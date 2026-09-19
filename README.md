# 化学方程式精确配平 API

基于 FastAPI 的化学方程式配平服务。服务端为反应物/生成物构建**元素守恒矩阵（含电荷行）**，用**有理数高斯消元**求零空间；仅当零空间为一维且基向量可化为全正时，返回最大公约数约简后的**最小正整数系数**及逐元素、电荷核算，否则返回 `422` 并说明原因（零空间维数或无法全正）。

另提供**化学计量接口**：在配平结果基础上，按内置原子量表（H、C、N、O、Na、Mg、Al、S、Cl、K、Ca、Fe、Cu、Mn）换算投料物质的量，以反应进度最小值判定（可并列的）限制反应物，返回理论反应进度、各生成物理论产量与各反应物剩余量。全程 `Decimal` 计算，输出按请求 `precision` 四舍五入。

以及**批次投料规划接口**：给定目标生成物的纯品质量与预计收率，反推理论反应进度，按配平系数求各反应物的化学计量摩尔数，再叠加每种原料的非负过量百分比、按纯度折算实际称量克数，返回配平方程、理论反应进度与逐原料投料规划。

以及**经验式推断接口**：由元素质量百分比反推经验式——先按内置原子量表换算摩尔比，再在 1～12 倍尺度内生成互质正整数下标候选，用候选式反算理论百分比并保留最大绝对偏差不超过 `tolerance` 者（并列时依次取总原子数较少、Hill 顺序字符串较小者）；提供实测摩尔质量时，按最接近的正整数倍给出分子式。

以及**氧化还原自动补全接口**：给定两侧已知物种（所在侧固定）与介质（`acidic`/`basic`），仅引入 H2O、H+（酸性）或 H2O、OH-（碱性）作为可出现在任一侧的辅助物种，求出元素与电荷均守恒的最小正整数净反应；多组最简解无法唯一判定或无法补全时返回 `422` 并说明原因。

以及**单反应平衡组成接口**：给出每个物种恰好一条初始摩尔数、体系体积与浓度平衡常数 Kc，复用精确配平得到化学计量数后，在所有物种摩尔数非负的反应进度区间内（同时允许正向 ξ>0 与逆向 ξ<0）用 `Decimal` 对 `ln Qc − ln Kc` 做有界求根（常规区域牛顿二次收敛、端点对数区域指数步跨越、异常退回纯二分），正确处理初始为零及端点趋于无穷的情形，返回反应方向、平衡反应进度、各物种平衡摩尔数与浓度、Qc 及对数残差。

## 快速启动（单条命令）

```bash
docker compose up --build
```

依赖随容器启动自动安装（`pip install -r requirements.txt` 后再拉起 uvicorn）。

- 服务地址：http://localhost:8000
- 交互式文档（Swagger UI）：http://localhost:8000/docs
- 健康检查：http://localhost:8000/health

## 配置

| 环境变量 | 默认值 | 说明 |
| --- | --- | --- |
| `PORT` | `8000` | 容器内服务监听端口 |

宿主机端口映射在 `docker-compose.yml` 的 `ports` 中调整（默认 `8000:8000`）。

## 分子式语法

| 语法 | 示例 |
| --- | --- |
| 元素与下标 | `H2O`、`Fe2O3` |
| 嵌套括号（`()` `[]` `{}` 可混用） | `Ca(OH)2`、`K4[Fe(CN)6]` |
| 点号结晶水（`.` `·` `*`） | `CuSO4.5H2O`、`Na2CO3·10H2O` |
| 插入符电荷（推荐，无歧义） | `SO4^2-`、`Fe^3+`、`NH4^+` |
| 纯符号尾电荷（每个 `+`/`-` 计 1） | `Na+`、`Cl-`、`Ca++`、`SO4--` |
| 电子（半反应） | `e-` |

> 注意：不带 `^` 时，符号前的数字属于前一个元素/基团——`NH4+` 解析为 NH₄⁺（电荷 +1）；要表示 Ca²⁺ 请写 `Ca^2+`。

## API

### `POST /api/v1/balance` — 配平

```bash
curl -s http://localhost:8000/api/v1/balance \
  -H 'Content-Type: application/json' \
  -d '{"reactants": ["KMnO4", "HCl"], "products": ["KCl", "MnCl2", "H2O", "Cl2"]}'
```

响应（`200`）：

```json
{
  "equation": "2KMnO4 + 16HCl -> 2KCl + 2MnCl2 + 8H2O + 5Cl2",
  "coefficients": [2, 16, 2, 2, 8, 5],
  "reactant_coefficients": [2, 16],
  "product_coefficients": [2, 2, 8, 5],
  "elements": [
    {"element": "Cl", "reactants": 16, "products": 16, "balanced": true},
    {"element": "H",  "reactants": 16, "products": 16, "balanced": true},
    {"element": "K",  "reactants": 2,  "products": 2,  "balanced": true},
    {"element": "Mn", "reactants": 2,  "products": 2,  "balanced": true},
    {"element": "O",  "reactants": 8,  "products": 8,  "balanced": true}
  ],
  "charge": {"reactants": 0, "products": 0, "balanced": true},
  "balanced": true
}
```

离子方程式（含电荷行守恒）：

```bash
curl -s http://localhost:8000/api/v1/balance \
  -H 'Content-Type: application/json' \
  -d '{"reactants": ["MnO4-", "H+", "e-"], "products": ["Mn^2+", "H2O"]}'
# -> MnO4- + 8H+ + 5e- -> Mn^2+ + 4H2O
```

无法配平时返回 `422`，并说明零空间维数或无法全正：

```json
{
  "detail": {
    "message": "零空间维数为 2（必须恰好为 1 才能唯一配平）",
    "null_space_dim": 2
  }
}
```

化学式语法错误返回 `400`。

### `POST /api/v1/verify` — 校验系数

```bash
curl -s http://localhost:8000/api/v1/verify \
  -H 'Content-Type: application/json' \
  -d '{
        "reactants": ["H2", "O2"],
        "products": ["H2O"],
        "reactant_coefficients": [2, 1],
        "product_coefficients": [2]
      }'
```

响应包含逐元素与电荷核算、`balanced` 结论，以及系数的 `gcd` 与是否已约简（`is_reduced`）。

### `POST /api/v1/redox/balance` — 氧化还原自动补全

给定两侧已知物种（**所在侧固定**，配平时不移动）与介质 `medium`（`acidic` 或 `basic`），仅引入辅助物种补齐守恒项，返回元素与电荷均守恒的**最小正整数净反应**：

- 酸性（`acidic`）：辅助物种为 H2O、H+；碱性（`basic`）：辅助物种为 H2O、OH-（等价于先按酸性配平、再把 H+ 按 `H+ + OH- -> H2O` 折算并约消两侧同种辅助物）。
- 已知物种中已有的辅助物种（如已把 H2O 写在某一侧）不会重复引入；辅助物种系数为零时不出现在结果中。
- `added_species` 列出实际新增的辅助物种及其所在侧与系数；`reactants`/`products` 为补全后的完整两侧物种列表。

```bash
curl -s http://localhost:8000/api/v1/redox/balance \
  -H 'Content-Type: application/json' \
  -d '{"reactants": ["MnO4-", "Fe^2+"], "products": ["Mn^2+", "Fe^3+"], "medium": "acidic"}'
```

响应（`200`）：

```json
{
  "medium": "acidic",
  "equation": "MnO4- + 5Fe^2+ + 8H+ -> Mn^2+ + 5Fe^3+ + 4H2O",
  "reactants": ["MnO4-", "Fe^2+", "H+"],
  "products": ["Mn^2+", "Fe^3+", "H2O"],
  "coefficients": [1, 5, 8, 1, 5, 4],
  "reactant_coefficients": [1, 5, 8],
  "product_coefficients": [1, 5, 4],
  "added_species": [
    {"formula": "H2O", "side": "products", "coefficient": 4},
    {"formula": "H+", "side": "reactants", "coefficient": 8}
  ],
  "elements": [
    {"element": "Fe", "reactants": 5, "products": 5, "balanced": true},
    {"element": "H",  "reactants": 8, "products": 8, "balanced": true},
    {"element": "Mn", "reactants": 1, "products": 1, "balanced": true},
    {"element": "O",  "reactants": 4, "products": 4, "balanced": true}
  ],
  "charge": {"reactants": 17, "products": 17, "balanced": true},
  "balanced": true
}
```

碱性介质（H+ 折算为 H2O/OH- 并约消）：

```bash
curl -s http://localhost:8000/api/v1/redox/balance \
  -H 'Content-Type: application/json' \
  -d '{"reactants": ["MnO4-", "SO3^2-"], "products": ["MnO2", "SO4^2-"], "medium": "basic"}'
# -> 2MnO4- + 3SO3^2- + H2O -> 2MnO2 + 3SO4^2- + 2OH-
```

错误返回：

| 状态码 | 场景 |
| --- | --- |
| `400` | 化学式无法解析 |
| `422` | 存在多组最简解无法唯一判定（`null_space_dim` > 1）；守恒方程组无非零解（`null_space_dim` = 0）；已知物种系数异号或为零，无法在不改变其所在侧的前提下补全；`medium` 不是 `acidic`/`basic` |

```json
{
  "detail": {
    "message": "零空间维数为 2（必须恰好为 1 才能唯一判定），存在多组最简解，无法唯一补全",
    "null_space_dim": 2
  }
}
```

### `POST /api/v1/stoichiometry` — 化学计量

在配平结果基础上，根据各反应物投料（单位仅支持 `mol` 或 `g`）计算限制反应物、理论反应进度、各生成物理论产量与各反应物剩余量。

- 内置原子量表（g/mol）：H 1.008、C 12.011、N 14.007、O 15.999、Na 22.98977、Mg 24.305、Al 26.982、S 32.06、Cl 35.45、K 39.0983、Ca 40.078、Fe 55.845、Cu 63.546、Mn 54.938（离子电荷的电子质量忽略不计）。
- 每种反应物必须恰好投料一次；以反应进度 ξ = n/ν 的最小值判定限制反应物（可并列多个）。
- 全程 `Decimal` 计算，所有小数字段按 `precision`（2～8，默认 4）四舍五入为**定长小数字符串**返回。

```bash
curl -s http://localhost:8000/api/v1/stoichiometry \
  -H 'Content-Type: application/json' \
  -d '{
        "reactants": ["H2", "O2"],
        "products": ["H2O"],
        "feed": [
          {"formula": "H2", "amount": 10,  "unit": "mol"},
          {"formula": "O2", "amount": 100, "unit": "g"}
        ],
        "precision": 4
      }'
```

响应（`200`）：

```json
{
  "equation": "2H2 + O2 -> 2H2O",
  "coefficients": [2, 1, 2],
  "reactant_coefficients": [2, 1],
  "product_coefficients": [2],
  "precision": 4,
  "extent": "3.1252",
  "limiting_reactants": ["O2"],
  "products": [
    {"formula": "H2O", "coefficient": 2, "theoretical_moles": "6.2504", "theoretical_grams": "112.6008"}
  ],
  "reactants": [
    {"formula": "H2", "coefficient": 2, "feed_moles": "10.0000", "feed_grams": "20.1600",
     "consumed_moles": "6.2504", "remaining_moles": "3.7496", "remaining_grams": "7.5592", "limiting": false},
    {"formula": "O2", "coefficient": 1, "feed_moles": "3.1252", "feed_grams": "100.0000",
     "consumed_moles": "3.1252", "remaining_moles": "0.0000", "remaining_grams": "0.0000", "limiting": true}
  ]
}
```

错误返回：

| 状态码 | 场景 |
| --- | --- |
| `400` | 化学式无法解析；缺失投料；同一物质重复投料；投料量非正数；投料物质不是反应物；元素无内置原子量 |
| `422` | 方程无法配平（零空间维数不为 1 或无法全正）；单位不是 `mol`/`g`；`precision` 超出 2～8 |

### `POST /api/v1/stoichiometry/plan` — 批次投料规划

由目标生成物的纯品质量与预计收率反推各原料称量：先配平方程，再由 `target_mass / yield_fraction` 得目标物理论产量、除以配平系数得**理论反应进度 ξ**；然后按配平系数求每种反应物的**化学计量摩尔数**，叠加该原料的**过量百分比**（`excess_percent`，非负）得到过量后摩尔数，最后按**纯度**（`purity`，(0,1] 质量分数）折算**实际称量克数**（`weigh_grams = 过量后摩尔数 × 摩尔质量 / 纯度`）。

- `target_product` 必须是 `products` 之一；`target_mass`（g）须为正数；`yield_fraction` 须在 (0,1] 内。
- 每种反应物必须在 `materials` 中恰好配置一条（纯度 + 过量百分比）；原子量表同化学计量接口。
- 全程 `Decimal` 计算，所有小数字段按 `precision`（2～8，默认 4）四舍五入为**定长小数字符串**返回。

```bash
curl -s http://localhost:8000/api/v1/stoichiometry/plan \
  -H 'Content-Type: application/json' \
  -d '{
        "reactants": ["H2", "O2"],
        "products": ["H2O"],
        "target_product": "H2O",
        "target_mass": "36.03",
        "yield_fraction": "0.9",
        "materials": [
          {"formula": "H2", "purity": "0.98", "excess_percent": "10"},
          {"formula": "O2", "purity": "1",    "excess_percent": "0"}
        ],
        "precision": 4
      }'
```

响应（`200`）：

```json
{
  "equation": "2H2 + O2 -> 2H2O",
  "coefficients": [2, 1, 2],
  "reactant_coefficients": [2, 1],
  "product_coefficients": [2],
  "precision": 4,
  "target_product": "H2O",
  "target_mass": "36.03",
  "yield_fraction": "0.9",
  "target_theoretical_moles": "2.2222",
  "target_theoretical_grams": "40.0333",
  "extent": "1.1111",
  "materials": [
    {"formula": "H2", "coefficient": 2, "purity": "0.98", "excess_percent": "10",
     "stoichiometric_moles": "2.2222", "moles_with_excess": "2.4444", "weigh_grams": "5.0286"},
    {"formula": "O2", "coefficient": 1, "purity": "1", "excess_percent": "0",
     "stoichiometric_moles": "1.1111", "moles_with_excess": "1.1111", "weigh_grams": "35.5533"}
  ]
}
```

错误返回：

| 状态码 | 场景 |
| --- | --- |
| `400` | 化学式无法解析；目标物不在生成物列表；`target_mass` 非正数；`yield_fraction` 不在 (0,1]；`purity` 不在 (0,1]；`excess_percent` 为负；原料配置缺失/重复/不是反应物；元素无内置原子量 |
| `422` | 方程无法配平（零空间维数不为 1 或无法全正）；`precision` 超出 2～8 |

### `POST /api/v1/formula/infer` — 经验式/分子式推断

由元素质量百分比反推经验式：按内置原子量表把百分比换算为摩尔比，在 **1～12 倍尺度**内四舍五入并约简为**互质正整数下标**候选；用候选式反算各元素理论百分比，保留**最大绝对偏差 ≤ `tolerance`** 的候选，偏差最小者胜出（并列时依次取总原子数较少、Hill 顺序字符串较小者）。提供 `molar_mass` 时，取**最接近实测值的正整数倍**作为分子式，其相对偏差（%）同样须 ≤ `tolerance`。

- 百分比须均为正数，且总和与 100 的偏差不超过 `tolerance`（默认 `0.5`，单位：百分点）。
- 元素须在内置原子量表内（同化学计量接口）。
- 所有小数字段按 `precision`（2～8，默认 4）四舍五入为**定长小数字符串**返回。

```bash
curl -s http://localhost:8000/api/v1/formula/infer \
  -H 'Content-Type: application/json' \
  -d '{
        "percentages": {"C": 40.0, "H": 6.7, "O": 53.3},
        "molar_mass": 180.0
      }'
```

响应（`200`）：

```json
{
  "empirical_formula": "CH2O",
  "subscripts": {"C": 1, "H": 2, "O": 1},
  "theoretical_percentages": {"C": "40.0020", "H": "6.7142", "O": "53.2838"},
  "deviations": {"C": "0.0020", "H": "0.0142", "O": "0.0162"},
  "max_deviation": "0.0162",
  "empirical_molar_mass": "30.0260",
  "molecular_formula": "C6H12O6",
  "multiplier": 6,
  "molecular_molar_mass": "180.1560",
  "molar_mass_relative_deviation": "0.0867"
}
```

不提供 `molar_mass` 时，分子式相关字段为 `null`：

```bash
curl -s http://localhost:8000/api/v1/formula/infer \
  -H 'Content-Type: application/json' \
  -d '{"percentages": {"N": 36.86, "O": 63.14}, "tolerance": 0.5}'
# -> empirical_formula: N2O3（尺度 2 的候选比尺度 1 的 NO2 偏差更小）
```

错误返回：

| 状态码 | 场景 |
| --- | --- |
| `400` | 百分比非正数；百分比总和与 100 的偏差超过 `tolerance`；元素无内置原子量；`molar_mass` 非正数 |
| `422` | 1～12 倍尺度内无候选式满足容差；摩尔质量无满足相对偏差容差的正整数倍；`tolerance` 为负；`precision` 超出 2～8 |

`422` 响应的 `detail` 附带诊断信息（最接近的候选式及其偏差，或摩尔质量倍数与相对偏差），例如：

```json
{
  "detail": {
    "message": "1～12 倍尺度内无候选式满足最大绝对偏差 ≤ tolerance=0.5",
    "tolerance": "0.5",
    "nearest_candidate": "CO",
    "nearest_max_deviation": "0.8811"
  }
}
```

### `POST /api/v1/equilibrium` — 单反应平衡组成

给定每个物种**恰好一条**初始摩尔数（允许为 `0`、不允许为负）、体系体积 `volume`（>0）与浓度平衡常数 `kc`（>0），计算单反应的平衡组成：

- 先复用精确配平得到最小正整数化学计量数；反应物的有符号计量数取 `−ν`、生成物取 `+ν`，组成按 `n_i(ξ) = n_i0 + ν_i·ξ` 随反应进度 ξ 变化。
- 由所有物种摩尔数非负确定可行区间 `[ξ_min, ξ_max]`：`ξ_max = min_{反应物} n_i0/ν_i`，`ξ_min = max_{生成物} −n_i0/ν_i`；区间覆盖正负两个方向（`ξ>0` 正向、`ξ<0` 逆向）。
- `f(ξ) = ln Qc − ln Kc` 在开区间内严格单调（`f'(ξ) = Σ ν_i²/n_i > 0`），端点处零摩尔物种使 `Qc → 0/∞`、`f → ∓∞`，故平衡根唯一且必在内部。全程 `Decimal` 计算（内部精度随 `Kc`/体积/投料量级自适应，至少 50 位）：常规区域用牛顿步（二次收敛），端点对数区域用指数步单步跨越数十个数量级，候选点始终严格落在括号区间内，数值异常时退回有界二分。
- 反应方向 `direction`：`forward`（Qc<Kc，正向）、`reverse`（Qc>Kc，逆向）、`equilibrium`（初始组成恰好满足 Qc=Kc，ξ=0）。初始缺反应物（Qc=0）只能逆向、初始缺生成物（Qc=∞）只能正向；若**两侧都存在零初始物种**，可行区间退化为单点，正向、逆向均无法启动，返回 `422`。
- 所有小数字段按 `precision`（2～8，默认 4）四舍五入为**定长小数字符串**返回；`residual` 为对数残差 `ln(Qc/Kc)`（平衡时为 `0.0000`）。

```bash
curl -s http://localhost:8000/api/v1/equilibrium \
  -H 'Content-Type: application/json' \
  -d '{
        "reactants": ["N2O4"],
        "products": ["NO2"],
        "initial_moles": [
          {"formula": "N2O4", "moles": 1},
          {"formula": "NO2",  "moles": 0}
        ],
        "volume": 2,
        "kc": "0.1",
        "precision": 4
      }'
```

响应（`200`）：

```json
{
  "equation": "N2O4 -> 2NO2",
  "coefficients": [1, 2],
  "reactant_coefficients": [1],
  "product_coefficients": [2],
  "precision": 4,
  "volume": "2",
  "kc": "0.1",
  "direction": "forward",
  "extent": "0.2000",
  "extent_min": "0.0000",
  "extent_max": "1.0000",
  "qc": "0.1000",
  "residual": "0.0000",
  "species": [
    {"formula": "N2O4", "side": "reactants", "coefficient": 1,
     "initial_moles": "1.0000", "equilibrium_moles": "0.8000", "concentration": "0.4000"},
    {"formula": "NO2",  "side": "products",  "coefficient": 2,
     "initial_moles": "0.0000", "equilibrium_moles": "0.4000", "concentration": "0.2000"}
  ]
}
```

逆向进行的例子（初始只有生成物、Qc=0）：

```bash
curl -s http://localhost:8000/api/v1/equilibrium \
  -H 'Content-Type: application/json' \
  -d '{
        "reactants": ["H2", "I2"],
        "products": ["HI"],
        "initial_moles": [
          {"formula": "H2", "moles": 0},
          {"formula": "I2", "moles": 0},
          {"formula": "HI", "moles": 2}
        ],
        "volume": 1,
        "kc": 49,
        "precision": 6
      }'
# H2 + I2 ⇌ 2HI：ξ = -2/9（逆向），H2=I2=0.222222，HI=1.555556，Qc=49
```

错误返回：

| 状态码 | 场景 |
| --- | --- |
| `400` | 化学式无法解析；初始摩尔数缺失/重复/为负/物种不在方程中；所有物种初始摩尔数均为 0；`volume` 或 `kc` 非正数 |
| `422` | 方程无法配平（零空间维数不为 1 或无法全正）；可行反应进度区间退化为单点（初始为零的反应物与生成物同时存在）；`precision` 超出 2～8 |

`422`（无可行区间）响应的 `detail` 附带诊断信息（零初始物种与区间端点）：

```json
{
  "detail": {
    "message": "可行反应进度区间退化为单点 ξ=0：反应物与生成物中均存在初始摩尔数为 0 的物种，正向与逆向都无法启动，不存在满足 Qc=Kc 的非平凡平衡组成",
    "zero_reactants": ["I2"],
    "zero_products": ["HI"],
    "extent_min": "0.0000",
    "extent_max": "0.0000"
  }
}
```

## 本地开发

```bash
pip install -r requirements-dev.txt
uvicorn app.main:app --reload --port 8000
pytest                 # 运行测试
```

## 项目结构

```
app/
  parser.py         # 分子式解析：嵌套括号、结晶水、离子电荷
  balancer.py       # 守恒矩阵（含电荷行）、有理数 RREF、零空间、整数化
  redox.py          # 氧化还原自动补全：H2O/H+/OH- 辅助列、介质折算、净反应组装
  stoichiometry.py  # 内置原子量表、投料换算、限制反应物、产量与剩余量、批次投料规划（Decimal）
  formula_infer.py  # 质量百分比 → 摩尔比 → 互质下标候选 → 经验式/分子式推断（Decimal）
  equilibrium.py    # 可行反应进度区间 → ln Qc-ln Kc 有界求根（牛顿/指数步/二分，Decimal）
  schemas.py        # 请求/响应模型
  main.py           # FastAPI 路由
tests/test_api.py / tests/test_redox.py / tests/test_stoichiometry.py / tests/test_stoichiometry_plan.py / tests/test_formula_infer.py / tests/test_equilibrium.py
Dockerfile / docker-compose.yml
```

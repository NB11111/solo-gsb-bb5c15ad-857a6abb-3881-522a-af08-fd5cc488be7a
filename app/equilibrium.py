"""单反应平衡组成：给定初始摩尔数、体积与 Kc，求平衡反应进度。

流程：
  1. 复用 balancer.balance_equation 得到精确配平的化学计量数；
     反应物的有符号计量数取 -ν_i、生成物取 +ν_i，平衡组成按
     n_i(ξ) = n_i0 + ν_i·ξ 随反应进度 ξ 变化（ξ>0 正向，ξ<0 逆向）；
  2. 校验每个物种恰好一条非负初始摩尔数，体积与 Kc 为正数；
  3. 由 n_i(ξ) ≥ 0 求可行区间 [ξ_min, ξ_max]：
       ξ_max = min_{ν_i<0} n_i0/(-ν_i)，ξ_min = max_{ν_i>0} -n_i0/ν_i；
     区间退化为单点（零初始的反应物与生成物同时存在，正向、逆向
     都无法启动）时无可平衡组成，抛 NoEquilibriumError（HTTP 422）；
  4. f(ξ) = ln Qc(ξ) - ln Kc = Σ ν_i·ln(n_i/V) - Δν·ln V - ln Kc
     在开区间内严格单调（f'(ξ) = Σ ν_i²/n_i > 0）；端点处某物种
     摩尔数为零使 Qc 趋于 0 或 +∞（f 分别趋于 -∞/+∞），故零点
     唯一且必在内部。常规区域用牛顿步（d=-f/f'，二次收敛），端点
     附近 f 由 P·ln δ 主导（δ 为到约束边界的距离，牛顿步退化为
     线性收敛）时改用指数步 δ→δ·exp(-|f|/(f'·δ))，单步可跨越
     数百个十进制量级；候选点始终严格落在括号区间 (ξ_min, ξ_max)
     内（初始为零的物种只在端点取零，不会出现 ln 0），数值异常时
     退回纯二分，直到 |f| 小于输出精度容差；
  5. 返回反应方向、平衡反应进度、各物种平衡摩尔数与浓度、Qc 及
     对数残差 ln(Qc/Kc)。

全程 Decimal 计算（内部精度按 |log Kc|、|log V| 与投料量级自适应，
至少 50 位），输出统一按 precision 四舍五入（ROUND_HALF_UP）为定长
小数字符串。输入数据本身不合法抛 EquilibriumError（HTTP 400）；
无正宽度可行区间抛 NoEquilibriumError（HTTP 422）；
FormulaError / BalanceError 由路由层按既有约定映射。
"""
from __future__ import annotations

from decimal import Decimal, localcontext
from typing import TYPE_CHECKING

from .balancer import balance_equation
from .stoichiometry import _round

if TYPE_CHECKING:  # 避免在逻辑层硬依赖 API 模型
    from .schemas import InitialMolesItem

# 二分停止策略：直至 |f(ξ)| 小于输出精度容差（或中点因有效数字耗尽与
# 端点重合）；迭代上限按内部精度自适应，另设硬上限防止病态输入无限循环。
_MAX_BISECTION_ITERATIONS = 20000
# 指数步单步最多跨越的十进制量级（过大会越过根，过小则需要多步）
_LOG_JUMP = 250
_LOG_JUMP_DECIMAL = Decimal(_LOG_JUMP)
# 牛顿相对步长低于该比例时判定进入端点对数区域，改用指数步
_NEWTON_MIN_FRACTION = Decimal("0.1")

_ZERO = Decimal(0)
_ONE = Decimal(1)


def _format_fixed(value: Decimal, precision: int) -> str:
    """定长小数格式化；对极大/极小值临时放大上下文精度，避免 quantize 溢出。"""
    with localcontext() as ctx:
        # quantize 为定点数所需的系数位数 ≈ 最高位指数 + precision
        ctx.prec = max(50, value.adjusted() + precision + 10)
        return _round(value + _ZERO, precision)


class EquilibriumError(ValueError):
    """初始摩尔数 / 体积 / Kc 等输入数据不合法（对应 HTTP 400）。"""


class NoEquilibriumError(ValueError):
    """不存在正宽度的可行反应进度区间（对应 HTTP 422）。"""

    def __init__(self, message: str, **context: object):
        super().__init__(message)
        self.message = message
        self.context = context


def compute_equilibrium(
    reactants: list[str],
    products: list[str],
    initial_moles: list[InitialMolesItem],
    volume: Decimal,
    kc: Decimal,
    precision: int,
) -> dict:
    """主流程：配平 → 投料与参数校验 → 可行区间 → Decimal 有界二分。"""
    balanced = balance_equation(reactants, products)
    r_coeffs = balanced["reactant_coefficients"]
    p_coeffs = balanced["product_coefficients"]

    formulas = list(reactants) + list(products)
    n_reactants = len(reactants)
    # 有符号化学计量数：反应物消耗（负）、生成物生成（正）
    nu = [Decimal(-c) for c in r_coeffs] + [Decimal(c) for c in p_coeffs]
    delta_nu = sum(nu)  # 气体计量数差 Δν = Σ生成物ν - Σ反应物ν

    # ---- 初始摩尔数校验：重复 / 非本方程物种 / 负数 / 缺失 ----
    n0: dict[str, Decimal] = {}
    for item in initial_moles:
        formula = item.formula.strip()
        if formula in n0:
            raise EquilibriumError(f"同一物种重复给出初始摩尔数: {formula!r}")
        if formula not in formulas:
            raise EquilibriumError(f"初始摩尔数物种 {formula!r} 不在反应方程式中")
        if not item.moles.is_finite() or item.moles < 0:
            raise EquilibriumError(
                f"初始摩尔数必须为非负数: {formula!r} 的 moles 为 {item.moles}"
            )
        n0[formula] = item.moles
    missing = [f for f in formulas if f not in n0]
    if missing:
        raise EquilibriumError(f"缺少物种初始摩尔数: {', '.join(missing)}")

    # ---- 体积与 Kc 校验 ----
    if not volume.is_finite() or volume <= 0:
        raise EquilibriumError(f"体积必须为正数: volume 为 {volume}")
    if not kc.is_finite() or kc <= 0:
        raise EquilibriumError(f"Kc 必须为正数: kc 为 {kc}")
    if all(n0[f] == 0 for f in formulas):
        raise EquilibriumError("至少需要一种物种的初始摩尔数为正数，不能全部为 0")

    initial = [n0[f] for f in formulas]

    # ---- 自适应内部精度 ----
    # 根可能距端点约 |Kc|^(-1/p)（零初始极限反应），有效数字需求随
    # |log Kc|、|log V| 及投料量级增长；再留 40 位保护位保证残差小。
    ln10_rough = Decimal(10).ln()

    def _magnitude(x: Decimal) -> int:
        return int((x.ln() / ln10_rough).copy_abs()) + 1

    magnitudes = [_magnitude(kc), _magnitude(volume)]
    magnitudes.extend(_magnitude(x) for x in initial if x > 0)
    inner_prec = max(50, max(magnitudes) + 40)

    with localcontext() as ctx:
        # 内部精度远高于输出精度（2～8 位），避免二分与对数的中间舍入误差
        ctx.prec = inner_prec
        ln10 = Decimal(10).ln()

        # ---- 可行区间 [xi_min, xi_max]：所有 n_i(ξ) ≥ 0 ----
        xi_max = None
        xi_min = None
        for n_i, v_i in zip(initial, nu):
            if v_i < 0:  # 反应物：n0 - c·ξ ≥ 0 ⇒ ξ ≤ n0/c
                bound = n_i / (-v_i)
                xi_max = bound if xi_max is None else min(xi_max, bound)
            else:  # 生成物：n0 + c·ξ ≥ 0 ⇒ ξ ≥ -n0/c
                bound = -n_i / v_i
                xi_min = bound if xi_min is None else max(xi_min, bound)
        # reactants/products 均非空（请求模型 min_length=1，配平系数全正），
        # 故两端点必为有限值且 xi_min ≤ 0 ≤ xi_max。
        if not xi_min < xi_max:
            zero_reactants = [
                form for form, v_i in zip(formulas, nu) if v_i < 0 and n0[form] == 0
            ]
            zero_products = [
                form for form, v_i in zip(formulas, nu) if v_i > 0 and n0[form] == 0
            ]
            raise NoEquilibriumError(
                "可行反应进度区间退化为单点 ξ=0：反应物与生成物中均存在"
                "初始摩尔数为 0 的物种，正向与逆向都无法启动，"
                "不存在满足 Qc=Kc 的非平凡平衡组成",
                zero_reactants=zero_reactants,
                zero_products=zero_products,
                extent_min=_round(xi_min, precision),
                extent_max=_round(xi_max, precision),
            )

        log_kc = kc.ln()
        log_volume = volume.ln()

        def f_and_df(xi: Decimal) -> tuple[Decimal, Decimal]:
            """f(ξ)=ln Qc-ln Kc 与 f'(ξ)=Σ ν_i²/n_i。要求 ξ 严格在区间内部。"""
            log_q = _ZERO
            df = _ZERO
            for n_i0, v_i in zip(initial, nu):
                n_i = n_i0 + v_i * xi
                log_q += v_i * n_i.ln()
                df += v_i * v_i / n_i
            return log_q - delta_nu * log_volume - log_kc, df

        def f(xi: Decimal) -> Decimal:
            """f(ξ)=ln Qc-ln Kc（仅取函数值，残差复核时使用）。"""
            return f_and_df(xi)[0]

        # ---- 反应方向判定（严格，不依赖二分近似与对数舍入） ----
        # 端点情形：某侧存在零初始摩尔物种时 Qc 严格为 0 或 +∞。
        zero_reactants = [
            form for form, v_i in zip(formulas, nu) if v_i < 0 and n0[form] == 0
        ]
        zero_products = [
            form for form, v_i in zip(formulas, nu) if v_i > 0 and n0[form] == 0
        ]
        if zero_reactants:
            # Qc=0（初始缺反应物，分子为 0）：只能逆向进行
            direction = "reverse"
        elif zero_products:
            # Qc=+∞（初始缺生成物，分母为 0）：只能正向进行
            direction = "forward"
        else:
            # Qc0 用整数次幂乘积精确计算（不经对数），与 Kc 直接比较，
            # 避免 Decimal.ln 的末位噪声把精确平衡态误判为有方向。
            qc0 = _ONE
            for n_i0, v_i in zip(initial, nu):
                qc0 *= (n_i0 / volume) ** int(v_i)
            if qc0 == kc:
                direction = "equilibrium"
            elif qc0 < kc:
                direction = "forward"
            else:
                direction = "reverse"

        # ---- 有界求根：f 严格递增，f(xi_min)→-∞、f(xi_max)→+∞ ----
        # 停止条件：|f(ξ)| 小于输出精度的 1/4 单位——此时残差四舍五入
        # 必为 0，Qc=Kc·e^res 与 Kc 的相对偏差远小于显示精度，各摩尔数
        # 同样已准确（端点附近的小量物种 n∝e^{f/p}，相对误差≈|f|/p）。
        residual_tol = Decimal(1).scaleb(-(precision + 2))
        if direction == "equilibrium":
            xi = _ZERO
            residual = _ZERO
        else:
            lo, hi = xi_min, xi_max
            xi = (lo + hi) / 2
            residual, slope = f_and_df(xi)
            # 迭代预算：纯二分每轮区间减半，极端 Kc 下需 ~3.32·量级轮；
            # 指数/牛顿加速通常远少于此。硬上限兜底。
            iterations = min(
                _MAX_BISECTION_ITERATIONS,
                max(100, (max(magnitudes) * 10) // 3 + 100),
            )
            # 步长策略（始终保证候选点严格位于括号 (lo, hi) 内）：
            #   * 常规区域：牛顿步 d=-f/f' 落在括号内部且相对朝根
            #     方向边界的距离不太小（≥ 10%）时直接采用，二次收敛；
            #   * 端点对数区域：朝根方向的边界处有物种摩尔数趋零，
            #     f≈∓P·ln(δ/δ*)（δ 为当前到该边界的距离）。此时
            #     牛顿步只移动 δ·|f|/P，比例虽大但模型应直接给出
            #     新距离 δ* = δ·exp(-|f|/(f'·δ))，候选点朝根移动
            #     δ-δ*（单步跨越 capped 在 _LOG_JUMP 个十进制量级）；
            #   * 数值异常时退回纯二分。
            for _ in range(iterations):
                if residual.copy_abs() <= residual_tol:
                    break
                newton = xi - residual / slope
                if residual > 0:  # 根在左侧，朝 lo 方向
                    d_boundary = xi - lo
                else:  # 根在右侧，朝 hi 方向
                    d_boundary = hi - xi
                newton_step = (newton - xi).copy_abs()
                if lo < newton < hi and newton_step >= _NEWTON_MIN_FRACTION * d_boundary:
                    candidate = newton
                else:
                    p_effective = slope * d_boundary  # 趋零物种的 |ν| 之和
                    log_shrink = residual.copy_abs() / p_effective
                    log_shrink = min(log_shrink, _LOG_JUMP_DECIMAL * ln10)
                    d_star = d_boundary * (-log_shrink).exp()
                    travel = d_boundary - d_star
                    if residual > 0:
                        candidate = xi - travel
                    else:
                        candidate = xi + travel
                    if not lo < candidate < hi:
                        candidate = (lo + hi) / 2
                    if candidate == lo or candidate == hi:
                        break  # 已达到当前有效数字极限
                value, candidate_slope = f_and_df(candidate)
                if value > 0:
                    hi = candidate
                elif value < 0:
                    lo = candidate
                else:
                    xi, residual, slope = candidate, value, candidate_slope
                    break
                xi, residual, slope = candidate, value, candidate_slope

        # ---- 组装平衡组成（防御性截断，抵消极端参数下的端点舍入） ----
        equilibrium_moles = [
            max(n_i0 + v_i * xi, _ZERO) for n_i0, v_i in zip(initial, nu)
        ]
        concentrations = [n_i / volume for n_i in equilibrium_moles]

        # 用收敛后的 ξ 重算残差与 Qc（Qc = Kc·exp(residual)），避免使用
        # 迭代过程中的陈旧值；退化到端点（有限精度下 n_j 为 0）时 Qc
        # 无法表示，按 Kc 报出、残差记 0。
        if all(n_i > 0 for n_i in equilibrium_moles):
            residual = f(xi)
            qc = (log_kc + residual).exp()
        else:
            residual = _ZERO
            qc = kc

        species_rows = [
            {
                "formula": formula,
                "side": "reactants" if i < n_reactants else "products",
                "coefficient": int(abs(v_i)),
                "initial_moles": _format_fixed(initial[i], precision),
                "equilibrium_moles": _format_fixed(equilibrium_moles[i], precision),
                "concentration": _format_fixed(concentrations[i], precision),
            }
            for i, (formula, v_i) in enumerate(zip(formulas, nu))
        ]

        return {
            "equation": balanced["equation"],
            "coefficients": balanced["coefficients"],
            "reactant_coefficients": r_coeffs,
            "product_coefficients": p_coeffs,
            "precision": precision,
            "volume": str(volume),
            "kc": str(kc),
            "direction": direction,
            "extent": _format_fixed(xi, precision),
            "extent_min": _format_fixed(xi_min, precision),
            "extent_max": _format_fixed(xi_max, precision),
            "qc": _format_fixed(qc, precision),
            "residual": _format_fixed(residual, precision),
            "species": species_rows,
        }

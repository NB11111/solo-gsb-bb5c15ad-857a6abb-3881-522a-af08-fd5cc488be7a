"""单反应平衡组成：精确配平 + Decimal 有界二分求平衡反应进度。

流程：
  1. 复用 balancer 的精确配平，得到最小正整数化学计量数
     （采用标准约定 ν_r<0、ν_p>0，n_i(ξ) = n_i0 + ν_i·ξ）；
  2. 校验每个物种恰好一条非负初始摩尔数、体积与 Kc 均为正；
  3. 由非负约束求可行反应进度区间 [ξ_min, ξ_max]
     （ξ_min = max_{生成物} -n0_i/ν_i ≤ 0，
      ξ_max = min_{反应物} n0_i/(-ν_i) ≥ 0）；
  4. 区间退化为一点（ξ_min == ξ_max，即两侧都存在初始为零的物种）
     时 Qc 无定义，抛 NoFeasibleIntervalError（对应 HTTP 422）；
  5. 在 [ξ_min, ξ_max] 上对 f(ξ) = ln Qc(ξ) - ln Kc 做有界二分：
     开区间内所有物种摩尔数严格为正，ln Qc 连续且严格单调递增
     （f'(ξ) = Σ ν_i²/n_i > 0），端点处按零摩尔数所在侧直接给出
     趋于 -∞/+∞ 的符号，因此初始为零与端点发散都无需取极限；
  6. 按根 ξ* 的符号给出反应方向（forward / reverse / at_equilibrium），
     返回各物种平衡摩尔数与浓度、Qc 与有符号残差。

全程 Decimal 计算（内部精度 200），输出按请求 precision 四舍五入
（ROUND_HALF_UP）为定长小数字符串；输入本身不合法抛 EquilibriumError
（对应 HTTP 400）。
"""
from __future__ import annotations

from decimal import Decimal, ROUND_HALF_UP, localcontext
from typing import TYPE_CHECKING

from .balancer import balance_equation
from .stoichiometry import _round

if TYPE_CHECKING:  # 避免在逻辑层硬依赖 API 模型
    from .schemas import InitialMolesItem

# 二分迭代上限：精度 200 时区间可收缩约 600 次到机器极限
_MAX_ITERATIONS = 600


class EquilibriumError(ValueError):
    """投料 / 体积 / Kc 等输入数据不合法（对应 HTTP 400）。"""


class NoFeasibleIntervalError(ValueError):
    """所有物种摩尔数非负的反应进度区间退化为一点，Qc 无定义（对应 HTTP 422）。"""

    def __init__(self, message: str, **context: object):
        super().__init__(message)
        self.message = message
        self.context = context


def _format_fixed(value: Decimal, precision: int, *, signed: bool = False) -> str:
    """按 precision 位小数四舍五入输出。

    常规值给定长小数字符串（负零归一化为正零）；定点形式会把非零值
    舍成零、或数值过大导致小数串过长时，退化为带 precision 位尾数
    小数的科学计数法（如 ``1.2346E-30``），避免输出 0.0000 式误导或
    超长数字串。signed=True 时正数也显式带 ``+``（用于有符号残差）。
    """
    quantum = Decimal(1).scaleb(-precision)
    try:
        rounded = value.quantize(quantum, rounding=ROUND_HALF_UP)
    except ArithmeticError:
        # 数值过大（如 1E200），定点表示所需位数超过上下文精度
        rounded = None

    def scientific(v: Decimal) -> str:
        exp = v.adjusted()
        mantissa = v.scaleb(-exp).quantize(quantum, rounding=ROUND_HALF_UP)
        # 尾数舍入进位时归一化到 [1, 10)（9.9999 舍成 10.0000）
        if abs(mantissa) >= 10:
            mantissa /= 10
            exp += 1
        sign = "+" if signed and v > 0 else ("-" if v < 0 else "")
        return f"{sign}{format(mantissa.copy_abs(), 'f')}E{exp:+03d}"

    if rounded is None:
        return scientific(value)
    if rounded == 0 and value == 0:
        text = format(rounded.copy_abs(), "f")
        return f"+{text}" if signed else text

    # 定点形式会把非零值舍成 0（如 1E-60 配 6 位小数），或数值过大使
    # 小数串过长时，退化为带 precision 位尾数小数的科学计数法
    if rounded == 0 or rounded.adjusted() > precision + 8:
        return scientific(value)
    text = format(rounded, "f")
    return f"+{text}" if signed and rounded > 0 else text


def _format_residual(value: Decimal, precision: int) -> str:
    """有符号残差按 precision 位定点小数输出（不足一位即显示 ±0.0000）。

    与 Qc/Kc 不同，残差是求解容差指标，极小（如 1E-180）正说明已充分
    收敛，故不走科学计数法，统一四舍五入为定长小数并显式带号；
    极端未收敛的大残差退化为科学计数法。
    """
    quantum = Decimal(1).scaleb(-precision)
    try:
        rounded = value.quantize(quantum, rounding=ROUND_HALF_UP)
    except ArithmeticError:
        return _format_fixed(value, precision, signed=True)
    if rounded == 0:
        rounded = rounded.copy_abs()
    text = format(rounded, "f")
    if rounded >= 0:
        return f"+{text}"
    return text


def _safe_round(value: Decimal, precision: int) -> str:
    """定长小数四舍五入；数值过大导致定点表示越上下文精度时退回科学计数法。"""
    try:
        return _round(value, precision)
    except ArithmeticError:
        return _format_fixed(value, precision)


def _moles_at(
    xi: Decimal,
    formulas: list[str],
    initial: dict[str, Decimal],
    stoich: list[int],
) -> list[Decimal]:
    """n_i(ξ) = n_i0 + ν_i·ξ（反应物 ν<0 消耗，生成物 ν>0 生成）。"""
    return [initial[f] + Decimal(stoich[i]) * xi for i, f in enumerate(formulas)]


def _ln_qc(
    moles: list[Decimal],
    stoich: list[int],
    volume: Decimal,
) -> Decimal:
    """ln Qc = Σ ν_i·ln(n_i/V)。调用方保证各 n_i 严格为正。

    常数项 ln V 可提出：Σ ν_i·ln(n_i/V) = Σ ν_i·ln n_i - (Σ ν_i)·ln V。
    """
    total = Decimal(0)
    for n, nu in zip(moles, stoich):
        total += Decimal(nu) * n.ln()
    return total - Decimal(sum(stoich)) * volume.ln()


def _bisect(
    lo: Decimal,
    hi: Decimal,
    formulas: list[str],
    initial: dict[str, Decimal],
    stoich: list[int],
    volume: Decimal,
    ln_kc: Decimal,
) -> tuple[Decimal, Decimal]:
    """在 [lo, hi] 上对 f(ξ) = ln Qc(ξ) - ln Kc 做有界二分。

    非退化可行区间内 ln Qc 严格单调递增，端点符号恒为
    f(lo)<0（某生成物在 lo 处摩尔数为 0，ln Qc -> -∞）、
    f(hi)>0（某反应物在 hi 处摩尔数为 0，ln Qc -> +∞），
    故端点不实际求值，仅在开区间内取中点——中点处各物种摩尔数严格为正，
    可安全取 ln。返回 (ξ, f(ξ))；区间收缩到机器极限（中点不再变化）
    或迭代上限时，返回历次中点中 |f| 最小者。
    """
    best_xi = (lo + hi) / 2
    best_f = None
    for _ in range(_MAX_ITERATIONS):
        mid = (lo + hi) / 2
        if mid == lo or mid == hi:
            break  # Decimal 精度耗尽，区间已不可再分
        f_mid = _ln_qc(
            _moles_at(mid, formulas, initial, stoich), stoich, volume
        ) - ln_kc
        if best_f is None or abs(f_mid) < abs(best_f):
            best_xi, best_f = mid, f_mid
        if f_mid == 0:
            return mid, Decimal(0)
        if f_mid < 0:
            lo = mid
        else:
            hi = mid
    return best_xi, best_f if best_f is not None else Decimal(0)


def compute_equilibrium(
    reactants: list[str],
    products: list[str],
    initial_moles: list[InitialMolesItem],
    volume: Decimal,
    kc: Decimal,
    precision: int,
) -> dict:
    """主流程：配平 → 投料校验 → 可行区间 → 二分 → 组装响应。

    FormulaError / BalanceError 继续向上抛，由路由层映射为 4xx。
    """
    # 物种列表不允许重复（同一侧重复或跨两侧出现都会使初始摩尔数归属不唯一）
    formulas = list(reactants) + list(products)
    if len(set(formulas)) != len(formulas):
        dup = sorted({f for f in formulas if formulas.count(f) > 1})
        raise EquilibriumError(
            "物种在反应物/生成物列表中重复出现，初始摩尔数归属不唯一: "
            + ", ".join(dup)
        )

    balanced = balance_equation(reactants, products)
    r_coeffs = balanced["reactant_coefficients"]
    p_coeffs = balanced["product_coefficients"]
    # 有符号化学计量数（IUPAC 约定）：反应物 -ν、生成物 +ν
    stoich = [-c for c in r_coeffs] + list(p_coeffs)

    # ---- 全局参数校验 ----
    if not volume.is_finite() or volume <= 0:
        raise EquilibriumError(f"体积必须为正数: volume 为 {volume}")
    if not kc.is_finite() or kc <= 0:
        raise EquilibriumError(f"平衡常数 Kc 必须为正数: kc 为 {kc}")

    # ---- 初始摩尔数校验：重复 / 非物种 / 非法值 / 缺失（每物种恰好一条） ----
    initial: dict[str, Decimal] = {}
    for item in initial_moles:
        formula = item.formula.strip()
        if formula in initial:
            raise EquilibriumError(f"同一物种重复给出初始摩尔数: {formula!r}")
        if formula not in formulas:
            raise EquilibriumError(f"初始摩尔数物种 {formula!r} 不在反应物/生成物列表中")
        if not item.moles.is_finite() or item.moles < 0:
            raise EquilibriumError(
                f"初始摩尔数必须为非负数: {formula!r} 的 moles 为 {item.moles}"
            )
        initial[formula] = item.moles
    missing = [f for f in formulas if f not in initial]
    if missing:
        raise EquilibriumError(f"缺少初始摩尔数: {', '.join(missing)}")

    n_reactants = len(reactants)

    with localcontext() as ctx:
        # 内部精度远高于输出精度；兼顾根距端点极近（极端 Kc）时
        # 二分需要的收缩次数（每约 3.32 次迭代换取一位十进制）
        ctx.prec = 200
        # 放开指数上限，允许极大 Kc 时 exp(lnQc) 还原 Qc
        ctx.Emax = 999999
        ctx.Emin = -999999

        # ---- 可行反应进度区间 [xi_min, xi_max]（xi_min ≤ 0 ≤ xi_max） ----
        # n_i(ξ) = n0_i + ν_i·ξ ≥ 0：
        # 反应物 ν<0 => ξ ≤ -n0/ν = n0/|ν|（正向极限，≥0）；
        # 生成物 ν>0 => ξ ≥ -n0/ν（逆向极限，≤0）。
        xi_min = max(
            -initial[products[i]] / Decimal(p_coeffs[i])
            for i in range(len(products))
        )
        xi_max = min(
            initial[reactants[i]] / Decimal(r_coeffs[i])
            for i in range(n_reactants)
        )

        if xi_min >= xi_max:
            # 退化（仅当 xi_min == xi_max == 0）：某反应物与某生成物初始
            # 均为零，反应向任何方向移动都会使其中一个摩尔数为负，
            # 唯一可行点处 Qc 为 0/0，无定义
            zero_reactants = [
                reactants[i]
                for i in range(n_reactants)
                if initial[reactants[i]] == 0
            ]
            zero_products = [p for p in products if initial[p] == 0]
            raise NoFeasibleIntervalError(
                "满足所有物种摩尔数非负的反应进度区间退化为一点，"
                "该点处 Qc 无定义（反应物与生成物中均存在初始摩尔数为零的物种），"
                "无法计算平衡组成",
                xi_min=_safe_round(xi_min, precision),
                xi_max=_safe_round(xi_max, precision),
                zero_reactants=zero_reactants,
                zero_products=zero_products,
            )

        ln_kc = kc.ln()

        # ---- 初始点 ξ=0 处若各物种均为正，可先判定是否恰好平衡；
        # 否则（存在初始为零的物种）f(0) 发散，直接二分——端点符号
        # 恒为左 -∞、右 +∞，二分区间始终夹住唯一根 ----
        moles_zero = _moles_at(Decimal(0), formulas, initial, stoich)
        any_zero = any(n == 0 for n in moles_zero)
        if not any_zero and _ln_qc(moles_zero, stoich, volume) == ln_kc:
            xi_star = Decimal(0)
            residual = Decimal(0)
            moles = moles_zero
        else:
            xi_star, residual = _bisect(
                xi_min, xi_max, formulas, initial, stoich, volume, ln_kc
            )
            moles = _moles_at(xi_star, formulas, initial, stoich)

        # ---- 反应方向：与输出的 extent 保持一致——|ξ| 不足半个输出单位
        # （即 extent 显示为 0.0000）时记为已平衡 ----
        equilibrium_tol = Decimal(1).scaleb(-precision) / 2
        if xi_star > equilibrium_tol:
            direction = "forward"
        elif xi_star < -equilibrium_tol:
            direction = "reverse"
        else:
            direction = "at_equilibrium"

        def qc_from(moles_: list[Decimal]) -> Decimal:
            product = Decimal(1)
            for n, nu in zip(moles_, stoich):
                product *= (n / volume) ** nu
            return product

        qc = qc_from(moles)

        species_rows = []
        for i, formula in enumerate(formulas):
            n_eq = moles[i]
            # 防御性截断：数学上严格非负，避免舍入噪声给出极小负数
            if n_eq < 0:
                n_eq = Decimal(0)
            species_rows.append(
                {
                    "formula": formula,
                    "side": "reactant" if i < n_reactants else "product",
                    "coefficient": abs(stoich[i]),
                    "initial_moles": _safe_round(initial[formula], precision),
                    "equilibrium_moles": _safe_round(n_eq, precision),
                    "concentration": _safe_round(n_eq / volume, precision),
                }
            )

        return {
            "equation": balanced["equation"],
            "coefficients": balanced["coefficients"],
            "reactant_coefficients": r_coeffs,
            "product_coefficients": p_coeffs,
            "precision": precision,
            "volume": _safe_round(volume, precision),
            "kc": _format_fixed(kc, precision),
            "direction": direction,
            "extent": _safe_round(xi_star, precision),
            "extent_min": _safe_round(xi_min, precision),
            "extent_max": _safe_round(xi_max, precision),
            "qc": _format_fixed(qc, precision),
            "residual": _format_residual(residual, precision),
            "species": species_rows,
        }

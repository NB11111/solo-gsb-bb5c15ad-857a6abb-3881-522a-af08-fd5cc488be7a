"""单反应平衡组成接口测试。"""
import pytest
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


def equilibrium(reactants, products, initial_moles, volume, kc, precision=None):
    payload = {
        "reactants": reactants,
        "products": products,
        "initial_moles": [
            {"formula": f, "moles": n} for f, n in initial_moles.items()
        ],
        "volume": str(volume),
        "kc": str(kc),
    }
    if precision is not None:
        payload["precision"] = precision
    resp = client.post("/api/v1/equilibrium", json=payload)
    assert resp.status_code == 200, resp.text
    return resp.json()


def equilibrium_err(reactants, products, initial_moles, volume, kc):
    resp = client.post(
        "/api/v1/equilibrium",
        json={
            "reactants": reactants,
            "products": products,
            "initial_moles": [
                {"formula": f, "moles": n} for f, n in initial_moles.items()
            ],
            "volume": str(volume),
            "kc": str(kc),
        },
    )
    assert resp.status_code in (400, 422), resp.text
    return resp


# ---------- 正常计算 ----------

def test_forward_zero_product_n2o4():
    # N2O4 ⇌ 2NO2，1 mol / 2 L，Kc=0.1：
    # (2ξ/2)²/((1-ξ)/2)=0.1 => 2ξ²+0.1ξ-0.1=0 => ξ=0.2
    data = equilibrium(
        ["N2O4"], ["NO2"], {"N2O4": 1, "NO2": 0}, volume=2, kc="0.1"
    )
    assert data["equation"] == "N2O4 -> 2NO2"
    assert data["coefficients"] == [1, 2]
    assert data["direction"] == "forward"
    assert data["extent"] == "0.2000"
    assert data["extent_min"] == "0.0000"
    assert data["extent_max"] == "1.0000"
    assert data["qc"] == "0.1000"
    assert data["residual"] == "0.0000"
    n2o4, no2 = data["species"]
    assert n2o4["side"] == "reactants"
    assert n2o4["initial_moles"] == "1.0000"
    assert n2o4["equilibrium_moles"] == "0.8000"
    assert n2o4["concentration"] == "0.4000"
    assert no2["side"] == "products"
    assert no2["coefficient"] == 2
    assert no2["initial_moles"] == "0.0000"
    assert no2["equilibrium_moles"] == "0.4000"
    assert no2["concentration"] == "0.2000"


def test_reverse_zero_reactants_hi():
    # H2 + I2 ⇌ 2HI，初始只有 2 mol HI，V=1，Kc=49：
    # (2-2ξ)²/ξ²=49 => ξ=-2/9（逆向），HI 余 14/9，H2=I2=2/9
    data = equilibrium(
        ["H2", "I2"], ["HI"], {"H2": 0, "I2": 0, "HI": 2}, volume=1, kc=49,
        precision=6,
    )
    assert data["direction"] == "reverse"
    assert data["extent"] == "-0.222222"
    assert data["extent_min"] == "-1.000000"
    assert data["extent_max"] == "0.000000"
    assert data["qc"] == "49.000000"
    assert data["residual"] == "0.000000"
    h2, i2, hi = data["species"]
    assert h2["equilibrium_moles"] == "0.222222"
    assert i2["equilibrium_moles"] == "0.222222"
    assert hi["equilibrium_moles"] == "1.555556"


def test_forward_all_positive_initial():
    # H2=I2=HI=1，V=1，Kc=49：(1+2ξ)²/((1-ξ)²)=49 => ξ=2/3
    data = equilibrium(
        ["H2", "I2"], ["HI"], {"H2": 1, "I2": 1, "HI": 1}, volume=1, kc=49,
        precision=6,
    )
    assert data["direction"] == "forward"
    assert data["extent"] == "0.666667"
    assert data["qc"] == "49.000000"
    assert data["residual"] == "0.000000"
    h2, i2, hi = data["species"]
    assert h2["equilibrium_moles"] == "0.333333"
    assert hi["equilibrium_moles"] == "2.333333"


def test_reverse_all_positive_initial():
    # H2=I2=HI=1，V=1，Kc=0.01：(1+2ξ)²/(1-ξ)²=0.01 => ξ=-3/7
    data = equilibrium(
        ["H2", "I2"], ["HI"], {"H2": 1, "I2": 1, "HI": 1}, volume=1, kc="0.01",
        precision=6,
    )
    assert data["direction"] == "reverse"
    assert data["extent"] == "-0.428571"
    assert data["qc"] == "0.010000"
    assert data["residual"] == "0.000000"


def test_initial_exactly_at_equilibrium():
    # H2=I2=0.1，HI=0.7，V=1：Qc=0.7²/(0.1·0.1)=49 恰为 Kc
    data = equilibrium(
        ["H2", "I2"], ["HI"],
        {"H2": "0.1", "I2": "0.1", "HI": "0.7"}, volume=1, kc=49,
    )
    assert data["direction"] == "equilibrium"
    assert data["extent"] == "0.0000"
    assert data["qc"] == "49.0000"
    assert data["residual"] == "0.0000"


def test_volume_affects_equilibrium_delta_nu_nonzero():
    # N2O4 ⇌ 2NO2，1 mol / 0.5 L，Kc=0.045：
    # 4ξ²/(0.5(1-ξ))=0.045，用返回组成反算 Qc 验证
    data = equilibrium(
        ["N2O4"], ["NO2"], {"N2O4": 1, "NO2": 0}, volume="0.5", kc="0.045",
        precision=8,
    )
    from decimal import Decimal

    xi = Decimal(data["extent"])
    V = Decimal("0.5")
    qc = (2 * xi / V) ** 2 / ((1 - xi) / V)
    assert abs(qc - Decimal("0.045")) < Decimal("0.0000005")


def test_nonunit_stoichiometry_nh3():
    # N2 + 3H2 ⇌ 2NH3，初始 1/3/0 mol，V=1，Kc=1
    data = equilibrium(
        ["N2", "H2"], ["NH3"], {"N2": 1, "H2": 3, "NH3": 0},
        volume=1, kc=1, precision=8,
    )
    assert data["coefficients"] == [1, 3, 2]
    assert data["direction"] == "forward"
    from decimal import Decimal

    rows = {s["formula"]: Decimal(s["concentration"]) for s in data["species"]}
    qc = rows["NH3"] ** 2 / (rows["N2"] * rows["H2"] ** 3)
    assert abs(qc - 1) < Decimal("0.0000005")


def test_concentration_uses_volume():
    # 同一组摩尔数、不同体积 => 浓度随 V 变，平衡进度也随 Δν≠0 变化
    data_2 = equilibrium(
        ["N2O4"], ["NO2"], {"N2O4": 1, "NO2": 0}, volume=2, kc="0.1",
        precision=6,
    )
    data_4 = equilibrium(
        ["N2O4"], ["NO2"], {"N2O4": 1, "NO2": 0}, volume=4, kc="0.1",
        precision=6,
    )
    assert data_2["extent"] == "0.200000"
    # V=4：(2ξ/4)²/((1-ξ)/4)=0.1 => ξ²+0.1ξ-0.1=0 => ξ≈0.270156
    assert data_4["extent"] == "0.270156"
    assert data_4["species"][1]["concentration"] == "0.135078"


def test_extreme_kc_near_completion():
    # Kc 极小、只有生成物：极限靠近 ξ=-0.5（HI 几乎全分解），不应报错
    data = equilibrium(
        ["H2", "I2"], ["HI"], {"H2": 0, "I2": 0, "HI": 1},
        volume=1, kc="1e-30", precision=6,
    )
    assert data["direction"] == "reverse"
    assert data["extent"] == "-0.500000"
    assert data["residual"] == "0.000000"


def test_precision_parameter_echoes():
    data = equilibrium(
        ["H2", "I2"], ["HI"], {"H2": 0, "I2": 0, "HI": 2}, volume=1, kc=49,
        precision=2,
    )
    assert data["precision"] == 2
    assert data["extent"] == "-0.22"
    data = equilibrium(
        ["H2", "I2"], ["HI"], {"H2": 0, "I2": 0, "HI": 2}, volume=1, kc=49,
        precision=8,
    )
    assert data["precision"] == 8
    assert data["extent"] == "-0.22222222"


# ---------- 4xx 错误 ----------

def test_singleton_feasible_interval_returns_422():
    resp = equilibrium_err(
        ["H2", "I2"], ["HI"], {"H2": 1, "I2": 0, "HI": 0}, volume=1, kc=49
    )
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert "可行" in detail["message"]
    assert detail["zero_reactants"] == ["I2"]
    assert detail["zero_products"] == ["HI"]
    assert detail["extent_min"] == detail["extent_max"] == "0.0000"


def test_all_zero_initial_moles_returns_400():
    resp = equilibrium_err(
        ["H2", "I2"], ["HI"], {"H2": 0, "I2": 0, "HI": 0}, volume=1, kc=49
    )
    assert resp.status_code == 400
    assert "至少" in resp.json()["detail"]


@pytest.mark.parametrize("bad", [-1, "-0.01"])
def test_negative_initial_moles_returns_400(bad):
    resp = equilibrium_err(
        ["H2", "I2"], ["HI"],
        {"H2": bad, "I2": 1, "HI": 0}, volume=1, kc=49,
    )
    assert resp.status_code == 400
    assert "非负" in resp.json()["detail"]


@pytest.mark.parametrize("bad_volume", [0, -1, "-2.5"])
def test_bad_volume_returns_400(bad_volume):
    resp = equilibrium_err(
        ["H2", "I2"], ["HI"], {"H2": 1, "I2": 1, "HI": 0},
        volume=bad_volume, kc=49,
    )
    assert resp.status_code == 400
    assert "体积" in resp.json()["detail"]


@pytest.mark.parametrize("bad_kc", [0, -10])
def test_bad_kc_returns_400(bad_kc):
    resp = equilibrium_err(
        ["H2", "I2"], ["HI"], {"H2": 1, "I2": 1, "HI": 0},
        volume=1, kc=bad_kc,
    )
    assert resp.status_code == 400
    assert "Kc" in resp.json()["detail"]


def test_missing_initial_moles_returns_400():
    resp = equilibrium_err(
        ["H2", "I2"], ["HI"], {"H2": 1, "I2": 1}, volume=1, kc=49
    )
    assert resp.status_code == 400
    assert "缺少" in resp.json()["detail"]
    assert "HI" in resp.json()["detail"]


def test_duplicate_initial_moles_returns_400():
    resp = client.post(
        "/api/v1/equilibrium",
        json={
            "reactants": ["H2", "I2"],
            "products": ["HI"],
            "initial_moles": [
                {"formula": "H2", "moles": 1},
                {"formula": "H2", "moles": 2},
                {"formula": "I2", "moles": 1},
                {"formula": "HI", "moles": 0},
            ],
            "volume": 1,
            "kc": 49,
        },
    )
    assert resp.status_code == 400
    assert "重复" in resp.json()["detail"]


def test_unknown_species_initial_moles_returns_400():
    resp = equilibrium_err(
        ["H2", "I2"], ["HI"], {"H2": 1, "I2": 1, "O2": 1}, volume=1, kc=49
    )
    assert resp.status_code == 400
    assert "不在" in resp.json()["detail"]


def test_bad_formula_returns_400():
    resp = equilibrium_err(
        ["H2O$"], ["H2"], {"H2O$": 1, "H2": 0}, volume=1, kc=1
    )
    assert resp.status_code == 400


def test_unbalanceable_returns_422():
    resp = equilibrium_err(
        ["H2", "O2"], ["H2O", "H2O2"],
        {"H2": 1, "O2": 1, "H2O": 0, "H2O2": 0}, volume=1, kc=1,
    )
    assert resp.status_code == 422
    assert resp.json()["detail"]["null_space_dim"] == 2


@pytest.mark.parametrize("precision", [1, 9])
def test_precision_out_of_range_returns_422(precision):
    resp = client.post(
        "/api/v1/equilibrium",
        json={
            "reactants": ["H2", "I2"],
            "products": ["HI"],
            "initial_moles": [
                {"formula": "H2", "moles": 1},
                {"formula": "I2", "moles": 1},
                {"formula": "HI", "moles": 0},
            ],
            "volume": 1,
            "kc": 49,
            "precision": precision,
        },
    )
    assert resp.status_code == 422


def test_min_length_validation_returns_422():
    resp = client.post(
        "/api/v1/equilibrium",
        json={
            "reactants": [],
            "products": ["HI"],
            "initial_moles": [{"formula": "HI", "moles": 1}],
            "volume": 1,
            "kc": 49,
        },
    )
    assert resp.status_code == 422

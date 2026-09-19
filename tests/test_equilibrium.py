"""单反应平衡组成接口测试。

解析/配平均复用现有模块，这里使用可解析的真实分子式：
  - N2O4 ⇌ 2NO2（Δn≠0，含初始为零、正/逆向）
  - NH4OCN ⇌ CO(NH2)2（1:1 异构化，Δn=0，体积无关）
  - N2 ⇌ 2N（Δn≠0，体积效应）
  - 2H2 + O2 ⇌ 2H2O（多物种、可行区间退化）
"""
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


def equilibrium(reactants, products, initial_moles, volume, kc, precision=None):
    payload = {
        "reactants": reactants,
        "products": products,
        "initial_moles": [
            {"formula": f, "moles": m} for f, m in initial_moles
        ],
        "volume": volume,
        "kc": kc,
    }
    if precision is not None:
        payload["precision"] = precision
    resp = client.post("/api/v1/equilibrium", json=payload)
    assert resp.status_code == 200, resp.text
    return resp.json()


# ---------- 正常计算 ----------

def test_forward_from_zero_product():
    # N2O4 ⇌ 2NO2，V=1 L，n0(N2O4)=1，Kc=0.04
    # 4ξ²/(1-ξ)=0.04 => ξ=(-0.04+√0.6416)/8≈0.0951249
    data = equilibrium(
        ["N2O4"], ["NO2"],
        [("N2O4", "1"), ("NO2", "0")],
        "1", "0.04",
    )
    assert data["equation"] == "N2O4 -> 2NO2"
    assert data["coefficients"] == [1, 2]
    assert data["reactant_coefficients"] == [1]
    assert data["product_coefficients"] == [2]
    assert data["precision"] == 4
    assert data["volume"] == "1.0000"
    assert data["kc"] == "0.0400"
    assert data["direction"] == "forward"
    assert data["extent"] == "0.0951"
    assert data["extent_min"] == "0.0000"
    assert data["extent_max"] == "1.0000"
    assert data["qc"] == "0.0400"
    assert data["residual"] == "+0.0000"

    n2o4, no2 = data["species"]
    assert n2o4["formula"] == "N2O4"
    assert n2o4["side"] == "reactant"
    assert n2o4["coefficient"] == 1
    assert n2o4["initial_moles"] == "1.0000"
    assert n2o4["equilibrium_moles"] == "0.9049"
    assert n2o4["concentration"] == "0.9049"
    assert no2["formula"] == "NO2"
    assert no2["side"] == "product"
    assert no2["coefficient"] == 2
    assert no2["initial_moles"] == "0.0000"
    assert no2["equilibrium_moles"] == "0.1902"
    assert no2["concentration"] == "0.1902"


def test_reverse_from_zero_reactant():
    # 初始只有生成物 NO2=1 mol，反应须逆向进行：
    # (1-2y)²/y=0.04 => y≈0.434113，ξ=-y，可行下界 -n0(NO2)/2=-0.5
    data = equilibrium(
        ["N2O4"], ["NO2"],
        [("N2O4", "0"), ("NO2", "1")],
        "1", "0.04", precision=6,
    )
    assert data["direction"] == "reverse"
    assert data["extent"] == "-0.434113"
    assert data["extent_min"] == "-0.500000"
    assert data["extent_max"] == "0.000000"
    assert data["qc"] == "0.040000"

    n2o4, no2 = data["species"]
    assert n2o4["equilibrium_moles"] == "0.434113"
    assert n2o4["concentration"] == "0.434113"
    assert no2["equilibrium_moles"] == "0.131774"


def test_reverse_with_both_sides_present():
    # 1:1 异构化 NH4OCN ⇌ CO(NH2)2，Kc=1；n0=1/3 => ξ=-1
    data = equilibrium(
        ["NH4OCN"], ["CO(NH2)2"],
        [("NH4OCN", "1"), ("CO(NH2)2", "3")],
        "1", "1",
    )
    assert data["direction"] == "reverse"
    assert data["extent"] == "-1.0000"
    assert data["extent_min"] == "-3.0000"
    assert data["extent_max"] == "1.0000"
    assert data["qc"] == "1.0000"
    react, prod = data["species"]
    assert react["equilibrium_moles"] == "2.0000"
    assert prod["equilibrium_moles"] == "2.0000"


def test_initially_at_equilibrium():
    # 等摩尔、Kc=1 => Qc=Kc，ξ=0，方向记为 at_equilibrium
    data = equilibrium(
        ["NH4OCN"], ["CO(NH2)2"],
        [("NH4OCN", "2"), ("CO(NH2)2", "2")],
        "1", "1",
    )
    assert data["direction"] == "at_equilibrium"
    assert data["extent"] == "0.0000"
    assert data["qc"] == "1.0000"
    assert data["residual"] == "+0.0000"
    assert all(s["equilibrium_moles"] == s["initial_moles"] for s in data["species"])


def test_delta_n_zero_volume_independent():
    # 1:1 反应与体积无关：V=10 与 V=1 都有 ξ=0.5
    for volume in ("1", "10"):
        data = equilibrium(
            ["NH4OCN"], ["CO(NH2)2"],
            [("NH4OCN", "1"), ("CO(NH2)2", "0")],
            volume, "1",
        )
        assert data["extent"] == "0.5000"
        assert data["qc"] == "1.0000"
        react, prod = data["species"]
        assert react["equilibrium_moles"] == "0.5000"
        assert prod["equilibrium_moles"] == "0.5000"
        assert react["concentration"] == prod["concentration"]


def test_volume_effect_delta_n_nonzero():
    # N2 ⇌ 2N，Kc=1：(2ξ/V)²/((1-ξ)/V)=1
    # V=4 => 4ξ²/(4(1-ξ))=... => ξ²/(1-ξ)=1 => ξ=(√5-1)/2≈0.618034
    data = equilibrium(
        ["N2"], ["N"],
        [("N2", "1"), ("N", "0")],
        "4", "1", precision=6,
    )
    assert data["extent"] == "0.618034"
    n2, n = data["species"]
    assert n2["equilibrium_moles"] == "0.381966"
    assert n2["concentration"] == "0.095492"
    assert n["equilibrium_moles"] == "1.236068"
    assert n["concentration"] == "0.309017"

    # V=1：4ξ²/(1-ξ)=1 => ξ=(√17-1)/8≈0.390388
    data = equilibrium(
        ["N2"], ["N"],
        [("N2", "1"), ("N", "0")],
        "1", "1", precision=6,
    )
    assert data["extent"] == "0.390388"


def test_multi_species_forward():
    # 2H2 + O2 ⇌ 2H2O，V=10 L，Kc=1000，初始 2/1/0
    data = equilibrium(
        ["H2", "O2"], ["H2O"],
        [("H2", "2"), ("O2", "1"), ("H2O", "0")],
        "10", "1000", precision=6,
    )
    assert data["coefficients"] == [2, 1, 2]
    assert data["direction"] == "forward"
    assert data["extent"] == "0.812420"
    assert data["extent_min"] == "0.000000"
    assert data["extent_max"] == "1.000000"
    h2, o2, h2o = data["species"]
    assert h2["coefficient"] == 2 and h2["side"] == "reactant"
    assert o2["coefficient"] == 1 and o2["side"] == "reactant"
    assert h2o["coefficient"] == 2 and h2o["side"] == "product"
    assert h2["equilibrium_moles"] == "0.375160"
    assert h2["concentration"] == "0.037516"
    assert o2["equilibrium_moles"] == "0.187580"
    assert o2["concentration"] == "0.018758"
    assert h2o["equilibrium_moles"] == "1.624840"
    assert h2o["concentration"] == "0.162484"
    assert data["qc"] == "1000.000000"


def test_precision_parameter():
    data = equilibrium(
        ["N2O4"], ["NO2"],
        [("N2O4", "1"), ("NO2", "0")],
        "1", "0.04", precision=8,
    )
    assert data["precision"] == 8
    assert data["extent"] == "0.09512492"
    assert data["qc"] == "0.04000000"


def test_multi_species_reverse():
    # 2H2 + O2 ⇌ 2H2O，初始只有 H2O=2 mol，V=1，Kc=2：
    # ξ=-0.5 时 H2O=1、H2=1、O2=0.5，Qc=1²/(1²·0.5)=2，恰为平衡
    data = equilibrium(
        ["H2", "O2"], ["H2O"],
        [("H2", "0"), ("O2", "0"), ("H2O", "2")],
        "1", "2",
    )
    assert data["direction"] == "reverse"
    assert data["extent"] == "-0.5000"
    assert data["extent_min"] == "-1.0000"
    assert data["extent_max"] == "0.0000"
    assert data["qc"] == "2.0000"
    h2, o2, h2o = data["species"]
    assert h2["equilibrium_moles"] == "1.0000"
    assert o2["equilibrium_moles"] == "0.5000"
    assert h2o["equilibrium_moles"] == "1.0000"


def test_extreme_kc_uses_scientific_notation():
    # Kc 极小，根极度靠近 ξ=0；Qc 定点形式会丢失有效数字，改科学计数法
    data = equilibrium(
        ["N2O4"], ["NO2"],
        [("N2O4", "1"), ("NO2", "0")],
        "1", "1E-50",
    )
    assert data["qc"] == "1.0000E-50"
    assert data["residual"] == "+0.0000"


# ---------- 4xx 错误 ----------

def post_raw(payload):
    return client.post("/api/v1/equilibrium", json=payload)


def base_payload(**overrides):
    payload = {
        "reactants": ["N2O4"],
        "products": ["NO2"],
        "initial_moles": [
            {"formula": "N2O4", "moles": "1"},
            {"formula": "NO2", "moles": "0"},
        ],
        "volume": "1",
        "kc": "0.04",
    }
    payload.update(overrides)
    return payload


def test_missing_initial_moles_returns_400():
    resp = post_raw(base_payload(
        initial_moles=[{"formula": "N2O4", "moles": "1"}],
    ))
    assert resp.status_code == 400
    assert "缺少" in resp.json()["detail"]
    assert "NO2" in resp.json()["detail"]


def test_duplicate_initial_moles_returns_400():
    resp = post_raw(base_payload(
        initial_moles=[
            {"formula": "N2O4", "moles": "1"},
            {"formula": "NO2", "moles": "0"},
            {"formula": "NO2", "moles": "0.5"},
        ],
    ))
    assert resp.status_code == 400
    assert "重复" in resp.json()["detail"]


def test_unknown_species_initial_moles_returns_400():
    resp = post_raw(base_payload(
        initial_moles=[
            {"formula": "N2O4", "moles": "1"},
            {"formula": "NO2", "moles": "0"},
            {"formula": "H2O", "moles": "1"},
        ],
    ))
    assert resp.status_code == 400
    assert "不在" in resp.json()["detail"]


def test_duplicate_species_in_equation_returns_400():
    resp = post_raw({
        "reactants": ["H2O"],
        "products": ["H2O"],
        "initial_moles": [{"formula": "H2O", "moles": "1"}],
        "volume": "1",
        "kc": "1",
    })
    assert resp.status_code == 400
    assert "重复" in resp.json()["detail"]


def test_negative_initial_moles_returns_400():
    resp = post_raw(base_payload(
        initial_moles=[
            {"formula": "N2O4", "moles": "-0.5"},
            {"formula": "NO2", "moles": "0"},
        ],
    ))
    assert resp.status_code == 400
    assert "非负" in resp.json()["detail"]


def test_zero_initial_moles_allowed():
    # 非负包含零：生成物初始为零是最常见场景，上面已覆盖；
    # 这里确认反应物初始为零（纯逆向）同样合法
    resp = post_raw(base_payload(
        initial_moles=[
            {"formula": "N2O4", "moles": "0"},
            {"formula": "NO2", "moles": "1"},
        ],
    ))
    assert resp.status_code == 200, resp.text
    assert resp.json()["direction"] == "reverse"


def test_invalid_volume_returns_400():
    for volume in ("0", "-1.5"):
        resp = post_raw(base_payload(volume=volume))
        assert resp.status_code == 400
        assert "体积" in resp.json()["detail"]


def test_invalid_kc_returns_400():
    for kc in ("0", "-1"):
        resp = post_raw(base_payload(kc=kc))
        assert resp.status_code == 400
        assert "Kc" in resp.json()["detail"]


def test_degenerate_interval_returns_422():
    # O2 与 H2O 初始均为零：可行区间退化为 ξ=0，Qc 为 0/0
    resp = post_raw({
        "reactants": ["H2", "O2"],
        "products": ["H2O"],
        "initial_moles": [
            {"formula": "H2", "moles": "2"},
            {"formula": "O2", "moles": "0"},
            {"formula": "H2O", "moles": "0"},
        ],
        "volume": "1",
        "kc": "1",
    })
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert "message" in detail
    assert detail["xi_min"] == "0.0000"
    assert detail["xi_max"] == "0.0000"
    assert "O2" in detail["zero_reactants"]
    assert "H2O" in detail["zero_products"]


def test_unbalanceable_returns_422():
    resp = post_raw({
        "reactants": ["H2", "O2"],
        "products": ["H2O", "H2O2"],
        "initial_moles": [
            {"formula": f, "moles": "1"}
            for f in ("H2", "O2", "H2O", "H2O2")
        ],
        "volume": "1",
        "kc": "1",
    })
    assert resp.status_code == 422
    assert resp.json()["detail"]["null_space_dim"] == 2


def test_bad_formula_returns_400():
    resp = post_raw({
        "reactants": ["N2O4$"],
        "products": ["NO2"],
        "initial_moles": [
            {"formula": "N2O4$", "moles": "1"},
            {"formula": "NO2", "moles": "0"},
        ],
        "volume": "1",
        "kc": "0.04",
    })
    assert resp.status_code == 400


def test_precision_out_of_range_returns_422():
    for precision in (1, 9):
        resp = post_raw(base_payload(precision=precision))
        assert resp.status_code == 422

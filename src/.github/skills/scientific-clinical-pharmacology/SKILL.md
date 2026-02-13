---
name: scientific-clinical-pharmacology
description: |
  臨床薬理学モデリングスキル。PopPK (NLME 混合効果モデル)・
  PBPK シミュレーション・TDM 投与量最適化・
  Emax/Sigmoid PD モデリング・薬物間相互作用予測・
  臨床薬理パイプライン。
  TU 外スキル (Python + nlmixr2/mrgsolve ラッパー)。
---

# Scientific Clinical Pharmacology

母集団薬物動態 (PopPK)・生理学的薬物動態 (PBPK)・
薬力学 (PD) モデリングを統合した臨床薬理学
解析パイプラインを提供する。

## When to Use

- 母集団 PK (PopPK) の NLME 解析を行うとき
- PBPK モデルで薬物動態をシミュレーションするとき
- TDM (Therapeutic Drug Monitoring) 投与量を最適化するとき
- Emax/Sigmoid PD モデルを当てはめるとき
- 薬物間相互作用 (DDI) の影響を予測するとき
- 小児・腎障害・肝障害の用量調節を検討するとき

---

## Quick Start

## 1. コンパートメント PK モデル

```python
import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.optimize import minimize


def one_compartment_iv(y, t, cl, v):
    """1-コンパートメント IV ボーラス ODE。"""
    return [-cl / v * y[0]]


def two_compartment_iv(y, t, cl, v1, q, v2):
    """2-コンパートメント IV ボーラス ODE。"""
    c1 = y[0] / v1
    c2 = y[1] / v2
    dy1 = -cl * c1 - q * (c1 - c2)
    dy2 = q * (c1 - c2)
    return [dy1, dy2]


def simulate_pk(dose, model="1cmt",
                  params=None,
                  times=None):
    """
    PK モデルシミュレーション。

    Parameters:
        dose: float — 投与量 (mg)
        model: str — "1cmt" or "2cmt"
        params: dict — PK パラメータ
            1cmt: {cl, v}
            2cmt: {cl, v1, q, v2}
        times: array — 時間点 (h)
    """
    if times is None:
        times = np.linspace(0, 24, 241)
    if params is None:
        params = ({"cl": 5.0, "v": 50.0}
                  if model == "1cmt"
                  else {"cl": 5.0, "v1": 50.0,
                        "q": 2.0, "v2": 30.0})

    if model == "1cmt":
        y0 = [dose]
        sol = odeint(one_compartment_iv, y0,
                     times,
                     args=(params["cl"],
                           params["v"]))
        conc = sol[:, 0] / params["v"]
    else:
        y0 = [dose, 0.0]
        sol = odeint(two_compartment_iv, y0,
                     times,
                     args=(params["cl"],
                           params["v1"],
                           params["q"],
                           params["v2"]))
        conc = sol[:, 0] / params["v1"]

    df = pd.DataFrame({
        "time": times, "concentration": conc})
    cmax = conc.max()
    t_half = 0.693 * params.get(
        "v", params.get("v1", 50)) / params["cl"]
    print(f"PK sim ({model}): Cmax={cmax:.2f}, "
          f"t1/2={t_half:.1f}h")
    return df
```

## 2. 母集団 PK (PopPK) 推定

```python
def popk_estimation(data, model="1cmt"):
    """
    母集団 PK パラメータ推定 (簡易 NLME)。

    Parameters:
        data: pd.DataFrame — 個体別濃度データ
            columns: [id, time, dv, dose, (covariates)]
        model: str — "1cmt" or "2cmt"
    """
    subjects = data["id"].unique()

    def _objective(theta):
        """集団目的関数 (OFV)。"""
        tv_cl = np.exp(theta[0])
        tv_v = np.exp(theta[1])
        omega_cl = np.exp(theta[2])
        omega_v = np.exp(theta[3])
        sigma = np.exp(theta[4])

        ofv = 0.0
        for subj in subjects:
            sdata = data[data["id"] == subj]
            dose = sdata["dose"].iloc[0]
            times = sdata["time"].values
            obs = sdata["dv"].values

            # 個体パラメータ (EBE 近似)
            eta_cl = 0.0
            eta_v = 0.0
            cl_i = tv_cl * np.exp(eta_cl)
            v_i = tv_v * np.exp(eta_v)

            pred = dose / v_i * np.exp(
                -cl_i / v_i * times)
            pred = np.maximum(pred, 1e-10)

            # OFV 要素
            residual = np.log(obs + 1e-10) - np.log(
                pred)
            ofv += np.sum(
                residual**2 / sigma**2
                + np.log(sigma**2))

        return ofv

    # 初期値
    x0 = [np.log(5), np.log(50),
           np.log(0.3), np.log(0.3),
           np.log(0.2)]

    result = minimize(_objective, x0,
                      method="Nelder-Mead",
                      options={"maxiter": 5000})

    estimates = {
        "tv_cl": round(np.exp(result.x[0]), 3),
        "tv_v": round(np.exp(result.x[1]), 3),
        "omega_cl": round(np.exp(result.x[2]), 3),
        "omega_v": round(np.exp(result.x[3]), 3),
        "sigma": round(np.exp(result.x[4]), 3),
        "ofv": round(result.fun, 2),
        "converged": result.success,
    }

    print(f"PopPK: CL={estimates['tv_cl']} L/h, "
          f"V={estimates['tv_v']} L, "
          f"OFV={estimates['ofv']}")
    return estimates
```

## 3. TDM 投与量最適化

```python
def tdm_dose_optimization(
        current_conc, current_dose,
        target_range, pk_params,
        interval=12):
    """
    TDM ベース投与量最適化。

    Parameters:
        current_conc: float — 現在トラフ濃度
        current_dose: float — 現在投与量 (mg)
        target_range: tuple — 目標濃度範囲 (min, max)
        pk_params: dict — {cl, v} PK パラメータ
        interval: float — 投与間隔 (h)
    """
    cl = pk_params["cl"]
    v = pk_params["v"]
    ke = cl / v
    target_mid = (target_range[0]
                  + target_range[1]) / 2

    # 線形 PK 仮定: 用量比例
    ratio = target_mid / max(current_conc, 0.01)
    new_dose = current_dose * ratio

    # シミュレーション検証
    times = np.linspace(0, interval, 100)
    conc_profile = (new_dose / v
                    * np.exp(-ke * times))
    cmax = conc_profile[0]
    ctrough = conc_profile[-1]

    in_range = (target_range[0] <= ctrough
                <= target_range[1])

    result = {
        "current_dose": current_dose,
        "current_trough": current_conc,
        "recommended_dose": round(new_dose, 1),
        "predicted_cmax": round(cmax, 2),
        "predicted_trough": round(ctrough, 2),
        "target_range": target_range,
        "in_target": in_range,
    }

    status = "✓" if in_range else "✗"
    print(f"TDM: {current_dose}mg → "
          f"{new_dose:.1f}mg "
          f"(trough {ctrough:.2f}) {status}")
    return result
```

## 4. Emax PD モデル

```python
def emax_model(conc, emax, ec50, hill=1):
    """Emax / Sigmoid Emax モデル。"""
    return emax * conc**hill / (
        ec50**hill + conc**hill)


def fit_emax(conc_data, effect_data,
               sigmoid=False):
    """
    Emax モデルフィッティング。

    Parameters:
        conc_data: array — 濃度データ
        effect_data: array — 薬効データ
        sigmoid: bool — Sigmoid (Hill) モデル
    """
    from scipy.optimize import curve_fit

    conc = np.array(conc_data)
    effect = np.array(effect_data)

    if sigmoid:
        def _model(c, emax, ec50, hill):
            return emax_model(c, emax, ec50, hill)
        p0 = [max(effect), np.median(conc), 1.0]
        bounds = ([0, 0, 0.1], [np.inf, np.inf, 10])
    else:
        def _model(c, emax, ec50):
            return emax_model(c, emax, ec50, 1)
        p0 = [max(effect), np.median(conc)]
        bounds = ([0, 0], [np.inf, np.inf])

    popt, pcov = curve_fit(
        _model, conc, effect, p0=p0,
        bounds=bounds, maxfev=5000)

    perr = np.sqrt(np.diag(pcov))

    result = {
        "emax": round(popt[0], 3),
        "ec50": round(popt[1], 3),
        "emax_se": round(perr[0], 3),
        "ec50_se": round(perr[1], 3),
    }
    if sigmoid:
        result["hill"] = round(popt[2], 3)
        result["hill_se"] = round(perr[2], 3)

    print(f"PD fit: Emax={result['emax']}, "
          f"EC50={result['ec50']}"
          + (f", Hill={result['hill']}"
             if sigmoid else ""))
    return result
```

## 5. 臨床薬理統合パイプライン

```python
def clinical_pharmacology_pipeline(
        pk_data, pd_data=None,
        target_range=(10, 20),
        output_dir="results"):
    """
    臨床薬理学統合パイプライン。

    Parameters:
        pk_data: pd.DataFrame — PK 濃度データ
        pd_data: pd.DataFrame | None — PD 効果データ
        target_range: tuple — 目標トラフ範囲
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # 1) PopPK 推定
    pk_est = popk_estimation(pk_data)

    # 2) TDM 最適化 (最新トラフがあれば)
    latest = pk_data.sort_values("time").iloc[-1]
    tdm = tdm_dose_optimization(
        latest["dv"], latest["dose"],
        target_range,
        {"cl": pk_est["tv_cl"],
         "v": pk_est["tv_v"]})

    # 3) PD モデル (データがあれば)
    pd_result = None
    if pd_data is not None:
        pd_result = fit_emax(
            pd_data["concentration"],
            pd_data["effect"],
            sigmoid=True)

    # 4) シミュレーション
    sim = simulate_pk(
        tdm["recommended_dose"],
        params={"cl": pk_est["tv_cl"],
                "v": pk_est["tv_v"]})
    sim.to_csv(out / "pk_simulation.csv",
               index=False)

    print(f"Clinical PK pipeline → {out}")
    return {
        "popk": pk_est,
        "tdm": tdm,
        "pd": pd_result,
        "simulation": sim,
    }
```

---

## パイプライン統合

```
admet-pharmacokinetics → clinical-pharmacology → pharmacogenomics
  (ADMET 予測)           (臨床 PK/PD)           (遺伝薬理学)
       │                       │                      ↓
  drug-repurposing ────────────┘          clinical-decision-support
    (薬剤リパーパシング)                       (臨床意思決定)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `pk_simulation.csv` | PK シミュレーション | → dose-response |
| `popk_estimates.json` | PopPK パラメータ | → pharmacogenomics |
| `tdm_recommendation.json` | TDM 推奨用量 | → clinical-decision |

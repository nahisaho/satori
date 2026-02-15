---
name: scientific-materials-characterization
description: |
  薄膜・材料キャラクタリゼーション解析のスキル。Thornton-Anders 構造ゾーンモデル（SZM）、
  XRD 結晶子サイズ解析（Scherrer 方程式）、Williamson-Hall プロット、多技法融合データ解析、
  PSP フレームワーク設計を行う際に使用。Scientific Skills Exp-13 で確立したパターン。
tu_tools:
  - key: materials_project
    name: Materials Project
    description: 材料特性データベース検索
---

# Scientific Materials Characterization

薄膜工学・材料科学における構造キャラクタリゼーション解析パイプラインスキル。
XRD、AFM、四探針法、UV-Vis などの測定データを統合し、Process-Structure-Property
（PSP）の因果連鎖を定量化する。

## When to Use

- 薄膜の成膜条件と膜特性の相関を解析したいとき
- XRD データから結晶子サイズ・格子歪を解析したいとき
- Thornton-Anders 構造ゾーンモデルにデータをマッピングしたいとき
- 多技法（XRD + AFM + 電気 + 光学）データを統合した解析
- PSP フレームワーク（Process → Structure → Property）

## Quick Start

## 1. 材料定義テンプレート

```python
# 材料物性定数辞書
MATERIAL_PROPERTIES = {
    "ZnO":  {"Tm_K": 2248, "Eg_eV": 3.37, "type": "TCO",  "crystal": "wurtzite"},
    "ITO":  {"Tm_K": 2200, "Eg_eV": 3.70, "type": "TCO",  "crystal": "cubic"},
    "Al2O3": {"Tm_K": 2345, "Eg_eV": 8.80, "type": "dielectric", "crystal": "corundum"},
    "HfO2": {"Tm_K": 3031, "Eg_eV": 5.80, "type": "dielectric", "crystal": "monoclinic"},
    "TiO2": {"Tm_K": 2116, "Eg_eV": 3.20, "type": "functional", "crystal": "rutile"},
    "SiO2": {"Tm_K": 1986, "Eg_eV": 8.90, "type": "dielectric", "crystal": "amorphous"},
}

# XRD 主要ピーク位置 (2θ, Cu-Kα)
XRD_PEAKS = {
    "ZnO":  {"(100)": 31.8, "(002)": 34.4, "(101)": 36.3},
    "ITO":  {"(222)": 30.6, "(400)": 35.5, "(440)": 51.0},
    "TiO2": {"(110)": 27.4, "(101)": 25.3, "(200)": 36.1},
}
```

## 2. 相同温度と構造ゾーンモデル（Thornton-Anders）

```python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def calculate_homologous_temperature(T_substrate_C, T_melt_K):
    """
    相同温度 T/Tm を算出する。
    T_substrate_C: 基板温度 (°C)
    T_melt_K: 融点 (K)
    """
    return (T_substrate_C + 273.15) / T_melt_K


def structure_zone_model(df, T_sub_col, material_col, materials_dict,
                          crystallite_col=None, figsize=(12, 8)):
    """
    Thornton-Anders 構造ゾーンモデルにデータをマッピングする。

    ゾーン境界:
    - Zone 1: T/Tm < 0.15 (繊維状 / アモルファス)
    - Zone T: 0.15 ≤ T/Tm < 0.30 (遷移帯 / 緻密微結晶)
    - Zone 2: 0.30 ≤ T/Tm < 0.50 (柱状結晶)
    - Zone 3: T/Tm ≥ 0.50 (等軸粒)
    """
    df = df.copy()

    # 相同温度を算出
    df["T_homologous"] = df.apply(
        lambda row: calculate_homologous_temperature(
            row[T_sub_col],
            materials_dict[row[material_col]]["Tm_K"]
        ), axis=1
    )

    # ゾーン分類
    def classify_zone(T_Tm):
        if T_Tm < 0.15:
            return "Zone 1"
        elif T_Tm < 0.30:
            return "Zone T"
        elif T_Tm < 0.50:
            return "Zone 2"
        else:
            return "Zone 3"

    df["Structure_Zone"] = df["T_homologous"].apply(classify_zone)

    # ──── 可視化 ────
    fig, ax = plt.subplots(figsize=figsize)

    zone_colors = {
        "Zone 1": "#FFE0B2",   # 薄橙
        "Zone T": "#C8E6C9",   # 薄緑
        "Zone 2": "#BBDEFB",   # 薄青
        "Zone 3": "#F8BBD0",   # 薄桃
    }

    # ゾーン背景
    boundaries = [0, 0.15, 0.30, 0.50, 0.80]
    zone_names = ["Zone 1\n(Fibrous)", "Zone T\n(Transition)",
                  "Zone 2\n(Columnar)", "Zone 3\n(Equiaxed)"]
    for i, (b_start, b_end) in enumerate(zip(boundaries[:-1], boundaries[1:])):
        rect = patches.Rectangle(
            (b_start, 0), b_end - b_start, 1,
            linewidth=0, facecolor=list(zone_colors.values())[i], alpha=0.3
        )
        ax.add_patch(rect)
        ax.text((b_start + b_end) / 2, 0.95, zone_names[i],
               ha="center", va="top", fontsize=9, fontweight="bold", alpha=0.6)

    # 材料別プロット
    unique_mats = df[material_col].unique()
    colors = plt.cm.Set2(np.linspace(0, 1, len(unique_mats)))

    for color, mat in zip(colors, unique_mats):
        mask = df[material_col] == mat
        subset = df[mask]

        # 正規化粒子エネルギー（出力 / 圧力の比のプロキシ）
        if "Power" in df.columns and "Working_Pressure" in df.columns:
            y_val = subset["Power"] / (subset["Working_Pressure"] * 100 + 1)
            y_val = (y_val - y_val.min()) / (y_val.max() - y_val.min() + 1e-10)
        else:
            y_val = np.random.uniform(0.3, 0.7, mask.sum())

        sizes = 50
        if crystallite_col and crystallite_col in df.columns:
            cs = subset[crystallite_col].values
            sizes = 30 + (cs - cs.min()) / (cs.max() - cs.min() + 1e-10) * 200

        ax.scatter(subset["T_homologous"], y_val,
                  s=sizes, c=[color], label=mat,
                  alpha=0.7, edgecolors="black", linewidth=0.5)

    ax.set_xlabel("Homologous Temperature T/Tₘ", fontsize=12)
    ax.set_ylabel("Normalized Particle Energy", fontsize=12)
    ax.set_title("Thornton-Anders Structure Zone Model", fontsize=14,
                fontweight="bold")
    ax.set_xlim(0, max(df["T_homologous"].max() * 1.1, 0.6))
    ax.set_ylim(0, 1)
    ax.legend(title="Material", bbox_to_anchor=(1.05, 1))

    # ゾーン境界線
    for b in [0.15, 0.30, 0.50]:
        ax.axvline(b, color="gray", linestyle="--", linewidth=1, alpha=0.5)

    plt.tight_layout()
    plt.savefig("figures/structure_zone_model.png", dpi=300, bbox_inches="tight")
    plt.close()

    # ゾーン統計
    zone_stats = df.groupby("Structure_Zone").agg(
        Count=("Structure_Zone", "size"),
        Mean_T_Tm=("T_homologous", "mean"),
    )
    if crystallite_col:
        zone_stats["Mean_Crystallite"] = df.groupby("Structure_Zone")[crystallite_col].mean()
    zone_stats.to_csv("results/structure_zone_statistics.csv")

    return df, zone_stats
```

## 3. XRD 結晶子サイズ解析（Scherrer 方程式）

```python
def scherrer_crystallite_size(beta_rad, two_theta_deg, K=0.9, wavelength_nm=0.15406):
    """
    Scherrer 方程式: D = Kλ / (β cos θ)

    Parameters:
        beta_rad: 半値全幅 FWHM (ラジアン)
        two_theta_deg: ブラッグ角 2θ (度)
        K: 形状因子 (通常 0.9)
        wavelength_nm: X 線波長 (nm, Cu-Kα = 0.15406)

    Returns:
        結晶子サイズ D (nm)
    """
    theta_rad = np.radians(two_theta_deg / 2)
    D = (K * wavelength_nm) / (beta_rad * np.cos(theta_rad))
    return D


def williamson_hall_analysis(two_theta_list, fwhm_list, wavelength_nm=0.15406):
    """
    Williamson-Hall プロット: β cos θ = Kλ/D + 4ε sin θ

    FWHM のブロードニングを結晶子サイズ効果と格子歪み効果に分離する。

    Returns:
        D: 結晶子サイズ (nm)
        epsilon: 格子歪み (%)
    """
    theta_rad = np.radians(np.array(two_theta_list) / 2)
    beta_rad = np.array(fwhm_list)

    x = 4 * np.sin(theta_rad)
    y = beta_rad * np.cos(theta_rad)

    # 線形回帰: y = intercept + slope * x
    from scipy.stats import linregress
    slope, intercept, r_value, p_value, std_err = linregress(x, y)

    D = (0.9 * wavelength_nm) / intercept if intercept > 0 else np.nan
    epsilon = slope * 100  # %

    return {
        "crystallite_size_nm": D,
        "lattice_strain_pct": epsilon,
        "r_squared": r_value ** 2,
        "p_value": p_value,
    }
```

## 4. 配向度 TC(hkl) の計算

```python
def texture_coefficient(I_measured, I_reference):
    """
    配向度 Texture Coefficient TC(hkl) を計算する。

    TC(hkl) = [I(hkl)/I₀(hkl)] / [(1/N) Σ I(hkl)/I₀(hkl)]

    TC = 1: ランダム配向
    TC > 1: 優先配向
    TC < 1: 劣位配向

    Parameters:
        I_measured: dict {(hkl): intensity} — 測定ピーク強度
        I_reference: dict {(hkl): intensity} — ICDD PDF カード参照強度
    """
    ratios = {}
    for hkl in I_measured:
        if hkl in I_reference and I_reference[hkl] > 0:
            ratios[hkl] = I_measured[hkl] / I_reference[hkl]

    N = len(ratios)
    if N == 0:
        return {}

    mean_ratio = np.mean(list(ratios.values()))

    tc = {hkl: ratio / mean_ratio for hkl, ratio in ratios.items()}
    return tc
```

## 5. 膜応力解析（Stoney 方程式）

```python
def stoney_film_stress(R_before_m, R_after_m, E_substrate_Pa,
                        nu_substrate, t_substrate_m, t_film_m):
    """
    Stoney 方程式: σ = (Es × ts²) / (6(1-νs) × tf) × (1/R_after - 1/R_before)

    曲率変化から膜応力を算出する。

    Parameters:
        R_before_m: 成膜前の曲率半径 (m)
        R_after_m: 成膜後の曲率半径 (m)
        E_substrate_Pa: 基板のヤング率 (Pa)
        nu_substrate: 基板のポアソン比
        t_substrate_m: 基板厚 (m)
        t_film_m: 膜厚 (m)

    Returns:
        膜応力 σ (Pa) — 正:引張、負:圧縮
    """
    curvature_change = (1 / R_after_m) - (1 / R_before_m)
    biaxial_modulus = E_substrate_Pa / (1 - nu_substrate)
    sigma = (biaxial_modulus * t_substrate_m**2 * curvature_change) / (6 * t_film_m)
    return sigma


# Si(100) 基板の標準定数
SI_SUBSTRATE = {
    "E_Pa": 130.2e9,        # ヤング率
    "nu": 0.279,            # ポアソン比
    "t_m": 525e-6,          # 標準ウエハ厚 525 μm
}
```

## 6. 多技法融合データ統合

```python
def merge_characterization_data(xrd_df, afm_df, electrical_df, optical_df,
                                  sample_id_col="Sample_ID"):
    """
    複数の測定手法のデータをサンプル ID で統合する。

    XRD → 結晶子サイズ, 配向度, 格子歪
    AFM → 表面粗さ Ra, Rq
    電気 → 比抵抗, シート抵抗
    光学 → 透過率, バンドギャップ
    """
    merged = xrd_df.copy()
    for df in [afm_df, electrical_df, optical_df]:
        if df is not None and sample_id_col in df.columns:
            merged = merged.merge(df, on=sample_id_col, how="left",
                                 suffixes=("", "_dup"))
            # 重複カラムを除去
            dup_cols = [c for c in merged.columns if c.endswith("_dup")]
            merged.drop(columns=dup_cols, inplace=True)

    return merged


def tauc_plot_bandgap(wavelength_nm, transmittance_pct, thickness_nm,
                       n_exponent=2, figsize=(8, 6)):
    """
    Tauc プロットからバンドギャップを推定する。

    (αhν)^(1/n) vs hν のプロットの線形外挿。
    n=2: 直接遷移, n=0.5: 間接遷移

    Parameters:
        wavelength_nm: 波長 (nm)
        transmittance_pct: 透過率 (%)
        thickness_nm: 膜厚 (nm)
        n_exponent: 遷移タイプ (2=直接, 0.5=間接)
    """
    T = transmittance_pct / 100
    alpha = -np.log(T + 1e-10) / (thickness_nm * 1e-7)  # cm⁻¹

    h = 6.626e-34   # Planck 定数 (J·s)
    c = 3e8          # 光速 (m/s)
    eV = 1.602e-19   # eV → J

    energy_eV = (h * c) / (wavelength_nm * 1e-9) / eV  # hν (eV)
    tauc = (alpha * energy_eV) ** (1 / n_exponent)

    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(energy_eV, tauc, "b-", linewidth=1.5)
    ax.set_xlabel("Photon Energy hν (eV)")
    ax.set_ylabel(f"(αhν)^(1/{n_exponent})")
    ax.set_title("Tauc Plot", fontweight="bold")
    plt.tight_layout()
    plt.savefig("figures/tauc_plot.png", dpi=300, bbox_inches="tight")
    plt.close()

    return energy_eV, tauc
```

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `materials_project` | Materials Project | 材料特性データベース検索 |

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/structure_zone_statistics.csv` | CSV |
| `results/xrd_analysis.csv` | CSV |
| `results/williamson_hall.csv` | CSV |
| `figures/structure_zone_model.png` | PNG |
| `figures/xrd_analysis.png` | PNG |
| `figures/tauc_plot.png` | PNG |

#### 参照実験

- **Exp-13**: 薄膜 PSP（Thornton-Anders SZM、XRD 結晶子サイズ、Stoney 膜応力）
- **Exp-12**: マテリアルサイエンスのプロセスデータ解析
- **Exp-11**: ARIM データポータルの材料データ活用

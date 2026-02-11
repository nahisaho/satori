---
name: scientific-data-simulation
description: |
  物理・化学・生物学に基づく合成データ生成のスキル。実験データが不足する場合に、
  ドメイン知識を反映したシミュレーションデータを生成する際に使用。
  Scientific Skills Exp-06, 07, 08, 09, 12, 13 で確立したパターン。
---

# Scientific Data Simulation & Generation

物理法則・化学モデル・生物学的知見に基づいて、現実的なシミュレーションデータを
生成するスキル。実データの不足を補い、ML パイプラインの開発・検証を可能にする。

## When to Use

- 実験データが未取得または不十分なとき
- ML パイプラインのプロトタイピングに合成データが必要なとき
- 特定の物理・化学法則に基づくデータ生成
- 既知の因果関係を組み込んだベンチマークデータの作成

## Quick Start

## 設計原則

1. **物理モデルベース**: 単純な正規分布ではなく、ドメインの因果関係を反映する
2. **ノイズの付加**: 実験の不確実性を模擬する適切なノイズレベル
3. **材料/群の多様性**: 複数の材料種・条件を生成し、群間差を反映
4. **範囲の現実性**: パラメータ範囲は実験装置・物理限界に基づく

## データ生成テンプレート

### 1. プロセスデータ生成（Exp-12, 13 パターン）

```python
import numpy as np
import pandas as pd

def generate_process_dataset(n_samples=500, seed=42):
    """
    物理ベースのプロセスデータを生成するテンプレート。
    因果関係: Process → Structure → Property
    """
    rng = np.random.default_rng(seed)

    # === Process パラメータ（独立変数）===
    temperature = rng.uniform(25, 500, n_samples)      # °C
    pressure = rng.uniform(0.1, 5.0, n_samples)        # Pa
    power = rng.uniform(50, 500, n_samples)             # W
    time = rng.uniform(5, 120, n_samples)               # min

    # === Structure 変数（Process に依存）===
    # 物理モデル：因果関係を組み込む
    dep_rate = (
        0.5
        + 0.02 * power                    # 出力依存
        - 0.005 * pressure ** 2           # 高圧での飽和
        + rng.normal(0, 0.5, n_samples)   # ノイズ
    )
    dep_rate = np.clip(dep_rate, 0.1, 30)

    thickness = dep_rate * time
    thickness = np.clip(thickness, 5, 2000)

    crystallite_size = (
        5
        + 0.1 * temperature               # 温度依存（アレニウス的）
        + 0.01 * time                      # 時間依存
        + rng.normal(0, 2, n_samples)      # ノイズ
    )
    crystallite_size = np.clip(crystallite_size, 2, 80)

    # === Property 変数（Structure に依存）===
    resistivity = (
        1e-2
        * np.exp(-0.005 * temperature)      # 温度活性化
        * (1 + 0.01 * pressure)
        * np.exp(rng.normal(0, 0.3, n_samples))  # 対数正規ノイズ
    )

    transmittance = (
        95
        - 0.02 * thickness                  # 膜厚依存
        + 0.05 * crystallite_size            # 結晶子サイズ依存
        + rng.normal(0, 1, n_samples)
    )
    transmittance = np.clip(transmittance, 40, 98)

    df = pd.DataFrame({
        "Temperature": temperature,
        "Pressure": pressure,
        "Power": power,
        "Time": time,
        "Deposition_Rate": dep_rate,
        "Thickness": thickness,
        "Crystallite_Size": crystallite_size,
        "Resistivity": resistivity,
        "Transmittance": transmittance,
    })

    return df
```

### 2. 多材料データ生成（Exp-13 パターン）

```python
def generate_multi_material_dataset(materials, n_per_material=100, seed=42):
    """
    複数材料の PSP データを生成する。
    materials: {"ZnO": {"Tm": 2248, "Eg": 3.3}, ...} の辞書
    """
    rng = np.random.default_rng(seed)
    all_data = []

    for mat_name, props in materials.items():
        n = n_per_material
        Tm = props["Tm"]  # 融点 (K)
        Eg = props["Eg"]  # バンドギャップ (eV)

        # Process
        Tsub = rng.uniform(25, 500, n)
        Pwork = rng.uniform(0.1, 5.0, n)
        Power = rng.uniform(50, 500, n)

        # Structure（材料依存の因果関係）
        T_homologous = (Tsub + 273.15) / Tm  # 相同温度
        crystallite = 5 + 80 * T_homologous + rng.normal(0, 3, n)
        crystallite = np.clip(crystallite, 2, 80)

        # Property（材料固有値 + プロセス依存変動）
        bandgap = Eg + rng.normal(0, 0.1, n)

        data = pd.DataFrame({
            "Material": mat_name,
            "Substrate_Temp": Tsub,
            "Working_Pressure": Pwork,
            "Power": Power,
            "T_homologous": T_homologous,
            "Crystallite_Size": crystallite,
            "Bandgap": bandgap,
        })
        all_data.append(data)

    return pd.concat(all_data, ignore_index=True)
```

### 3. 臨床試験データ生成（Exp-06 パターン）

```python
def generate_clinical_trial_data(n_total=500, effect_size=0.3, seed=42):
    """RCT シミュレーションデータを生成する。"""
    rng = np.random.default_rng(seed)
    n_per_arm = n_total // 2

    # 人口統計
    ages = np.concatenate([
        rng.normal(55, 12, n_per_arm),
        rng.normal(55, 12, n_per_arm),
    ])
    sex = rng.choice(["M", "F"], n_total)
    group = np.array(["Treatment"] * n_per_arm + ["Control"] * n_per_arm)

    # 主要エンドポイント
    baseline = rng.normal(100, 15, n_total)
    treatment_effect = np.where(group == "Treatment", effect_size * 15, 0)
    endpoint = baseline + treatment_effect + rng.normal(0, 10, n_total)

    # 生存時間（指数分布）
    survival_time = rng.exponential(
        np.where(group == "Treatment", 365 * 2, 365 * 1.5),
        n_total
    )
    event = rng.binomial(1, 0.7, n_total)

    return pd.DataFrame({
        "Patient_ID": range(1, n_total + 1),
        "Group": group,
        "Age": ages.astype(int),
        "Sex": sex,
        "Baseline": baseline,
        "Endpoint": endpoint,
        "Survival_Time": survival_time,
        "Event": event,
    })
```

### 4. スペクトルデータ生成（Exp-08, 11 パターン）

```python
def generate_spectrum(wavenumbers, peak_positions, peak_heights,
                       peak_widths, noise_level=0.02, seed=None):
    """
    ガウスピーク合成によるスペクトルを生成する。
    ラマン / IR / UV-Vis などに汎用。
    """
    rng = np.random.default_rng(seed)
    spectrum = np.zeros_like(wavenumbers, dtype=float)

    for pos, height, width in zip(peak_positions, peak_heights, peak_widths):
        spectrum += height * np.exp(-0.5 * ((wavenumbers - pos) / width) ** 2)

    # ノイズ付加
    spectrum += rng.normal(0, noise_level * spectrum.max(), len(wavenumbers))
    return spectrum


def generate_ecg_beat(t, hr=72):
    """合成 ECG 波形（PQRST パターン）を生成する。"""
    beat_duration = 60.0 / hr
    # P 波、QRS 群、T 波のガウス重ね合わせ
    p_wave = 0.1 * np.exp(-((t % beat_duration - 0.16) / 0.04) ** 2)
    qrs = 1.0 * np.exp(-((t % beat_duration - 0.25) / 0.01) ** 2)
    q_wave = -0.15 * np.exp(-((t % beat_duration - 0.22) / 0.015) ** 2)
    s_wave = -0.1 * np.exp(-((t % beat_duration - 0.28) / 0.015) ** 2)
    t_wave = 0.2 * np.exp(-((t % beat_duration - 0.40) / 0.05) ** 2)
    return p_wave + q_wave + qrs + s_wave + t_wave
```

## 品質チェックリスト

生成データの品質を保証するためのチェック項目：

- [ ] パラメータ範囲が物理的に現実的か
- [ ] 因果関係の方向が正しいか（温度↑ → 結晶子サイズ↑ など）
- [ ] ノイズレベルが実験の再現性に対応するか
- [ ] 外れ値の割合が現実的か（通常 1-5%）
- [ ] 変数間の相関構造が既知の物理法則と一致するか
- [ ] 群間差が効果量として検出可能な水準か

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `data/<dataset_name>.csv` | CSV |

#### 参照実験

- **Exp-06**: 臨床試験 RCT シミュレーション（500 名、2 群）
- **Exp-07**: メタボロミクス合成データ（100 サンプル × 200 代謝物）
- **Exp-08**: 合成 ECG/EEG 信号（PQRST、帯域合成）
- **Exp-09**: コドンバイアス反映ゲノム配列
- **Exp-12**: エッチングプロセスデータ（500 サンプル × 8 パラメータ）
- **Exp-13**: 薄膜成膜データ（600 サンプル × 6 材料 × PSP 3 階層）

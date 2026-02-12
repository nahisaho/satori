---
name: scientific-admet-pharmacokinetics
description: |
  ADMET 予測・薬物動態モデリングスキル。吸収(A)・分布(D)・代謝(M)・排泄(E)・毒性(T)の
  包括的予測パイプライン。DeepChem/ADMET-AI/PyTDC を活用した分子特性予測、
  PK/PD モデリング、ドラッグライクネス最適化、リード最適化戦略を提供。
  「ADMET 予測して」「薬物動態を評価して」「lead optimization して」で発火。
tu_tools:
  - key: pubchem
    name: PubChem
    description: 化合物・物質・生理活性アッセイデータベース
---

# Scientific ADMET & Pharmacokinetics

ADMET 予測と薬物動態モデリングのスキル。創薬初期段階における化合物の
薬物動態特性評価を包括的に支援する。

## When to Use

- 化合物の ADMET 特性を統合的に予測するとき
- リード化合物の最適化方針を策定するとき
- PK/PD パラメータを推定・モデリングするとき
- 化合物ライブラリのドラッグライクネスフィルタリングを行うとき
- 毒性リスク（hERG、Ames、肝毒性）を予測するとき

## Quick Start

### 1. ADMET 予測パイプライン概要

```
Input: SMILES / SDF
    ↓
Step 1: Absorption
  - Caco-2 透過性
  - HIA (Human Intestinal Absorption)
  - 経口バイオアベイラビリティ
  - Pgp 基質判定
    ↓
Step 2: Distribution
  - 血漿タンパク結合率 (PPB)
  - 血液脳関門 (BBB) 透過性
  - 分布容積 (VDss)
    ↓
Step 3: Metabolism
  - CYP 阻害/基質予測 (1A2, 2C9, 2C19, 2D6, 3A4)
  - 代謝安定性 (Half-life)
  - クリアランス予測
    ↓
Step 4: Excretion
  - 腎クリアランス
  - 総クリアランス
  - 排泄経路予測
    ↓
Step 5: Toxicity
  - hERG 阻害 (心毒性)
  - Ames 試験予測 (変異原性)
  - DILI (薬剤性肝障害)
  - LD50 急性毒性
    ↓
Output: ADMET Profile Card
```

---

## Phase 1: Absorption 予測

### Caco-2 透過性 & HIA

```python
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors

def predict_absorption_properties(smiles):
    """
    吸収関連パラメータの計算・予測。
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    # 物理化学パラメータ（吸収関連）
    properties = {
        "MW": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "HBA": Descriptors.NumHAcceptors(mol),
        "HBD": Descriptors.NumHDonors(mol),
        "RotBonds": Descriptors.NumRotatableBonds(mol),
    }

    # Rule-based absorption prediction
    properties["Lipinski_Violations"] = sum([
        properties["MW"] > 500,
        properties["LogP"] > 5,
        properties["HBA"] > 10,
        properties["HBD"] > 5,
    ])

    # TPSA-based intestinal absorption estimate
    # Ertl et al., J. Med. Chem. 2000
    if properties["TPSA"] < 140:
        properties["HIA_prediction"] = "High"
    else:
        properties["HIA_prediction"] = "Low"

    # BBB penetration estimate (TPSA < 90 Å²)
    properties["BBB_prediction"] = "Permeable" if properties["TPSA"] < 90 else "Non-permeable"

    return properties
```

---

## Phase 2: ドラッグライクネス評価

### 多重フィルタ評価

```python
def evaluate_druglikeness(smiles):
    """
    複数フィルタによるドラッグライクネス総合評価。
    """
    mol = Chem.MolFromSmiles(smiles)
    results = {}

    # Lipinski Rule of 5
    results["Lipinski"] = {
        "MW": Descriptors.MolWt(mol) <= 500,
        "LogP": Descriptors.MolLogP(mol) <= 5,
        "HBA": Descriptors.NumHAcceptors(mol) <= 10,
        "HBD": Descriptors.NumHDonors(mol) <= 5,
        "pass": None,  # True if ≤1 violation
    }
    violations = sum(1 for v in list(results["Lipinski"].values())[:4] if not v)
    results["Lipinski"]["pass"] = violations <= 1

    # Veber Rules (oral bioavailability)
    results["Veber"] = {
        "RotBonds": Descriptors.NumRotatableBonds(mol) <= 10,
        "TPSA": Descriptors.TPSA(mol) <= 140,
        "pass": None,
    }
    results["Veber"]["pass"] = all(
        v for k, v in results["Veber"].items() if k != "pass"
    )

    # QED (Quantitative Estimate of Drug-likeness)
    from rdkit.Chem import QED
    results["QED"] = QED.qed(mol)

    # Ghose Filter
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    mr = Descriptors.MolMR(mol)
    n_atoms = mol.GetNumAtoms()
    results["Ghose"] = {
        "pass": (160 <= mw <= 480 and -0.4 <= logp <= 5.6
                 and 40 <= mr <= 130 and 20 <= n_atoms <= 70)
    }

    return results
```

---

## Phase 3: 毒性予測

### 構造アラート & 毒性エンドポイント

```python
# PAINS (Pan-Assay Interference Compounds) フィルタ
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

def check_structural_alerts(smiles):
    """
    構造アラート (PAINS, Brenk) の検出。
    """
    mol = Chem.MolFromSmiles(smiles)
    alerts = []

    # PAINS フィルタ
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)
    if catalog.HasMatch(mol):
        entry = catalog.GetFirstMatch(mol)
        alerts.append({
            "type": "PAINS",
            "description": entry.GetDescription(),
            "severity": "Warning",
        })

    # Brenk フィルタ（毒性懸念構造）
    params_brenk = FilterCatalogParams()
    params_brenk.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    catalog_brenk = FilterCatalog(params_brenk)
    if catalog_brenk.HasMatch(mol):
        entry = catalog_brenk.GetFirstMatch(mol)
        alerts.append({
            "type": "Brenk",
            "description": entry.GetDescription(),
            "severity": "Warning",
        })

    return alerts
```

### 毒性リスクスコアカード

```markdown
## Toxicity Risk Assessment

| Endpoint       | Prediction | Confidence | Risk Level |
|----------------|------------|------------|------------|
| hERG IC50      |            |            | Low/Med/High |
| Ames           |            |            | Low/Med/High |
| DILI           |            |            | Low/Med/High |
| LD50 (rat)     |            |            | Low/Med/High |
| Skin Sensitiz. |            |            | Low/Med/High |
| Carcinogenicity|            |            | Low/Med/High |

### Structural Alerts
| Alert Type | Description | Severity |
|------------|-------------|----------|
```

---

## Phase 4: PK/PD モデリング

### コンパートメントモデル

```python
import numpy as np
from scipy.integrate import odeint

def one_compartment_pk(dose, ka, ke, vd, time_points):
    """
    1-コンパートメント PK モデル（経口投与）。
    dose: 投与量 (mg)
    ka: 吸収速度定数 (h⁻¹)
    ke: 消失速度定数 (h⁻¹)
    vd: 分布容積 (L)
    """
    def pk_ode(y, t, ka, ke):
        dAg_dt = -ka * y[0]          # 消化管
        dCp_dt = (ka * y[0] - ke * y[1] * vd) / vd  # 血漿中
        return [dAg_dt, dCp_dt]

    y0 = [dose, 0]
    solution = odeint(pk_ode, y0, time_points, args=(ka, ke))

    return {
        "time": time_points,
        "concentration": solution[:, 1],
        "Cmax": np.max(solution[:, 1]),
        "Tmax": time_points[np.argmax(solution[:, 1])],
        "AUC": np.trapz(solution[:, 1], time_points),
        "half_life": np.log(2) / ke,
    }
```

---

## ADMET Profile Card テンプレート

```markdown
# ADMET Profile: [COMPOUND NAME]

**SMILES**: [canonical SMILES]
**Molecular Formula**: [formula]
**Date**: [date]

## 1. Physicochemical Properties
| Property | Value | Optimal Range | Status |
|----------|-------|---------------|--------|
| MW       |       | 150-500       | ✓/✗   |
| LogP     |       | -0.4 to 5.6   | ✓/✗   |
| TPSA     |       | 20-140 Å²     | ✓/✗   |
| HBA      |       | ≤10           | ✓/✗   |
| HBD      |       | ≤5            | ✓/✗   |
| QED      |       | >0.5          | ✓/✗   |

## 2. Absorption
| Parameter | Value | Interpretation |
|-----------|-------|----------------|

## 3. Distribution
| Parameter | Value | Interpretation |
|-----------|-------|----------------|

## 4. Metabolism
| CYP Enzyme | Substrate? | Inhibitor? |
|------------|------------|------------|
| 1A2        |            |            |
| 2C9        |            |            |
| 2C19       |            |            |
| 2D6        |            |            |
| 3A4        |            |            |

## 5. Excretion
| Parameter | Value | Interpretation |
|-----------|-------|----------------|

## 6. Toxicity
| Endpoint | Prediction | Risk |
|----------|------------|------|

## 7. Druglikeness Summary
| Filter    | Pass? | Violations |
|-----------|-------|------------|
| Lipinski  |       |            |
| Veber     |       |            |
| Ghose     |       |            |
| PAINS     |       |            |
| Brenk     |       |            |

## 8. Recommendations
- [ ] Lead optimization priorities
- [ ] Key liabilities to address
- [ ] Suggested structural modifications
```

---

## Completeness Checklist

- [ ] 物理化学的性質: MW, LogP, TPSA, HBA, HBD, QED
- [ ] 吸収: Caco-2, HIA, Pgp 基質
- [ ] 分布: PPB, BBB, VDss
- [ ] 代謝: CYP 5 isoform (1A2/2C9/2C19/2D6/3A4)
- [ ] 排泄: クリアランス, Half-life
- [ ] 毒性: hERG, Ames, DILI ≥ 3 エンドポイント
- [ ] 構造アラート: PAINS + Brenk チェック
- [ ] ドラッグライクネス: Lipinski + Veber + QED

## Best Practices

1. **予測値は必ず信頼度を付与**: モデルの applicability domain を確認
2. **複数ツールで Cross-validate**: ADMET-AI, DeepChem, SwissADME の結果を比較
3. **既知化合物で Benchmark**: 同クラスの承認薬と比較して判断
4. **構造アラートは参考情報**: PAINS ヒットでも偽陽性の可能性を考慮
5. **PK パラメータは in vivo 補正**: allometric scaling で動物種間換算

## References

### Output Files

| ファイル | 形式 | 生成タイミング |
|---|---|---|
| `results/admet_profile.json` | ADMET プロファイル（JSON） | 全 5 段階評価完了時 |
| `results/admet_report.md` | ADMET 評価レポート（Markdown） | 全解析完了時 |
| `results/pk_model.json` | PK モデルパラメータ（JSON） | PK モデリング完了時 |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| ADMET-AI | `ADMETAI_predict_BBB_penetrance` | BBB 透過性予測 |
| ADMET-AI | `ADMETAI_predict_CYP_interactions` | CYP 相互作用予測 |
| ADMET-AI | `ADMETAI_predict_toxicity` | 毒性予測 |
| ADMET-AI | `ADMETAI_predict_bioavailability` | バイオアベイラビリティ予測 |
| PubChem | `PubChem_get_compound_properties_by_CID` | 化合物物性取得 |
| ChEMBL | `ChEMBL_get_molecule` | 分子情報取得 |
| ChEMBL | `ChEMBL_get_activity` | バイオアッセイデータ |

### 参照スキル

| スキル | 連携 |
|---|---|
| `scientific-drug-target-profiling` | ← ターゲットに結合する化合物の ADMET 評価 |
| `scientific-cheminformatics` | ← 分子記述子・構造情報の提供 |
| `scientific-drug-repurposing` | → ADMET 通過化合物のリポジショニング候補評価 |
| `scientific-clinical-decision-support` | → PK パラメータの臨床応用 |
| `scientific-academic-writing` | → 研究成果の論文化 |
| `scientific-regulatory-science` | → FDA 規制申請・承認履歴 |

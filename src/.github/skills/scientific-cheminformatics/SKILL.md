---
name: scientific-cheminformatics
description: |
  ケモインフォマティクス解析のスキル。RDKit を用いた分子記述子計算、Morgan フィンガープリント、
  Tanimoto 類似度、構造アラート検出、Lipinski Rule of 5 評価を行う際に使用。
  Scientific Skills Exp-02, 05 で確立したパターン。
---

# Scientific Cheminformatics Analysis

RDKit を用いた分子解析パイプラインスキル。SMILES → 分子記述子 → SAR 解析 →
毒性予測までの創薬ケモインフォマティクスワークフローを提供する。

## When to Use

- 化合物の物理化学的性質を算出したいとき
- SMILES 文字列から分子記述子を計算したいとき
- 化合物間の構造類似度を評価したいとき
- 構造活性相関（SAR）を解析したいとき
- 構造アラート（トキシコフォア）を検出したいとき
- Lipinski Rule of 5 / ドラッグライクネスを評価したいとき

## Quick Start

## 標準パイプライン

### 1. SMILES → 分子オブジェクト変換

```python
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, QED, Lipinski
from rdkit.Chem.Scaffolds import MurckoScaffold
import pandas as pd
import numpy as np

def smiles_to_mol(smiles):
    """SMILES 文字列から RDKit 分子オブジェクトを生成する。"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return mol
```

### 2. 分子記述子の一括計算

```python
def calculate_descriptors(smiles_list, names=None):
    """
    SMILES リストから主要な分子記述子を一括計算する。
    返値: DataFrame
    """
    records = []
    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue

        record = {
            "Name": names[i] if names else f"Mol_{i}",
            "SMILES": smi,
            "MW": Descriptors.MolWt(mol),
            "LogP": Descriptors.MolLogP(mol),
            "TPSA": Descriptors.TPSA(mol),
            "HBA": Descriptors.NumHAcceptors(mol),
            "HBD": Descriptors.NumHDonors(mol),
            "RotBonds": Descriptors.NumRotatableBonds(mol),
            "AromaticRings": Descriptors.NumAromaticRings(mol),
            "HeavyAtoms": mol.GetNumHeavyAtoms(),
            "QED": QED.qed(mol),
            "Fraction_CSP3": Descriptors.FractionCSP3(mol),
        }
        records.append(record)

    return pd.DataFrame(records)
```

### 3. Morgan フィンガープリント & Tanimoto 類似度

```python
from rdkit import DataStructs

def compute_fingerprints(smiles_list, radius=2, nBits=2048):
    """Morgan フィンガープリントを生成する。"""
    fps = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
            fps.append(fp)
    return fps

def tanimoto_similarity_matrix(fps, names=None):
    """Tanimoto 類似度行列を算出する。"""
    n = len(fps)
    sim_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            sim_matrix[i, j] = DataStructs.TanimotoSimilarity(fps[i], fps[j])

    if names is None:
        names = [f"Mol_{i}" for i in range(n)]

    sim_df = pd.DataFrame(sim_matrix, index=names, columns=names)
    sim_df.to_csv("results/tanimoto_similarity.csv")
    return sim_df
```

### 4. Lipinski Rule of 5 評価

```python
def lipinski_evaluation(desc_df):
    """Lipinski Rule of 5 の準拠チェック。"""
    desc_df = desc_df.copy()
    desc_df["Lipinski_MW"] = desc_df["MW"] <= 500
    desc_df["Lipinski_LogP"] = desc_df["LogP"] <= 5
    desc_df["Lipinski_HBA"] = desc_df["HBA"] <= 10
    desc_df["Lipinski_HBD"] = desc_df["HBD"] <= 5
    desc_df["Lipinski_Violations"] = 4 - (
        desc_df["Lipinski_MW"].astype(int) +
        desc_df["Lipinski_LogP"].astype(int) +
        desc_df["Lipinski_HBA"].astype(int) +
        desc_df["Lipinski_HBD"].astype(int)
    )
    desc_df["Lipinski_Pass"] = desc_df["Lipinski_Violations"] <= 1
    return desc_df
```

### 5. 構造アラート（トキシコフォア）検出（Exp-05）

```python
STRUCTURAL_ALERTS = {
    "Nitro": "[N+](=O)[O-]",
    "Epoxide": "C1OC1",
    "Aldehyde": "[CH]=O",
    "Michael_Acceptor": "C=CC(=O)",
    "Acyl_Halide": "C(=O)[F,Cl,Br,I]",
    "Aniline": "c1ccccc1N",
    "Hydrazine": "NN",
    "Sulfonate": "S(=O)(=O)[O-]",
}

def detect_structural_alerts(smiles_list, names=None, alerts=None):
    """SMARTS パターンによる構造アラートの検出。"""
    if alerts is None:
        alerts = STRUCTURAL_ALERTS

    results = []
    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue

        name = names[i] if names else f"Mol_{i}"
        for alert_name, smarts in alerts.items():
            pattern = Chem.MolFromSmarts(smarts)
            if mol.HasSubstructMatch(pattern):
                results.append({"Name": name, "SMILES": smi,
                               "Alert": alert_name, "SMARTS": smarts})

    return pd.DataFrame(results)
```

### 6. Murcko スキャフォールド解析

```python
def scaffold_analysis(smiles_list, names=None):
    """Murcko スキャフォールドの抽出と分類。"""
    scaffolds = []
    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smi)
        if mol:
            core = MurckoScaffold.GetScaffoldForMol(mol)
            scaffolds.append({
                "Name": names[i] if names else f"Mol_{i}",
                "SMILES": smi,
                "Scaffold": Chem.MolToSmiles(core),
            })
    return pd.DataFrame(scaffolds)
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/molecular_properties.csv` | CSV |
| `results/tanimoto_similarity.csv` | CSV |
| `results/structural_alerts.csv` | CSV |
| `figures/chemical_space_pca.png` | PNG |
| `figures/similarity_heatmap.png` | PNG |

#### 参照実験

- **Exp-02**: EGFR 阻害剤 SAR 解析（記述子、Tanimoto、MCS、Scaffold）
- **Exp-05**: 毒性予測（構造アラート、Morgan FP 分類モデル）

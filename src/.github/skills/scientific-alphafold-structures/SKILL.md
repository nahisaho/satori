---
name: scientific-alphafold-structures
description: |
  AlphaFold 構造予測スキル。AlphaFold Protein Structure Database
  REST API による予測構造取得・pLDDT 信頼度解析・PAE 残基間
  距離予測・構造カバレッジ分析。ToolUniverse 連携: alphafold。
tu_tools:
  - key: alphafold
    name: AlphaFold Database
    description: AlphaFold 予測構造・コンフィデンス・PAE 取得
---

# Scientific AlphaFold Structures

AlphaFold Protein Structure Database REST API を活用した
構造予測取得・信頼度解析パイプラインを提供する。

## When to Use

- UniProt ID から AlphaFold 予測構造を取得するとき
- pLDDT スコアで構造信頼度を評価するとき
- PAE (Predicted Aligned Error) で残基間予測精度を分析するとき
- 構造カバレッジ (モデル化割合) を確認するとき
- AlphaFold 構造を実験構造と比較するとき
- 大規模プロテオーム構造データをバッチ取得するとき

---

## Quick Start

## 1. AlphaFold 予測構造取得

```python
import requests
import pandas as pd
import numpy as np
from io import StringIO

AFDB_BASE = "https://alphafold.ebi.ac.uk/api"


def alphafold_get_prediction(uniprot_id, version=4):
    """
    AlphaFold DB — 予測構造取得。

    Parameters:
        uniprot_id: str — UniProt アクセッション (例: "P00533")
        version: int — AlphaFold モデルバージョン
    """
    url = f"{AFDB_BASE}/prediction/{uniprot_id}"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    entries = resp.json()

    if not entries:
        print(f"No AlphaFold prediction for {uniprot_id}")
        return None

    entry = entries[0] if isinstance(entries, list) else entries
    result = {
        "uniprot_id": uniprot_id,
        "entry_id": entry.get("entryId", ""),
        "gene": entry.get("gene", ""),
        "organism": entry.get("organismScientificName", ""),
        "tax_id": entry.get("taxId", ""),
        "sequence_length": entry.get("uniprotEnd", 0)
                           - entry.get("uniprotStart", 0) + 1,
        "model_url": entry.get("cifUrl", ""),
        "pdb_url": entry.get("pdbUrl", ""),
        "pae_url": entry.get("paeImageUrl", ""),
        "global_plddt": entry.get("globalMetricValue", None),
        "model_version": entry.get("latestVersion", version),
    }

    print(f"AlphaFold {uniprot_id}: pLDDT={result['global_plddt']}, "
          f"length={result['sequence_length']}")
    return result
```

## 2. pLDDT 信頼度プロファイル解析

```python
import biotite.structure.io.pdbx as pdbx
import biotite.structure as struc


def alphafold_plddt_profile(uniprot_id, output_dir="results"):
    """
    AlphaFold — pLDDT 残基別信頼度プロファイル。

    Parameters:
        uniprot_id: str — UniProt アクセッション
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # CIF 取得
    pred = alphafold_get_prediction(uniprot_id)
    if not pred:
        return pd.DataFrame()

    cif_url = pred["model_url"]
    resp = requests.get(cif_url, timeout=60)
    resp.raise_for_status()

    cif_path = output_dir / f"AF-{uniprot_id}.cif"
    cif_path.write_bytes(resp.content)

    # pLDDT = B-factor in AlphaFold structures
    pdbx_file = pdbx.CIFFile.read(str(cif_path))
    structure = pdbx.get_structure(pdbx_file, model=1)
    ca_mask = structure.atom_name == "CA"
    ca_atoms = structure[ca_mask]

    residues = []
    for i, atom in enumerate(ca_atoms):
        residues.append({
            "residue_index": i + 1,
            "residue_name": atom.res_name,
            "chain": atom.chain_id,
            "plddt": atom.b_factor,
        })

    df = pd.DataFrame(residues)

    # 信頼度カテゴリ
    df["confidence"] = pd.cut(
        df["plddt"],
        bins=[0, 50, 70, 90, 100],
        labels=["Very low", "Low", "Confident", "Very high"],
    )

    n_high = (df["plddt"] >= 70).sum()
    print(f"pLDDT profile {uniprot_id}: {len(df)} residues, "
          f"{n_high}/{len(df)} confident (≥70)")
    return df
```

## 3. PAE (Predicted Aligned Error) 解析

```python
def alphafold_pae_analysis(uniprot_id):
    """
    AlphaFold — PAE マトリクス解析。

    Parameters:
        uniprot_id: str — UniProt アクセッション
    """
    url = (f"https://alphafold.ebi.ac.uk/files/"
           f"AF-{uniprot_id}-F1-predicted_aligned_error_v4.json")
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    pae_data = data[0] if isinstance(data, list) else data
    pae_matrix = np.array(pae_data.get("predicted_aligned_error",
                                         pae_data.get("pae", [])))

    # ドメイン境界推定 (低 PAE ブロック検出)
    n_res = pae_matrix.shape[0]
    mean_pae = pae_matrix.mean()
    low_pae_mask = pae_matrix < mean_pae * 0.5

    # 行ごとの低 PAE 連続領域
    domain_scores = []
    for i in range(n_res):
        row = low_pae_mask[i]
        domain_scores.append(row.sum() / n_res)

    result = {
        "uniprot_id": uniprot_id,
        "pae_matrix_shape": pae_matrix.shape,
        "mean_pae": float(mean_pae),
        "median_pae": float(np.median(pae_matrix)),
        "min_pae": float(pae_matrix.min()),
        "max_pae": float(pae_matrix.max()),
        "domain_scores": domain_scores,
    }

    print(f"PAE {uniprot_id}: {n_res}x{n_res} matrix, "
          f"mean={mean_pae:.1f} Å")
    return result, pae_matrix
```

## 4. AlphaFold 統合パイプライン

```python
def alphafold_pipeline(uniprot_ids, output_dir="results"):
    """
    AlphaFold 構造解析統合パイプライン。

    Parameters:
        uniprot_ids: list[str] — UniProt IDs
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    all_predictions = []
    all_plddt = []

    for uid in uniprot_ids:
        # 1) 予測取得
        pred = alphafold_get_prediction(uid)
        if pred:
            all_predictions.append(pred)

        # 2) pLDDT プロファイル
        plddt = alphafold_plddt_profile(uid, output_dir=output_dir)
        if not plddt.empty:
            plddt["uniprot_id"] = uid
            all_plddt.append(plddt)

    # サマリ
    pred_df = pd.DataFrame(all_predictions)
    pred_df.to_csv(output_dir / "predictions.csv", index=False)

    if all_plddt:
        plddt_df = pd.concat(all_plddt, ignore_index=True)
        plddt_df.to_csv(output_dir / "plddt_profiles.csv", index=False)

    print(f"AlphaFold pipeline: {len(all_predictions)} structures")
    return {"predictions": pred_df}
```

---

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|---------|
| `alphafold` | AlphaFold Database | 予測構造・pLDDT・PAE 取得 |

## パイプライン統合

```
protein-structure-analysis → alphafold-structures → protein-design
  (PDB 実験構造)            (AlphaFold 予測)       (de novo 設計)
       │                          │                     ↓
  structural-proteomics ─────────┘              molecular-docking
  (EMDB/PDBe)               │                  (結合予測)
                              ↓
                    variant-effect-prediction
                    (構造ベース変異評価)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/predictions.csv` | 予測メタデータ | → protein-structure-analysis |
| `results/plddt_profiles.csv` | 残基別 pLDDT | → protein-design |
| `results/AF-*.cif` | 予測構造ファイル | → molecular-docking |

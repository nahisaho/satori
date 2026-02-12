---
name: scientific-protein-structure-analysis
description: |
  タンパク質構造解析スキル。PDB / AlphaFold DB / PDBe を活用した 3D 構造解析。
  構造アラインメント、結合部位検出、分子ドッキング準備、構造品質評価。
  ToolUniverse の Protein Structure Retrieval と claude-scientific-skills の
  PDB/AlphaFold スキルを統合。
  「タンパク質の構造を解析して」「PDB 構造を調べて」「ドッキング準備して」で発火。
tu_tools:
  - key: proteinsplus
    name: ProteinsPlus
    description: タンパク質結合部位検出・構造解析ツール群
---

# Scientific Protein Structure Analysis

タンパク質 3D 構造の包括的解析スキル。PDB、AlphaFold DB、PDBe から
構造データを取得し、構造品質評価・アラインメント・結合部位解析・
ドッキング準備までを一貫して実行する。

## When to Use

- タンパク質の 3D 構造を取得・評価したいとき
- AlphaFold 予測構造の信頼度 (pLDDT) を確認したいとき
- 結合部位・ポケットを特定したいとき
- 構造アラインメントで複数構造を比較したいとき
- 分子ドッキングの受容体準備を行うとき
- リガンド結合様式を解析したいとき

## Quick Start

### 1. 構造解析パイプライン

```
Input: Target protein (UniProt ID / Gene Symbol / PDB ID)
    ↓
Step 1: Structure Search
  - PDB 実験構造検索（X-ray, Cryo-EM, NMR）
  - AlphaFold 予測構造取得
  - 最適構造の選択
    ↓
Step 2: Quality Assessment
  - 分解能 / R-factor / Rfree
  - pLDDT (AlphaFold)
  - Ramachandran 分析
    ↓
Step 3: Structural Features
  - ドメインアーキテクチャ
  - 結合部位・ポケット検出
  - 翻訳後修飾サイト
    ↓
Step 4: Functional Analysis
  - リガンド結合解析
  - タンパク質間相互作用インターフェース
  - アロステリックサイト
    ↓
Output: Structure Report + Prepared Files
```

---

## Phase 1: Structure Retrieval

### PDB 構造取得パイプライン

```python
import requests

def search_pdb_structures(uniprot_id):
    """
    UniProt ID から PDB 構造を検索。
    RCSB PDB Search API v2 を使用。
    """
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": uniprot_id,
            }
        },
        "return_type": "entry",
        "request_options": {
            "sort": [{"sort_by": "rcsb_accession_info.deposit_date", "direction": "desc"}],
            "results_content_type": ["experimental"],
        }
    }

    response = requests.post(
        "https://search.rcsb.org/rcsbsearch/v2/query",
        json=query
    )
    return response.json().get("result_set", [])


def get_alphafold_structure(uniprot_id):
    """
    AlphaFold DB から予測構造を取得。
    """
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()[0]
        return {
            "pdb_url": data["pdbUrl"],
            "cif_url": data["cifUrl"],
            "pae_url": data.get("paeImageUrl"),
            "model_version": data.get("latestVersion"),
        }
    return None
```

### 構造選択基準

```markdown
## Structure Selection Priority

1. **X-ray crystallography** (Resolution < 2.5 Å) — Gold standard
2. **Cryo-EM** (Resolution < 3.5 Å) — Large complexes
3. **NMR** — Solution state dynamics
4. **AlphaFold** (pLDDT > 70) — No experimental structure available

### Selection Decision Tree
- Has PDB structure?
  - Yes → Resolution < 2.5 Å? → Use X-ray
  - Yes → Large complex? → Use Cryo-EM
  - No → AlphaFold pLDDT > 70? → Use AlphaFold
  - No → Homology model needed
```

---

## Phase 2: Quality Assessment

### 構造品質メトリクス

```python
def assess_structure_quality(pdb_id):
    """
    PDB 構造の品質メトリクスを取得。
    """
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    response = requests.get(url)
    data = response.json()

    quality = {
        "pdb_id": pdb_id,
        "method": data.get("exptl", [{}])[0].get("method"),
        "resolution": data.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0],
        "r_factor": data.get("refine", [{}])[0].get("ls_R_factor_R_work"),
        "r_free": data.get("refine", [{}])[0].get("ls_R_factor_R_free"),
        "deposit_date": data.get("rcsb_accession_info", {}).get("deposit_date"),
    }

    # 品質判定
    res = quality.get("resolution")
    if res:
        if res < 2.0:
            quality["quality_tier"] = "Excellent"
        elif res < 2.5:
            quality["quality_tier"] = "Good"
        elif res < 3.0:
            quality["quality_tier"] = "Moderate"
        else:
            quality["quality_tier"] = "Low"

    return quality
```

### AlphaFold 信頼度評価

```python
def assess_alphafold_confidence(pdb_file_path):
    """
    AlphaFold pLDDT スコアの解析。
    B-factor 列に pLDDT が格納されている。
    """
    from Bio.PDB import PDBParser
    import numpy as np

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("af", pdb_file_path)

    plddt_scores = []
    for atom in structure.get_atoms():
        if atom.name == "CA":
            plddt_scores.append(atom.bfactor)

    plddt_array = np.array(plddt_scores)

    return {
        "mean_plddt": np.mean(plddt_array),
        "median_plddt": np.median(plddt_array),
        "pct_very_high": np.sum(plddt_array > 90) / len(plddt_array) * 100,  # >90: very high
        "pct_confident": np.sum(plddt_array > 70) / len(plddt_array) * 100,   # >70: confident
        "pct_low": np.sum(plddt_array < 50) / len(plddt_array) * 100,         # <50: low/disordered
        "interpretation": {
            ">90": "Very high confidence (well modeled)",
            "70-90": "Confident (backbone reliable)",
            "50-70": "Low confidence (caution)",
            "<50": "Very low (likely disordered / use with care)",
        }
    }
```

---

## Phase 3: Binding Site Analysis

### ポケット検出

```python
def detect_binding_sites(pdb_id):
    """
    PDBe Arpeggio / fpocket を用いた結合部位検出。
    """
    # fpocket ベースのポケット検出
    # fpocket は Voronoi テッセレーションで溶媒アクセス可能なポケットを検出
    pocket_info = {
        "method": "fpocket / DoGSiteScorer",
        "pockets": [],
    }

    # PDBe Binding Sites API
    url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/binding_sites/{pdb_id}"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json().get(pdb_id.lower(), [])
        for site in data:
            pocket_info["pockets"].append({
                "site_id": site.get("site_id"),
                "details": site.get("details"),
                "residues": site.get("site_residues", []),
            })

    return pocket_info
```

---

## Report Template

```markdown
# Protein Structure Report: [PROTEIN NAME]

**UniProt**: [accession] | **Date**: [date]

## 1. Structure Summary
| Feature | Value |
|---------|-------|
| Best PDB | |
| Resolution | |
| Method | |
| AlphaFold available | |
| Mean pLDDT | |

## 2. Available Structures
| PDB ID | Method | Resolution | Ligands | Chains |
|--------|--------|------------|---------|--------|

## 3. Quality Assessment
### 3.1 Selected Structure
### 3.2 Ramachandran Statistics
### 3.3 AlphaFold Confidence Map

## 4. Domain Architecture
| Domain | Start | End | Pfam | Description |
|--------|-------|-----|------|-------------|

## 5. Binding Sites
| Site | Key Residues | Known Ligands | Druggability |
|------|-------------|---------------|-------------|

## 6. Structural Insights
### 6.1 Active Site
### 6.2 Allosteric Sites
### 6.3 PPI Interface

## 7. Files Generated
- [ ] Cleaned PDB (waters removed, single chain)
- [ ] AlphaFold structure
- [ ] Pocket analysis results
```

---

## Completeness Checklist

- [ ] PDB 検索: UniProt ID で全実験構造を取得
- [ ] AlphaFold 構造: 予測構造の取得と pLDDT 評価
- [ ] 品質評価: Resolution / R-factor / pLDDT
- [ ] ドメイン: InterPro / Pfam でドメインマッピング
- [ ] 結合部位: 少なくとも 1 ポケットの詳細解析
- [ ] リガンド: 共結晶リガンドの一覧

## Best Practices

1. **実験構造を優先**: AlphaFold は実験構造がない場合の補完
2. **分解能 < 2.5 Å を選択**: ドッキング用途では高分解能が必須
3. **pLDDT < 50 領域は信用しない**: ディスオーダー領域の可能性大
4. **複数構造の比較**: 異なるコンフォメーションをキャプチャ
5. **リガンド情報を活用**: 共結晶リガンドから pharmacophore を推定

## References

### Output Files

| ファイル | 形式 | 生成タイミング |
|---|---|---|
| `results/structure_report.md` | 構造解析レポート（Markdown） | 全解析完了時 |
| `results/structure_analysis.json` | 構造データ（JSON） | 品質評価完了時 |
| `results/binding_sites.json` | 結合部位データ（JSON） | ポケット解析完了時 |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| UniProt | `UniProt_get_entry_by_accession` | タンパク質エントリ取得 |
| InterPro | `InterPro_get_protein_domains` | ドメインアノテーション |
| InterPro | `InterProScan_scan_sequence` | 配列ドメインスキャン |
| BindingDB | `BindingDB_get_ligands_by_uniprot` | リガンド結合データ |
| Proteins API | `proteins_api_get_features` | タンパク質特徴情報 |
| AlphaMissense | `AlphaMissense_get_protein_scores` | 残基レベル病原性予測 |

### 参照スキル

| スキル | 連携 |
|---|---|
| `scientific-drug-target-profiling` | ← ターゲットタンパク質の構造解析依頼 |
| `scientific-bioinformatics` | ← 配列・ドメイン情報 |
| `scientific-sequence-analysis` | ← アミノ酸配列解析結果 |
| `scientific-protein-design` | → 構造に基づく de novo タンパク質設計 |
| `scientific-cheminformatics` | → 結合部位情報を分子ドッキングに利用 |
| `scientific-academic-writing` | → 研究成果の論文化 |

---
name: scientific-structural-proteomics
description: |
  構造プロテオミクス統合スキル。EMDB クライオ EM、PDBe 構造データ、
  Proteins API (UniProt)、Complex Portal 複合体、DeepGO 機能予測、
  EVE 変異影響評価の統合パイプライン。
---

# Scientific Structural Proteomics

6 つの構造・機能データベースを統合した
タンパク質構造・機能プロファイリングパイプラインを提供する。

## When to Use

- PDB エントリの詳細構造情報 (分解能, 実験手法, 二次構造) を取得するとき
- クライオ EM マップ情報を EMDB から取得するとき
- UniProt タンパク質のドメイン・変異・ゲノムマッピングを調べるとき
- タンパク質複合体のサブユニット構成を Complex Portal で調べるとき
- DeepGO で配列ベースの GO 機能予測を行うとき
- EVE で missense 変異の病原性スコアを取得するとき

---

## Quick Start

## 1. PDBe 構造データ取得

```python
import requests
import pandas as pd

PDBE_API = "https://www.ebi.ac.uk/pdbe/api"


def get_pdbe_entry(pdb_id):
    """
    PDBe エントリのサマリー・品質・二次構造を取得。

    Parameters:
        pdb_id: str — PDB ID (e.g., "1cbs")

    ToolUniverse:
        pdbe_get_entry_summary(pdb_id=pdb_id)
        pdbe_get_entry_quality(pdb_id=pdb_id)
        pdbe_get_entry_experiment(pdb_id=pdb_id)
        pdbe_get_entry_molecules(pdb_id=pdb_id)
        pdbe_get_entry_secondary_structure(pdb_id=pdb_id)
        pdbe_get_entry_assemblies(pdb_id=pdb_id)
        pdbe_get_entry_publications(pdb_id=pdb_id)
    """
    pdb = pdb_id.lower()

    # Summary
    resp_s = requests.get(f"{PDBE_API}/pdb/entry/summary/{pdb}")
    resp_s.raise_for_status()
    summary = resp_s.json().get(pdb, [{}])[0]

    # Quality
    resp_q = requests.get(f"{PDBE_API}/pdb/entry/quality/{pdb}")
    quality = (
        resp_q.json().get(pdb, [{}])[0] if resp_q.status_code == 200 else {}
    )

    # Experiment
    resp_e = requests.get(f"{PDBE_API}/pdb/entry/experiment/{pdb}")
    experiment = (
        resp_e.json().get(pdb, [{}])[0] if resp_e.status_code == 200 else {}
    )

    entry = {
        "pdb_id": pdb,
        "title": summary.get("title", ""),
        "method": experiment.get("experimental_method", ""),
        "resolution": experiment.get("resolution", ""),
        "release_date": summary.get("release_date", ""),
        "r_factor": quality.get("r_factor", ""),
        "r_free": quality.get("r_free", ""),
    }

    print(f"PDBe {pdb}: {entry['method']} @ {entry['resolution']}Å")
    return entry, summary, quality
```

## 2. EMDB クライオ EM 構造

```python
EMDB_API = "https://www.ebi.ac.uk/emdb/api"


def get_emdb_structure(emdb_id):
    """
    EMDB エントリの 3D EM マップ情報を取得。

    Parameters:
        emdb_id: str — EMDB ID (e.g., "EMD-1234")

    ToolUniverse:
        EMDB_get_structure(emdb_id=emdb_id)
        EMDB_get_map_info(emdb_id=emdb_id)
        EMDB_get_sample_info(emdb_id=emdb_id)
        EMDB_get_imaging_info(emdb_id=emdb_id)
        EMDB_get_validation(emdb_id=emdb_id)
        EMDB_get_publications(emdb_id=emdb_id)
        EMDB_search_structures(query=query)
    """
    url = f"{EMDB_API}/entry/{emdb_id}"
    resp = requests.get(url)
    resp.raise_for_status()
    data = resp.json()

    admin = data.get("admin", {})
    map_info = data.get("map", {})

    entry = {
        "emdb_id": emdb_id,
        "title": admin.get("title", ""),
        "resolution": map_info.get("resolution", {}).get("value", ""),
        "contour_level": map_info.get("contour_level", {}).get("value", ""),
        "deposition_date": admin.get("deposition_date", ""),
    }

    print(f"EMDB {emdb_id}: {entry['resolution']}Å resolution")
    return entry, data
```

## 3. Proteins API (UniProt) ドメイン・変異

```python
PROTEINS_API = "https://www.ebi.ac.uk/proteins/api"


def get_protein_features(uniprot_id):
    """
    Proteins API からタンパク質のドメイン・変異・修飾情報を取得。

    Parameters:
        uniprot_id: str — UniProt accession (e.g., "P04637")

    ToolUniverse:
        proteins_api_get_protein(accession=uniprot_id)
        proteins_api_get_features(accession=uniprot_id)
        proteins_api_get_variants(accession=uniprot_id)
        proteins_api_get_xrefs(accession=uniprot_id)
        proteins_api_get_genome_mappings(accession=uniprot_id)
    """
    headers = {"Accept": "application/json"}

    # Protein info
    resp_p = requests.get(
        f"{PROTEINS_API}/proteins/{uniprot_id}", headers=headers
    )
    resp_p.raise_for_status()
    protein = resp_p.json()

    # Features
    resp_f = requests.get(
        f"{PROTEINS_API}/features/{uniprot_id}", headers=headers
    )
    features = resp_f.json() if resp_f.status_code == 200 else {}

    feature_list = features.get("features", [])
    domains = [f for f in feature_list if f.get("type") == "DOMAIN"]
    variants = [f for f in feature_list if f.get("type") == "VARIANT"]

    print(f"Proteins API {uniprot_id}: "
          f"{len(domains)} domains, {len(variants)} variants")
    return protein, feature_list
```

## 4. Complex Portal 複合体

```python
def get_complex_portal(complex_id):
    """
    Complex Portal からタンパク質複合体の詳細を取得。

    Parameters:
        complex_id: str — Complex Portal ID (e.g., "CPX-1")

    ToolUniverse:
        ComplexPortal_get_complex(complex_id=complex_id)
        ComplexPortal_search_complexes(query=query)
    """
    url = (
        f"https://www.ebi.ac.uk/intact/complex-ws/complex/{complex_id}"
    )
    resp = requests.get(url, headers={"Accept": "application/json"})
    resp.raise_for_status()
    data = resp.json()

    participants = data.get("participants", [])
    subunits = []
    for p in participants:
        interactor = p.get("interactor", {})
        subunits.append({
            "identifier": interactor.get("identifier", ""),
            "name": interactor.get("name", ""),
            "stoichiometry": p.get("stoichiometry", ""),
        })

    print(f"Complex Portal {complex_id}: "
          f"{data.get('name', '?')}, {len(subunits)} subunits")
    return data, subunits
```

## 5. DeepGO 機能予測 & EVE 変異評価

```python
def predict_deepgo_function(sequence):
    """
    DeepGO でタンパク質配列から GO 機能予測。

    ToolUniverse:
        DeepGO_predict_function(sequence=sequence)
    """
    url = "https://deepgo.cbrc.kaust.edu.sa/deepgo/api/create"
    resp = requests.post(url, json={"data_format": "fasta", "data": sequence})
    resp.raise_for_status()
    data = resp.json()

    predictions = data.get("predictions", [])
    print(f"DeepGO: {len(predictions)} GO term predictions")
    return predictions


def get_eve_variant_score(gene_name, variant=None):
    """
    EVE (Evolutionary model of Variant Effect) による
    missense 変異の病原性スコアを取得。

    ToolUniverse:
        EVE_get_gene_info(gene_name=gene_name)
        EVE_get_variant_score(gene_name=gene_name, variant=variant)
    """
    url = f"https://evemodel.org/api/proteins/{gene_name}"

    if variant:
        url += f"/variants/{variant}"

    resp = requests.get(url)
    resp.raise_for_status()
    data = resp.json()

    if variant:
        score = data.get("eve_score", None)
        classification = data.get("eve_class", "")
        print(f"EVE {gene_name} {variant}: score={score}, class={classification}")
    else:
        print(f"EVE {gene_name}: found")

    return data
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/pdbe_entry.json` | JSON |
| `results/emdb_structure.json` | JSON |
| `results/protein_features.json` | JSON |
| `results/complex_portal.json` | JSON |
| `results/deepgo_predictions.json` | JSON |
| `results/eve_scores.json` | JSON |

### 利用可能ツール

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| PDBe | `pdbe_get_entry_summary` | エントリサマリー |
| PDBe | `pdbe_get_entry_quality` | 構造品質 |
| PDBe | `pdbe_get_entry_experiment` | 実験手法 |
| PDBe | `pdbe_get_entry_molecules` | 分子情報 |
| PDBe | `pdbe_get_entry_secondary_structure` | 二次構造 |
| PDBe | `pdbe_get_entry_assemblies` | 生物学的集合体 |
| PDBe | `pdbe_get_entry_publications` | 文献 |
| PDBe | `pdbe_get_entry_related_publications` | 関連文献 |
| PDBe | `pdbe_get_entry_observed_residues_ratio` | 観察残基 |
| PDBe | `pdbe_get_entry_status` | ステータス |
| EMDB | `EMDB_get_structure` | EM 構造 |
| EMDB | `EMDB_get_map_info` | マップ情報 |
| EMDB | `EMDB_get_sample_info` | サンプル情報 |
| EMDB | `EMDB_get_imaging_info` | イメージング |
| EMDB | `EMDB_get_validation` | バリデーション |
| EMDB | `EMDB_get_publications` | 文献 |
| EMDB | `EMDB_search_structures` | EM 構造検索 |
| Proteins API | `proteins_api_get_protein` | タンパク質情報 |
| Proteins API | `proteins_api_get_features` | ドメイン・特徴 |
| Proteins API | `proteins_api_get_variants` | 変異 |
| Proteins API | `proteins_api_get_xrefs` | 外部参照 |
| Proteins API | `proteins_api_get_genome_mappings` | ゲノムマッピング |
| Proteins API | `proteins_api_get_epitopes` | エピトープ |
| Proteins API | `proteins_api_get_proteomics` | プロテオミクス |
| Proteins API | `proteins_api_get_publications` | 文献 |
| Proteins API | `proteins_api_get_comments` | コメント |
| Proteins API | `proteins_api_search` | 検索 |
| Complex Portal | `ComplexPortal_get_complex` | 複合体詳細 |
| Complex Portal | `ComplexPortal_search_complexes` | 複合体検索 |
| DeepGO | `DeepGO_predict_function` | GO 機能予測 |
| EVE | `EVE_get_gene_info` | 遺伝子情報 |
| EVE | `EVE_get_variant_score` | 変異スコア |

### 参照スキル

| スキル | 関連 |
|---|---|
| `scientific-protein-structure-analysis` | タンパク質構造 |
| `scientific-proteomics-mass-spectrometry` | プロテオミクス |
| `scientific-protein-interaction-network` | PPI |
| `scientific-molecular-docking` | ドッキング |
| `scientific-variant-interpretation` | 変異解釈 |

### 依存パッケージ

`requests`, `pandas`

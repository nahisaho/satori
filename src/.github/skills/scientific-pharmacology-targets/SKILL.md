---
name: scientific-pharmacology-targets
description: |
  薬理学的ターゲットプロファイリングスキル。BindingDB 結合親和性、
  GPCRdb GPCR 構造-活性、GtoPdb 薬理学、BRENDA 酵素動態、
  Pharos 未解明ターゲット(TDL)の統合解析パイプライン。
---

# Scientific Pharmacology Targets

複数の薬理学データベース (BindingDB, GPCRdb, GtoPdb, BRENDA, Pharos) を
統合した包括的ターゲットプロファイリングパイプラインを提供する。

## When to Use

- 特定タンパク質の既知リガンド・結合親和性を調べるとき
- GPCR のリガンド・変異・構造情報を取得するとき
- 薬物-ターゲット相互作用のデータベース横断検索を行うとき
- 酵素阻害剤データ (BRENDA) を調べるとき
- 未解明ターゲット (Tdark/Tbio) のドラッガビリティを評価するとき

---

## Quick Start

## 1. BindingDB 結合親和性データ取得

```python
import requests
import pandas as pd


def get_bindingdb_ligands(uniprot_id, cutoff=None):
    """
    BindingDB から UniProt ID ベースのリガンド結合データを取得。

    Parameters:
        uniprot_id: str — UniProt accession (e.g., "P00533")
        cutoff: float | None — affinity cutoff nM

    ToolUniverse:
        BindingDB_get_ligands_by_uniprot(uniprot_id=uniprot_id)
        BindingDB_get_targets_by_compound(smiles=smiles)
    """
    url = "https://bindingdb.org/axis2/services/BDBService"
    params = {
        "uniprot": uniprot_id,
        "response": "json",
    }
    if cutoff:
        params["cutoff"] = cutoff

    resp = requests.get(f"{url}/getLigandsByUniprot", params=params)
    resp.raise_for_status()
    data = resp.json()

    ligands = data.get("getLigandsByUniprotResponse", {}).get("affinities", [])
    results = []
    for lig in ligands:
        results.append({
            "monomer_id": lig.get("monomerid", ""),
            "smiles": lig.get("smiles", ""),
            "affinity_type": lig.get("affinity_type", ""),
            "affinity_value_nm": lig.get("affinity", ""),
            "source": lig.get("source", ""),
        })

    df = pd.DataFrame(results)
    print(f"BindingDB '{uniprot_id}': {len(df)} ligands")
    return df
```

## 2. GPCRdb GPCR プロファイリング

```python
def get_gpcrdb_profile(protein_entry):
    """
    GPCRdb から GPCR のリガンド・変異・構造情報を取得。

    Parameters:
        protein_entry: str — GPCRdb entry name (e.g., "adrb2_human")

    ToolUniverse:
        GPCRdb_get_protein(entry_name=protein_entry)
        GPCRdb_get_ligands(entry_name=protein_entry)
        GPCRdb_get_mutations(entry_name=protein_entry)
        GPCRdb_get_structures(entry_name=protein_entry)
        GPCRdb_list_proteins()
    """
    base = "https://gpcrdb.org/services"

    # Protein info
    resp_p = requests.get(f"{base}/protein/{protein_entry}/")
    resp_p.raise_for_status()
    protein = resp_p.json()

    # Ligands
    resp_l = requests.get(f"{base}/ligands/{protein_entry}/")
    ligands = resp_l.json() if resp_l.status_code == 200 else []

    # Mutations
    resp_m = requests.get(f"{base}/mutants/{protein_entry}/")
    mutations = resp_m.json() if resp_m.status_code == 200 else []

    # Structures
    resp_s = requests.get(f"{base}/structure/protein/{protein_entry}/")
    structures = resp_s.json() if resp_s.status_code == 200 else []

    profile = {
        "entry_name": protein.get("entry_name", ""),
        "name": protein.get("name", ""),
        "family": protein.get("family", ""),
        "species": protein.get("species", ""),
        "num_ligands": len(ligands),
        "num_mutations": len(mutations),
        "num_structures": len(structures),
    }

    print(f"GPCRdb '{protein_entry}': "
          f"{profile['num_ligands']} ligands, "
          f"{profile['num_mutations']} mutations, "
          f"{profile['num_structures']} structures")
    return profile, ligands, mutations, structures
```

## 3. GtoPdb 薬理学データ

```python
def get_gtopdb_target_pharmacology(target_id):
    """
    Guide to PHARMACOLOGY (GtoPdb) から
    ターゲットの薬理学的相互作用データを取得。

    ToolUniverse:
        GtoPdb_get_target(target_id=target_id)
        GtoPdb_get_target_interactions(target_id=target_id)
        GtoPdb_get_ligand(ligand_id=ligand_id)
        GtoPdb_search_interactions(query=query)
    """
    base = "https://www.guidetopharmacology.org/services"

    # Target info
    resp_t = requests.get(f"{base}/targets/{target_id}")
    resp_t.raise_for_status()
    target = resp_t.json()

    # Interactions
    resp_i = requests.get(f"{base}/targets/{target_id}/interactions")
    interactions = resp_i.json() if resp_i.status_code == 200 else []

    results = []
    for ix in interactions:
        results.append({
            "ligand_id": ix.get("ligandId", ""),
            "ligand_name": ix.get("ligandName", ""),
            "type": ix.get("type", ""),
            "action": ix.get("action", ""),
            "affinity_type": ix.get("affinityType", ""),
            "affinity_median": ix.get("affinityMedian", ""),
            "approved": ix.get("approvedDrug", False),
        })

    df = pd.DataFrame(results)
    print(f"GtoPdb target {target_id} ({target.get('name', '')}): "
          f"{len(df)} interactions")
    return target, df
```

## 4. Pharos/TCRD 未解明ターゲット検索

```python
def search_pharos_targets(query, tdl=None):
    """
    Pharos / TCRD からターゲット情報を取得。
    Target Development Level (TDL) でフィルタ可能。

    Parameters:
        query: str — gene symbol or target name
        tdl: str | None — "Tclin", "Tchem", "Tbio", "Tdark"

    ToolUniverse:
        Pharos_search_targets(q=query)
        Pharos_get_target(q=query)
        Pharos_get_tdl_summary()
        Pharos_get_disease_targets(disease_name=disease_name)
    """
    url = "https://pharos-api.ncats.io/graphql"
    gql = """
    query TargetSearch($term: String!, $top: Int) {
        targets(filter: { term: $term }, top: $top) {
            targets {
                name
                sym
                uniprot { accession }
                tdl
                fam
                novelty
                jensenScore
                diseaseCount
                ligandCount
            }
            count
        }
    }
    """
    variables = {"term": query, "top": 20}
    resp = requests.post(url, json={"query": gql, "variables": variables})
    resp.raise_for_status()
    data = resp.json()["data"]["targets"]

    results = []
    for t in data["targets"]:
        results.append({
            "symbol": t.get("sym", ""),
            "name": t.get("name", ""),
            "uniprot": t.get("uniprot", {}).get("accession", ""),
            "tdl": t.get("tdl", ""),
            "family": t.get("fam", ""),
            "novelty": t.get("novelty", 0),
            "disease_count": t.get("diseaseCount", 0),
            "ligand_count": t.get("ligandCount", 0),
        })

    df = pd.DataFrame(results)
    if tdl:
        df = df[df["tdl"] == tdl]

    print(f"Pharos '{query}': {len(df)} targets"
          f"{f' (TDL={tdl})' if tdl else ''}")
    return df
```

## 5. 統合ターゲットプロファイリング

```python
def integrated_target_profile(uniprot_id, gene_symbol):
    """
    複数データベースを統合したターゲットプロファイル。

    ToolUniverse (横断):
        BindingDB_get_ligands_by_uniprot(uniprot_id)
        Pharos_get_target(q=gene_symbol)
        GtoPdb_get_targets() → GtoPdb_get_target_interactions()
    """
    profile = {
        "uniprot_id": uniprot_id,
        "gene_symbol": gene_symbol,
    }

    # BindingDB ligands
    try:
        bdb_df = get_bindingdb_ligands(uniprot_id)
        profile["bindingdb_ligand_count"] = len(bdb_df)
    except Exception:
        profile["bindingdb_ligand_count"] = 0

    # Pharos TDL
    try:
        pharos_df = search_pharos_targets(gene_symbol)
        if not pharos_df.empty:
            row = pharos_df.iloc[0]
            profile["tdl"] = row.get("tdl", "")
            profile["novelty"] = row.get("novelty", 0)
    except Exception:
        profile["tdl"] = "Unknown"

    print(f"Integrated profile {gene_symbol}: TDL={profile.get('tdl', '?')}, "
          f"BindingDB={profile['bindingdb_ligand_count']} ligands")
    return profile
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/bindingdb_ligands.csv` | CSV |
| `results/gpcrdb_profile.json` | JSON |
| `results/gtopdb_interactions.csv` | CSV |
| `results/pharos_targets.csv` | CSV |
| `results/integrated_target_profile.json` | JSON |

### 利用可能ツール

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| BindingDB | `BindingDB_get_ligands_by_uniprot` | UniProt→リガンド |
| BindingDB | `BindingDB_get_ligands_by_uniprots` | バッチ |
| BindingDB | `BindingDB_get_ligands_by_pdb` | PDB→リガンド |
| BindingDB | `BindingDB_get_targets_by_compound` | 化合物→ターゲット |
| GPCRdb | `GPCRdb_get_protein` | GPCR 詳細 |
| GPCRdb | `GPCRdb_get_ligands` | GPCR リガンド |
| GPCRdb | `GPCRdb_get_mutations` | GPCR 変異 |
| GPCRdb | `GPCRdb_get_structures` | GPCR 構造 |
| GPCRdb | `GPCRdb_list_proteins` | GPCR 一覧 |
| GtoPdb | `GtoPdb_get_target` | ターゲット情報 |
| GtoPdb | `GtoPdb_get_target_interactions` | 相互作用 |
| GtoPdb | `GtoPdb_get_ligand` | リガンド情報 |
| GtoPdb | `GtoPdb_get_targets` | ターゲット一覧 |
| GtoPdb | `GtoPdb_list_ligands` | リガンド一覧 |
| GtoPdb | `GtoPdb_get_disease` | 疾患関連 |
| GtoPdb | `GtoPdb_list_diseases` | 疾患一覧 |
| GtoPdb | `GtoPdb_search_interactions` | 相互作用検索 |
| BRENDA | `BRENDA_get_inhibitors` | 酵素阻害剤 |
| Pharos | `Pharos_search_targets` | ターゲット検索 |
| Pharos | `Pharos_get_target` | ターゲット詳細 |
| Pharos | `Pharos_get_tdl_summary` | TDL サマリー |
| Pharos | `Pharos_get_disease_targets` | 疾患→ターゲット |

### 参照スキル

| スキル | 関連 |
|---|---|
| `scientific-drug-target-interaction` | DTI 予測 |
| `scientific-compound-similarity` | 化合物類似性 |
| `scientific-compound-screening` | 化合物スクリーニング |
| `scientific-molecular-docking` | 分子ドッキング |
| `scientific-protein-interaction-network` | PPI ネットワーク |

### 依存パッケージ

`requests`, `pandas`

---
name: scientific-plant-biology
description: |
  植物バイオロジー統合スキル。Plant Reactome 代謝パスウェイ・
  TAIR Arabidopsis ゲノム情報・Phytozome 比較ゲノミクス・
  Ensembl Plants 種間オーソログ解析。
tu_tools:
  - key: tair
    name: TAIR
    description: シロイヌナズナゲノム・植物データ検索
---

# Scientific Plant Biology

Plant Reactome / TAIR / Phytozome / Ensembl Plants を活用した
植物ゲノム・代謝パスウェイ統合解析パイプラインを提供する。

## When to Use

- 植物代謝パスウェイ解析 (Plant Reactome) を実行するとき
- Arabidopsis thaliana の遺伝子・タンパク質情報を取得するとき
- 植物種間の比較ゲノミクス解析を行うとき
- 植物オーソログ・パラログを同定するとき
- 作物改良のための候補遺伝子を探索するとき
- 植物表現型データと遺伝子型を統合するとき

---

## Quick Start

## 1. Plant Reactome パスウェイ検索

```python
import requests
import pandas as pd
import json

PLANT_REACTOME = "https://plantreactome.gramene.org/ContentService"


def plant_reactome_search(query, species="Oryza sativa"):
    """
    Plant Reactome — 植物代謝/シグナルパスウェイ検索。

    Parameters:
        query: str — 検索クエリ (例: "photosynthesis")
        species: str — 種名
    """
    url = f"{PLANT_REACTOME}/search/query"
    params = {"query": query, "species": species, "cluster": True}
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for group in data.get("results", []):
        for entry in group.get("entries", []):
            results.append({
                "stId": entry.get("stId", ""),
                "name": entry.get("name", ""),
                "species": entry.get("species", ""),
                "type": entry.get("exactType", ""),
                "compartment": entry.get("compartmentNames", []),
            })

    df = pd.DataFrame(results)
    print(f"Plant Reactome: '{query}' → {len(df)} entries ({species})")
    return df


def plant_reactome_pathway_detail(pathway_id):
    """
    Plant Reactome パスウェイ詳細取得。

    Parameters:
        pathway_id: str — パスウェイ ID (例: "R-OSA-1119616")
    """
    url = f"{PLANT_REACTOME}/data/pathway/{pathway_id}/containedEvents"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    events = resp.json()

    steps = []
    for event in events:
        steps.append({
            "stId": event.get("stId", ""),
            "name": event.get("displayName", ""),
            "type": event.get("className", ""),
            "input_count": len(event.get("input", [])),
            "output_count": len(event.get("output", [])),
            "catalyst": event.get("catalystActivity", [{}])[0].get(
                "displayName", "") if event.get("catalystActivity") else "",
        })

    df = pd.DataFrame(steps)
    print(f"Pathway {pathway_id}: {len(df)} reaction steps")
    return df
```

## 2. TAIR Arabidopsis 遺伝子情報

```python
TAIR_BASE = "https://www.arabidopsis.org/api"


def tair_gene_search(gene_id=None, gene_name=None, keyword=None):
    """
    TAIR — Arabidopsis thaliana 遺伝子情報取得。

    Parameters:
        gene_id: str — AGI ID (例: "AT1G01010")
        gene_name: str — 遺伝子名 (例: "FLC")
        keyword: str — キーワード検索
    """
    if gene_id:
        url = f"{TAIR_BASE}/gene/{gene_id}"
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        return pd.DataFrame([{
            "agi_id": data.get("locus", ""),
            "name": data.get("name", ""),
            "description": data.get("description", ""),
            "chromosome": data.get("chromosome", ""),
            "start": data.get("start", ""),
            "end": data.get("end", ""),
            "strand": data.get("strand", ""),
            "gene_model_type": data.get("gene_model_type", ""),
        }])

    # キーワード検索
    search_term = gene_name or keyword or ""
    url = f"{TAIR_BASE}/search/gene"
    params = {"query": search_term, "limit": 50}
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for gene in data.get("results", []):
        results.append({
            "agi_id": gene.get("locus", ""),
            "name": gene.get("name", ""),
            "description": gene.get("description", ""),
            "chromosome": gene.get("chromosome", ""),
        })

    df = pd.DataFrame(results)
    print(f"TAIR: '{search_term}' → {len(df)} genes")
    return df


def tair_gene_expression(gene_id):
    """
    TAIR — 遺伝子発現パターン取得。

    Parameters:
        gene_id: str — AGI ID
    """
    url = f"{TAIR_BASE}/gene/{gene_id}/expression"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    tissues = []
    for expr in data.get("expression", []):
        tissues.append({
            "tissue": expr.get("tissue", ""),
            "stage": expr.get("developmental_stage", ""),
            "level": expr.get("expression_level", ""),
            "source": expr.get("source", ""),
        })

    df = pd.DataFrame(tissues)
    print(f"TAIR expression: {gene_id} → {len(df)} tissue records")
    return df
```

## 3. Ensembl Plants 種間比較

```python
ENSEMBL_PLANTS = "https://rest.ensembl.org"


def ensembl_plants_orthologs(gene_id, source_species="arabidopsis_thaliana",
                              target_species=None):
    """
    Ensembl Plants — 植物種間オーソログ検索。

    Parameters:
        gene_id: str — Ensembl Gene ID or symbol
        source_species: str — 起源種
        target_species: str — ターゲット種 (None = 全種)
    """
    url = f"{ENSEMBL_PLANTS}/homology/id/{gene_id}"
    params = {
        "type": "orthologues",
        "content-type": "application/json",
        "compara": "plants",
    }
    if target_species:
        params["target_species"] = target_species

    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    orthologs = []
    for homology in data.get("data", [{}])[0].get("homologies", []):
        target = homology.get("target", {})
        orthologs.append({
            "source_gene": gene_id,
            "source_species": source_species,
            "target_gene": target.get("id", ""),
            "target_species": target.get("species", ""),
            "target_protein": target.get("protein_id", ""),
            "identity": target.get("perc_id", 0),
            "dn_ds": homology.get("dn_ds", None),
            "type": homology.get("type", ""),
        })

    df = pd.DataFrame(orthologs)
    print(f"Ensembl Plants orthologs: {gene_id} → {len(df)} homologs")
    return df
```

## 4. Phytozome 比較ゲノミクス

```python
PHYTOZOME_BASE = "https://phytozome-next.jgi.doe.gov/api"


def phytozome_gene_family(gene_id, species="Athaliana"):
    """
    Phytozome — 遺伝子ファミリー・比較ゲノミクス。

    Parameters:
        gene_id: str — 遺伝子 ID
        species: str — 種略称
    """
    url = f"{PHYTOZOME_BASE}/search"
    params = {"query": gene_id, "organism": species}
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    families = []
    for hit in data.get("hits", []):
        families.append({
            "gene_id": hit.get("gene_id", ""),
            "family_id": hit.get("family_id", ""),
            "family_name": hit.get("family_name", ""),
            "species": hit.get("organism", ""),
            "annotation": hit.get("annotation", ""),
            "pfam_domains": hit.get("pfam", []),
        })

    df = pd.DataFrame(families)
    print(f"Phytozome: {gene_id} → {len(df)} family members")
    return df
```

## 5. 植物バイオロジー統合パイプライン

```python
def plant_biology_pipeline(gene_query, species="Oryza sativa",
                            output_dir="results"):
    """
    植物バイオロジー統合パイプライン。

    Parameters:
        gene_query: str — 遺伝子/パスウェイクエリ
        species: str — 対象種
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) Plant Reactome パスウェイ
    pathways = plant_reactome_search(gene_query, species=species)
    pathways.to_csv(output_dir / "plant_pathways.csv", index=False)

    # 2) TAIR (Arabidopsis ならば)
    tair_genes = tair_gene_search(keyword=gene_query)
    tair_genes.to_csv(output_dir / "tair_genes.csv", index=False)

    # 3) Ensembl Plants オーソログ
    if len(tair_genes) > 0:
        top_gene = tair_genes.iloc[0]["agi_id"]
        orthologs = ensembl_plants_orthologs(top_gene)
        orthologs.to_csv(output_dir / "orthologs.csv", index=False)
    else:
        orthologs = pd.DataFrame()

    print(f"Plant biology pipeline: {output_dir}")
    return {
        "pathways": pathways,
        "tair_genes": tair_genes,
        "orthologs": orthologs,
    }
```

---

## パイプライン統合

```
pathway-enrichment → plant-biology → environmental-ecology
  (KEGG/Reactome)   (PlantReactome)   (生態学/環境)
       │                   │                ↓
  gene-annotation ────────┘         marine-ecology
  (GO/InterPro)      │              (OBIS/WoRMS)
                      ↓
               comparative-genomics
               (Ensembl 比較)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/plant_pathways.csv` | Plant Reactome パスウェイ | → pathway-enrichment |
| `results/tair_genes.csv` | TAIR Arabidopsis 遺伝子 | → gene-annotation |
| `results/orthologs.csv` | 種間オーソログ | → comparative-genomics |

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `tair` | TAIR | シロイヌナズナゲノム・植物データ検索 |

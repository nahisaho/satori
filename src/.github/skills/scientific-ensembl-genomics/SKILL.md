---
name: scientific-ensembl-genomics
description: |
  Ensembl REST API ゲノミクススキル。遺伝子ルックアップ・配列取得・
  VEP (Variant Effect Predictor) バリアントアノテーション・
  クロスリファレンス・制御要素・系統樹・相同性検索・分類学統合パイプライン。
---

# Scientific Ensembl Genomics

Ensembl REST API (rest.ensembl.org) を活用したゲノミクスデータアクセス
パイプラインを提供する。遺伝子情報取得、VEP バリアント効果予測、
相同性検索、系統解析を統合。

## When to Use

- Ensembl Gene ID から遺伝子情報・座標を取得するとき
- VEP でバリアントの機能的影響を予測するとき (SIFT/PolyPhen/CADD)
- 遺伝子のオルソログ・パラログを検索するとき
- Ensembl ↔ UniProt / RefSeq / HGNC 間の ID 変換をするとき
- ゲノム領域の制御要素 (promoter/enhancer) を検索するとき
- 種間比較ゲノミクスデータを取得するとき

---

## Quick Start

## 1. 遺伝子ルックアップ

```python
import requests
import pandas as pd

ENSEMBL_REST = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}


def lookup_gene(gene_id, expand=True):
    """
    Ensembl 遺伝子情報取得。

    Parameters:
        gene_id: str — Ensembl Gene ID (例: "ENSG00000141510")
        expand: bool — トランスクリプト情報を含めるか

    ToolUniverse:
        ensembl_lookup_gene(gene_id=gene_id, species="homo_sapiens")
    """
    url = f"{ENSEMBL_REST}/lookup/id/{gene_id}"
    params = {"expand": 1 if expand else 0}
    resp = requests.get(url, params=params, headers=HEADERS)
    resp.raise_for_status()
    data = resp.json()

    info = {
        "id": data.get("id"),
        "display_name": data.get("display_name"),
        "biotype": data.get("biotype"),
        "species": data.get("species"),
        "assembly_name": data.get("assembly_name"),
        "seq_region_name": data.get("seq_region_name"),
        "start": data.get("start"),
        "end": data.get("end"),
        "strand": data.get("strand"),
        "description": data.get("description"),
    }

    if expand and "Transcript" in data:
        info["n_transcripts"] = len(data["Transcript"])
        info["canonical_transcript"] = data.get("canonical_transcript")

    print(f"Gene: {info['display_name']} ({info['id']}), "
          f"chr{info['seq_region_name']}:{info['start']}-{info['end']}")
    return info
```

## 2. 配列取得

```python
def get_sequence(seq_id, seq_type="genomic", species="homo_sapiens"):
    """
    Ensembl 配列取得 (DNA/cDNA/CDS/protein)。

    Parameters:
        seq_id: str — Ensembl ID (Gene/Transcript/Translation)
        seq_type: str — "genomic", "cdna", "cds", "protein"
        species: str — 生物種

    ToolUniverse:
        ensembl_get_sequence(id=seq_id, type=seq_type, species=species)
    """
    url = f"{ENSEMBL_REST}/sequence/id/{seq_id}"
    params = {"type": seq_type}
    resp = requests.get(url, params=params, headers=HEADERS)
    resp.raise_for_status()
    data = resp.json()

    result = {
        "id": data.get("id"),
        "seq_type": seq_type,
        "molecule": data.get("molecule"),
        "length": len(data.get("seq", "")),
        "sequence": data.get("seq"),
    }

    print(f"Sequence: {result['id']} ({seq_type}), {result['length']} bp/aa")
    return result
```

## 3. VEP (Variant Effect Predictor)

```python
def vep_region(species, chromosome, position, allele,
               sift=True, polyphen=True, cadd=False):
    """
    VEP によるバリアント効果予測。

    Parameters:
        species: str — "homo_sapiens"
        chromosome: str — 染色体番号
        position: int — ゲノム座標
        allele: str — 代替アレル (例: "T")
        sift: bool — SIFT 予測を含める
        polyphen: bool — PolyPhen 予測を含める
        cadd: bool — CADD スコアを含める

    ToolUniverse:
        ensembl_vep_region(
            species=species, region=f"{chromosome}:{position}:{position}",
            allele=allele, SIFT="b", PolyPhen="b"
        )
    """
    region = f"{chromosome}:{position}:{position}"
    url = f"{ENSEMBL_REST}/vep/{species}/region/{region}/{allele}"
    params = {}
    if sift:
        params["SIFT"] = "b"
    if polyphen:
        params["PolyPhen"] = "b"
    if cadd:
        params["CADD"] = 1

    resp = requests.get(url, params=params, headers=HEADERS)
    resp.raise_for_status()
    results = resp.json()

    consequences = []
    for r in results:
        for tc in r.get("transcript_consequences", []):
            cons = {
                "gene_symbol": tc.get("gene_symbol"),
                "gene_id": tc.get("gene_id"),
                "transcript_id": tc.get("transcript_id"),
                "consequence_terms": tc.get("consequence_terms", []),
                "impact": tc.get("impact"),
                "biotype": tc.get("biotype"),
                "amino_acids": tc.get("amino_acids"),
                "codons": tc.get("codons"),
            }
            if "sift_prediction" in tc:
                cons["sift"] = f"{tc['sift_prediction']}({tc.get('sift_score')})"
            if "polyphen_prediction" in tc:
                cons["polyphen"] = f"{tc['polyphen_prediction']}({tc.get('polyphen_score')})"
            consequences.append(cons)

    df = pd.DataFrame(consequences)
    print(f"VEP {chromosome}:{position} {allele}: "
          f"{len(df)} transcript consequences")
    return df
```

## 4. クロスリファレンス (ID 変換)

```python
def get_xrefs(ensembl_id, external_db=None):
    """
    Ensembl ID から外部 DB の ID を取得。

    Parameters:
        ensembl_id: str — Ensembl ID
        external_db: str — フィルタ DB 名 (例: "UniProt", "RefSeq", "HGNC")

    ToolUniverse:
        ensembl_get_xrefs(id=ensembl_id, external_db=external_db)
    """
    url = f"{ENSEMBL_REST}/xrefs/id/{ensembl_id}"
    params = {}
    if external_db:
        params["external_db"] = external_db

    resp = requests.get(url, params=params, headers=HEADERS)
    resp.raise_for_status()
    xrefs = resp.json()

    rows = []
    for x in xrefs:
        rows.append({
            "primary_id": x.get("primary_id"),
            "display_id": x.get("display_id"),
            "dbname": x.get("dbname"),
            "description": x.get("description", "")[:100],
        })

    df = pd.DataFrame(rows)
    print(f"Cross-references for {ensembl_id}: {len(df)} entries")
    return df
```

## 5. 相同性検索 (オルソログ/パラログ)

```python
def get_homology(species, gene_symbol, target_species=None,
                 homology_type="orthologues"):
    """
    遺伝子相同性検索。

    Parameters:
        species: str — ソース生物種
        gene_symbol: str — 遺伝子シンボル
        target_species: str — ターゲット生物種 (None で全種)
        homology_type: str — "orthologues", "paralogues", "all"

    ToolUniverse:
        ensembl_get_homology(
            species=species, symbol=gene_symbol,
            target_species=target_species, type=homology_type
        )
    """
    url = f"{ENSEMBL_REST}/homology/symbol/{species}/{gene_symbol}"
    params = {"type": homology_type}
    if target_species:
        params["target_species"] = target_species

    resp = requests.get(url, params=params, headers=HEADERS)
    resp.raise_for_status()

    homologies = resp.json().get("data", [{}])[0].get("homologies", [])

    rows = []
    for h in homologies:
        target = h.get("target", {})
        rows.append({
            "type": h.get("type"),
            "target_species": target.get("species"),
            "target_gene_id": target.get("id"),
            "target_symbol": target.get("protein_id"),
            "perc_id": target.get("perc_id"),
            "perc_pos": target.get("perc_pos"),
            "dn_ds": h.get("dn_ds"),
        })

    df = pd.DataFrame(rows)
    print(f"Homologs of {gene_symbol} ({species}): {len(df)} found")
    return df
```

## 6. 制御要素検索

```python
def get_regulatory_features(species, region):
    """
    ゲノム領域の制御要素 (promoter/enhancer/CTCF) 検索。

    Parameters:
        species: str — 生物種 (例: "homo_sapiens")
        region: str — ゲノム領域 (例: "7:140000000-140100000")

    ToolUniverse:
        ensembl_get_regulatory_features(region=region, species=species)
        ensembl_get_overlap_features(region=region)
    """
    url = f"{ENSEMBL_REST}/overlap/region/{species}/{region}"
    params = {"feature": "regulatory"}

    resp = requests.get(url, params=params, headers=HEADERS)
    resp.raise_for_status()
    features = resp.json()

    rows = []
    for f in features:
        rows.append({
            "id": f.get("id"),
            "feature_type": f.get("feature_type"),
            "start": f.get("start"),
            "end": f.get("end"),
            "strand": f.get("strand"),
            "description": f.get("description", ""),
        })

    df = pd.DataFrame(rows)
    print(f"Regulatory features in {region}: {len(df)}")
    return df
```

## 7. 遺伝子系統樹

```python
def get_gene_tree(gene_id, prune_species=None):
    """
    遺伝子ファミリーの系統樹取得。

    Parameters:
        gene_id: str — Ensembl Gene ID
        prune_species: list — 系統樹を制限する生物種リスト

    ToolUniverse:
        ensembl_get_genetree(id=gene_id, prune_species=species_list)
    """
    url = f"{ENSEMBL_REST}/genetree/member/id/{gene_id}"
    params = {"sequence": "none", "aligned": 0}
    if prune_species:
        params["prune_species"] = ";".join(prune_species)

    resp = requests.get(url, params=params, headers=HEADERS)
    resp.raise_for_status()
    tree = resp.json()

    result = {
        "tree_id": tree.get("tree", {}).get("id"),
        "type": tree.get("tree", {}).get("type"),
        "n_members": _count_leaves(tree.get("tree", {})),
    }

    print(f"Gene tree {result['tree_id']}: {result['n_members']} members")
    return tree


def _count_leaves(node):
    """系統樹リーフ数をカウント。"""
    if "children" not in node:
        return 1
    return sum(_count_leaves(c) for c in node["children"])
```

---

## パイプライン統合

```
bioinformatics ───→ ensembl-genomics ───→ variant-interpretation
  (Ensembl Gene ID)   (VEP アノテーション)    (ACMG/AMP 分類)
        │                    │                      ↓
genome-sequence-tools ──┘    │              variant-effect-prediction
  (BLAST/dbSNP)              │              (AlphaMissense/CADD)
                             ↓
                    regulatory-genomics → epigenomics-chromatin
                    (RegulomeDB/ReMap)    (ChIP-seq/ATAC-seq)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/ensembl_gene_info.json` | 遺伝子情報 | → bioinformatics |
| `results/vep_consequences.csv` | VEP バリアント効果 | → variant-interpretation |
| `results/homology_table.csv` | オルソログ/パラログ | → phylogenetics |
| `results/regulatory_features.csv` | 制御要素 | → regulatory-genomics |

## 利用可能ツール (ToolUniverse SMCP)

| ツール名 | 用途 |
|---------|------|
| `ensembl_lookup_gene` | 遺伝子ルックアップ |
| `ensembl_get_sequence` | 配列取得 |
| `ensembl_get_variants` | バリアント取得 |
| `ensembl_get_variation` | バリエーション詳細 |
| `ensembl_get_variation_phenotypes` | バリアント表現型 |
| `ensembl_vep_region` | VEP 効果予測 |
| `ensembl_get_xrefs` | クロスリファレンス |
| `ensembl_get_xrefs_by_name` | 名前ベース xref |
| `ensembl_get_regulatory_features` | 制御要素 |
| `ensembl_get_genetree` | 遺伝子系統樹 |
| `ensembl_get_homology` | 相同性検索 |
| `ensembl_get_alignment` | 配列アラインメント |
| `ensembl_get_taxonomy` | 分類学情報 |
| `ensembl_get_species` | 生物種一覧 |
| `ensembl_get_ontology_term` | GO オントロジー |
| `ensembl_get_overlap_features` | 領域オーバーラップ |

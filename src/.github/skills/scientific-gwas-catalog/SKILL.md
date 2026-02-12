---
name: scientific-gwas-catalog
description: |
  GWAS カタログスキル。NHGRI-EBI GWAS Catalog REST API によるゲノム
  ワイド関連研究メタデータ・関連シグナル・形質・遺伝子座検索。
  ToolUniverse 連携: gwas。
tu_tools:
  - key: gwas
    name: GWAS Catalog
    description: GWAS 関連シグナル・形質・遺伝子座検索
---

# Scientific GWAS Catalog

NHGRI-EBI GWAS Catalog REST API を活用した GWAS メタデータ
解析・遺伝子座レベル解釈パイプラインを提供する。

## When to Use

- GWAS Catalog から疾患/形質の関連バリアントを検索するとき
- 遺伝的関連シグナルのエフェクトサイズ・P値を取得するとき
- 特定遺伝子座の LD ブロック情報を解析するとき
- 多形質 PheWAS-like 解析を実施するとき
- GWAS サマリ統計量を下流解析に準備するとき
- 公開 GWAS データから PRS ウェイトを抽出するとき

---

## Quick Start

## 1. GWAS 関連シグナル検索

```python
import requests
import pandas as pd
import numpy as np

GWAS_BASE = "https://www.ebi.ac.uk/gwas/rest/api"


def gwas_search_associations(trait=None, gene=None, variant=None,
                               p_upper=5e-8, limit=100):
    """
    GWAS Catalog — 関連シグナル検索。

    Parameters:
        trait: str — 形質/疾患 EFO ID or 名前 (例: "EFO_0001645")
        gene: str — 遺伝子名 (例: "BRCA1")
        variant: str — rsID (例: "rs1234567")
        p_upper: float — P値上限
        limit: int — 最大結果数
    """
    if trait:
        url = f"{GWAS_BASE}/efoTraits/{trait}/associations"
    elif gene:
        url = f"{GWAS_BASE}/associations/search/findByGene"
    elif variant:
        url = f"{GWAS_BASE}/singleNucleotidePolymorphisms/{variant}/associations"
    else:
        url = f"{GWAS_BASE}/associations"

    params = {"size": limit}
    if gene:
        params["geneName"] = gene

    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    associations = data.get("_embedded", {}).get("associations", [])
    results = []
    for assoc in associations:
        p_value = assoc.get("pvalue", 1.0)
        if p_value and float(p_value) > p_upper:
            continue

        loci = assoc.get("loci", [{}])
        genes = []
        for locus in loci:
            for gene_info in locus.get("authorReportedGenes", []):
                genes.append(gene_info.get("geneName", ""))

        snps = []
        for snp_info in assoc.get("snps", []):
            snps.append(snp_info.get("rsId", ""))

        results.append({
            "association_id": assoc.get("associationId", ""),
            "p_value": float(p_value) if p_value else None,
            "p_value_mlog": assoc.get("pvalueMantissa", 0),
            "or_beta": assoc.get("orPerCopyNum", None),
            "beta_num": assoc.get("betaNum", None),
            "beta_direction": assoc.get("betaDirection", ""),
            "ci": assoc.get("range", ""),
            "risk_allele_freq": assoc.get("riskFrequency", ""),
            "snps": "; ".join(snps),
            "genes": "; ".join(genes),
            "trait": assoc.get("efoTraits", [{}])[0].get("trait", "")
                     if assoc.get("efoTraits") else "",
            "study_accession": assoc.get("study", {}).get(
                "accessionId", ""),
        })

    df = pd.DataFrame(results)
    print(f"GWAS associations: {len(df)} results "
          f"(trait={trait}, gene={gene}, p<{p_upper})")
    return df.sort_values("p_value") if not df.empty else df
```

## 2. GWAS 研究メタデータ検索

```python
def gwas_search_studies(query=None, efo_trait=None, limit=50):
    """
    GWAS Catalog — 研究メタデータ検索。

    Parameters:
        query: str — フリーテキスト検索
        efo_trait: str — EFO 形質 ID
        limit: int — 最大結果数
    """
    if efo_trait:
        url = f"{GWAS_BASE}/efoTraits/{efo_trait}/studies"
    else:
        url = f"{GWAS_BASE}/studies/search/findByDiseaseTrait"

    params = {"size": limit}
    if query:
        params["diseaseTrait"] = query

    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    studies = data.get("_embedded", {}).get("studies", [])
    results = []
    for s in studies:
        results.append({
            "accession": s.get("accessionId", ""),
            "title": s.get("title", ""),
            "pubmed_id": s.get("publicationInfo", {}).get(
                "pubmedId", ""),
            "author": s.get("publicationInfo", {}).get(
                "author", {}).get("fullname", ""),
            "journal": s.get("publicationInfo", {}).get(
                "publication", ""),
            "date": s.get("publicationInfo", {}).get(
                "publicationDate", ""),
            "initial_sample_size": s.get("initialSampleSize", ""),
            "replication_sample_size": s.get(
                "replicationSampleSize", ""),
            "ancestry": s.get("ancestries", []),
        })

    df = pd.DataFrame(results)
    print(f"GWAS studies: {len(df)} results")
    return df
```

## 3. GWAS 形質検索・PheWAS

```python
def gwas_phewas(variant_rsid, p_threshold=5e-8):
    """
    GWAS Catalog — バリアント PheWAS (形質横断検索)。

    Parameters:
        variant_rsid: str — rsID (例: "rs7903146")
        p_threshold: float — P値閾値
    """
    url = (f"{GWAS_BASE}/singleNucleotidePolymorphisms/"
           f"{variant_rsid}/associations")
    resp = requests.get(url, params={"size": 500}, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    associations = data.get("_embedded", {}).get("associations", [])
    results = []
    for assoc in associations:
        p_val = assoc.get("pvalue", 1.0)
        if p_val and float(p_val) > p_threshold:
            continue
        for trait in assoc.get("efoTraits", []):
            results.append({
                "variant": variant_rsid,
                "trait": trait.get("trait", ""),
                "efo_uri": trait.get("shortForm", ""),
                "p_value": float(p_val) if p_val else None,
                "or_beta": assoc.get("orPerCopyNum", None),
                "study": assoc.get("study", {}).get(
                    "accessionId", ""),
            })

    df = pd.DataFrame(results)
    if not df.empty:
        df = df.sort_values("p_value")
    print(f"PheWAS {variant_rsid}: {len(df)} trait associations")
    return df
```

## 4. GWAS 統合パイプライン

```python
def gwas_catalog_pipeline(trait_query, output_dir="results"):
    """
    GWAS Catalog 統合パイプライン。

    Parameters:
        trait_query: str — 形質/疾患名
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) 研究検索
    studies = gwas_search_studies(query=trait_query)
    studies.to_csv(output_dir / "gwas_studies.csv", index=False)

    # 2) 関連シグナル
    assocs = gwas_search_associations(gene=None, trait=None)
    assocs.to_csv(output_dir / "gwas_associations.csv", index=False)

    # 3) トップバリアントの PheWAS
    if not assocs.empty:
        top_snps = assocs["snps"].str.split("; ").explode().unique()[:5]
        phewas_all = []
        for rsid in top_snps:
            if rsid.startswith("rs"):
                phewas = gwas_phewas(rsid)
                phewas_all.append(phewas)
        if phewas_all:
            phewas_df = pd.concat(phewas_all, ignore_index=True)
            phewas_df.to_csv(output_dir / "phewas.csv", index=False)

    print(f"GWAS pipeline: {output_dir}")
    return {"studies": studies, "associations": assocs}
```

---

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|---------|
| `gwas` | GWAS Catalog | 関連シグナル・形質・研究メタデータ検索 |

## パイプライン統合

```
disease-research → gwas-catalog → variant-interpretation
  (DisGeNET/OMIM)  (GWAS Catalog)  (ACMG/AMP)
       │                 │                ↓
  population-genetics ──┘         variant-effect-prediction
  (Fst/PCA)           │          (CADD/SpliceAI)
                       ↓
                  precision-oncology
                  (臨床的意義判定)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/gwas_studies.csv` | GWAS 研究メタデータ | → literature-search |
| `results/gwas_associations.csv` | 関連シグナル | → variant-interpretation |
| `results/phewas.csv` | PheWAS 結果 | → disease-research |

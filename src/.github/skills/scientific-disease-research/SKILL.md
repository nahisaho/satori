---
name: scientific-disease-research
description: |
  疾患研究スキル。GWAS Catalog・Orphanet・OMIM・HPO・DisGeNET を統合し、
  疾患-遺伝子関連解析・希少疾患診断支援・表現型-遺伝型マッピング・
  疫学的特性評価を支援。
  「疾患と遺伝子の関連を調べて」「希少疾患を診断して」「GWAS 結果を解析して」で発火。
tu_tools:
  - key: disgenet
    name: DisGeNET
    description: 疾患-遺伝子関連スコア (GDA) データベース
---

# Scientific Disease Research

疾患研究のための統合解析スキル。ゲノムワイド関連解析（GWAS）、
希少疾患データベース、表現型オントロジーを横断活用し、
疾患メカニズムの解明と診断支援を行う。

## When to Use

- GWAS 結果の解釈と遺伝子-疾患関連の同定
- 希少疾患の表現型マッチング・診断候補絞り込み
- 遺伝的リスクスコア (PRS) の算出
- 疾患ネットワーク・共有メカニズム解析
- HPO (Human Phenotype Ontology) ベースの表現型解析
- DisGeNET での遺伝子-疾患関連スコアリング

## Quick Start

### 疾患研究パイプライン

```
Phase 1: Phenotype Characterization
  - HPO ベースの表現型プロファイリング
  - OMIM / Orphanet 疾患同定
  - 表現型クラスタリング
    ↓
Phase 2: Genetic Association
  - GWAS Catalog 検索 (EBI)
  - Significant loci 同定 (p < 5e-8)
  - LD (連鎖不平衡) ブロック解析
    ↓
Phase 3: Gene-Disease Mapping
  - DisGeNET GDA スコア算出
  - OMIM Morbid Map 参照
  - eQTL / sQTL 機能的アノテーション
    ↓
Phase 4: Rare Disease Diagnosis
  - Orphanet 疾患カタログ参照
  - HPO-based differential diagnosis
  - ACMG variant classification
    ↓
Phase 5: Network & Pathway Analysis
  - 疾患モジュール同定 (PPI ネットワーク上)
  - パスウェイ濃縮解析
  - 疾患間コモビディティ解析
    ↓
Phase 6: Report Generation
  - 疾患研究レポート (JSON + Markdown)
  - 遺伝子-疾患関連テーブル
  - 診断候補リスト
```

## Workflow

### 1. GWAS Catalog API

```python
import requests
import pandas as pd
import numpy as np

GWAS_API = "https://www.ebi.ac.uk/gwas/rest/api"

def search_gwas_associations(trait_keyword, p_threshold=5e-8):
    """GWAS Catalog で形質関連の SNP を検索"""
    resp = requests.get(
        f"{GWAS_API}/efoTraits/search/findBySearchQuery",
        params={"searchString": trait_keyword}
    )
    traits = resp.json().get("_embedded", {}).get("efoTraits", [])

    all_assocs = []
    for trait in traits[:5]:  # 上位 5 形質
        trait_uri = trait["_links"]["self"]["href"]
        assoc_resp = requests.get(
            f"{trait_uri}/associations",
            params={"size": 500}
        )
        for a in assoc_resp.json().get("_embedded", {}).get("associations", []):
            p_value = float(a.get("pvalue", 1))
            if p_value < p_threshold:
                for locus in a.get("loci", []):
                    for gene in locus.get("authorReportedGenes", []):
                        all_assocs.append({
                            "trait": trait.get("trait", ""),
                            "gene": gene.get("geneName", ""),
                            "rsid": a.get("snpInteraction", False),
                            "p_value": p_value,
                            "or_beta": a.get("orPerCopyNum", ""),
                            "risk_allele": "",
                            "study": a.get("study", {}).get("publicationInfo", {}).get("title", ""),
                        })

    df = pd.DataFrame(all_assocs).sort_values("p_value")
    print(f"GWAS associations for '{trait_keyword}': {len(df)}")
    return df

gwas_results = search_gwas_associations("type 2 diabetes")
```

### 2. DisGeNET 遺伝子-疾患関連

```python
DISGENET_API = "https://www.disgenet.org/api"
DISGENET_KEY = "YOUR_API_KEY"

def query_disgenet_gda(gene_symbol, source="ALL"):
    """DisGeNET で Gene-Disease Associations (GDA) を取得"""
    headers = {"Authorization": f"Bearer {DISGENET_KEY}"}
    resp = requests.get(
        f"{DISGENET_API}/gda/gene/{gene_symbol}",
        headers=headers,
        params={"source": source, "format": "json"}
    )
    data = resp.json()

    results = []
    for item in data:
        results.append({
            "gene": item.get("gene_symbol", ""),
            "disease": item.get("disease_name", ""),
            "disease_id": item.get("diseaseid", ""),
            "score": item.get("score", 0),
            "ei": item.get("ei", 0),  # Evidence Index
            "el": item.get("el", ""),  # Evidence Level
            "n_pmids": item.get("pmid_count", 0),
            "source": item.get("source", ""),
        })

    df = pd.DataFrame(results).sort_values("score", ascending=False)
    return df

# GDA Score 解釈:
# 0.0-0.3: Weak association
# 0.3-0.6: Moderate association
# 0.6-0.8: Strong association
# 0.8-1.0: Very strong / curated association
```

### 3. HPO (Human Phenotype Ontology) 表現型解析

```python
HPO_API = "https://hpo.jax.org/api"

def phenotype_matching(hpo_terms):
    """HPO タームリストから疾患候補をランキング"""

    # HPO term → disease mapping
    disease_scores = {}

    for hpo_id in hpo_terms:
        resp = requests.get(f"{HPO_API}/hpo/term/{hpo_id}/diseases")
        diseases = resp.json().get("diseases", [])

        for d in diseases:
            disease_name = d.get("diseaseName", "")
            disease_id = d.get("diseaseId", "")
            if disease_id not in disease_scores:
                disease_scores[disease_id] = {
                    "name": disease_name,
                    "matched_hpo": [],
                    "total_hpo": d.get("numberOfAnnotations", 0),
                }
            disease_scores[disease_id]["matched_hpo"].append(hpo_id)

    # Jaccard-like スコア
    results = []
    for did, info in disease_scores.items():
        matched = len(info["matched_hpo"])
        total_query = len(hpo_terms)
        total_disease = info["total_hpo"]
        # Harmonic mean of precision and recall
        precision = matched / total_query if total_query > 0 else 0
        recall = matched / total_disease if total_disease > 0 else 0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

        results.append({
            "disease_id": did,
            "disease_name": info["name"],
            "matched_hpo_count": matched,
            "total_disease_hpo": total_disease,
            "precision": round(precision, 3),
            "recall": round(recall, 3),
            "f1_score": round(f1, 3),
        })

    return pd.DataFrame(results).sort_values("f1_score", ascending=False)

# 例: 表現型ベースの診断
patient_hpo = ["HP:0001250", "HP:0001249", "HP:0000252"]  # seizures, ID, microcephaly
candidates = phenotype_matching(patient_hpo)
print("Top diagnostic candidates:")
print(candidates.head(10))
```

### 4. Orphanet 希少疾患検索

```python
ORPHANET_API = "https://api.orphacode.org"

def search_orphanet(query):
    """Orphanet で希少疾患を検索"""
    resp = requests.get(
        f"{ORPHANET_API}/EN/ClinicalEntity/ApproximateName/{query}",
        headers={"apiKey": "YOUR_ORPHANET_KEY"}
    )
    data = resp.json()

    results = []
    for item in data:
        results.append({
            "orpha_code": item.get("ORPHAcode", ""),
            "name": item.get("Preferred term", ""),
            "type": item.get("typology", ""),
            "prevalence": item.get("prevalence", ""),
            "inheritance": item.get("inheritance", ""),
            "age_of_onset": item.get("age_of_onset", ""),
        })

    return pd.DataFrame(results)


def get_disease_genes(orpha_code):
    """Orphanet 疾患に関連する遺伝子を取得"""
    resp = requests.get(
        f"{ORPHANET_API}/EN/ClinicalEntity/orphacode/{orpha_code}/Gene",
        headers={"apiKey": "YOUR_ORPHANET_KEY"}
    )
    data = resp.json()
    genes = []
    for g in data.get("DisorderGeneAssociationList", []):
        gene_info = g.get("Gene", {})
        genes.append({
            "gene_symbol": gene_info.get("Symbol", ""),
            "gene_name": gene_info.get("Name", ""),
            "association_type": g.get("DisorderGeneAssociationType", {}).get("Name", ""),
            "status": g.get("DisorderGeneAssociationStatus", {}).get("Name", ""),
        })
    return pd.DataFrame(genes)
```

### 5. Polygenic Risk Score (PRS) 算出

```python
def calculate_prs(gwas_summary_stats, individual_genotypes):
    """
    Polygenic Risk Score (PRS) by C+T (Clumping + Thresholding)

    PRS = Σ (beta_i × genotype_i) for i in selected SNPs
    """
    # LD Clumping (概念的)
    clumped = gwas_summary_stats[gwas_summary_stats["p_value"] < 5e-8].copy()

    # 効果量の方向統一
    clumped["beta"] = np.where(
        clumped["effect_allele"] == clumped["risk_allele"],
        clumped["beta"],
        -clumped["beta"]
    )

    # PRS 算出
    prs_scores = []
    for sample_id, geno in individual_genotypes.groupby("sample_id"):
        merged = clumped.merge(geno, on="rsid", how="inner")
        prs = (merged["beta"] * merged["dosage"]).sum()
        prs_scores.append({
            "sample_id": sample_id,
            "prs_raw": prs,
            "n_snps_used": len(merged),
        })

    prs_df = pd.DataFrame(prs_scores)

    # Z-score 標準化
    prs_df["prs_zscore"] = (prs_df["prs_raw"] - prs_df["prs_raw"].mean()) / prs_df["prs_raw"].std()

    # パーセンタイル
    prs_df["percentile"] = prs_df["prs_zscore"].rank(pct=True) * 100

    return prs_df
```

### 6. 疾患研究レポート生成

```python
import json

def generate_disease_report(disease_name, gwas_df, gda_df, hpo_df,
                            output_dir="results"):
    """疾患研究統合レポート"""
    report = {
        "disease": disease_name,
        "analysis_date": pd.Timestamp.now().isoformat(),
        "gwas_summary": {
            "total_associations": len(gwas_df),
            "genome_wide_significant": len(gwas_df[gwas_df["p_value"] < 5e-8]),
            "top_loci": gwas_df.nsmallest(10, "p_value").to_dict("records"),
        },
        "gene_disease_associations": {
            "total_genes": gda_df["gene"].nunique(),
            "strong_associations": len(gda_df[gda_df["score"] >= 0.6]),
            "top_genes": gda_df.nlargest(10, "score").to_dict("records"),
        },
        "phenotype_network": {
            "hpo_terms_used": len(hpo_df) if hpo_df is not None else 0,
        },
    }

    with open(f"{output_dir}/disease_research_report.json", "w") as f:
        json.dump(report, f, indent=2, default=str)

    md = f"# Disease Research Report: {disease_name}\n\n"
    md += f"## GWAS Summary\n\n"
    md += f"- Genome-wide significant loci: {report['gwas_summary']['genome_wide_significant']}\n\n"
    md += "| Gene | p-value | OR/Beta | Study |\n|---|---|---|---|\n"
    for locus in report["gwas_summary"]["top_loci"]:
        md += f"| {locus.get('gene', '')} | {locus.get('p_value', '')} | {locus.get('or_beta', '')} | {locus.get('study', '')[:50]} |\n"
    md += f"\n## Gene-Disease Associations (DisGeNET)\n\n"
    md += f"- Strong associations (score≥0.6): {report['gene_disease_associations']['strong_associations']}\n\n"
    md += "| Gene | Disease | GDA Score | Evidence |\n|---|---|---|---|\n"
    for g in report["gene_disease_associations"]["top_genes"]:
        md += f"| {g.get('gene', '')} | {g.get('disease', '')} | {g.get('score', '')} | {g.get('n_pmids', '')} PMIDs |\n"

    with open(f"{output_dir}/disease_research_report.md", "w") as f:
        f.write(md)

    return report
```

---

## Best Practices

1. **GWAS 有意水準**: ゲノムワイド有意水準は p < 5×10⁻⁸ を使用
2. **LD を考慮**: 隣接 SNP の連鎖不平衡を考慮し独立シグナルを同定
3. **DisGeNET スコア**: GDA score 0.6 以上を「強い関連」とする
4. **HPO の粒度**: 一般的すぎる HPO ターム（例: HP:0000001）は避ける
5. **PRS の限界**: PRS は集団レベルの予測であり個人の診断には限界がある
6. **希少疾患では遺伝子パネル優先**: GWAS より WES/WGS + ACMG 分類が有効
7. **交絡の考慮**: 集団層別化 (Population Stratification) を補正

## Completeness Checklist

- [ ] 疾患/形質の定義と HPO マッピング
- [ ] GWAS Catalog 検索・有意 loci 同定
- [ ] DisGeNET GDA スコア取得
- [ ] Orphanet 希少疾患情報参照
- [ ] HPO ベースの表現型マッチング
- [ ] パスウェイ・ネットワーク解析
- [ ] 疾患研究レポート（JSON + Markdown）生成

## References

### Output Files

| ファイル | 形式 | 生成タイミング |
|---|---|---|
| `results/disease_research_report.json` | 疾患研究レポート（JSON） | 解析完了時 |
| `results/disease_research_report.md` | 疾患研究レポート（Markdown） | レポート生成時 |
| `results/gwas_significant_loci.json` | GWAS 有意 loci（JSON） | GWAS 検索完了時 |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| OpenTargets | `OpenTargets_get_disease_id_description_by_name` | 疾患 ID・記述取得 |
| OpenTargets | `OpenTargets_get_associated_targets_by_disease_efoId` | 疾患関連ターゲット |
| EFO/OLS | `OSL_get_efo_id_by_disease_name` | EFO ID 解決 |
| HPO | `get_HPO_ID_by_phenotype` | 表現型→HPO マッピング |
| Monarch | `Monarch_get_gene_diseases` | 遺伝子-疾患関連 |
| ClinVar | `clinvar_search_variants` | 病原性バリアント検索 |
| ClinicalTrials | `search_clinical_trials` | 疾患関連臨床試験 |

### 参照スキル

| スキル | 連携 |
|---|---|
| `scientific-variant-interpretation` | → 候補バリアントの ACMG 分類 |
| `scientific-bioinformatics` | ← 発現データ・eQTL 解析 |
| `scientific-network-analysis` | ← 疾患モジュール・PPI ネットワーク |
| `scientific-meta-analysis` | ← GWAS メタアナリシス |
| `scientific-precision-oncology` | → がん疾患の精密医療連携 |
| `scientific-deep-research` | ← 疾患関連文献リサーチ |

---
name: scientific-cancer-genomics
description: |
  がんゲノミクスポータル統合スキル。COSMIC (体細胞変異カタログ)、
  cBioPortal (がんゲノミクスデータ解析)、DepMap (がん細胞依存性) の
  3 大がんゲノミクスデータベースを統合した変異プロファイリング、
  変異シグネチャー解析、遺伝子依存性 (essentiality) 評価、
  コピー数変化・がん種横断解析パイプライン。
  13 の ToolUniverse SMCP ツールと連携。
tu_tools:
  - key: cosmic
    name: COSMIC
    description: がん体細胞変異カタログ
  - key: cbioportal
    name: cBioPortal
    description: がんゲノミクスポータル
---

# Scientific Cancer Genomics

COSMIC / cBioPortal / DepMap の 3 大がんゲノミクスポータルを統合した
体細胞変異プロファイリング・機能解析パイプラインを提供する。

## When to Use

- がん関連遺伝子の体細胞変異をカタログ検索するとき
- cBioPortal でがん種横断の遺伝子変異頻度を調べるとき
- DepMap で遺伝子依存性 (essentiality) を評価するとき
- 変異シグネチャー解析 (SBS/DBS/ID) を行うとき
- コピー数変化 (CNA) のドライバー・パッセンジャー分類が必要なとき

---

## Quick Start

## 1. COSMIC 体細胞変異検索

```python
import pandas as pd
import numpy as np
import requests


def cosmic_search_mutations(gene, cancer_type=None, mutation_type=None):
    """
    COSMIC (Catalogue Of Somatic Mutations In Cancer) 変異検索。

    Parameters:
        gene: str — 遺伝子シンボル (e.g., "BRAF", "TP53")
        cancer_type: str — がん種フィルタ (e.g., "melanoma")
        mutation_type: str — 変異タイプ ("missense", "nonsense", "frameshift")
    """
    # ToolUniverse 経由: COSMIC_search_mutations, COSMIC_get_mutations_by_gene
    # COSMIC API は認証が必要 (Academic 無料)

    # Cancer Gene Census (CGC) チェック
    cgc_genes = {
        "TP53": {"role": "TSG", "tier": 1},
        "BRAF": {"role": "oncogene", "tier": 1},
        "KRAS": {"role": "oncogene", "tier": 1},
        "EGFR": {"role": "oncogene", "tier": 1},
        "PIK3CA": {"role": "oncogene", "tier": 1},
        "BRCA1": {"role": "TSG", "tier": 1},
        "BRCA2": {"role": "TSG", "tier": 1},
        "ALK": {"role": "oncogene", "tier": 1},
    }

    gene_info = cgc_genes.get(gene.upper(), {})

    result = {
        "gene": gene,
        "cgc_role": gene_info.get("role", "unknown"),
        "cgc_tier": gene_info.get("tier", None),
        "cancer_type_filter": cancer_type,
        "mutation_type_filter": mutation_type,
    }

    print(f"COSMIC query: {gene} "
          f"(CGC: {gene_info.get('role', 'N/A')}, "
          f"Tier {gene_info.get('tier', 'N/A')})")
    return result
```

## 2. cBioPortal がんゲノミクスデータ取得

```python
def cbioportal_query(genes, study_id=None, cancer_type=None,
                      data_types=None):
    """
    cBioPortal REST API によるがんゲノミクスデータ取得。

    Parameters:
        genes: list — 遺伝子シンボルリスト
        study_id: str — cBioPortal 研究 ID (e.g., "tcga_brca_pan_can_atlas_2018")
        cancer_type: str — がん種 (e.g., "Breast Cancer")
        data_types: list — ["mutations", "cna", "mrna", "methylation"]
    """
    base_url = "https://www.cbioportal.org/api"

    if data_types is None:
        data_types = ["mutations", "cna"]

    results = {}

    # 研究一覧取得
    if study_id is None:
        resp = requests.get(f"{base_url}/studies")
        studies = resp.json()
        if cancer_type:
            studies = [s for s in studies
                       if cancer_type.lower() in
                       s.get("cancerType", {}).get("name", "").lower()]
        print(f"cBioPortal: {len(studies)} studies for '{cancer_type}'")
        results["studies"] = pd.DataFrame([{
            "study_id": s["studyId"],
            "name": s["name"],
            "cancer_type": s.get("cancerType", {}).get("name", ""),
            "sample_count": s.get("allSampleCount", 0),
        } for s in studies[:20]])
    else:
        # 変異データ取得
        if "mutations" in data_types:
            url = f"{base_url}/molecular-profiles/{study_id}_mutations/mutations"
            params = {"projection": "DETAILED"}
            resp = requests.get(url, params=params)
            if resp.status_code == 200:
                mutations = resp.json()
                mut_df = pd.DataFrame([{
                    "gene": m.get("gene", {}).get("hugoGeneSymbol", ""),
                    "mutation": m.get("proteinChange", ""),
                    "mutation_type": m.get("mutationType", ""),
                    "chromosome": m.get("chr", ""),
                    "position": m.get("startPosition", ""),
                    "allele_freq": m.get("tumorAltCount", 0) /
                                   max(m.get("tumorRefCount", 1) +
                                       m.get("tumorAltCount", 1), 1),
                } for m in mutations
                    if m.get("gene", {}).get("hugoGeneSymbol", "") in genes])
                results["mutations"] = mut_df
                print(f"  Mutations: {len(mut_df)} found in {genes}")

    return results
```

## 3. DepMap 遺伝子依存性解析

```python
def depmap_gene_dependency(genes, cell_lineage=None):
    """
    DepMap (Cancer Dependency Map) 遺伝子依存性解析。

    Parameters:
        genes: list — 遺伝子シンボルリスト
        cell_lineage: str — 細胞系統フィルタ (e.g., "Lung", "Breast")
    """
    # ToolUniverse 経由:
    # DepMap_search_genes, DepMap_get_gene_dependencies
    # DepMap_get_cell_line, DepMap_get_cell_lines, DepMap_search_cell_lines

    # DepMap CRISPR (Chronos) dependency score:
    # negative = essential (依存), ~0 = non-essential
    # Common Essential: mean < -0.5 across 90% of lines
    # Selective Dependency: mean < -0.5 in specific lineages

    results = []
    for gene in genes:
        result = {
            "gene": gene,
            "cell_lineage": cell_lineage,
            "query_type": "CRISPR_dependency",
            # 実際のスコアは ToolUniverse 経由で取得
            "interpretation": (
                "Chronos score < 0: gene essentiality increases. "
                "score < -0.5: likely essential in this lineage. "
                "score ~ 0: non-essential."
            ),
        }
        results.append(result)

    df = pd.DataFrame(results)
    print(f"DepMap: queried {len(genes)} genes "
          f"(lineage: {cell_lineage or 'pan-cancer'})")
    return df
```

## 4. 変異シグネチャー解析

```python
def mutational_signature_analysis(mutations_df, genome="GRCh38",
                                    n_signatures=None):
    """
    体細胞変異シグネチャー解析 (COSMIC SBS signatures)。

    Parameters:
        mutations_df: DataFrame — columns: [chr, pos, ref, alt, sample]
        genome: str — 参照ゲノム
        n_signatures: int — 抽出シグネチャー数 (None=自動推定)
    """
    from itertools import product

    # 96 トリヌクレオチドコンテキスト
    bases = ["C", "T"]
    contexts = []
    for ref in bases:
        for alt in ["A", "C", "G", "T"]:
            if ref == alt:
                continue
            for five in "ACGT":
                for three in "ACGT":
                    contexts.append(f"{five}[{ref}>{alt}]{three}")

    # サンプルごとのカタログ構築
    samples = mutations_df["sample"].unique()
    catalog = pd.DataFrame(0, index=contexts, columns=samples)

    for _, row in mutations_df.iterrows():
        ref = row["ref"]
        alt = row["alt"]
        sample = row["sample"]
        context = row.get("trinucleotide_context", "N[N>N]N")
        if context in catalog.index:
            catalog.loc[context, sample] += 1

    print(f"Mutation catalog: {len(contexts)} contexts, "
          f"{len(samples)} samples, "
          f"{catalog.sum().sum():.0f} total mutations")

    # NMF 分解 (SigProfilerExtractor 代替)
    from sklearn.decomposition import NMF

    X = catalog.values.T  # samples × contexts
    if n_signatures is None:
        n_signatures = min(5, len(samples))

    model = NMF(n_components=n_signatures, random_state=42, max_iter=1000)
    W = model.fit_transform(X)  # exposure matrix
    H = model.components_  # signature profiles

    signatures = pd.DataFrame(H.T, index=contexts,
                               columns=[f"SBS_{i+1}" for i in range(n_signatures)])
    exposures = pd.DataFrame(W, index=samples,
                              columns=[f"SBS_{i+1}" for i in range(n_signatures)])

    print(f"Extracted {n_signatures} mutational signatures")
    return signatures, exposures
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/cosmic_mutations.csv` | CSV |
| `results/cbioportal_mutations.csv` | CSV |
| `results/depmap_dependencies.csv` | CSV |
| `results/mutation_signatures.csv` | CSV |
| `results/signature_exposures.csv` | CSV |
| `figures/mutation_spectrum.png` | PNG |
| `figures/signature_profiles.png` | PNG |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| COSMIC | `COSMIC_search_mutations` | 体細胞変異検索 |
| COSMIC | `COSMIC_get_mutations_by_gene` | 遺伝子別変異取得 |
| cBioPortal | `cBioPortal_get_cancer_studies` | がん研究一覧 |
| cBioPortal | `cBioPortal_get_mutations` | 変異データ取得 |
| cBioPortal | `cBioPortal_get_molecular_profiles` | 分子プロファイル |
| cBioPortal | `cBioPortal_get_patients` | 患者データ取得 |
| cBioPortal | `cBioPortal_get_sample_lists` | サンプルリスト |
| cBioPortal | `cBioPortal_get_samples` | サンプル詳細 |
| DepMap | `DepMap_get_gene_dependencies` | 遺伝子依存性スコア |
| DepMap | `DepMap_get_cell_line` | 細胞株情報 |
| DepMap | `DepMap_get_cell_lines` | 細胞株一覧 |
| DepMap | `DepMap_search_cell_lines` | 細胞株検索 |
| DepMap | `DepMap_search_genes` | 遺伝子検索 |

### 参照スキル

| スキル | 関連 |
|---|---|
| `scientific-precision-oncology` | 腫瘍プロファイル → 治療選択 |
| `scientific-variant-interpretation` | バリアント臨床解釈 |
| `scientific-variant-effect-prediction` | 計算病原性予測 |
| `scientific-disease-research` | GWAS → がんリスク |
| `scientific-drug-target-profiling` | 標的同定 → 依存性 |

### 依存パッケージ

`pandas`, `numpy`, `requests`, `scikit-learn`, `matplotlib`

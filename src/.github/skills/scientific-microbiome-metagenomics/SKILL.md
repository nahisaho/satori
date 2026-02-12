---
name: scientific-microbiome-metagenomics
description: |
  マイクロバイオーム・メタゲノミクス解析スキル。16S rRNA アンプリコン解析（DADA2）・
  ショットガンメタゲノム解析（MetaPhlAn / HUMAnN）・α/β 多様性・
  差次存在量解析（DESeq2 / ANCOM-BC）・機能的プロファイリング・
  組成データ解析（CoDA）パイプライン。
---

# Scientific Microbiome & Metagenomics

マイクロバイオーム解析の標準パイプラインを提供する。
16S rRNA アンプリコンおよびショットガンメタゲノムデータの
品質管理、分類学的プロファイリング、多様性評価、
差次存在量解析、機能的アノテーションを体系的に扱う。

## When to Use

- 16S rRNA アンプリコンシーケンスの解析が必要なとき
- ショットガンメタゲノムの分類学的・機能的プロファイリングを行うとき
- 群集の α / β 多様性を比較するとき
- 群間で差次存在量の微生物を同定するとき
- 組成データ（compositional data）の統計解析を行うとき

---

## Quick Start

## 1. 16S rRNA アンプリコン解析（DADA2）

```python
import numpy as np
import pandas as pd

def dada2_pipeline(fastq_dir, trim_left=20, trunc_len_f=240, trunc_len_r=200,
                    min_overlap=12):
    """
    DADA2 アンプリコン解析パイプライン。

    手順:
      1. filterAndTrim — 品質フィルタリング + プライマー除去
      2. learnErrors — エラーモデル学習
      3. dada — ASV（Amplicon Sequence Variant）推定
      4. mergePairs — ペアエンドマージ
      5. removeBimeraDenovo — キメラ除去
      6. assignTaxonomy — SILVA/GTDB による分類

    ASV vs OTU:
      ASV は 100% 配列同一性で分解（1 塩基差を区別）
      OTU は 97% 類似度でクラスタリング（旧来法）
    """
    import subprocess

    r_script = f"""
    library(dada2)

    path <- "{fastq_dir}"
    fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names=TRUE))
    fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names=TRUE))

    # Filter and trim
    filtFs <- file.path(path, "filtered", basename(fnFs))
    filtRs <- file.path(path, "filtered", basename(fnRs))
    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                          trimLeft={trim_left}, truncLen=c({trunc_len_f},{trunc_len_r}),
                          maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE)

    # Error learning
    errF <- learnErrors(filtFs, multithread=TRUE)
    errR <- learnErrors(filtRs, multithread=TRUE)

    # Denoise
    dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
    dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

    # Merge
    merged <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap={min_overlap})

    # ASV table
    seqtab <- makeSequenceTable(merged)
    seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus")

    # Taxonomy
    taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz")

    write.csv(seqtab.nochim, "results/asv_table.csv")
    write.csv(taxa, "results/taxonomy.csv")
    """

    with open("_dada2_pipeline.R", "w") as f:
        f.write(r_script)
    subprocess.run(["Rscript", "_dada2_pipeline.R"], check=True)

    asv_table = pd.read_csv("results/asv_table.csv", index_col=0)
    taxonomy = pd.read_csv("results/taxonomy.csv", index_col=0)
    print(f"  DADA2: {asv_table.shape[1]} ASVs from {asv_table.shape[0]} samples")
    return asv_table, taxonomy
```

## 2. ショットガン分類学的プロファイリング

```python
def shotgun_taxonomic_profiling(fastq_files, method="metaphlan"):
    """
    ショットガンメタゲノム分類学的プロファイリング。

    method:
      - "metaphlan": MetaPhlAn 4 — clade-specific marker 遺伝子ベース
      - "kraken2": Kraken2 — k-mer ベース（高速、メモリ大）
      - "sourmash": sourmash — MinHash ベース

    MetaPhlAn: 精度重視（微量種の検出に優れる）
    Kraken2: 速度重視（大規模データ向け）
    """
    import subprocess

    profiles = []
    for fq in fastq_files:
        sample = fq.split("/")[-1].replace(".fastq.gz", "")

        if method == "metaphlan":
            cmd = (f"metaphlan {fq} --input_type fastq "
                   f"--nproc 8 -o {sample}_profile.txt "
                   f"--bowtie2out {sample}.bt2out")
        elif method == "kraken2":
            cmd = (f"kraken2 --db kraken2_db --threads 8 "
                   f"--report {sample}_report.txt "
                   f"--output {sample}_kraken.txt {fq}")

        subprocess.run(cmd, shell=True, check=True)
        profile = pd.read_csv(f"{sample}_profile.txt", sep="\t",
                                comment="#", header=None)
        profile["sample"] = sample
        profiles.append(profile)

    merged = pd.concat(profiles, ignore_index=True)
    print(f"  Profiling ({method}): {len(fastq_files)} samples processed")
    return merged
```

## 3. α / β 多様性解析

```python
from scipy.spatial.distance import braycurtis, pdist, squareform
from scipy.stats import mannwhitneyu, kruskal
from skbio.diversity import alpha_diversity, beta_diversity

def alpha_diversity_analysis(asv_table, metadata, group_col,
                               metrics=None):
    """
    α 多様性（群集内多様性）解析。

    指標:
      - observed_features: 観察種数（Richness）
      - shannon: Shannon entropy H' = -Σ pᵢ ln(pᵢ)
      - simpson: Simpson index D = 1 - Σ pᵢ²
      - chao1: Chao1 推定種数 S_est = S_obs + f₁²/(2·f₂)
      - faith_pd: Faith's Phylogenetic Diversity（系統的多様性）
    """
    if metrics is None:
        metrics = ["observed_features", "shannon", "simpson", "chao1"]

    results = {}
    for metric in metrics:
        values = alpha_diversity(metric, asv_table.values, asv_table.index)
        results[metric] = values

    alpha_df = pd.DataFrame(results, index=asv_table.index)
    alpha_df = alpha_df.join(metadata[[group_col]])

    # 群間比較
    groups = alpha_df[group_col].unique()
    comparisons = {}
    for metric in metrics:
        if len(groups) == 2:
            g1 = alpha_df[alpha_df[group_col] == groups[0]][metric]
            g2 = alpha_df[alpha_df[group_col] == groups[1]][metric]
            stat, pval = mannwhitneyu(g1, g2)
        else:
            group_data = [alpha_df[alpha_df[group_col] == g][metric] for g in groups]
            stat, pval = kruskal(*group_data)
        comparisons[metric] = {"statistic": stat, "p_value": pval}

    print(f"  α diversity: {len(metrics)} indices computed for {len(alpha_df)} samples")
    return alpha_df, comparisons


def beta_diversity_analysis(asv_table, metadata, group_col,
                              metric="braycurtis", n_perms=999):
    """
    β 多様性（群集間距離）解析。

    距離指標:
      - braycurtis: Bray-Curtis dissimilarity
      - jaccard: Jaccard distance
      - unifrac: UniFrac（系統考慮、ツリー必要）
      - aitchison: Aitchison distance（CoDA 推奨）

    統計検定:
      - PERMANOVA (adonis2): 群間距離の有意差
      - PERMDISP: 分散均一性検定
    """
    dm = beta_diversity(metric, asv_table.values, asv_table.index)

    # PERMANOVA
    from skbio.stats.distance import permanova
    groups = metadata.loc[asv_table.index, group_col]
    permanova_result = permanova(dm, groups, permutations=n_perms)

    # PCoA
    from skbio.stats.ordination import pcoa
    pcoa_result = pcoa(dm)

    print(f"  β diversity ({metric}): PERMANOVA R²={permanova_result['test statistic']:.4f}, "
          f"p={permanova_result['p-value']:.4f}")
    return dm, pcoa_result, permanova_result
```

## 4. 差次存在量解析

```python
def differential_abundance(asv_table, metadata, group_col,
                             formula="~group", method="ancombc"):
    """
    差次存在量解析 — 群間で有意に異なる微生物の同定。

    method:
      - "ancombc": ANCOM-BC2 — バイアス補正・組成データ対応（推奨）
      - "deseq2": DESeq2 — 負の二項分布（RNA-seq 由来）
      - "aldex2": ALDEx2 — CLR 変換 + 効果量

    組成データの問題:
      相対存在量は合計=1 の制約がありスプリアス相関を生む。
      CLR 変換: clr(x) = log(xᵢ / geometric_mean(x))
    """
    import subprocess

    if method == "ancombc":
        r_script = f"""
        library(ANCOMBC)
        library(phyloseq)
        # ANCOM-BC2 analysis
        res <- ancombc2(data=ps, fix_formula="{formula}",
                         p_adj_method="holm", alpha=0.05)
        write.csv(res$res, "results/da_results.csv")
        """
        with open("_da_analysis.R", "w") as f:
            f.write(r_script)
        subprocess.run(["Rscript", "_da_analysis.R"], check=True)
        results = pd.read_csv("results/da_results.csv", index_col=0)

    elif method == "deseq2":
        r_script = f"""
        library(DESeq2)
        dds <- DESeqDataSetFromMatrix(countData=asv_counts,
                                       colData=sample_data,
                                       design={formula})
        dds <- DESeq(dds)
        res <- results(dds)
        write.csv(as.data.frame(res), "results/da_results.csv")
        """
        with open("_da_analysis.R", "w") as f:
            f.write(r_script)
        subprocess.run(["Rscript", "_da_analysis.R"], check=True)
        results = pd.read_csv("results/da_results.csv", index_col=0)

    n_sig = (results.get("padj", results.get("q_val", pd.Series())) < 0.05).sum()
    print(f"  DA ({method}): {n_sig} differentially abundant taxa")
    return results
```

## 5. 機能的プロファイリング

```python
def functional_profiling(fastq_files, method="humann"):
    """
    メタゲノム機能的プロファイリング。

    method:
      - "humann": HUMAnN 3 — UniRef90/MetaCyc パスウェイ
      - "picrust2": PICRUSt2 — 16S から機能予測

    HUMAnN 出力:
      1. Gene families (UniRef90/UniRef50)
      2. Pathway abundance (MetaCyc)
      3. Pathway coverage
    """
    import subprocess

    for fq in fastq_files:
        sample = fq.split("/")[-1].replace(".fastq.gz", "")
        cmd = (f"humann --input {fq} --output humann_results/{sample}/ "
               f"--threads 8 --nucleotide-database chocophlan "
               f"--protein-database uniref")
        subprocess.run(cmd, shell=True, check=True)

    # 結果のマージ
    subprocess.run("humann_join_tables -i humann_results/ -o results/pathway_abundance.tsv "
                    "--file_name pathabundance", shell=True, check=True)
    subprocess.run("humann_join_tables -i humann_results/ -o results/genefamilies.tsv "
                    "--file_name genefamilies", shell=True, check=True)

    pathways = pd.read_csv("results/pathway_abundance.tsv", sep="\t", index_col=0)
    print(f"  HUMAnN: {pathways.shape[0]} pathways across {pathways.shape[1]} samples")
    return pathways
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/asv_table.csv` | CSV |
| `results/taxonomy.csv` | CSV |
| `results/alpha_diversity.csv` | CSV |
| `results/beta_distance_matrix.csv` | CSV |
| `results/da_results.csv` | CSV |
| `results/pathway_abundance.tsv` | TSV |
| `figures/alpha_boxplot.png` | PNG |
| `figures/pcoa_plot.png` | PNG |
| `figures/barplot_taxonomy.png` | PNG |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| MGnify | `MGnify_search_studies` | メタゲノム研究検索 |
| MGnify | `MGnify_list_analyses` | メタゲノム解析一覧 |
| KEGG | `kegg_get_pathway_info` | 代謝パスウェイ情報 |
| KEGG | `kegg_search_pathway` | パスウェイ検索 |
| MetaCyc | `MetaCyc_search_pathways` | 代謝経路検索 |
| PubMed | `PubMed_search_articles` | マイクロバイオーム文献検索 |

### 参照スキル

| スキル | 連携内容 |
|---|---|
| [scientific-metabolomics](../scientific-metabolomics/SKILL.md) | 代謝物-微生物相関 |
| [scientific-network-analysis](../scientific-network-analysis/SKILL.md) | 微生物共起ネットワーク |
| [scientific-statistical-testing](../scientific-statistical-testing/SKILL.md) | 多重検定補正 |
| [scientific-multi-omics](../scientific-multi-omics/SKILL.md) | マルチオミクス統合 |
| [scientific-causal-inference](../scientific-causal-inference/SKILL.md) | 因果推論（微生物-表現型） |

#### 依存パッケージ

- scikit-bio, biom-format, qiime2, dada2 (R), ANCOM-BC (R), DESeq2 (R), HUMAnN, MetaPhlAn

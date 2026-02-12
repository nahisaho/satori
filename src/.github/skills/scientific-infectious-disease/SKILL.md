---
name: scientific-infectious-disease
description: |
  感染症ゲノミクス・疫学スキル。病原体ゲノム解析（SNP/系統樹）・
  AMR（薬剤耐性）遺伝子検出・分子疫学（MLST/cgMLST）・
  アウトブレイク調査トレーシング・疫学的 SIR/SEIR コンパートメントモデル・
  伝播ネットワーク推定パイプライン。
---

# Scientific Infectious Disease Genomics

病原体ゲノミクスと感染症疫学の統合解析パイプラインを提供する。
病原体配列タイピング、系統解析、薬剤耐性遺伝子検出、
アウトブレイク伝播推定、数理疫学モデルを体系的に扱う。

## When to Use

- 病原体の全ゲノムシーケンスデータの解析が必要なとき
- 薬剤耐性（AMR）遺伝子を検出・分類するとき
- 分子疫学タイピング（MLST, cgMLST, SNP）を行うとき
- アウトブレイクの伝播経路を推定するとき
- SIR / SEIR 等のコンパートメントモデルで感染拡大をシミュレーションするとき

---

## Quick Start

## 1. 病原体ゲノム前処理

```python
import numpy as np
import pandas as pd

def pathogen_qc_pipeline(fastq_r1, fastq_r2, reference_genome,
                           min_depth=30, min_coverage=0.95):
    """
    病原体 WGS 前処理パイプライン。

    手順:
      1. Fastp — read QC + adapter trimming
      2. BWA-MEM2 — リファレンスマッピング
      3. Samtools / Picard — dupmark + sort
      4. FreeBayes / GATK — variant calling
      5. カバレッジ / 深度 QC

    品質基準:
      - mean_depth ≥ min_depth (既定: 30x)
      - genome_coverage ≥ min_coverage (既定: 95%)
    """
    import subprocess

    cmds = [
        # QC + trimming
        f"fastp -i {fastq_r1} -I {fastq_r2} -o trim_R1.fq.gz -O trim_R2.fq.gz "
        f"--json qc_report.json",
        # Mapping
        f"bwa-mem2 mem -t 8 {reference_genome} trim_R1.fq.gz trim_R2.fq.gz | "
        f"samtools sort -@ 4 -o aligned.bam",
        # Mark duplicates
        f"samtools markdup aligned.bam dedup.bam",
        f"samtools index dedup.bam",
        # Variant calling
        f"freebayes -f {reference_genome} dedup.bam > variants.vcf",
        # Coverage stats
        f"samtools depth -a dedup.bam | awk '{{sum+=$3; n++}} END {{print sum/n}}'"
    ]

    for cmd in cmds:
        subprocess.run(cmd, shell=True, check=True)

    print(f"  Pipeline complete: variants.vcf generated")
    return "variants.vcf"
```

## 2. AMR 遺伝子検出

```python
def detect_amr_genes(assembly_fasta, database="resfinder"):
    """
    薬剤耐性（AMR）遺伝子の検出。

    データベース:
      - ResFinder: 後天性耐性遺伝子
      - CARD (RGI): 包括的 AMR データベース
      - AMRFinderPlus: NCBI 統合 AMR 検出

    結果カテゴリ:
      - 耐性遺伝子（acquired resistance genes）
      - 点変異（point mutations）
      - 耐性表現型予測
    """
    import subprocess
    import json

    if database == "resfinder":
        cmd = (f"python -m resfinder -ifa {assembly_fasta} "
               f"--acquired --point -o resfinder_results/")
        subprocess.run(cmd, shell=True, check=True)

        with open("resfinder_results/ResFinder_results_tab.txt") as f:
            lines = f.readlines()
        results = parse_resfinder_output(lines)

    elif database == "card":
        cmd = f"rgi main -i {assembly_fasta} -o rgi_results -t contig -a BLAST"
        subprocess.run(cmd, shell=True, check=True)
        results = pd.read_csv("rgi_results.txt", sep="\t")

    n_genes = len(results) if isinstance(results, list) else len(results)
    print(f"  AMR: {n_genes} resistance genes detected ({database})")
    return results


def parse_resfinder_output(lines):
    """ResFinder 出力をパースする。"""
    results = []
    for line in lines[1:]:
        fields = line.strip().split("\t")
        if len(fields) >= 6:
            results.append({
                "gene": fields[0],
                "identity": float(fields[1]),
                "coverage": float(fields[2]),
                "phenotype": fields[5] if len(fields) > 5 else "Unknown",
            })
    return results
```

## 3. 分子疫学タイピング

```python
def molecular_typing(assembly_fasta, organism, scheme="mlst"):
    """
    分子疫学タイピング。

    scheme:
      - "mlst": Multi-Locus Sequence Typing（7 遺伝子座）
      - "cgmlst": core genome MLST（数百〜数千遺伝子座）
      - "wgmlst": whole genome MLST

    MLST:
      各ハウスキーピング遺伝子座のアリル番号の組み合わせで
      Sequence Type（ST）を決定する。
    """
    import subprocess

    if scheme == "mlst":
        cmd = f"mlst {assembly_fasta} --scheme {organism}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        fields = result.stdout.strip().split("\t")
        typing = {
            "file": fields[0],
            "scheme": fields[1],
            "ST": fields[2],
            "alleles": fields[3:],
        }
    elif scheme == "cgmlst":
        cmd = f"chewbbaca AlleleCall -i {assembly_fasta} -g schema/ -o cgmlst_results/"
        subprocess.run(cmd, shell=True, check=True)
        typing = {"scheme": "cgMLST", "results_dir": "cgmlst_results/"}

    print(f"  Typing: ST={typing.get('ST', 'N/A')} ({scheme})")
    return typing
```

## 4. 系統解析・伝播推定

```python
def phylogenetic_analysis(alignment_fasta, method="iqtree", model="GTR+G"):
    """
    病原体系統解析パイプライン。

    method:
      - "iqtree": IQ-TREE 2 — 最尤法（ModelFinder 自動モデル選択）
      - "raxml": RAxML-NG — 最尤法
      - "beast": BEAST 2 — ベイズ系統年代学

    アウトブレイク推定:
      - SNP 距離行列 → 最小スパニングツリー
      - tMRCA (最近共通祖先時間) 推定
    """
    import subprocess
    from Bio import Phylo

    if method == "iqtree":
        cmd = (f"iqtree2 -s {alignment_fasta} -m {model} "
               f"-bb 1000 -alrt 1000 -nt AUTO")
        subprocess.run(cmd, shell=True, check=True)
        tree = Phylo.read(f"{alignment_fasta}.treefile", "newick")

    return tree


def transmission_network(snp_matrix, max_snp_distance=10):
    """
    SNP 距離ベースの伝播ネットワーク推定。

    基準:
      - 直接伝播: SNP 距離 ≤ max_snp_distance
      - 近縁クラスタ: SNP 距離 ≤ 2 × max_snp_distance

    アルゴリズム:
      1. ペアワイズ SNP 距離計算
      2. 閾値以下のペアをエッジとして接続
      3. 最小スパニングツリーで伝播方向推定
    """
    import networkx as nx

    G = nx.Graph()
    samples = snp_matrix.index.tolist()
    G.add_nodes_from(samples)

    for i, s1 in enumerate(samples):
        for j, s2 in enumerate(samples):
            if i < j:
                dist = snp_matrix.iloc[i, j]
                if dist <= max_snp_distance:
                    G.add_edge(s1, s2, weight=dist, snp_distance=dist)

    mst = nx.minimum_spanning_tree(G)
    clusters = list(nx.connected_components(G))

    print(f"  Transmission: {G.number_of_edges()} links, "
          f"{len(clusters)} clusters")
    return G, mst, clusters
```

## 5. SIR / SEIR コンパートメントモデル

```python
from scipy.integrate import odeint

def sir_model(y, t, beta, gamma, N):
    """
    SIR コンパートメントモデル。

    dS/dt = -β · S · I / N
    dI/dt =  β · S · I / N - γ · I
    dR/dt =  γ · I

    R₀ = β / γ (基本再生産数)
    """
    S, I, R = y
    dSdt = -beta * S * I / N
    dIdt = beta * S * I / N - gamma * I
    dRdt = gamma * I
    return [dSdt, dIdt, dRdt]


def seir_model(y, t, beta, sigma, gamma, N):
    """
    SEIR コンパートメントモデル（潜伏期あり）。

    dS/dt = -β · S · I / N
    dE/dt =  β · S · I / N - σ · E
    dI/dt =  σ · E - γ · I
    dR/dt =  γ · I

    σ: 潜伏期の逆数 (1/incubation_period)
    """
    S, E, I, R = y
    dSdt = -beta * S * I / N
    dEdt = beta * S * I / N - sigma * E
    dIdt = sigma * E - gamma * I
    dRdt = gamma * I
    return [dSdt, dEdt, dIdt, dRdt]


def run_epidemic_simulation(model="SIR", N=1e6, I0=10, R0=2.5,
                              gamma=1/10, sigma=1/5, days=180):
    """
    感染症拡大シミュレーション。

    Parameters:
        R0: 基本再生産数
        gamma: 回復率 (1/感染期間)
        sigma: 発症率 (1/潜伏期間、SEIR のみ)
        days: シミュレーション日数
    """
    beta = R0 * gamma
    t = np.linspace(0, days, days * 10)

    if model == "SIR":
        y0 = [N - I0, I0, 0]
        sol = odeint(sir_model, y0, t, args=(beta, gamma, N))
        df = pd.DataFrame(sol, columns=["S", "I", "R"])
    elif model == "SEIR":
        y0 = [N - I0, 0, I0, 0]
        sol = odeint(seir_model, y0, t, args=(beta, sigma, gamma, N))
        df = pd.DataFrame(sol, columns=["S", "E", "I", "R"])

    df["t"] = t
    peak_I = df["I"].max()
    peak_day = df.loc[df["I"].idxmax(), "t"]

    print(f"  {model}: R₀={R0:.1f}, peak infection={peak_I:.0f} at day {peak_day:.0f}")
    return df
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/amr_genes.csv` | CSV |
| `results/mlst_typing.json` | JSON |
| `results/snp_matrix.csv` | CSV |
| `results/transmission_network.json` | JSON |
| `results/epidemic_simulation.csv` | CSV |
| `figures/phylogenetic_tree.png` | PNG |
| `figures/transmission_network.png` | PNG |
| `figures/epidemic_curves.png` | PNG |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| EUHealthInfo | `euhealthinfo_search_infectious_diseases` | 感染症サーベイランスデータ |
| EUHealthInfo | `euhealthinfo_search_surveillance` | 疫学サーベイランス |
| CDC | `cdc_data_search_datasets` | CDC データセット検索 |
| CDC | `cdc_data_get_dataset` | CDC データ取得 |
| NCBI | `BLAST_nucleotide_search` | 病原体配列同定 |
| NCBI | `NCBI_get_sequence` | ゲノム配列取得 |
| PubMed | `PubMed_search_articles` | 感染症文献検索 |
| ClinicalTrials | `search_clinical_trials` | 感染症治療臨床試験 |

### 参照スキル

| スキル | 連携内容 |
|---|---|
| [scientific-sequence-analysis](../scientific-sequence-analysis/SKILL.md) | 配列アライメント・BLAST |
| [scientific-bioinformatics](../scientific-bioinformatics/SKILL.md) | ゲノムアノテーション |
| [scientific-network-analysis](../scientific-network-analysis/SKILL.md) | 伝播ネットワーク可視化 |
| [scientific-survival-clinical](../scientific-survival-clinical/SKILL.md) | 感染症アウトカム解析 |
| [scientific-bayesian-statistics](../scientific-bayesian-statistics/SKILL.md) | ベイズ系統年代学 |

#### 依存パッケージ

- biopython, ete3, scipy, networkx, subprocess (fastp, bwa-mem2, freebayes, iqtree2)

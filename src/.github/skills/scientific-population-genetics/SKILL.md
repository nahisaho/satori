---
name: scientific-population-genetics
description: |
  集団遺伝学解析スキル。アレル頻度解析・Hardy-Weinberg 平衡検定・
  集団構造解析（PCA / ADMIXTURE）・Fst 分化指標・選択圧検出（iHS / XP-EHH）・
  連鎖不平衡（LD）解析・GWAS Catalog / gnomAD データ統合パイプライン。
---

# Scientific Population Genetics

集団遺伝学に特化した解析パイプラインを提供する。
アレル頻度、集団構造、遺伝的分化、自然選択シグナル、
連鎖不平衡の解析を体系的に扱い、GWAS Catalog・gnomAD との統合を支援する。

## When to Use

- アレル頻度分布や Hardy-Weinberg 平衡を検定するとき
- 集団構造（PCA / ADMIXTURE / STRUCTURE）を解析するとき
- 集団間の遺伝的分化（Fst）を評価するとき
- 自然選択シグナル（iHS / Tajima's D / XP-EHH）を検出するとき
- GWAS 関連バリアントの集団遺伝学的解釈を行うとき

---

## Quick Start

## 1. QC・アレル頻度解析

```python
import numpy as np
import pandas as pd

def genotype_qc(plink_prefix, mind=0.02, geno=0.02, maf=0.01,
                  hwe_p=1e-6):
    """
    ジェノタイプ QC パイプライン（PLINK 2）。

    フィルタリング基準:
      - --mind: 個体ミッシング率 ≤ mind
      - --geno: SNP ミッシング率 ≤ geno
      - --maf: Minor Allele Frequency ≥ maf
      - --hwe: Hardy-Weinberg p ≥ hwe_p（コントロールのみ）

    追加 QC:
      - 性別不一致チェック
      - IBD 推定（近親者除外: π̂ > 0.25）
      - PCA アウトライアー除外
    """
    import subprocess

    # Step 1: ミッシング率フィルタ
    cmd = (f"plink2 --bfile {plink_prefix} "
           f"--mind {mind} --geno {geno} --maf {maf} "
           f"--hwe {hwe_p} "
           f"--make-bed --out {plink_prefix}_qc")
    subprocess.run(cmd, shell=True, check=True)

    # Step 2: IBD 推定（近親者検出）
    cmd = (f"plink2 --bfile {plink_prefix}_qc "
           f"--indep-pairwise 50 5 0.2 --out {plink_prefix}_prune")
    subprocess.run(cmd, shell=True, check=True)

    cmd = (f"plink2 --bfile {plink_prefix}_qc "
           f"--extract {plink_prefix}_prune.prune.in "
           f"--genome --out {plink_prefix}_ibd")
    subprocess.run(cmd, shell=True, check=True)

    return f"{plink_prefix}_qc"


def allele_frequency_stats(genotype_matrix, populations):
    """
    集団別アレル頻度統計。

    算出指標:
      - MAF: Minor Allele Frequency
      - Het: Observed heterozygosity = n_het / n_total
      - Expected Het (He): 2pq
      - HWE: Hardy-Weinberg 平衡検定 (χ² test)
        H₀: f(AA) = p², f(Aa) = 2pq, f(aa) = q²
    """
    from scipy.stats import chi2

    results = []
    for pop in populations["population"].unique():
        pop_idx = populations[populations["population"] == pop].index
        geno_pop = genotype_matrix.loc[pop_idx]

        for snp in geno_pop.columns:
            counts = geno_pop[snp].value_counts()
            n = counts.sum()
            n_0 = counts.get(0, 0)  # AA
            n_1 = counts.get(1, 0)  # Aa
            n_2 = counts.get(2, 0)  # aa

            p = (2 * n_0 + n_1) / (2 * n)
            q = 1 - p
            maf = min(p, q)

            # HWE test
            exp_0 = n * p**2
            exp_1 = n * 2*p*q
            exp_2 = n * q**2
            if exp_0 > 0 and exp_1 > 0 and exp_2 > 0:
                chi2_stat = ((n_0-exp_0)**2/exp_0 + (n_1-exp_1)**2/exp_1 +
                              (n_2-exp_2)**2/exp_2)
                hwe_p = 1 - chi2.cdf(chi2_stat, df=1)
            else:
                hwe_p = 1.0

            het_obs = n_1 / n
            het_exp = 2 * p * q

            results.append({
                "snp": snp, "population": pop,
                "MAF": round(maf, 4), "p": round(p, 4),
                "Het_obs": round(het_obs, 4), "Het_exp": round(het_exp, 4),
                "HWE_p": round(hwe_p, 6),
            })

    return pd.DataFrame(results)
```

## 2. 集団構造解析

```python
def population_structure(plink_prefix, n_components=10, method="pca"):
    """
    集団構造解析。

    method:
      - "pca": 主成分分析 — 集団間の遺伝的差異を 2D/3D で可視化
      - "admixture": ADMIXTURE — 各個体の祖先集団比率を推定
        K=2〜10 を試行し、CV error 最小の K を選択

    PCA on genotypes:
      X を (n_samples × n_snps) ジェノタイプ行列として
      共分散行列 C = XᵀX / n_snps の固有値分解
    """
    import subprocess

    if method == "pca":
        cmd = (f"plink2 --bfile {plink_prefix} "
               f"--pca {n_components} --out {plink_prefix}_pca")
        subprocess.run(cmd, shell=True, check=True)

        eigenvec = pd.read_csv(f"{plink_prefix}_pca.eigenvec", sep="\t")
        eigenval = pd.read_csv(f"{plink_prefix}_pca.eigenval", header=None)
        var_explained = eigenval[0] / eigenval[0].sum()

        print(f"  PCA: PC1={var_explained[0]:.3f}, PC2={var_explained[1]:.3f}")
        return eigenvec, var_explained

    elif method == "admixture":
        cv_errors = {}
        for K in range(2, 11):
            cmd = f"admixture --cv {plink_prefix}.bed {K}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            # CV error 抽出
            for line in result.stdout.split("\n"):
                if "CV error" in line:
                    cv_errors[K] = float(line.split(": ")[1])

        best_K = min(cv_errors, key=cv_errors.get)
        Q = pd.read_csv(f"{plink_prefix}.{best_K}.Q", sep=" ", header=None)

        print(f"  ADMIXTURE: best K={best_K} (CV error={cv_errors[best_K]:.4f})")
        return Q, cv_errors, best_K
```

## 3. 遺伝的分化（Fst）

```python
def calculate_fst(genotype_matrix, populations, method="weir_cockerham"):
    """
    集団間遺伝的分化指標 Fst 算出。

    Weir-Cockerham (1984) 推定量:
      F_ST = σ²_a / (σ²_a + σ²_b + σ²_w)
      σ²_a: 集団間分散
      σ²_b: 集団内個体間分散
      σ²_w: 個体内（アレル間）分散

    解釈:
      Fst = 0: 分化なし（パンミクシア）
      0 < Fst < 0.05: 低分化
      0.05 ≤ Fst < 0.15: 中程度の分化
      0.15 ≤ Fst < 0.25: 大きな分化
      Fst ≥ 0.25: 非常に大きな分化

    genome-wide Fst: 全 SNP の加重平均
    per-SNP Fst: 局所的適応シグナルの検出
    """
    pop_labels = populations["population"]
    unique_pops = pop_labels.unique()

    fst_per_snp = []
    for snp in genotype_matrix.columns:
        # 集団別アレル頻度
        pop_freqs = {}
        pop_sizes = {}
        for pop in unique_pops:
            idx = pop_labels[pop_labels == pop].index
            geno = genotype_matrix.loc[idx, snp].dropna()
            p = (2 * (geno == 0).sum() + (geno == 1).sum()) / (2 * len(geno))
            pop_freqs[pop] = p
            pop_sizes[pop] = len(geno)

        # Weir-Cockerham Fst
        n_pops = len(unique_pops)
        n_total = sum(pop_sizes.values())
        p_bar = sum(pop_freqs[p] * pop_sizes[p] for p in unique_pops) / n_total
        n_bar = n_total / n_pops

        MSP = sum(pop_sizes[p] * (pop_freqs[p] - p_bar)**2
                   for p in unique_pops) / (n_pops - 1)
        MSG = sum(pop_sizes[p] * pop_freqs[p] * (1 - pop_freqs[p])
                   for p in unique_pops) / (n_total - n_pops)

        nc = (n_total - sum(n**2 for n in pop_sizes.values()) / n_total) / (n_pops - 1)

        if (MSP + (nc - 1) * MSG) > 0:
            fst = (MSP - MSG) / (MSP + (nc - 1) * MSG)
        else:
            fst = 0

        fst_per_snp.append({"snp": snp, "Fst": max(fst, 0), "p_bar": p_bar})

    fst_df = pd.DataFrame(fst_per_snp)
    genome_fst = fst_df["Fst"].mean()

    print(f"  Fst: genome-wide={genome_fst:.4f}, "
          f"max per-SNP={fst_df['Fst'].max():.4f}")
    return fst_df, genome_fst
```

## 4. 自然選択シグナル検出

```python
def selection_scan(haplotype_matrix, positions, method="ihs"):
    """
    自然選択シグナルの検出。

    method:
      - "ihs": Integrated Haplotype Score — ローカル正の選択
        |iHS| > 2: 選択シグナル候補
      - "tajima_d": Tajima's D — 中立性検定
        D > 0: バランス選択 or 集団縮小
        D < 0: 正の選択 or 集団拡大  
        D ≈ 0: 中立進化
      - "xpehh": Cross-Population EHH — 集団間正の選択

    iHS:
      各 SNP について、派生アレル (derived) と祖先アレル (ancestral) の
      Extended Haplotype Homozygosity (EHH) を比較。
      iHS = ln(iHH_A / iHH_D)  → 標準化
    """
    if method == "tajima_d":
        # スライディングウィンドウ Tajima's D
        from allel import tajima_d
        import allel

        D_values = []
        window_size = 50000
        step = 10000

        for start in range(0, positions[-1], step):
            end = start + window_size
            mask = (positions >= start) & (positions < end)
            if mask.sum() > 5:
                ac = allel.AlleleCountsArray(
                    haplotype_matrix[:, mask].sum(axis=0).reshape(-1, 1))
                D = tajima_d(ac)
                D_values.append({"start": start, "end": end, "D": D,
                                  "n_snps": mask.sum()})

        df = pd.DataFrame(D_values)
        print(f"  Tajima's D: mean={df['D'].mean():.3f}, "
              f"range=[{df['D'].min():.3f}, {df['D'].max():.3f}]")
        return df

    elif method == "ihs":
        import allel
        ihs = allel.ihs(haplotype_matrix, positions)
        # 標準化
        ihs_std = (ihs - np.nanmean(ihs)) / np.nanstd(ihs)

        n_sig = np.sum(np.abs(ihs_std) > 2)
        print(f"  iHS: {n_sig} candidate regions (|iHS|>2)")
        return ihs_std
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/allele_frequencies.csv` | CSV |
| `results/pca_eigenvec.csv` | CSV |
| `results/admixture_Q.csv` | CSV |
| `results/fst_per_snp.csv` | CSV |
| `results/selection_scan.csv` | CSV |
| `figures/pca_populations.png` | PNG |
| `figures/admixture_barplot.png` | PNG |
| `figures/manhattan_fst.png` | PNG |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| gnomAD | `gnomad_get_variant` | バリアント集団頻度 |
| gnomAD | `gnomad_get_gene_constraints` | 遺伝子制約指標 |
| gnomAD | `gnomad_get_region` | 領域別バリアント |
| gnomAD | `gnomad_search_variants` | バリアント検索 |
| GWAS | `GWAS_search_associations_by_gene` | 遺伝子別 GWAS 関連 |
| GWAS | `gwas_search_studies` | GWAS 研究検索 |
| GWAS | `gwas_get_variants_for_trait` | 形質別バリアント |
| GWAS | `gwas_get_associations_for_snp` | SNP 別関連 |
| GWAS | `gwas_get_snps_for_gene` | 遺伝子近傍 SNP |

### 参照スキル

| スキル | 連携内容 |
|---|---|
| [scientific-variant-interpretation](../scientific-variant-interpretation/SKILL.md) | バリアント臨床解釈 |
| [scientific-bioinformatics](../scientific-bioinformatics/SKILL.md) | ゲノムアノテーション |
| [scientific-disease-research](../scientific-disease-research/SKILL.md) | 疾患-遺伝子関連 |
| [scientific-statistical-testing](../scientific-statistical-testing/SKILL.md) | 統計検定 |
| [scientific-pca-tsne](../scientific-pca-tsne/SKILL.md) | 次元削減 |

#### 依存パッケージ

- scikit-allel, plink2, admixture, pandas, numpy, scipy

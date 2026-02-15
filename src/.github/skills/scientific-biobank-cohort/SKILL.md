---
name: scientific-biobank-cohort
description: |
  バイオバンク・大規模コホートデータ解析スキル。UK Biobank /
  BBJ / All of Us 等の大規模コホートデータに対するフェノタイプ
  辞書検索・GWAS サマリー統計処理・PheWAS パイプライン。
tu_tools:
  - key: clinvar
    name: ClinVar
    description: バリアント臨床的意義データ検索
---

# Scientific Biobank Cohort

UK Biobank・バイオバンクジャパン (BBJ)・All of Us 等の大規模
コホートデータを活用したフェノタイプ辞書検索・GWAS サマリー
統計処理・PheWAS 解析パイプラインを提供する。

## When to Use

- バイオバンクのフェノタイプ辞書を検索するとき
- GWAS サマリー統計データを処理・可視化するとき
- PheWAS (Phenome-Wide Association Study) を実施するとき
- コホートの基本統計・人口統計特性を集計するとき
- バリアント-フェノタイプ関連を網羅的に検索するとき

---

## Quick Start

## 1. フェノタイプ辞書検索

```python
import pandas as pd
import numpy as np


def phenotype_dictionary(pheno_file, category=None,
                          keyword=None):
    """
    バイオバンク — フェノタイプ辞書検索。

    Parameters:
        pheno_file: str — フェノタイプ辞書 CSV パス
            (UK Biobank Data-Field listing 等)
        category: str — カテゴリフィルタ
        keyword: str — キーワードフィルタ
    """
    df = pd.read_csv(pheno_file)

    if category:
        df = df[df["Category"].str.contains(
            category, case=False, na=False)]
    if keyword:
        mask = (
            df["Field"].str.contains(
                keyword, case=False, na=False)
            | df["Description"].str.contains(
                keyword, case=False, na=False)
        )
        df = df[mask]

    print(f"Phenotype dict: {len(df)} fields matched")
    return df


def cohort_demographics(pheno_df, age_col="age",
                         sex_col="sex"):
    """
    バイオバンク — コホート人口統計サマリー。

    Parameters:
        pheno_df: DataFrame — 参加者フェノタイプデータ
        age_col: str — 年齢列名
        sex_col: str — 性別列名
    """
    summary = {
        "n_participants": len(pheno_df),
        "age_mean": pheno_df[age_col].mean(),
        "age_std": pheno_df[age_col].std(),
        "sex_distribution": (
            pheno_df[sex_col]
            .value_counts(normalize=True)
            .to_dict()
        ),
    }
    print(f"Cohort: n={summary['n_participants']}, "
          f"age={summary['age_mean']:.1f}±"
          f"{summary['age_std']:.1f}")
    return summary
```

## 2. GWAS サマリー統計処理

```python
def load_gwas_summary(sumstat_file, p_threshold=5e-8,
                       sep="\t"):
    """
    GWAS サマリー統計ファイル読み込み・フィルタリング。

    Parameters:
        sumstat_file: str — サマリー統計ファイルパス
            (TSV: CHR, POS, SNP, A1, A2, BETA, SE, P)
        p_threshold: float — P 値閾値
        sep: str — 区切り文字
    """
    df = pd.read_csv(sumstat_file, sep=sep)

    # 標準カラム名正規化
    col_map = {
        "chromosome": "CHR", "chr": "CHR",
        "position": "POS", "pos": "POS", "bp": "POS",
        "rsid": "SNP", "snp": "SNP", "variant_id": "SNP",
        "effect_allele": "A1", "a1": "A1",
        "other_allele": "A2", "a2": "A2",
        "beta": "BETA", "effect_size": "BETA",
        "se": "SE", "standard_error": "SE",
        "pval": "P", "p_value": "P", "pvalue": "P",
    }
    df.columns = [col_map.get(c.lower(), c)
                  for c in df.columns]

    # フィルタ
    sig = df[df["P"] < p_threshold].copy()
    sig.sort_values("P", inplace=True)

    print(f"GWAS summary: {len(df)} total, "
          f"{len(sig)} significant (P<{p_threshold})")
    return sig


def manhattan_data(gwas_df, chr_col="CHR",
                    pos_col="POS", p_col="P"):
    """
    Manhattan プロット用データ変換。

    Parameters:
        gwas_df: DataFrame — GWAS サマリー統計
        chr_col: str — 染色体列
        pos_col: str — 位置列
        p_col: str — P 値列
    """
    df = gwas_df.copy()
    df["-log10P"] = -np.log10(df[p_col])

    # 累積位置計算
    chr_lengths = (
        df.groupby(chr_col)[pos_col].max()
        .sort_index()
    )
    chr_offsets = chr_lengths.cumsum().shift(1).fillna(0)
    df["cumpos"] = df.apply(
        lambda r: r[pos_col] + chr_offsets.get(
            r[chr_col], 0),
        axis=1)

    print(f"Manhattan data: {len(df)} variants, "
          f"max -log10P={df['-log10P'].max():.1f}")
    return df
```

## 3. PheWAS (Phenome-Wide Association Study)

```python
def phewas_analysis(genotype_series, pheno_df,
                      pheno_cols=None,
                      p_threshold=0.05):
    """
    PheWAS — 1バリアントに対する多表現型アソシエーション。

    Parameters:
        genotype_series: Series — バリアント遺伝子型
            (0/1/2 コーディング)
        pheno_df: DataFrame — フェノタイプデータ
        pheno_cols: list — テスト対象表現型列
        p_threshold: float — Bonferroni 前閾値
    """
    from scipy import stats

    if pheno_cols is None:
        pheno_cols = [c for c in pheno_df.columns
                      if pheno_df[c].dtype in
                      [np.float64, np.int64]]

    results = []
    for col in pheno_cols:
        mask = pheno_df[col].notna()
        if mask.sum() < 50:
            continue
        geno = genotype_series[mask]
        pheno = pheno_df.loc[mask, col]

        # 数値 → 線形回帰 (簡易)
        slope, intercept, r, p, se = stats.linregress(
            geno, pheno)
        results.append({
            "phenotype": col,
            "beta": slope,
            "se": se,
            "p_value": p,
            "n": mask.sum(),
        })

    df = pd.DataFrame(results)
    n_tests = len(df)
    bonf = p_threshold / n_tests if n_tests > 0 else 0.05
    df["significant"] = df["p_value"] < bonf
    df.sort_values("p_value", inplace=True)

    n_sig = df["significant"].sum()
    print(f"PheWAS: {n_tests} phenotypes tested, "
          f"{n_sig} significant (Bonferroni)")
    return df
```

## 4. バイオバンク統合パイプライン

```python
def biobank_pipeline(sumstat_file, pheno_file=None,
                       output_dir="results"):
    """
    バイオバンク統合パイプライン。

    Parameters:
        sumstat_file: str — GWAS サマリー統計ファイル
        pheno_file: str — フェノタイプ辞書ファイル
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) GWAS サマリー統計読み込み
    gwas = load_gwas_summary(sumstat_file)
    gwas.to_csv(output_dir / "gwas_significant.csv",
                index=False)

    # 2) Manhattan プロットデータ
    manhattan = manhattan_data(gwas)
    manhattan.to_csv(
        output_dir / "manhattan_data.csv", index=False)

    # 3) フェノタイプ辞書検索 (利用可能な場合)
    if pheno_file:
        pheno_dict = phenotype_dictionary(pheno_file)
        pheno_dict.to_csv(
            output_dir / "phenotype_dict.csv",
            index=False)

    print(f"Biobank pipeline → {output_dir}")
    return {"gwas": gwas, "manhattan": manhattan}
```

---

## パイプライン統合

```
epidemiology-public-health → biobank-cohort → population-genetics
     (疫学デザイン)           (GWAS/PheWAS)    (集団遺伝解析)
           │                       │                 ↓
   mendelian-randomization ───────┘       rare-disease-genetics
     (因果推論)                            (Mendelian 解析)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/gwas_significant.csv` | Genome-wide significant SNP | → population-genetics |
| `results/manhattan_data.csv` | Manhattan プロットデータ | → GWAS 可視化 |
| `results/phenotype_dict.csv` | フェノタイプ辞書 | → PheWAS |

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `clinvar` | ClinVar | バリアント臨床的意義データ検索 |

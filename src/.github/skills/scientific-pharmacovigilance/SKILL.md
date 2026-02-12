---
name: scientific-pharmacovigilance
description: |
  ファーマコビジランス（医薬品安全性監視）スキル。FAERS/FDA 有害事象報告データベースを活用し、
  不均衡分析（PRR/ROR/IC）、MedDRA 階層構造、時系列トレンド、人口統計層別化を実施。
  市販後安全性シグナル検出と定量的評価を支援。
  「有害事象を分析して」「FAERS データを解析して」「安全性シグナルを検出して」で発火。
---

# Scientific Pharmacovigilance

医薬品安全性監視（ファーマコビジランス）のための解析スキル。
FDA FAERS（FDA Adverse Event Reporting System）を中心とした
市販後安全性データの統合的分析を支援する。

## When to Use

- 市販後有害事象データの不均衡分析（Signal Detection）
- 特定医薬品の安全性プロファイル評価
- MedDRA 階層を用いた有害事象の体系的分類
- 時系列トレンドによるシグナル推移の追跡
- 人口統計（年齢・性別・体重）別の層別化解析
- 薬物間の安全性比較

## Quick Start

### ファーマコビジランス解析パイプライン

```
Phase 1: Data Acquisition
  - FAERS / EudraVigilance / VigiBase からのデータ取得
  - 製品名・成分名の正規化
  - 重複報告の除去（Case ID ベース）
    ↓
Phase 2: MedDRA Coding & Hierarchy
  - LLT → PT → HLT → HLGT → SOC 階層マッピング
  - Preferred Term (PT) レベルでの集計
  - SMQ (Standardised MedDRA Query) 適用
    ↓
Phase 3: Disproportionality Analysis
  - PRR (Proportional Reporting Ratio) 算出
  - ROR (Reporting Odds Ratio) 算出
  - IC (Information Component, Bayesian) 算出
  - EBGM (Empirical Bayes Geometric Mean) 算出
  - シグナル閾値判定
    ↓
Phase 4: Temporal & Demographic Analysis
  - Time-to-Onset 分布解析
  - 四半期別報告トレンド
  - 年齢・性別・体重・適応症別層別化
  - Rechallenge / Dechallenge 分析
    ↓
Phase 5: Signal Evaluation & Reporting
  - シグナル優先順位付け
  - ケースナラティブ分析
  - 因果関係評価（WHO-UMC / Naranjo スケール）
  - 安全性シグナルレポート生成
    ↓
Phase 6: Risk-Benefit Assessment
  - NNH (Number Needed to Harm) 推定
  - リスク-ベネフィットバランス評価
  - REMS / RiskMAP 提言
    ↓
Phase 7: Regulatory Communication
  - PBRER / PSUR セクション対応
  - シグナルサマリー文書生成
  - 添付文書改訂案
```

## Workflow

### 1. FAERS データ取得・前処理

```python
import pandas as pd
import numpy as np
from scipy import stats

# === FAERS データ読み込み ===
# FAERS ASCII ファイル: DEMO, DRUG, REAC, OUTC, INDI, THER, RPSR
demo = pd.read_csv("faers/DEMO.txt", sep="$", dtype=str)
drug = pd.read_csv("faers/DRUG.txt", sep="$", dtype=str)
reac = pd.read_csv("faers/REAC.txt", sep="$", dtype=str)

# 重複除去: primaryid 最新版のみ保持
demo = demo.sort_values("fda_dt", ascending=False).drop_duplicates(subset=["caseid"], keep="first")

# Drug-Reaction ペア構築
merged = drug.merge(reac, on="primaryid", how="inner")
merged = merged.merge(demo[["primaryid", "age", "sex", "wt", "reporter_country"]], on="primaryid", how="left")

# 対象薬剤のフィルタ
target_drug = "ACETAMINOPHEN"
target_df = merged[merged["drugname"].str.upper().str.contains(target_drug)]
print(f"Target drug reports: {len(target_df):,}")
```

### 2. MedDRA 階層マッピング

```python
# === MedDRA 階層 ===
# LLT → PT → HLT → HLGT → SOC
meddra_pt = pd.read_csv("meddra/pt.asc", sep="$", header=None,
                         names=["pt_code", "pt_name", "null_field",
                                "pt_soc_code", "null2", "null3"])
meddra_soc = pd.read_csv("meddra/soc.asc", sep="$", header=None,
                          names=["soc_code", "soc_name", "soc_abbrev", "null1"])

# PT コードを Reaction テーブルに紐付け
merged["pt_name_clean"] = merged["pt"].str.strip().str.upper()

# SOC レベル集計
soc_summary = merged.merge(meddra_pt, left_on="pt_name_clean", right_on="pt_name")
soc_counts = soc_summary.groupby("pt_soc_code").size().reset_index(name="count")
soc_counts = soc_counts.merge(meddra_soc, left_on="pt_soc_code", right_on="soc_code")
soc_counts = soc_counts.sort_values("count", ascending=False)
print("Top SOC categories:")
print(soc_counts[["soc_name", "count"]].head(10).to_string(index=False))
```

### 3. 不均衡分析（Disproportionality Analysis）

```python
def disproportionality_analysis(drug_reactions_df, all_reactions_df, target_drug):
    """
    PRR, ROR, IC (Information Component) を算出。

    2x2 分割表:
                    Target AE   Other AE
    Target Drug      a           b
    Other Drugs      c           d
    """
    results = []

    target_pts = drug_reactions_df[
        drug_reactions_df["drugname"].str.upper().str.contains(target_drug)
    ]["pt"].value_counts()

    total_target = drug_reactions_df[
        drug_reactions_df["drugname"].str.upper().str.contains(target_drug)
    ]["primaryid"].nunique()

    total_all = all_reactions_df["primaryid"].nunique()
    total_other = total_all - total_target

    for pt, a in target_pts.items():
        # b: target drug, other AE
        b = total_target - a
        # c: other drugs, same AE
        c_total = all_reactions_df[all_reactions_df["pt"] == pt]["primaryid"].nunique()
        c = c_total - a
        # d: other drugs, other AE
        d = total_other - c

        # PRR
        if (a + b) > 0 and (c + d) > 0 and c > 0:
            prr = (a / (a + b)) / (c / (c + d))
            # PRR の 95% CI
            se_ln_prr = np.sqrt(1/a - 1/(a+b) + 1/c - 1/(c+d)) if a > 0 else np.inf
            prr_lower = np.exp(np.log(prr) - 1.96 * se_ln_prr) if prr > 0 else 0
        else:
            prr, prr_lower = 0, 0

        # ROR
        if b > 0 and c > 0 and d > 0:
            ror = (a * d) / (b * c)
            se_ln_ror = np.sqrt(1/a + 1/b + 1/c + 1/d) if a > 0 else np.inf
            ror_lower = np.exp(np.log(ror) - 1.96 * se_ln_ror) if ror > 0 else 0
        else:
            ror, ror_lower = 0, 0

        # IC (Information Component) - Bayesian
        expected = (a + b) * (a + c) / (a + b + c + d) if (a + b + c + d) > 0 else 1
        ic = np.log2((a + 0.5) / (expected + 0.5)) if expected > 0 else 0

        # Chi-square
        chi2, p_value = 0, 1
        if a + b + c + d > 0:
            try:
                chi2, p_value, _, _ = stats.chi2_contingency([[a, b], [c, d]])
            except ValueError:
                pass

        # Signal criteria
        is_signal = (prr >= 2 and chi2 >= 4 and a >= 3)
        is_signal_ror = (ror_lower > 1 and a >= 3)
        is_signal_ic = (ic > 0 and a >= 3)

        results.append({
            "pt": pt, "a": a, "b": b, "c": c, "d": d,
            "prr": round(prr, 3), "prr_lower": round(prr_lower, 3),
            "ror": round(ror, 3), "ror_lower": round(ror_lower, 3),
            "ic": round(ic, 3),
            "chi2": round(chi2, 3), "p_value": round(p_value, 6),
            "signal_prr": is_signal,
            "signal_ror": is_signal_ror,
            "signal_ic": is_signal_ic,
        })

    return pd.DataFrame(results).sort_values("prr", ascending=False)

signals = disproportionality_analysis(merged, merged, target_drug)
print(f"Signals detected (PRR≥2, χ²≥4, N≥3): {signals['signal_prr'].sum()}")
print(signals[signals["signal_prr"]].head(20))
```

### 4. 時系列トレンド分析

```python
import matplotlib.pyplot as plt

def temporal_trend_analysis(drug_df, target_drug, target_pt=None):
    """四半期別報告トレンドと Time-to-Onset 分布"""
    target = drug_df[drug_df["drugname"].str.upper().str.contains(target_drug)].copy()
    if target_pt:
        target = target[target["pt"].str.upper().str.contains(target_pt)]

    # 四半期別集計
    target["report_date"] = pd.to_datetime(target["fda_dt"], format="%Y%m%d", errors="coerce")
    target["quarter"] = target["report_date"].dt.to_period("Q")
    quarterly = target.groupby("quarter").size()

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # 報告トレンド
    quarterly.plot(kind="bar", ax=axes[0], color="#2196F3")
    axes[0].set_title(f"Quarterly Reports: {target_drug}")
    axes[0].set_ylabel("Report Count")
    axes[0].tick_params(axis="x", rotation=45)

    # Time-to-Onset
    if "event_dt" in target.columns and "start_dt" in target.columns:
        target["event_date"] = pd.to_datetime(target["event_dt"], format="%Y%m%d", errors="coerce")
        target["start_date"] = pd.to_datetime(target["start_dt"], format="%Y%m%d", errors="coerce")
        target["tto_days"] = (target["event_date"] - target["start_date"]).dt.days
        valid_tto = target["tto_days"].dropna()
        valid_tto = valid_tto[(valid_tto >= 0) & (valid_tto <= 365)]
        axes[1].hist(valid_tto, bins=30, color="#FF9800", edgecolor="black")
        axes[1].set_title("Time-to-Onset Distribution")
        axes[1].set_xlabel("Days")
        axes[1].set_ylabel("Frequency")

    plt.tight_layout()
    plt.savefig("figures/pv_temporal_trend.png", dpi=300, bbox_inches="tight")
    plt.show()

temporal_trend_analysis(merged, target_drug)
```

### 5. 人口統計層別化分析

```python
def demographic_stratification(drug_df, target_drug, target_pt):
    """年齢・性別・体重別の層別不均衡分析"""
    target = drug_df[
        (drug_df["drugname"].str.upper().str.contains(target_drug)) &
        (drug_df["pt"].str.upper().str.contains(target_pt))
    ].copy()

    # 年齢カテゴリ
    target["age_num"] = pd.to_numeric(target["age"], errors="coerce")
    target["age_group"] = pd.cut(target["age_num"],
                                  bins=[0, 18, 40, 65, 85, 120],
                                  labels=["<18", "18-40", "40-65", "65-85", "85+"])

    # 性別分布
    sex_dist = target["sex"].value_counts()

    # Weight-based analysis
    target["wt_num"] = pd.to_numeric(target["wt"], errors="coerce")

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Age distribution
    target["age_group"].value_counts().sort_index().plot(
        kind="bar", ax=axes[0], color="#4CAF50")
    axes[0].set_title(f"Age Distribution: {target_drug} + {target_pt}")
    axes[0].set_ylabel("Count")

    # Sex distribution
    sex_dist.plot(kind="pie", ax=axes[1], autopct="%1.1f%%",
                  colors=["#2196F3", "#E91E63", "#9E9E9E"])
    axes[1].set_title("Sex Distribution")

    # Outcome severity
    if "outc_cod" in target.columns:
        outcome_map = {"DE": "Death", "LT": "Life-Threatening",
                       "HO": "Hospitalization", "DS": "Disability",
                       "CA": "Congenital Anomaly", "OT": "Other"}
        target["outcome_label"] = target["outc_cod"].map(outcome_map)
        target["outcome_label"].value_counts().plot(
            kind="barh", ax=axes[2], color="#F44336")
        axes[2].set_title("Outcome Distribution")

    plt.tight_layout()
    plt.savefig("figures/pv_demographics.png", dpi=300, bbox_inches="tight")
    plt.show()

    return target

demographic_stratification(merged, target_drug, "HEPATOTOXICITY")
```

### 6. 因果関係評価

```python
def naranjo_assessment(case_data):
    """
    Naranjo Adverse Drug Reaction Probability Scale
    スコア: ≥9 = Definite, 5-8 = Probable, 1-4 = Possible, ≤0 = Doubtful
    """
    questions = [
        ("以前に確立された有害反応か", 1, 0, 0),
        ("薬剤投与後に出現したか", 2, -1, 0),
        ("中止で改善したか (Dechallenge)", 1, 0, 0),
        ("再投与で再発したか (Rechallenge)", 2, -1, 0),
        ("他の原因が除外できるか", -1, 2, 0),
        ("プラセボでも出現するか", -1, 1, 0),
        ("血中濃度は中毒域か", 1, 0, 0),
        ("用量依存性があるか", 1, 0, 0),
        ("類似薬で既往があるか", 1, 0, 0),
        ("客観的所見で確認されたか", 1, 0, 0),
    ]

    total_score = 0
    assessment = []
    for q, yes_score, no_score, dk_score in questions:
        # 実際のケースデータから判定
        answer = case_data.get(q, "dk")
        if answer == "yes":
            score = yes_score
        elif answer == "no":
            score = no_score
        else:
            score = dk_score
        total_score += score
        assessment.append({"question": q, "answer": answer, "score": score})

    if total_score >= 9:
        category = "Definite"
    elif total_score >= 5:
        category = "Probable"
    elif total_score >= 1:
        category = "Possible"
    else:
        category = "Doubtful"

    return {
        "total_score": total_score,
        "category": category,
        "details": assessment,
    }
```

### 7. EBGM (Empirical Bayes Geometric Mean)

```python
def calculate_ebgm(contingency_df, shrinkage_prior=0.5):
    """
    EBGM (Multi-item Gamma Poisson Shrinker) による
    ベイジアンシグナル検出。FDA OPIS で使用。

    EBGM = exp(E[log(λ)|N])
    EB05 = EBGM の下側 5% 信頼限界
    """
    results = []
    N_total = contingency_df[["a", "b", "c", "d"]].sum().sum()

    for _, row in contingency_df.iterrows():
        a = row["a"]
        E = ((row["a"] + row["b"]) * (row["a"] + row["c"])) / N_total
        if E > 0:
            # 簡易 EBGM (full GPS は EM で混合ガンマ推定)
            ebgm = (a + shrinkage_prior) / (E + shrinkage_prior)
            # EB05 近似 (Poisson 下限)
            from scipy.stats import poisson
            eb05 = poisson.ppf(0.05, a + shrinkage_prior) / (E + shrinkage_prior)
        else:
            ebgm, eb05 = 0, 0

        results.append({
            "pt": row["pt"],
            "observed": a,
            "expected": round(E, 3),
            "ebgm": round(ebgm, 3),
            "eb05": round(eb05, 3),
            "signal_ebgm": eb05 >= 2,  # EB05 ≥ 2 がシグナル基準
        })

    return pd.DataFrame(results).sort_values("ebgm", ascending=False)

ebgm_results = calculate_ebgm(signals)
print(f"EBGM signals (EB05≥2): {ebgm_results['signal_ebgm'].sum()}")
```

### 8. 安全性シグナルレポート生成

```python
import json

def generate_pv_report(target_drug, signals_df, ebgm_df, output_dir="results"):
    """安全性シグナルレポートの統合生成"""

    # 統合シグナル判定
    combined = signals_df.merge(ebgm_df[["pt", "ebgm", "eb05", "signal_ebgm"]], on="pt")
    combined["consensus_signal"] = (
        combined["signal_prr"] & combined["signal_ror"] & combined["signal_ebgm"]
    )

    report = {
        "drug": target_drug,
        "analysis_date": pd.Timestamp.now().isoformat(),
        "total_reports": int(signals_df["a"].sum()),
        "unique_pts_analyzed": len(signals_df),
        "signals_prr": int(signals_df["signal_prr"].sum()),
        "signals_ror": int(signals_df["signal_ror"].sum()),
        "signals_ic": int(signals_df["signal_ic"].sum()),
        "signals_ebgm": int(ebgm_df["signal_ebgm"].sum()),
        "consensus_signals": int(combined["consensus_signal"].sum()),
        "top_signals": combined[combined["consensus_signal"]].nlargest(20, "prr").to_dict("records"),
        "summary": {
            "method": "PRR + ROR + IC + EBGM consensus",
            "thresholds": {
                "prr": "≥2, χ²≥4, N≥3",
                "ror": "lower 95% CI > 1, N≥3",
                "ic": "IC > 0, N≥3",
                "ebgm": "EB05 ≥ 2",
            },
        },
    }

    with open(f"{output_dir}/pv_signal_report.json", "w") as f:
        json.dump(report, f, indent=2, default=str)

    # Markdown レポート
    md = f"# Pharmacovigilance Signal Report: {target_drug}\n\n"
    md += f"**Analysis Date**: {report['analysis_date']}\n\n"
    md += f"## Summary\n\n"
    md += f"| Metric | Value |\n|---|---|\n"
    md += f"| Total Reports | {report['total_reports']:,} |\n"
    md += f"| PTs Analyzed | {report['unique_pts_analyzed']} |\n"
    md += f"| PRR Signals | {report['signals_prr']} |\n"
    md += f"| ROR Signals | {report['signals_ror']} |\n"
    md += f"| IC Signals | {report['signals_ic']} |\n"
    md += f"| EBGM Signals | {report['signals_ebgm']} |\n"
    md += f"| **Consensus Signals** | **{report['consensus_signals']}** |\n\n"
    md += f"## Top Consensus Signals\n\n"
    md += "| PT | N | PRR | ROR | IC | EBGM |\n|---|---|---|---|---|---|\n"
    for s in report["top_signals"]:
        md += f"| {s['pt']} | {s['a']} | {s['prr']} | {s['ror']} | {s['ic']} | {s['ebgm']} |\n"

    with open(f"{output_dir}/pv_signal_report.md", "w") as f:
        f.write(md)

    return report
```

---

## Best Practices

1. **複数手法のコンセンサス**: PRR/ROR/IC/EBGM を併用し、2 手法以上で一致したシグナルを優先
2. **MedDRA PT レベルで集計**: LLT は粒度が細かすぎ、SOC は粗すぎる
3. **最小報告基準 N≥3**: 少数報告の偽陽性を防ぐ
4. **Weber 効果に注意**: 新薬発売直後は報告バイアスが大きい
5. **Notoriety Bias を考慮**: メディア報道後に報告が急増する可能性
6. **重複報告の除去**: FAERS は重複が多い（約 10-15%）。Case ID ベースで最新版を保持
7. **適応症の交絡**: indication-reaction の混同に注意

## Completeness Checklist

- [ ] FAERS データ取得と前処理完了
- [ ] MedDRA 階層マッピング適用
- [ ] PRR/ROR/IC/EBGM すべて算出
- [ ] シグナル閾値判定（コンセンサス方式）
- [ ] Time-to-Onset 分布分析
- [ ] 人口統計層別化（年齢・性別・転帰）
- [ ] 因果関係評価（Naranjo / WHO-UMC）
- [ ] 安全性シグナルレポート（JSON + Markdown）生成

## References

### Output Files

| ファイル | 形式 | 生成タイミング |
|---|---|---|
| `results/pv_signal_report.json` | シグナル検出結果（JSON） | 不均衡分析完了時 |
| `results/pv_signal_report.md` | シグナルレポート（Markdown） | レポート生成時 |
| `figures/pv_temporal_trend.png` | 時系列トレンド図 | トレンド分析時 |
| `figures/pv_demographics.png` | 人口統計分布図 | 層別化分析時 |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| FAERS | `FAERS_count_reactions_by_drug_event` | 有害事象カウント |
| FAERS | `FAERS_calculate_disproportionality` | PRR/ROR/IC 不均衡分析 |
| FAERS | `FAERS_stratify_by_demographics` | 年齢・性別・国別層別化 |
| FDA | `FDA_get_adverse_reactions_by_drug_name` | 添付文書副作用情報 |
| DailyMed | `DailyMed_search_spls` | 添付文書検索 |
| DailyMed | `DailyMed_parse_adverse_reactions` | 副作用テーブル抽出 |
| PharmGKB | `PharmGKB_get_dosing_guidelines` | PGx 用量ガイドライン |
| CPIC | `CPIC_get_guidelines` | CPIC ガイドライン取得 |

### 参照スキル

| スキル | 連携 |
|---|---|
| `scientific-survival-clinical` | ← 臨床安全性解析・有害事象グレーディング |
| `scientific-statistical-testing` | ← χ² 検定・多重比較補正 |
| `scientific-drug-target-profiling` | ← 標的プロファイル・薬理学的背景 |
| `scientific-admet-pharmacokinetics` | ← 毒性予測・代謝経路情報 |
| `scientific-clinical-decision-support` | → シグナル情報の臨床意思決定への反映 |
| `scientific-deep-research` | ← 安全性文献の深層リサーチ |
| `scientific-clinical-trials-analytics` | ← 臨床試験データベース照会 |
| `scientific-regulatory-science` | → FDA/FAERS 規制データ統合 |
| `scientific-pharmacogenomics` | ← PGx 代謝型別安全性評価 |

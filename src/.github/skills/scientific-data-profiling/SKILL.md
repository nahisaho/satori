---
name: scientific-data-profiling
description: |
  データプロファイリング・品質スキル。ydata-profiling 自動 EDA ・
  Great Expectations データバリデーション・データ品質スコア・
  型推論・相関検出・外れ値フラグ・データカタログ生成。
tu_tools:
  - key: biotools
    name: bio.tools
    description: データプロファイリングツール検索
---

# Scientific Data Profiling

データセットの包括的プロファイリング・品質評価・
自動 EDA レポートパイプラインを提供する。

## When to Use

- 新しいデータセットの全体像を素早く把握するとき
- データ品質スコアを算出して品質基準をチェックするとき
- ydata-profiling で自動 EDA レポートを生成するとき
- Great Expectations でデータバリデーションルールを定義するとき
- データカタログ (辞書) を自動生成するとき
- 相関・外れ値・欠損を一括診断するとき

---

## Quick Start

## 1. ydata-profiling 自動 EDA

```python
import numpy as np
import pandas as pd


def auto_profile_report(df, title="Data Profile Report",
                        minimal=False, output="profile_report.html"):
    """
    ydata-profiling 自動 EDA レポート。

    Parameters:
        df: pd.DataFrame — 入力データ
        title: str — レポートタイトル
        minimal: bool — 軽量モード
        output: str — 出力 HTML パス
    """
    from ydata_profiling import ProfileReport

    profile = ProfileReport(
        df, title=title, minimal=minimal,
        correlations={"pearson": {"calculate": True},
                      "spearman": {"calculate": True},
                      "kendall": {"calculate": True}},
        missing_diagrams={"bar": True, "matrix": True, "heatmap": True})

    profile.to_file(output)

    # サマリー抽出
    desc = profile.get_description()
    summary = {
        "n_rows": len(df),
        "n_cols": len(df.columns),
        "n_numeric": len(df.select_dtypes(include=[np.number]).columns),
        "n_categorical": len(df.select_dtypes(include=["object", "category"]).columns),
        "total_missing": int(df.isnull().sum().sum()),
        "missing_pct": float(df.isnull().sum().sum() / (len(df) * len(df.columns)) * 100),
        "n_duplicates": int(df.duplicated().sum()),
    }

    print(f"Profile Report → {output}")
    print(f"  {summary['n_rows']} rows × {summary['n_cols']} cols, "
          f"{summary['missing_pct']:.1f}% missing, "
          f"{summary['n_duplicates']} duplicates")
    return {"report_path": output, "summary": summary}
```

## 2. データ品質スコア

```python
def data_quality_score(df, rules=None):
    """
    データ品質スコア算出 (0-100)。

    Parameters:
        df: pd.DataFrame — 入力データ
        rules: dict | None — カスタムルール
    """
    scores = {}

    # 1. 完全性 (Completeness) — 非欠損率
    completeness = 1.0 - df.isnull().sum().sum() / (len(df) * len(df.columns))
    scores["completeness"] = completeness

    # 2. 一意性 (Uniqueness) — 非重複率
    uniqueness = 1.0 - df.duplicated().sum() / len(df) if len(df) > 0 else 1.0
    scores["uniqueness"] = uniqueness

    # 3. 一貫性 (Consistency) — 型一貫性
    type_consistent = 0
    for col in df.columns:
        non_null = df[col].dropna()
        if len(non_null) == 0:
            type_consistent += 1
            continue
        try:
            inferred = pd.api.types.infer_dtype(non_null, skipna=True)
            if inferred not in ["mixed", "mixed-integer"]:
                type_consistent += 1
        except Exception:
            pass
    consistency = type_consistent / len(df.columns) if len(df.columns) > 0 else 1.0
    scores["consistency"] = consistency

    # 4. 適時性 (Timeliness) — 日付カラムの新しさ
    date_cols = df.select_dtypes(include=["datetime64"]).columns
    if len(date_cols) > 0:
        max_date = df[date_cols[0]].max()
        freshness = 1.0  # Placeholder
        scores["timeliness"] = freshness
    else:
        scores["timeliness"] = 1.0

    # 5. 妥当性 (Validity) — 数値カラムの有限性
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    if len(numeric_cols) > 0:
        finite_rate = df[numeric_cols].apply(lambda x: np.isfinite(x.dropna()).mean()).mean()
        scores["validity"] = float(finite_rate)
    else:
        scores["validity"] = 1.0

    # 総合スコア
    weights = {"completeness": 0.3, "uniqueness": 0.2,
               "consistency": 0.2, "timeliness": 0.1, "validity": 0.2}
    total_score = sum(scores[k] * weights[k] for k in weights) * 100

    # カスタムルール
    rule_results = []
    if rules:
        for rule_name, rule_fn in rules.items():
            try:
                passed = rule_fn(df)
                rule_results.append({"rule": rule_name, "passed": passed})
            except Exception as e:
                rule_results.append({"rule": rule_name, "passed": False,
                                     "error": str(e)})

    print(f"Data Quality Score: {total_score:.1f}/100")
    for k, v in scores.items():
        print(f"  {k}: {v:.3f}")

    return {"total_score": total_score, "dimension_scores": scores,
            "rule_results": rule_results}
```

## 3. Great Expectations バリデーション

```python
def great_expectations_validate(df, expectations=None):
    """
    Great Expectations スタイルのデータバリデーション。

    Parameters:
        df: pd.DataFrame — 入力データ
        expectations: list[dict] | None — バリデーションルール
    """
    if expectations is None:
        expectations = _auto_generate_expectations(df)

    results = []
    for exp in expectations:
        exp_type = exp["type"]
        col = exp.get("column")
        kwargs = exp.get("kwargs", {})

        try:
            if exp_type == "expect_column_to_exist":
                success = col in df.columns
            elif exp_type == "expect_column_values_to_not_be_null":
                max_pct = kwargs.get("mostly", 1.0)
                non_null_pct = df[col].notnull().mean()
                success = non_null_pct >= max_pct
            elif exp_type == "expect_column_values_to_be_between":
                min_val, max_val = kwargs["min_value"], kwargs["max_value"]
                vals = df[col].dropna()
                success = bool((vals >= min_val).all() and (vals <= max_val).all())
            elif exp_type == "expect_column_values_to_be_unique":
                success = not df[col].duplicated().any()
            elif exp_type == "expect_column_values_to_be_in_set":
                valid_set = set(kwargs["value_set"])
                success = df[col].dropna().isin(valid_set).all()
            elif exp_type == "expect_table_row_count_to_be_between":
                success = kwargs["min_value"] <= len(df) <= kwargs["max_value"]
            else:
                success = None

            results.append({"expectation": exp_type, "column": col,
                            "success": success})
        except Exception as e:
            results.append({"expectation": exp_type, "column": col,
                            "success": False, "error": str(e)})

    results_df = pd.DataFrame(results)
    n_pass = results_df["success"].sum()
    n_total = len(results_df)
    print(f"Validation: {n_pass}/{n_total} expectations passed "
          f"({n_pass/n_total*100:.0f}%)")
    return results_df


def _auto_generate_expectations(df):
    """自動でバリデーションルールを推論。"""
    expectations = []
    for col in df.columns:
        expectations.append({"type": "expect_column_to_exist", "column": col})
        expectations.append({
            "type": "expect_column_values_to_not_be_null",
            "column": col,
            "kwargs": {"mostly": 0.9}})

        if df[col].dtype in [np.float64, np.int64]:
            q1, q3 = df[col].quantile([0.01, 0.99])
            iqr = q3 - q1
            expectations.append({
                "type": "expect_column_values_to_be_between",
                "column": col,
                "kwargs": {"min_value": float(q1 - 3 * iqr),
                           "max_value": float(q3 + 3 * iqr)}})
    return expectations
```

---

## パイプライン統合

```
[データ取得] → data-profiling → eda-correlation
                (品質診断)       (探索的解析)
                     │                ↓
          missing-data-analysis  anomaly-detection
            (欠損補完)             (異常検知)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `profile_report.html` | ydata-profiling レポート | → EDA |
| `quality_score.json` | データ品質スコア | → 品質管理 |
| `validation_results.csv` | バリデーション結果 | → データ修正 |

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `biotools` | bio.tools | データプロファイリングツール検索 |

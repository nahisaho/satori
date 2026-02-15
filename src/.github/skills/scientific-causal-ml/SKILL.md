---
name: scientific-causal-ml
description: |
  因果機械学習スキル。DoWhy 因果モデル・EconML CATE 推定・
  Double/Debiased ML・Causal Forest・メタラーナー (S/T/X)・
  異質的処置効果 (HTE)・因果特徴量選択。
tu_tools:
  - key: openml
    name: OpenML
    description: 因果推論 ML データセット参照
---

# Scientific Causal ML

機械学習ベースの因果推論パイプラインを提供し、
異質的処置効果 (HTE) の推定と因果特徴量発見を実現する。

## When to Use

- 処置効果が個人/サブグループで異なるとき (HTE 推定)
- Causal Forest で非パラメトリック因果効果を推定するとき
- Double ML で高次元データの処置効果を推定するとき
- メタラーナー (S/T/X-learner) で CATE を推定するとき
- DoWhy で因果モデルの同定・推定・反論をするとき
- 因果特徴量選択で重要な効果修飾因子を発見するとき

> **Note**: 統計的因果推論 (PSM/IPW/DID/RDD) は `scientific-causal-inference` を参照。

---

## Quick Start

## 1. DoWhy 因果モデル

```python
import numpy as np
import pandas as pd


def dowhy_causal_model(df, treatment, outcome, common_causes,
                       effect_modifiers=None, method="backdoor.linear_regression"):
    """
    DoWhy 因果推論パイプライン (同定→推定→反論)。

    Parameters:
        df: pd.DataFrame — 観測データ
        treatment: str — 処置変数
        outcome: str — 結果変数
        common_causes: list[str] — 共変量
        effect_modifiers: list[str] | None — 効果修飾因子
        method: str — 推定手法
    """
    import dowhy

    model = dowhy.CausalModel(
        data=df,
        treatment=treatment,
        outcome=outcome,
        common_causes=common_causes,
        effect_modifiers=effect_modifiers)

    # 同定
    estimand = model.identify_effect(proceed_when_unidentifiable=True)
    print(f"Identified estimand: {estimand.get_frontdoor_variables()}")

    # 推定
    estimate = model.estimate_effect(
        estimand, method_name=method)
    print(f"ATE = {estimate.value:.4f} (95% CI: [{estimate.get_confidence_intervals()[0]:.4f}, "
          f"{estimate.get_confidence_intervals()[1]:.4f}])")

    # 反論テスト
    refutations = {}
    for refuter_name in ["random_common_cause", "placebo_treatment_refuter",
                         "data_subset_refuter"]:
        try:
            refutation = model.refute_estimate(
                estimand, estimate, method_name=refuter_name)
            refutations[refuter_name] = {
                "new_effect": float(refutation.new_effect),
                "p_value": getattr(refutation, "refutation_result", {}).get("p_value", None)
            }
        except Exception:
            pass

    print(f"Refutation tests: {len(refutations)} passed")
    return {"model": model, "estimand": estimand,
            "estimate": estimate, "refutations": refutations}
```

## 2. EconML Double ML / Causal Forest

```python
def double_ml_estimate(df, treatment, outcome, features,
                       n_splits=5, model_type="linear"):
    """
    Double/Debiased ML による処置効果推定。

    Parameters:
        df: pd.DataFrame — データ
        treatment: str — 処置変数
        outcome: str — 結果変数
        features: list[str] — 共変量
        n_splits: int — クロスフィッティング分割数
        model_type: str — "linear" / "forest"
    """
    from econml.dml import LinearDML, CausalForestDML

    Y = df[outcome].values
    T = df[treatment].values
    X = df[features].values

    if model_type == "linear":
        est = LinearDML(cv=n_splits, random_state=42)
    else:
        est = CausalForestDML(
            n_estimators=200, cv=n_splits, random_state=42)

    est.fit(Y, T, X=X)

    ate = est.ate(X)
    ate_ci = est.ate_interval(X, alpha=0.05)

    # CATE (個人レベル)
    cate = est.effect(X)
    cate_ci = est.effect_interval(X, alpha=0.05)

    result_df = pd.DataFrame(X, columns=features)
    result_df["cate"] = cate
    result_df["cate_lower"] = cate_ci[0]
    result_df["cate_upper"] = cate_ci[1]

    print(f"Double ML ({model_type}): ATE={ate:.4f} "
          f"[{ate_ci[0]:.4f}, {ate_ci[1]:.4f}]")
    print(f"  CATE range: [{cate.min():.4f}, {cate.max():.4f}]")
    return {"ate": ate, "ate_ci": ate_ci,
            "cate_df": result_df, "model": est}


def causal_forest(df, treatment, outcome, features,
                  n_estimators=500):
    """
    Causal Forest — 非パラメトリック HTE 推定。

    Parameters:
        df: pd.DataFrame — データ
        treatment: str — 処置変数 (binary)
        outcome: str — 結果変数
        features: list[str] — 共変量
        n_estimators: int — 木の数
    """
    from econml.dml import CausalForestDML

    Y = df[outcome].values
    T = df[treatment].values
    X = df[features].values

    cf = CausalForestDML(
        n_estimators=n_estimators, random_state=42,
        min_samples_leaf=10)
    cf.fit(Y, T, X=X)

    cate = cf.effect(X)
    cate_ci = cf.effect_interval(X, alpha=0.05)

    # 特徴量重要度 (因果)
    importances = cf.feature_importances_
    feat_imp = pd.DataFrame({
        "feature": features,
        "causal_importance": importances
    }).sort_values("causal_importance", ascending=False)

    print(f"Causal Forest: {n_estimators} trees, "
          f"CATE median={np.median(cate):.4f}")
    print(f"  Top causal features: {feat_imp.head(5).to_dict('records')}")
    return {"cate": cate, "cate_ci": cate_ci,
            "feature_importance": feat_imp, "model": cf}
```

## 3. メタラーナー (S/T/X-Learner)

```python
def meta_learner(df, treatment, outcome, features,
                 learner_type="t", base_model=None):
    """
    メタラーナーによる CATE 推定。

    Parameters:
        df: pd.DataFrame — データ
        treatment: str — 処置変数 (binary 0/1)
        outcome: str — 結果変数
        features: list[str] — 共変量
        learner_type: str — "s" / "t" / "x"
        base_model: BaseEstimator | None — ベースモデル
    """
    from econml.metalearners import SLearner, TLearner, XLearner
    from sklearn.ensemble import GradientBoostingRegressor

    if base_model is None:
        base_model = GradientBoostingRegressor(
            n_estimators=200, max_depth=5, random_state=42)

    Y = df[outcome].values
    T = df[treatment].values
    X = df[features].values

    learners = {"s": SLearner, "t": TLearner, "x": XLearner}
    LearnerClass = learners[learner_type]

    if learner_type == "s":
        est = LearnerClass(overall_model=base_model)
    else:
        est = LearnerClass(models=base_model)

    est.fit(Y, T, X=X)
    cate = est.effect(X)

    result_df = pd.DataFrame(X, columns=features)
    result_df["cate"] = cate

    print(f"{learner_type.upper()}-Learner: "
          f"CATE mean={cate.mean():.4f}, std={cate.std():.4f}")
    return {"cate": cate, "cate_df": result_df, "model": est}
```

---

## パイプライン統合

```
causal-inference → causal-ml → feature-importance
   (統計的因果)    (因果 ML)     (特徴量解釈)
        │               │              ↓
  clinical-trial ───────┘      explainable-ai
   (臨床試験)                   (説明可能 AI)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `dowhy_causal_model.json` | DoWhy 因果モデル | → reporting |
| `cate_estimates.csv` | CATE 推定値 | → precision-medicine |
| `causal_feature_importance.csv` | 因果特徴量重要度 | → explainable-ai |

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `openml` | OpenML | 因果推論 ML データセット参照 |

---
name: scientific-explainable-ai
description: |
  説明可能 AI (XAI) スキル。SHAP・LIME・Captum・InterpretML を活用し、
  モデル予測の根拠説明・特徴量寄与分解・反実仮想説明・公平性監査を支援。
  「モデルの予測を説明して」「SHAP 値を計算して」「LIME で説明して」で発火。
tu_tools:
  - key: papers_with_code
    name: Papers with Code
    description: XAI 手法・ベンチマーク検索
---

# Scientific Explainable AI

説明可能 AI (Explainable AI / XAI) のための解析スキル。
ブラックボックスモデルの予測根拠を可視化・定量化し、
科学的知見の抽出と信頼性検証を支援する。

## When to Use

- ML/DL モデル予測の根拠説明
- SHAP 値による特徴量寄与の定量的分解
- LIME によるローカル説明生成
- 反実仮想説明（Counterfactual Explanation）
- モデル公平性・バイアス監査
- 規制対応の説明責任（FDA/EMA AI ガイドライン）

## Quick Start

### XAI パイプライン

```
Phase 1: Model Training
  - ベースモデルの学習
  - テストセット予測結果
    ↓
Phase 2: Global Explanation
  - SHAP summary plot (全体的特徴量重要度)
  - Permutation importance
  - Partial Dependence Plot (PDP)
    ↓
Phase 3: Local Explanation
  - SHAP force/waterfall plot (個別予測の説明)
  - LIME 説明 (局所線形近似)
  - Anchor 説明 (十分条件ルール)
    ↓
Phase 4: Interaction Analysis
  - SHAP interaction values
  - 特徴量間の依存関係
  - 非線形効果の可視化
    ↓
Phase 5: Counterfactual Analysis
  - What-if シナリオ
  - 反実仮想最近傍
  - Action recommendation
    ↓
Phase 6: Fairness & Audit
  - サブグループ間の公平性
  - バイアス検出
  - 説明一貫性検証
```

## Workflow

### 1. SHAP (SHapley Additive exPlanations)

```python
import shap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import train_test_split

# === モデル学習 ===
# 例: 分類モデル
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
model = GradientBoostingClassifier(n_estimators=200, random_state=42)
model.fit(X_train, y_train)

# === TreeSHAP (ツリーモデル用高速版) ===
explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(X_test)

# Global: Summary Plot
fig, ax = plt.subplots(figsize=(10, 8))
shap.summary_plot(shap_values, X_test, feature_names=feature_names,
                  show=False)
plt.tight_layout()
plt.savefig("figures/shap_summary.png", dpi=300, bbox_inches="tight")
plt.show()

# Global: Bar Plot (平均 |SHAP|)
shap.summary_plot(shap_values, X_test, feature_names=feature_names,
                  plot_type="bar", show=False)
plt.savefig("figures/shap_importance.png", dpi=300, bbox_inches="tight")
plt.show()

# Local: Waterfall Plot (個別サンプル)
shap.waterfall_plot(shap.Explanation(
    values=shap_values[0],
    base_values=explainer.expected_value,
    data=X_test.iloc[0],
    feature_names=feature_names,
))
plt.savefig("figures/shap_waterfall.png", dpi=300, bbox_inches="tight")
plt.show()
```

### 2. SHAP Interaction Values

```python
# === SHAP Interaction Effects ===
shap_interaction = explainer.shap_interaction_values(X_test)

# Interaction heatmap
interaction_matrix = np.abs(shap_interaction).mean(axis=0)
fig, ax = plt.subplots(figsize=(10, 8))
im = ax.imshow(interaction_matrix, cmap="YlOrRd")
ax.set_xticks(range(len(feature_names)))
ax.set_yticks(range(len(feature_names)))
ax.set_xticklabels(feature_names, rotation=45, ha="right")
ax.set_yticklabels(feature_names)
plt.colorbar(im, label="Mean |SHAP Interaction|")
ax.set_title("Feature Interaction Heatmap")
plt.tight_layout()
plt.savefig("figures/shap_interaction.png", dpi=300, bbox_inches="tight")
plt.show()

# Dependence plot with interaction
shap.dependence_plot(
    "feature_A", shap_values, X_test,
    interaction_index="feature_B",
    feature_names=feature_names,
)
```

### 3. LIME (Local Interpretable Model-agnostic Explanations)

```python
from lime.lime_tabular import LimeTabularExplainer

# === LIME Explainer ===
lime_explainer = LimeTabularExplainer(
    training_data=X_train.values,
    feature_names=feature_names,
    class_names=["Negative", "Positive"],
    mode="classification",
    discretize_continuous=True,
)

# 個別予測の説明
def explain_instance_lime(idx, num_features=10):
    """LIME で 1 インスタンスの予測を説明"""
    exp = lime_explainer.explain_instance(
        X_test.iloc[idx].values,
        model.predict_proba,
        num_features=num_features,
        num_samples=5000,
    )

    # 説明テーブル
    explanation_df = pd.DataFrame(
        exp.as_list(),
        columns=["Feature Rule", "Weight"]
    )
    explanation_df["Abs_Weight"] = explanation_df["Weight"].abs()
    explanation_df = explanation_df.sort_values("Abs_Weight", ascending=False)

    # R² (ローカル忠実度)
    local_fidelity = exp.score

    return {
        "prediction": model.predict_proba(X_test.iloc[idx:idx+1])[0],
        "local_fidelity_r2": round(local_fidelity, 4),
        "explanation": explanation_df,
    }

result = explain_instance_lime(0)
print(f"Local fidelity (R²): {result['local_fidelity_r2']}")
print(result["explanation"])
```

### 4. DeepSHAP / GradientSHAP (Deep Learning)

```python
import torch

# === DeepSHAP for Neural Networks ===
def deepshap_explain(model, background_data, test_data, feature_names):
    """Deep SHAP: DeepLIFT + SHAP"""
    background = torch.tensor(background_data[:100].values, dtype=torch.float32)
    test_tensor = torch.tensor(test_data.values, dtype=torch.float32)

    explainer = shap.DeepExplainer(model, background)
    shap_values = explainer.shap_values(test_tensor)

    return shap_values


# === GradientSHAP ===
def gradient_shap_explain(model, background_data, test_data):
    """Gradient SHAP: Integrated Gradients + SHAP"""
    explainer = shap.GradientExplainer(model, background_data)
    shap_values = explainer.shap_values(test_data, nsamples=200)
    return shap_values


# === Captum (PyTorch 公式 XAI) ===
from captum.attr import IntegratedGradients, LayerConductance, NeuronConductance

def captum_integrated_gradients(model, inputs, target_class=None):
    """Integrated Gradients for PyTorch model"""
    ig = IntegratedGradients(model)
    attributions = ig.attribute(inputs, target=target_class, n_steps=200)
    return attributions
```

### 5. 反実仮想説明 (Counterfactual Explanation)

```python
def counterfactual_explanation(model, instance, desired_class,
                                feature_names, feature_ranges,
                                n_cf=5, max_iter=1000):
    """
    反実仮想説明: 「何が変われば予測が変わるか」
    最小の特徴量変更で予測クラスを反転させるインスタンスを探索
    """
    from scipy.optimize import differential_evolution

    original_pred = model.predict(instance.reshape(1, -1))[0]

    def objective(x):
        # 予測クラスが desired_class であること
        pred_proba = model.predict_proba(x.reshape(1, -1))[0]
        class_loss = -pred_proba[desired_class]  # 目標クラスの確率を最大化

        # 変更量を最小化 (L1 距離)
        change_loss = np.sum(np.abs(x - instance) / (np.array(feature_ranges) + 1e-8))

        # スパース性ペナルティ
        sparsity_loss = np.sum(np.abs(x - instance) > 1e-3) * 0.1

        return class_loss + 0.5 * change_loss + sparsity_loss

    bounds = [(r[0], r[1]) for r in feature_ranges]
    result = differential_evolution(objective, bounds, maxiter=max_iter, seed=42)

    cf_instance = result.x
    cf_pred = model.predict(cf_instance.reshape(1, -1))[0]

    # 変更された特徴量
    changes = []
    for i, (orig, cf) in enumerate(zip(instance, cf_instance)):
        if abs(orig - cf) > 1e-3:
            changes.append({
                "feature": feature_names[i],
                "original": round(orig, 4),
                "counterfactual": round(cf, 4),
                "change": round(cf - orig, 4),
            })

    return {
        "original_prediction": int(original_pred),
        "counterfactual_prediction": int(cf_pred),
        "changes": changes,
        "n_features_changed": len(changes),
    }
```

### 6. 公平性監査

```python
def fairness_audit(model, X_test, y_test, sensitive_feature,
                   feature_names):
    """
    モデル公平性監査: サブグループ間の予測バイアスを検出
    """
    from sklearn.metrics import accuracy_score, recall_score, precision_score

    predictions = model.predict(X_test)
    proba = model.predict_proba(X_test)[:, 1]

    groups = X_test[sensitive_feature].unique()
    metrics = {}

    for group in groups:
        mask = X_test[sensitive_feature] == group
        metrics[group] = {
            "n": int(mask.sum()),
            "accuracy": round(accuracy_score(y_test[mask], predictions[mask]), 4),
            "recall": round(recall_score(y_test[mask], predictions[mask]), 4),
            "precision": round(precision_score(y_test[mask], predictions[mask], zero_division=0), 4),
            "positive_rate": round(predictions[mask].mean(), 4),
            "avg_score": round(proba[mask].mean(), 4),
        }

    # Demographic Parity Difference
    rates = [m["positive_rate"] for m in metrics.values()]
    dp_diff = max(rates) - min(rates)

    # Equalized Odds Difference
    recalls = [m["recall"] for m in metrics.values()]
    eo_diff = max(recalls) - min(recalls)

    audit_result = {
        "group_metrics": metrics,
        "demographic_parity_diff": round(dp_diff, 4),
        "equalized_odds_diff": round(eo_diff, 4),
        "fair_threshold": 0.1,  # < 0.1 is considered fair
        "is_fair_dp": dp_diff < 0.1,
        "is_fair_eo": eo_diff < 0.1,
    }

    return audit_result
```

---

## Best Practices

1. **グローバル + ローカルの併用**: SHAP summary (全体) + waterfall (個別) を両方報告
2. **SHAP > LIME**: SHAP は理論的保証 (一貫性・局所正確性) があり優先
3. **背景データの選択**: KernelSHAP では背景データ数が結果に影響 (100-1000 推奨)
4. **LIME の忠実度確認**: `exp.score` (R²) が 0.8 以上であることを確認
5. **反実仮想のアクション可能性**: 変更不可能な特徴（年齢・性別）は固定
6. **モデル非依存 vs モデル特化**: TreeSHAP はツリー用、DeepSHAP は NN 用
7. **規制準拠**: 医療 AI では FDA/EMA の説明可能性要件に対応

## Completeness Checklist

- [ ] SHAP summary plot (グローバル重要度)
- [ ] SHAP waterfall/force plot (ローカル説明)
- [ ] SHAP interaction values (特徴量相互作用)
- [ ] LIME 説明 (ローカル忠実度 R² 確認)
- [ ] PDP (部分依存プロット)
- [ ] 反実仮想説明（必要に応じて）
- [ ] 公平性監査（機密属性がある場合）
- [ ] 説明レポート生成

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `papers_with_code` | Papers with Code | XAI 手法・ベンチマーク検索 |

## References

### Output Files

| ファイル | 形式 | 生成タイミング |
|---|---|---|
| `results/xai_report.json` | XAI 結果サマリー（JSON） | 説明分析完了時 |
| `figures/shap_summary.png` | SHAP summary plot | グローバル説明時 |
| `figures/shap_waterfall.png` | SHAP waterfall plot | ローカル説明時 |
| `figures/shap_interaction.png` | 相互作用ヒートマップ | Interaction 分析時 |

### 参照スキル

| スキル | 連携 |
|---|---|
| `scientific-feature-importance` | ← Permutation/Tree importance との比較 |
| `scientific-ml-classification` | ← 分類モデルの説明 |
| `scientific-ml-regression` | ← 回帰モデルの説明 |
| `scientific-deep-learning` | ← DNN モデルの DeepSHAP/Captum 説明 |
| `scientific-graph-neural-networks` | ← GNN モデルの GNNExplainer 説明 |
| `scientific-clinical-decision-support` | → 臨床 AI の説明責任 |

---
name: scientific-multi-omics
description: |
  マルチオミクス統合解析スキル。ゲノム・トランスクリプトーム・プロテオーム・メタボローム
  データの統合手法（MOFA/SNF/DIABLO）、オミクス間相関解析、CCA/PLS 統合、
  パスウェイレベル統合、ネットワーク統合のテンプレートを提供。
---

# Scientific Multi-Omics Integration

複数のオミクスレイヤー（ゲノミクス、トランスクリプトミクス、プロテオミクス、
メタボロミクス）のデータを統合的に解析するためのスキル。各オミクスの個別解析は
それぞれの専門スキル（bioinformatics, metabolomics）に委ねつつ、本スキルは
**統合**に特化する。

## When to Use

- 同一サンプルから得られた複数オミクスデータを統合するとき
- オミクス間の相関構造を解明したいとき
- マルチオミクスバイオマーカー発見が必要なとき
- パスウェイレベルでの統合解析が必要なとき

---

## Quick Start

## 1. データ整合性チェック

```python
import numpy as np
import pandas as pd

def check_multiomics_alignment(omics_dict, sample_col="Sample_ID"):
    """
    マルチオミクスデータの整合性チェック。

    Parameters:
        omics_dict: {"transcriptomics": df, "proteomics": df, "metabolomics": df}
    """
    all_samples = {}
    for name, df in omics_dict.items():
        all_samples[name] = set(df[sample_col])

    # 共通サンプル
    common = set.intersection(*all_samples.values())
    report = {
        "n_omics_layers": len(omics_dict),
        "layers": list(omics_dict.keys()),
        "samples_per_layer": {k: len(v) for k, v in all_samples.items()},
        "common_samples": len(common),
        "features_per_layer": {k: df.shape[1] - 1 for k, df in omics_dict.items()},
    }

    # 共通サンプルでフィルタ
    aligned = {}
    for name, df in omics_dict.items():
        aligned[name] = df[df[sample_col].isin(common)].sort_values(sample_col).reset_index(drop=True)

    print(f"=== Multi-Omics Alignment ===")
    print(f"Common samples: {len(common)} / {max(report['samples_per_layer'].values())}")
    for k, v in report["features_per_layer"].items():
        print(f"  {k}: {v} features")

    return aligned, report
```

## 2. オミクス間相関解析

```python
from scipy.stats import spearmanr

def cross_omics_correlation(omics1_df, omics2_df, name1="Omics1", name2="Omics2",
                              top_n=50, method="spearman"):
    """
    2 つのオミクスレイヤー間の特徴量ペアワイズ相関を計算する。

    Parameters:
        omics1_df: DataFrame (samples × features), no sample_col
        omics2_df: DataFrame (samples × features), no sample_col
        top_n: 上位の強い相関ペアを返す数
    """
    features1 = omics1_df.columns.tolist()
    features2 = omics2_df.columns.tolist()

    correlations = []
    for f1 in features1:
        for f2 in features2:
            if method == "spearman":
                r, p = spearmanr(omics1_df[f1], omics2_df[f2])
            else:
                from scipy.stats import pearsonr
                r, p = pearsonr(omics1_df[f1], omics2_df[f2])
            correlations.append({
                f"{name1}_feature": f1,
                f"{name2}_feature": f2,
                "correlation": r,
                "p_value": p,
                "abs_correlation": abs(r),
            })

    corr_df = pd.DataFrame(correlations).sort_values("abs_correlation", ascending=False)

    return corr_df.head(top_n)
```

## 3. CCA (Canonical Correlation Analysis)

```python
from sklearn.cross_decomposition import CCA

def canonical_correlation_analysis(X1, X2, n_components=2):
    """
    正準相関分析: 2 つのオミクスデータ間の最大相関方向を見つける。

    Returns:
        cca_model, scores_X1, scores_X2, canonical_correlations
    """
    cca = CCA(n_components=n_components)
    scores_X1, scores_X2 = cca.fit_transform(X1, X2)

    # 正準相関係数
    canonical_corrs = []
    for i in range(n_components):
        r = np.corrcoef(scores_X1[:, i], scores_X2[:, i])[0, 1]
        canonical_corrs.append(r)

    return cca, scores_X1, scores_X2, canonical_corrs


def plot_cca_scores(scores_X1, scores_X2, labels, name1="Omics1",
                     name2="Omics2", figsize=(12, 5)):
    """CCA スコアプロットを描画する。"""
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    unique_labels = np.unique(labels)
    colors = plt.cm.Set1(np.linspace(0, 0.5, len(unique_labels)))

    for color, label in zip(colors, unique_labels):
        mask = labels == label
        ax1.scatter(scores_X1[mask, 0], scores_X1[mask, 1],
                   c=[color], label=label, alpha=0.7, edgecolors="black")
        ax2.scatter(scores_X2[mask, 0], scores_X2[mask, 1],
                   c=[color], label=label, alpha=0.7, edgecolors="black")

    ax1.set_title(f"CCA — {name1}", fontweight="bold")
    ax1.set_xlabel("CC1"); ax1.set_ylabel("CC2")
    ax2.set_title(f"CCA — {name2}", fontweight="bold")
    ax2.set_xlabel("CC1"); ax2.set_ylabel("CC2")
    ax1.legend(); ax2.legend()
    plt.tight_layout()
    plt.savefig("figures/cca_scores.png", dpi=300, bbox_inches="tight")
    plt.close()
```

## 4. SNF (Similarity Network Fusion)

```python
def similarity_network_fusion(omics_list, k_neighbors=20, n_iterations=20):
    """
    Similarity Network Fusion — 複数オミクスの類似度ネットワークを融合する。

    Wang et al., Nature Methods 2014

    Parameters:
        omics_list: list of np.arrays [(n_samples, p1), (n_samples, p2), ...]
        k_neighbors: KNN のK
        n_iterations: 融合反復回数

    Returns:
        fused_similarity: (n_samples, n_samples) 融合類似度行列
    """
    from sklearn.metrics import pairwise_distances

    n = omics_list[0].shape[0]
    n_views = len(omics_list)

    # 各ビューの類似度行列
    similarities = []
    for X in omics_list:
        D = pairwise_distances(X, metric="euclidean")
        mu = np.mean(np.sort(D, axis=1)[:, 1:k_neighbors+1], axis=1)
        S = np.exp(-D**2 / (mu[:, None] * mu[None, :] + 1e-10))
        np.fill_diagonal(S, 0)
        # 正規化
        S = S / (S.sum(axis=1, keepdims=True) + 1e-10)
        similarities.append(S)

    # KNN マスク
    knn_masks = []
    for X in omics_list:
        D = pairwise_distances(X, metric="euclidean")
        mask = np.zeros_like(D, dtype=bool)
        for i in range(n):
            nn = np.argsort(D[i])[:k_neighbors+1]
            mask[i, nn] = True
        knn_masks.append(mask)

    # 反復融合
    P = [s.copy() for s in similarities]
    for _ in range(n_iterations):
        P_new = []
        for v in range(n_views):
            # 他のビューの平均
            other_avg = np.mean([P[j] for j in range(n_views) if j != v], axis=0)
            # 局所構造の保持
            S_local = similarities[v] * knn_masks[v]
            S_local = S_local / (S_local.sum(axis=1, keepdims=True) + 1e-10)
            P_updated = S_local @ other_avg @ S_local.T
            P_updated = P_updated / (P_updated.sum(axis=1, keepdims=True) + 1e-10)
            P_new.append(P_updated)
        P = P_new

    fused = np.mean(P, axis=0)
    fused = (fused + fused.T) / 2

    return fused
```

## 5. パスウェイレベル統合

```python
def pathway_level_integration(omics_dict, pathway_mapping, pathway_col="Pathway",
                                feature_col="Feature"):
    """
    パスウェイレベルでオミクスデータを統合する。
    各パスウェイの活性スコアを PCA 第一主成分で算出。

    Parameters:
        omics_dict: {"transcriptomics": df, "proteomics": df}
        pathway_mapping: DataFrame (Feature, Pathway, Omics_Layer)
    """
    from sklearn.decomposition import PCA

    pathway_scores = {}
    pathways = pathway_mapping[pathway_col].unique()

    for pw in pathways:
        pw_features = pathway_mapping[pathway_mapping[pathway_col] == pw]
        combined_data = []

        for omics_name, df in omics_dict.items():
            features_in_omics = pw_features[pw_features["Omics_Layer"] == omics_name][feature_col]
            available = [f for f in features_in_omics if f in df.columns]
            if available:
                combined_data.append(df[available].values)

        if combined_data:
            X = np.hstack(combined_data)
            if X.shape[1] >= 2:
                pca = PCA(n_components=1)
                score = pca.fit_transform(X).ravel()
                pathway_scores[pw] = {
                    "score": score,
                    "explained_variance": pca.explained_variance_ratio_[0],
                    "n_features": X.shape[1],
                }

    return pathway_scores
```

## 6. マルチオミクスクラスタリング

```python
from sklearn.cluster import SpectralClustering

def multiomics_clustering(fused_similarity, n_clusters, labels_true=None):
    """融合類似度行列に基づくスペクトラルクラスタリング。"""
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

    clustering = SpectralClustering(n_clusters=n_clusters, affinity="precomputed",
                                    random_state=42)
    cluster_labels = clustering.fit_predict(fused_similarity)

    metrics = {"n_clusters": n_clusters}
    if labels_true is not None:
        metrics["ARI"] = adjusted_rand_score(labels_true, cluster_labels)
        metrics["NMI"] = normalized_mutual_info_score(labels_true, cluster_labels)

    return cluster_labels, metrics
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/cross_omics_correlation.csv` | CSV |
| `results/canonical_correlations.csv` | CSV |
| `results/pathway_activity_scores.csv` | CSV |
| `results/multiomics_clusters.csv` | CSV |
| `figures/cca_scores.png` | PNG |
| `figures/snf_heatmap.png` | PNG |
| `figures/multiomics_umap.png` | PNG |

#### 依存パッケージ

```
scikit-learn>=1.3
scipy>=1.10
```

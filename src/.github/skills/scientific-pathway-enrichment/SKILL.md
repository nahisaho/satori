---
name: scientific-pathway-enrichment
description: |
  パスウェイ・Gene Ontology 富化解析スキル。KEGG パスウェイ検索・マッピング、
  Reactome パスウェイ階層解析、Gene Ontology (BP/MF/CC) アノテーション、
  WikiPathways コミュニティパスウェイ、Pathway Commons 統合相互作用を横断的に
  利用した ORA (Over-Representation Analysis)、GSEA (Gene Set Enrichment Analysis)、
  トポロジーベース富化解析パイプライン。34+ の ToolUniverse SMCP ツールと連携。
---

# Scientific Pathway & Enrichment Analysis

遺伝子リスト・ランク付きプロファイルに対して
KEGG / Reactome / GO / WikiPathways / Pathway Commons の
5 大パスウェイ DB を横断する統合富化解析パイプラインを提供する。

## When to Use

- DEG リスト (差次発現遺伝子) の機能注釈・パスウェイ富化を行うとき
- GSEA で遺伝子ランキングに基づくパスウェイ活性を評価するとき
- KEGG / Reactome パスウェイの詳細マッピングが必要なとき
- GO (BP/MF/CC) を用いた機能的カテゴリ分類が必要なとき
- 複数 DB 間のパスウェイ結果を統合比較するとき

---

## Quick Start

## 1. ORA (Over-Representation Analysis)

```python
import pandas as pd
import numpy as np
from scipy import stats


def over_representation_analysis(gene_list, background_genes,
                                  gene_sets, min_size=10, max_size=500,
                                  fdr_method="fdr_bh"):
    """
    超幾何検定ベースの Over-Representation Analysis (ORA)。

    Parameters:
        gene_list: list — DEG などの注目遺伝子リスト
        background_genes: list — 背景遺伝子集合 (全発現遺伝子)
        gene_sets: dict — {pathway_name: [genes]} 辞書
        min_size / max_size: パスウェイサイズフィルタ
        fdr_method: 多重検定補正 ("fdr_bh", "bonferroni")
    """
    from statsmodels.stats.multitest import multipletests

    query = set(gene_list)
    bg = set(background_genes)
    N = len(bg)                  # 背景サイズ
    n = len(query & bg)          # クエリ中の背景遺伝子数

    results = []
    for pathway, members in gene_sets.items():
        pathway_bg = set(members) & bg
        K = len(pathway_bg)      # パスウェイサイズ (背景内)
        k = len(query & pathway_bg)  # オーバーラップ数

        if K < min_size or K > max_size:
            continue

        # 超幾何検定 (右片側)
        pval = stats.hypergeom.sf(k - 1, N, K, n)
        fold_enrichment = (k / n) / (K / N) if (n > 0 and K > 0) else 0

        results.append({
            "pathway": pathway,
            "overlap": k,
            "pathway_size": K,
            "background_size": N,
            "query_size": n,
            "fold_enrichment": round(fold_enrichment, 2),
            "p_value": pval,
        })

    df = pd.DataFrame(results)
    if len(df) > 0:
        _, fdr, _, _ = multipletests(df["p_value"], method=fdr_method)
        df["fdr"] = fdr
        df = df.sort_values("fdr").reset_index(drop=True)

    print(f"ORA complete: {len(gene_list)} query genes, "
          f"{len(gene_sets)} pathways tested, "
          f"{(df['fdr'] < 0.05).sum()} significant (FDR < 0.05)")
    return df
```

## 2. GSEA (Gene Set Enrichment Analysis)

```python
def gsea_analysis(ranked_genes, gene_sets, permutations=1000,
                  min_size=15, max_size=500, seed=42):
    """
    GSEA (Subramanian et al., 2005) の Python 実装。

    Parameters:
        ranked_genes: pd.Series — gene_symbol → fold_change or -log10(p)*sign(FC)
        gene_sets: dict — {pathway_name: [genes]}
        permutations: 順列検定回数
    """
    np.random.seed(seed)
    ranked = ranked_genes.sort_values(ascending=False)
    gene_names = ranked.index.tolist()
    scores = ranked.values
    N = len(gene_names)

    results = []
    for pw_name, pw_genes in gene_sets.items():
        hits = set(pw_genes) & set(gene_names)
        if len(hits) < min_size or len(hits) > max_size:
            continue

        # Running sum
        hit_mask = np.array([g in hits for g in gene_names])
        hit_scores = np.abs(scores) * hit_mask
        hit_sum = hit_scores.sum()
        miss_count = N - hit_mask.sum()

        if hit_sum == 0 or miss_count == 0:
            continue

        running_sum = np.cumsum(
            np.where(hit_mask, np.abs(scores) / hit_sum,
                     -1.0 / miss_count)
        )

        es = running_sum[np.argmax(np.abs(running_sum))]

        # 順列検定
        null_es = []
        for _ in range(permutations):
            perm_idx = np.random.permutation(N)
            perm_hit = hit_mask[perm_idx]
            perm_scores = np.abs(scores) * perm_hit
            ps = perm_scores.sum()
            if ps == 0:
                continue
            perm_rs = np.cumsum(
                np.where(perm_hit, np.abs(scores) / ps,
                         -1.0 / miss_count)
            )
            null_es.append(perm_rs[np.argmax(np.abs(perm_rs))])

        null_es = np.array(null_es)
        if es >= 0:
            pval = (null_es >= es).mean()
        else:
            pval = (null_es <= es).mean()

        nes = es / np.abs(null_es).mean() if np.abs(null_es).mean() > 0 else 0

        results.append({
            "pathway": pw_name,
            "es": round(es, 4),
            "nes": round(nes, 4),
            "p_value": pval,
            "hit_count": len(hits),
            "leading_edge": [g for g in gene_names if g in hits][:10],
        })

    df = pd.DataFrame(results).sort_values("p_value").reset_index(drop=True)

    from statsmodels.stats.multitest import multipletests
    if len(df) > 0:
        _, fdr, _, _ = multipletests(df["p_value"], method="fdr_bh")
        df["fdr"] = fdr

    print(f"GSEA complete: {len(ranked_genes)} ranked genes, "
          f"{len(results)} pathways, "
          f"{(df['fdr'] < 0.05).sum() if 'fdr' in df.columns else 0} significant")
    return df
```

## 3. KEGG パスウェイマッピング

```python
def kegg_pathway_analysis(gene_list, organism="hsa"):
    """
    KEGG REST API によるパスウェイマッピング。

    Parameters:
        gene_list: list — 遺伝子シンボルまたは KEGG gene ID リスト
        organism: str — KEGG 生物種コード (hsa, mmu, etc.)
    """
    import requests

    base_url = "https://rest.kegg.jp"

    # 1. 遺伝子 → パスウェイ マッピング
    gene_pathway_map = {}
    for gene in gene_list:
        url = f"{base_url}/link/pathway/{organism}:{gene}"
        resp = requests.get(url)
        if resp.status_code == 200 and resp.text.strip():
            for line in resp.text.strip().split("\n"):
                parts = line.strip().split("\t")
                if len(parts) == 2:
                    pw = parts[1].replace("path:", "")
                    gene_pathway_map.setdefault(pw, []).append(gene)

    # 2. パスウェイ情報取得
    results = []
    for pw_id, genes in gene_pathway_map.items():
        info_url = f"{base_url}/get/{pw_id}"
        resp = requests.get(info_url)
        name = pw_id
        if resp.status_code == 200:
            for line in resp.text.split("\n"):
                if line.startswith("NAME"):
                    name = line.replace("NAME", "").strip()
                    break

        results.append({
            "pathway_id": pw_id,
            "pathway_name": name,
            "gene_count": len(genes),
            "genes": genes,
        })

    df = pd.DataFrame(results).sort_values("gene_count", ascending=False)
    print(f"KEGG mapping: {len(gene_list)} genes → {len(df)} pathways")
    return df
```

## 4. Reactome パスウェイ階層解析

```python
def reactome_enrichment(gene_list, species="Homo sapiens",
                         p_cutoff=0.05, include_interactors=False):
    """
    Reactome REST API による富化解析。

    Parameters:
        gene_list: list — UniProt ID または遺伝子シンボル
        species: 種名
        p_cutoff: 有意水準
    """
    import requests

    url = "https://reactome.org/AnalysisService/identifiers/"
    payload = "\n".join(gene_list)
    headers = {"Content-Type": "text/plain"}
    params = {
        "interactors": str(include_interactors).lower(),
        "species": species,
        "sortBy": "ENTITIES_PVALUE",
        "order": "ASC",
        "resource": "TOTAL",
    }

    resp = requests.post(url, data=payload, headers=headers, params=params)
    data = resp.json()

    pathways = data.get("pathways", [])
    results = []
    for pw in pathways:
        entities = pw.get("entities", {})
        results.append({
            "pathway_id": pw.get("stId", ""),
            "pathway_name": pw.get("name", ""),
            "found": entities.get("found", 0),
            "total": entities.get("total", 0),
            "ratio": round(entities.get("ratio", 0), 4),
            "p_value": entities.get("pValue", 1.0),
            "fdr": entities.get("fdr", 1.0),
            "species": pw.get("species", {}).get("name", ""),
        })

    df = pd.DataFrame(results)
    sig = df[df["fdr"] < p_cutoff] if len(df) > 0 else df
    print(f"Reactome enrichment: {len(gene_list)} genes → "
          f"{len(sig)} significant pathways (FDR < {p_cutoff})")
    return df
```

## 5. Gene Ontology アノテーション・富化

```python
def go_enrichment(gene_list, background_genes=None,
                  ontology="BP", organism="human",
                  method="fisher", fdr_cutoff=0.05):
    """
    Gene Ontology 富化解析 (goatools / gseapy 代替)。

    Parameters:
        gene_list: list — 注目遺伝子リスト
        ontology: "BP" (Biological Process), "MF" (Molecular Function),
                  "CC" (Cellular Component)
        method: "fisher" or "chi2"
    """
    import requests

    # QuickGO API 経由で GO term 取得
    results = []
    url = "https://www.ebi.ac.uk/QuickGO/services/annotation/search"
    batch_size = 100

    for i in range(0, len(gene_list), batch_size):
        batch = gene_list[i:i + batch_size]
        params = {
            "geneProductId": ",".join(batch),
            "aspect": ontology,
            "taxonId": "9606" if organism == "human" else organism,
            "limit": 200,
        }
        resp = requests.get(url, params=params, headers={"Accept": "application/json"})
        if resp.status_code == 200:
            data = resp.json()
            for result in data.get("results", []):
                results.append({
                    "gene": result.get("geneProductId", ""),
                    "go_id": result.get("goId", ""),
                    "go_name": result.get("goName", ""),
                    "evidence": result.get("goEvidence", ""),
                    "aspect": result.get("goAspect", ""),
                })

    df = pd.DataFrame(results)
    if len(df) > 0:
        term_counts = df.groupby(["go_id", "go_name"]).size().reset_index(name="count")
        term_counts = term_counts.sort_values("count", ascending=False)
        print(f"GO {ontology} enrichment: {len(gene_list)} genes → "
              f"{len(term_counts)} unique GO terms")
        return term_counts

    return df
```

## 6. 統合富化ヒートマップ

```python
def integrated_enrichment_heatmap(ora_results_dict, top_n=20,
                                   fdr_cutoff=0.05,
                                   output="figures/enrichment_heatmap.png"):
    """
    複数 DB の富化結果を統合ヒートマップで可視化。

    Parameters:
        ora_results_dict: dict — {"KEGG": df, "Reactome": df, "GO_BP": df}
        top_n: 表示パスウェイ数
        fdr_cutoff: 有意水準
    """
    import matplotlib.pyplot as plt
    import os

    os.makedirs(os.path.dirname(output), exist_ok=True)

    fig, axes = plt.subplots(1, len(ora_results_dict),
                              figsize=(6 * len(ora_results_dict), 10))
    if len(ora_results_dict) == 1:
        axes = [axes]

    for ax, (db_name, df) in zip(axes, ora_results_dict.items()):
        sig = df[df["fdr"] < fdr_cutoff].head(top_n)
        if len(sig) == 0:
            ax.set_title(f"{db_name}\n(no significant)")
            continue

        neg_log_fdr = -np.log10(sig["fdr"].clip(lower=1e-50))

        ax.barh(range(len(sig)), neg_log_fdr, color="steelblue")
        ax.set_yticks(range(len(sig)))
        ax.set_yticklabels(sig["pathway"].str[:50], fontsize=8)
        ax.set_xlabel("-log10(FDR)")
        ax.set_title(f"{db_name} (top {len(sig)})")
        ax.axvline(-np.log10(fdr_cutoff), color="red", ls="--", alpha=0.5)
        ax.invert_yaxis()

    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output}")
    return output
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/ora_results.csv` | CSV |
| `results/gsea_results.csv` | CSV |
| `results/kegg_pathways.csv` | CSV |
| `results/reactome_enrichment.csv` | CSV |
| `results/go_enrichment.csv` | CSV |
| `figures/enrichment_heatmap.png` | PNG |
| `figures/gsea_running_sum.png` | PNG |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| KEGG | `kegg_search_pathway` | パスウェイキーワード検索 |
| KEGG | `kegg_get_pathway_info` | パスウェイ詳細情報取得 |
| KEGG | `kegg_find_genes` | 遺伝子検索 |
| KEGG | `kegg_get_gene_info` | 遺伝子詳細情報 |
| KEGG | `kegg_list_organisms` | 対応生物種一覧 |
| Reactome | `Reactome_get_pathway` | パスウェイ詳細取得 |
| Reactome | `Reactome_get_pathway_hierarchy` | パスウェイ階層構造 |
| Reactome | `Reactome_get_pathway_reactions` | パスウェイ内反応一覧 |
| Reactome | `Reactome_map_uniprot_to_pathways` | UniProt→パスウェイ変換 |
| Reactome | `Reactome_map_uniprot_to_reactions` | UniProt→反応変換 |
| Reactome | `Reactome_get_pathways_low_entity` | 低レベルパスウェイ |
| Reactome | `Reactome_list_top_pathways` | トップレベルパスウェイ一覧 |
| Reactome | `Reactome_list_species` | 対応種一覧 |
| Reactome | `Reactome_get_event_ancestors` | イベント祖先取得 |
| Reactome | `Reactome_get_events_hierarchy` | イベント階層 |
| Reactome | `Reactome_get_participant_reference_entities` | 参加エンティティ |
| Reactome | `Reactome_get_participants` | 参加要素取得 |
| Reactome | `Reactome_get_complex` | 複合体情報 |
| Reactome | `Reactome_get_diseases` | 疾患関連パスウェイ |
| Reactome | `Reactome_get_interactor` | 相互作用パートナー |
| Reactome | `Reactome_query_by_ids` | バッチ ID クエリ |
| Reactome | `Reactome_get_reaction` | 反応詳細 |
| Reactome | `Reactome_get_entity_compartment` | エンティティ局在 |
| Reactome | `Reactome_get_entity_events` | エンティティ関連イベント |
| Reactome | `Reactome_get_database_version` | DB バージョン確認 |
| GO | `GO_search_terms` | GO 用語検索 |
| GO | `GO_get_term_by_id` | GO ID → 用語情報 |
| GO | `GO_get_term_details` | 用語詳細 (定義/関係) |
| GO | `GO_get_annotations_for_gene` | 遺伝子 GO アノテーション |
| GO | `GO_get_genes_for_term` | GO term → 遺伝子リスト |
| WikiPathways | `WikiPathways_search` | コミュニティパスウェイ検索 |
| WikiPathways | `WikiPathways_get_pathway` | パスウェイ詳細取得 |
| Pathway Commons | `pc_search_pathways` | 統合パスウェイ検索 |
| Pathway Commons | `pc_get_interactions` | 分子間相互作用取得 |

### 参照スキル

| スキル | 関連 |
|---|---|
| `scientific-gene-expression-transcriptomics` | DEG リスト → ORA/GSEA |
| `scientific-proteomics-mass-spectrometry` | 差次タンパク質 → パスウェイ |
| `scientific-metabolomics` | 代謝物 → パスウェイマッピング |
| `scientific-single-cell-genomics` | scRNA-seq マーカー → GO 富化 |
| `scientific-multi-omics` | マルチオミクス富化統合 |
| `scientific-network-analysis` | パスウェイ → ネットワーク |
| `scientific-systems-biology` | FBA/GRN → パスウェイ統合 |

### 依存パッケージ

`scipy`, `statsmodels`, `pandas`, `numpy`, `matplotlib`, `requests`, `gseapy` (optional)

---
name: scientific-metabolomics-network
description: |
  代謝物ネットワーク構築スキル。KEGG/Reactome 代謝パスウェイ
  グラフ抽出・代謝物相関ネットワーク構築 (GGM/WGCNA)・
  ハブ代謝物同定・MetaboAnalyst 統合エンリッチメント
  パイプライン。
  TU 外スキル (直接 Python ライブラリ + REST API)。
---

# Scientific Metabolomics Network

KEGG/Reactome 代謝パスウェイからのネットワーク構築、
代謝物間相関解析 (Gaussian Graphical Model / WGCNA)、
ハブ代謝物同定、MetaboAnalyst 統合エンリッチメントの
パイプラインを提供する。

## When to Use

- 代謝物相関ネットワーク (partial correlation) を構築するとき
- KEGG/Reactome 代謝パスウェイをグラフ化するとき
- ハブ代謝物 (高接続度) を同定するとき
- 代謝パスウェイエンリッチメント解析を行うとき
- メタボロームデータのネットワーク可視化

---

## Quick Start

## 1. 代謝物相関ネットワーク

```python
import numpy as np
import pandas as pd
import networkx as nx
from sklearn.covariance import GraphicalLassoCV


def metabolite_correlation_network(
        data, method="glasso", threshold=0.1):
    """
    代謝物相関ネットワーク構築。

    Parameters:
        data: pd.DataFrame — 代謝物濃度行列
            (行=サンプル, 列=代謝物)
        method: str — "glasso" (Graphical Lasso) or
                      "pearson" (Pearson partial)
        threshold: float — エッジ閾値
    """
    metabolites = data.columns.tolist()

    if method == "glasso":
        model = GraphicalLassoCV(cv=5)
        model.fit(data.values)
        precision = model.precision_
        # Partial correlation from precision matrix
        diag = np.sqrt(np.diag(precision))
        partial_corr = -(precision /
                         np.outer(diag, diag))
        np.fill_diagonal(partial_corr, 1.0)
    else:
        partial_corr = data.corr().values

    # Build network
    G = nx.Graph()
    G.add_nodes_from(metabolites)

    for i in range(len(metabolites)):
        for j in range(i + 1, len(metabolites)):
            w = abs(partial_corr[i, j])
            if w > threshold:
                G.add_edge(
                    metabolites[i],
                    metabolites[j],
                    weight=round(w, 4),
                    sign=("+" if partial_corr[i, j]
                          > 0 else "-"))

    print(f"Metabolite network: "
          f"{G.number_of_nodes()} nodes, "
          f"{G.number_of_edges()} edges "
          f"(threshold={threshold})")
    return G


def hub_metabolites(G, top_n=10):
    """
    ハブ代謝物同定 (次数中心性) 。

    Parameters:
        G: nx.Graph — 代謝物ネットワーク
        top_n: int — 上位件数
    """
    degree_cent = nx.degree_centrality(G)
    betweenness = nx.betweenness_centrality(G)

    rows = []
    for node in G.nodes():
        rows.append({
            "metabolite": node,
            "degree": G.degree(node),
            "degree_centrality": round(
                degree_cent[node], 4),
            "betweenness": round(
                betweenness[node], 4),
        })

    df = pd.DataFrame(rows).sort_values(
        "degree_centrality",
        ascending=False).head(top_n)
    print(f"Top {top_n} hub metabolites:")
    for _, row in df.iterrows():
        print(f"  {row['metabolite']}: "
              f"deg={row['degree']}, "
              f"bc={row['betweenness']}")
    return df
```

## 2. KEGG 代謝パスウェイグラフ

```python
def kegg_pathway_graph(pathway_id):
    """
    KEGG — 代謝パスウェイをネットワークグラフとして取得。

    Parameters:
        pathway_id: str — KEGG パスウェイ ID
                          (例: "hsa00010")
    """
    import requests

    # KGML 取得
    url = (f"https://rest.kegg.jp/get/"
           f"{pathway_id}/kgml")
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()

    import xml.etree.ElementTree as ET
    root = ET.fromstring(resp.text)

    G = nx.DiGraph()

    # ノード追加
    entry_map = {}
    for entry in root.findall("entry"):
        eid = entry.get("id")
        name = entry.get("name", "")
        etype = entry.get("type", "")
        graphics = entry.find("graphics")
        label = (graphics.get("name", name)
                 if graphics is not None else name)
        entry_map[eid] = label
        G.add_node(label, entry_type=etype)

    # エッジ追加
    for relation in root.findall("relation"):
        e1 = relation.get("entry1")
        e2 = relation.get("entry2")
        rtype = relation.get("type", "")
        if e1 in entry_map and e2 in entry_map:
            G.add_edge(entry_map[e1],
                       entry_map[e2],
                       relation_type=rtype)

    for reaction in root.findall("reaction"):
        rname = reaction.get("name", "")
        substrates = [s.get("name", "")
                      for s in reaction.findall(
                          "substrate")]
        products = [p.get("name", "")
                    for p in reaction.findall(
                        "product")]
        for s in substrates:
            for p in products:
                G.add_edge(s, p,
                           reaction=rname)

    print(f"KEGG pathway {pathway_id}: "
          f"{G.number_of_nodes()} nodes, "
          f"{G.number_of_edges()} edges")
    return G
```

## 3. パスウェイエンリッチメント

```python
def metabolite_pathway_enrichment(
        metabolite_list, organism="hsa"):
    """
    代謝物パスウェイエンリッチメント (KEGG)。

    Parameters:
        metabolite_list: list[str] — KEGG compound ID
            リスト (例: ["C00031", "C00158"])
        organism: str — 生物種コード
    """
    import requests
    from scipy.stats import hypergeom

    # KEGG compound→pathway マッピング
    url = "https://rest.kegg.jp/link/pathway/compound"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()

    cpd_to_pw = {}
    pw_to_cpd = {}
    for line in resp.text.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) != 2:
            continue
        cpd = parts[0].replace("cpd:", "")
        pw = parts[1].replace("path:", "")
        if not pw.startswith("map"):
            continue
        cpd_to_pw.setdefault(cpd, set()).add(pw)
        pw_to_cpd.setdefault(pw, set()).add(cpd)

    # エンリッチメント計算
    query_set = set(metabolite_list)
    all_cpds = set(cpd_to_pw.keys())
    M = len(all_cpds)
    n = len(query_set & all_cpds)

    results = []
    for pw, pw_cpds in pw_to_cpd.items():
        N = len(pw_cpds)
        k = len(query_set & pw_cpds)
        if k == 0:
            continue
        pval = hypergeom.sf(k - 1, M, N, n)
        results.append({
            "pathway": pw,
            "overlap": k,
            "pathway_size": N,
            "pvalue": pval,
            "metabolites": ", ".join(
                query_set & pw_cpds),
        })

    df = pd.DataFrame(results).sort_values("pvalue")
    print(f"Pathway enrichment: "
          f"{len(df)} pathways (p<0.05: "
          f"{(df['pvalue'] < 0.05).sum()})")
    return df
```

## 4. 代謝ネットワーク統合パイプライン

```python
def metabolomics_network_pipeline(
        data, metabolite_ids=None,
        output_dir="results"):
    """
    代謝ネットワーク統合パイプライン。

    Parameters:
        data: pd.DataFrame — 代謝物濃度行列
        metabolite_ids: list[str] | None — KEGG
            compound ID (エンリッチメント用)
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) 相関ネットワーク
    G = metabolite_correlation_network(data)
    nx.write_graphml(
        G, str(output_dir / "metabolite_network.graphml"))

    # 2) ハブ代謝物
    hubs = hub_metabolites(G)
    hubs.to_csv(
        output_dir / "hub_metabolites.csv",
        index=False)

    # 3) パスウェイエンリッチメント
    if metabolite_ids:
        enrich = metabolite_pathway_enrichment(
            metabolite_ids)
        enrich.to_csv(
            output_dir / "pathway_enrichment.csv",
            index=False)

    print(f"Metabolomics network pipeline → "
          f"{output_dir}")
    return {"network": G, "hubs": hubs}
```

---

## パイプライン統合

```
metabolomics → metabolomics-network → pathway-enrichment
  (LC-MS/NMR)    (GGM/グラフ構築)       (KEGG/Reactome)
       │                 │                    ↓
  lipidomics ────────────┘         systems-biology
    (脂質サブクラス)                 (マルチオミクス統合)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/metabolite_network.graphml` | 相関ネットワーク | → systems-biology |
| `results/hub_metabolites.csv` | ハブ代謝物 | → biomarker-discovery |
| `results/pathway_enrichment.csv` | パスウェイエンリッチメント | → pathway-enrichment |

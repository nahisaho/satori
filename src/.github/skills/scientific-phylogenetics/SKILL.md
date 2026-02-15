---
name: scientific-phylogenetics
description: |
  系統解析スキル。ete3/ETE Toolkit による系統樹構築・可視化、
  scikit-bio 系統的多様性、配列アライメントベース進化解析、
  分子時計・分岐年代推定、祖先配列再構成パイプライン。
tu_tools:
  - key: ncbi_taxonomy
    name: NCBI Taxonomy
    description: 系統分類・分岐学データ検索
---

# Scientific Phylogenetics

ETE Toolkit / scikit-bio を中心とした
分子系統解析・進化生物学パイプラインを提供する。

## When to Use

- 分子系統樹を構築・可視化するとき (NJ/ML/ベイズ法)
- 多重配列アライメントから系統推定するとき
- 分岐年代推定 (分子時計) を行うとき
- 系統的多様性 (PD: Phylogenetic Diversity) を計算するとき
- 祖先配列再構成を行うとき
- 系統比較法 (PGLS 等) で形質進化を解析するとき

---

## Quick Start

## 1. ETE Toolkit 系統樹構築

```python
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace
import subprocess
import tempfile


def build_phylogenetic_tree(sequences_fasta, method="fasttree", model="GTR"):
    """
    配列アライメントから系統樹構築。

    Parameters:
        sequences_fasta: str — FASTA ファイルパス (アライン済み)
        method: str — "fasttree", "raxml", "iqtree"
        model: str — 進化モデル ("GTR", "JTT", "WAG", "LG")

    K-Dense: etetoolkit — Phylogenetics toolkit
    """
    commands = {
        "fasttree": ["fasttree", "-gtr", "-nt", sequences_fasta],
        "raxml": [
            "raxmlHPC", "-s", sequences_fasta, "-n", "tree",
            "-m", f"GTRGAMMA", "-p", "12345",
        ],
        "iqtree": [
            "iqtree2", "-s", sequences_fasta,
            "-m", model, "-bb", "1000", "--prefix", "iqtree_out",
        ],
    }

    cmd = commands.get(method, commands["fasttree"])
    result = subprocess.run(cmd, capture_output=True, text=True)

    if method == "fasttree":
        newick = result.stdout
    elif method == "iqtree":
        with open("iqtree_out.treefile", "r") as f:
            newick = f.read()
    else:
        newick = result.stdout

    tree = Tree(newick)
    print(f"Phylogenetic tree ({method}, {model}): "
          f"{len(tree)} leaves, {len(list(tree.traverse()))} total nodes")
    return tree


def visualize_tree(tree, output_file="phylogenetic_tree.png",
                   layout="rectangular", show_support=True):
    """
    ETE3 系統樹可視化。

    Parameters:
        tree: ete3.Tree — 系統樹オブジェクト
        output_file: str — 出力画像パス
        layout: str — "rectangular", "circular"
        show_support: bool — ブートストラップ値を表示
    """
    ts = TreeStyle()
    ts.mode = "c" if layout == "circular" else "r"
    ts.show_leaf_name = True
    ts.show_branch_length = True
    ts.show_branch_support = show_support
    ts.branch_vertical_margin = 10

    # Node styling
    for node in tree.traverse():
        nstyle = NodeStyle()
        if node.is_leaf():
            nstyle["fgcolor"] = "#2196F3"
            nstyle["size"] = 8
        else:
            nstyle["fgcolor"] = "#E91E63"
            nstyle["size"] = 5
            if show_support and node.support >= 0.9:
                nstyle["fgcolor"] = "#4CAF50"
        node.set_style(nstyle)

    tree.render(output_file, tree_style=ts, w=800, units="px")
    print(f"Tree rendered: {output_file} ({layout} layout)")
    return output_file
```

## 2. 多重配列アライメント

```python
from Bio import AlignIO, SeqIO
from Bio.Align.Applications import MafftCommandline, MuscleCommandline


def run_multiple_alignment(input_fasta, method="mafft", output_fasta=None):
    """
    多重配列アライメント。

    Parameters:
        input_fasta: str — 入力 FASTA パス
        method: str — "mafft", "muscle", "clustalw"
        output_fasta: str — 出力パス
    """
    if output_fasta is None:
        output_fasta = input_fasta.replace(".fasta", f"_aligned_{method}.fasta")

    if method == "mafft":
        cmd = f"mafft --auto {input_fasta} > {output_fasta}"
    elif method == "muscle":
        cmd = f"muscle -in {input_fasta} -out {output_fasta}"
    else:
        cmd = f"clustalw2 -INFILE={input_fasta} -OUTFILE={output_fasta}"

    subprocess.run(cmd, shell=True, check=True)

    alignment = AlignIO.read(output_fasta, "fasta")
    print(f"Alignment ({method}): {len(alignment)} sequences, "
          f"{alignment.get_alignment_length()} positions")
    return alignment
```

## 3. 系統的多様性 (Phylogenetic Diversity)

```python
import skbio
from skbio import TreeNode
from skbio.diversity import alpha_diversity, beta_diversity
import numpy as np


def calculate_phylogenetic_diversity(newick_string, sample_otus):
    """
    系統的多様性 (Faith's PD, UniFrac) 計算。

    Parameters:
        newick_string: str — Newick 形式系統樹
        sample_otus: dict — {sample_id: {otu_id: abundance}}

    K-Dense: scikit-bio — PD & UniFrac
    """
    tree = TreeNode.read([newick_string])

    # Prepare OTU table
    all_otus = sorted(set(
        otu for otus in sample_otus.values() for otu in otus
    ))
    sample_names = list(sample_otus.keys())
    otu_table = np.zeros((len(sample_names), len(all_otus)))
    for i, sample in enumerate(sample_names):
        for j, otu in enumerate(all_otus):
            otu_table[i, j] = sample_otus[sample].get(otu, 0)

    # Faith's PD (alpha diversity)
    pd_values = alpha_diversity("faith_pd", otu_table, ids=sample_names, tree=tree,
                                otu_ids=all_otus)
    print(f"Faith's PD: mean={pd_values.mean():.3f}, "
          f"range=[{pd_values.min():.3f}, {pd_values.max():.3f}]")

    # Weighted UniFrac (beta diversity)
    unifrac_dm = beta_diversity("weighted_unifrac", otu_table,
                                ids=sample_names, tree=tree, otu_ids=all_otus)
    print(f"Weighted UniFrac: mean distance = "
          f"{unifrac_dm.condensed_form().mean():.4f}")

    return {"faith_pd": pd_values, "unifrac": unifrac_dm}
```

## 4. 分子時計・分岐年代推定

```python
def estimate_divergence_times(tree, calibrations, rate_model="strict"):
    """
    分子時計による分岐年代推定。

    Parameters:
        tree: ete3.Tree — 系統樹
        calibrations: dict — {(taxon1, taxon2): (min_age, max_age)}
            e.g., {("human", "mouse"): (85, 95)}  # MYA
        rate_model: str — "strict" or "relaxed"
    """
    # Branch length to relative time conversion
    total_length = max(tree.get_distance(leaf) for leaf in tree.get_leaves())

    # Apply calibration
    for (t1, t2), (min_age, max_age) in calibrations.items():
        node1 = tree.search_nodes(name=t1)
        node2 = tree.search_nodes(name=t2)
        if node1 and node2:
            ancestor = tree.get_common_ancestor(node1[0], node2[0])
            dist = tree.get_distance(ancestor)
            calibration_age = (min_age + max_age) / 2
            rate = dist / calibration_age if calibration_age > 0 else 1
            print(f"Calibration {t1}-{t2}: {calibration_age} MYA, rate={rate:.6f}")

    # Estimate ages for all internal nodes
    node_ages = {}
    for node in tree.traverse("postorder"):
        if not node.is_leaf():
            dist = tree.get_distance(node)
            # Simple proportional dating
            estimated_age = (dist / total_length) * max(
                (min_age + max_age) / 2
                for (min_age, max_age) in calibrations.values()
            )
            node_ages[node.name or f"node_{id(node)}"] = estimated_age

    return node_ages
```

## 5. 祖先配列再構成

```python
def ancestral_sequence_reconstruction(alignment_file, tree_file, model="JTT"):
    """
    最尤法による祖先配列再構成。

    Parameters:
        alignment_file: str — アライメントファイルパス
        tree_file: str — 系統樹ファイルパス (Newick)
        model: str — アミノ酸置換モデル
    """
    # Using IQ-TREE for ASR
    cmd = [
        "iqtree2", "-s", alignment_file, "-te", tree_file,
        "-m", model, "-asr", "--prefix", "asr_output",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode == 0:
        # Parse ancestral sequences
        asr_file = "asr_output.state"
        ancestral_seqs = {}
        if os.path.exists(asr_file):
            import csv
            with open(asr_file) as f:
                reader = csv.reader(f, delimiter="\t")
                for row in reader:
                    if row and not row[0].startswith("#"):
                        node = row[0]
                        site = row[1]
                        state = row[2]
                        if node not in ancestral_seqs:
                            ancestral_seqs[node] = []
                        ancestral_seqs[node].append(state)

        print(f"ASR ({model}): {len(ancestral_seqs)} ancestral nodes reconstructed")
        return ancestral_seqs
    else:
        print(f"ASR failed: {result.stderr[:200]}")
        return None
```

---

## パイプライン出力

| 出力ファイル | 説明 | 連携先スキル |
|---|---|---|
| `results/phylogenetic_tree.nwk` | Newick 系統樹 | → infectious-disease, microbiome |
| `figures/phylogenetic_tree.png` | 系統樹可視化 | → publication-figures, presentation |
| `results/divergence_times.json` | 分岐年代推定 | → population-genetics, environmental-ecology |
| `results/ancestral_sequences.fasta` | 祖先配列 | → protein-design, sequence-analysis |
| `results/phylo_diversity.json` | 系統的多様性 | → microbiome-metagenomics |

## パイプライン統合

```
sequence-analysis ──→ phylogenetics ──→ infectious-disease
  (アライメント)      (系統樹構築)     (病原体系統解析)
                           │
                           ├──→ microbiome-metagenomics (UniFrac)
                           ├──→ population-genetics (分岐推定)
                           └──→ environmental-ecology (系統的多様性)
```

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `ncbi_taxonomy` | NCBI Taxonomy | 系統分類・分岐学データ検索 |

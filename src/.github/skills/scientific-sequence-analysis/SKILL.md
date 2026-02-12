---
name: scientific-sequence-analysis
description: |
  ゲノム配列解析スキル。コドン使用頻度（RSCU/CAI）、ペアワイズアラインメント
  （Needleman-Wunsch/Smith-Waterman）、系統解析（Jukes-Cantor/UPGMA/ブートストラップ）、
  ORF 探索、CpG 島検出、制限酵素マッピング、タンパク質特性（pI/GRAVY/疎水性プロファイル）
  の解析テンプレート。Scientific Skills Exp-09 で確立したパターン。
---

# Scientific Sequence Analysis

DNA / RNA / タンパク質の配列解析パイプライン。bioinformatics スキルが scRNA-seq や
バルク RNA-seq のオミクスレベル解析を扱うのに対し、本スキルは個々の配列レベル
（分子生物学レベル）の解析に特化する。

## When to Use

- DNA / タンパク質配列の特性解析が必要なとき
- コドン使用頻度・バイアスを分析するとき
- ペアワイズ/マルチプルアラインメント、系統樹を作成するとき
- ORF 探索、CpG 島検出、制限酵素マッピングが必要なとき

---

## Quick Start

## 1. 塩基組成解析

```python
from collections import Counter
import numpy as np
import pandas as pd

def sequence_composition(sequence, seq_type="dna"):
    """
    配列の塩基/アミノ酸組成を算出する。

    Returns:
        dict with counts, frequencies, GC content (DNA), MW estimate
    """
    seq = sequence.upper()
    counts = Counter(seq)
    total = len(seq)
    freq = {k: v / total * 100 for k, v in counts.items()}

    result = {"length": total, "composition": counts, "frequency_pct": freq}

    if seq_type == "dna":
        gc = (counts.get("G", 0) + counts.get("C", 0)) / total * 100
        at = (counts.get("A", 0) + counts.get("T", 0)) / total * 100
        result["gc_content_pct"] = gc
        result["at_content_pct"] = at
    elif seq_type == "rna":
        gc = (counts.get("G", 0) + counts.get("C", 0)) / total * 100
        result["gc_content_pct"] = gc

    return result
```

## 2. コドン使用頻度（RSCU / CAI）

```python
from Bio.Data import CodonTable

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

def compute_rscu(cds_sequence):
    """
    Relative Synonymous Codon Usage (RSCU) を算出する。

    RSCU = observed / expected
    expected = total for AA / number of synonymous codons
    RSCU = 1: 均等使用, >1: 優先使用, <1: 劣位使用
    """
    seq = cds_sequence.upper().replace("U", "T")
    codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3)]
    codon_counts = Counter(codons)

    # アミノ酸ごとの同義コドン集計
    aa_groups = {}
    for codon, aa in CODON_TABLE.items():
        if aa == "*":
            continue
        aa_groups.setdefault(aa, []).append(codon)

    rscu = {}
    for aa, synonymous in aa_groups.items():
        total = sum(codon_counts.get(c, 0) for c in synonymous)
        n_syn = len(synonymous)
        for codon in synonymous:
            observed = codon_counts.get(codon, 0)
            expected = total / n_syn if n_syn > 0 else 0
            rscu[codon] = observed / expected if expected > 0 else 0

    return rscu


def codon_adaptation_index(cds_sequence, reference_rscu):
    """
    Codon Adaptation Index (CAI) を算出する。

    CAI = exp(1/L * Σ ln(w_i))
    w_i = RSCU_i / RSCU_max (for synonymous family)
    """
    seq = cds_sequence.upper().replace("U", "T")
    codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3)]

    # 各同義ファミリーの最大 RSCU
    aa_groups = {}
    for codon, aa in CODON_TABLE.items():
        if aa == "*":
            continue
        aa_groups.setdefault(aa, []).append(codon)

    max_rscu = {}
    for aa, synonymous in aa_groups.items():
        max_val = max(reference_rscu.get(c, 0) for c in synonymous)
        for c in synonymous:
            max_rscu[c] = max_val

    w_values = []
    for codon in codons:
        if codon in CODON_TABLE and CODON_TABLE[codon] != "*":
            w = reference_rscu.get(codon, 0) / max_rscu.get(codon, 1)
            if w > 0:
                w_values.append(np.log(w))

    cai = np.exp(np.mean(w_values)) if w_values else 0
    return cai
```

## 3. ペアワイズアラインメント

```python
def needleman_wunsch(seq1, seq2, match=2, mismatch=-1, gap=-2):
    """
    Needleman-Wunsch グローバルアラインメント（動的計画法）。

    Returns:
        aligned_seq1, aligned_seq2, score
    """
    n, m = len(seq1), len(seq2)
    dp = np.zeros((n+1, m+1))
    traceback = np.zeros((n+1, m+1), dtype=int)  # 0=diag, 1=up, 2=left

    for i in range(1, n+1):
        dp[i][0] = i * gap
    for j in range(1, m+1):
        dp[0][j] = j * gap

    for i in range(1, n+1):
        for j in range(1, m+1):
            s = match if seq1[i-1] == seq2[j-1] else mismatch
            scores = [dp[i-1][j-1] + s, dp[i-1][j] + gap, dp[i][j-1] + gap]
            dp[i][j] = max(scores)
            traceback[i][j] = np.argmax(scores)

    # トレースバック
    a1, a2 = [], []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and traceback[i][j] == 0:
            a1.append(seq1[i-1]); a2.append(seq2[j-1]); i -= 1; j -= 1
        elif i > 0 and traceback[i][j] == 1:
            a1.append(seq1[i-1]); a2.append("-"); i -= 1
        else:
            a1.append("-"); a2.append(seq2[j-1]); j -= 1

    return "".join(reversed(a1)), "".join(reversed(a2)), dp[n][m]
```

## 4. 系統解析

```python
def jukes_cantor_distance(seq1, seq2):
    """Jukes-Cantor 距離: d = -3/4 ln(1 - 4p/3)"""
    aligned_len = min(len(seq1), len(seq2))
    mismatches = sum(1 for a, b in zip(seq1[:aligned_len], seq2[:aligned_len])
                     if a != b and a != "-" and b != "-")
    valid = sum(1 for a, b in zip(seq1[:aligned_len], seq2[:aligned_len])
                if a != "-" and b != "-")
    p = mismatches / valid if valid > 0 else 0

    if p >= 0.75:
        return float("inf")
    return -0.75 * np.log(1 - 4 * p / 3)


def upgma_tree(distance_matrix, names):
    """
    UPGMA (Unweighted Pair Group Method with Arithmetic Mean) 系統樹を構築する。

    Returns:
        Newick 形式の文字列
    """
    n = len(names)
    dm = distance_matrix.copy()
    clusters = {i: names[i] for i in range(n)}
    sizes = {i: 1 for i in range(n)}

    while len(clusters) > 1:
        keys = list(clusters.keys())
        min_dist = float("inf")
        merge_i, merge_j = 0, 0

        for a in range(len(keys)):
            for b in range(a+1, len(keys)):
                if dm[keys[a]][keys[b]] < min_dist:
                    min_dist = dm[keys[a]][keys[b]]
                    merge_i, merge_j = keys[a], keys[b]

        height = min_dist / 2
        new_name = f"({clusters[merge_i]}:{height:.4f},{clusters[merge_j]}:{height:.4f})"

        new_id = max(clusters.keys()) + 1
        clusters[new_id] = new_name
        sizes[new_id] = sizes[merge_i] + sizes[merge_j]

        # 新ノードへの距離更新（加重平均）
        new_row = {}
        for k in clusters.keys():
            if k != new_id:
                d = (sizes[merge_i] * dm[merge_i].get(k, 0) +
                     sizes[merge_j] * dm[merge_j].get(k, 0)) / sizes[new_id]
                new_row[k] = d

        dm[new_id] = new_row
        for k in new_row:
            dm.setdefault(k, {})[new_id] = new_row[k]

        del clusters[merge_i]
        del clusters[merge_j]

    return list(clusters.values())[0] + ";"
```

## 5. ORF 探索

```python
def find_orfs(sequence, min_length_aa=100):
    """
    全6リーディングフレームからORF（Open Reading Frame）を探索する。

    Returns:
        list of dict with frame, start, end, length_aa, protein_seq
    """
    seq = sequence.upper()
    reverse_comp = seq.translate(str.maketrans("ATGC", "TACG"))[::-1]

    orfs = []
    for strand, s in [("+", seq), ("-", reverse_comp)]:
        for frame in range(3):
            i = frame
            while i < len(s) - 2:
                codon = s[i:i+3]
                if codon == "ATG":
                    # ORF 開始 → 停止コドンまでスキャン
                    protein = []
                    j = i
                    while j < len(s) - 2:
                        c = s[j:j+3]
                        aa = CODON_TABLE.get(c, "X")
                        if aa == "*":
                            break
                        protein.append(aa)
                        j += 3
                    if len(protein) >= min_length_aa:
                        orfs.append({
                            "strand": strand,
                            "frame": frame + 1,
                            "start": i + 1,
                            "end": j + 3,
                            "length_aa": len(protein),
                            "protein_seq": "".join(protein[:50]) + "...",
                        })
                    i = j + 3
                else:
                    i += 3

    return sorted(orfs, key=lambda x: x["length_aa"], reverse=True)
```

## 6. CpG 島検出

```python
def detect_cpg_islands(sequence, window=200, step=1,
                        min_gc=0.50, min_obs_exp=0.60, min_length=200):
    """
    配列中の CpG island を検出する。

    判定基準（Gardiner-Garden & Frommer, 1987）:
    - GC含量 ≥ 50%
    - CpG observed/expected ≥ 0.60
    - 長さ ≥ 200 bp
    """
    seq = sequence.upper()
    islands = []
    in_island = False
    start = 0

    for i in range(0, len(seq) - window, step):
        w = seq[i:i+window]
        gc = (w.count("G") + w.count("C")) / window
        cpg_obs = w.count("CG") / window
        c_freq = w.count("C") / window
        g_freq = w.count("G") / window
        cpg_exp = c_freq * g_freq
        obs_exp = cpg_obs / cpg_exp if cpg_exp > 0 else 0

        if gc >= min_gc and obs_exp >= min_obs_exp:
            if not in_island:
                start = i
                in_island = True
        else:
            if in_island:
                length = i - start + window
                if length >= min_length:
                    islands.append({"start": start+1, "end": i+window,
                                   "length": length})
                in_island = False

    return islands
```

## 7. タンパク質特性

```python
AMINO_ACID_MW = {
    "A": 89.09, "R": 174.20, "N": 132.12, "D": 133.10, "C": 121.16,
    "E": 147.13, "Q": 146.15, "G": 75.03, "H": 155.16, "I": 131.17,
    "L": 131.17, "K": 146.19, "M": 149.21, "F": 165.19, "P": 115.13,
    "S": 105.09, "T": 119.12, "W": 204.23, "Y": 181.19, "V": 117.15,
}

KYTE_DOOLITTLE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "E": -3.5, "Q": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

def protein_properties(protein_seq):
    """タンパク質の物理化学的特性を算出する。"""
    seq = protein_seq.upper()
    mw = sum(AMINO_ACID_MW.get(aa, 0) for aa in seq) - 18.015 * (len(seq) - 1)
    gravy = np.mean([KYTE_DOOLITTLE.get(aa, 0) for aa in seq])

    return {
        "length_aa": len(seq),
        "molecular_weight_Da": mw,
        "gravy": gravy,
        "hydrophobic_pct": sum(1 for aa in seq if KYTE_DOOLITTLE.get(aa, 0) > 0) / len(seq) * 100,
    }
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/sequence_composition.csv` | CSV |
| `results/rscu_analysis.csv` | CSV |
| `results/orf_predictions.csv` | CSV |
| `results/cpg_islands.csv` | CSV |
| `figures/codon_usage_heatmap.png` | PNG |
| `figures/phylogenetic_tree.png` | PNG |
| `figures/hydrophobicity_profile.png` | PNG |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| BLAST | `BLAST_protein_search` | タンパク質相同性検索 |
| BLAST | `BLAST_nucleotide_search` | 核酸相同性検索 |
| UniProt | `UniProt_get_sequence_by_accession` | アミノ酸配列取得 |
| NCBI | `NCBI_get_sequence` | ヌクレオチド配列取得 |
| InterPro | `InterProScan_scan_sequence` | 配列ドメインスキャン |
| InterPro | `InterPro_get_protein_domains` | ドメインアノテーション |

#### 参照実験

- **Exp-09**: コドン使用頻度、ペアワイズアラインメント、系統解析、ORF 探索、CpG 島検出

---
name: scientific-crispr-design
description: |
  CRISPR gRNA 設計スキル。Cas9/Cas12a PAM 配列検索・
  オフターゲットスコアリング (CFD/MIT)・
  CRISPRscan/Rule Set 2 活性予測・検証プライマー設計・
  sgRNA スクリーニングライブラリ構築パイプライン。
  TU 外スキル (Python ライブラリ + ローカル解析)。
---

# Scientific CRISPR Design

CRISPR gRNA 設計・オフターゲット評価・活性予測を統合した
効率的なガイド RNA 選択パイプラインを提供する。

## When to Use

- CRISPR-Cas9/Cas12a の gRNA を設計するとき
- PAM 配列検索とガイド候補の列挙を行うとき
- オフターゲットスコア (CFD/MIT) で安全性を評価するとき
- gRNA 活性スコア (CRISPRscan/Rule Set 2) で効率を予測するとき
- CRISPR スクリーニングライブラリを構築するとき
- 検証用 PCR プライマーを設計するとき

---

## Quick Start

## 1. PAM 配列検索・gRNA 候補列挙

```python
import re
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


# PAM パターン定義
PAM_PATTERNS = {
    "SpCas9":  {"pam": "NGG", "guide_len": 20,
                "pam_side": "3prime"},
    "SaCas9":  {"pam": "NNGRRT", "guide_len": 21,
                "pam_side": "3prime"},
    "Cas12a":  {"pam": "TTTV", "guide_len": 23,
                "pam_side": "5prime"},
    "xCas9":   {"pam": "NG", "guide_len": 20,
                "pam_side": "3prime"},
}


def iupac_to_regex(pam):
    """IUPAC → 正規表現変換。"""
    iupac = {
        "N": "[ACGT]", "R": "[AG]", "Y": "[CT]",
        "S": "[GC]", "W": "[AT]", "K": "[GT]",
        "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
        "H": "[ACT]", "V": "[ACG]",
    }
    return "".join(iupac.get(c, c) for c in pam)


def find_grna_candidates(sequence, cas_type="SpCas9",
                           strand="both"):
    """
    gRNA 候補の列挙。

    Parameters:
        sequence: str — 標的 DNA 配列
        cas_type: str — Cas タイプ
        strand: str — "sense"/"antisense"/"both"
    """
    config = PAM_PATTERNS[cas_type]
    pam_re = iupac_to_regex(config["pam"])
    gl = config["guide_len"]
    side = config["pam_side"]
    seq = sequence.upper()

    candidates = []

    def _search_strand(s, s_name):
        for m in re.finditer(
                f"(?=({pam_re}))", s):
            pos = m.start()
            if side == "3prime":
                start = pos - gl
                if start < 0:
                    continue
                guide = s[start:pos]
            else:  # 5prime
                start = pos + len(config["pam"])
                end = start + gl
                if end > len(s):
                    continue
                guide = s[start:end]

            if len(guide) != gl:
                continue

            gc = (guide.count("G")
                  + guide.count("C")) / gl

            candidates.append({
                "guide": guide,
                "pam": m.group(1),
                "position": pos,
                "strand": s_name,
                "gc_content": round(gc, 3),
                "length": gl,
            })

    if strand in ("sense", "both"):
        _search_strand(seq, "+")
    if strand in ("antisense", "both"):
        rc = str(Seq(seq).reverse_complement())
        _search_strand(rc, "-")

    df = pd.DataFrame(candidates)

    # GC フィルタ (30-70%)
    if not df.empty:
        df = df[(df["gc_content"] >= 0.30)
                & (df["gc_content"] <= 0.70)]

    print(f"CRISPR {cas_type}: "
          f"{len(df)} gRNA candidates "
          f"(GC 30-70%)")
    return df.reset_index(drop=True)
```

## 2. オフターゲットスコアリング

```python
import numpy as np


# CFD スコア簡易実装 (Doench 2016)
def cfd_score(guide, off_target):
    """
    CFD (Cutting Frequency Determination) スコア。

    Parameters:
        guide: str — gRNA 配列 (20nt)
        off_target: str — オフターゲットサイト
    """
    # ポジション別ミスマッチペナルティ (簡易版)
    mm_penalty = {
        1: 0.0, 2: 0.0, 3: 0.014, 4: 0.0,
        5: 0.0, 6: 0.395, 7: 0.317, 8: 0.0,
        9: 0.389, 10: 0.079, 11: 0.445,
        12: 0.508, 13: 0.613, 14: 0.851,
        15: 0.732, 16: 0.828, 17: 0.615,
        18: 0.804, 19: 0.685, 20: 0.583,
    }

    score = 1.0
    for i in range(min(len(guide),
                       len(off_target))):
        if guide[i] != off_target[i]:
            pos = i + 1
            penalty = mm_penalty.get(pos, 0.5)
            score *= (1.0 - penalty)

    return round(score, 4)


def score_off_targets(guide, genome_fasta,
                        max_mismatches=4):
    """
    ゲノムワイドオフターゲットスコアリング。

    Parameters:
        guide: str — gRNA 配列
        genome_fasta: str — リファレンスゲノム
        max_mismatches: int — 最大ミスマッチ数
    """
    results = []
    gl = len(guide)
    guide_upper = guide.upper()

    for record in SeqIO.parse(
            genome_fasta, "fasta"):
        seq = str(record.seq).upper()
        for i in range(len(seq) - gl - 3):
            site = seq[i:i + gl]
            pam = seq[i + gl:i + gl + 3]
            if not re.match("[ACGT]GG", pam):
                continue

            mm = sum(1 for a, b in
                     zip(guide_upper, site)
                     if a != b)
            if mm <= max_mismatches:
                results.append({
                    "chrom": record.id,
                    "position": i,
                    "site": site,
                    "pam": pam,
                    "mismatches": mm,
                    "cfd_score": cfd_score(
                        guide_upper, site),
                })

    df = pd.DataFrame(results)
    df = df.sort_values("cfd_score",
                         ascending=False)
    print(f"Off-target: {len(df)} sites "
          f"(≤{max_mismatches} mm)")
    return df
```

## 3. gRNA 活性予測

```python
def rule_set2_score(guide_30mer):
    """
    Rule Set 2 活性スコア (Doench 2016 簡易版)。

    Parameters:
        guide_30mer: str — 30nt 配列
            (4nt upstream + 20nt guide + 3nt PAM
             + 3nt downstream)
    """
    seq = guide_30mer.upper()
    if len(seq) != 30:
        print(f"Warning: expected 30nt, "
              f"got {len(seq)}")
        return 0.0

    guide = seq[4:24]
    gc = (guide.count("G")
          + guide.count("C")) / 20

    # 位置重み付きスコア (簡易)
    score = 0.5

    # GC 最適範囲
    if 0.40 <= gc <= 0.70:
        score += 0.1
    elif gc < 0.30 or gc > 0.80:
        score -= 0.2

    # PAM 近傍優先塩基
    if guide[-1] == "G":
        score += 0.05
    if guide[-4] == "C":
        score += 0.03

    # ポリ T 回避 (Pol III 終結)
    if "TTTT" in guide:
        score -= 0.3

    return round(max(0, min(1, score)), 3)


def rank_grnas(candidates_df, genome_fasta=None):
    """
    gRNA 候補ランキング。

    Parameters:
        candidates_df: pd.DataFrame — gRNA 候補
        genome_fasta: str | None — オフタ解析用
    """
    df = candidates_df.copy()

    # 活性スコア (30mer が無い場合は guide のみ)
    df["activity_score"] = df["guide"].apply(
        lambda g: rule_set2_score(
            "AAAA" + g + "GGGNNN"
            if len(g) == 20
            else g.ljust(30, "N")))

    # オフターゲット (ゲノムがあれば)
    if genome_fasta:
        ot_scores = []
        for guide in df["guide"]:
            ot = score_off_targets(
                guide, genome_fasta, 3)
            specificity = (
                1.0 / (1.0 + len(ot))
                if not ot.empty else 1.0)
            ot_scores.append(round(specificity, 3))
        df["specificity"] = ot_scores
    else:
        df["specificity"] = 1.0

    # 総合スコア
    df["composite_score"] = (
        df["activity_score"] * 0.5
        + df["specificity"] * 0.3
        + df["gc_content"].clip(0.4, 0.6) * 0.2
    ).round(3)

    df = df.sort_values("composite_score",
                         ascending=False)
    print(f"gRNA ranking: top score = "
          f"{df['composite_score'].iloc[0]}")
    return df
```

## 4. sgRNA ライブラリ構築

```python
def build_sgrna_library(gene_list,
                          genome_fasta,
                          guides_per_gene=4,
                          cas_type="SpCas9"):
    """
    スクリーニング用 sgRNA ライブラリ構築。

    Parameters:
        gene_list: list[dict] — 遺伝子リスト
            [{"gene": "TP53", "sequence": "ATCG..."}]
        genome_fasta: str — リファレンスゲノム
        guides_per_gene: int — 遺伝子あたり gRNA 数
        cas_type: str — Cas タイプ
    """
    library = []

    for gene_info in gene_list:
        gene = gene_info["gene"]
        seq = gene_info["sequence"]

        candidates = find_grna_candidates(
            seq, cas_type)

        if candidates.empty:
            print(f"  {gene}: no candidates")
            continue

        ranked = rank_grnas(candidates)
        top = ranked.head(guides_per_gene)

        for _, row in top.iterrows():
            library.append({
                "gene": gene,
                "guide": row["guide"],
                "position": row["position"],
                "strand": row["strand"],
                "gc_content": row["gc_content"],
                "activity": row["activity_score"],
                "composite": row["composite_score"],
            })

    df = pd.DataFrame(library)
    n_genes = df["gene"].nunique()
    print(f"Library: {len(df)} sgRNAs for "
          f"{n_genes} genes")
    return df
```

---

## パイプライン統合

```
genome-sequence-tools → crispr-design → perturbation-analysis
  (ゲノム配列取得)       (gRNA 設計)     (摂動実験解析)
        │                     │                ↓
  variant-effect-prediction ─┘     functional-genomics
    (変異影響予測)                   (機能ゲノミクス)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `grna_candidates.csv` | gRNA 候補リスト | → ランキング |
| `off_target_report.csv` | オフターゲット評価 | → 安全性確認 |
| `sgrna_library.csv` | sgRNA ライブラリ | → perturbation-analysis |

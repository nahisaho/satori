---
name: scientific-drug-repurposing
description: |
  ドラッグリポジショニング（既存薬再創薬）スキル。ToolUniverse の Drug Repurposing
  パラダイムに準拠し、7 つの戦略（ターゲット型、化合物型、疾患駆動型、メカニズム型、
  ネットワーク型、表現型、構造型）で候補を体系的に探索。
  「ドラッグリポジショニングして」「既存薬の新規適応を探して」で発火。
---

# Scientific Drug Repurposing

既存薬・承認薬の新規適応症探索スキル。ToolUniverse（mims-harvard）の
Drug Repurposing 戦略を SATORI に統合し、多角的アプローチで
リポジショニング候補を体系的に評価する。

## When to Use

- 既存薬の新規適応症を探索するとき
- 疾患に対する既存治療候補を網羅的にスクリーニングするとき
- ネットワーク薬理学によるメカニズムベース探索を行うとき
- バーチャルスクリーニングで化合物ライブラリを評価するとき
- 緊急時（パンデミック等）の迅速な候補薬探索を行うとき

## Quick Start

### 7 つのリポジショニング戦略

```
Strategy 1: Target-Based（ターゲット型）
  Disease → Target genes → Known drugs for targets → Validation

Strategy 2: Compound-Based（化合物型）
  Drug → All known targets → Other diseases of targets → Score

Strategy 3: Disease-Driven（疾患駆動型）
  Disease → DEGs/GWAS targets → Pathway → Drugs in pathway

Strategy 4: Mechanism-Based（メカニズム型）
  Known MOA → Drugs with similar MOA → New indications

Strategy 5: Network-Based（ネットワーク型）
  Disease module → Proximity analysis → Proximal drugs

Strategy 6: Phenotype-Based（表現型型）
  Adverse events → Therapeutic potential → Repurpose

Strategy 7: Structure-Based（構造型）
  Active compound → Similar structures → Approved analogs
```

---

## Strategy 1: Target-Based Repurposing

### ワークフロー

```python
def target_based_repurposing(disease_name):
    """
    疾患 → ターゲット → 既存薬 のリポジショニングパイプライン。
    """
    pipeline = {
        "step1_disease": "EFO/MONDO IDで疾患を同定",
        "step2_targets": "Open Targets で上位 50 ターゲットを取得",
        "step3_drugs": "ターゲット毎に ChEMBL/DrugBank で既存薬検索",
        "step4_filter": "現在の適応症と異なる薬剤を候補に",
        "step5_validate": "文献エビデンス + 構造的妥当性で検証",
        "step6_score": "多基準スコアリングで最終ランキング",
    }
    return pipeline

# スコアリング関数
def repurposing_score(candidate):
    """
    リポジショニング候補のスコアリング。
    """
    score = 0.0
    weights = {
        "target_association": 0.25,  # Open Targets score
        "drug_approval_status": 0.20,  # FDA approved > Phase3 > Phase2
        "mechanism_relevance": 0.20,  # MOA と疾患の関連性
        "safety_profile": 0.15,       # 既知安全性プロファイル
        "literature_evidence": 0.10,  # PubMed エビデンス
        "structural_fit": 0.10,       # 分子ドッキングスコア
    }
    for key, weight in weights.items():
        score += candidate.get(key, 0) * weight
    return score
```

---

## Strategy 5: Network-Based Repurposing

### ネットワーク近接度解析

```python
import networkx as nx
import numpy as np

def network_proximity(drug_targets, disease_genes, ppi_network):
    """
    ネットワーク近接度法によるリポジショニング評価。
    Guney et al., Nature Communications 2016 に準拠。

    d(S,T) = 1/|T| * Σ min d(s,t) for s in S
    """
    G = ppi_network  # PPI ネットワーク (NetworkX Graph)
    distances = []

    for dt in drug_targets:
        if dt not in G:
            continue
        min_dist = float('inf')
        for dg in disease_genes:
            if dg not in G:
                continue
            try:
                d = nx.shortest_path_length(G, dt, dg)
                min_dist = min(min_dist, d)
            except nx.NetworkXNoPath:
                continue
        if min_dist < float('inf'):
            distances.append(min_dist)

    if not distances:
        return {"proximity": None, "significance": "N/A"}

    proximity = np.mean(distances)

    return {
        "proximity": proximity,
        "significance": "Close" if proximity < 2.0 else "Moderate" if proximity < 4.0 else "Distant",
        "n_drug_targets_in_network": len([t for t in drug_targets if t in G]),
        "n_disease_genes_in_network": len([g for g in disease_genes if g in G]),
    }
```

---

## Candidate Evaluation Matrix

### 多基準評価テンプレート

```markdown
## Repurposing Candidate Evaluation

### Candidate: [Drug Name]
| Criterion | Score (0-1) | Weight | Evidence |
|-----------|-------------|--------|----------|
| Target-disease association | | 0.25 | |
| FDA approval status | | 0.20 | |
| Mechanism relevance | | 0.20 | |
| Safety profile | | 0.15 | |
| Literature support | | 0.10 | |
| Structural compatibility | | 0.10 | |
| **Weighted Total** | | | |

### Evidence Summary
- Clinical trials: [NCT IDs]
- Case reports: [PMID]
- Mechanistic studies: [PMID]

### Risk Assessment
- Patent status: ___
- Formulation feasibility: ___
- Regulatory pathway: 505(b)(2) / Orphan / Standard
```

---

## Report Template

```markdown
# Drug Repurposing Report: [DISEASE / DRUG]

**Strategy**: [Target-Based / Network-Based / etc.]
**Date**: [date]

## 1. Background & Rationale

## 2. Methodology
### 2.1 Strategy Selection Justification
### 2.2 Data Sources
### 2.3 Scoring Criteria

## 3. Top Candidates
| Rank | Drug | Original Indication | Score | Key Evidence |
|------|------|---------------------|-------|-------------|
| 1 | | | | |
| 2 | | | | |
| 3 | | | | |

## 4. Detailed Candidate Profiles
### 4.1 Candidate 1
### 4.2 Candidate 2
### 4.3 Candidate 3

## 5. Validation Evidence

## 6. Recommendations
### 6.1 Immediate Next Steps
### 6.2 Clinical Trial Design Considerations
### 6.3 Regulatory Path

## 7. Limitations & Caveats

## 8. Data Sources
```

---

## Completeness Checklist

- [ ] 戦略選択: 7 戦略のうち少なくとも 2 つを適用
- [ ] データソース: Open Targets + ChEMBL + DrugBank を最低限使用
- [ ] 候補数: 上位 10 候補をリストアップ
- [ ] スコアリング: 多基準スコアで定量評価
- [ ] エビデンス: 各候補に文献エビデンスを付与
- [ ] 安全性: 既知の副作用プロファイルを確認
- [ ] 特許: 主要候補のパテント状況を確認

## Best Practices

1. **複数戦略を組み合わせ**: 単一戦略の偏りを回避
2. **承認薬を優先**: リスク最小化のため FDA 承認済み薬を上位に
3. **ネガティブエビデンスも記載**: 過去の失敗臨床試験を確認
4. **用量を検討**: リポジショニングでは用量が異なる場合がある
5. **ランダム化比較で検証**: 最終的な有効性は RCT で確認

## References

### Output Files

| ファイル | 形式 | 生成タイミング |
|---|---|---|
| `results/repurposing_candidates.json` | リポジショニング候補リスト（JSON） | 多基準スコアリング完了時 |
| `results/repurposing_report.md` | リポジショニング評価レポート（Markdown） | 全解析完了時 |
| `results/network_proximity.json` | ネットワーク近接スコア（JSON） | ネットワーク解析完了時 |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| ChEMBL | `ChEMBL_search_drugs` | 薬物検索 |
| ChEMBL | `ChEMBL_get_drug_mechanisms` | 薬物作用機序 |
| OpenTargets | `OpenTargets_get_associated_drugs_by_disease_efoId` | 疾患関連薬物 |
| DGIdb | `DGIdb_get_drug_gene_interactions` | 薬物-遺伝子相互作用 |
| ClinicalTrials | `search_clinical_trials` | 臨床試験マッチング |
| FAERS | `FAERS_count_reactions_by_drug_event` | 有害事象データ |
| PubMed | `PubMed_search_articles` | リポジショニング文献 |

### 参照スキル

| スキル | 連携 |
|---|---|
| `scientific-drug-target-profiling` | ← ターゲット情報の利用 |
| `scientific-admet-pharmacokinetics` | ← ADMET フィルタリング済み化合物 |
| `scientific-network-analysis` | ← 疾患-薬物ネットワーク構築 |
| `scientific-deep-research` | ← 文献エビデンス収集 |
| `scientific-clinical-decision-support` | → 候補薬の臨床意思決定 |
| `scientific-academic-writing` | → 研究成果の論文化 |

---
name: scientific-systematic-review
description: |
  PRISMA 2020 準拠系統的レビュースキル。マルチ DB 検索戦略立案
  (PubMed/Embase/Cochrane/Web of Science)、スクリーニングワークフロー
  (タイトル/抄録→全文)、品質評価 (RoB 2/ROBINS-I/NOS)、データ抽出
  テンプレート、PRISMA フロー図自動生成パイプライン。
---

# Scientific Systematic Review

PRISMA 2020 ガイドラインに準拠した
系統的レビュー・メタアナリシスの方法論パイプラインを提供する。

## When to Use

- 系統的レビューの検索戦略を設計するとき
- タイトル/抄録スクリーニングのワークフローが必要なとき
- バイアスリスク (RoB 2, ROBINS-I, NOS) 評価を行うとき
- PRISMA フロー図を生成するとき
- 系統的レビューのデータ抽出テーブルを作成するとき

---

## Quick Start

## 1. 検索戦略設計 (PICO → クエリ)

```python
import pandas as pd
import json


def design_search_strategy(pico, databases=None):
    """
    PICO フレームワークから検索戦略を設計。

    Parameters:
        pico: dict — {"P": "...", "I": "...", "C": "...", "O": "..."}
        databases: list — ["PubMed", "Embase", "Cochrane", "Web of Science"]
    """
    if databases is None:
        databases = ["PubMed", "Embase", "Cochrane"]

    strategy = {
        "pico": pico,
        "databases": databases,
        "search_blocks": [],
    }

    # P (Population) ブロック
    p_terms = pico.get("P", "").split(",")
    p_block = {
        "concept": "Population",
        "terms": [t.strip() for t in p_terms],
        "mesh_terms": [],  # 手動で MeSH を追加
        "boolean": "OR",
    }

    # I (Intervention) ブロック
    i_terms = pico.get("I", "").split(",")
    i_block = {
        "concept": "Intervention",
        "terms": [t.strip() for t in i_terms],
        "mesh_terms": [],
        "boolean": "OR",
    }

    # C (Comparison) ブロック
    c_terms = pico.get("C", "").split(",")
    c_block = {
        "concept": "Comparison",
        "terms": [t.strip() for t in c_terms if t.strip()],
        "boolean": "OR",
    }

    # O (Outcome) ブロック
    o_terms = pico.get("O", "").split(",")
    o_block = {
        "concept": "Outcome",
        "terms": [t.strip() for t in o_terms],
        "boolean": "OR",
    }

    strategy["search_blocks"] = [p_block, i_block]
    if c_block["terms"]:
        strategy["search_blocks"].append(c_block)
    if o_block["terms"]:
        strategy["search_blocks"].append(o_block)

    # PubMed クエリ生成
    pubmed_parts = []
    for block in strategy["search_blocks"]:
        terms = [f'"{t}"' for t in block["terms"]]
        mesh = [f'"{m}"[MeSH]' for m in block.get("mesh_terms", [])]
        all_terms = terms + mesh
        pubmed_parts.append(f"({' OR '.join(all_terms)})")

    strategy["pubmed_query"] = " AND ".join(pubmed_parts)

    print(f"Search strategy: {len(strategy['search_blocks'])} blocks, "
          f"{len(databases)} databases")
    print(f"PubMed query: {strategy['pubmed_query'][:200]}...")
    return strategy
```

## 2. スクリーニングワークフロー

```python
def screening_workflow(records_df, stage="title_abstract",
                        inclusion_criteria=None,
                        exclusion_criteria=None):
    """
    スクリーニングワークフロー管理。

    Parameters:
        records_df: DataFrame — columns: [id, title, abstract, source]
        stage: "title_abstract" or "fulltext"
        inclusion_criteria: list — 適格基準
        exclusion_criteria: list — 除外基準
    """
    if inclusion_criteria is None:
        inclusion_criteria = [
            "Published in English or Japanese",
            "Human subjects",
            "Original research (not review/editorial)",
        ]
    if exclusion_criteria is None:
        exclusion_criteria = [
            "Case reports (n < 5)",
            "Conference abstracts only",
            "Animal studies only",
        ]

    # 重複除去
    initial_count = len(records_df)
    records_df = records_df.drop_duplicates(subset=["title"], keep="first")
    duplicates_removed = initial_count - len(records_df)

    # スクリーニング結果テンプレート
    records_df["decision"] = "pending"
    records_df["excluded_reason"] = ""
    records_df["screener"] = ""

    result = {
        "stage": stage,
        "total_records": initial_count,
        "duplicates_removed": duplicates_removed,
        "unique_records": len(records_df),
        "inclusion_criteria": inclusion_criteria,
        "exclusion_criteria": exclusion_criteria,
    }

    print(f"Screening ({stage}): {initial_count} records → "
          f"{duplicates_removed} duplicates removed → "
          f"{len(records_df)} to screen")
    return records_df, result
```

## 3. バイアスリスク評価

```python
def risk_of_bias_assessment(studies_df, tool="RoB2"):
    """
    バイアスリスク評価。

    Parameters:
        studies_df: DataFrame — columns: [study_id, study_type, ...]
        tool: "RoB2" (RCT), "ROBINS-I" (非ランダム化), "NOS" (観察研究)
    """
    if tool == "RoB2":
        # Cochrane RoB 2 — 5 ドメイン
        domains = [
            "D1: Randomization process",
            "D2: Deviations from interventions",
            "D3: Missing outcome data",
            "D4: Measurement of the outcome",
            "D5: Selection of the reported result",
        ]
        levels = ["Low", "Some concerns", "High"]
    elif tool == "ROBINS-I":
        domains = [
            "D1: Confounding",
            "D2: Selection of participants",
            "D3: Classification of interventions",
            "D4: Deviations from intended interventions",
            "D5: Missing data",
            "D6: Measurement of outcomes",
            "D7: Selection of the reported result",
        ]
        levels = ["Low", "Moderate", "Serious", "Critical", "NI"]
    elif tool == "NOS":
        domains = [
            "Selection (0-4 stars)",
            "Comparability (0-2 stars)",
            "Outcome/Exposure (0-3 stars)",
        ]
        levels = ["0-3 (low quality)", "4-6 (moderate)", "7-9 (high quality)"]
    else:
        raise ValueError(f"Unknown tool: {tool}")

    # 評価テンプレート生成
    assessments = []
    for _, study in studies_df.iterrows():
        assessment = {
            "study_id": study.get("study_id", ""),
            "tool": tool,
        }
        for domain in domains:
            assessment[domain] = "pending"
        assessment["overall"] = "pending"
        assessments.append(assessment)

    df = pd.DataFrame(assessments)
    print(f"RoB assessment ({tool}): {len(df)} studies, "
          f"{len(domains)} domains")
    return df
```

## 4. PRISMA フロー図生成

```python
def generate_prisma_flowchart(counts, output="figures/prisma_flow.svg"):
    """
    PRISMA 2020 フロー図の自動生成。

    Parameters:
        counts: dict — {
            "databases": {"PubMed": 500, "Embase": 300, "Cochrane": 100},
            "other_sources": 20,
            "duplicates_removed": 150,
            "title_abstract_screened": 770,
            "title_abstract_excluded": 650,
            "fulltext_assessed": 120,
            "fulltext_excluded": {"not_relevant": 30, "wrong_design": 20, ...},
            "included_qualitative": 70,
            "included_quantitative": 50,
        }
    """
    import os
    os.makedirs(os.path.dirname(output), exist_ok=True)

    # Mermaid 形式で PRISMA フロー生成
    db_counts = counts.get("databases", {})
    total_db = sum(db_counts.values())
    other = counts.get("other_sources", 0)
    total = total_db + other
    dedup = counts.get("duplicates_removed", 0)
    screened = counts.get("title_abstract_screened", total - dedup)
    ta_excluded = counts.get("title_abstract_excluded", 0)
    ft_assessed = counts.get("fulltext_assessed", screened - ta_excluded)
    ft_excluded = counts.get("fulltext_excluded", {})
    ft_excluded_total = sum(ft_excluded.values()) if isinstance(ft_excluded, dict) else ft_excluded
    qualitative = counts.get("included_qualitative", ft_assessed - ft_excluded_total)
    quantitative = counts.get("included_quantitative", qualitative)

    mermaid = f"""flowchart TD
    A[Database検索<br>n={total_db}] --> C[重複除去後<br>n={total - dedup}]
    B[その他ソース<br>n={other}] --> C
    C --> D[タイトル/抄録スクリーニング<br>n={screened}]
    D --> E[除外<br>n={ta_excluded}]
    D --> F[全文評価<br>n={ft_assessed}]
    F --> G[除外<br>n={ft_excluded_total}]
    F --> H[質的統合<br>n={qualitative}]
    H --> I[量的統合 (メタアナリシス)<br>n={quantitative}]
"""

    # SVG として保存 (Mermaid CLI or fallback to text)
    mermaid_file = output.replace(".svg", ".mmd")
    with open(mermaid_file, "w") as f:
        f.write(mermaid)

    print(f"PRISMA flow: {total} identified → {qualitative} included")
    print(f"  Mermaid source: {mermaid_file}")
    return mermaid_file, counts
```

## 5. データ抽出テンプレート

```python
def create_extraction_template(study_type="RCT",
                                 custom_fields=None):
    """
    系統的レビュー用データ抽出テンプレート。

    Parameters:
        study_type: "RCT", "cohort", "cross-sectional", "case-control"
        custom_fields: list — 追加フィールド
    """
    base_fields = [
        "study_id", "first_author", "year", "country",
        "study_design", "sample_size", "population",
        "setting",
    ]

    if study_type == "RCT":
        type_fields = [
            "intervention", "comparator", "randomization_method",
            "blinding", "follow_up_duration",
            "primary_outcome", "primary_result",
            "secondary_outcomes", "adverse_events",
            "attrition_rate", "itt_analysis",
        ]
    elif study_type == "cohort":
        type_fields = [
            "exposure", "comparator", "follow_up_duration",
            "primary_outcome", "adjustment_variables",
            "effect_measure", "effect_estimate", "ci_95",
            "p_value", "loss_to_follow_up",
        ]
    else:
        type_fields = [
            "exposure", "outcome", "adjustment_variables",
            "effect_measure", "effect_estimate", "ci_95",
        ]

    all_fields = base_fields + type_fields
    if custom_fields:
        all_fields.extend(custom_fields)

    template = pd.DataFrame(columns=all_fields)
    print(f"Extraction template ({study_type}): {len(all_fields)} fields")
    return template
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/search_strategy.json` | JSON |
| `results/screening_records.csv` | CSV |
| `results/risk_of_bias.csv` | CSV |
| `results/data_extraction.csv` | CSV |
| `figures/prisma_flow.mmd` | Mermaid |
| `figures/prisma_flow.svg` | SVG |

### 利用可能ツール

> PubMed/EuropePMC ツールは `scientific-literature-search` スキルと共有。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| PubMed | `PubMed_search_articles` | 系統的検索 |
| PubMed | `PubMed_Guidelines_Search` | ガイドライン検索 |
| EuropePMC | `EuropePMC_search_articles` | 欧州文献検索 |

### 参照スキル

| スキル | 関連 |
|---|---|
| `scientific-literature-search` | マルチ DB 検索実行 |
| `scientific-meta-analysis` | 量的統合 (Forest/Funnel プロット) |
| `scientific-critical-review` | 品質評価・批判レビュー |
| `scientific-academic-writing` | レビュー論文執筆 |
| `scientific-scientific-schematics` | PRISMA 図作成 |

### 依存パッケージ

`pandas`, `json` (stdlib)

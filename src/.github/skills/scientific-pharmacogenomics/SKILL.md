---
name: scientific-pharmacogenomics
description: |
  ファーマコゲノミクス (薬理ゲノム学) 解析スキル。PharmGKB/ClinPGx による
  遺伝子-薬物相互作用照会、CPIC ガイドライン取得・解釈、Star アレル分類、
  代謝酵素表現型判定 (PM/IM/NM/RM/UM)、FDA 薬理ゲノムバイオマーカー、
  投与量レコメンデーション、PGx レポート生成を統合した
  個別化薬物療法支援パイプライン。
---

# Scientific Pharmacogenomics

遺伝子型に基づく薬物選択・用量調整を支援する
ファーマコゲノミクスパイプライン。CPIC/DPWG ガイドライン・PharmGKB・
FDA PGx バイオマーカーのデータベース照会と解釈を統合。

## When to Use

- 遺伝子型に基づく薬物選択・投与量調整が必要なとき
- CPIC/DPWG ガイドラインの系統的照会が必要なとき
- Star アレル分類・代謝酵素表現型 (PM/IM/NM/RM/UM) を判定するとき
- FDA 承認ファーマコゲノミクスバイオマーカーを確認するとき
- 個別化薬物療法レポートを生成するとき

---

## Quick Start

## 1. 遺伝子-薬物相互作用照会

```python
import pandas as pd
import json


def query_gene_drug_interactions(gene_symbol, data_source="PharmGKB"):
    """
    遺伝子-薬物相互作用の照会。

    主要 PGx 遺伝子:
    - CYP2D6: コデイン, タモキシフェン, トラマドール
    - CYP2C19: クロピドグレル, エスシタロプラム, オメプラゾール
    - CYP2C9: ワルファリン, フェニトイン
    - CYP3A5: タクロリムス
    - DPYD: フルオロウラシル, カペシタビン
    - TPMT/NUDT15: アザチオプリン, 6-MP
    - UGT1A1: イリノテカン
    - SLCO1B1: シンバスタチン
    - HLA-B: アバカビル (*57:01), カルバマゼピン (*15:02)
    - VKORC1: ワルファリン
    """
    print(f"  Querying {data_source} for gene: {gene_symbol}")

    # PharmGKB の主要フィールド
    query_fields = {
        "gene": gene_symbol,
        "clinical_annotations": "Tier 1A-3",
        "drug_labels": "FDA/EMA/PMDA/HCSC",
        "guideline_annotations": "CPIC/DPWG/CPNDS",
        "variant_annotations": "Clinical significance",
        "pathway_annotations": "PK/PD pathways",
    }

    return query_fields


def get_cpic_guidelines(gene_symbol=None, drug_name=None):
    """
    CPIC ガイドラインの取得。

    CPIC Level:
    - Level A: 処方変更を義務付ける強力なエビデンス
    - Level A/B: 処方変更を推奨するエビデンス
    - Level B: 考慮すべきエビデンス
    - Level C: 情報提供レベル
    - Level D: 有用なデータなし
    """
    print(f"  Querying CPIC guidelines:")
    if gene_symbol:
        print(f"    Gene: {gene_symbol}")
    if drug_name:
        print(f"    Drug: {drug_name}")

    return {"gene": gene_symbol, "drug": drug_name, "level": "pending"}
```

## 2. Star アレル分類・表現型判定

```python
import pandas as pd
import numpy as np


def star_allele_annotation(gene, genotype_variants):
    """
    Star アレル分類パイプライン。

    PharmCAT (Pharmacogenomics Clinical Annotation Tool) 互換の
    Star アレルコール → 表現型変換。

    Parameters:
        gene: 対象遺伝子 (e.g., "CYP2D6")
        genotype_variants: {rsID: genotype} 辞書
    """
    # CYP2D6 Star アレル例
    cyp2d6_alleles = {
        "*1": {"function": "Normal", "activity_score": 1.0},
        "*2": {"function": "Normal", "activity_score": 1.0},
        "*3": {"function": "No function", "activity_score": 0.0},
        "*4": {"function": "No function", "activity_score": 0.0},
        "*5": {"function": "No function (gene deletion)", "activity_score": 0.0},
        "*6": {"function": "No function", "activity_score": 0.0},
        "*9": {"function": "Decreased", "activity_score": 0.5},
        "*10": {"function": "Decreased", "activity_score": 0.25},
        "*17": {"function": "Decreased", "activity_score": 0.5},
        "*41": {"function": "Decreased", "activity_score": 0.5},
    }

    print(f"  Gene: {gene}")
    print(f"  Variants provided: {len(genotype_variants)}")
    print(f"  Reference: PharmCAT + PharmVar")

    return cyp2d6_alleles


def determine_metabolizer_phenotype(gene, diplotype, activity_scores):
    """
    Activity Score ベースの代謝酵素表現型判定。

    表現型分類:
    - PM (Poor Metabolizer): AS = 0
    - IM (Intermediate Metabolizer): 0 < AS < 1.25
    - NM (Normal Metabolizer): 1.25 ≤ AS ≤ 2.25
    - RM (Rapid Metabolizer): 2.25 < AS ≤ 3.0 (CYP2C19 のみ)
    - UM (Ultra-rapid Metabolizer): AS > 3.0 or gene duplication
    """
    as1, as2 = activity_scores
    total_as = as1 + as2

    if total_as == 0:
        phenotype = "PM (Poor Metabolizer)"
    elif total_as < 1.25:
        phenotype = "IM (Intermediate Metabolizer)"
    elif total_as <= 2.25:
        phenotype = "NM (Normal Metabolizer)"
    elif total_as <= 3.0:
        phenotype = "RM (Rapid Metabolizer)"
    else:
        phenotype = "UM (Ultra-rapid Metabolizer)"

    print(f"  Gene: {gene}")
    print(f"  Diplotype: {diplotype}")
    print(f"  Activity Score: {as1} + {as2} = {total_as}")
    print(f"  Phenotype: {phenotype}")

    # CPIC 投与量レコメンデーション
    if gene == "CYP2D6" and "PM" in phenotype:
        print("  ⚠ CPIC: コデイン禁忌 (モルフィンへの変換不可)")
        print("  ⚠ CPIC: タモキシフェン → 代替薬を推奨")

    return {"gene": gene, "diplotype": diplotype,
            "activity_score": total_as, "phenotype": phenotype}
```

## 3. FDA PGx バイオマーカー照会

```python
import pandas as pd


def query_fda_pgx_biomarkers(drug_name=None, gene_name=None,
                               biomarker_type=None):
    """
    FDA 承認ファーマコゲノミクスバイオマーカーの照会。

    FDA PGx Labeling Categories:
    - Required: テスト必須 (e.g., HLA-B*57:01 for abacavir)
    - Recommended: テスト推奨
    - Actionable: PGx 情報あり
    - Informative: 参考情報

    300+ の遺伝子-薬物ペアが FDA ラベルに記載。
    """
    print("  Querying FDA Pharmacogenomic Biomarkers:")
    if drug_name:
        print(f"    Drug: {drug_name}")
    if gene_name:
        print(f"    Gene: {gene_name}")
    if biomarker_type:
        print(f"    Type: {biomarker_type}")

    # FDA 主要バイオマーカー例
    key_biomarkers = [
        {"gene": "HLA-B*57:01", "drug": "Abacavir", "action": "Required",
         "recommendation": "HLA-B*57:01 陽性 → 禁忌 (過敏反応リスク)"},
        {"gene": "CYP2C19", "drug": "Clopidogrel", "action": "Actionable",
         "recommendation": "PM → 代替抗血小板薬 (ticagrelor/prasugrel)"},
        {"gene": "DPYD", "drug": "5-FU/Capecitabine", "action": "Recommended",
         "recommendation": "PM → 用量 50% 減量 or 代替薬"},
        {"gene": "UGT1A1*28", "drug": "Irinotecan", "action": "Recommended",
         "recommendation": "TA7/TA7 → 初回投与量減量"},
        {"gene": "TPMT/NUDT15", "drug": "Azathioprine", "action": "Recommended",
         "recommendation": "PM → 10%用量 or 代替薬"},
    ]

    return pd.DataFrame(key_biomarkers)


def pgx_dosing_recommendation(gene, phenotype, drug):
    """
    表現型に基づく CPIC 投与量レコメンデーション生成。
    """
    # CPIC 投与量テーブル例 (CYP2C19 × Clopidogrel)
    cpic_table = {
        ("CYP2C19", "UM", "Clopidogrel"): "標準用量 75 mg/day",
        ("CYP2C19", "RM", "Clopidogrel"): "標準用量 75 mg/day",
        ("CYP2C19", "NM", "Clopidogrel"): "標準用量 75 mg/day",
        ("CYP2C19", "IM", "Clopidogrel"): "代替抗血小板薬を推奨 (ticagrelor/prasugrel)",
        ("CYP2C19", "PM", "Clopidogrel"): "代替抗血小板薬を推奨 (ticagrelor/prasugrel)",
    }

    key = (gene, phenotype.split(" ")[0], drug)
    recommendation = cpic_table.get(key, "ガイドライン情報なし")

    result = {
        "gene": gene,
        "phenotype": phenotype,
        "drug": drug,
        "recommendation": recommendation,
        "source": "CPIC",
        "evidence_level": "Level A",
    }

    print(f"  PGx Dosing Recommendation:")
    print(f"    Gene: {gene} | Phenotype: {phenotype}")
    print(f"    Drug: {drug}")
    print(f"    Recommendation: {recommendation}")

    return result
```

## 4. PGx レポート生成

```python
import json
import pandas as pd
from datetime import datetime


def generate_pgx_report(patient_results, output_file="results/pgx_report.json"):
    """
    包括的 PGx レポート生成。

    含まれる情報:
    - 患者遺伝子型サマリ
    - Star アレル → 表現型マッピング
    - 各薬物の CPIC/DPWG レコメンデーション
    - FDA バイオマーカーステータス
    - アクショナブル所見ハイライト
    """
    import os
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    report = {
        "report_type": "Pharmacogenomics Report",
        "generated_at": datetime.now().isoformat(),
        "patient_id": patient_results.get("patient_id", "anonymous"),
        "genes_tested": [],
        "actionable_findings": [],
        "drug_recommendations": [],
    }

    for gene_result in patient_results.get("gene_results", []):
        gene_entry = {
            "gene": gene_result["gene"],
            "diplotype": gene_result["diplotype"],
            "phenotype": gene_result["phenotype"],
            "activity_score": gene_result.get("activity_score"),
        }
        report["genes_tested"].append(gene_entry)

        # アクショナブル所見の抽出
        if "PM" in gene_result["phenotype"] or "UM" in gene_result["phenotype"]:
            report["actionable_findings"].append({
                "gene": gene_result["gene"],
                "phenotype": gene_result["phenotype"],
                "clinical_significance": "Actionable — 投与量調整/代替薬検討を推奨",
            })

    with open(output_file, "w") as f:
        json.dump(report, f, ensure_ascii=False, indent=2)

    print(f"  PGx Report generated: {output_file}")
    print(f"  Genes tested: {len(report['genes_tested'])}")
    print(f"  Actionable findings: {len(report['actionable_findings'])}")

    return report
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/pgx_report.json` | JSON |
| `results/gene_drug_interactions.csv` | CSV |
| `results/star_allele_calls.csv` | CSV |
| `results/dosing_recommendations.csv` | CSV |
| `figures/pgx_phenotype_summary.png` | PNG |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| PharmGKB | `PharmGKB_search_drugs` | 薬物検索 (遺伝子関連) |
| PharmGKB | `PharmGKB_search_genes` | PGx 遺伝子検索 |
| PharmGKB | `PharmGKB_get_drug_details` | 薬物詳細 (クロスリファレンス) |
| PharmGKB | `PharmGKB_get_gene_details` | 遺伝子詳細 (アレル情報) |
| PharmGKB | `PharmGKB_get_clinical_annotations` | 遺伝子-薬物臨床アノテーション |
| PharmGKB | `PharmGKB_get_dosing_guidelines` | CPIC/DPWG 用量ガイドライン |
| PharmGKB | `PharmGKB_search_variants` | 遺伝的変異検索 |
| FDA | `fda_pharmacogenomic_biomarkers` | FDA PGx バイオマーカー一覧 |
| FDA | `FDA_get_pharmacogenomics_info_by_drug_name` | 薬物名で PGx 情報取得 |
| FDA | `FDA_get_drug_name_by_pharmacogenomics` | PGx から薬物名取得 |
| OpenTargets | `OpenTargets_drug_pharmacogenomics_data` | 薬物 PGx データ |

### 参照スキル

| スキル | 関連 |
|---|---|
| `scientific-variant-interpretation` | ACMG/AMP バリアント解釈 |
| `scientific-pharmacovigilance` | 市販後安全性監視 |
| `scientific-clinical-decision-support` | 臨床意思決定 |
| `scientific-precision-oncology` | 腫瘍 PGx (OncoKB) |
| `scientific-population-genetics` | 集団間アレル頻度差 |

### 依存パッケージ

`pandas`, `numpy`, `json`

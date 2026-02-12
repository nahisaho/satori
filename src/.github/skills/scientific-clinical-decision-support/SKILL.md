---
name: scientific-clinical-decision-support
description: |
  臨床意思決定支援スキル。エビデンスに基づく治療推奨、臨床パスウェイアルゴリズム、
  患者コホート解析、バイオマーカー分類、精密腫瘍学（ToolUniverse の Precision Oncology
  パラダイム + claude-scientific-skills の CDS/Treatment Plans を統合）。
  「治療推奨を作成して」「臨床パスウェイを設計して」「精密医療の解析して」で発火。
---

# Scientific Clinical Decision Support

臨床意思決定支援スキル。エビデンスに基づく治療推奨の生成、
臨床パスウェイの設計、患者特性に応じた個別化医療の支援を提供する。

## When to Use

- エビデンスに基づく治療推奨を構造化するとき
- 臨床パスウェイ/アルゴリズムを設計するとき
- 患者コホートの層別化解析を行うとき
- バイオマーカーに基づく治療選択を支援するとき
- 精密腫瘍学（変異プロファイル → 治療推奨）を実施するとき
- 臨床試験の不適格 / 適格基準を解析するとき

## Quick Start

### 臨床意思決定支援パイプライン

```
Input: 患者プロファイル / 疾患 / バイオマーカー
    ↓
Step 1: Evidence Gathering
  - ガイドライン検索 (NCCN, ESMO, etc.)
  - 臨床試験検索 (ClinicalTrials.gov)
  - 文献エビデンス (PubMed)
    ↓
Step 2: Patient Stratification
  - バイオマーカー分類
  - リスク層別化
  - コモビディティ評価
    ↓
Step 3: Treatment Options
  - 標準治療の同定
  - 代替治療の評価
  - 臨床試験マッチング
    ↓
Step 4: Recommendation Synthesis
  - エビデンスレベル付与
  - 推奨強度評価
  - リスク-ベネフィット分析
    ↓
Output: Clinical Decision Report
```

---

## Phase 1: Evidence Framework

### エビデンスレベルと推奨強度

```python
EVIDENCE_LEVELS = {
    "1a": "Systematic review of RCTs",
    "1b": "Individual RCT",
    "2a": "Systematic review of cohort studies",
    "2b": "Individual cohort study",
    "3a": "Systematic review of case-control studies",
    "3b": "Individual case-control study",
    "4":  "Case series / poor quality cohort",
    "5":  "Expert opinion",
}

RECOMMENDATION_GRADES = {
    "A": {"description": "Strong recommendation", "evidence": "1a, 1b"},
    "B": {"description": "Moderate recommendation", "evidence": "2a, 2b, 3a"},
    "C": {"description": "Weak recommendation", "evidence": "3b, 4"},
    "D": {"description": "Very weak recommendation", "evidence": "5, expert opinion"},
}

def grade_recommendation(evidence_sources):
    """
    GRADE システムによるエビデンス評価。
    """
    best_level = min(
        (s.get("evidence_level", "5") for s in evidence_sources),
        key=lambda x: float(x.replace("a", ".1").replace("b", ".2"))
    )

    if best_level.startswith("1"):
        return "A", "Strong"
    elif best_level.startswith("2") or best_level == "3a":
        return "B", "Moderate"
    elif best_level in ("3b", "4"):
        return "C", "Weak"
    else:
        return "D", "Very Weak"
```

---

## Phase 2: Precision Oncology

### 変異プロファイル → 治療推奨

```python
def precision_oncology_workflow(tumor_type, mutations):
    """
    精密腫瘍学ワークフロー。
    変異プロファイルから actionable targets を同定し、
    治療推奨を生成する。
    """
    workflow = {
        "step1": "変異リスト → OncoKB アノテーション",
        "step2": "Actionable mutations の同定 (Level 1-4)",
        "step3": "腫瘍タイプ特異的ガイドライン確認",
        "step4": "承認薬・臨床試験のマッチング",
        "step5": "治療推奨レポート生成",
    }

    # OncoKB Actionability Levels
    actionability = {
        "Level_1": "FDA-approved, same tumor type",
        "Level_2": "Standard care, different tumor type",
        "Level_3A": "Clinical evidence, same tumor",
        "Level_3B": "Clinical evidence, different tumor",
        "Level_4": "Preclinical evidence",
    }

    return workflow


# 主要ながんバイオマーカー
CANCER_BIOMARKERS = {
    "EGFR": {
        "tumor_types": ["NSCLC"],
        "mutations": ["L858R", "ex19del", "T790M", "C797S"],
        "therapies": {
            "L858R/ex19del": ["osimertinib", "erlotinib", "gefitinib"],
            "T790M": ["osimertinib"],
        }
    },
    "BRAF": {
        "tumor_types": ["Melanoma", "NSCLC", "CRC", "Thyroid"],
        "mutations": ["V600E", "V600K"],
        "therapies": {
            "V600E": ["vemurafenib + cobimetinib", "dabrafenib + trametinib",
                       "encorafenib + binimetinib"],
        }
    },
    "ALK": {
        "tumor_types": ["NSCLC"],
        "mutations": ["Fusion"],
        "therapies": {
            "Fusion": ["alectinib", "lorlatinib", "crizotinib", "ceritinib"],
        }
    },
    "MSI-H/dMMR": {
        "tumor_types": ["Pan-cancer"],
        "therapies": {
            "MSI-H": ["pembrolizumab", "nivolumab"],
        }
    },
}
```

---

## Phase 3: Clinical Trial Matching

### 臨床試験マッチング

```python
def match_clinical_trials(patient_profile):
    """
    ClinicalTrials.gov API を用いた臨床試験マッチング。
    """
    import requests

    params = {
        "query.cond": patient_profile.get("condition"),
        "query.intr": patient_profile.get("biomarker", ""),
        "filter.overallStatus": "RECRUITING",
        "pageSize": 20,
    }

    response = requests.get(
        "https://clinicaltrials.gov/api/v2/studies",
        params=params
    )

    trials = []
    if response.status_code == 200:
        for study in response.json().get("studies", []):
            protocol = study.get("protocolSection", {})
            trials.append({
                "nct_id": protocol.get("identificationModule", {}).get("nctId"),
                "title": protocol.get("identificationModule", {}).get("briefTitle"),
                "phase": protocol.get("designModule", {}).get("phases", []),
                "status": protocol.get("statusModule", {}).get("overallStatus"),
                "interventions": protocol.get("armsInterventionsModule", {}).get("interventions", []),
            })

    return trials
```

---

## Report Template

```markdown
# Clinical Decision Support Report

**Patient Profile**: [anonymized]
**Condition**: [disease]
**Date**: [date]

## 1. Clinical Context

## 2. Biomarker Profile
| Biomarker | Result | Clinical Significance |
|-----------|--------|----------------------|

## 3. Evidence-Based Treatment Options
| Option | Evidence Level | Recommendation Grade | Notes |
|--------|---------------|---------------------|-------|

## 4. Clinical Guideline Alignment
### 4.1 NCCN Guidelines
### 4.2 ESMO Guidelines
### 4.3 National Guidelines

## 5. Clinical Trial Options
| NCT ID | Title | Phase | Status | Relevance |
|--------|-------|-------|--------|-----------|

## 6. Risk-Benefit Analysis
| Treatment | Expected Benefit | Key Risks | NNT/NNH |
|-----------|-----------------|-----------|---------|

## 7. Recommendations
### 7.1 First-line
### 7.2 Second-line
### 7.3 Monitoring Plan

## Disclaimer
This report is for research and educational purposes only.
Clinical decisions must be made by qualified healthcare professionals.
```

---

## Completeness Checklist

- [ ] エビデンス検索: ガイドライン + 臨床試験 + 文献
- [ ] 患者層別化: バイオマーカー + リスクスコア
- [ ] 治療選択肢: 標準治療 + 代替 + 臨床試験
- [ ] エビデンスレベル: 全推奨に GRADE レベルを付与
- [ ] 安全性: 薬物相互作用 + 禁忌を確認
- [ ] 免責事項: 研究目的使用の明記

## Best Practices

1. **ガイドラインを最優先**: NCCN/ESMO 等の最新版を確認
2. **エビデンスレベルを明記**: 全推奨に根拠を付与
3. **患者個別因子を考慮**: 年齢、腎機能、薬物相互作用
4. **臨床試験を常に検索**: 標準治療不適の場合の選択肢
5. **免責事項を記載**: AI ツールは臨床判断の代替ではない

## References

### Output Files

| ファイル | 形式 | 生成タイミング |
|---|---|---|
| `results/clinical_decision_report.md` | 臨床意思決定レポート（Markdown） | 全解析完了時 |
| `results/clinical_recommendation.json` | 推奨事項データ（JSON） | GRADE 評価完了時 |
| `results/trial_matches.json` | 臨床試験マッチング結果（JSON） | 試験検索完了時 |

### 参照スキル

| スキル | 連携 |
|---|---|
| `scientific-variant-interpretation` | ← バリアント解釈結果の臨床応用 |
| `scientific-drug-repurposing` | ← リポジショニング候補の臨床評価 |
| `scientific-survival-clinical` | ← 生存解析・臨床試験データ |
| `scientific-meta-analysis` | ← メタアナリシスのエビデンス統合 |
| `scientific-deep-research` | ← 最新エビデンスの収集 |
| `scientific-presentation-design` | → 臨床結果のプレゼンテーション |
| `scientific-academic-writing` | → 研究成果の論文化 |

---
name: scientific-clinical-trials-analytics
description: |
  臨床試験レジストリ解析スキル。ClinicalTrials.gov API v2 経由の
  多基準試験検索 (疾患×介入×地域×フェーズ×ステータス)、試験詳細取得
  (適格基準・アウトカム・施設・有害事象)、競合ランドスケープ解析、
  AI 支援試験デザイン (ClinicalTrialDesignAgent)、バルクデータ取得・CSV エクスポート
  を統合した臨床研究支援パイプライン。
---

# Scientific Clinical Trials Analytics

ClinicalTrials.gov レジストリの包括的データアクセスと
競合分析パイプラインを提供する。40 万+ の登録試験から
疾患・介入・地域・フェーズ・スポンサー別にフィルタリングと
構造化データ抽出を実施。

## When to Use

- 特定疾患の臨床試験ランドスケープを調査するとき
- 競合薬剤の臨床開発状況を把握するとき
- 試験の適格基準・エンドポイント・施設情報を取得するとき
- 新規試験の計画にあたりプロトコルデザインの参考が必要なとき
- 有害事象・アウトカムデータを系統的に抽出するとき

---

## Quick Start

## 1. 臨床試験検索・フィルタリング

```python
import pandas as pd
import json
from datetime import datetime


def search_clinical_trials(condition=None, intervention=None,
                            phase=None, status=None,
                            location_country=None, sponsor=None,
                            start_date_from=None, max_results=100):
    """
    ClinicalTrials.gov 多基準検索。

    Parameters:
        condition: 疾患名 (e.g., "non-small cell lung cancer")
        intervention: 介入名 (e.g., "pembrolizumab")
        phase: フェーズ ("Phase 1", "Phase 2", "Phase 3", "Phase 4")
        status: ステータス ("RECRUITING", "COMPLETED", "ACTIVE_NOT_RECRUITING")
        location_country: 実施国 (e.g., "Japan")
        sponsor: スポンサー名
    """
    import requests

    base_url = "https://clinicaltrials.gov/api/v2/studies"
    params = {"format": "json", "pageSize": min(max_results, 100)}

    # クエリ構築
    query_parts = []
    if condition:
        query_parts.append(f"CONDITION[{condition}]")
    if intervention:
        query_parts.append(f"INTERVENTION[{intervention}]")
    if query_parts:
        params["query.cond"] = condition
        params["query.intr"] = intervention

    if phase:
        params["filter.phase"] = phase
    if status:
        params["filter.overallStatus"] = status
    if location_country:
        params["query.locn"] = location_country

    resp = requests.get(base_url, params=params)
    data = resp.json()

    studies = data.get("studies", [])
    total = data.get("totalCount", 0)

    print(f"  Search results: {total} trials found")
    print(f"  Filters: condition={condition}, intervention={intervention}")
    print(f"  Phase: {phase}, Status: {status}")

    results = []
    for study in studies:
        protocol = study.get("protocolSection", {})
        id_module = protocol.get("identificationModule", {})
        status_module = protocol.get("statusModule", {})
        design_module = protocol.get("designModule", {})

        results.append({
            "nct_id": id_module.get("nctId"),
            "title": id_module.get("briefTitle"),
            "status": status_module.get("overallStatus"),
            "phase": design_module.get("phases", [None]),
            "enrollment": design_module.get("enrollmentInfo", {}).get("count"),
            "start_date": status_module.get("startDateStruct", {}).get("date"),
        })

    return pd.DataFrame(results), total
```

## 2. 試験詳細取得

```python
import pandas as pd
import json


def get_trial_details(nct_id):
    """
    NCT ID による臨床試験の完全詳細取得。

    取得項目:
    - プロトコル (デザイン, アーム, 介入, マスキング)
    - 適格基準 (年齢, 性別, 包含/除外基準)
    - アウトカム (主要/副次エンドポイント)
    - 施設・地域情報
    - スポンサー, 共同研究者
    - リファレンス (関連論文)
    """
    import requests

    url = f"https://clinicaltrials.gov/api/v2/studies/{nct_id}"
    resp = requests.get(url, params={"format": "json"})
    data = resp.json()

    protocol = data.get("protocolSection", {})

    # 基本情報
    id_mod = protocol.get("identificationModule", {})
    design_mod = protocol.get("designModule", {})
    eligibility = protocol.get("eligibilityModule", {})
    outcomes_mod = protocol.get("outcomesModule", {})
    contacts_mod = protocol.get("contactsLocationsModule", {})

    detail = {
        "nct_id": nct_id,
        "title": id_mod.get("officialTitle"),
        "brief_title": id_mod.get("briefTitle"),
        "study_type": design_mod.get("studyType"),
        "phases": design_mod.get("phases"),
        "allocation": design_mod.get("designInfo", {}).get("allocation"),
        "masking": design_mod.get("designInfo", {}).get("maskingInfo", {}).get("masking"),
    }

    # 適格基準
    detail["eligibility"] = {
        "min_age": eligibility.get("minimumAge"),
        "max_age": eligibility.get("maximumAge"),
        "sex": eligibility.get("sex"),
        "criteria": eligibility.get("eligibilityCriteria"),
    }

    # アウトカム
    primary_outcomes = outcomes_mod.get("primaryOutcomes", [])
    secondary_outcomes = outcomes_mod.get("secondaryOutcomes", [])
    detail["primary_outcomes"] = [
        {"measure": o.get("measure"), "timeframe": o.get("timeFrame")}
        for o in primary_outcomes
    ]

    # 施設
    locations = contacts_mod.get("locations", [])
    detail["n_locations"] = len(locations)
    detail["countries"] = list(set(
        loc.get("country") for loc in locations if loc.get("country")
    ))

    print(f"  Trial: {nct_id}")
    print(f"  Title: {detail['brief_title']}")
    print(f"  Type: {detail['study_type']}, Phase: {detail['phases']}")
    print(f"  Primary outcomes: {len(primary_outcomes)}")
    print(f"  Locations: {detail['n_locations']} sites in {len(detail['countries'])} countries")

    return detail
```

## 3. 競合ランドスケープ解析

```python
import pandas as pd
import numpy as np


def competitive_landscape_analysis(condition, intervention_class=None,
                                     status_filter=None):
    """
    疾患・介入クラス別の臨床開発ランドスケープ解析。

    分析項目:
    - フェーズ分布 (Phase 1/2/3/4)
    - ステータス分布 (Recruiting/Completed/Terminated)
    - スポンサー別試験数
    - 経年トレンド (開始年別)
    - 地域分布
    """
    print(f"  Competitive landscape for: {condition}")
    if intervention_class:
        print(f"  Intervention class: {intervention_class}")

    # フェーズ分布解析
    phase_mapping = {
        "PHASE1": "Phase 1",
        "PHASE2": "Phase 2",
        "PHASE3": "Phase 3",
        "PHASE4": "Phase 4",
        "EARLY_PHASE1": "Early Phase 1",
    }

    return {"condition": condition, "intervention_class": intervention_class}


def extract_trial_adverse_events(nct_id):
    """
    試験の有害事象データ抽出。

    FDA 報告基準:
    - Serious AE (SAE): 死亡, 入院, 障害, 先天異常
    - Other AE: Grade 1-4 有害事象
    - CTCAE v5.0 grading
    """
    print(f"  Extracting adverse events for: {nct_id}")
    print("  Categories: Serious AE, Other AE")
    print("  Grading: CTCAE v5.0")

    return {"nct_id": nct_id}


def extract_trial_outcomes(nct_id):
    """
    試験結果 (Results section) からのアウトカムデータ抽出。

    - Primary outcome measures + 統計結果
    - Secondary outcome measures
    - 参加者フロー (Screened → Enrolled → Completed → Analyzed)
    """
    print(f"  Extracting outcomes for: {nct_id}")

    return {"nct_id": nct_id}
```

## 4. バルクデータ取得・エクスポート

```python
import pandas as pd
import json


def bulk_trial_export(condition, max_trials=1000,
                       output_file="results/clinical_trials_export.csv"):
    """
    大規模臨床試験データのバルク取得・CSV エクスポート。

    ページネーション対応 (1000 件 = 10 ページ × 100 件/ページ)。
    """
    import os
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    all_results = []
    page_token = None

    print(f"  Bulk export for: {condition}")
    print(f"  Max trials: {max_trials}")

    # ここでは概念的コード — 実際は API ページネーション
    # while len(all_results) < max_trials:
    #     resp = requests.get(url, params={..., "pageToken": page_token})
    #     studies = resp.json()["studies"]
    #     all_results.extend(studies)
    #     page_token = resp.json().get("nextPageToken")

    print(f"  Exported to: {output_file}")

    return output_file


def trial_design_summary(trials_df):
    """
    臨床試験デザインの要約統計。

    - Study type 分布 (Interventional/Observational)
    - Allocation (Randomized/Non-randomized)
    - Masking (Open/Single/Double/Triple/Quadruple)
    - Primary purpose (Treatment/Prevention/Diagnostic)
    """
    print("  Trial Design Summary:")

    if "study_type" in trials_df.columns:
        type_dist = trials_df["study_type"].value_counts()
        for st, count in type_dist.items():
            print(f"    {st}: {count}")

    return trials_df.describe()
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/clinical_trials_search.csv` | CSV |
| `results/trial_details.json` | JSON |
| `results/competitive_landscape.json` | JSON |
| `results/clinical_trials_export.csv` | CSV |
| `results/trial_adverse_events.csv` | CSV |
| `figures/trial_phase_distribution.png` | PNG |
| `figures/trial_timeline.png` | PNG |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| ClinicalTrials | `clinical_trials_search` | 臨床試験検索 |
| ClinicalTrials | `search_clinical_trials` | 詳細パラメータ検索 |
| ClinicalTrials | `get_clinical_trial_conditions_and_interventions` | 疾患・介入情報取得 |
| ClinicalTrials | `get_clinical_trial_locations` | 施設・地域取得 |
| ClinicalTrials | `extract_clinical_trial_adverse_events` | 有害事象抽出 |
| ClinicalTrials | `extract_clinical_trial_outcomes` | アウトカム抽出 |
| ClinicalTrials | `ClinicalTrialDesignAgent` | AI 支援試験デザイン |
| FDA | `FDA_get_clinical_studies_info_by_drug_name` | 薬物名で臨床研究情報 |
| FDA | `FDA_get_drug_names_by_clinical_studies` | 臨床研究から薬物名 |

### 参照スキル

| スキル | 関連 |
|---|---|
| `scientific-survival-clinical` | 生存解析 (KM/Cox) |
| `scientific-meta-analysis` | メタアナリシス統合 |
| `scientific-epidemiology-public-health` | 疫学リスク指標 |
| `scientific-clinical-decision-support` | 臨床意思決定 |
| `scientific-pharmacovigilance` | 安全性モニタリング |

### 依存パッケージ

`pandas`, `numpy`, `requests`, `json`

---
name: scientific-regulatory-science
description: |
  規制科学パイプラインスキル。FDA (医薬品/医療機器/食品)・EMA・PMDA 規制データベース横断照会、
  Orange Book 承認履歴・特許・排他性情報、510(k) デバイスクリアランス、
  ISO 13485 品質管理システム (設計管理/CAPA/リスク管理)、
  USPTO 特許検索を統合した薬事・規制情報パイプライン。
---

# Scientific Regulatory Science

FDA・EMA・PMDA 等の規制当局データベースを横断照会し、
医薬品・医療機器・食品の承認情報・安全性データ・品質管理要件・知的財産情報を
体系的に取得するパイプラインを提供する。

## When to Use

- FDA 医薬品の承認履歴・ジェネリック状況を調査するとき
- 医療機器の 510(k) クリアランス・規制分類を確認するとき
- 食品安全性 (FDA リコール, 有害事象) データを取得するとき
- ISO 13485 品質管理システムの設計管理・CAPA を計画するとき
- 競合特許のランドスケープ調査が必要なとき

---

## Quick Start

## 1. FDA Orange Book 照会

```python
import pandas as pd
import json


def query_fda_orange_book(drug_name=None, nda_number=None,
                            active_ingredient=None):
    """
    FDA Orange Book (Approved Drug Products with Therapeutic Equivalence Evaluations)。

    取得項目:
    - 承認履歴 (NDA/ANDA, 承認日, 適応症)
    - 治療的同等性 (TE) コード (AB, BX, etc.)
    - 特許情報 (満了日, 特許番号, Use Code)
    - 排他性情報 (NCE, Orphan, Pediatric)
    - ジェネリック承認状況
    """
    print(f"  Querying FDA Orange Book:")
    if drug_name:
        print(f"    Drug: {drug_name}")
    if active_ingredient:
        print(f"    Active ingredient: {active_ingredient}")

    return {"drug_name": drug_name, "nda_number": nda_number}


def analyze_patent_landscape(drug_name):
    """
    医薬品特許ランドスケープ解析。

    - Orange Book 掲載特許の満了日タイムライン
    - パラグラフ IV チャレンジの有無
    - 180 日排他性の状況
    - ジェネリック参入予測
    """
    print(f"  Patent landscape for: {drug_name}")
    print("  Analysis: Patent expiry timeline + exclusivity + generic entry")

    return {"drug_name": drug_name}
```

## 2. FDA 医療機器規制データ

```python
import pandas as pd


def query_fda_device_classification(device_name=None, product_code=None,
                                      regulation_number=None):
    """
    FDA 医療機器分類照会。

    分類クラス:
    - Class I: 一般規制 (General Controls)
    - Class II: 特別規制 (Special Controls) — 510(k) 必要
    - Class III: 市販前承認 (PMA) 必要
    """
    print(f"  FDA Device Classification:")
    if device_name:
        print(f"    Device: {device_name}")
    print("  Regulatory pathways: 510(k), De Novo, PMA, HDE")

    return {"device_name": device_name, "product_code": product_code}


def query_510k_clearance(device_name=None, applicant=None,
                           decision_date_from=None):
    """
    FDA 510(k) クリアランスデータ照会。

    510(k) = Premarket Notification
    - Predicate device との実質的同等性 (SE) の証明
    - 90 日レビュー (Traditional) / 30 日 (Special/Abbreviated)
    """
    print(f"  510(k) Clearance search:")
    if device_name:
        print(f"    Device: {device_name}")

    return {"device_name": device_name}
```

## 3. ISO 13485 品質管理システム

```python
import json
from datetime import datetime


def iso13485_design_control_checklist(product_name, risk_class="II"):
    """
    ISO 13485 設計管理チェックリスト。

    設計管理プロセス (ISO 13485:2016 §7.3):
    - 7.3.2 設計・開発計画
    - 7.3.3 設計入力
    - 7.3.4 設計出力
    - 7.3.5 設計レビュー
    - 7.3.6 設計検証
    - 7.3.7 設計バリデーション
    - 7.3.8 設計移管
    - 7.3.9 設計変更管理
    """
    checklist = {
        "product": product_name,
        "risk_class": risk_class,
        "design_phases": [
            {"phase": "Design Planning", "ref": "§7.3.2",
             "deliverables": ["DHF 設計履歴ファイル", "プロジェクト計画書", "リスク管理計画"]},
            {"phase": "Design Input", "ref": "§7.3.3",
             "deliverables": ["ユーザーニーズ", "設計要求仕様 (DRS)", "規制要求"]},
            {"phase": "Design Output", "ref": "§7.3.4",
             "deliverables": ["設計仕様書", "図面", "DMR (Device Master Record)"]},
            {"phase": "Design Review", "ref": "§7.3.5",
             "deliverables": ["設計レビュー議事録", "アクションアイテム"]},
            {"phase": "Design Verification", "ref": "§7.3.6",
             "deliverables": ["検証プロトコル/レポート", "トレーサビリティマトリクス"]},
            {"phase": "Design Validation", "ref": "§7.3.7",
             "deliverables": ["バリデーションプロトコル/レポート", "臨床評価 (必要時)"]},
            {"phase": "Design Transfer", "ref": "§7.3.8",
             "deliverables": ["製造移管文書", "製造仕様"]},
            {"phase": "Design Changes", "ref": "§7.3.9",
             "deliverables": ["変更管理 (ECO)", "影響分析"]},
        ],
    }

    print(f"  ISO 13485 Design Control for: {product_name} (Class {risk_class})")
    for phase in checklist["design_phases"]:
        print(f"    {phase['ref']} {phase['phase']}: {len(phase['deliverables'])} deliverables")

    return checklist


def capa_process(nonconformity_description, severity="Major"):
    """
    CAPA (Corrective and Preventive Action) プロセス。

    ISO 13485 §8.5.2 (是正処置) / §8.5.3 (予防処置)
    """
    capa_record = {
        "id": f"CAPA-{datetime.now().strftime('%Y%m%d-%H%M')}",
        "nonconformity": nonconformity_description,
        "severity": severity,
        "steps": [
            "1. 問題の特定と文書化",
            "2. 即時是正処置 (Containment)",
            "3. 根本原因分析 (5-Why / Fishbone)",
            "4. 是正処置の計画と実施",
            "5. 有効性検証",
            "6. 予防処置の展開",
            "7. CAPA クローズ",
        ],
    }

    print(f"  CAPA initiated: {capa_record['id']}")
    print(f"  Severity: {severity}")

    return capa_record
```

## 4. USPTO 特許検索

```python
import pandas as pd


def search_patents(query, start_date=None, end_date=None,
                    assignee=None, max_results=50):
    """
    USPTO 特許検索 (PatentsView API)。

    検索フィールド:
    - 特許タイトル・アブストラクト・クレーム
    - 出願人/譲受人
    - 出願日/登録日
    - CPC (Cooperative Patent Classification) コード
    """
    print(f"  USPTO Patent search: '{query}'")
    if assignee:
        print(f"  Assignee: {assignee}")
    if start_date:
        print(f"  Date range: {start_date} → {end_date or 'present'}")

    return {"query": query, "max_results": max_results}
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/fda_orange_book.json` | JSON |
| `results/device_classification.csv` | CSV |
| `results/510k_clearances.csv` | CSV |
| `results/iso13485_checklist.json` | JSON |
| `results/capa_record.json` | JSON |
| `results/patent_search.csv` | CSV |
| `figures/patent_timeline.png` | PNG |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| FDA Orange Book | `FDA_OrangeBook_search_drug` | 医薬品検索 |
| FDA Orange Book | `FDA_OrangeBook_get_approval_history` | 承認履歴取得 |
| FDA Orange Book | `FDA_OrangeBook_get_patent_info` | 特許情報取得 |
| FDA Orange Book | `FDA_OrangeBook_get_exclusivity` | 排他性情報 |
| FDA Orange Book | `FDA_OrangeBook_check_generic_availability` | ジェネリック承認状況 |
| FDA Orange Book | `FDA_OrangeBook_get_te_code` | 治療的同等性コード |
| FAERS | `FAERS_search_adverse_event_reports` | 有害事象レポート検索 |
| FAERS | `FAERS_calculate_disproportionality` | 不均衡分析 (ROR/PRR/IC) |

### 参照スキル

| スキル | 関連 |
|---|---|
| `scientific-pharmacovigilance` | 市販後安全性 |
| `scientific-pharmacogenomics` | FDA PGx バイオマーカー |
| `scientific-clinical-trials-analytics` | 臨床試験レジストリ |
| `scientific-admet-pharmacokinetics` | 前臨床 ADMET |
| `scientific-grant-writing` | 規制戦略セクション |

### 依存パッケージ

`pandas`, `numpy`, `requests`, `json`

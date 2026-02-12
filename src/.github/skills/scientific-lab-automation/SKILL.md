---
name: scientific-lab-automation
description: |
  実験室自動化・プロトコル管理スキル。PyLabRobot（液体ハンドリング）、
  Protocols.io（プロトコル共有）、Benchling/LabArchives（ELN/LIMS 統合）、
  Opentrons（ロボティクス）による実験自動化と再現性確保を支援。
  claude-scientific-skills の Lab Automation カテゴリを統合。
  「実験プロトコルを作成して」「液体ハンドリングを自動化して」で発火。
---

# Scientific Lab Automation

実験室自動化とプロトコル管理のスキル。液体ハンドリングロボットの制御、
プロトコルの構造化・共有、電子実験ノート（ELN）連携、
再現性のある実験ワークフローの構築を支援する。

## When to Use

- 液体ハンドリングプロトコルを自動化するとき
- 実験プロトコルを標準化・共有するとき
- ELN（電子実験ノート）や LIMS との連携を設定するとき
- 再現性のある実験ワークフローを設計するとき
- ハイスループットスクリーニング (HTS) を計画するとき

## Quick Start

### 実験自動化パイプライン

```
Step 1: Protocol Design
  - 手順の構造化（input/output/parameters）
  - SOP テンプレート作成
  - パラメータの定義
    ↓
Step 2: Automation Script
  - PyLabRobot / Opentrons Python API
  - ウェルプレートレイアウト定義
  - ピペッティングシーケンス
    ↓
Step 3: Validation
  - ドライラン（シミュレーション）
  - キャリブレーション確認
  - エッジケーステスト
    ↓
Step 4: Execution & Logging
  - 自動実行
  - データ自動記録
  - 異常検出・アラート
    ↓
Step 5: Documentation
  - ELN エントリ自動生成
  - Protocols.io 公開
  - データ管理（LIMS 連携）
```

---

## Phase 1: 液体ハンドリング自動化

### PyLabRobot プロトコル

```python
# PyLabRobot: Universal Python interface for liquid handling robots
# Supports: Hamilton STAR/Vantage, Opentrons OT-2/Flex, Tecan, etc.

from pylabrobot.liquid_handling import LiquidHandler
from pylabrobot.resources import Plate, Tip

async def serial_dilution_protocol(lh: LiquidHandler, source_well, dest_plate,
                                     dilution_factor=2, n_dilutions=8, volume_ul=100):
    """
    系列希釈プロトコル。
    PyLabRobot の統一 API でロボット非依存の記述。
    """
    # 希釈バッファーを全ウェルに分注
    buffer_volume = volume_ul * (1 - 1/dilution_factor)
    for i in range(n_dilutions):
        await lh.aspirate(
            resource=lh.deck.get_resource("buffer_trough"),
            volume=buffer_volume,
        )
        await lh.dispense(
            resource=dest_plate[0][i],
            volume=buffer_volume,
        )

    # 原液から系列希釈
    transfer_volume = volume_ul / dilution_factor
    await lh.aspirate(resource=source_well, volume=transfer_volume)
    await lh.dispense(resource=dest_plate[0][0], volume=transfer_volume)

    for i in range(n_dilutions - 1):
        await lh.aspirate(resource=dest_plate[0][i], volume=transfer_volume)
        await lh.dispense(resource=dest_plate[0][i+1], volume=transfer_volume)

    return {"status": "complete", "n_dilutions": n_dilutions}
```

### Opentrons OT-2 プロトコル

```python
from opentrons import protocol_api

metadata = {
    'protocolName': 'Serial Dilution',
    'author': 'SATORI',
    'apiLevel': '2.16',
}

def run(protocol: protocol_api.ProtocolContext):
    # ラビューアの定義
    tiprack = protocol.load_labware('opentrons_96_tiprack_300ul', 1)
    plate = protocol.load_labware('corning_96_wellplate_360ul_flat', 2)
    reservoir = protocol.load_labware('nest_12_reservoir_15ml', 3)
    pipette = protocol.load_instrument('p300_single_gen2', 'left', tip_racks=[tiprack])

    # 希釈バッファー分注
    pipette.distribute(
        100,
        reservoir['A1'],
        plate.columns()[1:12],
    )

    # 系列希釈
    pipette.transfer(
        100,
        plate.columns()[:11],
        plate.columns()[1:12],
        mix_after=(3, 50),
    )
```

---

## Phase 2: プロトコル構造化

### SOP テンプレート

```markdown
## Standard Operating Procedure

### Protocol: [Protocol Name]
**Version**: [X.Y] | **Date**: [YYYY-MM-DD] | **Author**: [name]

### 1. Purpose & Scope

### 2. Materials & Equipment
| Item | Catalog # | Vendor | Quantity |
|------|-----------|--------|----------|

### 3. Reagent Preparation
| Reagent | Final Conc | Stock Conc | Volume | Solvent |
|---------|------------|------------|--------|---------|

### 4. Procedure
| Step | Action | Parameters | Duration | Notes |
|------|--------|------------|----------|-------|
| 1 | | | | |
| 2 | | | | |

### 5. Quality Control
| QC Check | Acceptance Criteria | Method |
|----------|---------------------|--------|

### 6. Data Recording
- [ ] Raw data location: ___
- [ ] Analysis script: ___
- [ ] ELN entry: ___

### 7. Safety
| Hazard | PPE Required | SDS Reference |
|--------|-------------|---------------|

### 8. Revision History
| Version | Date | Changes | Author |
|---------|------|---------|--------|
```

---

## Phase 3: ELN / LIMS 連携

### Protocols.io 連携

```python
def create_protocol_io_entry(protocol_data):
    """
    Protocols.io API でプロトコルを作成。
    """
    import requests

    headers = {
        "Authorization": f"Bearer {PROTOCOLS_IO_TOKEN}",
        "Content-Type": "application/json",
    }

    payload = {
        "title": protocol_data["title"],
        "description": protocol_data["description"],
        "steps": [
            {"description": step, "duration": dur}
            for step, dur in zip(protocol_data["steps"], protocol_data["durations"])
        ],
        "materials": protocol_data.get("materials", []),
    }

    response = requests.post(
        "https://www.protocols.io/api/v4/protocols",
        headers=headers,
        json=payload,
    )
    return response.json()
```

---

## Report Template

```markdown
# Lab Automation Report: [Protocol Name]

**Robot**: [Hamilton STAR / Opentrons OT-2 / etc.]
**Date**: [date]

## 1. Protocol Overview
## 2. Deck Layout
## 3. Execution Log
| Step | Time | Action | Volume | Well | Status |
|------|------|--------|--------|------|--------|

## 4. QC Results
## 5. Data Files Generated
## 6. ELN Entry Reference
```

---

## Completeness Checklist

- [ ] プロトコル構造化: 全ステップが明確に定義
- [ ] パラメータ定義: ボリューム、温度、時間を指定
- [ ] バリデーション: ドライラン/シミュレーション実施
- [ ] QC 基準: 受入基準を事前定義
- [ ] データ管理: 自動記録の設定
- [ ] ドキュメント: SOP + ELN エントリ

## Best Practices

1. **ロボット非依存で設計**: PyLabRobot の統一 API で機種依存を排除
2. **ドライランを必ず実施**: 実液を使う前にシミュレーションで確認
3. **Dead volume を計算**: リザーバの dead volume を考慮した余裕量
4. **ピペッティング精度を検証**: 蛍光色素/重量法で実測値を確認
5. **バージョン管理**: プロトコルの変更を Git で追跡

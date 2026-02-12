---
name: scientific-protein-design
description: |
  タンパク質設計スキル。ESM タンパク質言語モデル、de novo 設計、指向性進化の
  計算的ガイド、安定性予測をカバー。ToolUniverse の Protein Therapeutic Design
  パラダイムと claude-scientific-skills の ESM/Adaptyv スキルを統合。
  「タンパク質を設計して」「ESM で評価して」「安定性を予測して」で発火。
---

# Scientific Protein Design

計算的タンパク質設計スキル。ESM (Evolutionary Scale Modeling) などの
タンパク質言語モデル、de novo バックボーン設計、配列最適化、
安定性・発現性予測を統合したワークフローを提供する。

## When to Use

- タンパク質の de novo 設計を行うとき
- ESM タンパク質言語モデルで配列を評価するとき
- 指向性進化のライブラリ設計をガイドするとき
- 変異体の安定性・機能への影響を予測するとき
- バインダー（抗体/ペプチド/ミニタンパク質）を設計するとき
- 酵素の活性部位を再設計するとき

## Quick Start

### タンパク質設計パイプライン

```
Phase 1: Target Definition
  - 設計目標の定義（バインダー/酵素/足場）
  - ターゲット構造の取得
  - 設計制約の設定
    ↓
Phase 2: Backbone Generation
  - De novo バックボーン生成（RFdiffusion 等）
  - トポロジー指定
  - 機能サイト保存
    ↓
Phase 3: Sequence Design
  - ProteinMPNN 配列設計
  - ESM スコアリング
  - 多様性確保
    ↓
Phase 4: In Silico Validation
  - ESMFold 構造予測
  - pLDDT / pTM 評価
  - 設計-予測一致度 (RMSD)
    ↓
Phase 5: Developability Assessment
  - 安定性予測
  - 発現予測
  - 凝集リスク
  - 最終ランキング
```

---

## Phase 1: ESM タンパク質言語モデル

### ESM-2 / ESM-1v による配列解析

```python
import torch
import esm

def load_esm_model(model_name="esm2_t33_650M_UR50D"):
    """
    ESM-2 モデルのロード。
    モデルサイズ: 8M, 35M, 150M, 650M, 3B, 15B
    推奨: 650M (精度と速度のバランス)
    """
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    model.eval()
    return model, alphabet, batch_converter


def compute_sequence_loglikelihood(model, alphabet, batch_converter, sequence):
    """
    配列の対数尤度を計算。設計品質の指標として使用。
    """
    data = [("protein", sequence)]
    batch_labels, batch_strs, batch_tokens = batch_converter(data)

    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[33])

    logits = results["logits"]
    log_probs = torch.log_softmax(logits, dim=-1)

    # 各位置の正解アミノ酸の対数確率を取得
    token_log_probs = []
    for i, aa in enumerate(sequence):
        aa_idx = alphabet.get_idx(aa)
        token_log_probs.append(log_probs[0, i+1, aa_idx].item())

    return {
        "sequence": sequence,
        "mean_log_likelihood": sum(token_log_probs) / len(token_log_probs),
        "total_log_likelihood": sum(token_log_probs),
        "per_position": token_log_probs,
    }
```

---

## Phase 2: 変異影響予測

### Zero-shot 変異スキャン

```python
def predict_mutation_effect(model, alphabet, batch_converter, wt_sequence, mutations):
    """
    ESM による zero-shot 変異影響予測。
    mutations: [("A", 42, "V"), ("G", 100, "D"), ...]  # (wt_aa, pos, mut_aa)

    出力: 各変異の ΔLL (log-likelihood ratio)
    正 = 有利な変異, 負 = 有害な変異
    """
    wt_data = [("wt", wt_sequence)]
    _, _, wt_tokens = batch_converter(wt_data)

    with torch.no_grad():
        wt_results = model(wt_tokens)

    wt_logits = wt_results["logits"]
    wt_log_probs = torch.log_softmax(wt_logits, dim=-1)

    mutation_effects = []
    for wt_aa, pos, mut_aa in mutations:
        wt_idx = alphabet.get_idx(wt_aa)
        mut_idx = alphabet.get_idx(mut_aa)
        # pos は 0-indexed、tokens は 1-indexed (BOS token)
        delta_ll = (wt_log_probs[0, pos+1, mut_idx] - wt_log_probs[0, pos+1, wt_idx]).item()

        mutation_effects.append({
            "mutation": f"{wt_aa}{pos+1}{mut_aa}",
            "delta_log_likelihood": delta_ll,
            "prediction": "Beneficial" if delta_ll > 0 else "Neutral" if delta_ll > -2 else "Deleterious",
        })

    return sorted(mutation_effects, key=lambda x: x["delta_log_likelihood"], reverse=True)
```

---

## Phase 3: De novo 設計ワークフロー

### 設計タスク別パイプライン

```markdown
## Task-Specific Pipelines

### Binder Design (バインダー設計)
1. ターゲット構造取得 → エピトープ定義
2. RFdiffusion でバックボーン生成 (≥5 candidates)
3. ProteinMPNN で配列設計 (≥8 sequences/backbone)
4. ESMFold で構造検証
5. Interface 品質評価 (pAE at interface)

### Scaffold Design (足場設計)
1. トポロジー指定 (α/β/mixed)
2. サイズ制約 (50-200 residues)
3. 機能サイト保存
4. 安定性最適化

### Enzyme Redesign (酵素再設計)
1. 活性部位ジオメトリ維持
2. 触媒残基保存
3. 基質アクセス確保
4. 折り畳みエネルギー最適化
```

---

## Phase 4: In Silico Validation

### ESMFold 検証基準

```python
def validate_design(designed_sequence, target_structure_path=None):
    """
    設計配列の in silico 検証。
    ESMFold で構造予測し、品質メトリクスを評価。
    """
    validation = {
        "sequence_length": len(designed_sequence),
        "esm_log_likelihood": None,  # mean LL
        "esmfold_plddt": None,       # mean pLDDT
        "esmfold_ptm": None,         # pTM score
        "rmsd_to_target": None,      # if target provided
        "pass_criteria": {},
    }

    # Pass/Fail 基準
    criteria = {
        "pLDDT": {"threshold": 70, "direction": ">"},
        "pTM": {"threshold": 0.5, "direction": ">"},
        "RMSD": {"threshold": 2.0, "direction": "<"},  # Å
    }

    return validation


# 発現系推奨
EXPRESSION_SYSTEMS = {
    "simple_scaffold": {"recommended": "E. coli", "alternative": "Insect cells"},
    "disulfide_containing": {"recommended": "Mammalian", "alternative": "Insect cells"},
    "glycosylated": {"recommended": "Mammalian (CHO/HEK)", "alternative": None},
    "toxic_protein": {"recommended": "Cell-free", "alternative": "Insect cells"},
    "large_complex": {"recommended": "Insect cells (baculovirus)", "alternative": "Mammalian"},
}
```

---

## Report Template

```markdown
# Protein Design Report: [PROJECT NAME]

**Design Type**: [Binder / Scaffold / Enzyme]
**Target**: [target protein]
**Date**: [date]

## 1. Design Objective

## 2. Target Analysis
### 2.1 Structure Used
### 2.2 Design Constraints
### 2.3 Key Residues

## 3. Backbone Candidates
| # | Method | Topology | Size | Score |
|---|--------|----------|------|-------|

## 4. Sequence Designs (Top 10)
| Rank | Backbone | Sequence | MPNN Score | ESM LL | pLDDT | pTM |
|------|----------|----------|------------|--------|-------|-----|

## 5. Validation Results
### 5.1 Structure Prediction
### 5.2 Design-Target RMSD
### 5.3 Interface Quality (if binder)

## 6. Developability
| Metric | Value | Status |
|--------|-------|--------|
| Aggregation risk | | |
| Isoelectric point | | |
| Expression prediction | | |

## 7. Final Candidates
### 7.1 Recommended for testing
### 7.2 Sequences (FASTA format)

## 8. Experimental Recommendations
### 8.1 Expression System
### 8.2 Purification Strategy
### 8.3 Characterization Assays
```

---

## Completeness Checklist

- [ ] ターゲット構造: PDB or AlphaFold を使用
- [ ] バックボーン: ≥5 候補を生成
- [ ] 配列設計: ≥8 配列/バックボーン
- [ ] ESM スコア: 全候補の LL を計算
- [ ] 構造検証: pLDDT >70 & pTM >0.5
- [ ] 発現系推奨: 設計タイプに応じた推奨
- [ ] 最終候補: ≥3 passing designs

## Best Practices

1. **多様性を確保**: 配列空間を広く探索（temperature parameter 調整）
2. **Negative design も考慮**: 望ましくない構造を排除
3. **保守的な変異から始める**: Consensus mutations を優先
4. **実験フィードバック**: 設計-テスト-学習サイクルを計画
5. **スケール意識**: 計算コストと精度のトレードオフを明示

## References

### Output Files

| ファイル | 形式 | 生成タイミング |
|---|---|---|
| `results/design_report.md` | 設計レポート（Markdown） | 全設計完了時 |
| `results/design_candidates.json` | 設計候補データ（JSON） | スクリーニング完了時 |
| `results/esm_scores.json` | ESM スコアデータ（JSON） | 変異スキャン完了時 |

### 参照スキル

| スキル | 連携 |
|---|---|
| `scientific-protein-structure-analysis` | ← ターゲット構造データの提供 |
| `scientific-sequence-analysis` | ← 配列・進化情報 |
| `scientific-lab-automation` | → 設計タンパク質の発現・精製プロトコル |
| `scientific-admet-pharmacokinetics` | → タンパク質治療薬の PK 評価 |
| `scientific-academic-writing` | → 研究成果の論文化 |

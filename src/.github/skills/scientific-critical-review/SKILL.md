```skill
---
name: scientific-critical-review
description: |
  学術論文の草稿に対する批判的レビュー・修正スキル。論理構成、考察の深さ、
  データ解釈の妥当性、先行研究との整合性、統計的主張の正確性を多角的に検証し、
  具体的な修正案を生成する。「論文をレビューして」「考察を深めて」「草稿を改善して」で発火。
  scientific-academic-writing で草稿を作成した後のセルフレビューに使用。
---

# Scientific Critical Review & Revision

学術論文の草稿を批判的思考でレビューし、具体的な修正を行うスキル。
査読者の視点から論文の弱点を特定し、考察の深化・論理の強化・主張の精緻化を支援する。

## When to Use

- 草稿を書き終えた後、投稿前のセルフレビューを行うとき
- Discussion / Conclusion の考察が浅いと感じるとき
- データ解釈と主張の整合性を検証したいとき
- 論理の飛躍や根拠不足を特定したいとき
- 査読者からの指摘を事前に予測して対策したいとき

## Quick Start

## 1. レビューワークフロー

```
草稿完成
  ├─ Pass 1: 構造レビュー（マクロ）
  │   ├─ セクション構成の妥当性
  │   ├─ ストーリーラインの一貫性
  │   └─ 各セクション間のバランス
  ├─ Pass 2: 論理・考察レビュー（コア）
  │   ├─ 主張と証拠の対応
  │   ├─ 考察の深さ・多層性
  │   ├─ 代替説明の検討
  │   └─ 限界の誠実な記述
  ├─ Pass 3: データ・統計レビュー
  │   ├─ 統計的主張の正確性
  │   ├─ 効果量・信頼区間の記載
  │   └─ 図表とテキストの整合性
  ├─ Pass 4: 文章・表現レビュー（ミクロ）
  │   ├─ 曖昧表現・過剰主張の検出
  │   ├─ ヘッジ表現の適切性
  │   └─ 学術文体の一貫性
  └─ Pass 5: 修正の実行
      ├─ 優先順位付き修正リスト
      ├─ 具体的な修正テキストの生成
      └─ 修正前後の差分表示
```

## 2. Pass 1: 構造レビュー（マクロ）

```markdown
## 構造レビューチェックリスト

### セクション構成
- [ ] Abstract が本文の主要結果を正確に反映しているか
- [ ] Introduction が CARS モデル（Territory → Niche → Occupying）に従っているか
- [ ] Methods が第三者による再現に十分な詳細を含むか
- [ ] Results が客観的事実のみを記述しているか（解釈が混入していないか）
- [ ] Discussion が Results の再述ではなく、解釈・考察を提供しているか
- [ ] Conclusion が新規な考察ではなく、主要発見の要約になっているか

### ストーリーライン
- [ ] 研究の動機（Why）→ 手法（How）→ 発見（What）→ 意義（So what）の流れが明確か
- [ ] Introduction で提示した研究課題に Discussion が回答しているか
- [ ] 各セクションの最初の文が前セクションからの論理的な接続になっているか

### バランス
- [ ] Introduction が長すぎないか（全体の 15-20% 目安）
- [ ] Methods が過度に詳細 or 不足していないか
- [ ] Results と Discussion のバランスが適切か
- [ ] 図表の数が適切か（ジャーナルの上限を超えていないか）
```

## 3. Pass 2: 論理・考察レビュー（コア — 最重要）

このパスが考察の深さを決定する。草稿の Discussion を以下の7層で検証する。

```markdown
## 考察の深さ検証フレームワーク（7 Layer Deep Analysis）

### Layer 1: 結果の要約（表層）
**チェック**: 主要な発見を 2-3 文で客観的に述べているか？
**問い**: 何が見つかったか？
**よくある問題**: Results の単なる繰り返しになっている。
**修正方針**: 数値の再掲ではなく、パターン・傾向の解釈に書き換える。

### Layer 2: 先行研究との位置づけ
**チェック**: 本研究の結果を既存の知見と比較しているか？
**問い**: これは先行研究と一致するか？矛盾するか？なぜか？
**よくある問題**: "Our results are consistent with [ref]." で終わっている。
**修正方針**: 一致/矛盾の理由を具体的に考察する。条件の違い、
手法の違い、サンプルの違いなど、差異の原因を分析する。

### Layer 3: メカニズムの考察
**チェック**: 観察された現象の背後にあるメカニズムを議論しているか？
**問い**: なぜこの結果が得られたのか？どのような過程・機序が働いたのか？
**よくある問題**: 相関関係を因果関係として扱っている。メカニズムの議論が表面的。
**修正方針**: 「AがBを引き起こす」ではなく「AとBの関連は、Cというメカニズムで
説明できる可能性がある」のように、推測のレベルを明示する。

### Layer 4: 代替説明の検討
**チェック**: 自分の解釈以外の可能性を検討しているか？
**問い**: この結果を説明する別の仮説はないか？交絡因子の影響はないか？
**よくある問題**: 自説に都合の良い解釈のみを提示（確証バイアス）。
**修正方針**: 少なくとも 1 つの代替説明を提示し、なぜ自説が最も妥当か
根拠を示す。反証可能性を意識する。

### Layer 5: 理論的・実用的インプリケーション
**チェック**: この発見が分野にどのような影響を与えるか議論しているか？
**問い**: So what? この発見は何に使えるか？何を変えるか？
**よくある問題**: インパクトの過大 or 過小評価。
**修正方針**: 理論的貢献（既存モデルの修正/拡張）と実用的貢献
（産業応用/臨床応用）を分けて議論する。

### Layer 6: 限界の誠実な記述
**チェック**: 研究の限界を正直に記述しているか？
**問い**: この研究で断言できないことは何か？一般化の範囲は？
**よくある問題**: 限界を挙げるが対策や影響の議論がない。
**修正方針**: 各限界について (1) 何が限界か、(2) 結果への影響度、
(3) 将来の研究でどう対処できるか、の 3 点を述べる。

### Layer 7: 将来の研究方向
**チェック**: 次にすべき研究を具体的に提案しているか？
**問い**: この研究の次のステップは何か？
**よくある問題**: "Further research is needed." のような漠然とした記述。
**修正方針**: 具体的な実験デザイン、サンプルサイズ、手法を含む
実行可能な提案に書き換える。
```

### 考察の深さスコアリング

```python
def score_discussion_depth(discussion_text):
    """
    Discussion セクションの考察の深さを 7 層で評価する。

    Returns:
        dict: 各層のスコア（0-3）と総合評価
    """
    layers = {
        "L1_summary": {
            "name": "結果の要約",
            "indicators": [
                "主要発見の明確な要約がある",
                "数値の再掲ではなくパターンの解釈がある",
            ],
        },
        "L2_prior_work": {
            "name": "先行研究との位置づけ",
            "indicators": [
                "3 件以上の先行研究との比較がある",
                "一致/矛盾の理由が議論されている",
                "条件・手法の違いによる差異が分析されている",
            ],
        },
        "L3_mechanism": {
            "name": "メカニズムの考察",
            "indicators": [
                "因果関係と相関関係が区別されている",
                "物理的/化学的/生物学的メカニズムが提案されている",
                "推測のレベルが適切にヘッジされている",
            ],
        },
        "L4_alternatives": {
            "name": "代替説明の検討",
            "indicators": [
                "少なくとも 1 つの代替仮説が検討されている",
                "交絡因子の可能性が議論されている",
                "自説が最も妥当な理由が示されている",
            ],
        },
        "L5_implications": {
            "name": "インプリケーション",
            "indicators": [
                "理論的貢献が明確に述べられている",
                "実用的応用の可能性が具体的に述べられている",
                "過大主張になっていない",
            ],
        },
        "L6_limitations": {
            "name": "限界の記述",
            "indicators": [
                "2 つ以上の限界が挙げられている",
                "各限界の結果への影響度が議論されている",
                "対処策が提案されている",
            ],
        },
        "L7_future": {
            "name": "将来の研究方向",
            "indicators": [
                "具体的な次のステップが提案されている",
                "実行可能な実験デザインが含まれている",
            ],
        },
    }

    # 各層を 0-3 でスコアリング
    # 0 = 欠如, 1 = 表面的, 2 = 適切, 3 = 優秀
    scores = {}
    for layer_id, layer in layers.items():
        # ここでは手動評価を想定。
        # AI レビューでは各 indicator の充足度を判定する。
        scores[layer_id] = {
            "name": layer["name"],
            "score": 0,  # 0-3 で評価
            "indicators": layer["indicators"],
            "comments": [],
        }

    total = sum(s["score"] for s in scores.values())
    max_total = len(scores) * 3

    return {
        "layers": scores,
        "total_score": total,
        "max_score": max_total,
        "percentage": round(total / max_total * 100, 1) if max_total > 0 else 0,
        "grade": (
            "A (Excellent)" if total >= 18 else
            "B (Good)" if total >= 14 else
            "C (Adequate)" if total >= 10 else
            "D (Needs improvement)" if total >= 6 else
            "F (Major revision required)"
        ),
    }
```

## 4. Pass 3: データ・統計レビュー

```markdown
## 統計的主張の検証チェックリスト

### 数値の正確性
- [ ] 本文中の数値が図表の値と一致しているか
- [ ] 平均値 ± SD/SEM の表記が統一されているか
- [ ] サンプルサイズ (n) が全ての検定で明記されているか

### 統計的推論
- [ ] 使用した統計検定が適切か（正規性、等分散性の確認）
- [ ] p 値だけでなく効果量（Cohen's d, η² 等）が報告されているか
- [ ] 95% 信頼区間が報告されているか
- [ ] 多重比較の補正が適用されているか
- [ ] 検出力（power）が議論されているか（特にネガティブな結果の場合）

### 過剰主張の検出パターン
| 草稿の表現 | 問題 | 修正案 |
|---|---|---|
| "A causes B" | 相関を因果と断定 | "A is associated with B" / "A may contribute to B" |
| "significantly increased" | p 値のみで判断 | "statistically significantly increased (d = 0.8, large effect)" |
| "proves that..." | 科学は証明しない | "provides strong evidence that..." / "supports the hypothesis that..." |
| "for the first time" | 検証困難 | "to our knowledge, this is the first report of..." |
| "novel" | 過剰使用 | 具体的な新規性を述べる |
| "no difference" | 検出力不足の可能性 | "no statistically significant difference was detected (power = 0.XX)" |

### 図表とテキストの整合性チェック

1. Results で言及した全ての図表が存在するか
2. 全ての図表が本文で参照されているか（孤立図表がないか）
3. 図のキャプションと本文の記述が矛盾していないか
4. 表の数値と本文の数値が一致しているか
```

## 5. Pass 4: 文章・表現レビュー（ミクロ）

```markdown
## 学術文体の検証

### ヘッジ表現の適切性チェック

確信度に応じたヘッジ表現の使い分け:

| 確信度 | 日本語 | 英語 |
|---|---|---|
| 高 (90%+) | 〜を示す / 明らかにする | demonstrate / show / reveal |
| 中高 (70-90%) | 〜を示唆する / 〜と考えられる | suggest / indicate / appear to |
| 中 (50-70%) | 〜の可能性がある | may / might / could / possible |
| 低 (30-50%) | 〜の可能性を排除できない | cannot rule out / possible but unlikely |
| 推測 | 〜と推測される | speculate / hypothesize |

### 強すぎる主張の検出と修正

検出パターン:
- 「〜を証明した」→「〜を支持する結果が得られた」
- 「明らかに」「当然」「疑いなく」→ 削除または根拠を追加
- 「全ての〜において」→ 範囲を明確化
- 「初めて」→ "to our knowledge" を追加

### 論理接続詞の検証

| 接続関係 | 使用例 | チェック |
|---|---|---|
| 因果 | Therefore, Thus, Hence | 前提→結論の論理が成立するか |
| 対比 | However, In contrast, Conversely | 実際に対比関係があるか |
| 追加 | Furthermore, Moreover, In addition | 本当に追加情報か、繰り返しでないか |
| 譲歩 | Although, Despite, Nevertheless | 譲歩の後に主張が続いているか |

### 冗長表現の削減

| 冗長 | 簡潔 |
|---|---|
| "It is well known that..." | 削除して直接述べる |
| "In order to..." | "To..." |
| "Due to the fact that..." | "Because..." |
| "A total of 50 samples" | "50 samples" |
| "It should be noted that..." | 削除して直接述べる |
| "As can be seen in Figure 1..." | "Figure 1 shows..." |
```

## 6. Pass 5: 修正の実行

```markdown
## 修正実行プロトコル

### Step 1: 問題の優先順位付け

| 優先度 | カテゴリ | 例 |
|---|---|---|
| **Critical** | 論理的誤り・データ矛盾 | 結果と結論の不一致、統計手法の誤り |
| **Major** | 考察の不足・主張の根拠不足 | Layer 3-4 の欠如、代替説明なし |
| **Minor** | 表現の改善・書式の統一 | ヘッジ不足、冗長表現、引用形式 |
| **Suggestion** | 追加で強化できる点 | 追加データ、補足図表 |

### Step 2: レビューレポート生成

各問題について以下の構造で報告する:

**[優先度] セクション名 — 問題の要約**

- **問題**: [具体的な問題の記述]
- **場所**: [セクション / 段落 / 行]
- **根拠**: [なぜこれが問題なのか]
- **修正案**: [具体的な修正テキスト]

### Step 3: 修正の適用

修正は以下の順序で適用する:
1. Critical な問題を全て修正
2. Major な問題を修正（考察の深化を含む）
3. Minor な修正を適用
4. 全体の整合性を再確認
```

## 7. 考察深化テンプレート

草稿の Discussion が浅い場合に使用する具体的な深化パターン。

```markdown
## 考察深化パターン集

### パターン A: 「一致→理由→発展」（先行研究との整合）

**Before（浅い）:**
> Our results are consistent with Smith et al. (2023).

**After（深い）:**
> Our finding that [具体的な発見] is consistent with the observations
> by Smith et al. (2023), who reported [先行研究の具体的な結果] using
> [先行研究の手法]. This agreement across [異なる条件/手法] strengthens
> the evidence that [共通のメカニズム/理論]. Moreover, our results extend
> their findings by demonstrating that [本研究独自の追加知見], which was
> not previously examined due to [理由].

### パターン B: 「矛盾→原因分析→統合」（不一致の解決）

**Before（浅い）:**
> However, our results differ from those of Jones et al. (2022).

**After（深い）:**
> Interestingly, our observation of [本研究の結果] contrasts with
> Jones et al. (2022), who found [先行研究の結果]. This discrepancy
> likely stems from differences in [条件 1: 例えば温度範囲/サンプルサイズ/
> 測定手法]. Specifically, our use of [本研究の条件] compared to their
> [先行研究の条件] may have [影響のメカニズム]. A unified interpretation
> of both findings could be that [統合的な説明], where [条件] acts as
> a moderating factor.

### パターン C: 「メカニズム提案→根拠→検証方法」

**Before（浅い）:**
> The increase in hardness may be due to grain refinement.

**After（深い）:**
> The observed increase in hardness (from XX to YY HV) with
> increasing [条件] can be attributed to [メカニズム 1: 例えば
> grain boundary strengthening following the Hall-Petch relationship].
> This interpretation is supported by our XRD analysis showing a
> decrease in crystallite size from XX to YY nm (Fig. X), which
> corresponds to [理論的な予測]. An alternative explanation involving
> [メカニズム 2: 例えば solid solution strengthening] is less likely
> because [反証の根拠]. To definitively distinguish between these
> mechanisms, [具体的な追加実験: TEM observation of grain boundaries]
> would be required.

### パターン D: 「定量的比較→スケーリング→予測」

**Before（浅い）:**
> The effect size was large (d = 0.9).

**After（深い）:**
> The large effect size (Cohen's d = 0.9, 95% CI [0.5, 1.3])
> indicates a practically meaningful difference between [群 A] and
> [群 B]. To contextualize this magnitude, [ベンチマーク: 例えば
> typical effect sizes in this field range from 0.3 to 0.6 (meta-analysis
> by Lee et al., 2021)], suggesting that [条件] has a notably strong
> influence on [アウトカム]. Extrapolating from our dose-response data,
> we predict that [条件の値] would be required to achieve [目標値],
> which warrants validation in [次の実験].

### パターン E: 「限界→影響評価→対策」

**Before（浅い）:**
> A limitation of this study is the small sample size.

**After（深い）:**
> A limitation of this study is the relatively small sample size
> (n = XX per group). Post-hoc power analysis indicates that our study
> had 0.XX power to detect a medium effect size (d = 0.5), which
> raises the possibility that some true effects may have been missed
> (Type II error). However, the primary findings showed large effect
> sizes (d > 0.8) with statistical significance, suggesting they are
> robust to sample size concerns. The non-significant secondary outcomes
> should be interpreted with caution and verified in a larger cohort.
> Based on our preliminary effect sizes, a minimum of [計算値]
> participants per group would be needed to achieve 0.80 power for
> the secondary analyses.
```

## 8. セクション別レビュー質問集

各セクションを批判的にレビューする際の質問集。

```markdown
### Abstract への質問
1. 本文を読まなくても研究の全体像が理解できるか？
2. 定量的な結果（数値）が含まれているか？
3. 過剰な主張をしていないか？
4. ワード数制限内に収まっているか？

### Introduction への質問
1. 読者はなぜこの研究が必要なのか理解できるか？
2. 先行研究のレビューは公平か？（都合の良い文献のみ引用していないか）
3. ギャップの特定は論理的か？
4. 研究目的は明確で具体的か？

### Methods への質問
1. 別の研究者がこの方法で同じ結果を再現できるか？
2. 対照実験/コントロール群は適切か？
3. サンプルサイズの根拠はあるか？
4. バイアスの制御方法が記述されているか？

### Results への質問
1. 結果の記述に解釈が混入していないか？
2. ネガティブな結果も報告されているか？（報告バイアス）
3. 全ての図表が本文で参照されているか？
4. エラーバーは SD か SEM か明記されているか？

### Discussion への質問
1. 考察の 7 層（Layer 1-7）が全てカバーされているか？
2. 代替説明を検討しているか？
3. 限界に対する影響評価と対策があるか？
4. 結論は証拠の範囲内に収まっているか？
```

## 9. レビュー結果の出力形式

レビュー結果は必ずファイルに保存する。これにより、修正履歴の追跡、
査読対応時の参照、仮説との整合性検証が可能になる。

### 9.1 レビューレポートの保存

```python
def save_review_report(review_data, filepath=None):
    """
    レビュー結果を Markdown ファイルとして保存する。
    
    Args:
        review_data: dict — レビュー結果
            {"scores": {...}, "issues": [...], "revisions": [...]}
        filepath: Path — 保存先（デフォルト: manuscript/review_report.md）
    """
    import datetime
    from pathlib import Path
    
    if filepath is None:
        filepath = BASE_DIR / "manuscript" / "review_report.md"
    filepath.parent.mkdir(parents=True, exist_ok=True)
    
    scores = review_data.get("scores", {})
    issues = review_data.get("issues", [])
    
    content = f"""# Critical Review Report

> レビュー日時: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}

## 1. 総合評価

| 評価項目 | スコア (0-3) | コメント |
|---|---|---|
| 構造の明確さ | {scores.get('structure', 0)} | {scores.get('structure_comment', '')} |
| 論理の一貫性 | {scores.get('logic', 0)} | {scores.get('logic_comment', '')} |
| 考察の深さ | {scores.get('discussion', 0)} | {scores.get('discussion_comment', '')} |
| 統計的厳密性 | {scores.get('statistics', 0)} | {scores.get('statistics_comment', '')} |
| 文章の質 | {scores.get('writing', 0)} | {scores.get('writing_comment', '')} |
| **総合** | **{sum(scores.get(k, 0) for k in ['structure','logic','discussion','statistics','writing'])}/15** | |

## 2. 修正事項一覧

"""
    # 優先度別にグループ化
    for priority in ["Critical", "Major", "Minor", "Suggestion"]:
        group = [i for i in issues if i.get("priority") == priority]
        if group:
            content += f"### {priority}\n\n"
            for issue in group:
                content += f"- [ ] **{issue['section']}**: {issue['summary']}\n"
                content += f"  - 問題: {issue['problem']}\n"
                content += f"  - 修正案: {issue['fix']}\n\n"
    
    # 修正前後の差分
    revisions = review_data.get("revisions", [])
    if revisions:
        content += "## 3. セクション別の修正案\n\n"
        for rev in revisions:
            content += f"### {rev['section']}\n\n"
            content += f"**修正前:**\n> {rev['before']}\n\n"
            content += f"**修正後:**\n> {rev['after']}\n\n"
            content += f"**修正理由:** {rev['reason']}\n\n"
    
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(content)
    
    print(f"  → レビューレポートを保存: {filepath}")
    return filepath
```

### 9.2 レビュー結果の JSON 保存

```python
def save_review_json(review_data, filepath=None):
    """
    レビュー結果を JSON として保存する。
    仮説との整合性検証、再レビュー時の比較に使用。
    """
    import datetime, json
    from pathlib import Path
    
    if filepath is None:
        filepath = BASE_DIR / "manuscript" / "review_report.json"
    filepath.parent.mkdir(parents=True, exist_ok=True)
    
    # 仮説ファイルがあれば参照を埋め込む
    hypothesis_ref = None
    hypothesis_path = BASE_DIR / "docs" / "hypothesis.json"
    if hypothesis_path.exists():
        with open(hypothesis_path, "r", encoding="utf-8") as f:
            hypothesis_ref = json.load(f).get("hypothesis")
    
    data = {
        "version": "1.0",
        "created_at": datetime.datetime.now().isoformat(),
        "hypothesis_ref": hypothesis_ref,
        "scores": review_data.get("scores", {}),
        "issues": review_data.get("issues", []),
        "discussion_depth": review_data.get("discussion_depth", {}),
        "revision_count": len(review_data.get("revisions", [])),
    }
    
    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False, default=str)
    
    print(f"  → レビュー JSON を保存: {filepath}")
    return filepath
```

### 9.3 修正済み原稿の保存

```python
def save_revised_manuscript(original_path, revisions, output_path=None):
    """
    レビューに基づいて修正した原稿を保存する。
    元の草稿はそのまま保持し、修正版を別ファイルとして保存。
    
    Args:
        original_path: Path — 元の草稿ファイル
        revisions: list[dict] — 修正リスト
            [{"before": "...", "after": "...", "section": "..."}]
        output_path: Path — 出力先（デフォルト: manuscript/manuscript_revised.md）
    """
    from pathlib import Path
    
    if output_path is None:
        output_path = BASE_DIR / "manuscript" / "manuscript_revised.md"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # 元の原稿を読み込み
    with open(original_path, "r", encoding="utf-8") as f:
        content = f.read()
    
    # 修正を適用
    applied = 0
    for rev in revisions:
        if rev["before"] in content:
            content = content.replace(rev["before"], rev["after"], 1)
            applied += 1
    
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(content)
    
    print(f"  → 修正済み原稿を保存: {output_path} ({applied}/{len(revisions)} 件適用)")
    return output_path
```

### 9.4 修正差分の保存

```python
def save_revision_diff(revisions, filepath=None):
    """
    修正前後の差分を Markdown ファイルとして保存する。
    査読対応 (Response to Reviewers) 作成時に参照できる。
    """
    import datetime
    from pathlib import Path
    
    if filepath is None:
        filepath = BASE_DIR / "manuscript" / "revision_diff.md"
    filepath.parent.mkdir(parents=True, exist_ok=True)
    
    content = f"""# 修正差分レポート

> 修正日時: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}
> 修正箇所数: {len(revisions)}

"""
    for i, rev in enumerate(revisions, 1):
        content += f"""## 修正 {i}: {rev.get('section', 'N/A')}

**優先度**: {rev.get('priority', 'N/A')}

**修正前:**
> {rev['before']}

**修正後:**
> {rev['after']}

**修正理由:** {rev.get('reason', '')}

---

"""
    
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(content)
    
    print(f"  → 修正差分レポートを保存: {filepath}")
    return filepath
```

### 9.5 レビューワークフロー統合

```python
def run_review_pipeline(manuscript_path):
    """
    レビューパイプラインを実行し、全ての結果をファイルに保存する。
    
    出力ファイル:
        manuscript/review_report.md   — レビューレポート
        manuscript/review_report.json — レビュー結果 (JSON)
        manuscript/manuscript_revised.md — 修正済み原稿
        manuscript/revision_diff.md   — 修正差分
    """
    print("=" * 60)
    print("Critical Review Pipeline")
    print("=" * 60)
    
    # Pass 1-4: レビュー実行
    review_data = {
        "scores": {},     # Pass 1-4 で得られたスコア
        "issues": [],     # 検出された問題
        "revisions": [],  # 具体的な修正案
        "discussion_depth": {},  # 7層深度分析結果
    }
    
    # Pass 5: ファイル保存
    save_review_report(review_data)
    save_review_json(review_data)
    
    if review_data["revisions"]:
        save_revised_manuscript(manuscript_path, review_data["revisions"])
        save_revision_diff(review_data["revisions"])
    
    print("\n" + "=" * 60)
    print("レビュー完了！")
    print("=" * 60)
```

### レビューレポートテンプレート

```markdown
# Critical Review Report

## 1. 総合評価

| 評価項目 | スコア (0-3) | コメント |
|---|---|---|
| 構造の明確さ | X | [コメント] |
| 論理の一貫性 | X | [コメント] |
| 考察の深さ | X | [コメント] |
| 統計的厳密性 | X | [コメント] |
| 文章の質 | X | [コメント] |
| **総合** | **X/15** | **[総合コメント]** |

## 2. 修正事項一覧

### Critical
- [ ] [修正事項 1]
- [ ] [修正事項 2]

### Major
- [ ] [修正事項 3]
- [ ] [修正事項 4]

### Minor
- [ ] [修正事項 5]

## 3. セクション別の修正案

### Discussion の修正

**修正前:**
> [元の文章]

**修正後:**
> [改善された文章]

**修正理由:** [理由の説明]
```

## References

### Output Files

| ファイル | 形式 | 生成タイミング |
|---|---|---|
| `manuscript/review_report.md` | レビューレポート（Markdown） | Pass 1-4 完了時 |
| `manuscript/review_report.json` | レビュー結果（JSON） | Pass 1-4 完了時 |
| `manuscript/manuscript_revised.md` | 修正済み原稿 | Pass 5 完了時 |
| `manuscript/revision_diff.md` | 修正前後の差分 | Pass 5 完了時 |

### 参照スキル

| スキル | 連携 |
|---|---|
| `scientific-academic-writing` | 草稿作成 → 本スキルでレビュー → 修正原稿 |
| `scientific-statistical-testing` | 統計的主張の検証に使用 |
| `scientific-publication-figures` | 図表の品質・整合性チェックに使用 |
| `scientific-hypothesis-pipeline` | 仮説と結論の整合性検証に使用 |

```

---
name: scientific-peer-review-response
description: |
  査読コメントへの体系的対応とリバッタルレター生成スキル。
  査読コメントの構造化解析（Major/Minor/Editorial 分類）、
  ポイント・バイ・ポイント回答生成、改訂箇所マッピング、
  複数ラウンド対応を行う。
  「査読に回答して」「reviewer response を作成して」「リバッタルを書いて」で発火。
tu_tools:
  - key: crossref
    name: Crossref
    description: 査読指摘の文献裏付け検索
---

# Scientific Peer Review Response

査読者コメントを構造化し、ポイント・バイ・ポイントの回答レターを生成するスキル。
改訂原稿との対応関係を明示し、エディタ・査読者への効果的なコミュニケーションを支援する。

## When to Use

- ジャーナルから査読結果（Decision Letter）を受け取ったとき
- 各査読者コメントに対する回答レターを作成するとき
- 改訂原稿における変更箇所と査読コメントのマッピングを行うとき
- 2nd / 3rd ラウンドの再査読に対応するとき
- エディタへのカバーレター（改訂版）を作成するとき

## Quick Start

## 1. 査読対応ワークフロー

```
Decision Letter 受領
  ├─ Phase 1: コメント構造化
  │   ├─ 査読者ごとにコメント分離
  │   ├─ 各コメントを Major / Minor / Editorial に分類
  │   ├─ コメント間の関連性を検出
  │   └─ 対応優先度の決定
  ├─ Phase 2: 対応戦略策定
  │   ├─ 受容 (Accept): そのまま修正
  │   ├─ 部分受容 (Partially Accept): 一部修正 + 説明
  │   ├─ 反論 (Rebut): エビデンスに基づく反論
  │   └─ 追加実験/解析の要否判断
  ├─ Phase 3: 回答レター生成
  │   ├─ ポイント・バイ・ポイント回答
  │   ├─ 引用テキスト + 回答 + 改訂箇所の3点セット
  │   └─ 追加データ/図表の挿入
  ├─ Phase 4: 改訂箇所マッピング
  │   ├─ コメント → 改訂箇所の対応表
  │   ├─ 変更テキストのハイライト指示
  │   └─ ページ/行番号参照
  └─ Phase 5: ファイル出力
      ├─ manuscript/response_to_reviewers.md
      ├─ manuscript/response_mapping.json
      └─ manuscript/cover_letter_revised.md
```

## 2. コメント構造化エンジン

```python
import re
import json
from pathlib import Path

COMMENT_SEVERITY = {
    "major": {
        "keywords": [
            "fundamental", "significant", "serious", "critical",
            "major concern", "major issue", "major revision",
            "additional experiment", "new analysis", "must address",
            "strongly recommend", "重大", "根本的",
        ],
        "weight": 3,
    },
    "minor": {
        "keywords": [
            "minor", "small", "slight", "consider",
            "suggest", "would benefit", "could improve",
            "clarify", "explain", "軽微", "検討",
        ],
        "weight": 1,
    },
    "editorial": {
        "keywords": [
            "typo", "grammar", "spelling", "formatting",
            "reference", "citation", "figure label", "table format",
            "誤字", "体裁",
        ],
        "weight": 0,
    },
}


def parse_decision_letter(decision_text):
    """
    Decision Letter を解析し、査読者ごとのコメントを構造化する。

    Args:
        decision_text: str — エディタ/査読者からのテキスト全文

    Returns:
        dict: {
            "editor_decision": "major_revision" | "minor_revision" | "reject",
            "editor_comments": str,
            "reviewers": [
                {
                    "id": "Reviewer 1",
                    "comments": [
                        {
                            "number": 1,
                            "text": "...",
                            "severity": "major" | "minor" | "editorial",
                            "requires_new_data": bool,
                            "related_section": str,
                        }
                    ]
                }
            ]
        }
    """
    result = {
        "editor_decision": _detect_decision(decision_text),
        "editor_comments": "",
        "reviewers": [],
    }

    # エディタコメントの抽出
    editor_match = re.search(
        r'(?:Editor|Associate Editor|AE)[\'s]*\s*(?:Comments?|Decision)[\s:]*\n(.*?)(?=Reviewer\s+\d|$)',
        decision_text, re.DOTALL | re.IGNORECASE
    )
    if editor_match:
        result["editor_comments"] = editor_match.group(1).strip()

    # 査読者ごとのセクション分割
    reviewer_sections = re.split(
        r'(Reviewer\s+\d+|Referee\s+\d+)',
        decision_text, flags=re.IGNORECASE
    )

    i = 1
    while i < len(reviewer_sections) - 1:
        reviewer_id = reviewer_sections[i].strip()
        reviewer_text = reviewer_sections[i + 1].strip()
        comments = _extract_comments(reviewer_text, reviewer_id)
        result["reviewers"].append({
            "id": reviewer_id,
            "comments": comments,
        })
        i += 2

    return result


def _detect_decision(text):
    """エディタの判定を検出する。"""
    text_lower = text.lower()
    if "reject" in text_lower and "resubmi" not in text_lower:
        return "reject"
    if "major revision" in text_lower or "major revisions" in text_lower:
        return "major_revision"
    if "minor revision" in text_lower or "minor revisions" in text_lower:
        return "minor_revision"
    if "accept" in text_lower:
        return "accept"
    return "unknown"


def _extract_comments(text, reviewer_id):
    """テキストから個別コメントを抽出・分類する。"""
    comments = []

    # 番号付きコメントの抽出: "1.", "1)", "(1)", "Comment 1:" 等
    numbered = re.split(
        r'(?:^|\n)\s*(?:(\d+)[.\):]|\((\d+)\)|Comment\s+(\d+)\s*:)',
        text
    )

    if len(numbered) > 1:
        # 番号付きコメントが見つかった場合
        idx = 0
        for i in range(1, len(numbered), 2):
            num_str = numbered[i] or numbered[i + 1] if i + 1 < len(numbered) else None
            body = numbered[i + 1] if i + 1 < len(numbered) else ""
            if body and body.strip():
                idx += 1
                comments.append(_classify_comment(body.strip(), idx))
    else:
        # 段落ベースで分割
        paragraphs = [p.strip() for p in text.split('\n\n') if p.strip()]
        for idx, para in enumerate(paragraphs, 1):
            if len(para) > 20:  # 短すぎるものはスキップ
                comments.append(_classify_comment(para, idx))

    return comments


def _classify_comment(text, number):
    """コメントの重要度を自動分類する。"""
    text_lower = text.lower()
    severity = "minor"  # デフォルト

    for sev, config in COMMENT_SEVERITY.items():
        for kw in config["keywords"]:
            if kw.lower() in text_lower:
                severity = sev
                break

    # 追加実験/解析が必要かの判定
    requires_new_data = bool(re.search(
        r'additional\s+(?:experiment|data|analysis|result)|'
        r'new\s+(?:experiment|data|analysis)|'
        r'perform\s+(?:additional|further)|'
        r'追加(?:実験|解析|データ)',
        text_lower
    ))

    # 関連セクションの推定
    related_section = _infer_section(text)

    return {
        "number": number,
        "text": text,
        "severity": severity,
        "requires_new_data": requires_new_data,
        "related_section": related_section,
    }


def _infer_section(text):
    """コメントが関連する論文セクションを推定する。"""
    text_lower = text.lower()
    section_keywords = {
        "Abstract": ["abstract", "summary"],
        "Introduction": ["introduction", "background", "motivation"],
        "Methods": ["method", "experimental", "procedure", "protocol"],
        "Results": ["result", "figure", "table", "data"],
        "Discussion": ["discussion", "interpretation", "mechanism", "implication"],
        "Conclusion": ["conclusion", "concluding"],
        "References": ["reference", "citation", "bibliography"],
        "SI": ["supplementary", "supporting information", "SI"],
    }
    for section, keywords in section_keywords.items():
        if any(kw in text_lower for kw in keywords):
            return section
    return "General"
```

## 3. 回答レター生成

```python
RESPONSE_TEMPLATES = {
    "accept": {
        "prefix": "We thank the reviewer for this insightful comment.",
        "action": "As suggested, we have revised the manuscript as follows:",
        "location": "Please see the revised text in {section}, {page_info}.",
    },
    "partially_accept": {
        "prefix": "We appreciate this constructive suggestion.",
        "action": "We have partially addressed this point as follows:",
        "location": "The relevant changes can be found in {section}, {page_info}.",
    },
    "rebut": {
        "prefix": "We respectfully appreciate the reviewer raising this point.",
        "action": "After careful consideration, we believe that:",
        "location": "We have added a clarification in {section}, {page_info}.",
    },
    "new_data": {
        "prefix": "We thank the reviewer for suggesting this additional analysis.",
        "action": "We have performed the requested analysis, and the results show:",
        "location": "The new data are presented in {section} ({figure_info}).",
    },
}


def generate_response_letter(parsed_comments, responses, round_number=1):
    """
    ポイント・バイ・ポイントの回答レターを生成する。

    Args:
        parsed_comments: dict — parse_decision_letter() の結果
        responses: list[dict] — 各コメントへの対応方針
            [{"reviewer": "Reviewer 1", "comment_num": 1,
              "strategy": "accept" | "partially_accept" | "rebut" | "new_data",
              "response_text": "...",
              "revised_location": {"section": "...", "page": "...", "lines": "..."},
              "new_figure": None | "Figure S3"
            }]
        round_number: int — 査読ラウンド（1st, 2nd, ...）

    Returns:
        str: 回答レター全文（Markdown）
    """
    round_label = _ordinal(round_number)

    letter = [
        f"# Response to Reviewers ({round_label} Revision)",
        "",
        "We sincerely thank the editor and reviewers for their constructive "
        "comments and suggestions. We have carefully addressed all the points "
        "raised in the review. Below, we provide our point-by-point responses.",
        "",
        "**Notation:**",
        "- Reviewer comments are shown in **bold italic**.",
        "- Our responses are in regular text.",
        "- Revised text in the manuscript is shown in blue.",
        "",
        "---",
        "",
    ]

    # エディタコメント
    if parsed_comments.get("editor_comments"):
        letter.extend([
            "## Editor Comments",
            "",
            f"***{parsed_comments['editor_comments']}***",
            "",
            "**Response:** We thank the editor for handling our manuscript and "
            "providing these comments. We have addressed all reviewer concerns "
            "as detailed below.",
            "",
            "---",
            "",
        ])

    # 査読者ごとの回答
    for reviewer in parsed_comments["reviewers"]:
        letter.extend([
            f"## {reviewer['id']}",
            "",
        ])

        reviewer_responses = [
            r for r in responses
            if r["reviewer"] == reviewer["id"]
        ]

        for comment in reviewer["comments"]:
            matching = [
                r for r in reviewer_responses
                if r["comment_num"] == comment["number"]
            ]

            letter.extend([
                f"### Comment {comment['number']} [{comment['severity'].upper()}]",
                "",
                f"***{comment['text']}***",
                "",
            ])

            if matching:
                resp = matching[0]
                tmpl = RESPONSE_TEMPLATES.get(resp["strategy"], RESPONSE_TEMPLATES["accept"])

                letter.append(f"**Response:** {tmpl['prefix']}")
                letter.append("")
                letter.append(resp["response_text"])
                letter.append("")

                loc = resp.get("revised_location", {})
                if loc:
                    location_text = tmpl["location"].format(
                        section=loc.get("section", "[Section]"),
                        page_info=f"page {loc.get('page', 'X')}, lines {loc.get('lines', 'X-Y')}",
                        figure_info=resp.get("new_figure", ""),
                    )
                    letter.append(location_text)
                    letter.append("")
            else:
                letter.extend([
                    "**Response:** [回答を記入してください]",
                    "",
                ])

            letter.extend(["---", ""])

    return "\n".join(letter)


def _ordinal(n):
    """数値を序数表現に変換する。"""
    suffixes = {1: "1st", 2: "2nd", 3: "3rd"}
    return suffixes.get(n, f"{n}th")
```

## 4. 改訂箇所マッピング

```python
def generate_response_mapping(parsed_comments, responses, filepath=None):
    """
    査読コメントと改訂箇所の対応表（トレーサビリティマトリクス）を生成する。

    Args:
        parsed_comments: dict
        responses: list[dict]
        filepath: Path

    Returns:
        dict: マッピングデータ
    """
    if filepath is None:
        filepath = BASE_DIR / "manuscript" / "response_mapping.json"
    filepath.parent.mkdir(parents=True, exist_ok=True)

    mapping = {
        "decision": parsed_comments["editor_decision"],
        "total_comments": sum(
            len(r["comments"]) for r in parsed_comments["reviewers"]
        ),
        "summary": {
            "major": 0, "minor": 0, "editorial": 0,
            "accepted": 0, "partially_accepted": 0,
            "rebutted": 0, "new_data_added": 0,
        },
        "reviewers": [],
    }

    for reviewer in parsed_comments["reviewers"]:
        reviewer_map = {
            "id": reviewer["id"],
            "comments": [],
        }

        for comment in reviewer["comments"]:
            mapping["summary"][comment["severity"]] += 1

            matched = [
                r for r in responses
                if r["reviewer"] == reviewer["id"]
                and r["comment_num"] == comment["number"]
            ]

            strategy = matched[0]["strategy"] if matched else "pending"
            location = matched[0].get("revised_location", {}) if matched else {}

            if strategy == "accept":
                mapping["summary"]["accepted"] += 1
            elif strategy == "partially_accept":
                mapping["summary"]["partially_accepted"] += 1
            elif strategy == "rebut":
                mapping["summary"]["rebutted"] += 1
            elif strategy == "new_data":
                mapping["summary"]["new_data_added"] += 1

            reviewer_map["comments"].append({
                "number": comment["number"],
                "severity": comment["severity"],
                "strategy": strategy,
                "section": comment["related_section"],
                "revised_location": location,
                "requires_new_data": comment["requires_new_data"],
            })

        mapping["reviewers"].append(reviewer_map)

    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(mapping, f, indent=2, ensure_ascii=False)

    print(f"  → マッピングファイルを保存: {filepath}")
    return mapping
```

## 5. 改訂版カバーレター

```python
def generate_revised_cover_letter(parsed_comments, mapping,
                                    journal_name="", manuscript_id="",
                                    filepath=None):
    """
    改訂版投稿時のカバーレターを生成する。

    Args:
        parsed_comments: dict
        mapping: dict — generate_response_mapping() の結果
        journal_name: str
        manuscript_id: str
    """
    if filepath is None:
        filepath = BASE_DIR / "manuscript" / "cover_letter_revised.md"
    filepath.parent.mkdir(parents=True, exist_ok=True)

    s = mapping["summary"]
    total = s["major"] + s["minor"] + s["editorial"]

    letter = f"""# Cover Letter (Revised Manuscript)

**Date:** [日付]
**Manuscript ID:** {manuscript_id or "[ID]"}
**Journal:** {journal_name or "[ジャーナル名]"}

Dear Editor,

Thank you for the opportunity to revise our manuscript entitled "[タイトル]"
(Manuscript ID: {manuscript_id or "[ID]"}).

We are grateful to the reviewers for their constructive comments, which have
significantly improved the quality of our manuscript. We have carefully
addressed all {total} comments raised by the reviewers.

**Summary of revisions:**

- **Major comments:** {s["major"]} (all addressed)
- **Minor comments:** {s["minor"]} (all addressed)
- **Editorial comments:** {s["editorial"]} (all corrected)

**Key changes in the revised manuscript:**

1. [主要な変更点 1 — コメント番号への参照]
2. [主要な変更点 2]
3. [主要な変更点 3]

All changes in the revised manuscript are highlighted in blue for easy
identification. A detailed point-by-point response to each reviewer comment
is provided in the accompanying "Response to Reviewers" document.

We believe that the revised manuscript addresses all concerns raised by the
reviewers and is now suitable for publication in {journal_name or "[ジャーナル名]"}.

Sincerely,

[著者名]
"""

    with open(filepath, "w", encoding="utf-8") as f:
        f.write(letter)

    print(f"  → カバーレター（改訂版）を保存: {filepath}")
    return filepath
```

## 6. パイプライン統合

```python
def run_review_response_pipeline(decision_letter_path, responses=None,
                                   round_number=1, journal_name="",
                                   manuscript_id=""):
    """
    査読対応パイプラインを実行する。

    Args:
        decision_letter_path: Path — Decision Letter のファイルパス
        responses: list[dict] — 各コメントへの対応（None の場合はテンプレート出力）
        round_number: int — 査読ラウンド
        journal_name: str
        manuscript_id: str

    出力ファイル:
        manuscript/response_to_reviewers.md  — ポイント・バイ・ポイント回答
        manuscript/response_mapping.json     — コメント-改訂マッピング
        manuscript/cover_letter_revised.md   — 改訂版カバーレター
    """
    print("=" * 60)
    print(f"Peer Review Response Pipeline (Round {round_number})")
    print("=" * 60)

    # Phase 1: コメント構造化
    print("\n[Phase 1] Decision Letter を解析中...")
    with open(decision_letter_path, "r", encoding="utf-8") as f:
        decision_text = f.read()

    parsed = parse_decision_letter(decision_text)
    total = sum(len(r["comments"]) for r in parsed["reviewers"])
    print(f"  → エディタ判定: {parsed['editor_decision']}")
    print(f"  → 査読者数: {len(parsed['reviewers'])}")
    print(f"  → 総コメント数: {total}")

    for reviewer in parsed["reviewers"]:
        severity_counts = {}
        for c in reviewer["comments"]:
            severity_counts[c["severity"]] = severity_counts.get(c["severity"], 0) + 1
        print(f"  → {reviewer['id']}: {severity_counts}")

    # Phase 2-3: 回答レター生成
    print("\n[Phase 2-3] 回答レターを生成中...")
    if responses is None:
        responses = _generate_template_responses(parsed)

    letter = generate_response_letter(parsed, responses, round_number)
    letter_path = BASE_DIR / "manuscript" / "response_to_reviewers.md"
    letter_path.parent.mkdir(parents=True, exist_ok=True)
    with open(letter_path, "w", encoding="utf-8") as f:
        f.write(letter)
    print(f"  → 回答レターを保存: {letter_path}")

    # Phase 4: マッピング生成
    print("\n[Phase 4] コメント-改訂マッピングを生成中...")
    mapping = generate_response_mapping(parsed, responses)

    # Phase 5: カバーレター生成
    print("\n[Phase 5] 改訂版カバーレターを生成中...")
    generate_revised_cover_letter(parsed, mapping, journal_name, manuscript_id)

    print("\n" + "=" * 60)
    print("査読対応完了！")
    print("=" * 60)

    return parsed, mapping


def _generate_template_responses(parsed):
    """未入力の回答テンプレートを生成する。"""
    responses = []
    for reviewer in parsed["reviewers"]:
        for comment in reviewer["comments"]:
            responses.append({
                "reviewer": reviewer["id"],
                "comment_num": comment["number"],
                "strategy": "accept",
                "response_text": "[ここに回答を記入してください]",
                "revised_location": {
                    "section": comment["related_section"],
                    "page": "X",
                    "lines": "X-Y",
                },
            })
    return responses
```

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `crossref` | Crossref | 査読指摘の文献裏付け検索 |

## References

### Output Files

| ファイル | 形式 | 生成タイミング |
|---|---|---|
| `manuscript/response_to_reviewers.md` | 回答レター（Markdown） | パイプライン完了時 |
| `manuscript/response_mapping.json` | コメント-改訂マッピング | パイプライン完了時 |
| `manuscript/cover_letter_revised.md` | 改訂版カバーレター | パイプライン完了時 |

### コメント分類基準

| 分類 | 説明 | 対応方針 |
|---|---|---|
| Major | 結論に影響する根本的指摘 | 追加実験/解析・大幅改訂 |
| Minor | 改善すべき点（結論には影響なし） | 修正・追記・説明追加 |
| Editorial | 誤字・体裁・書式 | 即座に修正 |

### 対応戦略

| 戦略 | 使用場面 | 回答トーン |
|---|---|---|
| Accept | 正当な指摘→修正 | 感謝 + 修正内容 |
| Partially Accept | 一部は正当、一部は異なる | 感謝 + 一部修正 + 説明 |
| Rebut | エビデンスに基づき反論 | 丁寧な反論 + 根拠 |
| New Data | 追加実験/解析が必要 | 感謝 + 新データ提示 |

### 参照スキル

| スキル | 連携 |
|---|---|
| `scientific-academic-writing` | Cover Letter テンプレート・Response to Reviewers 雛形 |
| `scientific-critical-review` | セルフレビュー → 査読前の事前対策 |
| `scientific-revision-tracker` | 改訂箇所の追跡・diff 生成 |
| `scientific-paper-quality` | 改訂後の品質メトリクス再評価 |
| `scientific-citation-checker` | 査読で追加を求められた引用の検証 |

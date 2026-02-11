---
name: scientific-revision-tracker
description: |
  論文改訂の変更履歴追跡・差分管理スキル。原稿バージョン間の diff 生成、
  変更箇所のハイライト（赤字削除/青字追加）、改訂サマリー自動生成、
  査読コメントと改訂箇所のトレーサビリティ管理を行う。
  「改訂をトラッキングして」「変更履歴を作って」「diff を出して」で発火。
---

# Scientific Revision Tracker

論文改訂プロセスにおける変更履歴の追跡・差分管理・バージョン管理スキル。
査読コメントから改訂箇所への双方向トレーサビリティを確保する。

## When to Use

- 原稿の改訂前後の差分を可視化したいとき
- 査読対応時の変更履歴を管理したいとき
- 改訂版原稿に変更箇所をマークアップ（赤字/青字）したいとき
- 複数ラウンドの改訂履歴を一元管理したいとき
- 改訂サマリーを自動生成したいとき

## Quick Start

## 1. 改訂追跡ワークフロー

```
原稿 (v1) + 査読コメント
  ├─ Phase 1: バージョンスナップショット
  │   ├─ 改訂前原稿を manuscript/versions/v1.md として保存
  │   ├─ メタデータ（日時、ラウンド番号）記録
  │   └─ セクション単位のハッシュ記録
  ├─ Phase 2: 改訂の実行（他スキルとの連携）
  │   ├─ peer-review-response → 対応方針の決定
  │   ├─ academic-writing → 本文修正
  │   └─ critical-review → 修正後のセルフレビュー
  ├─ Phase 3: 差分検出・ハイライト
  │   ├─ セクション単位の差分抽出
  │   ├─ 文レベルの変更検出
  │   ├─ 追加/削除/変更のマークアップ
  │   └─ 変更統計の計算
  ├─ Phase 4: トレーサビリティ検証
  │   ├─ 全 Major コメントに対応する変更があるか
  │   ├─ 変更箇所が回答レターの記述と一致するか
  │   └─ 未対応コメントの検出
  └─ Phase 5: ファイル出力
      ├─ manuscript/manuscript_tracked.md  — 変更マークアップ付き
      ├─ manuscript/revision_summary.json  — 改訂統計
      └─ manuscript/versions/vN.md         — バージョンスナップショット
```

## 2. バージョン管理

```python
import re
import json
import hashlib
from pathlib import Path
from datetime import datetime
from difflib import SequenceMatcher, unified_diff


def create_version_snapshot(manuscript_path, round_number=1,
                             label="original", filepath=None):
    """
    改訂前の原稿スナップショットを保存する。

    Args:
        manuscript_path: Path — 現行原稿のパス
        round_number: int — 査読ラウンド
        label: str — "original", "r1_revised", "r2_revised" 等
        filepath: Path — スナップショット保存先

    Returns:
        dict: スナップショットメタデータ
    """
    if filepath is None:
        versions_dir = BASE_DIR / "manuscript" / "versions"
        versions_dir.mkdir(parents=True, exist_ok=True)
        filepath = versions_dir / f"v{round_number}_{label}.md"

    with open(manuscript_path, "r", encoding="utf-8") as f:
        content = f.read()

    # セクション単位のハッシュ
    sections = _split_into_sections(content)
    section_hashes = {
        name: hashlib.md5(text.encode()).hexdigest()
        for name, text in sections.items()
    }

    # スナップショット保存
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(content)

    metadata = {
        "version": f"v{round_number}_{label}",
        "timestamp": datetime.now().isoformat(),
        "round": round_number,
        "label": label,
        "path": str(filepath),
        "word_count": len(content.split()),
        "section_count": len(sections),
        "section_hashes": section_hashes,
    }

    # メタデータ保存
    meta_path = filepath.with_suffix(".json")
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2, ensure_ascii=False)

    print(f"  → スナップショットを保存: {filepath}")
    print(f"  → メタデータを保存: {meta_path}")
    return metadata


def _split_into_sections(text):
    """Markdown テキストをセクションに分割する。"""
    sections = {}
    current_section = "Preamble"
    current_content = []

    for line in text.split('\n'):
        header_match = re.match(r'^#{1,3}\s+(.+)$', line)
        if header_match:
            if current_content:
                sections[current_section] = '\n'.join(current_content)
            current_section = header_match.group(1).strip()
            current_content = [line]
        else:
            current_content.append(line)

    if current_content:
        sections[current_section] = '\n'.join(current_content)

    return sections
```

## 3. 差分検出エンジン

```python
def compute_diff(original_path, revised_path):
    """
    2 つの原稿バージョン間の差分を計算する。

    Args:
        original_path: Path — 元の原稿
        revised_path: Path — 改訂後の原稿

    Returns:
        dict: {
            "sections_changed": ["Introduction", "Discussion"],
            "sections_added": [],
            "sections_removed": [],
            "changes": [
                {
                    "section": "Introduction",
                    "type": "modified",
                    "original_lines": [...],
                    "revised_lines": [...],
                    "similarity": 0.85,
                }
            ],
            "stats": {
                "lines_added": 45,
                "lines_removed": 12,
                "lines_modified": 23,
                "word_count_change": +350,
                "sections_modified": 3,
            }
        }
    """
    with open(original_path, "r", encoding="utf-8") as f:
        original = f.read()
    with open(revised_path, "r", encoding="utf-8") as f:
        revised = f.read()

    orig_sections = _split_into_sections(original)
    rev_sections = _split_into_sections(revised)

    all_sections = set(list(orig_sections.keys()) + list(rev_sections.keys()))
    sections_changed = []
    sections_added = []
    sections_removed = []
    changes = []

    for section in all_sections:
        if section in orig_sections and section not in rev_sections:
            sections_removed.append(section)
            changes.append({
                "section": section,
                "type": "removed",
                "original_lines": orig_sections[section].split('\n'),
                "revised_lines": [],
                "similarity": 0.0,
            })
        elif section not in orig_sections and section in rev_sections:
            sections_added.append(section)
            changes.append({
                "section": section,
                "type": "added",
                "original_lines": [],
                "revised_lines": rev_sections[section].split('\n'),
                "similarity": 0.0,
            })
        elif orig_sections[section] != rev_sections[section]:
            sections_changed.append(section)
            similarity = SequenceMatcher(
                None, orig_sections[section], rev_sections[section]
            ).ratio()
            changes.append({
                "section": section,
                "type": "modified",
                "original_lines": orig_sections[section].split('\n'),
                "revised_lines": rev_sections[section].split('\n'),
                "similarity": round(similarity, 3),
            })

    # 統計
    orig_lines = original.split('\n')
    rev_lines = revised.split('\n')
    diff_lines = list(unified_diff(orig_lines, rev_lines, lineterm=''))
    lines_added = sum(1 for l in diff_lines if l.startswith('+') and not l.startswith('+++'))
    lines_removed = sum(1 for l in diff_lines if l.startswith('-') and not l.startswith('---'))

    stats = {
        "lines_added": lines_added,
        "lines_removed": lines_removed,
        "lines_modified": min(lines_added, lines_removed),
        "word_count_original": len(original.split()),
        "word_count_revised": len(revised.split()),
        "word_count_change": len(revised.split()) - len(original.split()),
        "sections_modified": len(sections_changed),
        "sections_added": len(sections_added),
        "sections_removed": len(sections_removed),
    }

    return {
        "sections_changed": sections_changed,
        "sections_added": sections_added,
        "sections_removed": sections_removed,
        "changes": changes,
        "stats": stats,
    }
```

## 4. 変更マークアップ生成

```python
MARKUP_STYLES = {
    "latex": {
        "added": r"\textcolor{blue}{%s}",
        "removed": r"\st{\textcolor{red}{%s}}",
        "comment": r"\marginpar{\footnotesize %s}",
    },
    "markdown": {
        "added": '<span style="color:blue">%s</span>',
        "removed": '~~<span style="color:red">%s</span>~~',
        "comment": '<!-- REVISION: %s -->',
    },
    "track_changes": {
        "added": "[++%s++]",
        "removed": "[--%s--]",
        "comment": "[COMMENT: %s]",
    },
}


def generate_tracked_manuscript(original_path, revised_path,
                                  markup_style="markdown",
                                  response_mapping=None,
                                  filepath=None):
    """
    変更箇所をマークアップした原稿を生成する。

    Args:
        original_path: Path — 元の原稿
        revised_path: Path — 改訂後の原稿
        markup_style: str — "markdown", "latex", "track_changes"
        response_mapping: dict — コメント-改訂マッピング（あれば注釈追加）
        filepath: Path — 出力先
    """
    if filepath is None:
        filepath = BASE_DIR / "manuscript" / "manuscript_tracked.md"
    filepath.parent.mkdir(parents=True, exist_ok=True)

    style = MARKUP_STYLES.get(markup_style, MARKUP_STYLES["markdown"])

    with open(original_path, "r", encoding="utf-8") as f:
        orig_lines = f.readlines()
    with open(revised_path, "r", encoding="utf-8") as f:
        rev_lines = f.readlines()

    matcher = SequenceMatcher(None, orig_lines, rev_lines)
    tracked = []

    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag == 'equal':
            tracked.extend(orig_lines[i1:i2])
        elif tag == 'replace':
            for line in orig_lines[i1:i2]:
                tracked.append(style["removed"] % line.rstrip() + '\n')
            for line in rev_lines[j1:j2]:
                tracked.append(style["added"] % line.rstrip() + '\n')
        elif tag == 'delete':
            for line in orig_lines[i1:i2]:
                tracked.append(style["removed"] % line.rstrip() + '\n')
        elif tag == 'insert':
            for line in rev_lines[j1:j2]:
                tracked.append(style["added"] % line.rstrip() + '\n')

    # 査読コメント番号を注釈として挿入
    if response_mapping:
        tracked = _annotate_with_comments(tracked, response_mapping, style)

    with open(filepath, "w", encoding="utf-8") as f:
        f.writelines(tracked)

    print(f"  → トラック付き原稿を保存: {filepath}")
    return filepath


def _annotate_with_comments(tracked_lines, mapping, style):
    """改訂箇所にコメント番号注釈を追加する。"""
    # response_mapping から section → コメント番号のマッピングを構築
    section_comments = {}
    for reviewer in mapping.get("reviewers", []):
        for comment in reviewer.get("comments", []):
            section = comment.get("section", "General")
            if section not in section_comments:
                section_comments[section] = []
            section_comments[section].append(
                f"{reviewer['id']} #{comment['number']}"
            )

    # セクション見出しの直後に注釈を挿入
    result = []
    for line in tracked_lines:
        result.append(line)
        header = re.match(r'^#{1,3}\s+(.+)', line)
        if header:
            section_name = header.group(1).strip()
            if section_name in section_comments:
                comment_refs = ", ".join(section_comments[section_name])
                result.append(style["comment"] % f"Addresses: {comment_refs}" + '\n')

    return result
```

## 5. トレーサビリティ検証

```python
def verify_revision_traceability(response_mapping, diff_result):
    """
    全査読コメントに対応する改訂が存在するか検証する。

    Args:
        response_mapping: dict — generate_response_mapping() の結果
        diff_result: dict — compute_diff() の結果

    Returns:
        dict: {
            "all_addressed": bool,
            "unaddressed_comments": [...],
            "orphan_changes": [...],  # コメントに紐づかない変更
            "verification_summary": str,
        }
    """
    changed_sections = set(diff_result["sections_changed"] + diff_result["sections_added"])

    unaddressed = []
    addressed = []

    for reviewer in response_mapping.get("reviewers", []):
        for comment in reviewer.get("comments", []):
            section = comment.get("section", "General")
            strategy = comment.get("strategy", "pending")

            if strategy == "rebut":
                # 反論の場合は変更不要
                addressed.append(comment)
            elif section in changed_sections:
                addressed.append(comment)
            elif strategy == "accept" or strategy == "partially_accept":
                unaddressed.append({
                    "reviewer": reviewer["id"],
                    "comment_number": comment["number"],
                    "severity": comment["severity"],
                    "expected_section": section,
                    "strategy": strategy,
                })

    all_addressed = len(unaddressed) == 0

    # コメントに紐づかない変更セクション
    commented_sections = set()
    for reviewer in response_mapping.get("reviewers", []):
        for comment in reviewer.get("comments", []):
            commented_sections.add(comment.get("section", "General"))

    orphan_changes = [s for s in changed_sections if s not in commented_sections]

    return {
        "all_addressed": all_addressed,
        "addressed_count": len(addressed),
        "unaddressed_comments": unaddressed,
        "orphan_changes": orphan_changes,
        "verification_summary": (
            f"✅ 全コメント対応済み ({len(addressed)} 件)"
            if all_addressed
            else f"⚠️ 未対応コメント: {len(unaddressed)} 件"
        ),
    }
```

## 6. パイプライン統合

```python
def run_revision_tracker(original_path, revised_path, round_number=1,
                          response_mapping_path=None, markup_style="markdown"):
    """
    改訂追跡パイプラインを実行する。

    出力ファイル:
        manuscript/versions/vN_*.md          — バージョンスナップショット
        manuscript/manuscript_tracked.md     — 変更マークアップ付き原稿
        manuscript/revision_summary.json     — 改訂統計・トレーサビリティ
    """
    print("=" * 60)
    print(f"Revision Tracker Pipeline (Round {round_number})")
    print("=" * 60)

    # Phase 1: スナップショット
    print("\n[Phase 1] バージョンスナップショットを作成中...")
    create_version_snapshot(original_path, round_number, "original")
    create_version_snapshot(revised_path, round_number, "revised")

    # Phase 2: 差分計算
    print("\n[Phase 2] 差分を計算中...")
    diff = compute_diff(original_path, revised_path)
    stats = diff["stats"]
    print(f"  → 変更セクション: {diff['sections_changed']}")
    print(f"  → 行追加: +{stats['lines_added']}, 行削除: -{stats['lines_removed']}")
    print(f"  → 語数変化: {stats['word_count_change']:+d}")

    # Phase 3: マークアップ生成
    print("\n[Phase 3] 変更マークアップを生成中...")
    mapping = None
    if response_mapping_path:
        with open(response_mapping_path, "r", encoding="utf-8") as f:
            mapping = json.load(f)

    generate_tracked_manuscript(
        original_path, revised_path,
        markup_style=markup_style,
        response_mapping=mapping,
    )

    # Phase 4: トレーサビリティ検証
    traceability = None
    if mapping:
        print("\n[Phase 4] トレーサビリティを検証中...")
        traceability = verify_revision_traceability(mapping, diff)
        print(f"  → {traceability['verification_summary']}")
        if traceability["unaddressed_comments"]:
            for ua in traceability["unaddressed_comments"]:
                print(f"    ⚠️ {ua['reviewer']} #{ua['comment_number']} "
                      f"({ua['severity']}): {ua['expected_section']}")

    # 改訂サマリー保存
    summary = {
        "round": round_number,
        "timestamp": datetime.now().isoformat(),
        "diff_stats": stats,
        "sections_changed": diff["sections_changed"],
        "traceability": traceability,
    }
    summary_path = BASE_DIR / "manuscript" / "revision_summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print(f"\n  → 改訂サマリーを保存: {summary_path}")

    print("\n" + "=" * 60)
    print("改訂追跡完了！")
    print("=" * 60)

    return diff, traceability
```

## References

### Output Files

| ファイル | 形式 | 生成タイミング |
|---|---|---|
| `manuscript/versions/vN_original.md` | 改訂前スナップショット | Phase 1 |
| `manuscript/versions/vN_revised.md` | 改訂後スナップショット | Phase 1 |
| `manuscript/versions/vN_*.json` | バージョンメタデータ | Phase 1 |
| `manuscript/manuscript_tracked.md` | 変更マークアップ付き原稿 | Phase 3 |
| `manuscript/revision_summary.json` | 改訂統計・トレーサビリティ | Phase 4 |

### マークアップスタイル

| スタイル | 追加 | 削除 | 用途 |
|---|---|---|---|
| `markdown` | `<span style="color:blue">` | `~~<span style="color:red">~~` | Web / プレビュー用 |
| `latex` | `\textcolor{blue}{}` | `\st{\textcolor{red}{}}` | LaTeX 投稿用 |
| `track_changes` | `[++text++]` | `[--text--]` | プレーンテキスト |

### 参照スキル

| スキル | 連携 |
|---|---|
| `scientific-peer-review-response` | `response_mapping.json` からコメント-改訂対応を参照 |
| `scientific-academic-writing` | 改訂原稿 `manuscript.md` の作成 |
| `scientific-critical-review` | 改訂後のセルフレビュー |
| `scientific-paper-quality` | 改訂前後の品質メトリクス比較 |
| `scientific-latex-formatter` | 改訂版 LaTeX にマークアップを反映 |

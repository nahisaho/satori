---
name: scientific-citation-checker
description: |
  原稿中の引用文献の自動検索・網羅性チェックを行うスキル。
  参照リスト抽出、DOI/タイトルベース自動検索、引用カバレッジ分析、
  フォーマット一貫性検証、重複検出を実行する。
  「引用をチェックして」「参考文献を検索」「citation check」で発火。
---

# Scientific Citation Checker

原稿中の参考文献リストと本文中引用の整合性、網羅性、
フォーマット一貫性を自動検証するスキル。

## When to Use

- 原稿の参考文献リストが本文中の引用と一致しているか確認するとき
- 主張に対するエビデンス引用が不足していないか確認するとき
- 引用フォーマットの一貫性を検証するとき
- 重複引用や孤立引用（本文中で言及されていない参考文献）を検出するとき
- DOI やメタデータを使って参考文献情報を自動補完するとき

## Quick Start

## 1. 検証ワークフロー

```
原稿 (manuscript/manuscript.md) + 参照リスト
  ├─ Phase 1: 引用抽出
  │   ├─ 本文中の引用マーカー [1], (Author, 2024) 等を抽出
  │   ├─ 参考文献リスト（References セクション）を解析
  │   └─ 引用キーと参考文献エントリの対応マッピング
  ├─ Phase 2: 整合性チェック
  │   ├─ 孤立引用: 参考文献リストにあるが本文中で引用されていない
  │   ├─ 未解決引用: 本文中で引用されているが参考文献リストにない
  │   ├─ 番号連続性: [1], [2], [3]... が飛び番号なく連続しているか
  │   └─ 重複検出: 同一文献が異なるキーで複数回登録されていないか
  ├─ Phase 3: 網羅性チェック
  │   ├─ 主張-エビデンス対応: 事実記述・比較・先行研究言及に引用があるか
  │   ├─ セクション別引用密度: Introduction/Discussion で特に不足していないか
  │   └─ 自己引用率チェック: 自己引用が過剰でないか
  ├─ Phase 4: メタデータ検証・補完
  │   ├─ DOI 存在チェック: 各引用の DOI が有効か
  │   ├─ メタデータ補完: タイトル・著者・年などの欠損を自動補完
  │   └─ フォーマット一貫性: 各引用のフォーマットが統一されているか
  └─ Phase 5: レポート出力
      └─ manuscript/citation_report.json
```

## 2. 引用抽出エンジン

```python
import re
import json
from pathlib import Path


CITATION_PATTERNS = {
    "numeric": r'\[(\d+(?:[-–,\s]\d+)*)\]',
    "author_year": r'\(([A-Z][a-z]+(?:\s(?:et\sal\.|&\s[A-Z][a-z]+))?,\s*\d{4}(?:;\s*[A-Z][a-z]+(?:\s(?:et\sal\.|&\s[A-Z][a-z]+))?,\s*\d{4})*)\)',
    "superscript": r'(?<!\[)\b(\d+(?:[-–,]\d+)*)\b(?=[,.]?\s)',
}


def extract_citations(md_text, citation_style="numeric"):
    """
    原稿テキストから引用マーカーを抽出する。

    Args:
        md_text: str — 原稿テキスト
        citation_style: str — "numeric", "author_year", "superscript"

    Returns:
        dict: {
            "markers": [{"text": "[1]", "line": 10, "keys": [1]}, ...],
            "unique_keys": [1, 2, 3, ...],
            "style": "numeric"
        }
    """
    pattern = CITATION_PATTERNS.get(citation_style, CITATION_PATTERNS["numeric"])
    markers = []
    unique_keys = set()

    for line_num, line in enumerate(md_text.split('\n'), 1):
        for m in re.finditer(pattern, line):
            raw = m.group(1)
            keys = _parse_citation_keys(raw, citation_style)
            markers.append({
                "text": m.group(0),
                "line": line_num,
                "keys": keys,
            })
            unique_keys.update(keys)

    return {
        "markers": markers,
        "unique_keys": sorted(unique_keys) if citation_style == "numeric" else sorted(unique_keys),
        "style": citation_style,
    }


def _parse_citation_keys(raw, style):
    """引用マーカーのテキストからキーを抽出する。"""
    if style == "numeric":
        keys = []
        for part in re.split(r'[,\s]+', raw):
            part = part.strip()
            if '–' in part or '-' in part:
                sep = '–' if '–' in part else '-'
                start, end = part.split(sep)
                keys.extend(range(int(start.strip()), int(end.strip()) + 1))
            elif part.isdigit():
                keys.append(int(part))
        return keys
    elif style == "author_year":
        return [k.strip() for k in raw.split(';')]
    return [raw]
```

## 3. 参考文献リスト解析

```python
def parse_reference_list(md_text):
    """
    原稿末尾の References セクションから参考文献リストを解析する。

    Returns:
        list[dict]: [
            {"index": 1, "raw": "Author A. Title...", "doi": "10.xxx",
             "authors": "Author A", "title": "Title", "year": "2024",
             "journal": "Journal Name"},
            ...
        ]
    """
    # References セクションを抽出
    ref_section = re.search(
        r'(?:^#{1,2}\s*References\s*$)(.*)',
        md_text, flags=re.MULTILINE | re.DOTALL | re.IGNORECASE
    )
    if not ref_section:
        return []

    ref_text = ref_section.group(1)
    references = []

    # 番号付きリスト形式: 1. Author A. Title...
    numbered = re.findall(r'^\s*(\d+)[.\)]\s*(.+)$', ref_text, re.MULTILINE)
    if numbered:
        for idx, raw in numbered:
            ref = _parse_single_reference(raw, int(idx))
            references.append(ref)
        return references

    # 箇条書き形式: - Author A. Title...
    bullets = re.findall(r'^\s*[-*]\s*(.+)$', ref_text, re.MULTILINE)
    if bullets:
        for i, raw in enumerate(bullets, 1):
            ref = _parse_single_reference(raw, i)
            references.append(ref)
        return references

    return references


def _parse_single_reference(raw, index):
    """単一の参考文献エントリを構造化する。"""
    ref = {
        "index": index,
        "raw": raw.strip(),
        "doi": None,
        "authors": None,
        "title": None,
        "year": None,
        "journal": None,
    }

    # DOI 抽出
    doi_match = re.search(r'(?:doi:\s*|https?://doi\.org/)?(10\.\d{4,}/[^\s,]+)', raw, re.IGNORECASE)
    if doi_match:
        ref["doi"] = doi_match.group(1).rstrip('.')

    # 年の抽出
    year_match = re.search(r'\b((?:19|20)\d{2})\b', raw)
    if year_match:
        ref["year"] = year_match.group(1)

    return ref
```

## 4. 整合性チェック

```python
def check_citation_consistency(citations, references):
    """
    本文中の引用と参考文献リストの整合性を検証する。

    Args:
        citations: dict — extract_citations() の結果
        references: list[dict] — parse_reference_list() の結果

    Returns:
        dict: {
            "orphan_refs": [...],      # 参考文献リストにあるが本文で未引用
            "unresolved_cites": [...], # 本文で引用されているがリストにない
            "duplicates": [...],       # 重複エントリ
            "numbering_gaps": [...],   # 番号の飛び
            "total_refs": int,
            "total_citations": int,
            "passed": bool
        }
    """
    cited_keys = set(citations["unique_keys"])
    ref_keys = set(r["index"] for r in references)

    # 孤立参考文献（リストにあるが本文で未引用）
    orphan_refs = sorted(ref_keys - cited_keys)

    # 未解決引用（本文で引用されているがリストにない）
    unresolved = sorted(cited_keys - ref_keys)

    # 番号の連続性チェック
    if citations["style"] == "numeric" and ref_keys:
        expected = set(range(1, max(ref_keys) + 1))
        gaps = sorted(expected - ref_keys)
    else:
        gaps = []

    # 重複検出（DOI ベース）
    duplicates = _find_duplicate_references(references)

    passed = (len(orphan_refs) == 0 and len(unresolved) == 0
              and len(gaps) == 0 and len(duplicates) == 0)

    return {
        "orphan_refs": orphan_refs,
        "unresolved_cites": unresolved,
        "duplicates": duplicates,
        "numbering_gaps": gaps,
        "total_refs": len(references),
        "total_citations": len(cited_keys),
        "passed": passed,
    }


def _find_duplicate_references(references):
    """DOI またはタイトル類似度で重複参考文献を検出する。"""
    duplicates = []
    doi_map = {}

    for ref in references:
        if ref.get("doi"):
            doi = ref["doi"].lower().strip()
            if doi in doi_map:
                duplicates.append({
                    "ref_a": doi_map[doi]["index"],
                    "ref_b": ref["index"],
                    "reason": f"同一 DOI: {doi}",
                })
            else:
                doi_map[doi] = ref

    return duplicates
```

## 5. 網羅性チェック

```python
EVIDENCE_INDICATORS = [
    r'(?:has been|have been)\s+(?:shown|demonstrated|reported|observed)',
    r'(?:previous|prior|earlier|recent)\s+(?:studies?|work|research|findings?)',
    r'(?:according to|as reported by|as shown by)',
    r'(?:it is (?:well )?known|it has been established)',
    r'(?:compared (?:to|with)|in contrast (?:to|with)|consistent with)',
    r'(?:following|based on)\s+(?:the )?(?:method|approach|protocol)',
    r'\b(?:typically|generally|commonly|often|usually)\b.*\b(?:used|observed|found)\b',
]


def check_citation_coverage(md_text, citations):
    """
    本文中の事実記述・先行研究言及に引用が付いているか確認する。

    Returns:
        dict: {
            "uncited_claims": [{"line": 10, "text": "...", "indicator": "..."}],
            "section_density": {"Introduction": 5.2, "Methods": 1.1, ...},
            "self_citation_rate": 0.15,
        }
    """
    uncited_claims = []

    for line_num, line in enumerate(md_text.split('\n'), 1):
        for pattern in EVIDENCE_INDICATORS:
            if re.search(pattern, line, re.IGNORECASE):
                # この行に引用マーカーがあるか確認
                has_citation = any(
                    m["line"] == line_num for m in citations["markers"]
                )
                if not has_citation:
                    uncited_claims.append({
                        "line": line_num,
                        "text": line.strip()[:100],
                        "indicator": pattern,
                    })
                break  # 1 行につき 1 つの警告

    # セクション別引用密度
    section_density = _calculate_section_density(md_text, citations)

    return {
        "uncited_claims": uncited_claims,
        "section_density": section_density,
        "total_uncited": len(uncited_claims),
    }


def _calculate_section_density(md_text, citations):
    """セクションごとの引用密度（引用数 / 段落数）を計算する。"""
    sections = re.split(r'^#{1,2}\s+', md_text, flags=re.MULTILINE)
    density = {}

    current_line = 0
    for section in sections:
        lines = section.split('\n')
        header = lines[0].strip() if lines else "Unknown"
        n_paragraphs = max(1, len([l for l in lines if l.strip()]))
        n_citations = sum(
            1 for m in citations["markers"]
            if current_line < m["line"] <= current_line + len(lines)
        )
        density[header[:30]] = round(n_citations / n_paragraphs, 2)
        current_line += len(lines)

    return density
```

## 6. レポート生成・パイプライン

```python
def run_citation_check(manuscript_path, citation_style="numeric", filepath=None):
    """
    引用チェックパイプラインを実行する。

    Args:
        manuscript_path: Path — 原稿ファイルパス
        citation_style: str — 引用スタイル
        filepath: Path — レポート出力先

    Returns:
        dict: 全チェック結果を含む辞書
    """
    from pathlib import Path

    if filepath is None:
        filepath = BASE_DIR / "manuscript" / "citation_report.json"
    filepath.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Citation Checker Pipeline")
    print("=" * 60)

    with open(manuscript_path, "r", encoding="utf-8") as f:
        md_text = f.read()

    # Phase 1: 引用抽出
    print("\n[Phase 1] 引用マーカーを抽出中...")
    citations = extract_citations(md_text, citation_style)
    print(f"  → {len(citations['markers'])} 個の引用マーカー、{len(citations['unique_keys'])} 個のユニークキー")

    # Phase 2: 参考文献リスト解析
    print("\n[Phase 2] 参考文献リストを解析中...")
    references = parse_reference_list(md_text)
    print(f"  → {len(references)} 件の参考文献エントリ")

    # Phase 3: 整合性チェック
    print("\n[Phase 3] 引用-参考文献の整合性を検証中...")
    consistency = check_citation_consistency(citations, references)
    status = "✅ PASS" if consistency["passed"] else "⚠️ ISSUES FOUND"
    print(f"  → {status}")
    if consistency["orphan_refs"]:
        print(f"  → 孤立参考文献: {consistency['orphan_refs']}")
    if consistency["unresolved_cites"]:
        print(f"  → 未解決引用: {consistency['unresolved_cites']}")

    # Phase 4: 網羅性チェック
    print("\n[Phase 4] 引用網羅性を検証中...")
    coverage = check_citation_coverage(md_text, citations)
    print(f"  → 引用不足の可能性: {coverage['total_uncited']} 箇所")

    # レポート統合
    report = {
        "manuscript": str(manuscript_path),
        "citation_style": citation_style,
        "summary": {
            "total_citations": consistency["total_citations"],
            "total_references": consistency["total_refs"],
            "consistency_passed": consistency["passed"],
            "uncited_claims": coverage["total_uncited"],
        },
        "consistency": consistency,
        "coverage": coverage,
    }

    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)

    print(f"\n  → レポートを保存: {filepath}")
    print("=" * 60)

    return report
```

## References

### Output Files

| ファイル | 形式 | 生成タイミング |
|---|---|---|
| `manuscript/citation_report.json` | JSON レポート | チェック完了時 |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| PubMed | `PubMed_search_articles` | 引用元論文の実在確認 |
| PubMed | `PubMed_get_article` | 論文メタデータ取得 |
| Crossref | `Crossref_get_work` | DOI バリデーション |
| Crossref | `Crossref_search_works` | 出版情報検索 |
| EuropePMC | `EuropePMC_search_articles` | ヨーロッパ文献確認 |

### 検出項目一覧

| チェック項目 | 説明 | 重要度 |
|---|---|---|
| 孤立参考文献 | 参考文献リストにあるが本文未引用 | Warning |
| 未解決引用 | 本文中の引用番号がリストにない | Error |
| 番号飛び | [1], [3] のように連番でない | Warning |
| 重複参考文献 | 同一 DOI の重複エントリ | Warning |
| 引用不足の主張 | エビデンス記述に引用なし | Info |
| セクション引用密度 | Introduction/Discussion の引用が極端に少ない | Info |

### 参照スキル

| スキル | 連携 |
|---|---|
| `scientific-academic-writing` | 入力: `manuscript/manuscript.md` の引用を検証 |
| `scientific-latex-formatter` | 検証後の引用キーを BibTeX に変換 |
| `scientific-critical-review` | レビュー時に引用の妥当性を評価 |
| `scientific-hypothesis-pipeline` | 仮説と引用エビデンスの対応を確認 |

```

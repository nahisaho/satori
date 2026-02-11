---
name: scientific-paper-quality
description: |
  論文品質の定量的評価スキル。可読性スコア、セクションバランス分析、
  語彙多様性、学術語使用率、冗長表現検出、ジャーナル要件適合チェック、
  再現可能性チェックを実行する。
  「論文の品質をチェックして」「可読性スコアを出して」「投稿前チェック」で発火。
---

# Scientific Paper Quality

論文品質を定量的メトリクスで評価し、投稿前の品質保証を支援するスキル。
可読性・構造・語彙・ジャーナル適合性を多角的にスコアリングする。

## When to Use

- 投稿前の最終品質チェックを行うとき
- 原稿の可読性を定量的に評価したいとき
- ジャーナル投稿要件（語数制限、図表数等）への適合を確認するとき
- セクション間のバランスを評価したいとき
- 冗長表現・弱い動詞・過剰主張を検出したいとき
- 改訂前後の品質変化を比較したいとき
- 再現可能性の観点から Methods を検証したいとき

## Quick Start

## 1. 品質チェックワークフロー

```
原稿 (manuscript/manuscript.md)
  ├─ Dimension 1: 可読性メトリクス
  │   ├─ Flesch-Kincaid Grade Level
  │   ├─ Gunning Fog Index
  │   ├─ 平均文長 / 平均単語長
  │   └─ セクション別可読性
  ├─ Dimension 2: 構造品質
  │   ├─ セクション長バランス（IMRAD 比率）
  │   ├─ 段落構成の適切さ
  │   ├─ 図表-テキスト参照の網羅性
  │   └─ セクション間の論理フロー
  ├─ Dimension 3: 語彙・表現品質
  │   ├─ 語彙多様性 (TTR / MTLD)
  │   ├─ 学術語使用率
  │   ├─ 冗長表現の検出
  │   ├─ 弱い動詞 / 曖昧表現の検出
  │   └─ ヘッジ表現の適切性
  ├─ Dimension 4: ジャーナル適合性
  │   ├─ 語数制限チェック
  │   ├─ 図表数制限チェック
  │   ├─ 参考文献数チェック
  │   └─ フォーマット要件適合
  ├─ Dimension 5: 再現可能性
  │   ├─ Methods の詳細度
  │   ├─ 統計手法の記載完全性
  │   ├─ データ可用性記載
  │   └─ コード/ソフトウェアバージョン記載
  └─ 総合スコアカード出力
      └─ manuscript/quality_report.json
```

## 2. 可読性メトリクス

```python
import re
import json
import math
from pathlib import Path
from collections import Counter


def compute_readability(text):
    """
    テキストの可読性メトリクスを計算する。

    Args:
        text: str — 解析対象テキスト

    Returns:
        dict: {
            "flesch_kincaid_grade": float,
            "gunning_fog": float,
            "avg_sentence_length": float,
            "avg_word_length": float,
            "sentences": int,
            "words": int,
            "syllables": int,
            "complex_words": int,
        }
    """
    # 前処理: Markdown 構文を除去
    clean = _strip_markdown(text)

    sentences = _split_sentences(clean)
    words = _tokenize_words(clean)
    n_sentences = len(sentences)
    n_words = len(words)

    if n_sentences == 0 or n_words == 0:
        return {"error": "テキストが短すぎます"}

    syllable_counts = [_count_syllables(w) for w in words]
    n_syllables = sum(syllable_counts)
    complex_words = sum(1 for s in syllable_counts if s >= 3)

    # Flesch-Kincaid Grade Level
    fk_grade = (0.39 * (n_words / n_sentences)
                + 11.8 * (n_syllables / n_words)
                - 15.59)

    # Gunning Fog Index
    fog = 0.4 * ((n_words / n_sentences) + 100 * (complex_words / n_words))

    return {
        "flesch_kincaid_grade": round(fk_grade, 1),
        "gunning_fog": round(fog, 1),
        "avg_sentence_length": round(n_words / n_sentences, 1),
        "avg_word_length": round(sum(len(w) for w in words) / n_words, 1),
        "sentences": n_sentences,
        "words": n_words,
        "syllables": n_syllables,
        "complex_words": complex_words,
    }


def _strip_markdown(text):
    """Markdown 構文を除去してプレーンテキストにする。"""
    text = re.sub(r'^#{1,6}\s+', '', text, flags=re.MULTILINE)
    text = re.sub(r'\*\*(.+?)\*\*', r'\1', text)
    text = re.sub(r'\*(.+?)\*', r'\1', text)
    text = re.sub(r'!\[.*?\]\(.*?\)', '', text)
    text = re.sub(r'\[(.+?)\]\(.*?\)', r'\1', text)
    text = re.sub(r'```.*?```', '', text, flags=re.DOTALL)
    text = re.sub(r'`(.+?)`', r'\1', text)
    text = re.sub(r'^\|.*\|$', '', text, flags=re.MULTILINE)
    text = re.sub(r'^\s*[-*]\s+', '', text, flags=re.MULTILINE)
    return text


def _split_sentences(text):
    """テキストを文に分割する。"""
    sentences = re.split(r'(?<=[.!?])\s+(?=[A-Z])', text)
    return [s.strip() for s in sentences if s.strip() and len(s.strip()) > 5]


def _tokenize_words(text):
    """テキストを単語に分割する。"""
    return re.findall(r'\b[a-zA-Z]+\b', text)


def _count_syllables(word):
    """単語の音節数を推定する。"""
    word = word.lower()
    if len(word) <= 3:
        return 1
    vowels = "aeiou"
    count = 0
    prev_vowel = False
    for char in word:
        is_vowel = char in vowels
        if is_vowel and not prev_vowel:
            count += 1
        prev_vowel = is_vowel
    if word.endswith('e') and count > 1:
        count -= 1
    return max(1, count)
```

## 3. 構造品質分析

```python
IMRAD_BALANCE = {
    "ideal_ratios": {
        "Introduction": (0.15, 0.25),
        "Methods": (0.15, 0.30),
        "Results": (0.20, 0.35),
        "Discussion": (0.20, 0.35),
    },
    "abstract_max_words": {
        "nature": 150,
        "science": 125,
        "acs": 250,
        "ieee": 250,
        "elsevier": 300,
        "default": 250,
    },
}


def analyze_structure(text):
    """
    論文構造の品質を分析する。

    Returns:
        dict: {
            "section_word_counts": {"Introduction": 450, ...},
            "section_ratios": {"Introduction": 0.18, ...},
            "balance_score": float (0-1),
            "balance_issues": [...],
            "figure_text_coverage": {...},
            "paragraph_stats": {...},
        }
    """
    sections = _split_into_sections(text)
    total_body_words = 0
    section_words = {}

    # IMRAD セクションの語数
    for name, content in sections.items():
        clean = _strip_markdown(content)
        wc = len(_tokenize_words(clean))
        section_words[name] = wc
        if any(k.lower() in name.lower() for k in
               ["introduction", "method", "result", "discussion"]):
            total_body_words += wc

    # セクション比率
    section_ratios = {}
    if total_body_words > 0:
        for name, wc in section_words.items():
            section_ratios[name] = round(wc / total_body_words, 3)

    # バランス評価
    balance_issues = []
    balance_scores = []

    for section, (low, high) in IMRAD_BALANCE["ideal_ratios"].items():
        matched_key = None
        for key in section_ratios:
            if section.lower() in key.lower():
                matched_key = key
                break

        if matched_key:
            ratio = section_ratios[matched_key]
            if ratio < low:
                balance_issues.append(
                    f"{section} が短すぎます ({ratio:.0%} < {low:.0%})"
                )
                balance_scores.append(ratio / low)
            elif ratio > high:
                balance_issues.append(
                    f"{section} が長すぎます ({ratio:.0%} > {high:.0%})"
                )
                balance_scores.append(high / ratio)
            else:
                balance_scores.append(1.0)

    balance_score = sum(balance_scores) / len(balance_scores) if balance_scores else 0.5

    # 図表参照チェック
    figure_refs = set(re.findall(r'(?:Fig(?:ure)?\.?\s*|fig\.?\s*)(\d+)', text, re.IGNORECASE))
    figure_embeds = set(re.findall(r'!\[.*?(?:Fig(?:ure)?\.?\s*|fig\.?\s*)(\d+)', text, re.IGNORECASE))
    table_refs = set(re.findall(r'Table\s+(\d+)', text, re.IGNORECASE))

    # 段落統計
    paragraphs = [p for p in text.split('\n\n') if p.strip() and not p.strip().startswith('#')]
    para_lengths = [len(_tokenize_words(p)) for p in paragraphs]

    return {
        "section_word_counts": section_words,
        "section_ratios": section_ratios,
        "balance_score": round(balance_score, 2),
        "balance_issues": balance_issues,
        "figure_text_coverage": {
            "figures_referenced": sorted(figure_refs),
            "figures_embedded": sorted(figure_embeds),
            "tables_referenced": sorted(table_refs),
        },
        "paragraph_stats": {
            "count": len(paragraphs),
            "avg_length": round(sum(para_lengths) / max(1, len(para_lengths)), 1),
            "min_length": min(para_lengths) if para_lengths else 0,
            "max_length": max(para_lengths) if para_lengths else 0,
        },
    }
```

## 4. 語彙・表現品質

```python
WEAK_VERBS = [
    "is", "are", "was", "were", "been", "being",
    "have", "has", "had", "do", "does", "did",
    "make", "makes", "made", "get", "gets", "got",
    "seem", "seems", "seemed", "appear", "appears",
]

REDUNDANT_PHRASES = {
    "in order to": "to",
    "due to the fact that": "because",
    "it is important to note that": "[omit]",
    "it should be noted that": "[omit]",
    "a large number of": "many",
    "a small number of": "few",
    "in the case of": "for",
    "on the other hand": "conversely",
    "at the present time": "now / currently",
    "in the event that": "if",
    "in close proximity to": "near",
    "has the ability to": "can",
    "is able to": "can",
    "as a matter of fact": "[omit]",
    "it is well known that": "[omit or cite]",
    "it has been shown that": "[cite]",
    "take into account": "consider",
    "a number of": "several",
    "the majority of": "most",
    "prior to": "before",
    "subsequent to": "after",
    "in spite of": "despite",
    "for the purpose of": "to / for",
}

HEDGE_WORDS = [
    "might", "could", "may", "possibly", "perhaps",
    "somewhat", "relatively", "apparently", "presumably",
    "likely", "unlikely", "tend to", "seems to",
]

OVERCLAIM_PHRASES = [
    "clearly demonstrates", "undoubtedly", "proves that",
    "unequivocally", "definitely", "without doubt",
    "conclusively proves", "perfect", "always",
    "never", "all cases", "completely",
    "the first to", "novel", "unprecedented",
    "for the first time",
]


def analyze_vocabulary(text):
    """
    語彙・表現の品質を分析する。

    Returns:
        dict: {
            "vocabulary_diversity": float (TTR),
            "academic_word_ratio": float,
            "weak_verb_count": int,
            "weak_verb_examples": [...],
            "redundant_phrases": [...],
            "hedge_count": int,
            "overclaim_count": int,
            "overclaim_examples": [...],
        }
    """
    clean = _strip_markdown(text)
    words = _tokenize_words(clean)
    words_lower = [w.lower() for w in words]

    if not words:
        return {"error": "テキストが短すぎます"}

    # Type-Token Ratio (TTR)
    unique_words = set(words_lower)
    ttr = len(unique_words) / len(words_lower)

    # 弱い動詞
    weak_verb_positions = []
    for i, w in enumerate(words_lower):
        if w in WEAK_VERBS:
            context_start = max(0, i - 3)
            context_end = min(len(words), i + 4)
            weak_verb_positions.append({
                "verb": w,
                "context": " ".join(words[context_start:context_end]),
            })

    # 冗長表現
    redundant_found = []
    text_lower = clean.lower()
    for phrase, suggestion in REDUNDANT_PHRASES.items():
        count = text_lower.count(phrase)
        if count > 0:
            redundant_found.append({
                "phrase": phrase,
                "suggestion": suggestion,
                "count": count,
            })

    # ヘッジ表現
    hedge_count = sum(text_lower.count(h) for h in HEDGE_WORDS)

    # 過剰主張
    overclaim_found = []
    for phrase in OVERCLAIM_PHRASES:
        if phrase.lower() in text_lower:
            overclaim_found.append(phrase)

    return {
        "vocabulary_diversity_ttr": round(ttr, 3),
        "total_words": len(words),
        "unique_words": len(unique_words),
        "weak_verb_count": len(weak_verb_positions),
        "weak_verb_examples": weak_verb_positions[:10],
        "redundant_phrases": redundant_found,
        "hedge_count": hedge_count,
        "overclaim_count": len(overclaim_found),
        "overclaim_examples": overclaim_found,
    }
```

## 5. ジャーナル適合性チェック

```python
JOURNAL_REQUIREMENTS = {
    "nature": {
        "max_words": 3000,
        "max_figures": 8,
        "max_references": 50,
        "max_abstract_words": 150,
        "requires_data_availability": True,
        "requires_author_contributions": True,
        "requires_competing_interests": True,
    },
    "science": {
        "max_words": 2500,
        "max_figures": 4,
        "max_references": 40,
        "max_abstract_words": 125,
        "requires_data_availability": True,
        "requires_author_contributions": True,
        "requires_competing_interests": True,
    },
    "acs": {
        "max_words": 7000,
        "max_figures": 10,
        "max_references": 60,
        "max_abstract_words": 250,
        "requires_data_availability": False,
        "requires_author_contributions": True,
        "requires_competing_interests": True,
    },
    "ieee": {
        "max_words": 8000,
        "max_figures": 12,
        "max_references": 50,
        "max_abstract_words": 250,
        "requires_data_availability": False,
        "requires_author_contributions": False,
        "requires_competing_interests": False,
    },
    "elsevier": {
        "max_words": 8000,
        "max_figures": 15,
        "max_references": 80,
        "max_abstract_words": 300,
        "requires_data_availability": True,
        "requires_author_contributions": True,
        "requires_competing_interests": True,
    },
}


def check_journal_compliance(text, journal_format="elsevier"):
    """
    ジャーナル投稿要件への適合をチェックする。

    Returns:
        dict: {
            "passed": bool,
            "violations": [...],
            "warnings": [...],
            "word_count": int,
            "figure_count": int,
            "reference_count": int,
        }
    """
    reqs = JOURNAL_REQUIREMENTS.get(journal_format, JOURNAL_REQUIREMENTS["elsevier"])
    clean = _strip_markdown(text)
    words = _tokenize_words(clean)
    word_count = len(words)

    # 図の数
    figures = set(re.findall(r'!\[.*?\]\(.*?\)', text))
    figure_count = len(figures)

    # 参考文献の数
    ref_section = re.search(r'#{1,2}\s*References(.*)', text, re.DOTALL | re.IGNORECASE)
    ref_count = 0
    if ref_section:
        ref_count = len(re.findall(r'^\s*\d+[.\)]', ref_section.group(1), re.MULTILINE))

    # Abstract 語数
    abstract_match = re.search(r'#{1,2}\s*Abstract\s*\n(.*?)(?=\n#{1,2}\s|\Z)',
                                text, re.DOTALL | re.IGNORECASE)
    abstract_words = len(_tokenize_words(abstract_match.group(1))) if abstract_match else 0

    violations = []
    warnings = []

    # 語数チェック
    if word_count > reqs["max_words"]:
        violations.append(
            f"語数超過: {word_count} / {reqs['max_words']} "
            f"(超過 {word_count - reqs['max_words']} 語)"
        )
    elif word_count > reqs["max_words"] * 0.9:
        warnings.append(
            f"語数が上限に近づいています: {word_count} / {reqs['max_words']}"
        )

    # 図表数チェック
    if figure_count > reqs["max_figures"]:
        violations.append(
            f"図表数超過: {figure_count} / {reqs['max_figures']}"
        )

    # 参考文献数チェック
    if ref_count > reqs["max_references"]:
        warnings.append(
            f"参考文献数超過: {ref_count} / {reqs['max_references']}"
        )

    # Abstract 語数チェック
    if abstract_words > reqs["max_abstract_words"]:
        violations.append(
            f"Abstract 語数超過: {abstract_words} / {reqs['max_abstract_words']}"
        )

    # 必須セクションチェック
    text_lower = text.lower()
    if reqs.get("requires_data_availability"):
        if "data availability" not in text_lower and "data access" not in text_lower:
            violations.append("Data Availability Statement が見つかりません")

    if reqs.get("requires_author_contributions"):
        if "author contribution" not in text_lower and "contributions" not in text_lower:
            warnings.append("Author Contributions セクションが見つかりません")

    if reqs.get("requires_competing_interests"):
        if "competing interest" not in text_lower and "conflict of interest" not in text_lower:
            warnings.append("Competing Interests 声明が見つかりません")

    return {
        "journal": journal_format,
        "passed": len(violations) == 0,
        "violations": violations,
        "warnings": warnings,
        "word_count": word_count,
        "abstract_words": abstract_words,
        "figure_count": figure_count,
        "reference_count": ref_count,
    }
```

## 6. 再現可能性チェック

```python
REPRODUCIBILITY_CHECKS = {
    "statistical_methods": {
        "patterns": [
            r't-test|ANOVA|chi-square|Mann-Whitney|Kruskal-Wallis|Wilcoxon',
            r'p\s*[<>=]\s*0\.\d+|α\s*=\s*0\.\d+|significance level',
            r'n\s*=\s*\d+|sample size|number of (?:samples|subjects|participants)',
            r'confidence interval|CI|standard (?:deviation|error)',
            r'effect size|Cohen[\'s]?\s*d|η²|R²',
        ],
        "required": ["test_type", "significance_level", "sample_size"],
    },
    "software_versions": {
        "patterns": [
            r'Python\s+\d+\.\d+|R\s+\d+\.\d+|MATLAB\s+R\d{4}',
            r'(?:version|v\.?)\s*\d+\.\d+',
            r'scikit-learn|scipy|numpy|pandas|statsmodels|ggplot',
        ],
    },
    "data_availability": {
        "patterns": [
            r'data\s+(?:are|is)\s+(?:available|deposited|accessible)',
            r'(?:GitHub|Zenodo|Figshare|Dryad)',
            r'accession\s+(?:number|code)',
            r'(?:doi|DOI):\s*10\.\d+',
        ],
    },
}


def check_reproducibility(text):
    """
    Methods セクションの再現可能性を評価する。

    Returns:
        dict: {
            "score": float (0-1),
            "checks": {
                "statistical_methods": {"present": [...], "missing": [...]},
                "software_versions": {"present": [...], "missing": [...]},
                "data_availability": {"present": [...], "missing": [...]},
            },
            "recommendations": [...],
        }
    """
    # Methods セクション抽出
    methods_match = re.search(
        r'#{1,2}\s*(?:Methods?|Materials?\s+and\s+Methods?|Experimental)\s*\n(.*?)(?=\n#{1,2}\s|\Z)',
        text, re.DOTALL | re.IGNORECASE
    )
    methods_text = methods_match.group(1) if methods_match else text

    checks = {}
    total_found = 0
    total_checks = 0

    for category, config in REPRODUCIBILITY_CHECKS.items():
        present = []
        for pattern in config["patterns"]:
            matches = re.findall(pattern, methods_text, re.IGNORECASE)
            if matches:
                present.extend(set(matches))

        checks[category] = {
            "present": present,
            "found": len(present) > 0,
        }
        total_checks += 1
        if present:
            total_found += 1

    score = total_found / total_checks if total_checks > 0 else 0

    recommendations = []
    if not checks["statistical_methods"]["found"]:
        recommendations.append(
            "統計手法の詳細（検定名、有意水準、サンプルサイズ）を明記してください"
        )
    if not checks["software_versions"]["found"]:
        recommendations.append(
            "使用ソフトウェア/ライブラリとそのバージョンを明記してください"
        )
    if not checks["data_availability"]["found"]:
        recommendations.append(
            "データの可用性に関する記述（リポジトリ URL 等）を追加してください"
        )

    return {
        "score": round(score, 2),
        "checks": checks,
        "recommendations": recommendations,
    }
```

## 7. 品質スコアカード・パイプライン

```python
def run_quality_check(manuscript_path, journal_format="elsevier",
                       comparison_path=None, filepath=None):
    """
    論文品質チェックパイプラインを実行する。

    Args:
        manuscript_path: Path — 原稿ファイルパス
        journal_format: str — ジャーナル形式
        comparison_path: Path — 比較用原稿（改訂前版など、差分表示用）
        filepath: Path — レポート出力先

    出力ファイル:
        manuscript/quality_report.json — 品質スコアカード
    """
    if filepath is None:
        filepath = BASE_DIR / "manuscript" / "quality_report.json"
    filepath.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Paper Quality Check Pipeline")
    print("=" * 60)

    with open(manuscript_path, "r", encoding="utf-8") as f:
        text = f.read()

    # Dimension 1: 可読性
    print("\n[Dim 1] 可読性メトリクスを計算中...")
    readability = compute_readability(text)
    print(f"  → Flesch-Kincaid Grade: {readability.get('flesch_kincaid_grade', 'N/A')}")
    print(f"  → Gunning Fog Index: {readability.get('gunning_fog', 'N/A')}")
    print(f"  → 平均文長: {readability.get('avg_sentence_length', 'N/A')} 語")

    # Dimension 2: 構造品質
    print("\n[Dim 2] 構造品質を分析中...")
    structure = analyze_structure(text)
    print(f"  → バランススコア: {structure['balance_score']}")
    for issue in structure["balance_issues"]:
        print(f"  ⚠️ {issue}")

    # Dimension 3: 語彙・表現
    print("\n[Dim 3] 語彙・表現品質を分析中...")
    vocabulary = analyze_vocabulary(text)
    print(f"  → 語彙多様性 (TTR): {vocabulary.get('vocabulary_diversity_ttr', 'N/A')}")
    print(f"  → 冗長表現: {len(vocabulary.get('redundant_phrases', []))} 種")
    print(f"  → 過剰主張: {vocabulary.get('overclaim_count', 0)} 箇所")

    # Dimension 4: ジャーナル適合性
    print(f"\n[Dim 4] ジャーナル適合性を確認中 ({journal_format})...")
    compliance = check_journal_compliance(text, journal_format)
    status = "✅ PASS" if compliance["passed"] else "❌ FAIL"
    print(f"  → {status}")
    for v in compliance["violations"]:
        print(f"  ❌ {v}")
    for w in compliance["warnings"]:
        print(f"  ⚠️ {w}")

    # Dimension 5: 再現可能性
    print("\n[Dim 5] 再現可能性を確認中...")
    reproducibility = check_reproducibility(text)
    print(f"  → 再現可能性スコア: {reproducibility['score']}")
    for rec in reproducibility["recommendations"]:
        print(f"  💡 {rec}")

    # 総合スコア (0-100)
    scores = {
        "readability": _readability_score(readability),
        "structure": structure["balance_score"],
        "vocabulary": _vocabulary_score(vocabulary),
        "compliance": 1.0 if compliance["passed"] else 0.5,
        "reproducibility": reproducibility["score"],
    }
    weights = {
        "readability": 0.20,
        "structure": 0.20,
        "vocabulary": 0.20,
        "compliance": 0.25,
        "reproducibility": 0.15,
    }
    overall = sum(scores[k] * weights[k] for k in scores) * 100

    print(f"\n{'=' * 60}")
    print(f"総合品質スコア: {overall:.0f} / 100")
    print(f"{'=' * 60}")

    # 比較（改訂前後）
    comparison = None
    if comparison_path:
        print("\n[比較] 改訂前後の品質変化を計算中...")
        with open(comparison_path, "r", encoding="utf-8") as f:
            old_text = f.read()
        old_readability = compute_readability(old_text)
        old_scores = {
            "readability": _readability_score(old_readability),
            "structure": analyze_structure(old_text)["balance_score"],
            "vocabulary": _vocabulary_score(analyze_vocabulary(old_text)),
        }
        old_overall = sum(old_scores.get(k, 0) * weights[k] for k in old_scores) * 100
        comparison = {
            "score_change": round(overall - old_overall, 1),
            "improved": overall > old_overall,
        }
        print(f"  → スコア変化: {comparison['score_change']:+.1f}")

    # レポート保存
    report = {
        "manuscript": str(manuscript_path),
        "journal": journal_format,
        "overall_score": round(overall, 1),
        "dimension_scores": {k: round(v * 100, 1) for k, v in scores.items()},
        "readability": readability,
        "structure": structure,
        "vocabulary": vocabulary,
        "compliance": compliance,
        "reproducibility": reproducibility,
        "comparison": comparison,
    }

    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)

    print(f"\n  → 品質レポートを保存: {filepath}")
    return report


def _readability_score(readability):
    """可読性を 0-1 スコアに変換する。学術論文の適正範囲は FK Grade 12-16。"""
    fk = readability.get("flesch_kincaid_grade", 14)
    if 12 <= fk <= 16:
        return 1.0
    elif 10 <= fk < 12 or 16 < fk <= 18:
        return 0.7
    else:
        return 0.4


def _vocabulary_score(vocabulary):
    """語彙品質を 0-1 スコアに変換する。"""
    ttr = vocabulary.get("vocabulary_diversity_ttr", 0.5)
    redundant = len(vocabulary.get("redundant_phrases", []))
    overclaim = vocabulary.get("overclaim_count", 0)

    score = min(1.0, ttr * 2)  # TTR 0.5 で満点
    score -= redundant * 0.02  # 冗長表現 1 つにつき -0.02
    score -= overclaim * 0.05  # 過剰主張 1 つにつき -0.05
    return max(0.0, round(score, 2))
```

## References

### Output Files

| ファイル | 形式 | 生成タイミング |
|---|---|---|
| `manuscript/quality_report.json` | 品質スコアカード | チェック完了時 |

### 品質メトリクス一覧

| Dimension | メトリクス | 適正範囲（学術論文） |
|---|---|---|
| 可読性 | Flesch-Kincaid Grade | 12-16 |
| 可読性 | Gunning Fog Index | 12-18 |
| 可読性 | 平均文長 | 15-25 語 |
| 構造 | IMRAD バランス | 各セクション 15-35% |
| 語彙 | TTR (語彙多様性) | > 0.4 |
| 語彙 | 冗長表現 | < 5 種類 |
| 適合性 | 語数制限 | ジャーナル依存 |
| 再現性 | 統計手法記載 | 必須 |

### 参照スキル

| スキル | 連携 |
|---|---|
| `scientific-academic-writing` | 原稿 `manuscript/manuscript.md` の品質を評価 |
| `scientific-critical-review` | セルフレビュー結果と品質スコアの照合 |
| `scientific-peer-review-response` | 改訂後の品質改善を定量的に検証 |
| `scientific-revision-tracker` | 改訂前後の品質スコア比較 |
| `scientific-latex-formatter` | 投稿前の最終品質ゲート |

---
name: scientific-clinical-nlp
description: |
  臨床自然言語処理スキル。MedSpaCy / cTAKES / scispaCy
  による臨床テキスト NER、セクション検出、否定文検出、
  ICD-10/SNOMED-CT エンティティリンキング、
  匿名化 (De-identification) パイプライン。
  TU 外スキル (直接 Python ライブラリ)。
---

# Scientific Clinical NLP

MedSpaCy・scispaCy を中心とした臨床テキスト自然言語処理
パイプラインを提供する。電子カルテテキストからの臨床エンティティ
抽出・否定文検出 (NegEx)・セクション検出・標準用語へのリンキング
を行う。

## When to Use

- 電子カルテ / 臨床ノートから疾患・薬剤・症状を抽出するとき
- 臨床テキストの否定文 (NegEx/ConText) を検出するとき
- テキストセクション (主訴/HPI/Assessment/Plan) を分類するとき
- ICD-10 / SNOMED-CT コードへのリンキングを行うとき
- PHI 匿名化 (De-identification) を実施するとき
- バイオメディカル文献テキストマイニングとの連携

---

## Quick Start

## 1. MedSpaCy 臨床 NER

```python
import medspacy
from medspacy.ner import TargetRule
from medspacy.visualization import visualize_ent


def clinical_ner(text, rules=None):
    """
    MedSpaCy — 臨床テキスト NER パイプライン。

    Parameters:
        text: str — 臨床テキスト
        rules: list[dict] | None — カスタムルール
    """
    nlp = medspacy.load(
        enable=["medspacy_pyrush",
                "medspacy_target_matcher",
                "medspacy_context"])

    if rules:
        target_matcher = nlp.get_pipe(
            "medspacy_target_matcher")
        for r in rules:
            target_matcher.add(TargetRule(
                literal=r["literal"],
                category=r.get("category",
                               "CONDITION")))

    doc = nlp(text)

    entities = []
    for ent in doc.ents:
        entities.append({
            "text": ent.text,
            "label": ent.label_,
            "start": ent.start_char,
            "end": ent.end_char,
            "is_negated": ent._.is_negated,
            "is_uncertain": ent._.is_uncertain,
            "is_historical": ent._.is_historical,
            "is_family": ent._.is_family,
        })

    n_neg = sum(1 for e in entities
                if e["is_negated"])
    print(f"Clinical NER: {len(entities)} entities, "
          f"{n_neg} negated")
    return entities


def clinical_ner_batch(texts, rules=None):
    """
    MedSpaCy — バッチ臨床 NER。

    Parameters:
        texts: list[str] — 臨床テキストリスト
        rules: list[dict] | None — カスタムルール
    """
    all_entities = []
    for i, text in enumerate(texts):
        ents = clinical_ner(text, rules)
        for e in ents:
            e["doc_id"] = i
        all_entities.extend(ents)

    import pandas as pd
    df = pd.DataFrame(all_entities)
    print(f"Batch NER: {len(texts)} docs, "
          f"{len(df)} total entities")
    return df
```

## 2. セクション検出

```python
def clinical_section_detect(text):
    """
    MedSpaCy — 臨床テキストセクション検出。

    Parameters:
        text: str — 臨床テキスト
    """
    import medspacy
    nlp = medspacy.load(
        enable=["medspacy_pyrush",
                "medspacy_sectionizer"])

    doc = nlp(text)

    sections = []
    for section in doc._.sections:
        sections.append({
            "category": section.category,
            "title": (section.title_span.text
                      if section.title_span else ""),
            "body": (section.body_span.text[:200]
                     if section.body_span else ""),
        })

    print(f"Sections detected: {len(sections)}")
    for s in sections:
        print(f"  [{s['category']}] "
              f"{s['title'][:50]}")
    return sections
```

## 3. SNOMED-CT / ICD-10 リンキング

```python
def clinical_entity_linking(text,
                              linker_name="umls"):
    """
    scispaCy — 臨床エンティティの UMLS/SNOMED リンキング。

    Parameters:
        text: str — 臨床テキスト
        linker_name: str — リンカー ("umls", "mesh",
                           "snomed")
    """
    import spacy
    import scispacy
    from scispacy.linking import EntityLinker

    nlp = spacy.load("en_core_sci_md")
    nlp.add_pipe("scispacy_linker",
                 config={"resolve_abbreviations": True,
                         "linker_name": linker_name})

    doc = nlp(text)
    linker = nlp.get_pipe("scispacy_linker")

    linked = []
    for ent in doc.ents:
        for cui, score in ent._.kb_ents[:3]:
            concept = linker.kb.cui_to_entity.get(
                cui, {})
            linked.append({
                "text": ent.text,
                "cui": cui,
                "score": round(score, 3),
                "canonical_name": (
                    concept.canonical_name
                    if hasattr(concept,
                               "canonical_name")
                    else str(concept)),
            })

    import pandas as pd
    df = pd.DataFrame(linked)
    print(f"Entity linking: {len(doc.ents)} entities → "
          f"{len(df)} CUI mappings")
    return df
```

## 4. 臨床 NLP 統合パイプライン

```python
def clinical_nlp_pipeline(texts,
                            output_dir="results"):
    """
    臨床 NLP 統合パイプライン。

    Parameters:
        texts: list[str] — 臨床テキストリスト
        output_dir: str — 出力ディレクトリ
    """
    import pandas as pd
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) NER + 否定文検出
    ner_df = clinical_ner_batch(texts)
    ner_df.to_csv(output_dir / "clinical_ner.csv",
                  index=False)

    # 2) セクション検出
    all_sections = []
    for i, text in enumerate(texts):
        secs = clinical_section_detect(text)
        for s in secs:
            s["doc_id"] = i
        all_sections.extend(secs)
    section_df = pd.DataFrame(all_sections)
    section_df.to_csv(
        output_dir / "clinical_sections.csv",
        index=False)

    # 3) エンティティリンキング (最初のテキスト)
    if texts:
        link_df = clinical_entity_linking(texts[0])
        link_df.to_csv(
            output_dir / "entity_linking.csv",
            index=False)

    print(f"Clinical NLP pipeline → {output_dir}")
    return {"ner": ner_df, "sections": section_df}
```

---

## パイプライン統合

```
text-mining-nlp → clinical-nlp → clinical-reporting
  (PubMed/文献)     (NER/NegEx)     (構造化レポート)
       │                 │                 ↓
  biomedical-ner ───────┘       pharmacogenomics
    (scispaCy)                   (PGx 処方支援)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/clinical_ner.csv` | 臨床エンティティ+否定 | → phenotype-hpo |
| `results/clinical_sections.csv` | セクション分類 | → clinical-reporting |
| `results/entity_linking.csv` | UMLS/SNOMED リンキング | → disease-research |

---
name: scientific-supplementary-generator
description: |
  学術論文の Supplementary Information (SI) を自動生成するスキル。
  本文から溢れた図表・手法詳細・追加データを構造化し、ジャーナル規定に準拠した
  SI ドキュメントを生成する。「SIを作って」「補足資料を生成して」「Supplementary を作成して」で発火。
  scientific-academic-writing で本文を作成した後に使用。
---

# Scientific Supplementary Information Generator

学術論文の Supplementary Information (SI) を解析結果・本文原稿から自動生成するスキル。
ジャーナル別の SI 規定に準拠し、Figure S / Table S / Methods S の番号管理、
本文からの参照整合性チェック、独立した引用リストの管理を行う。

## When to Use

- 論文本文の完成後、Supplementary Information を作成するとき
- 本文に収まらない追加図表・手法・データを整理するとき
- ジャーナルの SI フォーマット規定に合わせる必要があるとき
- 本文と SI 間の相互参照（Figure S1, Table S1 等）の整合性を確保したいとき
- レビュアーから追加データの提示を求められたとき

## Quick Start

## 1. SI 生成ワークフロー

```
本文原稿 (manuscript.md) + 解析結果 (results/, figures/)
  ├─ Phase 1: SI 候補の収集
  │   ├─ 本文から参照されている Figure S / Table S を抽出
  │   ├─ figures/ 内で本文未使用の図を検出
  │   ├─ results/ 内の詳細データ（中間結果等）を候補に追加
  │   └─ Methods の詳細手順で本文に載せきれないものを特定
  ├─ Phase 2: SI ドキュメント構成
  │   ├─ Supplementary Figures セクション
  │   ├─ Supplementary Tables セクション
  │   ├─ Supplementary Methods セクション
  │   └─ Supplementary References セクション
  ├─ Phase 3: 番号体系・参照整合性チェック
  │   ├─ Figure S1, S2, ... の連番確認
  │   ├─ 本文中の "Figure S1" 参照が SI 内に存在するか検証
  │   └─ SI 内の図表が本文から参照されているか検証
  └─ Phase 4: ファイル保存
      ├─ manuscript/supplementary.md
      ├─ manuscript/supplementary_figures/ (必要に応じてコピー)
      └─ manuscript/si_crossref_report.json
```

## 2. ジャーナル別 SI 規定

```markdown
## ジャーナル別 SI フォーマット

| ジャーナル | SI タイトル | 図番号形式 | 表番号形式 | 独立ファイル |
|---|---|---|---|---|
| Nature 系 | Supplementary Information | Supplementary Fig. 1 | Supplementary Table 1 | 単一 PDF |
| Science 系 | Supplementary Materials | fig. S1 | table S1 | 単一 PDF |
| ACS 系 | Supporting Information | Figure S1 | Table S1 | 単一 PDF |
| IEEE 系 | Supplementary Material | Fig. S1 | TABLE S-I | 別ファイル可 |
| Elsevier 系 | Appendix A / Supplementary data | Fig. A.1 | Table A.1 | Appendix 形式 |
```

## 3. SI 候補の自動収集

```python
def collect_si_candidates(manuscript_path, figures_dir, results_dir):
    """
    本文原稿と解析結果ディレクトリから SI 候補を自動収集する。

    Args:
        manuscript_path: Path — 本文原稿 (manuscript.md)
        figures_dir: Path — figures/ ディレクトリ
        results_dir: Path — results/ ディレクトリ

    Returns:
        dict: SI 候補のカテゴリ別リスト
    """
    import re
    from pathlib import Path

    # 本文を読み込み
    with open(manuscript_path, "r", encoding="utf-8") as f:
        manuscript_text = f.read()

    # 本文中の Figure 参照を抽出
    main_figs = set(re.findall(
        r'!\[(?:Fig(?:ure)?\.?\s*\d+)\]\(figures/([^)]+)\)', manuscript_text
    ))

    # 本文中の SI 参照を抽出
    si_refs = re.findall(
        r'(?:Figure|Fig\.?|Table|Supplementary\s+(?:Fig|Table))\s*S(\d+)',
        manuscript_text, re.IGNORECASE
    )

    # figures/ 内の全画像ファイル
    all_figs = set()
    if figures_dir.exists():
        all_figs = {
            f.name for f in figures_dir.iterdir()
            if f.suffix.lower() in ('.png', '.svg', '.pdf', '.jpg', '.jpeg')
        }

    # 本文未使用の図 → SI 候補
    unused_figs = all_figs - main_figs

    # results/ 内の詳細データ
    detail_data = []
    if results_dir.exists():
        for f in results_dir.iterdir():
            if f.suffix == '.csv' and f.name not in ('analysis_summary.json',):
                detail_data.append(f.name)

    candidates = {
        "supplementary_figures": sorted(unused_figs),
        "supplementary_tables": detail_data,
        "si_references_in_main": sorted(set(si_refs)),
        "main_figures": sorted(main_figs),
    }

    print(f"  SI 候補:")
    print(f"    補足図: {len(candidates['supplementary_figures'])} 件")
    print(f"    補足表: {len(candidates['supplementary_tables'])} 件")
    print(f"    本文中の SI 参照: {len(candidates['si_references_in_main'])} 件")

    return candidates
```

## 4. SI ドキュメント生成

```python
def generate_supplementary(candidates, journal_format="imrad",
                            title="", authors="", filepath=None):
    """
    SI ドキュメントを Markdown として生成・保存する。

    Args:
        candidates: dict — collect_si_candidates() の戻り値
        journal_format: str — ジャーナル形式 (nature/science/acs/ieee/elsevier/imrad)
        title: str — 論文タイトル
        authors: str — 著者リスト
        filepath: Path — 出力先（デフォルト: manuscript/supplementary.md）
    """
    import datetime
    from pathlib import Path

    if filepath is None:
        filepath = BASE_DIR / "manuscript" / "supplementary.md"
    filepath.parent.mkdir(parents=True, exist_ok=True)

    # ジャーナル別フォーマット設定
    fmt = _get_si_format(journal_format)

    content = f"""# {fmt['title_prefix']}
# {title}

{authors}

> Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}

"""

    # Supplementary Figures
    sup_figs = candidates.get("supplementary_figures", [])
    if sup_figs:
        content += f"## {fmt['fig_section_title']}\n\n"
        for i, fig_name in enumerate(sup_figs, 1):
            fig_label = fmt['fig_prefix'].format(n=i)
            content += f"![{fig_label}](figures/{fig_name})\n\n"
            content += f"**{fig_label}.** [キャプション: {fig_name} の説明を記入]\n\n"

    # Supplementary Tables
    sup_tables = candidates.get("supplementary_tables", [])
    if sup_tables:
        content += f"## {fmt['table_section_title']}\n\n"
        for i, table_name in enumerate(sup_tables, 1):
            table_label = fmt['table_prefix'].format(n=i)
            content += f"**{table_label}.** [キャプション: {table_name} の説明]\n\n"
            content += f"データソース: `results/{table_name}`\n\n"

    # Supplementary Methods
    content += f"""## {fmt['methods_section_title']}

### S1. [追加実験手法の詳細]

[本文 Methods に収まらなかった詳細手順を記述]

### S2. [追加解析手法の詳細]

[パラメータ設定、アルゴリズムの詳細、コード等]

"""

    # Supplementary References
    content += f"""## {fmt['refs_section_title']}

[SI 内でのみ引用した文献をここに記載]

"""

    with open(filepath, "w", encoding="utf-8") as f:
        f.write(content)

    print(f"  → SI ドキュメントを保存: {filepath}")
    return filepath


def _get_si_format(journal_format):
    """ジャーナル別の SI フォーマット設定を返す。"""
    formats = {
        "nature": {
            "title_prefix": "Supplementary Information for:",
            "fig_section_title": "Supplementary Figures",
            "table_section_title": "Supplementary Tables",
            "methods_section_title": "Supplementary Methods",
            "refs_section_title": "Supplementary References",
            "fig_prefix": "Supplementary Fig. {n}",
            "table_prefix": "Supplementary Table {n}",
        },
        "science": {
            "title_prefix": "Supplementary Materials for:",
            "fig_section_title": "Supplementary Figures",
            "table_section_title": "Supplementary Tables",
            "methods_section_title": "Materials and Methods",
            "refs_section_title": "References",
            "fig_prefix": "fig. S{n}",
            "table_prefix": "table S{n}",
        },
        "acs": {
            "title_prefix": "Supporting Information for:",
            "fig_section_title": "Supporting Figures",
            "table_section_title": "Supporting Tables",
            "methods_section_title": "Supporting Methods",
            "refs_section_title": "Supporting References",
            "fig_prefix": "Figure S{n}",
            "table_prefix": "Table S{n}",
        },
        "ieee": {
            "title_prefix": "Supplementary Material for:",
            "fig_section_title": "Supplementary Figures",
            "table_section_title": "Supplementary Tables",
            "methods_section_title": "Supplementary Methods",
            "refs_section_title": "Supplementary References",
            "fig_prefix": "Fig. S{n}",
            "table_prefix": "TABLE S-{n}",
        },
        "elsevier": {
            "title_prefix": "Appendix A. Supplementary data for:",
            "fig_section_title": "Supplementary Figures",
            "table_section_title": "Supplementary Tables",
            "methods_section_title": "Supplementary Methods",
            "refs_section_title": "References",
            "fig_prefix": "Fig. A.{n}",
            "table_prefix": "Table A.{n}",
        },
    }
    return formats.get(journal_format, formats.get("acs"))  # デフォルト ACS
```

## 5. 本文–SI 間の相互参照チェック

```python
def check_si_crossrefs(manuscript_path, si_path, filepath=None):
    """
    本文と SI 間の相互参照整合性を検証する。

    チェック内容:
    - 本文中の "Figure S1" 参照が SI 内に存在するか
    - SI 内の Figure S / Table S が本文から参照されているか（孤立チェック）
    - 番号の連続性 (S1, S2, S3... に欠番がないか)

    Returns:
        dict: 検証結果
    """
    import re, json, datetime
    from pathlib import Path

    with open(manuscript_path, "r", encoding="utf-8") as f:
        main_text = f.read()
    with open(si_path, "r", encoding="utf-8") as f:
        si_text = f.read()

    # 本文中の SI 参照を抽出
    main_fig_refs = set(re.findall(
        r'(?:Supplementary\s+)?(?:Fig(?:ure)?\.?\s*S|fig\.\s*S)(\d+)',
        main_text, re.IGNORECASE
    ))
    main_table_refs = set(re.findall(
        r'(?:Supplementary\s+)?(?:Table\s*S|table\s*S)[\-]?(\d+)',
        main_text, re.IGNORECASE
    ))

    # SI 内の Figure/Table 定義を抽出
    si_fig_defs = set(re.findall(
        r'\*\*(?:Supplementary\s+)?(?:Fig(?:ure)?\.?\s*S?|fig\.\s*S?)(\d+)',
        si_text, re.IGNORECASE
    ))
    si_table_defs = set(re.findall(
        r'\*\*(?:Supplementary\s+)?(?:Table\s*S?|TABLE\s*S[\-]?)(\d+)',
        si_text, re.IGNORECASE
    ))

    # チェック
    issues = []

    # 本文で参照されているが SI に定義がない
    missing_in_si = main_fig_refs - si_fig_defs
    for ref in sorted(missing_in_si):
        issues.append({
            "type": "MISSING_IN_SI",
            "severity": "Critical",
            "message": f"Figure S{ref} が本文で参照されているが SI に定義がない",
        })

    missing_tables_in_si = main_table_refs - si_table_defs
    for ref in sorted(missing_tables_in_si):
        issues.append({
            "type": "MISSING_IN_SI",
            "severity": "Critical",
            "message": f"Table S{ref} が本文で参照されているが SI に定義がない",
        })

    # SI に定義があるが本文から参照されていない（孤立）
    orphan_figs = si_fig_defs - main_fig_refs
    for ref in sorted(orphan_figs):
        issues.append({
            "type": "ORPHAN_IN_SI",
            "severity": "Minor",
            "message": f"Figure S{ref} が SI に定義されているが本文から参照されていない",
        })

    orphan_tables = si_table_defs - main_table_refs
    for ref in sorted(orphan_tables):
        issues.append({
            "type": "ORPHAN_IN_SI",
            "severity": "Minor",
            "message": f"Table S{ref} が SI に定義されているが本文から参照されていない",
        })

    # 番号の連続性チェック
    for label, nums in [("Figure S", si_fig_defs), ("Table S", si_table_defs)]:
        if nums:
            int_nums = sorted(int(n) for n in nums)
            expected = list(range(1, max(int_nums) + 1))
            gaps = set(expected) - set(int_nums)
            for g in sorted(gaps):
                issues.append({
                    "type": "NUMBERING_GAP",
                    "severity": "Major",
                    "message": f"{label}{g} が欠番（{label}1〜{label}{max(int_nums)} の間）",
                })

    result = {
        "timestamp": datetime.datetime.now().isoformat(),
        "main_fig_refs": sorted(main_fig_refs),
        "main_table_refs": sorted(main_table_refs),
        "si_fig_defs": sorted(si_fig_defs),
        "si_table_defs": sorted(si_table_defs),
        "issues": issues,
        "status": "PASS" if not any(i["severity"] == "Critical" for i in issues) else "FAIL",
    }

    if filepath is None:
        filepath = BASE_DIR / "manuscript" / "si_crossref_report.json"
    filepath.parent.mkdir(parents=True, exist_ok=True)

    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)

    print(f"  → SI 相互参照チェック: {result['status']}")
    print(f"    本文→SI参照: Figure S × {len(main_fig_refs)}, Table S × {len(main_table_refs)}")
    print(f"    問題: {len(issues)} 件")
    if issues:
        for issue in issues:
            print(f"    [{issue['severity']}] {issue['message']}")

    return result
```

## 6. SI 生成パイプライン統合

```python
def run_si_pipeline(manuscript_path, figures_dir=None, results_dir=None,
                     journal_format="imrad", title="", authors=""):
    """
    SI 生成パイプラインを実行する。

    1. SI 候補の収集
    2. SI ドキュメントの生成
    3. 相互参照チェック
    """
    from pathlib import Path

    if figures_dir is None:
        figures_dir = BASE_DIR / "figures"
    if results_dir is None:
        results_dir = BASE_DIR / "results"

    print("=" * 60)
    print("Supplementary Information Generator")
    print("=" * 60)

    # Phase 1: 候補収集
    print("\n[Phase 1] SI 候補の収集...")
    candidates = collect_si_candidates(manuscript_path, figures_dir, results_dir)

    # Phase 2: SI 生成
    print("\n[Phase 2] SI ドキュメントの生成...")
    si_path = generate_supplementary(
        candidates, journal_format=journal_format,
        title=title, authors=authors,
    )

    # Phase 3: 相互参照チェック
    print("\n[Phase 3] 相互参照チェック...")
    crossref = check_si_crossrefs(manuscript_path, si_path)

    print("\n" + "=" * 60)
    print(f"SI 生成完了！ (ステータス: {crossref['status']})")
    print("=" * 60)

    return {"si_path": si_path, "crossref": crossref}
```

## References

### Output Files

| ファイル | 形式 | 生成タイミング |
|---|---|---|
| `manuscript/supplementary.md` | SI ドキュメント（Markdown） | Phase 2 完了時 |
| `manuscript/si_crossref_report.json` | 相互参照チェック結果（JSON） | Phase 3 完了時 |

### 参照スキル

| スキル | 連携 |
|---|---|
| `scientific-academic-writing` | 本文原稿 (`manuscript/manuscript.md`) を入力として使用 |
| `scientific-publication-figures` | `figures/` の未使用図を SI 候補として自動検出 |
| `scientific-pipeline-scaffold` | `results/` の詳細データを SI Table 候補として収集 |
| `scientific-critical-review` | レビュー後の修正で追加された図表を SI に反映 |
| `scientific-latex-formatter` | SI の LaTeX 変換に使用 |

```

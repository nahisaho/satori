---
name: scientific-reproducible-reporting
description: |
  再現可能レポーティングスキル。Quarto 科学文書・
  Jupyter Book 多章構成・Papermill パラメトリック実行・
  nbconvert 自動変換・Sphinx-Gallery コード例ドキュメント。
tu_tools:
  - key: biotools
    name: bio.tools
    description: 再現可能レポーティングツール検索
---

# Scientific Reproducible Reporting

再現可能な科学レポート・文書生成パイプラインを提供し、
コード → 実行 → 文書化の自動化を実現する。

## When to Use

- Quarto で再現可能な科学文書を作成するとき
- Jupyter Book で多章構成の文書を構築するとき
- Papermill でパラメトリック実行を自動化するとき
- nbconvert でノートブックを各種形式に変換するとき
- CI/CD で解析レポートを自動生成するとき
- 複数パラメータセットで解析を繰り返し実行するとき

---

## Quick Start

## 1. Quarto 科学文書テンプレート

```python
import os


def generate_quarto_document(title="Scientific Analysis Report",
                             author="SATORI",
                             format_type="html",
                             output_dir="quarto_project"):
    """
    Quarto 科学文書テンプレート生成。

    Parameters:
        title: str — ドキュメントタイトル
        author: str — 著者名
        format_type: str — "html" / "pdf" / "docx" / "revealjs"
        output_dir: str — 出力ディレクトリ
    """
    os.makedirs(output_dir, exist_ok=True)

    # _quarto.yml
    quarto_config = f"""project:
  type: default
  output-dir: _output

format:
  {format_type}:
    toc: true
    toc-depth: 3
    number-sections: true
    code-fold: true
    code-tools: true
    theme: cosmo

execute:
  echo: true
  warning: false
  cache: true

bibliography: references.bib
csl: nature.csl
"""

    # メインドキュメント
    main_qmd = f"""---
title: "{title}"
author: "{author}"
date: today
format:
  {format_type}:
    code-fold: true
    code-tools: true
jupyter: python3
---

## はじめに

このレポートは SATORI スキルを用いた再現可能な科学解析文書です。

```{{python}}
#| label: setup
#| echo: false

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

# パラメータ (Papermill 互換)
n_samples = 1000
random_seed = 42
```

## データ概要

```{{python}}
#| label: data-summary
#| tbl-cap: "データセット概要"

np.random.seed(random_seed)
df = pd.DataFrame({{
    "x": np.random.randn(n_samples),
    "y": np.random.randn(n_samples),
    "group": np.random.choice(["A", "B", "C"], n_samples)
}})
df.describe()
```

## 可視化

```{{python}}
#| label: fig-scatter
#| fig-cap: "散布図"

fig, ax = plt.subplots(figsize=(8, 6))
for g, sub in df.groupby("group"):
    ax.scatter(sub["x"], sub["y"], label=g, alpha=0.6, s=20)
ax.legend()
ax.set_xlabel("X")
ax.set_ylabel("Y")
plt.show()
```

## 結論

解析結果のサマリーを記載する。

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `biotools` | bio.tools | 再現可能レポーティングツール検索 |

## References
"""

    # references.bib (空テンプレート)
    bib_template = """@article{example2024,
  title={Example Reference},
  author={Author, A.},
  journal={Journal},
  year={2024}
}
"""

    with open(os.path.join(output_dir, "_quarto.yml"), "w") as f:
        f.write(quarto_config)
    with open(os.path.join(output_dir, "report.qmd"), "w") as f:
        f.write(main_qmd)
    with open(os.path.join(output_dir, "references.bib"), "w") as f:
        f.write(bib_template)

    print(f"Quarto project → {output_dir}/")
    print(f"  Build: cd {output_dir} && quarto render report.qmd")
    return output_dir
```

## 2. Papermill パラメトリック実行

```python
def papermill_parametric_run(template_notebook, output_dir,
                             parameter_sets, kernel="python3"):
    """
    Papermill パラメトリック実行 — 複数パラメータセットで自動実行。

    Parameters:
        template_notebook: str — テンプレートノートブックパス
        output_dir: str — 出力ディレクトリ
        parameter_sets: list[dict] — パラメータセットのリスト
        kernel: str — カーネル名
    """
    import papermill as pm

    os.makedirs(output_dir, exist_ok=True)
    results = []

    for i, params in enumerate(parameter_sets):
        output_path = os.path.join(output_dir, f"run_{i:03d}.ipynb")
        try:
            pm.execute_notebook(
                template_notebook,
                output_path,
                parameters=params,
                kernel_name=kernel)
            results.append({
                "run": i, "params": params,
                "output": output_path, "status": "success"})
        except Exception as e:
            results.append({
                "run": i, "params": params,
                "output": output_path, "status": f"error: {str(e)}"})

    import pandas as pd
    results_df = pd.DataFrame(results)
    n_success = (results_df["status"] == "success").sum()
    print(f"Papermill: {n_success}/{len(parameter_sets)} runs succeeded")
    return results_df
```

## 3. Jupyter Book 多章構成

```python
def generate_jupyter_book(title="Scientific Analysis Book",
                          chapters=None,
                          output_dir="jupyter_book"):
    """
    Jupyter Book プロジェクトテンプレート生成。

    Parameters:
        title: str — 書籍タイトル
        chapters: list[dict] | None — 章情報 [{"title": ..., "file": ...}]
        output_dir: str — 出力ディレクトリ
    """
    os.makedirs(output_dir, exist_ok=True)

    if chapters is None:
        chapters = [
            {"title": "Introduction", "file": "intro"},
            {"title": "Data Loading", "file": "ch01_data"},
            {"title": "Exploratory Analysis", "file": "ch02_eda"},
            {"title": "Modeling", "file": "ch03_model"},
            {"title": "Results", "file": "ch04_results"},
        ]

    # _config.yml
    config = f"""title: "{title}"
author: SATORI
execute:
  execute_notebooks: auto
  timeout: 600
repository:
  url: ""
launch_buttons:
  binderhub_url: ""
sphinx:
  extra_extensions:
    - sphinx_proof
"""

    # _toc.yml
    toc_entries = "\n".join(
        [f"  - file: {ch['file']}" for ch in chapters])
    toc = f"""format: jb-book
root: intro
chapters:
{toc_entries}
"""

    with open(os.path.join(output_dir, "_config.yml"), "w") as f:
        f.write(config)
    with open(os.path.join(output_dir, "_toc.yml"), "w") as f:
        f.write(toc)

    # 各章テンプレート
    for ch in chapters:
        filepath = os.path.join(output_dir, f"{ch['file']}.md")
        if not os.path.exists(filepath):
            content = f"# {ch['title']}\n\nThis chapter covers {ch['title'].lower()}.\n"
            with open(filepath, "w") as f:
                f.write(content)

    print(f"Jupyter Book → {output_dir}/")
    print(f"  Build: jupyter-book build {output_dir}")
    return output_dir
```

## 4. nbconvert 自動変換

```python
def batch_convert_notebooks(notebook_dir, output_format="html",
                            output_dir=None, execute=True):
    """
    ノートブック一括変換。

    Parameters:
        notebook_dir: str — ノートブックディレクトリ
        output_format: str — "html" / "pdf" / "markdown" / "script"
        output_dir: str | None — 出力先 (None=同ディレクトリ)
        execute: bool — 実行後に変換
    """
    import subprocess
    import glob

    notebooks = sorted(glob.glob(os.path.join(notebook_dir, "*.ipynb")))
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    results = []
    for nb_path in notebooks:
        cmd = ["jupyter", "nbconvert", f"--to={output_format}"]
        if execute:
            cmd.append("--execute")
        if output_dir:
            cmd.extend(["--output-dir", output_dir])
        cmd.append(nb_path)

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            results.append({"notebook": nb_path, "status": "success"})
        except subprocess.CalledProcessError as e:
            results.append({"notebook": nb_path, "status": f"error: {e.stderr[:100]}"})

    import pandas as pd
    results_df = pd.DataFrame(results)
    n_ok = (results_df["status"] == "success").sum()
    print(f"nbconvert ({output_format}): {n_ok}/{len(notebooks)} converted")
    return results_df
```

---

## パイプライン統合

```
[解析完了] → reproducible-reporting → presentation-design
              (レポート自動生成)         (プレゼン作成)
                     │                        ↓
            interactive-dashboard      academic-writing
              (ダッシュボード)           (論文執筆)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `quarto_project/` | Quarto プロジェクト | → quarto render |
| `papermill_runs/` | パラメトリック実行結果 | → 集計 |
| `jupyter_book/` | Jupyter Book プロジェクト | → jb build |

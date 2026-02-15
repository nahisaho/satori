---
name: scientific-interactive-dashboard
description: |
  インタラクティブダッシュボードスキル。
  Streamlit / Dash / Panel / Voilà による
  科学データダッシュボード構築・リアルタイムパラメータ探索 UI ・
  ウィジェット連動・データアップロード・解析パイプライン UI 化。
tu_tools:
  - key: biotools
    name: bio.tools
    description: インタラクティブ可視化ツール検索
---

# Scientific Interactive Dashboard

科学データのインタラクティブダッシュボードを構築し、
パラメータ探索・結果共有・リアルタイム解析を実現する。

## When to Use

- Streamlit で迅速にデータ探索ダッシュボードを構築するとき
- Dash でカスタマイズ性の高い解析 UI を作成するとき
- Panel / Voilà で Jupyter ノートブックをダッシュボード化するとき
- パラメータスライダー + リアルタイム更新の UI を実装するとき
- 複数人で解析結果を共有するとき
- 非プログラマーに解析ツールを提供するとき

---

## Quick Start

## 1. Streamlit 科学データダッシュボード

```python
def generate_streamlit_dashboard(output_path="dashboard_app.py"):
    """
    Streamlit ダッシュボードテンプレート生成。

    Parameters:
        output_path: str — 出力 Python ファイル
    """
    code = '''
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px


st.set_page_config(page_title="Scientific Data Dashboard",
                   layout="wide", page_icon="🔬")

st.title("🔬 Scientific Data Dashboard")

# --- サイドバー: データアップロード & パラメータ ---
st.sidebar.header("Settings")

uploaded_file = st.sidebar.file_uploader(
    "Upload CSV / Excel", type=["csv", "xlsx"])

if uploaded_file is not None:
    if uploaded_file.name.endswith(".csv"):
        df = pd.read_csv(uploaded_file)
    else:
        df = pd.read_excel(uploaded_file)
else:
    # デモデータ
    np.random.seed(42)
    n = 500
    df = pd.DataFrame({
        "x": np.random.randn(n),
        "y": np.random.randn(n),
        "z": np.random.randn(n),
        "category": np.random.choice(["A", "B", "C"], n),
        "value": np.random.exponential(2, n)
    })
    st.sidebar.info("Demo data loaded (upload your own CSV)")

# --- データ概要 ---
col1, col2, col3 = st.columns(3)
col1.metric("Rows", len(df))
col2.metric("Columns", len(df.columns))
col3.metric("Missing", int(df.isnull().sum().sum()))

# --- タブ ---
tab1, tab2, tab3, tab4 = st.tabs(
    ["📊 Explorer", "📈 Distribution", "🔗 Correlation", "📋 Data"])

numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
cat_cols = df.select_dtypes(exclude=[np.number]).columns.tolist()

with tab1:
    st.subheader("Interactive Explorer")
    c1, c2 = st.columns(2)
    x_col = c1.selectbox("X axis", numeric_cols, index=0)
    y_col = c2.selectbox("Y axis", numeric_cols,
                          index=min(1, len(numeric_cols)-1))
    color_col = st.selectbox("Color", [None] + cat_cols + numeric_cols)

    fig = px.scatter(df, x=x_col, y=y_col, color=color_col,
                     opacity=0.7, title=f"{x_col} vs {y_col}")
    st.plotly_chart(fig, use_container_width=True)

with tab2:
    st.subheader("Distribution Analysis")
    dist_col = st.selectbox("Column", numeric_cols, key="dist")
    n_bins = st.slider("Bins", 10, 100, 30)
    fig2 = px.histogram(df, x=dist_col, nbins=n_bins,
                        marginal="box", title=f"Distribution: {dist_col}")
    st.plotly_chart(fig2, use_container_width=True)

with tab3:
    st.subheader("Correlation Matrix")
    corr = df[numeric_cols].corr()
    fig3 = px.imshow(corr, text_auto=".2f", color_continuous_scale="RdBu_r",
                     title="Correlation Heatmap")
    st.plotly_chart(fig3, use_container_width=True)

with tab4:
    st.subheader("Raw Data")
    st.dataframe(df, use_container_width=True)
    csv = df.to_csv(index=False)
    st.download_button("Download CSV", csv, "data.csv", "text/csv")
'''

    with open(output_path, "w") as f:
        f.write(code)

    print(f"Streamlit dashboard → {output_path}")
    print(f"  Run: streamlit run {output_path}")
    return output_path
```

## 2. Dash コールバックダッシュボード

```python
def generate_dash_dashboard(output_path="dash_app.py"):
    """
    Dash ダッシュボードテンプレート生成。

    Parameters:
        output_path: str — 出力 Python ファイル
    """
    code = '''
from dash import Dash, html, dcc, Input, Output, dash_table
import pandas as pd
import numpy as np
import plotly.express as px

app = Dash(__name__)

# デモデータ
np.random.seed(42)
n = 500
df = pd.DataFrame({
    "x": np.random.randn(n),
    "y": np.random.randn(n),
    "z": np.random.randn(n),
    "group": np.random.choice(["Control", "Treatment A", "Treatment B"], n),
    "response": np.random.exponential(2, n)
})

numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()

app.layout = html.Div([
    html.H1("Scientific Data Dashboard", style={"textAlign": "center"}),

    html.Div([
        html.Div([
            html.Label("X Axis"),
            dcc.Dropdown(id="x-col", options=numeric_cols,
                         value=numeric_cols[0])
        ], style={"width": "30%", "display": "inline-block"}),
        html.Div([
            html.Label("Y Axis"),
            dcc.Dropdown(id="y-col", options=numeric_cols,
                         value=numeric_cols[1])
        ], style={"width": "30%", "display": "inline-block"}),
        html.Div([
            html.Label("Color"),
            dcc.Dropdown(id="color-col",
                         options=df.columns.tolist(),
                         value="group")
        ], style={"width": "30%", "display": "inline-block"}),
    ], style={"padding": "20px"}),

    html.Div([
        html.Div([dcc.Graph(id="scatter-plot")],
                 style={"width": "50%", "display": "inline-block"}),
        html.Div([dcc.Graph(id="histogram")],
                 style={"width": "50%", "display": "inline-block"}),
    ]),

    html.Div([
        html.H3("Summary Statistics"),
        dash_table.DataTable(
            id="summary-table",
            columns=[{"name": c, "id": c}
                     for c in ["stat"] + numeric_cols],
            style_table={"overflowX": "auto"})
    ], style={"padding": "20px"})
])

@app.callback(
    [Output("scatter-plot", "figure"),
     Output("histogram", "figure"),
     Output("summary-table", "data")],
    [Input("x-col", "value"),
     Input("y-col", "value"),
     Input("color-col", "value")]
)
def update_plots(x_col, y_col, color_col):
    fig1 = px.scatter(df, x=x_col, y=y_col, color=color_col,
                      opacity=0.7, title=f"{x_col} vs {y_col}")
    fig2 = px.histogram(df, x=x_col, color=color_col,
                        marginal="box", barmode="overlay", opacity=0.7)
    stats = df[numeric_cols].describe().reset_index()
    stats.columns = ["stat"] + numeric_cols
    return fig1, fig2, stats.to_dict("records")

if __name__ == "__main__":
    app.run(debug=True, port=8050)
'''

    with open(output_path, "w") as f:
        f.write(code)

    print(f"Dash dashboard → {output_path}")
    print(f"  Run: python {output_path}")
    return output_path
```

## 3. Panel ダッシュボード

```python
def generate_panel_dashboard(output_path="panel_app.py"):
    """
    Panel ダッシュボードテンプレート生成。

    Parameters:
        output_path: str — 出力 Python ファイル
    """
    code = '''
import panel as pn
import pandas as pd
import numpy as np
import plotly.express as px

pn.extension("plotly")

# デモデータ
np.random.seed(42)
n = 500
df = pd.DataFrame({
    "x": np.random.randn(n),
    "y": np.random.randn(n),
    "z": np.random.randn(n),
    "group": np.random.choice(["A", "B", "C"], n),
    "value": np.random.exponential(2, n)
})

numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()

# ウィジェット
x_select = pn.widgets.Select(name="X Axis", options=numeric_cols, value="x")
y_select = pn.widgets.Select(name="Y Axis", options=numeric_cols, value="y")
n_bins = pn.widgets.IntSlider(name="Histogram Bins", start=10, end=100, value=30)


@pn.depends(x_select, y_select)
def scatter_plot(x_col, y_col):
    fig = px.scatter(df, x=x_col, y=y_col, color="group",
                     opacity=0.7, title=f"{x_col} vs {y_col}")
    return fig


@pn.depends(x_select, n_bins)
def hist_plot(x_col, bins):
    fig = px.histogram(df, x=x_col, nbins=bins, color="group",
                       barmode="overlay", opacity=0.7)
    return fig


dashboard = pn.template.FastListTemplate(
    title="Scientific Data Dashboard",
    sidebar=[x_select, y_select, n_bins],
    main=[
        pn.Row(pn.pane.Plotly(scatter_plot, sizing_mode="stretch_width"),
               pn.pane.Plotly(hist_plot, sizing_mode="stretch_width")),
        pn.pane.DataFrame(df.describe().T, sizing_mode="stretch_width")
    ]
)

dashboard.servable()
'''

    with open(output_path, "w") as f:
        f.write(code)

    print(f"Panel dashboard → {output_path}")
    print(f"  Run: panel serve {output_path}")
    return output_path
```

## 4. ダッシュボード比較ガイド

```python
def compare_dashboard_frameworks():
    """
    Streamlit / Dash / Panel / Voilà 比較表を出力。
    """
    comparison = pd.DataFrame({
        "Framework": ["Streamlit", "Dash", "Panel", "Voilà"],
        "Ease_of_Use": ["★★★★★", "★★★☆☆", "★★★★☆", "★★★★★"],
        "Customization": ["★★★☆☆", "★★★★★", "★★★★☆", "★★☆☆☆"],
        "Interactivity": ["★★★★☆", "★★★★★", "★★★★★", "★★★☆☆"],
        "Performance": ["★★★☆☆", "★★★★★", "★★★★☆", "★★★☆☆"],
        "Deployment": ["Streamlit Cloud", "Heroku/AWS", "Any ASGI", "Binder/Hub"],
        "Best_For": [
            "Rapid prototyping, data exploration",
            "Production apps, complex callbacks",
            "Jupyter integration, scientific viz",
            "Notebook → dashboard conversion"
        ]
    })

    print("=== Dashboard Framework Comparison ===")
    print(comparison.to_string(index=False))
    return comparison
```

---

## パイプライン統合

```
advanced-visualization → interactive-dashboard → presentation-design
    (高度可視化)           (ダッシュボード)          (プレゼン)
         │                       │                      ↓
  missing-data-analysis ────────┘            scientific-schematics
    (欠損値解析)                                (図式デザイン)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `dashboard_app.py` | Streamlit ダッシュボード | → deployment |
| `dash_app.py` | Dash ダッシュボード | → deployment |
| `panel_app.py` | Panel ダッシュボード | → deployment |
| `framework_comparison.csv` | フレームワーク比較 | → 選択指針 |

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `biotools` | bio.tools | インタラクティブ可視化ツール検索 |

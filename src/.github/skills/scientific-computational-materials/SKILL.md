---
name: scientific-computational-materials
description: |
  計算材料科学スキル。pymatgen による結晶構造操作・対称性解析、
  Materials Project API による材料データベース照会、
  相図計算 (凸包解析)、電子バンド構造・状態密度 (DOS) 可視化、
  VASP/Quantum ESPRESSO 入出力、高スループットスクリーニングパイプライン。
---

# Scientific Computational Materials

無機結晶材料を中心に、結晶構造操作・物性データベース照会・
第一原理計算入出力・相安定性解析・電子構造可視化を提供する
計算材料科学パイプライン。

## When to Use

- 結晶構造の生成・変換・対称性解析を行うとき
- Materials Project API で材料物性データを照会するとき
- 相図 (凸包) による安定性解析が必要なとき
- DFT 計算 (VASP/QE) の入力ファイル生成・出力解析をするとき
- バンド構造・状態密度 (DOS) を可視化するとき
- 高スループット材料スクリーニングを行うとき

---

## Quick Start

## 1. 結晶構造操作

```python
from pymatgen.core import Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def create_crystal_structure(lattice_params, species, coords,
                               output_file="results/structure.cif"):
    """
    pymatgen による結晶構造の生成と解析。

    Parameters:
    - lattice_params: (a, b, c, alpha, beta, gamma) or 3×3 matrix
    - species: 元素リスト ["Si", "O", ...]
    - coords: 分率座標 [[x,y,z], ...]
    """
    lattice = Lattice.from_parameters(*lattice_params)
    structure = Structure(lattice, species, coords)

    sga = SpacegroupAnalyzer(structure, symprec=0.1)

    print(f"  Crystal structure:")
    print(f"    Formula: {structure.composition.reduced_formula}")
    print(f"    Space group: {sga.get_space_group_symbol()} ({sga.get_space_group_number()})")
    print(f"    Crystal system: {sga.get_crystal_system()}")
    print(f"    Lattice: a={lattice.a:.3f}, b={lattice.b:.3f}, c={lattice.c:.3f}")
    print(f"    Volume: {structure.volume:.3f} Å³")
    print(f"    Density: {structure.density:.3f} g/cm³")

    # 慣用セルに変換
    conventional = sga.get_conventional_standard_structure()
    primitive = sga.get_primitive_standard_structure()

    print(f"    Conventional cell: {len(conventional)} atoms")
    print(f"    Primitive cell: {len(primitive)} atoms")

    # CIF 出力
    structure.to(filename=output_file)

    return structure, sga


def analyze_structure_symmetry(structure, symprec=0.01):
    """
    結晶対称性の詳細解析。

    - 空間群・点群
    - Wyckoff 位置
    - サイト対称性
    """
    sga = SpacegroupAnalyzer(structure, symprec=symprec)

    symmetry_info = {
        "space_group_symbol": sga.get_space_group_symbol(),
        "space_group_number": sga.get_space_group_number(),
        "crystal_system": sga.get_crystal_system(),
        "point_group": sga.get_point_group_symbol(),
        "hall_symbol": sga.get_hall(),
        "symmetry_operations": len(sga.get_symmetry_operations()),
    }

    # Wyckoff 位置
    sym_struct = sga.get_symmetrized_structure()
    wyckoff_sites = sym_struct.wyckoff_symbols

    print(f"  Symmetry analysis:")
    for k, v in symmetry_info.items():
        print(f"    {k}: {v}")
    print(f"    Wyckoff sites: {wyckoff_sites}")

    return symmetry_info
```

## 2. Materials Project API 照会

```python
from mp_api.client import MPRester
import pandas as pd


def query_materials_project(formula=None, elements=None,
                              band_gap_range=None,
                              e_above_hull_max=0.025):
    """
    Materials Project API による材料データベース照会。

    検索条件:
    - 化学式 (formula) — exact or reduced
    - 構成元素 (elements) — 含む/排除
    - バンドギャップ範囲 (band_gap_range) — eV
    - 凸包上エネルギー (e_above_hull) — 安定性指標
    """
    with MPRester() as mpr:
        criteria = {}
        if formula:
            criteria["formula"] = formula
        if elements:
            criteria["elements"] = elements

        docs = mpr.materials.summary.search(
            **criteria,
            energy_above_hull=(0, e_above_hull_max) if e_above_hull_max else None,
            band_gap=band_gap_range,
            fields=[
                "material_id", "formula_pretty", "volume",
                "density", "band_gap", "energy_above_hull",
                "formation_energy_per_atom", "is_stable",
                "symmetry", "nsites",
            ],
        )

    results = []
    for doc in docs:
        results.append({
            "mp_id": doc.material_id,
            "formula": doc.formula_pretty,
            "space_group": doc.symmetry.symbol if doc.symmetry else None,
            "band_gap_eV": doc.band_gap,
            "e_above_hull_eV": doc.energy_above_hull,
            "formation_energy_eV": doc.formation_energy_per_atom,
            "density_g_cm3": doc.density,
            "nsites": doc.nsites,
            "is_stable": doc.is_stable,
        })

    df = pd.DataFrame(results)
    print(f"  Materials Project query:")
    print(f"    Found: {len(df)} materials")
    if len(df) > 0:
        print(f"    Stable: {df['is_stable'].sum()}")
        print(f"    Band gap range: {df['band_gap_eV'].min():.2f}–{df['band_gap_eV'].max():.2f} eV")

    return df
```

## 3. 相図・凸包解析

```python
import numpy as np
import pandas as pd


def compute_phase_diagram(system_elements, output_file="figures/phase_diagram.png"):
    """
    相図 (凸包) 計算。

    凸包 (Convex Hull):
    - 安定相: 凸包上の点 (e_above_hull = 0)
    - 準安定相: 凸包上方の点 (e_above_hull > 0)
    - 分解反応: 凸包上の隣接安定相への分解
    """
    from mp_api.client import MPRester
    from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter

    with MPRester() as mpr:
        entries = mpr.get_entries_in_chemsys(system_elements)

    pd_obj = PhaseDiagram(entries)

    print(f"  Phase diagram: {'-'.join(system_elements)}")
    print(f"    Total entries: {len(entries)}")
    print(f"    Stable phases: {len(pd_obj.stable_entries)}")

    for entry in pd_obj.stable_entries:
        formula = entry.composition.reduced_formula
        e_form = pd_obj.get_form_energy_per_atom(entry)
        print(f"    {formula}: ΔHf = {e_form:.4f} eV/atom")

    # 可視化
    plotter = PDPlotter(pd_obj)
    plotter.get_plot().savefig(output_file, dpi=300, bbox_inches="tight")

    return pd_obj
```

## 4. 電子バンド構造・DOS

```python
import numpy as np


def plot_band_structure(material_id, output_file="figures/band_structure.png"):
    """
    電子バンド構造の取得と可視化。

    - 高対称 k-path (Setyawan-Curtarolo 規約)
    - バンドギャップ判定 (直接/間接)
    - フェルミレベル基準
    """
    from mp_api.client import MPRester
    from pymatgen.electronic_structure.plotter import BSPlotter

    with MPRester() as mpr:
        bs = mpr.get_bandstructure_by_material_id(material_id)

    if bs is None:
        print(f"  No band structure available for {material_id}")
        return None

    gap = bs.get_band_gap()
    print(f"  Band structure: {material_id}")
    print(f"    Band gap: {gap['energy']:.3f} eV")
    print(f"    Direct: {gap['direct']}")
    print(f"    Transition: {gap['transition']}")

    plotter = BSPlotter(bs)
    plotter.get_plot().savefig(output_file, dpi=300, bbox_inches="tight")

    return bs


def plot_density_of_states(material_id, output_file="figures/dos.png"):
    """
    電子状態密度 (DOS) の取得と可視化。

    - Total DOS + projected DOS (元素分解)
    - スピン偏極 (該当時)
    """
    from mp_api.client import MPRester
    from pymatgen.electronic_structure.plotter import DosPlotter

    with MPRester() as mpr:
        dos = mpr.get_dos_by_material_id(material_id)

    if dos is None:
        print(f"  No DOS available for {material_id}")
        return None

    print(f"  DOS: {material_id}")
    print(f"    Efermi: {dos.efermi:.3f} eV")

    plotter = DosPlotter()
    plotter.add_dos("Total", dos)
    plotter.get_plot().savefig(output_file, dpi=300, bbox_inches="tight")

    return dos
```

## 5. VASP/Quantum ESPRESSO 入出力

```python
from pymatgen.core import Structure


def generate_vasp_inputs(structure, output_dir="vasp_inputs",
                          calculation_type="relaxation"):
    """
    VASP 入力ファイル生成。

    計算タイプ:
    - relaxation: 構造緩和 (ISIF=3)
    - static: 静的計算 (NSW=0)
    - band: バンド構造 (ICHARG=11)
    - dos: 状態密度 (LORBIT=11)
    """
    from pymatgen.io.vasp.sets import (
        MPRelaxSet, MPStaticSet,
    )
    import os
    os.makedirs(output_dir, exist_ok=True)

    if calculation_type == "relaxation":
        vis = MPRelaxSet(structure)
    elif calculation_type == "static":
        vis = MPStaticSet(structure)
    else:
        vis = MPRelaxSet(structure)

    vis.write_input(output_dir)

    print(f"  VASP inputs generated: {output_dir}")
    print(f"    Calculation: {calculation_type}")
    print(f"    Files: INCAR, POSCAR, POTCAR, KPOINTS")

    return output_dir


def parse_vasp_output(vasprun_file="vasprun.xml"):
    """
    VASP 出力 (vasprun.xml) の解析。
    """
    from pymatgen.io.vasp.outputs import Vasprun

    vr = Vasprun(vasprun_file, parse_dos=True, parse_eigen=True)

    print(f"  VASP output: {vasprun_file}")
    print(f"    Final energy: {vr.final_energy:.6f} eV")
    print(f"    Converged: {vr.converged}")
    print(f"    Band gap: {vr.get_band_structure().get_band_gap()['energy']:.3f} eV")

    return vr
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/structure.cif` | CIF |
| `results/materials_query.csv` | CSV |
| `figures/phase_diagram.png` | PNG |
| `figures/band_structure.png` | PNG |
| `figures/dos.png` | PNG |
| `vasp_inputs/` | VASP input set |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

なし — Materials Project API (mp-api) を直接利用。

### 参照スキル

| スキル | 関連 |
|---|---|
| `scientific-quantum-computing` | 量子計算・VQE・量子化学 |
| `scientific-cheminformatics` | 分子記述子・構造解析 |
| `scientific-publication-figures` | 構造・相図可視化 |
| `scientific-materials-characterization` | XRD・SEM・実験材料特性 |

### 依存パッケージ

`pymatgen`, `mp-api`, `pandas`, `numpy`, `matplotlib`

---
name: scientific-md-simulation
description: |
  分子動力学シミュレーション解析スキル。MDAnalysis によるトラジェクトリ解析・
  RMSD/RMSF/Rg 時系列指標・水素結合解析・二次構造変化追跡・
  OpenFF Toolkit力場パラメータ化・溶媒和自由エネルギー推定パイプライン。
tu_tools:
  - key: pdb
    name: PDB
    description: 分子構造データベース参照
---

# Scientific MD Simulation

MDAnalysis と OpenFF Toolkit を活用した分子動力学 (MD) シミュレーション
解析パイプラインを提供する。トラジェクトリ読込みから構造指標計算、
水素結合解析、力場パラメータ化まで統合。

## When to Use

- MD トラジェクトリ (DCD/XTC/TRR/GRO) を解析するとき
- RMSD/RMSF/Radius of Gyration 時系列を計算するとき
- 水素結合パターンを解析するとき
- OpenFF で小分子の力場パラメータを自動生成するとき
- タンパク質-リガンド複合体の安定性を評価するとき
- 溶媒接触表面積 (SASA) を計算するとき

---

## Quick Start

## 1. トラジェクトリ読込み & 基本情報

```python
import MDAnalysis as mda
import numpy as np
import pandas as pd


def load_trajectory(topology, trajectory):
    """
    MD トラジェクトリ読込み。

    Parameters:
        topology: str — トポロジファイル (PSF/PDB/GRO/PRMTOP)
        trajectory: str — トラジェクトリファイル (DCD/XTC/TRR)

    K-Dense: mdanalysis
    """
    u = mda.Universe(topology, trajectory)

    info = {
        "n_atoms": u.atoms.n_atoms,
        "n_residues": u.residues.n_residues,
        "n_segments": u.segments.n_segments,
        "n_frames": u.trajectory.n_frames,
        "dt_ps": u.trajectory.dt,
        "total_time_ns": u.trajectory.n_frames * u.trajectory.dt / 1000,
        "topology_format": topology.split(".")[-1],
        "trajectory_format": trajectory.split(".")[-1],
    }

    print(f"Loaded: {info['n_atoms']} atoms, {info['n_frames']} frames, "
          f"{info['total_time_ns']:.1f} ns")
    return u, info
```

## 2. RMSD 解析

```python
from MDAnalysis.analysis.rms import RMSD


def compute_rmsd(universe, selection="backbone", ref_frame=0):
    """
    RMSD (Root Mean Square Deviation) 計算。

    Parameters:
        universe: mda.Universe — MD ユニバース
        selection: str — 原子選択文字列
        ref_frame: int — 参照フレーム
    """
    R = RMSD(universe, universe, select=selection, ref_frame=ref_frame)
    R.run()

    df = pd.DataFrame({
        "frame": R.results.rmsd[:, 0].astype(int),
        "time_ps": R.results.rmsd[:, 1],
        "rmsd_A": R.results.rmsd[:, 2],
    })
    df["time_ns"] = df["time_ps"] / 1000

    print(f"RMSD ({selection}): mean={df['rmsd_A'].mean():.2f} Å, "
          f"max={df['rmsd_A'].max():.2f} Å")
    return df
```

## 3. RMSF 解析

```python
from MDAnalysis.analysis.rms import RMSF as RMSFAnalysis


def compute_rmsf(universe, selection="name CA"):
    """
    RMSF (Root Mean Square Fluctuation) per residue。

    Parameters:
        universe: mda.Universe — MD ユニバース
        selection: str — 原子選択 (通常 Cα)
    """
    atoms = universe.select_atoms(selection)
    R = RMSFAnalysis(atoms).run()

    df = pd.DataFrame({
        "resid": atoms.resids,
        "resname": atoms.resnames,
        "rmsf_A": R.results.rmsf,
    })

    # 柔軟領域の同定
    threshold = df["rmsf_A"].mean() + 2 * df["rmsf_A"].std()
    df["flexible"] = df["rmsf_A"] > threshold

    print(f"RMSF: mean={df['rmsf_A'].mean():.2f} Å, "
          f"flexible residues={df['flexible'].sum()}")
    return df
```

## 4. Radius of Gyration

```python
def compute_radius_of_gyration(universe, selection="protein"):
    """
    Radius of Gyration (Rg) 時系列計算。

    Parameters:
        universe: mda.Universe — MD ユニバース
        selection: str — 原子選択
    """
    protein = universe.select_atoms(selection)
    rg_data = []

    for ts in universe.trajectory:
        rg = protein.radius_of_gyration()
        rg_data.append({
            "frame": ts.frame,
            "time_ns": ts.time / 1000,
            "rg_A": rg,
        })

    df = pd.DataFrame(rg_data)
    print(f"Rg: mean={df['rg_A'].mean():.2f} Å, "
          f"std={df['rg_A'].std():.2f} Å")
    return df
```

## 5. 水素結合解析

```python
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis


def hydrogen_bond_analysis(universe, donor_sel="protein", acceptor_sel="protein",
                           d_a_cutoff=3.0, angle_cutoff=150):
    """
    水素結合解析。

    Parameters:
        universe: mda.Universe — MD ユニバース
        donor_sel: str — ドナー選択
        acceptor_sel: str — アクセプター選択
        d_a_cutoff: float — D-A 距離閾値 (Å)
        angle_cutoff: float — D-H-A 角度閾値 (°)
    """
    hbonds = HydrogenBondAnalysis(
        universe,
        donors_sel=donor_sel,
        acceptors_sel=acceptor_sel,
        d_a_cutoff=d_a_cutoff,
        d_h_a_angle_cutoff=angle_cutoff,
    )
    hbonds.run()

    # フレームあたり水素結合数
    counts = hbonds.count_by_time()
    df_counts = pd.DataFrame({
        "time_ps": counts[:, 0],
        "n_hbonds": counts[:, 1].astype(int),
    })

    print(f"H-bonds: mean={df_counts['n_hbonds'].mean():.1f}/frame, "
          f"total unique={len(hbonds.results.hbonds)}")
    return hbonds, df_counts
```

## 6. SASA (溶媒接触表面積)

```python
from MDAnalysis.analysis.sasa import SASA


def compute_sasa(universe, selection="protein"):
    """
    Solvent Accessible Surface Area (SASA) 計算。

    Parameters:
        universe: mda.Universe
        selection: str — 原子選択
    """
    atoms = universe.select_atoms(selection)
    sasa = SASA(atoms)
    sasa.run()

    df = pd.DataFrame({
        "frame": range(len(sasa.results.area)),
        "sasa_A2": sasa.results.area,
    })

    print(f"SASA: mean={df['sasa_A2'].mean():.1f} Å², "
          f"std={df['sasa_A2'].std():.1f} Å²")
    return df
```

## 7. OpenFF 力場パラメータ化

```python
def parameterize_with_openff(smiles, force_field="openff-2.1.0"):
    """
    OpenFF Toolkit でリガンドの力場パラメータ自動生成。

    Parameters:
        smiles: str — リガンド SMILES
        force_field: str — OpenFF 力場バージョン

    K-Dense: openff
    """
    from openff.toolkit import Molecule, ForceField
    from openff.interchange import Interchange

    mol = Molecule.from_smiles(smiles)
    mol.generate_conformers(n_conformers=1)

    ff = ForceField(f"{force_field}.offxml")
    topology = mol.to_topology()
    interchange = Interchange.from_smirnoff(ff, topology)

    result = {
        "smiles": smiles,
        "force_field": force_field,
        "n_atoms": mol.n_atoms,
        "n_bonds": mol.n_bonds,
        "n_conformers": mol.n_conformers,
        "partial_charges": mol.partial_charges,
    }

    print(f"OpenFF parameterized: {smiles} ({mol.n_atoms} atoms, "
          f"FF={force_field})")
    return interchange, result
```

## 8. 統合 MD 解析パイプライン

```python
def md_analysis_pipeline(topology, trajectory, selection="protein"):
    """
    MD トラジェクトリ統合解析パイプライン。

    Pipeline:
        load → RMSD → RMSF → Rg → H-bonds → summary
    """
    u, info = load_trajectory(topology, trajectory)

    rmsd = compute_rmsd(u, "backbone")
    rmsf = compute_rmsf(u, "name CA")
    rg = compute_radius_of_gyration(u, selection)
    hbonds, hb_counts = hydrogen_bond_analysis(u)

    summary = {
        "system": info,
        "rmsd_mean_A": round(rmsd["rmsd_A"].mean(), 2),
        "rmsd_std_A": round(rmsd["rmsd_A"].std(), 2),
        "rmsf_mean_A": round(rmsf["rmsf_A"].mean(), 2),
        "n_flexible_residues": int(rmsf["flexible"].sum()),
        "rg_mean_A": round(rg["rg_A"].mean(), 2),
        "hbonds_mean_per_frame": round(hb_counts["n_hbonds"].mean(), 1),
    }

    print(f"\n=== MD Analysis Summary ===")
    for k, v in summary.items():
        if isinstance(v, dict):
            continue
        print(f"  {k}: {v}")

    return summary
```

---

## パイプライン統合

```
molecular-docking ──→ md-simulation ──→ admet-pharmacokinetics
  (ドッキングポーズ)    (MD 安定性評価)     (PK パラメータ推定)
        │                    │                      ↓
protein-structure ──┘        │              drug-target-profiling
  (PDB/AlphaFold)            ↓                (候補評価)
                    computational-materials
                    (pymatgen/VASP 連携)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/rmsd_timeseries.csv` | RMSD 時系列 | → publication-figures |
| `results/rmsf_per_residue.csv` | RMSF 残基別 | → protein-structure-analysis |
| `results/hbond_analysis.csv` | 水素結合解析 | → molecular-docking |
| `results/md_summary.json` | 統合サマリ | → admet-pharmacokinetics |

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `pdb` | PDB | 分子構造データベース参照 |

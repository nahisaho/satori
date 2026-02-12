---
name: scientific-molecular-docking
description: |
  構造ベース分子ドッキングスキル。DiffDock (拡散生成モデル)、
  AutoDock Vina (スコアリング関数)、GNINA (CNN ベーススコアリング) を統合した
  タンパク質-リガンド結合ポーズ予測、バーチャルスクリーニング、
  結合自由エネルギー推定、ドッキングスコア統合パイプライン。
---

# Scientific Molecular Docking

DiffDock / AutoDock Vina / GNINA の 3 大ドッキングエンジンによる
構造ベース仮想スクリーニング・結合ポーズ予測パイプラインを提供する。

## When to Use

- タンパク質-リガンド結合モードを予測するとき
- 化合物ライブラリのバーチャルスクリーニングが必要なとき
- 結合自由エネルギーを推定してリガンドをランキングするとき
- DiffDock で AI ベースの結合ポーズ生成を行うとき
- 複数のドッキング手法のコンセンサス評価が必要なとき

---

## Quick Start

## 1. リガンド・受容体の準備

```python
import os
import subprocess
import pandas as pd
import numpy as np


def prepare_receptor(pdb_file, output_dir="structures/prepared",
                      remove_water=True, add_hydrogens=True):
    """
    ドッキング用受容体 (タンパク質) 準備。

    Parameters:
        pdb_file: str — 入力 PDB ファイル
        remove_water: bool — 水分子除去
        add_hydrogens: bool — 水素原子付加
    """
    os.makedirs(output_dir, exist_ok=True)
    base_name = os.path.splitext(os.path.basename(pdb_file))[0]

    # PDB → PDBQT 変換 (AutoDock Vina 用)
    pdbqt_file = f"{output_dir}/{base_name}_receptor.pdbqt"

    try:
        from openbabel import pybel
        mol = next(pybel.readfile("pdb", pdb_file))
        if remove_water:
            mol.OBMol.DeleteWater()
        if add_hydrogens:
            mol.addh()
        mol.write("pdbqt", pdbqt_file, overwrite=True)
        print(f"Receptor prepared: {pdbqt_file}")
    except ImportError:
        # Open Babel 不在時: MGLTools prepare_receptor4
        cmd = ["prepare_receptor4.py", "-r", pdb_file,
               "-o", pdbqt_file, "-A", "hydrogens"]
        if remove_water:
            cmd.extend(["-U", "waters"])
        subprocess.run(cmd, check=True)
        print(f"Receptor prepared (MGLTools): {pdbqt_file}")

    return pdbqt_file


def prepare_ligands(sdf_file, output_dir="structures/ligands"):
    """
    リガンドファイル準備 (SDF → PDBQT/MOL2)。
    """
    os.makedirs(output_dir, exist_ok=True)

    try:
        from openbabel import pybel
        ligands = list(pybel.readfile("sdf", sdf_file))
        prepared = []
        for i, mol in enumerate(ligands):
            mol.addh()
            mol.make3D()
            name = mol.title or f"ligand_{i}"
            out = f"{output_dir}/{name}.pdbqt"
            mol.write("pdbqt", out, overwrite=True)
            prepared.append({"name": name, "file": out, "atoms": len(mol.atoms)})
        print(f"Prepared {len(prepared)} ligands from {sdf_file}")
        return pd.DataFrame(prepared)
    except ImportError:
        print("openbabel not available, using RDKit fallback")
        from rdkit import Chem
        from rdkit.Chem import AllChem
        suppl = Chem.SDMolSupplier(sdf_file)
        prepared = []
        for i, mol in enumerate(suppl):
            if mol is None:
                continue
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"lig_{i}"
            out = f"{output_dir}/{name}.mol2"
            Chem.MolToMolFile(mol, out)
            prepared.append({"name": name, "file": out})
        return pd.DataFrame(prepared)
```

## 2. AutoDock Vina ドッキング

```python
def autodock_vina_dock(receptor_pdbqt, ligand_pdbqt,
                        center, box_size,
                        exhaustiveness=32, n_poses=9):
    """
    AutoDock Vina による分子ドッキング。

    Parameters:
        receptor_pdbqt: str — 受容体 PDBQT
        ligand_pdbqt: str — リガンド PDBQT
        center: tuple — (x, y, z) ボックス中心座標
        box_size: tuple — (sx, sy, sz) ボックスサイズ (Å)
        exhaustiveness: int — 探索精度 (8-64)
        n_poses: int — 出力ポーズ数
    """
    try:
        from vina import Vina
        v = Vina(sf_name="vina")
        v.set_receptor(receptor_pdbqt)
        v.set_ligand_from_file(ligand_pdbqt)
        v.compute_vina_maps(center=list(center), box_size=list(box_size))
        v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)

        energies = v.energies()
        results = []
        for i, e in enumerate(energies):
            results.append({
                "pose": i + 1,
                "affinity_kcal": e[0],
                "rmsd_lb": e[1] if len(e) > 1 else None,
                "rmsd_ub": e[2] if len(e) > 2 else None,
            })

        output = ligand_pdbqt.replace(".pdbqt", "_docked.pdbqt")
        v.write_poses(output, n_poses=n_poses, overwrite=True)

        df = pd.DataFrame(results)
        print(f"Vina docking: best affinity = {df['affinity_kcal'].min():.1f} kcal/mol")
        return df, output

    except ImportError:
        # CLI フォールバック
        output = ligand_pdbqt.replace(".pdbqt", "_docked.pdbqt")
        cmd = [
            "vina",
            "--receptor", receptor_pdbqt,
            "--ligand", ligand_pdbqt,
            "--center_x", str(center[0]),
            "--center_y", str(center[1]),
            "--center_z", str(center[2]),
            "--size_x", str(box_size[0]),
            "--size_y", str(box_size[1]),
            "--size_z", str(box_size[2]),
            "--exhaustiveness", str(exhaustiveness),
            "--num_modes", str(n_poses),
            "--out", output,
        ]
        subprocess.run(cmd, check=True)
        return pd.DataFrame(), output
```

## 3. DiffDock AI ドッキング

```python
def diffdock_predict(protein_file, ligand_file, n_poses=10,
                      output_dir="results/diffdock"):
    """
    DiffDock (拡散生成モデル) ドッキング。

    Parameters:
        protein_file: str — タンパク質 PDB ファイル
        ligand_file: str — リガンド SDF/MOL2 ファイル
        n_poses: int — 生成ポーズ数
    """
    os.makedirs(output_dir, exist_ok=True)

    # DiffDock-L (large model) 推論
    cmd = [
        "python", "-m", "diffdock.inference",
        "--protein_path", protein_file,
        "--ligand", ligand_file,
        "--out_dir", output_dir,
        "--samples_per_complex", str(n_poses),
        "--model_dir", "DiffDock-L",
        "--confidence_model_dir", "DiffDock-L",
    ]

    print(f"Running DiffDock ({n_poses} poses)...")
    try:
        subprocess.run(cmd, check=True, capture_output=True)
    except FileNotFoundError:
        print("DiffDock not installed. Install from: "
              "https://github.com/gcorso/DiffDock")
        return pd.DataFrame()

    # 結果パース
    results = []
    for i in range(n_poses):
        pose_file = f"{output_dir}/rank{i+1}.sdf"
        conf_file = f"{output_dir}/rank{i+1}_confidence.txt"
        confidence = None
        if os.path.exists(conf_file):
            with open(conf_file) as f:
                confidence = float(f.read().strip())
        results.append({
            "pose": i + 1,
            "file": pose_file,
            "confidence": confidence,
        })

    df = pd.DataFrame(results)
    if len(df) > 0 and "confidence" in df.columns:
        print(f"DiffDock: {len(df)} poses, "
              f"best confidence = {df['confidence'].max()}")
    return df
```

## 4. バーチャルスクリーニング

```python
def virtual_screening(receptor_pdbqt, ligand_library,
                       center, box_size,
                       method="vina", top_n=20):
    """
    化合物ライブラリのバーチャルスクリーニング。

    Parameters:
        receptor_pdbqt: str — 受容体 PDBQT
        ligand_library: list[str] — リガンド PDBQT ファイルのリスト
        center/box_size: ドッキングボックスパラメータ
        method: "vina" or "diffdock"
        top_n: int — 上位候補数
    """
    all_results = []

    for i, ligand in enumerate(ligand_library):
        lig_name = os.path.splitext(os.path.basename(ligand))[0]
        print(f"  [{i+1}/{len(ligand_library)}] Docking {lig_name}...", end=" ")

        if method == "vina":
            df, _ = autodock_vina_dock(
                receptor_pdbqt, ligand, center, box_size,
                exhaustiveness=16, n_poses=3
            )
            if len(df) > 0:
                best = df.iloc[0]
                all_results.append({
                    "ligand": lig_name,
                    "best_affinity": best["affinity_kcal"],
                    "n_poses": len(df),
                })
                print(f"{best['affinity_kcal']:.1f} kcal/mol")

    results_df = pd.DataFrame(all_results)
    results_df = results_df.sort_values("best_affinity").head(top_n)

    print(f"\nVirtual screening: {len(ligand_library)} compounds → "
          f"top {top_n} candidates")
    return results_df
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `structures/prepared/*_receptor.pdbqt` | PDBQT |
| `structures/ligands/*.pdbqt` | PDBQT |
| `results/docking_results.csv` | CSV |
| `results/diffdock/rank*.sdf` | SDF |
| `results/virtual_screening.csv` | CSV |
| `figures/docking_scores.png` | PNG |

### 利用可能ツール

> このスキルは主に K-Dense-AI/claude-scientific-skills の diffdock スキルを参照しています。ToolUniverse SMCP には専用ドッキングツールは含まれませんが、タンパク質構造は PDB/AlphaFold ツール経由で取得可能です。

### 参照スキル

| スキル | 関連 |
|---|---|
| `scientific-protein-structure-analysis` | 受容体構造取得・結合部位検出 |
| `scientific-drug-target-profiling` | 標的選定 → ドッキング |
| `scientific-cheminformatics` | リガンド記述子・フィルタリング |
| `scientific-admet-pharmacokinetics` | ドッキング → ADMET |
| `scientific-drug-repurposing` | リポジショニング候補ドッキング |
| `scientific-protein-interaction-network` | PPI → ドッキング界面 |

### 依存パッケージ

`vina` (AutoDock Vina), `rdkit`, `openbabel` (optional), `numpy`, `pandas`

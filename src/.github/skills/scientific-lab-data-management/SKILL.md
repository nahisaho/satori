---
name: scientific-lab-data-management
description: |
  ラボデータ管理スキル。Benchling (ELN/DNA 設計/レジストリ)、
  DNAnexus (ゲノミクス PaaS)、LatchBio (ワークフロー)、
  OMERO (バイオイメージング)、Protocols.io (プロトコル共有)
  を統合したウェット・ドライラボデータ管理パイプライン。
---

# Scientific Lab Data Management

ウェットラボ実験管理からゲノミクスデータ処理まで、
ラボデータの生成・記録・解析・共有を統合管理するパイプライン。

## When to Use

- 電子実験ノート (ELN) でプロトコル・結果を記録するとき
- DNA 配列設計・クローニング計画を管理するとき
- ゲノミクス大規模データを PaaS 上で解析するとき
- バイオイメージングデータを構造化管理するとき
- 実験プロトコルを共有・再現するとき

---

## Quick Start

## 1. Benchling ELN / DNA 設計

```python
import json
import requests


class BenchlingClient:
    """
    Benchling API クライアント。

    Benchling 機能:
    - ELN (Electronic Lab Notebook): 実験記録
    - Molecular Biology: DNA 配列設計, プライマー設計, クローニング
    - Registry: サンプル・試薬レジストリ
    - Inventory: 在庫管理
    """

    def __init__(self, api_key, tenant_url):
        self.base_url = f"https://{tenant_url}/api/v2"
        self.headers = {
            "Authorization": f"Basic {api_key}",
            "Content-Type": "application/json",
        }

    def create_dna_sequence(self, name, bases, folder_id,
                              annotations=None):
        """
        DNA 配列の登録。

        Parameters:
        - name: 配列名
        - bases: 塩基配列 (ATCG)
        - folder_id: 保存先フォルダ
        - annotations: アノテーション [{name, start, end, type, strand}]
        """
        payload = {
            "name": name,
            "bases": bases,
            "folderId": folder_id,
            "isCircular": False,
            "annotations": annotations or [],
        }

        print(f"  Benchling DNA sequence: {name}")
        print(f"    Length: {len(bases)} bp")
        if annotations:
            print(f"    Annotations: {len(annotations)}")

        return payload

    def search_registry(self, query, schema_id=None, page_size=50):
        """
        Benchling Registry 検索。

        レジストリエンティティ:
        - プラスミド, 菌株, 抗体, 細胞株, 化合物
        """
        params = {
            "query": query,
            "pageSize": page_size,
        }
        if schema_id:
            params["schemaId"] = schema_id

        print(f"  Benchling registry search: '{query}'")

        return params

    def create_entry(self, name, folder_id, template_id=None):
        """
        ELN エントリ (実験ノート) 作成。
        """
        payload = {
            "name": name,
            "folderId": folder_id,
        }
        if template_id:
            payload["entryTemplateId"] = template_id

        print(f"  Benchling ELN entry: {name}")

        return payload
```

## 2. DNAnexus ゲノミクス PaaS

```python
import json


class DNAnexusClient:
    """
    DNAnexus Platform API クライアント。

    DNAnexus 機能:
    - データストレージ: FASTQ, BAM, VCF 等の大規模ファイル
    - ワークフロー実行: WDL/CWL/Applet ベース
    - コンプライアンス: HIPAA, GxP, FedRAMP
    - コラボレーション: プロジェクト単位のアクセス管理
    """

    def __init__(self, token):
        self.token = token
        self.base_url = "https://api.dnanexus.com"

    def upload_file(self, local_path, project_id, folder="/"):
        """
        ファイルアップロード。

        対応形式: FASTQ(.gz), BAM, CRAM, VCF, BED, etc.
        """
        print(f"  DNAnexus upload: {local_path}")
        print(f"    Project: {project_id}")
        print(f"    Destination: {folder}")

        return {"local_path": local_path, "project_id": project_id}

    def run_workflow(self, workflow_id, project_id, inputs):
        """
        ワークフロー実行。

        ワークフロー例:
        - GATK Best Practices (germline/somatic)
        - RNA-STAR alignment + featureCounts
        - DeepVariant caller
        - Structural variant calling
        """
        print(f"  DNAnexus workflow: {workflow_id}")
        print(f"    Project: {project_id}")
        print(f"    Inputs: {len(inputs)} parameters")

        return {
            "workflow_id": workflow_id,
            "project_id": project_id,
            "inputs": inputs,
        }

    def list_project_files(self, project_id, folder="/", name_glob=None):
        """
        プロジェクト内ファイル一覧。
        """
        params = {"folder": folder}
        if name_glob:
            params["name"] = {"glob": name_glob}

        print(f"  DNAnexus list: {project_id}{folder}")

        return params
```

## 3. OMERO バイオイメージング管理

```python
import json


class OMEROClient:
    """
    OMERO (Open Microscopy Environment Remote Objects) クライアント。

    OMERO 機能:
    - 画像データ管理: 150+ 画像フォーマット (Bio-Formats)
    - メタデータ: Key-Value, タグ, ROI
    - 解析統合: ImageJ/Fiji, CellProfiler, napari
    - アクセス制御: プロジェクト/グループ権限
    """

    def __init__(self, host, port=4064):
        self.host = host
        self.port = port

    def import_images(self, file_paths, dataset_id):
        """
        画像インポート。

        対応フォーマット (Bio-Formats):
        - OME-TIFF, ND2 (Nikon), CZI (Zeiss), LIF (Leica)
        - VSI (Olympus), SVS (Aperio), DICOM
        """
        print(f"  OMERO import: {len(file_paths)} images → Dataset {dataset_id}")

        return {"files": file_paths, "dataset_id": dataset_id}

    def create_roi(self, image_id, shapes):
        """
        ROI (Region of Interest) 作成。

        Shape タイプ:
        - Rectangle, Ellipse, Polygon
        - Line, Polyline, Point
        - Mask (binary mask)
        """
        print(f"  OMERO ROI: Image {image_id}, {len(shapes)} shapes")

        return {"image_id": image_id, "shapes": shapes}

    def query_images(self, project=None, dataset=None,
                      key_value_pairs=None):
        """
        画像検索 (メタデータベース)。

        フィルタ:
        - プロジェクト/データセット階層
        - Key-Value annotation
        - タグ
        - 取得日, 機器名
        """
        print(f"  OMERO query:")
        if project:
            print(f"    Project: {project}")
        if key_value_pairs:
            print(f"    Key-Value: {key_value_pairs}")

        return {"project": project, "dataset": dataset}
```

## 4. Protocols.io プロトコル共有

```python
import json


def create_protocol(title, description, steps, reagents=None,
                      doi_prefix="dx.doi.org/10.17504"):
    """
    Protocols.io プロトコル作成。

    Protocols.io:
    - DOI 付与による引用可能なプロトコル
    - バージョン管理
    - フォーク・改変・派生
    - JOVE, Nature Protocol Exchange 連携
    """
    protocol = {
        "title": title,
        "description": description,
        "steps": [],
        "reagents": reagents or [],
    }

    for i, step in enumerate(steps, 1):
        protocol["steps"].append({
            "step_number": i,
            "description": step.get("description", ""),
            "duration": step.get("duration"),
            "temperature": step.get("temperature"),
            "critical_step": step.get("critical", False),
            "expected_result": step.get("expected_result"),
        })

    print(f"  Protocol: {title}")
    print(f"    Steps: {len(steps)}")
    if reagents:
        print(f"    Reagents: {len(reagents)}")
    print(f"    DOI: {doi_prefix}/protocols.io...")

    return protocol


def fork_protocol(original_protocol_id, modifications):
    """
    既存プロトコルのフォークと改変。

    - 変更点の追跡
    - 元プロトコルへのリンク
    - バージョン番号の自動付与
    """
    print(f"  Forking protocol: {original_protocol_id}")
    print(f"    Modifications: {len(modifications)}")

    return {
        "forked_from": original_protocol_id,
        "modifications": modifications,
    }
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/benchling_sequences.json` | JSON |
| `results/benchling_registry.json` | JSON |
| `results/dnanexus_workflow_output.json` | JSON |
| `results/omero_image_metadata.json` | JSON |
| `results/protocol.json` | JSON |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

なし — 各プラットフォームの REST API を直接利用。

### 参照スキル

| スキル | 関連 |
|---|---|
| `scientific-bioinformatics` | ゲノミクスデータ解析 |
| `scientific-image-analysis` | 顯微鏡画像解析 |
| `scientific-gene-expression-transcriptomics` | RNA-seq データ管理 |
| `scientific-single-cell-genomics` | scRNA-seq データ管理 |
| `scientific-data-preprocessing` | データ前処理 |

### 依存パッケージ

`requests`, `json`, `pandas` (各プラットフォーム SDK: `benchling-sdk`, `dxpy`, `omero-py`)

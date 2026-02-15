---
name: scientific-spectral-signal
description: |
  分光スペクトルおよび生体信号の前処理・解析スキル。ベースライン補正、フィルタリング、
  ピーク検出、帯域パワー解析を行う際に使用。
  Scientific Skills Exp-08（ECG/EEG）、Exp-11（ラマン分光）で確立したパターン。
tu_tools:
  - key: biotools
    name: bio.tools
    description: スペクトル信号処理ツール検索
---

# Scientific Spectral & Signal Processing

分光スペクトル（ラマン、IR、UV-Vis など）および生体信号（ECG, EEG）の
前処理・解析パイプラインスキル。ノイズ除去、ベースライン補正、ピーク検出、
周波数解析の標準ワークフローを提供する。

## When to Use

- スペクトルデータの前処理（ベースライン補正、平滑化、正規化）
- 時系列信号のフィルタリング（バンドパス、ノッチ）
- ピーク検出と定量化
- 周波数帯域パワー解析（EEG δ/θ/α/β/γ）
- スペクトル類似度解析と分類

## Quick Start

## スペクトル前処理パイプライン

### 1. ALS ベースライン補正（Exp-11）

```python
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

def baseline_correction_als(y, lam=1e5, p=0.01, niter=10):
    """
    Asymmetric Least Squares (ALS) によるベースライン推定。
    蛍光バックグラウンドの除去に使用。
    lam: 平滑化パラメータ（大きいほど滑らか）
    p: 非対称性パラメータ（小さいほどベースライン追従が強い）
    """
    L = len(y)
    D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L - 2))
    w = np.ones(L)
    for _ in range(niter):
        W = sparse.diags(w)
        Z = W + lam * D.dot(D.transpose())
        z = spsolve(Z, w * y)
        w = p * (y > z) + (1 - p) * (y < z)
    return z

def remove_baseline(spectrum, wavenumbers, lam=1e5, p=0.01):
    """ベースライン補正済みスペクトルを返す。"""
    baseline = baseline_correction_als(spectrum, lam=lam, p=p)
    return spectrum - baseline
```

### 2. Savitzky-Golay 平滑化

```python
from scipy.signal import savgol_filter

def smooth_spectrum(spectrum, window_length=11, polyorder=3):
    """Savitzky-Golay フィルタでスペクトルを平滑化する。"""
    return savgol_filter(spectrum, window_length=window_length, polyorder=polyorder)
```

### 3. 正規化

```python
def normalize_spectrum(spectrum, method="minmax"):
    """
    スペクトルを正規化する。
    method: 'minmax', 'snv' (Standard Normal Variate), 'area'
    """
    if method == "minmax":
        return (spectrum - spectrum.min()) / (spectrum.max() - spectrum.min() + 1e-10)
    elif method == "snv":
        return (spectrum - spectrum.mean()) / (spectrum.std() + 1e-10)
    elif method == "area":
        return spectrum / (np.trapz(np.abs(spectrum)) + 1e-10)
    else:
        raise ValueError(f"Unknown method: {method}")
```

### 4. ピーク検出

```python
from scipy.signal import find_peaks

def detect_peaks(spectrum, wavenumbers=None, height=None,
                 prominence=0.05, distance=10, width=None):
    """
    スペクトル中のピークを検出する。
    返値: ピーク位置インデックス、ピーク属性辞書
    """
    peaks, properties = find_peaks(
        spectrum, height=height, prominence=prominence,
        distance=distance, width=width
    )

    if wavenumbers is not None:
        peak_positions = wavenumbers[peaks]
    else:
        peak_positions = peaks

    return peaks, peak_positions, properties
```

## 生体信号処理パイプライン

### 5. バンドパス/ノッチフィルタ（Exp-08）

```python
from scipy.signal import butter, sosfilt, iirnotch

def bandpass_filter(signal, fs, lowcut, highcut, order=4):
    """Butterworth バンドパスフィルタ。"""
    nyq = 0.5 * fs
    sos = butter(order, [lowcut / nyq, highcut / nyq], btype="band", output="sos")
    return sosfilt(sos, signal)

def notch_filter(signal, fs, freq=50.0, Q=30.0):
    """商用電源ノイズ除去用ノッチフィルタ（50/60 Hz）。"""
    b, a = iirnotch(freq, Q, fs)
    from scipy.signal import filtfilt
    return filtfilt(b, a, signal)
```

### 6. 周波数帯域パワー解析（EEG）

```python
from scipy.signal import welch

EEG_BANDS = {
    "delta": (0.5, 4),
    "theta": (4, 8),
    "alpha": (8, 13),
    "beta":  (13, 30),
    "gamma": (30, 100),
}

def band_power(signal, fs, band, method="welch"):
    """指定周波数帯域のパワーを算出する。"""
    freqs, psd = welch(signal, fs=fs, nperseg=min(len(signal), 256))
    low, high = band
    idx = np.logical_and(freqs >= low, freqs <= high)
    return np.trapz(psd[idx], freqs[idx])

def eeg_band_powers(signal, fs, bands=None):
    """EEG の全帯域パワーを一括計算する。"""
    if bands is None:
        bands = EEG_BANDS
    powers = {}
    for name, (low, high) in bands.items():
        powers[name] = band_power(signal, fs, (low, high))
    total = sum(powers.values())
    relative = {k: v / total for k, v in powers.items()}
    return powers, relative
```

### 7. R 波検出 & HRV 解析（ECG）

```python
def detect_r_peaks(ecg_signal, fs, height_factor=0.5, distance_ms=300):
    """ECG 信号から R 波を検出する。"""
    min_distance = int(distance_ms * fs / 1000)
    threshold = height_factor * np.max(ecg_signal)
    peaks, _ = find_peaks(ecg_signal, height=threshold, distance=min_distance)
    return peaks

def hrv_analysis(r_peaks, fs):
    """R-R 間隔から HRV 指標を算出する。"""
    rr_intervals = np.diff(r_peaks) / fs * 1000  # ms
    return {
        "mean_RR": np.mean(rr_intervals),
        "SDNN": np.std(rr_intervals, ddof=1),
        "RMSSD": np.sqrt(np.mean(np.diff(rr_intervals) ** 2)),
        "pNN50": np.sum(np.abs(np.diff(rr_intervals)) > 50) / len(rr_intervals) * 100,
        "mean_HR": 60000 / np.mean(rr_intervals),
    }
```

### 8. スペクトル類似度解析（Exp-11）

```python
from scipy.spatial.distance import cosine
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

def spectral_similarity_matrix(spectra_dict, method="cosine"):
    """
    スペクトル間の類似度行列を算出する。
    spectra_dict: {name: spectrum_array}
    """
    names = list(spectra_dict.keys())
    n = len(names)
    sim_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            if method == "cosine":
                sim_matrix[i, j] = 1 - cosine(spectra_dict[names[i]],
                                               spectra_dict[names[j]])
            elif method == "pearson":
                sim_matrix[i, j] = np.corrcoef(spectra_dict[names[i]],
                                                spectra_dict[names[j]])[0, 1]

    return pd.DataFrame(sim_matrix, index=names, columns=names)
```

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `biotools` | bio.tools | スペクトル信号処理ツール検索 |

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/peak_detection_results.csv` | CSV |
| `results/hrv_metrics.csv` | CSV |
| `results/spectral_similarity.csv` | CSV |
| `figures/spectrum_processed.png` | PNG |
| `figures/peak_detection.png` | PNG |
| `figures/eeg_band_powers.png` | PNG |

#### 参照実験

- **Exp-08**: ECG/EEG 生体信号処理（フィルタ、R 波検出、HRV、帯域パワー）
- **Exp-11**: ラマン分光（ALS ベースライン、ピーク検出、類似度行列、クラスタリング）

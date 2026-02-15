---
name: scientific-quantum-computing
description: |
  量子計算スキル。Qiskit・Cirq・PennyLane・QuTiP を活用し、
  量子回路設計・シミュレーション・変分量子アルゴリズム（VQE/QAOA）・
  量子化学計算・量子機械学習を支援。
  「量子回路を設計して」「VQE で基底エネルギーを求めて」「量子シミュレーションして」で発火。
tu_tools:
  - key: papers_with_code
    name: Papers with Code
    description: 量子計算論文・ベンチマーク検索
---

# Scientific Quantum Computing

量子計算のための解析・シミュレーションスキル。
主要量子フレームワーク（Qiskit, Cirq, PennyLane, QuTiP）を活用し、
量子アルゴリズム設計から実機実行までを支援する。

## When to Use

- 量子回路の設計・シミュレーション
- 変分量子固有値ソルバー (VQE) による分子エネルギー計算
- QAOA (Quantum Approximate Optimization Algorithm) による組合せ最適化
- 量子機械学習 (QML) のパラメトリック回路設計
- 量子ノイズモデリング・エラー緩和
- 量子化学計算（分子ハミルトニアン構築）

## Quick Start

### 量子計算パイプライン

```
Phase 1: Problem Formulation
  - 問題のハミルトニアン / コスト関数定義
  - 量子ビット数・エンコーディング設計
  - 古典/量子ハイブリッド戦略選定
    ↓
Phase 2: Circuit Design
  - Ansatz 選択 (EfficientSU2, UCCSD, HEA, etc.)
  - パラメトリック回路構築
  - ゲート分解・トランスパイル
    ↓
Phase 3: Simulation / Execution
  - ステートベクトルシミュレータ
  - ノイズモデル付きシミュレーション
  - 実機バックエンド (IBM Quantum, Google QCS)
    ↓
Phase 4: Optimization
  - 古典オプティマイザ (COBYLA, SPSA, L-BFGS-B)
  - Gradient 計算 (Parameter Shift Rule)
  - 収束判定
    ↓
Phase 5: Error Mitigation
  - Zero Noise Extrapolation (ZNE)
  - Probabilistic Error Cancellation (PEC)
  - Measurement Error Mitigation
    ↓
Phase 6: Analysis & Reporting
  - エネルギー収束プロット
  - 量子状態トモグラフィー
  - 結果の古典ベンチマーク比較
```

## Workflow

### 1. Qiskit: 量子回路設計

```python
from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from qiskit.visualization import plot_histogram, plot_bloch_multivector

# === 基本量子回路 ===
qc = QuantumCircuit(2, 2)

# Bell State 生成
qc.h(0)           # Hadamard on qubit 0
qc.cx(0, 1)       # CNOT: qubit 0 → qubit 1
qc.measure([0, 1], [0, 1])

# シミュレーション
simulator = AerSimulator()
compiled = transpile(qc, simulator)
result = simulator.run(compiled, shots=10000).result()
counts = result.get_counts()
print(f"Bell state measurement: {counts}")  # {'00': ~5000, '11': ~5000}

# 回路可視化
print(qc.draw("text"))
```

### 2. VQE (Variational Quantum Eigensolver)

```python
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_algorithms import VQE
from qiskit_algorithms.optimizers import COBYLA, SPSA
from qiskit.circuit.library import EfficientSU2
from qiskit.primitives import Estimator

# === H2 分子の基底エネルギー ===
# Step 1: 分子ハミルトニアン構築
driver = PySCFDriver(atom="H 0 0 0; H 0 0 0.735", basis="sto-3g")
problem = driver.run()

# Step 2: Jordan-Wigner マッピング (フェルミオン → qubit)
mapper = JordanWignerMapper()
qubit_op = mapper.map(problem.second_q_ops()[0])
print(f"Qubit Hamiltonian terms: {len(qubit_op)}")
print(f"Number of qubits: {qubit_op.num_qubits}")

# Step 3: Ansatz (パラメトリック回路)
ansatz = EfficientSU2(num_qubits=qubit_op.num_qubits, reps=2,
                       entanglement="linear")
print(f"Ansatz parameters: {ansatz.num_parameters}")

# Step 4: VQE 実行
estimator = Estimator()
optimizer = COBYLA(maxiter=500)
vqe = VQE(estimator=estimator, ansatz=ansatz, optimizer=optimizer)

result = vqe.compute_minimum_eigenvalue(qubit_op)
print(f"VQE ground state energy: {result.eigenvalue:.6f} Ha")
print(f"Exact (FCI): -1.137276 Ha")
print(f"Chemical accuracy (1.6 mHa): {abs(result.eigenvalue - (-1.137276)) < 0.0016}")
```

### 3. PennyLane: 量子機械学習

```python
import pennylane as qml
from pennylane import numpy as pnp

# === 量子カーネル分類器 ===
n_qubits = 4
dev = qml.device("default.qubit", wires=n_qubits)

@qml.qnode(dev)
def quantum_kernel_circuit(x1, x2):
    """Quantum kernel: |<φ(x1)|φ(x2)>|²"""
    # Feature map for x1
    for i in range(n_qubits):
        qml.RY(x1[i], wires=i)
    for i in range(n_qubits - 1):
        qml.CNOT(wires=[i, i + 1])
    for i in range(n_qubits):
        qml.RZ(x1[i] ** 2, wires=i)

    # Adjoint feature map for x2
    qml.adjoint(qml.RZ)(x2[0] ** 2, wires=0)
    for i in reversed(range(n_qubits - 1)):
        qml.CNOT(wires=[i, i + 1])
    for i in reversed(range(n_qubits)):
        qml.adjoint(qml.RY)(x2[i], wires=i)

    return qml.probs(wires=range(n_qubits))


# === Variational Quantum Classifier ===
@qml.qnode(dev)
def variational_classifier(weights, x):
    """パラメトリック量子分類回路"""
    # Data encoding (Angle Encoding)
    for i in range(n_qubits):
        qml.RX(x[i], wires=i)

    # Variational layers
    n_layers = weights.shape[0]
    for layer in range(n_layers):
        for i in range(n_qubits):
            qml.RY(weights[layer, i, 0], wires=i)
            qml.RZ(weights[layer, i, 1], wires=i)
        for i in range(n_qubits - 1):
            qml.CNOT(wires=[i, i + 1])

    return qml.expval(qml.PauliZ(0))


def train_vqc(X_train, y_train, n_layers=3, epochs=100, lr=0.01):
    """VQC トレーニング"""
    weights = pnp.random.randn(n_layers, n_qubits, 2, requires_grad=True)
    bias = pnp.array(0.0, requires_grad=True)
    opt = qml.AdamOptimizer(stepsize=lr)

    for epoch in range(epochs):
        def cost_fn(w, b):
            predictions = [variational_classifier(w, x) + b for x in X_train]
            loss = pnp.mean((pnp.array(predictions) - y_train) ** 2)
            return loss

        weights, bias = opt.step(cost_fn, weights, bias)

        if (epoch + 1) % 20 == 0:
            loss = cost_fn(weights, bias)
            print(f"Epoch {epoch+1}: Loss = {loss:.4f}")

    return weights, bias
```

### 4. Cirq: Google Quantum AI

```python
import cirq
import numpy as np

# === Cirq 基本回路 ===
qubits = cirq.LineQubit.range(3)

circuit = cirq.Circuit([
    cirq.H(qubits[0]),
    cirq.CNOT(qubits[0], qubits[1]),
    cirq.CNOT(qubits[1], qubits[2]),
    cirq.measure(*qubits, key="result"),
])
print(circuit)

# ノイズモデル付きシミュレーション
noise_model = cirq.ConstantQubitNoiseModel(
    qubit_noise_gate=cirq.DepolarizingChannel(p=0.01)
)
noisy_simulator = cirq.DensityMatrixSimulator(noise=noise_model)
result = noisy_simulator.simulate(circuit[:-1])  # without measurement
print(f"State fidelity with noise: {cirq.fidelity(result.final_density_matrix, np.array([[0.5, 0, 0, 0.5], [0, 0, 0, 0], [0, 0, 0, 0], [0.5, 0, 0, 0.5]])):.4f}")
```

### 5. QAOA (Quantum Approximate Optimization Algorithm)

```python
from qiskit_algorithms import QAOA
from qiskit_algorithms.optimizers import COBYLA
from qiskit.quantum_info import SparsePauliOp
from qiskit.primitives import Sampler

def solve_maxcut_qaoa(adjacency_matrix, p=2):
    """
    QAOA で MaxCut 問題を解く

    C(z) = Σ_{(i,j)∈E} (1 - z_i z_j) / 2
    """
    n = len(adjacency_matrix)

    # コスト Hamilton 構築
    pauli_list = []
    for i in range(n):
        for j in range(i + 1, n):
            if adjacency_matrix[i][j] != 0:
                # 0.5 * (I - Z_i Z_j)
                z_str = ["I"] * n
                z_str[i] = "Z"
                z_str[j] = "Z"
                pauli_list.append(("".join(z_str), -0.5 * adjacency_matrix[i][j]))
                pauli_list.append(("I" * n, 0.5 * adjacency_matrix[i][j]))

    cost_op = SparsePauliOp.from_list(pauli_list).simplify()

    # QAOA
    sampler = Sampler()
    optimizer = COBYLA(maxiter=300)
    qaoa = QAOA(sampler=sampler, optimizer=optimizer, reps=p)

    result = qaoa.compute_minimum_eigenvalue(cost_op)
    return result
```

### 6. QuTiP: 量子系ダイナミクス

```python
import qutip

# === 量子系の時間発展 ===
# Jaynes-Cummings モデル (光-物質相互作用)
N = 20  # Fock 状態数
wc = 1.0  # Cavity frequency
wa = 1.0  # Atom frequency
g = 0.05  # Coupling strength

# ハミルトニアン
a = qutip.tensor(qutip.destroy(N), qutip.qeye(2))
sm = qutip.tensor(qutip.qeye(N), qutip.destroy(2))
H = wc * a.dag() * a + wa * sm.dag() * sm + g * (a.dag() * sm + a * sm.dag())

# 初期状態: |1 photon, excited atom>
psi0 = qutip.tensor(qutip.basis(N, 1), qutip.basis(2, 0))

# 時間発展
tlist = np.linspace(0, 25, 1000)
result = qutip.mesolve(H, psi0, tlist, [], [a.dag() * a, sm.dag() * sm])

# Rabi 振動プロット
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(tlist, result.expect[0], label="Cavity photon number")
ax.plot(tlist, result.expect[1], label="Atom excitation")
ax.set_xlabel("Time")
ax.set_ylabel("Expectation value")
ax.set_title("Jaynes-Cummings Model: Rabi Oscillations")
ax.legend()
plt.savefig("figures/quantum_rabi.png", dpi=300, bbox_inches="tight")
plt.show()
```

---

## Best Practices

1. **Ansatz 深さとノイズのトレードオフ**: 深い回路は表現力が高いがノイズに弱い (p ≤ 5 推奨)
2. **Parameter Shift Rule**: ゲート勾配の正確な解析的計算に使用
3. **Barren Plateau 回避**: ランダム初期化を避け、問題特化 ansatz を選択
4. **エラー緩和必須**: 実機結果は常に ZNE / PEC で補正
5. **古典ベンチマーク**: VQE 結果は必ず FCI / CCSD(T) と比較
6. **ノイズモデルの検証**: 実機ノイズキャリブレーションデータを使用

## Completeness Checklist

- [ ] 問題のハミルトニアン / コスト関数定義
- [ ] 量子ビットエンコーディング・ansatz 設計
- [ ] シミュレータによる検証 (ステートベクトル)
- [ ] ノイズモデル付きシミュレーション
- [ ] オプティマイザ収束確認
- [ ] 古典手法とのベンチマーク比較
- [ ] エラー緩和手法適用（実機の場合）
- [ ] 結果レポート・収束プロット生成

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `papers_with_code` | Papers with Code | 量子計算論文・ベンチマーク検索 |

## References

### Output Files

| ファイル | 形式 | 生成タイミング |
|---|---|---|
| `results/quantum_result.json` | 量子計算結果（JSON） | 計算完了時 |
| `figures/quantum_convergence.png` | 収束プロット | VQE/QAOA 完了時 |
| `figures/quantum_rabi.png` | 量子ダイナミクス図 | QuTiP シミュレーション時 |

### 参照スキル

| スキル | 連携 |
|---|---|
| `scientific-cheminformatics` | ← 分子構造 → ハミルトニアン |
| `scientific-process-optimization` | ← QAOA による組合せ最適化 |
| `scientific-statistical-testing` | ← 量子測定結果の統計解析 |
| `scientific-deep-learning` | → 量子-古典ハイブリッド ML |
| `scientific-bayesian-statistics` | ← ベイズ最適化による VQE パラメータ最適化 |

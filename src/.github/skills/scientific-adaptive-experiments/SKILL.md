---
name: scientific-adaptive-experiments
description: |
  適応的実験計画スキル。多腕バンディット (Thompson Sampling/UCB)・
  ベイズ適応設計・逐次検定 (SPRT)・
  Response-Adaptive Randomization・早期停止規則。
---

# Scientific Adaptive Experiments

実験中のデータに基づいて動的に実験計画を修正する
適応的実験設計パイプラインを提供する。

## When to Use

- A/B テストの適応的割り当てを行うとき
- 多腕バンディットで探索と活用を最適化するとき
- 逐次解析で早期停止判定を行うとき
- ベイズ適応設計で用量探索するとき
- Response-Adaptive Randomization で臨床試験を設計するとき

---

## Quick Start

## 1. 多腕バンディット

```python
import numpy as np
import pandas as pd
from typing import List, Dict


class ThompsonSamplingBandit:
    """
    Thompson Sampling ベータ-ベルヌーイバンディット。
    """

    def __init__(self, n_arms, prior_alpha=1.0, prior_beta=1.0):
        """
        Parameters:
            n_arms: int — アーム数
            prior_alpha: float — Beta 事前分布の α
            prior_beta: float — Beta 事前分布の β
        """
        self.n_arms = n_arms
        self.alphas = np.full(n_arms, prior_alpha)
        self.betas = np.full(n_arms, prior_beta)
        self.history = []

    def select_arm(self):
        samples = np.array([
            np.random.beta(a, b)
            for a, b in zip(self.alphas, self.betas)])
        return int(np.argmax(samples))

    def update(self, arm, reward):
        if reward > 0:
            self.alphas[arm] += 1
        else:
            self.betas[arm] += 1
        self.history.append({"arm": arm, "reward": reward})

    def get_summary(self):
        records = []
        for i in range(self.n_arms):
            a, b = self.alphas[i], self.betas[i]
            records.append({
                "arm": i,
                "alpha": a, "beta": b,
                "mean": a / (a + b),
                "n_pulls": int(a + b - 2),
                "95%_lower": float(np.percentile(
                    np.random.beta(a, b, 10000), 2.5)),
            })
        return pd.DataFrame(records)


def run_bandit_experiment(true_rates, n_rounds=1000,
                          algorithm="thompson", seed=42):
    """
    バンディット実験シミュレーション。

    Parameters:
        true_rates: list[float] — 各アームの真の成功率
        n_rounds: int — ラウンド数
        algorithm: str — "thompson" / "ucb" / "epsilon_greedy"
        seed: int — 乱数シード
    """
    rng = np.random.default_rng(seed)
    n_arms = len(true_rates)

    if algorithm == "thompson":
        bandit = ThompsonSamplingBandit(n_arms)
        for t in range(n_rounds):
            arm = bandit.select_arm()
            reward = int(rng.random() < true_rates[arm])
            bandit.update(arm, reward)
    elif algorithm == "ucb":
        counts = np.zeros(n_arms)
        values = np.zeros(n_arms)
        history = []
        for t in range(n_rounds):
            if t < n_arms:
                arm = t
            else:
                ucb = values + np.sqrt(2 * np.log(t) / (counts + 1e-10))
                arm = int(np.argmax(ucb))
            reward = int(rng.random() < true_rates[arm])
            counts[arm] += 1
            values[arm] += (reward - values[arm]) / counts[arm]
            history.append({"arm": arm, "reward": reward})
        bandit = type("UCB", (), {
            "history": history,
            "get_summary": lambda self: pd.DataFrame([
                {"arm": i, "n_pulls": int(counts[i]),
                 "mean": values[i]} for i in range(n_arms)])
        })()
    else:  # epsilon_greedy
        epsilon = 0.1
        counts = np.zeros(n_arms)
        values = np.zeros(n_arms)
        history = []
        for t in range(n_rounds):
            if rng.random() < epsilon:
                arm = rng.integers(n_arms)
            else:
                arm = int(np.argmax(values))
            reward = int(rng.random() < true_rates[arm])
            counts[arm] += 1
            values[arm] += (reward - values[arm]) / counts[arm]
            history.append({"arm": arm, "reward": reward})
        bandit = type("EG", (), {
            "history": history,
            "get_summary": lambda self: pd.DataFrame([
                {"arm": i, "n_pulls": int(counts[i]),
                 "mean": values[i]} for i in range(n_arms)])
        })()

    summary = bandit.get_summary()
    best_arm = summary.loc[summary["mean"].idxmax(), "arm"]
    print(f"Bandit ({algorithm}, {n_rounds} rounds):")
    print(f"  Best arm: {best_arm} (est. rate={summary.loc[summary['arm']==best_arm, 'mean'].values[0]:.3f})")
    print(f"  True best: {np.argmax(true_rates)} (rate={max(true_rates):.3f})")
    return bandit, summary
```

## 2. 逐次検定 (SPRT)

```python
def sequential_probability_ratio_test(data_stream,
                                       h0_rate=0.5,
                                       h1_rate=0.6,
                                       alpha=0.05, beta=0.2):
    """
    Wald の逐次確率比検定 (SPRT)。

    Parameters:
        data_stream: iterable — 逐次観測データ (0/1)
        h0_rate: float — 帰無仮説下の成功率
        h1_rate: float — 対立仮説下の成功率
        alpha: float — 第 1 種の過誤確率
        beta: float — 第 2 種の過誤確率
    """
    A = np.log((1 - beta) / alpha)  # H1 採択境界
    B = np.log(beta / (1 - alpha))  # H0 採択境界

    log_lr = 0.0
    history = []

    for i, x in enumerate(data_stream):
        # 対数尤度比の更新
        if x == 1:
            log_lr += np.log(h1_rate / h0_rate)
        else:
            log_lr += np.log((1 - h1_rate) / (1 - h0_rate))

        decision = None
        if log_lr >= A:
            decision = "reject_H0"
        elif log_lr <= B:
            decision = "accept_H0"

        history.append({
            "step": i + 1, "observation": x,
            "log_lr": log_lr,
            "upper_bound": A, "lower_bound": B,
            "decision": decision,
        })

        if decision is not None:
            print(f"SPRT: {decision} at step {i+1} "
                  f"(log_LR={log_lr:.3f})")
            break

    df = pd.DataFrame(history)
    if df["decision"].iloc[-1] is None:
        print(f"SPRT: Inconclusive after {len(df)} observations")
    return df
```

## 3. ベイズ適応用量探索

```python
def bayesian_adaptive_dose_finding(dose_levels, n_patients=30,
                                    target_toxicity=0.33,
                                    prior_alpha=1.0,
                                    prior_beta=1.0, seed=42):
    """
    ベイズ適応用量探索 (CRM 簡易版)。

    Parameters:
        dose_levels: list[float] — 用量レベル
        n_patients: int — 患者数
        target_toxicity: float — 目標毒性率
        prior_alpha: float — Beta 事前分布 α
        prior_beta: float — Beta 事前分布 β
        seed: int — 乱数シード
    """
    rng = np.random.default_rng(seed)
    n_doses = len(dose_levels)
    alphas = np.full(n_doses, prior_alpha)
    betas = np.full(n_doses, prior_beta)
    history = []

    # 真の毒性率 (シミュレーション用)
    true_tox = np.linspace(0.05, 0.60, n_doses)

    current_dose = 0
    for patient in range(n_patients):
        toxicity = int(rng.random() < true_tox[current_dose])

        if toxicity:
            alphas[current_dose] += 1
        else:
            betas[current_dose] += 1

        history.append({
            "patient": patient + 1,
            "dose_level": current_dose,
            "dose_value": dose_levels[current_dose],
            "toxicity": toxicity,
        })

        # 次の用量選択 (最も target に近い推定毒性率)
        means = alphas / (alphas + betas)
        distances = np.abs(means - target_toxicity)
        current_dose = int(np.argmin(distances))

    # MTD 推定
    final_means = alphas / (alphas + betas)
    mtd_idx = int(np.argmin(np.abs(final_means - target_toxicity)))

    summary = pd.DataFrame({
        "dose_level": range(n_doses),
        "dose_value": dose_levels,
        "est_toxicity": final_means,
        "n_treated": [sum(1 for h in history
                         if h["dose_level"] == i)
                      for i in range(n_doses)],
    })

    print(f"Bayesian adaptive: MTD = dose {mtd_idx} "
          f"(value={dose_levels[mtd_idx]}, "
          f"est_tox={final_means[mtd_idx]:.3f})")
    return summary, pd.DataFrame(history)
```

---

## パイプライン統合

```
[実験目的] → adaptive-experiments → statistical-testing
               (適応設計)              (最終解析)
                    │
              doe ← statistical-simulation
          (古典計画)    (検出力分析)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `bandit_summary.csv` | バンディット結果 | → 最適アーム |
| `sprt_history.csv` | SPRT 検定履歴 | → 判定結果 |
| `dose_finding.csv` | 用量探索結果 | → MTD 推定 |

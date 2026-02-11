# IEEE Transactions テンプレート

> 対象: IEEE Transactions on XXX, IEEE Access, IEEE Journal of XXX 等
> 構成: Abstract → Index Terms → I. Introduction → II–IV. Body → V. Conclusion
> 特徴: ローマ数字セクション番号。TABLE I (大文字ローマ数字)。Fig. 1(a) 形式。
> 語数目安: Regular Paper 6,000–8,000 words / Short Paper/Letter ≤ 4,000 words

---

# [タイトル — Title Case。略語はフルスペル後に括弧で定義]

Author A, *Member, IEEE*, Author B, *Senior Member, IEEE*,
and Author C, *Fellow, IEEE*

---

## Abstract

<!-- 最大 250 words。段落分けなし。一人称を避ける。 -->

[背景と課題を述べる。]
[提案手法を述べる。]
[実験結果を定量的に示す。]
[結論と意義を述べる。]

---

## Index Terms

<!-- IEEE Thesaurus から選択。アルファベット順。各語を大文字開始。 -->

Index Term 1, Index Term 2, Index Term 3, Index Term 4.

---

## I. Introduction

<!-- 参照は角括弧 [1]。本文と括弧の間にスペースなし。 -->

[研究分野の背景と重要性] [1], [2].

[先行研究のレビュー] [3]–[5].

[先行研究の課題・限界を述べる].

[本論文の貢献を明示する]:
1) [貢献 1]
2) [貢献 2]
3) [貢献 3]

The remainder of this paper is organized as follows.
Section II describes [内容]. Section III presents [内容].
Section IV discusses [内容]. Finally, Section V concludes the paper.

---

## II. [Related Work / System Model / Problem Formulation]

### *A. [Subsection タイトル]*

[関連研究のレビューまたはシステムモデルの記述]

### *B. [Subsection タイトル]*

[問題定式化または理論的背景]

The objective function can be expressed as
$$
\min_{x} f(x) = \sum_{i=1}^{N} g_i(x_i)
$$
subject to $h(x) \leq 0$.

---

## III. [Proposed Method / Methodology]

### *A. [Subsection タイトル]*

[提案手法の詳細な説明]

**Algorithm 1** [アルゴリズム名]

```
Input: [入力パラメータ]
Output: [出力]
1: Initialize parameters
2: for t = 1 to T do
3:    Compute [ステップの説明]
4:    Update [パラメータの更新]
5: end for
6: return [結果]
```

### *B. [Subsection タイトル]*

[手法の詳細の続き]

---

## IV. [Experimental Results / Performance Evaluation]

### *A. Experimental Setup*

[データセット、ハードウェア、ソフトウェア環境の詳細]

The experiments were conducted on [hardware spec].
The proposed method was implemented in [language/framework].

### *B. [実験 1 のタイトル]*

[結果の記述。Fig. 1(a) を参照。]

As shown in Fig. 1, the proposed method achieves [結果].
TABLE I summarizes the quantitative comparison.

### *C. [実験 2 のタイトル]*

[結果の記述。Fig. 2 を参照。]

### *D. Computational Complexity*

[計算量の分析。O(n log n) 等の表記を含む。]

---

## V. Conclusion

[主要な貢献を再掲する。]
[結果のハイライトを述べる。]
[今後の課題を 1-2 文で述べる。]

---

## Appendix

<!-- 必要に応じて。定理の証明、追加の数式展開等。 -->

### A. [Proof of Theorem 1]

[証明の詳細]

---

## Acknowledgment

<!-- "Acknowledgments" ではなく単数形 "Acknowledgment" -->

The authors would like to thank [謝辞の対象] for [貢献内容].

---

## References

<!-- 角括弧番号順。著者名はイニシャル先行。ジャーナル名はイタリック。 -->

[1] A. B. Author, C. D. Author, and E. F. Author, "Title of article,"
    *IEEE Trans. XXX*, vol. XX, no. Y, pp. 123–135, Mon. Year.

[2] A. B. Author and C. D. Author, "Title of conference paper," in
    *Proc. IEEE Conf. XXX*, City, Country, Year, pp. 100–105.

[3] A. B. Author, *Title of Book*, Xth ed. City, State: Publisher, Year.

---

## Figure Captions

**Fig. 1.** [タイトル。] (a) [パネル a の説明。] (b) [パネル b の説明。]

**Fig. 2.** [タイトル。] [説明]

---

## Tables

**TABLE I**
[タイトル — 大文字]

| Method | Metric A | Metric B | Metric C |
|---|---|---|---|
| Baseline [3] | XX.X | YY.Y | ZZ.Z |
| Method A [5] | XX.X | YY.Y | ZZ.Z |
| **Proposed** | **XX.X** | **YY.Y** | **ZZ.Z** |

<!-- ベスト値は太字で強調 -->

---

## Author Biography

<!-- IEEE Transactions では著者略歴が必要 -->

**Author A** (M'XX–SM'YY) received the B.S. degree in [分野] from
[大学], [国], in [年], and the Ph.D. degree from [大学] in [年].
[現在の所属と研究分野を記述。]

**Author B** [略歴]

---

## IEEE 数式表記規約

```
- スカラー: イタリック (x, y, α)
- ベクトル: 太字小文字 (𝐱, 𝐲)
- 行列: 太字大文字 (𝐀, 𝐁)
- 集合: カリグラフィック (𝒜, ℬ)
- 確率: Pr{·} or P(·)
- 期待値: E[·]
- 数式番号: 右寄せ (1), (2), ...
- 参照: as shown in (1) or from (1)
```

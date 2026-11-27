---
title: 固有値と対角化
description: 固有値・固有ベクトル・特性多項式・対角化
notation:
  - symbol: '\Phi_A(t)'
    meaning: 行列 A の特性多項式
---

## 固有値と固有ベクトル

:::definition[固有値]{#def-linalg-eigenvalue}
$V$ を体 $K$ 上のベクトル空間、$f\colon V\to V$ を線形写像とする。$f(v)=\lambda v$ を満たす $v\neq 0$ が存在するとき、$\lambda\in K$ を $f$ の :term[固有値]{en="eigenvalue" yomi="こゆうち"} という。
:::

$$
\Phi_A(t) = \det(tI - A) \label{eq-linalg-charpoly}
$$

:::theorem[Cayley–Hamilton]{#thm-linalg-cayley-hamilton}
$\Phi_A(A) = O$ である。
:::

:::proof
:eqref[eq-linalg-charpoly] を用いる。詳細は :cite[saito-linalg]{loc="第4章"} を参照。
:::

:ref[def-linalg-eigenvalue] と :ref[thm-linalg-cayley-hamilton]{full} を参照。時刻 12:30 と比 a:b。

```tikzcd caption="可換図式"
V \arrow[r, "f"] \arrow[d, "\varphi"'] \arrow[dr, "h"] & W \arrow[d, "\psi"] \\
V' \arrow[r, "g"'] & W'
```

## 演習

::::exercise{level=basic}
$A=\begin{pmatrix}2&1\\1&2\end{pmatrix}$ の固有値を求めよ。
:::solution
$\Phi_A(t)=(t-1)(t-3)$ より $1, 3$。
:::
::::

## 確認テスト

::::quiz
$2\times 2$ 行列 $A$ の固有値の和に等しいのはどれか。

- [ ] $\det A$
- [x] $\tr A$
- [ ] $\rank A$

:::explanation
特性多項式の $t$ の係数を比較する。
:::
::::

::::quiz{answer="3"}
上の演習の行列の最大固有値はいくつか。
::::

::interactive[行列の作用]{name="matrix-action"}

::::quiz{answer="1"}
ダミー
::::

---
title: 数列と級数
description: ε–N 論法による数列の極限と無限大への発散、単調収束・Bolzano–Weierstrass・Cauchy 完備性と上限の公理の同値性、上極限と下極限、正項級数と各種判定法、Abel の変形と Dirichlet の判定法、絶対収束と再配列、級数の積を扱う。
notation:
  - symbol: 'a_n\to a'
    meaning: 数列の極限が a であること
  - symbol: 'a_n\to+\infty'
    meaning: 数列が正の無限大に発散すること
  - symbol: 'a_{n_k}'
    meaning: 添字を狭義増加に選んだ部分列の第 k 項
  - symbol: '\limsup_{n\to\infty}a_n'
    meaning: 数列の上極限
  - symbol: '\liminf_{n\to\infty}a_n'
    meaning: 数列の下極限
  - symbol: 's_n=\sum_{k=1}^n a_k'
    meaning: 級数の第 n 部分和
  - symbol: 'e'
    meaning: 複利の極限として定める定数（Napier の数）
---

## この章の目標

- 誤差を指定して添字を選び、数列の収束と無限大への発散、四則演算の極限法則、基本的な極限を証明する。
- 単調収束定理、Bolzano–Weierstrass の定理、Cauchy 完備性を証明し、これらが順序体において上限の公理と同値であることを示す。
- 上極限と下極限によって振動する数列を記述し、その特徴付けと演算に関する性質を導く。
- 正項級数の比較・凝集・比・べき根の判定法、Abel の変形による Dirichlet の判定法と Leibniz の判定法を証明し、使い分ける。
- 絶対収束と条件収束を区別し、再配列と級数の積（Cauchy 積）について何が成り立つかを証明する。

前章では、集合の上限を求めるのに、その境界にどこまで近づけるかを調べた。
数列では項に順序がある。極限を調べるには、いったん境界に近づいた後も近くにとどまるかを確かめる必要がある。
以下、数列の項は実数とし、添字は特に断らない限り $1$ から始める。

## 数列の極限

### 誤差と添字の順序

:::definition[数列の収束]{#def-anl-sequence-limit}
写像 $a:\N\to\R$ を実数列といい、$a(n)=a_n$ と書く。ある $L\in\R$ について

$$
\forall\eps>0\quad\exists N\in\N\quad\forall n\ge N
\quad\abs{a_n-L}<\eps \label{eq-anl-sequence-limit}
$$

が成り立つとき、$(a_n)$ は $L$ に :term[収束]{en="convergence" yomi="しゅうそく"} するという。
$L$ をその :term[極限]{en="limit" yomi="きょくげん"} と呼び、$a_n\to L$ または $\lim_{n\to\infty}a_n=L$ と書く。
収束しない数列は :term[発散]{en="divergence" yomi="はっさん"} するという。
:::

式 :eqref[eq-anl-sequence-limit] では $N$ は $\eps$ に応じて選べるが、一つの $N$ がそれ以降のすべての添字に通用しなければならない。
従って「$a_n$ が $L$ に収束しない」とは、ある $\eps_0>0$ に対して、どの $N$ の後にも $\abs{a_n-L}\ge\eps_0$ となる $n\ge N$ があることをいう。
また、$\abs{a_n-L}<\eps$ を $\abs{a_n-L}\le\eps$ や $\abs{a_n-L}<C\eps$（$C>0$ は $\eps$ によらない定数）に替えても同じ条件になる。$\eps$ の代わりに $\eps/2$ や $\eps/C$ を指定すればよいからである。

:::example[有理式の誤差から添字を決める]
$a_n=(3n-2)/(2n+1)$ の極限を $3/2$ と予想する。差を通分すると

$$
\abs{a_n-\frac32}=\frac{7}{4n+2}<\frac{7}{4n}.
$$

従って $N=\lfloor7/(4\eps)\rfloor+1$ とすれば、$n\ge N$ で右辺は $\eps$ 未満となる。
例えば $\eps=1/100$ なら $N=176$ とできる。定義は $N$ の最小性を要求しない。
:::

:::proposition[極限の一意性と有界性]{#prop-anl-limit-basics}
収束列の極限は一意であり、収束列は有界である。数列が有界であるとは、ある $M$ に対してすべての $n$ で $\abs{a_n}\le M$ となることである。
:::

:::proof
$a_n\to a$ かつ $a_n\to b$ で $a\ne b$ とする。$\eps=\abs{a-b}/2$ を指定し、両方の収束に対応する添字の大きい方を $N$ とする。
すると $\abs{a-b}\le\abs{a-a_N}+\abs{a_N-b}<2\eps=\abs{a-b}$ となり矛盾する。

次に $a_n\to a$ で誤差を $1$ と指定すると、ある $N$ 以降で $\abs{a_n}\le\abs{a}+\abs{a_n-a}<\abs{a}+1$ となる。
この値と、有限個の初めの項の絶対値 $\abs{a_j}$（$j<N$）の最大をとれば、すべての項を抑えられる。
:::

### 演算と不等式

:::theorem[極限の四則演算]{#thm-anl-limit-laws}
$a_n\to a,\ b_n\to b$ とする。このとき $a_n+b_n\to a+b$、$a_n-b_n\to a-b$、$a_nb_n\to ab$ である。
さらに $b\ne0$ ならば十分後のすべての $n$ で $b_n\ne0$ となり、その範囲で $a_n/b_n\to a/b$ である。
:::

:::proof
$\eps>0$ とする。和と差については、各列に誤差 $\eps/2$ を指定する。対応する二つの添字の最大以降で、和や差の誤差は三角不等式により $\eps$ 未満となる。

積については :ref[prop-anl-limit-basics] により $\abs{a_n}\le M$ となる $M\ge1$ を選ぶ。
次のように差を二つに分ける。

$$
\abs{a_nb_n-ab}\le\abs{a_n}\abs{b_n-b}+\abs{b}\abs{a_n-a}\le M\abs{b_n-b}+\abs{b}\abs{a_n-a}.
$$

そこで $b_n$ の誤差を $\eps/(2M)$ 未満、$a_n$ の誤差を $\eps/(2(\abs{b}+1))$ 未満とすればよい。
両方の条件は、対応する添字の最大をとると同時に成り立ち、そのとき右辺は $\eps/2+\eps/2=\eps$ 未満である。

$b\ne0$ の場合には $\abs{b_n-b}<\abs{b}/2$ となる範囲で $\abs{b_n}\ge\abs{b}-\abs{b_n-b}>\abs{b}/2$ である。従って

$$
\abs{\frac1{b_n}-\frac1b}
=\frac{\abs{b_n-b}}{\abs{b_n}\abs{b}}
\le\frac{2}{\abs{b}^2}\abs{b_n-b}.
$$

$\abs{b_n-b}<\min\{\abs{b}/2,\eps\abs{b}^2/2\}$ と選べば逆数の収束が従う。
最後に積の法則を $a_n$ と $1/b_n$ に適用する。初めの有限個の項は収束に影響しない。
:::

定数列 $c$ は $c$ に収束するので、特に $ca_n\to ca$ である。また三角不等式 $\bigl|\abs{a_n}-\abs{a}\bigr|\le\abs{a_n-a}$ より、$a_n\to a$ なら $\abs{a_n}\to\abs{a}$ である。

:::proposition[順序と挟み撃ち]{#prop-anl-squeeze}
十分後で $a_n\le b_n$ かつ $a_n\to a,\ b_n\to b$ ならば $a\le b$ である。
また十分後で $u_n\le x_n\le v_n$ かつ $u_n,v_n\to L$ ならば $x_n\to L$ である。後者を :term[挟み撃ちの定理]{en="squeeze theorem" yomi="はさみうちのていり"} という。
:::

:::proof
もし $a>b$ ならば $\eps=(a-b)/2$ とおくと、十分後で $a_n>a-\eps=b+\eps>b_n$ となり仮定に反する。
後者では、任意の $\eps>0$ をとる。三つの不等式に必要な添字をまとめて大きくすれば
$L-\eps<u_n\le x_n\le v_n<L+\eps$ となり、定義から従う。
:::

前半で $a_n<b_n$ を仮定しても、結論は $a<b$ にはならない。例えば $0<1/n$ であるが、両辺の極限はともに $0$ である。極限をとると狭義の不等号は等号付きに弱まる。

:::example[符号が変わる列]
$x_n=(-1)^n/n$ は $-1/n\le x_n\le1/n$ を満たす。$N=\lfloor1/\eps\rfloor+1$ とすれば $n\ge N$ で $\abs{x_n}<\eps$ となり、極限は $0$ である。
これに対して $(-1)^n$ では、どこまで後に進んでも距離 $2$ の二項が残る。
もし $L$ に収束すると仮定すれば、誤差 $1$ に対する $N$ の後の隣り合う二項は $\abs{(-1)^n-(-1)^{n+1}}\le\abs{(-1)^n-L}+\abs{L-(-1)^{n+1}}<2$ を満たし、距離が $2$ であることに矛盾する。
:::

### 無限大への発散

:::definition[無限大への発散]{#def-anl-diverge-infinity}
数列 $(a_n)$ が

$$
\forall M\in\R\quad\exists N\in\N\quad\forall n\ge N\quad a_n>M
$$

を満たすとき、$(a_n)$ は正の無限大に発散するといい、$a_n\to+\infty$ または $\lim_{n\to\infty}a_n=+\infty$ と書く。
$-a_n\to+\infty$ のとき $a_n\to-\infty$ と書く。
:::

$a_n\to+\infty$ は収束ではなく発散の一種である。$\lim a_n=+\infty$ と書くのは、拡張実数 $\overline{\R}$（:ref[def-anl-extended-reals]）の中での極限とみなす便宜による。
発散する数列には、$(-1)^n$ のように振動して $\pm\infty$ のどちらにも向かわないものもある。

:::proposition[無限大への発散の性質]{#prop-anl-infinite-limits}
1. $a_n\to+\infty$ で $(b_n)$ が下に有界ならば $a_n+b_n\to+\infty$ である。
2. $a_n\to+\infty$ で、ある $c>0$ について十分後で $b_n\ge c$ ならば $a_nb_n\to+\infty$ である。
3. 十分後で $a_n\le b_n$ かつ $a_n\to+\infty$ ならば $b_n\to+\infty$ である。
4. すべての $n$ で $a_n>0$ とする。$a_n\to+\infty$ であることと $1/a_n\to0$ であることは同値である。
:::

:::proof
1. $b_n\ge-K$ とする。$M$ に対して $a_n>M+K$ となる範囲で $a_n+b_n>M$ である。
2. $M$ に対して、$a_n>\max\{M,0\}/c$ かつ $b_n\ge c$ となる範囲で $a_nb_n\ge a_nc>M$ である（$a_n>0$ に注意）。
3. $a_n>M$ となる範囲で $b_n\ge a_n>M$ である。
4. $a_n\to+\infty$ とし $\eps>0$ をとる。$a_n>1/\eps$ となる範囲で $0<1/a_n<\eps$ である。
逆に $1/a_n\to0$ とし $M$ をとる。$M\le0$ なら全項が $M$ を超える。$M>0$ なら $1/a_n<1/M$ となる範囲で $a_n>M$ である。
:::

### 基本的な極限

以下の極限は、級数の判定法や後の章で繰り返し用いる。証明には :ref[prop-anl-bernoulli] と :ref[prop-anl-binomial] を使う。

:::proposition[基本的な極限]{#prop-anl-standard-limits}
1. $\abs r<1$ ならば $r^n\to0$ である。$r>1$ ならば $r^n\to+\infty$ である。
2. $a>0$ ならば $\sqrt[n]a\to1$ である。
3. $\sqrt[n]n\to1$ である。
4. $\abs r<1$ と $k\in\N$ に対して $n^kr^n\to0$ である。
5. 任意の $a\in\R$ に対して $a^n/n!\to0$ である。
:::

:::proof
1. $r=0$ なら明らかである。$0<\abs r<1$ なら $h=1/\abs r-1>0$ とおくと、Bernoulli の不等式から

$$
\abs r^n=\frac1{(1+h)^n}\le\frac1{1+nh}<\frac1{nh}
$$

であり、$n>1/(h\eps)$ で $\abs{r^n}<\eps$ となる。$r>1$ なら $r=1+h$ として $r^n\ge1+nh$ であるから、$n>M/h$ で $r^n>M$ である。

2. $a>1$ とする。$x_n=\sqrt[n]a-1$ は正である（$\sqrt[n]a\le1$ なら $a\le1$ となるため）。Bernoulli の不等式により $a=(1+x_n)^n\ge1+nx_n$ なので $0<x_n\le(a-1)/n$ であり、$x_n\to0$ である。
$a=1$ なら全項が $1$ である。$0<a<1$ なら $\sqrt[n]a\cdot\sqrt[n]{1/a}=\sqrt[n]1=1$ より $\sqrt[n]a=1/\sqrt[n]{1/a}\to1$ である（:ref[thm-anl-limit-laws]）。

3. $x_n=\sqrt[n]n-1\ge0$ とおく。$n\ge2$ のとき二項定理で $k=2$ の項だけを残すと

$$
n=(1+x_n)^n\ge\binom n2x_n^2=\frac{n(n-1)}2x_n^2
$$

となり、$0\le x_n\le\sqrt{2/(n-1)}$ である。$\eps>0$ に対し $n>1+2/\eps^2$ なら $x_n<\eps$ となる。

4. $\abs r=1/(1+h)$、$h>0$ と書く。$n\ge2(k+1)$ なら $n-k\ge n/2$ なので、二項定理の $k+1$ 次の項だけを残して

$$
(1+h)^n\ge\binom n{k+1}h^{k+1}=\frac{n(n-1)\cdots(n-k)}{(k+1)!}h^{k+1}\ge\frac{(n/2)^{k+1}}{(k+1)!}h^{k+1}
$$

を得る。従って $0\le n^k\abs r^n\le C/n$、$C=2^{k+1}(k+1)!/h^{k+1}$ であり、右辺は $0$ に収束する。

5. 整数 $K>2\abs a$ をとる。$n>K$ ならば

$$
\frac{\abs a^n}{n!}=\frac{\abs a^K}{K!}\prod_{j=K+1}^n\frac{\abs a}{j}\le\frac{\abs a^K}{K!}\left(\frac12\right)^{n-K}
$$

であり、1 により右辺は $0$ に収束する。挟み撃ちで結論を得る。
:::

4 は「指数関数的な減衰はどんな多項式的な増大にも勝つ」ことを、5 は「階乗はどんな指数的な増大にも勝つ」ことを述べている。

## 単調収束定理

### 上限を数列の極限とみる

:::definition[単調な数列]
すべての $n$ で $a_n\le a_{n+1}$ を満たす数列を :term[単調増加列]{en="nondecreasing sequence" yomi="たんちょうぞうかれつ"}、$a_n\ge a_{n+1}$ を満たす数列を :term[単調減少列]{en="nonincreasing sequence" yomi="たんちょうげんしょうれつ"} という。
ここでは等号を許す。等号を許さない場合は「狭義」と明記する。いずれかであるとき単調列という。
:::

:::theorem[単調収束定理]{#thm-anl-monotone-convergence}
上に有界な単調増加列は収束し、その極限は項全体の上限である。下に有界な単調減少列の極限は項全体の下限である。
:::

:::proof
増加列の場合、:ref[ax-anl-supremum] により $s=\sup\{a_n\mid n\in\N\}$ が存在する。
$\eps>0$ に対し :ref[thm-anl-bound-epsilon] で $s-\eps<a_N$ となる $N$ を選ぶ。
$n\ge N$ なら $s-\eps<a_N\le a_n\le s$ なので $\abs{a_n-s}<\eps$ である。
減少列に対しては増加列 $(-a_n)$ に適用し、:ref[prop-anl-infimum] により符号を戻せばよい。
:::

有界性を外すと、上に有界でない増加列は任意の $M$ に対して $a_N>M$ を満たす一項をもち、$n\ge N$ で $a_n\ge a_N>M$ が続く。すなわち $a_n\to+\infty$ である。
従って単調増加列は、拡張実数の意味で常に極限 $\sup_n a_n\in(-\infty,+\infty]$ をもつ。単調減少列も同様に $\inf_n a_n$ に向かう。

単調収束定理は極限の値を知らなくても存在を保証する。漸化式で定まる数列の極限を求めるときには、まず単調性と有界性で存在を示し、次に漸化式の極限をとって値を決める。

:::example[√2 の Newton 近似]{#exm-anl-newton-sqrt2}
$x_1=2$、$x_{n+1}=\dfrac12\left(x_n+\dfrac2{x_n}\right)$ とおく。帰納的に $x_n>0$ である。
$\sqrt2$ との差について

$$
x_{n+1}-\sqrt2=\frac{x_n^2-2\sqrt2x_n+2}{2x_n}=\frac{(x_n-\sqrt2)^2}{2x_n}\ge0,\qquad
x_{n+1}-x_n=\frac{2-x_n^2}{2x_n}
$$

が成り立つ。前者から $n\ge2$ で $x_n\ge\sqrt2$（$x_1=2$ も同様）、従って後者から $x_{n+1}\le x_n$ である。
単調減少で $\sqrt2$ で下から抑えられるので、極限 $L\ge\sqrt2>0$ が存在する。漸化式で極限をとると $L=(L+2/L)/2$、すなわち $L^2=2$ であり、$L=\sqrt2$ である。
さらに前者の等式と $x_n\ge\sqrt2$ から $0\le x_{n+1}-\sqrt2\le(x_n-\sqrt2)^2/(2\sqrt2)$ となり、誤差は一段ごとにおよそ二乗される。
実際 $x_2=3/2$、$x_3=17/12$、$x_4=577/408$ であり、$x_4-\sqrt2$ は $2\times10^{-6}$ 程度である。
:::

### 定数 e の存在

:::proposition[複利の極限]{#prop-anl-e-limit}
$u_n=(1+1/n)^n$ は狭義単調増加かつ上に有界である。その極限を $e$ と書く。
:::

:::proof
二項定理で $(1+1/n)^n$ を展開し、$\binom nk/n^k$ を積で書くと

$$
u_n=\sum_{k=0}^n\frac{1}{k!}
 \prod_{j=0}^{k-1}\left(1-\frac{j}{n}\right)
 \label{eq-anl-e-binomial}
$$

となる。ここで因子がない積は $1$ とする。$k\le n$ の各因子は正であり、$n$ を $n+1$ に替えても減少しない。
さらに $k=n+1$ の正の項が加わるので $u_n<u_{n+1}$ である。

各因子は $1$ 以下である。また $k\ge1$ なら $k!\ge2^{k-1}$ であるから（$k!=\prod_{j=2}^kj$ の各因子が $2$ 以上）

$$
u_n\le1+\sum_{k=1}^n 2^{1-k}=3-2^{1-n}<3.
$$

有限和の等号は和を半分にしたものを元の和から引いて得られる。
よって :ref[thm-anl-monotone-convergence] により極限 $e$ が存在する。
:::

式 :eqref[eq-anl-e-binomial] は極限の存在を指数関数に頼らずに示している。
例えば $u_1=2,\ u_2=9/4$ なので $9/4<e\le3$ ともわかる。
この定数と指数関数の関係は [初等関数の厳密な構成](/math/analysis/elementary-functions/) で扱う。

:::proposition[e の級数表示と誤差]{#prop-anl-e-series}
$s_n=\sum_{k=0}^n1/k!$ とおくと $s_n\to e$ であり、すべての $n\in\N$ について

$$
0<e-s_n<\frac1{n!\,n}.
$$
:::

:::proof
$s_n$ は狭義増加であり、上と同じ評価で $s_n<3$ なので、極限 $\sigma$ が存在する。
:eqref[eq-anl-e-binomial] の各積は $1$ 以下なので $u_n\le s_n$、極限をとって $e\le\sigma$ である。
逆に $m$ を固定すると、$n\ge m$ では $k>m$ の項を捨てて

$$
u_n\ge\sum_{k=0}^m\frac1{k!}\prod_{j=0}^{k-1}\left(1-\frac jn\right)
$$

である。右辺は $m+1$ 個の項の和で、各積は $n\to\infty$ で $1$ に収束するので、右辺は $s_m$ に収束する。:ref[prop-anl-squeeze] により $e\ge s_m$ がすべての $m$ で成り立ち、$e\ge\sigma$ を得る。従って $e=\sigma$ である。

誤差を評価する。$m>n$ なら、$k\ge n+1$ で $k!\ge(n+1)!\,(n+1)^{k-n-1}$ であるから

$$
s_m-s_n=\sum_{k=n+1}^m\frac1{k!}\le\frac1{(n+1)!}\sum_{i=0}^{m-n-1}\frac1{(n+1)^i}<\frac1{(n+1)!}\cdot\frac{n+1}{n}
$$

となる。$m\to\infty$ として $e-s_n\le\dfrac1{n!\,n}$ を得る。これを $n+1$ に適用すると

$$
e-s_n=\frac1{(n+1)!}+(e-s_{n+1})\le\frac1{(n+1)!}\left(1+\frac1{n+1}\right)<\frac1{(n+1)!}\cdot\frac{n+1}n=\frac1{n!\,n}
$$

となる。最後の不等号は $n(n+2)<(n+1)^2$ による。$e-s_n>0$ は $s_m$ が狭義増加であることから従う。
:::

この誤差評価は収束が非常に速いことを示す。例えば $n=10$ で誤差は $1/(10!\cdot10)<3\times10^{-8}$ である。

:::theorem[e の無理性]{#thm-anl-e-irrational}
$e$ は無理数である。
:::

:::proof
$e=p/q$（$p\in\Z$、$q\in\N$）と仮定する。$q!\,e=(q-1)!\,p$ と $q!\,s_q=\sum_{k=0}^qq!/k!$ はともに整数であるから、$x=q!\,(e-s_q)$ も整数である。
ところが :ref[prop-anl-e-series] より $0<x<q!/(q!\,q)=1/q\le1$ となり、$0$ と $1$ の間に整数はないことに反する。
:::

## Bolzano–Weierstrass の定理

### 有界な数列から部分列を選ぶ

:::definition[部分列]{#def-anl-subsequence}
自然数 $n_1<n_2<\cdots$ を選んで得る数列 $(a_{n_k})_{k\in\N}$ を $(a_n)$ の :term[部分列]{en="subsequence" yomi="ぶぶんれつ"} という。
部分列の極限として得られる実数を元の数列の :term[集積値]{en="subsequential limit" yomi="しゅうせきち"} という。
:::

$n_1\ge1$ であり、$n_k\ge k$ なら $n_{k+1}\ge n_k+1\ge k+1$ であるから、帰納法ですべての $k$ について $n_k\ge k$ である。
これを用いると、$a_n\to L$ ならばすべての部分列も $L$ に収束する。
実際、誤差 $\eps>0$ に対応する元の添字を $N$ とすれば、$k\ge N$ のとき $n_k\ge k\ge N$ である。

:::proposition[集積値の特徴付け]{#prop-anl-cluster-value}
$c\in\R$ が $(a_n)$ の集積値であるための必要十分条件は、任意の $\eps>0$ と任意の $N$ に対して $\abs{a_n-c}<\eps$ となる $n\ge N$ が存在することである。
:::

:::proof
$a_{n_k}\to c$ とする。$\eps,N$ に対して $k\ge N$ かつ $\abs{a_{n_k}-c}<\eps$ となる $k$ をとれば、$n=n_k\ge k\ge N$ が条件を満たす。
逆に条件を仮定する。$\eps=1,N=1$ として $n_1$ を選び、$n_{j-1}$ まで選んだら $\eps=1/j$、$N=n_{j-1}+1$ として $n_j$ を選ぶ。
$n_1<n_2<\cdots$ かつ $\abs{a_{n_j}-c}<1/j$ なので、部分列は $c$ に収束する。
:::

:::lemma[単調部分列の存在]{#lem-anl-monotone-subsequence}
任意の実数列は単調な部分列をもつ。
:::

:::proof
添字 $m$ が「すべての $n>m$ について $a_m\ge a_n$」を満たすとき、$m$ を山頂と呼ぶ。

山頂が無限にあれば、それらを小さい順に $m_1<m_2<\cdots$ と並べると、山頂の定義により $a_{m_1}\ge a_{m_2}\ge\cdots$ となり、単調減少な部分列を得る。
山頂が有限個なら、それらより大きい添字 $n_1$ をとる。$n_1$ は山頂でないので $a_{n_2}>a_{n_1}$ となる $n_2>n_1$ がある。$n_2$ も山頂でないので同様に $n_3>n_2$ がとれる。
これを繰り返して狭義単調増加な部分列を得る。
:::

この補題の証明は数の大小しか使わないので、任意の順序体で成り立つ。

:::theorem[Bolzano–Weierstrass の定理]{#thm-anl-bolzano-weierstrass}
有界な実数列は、少なくとも一つの収束部分列をもつ。
:::

:::proof
:ref[lem-anl-monotone-subsequence] により単調な部分列がとれる。それは有界であるから、:ref[thm-anl-monotone-convergence] により収束する。
:::

この短い証明のほかに、区間の二分による証明もある。こちらは部分列の極限の位置を具体的に絞り込む。

:::proof[区間の二分による別証明]
すべての項を含む有界閉区間 $I_1=[u,v]$ をとる。中央で二つの閉区間に分けると、少なくとも一方には無限個の添字に対応する項が入る。
実際、両方に有限個の添字しか入らなければ、それらを合わせても有限個である。値が同じ項も添字ごとに数える。

このような側を $I_2$ とし、同じ操作を繰り返す。各 $I_k$ には無限個の添字があるので、$a_{n_k}\in I_k$ を満たす $n_k$ を一つ前より大きく選べる。
区間は入れ子で長さは $(v-u)2^{1-k}$ であるから、:ref[thm-anl-nested-intervals] によって共通点 $L$ がある。

$L,a_{n_k}\in I_k$ なので $\abs{a_{n_k}-L}\le(v-u)2^{1-k}\le(v-u)/k$ である。
任意の $\eps>0$ に対して $K=\lfloor(v-u)/\eps\rfloor+1$ とすれば、$k\ge K$ でこの誤差は $\eps$ 未満になる。
:::

:::example[全体の極限と部分列の極限]
$x_n=(-1)^n+1/n$ の偶数番目の項は $1+1/(2k)$、奇数番目の項は $-1+1/(2k-1)$ であるから、それぞれ $1,-1$ に収束する。
元の列が収束するなら両方の極限は一致するはずなので、元の列は発散する。
それでもすべての項は $[-2,2]$ に入り、上の定理の仮定を満たす。
:::

有界性は外せない。$a_n=n$ のどの部分列も $a_{n_k}=n_k\ge k$ なので有界でなく、:ref[prop-anl-limit-basics] により収束しない。

## Cauchy 列と完備性

### 極限の候補を使わない条件

:::definition[Cauchy 列]{#def-anl-cauchy}
$(a_n)$ が :term[Cauchy 列]{en="Cauchy sequence" yomi="こーしーれつ"} であるとは

$$
\forall\eps>0\quad\exists N\in\N\quad\forall m,n\ge N
\quad\abs{a_m-a_n}<\eps
$$

となることをいう。隣り合う二項だけでなく、$N$ 以降のどの二項を選んでもよい。
:::

隣り合う二項の差だけが $0$ に近づいても Cauchy 列とは限らない。例えば $H_n=\sum_{k=1}^n1/k$ では $H_{n+1}-H_n=1/(n+1)\to0$ であるが、後で示すように $H_{2n}-H_n\ge1/2$ であり、Cauchy 列でない。

:::theorem[実数の Cauchy 完備性]{#thm-anl-cauchy-complete}
実数列が収束するための必要十分条件は、Cauchy 列であることである。
:::

:::proof
$a_n\to L$ ならば任意の $\eps>0$ に対し、ある $N$ 以降で $\abs{a_n-L}<\eps/2$ である。
従って $m,n\ge N$ なら $\abs{a_m-a_n}\le\abs{a_m-L}+\abs{a_n-L}<\eps$ となる。

逆に $(a_n)$ を Cauchy 列とする。まず誤差 $1$ に対応する $N$ をとると、$m\ge N$ では $\abs{a_m}\le\abs{a_N}+1$ である。
初めの有限個の項も含めた上界をとれば列は有界なので、:ref[thm-anl-bolzano-weierstrass] によってある部分列 $a_{n_k}$ が $L$ に収束する。

任意の $\eps>0$ に対して Cauchy 条件を誤差 $\eps/2$ で満たす $N$ をとる。
さらに $n_k\ge N$ かつ $\abs{a_{n_k}-L}<\eps/2$ となる $k$ を一つ固定する。
するとすべての $n\ge N$ で $\abs{a_n-L}\le\abs{a_n-a_{n_k}}+\abs{a_{n_k}-L}<\eps$ となる。
:::

:::example[有理数だけを許す場合]
$q_n=\lfloor 2^n\sqrt{2}\rfloor/2^n$ は有理数列であり、$0\le\sqrt{2}-q_n<2^{-n}$ である。
従って $\R$ では $\sqrt{2}$ に収束し、特に Cauchy 列である。しかし $\Q$ の中では極限をもたない。
もし有理数に収束すれば $\R$ でも同じ数に収束するので、極限の一意性と :ref[thm-anl-sqrt-two] に反する。
:::

:::example[縮小的な数列]{#exm-anl-contractive}
定数 $0<c<1$ について、すべての $n$ で $\abs{a_{n+2}-a_{n+1}}\le c\abs{a_{n+1}-a_n}$ が成り立つとする。
帰納的に $\abs{a_{k+1}-a_k}\le c^{k-1}d$（$d=\abs{a_2-a_1}$）であるから、$m>n$ なら

$$
\abs{a_m-a_n}\le\sum_{k=n}^{m-1}\abs{a_{k+1}-a_k}\le d\sum_{k=n}^{m-1}c^{k-1}<\frac{dc^{n-1}}{1-c}
$$

となり、:ref[prop-anl-standard-limits] の 1 から $(a_n)$ は Cauchy 列、従って収束する。

例えば $a_1=1$、$a_{n+1}=1+1/a_n$ とすると帰納的に $a_n\ge1$ であり、$a_{n}a_{n-1}=a_{n-1}+1\ge2$ から

$$
\abs{a_{n+1}-a_n}=\abs{\frac1{a_n}-\frac1{a_{n-1}}}=\frac{\abs{a_n-a_{n-1}}}{a_na_{n-1}}\le\frac12\abs{a_n-a_{n-1}}
$$

である。極限 $L\ge1$ は $L=1+1/L$ を満たすので $L^2-L-1=0$、$L=(1+\sqrt5)/2$ である。
この数列は単調ではない（$1,2,3/2,5/3,\dots$）ので、単調収束定理は直接には使えない。
:::

### 完備性の各形式から上限の公理へ

[実数の連続性](/math/analysis/real-numbers/) では上限の公理から出発した。ここで証明の出発点を替え、順序体の中で上限を取り戻せることを確かめる。
以下では順序体 $K$ の中で数列の収束、部分列、Cauchy 列を同じ式で定義する（誤差 $\eps$ は $K$ の正の元を動く）。

:::remark[順序体で成り立つ性質]
:ref[prop-anl-limit-basics]、:ref[thm-anl-limit-laws]、:ref[prop-anl-squeeze]、部分列が同じ極限に収束すること、:ref[lem-anl-monotone-subsequence]、および :ref[thm-anl-cauchy-complete] の証明の前半（収束列は Cauchy 列）と「Cauchy 列は有界」の部分は、順序体の公理と :ref[thm-anl-bound-epsilon] だけを用いて証明した。
従ってこれらは任意の順序体で成り立ち、以下の証明で用いても循環論法にはならない。
これに対し、上限・単調収束・Bolzano–Weierstrass・Archimedes の性質は、以下で仮定として明示したときにだけ用いる。
:::

:::proposition[完備性の同値な条件]{#prop-anl-completeness-equivalence}
順序体 $K$ について、次の四条件は互いに同値である。

1. 上限の公理：空でない上に有界な部分集合は上限をもつ。
2. 単調収束：上に有界な単調増加列は収束する。
3. Bolzano–Weierstrass の性質：有界な数列は収束部分列をもつ。
4. Cauchy 完備性と Archimedes の性質：Cauchy 列は収束し、かつ $\N$ は上に有界でない。
:::

:::proof
**1 ⇒ 2。** :ref[thm-anl-monotone-convergence] の証明は上限の公理と :ref[thm-anl-bound-epsilon] しか使わない。

**2 ⇒ 3。** 有界な数列は :ref[lem-anl-monotone-subsequence] により単調な部分列をもち、その部分列も有界である。
増加なら 2 により収束し、減少なら符号を替えた増加列に 2 を適用すればよい。

**3 ⇒ 4。** まず Archimedes の性質を示す。$\N$ が上に有界なら、数列 $x_n=n$ は $1\le x_n\le u$ で有界であり、3 により収束部分列 $n_k\to L$ をもつ。
収束列は Cauchy 列なので、ある $K$ 以降で $\abs{n_{k+1}-n_k}<1/2$ となる。しかし $n_{k+1}>n_k$ は自然数なので :ref[prop-anl-natural-basics] の 5 により $n_{k+1}-n_k\ge1$ であり、矛盾する。
次に Cauchy 列は有界なので、3 により収束部分列をもつ。:ref[thm-anl-cauchy-complete] の証明の後半は三角不等式しか使わないので、Cauchy 列そのものが同じ極限に収束する。

**4 ⇒ 1。** 入れ子の有界閉区間 $[l_n,u_n]$ で、各 $\eps>0$ に対し長さが $\eps$ 未満になるものをとる。
$\eps>0$ に対し $u_N-l_N<\eps$ と選べば、$m,n\ge N$ では $l_m,l_n\in[l_N,u_N]$ なので $\abs{l_m-l_n}<\eps$ である。従って $(l_n)$ は Cauchy 列であり、極限 $s$ をもつ。
各 $k$ について $n\ge k$ なら $l_k\le l_n\le u_k$ であり、:ref[prop-anl-squeeze] で極限をとると $l_k\le s\le u_k$ となる。従って共通点がある。
Archimedes の性質も仮定されているので、:ref[prop-anl-bisection-sup] により上限の公理が従う。
:::

[実数の連続性](/math/analysis/real-numbers/) の :ref[cor-anl-completeness-forms] と合わせると、上限の公理、Dedekind の性質、区間縮小法と Archimedes の性質、単調収束、Bolzano–Weierstrass の性質、Cauchy 完備性と Archimedes の性質の六つは、順序体においてすべて同値である。

## 上極限と下極限

### 後ろに残る項の上限と下限

:::definition[上極限と下極限]{#def-anl-limit-extrema}
$A_n=\{a_k\mid k\ge n\}$ とし、拡張実数 $\overline{\R}$ の中で（:ref[prop-anl-extended-sup]）$u_n=\sup A_n,\ l_n=\inf A_n$ とおく。

$$
\limsup_{n\to\infty}a_n=\inf_n u_n,\qquad
\liminf_{n\to\infty}a_n=\sup_n l_n .
$$

これらをそれぞれ :term[上極限]{en="limit superior" yomi="じょうきょくげん"}、:term[下極限]{en="limit inferior" yomi="かきょくげん"} という。
:::

$A_{n+1}\subset A_n$ より $u_n$ は単調減少、$l_n$ は単調増加である。
$A_n$ と $A_1$ は有限個の項しか違わないので、$(a_n)$ が上に有界でなければすべての $u_n$ が $+\infty$ であり、上極限は $+\infty$ である。
上に有界なら $u_n$ はすべて実数であり、単調収束定理と無限大への発散の議論により、$u_n$ は拡張実数の意味で上極限に収束する。
このとき上極限が $-\infty$ となるのは、任意の $M\in\R$ に対して $u_n<M$ となる $n$ があるとき、すなわち $u_n\to-\infty$ のときであり、$a_n\le u_n$ から $a_n\to-\infty$ となる。

:::proposition[集積値と収束の判定]{#prop-anl-limsup-cluster}
有界列について $l=\liminf a_n,\ u=\limsup a_n$ とおく。
集積値の最小値は $l$、最大値は $u$ である。また列が $L$ に収束することと $l=u=L$ であることは同値である。
:::

:::proof
$k\ge n$ なら $l_n\le a_k\le u_n$ であるから、どの収束部分列の極限 $b$ も各 $n$ に対して $l_n\le b\le u_n$ を満たす（部分列の十分後の項は $k\ge n$ を満たす）。
さらに $n$ の極限をとると $l\le b\le u$ となる。

$u$ 自身を実現するため、$n_0=0$ とし、各 $j\ge1$ で $m_j=\max\{j,n_{j-1}+1\}$ とおく。
上限の特徴付けにより、$k\ge m_j$ のうち

$$
u_{m_j}-1/j<a_k\le u_{m_j}
$$

となるものを $n_j$ に選べる。$n_j\ge m_j>n_{j-1}$ なので部分列である。$m_j\ge j$ なので $u\le u_{m_j}\le u_j$ から $u_{m_j}\to u$ であり、挟み撃ちにより $a_{n_j}\to u$ である。
同じ手順で $l_{m_j}\le a_k<l_{m_j}+1/j$ を満たす項を選ぶと $l$ に収束する部分列を得る。

列が $L$ に収束すればすべての部分列の極限が $L$ なので $l=u=L$ である。
逆に $l=u=L$ なら $l_n\le a_n\le u_n$ に挟み撃ちを用いて $a_n\to L$ となる。
:::

上に有界でない数列では、$a_{n_j}>j$ かつ $n_j>n_{j-1}$ となる添字を順に選べる（各 $A_n$ が上に有界でないため）。従って $+\infty$ に発散する部分列があり、上極限 $+\infty$ はこの意味で「最大の集積値」である。

例えば $a_n=(-1)^n+1/n$ では、偶数番目と奇数番目の項の極限はそれぞれ $1,-1$ である。
またすべての $k\ge n$ で $-1\le a_k\le1+1/n$ であるから、集積値は $[-1,1]$ に入る。
従って上極限は $1$、下極限は $-1$ である。

### 誤差による特徴付けと演算

:::theorem[上極限の ε による特徴付け]{#thm-anl-limsup-epsilon}
実数列 $(a_n)$ と $u\in\R$ について、$u=\limsup a_n$ であるための必要十分条件は、各 $\eps>0$ に対して次の二つが成り立つことである。

1. 十分後のすべての $n$ で $a_n<u+\eps$ である。
2. $a_n>u-\eps$ となる $n$ が無限に多くある。

下極限についても、不等号の向きを逆にした同様の特徴付けが成り立つ。
:::

:::proof
$u=\limsup a_n$ とする。$u_n\to u$ なので $u_N<u+\eps$ となる $N$ があり、$n\ge N$ で $a_n\le u_N<u+\eps$ である。
また各 $N$ について $u_N\ge u>u-\eps$ なので、上限の特徴付けにより $a_n>u-\eps$ となる $n\ge N$ がある。$N$ は任意なので、そのような $n$ は無限にある。

逆に 1, 2 を仮定する。1 から、ある $N$ 以降で $u+\eps$ は $A_N$ の上界なので $u_N\le u+\eps$、従って $\limsup a_n\le u+\eps$ である。
2 から、各 $n$ について $a_k>u-\eps$ となる $k\ge n$ があるので $u_n>u-\eps$、従って $\limsup a_n\ge u-\eps$ である。
$\eps$ は任意なので $\limsup a_n=u$ である。下極限は $(-a_n)$ に適用して符号を戻せばよい（次の命題の 2）。
:::

:::proposition[上極限・下極限の性質]{#prop-anl-limsup-properties}
1. $\liminf a_n\le\limsup a_n$ である。
2. $\limsup(-a_n)=-\liminf a_n$ である。
3. 十分後で $a_n\le b_n$ ならば $\limsup a_n\le\limsup b_n$、$\liminf a_n\le\liminf b_n$ である。
4. $(a_n),(b_n)$ が有界ならば
$\liminf a_n+\liminf b_n\le\liminf(a_n+b_n)\le\limsup(a_n+b_n)\le\limsup a_n+\limsup b_n$ である。
5. $b_n\to b\in\R$ かつ $(a_n)$ が有界ならば $\limsup(a_n+b_n)=\limsup a_n+b$ である。
:::

:::proof
1. 任意の $n,m$ について $p=\max\{n,m\}$ とおけば $l_n\le l_p\le u_p\le u_m$ である。:ref[prop-anl-sup-properties] の 4 により $\sup_nl_n\le\inf_mu_m$ である。
2. $\sup\{-a_k\mid k\ge n\}=-\inf A_n$ であり、$\inf_n(-l_n)=-\sup_nl_n$ である。
3. $a_n\le b_n$ が $n\ge N$ で成り立てば、$n\ge N$ で $\sup_{k\ge n}a_k\le\sup_{k\ge n}b_k$ である。減少列 $u_n$ の下限は $n\ge N$ の部分だけで決まるので結論を得る。下極限も同様である。
4. $k\ge n$ で $a_k+b_k\le\sup_{j\ge n}a_j+\sup_{j\ge n}b_j$ なので、左辺の $k\ge n$ についての上限も右辺以下である。$n\to\infty$ として :ref[thm-anl-limit-laws] を用いれば右端の不等式を得る。左端は 2 を用いて符号を替えればよく、中央は 1 である。
5. 4 と $\limsup b_n=b$（:ref[prop-anl-limsup-cluster]）から $\le$ が従う。逆向きは $a_n=(a_n+b_n)+(-b_n)$ に同じ不等式を用いて $\limsup a_n\le\limsup(a_n+b_n)-b$ を得る。
:::

4 の不等式は一般に等号にならない。$a_n=(-1)^n$、$b_n=-(-1)^n$ では $\limsup(a_n+b_n)=0$ であるが、$\limsup a_n+\limsup b_n=2$ である。

## 級数の収束と判定法

### 部分和と剰余

:::definition[級数と部分和]{#def-anl-series}
列 $(a_n)$ から $s_0=0,\ s_n=\sum_{k=1}^n a_k$ をつくり、$s_n$ を第 $n$ :term[部分和]{en="partial sum" yomi="ぶぶんわ"} という。
:term[級数]{en="series" yomi="きゅうすう"} $\sum_{n=1}^\infty a_n$ が収束するとは、部分和の列が有限な極限 $s$ をもつことであり、その極限を級数の和という。
部分和が収束しないとき級数は発散するといい、$s_n\to\pm\infty$ のとき $\sum a_n=\pm\infty$ と書く。
:::

:ref[thm-anl-cauchy-complete] を部分和に適用すると、収束の必要十分条件は、各 $\eps>0$ についてある $N$ があり

$$
m>n\ge N\quad\Longrightarrow\quad
\left|\sum_{k=n+1}^m a_k\right|<\eps \label{eq-anl-series-cauchy}
$$

となることである。特に $m=n+1$ とすれば、級数が収束するためには $a_n\to0$ が必要である。

:::proposition[級数の基本的な性質]{#prop-anl-series-basics}
1. $\sum a_n=s$、$\sum b_n=t$ が収束すれば、実数 $\alpha,\beta$ に対し $\sum(\alpha a_n+\beta b_n)=\alpha s+\beta t$ である。
2. 有限個の項を変更・追加・削除しても収束・発散は変わらない。収束する級数の剰余 $r_n=\sum_{k=n+1}^\infty a_k$ は収束し、$r_n=s-s_n\to0$ である。
3. 収束する級数の項を順序を変えずに括弧でまとめた級数は、同じ和に収束する。
:::

:::proof
1. 部分和について $\sum_{k=1}^n(\alpha a_k+\beta b_k)=\alpha s_n+\beta t_n$ であり、:ref[thm-anl-limit-laws] を適用する。
2. 有限個の項の変更は、ある番号以降の部分和に一定の数を加えるだけである。剰余の部分和は $\sum_{k=n+1}^ma_k=s_m-s_n$ なので $m\to\infty$ で $s-s_n$ に収束し、これは $n\to\infty$ で $0$ に近づく。
3. 括弧でまとめた級数の部分和は、元の部分和の部分列 $s_{n_1},s_{n_2},\dots$ である。収束列の部分列は同じ極限に収束する。
:::

3 の逆は成り立たない。$(1-1)+(1-1)+\cdots=0$ であるが、括弧を外した $1-1+1-1+\cdots$ の部分和は $1,0,1,0,\dots$ と振動する。
ただし項がすべて非負なら、部分和は単調なので、部分列の収束から全体の収束が従い、括弧を自由に外せる。

:::example[幾何級数]
初項 $1$、公比 $r$ の幾何級数 $\sum_{n=0}^\infty r^n$ を考える。$r\ne1$ なら、$1+r+\cdots+r^{n-1}$ に $1-r$ を掛けると中間項が消え、部分和は $(1-r^n)/(1-r)$ となる。
$\abs r<1$ なら :ref[prop-anl-standard-limits] の 1 により $r^n\to0$ なので、級数は $1/(1-r)$ に収束する。
$\abs r\ge1$ なら $\abs{r^n}\ge1$ なので一般項が $0$ に収束せず、級数は発散する。
:::

:::example[調和級数]
$a_n=1/n$ は $0$ に収束するが、$\sum_{n=1}^\infty 1/n$ は発散する。
実際、任意の $j\ge0$ に対して

$$
\sum_{k=2^j+1}^{2^{j+1}}\frac1k
\ge 2^j\frac1{2^{j+1}}=\frac12 .
$$

従って $s_{2^m}\ge1+m/2$ となる。部分和は増加するので $s_n\to+\infty$ である。
一般項が $0$ に収束することは収束の必要条件にすぎない。
:::

### 正項級数と比較

項がすべて非負の級数を正項級数という。

:::proposition[正項級数の収束]{#prop-anl-positive-series}
正項級数は、部分和が上に有界ならば収束し、そうでなければ $+\infty$ に発散する。
:::

:::proof
$a_n\ge0$ なので $s_{n+1}=s_n+a_{n+1}\ge s_n$ であり、部分和は単調増加である。:ref[thm-anl-monotone-convergence] とその後の注意（上に有界でない増加列は $+\infty$ に発散する）による。
:::

:::theorem[比較判定法]{#thm-anl-series-comparison}
十分大きな $n$ で $0\le a_n\le b_n$ ならば、$\sum b_n$ の収束から $\sum a_n$ の収束が従う。
また $\sum a_n$ が $+\infty$ に発散すれば、$\sum b_n$ も $+\infty$ に発散する。
:::

:::proof
有限個の初項を除き、不等式が全項で成り立つ場合を考える。除いた項の和は後の部分和に定数を加えるだけである。
$a_n$ の部分和は増加し、$b_n$ の対応する部分和以下である。
後者が収束すれば有界なので、前者は :ref[thm-anl-monotone-convergence] によって収束する。
前者がいくらでも大きくなる場合には、その上にある後者もいくらでも大きくなる。
:::

符号をもつ項を比較する際は、$\sum\abs{a_n}$ の収束を示せばよい。
実際、$\bigl|\sum_{k=n+1}^m a_k\bigr|\le\sum_{k=n+1}^m\abs{a_k}$ より、絶対値の級数が :eqref[eq-anl-series-cauchy] を満たせば元の級数も満たす。

:::corollary[極限による比較]{#cor-anl-limit-comparison}
$a_n\ge0$、$b_n>0$ で $a_n/b_n\to c$ とする。$0<c<\infty$ ならば、$\sum a_n$ と $\sum b_n$ は同時に収束または発散する。
$c=0$ ならば $\sum b_n$ の収束から $\sum a_n$ の収束が従う。
:::

:::proof
$0<c<\infty$ なら、十分後で $c/2<a_n/b_n<2c$、すなわち $(c/2)b_n<a_n<2cb_n$ であり、:ref[thm-anl-series-comparison] を両方向に用いる。
$c=0$ なら十分後で $a_n\le b_n$ である。
:::

例えば $\sum(2n+1)/(n^3+1)$ は $b_n=1/n^2$ と比べて $a_n/b_n\to2$ であり、次の例により収束する。

:::theorem[Cauchy の凝集判定法]{#thm-anl-condensation}
$a_n\ge0$ が単調減少ならば、$\sum_{n=1}^\infty a_n$ が収束するための必要十分条件は、$\sum_{k=0}^\infty2^ka_{2^k}$ が収束することである。
:::

:::proof
$s_n=\sum_{j=1}^na_j$、$t_k=\sum_{i=0}^k2^ia_{2^i}$ とおく。いずれも単調増加である。
$n<2^{k+1}$ ならば、項を $2^i\le j<2^{i+1}$ ごとにまとめ、各組で $a_j\le a_{2^i}$ を用いて

$$
s_n\le a_1+(a_2+a_3)+\cdots+(a_{2^k}+\cdots+a_{2^{k+1}-1})\le a_1+2a_2+\cdots+2^ka_{2^k}=t_k
$$

を得る。逆に $n\ge2^k$ ならば、$2^{i-1}<j\le2^i$ ごとにまとめ、各組で $a_j\ge a_{2^i}$ を用いて

$$
s_n\ge a_1+a_2+(a_3+a_4)+\cdots+(a_{2^{k-1}+1}+\cdots+a_{2^k})\ge\tfrac12a_1+a_2+2a_4+\cdots+2^{k-1}a_{2^k}=\tfrac12t_k
$$

を得る。従って $(s_n)$ と $(t_k)$ は一方が有界なら他方も有界であり、正項級数の性質から結論が従う。
:::

:::example[p 級数]
自然数 $p$ に対して $\sum1/n^p$ を考える。$a_n=1/n^p$ は単調減少であり

$$
2^ka_{2^k}=2^k\cdot2^{-kp}=\left(2^{1-p}\right)^k
$$

である。$p\ge2$ なら公比 $2^{1-p}\le1/2$ の幾何級数なので収束し、$p=1$ なら各項が $1$ なので発散する。
従って $\sum1/n^p$ は $p\ge2$ で収束し、$p=1$ で発散する。
実数べき $n^p$ を [初等関数の厳密な構成](/math/analysis/elementary-functions/) で定義した後は、同じ計算で「$p>1$ で収束、$p\le1$ で発散」がわかる。
また $1/\sqrt n\ge1/n$ なので、$\sum1/\sqrt n$ は比較判定法により発散する。
:::

### 比、べき根

:::theorem[比判定法]{#thm-anl-ratio-test}
十分後で $a_n\ne0$ とする。$r_n=\abs{a_{n+1}/a_n}$ に対し $\limsup r_n<1$ ならば $\sum\abs{a_n}$ は収束する。
十分後で $r_n\ge1$ ならば（特に $\liminf r_n>1$ ならば）$\sum a_n$ は発散する。
:::

:::proof
前者では $\limsup r_n<q<1$ をとると、:ref[thm-anl-limsup-epsilon] によりある $N$ 以降で $r_n<q$ となる。
従って $j\ge0$ で $\abs{a_{N+j}}\le\abs{a_N}q^j$ である。幾何級数との比較により結論を得る。
後者では、ある $N$ 以降で $\abs{a_{n+1}}\ge\abs{a_n}$ となるので $\abs{a_n}\ge\abs{a_N}>0$ である。一般項が $0$ に収束しないため級数は発散する。
$\liminf r_n>1$ なら、下極限の特徴付けにより十分後で $r_n>1$ である。
:::

次の判定法で用いる $\sqrt[n]{t}$ は、:ref[thm-anl-nth-root] で存在を示した非負の $n$ 乗根である。

:::theorem[べき根判定法]{#thm-anl-root-test}
$\rho=\limsup_{n\to\infty}\sqrt[n]{\abs{a_n}}$ とする。$\rho<1$ ならば $\sum\abs{a_n}$ は収束し、$\rho>1$ ならば $\sum a_n$ は発散する。
:::

:::proof
$\rho<1$ ならば $\rho<q<1$ を選ぶ。十分大きな $n$ で $\sqrt[n]{\abs{a_n}}<q$ となるから、$\abs{a_n}<q^n$ を比較判定法に用いればよい。
$\rho>1$ のときは $1<q<\rho$ を選ぶ（$\rho=+\infty$ なら $q=2$ とする）。:ref[thm-anl-limsup-epsilon] の 2（$\rho=+\infty$ なら上極限の定義）により、無限に多くの $n$ で $\sqrt[n]{\abs{a_n}}>q$ となる。
該当する項は $\abs{a_n}>q^n>1$ なので一般項は $0$ に収束しない。
:::

二つの判定法の強さは次の不等式で比較できる。

:::proposition[比とべき根の比較]{#prop-anl-ratio-root}
正の数列 $(c_n)$ について

$$
\liminf_{n\to\infty}\frac{c_{n+1}}{c_n}\le\liminf_{n\to\infty}\sqrt[n]{c_n}\le\limsup_{n\to\infty}\sqrt[n]{c_n}\le\limsup_{n\to\infty}\frac{c_{n+1}}{c_n}.
$$
:::

:::proof
中央は :ref[prop-anl-limsup-properties] の 1 である。右端を示す。$\beta=\limsup c_{n+1}/c_n$ とし、$\beta=+\infty$ なら明らかなので $\beta<\infty$ とする。
$\beta'>\beta$ を任意にとると、ある $N$ 以降で $c_{n+1}/c_n<\beta'$ となり、$n>N$ で $c_n\le c_N\beta'^{\,n-N}$ である。従って

$$
\sqrt[n]{c_n}\le\sqrt[n]{c_N\beta'^{-N}}\cdot\beta'.
$$

右辺は :ref[prop-anl-standard-limits] の 2 により $\beta'$ に収束するので、:ref[prop-anl-limsup-properties] の 3 により $\limsup\sqrt[n]{c_n}\le\beta'$ である。$\beta'>\beta$ は任意なので右端の不等式を得る。
左端は、$\alpha=\liminf c_{n+1}/c_n>0$ の場合に $0<\alpha'<\alpha$ をとって同様に $c_n\ge c_N\alpha'^{\,n-N}$ から従う（$\alpha=0$ なら明らか）。
:::

従って $\limsup r_n<1$ ならべき根判定法の $\rho<1$ も成り立ち、$\liminf r_n>1$ なら $\rho>1$ も成り立つ。比判定法のこの二つの形で判定できる級数は、べき根判定法でも判定できる。比の極限が存在すればべき根の極限も同じ値で存在する。
ただし「十分後で $r_n\ge1$」という発散の条件は別である。$a_n=1$ では $r_n=1$ により発散がわかるが、$\rho=1$ なのでべき根判定法では判定できない。
一方、べき根判定法の方が強い例は次の例と演習問題で見る。

:::example[判定法の選択]
$\sum n/2^n$ では隣り合う項の比が $(n+1)/(2n)\to1/2$ なので比判定法を使える。
これに対し、$a_n=2^{-n}$（$n$ が偶数）、$a_n=0$（$n$ が奇数）とすると、比の分母が何度も $0$ になる。
べき根の方は全項で定義できる。偶数番目では $\sqrt[n]{a_n}=1/2$、奇数番目では $0$ なので、その上極限は $1/2$ である。従ってこの級数も収束する。

判定値が $1$ になるときは別の議論が要る。例えば $1/n$ と $1/(n(n+1))$ はともに隣り合う項の比が $1$ に近づく。
前者の級数は発散するが、後者は $1/n-1/(n+1)$ と書けるため、第 $n$ 部分和が $1-1/(n+1)$ となって収束する。
:::

### Abel の変形と Dirichlet の判定法

符号が変わる級数では、部分和の打ち消し合いを利用する。そのための基本的な道具は、積分の部分積分に相当する次の恒等式である。

:::lemma[Abel の変形]{#lem-anl-abel-summation}
$A_k=\sum_{j=1}^ka_j$（$A_0=0$）とおく。$1\le m\le n$ について

$$
\sum_{k=m}^na_kb_k=A_nb_n-A_{m-1}b_m+\sum_{k=m}^{n-1}A_k(b_k-b_{k+1})
\label{eq-anl-abel-summation}
$$

が成り立つ。
:::

:::proof
$a_k=A_k-A_{k-1}$ を代入し、第二の和の添字を一つずらすと

$$
\sum_{k=m}^n(A_k-A_{k-1})b_k=\sum_{k=m}^nA_kb_k-\sum_{k=m-1}^{n-1}A_kb_{k+1}
=A_nb_n-A_{m-1}b_m+\sum_{k=m}^{n-1}A_k(b_k-b_{k+1})
$$

となる。
:::

:::theorem[Dirichlet の判定法]{#thm-anl-dirichlet-test}
$\sum a_n$ の部分和 $A_n$ が有界（$\abs{A_n}\le M$）で、$(b_n)$ が単調減少して $0$ に収束するならば、$\sum a_nb_n$ は収束する。さらに $m\le n$ について

$$
\left|\sum_{k=m}^na_kb_k\right|\le2Mb_m .
$$
:::

:::proof
$b_k\ge0$ かつ $b_k-b_{k+1}\ge0$ であるから、:eqref[eq-anl-abel-summation] の各項を絶対値で抑えて

$$
\left|\sum_{k=m}^na_kb_k\right|\le Mb_n+Mb_m+M\sum_{k=m}^{n-1}(b_k-b_{k+1})=Mb_n+Mb_m+M(b_m-b_n)=2Mb_m
$$

を得る。$b_m\to0$ なので、:eqref[eq-anl-series-cauchy] により級数は収束する。
:::

:::theorem[Abel の判定法]{#thm-anl-abel-test}
$\sum a_n$ が収束し、$(b_n)$ が単調かつ有界ならば、$\sum a_nb_n$ は収束する。
:::

:::proof
単調収束定理により $b_n\to b$ である。$b_n$ が減少なら $b_n-b$ は単調減少して $0$ に収束し、増加なら $b-b_n$ がそうである。
前者の場合 $a_nb_n=a_n(b_n-b)+ba_n$ と分ける。$\sum a_n$ の部分和は収束するので有界であり、:ref[thm-anl-dirichlet-test] により $\sum a_n(b_n-b)$ は収束する。$\sum ba_n$ も収束するので和も収束する。
後者の場合は $a_nb_n=-a_n(b-b_n)+ba_n$ と分ければ同様である。
:::

:::example[周期的な符号]
$\eps_n$ を $+1,+1,-1,-1$ の繰り返しとする。部分和 $\sum_{k=1}^n\eps_k$ は $1,2,1,0$ を繰り返すので $\abs{\cdot}\le2$ であり、$1/n$ は単調減少して $0$ に収束する。
Dirichlet の判定法により $\sum\eps_n/n=1+\tfrac12-\tfrac13-\tfrac14+\tfrac15+\cdots$ は収束する。絶対値をとると調和級数となるので発散する。
:::

### 交代級数

:::theorem[Leibniz の判定法]{#thm-anl-leibniz}
$b_n\ge0$ が単調減少して $0$ に収束するとき、:term[交代級数]{en="alternating series" yomi="こうたいきゅうすう"} $\sum_{n=1}^\infty(-1)^{n-1}b_n$ は収束する。
その和 $s$ と部分和 $s_n$ の間には $\abs{s-s_n}\le b_{n+1}$ が成り立つ。
:::

:::proof
$s_{2k+2}-s_{2k}=b_{2k+1}-b_{2k+2}\ge0$、$s_{2k+1}-s_{2k-1}=-b_{2k}+b_{2k+1}\le0$ である。
しかも $s_{2k}\le s_{2k}+b_{2k+1}=s_{2k+1}\le s_1$ かつ $s_{2k}\ge s_0=0$ である。よって偶数番目の部分和は上に有界な増加列、奇数番目の部分和は下に有界な減少列であり、それぞれ収束する。
$s_{2k+1}-s_{2k}=b_{2k+1}\to0$ なので、両者の極限は一致する。両方の列の誤差が $\eps$ 未満になる段階を選べば、元の部分和もそれ以降は誤差 $\eps$ 未満となる。

共通の極限 $s$ は各偶数部分和以上、各奇数部分和以下である。
従って $s$ は $s_n$ と $s_{n+1}$ の間にあり、この区間の長さが $b_{n+1}$ である。
:::

収束だけなら、部分和が $0$ と $1$ を交互にとる $\sum(-1)^{n-1}$ に :ref[thm-anl-dirichlet-test] を適用しても得られる。Leibniz の判定法の利点は誤差評価が具体的なことである。

例えば $1-1/2+1/4-1/8+\cdots$ は、公比 $-1/2$ の幾何級数で、和は $2/3$ である。
下図では、部分和が和をまたぐたびに、次の一項が残りの誤差の上界を与えている。

```tikz caption="交代する幾何級数の部分和。偶数番目の点は下から、奇数番目の点は上から同じ高さに近づく。"
\begin{tikzpicture}[x=0.8cm,y=4cm]
\draw[->] (0,0.42) -- (7.65,0.42) node[right] {$n$};
\draw[->] (0,0.42) -- (0,1.09) node[above] {$s_n$};
\draw[dashed] (0,{2/3}) -- (7.4,{2/3});
\foreach \y/\v in {1/1,0.5/{1/2},{2/3}/{2/3}}
  \node[left] at (0,\y) {$\v$};
\draw (1,1) -- (2,0.5) -- (3,0.75) -- (4,0.625)
  -- (5,0.6875) -- (6,0.65625) -- (7,0.671875);
\foreach \x/\y in {1/1,2/0.5,3/0.75,4/0.625,5/0.6875,6/0.65625,7/0.671875}{
  \fill (\x,\y) circle (1.5pt);
  \node[below] at (\x,0.42) {$\x$};
}
\end{tikzpicture}
```

## 絶対収束と条件収束・再配列

### 絶対値を付けても収束する場合

:::definition[絶対収束・条件収束・再配列]{#def-anl-absolute-convergence}
$\sum\abs{a_n}$ が収束するとき $\sum a_n$ は :term[絶対収束]{en="absolute convergence" yomi="ぜったいしゅうそく"} するという。
$\sum a_n$ が収束し $\sum\abs{a_n}$ が発散するとき :term[条件収束]{en="conditional convergence" yomi="じょうけんしゅうそく"} するという。
全単射 $\pi:\N\to\N$ によって得られる $\sum_{k=1}^\infty a_{\pi(k)}$ を :term[再配列]{en="rearrangement" yomi="さいはいれつ"} という。
:::

絶対収束から通常の収束が従うことは、部分和の Cauchy 条件で示した。その和については $\bigl|\sum a_n\bigr|\le\sum\abs{a_n}$ が成り立つ（有限和の三角不等式で極限をとる）。
一方、交代調和級数 $\sum(-1)^{n-1}/n$ は Leibniz の判定法で収束し、絶対値を付けた調和級数は発散する。従ってこれは条件収束の例である。

:::theorem[絶対収束する級数の再配列]{#thm-anl-absolute-rearrangement}
絶対収束する級数の任意の再配列は絶対収束し、元の級数と同じ和をもつ。
:::

:::proof
元の和を $s$、再配列の第 $m$ 部分和を $t_m$ とする。$\eps>0$ に対して $\sum_{n>K}\abs{a_n}<\eps$ となる $K$ を選ぶ（剰余が $0$ に収束するため）。
$M=\max\{\pi^{-1}(1),\dots,\pi^{-1}(K)\}$ とすれば $m\ge M$ で $t_m$ には $a_1,\dots,a_K$ がすべて含まれる。

$m\ge M$ を固定し、$J\ge\max\{\pi(1),\dots,\pi(m),K\}$ とする。
$s_J-t_m$ は、添字が $\{1,\dots,J\}$ に属し $\{\pi(1),\dots,\pi(m)\}$ に属さない項の和であり、それらの添字はすべて $K$ より大きい。従って

$$
\abs{s_J-t_m}\le\sum_{n=K+1}^J\abs{a_n}\le\sum_{n>K}\abs{a_n}<\eps .
$$

$J\to\infty$ とすれば $\abs{s-t_m}\le\sum_{n>K}\abs{a_n}<\eps$ となる。
再配列の絶対値の部分和は、元の絶対値の級数の部分和 $\sum_{n\le\max\pi(\{1,\dots,m\})}\abs{a_n}$ 以下なので有界であり、再配列も絶対収束する。
:::

### 正項と負項を別々に選ぶ

:::theorem[Riemann の再配列定理]{#thm-anl-riemann-rearrangement}
実数項の級数 $\sum a_n$ が条件収束するとする。任意に指定した $T\in\R$ に収束する再配列が存在する。
さらに $+\infty$、$-\infty$ に発散する再配列、および上にも下にも限りなく振動する再配列が存在する。
:::

:::proof
$a_n^+=\max\{a_n,0\},\ a_n^-=\max\{-a_n,0\}$ とおくと $a_n^\pm=(\abs{a_n}\pm a_n)/2$ である。
もし $\sum a_n^+$ が収束すれば $\sum\abs{a_n}=2\sum a_n^+-\sum a_n$ も収束してしまう。従って $\sum a_n^+=+\infty$ であり、同様に $\sum a_n^-=+\infty$ である。
また $a_n\to0$ である。従って正項と負項はそれぞれ無限にあり、元の順に正項を $p_1,p_2,\dots$、負項を $-q_1,-q_2,\dots$（$q_j>0$）と並べれば $p_i,q_j\to0$ である。$\sum p_i=\sum q_j=+\infty$ なので、有限個を除いた残りの和もやはり $+\infty$ である。

**有限の目標値。** 部分和 $0$ から出発し、次の段階を交互に行う。
上昇段階では、まだ使っていない正項を元の順に加え、部分和が初めて $T$ より大きくなったところで止める（最初の段階で既に $0>T$ なら何も加えない）。
下降段階では、まだ使っていない負項を元の順に加え、部分和が初めて $T$ より小さくなったところで止める。
各段階は、残りの正項・負項の和が $+\infty$ なので有限回で終わる。最初を除く各段階の開始時には部分和が目標の反対側にあるので、少なくとも一項を使う。
従ってすべての正項と負項はいずれ使われる。零の項があれば、上昇段階と下降段階の各組の後に、未使用の零項を元の順に一つずつ挿入する。こうして全項をちょうど一度ずつ使う再配列が得られる。

上昇段階の終わりの部分和 $\sigma$ は、最後に加えた正項を $p$ とすると、加える前は $T$ 以下だったので $T<\sigma\le T+p$ である。
同様に下降段階の終わりの部分和は、最後に加えた負項を $-q$ として $T-q\le\sigma<T$ を満たす。
上昇段階の途中の部分和は、その段階の開始値（直前の下降段階の終わりの値 $\ge T-q'$）から単調に増加して終わりの値（$\le T+p$）に至るので、すべて $[T-q',T+p]$ に入る。下降段階も同様である。零項の挿入は部分和を変えない。

$\eps>0$ を指定する。$i\ge I$ で $p_i<\eps$、$j\ge J$ で $q_j<\eps$ となる $I,J$ をとる。
$p_1,\dots,p_I$ と $q_1,\dots,q_J$ がすべて使われた段階より二段階後以降では、上の評価に現れる $p,q,q'$ はすべて $\eps$ 未満であり、部分和はすべて $(T-\eps,T+\eps)$ に入る。これで和 $T$ が得られた。

**無限大の目標値。** 第 $k$ 段階では、未使用の正項を少なくとも一項、部分和が $k+q_k$ を超えるまで加え、次に負項 $-q_k$ を一つ加える。
この直後の部分和は $k$ より大きい。次の段階で正項を加えている間はさらに増加し、その段階の終わりの負項を加えても $k+1$ より大きいので、第 $k$ 段階の終わり以降の部分和はすべて $k$ より大きい。従って部分和は $+\infty$ に発散する。
各段階で両符号の項を一つ以上使うので取り残しはなく、零項も先ほどと同じく挿入できる。級数全体の符号を替えてこの操作をし、符号を戻せば $-\infty$ に向かう再配列を得る。

**振動。** 第 $k$ 段階で未使用の正項を加えて部分和を $k$ より大きくし、次に未使用の負項を加えて $-k$ より小さくする。
同じ発散性により各段階は有限回で終わり、各段階で両符号の項を使うので全項を一度ずつ使える。部分和の上極限は $+\infty$、下極限は $-\infty$ である。
:::

この証明では各符号の項を元の順で使うことで全単射を作った。単に和を目標に近づける項だけ選ぶと、元の項をすべて使ったことにはならない。
数列と級数のこれらの定理については :cite[rudin-pma] も参照されたい。

### 級数の積

二つの級数 $\sum a_n$、$\sum b_n$ の積を考える。積 $a_jb_k$ をすべて足し合わせるには順序を決めなければならない。
べき級数の積と相性がよいのは、$j+k$ が一定の項をまとめる順序である。この項では添字を $0$ から始める。

:::definition[Cauchy 積]{#def-anl-cauchy-product}
級数 $\sum_{n=0}^\infty a_n$、$\sum_{n=0}^\infty b_n$ に対し

$$
c_n=\sum_{k=0}^na_kb_{n-k}=a_0b_n+a_1b_{n-1}+\cdots+a_nb_0
$$

とおき、$\sum_{n=0}^\infty c_n$ を二つの級数の :term[Cauchy 積]{en="Cauchy product" yomi="こーしーせき"} という。
:::

:::theorem[Mertens の定理]{#thm-anl-cauchy-product}
$\sum a_n=A$ が絶対収束し、$\sum b_n=B$ が収束するならば、Cauchy 積 $\sum c_n$ は $AB$ に収束する。
さらに両方が絶対収束するならば、Cauchy 積も絶対収束する。
:::

:::proof
部分和を $A_n,B_n,C_n$ とし、$\beta_n=B_n-B$ とおく。和の順序を入れ替えると

$$
C_n=\sum_{m=0}^n\sum_{k=0}^ma_kb_{m-k}=\sum_{k=0}^na_kB_{n-k}=A_nB+\gamma_n,\qquad
\gamma_n=\sum_{k=0}^na_k\beta_{n-k}
$$

となる。$A_nB\to AB$ なので、$\gamma_n\to0$ を示せばよい。

$\alpha=\sum\abs{a_n}$ とおき、$\eps>0$ をとる。$\beta_n\to0$ なので、$n\ge N$ で $\abs{\beta_n}\le\eps$ となる $N$ がある。$n\ge N$ について

$$
\abs{\gamma_n}\le\left|\sum_{i=0}^{N}\beta_ia_{n-i}\right|+\sum_{i=N+1}^n\abs{\beta_i}\abs{a_{n-i}}
\le\left|\sum_{i=0}^{N}\beta_ia_{n-i}\right|+\eps\alpha
$$

である（$i=n-k$ と添字を替えた）。$N$ を固定すると、第一項は $N+1$ 個の項の和で $a_n\to0$ より $0$ に収束する。
従って :ref[prop-anl-limsup-properties] の 3 により $\limsup\abs{\gamma_n}\le\eps\alpha$ であり、$\eps$ は任意なので $\limsup\abs{\gamma_n}\le0$ である。
$\abs{\gamma_n}\ge0$ より $\liminf\abs{\gamma_n}\ge0$ でもあるから、:ref[prop-anl-limsup-cluster] により $\abs{\gamma_n}\to0$、すなわち $\gamma_n\to0$ である。

両方が絶対収束するなら

$$
\sum_{n=0}^N\abs{c_n}\le\sum_{n=0}^N\sum_{k=0}^n\abs{a_k}\abs{b_{n-k}}\le\left(\sum_{k=0}^N\abs{a_k}\right)\left(\sum_{j=0}^N\abs{b_j}\right)\le\sum\abs{a_n}\cdot\sum\abs{b_n}
$$

であり、部分和が有界なので $\sum\abs{c_n}$ は収束する。
:::

:::example[絶対収束の仮定が必要な例]
$a_n=b_n=(-1)^n/\sqrt{n+1}$（$n\ge0$）とすると、Leibniz の判定法により両方とも収束する（条件収束）。Cauchy 積の項は

$$
c_n=(-1)^n\sum_{k=0}^n\frac1{\sqrt{(k+1)(n-k+1)}}
$$

である。相加平均と相乗平均の不等式から $(k+1)(n-k+1)\le\bigl((n+2)/2\bigr)^2$ なので、各項は $2/(n+2)$ 以上であり、$\abs{c_n}\ge2(n+1)/(n+2)\ge1$ となる。
一般項が $0$ に収束しないので $\sum c_n$ は発散する。
:::

:::example[幾何級数の積]
$\abs x<1$ のとき $\sum x^n=1/(1-x)$ は絶対収束する。自分自身との Cauchy 積は $c_n=\sum_{k=0}^nx^kx^{n-k}=(n+1)x^n$ であるから、Mertens の定理により

$$
\sum_{n=0}^\infty(n+1)x^n=\frac1{(1-x)^2}\qquad(\abs x<1)
$$

である。例えば $x=1/2$ とすると $\sum_{n=0}^\infty(n+1)/2^n=4$ である。
指数関数の加法定理 $\exp(x)\exp(y)=\exp(x+y)$ も、[初等関数の厳密な構成](/math/analysis/elementary-functions/) でこの定理から導かれる。
:::

## 演習問題

::::exercise{level=basic}
$a_n=(2n+5)/(n+2)$ の極限を求め、誤差 $1/100$ 未満を保証する最小の $N$ を求めよ。
:::solution
$a_n=2+1/(n+2)$ である。$\eps>0$ に対して $N=\lfloor1/\eps\rfloor+1$ とすれば、$n\ge N$ で $1/(n+2)<1/N<\eps$ なので極限は $2$ である。
指定された誤差については正確に $1/(n+2)<1/100$、すなわち $n+2>100$ を解けばよい。$n>98$ なので、$n\ge N$ のすべてで条件が成り立つ最小の $N$ は $99$ であり、$n=98$ では誤差がちょうど $1/100$ になる。
:::
::::

::::exercise{level=basic}
$x_1=0,\ x_{n+1}=(x_n+3)/2$ で定まる列について、単調性と有界性を示し、極限を求めよ。
:::solution
$x_1<3$ で、$x_n<3$ なら $x_{n+1}<(3+3)/2=3$ なので、帰納法で全項が $3$ 未満である。
従って $x_{n+1}-x_n=(3-x_n)/2>0$ となり単調増加である。単調収束定理により極限 $L$ が存在する。
漸化式で極限をとると $L=(L+3)/2$、従って $L=3$ である。
実際 $3-x_{n+1}=(3-x_n)/2$ から誤差は $3-x_n=3/2^{n-1}$ ともわかる。
:::
::::

::::exercise{level=basic}
$a_{2k}=2-1/k,\ a_{2k-1}=-1+1/k$ とする。上極限と下極限を求めよ。
:::solution
偶数番目の部分列は $2$ に、奇数番目の部分列は $-1$ に収束するので、$2$ と $-1$ は集積値である。
一方、すべての項は $[-1,2]$ に入るので、どの集積値もこの区間に入る（:ref[prop-anl-squeeze]）。
:ref[prop-anl-limsup-cluster] より、上極限は集積値の最大値 $2$、下極限は最小値 $-1$ である。
:::
::::

::::exercise{level=basic}
$\sum_{n=1}^\infty 1/(n(n+2))$ の和と、第 $N$ 部分和を使うときの誤差を求めよ。
:::solution
$1/(n(n+2))=(1/n-1/(n+2))/2$ である。これを $n=1,\dots,N$ について足すと、初めの二項と最後の二項を残して打ち消し合い、

$$
s_N=\frac12\left(1+\frac12-\frac1{N+1}-\frac1{N+2}\right)=\frac34-\frac{1}{2(N+1)}-\frac{1}{2(N+2)}
$$

を得る。従って和は $3/4$、誤差は $3/4-s_N=1/(2(N+1))+1/(2(N+2))$ である。
:::
::::

::::exercise{level=basic}
$\sum_{n=1}^\infty n/2^n$ の和を求めよ。
:::solution
本文の幾何級数の積の例で $x=1/2$ とすると $\sum_{n=0}^\infty(n+1)/2^n=4$ である。添字を $m=n+1$ と替えると

$$
\sum_{m=1}^\infty\frac m{2^{m-1}}=4,\qquad\text{従って}\qquad\sum_{m=1}^\infty\frac m{2^m}=\frac12\cdot4=2
$$

である（:ref[prop-anl-series-basics] の 1）。
:::
::::

::::exercise{level=advanced}
非負項の級数 $\sum b_n$ が収束し、$(b_n)$ が単調減少ならば $nb_n\to0$ であることを示せ。
:::hint
$\lfloor n/2\rfloor$ より後の項だけを足して $b_n$ と比較する。
:::
:::solution
$n\ge2$ について $m=\lfloor n/2\rfloor$ とすると $n-m\ge n/2$ であり、$m<k\le n$ で $b_k\ge b_n$ となるから

$$
0\le nb_n\le2(n-m)b_n\le2\sum_{k=m+1}^n b_k .
$$

$\eps>0$ に対して部分和の Cauchy 条件で、$n>m\ge M$ なら右の和が $\eps/2$ 未満になる添字 $M$ を選ぶ。
$n\ge2M$ なら $m\ge M$ かつ $n>m$ なので $nb_n<\eps$ となる。
:::
::::

::::exercise{level=advanced}
交代調和級数の和を $s$ とする。正項を二つ、負項を一つずつ元の順に並べた再配列 $1+\tfrac13-\tfrac12+\tfrac15+\tfrac17-\tfrac14+\cdots$ の和を $s$ で表せ。
:::solution
$H_m=\sum_{j=1}^m1/j$ とおく。元の級数の第 $2k$ 部分和は、奇数の逆数の和から偶数の逆数の和を引いたもので

$$
\sum_{j=1}^{2k}\frac{(-1)^{j-1}}j=H_{2k}-2\cdot\frac12H_k=H_{2k}-H_k
$$

となるので $H_{2k}-H_k\to s$ である。
再配列の第 $3k$ 部分和では、正項を $2k$ 個（分母 $1,3,\dots,4k-1$）、負項を $k$ 個（分母 $2,4,\dots,2k$）使っている。分母が $4k$ 以下の奇数の逆数の和は $H_{4k}-\tfrac12H_{2k}$ であるから

$$
t_{3k}=H_{4k}-\tfrac12H_{2k}-\tfrac12H_k
=(H_{4k}-H_{2k})+\tfrac12(H_{2k}-H_k)\to s+\tfrac12s=\tfrac32s .
$$

次の二つの部分和は $t_{3k+1}=t_{3k}+1/(4k+1)$、$t_{3k+2}=t_{3k+1}+1/(4k+3)$ であり、加えた項は $0$ に近づく。
従って全体の和も $3s/2$ である。しかも Leibniz の判定法の証明から $s\ge s_2=1-1/2>0$ なので、和は実際に変化した。
:::
::::

::::exercise{level=advanced}
$a_n\to a$ ならば $\sigma_n=(a_1+\cdots+a_n)/n\to a$ であることを示せ（Cesàro 平均）。また逆が成り立たない例を挙げよ。
:::solution
$\eps>0$ に対し、$n\ge N$ で $\abs{a_n-a}<\eps/2$ となる $N$ をとり、$C=\sum_{k=1}^{N}\abs{a_k-a}$ とおく。$n>N$ なら

$$
\abs{\sigma_n-a}=\left|\frac1n\sum_{k=1}^n(a_k-a)\right|\le\frac Cn+\frac{n-N}n\cdot\frac\eps2<\frac Cn+\frac\eps2
$$

である。$n>2C/\eps$ も満たせば右辺は $\eps$ 未満となるので、$\sigma_n\to a$ である。

逆の反例として $a_n=(-1)^n$ をとると、$a_1+\cdots+a_n$ は $-1$ と $0$ を交互にとるので $\abs{\sigma_n}\le1/n\to0$ であるが、$(a_n)$ は収束しない。
:::
::::

::::exercise{level=advanced}
$a_n=2^{-n+(-1)^n}$ とする。比判定法では $\sum a_n$ の収束が判定できず、べき根判定法では判定できることを示せ。
:::solution
$a_{n+1}/a_n=2^{-1+(-1)^{n+1}-(-1)^n}=2^{-1-2(-1)^n}$ は、$n$ が偶数なら $2^{-3}=1/8$、奇数なら $2^{1}=2$ である。
従って $\limsup a_{n+1}/a_n=2>1$、$\liminf a_{n+1}/a_n=1/8<1$ であり、:ref[thm-anl-ratio-test] のどちらの条件も満たさない。

一方 $\sqrt[n]{a_n}=\tfrac12\cdot\sqrt[n]{2^{(-1)^n}}$ であり、$n$ が偶数なら $\tfrac12\sqrt[n]2$、奇数なら $\tfrac12/\sqrt[n]2$ である。
:ref[prop-anl-standard-limits] の 2 により $\sqrt[n]2\to1$ なので、どちらの部分列も $1/2$ に収束し、$\sqrt[n]{a_n}\to1/2<1$ である。べき根判定法により級数は収束する。
これは :ref[prop-anl-ratio-root] の両端の不等式が狭義になりうることも示している（$1/8<1/2=1/2<2$）。
:::
::::

::::exercise{level=advanced}
$a_n\ge0$ で $\sum a_n$ が収束するならば、$\sum\sqrt{a_n}/n$ も収束することを示せ。
:::hint
相加平均と相乗平均の不等式を用いる。
:::
:::solution
$x,y\ge0$ について $\sqrt{xy}\le(x+y)/2$ である（[実数の連続性](/math/analysis/real-numbers/) の演習問題）。$x=a_n$、$y=1/n^2$ とすると

$$
0\le\frac{\sqrt{a_n}}n=\sqrt{a_n\cdot\frac1{n^2}}\le\frac12\left(a_n+\frac1{n^2}\right)
$$

である。右辺は収束する二つの級数 $\sum a_n$ と $\sum1/n^2$（凝集判定法の例で $p=2$）の和の半分なので収束し、比較判定法により左辺の級数も収束する。
:::
::::

## 確認テスト

::::quiz
実数列について、収束を保証する条件をすべて選べ。

- [ ] 有界である。
- [x] 上に有界で単調増加である。
- [x] Cauchy 列である。
- [ ] 収束部分列が一つある。

:::explanation
(B) は :ref[thm-anl-monotone-convergence]、(C) は :ref[thm-anl-cauchy-complete] による。
$(-1)^n$ は有界で収束部分列（偶数番目の項は定数 $1$）ももつが、二つの集積値 $1,-1$ をもつため収束しない。従って (A) と (D) は収束を保証しない。
:::
::::

::::quiz{answer="4"}
$\sum_{n=0}^\infty (3/4)^n$ の和を求めよ。

:::explanation
初項が $1$、公比が $3/4$ の幾何級数であり、$\abs{3/4}<1$ なので収束する。部分和は $(1-(3/4)^n)/(1-3/4)$ であり、和は $1/(1-3/4)=4$ である。
:::
::::

::::quiz{answer="99"}
交代調和級数の和を第 $N$ 部分和で近似する。Leibniz の誤差評価 $1/(N+1)$ が $1/100$ **以下**となる最小の $N$ を求めよ。

:::explanation
:ref[thm-anl-leibniz] の評価は $\abs{s-s_N}\le b_{N+1}=1/(N+1)$ である。$1/(N+1)\le1/100$ は $N+1\ge100$ と同値なので、最小値は $99$ である。
:::
::::

::::quiz
$a_n\ne0$ で $\lim\abs{a_{n+1}/a_n}=1$ のとき、必ず正しいものはどれか。

- [ ] 級数は絶対収束する。
- [ ] 級数は発散する。
- [x] この比の極限だけでは級数の収束は決まらない。

:::explanation
$a_n=1/n$ と $a_n=1/(n(n+1))$ はともに比の極限が $1$ である。調和級数は発散し、後者の部分和 $1-1/(n+1)$ は $1$ に収束する。
従って比の極限が $1$ という情報だけからは、どちらの結論も導けない。:ref[prop-anl-ratio-root] により、べき根判定法でもこの場合は判定できない。
:::
::::

::::quiz
条件収束する実数項級数について正しいものをすべて選べ。

- [x] 正項だけを足した部分和は上に有界でない。
- [ ] すべての再配列が同じ和をもつ。
- [x] 和が $0$ になる再配列がある。

:::explanation
(A)：正項の部分和は $\sum a_n^+$ の部分和であり、もし有界なら $\sum\abs{a_n}=2\sum a_n^+-\sum a_n$ が収束して条件収束であることに反する。
(C)：:ref[thm-anl-riemann-rearrangement] では任意の実数を指定できるから、$0$ も他の値も再配列の和になる。従って (B) は誤りである。
:::
::::

::::quiz{answer="9"}
$s_n=\sum_{k=0}^n1/k!$ で $e$ を近似する。:ref[prop-anl-e-series] の誤差の上界 $1/(n!\,n)$ が $10^{-6}$ 未満となる最小の $n$ を求めよ。

:::explanation
$n!\,n$ は $n$ について増加し、$n=8$ では $8!\cdot8=322560<10^6$、$n=9$ では $9!\cdot9=3265920>10^6$ である。
従って $1/(n!\,n)<10^{-6}$ となる最小の $n$ は $9$ である。実際 $s_9=2.7182815\cdots$ であり、$e=2.7182818\cdots$ との差は約 $3\times10^{-7}$ である。
:::
::::

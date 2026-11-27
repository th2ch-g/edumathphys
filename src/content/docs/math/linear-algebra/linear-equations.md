---
title: 連立一次方程式と行列
description: 行列の演算、行基本変形と簡約行階段形、連立一次方程式の解法、階数、正則行列と逆行列
notation:
  - symbol: 'M_{m,n}(K)'
    meaning: 体 K の元を成分とする m×n 行列全体
  - symbol: 'M_n(K)'
    meaning: K 上の n 次正方行列全体
  - symbol: 'I_n'
    meaning: n 次単位行列
  - symbol: '\transp{A}'
    meaning: 行列 A の転置
  - symbol: '\rank A'
    meaning: 行列 A の階数
  - symbol: '(A\mid b)'
    meaning: 係数行列 A に右辺 b を付け加えた拡大係数行列
  - symbol: 'E_{ij}'
    meaning: 指定したサイズで第 (i,j) 成分だけが 1 の行列
  - symbol: 'A^{-1}'
    meaning: 正則行列 A の逆行列
---

## この章の目標

- 行列の和・スカラー倍・積を定義し、その代数的性質を理解する。
- 行基本変形によって任意の行列を簡約行階段形に変形でき、しかもその結果が一意であることを証明する。
- 連立一次方程式の解の存在と自由度を、行列の階数によって完全に記述する。
- 正則行列のさまざまな同値な特徴付けを証明し、逆行列を計算できるようになる。

この章では、係数を主に $K=\R$（実数体）または $K=\C$（複素数体）にとる。計算には加減乗除だけを使うため、以下の結果は、仮定を明記した箇所を除いて任意の**体** $K$ に対して成り立つ。体の公理は次章 [ベクトル空間・基底・次元](/math/linear-algebra/vector-spaces/) で述べる。

## 行列とその演算

### 行列の定義

複数の一次式を扱うときは、各式で変数を同じ順に並べると、変数名を繰り返さずに係数だけで式を記録できる。例えば $x_1-2x_2+3x_3$ と $4x_1+x_3$ の係数を並べると、次の表が得られる。

$$
\begin{pmatrix}1&-2&3\\4&0&1\end{pmatrix}
$$

行番号は式を、列番号は変数を指定する。第 2 式の $x_2$ の係数が $0$ であることも、表には明示する。

:::definition[行列]{#def-linalg-matrix}
$m, n$ を正の整数とする。$K$ の元 $a_{ij}$（$1\le i\le m$, $1\le j\le n$）を縦 $m$ 個、横 $n$ 個の長方形に並べたもの

$$
A = \begin{pmatrix}
a_{11} & a_{12} & \cdots & a_{1n}\\
a_{21} & a_{22} & \cdots & a_{2n}\\
\vdots & \vdots & \ddots & \vdots\\
a_{m1} & a_{m2} & \cdots & a_{mn}
\end{pmatrix}
$$

を $K$ 上の $m\times n$ :term[行列]{en="matrix" yomi="ぎょうれつ"} という。$a_{ij}$ を $A$ の $(i,j)$ :term[成分]{en="entry" yomi="せいぶん"} と呼び、$A=(a_{ij})$ と略記する。$K$ 上の $m\times n$ 行列全体の集合を $M_{m,n}(K)$ で表す。$m=n$ のとき $A$ を $n$ 次 :term[正方行列]{en="square matrix" yomi="せいほうぎょうれつ"} といい、$M_n(K)\defeq M_{n,n}(K)$ と書く。
:::

横に並んだ $(a_{i1},\dots,a_{in})$ を $A$ の第 $i$ **行**、縦に並んだ成分を第 $j$ **列**という。$m\times 1$ 行列を :term[列ベクトル]{en="column vector" yomi="れつべくとる"}、$1\times n$ 行列を :term[行ベクトル]{en="row vector" yomi="ぎょうべくとる"} と呼ぶ。本サイトでは $K^n$ を $n$ 次列ベクトル全体の集合とし、その元を

$$
x = \begin{pmatrix} x_1\\ \vdots\\ x_n\end{pmatrix}
$$

と書く。2 つの行列 $A=(a_{ij})$, $B=(b_{ij})$ は、サイズが等しく、かつすべての $i,j$ について $a_{ij}=b_{ij}$ であるときに限り**等しい**といい、$A=B$ と書く。

:::definition[転置行列]{#def-linalg-transpose}
$A=(a_{ij})\in M_{m,n}(K)$ に対し、$(j,i)$ 成分が $a_{ij}$ である $n\times m$ 行列を $A$ の :term[転置行列]{en="transpose" yomi="てんちぎょうれつ"} といい、$\transp{A}$ と書く。例えば $\transp{(x_1,\dots,x_n)}$ は行 $(x_1,\dots,x_n)$ を転置した列ベクトルである。
:::

すべての成分が $0$ である行列を :term[零行列]{en="zero matrix" yomi="れいぎょうれつ"} といい $O$（サイズを明示するときは $O_{m,n}$）で表す。正方行列 $A=(a_{ij})$ の成分 $a_{11},a_{22},\dots,a_{nn}$ を**対角成分**という。対角成分がすべて $1$ でそれ以外の成分がすべて $0$ である $n$ 次正方行列を $n$ 次 :term[単位行列]{en="identity matrix" yomi="たんいぎょうれつ"} といい $I_n$ と書く。Kronecker のデルタ

$$
\delta_{ij} = \begin{cases} 1 & (i=j)\\ 0 & (i\ne j)\end{cases}
$$

を用いれば $I_n=(\delta_{ij})$ である。$I_n$ の第 $j$ 列を $e_j$ と書き、**基本ベクトル**と呼ぶ。

### 特別な形の行列

成分の多くが $0$ と分かっている行列では、計算に必要な項を先に絞れる。正方行列について、$i>j$ で $a_{ij}=0$ であるものを :term[上三角行列]{en="upper triangular matrix" yomi="うえさんかくぎょうれつ"}、$i<j$ で $a_{ij}=0$ であるものを :term[下三角行列]{en="lower triangular matrix" yomi="したさんかくぎょうれつ"} という。上下の両方の条件を満たすものが :term[対角行列]{en="diagonal matrix" yomi="たいかくぎょうれつ"} であり、対角成分を並べて $\diag(d_1,\dots,d_n)$ と書く。

:::example[成分と行列の形]
$\begin{pmatrix}2&-1&3\\0&0&5\\0&0&4\end{pmatrix}$ は上三角行列である。対角成分が $0$ でも、対角線より下がすべて $0$ ならば上三角行列である。$\diag(2,0,4)$ は上三角行列でも下三角行列でもあり、$I_3=\diag(1,1,1)$ は単位行列である。
:::

指定したサイズの行列で、$(i,j)$ 成分だけが $1$ で他の成分が $0$ のものを :term[行列単位]{en="matrix unit" yomi="ぎょうれつたんい"} $E_{ij}$ と呼ぶ。例えば $2\times 3$ 行列の $E_{23}$ は、第 1 行が零行で第 2 行が $(0,0,1)$ である。

### 和・スカラー倍・積

:::definition[行列の和とスカラー倍]{#def-linalg-matrix-sum}
$A=(a_{ij}), B=(b_{ij})\in M_{m,n}(K)$ と $c\in K$ に対し、和 $A+B$ とスカラー倍 $cA$ を

$$
A+B \defeq (a_{ij}+b_{ij}),\qquad cA \defeq (c\,a_{ij})
$$

で定める。また $-A\defeq(-1)A$, $A-B\defeq A+(-B)$ とする。
:::

和とスカラー倍の法則は、行列の等しさを成分ごとに判定すれば証明できる。例えば $(A+B)+C$ と $A+(B+C)$ の $(i,j)$ 成分はそれぞれ $(a_{ij}+b_{ij})+c_{ij}$ と $a_{ij}+(b_{ij}+c_{ij})$ であり、$K$ の加法の結合法則で一致する。各成分で次の等式を用いると、残りの行列の等式が得られる。

$$
\begin{aligned}
a_{ij}+b_{ij}&=b_{ij}+a_{ij},&
a_{ij}+0&=a_{ij},&
a_{ij}+(-a_{ij})&=0,\\
s(a_{ij}+b_{ij})&=sa_{ij}+sb_{ij},&
(s+t)a_{ij}&=sa_{ij}+ta_{ij},\\
(st)a_{ij}&=s(ta_{ij}),&
1a_{ij}&=a_{ij}.
\end{aligned}
$$

成分ごとに計算するため、任意の $A\in M_{m,n}(K)$ は $\sum_{i=1}^m\sum_{j=1}^n a_{ij}E_{ij}$ と書ける。右辺の $(i,j)$ 成分に寄与する項は $a_{ij}E_{ij}$ だけである。

:::example[二つの記録の差]
$A=\begin{pmatrix}1&-2&3\\4&0&1\end{pmatrix}$、$B=\begin{pmatrix}0&1&2\\-1&3&0\end{pmatrix}$ とする。例えば第 $(2,1)$ 成分は $2\cdot4-(-1)=9$ なので、$2A-B=\begin{pmatrix}2&-5&4\\9&-3&2\end{pmatrix}$ である。$B$ のサイズが $3\times2$ ならば、同じ成分数でも $A$ と $B$ の差は定義されない。
:::

:::definition[行列の積]{#def-linalg-matrix-product}
$A=(a_{ij})\in M_{l,m}(K)$, $B=(b_{jk})\in M_{m,n}(K)$ に対し、積 $AB\in M_{l,n}(K)$ を

$$
(AB)_{ik} \defeq \sum_{j=1}^{m} a_{ij}b_{jk}\qquad(1\le i\le l,\ 1\le k\le n) \label{eq-linalg-product}
$$

で定める。積 $AB$ は、$A$ の列の個数と $B$ の行の個数が等しいときに限り定義される。
:::

式 :eqref[eq-linalg-product] は「$A$ の第 $i$ 行と $B$ の第 $k$ 列の対応する成分の積の和」である。ここで、各ベクトルを定数倍して足したものを**一次結合**という。特に $A\in M_{m,n}(K)$ と $x\in K^n$ の積 $Ax\in K^m$ は

$$
Ax = \begin{pmatrix} a_{11}x_1+\cdots+a_{1n}x_n\\ \vdots\\ a_{m1}x_1+\cdots+a_{mn}x_n\end{pmatrix}
= x_1 a_1 + x_2 a_2 + \cdots + x_n a_n
$$

と書ける。ここで $a_j$ は $A$ の第 $j$ 列である。すなわち **$Ax$ は $A$ の列ベクトルの $x_j$ を係数とする一次結合**である。同様に、$AB$ の第 $k$ 列は $A$ に $B$ の第 $k$ 列を掛けたもの $A b_k$ に等しい。

:::proposition[積の性質]{#prop-linalg-product-laws}
行列のサイズは、以下の演算が定義される範囲で任意とする。$c\in K$ とする。

1. （結合法則）$(AB)C = A(BC)$
2. （分配法則）$A(B+C) = AB+AC$, $(A+B)C = AC+BC$
3. $(cA)B = A(cB) = c(AB)$
4. $A\in M_{m,n}(K)$ ならば $I_mA = A = AI_n$
:::

:::proof
1. $A\in M_{k,l}(K)$, $B\in M_{l,m}(K)$, $C\in M_{m,n}(K)$ とする。$(i,q)$ 成分を比較すると

$$
((AB)C)_{iq} = \sum_{r=1}^{m}\Bigl(\sum_{j=1}^{l} a_{ij}b_{jr}\Bigr)c_{rq}
= \sum_{j=1}^{l}\sum_{r=1}^{m} a_{ij}b_{jr}c_{rq}
= \sum_{j=1}^{l} a_{ij}\Bigl(\sum_{r=1}^{m} b_{jr}c_{rq}\Bigr) = (A(BC))_{iq}.
$$

有限和の順序交換と、体における分配法則・結合法則を用いた。

2. $(A(B+C))_{ik}=\sum_j a_{ij}(b_{jk}+c_{jk})=\sum_j a_{ij}b_{jk}+\sum_j a_{ij}c_{jk}=(AB+AC)_{ik}$ である。また、$((A+B)C)_{ik}=\sum_j(a_{ij}+b_{ij})c_{jk}=(AC+BC)_{ik}$ である。
3. $((cA)B)_{ik}=\sum_j(ca_{ij})b_{jk}=c\sum_j a_{ij}b_{jk}$ であり、$(A(cB))_{ik}=\sum_j a_{ij}(cb_{jk})$ も $K$ の乗法の交換法則と結合法則により同じ値である。
4. $(I_mA)_{ik}=\sum_j\delta_{ij}a_{jk}=a_{ik}$、$(AI_n)_{ik}=\sum_j a_{ij}\delta_{jk}=a_{ik}$ である。
:::

積の定義から $Ae_j$ は $A$ の第 $j$ 列、$\transp{e_i}A$ は $A$ の第 $i$ 行である。さらに $E_{ij}E_{kl}$ の成分の計算では、第 1 因子の第 $j$ 列と第 2 因子の第 $k$ 行が一致する添字でだけ積が残る。したがって、サイズが合うとき $E_{ij}E_{kl}=\delta_{jk}E_{il}$ である。

:::example[二段階の計算]
$x=\transp{(x_1,x_2,x_3)}$ から $y=\transp{(x_1+2x_2,\ x_2-x_3)}$ を作り、次に $z=\transp{(2y_1+y_2,\ y_1-y_2)}$ を作る。二つの計算を表す行列をそれぞれ $B,A$ とすれば、$z=A(Bx)=(AB)x$ であり、

$$
AB=\begin{pmatrix}2&1\\1&-1\end{pmatrix}
\begin{pmatrix}1&2&0\\0&1&-1\end{pmatrix}
=\begin{pmatrix}2&5&-1\\1&1&1\end{pmatrix}.
$$

$x=\transp{(1,0,2)}$ のとき、$y=\transp{(1,-2)}$、$z=\transp{(0,3)}$ であり、まとめた行列でも $(AB)x=\transp{(0,3)}$ となる。右側の因子が先に $x$ に作用するため、この例では $B$ が先である。
:::

:::proposition[三角行列の積]{#prop-linalg-triangular-product}
同じサイズの上三角行列の積は上三角行列であり、その対角成分は元の対応する対角成分の積である。下三角行列どうしにも同じ性質が成り立つ。
:::

:::proof
$A,B$ が上三角行列なら、積の第 $(i,j)$ 成分で $a_{ik}b_{kj}$ が $0$ でないためには $i\le k\le j$ が必要である。$i>j$ なら該当する $k$ がなく、成分は $0$ である。$i=j$ なら $k=i$ の項だけが残り、$(AB)_{ii}=a_{ii}b_{ii}$ となる。下三角の場合は条件が $i\ge k\ge j$ となり、$i<j$ で $0$、$i=j$ で $a_{ii}b_{ii}$ となる。
:::

:::remark
行列の積は一般に**可換でない**。例えば

$$
P=\begin{pmatrix}1&2\\3&4\end{pmatrix},\quad Q=\begin{pmatrix}0&1\\1&0\end{pmatrix}
\quad\Longrightarrow\quad
PQ=\begin{pmatrix}2&1\\4&3\end{pmatrix}\ne\begin{pmatrix}3&4\\1&2\end{pmatrix}=QP.
$$

また、$A\ne O$, $B\ne O$ でも $AB=O$ となることがある（例：$A=B=\begin{pmatrix}0&1\\0&0\end{pmatrix}$）。したがって、数の場合の「$ab=0$ ならば $a=0$ または $b=0$」は行列では成り立たない。
:::

### 転置とブロック分割

:::proposition{#prop-linalg-transpose-laws}
行列のサイズは以下の演算が定義されるものとし、$c\in K$ とする。このとき $(A+B)^\top = \transp{A}+\transp{B}$、$(cA)^\top = c\,\transp{A}$、$(\transp{A})^\top = A$、$(AB)^\top = \transp{B}\,\transp{A}$ が成り立つ。
:::

:::proof
$A,B\in M_{m,n}(K)$ なら、転置の定義により

$$
\begin{aligned}
((A+B)^\top)_{ji}&=a_{ij}+b_{ij}=(\transp{A}+\transp{B})_{ji},\\
((cA)^\top)_{ji}&=ca_{ij}=(c\transp{A})_{ji},\\
((\transp{A})^\top)_{ij}&=(\transp{A})_{ji}=a_{ij}.
\end{aligned}
$$

積については $A\in M_{l,m}(K)$, $B\in M_{m,n}(K)$ とする。このとき $((AB)^\top)_{ki}=(AB)_{ik}=\sum_j a_{ij}b_{jk}=\sum_j(\transp{B})_{kj}(\transp{A})_{ji}=(\transp{B}\transp{A})_{ki}$ である。
:::

行列を縦横の区切り線で小さな行列（**ブロック**）に分けて扱うと便利なことが多い。例えば $A=\begin{pmatrix}A_{11}&A_{12}\\A_{21}&A_{22}\end{pmatrix}$, $B=\begin{pmatrix}B_{11}&B_{12}\\B_{21}&B_{22}\end{pmatrix}$ と分割し、各積 $A_{ij}B_{jk}$ が定義されるように区切り方が整合していれば、

$$
AB = \begin{pmatrix}A_{11}B_{11}+A_{12}B_{21} & A_{11}B_{12}+A_{12}B_{22}\\ A_{21}B_{11}+A_{22}B_{21} & A_{21}B_{12}+A_{22}B_{22}\end{pmatrix}
$$

が成り立つ。実際、$A$ の列と $B$ の行を同じ位置で区切り、第 1 区間を $J_1$、第 2 区間を $J_2$ とする。第 $(p,q)$ 成分の和を $\sum_{j\in J_1}a_{pj}b_{jq}+\sum_{j\in J_2}a_{pj}b_{jq}$ と分けると、$p$ が第 $i$ 行ブロック、$q$ が第 $k$ 列ブロックにあるとき、これは $A_{i1}B_{1k}+A_{i2}B_{2k}$ の対応する成分である。区間を増やした場合も、区間ごとに和をとれば同じ規則が成り立つ。

:::example[幅の異なるブロック]
$D=\begin{pmatrix}1&2&3\\0&1&4\end{pmatrix}$ を最初の 2 列と最後の 1 列に分け、$F=\begin{pmatrix}1&0\\2&1\\-1&2\end{pmatrix}$ を最初の 2 行と最後の 1 行に分ける。対応する積を足すと

$$
\begin{aligned}
DF&=\begin{pmatrix}1&2\\0&1\end{pmatrix}\begin{pmatrix}1&0\\2&1\end{pmatrix}
+\begin{pmatrix}3\\4\end{pmatrix}\begin{pmatrix}-1&2\end{pmatrix}\\
&=\begin{pmatrix}5&2\\2&1\end{pmatrix}
+\begin{pmatrix}-3&6\\-4&8\end{pmatrix}
=\begin{pmatrix}2&8\\-2&9\end{pmatrix}.
\end{aligned}
$$

第 2 項は列ベクトルと行ベクトルの積なので $2\times2$ 行列である。行ベクトルと列ベクトルを逆順に掛けると $1\times1$ 行列になる。
:::

複素行列では、転置と成分ごとの複素共役をともに取ったものを $\adj{A}=(\conj{a_{ji}})$ と書き、**共役転置**という。例えば $A=\begin{pmatrix}1&i\\2-i&0\end{pmatrix}$ なら $\adj{A}=\begin{pmatrix}1&2+i\\-i&0\end{pmatrix}$ である。和と積の複素共役は成分ごとの和と積の共役に等しいので、:ref[prop-linalg-transpose-laws] と合わせると $(AB)^*=\adj{B}\adj{A}$、$(cA)^*=\conj{c}\adj{A}$、$(\adj{A})^*=A$ が従う。複素内積との関係は [内積空間](/math/linear-algebra/inner-product-spaces/) で扱う。

### 連立一次方程式の行列表示

$n$ 個の未知数 $x_1,\dots,x_n$ に関する $m$ 本の方程式

$$
\begin{cases}
a_{11}x_1 + \cdots + a_{1n}x_n = b_1\\
\qquad\vdots\\
a_{m1}x_1 + \cdots + a_{mn}x_n = b_m
\end{cases}
$$

を :term[連立一次方程式]{en="system of linear equations" yomi="れんりついちじほうていしき"} という。:term[係数行列]{en="coefficient matrix" yomi="けいすうぎょうれつ"} を $A=(a_{ij})\in M_{m,n}(K)$、右辺の列を $b=\transp{(b_1,\dots,b_m)}$ とおけば、この方程式は

$$
Ax = b \label{eq-linalg-system}
$$

と書ける。$b=0$ のとき :term[斉次]{en="homogeneous" yomi="せいじ"}（同次）であるという。$A$ の右に $b$ を付け加えた $m\times(n+1)$ 行列 $(A\mid b)$ を :term[拡大係数行列]{en="augmented matrix" yomi="かくだいけいすうぎょうれつ"} と呼ぶ。

## 基本変形と行階段形

### 行基本変形と基本行列

連立一次方程式を解く消去法の操作を、行列の言葉で定式化する。

:::definition[行基本変形]{#def-linalg-row-operations}
行列に対する次の 3 種類の操作を :term[行基本変形]{en="elementary row operation" yomi="ぎょうきほんへんけい"} という。

1. 第 $i$ 行と第 $j$ 行を入れ替える（$i\ne j$）。
2. 第 $i$ 行を $c$ 倍する（$c\in K$, $c\ne 0$）。
3. 第 $i$ 行に第 $j$ 行の $c$ 倍を加える（$i\ne j$, $c\in K$）。

行列 $A$ に有限回の行基本変形を施して $B$ が得られるとき、$A$ と $B$ は :term[行同値]{en="row equivalent" yomi="ぎょうどうち"} であるという。変形を 0 回とすれば $A$ は自分自身と行同値である。各変形を逆順に取り消せば $B$ から $A$ に戻り、変形の列をつなげば推移性が成り立つ。したがって行同値は同値関係である。
:::

それぞれの行基本変形は、単位行列に同じ変形を施した行列を左から掛けることで実現できる。ここで、正方行列 $E$ に対して $EE'=E'E=I$ を満たす行列 $E'$ があるとき、$E$ は**正則**であるといい、$E'$ を逆行列 $E^{-1}$ と書く。逆行列の一意性と積の扱いは :ref[prop-linalg-inverse-basic] で証明する。以下、基本行列の積に対応する逆操作は、各行操作を逆順に取り消して得られる。

:::definition[基本行列]{#def-linalg-elementary-matrix}
$I_m$ に上の 1, 2, 3 の変形を施した行列をそれぞれ $P_{ij}$, $Q_i(c)$, $R_{ij}(c)$ と書き、これらを総称して :term[基本行列]{en="elementary matrix" yomi="きほんぎょうれつ"} という。
:::

例えば $m=3$ のとき

$$
P_{12}=\begin{pmatrix}0&1&0\\1&0&0\\0&0&1\end{pmatrix},\quad
Q_2(c)=\begin{pmatrix}1&0&0\\0&c&0\\0&0&1\end{pmatrix},\quad
R_{13}(c)=\begin{pmatrix}1&0&c\\0&1&0\\0&0&1\end{pmatrix}
$$

である。

:::lemma{#lem-linalg-elementary-left-mult}
$A\in M_{m,n}(K)$ に行基本変形を施した結果は、対応する $m$ 次基本行列 $E$ を用いて $EA$ と表される。さらに基本行列は正則であり（後述 :ref[def-linalg-invertible]）、

$$
P_{ij}P_{ij}=I_m,\qquad Q_i(c)\,Q_i(c^{-1}) = I_m,\qquad R_{ij}(c)\,R_{ij}(-c) = I_m
$$

が成り立つ（左右を入れ替えても同様）。
:::

:::proof
$E=(e_{kl})$ を基本行列とすると、$(EA)$ の第 $k$ 行は $\sum_l e_{kl}\cdot(A\text{ の第 }l\text{ 行})$ である。$P_{ij}$ の第 $i,j$ 行はそれぞれ $\transp{e_j},\transp{e_i}$ だから、$P_{ij}A$ では第 $i,j$ 行が入れ替わる。$Q_i(c)$ の第 $i$ 行は $c\transp{e_i}$ だから、$Q_i(c)A$ では第 $i$ 行が $c$ 倍される。$R_{ij}(c)$ の第 $i$ 行は $\transp{e_i}+c\transp{e_j}$ だから、$R_{ij}(c)A$ では第 $i$ 行に第 $j$ 行の $c$ 倍が加わる。いずれの場合も、指定した行以外では $E$ の行が $I_m$ の対応する行に等しいため、$A$ の行がそのまま残る。

入れ替えは同じ入れ替えで戻る。$c$ 倍は $c^{-1}$ 倍で、$c$ 倍の加算は $-c$ 倍の加算で戻る。$I_m$ にこれらの操作と逆操作を続けて施すと、どちらの順でも $I_m$ が得られる。この事実を左からの積で表せば、主張の 3 つの等式と逆順の等式が得られる。
:::

:::example[行操作をまとめる]
$A=\begin{pmatrix}1&2\\3&5\end{pmatrix}$ に対し、第 2 行から第 1 行の 3 倍を引いてから、第 1 行と第 2 行を交換する。まとめた変形を表す行列は

$$
E=P_{12}R_{21}(-3)
=\begin{pmatrix}0&1\\1&0\end{pmatrix}
\begin{pmatrix}1&0\\-3&1\end{pmatrix}
=\begin{pmatrix}-3&1\\1&0\end{pmatrix}
$$

であり、$EA=\begin{pmatrix}0&-1\\1&2\end{pmatrix}$ となる。元に戻すには入れ替えを先に取り消すので、逆の積は $R_{21}(3)P_{12}$ である。
:::

### 列基本変形

行を混ぜる操作は方程式を取り替え、列を混ぜる操作は変数の取り方を変える。この違いを記録するため、列操作に対応する行列は右から掛ける。

:::proposition[列基本変形と右からの積]{#prop-linalg-elementary-right-mult}
$A\in M_{m,n}(K)$ の第 $i,j$ 列を入れ替える操作、第 $i$ 列を $c\ne0$ 倍する操作、第 $i$ 列に第 $j$ 列の $c$ 倍を加える操作は、それぞれ $n$ 次基本行列 $P_{ij}$、$Q_i(c)$、$R_{ji}(c)$ を右から掛けることに対応する。
:::

:::proof
転置した行列の対応する行に操作を施すと $E\transp{A}$ となる。これを転置し直せば :ref[prop-linalg-transpose-laws] により $A\transp{E}$ である。$P_{ij}$ と $Q_i(c)$ の転置はそれ自身である。$R_{ij}(c)=I_n+cE_{ij}$ の転置は $I_n+cE_{ji}=R_{ji}(c)$ なので、加算では添字が逆になる。
:::

右から基本行列の積 $F$ を掛けた係数行列を $AF$ とする。$(AF)y=b$ の解 $y$ に対して $x=Fy$ は $Ax=b$ の解である。逆に、各列操作が可逆なので $Ax=b$ の解から $y=F^{-1}x$ を得る。例えば $A=(1\ \ 2)$ の列を入れ替えると $(AF)y=y_1\cdot2+y_2$ となる。$b=1$ のとき $y=\transp{(0,1)}$ を選べば、元の解は $x=\transp{(1,0)}$ である。

### 行階段形と簡約行階段形

:::definition[行階段形・簡約行階段形]{#def-linalg-echelon}
行列の各行について、左から見て最初に現れる $0$ でない成分をその行の :term[主成分]{en="pivot" yomi="しゅせいぶん"}（ピボット）という。行列 $R$ が次の条件を満たすとき :term[行階段形]{en="row echelon form" yomi="ぎょうかいだんけい"} であるという。

1. 零行（成分がすべて $0$ の行）は、零でない行よりも下にある。
2. 零でない行の主成分の位置（列番号）は、下の行ほど真に右にある。

さらに次の条件も満たすとき、$R$ は :term[簡約行階段形]{en="reduced row echelon form" yomi="かんやくぎょうかいだんけい"} であるという。

3. 各主成分は $1$ である。
4. 主成分を含む列において、主成分以外の成分はすべて $0$ である。

主成分を含む列を**主列**という。
:::

例えば次の行列は簡約行階段形である（主成分を $\boxed{1}$ で示した）。

$$
\begin{pmatrix}
\boxed{1} & 2 & 0 & 2 & 3\\
0 & 0 & \boxed{1} & 1 & 2\\
0 & 0 & 0 & 0 & 0
\end{pmatrix}
$$

:::theorem[Gauss–Jordan の消去法]{#thm-linalg-gauss-jordan}
任意の行列は、有限回の行基本変形によって簡約行階段形に変形できる。
:::

:::proof
$A\in M_{m,n}(K)$ に対し、左から順に列を処理する。「処理済みの最初の $k$ 列からなる部分は簡約行階段形である」ことを $k=0,1,\dots,n$ について帰納的に示す。$k=0$ では処理済みの列がないため、この条件を満たす。

$k$ 列目まで処理が済んでおり、最初の $k$ 列の部分に主成分が $r$ 個あるとする。このとき第 $r+1$ 行以下の最初の $k$ 成分はすべて $0$ である。第 $k+1$ 列について、第 $r+1$ 行以下に $0$ でない成分がなければ、最初の $k+1$ 列の部分はすでに簡約行階段形である。そうでなければ、その成分を含む行を第 $r+1$ 行に移し、その成分の逆数倍をして $(r+1,k+1)$ 成分を $1$ にする。

続いて $i\ne r+1$ の各行について、第 $k+1$ 列の成分を $\alpha_i$ とし、第 $i$ 行から第 $r+1$ 行の $\alpha_i$ 倍を引く。これで第 $k+1$ 列は $e_{r+1}$ になる。この操作は、第 $r+1$ 行の最初の $k$ 成分がすべて $0$ であることから、最初の $k$ 列を変えない。従って最初の $k+1$ 列の部分が簡約行階段形になる。$k=n$ まで繰り返せばよい。
:::

この証明の手順は、そのまま計算のアルゴリズムになる。簡約行階段形は変形の手順によらず一意に定まる。

:::theorem[簡約行階段形の一意性]{#thm-linalg-rref-unique}
同じサイズの簡約行階段形行列 $R,R'$ について、$Rx=0$ と $R'x=0$ の解集合が同じならば $R=R'$ である。特に、各行列 $A$ に行同値な簡約行階段形はただ 1 つに定まる。これを $A$ の簡約行階段形と呼ぶ。
:::

証明には次の事実を用いる。

:::lemma{#lem-linalg-row-equiv-solutions}
$A$ と $B$ が行同値ならば、斉次方程式 $Ax=0$ と $Bx=0$ の解の集合は一致する。
:::

:::proof
:ref[lem-linalg-elementary-left-mult] より、ある基本行列の積 $E=E_k\cdots E_1$ によって $B=EA$ と書ける。$Ax=0$ ならば $Bx=EAx=0$ である。逆に $E$ は正則で $A=E^{-1}B$ だから、$Bx=0$ ならば $Ax=0$ である。
:::

:::proof[:ref[thm-linalg-rref-unique] の証明]
$R, R'\in M_{m,n}(K)$ を同じ斉次解をもつ簡約行階段形とし、$n$ に関する帰納法で $R=R'$ を示す。$n=1$ のとき、$1$ 列の簡約行階段形は $0$ か $e_1$ のどちらかである。前者の斉次解は $K$ 全体、後者の斉次解は $\{0\}$ なので、$R=R'$ である。

$n\ge 2$ とする。最後の列を除いた行列 $R_0, R_0'$ も簡約行階段形である。最後の列が主列の場合はその行が零行になるが、その下には元から零行しかないためである。$R_0u=0$ と $R\begin{pmatrix}u\\0\end{pmatrix}=0$ は同値であり、$R'_0$ と $R'$ にも同じことがいえる。従って $R_0$ と $R_0'$ は同じ斉次解をもち、帰納法の仮定より $R_0=R_0'$ となる。

$R\ne R'$ と仮定して矛盾を導く。$Rx=0$ を満たす任意の $x$ は、仮定により $R'x=0$ も満たす。$R-R'$ は最後の列 $d\ne 0$ 以外が零なので $(R-R')x=x_nd=0$ であり、すべての解は $x_n=0$ を満たす。

もし $R$ の第 $n$ 列が主列でなければ、$x_n=1$ とし、他の非主列に対応する変数を $0$ として、主列に対応する変数を決めれば $Rx=0$ の解が得られる（主列は単位ベクトルなので、各主成分の行の方程式がちょうど 1 つの主列変数を決める）。これは $x_n=0$ に反する。よって第 $n$ 列は $R$ の主列であり、同様に $R'$ の主列でもある。$R_0=R_0'$ の零でない行の個数を $r$ とすれば、$R$ と $R'$ の最後の列はともに $e_{r+1}$ となり、$R=R'$ が従う。これは仮定に反する。

最後に、$A$ と行同値な二つの簡約行階段形は、:ref[lem-linalg-row-equiv-solutions] によりいずれも $A$ と同じ斉次解をもつため等しい。存在は :ref[thm-linalg-gauss-jordan] による。
:::

### 前進消去と後退代入

方程式を一度解く場合は、主成分の下だけを消去して行階段形にする方法もある。これを :term[前進消去]{en="forward elimination" yomi="ぜんしんしょうきょ"} と呼ぶ。:ref[thm-linalg-gauss-jordan] の手順から主成分の正規化と上側の消去を省き、主成分を $p\ne0$ として下の各行から主成分の行の $a_{ij}/p$ 倍を引く。採用した行の左側は処理済みの列ですべて $0$ だから、その部分の階段形は保たれる。

行階段形になった係数行列で、主列に対応する未知数を**主変数**、それ以外を :term[自由変数]{en="free variable" yomi="じゆうへんすう"} という。最下段の零でない方程式から順に主変数を求める操作を :term[後退代入]{en="back substitution" yomi="こうたいだいにゅう"} という。第 $i$ 行の主列を $j_i$ とし、右辺を $c_i$ とすると、

$$
x_{j_i}=\frac{c_i-\sum_{k>j_i}a_{ik}x_k}{a_{i,j_i}}.
$$

$j_i$ より右にある主変数は下段で決まる。したがって自由変数を先に決めておけば、この式を下から上へ適用できる。矛盾する零行がなければ、逆にどの解もこの順で値を復元できる。

行列そのものを簡約行階段形にしたい場合は、前進消去の後、主成分を下から順に $1$ に正規化して上側を消す。主成分を含む行はそれより左の成分が $0$ なので、左にある主成分は変わらない。さらに右の主列はすでに簡約されており、いま使う行のそれらの成分は $0$ なので、処理済みの主列も変わらない。

:::example[零を避ける行交換と後退代入]
実数の範囲で、拡大係数行列を次のように変形する。

$$
\begin{aligned}
\left(\begin{array}{ccc|c}0&2&1&5\\1&1&0&3\\2&3&1&9\end{array}\right)
&\longrightarrow
\left(\begin{array}{ccc|c}1&1&0&3\\0&2&1&5\\0&1&1&3\end{array}\right)\\
&\longrightarrow
\left(\begin{array}{ccc|c}1&1&0&3\\0&2&1&5\\0&0&\tfrac12&\tfrac12\end{array}\right).
\end{aligned}
$$

最初の矢印では第 1 行と第 2 行を交換し、新しい第 1 行の 2 倍を第 3 行から引いた。次に第 2 行の $1/2$ 倍を第 3 行から引いた。最下段から $x_3=1$、第 2 行から $x_2=(5-1)/2=2$、第 1 行から $x_1=3-2=1$ である。右辺の列まで一緒に更新するため、この解は元の連立方程式も満たす。
:::

## 連立一次方程式の解法

### 行基本変形と解の集合

前進消去や簡約化で変えるのは方程式の書き方である。そこで、変形後の方程式が元の方程式と同じ解をもつことを示しておく。

:::proposition{#prop-linalg-row-ops-preserve-solutions}
拡大係数行列 $(A\mid b)$ に行基本変形を施して $(A'\mid b')$ が得られたとすると、$Ax=b$ と $A'x=b'$ の解の集合は一致する。
:::

:::proof
ある正則行列 $E$（基本行列の積）により $A'=EA$, $b'=Eb$ と書ける。$Ax=b$ ならば $A'x=EAx=Eb=b'$ であり、逆も $E^{-1}$ を掛ければよい。
:::

したがって、$(A\mid b)$ を簡約行階段形 $(R\mid c)$ に変形してから解けばよい。主変数と自由変数の区別には、係数側 $R$ の主列を用いる。

:::theorem[解の存在と表示]{#thm-linalg-solve-rref}
$A\in M_{m,n}(K)$、$b\in K^m$ とする。$(A\mid b)$ の簡約行階段形を $(R\mid c)$ とし、$R$ の零でない行の個数を $r$ とする。

1. $Ax=b$ が解をもつための必要十分条件は、$(R\mid c)$ の最後の列が主列でないこと、すなわち $c_{r+1}=\cdots=c_m=0$ である。
2. 解をもつとき、自由変数（$n-r$ 個）に任意の値を与えると、主変数は定数項を含みうる一次式として一意に定まる。こうしてすべての解が得られる。
:::

:::proof
先に、$R$ 自身が簡約行階段形であることを確認する。最後の列を削除しても、$R$ 内にある主成分の位置や値は変わらない。削除した列に主成分があった場合、その行は最後の非零行であり、削除後は零行になる。従っていずれの場合も簡約行階段形の条件が保たれる。

1. 最後の列が主列であることは、$R$ が零行である行に $c$ の非零成分があることと同値である。その行は $(0,\dots,0\mid 1)$ であり、方程式 $0=1$ を与えるので解はない。最後の列が主列でなければ、第 $r+1$ 行以下は $0=0$ である。また、第 $i$ 行（$i\le r$）の方程式は、主成分の列に対応する主変数 $x_{j_i}$ について

$$
x_{j_i} = c_i - \sum_{k:\,\text{自由変数}} r_{ik}x_k
$$

と解かれた形をしている。例えばすべての自由変数を $0$ とすれば、残る式を満たす解が得られる。

2. 自由変数を決めて上式で主変数を定めると、第 $1$ 行から第 $r$ 行までの式をすべて満たす。残りの行は $0=0$ なので追加の条件はない。逆に任意の解の主変数はこの式を満たすため、その解の自由変数を選べば元の解が再現される。
:::

:::example[自由変数をもつ連立方程式]{#exm-linalg-solve}
次の連立方程式を解く。

$$
\begin{cases}
x_1 + 2x_2 - x_3 + x_4 = 1\\
2x_1 + 4x_2 - x_3 + 3x_4 = 4\\
-x_1 - 2x_2 + 2x_3 = 1
\end{cases}
$$

拡大係数行列に行基本変形を施す。第 2 行から第 1 行の 2 倍を引き、第 3 行に第 1 行を加え、さらに第 3 行から第 2 行を引き、第 1 行に第 2 行を加えると

$$
\left(\begin{array}{cccc|c}
1&2&-1&1&1\\ 2&4&-1&3&4\\ -1&-2&2&0&1
\end{array}\right)
\longrightarrow
\left(\begin{array}{cccc|c}
1&2&0&2&3\\ 0&0&1&1&2\\ 0&0&0&0&0
\end{array}\right)
$$

となる。最後の列は主列でないので解が存在する。主変数は $x_1, x_3$、自由変数は $x_2, x_4$ である。$x_2=s$, $x_4=t$（$s,t\in K$ は任意）とおけば

$$
x = \begin{pmatrix}3-2s-2t\\ s\\ 2-t\\ t\end{pmatrix}
= \begin{pmatrix}3\\0\\2\\0\end{pmatrix} + s\begin{pmatrix}-2\\1\\0\\0\end{pmatrix} + t\begin{pmatrix}-2\\0\\-1\\1\end{pmatrix}.
$$
:::

:::example[係数の消去と右辺の条件]
実数 $a,b$ に対し

$$
\begin{cases}
x+y=1,\\
x+ay=b
\end{cases}
\quad\longleftrightarrow\quad
\left(\begin{array}{cc|c}1&1&1\\1&a&b\end{array}\right)
\longrightarrow
\left(\begin{array}{cc|c}1&1&1\\0&a-1&b-1\end{array}\right)
$$

と変形する。$a\ne1$ のときは $y=(b-1)/(a-1)$、$x=1-y=(a-b)/(a-1)$ に一意に定まる。$a=1$ の場合は、最後の式は $0=b-1$ になる。したがって $b=1$ ならば $(x,y)=(1-t,t)$ が全解になり、$b\ne1$ ならば解は存在しない。

分母 $a-1$ が $0$ になる場合を先に分けると、除算で解を失わずにすむ。分岐点では、係数が消えて右辺だけが残るかどうかを調べる。
:::

:::example[複素係数の消去]
$x+iy=2$、$ix+y=1+i$ を $\C$ 上で解く。第 2 式から第 1 式の $i$ 倍を引くと $(1-i^2)y=1-i$ なので、$y=(1-i)/2$ である。第 1 式に戻すと $x=2-i(1-i)/2=(3-i)/2$ となる。実際、

$$
x+iy=\frac{3-i+i+1}{2}=2,\qquad
ix+y=\frac{3i+1+1-i}{2}=1+i.
$$

消去で用いるのは $i$ 倍や $2$ による除算など、実係数の場合と同じ種類の演算である。
:::

### 同じ係数行列で複数の右辺を解く

$B=(b_1\ \cdots\ b_p)$ に対し、$AX=B$ は $p$ 本の連立方程式 $Ax_j=b_j$ をまとめたものである。係数行列 $A$ を簡約化する行操作を $(A\mid B)$ の全列に一度施せば、各右辺を別々に処理するときと同じ結果が得られる。実際、基本行列の積 $E$ に対して $E(A\mid B)=(EA\mid EB)$ であり、$EB$ の第 $j$ 列は $Eb_j$ である。ここで必要なのは左側 $EA$ の簡約である。右辺の列を使って追加の消去をする前に、各列の矛盾を独立に判定できる。

:::example[同じ左辺でも解の存在は右辺で変わる]
実数の範囲で

$$
\left(\begin{array}{cc|cc}
1&2&3&3\\
2&4&6&7
\end{array}\right)
\longrightarrow
\left(\begin{array}{cc|cc}
1&2&3&3\\
0&0&0&1
\end{array}\right)
$$

とする。操作は第 2 行から第 1 行の 2 倍を引くことだけである。第 1 の右辺では $(x_1,x_2)=(3-2t,t)$ が全解である。第 2 の右辺では $0=1$ が出るので解はない。$AX=B$ として二つを同時に要求すると、第 2 列を作れないため $X$ は存在しない。
:::

## 階数と解の構造

### 階数

係数の消去で零行が現れると、その行の右辺を見て系が矛盾するかどうかを判定できる。一方、残る係数の主成分を数えると、可解な場合に何個の主変数が決まるかが分かる。

:::definition[階数]{#def-linalg-rank}
行列 $A$ の簡約行階段形の零でない行の個数（すなわち主成分の個数）を $A$ の :term[階数]{en="rank" yomi="かいすう"} といい、$\rank A$ で表す。
:::

:ref[thm-linalg-rref-unique] によって簡約行階段形は一意なので、階数は消去の手順によらず定まる。行同値な行列は同じ簡約行階段形をもつので、階数も等しい。行階段形から主成分を下から順に正規化して上側を消す操作では主成分の位置が変わらない。従って、前進消去を終えた時点でも零でない行を数えれば階数が得られる。$A\in M_{m,n}(K)$ の各主成分は異なる行と列に位置するので、$\rank A\le\min\{m,n\}$ である。階数と行・列の張る空間との関係は [線形写像と表現行列](/math/linear-algebra/linear-maps/) で扱う。

:::theorem[解の存在と自由度]{#thm-linalg-rank-criterion}
$A\in M_{m,n}(K)$, $b\in K^m$ とする。

1. $Ax=b$ が解をもつための必要十分条件は $\rank A = \rank (A\mid b)$ である。
2. 解をもつとき、一般解は $n-\rank A$ 個の任意定数を含む。従って $Ax=b$ がちょうど 1 つの解をもつための必要十分条件は $\rank (A\mid b)=\rank A=n$ である。
:::

:::proof
:ref[thm-linalg-solve-rref] の証明により、$(A\mid b)$ の簡約行階段形 $(R\mid c)$ の係数部分 $R$ は簡約行階段形である。同じ行操作を $A$ に施せば $R$ が得られるため、$R$ は $A$ の簡約行階段形である。よって $\rank(A\mid b)$ が $\rank A$ に等しいことは、最後の列が主列でないことと同値である。これと :ref[thm-linalg-solve-rref] により 1 が従う。

解があるときは $n$ 個の変数から $\rank A$ 個の主変数を除いたものが自由変数である。各自由変数の値は解の成分にそのまま現れるので、自由変数があれば $0$ と $1$ を選ぶことで異なる解が得られる。自由変数がなければ主変数の値がすべて決まり、解はただ 1 つになる。これで 2 が従う。
:::

例えば $A=\transp{(1,0)}$ と $b=\transp{(0,1)}$ では、未知数は $1$ 個で $\rank A=1$ である。しかし第 2 式が $0=1$ であるため、$\rank(A\mid b)=2$ となり解はない。このように $\rank A=n$ から分かるのは、解が高々一つであることまでである。

:::corollary{#cor-linalg-homogeneous-nontrivial}
斉次方程式 $Ax=0$（$A\in M_{m,n}(K)$）が自明でない解 $x\ne 0$ をもつ $\iff$ $\rank A<n$。特に $m<n$（方程式の本数が未知数の個数より少ない）ならば、自明でない解が存在する。
:::

:::proof
斉次方程式は常に $x=0$ を解にもつ。:ref[thm-linalg-rank-criterion] より、解が $x=0$ のみであることは $\rank A=n$ と同値である。$m<n$ ならば $\rank A\le m<n$ である。
:::

### 解の構造

:::proposition[一般解 = 特殊解 + 斉次方程式の一般解]{#prop-linalg-solution-structure}
$x_0$ を $Ax=b$ の 1 つの解（**特殊解**）とする。このとき $Ax=b$ の解全体は

$$
\{\, x_0 + y \mid Ay = 0 \,\}
$$

に等しい。
:::

:::proof
$Ay=0$ ならば $A(x_0+y)=Ax_0+Ay=b$ である。逆に $Ax=b$ ならば $y\defeq x-x_0$ は $Ay=b-b=0$ を満たし、$x=x_0+y$ である。
:::

:ref[exm-linalg-solve] では $x_0=\transp{(3,0,2,0)}$ が特殊解であり、$s, t$ を含む部分が斉次方程式の一般解である。

この表示を一般に作るには、$r=\rank A$ とし、自由変数をすべて $0$ にした特殊解を $x_0$ とする。次に斉次方程式の解で、第 $j$ 自由変数を $1$、残りをすべて $0$ としたものを $v_j$（$j=1,\dots,n-r$）とする。和 $\sum_jt_jv_j$ は斉次方程式を満たし、その自由変数は $t_1,\dots,t_{n-r}$ である。従って、これが指定した自由変数をもつ唯一の斉次解であり、非斉次系の全解は $x_0+\sum_jt_jv_j$ と書ける。自由変数がなければ空の和を $0$ とする。この $v_j$ の組は [ベクトル空間・基底・次元](/math/linear-algebra/vector-spaces/) で解空間の基底として扱う。

:::example[矛盾と自由度を分けて読む]
$A=\begin{pmatrix}1&0&1\\0&1&-1\\1&1&0\end{pmatrix}$、$b=\transp{(2,3,h)}$ とする。第 3 行から第 1 行と第 2 行を引くと、係数側の第 3 行は零行になり、右辺は $h-5$ になる。よって $\rank A=2$ であり、$h=5$ ならば $\rank(A\mid b)=2$、$h\ne5$ ならば $\rank(A\mid b)=3$ である。

$h=5$ のとき、$x_3=0$ とした特殊解は $\transp{(2,3,0)}$、斉次系で $x_3=1$ とした解は $\transp{(-1,1,1)}$ である。したがって

$$
x=\begin{pmatrix}2\\3\\0\end{pmatrix}
+t\begin{pmatrix}-1\\1\\1\end{pmatrix},\qquad t\in K.
$$

もし $h\ne5$ であれば、自由変数をどう選んでも $0=h-5$ は満たせない。変数の自由度を数える前に、まず右辺の条件を判定する。
:::

### 解が存在しないことの証明

右辺に対する制約は、消去で零行を作った行の組み合わせから求められる。行操作を単位行列にも施し、その組み合わせを記録しておくと、元の方程式からどのように矛盾が出たかを示せる。

:::proposition[可解条件を与える行の組み合わせ]{#prop-linalg-left-annihilator-criterion}
$A\in M_{m,n}(K)$、$b\in K^m$ とする。$Ax=b$ に解があることと、すべての $y\in K^m$ について $\transp{y}A=0\Rightarrow\transp{y}b=0$ が成り立つことは同値である。
:::

:::proof
$Ax=b$ に解があれば、$\transp{y}b=(\transp{y}A)x=0$ である。逆を示すため、$E$ を行操作の基本行列の積として $EA=R$ を簡約行階段形にする。$r=\rank A$ とおくと、$E$ の第 $r+1$ 行以下をそれぞれ $\transp{y}$ としたとき $\transp{y}A=0$ である。仮定から、それらの行では $\transp{y}b=0$ でもある。すなわち $Eb$ の第 $r+1$ 成分以降は $0$ なので、主変数を決める式に戻って $Rx=Eb$ を解ける。$E$ は可逆なので $Ax=b$ にも解がある。
:::

:::example[右辺に課される二つの条件]
$A=\begin{pmatrix}1&2\\2&4\\-1&-2\end{pmatrix}$ とする。第 2 行から第 1 行の 2 倍を引き、第 3 行に第 1 行を加えると後ろの 2 行が零になる。したがって $b=\transp{(b_1,b_2,b_3)}$ に対する可解条件は

$$
b_2-2b_1=0,\qquad b_3+b_1=0
$$

である。この場合の $E$ の後ろの 2 行は $(-2,1,0)$、$(1,0,1)$ であり、それぞれを $A$ の左から掛けると零行になる。例えば $b=\transp{(1,2,0)}$ では $(1,0,1)A=0$ だが $(1,0,1)b=1$ なので、解がないことをこの 1 行で示せる。
:::

### 主列から得られる階数分解

簡約行階段形の主列は単位ベクトルなので、残りの列がそれらのどの組み合わせでできているかを読み取れる。行操作は列どうしの関係を保つので、その係数を元の行列にも使える。

:::theorem[主列による階数分解]{#thm-linalg-rank-factorization}
$A\in M_{m,n}(K)$ の簡約行階段形を $R$ とし、その主列番号を $j_1<\cdots<j_r$ とする。$r=\rank A>0$ のとき、$C=(a_{j_1}\ \cdots\ a_{j_r})$ を $A$ の対応する列で作り、$D$ を $R$ の最初の $r$ 行とする。このとき

$$
A=CD,\qquad \rank C=\rank D=r.
$$

$A=FG$ と書ける任意の分解について、$F$ の列数（$G$ の行数）は $r$ 以上である。
:::

:::proof
$R=EA$ と書くと $EC=\begin{pmatrix}I_r\\O\end{pmatrix}$ である。従って $ECD=\begin{pmatrix}D\\O\end{pmatrix}=R=EA$ であり、$E$ の逆行列を左から掛けると $CD=A$ となる。$EC$ の主列は $r$ 本あり、$D$ は簡約行階段形で零でない行が $r$ 本あるので、どちらの階数も $r$ である。

$A=FG$ と書け、$F$ の列数が $k<r$ だと仮定する。$G$ の第 $j_1,\dots,j_r$ 列を集めた $k\times r$ 行列を $G_0$ とすると $C=FG_0$ である。:ref[cor-linalg-homogeneous-nontrivial] により $G_0u=0$ を満たす $u\ne0$ があるので $Cu=0$ になる。しかし $ECu=\transp{(u_1,\dots,u_r,0,\dots,0)}$ なので $Cu=0$ は $u=0$ を強制する。矛盾から $k\ge r$ である。
:::

零行列は階数 $0$ である。空の行・列を許す記法を採用すれば、階数分解は内側のサイズが $0$ の積としても書ける。その積の全成分は空の和なので $0$ である。

:::example[元の行列の列を選ぶ]{#exm-linalg-rank-factorization}
$A=\begin{pmatrix}1&2&3\\2&4&6\\1&3&5\end{pmatrix}$ とする。第 2 行から第 1 行の 2 倍を引き、第 3 行から第 1 行を引く。零行を末尾に移し、第 1 行から新しい第 2 行の 2 倍を引くと

$$
R=\begin{pmatrix}1&0&-1\\0&1&2\\0&0&0\end{pmatrix},\qquad
A=\begin{pmatrix}1&2\\2&4\\1&3\end{pmatrix}
\begin{pmatrix}1&0&-1\\0&1&2\end{pmatrix}
$$

を得る。右側の因子の第 3 列から、元の行列でも $a_3=-a_1+2a_2$ と分かる。一方、簡約行階段形の主列 $e_1,e_2$ を左側の因子に使えば積は $R$ になる。元の列を復元するには、主列番号だけを $R$ から読み取り、列そのものは $A$ から取る。
:::

:::corollary[積・転置・和の階数]{#cor-linalg-rank-calculus}
行列のサイズがそれぞれの演算に適合するとき、

$$
\rank(AB)\le\min\{\rank A,\rank B\},\quad
\rank(\transp{A})=\rank A,\quad
\rank(A+B)\le\rank A+\rank B
$$

が成り立つ。さらに、正則行列を左右から掛けても階数は変わらない。
:::

:::proof
$A=CD$ が階数分解ならば $AB=C(DB)$ なので、分解の内側のサイズの下限より $\rank(AB)\le\rank A$ である。$B$ の階数分解を使って $AB=(AC')D'$ と書けば $\rank(AB)\le\rank B$ を得る。

$\transp{A}=\transp{D}\transp{C}$ から $\rank(\transp{A})\le\rank A$ である。$A$ と $\transp{A}$ を取り替えると逆向きの不等式が得られる。和は $A+B=(C\mid C')\begin{pmatrix}D\\D'\end{pmatrix}$ と書けるので、同じ下限を使えばよい。各議論で零行列が現れる場合は、零との積が零であること、零の転置が零であること、零を足しても変わらないことから該当する結論が従う。

正則な $P,Q$ に対しては積の不等式から $\rank(PAQ)\le\rank A$ である。逆に $A=P^{-1}(PAQ)Q^{-1}$ なので $\rank A\le\rank(PAQ)$ である。
:::

### 行と列の両方を用いる標準形

:::proposition[階数による標準形]{#prop-linalg-elementary-rank-form}
$A\in M_{m,n}(K)$ の階数を $r$ とする。行基本変形と列基本変形を組み合わせると、左上に $I_r$ があり、他の成分がすべて $0$ である $m\times n$ 行列に変形できる。従って、同じサイズの二つの行列が行・列基本変形で互いに移り合うことと、階数が等しいことは同値である。
:::

:::proof
行基本変形で簡約行階段形にし、主列を左側に元の順で集める。主列は上から $r$ 個の基本ベクトルなので、変形後は $\begin{pmatrix}I_r&F\\O&O\end{pmatrix}$ である。右側の各列から、左側の第 $i$ 列の $F_{ij}$ 倍を $i=1,\dots,r$ について順に引くと、その列の成分がすべて消える。$r=0$ なら行列は最初から零行列である。

階数が等しければ、両方を同じ標準形に変形し、片方の操作を逆順に戻して互いに移せる。逆に、行・列基本変形は左・右から可逆な基本行列を掛ける操作なので、:ref[cor-linalg-rank-calculus] により階数を保つ。
:::

例えば $R=\begin{pmatrix}1&0&-1\\0&1&2\\0&0&0\end{pmatrix}$ では、第 3 列に第 1 列を足し、第 2 列の 2 倍を引けば $\diag(1,1,0)$ になる。行だけで簡約化した時点では第 3 列の関係が残り、列操作を加えると階数だけが残る。

:::corollary[Sylvester の階数不等式]{#cor-linalg-rank-sylvester}
$A\in M_{m,n}(K)$、$B\in M_{n,p}(K)$ ならば

$$
\rank(AB)\ge \rank A+\rank B-n
$$

である。
:::

:::proof
$r=\rank A$ とする。$r=0$ または $\rank B=0$ なら右辺は $0$ 以下なので成り立つ。$r=n$ なら、階数標準形を用いて $PAQ=\begin{pmatrix}I_n\\O\end{pmatrix}$ と書ける。よって $PAB=\begin{pmatrix}Q^{-1}B\\O\end{pmatrix}$ であり、零行の追加は階数を変えないので $\rank(AB)=\rank B$ である。

$0<r<n$ のときも :ref[prop-linalg-elementary-rank-form] により、可逆な $P,Q$ で $PAQ$ を階数標準形にできる。$Q^{-1}B$ を上 $r$ 行の $B_1$ と下 $n-r$ 行の $B_2$ に分ける。このとき

$$
PAB=(PAQ)(Q^{-1}B)=\begin{pmatrix}B_1\\O\end{pmatrix},
\qquad
Q^{-1}B=\begin{pmatrix}B_1\\O\end{pmatrix}+\begin{pmatrix}O\\B_2\end{pmatrix}.
$$

ここで、零行を移動して末尾に置いてから消去すれば、その零行は主成分を増やさない。従って $\rank\begin{pmatrix}O\\B_2\end{pmatrix}=\rank B_2\le n-r$ である。和の階数の上界と、可逆行列を掛けても階数が変わらないこと（:ref[cor-linalg-rank-calculus]）を使えば、$\rank B\le\rank(AB)+(n-r)$ を得る。
:::

例えば $A=\diag(1,1,0)$、$B=\diag(0,1,1)$ の積は $\diag(0,1,0)$ である。両方の階数が $2$ で、積の階数は $1$ なので下界 $2+2-3=1$ に達する。

## 正則行列と逆行列

$Ax=b$ をどの右辺についても一意に解けるなら、$b$ から $x$ を復元する計算も一つの行列で表せる。その行列は $A$ を掛ける操作を元に戻す。

:::definition[正則行列・逆行列]{#def-linalg-invertible}
$A\in M_n(K)$ に対し、$AB=BA=I_n$ を満たす $B\in M_n(K)$ が存在するとき、$A$ は :term[正則]{en="invertible" yomi="せいそく"} であるという。このような $B$ はただ 1 つであり（下の命題）、これを $A$ の :term[逆行列]{en="inverse matrix" yomi="ぎゃくぎょうれつ"} といい $A^{-1}$ と書く。
:::

:::proposition{#prop-linalg-inverse-basic}
1. 逆行列は存在すれば一意である。
2. $A, B\in M_n(K)$ が正則ならば $AB$ も正則で、$(AB)^{-1}=B^{-1}A^{-1}$。
3. $A$ が正則ならば $\transp{A}$ も正則で、$(\transp{A})^{-1}=(A^{-1})^\top$。
:::

:::proof
1. $B, B'$ がともに逆行列ならば $B = BI_n = B(AB') = (BA)B' = I_nB' = B'$。
2. $(AB)(B^{-1}A^{-1}) = A(BB^{-1})A^{-1} = AA^{-1} = I_n$、逆順も同様。
3. :ref[prop-linalg-transpose-laws] より $\transp{A}(A^{-1})^\top = (A^{-1}A)^\top = I_n$、逆順も同様。
:::

次の定理は、正方行列の正則性を確かめるには片側の逆行列だけで十分であることを含む、基本的な結果である。

:::theorem[正則性の同値条件]{#thm-linalg-invertible-equiv}
$A\in M_n(K)$ について、次の条件は互いに同値である。

1. $A$ は正則である。
2. $\rank A = n$。
3. $A$ の簡約行階段形は $I_n$ である。
4. $Ax=0$ の解は $x=0$ のみである。
5. $A$ は基本行列の積で表される。
6. $BA=I_n$ を満たす $B\in M_n(K)$ が存在する。
7. $AB=I_n$ を満たす $B\in M_n(K)$ が存在する。
8. すべての $b\in K^n$ に対して $Ax=b$ は解をもつ。
:::

:::proof
(1)$\Rightarrow$(6), (1)$\Rightarrow$(7) は定義から明らかである。

(6)$\Rightarrow$(4)：$Ax=0$ ならば $x = BAx = 0$。

(4)$\Rightarrow$(2)：:ref[cor-linalg-homogeneous-nontrivial] による。

(2)$\Rightarrow$(3)：$n$ 次正方行列の簡約行階段形が $n$ 個の主成分をもてば、すべての列が主列であり、簡約行階段形の定義から $I_n$ に等しい。

(3)$\Rightarrow$(5)：基本行列 $E_1,\dots,E_k$ により $E_k\cdots E_1A=I_n$ と書ける。各 $E_i$ は正則なので $A = E_1^{-1}\cdots E_k^{-1}$ であり、:ref[lem-linalg-elementary-left-mult] より基本行列の逆行列はまた基本行列である。

(5)$\Rightarrow$(1)：基本行列は正則であり、正則行列の積は正則である（:ref[prop-linalg-inverse-basic]）。

(7)$\Rightarrow$(1)：$AB=I_n$ ならば $\transp{B}\,\transp{A}=I_n$ なので、$\transp{A}$ に対して (6) が成り立つ。すでに示した (6)$\Rightarrow$(1) により $\transp{A}$ は正則であり、:ref[prop-linalg-inverse-basic] の 3 から $A=(\transp{A})^\top$ も正則である。

(1)$\Rightarrow$(8)：$x=A^{-1}b$ は $Ax=b$ を満たす。(8)$\Rightarrow$(7)：右辺 $e_j$ に対する解を $u_j$ とする。$B=(u_1\ \cdots\ u_n)$ とおけば、$AB$ の第 $j$ 列は $Au_j=e_j$ である。
:::

:::remark
(6) や (7) の $B$ は自動的に $A^{-1}$ に等しい。(6) なら $B=(BA)A^{-1}=A^{-1}$、(7) なら $B=A^{-1}(AB)=A^{-1}$ である。長方形の行列では、左と右で単位行列のサイズが異なる。以下でそれぞれの存在条件を調べる。
:::

:::proposition[正則な部分行列と階数]{#prop-linalg-rank-invertible-submatrix}
行列 $A$ の階数が $r>0$ ならば、行と列をそれぞれ $r$ 本ずつ選んで得られる正則な $r$ 次正方行列がある。さらに、正則な正方部分行列のサイズは $r$ 以下である。階数 $0$ の行列には正則な正方部分行列がない。
:::

:::proof
行や列の選択は、単位行列から選んだ行や列を並べた行列を左右から掛ける操作である。例えば第 $i$ 行の選択には左から $\transp{e_i}$ を、第 $j$ 列の選択には右から $e_j$ を掛ける。従って正方部分行列の階数は、:ref[cor-linalg-rank-calculus] により $\rank A$ 以下である。正則ならそのサイズと階数は等しいので、サイズの上界が従う。

一方 $r>0$ ならば、:ref[thm-linalg-rank-factorization] の列の選び方で階数 $r$ の $m\times r$ 行列 $C$ を取れる。$\transp{C}$ の階数も $r$ なので、同じ選び方で $\transp{C}$ から $r$ 本の列を選べば、階数 $r$ の $r$ 次正方行列が得られる。これを転置した行列は $C$ から $r$ 本の行を選んだものであり、:ref[thm-linalg-invertible-equiv] により正則である。
:::

例えば :ref[exm-linalg-rank-factorization] の行列の階数は $2$ である。第 1、3 行と第 1、2 列を選ぶと $\begin{pmatrix}1&2\\1&3\end{pmatrix}$ が得られ、第 2 行から第 1 行を引くと主成分が二つ現れる。

### 長方形の行列の片側逆

:::proposition[左逆行列と右逆行列]{#prop-linalg-one-sided-inverses}
$A\in M_{m,n}(K)$ とする。$LA=I_n$ を満たす $L\in M_{n,m}(K)$ を**左逆行列**、$AR=I_m$ を満たす $R\in M_{n,m}(K)$ を**右逆行列**という。左逆行列が存在することは $\rank A=n$ と、右逆行列が存在することは $\rank A=m$ と同値である。
:::

:::proof
$LA=I_n$ なら $\rank A\ge\rank(LA)=n$ であり、$\rank A\le n$ と合わせて $\rank A=n$ である。逆に $\rank A=n$ なら、行操作を表す $E$ により $EA=\begin{pmatrix}I_n\\O\end{pmatrix}$ とできる。$E$ の上から $n$ 行を $L$ とすれば $LA=I_n$ になる。

右逆行列については $AR=I_m$ を転置して $\transp{R}\transp{A}=I_m$ とする。これは $\transp{A}$ の左逆行列の条件である。上の結果と $\rank(\transp{A})=\rank A$ から結論が従う。
:::

:::example[右逆行列が複数ある場合]
$A=\begin{pmatrix}1&0&1\\0&1&1\end{pmatrix}$ は階数 $2$ なので右逆行列をもつ。各右辺 $e_1,e_2$ に対して第 3 成分を自由に選ぶと、右逆行列全体は

$$
R=\begin{pmatrix}1-s&-t\\-s&1-t\\s&t\end{pmatrix},
\qquad s,t\in K
$$

となる。$AR$ の第 1 行は $R$ の第 1 行と第 3 行の和、第 2 行は第 2 行と第 3 行の和なので、$AR=I_2$ である。一方 $A$ の列数は $3$ なので、左逆行列の条件 $\rank A=3$ は満たせない。
:::

### 逆行列の計算

:ref[thm-linalg-invertible-equiv] の証明から、次の計算法が得られる。$n\times 2n$ 行列 $(A\mid I_n)$ に行基本変形を施して左半分を簡約化する。左半分が $I_n$ になれば $A$ は正則で、そのとき右半分が $A^{-1}$ である。左半分が $I_n$ にならなければ $A$ は正則でない。

実際、行基本変形の列に対応する基本行列の積を $E$ とすれば $E(A\mid I_n)=(EA\mid E)$ であり、$EA=I_n$ ならば $E=A^{-1}$ である（:ref[thm-linalg-invertible-equiv] の (6)）。

:::example{#exm-linalg-inverse}
$A=\begin{pmatrix}1&2&0\\0&1&1\\1&0&1\end{pmatrix}$ の逆行列を求める。以下では $K=\R,\C$ のように $3\ne0$ である場合を扱う。

$$
\left(\begin{array}{ccc|ccc}
1&2&0&1&0&0\\ 0&1&1&0&1&0\\ 1&0&1&0&0&1
\end{array}\right)
\to
\left(\begin{array}{ccc|ccc}
1&2&0&1&0&0\\ 0&1&1&0&1&0\\ 0&-2&1&-1&0&1
\end{array}\right)
\to
\left(\begin{array}{ccc|ccc}
1&0&-2&1&-2&0\\ 0&1&1&0&1&0\\ 0&0&3&-1&2&1
\end{array}\right)
$$

第 3 行を $1/3$ 倍し、第 1 行に第 3 行の 2 倍を加え、第 2 行から第 3 行を引くと

$$
\left(\begin{array}{ccc|ccc}
1&0&0&\tfrac13&-\tfrac23&\tfrac23\\[2pt] 0&1&0&\tfrac13&\tfrac13&-\tfrac13\\[2pt] 0&0&1&-\tfrac13&\tfrac23&\tfrac13
\end{array}\right),
\qquad
A^{-1} = \frac13\begin{pmatrix}1&-2&2\\1&1&-1\\-1&2&1\end{pmatrix}.
$$

例えば $AA^{-1}$ の第 1 行は $\frac13(1,-2,2)+\frac23(1,1,-1)=(1,0,0)$ である。第 2 行は $\frac13(1,1,-1)+\frac13(-1,2,1)=(0,1,0)$、第 3 行は $\frac13(1,-2,2)+\frac13(-1,2,1)=(0,0,1)$ である。$3=0$ の体では、途中の左側の第 3 行が零行になり、階数が $2$ なので逆行列は存在しない。
:::

:::example[二次正方行列の公式]
$A=\begin{pmatrix}a&b\\c&d\end{pmatrix}$ とする。$\Delta=ad-bc$ とおくと、

$$
\begin{pmatrix}a&b\\c&d\end{pmatrix}\begin{pmatrix}d&-b\\-c&a\end{pmatrix}
=\begin{pmatrix}\Delta&0\\0&\Delta\end{pmatrix}.
$$

従って $\Delta\ne0$ なら、右側の因子を $\Delta^{-1}$ 倍したものが $A^{-1}$ である。$\Delta=0$ で第 1 行が零でなければ、$x=\transp{(-b,a)}\ne0$ は $Ax=\transp{(0,\Delta)}=0$ を満たす。第 1 行が零なら $\rank A\le1$ である。いずれの場合も :ref[thm-linalg-invertible-equiv] により逆行列は存在しない。例えば $\begin{pmatrix}2&1\\3&2\end{pmatrix}$ では $\Delta=1$ なので、逆行列は $\begin{pmatrix}2&-1\\-3&2\end{pmatrix}$ である。この $\Delta$ は [行列式](/math/linear-algebra/determinants/) の章で扱う二次の行列式である。
:::

### 三角行列の逆

:::proposition[三角行列の正則性]{#prop-linalg-triangular-inverse}
上三角行列が正則であるための必要十分条件は、対角成分がすべて $0$ でないことである。正則な上三角行列の逆行列は上三角行列であり、対角成分は元の対角成分の逆数である。下三角行列にも同じ結論が成り立つ。
:::

:::proof
$U$ を上三角行列とする。すべての対角成分が $0$ でなければ、任意の $b$ について $Ux=b$ は後退代入で一意に解ける。従って $U$ は正則である。

逆に $0$ である対角成分があり、その添字の最小を $k$ とする。$x_k=1$、$x_j=0$（$j>k$）とおけば、第 $k$ 行以下の斉次方程式は満たされる。第 $k-1$ 行から第 1 行までは対角成分が $0$ でないので、順に $x_{k-1},\dots,x_1$ を決められる。こうして $Ux=0$ の解 $x\ne0$ が得られるため、$U$ は正則でない。

正則な場合、逆行列の第 $j$ 列は $Ux=e_j$ の解である。第 $n$ 行から第 $j+1$ 行までは右辺が $0$ なので $x_n=\cdots=x_{j+1}=0$ であり、第 $j$ 行で $x_j=1/u_{jj}$ を得る。これが逆行列の三角性と対角成分を与える。下三角の場合は転置して上三角の場合を適用し、:ref[prop-linalg-inverse-basic] の 3 を使う。
:::

:::example[三角系で逆行列の一列を求める]
実行列 $U=\begin{pmatrix}1&2&1\\0&2&1\\0&0&1\end{pmatrix}$ に対する $Ux=e_3$ では、$x_3=1$、$2x_2+1=0$、$x_1+2x_2+1=0$ の順になる。従って $U^{-1}$ の第 3 列は $\transp{(0,-1/2,1)}$ である。同じ代入を $e_1,e_2$ に行うと第 1 列は $\transp{(1,0,0)}$、第 2 列は $\transp{(-1,1/2,0)}$ となる。
:::

### 行交換を含む LU 分解

同じ係数行列で右辺を何度も取り替えるときは、前進消去で使う交換と加算を行列として保存できる。各行・各列に $1$ が一つずつあり、他の成分が $0$ である正方行列を :term[置換行列]{en="permutation matrix" yomi="ちかんぎょうれつ"} という。左から掛けると各行を一度ずつ選んで並べ替える。単位行列の行を交換しても置換行列の条件が保たれるため、行交換の積は置換行列である。置換行列の転置を掛けると元の順に戻る。実際、異なる行には同じ位置に $1$ が立たず、各行には $1$ が一つあるから $P\transp{P}=I$ であり、列についても同様に $\transp{P}P=I$ である。

:::theorem[行交換を含む LU 分解]{#thm-linalg-plu-factorization}
正則な $A\in M_n(K)$ に対し、置換行列 $P$、対角成分がすべて $1$ の下三角行列 $L$、正則な上三角行列 $U$ が存在して $PA=LU$ となる。$P$ を固定したとき、この条件を満たす $L,U$ は一意である。
:::

:::proof
$n=1$ のとき $P=L=I_1$、$U=A$ でよい。$n\ge2$ とし、より小さいサイズでの存在を仮定する。第 1 列が零なら $Ae_1=0$ となるので、正則性から第 1 列には $0$ でない成分がある。行交換 $P_0$ でその行を先頭に移し、

$$
P_0A=\begin{pmatrix}p&u\\v&B\end{pmatrix},
\qquad
S=B-vp^{-1}u
$$

と書く。ここで $p\ne0$ はスカラー、$u$ は行、$v$ は列である。$Sz=0$ を満たす $z\ne0$ があれば、$P_0A\begin{pmatrix}-p^{-1}uz\\z\end{pmatrix}=0$ となって $P_0A$ の正則性に反する。従って $S$ は正則であり、帰納法の仮定により $P_1S=L_1U_1$ と分解できる。

残りの行に $P_1$ を適用し、$P=\begin{pmatrix}1&0\\0&P_1\end{pmatrix}P_0$ とすると、ブロック積により

$$
PA=\begin{pmatrix}p&u\\P_1v&P_1B\end{pmatrix}
=
\underbrace{\begin{pmatrix}1&0\\P_1vp^{-1}&L_1\end{pmatrix}}_{L}
\underbrace{\begin{pmatrix}p&u\\0&U_1\end{pmatrix}}_{U}
$$

となる。右下のブロックは $P_1vp^{-1}u+L_1U_1=P_1B$ である。$L$ は対角成分が $1$ の下三角行列であり、$U$ は対角成分が $0$ でない上三角行列だから正則である。

同じ $P$ で $LU=L'U'$ なら $(L')^{-1}L=U'U^{-1}$ である。左辺は対角成分が $1$ の下三角行列、右辺は上三角行列である（:ref[prop-linalg-triangular-product]、:ref[prop-linalg-triangular-inverse]）。従って両辺は単位行列であり、$L=L'$、$U=U'$ である。
:::

:::example[分解してから右辺を代入する]
前進消去の例で用いた $A=\begin{pmatrix}0&2&1\\1&1&0\\2&3&1\end{pmatrix}$ では、第 1 行と第 2 行を交換する置換行列を $P$ とすると、

$$
PA=\begin{pmatrix}1&1&0\\0&2&1\\2&3&1\end{pmatrix},\qquad
L=\begin{pmatrix}1&0&0\\0&1&0\\2&\tfrac12&1\end{pmatrix},\qquad
U=\begin{pmatrix}1&1&0\\0&2&1\\0&0&\tfrac12\end{pmatrix}.
$$

$LU$ の第 3 行は $2(1,1,0)+\frac12(0,2,1)+(0,0,\frac12)=(2,3,1)$ であり、第 1、2 行は $U$ の対応する行である。従って $PA=LU$ である。

$b=\transp{(5,3,9)}$ に対し、まず $Ly=Pb=\transp{(3,5,9)}$ を上から順に解く。この手順を :term[前進代入]{en="forward substitution" yomi="ぜんしんだいにゅう"} という。ここでは $y_1=3$、$y_2=5$、$y_3=9-2\cdot3-5/2=1/2$ である。次に $Ux=y$ を後退代入で解けば $x=\transp{(1,2,1)}$ を得る。別の $b$ に対しても、同じ $P,L,U$ を使ってまず $Ly=Pb$ を解き、その解を $Ux=y$ の右辺に使う。
:::

:::proposition[行交換なしで LU 分解できる条件]{#prop-linalg-lu-leading-blocks}
正則な $A\in M_n(K)$ について、$A=LU$ と分解できることは、すべての $k=1,\dots,n$ について $A$ の先頭 $k$ 行・$k$ 列からなる正方ブロックが正則であることと同値である。ここで $L$ は対角成分が $1$ の下三角行列、$U$ は正則な上三角行列とする。
:::

:::proof
$A=LU$ と書ければ、先頭の $k\times k$ ブロックは $L$ と $U$ の対応するブロックどうしの積になる。積の各成分で $k$ より大きい添字は、$L$ の上三角側の零成分のため寄与しない。二つの三角ブロックの対角成分は非零なので、その積は正則である。

逆に、先頭のすべての正方ブロックが正則であるとする。まず $A$ の第 1 成分は非零なので、第 1 列の主成分として使える。行交換を行わずに最初の $k-1$ 列を前進消去できたとする。第 $k$ 行の先頭 $k-1$ 成分は $0$ である。この第 $k$ 行の第 $k$ 成分が $0$ なら、現在の先頭 $k\times k$ ブロックは最終行が零なので正則でない。ところが元の先頭 $k\times k$ ブロックに同じ消去を施して現在のブロックが得られ、途中で上の $k$ 行に外から足された行はない。元のブロックは正則であり、行基本変形で正則性は保たれるから矛盾である。よって第 $k$ 成分を主成分として使える。

こうして行交換なしで上三角行列 $U$ に到達する。各消去は対角成分が $1$ の下三角基本行列を左から掛ける操作である。その積を $E$ とすれば $EA=U$ であり、:ref[prop-linalg-triangular-product] と :ref[prop-linalg-triangular-inverse] により $L=E^{-1}$ も対角成分が $1$ の下三角行列である。$U$ の対角成分はすべて主成分なので非零であり、$A=LU$ が求める分解である。
:::

上の例では $A$ の先頭 $1\times 1$ ブロックが $(0)$ なので、そのままでは $A=LU$ に分解できない。交換後の $PA$ で分解できたのは、先頭の各正方ブロックが正則になったためである。

### ブロック消去と Schur 補行列

変数の組をまとめて消去できると、行列を分割して方程式を解ける。ここでは左上の正方ブロックが正則な場合を扱う。

:::proposition[Schur 補行列による可解条件]{#prop-linalg-schur-elimination}
$M=\begin{pmatrix}A&B\\C&D\end{pmatrix}$ を正方行列とし、左上の正方ブロック $A$ は正則とする。$S=D-CA^{-1}B$ を $A$ に関する :term[Schur 補行列]{en="Schur complement" yomi="しゅあほぎょうれつ"} と呼ぶ。$M\begin{pmatrix}x\\y\end{pmatrix}=\begin{pmatrix}b\\c\end{pmatrix}$ の全解は

$$
Sy=c-CA^{-1}b,\qquad x=A^{-1}(b-By)
$$

で与えられる。特に $M$ が正則であることと $S$ が正則であることは同値であり、正則な場合は

$$
M^{-1}=
\begin{pmatrix}
A^{-1}+A^{-1}BS^{-1}CA^{-1}&-A^{-1}BS^{-1}\\
-S^{-1}CA^{-1}&S^{-1}
\end{pmatrix}
$$

である。
:::

:::proof
上側の式 $Ax+By=b$ は $x=A^{-1}(b-By)$ と同値である。これを下側の式 $Cx+Dy=c$ に代入すると $CA^{-1}b+(D-CA^{-1}B)y=c$ となる。逆にこの $y$ とそれに対応する $x$ は両式を満たす。

斉次系では $Sy=0$ の解 $y$ に対し $x=-A^{-1}By$ が一意に決まる。従って $M$ の斉次系に零以外の解がないことと、$S$ の斉次系に零以外の解がないことは同値である。二つとも正方行列なので :ref[thm-linalg-invertible-equiv] が適用できる。$S$ が正則なら上の式から $y=S^{-1}c-S^{-1}CA^{-1}b$ が得られ、その式を $x=A^{-1}b-A^{-1}By$ に代入して $b,c$ の係数を集めると、表示された逆行列になる。
:::

:::example[二つの変数をまとめて消去する]
$M=\begin{pmatrix}2&0&1\\0&1&2\\1&1&4\end{pmatrix}$ の左上 $2\times2$ ブロックを $A=\diag(2,1)$ とする。実数の範囲では $A^{-1}=\diag(1/2,1)$ であり、

$$
S=4-(1\ \ 1)\begin{pmatrix}\tfrac12&0\\0&1\end{pmatrix}
\begin{pmatrix}1\\2\end{pmatrix}
=4-\tfrac52=\tfrac32.
$$

右辺を $\transp{(3,4,7)}$ とすれば、$Sy=7-(3/2+4)=3/2$ から $y=1$ である。上の二つの値は $x=\diag(1/2,1)\transp{(3-1,4-2)}=\transp{(1,2)}$ となる。
:::

:::corollary[ブロック消去と階数]{#cor-linalg-schur-rank}
左上の正方ブロック $A$ が正則なとき、$M=\begin{pmatrix}A&B\\C&D\end{pmatrix}$ の階数は $\rank A+\rank(D-CA^{-1}B)$ である。全体の行数と列数は異なっていてもよい。
:::

:::proof
対角ブロックに適切なサイズの単位行列を用いると

$$
\begin{pmatrix}I&0\\-CA^{-1}&I\end{pmatrix}
\begin{pmatrix}A&B\\C&D\end{pmatrix}
\begin{pmatrix}I&-A^{-1}B\\0&I\end{pmatrix}
=\begin{pmatrix}A&0\\0&D-CA^{-1}B\end{pmatrix}.
$$

左右の因子は、それぞれ非対角ブロックの符号を変えた行列を逆行列にもつ。従って :ref[cor-linalg-rank-calculus] により右辺の階数は $\rank M$ に等しい。$A$ 側と $D-CA^{-1}B$ 側を別々に簡約化すると、$A$ は単位行列になり、右下の主成分はすべてその右にある。零行は末尾に集まるので、全体の主成分の個数は両ブロックの階数の和である。
:::

例えば実行列 $\begin{pmatrix}2&1\\4&2\end{pmatrix}$ の左上ブロックを $(2)$ とすると補行列は $(2-4\cdot\frac12\cdot1)=(0)$ であり、全体の階数は $1+0=1$ である。

## 演習問題

::::exercise{level=basic}
$K=\R$ または $\C$ とする。$P=\begin{pmatrix}1&2\\3&4\end{pmatrix}$ と可換な（$PX=XP$ を満たす）2 次正方行列 $X$ をすべて求めよ。
:::solution
$X=\begin{pmatrix}a&b\\c&d\end{pmatrix}$ とおくと

$$
PX=\begin{pmatrix}a+2c & b+2d\\ 3a+4c & 3b+4d\end{pmatrix},\qquad
XP=\begin{pmatrix}a+3b & 2a+4b\\ c+3d & 2c+4d\end{pmatrix}.
$$

成分を比較すると $2c=3b$, $b+2d=2a+4b$, $3a+4c=c+3d$, $3b+4d=2c+4d$ である。第 1 式と第 4 式は同じ条件 $c=\tfrac32 b$ であり、第 2 式は $d = a + \tfrac32 b$、第 3 式に代入すると $3a + 3c = 3d$ すなわち $d=a+c=a+\tfrac32b$ となって第 2 式と一致する。よって $a, b$ を任意として

$$
X=\begin{pmatrix}a & b\\ \tfrac32 b & a+\tfrac32 b\end{pmatrix} = \Bigl(a-\tfrac12 b\Bigr) I_2 + \frac{b}{2}P\qquad(a,b\in K).
$$

また $I_2$ と $P$ はどちらも $P$ と可換なので、その一次結合はすべて条件を満たす。
:::
::::

::::exercise{level=basic}
次の連立一次方程式を解け。

$$
\begin{cases}
x + y + z = 1\\
x + 2y + 3z = 2\\
2x + 3y + 4z = 3
\end{cases}
$$
:::solution
拡大係数行列を簡約化する。第 2 行から第 1 行を、第 3 行から第 1 行の 2 倍を引き、さらに第 3 行から第 2 行を、第 1 行から第 2 行を引くと

$$
\left(\begin{array}{ccc|c}1&1&1&1\\1&2&3&2\\2&3&4&3\end{array}\right)
\to
\left(\begin{array}{ccc|c}1&0&-1&0\\0&1&2&1\\0&0&0&0\end{array}\right).
$$

$\rank A=\rank(A\mid b)=2$ なので解が存在し、自由変数は $z$ の 1 個である。$z=t$ とおいて

$$
(x,y,z) = (t,\ 1-2t,\ t)\qquad(t\in K).
$$
:::
::::

::::exercise{level=basic}
$a\in\R$ をパラメータとして、行列 $C=\begin{pmatrix}1&1&1\\1&2&a\\1&a&4\end{pmatrix}$ の階数を求めよ。
:::solution
第 2 行・第 3 行から第 1 行を引き、さらに第 3 行から第 2 行の $(a-1)$ 倍を引くと

$$
\begin{pmatrix}1&1&1\\0&1&a-1\\0&a-1&3\end{pmatrix}
\to
\begin{pmatrix}1&1&1\\0&1&a-1\\0&0&3-(a-1)^2\end{pmatrix}.
$$

$3-(a-1)^2=0 \iff a=1\pm\sqrt3$ である。よって $a=1\pm\sqrt3$ のとき $\rank C=2$、それ以外のとき $\rank C=3$。
:::
::::

::::exercise{level=basic}
$K=\R$ または $\C$ とし、$B=\begin{pmatrix}2&1&1\\1&2&1\\1&1&2\end{pmatrix}$ の逆行列を求めよ。
:::hint
$J$ をすべての成分が $1$ の 3 次正方行列とすると $B=I_3+J$, $J^2=3J$ である。$B^{-1}=\alpha I_3+\beta J$ の形を仮定してもよい。
:::
:::solution
ヒントの記法で $(I_3+J)(\alpha I_3+\beta J) = \alpha I_3 + (\alpha+\beta+3\beta)J$ である。これが $I_3$ となるには $\alpha=1$, $\alpha+4\beta=0$ すなわち $\beta=-\tfrac14$ であればよい。$I_3+J$ と $\alpha I_3+\beta J$ は可換なので、:ref[thm-linalg-invertible-equiv] を持ち出すまでもなく両側の逆行列であり、

$$
B^{-1} = I_3 - \frac14 J = \frac14\begin{pmatrix}3&-1&-1\\-1&3&-1\\-1&-1&3\end{pmatrix}.
$$

もちろん $(B\mid I_3)$ の簡約化でも同じ結果が得られる。
:::
::::

::::exercise{level=advanced}
$N\in M_n(K)$ がある正の整数 $k$ について $N^k=O$ を満たす（**冪零**である）とき、$I_n-N$ は正則であって

$$
(I_n-N)^{-1} = I_n + N + N^2 + \cdots + N^{k-1}
$$

となることを示せ。また $N=\begin{pmatrix}0&1&2\\0&0&3\\0&0&0\end{pmatrix}$ のとき $(I_3-N)^{-1}$ を求めよ。
:::solution
$S\defeq I_n+N+\cdots+N^{k-1}$ とおくと、

$$
(I_n-N)S = S - NS = (I_n+N+\cdots+N^{k-1}) - (N+N^2+\cdots+N^{k}) = I_n - N^k = I_n
$$

であり、$S(I_n-N)=I_n$ も同様である。よって $I_n-N$ は正則で $S$ がその逆行列である。

与えられた $N$ については $N^2=\begin{pmatrix}0&0&3\\0&0&0\\0&0&0\end{pmatrix}$, $N^3=O$ なので

$$
(I_3-N)^{-1} = I_3+N+N^2 = \begin{pmatrix}1&1&5\\0&1&3\\0&0&1\end{pmatrix}.
$$
:::
::::

::::exercise{level=advanced}
正方行列 $A=(a_{ij})$ の対角成分の和 $\tr A \defeq \sum_i a_{ii}$ を :term[トレース]{en="trace" yomi="とれーす"} という。

1. $A\in M_{m,n}(K)$, $B\in M_{n,m}(K)$ に対し $\tr(AB)=\tr(BA)$ を示せ。
2. $K=\R$ または $\C$ のとき、$AB-BA=I_n$ を満たす $A, B\in M_n(K)$ は存在しないことを示せ。
3. $m\ne n$ のとき、$A\in M_{m,n}(K)$ が $BA=I_n$ と $AC=I_m$ を同時に満たす $B, C$ をもつことはないことを示せ。
:::solution
1. $\tr(AB)=\sum_{i=1}^m\sum_{j=1}^n a_{ij}b_{ji} = \sum_{j=1}^n\sum_{i=1}^m b_{ji}a_{ij} = \tr(BA)$。

2. 存在したとすると、1 より $\tr(AB-BA)=\tr(AB)-\tr(BA)=0$ だが $\R$ または $\C$ では $\tr I_n = n\ne 0$ なので矛盾する。

3. :ref[prop-linalg-one-sided-inverses] により $BA=I_n$ から $\rank A=n$、$AC=I_m$ から $\rank A=m$ となる。階数は行数と列数で上から抑えられる整数なので、任意の体で $m=n$ を得る。
:::
::::

::::exercise{level=basic}
実数の範囲で次の二つの式を同じ係数行列を用いて解け。

$$
\begin{pmatrix}0&1&1\\1&2&0\\2&4&0\end{pmatrix}x
=\begin{pmatrix}3\\1\\2\end{pmatrix},
\qquad
\begin{pmatrix}0&1&1\\1&2&0\\2&4&0\end{pmatrix}x
=\begin{pmatrix}3\\1\\4\end{pmatrix}.
$$

解が存在しない場合は、元の方程式から矛盾を作る行の組み合わせを示せ。
:::solution
第 3 行から第 2 行の 2 倍を引き、第 1 行と第 2 行を交換する。二つの右辺をまとめた行列は

$$
\left(\begin{array}{ccc|cc}0&1&1&3&3\\1&2&0&1&1\\2&4&0&2&4\end{array}\right)
\longrightarrow
\left(\begin{array}{ccc|cc}1&2&0&1&1\\0&1&1&3&3\\0&0&0&0&2\end{array}\right)
$$

となる。第 1 の右辺では $x_3=t$ とおけば $x_2=3-t$、$x_1=1-2x_2=-5+2t$ となる。各 $t\in\R$ で全行を満たし、どの解もこの形である。第 2 の右辺では、元の第 3 式から第 2 式の 2 倍を引くと $0=2$ となるため解は存在しない。
:::
::::

::::exercise{level=basic}
$A=\begin{pmatrix}1&2\\0&1\\2&0\end{pmatrix}$ の左逆行列をすべて求めよ。右逆行列の有無も答えよ。
:::solution
$LA=I_2$ の両辺を転置すると $\transp{A}\transp{L}=I_2$ である。$L$ の第 1 行を $(a,b,s)$ とすると $a+2s=1$、$2a+b=0$ なので $a=1-2s$、$b=-2+4s$ を得る。第 2 行を $(c,d,t)$ とすると $c+2t=0$、$2c+d=1$ なので $c=-2t$、$d=1+4t$ を得る。従って

$$
L=\begin{pmatrix}1-2s&-2+4s&s\\-2t&1+4t&t\end{pmatrix},
\qquad s,t\in K
$$

が全解である。これを $A$ に掛けると、各行の条件から $I_2$ になる。$A$ の階数は列数 $2$ を超えられないので、行数 $3$ が必要な右逆行列は存在しない。
:::
::::

::::exercise{level=advanced}
$A\in M_{m,n}(K)$ とする。$A$ に解がある右辺全体と $A'$ に解がある右辺全体が一致することと、$A$ と $A'\in M_{m,n}(K)$ が列基本変形で互いに移り合うことが同値であることを示せ。
:::solution
$A'=AF$ と書け、$F$ が列基本変形を表す基本行列の積であるとする。このとき $x=Fy$ と $y=F^{-1}x$ により、$Ax=b$ の可解性と $A'y=b$ の可解性は一致する。

逆に同じ右辺に対して解をもつとする。行ベクトル $v$ が $vA=0$ を満たすなら、$A'$ の各列 $a'_j$ は $A$ を係数行列とする方程式の可解な右辺だから、ある $x_j$ によって $a'_j=Ax_j$ と書ける。よって $va'_j=0$ であり $vA'=0$ である。$A,A'$ を取り替えれば逆が成り立つので、$\transp{A}$ と $\transp{A'}$ の斉次方程式は同じ解をもつ。

行操作は斉次解を保つので、$\transp{A}$ と $\transp{A'}$ の簡約行階段形も同じ斉次解をもつ。:ref[thm-linalg-rref-unique] により二つの簡約行階段形は等しく、$\transp{A}$ と $\transp{A'}$ は行同値である。各操作を転置すると、$A$ と $A'$ は列基本変形で互いに移り合う。
:::
::::

::::exercise{level=advanced}
実数の範囲で $M=\begin{pmatrix}2&1&1\\0&1&2\\1&1&4\end{pmatrix}$ の左上を $2\times2$ ブロックとして、Schur 補行列を用い $M^{-1}$ を求めよ。
:::solution
$A=\begin{pmatrix}2&1\\0&1\end{pmatrix}$、$B=\transp{(1,2)}$、$C=(1\ \ 1)$、$D=(4)$ とする。$A^{-1}=\begin{pmatrix}1/2&-1/2\\0&1\end{pmatrix}$ なので $A^{-1}B=\transp{(-1/2,2)}$、$CA^{-1}=(1/2\ \ 1/2)$、$S=4-(-1/2+2)=5/2$ である。

右辺を $\transp{(b_1,b_2,c)}$ とすれば、:ref[prop-linalg-schur-elimination] により

$$
\begin{aligned}
y&=\tfrac25(c-\tfrac12b_1-\tfrac12b_2),\\
x_1&=\tfrac12(b_1-b_2)+\tfrac12y=\tfrac25b_1-\tfrac35b_2+\tfrac15c,\\
x_2&=b_2-2y=\tfrac25b_1+\tfrac75b_2-\tfrac45c.
\end{aligned}
$$

よって $M^{-1}=\tfrac15\begin{pmatrix}2&-3&1\\2&7&-4\\-1&-1&2\end{pmatrix}$ である。
:::
::::

## 確認テスト

::::quiz
$K=\R$ とする。次のうち、行基本変形**でない**ものはどれか。

- [ ] 第 1 行と第 3 行を入れ替える
- [x] 第 2 行を $0$ 倍する
- [ ] 第 2 行に第 1 行の $-5$ 倍を加える
- [ ] 第 3 行を $\tfrac12$ 倍する

:::explanation
:ref[def-linalg-row-operations] により、異なる行の交換、他の行の定数倍の加算、非零定数による行のスカラー倍はすべて行基本変形である。$0$ 倍だけは逆操作が存在しない。例えば $I_3$ の第 2 行を $0$ 倍すると階数が $3$ から $2$ に変わるので、得られた行列は $I_3$ と行同値でない。
:::
::::

::::quiz{answer="2" tol="1e-6"}
実行列 $\begin{pmatrix}1&2&3\\2&4&6\\1&0&1\end{pmatrix}$ の階数を求めよ。

:::explanation
実数の範囲で、第 2 行から第 1 行の 2 倍を引いて零行にし、第 3 行から第 1 行を引いて $(0,-2,-2)$ を得る。零行を末尾に移すと、主成分の列が第 1 列と第 2 列の行階段形になる。したがって階数は $2$ である。
:::
::::

::::quiz
$A\in M_{3,5}(\R)$ とする。正しいものをすべて選べ。

- [x] 斉次方程式 $Ax=0$ は自明でない解をもつ
- [ ] 任意の $b\in\R^3$ に対し $Ax=b$ は解をもつ
- [x] $\rank A\le 3$
- [ ] $Ax=0$ の一般解はちょうど 2 個の任意定数を含む

:::explanation
$\rank A\le\min\{3,5\}=3<5$ なので :ref[cor-linalg-homogeneous-nontrivial] により自明でない解をもつ。$\rank A$ が 3 未満のこともある（例えば $A=O$）ので、$Ax=b$ が常に解けるとは限らず、任意定数の個数 $5-\rank A$ も 2 とは限らない。
:::
::::

::::quiz{answer="2" tol="1e-6"}
:ref[exm-linalg-solve] の連立方程式の一般解に含まれる任意定数の個数はいくつか。

:::explanation
未知数は 4 個、係数行列の階数は 2 なので、:ref[thm-linalg-rank-criterion] より $4-2=2$ 個である。
:::
::::

::::quiz
$A, B\in M_n(K)$ が $AB=I_n$ を満たすとする。正しいものを選べ。

- [ ] $BA=I_n$ とは限らない
- [x] $BA=I_n$ が成り立ち、$B=A^{-1}$ である
- [ ] $A=B$ である
- [ ] $A$ と $B$ はともに対称行列である

:::explanation
:ref[thm-linalg-invertible-equiv] の (7)$\Rightarrow$(1) により $A$ は正則である。$AB=I_n$ の左から $A^{-1}$ を掛ければ $B=A^{-1}$ となり、従って $BA=I_n$ である。また、実数の範囲で $A=\begin{pmatrix}1&1\\0&1\end{pmatrix}$、$B=\begin{pmatrix}1&-1\\0&1\end{pmatrix}$ は $AB=I_2$ を満たすが、$A$ と $B$ は異なり、どちらも対称ではない。
:::
::::

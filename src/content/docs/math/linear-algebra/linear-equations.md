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
  - symbol: 'A^{-1}'
    meaning: 正則行列 A の逆行列
---

## この章の目標

- 行列の和・スカラー倍・積を定義し、その代数的性質を理解する。
- 行基本変形によって任意の行列を簡約行階段形に変形でき、しかもその結果が一意であることを証明する。
- 連立一次方程式の解の存在と自由度を、行列の階数によって完全に記述する。
- 正則行列のさまざまな同値な特徴付けを証明し、逆行列を計算できるようになる。

この章では、係数を一般の**体**（加減乗除が自由にできる数の体系）$K$ にとる。具体的には $K=\R$（実数体）または $K=\C$（複素数体）と思って読んでよい。体の公理は次章 [ベクトル空間・基底・次元](/math/linear-algebra/vector-spaces/) で正確に述べる。

## 行列とその演算

### 行列の定義

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

横に並んだ $(a_{i1},\dots,a_{in})$ を $A$ の第 $i$ **行**、縦に並んだ成分を第 $j$ **列**という。$m\times 1$ 行列を**列ベクトル**、$1\times n$ 行列を**行ベクトル**と呼ぶ。本サイトでは $K^n$ を $n$ 次列ベクトル全体の集合とし、その元を

$$
x = \begin{pmatrix} x_1\\ \vdots\\ x_n\end{pmatrix}
$$

と書く。2 つの行列 $A=(a_{ij})$, $B=(b_{ij})$ は、サイズが等しく、かつすべての $i,j$ について $a_{ij}=b_{ij}$ であるときに限り**等しい**といい、$A=B$ と書く。

すべての成分が $0$ である行列を**零行列**といい $O$（サイズを明示するときは $O_{m,n}$）で表す。正方行列 $A=(a_{ij})$ の成分 $a_{11},a_{22},\dots,a_{nn}$ を**対角成分**という。対角成分がすべて $1$ でそれ以外の成分がすべて $0$ である $n$ 次正方行列を $n$ 次 :term[単位行列]{en="identity matrix" yomi="たんいぎょうれつ"} といい $I_n$ と書く。Kronecker のデルタ

$$
\delta_{ij} = \begin{cases} 1 & (i=j)\\ 0 & (i\ne j)\end{cases}
$$

を用いれば $I_n=(\delta_{ij})$ である。$I_n$ の第 $j$ 列を $e_j$ と書き、**基本ベクトル**と呼ぶ。

### 和・スカラー倍・積

:::definition[行列の和とスカラー倍]{#def-linalg-matrix-sum}
$A=(a_{ij}), B=(b_{ij})\in M_{m,n}(K)$ と $c\in K$ に対し、和 $A+B$ とスカラー倍 $cA$ を

$$
A+B \defeq (a_{ij}+b_{ij}),\qquad cA \defeq (c\,a_{ij})
$$

で定める。また $-A\defeq(-1)A$, $A-B\defeq A+(-B)$ とする。
:::

和とスカラー倍については、数の場合と同様の計算法則（結合法則・交換法則・分配法則など）が成り立つ。これらは成分ごとの等式に帰着するので証明は容易である。

:::definition[行列の積]{#def-linalg-matrix-product}
$A=(a_{ij})\in M_{l,m}(K)$, $B=(b_{jk})\in M_{m,n}(K)$ に対し、積 $AB\in M_{l,n}(K)$ を

$$
(AB)_{ik} \defeq \sum_{j=1}^{m} a_{ij}b_{jk}\qquad(1\le i\le l,\ 1\le k\le n) \label{eq-linalg-product}
$$

で定める。積 $AB$ は、$A$ の列の個数と $B$ の行の個数が等しいときに限り定義される。
:::

式 :eqref[eq-linalg-product] は「$A$ の第 $i$ 行と $B$ の第 $k$ 列の対応する成分の積の和」である。特に $A\in M_{m,n}(K)$ と $x\in K^n$ の積 $Ax\in K^m$ は

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
1 のみ示す（他は成分を比較すれば直ちに従う）。$A\in M_{k,l}(K)$, $B\in M_{l,m}(K)$, $C\in M_{m,n}(K)$ とする。$(i,q)$ 成分を比較すると

$$
((AB)C)_{iq} = \sum_{r=1}^{m}\Bigl(\sum_{j=1}^{l} a_{ij}b_{jr}\Bigr)c_{rq}
= \sum_{j=1}^{l}\sum_{r=1}^{m} a_{ij}b_{jr}c_{rq}
= \sum_{j=1}^{l} a_{ij}\Bigl(\sum_{r=1}^{m} b_{jr}c_{rq}\Bigr) = (A(BC))_{iq}.
$$

有限和の順序交換と、体における分配法則・結合法則を用いた。
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

:::definition[転置行列]{#def-linalg-transpose}
$A=(a_{ij})\in M_{m,n}(K)$ に対し、$(j,i)$ 成分が $a_{ij}$ である $n\times m$ 行列を $A$ の :term[転置行列]{en="transpose" yomi="てんちぎょうれつ"} といい、$\transp{A}$ と書く。
:::

:::proposition{#prop-linalg-transpose-laws}
$(A+B)^\top = \transp{A}+\transp{B}$, $(cA)^\top = c\,\transp{A}$, $(\transp{A})^\top = A$, $(AB)^\top = \transp{B}\,\transp{A}$ が成り立つ。
:::

:::proof
最後の式のみ示す。$A\in M_{l,m}(K)$, $B\in M_{m,n}(K)$ とすると、$((AB)^\top)_{ki} = (AB)_{ik} = \sum_j a_{ij}b_{jk} = \sum_j (\transp{B})_{kj}(\transp{A})_{ji} = (\transp{B}\transp{A})_{ki}$ である。
:::

行列を縦横の区切り線で小さな行列（**ブロック**）に分けて扱うと便利なことが多い。例えば $A=\begin{pmatrix}A_{11}&A_{12}\\A_{21}&A_{22}\end{pmatrix}$, $B=\begin{pmatrix}B_{11}&B_{12}\\B_{21}&B_{22}\end{pmatrix}$ と分割し、各積 $A_{ij}B_{jk}$ が定義されるように区切り方が整合していれば、

$$
AB = \begin{pmatrix}A_{11}B_{11}+A_{12}B_{21} & A_{11}B_{12}+A_{12}B_{22}\\ A_{21}B_{11}+A_{22}B_{21} & A_{21}B_{12}+A_{22}B_{22}\end{pmatrix}
$$

が成り立つ。これは積の定義 :eqref[eq-linalg-product] の和を区切りに沿って分けたものにすぎない。

### 連立一次方程式の行列表示

$n$ 個の未知数 $x_1,\dots,x_n$ に関する $m$ 本の方程式

$$
\begin{cases}
a_{11}x_1 + \cdots + a_{1n}x_n = b_1\\
\qquad\vdots\\
a_{m1}x_1 + \cdots + a_{mn}x_n = b_m
\end{cases}
$$

を**連立一次方程式**という。$A=(a_{ij})\in M_{m,n}(K)$ を**係数行列**、$b=\transp{(b_1,\dots,b_m)}$ とおけば、この方程式は

$$
Ax = b \label{eq-linalg-system}
$$

と書ける。$b=0$ のとき**斉次**（同次）であるという。$A$ の右に $b$ を付け加えた $m\times(n+1)$ 行列 $(A\mid b)$ を :term[拡大係数行列]{en="augmented matrix" yomi="かくだいけいすうぎょうれつ"} と呼ぶ。

## 基本変形と行階段形

### 行基本変形と基本行列

連立一次方程式を解く消去法の操作を、行列の言葉で定式化する。

:::definition[行基本変形]{#def-linalg-row-operations}
行列に対する次の 3 種類の操作を :term[行基本変形]{en="elementary row operation" yomi="ぎょうきほんへんけい"} という。

1. 第 $i$ 行と第 $j$ 行を入れ替える（$i\ne j$）。
2. 第 $i$ 行を $c$ 倍する（$c\in K$, $c\ne 0$）。
3. 第 $i$ 行に第 $j$ 行の $c$ 倍を加える（$i\ne j$, $c\in K$）。

行列 $A$ に有限回の行基本変形を施して $B$ が得られるとき、$A$ と $B$ は**行同値**であるという。
:::

それぞれの行基本変形は、単位行列に同じ変形を施した行列を左から掛けることで実現できる。

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
$E=(e_{kl})$ を基本行列とすると、$(EA)$ の第 $k$ 行は $\sum_l e_{kl}\cdot(A\text{ の第 } l\text{ 行})$ である。$E=R_{ij}(c)$ の場合、$k\ne i$ なら $E$ の第 $k$ 行は $I_m$ の第 $k$ 行なので $EA$ の第 $k$ 行は $A$ の第 $k$ 行のままであり、$k=i$ なら $E$ の第 $i$ 行は $\transp{e_i}+c\,\transp{e_j}$ なので $EA$ の第 $i$ 行は「$A$ の第 $i$ 行 $+\,c\times$ 第 $j$ 行」となる。他の 2 種類も同様である。後半は、変形を元に戻す変形が同じ種類の基本変形であることから、前半を $I_m$ に適用すれば従う。
:::

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
$A\in M_{m,n}(K)$ に対し、列の個数 $n$ に関する帰納法で示す。より強く「$A$ の最初の $k$ 列からなる部分が簡約行階段形になるように変形できる」ことを $k=1,\dots,n$ について順に示せばよい。

$k$ 列目まで処理が済んでおり、最初の $k$ 列の部分に主成分が $r$ 個あるとする（このとき第 $r+1$ 行以下の最初の $k$ 成分はすべて $0$ である）。第 $k+1$ 列について、第 $r+1$ 行以下に $0$ でない成分がなければ、最初の $k+1$ 列の部分はすでに簡約行階段形である。そうでなければ、その成分を含む行を第 $r+1$ 行と入れ替え（変形 1）、第 $r+1$ 行をその成分の逆数倍して $(r+1,k+1)$ 成分を $1$ にする（変形 2）。さらに $i\ne r+1$ の各行から第 $r+1$ 行の $a_{i,k+1}$ 倍を引けば（変形 3）、第 $k+1$ 列は $e_{r+1}$ になる。この操作は、第 $r+1$ 行の最初の $k$ 成分がすべて $0$ であることから、最初の $k$ 列を変えない。よって最初の $k+1$ 列の部分が簡約行階段形になる。$k=n$ まで繰り返せばよい。
:::

この証明の手順は、そのまま計算のアルゴリズムになる。簡約行階段形は変形の手順によらず一意に定まる。

:::theorem[簡約行階段形の一意性]{#thm-linalg-rref-unique}
行同値な 2 つの簡約行階段形行列は等しい。したがって、各行列 $A$ に行同値な簡約行階段形はただ 1 つに定まる。これを $A$ の簡約行階段形と呼ぶ。
:::

証明には次の事実を用いる。

:::lemma{#lem-linalg-row-equiv-solutions}
$A$ と $B$ が行同値ならば、斉次方程式 $Ax=0$ と $Bx=0$ の解の集合は一致する。
:::

:::proof
:ref[lem-linalg-elementary-left-mult] より、ある基本行列の積 $E=E_k\cdots E_1$ によって $B=EA$ と書ける。$Ax=0$ ならば $Bx=EAx=0$ である。逆に $E$ は正則で $A=E^{-1}B$ だから、$Bx=0$ ならば $Ax=0$ である。
:::

:::proof[:ref[thm-linalg-rref-unique] の証明]
$R, R'\in M_{m,n}(K)$ を行同値な簡約行階段形とし、$n$ に関する帰納法で $R=R'$ を示す。$n=1$ のとき、$1$ 列の簡約行階段形は $0$ か $e_1$ のどちらかであり、零行列は零行列としか行同値でないから $R=R'$ である。

$n\ge 2$ とする。$R, R'$ から最後の列を除いた行列 $R_0, R_0'$ もまた簡約行階段形であり、$R$ を $R'$ に変える行基本変形によって $R_0$ は $R_0'$ に変わるので、$R_0$ と $R_0'$ は行同値である。帰納法の仮定より $R_0=R_0'$ となる。したがって $R$ と $R'$ は最後の列だけが異なりうる。

$R\ne R'$ と仮定して矛盾を導く。$Rx=0$ を満たす任意の $x$ は、:ref[lem-linalg-row-equiv-solutions] により $R'x=0$ も満たすから $(R-R')x=0$ である。$R-R'$ は最後の列 $d\ne 0$ 以外が零なので $(R-R')x = x_n d$ であり、したがって $x_n=0$ である。つまり、$Rx=0$ の解はすべて $x_n=0$ を満たす。

もし $R$ の第 $n$ 列が主列でなければ、$x_n=1$ とし、他の非主列に対応する変数を $0$ として、主列に対応する変数を決めれば $Rx=0$ の解が得られる（主列は単位ベクトルなので、各主成分の行の方程式がちょうど 1 つの主列変数を決める）。これは $x_n=0$ に反する。よって第 $n$ 列は $R$ の主列であり、同様に $R'$ の主列でもある。$R_0=R_0'$ の零でない行の個数を $r$ とすれば、$R$ と $R'$ の最後の列はともに $e_{r+1}$ となり、$R=R'$ が従う。これは仮定に反する。
:::

## 連立一次方程式の解法

### 行基本変形と解の集合

:::proposition{#prop-linalg-row-ops-preserve-solutions}
拡大係数行列 $(A\mid b)$ に行基本変形を施して $(A'\mid b')$ が得られたとすると、$Ax=b$ と $A'x=b'$ の解の集合は一致する。
:::

:::proof
ある正則行列 $E$（基本行列の積）により $A'=EA$, $b'=Eb$ と書ける。$Ax=b$ ならば $A'x=EAx=Eb=b'$ であり、逆も $E^{-1}$ を掛ければよい。
:::

したがって、$(A\mid b)$ を簡約行階段形 $(R\mid c)$ に変形してから解けばよい。$R$ の主列に対応する変数を**主変数**、それ以外の変数を :term[自由変数]{en="free variable" yomi="じゆうへんすう"} と呼ぶ。

:::theorem[解の存在と表示]{#thm-linalg-solve-rref}
$(A\mid b)$ の簡約行階段形を $(R\mid c)$ とし、$R$ の零でない行の個数を $r$ とする。

1. $Ax=b$ が解をもつための必要十分条件は、$(R\mid c)$ の最後の列が主列でないこと、すなわち $c_{r+1}=\cdots=c_m=0$ である。
2. 解をもつとき、自由変数（$n-r$ 個）に任意の値を与えると、主変数はそれらの一次式として一意に定まる。こうしてすべての解が得られる。
:::

:::proof
1. 最後の列が主列であれば、ある行が $(0,\dots,0\mid 1)$ となり、方程式 $0=1$ を含むので解はない。主列でなければ、$R$ の第 $i$ 行（$i\le r$）の方程式は、主成分の列に対応する主変数 $x_{j_i}$ について

$$
x_{j_i} = c_i - \sum_{k:\,\text{自由変数}} r_{ik}x_k
$$

と解かれた形をしており、第 $r+1$ 行以下は $0=0$ である。

2. 上の式から明らかである。主列は単位ベクトルなので、各主変数はちょうど 1 つの式で決まる。
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

## 階数と解の構造

### 階数

:::definition[階数]{#def-linalg-rank}
行列 $A$ の簡約行階段形の零でない行の個数（すなわち主成分の個数）を $A$ の :term[階数]{en="rank" yomi="かいすう"} といい、$\rank A$ で表す。
:::

:ref[thm-linalg-rref-unique] によって簡約行階段形は一意なので、階数は well-defined である。また行同値な行列は同じ簡約行階段形をもつから、同じ階数をもつ。$A\in M_{m,n}(K)$ ならば明らかに $\rank A\le\min\{m,n\}$ である。階数のより概念的な意味（列ベクトルの張る空間の次元）は第 3 章で明らかになる。

:::theorem[解の存在と自由度]{#thm-linalg-rank-criterion}
$A\in M_{m,n}(K)$, $b\in K^m$ とする。

1. $Ax=b$ が解をもつ $\iff$ $\rank A = \rank (A\mid b)$。
2. 解をもつとき、一般解は $n-\rank A$ 個の任意定数を含む。特に、解がただ 1 つである $\iff$ $\rank A=n$。
:::

:::proof
$(A\mid b)$ の簡約行階段形 $(R\mid c)$ において、$R$ は $A$ の簡約行階段形である（$(A\mid b)$ を簡約化する操作は、最初の $n$ 列について $A$ を簡約化する操作になっている）。よって $\rank A$ は $R$ の主成分の個数、$\rank(A\mid b)$ は $(R\mid c)$ の主成分の個数であり、両者が等しいことは最後の列が主列でないことと同値である。あとは :ref[thm-linalg-solve-rref] から従う。
:::

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

## 正則行列と逆行列

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
:::

:::proof
(1)$\Rightarrow$(6), (1)$\Rightarrow$(7) は定義から明らかである。

(6)$\Rightarrow$(4)：$Ax=0$ ならば $x = BAx = 0$。

(4)$\Rightarrow$(2)：:ref[cor-linalg-homogeneous-nontrivial] による。

(2)$\Rightarrow$(3)：$n$ 次正方行列の簡約行階段形が $n$ 個の主成分をもてば、すべての列が主列であり、簡約行階段形の定義から $I_n$ に等しい。

(3)$\Rightarrow$(5)：基本行列 $E_1,\dots,E_k$ により $E_k\cdots E_1A=I_n$ と書ける。各 $E_i$ は正則なので $A = E_1^{-1}\cdots E_k^{-1}$ であり、:ref[lem-linalg-elementary-left-mult] より基本行列の逆行列はまた基本行列である。

(5)$\Rightarrow$(1)：基本行列は正則であり、正則行列の積は正則である（:ref[prop-linalg-inverse-basic]）。

(7)$\Rightarrow$(1)：$AB=I_n$ ならば $\transp{B}\,\transp{A}=I_n$ なので、$\transp{A}$ に対して (6) が成り立つ。すでに示した (6)$\Rightarrow$(1) により $\transp{A}$ は正則であり、:ref[prop-linalg-inverse-basic] の 3 から $A=(\transp{A})^\top$ も正則である。
:::

:::remark
(6) や (7) の $B$ は自動的に $A^{-1}$ に等しい。実際 $A$ が正則なら $B = BAA^{-1} = A^{-1}$（(6) の場合）である。正方でない行列については、片側の逆行列が存在しても両側の逆行列は存在しない（演習で確かめよ）。
:::

### 逆行列の計算

:ref[thm-linalg-invertible-equiv] の証明から、次の計算法が得られる。$n\times 2n$ 行列 $(A\mid I_n)$ に行基本変形を施して左半分を簡約化する。左半分が $I_n$ になれば $A$ は正則で、そのとき右半分が $A^{-1}$ である。左半分が $I_n$ にならなければ $A$ は正則でない。

実際、行基本変形の列に対応する基本行列の積を $E$ とすれば $E(A\mid I_n)=(EA\mid E)$ であり、$EA=I_n$ ならば $E=A^{-1}$ である（:ref[thm-linalg-invertible-equiv] の (6)）。

:::example{#exm-linalg-inverse}
$A=\begin{pmatrix}1&2&0\\0&1&1\\1&0&1\end{pmatrix}$ の逆行列を求める。

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

$AA^{-1}=I_3$ を直接確かめることができる。
:::

## 演習問題

::::exercise{level=basic}
$P=\begin{pmatrix}1&2\\3&4\end{pmatrix}$ と可換な（$PX=XP$ を満たす）2 次正方行列 $X$ をすべて求めよ。
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

すなわち $P$ と可換な行列は $I_2$ と $P$ の一次結合に限る。
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
$B=\begin{pmatrix}2&1&1\\1&2&1\\1&1&2\end{pmatrix}$ の逆行列を求めよ。
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
正方行列 $A=(a_{ij})$ の対角成分の和 $\tr A \defeq \sum_i a_{ii}$ を**トレース**という。

1. $A\in M_{m,n}(K)$, $B\in M_{n,m}(K)$ に対し $\tr(AB)=\tr(BA)$ を示せ。
2. $K=\R$ または $\C$ のとき、$AB-BA=I_n$ を満たす $A, B\in M_n(K)$ は存在しないことを示せ。
3. $m\ne n$ のとき、$A\in M_{m,n}(K)$ が $BA=I_n$ と $AC=I_m$ を同時に満たす $B, C$ をもつことはないことを示せ。
:::solution
1. $\tr(AB)=\sum_{i=1}^m\sum_{j=1}^n a_{ij}b_{ji} = \sum_{j=1}^n\sum_{i=1}^m b_{ji}a_{ij} = \tr(BA)$。

2. 存在したとすると、1 より $\tr(AB-BA)=\tr(AB)-\tr(BA)=0$ だが $\tr I_n = n\ne 0$ であり矛盾する（$\R$, $\C$ では $n\ne 0$ である。標数 $p$ の体では $p\mid n$ のとき反例がある）。

3. $BA=I_n$, $AC=I_m$ とすると $B = B(AC) = (BA)C = C$ である。1 より $n=\tr I_n=\tr(BA)=\tr(AB)=\tr(AC)=\tr I_m=m$ となり、$m\ne n$ に反する（ここでも $K=\R,\C$ を仮定した。一般の体では第 3 章の次元論を用いて示せる）。
:::
::::

## 確認テスト

::::quiz
次のうち、行基本変形**でない**ものはどれか。

- [ ] 第 1 行と第 3 行を入れ替える
- [x] 第 2 行を $0$ 倍する
- [ ] 第 2 行に第 1 行の $-5$ 倍を加える
- [ ] 第 3 行を $\tfrac12$ 倍する

:::explanation
行のスカラー倍は $c\ne0$ の場合に限る。$0$ 倍は元に戻せない（逆の変形が存在しない）ので、行同値性が保たれない。
:::
::::

::::quiz{answer="2"}
行列 $\begin{pmatrix}1&2&3\\2&4&6\\1&0&1\end{pmatrix}$ の階数を求めよ。

:::explanation
第 2 行は第 1 行の 2 倍なので、第 2 行から第 1 行の 2 倍を引くと零行になる。残る 2 行 $(1,2,3)$, $(1,0,1)$ は簡約化すると主成分を 2 つもつ。よって階数は 2。
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

::::quiz{answer="2"}
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
:ref[thm-linalg-invertible-equiv] の (7)$\Rightarrow$(1) により $A$ は正則であり、$B=A^{-1}AB=A^{-1}$ となる。正方行列に対しては片側の逆行列が両側の逆行列になる。
:::
::::

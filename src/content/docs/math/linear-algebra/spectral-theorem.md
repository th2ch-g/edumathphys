---
title: 正規作用素とスペクトル定理・二次形式
description: 自己随伴・ユニタリ・正規作用素を内積で特徴づけ、正規作用素のスペクトル定理と実対称行列の直交対角化を証明し、二次形式の主軸変換と慣性法則、正定値行列の判定・平方根・極分解へ進む。
notation:
  - symbol: '\adj{T}'
    meaning: 内積に関する T の随伴作用素
  - symbol: 'E_\lambda(T)'
    meaning: 固有値 λ に対する T の固有空間
  - symbol: 'P_\lambda'
    meaning: 正規作用素の固有空間 E_λ(T) への直交射影（スペクトル射影）
  - symbol: '\mathrm{U}(n),\ \mathrm{O}(n)'
    meaning: n 次ユニタリ群と n 次直交群
  - symbol: '\SU(n),\ \SO(n)'
    meaning: 行列式 1 の元からなる特殊ユニタリ群と特殊直交群
  - symbol: 'R_\theta'
    meaning: 平面の角 θ の回転行列
  - symbol: 'H>0,\ H\ge0'
    meaning: Hermite 行列 H が正定値・半正定値であること
  - symbol: '\Delta_k'
    meaning: 左上 k 次の首座小行列式
  - symbol: 'H^{1/2}'
    meaning: 半正定値行列 H の半正定値平方根
---

## この章の目標

- 自己随伴・ユニタリ・正規作用素を定義し、$\inner{Tv}{v}$ やノルムによる特徴づけ、ユニタリ群・直交群の基本性質を証明する。
- Schur の三角化を証明し、複素内積空間の作用素が正規直交固有基底をもつことと正規であることの同値性（スペクトル定理）を導く。
- 実対称行列の直交対角化、実正規行列と直交行列の標準形を導き、3 次元の回転の軸と角を求める。
- 二次形式を主軸変換と平方完成で標準形に直し、Sylvester の慣性法則で分類して、二次曲線・二次曲面の分類に応用する。
- 正定値性を首座小行列式・Cholesky 分解で判定し、半正定値平方根・極分解・二つの形式の同時対角化を構成する。

## 自己随伴作用素とユニタリ作用素

[固有値と対角化](/math/linear-algebra/eigenvalues/) では、対角化に用いる基底変換に何の制限も付けなかった。
そのため基底変換は長さや角度を自由に変えてしまう。
内積が与えられた空間では、正規直交基底どうしの変換だけで対角化できるかが自然な問題になる。
この章の主結果は、その可否が $T$ と随伴 $\adj{T}$ の可換性という代数的条件だけで決まるということである。

この章を通じて、特に断らない限り $K=\R$ または $\C$ とし、$V$ を $K$ 上の $n$ 次元内積空間（$n\ge1$）、$I=\id_V$ とする。
内積は第 1 引数について線形である（:ref[def-linalg-inner-product]）。
$T\in\End(V)$ の随伴 $\adj{T}$ は :ref[thm-linalg-adjoint] により一意に存在し、すべての $x,y$ について $\inner{Tx}{y}=\inner{x}{\adj{T}y}$ を満たす。
行列 $A\in M_n(K)$ は、標準内積をもつ $K^n$ 上の作用素 $x\mapsto Ax$ と同一視する。
その随伴は共役転置 $\adj{A}$（実行列では転置 $\transp{A}$）で表される（:ref[prop-linalg-adjoint-matrix]）。

:::definition[自己随伴・ユニタリ・正規]{#def-linalg-normal}
$T\in\End(V)$ とする。

1. $\adj{T}=T$ を満たすとき、$T$ を :term[自己随伴作用素]{en="self-adjoint operator" yomi="じこずいはんさようそ"} という。
2. $\adj{T}T=T\adj{T}=I$ を満たすとき、$T$ を :term[ユニタリ作用素]{en="unitary operator" yomi="ゆにたりさようそ"} という。$K=\R$ の場合は :term[直交変換]{en="orthogonal transformation" yomi="ちょっこうへんかん"} ともいう。
3. $\adj{T}T=T\adj{T}$ を満たすとき、$T$ を :term[正規作用素]{en="normal operator" yomi="せいきさようそ"} という。

行列 $A\in M_n(\C)$ については、$\adj{A}=A$ のとき :term[Hermite 行列]{en="Hermitian matrix" yomi="えるみーとぎょうれつ"}、$\adj{A}A=I$ のとき :term[ユニタリ行列]{en="unitary matrix" yomi="ゆにたりぎょうれつ"}、$\adj{A}A=A\adj{A}$ のとき :term[正規行列]{en="normal matrix" yomi="せいきぎょうれつ"} という。
実行列では、$\transp{A}=A$ のとき :term[対称行列]{en="symmetric matrix" yomi="たいしょうぎょうれつ"}、$\transp{A}A=I$ のとき :term[直交行列]{en="orthogonal matrix" yomi="ちょっこうぎょうれつ"}、$\transp{A}=-A$ のとき交代行列という。
:::

自己随伴作用素とユニタリ作用素はいずれも正規である。
正規直交基底 $\mathcal B$ に関する $T$ の表現行列を $A$ とすれば、$\adj{T}$ の表現行列は $\adj{A}$ である（:ref[prop-linalg-adjoint-matrix]）。
表現行列をとる対応は積を保つ（:ref[prop-linalg-composition-matrix]）ので、$T$ が自己随伴・ユニタリ・正規であることは、$A$ がそれぞれ Hermite（実対称）・ユニタリ（直交）・正規行列であることと同値である。
正規直交でない基底では、この対応は一般に成り立たない。
例えば自己随伴な作用素でも、斜交座標での表現行列は対称とは限らない。

### 二次式による判定

作用素そのものの代わりに、一変数の量 $v\mapsto\inner{Tv}{v}$ だけを調べて済む場合がある。
複素空間ではこの量が作用素を完全に決めるが、実空間では決めない。

:::proposition[内積 $\inner{Tv}{v}$ による判定]{#prop-linalg-quadratic-zero}
$T\in\End(V)$ とする。

1. $K=\C$ とする。すべての $v\in V$ で $\inner{Tv}{v}=0$ ならば $T=0$ である。
2. $K=\C$ とする。$T$ が自己随伴であるための必要十分条件は、すべての $v\in V$ で $\inner{Tv}{v}\in\R$ となることである。
3. $K=\R,\C$ のいずれでも、$T$ が自己随伴で、すべての $v$ で $\inner{Tv}{v}=0$ ならば $T=0$ である。
:::

:::proof
$Q(v)=\inner{Tv}{v}$ とする。$Q(u+i^kw)$ を展開すると、二つの対角項 $Q(u)+Q(w)$ と交差項 $i^{-k}\inner{Tu}{w}+i^k\inner{Tw}{u}$ を得る。$i^k$ を掛けて四つの $k$ で足せば $\sum i^k=\sum i^{2k}=0$ より
$$
\inner{Tu}{w}=\frac14\sum_{k=0}^3 i^k Q(u+i^kw). \label{eq-linalg-operator-polarization}
$$
したがって $Q=0$ なら $w=Tu$ として $Tu=0$ であり、1 が従う。$Q(v)-\overline{Q(v)}=\inner{(T-\adj T)v}{v}$ だから、1 を差に適用すれば $Q$ が実数値であることと $T=\adj T$ の同値性も従う。

3 の実の場合は自己随伴性と対称性から $\inner{Tu}{w}=\inner{Tw}{u}$ であり、展開して $Q(u+w)-Q(u-w)=4\inner{Tu}{w}$ となる。再び $w=Tu$ とすればよい。複素の場合は 1 に含まれる。
:::

:::example[実空間では二次式が作用素を決めない]{#exm-linalg-rotation-quadratic}
$\R^2$ の $90^\circ$ 回転

$$
R=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\qquad \inner{Rv}{v}=-v_2v_1+v_1v_2=0
$$

は零ではないが、すべての $v$ で $\inner{Rv}{v}=0$ を満たす。
したがって 1 は $K=\R$ では成り立たない。
$R$ は自己随伴でない（$\transp{R}=-R$）ので 3 とは矛盾しない。
また $K=\R$ では $\inner{Tv}{v}$ は常に実数だから、2 の条件は何も制限しない。
:::

:::proposition[自己随伴作用素とユニタリ作用素の固有値]{#prop-linalg-selfadjoint-eigen}
1. $T$ が自己随伴ならば、$T$ の固有値はすべて実数である。さらに $\Phi_T$ は $\C$ 上の根がすべて実数であり、$\R$ 上で一次式の積に分解する。
2. $T$ がユニタリならば、$T$ の固有値 $\lambda$ はすべて $\abs\lambda=1$ を満たす。
:::

:::proof
1. $Tv=\lambda v$、$v\ne0$ とすると

$$
\lambda\norm{v}^2=\inner{\lambda v}{v}=\inner{Tv}{v}=\inner{v}{Tv}=\inner{v}{\lambda v}=\bar\lambda\norm{v}^2
$$

であり、$\norm{v}\ne0$ より $\lambda=\bar\lambda$ である。
後半を示すため、正規直交基底をとって $T$ の表現行列を $A$ とする。
$A$ は Hermite 行列（$K=\R$ なら実対称行列）であり、複素行列とみれば $\C^n$ 上の標準内積について自己随伴である。
$\Phi_T=\Phi_A$ の複素根 $\lambda$ は $\C^n$ 上の $A$ の固有値なので、前半により実数である。
:ref[thm-linalg-triangularization] により $\Phi_A$ は $\C$ 上で $\prod_{j=1}^n(t-\lambda_j)$ と分解し、各 $\lambda_j$ は実数だから、これは $\R$ 上の分解でもある。

2. $Tv=\lambda v$、$v\ne0$ ならば $\norm{v}^2=\inner{\adj{T}Tv}{v}=\inner{Tv}{Tv}=\abs\lambda^2\norm{v}^2$ なので $\abs\lambda=1$ である。
:::

したがって実対称行列の固有値は、複素数上で求めてもすべて実数である。

:::example[鏡映と回転]
$\diag(1,-1)$ は自己随伴であり、直交行列でもある。
これに対し $R=\begin{pmatrix}0&-1\\1&0\end{pmatrix}$ は $\transp{R}R=I$ を満たすが、$\transp{R}=-R$ である。
複素数まで許せば固有値は $i,-i$ であり、どちらの絶対値も $1$ である。
また $\diag(2,-3)$ は自己随伴だが、固有値の絶対値が $1$ でないのでユニタリではない。
:::

### ユニタリ作用素と直交行列

:::proposition[長さを保つ変換]{#prop-linalg-unitary-equiv}
$U\in\End(V)$ について、次は同値である。

1. $\adj{U}U=I$。
2. $U\adj{U}=I$。
3. すべての $x,y\in V$ について $\inner{Ux}{Uy}=\inner{x}{y}$。
4. すべての $x\in V$ について $\norm{Ux}=\norm{x}$。
5. $U$ は任意の正規直交基底を正規直交基底に移す。
6. $U$ はある正規直交基底を正規直交基底に移す。

特に、行列 $A\in M_n(K)$ がユニタリ（直交）行列であることは、$A$ の列が $K^n$ の正規直交基底をなすことと同値であり、行が正規直交基底をなすこととも同値である。
:::

:::proof
$\adj U U=I$ なら $U$ は単射であり、有限次元の自己写像だから可逆で $\adj U=U^{-1}$ である。$U\adj U=I$ からも $\adj U$ に同じ議論を用いれば逆を得る。さらに $\inner{Ux}{Uy}=\inner{\adj U Ux}{y}$ であるから、1 は 3 を与える。逆に 3 が成り立てば $y=(\adj U U-I)x$ と取ってこの差のノルムが零になる。

3 に $y=x$ を入れれば 4 を得る。逆に線形写像 $U$ が長さを保てば、:ref[prop-linalg-polarization] の右辺に現れる各和・差の長さを保つので内積も保つ。3 ならどの正規直交基底の像も正規直交する $n$ 本であり、一次独立だから基底である。よって 5、したがって 6 を得る。6 の基底で $x=\sum a_j e_j,y=\sum b_j e_j$ と書けば、像の内積も元の内積も $\sum_j a_j\bar b_j$ に等しいため 3 に戻る。

行列では $(\adj A A)_{ij}=\inner{a_j}{a_i}$ が第 $j,i$ 列の内積である。したがって 1 は列の正規直交性、2 は同様に行の正規直交性と同値である。
:::

4 から、ユニタリ作用素は距離 $\norm{x-y}$ も保つ。
逆に、原点を固定して距離を保つ $\R^n$ の写像は自動的に線形で直交行列で表される（演習問題を参照）。

:::proposition[ユニタリ群と直交群]{#prop-linalg-unitary-group}
$n$ 次ユニタリ行列全体を $\mathrm{U}(n)$、$n$ 次実直交行列全体を $\mathrm{O}(n)$ と書く。

1. $\mathrm{U}(n)$ と $\mathrm{O}(n)$ は行列の積について群をなす。すなわち単位行列を含み、積と逆行列について閉じている。
2. $A\in\mathrm{U}(n)$ ならば $\abs{\det A}=1$ であり、$A\in\mathrm{O}(n)$ ならば $\det A=\pm1$ である。
3. $\SU(n)=\{A\in\mathrm{U}(n)\mid\det A=1\}$ と $\SO(n)=\{A\in\mathrm{O}(n)\mid\det A=1\}$ も群をなす。
4. $V$ の二つの正規直交基底の間の基底変換行列はユニタリ（$K=\R$ では直交）行列である。
5. $U$ がユニタリ行列ならば、$A$ が Hermite・ユニタリ・正規であることと、$\adj{U}AU$ がそれぞれ同じ性質をもつことは同値である。
:::

:::proof
$\adj I I=I$、$\adj{(AB)}AB=\adj B(\adj A A)B=I$、また $A^{-1}=\adj A$ も $A\adj A=I$ を満たす。よって単位元・積・逆元を含む。実行列どうしではすべての演算が実数内で行える。Leibniz の公式と :ref[prop-linalg-det-transpose] から $\det\adj A=\overline{\det A}$ なので、乗法性により $1=\det(\adj A A)=|\det A|^2$ である。実数なら値は $\pm1$ である。乗法性から、行列式 $1$ の元の積と逆元も行列式 $1$ となる。

二つの正規直交基底の間の変換は :ref[cor-linalg-coordinate-isometry] より列が正規直交し、ユニタリである。最後に $B=\adj U AU$ の随伴は $\adj B=\adj U\adj A U$ である。写像 $X\mapsto\adj U XU$ は $U\adj U=I$ により積と単位元を保ち、逆写像 $X\mapsto UX\adj U$ をもつ。したがって三種類の定義等式がすべて両方向に保たれる。
:::

$\mathrm{U}(n)$ を :term[ユニタリ群]{en="unitary group" yomi="ゆにたりぐん"}、$\mathrm{O}(n)$ を :term[直交群]{en="orthogonal group" yomi="ちょっこうぐん"} という。
また $\SU(n)$ を :term[特殊ユニタリ群]{en="special unitary group" yomi="とくしゅゆにたりぐん"}、$\SO(n)$ を :term[特殊直交群]{en="special orthogonal group" yomi="とくしゅちょっこうぐん"} という。
ユニタリ行列 $U$ によって $B=\adj{U}AU$ となるとき、$A$ と $B$ は :term[ユニタリ同値]{en="unitarily equivalent" yomi="ゆにたりどうち"}（実の場合は直交同値）であるという。
ユニタリ同値は相似の特別な場合であり、同じ作用素を二つの正規直交基底で表した行列どうしの関係である。

:::example[2 次の直交行列]{#exm-linalg-o2}
$A\in\mathrm{O}(2)$ の第 1 列は単位ベクトルなので $(\cos\theta,\sin\theta)$ と書ける。
第 2 列はこれに直交する単位ベクトルなので $\pm(-\sin\theta,\cos\theta)$ である。
したがって

$$
A=R_\theta=\begin{pmatrix}\cos\theta&-\sin\theta\\\sin\theta&\cos\theta\end{pmatrix}
\quad\text{または}\quad
A=\begin{pmatrix}\cos\theta&\sin\theta\\\sin\theta&-\cos\theta\end{pmatrix}
$$

である。
前者は角 $\theta$ の回転で $\det=1$、すなわち $\SO(2)=\{R_\theta\}$ である。
後者は $\det=-1$ であり、$(\cos\frac\theta2,\sin\frac\theta2)$ を固定し、それに直交する $(-\sin\frac\theta2,\cos\frac\theta2)$ を $-1$ 倍する。
実際、加法定理から

$$
\begin{pmatrix}\cos\theta&\sin\theta\\\sin\theta&-\cos\theta\end{pmatrix}\begin{pmatrix}\cos\frac\theta2\\\sin\frac\theta2\end{pmatrix}
=\begin{pmatrix}\cos(\theta-\frac\theta2)\\\sin(\theta-\frac\theta2)\end{pmatrix}
$$

であり、もう一方も同様に確かめられる。
したがって後者は角 $\theta/2$ の直線に関する鏡映である。
:::

### 正規作用素の基本性質

:::proposition[正規作用素の特徴づけと固有空間]{#prop-linalg-normal-basic}
$T\in\End(V)$ とする。

1. $T$ が正規であるための必要十分条件は、すべての $x\in V$ で $\norm{Tx}=\norm{\adj{T}x}$ となることである。
2. $T$ が正規ならば、任意の $\lambda\in K$ について $T-\lambda I$ も正規であり、$E_\lambda(T)=E_{\bar\lambda}(\adj{T})$ である。特に $\Ker T=\Ker\adj{T}$ である。
3. $T$ が正規ならば、相異なる固有値に属する固有ベクトルは直交する。
:::

:::proof
1. 任意の $x$ について

$$
\norm{Tx}^2-\norm{\adj{T}x}^2=\inner{\adj{T}Tx}{x}-\inner{T\adj{T}x}{x}=\inner{(\adj{T}T-T\adj{T})x}{x}
$$

である。
$\adj{T}T-T\adj{T}$ は自己随伴である（:ref[prop-linalg-adjoint-properties]）。
したがって :ref[prop-linalg-quadratic-zero] の 3 により、右辺がすべての $x$ で零であることは $\adj{T}T=T\adj{T}$ と同値である。
この議論は $K=\R$ でも有効である。

2. :ref[prop-linalg-adjoint-properties] より $\adj{(T-\lambda I)}=\adj{T}-\bar\lambda I$ である。
展開すれば $(T-\lambda I)(\adj{T}-\bar\lambda I)$ と $(\adj{T}-\bar\lambda I)(T-\lambda I)$ はどちらも $T\adj{T}-\bar\lambda T-\lambda\adj{T}+\abs\lambda^2I$ に等しいので、$T-\lambda I$ は正規である。
1 を $T-\lambda I$ に適用すると $\norm{(T-\lambda I)v}=\norm{(\adj{T}-\bar\lambda I)v}$ となるから、$Tv=\lambda v$ と $\adj{T}v=\bar\lambda v$ は同値である。

3. $Tv=\lambda v$、$Tw=\mu w$、$\lambda\ne\mu$ とする。
2 より $\adj{T}w=\bar\mu w$ なので

$$
\lambda\inner{v}{w}=\inner{Tv}{w}=\inner{v}{\adj{T}w}=\inner{v}{\bar\mu w}=\mu\inner{v}{w}
$$

である。
$\lambda\ne\mu$ より $\inner{v}{w}=0$ である。
:::

:::example[正規性の有無と核]
$c\in\C$ に対するスカラー作用素 $cI$ は常に正規であり、$c=2i$ なら自己随伴でもユニタリでもない。
$\begin{pmatrix}1&-1\\1&1\end{pmatrix}=\sqrt2R_{\pi/4}$ も正規だが、自己随伴でもユニタリでもない。
一方、$N=\begin{pmatrix}0&1\\0&0\end{pmatrix}$ では $\adj{N}N=\diag(0,1)$、$N\adj{N}=\diag(1,0)$ なので正規でない。
実際 $\Ker N=\spn\{e_1\}$ だが $\Ker\adj{N}=\spn\{e_2\}$ であり、2 の結論が成り立たない。
:::

### 随伴で保たれる部分空間

作用素 $T$ の不変部分空間 $W$（:ref[def-linalg-invariant-subspace]）があっても、直交補空間 $W^\perp$ が不変とは限らない。
:ref[prop-linalg-invariant-adjoint] により

$$
T(W)\subset W \quad\Longleftrightarrow\quad \adj{T}(W^\perp)\subset W^\perp
$$

であるから、$W^\perp$ の $T$ 不変性は $W$ の $\adj{T}$ 不変性と同値である。
$W$ と $W^\perp$ がともに $T$ 不変なら、両者に合わせた正規直交基底で $T$ はブロック対角になる（:ref[prop-linalg-invariant-block]）。
正規作用素ではこれが自動的に成り立つ。

:::proposition[正規作用素の不変部分空間]{#prop-linalg-normal-invariant}
$K=\R$ または $\C$ とし、$T\in\End(V)$ を正規作用素、$W$ を $T$ 不変部分空間とする。このとき次が成り立つ。

1. $W^\perp$ も $T$ 不変である。
2. $W$ は $\adj{T}$ 不変であり、$\adj{(T|_W)}=\adj{T}|_W$ である。
3. $T|_W$ と $T|_{W^\perp}$ はそれぞれ $W$ と $W^\perp$ 上の正規作用素である。
:::

:::proof
$W$ の正規直交基底を :ref[cor-linalg-orthonormal-basis] により全体の正規直交基底に延長する。追加した列は :ref[prop-linalg-complement-properties] により $W^\perp$ の基底をなす。不変性と :ref[prop-linalg-invariant-block] から、この基底で $T$ の行列は $M=\begin{pmatrix}A&B\\0&D\end{pmatrix}$ である。随伴の行列は $\adj M$ だから、$\adj M M=M\adj M$ の左上ブロックは
$$
\adj A A=A\adj A+B\adj B.
$$
:ref[lem-linalg-trace-commute] を使ってトレースを取ると $\sum_{ij}|b_{ij}|^2=\tr(B\adj B)=0$ となり $B=0$ である。よって $M$ と $\adj M$ は同じ二空間を保つ。$W$ 上で $\adj T$ の行列は $\adj A$ なので制限の随伴に等しい。さらに対角ブロックの比較から $\adj A A=A\adj A,\ \adj D D=D\adj D$ を得る。零次元の成分はそのまま省けばよい。
:::

この証明は複素数の固有値を使わないので、実内積空間でもそのまま成り立つ。
これが後の実正規行列の標準形（:ref[prop-linalg-real-normal-form]）の証明を支える。
また、$W$ と $W^\perp$ がともに $T$ 不変であることは、直交射影 $P_W$ が $T$ と可換であることと同値である。
実際、両者が不変なら $x=w+w'$（$w\in W,\ w'\in W^\perp$）について $P_WTx=Tw=TP_Wx$ である。
逆に $P_WT=TP_W$ なら、$T$ は $P_W$ の像 $W$ と核 $W^\perp$ をそれぞれ保つ。

:::example[同じ部分空間を保つ二つの変換]
$W=\spn\{e_1\}$ とし、$A=\begin{pmatrix}1&2\\0&3\end{pmatrix}$ と $D=\diag(1,3)$ を比べる。
ともに $W$ を保つが、$A$ は $e_2$ を $(2,3)$ に送るので $W^\perp=\spn\{e_2\}$ を保たない。
:ref[prop-linalg-normal-invariant] より $A$ は正規でありえず、実際 $\adj{A}A$ と $A\adj{A}$ の $(1,1)$ 成分は $1$ と $5$ で異なる。
射影との交換子は

$$
P_W A-AP_W=\begin{pmatrix}0&2\\0&0\end{pmatrix}
$$

であり、証明中のブロック $B$ に当たる右上成分が残っている。
$D$ の場合はこの差が零なので $W^\perp$ も独立して調べられる。
:::

## 正規作用素のスペクトル定理

:ref[thm-linalg-triangularization] では、特性多項式が一次式の積に分解すれば上三角化できることを示した。
内積空間では、基底を正規直交基底に制限してもなお上三角化できる。
これがスペクトル定理への出発点である。

:::theorem[Schur の三角化]{#thm-linalg-schur}
$T\in\End(V)$ の特性多項式が $K$ 上で $\Phi_T=\prod_{i=1}^n(t-\lambda_i)$ と一次式の積に分解するとする（$K=\C$ なら常に成り立つ。:ref[thm-linalg-triangularization]）。
このとき $V$ の正規直交基底で、$T$ の表現行列が対角成分 $\lambda_1,\ldots,\lambda_n$ をこの順にもつ上三角行列になるものが存在する。
行列で言えば、$\Phi_A$ が $K$ 上で分解する $A\in M_n(K)$ に対し、ユニタリ行列（$K=\R$ では直交行列）$U$ で $\adj{U}AU$ が上三角になるものが存在し、対角成分の並び順は任意に指定できる。
:::

:::proof
指定順の最初の固有値 $\lambda_1$ に対し、単位固有ベクトル $q_1$ を選ぶ。$W=Kq_1$ は不変である。$q_1$ を正規直交基底に延長すると、:ref[prop-linalg-invariant-block] より行列は $\begin{pmatrix}\lambda_1&b\\0&C\end{pmatrix}$ になる。$C$ は $W^\perp$ 上の作用素 $S=P_{W^\perp}T|_{W^\perp}$ の行列である。実際、残りの基底ベクトルの像から $q_1$ 成分を除けば $C$ の各列が残る。この $S$ を :term[圧縮]{en="compression" yomi="あっしゅく"} という。

:ref[cor-linalg-invariant-charpoly] より $\Phi_T=(t-\lambda_1)\Phi_C$ であり、:ref[lem-linalg-polynomial-division] により約分すると $\Phi_C=\prod_{j=2}^n(t-\lambda_j)$ である。$n=1$ ならここで終了する。$n>1$ なら次元に関する帰納法で $W^\perp$ に $S$ の上三角形を与える正規直交基底 $q_2,\ldots,q_n$ を、指定順で選べる。$Tq_j=\inner{Tq_j}{q_1}q_1+Sq_j$ なので、全基底でも第 $j$ 列は第 $j$ 行までに収まり、対角成分は $\lambda_j$ である。行列の場合、この基底を列に並べた $U$ はユニタリで $\adj U A U$ が求めた上三角行列となる。
:::

:::example[Schur 形が最良である場合]
$A=\begin{pmatrix}1&1\\-1&3\end{pmatrix}$ では $\Phi_A=(t-1)(t-3)+1=(t-2)^2$ である。
固有値 $2$ の固有ベクトルは $(1,1)$ の倍数だけなので、$A$ は対角化可能でない（:ref[thm-linalg-diagonalization-criterion]）。
$q_1=(1,1)/\sqrt2$ を $q_2=(1,-1)/\sqrt2$ で正規直交基底に延長し、$U=(q_1\ q_2)$ とすると

$$
\adj{U}AU=\begin{pmatrix}2&-2\\0&2\end{pmatrix}
$$

となる。
右上の $-2$ は正規直交基底をどう選んでも消せない。
実際、ユニタリ同値は $\tr(\adj{A}A)=\sum_{i,j}\abs{a_{ij}}^2$ を保つ（:ref[lem-linalg-trace-commute]）ので、上三角形の右上成分の絶対値の二乗は $\tr(\adj{A}A)-\abs{2}^2-\abs{2}^2=12-8=4$ に固定されている。
:::

:::example[圧縮が必要な理由]
$T=\begin{pmatrix}1&2\\0&3\end{pmatrix}$ では $e_1$ は固有ベクトルだが、$Te_2=2e_1+3e_2$ は $e_1^\perp$ に入らない。
したがって $T$ を $e_1^\perp$ に制限することはできず、証明では $e_1$ 成分を捨てる圧縮 $e_2\mapsto3e_2$ を用いた。
この例では $\adj{T}T$ と $T\adj{T}$ の左上成分が $1,5$ なので $T$ は正規でない。
:::

$K=\R$ では仮定は本質的である。
例えば回転 $R_\theta$（$0<\theta<\pi$）は実固有値をもたないので、どの実正規直交基底でも上三角にならない。
実行列一般については、次の準三角形まで簡約できる。

:::theorem[実 Schur 分解]{#thm-linalg-real-schur}
任意の $A\in M_n(\R)$ に対し、直交行列 $Q$ で $\transp{Q}AQ$ がブロック上三角行列になり、各対角ブロックが $1\times1$ または $2\times2$ であるものが存在する。
$1\times1$ ブロックは $A$ の実固有値であり、$2\times2$ ブロックの固有値は $A$ の互いに共役な非実固有値の対である。
:::

:::proof
次元が零なら空の列でよい。実固有値があればその実固有直線を $W$ とする。なければ複素固有値 $a+ib$（$b\ne0$）と固有ベクトル $u+iv$（$u,v$ は実）を選ぶ。実部と虚部をとると $Au=au-bv,\ Av=bu+av$ なので $W=\spn\{u,v\}$ は不変である。もし一次元ならその上での $A$ は実スカラー倍となり、$u+iv$ の固有値が非実であることに反する。よってこの場合は二次元であり、基底 $(u,v)$ での制限の特性多項式は $(t-a)^2+b^2$ である。

$W$ の正規直交基底を全体に延長する。:ref[prop-linalg-invariant-block] により行列は $\begin{pmatrix}A_1&B\\0&C\end{pmatrix}$ となり、$A_1$ は上で調べた制限の行列である。$C$ に次元についての帰納法を適用し、残りの基底の中だけで直交変換を施せばよい。:ref[cor-linalg-invariant-charpoly] によりこの操作は重複度も含めて固有値をすべて振り分ける。
:::

### スペクトル定理

正規性はユニタリ同値で保たれる（:ref[prop-linalg-unitary-group] の 5）。
したがって Schur の三角化で得られる上三角行列も正規である。
ところが正規な上三角行列は対角行列しかない。

:::lemma[正規な上三角行列]{#lem-linalg-normal-triangular}
上三角行列 $B\in M_n(\C)$ が正規ならば、$B$ は対角行列である。
:::

:::proof
$n$ に関する帰納法で示す。
$n=1$ は明らかである。
$B$ の第 1 列は $(b_{11},0,\ldots,0)$ なので、$\adj{B}B$ と $B\adj{B}$ の $(1,1)$ 成分はそれぞれ

$$
(\adj{B}B)_{11}=\sum_{i=1}^n\abs{b_{i1}}^2=\abs{b_{11}}^2,\qquad
(B\adj{B})_{11}=\sum_{j=1}^n\abs{b_{1j}}^2
$$

である。
両者が等しいので $b_{1j}=0$（$j\ge2$）となり、$B=\diag(b_{11},B')$ とブロック対角に書ける（$B'\in M_{n-1}(\C)$ は上三角）。
ブロックごとに積をとると

$$
\adj{B}B=\diag(\abs{b_{11}}^2,\adj{B'}B'),\qquad B\adj{B}=\diag(\abs{b_{11}}^2,B'\adj{B'})
$$

であり、$\adj{B}B=B\adj{B}$ から $\adj{B'}B'=B'\adj{B'}$ を得る。
よって $B'$ は正規な上三角行列であり、帰納法の仮定により対角行列である。
:::

:::theorem[複素スペクトル定理]{#thm-linalg-spectral}
$K=\C$ とする。
$T\in\End(V)$ が正規であるための必要十分条件は、$V$ が $T$ の固有ベクトルからなる正規直交基底をもつことである。
このとき $T$ の相異なる固有値を $\lambda_1,\ldots,\lambda_r$、固有空間 $E_{\lambda_j}(T)$ への直交射影を $P_j$ とすれば、

$$
V=E_{\lambda_1}(T)\oplus\cdots\oplus E_{\lambda_r}(T)\quad(\text{互いに直交}),\qquad
I=\sum_{j=1}^rP_j,\qquad
T=\sum_{j=1}^r\lambda_jP_j,\qquad
P_iP_j=\delta_{ij}P_j \label{eq-linalg-spectral-decomposition}
$$

が成り立つ。
行列で言えば、$A\in M_n(\C)$ が正規であることと、ユニタリ行列 $U$ と対角行列 $D$ で $A=UD\adj{U}$ と書けることは同値である。
:::

:::proof
正規なら :ref[thm-linalg-schur] の正規直交基底で上三角行列 $B$ を得る。:ref[prop-linalg-adjoint-matrix] により $\adj B B=B\adj B$ だから :ref[lem-linalg-normal-triangular] で $B$ は対角行列である。逆に正規直交固有基底があれば $T,\adj T$ を表す対角行列 $D,\bar D$ は可換なので $T$ は正規である。

この基底を固有値ごとに分ける。$x=\sum_e c_e e$ に $T-\lambda_j I$ を作用させると $e$ の係数は $(\lambda_e-\lambda_j)c_e$ だから、値が零になるのは $\lambda_e\ne\lambda_j$ の係数がすべて零のときである。したがって各組が対応する固有空間の基底である。$P_j$ はその組の座標だけを残すので $\sum_jP_j=I$、$P_iP_j=\delta_{ij}P_j$ である。各組で $T$ が $\lambda_j$ 倍になることから $T=\sum_j\lambda_jP_j$ を得る。基底の列を $U$ に並べれば $\adj U U=I$ かつ $AU=UD$ であり、行列表示の主張も従う。
:::

:::remark
:ref[prop-linalg-normal-invariant] を用いる別証明もある。
$K=\C$ なので $T$ は単位固有ベクトル $v$ をもち、$\spn\{v\}$ は $T$ 不変である。
よって $v^\perp$ も $T$ 不変で、$T|_{v^\perp}$ は正規である。
次元に関する帰納法で $v^\perp$ に固有ベクトルからなる正規直交基底をとり、$v$ を加えればよい。
この証明は Schur の三角化を経由しないが、本質的には同じ「固有ベクトルを一本取って直交補空間へ降りる」操作である。
:::

式 :eqref[eq-linalg-spectral-decomposition] の右辺を $T$ の :term[スペクトル分解]{en="spectral decomposition" yomi="すぺくとるぶんかい"}、$P_j$ を :term[スペクトル射影]{en="spectral projection" yomi="すぺくとるしゃえい"} という。
固有値が重複する場合、固有ベクトルからなる正規直交基底は一意でないが、スペクトル射影は $T$ だけで決まる。

:::proposition[スペクトル分解の一意性と多項式表示]{#prop-linalg-spectral-projections}
$K=\C$、$T$ を正規作用素とし、記号を :ref[thm-linalg-spectral] のとおりとする。

1. $p_j(t)=\prod_{k\ne j}\dfrac{t-\lambda_k}{\lambda_j-\lambda_k}$ とおけば $P_j=p_j(T)$ である。
2. 相異なる $\mu_1,\ldots,\mu_s\in\C$ と、零でない直交射影 $Q_1,\ldots,Q_s$ が $Q_iQ_k=0$（$i\ne k$）、$\sum_iQ_i=I$、$T=\sum_i\mu_iQ_i$ を満たすならば、$s=r$ であり、番号を付け替えると $\mu_i=\lambda_i$、$Q_i=P_i$ である。
3. $\adj{T}=\sum_j\bar\lambda_jP_j$ であり、特に $\adj{T}$ は $T$ の多項式として表される。
:::

:::proof
:ref[lem-linalg-polynomial-substitution] より $p_j(T)$ は $E_{\lambda_k}$ 上で $p_j(\lambda_k)=\delta_{jk}$ 倍である。固有空間の直和が全体を覆うので $p_j(T)=P_j$ となる。

別の表示については $TQ_i=\mu_iQ_i$ より $\Img Q_i\subset E_{\mu_i}$ で、$Q_i\ne0$ だから各 $\mu_i$ は固有値である。逆に $Tx=\lambda x$ なら $Q_iT=\mu_iQ_i$ から $(\lambda-\mu_i)Q_ix=0$ となる。$x=\sum_iQ_ix$ より、$\lambda$ が各 $\mu_i$ と異なるなら $x=0$ であり、$\lambda=\mu_j$ なら $x=Q_jx$ である。よって固有値集合は一致し、$\Img Q_j=E_{\mu_j}$ である。その空間への直交射影は一意なので $Q_j$ も決まる。

各 $P_j$ は自己随伴だから、分解の随伴を取ると $\adj T=\sum_j\bar\lambda_jP_j$ となる。$P_j=p_j(T)$ を代入すれば多項式表示も得られる。
:::

:::corollary[正規作用素の種類と固有値]{#cor-linalg-normal-classes}
$K=\C$、$T$ を正規作用素とする。

1. $T$ が自己随伴 ⇔ $T$ の固有値はすべて実数。
2. $T$ がユニタリ ⇔ $T$ の固有値はすべて絶対値 $1$。
3. $\adj{T}=-T$ ⇔ $T$ の固有値はすべて純虚数または $0$。
4. $T$ が直交射影 ⇔ $T$ の固有値はすべて $0$ または $1$。
:::

:::proof
スペクトル分解 $T=\sum_j\lambda_jP_j$ と :ref[prop-linalg-spectral-projections] の 3 を用いる。
$E_{\lambda_j}(T)$ の元 $x\ne0$ に作用させると $Tx=\lambda_jx$、$\adj{T}x=\bar\lambda_jx$、$\adj{T}Tx=\abs{\lambda_j}^2x$、$T^2x=\lambda_j^2x$ である。
したがって $T=\adj{T}$ なら $\lambda_j=\bar\lambda_j$、$\adj{T}T=I$ なら $\abs{\lambda_j}=1$、$\adj{T}=-T$ なら $\bar\lambda_j=-\lambda_j$ である。
$T$ が直交射影なら $T^2=T$ より $\lambda_j^2=\lambda_j$ である。
逆に固有値がそれぞれの条件を満たせば、$\adj{T}=\sum_j\bar\lambda_jP_j$、$\adj{T}T=\sum_j\abs{\lambda_j}^2P_j$、$T^2=\sum_j\lambda_j^2P_j$ を $T$ や $I=\sum_jP_j$ と比べて結論を得る。
4 では $T^2=T=\adj{T}$ となるので :ref[prop-linalg-orthogonal-projection-char] を用いる。
:::

正規性の仮定は外せない。
例えば $\begin{pmatrix}1&2\\0&3\end{pmatrix}$ の固有値は実数だが、自己随伴ではない。
2 から、ユニタリ行列 $U$ はユニタリ行列 $W$ により $\adj{W}UW=\diag(e^{i\theta_1},\ldots,e^{i\theta_n})$ と対角化される。

:::example[射影を固有値から計算する]
$H=\begin{pmatrix}0&-i\\i&0\end{pmatrix}$ は Hermite 行列である。
$H^2=I$ でトレースが零なので、:ref[cor-linalg-trace-det-eigenvalues] より固有値は $1,-1$ である。
:ref[prop-linalg-spectral-projections] の 1 から

$$
P_{1}=\frac{H-(-1)I}{1-(-1)}=\frac{I+H}2,\qquad P_{-1}=\frac{I-H}2
$$

であり、$P_1+P_{-1}=I$、$P_1P_{-1}=(I-H^2)/4=0$、$P_1-P_{-1}=H$ を直接確かめられる。
:::

### 可換性と共通固有基底

:::theorem[同時対角化]{#thm-linalg-simultaneous}
$K=\C$ とし、$T_1,\ldots,T_k\in\End(V)$ を正規作用素とする。
これらすべてを同時に対角化する正規直交基底が存在するための必要十分条件は、どの二つも可換であることである。
:::

:::proof
共通の正規直交基底で対角行列になれば、対角行列どうしは可換なので $T_iT_j=T_jT_i$ である。

逆を $\dim V$ に関する帰納法で示す。
$\dim V=1$ なら任意の単位ベクトルでよい。
すべての $T_i$ がスカラー倍作用素なら、任意の正規直交基底でよい。
そうでないとき、例えば $T_1$ がスカラーでないとすると、:ref[thm-linalg-spectral] により $T_1$ は二つ以上の相異なる固有値をもち、各固有空間 $E=E_\lambda(T_1)$ は $\dim E<\dim V$ を満たす。
$x\in E$ なら $T_1T_jx=T_jT_1x=\lambda T_jx$ なので $E$ は各 $T_j$ で不変である。
:ref[prop-linalg-normal-invariant] より $T_j|_E$ は $E$ 上の正規作用素であり、これらは互いに可換である。
帰納法の仮定により $E$ にはすべての $T_j|_E$ を対角化する正規直交基底がある。
$V$ は $T_1$ の固有空間の直交直和なので、各固有空間の基底を合わせれば $V$ の正規直交基底になり、すべての $T_j$ を対角化する。
:::

:::example[重複固有値の中での基底選択]
$A=\diag(2,2,3)$ と、最初の二座標を入れ替え第 3 座標を固定する置換行列 $B$ は可換である。
$A$ の固有値 $2$ の平面で $(1,1,0)/\sqrt2,(1,-1,0)/\sqrt2$ を選ぶと、$B$ はそれぞれを $1,-1$ 倍する。
残りに $e_3$ を加えれば共通の固有ベクトルからなる正規直交基底である。
$A$ の固有ベクトルからなる正規直交基底を任意に選んでも（例えば $e_1,e_2,e_3$）、それが $B$ を対角化するとは限らない。
:::

### 固有基底を求めずに関数を作用させる

スペクトル射影が分かれば、固有値が重複していても、固有空間の中で個々の軸を選ぶ必要はない。
:eqref[eq-linalg-spectral-decomposition] と $P_iP_j=\delta_{ij}P_j$ から、帰納的に $T^m=\sum_j\lambda_j^mP_j$ であり、多項式 $p$ について

$$
p(T)=\sum_jp(\lambda_j)P_j
$$

となる（定数項では $I=\sum_jP_j$ を用いる）。
この右辺は多項式に限らず意味をもつ。

:::definition[正規作用素の関数]{#def-linalg-normal-function}
$K=\C$、$T=\sum_j\lambda_jP_j$ を正規作用素のスペクトル分解とする。
固有値の集合上の複素数値関数 $f$ に対し $f(T)=\sum_jf(\lambda_j)P_j$ と定める。
:::

$f(T)$ は多項式 $q$ で $q(\lambda_j)=f(\lambda_j)$（$j=1,\ldots,r$）となるもの、例えば $q=\sum_jf(\lambda_j)p_j$ によって $f(T)=q(T)$ と書ける（:ref[prop-linalg-spectral-projections] の 1）。
したがって $f(T)$ は $T$ と可換なすべての作用素と可換である。
また $P_iP_j=\delta_{ij}P_j$ から $(fg)(T)=f(T)g(T)$、$(f+g)(T)=f(T)+g(T)$ が、:ref[prop-linalg-spectral-projections] の 3 と同様に $\bar f(T)=\adj{f(T)}$ が成り立つ。
特に $f(T)$ は正規であり、その固有値は $f(\lambda_j)$ である。
指数関数などの一般の行列関数は [行列解析](/math/linear-algebra/matrix-analysis/) で扱う。

:::example[一つの固有値が重複する空間]
$\R^3$ で $u=(1,1,1)$ とし、$A=2I+u\transp{u}$ を考える。
$\transp{u}u=3$ より $Au=5u$ であり、$x\perp u$ なら $Ax=2x$ である。
したがって固有値は $5$（重複度 $1$）と $2$（重複度 $2$）で、射影は

$$
P_5=u\transp{u}/3,\qquad P_2=I-P_5
$$

である。
$f(t)=1/t$ を作用させて $A^{-1}=P_2/2+P_5/5$ を得る。
$b=(3,1,2)$ の $u$ 方向の成分は $(\transp{u}b/3)u=(2,2,2)$、直交成分は $(1,-1,0)$ なので、$Ax=b$ の解は

$$
x=\tfrac12(1,-1,0)+\tfrac15(2,2,2)
=(9/10,-1/10,2/5)
$$

である。
固有値 $2$ の平面上の正規直交基底には無数の選択肢があるが、$P_2$ や解 $x$ はどの選択にも依存しない。
:::

## 実対称行列の直交対角化

実係数では、回転のように実固有値をもたない作用素がある。
自己随伴性があれば、:ref[prop-linalg-selfadjoint-eigen] により特性多項式は実数の範囲で分解するので、実空間の中だけで直交する軸を選べる。

:::theorem[実スペクトル定理]{#thm-linalg-real-spectral}
$K=\R$ とする。
$T\in\End(V)$ が自己随伴であるための必要十分条件は、$V$ が $T$ の固有ベクトルからなる正規直交基底をもつことである。
このとき :ref[thm-linalg-spectral] と同じ形のスペクトル分解 $T=\sum_j\lambda_jP_j$（$\lambda_j\in\R$）が成り立つ。
行列で言えば、$A\in M_n(\R)$ が対称であることと、直交行列 $Q$ と実対角行列 $D$ で $A=QD\transp{Q}$ と書けることは同値である。
:::

:::proof
$T$ が自己随伴ならば、:ref[prop-linalg-selfadjoint-eigen] の 1 により $\Phi_T$ は $\R$ 上で一次式の積に分解する。
したがって :ref[thm-linalg-schur] が $K=\R$ でも適用でき、表現行列 $B$ が上三角になる正規直交基底がとれる。
$T$ は自己随伴なので $B$ は対称行列でもある（:ref[prop-linalg-adjoint-matrix]）。
上三角かつ対称な行列は対角行列であるから、この基底は固有ベクトルからなる。

逆に固有ベクトルからなる正規直交基底があれば、そこでの表現行列 $D$ は実対角行列で $\transp{D}=D$ なので、$T$ は自己随伴である。
スペクトル分解は :ref[thm-linalg-spectral] の証明の後半と同じ議論で得られる。
行列の主張は、固有ベクトルからなる正規直交基底を列に並べて $Q$ とすればよい。
逆に $A=QD\transp{Q}$ なら $\transp{A}=Q\transp{D}\transp{Q}=A$ である。
:::

$K=\C$ の自己随伴作用素は正規なので、:ref[thm-linalg-spectral] と :ref[cor-linalg-normal-classes] により、実固有値をもつ正規直交固有基底をもつ。
これに対し $K=\R$ では、自己随伴でない正規作用素（例えば回転）は一般に実の固有基底をもたない。
実数上では直交行列による対角化と自己随伴性が同値である。例えば $\begin{pmatrix}1&2\\0&3\end{pmatrix}$ は異なる二実固有値をもつので対角化可能だが、転置と一致しないため直交対角化はできない。

:::example[二次元で軸を求める]
$A=\begin{pmatrix}2&1\\1&2\end{pmatrix}$ の特性多項式は $(t-3)(t-1)$ である。
$q_1=(1,1)/\sqrt2,\ q_2=(1,-1)/\sqrt2$ とすれば、固有値はこの順に $3,1$ である。
新座標 $u=(x+y)/\sqrt2,\ v=(x-y)/\sqrt2$ を用いると $2x^2+2xy+2y^2=3u^2+v^2$ になる。
:::

:::example[重複する固有値をもつ 3 次の対称行列]{#exm-linalg-symmetric-3x3}
$A=\begin{pmatrix}2&1&1\\1&2&1\\1&1&2\end{pmatrix}=I+J$ とする。ここで $J$ は全成分が $1$ の行列である。
$J(1,1,1)=3(1,1,1)$ であり、$x_1+x_2+x_3=0$ なら $Jx=0$ である。
よって $A$ の固有値は $4$（固有空間 $\spn\{(1,1,1)\}$）と $1$（固有空間は平面 $x_1+x_2+x_3=0$）である。
固有値 $1$ の固有空間の基底 $(1,-1,0),(1,0,-1)$ に :ref[thm-linalg-gram-schmidt] を適用すると、

$$
q_2=\frac{(1,-1,0)}{\sqrt2},\qquad
(1,0,-1)-\frac12(1,-1,0)=\Bigl(\frac12,\frac12,-1\Bigr),\qquad
q_3=\frac{(1,1,-2)}{\sqrt6}
$$

を得る。
$q_1=(1,1,1)/\sqrt3$ と合わせて $Q=(q_1\ q_2\ q_3)$ とすれば、$\transp{Q}AQ=\diag(4,1,1)$ である。
異なる固有値の固有空間は自動的に直交するので（:ref[prop-linalg-normal-basic] の 3）、直交化が必要なのは同じ固有空間の内部だけである。
:::

### Rayleigh 商と固有値の極値

:::proposition[Rayleigh 商の値域]{#prop-linalg-rayleigh}
$T$ を自己随伴作用素とし、固有値を重複度込みで $\lambda_1\le\cdots\le\lambda_n$ と並べる。
$x\ne0$ に対し

$$
R_T(x)=\frac{\inner{Tx}{x}}{\norm{x}^2}
$$

を :term[Rayleigh 商]{en="Rayleigh quotient" yomi="れいりーしょう"} という。
このとき $\lambda_1\le R_T(x)\le\lambda_n$ であり、$R_T(x)=\lambda_n$ となるのは $x\in E_{\lambda_n}(T)$ のとき、$R_T(x)=\lambda_1$ となるのは $x\in E_{\lambda_1}(T)$ のときに限る。
特に $\lambda_n=\max_{x\ne0}R_T(x)$、$\lambda_1=\min_{x\ne0}R_T(x)$ である。
:::

:::proof
固有ベクトルからなる正規直交基底 $(e_j)$、$Te_j=\lambda_je_j$ をとり、$x=\sum_jc_je_j$ とする。
:ref[prop-linalg-orthonormal-expansion] より

$$
R_T(x)=\frac{\sum_j\lambda_j\abs{c_j}^2}{\sum_j\abs{c_j}^2}
$$

であり、これは固有値の重み $\abs{c_j}^2/\sum_k\abs{c_k}^2$ による平均なので、$\lambda_1$ 以上 $\lambda_n$ 以下である。
$\lambda_n-R_T(x)=\sum_j(\lambda_n-\lambda_j)\abs{c_j}^2/\norm{x}^2$ の各項は非負なので、これが零であることは $\lambda_j<\lambda_n$ となるすべての $j$ で $c_j=0$ であること、すなわち $x\in E_{\lambda_n}(T)$ と同値である。
最小値についても同様である。
:::

中間の固有値 $\lambda_k$ も、部分空間の上で Rayleigh 商の最大値・最小値をとる操作を組み合わせて特徴づけられる（Courant–Fischer の変分原理）。
これは Hermite 行列の固有値の摂動評価（Weyl の不等式）の基礎であり、[行列解析](/math/linear-algebra/matrix-analysis/) で証明する。
行列の摂動論における応用は :cite[horn-johnson] でも扱われる。

対ごとに可換な実対称行列も、一つの直交行列で同時に対角化できる。
証明は :ref[thm-linalg-simultaneous] と同じ帰納法であり、:ref[thm-linalg-spectral] の代わりに実スペクトル定理を、固有空間への制限の自己随伴性には :ref[prop-linalg-normal-invariant] を用いればよい。

### 実正規作用素の一次元・二次元分解

実正規行列が実数上で対角化できないとき、どこまで実の正規直交基底で簡単にできるだろうか。
複素固有値を共役な対としてまとめると、各対を実平面上の回転と伸縮として表せる。

:::proposition[実正規行列の直交標準形]{#prop-linalg-real-normal-form}
$A\in M_n(\R)$ が $\transp{A}A=A\transp{A}$ を満たすとする。
このとき直交行列 $Q$ により、$\transp{Q}AQ$ は実数からなる $1\times1$ ブロックと

$$
\begin{pmatrix}a&-b\\ b&a\end{pmatrix},\qquad a,b\in\R,\ b>0
$$

の形の $2\times2$ ブロックとのブロック対角行列になる。
$1\times1$ ブロックは $A$ の実固有値であり、$2\times2$ ブロックは $A$ の非実固有値の対 $a\pm ib$ に対応する。
逆に、この形の行列と直交同値な行列は正規である。
:::

:::proof
複素数上で実行列 $A$ を考える。実固有値に対しては行基本変形を実数上で行えるので、実の単位固有ベクトルを選べる。張る直線とその直交補は :ref[prop-linalg-normal-invariant] により不変で、直交補への制限も正規である。

非実固有値は共役対で現れるので $\lambda=a+ib,\ b>0$ を選ぶ。固有ベクトル $z=u+iv$ と $\bar z$ は異なる固有値に属し、複素化した行列も正規だから直交する。よって
$$
0=\inner{z}{\bar z}=\sum_j z_j^2=\norm u^2-\norm v^2+2i\,\transp u v.
$$
この実部・虚部から $u,v$ は直交し、長さは等しく非零である。$Az=\lambda z$ の実部と虚部は $Au=au-bv,\ Av=bu+av$ であるから、$(u,-v)$ をそれぞれ規格化した順序では制限の行列が所定の二次ブロックになる。その平面の直交補も不変で制限は正規である。

いずれの場合も直交補に同じ操作を繰り返すと次元が下がり、零次元で終了する。二次ブロックの特性多項式は $(t-a)^2+b^2$ であり、各ブロックの特性多項式の積が全体の特性多項式である。逆に二次ブロック $C$ は $\transp C C=C\transp C=(a^2+b^2)I$ を満たし、一次ブロックも正規である。直和と直交相似はこの等式を保つ。
:::

$A$ が実対称なら非実固有値がないので $2\times2$ ブロックは現れず、これは実スペクトル定理に一致する。
$A$ が交代行列（$\transp{A}=-A$）なら固有値は純虚数または $0$ なので、$1\times1$ ブロックは $0$、$2\times2$ ブロックは $a=0$ の $\begin{pmatrix}0&-b\\b&0\end{pmatrix}$ である。

:::corollary[直交行列の標準形]{#cor-linalg-orthogonal-normal-form}
$A\in\mathrm{O}(n)$ ならば、$Q\in\mathrm{O}(n)$ で

$$
\transp{Q}AQ=\diag(I_p,\,-I_q,\,R_{\theta_1},\ldots,R_{\theta_m}),\qquad 0<\theta_j<\pi
$$

となるものが存在する。
ここで $R_\theta$ は角 $\theta$ の回転行列であり、$\det A=(-1)^q$ である。
:::

:::proof
直交行列は正規なので :ref[prop-linalg-real-normal-form] が適用できる。
$\transp{Q}AQ$ は直交行列であり（:ref[prop-linalg-unitary-group] の 1）、ブロック対角行列が直交行列であることは各ブロックが直交行列であることと同値である。
$1\times1$ ブロック $\lambda$ は $\lambda^2=1$ を満たすので $\pm1$ である。
$2\times2$ ブロック $\begin{pmatrix}a&-b\\b&a\end{pmatrix}$ の列の長さが $1$ であることから $a^2+b^2=1$ であり、$b>0$ より $a=\cos\theta,\ b=\sin\theta$（$0<\theta<\pi$）と書ける。
ブロックの順序は基底の並べ替え（置換行列は直交行列）で調整できる。
$\det R_\theta=1$ なので $\det A=(-1)^q$ である。
:::

### 3 次元の回転

:::theorem[3 次元の回転の軸]{#thm-linalg-rotation-axis}
$A\in\SO(3)$ ならば、$Q\in\SO(3)$ と $0\le\theta\le\pi$ で

$$
\transp{Q}AQ=\begin{pmatrix}1&0&0\\0&\cos\theta&-\sin\theta\\0&\sin\theta&\cos\theta\end{pmatrix}
$$

となるものが存在する。
特に $Q$ の第 1 列 $u$ は $Au=u$ を満たし、$A$ は $u$ を軸とする角 $\theta$ の回転である。
角は $\tr A=1+2\cos\theta$ で定まり、$A\ne I$ なら回転軸 $\spn\{u\}=\Ker(A-I)$ も一意である。
:::

:::proof
:ref[cor-linalg-orthogonal-normal-form] で $n=3$ とすると、$p+q+2m=3$ かつ $(-1)^q=\det A=1$ である。
$m=1$ なら $p+q=1$ で $q$ は偶数なので $p=1,q=0$ であり、$\diag(1,R_\theta)$ の形になる。
$m=0$ なら $q\in\{0,2\}$ であり、$q=0$ は $A=I$（$\theta=0$）、$q=2$ は $\diag(1,-1,-1)=\diag(1,R_\pi)$（$\theta=\pi$）である。
得られた $Q$ が $\det Q=-1$ なら、第 1 列 $u$ を $-u$ に替える。
これは $Q$ に右から $D=\diag(-1,1,1)$ を掛けることであり、$D\diag(1,R_\theta)D=\diag(1,R_\theta)$ なので標準形は変わらず、$\det(QD)=1$ となる。

トレースは相似で不変なので（:ref[prop-linalg-similarity-invariants]）、$\tr A=1+2\cos\theta$ である。
$0\le\theta\le\pi$ では $\cos$ は単射なので $\theta$ は $A$ から決まる。
$A\ne I$ なら $\theta>0$ であり、標準形から $\Ker(A-I)$ は 1 次元で $u$ により張られる。
:::

角 $\theta$ の符号は、軸の向き $u$ と、$Q=(u\ w_1\ w_2)$ が $\det Q=1$ を満たすという向き付けの規約で決まる。
軸の向きを逆にすると、同じ回転が角 $-\theta$ の回転として記述される。

:::example[回転の軸と角を求める]{#exm-linalg-rotation-axis}
次の行列を考える。

$$
A=\frac13\begin{pmatrix}2&-1&2\\2&2&-1\\-1&2&2\end{pmatrix}
$$

列は互いに直交する単位ベクトルであり、$\det A=1$ なので $A\in\SO(3)$ である。
各行の和が $1$ なので $A(1,1,1)=(1,1,1)$ であり、軸は $u=(1,1,1)/\sqrt3$ である。
$\tr A=2=1+2\cos\theta$ より $\cos\theta=1/2$、$\theta=\pi/3$ である。
向きを確かめるため、$w_1=(1,-1,0)/\sqrt2$、$w_2=(1,1,-2)/\sqrt6$ とおく。
$Q=(u\ w_1\ w_2)$ は直交行列で $\det Q=1$ である。
$Aw_1=(1,0,-1)/\sqrt2$ より

$$
\inner{Aw_1}{w_1}=\frac12=\cos\frac\pi3,\qquad \inner{Aw_1}{w_2}=\frac{1+2}{\sqrt{12}}=\frac{\sqrt3}2=\sin\frac\pi3
$$

であり、$\transp{Q}AQ=\diag(1,R_{\pi/3})$ となる。
:::

## 二次形式と Sylvester の慣性法則

$2x^2+2xy+2y^2$ のような二次式は、座標を回転すれば $3u^2+v^2$ のように交差項のない形になった。
固有値を決めるためには軸の長さを保つ必要があったが、二次式の符号の様子だけが知りたいときは、軸を伸縮して係数の絶対値を $1$ にしてよい。
この場合に何が不変量として残るかを調べる。
この節の形式の理論は内積を必要とせず、$V$ は $K$ 上の有限次元ベクトル空間であればよい。
内積は主軸変換のときにだけ用いる。

### 双線形形式・Hermite 形式と表現行列

:::definition[対称双線形形式・Hermite 形式・二次形式]{#def-linalg-hermitian-form}
$V$ を $K$ 上の $n$ 次元ベクトル空間とする。

1. $K=\R$ のとき、写像 $b:V\times V\to\R$ が各引数について線形で $b(x,y)=b(y,x)$ を満たすとき、$b$ を :term[対称双線形形式]{en="symmetric bilinear form" yomi="たいしょうそうせんけいけいしき"} という。
2. $K=\C$ のとき、写像 $b:V\times V\to\C$ が第 1 引数について線形、第 2 引数について共役線形（:term[半双線形形式]{en="sesquilinear form" yomi="はんそうせんけいけいしき"}）で、$b(y,x)=\conj{b(x,y)}$ を満たすとき、$b$ を :term[Hermite 形式]{en="Hermitian form" yomi="えるみーとけいしき"} という。

いずれの場合も $q(x)=b(x,x)$ は実数値であり、これを $b$ に付随する :term[二次形式]{en="quadratic form" yomi="にじけいしき"} という。
基底 $(v_1,\ldots,v_n)$ に対し $H_{ij}=b(v_j,v_i)$ を成分とする行列 $H$ を $b$ の表現行列という。
:::

内積は正定値な対称双線形形式（$K=\R$）または正定値な Hermite 形式（$K=\C$）にほかならない。
$x=\sum_jc_jv_j$、$y=\sum_id_iv_i$ とすれば

$$
b(x,y)=\sum_{i,j}\bar d_ic_jb(v_j,v_i)=\adj{d}Hc,\qquad q(x)=\adj{c}Hc
$$

であり、$b(v_i,v_j)=\conj{b(v_j,v_i)}$ から $\adj{H}=H$ である。
すなわち表現行列は Hermite 行列（$K=\R$ では対称行列）である。
逆に Hermite 行列 $H$ を与えれば、$b(x,y)=\adj{d}Hc$ は Hermite 形式を定める。
新しい基底 $v'_k=\sum_js_{jk}v_j$ に移ると、新しい座標 $c'$ と古い座標の関係は $c=Sc'$ なので

$$
b(x,y)=\adj{(Sd')}H(Sc')=\adj{d'}(\adj{S}HS)c'
$$

となり、表現行列は $\adj{S}HS$ に変わる。

:::definition[合同]{#def-linalg-congruence}
Hermite 行列（実対称行列）$H,H'$ について、正則行列 $S$（$K=\R$ では実の正則行列）で $H'=\adj{S}HS$ となるものが存在するとき、$H$ と $H'$ は :term[合同]{en="congruent" yomi="ごうどう"} であるという。
:::

合同は同値関係であり、同じ形式を異なる基底で表した行列どうしの関係にほかならない。
相似 $S^{-1}HS$ とは異なり、合同は固有値を保たない（例えば $\adj{S}HS$ で $S=cI$ とすれば行列は $\abs{c}^2$ 倍になる）。
ただし $S$ がユニタリなら $\adj{S}=S^{-1}$ なので、合同と相似は一致する。

形式 $b$ は二次形式 $q$ だけから復元できる。
$K=\R$ では $q(x+y)-q(x-y)=4b(x,y)$ であり、$K=\C$ では :ref[prop-linalg-quadratic-zero] の証明と同じ計算により

$$
b(x,y)=\frac14\sum_{k=0}^3i^k\,q(x+i^ky)
$$

である。
これらは内積の偏極恒等式（:ref[prop-linalg-polarization]）の一般化であり、証明には正定値性を使わない。
なお、複素空間の双線形形式 $\transp{x}Ax$ は実数値とは限らないので、以下の符号による分類の対象は実対称形式と複素 Hermite 形式である。

:::proposition[形式と自己随伴作用素の対応]{#prop-linalg-form-operator}
$V$ を内積空間とする。
$V$ 上の Hermite 形式（$K=\R$ では対称双線形形式）$b$ に対し、$b(x,y)=\inner{Tx}{y}$ がすべての $x,y$ で成り立つ自己随伴作用素 $T$ がただ一つ存在する。
逆に自己随伴作用素 $T$ に対し $\inner{Tx}{y}$ は Hermite 形式である。
正規直交基底に関する $b$ の表現行列は、同じ基底に関する $T$ の表現行列に等しい。
:::

:::proof
$x$ を固定し、線形汎関数 $y\mapsto\overline{b(x,y)}$ に :ref[thm-linalg-riesz-finite] を適用する。共役を取ると $b(x,y)=\inner{Tx}{y}$ となるベクトル $Tx$ がただ一つある。第 1 引数の線形性と表示の一意性から $T$ は線形で、$\inner{Tx}{y}=\overline{b(y,x)}=\inner{x}{Ty}$ から自己随伴である。逆にこの等式は自己随伴な $T$ が定める形式の共役対称性も与える。正規直交基底では $Te_j$ の第 $i$ 成分が $\inner{Te_j}{e_i}=b(e_j,e_i)$ なので、二つの表現行列は一致する。
:::

:::theorem[主軸定理]{#thm-linalg-principal-axes}
$V$ を内積空間、$b$ を Hermite 形式（$K=\R$ では対称双線形形式）とする。
$V$ の正規直交基底 $(e_1,\ldots,e_n)$ と実数 $\lambda_1,\ldots,\lambda_n$ で

$$
b(x,x)=\sum_{j=1}^n\lambda_j\abs{\inner{x}{e_j}}^2\qquad(x\in V)
$$

となるものが存在する。
$\lambda_j$ は、任意の正規直交基底に関する $b$ の表現行列の固有値（重複度込み）である。
特に実対称行列 $A$ に対し、直交行列 $Q$ による座標変換 $x=Qy$ で

$$
\transp{x}Ax=\lambda_1y_1^2+\cdots+\lambda_ny_n^2
$$

とできる。
:::

:::proof
:ref[prop-linalg-form-operator] の自己随伴作用素 $T$ に :ref[thm-linalg-spectral]（$K=\C$）または :ref[thm-linalg-real-spectral]（$K=\R$）を適用し、固有ベクトルからなる正規直交基底 $(e_j)$、$Te_j=\lambda_je_j$ をとる。
$x=\sum_jc_je_j$、$c_j=\inner{x}{e_j}$ とすれば $b(x,x)=\inner{Tx}{x}=\sum_j\lambda_j\abs{c_j}^2$ である。
任意の正規直交基底に関する $b$ の表現行列は $T$ の表現行列なので、その固有値は $T$ の固有値 $\lambda_j$ である。
行列の主張は、$\R^n$ の標準内積と $b(x,y)=\transp{y}Ax$ に適用し、$Q=(e_1\ \cdots\ e_n)$ とすればよい。
:::

正規直交基底による座標変換で二次形式を対角化することを :term[主軸変換]{en="principal axis transformation" yomi="しゅじくへんかん"} といい、$e_j$ の方向を主軸という。

### Sylvester の慣性法則

:::theorem[Sylvester の慣性法則]{#thm-linalg-inertia}
任意の Hermite 行列（実対称行列）$H\in M_n(K)$ は、ある $p,q,r\ge0$（$p+q+r=n$）について $\diag(I_p,-I_q,O_r)$ と合同である。
整数 $p,q,r$ は $H$ によって一意に定まる。
さらに $p$ と $q$ は $H$ の正の固有値と負の固有値の個数（重複度込み）に等しく、$r=n-\rank H$ である。
:::

:::proof
正規直交固有基底に移り、実固有値 $\lambda_j$ を正・負・零の順に並べる。非零固有値の基底ベクトルだけを $|\lambda_j|^{-1/2}$ 倍すると、形式の対応する係数が $+1$ または $-1$ になり、所定の標準形を得る。

正定値に制限できる部分空間の最大次元は合同不変である。実際、$H'=\adj S HS$ では $x$ での $H'$ の値が $Sx$ での $H$ の値に等しく、$S$ と $S^{-1}$ がそのような部分空間を同じ次元で対応させる。標準形で先頭 $p$ 本の座標軸が張る空間は正定値である。一方、次元が $p$ より大きい部分空間では、先頭 $p$ 座標への射影に非自明な核があり、その核上の値は $-\sum_{p<j\le p+q}|x_j|^2\le0$ である。よって最大次元は $p$ に等しい。$-H$ に同じ議論を用いれば $q$ も決まり、$r=n-p-q$ も決まる。構成から $p,q$ はそれぞれ正・負の固有値の個数であり、$r$ は零固有空間の次元 $n-\rank H$ である。
:::

三つ組 $(p,q,r)$ を $H$（または形式 $b$）の :term[慣性]{en="inertia" yomi="かんせい"} といい、$p$ を正の慣性指数、$q$ を負の慣性指数という。
$(p,q)$ あるいは差 $p-q$ を :term[符号数]{en="signature" yomi="ふごうすう"} と呼ぶ流儀もある。
慣性は形式 $b$ そのものの不変量であり、どの基底で表現行列を計算しても同じ値になる。

:::corollary[合同による分類]{#cor-linalg-congruence-classes}
二つの $n$ 次 Hermite 行列（実対称行列）が合同であるための必要十分条件は、慣性が一致することである。
したがって合同類は $p+q+r=n$ を満たす非負整数の組 $(p,q,r)$ と一対一に対応し、その個数は $\binom{n+2}2$ である。
:::

:::proof
慣性が等しければ、どちらも同じ $\diag(I_p,-I_q,O_r)$ と合同であり、合同は同値関係なので両者は合同である。
逆は :ref[thm-linalg-inertia] の一意性である。
組の個数は、$n$ 個の玉と 2 本の仕切りの並べ方の数 $\binom{n+2}2$ に等しい。
:::

### 平方完成による標準形

慣性を求めるのに固有値は要らない。
:term[平方完成]{en="completing the square" yomi="へいほうかんせい"} を繰り返す Lagrange の方法は、正則な座標変換だけで形式を対角化する。
$K=\R$ で述べる。

1. ある対角成分、例えば $h_{11}$ が零でなければ、$x_1$ を含む項をまとめて $h_{11}\bigl(x_1+\sum_{j\ge2}\frac{h_{1j}}{h_{11}}x_j\bigr)^2$ と残りの $x_2,\ldots,x_n$ だけの二次形式の和に書く。新しい変数 $y_1=x_1+\sum_{j\ge2}(h_{1j}/h_{11})x_j$、$y_j=x_j$（$j\ge2$）への変換は正則である。
2. 対角成分がすべて零でも非零の非対角成分があれば、添字を替えて $h_{12}\ne0$ とする。このとき、$x_1=y_1+y_2,\ x_2=y_1-y_2$ と置くと $2h_{12}x_1x_2=2h_{12}(y_1^2-y_2^2)$ となり、対角成分が現れる。

残りの変数について同じ操作を繰り返せば、正則な変換で $\sum_jc_jy_j^2$ の形になる。
残りの行列が零ならそこで止める。$c_j\ne0$ の座標だけを $\abs{c_j}^{-1/2}$ 倍し、符号順に並べれば Sylvester の標準形になる。

:::example[平方完成で慣性を読む]
実二次形式を一変数ずつ平方完成する。

$$
q(x,y,z)=x^2+4xy+3y^2=(x+2y)^2-y^2
$$

$u=x+2y,\ v=y,\ w=z$ は可逆な座標変換なので慣性は $(1,1,1)$ である。
$q=0$ を満たすベクトルは核の方向だけではない。
例えば $(x,y,z)=(-1,1,0)$ は $q=0$ を満たすが、表現行列を掛けると $(1,1,0)\ne0$ になる。
:::

:::example[対角成分が消えているとき]
$q(x,y)=2xy$ を座標軸の上だけで調べると値はどちらも零である。
しかし $x=u+v,\ y=u-v$ とすれば $q=2u^2-2v^2$ であり、慣性は $(1,1,0)$ である。
主軸変換 $u'=(x+y)/\sqrt2,\ v'=(x-y)/\sqrt2$ でも $q=u'^2-v'^2$ となり、どちらの方法でも同じ慣性を得る。
対角成分の符号だけでは判定できない。
:::

### 零方向と形式の核を区別する

形式 $b$ の :term[根基]{en="radical" yomi="こんき"} とは、すべての $y$ に対して $b(x,y)=0$ となる $x$ の全体である。
座標で書けば $\adj{d}Hc=0$ がすべての $d$ で成り立つことなので、$Hc=0$ と同値であり、根基の次元は $n-\rank H=r$ である。
根基が $\{0\}$、すなわち $\det H\ne0$ のとき、$b$ は :term[非退化]{en="nondegenerate" yomi="ひたいか"} であるという。
一方、$q(x)=b(x,x)=0$ という条件は相手を $y=x$ に限定したものなので、一般に根基に属するという条件より弱い。

:::example[不定値形式の零集合]
$q(x,y,z)=x^2-y^2$ の標準形はすでに $\diag(1,-1,0)$ で、根基は $z$ 軸である。
ところが $q=0$ となる点は

$$
\{x=y\}\cup\{x=-y\}
$$

という二つの平面の和集合になる。
例えば $(1,1,0)$ と $(1,-1,0)$ はともにこの集合に属するが、和 $(2,0,0)$ の値は $4$ である。
したがって零集合は部分空間ですらない。
半正定値の標準形では非負の平方の和が零になるのは各平方が零のときだけなので、半正定値形式では零集合と根基が一致する（演習問題を参照）。
:::

### 二次曲線と二次曲面の分類

主軸変換の古典的な応用として、$\R^n$ の二次曲面

$$
\transp{x}Ax+2\transp{b}x+c=0,\qquad A=\transp{A}\ne O,\ b\in\R^n,\ c\in\R
$$

を、距離を保つ座標変換 $x=Qy+x_0$（$Q\in\mathrm{O}(n)$、$x_0\in\R^n$）で標準形に直す。
このような変換を :term[合同変換]{en="rigid motion" yomi="ごうどうへんかん"}（運動）といい、図形の形と大きさを変えない。
行列の合同と同じ言葉だが別の概念なので注意する。

:::theorem[二次曲面の標準形]{#thm-linalg-quadric-normal-form}
$k=\rank A$ とし、$A$ の非零固有値を $\lambda_1,\ldots,\lambda_k$ とする。
合同変換により、上の方程式は次のいずれかの形になる。

1. $\lambda_1y_1^2+\cdots+\lambda_ky_k^2+c'=0$（$c'\in\R$）。
2. $\lambda_1y_1^2+\cdots+\lambda_ky_k^2+2\beta y_{k+1}=0$（$k<n$、$\beta>0$）。
:::

:::proof
主軸変換 $x=Qz$ で二次部分を $\sum_{j\le k}\lambda_j z_j^2$ にし、$b'=\transp Q b$ とする。$j\le k$ で $w_j=z_j+b'_j/\lambda_j$、残りでは $w_j=z_j$ とおくと、平方完成により
$$
\sum_{j\le k}\lambda_jw_j^2+2\sum_{j>k}b'_jw_j+c''=0,\qquad
c''=c-\sum_{j\le k}(b'_j)^2/\lambda_j
$$
を得る。残った係数列 $d=(b'_{k+1},\ldots,b'_n)$ が零なら第 1 の形である。非零なら核の座標内で $d/\norm d$ を先頭とする正規直交基底を選ぶ。その座標を $v$ とすれば一次の部分は $2\norm d\,v_{k+1}$ となる。$y_{k+1}=v_{k+1}+c''/(2\norm d)$ の平行移動で定数を吸収すると第 2 の形になる。いずれの変換も直交変換または平行移動である。
:::

$n=2$ の場合（二次曲線）を、$\det A=\lambda_1\lambda_2$ の符号で整理する。

- $\det A>0$（楕円型）：$\lambda_1y_1^2+\lambda_2y_2^2+c'=0$ で $\lambda_1,\lambda_2$ は同符号である。$c'$ が $\lambda_j$ と逆符号なら :term[楕円]{en="ellipse" yomi="だえん"}、$c'=0$ なら 1 点、同符号なら空集合である。
- $\det A<0$（双曲型）：$\lambda_1,\lambda_2$ は異符号である。$c'\ne0$ なら :term[双曲線]{en="hyperbola" yomi="そうきょくせん"}、$c'=0$ なら交わる 2 直線である。
- $\det A=0$（放物型）：$k=1$ である。2 の形なら $\lambda_1y_1^2+2\beta y_2=0$ で :term[放物線]{en="parabola" yomi="ほうぶつせん"} である。1 の形なら $\lambda_1y_1^2+c'=0$ で、平行な 2 直線、1 本の直線（重複）、空集合のいずれかである。

:::example[楕円の主軸]{#exm-linalg-ellipse}
$5x^2-4xy+8y^2=36$ を考える。
$A=\begin{pmatrix}5&-2\\-2&8\end{pmatrix}$ の特性多項式は $t^2-13t+36=(t-4)(t-9)$ であり、$A(2,1)=4(2,1)$、$A(1,-2)=9(1,-2)$ である。
$u=(2x+y)/\sqrt5,\ v=(x-2y)/\sqrt5$ とすれば方程式は

$$
4u^2+9v^2=36,\qquad\text{すなわち}\qquad \frac{u^2}{9}+\frac{v^2}{4}=1
$$

となる。
これは $(2,1)/\sqrt5$ 方向の半軸の長さが $3$、$(1,-2)/\sqrt5$ 方向の半軸の長さが $2$ の楕円である。
:::

:::example[放物線の標準形]
$x^2+2xy+y^2+\sqrt2x-\sqrt2y=0$ を考える。
$A=\begin{pmatrix}1&1\\1&1\end{pmatrix}$ の固有値は $2$（方向 $(1,1)/\sqrt2$）と $0$（方向 $(1,-1)/\sqrt2$）で、$\det A=0$ である。
$u=(x+y)/\sqrt2,\ v=(x-y)/\sqrt2$ と置けば $(x+y)^2=2u^2$、$\sqrt2(x-y)=2v$ なので、方程式は $2u^2+2v=0$、すなわち $v=-u^2$ となる。
これは放物線である。
:::

$n=3$ では、$c'\ne0$ の 1 の形 $\sum_j\mu_jy_j^2=1$（$\mu_j=\lambda_j/(-c')$）は、$\diag(\mu_j)$ の慣性が $(3,0,0)$ なら :term[楕円面]{en="ellipsoid" yomi="だえんめん"}、$(2,1,0)$ なら :term[一葉双曲面]{en="hyperboloid of one sheet" yomi="いちようそうきょくめん"}、$(1,2,0)$ なら :term[二葉双曲面]{en="hyperboloid of two sheets" yomi="にようそうきょくめん"}、$(0,3,0)$ なら空集合である。
$c'=0$ で慣性が $(2,1,0)$ または $(1,2,0)$ なら二次錐面である。
2 の形で $k=2$ なら、$\lambda_1,\lambda_2$ が同符号のとき楕円放物面、異符号のとき双曲放物面である。

:::example[二次曲面の判定]
$2xy+2yz+2zx=1$ の行列 $\begin{pmatrix}0&1&1\\1&0&1\\1&1&0\end{pmatrix}=J-I$ の固有値は $3-1=2$ と $0-1=-1$（重複度 $2$）である（:ref[exm-linalg-symmetric-3x3] の $J$ を用いた）。
主軸変換で $2u^2-v^2-w^2=1$ となり、慣性 $(1,2,0)$ なので二葉双曲面である。
$u$ 軸は $(1,1,1)$ の方向であり、二つの葉は $u\ge1/\sqrt2$ と $u\le-1/\sqrt2$ にある。
:::

## 正定値行列

エネルギーや距離の二乗に使う形式は、どの非零方向でも正でなければならない。
慣性の言葉では $p=n$ の場合であるが、固有値をすべて求める手間を省き、行列式だけで判定する方法がある。

:::definition[正定値と半正定値]{#def-linalg-positive-definite}
Hermite 行列 $H\in M_n(K)$（$K=\R$ では対称行列）について、すべての $x\ne0$ で $\adj{x}Hx>0$ となるとき :term[正定値]{en="positive definite" yomi="せいていち"}、すべての $x$ で $\adj{x}Hx\ge0$ となるとき :term[半正定値]{en="positive semidefinite" yomi="はんせいていち"} といい、それぞれ $H>0,\ H\ge0$ と書く。
$-H$ が正定値（半正定値）のとき $H$ を負定値（半負定値）といい、$\adj{x}Hx$ が正の値も負の値もとるとき :term[不定値]{en="indefinite" yomi="ふていち"} という。
自己随伴作用素 $T$ についても、$\inner{Tx}{x}$ を用いて同様に定める。
:::

:ref[thm-linalg-inertia] の言葉では、$H>0$ は慣性 $(n,0,0)$、$H\ge0$ は $q=0$ と同値である。
$H$ の左上 $k\times k$ 部分行列を $H_{[k]}$ とし、その行列式 $\Delta_k=\det H_{[k]}$ を左上 $k$ 次の :term[首座小行列式]{en="leading principal minor" yomi="しゅざしょうぎょうれつしき"} という（$\Delta_0=1$ とおく）。
より一般に、添字の集合 $J\subset\{1,\ldots,n\}$ に対し、$J$ に属する行と列だけを残した部分行列 $H_J$ の行列式を :term[主小行列式]{en="principal minor" yomi="しゅしょうぎょうれつしき"} という。

### 首座小行列式による判定と Cholesky 分解

:::theorem[正定値性の判定と Cholesky 分解]{#thm-linalg-positive-criteria}
Hermite 行列 $H\in M_n(K)$ について次は同値である。

1. $H>0$。
2. $H$ のすべての固有値が正である。
3. $\Delta_k>0$（$k=1,\ldots,n$）。
4. 対角成分が正の実数である上三角行列 $R\in M_n(K)$ により $H=\adj{R}R$ と書ける。

4 の $R$ は一意であり、$H=\adj{R}R$ を $H$ の :term[Cholesky 分解]{en="Cholesky decomposition" yomi="これすきーぶんかい"} という。
:::

:::proof
正規直交固有基底で形式は $\sum_j\lambda_j|c_j|^2$ だから、各座標軸での評価と和の符号により 1 と 2 は同値である。正定値な形式を最初の $k$ 座標の部分空間に制限しても正定値である。その行列式は正の固有値の積だから $\Delta_k>0$ となる。

逆に全首座小行列式が正とする。$H=\begin{pmatrix}a&\adj b\\ b&C\end{pmatrix}$ と書くと $a=\Delta_1>0$ である。$S=C-b\adj b/a$ とおき、第 1 行の倍数を各行から引けば、左下ブロックが零、右下が $S$ になる。同じ行操作がどの首座部分にも使えるため $\Delta_{j+1}(H)=a\Delta_j(S)$ であり、$S$ も条件 3 を満たす。一次元では $R=(\sqrt a)$ である。次元について帰納し、$S=\adj{R_0}R_0$ と分解すれば
$$
R=\begin{pmatrix}\sqrt a&\adj b/\sqrt a\\0&R_0\end{pmatrix},
\qquad \adj R R=\begin{pmatrix}a&\adj b\\b&b\adj b/a+S\end{pmatrix}=H
$$
となり 4 を得る。逆に 4 の $R$ は正則なので $x\ne0$ に対し $\adj x Hx=\norm{Rx}^2>0$ である。

一意性も同じ消去から従う。正の対角を指定すると第 1 行の先頭は $\sqrt a$、残りは $\adj b/\sqrt a$ に決まる。残りの三角因子は $S$ の分解だから帰納法で一意である。
:::

:::example[パラメータを含む判定]
$H=\begin{pmatrix}2&1\\1&t\end{pmatrix}$ は $\Delta_1=2,\ \Delta_2=2t-1$ より $t>1/2$ で正定値である。
実際、形式は $2\abs{x+y/2}^2+(t-1/2)\abs{y}^2$ であり、$t=1/2$ では $(1,-2)$ 方向で値が零になる。
Cholesky 分解は $R=\begin{pmatrix}\sqrt2&1/\sqrt2\\0&\sqrt{t-1/2}\end{pmatrix}$ である。
:::

:::corollary[負定値性の判定]{#cor-linalg-negative-definite}
Hermite 行列 $H$ が負定値であるための必要十分条件は、$(-1)^k\Delta_k>0$（$k=1,\ldots,n$）、すなわち $\Delta_1<0,\ \Delta_2>0,\ \Delta_3<0,\ldots$ と符号が交互に変わることである。
:::

:::proof
$H$ が負定値であることは $-H>0$ と同値である。
$(-H)_{[k]}=-H_{[k]}$ の行列式は :ref[prop-linalg-det-operations] の 2 により $(-1)^k\Delta_k$ である。
:ref[thm-linalg-positive-criteria] を $-H$ に適用すればよい。
:::

首座小行列式がすべて零でなければ、正定値・負定値でない場合も含めて慣性まで読み取れる。

:::theorem[Jacobi の符号法則]{#thm-linalg-jacobi-signature}
Hermite 行列 $H\in M_n(K)$ の首座小行列式 $\Delta_1,\ldots,\Delta_n$ がすべて零でないとする。
このとき対角成分が $1$ の下三角行列 $L$ により

$$
H=LD\adj{L},\qquad D=\diag\Bigl(\frac{\Delta_1}{\Delta_0},\frac{\Delta_2}{\Delta_1},\ldots,\frac{\Delta_n}{\Delta_{n-1}}\Bigr)
$$

と書ける。
したがって $H$ の負の慣性指数 $q$ は数列 $\Delta_0=1,\Delta_1,\ldots,\Delta_n$ の符号の変化の回数に等しく、$p=n-q$、$r=0$ である。
:::

:::proof
$H=\begin{pmatrix}a&\adj b\\ b&C\end{pmatrix}$ で $a=\Delta_1\ne0$ とし、$E=\begin{pmatrix}1&0\\-b/a&I\end{pmatrix}$ とおく。直接掛けると $EH\adj E=\diag(a,S),\ S=C-b\adj b/a$ である。この消去を左上部分に制限しても成り立つから $\Delta_j(S)=\Delta_{j+1}(H)/a$ であり、残りの首座小行列式も非零である。一次元なら $L=1,D=a$ である。帰納的に $S=L_0D_0\adj{L_0}$ とすれば $L=E^{-1}\diag(1,L_0)$ は単位下三角で、$H=L\diag(a,D_0)\adj L$ を得る。$S$ の相続く首座小行列式の比は $\Delta_{j+1}(H)/\Delta_j(H)$ に等しいため、所定の対角が得られる。合同で慣性が保たれるので、負の対角成分を数えればよい。これは相続く $\Delta_j$ の符号が変わる回数である。
:::

正定値（すべて $\Delta_k>0$）と負定値（交互の符号）は、この定理の特別な場合とみなせる。
多変数関数の極値判定では、臨界点での Hesse 行列にこれらの判定法を適用する（[多変数関数の微分](/math/analysis/multivariable-differentiation/) を参照）。

:::example[首座小行列式から慣性を読む]{#exm-linalg-jacobi}
$H=\begin{pmatrix}1&2&0\\2&1&1\\0&1&1\end{pmatrix}$ では

$$
\Delta_1=1,\qquad \Delta_2=1-4=-3,\qquad \Delta_3=1\cdot(1-1)-2\cdot(2-0)=-4
$$

である。
数列 $1,1,-3,-4$ の符号の変化は 1 回なので、慣性は $(2,1,0)$ である。
実際、$\Phi_H=(t-1)(t^2-2t-4)$ であり、固有値は $1,\ 1\pm\sqrt5$ で、正が二つ、負が一つである。
:::

いずれかの $\Delta_k$ が零になると Jacobi の符号法則は使えない。
例えば $\diag(0,-1)$ と $\diag(0,1)$ はともに $\Delta_1=\Delta_2=0$ であるが、前者は半正定値でなく、後者は半正定値である。
半正定値性を首座小行列式だけで判定することはできず、主小行列式による必要十分条件では首座以外も含める。
その証明のために、次の展開公式を用意する。

:::lemma[主小行列式による展開]{#lem-linalg-principal-minor-expansion}
$B\in M_k(K)$ と $\eps\in K$ に対し

$$
\det(B+\eps I_k)=\sum_{j=0}^{k}\eps^{k-j}E_j(B)
$$

が成り立つ。
ここで $E_j(B)$ は $B$ の $j$ 次主小行列式すべての和であり、$E_0(B)=1$ とする。
:::

:::proof
第 $i$ 列を $b_i+\eps e_i$ と書き、多重線形性で展開する。$b_i$ を選んだ添字集合を $J$ とすると、他の列からは $\eps^{k-|J|}$ が出る。残った単位列ごとに余因子展開すると、対角位置の $1$ だけが寄与し、行と列を同じ位置で消すので符号は正である。残る行列式は $\det B_J$（$J=\varnothing$ では $1$）となる。$|J|=j$ ごとに集めれば公式を得る。
:::

同じ計算を多項式環 $K[t]$ 上の行列で行えば（:ref[prop-linalg-det-ring]）、$\Phi_B(t)=\det(tI-B)=\sum_{j=0}^k(-1)^jE_j(B)t^{k-j}$ も得られる。
これは :ref[prop-linalg-charpoly-coefficients] の係数の表示を完成させるものである。

:::proposition[半正定値性の判定]{#prop-linalg-psd-criteria}
Hermite 行列 $H\in M_n(K)$ について次は同値である。

1. $H\ge0$。
2. $H$ のすべての固有値が非負である。
3. $H$ のすべての主小行列式が非負である。
4. ある $B\in M_n(K)$ により $H=\adj{B}B$ と書ける。
5. ある内積空間のベクトル $x_1,\ldots,x_n$ により $H_{ij}=\inner{x_j}{x_i}$ と書ける（$H$ は Gram 行列である）。

さらに、$H>0$ であることは、5 で $x_1,\ldots,x_n$ を一次独立にとれることと同値であり、4 で $B$ を正則にとれることとも同値である。
:::

:::proof
固有基底で $\adj x Hx=\sum_j\lambda_j|c_j|^2$ だから 1 と 2 は同値である。1 なら任意の座標部分空間への制限 $H_J$ も半正定値であり、その行列式は非負固有値の積なので 3 が従う。逆に 3 なら :ref[lem-linalg-principal-minor-expansion] よりどの $k$ でも
$$
\det(H_{[k]}+\eps I_k)=\sum_{j=0}^k\eps^{k-j}E_j(H_{[k]})\ge\eps^k>0\qquad(\eps>0).
$$
したがって $H+\eps I>0$ であり、任意の非零 $x$ に対し $\adj x Hx>-\eps\norm x^2$ となる。$\eps\downarrow0$ とすれば 1 を得る。零ベクトルでの値は零である。

2 なら $H=UD\adj U$ で $B=\diag(\sqrt{\lambda_j})\adj U$ とおけば $\adj B B=H$ である。4 の行列 $B$ の列を $x_j$ と取れば $(\adj B B)_{ij}=\inner{x_j}{x_i}$ となり 5 が従う。5 なら $\adj c Hc=\norm{\sum_j c_j x_j}^2\ge0$ で 1 に戻る。この式で非零の $c$ に対し常に正となることは列の独立性と同値で、正方行列 $B$ の場合は正則性と同値である。
:::

5 は、内積空間の Gram 行列（[内積空間](/math/linear-algebra/inner-product-spaces/) を参照）がちょうど半正定値行列であることを述べている。

### 正の平方根と極分解

固有値ごとに非負の平方根を選ぶ。

:::theorem[半正定値平方根と極分解]{#thm-linalg-polar}
1. 半正定値行列 $H$ に対し、$B\ge0$ かつ $B^2=H$ を満たす行列 $B$ がただ一つ存在する。これを $H^{1/2}$ と書き、$H$ の半正定値平方根という。$H^{1/2}$ は $H$ の多項式で表され、$H>0$ なら $H^{1/2}>0$ である。
2. 任意の $A\in M_n(K)$ は、$P=(\adj{A}A)^{1/2}$ とユニタリ行列（$K=\R$ では直交行列）$U$ により $A=UP$ と書ける。これを $A$ の :term[極分解]{en="polar decomposition" yomi="きょくぶんかい"} という。$P$ は常に一意であり、$A$ が正則なら $U$ も一意である。
:::

:::proof
スペクトル分解 $H=\sum_j\lambda_jP_j$ で $\lambda_j\ge0$ だから、$B=\sum_j\sqrt{\lambda_j}P_j$ は半正定値で $B^2=H$ である。各射影は :ref[prop-linalg-spectral-projections] の補間式で $H$ の多項式となる。この式は実固有値には実係数を与える。$H>0$ なら $\sqrt{\lambda_j}>0$ である。別の半正定値平方根 $B'$ は自身の二乗 $H$ と可換なので、各 $E_{\lambda_j}(H)$ を保つ。その上で $B'$ は自己随伴で固有基底をもち、各固有値 $\mu$ は $\mu\ge0,\ \mu^2=\lambda_j$ を満たす。よってこの空間上で $B'=\sqrt{\lambda_j}I$ となり、全体で $B'=B$ である。

$\adj A A\ge0$ なので $P=(\adj A A)^{1/2}$ が存在する。$\norm{Px}^2=\inner{P^2x}{x}=\norm{Ax}^2$ より $\Ker P=\Ker A$ であり、$Px\mapsto Ax$ は代表元によらず $\Img P$ から $\Img A$ への線形等長同型を与える。同じ次元の両像の直交補で正規直交基底を選び、一方を他方に送って延長すれば、全空間のユニタリ作用素 $U$ が得られ、$A=UP$ である。逆にそのような分解では $\adj A A=P^2$ だから $P$ は一意である。$A$ が正則なら $U=AP^{-1}$、特異なら直交補上で選んだ対応を符号反転して別の $U$ を得る。
:::

複素数 $z=e^{i\theta}r$（$r\ge0$）の極形式との類比で、$U$ は位相、$P$ は絶対値に当たる。
極分解は、任意の線形変換が「ある直交軸に沿った伸縮のあとに等長変換を施したもの」であることを述べている。
$P$ を対角化すれば特異値分解が得られる（[行列解析](/math/linear-algebra/matrix-analysis/) で扱う）。

:::example[伸縮してから回す]
入力の軸を別々の倍率で伸ばしてから回す例を考える。

$$
A=\begin{pmatrix}0&-2\\3&0\end{pmatrix},\qquad \adj{A}A=\begin{pmatrix}9&0\\0&4\end{pmatrix}.
$$

$P=\diag(3,2)$ だから、$U=AP^{-1}=\begin{pmatrix}0&-1\\1&0\end{pmatrix}$ は角 $\pi/2$ の回転行列である。
特異な例 $A=\diag(2,0)$ では $P=A$ で、$U=I$ と $U=\diag(1,-1)$ のどちらでも同じ積 $UP=A$ を得る。
:::

### ブロック消去で正定値性を調べる

$H$ が大きくても、その一部の正定値性が既知なら、残りの小さい問題に帰着できる。

:::proposition[Schur 補行列]{#prop-linalg-schur-complement}
Hermite 行列 $H=\begin{pmatrix}B&C\\\adj{C}&D\end{pmatrix}$ において $B>0$ とし、$S=D-\adj{C}B^{-1}C$ とおく。
$S$ を $B$ に関する :term[Schur 補行列]{en="Schur complement" yomi="しゅーるほぎょうれつ"} という。このとき

1. $H$ は $\diag(B,S)$ と合同であり、$\det H=\det B\det S$ である。
2. $H>0\iff S>0$、$H\ge0\iff S\ge0$ である。
3. $H$ の慣性は $B$ の慣性 $(\dim B,0,0)$ と $S$ の慣性の和である。
:::

:::proof
$E=\begin{pmatrix}I&B^{-1}C\\0&I\end{pmatrix}$ とする。$\adj{(B^{-1})}=B^{-1}$ なので、ブロック積を取ると
$$
\adj E\diag(B,S)E
=\begin{pmatrix}B&C\\\adj C&\adj C B^{-1} C+S\end{pmatrix}=H.
$$
$\det E=1$ だから行列式を取れば積の公式を得る。:ref[thm-linalg-inertia] より合同で慣性は変わらない。各ブロックを標準形にする基底を並べると慣性は成分ごとの和になり、$B$ は正の方向だけを加える。したがって正定値性および半正定値性は $S$ の対応する条件とそれぞれ同値である。
:::

:::example[三変数の形式から二変数を消す]
$B=\diag(2,3),\ C=\transp{(2,3)},\ D=(t)$ とする。
$B^{-1}C=\transp{(1,1)}$ より $S=t-\transp{C}B^{-1}C=t-5$ であり、$\det H=6(t-5)$ である。
$H$ の形式は

$$
2(x+z)^2+3(y+z)^2+(t-5)z^2
$$

となる。
$t>5$ なら正定値、$t=5$ なら慣性 $(2,0,1)$ で根基は $(-1,-1,1)$ の張る直線、$t<5$ なら慣性 $(2,1,0)$ である。
固有値をすべて計算せずに慣性までわかる。
:::

### 二つの形式の同時対角化

$H>0$ と Hermite 行列 $B$ に対し、$Bx=\lambda Hx$ を満たす $\lambda$ と $x\ne0$ を求める問題を一般化固有値問題という。
力学の微小振動では、$H$ が運動エネルギー、$B$ がポテンシャルエネルギーの二次形式に当たる（[微小振動](/physics/analytical-mechanics/small-oscillations/) を参照）。

:::theorem[二つの Hermite 形式の同時対角化]{#thm-linalg-simultaneous-forms}
$H,B\in M_n(K)$ を Hermite 行列（実対称行列）とし、$H>0$ とする。
正則行列 $S$ で

$$
\adj{S}HS=I,\qquad \adj{S}BS=\diag(\mu_1,\ldots,\mu_n),\qquad \mu_j\in\R
$$

となるものが存在する。
$\mu_1,\ldots,\mu_n$ は方程式 $\det(tH-B)=0$ の根（重複度込み）であり、$S$ の第 $j$ 列 $s_j$ は $Bs_j=\mu_jHs_j$ を満たす。
:::

:::proof
Cholesky 分解 $H=\adj{R}R$（:ref[thm-linalg-positive-criteria]）をとる。
$C=\adj{(R^{-1})}BR^{-1}$ は Hermite 行列なので、ユニタリ行列 $Q$ で $\adj{Q}CQ=D=\diag(\mu_1,\ldots,\mu_n)$、$\mu_j\in\R$ とできる。
$S=R^{-1}Q$ とおけば

$$
\adj{S}HS=\adj{Q}\adj{(R^{-1})}\adj{R}RR^{-1}Q=I,\qquad \adj{S}BS=\adj{Q}CQ=D
$$

である。
$\adj{S}HS=I$ より $\adj{S}=S^{-1}H^{-1}$ であり、$\adj{S}BS=D$ に左から $HS$ を掛けると $BS=HSD$、すなわち $Bs_j=\mu_jHs_j$ である。
また $tH-B=\adj{(S^{-1})}(tI-D)S^{-1}$ なので、$K[t]$ 上の行列式の乗法性（:ref[prop-linalg-det-ring]）により

$$
\det(tH-B)=\abs{\det S}^{-2}\prod_{j=1}^n(t-\mu_j)
$$

である。
:::

$\adj{s_i}Hs_j=\delta_{ij}$ なので、一般化固有ベクトルは標準内積ではなく $H$ による内積 $\inner{x}{y}_H=\adj{y}Hx$ について正規直交である。
二つの形式の可換性は仮定していないことに注意する。
一方、どちらの形式も定値でない場合には同時対角化はできないことがある。
例えば $H=\diag(1,-1)$、$B=\begin{pmatrix}0&1\\1&0\end{pmatrix}$ に対し、正則な $S$ で $\adj{S}HS$ と $\adj{S}BS$ がともに実対角行列 $D_1,D_2$ になったとすると、$\det(tH-B)=\abs{\det S}^{-2}\det(tD_1-D_2)$ は実係数の一次以下の式の積である。
ところが $\det(tH-B)=-t^2-1$ は実根をもたない 2 次式なので、矛盾する。

:::example[通常の内積では直交しない固有方向]
$H=\diag(4,1),\ B=\begin{pmatrix}0&2\\2&0\end{pmatrix}$ とする。
$H^{1/2}=\diag(2,1)$ であり、$H^{-1/2}BH^{-1/2}=\begin{pmatrix}0&1\\1&0\end{pmatrix}$ の固有値は $1,-1$ である。
固有ベクトル $(1,1),(1,-1)$ を $H^{-1/2}$ で戻すと、元の方向は $(1,2),(1,-2)$ に比例する。
標準内積では $1-4=-3$ で直交しないが、$H$ による内積では $4-4=0$ である。
どちらも $H$ での長さは $\sqrt8$ なので、これで割れば $S$ の列が得られる。
この計算は Cholesky 因子 $R$ の代わりに平方根 $H^{1/2}$ を用いても同様に行える。
:::

### 極分解における左右の順序

$A=UP$ の $P$ と $U$ の順を入れ替えるだけでは一般に $A$ に戻らない。
$Q=UP\adj{U}$ とおくと、$Q$ は $P$ と同じ固有値をもつので半正定値であり、$Q^2=U\adj{A}A\adj{U}=A\adj{A}$、$QU=UP=A$ である。
平方根の一意性から $Q=(A\adj{A})^{1/2}$ であり、$A=QU$ を左極分解という。
$P=Q$ は $\adj{A}A=A\adj{A}$、すなわち $A$ の正規性と同値である。
例えば $A=\begin{pmatrix}0&-2\\3&0\end{pmatrix}$ では $P=\diag(3,2)$、$Q=\diag(2,3)$ であり、$UPe_1=QUe_1=(0,3)$ だが $PUe_1=(0,2)$ となる。

## 演習問題

::::exercise{level=basic}
$H=\begin{pmatrix}0&-i\\i&0\end{pmatrix}$ の固有ベクトルからなる正規直交基底を求め、スペクトル分解を書け。

:::solution
$H(x,y)=(-iy,ix)$ である。
固有値 $1$ では $-iy=x$、すなわち $y=ix$ であり、固有値 $-1$ では $y=-ix$ である。
したがって $q_+=(1,i)/\sqrt2,\ q_-=(1,-i)/\sqrt2$ が求める基底である。
実際 $\inner{q_+}{q_-}=\frac12(1\cdot1+i\cdot\conj{(-i)})=\frac12(1+i\cdot i)=0$ である。
スペクトル分解は

$$
H=q_+\adj{q_+}-q_-\adj{q_-},\qquad q_\pm\adj{q_\pm}=\frac12\begin{pmatrix}1&\mp i\\\pm i&1\end{pmatrix}=\frac{I\pm H}2
$$

であり、本文の例と一致する。
:::
::::

::::exercise{level=basic}
次の実対称行列の慣性と根基を求めよ。

$$
A=\begin{pmatrix}4&2\\2&1\end{pmatrix}
$$

:::solution
形式は $4x^2+4xy+y^2=(2x+y)^2$ である。
$u=2x+y,\ v=x$ は正則な座標変換（行列式 $-1$）なので、形式は $u^2+0\cdot v^2$ となり、慣性は $(1,0,1)$ である。
根基は $A\transp{(x,y)}=0$ の解、すなわち $2x+y=0$ で、$(1,-2)$ の張る直線である。
固有値で確かめると $\tr A=5$、$\det A=0$ より固有値は $5,0$ であり、正が一つ、零が一つである。
:::
::::

::::exercise{level=basic}
$\begin{pmatrix}4&2\\2&5\end{pmatrix}$ の Cholesky 分解を求めよ。

:::solution
$R=\begin{pmatrix}r_{11}&r_{12}\\0&r_{22}\end{pmatrix}$ とおくと $\transp{R}R=\begin{pmatrix}r_{11}^2&r_{11}r_{12}\\r_{11}r_{12}&r_{12}^2+r_{22}^2\end{pmatrix}$ である。
$r_{11}^2=4$ と $r_{11}>0$ から $r_{11}=2$、$2r_{12}=2$ から $r_{12}=1$ である。
残りは $r_{22}^2=5-1=4$ なので、$r_{22}>0$ を選べば

$$
R=\begin{pmatrix}2&1\\0&2\end{pmatrix}
$$

となる。
首座小行列式 $4,\ 16$ はともに正であり、:ref[thm-linalg-positive-criteria] と整合する。
:::
::::

::::exercise{level=basic}
次の正定値行列の半正定値平方根を求めよ。

$$
H=\begin{pmatrix}5&4\\4&5\end{pmatrix}
$$

:::solution
$H(1,1)=9(1,1)$、$H(1,-1)=(1,-1)$ なので、固有値は $9,1$ であり、スペクトル射影は $P_9=\frac12\begin{pmatrix}1&1\\1&1\end{pmatrix}$、$P_1=\frac12\begin{pmatrix}1&-1\\-1&1\end{pmatrix}$ である。
したがって

$$
H^{1/2}=3P_9+1\cdot P_1=\frac12\begin{pmatrix}3+1&3-1\\3-1&3+1\end{pmatrix}=\begin{pmatrix}2&1\\1&2\end{pmatrix}
$$

となる。
$\begin{pmatrix}2&1\\1&2\end{pmatrix}^2=\begin{pmatrix}5&4\\4&5\end{pmatrix}$ であり、固有値 $3,1$ は正なので、一意性によりこれが求める平方根である。
:::
::::

::::exercise{level=basic}
二次曲線 $3x^2+2xy+3y^2=8$ を主軸変換で標準形に直し、曲線の種類と主軸の方向・半軸の長さを求めよ。

:::solution
$A=\begin{pmatrix}3&1\\1&3\end{pmatrix}$ の固有値は $4$（固有ベクトル $(1,1)$）と $2$（固有ベクトル $(1,-1)$）である。
$u=(x+y)/\sqrt2,\ v=(x-y)/\sqrt2$ とすれば $3x^2+2xy+3y^2=4u^2+2v^2$ なので、方程式は

$$
\frac{u^2}{2}+\frac{v^2}{4}=1
$$

となる。
二つの固有値がともに正で右辺が $8>0$ なので楕円であり、$(1,1)/\sqrt2$ 方向の半軸は $\sqrt2$、$(1,-1)/\sqrt2$ 方向の半軸は $2$ である。
:::
::::

::::exercise{level=basic}
次の行列が $\SO(3)$ に属することを確かめ、回転軸と回転角を求めよ。

$$
A=\frac13\begin{pmatrix}-1&2&2\\2&-1&2\\2&2&-1\end{pmatrix}
$$

:::solution
$A=\frac23J-I$（$J$ は全成分 $1$ の行列）と書ける。
$J^2=3J$ より $A^2=\frac49\cdot3J-\frac43J+I=I$ であり、$A$ は対称なので $\transp{A}A=A^2=I$ である。
$J$ の固有値は $3$（方向 $(1,1,1)$）と $0$（重複度 $2$）なので、$A$ の固有値は $1$ と $-1$（重複度 $2$）であり、$\det A=1$ である。
よって $A\in\SO(3)$ である。
回転軸は $\Ker(A-I)=\spn\{(1,1,1)\}$ である。
$\tr A=-1=1+2\cos\theta$ より $\cos\theta=-1$、$\theta=\pi$ であり、$A$ は軸 $(1,1,1)$ のまわりの半回転である。
:::
::::

::::exercise{level=basic}
$H=\begin{pmatrix}2&1&0\\1&-1&1\\0&1&3\end{pmatrix}$ の慣性を、固有値を求めずに決定せよ。

:::solution
首座小行列式は

$$
\Delta_1=2,\qquad \Delta_2=-2-1=-3,\qquad \Delta_3=2(-3-1)-1\cdot(3-0)=-11
$$

であり、いずれも零でない。
数列 $1,2,-3,-11$ の符号の変化は $2\to-3$ の 1 回だけなので、:ref[thm-linalg-jacobi-signature] より慣性は $(2,1,0)$ である。
:::
::::

::::exercise{level=advanced}
任意の複素正方行列 $A$ は Hermite 行列 $H,K$ により $A=H+iK$ と一意に書けることを示せ。
さらに、$A$ が正規であることと $HK=KH$ が同値であることを示せ。

:::solution
$A=H+iK$ と書けたとすると、随伴をとって $\adj{A}=H-iK$ である。
二つの式を加減すれば

$$
H=\frac{A+\adj{A}}2,\qquad K=\frac{A-\adj{A}}{2i}
$$

が強制されるので一意である。
逆にこの式で $H,K$ を定めれば、どちらも自己随伴であり（$\adj{K}=\frac{\adj{A}-A}{-2i}=K$）、$H+iK=A$ となる。
積を展開すると

$$
\adj{A}A-A\adj{A}=(H-iK)(H+iK)-(H+iK)(H-iK)=2i(HK-KH)
$$

なので、$\adj{A}A=A\adj{A}$ と $HK=KH$ は同値である。
:::
::::

::::exercise{level=advanced}
$H=\begin{pmatrix}1&1\\1&2\end{pmatrix}$、$B=\begin{pmatrix}1&1\\1&0\end{pmatrix}$ について、$\transp{S}HS=I$ かつ $\transp{S}BS$ が対角行列となる正則行列 $S$ を一つ求めよ。

:::solution
$\Delta_1=1,\ \Delta_2=1$ より $H>0$ なので :ref[thm-linalg-simultaneous-forms] が使える。
$\det(tH-B)=\det\begin{pmatrix}t-1&t-1\\t-1&2t\end{pmatrix}=(t-1)(2t-(t-1))=(t-1)(t+1)$ であり、$\mu=1,-1$ である。
$\mu=1$：$(B-H)s=0$ は $\begin{pmatrix}0&0\\0&-2\end{pmatrix}s=0$ なので $s=(1,0)$ であり、$\transp{s}Hs=1$ である。
$\mu=-1$：$(B+H)s=0$ は $\begin{pmatrix}2&2\\2&2\end{pmatrix}s=0$ なので $s=(1,-1)$ であり、$\transp{s}Hs=1-2+2=1$ である。
したがって

$$
S=\begin{pmatrix}1&1\\0&-1\end{pmatrix},\qquad \transp{S}HS=I,\qquad \transp{S}BS=\diag(1,-1)
$$

である。
二つの列は $H$ による内積で直交するが、標準内積では $(1,0)\cdot(1,-1)=1\ne0$ である。
:::
::::

::::exercise{level=advanced}
$K=\C$ とし、$T\in\End(V)$ のすべての不変部分空間 $W$ について $W^\perp$ も $T$ 不変であるとする。
このとき $T$ は正規であることを示せ（:ref[prop-linalg-normal-invariant] の逆）。

:::solution
$\dim V$ に関する帰納法で、$V$ が $T$ の固有ベクトルからなる正規直交基底をもつことを示す。
$\dim V=1$ なら明らかである。
$K=\C$ なので $T$ は単位固有ベクトル $v$ をもち、$\spn\{v\}$ は $T$ 不変だから、仮定により $U=v^\perp$ も $T$ 不変である。
$T|_U$ も同じ性質をもつことを確かめる。
$W\subset U$ が $T|_U$ 不変なら $W$ は $T$ 不変なので、$W^\perp$（$V$ での直交補空間）は $T$ 不変であり、$U$ における $W$ の直交補空間 $W^\perp\cap U$ も二つの不変部分空間の共通部分として $T$ 不変である。
帰納法の仮定から $U$ は $T|_U$ の固有ベクトルからなる正規直交基底をもち、それに $v$ を加えれば $V$ の正規直交基底になる。
よって :ref[thm-linalg-spectral] により $T$ は正規である。
:::
::::

::::exercise{level=advanced}
半正定値行列 $H$ について、$\adj{x}Hx=0$ と $Hx=0$ が同値であることを示せ。
また不定値の場合の反例を挙げよ。

:::solution
$P=H^{1/2}$（:ref[thm-linalg-polar]）は自己随伴なので $\adj{x}Hx=\adj{x}P^2x=\norm{Px}^2$ である。
これが零なら $Px=0$ であり、$Hx=P(Px)=0$ である。
逆向きは $Hx=0$ なら $\adj{x}Hx=0$ となることから明らかである。
不定値の $H=\diag(1,-1)$、$x=(1,1)$ では $\adj{x}Hx=1-1=0$ だが $Hx=(1,-1)\ne0$ である。
この結果から、半正定値形式では二次形式の零集合と根基が一致する。
:::
::::

::::exercise{level=advanced}
Hermite 行列 $A,B$ について、ユニタリ同値であることと固有値の多重集合が等しいことが同値であることを示せ。
また $A=\diag(1,4)$、$B=\diag(1,9)$ を用いて、ユニタリ同値と合同の違いを説明せよ。

:::solution
固有値を同じ順に並べた対角行列 $D$ によって $A=UD\adj{U}$、$B=WD\adj{W}$（$U,W$ ユニタリ）とできるなら、$S=U\adj{W}$ はユニタリで

$$
\adj{S}AS=W\adj{U}AU\adj{W}=WD\adj{W}=B
$$

である。
逆にユニタリ同値なら相似なので特性多項式が等しく（:ref[prop-linalg-charpoly-similarity]）、固有値の多重集合は等しい。

例の $A,B$ は固有値 $4$ と $9$ が異なるのでユニタリ同値ではない。
しかし慣性はともに $(2,0,0)$ なので合同である（:ref[cor-linalg-congruence-classes]）。
実際、$C=\diag(1,3/2)$ とすれば $\transp{C}AC=\diag(1,9)=B$ である。
$C$ は第 2 座標を伸縮するので、合同は符号の個数を保ちつつ、単位ベクトル上の形式の値を変えてしまう。
:::
::::

::::exercise{level=advanced}
写像 $f:\R^n\to\R^n$ が $f(0)=0$ とすべての $x,y$ について $\norm{f(x)-f(y)}=\norm{x-y}$ を満たすとする。
$f$ は線形であり、ある $A\in\mathrm{O}(n)$ により $f(x)=Ax$ と書けることを示せ。

:::solution
$y=0$ とすれば $\norm{f(x)}=\norm{x}$ である。
$\norm{f(x)-f(y)}^2=\norm{f(x)}^2-2\inner{f(x)}{f(y)}+\norm{f(y)}^2$ と $\norm{x-y}^2$ を比べると、$\inner{f(x)}{f(y)}=\inner{x}{y}$ がすべての $x,y$ で成り立つ。
特に $a_j=f(e_j)$ は正規直交基底をなす。
任意の $x$ について $f(x)=\sum_j\inner{f(x)}{a_j}a_j=\sum_j\inner{x}{e_j}a_j=\sum_jx_ja_j$ である（:ref[prop-linalg-orthonormal-expansion]）。
よって $A=(a_1\ \cdots\ a_n)$ とすれば $f(x)=Ax$ であり、$A$ の列は正規直交基底なので $A\in\mathrm{O}(n)$ である（:ref[prop-linalg-unitary-equiv]）。
:::
::::

::::exercise{level=advanced}
$K=\C$ とし、$A\in M_n(\C)$ を正規行列とする。

1. $A$ が冪零（ある $m$ で $A^m=O$）ならば $A=O$ であることを示せ。
2. $AB=BA$ を満たす任意の $B\in M_n(\C)$ について $\adj{A}B=B\adj{A}$ であることを示せ。
3. 2 は $A$ が正規でないと成り立たないことを例で示せ。

:::solution
1. スペクトル分解 $A=\sum_j\lambda_jP_j$ より $A^m=\sum_j\lambda_j^mP_j=O$ である。
$E_{\lambda_j}(A)$ の元に作用させれば $\lambda_j^m=0$、すなわち $\lambda_j=0$ なので、$A=O$ である。

2. :ref[prop-linalg-spectral-projections] の 3 により $\adj{A}=q(A)$ となる多項式 $q$ がある。
$B$ は $A$ と可換なので $A$ の冪とも可換であり、その一次結合 $q(A)=\adj{A}$ とも可換である。

3. $A=\begin{pmatrix}0&1\\0&0\end{pmatrix}$、$B=A$ とすれば $AB=BA$ だが、$\adj{A}B=\diag(0,1)$、$B\adj{A}=\diag(1,0)$ で異なる。
:::
::::

## 確認テスト

::::quiz
複素正規行列で常に成り立つ性質をすべて選べ。

- [ ] すべての固有値が実数である
- [x] 異なる固有値に属する固有ベクトルは直交する
- [ ] すべての固有値の絶対値が $1$ である
- [x] $\Ker A=\Ker\adj{A}$ である

:::explanation
正解は (B) と (D) である。
:ref[prop-linalg-normal-basic] の 3 により異なる固有値の固有ベクトルは直交し、2 により $\Ker A=E_0(A)=E_0(\adj{A})=\Ker\adj{A}$ である。
一方、スカラー行列 $2iI$ は正規だが、固有値 $2i$ は実数でも絶対値 $1$ でもないので (A) と (C) は誤りである。
:::
::::

::::quiz{answer="2" tol="1e-6"}
$q(x,y,z)=x^2-y^2+3z^2$ の正の慣性指数を求めよ。

:::explanation
$z$ を $z'=\sqrt3z$ に取り直すと $q=x^2-y^2+z'^2$ であり、正の平方が二つ、負の平方が一つである。
慣性は $(2,1,0)$ なので、正の慣性指数は $2$ である。
:::
::::

::::quiz
Hermite 行列 $H$ が半正定値であるための必要十分条件を選べ。

- [ ] 首座小行列式がすべて非負である
- [x] 主小行列式がすべて非負である
- [ ] 行列式が非負である

:::explanation
正解は (B) である（:ref[prop-linalg-psd-criteria]）。
全主小行列式が非負なら :ref[lem-linalg-principal-minor-expansion] により $H+\eps I$ の首座小行列式はすべて正なので $H+\eps I>0$ であり、$\eps\to0$ として $H\ge0$ を得る。
逆は座標部分空間への制限から従う。
$\diag(0,-1)$ は首座小行列式 $0,0$ と行列式 $0$ がすべて非負だが、$e_2$ で形式の値が $-1$ なので半正定値でない。
よって (A) と (C) は十分条件でない。
:::
::::

::::quiz{answer="3" tol="1e-6"}
$\begin{pmatrix}2&1\\1&2\end{pmatrix}$ の Rayleigh 商 $\inner{Ax}{x}/\norm{x}^2$（$x\ne0$）の最大値を求めよ。

:::explanation
固有値は $3$（方向 $(1,1)$）と $1$（方向 $(1,-1)$）であり、:ref[prop-linalg-rayleigh] により Rayleigh 商は固有値の重み付き平均である。
最大値は最大固有値 $3$ で、$x$ が $(1,1)$ の倍数のときに達成される。
:::
::::

::::quiz{answer="120" tol="1e-6"}
$A\in\SO(3)$ が $\tr A=0$ を満たすとき、回転角 $\theta$（$0\le\theta\le\pi$）を度数法で求めよ。

:::explanation
:ref[thm-linalg-rotation-axis] により $\tr A=1+2\cos\theta$ である。
$1+2\cos\theta=0$ より $\cos\theta=-1/2$ であり、$0\le\theta\le\pi$ では $\theta=2\pi/3$、すなわち $120^\circ$ である。
例えば座標の巡回置換行列はこの回転で、軸は $(1,1,1)$ である。
:::
::::

::::quiz
実対称行列 $A$ について常に成り立つものをすべて選べ。

- [x] 固有値はすべて実数である
- [x] 直交行列で対角化できる
- [ ] 対角成分がすべて正なら正定値である
- [ ] 正則である

:::explanation
正解は (A) と (B) である。
(A) は :ref[prop-linalg-selfadjoint-eigen]、(B) は :ref[thm-linalg-real-spectral] による。
(C) は誤りで、$\begin{pmatrix}1&2\\2&1\end{pmatrix}$ は対角成分が正だが $\Delta_2=-3<0$ であり、$(1,-1)$ で形式の値は $-2$ である。
(D) も誤りで、零行列は実対称だが正則でない。
:::
::::

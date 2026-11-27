---
title: 固有値と対角化
description: 固有値・固有空間と不変部分空間、特性多項式と二種類の重複度、対角化可能性の判定とその応用、Cayley–Hamilton の定理、三角化と同時三角化
notation:
  - symbol: 'E_\lambda(T)'
    meaning: 線形変換 T の固有値 λ に対する固有空間
  - symbol: '\Phi_A(t)'
    meaning: 正方行列 A の特性多項式 det(tI − A)
  - symbol: '\Phi_T(t)'
    meaning: 線形変換 T の特性多項式
  - symbol: 'm_\lambda'
    meaning: 固有値 λ の代数的重複度
  - symbol: 'g_\lambda'
    meaning: 固有値 λ の幾何学的重複度
  - symbol: 'T|_W'
    meaning: 不変部分空間 W への T の制限
  - symbol: 'p(A)'
    meaning: 多項式 p の変数に正方行列 A を代入した行列
  - symbol: '\diag(\lambda_1,\ldots,\lambda_n)'
    meaning: 対角成分が λ1, …, λn の対角行列
---

## この章の目標

- 固有値・固有ベクトル・固有空間と不変部分空間を定義し、特性多項式を用いて計算できるようになる。
- 固有値の代数的重複度と幾何学的重複度を区別し、対角化可能性を判定する。
- 対角化を用いて行列の冪や線形漸化式を計算する。
- Cayley–Hamilton の定理を証明し、行列の多項式の計算に用いる。
- 特性多項式が一次式の積に分解すれば三角化できることを示し、トレース・行列式と固有値の関係や、可換な行列の同時三角化を導く。

[線形写像と表現行列](/math/linear-algebra/linear-maps/) で見たように、線形変換の表現行列は基底を取り替えると相似な行列に移る（:ref[prop-linalg-similarity-equivalence]）。そこで、表現行列がなるべく簡単になる基底を探したい。最も理想的なのは、各基底ベクトルが自分自身のスカラー倍に移る基底であり、そのとき表現行列は対角行列になる。対角行列ならば冪 $A^k$ の計算は各対角成分の $k$ 乗で済み、漸化式 $x_{k+1}=Ax_k$ や微分方程式 $x'=Ax$ の解析も各座標ごとの 1 次元の問題に分解される。本章では、このような基底がいつ存在するかを調べる。

本章の定義と定理は、特に断らない限り任意の体 $K$ 上で成り立つ。代数学の基本定理を用いる箇所（$K=\C$）や、$K=\R$ または $\C$ を仮定する例ではその旨を記す。$V$ は $K$ 上の $n$ 次元ベクトル空間（$n\ge1$）とする。ただし固有値の定義自体は無限次元の空間でも意味をもつ。

## 固有値と固有ベクトル

### 定義

:::definition[固有値・固有ベクトル・固有空間]{#def-linalg-eigenvalue}
$V$ を $K$ 上のベクトル空間、$T:V\to V$ を線形変換とし、$\lambda\in K$ とする。$T(v)=\lambda v$ を満たす $v\ne0$ が存在するとき、$\lambda$ を $T$ の :term[固有値]{en="eigenvalue" yomi="こゆうち"} という。このとき、$T(v)=\lambda v$ を満たす $0$ でないベクトル $v$ を、固有値 $\lambda$ に属する :term[固有ベクトル]{en="eigenvector" yomi="こゆうべくとる"} という。固有値 $\lambda$ に対して

$$
E_\lambda(T)\defeq\Ker(T-\lambda\id_V)=\{v\in V\mid T(v)=\lambda v\}
$$

を $\lambda$ に対する :term[固有空間]{en="eigenspace" yomi="こゆうくうかん"} という。行列 $A\in M_n(K)$ については $T=L_A$ として同じ用語を用い、$E_\lambda(A)\defeq\Ker L_{A-\lambda I_n}=\{x\in K^n\mid Ax=\lambda x\}$ と書く。
:::

零ベクトルはどの $\lambda$ についても $T(0)=\lambda0$ を満たすので、固有値の定義では**零でない解の存在**が本質的である。一方、固有空間は線形写像 $T-\lambda\id_V$ の核なので部分空間であり、零ベクトルを含む。$E_\lambda(T)$ は固有値 $\lambda$ に属する固有ベクトル全体に $0$ を付け加えたものである。

$0$ も固有値になりうる。$0$ が $T$ の固有値であることは $\Ker T\ne\{0\}$、すなわち $T$ が単射でないことと同値である（:ref[prop-linalg-kernel-injective]）。有限次元ならば、これは $T$ が同型写像でないことと同値である（:ref[cor-linalg-equal-dimension-map]）。

:::proposition[固有値と表現行列]{#prop-linalg-eigen-coordinates}
$V$ を有限次元とし、$\mathcal B$ を $V$ の基底、$A=[T]_{\mathcal B}$ とする。$v\in V$ について、$T(v)=\lambda v$ であることと $A[v]_{\mathcal B}=\lambda[v]_{\mathcal B}$ であることは同値である。特に $T$ と $A$ の固有値は一致し、座標写像は $E_\lambda(T)$ から $E_\lambda(A)$ への同型写像を与える。
:::

:::proof
:eqref[eq-linalg-coordinate-action] より $[T(v)]_{\mathcal B}=A[v]_{\mathcal B}$ であり、座標写像 $\varphi_{\mathcal B}$ は線形かつ単射なので、$T(v)=\lambda v\iff[T(v)]_{\mathcal B}=\lambda[v]_{\mathcal B}\iff A[v]_{\mathcal B}=\lambda[v]_{\mathcal B}$ である。よって $\varphi_{\mathcal B}(E_\lambda(T))=E_\lambda(A)$ であり、$\varphi_{\mathcal B}$ は同型写像なので（:ref[prop-linalg-coordinate-isomorphism]）その制限も同型写像である。$v\ne0\iff[v]_{\mathcal B}\ne0$ なので固有値も一致する。
:::

したがって抽象的な空間の線形変換でも、基底を選べば行列の固有値問題、すなわち斉次方程式 $(A-\lambda I_n)x=0$ が自明でない解をもつかという問題に帰着する。

:::example[2 つの固有方向]{#exm-linalg-eigen-2x2}
$A=\begin{pmatrix}2&1\\1&2\end{pmatrix}$ について

$$
A\begin{pmatrix}1\\1\end{pmatrix}=3\begin{pmatrix}1\\1\end{pmatrix},
\qquad
A\begin{pmatrix}1\\-1\end{pmatrix}=\begin{pmatrix}1\\-1\end{pmatrix}
$$

なので、$3$ と $1$ は $A$ の固有値であり、$\transp{(1,1)}$, $\transp{(1,-1)}$ はそれぞれに属する固有ベクトルである。$L_A$ は直線 $y=x$ 上では $3$ 倍、直線 $y=-x$ 上では $1$ 倍として作用する。一方、例えば $v=\transp{(1,0)}$ はこれらの直線上になく、$Av=\transp{(2,1)}$ は $v$ のスカラー倍でない。
:::

::interactive[2 次正方行列による単位円の像と固有方向]{name="matrix-action"}

初期値は上の $A$ である。単位円上の点 $v$ を動かすと、その像 $Av$ は楕円上を動く。$v$ を直線 $y=x$ 上に置くと $Av$ は同じ直線上で原点から $3$ 倍の距離にあり、$y=-x$ 上に置くと $v$ と $Av$ は重なる。また $a=d=0$, $b=-1$, $c=1$ とすると $Av=\transp{(-y,x)}$ となり（$90^\circ$ 回転）、実平面ではどの $0$ でないベクトルも自分と同じ直線上に移らない。すなわち実固有値は存在しない（:ref[exm-linalg-rotation-eigen]）。

:::example[微分作用素の固有値]{#exm-linalg-eigen-derivative}
1. $V=C^\infty(\R)$（何回でも微分可能な実数値関数全体、実ベクトル空間）上の微分 $D:y\mapsto y'$ を考える。任意の $\lambda\in\R$ について $D(e^{\lambda x})=\lambda e^{\lambda x}$ なので、すべての実数が $D$ の固有値である。さらに $y'=\lambda y$ ならば $(ye^{-\lambda x})'=(y'-\lambda y)e^{-\lambda x}=0$ なので $ye^{-\lambda x}$ は定数であり、$E_\lambda(D)=\langle e^{\lambda x}\rangle$ である。無限次元空間上の線形変換は、無限個の固有値をもちうる。
2. 一方、$K=\R$ または $\C$ とし、$D$ を $K[t]_{\le n}$ 上の微分とすると、$D$ の固有値は $0$ だけである。実際、$p\ne0$ で $Dp=\lambda p$ ならば、$\deg Dp<\deg p$（$\deg p\ge1$ のとき）または $Dp=0$（$p$ が定数のとき）なので、$\lambda\ne0$ では両辺の次数が合わない。$E_0(D)$ は定数多項式全体である。
3. :ref[exm-linalg-projection] の射影 $p:V\to V$（$V=W_1\oplus W_2$, $W_1,W_2\ne\{0\}$）の固有値は $1$ と $0$ であり、$E_1(p)=W_1$, $E_0(p)=W_2$ である。実際 $p(v)=\lambda v$ ならば $\lambda v=p(v)=p^2(v)=\lambda^2v$ より $\lambda=0$ または $1$ である。
:::

### 異なる固有値に属する固有ベクトル

:::theorem[異なる固有値に属する固有ベクトルの一次独立性]{#thm-linalg-eigenvectors-independent}
$\lambda_1,\ldots,\lambda_r$ を $T$ の相異なる固有値とし、それぞれに属する固有ベクトル $v_i\in E_{\lambda_i}(T)$ を 1 つずつ選ぶ。このとき $v_1,\ldots,v_r$ は一次独立である。
:::

:::proof
$r$ に関する帰納法で示す。$r=1$ のときは $v_1\ne0$ より一次独立である。$r-1$ 個のとき主張が成り立つとし、$\sum_{i=1}^r c_i v_i=0$ とする。両辺に $T-\lambda_r\id_V$ を施すと、$(T-\lambda_r\id_V)v_i=(\lambda_i-\lambda_r)v_i$ なので第 $r$ 項が消え、

$$
\sum_{i=1}^{r-1}c_i(\lambda_i-\lambda_r)v_i=0
$$

を得る。帰納法の仮定より各係数 $c_i(\lambda_i-\lambda_r)$ は $0$ であり、$\lambda_i\ne\lambda_r$ なので $c_i=0$（$i<r$）である。元の式は $c_rv_r=0$ となり、$v_r\ne0$ より $c_r=0$ である。
:::

:::proposition[固有空間の和は直和]{#prop-linalg-eigenspaces-direct-sum}
$\lambda_1,\ldots,\lambda_r$ を $T$ の相異なる固有値とすると、$E_{\lambda_1}(T)+\cdots+E_{\lambda_r}(T)$ は直和である。特に $V$ が $n$ 次元ならば $\sum_i\dim E_{\lambda_i}(T)\le n$ であり、$T$ の相異なる固有値は高々 $n$ 個である。
:::

:::proof
$u_i\in E_{\lambda_i}(T)$ が $u_1+\cdots+u_r=0$ を満たすとき、すべての $u_i$ が $0$ であることを示せばよい（表示の一意性は差をとれば零の表示の一意性に帰着する。:ref[def-linalg-direct-sum]）。$0$ でない $u_i$ があったとすると、$0$ でないものだけを集めた和が $0$ となるが、それらは相異なる固有値に属する固有ベクトルなので、:ref[thm-linalg-eigenvectors-independent] に反する。後半は :ref[prop-linalg-direct-sum-criterion] の 2 と :ref[cor-linalg-subspace-dimension] による。
:::

### 不変部分空間

固有ベクトル $v$ の張る直線 $\langle v\rangle$ は、$T$ によってそれ自身の中に移される。この性質を一般化する。

:::definition[不変部分空間]{#def-linalg-invariant-subspace}
$T\in\End(V)$ とする。$V$ の部分空間 $W$ が $T(W)\subset W$ を満たすとき、$W$ を $T$ の :term[不変部分空間]{en="invariant subspace" yomi="ふへんぶぶんくうかん"}（$T$ 不変部分空間）という。このとき $T$ の $W$ への制限 $T|_W:W\to W$ は $W$ 上の線形変換である。
:::

$\{0\}$, $V$, $\Ker T$, $\Img T$, 各固有空間 $E_\lambda(T)$ は $T$ 不変である（例えば $v\in E_\lambda(T)$ ならば $T(v)=\lambda v\in E_\lambda(T)$）。また $1$ 次元の部分空間 $\langle v\rangle$ が $T$ 不変であることと、$v$ が $T$ の固有ベクトルであることは同値である。より一般に、$S\in\End(V)$ が $T$ と可換（$ST=TS$）ならば、$\Ker S$, $\Img S$ および $S$ の固有空間は $T$ 不変である。例えば $v\in E_\mu(S)$ ならば $S(T(v))=T(S(v))=\mu T(v)$ より $T(v)\in E_\mu(S)$ である。

:::proposition[不変部分空間とブロック三角行列]{#prop-linalg-invariant-block}
$W$ を $T$ 不変部分空間とし、$W$ の基底 $(w_1,\dots,w_k)$ を $V$ の基底 $\mathcal B=(w_1,\dots,w_k,v_{k+1},\dots,v_n)$ に延長する。このとき

$$
[T]_{\mathcal B}=\begin{pmatrix}A_1&H\\O&A_2\end{pmatrix},\qquad A_1=[T|_W]_{(w_1,\dots,w_k)}\in M_k(K),\ A_2\in M_{n-k}(K)
$$

の形になる。逆に、ある基底でこの形になれば、最初の $k$ 個の基底ベクトルの張る部分空間は $T$ 不変である。さらに $V=W_1\oplus W_2$ で $W_1,W_2$ がともに $T$ 不変ならば、それぞれの基底を並べた基底で $[T]$ はブロック対角行列 $\begin{pmatrix}A_1&O\\O&A_2\end{pmatrix}$ になる。
:::

:::proof
$j\le k$ ならば $T(w_j)\in W=\langle w_1,\dots,w_k\rangle$ なので、$[T]_{\mathcal B}$ の第 $j$ 列の第 $k+1$ 成分以下は $0$ であり、上の $k$ 成分は $T(w_j)=T|_W(w_j)$ を $(w_1,\dots,w_k)$ で表した座標である。逆に第 $j$ 列（$j\le k$）の第 $k+1$ 成分以下が $0$ ならば $T(w_j)\in W$ であり、:ref[prop-linalg-image-preimage] より $T(W)=\langle T(w_1),\dots,T(w_k)\rangle\subset W$ である。直和の場合は、$W_2$ の基底ベクトルの像が $W_2$ に属するので右上のブロックも $O$ になる。
:::

## 特性多項式

### 多項式についての準備

固有値を求める方程式は多項式の方程式になる。そこで、以下で用いる多項式の基本性質をまとめて証明しておく。$K[t]$ の元 $p=\sum_{k=0}^dc_kt^k$（$c_d\ne0$）について、$d$ を $p$ の**次数**（$\deg p$）、$c_d$ を**最高次係数**といい、最高次係数が $1$ の多項式を**モニック**という。零多項式の次数は定めない。多項式の積は $(\sum_ia_it^i)(\sum_jb_jt^j)\defeq\sum_k(\sum_{i+j=k}a_ib_j)t^k$ で定義され、$K[t]$ は可換環である（:ref[prop-linalg-det-ring] の直前の説明）。多項式 $p$ は記号の列（係数の列）であって、関数 $K\to K$ とは区別する。

:::lemma[行列の多項式]{#lem-linalg-polynomial-substitution}
$A\in M_n(K)$ と $p=\sum_kc_kt^k\in K[t]$ に対し、$A^0\defeq I_n$ として $p(A)\defeq\sum_kc_kA^k\in M_n(K)$ と定める（$n=1$ のときは、スカラー $c\in K$ に対する値 $p(c)$ である）。$p,q\in K[t]$, $c\in K$ について次が成り立つ。

1. $(p+q)(A)=p(A)+q(A)$, $(cp)(A)=c\,p(A)$, $(pq)(A)=p(A)q(A)$。特に $p(A)$ と $q(A)$ は可換である。
2. $Av=\lambda v$ ならば $p(A)v=p(\lambda)v$。
3. $P$ が正則ならば $p(P^{-1}AP)=P^{-1}p(A)P$。
4. $T\in\End(V)$ についても $T^0\defeq\id_V$ として $p(T)\defeq\sum_kc_kT^k$ と定めると、1 が同様に成り立ち、任意の基底 $\mathcal B$ について $[p(T)]_{\mathcal B}=p([T]_{\mathcal B})$ である。
:::

:::proof
1. 和とスカラー倍は定義から明らかである。$p=\sum_ia_it^i$, $q=\sum_jb_jt^j$ とすると、行列の積の分配法則と $A^iA^j=A^{i+j}$（結合法則）より

$$
p(A)q(A)=\sum_{i,j}a_ib_jA^{i+j}=\sum_k\Bigl(\sum_{i+j=k}a_ib_j\Bigr)A^k=(pq)(A)
$$

である。$pq=qp$ なので $p(A)q(A)=(pq)(A)=(qp)(A)=q(A)p(A)$。

2. $A^kv=\lambda^kv$ が $k$ に関する帰納法で成り立つ（$A^{k+1}v=A(\lambda^kv)=\lambda^{k+1}v$）ので、$p(A)v=\sum_kc_k\lambda^kv=p(\lambda)v$。

3. :ref[prop-linalg-similarity-invariants] の 3 より $(P^{-1}AP)^k=P^{-1}A^kP$ なので、一次結合をとればよい。

4. 1 の証明は $\End(V)$ でもそのまま成り立つ。:ref[cor-linalg-invertible-map-matrix] より $[T^k]_{\mathcal B}=[T]_{\mathcal B}^k$ であり、:ref[prop-linalg-hom-matrix-iso] より $T\mapsto[T]_{\mathcal B}$ は線形なので、$[p(T)]_{\mathcal B}=\sum_kc_k[T]_{\mathcal B}^k=p([T]_{\mathcal B})$。
:::

:::lemma[多項式の割り算]{#lem-linalg-polynomial-division}
1. $p,q\in K[t]$ が零でなければ $pq\ne0$ であり、$\deg(pq)=\deg p+\deg q$ である。したがって $p\ne0$ かつ $pq=pr$ ならば $q=r$（約分できる）。
2. $g\ne0$ とする。任意の $p\in K[t]$ に対し、$p=qg+r$ かつ「$r=0$ または $\deg r<\deg g$」を満たす $q,r\in K[t]$ がただ 1 組存在する。$r$ を $p$ を $g$ で割った**余り**という。
:::

:::proof
1. 最高次係数をそれぞれ $a,b$ とすると、$pq$ の $t^{\deg p+\deg q}$ の係数は $ab\ne0$ であり、それより高い次数の係数は $0$ である。$pq=pr$ ならば $p(q-r)=0$ なので $q-r=0$。

2. 存在：$p$ の次数に関する帰納法で示す。$p=0$ または $\deg p<\deg g$ ならば $q=0$, $r=p$ とすればよい。$\deg p=d\ge e=\deg g$ とし、最高次係数をそれぞれ $a,b$ とすると、$p_1\defeq p-\frac{a}{b}t^{d-e}g$ は $t^d$ の係数が $0$ なので、$p_1=0$ または $\deg p_1<d$ である。帰納法の仮定より $p_1=q_1g+r$ と書けるので、$p=(q_1+\frac abt^{d-e})g+r$。一意性：$qg+r=q'g+r'$ ならば $(q-q')g=r'-r$ である。$q\ne q'$ ならば 1 より左辺の次数は $\deg g$ 以上だが、右辺は $0$ か次数 $\deg g$ 未満なので矛盾する。よって $q=q'$, $r=r'$。
:::

:::lemma[因数定理と根の個数]{#lem-linalg-factor-theorem}
$p\in K[t]$, $c\in K$ とする。

1. （因数定理）$p(c)=0$ であるための必要十分条件は、$p=(t-c)q$ となる $q\in K[t]$ が存在することである。
2. $p\ne0$, $\deg p=d$ ならば、$p(c)=0$ となる $c\in K$（$p$ の**根**）は高々 $d$ 個である。
3. $p\ne0$ ならば、$p=(t-c)^mq$ かつ $q(c)\ne0$ を満たす整数 $m\ge0$ と $q\in K[t]$ がただ 1 組存在する。$m$ は $(t-c)^m$ が $p$ を割り切るような最大の整数であり、これを根 $c$ の**重複度**という。
4. $p=a\prod_{i=1}^d(t-c_i)$（$a\ne0$）と一次式の積に分解されるとき、$p$ の根は $c_1,\dots,c_d$ のいずれかであり、根 $c$ の重複度は $c_i=c$ となる $i$ の個数に等しい。
5. $fg=\prod_{i=1}^d(t-c_i)$ ならば、$f$ はある部分集合 $I\subset\{1,\dots,d\}$ と $a\in K$ を用いて $f=a\prod_{i\in I}(t-c_i)$ と書ける。すなわち、一次式の積に分解する多項式の約数も一次式の積に分解する。
:::

:::proof
1. :ref[lem-linalg-polynomial-division] より $p=(t-c)q+r$, $r$ は定数（$r=0$ または次数 $0$）と書ける。:ref[lem-linalg-polynomial-substitution] の 1（$n=1$）より両辺に $c$ を代入すると $p(c)=0\cdot q(c)+r=r$ なので、$p(c)=0\iff r=0$。

2. $d$ に関する帰納法。$d=0$ ならば $p$ は $0$ でない定数で根をもたない。$d\ge1$ で根 $c$ があれば、1 より $p=(t-c)q$, $\deg q=d-1$ である。$c'\ne c$ が根ならば $0=p(c')=(c'-c)q(c')$ より $q(c')=0$。よって $p$ の根は $c$ と $q$ の根（高々 $d-1$ 個）に限られる。

3. 存在：$p(c)\ne0$ ならば $m=0$, $q=p$。そうでなければ 1 により $(t-c)$ をくくり出す操作を繰り返す。各段階で次数が下がるので有限回で $q(c)\ne0$ となる。一意性：$(t-c)^mq=(t-c)^{m'}q'$（$q(c),q'(c)\ne0$）で $m<m'$ とすると、:ref[lem-linalg-polynomial-division] の 1 により $(t-c)^m$ を約分して $q=(t-c)^{m'-m}q'$、したがって $q(c)=0$ となり矛盾する。最後に $(t-c)^{m+1}$ が $p$ を割り切るとすると、$p=(t-c)^{m+1}s$ と $p=(t-c)^mq$ から同様に約分して $q=(t-c)s$、$q(c)=0$ となり矛盾する。$(t-c)^m$ は明らかに $p$ を割り切る。

4. $p(c)=a\prod_i(c-c_i)$ であり、体では積が $0$ になるのはいずれかの因子が $0$ のときに限るので、根は $c_i$ のいずれかである。$c_i=c$ となる $i$ の個数を $k$ とすれば $p=(t-c)^k\cdot a\prod_{c_i\ne c}(t-c_i)$ であり、後者は $c$ で $0$ にならないので、3 の一意性から重複度は $k$ である。

5. $d$ に関する帰納法。$d=0$ ならば $fg=1$ で、1 より $f,g$ は定数である。$d\ge1$ とする。$f(c_d)g(c_d)=0$ なので $f(c_d)=0$ または $g(c_d)=0$ である。$f(c_d)=0$ ならば $f=(t-c_d)f_1$ と書け、$(t-c_d)$ を約分して $f_1g=\prod_{i<d}(t-c_i)$ となるので、帰納法の仮定を $f_1$ に適用すればよい。$g(c_d)=0$ ならば同様に $g=(t-c_d)g_1$ として $fg_1=\prod_{i<d}(t-c_i)$ に帰納法の仮定を適用する。
:::

:::remark
$K=\R$ や $\C$ のような無限体では、2 により、値がすべての $c\in K$ で一致する 2 つの多項式は等しい（差が無限個の根をもつ）。しかし有限体では事情が異なり、例えば $\F_2$ 上の $t^2+t$ は零多項式でないが、$t=0,1$ のどちらでも値 $0$ をとる。多項式と多項式関数を区別するのはこのためである。
:::

### 特性多項式

$A\in M_n(K)$ に対し、$\lambda$ が固有値であることは $(\lambda I_n-A)x=0$ が自明でない解をもつこと、すなわち $\det(\lambda I_n-A)=0$ と同値である（:ref[cor-linalg-det-invertible]）。そこで $\lambda$ を変数 $t$ に置き換えた多項式を考える。

:::definition[特性多項式]{#def-linalg-characteristic-polynomial}
$A\in M_n(K)$ に対し、$tI_n-A$ を多項式環 $K[t]$ の元を成分とする行列とみて、その行列式（:ref[prop-linalg-det-ring]）

$$
\Phi_A(t)\defeq\det(tI_n-A)\in K[t] \label{eq-linalg-charpoly}
$$

を $A$ の :term[特性多項式]{en="characteristic polynomial" yomi="とくせいたこうしき"}（固有多項式）という。方程式 $\Phi_A(t)=0$ を**特性方程式**（固有方程式）という。
:::

:::proposition[特性多項式の係数]{#prop-linalg-charpoly-coefficients}
$\Phi_A$ は $n$ 次のモニック多項式であり、

$$
\Phi_A(t)=t^n-(\tr A)\,t^{n-1}+\cdots+(-1)^n\det A
$$

である。すなわち $t^{n-1}$ の係数は $-\tr A$、定数項は $(-1)^n\det A$ である。特に $n=2$ ならば $\Phi_A(t)=t^2-(\tr A)t+\det A$ である。
:::

:::proof
$tI_n-A$ の $(i,j)$ 成分は $t\delta_{ij}-a_{ij}$ であり、$t$ を含むのは対角成分だけである。Leibniz の公式 :eqref[eq-linalg-leibniz] の各項は各列から 1 つずつ成分を選んだ積である。$\sigma\ne\id$ ならば $\sigma(j)\ne j$ となる $j$ が少なくとも 2 つある（1 つだけでは $\sigma$ が全単射にならない）ので、その項は対角成分を高々 $n-2$ 個しか含まず、次数は $n-2$ 以下である。したがって $t^n$, $t^{n-1}$ の係数は $\sigma=\id$ の項 $\prod_i(t-a_{ii})$ だけから決まり、それぞれ $1$ と $-\sum_ia_{ii}=-\tr A$ である。定数項は $t=0$ を代入した値 $\Phi_A(0)=\det(-A)=(-1)^n\det A$ である（:ref[prop-linalg-det-ring] の 5 と :ref[prop-linalg-det-operations] の 2）。
:::

:::proposition[特性多項式の根と相似不変性]{#prop-linalg-charpoly-similarity}
1. $\lambda\in K$ が $A$ の固有値であるための必要十分条件は $\Phi_A(\lambda)=0$ である。
2. $B=P^{-1}AP$（$P$ は正則）ならば $\Phi_B=\Phi_A$ である。さらに $\lambda\in K$ について $L_P$ は $E_\lambda(B)$ から $E_\lambda(A)$ への同型写像を与える。
:::

:::proof
1. :ref[prop-linalg-det-ring] の 5（$t=\lambda$ の代入）より $\Phi_A(\lambda)=\det(\lambda I_n-A)$ である。$(\lambda I_n-A)x=0$ が自明でない解をもつことは :ref[cor-linalg-det-invertible] より $\det(\lambda I_n-A)=0$ と同値である。

2. $K[t]$ 上の行列として $tI_n-B=P^{-1}(tI_n-A)P$ である。:ref[prop-linalg-det-ring] の 2（積公式）より $\Phi_B=\det(P^{-1})\Phi_A\det P=\Phi_A$。また $Bx=\lambda x\iff A(Px)=\lambda Px$ なので $L_P(E_\lambda(B))=E_\lambda(A)$ であり、$L_P$ は同型写像である。
:::

:::definition[線形変換の特性多項式]{#def-linalg-operator-charpoly}
$T\in\End(V)$ に対し、$V$ の基底 $\mathcal B$ を 1 つ選んで $\Phi_T\defeq\Phi_{[T]_{\mathcal B}}$ と定め、$T$ の特性多項式という。:ref[prop-linalg-similarity-equivalence] と :ref[prop-linalg-charpoly-similarity] の 2 により、これは基底の選び方によらない。:ref[prop-linalg-eigen-coordinates] より、$T$ の固有値は $\Phi_T$ の $K$ における根と一致する。
:::

:ref[lem-linalg-factor-theorem] の 2 より、$T$ の固有値は高々 $n$ 個である（:ref[prop-linalg-eigenspaces-direct-sum] の別証明）。

:::corollary[不変部分空間と特性多項式]{#cor-linalg-invariant-charpoly}
:ref[prop-linalg-invariant-block] の状況で $\Phi_T=\Phi_{T|_W}\cdot\Phi_{A_2}$ が成り立つ。特に $\Phi_{T|_W}$ は $\Phi_T$ を割り切る。
:::

:::proof
$tI_n-[T]_{\mathcal B}$ は対角ブロック $tI_k-A_1$, $tI_{n-k}-A_2$ をもつ $K[t]$ 上のブロック上三角行列なので、:ref[prop-linalg-det-ring] の 2（$K[t]$ 上の :ref[prop-linalg-det-block-triangular]）より $\Phi_T=\det(tI_k-A_1)\det(tI_{n-k}-A_2)=\Phi_{T|_W}\Phi_{A_2}$ である。
:::

:::example[回転と係数体]{#exm-linalg-rotation-eigen}
$R=\begin{pmatrix}0&-1\\1&0\end{pmatrix}$（$90^\circ$ 回転、:ref[exm-linalg-rotation]）について

$$
\Phi_R(t)=\det\begin{pmatrix}t&1\\-1&t\end{pmatrix}=t^2+1
$$

である。$t^2+1=0$ は実数解をもたないので、$R$ は $\R$ 上では固有値をもたない。一方、同じ行列を $\C^2$ に作用させると固有値 $i,-i$ をもつ。$Rx=\lambda x$ の第 1 成分は $-y=\lambda x$ であり、第 2 成分 $x=\lambda y$ は $\lambda^2=-1$ のもとで第 1 成分から従うので、

$$
E_i(R)=\bigl\langle\transp{(1,-i)}\bigr\rangle,\qquad
E_{-i}(R)=\bigl\langle\transp{(1,i)}\bigr\rangle
$$

である。固有値の有無は、変換をどの体の上で考えるかによって変わる。
:::

:::proposition[三角行列の固有値]{#prop-linalg-triangular-eigenvalues}
$A$ が上三角行列または下三角行列ならば $\Phi_A(t)=\prod_{i=1}^n(t-a_{ii})$ であり、$A$ の固有値は対角成分 $a_{11},\dots,a_{nn}$ に限る。
:::

:::proof
$tI_n-A$ も（$K[t]$ 上の）三角行列で、対角成分は $t-a_{ii}$ なので、:ref[prop-linalg-det-ring] の 2（:ref[prop-linalg-det-block-triangular]）から従う。根については :ref[lem-linalg-factor-theorem] の 4 による。
:::

:::proposition[実行列の複素固有値]{#prop-linalg-real-conjugate}
$A\in M_n(\R)$ を $\C^n$ に作用する複素行列とみる。$\lambda\in\C$ が $A$ の固有値で $v\in\C^n$ がそれに属する固有ベクトルならば、$\conj{\lambda}$ も固有値で、$\conj{v}$（成分ごとの複素共役）がそれに属する固有ベクトルである。さらに $\lambda$ と $\conj\lambda$ の代数的重複度（$\Phi_A$ の根としての重複度）は等しい。
:::

:::proof
$Av=\lambda v$ の両辺の複素共役をとると、$A$ の成分は実数なので $A\conj{v}=\conj{\lambda}\,\conj{v}$ であり、$\conj v\ne0$ である。重複度については、$\Phi_A$ は実係数なので、多項式 $q$ の係数の共役をとったものを $\conj q$ と書くと $\conj{\Phi_A}=\Phi_A$ である。$\Phi_A=(t-\lambda)^mq$, $q(\lambda)\ne0$ ならば、係数の共役をとって $\Phi_A=(t-\conj\lambda)^m\conj q$, $\conj q(\conj\lambda)=\conj{q(\lambda)}\ne0$ を得るので、:ref[lem-linalg-factor-theorem] の 3 の一意性より $\conj\lambda$ の重複度も $m$ である。
:::

### 二種類の重複度

:::definition[代数的重複度と幾何学的重複度]{#def-linalg-eigenvalue-multiplicities}
$\lambda$ を $T$ の固有値とする。$\lambda$ の $\Phi_T$ の根としての重複度（:ref[lem-linalg-factor-theorem] の 3）、すなわち $(t-\lambda)^m$ が $\Phi_T$ を割り切る最大の整数 $m$ を、$\lambda$ の :term[代数的重複度]{en="algebraic multiplicity" yomi="だいすうてきちょうふくど"} といい $m_\lambda$ と書く。一方、$g_\lambda\defeq\dim E_\lambda(T)$ を $\lambda$ の :term[幾何学的重複度]{en="geometric multiplicity" yomi="きかがくてきちょうふくど"} という。行列についても同様に定める。
:::

:::proposition[重複度の不等式]{#prop-linalg-multiplicity-bound}
任意の固有値 $\lambda$ に対して $1\le g_\lambda\le m_\lambda$ が成り立つ。
:::

:::proof
固有ベクトルが存在するので $g\defeq g_\lambda\ge1$ である。$W=E_\lambda(T)$ は $T$ 不変で、$T|_W=\lambda\id_W$ である。:ref[prop-linalg-invariant-block] と :ref[cor-linalg-invariant-charpoly] により、$W$ の基底を延長した基底 $\mathcal B$ で

$$
[T]_{\mathcal B}=\begin{pmatrix}\lambda I_g&H\\O&B\end{pmatrix},
\qquad
\Phi_T(t)=(t-\lambda)^g\,\Phi_B(t)
$$

となる（$g=n$ のときは $\Phi_B=1$ と読む）。よって $(t-\lambda)^g$ は $\Phi_T$ を割り切り、$m_\lambda$ の最大性から $g\le m_\lambda$ である。
:::

:::example[同じ特性多項式、異なる固有空間]{#exm-linalg-multiplicities}
2 つの行列

$$
A=\begin{pmatrix}2&0\\0&2\end{pmatrix},\qquad
B=\begin{pmatrix}2&1\\0&2\end{pmatrix}
$$

はどちらも特性多項式 $(t-2)^2$ をもち、$m_2=2$ である。$A-2I_2=O$ なので $E_2(A)=K^2$, $g_2=2$ である。一方 $(B-2I_2)\transp{(x,y)}=\transp{(y,0)}$ なので $E_2(B)=\langle e_1\rangle$, $g_2=1$ である。:ref[prop-linalg-charpoly-similarity] の 2 より相似な行列の固有空間の次元は等しいので、$A$ と $B$ は相似でない。特性多項式は相似の不変量だが、相似類を完全には決定しない。
:::

## 対角化可能性の判定

### 対角化可能性の必要十分条件

:::definition[対角化可能]{#def-linalg-diagonalizable}
$A\in M_n(K)$ について、$P^{-1}AP$ が対角行列となる正則行列 $P\in M_n(K)$ が存在するとき、$A$ は $K$ 上 :term[対角化可能]{en="diagonalizable" yomi="たいかくかかのう"} であるといい、そのような $P$ を求めることを $A$ の :term[対角化]{en="diagonalization" yomi="たいかくか"} という。線形変換 $T\in\End(V)$ については、ある基底 $\mathcal B$ に関する表現行列 $[T]_{\mathcal B}$ が対角行列になるとき、$T$ は対角化可能であるという。
:::

:ref[prop-linalg-similarity-equivalence] により、$A$ が対角化可能であることと $L_A$ が対角化可能であることは同値である。

:::theorem[対角化の必要十分条件]{#thm-linalg-diagonalization-criterion}
$T\in\End(V)$ の $K$ における相異なる固有値の全体を $\lambda_1,\ldots,\lambda_r$ とする。次の条件は互いに同値である。

1. $T$ は対角化可能である。
2. $V$ は $T$ の固有ベクトルからなる基底をもつ。
3. $V=E_{\lambda_1}(T)\oplus\cdots\oplus E_{\lambda_r}(T)$。
4. $\Phi_T$ は $K$ 上で一次式の積に分解し、各固有値 $\lambda_i$ について $g_{\lambda_i}=m_{\lambda_i}$ である。
5. $\sum_{i=1}^rg_{\lambda_i}=n$。
:::

:::proof
(1)$\Leftrightarrow$(2)：基底 $\mathcal B=(v_1,\dots,v_n)$ に関する表現行列が対角行列 $\diag(\mu_1,\dots,\mu_n)$ であることは、第 $j$ 列が $\mu_je_j$、すなわち $T(v_j)=\mu_jv_j$（$1\le j\le n$）であることと同値である。

(2)$\Rightarrow$(3)：各基底ベクトルはいずれかの固有空間に属するので、$V=\sum_iE_{\lambda_i}(T)$ であり、:ref[prop-linalg-eigenspaces-direct-sum] よりこの和は直和である。

(3)$\Rightarrow$(5)：:ref[prop-linalg-direct-sum-criterion] の 2 による。

(5)$\Rightarrow$(2)：:ref[prop-linalg-eigenspaces-direct-sum] より $W\defeq\bigoplus_iE_{\lambda_i}(T)$ の次元は $\sum_ig_{\lambda_i}=n$ なので、:ref[cor-linalg-subspace-dimension] より $W=V$ である。各 $E_{\lambda_i}(T)$ の基底を並べたものは $W=V$ の基底であり（:ref[prop-linalg-direct-sum-criterion] の 2 の証明）、固有ベクトルからなる。

(2)$\Rightarrow$(4)：(1) より、ある基底で $[T]=\diag(\mu_1,\dots,\mu_n)$ であり、各 $\mu_j$ は固有値である。$\lambda_i$ が $\mu_1,\dots,\mu_n$ の中に $k_i$ 回現れるとすると、:ref[prop-linalg-triangular-eigenvalues] より $\Phi_T=\prod_j(t-\mu_j)$ は一次式の積に分解し、:ref[lem-linalg-factor-theorem] の 4 より $m_{\lambda_i}=k_i$ である。一方、$\lambda_i$ に対応する $k_i$ 個の基底ベクトルは $E_{\lambda_i}(T)$ の一次独立な元なので $g_{\lambda_i}\ge k_i=m_{\lambda_i}$ であり、:ref[prop-linalg-multiplicity-bound] と合わせて $g_{\lambda_i}=m_{\lambda_i}$ である。

(4)$\Rightarrow$(5)：$\Phi_T=\prod_{j=1}^n(t-c_j)$ と分解すると、:ref[lem-linalg-factor-theorem] の 4 より根 $c_j$ は $T$ の固有値 $\lambda_1,\dots,\lambda_r$ のいずれかであり、$\lambda_i$ の重複度 $m_{\lambda_i}$ は $c_j=\lambda_i$ となる $j$ の個数である。よって $\sum_ig_{\lambda_i}=\sum_im_{\lambda_i}=n$。
:::

:::corollary[相異なる固有値が n 個あれば対角化可能]{#cor-linalg-distinct-eigenvalues-diagonalizable}
$T$ が $K$ に $n$ 個の相異なる固有値をもてば、$T$ は対角化可能である。
:::

:::proof
各固有値について $g_{\lambda_i}\ge1$ なので $\sum_ig_{\lambda_i}\ge n$ であり、:ref[prop-linalg-eigenspaces-direct-sum] より $\le n$ なので等号が成り立つ。:ref[thm-linalg-diagonalization-criterion] の (5) による。
:::

この系の逆は成り立たない（$I_n$ は固有値 $1$ のみをもつが対角行列である）。対角化可能性の判定と対角化は、次の手順で行う。

1. $\Phi_A(t)$ を計算し、$K$ における根（固有値）$\lambda_1,\dots,\lambda_r$ とその重複度を求める。$\Phi_A$ が $K$ 上で一次式の積に分解しなければ、$K$ 上では対角化できない。
2. 各 $\lambda_i$ について斉次方程式 $(A-\lambda_iI_n)x=0$ を解き、$E_{\lambda_i}(A)$ の基底を求める。:ref[cor-linalg-solution-space-dim] より $g_{\lambda_i}=n-\rank(A-\lambda_iI_n)$ である。単根（$m_{\lambda_i}=1$）では :ref[prop-linalg-multiplicity-bound] より必ず $g_{\lambda_i}=1$ なので、確認が必要なのは重根だけである。
3. すべての $i$ で $g_{\lambda_i}=m_{\lambda_i}$ ならば、各固有空間の基底を並べた行列 $P$ によって $P^{-1}AP$ は対角行列になる。対角成分には、$P$ の列の順に対応する固有値が並ぶ。

:::example[重根をもつが対角化可能な行列]{#exm-linalg-diagonalization-3x3}
$A=\begin{pmatrix}2&1&1\\1&2&1\\1&1&2\end{pmatrix}$ を対角化する。$J$ をすべての成分が $1$ の 3 次正方行列とすると $A=I_3+J$ であり、

$$
\Phi_A(t)=\det\begin{pmatrix}t-2&-1&-1\\-1&t-2&-1\\-1&-1&t-2\end{pmatrix}=(t-1)^2(t-4)
$$

である（第 2, 3 列を第 1 列に加えると第 1 列が $(t-4)\transp{(1,1,1)}$ になり、$t-4$ をくくり出した後、第 2, 3 行から第 1 行を引けば上三角行列になる）。$A-I_3=J$ の階数は $1$ なので $g_1=3-1=2=m_1$ であり、$E_1(A)=\{x\mid x_1+x_2+x_3=0\}$ の基底として $\transp{(1,-1,0)},\transp{(1,0,-1)}$ がとれる。$E_4(A)$ は $\transp{(1,1,1)}$ で張られる。よって

$$
P=\begin{pmatrix}1&1&1\\-1&0&1\\0&-1&1\end{pmatrix},\qquad P^{-1}AP=\diag(1,1,4)
$$

である（$\det P=3\ne0$）。
:::

### 対角化の応用

$A=PDP^{-1}$, $D=\diag(\lambda_1,\dots,\lambda_n)$ ならば、:ref[lem-linalg-polynomial-substitution] の 3 より、任意の多項式 $p$ について $p(A)=P\,\diag(p(\lambda_1),\dots,p(\lambda_n))\,P^{-1}$ である。特に $A^k=PD^kP^{-1}$ であり、$D^k=\diag(\lambda_1^k,\dots,\lambda_n^k)$ である。

:::example[行列の冪]{#exm-linalg-matrix-power}
:ref[exm-linalg-eigen-2x2] の $A=\begin{pmatrix}2&1\\1&2\end{pmatrix}$ の特性多項式は $(t-2)^2-1=(t-3)(t-1)$ である。固有値 $3,1$ に属する固有ベクトルをこの順に並べて

$$
P=\begin{pmatrix}1&1\\1&-1\end{pmatrix},\qquad
D=\diag(3,1)
$$

とおくと $AP=PD$ であり、$P^2=2I_2$ より $P^{-1}=\tfrac12P$ である。したがって任意の $k\ge0$ について

$$
A^k=PD^kP^{-1}
=\frac12\begin{pmatrix}3^k+1&3^k-1\\3^k-1&3^k+1\end{pmatrix}
$$

となる。例えば $A^3e_1=\transp{(14,13)}$ である。これは $e_1=\tfrac12\transp{(1,1)}+\tfrac12\transp{(1,-1)}$ と固有ベクトルの和に分解し、第 1 の成分だけを $3^3$ 倍したものである。
:::

:::example[Fibonacci 数列]{#exm-linalg-fibonacci}
$F_0=0$, $F_1=1$, $F_{k+2}=F_{k+1}+F_k$ で定まる Fibonacci 数列の一般項を求める（$K=\R$）。$x_k\defeq\transp{(F_{k+1},F_k)}$ とおくと

$$
x_{k+1}=Ax_k,\qquad A=\begin{pmatrix}1&1\\1&0\end{pmatrix}
$$

なので $x_k=A^kx_0$, $x_0=\transp{(1,0)}$ である。$\Phi_A(t)=t^2-t-1$ の根は $\varphi=\frac{1+\sqrt5}{2}$, $\psi=\frac{1-\sqrt5}{2}$ で、相異なるので $A$ は対角化可能である。$\varphi^2=\varphi+1$ より $A\transp{(\varphi,1)}=\transp{(\varphi+1,\varphi)}=\varphi\transp{(\varphi,1)}$ であり、同様に $\transp{(\psi,1)}$ は $\psi$ に属する固有ベクトルである。$x_0=\frac{1}{\sqrt5}\bigl(\transp{(\varphi,1)}-\transp{(\psi,1)}\bigr)$（$\varphi-\psi=\sqrt5$）と分解すれば

$$
x_k=\frac{1}{\sqrt5}\Bigl(\varphi^k\begin{pmatrix}\varphi\\1\end{pmatrix}-\psi^k\begin{pmatrix}\psi\\1\end{pmatrix}\Bigr),\qquad
F_k=\frac{\varphi^k-\psi^k}{\sqrt5}
$$

を得る（**Binet の公式**）。$\abs\psi<1$ なので、$k$ が大きいとき $F_k\approx\varphi^k/\sqrt5$ であり、$F_{k+1}/F_k\to\varphi$（黄金比）である。一般に、漸化式 $x_{k+1}=Ax_k$ の長期的な振る舞いは、絶対値最大の固有値とその固有ベクトルによって支配される。
:::

同様に、定数係数の線形微分方程式系 $x'(s)=Ax(s)$ は、$A=PDP^{-1}$ ならば $y=P^{-1}x$ とおくことで $y_i'=\lambda_iy_i$ という独立な方程式に分解され、$y_i(s)=e^{\lambda_is}y_i(0)$ と解ける。対角化できない場合も含めた一般論は行列の指数関数を用いて [行列解析](/math/linear-algebra/matrix-analysis/) で扱う。

## Cayley–Hamilton の定理

### 定理と証明

行列の多項式 $p(A)$ は :ref[lem-linalg-polynomial-substitution] で定義した。$n$ 次正方行列全体 $M_n(K)$ は $n^2$ 次元なので、$n^2+1$ 個の行列 $I_n,A,\dots,A^{n^2}$ は一次従属であり、$p(A)=O$ となる $0$ でない多項式 $p$（次数 $n^2$ 以下）が必ず存在する。実は、次数 $n$ の多項式 $\Phi_A$ がそのような多項式になっている。

:::theorem[Cayley–Hamilton の定理]{#thm-linalg-cayley-hamilton}
任意の $A\in M_n(K)$ について、その特性多項式に $A$ を代入すると零行列になる。すなわち

$$
\Phi_A(A)=O. \label{eq-linalg-cayley-hamilton}
$$

同様に、任意の $T\in\End(V)$ について $\Phi_T(T)=0$ である。
:::

:::proof
$\Phi_A(t)=\sum_{k=0}^n c_kt^k$（$c_n=1$）とする。$K[t]$ を成分とする行列 $tI_n-A$ の余因子行列 $\operatorname{adj}(tI_n-A)$ を考える。その成分は $tI_n-A$ の $(n-1)$ 次小行列式（の符号を変えたもの）であり、Leibniz の公式の各項は $t$ を含む成分を高々 $n-1$ 個しか含まないので、次数 $n-1$ 以下の多項式である。そこで、各成分の $t^k$ の係数を並べた行列を $B_k\in M_n(K)$ とすれば

$$
\operatorname{adj}(tI_n-A)=\sum_{k=0}^{n-1}t^kB_k
$$

と書ける。:ref[prop-linalg-det-ring] の 3（$K[t]$ 上の :eqref[eq-linalg-adjugate]）より

$$
(tI_n-A)\sum_{k=0}^{n-1}t^kB_k=\Phi_A(t)\,I_n=\sum_{k=0}^nc_kt^kI_n
$$

が $M_n(K[t])$ において成り立つ。左辺の $(i,j)$ 成分は $\sum_kt^{k+1}(B_k)_{ij}-\sum_kt^k(AB_k)_{ij}$ である。2 つの多項式が等しいことは係数がすべて等しいことなので、$B_{-1}\defeq O$, $B_n\defeq O$ とおいて $t^k$ の係数を比較すると

$$
B_{k-1}-AB_k=c_k I_n\qquad(0\le k\le n)
$$

を得る。第 $k$ 式に左から $A^k$ を掛け、$k=0$ から $n$ まで加えると

$$
\Phi_A(A)=\sum_{k=0}^n c_k A^k
=\sum_{k=0}^n A^kB_{k-1}-\sum_{k=0}^n A^{k+1}B_k
$$

となる。第 1 の和は $k=j+1$ と書き直すと $\sum_{j=-1}^{n-1}A^{j+1}B_j$、第 2 の和は $\sum_{j=0}^nA^{j+1}B_j$ であり、$B_{-1}=B_n=O$ なので両者は同じ $\sum_{j=0}^{n-1}A^{j+1}B_j$ である。よって $\Phi_A(A)=O$ である。この計算では $A$ と $B_j$ の積の順序を入れ替えていないことに注意する。

線形変換については、基底 $\mathcal B$ をとって $A=[T]_{\mathcal B}$ とすれば、:ref[lem-linalg-polynomial-substitution] の 4 より $[\Phi_T(T)]_{\mathcal B}=\Phi_A(A)=O$ なので、$\Phi_T(T)=0$ である。
:::

:::remark
「$\Phi_A(t)=\det(tI_n-A)$ に $t=A$ を代入すれば $\det(AI_n-A)=\det O=0$ だから明らか」という議論は**誤り**である。$\Phi_A(A)$ は $n$ 次正方行列であり、スカラー $\det(A-A)=0$ とは別物である。$\det(tI_n-A)$ を展開してから $t$ に $A$ を代入することと、行列 $tI_n-A$ の成分の $t$ に $A$ を代入してから行列式をとることは、意味がまったく異なる（後者は成分が行列の「行列」になってしまう）。例えば $n=2$ では、定理の主張は

$$
A^2-(\tr A)A+(\det A)I_2=O
$$

という行列の等式であり、成分を計算して直接確かめることができる。
:::

### Cayley–Hamilton の定理の応用

:::corollary[行列の多項式の次数の引き下げ]{#cor-linalg-ch-reduction}
$A\in M_n(K)$ とする。

1. 任意の $p\in K[t]$ について、$p$ を $\Phi_A$ で割った余りを $r$ とすると $p(A)=r(A)$ である。特に $p(A)$ は $I_n,A,\dots,A^{n-1}$ の一次結合である。
2. $\Phi_A(t)=t^n+c_{n-1}t^{n-1}+\cdots+c_1t+c_0$ とする。$A$ が正則であることは $c_0\ne0$ と同値であり、このとき

$$
A^{-1}=-\frac{1}{c_0}\bigl(A^{n-1}+c_{n-1}A^{n-2}+\cdots+c_1I_n\bigr)
$$

である。特に $A^{-1}$ は $A$ の多項式である。
:::

:::proof
1. $p=q\Phi_A+r$ ならば、:ref[lem-linalg-polynomial-substitution] の 1 と Cayley–Hamilton の定理より $p(A)=q(A)\Phi_A(A)+r(A)=r(A)$。
2. $c_0=(-1)^n\det A$（:ref[prop-linalg-charpoly-coefficients]）なので、前半は :ref[cor-linalg-det-invertible] による。$\Phi_A(A)=O$ を $A(A^{n-1}+c_{n-1}A^{n-2}+\cdots+c_1I_n)=-c_0I_n$ と書き直し、$-c_0$ で割ればよい（左右の因子は可換である）。
:::

:::example[対角化できない行列の冪]{#exm-linalg-ch-power}
$B=\begin{pmatrix}2&1\\0&2\end{pmatrix}$ は :ref[exm-linalg-multiplicities] で見たように対角化できない。$\Phi_B(t)=(t-2)^2$ なので、定理より $N\defeq B-2I_2$ は $N^2=O$ を満たす。$2I_2$ と $N$ は可換なので二項展開が使え、$k\ge1$ に対し $N$ を 2 つ以上含む項が消えて

$$
B^k=(2I_2+N)^k=2^k I_2+k2^{k-1}N
=\begin{pmatrix}2^k&k2^{k-1}\\0&2^k\end{pmatrix}
$$

となる。冪の計算に対角化は必須ではない。別の方法として、$t^k$ を $(t-2)^2$ で割った余り $r(t)=at+b$ を求めてもよい。$t^k=q(t)(t-2)^2+at+b$ に $t=2$ を代入し、両辺を微分してから $t=2$ を代入すると $2^k=2a+b$, $k2^{k-1}=a$ となり、$B^k=aB+bI_2$ から同じ結果を得る。
:::

Cayley–Hamilton の定理により $\Phi_A$ は「$A$ を代入すると $0$ になる多項式」の 1 つであるが、そのような多項式のうち次数が最小のものはより低い次数でありうる（例えば $I_n$ は $t-1$ を代入すると $O$ になる）。この最小のものを最小多項式といい、対角化可能性の別の判定法を与える（[最小多項式と Jordan 標準形](/math/linear-algebra/jordan-form/)）。

## 三角化

### 三角化可能性

対角化はいつも可能とは限らない。しかし、上三角行列までなら、特性多項式が一次式の積に分解しさえすれば常に到達できる。

:::definition[三角化]{#def-linalg-triangularization}
$A\in M_n(K)$ について、$P^{-1}AP$ が上三角行列となる正則行列 $P\in M_n(K)$ が存在するとき、$A$ は $K$ 上**三角化可能**であるといい、そのような $P$ を求めることを :term[三角化]{en="triangularization" yomi="さんかくか"} という。線形変換 $T$ については、表現行列が上三角行列となる基底が存在するとき三角化可能という。
:::

基底 $(v_1,\ldots,v_n)$ に関する $T$ の表現行列が上三角であることは、各 $j$ について

$$
T(v_j)\in\langle v_1,\ldots,v_j\rangle
$$

が成り立つこと、すなわち部分空間の列 $V_j\defeq\langle v_1,\dots,v_j\rangle$（$\{0\}=V_0\subset V_1\subset\cdots\subset V_n=V$, $\dim V_j=j$）がすべて $T$ 不変であることと同値である。実際、第 $j$ 列の第 $j+1$ 成分以下が $0$ であることが $T(v_j)\in V_j$ を意味する。

:::theorem[三角化と一次式への分解]{#thm-linalg-triangularization}
$A\in M_n(K)$ が $K$ 上で三角化可能であるための必要十分条件は、$\Phi_A$ が $K$ 上で一次式の積に分解することである。特に、複素正方行列はすべて $\C$ 上で三角化可能である。線形変換についても同様である。
:::

:::proof
**必要性**：$P^{-1}AP=U$ が対角成分 $d_1,\dots,d_n$ の上三角行列ならば、:ref[prop-linalg-charpoly-similarity] と :ref[prop-linalg-triangular-eigenvalues] より $\Phi_A=\Phi_U=\prod_j(t-d_j)$ である。

**十分性**：$\Phi_A=\prod_{i=1}^n(t-\lambda_i)$ と仮定し、$n$ に関する帰納法で示す。$n=1$ の行列はすでに上三角である。$n\ge2$ とする。$\lambda\defeq\lambda_1$ は $\Phi_A$ の根なので固有値であり（:ref[prop-linalg-charpoly-similarity]）、固有ベクトル $v_1\in K^n$ がとれる。$\langle v_1\rangle$ は $L_A$ 不変なので、$v_1$ を $K^n$ の基底に延長し、それを列に並べた正則行列を $Q_0$ とすると、:ref[prop-linalg-invariant-block] と :ref[cor-linalg-invariant-charpoly] より

$$
Q_0^{-1}AQ_0=\begin{pmatrix}\lambda&h\\O&B\end{pmatrix}\qquad(h\in M_{1,n-1}(K),\ B\in M_{n-1}(K)),
\qquad
\Phi_A=(t-\lambda)\Phi_B
$$

となる。$(t-\lambda)\Phi_B=(t-\lambda)\prod_{i\ge2}(t-\lambda_i)$ から :ref[lem-linalg-polynomial-division] の 1 により $(t-\lambda)$ を約分して $\Phi_B=\prod_{i\ge2}(t-\lambda_i)$ を得る。帰納法の仮定により、$Q^{-1}BQ$ が上三角となる正則行列 $Q\in M_{n-1}(K)$ がある。$S\defeq\begin{pmatrix}1&O\\O&Q\end{pmatrix}$ とおくと $S^{-1}=\begin{pmatrix}1&O\\O&Q^{-1}\end{pmatrix}$ であり、ブロックの積を計算すると

$$
S^{-1}(Q_0^{-1}AQ_0)S=\begin{pmatrix}\lambda&hQ\\O&Q^{-1}BQ\end{pmatrix}
$$

は上三角行列である。$P\defeq Q_0S$ とすればよい。線形変換については、表現行列に適用し :ref[prop-linalg-similarity-equivalence] の 3 を用いればよい。

$K=\C$ の場合：**代数学の基本定理**（次数 $1$ 以上の複素係数多項式は複素数の根をもつ）により、$\Phi_A$ は根 $c$ をもち、:ref[lem-linalg-factor-theorem] の 1 で $\Phi_A=(t-c)q$ と書ける。$q$ に同じ議論を繰り返せば、$\Phi_A$ は一次式の積に分解する。代数学の基本定理は複素解析の定理であり、その証明は [複素解析](/math/complex-analysis/) の Liouville の定理の節（または :cite[ahlfors-complex]）を参照されたい。
:::

この証明では、最初の固有ベクトル以外の基底ベクトルを固有ベクトルに選ぶ必要はない。右上の $h$ が $0$ でなくても、右下のブロックを三角化すればよいからである。:ref[exm-linalg-multiplicities] の $B=\begin{pmatrix}2&1\\0&2\end{pmatrix}$ は上三角だが対角化できない。対角化可能性は三角化可能性より真に強い条件である。$K=\R$ では、:ref[exm-linalg-rotation-eigen] の回転行列のように、固有値をもたず三角化もできない行列がある。

:::example[多項式の平行移動]{#exm-linalg-shift-triangular}
$K=\R$ または $\C$ とし、$K[t]_{\le2}$ 上の線形変換 $T(p)(t)\defeq p(t+1)$ を考える（代入は和と定数倍を保つので線形である）。基底 $(1,t,t^2)$ の像は $1$, $1+t$, $1+2t+t^2$ なので

$$
[T]=\begin{pmatrix}1&1&1\\0&1&2\\0&0&1\end{pmatrix}
$$

である。各 $k$ について次数 $k$ 以下の多項式の空間は $T$ 不変であり、これが行列の上三角性に対応する。:ref[prop-linalg-triangular-eigenvalues] より $\Phi_T(s)=(s-1)^3$ であり、固有値は $1$ だけで $m_1=3$ である。$p=a+bt+ct^2$ が $T(p)=p$ を満たす条件は $(b+c)+2ct=0$ であり、$b=c=0$ となる。よって $E_1(T)$ は定数多項式全体で $g_1=1<3$ であり、$T$ は三角化可能だが対角化できない。
:::

### 固有値とトレース・行列式

:::corollary[固有値の和と積]{#cor-linalg-trace-det-eigenvalues}
$\Phi_A(t)=\prod_{i=1}^n(t-\lambda_i)$ と分解されるとする（$\lambda_i$ は重複度の分だけ繰り返して並べた固有値）。このとき

$$
\tr A=\sum_{i=1}^n\lambda_i,\qquad \det A=\prod_{i=1}^n\lambda_i
$$

である。特に $K=\C$ ではすべての $A\in M_n(\C)$ についてこれが成り立ち、実行列についても複素固有値を用いて成り立つ。
:::

:::proof
$\prod_i(t-\lambda_i)$ を展開すると、$t^{n-1}$ の係数は $-\sum_i\lambda_i$、定数項は $(-1)^n\prod_i\lambda_i$ である。:ref[prop-linalg-charpoly-coefficients] と係数を比較すればよい。実行列 $A\in M_n(\R)\subset M_n(\C)$ の特性多項式は $\R$ 上で計算しても $\C$ 上で計算しても同じ多項式である。
:::

:::corollary[多項式による固有値の変換]{#cor-linalg-spectral-mapping}
$\Phi_A(t)=\prod_{i=1}^n(t-\lambda_i)$ と分解されるとする。任意の $p\in K[t]$ について

$$
\Phi_{p(A)}(t)=\prod_{i=1}^n\bigl(t-p(\lambda_i)\bigr)
$$

である。特に $p(A)$ の固有値は $p(\lambda_1),\dots,p(\lambda_n)$ であり、$A^k$ の固有値は $\lambda_i^k$ である。また $A$ が正則ならば、$A^{-1}$ の固有値は $\lambda_i^{-1}$ である（重複度も込めて）。
:::

:::proof
:ref[thm-linalg-triangularization] により $P^{-1}AP=U$ が上三角行列となる $P$ があり、:ref[prop-linalg-triangular-eigenvalues] と :ref[lem-linalg-factor-theorem] の 4 より $U$ の対角成分は $\lambda_1,\dots,\lambda_n$ を並べ替えたものである。上三角行列の和・スカラー倍・積は上三角行列であり、積の対角成分は $(UU')_{ii}=\sum_ku_{ik}u'_{ki}=u_{ii}u'_{ii}$（$k<i$ では $u_{ik}=0$、$k>i$ では $u'_{ki}=0$）である。よって $p(U)$ は対角成分 $p(u_{ii})$ の上三角行列であり、:ref[lem-linalg-polynomial-substitution] の 3 より $p(A)=Pp(U)P^{-1}$ なので、$\Phi_{p(A)}=\Phi_{p(U)}=\prod_i(t-p(\lambda_i))$ である。

$A$ が正則ならば、:ref[cor-linalg-ch-reduction] の 2 より $A^{-1}=q(A)$ となる多項式 $q$ がある。$q(U)U=P^{-1}q(A)AP=I_n$ の対角成分を比べると $q(\lambda_i)\lambda_i=1$ なので、$A^{-1}=q(A)$ の固有値は $q(\lambda_i)=\lambda_i^{-1}$ である。
:::

:::remark
Cayley–Hamilton の定理は、$K=\C$ では三角化を用いて次のようにも証明できる。$A=PUP^{-1}$, $U$ は対角成分 $\lambda_1,\dots,\lambda_n$ の上三角行列とする。$W_j\defeq\langle e_1,\dots,e_j\rangle$（$W_0=\{0\}$）とおくと、$U$ の上三角性から $UW_j\subset W_j$ であり、$(U-\lambda_jI_n)e_j\in W_{j-1}$ なので $(U-\lambda_jI_n)W_j\subset W_{j-1}$ である。したがって

$$
\Phi_A(U)=(U-\lambda_1I_n)(U-\lambda_2I_n)\cdots(U-\lambda_nI_n)
$$

を右の因子から順に $W_n=K^n$ に施すと、$W_n\to W_{n-1}\to\cdots\to W_0=\{0\}$ となるので $\Phi_A(U)=O$ であり、$\Phi_A(A)=P\Phi_A(U)P^{-1}=O$ である。本文の証明は余因子行列を用いることで、代数学の基本定理を使わずに任意の体で成り立つ。
:::

### 同時三角化

:::theorem[可換な行列の同時三角化]{#thm-linalg-simultaneous-triangularization}
$A,B\in M_n(K)$ が $AB=BA$ を満たし、$\Phi_A$ と $\Phi_B$ がともに $K$ 上で一次式の積に分解するとする。このとき、$P^{-1}AP$ と $P^{-1}BP$ がともに上三角行列となる正則行列 $P$ が存在する。特に、可換な複素正方行列は同時に三角化できる。
:::

:::proof
まず $A,B$ に**共通の固有ベクトル**があることを示す。$\lambda$ を $A$ の固有値とすると、$B$ は $A$ と可換なので $E\defeq E_\lambda(A)$ は $L_B$ 不変である（:ref[def-linalg-invariant-subspace] の直後の説明）。:ref[cor-linalg-invariant-charpoly] より、ある正方行列 $B'$ を用いて $\Phi_B=\Phi_{L_B|_E}\cdot\Phi_{B'}$ と書け、:ref[lem-linalg-factor-theorem] の 5 より $\Phi_{L_B|_E}$ も一次式の積に分解する。$E\ne\{0\}$ なので $\Phi_{L_B|_E}$ は次数 $1$ 以上であり、根 $\mu$ をもつ。よって $L_B|_E$ は固有ベクトル $w\in E$ をもち、$Aw=\lambda w$, $Bw=\mu w$ である。

$n$ に関する帰納法で定理を示す。$n=1$ は自明である。$n\ge2$ のとき、共通の固有ベクトル $w$ を $K^n$ の基底に延長して列に並べた行列を $Q_0$ とすると、:ref[prop-linalg-invariant-block] より

$$
Q_0^{-1}AQ_0=\begin{pmatrix}\lambda&h\\O&A_2\end{pmatrix},\qquad Q_0^{-1}BQ_0=\begin{pmatrix}\mu&h'\\O&B_2\end{pmatrix}
$$

となる。この 2 つの行列も可換であり、積の右下ブロックを比較すると $A_2B_2=B_2A_2$ である。また :ref[cor-linalg-invariant-charpoly] より $\Phi_A=(t-\lambda)\Phi_{A_2}$ であり、:ref[lem-linalg-factor-theorem] の 5 より $\Phi_{A_2}$ も一次式の積に分解し、$\Phi_{B_2}$ も同様である。帰納法の仮定により $Q^{-1}A_2Q$, $Q^{-1}B_2Q$ がともに上三角となる $Q$ があり、:ref[thm-linalg-triangularization] の証明と同じく $P=Q_0\begin{pmatrix}1&O\\O&Q\end{pmatrix}$ とすればよい。
:::

同時三角化を用いると、例えば可換な $A,B$ について $A+B$ や $AB$ の固有値は、適当な順序で並べた $A$ と $B$ の固有値の和や積 $\lambda_i+\mu_i$, $\lambda_i\mu_i$ であることが分かる（両者を同時に上三角にすれば、対角成分どうしが足され、掛けられる）。可換でない行列ではこれは一般に成り立たない。

## 演習問題

::::exercise{level=basic}
$A=\begin{pmatrix}4&-2\\1&1\end{pmatrix}$ の固有値と各固有空間を求め、$A$ を対角化する正則行列 $P$ を 1 つ与えよ。
:::solution
$\Phi_A(t)=(t-4)(t-1)+2=t^2-5t+6=(t-2)(t-3)$ である。$\lambda=2$ では $(A-2I_2)x=0$ が $2x_1-2x_2=0$ に、$\lambda=3$ では $(A-3I_2)x=0$ が $x_1-2x_2=0$ に帰着するので

$$
E_2(A)=\bigl\langle\transp{(1,1)}\bigr\rangle,\qquad
E_3(A)=\bigl\langle\transp{(2,1)}\bigr\rangle
$$

である。$P=\begin{pmatrix}1&2\\1&1\end{pmatrix}$ は $\det P=-1\ne0$ なので正則であり、列の順序に従って $AP=P\diag(2,3)$、すなわち $P^{-1}AP=\diag(2,3)$ である。トレース $4+1=5=2+3$、行列式 $4+2=6=2\cdot3$ も :ref[cor-linalg-trace-det-eigenvalues] と整合する。
:::
::::

::::exercise{level=basic}
$a\in K$ とし、$A=\begin{pmatrix}2&a&0\\0&2&0\\0&0&3\end{pmatrix}$ を考える。各固有値の 2 種類の重複度を求め、$A$ が対角化可能となる $a$ を決定せよ。
:::solution
上三角行列なので $\Phi_A(t)=(t-2)^2(t-3)$ であり、$m_2=2$, $m_3=1$ である。$(A-2I_3)x=0$ は $ax_2=0$, $x_3=0$ なので、$a=0$ ならば $g_2=2$、$a\ne0$ ならば $g_2=1$ である。単根 $3$ については常に $g_3=1$ である。:ref[thm-linalg-diagonalization-criterion] より、対角化可能であるのは $a=0$ の場合に限り、そのとき $A$ はすでに対角行列である。
:::
::::

::::exercise{level=basic}
$A=\begin{pmatrix}0&0&1\\1&0&1\\0&1&0\end{pmatrix}$ に対し、Cayley–Hamilton の定理を用いて $A^{-1}$ と $A^5$ を $I_3,A,A^2$ の一次結合で表せ。
:::solution
$tI_3-A$ の第 1 行に沿って展開すると

$$
\Phi_A(t)=\det\begin{pmatrix}t&0&-1\\-1&t&-1\\0&-1&t\end{pmatrix}=t(t^2-1)+(-1)\cdot(1-0)=t^3-t-1
$$

である。Cayley–Hamilton の定理より $A^3=A+I_3$ であり、これを $A(A^2-I_3)=I_3$ と書けば、$A^2-I_3$ は $A$ と可換なので $A^{-1}=A^2-I_3$ である。また $A^5=A^2A^3=A^3+A^2=A^2+A+I_3$ である。:ref[cor-linalg-ch-reduction] の 2 の公式（$c_0=-1$, $c_1=-1$, $c_2=0$）からも $A^{-1}=A^2-I_3$ を得る。
:::
::::

::::exercise{level=basic}
$A=\begin{pmatrix}5&-6&-6\\-1&4&2\\3&-6&-4\end{pmatrix}$ を対角化し、$A^k$（$k\ge0$）を求めよ。
:::solution
特性多項式を計算すると $\Phi_A(t)=t^3-5t^2+8t-4=(t-1)(t-2)^2$ である（$\tr A=5$, $\det A=4$ と整合する）。$A-2I_3=\begin{pmatrix}3&-6&-6\\-1&2&2\\3&-6&-6\end{pmatrix}$ の各行は $(1,-2,-2)$ の定数倍なので階数は $1$ であり、$g_2=2=m_2$ である。$E_2(A)=\{x\mid x_1-2x_2-2x_3=0\}$ の基底として $\transp{(2,1,0)},\transp{(2,0,1)}$ がとれる。$A-I_3=\begin{pmatrix}4&-6&-6\\-1&3&2\\3&-6&-5\end{pmatrix}$ を簡約化すると $E_1(A)=\bigl\langle\transp{(3,-1,3)}\bigr\rangle$ を得る。よって

$$
P=\begin{pmatrix}3&2&2\\-1&1&0\\3&0&1\end{pmatrix},\qquad P^{-1}=\begin{pmatrix}-1&2&2\\-1&3&2\\3&-6&-5\end{pmatrix},\qquad P^{-1}AP=\diag(1,2,2)
$$

である（$\det P=-1$）。$A^k=P\diag(1,2^k,2^k)P^{-1}$ を成分ごとに計算してもよいが、

$$
A^k=P\diag(1,0,0)P^{-1}+2^kP\diag(0,1,1)P^{-1}
$$

と分けるほうが見通しがよい。$E\defeq P\diag(1,0,0)P^{-1}$ とおくと $P\diag(0,1,1)P^{-1}=I_3-E$ である。$k=1$ の式 $A=E+2(I_3-E)$ から $E=2I_3-A$ が分かるので

$$
E=\begin{pmatrix}-3&6&6\\1&-2&-2\\-3&6&6\end{pmatrix},\qquad
A^k=E+2^k(I_3-E)=(2I_3-A)+2^k(A-I_3)
$$

である。$E$ は $P$ の第 1 列 $\transp{(3,-1,3)}$ と $P^{-1}$ の第 1 行 $(-1,2,2)$ の積であり、$E_2(A)$ に沿った $E_1(A)$ への射影（:ref[exm-linalg-projection]）である。$k=0$ で $I_3$ になることも確かめられる。
:::
::::

::::exercise{level=basic}
$A\in M_3(\C)$ の固有値が $1,-1,2$ であるとする。$\det A$, $\tr A$, $\det(A^2+A+I_3)$, $\tr(A^{-1})$ を求めよ。
:::solution
固有値が相異なる 3 個なので $\Phi_A(t)=(t-1)(t+1)(t-2)$ である。:ref[cor-linalg-trace-det-eigenvalues] より $\det A=1\cdot(-1)\cdot2=-2$, $\tr A=1-1+2=2$。:ref[cor-linalg-spectral-mapping] より $A^2+A+I_3$ の固有値は $\lambda^2+\lambda+1$ に $\lambda=1,-1,2$ を代入した $3,1,7$ なので、$\det(A^2+A+I_3)=21$。$A^{-1}$ の固有値は $1,-1,\tfrac12$ なので $\tr(A^{-1})=\tfrac12$。
:::
::::

::::exercise{level=advanced}
$\alpha,\beta\in K$ を相異なる元とし、線形変換 $T$ が $(T-\alpha\id_V)(T-\beta\id_V)=0$ を満たすとする。任意の $v\in V$ を $\Ker(T-\alpha\id_V)$ の元と $\Ker(T-\beta\id_V)$ の元の和に分解し、$T$ が対角化可能であることを示せ。
:::hint
$v$ に $T-\beta\id_V$ を施すと、どちらの核の元が得られるか。
:::
:::solution
$\alpha-\beta\ne0$ なので

$$
v_\alpha\defeq\frac{1}{\alpha-\beta}(T-\beta\id_V)v,
\qquad
v_\beta\defeq\frac{1}{\beta-\alpha}(T-\alpha\id_V)v
$$

とおける。$v_\alpha+v_\beta=\frac{1}{\alpha-\beta}\bigl((T-\beta\id_V)-(T-\alpha\id_V)\bigr)v=v$ である。仮定より $(T-\alpha\id_V)v_\alpha=0$ である。また 2 つの因子は $T$ の多項式なので可換であり（:ref[lem-linalg-polynomial-substitution]）、$(T-\beta\id_V)v_\beta=0$ でもある。両方の核に属する $w$ は $\alpha w=T(w)=\beta w$ を満たすので $w=0$ である。よって $V=\Ker(T-\alpha\id_V)\oplus\Ker(T-\beta\id_V)$ であり、各核の基底を並べれば固有ベクトルからなる基底を得る（片方の核が $\{0\}$ ならば、もう片方が $V$ 全体であり $T$ はスカラー倍である）。:ref[thm-linalg-diagonalization-criterion] より $T$ は対角化可能である。例えば射影（$T^2=T$, $\alpha=1$, $\beta=0$）や対合（$T^2=\id_V$, $\alpha=1$, $\beta=-1$、$1+1\ne0$ のとき）は対角化可能である。
:::
::::

::::exercise{level=advanced}
$A\in M_n(\R)$ が $A^2=-I_n$ を満たすとする。

1. $n$ は偶数であることを示せ。
2. $A$ は $\R$ 上で三角化できないことを示せ。
3. $A$ を $\C^n$ に作用する行列とみると対角化可能であることを示し、固有値 $i$ と $-i$ の重複度を求めよ。
:::solution
1. 行列式をとると $(\det A)^2=\det(-I_n)=(-1)^n$ である。$\det A$ は実数なので左辺は $0$ 以上であり、$(-1)^n=1$、すなわち $n$ は偶数である。

2. 実固有値 $\lambda$ と実固有ベクトル $v$ があれば、$-v=A^2v=\lambda^2v$ より $\lambda^2=-1$ となり矛盾する。$\R$ 上で $P^{-1}AP=U$ が上三角ならば、$Ue_1=u_{11}e_1$ より $Pe_1$ は実固有ベクトルになるので、三角化できない。

3. $(A-iI_n)(A+iI_n)=A^2+I_n=O$ なので、前問を $\alpha=i$, $\beta=-i$ に適用すれば $\C^n=E_i(A)\oplus E_{-i}(A)$ であり、$A$ は対角化可能である。:ref[prop-linalg-real-conjugate] より $v\mapsto\conj v$ は $E_i(A)$ から $E_{-i}(A)$ への全単射であり、一次独立なベクトルを一次独立なベクトルに移す（$\sum c_j\conj{v_j}=0\Rightarrow\sum\conj{c_j}v_j=0$）ので、$\dim E_i(A)=\dim E_{-i}(A)=n/2$ である。対角化可能なので、代数的重複度もともに $n/2$ である。
:::
::::

::::exercise{level=advanced}
1. $A,B\in M_n(K)$ に対し $\Phi_{AB}=\Phi_{BA}$ を示せ。
2. $A\in M_{m,n}(K)$, $B\in M_{n,m}(K)$ に対し $t^n\Phi_{AB}(t)=t^m\Phi_{BA}(t)$ を示せ。特に $AB$ と $BA$ の $0$ でない固有値は（重複度も込めて）一致する。
:::hint
$K[t]$ 上のブロック行列 $\begin{pmatrix}tI_m&A\\B&I_n\end{pmatrix}$ に、左右から適当なブロック三角行列を掛ける。
:::
:::solution
2 を示せば、$m=n$ の場合に $K[t]$ で $t^n$ を約分して（:ref[lem-linalg-polynomial-division] の 1）1 を得る。$M\defeq\begin{pmatrix}tI_m&A\\B&I_n\end{pmatrix}$, $N\defeq\begin{pmatrix}I_m&O\\-B&tI_n\end{pmatrix}$ を $K[t]$ 上の $m+n$ 次正方行列とする。ブロックの積を計算すると

$$
NM=\begin{pmatrix}tI_m&A\\O&tI_n-BA\end{pmatrix},\qquad
MN=\begin{pmatrix}tI_m-AB&tA\\O&tI_n\end{pmatrix}
$$

である。:ref[prop-linalg-det-ring] の 2 より、$\det N=t^n$ であり、$\det(NM)=t^m\Phi_{BA}(t)$, $\det(MN)=\Phi_{AB}(t)\,t^n$ である。積公式より $\det(NM)=\det N\det M=\det(MN)$ なので $t^n\Phi_{AB}=t^m\Phi_{BA}$。

$\lambda\ne0$ について、$t^n\Phi_{AB}=t^m\Phi_{BA}$ における $t-\lambda$ の重複度は左右で等しく、$t^n$, $t^m$ は $\lambda$ を根にもたないので、$\lambda$ の $\Phi_{AB}$ と $\Phi_{BA}$ における重複度は等しい（:ref[lem-linalg-factor-theorem] の 3）。例えば $u,v\in K^n$ に対し $u\transp{v}\in M_n(K)$ と $\transp{v}u\in M_1(K)$ を考えれば、$\Phi_{u\transp v}(t)=t^{n-1}(t-\transp{v}u)$ である。
:::
::::

::::exercise{level=advanced}
モニック多項式 $p(t)=t^n+c_{n-1}t^{n-1}+\cdots+c_1t+c_0\in K[t]$（$n\ge1$）に対し、

$$
C_p\defeq\begin{pmatrix}0&0&\cdots&0&-c_0\\1&0&\cdots&0&-c_1\\0&1&\cdots&0&-c_2\\ \vdots&&\ddots&&\vdots\\0&0&\cdots&1&-c_{n-1}\end{pmatrix}
$$

（$p$ の**同伴行列**）とおく。

1. $\Phi_{C_p}=p$ を示せ。したがって任意のモニック多項式はある行列の特性多項式である。
2. $C_p$ の各固有値の幾何学的重複度は $1$ であることを示せ。したがって $C_p$ が対角化可能であるための必要十分条件は、$p$ が $K$ 上で相異なる $n$ 個の一次式の積に分解することである。
:::solution
1. $n$ に関する帰納法で示す。$n=1$ ならば $C_p=(-c_0)$ で $\Phi_{C_p}=t+c_0=p$。$n\ge2$ とする。

$$
tI_n-C_p=\begin{pmatrix}t&0&\cdots&0&c_0\\-1&t&\cdots&0&c_1\\ \vdots&\ddots&\ddots&&\vdots\\0&\cdots&-1&t&c_{n-2}\\0&\cdots&0&-1&t+c_{n-1}\end{pmatrix}
$$

を第 1 行に沿って展開する。$(1,1)$ 成分 $t$ の余因子は、$p_1(t)\defeq t^{n-1}+c_{n-1}t^{n-2}+\cdots+c_1$ の同伴行列 $C_{p_1}$ について $\det(tI_{n-1}-C_{p_1})$ であり、帰納法の仮定より $p_1(t)$ に等しい。$(1,n)$ 成分 $c_0$ の余因子は $(-1)^{1+n}$ と、第 1 行と第 $n$ 列を除いた $(n-1)$ 次行列の行列式の積である。後者は対角成分がすべて $-1$ の上三角行列なので行列式は $(-1)^{n-1}$ である。よって

$$
\Phi_{C_p}(t)=t\,p_1(t)+c_0(-1)^{1+n}(-1)^{n-1}=t\,p_1(t)+c_0=p(t).
$$

2. 固有値 $\lambda$ について、$C_p-\lambda I_n$ の第 2–$n$ 行・第 1–$(n-1)$ 列からなる $n-1$ 次小行列式は、対角成分がすべて $1$ の上三角行列の行列式なので $1\ne0$ である。:ref[thm-linalg-rank-minor] より $\rank(C_p-\lambda I_n)\ge n-1$ であり、$g_\lambda=n-\rank(C_p-\lambda I_n)\le1$。固有値なので $g_\lambda=1$ である。:ref[thm-linalg-diagonalization-criterion] より、対角化可能であることは「$p$ が一次式の積に分解し、すべての根で $m_\lambda=g_\lambda=1$」、すなわち $p$ が相異なる $n$ 個の一次式の積に分解することと同値である。
:::
::::

::::exercise{level=advanced}
$A,B\in M_n(K)$ がともに対角化可能で $AB=BA$ を満たすならば、$P^{-1}AP$ と $P^{-1}BP$ がともに対角行列となる正則行列 $P$ が存在する（**同時対角化**）ことを示せ。逆に、同時対角化できる 2 つの行列は可換であることも示せ。
:::hint
まず「$B$ が対角化可能で $W$ が $L_B$ 不変ならば、$W$ は $B$ の固有ベクトルからなる基底をもつ」ことを示す。
:::
:::solution
**補題**：$B$ が対角化可能で、部分空間 $W$ が $L_B$ 不変ならば、$W=\bigoplus_j(W\cap E_{\mu_j}(B))$ である（$\mu_1,\dots,\mu_s$ は $B$ の相異なる固有値）。実際、$w\in W$ を $w=\sum_ju_j$（$u_j\in E_{\mu_j}(B)$）と分解する（:ref[thm-linalg-diagonalization-criterion] の (3)）。$q_i(t)\defeq\prod_{j\ne i}(t-\mu_j)$ とおくと、:ref[lem-linalg-polynomial-substitution] の 2 より $q_i(B)u_j=q_i(\mu_j)u_j$ であり、$q_i(\mu_j)=0$（$j\ne i$）なので $q_i(B)w=q_i(\mu_i)u_i$ である。$W$ は $L_B$ 不変なので $q_i(B)w\in W$ であり、$q_i(\mu_i)\ne0$ より $u_i\in W$ である。よって $W=\sum_j(W\cap E_{\mu_j}(B))$ であり、この和は直和である（:ref[prop-linalg-eigenspaces-direct-sum]）。

**本題**：$A$ の相異なる固有値を $\lambda_1,\dots,\lambda_r$ とすると $K^n=\bigoplus_iE_{\lambda_i}(A)$ である。$AB=BA$ より各 $E_{\lambda_i}(A)$ は $L_B$ 不変なので、補題により $E_{\lambda_i}(A)$ は $B$ の固有ベクトルからなる基底をもつ。その元は $A$ の固有ベクトルでもある。これらの基底を $i$ について並べると $K^n$ の基底になり、$A$ と $B$ の共通の固有ベクトルからなる。それを列に並べた行列を $P$ とすれば、$P^{-1}AP$ と $P^{-1}BP$ はともに対角行列である。

逆に $P^{-1}AP=D$, $P^{-1}BP=D'$ が対角行列ならば、対角行列どうしは可換なので $AB=PDD'P^{-1}=PD'DP^{-1}=BA$ である。
:::
::::

::::exercise{level=advanced}
$A\in M_n(\C)$ について、次の 3 条件が同値であることを示せ。

1. $A$ は**冪零**である（ある $k\ge1$ について $A^k=O$）。
2. $\Phi_A(t)=t^n$。
3. $A$ の固有値は $0$ だけである。

また、冪零行列が対角化可能ならば $A=O$ であることを示せ。
:::solution
(1)$\Rightarrow$(3)：$Av=\lambda v$, $v\ne0$ ならば $0=A^kv=\lambda^kv$ より $\lambda=0$。

(3)$\Rightarrow$(2)：$K=\C$ では $\Phi_A$ は一次式の積に分解し（代数学の基本定理）、その根は固有値なので、:ref[lem-linalg-factor-theorem] の 4 より $\Phi_A=\prod_{i=1}^n(t-0)=t^n$。

(2)$\Rightarrow$(1)：Cayley–Hamilton の定理より $A^n=\Phi_A(A)=O$。

最後に、冪零な $A$ が対角化可能ならば $P^{-1}AP=\diag(\lambda_1,\dots,\lambda_n)$ の対角成分は固有値なのですべて $0$ であり、$A=POP^{-1}=O$ である。例えば $\begin{pmatrix}0&1\\0&0\end{pmatrix}$ は冪零だが $O$ でないので、対角化できない。
:::
::::

## 確認テスト

::::quiz
$0$ が正方行列 $A\in M_n(K)$ の固有値であることと同値な条件をすべて選べ。

- [x] $Ax=0$ が自明でない解をもつ
- [ ] $A=O$
- [x] $\det A=0$
- [x] $\rank A<n$

:::explanation
$0$ が固有値であることは定義から $Ax=0\cdot x=0$ が自明でない解をもつことであり、:ref[cor-linalg-det-invertible] より $\det A=0$ と、:ref[thm-linalg-invertible-equiv] より $\rank A<n$ と同値である。$A=O$ は十分条件だが必要条件ではない（例えば $\diag(0,1)$）。
:::
::::

::::quiz{answer="2"}
$4$ 次正方行列 $A$ が $\rank(A-3I_4)=2$ を満たすとき、固有値 $3$ の幾何学的重複度はいくつか。

:::explanation
幾何学的重複度は $E_3(A)=\Ker L_{A-3I_4}$ の次元であり、:ref[cor-linalg-solution-space-dim] より $4-\rank(A-3I_4)=4-2=2$ である。
:::
::::

::::quiz
$3$ 次実正方行列 $A$ の特性多項式が $(t-1)^2(t+2)$ である。$A$ が $\R$ 上対角化可能であるための必要十分条件はどれか。

- [ ] $E_1(A)$ が 1 次元である
- [x] $E_1(A)$ が 2 次元である
- [ ] $\det A\ne0$ である
- [ ] $E_{-2}(A)$ が 1 次元である

:::explanation
単根 $-2$ の固有空間は :ref[prop-linalg-multiplicity-bound] により常に 1 次元なので、条件にならない。重根 $1$ の代数的重複度は $2$ なので、:ref[thm-linalg-diagonalization-criterion] より $\dim E_1(A)=2$ が必要十分である。$\det A=1\cdot1\cdot(-2)=-2\ne0$ は常に成り立つ。
:::
::::

::::quiz{answer="5"}
$2$ 次正方行列 $A$ の特性多項式が $\Phi_A(t)=t^2-5t+6$ である。Cayley–Hamilton の定理から得られる等式 $A^2=cA-6I_2$ の係数 $c$ を答えよ。

:::explanation
$\Phi_A(A)=A^2-5A+6I_2=O$ を移項すれば $A^2=5A-6I_2$ なので $c=5$ である。$c=\tr A$ でもある（:ref[prop-linalg-charpoly-coefficients]）。
:::
::::

::::quiz{answer="-2"}
$3$ 次正方行列 $A$ の固有値が重複度も込めて $1,2,2$ であるとき、$\det(A-3I_3)$ を求めよ。

:::explanation
:ref[cor-linalg-spectral-mapping] より $A-3I_3$ の固有値は $1-3,2-3,2-3=-2,-1,-1$ であり、:ref[cor-linalg-trace-det-eigenvalues] より $\det(A-3I_3)=(-2)(-1)(-1)=-2$ である。あるいは $\det(A-3I_3)=(-1)^3\Phi_A(3)=-(3-1)(3-2)^2=-2$ と計算してもよい。
:::
::::

::::quiz
すべての複素正方行列について成り立つものをすべて選べ。

- [x] 上三角行列と相似である
- [ ] 対角行列と相似である
- [x] 自分自身の特性多項式に代入すると零行列になる
- [x] トレースは固有値（重複度込み）の和に等しい

:::explanation
複素数上では特性多項式が一次式の積に分解するので、:ref[thm-linalg-triangularization] により三角化できる。:ref[thm-linalg-cayley-hamilton] は対角化可能性によらず成り立ち、:ref[cor-linalg-trace-det-eigenvalues] もすべての複素行列について成り立つ。対角化については、$\begin{pmatrix}0&1\\0&0\end{pmatrix}$ のように幾何学的重複度が代数的重複度より小さい行列が反例である。
:::
::::

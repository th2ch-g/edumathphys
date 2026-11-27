---
title: 線形写像と表現行列
description: 線形写像の核と像、次元定理、基底を通した線形写像と行列の対応（表現行列）、基底変換・行列の同値と相似
notation:
  - symbol: '\Ker f'
    meaning: 線形写像 f の核
  - symbol: '\Img f'
    meaning: 線形写像 f の像
  - symbol: '\rank f'
    meaning: 線形写像 f の階数（像の次元）
  - symbol: 'V\cong W'
    meaning: V と W が同型であること
  - symbol: '\Hom(V,W)'
    meaning: V から W への線形写像全体のなすベクトル空間
  - symbol: '\End(V)'
    meaning: V 上の線形変換全体 Hom(V,V)
  - symbol: 'L_A'
    meaning: 行列 A を左から掛ける線形写像 x ↦ Ax
  - symbol: '[v]_{\mathcal B}'
    meaning: 基底 B に関する v の座標ベクトル
  - symbol: '[f]_{\mathcal C\leftarrow\mathcal B}'
    meaning: 定義域の基底 B と終域の基底 C に関する f の表現行列
  - symbol: 'P_{\mathcal B\leftarrow\mathcal B^{\prime}}'
    meaning: 基底 B′ に関する座標を基底 B に関する座標に変える基底変換行列
  - symbol: '\tr A'
    meaning: 正方行列 A のトレース（対角成分の和）
---

## この章の目標

- 線形写像を定義し、基底上の値を指定すれば線形写像がただ 1 つに定まることを理解する。
- 核と像によって単射性・全射性を判定し、次元定理（階数と退化次数の和が定義域の次元に等しいこと）を証明する。
- 基底を選ぶことで線形写像と行列が 1 対 1 に対応し（表現行列）、合成が行列の積に対応することを示す。
- 行列の階数が列空間の次元にも行空間の次元にも等しいことを示し、階数標準形を導く。
- 基底変換による表現行列の変化を計算し、相似な行列が共有する量（階数・トレース）を理解する。

前章 [ベクトル空間・基底・次元](/math/linear-algebra/vector-spaces/) ではベクトル空間そのものを調べた。本章では、ベクトル空間の間の「和とスカラー倍を保つ写像」である線形写像を調べる。行列 $A$ を掛ける操作 $x\mapsto Ax$ は和とスカラー倍を保つ（:ref[prop-linalg-product-laws]）。この 2 つの性質だけを取り出したものが線形写像であり、多項式の微分や行列の転置のように数ベクトルでない対象に対する操作も、同じ枠組みで扱える。逆に、有限次元空間に基底を選べば、どの線形写像も行列で表される。この対応によって、行列の計算は線形写像の性質として、線形写像の問題は行列の計算として理解できるようになる。

本章の定義と定理はすべて任意の体 $K$ 上で成り立つ（例で $K=\R$ または $\C$ を仮定する箇所はその旨を記す）。

## 線形写像の定義と例

### 線形写像

:::definition[線形写像]{#def-linalg-linear-map}
$V, W$ を同じ体 $K$ 上のベクトル空間とする。写像 $f:V\to W$ が、任意の $u,v\in V$ と $a,b\in K$ について

$$
f(au+bv)=af(u)+bf(v)
$$

を満たすとき、$f$ を（$K$ 上の）:term[線形写像]{en="linear map" yomi="せんけいしゃぞう"} という。$V=W$ のとき、すなわち線形写像 $f:V\to V$ を $V$ 上の :term[線形変換]{en="linear transformation" yomi="せんけいへんかん"}（線形作用素）という。
:::

:::proposition[線形写像の基本性質]{#prop-linalg-linear-map-basic}
$f:V\to W$ を写像とする。

1. $f$ が線形であるための必要十分条件は、**加法性** $f(u+v)=f(u)+f(v)$（$u,v\in V$）と**斉次性** $f(av)=af(v)$（$a\in K$, $v\in V$）がともに成り立つことである。
2. $f$ が線形ならば $f(0)=0$, $f(-v)=-f(v)$ であり、任意の $v_1,\dots,v_k\in V$ と $a_1,\dots,a_k\in K$ について

$$
f\Bigl(\sum_{j=1}^k a_jv_j\Bigr)=\sum_{j=1}^k a_jf(v_j) \label{eq-linalg-linear-combination}
$$

が成り立つ。すなわち線形写像は一次結合を、同じ係数の一次結合に移す。
:::

:::proof
1. $f$ が線形ならば、$a=b=1$ とおいて加法性を得る。また $b=0$ とおくと、:ref[prop-linalg-vs-basic] より $0\cdot f(v)=0$ なので $f(au)=af(u)$ となり、斉次性を得る。逆に 2 条件が成り立てば $f(au+bv)=f(au)+f(bv)=af(u)+bf(v)$ である。

2. 斉次性で $a=0$ とすれば $f(0)=f(0\cdot0)=0\cdot f(0)=0$、$a=-1$ とすれば :ref[prop-linalg-vs-basic] より $f(-v)=f((-1)v)=(-1)f(v)=-f(v)$ である。最後の等式は $k$ に関する帰納法で示す。$k=1$ のときは斉次性そのものである。$k-1$ 個で成り立つとすると、加法性と斉次性から

$$
f\Bigl(\sum_{j=1}^k a_jv_j\Bigr)=f\Bigl(\sum_{j=1}^{k-1}a_jv_j\Bigr)+f(a_kv_k)=\sum_{j=1}^{k-1}a_jf(v_j)+a_kf(v_k)
$$

となる。$k=0$（空和）の場合も、両辺が $0$ なので $f(0)=0$ により成り立つ。
:::

### 線形写像の例

:::example[線形写像の例]{#exm-linalg-linear-maps}
1. **零写像・恒等写像・スカラー倍**：すべての $v$ を $0$ に移す**零写像** $0:V\to W$ と、恒等写像 $\id_V:V\to V$ は線形である。$c\in K$ に対し $c\,\id_V:v\mapsto cv$ も線形である。実際 $c(au+bv)=(ca)u+(cb)v=a(cu)+b(cv)$ であり、ここで体の乗法の可換性 $ca=ac$ を用いた。
2. **行列の定める写像**：$A\in M_{m,n}(K)$ に対し、$L_A:K^n\to K^m$ を $L_A(x)\defeq Ax$ で定める。:ref[prop-linalg-product-laws] より $A(au+bv)=a\,Au+b\,Av$ なので $L_A$ は線形である。例えば

$$
A=\begin{pmatrix}1&0&1\\0&1&1\end{pmatrix}
\quad\Longrightarrow\quad
L_A\begin{pmatrix}x\\y\\z\end{pmatrix}=\begin{pmatrix}x+z\\y+z\end{pmatrix}.
$$

3. **転置とトレース**：転置 $M_{m,n}(K)\to M_{n,m}(K)$, $X\mapsto\transp{X}$ は :ref[prop-linalg-transpose-laws] により線形である。正方行列 $X=(x_{ij})\in M_n(K)$ に対角成分の和

$$
\tr X\defeq\sum_{i=1}^n x_{ii}
$$

を対応させる写像 $\tr:M_n(K)\to K$ も線形である。$\tr X$ を $X$ の :term[トレース]{en="trace" yomi="とれーす"}（跡）という。
4. **多項式の微分・代入・掛け算**：$K[t]$ 上の（形式的な）**微分** $D$ を

$$
D\Bigl(\sum_{k=0}^d c_kt^k\Bigr)\defeq\sum_{k=1}^d kc_kt^{k-1}
$$

で定める（$kc_k$ は $c_k$ を $k$ 回加えたもの）。$D$ は「$k$ 次の係数を $k$ 倍して $k-1$ 次の係数にする」操作であり、多項式の和とスカラー倍は係数ごとに計算されるから、$D$ は線形変換である。$D$ は $K[t]_{\le n}$ を $K[t]_{\le n-1}$ に移す。また $c\in K$ での**代入** $\mathrm{ev}_c:K[t]\to K$, $p\mapsto p(c)$ と、**$t$ を掛ける写像** $M_t:K[t]\to K[t]$, $p\mapsto tp$ も線形である。
5. **関数空間上の写像**：$K=\R$ とする。連続関数にその積分を対応させる $C([a,b])\to\R$, $g\mapsto\int_a^b g(s)\,\dd s$ や、$C^1$ 級関数に導関数を対応させる $C^1([a,b])\to C([a,b])$, $g\mapsto g'$ は、積分と微分の線形性（[解析学](/math/analysis/)）により線形写像である。
:::

:::example[線形でない写像]{#exm-linalg-nonlinear-maps}
1. $A\in M_{m,n}(K)$, $b\in K^m$, $b\ne0$ とする。$x\mapsto Ax+b$ は $0$ を $b\ne0$ に移すので、:ref[prop-linalg-linear-map-basic] の 2 により線形でない。このような写像を**アフィン写像**という。
2. $\R\to\R$, $x\mapsto x^2$ は $f(2\cdot1)=4\ne 2=2f(1)$ なので線形でない。
3. 複素共役 $\C\to\C$, $z\mapsto\conj{z}$ は和を保ち、実数 $r$ に対して $\conj{rz}=r\conj{z}$ を満たすので、$\C$ を $\R$ 上のベクトル空間とみれば線形である。しかし $\conj{i\cdot1}=-i\ne i=i\cdot\conj{1}$ なので $\C$ 上線形ではない。線形性は、スカラーの体を指定して初めて意味をもつ。
:::

### 数ベクトル空間の間の線形写像

$K^n$ から $K^m$ への線形写像は、すべて行列を掛ける写像である。

:::proposition[数ベクトル空間の間の線形写像]{#prop-linalg-matrix-map}
1. 任意の線形写像 $f:K^n\to K^m$ に対し、$f=L_A$ となる $A\in M_{m,n}(K)$ がただ 1 つ存在する。それは第 $j$ 列が $f(e_j)$ の行列 $A=(f(e_1)\ \cdots\ f(e_n))$ である。
2. 写像の和とスカラー倍を値ごとに $(f+g)(x)\defeq f(x)+g(x)$, $(cf)(x)\defeq cf(x)$ で定めると、$A,B\in M_{m,n}(K)$, $c\in K$ に対して $L_{A+B}=L_A+L_B$, $L_{cA}=cL_A$ である。また $A\in M_{l,m}(K)$, $B\in M_{m,n}(K)$ に対して $L_A\circ L_B=L_{AB}$ である。
:::

:::proof
1. $x=\transp{(x_1,\dots,x_n)}=\sum_jx_je_j$ と書けば、:eqref[eq-linalg-linear-combination] より $f(x)=\sum_jx_jf(e_j)$ である。一方、行列とベクトルの積は列の一次結合であるから（[連立一次方程式と行列](/math/linear-algebra/linear-equations/) の積の定義の直後の式）、$A=(f(e_1)\ \cdots\ f(e_n))$ に対して $Ax=\sum_jx_jf(e_j)$ である。よって $f=L_A$ である。$L_A=L_B$ ならば $Ae_j=Be_j$ であり、これは $A$ と $B$ の第 $j$ 列が等しいことを意味するから $A=B$ である。
2. :ref[prop-linalg-product-laws] より $(A+B)x=Ax+Bx$, $(cA)x=c(Ax)$, $(AB)x=A(Bx)$ である。
:::

:::example[平面の回転]{#exm-linalg-rotation}
$\theta\in\R$ に対し

$$
R_\theta\defeq\begin{pmatrix}\cos\theta&-\sin\theta\\ \sin\theta&\cos\theta\end{pmatrix}
$$

とおく。$x\in\R^2$ を極座標で $x=\transp{(r\cos\varphi,\ r\sin\varphi)}$ と書くと、三角関数の加法定理により

$$
R_\theta x=\begin{pmatrix}r(\cos\theta\cos\varphi-\sin\theta\sin\varphi)\\ r(\sin\theta\cos\varphi+\cos\theta\sin\varphi)\end{pmatrix}=\begin{pmatrix}r\cos(\varphi+\theta)\\ r\sin(\varphi+\theta)\end{pmatrix}
$$

となる。したがって $L_{R_\theta}$ は原点のまわりの角 $\theta$ の回転であり、回転は線形写像である。$R_\theta$ の第 1 列・第 2 列は、それぞれ $e_1, e_2$ を角 $\theta$ 回転したものである（:ref[prop-linalg-matrix-map] の 1）。さらに積を計算して加法定理を用いると $R_\alpha R_\beta=R_{\alpha+\beta}$ となる。:ref[prop-linalg-matrix-map] の 2 により、これは「角 $\beta$ の回転の後に角 $\alpha$ の回転を行うと角 $\alpha+\beta$ の回転になる」ことの行列による表現である。
:::

### 基底の像から写像を作る

有限次元空間上の線形写像は、基底の行き先だけで決まり、しかも基底の行き先は自由に選べる。

:::proposition[基底上で指定した写像の延長]{#prop-linalg-map-on-basis}
$\mathcal B=(v_1,\ldots,v_n)$ を $V$ の基底とする。任意の $w_1,\ldots,w_n\in W$ に対し、$f(v_j)=w_j$（$1\le j\le n$）を満たす線形写像 $f:V\to W$ がただ 1 つ存在する。
:::

:::proof
:ref[def-linalg-basis] と :ref[prop-linalg-independence-unique] より、各 $v\in V$ はただ一通りに $v=\sum_j x_jv_j$ と書ける。そこで

$$
f(v)\defeq \sum_{j=1}^n x_jw_j
$$

とおく。係数 $x_j$ が $v$ から一意に定まるので、これは写像 $f:V\to W$ を定める。$u=\sum_j y_jv_j$ とすれば $au+bv=\sum_j(ay_j+bx_j)v_j$ であり、これが $au+bv$ の（唯一の）表示であるから

$$
f(au+bv)=\sum_j(ay_j+bx_j)w_j=a\sum_jy_jw_j+b\sum_jx_jw_j=af(u)+bf(v)
$$

となる。よって $f$ は線形であり、$v_j$ の表示は $x_j=1$, 他の係数 $0$ なので $f(v_j)=w_j$ である。逆に $g$ が $g(v_j)=w_j$ を満たす線形写像ならば、:eqref[eq-linalg-linear-combination] より $g(\sum_jx_jv_j)=\sum_jx_jw_j=f(\sum_jx_jv_j)$ となり $g=f$ である。
:::

例えば $f:K^2\to K$ に $f(e_1)=2$, $f(e_2)=-1$ を指定すると、$f(\transp{(s,t)})=2s-t$ と決まる。基底でない（一次従属な）ベクトルに値を指定すると、矛盾が起こりうる。例えば $f(e_1)=2$ と指定した後に $f(2e_1)=5$ を要求すると、線形性が要求する値 $4$ と矛盾する。また生成系でないベクトルに値を指定しただけでは、線形写像は一意に定まらない。

:::example[射影]{#exm-linalg-projection}
$V=W_1\oplus W_2$ を直和分解とする（:ref[def-linalg-direct-sum]）。各 $v\in V$ は $v=w_1+w_2$（$w_i\in W_i$）とただ一通りに書けるので、$p(v)\defeq w_1$ と定めることができる。$v'=w_1'+w_2'$ も同様に分解すれば $av+bv'=(aw_1+bw_1')+(aw_2+bw_2')$ が $av+bv'$ の分解であるから、$p(av+bv')=aw_1+bw_1'=ap(v)+bp(v')$ となり、$p$ は線形変換である。$p$ を $W_2$ に沿った $W_1$ への :term[射影]{en="projection" yomi="しゃえい"} という。定義から

$$
p\circ p=p,\qquad p(V)=W_1,\qquad \{v\mid p(v)=0\}=W_2
$$

であり、$\id_V-p$ は $W_1$ に沿った $W_2$ への射影である。射影は $W_1$ だけでは決まらず、補空間 $W_2$ の選び方に依存する。例えば $\R^2=\langle e_1\rangle\oplus\langle e_2\rangle$ に関する射影は $\transp{(x,y)}\mapsto\transp{(x,0)}$ であるが、$\R^2=\langle e_1\rangle\oplus\langle\transp{(1,1)}\rangle$ に関する射影は、$\transp{(x,y)}=(x-y)e_1+y\transp{(1,1)}$ より $\transp{(x,y)}\mapsto\transp{(x-y,0)}$ である。
:::

### 線形写像の和・スカラー倍・合成

:::definition[線形写像の空間]{#def-linalg-hom}
$V$ から $W$ への線形写像全体の集合を $\Hom(V,W)$ と書く（$K$ を明示するときは $\Hom_K(V,W)$）。$f,g\in\Hom(V,W)$, $c\in K$ に対し、和とスカラー倍を値ごとに

$$
(f+g)(v)\defeq f(v)+g(v),\qquad (cf)(v)\defeq c\,f(v)
$$

で定める。また $\End(V)\defeq\Hom(V,V)$ と書き、$V$ 上の線形変換全体を表す。$T\in\End(V)$ と $k\ge1$ に対し $T^k\defeq T\circ\cdots\circ T$（$k$ 個の合成）、$T^0\defeq\id_V$ とする。
:::

:::proposition[Hom の構造]{#prop-linalg-hom-structure}
1. $\Hom(V,W)$ は上の演算によって $K$ 上のベクトル空間になる。
2. $f\in\Hom(U,V)$, $g\in\Hom(V,W)$ ならば $g\circ f\in\Hom(U,W)$ である。
3. 合成は各変数について線形である：$(g_1+g_2)\circ f=g_1\circ f+g_2\circ f$, $g\circ(f_1+f_2)=g\circ f_1+g\circ f_2$, $(cg)\circ f=c(g\circ f)=g\circ(cf)$。
:::

:::proof
1. $V$ から $W$ への写像全体 $\mathcal F(V,W)$ は、値ごとの和とスカラー倍によってベクトル空間になる（公理は $W$ での公理から値ごとに従う。:ref[exm-linalg-vector-spaces] の 4 と同じである）。$\Hom(V,W)$ がその部分空間であることを示せばよい。零写像は線形であり、$f,g$ が線形ならば

$$
(f+g)(au+bv)=af(u)+bf(v)+ag(u)+bg(v)=a(f+g)(u)+b(f+g)(v)
$$

となり、$cf$ についても同様である。

2. $g(f(au+bv))=g(af(u)+bf(v))=ag(f(u))+bg(f(v))$。

3. 第 1 式と $(cg)\circ f=c(g\circ f)$ は値ごとの定義から直ちに従う。第 2 式は $g(f_1(u)+f_2(u))=g(f_1(u))+g(f_2(u))$ と $g$ の加法性から、$g\circ(cf)=c(g\circ f)$ は $g$ の斉次性から従う。
:::

したがって $\End(V)$ には和・スカラー倍・合成（積）が定義され、合成は結合法則を満たすが、一般に可換ではない。

:::remark
$K=\R$ または $\C$ とし、$K[t]$ 上の微分 $D$ と $M_t$ を考える。$p=\sum_kc_kt^k$ に対し

$$
D(M_tp)=\sum_k(k+1)c_kt^k=\sum_kc_kt^k+\sum_kkc_kt^k=p+M_t(Dp)
$$

なので $D\circ M_t-M_t\circ D=\id$ である。一方、有限次元空間ではこのような関係は成り立たない（[連立一次方程式と行列](/math/linear-algebra/linear-equations/) の演習と、後述の :ref[lem-linalg-trace-commute] による）。量子力学の正準交換関係 $[\hat x,\hat p]=i\hbar$ を満たす作用素が無限次元空間に作用しなければならないのは、この理由による（[量子力学](/physics/quantum-mechanics/)）。
:::

## 核と像

### 核と像の定義

:::definition[核と像]{#def-linalg-kernel-image}
線形写像 $f:V\to W$ に対し、

$$
\Ker f\defeq\{v\in V\mid f(v)=0\},\qquad
\Img f\defeq\{f(v)\mid v\in V\}
$$

をそれぞれ $f$ の :term[核]{en="kernel" yomi="かく"}、:term[像]{en="image" yomi="ぞう"} という。
:::

:::proposition[核と像による判定]{#prop-linalg-kernel-injective}
$\Ker f$ は $V$ の部分空間であり、$\Img f$ は $W$ の部分空間である。また

$$
f\text{ が単射}\ \iff\ \Ker f=\{0\},\qquad
f\text{ が全射}\ \iff\ \Img f=W
$$

が成り立つ。
:::

:::proof
$f(0)=0$ より $0\in\Ker f$, $0\in\Img f$ である。$u,v\in\Ker f$ ならば $f(au+bv)=af(u)+bf(v)=0$ なので $au+bv\in\Ker f$ である。また像の任意の 2 元 $f(u),f(v)$ について $af(u)+bf(v)=f(au+bv)\in\Img f$ である。これで :ref[def-linalg-subspace] の条件が確かめられた。

$f$ が単射ならば、$f(v)=0=f(0)$ から $v=0$ である。逆に $\Ker f=\{0\}$ とし、$f(u)=f(v)$ とすると $f(u-v)=f(u)-f(v)=0$ より $u-v\in\Ker f=\{0\}$、すなわち $u=v$ である。全射についての同値は全射の定義そのものである。
:::

核は「$f$ で潰れる方向」、像は「$f$ で到達できる範囲」を表す。単射性を調べるには、任意の 2 点 $u,v$ を比べる必要はなく、$0$ に移る点だけを調べればよい。これが線形写像の大きな利点である。

:::proposition[部分空間の像と逆像]{#prop-linalg-image-preimage}
$f:V\to W$ を線形写像とする。

1. $V$ の部分空間 $U$ の像 $f(U)=\{f(u)\mid u\in U\}$ は $W$ の部分空間である。特に $S\subset V$ に対して $f(\langle S\rangle)=\langle f(S)\rangle$ であり、$V=\langle S\rangle$ ならば $\Img f=\langle f(S)\rangle$ である。
2. $W$ の部分空間 $W'$ の逆像 $f^{-1}(W')=\{v\in V\mid f(v)\in W'\}$ は $\Ker f$ を含む $V$ の部分空間である。
:::

:::proof
1. $f(U)$ が部分空間であることは :ref[prop-linalg-kernel-injective] の像についての証明と同じである。$\langle S\rangle$ の元は $\sum_ic_is_i$（$s_i\in S$）の形であり、その像 $\sum_ic_if(s_i)$ は $f(S)$ の元の一次結合である。逆に $f(S)$ の元の一次結合 $\sum_ic_if(s_i)$ は $f(\sum_ic_is_i)\in f(\langle S\rangle)$ である。
2. $f(0)=0\in W'$ より $0\in f^{-1}(W')$。$f(u),f(v)\in W'$ ならば $f(au+bv)=af(u)+bf(v)\in W'$ である。$\Ker f=f^{-1}(\{0\})\subset f^{-1}(W')$ も明らかである。
:::

:::example[核と像の計算]{#exm-linalg-kernel-image}
:ref[exm-linalg-linear-maps] の $f=L_A:K^3\to K^2$, $f(\transp{(x,y,z)})=\transp{(x+z,y+z)}$ を考える。$f(\transp{(x,y,z)})=0$ は $x=y=-z$ と同値なので

$$
\Ker f=\bigl\langle\transp{(-1,-1,1)}\bigr\rangle
$$

である。一方、任意の $\transp{(s,t)}$ は $\transp{(s,t,0)}$ の像なので $\Img f=K^2$ である。特に $\transp{(2,-1)}$ に移る点をすべて求めると

$$
\begin{pmatrix}2\\-1\\0\end{pmatrix}+z\begin{pmatrix}-1\\-1\\1\end{pmatrix}\qquad(z\in K)
$$

となる。どの値を指定しても、その値に移る点の集合は核の方向に平行な直線である。
:::

この最後の観察は一般に成り立つ。

:::proposition[逆像の構造]{#prop-linalg-fiber}
$f:V\to W$ を線形写像、$w\in W$ とする。$w\notin\Img f$ ならば $f(v)=w$ を満たす $v$ は存在しない。$w\in\Img f$ ならば、$f(v_0)=w$ を満たす $v_0$ を 1 つ選ぶと

$$
\{v\in V\mid f(v)=w\}=v_0+\Ker f\defeq\{v_0+u\mid u\in\Ker f\}
$$

である。
:::

:::proof
前半は像の定義そのものである。$u\in\Ker f$ ならば $f(v_0+u)=w+0=w$ である。逆に $f(v)=w$ ならば $f(v-v_0)=0$ なので $u\defeq v-v_0\in\Ker f$ であり、$v=v_0+u$ である。
:::

$f=L_A$ に適用すれば、これは :ref[prop-linalg-solution-structure]（一般解 = 特殊解 + 斉次方程式の一般解）にほかならない。線形微分方程式など、「線形写像 $f$ について $f(v)=w$ を解く」形の問題はすべてこの構造をもつ。

### 単射・全射と基底

:::proposition[単射・全射と基底の像]{#prop-linalg-injective-basis}
$\mathcal B=(v_1,\dots,v_n)$ を $V$ の基底とし、$f:V\to W$ を線形写像とする。

1. $f$ が単射 $\iff$ $f(v_1),\dots,f(v_n)$ が一次独立。
2. $f$ が全射 $\iff$ $f(v_1),\dots,f(v_n)$ が $W$ を生成する。
3. $f$ が全単射 $\iff$ $(f(v_1),\dots,f(v_n))$ が $W$ の基底。

また $f$ が単射ならば、$V$ の任意の一次独立なベクトル $u_1,\dots,u_k$ について $f(u_1),\dots,f(u_k)$ は一次独立である。
:::

:::proof
最後の主張から示す。$f$ が単射で $\sum_ic_if(u_i)=0$ ならば $f(\sum_ic_iu_i)=0$ より $\sum_ic_iu_i\in\Ker f=\{0\}$ であり、$u_i$ の一次独立性から $c_i=0$ である。これで 1 の $\Rightarrow$ も示された。

1 の $\Leftarrow$：$v=\sum_jx_jv_j\in\Ker f$ とすると $\sum_jx_jf(v_j)=0$ であり、仮定より $x_j=0$、すなわち $v=0$ である。:ref[prop-linalg-kernel-injective] より $f$ は単射である。

2 は :ref[prop-linalg-image-preimage] の 1 より $\Img f=\langle f(v_1),\dots,f(v_n)\rangle$ であることから従う。3 は 1 と 2 を合わせたものである。
:::

### 同型写像

:::definition[同型写像]{#def-linalg-isomorphism}
全単射な線形写像 $f:V\to W$ を :term[同型写像]{en="isomorphism" yomi="どうけいしゃぞう"} という。$V$ から $W$ への同型写像が存在するとき、$V$ と $W$ は（$K$ 上）**同型**であるといい、$V\cong W$ と書く。$V$ から $V$ 自身への同型写像を $V$ の**自己同型**といい、その全体を $\GL(V)$ と書く。
:::

:::proposition[同型写像の性質]{#prop-linalg-isomorphism-inverse}
1. $f:V\to W$ が同型写像ならば、逆写像 $f^{-1}:W\to V$ も線形であり、同型写像である。
2. $f:U\to V$, $g:V\to W$ が同型写像ならば、$g\circ f$ も同型写像であり $(g\circ f)^{-1}=f^{-1}\circ g^{-1}$ である。
3. 同型 $\cong$ は同値関係である。すなわち $V\cong V$、$V\cong W\Rightarrow W\cong V$、$U\cong V$ かつ $V\cong W\Rightarrow U\cong W$ が成り立つ。
:::

:::proof
1. $y,z\in W$ をとる。$f$ は全射なので $y=f(u)$, $z=f(v)$ となる $u,v\in V$ があり、$f$ は単射なのでこれらは一意で $u=f^{-1}(y)$, $v=f^{-1}(z)$ である。$f(au+bv)=ay+bz$ なので

$$
f^{-1}(ay+bz)=au+bv=af^{-1}(y)+bf^{-1}(z)
$$

となる。$f^{-1}$ は全単射（逆写像が $f$）なので同型写像である。

2. 全単射の合成は全単射であり（[集合と論理](/math/sets-logic/)）、:ref[prop-linalg-hom-structure] の 2 より $g\circ f$ は線形である。$(f^{-1}\circ g^{-1})\circ(g\circ f)=\id_U$, $(g\circ f)\circ(f^{-1}\circ g^{-1})=\id_W$ である。

3. $\id_V$ は同型写像であり、対称律は 1、推移律は 2 による。
:::

特に $\GL(V)$ は合成と逆写像について閉じている（[群論](/math/group-theory/) の言葉では、合成について群をなす）。同型写像は一次独立性・生成性・基底・部分空間・次元をすべて保つ（:ref[prop-linalg-injective-basis]、:ref[prop-linalg-image-preimage]）。したがって線形代数の性質に関する限り、同型な空間は区別する必要がない。

## 次元定理

### 次元定理とその証明

:::theorem[次元定理]{#thm-linalg-rank-nullity}
$V$ を有限次元とし、$f:V\to W$ を線形写像とする。このとき $\Img f$ は有限次元であり、

$$
\dim V=\dim\Ker f+\dim\Img f \label{eq-linalg-rank-nullity}
$$

が成り立つ。$\dim\Img f$ を $f$ の**階数**といい $\rank f$ と書き、$\dim\Ker f$ を $f$ の :term[退化次数]{en="nullity" yomi="たいかじすう"} という。
:::

:::proof
$n=\dim V$ とおく。:ref[cor-linalg-subspace-dimension] より $\Ker f$ は有限次元で、$k\defeq\dim\Ker f\le n$ である。$\Ker f$ の基底 $(u_1,\ldots,u_k)$（$k=0$ のときは空の組）を :ref[thm-linalg-basis-extension] により $V$ の基底

$$
(u_1,\ldots,u_k,v_1,\ldots,v_{n-k})
$$

に延長する。$(f(v_1),\ldots,f(v_{n-k}))$ が $\Img f$ の基底であることを示す。

**生成すること**：$f(v)\in\Img f$ をとり、$v=\sum_ia_iu_i+\sum_jb_jv_j$ と展開すると、$f(u_i)=0$ より $f(v)=\sum_jb_jf(v_j)$ である。

**一次独立であること**：$\sum_jb_jf(v_j)=0$ とすると $f(\sum_jb_jv_j)=0$ なので $\sum_jb_jv_j\in\Ker f$ であり、$\sum_jb_jv_j=\sum_ia_iu_i$ と書ける。移項すると $\sum_ia_iu_i-\sum_jb_jv_j=0$ は $V$ の基底の一次関係なので、すべての $b_j$ が $0$ である。

よって $\Img f$ は $n-k$ 個のベクトルからなる基底をもち、$\dim\Img f=n-k$ である。$n-k=0$ のときは、空の組が $\Img f=\{0\}$ の基底である。
:::

次元定理は「入力の自由度 $\dim V$ が、潰される自由度 $\dim\Ker f$ と出力に残る自由度 $\dim\Img f$ に分配される」ことを述べている。証明の中の部分空間 $U=\langle v_1,\dots,v_{n-k}\rangle$ は $\Ker f$ の補空間である。より一般に次が成り立つ。

:::proposition[核の補空間と像]{#prop-linalg-complement-image}
$f:V\to W$ を線形写像とし、$U$ を $\Ker f$ の任意の補空間とする（$V=\Ker f\oplus U$）。このとき $f$ の $U$ への制限 $f|_U:U\to\Img f$ は同型写像である。
:::

:::proof
$f|_U$ の核は $U\cap\Ker f=\{0\}$ なので単射である。また $v\in V$ を $v=u_0+u$（$u_0\in\Ker f$, $u\in U$）と分解すれば $f(v)=f(u)$ なので、$f(U)=\Img f$ であり全射である。
:::

:::example[2 点での値]{#exm-linalg-evaluation}
$E:K[t]_{\le 2}\to K^2$ を $E(p)\defeq\transp{(p(0),p(1))}$ で定める。代入は線形なので $E$ は線形である。$p=a+bt+ct^2$ と書けば、核の条件は

$$
a=0,\qquad a+b+c=0
$$

となるので、$\Ker E=\langle t^2-t\rangle$, $\dim\Ker E=1$ である。$\dim K[t]_{\le2}=3$ だから、次元定理より $\dim\Img E=2$ であり、$\Img E=K^2$、すなわち $E$ は全射である。

実際、値 $s_0, s_1$ を指定すると、1 次式 $p_0(t)=s_0+(s_1-s_0)t$ がその値をとる。:ref[prop-linalg-fiber] より、同じ値をとる 2 次以下の多項式は、ちょうど $p_0(t)+c(t^2-t)$（$c\in K$）である。
:::

### 次元と単射・全射

:::corollary[同じ有限次元の空間の間の線形写像]{#cor-linalg-equal-dimension-map}
$\dim V=\dim W=n<\infty$ のとき、線形写像 $f:V\to W$ について、単射であること、全射であること、同型写像であることは互いに同値である。
:::

:::proof
:eqref[eq-linalg-rank-nullity] より、$\dim\Ker f=0$ と $\dim\Img f=n$ は同値である。前者は $\Ker f=\{0\}$ と同値であり、後者は :ref[cor-linalg-subspace-dimension] より $\Img f=W$ と同値である。よって :ref[prop-linalg-kernel-injective] から単射と全射は同値であり、どちらかが成り立てば両方が成り立つので同型写像である。
:::

:::corollary[次元の大小と単射・全射]{#cor-linalg-map-dimension-bounds}
$V, W$ を有限次元とし、$f:V\to W$ を線形写像とする。

1. $\rank f\le\min\{\dim V,\dim W\}$ である。
2. $\dim V>\dim W$ ならば $f$ は単射でない。
3. $\dim V<\dim W$ ならば $f$ は全射でない。
:::

:::proof
1. :eqref[eq-linalg-rank-nullity] より $\rank f\le\dim V$、:ref[cor-linalg-subspace-dimension] より $\rank f\le\dim W$。2. $\dim\Ker f=\dim V-\rank f\ge\dim V-\dim W>0$。3. $\dim\Img f\le\dim V<\dim W$ なので $\Img f\ne W$。
:::

$f=L_A$（$A\in M_{m,n}(K)$, $m<n$）に 2 を適用すると、「方程式の本数が未知数の個数より少ない斉次方程式は自明でない解をもつ」こと（:ref[cor-linalg-homogeneous-nontrivial]）が再び得られる。ただし次元論そのものが :ref[cor-linalg-homogeneous-nontrivial] に基づいて構築されているので（:ref[lem-linalg-steinitz]）、これは独立な別証明ではなく、同じ事実の言い換えである。

:::example[Lagrange 補間]{#exm-linalg-lagrange-interpolation}
$c_0,c_1,\dots,c_n\in K$ を相異なる元とし、

$$
E:K[t]_{\le n}\to K^{n+1},\qquad E(p)\defeq\transp{\bigl(p(c_0),\dots,p(c_n)\bigr)}
$$

を考える。$i=0,\dots,n$ に対し

$$
\ell_i(t)\defeq\prod_{j\ne i}\frac{t-c_j}{c_i-c_j}
$$

とおく（$c_i\ne c_j$ なので分母は $0$ でない。$n=0$ のときは $\ell_0=1$）。$\ell_i$ は $n$ 次多項式であり、$\ell_i(c_k)=\delta_{ik}$ を満たす。よって任意の $y=\transp{(y_0,\dots,y_n)}$ について $E(\sum_iy_i\ell_i)=y$ となり、$E$ は全射である。$\dim K[t]_{\le n}=n+1=\dim K^{n+1}$ なので、:ref[cor-linalg-equal-dimension-map] により $E$ は同型写像である。これから次が従う。

1. （**Lagrange 補間**）任意の $y_0,\dots,y_n\in K$ に対し、$p(c_i)=y_i$（$0\le i\le n$）を満たす $n$ 次以下の多項式 $p$ がただ 1 つ存在し、それは $p=\sum_iy_i\ell_i$ である。
2. 次数 $n$ 以下の多項式が相異なる $n+1$ 個の点で $0$ になるならば、それは零多項式である（$E$ の単射性）。
3. $(\ell_0,\dots,\ell_n)$ は $K[t]_{\le n}$ の基底であり、この基底に関する $p$ の座標は $\transp{(p(c_0),\dots,p(c_n))}$ である。

存在（全射性）を具体的な構成で示せば、一意性（単射性）は次元の比較から自動的に従う。これが :ref[cor-linalg-equal-dimension-map] の典型的な使い方である。
:::

:::example[無限次元では単射と全射は同値でない]{#exm-linalg-infinite-dim-maps}
$K=\R$ または $\C$ とし、$K[t]$ 上の微分 $D$ と、$J(\sum_kc_kt^k)\defeq\sum_k\frac{c_k}{k+1}t^{k+1}$ で定まる線形変換 $J$（積分）を考える。

- $D\circ J=\id$ なので $D$ は全射である。しかし $D(1)=0$ なので $D$ は単射でない。実際 $\Ker D=K$（定数多項式全体）である（$kc_k=0$, $k\ge1$ から $c_k=0$）。
- $J$ は単射（$c_k/(k+1)=0\Rightarrow c_k=0$）だが、$J$ の像の元は定数項が $0$ なので全射でない。$M_t$ も同様に単射だが全射でない。
- $D\circ J=\id$ であるにもかかわらず $J\circ D\ne\id$（$J(D(1))=0$）である。

有限次元では、:ref[cor-linalg-equal-dimension-map] と後述の :ref[cor-linalg-invertible-map-matrix] により、$S\circ T=\id_V$ ならば $T\circ S=\id_V$ である（正方行列について :ref[thm-linalg-invertible-equiv] の (6), (7) が成り立つことに対応する）。無限次元ではこれが成り立たない。
:::

## 表現行列

### 座標写像

:::definition[座標写像]{#def-linalg-coordinate-map}
$\mathcal B=(v_1,\ldots,v_n)$ を $V$ の基底とする。$v\in V$ を $v=\sum_jx_jv_j$ と表したときの座標（[ベクトル空間・基底・次元](/math/linear-algebra/vector-spaces/) の座標の節）を

$$
[v]_{\mathcal B}\defeq\begin{pmatrix}x_1\\ \vdots\\ x_n\end{pmatrix}\in K^n
$$

と書き、写像 $\varphi_{\mathcal B}:V\to K^n$, $v\mapsto[v]_{\mathcal B}$ を基底 $\mathcal B$ に関する :term[座標写像]{en="coordinate map" yomi="ざひょうしゃぞう"} という。
:::

:::proposition[座標写像は同型]{#prop-linalg-coordinate-isomorphism}
座標写像 $\varphi_{\mathcal B}:V\to K^n$ は同型写像であり、$\varphi_{\mathcal B}(v_j)=e_j$、逆写像は $\varphi_{\mathcal B}^{-1}(x)=\sum_jx_jv_j$ である。
:::

:::proof
$\psi(x)\defeq\sum_jx_jv_j$ で $\psi:K^n\to V$ を定めると、$\psi$ は線形であり（:ref[prop-linalg-map-on-basis] で $e_j\mapsto v_j$ と指定した写像）、標準基底を基底 $\mathcal B$ に移すので、:ref[prop-linalg-injective-basis] の 3 により同型写像である。座標の定義から $\varphi_{\mathcal B}$ は $\psi$ の逆写像であり、:ref[prop-linalg-isomorphism-inverse] により同型写像である。
:::

:::theorem[有限次元空間の同型類]{#thm-linalg-iso-dimension}
$V$ を有限次元とする。$V\cong W$ であるための必要十分条件は、$W$ が有限次元で $\dim W=\dim V$ となることである。特に $\dim V=n$ ならば $V\cong K^n$ である。
:::

:::proof
$f:V\to W$ を同型写像とし、$(v_1,\dots,v_n)$ を $V$ の基底とすると、:ref[prop-linalg-injective-basis] の 3 により $(f(v_1),\dots,f(v_n))$ は $W$ の基底なので、$\dim W=n$ である。逆に $\dim V=\dim W=n$ ならば、それぞれの基底 $\mathcal B,\mathcal C$ をとって $\varphi_{\mathcal C}^{-1}\circ\varphi_{\mathcal B}:V\to K^n\to W$ が同型写像である（:ref[prop-linalg-coordinate-isomorphism], :ref[prop-linalg-isomorphism-inverse]）。
:::

:::remark
有限次元空間は次元だけで同型類が決まる。しかし同型写像 $V\to K^n$ は基底の選び方に依存し、一般に「標準的な」選び方は存在しない。例えば $V=\{x\in\R^3\mid x_1+x_2+x_3=0\}$ は $\R^2$ と同型だが、$V$ の基底を 1 組選ばなければ具体的な同型写像は書けない。基底によらない議論と、基底を選んで行列で計算する議論とを行き来することが、線形代数の基本的な方法である。
:::

### 表現行列

:::definition[表現行列]{#def-linalg-representation-matrix}
$f:V\to W$ を線形写像とし、$\mathcal B=(v_1,\ldots,v_n)$ と $\mathcal C=(w_1,\ldots,w_m)$ をそれぞれ $V$ と $W$ の基底とする。$f(v_j)=\sum_{i=1}^m a_{ij}w_i$ と展開したときの係数を並べた $m\times n$ 行列

$$
[f]_{\mathcal C\leftarrow\mathcal B}
\defeq
\bigl([f(v_1)]_{\mathcal C}\ \cdots\ [f(v_n)]_{\mathcal C}\bigr)
=(a_{ij})
$$

を、基底 $\mathcal B,\mathcal C$ に関する $f$ の :term[表現行列]{en="representation matrix" yomi="ひょうげんぎょうれつ"} という。$V=W$ で同じ基底を用いるときは $[f]_{\mathcal B}\defeq[f]_{\mathcal B\leftarrow\mathcal B}$ と書く。
:::

表現行列の**第 $j$ 列は第 $j$ 基底ベクトルの像の座標**である。矢印の右側が入力の基底、左側が出力の基底であり、列の個数は入力の次元に、行の個数は出力の次元に等しい。例えば $L_A:K^n\to K^m$ の標準基底に関する表現行列は $A$ 自身である（:ref[prop-linalg-matrix-map]）。

:::proposition[座標での作用]{#prop-linalg-coordinate-action}
$A=[f]_{\mathcal C\leftarrow\mathcal B}$ は、すべての $v\in V$ について

$$
[f(v)]_{\mathcal C}=A[v]_{\mathcal B} \label{eq-linalg-coordinate-action}
$$

を満たす。逆に、この式をすべての $v\in V$ について満たす $m\times n$ 行列は $A$ に限る。
:::

:::proof
$v=\sum_j x_jv_j$ ならば、線形性により

$$
f(v)=\sum_jx_jf(v_j)=\sum_jx_j\sum_ia_{ij}w_i=\sum_i\Bigl(\sum_ja_{ij}x_j\Bigr)w_i
$$

である。括弧内は $A[v]_{\mathcal B}$ の第 $i$ 成分なので式が従う。逆に行列 $A'$ がこの式を満たすとき、$v=v_j$ を代入すると $[v_j]_{\mathcal B}=e_j$ なので、$A'e_j$（$A'$ の第 $j$ 列）$=[f(v_j)]_{\mathcal C}$（$A$ の第 $j$ 列）である。
:::

:eqref[eq-linalg-coordinate-action] は $\varphi_{\mathcal C}\circ f=L_A\circ\varphi_{\mathcal B}$ と書ける。すなわち次の図式は**可換**である（どちらの経路をたどっても同じ写像になる）。

$$
\begin{array}{ccc}
V & \xrightarrow{\quad f\quad} & W\\
{\scriptstyle\varphi_{\mathcal B}}\big\downarrow & & \big\downarrow{\scriptstyle\varphi_{\mathcal C}}\\
K^n & \xrightarrow[\quad L_A\quad]{} & K^m
\end{array}
$$

座標写像は同型なので、$f$ についての問題は $L_A$ についての問題、すなわち行列の計算に翻訳できる。

:::example[微分の表現行列]{#exm-linalg-derivative-matrix}
$D:K[t]_{\le2}\to K[t]_{\le1}$ について、$\mathcal B=(1,t,t^2)$, $\mathcal C=(1,t)$ とする。$D(1)=0$, $D(t)=1$, $D(t^2)=2t$ なので

$$
[D]_{\mathcal C\leftarrow\mathcal B}
=\begin{pmatrix}0&1&0\\0&0&2\end{pmatrix}.
$$

例えば $p=3-t+2t^2$ の座標 $\transp{(3,-1,2)}$ にこの行列を掛けると $\transp{(-1,4)}$ を得る。出力の基底に戻すと $Dp=-1+4t$ である。
:::

:::example[転置の表現行列]{#exm-linalg-transpose-matrix}
$\tau:M_2(K)\to M_2(K)$, $X\mapsto\transp{X}$ を、基底 $\mathcal E=(E_{11},E_{12},E_{21},E_{22})$（行列単位、:ref[exm-linalg-dimensions]）で表す。$\tau(E_{11})=E_{11}$, $\tau(E_{12})=E_{21}$, $\tau(E_{21})=E_{12}$, $\tau(E_{22})=E_{22}$ なので

$$
[\tau]_{\mathcal E}=\begin{pmatrix}1&0&0&0\\0&0&1&0\\0&1&0&0\\0&0&0&1\end{pmatrix}.
$$

「行列に作用する線形変換」も、$M_2(K)\cong K^4$ を通して $4$ 次正方行列で表される。
:::

### 合成と積

:::proposition[合成の表現行列]{#prop-linalg-composition-matrix}
線形写像 $f:U\to V$, $g:V\to W$ と、$U,V,W$ の基底 $\mathcal A,\mathcal B,\mathcal C$ に対し

$$
[g\circ f]_{\mathcal C\leftarrow\mathcal A}
=[g]_{\mathcal C\leftarrow\mathcal B}
 [f]_{\mathcal B\leftarrow\mathcal A}. \label{eq-linalg-composition-matrix}
$$
:::

:::proof
:ref[prop-linalg-hom-structure] の 2 より $g\circ f$ は線形である。任意の $u\in U$ に対し、:eqref[eq-linalg-coordinate-action] を 2 回用いると

$$
[g(f(u))]_{\mathcal C}
=[g]_{\mathcal C\leftarrow\mathcal B}[f(u)]_{\mathcal B}
=[g]_{\mathcal C\leftarrow\mathcal B}[f]_{\mathcal B\leftarrow\mathcal A}[u]_{\mathcal A}
$$

となる。:ref[prop-linalg-coordinate-action] の一意性から結論を得る。
:::

行列の積の定義（:eqref[eq-linalg-product]）は、まさにこの式が成り立つように定められている。行列の積が結合的であるのは写像の合成が結合的であることの反映であり、可換でないのは写像の合成が順序に依存することの反映である。

例えば :ref[exm-linalg-derivative-matrix] の微分の後に代入 $\mathrm{ev}_1:K[t]_{\le1}\to K$ を施すとする。基底 $(1,t)$ と $K$ の基底 $(1)$ に関する $\mathrm{ev}_1$ の表現行列は $(1\ \ 1)$ なので、合成の表現行列は

$$
\begin{pmatrix}1&1\end{pmatrix}
\begin{pmatrix}0&1&0\\0&0&2\end{pmatrix}
=\begin{pmatrix}0&1&2\end{pmatrix}
$$

である。実際、$p=a+bt+ct^2$ に対し $p'(1)=b+2c$ となる。

:::proposition[線形写像の空間と行列の空間]{#prop-linalg-hom-matrix-iso}
$\dim V=n$, $\dim W=m$ とし、基底 $\mathcal B,\mathcal C$ を固定する。写像

$$
\Phi:\Hom(V,W)\to M_{m,n}(K),\qquad f\mapsto[f]_{\mathcal C\leftarrow\mathcal B}
$$

は同型写像である。特に $\dim\Hom(V,W)=mn$ である。
:::

:::proof
$[(f+g)(v_j)]_{\mathcal C}=[f(v_j)]_{\mathcal C}+[g(v_j)]_{\mathcal C}$ と $[(cf)(v_j)]_{\mathcal C}=c[f(v_j)]_{\mathcal C}$ は座標写像の線形性による。これは $\Phi(f+g)=\Phi(f)+\Phi(g)$, $\Phi(cf)=c\Phi(f)$ を列ごとに述べたものなので、$\Phi$ は線形である。任意の $A=(a_{ij})\in M_{m,n}(K)$ に対し、:ref[prop-linalg-map-on-basis] により $f(v_j)=\sum_ia_{ij}w_i$（$1\le j\le n$）を満たす線形写像 $f$ がただ 1 つ存在する。これは $\Phi(f)=A$ となる $f$ がただ 1 つ存在すること、すなわち $\Phi$ が全単射であることを意味する。最後の主張は :ref[thm-linalg-iso-dimension] と $\dim M_{m,n}(K)=mn$（:ref[exm-linalg-dimensions]）による。
:::

:::corollary[同型写像と正則行列]{#cor-linalg-invertible-map-matrix}
$\dim V=\dim W=n$ とし、$\mathcal B,\mathcal C$ をそれぞれ $V,W$ の基底とする。線形写像 $f:V\to W$ が同型写像であるための必要十分条件は $A=[f]_{\mathcal C\leftarrow\mathcal B}$ が正則であることであり、このとき $[f^{-1}]_{\mathcal B\leftarrow\mathcal C}=A^{-1}$ である。特に $L_A:K^n\to K^n$ が同型写像であることと $A$ が正則であることは同値である。また $T\in\End(V)$ と $\mathcal B$ について $[\id_V]_{\mathcal B}=I_n$, $[S\circ T]_{\mathcal B}=[S]_{\mathcal B}[T]_{\mathcal B}$ が成り立つ。
:::

:::proof
最後の主張は :ref[prop-linalg-composition-matrix] と $\id_V(v_j)=v_j$ による。$f$ が同型写像ならば、:eqref[eq-linalg-composition-matrix] より

$$
[f^{-1}]_{\mathcal B\leftarrow\mathcal C}A=[f^{-1}\circ f]_{\mathcal B\leftarrow\mathcal B}=[\id_V]_{\mathcal B}=I_n
$$

であり、同様に $A[f^{-1}]_{\mathcal B\leftarrow\mathcal C}=I_n$ である。逆に $A$ が正則ならば、:ref[prop-linalg-hom-matrix-iso] により $[g]_{\mathcal B\leftarrow\mathcal C}=A^{-1}$ となる線形写像 $g:W\to V$ がある。$[g\circ f]_{\mathcal B}=A^{-1}A=I_n=[\id_V]_{\mathcal B}$ なので、:ref[prop-linalg-hom-matrix-iso] の単射性から $g\circ f=\id_V$ であり、同様に $f\circ g=\id_W$ である。よって $f$ は全単射である。
:::

この系により、$\End(V)$ における和・スカラー倍・合成は、基底を 1 つ固定すれば $M_n(K)$ における和・スカラー倍・積とそっくり対応する。線形変換について成り立つことは行列について成り立ち、逆も成り立つ。

### 階数と列空間・行空間

:::theorem[行階数と列階数の一致]{#thm-linalg-row-column-rank}
$A\in M_{m,n}(K)$ の列が張る $K^m$ の部分空間を $A$ の :term[列空間]{en="column space" yomi="れつくうかん"}、行を転置した列が張る $K^n$ の部分空間を $A$ の :term[行空間]{en="row space" yomi="ぎょうくうかん"} という。このとき

$$
\dim\Img L_A=\dim(\text{列空間})=\rank A=\dim(\text{行空間})=\rank \transp{A}
$$

である。ここで $\rank A$ は :ref[def-linalg-rank] で定めた簡約行階段形の主成分の個数である。
:::

:::proof
$A$ の列を $a_1,\ldots,a_n$ とすると $Ax=\sum_j x_ja_j$ であり、$x$ の成分は任意に選べるので、$\Img L_A$ は列空間に等しい。

$A$ の簡約行階段形を $R$ とし、$R=EA$（$E$ は基本行列の積で正則、:ref[lem-linalg-elementary-left-mult]）と書く。$R$ の第 $j$ 列は $Ea_j$ なので、$L_E$ は $A$ の列空間を $R$ の列空間に移す。$L_E$ は単射なので（:ref[cor-linalg-invertible-map-matrix]）、この制限は同型写像であり、2 つの列空間の次元は等しい。$r=\rank A$ とおく。$R$ の主列は $e_1,\ldots,e_r$ であり、$R$ のどの列も第 $r+1$ 成分以下は $0$ なので $\langle e_1,\dots,e_r\rangle$ に属する。よって $(e_1,\dots,e_r)$ は $R$ の列空間の基底であり、$A$ の列空間も $r$ 次元である。

次に行空間を考える。$R=EA$ より $R$ の各行は $A$ の行の一次結合であり、$A=E^{-1}R$ より逆も成り立つ。したがって $A$ と $R$ の行空間は等しい。$R$ の零でない $r$ 個の行は一次独立である。実際、第 $i$ 行の主成分のある列（第 $j_i$ 列）では、第 $i$ 行の成分が $1$ で他の行の成分は $0$ なので、これらの行の一次結合が $0$ ならば第 $j_i$ 成分を見て第 $i$ 行の係数が $0$ となる。よって行空間の次元は $r$ である。

最後に、$\transp{A}$ の列空間は $A$ の行空間そのものなので、すでに示した「列空間の次元 = 階数」を $\transp{A}$ に適用すれば $\rank\transp{A}=\dim(A\text{ の行空間})=r$ を得る。
:::

この証明から、$A$ の列空間の基底は、$R$ の主列と**同じ番号の列を元の $A$ から**取り出せば得られる（$L_{E^{-1}}$ が $R$ の主列 $e_1,\dots,e_r$ を $A$ の対応する列に戻すから）。一方、行空間の基底としては $R$ 自身の零でない行（の転置）を使える。$R$ の列空間と $A$ の列空間は一般に異なるので、列空間の基底に $R$ の列を使ってはならない。

:::corollary[解空間の次元]{#cor-linalg-solution-space-dim}
$A\in M_{m,n}(K)$ に対し、斉次方程式 $Ax=0$ の解空間 $\Ker L_A$ の次元は $n-\rank A$ である。
:::

:::proof
$L_A$ に次元定理を適用し、:ref[thm-linalg-row-column-rank] の $\dim\Img L_A=\rank A$ を用いる。
:::

[連立一次方程式と行列](/math/linear-algebra/linear-equations/) の :ref[thm-linalg-rank-criterion] で「一般解は $n-\rank A$ 個の任意定数を含む」と述べたことの正確な意味は、解空間が $n-\rank A$ 次元のベクトル空間であるということである。

:::proposition[線形写像の階数と表現行列の階数]{#prop-linalg-rank-map-matrix}
有限次元空間の間の線形写像 $f:V\to W$ と基底 $\mathcal B,\mathcal C$ について、$\rank f=\rank[f]_{\mathcal C\leftarrow\mathcal B}$ である。
:::

:::proof
$A=[f]_{\mathcal C\leftarrow\mathcal B}$ とおく。:eqref[eq-linalg-coordinate-action] と $\varphi_{\mathcal B}(V)=K^n$ より $\varphi_{\mathcal C}(\Img f)=\{A[v]_{\mathcal B}\mid v\in V\}=\Img L_A$ である。$\varphi_{\mathcal C}$ は同型なので $\dim\Img f=\dim\Img L_A=\rank A$（:ref[thm-linalg-row-column-rank]）。
:::

:::proposition[積の階数]{#prop-linalg-rank-product}
1. $A\in M_{l,m}(K)$, $B\in M_{m,n}(K)$ ならば $\rank(AB)\le\min\{\rank A,\rank B\}$。
2. $A\in M_{m,n}(K)$ と正則行列 $Q\in M_m(K)$, $P\in M_n(K)$ に対し $\rank(QAP)=\rank A$。
:::

:::proof
まず一般に、線形写像 $g$ と部分空間 $U$ について $\dim g(U)\le\dim U$ であることに注意する（$g|_U$ に次元定理を適用すれば $\dim g(U)=\dim U-\dim(U\cap\Ker g)$）。

1. $\Img L_{AB}=L_A(\Img L_B)$ である。これは $\Img L_A$ に含まれるので $\rank(AB)\le\rank A$。また上の注意から $\dim L_A(\Img L_B)\le\dim\Img L_B=\rank B$。
2. 1 より $\rank(QAP)\le\rank A$ であり、$A=Q^{-1}(QAP)P^{-1}$ に再び 1 を用いれば $\rank A\le\rank(QAP)$。
:::

## 基底変換と行列の相似

### 基底変換行列

同じベクトル空間に 2 つの基底があるとき、同じベクトルの座標は基底ごとに異なる。その関係を表すのが基底変換行列である。

:::definition[基底変換行列]{#def-linalg-change-of-basis}
$V$ の 2 つの基底 $\mathcal B=(v_1,\ldots,v_n)$, $\mathcal B'=(v'_1,\ldots,v'_n)$ に対し、恒等写像の表現行列

$$
P_{\mathcal B\leftarrow\mathcal B'}
\defeq[\id_V]_{\mathcal B\leftarrow\mathcal B'}
=\bigl([v'_1]_{\mathcal B}\ \cdots\ [v'_n]_{\mathcal B}\bigr)
$$

を、$\mathcal B'$ から $\mathcal B$ への :term[基底変換行列]{en="change-of-basis matrix" yomi="きていへんかんぎょうれつ"} という。すなわち $v'_j=\sum_ip_{ij}v_i$ である。
:::

:::proposition[基底変換行列の性質]{#prop-linalg-change-of-basis-matrix}
$\mathcal B,\mathcal B',\mathcal B''$ を $V$ の基底とする。

1. すべての $v\in V$ について $[v]_{\mathcal B}=P_{\mathcal B\leftarrow\mathcal B'}[v]_{\mathcal B'}$。
2. $P_{\mathcal B\leftarrow\mathcal B''}=P_{\mathcal B\leftarrow\mathcal B'}P_{\mathcal B'\leftarrow\mathcal B''}$, $P_{\mathcal B\leftarrow\mathcal B}=I_n$。
3. $P_{\mathcal B\leftarrow\mathcal B'}$ は正則であり、$(P_{\mathcal B\leftarrow\mathcal B'})^{-1}=P_{\mathcal B'\leftarrow\mathcal B}$。
4. 逆に、任意の正則行列 $P=(p_{ij})\in M_n(K)$ に対し、$v'_j\defeq\sum_ip_{ij}v_i$ とおけば $\mathcal B'=(v'_1,\dots,v'_n)$ は $V$ の基底であり、$P_{\mathcal B\leftarrow\mathcal B'}=P$ である。
:::

:::proof
1 は :eqref[eq-linalg-coordinate-action] を $\id_V$ に適用したもの、2 は :eqref[eq-linalg-composition-matrix] を $\id_V\circ\id_V=\id_V$ に適用したものである。3 は 2 で $\mathcal B''=\mathcal B$ とした式と、$\mathcal B,\mathcal B'$ の役割を入れ替えた式から従う。

4. $v'_j=\varphi_{\mathcal B}^{-1}(p_j)$（$p_j$ は $P$ の第 $j$ 列）である。$P$ は正則なので $\rank P=n$（:ref[thm-linalg-invertible-equiv]）であり、:ref[prop-linalg-independence-rank] と :ref[cor-linalg-n-vectors-basis] により $(p_1,\dots,p_n)$ は $K^n$ の基底である。同型写像 $\varphi_{\mathcal B}^{-1}$ は基底を基底に移すので（:ref[prop-linalg-injective-basis]）、$\mathcal B'$ は $V$ の基底である。$[v'_j]_{\mathcal B}=p_j$ なので $P_{\mathcal B\leftarrow\mathcal B'}=P$ である。
:::

:::remark
$P_{\mathcal B\leftarrow\mathcal B'}$ の列は**新しい基底 $\mathcal B'$ のベクトルを古い基底 $\mathcal B$ で表した座標**であるが、この行列が変換するのは**新しい座標から古い座標へ**である。新しい座標を求めるには逆行列 $P_{\mathcal B'\leftarrow\mathcal B}=(P_{\mathcal B\leftarrow\mathcal B'})^{-1}$ が必要である。向きを取り違えやすいので、記号の矢印で確認するとよい。
:::

:::example[座標の変換]{#exm-linalg-coordinate-change}
$\R^2$ の標準基底 $\mathcal E=(e_1,e_2)$ と、基底 $\mathcal B'=(\transp{(1,1)},\transp{(1,2)})$ を考える。

$$
P=P_{\mathcal E\leftarrow\mathcal B'}=\begin{pmatrix}1&1\\1&2\end{pmatrix},\qquad
P^{-1}=\begin{pmatrix}2&-1\\-1&1\end{pmatrix}
$$

である。$v=\transp{(3,5)}$ の $\mathcal B'$ に関する座標は $[v]_{\mathcal B'}=P^{-1}\transp{(3,5)}=\transp{(1,2)}$ であり、実際 $1\cdot\transp{(1,1)}+2\cdot\transp{(1,2)}=\transp{(3,5)}$ である。
:::

### 表現行列の基底変換

:::theorem[表現行列の基底変換]{#thm-linalg-change-of-basis}
$f:V\to W$ を線形写像とし、$A=[f]_{\mathcal C\leftarrow\mathcal B}$ とする。$V$ の別の基底 $\mathcal B'$ と $W$ の別の基底 $\mathcal C'$ に対し、$P=P_{\mathcal B\leftarrow\mathcal B'}$, $Q=P_{\mathcal C\leftarrow\mathcal C'}$ とおけば

$$
[f]_{\mathcal C'\leftarrow\mathcal B'}=Q^{-1}AP \label{eq-linalg-change-of-basis}
$$

が成り立つ。
:::

:::proof
$f=\id_W\circ f\circ\id_V$ に :eqref[eq-linalg-composition-matrix] を適用すると

$$
[f]_{\mathcal C'\leftarrow\mathcal B'}=[\id_W]_{\mathcal C'\leftarrow\mathcal C}\,[f]_{\mathcal C\leftarrow\mathcal B}\,[\id_V]_{\mathcal B\leftarrow\mathcal B'}=P_{\mathcal C'\leftarrow\mathcal C}AP
$$

であり、:ref[prop-linalg-change-of-basis-matrix] の 3 より $P_{\mathcal C'\leftarrow\mathcal C}=Q^{-1}$ である。
:::

特に $f=L_A:K^n\to K^m$ で、$\mathcal B,\mathcal C$ が標準基底のとき、新しい基底のベクトルを列に並べた行列を $P$, $Q$ とすれば $[L_A]_{\mathcal C'\leftarrow\mathcal B'}=Q^{-1}AP$ である。

:::example[入力と出力で異なる基底変換]{#exm-linalg-change-of-basis-evaluation}
:ref[exm-linalg-evaluation] の $E(p)=\transp{(p(0),p(1))}$ を考える。入力側の基底 $(1,t,t^2)$ と出力側の標準基底に関する表現行列は

$$
A=\begin{pmatrix}1&0&0\\1&1&1\end{pmatrix}
$$

である。入力側を $(1,\ t,\ t^2-t)$ に、出力側を $(\transp{(1,1)},\transp{(0,1)})$ に変える。新しい基底を古い基底の座標で表して列に並べると

$$
P=\begin{pmatrix}1&0&0\\0&1&-1\\0&0&1\end{pmatrix},\qquad
Q=\begin{pmatrix}1&0\\1&1\end{pmatrix},\qquad
Q^{-1}AP=\begin{pmatrix}1&0&0\\0&1&0\end{pmatrix}
$$

となる。この結果は写像の値からも読み取れる。$E(1)=\transp{(1,1)}$ と $E(t)=\transp{(0,1)}$ は順に出力側の新しい基底ベクトルであり、$t^2-t$ は核に属するので第 3 列は $0$ である。
:::

### 行列の同値と階数標準形

上の例のように、入力と出力の基底を独立に選べるならば、表現行列は非常に簡単な形にできる。

:::definition[行列の同値]{#def-linalg-matrix-equivalence}
$A,B\in M_{m,n}(K)$ について、正則行列 $Q\in M_m(K)$, $P\in M_n(K)$ で $B=Q^{-1}AP$ となるものが存在するとき、$A$ と $B$ は**同値**であるという。
:::

同値は同値関係である（$P,Q$ を単位行列にとれば反射律、$A=QBP^{-1}=(Q^{-1})^{-1}BP^{-1}$ から対称律、$C=Q'^{-1}BP'$ ならば $C=(QQ')^{-1}A(PP')$ から推移律）。:ref[thm-linalg-change-of-basis] と :ref[prop-linalg-change-of-basis-matrix] の 4 により、$A,B$ が同値であることは、同じ線形写像 $L_A$ を適当な基底で表したものが $B$ であることと同値である。

:::theorem[階数標準形]{#thm-linalg-rank-normal-form}
$\dim V=n$, $\dim W=m$ とし、$f:V\to W$ を $\rank f=r$ の線形写像とする。このとき $V$ の基底 $\mathcal B$ と $W$ の基底 $\mathcal C$ をうまく選べば

$$
[f]_{\mathcal C\leftarrow\mathcal B}=\begin{pmatrix}I_r&O\\O&O\end{pmatrix}\qquad(m\times n\text{ 行列}) \label{eq-linalg-rank-normal-form}
$$

となる（$r=0$ のときは零行列、$r=m$ や $r=n$ のときは対応する零ブロックを省く）。
:::

:::proof
次元定理より $\dim\Ker f=n-r$ である。$\Ker f$ の基底 $(u_1,\dots,u_{n-r})$ を $V$ の基底に延長し、付け加えたベクトルを**先頭に**並べて $\mathcal B=(v_1,\dots,v_r,u_1,\dots,u_{n-r})$ とする。:ref[thm-linalg-rank-nullity] の証明により $(f(v_1),\dots,f(v_r))$ は $\Img f$ の基底である。これを $W$ の基底 $\mathcal C=(f(v_1),\dots,f(v_r),w_{r+1},\dots,w_m)$ に延長する。このとき $f(v_j)$ は $\mathcal C$ の第 $j$ ベクトル（$j\le r$）、$f(u_k)=0$ なので、表現行列は :eqref[eq-linalg-rank-normal-form] となる。
:::

:::corollary[同値類は階数で決まる]{#cor-linalg-equivalence-rank}
$A\in M_{m,n}(K)$ が $\rank A=r$ を満たすならば、$A$ は :eqref[eq-linalg-rank-normal-form] の行列 $\begin{pmatrix}I_r&O\\O&O\end{pmatrix}$ と同値である。したがって $A,B\in M_{m,n}(K)$ が同値であるための必要十分条件は $\rank A=\rank B$ である。
:::

:::proof
$L_A$ に :ref[thm-linalg-rank-normal-form] を適用すると、$K^n, K^m$ の基底 $\mathcal B',\mathcal C'$ で $[L_A]_{\mathcal C'\leftarrow\mathcal B'}$ が標準形になるものがある（$\rank L_A=\rank A$ は :ref[thm-linalg-row-column-rank] による）。:ref[thm-linalg-change-of-basis] より、これは標準基底に関する表現行列 $A$ と同値である。同値な行列の階数が等しいことは :ref[prop-linalg-rank-product] の 2 による。逆に階数が等しければ、両者は同じ標準形と同値なので、推移律により互いに同値である。
:::

:::example[階数標準形の計算]{#exm-linalg-rank-normal-form}
$A=\begin{pmatrix}1&2&1\\2&4&3\end{pmatrix}$ を考える。行基本変形（第 2 行から第 1 行の 2 倍を引き、次に第 1 行から第 2 行を引く）で

$$
A\longrightarrow\begin{pmatrix}1&2&1\\0&0&1\end{pmatrix}\longrightarrow\begin{pmatrix}1&2&0\\0&0&1\end{pmatrix}=EA,\qquad E=\begin{pmatrix}3&-1\\-2&1\end{pmatrix}
$$

となる。さらに列基本変形（第 2 列から第 1 列の 2 倍を引き、次に第 2 列と第 3 列を入れ替える）を施す。列基本変形は、:ref[lem-linalg-elementary-left-mult] を転置して分かるとおり、基本行列を**右から**掛けることで実現される。

$$
EA\longrightarrow\begin{pmatrix}1&0&0\\0&0&1\end{pmatrix}\longrightarrow\begin{pmatrix}1&0&0\\0&1&0\end{pmatrix}=EAP,\qquad P=\begin{pmatrix}1&0&-2\\0&0&1\\0&1&0\end{pmatrix}.
$$

よって $Q\defeq E^{-1}=\begin{pmatrix}1&1\\2&3\end{pmatrix}$ とおけば $Q^{-1}AP=\begin{pmatrix}1&0&0\\0&1&0\end{pmatrix}$ であり、$\rank A=2$ である。一般に、行基本変形で簡約行階段形にした後、列の入れ替えで主列を左に集め、列基本変形で主成分の右の成分を消せば、階数標準形が得られる。
:::

### 行列の相似

線形変換 $T:V\to V$ を表すときは、入力と出力に**同じ基底**を用いるのが自然である（例えば $T$ を繰り返し施すとき）。このとき :eqref[eq-linalg-change-of-basis] で $Q=P$ となる。

:::definition[相似]{#def-linalg-similarity}
$A,B\in M_n(K)$ に対して、正則行列 $P\in M_n(K)$ で $B=P^{-1}AP$ となるものが存在するとき、$A$ と $B$ は :term[相似]{en="similar" yomi="そうじ"} であるという。
:::

:::proposition[相似と線形変換]{#prop-linalg-similarity-equivalence}
1. 相似は $M_n(K)$ 上の同値関係である。
2. $T\in\End(V)$ と $V$ の基底 $\mathcal B,\mathcal B'$ について、$P=P_{\mathcal B\leftarrow\mathcal B'}$ とおけば $[T]_{\mathcal B'}=P^{-1}[T]_{\mathcal B}P$ である。
3. 逆に $B=P^{-1}[T]_{\mathcal B}P$（$P$ は正則）ならば、$[T]_{\mathcal B'}=B$ となる基底 $\mathcal B'$ が存在する。

したがって、2 つの正方行列が相似であることと、それらが同じ線形変換を（一般に異なる）基底で表したものであることは同値である。
:::

:::proof
1. $A=I^{-1}AI$。$B=P^{-1}AP$ ならば $A=PBP^{-1}=(P^{-1})^{-1}BP^{-1}$。$B=P^{-1}AP$, $C=Q^{-1}BQ$ ならば $C=(PQ)^{-1}A(PQ)$。
2. :ref[thm-linalg-change-of-basis] で $V=W$, $\mathcal C=\mathcal B$, $\mathcal C'=\mathcal B'$ とすればよい。
3. :ref[prop-linalg-change-of-basis-matrix] の 4 により $P_{\mathcal B\leftarrow\mathcal B'}=P$ となる基底 $\mathcal B'$ があり、2 から $[T]_{\mathcal B'}=P^{-1}[T]_{\mathcal B}P=B$ である。
:::

行列の同値類は階数だけで決まった（:ref[cor-linalg-equivalence-rank]）が、相似による分類ははるかに精密で難しい。相似で変わらない量を調べるために、トレースの基本性質を示す。

:::lemma[トレースの巡回性]{#lem-linalg-trace-commute}
$A\in M_{m,n}(K)$, $B\in M_{n,m}(K)$ に対し $\tr(AB)=\tr(BA)$ である。
:::

:::proof
$AB\in M_m(K)$, $BA\in M_n(K)$ であり、

$$
\tr(AB)=\sum_{i=1}^m\sum_{j=1}^n a_{ij}b_{ji}=\sum_{j=1}^n\sum_{i=1}^m b_{ji}a_{ij}=\tr(BA).
$$
:::

:::proposition[相似で不変な量]{#prop-linalg-similarity-invariants}
$B=P^{-1}AP$ とすると、次が成り立つ。

1. $\rank B=\rank A$。
2. $\tr B=\tr A$。
3. 任意の $k\ge0$ について $B^k=P^{-1}A^kP$。
4. $B$ が正則であることと $A$ が正則であることは同値である。
:::

:::proof
1 は :ref[prop-linalg-rank-product] の 2 による。2 は :ref[lem-linalg-trace-commute] より $\tr(P^{-1}(AP))=\tr((AP)P^{-1})=\tr A$。3 は $k$ に関する帰納法で、$B^{k+1}=(P^{-1}A^kP)(P^{-1}AP)=P^{-1}A^{k+1}P$。4 は 1 と :ref[thm-linalg-invertible-equiv] による。
:::

:::definition[線形変換のトレースと階数]{#def-linalg-operator-trace}
有限次元空間 $V$ 上の線形変換 $T$ に対し、$V$ の基底 $\mathcal B$ を 1 つ選んで

$$
\tr T\defeq\tr[T]_{\mathcal B}
$$

と定め、$T$ の**トレース**という。:ref[prop-linalg-similarity-equivalence] の 2 と :ref[prop-linalg-similarity-invariants] の 2 により、これは基底の選び方によらない。同様に、$\rank[T]_{\mathcal B}$ も基底によらず $\rank T$ に等しい（:ref[prop-linalg-rank-map-matrix]）。
:::

行列式や特性多項式も相似で不変であり、基底によらない線形変換の量を与える（[行列式](/math/linear-algebra/determinants/)、[固有値と対角化](/math/linear-algebra/eigenvalues/)）。

:::example[成分の交換]{#exm-linalg-swap-diagonal}
$T(\transp{(x,y)})=\transp{(y,x)}$ は、標準基底では $A=\begin{pmatrix}0&1\\1&0\end{pmatrix}$ で表される。新しい基底を $u=\transp{(1,1)}$, $v=\transp{(1,-1)}$ の順に選ぶと

$$
P=\begin{pmatrix}1&1\\1&-1\end{pmatrix},\qquad
P^{-1}=\frac12\begin{pmatrix}1&1\\1&-1\end{pmatrix},\qquad
P^{-1}AP=\begin{pmatrix}1&0\\0&-1\end{pmatrix}
$$

となる。実際 $T(u)=u$, $T(v)=-v$ なので、この基底では第 2 座標の符号を変えるだけの変換である。トレースは $0$、階数は $2$ で、どちらの表示でも一致している。このように線形変換をなるべく簡単な行列で表す基底の選び方を、[固有値と対角化](/math/linear-algebra/eigenvalues/) で調べる。
:::

:::remark
階数とトレースが等しくても相似とは限らない。例えば $I_2$ と $N=\begin{pmatrix}1&1\\0&1\end{pmatrix}$ はどちらも階数 $2$、トレース $2$ であるが、任意の正則行列 $P$ について $P^{-1}I_2P=I_2\ne N$ なので相似でない。相似による完全な分類は [最小多項式と Jordan 標準形](/math/linear-algebra/jordan-form/) で与えられる。
:::

## 演習問題

::::exercise{level=basic}
次の写像が線形写像かどうかを判定せよ。$A\in M_n(K)$ は固定した行列とする。

1. $f:\R^2\to\R^2$, $f(\transp{(x,y)})=\transp{(x+y,\ xy)}$
2. $g:M_n(K)\to M_n(K)$, $g(X)=AX-XA$
3. $h:K[t]\to K[t]$, $h(p)(t)=p(t^2)$
4. $k:\R^2\to\R$, $k(\transp{(x,y)})=\abs{x}+y$
:::solution
1. 線形でない。$f(\transp{(1,1)})=\transp{(2,1)}$ だが $f(\transp{(2,2)})=\transp{(4,4)}\ne2f(\transp{(1,1)})=\transp{(4,2)}$。
2. 線形である。$g(aX+bY)=A(aX+bY)-(aX+bY)A=a(AX-XA)+b(AY-YA)$（:ref[prop-linalg-product-laws]）。
3. 線形である。$p=\sum_kc_kt^k$ ならば $h(p)=\sum_kc_kt^{2k}$ であり、$h(p)$ の各係数は $p$ の係数そのもの（または $0$）なので、係数ごとに計算される和とスカラー倍を保つ。
4. 線形でない。$k(\transp{(-1,0)})=1$ だが $-k(\transp{(1,0)})=-1$ であり、斉次性 $k(-v)=-k(v)$ が成り立たない。
:::
::::

::::exercise{level=basic}
$T:K^2\to K^2$ は線形で、$\transp{(1,1)}$ を $\transp{(3,1)}$ に、$\transp{(1,0)}$ を $\transp{(1,-1)}$ に移す。$T$ の標準基底に関する表現行列と $T(\transp{(2,3)})$ を求めよ。
:::solution
$e_1=\transp{(1,0)}$ の像は $\transp{(1,-1)}$ である。$e_2=\transp{(1,1)}-e_1$ なので、線形性より $T(e_2)=\transp{(3,1)}-\transp{(1,-1)}=\transp{(2,2)}$ である。:ref[prop-linalg-matrix-map] より

$$
A=\begin{pmatrix}1&2\\-1&2\end{pmatrix},\qquad
T\begin{pmatrix}2\\3\end{pmatrix}=A\begin{pmatrix}2\\3\end{pmatrix}=\begin{pmatrix}8\\4\end{pmatrix}.
$$

$\transp{(1,1)},\transp{(1,0)}$ は $K^2$ の基底なので、:ref[prop-linalg-map-on-basis] によりこのような $T$ はただ 1 つであり、答えも一意である。
:::
::::

::::exercise{level=basic}
$A=\begin{pmatrix}1&2&3\\2&4&6\\0&1&1\end{pmatrix}$ の列空間、行空間、$\Ker L_A$ の基底をそれぞれ 1 組求めよ。
:::solution
第 2 行から第 1 行の 2 倍を引き、第 2 行と第 3 行を入れ替え、さらに第 1 行から新しい第 2 行の 2 倍を引くと

$$
R=\begin{pmatrix}1&0&1\\0&1&1\\0&0&0\end{pmatrix}
$$

となる。主列は第 1, 2 列なので、列空間の基底には元の行列の第 1, 2 列 $\transp{(1,2,0)},\transp{(2,4,1)}$ をとる（$R$ の列 $e_1,e_2$ ではない）。行空間の基底には $R$ の零でない行の転置 $\transp{(1,0,1)},\transp{(0,1,1)}$ をとる。$Rx=0$ は $x_1=-x_3$, $x_2=-x_3$ と同値なので、核の基底は $\transp{(-1,-1,1)}$ である。:ref[cor-linalg-solution-space-dim] のとおり $\dim\Ker L_A=3-2=1$ である。
:::
::::

::::exercise{level=basic}
$D:K[t]_{\le2}\to K[t]_{\le1}$ を微分とする。入力側に基底 $(1,\ t-1,\ (t-1)^2)$、出力側に基底 $(1,t)$ を用いるときの表現行列を求めよ。さらに出力側の基底を $(1,\ t-1)$ に変えるとどうなるか。
:::solution
入力基底の像は $D(1)=0$, $D(t-1)=1$, $D((t-1)^2)=2(t-1)=2t-2$ である。$2t-2$ の座標は、基底 $(1,t)$ では $\transp{(-2,2)}$、基底 $(1,t-1)$ では $\transp{(0,2)}$ である。したがって求める表現行列はそれぞれ

$$
\begin{pmatrix}0&1&-2\\0&0&2\end{pmatrix},\qquad
\begin{pmatrix}0&1&0\\0&0&2\end{pmatrix}
$$

である。第 3 列はどちらも同じ多項式 $2(t-1)$ を表している。後者は :ref[exm-linalg-derivative-matrix] の行列と同じ形であり、「$t-1$ を変数とみて微分する」ことを表している。
:::
::::

::::exercise{level=basic}
$c_0=0$, $c_1=1$, $c_2=2$ について、:ref[exm-linalg-lagrange-interpolation] の多項式 $\ell_0,\ell_1,\ell_2$ を求めよ。また $p(0)=1$, $p(1)=3$, $p(2)=7$ を満たす 2 次以下の多項式 $p$ と、基底 $(\ell_0,\ell_1,\ell_2)$ に関する $t^2$ の座標を求めよ。
:::solution
定義から

$$
\ell_0=\frac{(t-1)(t-2)}{(0-1)(0-2)}=\frac{t^2-3t+2}{2},\quad
\ell_1=\frac{t(t-2)}{1\cdot(1-2)}=-t^2+2t,\quad
\ell_2=\frac{t(t-1)}{2\cdot1}=\frac{t^2-t}{2}.
$$

$p=1\cdot\ell_0+3\ell_1+7\ell_2=\tfrac12(t^2-3t+2)+(-3t^2+6t)+\tfrac72(t^2-t)=t^2+t+1$ であり、実際 $p(0)=1$, $p(1)=3$, $p(2)=7$ である。このような 2 次以下の多項式はこれ 1 つに限る。$t^2$ の座標は $\transp{(0^2,1^2,2^2)}=\transp{(0,1,4)}$ であり、$\ell_1+4\ell_2=-t^2+2t+2t^2-2t=t^2$ と確かめられる。
:::
::::

::::exercise{level=basic}
$A=\begin{pmatrix}3&-1\\2&0\end{pmatrix}$ の定める線形変換 $L_A$ を、基底 $\mathcal B'=(\transp{(1,1)},\transp{(1,2)})$ で表せ。それを用いて $A^k$（$k\ge0$）を求めよ。
:::solution
$P=P_{\mathcal E\leftarrow\mathcal B'}=\begin{pmatrix}1&1\\1&2\end{pmatrix}$, $P^{-1}=\begin{pmatrix}2&-1\\-1&1\end{pmatrix}$ である（:ref[exm-linalg-coordinate-change]）。$AP=\begin{pmatrix}2&1\\2&2\end{pmatrix}$ より

$$
[L_A]_{\mathcal B'}=P^{-1}AP=\begin{pmatrix}2&0\\0&1\end{pmatrix}.
$$

実際 $A\transp{(1,1)}=2\transp{(1,1)}$, $A\transp{(1,2)}=\transp{(1,2)}$ である。:ref[prop-linalg-similarity-invariants] の 3 より $A^k=P\diag(2^k,1)P^{-1}$ であり、

$$
A^k=\begin{pmatrix}2^k&1\\2^k&2\end{pmatrix}\begin{pmatrix}2&-1\\-1&1\end{pmatrix}
=\begin{pmatrix}2^{k+1}-1&1-2^k\\2^{k+1}-2&2-2^k\end{pmatrix}.
$$

$k=0$ で $I_2$、$k=1$ で $A$ に一致する。
:::
::::

::::exercise{level=basic}
トレース $\tr:M_n(K)\to K$ について、$\dim\Ker\tr=n^2-1$ を示し、$\Ker\tr$ の基底を 1 組求めよ。また任意の $X,Y\in M_n(K)$ について $XY-YX\in\Ker\tr$ を示せ。
:::solution
$\tr$ は線形であり、$\tr(cE_{11})=c$ なので全射である。次元定理より $\dim\Ker\tr=n^2-1$ である。

$E_{ij}$（$i\ne j$、$n^2-n$ 個）と $E_{ii}-E_{nn}$（$1\le i\le n-1$、$n-1$ 個）はいずれもトレースが $0$ である。これら計 $n^2-1$ 個が一次独立であることを示す。$\sum_{i\ne j}a_{ij}E_{ij}+\sum_{i<n}b_i(E_{ii}-E_{nn})=O$ とすると、$(i,j)$ 成分（$i\ne j$）から $a_{ij}=0$、$(i,i)$ 成分（$i<n$）から $b_i=0$ である。よって :ref[cor-linalg-n-vectors-basis] により、これらは $\Ker\tr$ の基底である。最後の主張は :ref[lem-linalg-trace-commute] より $\tr(XY-YX)=\tr(XY)-\tr(YX)=0$。
:::
::::

::::exercise{level=advanced}
有限次元空間の間の線形写像 $f:U\to V$, $g:V\to W$ について、

$$
\rank(g\circ f)=\rank f-\dim(\Img f\cap\Ker g)
$$

を示せ。これを用いて、$A\in M_{l,m}(K)$, $B\in M_{m,n}(K)$ について **Sylvester の不等式** $\rank A+\rank B-m\le\rank(AB)$ を証明せよ。
:::hint
$g$ を $\Img f$ に制限した写像に次元定理を適用する。
:::
:::solution
$h\defeq g|_{\Img f}:\Img f\to W$ とおく。$h(y)=0$ となるのは $y\in\Img f\cap\Ker g$ のときなので $\Ker h=\Img f\cap\Ker g$ である。また $\Img h=g(\Img f)=\Img(g\circ f)$ である。$h$ に次元定理を適用すると

$$
\rank f=\dim\Img f=\dim\Ker h+\dim\Img h=\dim(\Img f\cap\Ker g)+\rank(g\circ f)
$$

となる。$f=L_B:K^n\to K^m$, $g=L_A:K^m\to K^l$ とすると、$\dim(\Img f\cap\Ker g)\le\dim\Ker L_A=m-\rank A$（:ref[cor-linalg-solution-space-dim]）なので

$$
\rank(AB)=\rank B-\dim(\Img L_B\cap\Ker L_A)\ge\rank B-(m-\rank A).
$$

等号は $\Ker L_A\subset\Img L_B$ のときに成り立つ。
:::
::::

::::exercise{level=advanced}
有限次元の $V$ 上の線形変換 $H$ が $H^2=H$ を満たすとする。

1. $V=\Img H\oplus\Ker H$ を示し、$H$ はこの直和分解に関する $\Img H$ への射影（:ref[exm-linalg-projection]）であることを示せ。
2. この分解に沿った基底での表現行列を求め、$\tr H=(\rank H)\cdot1$ を示せ。
:::solution
1. 各 $v$ を $v=H(v)+(v-H(v))$ と書く。第 1 項は $\Img H$ に属し、第 2 項は $H(v-H(v))=H(v)-H^2(v)=0$ より $\Ker H$ に属する。よって $V=\Img H+\Ker H$。また $y=H(x)\in\Ker H$ ならば $y=H(x)=H^2(x)=H(y)=0$ である。:ref[prop-linalg-direct-sum-criterion] より和は直和である。$y\in\Img H$ ならば $H(y)=y$（上の計算）、$u\in\Ker H$ ならば $H(u)=0$ なので、$H(y+u)=y$ となり、$H$ は $\Ker H$ に沿った $\Img H$ への射影である。

2. $\Img H$ の基底の後に $\Ker H$ の基底を並べれば $V$ の基底 $\mathcal B$ を得る。上で見たとおり

$$
[H]_{\mathcal B}=\begin{pmatrix}I_r&O\\O&O\end{pmatrix},\qquad r=\rank H
$$

である（$r=0$ または $r=n$ なら対応するブロックを省く）。よって $\tr H=\tr[H]_{\mathcal B}=r\cdot1$（$1$ を $r$ 回加えたもの）である。
:::
::::

::::exercise{level=advanced}
$K$ を任意の体とし、$A\in M_{m,n}(K)$ とする。

1. $BA=I_n$ を満たす $B\in M_{n,m}(K)$ が存在する $\iff$ $\rank A=n$ を示せ。
2. $AC=I_m$ を満たす $C\in M_{n,m}(K)$ が存在する $\iff$ $\rank A=m$ を示せ。
3. $m\ne n$ ならば、$BA=I_n$ と $AC=I_m$ を同時に満たす $B,C$ は存在しないことを示せ。
:::solution
1. $BA=I_n$ ならば $L_B\circ L_A=\id$ なので $L_A$ は単射であり、:ref[cor-linalg-solution-space-dim] から $n-\rank A=0$。逆に $\rank A=n$ ならば $\Ker L_A=\{0\}$ で $L_A$ は単射なので、:ref[prop-linalg-injective-basis] より $Ae_1,\dots,Ae_n$ は一次独立である。これを $K^m$ の基底 $(Ae_1,\dots,Ae_n,w_{n+1},\dots,w_m)$ に延長し（:ref[thm-linalg-basis-extension]）、:ref[prop-linalg-map-on-basis] により $g(Ae_j)=e_j$, $g(w_k)=0$ となる線形写像 $g:K^m\to K^n$ をとる。$g=L_B$（:ref[prop-linalg-matrix-map]）とすれば $BAe_j=e_j$ がすべての $j$ で成り立つので $BA=I_n$ である。

2. $AC=I_m$ ならば $L_A\circ L_C=\id$ なので $L_A$ は全射であり、$\rank A=\dim\Img L_A=m$。逆に $\rank A=m$ ならば $L_A$ は全射なので、各 $i$ について $Ac_i=e_i$ となる $c_i\in K^n$ がある。$C=(c_1\ \cdots\ c_m)$ とおけば $AC=I_m$。

3. 両方が存在すれば 1, 2 より $n=\rank A=m$ である。これで [連立一次方程式と行列](/math/linear-algebra/linear-equations/) の最後の演習の 3 が、トレースを使わずに任意の体で示された。なお、片側逆行列は一意とは限らない。例えば $A=\transp{(1,0)}\in M_{2,1}(K)$ に対し、任意の $b\in K$ について $B=(1\ \ b)$ が $BA=I_1$ を満たす。
:::
::::

::::exercise{level=advanced}
$\dim V=n$ とし、$W$ を $V$ の $k$ 次元部分空間とする。$\mathcal S_W\defeq\{T\in\End(V)\mid T(W)\subset W\}$ は $\End(V)$ の部分空間であり、$\dim\mathcal S_W=n^2-k(n-k)$ であることを示せ。
:::solution
$W$ の基底 $(w_1,\dots,w_k)$ を $V$ の基底 $\mathcal B=(w_1,\dots,w_k,v_{k+1},\dots,v_n)$ に延長する。$T(W)\subset W$ は $T(w_j)\in W$（$1\le j\le k$）と同値である（$W=\langle w_1,\dots,w_k\rangle$ と :ref[prop-linalg-image-preimage]）。$T(w_j)\in W$ は、$[T(w_j)]_{\mathcal B}$ の第 $k+1$ 成分以下が $0$ であること、すなわち $[T]_{\mathcal B}$ の第 $j$ 列の第 $k+1$ 行以下が $0$ であることと同値である。よって

$$
T\in\mathcal S_W\iff[T]_{\mathcal B}=\begin{pmatrix}X&Y\\O&Z\end{pmatrix}\quad(X\in M_k(K),\ Z\in M_{n-k}(K))
$$

である。このような行列全体は、左下の $(n-k)\times k$ ブロックの成分がすべて $0$ という条件で定まる $M_n(K)$ の部分空間であり、行列単位のうち左下ブロック以外の $n^2-k(n-k)$ 個が基底をなす。:ref[prop-linalg-hom-matrix-iso] の同型 $T\mapsto[T]_{\mathcal B}$ により、$\mathcal S_W$ はこの部分空間に対応するので、主張が従う。
:::
::::

::::exercise{level=advanced}
$V$ を有限次元とし、$T\in\End(V)$ がすべての $S\in\End(V)$ と可換（$TS=ST$）であるとする。$T=c\,\id_V$ となる $c\in K$ が存在することを示せ。
:::hint
$v\ne0$ に対し、$\langle v\rangle$ への射影 $S$ を考える。
:::
:::solution
$v\ne0$ を任意にとる。:ref[cor-linalg-complement-exists] により $V=\langle v\rangle\oplus U$ となる $U$ があり、この分解に関する $\langle v\rangle$ への射影を $S$ とする（:ref[exm-linalg-projection]）。$S(v)=v$ なので $T(v)=T(S(v))=S(T(v))\in\Img S=\langle v\rangle$ である。よって各 $v\ne0$ に対し $T(v)=\lambda_vv$ となる $\lambda_v\in K$ がある。

$v,w\ne0$ について $\lambda_v=\lambda_w$ を示す。$v,w$ が一次独立ならば、$v+w\ne0$ であり

$$
\lambda_{v+w}v+\lambda_{v+w}w=T(v+w)=\lambda_vv+\lambda_ww
$$

から $\lambda_v=\lambda_{v+w}=\lambda_w$ である。一次従属ならば $w=av$（$a\ne0$）で、$T(w)=aT(v)=\lambda_vw$ より $\lambda_w=\lambda_v$ である。この共通の値を $c$ とすれば $T=c\,\id_V$ である（$V=\{0\}$ のときは任意の $c$ でよい）。行列でいえば、すべての $n$ 次正方行列と可換な行列はスカラー行列 $cI_n$ に限る。
:::
::::

## 確認テスト

::::quiz
線形写像 $f:V\to W$（$V,W$ は有限次元）が単射であることと同値な条件を選べ。

- [ ] $\Img f=W$
- [x] $\Ker f=\{0\}$
- [ ] $f(0)=0$
- [ ] $\dim V\le\dim W$

:::explanation
$f(u)=f(v)$ は $u-v\in\Ker f$ と同値なので、単射性は $\Ker f=\{0\}$ と同値である（:ref[prop-linalg-kernel-injective]）。$\Img f=W$ は全射性の条件であり、$f(0)=0$ はすべての線形写像が満たす。$\dim V\le\dim W$ は単射であるための必要条件（:ref[cor-linalg-map-dimension-bounds]）だが十分条件ではない（零写像が反例）。
:::
::::

::::quiz{answer="3"}
$\dim V=7$ で、線形写像 $f:V\to K^5$ の像の次元が $4$ である。核の次元はいくつか。

:::explanation
次元定理 :eqref[eq-linalg-rank-nullity] より $\dim\Ker f=7-4=3$ である。終域の次元 $5$ は関係しない。
:::
::::

::::quiz
基底変換行列 $P_{\mathcal B\leftarrow\mathcal B'}$ の第 $j$ 列は何か。

- [ ] 古い基底 $\mathcal B$ の第 $j$ ベクトルを新しい基底 $\mathcal B'$ で表した座標
- [ ] 古い基底 $\mathcal B$ の第 $j$ ベクトルを古い基底 $\mathcal B$ で表した座標
- [x] 新しい基底 $\mathcal B'$ の第 $j$ ベクトルを古い基底 $\mathcal B$ で表した座標

:::explanation
$P_{\mathcal B\leftarrow\mathcal B'}=[\id_V]_{\mathcal B\leftarrow\mathcal B'}$ であり、表現行列の第 $j$ 列は入力基底の第 $j$ ベクトル $v'_j$ の像 $\id_V(v'_j)=v'_j$ を出力基底 $\mathcal B$ で表した座標である。この行列は $\mathcal B'$ に関する座標を $\mathcal B$ に関する座標に変換する（:ref[prop-linalg-change-of-basis-matrix]）。
:::
::::

::::quiz{answer="2"}
$4\times 5$ 行列 $A$ の行空間が $2$ 次元である。$\rank\transp{A}$ はいくつか。

:::explanation
$\transp{A}$ の列空間は $A$ の行空間であり、:ref[thm-linalg-row-column-rank] よりその次元 $2$ が $\rank\transp{A}$ に等しい（$\rank A$ も $2$ である）。
:::
::::

::::quiz{answer="12"}
$\dim_\R\Hom(\R^3,\R^4)$ はいくつか。

:::explanation
:ref[prop-linalg-hom-matrix-iso] より $\Hom(\R^3,\R^4)\cong M_{4,3}(\R)$ であり、次元は $4\cdot3=12$ である。
:::
::::

::::quiz
$A,B\in M_n(K)$ が相似であるとき、必ず一致するものをすべて選べ。

- [x] 階数
- [x] トレース
- [ ] $(1,1)$ 成分
- [x] 正則であるかどうか
- [ ] 対称行列であるかどうか

:::explanation
階数・トレース・正則性は :ref[prop-linalg-similarity-invariants] により相似で不変である。$(1,1)$ 成分は不変でない：:ref[exm-linalg-swap-diagonal] で $\begin{pmatrix}0&1\\1&0\end{pmatrix}$ は $\diag(1,-1)$ と相似である。対称性も不変でない：$A=\diag(1,2)$, $P=\begin{pmatrix}1&1\\0&1\end{pmatrix}$ とすると $P^{-1}AP=\begin{pmatrix}1&-1\\0&2\end{pmatrix}$ は対称でない。
:::
::::

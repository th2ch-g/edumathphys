---
title: ベクトル空間・基底・次元
description: 体とベクトル空間の公理、部分空間、一次独立性、基底と次元、和空間と直和
notation:
  - symbol: '\langle S\rangle'
    meaning: 集合 S で生成される部分空間（S の張る空間）
  - symbol: '\dim V'
    meaning: ベクトル空間 V の次元
  - symbol: 'W_1+W_2'
    meaning: 部分空間の和
  - symbol: 'W_1\oplus W_2'
    meaning: 部分空間の直和
  - symbol: 'K[x]_{\le n}'
    meaning: 次数 n 以下の多項式全体
---

## この章の目標

- 体とベクトル空間を公理によって定義し、数ベクトル以外の多様な例を同じ枠組みで扱えるようにする。
- 一次独立性と生成系の概念を理解し、基底の存在と「基底の個数は一定である」ことを証明する。
- 次元の性質（部分空間の次元・基底の延長・次元公式）を自在に使えるようになる。
- 和空間と直和を理解する。

前章 [連立一次方程式と行列](/math/linear-algebra/linear-equations/) の結果、特に :ref[cor-linalg-homogeneous-nontrivial] が、本章の次元論の要となる。

## 体とベクトル空間の公理

### 体

:::definition[体]{#def-linalg-field}
集合 $K$ に 2 つの演算、加法 $+$ と乗法 $\cdot$ が定義されていて、次の条件を満たすとき、$K$ を :term[体]{en="field" yomi="たい"} という。

1. 加法について：結合法則 $(a+b)+c=a+(b+c)$、交換法則 $a+b=b+a$ が成り立ち、零元 $0$（$a+0=a$）と各 $a$ の加法逆元 $-a$（$a+(-a)=0$）が存在する。
2. 乗法について：結合法則と交換法則が成り立ち、単位元 $1\ne 0$（$1\cdot a=a$）と、$0$ でない各 $a$ の乗法逆元 $a^{-1}$（$a\cdot a^{-1}=1$）が存在する。
3. 分配法則 $a\cdot(b+c)=a\cdot b+a\cdot c$ が成り立つ。
:::

:::example{#exm-linalg-fields}
有理数体 $\Q$、実数体 $\R$、複素数体 $\C$ は体である。整数全体 $\Z$ は $2$ の乗法逆元が $\Z$ 内に存在しないので体ではない。

$p$ を素数とするとき、$\{0,1,\dots,p-1\}$ に $p$ で割った余りによる加法と乗法を定めたもの $\F_p$ は体になる（乗法逆元の存在に $p$ が素数であることを用いる。[集合と論理](/math/sets-logic/) の同値関係の章や [群論](/math/group-theory/) を参照）。特に $\F_2=\{0,1\}$ では $1+1=0$ である。
:::

体の公理から、数の場合と同じ計算規則（$0\cdot a=0$、$(-1)a=-a$、$ab=0\Rightarrow a=0$ または $b=0$ など）が導かれる。前章の議論はすべて体の公理のみを用いていたので、任意の体 $K$ 上で成り立つ。

### ベクトル空間

:::definition[ベクトル空間]{#def-linalg-vector-space}
$K$ を体とする。集合 $V$ に、加法 $V\times V\to V$, $(u,v)\mapsto u+v$ と、スカラー倍 $K\times V\to V$, $(c,v)\mapsto cv$ が定義されていて、任意の $u,v,w\in V$, $a,b\in K$ について次が成り立つとき、$V$ を $K$ 上の :term[ベクトル空間]{en="vector space" yomi="べくとるくうかん"}（線形空間）といい、その元を**ベクトル**と呼ぶ。

1. $(u+v)+w=u+(v+w)$
2. $u+v=v+u$
3. ある元 $0\in V$（**零ベクトル**）が存在して、すべての $v$ について $v+0=v$
4. 各 $v$ に対し $v+v'=0$ となる $v'\in V$ が存在する（これを $-v$ と書く）
5. $a(u+v)=au+av$
6. $(a+b)v=av+bv$
7. $(ab)v=a(bv)$
8. $1v=v$

$K=\R$ のとき**実ベクトル空間**、$K=\C$ のとき**複素ベクトル空間**という。
:::

:::proposition{#prop-linalg-vs-basic}
$V$ を $K$ 上のベクトル空間とし、$v\in V$, $c\in K$ とする。

1. 零ベクトルと各 $v$ の逆元 $-v$ は一意である。
2. $0v=0$, $c0=0$, $(-1)v=-v$。
3. $cv=0$ ならば $c=0$ または $v=0$。
:::

:::proof
1. $0, 0'$ がともに零ベクトルなら $0=0+0'=0'$。$v', v''$ がともに $v$ の逆元なら $v'=v'+(v+v'')=(v'+v)+v''=v''$。

2. $0v=(0+0)v=0v+0v$ の両辺に $-(0v)$ を加えて $0=0v$。$c0=0$ も同様。$v+(-1)v=(1+(-1))v=0v=0$ なので $(-1)v=-v$。

3. $c\ne0$ ならば $v=1v=(c^{-1}c)v=c^{-1}(cv)=c^{-1}0=0$。
:::

:::example[ベクトル空間の例]{#exm-linalg-vector-spaces}
以下はいずれも、自然な加法とスカラー倍について $K$ 上のベクトル空間である。

1. 数ベクトル空間 $K^n$。
2. 行列の空間 $M_{m,n}(K)$。
3. $K$ 係数の多項式全体 $K[x]$、および次数が $n$ 以下の多項式全体 $K[x]_{\le n}$（零多項式も含める）。
4. 集合 $X$ から $K$ への写像全体 $\mathcal F(X,K)$。和とスカラー倍は $(f+g)(t)=f(t)+g(t)$, $(cf)(t)=cf(t)$ で定める。$X=\{1,\dots,n\}$ とすれば $K^n$ と同じものである。
5. 区間 $[a,b]$ 上の実数値連続関数全体 $C([a,b])$（実ベクトル空間）。
6. $A\in M_{m,n}(K)$ に対し、斉次方程式 $Ax=0$ の解全体 $\{x\in K^n\mid Ax=0\}$。
7. $\C$ は $\C$ 上のベクトル空間であると同時に、スカラー倍を実数に制限すれば $\R$ 上のベクトル空間でもある。
:::

## 部分空間

:::definition[部分空間]{#def-linalg-subspace}
$V$ を $K$ 上のベクトル空間とする。部分集合 $W\subset V$ が次を満たすとき、$W$ を $V$ の :term[部分空間]{en="subspace" yomi="ぶぶんくうかん"} という。

1. $0\in W$
2. $u,v\in W\Rightarrow u+v\in W$
3. $c\in K$, $v\in W\Rightarrow cv\in W$
:::

部分空間は $V$ の演算を制限することで、それ自身ベクトル空間になる（公理 1, 2, 5–8 は $V$ で成り立つので $W$ でも成り立ち、公理 4 は $-v=(-1)v\in W$ による）。$\{0\}$ と $V$ 自身は常に部分空間であり、これらを**自明な部分空間**という。

:::example{#exm-linalg-subspaces}
1. $A\in M_{m,n}(K)$ に対し $\{x\in K^n\mid Ax=0\}$ は $K^n$ の部分空間である。これを $A$ の**解空間**という。一方、$b\ne0$ のとき $\{x\mid Ax=b\}$ は $0$ を含まないので部分空間ではない。
2. $K[x]_{\le n}$ は $K[x]$ の部分空間である。
3. $\R^2$ において、原点を通る直線 $\{t(a,b)^\top\mid t\in\R\}$ は部分空間だが、原点を通らない直線や、2 本の直線の和集合は部分空間でない。
:::

:::proposition{#prop-linalg-subspace-intersection}
$\{W_\lambda\}_{\lambda\in\Lambda}$ を $V$ の部分空間の族とすると、共通部分 $\bigcap_\lambda W_\lambda$ も $V$ の部分空間である。
:::

:::proof
各 $W_\lambda$ が条件 1–3 を満たすので、共通部分も満たす。
:::

## 一次独立と生成系

### 一次結合と生成される部分空間

:::definition[一次結合・生成系]{#def-linalg-span}
$v_1,\dots,v_k\in V$ と $c_1,\dots,c_k\in K$ に対し、$c_1v_1+\cdots+c_kv_k$ の形のベクトルを $v_1,\dots,v_k$ の :term[一次結合]{en="linear combination" yomi="いちじけつごう"} という。部分集合 $S\subset V$ の元の一次結合全体を $\langle S\rangle$（あるいは $\spn S$）と書き、$S$ で**生成される部分空間**（$S$ の張る空間）という。ただし $\langle\varnothing\rangle=\{0\}$ とする。$\langle S\rangle=V$ のとき、$S$ は $V$ を**生成する**、あるいは $S$ は $V$ の :term[生成系]{en="spanning set" yomi="せいせいけい"} であるという。
:::

:::proposition{#prop-linalg-span-smallest}
$\langle S\rangle$ は $S$ を含む $V$ の部分空間のうちで最小のものである。すなわち、$\langle S\rangle$ は $S$ を含む部分空間すべての共通部分に等しい。
:::

:::proof
一次結合の和やスカラー倍はまた一次結合なので、$\langle S\rangle$ は部分空間であり、明らかに $S$ を含む。逆に $S$ を含む部分空間 $W$ は、条件 2, 3 によって $S$ の元の一次結合をすべて含むので $\langle S\rangle\subset W$ である。
:::

### 一次独立性

:::definition[一次独立・一次従属]{#def-linalg-linear-independence}
$v_1,\dots,v_k\in V$ が :term[一次独立]{en="linearly independent" yomi="いちじどくりつ"} であるとは、

$$
c_1v_1+\cdots+c_kv_k=0 \implies c_1=\cdots=c_k=0
$$

が成り立つことをいう。一次独立でないとき :term[一次従属]{en="linearly dependent" yomi="いちじじゅうぞく"} という。無限集合 $S\subset V$ は、その任意の有限個の相異なる元が一次独立であるとき一次独立であるという。
:::

一次従属とは、少なくとも 1 つは $0$ でない係数をもつ自明でない一次関係が存在することである。$v_1,\dots,v_k$ が一次従属であることと、ある $v_i$ が他のベクトルの一次結合で書けることは同値である（$c_i\ne0$ ならば $v_i=-c_i^{-1}\sum_{j\ne i}c_jv_j$）。

:::proposition[一次結合の表示の一意性]{#prop-linalg-independence-unique}
$v_1,\dots,v_k$ が一次独立であることは、$\langle v_1,\dots,v_k\rangle$ の各元が $v_1,\dots,v_k$ の一次結合として**ただ一通りに**表されることと同値である。
:::

:::proof
$\sum c_iv_i=\sum d_iv_i$ ならば $\sum(c_i-d_i)v_i=0$ である。一次独立なら $c_i=d_i$ となる。逆に表示が一意なら、$\sum c_iv_i=0=\sum 0\,v_i$ から $c_i=0$ である。
:::

:::lemma{#lem-linalg-add-independent}
$v_1,\dots,v_k$ が一次独立で、$v\notin\langle v_1,\dots,v_k\rangle$ ならば、$v_1,\dots,v_k,v$ も一次独立である。
:::

:::proof
$c_1v_1+\cdots+c_kv_k+cv=0$ とする。$c\ne0$ ならば $v=-c^{-1}\sum c_iv_i\in\langle v_1,\dots,v_k\rangle$ となり仮定に反する。よって $c=0$ であり、$v_1,\dots,v_k$ の一次独立性から $c_1=\cdots=c_k=0$ である。
:::

数ベクトルの一次独立性は、階数によって判定できる。

:::proposition{#prop-linalg-independence-rank}
$a_1,\dots,a_n\in K^m$ を並べた行列を $A=(a_1\ \cdots\ a_n)\in M_{m,n}(K)$ とする。$a_1,\dots,a_n$ が一次独立であるための必要十分条件は $\rank A=n$ である。
:::

:::proof
$c=\transp{(c_1,\dots,c_n)}$ とおけば $c_1a_1+\cdots+c_na_n=Ac$ である。よって一次独立性は「$Ac=0\Rightarrow c=0$」と同値であり、:ref[cor-linalg-homogeneous-nontrivial] より $\rank A=n$ と同値である。
:::

次の補題が次元論の基礎である。

:::lemma[Steinitz の交換補題の弱形]{#lem-linalg-steinitz}
$V$ が $m$ 個のベクトル $u_1,\dots,u_m$ で生成されるならば、$V$ の $m$ 個より多いベクトルは一次従属である。
:::

:::proof
$v_1,\dots,v_n\in V$, $n>m$ とする。各 $v_j$ は $u_i$ の一次結合なので $v_j=\sum_{i=1}^m a_{ij}u_i$ と書ける。$A=(a_{ij})\in M_{m,n}(K)$ は方程式の本数 $m$ が未知数の個数 $n$ より少ないので、:ref[cor-linalg-homogeneous-nontrivial] より $Ac=0$ を満たす $c\ne0$ が存在する。このとき

$$
\sum_{j=1}^n c_jv_j = \sum_{j=1}^n c_j\sum_{i=1}^m a_{ij}u_i = \sum_{i=1}^m\Bigl(\sum_{j=1}^n a_{ij}c_j\Bigr)u_i = \sum_{i=1}^m (Ac)_i\,u_i = 0
$$

となり、$v_1,\dots,v_n$ は一次従属である。
:::

## 基底と次元

:::definition[基底]{#def-linalg-basis}
$V$ のベクトルの組 $(v_1,\dots,v_n)$ が一次独立であり、かつ $V$ を生成するとき、これを $V$ の :term[基底]{en="basis" yomi="きてい"} という。有限個のベクトルで生成されるベクトル空間を**有限生成**であるという。
:::

基底の定義では順序を考えた組を用いる（座標を考えるときに順序が必要になるため）。$K^n$ の基本ベクトル $(e_1,\dots,e_n)$ は基底であり、これを**標準基底**という。

:::theorem[基底の存在]{#thm-linalg-basis-exists}
$V\ne\{0\}$ が有限生成ならば、$V$ の任意の有限生成系から基底を選び出すことができる。特に $V$ は基底をもつ。
:::

:::proof
$S=\{u_1,\dots,u_m\}$ を生成系とする。$S$ の部分集合で一次独立なものの中から、元の個数が最大のもの $\{v_1,\dots,v_n\}$ を 1 つ選ぶ（$V\ne\{0\}$ なので $S$ は零でない元を含み、$n\ge1$ である）。任意の $u\in S$ は $\langle v_1,\dots,v_n\rangle$ に属する。実際、属さなければ :ref[lem-linalg-add-independent] により $\{v_1,\dots,v_n,u\}$ が一次独立となって、$n$ の最大性に反する。よって $V=\langle S\rangle\subset\langle v_1,\dots,v_n\rangle$ となり、$(v_1,\dots,v_n)$ は基底である。
:::

:::theorem[次元の不変性]{#thm-linalg-dimension-invariance}
有限生成なベクトル空間 $V$ の基底を構成するベクトルの個数は、基底の取り方によらず一定である。
:::

:::proof
$(v_1,\dots,v_n)$, $(w_1,\dots,w_m)$ をともに基底とする。$V$ は $v_1,\dots,v_n$ で生成され、$w_1,\dots,w_m$ は一次独立なので、:ref[lem-linalg-steinitz] より $m\le n$ である。役割を入れ替えれば $n\le m$ である。
:::

:::definition[次元]{#def-linalg-dimension}
有限生成なベクトル空間 $V$ の基底の個数を $V$ の :term[次元]{en="dimension" yomi="じげん"} といい、$\dim V$（体を明示するときは $\dim_K V$）と書く。$\dim\{0\}=0$ と定める。有限生成でないベクトル空間は**無限次元**であるといい、$\dim V=\infty$ と書く。
:::

:::example{#exm-linalg-dimensions}
1. $\dim K^n=n$（標準基底）。
2. $\dim M_{m,n}(K)=mn$。$(i,j)$ 成分のみが $1$ で他が $0$ の行列 $E_{ij}$ 全体が基底をなす（**行列単位**）。
3. $\dim K[x]_{\le n}=n+1$。$(1,x,\dots,x^n)$ が基底である。
4. $K[x]$ は無限次元である。実際、有限個の多項式の次数の最大値を $d$ とすると、それらの一次結合の次数は $d$ 以下なので $x^{d+1}$ を表せない。
5. $\dim_\C\C=1$ だが $\dim_\R\C=2$（基底 $(1,i)$）。次元はスカラーの体に依存する。
:::

:::remark
無限次元のベクトル空間も、選択公理（Zorn の補題）を用いれば基底をもつことが証明できる（[集合と論理](/math/sets-logic/) の選択公理の章で扱う）。ただし、例えば $C([a,b])$ の基底を具体的に書き下すことはできない。無限次元空間では、本章の意味の基底よりも、関数解析で扱う位相的な「完全正規直交系」のほうが有用であることが多い。
:::

### 次元の基本性質

:::theorem[基底の延長]{#thm-linalg-basis-extension}
$\dim V=n$ とし、$v_1,\dots,v_k\in V$ が一次独立であるとする。このとき $k\le n$ であり、$v_{k+1},\dots,v_n\in V$ を補って $(v_1,\dots,v_n)$ を $V$ の基底にできる。
:::

:::proof
$k\le n$ は :ref[lem-linalg-steinitz] による。$\langle v_1,\dots,v_k\rangle\ne V$ ならば、$v_{k+1}\in V\setminus\langle v_1,\dots,v_k\rangle$ をとると :ref[lem-linalg-add-independent] より $v_1,\dots,v_{k+1}$ は一次独立である。これを繰り返すと、一次独立なベクトルの個数は $n$ を超えられないので、有限回で $\langle v_1,\dots,v_l\rangle=V$ となり基底が得られる。このとき基底の個数は $n$ なので $l=n$ である。
:::

:::corollary{#cor-linalg-n-vectors-basis}
$\dim V=n$ のとき、$V$ の $n$ 個のベクトル $v_1,\dots,v_n$ について次は同値である。

1. $(v_1,\dots,v_n)$ は基底である。
2. $v_1,\dots,v_n$ は一次独立である。
3. $v_1,\dots,v_n$ は $V$ を生成する。
:::

:::proof
(2)$\Rightarrow$(1)：:ref[thm-linalg-basis-extension] で補うべきベクトルは $n-n=0$ 個である。(3)$\Rightarrow$(1)：:ref[thm-linalg-basis-exists] により $\{v_1,\dots,v_n\}$ から基底を選び出せるが、基底の個数は $n$ なので全部が選ばれる。(1)$\Rightarrow$(2), (3) は定義から明らかである。
:::

:::corollary{#cor-linalg-subspace-dimension}
$\dim V=n$ とし、$W$ を $V$ の部分空間とする。このとき $W$ は有限生成で $\dim W\le n$ であり、等号が成り立つのは $W=V$ のときに限る。
:::

:::proof
$W$ の一次独立なベクトルは高々 $n$ 個なので、:ref[lem-linalg-add-independent] を用いて $W$ の中で一次独立なベクトルを付け加えていく操作は高々 $n$ 回で止まり、そのとき $W$ の基底が得られる。よって $\dim W\le n$ である。$\dim W=n$ ならば $W$ の基底は $V$ の $n$ 個の一次独立なベクトルなので、:ref[cor-linalg-n-vectors-basis] により $V$ の基底でもあり、$W=V$ となる。
:::

### 座標

$\mathcal B=(v_1,\dots,v_n)$ を $V$ の基底とすると、:ref[prop-linalg-independence-unique] により各 $v\in V$ は $v=x_1v_1+\cdots+x_nv_n$ とただ一通りに表される。数ベクトル $\transp{(x_1,\dots,x_n)}\in K^n$ を、基底 $\mathcal B$ に関する $v$ の :term[座標]{en="coordinates" yomi="ざひょう"} という。これにより、$n$ 次元ベクトル空間の計算は $K^n$ の計算に帰着される。この対応の意味は次章で「同型写像」として明確になる。

:::example{#exm-linalg-coordinates}
$\R^3$ のベクトル $u_1=\transp{(1,1,0)}$, $u_2=\transp{(1,0,1)}$, $u_3=\transp{(0,1,1)}$ を並べた行列 $U=(u_1\ u_2\ u_3)$ は簡約化すると $I_3$ になるので、:ref[prop-linalg-independence-rank] と :ref[cor-linalg-n-vectors-basis] から $(u_1,u_2,u_3)$ は $\R^3$ の基底である。$v=\transp{(1,2,3)}$ の座標 $x$ は $Ux=v$ を解いて $x=\transp{(0,1,2)}$、すなわち $v=u_2+2u_3$ である。

一方、同じ 3 つのベクトルを $\F_2^3$ の元とみると、$u_1+u_2+u_3=\transp{(2,2,2)}=0$ となり一次従属である。一次独立性は係数体に依存する。
:::

## 和空間と直和

:::definition[和空間]{#def-linalg-sum-space}
$V$ の部分空間 $W_1, W_2$ に対し、

$$
W_1+W_2 \defeq \{\,w_1+w_2 \mid w_1\in W_1,\ w_2\in W_2\,\}
$$

を $W_1$ と $W_2$ の**和空間**という。これは $W_1\cup W_2$ で生成される部分空間 $\langle W_1\cup W_2\rangle$ に等しい。
:::

:::theorem[次元公式]{#thm-linalg-dimension-formula}
$W_1, W_2$ を有限次元ベクトル空間 $V$ の部分空間とすると

$$
\dim(W_1+W_2) = \dim W_1 + \dim W_2 - \dim(W_1\cap W_2). \label{eq-linalg-dim-formula}
$$
:::

:::proof
$\dim(W_1\cap W_2)=r$, $\dim W_1=r+s$, $\dim W_2=r+t$ とする。$W_1\cap W_2$ の基底 $(u_1,\dots,u_r)$ を :ref[thm-linalg-basis-extension] により $W_1$ の基底 $(u_1,\dots,u_r,v_1,\dots,v_s)$ と $W_2$ の基底 $(u_1,\dots,u_r,w_1,\dots,w_t)$ に延長する。$u_1,\dots,u_r,v_1,\dots,v_s,w_1,\dots,w_t$ が $W_1+W_2$ の基底であることを示せば、$\dim(W_1+W_2)=r+s+t$ となって主張が従う。

これらが $W_1+W_2$ を生成することは明らかである。一次独立性を示す。

$$
\sum_i a_iu_i + \sum_j b_jv_j + \sum_k c_kw_k = 0
$$

とすると、$\sum_k c_kw_k = -\sum_i a_iu_i-\sum_j b_jv_j$ の左辺は $W_2$ に、右辺は $W_1$ に属するので、これは $W_1\cap W_2$ の元である。よって $\sum_k c_kw_k=\sum_i d_iu_i$ と書けるが、$(u_1,\dots,u_r,w_1,\dots,w_t)$ は一次独立なので $c_k=0$（かつ $d_i=0$）である。すると $\sum a_iu_i+\sum b_jv_j=0$ となり、$W_1$ の基底の一次独立性から $a_i=b_j=0$ である。
:::

:::definition[直和]{#def-linalg-direct-sum}
$V$ の部分空間 $W_1,\dots,W_k$ について、和空間 $W=W_1+\cdots+W_k$ の各元が $w=w_1+\cdots+w_k$（$w_i\in W_i$）と**ただ一通りに**表されるとき、$W$ は $W_1,\dots,W_k$ の :term[直和]{en="direct sum" yomi="ちょくわ"} であるといい、

$$
W = W_1\oplus\cdots\oplus W_k
$$

と書く。$V=W_1\oplus W_2$ のとき、$W_2$ を $W_1$ の**補空間**という。
:::

:::proposition[直和の判定]{#prop-linalg-direct-sum-criterion}
1. $W_1+W_2$ が直和であるための必要十分条件は $W_1\cap W_2=\{0\}$ である。
2. $V$ が有限次元のとき、$W_1+\cdots+W_k$ が直和であるための必要十分条件は

$$
\dim(W_1+\cdots+W_k) = \dim W_1+\cdots+\dim W_k.
$$
:::

:::proof
1. 表示の一意性は「$w_1+w_2=0$（$w_i\in W_i$）ならば $w_1=w_2=0$」と同値である（差をとればよい）。$w_1+w_2=0$ ならば $w_1=-w_2\in W_1\cap W_2$ なので、$W_1\cap W_2=\{0\}$ ならば一意性が成り立つ。逆に $0\ne w\in W_1\cap W_2$ があれば $w+(-w)=0$ が自明でない表示を与える。

2. 各 $W_i$ の基底を並べた組を $\mathcal B$ とすると、$\mathcal B$ は $W_1+\cdots+W_k$ を生成する。直和であることは、$\sum w_i=0\Rightarrow$ すべての $w_i=0$ と同値であり、これは $\mathcal B$ が一次独立であることと同値である（各 $w_i$ を $W_i$ の基底で展開すればよい）。$\mathcal B$ は $\sum\dim W_i$ 個のベクトルからなる生成系なので、それが一次独立であることは、:ref[cor-linalg-n-vectors-basis] と :ref[thm-linalg-basis-exists] により $\dim(W_1+\cdots+W_k)=\sum\dim W_i$ と同値である。
:::

:::corollary[補空間の存在]{#cor-linalg-complement-exists}
有限次元ベクトル空間 $V$ の任意の部分空間 $W$ は補空間をもつ。
:::

:::proof
$W$ の基底 $(w_1,\dots,w_k)$ を $V$ の基底 $(w_1,\dots,w_k,v_{k+1},\dots,v_n)$ に延長し、$W'\defeq\langle v_{k+1},\dots,v_n\rangle$ とおけば $V=W\oplus W'$ である。
:::

:::remark
補空間は一意ではない。例えば $\R^2$ において $W=\langle e_1\rangle$ の補空間は、$e_1$ に平行でない任意の直線 $\langle v\rangle$ である。「直交補空間」という特別な選び方は、内積を導入する第 6 章で扱う。
:::

## 演習問題

::::exercise{level=basic}
次の $\R^3$ の部分集合のうち、部分空間であるものを選び、理由を述べよ。

1. $W_1=\{x\in\R^3\mid x_1+2x_2-x_3=0\}$
2. $W_2=\{x\in\R^3\mid x_1x_2=0\}$
3. $W_3=\{x\in\R^3\mid x_1=x_2=x_3\}$
4. $W_4=\{x\in\R^3\mid x_1+x_2+x_3=1\}$
:::solution
1. 部分空間である。$1\times3$ 行列 $(1\ 2\ {-1})$ の解空間だからである（:ref[exm-linalg-subspaces]）。
2. 部分空間でない。$e_1, e_2\in W_2$ だが $e_1+e_2=\transp{(1,1,0)}\notin W_2$。
3. 部分空間である。$x_1-x_2=0$, $x_2-x_3=0$ の解空間であり、$\langle\transp{(1,1,1)}\rangle$ に等しい。
4. 部分空間でない。$0\notin W_4$。
:::
::::

::::exercise{level=basic}
$\R^4$ のベクトル

$$
v_1=\begin{pmatrix}1\\2\\0\\1\end{pmatrix},\
v_2=\begin{pmatrix}2\\1\\1\\0\end{pmatrix},\
v_3=\begin{pmatrix}0\\3\\-1\\2\end{pmatrix},\
v_4=\begin{pmatrix}1\\-1\\1\\-1\end{pmatrix}
$$

で生成される部分空間 $W$ の次元と基底を 1 組求めよ。また $v_3, v_4$ を基底の一次結合で表せ。
:::solution
$A=(v_1\ v_2\ v_3\ v_4)$ を簡約化すると

$$
\begin{pmatrix}1&0&2&-1\\0&1&-1&1\\0&0&0&0\\0&0&0&0\end{pmatrix}
$$

となる。主列は第 1, 2 列なので、:ref[prop-linalg-independence-rank] より $v_1, v_2$ は一次独立である。また簡約行階段形の第 3, 4 列の成分から $v_3=2v_1-v_2$, $v_4=-v_1+v_2$ が読み取れる（行基本変形は列の間の一次関係を保つ。これは :ref[lem-linalg-row-equiv-solutions] による）。したがって $W=\langle v_1,v_2\rangle$, $\dim W=2$ で、$(v_1,v_2)$ が基底である。
:::
::::

::::exercise{level=basic}
$\mathrm{Sym}_n(K)=\{A\in M_n(K)\mid\transp{A}=A\}$（**対称行列**全体）と $\mathrm{Alt}_n(K)=\{A\in M_n(K)\mid\transp{A}=-A\}$（**交代行列**全体）を考える。$K=\R$ とする。

1. これらは $M_n(\R)$ の部分空間であることを示し、次元を求めよ。
2. $M_n(\R)=\mathrm{Sym}_n(\R)\oplus\mathrm{Alt}_n(\R)$ を示せ。
:::solution
1. 転置は和とスカラー倍を保つ（:ref[prop-linalg-transpose-laws]）ので、どちらも部分空間である。$i\le j$ について $E_{ij}+E_{ji}$（$i=j$ のときは $E_{ii}$）が $\mathrm{Sym}_n$ の基底をなすので $\dim\mathrm{Sym}_n=n(n+1)/2$。$i<j$ について $E_{ij}-E_{ji}$ が $\mathrm{Alt}_n$ の基底をなすので $\dim\mathrm{Alt}_n=n(n-1)/2$（交代行列の対角成分は $a_{ii}=-a_{ii}$ より $0$）。

2. 任意の $A$ は $A=\tfrac12(A+\transp{A})+\tfrac12(A-\transp{A})$ と対称部分と交代部分に分解されるので $M_n=\mathrm{Sym}_n+\mathrm{Alt}_n$。また $A\in\mathrm{Sym}_n\cap\mathrm{Alt}_n$ ならば $A=\transp{A}=-A$ より $A=O$。:ref[prop-linalg-direct-sum-criterion] により直和である。次元も $n(n+1)/2+n(n-1)/2=n^2$ と整合する。
:::
::::

::::exercise{level=basic}
$\R^4$ の部分空間

$$
W_1=\{x\mid x_1+x_2+x_3+x_4=0\},\qquad W_2=\{x\mid x_1=x_2,\ x_3=x_4\}
$$

について、$\dim W_1$, $\dim W_2$, $\dim(W_1\cap W_2)$, $\dim(W_1+W_2)$ を求めよ。
:::solution
$W_1$ は階数 1 の $1\times4$ 行列の解空間なので $\dim W_1=4-1=3$。$W_2$ は階数 2 の行列の解空間で $\dim W_2=2$（基底 $\transp{(1,1,0,0)}, \transp{(0,0,1,1)}$）。$W_1\cap W_2$ の元は $x=\transp{(s,s,t,t)}$ で $2s+2t=0$ を満たすものなので、$W_1\cap W_2=\langle\transp{(1,1,-1,-1)}\rangle$、$\dim=1$。次元公式 :eqref[eq-linalg-dim-formula] より $\dim(W_1+W_2)=3+2-1=4$、すなわち $W_1+W_2=\R^4$ である。
:::
::::

::::exercise{level=advanced}
$V$ の部分空間 $W_1, W_2$ について、$W_1\cup W_2$ が部分空間であるための必要十分条件は、$W_1\subset W_2$ または $W_2\subset W_1$ であることを示せ。
:::solution
十分性は明らか。必要性を示すため、$W_1\not\subset W_2$ かつ $W_2\not\subset W_1$ と仮定し、$w_1\in W_1\setminus W_2$, $w_2\in W_2\setminus W_1$ をとる。$W_1\cup W_2$ が部分空間なら $w_1+w_2\in W_1\cup W_2$ である。$w_1+w_2\in W_1$ ならば $w_2=(w_1+w_2)-w_1\in W_1$ となり矛盾、$w_1+w_2\in W_2$ ならば $w_1\in W_2$ となり矛盾。よって $W_1\cup W_2$ は部分空間でない。
:::
::::

::::exercise{level=advanced}
相異なる実数 $a_1,\dots,a_k$ に対し、関数 $f_i(t)=e^{a_it}$ は実ベクトル空間 $C(\R)$ において一次独立であることを示せ（微分の基本性質を用いてよい）。
:::hint
$k$ に関する帰納法。$\sum c_ie^{a_it}=0$ に $e^{-a_kt}$ を掛けてから微分する。
:::
:::solution
$k=1$ のときは $e^{a_1t}\ne0$ より明らか。$k-1$ 個のとき成り立つとし、$\sum_{i=1}^k c_ie^{a_it}=0$（すべての $t$ で）とする。両辺に $e^{-a_kt}$ を掛けると $\sum_{i=1}^{k-1}c_ie^{(a_i-a_k)t}+c_k=0$ であり、これを $t$ で微分すると

$$
\sum_{i=1}^{k-1} c_i(a_i-a_k)e^{(a_i-a_k)t}=0.
$$

$a_i-a_k$（$i<k$）は相異なるので、帰納法の仮定から $c_i(a_i-a_k)=0$、$a_i\ne a_k$ より $c_i=0$（$i<k$）である。すると $c_k=0$ も従う。
:::
::::

## 確認テスト

::::quiz
次のうち $\R^2$ の部分空間であるものはどれか。

- [ ] $\{(x,y)^\top\mid y=x+1\}$
- [x] $\{(x,y)^\top\mid 3x-2y=0\}$
- [ ] $\{(x,y)^\top\mid x\ge0\}$
- [ ] $\{(x,y)^\top\mid xy\ge0\}$

:::explanation
原点を通る直線のみが部分空間である。1 番目は原点を含まない。3 番目はスカラー倍（$-1$ 倍）で閉じていない。4 番目は和で閉じていない（$(1,0)+(0,-1)$ など）。
:::
::::

::::quiz{answer="6"}
3 次実対称行列全体 $\mathrm{Sym}_3(\R)$ の次元はいくつか。

:::explanation
$n(n+1)/2=3\cdot4/2=6$。対角成分 3 個と上三角部分 3 個を自由に選べる。
:::
::::

::::quiz{answer="2"}
$\R^5$ の部分空間 $W_1, W_2$ が $\dim W_1=3$, $\dim W_2=4$ を満たすとき、$\dim(W_1\cap W_2)$ のとりうる最小値はいくつか。

:::explanation
次元公式より $\dim(W_1\cap W_2)=3+4-\dim(W_1+W_2)\ge 7-5=2$。等号は $W_1+W_2=\R^5$ のときに成り立つ。
:::
::::

::::quiz
$\dim V=n$ とする。正しいものをすべて選べ。

- [x] $n+1$ 個のベクトルは必ず一次従属である
- [x] $n$ 個の一次独立なベクトルは基底をなす
- [ ] $n-1$ 個のベクトルが一次独立なら、それらは $V$ を生成する
- [x] 任意の部分空間は補空間をもつ

:::explanation
1 は :ref[lem-linalg-steinitz]、2 は :ref[cor-linalg-n-vectors-basis]、4 は :ref[cor-linalg-complement-exists] による。$n-1$ 個のベクトルが張る空間の次元は高々 $n-1$ なので $V$ を生成しない。
:::
::::

::::quiz
$\C$ を $\R$ 上のベクトル空間とみたとき、次元はいくつか。

- [ ] 1
- [x] 2
- [ ] 無限大
- [ ] 体によらず定まらない

:::explanation
$(1, i)$ が $\R$ 上の基底である。$\C$ 上では $\dim_\C\C=1$ であり、次元は係数体に依存する。
:::
::::

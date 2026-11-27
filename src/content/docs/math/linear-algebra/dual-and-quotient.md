---
title: 双対空間と商空間
description: 線形汎関数の空間と双対基底、零化空間、双対写像と転置、二重双対への自然な同型、双線形形式による双対性、商空間とその普遍性、同型定理と完全列を扱う。
notation:
  - symbol: 'V^\vee'
    meaning: V 上の線形汎関数全体（双対空間）
  - symbol: 'e^1,\ldots,e^n'
    meaning: 基底 e_1, …, e_n の双対基底
  - symbol: 'U^\circ'
    meaning: 部分空間 U の零化空間
  - symbol: 'F_\circ'
    meaning: 汎関数の集合 F が共通に消すベクトル全体
  - symbol: 'T^\vee'
    meaning: 線形写像 T の双対写像
  - symbol: '\iota_V'
    meaning: V から二重双対への評価写像
  - symbol: 'V/U'
    meaning: 部分空間 U による商空間
  - symbol: '\pi_U'
    meaning: 商空間への標準射影
  - symbol: '\operatorname{Coker}T'
    meaning: 線形写像 T の余核 W/Img T
---

## この章の目標

- 双対空間と双対基底を定義し、汎関数の係数と基底変換の関係（転置と逆転置）を導く。
- 零化空間を通じて部分空間と連立一次方程式の系を対応させ、核と像の双対関係から行階数と列階数の一致を再証明する。
- 評価写像による有限次元空間と二重双対の自然な同型を示し、非退化な双線形形式による双対性を理解する。
- 商空間を構成して普遍性を証明し、商空間と双対空間の関係、誘導写像の行列を調べる。
- 第一・第二・第三同型定理、部分空間の対応、短完全列の分裂と双対を証明し、次元の計算に使えるようにする。

## 双対空間と双対基底

連立方程式の一行 $a_1x_1+\cdots+a_nx_n$ は、列ベクトルから一つの数を取り出す線形な規則である。
同じ未知数に対する行を足したり定数倍したりできるので、こうした規則もまたベクトル空間をなす。
この空間を元の空間と対等な対象として扱うと、部分空間を方程式で記述すること、行列の転置、連立方程式の可解条件などが一つの枠組みで理解できる。

本章の係数体 $K$ は任意とする。
次元の等式や基底の本数に触れるときは有限次元を仮定し、そのつど明記する。

:::definition[双対空間]{#def-linalg-dual-space}
$K$ 上のベクトル空間 $V$ に対し、線形写像 $f:V\to K$ を $V$ 上の :term[線形汎関数]{en="linear functional" yomi="せんけいはんかんすう"} という。
線形汎関数全体 $V^\vee\defeq\Hom(V,K)$ を、:ref[def-linalg-hom] の和とスカラー倍
$(af+bg)(v)=af(v)+bg(v)$ によりベクトル空間とみなしたものを $V$ の :term[双対空間]{en="dual space" yomi="そうついくうかん"} という。
:::

$V^\vee$ を $V^*$ と書く文献も多い。
本サイトでは内積から定まる随伴 $\adj{T}$ と区別するため、代数的な双対には $\vee$ を用いる。
複素数上でも $f(iv)=if(v)$ であり、汎関数は複素線形で、引数に共役は付かない。

:::example[数ベクトル空間の双対]
$V=K^n$（列ベクトル）とする。
行ベクトル $a=(a_1,\ldots,a_n)$ は $f_a(x)=ax=\sum_ia_ix_i$ により汎関数を定める。
逆に任意の $f\in(K^n)^\vee$ は $a_i=f(e_i)$ とおけば $f(x)=\sum_ix_if(e_i)=f_a(x)$ となる。
$a\mapsto f_a$ は線形かつ全単射なので、$(K^n)^\vee$ は行ベクトルの空間 $K^{1\times n}$ と同一視できる。
「列ベクトルの双対は行ベクトル」というこの同一視が、以下で転置が現れる理由である。
:::

:::theorem[双対基底]{#thm-linalg-dual-basis}
$e_1,\ldots,e_n$ が $V$ の基底なら、$e^i(e_j)=\delta_{ij}$ を満たす汎関数 $e^1,\ldots,e^n$ が一意に存在し、$V^\vee$ の基底となる。
これを $(e_j)$ の :term[双対基底]{en="dual basis" yomi="そうついきてい"} という。
特に $\dim V^\vee=\dim V$ である。
さらに任意の $v\in V$、$f\in V^\vee$ について

$$
v=\sum_{i=1}^n e^i(v)\,e_i,\qquad f=\sum_{i=1}^n f(e_i)\,e^i \label{eq-linalg-dual-expansion}
$$

が成り立つ。
:::

:::proof
$v=\sum_j x_je_j$ の表示は一意なので、$e^i(v)=x_i$ により写像 $e^i:V\to K$ が定まる。
$v+w$ の座標は座標の和、$av$ の座標は座標の $a$ 倍だから $e^i$ は線形であり、$e^i(e_j)=\delta_{ij}$ を満たす。
線形写像は基底上の値で決まる（:ref[prop-linalg-map-on-basis]）ので、この条件を満たす汎関数は $e^i$ のみである。

$v=\sum_jx_je_j$ なら $x_i=e^i(v)$ なので、:eqref[eq-linalg-dual-expansion] の第一式は座標の定義そのものである。
任意の $f$ に対して

$$
f(v)=\sum_j x_j f(e_j)=\sum_j f(e_j)e^j(v)=\Bigl(\sum_j f(e_j)e^j\Bigr)(v)
$$

であるから第二式が成り立ち、$e^1,\ldots,e^n$ は $V^\vee$ を生成する。
$\sum_ja_je^j=0$ を $e_i$ に適用すると $a_i=0$ となるので一次独立でもある。
:::

双対基底 $e^i$ は「第 $i$ 座標を読む汎関数」である。
注意すべきは、$e^i$ は $e_i$ だけでなく基底全体に依存することである。
$e_1$ を固定して $e_2$ だけを取り替えても $e^1$ は変わりうる（次の例）。

:::example[長さが 1 でない直交基底の座標を読む]
$\R^2$ の基底 $b_1=(1,1),\ b_2=(1,-1)$ に対し $(x,y)=\frac{x+y}{2}b_1+\frac{x-y}{2}b_2$ である。
したがって新しい係数を読み取る二つの汎関数は

$$
b^1(x,y)=\frac{x+y}{2},\qquad b^2(x,y)=\frac{x-y}{2}
$$

である。
$f(x,y)=2x+3y$ は $f(b_1)=5,\ f(b_2)=-1$ なので :eqref[eq-linalg-dual-expansion] より $f=5b^1-b^2$ と書ける。
一方、基底 $b_1,\ e_2=(0,1)$ に対しては $(x,y)=xb_1+(y-x)e_2$ なので、$b_1$ に対応する双対基底元は $x$ であり、$b^1$ とは異なる。
:::

### 基底変換と双対基底

基底を取り替えたとき双対基底がどう変わるかを、行列で正確に述べておく。
:ref[def-linalg-change-of-basis] と同じく、新しい基底ベクトルの旧座標を列に並べた行列を用いる。

:::proposition[双対基底の変換則]{#prop-linalg-dual-basis-change}
$V$ の基底 $(e_i)$ と $(b_j)$ が $b_j=\sum_iS_{ij}e_i$ で結ばれているとする（$S=P_{\mathcal E\leftarrow\mathcal B}$ は正則）。
このとき双対基底は

$$
b^j=\sum_i (S^{-1})_{ji}\,e^i
$$

で結ばれる。
すなわち新しい双対基底元の旧双対基底での係数を列に並べた行列は $\transp{(S^{-1})}$ である。
また汎関数 $f$ の係数ベクトル $a=(f(e_i))_i$、$a'=(f(b_j))_j$ の間には $a'=\transp{S}a$ が成り立つ。
:::

:::proof
$b^j=\sum_iL_{ji}e^i$ と書くと、$b^j(b_k)=\sum_iL_{ji}e^i\bigl(\sum_lS_{lk}e_l\bigr)=\sum_iL_{ji}S_{ik}=(LS)_{jk}$ である。
これが $\delta_{jk}$ となる条件は $LS=I$、すなわち $L=S^{-1}$ である。
後半は $a'_j=f(b_j)=\sum_iS_{ij}f(e_i)=(\transp{S}a)_j$ による。
:::

ベクトルの座標は $x'=S^{-1}x$ と「逆」に変わり、汎関数の係数は $a'=\transp{S}a$ と変わる。
したがって値 $f(v)=\transp{a}x=\transp{(a')}x'$ は基底によらない。
物理で言う反変成分と共変成分の区別は、この変換則の違いである。

:::example[基底変換の具体例]
$S=\begin{pmatrix}1&1\\0&1\end{pmatrix}$、つまり $b_1=e_1,\ b_2=e_1+e_2$ とする。
$S^{-1}=\begin{pmatrix}1&-1\\0&1\end{pmatrix}$ の行を読むと $b^1=e^1-e^2,\ b^2=e^2$、すなわち $b^1(x,y)=x-y,\ b^2(x,y)=y$ である。
実際 $b^1(b_2)=1-1=0$ である。
汎関数 $f=e^1$（第一座標）の新しい係数は $\transp{S}(1,0)=(1,1)$ であり、$f=b^1+b^2$ と確かめられる。
$b_1=e_1$ は変わらないのに $b^1\ne e^1$ である点に注意する。
:::

:::example[Lagrange 補間と数値積分の重み]{#exm-linalg-lagrange-dual}
$V=K[t]_{\le n}$ とし、相異なる $a_0,\ldots,a_n\in K$ を取る。

$$
\ell_i(t)=\prod_{j\ne i}\frac{t-a_j}{a_i-a_j}\qquad(i=0,\ldots,n)
$$

は $\ell_i(a_k)=\delta_{ik}$ を満たす。
評価汎関数 $\mathrm{ev}_{a_k}(p)=p(a_k)$ を $\sum_ic_i\ell_i=0$ に適用すると $c_k=0$ なので $\ell_0,\ldots,\ell_n$ は一次独立であり、$\dim V=n+1$ より基底である。
その双対基底がちょうど $\mathrm{ev}_{a_0},\ldots,\mathrm{ev}_{a_n}$ である。
:eqref[eq-linalg-dual-expansion] の第一式は Lagrange の補間公式 $p=\sum_ip(a_i)\ell_i$ にほかならない。

第二式からは数値積分の公式が得られる。
$K=\R$、$n=2$、$a=-1,0,1$ とし、$I(p)=\int_{-1}^1p(t)\,\dd t$ とおく。
$\ell_{-1}=t(t-1)/2,\ \ell_0=1-t^2,\ \ell_1=t(t+1)/2$ の積分はそれぞれ $1/3,\ 4/3,\ 1/3$ なので

$$
\int_{-1}^1p(t)\,\dd t=\frac13p(-1)+\frac43p(0)+\frac13p(1)\qquad(\deg p\le2)
$$

である（Simpson の公式）。
$t^3$ では両辺とも $0$ なので、実は三次以下で正確である。
:::

### 零化空間

部分空間を「どの汎関数で消えるか」で記述する。

:::definition[零化空間]{#def-linalg-annihilator}
部分集合 $U\subset V$ に対し、$U$ の :term[零化空間]{en="annihilator" yomi="れいかくうかん"} を次で定める。

$$
U^\circ=\{f\in V^\vee:f(u)=0\ (\forall u\in U)\}
$$

逆に $F\subset V^\vee$ に対し、$F$ のすべての元で消えるベクトル全体を次で表す。

$$
F_\circ=\{v\in V:f(v)=0\ (\forall f\in F)\}
$$

:::

$U^\circ$ は $V^\vee$ の部分空間、$F_\circ$ は $V$ の部分空間である（消える条件は線形な条件だから）。
$F_\circ=\bigcap_{f\in F}\Ker f$ であり、$F$ を方程式系とみればその解空間である。
また $U^\circ=(\spn U)^\circ$ である。

:::proposition[零化空間の次元と復元]{#prop-linalg-annihilator}
$\dim V=n$ とし、$U\subset V$、$F\subset V^\vee$ を部分空間とする。
1. $\dim U^\circ=n-\dim U$ であり、$(U^\circ)_\circ=U$ である。
2. $\dim F_\circ=n-\dim F$ であり、$(F_\circ)^\circ=F$ である。

したがって $U\mapsto U^\circ$ と $F\mapsto F_\circ$ は、$V$ の部分空間と $V^\vee$ の部分空間の間の互いに逆な全単射であり、包含の向きを逆にする。
:::

:::proof
1. $U$ の基底 $e_1,\ldots,e_k$ を :ref[thm-linalg-basis-extension] により $V$ の基底 $e_1,\ldots,e_n$ に延長し、双対基底を $e^i$ とする。
$f=\sum_ia_ie^i$ が $U$ を消す条件は $f(e_j)=a_j=0\ (j\le k)$ だから、$U^\circ=\spn\{e^{k+1},\ldots,e^n\}$ であり、$\dim U^\circ=n-k$ である。
$v=\sum_jx_je_j$ がこれらすべてで消える条件は $x_{k+1}=\cdots=x_n=0$、すなわち $v\in U$ である。

2. $F$ の基底 $f_1,\ldots,f_r$ を取り、写像 $\Phi:V\to K^r,\ v\mapsto(f_1(v),\ldots,f_r(v))$ を考える。
$F_\circ=\Ker\Phi$ である。
$\Phi$ が全射でないとすると、像は $K^r$ の真の部分空間であり、1 を $K^r$ に適用すれば像を消す非零の汎関数 $y\mapsto\sum_ic_iy_i$ が存在する。
すると $\sum_ic_if_i(v)=0$ がすべての $v$ で成り立ち、$\sum_ic_if_i=0$ は $f_i$ の独立性に反する。
よって $\Phi$ は全射で、:ref[thm-linalg-rank-nullity] より $\dim F_\circ=n-r$ である。
定義から $F\subset(F_\circ)^\circ$ であり、1 より $\dim(F_\circ)^\circ=n-(n-r)=r=\dim F$ なので両者は一致する。

最後の主張は 1, 2 から直ちに従う。
包含の逆転は定義から明らかである（消すべきベクトルが増えれば条件が増える）。
:::

:::example[部分空間を方程式で指定する]
$\R^3$ で $U=\spn\{(1,1,0),(0,1,1)\}$ とする。
行 $(a,b,c)$ が $U$ を消す条件は $a+b=b+c=0$ だから、$U^\circ$ は行 $(1,-1,1)$ が張る一次元空間である。
したがって $U=(U^\circ)_\circ$ はただ一つの方程式 $x-y+z=0$ の解空間である。
ここでは内積や角度を選ぶ必要はない。
:::

:::proposition[和と共通部分の零化空間]{#prop-linalg-annihilator-lattice}
$U,W\subset V$ を部分空間とする。
1. $U\subset W$ なら $W^\circ\subset U^\circ$ である。
2. $(U+W)^\circ=U^\circ\cap W^\circ$ である。
3. $V$ が有限次元なら $(U\cap W)^\circ=U^\circ+W^\circ$ である。
:::

:::proof
1 は定義から明らかである。
2. $f$ が $U+W$ を消すことと、$U$ と $W$ の両方を消すことは同値である（$f(u+w)=f(u)+f(w)$）。

3. $f\in U^\circ,\ g\in W^\circ$ なら $f+g$ は $U\cap W$ を消すので、$U^\circ+W^\circ\subset(U\cap W)^\circ$ である。
逆向きは次元で示す。
2 を $V^\vee$ の部分空間 $U^\circ,W^\circ$ と $F\mapsto F_\circ$ について同様に考えると $(U^\circ+W^\circ)_\circ=(U^\circ)_\circ\cap(W^\circ)_\circ=U\cap W$ である（:ref[prop-linalg-annihilator] の 1）。
両辺の零化空間を取り、:ref[prop-linalg-annihilator] の 2 を使えば $U^\circ+W^\circ=(U\cap W)^\circ$ を得る。
:::

### 超平面と汎関数の一次結合

$n$ 次元空間の $n-1$ 次元部分空間を :term[超平面]{en="hyperplane" yomi="ちょうへいめん"} という（無限次元では余次元 $1$、すなわち $\dim V/H=1$ の部分空間をそう呼ぶ）。

:::proposition[超平面と非零汎関数]{#prop-linalg-hyperplane}
$\dim V=n\ge1$ とする。
非零の $f\in V^\vee$ の核は超平面であり、逆に任意の超平面 $H$ は $H=\Ker f$ と書ける。
さらに $\Ker f=\Ker g$（$f,g\ne0$）なら $g=cf$ となる $c\ne0$ がある。
:::

:::proof
$f\ne0$ なら像は $K$ の非零部分空間なので $K$ 全体であり、次元定理より $\dim\Ker f=n-1$ である。
超平面 $H$ に対しては :ref[prop-linalg-annihilator] より $\dim H^\circ=1$ なので、非零の $f\in H^\circ$ を取れば $H\subset\Ker f$、次元が等しいので $H=\Ker f$ である。
$\Ker f=\Ker g=H$ なら $f,g$ はともに一次元空間 $H^\circ$ の非零元なので互いに非零定数倍である。
:::

最後の主張は次の命題の特別な場合である。
この命題は無限次元でも成り立ち、Lagrange の未定乗数法や関数解析で繰り返し使われる。

:::proposition[汎関数が一次結合になる条件]{#prop-linalg-functional-combination}
$V$ を任意のベクトル空間とし、$f_1,\ldots,f_k,g\in V^\vee$ とする。
$g\in\spn\{f_1,\ldots,f_k\}$ であるための必要十分条件は

$$
\Ker f_1\cap\cdots\cap\Ker f_k\subset\Ker g
$$

が成り立つことである。
:::

:::proof
必要性は明らかである。
十分性を示す。
$\Phi:V\to K^k,\ \Phi(v)=(f_1(v),\ldots,f_k(v))$ とおくと、仮定は $\Ker\Phi\subset\Ker g$ である。
$\Img\Phi$ 上の写像 $h$ を $h(\Phi(v))=g(v)$ で定める。
$\Phi(v)=\Phi(v')$ なら $v-v'\in\Ker\Phi\subset\Ker g$ なので値は代表によらず、$h$ は線形である。
$\Img\Phi$ は有限次元空間 $K^k$ の部分空間なので、その基底を $K^k$ の基底に延長し、追加した基底上で $0$ と定めれば $h$ を $K^k$ 全体の線形汎関数に延長できる。
その延長は $h(y)=\sum_ic_iy_i$ の形なので、$g(v)=h(\Phi(v))=\sum_ic_if_i(v)$ である。
:::

:::corollary[方程式の本数と余次元]{#cor-linalg-codimension}
$\dim V=n$、$U\subset V$ を部分空間とし、$r=n-\dim U$ とする。
$U^\circ$ の基底 $f_1,\ldots,f_r$ に対して $U=\Ker f_1\cap\cdots\cap\Ker f_r$ である。
また $U$ を $k$ 本の一次方程式の共通解として表すには $k\ge r$ が必要である。
:::

:::proof
前半は $(U^\circ)_\circ=U$ そのものである。
$U=\bigcap_{i=1}^k\Ker g_i$ なら、$U^\circ$ の任意の元 $f$ は $\bigcap\Ker g_i\subset\Ker f$ を満たすので :ref[prop-linalg-functional-combination] より $g_i$ の一次結合である。
よって $r=\dim U^\circ\le k$ である。
:::

独立な条件を数えることと、消した方向の数を数えることが双対になっている。

:::example[独立な制約だけを取り出す]
$\R^4$ の未知数 $(x,y,z,w)$ に $x-y=0,\ z-w=0,\ x-y+z-w=0$ を課す。
最後の行は初めの二行の和なので、共通の解は $(s,s,t,t)$ の形の二次元空間 $U$ である。
$U^\circ$ は $(1,-1,0,0),(0,0,1,-1)$ の張る二次元空間であり、:ref[cor-linalg-codimension] により二本の方程式が必要かつ十分である。
例えば $2x-2y-z+w=0$ を追加しても、対応する行が $U^\circ$ に入るので解集合は変わらない。
:::

## 双対写像と転置

写像 $T:V\to W$ の出力を汎関数 $g$ で測った値 $g(Tv)$ は、入力 $v$ の汎関数とみなせる。
こうして汎関数を運ぶ向きは $T$ の矢印と逆になる。

:::definition[双対写像]{#def-linalg-dual-map}
線形写像 $T:V\to W$ の :term[双対写像]{en="dual map" yomi="そうついしゃぞう"}（転置写像ともいう）を

$$
T^\vee:W^\vee\to V^\vee,\qquad T^\vee g=g\circ T
$$

で定める。
$g\circ T$ は線形写像の合成なので $V^\vee$ の元であり、$(ag+bh)\circ T=a(g\circ T)+b(h\circ T)$ より $T^\vee$ は線形である。
:::

:::proposition[双対写像の性質と転置行列]{#prop-linalg-dual-transpose}
1. $(aS+bT)^\vee=aS^\vee+bT^\vee$、$(ST)^\vee=T^\vee S^\vee$、$\id_V^\vee=\id_{V^\vee}$ が成り立つ（$ST$ は合成が定義される場合）。
2. $T$ が同型なら $T^\vee$ も同型で $(T^\vee)^{-1}=(T^{-1})^\vee$ である。
3. 有限次元で、基底 $(e_j)$, $(w_i)$ に関する $T$ の表現行列が $A$ なら、双対基底 $(w^i)$, $(e^j)$ に関する $T^\vee$ の表現行列は $\transp{A}$ である。
:::

:::proof
1. $g\circ(aS+bT)=a(g\circ S)+b(g\circ T)$、$g\circ(S\circ T)=(g\circ S)\circ T$、$g\circ\id=g$ による。
2. 1 より $T^\vee(T^{-1})^\vee=(T^{-1}T)^\vee=\id$、$(T^{-1})^\vee T^\vee=(TT^{-1})^\vee=\id$ である。
3. $Te_j=\sum_ia_{ij}w_i$ に $w^i$ を適用すると $(T^\vee w^i)(e_j)=w^i(Te_j)=a_{ij}$ である。
:eqref[eq-linalg-dual-expansion] より $T^\vee w^i=\sum_ja_{ij}e^j$ なので、$T^\vee$ の表現行列の第 $i$ 列は $(a_{i1},\ldots,a_{in})$、すなわち $A$ の第 $i$ 行を立てたものである。
:::

数ベクトル空間で汎関数を行ベクトルとみると、$T^\vee$ は行ベクトル $y$ を $yA$ に移す写像である。
これを列ベクトルで書いたものが $\transp{A}\transp{y}$ である。
転置の公式 $\transp{(AB)}=\transp{B}\transp{A}$（:ref[prop-linalg-transpose-laws]）は 1 の $(ST)^\vee=T^\vee S^\vee$ の行列版にほかならない。

:::example[複素数でも転置する]
$T:\C^2\to\C,\ T(x,y)=x+iy$ の行列は $(1\ \ i)$ である。
標準基底で $T^\vee$ の行列は $\transp{(1,i)}$ である。
一方、標準内積についての随伴写像（:ref[thm-linalg-adjoint]）の行列は共役転置 $\transp{(1,-i)}$ である。
例えば $g(z)=iz$ を引き戻すと $T^\vee g(x,y)=i(x+iy)=ix-y$ となる。
双対写像は内積なしに定まり、共役を含まない。
:::

内積空間では、双対写像と随伴写像は Riesz の表現（:ref[thm-linalg-riesz-finite]）によって結ばれる。

:::proposition[随伴と双対写像]{#prop-linalg-adjoint-dual}
$V,W$ を有限次元内積空間とし、$R_V:V\to V^\vee,\ R_V(z)=\inner{\cdot}{z}$ とおく（$W$ も同様）。
$R_V$ は全単射で、$R_V(az+bz')=\conj{a}R_V(z)+\conj{b}R_V(z')$（共役線形）を満たす。
線形写像 $T:V\to W$ に対して

$$
\adj{T}=R_V^{-1}\circ T^\vee\circ R_W
$$

が成り立つ。
:::

:::proof
第一引数についての線形性から $R_V(z)\in V^\vee$ であり、第二引数についての共役線形性から $R_V$ は共役線形である。
:ref[thm-linalg-riesz-finite] は $R_V$ が全単射であることを述べている。
$y\in W$ に対し、$T^\vee(R_Wy)$ は $x\mapsto\inner{Tx}{y}=\inner{x}{\adj{T}y}$ という汎関数なので $R_V(\adj{T}y)$ に等しい。
:::

共役線形な $R$ を二度挟むため、正規直交基底では $T^\vee$ の行列 $\transp{A}$ が $\adj{T}$ の行列 $\adj{A}=\transp{\conj{A}}$ に変わる。
実内積空間では $R_V$ は線形な同型であり、$\adj{T}$ と $T^\vee$ は同一視できる。

:::theorem[核と像の双対関係]{#thm-linalg-dual-kernel}
線形写像 $T:V\to W$ について次が成り立つ。

$$
\Ker T^\vee=(\Img T)^\circ,\qquad \Img T^\vee\subset(\Ker T)^\circ .
$$

$V,W$ が有限次元なら第二式は等号 $\Img T^\vee=(\Ker T)^\circ$ となり、さらに $T$ が単射であることと $T^\vee$ が全射であること、$T$ が全射であることと $T^\vee$ が単射であることはそれぞれ同値である。
:::

:::proof
$T^\vee g=g\circ T=0$ は、$g$ が $T$ の像を消すという条件そのものである。
また $g\circ T$ は $\Ker T$ を必ず消すので第二式の包含が成り立つ。

有限次元で逆の包含を示す。
$f$ が $\Ker T$ を消すとする。
$\Img T$ 上の写像 $g_0$ を $g_0(Tv)=f(v)$ で定める。
$Tv=Tv'$ なら $v-v'\in\Ker T$ なので $f(v)=f(v')$ であり、$g_0$ は矛盾なく定まって線形である。
$\Img T$ の基底を $W$ の基底に延長し、追加した基底上で $0$ と置いて $g_0$ を $g\in W^\vee$ に延長すれば、$T^\vee g=g\circ T=f$ である。

単射・全射の同値性は次のように従う。
$T^\vee$ が単射 $\iff(\Img T)^\circ=0\iff\Img T=W$（:ref[prop-linalg-annihilator] の次元公式）である。
$T^\vee$ が全射 $\iff(\Ker T)^\circ=V^\vee\iff\Ker T=0$ である。
:::

:::corollary[双対写像の階数]{#cor-linalg-dual-rank}
有限次元空間の間の線形写像 $T$ について $\rank T^\vee=\rank T$ である。
特に任意の行列 $A$ について $\rank\transp{A}=\rank A$、すなわち行階数と列階数は等しい。
:::

:::proof
$\dim V=n$ とすると、:ref[thm-linalg-dual-kernel] と :ref[prop-linalg-annihilator] より
$\rank T^\vee=\dim(\Ker T)^\circ=n-\dim\Ker T$ であり、次元定理 :ref[thm-linalg-rank-nullity] よりこれは $\rank T$ に等しい。
行列については :ref[prop-linalg-dual-transpose] の 3 により $T^\vee$ の行列が $\transp{A}$ であることを使えばよい。
:::

この証明は行基本変形を使わない、:ref[thm-linalg-row-column-rank] の別証明になっている。

### 直和の双対と射影の引き戻し

:::proposition[直和の双対]{#prop-linalg-dual-direct-sum}
$V=U\oplus W$ とし、$p_U,p_W$ を直和分解に沿った射影とする。
写像 $V^\vee\to U^\vee\times W^\vee,\ f\mapsto(f|_U,f|_W)$ は同型であり、その逆は $(g,h)\mapsto g\circ p_U+h\circ p_W$ である。
この同型により $W^\circ$ は $U^\vee$ と、$U^\circ$ は $W^\vee$ と同一視される。
:::

:::proof
$v=p_Uv+p_Wv$ なので $f=f|_U\circ p_U+f|_W\circ p_W$ であり、二つの写像の合成 $V^\vee\to U^\vee\times W^\vee\to V^\vee$ は恒等写像である。
逆向きの合成では $(g\circ p_U+h\circ p_W)|_U=g$（$p_U|_U=\id_U,\ p_W|_U=0$）などから恒等写像を得る。
$f\in W^\circ$ は $f|_W=0$ と同値なので $(g,0)$ の形の元に対応する。
:::

射影 $p_U:V\to U$ の双対 $p_U^\vee:U^\vee\to V^\vee$ は $g\mapsto g\circ p_U$ で、その像はちょうど $W^\circ$ である。
同じ $U$ でも補空間 $W$ を取り替えれば $p_U$ も $W^\circ$ も変わることに注意する。

:::example[差分の転置と境界での係数]
三点の値から隣り合う値の差を作る写像 $T:\R^3\to\R^2$ を考える。

$$
T(x_1,x_2,x_3)=(x_2-x_1,x_3-x_2)
$$

出力に係数 $(a,b)$ を掛けて測った量を整理すると

$$
a(x_2-x_1)+b(x_3-x_2)=-ax_1+(a-b)x_2+bx_3
$$

となる。
したがって $T^\vee$ の像は係数 $(-a,a-b,b)$ の行全体であり、その成分和は零である。
逆に $(r,s,t)$ の和が零なら $a=-r,\ b=t$ とすれば中間の係数も一致する。
$\Ker T$ は定数列 $(c,c,c)$ の全体であり、:ref[thm-linalg-dual-kernel] の $\Img T^\vee=(\Ker T)^\circ$ が「係数和が零」という条件として現れたのである。
これは部分和分（離散的な部分積分）の公式の有限次元版である。
:::

### 方程式の整合条件を出力側で調べる

:::corollary[Fredholm の択一定理（有限次元）]{#cor-linalg-fredholm-alternative}
$A\in K^{m\times n}$、$b\in K^m$ とする。
$Ax=b$ が解をもつための必要十分条件は、$\transp{A}y=0$ を満たすすべての $y\in K^m$ について $\transp{y}b=0$ となることである。
特に「$Ax=b$ がすべての $b$ に対して解をもつ」ことと「$\transp{A}y=0$ が自明な解しかもたない」ことは同値である。
:::

:::proof
$y\in K^m$ を汎関数 $g(w)=\transp{y}w$ とみると、$\transp{A}y=0$ は $g\circ A=0$、すなわち $g\in\Ker A^\vee=(\Img A)^\circ$ と同値である。
可解性 $b\in\Img A$ は :ref[prop-linalg-annihilator] の $((\Img A)^\circ)_\circ=\Img A$ により、$(\Img A)^\circ$ のすべての元が $b$ を消すことと同値である。
後半は $b$ を任意にした場合で、$(\Img A)^\circ=0$ と $\Img A=K^m$ の同値性である。
:::

行列でいえば、$A$ の行の任意の一次関係（$\transp{y}A=0$）に右辺も従うこと（$\transp{y}b=0$）が必要であり、それで十分でもある。
これは :ref[thm-linalg-rank-criterion] の拡大係数行列による判定を、消去操作とは別の側から見たものである。
必要な独立条件の数は $\dim(\Img A)^\circ=m-\rank A$ である。

:::example[左側の核が与える測定値の制約]
二つの実数 $x,y$ から三つの測定値 $(x+y,x-y,2x)$ を作る。
測定行列 $A=\begin{pmatrix}1&1\\1&-1\\2&0\end{pmatrix}$ では、最初の二行の和が最後の行である。
階数は $2$ なので $\Ker\transp{A}$ は一次元で、$(1,1,-1)$ が張る。
したがってデータ $b$ に解がある条件は $b_1+b_2-b_3=0$ である。

$b=(4,2,6)$ ではこの条件が成立し、最初の二式を解くと $(x,y)=(3,1)$ で第三式も満たす。
$b=(4,2,5)$ では $b_1+b_2-b_3=1\ne0$ であり、この一つの数だけで解がないことが証明できる。
:::

## 二重双対

$V$ と $V^\vee$ は有限次元なら同じ次元をもつので同型である。
しかし同型写像を一つ指定するには、基底や内積などの選択が要る（下の例を見よ）。
一方、$V^\vee$ の双対 $V^{\vee\vee}=(V^\vee)^\vee$ への写像は、ベクトルを汎関数に「入力する」だけで何の選択もなしに作れる。

:::theorem[評価写像と自然性]{#thm-linalg-double-dual}
$V^{\vee\vee}=(V^\vee)^\vee$ を $V$ の :term[二重双対]{en="double dual" yomi="にじゅうそうつい"} と呼ぶ。
:term[評価写像]{en="evaluation map" yomi="ひょうかしゃぞう"}

$$
\iota_V:V\longrightarrow V^{\vee\vee},\qquad \iota_V(v)(f)=f(v)
$$

は線形であり、$V$ が有限次元なら同型である。
さらに任意の線形写像 $T:V\to W$ に対して

$$
T^{\vee\vee}\circ\iota_V=\iota_W\circ T \label{eq-linalg-double-dual-natural}
$$

が成り立つ。
:::

:::proof
各 $v$ について $f\mapsto f(v)$ は $f$ について線形なので $\iota_V(v)\in V^{\vee\vee}$ である。
$\iota_V(av+bw)(f)=af(v)+bf(w)$ がすべての $f$ で成り立つので $\iota_V$ は線形である。

有限次元とする。
$v\ne0$ なら $v$ を基底 $v=e_1,e_2,\ldots,e_n$ に延長でき、双対基底の $e^1$ は $e^1(v)=1$ を満たす。
よって $\iota_V(v)\ne0$ であり $\iota_V$ は単射である。
:ref[thm-linalg-dual-basis] を二度使うと $\dim V^{\vee\vee}=\dim V^\vee=\dim V$ なので、:ref[cor-linalg-equal-dimension-map] より全射でもある。

最後に $g\in W^\vee$ に対し

$$
\bigl(T^{\vee\vee}\iota_V(v)\bigr)(g)=\iota_V(v)(T^\vee g)=(g\circ T)(v)=g(Tv)=\iota_W(Tv)(g)
$$

である。
:::

:eqref[eq-linalg-double-dual-natural] は、基底を選ばずに定めた対応がすべての線形写像と両立することを表す。
これを $\iota$ の :term[自然性]{en="naturality" yomi="しぜんせい"} という。
有限次元では $\iota_V$ によって $V=V^{\vee\vee}$ と同一視することが多く、この同一視の下で $T^{\vee\vee}=T$ となる。
行列でいえば $\transp{(\transp{A})}=A$ である。

:::example[二重双対の元をベクトルに戻す]
$V=\R^2$ で $L(f)=2f(e_1)-3f(e_2)$ とする。
$f$ は線形なので $L(f)=f(2e_1-3e_2)$ であり、$\iota_V^{-1}L=(2,-3)$ である。
一般に有限次元では $\iota_V^{-1}L=\sum_iL(e^i)e_i$ である。
実際、$f=\sum_if(e_i)e^i$ を代入すると $L(f)=\sum_if(e_i)L(e^i)=f\bigl(\sum_iL(e^i)e_i\bigr)$ となる。
:::

### 自然な同型と基底に依存する同型

基底 $(e_i)$ を一つ選ぶと、$e_i\mapsto e^i$ を線形に延長して同型 $\phi:V\to V^\vee$ が得られる。
しかしこれは基底の選び方に依存する。
実一次元空間 $V=\R e$ で、基底 $e$ に対する規則は $e\mapsto e^1$（$e^1(e)=1$）である。
新しい基底 $b=2e$ を取ると $b^1=e^1/2$ であり、規則 $b\mapsto b^1$ を線形に延長すると $e=b/2\mapsto b^1/2=e^1/4$ となる。
二つの規則は $e$ の行き先がすでに異なる。

これに対して $\iota_V$ はどの基底で計算しても同じ写像である。
実際、任意の基底 $(e_i)$ に対して $\iota_V(e_i)(e^j)=e^j(e_i)=\delta_{ji}$ なので、$\iota_V$ は基底 $(e_i)$ を「双対基底の双対基底」に送る。
「基底を選ぶと得られる」写像のうち、選び方によらないものがこれである。

:::example[内積を選んだ場合との違い]
複素内積空間で $z$ に $R_V(z)=\inner{\cdot}{z}$ を対応させると、$R_V(iz)=-iR_V(z)$ である（:ref[prop-linalg-adjoint-dual]）。
よって内積が与える $V\to V^\vee$ は共役線形であり、内積の選び方にも依存する。
これに対し $\iota_V(iz)(f)=f(iz)=if(z)$ なので、評価による $V\to V^{\vee\vee}$ は複素線形で、追加の構造を必要としない。
:::

零化空間を二重双対の中で書けば、:ref[prop-linalg-annihilator] の $(U^\circ)_\circ=U$ は $(U^\circ)^\circ=\iota_V(U)$ と同値である。
左辺は $V^{\vee\vee}$ の部分空間であり、$V$ の部分空間と等号で結ぶにはこの同一視を挟む必要がある。

### 双線形形式による双対性

評価 $(f,v)\mapsto f(v)$ は二変数のそれぞれについて線形である。
この性質をもつ一般の対を考えると、実内積や行列のトレースによる形式も同じ枠組みで扱える。複素内積では第 2 引数を固定した場合にだけ複素線形な汎関数が得られる。

:::definition[双線形形式と非退化性]{#def-linalg-bilinear-form}
$V,W$ を $K$ 上の空間とする。
写像 $B:V\times W\to K$ が、$w$ を固定するごとに $v$ について線形、$v$ を固定するごとに $w$ について線形であるとき、:term[双線形形式]{en="bilinear form" yomi="そうせんけいけいしき"} という。
$B$ は二つの線形写像

$$
B_L:V\to W^\vee,\ B_L(v)=B(v,\cdot),\qquad B_R:W\to V^\vee,\ B_R(w)=B(\cdot,w)
$$

を定める。
$\Ker B_L$ を左核、$\Ker B_R$ を右核という。
左核と右核がともに $0$ のとき $B$ は :term[非退化]{en="nondegenerate" yomi="ひたいか"} であるという。
:::

$V$ の基底 $(e_i)$ と $W$ の基底 $(f_j)$ を選ぶと、$M=(B(e_i,f_j))_{i,j}$ を $B$ の行列（Gram 行列）といい、$B(v,w)=\transp{x}My$（$x,y$ は座標）である。
基底を $e'=eS$、$f'=fT$（$S,T$ は基底変換行列）に取り替えると、行列は $\transp{S}MT$ に変わる。
$V=W$ で同じ基底を使えば $M\mapsto\transp{S}MS$ であり、実対称形式では [二次形式の章](/math/linear-algebra/spectral-theorem/) の合同変換と一致する。複素数上では双線形形式には転置を用い、Hermite 形式の合同変換には共役転置を用いる。

:::proposition[完全な対]{#prop-linalg-perfect-pairing}
$V,W$ を有限次元とし、$B:V\times W\to K$ を双線形形式、$M$ をその行列とする。
1. $B_R=B_L^\vee\circ\iota_W$ であり、$B_L$ と $B_R$ の表現行列は双対基底について $\transp{M}$ と $M$ である。
2. $\rank B_L=\rank B_R=\rank M$ である。
3. $B$ が非退化であるための必要十分条件は、$\dim V=\dim W$ かつ $M$ が正則であることである。このとき $B_L,B_R$ はともに同型である。
4. $\dim V=\dim W$ なら、左核が $0$ であることと右核が $0$ であることは同値である。
:::

:::proof
1. $w\in W,\ v\in V$ に対し $\bigl(B_L^\vee(\iota_Ww)\bigr)(v)=\iota_W(w)(B_Lv)=B(v,w)=B_R(w)(v)$ である。
$B_L(e_i)=\sum_jB(e_i,f_j)f^j$ なので、$B_L$ の表現行列の第 $i$ 列は $M$ の第 $i$ 行であり、行列は $\transp{M}$ である。
同様に $B_R(f_j)=\sum_iB(e_i,f_j)e^i$ なので $B_R$ の行列は $M$ である。

2. 1 と :ref[cor-linalg-dual-rank] による。

3. $B_L$ が単射 $\iff\rank M=\dim V$、$B_R$ が単射 $\iff\rank M=\dim W$ である。
両方が成り立つことは $\dim V=\dim W=\rank M$、すなわち $M$ が正方かつ正則であることと同値である。
このとき $B_L,B_R$ は同じ次元の空間の間の単射なので同型である。

4. 3 の証明から、$\dim V=\dim W$ ならどちらの条件も $\rank M=\dim V$ と同値である。
:::

非退化な $B$ があれば、$V$ の基底 $(e_i)$ に対して $B(e_i,f_j)=\delta_{ij}$ を満たす $W$ の基底 $(f_j)$ がただ一つ存在する。
実際 $f_j=B_R^{-1}(e^j)$ とすればよい。
これを $B$ に関する双対基底という。

:::example[逆格子]
$\R^3$ の標準内積 $B(v,w)=v\cdot w$ は非退化である。
結晶の基本並進ベクトル $a_1=(1,0,0),\ a_2=(1,1,0),\ a_3=(0,1,1)$ に対し、$a_i\cdot b_j=\delta_{ij}$ を満たす $b_j$ を固体物理では逆格子ベクトル（の $2\pi$ 倍を除いたもの）という。
$a_i$ を列とする行列を $A$ とすれば、$b_j$ を列とする行列は $(\transp{A})^{-1}$ であり、計算すると

$$
b_1=(1,-1,1),\qquad b_2=(0,1,-1),\qquad b_3=(0,0,1)
$$

となる。
三次元では $b_1=\frac{a_2\times a_3}{a_1\cdot(a_2\times a_3)}$ などの公式もあり、ここでは $a_1\cdot(a_2\times a_3)=1$ なので $b_1=a_2\times a_3=(1,-1,1)$ と一致する。
:::

:::example[トレースによる対]
$M_n(K)$ 上の $B(A,X)=\tr(AX)$ は双線形で、行列単位 $E_{ij}$ について $\tr(E_{ij}E_{kl})=\delta_{jk}\delta_{il}$ である。
したがって基底 $(E_{ij})$ に関する $B$ の双対基底は $(E_{ji})$ であり、$B$ は非退化である。
:ref[prop-linalg-perfect-pairing] より、$M_n(K)$ 上の任意の線形汎関数 $\varphi$ は、ただ一つの $A$ を用いて $\varphi(X)=\tr(AX)$ と書ける。
:::

### 有限次元という仮定はどこに効くか

双対基底の証明では $f=\sum_jf(e_j)e^j$ と有限和で書いた。
無限基底に対しては、ベクトルを作る一次結合は依然として有限和だが、汎関数が非零の値を取る基底ベクトルの個数には制限がない。

:::example[有限台の列とすべての列]
$V=K^{(\N)}$ を、有限個の成分だけが非零である数列の空間とする。
各座標の単位列 $e_j$ は基底である。
任意の数列 $a=(a_j)_{j\ge1}$ は

$$
f_a(x)=\sum_{j\ge1}a_jx_j
$$

によって汎関数を与える。
$x$ が有限台なので毎回有限和となり、収束の仮定は要らない。
逆に任意の $f$ は $a_j=f(e_j)$ によりこの形である。

したがって $V^\vee$ は有限台に制限しない全数列の空間 $K^\N$ と同型である。
$a=(1,1,\ldots)$ の汎関数は単位列すべてに $1$ を返すので、有限個の $e^j$ の一次結合にはならない。
$e^j$ たちは一次独立だが $V^\vee$ を生成しない。
:::

同じ例で評価写像が全射でないことも確認できる。
$K^\N$ の部分空間 $F=K^{(\N)}+K\mathbf1$（$\mathbf1$ は定数列 $1$）上で

$$
L_0(b+c\mathbf1)=c\qquad(b\in K^{(\N)},\ c\in K)
$$

と定める。
有限台の定数列は零だけなのでこの表示は一意で、$L_0$ は線形である。
選択公理（Zorn の補題）を仮定すれば、任意の空間の一次独立系は基底に延長できる。
これを用いて $F$ の基底を $K^\N$ の基底に延長し、追加した基底上で値を零と定めて $L_0$ を $L\in(K^\N)^\vee=V^{\vee\vee}$ に延長する。

もし $L=\iota_V(x)$（$x\in V$）なら、各 $e_j\in K^{(\N)}\subset K^\N$ を汎関数とみて評価すると $x_j=e_j(x)=L(e_j)=L_0(e_j)=0$ となる。
すると $x=0$ だが $L(\mathbf1)=1$ であり矛盾する。
したがってこの無限次元空間は代数的二重双対に真に含まれる。
（選択公理の下で、一般の $V$ に対しても $\iota_V$ は単射である。$v\ne0$ を基底に延長すれば $f(v)=1$ となる $f$ を作れるからである。）

## 商空間

連立方程式 $Tx=b$ の解は、:ref[prop-linalg-solution-structure] により $\Ker T$ の元を足す自由をもつ。
出力が区別できない入力を一つにまとめると、その自由度を取り除いた空間が得られる。
商空間は、部分空間方向の違いを「無視する」操作を、補空間を選ばずに実現する。

:::definition[商空間と標準射影]{#def-linalg-quotient}
$U\subset V$ を部分空間とする。
$x-y\in U$ のとき $x\sim y$ と定める。
$0\in U$ から反射律、$-U=U$ から対称律、和についての閉性から推移律が従い、$\sim$ は同値関係である。
$x$ の同値類 $[x]=x+U=\{x+u:u\in U\}$ を $U$ を法とする :term[剰余類]{en="coset" yomi="じょうよるい"}、剰余類全体の集合を :term[商空間]{en="quotient space" yomi="しょうくうかん"} $V/U$ と呼ぶ。
写像 $\pi_U:V\to V/U,\ \pi_U(x)=[x]$ を :term[標準射影]{en="canonical projection" yomi="ひょうじゅんしゃえい"} という。
:::

幾何的には、$[x]$ は $U$ を $x$ だけ平行移動したアフィン部分空間であり、$V/U$ は $U$ に平行な「平行な平面の族」全体である。
$[x]=[y]$ と $x-y\in U$ は同値である。

:::proposition[商の演算、基底、次元]{#prop-linalg-quotient-structure}
$[x]+[y]=[x+y],\ a[x]=[ax]$ により $V/U$ はベクトル空間になり、$\pi_U$ は全射線形で $\Ker\pi_U=U$ である。
$V$ が有限次元で、$U$ の基底 $u_1,\ldots,u_k$ を $V$ の基底 $u_1,\ldots,u_k,v_1,\ldots,v_m$ に延長すれば、$[v_1],\ldots,[v_m]$ は $V/U$ の基底である。
したがって

$$
\dim V/U=\dim V-\dim U .
$$

:::

:::proof
演算が代表元によらないことを確かめる。
$[x]=[x'],\ [y]=[y']$ とし $x'=x+u,\ y'=y+w\ (u,w\in U)$ と書くと、$(x'+y')-(x+y)=u+w\in U$、$ax'-ax=au\in U$ である。
よって和とスカラー倍は定まる。
ベクトル空間の公理は、代表元で成り立つ等式（結合法則など）に $\pi_U$ を施せば従う。
零元は $[0]=U$、$[x]$ の逆元は $[-x]$ である。
$\pi_U$ の線形性は演算の定義そのものであり、全射は明らかで、$[x]=[0]\iff x\in U$ より $\Ker\pi_U=U$ である。

有限次元の場合、$x=\sum a_iu_i+\sum b_jv_j$ なら $[x]=\sum b_j[v_j]$ なので $[v_j]$ が生成する。
$\sum b_j[v_j]=0$ なら $\sum b_jv_j\in U$ であり、$\sum b_jv_j=\sum a_iu_i$ と書けるので、延長した基底の独立性から $b_j=0$ である。
:::

:::example[同値類を座標で読む]
$V=\R^2,\ U=\spn\{(1,1)\}$ とする。
$(3,1),(2,0),(0,-2)$ はすべて同じ類であり、違いは $U$ 方向だけである。
各類には $(d,0)$ の形の元がただ一つあり、$d=x-y$ が商の座標になる（$V/U\to\R,\ [(x,y)]\mapsto x-y$ は同型）。
$U$ の元の類は零元一個にまとまるのであって、$U$ の元が商から単に取り除かれるわけではない。
:::

:::example[多項式の商]
$V=K[t]_{\le3},\ U=\spn\{t^2,t^3\}$ とする。
$[a+bt+ct^2+dt^3]=[a+bt]$ なので $[1],[t]$ が商の基底である。
例えば $[1+t^2]+[2t-t^2]=[1+2t]$ である。
ここで定めたのは和とスカラー倍であり、一般の部分空間を法とする商に積が自動的に入るわけではない。
（$U$ が積について「イデアル」になっている場合、例えば $K[t]$ を $t^2K[t]$ で割る場合には積も定まる。）
:::

:::theorem[商空間の普遍性]{#thm-linalg-quotient-universal}
$U\subset V$ を部分空間、$T:V\to W$ を線形写像とする。
$T=\bar T\circ\pi_U$ を満たす線形写像 $\bar T:V/U\to W$ が存在するための必要十分条件は $U\subset\Ker T$ である。
このとき $\bar T$ は一意で、$\Img\bar T=\Img T$、$\Ker\bar T=\pi_U(\Ker T)=\Ker T/U$ である。
:::

:::proof
分解できれば $u\in U$ に対して $Tu=\bar T[u]=\bar T[0]=0$ となる。
逆に $U\subset\Ker T$ なら、$[x]=[y]$ のとき $T(x-y)=0$ なので $\bar T[x]=Tx$ は代表元によらず定まる。
$\bar T([x]+[y])=\bar T[x+y]=T(x+y)=Tx+Ty$ などから線形である。
すべての類が $\pi_U(x)$ の形なので、$\bar T\pi_U=T$ という条件が $\bar T$ を一意に決める。
像は明らかに等しく、$\bar T[x]=0\iff x\in\Ker T$ より $\Ker\bar T=\{[x]:x\in\Ker T\}$ である。
:::

```tikzcd caption="商空間の普遍性：U⊂Ker T なら T は商を通って一意に分解する"
V \arrow[r, "T"] \arrow[d, "\pi_U"'] & W \\
V/U \arrow[ur, dashed, "\bar T"'] &
```

### 補空間と商

:::proposition[補空間による商の実現]{#prop-linalg-complement-quotient}
$V=U\oplus C$ なら、$\pi_U|_C:C\to V/U$ は同型である。
逆に有限次元では補空間が存在するので（:ref[cor-linalg-complement-exists]）、$V/U$ は任意の補空間と同型である。
:::

:::proof
$\pi_U|_C$ の核は $C\cap U=0$ なので単射である。
任意の $x\in V$ を $x=u+c$ と分解すれば $[x]=[c]$ なので全射である。
:::

補空間 $C$ 自体は一意でない（$\R^2$ で $U$ を横軸とすれば、横軸以外のどの直線も補空間である）。
これに対し $V/U$ は $U$ だけから定まる。
補空間を一つ選ぶことは、各類から代表元を「線形に」一つずつ選ぶことに等しい。

### 商空間と双対空間

:::proposition[商の双対と部分空間の双対]{#prop-linalg-quotient-dual}
$U\subset V$ を部分空間とし、$i:U\to V$ を包含写像とする。
1. $\pi_U^\vee:(V/U)^\vee\to V^\vee$ は単射で、その像は $U^\circ$ である。したがって $(V/U)^\vee\simeq U^\circ$ である。
2. $V$ が有限次元なら、制限 $i^\vee:V^\vee\to U^\vee,\ f\mapsto f|_U$ は全射で核は $U^\circ$ である。したがって $V^\vee/U^\circ\simeq U^\vee$ である。
:::

:::proof
1. $\pi_U$ は全射なので、$h\circ\pi_U=0$ なら $h=0$ であり、$\pi_U^\vee$ は単射である。
$h\circ\pi_U$ は $U$ を消すので像は $U^\circ$ に含まれる。
逆に $f\in U^\circ$ なら :ref[thm-linalg-quotient-universal] より $f=\bar f\circ\pi_U$ と分解でき、$f=\pi_U^\vee\bar f$ である。

2. $U$ 上の汎関数は、$U$ の基底を $V$ の基底に延長し追加した基底上で $0$ と定めれば $V$ 上に延長できるので、$i^\vee$ は全射である。
$f|_U=0$ と $f\in U^\circ$ は同値である。
したがって $[f]\mapsto f|_U$ は代表元によらず全射であり、零に移る類は $U^\circ$ 自身だけである。よって $V^\vee/U^\circ\simeq U^\vee$ である。
:::

次元を比べると、1 は $\dim V-\dim U=\dim U^\circ$、2 は $\dim V-\dim U^\circ=\dim U$ となり、:ref[prop-linalg-annihilator] と整合する。
双対を取ると「部分空間」と「商空間」の役割が入れ替わる。

### 計量があれば最短の代表元を選べる

商 $V/U$ にベクトルの長さは自動的には定まらない。
有限次元の $V$ に内積を与えたときは、各類 $[v]$ の中で長さ最小の元として $v-P_Uv$ を選べる（$P_U$ は直交射影）。
実際、$v+u=(v-P_Uv)+(P_Uv+u)$ は直交和なので $\norm{v+u}^2=\norm{v-P_Uv}^2+\norm{P_Uv+u}^2$ であり、第二項を消す $u=-P_Uv$ でだけ最小になる。
また代表元を $v+w\ (w\in U)$ に変えると $P_U(v+w)=P_Uv+w$ なので、差 $v-P_Uv$ は変わらない。

したがって

$$
[v]\longmapsto (I-P_U)v
$$

は $V/U$ から $U^\perp$ への線形同型であり、:ref[prop-linalg-complement-quotient] で $C=U^\perp$ と選んだ場合の逆写像である。
ここで用いた最良近似の存在と一意性は :ref[thm-linalg-projection] の内容である。
$\norm{[v]}\defeq\min_{u\in U}\norm{v+u}=\norm{(I-P_U)v}$ を商ノルムという。

:::example[同じ類の最短代表]
標準内積の $\R^2$、$U=\spn\{(1,1)\}$ で、$(3,1)$ の類の最短代表は $(3,1)-(2,2)=(1,-1)$ である。
長さは $\sqrt2$ で、同じ類の代表 $(2,0)$ の長さ $2$ より短い。
座標の一方を零にする代表は見つけやすいが、内積を変えれば最短になる代表も変わる。
:::

### 商空間上の誘導写像

:::proposition[不変部分空間と誘導作用素]{#prop-linalg-induced-operator}
$T\in\End(V)$ が $T(U)\subset U$ を満たすとする。
1. $\widetilde T[v]=[Tv]$ により $\widetilde T\in\End(V/U)$ が定まり、$\pi_U T=\widetilde T\pi_U$ である。
2. 任意の多項式 $p$ について $p(\widetilde T)=\widetilde{p(T)}$ である。特に $p(T)=0$ なら $p(\widetilde T)=0$ かつ $p(T|_U)=0$ である。
3. $V$ が有限次元なら、$U$ の基底 $u_1,\ldots,u_k$ を $V$ の基底 $u_1,\ldots,u_k,v_1,\ldots,v_m$ に延長すると、$T$ の表現行列は $\begin{pmatrix}B&C\\0&D\end{pmatrix}$ の形になり、$B$ は $T|_U$ の、$D$ は基底 $[v_1],\ldots,[v_m]$ に関する $\widetilde T$ の表現行列である。
4. 有限次元では $\Phi_T=\Phi_{T|_U}\Phi_{\widetilde T}$、$\det T=\det(T|_U)\det\widetilde T$ である。
:::

:::proof
1. $\pi_UT:V\to V/U$ は $U$ を消す（$Tu\in U$）ので、:ref[thm-linalg-quotient-universal] により $V/U$ を通って分解する。
2. 不変部分空間を保つ写像について $\widetilde I=I$ であり、$\widetilde{ST}=\widetilde S\widetilde T$、$\widetilde{S+T}=\widetilde S+\widetilde T$ は代表元で計算すれば明らかであり、これを :ref[lem-linalg-polynomial-substitution] の多項式表示の各項に適用すればよい。
$p(T)=0$ なら $p(\widetilde T)=\widetilde0=0$ であり、$p(T|_U)=p(T)|_U=0$ である。
3. $Tu_j\in U$ なので最初の $k$ 列の下側は零である。
$Tv_j=\sum_ic_{ij}u_i+\sum_ld_{lj}v_l$ に $\pi_U$ を施すと $\widetilde T[v_j]=\sum_ld_{lj}[v_l]$ となり、$D=(d_{lj})$ が $\widetilde T$ の行列である。
4. まず :ref[prop-linalg-invariant-block] のブロック表示に :ref[prop-linalg-det-block-triangular] を適用すると $\det T=\det B\det D$ である。同じ基底で $tI-T$ は対角ブロックが $tI-B,tI-D$ の上三角ブロック行列になる。:ref[prop-linalg-det-ring] により $K[t]$ でもこの計算ができるので、:ref[cor-linalg-invariant-charpoly] の第二因子は商上の特性多項式 $\Phi_{\widetilde T}$ である。
:::

有限次元では 2 と :ref[prop-linalg-minimal-polynomial-divides] から、$T|_U$ と $\widetilde T$ の最小多項式はともに $T$ の最小多項式を割り切る。
特性多項式には非対角ブロック $C$ が現れないが、$C$ を無視して $V$ 上の作用まで直和とみなすことはできない。

:::example[制限と商だけでは得られない情報]
$T=\begin{pmatrix}1&1\\0&1\end{pmatrix}$、$U=\spn\{e_1\}$ とする。
制限 $T|_U$ も誘導作用素 $\widetilde T$ も一次元の恒等写像である。
しかし $(T-I)e_2=e_1$ なので $T\ne I$ であり、$U$ の $T$ 不変な補空間 $L$ も存在しない。
もし存在すれば、$\pi_U|_L:L\to V/U$ は $T|_L$ と $\widetilde T=\id$ を結ぶ同型なので $T|_L=\id$ となり、$V=U\oplus L$ 全体で $T=I$ となってしまうからである。
このずれを鎖の長さで記述するのが前章の Jordan 標準形（:ref[thm-linalg-jordan]）である。
:::

:::example[商による三角化の別証明]
商空間を使うと、:ref[thm-linalg-triangularization] の十分性を次のように示せる。
$\dim V=n$ とし、$\Phi_T$ が一次式の積に分解するとする。
$n=1$ なら明らかなので $n\ge2$ とし、$n-1$ 次元以下で主張が成り立つとする。
$\Phi_T$ の根 $\lambda$ は固有値なので固有ベクトル $v_1$ があり、$U=\spn\{v_1\}$ は不変である。
:ref[prop-linalg-induced-operator] の 4 より $\Phi_T=(t-\lambda)\Phi_{\widetilde T}$ なので、$\Phi_{\widetilde T}$ も一次式の積に分解する。
帰納法の仮定から、$V/U$ の基底 $[v_2],\ldots,[v_n]$ で $\widetilde T[v_j]\in\spn\{[v_2],\ldots,[v_j]\}$ となるものがある。
代表 $v_2,\ldots,v_n$ を選ぶと $v_1,\ldots,v_n$ は $V$ の基底である。
実際 $\sum_ic_iv_i=0$ に $\pi_U$ を施すと $\sum_{i\ge2}c_i[v_i]=0$ より $c_2=\cdots=c_n=0$、残る $c_1v_1=0$ より $c_1=0$ であり、$n$ 本の独立なベクトルは基底である（:ref[cor-linalg-n-vectors-basis]）。
さらに $Tv_j\in\spn\{v_1,\ldots,v_j\}$ が成り立つ。
実際 $[Tv_j]\in\spn\{[v_2],\ldots,[v_j]\}$ は $Tv_j-\sum_{i=2}^jc_iv_i\in U$ を意味する。
したがってこの基底で $T$ は上三角である。
:::

## 同型定理

商の普遍性が保証する写像の核と像を調べれば、異なる手順で構成した空間が同型になることを確かめられる。
まず写像が失った情報だけを入力から取り除く。

:::theorem[第一同型定理]{#thm-linalg-first-isomorphism}
線形写像 $T:V\to W$ に対し、

$$
\bar T:V/\Ker T\to\Img T,\qquad [v]\mapsto Tv
$$

は矛盾なく定まる線形同型であり、$T=\iota\circ\bar T\circ\pi_{\Ker T}$（$\iota:\Img T\to W$ は包含）と分解する。
:::

:::proof
:ref[thm-linalg-quotient-universal] を $U=\Ker T$ に適用すると、$\bar T$ は代表元によらず線形で、$\Ker\bar T=\Ker T/\Ker T=0$、$\Img\bar T=\Img T$ である。
よって値域を $\Img T$ とすれば全単射である。
:::

任意の線形写像は「全射 $\pi$」「同型 $\bar T$」「単射 $\iota$」の合成に分解される。
$V$ が有限次元なら :ref[prop-linalg-quotient-structure] と合わせて $\dim V-\dim\Ker T=\dim\Img T$ を得る。
これは次元定理 :ref[thm-linalg-rank-nullity] の、基底の延長を商に押し込めた証明である。

:::example[微分で失う自由度]
$D:\R[t]_{\le3}\to\R[t]_{\le2},\ p\mapsto p'$ の核は定数であり、二次以下の多項式には三次以下の原始関数があるので $D$ は全射である。
第一同型定理より、定数を法とする商から $\R[t]_{\le2}$ への同型 $[p]\mapsto p'$ を得る。
例えば $3t^2-2$ に対応する類は $[t^3-2t]$ であり、「不定積分は定数の差を除いて定まる」という言明の正確な形がこれである。
:::

:::theorem[第二・第三同型定理]{#thm-linalg-quotient-isomorphisms}
1. 部分空間 $U,W\subset V$ に対し、$u+(U\cap W)\mapsto u+W$ は同型 $U/(U\cap W)\simeq(U+W)/W$ を与える。
2. $U\subset W\subset V$ なら、$W/U$ は $V/U$ の部分空間であり、$(v+U)+W/U\mapsto v+W$ は同型 $(V/U)/(W/U)\simeq V/W$ を与える。
:::

:::proof
1. $\phi:U\to(U+W)/W,\ u\mapsto u+W$ は線形である。
任意の元 $(u+w)+W$ は $u+W=\phi(u)$ に等しいので全射で、$\phi(u)=0\iff u\in W$ より核は $U\cap W$ である。
第一同型定理を適用すればよい。

2. $W/U=\pi_U(W)$ は部分空間の像なので部分空間である。
$\psi:V/U\to V/W,\ v+U\mapsto v+W$ は、$U\subset W$ より代表元の差が $U$ にあれば $W$ にもあるので定義でき（:ref[thm-linalg-quotient-universal] を $\pi_W$ に適用）、線形かつ全射である。
$\psi(v+U)=0\iff v\in W$ なので核は $W/U$ である。
再び第一同型定理を適用する。
:::

:::corollary[次元公式の別証明]
有限次元の部分空間 $U,W$ について $\dim(U+W)+\dim(U\cap W)=\dim U+\dim W$ である。
:::

:::proof
第二同型定理の両辺の次元を :ref[prop-linalg-quotient-structure] で計算すると $\dim U-\dim(U\cap W)=\dim(U+W)-\dim W$ である。
これは :ref[thm-linalg-dimension-formula] にほかならない。
:::

:::example[段階的に方向を消す]
$V=\R^3,\ U=\spn\{e_1\},\ W=\spn\{e_1,e_2\}$ とする。
最初に $U$ を法とすると座標 $(y,z)$ が残り、次に $W/U$（$y$ 軸の像）を法とすると $z$ が残る。
第三同型定理の同型は、こうして得た類を、初めから $W$ を法とした類 $(0,0,z)+W$ へ送る。
:::

:::example[第二同型定理で像に残る方向を特定する]
$\R^3$ で $U=\spn\{e_1,e_2\},\ W=\spn\{e_2,e_1+e_3\}$ とする。
$W$ の元 $ae_2+b(e_1+e_3)$ が $U$ に入るのは $e_3$ 成分が消える $b=0$ のときなので $U\cap W=\spn\{e_2\}$ である。
一方 $U+W$ は $e_1,e_2$ と $(e_1+e_3)-e_1=e_3$ を含み、$\R^3$ 全体になる。

$U/(U\cap W)$ は $e_1+(U\cap W)$ を基底にもつ一次元空間である。
$(U+W)/W$ では $e_2+W=0$、$(e_1+e_3)+W=0$ なので $e_1+W=-e_3+W$ であり、これは非零の類である。
第二同型定理の写像は $e_1+(U\cap W)$ をこの類に送る。
例えば $2e_1+3e_2$ の右側の類は $-2e_3$ でも代表できる。
$-2e_3\notin U$ であるが、$2e_1+3e_2-(-2e_3)=2(e_1+e_3)+3e_2\in W$ なので商では同じ元である。
:::

### 部分空間の対応

:::proposition[対応定理]{#prop-linalg-correspondence}
$U\subset V$ を部分空間とする。
$U$ を含む $V$ の部分空間 $L$ と、$V/U$ の部分空間 $M$ は

$$
L\longmapsto L/U=\pi_U(L),\qquad M\longmapsto\pi_U^{-1}(M)
$$

により一対一に対応し、この対応は包含関係と和・共通部分を保つ。
有限次元では $\dim L/U=\dim L-\dim U$ である。
:::

:::proof
$\pi_U^{-1}(M)$ は部分空間の逆像なので部分空間で、$U=\pi_U^{-1}(0)$ を含む。
$\pi_U$ は全射なので $\pi_U(\pi_U^{-1}(M))=M$ である。
他方向では、$v\in\pi_U^{-1}(\pi_U(L))\iff v+U=l+U\ (\exists l\in L)\iff v\in L+U=L$ なので $\pi_U^{-1}(\pi_U(L))=L$ である。
包含を保つことは両写像の定義から明らかで、順写像と逆写像の両方が包含を保つため、最小上界（和）と最大下界（共通部分）も対応する。
次元は :ref[prop-linalg-quotient-structure] を $L$ に適用すればよい。
:::

### 余核と指数

:::definition[余核]{#def-linalg-cokernel}
線形写像 $T:V\to W$ に対し $\operatorname{Coker}T\defeq W/\Img T$ を $T$ の :term[余核]{en="cokernel" yomi="よかく"} という。
:::

核が「単射からのずれ」を測るのに対し、余核は「全射からのずれ」を測る。
$\operatorname{Coker}T=0$ と $T$ の全射性は同値である。
有限次元では $\dim\operatorname{Coker}T=\dim W-\rank T$ であり、:ref[prop-linalg-quotient-dual] の 1 より

$$
(\operatorname{Coker}T)^\vee\simeq(\Img T)^\circ=\Ker T^\vee
$$

である。
すなわち :ref[cor-linalg-fredholm-alternative] は「余核の双対は双対写像の核」と言い換えられる。

:::proposition[指数公式]{#prop-linalg-index}
有限次元空間の間の線形写像 $T:V\to W$ について

$$
\dim\Ker T-\dim\operatorname{Coker}T=\dim V-\dim W .
$$

特に左辺（$T$ の指数）は $T$ によらず、$V,W$ の次元だけで決まる。
:::

:::proof
次元定理より $\dim\Ker T=\dim V-\rank T$、定義より $\dim\operatorname{Coker}T=\dim W-\rank T$ である。
差を取ればよい。
:::

$\dim V=\dim W$ なら指数は $0$ で、「単射 $\iff$ 全射」（:ref[cor-linalg-equal-dimension-map]）は指数公式の特別な場合である。
無限次元では両辺が意味をもたないことがあり、Banach 空間の間の有界線形写像では、像が閉じており、核と余核がともに有限次元であるものを Fredholm 作用素という。この場合は左辺だけで指数を定義できる。

### 短完全列

:::definition[完全列]{#def-linalg-exact-sequence}
ベクトル空間と線形写像の列 $\cdots\to V_{i-1}\xrightarrow{T_{i-1}}V_i\xrightarrow{T_i}V_{i+1}\to\cdots$ が各 $V_i$ で $\Img T_{i-1}=\Ker T_i$ を満たすとき、:term[完全列]{en="exact sequence" yomi="かんぜんれつ"} という。
特に

$$
0\longrightarrow U\xrightarrow{\ i\ }V\xrightarrow{\ p\ }W\longrightarrow0
$$

が完全であること、すなわち $i$ が単射、$p$ が全射、$\Img i=\Ker p$ であることを :term[短完全列]{en="short exact sequence" yomi="たんかんぜんれつ"} という。
:::

部分空間 $U\subset V$ に対し、包含 $U\to V$ と射影 $V\to V/U$ は短完全列をなす。
逆に任意の短完全列では、第一同型定理より $p$ は同型 $V/i(U)\simeq W$ を誘導するので、すべての短完全列は本質的にこの形である。

:::proposition[短完全列の分裂]{#prop-linalg-split-exact}
$0\to U\xrightarrow{i}V\xrightarrow{p}W\to0$ を短完全列とし、$W$ は有限次元とする。
このとき $p\circ s=\id_W$ を満たす線形写像 $s:W\to V$（切断）が存在し、任意の切断について $V=i(U)\oplus s(W)$ である。
特に $V$ が有限次元なら $\dim V=\dim U+\dim W$ である。
:::

:::proof
$W$ の基底 $w_1,\ldots,w_m$ の各々について $p$ による原像 $v_j$ を一つ選び（$p$ は全射）、$s(w_j)=v_j$ として線形に延長すれば $ps=\id_W$ である。
任意の $v$ について $v=(v-sp(v))+sp(v)$ と書くと、$p(v-sp(v))=p(v)-p(v)=0$ より第一項は $\Ker p=i(U)$ に、第二項は $s(W)$ に入る。
$x\in i(U)\cap s(W)$ を $x=s(w)$ と書くと $0=p(x)=ps(w)=w$ より $x=0$ である。
次元の式は $i,s$ が単射であることから従う。
:::

一つの切断 $s_0$ を固定する。別の切断 $s$ との差は $p(s-s_0)=0$ を満たすので $i(U)$ に値をとる。$i:U\to i(U)$ が同型であることから、$s=s_0+iL$ を満たす線形写像 $L:W\to U$ が一意にある。逆にどの $L$ に対しても $p(s_0+iL)=\id_W$ である。これで補空間を選ぶ自由が具体的にわかる。特に $U$ または $W$ が零なら切断は一意である。
一方、商 $V/i(U)$ と $W$ の同型は $p$ だけで決まる。

:::example[切断の違いを一つのパラメータで表す]
$p:\R^2\to\R,\ p(x,y)=y$ とする。
線形な切断は $s(1)=(a,1)$ で決まり、$s_a(t)=(at,t)$ と書ける（$ps(1)=1$ が必要十分）。
$\Ker p$ はつねに横軸で、補空間 $s_a(\R)$ の傾きだけが選択に依存する。
:::

:::proposition[双対列の完全性]{#prop-linalg-dual-exact}
有限次元空間の短完全列 $0\to U\xrightarrow{i}V\xrightarrow{p}W\to0$ の双対

$$
0\longrightarrow W^\vee\xrightarrow{\ p^\vee\ }V^\vee\xrightarrow{\ i^\vee\ }U^\vee\longrightarrow0
$$

も短完全列である。
:::

:::proof
:ref[thm-linalg-dual-kernel] により、$p$ が全射なので $p^\vee$ は単射、$i$ が単射なので $i^\vee$ は全射である。
中央では $\Img p^\vee=(\Ker p)^\circ=(\Img i)^\circ=\Ker i^\vee$ である。
:::

$i$ を包含 $U\to V$、$p$ を $\pi_U$ とすれば、この命題は :ref[prop-linalg-quotient-dual] の二つの同型を一つにまとめたものである。

:::proposition[完全列の次元の交代和]{#prop-linalg-euler-characteristic}
有限次元空間の完全列 $0\to V_1\xrightarrow{T_1}V_2\xrightarrow{T_2}\cdots\xrightarrow{T_{k-1}}V_k\to0$ について

$$
\sum_{i=1}^k(-1)^i\dim V_i=0 .
$$

:::

:::proof
各 $i$ について次元定理より $\dim V_i=\dim\Ker T_i+\dim\Img T_i$ である（$T_k=0$ とする）。
完全性から $\dim\Ker T_i=\dim\Img T_{i-1}$（$\Img T_0=0$ とする）なので、$r_i=\dim\Img T_i$ とおくと $\dim V_i=r_{i-1}+r_i$、$r_0=r_k=0$ である。
交代和を取ると $\sum_i(-1)^i(r_{i-1}+r_i)$ は隣り合う項が打ち消し合い、$-r_0+(-1)^kr_k=0$ となる。
:::

指数公式 :ref[prop-linalg-index] は、完全列 $0\to\Ker T\to V\xrightarrow{T}W\to\operatorname{Coker}T\to0$ にこの命題を適用したものである。

### 二つの部分空間を法として写像を移す

:::proposition[商から商への誘導写像]{#prop-linalg-induced-quotient-map}
$A:V\to W$ を線形写像とし、$U\subset V,\ Z\subset W$ を部分空間とする。
$\widetilde A[x]_U=[Ax]_Z$ により線形写像 $\widetilde A:V/U\to W/Z$ が定まるための必要十分条件は $A(U)\subset Z$ である。
このとき

$$
\Ker\widetilde A=A^{-1}(Z)/U,\qquad
\Img\widetilde A=(\Img A+Z)/Z
$$

であり、$\widetilde A$ が単射 $\iff A^{-1}(Z)=U$、全射 $\iff\Img A+Z=W$ である。
:::

:::proof
$\pi_Z\circ A:V\to W/Z$ に :ref[thm-linalg-quotient-universal] を適用すると、分解できる条件は $U\subset\Ker(\pi_ZA)=A^{-1}(Z)$、すなわち $A(U)\subset Z$ である。
定理の後半から $\Ker\widetilde A=\Ker(\pi_ZA)/U=A^{-1}(Z)/U$、$\Img\widetilde A=\pi_Z(\Img A)=(\Img A+Z)/Z$ である。
単射・全射の条件は :ref[prop-linalg-correspondence] により従う。
:::

:::example[入力と出力で異なる方向を無視する]
$A:\R^3\to\R^3,\ A(x,y,z)=(x+y,z,0)$、$U=Z=\spn\{e_1\}$ とする。
$A(U)\subset Z$ なので商には写像が定まる。
入力の商座標を $(y,z)$、出力の商座標を第二、第三成分の順に取ると、$\widetilde A(y,z)=(z,0)$ である。

元の座標では $A^{-1}(Z)$ は $z=0$ の平面であり、$U$ を法とすると $y$ 方向だけが残る。
これは確かに $\Ker\widetilde A$ である。
また $\Img A$ は第一、第二座標平面なので、$Z$ を法とした像は第二座標の直線となる。
商を先に作っても、元の空間で逆像と像を計算してから商を取っても同じ結果に至る。
:::

:::example[同じ写像を誘導する条件]
$A,B:V\to W$ がともに $A(U),B(U)\subset Z$ を満たすとき、同じ写像 $V/U\to W/Z$ を誘導するための必要十分条件は $\Img(A-B)\subset Z$ である。
実際、$[Av]_Z=[Bv]_Z$ がすべての $v$ で成り立つことは $(A-B)v\in Z$ と同値である。
例えば $V=W=\R^2,\ U=Z=\spn\{e_1\}$ で

$$
A=\begin{pmatrix}1&2\\0&3\end{pmatrix},\qquad B=\begin{pmatrix}4&-1\\0&3\end{pmatrix}
$$

はどちらも $U$ を $Z$ に移し、$A-B$ の第二行が零なので差の像は $Z$ に入る。
商で残る第二座標への作用はともに $y\mapsto3y$ であり、元の行列が異なっていても商上の作用は区別できない。
:::

## 演習問題

::::exercise{level=basic}
$\R^2$ の基底 $b_1=(2,1),\ b_2=(1,1)$ の双対基底を求めよ。

:::solution
$b_j$ を列とする行列 $S=\begin{pmatrix}2&1\\1&1\end{pmatrix}$ の行列式は $1$ で、$S^{-1}=\begin{pmatrix}1&-1\\-1&2\end{pmatrix}$ である。
:ref[prop-linalg-dual-basis-change] より双対基底は $S^{-1}$ の行で与えられ、$b^1(x,y)=x-y,\ b^2(x,y)=-x+2y$ である。
検算：$b^1(b_1)=1,\ b^1(b_2)=0,\ b^2(b_1)=0,\ b^2(b_2)=1$。
:::
::::

::::exercise{level=basic}
$T:\R^3\to\R^2,\ T(x,y,z)=(x+z,y-z)$ に対し、$g(u,v)=2u+3v$ の双対写像による像を求めよ。

:::solution
定義から

$$
(T^\vee g)(x,y,z)=g(T(x,y,z))=2(x+z)+3(y-z)=2x+3y-z
$$

である。
係数列 $(2,3,-1)$ は、$T$ の行列 $\begin{pmatrix}1&0&1\\0&1&-1\end{pmatrix}$ の転置を $\transp{(2,3)}$ に掛けたものに等しい（:ref[prop-linalg-dual-transpose]）。
:::
::::

::::exercise{level=basic}
$\R^3$ の $U=\spn\{(1,2,3)\}$ の零化空間の基底を求め、$U$ を二本の一次方程式の共通解として表せ。

:::solution
$(a,b,c)$ が $U$ を消す条件は $a+2b+3c=0$ であり、解は $b(-2,1,0)+c(-3,0,1)$ である。
$b,c$ は自由で表示は一意なので、$(-2,1,0),(-3,0,1)$ が $U^\circ$ の基底である（$\dim U^\circ=3-1=2$ と整合する）。
$(U^\circ)_\circ=U$ より、$U$ は $-2x+y=0,\ -3x+z=0$ の共通解であり、実際 $(1,2,3)$ の定数倍全体である。
:::
::::

::::exercise{level=basic}
$V=\R[t]_{\le3}$、$U=\{p\in V:p(1)=p(-1)=0\}$ のとき $\dim U$ と $V/U$ の基底を求めよ。

:::solution
$\mathrm{ev}_1,\mathrm{ev}_{-1}$ は一次独立な汎関数（一次関係の係数を $a,b$ とすれば、$1$ で $a+b=0$、$t$ で $a-b=0$ となり、実数上では $a=b=0$）なので、:ref[prop-linalg-annihilator] より $\dim U=4-2=2$ である。
$p\in V$ を $t^2-1$ で割ると（:ref[lem-linalg-polynomial-division]）$p=(t^2-1)q+a+bt$ と書け、$(t^2-1)q\in U$ なので $[p]=[a+bt]$ である。
よって $[1],[t]$ が $V/U$ を生成し、$\dim V/U=2$ なので基底である。
写像 $[p]\mapsto(p(1),p(-1))$ が $V/U\simeq\R^2$ を与えることも第一同型定理からわかる。
:::
::::

::::exercise{level=basic}
$M_n(K)$ 上の線形汎関数はすべて、ただ一つの $A\in M_n(K)$ により $X\mapsto\tr(AX)$ と書けることを、双対基底を用いて直接示せ。
さらに $n\ge2$ として $\varphi(X)=\tr X$ と $\psi(X)=x_{12}$ に対応する $A$ を求めよ。

:::solution
行列単位 $E_{ij}$ は $M_n(K)$ の基底であり、その双対基底 $E^{ij}$ は $X$ の $(i,j)$ 成分 $x_{ij}$ を取り出す汎関数である。
$\tr(E_{ji}X)=\sum_k(E_{ji}X)_{kk}=x_{ij}$ なので $E^{ij}(X)=\tr(E_{ji}X)$ である。
任意の $\varphi$ は :eqref[eq-linalg-dual-expansion] より $\varphi=\sum_{i,j}\varphi(E_{ij})E^{ij}$ なので

$$
\varphi(X)=\tr(AX),\qquad A=\sum_{i,j}\varphi(E_{ij})E_{ji}
$$

と書ける。
一意性は、$\tr(AX)=0$ が全 $X$ で成り立てば $X=E_{ji}$ として $a_{ij}=0$ となることによる。
$\varphi=\tr$ では $\varphi(E_{ij})=\delta_{ij}$ なので $A=I$、$\psi$ では $\psi(E_{ij})=\delta_{i1}\delta_{j2}$ なので $A=E_{21}$ である。
:::
::::

::::exercise{level=basic}
$T(x,y)=(y,x)$ は $\R^2$ の $U=\spn\{(1,0)\}$ を法とする商に作用素を誘導するか。
代表元でも確かめよ。

:::solution
$T(1,0)=(0,1)\notin U$ なので、:ref[prop-linalg-induced-operator] の不変性の条件を満たさない。
実際 $(0,0)$ と $(1,0)$ は同じ類だが、像 $(0,0)$ と $(0,1)$ は異なる類であり、$[v]\mapsto[Tv]$ は代表元によって値が変わる。
:::
::::

::::exercise{level=advanced}
$V$ を任意のベクトル空間、$f_1,\ldots,f_k\in V^\vee$ とする。
$f_1,\ldots,f_k$ が一次独立であるための必要十分条件は、$f_i(v_j)=\delta_{ij}$ を満たす $v_1,\ldots,v_k\in V$ が存在することであることを示せ。

:::hint
$\Phi(v)=(f_1(v),\ldots,f_k(v))$ の全射性を調べよ。
:::

:::solution
$v_j$ が存在すれば、$\sum_ic_if_i=0$ を $v_j$ に適用して $c_j=0$ となるので独立である。
逆に独立とし、$\Phi:V\to K^k$ を上のように定める。
$\Img\Phi\ne K^k$ なら、:ref[prop-linalg-annihilator] を $K^k$ に適用して $\Img\Phi$ を消す非零の汎関数 $y\mapsto\sum_ic_iy_i$ が得られ、$\sum_ic_if_i=0$ となって独立性に反する。
よって $\Phi$ は全射であり、$\Phi(v_j)=e_j$ を満たす $v_j$ を選べばよい。
（$V$ が無限次元でも $K^k$ は有限次元なので議論はそのまま通る。）
:::
::::

::::exercise{level=advanced}
有限次元空間 $V$ の部分空間 $U,W$ に対し、$\dim(U\cap W)^\circ$ と $\dim(U+W)^\circ$ を $\dim V,\dim U,\dim W,\dim(U\cap W)$ で表せ。
また $U^\circ\cap W^\circ=0$ となる条件を述べよ。

:::solution
:ref[prop-linalg-annihilator] と :ref[thm-linalg-dimension-formula] より

$$
\dim(U\cap W)^\circ=\dim V-\dim(U\cap W),\qquad
\dim(U+W)^\circ=\dim V-\dim U-\dim W+\dim(U\cap W)
$$

である。
:ref[prop-linalg-annihilator-lattice] より $U^\circ\cap W^\circ=(U+W)^\circ$ なので、これが $0$ となるのは $U+W=V$ のときに限る。
:::
::::

::::exercise{level=advanced}
有限次元の $V$ と部分空間 $U$ について、包含 $i:U\to V$ と射影 $\pi:V\to V/U$ に対し $\Img\pi^\vee=\Ker i^\vee$ を直接示し、$i^\vee$ が全射であることを説明せよ。

:::solution
$i^\vee f=f\circ i=f|_U$ なので $\Ker i^\vee=U^\circ$ である。
$\pi^\vee h=h\circ\pi$ は $U$ を消すので $\Img\pi^\vee\subset U^\circ$ であり、逆に $f\in U^\circ$ は商の普遍性により $f=\bar f\circ\pi=\pi^\vee\bar f$ と書けるので $U^\circ\subset\Img\pi^\vee$ である。
$U$ 上の汎関数 $g$ は、$U$ の基底を $V$ の基底に延長し、追加した基底上で値を $0$ とすれば $V$ 上の汎関数 $f$ に延長でき、$i^\vee f=g$ となるので $i^\vee$ は全射である。
:::
::::

::::exercise{level=advanced}
$T\in\End(V)$、$\dim V<\infty$、$T(U)\subset U$ とする。
$T$ が可逆であるための必要十分条件は、$T|_U$ と誘導作用素 $\widetilde T\in\End(V/U)$ がともに可逆であることを示せ。
可逆なとき $(T^{-1})(U)\subset U$ も示せ。

:::solution
:ref[prop-linalg-induced-operator] の 4 より $\det T=\det(T|_U)\det\widetilde T$ なので、$\det T\ne0$ と両因子が非零であることは同値である（:ref[prop-linalg-operator-det]）。
可逆なとき $T|_U:U\to U$ は単射なので全単射であり、任意の $u\in U$ は $u=Tu'$（$u'\in U$）と書ける。
よって $T^{-1}u=u'\in U$ である。
:::
::::

::::exercise{level=advanced}
有限次元空間の線形写像 $T:V\to W$、$S:W\to Z$ について

$$
\rank(ST)=\rank T-\dim(\Img T\cap\Ker S)
$$

を示し、Sylvester の階数不等式 $\rank T+\rank S\le\rank(ST)+\dim W$ を導け。

:::solution
$\Img(ST)=S(\Img T)$ である。
制限 $S|_{\Img T}:\Img T\to Z$ の核は $\Img T\cap\Ker S$、像は $S(\Img T)$ なので、次元定理（あるいは第一同型定理）より
$\rank(ST)=\dim\Img T-\dim(\Img T\cap\Ker S)$ である。
$\dim(\Img T\cap\Ker S)\le\dim\Ker S=\dim W-\rank S$ を代入すると $\rank(ST)\ge\rank T+\rank S-\dim W$ を得る。
等号は $\Ker S\subset\Img T$ のときに成り立つ。
:::
::::

::::exercise{level=advanced}
$\dim V=n$ とし、$V$ 上の双線形形式 $B:V\times V\to K$ が非退化であるとする。
部分空間 $U$ に対し $U^{\perp_B}=\{w\in V:B(u,w)=0\ (\forall u\in U)\}$ とおくとき、$\dim U^{\perp_B}=n-\dim U$ を示せ。
さらに $V=\R^2$、$B(v,w)=v_1w_2-v_2w_1$ について、任意の一次元部分空間 $U$ で $U^{\perp_B}=U$ となることを確かめよ。

:::solution
$B_R:V\to V^\vee,\ w\mapsto B(\cdot,w)$ は :ref[prop-linalg-perfect-pairing] より同型である。
$w\in U^{\perp_B}\iff B_R(w)\in U^\circ$ なので $U^{\perp_B}=B_R^{-1}(U^\circ)$ であり、$\dim U^{\perp_B}=\dim U^\circ=n-\dim U$ である。
例の $B$ は行列 $\begin{pmatrix}0&1\\-1&0\end{pmatrix}$（正則）をもつので非退化である。
$U=\spn\{u\}$ なら $B(u,u)=u_1u_2-u_2u_1=0$ より $U\subset U^{\perp_B}$ で、次元がともに $1$ なので等しい。
内積と異なり、非退化な双線形形式では $U\cap U^{\perp_B}=0$ とは限らない。
:::
::::

## 確認テスト

::::quiz
$T:V\to W$ から得る双対写像 $T^\vee$ の向きと作用を選べ。

- [ ] $V^\vee\to W^\vee,\ f\mapsto f\circ T$
- [x] $W^\vee\to V^\vee,\ g\mapsto g\circ T$
- [ ] $W\to V,\ w\mapsto T^{-1}w$

:::explanation
$g\in W^\vee$ は $W$ 上で値を測る。
先に $T$ を施してから $g$ を適用すれば $V$ 上の汎関数 $g\circ T$ となるので、$T^\vee:W^\vee\to V^\vee$ である。
$T$ の可逆性は必要なく、$f\circ T$ は $f\in V^\vee$ に対しては定義されない（合成の順序が合わない）。
:::
::::

::::quiz{answer="3" tol="1e-6"}
$\dim V=7,\ \dim U=4$ のとき $\dim(V/U)^\vee$ はいくつか。

:::explanation
$\dim V/U=7-4=3$ であり、有限次元空間とその双対の次元は等しいので $3$ である。
$(V/U)^\vee\simeq U^\circ$（:ref[prop-linalg-quotient-dual]）と $\dim U^\circ=7-4$ からも同じ答えが出る。
:::
::::

::::quiz
線形写像 $T:V\to W$ が $T=\bar T\circ\pi_U$ と商 $V/U$ を通って分解する条件を選べ。

- [x] $U\subset\Ker T$
- [ ] $\Ker T\subset U$
- [ ] $T$ が全射である

:::explanation
代表元の差 $u\in U$ が出力に影響しないこと、すなわちすべての $u\in U$ で $Tu=0$ が必要十分である（:ref[thm-linalg-quotient-universal]）。
逆の包含 $\Ker T\subset U$ だけでは分解は保証されない。例えば $T=\id_V$ と $U=V\ne0$ なら、$T$ は全射でこの包含も満たすが、零空間 $V/U$ を通る写像は零写像だけである。分解が存在したうえで $\bar T$ が単射になるのは $U=\Ker T$ のときである。
:::
::::

::::quiz{answer="7" tol="1e-6"}
$V=\R^2,\ v=(2,1),\ f(x,y)=3x+y$ とするとき $\iota_V(v)(f)$ を求めよ。

:::explanation
評価写像の定義から $\iota_V(v)(f)=f(v)=3\cdot2+1=7$ である。
:::
::::

::::quiz
実有限次元空間 $V$（$\dim V\ge1$）とその部分空間 $U$ について、基底の選び方によらず定まる線形同型をすべて選べ。

- [ ] $V\to V^\vee,\ e_i\mapsto e^i$
- [x] $V\to V^{\vee\vee},\ v\mapsto(f\mapsto f(v))$
- [x] $(V/U)^\vee\to U^\circ,\ h\mapsto h\circ\pi_U$

:::explanation
$e_i\mapsto e^i$ は基底に依存する（本文の一次元の例で、基底 $e$ と $2e$ で行き先が異なる）。
評価写像 $\iota_V$ は定義に基底を使わず、有限次元では同型である（:ref[thm-linalg-double-dual]）。
$\pi_U^\vee$ も定義に選択を含まず、:ref[prop-linalg-quotient-dual] より $U^\circ$ への同型である。
:::
::::

::::quiz{answer="3" tol="1e-6"}
$A=\begin{pmatrix}1&2&3\\2&4&6\end{pmatrix}$ を $\R^3\to\R^2$ の写像とみるとき、$\dim\Ker A+\dim\operatorname{Coker}A$ を求めよ。

:::explanation
第二行は第一行の $2$ 倍なので $\rank A=1$ である。
次元定理より $\dim\Ker A=3-1=2$、余核の定義より $\dim\operatorname{Coker}A=2-1=1$ なので、和は $3$ である。
指数公式 :ref[prop-linalg-index] による検算：$\dim\Ker A-\dim\operatorname{Coker}A=2-1=1=3-2$ である。
余核が一次元であることは、$Ax=b$ の可解条件が一本（$2b_1-b_2=0$）であることに対応する。
:::
::::

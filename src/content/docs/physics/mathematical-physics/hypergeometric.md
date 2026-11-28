---
title: "超幾何関数とその他の特殊関数"
description: "超幾何方程式の特異点と接続から、合流型の関数、Airy 関数、楕円関数、Mathieu 関数を構成する。"
notation:
  - symbol: '(a)_n'
    meaning: "上昇階乗"
  - symbol: '{}_2F_1(a,b;c;z)'
    meaning: "Gauss の超幾何関数"
  - symbol: 'M(a,b,z),U(a,b,z)'
    meaning: "Kummer の関数"
  - symbol: 'K(k)'
    meaning: "法 k を用いた第 1 種完全楕円積分"
---

## この章の目標

- 微分方程式の漸化式から超幾何級数を導く。
- 特異点を移動した局所解を接続する。
- 合流型の方程式から物理で使う特殊関数を構成する。
- 非線形振子と周期係数の方程式で適切な関数を選ぶ。

複素平面での局所解の接続については [複素領域の常微分方程式](/math/ode/complex-domain/) を参照する。正の実引数の冪は実対数を用い、複素引数では連続な枝を指定する。

## Gauss の超幾何関数と超幾何級数

### 係数比で関数を定義する

$(a)_0=1$、$(a)_n=a(a+1)\cdots(a+n-1)$ とする。この記号を :term[ポッホハマー記号]{en="Pochhammer symbol" yomi="ぽっほはまーきごう"} という。$c$ が非正整数でないとき
$$
{}_2F_1(a,b;c;z)=\sum_{n\ge0}\frac{(a)_n(b)_n}{(c)_n n!}z^n
\label{eq-mp-hypergeometric-series}
$$
を :term[超幾何関数]{en="hypergeometric function" yomi="ちょうきかかんすう"} とする。$a$ または $b$ が非正整数なら有限、その他は係数比から $|z|<1$ で収束する。$z(1-z)y''+[c-(a+b+1)z]y'-aby=0$ に代入すると $(n+1)(n+c)d_{n+1}=(n+a)(n+b)d_n$ であり、$d_0=1$ で上の級数になる。有限特異点の候補 $0,1$ での決定根はそれぞれ $(0,1-c),(0,c-a-b)$ である。無限遠では $y\sim z^{-r}$ として最高次から $-(r-a)(r-b)=0$ を得る。例えば ${}_2F_1(1,1;2;z)=-\log(1-z)/z$ は係数 $1/(n+1)$ を積分することで従う。

### 分子が負整数の級数

Gauss 級数の項の比は $(a+n)(b+n)/((c+n)(n+1))$ である。$c$ が零以下の整数でなければ、$a=-m$ で有限和になる。分母にも零が生じる場合は、どのパラメータを先に極限に取るかを指定する。

:::example[終端する Gauss 級数（1）]
$a=-m,b=2,c=1$ を使う。

$$
m=1
$$

級数を多項式にする。

$$
\begin{aligned}
F(-m,2;1;x)\big|_{m=1} &= 1 - 2 x
\end{aligned}
$$

一階項の係数は $c-(a+b+1)x$ である。

$$
\begin{aligned}
c-(a+b+1)x\big|_{m=1} &= 1 - 2 x
\end{aligned}
$$

零階項の係数を読む。

$$
\begin{aligned}
-ab\big|_{m=1} &= 2
\end{aligned}
$$

$x=\frac{1}{4}$ で値と導関数を評価する。

$$
\begin{aligned}
F(-m,2;1;x)\big|_{x=\frac{1}{4}} &= \frac{1}{2} \\
\frac{\partial^1 (F(-m,2;1;x))}{\partial x^1} &= -2 \\
\left.\frac{\partial^1 (F(-m,2;1;x))}{\partial x^1}\right|_{x=\frac{1}{4}} &= -2
\end{aligned}
$$

$x=\frac{1}{4}$ で値と導関数を評価する。

$$
\begin{aligned}
c-(a+b+1)x\big|_{x=\frac{1}{4}} &= \frac{1}{2} \\
\frac{\partial^1 (c-(a+b+1)x)}{\partial x^1} &= -2 \\
\left.\frac{\partial^1 (c-(a+b+1)x)}{\partial x^1}\right|_{x=\frac{1}{4}} &= -2
\end{aligned}
$$
:::

:::example[終端する Gauss 級数（2）]
$a=-m,b=2,c=1$ を使う。次数を一つ上げる。

$$
m=2
$$

級数を多項式にする。

$$
\begin{aligned}
F(-m,2;1;x)\big|_{m=2} &= \left(x - 1\right) \left(3 x - 1\right)
\end{aligned}
$$

一階項の係数は $c-(a+b+1)x$ である。

$$
\begin{aligned}
c-(a+b+1)x\big|_{m=2} &= 1 - x
\end{aligned}
$$

零階項の係数を読む。

$$
\begin{aligned}
-ab\big|_{m=2} &= 4
\end{aligned}
$$

$x=\frac{1}{4}$ で値と導関数を評価する。

$$
\begin{aligned}
F(-m,2;1;x)\big|_{x=\frac{1}{4}} &= \frac{3}{16} \\
\frac{\partial^1 (F(-m,2;1;x))}{\partial x^1} &= 6 x - 4 \\
\left.\frac{\partial^1 (F(-m,2;1;x))}{\partial x^1}\right|_{x=\frac{1}{4}} &= - \frac{5}{2}
\end{aligned}
$$

$x=\frac{1}{4}$ で値と導関数を評価する。

$$
\begin{aligned}
c-(a+b+1)x\big|_{x=\frac{1}{4}} &= \frac{3}{4} \\
\frac{\partial^1 (c-(a+b+1)x)}{\partial x^1} &= -1 \\
\left.\frac{\partial^1 (c-(a+b+1)x)}{\partial x^1}\right|_{x=\frac{1}{4}} &= -1
\end{aligned}
$$
:::

:::example[終端する Gauss 級数（3）]
$a=-m,b=2,c=1$ を使う。さらに次数を上げて、低次の結果と照合する。

$$
m=3
$$

級数を多項式にする。

$$
\begin{aligned}
F(-m,2;1;x)\big|_{m=3} &= \left(1 - 4 x\right) \left(x - 1\right)^{2}
\end{aligned}
$$

一階項の係数は $c-(a+b+1)x$ である。

$$
\begin{aligned}
c-(a+b+1)x\big|_{m=3} &= 1
\end{aligned}
$$

零階項の係数を読む。

$$
\begin{aligned}
-ab\big|_{m=3} &= 6
\end{aligned}
$$

$x=\frac{1}{4}$ で値と導関数を評価する。

$$
\begin{aligned}
F(-m,2;1;x)\big|_{x=\frac{1}{4}} &= 0 \\
\frac{\partial^1 (F(-m,2;1;x))}{\partial x^1} &= 6 \left(1 - 2 x\right) \left(x - 1\right) \\
\left.\frac{\partial^1 (F(-m,2;1;x))}{\partial x^1}\right|_{x=\frac{1}{4}} &= - \frac{9}{4}
\end{aligned}
$$
:::

:::example[終端する Gauss 級数（4）]
$a=-m,b=2,c=1$ を使う。四次の結果まで具体化する。

$$
m=4
$$

級数を多項式にする。

$$
\begin{aligned}
F(-m,2;1;x)\big|_{m=4} &= \left(x - 1\right)^{3} \left(5 x - 1\right)
\end{aligned}
$$

一階項の係数は $c-(a+b+1)x$ である。

$$
\begin{aligned}
c-(a+b+1)x\big|_{m=4} &= x + 1
\end{aligned}
$$

零階項の係数を読む。

$$
\begin{aligned}
-ab\big|_{m=4} &= 8
\end{aligned}
$$

$x=\frac{1}{4}$ で値と導関数を評価する。

$$
\begin{aligned}
F(-m,2;1;x)\big|_{x=\frac{1}{4}} &= - \frac{27}{256} \\
\frac{\partial^1 (F(-m,2;1;x))}{\partial x^1} &= \left(x - 1\right)^{2} \left(20 x - 8\right) \\
\left.\frac{\partial^1 (F(-m,2;1;x))}{\partial x^1}\right|_{x=\frac{1}{4}} &= - \frac{27}{16}
\end{aligned}
$$

$x=\frac{1}{4}$ で値と導関数を評価する。

$$
\begin{aligned}
c-(a+b+1)x\big|_{x=\frac{1}{4}} &= \frac{5}{4} \\
\frac{\partial^1 (c-(a+b+1)x)}{\partial x^1} &= 1 \\
\left.\frac{\partial^1 (c-(a+b+1)x)}{\partial x^1}\right|_{x=\frac{1}{4}} &= 1
\end{aligned}
$$
:::

### 一般化超幾何級数の収束域

分子パラメータを $a_1,\ldots,a_p$、分母パラメータを $b_1,\ldots,b_q$ とし、$b_j$ は零以下の整数を避ける。
$$
{}_pF_q(\mathbf a;\mathbf b;z)=\sum_{n\ge0}
 \frac{\prod_{i=1}^p(a_i)_n}{\prod_{j=1}^q(b_j)_n}\frac{z^n}{n!}
$$
を :term[一般化超幾何級数]{en="generalized hypergeometric series" yomi="いっぱんかちょうきかきゅうすう"} という。分子に非正整数があれば有限次数で終わる。終端しない場合、隣接項の比は $z\,\prod_i(n+a_i)/[(n+1)\prod_j(n+b_j)]$ で、その大きさは $|z|n^{p-q-1}$ に漸近する。従って $p\le q$ なら全関数、$p=q+1$ なら収束半径 $1$、$p>q+1$ なら $z\ne0$ で発散する。

$\vartheta=z\partial_z$ と書くと、$\vartheta z^n=nz^n$ である。係数の比を次数ごとに書けば
$$
\left[\vartheta\prod_j(\vartheta+b_j-1)-z\prod_i(\vartheta+a_i)\right]y=0
$$
となる。実際、第 1 項で $c_n n\prod_j(n+b_j-1)$、第 2 項で $c_{n-1}\prod_i(n-1+a_i)$ が $z^n$ の係数になり、再帰で一致する。$z^0$ の第 1 項は $\vartheta$ で消える。

例えば $\,{}_0F_0(z)=e^z$ であり、$\,{}_1F_0(1;;z)$ は幾何級数 $1/(1-z)$ になる。後者の引数が $-1/2$ なら項が $(-1/2)^n$、和が $2/3$ である。引数が $1/2$ なら和は $2$ になる。さらに $\,{}_2F_0(1,1;;z)$ は係数が $n!$ なので任意の $z\ne0$ で発散する。この最後の係数の増加は [漸近展開と級数](/physics/mathematical-physics/asymptotics/) で扱う Borel 和の出発点である。

### 収束円周上の判定

$p=q+1$ で級数が終端しないとし、$d=\sum_j b_j-\sum_i a_i$ と置く。$z^n$ の係数 $c_n$ は、Gamma の漸化式と Stirling の公式から $c_n=Cn^{-d-1}(1+O(n^{-1}))$、$C=\prod_j\Gamma(b_j)/\prod_i\Gamma(a_i)\ne0$ となる。従って $|z|=1$ では $\Re d>0$ のとき、かつこのときに限り絶対収束する。$z\ne1$ の場合、$\Re d>-1$ なら $c_n\to0$ で、係数比から $c_{n+1}-c_n=O(n^{-\Re d-2})$ である。幾何和 $\sum_{j\le n}z^j$ の有界性と部分和分によって収束する。$\Re d\le-1$ では一般項が零に近づかない。

$z=1$ では $\Re d\le0$ のときも収束しないことを、区間 $N\le n\le tN$（固定した $t>1$）の部分和で確かめられる。係数の漸近式をこの区間で用いると、その和を $CN^{-d}$ で割った値は、Riemann 和として $\int_1^t u^{-d-1}\dd u$ に収束する。$t$ を $1$ に十分近く選べば積分は非零である。$\Re d\le0$ では $|N^{-d}|\ge1$ なので Cauchy の収束条件が破れる。例えば $F(1,1;2;z)$ は $d=0$ で、$z=1$ では調和級数、$z=-1$ では収束する交代級数になる。

## 超幾何関数の接続公式

### 別の特異点で同じ解を表す

$\Re c>\Re b>0$ なら、二項級数と式 :eqref[eq-mp-beta-gamma] により
$$
{}_2F_1(a,b;c;z)=
\frac{\Gamma(c)}{\Gamma(b)\Gamma(c-b)}\int_0^1 t^{b-1}(1-t)^{c-b-1}(1-zt)^{-a}\dd t.
$$
$|z|<1$ のコンパクト集合で級数は一様収束し、主枝で切断 $[1,\infty)$ を避ければ積分で接続できる。$\Re(c-a-b)>0$ で実軸上を $z\uparrow1$ とすると、優収束と Beta 積分から Gauss の値 $A=\Gamma(c)\Gamma(c-a-b)/[\Gamma(c-a)\Gamma(c-b)]$ になる。

$d=c-a-b\notin\Z$ とすると、$s=1-z$ と置き、$0<|s|<1$、$-\pi<\arg s<\pi$ で方程式を解き直して
$$
\begin{aligned}
{}_2F_1(a,b;c;z)
&=A\,{}_2F_1(a,b;1-d;s)\\
&\quad+B s^d\,{}_2F_1(c-a,c-b;1+d;s),\\
B&=\frac{\Gamma(c)\Gamma(-d)}{\Gamma(a)\Gamma(b)} .
\end{aligned}
\label{eq-mp-hypergeometric-connection}
$$
二つの局所解は決定根が異なるので独立である。$\Re d>0$ では $s\to0$ で第 2 項が消え、定数 $A$ は上の値になる。他方、$y=(1-z)^d u$ を代入すると $(c-a,c-b;c)$ の同じ型の式が得られ、原点で両者は $1$ の解析解なので $F(a,b;c;z)=(1-z)^d F(c-a,c-b;c;z)$ である。$\Re d<0$ でこれを $s^d$ で割って Gauss の値を取ると $B$ が求まる。局所解の係数は任意の共通正則点で値と微分を合わせる二元連立であり、パラメータに有理型なので、別々に求めた係数を同じ領域へ解析接続できる。整数差では両項の極が相殺する極限を取り、$s^d$ の指数の微分から対数が生じ得る。例えば $F(1,1;2;z)$ は $d=0$ で実際に対数を持つ。

### Pfaff の変換で負の引数を移す

Euler の積分で $t=1-v$ とすれば $1-zt=(1-z)[1-zv/(z-1)]$ である。従って
$$
{}_2F_1(a,b;c;z)=(1-z)^{-a}\,
{}_2F_1(a,c-b;c;z/(z-1)).
$$
積分で証明した範囲から一意接続すれば枝の一致する領域で成り立つ。負の実軸では $z/(z-1)\in(0,1)$ となるので、元の級数が使えない大きさの負の引数も、収束する側から評価できる。例えば $F(1,1;2;-2)=\log3/2$ であり、右辺でも $(1/3)F(1,1;2;2/3)=(1/3)(3\log3/2)$ と一致する。

### 無限遠の二つの冪への接続

$b-a\notin\Z$ とする。Pfaff の変換の右辺の引数 $w=z/(z-1)$ は、無限遠で $1$ に近づく。式 :eqref[eq-mp-hypergeometric-connection] を $(a,c-b;c;w)$ に適用すると、その指数差は $b-a$ となる。$1-w=1/(1-z)$ であり、両方の級数にもう一度 Pfaff の変換を適用すると、引数は $(1-w)/(-w)=1/z$ になる。まず $z<-1$ の実軸上で冪をすべて正の底で合わせると
$$
\begin{aligned}
 {}_2F_1(a,b;c;z)
 &=\frac{\Gamma(c)\Gamma(b-a)}{\Gamma(b)\Gamma(c-a)}
   (-z)^{-a}{}_2F_1(a,1-c+a;1-b+a;1/z)\\
 &\quad+\frac{\Gamma(c)\Gamma(a-b)}{\Gamma(a)\Gamma(c-b)}
   (-z)^{-b}{}_2F_1(b,1-c+b;1-a+b;1/z)
\end{aligned}
$$
となる。ここで $c$ は非正整数を避け、分母の Gamma の極は逆 Gamma の零として読む。この等式は $|z|>1$、$|\arg(-z)|<\pi$ に主枝で延びる。二つの冪は無限遠の決定根 $(a,b)$ に対応する。例えば $F(1/2,1;1;z)=(1-z)^{-1/2}$ では第 2 係数が $1/\Gamma(0)=0$ を含んで消え、第 1 冪だけが残る。$b-a$ が整数の場合は右辺を和のまま極限に取ると、原点側で接続した関数が得られる。

### 分岐点での対数

接続点の指数差が整数でなければ接続係数は Gamma の商で書ける。整数差に近づけると、個々の係数が発散して和が有限になることがある。二つの冪の差を指数差で割った極限から対数が生じる。複素領域へ延ばすときは対数の切断を固定する。

:::example[Euler 積分からの接続（1）]
$F(1,m;m+1;x)=m\int_0^1 t^{m-1}/(1-xt)\dd t$ として $v=xt$ を入れる。

$$
m=1
$$

分子を分母で割ると有限和と対数になる。

$$
\begin{aligned}
F(1,m;m+1;x)\big|_{m=1} &= - \frac{\log{\left(1 - x \right)}}{x}
\end{aligned}
$$

分岐点付近では $x^{-m}\to1$ である。

$$
\begin{aligned}
[\log(1-x)]F\big|_{m=1} &= -1
\end{aligned}
$$

対数を差し引いた極限は調和和になる。

$$
\begin{aligned}
\lim_{x\to1^-}(F+m\log(1-x))\big|_{m=1} &= 0
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
F(1,m;m+1;x)\big|_{x=\frac{1}{2}} &= \log{\left(4 \right)}
\end{aligned}
$$
:::

:::example[Euler 積分からの接続（2）]
$F(1,m;m+1;x)=m\int_0^1 t^{m-1}/(1-xt)\dd t$ として $v=xt$ を入れる。次数を一つ上げる。

$$
m=2
$$

分子を分母で割ると有限和と対数になる。

$$
\begin{aligned}
F(1,m;m+1;x)\big|_{m=2} &= \frac{- 2 x - 2 \log{\left(1 - x \right)}}{x^{2}}
\end{aligned}
$$

分岐点付近では $x^{-m}\to1$ である。

$$
\begin{aligned}
[\log(1-x)]F\big|_{m=2} &= -2
\end{aligned}
$$

対数を差し引いた極限は調和和になる。

$$
\begin{aligned}
\lim_{x\to1^-}(F+m\log(1-x))\big|_{m=2} &= -2
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
F(1,m;m+1;x)\big|_{x=\frac{1}{2}} &= -4 + \log{\left(256 \right)}
\end{aligned}
$$
:::

:::example[Euler 積分からの接続（3）]
$F(1,m;m+1;x)=m\int_0^1 t^{m-1}/(1-xt)\dd t$ として $v=xt$ を入れる。さらに次数を上げて、低次の結果と照合する。

$$
m=3
$$

分子を分母で割ると有限和と対数になる。

$$
\begin{aligned}
F(1,m;m+1;x)\big|_{m=3} &= \frac{- \frac{3 x^{2}}{2} - 3 x - 3 \log{\left(1 - x \right)}}{x^{3}}
\end{aligned}
$$

分岐点付近では $x^{-m}\to1$ である。

$$
\begin{aligned}
[\log(1-x)]F\big|_{m=3} &= -3
\end{aligned}
$$

対数を差し引いた極限は調和和になる。

$$
\begin{aligned}
\lim_{x\to1^-}(F+m\log(1-x))\big|_{m=3} &= - \frac{9}{2}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
F(1,m;m+1;x)\big|_{x=\frac{1}{2}} &= -15 + \log{\left(16777216 \right)}
\end{aligned}
$$
:::

:::example[Euler 積分からの接続（4）]
$F(1,m;m+1;x)=m\int_0^1 t^{m-1}/(1-xt)\dd t$ として $v=xt$ を入れる。四次の結果まで具体化する。

$$
m=4
$$

分子を分母で割ると有限和と対数になる。

$$
\begin{aligned}
F(1,m;m+1;x)\big|_{m=4} &= \frac{- \frac{4 x^{3}}{3} - 2 x^{2} - 4 x - 4 \log{\left(1 - x \right)}}{x^{4}}
\end{aligned}
$$

分岐点付近では $x^{-m}\to1$ である。

$$
\begin{aligned}
[\log(1-x)]F\big|_{m=4} &= -4
\end{aligned}
$$

対数を差し引いた極限は調和和になる。

$$
\begin{aligned}
\lim_{x\to1^-}(F+m\log(1-x))\big|_{m=4} &= - \frac{22}{3}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
F(1,m;m+1;x)\big|_{x=\frac{1}{2}} &= - \frac{128}{3} + \log{\left(18446744073709551616 \right)}
\end{aligned}
$$
:::

## 合流型超幾何関数（Kummer 関数）

### 特異点を合流させる

Gauss 級数で $b$ を正の実数として $+\infty$ に送り、引数を $z/b$ とすると $(b)_n/b^n\to1$ である。項別極限は $|z|\le R$ で正当化できる。実際、十分大きい $n,b$ では項の比の絶対値が $2R(1/n+1/b)\le1/2$ になる。先頭が有限で止まる場合は有限和の極限でよい。$c$ が非正整数でないとき、こうして得られる $M(a,c,z)=\sum (a)_n z^n/[(c)_n n!]$ は全平面で収束し、$zy''+(c-z)y'-ay=0$ を満たす。同じ方程式は元の Gauss の方程式を拡大して係数の極限を取っても得られる。$z=0$ のもう一つの根は $1-c$ だから $c\notin\Z$ なら $z^{1-c}M(a-c+1,2-c,z)$ が独立解となる。

$\Re a>0,\Re z>0$ で
$$
U(a,c,z)=\Gamma(a)^{-1}\int_0^\infty e^{-zt}t^{a-1}(1+t)^{c-a-1}\dd t
$$
と定める。$\Gamma(a)$ を外した被積分関数に $z\partial_z^2+(c-z)\partial_z-a$ を作用させると $-\partial_t[e^{-zt}t^a(1+t)^{c-a}]$ になるので端点で消える。実軸上で $z\to+\infty$ とし $u=zt$ と置くと、$\Re a>0$ のもとで優収束により $z^aU\to1$ となる。これは指数の枝と区別する代数的な規格化である。$\Re a>0$ なら減衰し、例えば $a=0$ への延長は $U(0,c,z)=1$ である。例えば $U(1,2,z)=1/z$、$a\notin\{0,-1,-2,\ldots\}$ なら $M(a,a,z)=e^z$ である。$M(-n,\alpha+1,z)$ は係数比較で $L_n^{(\alpha)}(z)n!/(\alpha+1)_n$ に等しい。

### 原点の二解に減衰側を合わせる

$c\notin\Z$ で $U=C_0 M(a,c,z)+C_1 z^{1-c}M(a-c+1,2-c,z)$ とする。$\Re a>0,\Re c<1$ では $z\downarrow0$ に積分の極限があり、$t=v/(1-v)$ を使うと $C_0=\Gamma(1-c)/\Gamma(a-c+1)$ となる。$\Re c>1$ では $z^{c-1}U$ の極限を取る。$u=zt$ と置いた積分で $u>0$ の点ごとの極限は $\Gamma(a)^{-1}e^{-u}u^{c-2}$ である。$t$ の有界部分は掛けた $z^{c-1}$ によって消え、$t\ge1$ では $(1+t)^{c-a-1}$ を $t^{c-a-1}$ の定数倍で挟めるので優収束が使える。従って $C_1=\Gamma(c-1)/\Gamma(a)$ である。重なりの正則点で値と微分を合わせる式に戻れば両係数はパラメータに有理型なので、最初の仮定から接続できる。整数 $c$ では二項全体の極限を取り、個々の Gamma 因子の極で式を打ち切らない。

### 無限遠での合流型の展開

積分内の $(1+t)^{c-a-1}$ を $t=0$ で展開する。$0\le t\le1/2$ で $N-1$ 次まで残すと剰余は $O(t^N)$ であり、$t>1/2$ での積分は $z$ に関して指数的に小さい。従って $\Re a>0$、$z\to+\infty$ のとき
$$
U(a,c,z)=z^{-a}\left[\sum_{j=0}^{N-1}
 \frac{(-1)^j(a)_j(a-c+1)_j}{j!\,z^j}
 +O(z^{-N})\right].
$$
係数は $t^j$ の二項係数 $(-1)^j(a-c+1)_j/j!$ に積分の Gamma 比 $(a)_j$ を掛けたものである。例えば $U(1,2,z)=1/z$ では $(a-c+1)_j=(0)_j$ により補正項は全て零である。一方 $U(1,1,z)$ では係数が $(-1)^j j!$ なので級数は発散する。これは漸近展開の各有限打ち切りの誤差評価とは両立する。

規格化を任意の $a$ へ延ばすには、同じ積分を $t^a/a$ で部分積分して得る恒等式
$$
 U(a,c,z)=z U(a+1,c+1,z)+(a+1-c)U(a+1,c,z)
$$
を用いる。十分回数だけ右へずらせば全積分に $\Re a>0$ が適用できる。右辺は接続公式の解析接続と一致し、第一項が先頭の $z^{-a}$ を与え、第二項は一段小さい。この操作によって $U(0,c,z)=1$ も分かる。実際、$U(1,c,z)$ の積分を使うと $z U(1,c+1,z)+(1-c)U(1,c,z)$ は $-\partial_t[e^{-zt}(1+t)^{c-1}]$ の積分になり、値は $1$ である。

### Laguerre の規格化と合流

係数を比較すると $M(-m;b;x)=m!L_m^{(b-1)}(x)/(b)_m$ になる。同じ方程式と同じ原点値を使って照合することもできる。動径問題ではこの多項式に指数と冪を掛け、無限遠と原点の条件を組み込む。

:::example[Kummer 多項式の規格化（1）]
$b=2$ で積分ノルムを比べる。

$$
m=1
$$

原点で $1$ になる解を作る。

$$
\begin{aligned}
M(-m,2,x)\big|_{m=1} &= 1 - \frac{x}{2}
\end{aligned}
$$

Laguerre の原点値を求める。

$$
\begin{aligned}
L_m^{(1)}(0)\big|_{m=1} &= 2
\end{aligned}
$$

二乗積分には倍率の二乗が入る。

$$
\begin{aligned}
\int_0^\infty |M(-m,2,x)|^2 x e^{-x}\dd x\big|_{m=1} &= \frac{1}{2}
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
M(-m,2,x)\big|_{x=0} &= 1 \\
\frac{\partial^1 (M(-m,2,x))}{\partial x^1} &= - \frac{1}{2} \\
\left.\frac{\partial^1 (M(-m,2,x))}{\partial x^1}\right|_{x=0} &= - \frac{1}{2} \\
\frac{\partial^2 (M(-m,2,x))}{\partial x^2} &= 0 \\
\left.\frac{\partial^2 (M(-m,2,x))}{\partial x^2}\right|_{x=0} &= 0
\end{aligned}
$$
:::

:::example[Kummer 多項式の規格化（2）]
$b=2$ で積分ノルムを比べる。次数を一つ上げる。

$$
m=2
$$

原点で $1$ になる解を作る。

$$
\begin{aligned}
M(-m,2,x)\big|_{m=2} &= \frac{x^{2}}{6} - x + 1
\end{aligned}
$$

Laguerre の原点値を求める。

$$
\begin{aligned}
L_m^{(1)}(0)\big|_{m=2} &= 3
\end{aligned}
$$

二乗積分には倍率の二乗が入る。

$$
\begin{aligned}
\int_0^\infty |M(-m,2,x)|^2 x e^{-x}\dd x\big|_{m=2} &= \frac{1}{3}
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
M(-m,2,x)\big|_{x=0} &= 1 \\
\frac{\partial^1 (M(-m,2,x))}{\partial x^1} &= \frac{x}{3} - 1 \\
\left.\frac{\partial^1 (M(-m,2,x))}{\partial x^1}\right|_{x=0} &= -1 \\
\frac{\partial^2 (M(-m,2,x))}{\partial x^2} &= \frac{1}{3} \\
\left.\frac{\partial^2 (M(-m,2,x))}{\partial x^2}\right|_{x=0} &= \frac{1}{3}
\end{aligned}
$$
:::

:::example[Kummer 多項式の規格化（3）]
$b=2$ で積分ノルムを比べる。さらに次数を上げて、低次の結果と照合する。

$$
m=3
$$

原点で $1$ になる解を作る。

$$
\begin{aligned}
M(-m,2,x)\big|_{m=3} &= - \frac{x^{3}}{24} + \frac{x^{2}}{2} - \frac{3 x}{2} + 1
\end{aligned}
$$

Laguerre の原点値を求める。

$$
\begin{aligned}
L_m^{(1)}(0)\big|_{m=3} &= 4
\end{aligned}
$$

二乗積分には倍率の二乗が入る。

$$
\begin{aligned}
\int_0^\infty |M(-m,2,x)|^2 x e^{-x}\dd x\big|_{m=3} &= \frac{1}{4}
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
M(-m,2,x)\big|_{x=0} &= 1 \\
\frac{\partial^1 (M(-m,2,x))}{\partial x^1} &= - \frac{x^{2}}{8} + x - \frac{3}{2} \\
\left.\frac{\partial^1 (M(-m,2,x))}{\partial x^1}\right|_{x=0} &= - \frac{3}{2} \\
\frac{\partial^2 (M(-m,2,x))}{\partial x^2} &= 1 - \frac{x}{4} \\
\left.\frac{\partial^2 (M(-m,2,x))}{\partial x^2}\right|_{x=0} &= 1
\end{aligned}
$$
:::

:::example[Kummer 多項式の規格化（4）]
$b=2$ で積分ノルムを比べる。四次の結果まで具体化する。

$$
m=4
$$

原点で $1$ になる解を作る。

$$
\begin{aligned}
M(-m,2,x)\big|_{m=4} &= \frac{x^{4}}{120} - \frac{x^{3}}{6} + x^{2} - 2 x + 1
\end{aligned}
$$

Laguerre の原点値を求める。

$$
\begin{aligned}
L_m^{(1)}(0)\big|_{m=4} &= 5
\end{aligned}
$$

二乗積分には倍率の二乗が入る。

$$
\begin{aligned}
\int_0^\infty |M(-m,2,x)|^2 x e^{-x}\dd x\big|_{m=4} &= \frac{1}{5}
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
M(-m,2,x)\big|_{x=0} &= 1 \\
\frac{\partial^1 (M(-m,2,x))}{\partial x^1} &= \frac{x^{3}}{30} - \frac{x^{2}}{2} + 2 x - 2 \\
\left.\frac{\partial^1 (M(-m,2,x))}{\partial x^1}\right|_{x=0} &= -2 \\
\frac{\partial^2 (M(-m,2,x))}{\partial x^2} &= \frac{x^{2}}{10} - x + 2 \\
\left.\frac{\partial^2 (M(-m,2,x))}{\partial x^2}\right|_{x=0} &= 2
\end{aligned}
$$
:::

:::example[次数 1 の方程式への代入]
級数で求めた関数が元の方程式を満たすか、各微分項を分けて確かめる。係数を読み、試行関数を微分して方程式へ戻す。

この例の係数を左から二階、一階、零階の順に書く。右辺がある場合は三つの項の和と比較する。

$$
\left( x, \  2 - x, \  1\right)
$$

関数とその二つの導関数を並べる。特に一つ下の階数の導関数には、一階項の係数を掛ける。

$$
\begin{aligned}u&=1 - \frac{x}{2}\\u^{\prime}&=- \frac{1}{2}\\u^{\prime\prime}&=0\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=0+\left(\frac{x}{2} - 1\right)+\left(1 - \frac{x}{2}\right)\\&=0\end{aligned}
$$

$x$ が $0$ のところで値を評価する。変換した方程式を使う場合も、初期値や境界値は元の変数で照合する。

$$
\begin{aligned}u&=1\\ u^\prime&=- \frac{1}{2}\\ Lu&=0\end{aligned}
$$
:::
:::example[次数 2 の方程式への代入]
級数で求めた関数が元の方程式を満たすか、各微分項を分けて確かめる。係数を読み、試行関数を微分して方程式へ戻す。

この例の係数を左から二階、一階、零階の順に書く。右辺がある場合は三つの項の和と比較する。

$$
\left( x, \  2 - x, \  2\right)
$$

関数とその二つの導関数を並べる。特に一つ下の階数の導関数には、一階項の係数を掛ける。

$$
\begin{aligned}u&=\frac{x^{2}}{6} - x + 1\\u^{\prime}&=\frac{x}{3} - 1\\u^{\prime\prime}&=\frac{1}{3}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=\frac{x}{3}+\left(\left(2 - x\right) \left(\frac{x}{3} - 1\right)\right)+\left(\frac{x^{2}}{3} - 2 x + 2\right)\\&=0\end{aligned}
$$

$x$ が $0$ のところで値を評価する。変換した方程式を使う場合も、初期値や境界値は元の変数で照合する。

$$
\begin{aligned}u&=1\\ u^\prime&=-1\\ Lu&=0\end{aligned}
$$
:::
:::example[次数 3 の方程式への代入]
級数で求めた関数が元の方程式を満たすか、各微分項を分けて確かめる。係数を読み、試行関数を微分して方程式へ戻す。

この例の係数を左から二階、一階、零階の順に書く。右辺がある場合は三つの項の和と比較する。

$$
\left( x, \  2 - x, \  3\right)
$$

関数とその二つの導関数を並べる。特に一つ下の階数の導関数には、一階項の係数を掛ける。

$$
\begin{aligned}u&=- \frac{x^{3}}{24} + \frac{x^{2}}{2} - \frac{3 x}{2} + 1\\u^{\prime}&=- \frac{x^{2}}{8} + x - \frac{3}{2}\\u^{\prime\prime}&=1 - \frac{x}{4}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=x \left(1 - \frac{x}{4}\right)+\left(\left(2 - x\right) \left(- \frac{x^{2}}{8} + x - \frac{3}{2}\right)\right)+\left(- \frac{x^{3}}{8} + \frac{3 x^{2}}{2} - \frac{9 x}{2} + 3\right)\\&=0\end{aligned}
$$

$x$ が $0$ のところで値を評価する。変換した方程式を使う場合も、初期値や境界値は元の変数で照合する。

$$
\begin{aligned}u&=1\\ u^\prime&=- \frac{3}{2}\\ Lu&=0\end{aligned}
$$
:::
:::example[次数 4 の方程式への代入]
級数で求めた関数が元の方程式を満たすか、各微分項を分けて確かめる。係数を読み、試行関数を微分して方程式へ戻す。

この例の係数を左から二階、一階、零階の順に書く。右辺がある場合は三つの項の和と比較する。

$$
\left( x, \  2 - x, \  4\right)
$$

関数とその二つの導関数を並べる。特に一つ下の階数の導関数には、一階項の係数を掛ける。

$$
\begin{aligned}u&=\frac{x^{4}}{120} - \frac{x^{3}}{6} + x^{2} - 2 x + 1\\u^{\prime}&=\frac{x^{3}}{30} - \frac{x^{2}}{2} + 2 x - 2\\u^{\prime\prime}&=\frac{x^{2}}{10} - x + 2\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=x \left(\frac{x^{2}}{10} - x + 2\right)+\left(\left(2 - x\right) \left(\frac{x^{3}}{30} - \frac{x^{2}}{2} + 2 x - 2\right)\right)+\left(\frac{x^{4}}{30} - \frac{2 x^{3}}{3} + 4 x^{2} - 8 x + 4\right)\\&=0\end{aligned}
$$

$x$ が $0$ のところで値を評価する。変換した方程式を使う場合も、初期値や境界値は元の変数で照合する。

$$
\begin{aligned}u&=1\\ u^\prime&=-2\\ Lu&=0\end{aligned}
$$
:::

## Whittaker 関数と放物柱関数

### 一階項を消して形を整える

$y=e^{-z/2}z^{\mu+1/2}u$ を Kummer 方程式の解 $u$ から作る。$a=\mu-\kappa+1/2,c=1+2\mu$ とすると $y''+[-1/4+\kappa/z+(1/4-\mu^2)/z^2]y=0$ になる。実際、一階項の係数は $(2\mu+1)/z-1$ と一致し、残る係数から $a$ が定まる。$u=M,U$ を選んだものをそれぞれ $M_{\kappa,\mu},W_{\kappa,\mu}$ とする。例えば $W_{0,1/2}=e^{-z/2}$ は $\kappa=0,\mu=1/2$ で $U(1,2,z)=1/z$ から従う。

$y''+(\nu+1/2-x^2/4)y=0$ の放物柱関数の一つを $D_\nu(x)$ とする。正の $x$ では $D_\nu(x)=2^{\nu/2}e^{-x^2/4}U(-\nu/2,1/2,x^2/2)$ と選び、$+\infty$ で $x^\nu e^{-x^2/4}$ と規格化する。$y=e^{-x^2/4}u(x^2/2)$ と置くと Kummer の係数が上の値になる。最初は $\Re\nu<0$ で積分表示を使い、他はパラメータで接続する。整数 $n\ge0$ では $2^{-n/2}e^{-x^2/4}H_n(x/\sqrt2)$ が同じ式と規格化を持つ。実次数の減衰解は定数倍を除いて一つである。実際、十分大きい $x$ では $D_\nu$ は先頭の漸近形により非零であり、第 2 解 $D_\nu(x)\int^x D_\nu(t)^{-2}\dd t$ は $\int^x t^{-2\nu}e^{t^2/2}\dd t\sim x^{-2\nu-1}e^{x^2/2}$ から増大する。この積分の漸近形は導関数の比から得られる。従って整数次数で示した Hermite の解は同じ規格化により $D_n$ と一致する。

### 一階微分を除いた方程式

一般の $y^{\prime\prime}+P y^\prime+Qy=0$ で $y=e^{-\int P/2}u$ と置くと、一階項は $2(-P/2)+P=0$ で消える。零階項は $Q-P^\prime/2-P^2/4$ になる。放物柱方程式では Gaussian を外した後が Hermite 型になる。

:::example[放物柱関数の整数次数（1）]
減衰解を Hermite 多項式で表す。

$$
m=1
$$

引数を $\sqrt2$ で割って規格化を合わせる。

$$
\begin{aligned}
D_m(x)\big|_{m=1} &= x e^{- \frac{x^{2}}{4}}
\end{aligned}
$$

Gaussian を外した多項式の最高次は $1$ である。

$$
\begin{aligned}
e^{x^2/4}D_m(x)\big|_{m=1} &= x
\end{aligned}
$$

方程式で掛かる零階項を読む。

$$
\begin{aligned}
V_{\rm eff}(x)\big|_{m=1} &= \frac{3}{2} - \frac{x^{2}}{4}
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
D_m(x)\big|_{x=0} &= 0 \\
\frac{\partial^1 (D_m(x))}{\partial x^1} &= \frac{\left(2 - x^{2}\right) e^{- \frac{x^{2}}{4}}}{2} \\
\left.\frac{\partial^1 (D_m(x))}{\partial x^1}\right|_{x=0} &= 1 \\
\frac{\partial^2 (D_m(x))}{\partial x^2} &= \frac{x \left(x^{2} - 6\right) e^{- \frac{x^{2}}{4}}}{4} \\
\left.\frac{\partial^2 (D_m(x))}{\partial x^2}\right|_{x=0} &= 0
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
e^{x^2/4}D_m(x)\big|_{x=0} &= 0 \\
\frac{\partial^1 (e^{x^2/4}D_m(x))}{\partial x^1} &= 1 \\
\left.\frac{\partial^1 (e^{x^2/4}D_m(x))}{\partial x^1}\right|_{x=0} &= 1 \\
\frac{\partial^2 (e^{x^2/4}D_m(x))}{\partial x^2} &= 0 \\
\left.\frac{\partial^2 (e^{x^2/4}D_m(x))}{\partial x^2}\right|_{x=0} &= 0
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
V_{\rm eff}(x)\big|_{x=0} &= \frac{3}{2} \\
\frac{\partial^1 (V_{\rm eff}(x))}{\partial x^1} &= - \frac{x}{2} \\
\left.\frac{\partial^1 (V_{\rm eff}(x))}{\partial x^1}\right|_{x=0} &= 0 \\
\frac{\partial^2 (V_{\rm eff}(x))}{\partial x^2} &= - \frac{1}{2} \\
\left.\frac{\partial^2 (V_{\rm eff}(x))}{\partial x^2}\right|_{x=0} &= - \frac{1}{2}
\end{aligned}
$$
:::

:::example[放物柱関数の整数次数（2）]
減衰解を Hermite 多項式で表す。次数を一つ上げる。

$$
m=2
$$

引数を $\sqrt2$ で割って規格化を合わせる。

$$
\begin{aligned}
D_m(x)\big|_{m=2} &= \left(x^{2} - 1\right) e^{- \frac{x^{2}}{4}}
\end{aligned}
$$

Gaussian を外した多項式の最高次は $1$ である。

$$
\begin{aligned}
e^{x^2/4}D_m(x)\big|_{m=2} &= x^{2} - 1
\end{aligned}
$$

方程式で掛かる零階項を読む。

$$
\begin{aligned}
V_{\rm eff}(x)\big|_{m=2} &= \frac{5}{2} - \frac{x^{2}}{4}
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
D_m(x)\big|_{x=0} &= -1 \\
\frac{\partial^1 (D_m(x))}{\partial x^1} &= \frac{x \left(5 - x^{2}\right) e^{- \frac{x^{2}}{4}}}{2} \\
\left.\frac{\partial^1 (D_m(x))}{\partial x^1}\right|_{x=0} &= 0 \\
\frac{\partial^2 (D_m(x))}{\partial x^2} &= \frac{\left(x^{4} - 11 x^{2} + 10\right) e^{- \frac{x^{2}}{4}}}{4} \\
\left.\frac{\partial^2 (D_m(x))}{\partial x^2}\right|_{x=0} &= \frac{5}{2}
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
e^{x^2/4}D_m(x)\big|_{x=0} &= -1 \\
\frac{\partial^1 (e^{x^2/4}D_m(x))}{\partial x^1} &= 2 x \\
\left.\frac{\partial^1 (e^{x^2/4}D_m(x))}{\partial x^1}\right|_{x=0} &= 0 \\
\frac{\partial^2 (e^{x^2/4}D_m(x))}{\partial x^2} &= 2 \\
\left.\frac{\partial^2 (e^{x^2/4}D_m(x))}{\partial x^2}\right|_{x=0} &= 2
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
V_{\rm eff}(x)\big|_{x=0} &= \frac{5}{2} \\
\frac{\partial^1 (V_{\rm eff}(x))}{\partial x^1} &= - \frac{x}{2} \\
\left.\frac{\partial^1 (V_{\rm eff}(x))}{\partial x^1}\right|_{x=0} &= 0 \\
\frac{\partial^2 (V_{\rm eff}(x))}{\partial x^2} &= - \frac{1}{2} \\
\left.\frac{\partial^2 (V_{\rm eff}(x))}{\partial x^2}\right|_{x=0} &= - \frac{1}{2}
\end{aligned}
$$
:::

:::example[放物柱関数の整数次数（3）]
減衰解を Hermite 多項式で表す。さらに次数を上げて、低次の結果と照合する。

$$
m=3
$$

引数を $\sqrt2$ で割って規格化を合わせる。

$$
\begin{aligned}
D_m(x)\big|_{m=3} &= x \left(x^{2} - 3\right) e^{- \frac{x^{2}}{4}}
\end{aligned}
$$

Gaussian を外した多項式の最高次は $1$ である。

$$
\begin{aligned}
e^{x^2/4}D_m(x)\big|_{m=3} &= x \left(x^{2} - 3\right)
\end{aligned}
$$

方程式で掛かる零階項を読む。

$$
\begin{aligned}
V_{\rm eff}(x)\big|_{m=3} &= \frac{7}{2} - \frac{x^{2}}{4}
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
D_m(x)\big|_{x=0} &= 0 \\
\frac{\partial^1 (D_m(x))}{\partial x^1} &= \frac{\left(- x^{4} + 9 x^{2} - 6\right) e^{- \frac{x^{2}}{4}}}{2} \\
\left.\frac{\partial^1 (D_m(x))}{\partial x^1}\right|_{x=0} &= -3 \\
\frac{\partial^2 (D_m(x))}{\partial x^2} &= \frac{x \left(x^{4} - 17 x^{2} + 42\right) e^{- \frac{x^{2}}{4}}}{4} \\
\left.\frac{\partial^2 (D_m(x))}{\partial x^2}\right|_{x=0} &= 0
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
e^{x^2/4}D_m(x)\big|_{x=0} &= 0 \\
\frac{\partial^1 (e^{x^2/4}D_m(x))}{\partial x^1} &= 3 x^{2} - 3 \\
\left.\frac{\partial^1 (e^{x^2/4}D_m(x))}{\partial x^1}\right|_{x=0} &= -3 \\
\frac{\partial^2 (e^{x^2/4}D_m(x))}{\partial x^2} &= 6 x \\
\left.\frac{\partial^2 (e^{x^2/4}D_m(x))}{\partial x^2}\right|_{x=0} &= 0
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
V_{\rm eff}(x)\big|_{x=0} &= \frac{7}{2} \\
\frac{\partial^1 (V_{\rm eff}(x))}{\partial x^1} &= - \frac{x}{2} \\
\left.\frac{\partial^1 (V_{\rm eff}(x))}{\partial x^1}\right|_{x=0} &= 0 \\
\frac{\partial^2 (V_{\rm eff}(x))}{\partial x^2} &= - \frac{1}{2} \\
\left.\frac{\partial^2 (V_{\rm eff}(x))}{\partial x^2}\right|_{x=0} &= - \frac{1}{2}
\end{aligned}
$$
:::

:::example[放物柱関数の整数次数（4）]
減衰解を Hermite 多項式で表す。四次の結果まで具体化する。

$$
m=4
$$

引数を $\sqrt2$ で割って規格化を合わせる。

$$
\begin{aligned}
D_m(x)\big|_{m=4} &= \left(x^{4} - 6 x^{2} + 3\right) e^{- \frac{x^{2}}{4}}
\end{aligned}
$$

Gaussian を外した多項式の最高次は $1$ である。

$$
\begin{aligned}
e^{x^2/4}D_m(x)\big|_{m=4} &= x^{4} - 6 x^{2} + 3
\end{aligned}
$$

方程式で掛かる零階項を読む。

$$
\begin{aligned}
V_{\rm eff}(x)\big|_{m=4} &= \frac{9}{2} - \frac{x^{2}}{4}
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
D_m(x)\big|_{x=0} &= 3 \\
\frac{\partial^1 (D_m(x))}{\partial x^1} &= \frac{x \left(- x^{4} + 14 x^{2} - 27\right) e^{- \frac{x^{2}}{4}}}{2} \\
\left.\frac{\partial^1 (D_m(x))}{\partial x^1}\right|_{x=0} &= 0 \\
\frac{\partial^2 (D_m(x))}{\partial x^2} &= \frac{\left(x^{6} - 24 x^{4} + 111 x^{2} - 54\right) e^{- \frac{x^{2}}{4}}}{4} \\
\left.\frac{\partial^2 (D_m(x))}{\partial x^2}\right|_{x=0} &= - \frac{27}{2}
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
e^{x^2/4}D_m(x)\big|_{x=0} &= 3 \\
\frac{\partial^1 (e^{x^2/4}D_m(x))}{\partial x^1} &= 4 x \left(x^{2} - 3\right) \\
\left.\frac{\partial^1 (e^{x^2/4}D_m(x))}{\partial x^1}\right|_{x=0} &= 0 \\
\frac{\partial^2 (e^{x^2/4}D_m(x))}{\partial x^2} &= 12 x^{2} - 12 \\
\left.\frac{\partial^2 (e^{x^2/4}D_m(x))}{\partial x^2}\right|_{x=0} &= -12
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
V_{\rm eff}(x)\big|_{x=0} &= \frac{9}{2} \\
\frac{\partial^1 (V_{\rm eff}(x))}{\partial x^1} &= - \frac{x}{2} \\
\left.\frac{\partial^1 (V_{\rm eff}(x))}{\partial x^1}\right|_{x=0} &= 0 \\
\frac{\partial^2 (V_{\rm eff}(x))}{\partial x^2} &= - \frac{1}{2} \\
\left.\frac{\partial^2 (V_{\rm eff}(x))}{\partial x^2}\right|_{x=0} &= - \frac{1}{2}
\end{aligned}
$$
:::

:::example[次数 1 の方程式への代入]
級数で求めた関数が元の方程式を満たすか、各微分項を分けて確かめる。係数を読み、試行関数を微分して方程式へ戻す。

この例の係数を左から二階、一階、零階の順に書く。右辺がある場合は三つの項の和と比較する。

$$
\left( 1, \  0, \  \frac{3}{2} - \frac{x^{2}}{4}\right)
$$

関数とその二つの導関数を並べる。特に一つ下の階数の導関数には、一階項の係数を掛ける。

$$
\begin{aligned}u&=x e^{- \frac{x^{2}}{4}}\\u^{\prime}&=- \frac{x^{2} e^{- \frac{x^{2}}{4}}}{2} + e^{- \frac{x^{2}}{4}}\\u^{\prime\prime}&=\frac{x \left(x^{2} - 6\right) e^{- \frac{x^{2}}{4}}}{4}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=\frac{x \left(x^{2} - 6\right) e^{- \frac{x^{2}}{4}}}{4}+\left(0\right)+\left(x \left(\frac{3}{2} - \frac{x^{2}}{4}\right) e^{- \frac{x^{2}}{4}}\right)\\&=0\end{aligned}
$$

$x$ が $0$ のところで値を評価する。変換した方程式を使う場合も、初期値や境界値は元の変数で照合する。

$$
\begin{aligned}u&=0\\ u^\prime&=1\\ Lu&=0\end{aligned}
$$
:::
:::example[次数 2 の方程式への代入]
級数で求めた関数が元の方程式を満たすか、各微分項を分けて確かめる。係数を読み、試行関数を微分して方程式へ戻す。

この例の係数を左から二階、一階、零階の順に書く。右辺がある場合は三つの項の和と比較する。

$$
\left( 1, \  0, \  \frac{5}{2} - \frac{x^{2}}{4}\right)
$$

関数とその二つの導関数を並べる。特に一つ下の階数の導関数には、一階項の係数を掛ける。

$$
\begin{aligned}u&=\frac{\left(2 x^{2} - 2\right) e^{- \frac{x^{2}}{4}}}{2}\\u^{\prime}&=- \frac{x \left(2 x^{2} - 2\right) e^{- \frac{x^{2}}{4}}}{4} + 2 x e^{- \frac{x^{2}}{4}}\\u^{\prime\prime}&=\left(- 2 x^{2} + \frac{\left(x^{2} - 2\right) \left(x^{2} - 1\right)}{4} + 2\right) e^{- \frac{x^{2}}{4}}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=\left(- 2 x^{2} + \frac{\left(x^{2} - 2\right) \left(x^{2} - 1\right)}{4} + 2\right) e^{- \frac{x^{2}}{4}}+\left(0\right)+\left(\frac{\left(\frac{5}{2} - \frac{x^{2}}{4}\right) \left(2 x^{2} - 2\right) e^{- \frac{x^{2}}{4}}}{2}\right)\\&=0\end{aligned}
$$

$x$ が $0$ のところで値を評価する。変換した方程式を使う場合も、初期値や境界値は元の変数で照合する。

$$
\begin{aligned}u&=-1\\ u^\prime&=0\\ Lu&=0\end{aligned}
$$
:::
:::example[次数 3 の方程式への代入]
級数で求めた関数が元の方程式を満たすか、各微分項を分けて確かめる。係数を読み、試行関数を微分して方程式へ戻す。

この例の係数を左から二階、一階、零階の順に書く。右辺がある場合は三つの項の和と比較する。

$$
\left( 1, \  0, \  \frac{7}{2} - \frac{x^{2}}{4}\right)
$$

関数とその二つの導関数を並べる。特に一つ下の階数の導関数には、一階項の係数を掛ける。

$$
\begin{aligned}u&=\frac{\sqrt{2} \left(2 \sqrt{2} x^{3} - 6 \sqrt{2} x\right) e^{- \frac{x^{2}}{4}}}{4}\\u^{\prime}&=- \frac{\sqrt{2} x \left(2 \sqrt{2} x^{3} - 6 \sqrt{2} x\right) e^{- \frac{x^{2}}{4}}}{8} + \frac{\sqrt{2} \left(6 \sqrt{2} x^{2} - 6 \sqrt{2}\right) e^{- \frac{x^{2}}{4}}}{4}\\u^{\prime\prime}&=x \left(- 3 x^{2} + \frac{\left(x^{2} - 3\right) \left(x^{2} - 2\right)}{4} + 9\right) e^{- \frac{x^{2}}{4}}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=x \left(- 3 x^{2} + \frac{\left(x^{2} - 3\right) \left(x^{2} - 2\right)}{4} + 9\right) e^{- \frac{x^{2}}{4}}+\left(0\right)+\left(\frac{\sqrt{2} \left(\frac{7}{2} - \frac{x^{2}}{4}\right) \left(2 \sqrt{2} x^{3} - 6 \sqrt{2} x\right) e^{- \frac{x^{2}}{4}}}{4}\right)\\&=0\end{aligned}
$$

$x$ が $0$ のところで値を評価する。変換した方程式を使う場合も、初期値や境界値は元の変数で照合する。

$$
\begin{aligned}u&=0\\ u^\prime&=-3\\ Lu&=0\end{aligned}
$$
:::
:::example[次数 4 の方程式への代入]
級数で求めた関数が元の方程式を満たすか、各微分項を分けて確かめる。係数を読み、試行関数を微分して方程式へ戻す。

この例の係数を左から二階、一階、零階の順に書く。右辺がある場合は三つの項の和と比較する。

$$
\left( 1, \  0, \  \frac{9}{2} - \frac{x^{2}}{4}\right)
$$

関数とその二つの導関数を並べる。特に一つ下の階数の導関数には、一階項の係数を掛ける。

$$
\begin{aligned}u&=\frac{\left(4 x^{4} - 24 x^{2} + 12\right) e^{- \frac{x^{2}}{4}}}{4}\\u^{\prime}&=- \frac{x \left(4 x^{4} - 24 x^{2} + 12\right) e^{- \frac{x^{2}}{4}}}{8} + \frac{\left(16 x^{3} - 48 x\right) e^{- \frac{x^{2}}{4}}}{4}\\u^{\prime\prime}&=\left(- 4 x^{2} \left(x^{2} - 3\right) + 12 x^{2} + \frac{\left(x^{2} - 2\right) \left(x^{4} - 6 x^{2} + 3\right)}{4} - 12\right) e^{- \frac{x^{2}}{4}}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=\left(- 4 x^{2} \left(x^{2} - 3\right) + 12 x^{2} + \frac{\left(x^{2} - 2\right) \left(x^{4} - 6 x^{2} + 3\right)}{4} - 12\right) e^{- \frac{x^{2}}{4}}+\left(0\right)+\left(\frac{\left(\frac{9}{2} - \frac{x^{2}}{4}\right) \left(4 x^{4} - 24 x^{2} + 12\right) e^{- \frac{x^{2}}{4}}}{4}\right)\\&=0\end{aligned}
$$

$x$ が $0$ のところで値を評価する。変換した方程式を使う場合も、初期値や境界値は元の変数で照合する。

$$
\begin{aligned}u&=3\\ u^\prime&=0\\ Lu&=0\end{aligned}
$$
:::

## Airy 関数

### 線形な転回点の近く

$y''-xy=0$ の解は [級数解法](/physics/mathematical-physics/ode-series/) で二つ構成した。Airy 関数 $\mathrm{Ai}(x)$ は輪郭 $\mathcal C$ を $\arg t=-\pi/3$ から $\pi/3$ へ向けて
$$
\mathrm{Ai}(x)=\frac1{2\pi i}\int_{\mathcal C}e^{t^3/3-xt}\dd t
$$
と取る。両端で $t^3=-|t|^3$ のため積分内で二回微分でき、$y''-xy$ は $t$ 微分の積分となって消える。$x=0$ で二本の半直線を使い $v=|t|^3/3$ と置くと、反射公式を用いて $\mathrm{Ai}(0)=1/[3^{2/3}\Gamma(2/3)]$、$\mathrm{Ai}'(0)=-1/[3^{1/3}\Gamma(1/3)]$ となる。第 2 解 $\mathrm{Bi}$ は $\mathrm{Bi}(0)=\sqrt3\,\mathrm{Ai}(0)$、$\mathrm{Bi}'(0)=-\sqrt3\,\mathrm{Ai}'(0)$ で定める。この初期値により Wronskian は $1/\pi$ である。$x>0$ の $\mathrm{Ai}$ の減衰主項は $e^{-2x^{3/2}/3}/(2\sqrt\pi x^{1/4})$ である。輪郭を $t=\sqrt x+iv$ に変形すると指数は $-2x^{3/2}/3-\sqrt x\,v^2-iv^3/3$ となる。高さ $\pm R$ の水平線で二本の半直線とこの縦線を結ぶ。水平線上の $t=\sigma\pm iR$ では $\sqrt x\le\sigma\le R/\sqrt3$ なので $\Re(t^3/3-xt)\le-8\sigma R^2/9$ である。従って接続線分では $\Re(t^3/3-xt)$ が減衰するため Cauchy の定理が使え、$v=x^{-1/4}u$ の Gaussian 極限で主項を得る。

### 転回点の反対側で振動が出る

$x=-X$、$X>0$ とすると Airy の輪郭は、絶対収束する傾いた半直線を虚軸へ近づけた振動積分で
$$
\mathrm{Ai}(-X)=\frac1{2\pi}\int_\R e^{i(Xv-v^3/3)}\dd v
$$
と書ける。例えば虚軸に平行な $\Re t=\eta>0$ を経由すれば絶対値に $e^{-\eta v^2}$ があり、$\eta\downarrow0$ では無限遠の積分を位相の導関数で部分積分して同じ極限を得る。$v=\sqrt X\,u$ と置くと大きいパラメータは $X^{3/2}$、相は $u-u^3/3$ になる。停留点 $u=\pm1$ での相は $\pm2/3$、二階微分は $\mp2$ である。二点の寄与を足して
$$
\mathrm{Ai}(-X)=\frac{X^{-1/4}}{\sqrt\pi}
\cos(2X^{3/2}/3-\pi/4)+O(X^{-1})
$$
を得る。停留点の近くを除いた尾は、位相の導関数が $u^2$ 程度で増えるので一回の部分積分でこの誤差に入る。正側の指数減衰解が負側では二つの進行波の和に接続する。

### 増大解の規格化

$W(\mathrm{Ai},\mathrm{Bi})=1/\pi$ と正の実軸の $\mathrm{Ai}$ の漸近形を用いる。十分大きい $x$ では $\mathrm{Ai}$ は非零なので、適当な $x_0$ に対して
$$
 \frac{\mathrm{Bi}(x)}{\mathrm{Ai}(x)}
 =\frac{\mathrm{Bi}(x_0)}{\mathrm{Ai}(x_0)}
 +\frac1\pi\int_{x_0}^x\frac{\dd t}{\mathrm{Ai}(t)^2}
$$
となる。積分の被積分関数は $4\pi\sqrt t\,e^{4t^{3/2}/3}$ に漸近する。指数 $e^{4t^{3/2}/3}$ の導関数は $2\sqrt t\,e^{4t^{3/2}/3}$ なので、導関数の比の極限から積分は $2\pi e^{4x^{3/2}/3}$ に漸近する。従って $\mathrm{Bi}(x)\sim x^{-1/4}e^{2x^{3/2}/3}/\sqrt\pi$ である。正側の転回点から離れた領域では、同じ方程式の二つの枝が逆向きの指数を持つ。

### 複素回転で振動する第 2 解を得る

$\omega=e^{2\pi i/3}$ とすると、$\mathrm{Ai}(\omega x)$ と $\mathrm{Ai}(\omega^{-1}x)$ も Airy 方程式を満たす。実際、二階微分による係数 $\omega^{\pm2}$ に方程式から $\omega^{\pm1}x$ が掛かり、$\omega^{\pm3}=1$ になる。従って
$$
 \mathrm{Bi}(x)
 =e^{-i\pi/6}\mathrm{Ai}(\omega^{-1}x)
 +e^{i\pi/6}\mathrm{Ai}(\omega x)
$$
である。右辺を $x=0$ で評価すると係数の和は $\sqrt3$、一回微分後の係数の和は $-\sqrt3$ であり、先に定めた初期値と一致するからである。

$\mathrm{Ai}(z)$ の縦線による評価は、主枝の $|\arg z|\le\pi/3$ にも延びる。$t=\sqrt z+iv$ とすれば $\Re\sqrt z>0$ で、横につなぐ線分では三次項の実部が負の二次の高さに比例し、$-zt$ の実部は高さの一次程度である。従って線分を無限遠へ送ることができる。$z=R e^{i\vartheta}$、$v=R^{-1/4}u$ と置くと、残る積分の因子は $e^{-e^{i\vartheta/2}u^2}e^{-iR^{-3/4}u^3/3}$ である。第 2 因子の一次の補正は奇関数なので積分から消え、その次の剰余は Gaussian に $O(R^{-3/2}|u|^6)$ を掛けたものに抑えられる。複素 Gaussian の積分は正の実係数から $\Re(e^{i\vartheta/2})>0$ の範囲へ解析接続して求められる。このため
$$
 \mathrm{Ai}(z)=\frac{z^{-1/4}}{2\sqrt\pi}e^{-2z^{3/2}/3}
 \left(1+O(R^{-3/2})\right)
$$
となる。$x=-X$ では回転後の二点の偏角は $\mp\pi/3$ である。上式を二項に戻して足すと
$$
 \mathrm{Bi}(-X)=\frac{X^{-1/4}}{\sqrt\pi}
 \cos(2X^{3/2}/3+\pi/4)+O(X^{-7/4})
 \qquad (X\to+\infty)
$$
を得る。$\mathrm{Ai}(-X)$ の位相 $-\pi/4$ と合わせると、転回点の左では二解の振動が四分の一周期ずれている。

### Airy の二系列

Airy 方程式 $u^{\prime\prime}=xu$ では $a_2=0$、$a_{n+3}=a_n/((n+3)(n+2))$ となる。定数から始まる系列と一次から始まる系列はいずれも全平面で収束する。従って原点で有限という条件だけで遠方の減衰解は選べない。

:::example[Airy 級数の係数（1）]
初期値 $(1,0)$ と $(0,1)$ の二系列を比較する。

$$
m=1
$$

第一の系列に現れる冪は三の倍数である。

$$
\begin{aligned}
a_{3m}\big|_{m=1} &= \frac{1}{6}
\end{aligned}
$$

第二の系列の指数は一つずれる。

$$
\begin{aligned}
b_{3m+1}\big|_{m=1} &= \frac{1}{12}
\end{aligned}
$$

二階微分で前の次数へ戻る係数を求める。

$$
\begin{aligned}
(3m)(3m-1)a_{3m}\big|_{m=1} &= 1
\end{aligned}
$$
:::

:::example[Airy 級数の係数（2）]
初期値 $(1,0)$ と $(0,1)$ の二系列を比較する。次数を一つ上げる。

$$
m=2
$$

第一の系列に現れる冪は三の倍数である。

$$
\begin{aligned}
a_{3m}\big|_{m=2} &= \frac{1}{180}
\end{aligned}
$$

第二の系列の指数は一つずれる。

$$
\begin{aligned}
b_{3m+1}\big|_{m=2} &= \frac{1}{504}
\end{aligned}
$$

二階微分で前の次数へ戻る係数を求める。

$$
\begin{aligned}
(3m)(3m-1)a_{3m}\big|_{m=2} &= \frac{1}{6}
\end{aligned}
$$
:::

:::example[Airy 級数の係数（3）]
初期値 $(1,0)$ と $(0,1)$ の二系列を比較する。さらに次数を上げて、低次の結果と照合する。

$$
m=3
$$

第一の系列に現れる冪は三の倍数である。

$$
\begin{aligned}
a_{3m}\big|_{m=3} &= \frac{1}{12960}
\end{aligned}
$$

第二の系列の指数は一つずれる。

$$
\begin{aligned}
b_{3m+1}\big|_{m=3} &= \frac{1}{45360}
\end{aligned}
$$

二階微分で前の次数へ戻る係数を求める。

$$
\begin{aligned}
(3m)(3m-1)a_{3m}\big|_{m=3} &= \frac{1}{180}
\end{aligned}
$$
:::

:::example[Airy 級数の係数（4）]
初期値 $(1,0)$ と $(0,1)$ の二系列を比較する。四次の結果まで具体化する。

$$
m=4
$$

第一の系列に現れる冪は三の倍数である。

$$
\begin{aligned}
a_{3m}\big|_{m=4} &= \frac{1}{1710720}
\end{aligned}
$$

第二の系列の指数は一つずれる。

$$
\begin{aligned}
b_{3m+1}\big|_{m=4} &= \frac{1}{7076160}
\end{aligned}
$$

二階微分で前の次数へ戻る係数を求める。

$$
\begin{aligned}
(3m)(3m-1)a_{3m}\big|_{m=4} &= \frac{1}{12960}
\end{aligned}
$$
:::

## 楕円積分と Jacobi の楕円関数

### 振子の周期を積分する

$F(\phi,k)=\int_0^\phi(1-k^2\sin^2\theta)^{-1/2}\dd\theta$ を :term[第 1 種不完全楕円積分]{en="incomplete elliptic integral of the first kind" yomi="だいいっしゅふかんぜんだえんせきぶん"}、$K(k)=F(\pi/2,k)$ を完全楕円積分とする。第 2 種は逆平方根を平方根に替えた $E(\phi,k)$ である。$|k|<1$ で二項級数を積分し、$\int_0^{\pi/2}\sin^{2n}\theta\dd\theta=(\pi/2)(1/2)_n/n!$ を Beta 積分から得ると $K=(\pi/2)F(1/2,1/2;1;k^2)$、$E(\pi/2,k)=(\pi/2)F(-1/2,1/2;1;k^2)$ となる。ここで右辺の四引数の $F$ は ${}_2F_1$ である。

実 $0\le k<1$ で $u=F(\phi,k)$ を逆に解き、$\operatorname{sn}u=\sin\phi$、$\operatorname{cn}u=\cos\phi$、$\operatorname{dn}u=\sqrt{1-k^2\sin^2\phi}$ とする。$\phi'(u)=\operatorname{dn}u$ だから $\operatorname{sn}'=\operatorname{cn}\operatorname{dn}$、$\operatorname{cn}'=-\operatorname{sn}\operatorname{dn}$、$\operatorname{dn}'=-k^2\operatorname{sn}\operatorname{cn}$ である。例えば $k=0$ では $\sin u,\cos u,1$ に戻る。長さ $\ell$ の単振子では $\dot\theta^2=(4g/\ell)(k^2-\sin^2(\theta/2))$、$k=\sin(\theta_{\max}/2)$ をエネルギー保存から得る。$\sin(\theta/2)=k\sin\phi$ を正の四分周期に入れると周期は $4\sqrt{\ell/g}K(k)$（$\mathrm s$）である。小振幅の相対補正は $k^2/4+O(k^4)$ である。

### 楕円積分の第三種と振幅の区別

第 3 種の不完全楕円積分は
$$
\Pi(n;\phi,k)=\int_0^\phi\frac{\dd\theta}
{(1-n\sin^2\theta)\sqrt{1-k^2\sin^2\theta}}
$$
である。例えば実 $0\le k<1$、$n<1$、$0\le\phi\le\pi/2$ では積分に極がない。$\phi=\pi/2$ を第 3 種の完全積分とする。$n=0$ なら第 1 種に一致する。さらに $k=0$ で完全にすると $v=\tan\theta$ により $\int_0^\infty[1+(1-n)v^2]^{-1}\dd v=\pi/(2\sqrt{1-n})$ である。例えば $n=3/4$ で $\pi$ になる。$n\uparrow1$ で発散するのは分母のもう一つの因子が $\theta=\pi/2$ で零に近づくためである。法 $k$ を二乗した母数 $m=k^2$ で表を引く場合は、係数の $k^2$ と $m$ を同じものとして照合する必要がある。

### 実周期と振子の時刻表示

$0\le k<1$ では振幅 $\phi$ に関する被積分関数が正かつ $\pi$ 周期で、$\int_0^\pi(1-k^2\sin^2\theta)^{-1/2}\dd\theta=2K$ である。従って同じ法で $F(\phi+\pi,k)=F(\phi,k)+2K$、逆振幅は $\phi(u+2K)=\phi(u)+\pi$ である。よって $\operatorname{sn},\operatorname{cn}$ は $2K$ の平行移動で符号が変わり $4K$ 周期となる。$\operatorname{dn}$ は $2K$ 周期である。$k=0$ では $\operatorname{dn}$ が定数なので任意の平行移動で不変である。

初期条件 $\theta(0)=0$、$\dot\theta(0)=2k\sqrt{g/\ell}$ の振子では
$$
\theta(t)=2\arcsin\bigl(k\,\operatorname{sn}(\sqrt{g/\ell}\,t,k)\bigr)
$$
となる。実際、微分により $\dot\theta=2k\sqrt{g/\ell}\,\operatorname{cn}$ となる。分母の $\sqrt{1-k^2\operatorname{sn}^2}$ が $\operatorname{dn}$ と等しく相殺されるためである。さらに $\dd\operatorname{cn}/\dd u=-\operatorname{sn}\operatorname{dn}$ と $\sin\theta=2k\operatorname{sn}\operatorname{dn}$ を使えば $\ddot\theta=-(g/\ell)\sin\theta$ が直接確かめられる。この式から、速度が転回点で符号を変えることも $\operatorname{cn}$ の零点で分かる。

### Jacobi の加法公式

同じ法 $k$ を固定し、$s,c,d$ で $u$ での $\operatorname{sn},\operatorname{cn},\operatorname{dn}$ を、$S,C,D$ で $v$ での値を表す。$\Delta=1-k^2s^2S^2>0$ と置けば、$u+v$ での三つの値は
$$
 \widetilde s=\frac{sCD+Scd}{\Delta},\qquad
 \widetilde c=\frac{cC-sSdD}{\Delta},\qquad
 \widetilde d=\frac{dD-k^2sScC}{\Delta}
$$
である。これを示すため、まず $c^2=1-s^2$、$d^2=1-k^2s^2$ と $S,C,D$ に対する同じ関係を使う。分子をそれぞれ $N_s,N_c,N_d$ とすると、二乗の交差項が相殺し、
$$
\begin{aligned}
N_s^2+N_c^2&=(s^2D^2+c^2)(C^2+S^2d^2)=\Delta^2,\\
N_d^2+k^2N_s^2&=(d^2+k^2s^2 C^2)(D^2+k^2 S^2 c^2)=\Delta^2
\end{aligned}
$$
となる。さらに $u$ で微分すると $\Delta'=-2k^2sS^2cd$、$N_s'=cdCD-sS(d^2+k^2c^2)$ である。$N_s'\Delta-N_s\Delta'$ を展開して $c^2,d^2,C^2,D^2$ を消去すると $N_cN_d$ になる。従って $\widetilde s'=\widetilde c\widetilde d$ である。上の二乗の関係を微分して $\widetilde c'=-\widetilde s\widetilde d$、$\widetilde d'=-k^2\widetilde s\widetilde c$ を得る。$\widetilde c$ の零点にも等式は連続に延びる。実際、$\widetilde c$ は解析的であり、$u=0$ で値 $C$ が零の場合でも導関数 $-SD$ は非零である。$u=0$ の三つの値は $(S,C,D)$ なので、一次系の初期値の一意性により $u+v$ の値と一致する。

例えば $v=u$ とすると $\operatorname{sn}(2u)=2scd/(1-k^2s^4)$ である。$k=0$ を入れれば三式は $\sin,\cos,1$ の加法公式になり、法を $1$ に近づけると差が現れる。実際、$k=1$ の積分を $|\phi|<\pi/2$ で行うと $u=\operatorname{artanh}(\sin\phi)$ だから、その逆関数は $(\tanh u,\operatorname{sech}u,\operatorname{sech}u)$ である。この極限では元の有限周期が失われる。$K(k)$ の被積分関数は $k\uparrow1$ で $\theta=\pi/2$ の近くの $1/\cos\theta$ に単調収束するため、$K(k)$ は無限大に発散する。

### 完全楕円積分の小振幅計算

二項級数を $|k|$ が $1$ 未満の閉区間に収まる範囲で展開すると、角度について一様収束する。各項の $\sin^{2m}\theta$ の積分にも中央二項係数が出る。振子の周期では $k=\sin(\theta_0/2)$ を入れるため、法の冪を得てから開き角に戻す。

:::example[完全楕円積分の低次係数（1）]
変数を $x=k^2$ とする。

$$
m=1
$$

積分前の二項級数の係数を出す。

$$
\begin{aligned}
c_m\big|_{m=1} &= \frac{1}{2}
\end{aligned}
$$

正弦の偶数冪を端点まで積分する。

$$
\begin{aligned}
\int_0^{\pi/2}\sin^{2m}\theta\,\dd\theta\big|_{m=1} &= \frac{\pi}{4}
\end{aligned}
$$

二つの値を掛けて係数を得る。

$$
\begin{aligned}
[x^m]K(\sqrt{x})\big|_{m=1} &= \frac{\pi}{8}
\end{aligned}
$$
:::

:::example[完全楕円積分の低次係数（2）]
変数を $x=k^2$ とする。次数を一つ上げる。

$$
m=2
$$

積分前の二項級数の係数を出す。

$$
\begin{aligned}
c_m\big|_{m=2} &= \frac{3}{8}
\end{aligned}
$$

正弦の偶数冪を端点まで積分する。

$$
\begin{aligned}
\int_0^{\pi/2}\sin^{2m}\theta\,\dd\theta\big|_{m=2} &= \frac{3 \pi}{16}
\end{aligned}
$$

二つの値を掛けて係数を得る。

$$
\begin{aligned}
[x^m]K(\sqrt{x})\big|_{m=2} &= \frac{9 \pi}{128}
\end{aligned}
$$
:::

:::example[完全楕円積分の低次係数（3）]
変数を $x=k^2$ とする。さらに次数を上げて、低次の結果と照合する。

$$
m=3
$$

積分前の二項級数の係数を出す。

$$
\begin{aligned}
c_m\big|_{m=3} &= \frac{5}{16}
\end{aligned}
$$

正弦の偶数冪を端点まで積分する。

$$
\begin{aligned}
\int_0^{\pi/2}\sin^{2m}\theta\,\dd\theta\big|_{m=3} &= \frac{5 \pi}{32}
\end{aligned}
$$

二つの値を掛けて係数を得る。

$$
\begin{aligned}
[x^m]K(\sqrt{x})\big|_{m=3} &= \frac{25 \pi}{512}
\end{aligned}
$$
:::

:::example[完全楕円積分の低次係数（4）]
変数を $x=k^2$ とする。四次の結果まで具体化する。

$$
m=4
$$

積分前の二項級数の係数を出す。

$$
\begin{aligned}
c_m\big|_{m=4} &= \frac{35}{128}
\end{aligned}
$$

正弦の偶数冪を端点まで積分する。

$$
\begin{aligned}
\int_0^{\pi/2}\sin^{2m}\theta\,\dd\theta\big|_{m=4} &= \frac{35 \pi}{256}
\end{aligned}
$$

二つの値を掛けて係数を得る。

$$
\begin{aligned}
[x^m]K(\sqrt{x})\big|_{m=4} &= \frac{1225 \pi}{32768}
\end{aligned}
$$
:::

## Mathieu 関数

### 周期係数が周波数を結び付ける

:term[マチュー方程式]{en="Mathieu equation" yomi="まちゅーほうていしき"} は $y''+(a-2q\cos2x)y=0$ である。$2\pi$ 周期解を $y=\sum c_ne^{inx}$ とすれば $(a-n^2)c_n=q(c_{n-2}+c_{n+2})$ となる。偶奇の添字は混ざらず、$c_{-n}=\pm c_n$ でも分かれる。$q=0$ では $a=n^2$ が周期固有値である。$n=1$ の二重固有値では、$\cos x$ と $\sin x$ を正規化して $2q\cos2x$ を射影すると対角が $q,-q$、非対角が零になる。実際 $2\cos2x\cos x=\cos3x+\cos x$、$2\cos2x\sin x=\sin3x-\sin x$ である。従って $a_1(q)=1+q+O(q^2)$、$b_1(q)=1-q+O(q^2)$ である。各反射対称部分空間で固有値は孤立で単純であり、補空間の逆作用素は残りの固有値との差で有界である。そこへ級数を反復すると固有ベクトルの修正が $O(q)$、固有値の残りが $O(q^2)$ となる。一般の解の周期から周期への増幅を調べるには基本解行列を長さ $\pi$ だけ進めた行列を作る。Wronskian が一定なのでその行列式は $1$ で、実 $a,q$ でトレースの絶対値が $2$ を超えると一つの乗数の絶対値が $1$ を超える。

### 周波数の隣接結合

Mathieu 方程式の作用素 $-\partial_x^2+2q\cos 2x$ は、指数基底で次数が $\pm2$ ずれた成分を結ぶ。このため偶数と奇数は混ざらない。余弦基底の $m=1$ では負側の隣接項が $\cos(-x)=\cos x$ となるので、作用のうち一次の項が元のモードにも残る。

:::example[余弦成分から出る側帯（1）]
余弦基底に作用素を掛ける。

$$
m=1
$$

二階微分では同じモードに次数の二乗が掛かる。

$$
\begin{aligned}
-\partial_x^2\cos(mx)\big|_{m=1} &= \cos{\left(x \right)}
\end{aligned}
$$

ポテンシャルの積を積和公式で二つにする。

$$
\begin{aligned}
2q\cos(2x)\cos(mx)\big|_{m=1} &= q \left(\cos{\left(x \right)} + \cos{\left(3 x \right)}\right)
\end{aligned}
$$

二つの作用を足す。

$$
\begin{aligned}
L\cos(mx)\big|_{m=1} &= q \left(\cos{\left(x \right)} + \cos{\left(3 x \right)}\right) + \cos{\left(x \right)}
\end{aligned}
$$

$x=\frac{\pi}{6}$ で値と導関数を評価する。

$$
\begin{aligned}
-\partial_x^2\cos(mx)\big|_{x=\frac{\pi}{6}} &= \frac{\sqrt{3}}{2} \\
\frac{\partial^1 (-\partial_x^2\cos(mx))}{\partial x^1} &= - \sin{\left(x \right)} \\
\left.\frac{\partial^1 (-\partial_x^2\cos(mx))}{\partial x^1}\right|_{x=\frac{\pi}{6}} &= - \frac{1}{2} \\
\frac{\partial^2 (-\partial_x^2\cos(mx))}{\partial x^2} &= - \cos{\left(x \right)} \\
\left.\frac{\partial^2 (-\partial_x^2\cos(mx))}{\partial x^2}\right|_{x=\frac{\pi}{6}} &= - \frac{\sqrt{3}}{2}
\end{aligned}
$$

$x=\frac{\pi}{6}$ で値と導関数を評価する。

$$
\begin{aligned}
2q\cos(2x)\cos(mx)\big|_{x=\frac{\pi}{6}} &= \frac{\sqrt{3} q}{2} \\
\frac{\partial^1 (2q\cos(2x)\cos(mx))}{\partial x^1} &= - q \left(\sin{\left(x \right)} + 3 \sin{\left(3 x \right)}\right) \\
\left.\frac{\partial^1 (2q\cos(2x)\cos(mx))}{\partial x^1}\right|_{x=\frac{\pi}{6}} &= - \frac{7 q}{2} \\
\frac{\partial^2 (2q\cos(2x)\cos(mx))}{\partial x^2} &= - q \left(\cos{\left(x \right)} + 9 \cos{\left(3 x \right)}\right) \\
\left.\frac{\partial^2 (2q\cos(2x)\cos(mx))}{\partial x^2}\right|_{x=\frac{\pi}{6}} &= - \frac{\sqrt{3} q}{2}
\end{aligned}
$$

$x=\frac{\pi}{6}$ で値と導関数を評価する。

$$
\begin{aligned}
L\cos(mx)\big|_{x=\frac{\pi}{6}} &= \frac{\sqrt{3} \left(q + 1\right)}{2} \\
\frac{\partial^1 (L\cos(mx))}{\partial x^1} &= - q \sin{\left(x \right)} - 3 q \sin{\left(3 x \right)} - \sin{\left(x \right)} \\
\left.\frac{\partial^1 (L\cos(mx))}{\partial x^1}\right|_{x=\frac{\pi}{6}} &= - \frac{7 q}{2} - \frac{1}{2} \\
\frac{\partial^2 (L\cos(mx))}{\partial x^2} &= - q \cos{\left(x \right)} - 9 q \cos{\left(3 x \right)} - \cos{\left(x \right)} \\
\left.\frac{\partial^2 (L\cos(mx))}{\partial x^2}\right|_{x=\frac{\pi}{6}} &= \frac{\sqrt{3} \left(- q - 1\right)}{2}
\end{aligned}
$$
:::

:::example[余弦成分から出る側帯（2）]
余弦基底に作用素を掛ける。次数を一つ上げる。

$$
m=2
$$

二階微分では同じモードに次数の二乗が掛かる。

$$
\begin{aligned}
-\partial_x^2\cos(mx)\big|_{m=2} &= 4 \cos{\left(2 x \right)}
\end{aligned}
$$

ポテンシャルの積を積和公式で二つにする。

$$
\begin{aligned}
2q\cos(2x)\cos(mx)\big|_{m=2} &= q \left(\cos{\left(4 x \right)} + 1\right)
\end{aligned}
$$

二つの作用を足す。

$$
\begin{aligned}
L\cos(mx)\big|_{m=2} &= q \cos{\left(4 x \right)} + q + 4 \cos{\left(2 x \right)}
\end{aligned}
$$

$x=\frac{\pi}{6}$ で値と導関数を評価する。

$$
\begin{aligned}
-\partial_x^2\cos(mx)\big|_{x=\frac{\pi}{6}} &= 2 \\
\frac{\partial^1 (-\partial_x^2\cos(mx))}{\partial x^1} &= - 8 \sin{\left(2 x \right)} \\
\left.\frac{\partial^1 (-\partial_x^2\cos(mx))}{\partial x^1}\right|_{x=\frac{\pi}{6}} &= - 4 \sqrt{3} \\
\frac{\partial^2 (-\partial_x^2\cos(mx))}{\partial x^2} &= - 16 \cos{\left(2 x \right)} \\
\left.\frac{\partial^2 (-\partial_x^2\cos(mx))}{\partial x^2}\right|_{x=\frac{\pi}{6}} &= -8
\end{aligned}
$$

$x=\frac{\pi}{6}$ で値と導関数を評価する。

$$
\begin{aligned}
2q\cos(2x)\cos(mx)\big|_{x=\frac{\pi}{6}} &= \frac{q}{2} \\
\frac{\partial^1 (2q\cos(2x)\cos(mx))}{\partial x^1} &= - 4 q \sin{\left(4 x \right)} \\
\left.\frac{\partial^1 (2q\cos(2x)\cos(mx))}{\partial x^1}\right|_{x=\frac{\pi}{6}} &= - 2 \sqrt{3} q \\
\frac{\partial^2 (2q\cos(2x)\cos(mx))}{\partial x^2} &= - 16 q \cos{\left(4 x \right)} \\
\left.\frac{\partial^2 (2q\cos(2x)\cos(mx))}{\partial x^2}\right|_{x=\frac{\pi}{6}} &= 8 q
\end{aligned}
$$

$x=\frac{\pi}{6}$ で値と導関数を評価する。

$$
\begin{aligned}
L\cos(mx)\big|_{x=\frac{\pi}{6}} &= \frac{q}{2} + 2 \\
\frac{\partial^1 (L\cos(mx))}{\partial x^1} &= - 4 q \sin{\left(4 x \right)} - 8 \sin{\left(2 x \right)} \\
\left.\frac{\partial^1 (L\cos(mx))}{\partial x^1}\right|_{x=\frac{\pi}{6}} &= 2 \sqrt{3} \left(- q - 2\right) \\
\frac{\partial^2 (L\cos(mx))}{\partial x^2} &= - 16 q \cos{\left(4 x \right)} - 16 \cos{\left(2 x \right)} \\
\left.\frac{\partial^2 (L\cos(mx))}{\partial x^2}\right|_{x=\frac{\pi}{6}} &= 8 q - 8
\end{aligned}
$$
:::

:::example[余弦成分から出る側帯（3）]
余弦基底に作用素を掛ける。さらに次数を上げて、低次の結果と照合する。

$$
m=3
$$

二階微分では同じモードに次数の二乗が掛かる。

$$
\begin{aligned}
-\partial_x^2\cos(mx)\big|_{m=3} &= 9 \cos{\left(3 x \right)}
\end{aligned}
$$

ポテンシャルの積を積和公式で二つにする。

$$
\begin{aligned}
2q\cos(2x)\cos(mx)\big|_{m=3} &= q \left(\cos{\left(x \right)} + \cos{\left(5 x \right)}\right)
\end{aligned}
$$

二つの作用を足す。

$$
\begin{aligned}
L\cos(mx)\big|_{m=3} &= q \left(\cos{\left(x \right)} + \cos{\left(5 x \right)}\right) + 9 \cos{\left(3 x \right)}
\end{aligned}
$$

$x=\frac{\pi}{6}$ で値と導関数を評価する。

$$
\begin{aligned}
-\partial_x^2\cos(mx)\big|_{x=\frac{\pi}{6}} &= 0 \\
\frac{\partial^1 (-\partial_x^2\cos(mx))}{\partial x^1} &= - 27 \sin{\left(3 x \right)} \\
\left.\frac{\partial^1 (-\partial_x^2\cos(mx))}{\partial x^1}\right|_{x=\frac{\pi}{6}} &= -27 \\
\frac{\partial^2 (-\partial_x^2\cos(mx))}{\partial x^2} &= - 81 \cos{\left(3 x \right)} \\
\left.\frac{\partial^2 (-\partial_x^2\cos(mx))}{\partial x^2}\right|_{x=\frac{\pi}{6}} &= 0
\end{aligned}
$$

$x=\frac{\pi}{6}$ で値と導関数を評価する。

$$
\begin{aligned}
2q\cos(2x)\cos(mx)\big|_{x=\frac{\pi}{6}} &= 0 \\
\frac{\partial^1 (2q\cos(2x)\cos(mx))}{\partial x^1} &= - q \left(\sin{\left(x \right)} + 5 \sin{\left(5 x \right)}\right) \\
\left.\frac{\partial^1 (2q\cos(2x)\cos(mx))}{\partial x^1}\right|_{x=\frac{\pi}{6}} &= - 3 q \\
\frac{\partial^2 (2q\cos(2x)\cos(mx))}{\partial x^2} &= - q \left(\cos{\left(x \right)} + 25 \cos{\left(5 x \right)}\right) \\
\left.\frac{\partial^2 (2q\cos(2x)\cos(mx))}{\partial x^2}\right|_{x=\frac{\pi}{6}} &= 12 \sqrt{3} q
\end{aligned}
$$

$x=\frac{\pi}{6}$ で値と導関数を評価する。

$$
\begin{aligned}
L\cos(mx)\big|_{x=\frac{\pi}{6}} &= 0 \\
\frac{\partial^1 (L\cos(mx))}{\partial x^1} &= - q \sin{\left(x \right)} - 5 q \sin{\left(5 x \right)} - 27 \sin{\left(3 x \right)} \\
\left.\frac{\partial^1 (L\cos(mx))}{\partial x^1}\right|_{x=\frac{\pi}{6}} &= - 3 q - 27 \\
\frac{\partial^2 (L\cos(mx))}{\partial x^2} &= - q \cos{\left(x \right)} - 25 q \cos{\left(5 x \right)} - 81 \cos{\left(3 x \right)} \\
\left.\frac{\partial^2 (L\cos(mx))}{\partial x^2}\right|_{x=\frac{\pi}{6}} &= 12 \sqrt{3} q
\end{aligned}
$$
:::

:::example[余弦成分から出る側帯（4）]
余弦基底に作用素を掛ける。四次の結果まで具体化する。

$$
m=4
$$

二階微分では同じモードに次数の二乗が掛かる。

$$
\begin{aligned}
-\partial_x^2\cos(mx)\big|_{m=4} &= 16 \cos{\left(4 x \right)}
\end{aligned}
$$

ポテンシャルの積を積和公式で二つにする。

$$
\begin{aligned}
2q\cos(2x)\cos(mx)\big|_{m=4} &= q \left(\cos{\left(2 x \right)} + \cos{\left(6 x \right)}\right)
\end{aligned}
$$

二つの作用を足す。

$$
\begin{aligned}
L\cos(mx)\big|_{m=4} &= q \left(\cos{\left(2 x \right)} + \cos{\left(6 x \right)}\right) + 16 \cos{\left(4 x \right)}
\end{aligned}
$$

$x=\frac{\pi}{6}$ で値と導関数を評価する。

$$
\begin{aligned}
-\partial_x^2\cos(mx)\big|_{x=\frac{\pi}{6}} &= -8 \\
\frac{\partial^1 (-\partial_x^2\cos(mx))}{\partial x^1} &= - 64 \sin{\left(4 x \right)} \\
\left.\frac{\partial^1 (-\partial_x^2\cos(mx))}{\partial x^1}\right|_{x=\frac{\pi}{6}} &= - 32 \sqrt{3} \\
\frac{\partial^2 (-\partial_x^2\cos(mx))}{\partial x^2} &= - 256 \cos{\left(4 x \right)} \\
\left.\frac{\partial^2 (-\partial_x^2\cos(mx))}{\partial x^2}\right|_{x=\frac{\pi}{6}} &= 128
\end{aligned}
$$

$x=\frac{\pi}{6}$ で値と導関数を評価する。

$$
\begin{aligned}
2q\cos(2x)\cos(mx)\big|_{x=\frac{\pi}{6}} &= - \frac{q}{2} \\
\frac{\partial^1 (2q\cos(2x)\cos(mx))}{\partial x^1} &= - 2 q \left(\sin{\left(2 x \right)} + 3 \sin{\left(6 x \right)}\right) \\
\left.\frac{\partial^1 (2q\cos(2x)\cos(mx))}{\partial x^1}\right|_{x=\frac{\pi}{6}} &= - \sqrt{3} q \\
\frac{\partial^2 (2q\cos(2x)\cos(mx))}{\partial x^2} &= - 4 q \left(\cos{\left(2 x \right)} + 9 \cos{\left(6 x \right)}\right) \\
\left.\frac{\partial^2 (2q\cos(2x)\cos(mx))}{\partial x^2}\right|_{x=\frac{\pi}{6}} &= 34 q
\end{aligned}
$$

$x=\frac{\pi}{6}$ で値と導関数を評価する。

$$
\begin{aligned}
L\cos(mx)\big|_{x=\frac{\pi}{6}} &= - \frac{q}{2} - 8 \\
\frac{\partial^1 (L\cos(mx))}{\partial x^1} &= - 2 q \left(\sin{\left(2 x \right)} + 3 \sin{\left(6 x \right)}\right) - 64 \sin{\left(4 x \right)} \\
\left.\frac{\partial^1 (L\cos(mx))}{\partial x^1}\right|_{x=\frac{\pi}{6}} &= \sqrt{3} \left(- q - 32\right) \\
\frac{\partial^2 (L\cos(mx))}{\partial x^2} &= - 4 q \left(\cos{\left(2 x \right)} + 9 \cos{\left(6 x \right)}\right) - 256 \cos{\left(4 x \right)} \\
\left.\frac{\partial^2 (L\cos(mx))}{\partial x^2}\right|_{x=\frac{\pi}{6}} &= 34 q + 128
\end{aligned}
$$
:::

### Mathieu の最低特性値のずれ

周期 $2\pi$ で偶関数の空間に制限し、$-u''+2q\cos(2x)u=a(q)u$ の $q=0$ での最低固有関数を定数 $1$ に選ぶ。平均を常に $1$ に固定し $u=1+qv+q^2w+\cdots$、$a=a_1q+a_2q^2+\cdots$ とおく。一階は $-v''+2\cos 2x=a_1$ である。周期積分で $a_1=0$、平均零で $v=-\cos 2x/2$ になる。二階を周期積分すると $a_2$ は $2\cos(2x)v$ の平均であり、$a_2=-1/2$ である。従って最低特性値は $a(q)=-q^2/2+O(q^3)$ となる。

この展開を正当化するには、平均零の偶周期関数上で $-\partial_x^2$ の逆が有界であることを使う。Fourier 係数ごとに $1/n^2$ が掛かるため、平均零の $L^2$ から二階微分が $L^2$ の空間へ有界に作用する。方程式の平均と平均零への射影を別々に取り、未知の平均零部分と $a$ で微分すると、$q=0$ の導関数はこの微分作用素と符号付きのスカラー恒等写像からなる可逆なブロック行列である。陰関数定理により十分小さい実 $q$ で滑らかな解枝が存在し、二階までの剰余が得られる。他の固有値が近づかないことは、平均零の偶周期関数に $\int |u'|^2\ge\int |u|^2$ が成り立つことから示せる。平均が零であれば、互いに直交する Fourier モードの波数が $1$ 以上だからである。二つの独立な固有関数がともに $a<1-2|q|$ に属すると、その線形結合に非零の平均零関数を取れる。しかしこの部分空間での二次形式は $(1-2|q|)\|u\|^2$ 以上であり、直交した二固有関数の固有値の重み付き平均と矛盾する。

### Floquet 乗数と周期解

実パラメータ $a,q$ に対して $(y,y')$ を一周期 $\pi$ だけ送る実行列を $T$、そのトレースを $\Delta$ とする。周期性と初期値の一意性から、$T$ の固有ベクトルに対応する解は $y(x+\pi)=\rho y(x)$ を満たす。Abel の公式 :eqref[eq-mp-abel-wronskian] により $\det T=1$ で、その固有値は $\rho^2-\Delta\rho+1=0$ の二根である。$|\Delta|<2$ なら相異なる共役な単位円上の根なので全解が有界である。$|\Delta|>2$ なら相異なる実根で、一方は絶対値が $1$ を超える。従ってその枝は周期を追うと指数的に増大する。$\Delta=\pm2$ では $T=\pm I$ なら全解が有界で、そうでなければ Jordan 形から一方の枝に周期数に比例する増大が現れる。それぞれのケースで周期の途中までの発展は有界な行列を掛けるだけなので、周期ごとの判定を連続な $x$ へ移せる。

$\rho=1$ の枝は $\pi$ 周期、$\rho=-1$ の枝は半周期で符号が変わる $2\pi$ 周期である。反射 $x\mapsto-x$ は同じ方程式を保つので、周期解は偶解と奇解へ分けられる。偶解を $\mathrm{ce}_n(x,q)$、奇解を $\mathrm{se}_n(x,q)$ と書き、$q=0$ に連続な符号を取り、区間 $[0,2\pi]$ で二乗積分を $\pi$ に規格化する。$q=0$ の $n\ge1$ では余弦と正弦、$n=0$ の偶解では $1/\sqrt2$ に一致する。対応する特性値を $a_n(q),b_n(q)$ と書く。例えば $q=0$ の $\Delta$ は $a>0$ で $2\cos(\pi\sqrt a)$、$a<0$ で $2\cosh(\pi\sqrt{-a})$ である。$a=0$ では $1,x$ が解なので $T$ の対角は $1$、右上は $\pi$ となり、重複乗数 $1$ があっても全解は有界ではない。

## 演習問題

::::exercise{level=basic}
${}_2F_1(-2,3;1;z)$ を求めよ。
:::solution
三次以降は $(-2)_n=0$ で、係数は $1,-6,6$ なので $1-6z+6z^2$ である。
:::
::::

::::exercise{level=basic}
$M(-1,3,z)$ と $U(1,2,z)$ を求めよ。
:::solution
$M$ は $1-z/3$ で級数が止まる。$U$ の積分では $(1+t)^0=1$、$\Gamma(1)=1$ なので $1/z$ である。
:::
::::

::::exercise{level=basic}
Mathieu 方程式で $q=0$ の周期解 $y=\cos3x$ の $a$ を求めよ。
:::solution
$y''=-9y$ を代入すると $a=9$ が必要で、実際に周期解である。
:::
::::

::::exercise{level=basic}
$D_0(x)$ の方程式と減衰規格化を確認せよ。
:::solution
$D_0=e^{-x^2/4}$ の二階微分は $(x^2/4-1/2)D_0$ で、規格化係数は $1$ である。
:::
::::

::::exercise{level=advanced}
$F(1,1;2;z)$ の $z=1$ 近くの対数の係数を求めよ。
:::solution
$F=-\log(1-z)/z$ である。$s=1-z$ とすると $F=-\log s+O(s|\log s|)$ なので先頭の係数は $-1$ である。
:::
::::

::::exercise{level=advanced}
$\operatorname{sn}(u,k)$ が $s''=-(1+k^2)s+2k^2s^3$ を満たすことを示せ。
:::solution
$s'=\operatorname{cn}\operatorname{dn}$ を微分すると $-s\operatorname{dn}^2-k^2s\operatorname{cn}^2$ である。$\operatorname{dn}^2=1-k^2s^2$ と $\operatorname{cn}^2=1-s^2$ を入れる。
:::
::::

::::exercise{level=advanced}
$K(k)$ の $k^4$ の係数を求めよ。
:::solution
二項級数の第 2 項の後の係数は $3/8$ であり、$\int_0^{\pi/2}\sin^4\theta\dd\theta=3\pi/16$ なので $9\pi/128$ である。
:::
::::

::::exercise{level=advanced}
放物柱方程式で $\nu=1$ の正規化した減衰解を直接求めよ。
:::solution
$y=xe^{-x^2/4}$ の二階微分は $(x^2/4-3/2)y$ である。無限遠で $xe^{-x^2/4}$ に一致するので $D_1$ である。
:::
::::

## 確認テスト

::::quiz{answer="-6" tol="1e-9"}
${}_2F_1(-2,3;1;z)$ の一次係数を答えよ。

:::explanation
$n=1$ の係数は $(-2)\cdot3/(1\cdot1)=-6$ である。
:::
::::

::::quiz
Gauss の方程式の三つの特異点の候補はどれか。

- [x] $0,1,\infty$
- [ ] $0,-1,\infty$
- [ ] $0,1,2$

:::explanation
有限点では最高階係数 $z(1-z)$ の零点が候補であり、無限遠は $1/z$ を局所座標に調べる。
:::
::::

::::quiz{answer="1" tol="1e-9"}
$\operatorname{dn}(0,k)$ の値を求めよ。

:::explanation
$u=0$ なら逆振幅は $\phi=0$ だから、$\sqrt{1-k^2\sin^2\phi}=1$ である。
:::
::::

::::quiz
Kummer 方程式の無限遠の性質はどれか。

- [ ] 常に正則点
- [ ] 確定特異点
- [x] 不確定特異点

:::explanation
正規形の一階係数は $c/z-1$ で、無限遠の確定特異性に必要な $O(1/z)$ にならない。
:::
::::

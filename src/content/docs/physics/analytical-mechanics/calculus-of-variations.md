---
title: 変分法
description: 汎関数の第一変分から Euler–Lagrange 方程式を導き、正則性・自然境界条件・第二変分と十分条件を調べ、最速降下線・懸垂線・等周問題・各点での拘束を伴う変分問題を解く。
notation:
  - symbol: 'J[y]'
    meaning: 関数 y に実数を対応させる汎関数
  - symbol: '\norm{y}_0,\ \norm{y}_1'
    meaning: 関数の最大値ノルムと C¹ ノルム
  - symbol: '\eta'
    meaning: 曲線を変える方向を指定する関数
  - symbol: '\delta J[y;\eta]'
    meaning: 関数 y における方向 η への第一変分
  - symbol: '\delta^2 J[y;\eta]'
    meaning: 関数 y における方向 η への第二変分
  - symbol: 'F_p'
    meaning: 被積分関数 F の傾き引数 p に関する偏導関数
  - symbol: '\lambda'
    meaning: 拘束を課すための Lagrange 乗数
---

## この章の目標

- 関数の近さを測るノルムを導入して汎関数の極値を定義し、端点を固定した変形に対する第一変分を計算する。
- 変分法の基本補題と du Bois-Reymond の補題を証明し、Euler–Lagrange 方程式、その解の正則性、自然境界条件、高階導関数を含む場合を導く。
- 第二変分と凸性を使って、停留曲線が実際に極小を与えるかを判定する。
- 最速降下線を求めてその最速性を証明し、懸垂面の解の個数を調べる。
- 積分で課す拘束と各点で課す拘束に対して Lagrange 乗数の存在を証明し、等周問題、鎖の形、球面上の測地線に応用する。

重力の下で、二点を結ぶどの曲線に沿って滑れば最も早く着くか。
短い経路では移動距離を節約できるが、初めに深く降りる経路では早く加速できる。
この二つの効果を比較するには、曲線を少しずつ変えたときの所要時間を計算したい。
そのための微分法が :term[変分法]{en="calculus of variations" yomi="へんぶんほう"} である。
本章の内容は数学としての変分法の基礎であり、次章 [Lagrange 形式](/physics/analytical-mechanics/lagrangian/) で力学の運動方程式を作用の停留条件として書き直すために使う。
数学的な扱いは :cite[courant-hilbert]、力学への応用は :cite[arnold-mmcm] にも述べられている。

## 汎関数と変分

### 関数を入力とする量

:::definition[汎関数と許容関数]{#def-am-functional}
指定された条件を満たす関数全体を $\mathcal A$ とする。$y\in\mathcal A$ に数 $J[y]$ を対応させる写像を
:term[汎関数]{en="functional" yomi="はんかんすう"} といい、$\mathcal A$ の元を**許容関数**と呼ぶ。
本章の基本例は、定数 $a<b$、指定した端点値 $A,B$ に対して

$$
J[y]=\int_a^b F(x,y(x),y'(x))\,\dd x,\qquad y(a)=A,\quad y(b)=B
\label{eq-am-functional}
$$

と書けるものである。プライムは独立変数 $x$ に関する微分を表す。特に断らない限り、許容関数は $[a,b]$ 上の $C^1$ 級関数で端点条件を満たすものとする。
:::

例えば平面上のグラフの長さは $F(x,y,p)=\sqrt{1+p^2}$ で得られる。
積分の中の $y'(x)$ は曲線に沿って決まるが、$F$ 自体を偏微分するときには、その三つの引数 $(x,y,p)$ は独立に扱う。
以下、$F_y=\partial F/\partial y$、$F_p=\partial F/\partial p$、$F_{pp}=\partial^2F/\partial p^2$ などと書き、式を曲線上で評価するときに $p=y'(x)$ を入れる。

本章で扱う汎関数の例を挙げておく。

1. 曲線の長さ $\int_a^b\sqrt{1+(y')^2}\,\dd x$。
2. 重力の下での降下時間 $\int_0^X\sqrt{1+(y')^2}/\sqrt{2gy}\,\dd x$（最速降下線）。
3. 回転面の面積 $2\pi\int_a^b y\sqrt{1+(y')^2}\,\dd x$（懸垂面）。
4. 力学系の作用 $\int_{t_0}^{t_1}L(q,\dot q,t)\,\dd t$（次章）。

### 関数の近さと極値

関数を少し動かしたときに汎関数がどう変わるかを調べるには、関数どうしの近さを決めておく必要がある。
$[a,b]$ 上の連続関数 $u$ と $C^1$ 級関数 $u$ に対し

$$
\norm{u}_0=\max_{a\le x\le b}\abs{u(x)},\qquad
\norm{u}_1=\max_{a\le x\le b}\abs{u(x)}+\max_{a\le x\le b}\abs{u'(x)}
$$

とおく。前者を最大値ノルム、後者を $C^1$ ノルムという。$\norm{u}_1$ が小さければ $u$ も $u'$ も小さい。

:::definition[強い極値と弱い極値]{#def-am-local-extremum}
許容関数 $y$ について、正の数 $r$ があって、$\norm{\tilde y-y}_0<r$ を満たすすべての許容関数 $\tilde y$ に対し $J[\tilde y]\ge J[y]$ となるとき、$y$ は $J$ の**強い極小**を与えるという。
$\norm{\tilde y-y}_1<r$ を満たすすべての許容関数について同じ不等式が成り立つとき、**弱い極小**を与えるという。極大も同様に定める。
:::

$\norm{\tilde y-y}_1<r$ ならば $\norm{\tilde y-y}_0<r$ なので、強い極小は弱い極小である。
逆は成り立たない。傾きは大きく違うが値は近い曲線と比べるかどうかが両者の差である（演習問題を参照）。

関数の変数が無限個あることから、有限次元にはない現象も起こる。最小値が存在しないことがあるのはその一つである。

:::example[最小値をとらない汎関数]{#exm-am-weierstrass}
$J[y]=\int_{-1}^1x^2(y')^2\,\dd x$ を、$y(-1)=-1$、$y(1)=1$ を満たす $C^1$ 級関数の上で考える。
$J\ge0$ である。$\eps>0$ に対し $c_\eps=\arctan(1/\eps)$、$y_\eps(x)=\arctan(x/\eps)/c_\eps$ とおくと端点条件を満たし、$y_\eps'=\eps/(c_\eps(\eps^2+x^2))$ である。$x^2\le\eps^2+x^2$ を使うと

$$
J[y_\eps]\le\frac{1}{c_\eps^2}\int_{-1}^1\frac{\eps^2}{\eps^2+x^2}\,\dd x=\frac{2\eps}{c_\eps}\to0\qquad(\eps\to0)
$$

となる。したがって下限は $0$ である。一方 $J[y]=0$ なら $x\ne0$ で $y'=0$ であり、$y$ は $[-1,0)$ と $(0,1]$ でそれぞれ定数 $-1$、$1$ となって $x=0$ で連続でない。よって下限 $0$ をとる許容関数は存在しない。
:::

このため変分問題では、停留曲線を求めることと、それが実際に最小を与えることを区別して扱う。

### 端点を固定した曲線の変形

比較したい曲線を $y$ とし、$y_\eps(x)=y(x)+\eps\eta(x)$ という一変数の族を作る。
端点を固定するため、$\eta\in C^1([a,b])$ は $\eta(a)=\eta(b)=0$ を満たすとする。この条件の下で $\eta$ の形は自由に選べる。このような $\eta$ を**許容方向**と呼ぶ。
前章の仮想変位（:ref[def-am-virtual-displacement]）が一つの時刻での配位を動かしたのに対して、ここでは区間全体の曲線を一度に変える。

:::definition[第一変分と停留]{#def-am-first-variation}
$\eps$ を実数として $J[y_\eps]$ を $\eps$ で微分した値

$$
\delta J[y;\eta]=\left.\dv{}{\eps}J[y+\eps\eta]\right|_{\eps=0}
$$

を $\eta$ 方向の :term[第一変分]{en="first variation" yomi="だいいちへんぶん"} という。
すべての許容方向 $\eta$ に対してこれが零のとき、$y$ は $J$ の :term[停留曲線]{en="stationary curve" yomi="ていりゅうきょくせん"}（:term[極値曲線]{en="extremal" yomi="きょくちきょくせん"}）である。
:::

$y$ が弱い極小を与えるとする。固定した $\eta$ に対して $\norm{\eps\eta}_1=\abs{\eps}\norm{\eta}_1$ なので、$\abs\eps$ が小さければ $y+\eps\eta$ は $y$ の $C^1$ 近傍にある。
したがって一変数関数 $\eps\mapsto J[y+\eps\eta]$ は $\eps=0$ で極小をとり、微分可能ならその微分は零である。
すなわち**弱い極値（したがって強い極値）を与える曲線は停留曲線である**。
逆向きの判断には、第一変分より先の項や、他の曲線との直接比較が必要になる。

:::proposition[積分汎関数の第一変分]{#prop-am-first-variation}
$y,\eta\in C^1([a,b])$ とし、グラフ $\{(x,y(x),y'(x))\}$ を含む開集合で $F$ は $C^1$ 級とする。このとき

$$
\delta J[y;\eta]=\int_a^b \big(F_y\eta+F_p\eta'\big)\,\dd x.
\label{eq-am-first-variation}
$$
:::

:::proof
グラフはコンパクトなので、ある $\rho>0$ について、グラフからの距離が $\rho$ 以下の点全体 $K$ も $F$ が $C^1$ 級である開集合に含まれる。
$\abs\eps\,\norm{\eta}_1\le\rho$ なら $(x,y+\eps\eta,y'+\eps\eta')\in K$ である。$\eps$ に依存する被積分関数 $F(x,y+\eps\eta,y'+\eps\eta')$ を $\eps$ で偏微分すると、連鎖律により、変形後の引数で評価した $F_y\eta+F_p\eta'$ となる。
これは $(x,\eps)$ の連続関数であり、コンパクト集合上で一様連続である。したがって積分記号の下で微分でき（[Riemann 積分](/math/analysis/riemann-integral/)）、$\eps=0$ を代入して主張を得る。
:::

:::example[微分してから積分する]
$J[y]=\tfrac12\int_0^1 (y')^2\,\dd x$、$y(0)=0,\ y(1)=1$ とする。
直線 $y=x$ のまわりでは、平方を展開して

$$
J[x+\eps\eta]=\tfrac12+\eps\int_0^1\eta'\,\dd x
+\tfrac{\eps^2}{2}\int_0^1(\eta')^2\,\dd x
$$

となる。一次の係数は $\eta(1)-\eta(0)=0$ である。さらに二次の項は非負であり、
$\eta$ が恒等的に零でない限り正になる。ここでは、直線が唯一の最小解であることまで結論できる。
:::

## Euler–Lagrange 方程式

### 積分から各点の条件へ

第一変分には $\eta$ と $\eta'$ が現れるが、この二つは独立には選べない。
そこで $\eta'$ の項を部分積分し、一種類の任意関数 $\eta$ だけが残る形に直す。
その積分を各点の方程式に読み替えるのが次の補題である。

:::lemma[変分法の基本補題]{#lem-am-fundamental-variation}
$h\in C([a,b])$ とし、$k\ge1$ を整数とする。$[a,b]$ の内部の閉区間の外で零になるすべての $\eta\in C^k([a,b])$ に対して

$$
\int_a^b h(x)\eta(x)\,\dd x=0
$$

ならば $h$ は区間上で恒等的に零である。特に、両端で零になるすべての $\eta\in C^1([a,b])$ についてこの積分が零なら $h=0$ である。
:::

:::proof
ある内点で $h>0$ と仮定する。連続性により、その点を含む小区間 $[c,d]\subset(a,b)$ と定数 $\kappa>0$ を選んで $[c,d]$ 上で $h\ge\kappa$ とできる。
この小区間内では $\eta(x)=\big((x-c)(d-x)\big)^{k+1}$、外では $\eta(x)=0$ とおく。
$c,d$ では値と $k$ 階までの導関数がすべて零なので、つないだ関数は $C^k$ 級である。しかも

$$
\int_a^b h\eta\,\dd x\ge \kappa\int_c^d\big((x-c)(d-x)\big)^{k+1}\,\dd x>0
$$

となり仮定に反する。$h<0$ の内点がある場合も $-h$ に同じ議論を適用できる。
よって内点では $h=0$ であり、両端での値も連続性から零になる。
最後の主張は、この $\eta$ が両端で零になる $C^1$ 級関数であることから従う。
:::

$y$ が $C^2$ 級でなければ $F_p$ を $x$ で微分できるとは限らず、上の部分積分が使えない。その場合には、逆に $\eta$ の項を部分積分して $\eta'$ だけを残す。

:::lemma[du Bois-Reymond の補題]{#lem-am-du-bois-reymond}
$g,h\in C([a,b])$ とする。両端で零になるすべての $\eta\in C^1([a,b])$ に対して

$$
\int_a^b\big(g(x)\eta(x)+h(x)\eta'(x)\big)\,\dd x=0
$$

ならば、$h$ は $C^1$ 級で $h'=g$ である。特に $g=0$ なら $h$ は定数である。
:::

:::proof
まず $g=0$ の場合を示す。$h$ の平均値を $c=\frac{1}{b-a}\int_a^bh\,\dd x$ とし、$\eta(x)=\int_a^x(h(t)-c)\,\dd t$ とおく。
$\eta$ は $C^1$ 級で $\eta(a)=0$、$\eta(b)=\int_a^bh\,\dd t-c(b-a)=0$ を満たす。仮定から

$$
0=\int_a^bh\eta'\,\dd x=\int_a^bh(h-c)\,\dd x=\int_a^b(h-c)^2\,\dd x+c\int_a^b(h-c)\,\dd x=\int_a^b(h-c)^2\,\dd x
$$

である。被積分関数は連続で非負なので $h=c$ が恒等的に成り立つ。

一般の場合、$G(x)=\int_a^xg(t)\,\dd t$ とおけば $G$ は $C^1$ 級で $G'=g$ である。$\eta$ は両端で零なので部分積分により $\int_a^bg\eta\,\dd x=-\int_a^bG\eta'\,\dd x$ となり、仮定は $\int_a^b(h-G)\eta'\,\dd x=0$ と同値である。
前半から $h-G$ は定数であり、$h=G+\text{定数}$ は $C^1$ 級で $h'=g$ となる。
:::

### 一つの未知関数の方程式

:::theorem[Euler–Lagrange 方程式]{#thm-am-euler-lagrange}
固定端の汎関数 :eqref[eq-am-functional] で $F$ は $C^2$ 級とし、$y\in C^1([a,b])$ を許容関数とする。
$y$ が停留曲線であるための必要十分条件は、関数 $x\mapsto F_p(x,y(x),y'(x))$ が $C^1$ 級で、各点で $\dd F_p/\dd x=F_y$ が成り立つことである。
特に $y$ が $C^2$ 級なら、この条件は

$$
\dv{}{x}\pdv{F}{p}-\pdv{F}{y}=0,\qquad p=y'(x)
\label{eq-am-euler-lagrange}
$$

と書ける。これを :term[Euler–Lagrange 方程式]{en="Euler–Lagrange equation" yomi="おいらーらぐらんじゅほうていしき"} という。
:::

:::proof
:ref[prop-am-first-variation] により、停留条件はすべての許容方向 $\eta$ について $\int_a^b(F_y\eta+F_p\eta')\,\dd x=0$ となることである。
ここで $F_y,F_p$ は曲線上で評価した $x$ の連続関数である。:ref[lem-am-du-bois-reymond] により、これから $F_p$ が $C^1$ 級で $\dd F_p/\dd x=F_y$ が従う。
逆にこの条件が成り立てば、部分積分により

$$
\delta J[y;\eta]=[F_p\eta]_a^b+\int_a^b\left(F_y-\dv{}{x}F_p\right)\eta\,\dd x
$$

の境界項は $\eta(a)=\eta(b)=0$ で消え、積分も零である。
:::

全微分 $\dd F_p/\dd x$ には、$F_p$ に陽に現れる $x$ だけでなく $y,y'$ の変化も寄与する。$y$ が $C^2$ 級なら具体的には

$$
\dv{F_p}{x}=F_{px}+F_{py}y'+F_{pp}y''.
$$

したがって $F_{pp}\ne0$ の領域では、方程式を $y''$ について解ける。つまり Euler–Lagrange 方程式は二階の常微分方程式であり、一般解は二つの任意定数を含む。二つの端点条件がそれを決める。

:::proposition[停留曲線の正則性]{#prop-am-regularity}
$F$ は $C^2$ 級、$y\in C^1([a,b])$ は停留曲線とする。点 $x_0$ で $F_{pp}(x_0,y(x_0),y'(x_0))\ne0$ ならば、$y$ は $x_0$ の近くで $C^2$ 級である。
:::

:::proof
:ref[thm-am-euler-lagrange] により $\varphi(x)=F_p(x,y(x),y'(x))$ は $C^1$ 級である。
$\Phi(x,p)=F_p(x,y(x),p)-\varphi(x)$ とおくと、$F_p$ は $C^1$ 級、$y$ と $\varphi$ も $C^1$ 級なので、$\Phi$ は $(x,p)$ について $C^1$ 級である。
$\Phi(x,y'(x))=0$ がすべての $x$ で成り立ち、$\partial\Phi/\partial p=F_{pp}$ は $(x_0,y'(x_0))$ で零でない。
陰関数定理により、$(x_0,y'(x_0))$ のある近傍で $\Phi(x,p)=0$ の解は $p=P(x)$ の形にただ一つ書け、$P$ は $C^1$ 級である。
$y'$ は連続なので $x_0$ の近くでは $(x,y'(x))$ がこの近傍にとどまり、一意性から $y'=P$ となる。よって $y'$ は $C^1$ 級である。
:::

:::remark[角をもつ最小解]
$F_{pp}=0$ となる点では、停留曲線が滑らかとは限らない。
例えば $J[y]=\int_{-1}^1y^2(1-y')^2\,\dd x$、$y(-1)=0$、$y(1)=1$ では $J\ge0$ であり、角をもつ関数 $y=\max(x,0)$ で $J=0$ となる。
しかし $J[y]=0$ を満たす $C^1$ 級関数は存在しない。実際、そのような $y$ では各点で $y=0$ か $y'=1$ である。零点の上限を $x_1$ とすると $(x_1,1]$ で $y'=1$ なので $y=x-x_1$、連続性から $y'(x_1)=1$ である。よって $x_1$ の左側で $y<0$、したがって $y'=1$ が続き、$y=x-x_1$ が $[-1,1]$ 全体で成り立つ。これは二つの端点条件を同時に満たさない。
この例では角の点で $F_{pp}=2y^2=0$ であり、:ref[prop-am-regularity] と矛盾しない。角をもつ曲線まで許容関数を広げる理論（Weierstrass–Erdmann の角条件）は本サイトでは扱わない。
:::

### 複数の関数と自然境界条件

未知関数を $y^1,\ldots,y^n$ とし、それぞれの傾き引数を $p^1,\ldots,p^n$ とする。
すべての $y^i$ を $C^2$ 級とし、全成分の端点を固定して $y^i_\eps=y^i+\eps\eta^i$ とおくと、

$$
\begin{aligned}
\delta J
&=\int_a^b\sum_{i=1}^n
\left(\pdv{F}{y^i}\eta^i+\pdv{F}{p^i}(\eta^i)'\right)\dd x\\
&=\int_a^b\sum_{i=1}^n
\left(\pdv{F}{y^i}-\dv{}{x}\pdv{F}{p^i}\right)\eta^i\,\dd x .
\end{aligned}
$$

成分は独立に変えられるので、第 $i$ 成分だけを変えて他の $\eta^j$ を零にできる。
したがって基本補題から

$$
\dv{}{x}\pdv{F}{p^i}-\pdv{F}{y^i}=0,\qquad
p^i=(y^i)',\qquad i=1,\ldots,n
\label{eq-am-vector-variation}
$$

を得る。逆にこれらの式が成り立てば第一変分は零である。
時間を独立変数とする力学では、$y^i$ は系の一般化座標 $q^i$、$p^i$ に当たる引数は一般化速度 $\dot q^i$ である。
ここでの $p^i$ は傾きを指すため、一般化運動量とは区別する。

固定端を外すと何が変わるだろうか。例えば $x=b$ を固定したまま終点の**高さ**だけを自由にする。

:::proposition[自然境界条件]{#prop-am-natural-boundary}
$F$ は $C^2$ 級、$y\in C^2([a,b])$ とし、許容関数を $y(a)=A$ だけを課した $C^1$ 級関数全体とする。
$y$ が停留する（$\eta(a)=0$ を満たすすべての $\eta\in C^1$ で $\delta J[y;\eta]=0$）ための必要十分条件は、:eqref[eq-am-euler-lagrange] と

$$
F_p\big(b,y(b),y'(b)\big)=0
\label{eq-am-natural-boundary}
$$

が成り立つことである。
:::

:::proof
部分積分により

$$
\delta J[y;\eta]=F_p(b)\eta(b)+\int_a^b\left(F_y-\dv{}{x}F_p\right)\eta\,\dd x
$$

である。$\eta(a)=0$ のみが課されている。まず $\eta(b)=0$ も満たす $\eta$ に限ると、基本補題から :eqref[eq-am-euler-lagrange] が従う。
すると積分は常に零なので、$\eta(x)=(x-a)/(b-a)$ を選べば $F_p(b)=0$ を得る。逆は上の式から明らかである。
:::

:eqref[eq-am-natural-boundary] を :term[自然境界条件]{en="natural boundary condition" yomi="しぜんきょうかいじょうけん"} という。
許容関数に課さなかった端点条件を、停留性が自動的に要求するという意味である。
複数の関数の場合は、自由な成分ごとに $F_{p^i}(b)=0$ となる。終点が曲線上を動く場合の条件は演習問題で導く。

### 高階導関数を含む場合

弾性はりのたわみのように、被積分関数が二階導関数を含む問題もある。

:::proposition[Euler–Poisson 方程式]{#prop-am-euler-poisson}
$J[y]=\int_a^bF(x,y,y',y'')\,\dd x$ で、$F(x,y,p,r)$ は $C^3$ 級、$y$ は $C^4$ 級とし、両端で $y$ と $y'$ の値を固定する。$y$ が停留する（両端で $\eta=\eta'=0$ となるすべての $\eta\in C^2$ で第一変分が零）ための必要十分条件は

$$
F_y-\dv{}{x}F_p+\dvn{}{x}{2}F_r=0
\label{eq-am-euler-poisson}
$$

である。ここで $F_y,F_p,F_r$ は $(x,y,y',y'')$ で評価する。
:::

:::proof
:ref[prop-am-first-variation] と同じ理由で $\delta J=\int_a^b(F_y\eta+F_p\eta'+F_r\eta'')\,\dd x$ である。
仮定から $x\mapsto F_r$ は $C^2$ 級、$x\mapsto F_p$ は $C^1$ 級である。第二項を一回、第三項を二回部分積分すると

$$
\delta J=\Big[F_p\eta+F_r\eta'-\Big(\dv{}{x}F_r\Big)\eta\Big]_a^b+\int_a^b\Big(F_y-\dv{}{x}F_p+\dvn{}{x}{2}F_r\Big)\eta\,\dd x
$$

となり、境界項は $\eta=\eta'=0$ で消える。:ref[lem-am-fundamental-variation] を $k=2$ で使えば必要性が、上の式から十分性が従う。
:::

例えば曲げ剛性 $EI$ のはりに、単位長さ当たり $w(x)$ の荷重がかかるとき、エネルギーは $\int(\tfrac12EI(y'')^2-wy)\,\dd x$ であり、:eqref[eq-am-euler-poisson] は $EI\,y''''=w$ を与える。

### 独立変数または従属変数が陽に現れないとき

:::proposition[Beltrami の恒等式]{#prop-am-beltrami}
$F$ は $C^2$ 級とする。$C^2$ 級の停留曲線に沿って次式が成り立つ。

$$
\dv{}{x}(F-y'F_p)=F_x.
$$

特に $F$ が $x$ に陽に依存しなければ

$$
F-y'F_p=C \quad (C\text{ は定数})
\label{eq-am-beltrami}
$$

となる。この等式を **Beltrami の恒等式**と呼ぶ。
:::

:::proof
曲線に沿って微分すると

$$
\dv{}{x}(F-y'F_p)
=F_x+F_y y'+F_p y''-y''F_p-y'\dv{F_p}{x}
=F_x+y'\left(F_y-\dv{F_p}{x}\right).
$$

右辺の括弧は :eqref[eq-am-euler-lagrange] で零である。
:::

もう一つの簡約は、$F$ に従属変数 $y$ が陽に現れない場合である。
この場合は $F_y=0$ なので Euler–Lagrange 方程式から $F_p$ が一定になる。
どちらの量が一定かは、陽に現れない変数が $x$ か $y$ かで判定できる。
どちらの場合も、二階の方程式が一階の方程式に下がる。力学ではこれらがエネルギーと運動量の保存に当たる（[対称性と保存則](/physics/analytical-mechanics/symmetries/)）。

:::example[平面上の最短線]
$F=\sqrt{1+p^2}$ は $y$ を含まないので $F_p=p/\sqrt{1+p^2}$ が一定である。$p\mapsto p/\sqrt{1+p^2}$ は狭義単調増加なので $y'$ が一定であり、停留曲線は直線である。
それが実際に最短であることは、後の :ref[thm-am-convex-sufficient] で示す。
:::

:::example[Fermat の原理と屈折]{#exm-am-snell}
屈折率が高さだけに依存する媒質 $n(y)$ の中では、光速は $c/n(y)$ である。光は所要時間

$$
T[y]=\frac1c\int_a^bn(y)\sqrt{1+(y')^2}\,\dd x
$$

を停留させる経路を進む（:term[Fermat の原理]{en="Fermat's principle" yomi="ふぇるまーのげんり"}）。
被積分関数は $x$ を含まないので、Beltrami の恒等式から

$$
F-y'F_p=\frac{n(y)}{\sqrt{1+(y')^2}}=\text{定数}
$$

である。光線が層の法線（$y$ 軸）となす角を $\varphi$ とすると $\sin\varphi=1/\sqrt{1+(y')^2}$ なので、$n(y)\sin\varphi$ が経路に沿って一定となる。これは Snell の屈折の法則である。
:::

### 停留曲線は極小か

停留曲線が極小を与えるかを調べるには、一変数関数 $\eps\mapsto J[y+\eps\eta]$ の二階微分を見る。
$F$ が $C^2$ 級なら、:ref[prop-am-first-variation] と同じ議論で積分記号の下で二回微分でき、

$$
\delta^2J[y;\eta]=\left.\dvn{}{\eps}{2}J[y+\eps\eta]\right|_{\eps=0}
=\int_a^b\big(F_{yy}\eta^2+2F_{yp}\eta\eta'+F_{pp}(\eta')^2\big)\,\dd x
\label{eq-am-second-variation}
$$

となる。これを :term[第二変分]{en="second variation" yomi="だいにへんぶん"} という。
$y$ が弱い極小を与えるなら、$\eps\mapsto J[y+\eps\eta]$ は $\eps=0$ で極小をとるので、すべての許容方向で $\delta^2J[y;\eta]\ge0$ である。

:::proposition[Legendre の条件]{#prop-am-legendre}
$F$ は $C^2$ 級とする。$y\in C^1([a,b])$ が固定端の汎関数 :eqref[eq-am-functional] の弱い極小を与えるならば、すべての $x\in[a,b]$ で

$$
F_{pp}\big(x,y(x),y'(x)\big)\ge0
$$

である。
:::

:::proof
ある内点 $x_0$ で $F_{pp}<0$ と仮定する。連続性により、$\kappa>0$ と $\delta_0>0$ があって $[x_0-\delta_0,x_0+\delta_0]\subset(a,b)$ 上で $F_{pp}\le-\kappa$ となる。
曲線上での $\abs{F_{yy}},\abs{F_{yp}}$ の最大値を $C$ とする。
$\phi(s)=(1-s^2)^2$（$\abs s\le1$）、$\phi(s)=0$（$\abs s>1$）は $C^1$ 級である。$0<\delta\le\delta_0$ に対し $\eta_\delta(x)=\phi((x-x_0)/\delta)$ とおくと $\eta_\delta'(x)=\phi'((x-x_0)/\delta)/\delta$ であり、$s=(x-x_0)/\delta$ と置換して

$$
\delta^2J[y;\eta_\delta]\le C\delta\int_{-1}^1\phi^2\,\dd s+2C\int_{-1}^1\abs{\phi\phi'}\,\dd s-\frac{\kappa}{\delta}\int_{-1}^1(\phi')^2\,\dd s
$$

を得る。最後の積分は正なので、$\delta$ を十分小さくすると右辺は負になり、$\delta^2J\ge0$ に反する。
したがって内点で $F_{pp}\ge0$ であり、端点でも連続性から同じである。
:::

直観的には、振幅が小さくても急に振動する変形では $\eta'$ の項が支配的になり、その係数 $F_{pp}$ の符号が効く。
Legendre の条件は必要条件であって十分条件ではない（演習問題の $\int((y')^2-y^2)\,\dd x$ では $F_{pp}=2>0$ なのに極小にならない）。
十分条件として最も使いやすいのは凸性である。

:::theorem[凸性による十分条件]{#thm-am-convex-sufficient}
$F$ は $C^1$ 級で、各 $x\in[a,b]$ について $(y,p)\mapsto F(x,y,p)$ が $\R^2$ 上の凸関数であるとする。
固定端の停留曲線 $y_*\in C^1([a,b])$ は、同じ端点条件を満たすすべての許容関数 $y$ に対して $J[y]\ge J[y_*]$ を満たす。
さらに各 $x$ について $(y,p)\mapsto F$ が狭義凸ならば、最小を与える許容関数は $y_*$ だけである。
:::

:::proof
微分可能な凸関数 $\varphi$ は接平面より上にある：$\varphi(z+h)\ge\varphi(z)+\grad\varphi(z)\cdot h$。
実際、$g(s)=\varphi(z+sh)$ は凸なので差分商 $(g(s)-g(0))/s$ は $s>0$ について単調非減少であり、$g(1)-g(0)\ge\lim_{s\downarrow0}(g(s)-g(0))/s=g'(0)$ である。
$\eta=y-y_*$ とおき、各 $x$ でこの不等式を $z=(y_*,y_*')$、$h=(\eta,\eta')$ に適用して積分すると

$$
J[y]-J[y_*]\ge\int_a^b\big(F_y\eta+F_p\eta'\big)\,\dd x=\delta J[y_*;\eta]=0
$$

となる。狭義凸なら、$h\ne0$ である $x$ では不等式が狭義である。$\eta$ が恒等的に零でなければそのような $x$ があり、両辺の差は連続で非負なので、その近傍で正となって積分も正になる。
:::

例えば $F=\sqrt{1+p^2}$ は $F_{pp}=(1+p^2)^{-3/2}>0$ より $p$ について狭義凸で、$y$ には依存しないから $(y,p)$ について凸である。よって直線は端点を結ぶ曲線のうちで最短である。
等号が成り立つには各点で $\sqrt{1+(y_*'+\eta')^2}=\sqrt{1+(y_*')^2}+F_p\eta'$、すなわち $p$ についての狭義凸性から $\eta'=0$ が必要で、$\eta(a)=0$ より $\eta=0$ である。したがって最短線は直線に限る。

凸性が成り立たない場合の典型例として、調和振動子に対応する汎関数を調べる。そのための不等式を準備する。

:::lemma[Wirtinger の不等式]{#lem-am-wirtinger}
$\ell>0$ とし、$\eta\in C^1([0,\ell])$ は $\eta(0)=\eta(\ell)=0$ を満たすとする。このとき

$$
\int_0^\ell(\eta')^2\,\dd x\ge\left(\frac{\pi}{\ell}\right)^2\int_0^\ell\eta^2\,\dd x
\label{eq-am-wirtinger}
$$

であり、等号は $\eta(x)=c\sin(\pi x/\ell)$（$c$ は定数）のときに限る。
:::

:::proof
$k=\pi/\ell$、$w(x)=k\cot(kx)$（$0<x<\ell$）とおく。$w'=-k^2/\sin^2(kx)=-k^2-w^2$ を使うと、

$$
(\eta'-w\eta)^2+\dv{}{x}(w\eta^2)
=(\eta')^2-2w\eta\eta'+w^2\eta^2+(-k^2-w^2)\eta^2+2w\eta\eta'
=(\eta')^2-k^2\eta^2
$$

が成り立つ。$0<\delta<\ell/2$ として $[\delta,\ell-\delta]$ で積分すると

$$
\int_\delta^{\ell-\delta}\big((\eta')^2-k^2\eta^2\big)\,\dd x=\int_\delta^{\ell-\delta}(\eta'-w\eta)^2\,\dd x+\big[w\eta^2\big]_\delta^{\ell-\delta}.
$$

$\eta(0)=0$ と平均値の定理から $\abs{\eta(\delta)}\le\norm{\eta'}_0\delta$ であり、$\abs{w(\delta)}\le1/\delta$（$\tan u\ge u$ による）なので $\abs{w(\delta)\eta(\delta)^2}\le\norm{\eta'}_0^2\delta\to0$ である。$x=\ell-\delta$ でも同様である。
$\delta\to0$ とすれば $\int_0^\ell((\eta')^2-k^2\eta^2)\,\dd x=\int_0^\ell(\eta'-w\eta)^2\,\dd x\ge0$ を得る（右辺は非負の関数の広義積分である）。
等号が成り立てば $(0,\ell)$ で $\eta'=w\eta$ であり、$\big(\eta/\sin(kx)\big)'=(\eta'-w\eta)/\sin(kx)=0$ から $\eta=c\sin(kx)$ となる。逆にこの関数で等号が成り立つことは直接計算でわかる。
:::

:::example[区間の長さで変わる極小性]{#exm-am-jacobi-example}
$J[y]=\tfrac12\int_0^\ell((y')^2-y^2)\,\dd x$、$y(0)=y(\ell)=0$ とする。Euler–Lagrange 方程式は $y''+y=0$ で、$y=0$ は停留曲線である。
被積分関数は二次式なので、任意の許容関数 $\eta$ について展開に余りがなく $J[\eta]-J[0]=J[\eta]$ である。

1. $\ell<\pi$ なら、:ref[lem-am-wirtinger] により $J[\eta]\ge\tfrac12(\pi^2/\ell^2-1)\int_0^\ell\eta^2\,\dd x$ であり、$\eta\ne0$ なら正である。よって $y=0$ は唯一の最小解である。
2. $\ell=\pi$ なら $J\ge0$ だが、$\eta=c\sin x$ で $J=0$ となる。最小解は一意でない。
3. $\ell>\pi$ なら $\eta=\sin(\pi x/\ell)$ で $J[\eta]=\tfrac{\ell}{4}(\pi^2/\ell^2-1)<0$ であり、$y=0$ は極小でない。

$F_{pp}=1>0$ なので Legendre の条件はどの $\ell$ でも満たされる。極小性が失われる境目 $\ell=\pi$ は、$x=0$ で零になる方程式の解 $\sin x$ が再び零になる点である。
このような点を :term[共役点]{en="conjugate point" yomi="きょうやくてん"} と呼び、一般の汎関数でも共役点の有無が弱い極小性を左右する（Jacobi の条件）。一般論には立ち入らないが、この例は次章で調和振動子の作用が極小になる条件として再び現れる。
:::

## 最速降下線と懸垂線

### 重力の下で二点を結ぶ時間

鉛直面に水平右向きの $x$ 軸と鉛直**下向き**の $y$ 軸をとる。
質点を原点で静かに放し、一様な重力（加速度 $g>0$）の下、固定した滑らかな針金に沿って滑らせる。
摩擦、空気抵抗、質点の回転は無視する。垂直抗力は仕事をしないから、エネルギー保存より
$\tfrac12mv^2=mgy$、すなわち $v=\sqrt{2gy}$ である。
到達点を $(X,Y)$、$X,Y>0$ とし、途中では $y>0$ を保つ。所要時間が最小の曲線を :term[最速降下線]{en="brachistochrone" yomi="さいそくこうかせん"} という。
この問題は 1696 年に Johann Bernoulli が提出し、変分法誕生のきっかけとなった。

グラフ $y=y(x)$ の微小な長さは $\dd s=\sqrt{1+(y')^2}\,\dd x$ である。
したがって時間は

$$
T[y]=\int_0^X\frac{\sqrt{1+(y')^2}}{\sqrt{2gy}}\,\dd x.
\label{eq-am-descent-time}
$$

出発点では $v=0$ なので、積分は内部から端点に近づける広義積分とする。
被積分関数 $F=\sqrt{1+p^2}/\sqrt{2gy}$ は $y>0$ で滑らかであり、$F_{pp}=1/(\sqrt{2gy}(1+p^2)^{3/2})>0$ なので、$C^1$ 級の停留曲線は :ref[prop-am-regularity] により $C^2$ 級である。
端点を除いた区間で :eqref[eq-am-beltrami] を使うと、

$$
F-pF_p=\frac{1}{\sqrt{2gy}\sqrt{1+p^2}},
\qquad y(1+p^2)=2a
$$

となる。$a>0$ は長さの次元をもつ積分定数である。

### サイクロイドと端点条件

$y=a(1-\cos\theta)=2a\sin^2(\theta/2)$ と置き、$0<\theta<2\pi$ とする。
このとき $1+p^2=2a/y=1/\sin^2(\theta/2)$、すなわち $p^2=\cot^2(\theta/2)$ である。
初めは降下し、最深点を越えれば上昇する枝を一続きに記述するには $p=\cot(\theta/2)$ と選べばよい。
$\dd y/\dd\theta=a\sin\theta$ と $p=\dd y/\dd x$ を組み合わせると

$$
\dv{x}{\theta}=\frac{a\sin\theta}{\cot(\theta/2)}=2a\sin^2(\theta/2)=a(1-\cos\theta)
$$

を得る。最深点 $\theta=\pi$ での値も連続に定まる。
原点を通る条件を入れて積分すれば

$$
x=a(\theta-\sin\theta),\qquad y=a(1-\cos\theta).
\label{eq-am-cycloid}
$$

これは半径 $a$ の円を直線上で転がすとき、円周上の一点が描く :term[サイクロイド]{en="cycloid" yomi="さいくろいど"} である。
始点での接線は鉛直である。$\theta=\pi$ で接線は水平となり、$y=2a$ に達する。

終点に対応する $\theta_1$ と尺度 $a$ は

$$
\frac{X}{Y}=\frac{\theta_1-\sin\theta_1}{1-\cos\theta_1},\qquad
a=\frac{Y}{1-\cos\theta_1}
\label{eq-am-cycloid-endpoint}
$$

によって一意に決まる。実際、比を $h(\theta)$ と書くと、
$h'(\theta)$ の分母は $(1-\cos\theta)^2$、分子は $2(1-\cos\theta)-\theta\sin\theta$ である。
$0<\theta<\pi$ では $u=\theta/2$ として $\sin u-u\cos u>0$ である。これは左辺が $u=0$ で零、その導関数が $u\sin u>0$ であることに従う。
よって分子 $4\sin u(\sin u-u\cos u)$ は正である。$\pi\le\theta<2\pi$ でも $\sin\theta\le0$ から正である。
さらに $\theta\downarrow0$ で $h(\theta)\sim\theta/3$、$\theta\uparrow2\pi$ で $h(\theta)\to\infty$ なので、$h$ は $(0,2\pi)$ から $(0,\infty)$ への狭義単調増加な全単射である。

所要時間もこのパラメータで計算できる。$0<\theta<2\pi$ では

$$
\dv{s}{\theta}=a\sqrt{(1-\cos\theta)^2+\sin^2\theta}=2a\sin(\theta/2),\qquad
v=2\sqrt{ga}\sin(\theta/2),\qquad
T=\sqrt{\frac ag}\,\theta_1.
$$

出発直後は $v$ が小さいが、$\dd s/\dd\theta$ も同じ割合で小さく、時間は有限である。
$X/Y>\pi/2$ なら $\theta_1>\pi$ なので、最速の経路は終点の高さより一度深く降りる。

```tikz caption="落差に比べて水平距離が大きいと、最速降下線は終点より深く降りる。破線は同じ端点を結ぶ直線である。"
\begin{tikzpicture}[>=Stealth]
\draw[->] (0,0) -- (6.2,0) node[right] {$x$};
\draw[->] (0,0) -- (0,-2.5) node[below] {$y$};
\coordinate (B) at ({1.5*pi+1},-1);
\draw[dashed] (0,0) -- (B);
\draw[thick,domain=0:4.712389,samples=70,variable=\t] plot ({\t-sin(\t r)},{cos(\t r)-1});
\fill (0,0) circle (1.5pt) node[above left] {$O$};
\fill (3.14159,-2) circle (1.5pt) node[below] {$\theta=\pi$};
\fill (B) circle (1.5pt) node[right] {$(X,Y)$};
\end{tikzpicture}
```

### 他の経路との比較

停留条件だけから最速とはいえない。しかも降下時間の被積分関数は $(y,p)$ について凸でないので、:ref[thm-am-convex-sufficient] も使えない。
そこで、原点から各点までサイクロイドをたどる時間を関数として用いる。

:::theorem[サイクロイドの最速性]{#thm-am-brachistochrone}
$X,Y>0$ とし、$\theta_1,a$ を :eqref[eq-am-cycloid-endpoint] で定める。原点と $(X,Y)$ を結び、原点以外では $y>0$ にある区分的に $C^1$ 級の任意の曲線 $\gamma$ について

$$
\int_\gamma\frac{\dd s}{\sqrt{2gy}}\ge\sqrt{\frac ag}\,\theta_1
$$

であり、等号はサイクロイド :eqref[eq-am-cycloid] で成り立つ。曲線はグラフでなくてもよく、途中で水平方向に逆行してもよい。
:::

:::proof
写像 $\Psi:(a,\theta)\mapsto(x,y)=(a(\theta-\sin\theta),a(1-\cos\theta))$ を $(0,\infty)\times(0,2\pi)$ 上で考える。
上で示したことから、第一象限の各点 $(x,y)$ に対し $h(\theta)=x/y$ となる $\theta$ がただ一つあり、$a=y/(1-\cos\theta)$ と合わせて $\Psi$ は第一象限への全単射である。
その Jacobi 行列式は

$$
\det\begin{pmatrix}\theta-\sin\theta&a(1-\cos\theta)\\1-\cos\theta&a\sin\theta\end{pmatrix}
=-a\big(2(1-\cos\theta)-\theta\sin\theta\big)
$$

で、括弧内は $h'$ の分子なので零でない。逆関数定理により逆写像は滑らかである。
したがって $U(x,y)=\sqrt{a/g}\,\theta$ は第一象限で滑らかな関数である。

$U$ の微分を求める。全微分 $\dd x=(\theta-\sin\theta)\dd a+a(1-\cos\theta)\dd\theta$、$\dd y=(1-\cos\theta)\dd a+a\sin\theta\,\dd\theta$ と
$\cot(\theta/2)(1-\cos\theta)=\sin\theta$、$\cot(\theta/2)\sin\theta=1+\cos\theta$ から、

$$
\dd x+\cot(\theta/2)\,\dd y=\theta\,\dd a+2a\,\dd\theta
$$

である。これを $2\sqrt{ag}$ で割ったものは $\dd U=\frac{\theta}{2\sqrt{ag}}\dd a+\sqrt{\frac ag}\dd\theta$ に一致する。よって

$$
\grad U=\frac{1}{2\sqrt{ag}}\big(1,\cot(\theta/2)\big),\qquad
\norm{\grad U}^2=\frac{1+\cot^2(\theta/2)}{4ag}=\frac{1}{4ag\sin^2(\theta/2)}=\frac{1}{2gy}
$$

を得る。

左半平面では $U(x,y)=U(-x,y)$ と定める。$x\downarrow0$、$y>0$ のときは $h(\theta)\to0$ より $\theta\to0$、$a=y/(1-\cos\theta)\sim2y/\theta^2$ となり、
$U\to\sqrt{2y/g}$、$U_x=1/(2\sqrt{ag})\to0$、$U_y\to1/\sqrt{2gy}$ である。左側からの極限も同じなので、$U$ は上半平面 $y>0$ 全体で $C^1$ 級であり、そこで $\norm{\grad U}^2=1/(2gy)$ が成り立つ。
また原点では $U\to0$ である。実際、$U^2/(\abs{x}+y)$ は
$\theta^2/[g(\theta-\sin\theta+1-\cos\theta)]$ に等しく、$0<\theta<2\pi$ で有界（$\theta\to0$ の極限は $2/g$）である。

$\gamma$ を $[0,1]$ 上の区分的に $C^1$ 級の写像とし、$\gamma(0)$ が原点とする。$0<t_0<1$ について、Cauchy–Schwarz の不等式より

$$
U(X,Y)-U(\gamma(t_0))=\int_{t_0}^1\grad U\cdot\gamma'\,\dd t
\le\int_{t_0}^1\norm{\grad U}\,\abs{\gamma'}\,\dd t
=\int_{t_0}^1\frac{\abs{\gamma'}\,\dd t}{\sqrt{2gy}}
$$

である。$t_0\downarrow0$ とすると左辺は $U(X,Y)=\sqrt{a/g}\,\theta_1$ に、右辺は降下時間に収束する。
サイクロイドでは接ベクトル $(\dd x/\dd\theta,\dd y/\dd\theta)$ が $(1,\cot(\theta/2))$ に正の数を掛けたものなので Cauchy–Schwarz の等号が成り立ち、降下時間は $\sqrt{a/g}\,\theta_1$ に等しい。
:::

この証明で使った $U$ は「各点に最速で到達する時間」であり、$\norm{\grad U}^2=1/(2gy)$ は光学の波面に対するアイコナール方程式と同じ形をしている。
この見方は [Hamilton–Jacobi 理論](/physics/analytical-mechanics/hamilton-jacobi/) で一般化される。

### 回転面の面積と懸垂線

ここから $y$ は回転軸からの距離を表す。平行な二つの円環を共通の中心軸に垂直に置く。
軸方向の座標を $x$ とし、正の半径 $y(x)$ のグラフを $x$ 軸のまわりに回して膜を作る。
圧力差はなく、表面張力は一定とする。このとき膜のエネルギーは面積に比例するので、平衡形は面積の停留曲面である。
幅 $\dd s$ の帯の面積は円周 $2\pi y$ と $\dd s$ の積なので、面積は

$$
A[y]=2\pi\int_{x_1}^{x_2} y\sqrt{1+(y')^2}\,\dd x .
$$

$y>0$ では $F_{pp}=y(1+p^2)^{-3/2}>0$ なので、停留曲線は $C^2$ 級である。
固定した円環では両端の半径が指定される。定数 $2\pi$ を除き :eqref[eq-am-beltrami] を用いると

$$
\frac{y}{\sqrt{1+(y')^2}}=a,\qquad a>0.
$$

この式を $y'$ について解くと符号の選び方が残り、水平接線の点で扱いにくい。
そこで元の Euler–Lagrange 方程式に戻り、
$F_y=\sqrt{1+p^2}$、$F_p=yp/\sqrt{1+p^2}$ を代入して分母を払うと

$$
yy''=1+(y')^2=\frac{y^2}{a^2}.
$$

$y>0$ なので $y''=y/a^2$ である。一般解 $y=C e^{x/a}+D e^{-x/a}$ を第一積分 $y^2=a^2(1+(y')^2)$ に代入すると $4CD=a^2$ を得る。
正の半径には $C,D>0$ が必要なので、水平移動を表す定数 $b$ で書き直せば

$$
y(x)=a\cosh\frac{x-b}{a}.
\label{eq-am-catenary}
$$

この曲線を :term[懸垂線]{en="catenary" yomi="けんすいせん"}、これを回転してできる面を :term[懸垂面]{en="catenoid" yomi="けんすいめん"} という。
例えば $x=\pm d$ に同半径 $R$ の円環があれば、$\cosh((d-b)/a)=\cosh((-d-b)/a)$ から $b=0$、$R=a\cosh(d/a)$ となる。
解の個数を調べるには $u=d/a>0$ とおくとよい。端点方程式は $R/d=\cosh u/u$ となり、右辺は $u\downarrow0$ と $u\to\infty$ の両方で発散する。
その導関数 $(u\sinh u-\cosh u)/u^2$ の符号は $u\tanh u-1$ の符号で決まる。$u\tanh u$ は零から狭義単調に増加して無限大に向かうので、右辺にはただ一つの極小値がある。
数値的には極小点は $u_*\simeq1.1997$、極小値は $\cosh u_*/u_*\simeq1.5089$ である。
したがって $R/d>1.5089$ なら、首の太さが異なる二つの懸垂面が同じ円環を結ぶ。等しいときは一つ、下回るときは滑らかな正の半径をもつ解はない。
後者では、二つの円環を離していくと石鹸膜が首でくびれて切れ、各円環に平らな円板の膜が残る。
二つの停留解が得られる場合、面積が最小かどうかを調べるにはさらに第二変分や他の曲面との比較が必要になる。

ここまでは回転対称な変形を用いた。対称性を崩す変形についても第一変分が消えることを確かめよう。
方位角を $\varphi$ とし、半径を $r(x,\varphi)$ と書けば、曲面の位置は $(x,r\cos\varphi,r\sin\varphi)$ である。
二つの接ベクトル $(1,r_x\cos\varphi,r_x\sin\varphi)$ と $(0,r_\varphi\cos\varphi-r\sin\varphi,r_\varphi\sin\varphi+r\cos\varphi)$ の外積の大きさから、面積要素は $\sqrt{r^2(1+r_x^2)+r_\varphi^2}\,\dd x\,\dd\varphi$ となる。
回転対称な解 $r=y(x)$ では $r_\varphi=0$ なので、$r=y+\eps\eta(x,\varphi)$ の一次まででは $r_\varphi^2=\eps^2\eta_\varphi^2$ の寄与が消える。
残る第一変分は、各 $\varphi$ で上と同じ一次元の第一変分をとり、それを角度で積分したものである。両端の円環を固定すると各角度で零であり、懸垂面はこれらの変形に対しても停留する。

## 拘束条件付き変分問題

### 積分を一つ指定する場合

汎関数 $J[y]$ を、端点条件に加えて

$$
K[y]=\int_a^b G(x,y,y')\,\dd x=K_0
$$

を保ったまま停留させたい。曲線の長さを指定する場合を典型として、この種の問題を :term[等周問題]{en="isoperimetric problem" yomi="とうしゅうもんだい"} という。
ここで $y$ が**拘束の下で停留する**とは、$y$ を通り拘束を満たす許容関数の $C^1$ 級の一変数族 $y_s$（$y_0=y$、$(s,x)\mapsto y_s(x)$ とその $x$ 微分が $s$ について $C^1$ 級）すべてについて、$\dd J[y_s]/\dd s|_{s=0}=0$ となることとする。

:::proposition[等周問題の乗数]{#prop-am-isoperimetric}
$F,G$ は $C^2$ 級とする。許容関数 $y$ が拘束の下で停留し、端点で零となるある方向 $\zeta$ で $\delta K[y;\zeta]\ne0$ とする。
このとき一つの実定数 $\lambda$ が存在し、$y$ は $J+\lambda(K-K_0)$ を、固定端のすべての方向について停留させる。
この $\lambda$ を **Lagrange 乗数** という（:ref[prop-am-ideal-multiplier] と同じ用語である）。
:::

:::proof
任意の固定端の方向 $\eta$ を選び、二変数の曲線族 $y+s\eta+u\zeta$ を作る。
:ref[prop-am-first-variation] の議論から、$k(s,u)=K[y+s\eta+u\zeta]$ と $j(s,u)=J[y+s\eta+u\zeta]$ は原点の近くで $C^1$ 級であり、偏導関数は第一変分で与えられる。
$(s,u)=(0,0)$ で $k=K_0$、しかも $u$ 偏微分は $\delta K[y;\zeta]\ne0$ である。
陰関数定理により $u(0)=0$、$k(s,u(s))=K_0$ となる $C^1$ 級関数 $u(s)$ が局所的に存在する。
族 $y+s\eta+u(s)\zeta$ は拘束を満たすので、停留の仮定から $\dd j(s,u(s))/\dd s=0$（$s=0$）である。拘束と合わせて $s=0$ で微分すると

$$
\delta K[\eta]+u'(0)\delta K[\zeta]=0,\qquad
\delta J[\eta]+u'(0)\delta J[\zeta]=0
$$

となる。ここでは基点 $y$ を省略した。$u'(0)$ を消去し、$\lambda=-\delta J[\zeta]/\delta K[\zeta]$ とおけば
$\delta J[\eta]+\lambda\delta K[\eta]=0$ を得る。$\lambda$ は $\eta$ の選び方に依存しない。
:::

したがって被積分関数を $F+\lambda G$ に置き換えた Euler–Lagrange 方程式を解き、元の端点条件と $K[y]=K_0$ を連立すればよい。未知数は方程式の二つの積分定数と $\lambda$ の三つで、条件も三つである。
もし $\delta K$ 自体がすべての方向で零なら、上の割り算はできない。そのような候補（$K$ の停留曲線）は別に調べる必要がある。

拘束が $m$ 個の積分 $K_\alpha[y]=K_{\alpha,0}$ の場合も同様である。方向 $\zeta_1,\ldots,\zeta_m$ で行列 $A_{\alpha\beta}=\delta K_\alpha[\zeta_\beta]$ が正則になるとする。
族 $y+s\eta+\sum_\beta u_\beta\zeta_\beta$ に陰関数定理を適用して $u(s)\in\R^m$ を定めると、$\delta K_\alpha[\eta]+\sum_\beta A_{\alpha\beta}u'_\beta(0)=0$ と $\delta J[\eta]+\sum_\beta\delta J[\zeta_\beta]u'_\beta(0)=0$ を得る。
前者から $u'(0)=-A^{-1}\delta K[\eta]$ を後者に代入し、$\lambda_\alpha=-\sum_\beta\delta J[\zeta_\beta](A^{-1})_{\beta\alpha}$ とおけば、$\delta J[\eta]+\sum_\alpha\lambda_\alpha\delta K_\alpha[\eta]=0$ となる。

:::example[Dido の問題]{#exm-am-dido}
長さ $\ell$ の紐の両端を直線（海岸）上の二点 $(\pm d,0)$ に固定し、紐と直線の間の面積を最大にする。$2d<\ell<\pi d$ とし、紐をグラフ $y(x)\ge0$ で表す。
面積 $J=\int_{-d}^dy\,\dd x$ を、長さ $K=\int_{-d}^d\sqrt{1+(y')^2}\,\dd x=\ell$ の下で停留させる。
$K$ の停留曲線は直線 $y=0$ だけで、その長さは $2d\ne\ell$ なので、乗数を使える。$F+\lambda G=y+\lambda\sqrt{1+p^2}$ は $x$ を含まないから、Beltrami の恒等式より

$$
y+\frac{\lambda}{\sqrt{1+(y')^2}}=C .
$$

$\lambda=0$ なら $y$ は定数となり端点条件から $y=0$ で長さが合わない。よって $\lambda\ne0$ で、$y-C=-\lambda\cos\psi$、$y'=\tan\psi$（$\abs\psi<\pi/2$）と置ける。すると $\dd x=\dd y/\tan\psi=\lambda\cos\psi\,\dd\psi$ から $x-x_0=\lambda\sin\psi$ であり、

$$
(x-x_0)^2+(y-C)^2=\lambda^2 .
$$

停留曲線は円弧である。対称性から $x_0=0$ であり、半径 $\rho=\abs\lambda$、弦の半分 $d$ の円弧の中心角の半分を $\alpha$ とすると $d=\rho\sin\alpha$、$\ell=2\rho\alpha$ である。
$\ell/(2d)=\alpha/\sin\alpha$ は $0<\alpha<\pi/2$ で $1$ から $\pi/2$ まで狭義単調に増えるので、条件 $2d<\ell<\pi d$ の下で $\alpha$ はただ一つに決まる。
この円弧が実際に面積を最大にすることは、閉曲線についての等周不等式（[平面曲線](/math/curves-surfaces/plane-curves/)）から従う。
:::

### 長さを固定して鎖をつるす

単位長さ当たりの質量 $\rho$ が一定で、曲げに抵抗しない伸びない鎖を考える。
両端を $(-d,0),(d,0)$ に止め、鉛直**上向き**を $y$ とする。長さ $\ell>2d$ として、滑らかなグラフ状のつり合いを求める。
重力による位置エネルギーと長さは

$$
V[y]=\rho g\int_{-d}^d y\sqrt{1+(y')^2}\,\dd x,\qquad
K[y]=\int_{-d}^d\sqrt{1+(y')^2}\,\dd x .
$$

長さそのものが停留するグラフでは $y'/\sqrt{1+(y')^2}$ が一定、したがって直線になる。
端点を結ぶ直線の長さは $2d<\ell$ なので、解の候補はこの直線ではなく、:ref[prop-am-isoperimetric] の仮定を満たす。
$V+\lambda(K-\ell)$ の被積分関数 $(\rho gy+\lambda)\sqrt{1+p^2}$ に Beltrami の恒等式を適用すると

$$
\frac{\rho g y+\lambda}{\sqrt{1+(y')^2}}=C
$$

となる。$C=0$ なら $y$ は定数 $-\lambda/(\rho g)$ であり、端点条件から $y=0$ となって長さの条件に反する。
したがって $C\ne0$ である。$z=y+\lambda/(\rho g)$、$a=C/(\rho g)\ne0$ とおくと $z/\sqrt{1+(z')^2}=a$ であり、$z$ は $a$ と同じ符号をもち零にならない。
元の Euler–Lagrange 方程式から回転面と同じ計算で $zz''=1+(z')^2=z^2/a^2$ を得るので、$z=Ce^{x/a}+De^{-x/a}$、$4CD=a^2$ である。$C,D$ は $z$ と同符号なので $z=a\cosh((x-b)/a)$ と書ける。
同じ高さの端点から $b=0$ となり、

$$
y(x)=a\big(\cosh(x/a)-\cosh(d/a)\big),\qquad
\ell=2\abs{a}\sinh(d/\abs{a}).
\label{eq-am-hanging-chain}
$$

長さの式では $\sqrt{1+\sinh^2(x/a)}=\cosh(x/a)$ を積分した。
$u=d/\abs a>0$ とすると $\ell/(2d)=\sinh u/u$ である。この比は $u\downarrow0$ で $1$ に近づき、無限大まで狭義単調に増える。
実際、導関数の分子 $u\cosh u-\sinh u$ は $u=0$ で零であり、その導関数は $u\sinh u>0$ である。
よって $\abs a$ はただ一つに決まり、停留曲線は $a>0$ と $a<0$ の二つである。

$a>0$ は下にたわむ懸垂線、$a<0$ はそれを上下に反転したアーチである。
両者の位置エネルギーは符号だけが異なり、$a>0$ では $y\le0$ なので $V<0$、アーチでは $V>0$ である。したがってアーチは最小を与えない。
物理的には、アーチの形を保つには鎖が圧縮力を支えなければならない。次の注意のように、$a$ は鎖の張力の水平成分を $\rho g$ で割ったものであり、$a<0$ は圧縮を意味する。
糸や鎖は圧縮に耐えないので、実現するのは $a>0$ の懸垂線である。

:::remark[力のつり合いによる確認]
最下点 $x=0$ から $x$ までの鎖の部分に働く力を考える。張力は鎖の接線方向に働き、水平方向には外力がないので、張力の水平成分 $H$ は鎖に沿って一定である。
鉛直方向のつり合いから、点 $x$ での張力の鉛直成分はこの部分の重さ $\rho g s(x)$ に等しい。ここで $s(x)$ は最下点からの弧長である。
よって $y'=\rho gs/H$ であり、微分して $y''=(\rho g/H)\sqrt{1+(y')^2}$ を得る。$w=y'$ についてこれを解けば $w=\sinh(\rho gx/H)$、したがって $y=(H/\rho g)\cosh(\rho gx/H)+$ 定数となり、$a=H/(\rho g)$ である。
:::

### 各点で等式を課す場合

積分一つを指定する場合と異なり、曲線の全点で $H(x,y^1,\ldots,y^n)=0$ を課す場合がある。
この場合は各点で反力に相当する量が変わるため、乗数も一般に関数 $\lambda(x)$ となる。
以下、$E_i=\dd F_{p^i}/\dd x-F_{y^i}$ と書く。これは :eqref[eq-am-vector-variation] の左辺である。

:::theorem[各点の拘束と乗数関数]{#thm-am-pointwise-constraint}
$F$ と $H$ は $C^2$ 級とする。$C^2$ 級の曲線 $y=(y^1,\ldots,y^n)$ が拘束 $H(x,y(x))=0$ を満たし、曲線上で $\nabla_yH=(H_{y^1},\ldots,H_{y^n})\ne0$ とする。
$y$ が、端点が同じで拘束を満たす $C^2$ 級の曲線の中で停留する（拘束を満たす族について第一変分が零）ならば、連続関数 $\lambda(x)$ があって

$$
E_i=\lambda(x)H_{y^i},\qquad i=1,\ldots,n
\label{eq-am-pointwise-multiplier}
$$

が $(a,b)$ の各点で成り立つ。$H$ は傾きを含まないので、これは $F+\lambda(x)H$ に対する Euler–Lagrange 方程式にほかならない。
:::

:::proof
内点 $x_0$ を一つとる。番号を付け替えて $H_{y^n}(x_0,y(x_0))\ne0$ としてよい。
陰関数定理により、点 $(x_0,y^1(x_0),\ldots,y^{n-1}(x_0))$ の近傍 $W$ と $C^2$ 級関数 $g$ があって、$(x_0,y(x_0))$ の近くでは $H(x,y)=0$ と $y^n=g(x,y^1,\ldots,y^{n-1})$ が同値になる。
$x_0$ を内部に含む閉区間 $I\subset(a,b)$ を十分短くとり、$I$ 上で $(x,y^1(x),\ldots,y^{n-1}(x))\in W$ かつ $H_{y^n}\ne0$ とする。このとき $I$ 上で $y^n=g(x,y^1,\ldots,y^{n-1})$ である。

$I$ の内部の閉区間の外で零になる $C^2$ 級関数 $\eta^1,\ldots,\eta^{n-1}$ を任意に選ぶ。$I$ 上で

$$
\tilde y^i=y^i+\eps\eta^i\ (i<n),\qquad \tilde y^n=g(x,\tilde y^1,\ldots,\tilde y^{n-1})
$$

とおき、$I$ の外では $\tilde y=y$ とする。$\abs\eps$ が小さければ $\tilde y$ は $W$ に収まり、拘束を満たし、$\eta^i$ が零になる $I$ の端の近くで $y$ と一致するので $C^2$ 級につながる。端点も同じである。
$\eps$ で微分すると、$\eta^n=\partial\tilde y^n/\partial\eps|_{\eps=0}=\sum_{i<n}g_{y^i}\eta^i=-\sum_{i<n}(H_{y^i}/H_{y^n})\eta^i$ である。最後の等号は $H(x,y^1,\ldots,y^{n-1},g)=0$ を微分して得られる。
$\eta^n$ も $I$ の内部の閉区間の外で零なので、部分積分の境界項は消え、停留条件は

$$
0=\int_I\sum_{i=1}^nE_i\eta^i\,\dd x=\int_I\sum_{i=1}^{n-1}\left(E_i-E_n\frac{H_{y^i}}{H_{y^n}}\right)\eta^i\,\dd x
$$

となる。各 $\eta^i$ を独立に選び :ref[lem-am-fundamental-variation]（$k=2$）を使えば、$I$ の内部で $E_i-E_nH_{y^i}/H_{y^n}=0$ $(i<n)$ である。
$\lambda=E_n/H_{y^n}$ とおけば、$i<n$ についても $i=n$ についても :eqref[eq-am-pointwise-multiplier] が成り立つ。

最後に、局所的に得た $\lambda$ をつなぐ。:eqref[eq-am-pointwise-multiplier] が成り立てば、$\nabla_yH\ne0$ なので $\lambda=\sum_iE_iH_{y^i}/\abs{\nabla_yH}^2$ であり、$\lambda$ はどの成分で解いたかによらずこの式で決まる。
よって各点の近くで得た $\lambda$ は重なりで一致し、$(a,b)$ 上の連続関数を定める。
:::

$\lambda$ の変分をとれば元の拘束 $H=0$ が戻る。したがって $F+\lambda H$ を $(y,\lambda)$ の汎関数とみなして変分すれば、方程式と拘束を一度に得る。
:eqref[eq-am-pointwise-multiplier] の $\lambda H_{y^i}$ は、次章で拘束力の一般化力として現れる。

:::example[球面上の測地線]{#exm-am-sphere-geodesic}
単位球面上の曲線 $\vb{r}(t)$ に対し、$J=\int_{t_0}^{t_1}\tfrac12\abs{\dot{\vb{r}}}^2\,\dd t$ を拘束 $H=\tfrac12(\abs{\vb{r}}^2-1)=0$ の下で停留させる。
$\nabla H=\vb{r}\ne0$ なので :ref[thm-am-pointwise-constraint] を使えて、$E_i=\ddot r_i$ から

$$
\ddot{\vb{r}}=\lambda(t)\vb{r}
$$

を得る。$\vb{r}\times\ddot{\vb{r}}=0$ より $\dv{}{t}(\vb{r}\times\dot{\vb{r}})=0$ であり、$\vb{h}=\vb{r}\times\dot{\vb{r}}$ は一定である。
$\vb{h}\ne0$ なら $\vb{r}\cdot\vb{h}=0$ が常に成り立つので、曲線は原点を通り $\vb{h}$ に垂直な平面と球面の交わり、すなわち大円上にある。
さらに $\abs{\vb{r}}^2=1$ を二回微分すると $\vb{r}\cdot\ddot{\vb{r}}+\abs{\dot{\vb{r}}}^2=0$、したがって $\lambda=-\abs{\dot{\vb{r}}}^2$ である。
$\dv{}{t}\abs{\dot{\vb{r}}}^2=2\dot{\vb{r}}\cdot\lambda\vb{r}=0$ なので速さは一定である。
力学的には、球面に拘束された自由粒子が一定の速さで大円上を動き、$\lambda\vb{r}$ がその向心力を与えている。
:::

## 演習問題

::::exercise{level=basic}
$J[y]=\int_0^1((y')^2/2+2y)\,\dd x$、$y(0)=y(1)=0$ の停留曲線を求め、それが唯一の最小解であることを示せ。
:::solution
$F_p=y'$、$F_y=2$ なので $y''=2$ である。二度積分して端点条件を入れると $y_*(x)=x^2-x$ を得る。
任意の許容曲線を $y_*+\eta$ と書き、差を展開すると

$$
J[y_*+\eta]-J[y_*]=\int_0^1(y_*'\eta'+2\eta)\,\dd x+
\tfrac12\int_0^1(\eta')^2\,\dd x.
$$

初めの積分は部分積分と $y_*''=2$ により零である。残る項は非負で、零なら $\eta'=0$、固定端より $\eta=0$ である。
（$F$ は $(y,p)$ について凸で $p$ について狭義凸なので、:ref[thm-am-convex-sufficient] からも同じ結論が得られる。）
:::
::::

::::exercise{level=basic}
長さ $J[y]=\int_0^b\sqrt{1+(y')^2}\,\dd x$ に対し、始点を $(0,A)$、終点を直線 $x=b$ 上で自由とする。停留曲線と最小の長さを求めよ。
:::solution
内点では $F_p=y'/\sqrt{1+(y')^2}$ が一定になる。終端での自然境界条件 :eqref[eq-am-natural-boundary] は $F_p(b)=0$ だから、この定数は零である。
よって $y'=0$、$y=A$ となる。その長さは $b$ であり、任意の競合曲線でも $\sqrt{1+(y')^2}\ge1$ から長さは $b$ 以上である。
:::
::::

::::exercise{level=basic}
最速降下線の終点を $(\pi R,2R)$（$R>0$）とする。曲線の尺度、到達時間、同じ端点を結ぶ直線を滑った場合の所要時間を求めよ。
:::solution
:eqref[eq-am-cycloid-endpoint] の比は $\pi R/(2R)=\pi/2=h(\pi)$ であり、$\theta_1=\pi$、$a=2R/(1-\cos\pi)=R$ となる。
従って最速時間は $\pi\sqrt{R/g}$ である。直線の長さは $s=R\sqrt{\pi^2+4}$、斜面方向の加速度は $g(2R/s)$ と一定なので、

$$
s=\tfrac12 g(2R/s)t^2,\qquad t=\sqrt{\pi^2+4}\sqrt{R/g}\simeq3.72\sqrt{R/g}.
$$

直線のほうが $\pi\simeq3.14$ より大きい係数をもち、時間を要する。:ref[thm-am-brachistochrone] と整合する。
:::
::::

::::exercise{level=basic}
一様な荷重を受け、両端を固定したはりのたわみ $y(x)$ が $J[y]=\int_0^1(\tfrac12(y'')^2-y)\,\dd x$ を停留させるとする。端点条件 $y(0)=y'(0)=y(1)=y'(1)=0$ の下で $y$ を求め、最大たわみを計算せよ。
:::solution
$F=\tfrac12r^2-y$ なので $F_y=-1$、$F_p=0$、$F_r=y''$ である。:eqref[eq-am-euler-poisson] は $-1+y''''=0$、すなわち $y''''=1$ となる。
一般解は $y=x^4/24+c_3x^3+c_2x^2+c_1x+c_0$ である。$y(0)=y'(0)=0$ から $c_0=c_1=0$、$y(1)=y'(1)=0$ から

$$
\tfrac1{24}+c_3+c_2=0,\qquad \tfrac16+3c_3+2c_2=0
$$

で、$c_3=-\tfrac1{12}$、$c_2=\tfrac1{24}$ を得る。よって $y=x^2(1-x)^2/24$ である。
対称性から最大は $x=1/2$ で、その値は $(1/4)(1/4)/24=1/384$ である。
:::
::::

::::exercise{level=basic}
屈折率が $n(y)^2=n_0^2(1+\kappa y)$（$n_0,\kappa>0$）で与えられる層状媒質の中の光線は放物線であることを示せ。
:::solution
$F=n(y)\sqrt{1+p^2}$ では $F_y=n'(y)\sqrt{1+p^2}$、$F_p=n(y)p/\sqrt{1+p^2}$ である。$\dd F_p/\dd x=n'(y')^2/\sqrt{1+(y')^2}+ny''/(1+(y')^2)^{3/2}$ なので、Euler–Lagrange 方程式は

$$
\frac{n\,y''}{(1+(y')^2)^{3/2}}=\frac{n'}{\sqrt{1+(y')^2}},\qquad\text{すなわち}\qquad y''=\frac{n'}{n}\big(1+(y')^2\big)
$$

となる。$n^2=n_0^2(1+\kappa y)$ から $n'/n=(n^2)'/(2n^2)=\kappa/(2(1+\kappa y))$ である。
一方 :ref[exm-am-snell] の第一積分 $n/\sqrt{1+(y')^2}=C$（$C>0$）から $1+(y')^2=n_0^2(1+\kappa y)/C^2$ である。代入すると

$$
y''=\frac{\kappa n_0^2}{2C^2}
$$

は定数であり、光線は放物線である。
これは地表付近で屈折率が高さとともに変わる大気中の光線の曲がり（蜃気楼）の簡単なモデルである。
:::
::::

::::exercise{level=advanced}
$J[y]=\tfrac12\int_0^{2\pi}((y')^2-y^2)\,\dd x$ で端点をともに零とする。$y=0$ が停留するが、局所的な最大にも最小にもならないことを示せ。
:::solution
方程式は $y''+y=0$ なので $y=0$ は停留する。$\eta_n(x)=\sin(nx/2)$（正整数 $n$）は固定端の条件を満たす。二乗した正弦と余弦の $[0,2\pi]$ 上の積分はいずれも $\pi$ であり、

$$
J[\eps\eta_n]=\tfrac{\pi\eps^2}{2}(n^2/4-1).
$$

$n=1$ では負、$n=3$ では正となる。$\norm{\eps\eta_n}_1=\abs\eps(1+n/2)$ はいくらでも小さくできるので、$C^1$ ノルムで測った原点の任意の近傍に両方の符号の値がある。
:ref[exm-am-jacobi-example] でいえば、区間の長さ $2\pi$ が共役点 $\pi$ を越えている場合である。
:::
::::

::::exercise{level=advanced}
$y(0)=y(1)=0$ と $\int_0^1y\,\dd x=A$ の下で、$J[y]=\tfrac12\int_0^1(y')^2\,\dd x$ を最小にせよ。
:::solution
$G=y$ に対し $\zeta=x(1-x)$ なら $\delta K[\zeta]=\int_0^1\zeta\,\dd x=1/6\ne0$ なので乗数を使える。
$F+\lambda G=(y')^2/2+\lambda y$ から $y''=\lambda$、両端条件から $y=\lambda(x^2-x)/2$ を得る。
積分条件は $-\lambda/12=A$ である。よって

$$
\lambda=-12A,\qquad y_*=6Ax(1-x),\qquad J[y_*]=6A^2.
$$

他の許容曲線を $y_*+\eta$ とおくと $\int_0^1\eta\,\dd x=0$ である。
交差項 $\int_0^1 y_*'\eta'\,\dd x=-\int_0^1y_*''\eta\,\dd x=-\lambda\int_0^1\eta\,\dd x$ は消え、差は $\tfrac12\int_0^1(\eta')^2\,\dd x$ となる。したがって最小解は一意である。
:::
::::

::::exercise{level=advanced}
終点が曲線 $y=\psi(x)$ 上を動く問題を考える。$y(a)=A$ を固定し、$J[y]=\int_a^{b}F(x,y,y')\,\dd x$ の終点 $b$ と $y(b)=\psi(b)$ を自由にする。$y$ が停留するならば、Euler–Lagrange 方程式に加えて終点で

$$
F+(\psi'-y')F_p=0
$$

（:term[横断性条件]{en="transversality condition" yomi="おうだんせいじょうけん"}）が成り立つことを示せ。さらに、点から曲線までの最短線は曲線と直交することを導け。$F$、$\psi$、$y$ は十分滑らかとし、$y$ は $b$ の少し先まで $C^2$ 級に延長しておく。
:::solution
まず終点を固定した変形から Euler–Lagrange 方程式が従う。次に終点を動かす族を作る。
$\chi$ を $C^2$ 級で $\chi(a)=0$、$\chi(b)=1$ を満たす関数とし、$y_\eps(x)=y(x)+c\chi(x)$、終点を $b+\eps$ とする。
終点条件 $y(b+\eps)+c\chi(b+\eps)=\psi(b+\eps)$ の左辺から右辺を引いたものは $(\eps,c)=(0,0)$ で零、$c$ 偏微分は $\chi(b)=1\ne0$ なので、陰関数定理により $c=c(\eps)$、$c(0)=0$ が定まる。$\eps=0$ で微分すると

$$
y'(b)+c'(0)=\psi'(b),\qquad c'(0)=\psi'(b)-y'(b).
$$

$j(\eps)=\int_a^{b+\eps}F(x,y_\eps,y_\eps')\,\dd x$ を微分すると、$\eta=c'(0)\chi$ として

$$
j'(0)=F\big|_{x=b}+\int_a^b(F_y\eta+F_p\eta')\,\dd x
=F\big|_{x=b}+F_p\big|_{x=b}\eta(b)+\int_a^b\Big(F_y-\dv{}{x}F_p\Big)\eta\,\dd x .
$$

積分は Euler–Lagrange 方程式で零、$\eta(b)=c'(0)=\psi'-y'$ なので、停留条件 $j'(0)=0$ は $F+(\psi'-y')F_p=0$ となる。

$F=\sqrt{1+p^2}$ なら $F_p=y'/\sqrt{1+(y')^2}$ で、条件に $\sqrt{1+(y')^2}$ を掛けると $1+(y')^2+(\psi'-y')y'=0$、すなわち $1+\psi'y'=0$ である。
これは直線の方向 $(1,y')$ と曲線の接線 $(1,\psi')$ の内積が零、つまり直交することを意味する。
:::
::::

::::exercise{level=advanced}
サイクロイド :eqref[eq-am-cycloid] に沿って、$\theta=\theta_0\in(0,\pi)$ の点から初速度零で滑らせる。最深点 $\theta=\pi$ に着くまでの時間が $\theta_0$ によらず $\pi\sqrt{a/g}$ であることを示せ（Huygens の等時性）。
:::solution
エネルギー保存から $v^2=2g(y-y_0)=2ga(\cos\theta_0-\cos\theta)$ である。$\dd s=2a\sin(\theta/2)\,\dd\theta$ と $\cos\theta=2\cos^2(\theta/2)-1$ を使い、$u=\cos(\theta/2)$、$u_0=\cos(\theta_0/2)$ とおくと $\dd u=-\tfrac12\sin(\theta/2)\,\dd\theta$ なので

$$
t=\int_{\theta_0}^\pi\frac{2a\sin(\theta/2)\,\dd\theta}{\sqrt{4ga(u_0^2-u^2)}}
=2\sqrt{\frac ag}\int_0^{u_0}\frac{\dd u}{\sqrt{u_0^2-u^2}}
=2\sqrt{\frac ag}\cdot\frac\pi2=\pi\sqrt{\frac ag}.
$$

最後の積分は $u=u_0\sin\sigma$ で $\int_0^{\pi/2}\dd\sigma$ となる。結果は $\theta_0$ を含まない。
:::
::::

::::exercise{level=advanced}
$J[y]=\int_0^1(y')^3\,\dd x$、$y(0)=0$、$y(1)=1$ とする。$y=x$ は弱い極小を与えるが強い極小を与えないことを示せ。強い極小の比較には、角をもつ区分的に $C^1$ 級の関数を用いてよい（角を少し丸めても $J$ の値はわずかしか変わらない）。
:::solution
弱い極小：$\eta$ を許容方向とすると $\int_0^1\eta'\,\dd x=0$ であり、

$$
J[x+\eta]=\int_0^1(1+\eta')^3\,\dd x=1+3\int_0^1\eta'\,\dd x+\int_0^1(\eta')^2(3+\eta')\,\dd x=1+\int_0^1(\eta')^2(3+\eta')\,\dd x .
$$

$\norm{\eta}_1<3$ なら $3+\eta'>0$ なので $J[x+\eta]\ge J[x]=1$ である。

強い極小でないこと：整数 $N\ge2$ に対し $\eps=1/N^2$ とし、$[0,\eps]$ では傾き $-N$、$[\eps,1]$ では傾き $s=(1+N\eps)/(1-\eps)$ の折れ線 $y_N$ を作る。端点条件を満たし、

$$
J[y_N]=-N^3\eps+s^3(1-\eps)=-N+\frac{(1+1/N)^3}{(1-1/N^2)^2}\to-\infty .
$$

一方 $\abs{y_N(x)-x}$ は $x=\eps$ で最大となり、その値は $(N+1)\eps\to0$ である。よって最大値ノルムで $y=x$ のいくらでも近くに $J$ の値が小さい関数がある。
差 $y_N-x$ の傾きは $-N-1$ と大きく、$C^1$ ノルムでは近くないことが二つの極小の違いである。
:::
::::

## 確認テスト

::::quiz
固定した両端の間で第一変分がすべての方向で零になった。何が結論できるか。

- [ ] 汎関数は必ず最小値をとる
- [x] 曲線は停留している
- [ ] 被積分関数が曲線の各点で零になる

:::explanation
停留は第一変分が零になることそのものである。極小かどうかの判定には、第二変分（:ref[prop-am-legendre]、:ref[exm-am-jacobi-example]）や凸性（:ref[thm-am-convex-sufficient]）、競合曲線との直接比較が必要になる。:ref[exm-am-weierstrass] のように最小値が存在しないこともある。
:::
::::

::::quiz{answer="2"}
$F=(y')^2/2+2y$ の滑らかな停留曲線で $y''$ はいくらか。

:::explanation
$F_p=y'$、$F_y=2$ から、$\dd F_p/\dd x=F_y$ は $y''=2$ となる。
:::
::::

::::quiz
Beltrami の恒等式から $F-y'F_p$ が一定となるための十分条件はどれか。曲線は $C^2$ 級の停留曲線とする。

- [ ] $F$ が $y'$ に陽に依存しない
- [ ] $F$ が $y$ に陽に依存しない
- [x] $F$ が $x$ に陽に依存しない

:::explanation
曲線に沿った微分は $F_x$ に等しいので、$F_x=0$ なら一定である。$F_y=0$ から得られる定数は $F_p$ である。
:::
::::

::::quiz
$F$ が $C^2$ 級で、$C^1$ 級の曲線 $y$ が固定端の汎関数の弱い極小を与える。必ず成り立つものをすべて選べ。

- [x] $F_p(x,y,y')$ は $x$ について $C^1$ 級で $\dd F_p/\dd x=F_y$
- [x] $F_{pp}(x,y,y')\ge0$
- [ ] $F$ は $(y,p)$ について凸である
- [ ] $y$ は強い極小も与える

:::explanation
弱い極小は停留曲線なので、:ref[thm-am-euler-lagrange] により一つ目が成り立つ。二つ目は :ref[prop-am-legendre] である。凸性は十分条件であって必要条件ではない。強い極小は弱い極小より強い要求であり、$\int(y')^3\,\dd x$ の演習のように成り立たないことがある。
:::
::::

::::quiz{answer="3.1416" tol="0.0005"}
$J[y]=\tfrac12\int_0^\ell((y')^2-y^2)\,\dd x$、$y(0)=y(\ell)=0$ で、$y=0$ が唯一の最小解となる $\ell$ の上限はいくらか。小数第 4 位まで答えよ。

:::explanation
:ref[exm-am-jacobi-example] により、$\ell<\pi$ なら Wirtinger の不等式から唯一の最小解、$\ell=\pi$ では $c\sin x$ も最小解、$\ell>\pi$ では極小でない。上限は $\pi\simeq3.1416$ である。
:::
::::

::::quiz{answer="1.571" tol="0.0005"}
最速降下線がちょうど終点で最深点に達するとき、水平距離と落差の比 $X/Y$ はいくらか。小数第 3 位まで答えよ。

:::explanation
$\theta_1=\pi$ を :eqref[eq-am-cycloid] に入れると $X=\pi a,\ Y=2a$ であり、比は $\pi/2\simeq1.571$ である。
:::
::::

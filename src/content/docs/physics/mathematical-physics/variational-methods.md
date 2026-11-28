---
title: "変分法"
description: "積分の停留条件と境界・拘束・安定性を導き、固有値を有限個の試行関数で近似する。"
notation:
  - symbol: '\delta J'
    meaning: "汎関数 J の第一変分"
  - symbol: '\delta^2J'
    meaning: "汎関数 J の第二変分"
  - symbol: '\mathcal R[u]'
    meaning: "Rayleigh 商"
---

## この章の目標

- 場や高階導関数を含む汎関数の停留方程式を導く。
- 拘束と自由な端点から積分乗数と境界条件を求める。
- 第二変分で停留解の安定性を判定する。
- Rayleigh–Ritz 法で固有値の上界を計算する。

## 多重積分の変分問題

### 場を近くの場に置き換える

滑らかな有界領域 $\Omega$ で $J[u]=\int_\Omega F(x,u,\grad u)\dd x$ とする。許される変分 $\eta$ に対して $u+\epsilon\eta$ を入れ、$\epsilon=0$ での微分を第一変分という。滑らかな $u,\eta,F$ について連鎖律と発散定理から
$$
\delta J[u;\eta]
=\int_\Omega(F_u-\divg F_{\grad u})\eta\,\dd x+
\int_{\partial\Omega}(F_{\grad u}\cdot n)\eta\,\dd S
\label{eq-mp-field-variation}
$$
である。内部に台を持つ $\eta$ を全て選べるなら、停留解は $F_u-\divg F_{\grad u}=0$ を満たす。これには :ref[lem-am-fundamental-variation] の議論と同じく、連続な係数が正または負の領域に非負の試験関数を取る。多成分の場では一度に一成分だけ $\eta$ を入れれば成分ごとの方程式を得る。

:::example[Poisson 方程式から最小化へ]
$F=|\grad u|^2/2-fu$ なら方程式は $-\Delta u=f$ である。Dirichlet を固定すれば $\eta=0$ が境界で成り立つ。解の周りで $J[u+\eta]-J[u]=\int|\grad\eta|^2/2$ であるから、連結領域で解は唯一の最小解である。例えば単位円板で $f=4$、外周で零なら $u=1-r^2$ であり、二次関数の Laplacian は $-4$ である。
:::

面積 $\int_\Omega\sqrt{1+|\grad u|^2}$ の停留グラフは $\divg(\grad u/\sqrt{1+|\grad u|^2})=0$ に従う。小傾斜 $|\grad u|\ll1$ では主項が $\Delta u=0$ となる。これは面積密度を $1+|\grad u|^2/2+O(|\grad u|^4)$ と近似することに対応する。

### 多成分の場と保存則

$u=(u^1,\ldots,u^m)$ に対して $p_i^\alpha=\partial_i u^\alpha$ と書く。停留条件は各 $\alpha$ で $F_{u^\alpha}-\sum_i\partial_i F_{p_i^\alpha}=0$ である。特に積分密度が $u^1$ に陽に依存しなければ、$j_i=F_{p_i^1}$ の発散は零となる。例えば時間 $t$ と空間 $x$ の場で
$$
F=\frac12(u_t^2-c^2u_x^2)
$$
を取れば、停留式は波動方程式 $u_{tt}-c^2u_{xx}=0$ であり、$j=(u_t,-c^2u_x)$ が発散零になる。この作用は鞍点を持つ。例えば $\eta(t,x)=\chi(t)\psi(x)\sin(Nt)$ を内部に台を持つ試験関数に取れば、大きい $N$ で $\int\eta_t^2$ が $\int c^2\eta_x^2$ を上回る。空間の振動を速くすれば二次の差は負になる。

### 厳密な二次の差として評価する

Dirichlet 問題で $J[u]=\int_\Omega(\frac12|\grad u|^2-fu)$ を用いる。解 $u_0$ と同じ境界値を持つ $u_0+\eta$ を代入して引くと、一次の差は $\int(\grad u_0\cdot\grad\eta-f\eta)=0$ になる。境界項は $\eta=0$ で消え、内部項は $-\Delta u_0=f$ で消える。そのため残る差は $\frac12\|\grad\eta\|_2^2$ であり、停留に加えて最小性も得られる。連結領域で差が零なら $\eta$ は定数となり、境界条件でその定数は零になる。

:::example[単位正方形での正弦解（1）]
右辺を正弦モードに選び、解とエネルギーを求める。

$$
m=1
$$

両辺で消える正弦の積を使う。

$$
\begin{aligned}
f(x,y)\big|_{m=1} &= \sin{\left(\pi x \right)} \sin{\left(\pi y \right)}
\end{aligned}
$$

二階微分で掛かる固有値の逆を掛ける。

$$
\begin{aligned}
u_0(x,y)\big|_{m=1} &= \frac{\sin{\left(\pi x \right)} \sin{\left(\pi y \right)}}{2 \pi^{2}}
\end{aligned}
$$

解でのエネルギーは部分積分により $-\int fu_0/2$ になる。

$$
\begin{aligned}
J[u_0]\big|_{m=1} &= - \frac{1}{16 \pi^{2}}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
f(x,y)\big|_{x=\frac{1}{3}} &= \frac{\sqrt{3} \sin{\left(\pi y \right)}}{2} \\
\frac{\partial^1 (f(x,y))}{\partial x^1} &= \pi \sin{\left(\pi y \right)} \cos{\left(\pi x \right)} \\
\left.\frac{\partial^1 (f(x,y))}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{\pi \sin{\left(\pi y \right)}}{2}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
u_0(x,y)\big|_{x=\frac{1}{3}} &= \frac{\sqrt{3} \sin{\left(\pi y \right)}}{4 \pi^{2}} \\
\frac{\partial^1 (u_0(x,y))}{\partial x^1} &= \frac{\sin{\left(\pi y \right)} \cos{\left(\pi x \right)}}{2 \pi} \\
\left.\frac{\partial^1 (u_0(x,y))}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{\sin{\left(\pi y \right)}}{4 \pi}
\end{aligned}
$$
:::

:::example[単位正方形での正弦解（2）]
右辺を正弦モードに選び、解とエネルギーを求める。

$$
m=2
$$

両辺で消える正弦の積を使う。

$$
\begin{aligned}
f(x,y)\big|_{m=2} &= \sin{\left(2 \pi x \right)} \sin{\left(\pi y \right)}
\end{aligned}
$$

二階微分で掛かる固有値の逆を掛ける。

$$
\begin{aligned}
u_0(x,y)\big|_{m=2} &= \frac{\sin{\left(2 \pi x \right)} \sin{\left(\pi y \right)}}{5 \pi^{2}}
\end{aligned}
$$

解でのエネルギーは部分積分により $-\int fu_0/2$ になる。

$$
\begin{aligned}
J[u_0]\big|_{m=2} &= - \frac{1}{40 \pi^{2}}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
f(x,y)\big|_{x=\frac{1}{3}} &= \frac{\sqrt{3} \sin{\left(\pi y \right)}}{2} \\
\frac{\partial^1 (f(x,y))}{\partial x^1} &= 2 \pi \sin{\left(\pi y \right)} \cos{\left(2 \pi x \right)} \\
\left.\frac{\partial^1 (f(x,y))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - \pi \sin{\left(\pi y \right)}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
u_0(x,y)\big|_{x=\frac{1}{3}} &= \frac{\sqrt{3} \sin{\left(\pi y \right)}}{10 \pi^{2}} \\
\frac{\partial^1 (u_0(x,y))}{\partial x^1} &= \frac{2 \sin{\left(\pi y \right)} \cos{\left(2 \pi x \right)}}{5 \pi} \\
\left.\frac{\partial^1 (u_0(x,y))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - \frac{\sin{\left(\pi y \right)}}{5 \pi}
\end{aligned}
$$
:::

:::example[単位正方形での正弦解（3）]
右辺を正弦モードに選び、解とエネルギーを求める。

$$
m=3
$$

両辺で消える正弦の積を使う。

$$
\begin{aligned}
f(x,y)\big|_{m=3} &= \sin{\left(3 \pi x \right)} \sin{\left(\pi y \right)}
\end{aligned}
$$

二階微分で掛かる固有値の逆を掛ける。

$$
\begin{aligned}
u_0(x,y)\big|_{m=3} &= \frac{\sin{\left(3 \pi x \right)} \sin{\left(\pi y \right)}}{10 \pi^{2}}
\end{aligned}
$$

解でのエネルギーは部分積分により $-\int fu_0/2$ になる。

$$
\begin{aligned}
J[u_0]\big|_{m=3} &= - \frac{1}{80 \pi^{2}}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
f(x,y)\big|_{x=\frac{1}{3}} &= 0 \\
\frac{\partial^1 (f(x,y))}{\partial x^1} &= 3 \pi \sin{\left(\pi y \right)} \cos{\left(3 \pi x \right)} \\
\left.\frac{\partial^1 (f(x,y))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - 3 \pi \sin{\left(\pi y \right)}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
u_0(x,y)\big|_{x=\frac{1}{3}} &= 0 \\
\frac{\partial^1 (u_0(x,y))}{\partial x^1} &= \frac{3 \sin{\left(\pi y \right)} \cos{\left(3 \pi x \right)}}{10 \pi} \\
\left.\frac{\partial^1 (u_0(x,y))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - \frac{3 \sin{\left(\pi y \right)}}{10 \pi}
\end{aligned}
$$
:::

:::example[単位正方形での正弦解（4）]
右辺を正弦モードに選び、解とエネルギーを求める。

$$
m=4
$$

両辺で消える正弦の積を使う。

$$
\begin{aligned}
f(x,y)\big|_{m=4} &= \sin{\left(4 \pi x \right)} \sin{\left(\pi y \right)}
\end{aligned}
$$

二階微分で掛かる固有値の逆を掛ける。

$$
\begin{aligned}
u_0(x,y)\big|_{m=4} &= \frac{\sin{\left(4 \pi x \right)} \sin{\left(\pi y \right)}}{17 \pi^{2}}
\end{aligned}
$$

解でのエネルギーは部分積分により $-\int fu_0/2$ になる。

$$
\begin{aligned}
J[u_0]\big|_{m=4} &= - \frac{1}{136 \pi^{2}}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
f(x,y)\big|_{x=\frac{1}{3}} &= - \frac{\sqrt{3} \sin{\left(\pi y \right)}}{2} \\
\frac{\partial^1 (f(x,y))}{\partial x^1} &= 4 \pi \sin{\left(\pi y \right)} \cos{\left(4 \pi x \right)} \\
\left.\frac{\partial^1 (f(x,y))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - 2 \pi \sin{\left(\pi y \right)}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
u_0(x,y)\big|_{x=\frac{1}{3}} &= - \frac{\sqrt{3} \sin{\left(\pi y \right)}}{34 \pi^{2}} \\
\frac{\partial^1 (u_0(x,y))}{\partial x^1} &= \frac{4 \sin{\left(\pi y \right)} \cos{\left(4 \pi x \right)}}{17 \pi} \\
\left.\frac{\partial^1 (u_0(x,y))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - \frac{2 \sin{\left(\pi y \right)}}{17 \pi}
\end{aligned}
$$
:::

### 面積汎関数と小傾斜の極限

グラフ $z=u(x,y)$ の面積は $J[u]=\int_\Omega\sqrt{1+|\grad u|^2}\dd x\dd y$ である。二つの接ベクトル $(1,0,u_x),(0,1,u_y)$ の外積の長さが平方根になる。値を境界で固定した変分 $\eta$ を入れ、$\epsilon$ で微分すると $\int\nabla u\cdot\nabla\eta/\sqrt{1+|\nabla u|^2}$ を得る。部分積分により内部の方程式は $\divg(\nabla u/\sqrt{1+|\nabla u|^2})=0$ である。例えば $u=ax+by+d$ では微分の内側が定数なので、固定境界の平面は停留する。

この平面は同じ境界値のグラフの中で最小でもある。$F(p)=(1+|p|^2)^{1/2}$ とすると、Hessian の二次形式は
$$
q^\mathsf T D^2F(p)q
=\frac{|q|^2(1+|p|^2)-(p\cdot q)^2}{(1+|p|^2)^{3/2}}
\ge\frac{|q|^2}{(1+|p|^2)^{3/2}}
$$
であり、$q\ne0$ なら正である。従って点ごとに $F(p+\nabla\eta)\ge F(p)+DF(p)\nabla\eta$ で、平面では右辺の追加分は境界で $\eta=0$ のため積分が零になる。$|\nabla u|\ll1$ の範囲で二次まで展開すると $F=1+\frac12|\nabla u|^2+O(|\nabla u|^4)$ となり、元の非線形方程式に対する小傾斜の近似として Laplace 方程式が出る。境界で傾斜が大きい箇所を含むと、この一様な近似は使えない。
:::example[傾いた平面での単位面積と曲率]
平面 $u=2x+y$ について、傾きのある向きとそれに直交する向きの第二変分の係数を比較する。

面積の被積分関数を計算する。単位正方形なので積分値は同じ数になる。

$$
\begin{aligned}
|\nabla u|^2 &= 5 \\
J[u] &= \sqrt{6}
\end{aligned}
$$

勾配の向きと垂直な向きに対して、単位ベクトルの Hessian の値を得る。

$$
\begin{aligned}
q_\parallel^\mathsf T D^2F q_\parallel &= \frac{\sqrt{6}}{36} \\
q_\perp^\mathsf T D^2F q_\perp &= \frac{\sqrt{6}}{6}
\end{aligned}
$$
:::

## 高階導関数を含む変分問題

### 何回も部分積分する

$J[y]=\int_a^b F(x,y,y',\ldots,y^{(m)})\dd x$ で $F_j=\partial F/\partial y^{(j)}$ とする。$\delta J=\int\sum_{j=0}^m F_j\eta^{(j)}$ を $j$ 回ずつ部分積分すると、内部係数は $\sum_{j=0}^m(-\partial_x)^jF_j$ になる。境界には $j\ge1$ について $\sum_{k=0}^{j-1}[(-1)^kF_j^{(k)}\eta^{(j-1-k)}]_a^b$ が残る。各段で一つ微分を移すごとに $-1$ が掛かる。この式は :ref[prop-am-euler-poisson] の一変数の議論を全階数でまとめたものである。

例えば $J[y]=\int_0^1[(y'')^2/2-y]\dd x$ で両端の $y,y'$ を零とすると $y''''=1$ である。$y=x^2(1-x)^2/24$ を代入すると $y''''=1$ かつ四つの境界条件を満たす。梁では $(y'')^2$ の係数に曲げ剛性 $EI$（$\mathrm{N\,m^2}$）を掛け、$y$ を $\mathrm m$ で測れば積分はエネルギーの単位を持つ。ここでの式は係数を無次元化している。

### 板の曲げと二重 Laplacian

平面領域の単純な曲げモデル $J[u]=\int_\Omega[B(\Delta u)^2/2-fu]\dd x$ で $B$ を正の定数とする。$\int(\Delta u)(\Delta\eta)$ に Green の恒等式を二回使うと、内部は $\int\eta\Delta^2u$、境界は $\int_{\partial\Omega}[(\Delta u)\partial_n\eta-\eta\partial_n(\Delta u)]$ になる。固定された板の端で $u,\partial_nu$ を指定すれば二つの変分の跡はともに零であり、内部条件は $B\Delta^2u=f$ である。例えば単位円板で $B=1$、$u=(1-r^2)^2$ は境界で $u=\partial_nu=0$、$\Delta u=-8+16r^2$、$\Delta^2u=64$ である。したがって一様源 $f=64$ に対する停留解になる。一般の材料では曲げエネルギーの形は構成則に依存する。

### 曲げと張力の境界項を合わせる

$J[y]=\int_a^b [B(y'')^2/2+T(y')^2/2-fy]\dd x$ とする。$B,T$ は正の定数である。$\delta J$ を部分積分すると内部は $B y^{(4)}-Ty''-f$ となり、端には $[B y''\eta'+(Ty'-B y''')\eta]_a^b$ が残る。自由端で二つの変分を別々に指定できるなら、二つの係数をそれぞれ零にする。単純支持端で変位だけを固定すれば $\eta=0$ となるため、二階微分だけが自然条件で零になる。

:::example[片持ち梁の端条件と先端の変位]
$0\le x\le L$ で $T=0$、$B$ と荷重 $f$ が定数、左端を固定し右端を自由とする。条件は $y(0)=y'(0)=y''(L)=y'''(L)=0$ である。$By^{(4)}=f$ を右から二回積分すると $y''(x)=f(L-x)^2/(2B)$、さらに左から二回積分すると
$$
y(x)=\frac{f}{24B}x^2(x^2-4Lx+6L^2)
$$
になる。先端の変位は $fL^4/(8B)$ である。分布荷重 $f$ の単位は $\mathrm{N\,m^{-1}}$、曲げ剛性 $B$ は $\mathrm{N\,m^2}$ なので、この値の単位は $\mathrm m$ である。
:::

:::example[単純支持の曲げモード（1）]
単位区間で曲げ剛性と張力を $1$ とする。

$$
m=1
$$

変位と曲率が両端で零になる試行関数を作る。

$$
\begin{aligned}
y(x)\big|_{m=1} &= \sin{\left(\pi x \right)}
\end{aligned}
$$

内部作用素を作用させると一つの固有値が掛かる。

$$
\begin{aligned}
(D^4-D^2)y\big|_{m=1} &= \pi^{2} \left(1 + \pi^{2}\right) \sin{\left(\pi x \right)}
\end{aligned}
$$

弾性エネルギーは曲率の二乗と傾斜の二乗の積分の和から得る。

$$
\begin{aligned}
J_{\rm elastic}\big|_{m=1} &= \frac{\pi^{2} \left(1 + \pi^{2}\right)}{4}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
y(x)\big|_{x=\frac{1}{3}} &= \frac{\sqrt{3}}{2} \\
\frac{\partial^1 (y(x))}{\partial x^1} &= \pi \cos{\left(\pi x \right)} \\
\left.\frac{\partial^1 (y(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{\pi}{2}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
(D^4-D^2)y\big|_{x=\frac{1}{3}} &= \frac{\sqrt{3} \pi^{2} \left(1 + \pi^{2}\right)}{2} \\
\frac{\partial^1 ((D^4-D^2)y)}{\partial x^1} &= \pi^{3} \left(1 + \pi^{2}\right) \cos{\left(\pi x \right)} \\
\left.\frac{\partial^1 ((D^4-D^2)y)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{\pi^{3} \left(1 + \pi^{2}\right)}{2}
\end{aligned}
$$
:::

:::example[単純支持の曲げモード（2）]
単位区間で曲げ剛性と張力を $1$ とする。

$$
m=2
$$

変位と曲率が両端で零になる試行関数を作る。

$$
\begin{aligned}
y(x)\big|_{m=2} &= \sin{\left(2 \pi x \right)}
\end{aligned}
$$

内部作用素を作用させると一つの固有値が掛かる。

$$
\begin{aligned}
(D^4-D^2)y\big|_{m=2} &= 4 \pi^{2} \left(1 + 4 \pi^{2}\right) \sin{\left(2 \pi x \right)}
\end{aligned}
$$

弾性エネルギーは曲率の二乗と傾斜の二乗の積分の和から得る。

$$
\begin{aligned}
J_{\rm elastic}\big|_{m=2} &= \pi^{2} + 4 \pi^{4}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
y(x)\big|_{x=\frac{1}{3}} &= \frac{\sqrt{3}}{2} \\
\frac{\partial^1 (y(x))}{\partial x^1} &= 2 \pi \cos{\left(2 \pi x \right)} \\
\left.\frac{\partial^1 (y(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - \pi
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
(D^4-D^2)y\big|_{x=\frac{1}{3}} &= \sqrt{3} \pi^{2} \left(2 + 8 \pi^{2}\right) \\
\frac{\partial^1 ((D^4-D^2)y)}{\partial x^1} &= 8 \pi^{3} \left(1 + 4 \pi^{2}\right) \cos{\left(2 \pi x \right)} \\
\left.\frac{\partial^1 ((D^4-D^2)y)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \pi^{3} \left(- 16 \pi^{2} - 4\right)
\end{aligned}
$$
:::

:::example[単純支持の曲げモード（3）]
単位区間で曲げ剛性と張力を $1$ とする。

$$
m=3
$$

変位と曲率が両端で零になる試行関数を作る。

$$
\begin{aligned}
y(x)\big|_{m=3} &= \sin{\left(3 \pi x \right)}
\end{aligned}
$$

内部作用素を作用させると一つの固有値が掛かる。

$$
\begin{aligned}
(D^4-D^2)y\big|_{m=3} &= 9 \pi^{2} \left(1 + 9 \pi^{2}\right) \sin{\left(3 \pi x \right)}
\end{aligned}
$$

弾性エネルギーは曲率の二乗と傾斜の二乗の積分の和から得る。

$$
\begin{aligned}
J_{\rm elastic}\big|_{m=3} &= \frac{9 \pi^{2} \left(1 + 9 \pi^{2}\right)}{4}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
y(x)\big|_{x=\frac{1}{3}} &= 0 \\
\frac{\partial^1 (y(x))}{\partial x^1} &= 3 \pi \cos{\left(3 \pi x \right)} \\
\left.\frac{\partial^1 (y(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - 3 \pi
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
(D^4-D^2)y\big|_{x=\frac{1}{3}} &= 0 \\
\frac{\partial^1 ((D^4-D^2)y)}{\partial x^1} &= 27 \pi^{3} \left(1 + 9 \pi^{2}\right) \cos{\left(3 \pi x \right)} \\
\left.\frac{\partial^1 ((D^4-D^2)y)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \pi^{3} \left(- 243 \pi^{2} - 27\right)
\end{aligned}
$$
:::

:::example[単純支持の曲げモード（4）]
単位区間で曲げ剛性と張力を $1$ とする。

$$
m=4
$$

変位と曲率が両端で零になる試行関数を作る。

$$
\begin{aligned}
y(x)\big|_{m=4} &= \sin{\left(4 \pi x \right)}
\end{aligned}
$$

内部作用素を作用させると一つの固有値が掛かる。

$$
\begin{aligned}
(D^4-D^2)y\big|_{m=4} &= 16 \pi^{2} \left(1 + 16 \pi^{2}\right) \sin{\left(4 \pi x \right)}
\end{aligned}
$$

弾性エネルギーは曲率の二乗と傾斜の二乗の積分の和から得る。

$$
\begin{aligned}
J_{\rm elastic}\big|_{m=4} &= \pi^{2} \left(4 + 64 \pi^{2}\right)
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
y(x)\big|_{x=\frac{1}{3}} &= - \frac{\sqrt{3}}{2} \\
\frac{\partial^1 (y(x))}{\partial x^1} &= 4 \pi \cos{\left(4 \pi x \right)} \\
\left.\frac{\partial^1 (y(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - 2 \pi
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
(D^4-D^2)y\big|_{x=\frac{1}{3}} &= \sqrt{3} \pi^{2} \left(- 128 \pi^{2} - 8\right) \\
\frac{\partial^1 ((D^4-D^2)y)}{\partial x^1} &= 64 \pi^{3} \left(1 + 16 \pi^{2}\right) \cos{\left(4 \pi x \right)} \\
\left.\frac{\partial^1 ((D^4-D^2)y)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \pi^{3} \left(- 512 \pi^{2} - 32\right)
\end{aligned}
$$
:::

### 梁の端条件の数を確かめる

単純支持の梁に一様荷重 $f$ がかかるとし、剛性を $B=1$、張力を零とする。方程式 $u^{(4)}=f$ は四階であり、左右の変位と左右の曲率を零にして四条件を課す。$f=24$ なら $u=x^4+Ax^3+Bx^2+Cx+D$ と書ける。左の条件で $B=D=0$、右の曲率で $A=-2$、右の変位で $C=1$ となる。すると $u=x^4-2x^3+x$ は中央で $5/16$、両端で零である。変分 $\eta$ が両端で零なら、交差項 $\int u''\eta''-24\eta$ は内部方程式と $u''=0$ の境界条件で消える。エネルギー差は $\frac12\int(\eta'')^2$ であり、差が零なら $\eta$ は一次になる。端点で零なのでこの一次関数は零である。
:::example[曲げモーメントの符号を読む]
上の単純支持梁では、境界項で $\eta'$ に掛かる量を $M=Bu''$ と取る。この符号規約で $M=12x(x-1)$ だから、両端では零、中央では $-3$ である。中央変位は $5/16$、傾きは零である。荷重と比較するにはモーメントをさらに二回微分し、$M''=24$ を得る。
:::

## 等周問題と拘束付き変分

### 許される変位の方向を限定する

$C[u]=c$ の拘束で $\delta C$ が零汎関数でないとする。$\delta C[\zeta]\ne0$ となる方向 $\zeta$ を取る。任意の $\eta$ に対し $C[u+s\eta+t(s)\zeta]=c$ を陰関数定理で解けば $t'(0)=-\delta C[\eta]/\delta C[\zeta]$ である。従って拘束付き停留性は $\delta J[\eta]=\lambda\delta C[\eta]$ を与え、$\lambda=\delta J[\zeta]/\delta C[\zeta]$ は $\eta$ によらない。これが等周問題で用いる定数の Lagrange 乗数である。複数拘束なら同じ議論を可逆な拘束微分行列を持つ方向の組で行う。拘束の微分が独立でないとこの導出は使えない。

例えば $y(0)=y(1)=0$、$\int_0^1 y=A$ を制約に $\int_0^1(y')^2/2$ を最小化する。$-y''=\lambda$ より $y=(\lambda/2)x(1-x)$ である。積分は $\lambda/12$ なので $\lambda=12A$ である。許される差 $\eta$ は $\int\eta=0$ で、一次変分は消え、残りが $\int(\eta')^2/2\ge0$ なので最小性も従う。

### 周長で囲まれた面積の上限

長さ $L$ の単純な閉曲線を弧長 $s$ で表示し、$x(s),y(s)$ の平均を平行移動で零にする。向き付き面積は Green の定理で $A=\tfrac12\int_0^L(xy'-yx')\dd s$ である。Cauchy–Schwarz と $x'^2+y'^2=1$ によって
$$
|A|\le\frac12\left(\int_0^L(x^2+y^2)\dd s\right)^{1/2}\sqrt L.
$$
平均零の周期関数 $f$ の Fourier 級数では定数項が無く、導関数の各重みは少なくとも $(2\pi/L)^2$ である。従って $\int f^2\le (L/(2\pi))^2\int(f')^2$ である。$f=x,y$ に適用すると $|A|\le L^2/(4\pi)$ となる。円では $r=L/(2\pi)$ で面積がその上限に一致する。等号の場合には Fourier の一次モードだけがあり、$x,y$ の組は $a\cos(2\pi s/L)+b\sin(2\pi s/L)$（$a,b$ は平面ベクトル）である。弧長条件を代入すると $|a|=|b|=L/(2\pi)$、$a\cdot b=0$ となるので円である。例えば $L=2\pi\,\mathrm m$ なら最大面積は $\pi\,\mathrm{m^2}$ である。この場合は停留条件だけでなく大域的な最適性も判定できる。

### 積分拘束を満たす最小関数

値を固定した端点の下で $\int y=C$ とする。この拘束の微分は $\int\eta$ である。停留条件がその核で成り立てば、任意の $\eta$ から単位積分を持つ試験関数の倍数を引いて核に入れられる。従って $J^\prime$ は拘束の微分の定数倍になる。次の二次汎関数では、拘束を保つ差について交差項が消えるので解は最小である。

:::example[面積を指定した放物線（1）]
単位区間で $y(0)=y(1)=0$、$\int y=m$、$J=\frac12\int (y^\prime)^2$ とする。

$$
m=1
$$

定数の二階微分を積分し、両端条件と面積から倍率を決める。

$$
\begin{aligned}
y(x)\big|_{m=1} &= 6 x \left(1 - x\right)
\end{aligned}
$$

乗数は $-y^{\prime\prime}$ として読む。これは $J-\lambda\int y$ の規約である。

$$
\begin{aligned}
\lambda\big|_{m=1} &= 12
\end{aligned}
$$

最小エネルギーは $\lambda\int y/2$ になる。

$$
\begin{aligned}
J[y]\big|_{m=1} &= 6
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
y(x)\big|_{x=\frac{1}{3}} &= \frac{4}{3} \\
\frac{\partial^1 (y(x))}{\partial x^1} &= 6 - 12 x \\
\left.\frac{\partial^1 (y(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= 2
\end{aligned}
$$
:::

:::example[面積を指定した放物線（2）]
単位区間で $y(0)=y(1)=0$、$\int y=m$、$J=\frac12\int (y^\prime)^2$ とする。

$$
m=2
$$

定数の二階微分を積分し、両端条件と面積から倍率を決める。

$$
\begin{aligned}
y(x)\big|_{m=2} &= 12 x \left(1 - x\right)
\end{aligned}
$$

乗数は $-y^{\prime\prime}$ として読む。これは $J-\lambda\int y$ の規約である。

$$
\begin{aligned}
\lambda\big|_{m=2} &= 24
\end{aligned}
$$

最小エネルギーは $\lambda\int y/2$ になる。

$$
\begin{aligned}
J[y]\big|_{m=2} &= 24
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
y(x)\big|_{x=\frac{1}{3}} &= \frac{8}{3} \\
\frac{\partial^1 (y(x))}{\partial x^1} &= 12 - 24 x \\
\left.\frac{\partial^1 (y(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= 4
\end{aligned}
$$
:::

:::example[面積を指定した放物線（3）]
単位区間で $y(0)=y(1)=0$、$\int y=m$、$J=\frac12\int (y^\prime)^2$ とする。

$$
m=3
$$

定数の二階微分を積分し、両端条件と面積から倍率を決める。

$$
\begin{aligned}
y(x)\big|_{m=3} &= 18 x \left(1 - x\right)
\end{aligned}
$$

乗数は $-y^{\prime\prime}$ として読む。これは $J-\lambda\int y$ の規約である。

$$
\begin{aligned}
\lambda\big|_{m=3} &= 36
\end{aligned}
$$

最小エネルギーは $\lambda\int y/2$ になる。

$$
\begin{aligned}
J[y]\big|_{m=3} &= 54
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
y(x)\big|_{x=\frac{1}{3}} &= 4 \\
\frac{\partial^1 (y(x))}{\partial x^1} &= 18 - 36 x \\
\left.\frac{\partial^1 (y(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= 6
\end{aligned}
$$
:::

:::example[面積を指定した放物線（4）]
単位区間で $y(0)=y(1)=0$、$\int y=m$、$J=\frac12\int (y^\prime)^2$ とする。

$$
m=4
$$

定数の二階微分を積分し、両端条件と面積から倍率を決める。

$$
\begin{aligned}
y(x)\big|_{m=4} &= 24 x \left(1 - x\right)
\end{aligned}
$$

乗数は $-y^{\prime\prime}$ として読む。これは $J-\lambda\int y$ の規約である。

$$
\begin{aligned}
\lambda\big|_{m=4} &= 48
\end{aligned}
$$

最小エネルギーは $\lambda\int y/2$ になる。

$$
\begin{aligned}
J[y]\big|_{m=4} &= 96
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
y(x)\big|_{x=\frac{1}{3}} &= \frac{16}{3} \\
\frac{\partial^1 (y(x))}{\partial x^1} &= 24 - 48 x \\
\left.\frac{\partial^1 (y(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= 8
\end{aligned}
$$
:::

### ベクトル値の拘束と独立性

拘束が $C_1[u],C_2[u]$ の二つで、線形化した写像 $\eta\mapsto(C_1'[u]\eta,C_2'[u]\eta)$ が全射とする。$\eta_1,\eta_2$ をそれぞれ標準基底へ送る変分として取れる。任意の $\eta$ から $(C_1'[u]\eta)\eta_1+(C_2'[u]\eta)\eta_2$ を引くと、線形化された拘束の核に入る。線形な拘束ならこの核の変分は有限の大きさでも許される。非線形で $C^1$ の拘束の場合には、$u+\epsilon\eta+s_1\eta_1+s_2\eta_2$ を使い、二つの $s_j$ に関する導関数が単位行列なので陰関数定理により拘束を保つ曲線を作れる。核に沿っては $s_j'(0)=0$ である。従って停留条件は $J'[u]=\lambda_1 C_1'[u]+\lambda_2 C_2'[u]$ になる。拘束が従属ならこの手順で二つの独立な乗数は定まらない。

例えば単位区間の両端を固定し、$\int y=1$ かつ $\int xy=1/4$ を課す。作用 $\frac12\int(y')^2$ の停留解は $-y''=\lambda_1+\lambda_2 x$ を満たす。境界値を先に合わせると $y=\lambda_1 x(1-x)/2+\lambda_2 x(1-x^2)/6$ である。この二基底の二つのモーメントを積分して連立系を解けばよい。差の二乗エネルギーが非負になるため、一解が得られれば一意の最小解である。
:::example[二つのモーメントを固定したエネルギー]
右辺の二つの係数を別々の基底にし、線形な拘束を合わせる。

二つの Dirichlet 解を並べる。

$$
\begin{aligned}
\phi_1 &= \frac{x \left(1 - x\right)}{2} \\
\phi_2 &= \frac{x \left(1 - x^{2}\right)}{6}
\end{aligned}
$$

それぞれのモーメントを列に並べる。

$$
\begin{aligned}
B &= \begin{pmatrix}\frac{1}{12} & \frac{1}{24}\\\frac{1}{24} & \frac{1}{45}\end{pmatrix}
\end{aligned}
$$

解くと、二つの乗数と試行関数が求まる。

$$
\begin{aligned}
(\lambda_1,\lambda_2)^{\mathsf T} &= \begin{pmatrix}102\\-180\end{pmatrix} \\
y &= 51 x \left(1 - x\right) - 30 x \left(1 - x^{2}\right)
\end{aligned}
$$

結果の面積と一次モーメントを積分し直す。

$$
\begin{aligned}
\int_0^1 y\,\dd x &= 1 \\
\int_0^1xy\,\dd x &= \frac{1}{4}
\end{aligned}
$$
:::

### 点ごとの拘束

積分値を固定する拘束に対しては定数の乗数が現れる。これに対し、各点で場の値を制限すると乗数は位置に依存し得る。例えば $u:\Omega\to\R^m$ に $|u(x)|^2=1$ を課し、固定境界で $J=\frac12\int|\grad u|^2$ を停留させる。接方向 $u\cdot\eta=0$ に対して $(u+\epsilon\eta)/|u+\epsilon\eta|$ を取れば拘束を保ち、その微分は $\eta$ である。従って $-\Delta u$ の接成分が零となり、$-\Delta u=\lambda(x)u$ と書ける。$|u|^2=1$ を二回微分すると $u\cdot\Delta u=-|\grad u|^2$ だから、$\lambda=|\grad u|^2$ である。例えば区間上の $u(x)=(\cos ax,\sin ax)$ では $-u''=a^2 u$ であり、乗数も $a^2$ になる。

## 自然境界条件と横断条件

### 境界を動かすと何が残るか

式 :eqref[eq-mp-field-variation] で境界値が自由なら $F_{\grad u}\cdot n=0$ が必要であり、これを :term[自然境界条件]{en="natural boundary condition" yomi="しぜんきょうかいじょうけん"} という。境界値が自由な部分で任意の滑らかな跡を持つ変分は、法線方向にカットオフを掛けて構成できる。例えば $F=|\grad u|^2/2-fu$ なら $\partial_nu=0$ である。高階の前節の式で $m=2$ なら境界項は $[(F_1-F_2')\eta+F_2\eta']$ であり、$\eta,\eta'$ を独立に変えられる端点で $F_2=F_1-F_2'=0$ となる。

一変数の積分で上端が動き、端点の全変化を $(\delta b,\delta y_b)$ とする。固定座標での変分は $\eta(b)=\delta y_b-y'(b)\delta b$ であり、積分区間の変化は $F(b)\delta b$ を加える。従って停留解の端点項は $(F-y'F_{y'})\delta b+F_{y'}\delta y_b$ になる。端点が曲線 $y_b=\gamma(b)$ に沿うなら $\delta y_b=\gamma'(b)\delta b$ なので $F+(\gamma'-y')F_{y'}=0$ である。長さ $F=\sqrt{1+(y')^2}$ では $1+\gamma' y'=0$ となり、二つの接線は直交する。例えば $(0,0)$ から $y=2-x$ までの最短直線は方向 $(1,1)$ で、交点は $(1,1)$、長さは $\sqrt2$ である。

:::example[境界にもエネルギーがある場合]
体積の積分密度 $|\grad u|^2/2-fu$ に境界積分 $\int_{\partial\Omega}(\alpha u^2/2-gu)\dd S$ を加える。境界で $u$ が自由なら、境界項は $(\partial_nu+\alpha u-g)\eta$ になる。従って $\partial_nu+\alpha u=g$ が自然に出る。$\alpha\ge0$ のとき停留解からのエネルギー差は $\frac12\int_\Omega|\grad\eta|^2+\frac12\int_{\partial\Omega}\alpha\eta^2$ であり、連結領域で $\alpha>0$ の部分の面積が正なら零方向はない。
:::

区分的に滑らかな曲線の折点 $(x_0,y_0)$ を動かせる場合には、区間をその点で分けて同じ端点変分を使う。左側では上端、右側では下端なので、二つの端点項には逆符号が付く。$\delta x_0,\delta y_0$ が独立なら、$F_{y'}$ と $F-y'F_{y'}$ は折点で左右の値が一致する。これを :term[ワイエルシュトラス・エルドマンの角条件]{en="Weierstrass–Erdmann corner conditions" yomi="わいえるしゅとらすえるどまんのかどじょうけん"} という。例えば $F=(y')^2/2$ では第一条件だけで傾きの跳びが消える。

### 動く端点の速度と横断条件

終点 $(b,y(b))$ を曲線 $Y=\psi(X)$ の上で動かす。変分の端での値は $\eta(b)=(\psi^\prime(b)-y^\prime(b))\delta b$ であり、積分区間の変化からさらに $F(b)\delta b$ が出る。従って横断条件は $F+F_{y^\prime}(\psi^\prime-y^\prime)=0$ である。長さの場合は、直線解の接線と端点曲線の接線の内積が零という条件になる。

:::example[傾いた直線への最近点（1）]
原点から $Y=m-X$ への距離を最小にする。

$$
m=1
$$

終点を $(x,m-x)$ とし、距離の二乗を作る。

$$
\begin{aligned}
d^2(x)\big|_{m=1} &= x^{2} + \left(x - 1\right)^{2}
\end{aligned}
$$

終点変数を微分して停留条件を解く。

$$
\begin{aligned}
x_*\big|_{m=1} &= \frac{1}{2}
\end{aligned}
$$

この二次式は凸なので距離の最小値が求まる。

$$
\begin{aligned}
d(x_*)\big|_{m=1} &= \frac{\sqrt{2}}{2}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
d^2(x)\big|_{x=\frac{1}{3}} &= \frac{5}{9} \\
\frac{\partial^1 (d^2(x))}{\partial x^1} &= 4 x - 2 \\
\left.\frac{\partial^1 (d^2(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - \frac{2}{3}
\end{aligned}
$$
:::

:::example[傾いた直線への最近点（2）]
原点から $Y=m-X$ への距離を最小にする。

$$
m=2
$$

終点を $(x,m-x)$ とし、距離の二乗を作る。

$$
\begin{aligned}
d^2(x)\big|_{m=2} &= x^{2} + \left(x - 2\right)^{2}
\end{aligned}
$$

終点変数を微分して停留条件を解く。

$$
\begin{aligned}
x_*\big|_{m=2} &= 1
\end{aligned}
$$

この二次式は凸なので距離の最小値が求まる。

$$
\begin{aligned}
d(x_*)\big|_{m=2} &= \sqrt{2}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
d^2(x)\big|_{x=\frac{1}{3}} &= \frac{26}{9} \\
\frac{\partial^1 (d^2(x))}{\partial x^1} &= 4 x - 4 \\
\left.\frac{\partial^1 (d^2(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - \frac{8}{3}
\end{aligned}
$$
:::

:::example[傾いた直線への最近点（3）]
原点から $Y=m-X$ への距離を最小にする。

$$
m=3
$$

終点を $(x,m-x)$ とし、距離の二乗を作る。

$$
\begin{aligned}
d^2(x)\big|_{m=3} &= x^{2} + \left(x - 3\right)^{2}
\end{aligned}
$$

終点変数を微分して停留条件を解く。

$$
\begin{aligned}
x_*\big|_{m=3} &= \frac{3}{2}
\end{aligned}
$$

この二次式は凸なので距離の最小値が求まる。

$$
\begin{aligned}
d(x_*)\big|_{m=3} &= \frac{3 \sqrt{2}}{2}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
d^2(x)\big|_{x=\frac{1}{3}} &= \frac{65}{9} \\
\frac{\partial^1 (d^2(x))}{\partial x^1} &= 4 x - 6 \\
\left.\frac{\partial^1 (d^2(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - \frac{14}{3}
\end{aligned}
$$
:::

:::example[傾いた直線への最近点（4）]
原点から $Y=m-X$ への距離を最小にする。

$$
m=4
$$

終点を $(x,m-x)$ とし、距離の二乗を作る。

$$
\begin{aligned}
d^2(x)\big|_{m=4} &= x^{2} + \left(x - 4\right)^{2}
\end{aligned}
$$

終点変数を微分して停留条件を解く。

$$
\begin{aligned}
x_*\big|_{m=4} &= 2
\end{aligned}
$$

この二次式は凸なので距離の最小値が求まる。

$$
\begin{aligned}
d(x_*)\big|_{m=4} &= 2 \sqrt{2}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
d^2(x)\big|_{x=\frac{1}{3}} &= \frac{122}{9} \\
\frac{\partial^1 (d^2(x))}{\partial x^1} &= 4 x - 8 \\
\left.\frac{\partial^1 (d^2(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - \frac{20}{3}
\end{aligned}
$$
:::

### 端点の移動がない場合とのつながり

固定端なら $\delta b=\eta(b)=0$ であり、境界項から新しい条件は出ない。位置の $x$ 成分を固定して高さだけ自由にすれば $\delta b=0$、$\eta(b)$ は任意なので $F_{y'}(b)=0$ になる。これに対して高さそのものを固定して $b$ を動かすと、横断条件には $\psi'=0$ を入れて $F-y'F_{y'}=0$ となる。例えば $F=\sqrt{1+(y')^2}$ ではこの最後の差は $1/\sqrt{1+(y')^2}>0$ であり、有限傾斜のグラフでは満たせない。幾何学的には水平な終点曲線へ直交する直線が垂直になるためである。従って端点を自由にした問題では、許容関数の表現で垂直な線分を扱えるかも確認する必要がある。

## 第 2 変分と Legendre 条件

### 停留解の近くの二次形式

$F(x,y,y')$ の第二変分は
$$
\delta^2J[u;\eta]=\int_a^b
\bigl(F_{yy}\eta^2+2F_{yy'}\eta\eta'+F_{y'y'}(\eta')^2\bigr)\dd x
$$
である。これは $\epsilon$ の二階微分であり、Taylor の二次項の係数はその半分である。弱い局所最小ではこの二次形式は全変分で非負となる。とくに $F_{y'y'}\ge0$ が Legendre の必要条件である。もし内部の一点で負なら、その周りの小区間で負の上限を取る。台をその区間に収めて $\eta=\psi(x)\sin(Nx)$ と置くと、最高階の項は $-cN^2\int\psi^2\cos^2(Nx)$ 以下、他の項は $O(N)$ 以下である。振動積分の消失により $\int\psi^2\cos^2(Nx)\to\int\psi^2/2>0$ なので矛盾する。

この条件だけでは最小性は出ない。$J=\frac12\int_0^1((u')^2-\alpha u^2)$ の零解では、$F_{y'y'}=1$ であるが、$\eta=\sin\pi x$ に対して第二変分は $(\pi^2-\alpha)/2$ である。$\alpha>\pi^2$ なら負になる。逆に $\alpha<\pi^2$ では正弦基底で各係数が正で、零解は最小である。この例では作用そのものが二次形式なので高次の剰余は無い。

### 第二変分を Jacobi 方程式で調べる

固定端で交差項を部分積分すると第二変分は $\int(P(\eta')^2+Q\eta^2)$、$P=F_{y'y'}$、$Q=F_{yy}-(F_{yy'})'$ である。$P>0$ とし、区間全体で正の解 $v$ が $-(Pv')'+Qv=0$ を満たすとする。$\eta=v h$ と書くと
$$
P(\eta')^2+Q\eta^2=Pv^2(h')^2+(Pvv'h^2)' .
$$
右辺を微分して交差項を比較し、$(Pv')'=Qv$ を使えば分かる。固定端では境界項が零であり、第二変分は非負、非零の $\eta$ なら正である。例えば $P=1,Q=-\alpha$、$[0,1]$、$0\le\alpha<\pi^2$ なら $v=\cos[\sqrt\alpha(x-1/2)]>0$ と取れる。$\alpha=\pi^2$ ではこの $v$ が端で零になり、実際 $\eta=\sin(\pi x)$ が第二変分の零方向になる。端での正性の仮定が効いている。

### 第一変分が零でも最小とは限らない

$J[y]=\frac12\int[(y^\prime)^2-\mu y^2]$ とし、両端で零にする。零解はどの $\mu$ でも停留する。最高階微分の係数は正なので弱い Legendre 条件を満たすが、低い周波数では $-\mu y^2$ の影響が大きい。正弦基底へ分解すると二次形式は $(m^2\pi^2-\mu)$ の符号で決まる。従って微分項の局所条件に加え、領域全体での低いモードを調べる必要がある。

:::example[負になるモードと零になるモード（1）]
単位区間で $\mu=4\pi^2$ とする。

$$
m=1
$$

端点で零になる変分を選ぶ。

$$
\begin{aligned}
\eta(x)\big|_{m=1} &= \sin{\left(\pi x \right)}
\end{aligned}
$$

二乗積分を評価する。

$$
\begin{aligned}
\int_0^1\eta^2\dd x &= \frac{1}{2}
\end{aligned}
$$

第二変分は各モードで次の値になる。

$$
\begin{aligned}
\delta^2 J[\eta]\big|_{m=1} &= - \frac{3 \pi^{2}}{2}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
\eta(x)\big|_{x=\frac{1}{3}} &= \frac{\sqrt{3}}{2} \\
\frac{\partial^1 (\eta(x))}{\partial x^1} &= \pi \cos{\left(\pi x \right)} \\
\left.\frac{\partial^1 (\eta(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{\pi}{2}
\end{aligned}
$$
:::

:::example[負になるモードと零になるモード（2）]
単位区間で $\mu=4\pi^2$ とする。

$$
m=2
$$

端点で零になる変分を選ぶ。

$$
\begin{aligned}
\eta(x)\big|_{m=2} &= \sin{\left(2 \pi x \right)}
\end{aligned}
$$

二乗積分を評価する。

$$
\begin{aligned}
\int_0^1\eta^2\dd x &= \frac{1}{2}
\end{aligned}
$$

第二変分は各モードで次の値になる。

$$
\begin{aligned}
\delta^2 J[\eta]\big|_{m=2} &= 0
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
\eta(x)\big|_{x=\frac{1}{3}} &= \frac{\sqrt{3}}{2} \\
\frac{\partial^1 (\eta(x))}{\partial x^1} &= 2 \pi \cos{\left(2 \pi x \right)} \\
\left.\frac{\partial^1 (\eta(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - \pi
\end{aligned}
$$
:::

:::example[負になるモードと零になるモード（3）]
単位区間で $\mu=4\pi^2$ とする。

$$
m=3
$$

端点で零になる変分を選ぶ。

$$
\begin{aligned}
\eta(x)\big|_{m=3} &= \sin{\left(3 \pi x \right)}
\end{aligned}
$$

二乗積分を評価する。

$$
\begin{aligned}
\int_0^1\eta^2\dd x &= \frac{1}{2}
\end{aligned}
$$

第二変分は各モードで次の値になる。

$$
\begin{aligned}
\delta^2 J[\eta]\big|_{m=3} &= \frac{5 \pi^{2}}{2}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
\eta(x)\big|_{x=\frac{1}{3}} &= 0 \\
\frac{\partial^1 (\eta(x))}{\partial x^1} &= 3 \pi \cos{\left(3 \pi x \right)} \\
\left.\frac{\partial^1 (\eta(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - 3 \pi
\end{aligned}
$$
:::

:::example[負になるモードと零になるモード（4）]
単位区間で $\mu=4\pi^2$ とする。

$$
m=4
$$

端点で零になる変分を選ぶ。

$$
\begin{aligned}
\eta(x)\big|_{m=4} &= \sin{\left(4 \pi x \right)}
\end{aligned}
$$

二乗積分を評価する。

$$
\begin{aligned}
\int_0^1\eta^2\dd x &= \frac{1}{2}
\end{aligned}
$$

第二変分は各モードで次の値になる。

$$
\begin{aligned}
\delta^2 J[\eta]\big|_{m=4} &= 6 \pi^{2}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
\eta(x)\big|_{x=\frac{1}{3}} &= - \frac{\sqrt{3}}{2} \\
\frac{\partial^1 (\eta(x))}{\partial x^1} &= 4 \pi \cos{\left(4 \pi x \right)} \\
\left.\frac{\partial^1 (\eta(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - 2 \pi
\end{aligned}
$$
:::

### 正の Jacobi 解による二次形式の判定

:::proposition[正の Jacobi 解と二次形式]{#prop-mp-jacobi-positive}
$P>0$ を $C^1$、$Q$ を連続な実関数とし、閉区間 $[a,b]$ で正の解 $z$ が $-(Pz')'+Qz=0$ を満たすとする。両端が零の非零な実変分 $\eta\in C^1([a,b])$ に対し、$\int(P(\eta')^2+Q\eta^2)>0$ である。
:::

:::proof
$\eta=zv$ として被積分関数を展開する。$P(\eta')^2=Pz^2(v')^2+2Pzz'vv'+P(z')^2v^2$ である。一方 $(Pzz')'=P(z')^2+z(Pz')'=P(z')^2+Qz^2$ なので、後ろの二項と $Q\eta^2$ を足したものは $(Pzz'v^2)'$ である。$z>0$ だから境界の $v$ は零であり、この全微分の積分は消える。残った $\int Pz^2(v')^2$ は非負で、零なら $v$ は定数になる。境界値でその定数は零となり、$\eta\ne0$ に反する。
:::

例えば $P=1,Q=-\omega^2$ を $[0,1]$ で考える。$\omega<\pi$ なら $z(x)=\cos(\omega(x-1/2))$ は両端を含めて正である。$z$ は二回微分で $-\omega^2 z$ になるので、上の仮定を満たす。$\omega=\pi$ では $z=\sin(\pi x)$ が原点および右端で零になり、$\eta=z$ と取ると二次形式も零になる。さらに $\omega>\pi$ で同じ $\eta=\sin(\pi x)$ を取ると値は $(\pi^2-\omega^2)/2<0$ になる。したがって $P>0$ のままでも、区間全体で変分を抑えられるかどうかが変化する。

### 共役点と弱い局所最小

$P>0$ とし、$-(Pv')'+Qv=0$ の解を $v(a)=0$、$P(a)v'(a)=1$ で定める。$a$ より右で最初に $v(\xi)=0$ となる点を $a$ に対する :term[共役点]{en="conjugate point" yomi="きょうやくてん"} という。$a<\xi<b$ なら、$\eta_0=v$ を $[a,\xi]$ に置き、$[\xi,b]$ で零とする。この連続な区分的滑らかな変分では二次形式は端点項 $[Pvv']_a^\xi=0$ に等しい。しかし二次形式の双線形部分は、固定端の滑らかな $\psi$ に対して
$$
B(\eta_0,\psi)=\int_a^\xi(Pv'\psi'+Qv\psi)\dd x=P(\xi)v'(\xi)\psi(\xi)
$$
である。解の一意性から $v'(\xi)\ne0$ なので、$\psi(\xi)\ne0$ を選べば $B$ は零でない。そこで $\eta_0+s\psi$ の二次形式は、小さい $s$ の符号を選ぶと負になる。固定端を保って折点の近くを滑らかにすれば、$H^1$ ノルムの近似と二次形式の連続性から負の滑らかな変分も得られる。従って弱い局所最小には区間内部に共役点がないことが必要である。$\xi=b$ なら $\eta=v$ が零方向になる。

逆に $v$ が $(a,b]$ で正とする。初期値 $w(a)=1$、$P(a)w'(a)=0$ の Jacobi 解を取り、$z=v+\epsilon w$ と置く。$a$ の近くでは $w>0$、そこから離れた閉区間では $v$ に正の下限があるので、十分小さい $\epsilon>0$ で $z$ は閉区間全体で正である。:ref[prop-mp-jacobi-positive] の変形で、第二変分は $\int Pz^2[(\eta/z)']^2$ に等しい。固定端の Poincaré 不等式と $z,z',1/z$ の有界性により、この式は $c\|\eta\|_{H^1}^2$（$c>0$）以上である。さらに $F$ が $C^2$ なら、$\|\eta\|_{C^1}\to0$ で Taylor の剰余は $o(\|\eta\|_{H^1}^2)$ になる。実際、剰余の係数は $F$ の二階微分の一様連続性で一様に零へ近づく。従って上の条件は固定端の厳密な弱局所最小を与える。

### 強い局所最小と Weierstrass の余剰関数

関数と導関数をともに小さく変える弱い局所最小に対し、強い局所最小では関数の一様ノルムだけを小さくした曲線と比較する。従って傾きの大きな変化を短区間に閉じ込めた変分も許される。固定端の $C^2$ 停留解の傾きを $p=y'(x_0)$ とする。強い局所最小には任意の許容傾き $q$ に対して
$$
E(x_0,y(x_0),p,q)
=F(x_0,y(x_0),q)-F(x_0,y(x_0),p)
-(q-p)F_{y'}(x_0,y(x_0),p)\ge0
$$
が必要である。左辺を :term[ワイエルシュトラスの余剰関数]{en="Weierstrass excess function" yomi="わいえるしゅとらすのよじょうかんすう"} と呼ぶ。

実際、内部の $x_0$ から長さ $h$ だけ傾きを $q$ にして、次の固定幅の区間で元の曲線へ滑らかに戻す。曲線の差 $\eta_h$ の大きさは $O(h)$、戻す区間での導関数の差も $O(h)$ である。$\delta J[y;\eta_h]=0$ は区分的に滑らかな変分にも部分積分で成り立つ。そこで $J[y+\eta_h]-J[y]$ からこの零の一次変分を引く。戻す区間の剰余は Taylor 展開で $O(h^2)$、短い区間の剰余は連続性で $hE+o(h)$ となる。従って $E<0$ なら、一様には任意に小さい変形でも $J$ が減る。滑らかな曲線だけを比較する場合も、折点の近くを長さ $o(h)$ の範囲で、傾きを有界に保って丸めれば同じ極限を得る。$q$ を $p$ に近づけると $E=(q-p)^2F_{y'y'}/2+o((q-p)^2)$ なので、Legendre の条件も従う。

:::example[弱い最小だが強い最小でない場合]
$F=(y')^2-(y')^4$ で、両端を零に固定する。$y=0$ の周りで $\|y'\|_\infty<1$ なら被積分関数は非負であり、零になる曲線は $y=0$ だけである。従って弱い意味では厳密な最小である。しかし $p=0$ では $E(0,q)=q^2-q^4$ だから、例えば $q=2$ で $-12<0$ となる。短い区間だけ傾きを変えれば、曲線の一様な差を小さく保ったまま作用を下げられる。
:::

## 固有値問題の変分原理と Rayleigh–Ritz 法

### 有限次元の行列に圧縮する

正則 Dirichlet 問題の Rayleigh 商と最小最大原理 :eqref[eq-mp-sl-minimax] を使う。試行関数の空間は、両端の値が零で二乗可積分な弱導関数を持つ $H^1_0(a,b)$ とする。区間の内部にコンパクトな台を持つ滑らかな関数は作用素の定義域に含まれ、$H^1_0$ で稠密である。従って最小最大原理の直交性の議論は、内積と形式の連続性によってこの空間にも移る。具体的には、最初の $j-1$ 固有関数に直交する $v$ を作用素の定義域から近似する際、近似からその $j-1$ 成分を引けば直交性を保てる。

試行基底を $v_1,\ldots,v_N$ とすると、二次形式とノルムはそれぞれ行列 $A_{ij}=\int (p v_j'\overline{v_i'}+qv_j\overline{v_i})$、$B_{ij}=\int wv_j\overline{v_i}$ になる。$B$ は独立な試行関数に対して正定値である。商を係数で停留させると $Ac=\lambda Bc$ となる。これは $B^{-1/2}AB^{-1/2}$ の通常の固有値問題であり、最小最大原理の部分空間を小さくしたので第 $j$ Ritz 値は真の第 $j$ 固有値以上である。試行空間を包含で広げれば上界は非増加である。Robin 条件 $p\partial_nu+\alpha u=0$ を課す端では、試行関数の端値を自由にし、分子に $\alpha |u|^2$ をその端で加える。部分積分するとこの項が端での流束と相殺する。

$-u''$ の単位区間 Dirichlet で $v_1=x(1-x),v_2=x(1-x)(2x-1)$ を使う。中点について一方が偶、他方が奇なので二行列の非対角は零である。$v_1$ の商は $10$、$v_2$ の分子は $1/5$、分母は $1/210$ なので商は $42$ である。第一値は第一の試行空間から改善せず、偶の近似には別の偶関数を足す必要がある。

### 近似が固有値へ収束する条件

固定した $j$ に対し、試行空間の増大列がエネルギー形式のノルムで $u_1,\ldots,u_j$ を近似できるとする。例えば $q+Cw>0$ となる $C$ を足した分子の平方根をこのノルムに取れる。$v_i^{(N)}\to u_i$ なら各内積の連続性により、$j$ 個の試行関数の Gram 行列は $I_j$ へ、形式の行列は $\diag(\lambda_1,\ldots,\lambda_j)$ へ収束する。行列の次元は固定なので係数の単位球上で二次形式は一様収束する。分母の最小固有値は $1$ に近づいて零から離れる。従ってこの $j$ 次元空間上の最大商は $\lambda_j$ に収束し、上界と最小最大原理で第 $j$ Ritz 値もそこへ収束する。逆に例えば中点について偶の関数だけを選ぶと、偶関数の固有値しか近似できず、全試行空間での第 $2$ 固有値に関する収束は出ない。

### 有限次元の商と質量行列

試行基底 $\phi_j$ に対し $A_{ij}=a(\phi_j,\phi_i)$、$B_{ij}=\langle\phi_j,\phi_i\rangle$ とする。有限次元の商 $c^*Ac/(c^*Bc)$ の停留条件は $Ac=\lambda Bc$ である。基底を非直交にしても $B$ を含めれば同じ試行空間を記述できる。$B$ は正定値なので $B^{-1/2}AB^{-1/2}$ を対称固有値問題にできる。同じ空間の基底を交換しても固有値は変わらない。

:::example[一方向の試行空間での Ritz 計算（1）]
単位区間の Dirichlet 作用素に $x^m(1-x)^2$ を試す。

$$
m=1
$$

基底関数を作り、両端の値を確認する。

$$
\begin{aligned}
\phi(x)\big|_{m=1} &= x \left(x - 1\right)^{2}
\end{aligned}
$$

質量行列の成分を Beta 積分で求める。

$$
\begin{aligned}
B_{11}\big|_{m=1} &= \frac{1}{105}
\end{aligned}
$$

剛性行列では積の微分を二乗する。

$$
\begin{aligned}
A_{11}\big|_{m=1} &= \frac{2}{15}
\end{aligned}
$$

分子を分母で割って一試行空間での上界を求める。

$$
\begin{aligned}
\lambda_{\rm Ritz}\big|_{m=1} &= 14
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
\phi(x)\big|_{x=\frac{1}{3}} &= \frac{4}{27} \\
\frac{\partial^1 (\phi(x))}{\partial x^1} &= \left(x - 1\right) \left(3 x - 1\right) \\
\left.\frac{\partial^1 (\phi(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= 0
\end{aligned}
$$
:::

:::example[一方向の試行空間での Ritz 計算（2）]
単位区間の Dirichlet 作用素に $x^m(1-x)^2$ を試す。

$$
m=2
$$

基底関数を作り、両端の値を確認する。

$$
\begin{aligned}
\phi(x)\big|_{m=2} &= x^{2} \left(x - 1\right)^{2}
\end{aligned}
$$

質量行列の成分を Beta 積分で求める。

$$
\begin{aligned}
B_{11}\big|_{m=2} &= \frac{1}{630}
\end{aligned}
$$

剛性行列では積の微分を二乗する。

$$
\begin{aligned}
A_{11}\big|_{m=2} &= \frac{2}{105}
\end{aligned}
$$

分子を分母で割って一試行空間での上界を求める。

$$
\begin{aligned}
\lambda_{\rm Ritz}\big|_{m=2} &= 12
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
\phi(x)\big|_{x=\frac{1}{3}} &= \frac{4}{81} \\
\frac{\partial^1 (\phi(x))}{\partial x^1} &= 2 x \left(x - 1\right) \left(2 x - 1\right) \\
\left.\frac{\partial^1 (\phi(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{4}{27}
\end{aligned}
$$
:::

:::example[一方向の試行空間での Ritz 計算（3）]
単位区間の Dirichlet 作用素に $x^m(1-x)^2$ を試す。

$$
m=3
$$

基底関数を作り、両端の値を確認する。

$$
\begin{aligned}
\phi(x)\big|_{m=3} &= x^{3} \left(x - 1\right)^{2}
\end{aligned}
$$

質量行列の成分を Beta 積分で求める。

$$
\begin{aligned}
B_{11}\big|_{m=3} &= \frac{1}{2310}
\end{aligned}
$$

剛性行列では積の微分を二乗する。

$$
\begin{aligned}
A_{11}\big|_{m=3} &= \frac{2}{315}
\end{aligned}
$$

分子を分母で割って一試行空間での上界を求める。

$$
\begin{aligned}
\lambda_{\rm Ritz}\big|_{m=3} &= \frac{44}{3}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
\phi(x)\big|_{x=\frac{1}{3}} &= \frac{4}{243} \\
\frac{\partial^1 (\phi(x))}{\partial x^1} &= x^{2} \left(x - 1\right) \left(5 x - 3\right) \\
\left.\frac{\partial^1 (\phi(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{8}{81}
\end{aligned}
$$
:::

:::example[一方向の試行空間での Ritz 計算（4）]
単位区間の Dirichlet 作用素に $x^m(1-x)^2$ を試す。

$$
m=4
$$

基底関数を作り、両端の値を確認する。

$$
\begin{aligned}
\phi(x)\big|_{m=4} &= x^{4} \left(x - 1\right)^{2}
\end{aligned}
$$

質量行列の成分を Beta 積分で求める。

$$
\begin{aligned}
B_{11}\big|_{m=4} &= \frac{1}{6435}
\end{aligned}
$$

剛性行列では積の微分を二乗する。

$$
\begin{aligned}
A_{11}\big|_{m=4} &= \frac{2}{693}
\end{aligned}
$$

分子を分母で割って一試行空間での上界を求める。

$$
\begin{aligned}
\lambda_{\rm Ritz}\big|_{m=4} &= \frac{130}{7}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
\phi(x)\big|_{x=\frac{1}{3}} &= \frac{4}{729} \\
\frac{\partial^1 (\phi(x))}{\partial x^1} &= 2 x^{3} \left(x - 1\right) \left(3 x - 2\right) \\
\left.\frac{\partial^1 (\phi(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{4}{81}
\end{aligned}
$$
:::

:::example[非直交な二多項式の試行空間]
一次の因子が同じ対称・反対称成分を含むので、基底の内積が零にならない場合を扱う。

各関数が端点で零であることを確かめ、二つを列の順に固定する。

$$
\begin{aligned}
\phi_{1}(x) &= x \left(1 - x\right) \\
\phi_{2}(x) &= x^{2} \left(1 - x\right)
\end{aligned}
$$

微分して、積を積分すると剛性行列を得る。

$$
\begin{aligned}
\phi_{1}^\prime &= 1 - 2 x \\
\phi_{2}^\prime &= - x^{2} + 2 x \left(1 - x\right) \\
A &= \begin{pmatrix}\frac{1}{3} & \frac{1}{6}\\\frac{1}{6} & \frac{2}{15}\end{pmatrix}
\end{aligned}
$$

関数自体の積を積分して、質量行列を作る。ここで零でない対角外成分も残して計算する。

$$
\begin{aligned}
B &= \begin{pmatrix}\frac{1}{30} & \frac{1}{60}\\\frac{1}{60} & \frac{1}{105}\end{pmatrix} \\
\det B &= \frac{1}{25200}
\end{aligned}
$$

一般化固有方程式を行列式にして根を取る。第一の根はこの試行空間での最小の Rayleigh 商である。

$$
\begin{aligned}
\det(tB-A) &= \frac{\left(t - 42\right) \left(t - 10\right)}{25200} \\
(\lambda_j)_j &= \left( 10, \  42\right)
\end{aligned}
$$
:::


:::example[対称な試行空間を増やす]
反転に対して偶な二多項式を選ぶ。特に第一固有値の上界を、一試行関数の上界と比べる。

各関数が端点で零であることを確かめ、二つを列の順に固定する。

$$
\begin{aligned}
\phi_{1}(x) &= x \left(1 - x\right) \\
\phi_{2}(x) &= x^{2} \left(1 - x\right)^{2}
\end{aligned}
$$

微分して、積を積分すると剛性行列を得る。

$$
\begin{aligned}
\phi_{1}^\prime &= 1 - 2 x \\
\phi_{2}^\prime &= x^{2} \left(2 x - 2\right) + 2 x \left(1 - x\right)^{2} \\
A &= \begin{pmatrix}\frac{1}{3} & \frac{1}{15}\\\frac{1}{15} & \frac{2}{105}\end{pmatrix}
\end{aligned}
$$

関数自体の積を積分して、質量行列を作る。ここで零でない対角外成分も残して計算する。

$$
\begin{aligned}
B &= \begin{pmatrix}\frac{1}{30} & \frac{1}{140}\\\frac{1}{140} & \frac{1}{630}\end{pmatrix} \\
\det B &= \frac{1}{529200}
\end{aligned}
$$

一般化固有方程式を行列式にして根を取る。第一の根はこの試行空間での最小の Rayleigh 商である。

$$
\begin{aligned}
\det(tB-A) &= \frac{t^{2} - 112 t + 1008}{529200} \\
(\lambda_j)_j &= \left( 56 - 4 \sqrt{133}, \  4 \sqrt{133} + 56\right)
\end{aligned}
$$
:::


:::example[偶奇に分けて同じ空間を解く]
最初の空間と同じ空間を、反転に関する偶・奇の基底で表す。係数行列の対角外成分が両方とも消える。

各関数が端点で零であることを確かめ、二つを列の順に固定する。

$$
\begin{aligned}
\phi_{1}(x) &= x \left(1 - x\right) \\
\phi_{2}(x) &= x \left(1 - x\right) \left(2 x - 1\right)
\end{aligned}
$$

微分して、積を積分すると剛性行列を得る。

$$
\begin{aligned}
\phi_{1}^\prime &= 1 - 2 x \\
\phi_{2}^\prime &= 2 x \left(1 - x\right) - x \left(2 x - 1\right) + \left(1 - x\right) \left(2 x - 1\right) \\
A &= \begin{pmatrix}\frac{1}{3} & 0\\0 & \frac{1}{5}\end{pmatrix}
\end{aligned}
$$

関数自体の積を積分して、質量行列を作る。ここで零でない対角外成分も残して計算する。

$$
\begin{aligned}
B &= \begin{pmatrix}\frac{1}{30} & 0\\0 & \frac{1}{210}\end{pmatrix} \\
\det B &= \frac{1}{6300}
\end{aligned}
$$

一般化固有方程式を行列式にして根を取る。第一の根はこの試行空間での最小の Rayleigh 商である。

$$
\begin{aligned}
\det(tB-A) &= \frac{\left(t - 42\right) \left(t - 10\right)}{6300} \\
(\lambda_j)_j &= \left( 10, \  42\right)
\end{aligned}
$$
:::

### 試行空間の誤差を残差で測る

規格化された試行関数 $u$ が作用素 $L$ の領域に属するとし、$\rho=\inner{Lu}{u}$、$r=Lu-\rho u$ と置く。離散スペクトルなら $u=\sum_j c_j\phi_j$ であり、
$$
\|r\|^2=\sum_j(\lambda_j-\rho)^2 |c_j|^2
$$
である。よって $\rho$ と全固有値との距離は $\|r\|$ 以下になる。実際、各因子がこの距離の二乗以上であり、係数の二乗和が $1$ だからである。近い固有値を識別するには、例えば固有値の間隔について別の情報を使う。自然境界条件を変分の形式で満たすだけの試行関数では、$Lu$ を作用素の領域で定義できるかを先に確認する。

単位区間で $L=-D^2$、$u=\sqrt{30}\,x(1-x)$ とすると $\rho=10$、$Lu=2\sqrt{30}$ である。規格化により $\|r\|^2=\|Lu\|^2-\rho^2=120-100=20$ となる。区間 $[10-\sqrt{20},10+\sqrt{20}]$ には第一固有値 $\pi^2$ が入り、第二固有値 $4\pi^2$ は入らない。最小最大原理が与える $10$ は片側の上界であり、残差からはスペクトルまでの距離も評価できる。

## 演習問題

::::exercise{level=basic}
$J[u]=\int_\Omega(|\grad u|^2+2fu)$ の停留方程式を求めよ。
:::solution
$F_u=2f,F_{\grad u}=2\grad u$ なので $\delta J$ の内部係数は $2f-2\Delta u$ である。よって $\Delta u=f$ となる。
:::
::::

::::exercise{level=basic}
$F(y'')=(y'')^2/2$ の内部停留条件を求めよ。
:::solution
$F_2=y''$ であり、二回部分積分すると $F_2''=y''''$ が残る。従って $y''''=0$ である。
:::
::::

::::exercise{level=basic}
$y(0)=0$、$y(1)$ が自由な $J=\int_0^1 (y')^2/2$ の停留解を求めよ。
:::solution
内部は $-y''=0$ で、自由端では $y'(1)=0$ である。線形解の傾きが零、固定端で定数も零なので $y=0$ である。
:::
::::

::::exercise{level=basic}
$\int_0^1 y=1$ を課した本文の等周問題の $\lambda$ と $y(1/2)$ を求めよ。
:::solution
$\lambda=12$ で $y=6x(1-x)$ となる。$x=1/2$ では $y=3/2$ である。
:::
::::

::::exercise{level=advanced}
自由端を持つ $J=\int [(y'')^2/2+q(x)y]\dd x$ の端条件を求めよ。
:::solution
境界項は $[y''\eta'-y'''\eta]$ である。端点付近の一次式にカットオフを掛ければ $\eta,\eta'$ の端値を独立に指定できるため、$y''=y'''=0$ が必要である。
:::
::::

::::exercise{level=advanced}
$J=\int_0^1[(y')^2-16y^2]/2$ の零解が最小でないことを示せ。
:::solution
$y=\epsilon\sin\pi x$ は端点を保つ。値は $\epsilon^2(\pi^2-16)/4<0$ で、$\epsilon$ を任意に小さくできるので零解は局所最小でもない。
:::
::::

::::exercise{level=advanced}
$v=x(1-x)(2x-1)$ の単位区間 Dirichlet の Ritz 値を計算せよ。
:::solution
$v=-2x^3+3x^2-x$ である。二乗を展開して積分すると $\int_0^1v^2=1/210$、導関数の二乗は $\int_0^1(v')^2=1/5$ である。比は $42$ となる。
:::
::::

::::exercise{level=advanced}
グラフの長さの端点が $y=2-x$ 上を動く場合に $(1,1)$ が最短であることを、横断条件以外でも確かめよ。
:::solution
点 $(x,2-x)$ まで任意の曲線の長さは直線距離以上である。距離の二乗は $x^2+(2-x)^2=2(x-1)^2+2$ だから $x=1$ で唯一最小になる。その点への直線の傾きは $1$、境界の傾きは $-1$ で横断条件とも一致する。
:::
::::

## 確認テスト

::::quiz{answer="10" tol="1e-9"}
$v=x(1-x)$ の単位区間 Dirichlet の Rayleigh 商を答えよ。

:::explanation
分子は $1/3$、分母は $1/30$ なので $(1/3)/(1/30)=10$ である。
:::
::::

::::quiz
端点が自由で $F=(y')^2/2$ の場合、自然境界条件はどれか。

- [x] $y'=0$
- [ ] $y=0$
- [ ] $y''=1$

:::explanation
端点に $y'\eta$ が残る。$\eta$ が任意なので係数の $y'$ が零である。
:::
::::

::::quiz{answer="12" tol="1e-9"}
$\int_0^1y=1$ の本文の等周問題で $\lambda$ を答えよ。

:::explanation
$y=(\lambda/2)x(1-x)$ の積分は $\lambda/12$ だから $\lambda=12$ が必要である。
:::
::::

::::quiz
試行空間を包含の関係で広げると第一 Ritz 値はどうなるか。

- [x] 増えない
- [ ] 必ず増える
- [ ] 行列の対称性を失う

:::explanation
第一値は商の最小値であり、以前の試行関数も新しい空間に含まれるので以前より大きくならない。
:::
::::

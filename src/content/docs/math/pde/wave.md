---
title: 波動方程式
description: 一次元から三次元までの解の表示を導き、有限伝播と有界領域での振動を調べる。
notation:
  - symbol: 'c'
    meaning: 波の伝播速度
  - symbol: 'M_r'
    meaning: 半径 r の球面平均
  - symbol: 'S(t)'
    meaning: 初速度から時刻 t の変位への解作用素
---

## この章の目標

- 特性変数と反射から一次元の波を解く。
- 球面平均と次元の降下から三次元と二次元の公式を導く。
- Huygens の原理と有限伝播速度を区別する。
- 外力のある波を Duhamel の公式と固有関数で表す。

## d'Alembert の公式と反射法

一定速度 $c>0$ とし、$u_{tt}-c^2\Delta u=f$ を :term[波動方程式]{en="wave equation" yomi="はどうほうていしき"} という。$g=u(0,\cdot)$、$h=u_t(0,\cdot)$ をそれぞれ初期変位と初速度という。一定張力 $\tau$ と線密度 $\mu$ の弦では、小さな傾きに対する横方向の力が長さ当たり $\tau u_{xx}$ となる。従って $\mu u_{tt}=\tau u_{xx}$、$c=\sqrt{\tau/\mu}$ である。例えば $\tau=100\,\mathrm N$、$\mu=10^{-2}\,\mathrm{kg\,m^{-1}}$ なら $c=100\,\mathrm{m\,s^{-1}}$ である。

:::theorem[d'Alembert の公式]{#thm-pde-dalembert}
全直線で $g\in C^2$、$h\in C^1$ とする。斉次問題の一意な古典解は

$$
u(t,x)=\tfrac12(g(x-ct)+g(x+ct))+\frac1{2c}\int_{x-ct}^{x+ct}h(y)\,\dd y .
\label{eq-pde-dalembert}
$$

:::

:::proof
$\xi=x-ct$、$\eta=x+ct$ とおくと作用素は $-4c^2\partial_{\xi\eta}$ であり、解は $A(\xi)+B(\eta)$ である。初期条件は $A+B=g$、$-cA'+cB'=h$ となる。$A'=(g'-h/c)/2$、$B'=(g'+h/c)/2$ を積分すると表示を得る。積分定数は和だけで現れ、$g$ により決まる。
:::

半直線の斉次 Dirichlet 条件では $g,h$ を奇延長し、斉次 Neumann 条件では偶延長する。式は延長の対称性を保つので、境界で変位またはその空間微分が消える。古典解のためには、奇延長では $g(0)=g''(0)=h(0)=0$、偶延長では $g'(0)=h'(0)=0$ を要求する。有限区間では両端で反射を繰り返す。同じ種類の両端条件なら $2L$ 周期、異なる種類なら $2L$ の平行移動で符号が反転して $4L$ 周期となる。

:::example[二つの向きに進む成分]
$c=2$、$g=\cos x$、$h=0$ なら $u=\cos x\cos2t$ である。$g=0$、$h=\cos x$ なら積分により $u=\cos x\sin2t/2$ となる。従って同じ空間モードでも初期速度によって位相が変わる。
:::

:::example[反射で変わる境界値]
$c=1$、$g=0$、$h$ を $(1,2)$ に台を持つ非負関数とする。$x=0$ では奇延長なら積分が打ち消され $u(t,0)=0$ となる。偶延長なら $u(t,0)=\int_0^t h$ となり、$t>2$ では $h$ の積分に等しい。後者で消えるのは $u_x(t,0)$ である。
:::

### 特性座標で初期線を横切る二方向

$\xi=x-ct$、$\eta=x+ct$ と置くと $\partial_t=-c\partial_\xi+c\partial_\eta$、$\partial_x=\partial_\xi+\partial_\eta$ である。従って $u_{tt}-c^2u_{xx}=-4c^2u_{\xi\eta}$ であり、混合微分が零なら $u=F(\xi)+G(\eta)$ と書ける。初期線 $\xi=\eta=x$ で $F+G=g$、$-cF'+cG'=h$ なので
$$
F'=(g'-h/c)/2,\qquad G'=(g'+h/c)/2 .
$$
これを積分すると d'Alembert の公式であり、$F,G$ の積分定数は和だけが決まる。初期値の一つが欠けると、二方向への分解は決まらない。特に $h=-cg'$ を指定すると $G'=0$ なので純粋な右向き波になる。

:::example[両端が異なる型の反射]
区間 $(0,L)$ で $u(t,0)=0$、$u_x(t,L)=0$ とする。初期データを左端の周りで奇、右端の周りで偶に延長する。二つの鏡映の合成は $2L$ の平行移動であり、符号が一回反転するため $g(x+2L)=-g(x)$、$g(x+4L)=g(x)$ となる。$h$ も同じである。全直線の解は d'Alembert の式から両方の対称性を保存する。$L=\pi$ で $g(x)=\sin(x/2)$、$h=0$ とすれば $u=\cos(ct/2)\sin(x/2)$ である。波形が符号を反転する移動距離は $2\pi$ であり、往復一回で奇延長の反転が一度現れる。
:::

:::example[左右の波を足して四次の初期データを運ぶ]
$c=2$、$g(x)=x^4$、$h(x)=3x^2$、$P=4\partial_{xx}$ とする。二方向の表示を Taylor 係数からも得る。

時間二階の方程式なので偶数階と奇数階を分ける。$g$ からは $g,Pg,P^2g,\ldots$、$h$ からは $h,Ph,P^2h,\ldots$ を得る。空間の多項式に定数係数の二階作用素を掛けると次数が少なくとも二つ下がる。従ってこの場合には両方の列が有限回で零になり、級数の収束を別に仮定する必要はない。

$$
\begin{aligned}
P^{0}g&=x^{4}\\
P^{1}g&=48 x^{2}\\
P^{2}g&=384\\
P^{0}h&=3 x^{2}\\
P^{1}h&=24\\
P^{3}g&=0
\end{aligned}
$$

$P^j g$ を $t^{2j}/(2j)!$、$P^j h$ を $t^{2j+1}/(2j+1)!$ と掛けて足す。時間で二回微分すると各項の添字が一つ下がり、空間作用素を一回掛けた列と一致する。時刻零の変位と速度は別々に確認する。

$$
\begin{aligned}
u(t,x)&=16 t^{4} + 4 t^{3} + 24 t^{2} x^{2} + 3 t x^{2} + x^{4}\\
u_{tt}&=\frac{\partial^{2}}{\partial t^{2}} \left(16 t^{4} + 4 t^{3} + 24 t^{2} x^{2} + 3 t x^{2} + x^{4}\right)=24 \left(8 t^{2} + t + 2 x^{2}\right)\\
Pu&=24 \left(8 t^{2} + t + 2 x^{2}\right)\\
u(0,x)&=x^{4},\qquad t=0\\
u_t(0,x)&=3 x^{2},\qquad t=0
\end{aligned}
$$

空間の導関数を測ると、時間の多項式の各係数が局所でどう効くかが分かる。初期関数との差は単に平行移動で置き換えず、以下では導関数と時刻を代入して計算する。

$$
\begin{aligned}
\partial_{x} u&=\frac{\partial}{\partial x} \left(16 t^{4} + 4 t^{3} + 24 t^{2} x^{2} + 3 t x^{2} + x^{4}\right)=2 x \left(24 t^{2} + 3 t + 2 x^{2}\right)\\
g&=1,\qquad x=1\\
u&=10,\qquad t=\frac{1}{2},x=1\\
u_t&=38,\qquad t=\frac{1}{2},x=1\\
h&=3,\qquad x=1
\end{aligned}
$$

d’Alembert の変位は $((x+2t)^4+(x-2t)^4)/2$ であり、奇数の交差項が相殺する。速度の積分は $((x+2t)^3-(x-2t)^3)/4$ で、$t$ の奇数べきだけが残る。両者を足すと上の多項式に一致する。例えば $x=1,t=1/2$ では足し合わせる端点は $0,2$ である。
:::

### 二つの端点を観測することと区間を平均すること

$d$ を測定点から初期速度の支持までの最短距離、$D$ を最長距離とする。$g=0$、$h\in L^1(\R)$ がコンパクト台を持つなら
$$
u(t,x)=\frac{1}{2c}\int_{x-ct}^{x+ct}h(y)\dd y
$$
であり、$ct<d$ では零、$ct>D$ では一定 $\int h/(2c)$ になる。速度の総量が零なら通過後の変位も零になる。これに対し初期変位の二項は端点を読むので、両端が支持を通り過ぎれば総量にかかわらず零になる。有限伝播は信号が届く前の値を決めるが、信号が通過した後の値を零とするにはデータの種類も必要である。

:::example[符号の異なる二つの初速度]
$c=2$、$h=\mathbf1_{[1,2]}-\mathbf1_{[2,3]}$ とし、$x=0$ で観測する。$0<t<1/2$ では積分範囲はまだ台に届かない。$1/2<t<1$ では正の区間のうち長さ $2t-1$ が入るので $u=(2t-1)/4$ である。$1<t<3/2$ では正の区間全体と負の区間の長さ $2t-2$ が入るので $u=(3-2t)/4$ である。$t>3/2$ では全積分が零なので $u=0$。$t=3/4$ と $5/4$ では同じ $1/8$ だが時間微分の符号は反対である。段差を平滑化して式から極限を取れば、この応答は連続の弱解として正当化できる。
:::

## 球面平均と Kirchhoff の公式

滑らかな $g$ の球面平均を $M_rg(x)=\sigma_{n-1}^{-1}\int_{\mathbb S^{n-1}}g(x+r\theta)\,\dd S_\theta$ とする。:ref[thm-anl-divergence] により $\sigma_{n-1}r^{n-1}\partial_rM_rg=\int_{B_r(x)}\Delta g$ である。両辺を半径で微分して面積で割ると

$$
\partial_{rr}M_rg+\frac{n-1}{r}\partial_rM_rg=\Delta_xM_rg.
\label{eq-pde-euler-poisson-darboux}
$$

右辺では積分の下で $x$ を微分した。これを :term[Euler–Poisson–Darboux 方程式]{en="Euler–Poisson–Darboux equation" yomi="おいらーぽあそんだるぶーほうていしき"} と呼ぶ。

:::theorem[Kirchhoff の公式]{#thm-pde-kirchhoff}
三次元全空間の滑らかなデータには

$$
u(t,x)=\partial_t(tM_{ct}g(x))+tM_{ct}h(x)
\label{eq-pde-kirchhoff}
$$

が斉次問題の古典解を与える。
:::

:::proof
$v(t,x)=tM_{ct}g(x)$ の二階時間微分は $2c\partial_rM_rg+tc^2\partial_{rr}M_rg$ である。$r=ct$、$n=3$ の :eqref[eq-pde-euler-poisson-darboux] によりこれは $c^2\Delta_xv$ である。球面平均は半径について偶なので $v(0)=v_{tt}(0)=0$、$v_t(0)=g$ である。よって $v_t$ は初期値 $(g,0)$ の解であり、$g$ を $h$ に替えた $v$ は初期値 $(0,h)$ の解である。後の円錐エネルギー評価で一意性を示す。
:::

:::example[面積を使う表示]
$tM_{ct}h=(4\pi c^2t)^{-1}\int_{|y-x|=ct}h(y)\,\dd S_y$ である。$c=1$、$h=e^{-|y|^2}$、$g=0$ のとき、中心での値は $te^{-t^2}$ である。球面積 $4\pi t^2$ と積分の前の $1/(4\pi t)$ が合わさって $t$ になる。
:::

:::example[多項式の平均]
$g=|x|^2$、$h=0$ なら $M_{ct}g=|x|^2+c^2t^2$ である。従って $u=|x|^2+3c^2t^2$ となり、$u_{tt}=6c^2=c^2\Delta u$ である。
:::

### 球面上の規格化と分布核

三次元の初速度作用素を $S(t)h=tM_{ct}h$ とする。測度 $E_t=(4\pi c^2t)^{-1}\dd S_{|y|=ct}$ と書けば $E_t*h=S(t)h$ である。総質量は $t$ なので $E_t\to0$、$E_t/t\rightharpoonup\delta_0$ となる。$t<0$ で零に延長すると $(\partial_t^2-c^2\Delta)E=\delta_{(0,0)}$ である。実際、$t>\eps$ で試験関数と部分積分すると下端に $\langle E'_\eps,\phi(\eps)\rangle-\langle E_\eps,\phi_t(\eps)\rangle$ が残る。$E'_\eps=\partial_\eps(\eps M_{c\eps})$ は試験関数に対して恒等作用素に収束し、第二項は零に収束する。この因果的な基本解は外力の時空間畳み込みにそのまま使える。

:::example[四次多項式の球面平均]
$h(y)=y_1^4$ とする。単位球面で $z=\theta_1$ の面素は角度を積分すると $2\pi\dd z$ なので、$\theta_1^2$ と $\theta_1^4$ の平均はそれぞれ $1/3,1/5$ である。奇数の項は消え、
$$
S(t)h(x)=t x_1^4+2c^2t^3x_1^2+\frac{c^4t^5}{5}
$$
となる。時間二階微分は $12c^2 t x_1^2+4c^4t^3$ であり、空間 Laplacian の $c^2$ 倍に一致する。中心で $t^5$ の項が現れるのは、球面上で四次の初速度を測るためである。
:::

:::example[球面の偶数モーメントで混合項を求める]
三次元で $c=1$、$g(x,y,z)=x^2y^2$、$h=x^2+2z^2$、$P=\Delta$ とする。半径方向の平均では座標の混合が現れる。

時間二階の方程式なので偶数階と奇数階を分ける。$g$ からは $g,Pg,P^2g,\ldots$、$h$ からは $h,Ph,P^2h,\ldots$ を得る。空間の多項式に定数係数の二階作用素を掛けると次数が少なくとも二つ下がる。従ってこの場合には両方の列が有限回で零になり、級数の収束を別に仮定する必要はない。

$$
\begin{aligned}
P^{0}g&=x^{2} y^{2}\\
P^{1}g&=2 \left(x^{2} + y^{2}\right)\\
P^{2}g&=8\\
P^{0}h&=x^{2} + 2 z^{2}\\
P^{1}h&=6\\
P^{3}g&=0
\end{aligned}
$$

$P^j g$ を $t^{2j}/(2j)!$、$P^j h$ を $t^{2j+1}/(2j+1)!$ と掛けて足す。時間で二回微分すると各項の添字が一つ下がり、空間作用素を一回掛けた列と一致する。時刻零の変位と速度は別々に確認する。

$$
\begin{aligned}
u(t,x)&=\frac{t^{4}}{3} + t^{3} + t^{2} \left(x^{2} + y^{2}\right) + t \left(x^{2} + 2 z^{2}\right) + x^{2} y^{2}\\
u_{tt}&=\frac{\partial^{2}}{\partial t^{2}} \left(\frac{t^{4}}{3} + t^{3} + t^{2} \left(x^{2} + y^{2}\right) + t \left(x^{2} + 2 z^{2}\right) + x^{2} y^{2}\right)=2 \left(2 t^{2} + 3 t + x^{2} + y^{2}\right)\\
Pu&=2 \left(2 t^{2} + 3 t + x^{2} + y^{2}\right)\\
u(0,x)&=x^{2} y^{2},\qquad t=0\\
u_t(0,x)&=x^{2} + 2 z^{2},\qquad t=0
\end{aligned}
$$

空間の導関数を測ると、時間の多項式の各係数が局所でどう効くかが分かる。初期関数との差は単に平行移動で置き換えず、以下では導関数と時刻を代入して計算する。

$$
\begin{aligned}
\partial_{x} u&=\frac{\partial}{\partial x} \left(\frac{t^{4}}{3} + t^{3} + t^{2} \left(x^{2} + y^{2}\right) + t \left(x^{2} + 2 z^{2}\right) + x^{2} y^{2}\right)=2 x \left(t^{2} + t + y^{2}\right)\\
\partial_{y} u&=\frac{\partial}{\partial y} \left(\frac{t^{4}}{3} + t^{3} + t^{2} \left(x^{2} + y^{2}\right) + t \left(x^{2} + 2 z^{2}\right) + x^{2} y^{2}\right)=2 y \left(t^{2} + x^{2}\right)\\
\partial_{z} u&=\frac{\partial}{\partial z} \left(\frac{t^{4}}{3} + t^{3} + t^{2} \left(x^{2} + y^{2}\right) + t \left(x^{2} + 2 z^{2}\right) + x^{2} y^{2}\right)=4 t z\\
g&=4,\qquad x=1,y=2,z=0\\
u&=\frac{118}{3},\qquad t=2,x=1,y=2,z=0\\
u_t&=\frac{131}{3},\qquad t=2,x=1,y=2,z=0\\
h&=1,\qquad x=1,y=2,z=0
\end{aligned}
$$

球面上で $\theta_i^2$ の平均は $1/3$、$\theta_i^2\theta_j^2$（$i\ne j$）は $1/15$ である。後者は $(\sum\theta_i^2)^2=1$ を平均し、各 $\theta_i^4$ の平均が $1/5$ であることから得られる。従って Kirchhoff の公式を使うとき、初期変位の $t^4$ の係数は $(1/15)$ を $\partial_t(t\cdot)$ で処理して $1/3$ になる。ここで独立な一次元の平均の積を用いると $1/9$ になり、球面での座標の従属性を見落とす。
:::

### Euler–Poisson–Darboux の式の原点条件

球面平均 $M_rg$ は半径 $r$ の偶関数に延長でき、$M_0g=g$、$\partial_rM_rg|_{r=0}=0$ である。後者は各方向 $\theta$ の寄与を $-\theta$ と対にすると従う。球面平均を計算して $r$ で微分するだけでなく、$x$ についての微分も積分の中に入る。従って $M_r\Delta g=\Delta_x M_rg$ であり、球面平均は
$$
\partial_{rr} M_rg+\frac{n-1}{r}\partial_r M_rg=\Delta_x M_rg
$$
という方程式を満たす。$n=3$ では $W=rM_rg$ に対し $W_{rr}=r(M_{rr}+2M_r/r)=\Delta_x W$。$r=ct$ と置けば波動方程式になり、$W(0)=0$、$W_r(0)=g$ である。従って初速度の解は $W(ct)/c$ である。初期変位の解を得るにはさらに時間微分し、$W_{rr}(0)=0$ を使う。ここで半径因子を掛けずに $M_{ct}$ だけを用いると、一次の項 $2M_r/r$ が残って波動方程式にはならない。

:::example[球面平均で三角関数をそのまま積分する]
$g(y)=\cos(k y_1)$、$n=3$ とする。方位角を積分すると $\theta_1$ は $[-1,1]$ に密度 $1/2$ を持つため
$$
M_{ct}g(x)=\frac12\int_{-1}^1 \cos(k(x_1+ctz))\dd z
=\cos(kx_1)\frac{\sin(kct)}{kct}.
$$
$k=0$ は極限で読む。従って初速度の作用素は $\sin(kct)/(kc)$、変位の作用素は $\cos(kct)$ を掛ける。例えば $k=2,c=3$ で $t=\pi/12$ なら初速度からの係数は $1/6$、同じ空間形の初期変位からの係数は零になる。これは Fourier 変換で波の周波数 $kc=6$ を得る計算と一致する。
:::

## Hadamard の降下法と Poisson の公式

:::theorem[二次元の Poisson の公式]{#thm-pde-wave-poisson}
二次元全空間で初期変位 $g$ と初速度 $h$ を滑らかとする。初期変位 $0$ の場合の解を $S(t)h$ と書くと、一般の斉次問題の解は

$$
S(t)h(x)=\frac1{2\pi c}\int_{|y-x|<ct}\frac{h(y)}{\sqrt{c^2t^2-|y-x|^2}}\,\dd y,\qquad
u=\partial_tS(t)g+S(t)h .
\label{eq-pde-wave-poisson}
$$

:::

:::proof
三次元のデータを第 3 座標から独立に取ると、公式から解も第 3 座標に依存しない。従って三次元の式は二次元の PDE を満たす。この手続きを :term[降下法]{en="method of descent" yomi="こうかほう"} という。半径 $R=ct$ の球面の上下半分は高さ $\pm\sqrt{R^2-|y-x|^2}$ のグラフであり、各面素は $R/\sqrt{R^2-|y-x|^2}\,\dd y$ である。二つの面の寄与に Kirchhoff の面積係数 $1/(4\pi c^2t)$ を掛けると表示を得る。$y=x+ctz$ と置き換えれば時間微分では積分領域が固定され、重み $(1-|z|^2)^{-1/2}$ は可積分である。
:::

:::example[初速度の規格化]
$h=2$ なら $S(t)h=(2/c)\int_0^{ct}r/\sqrt{c^2t^2-r^2}\,\dd r=2t$ である。$h(y)=|y|^2$ なら同じ中心積分で $S(t)h(0)=2c^2t^3/3$ になる。後者では $u=t|x|^2+2c^2t^3/3$ と書け、$u_{tt}=4c^2t=c^2\Delta u$ と確認できる。
:::

:::example[角度の平均で消える初速度]
二次元で $h(y)=y_1$、$g=0$ とする。中心を $x$ とした Poisson 積分で $y=x+z$ とおくと、$z_1$ の項だけが反射に関して奇になる。残る定数 $x_1$ には先の規格化を適用できるので $u=tx_1$ である。初期データの奇偶によっては、正の核の領域内でも正負の値が打ち消し合う。
:::

### 円盤核の端点の特異性を扱う

二次元の公式では被積分核が円盤の縁で発散する。しかし $\int_0^{ct}r((ct)^2-r^2)^{-1/2}\dd r=ct$ なので可積分である。初期変位に対して時間微分するときは、先に $y=x+ctz$ と変数変換し、
$$
S_2(t)g(x)=\frac{t}{2\pi}\int_{|z|<1}\frac{g(x+ctz)}{\sqrt{1-|z|^2}}\dd z
$$
と書く。積分領域も特異因子も時間に依存しなくなるので、微分は $t$ と $g$ の引数だけに作用する。$g$ が $C^1$ ならコンパクトな範囲で導関数が有界であり、優収束により微分を積分の中へ入れられる。元の式のまま積分端の値を代入すると無限大が現れ、誤った式になる。

:::example[円盤の中心における段差の応答]
$c=1$、$g=0$、$h=\mathbf1_{B_1}$ とする。中心の値は $u(t,0)=\int_0^{\min(t,1)}r(t^2-r^2)^{-1/2}\dd r$ なので、$t\le1$ では $t$、$t>1$ では $t-\sqrt{t^2-1}$ である。例えば $t=5/4$ では $1/2$ である。$t\to\infty$ で有理化すると $u=1/(t+\sqrt{t^2-1})\sim(2t)^{-1}$ となる。三次元で同じ半径の台を選ぶと $t>1$ に中心は直ちに零になり、両者の差は核が支持される位置の違いによる。
:::

:::example[降下した核で半径四次の寄与を測る]
二次元で $c=1$、$g=0$、$h=(x^2+y^2)^2$、$P=\Delta$ とする。円盤全体からの初速度の寄与を中心でも求める。

時間二階の方程式なので偶数階と奇数階を分ける。$g$ からは $g,Pg,P^2g,\ldots$、$h$ からは $h,Ph,P^2h,\ldots$ を得る。空間の多項式に定数係数の二階作用素を掛けると次数が少なくとも二つ下がる。従ってこの場合には両方の列が有限回で零になり、級数の収束を別に仮定する必要はない。

$$
\begin{aligned}
P^{0}h&=\left(x^{2} + y^{2}\right)^{2}\\
P^{1}h&=16 \left(x^{2} + y^{2}\right)\\
P^{2}h&=64\\
P^{3}g&=0
\end{aligned}
$$

$P^j g$ を $t^{2j}/(2j)!$、$P^j h$ を $t^{2j+1}/(2j+1)!$ と掛けて足す。時間で二回微分すると各項の添字が一つ下がり、空間作用素を一回掛けた列と一致する。時刻零の変位と速度は別々に確認する。

$$
\begin{aligned}
u(t,x)&=\frac{8 t^{5}}{15} + \frac{8 t^{3} \left(x^{2} + y^{2}\right)}{3} + t \left(x^{2} + y^{2}\right)^{2}\\
u_{tt}&=\frac{\partial^{2}}{\partial t^{2}} \left(\frac{8 t^{5}}{15} + \frac{8 t^{3} \left(x^{2} + y^{2}\right)}{3} + t \left(x^{2} + y^{2}\right)^{2}\right)=\frac{16 t \left(2 t^{2} + 3 x^{2} + 3 y^{2}\right)}{3}\\
Pu&=\frac{16 t \left(2 t^{2} + 3 x^{2} + 3 y^{2}\right)}{3}\\
u(0,x)&=0,\qquad t=0\\
u_t(0,x)&=\left(x^{2} + y^{2}\right)^{2},\qquad t=0
\end{aligned}
$$

空間の導関数を測ると、時間の多項式の各係数が局所でどう効くかが分かる。初期関数との差は単に平行移動で置き換えず、以下では導関数と時刻を代入して計算する。

$$
\begin{aligned}
\partial_{x} u&=\frac{\partial}{\partial x} \left(\frac{8 t^{5}}{15} + \frac{8 t^{3} \left(x^{2} + y^{2}\right)}{3} + t \left(x^{2} + y^{2}\right)^{2}\right)=\frac{4 t x \left(4 t^{2} + 3 x^{2} + 3 y^{2}\right)}{3}\\
\partial_{y} u&=\frac{\partial}{\partial y} \left(\frac{8 t^{5}}{15} + \frac{8 t^{3} \left(x^{2} + y^{2}\right)}{3} + t \left(x^{2} + y^{2}\right)^{2}\right)=\frac{4 t y \left(4 t^{2} + 3 x^{2} + 3 y^{2}\right)}{3}\\
g&=0,\qquad x=0,y=0\\
u&=\frac{8}{15},\qquad t=1,x=0,y=0\\
u_t&=\frac{8}{3},\qquad t=1,x=0,y=0\\
h&=0,\qquad x=0,y=0
\end{aligned}
$$

中心で Poisson の式を極座標にすると $\int_0^t r^5(t^2-r^2)^{-1/2}\dd r=t^5\int_0^{\pi/2}\sin^5\theta\dd\theta=8t^5/15$ である。積分は $r=t\sin\theta$ と置いた後、$\sin^5\theta=(1-\cos^2\theta)^2\sin\theta$ として求められる。$r$ を端の値 $t$ に固定しただけの球面平均とは係数が異なる。
:::

### 一次元への降下でも規格化を確認する

二次元の Poisson 表示で初速度が $y_1$ だけに依存するとする。$x=0$、固定した $y_1$ に対し $A^2=(ct)^2-y_1^2>0$ であり、もう一方の変数を積分すると
$$
\int_{-A}^A\frac{\dd y_2}{\sqrt{A^2-y_2^2}}=\pi .
$$
これは $y_2=A\sin\theta$ の置換から従う。従って二次元の初速度作用素は $1/(2\pi c)$ に $\pi$ が掛かって、一次元の $1/(2c)\int_{-ct}^{ct} h(y_1)\dd y_1$ に一致する。次元を落とすたびに積分核の形が変わるが、定数初速度 $h=1$ への作用はどの次元でも $t$ である。

## Huygens の原理

:term[Huygens の原理]{en="Huygens principle" yomi="ほいへんすのげんり"} とは、自由波の初期データが後ろ向き特性円錐と初期面の交わりの近傍で消えていれば、その観測点で波が零になる性質である。三次元では Kirchhoff の公式がこの性質を与える。初期変位の項を展開すると $M_{ct}g+ct\,\sigma_2^{-1}\int_{\mathbb S^2}\theta\cdot Dg(x+ct\theta)\,\dd S_\theta$ であり、球面上の値と導関数だけが現れる。一方、二次元では核が円盤内でも正である。有限伝播という性質だけでは円盤内の過去の値が残るかどうかは決まらない。

:::example[同じ台を持つデータを次元ごとに与える]
$c=1$、$g=0$、$h\ge0$ を単位球の内部に台を持つ積分が正の関数とする。三次元の中心では $t>1$ に球面が台と交わらず $u(t,0)=0$ になる。二次元では $t>1$ に $u(t,0)>0$ である。さらに円盤の台上で $(t^2-|y|^2)^{-1/2}=t^{-1}+O(t^{-3})$ なので、$t\to\infty$ の値は $(2\pi t)^{-1}\int h+O(t^{-3})$ である。
:::

### 一般の次元での違い

奇数次元 $n=2m+1\ge3$ では、三次元の公式を球面平均の微分によって拡張できる。$D=t^{-1}\partial_t$、$P_q=\partial_{tt}-q t^{-1}\partial_t$ と置くと $P_{q-2}D=DP_q$ である。実際、いずれも $v'''/t-qv''/t^2+qv'/t^3$ を与える。式 :eqref[eq-pde-euler-poisson-darboux] から $v=t^{n-2}M_{ct}h$ は $P_{n-3}v=c^2\Delta v$ を満たす。従って $D$ を $m-1$ 回作用させれば $P_0=\partial_{tt}$ を含む波動方程式になる。$M_{ct}h=h+O(t^2)$ は滑らかな偶関数なので、これらの微分の後に残る一次の項は $(2m-1)!!\,t h$ である。ここで $(2m-1)!!=(2m-1)(2m-3)\dotsm1$ である。よって初速度の作用素は

$$
S_n(t)h=\frac{1}{(2m-1)!!}D^{m-1}\big(t^{2m-1}M_{ct}h\big).
\label{eq-pde-wave-odd-dimension}
$$

初期変位には $\partial_t S_n(t)$ が対応する。いずれの微分も球面上の初期値と導関数のみを使うので、奇数次元 $n\ge3$ には Huygens の原理が成り立つ。例えば五次元では $S_5(t)h=tM_{ct}h+(t^2/3)\partial_t M_{ct}h$ である。$h=|x|^2$ を入れると $u=t|x|^2+5c^2t^3/3$ となり、$u_{tt}=10c^2t=c^2\Delta u$ と確かめられる。

偶数次元 $n=2m$ では、$n+1$ 次元から降下させる。半径 $ct$ の球面の上下を $n$ 次元の球へ射影し、各面素 $ct/\sqrt{c^2t^2-|x-y|^2}$ を球面積 $\sigma_{2m}(ct)^{2m}$ で割る。これを :eqref[eq-pde-wave-odd-dimension] に入れると
$$
S_{2m}(t)h(x)=
\frac{2}{\sigma_{2m}c^{2m-1}(2m-1)!!}\,
D^{m-1}\int_{|y-x|<ct}\frac{h(y)}{\sqrt{c^2t^2-|x-y|^2}}\,\dd y
$$
となる。時間微分は $y=x+ctz$ と置換してから取れる。$m=1$ では $\sigma_2=4\pi$ なので、二次元の Poisson の公式に戻る。

$h$ の台が半径 $R$ の球に収まり、中心で $t>2R/c$ とする。被積分核は $t^{-1}$ を掛けた $t^{-2}$ の収束べき級数であり、最初の係数は $1/c$ である。固定した台の上でこの級数を項別に微分できる。$D^{m-1}t^{-1}=(-1)^{m-1}(2m-3)!!t^{-(2m-1)}$ であり、$m=1$ では $(-1)!!=1$ と読む。従って $\int h\ne0$ なら中心での主要項は非零であり、十分遅い時刻にも波が残る。一次元でも、正の初速度を持つ台付きデータは :eqref[eq-pde-dalembert] の区間積分により通過後に定数の尾を残す。

:::example[固定点を通過した後の一次元の値]
$c=1$、$g=0$、$h=\mathbf1_{[-1,1]}$ とする。積分公式はこの有界データでも連続な弱解を与える。実際、段差を平滑化した古典解の表示は局所一様に収束し、分布の式を保つ。$x=0$ では $u(t,0)=\min(t,1)$ である。$t>1$ で両端 $-t,t$ の初速度は $0$ でも、二点の間の積分は非零である。
:::

### 球面値だけを比較するときの注意

三次元で半径 $R=ct$ の球面上の初期変位の値が零でも、法線微分は零とは限らない。例えば $g(y)=|y|^2-R^2$、$h=0$ とする。中心の変位作用素は $\partial_t[t((ct)^2-R^2)]$ なので、$t=R/c$ で $2R^2$ である。一方、球面の近傍で $g$ が恒等的に零なら法線微分も零になる。従って Huygens の主張は基本解の分布としての台についてのものであり、変位作用素に含まれる面測度の微分も同じ球面に台を持つ。

:::example[三次元の観測窓]
$c=2$、源の台が中心から距離 $3$ の点を中心とする半径 $1/2$ の球内にあるとする。観測点を原点とすれば、初速度だけの自由波が観測され得るのは $ct\in[5/2,7/2]$、すなわち $t\in[5/4,7/4]$ の間である。初速度が非負で、台の中で正な開集合を持つなら、その開集合を横切る球面については値は正である。二次元の式で同じ配置なら、$ct>7/2$ 以後も核は台全体で正なので、初速度が非負なら値が残る。
:::

:::example[球面からの作用と球面値との差]
三次元で $c=1$ とし、$g=(x^2+y^2+z^2)^2-16$、$h=0$ とする。半径 $2$ の球面における変位の値は零であり、時刻 $2$ の中心を調べる。

時間二階の方程式なので偶数階と奇数階を分ける。$g$ からは $g,Pg,P^2g,\ldots$、$h$ からは $h,Ph,P^2h,\ldots$ を得る。空間の多項式に定数係数の二階作用素を掛けると次数が少なくとも二つ下がる。従ってこの場合には両方の列が有限回で零になり、級数の収束を別に仮定する必要はない。

$$
\begin{aligned}
P^{0}g&=\left(x^{2} + y^{2} + z^{2}\right)^{2} - 16\\
P^{1}g&=20 \left(x^{2} + y^{2} + z^{2}\right)\\
P^{2}g&=120\\
P^{3}g&=0
\end{aligned}
$$

$P^j g$ を $t^{2j}/(2j)!$、$P^j h$ を $t^{2j+1}/(2j+1)!$ と掛けて足す。時間で二回微分すると各項の添字が一つ下がり、空間作用素を一回掛けた列と一致する。時刻零の変位と速度は別々に確認する。

$$
\begin{aligned}
u(t,x)&=5 t^{4} + 10 t^{2} \left(x^{2} + y^{2} + z^{2}\right) + \left(x^{2} + y^{2} + z^{2}\right)^{2} - 16\\
u_{tt}&=\frac{\partial^{2}}{\partial t^{2}} \left(5 t^{4} + 10 t^{2} \left(x^{2} + y^{2} + z^{2}\right) + \left(x^{2} + y^{2} + z^{2}\right)^{2} - 16\right)=20 \left(3 t^{2} + x^{2} + y^{2} + z^{2}\right)\\
Pu&=20 \left(3 t^{2} + x^{2} + y^{2} + z^{2}\right)\\
u(0,x)&=\left(x^{2} + y^{2} + z^{2}\right)^{2} - 16,\qquad t=0\\
u_t(0,x)&=0,\qquad t=0
\end{aligned}
$$

空間の導関数を測ると、時間の多項式の各係数が局所でどう効くかが分かる。初期関数との差は単に平行移動で置き換えず、以下では導関数と時刻を代入して計算する。

$$
\begin{aligned}
\partial_{x} u&=\frac{\partial}{\partial x} \left(5 t^{4} + 10 t^{2} \left(x^{2} + y^{2} + z^{2}\right) + \left(x^{2} + y^{2} + z^{2}\right)^{2} - 16\right)=4 x \left(5 t^{2} + x^{2} + y^{2} + z^{2}\right)\\
\partial_{y} u&=\frac{\partial}{\partial y} \left(5 t^{4} + 10 t^{2} \left(x^{2} + y^{2} + z^{2}\right) + \left(x^{2} + y^{2} + z^{2}\right)^{2} - 16\right)=4 y \left(5 t^{2} + x^{2} + y^{2} + z^{2}\right)\\
\partial_{z} u&=\frac{\partial}{\partial z} \left(5 t^{4} + 10 t^{2} \left(x^{2} + y^{2} + z^{2}\right) + \left(x^{2} + y^{2} + z^{2}\right)^{2} - 16\right)=4 z \left(5 t^{2} + x^{2} + y^{2} + z^{2}\right)\\
g&=-16,\qquad x=0,y=0,z=0\\
u&=64,\qquad t=2,x=0,y=0,z=0\\
u_t&=160,\qquad t=2,x=0,y=0,z=0\\
h&=0,\qquad x=0,y=0,z=0
\end{aligned}
$$

球面上で $g=0$ でも、半径微分は $4r^3=32$ である。Kirchhoff の変位項のうち $t$ を掛けた半径微分が中心で残り、$2\cdot32=64$ になる。従って Huygens の性質から消えるのは球面近傍にデータがない場合であり、球面の値が零という事実だけでは結論できない。
:::

### 初期速度を見て尾と打消しを区別する

二次元で $h$ は原点から半径 $R$ の中にコンパクトな台を持つとする。$g=0$ で $t>2R/c$ なら展開 $(1-\eta)^{-1/2}=1+\eta/2+O(\eta^2)$ を $\eta=|y|^2/(c^2t^2)$ に一様に使える。従って中心の値は
$$
u(t,0)=\frac{1}{2\pi c^2 t}\int h
+\frac{1}{4\pi c^4 t^3}\int |y|^2 h+O(t^{-5})
$$
となる。ここでは Poisson の係数 $1/(2\pi c)$ と平方根からさらに出る $1/(ct)$ をともに用いた。$\int h=0$ なら第一項が相殺するが、二次モーメントまで零になるとは限らない。例えば $h=\mathbf1_{B_1}-\tfrac14\mathbf1_{B_2}$ は総量が零であり、$\int |y|^2h=\pi/2-(1/4)8\pi=-3\pi/2$。よって $c=1$ で中心の尾は $-3/(8t^3)$ を主要項に持つ。三次元以上の奇数次元で台を球面が通過した後に厳密に消える性質と、二次元で低次モーメントが打ち消されて速く減る現象は区別される。

## エネルギー保存と依存領域

:::theorem[円錐上のエネルギー評価]{#thm-pde-wave-cone}
斉次の滑らかな波に対し $e=(u_t^2+c^2|Du|^2)/2$ とすると $\partial_te=c^2\divg(u_tDu)$ である。$B_{c(T-t)}(x_0)$ 上の $e$ の積分は $0<t<T$ で非増加となる。
:::

:::proof
方程式に $u_t$ を掛けて積の微分を整理すると局所式を得る。球の境界は法線方向に速度 $-c$ で動くため、積分の時間微分は境界上の

$$
c^2u_t\partial_\nu u-ce
=-\tfrac c2(u_t-c\partial_\nu u)^2-\tfrac{c^3}2|D_{\rm tan}u|^2
\label{eq-pde-wave-cone-flux}
$$

の積分である。ここで $D_{\rm tan}u$ は勾配の接成分である。右辺は非正なので所望の単調性が従う。
:::

二解の初期変位と初速度が円錐の底の球の内部で等しければ、差のエネルギーが円錐内で $0$ になる。したがって差の時間微分も各内点で $0$ となり、初期の差が $0$ なので頂点でも $0$ である。この議論は有限の球だけを使うため、全空間での一意性には無限遠の増大条件が要らない。固定領域の全エネルギーの時間微分は境界流束 $\int c^2u_t\partial_\nu u$ であり、斉次 Dirichlet または斉次 Neumann 条件で消える。外力があれば内部に $\int f u_t$ が加わる。

円錐内に外力がある場合も、縮む球 $B(t)=B_{c(T-t)}(x_0)$ 上のエネルギーを $E(t)$ とすれば、境界項が非正なので $E'(t)\le \norm{f(t)}_{L^2(B(t))}\sqrt{2E(t)}$ である。$\sqrt{2E(t)+\eps}$ を微分して評価し、積分後に $\eps\downarrow0$ とすると
$$
\sqrt{2E(t)}\le\sqrt{2E(0)}+\int_0^t\norm{f(s)}_{L^2(B(s))}\,\dd s
$$
を得る。したがって初期値に加えて円錐内の外力も同じなら、二解は頂点で一致する。全空間の有限エネルギーの斉次解では、$E(t)$ を時刻 $t$ の半径 $R$ の球で取って $R\to\infty$ とすると、全エネルギーが初期値以下になる。時間を反転した方程式に同じ評価を使えば逆向きも従い、全エネルギーは保存する。

:::example[観測点までの距離]
$c=2$ の一次元の波の初期台が $[-1,1]$ なら $x=7$ で $0<t<3$ の値は $0$ である。円錐が初めて初期台に達するには距離 $6$ をさかのぼる必要があるからである。
:::

:::example[定在波のエネルギー]
$(0,\pi)$ で $c=1$、$u=\cos(2t)\sin2x$ とする。運動項は $\pi\sin^2(2t)$、空間項は $\pi\cos^2(2t)$ で、時間に依存するが和は $\pi$ である。
:::

### 移動境界での流束の平方

$e=(u_t^2+c^2|Du|^2)/2$ とする。球が内向き速さ $c$ で縮むとき、その境界の寄与は $c^2u_t\partial_\nu u-ce$ になる。勾配を法線成分と接成分 $D_\tau u$ に分ければ
$$
c^2u_t\partial_\nu u-ce
=-\frac c2(u_t-c\partial_\nu u)^2-\frac{c^3}2|D_\tau u|^2 .
$$
従って円錐での単調性は、波の速度と境界の縮み方が同じであることから出る。縮む速度を $v>c$ にするとさらに $-(v-c)e$ が加わり、$v<c$ では符号を保証できない。一次元で $u=F(x-ct)$ の左側の端点では $\nu=-1$、$u_t=c\partial_\nu u$ なので、上の境界項は零になる。右端では $\nu=1$ なので $-2c^3(F')^2$ となる。左側の特性境界と並走するエネルギーと、右端から出るエネルギーの違いがこの平方に現れる。

:::example[初期速度の局所エネルギーから波の動きを読む]
一次元で $c=1$、$g=0$、$h=x^4$ とし、$P=\partial_{xx}$ を使う。多項式は全空間で有限エネルギーではないが、各後退円錐上の積分は有限である。

時間二階の方程式なので偶数階と奇数階を分ける。$g$ からは $g,Pg,P^2g,\ldots$、$h$ からは $h,Ph,P^2h,\ldots$ を得る。空間の多項式に定数係数の二階作用素を掛けると次数が少なくとも二つ下がる。従ってこの場合には両方の列が有限回で零になり、級数の収束を別に仮定する必要はない。

$$
\begin{aligned}
P^{0}h&=x^{4}\\
P^{1}h&=12 x^{2}\\
P^{2}h&=24\\
P^{3}g&=0
\end{aligned}
$$

$P^j g$ を $t^{2j}/(2j)!$、$P^j h$ を $t^{2j+1}/(2j+1)!$ と掛けて足す。時間で二回微分すると各項の添字が一つ下がり、空間作用素を一回掛けた列と一致する。時刻零の変位と速度は別々に確認する。

$$
\begin{aligned}
u(t,x)&=\frac{t^{5}}{5} + 2 t^{3} x^{2} + t x^{4}\\
u_{tt}&=\frac{\partial^{2}}{\partial t^{2}} \left(\frac{t^{5}}{5} + 2 t^{3} x^{2} + t x^{4}\right)=4 t \left(t^{2} + 3 x^{2}\right)\\
Pu&=4 t \left(t^{2} + 3 x^{2}\right)\\
u(0,x)&=0,\qquad t=0\\
u_t(0,x)&=x^{4},\qquad t=0
\end{aligned}
$$

空間の導関数を測ると、時間の多項式の各係数が局所でどう効くかが分かる。初期関数との差は単に平行移動で置き換えず、以下では導関数と時刻を代入して計算する。

$$
\begin{aligned}
\partial_{x} u&=\frac{\partial}{\partial x} \left(\frac{t^{5}}{5} + 2 t^{3} x^{2} + t x^{4}\right)=4 t x \left(t^{2} + x^{2}\right)\\
g&=0,\qquad x=1\\
u&=\frac{121}{160},\qquad t=\frac{1}{2},x=1\\
u_t&=\frac{41}{16},\qquad t=\frac{1}{2},x=1\\
h&=1,\qquad x=1
\end{aligned}
$$

初期断面 $[x-t,x+t]$ 上の $h$ が与えられれば、その区間を外れた場所で値をどう変えても観測値は変わらない。例えば計算した点の区間は $[1/2,3/2]$ である。この区間の近傍だけに台を持つ滑らかな切断を初期データに掛けても、元の問題と同じ値が出る。エネルギーの議論は無限遠での増大よりも、この局所領域での差を制御している。
:::

### 有界領域で境界を通るエネルギー

滑らかな解 $u$ に対し $u_{tt}-c^2\Delta u=f$ に $u_t$ を掛けて積分すると
$$
\frac{\dd}{\dd t}\int_\Omega\tfrac12(u_t^2+c^2|Du|^2)
=\int_\Omega fu_t+c^2\int_{\partial\Omega}u_t\partial_\nu u .
$$
固定端の斉次 Dirichlet 条件では $u_t=0$、自由端の斉次 Neumann では $\partial_\nu u=0$ なので境界からの仕事は零である。一方、線形な境界ブレーキ $\partial_\nu u+\beta u_t=0$（$\beta\ge0$）では $-c^2\int_{\partial\Omega}\beta u_t^2$ となり、エネルギーが減る。この条件が変位を固定しないことは、空間的に一定な変位が許されることからも分かる。外力がなく $\beta$ を定数とすれば、積分された散逸は初期からその時刻までのエネルギーの減少に等しい。

:::example[一次元の非反射境界]
区間の右端 $x=L$ で $u_x+u_t/c=0$ を課す。$u=F(x-ct)+G(x+ct)$ に入れると $2G'=0$ なので、右から領域へ戻る波の変動が排除される。左端では外法線が逆で $-u_x+u_t/c=0$ とすれば $F'=0$ が入流条件になる。$c=2$、$F'(L-2t)=3$ の瞬間の右端では $u_t=-6$、$u_x=3$ なので境界寄与は $4(-6)3=-72$ である。境界を出る右向き波のエネルギーは再び左向きに領域へ戻らない。
:::

## 非斉次方程式と Duhamel の原理

:::proposition[波動方程式の Duhamel 表示]{#prop-pde-wave-duhamel}
滑らかなデータと外力に対し、

$$
u(t)=\partial_tS(t)g+S(t)h+\int_0^tS(t-s)f(s)\,\dd s
\label{eq-pde-wave-duhamel}
$$

は非斉次問題を解く。
:::

:::proof
積分項を $v$ とする。時間微分では $S(0)=0$ なので端点項は消える。二度目には $S'(0)=I$ により $f(t)$ が現れ、残った積分は $c^2\Delta v$ になる。正則なデータに対して既出の積分公式の領域を固定すれば微分と積分を交換できる。また $v(0)=v_t(0)=0$ である。
:::

:::example[一定入力と共鳴]
全空間の $f=6$ と零初期値では $\int_0^t6(t-s)\,\dd s=3t^2$ となる。一方、$(0,\pi)$ の固定端で $c=1$、$f=\cos t\sin x$ ならモードの係数は $a''+a=\cos t$ を満たす。零初期値の解 $a=t\sin t/2$ は二度微分で確認できる。周波数が一致するため応答の振幅は時間に比例する。
:::

:::example[時間だけの衝撃入力を幅で近似する]
$(0,\pi)$ の $\sin2x$ モードに $f_\delta(t)=\delta^{-1}\mathbf1_{[s,s+\delta]}(t)$ を係数として加える。$t>s+\delta$ では応答の係数は $(4\delta)^{-1}(\cos(2(t-s-\delta))-\cos(2(t-s)))$ である。$\delta\downarrow0$ で $\sin(2(t-s))/2$ に収束する。時刻 $s$ で変位が連続のまま、初速度の係数が $1$ 増える理想化に対応する。
:::

### 半直線の境界入力

$x>0$、$u(0,x)=u_t(0,x)=0$、$u(t,0)=b(t)$ とし、$b$ を負時刻で零に延ばして $C^2$ であるとする。$u(t,x)=b(t-x/c)$ は $u_{tt}-c^2u_{xx}=0$ を満たし、初期値・境界値を実現する。一意性は零境界との差を奇延長し、依存領域の評価を使えば従う。零境界と整合する初期値・外力では、斉次 Dirichlet の解を先に作り、それに境界入力の解を足せる。一般の適合データでは、まず境界値を実現する滑らかな延長を引き、残った方程式の源項と初期値を合わせて調べる。速度 $c$ による遅れを伴うので、境界を時刻 $s$ に変えても点 $x$ へは $s+x/c$ まで届かない。

:::example[立ち上がりの時刻を遅らせた境界値]
$c=3$ とし、$b(t)=(t-1)_+^3$ とする。信号を $x=6$ で測れば $u(t,6)=(t-3)_+^3$ であり、$t=7/2$ で $1/8$ である。境界から領域へ入るエネルギーの率は $-c^2u_tu_x|_{x=0}=c(b')^2$。従って $1<t<2$ の間に入力したエネルギーは $\int_1^2 3\cdot9(t-1)^4\dd t=27/5$ である。時刻 $2$ での分布を積分しても $u_t=3(1-x/3)^2$、$u_x=-u_t/3$ が $0<x<3$ で成り立ち、同じ値になる。
:::

:::example[強制振動の時間因子を初期値と合わせる]
$(0,\pi)$ で $c=1$、零境界、変位と速度は初期に零、外力は $4\cos(2t)\sin(2x)$ とする。共鳴する振動子の解を検算する。

まず初期値を固有関数ごとに分け、対応する時間方程式を解く。外力がある場合も、同じ固有関数に属する外力だけを時間方程式に入れる。以下では得られた表示の微分を取って、式の各項の関係を確かめる。

$$
\begin{aligned}
u&=t \sin{\left(2 t \right)} \sin{\left(2 x \right)}\\
u_{tt}&=\frac{\partial^{2}}{\partial t^{2}} t \sin{\left(2 t \right)} \sin{\left(2 x \right)}=- 4 \left(t \sin{\left(2 t \right)} - \cos{\left(2 t \right)}\right) \sin{\left(2 x \right)}\\
u_{xx}&=\frac{\partial^{2}}{\partial x^{2}} t \sin{\left(2 t \right)} \sin{\left(2 x \right)}=- 4 t \sin{\left(2 t \right)} \sin{\left(2 x \right)}\\
f&=4 t \sin{\left(2 t \right)} \sin{\left(2 x \right)} + 4 \left(- t \sin{\left(2 t \right)} + \cos{\left(2 t \right)}\right) \sin{\left(2 x \right)}=4 \sin{\left(2 x \right)} \cos{\left(2 t \right)}
\end{aligned}
$$

境界条件を各モードで確かめてから、時刻零のデータに戻る。端点値がすべて零でも時間の係数は任意に選べない。時刻を零に戻して初めて各係数の初期条件が確認できる。

$$
\begin{aligned}
u(0,x)&=0,\qquad t=0\\
u_t(0,x)&=0,\qquad t=0\\
u(t,0)&=0,\qquad x=0\\
u(t,L)&=0,\qquad x=\pi
\end{aligned}
$$

エネルギーを空間積分で測る。異なるモードの積の積分は零であり、一つずつ積分して足し合わせられる。その時間微分を取り、解の式を使って積分する方法と比較する。

$$
\begin{aligned}
e(t,x)&=2 t^{2} \sin^{2}{\left(2 t \right)} \cos^{2}{\left(2 x \right)} + \frac{\left(2 t \sin{\left(2 x \right)} \cos{\left(2 t \right)} + \sin{\left(2 t \right)} \sin{\left(2 x \right)}\right)^{2}}{2}\\
\mathcal E(t)&=\frac{\pi \left(4 t^{2} + 2 t \sin{\left(4 t \right)} - \cos^{2}{\left(2 t \right)} + 1\right)}{4}\\
\mathcal E^{\prime}&=\frac{d}{d t} \frac{\pi \left(4 t^{2} + 2 t \sin{\left(4 t \right)} - \cos^{2}{\left(2 t \right)} + 1\right)}{4}=\pi \left(2 t \cos{\left(4 t \right)} + 2 t + \sin{\left(4 t \right)}\right)\\
\text{入力と散逸の和}&=\int\limits_{0}^{\pi} 4 \left(2 t \sin{\left(2 x \right)} \cos{\left(2 t \right)} + \sin{\left(2 t \right)} \sin{\left(2 x \right)}\right) \sin{\left(2 x \right)} \cos{\left(2 t \right)}\, dx=\pi \left(2 t \cos{\left(4 t \right)} + 2 t + \sin{\left(4 t \right)}\right)
\end{aligned}
$$

最後に一つの観測点で変位または温度と勾配を計算する。点で値が相殺されることと、区間全体のエネルギーが零になることは区別する。

$$
\begin{aligned}
u&=\frac{\pi}{4},\qquad t=\frac{\pi}{4},x=\frac{\pi}{4}\\
u_x&=0,\qquad t=\frac{\pi}{4},x=\frac{\pi}{4}\\
\mathcal E&=\frac{\pi \left(1 + \frac{\pi^{2}}{4}\right)}{4}=\frac{\pi \left(4 + \pi^{2}\right)}{16},\qquad t=\frac{\pi}{4}
\end{aligned}
$$

係数 $t\sin(2t)$ を二回微分すると $4\cos(2t)-4t\sin(2t)$ であり、復元力の項が後者を消す。外力が一つの固有関数だけを含むため、他の係数は零である。外力の入力とエネルギーの微分は一致するが、外力の符号が変わるので瞬間的にエネルギーが減る時刻もある。
:::

### 時空間の入力を円錐で積分する

一次元の斉次初期値に Duhamel を使うと、外力応答は
$$
u(t,x)=\frac1{2c}\int_0^t\int_{x-c(t-s)}^{x+c(t-s)} f(s,y)\dd y\dd s
$$
である。ここでは初速度作用素が区間積分なので、時間方向の積分と合わせ過去の円錐の内部全体を積分する。$f$ が連続で微分可能な場合は、内側の積分をまず両端について微分して $u_t$ を求める。もう一回微分すると $s=t$ の上端から $f(t,x)$ が一回現れ、二つの端点の位置微分は $c^2u_{xx}$ と一致する。$f\in L^1_{\rm loc}$ でもこの式が局所で意味を持つ場合は、滑らかな近似と局所積分により分布の式が得られる。任意の観測点に現れるのはその円錐内にある外力だけである。

:::example[時間だけでなく空間にも幅がある入力]
$c=1$、$f(s,y)=\mathbf1_{[0,1]}(s)\mathbf1_{[-1,1]}(y)$、$(t,x)=(2,0)$ とする。$0\le s\le1$ では円錐断面の半径 $2-s\ge1$ なので、入力の全幅 $2$ が入る。従って $u(2,0)=\tfrac12\int_0^1 2\dd s=1$ である。$t=1$ なら断面の半径は $1-s\le1$ で、$u(1,0)=\int_0^1(1-s)\dd s=1/2$。$t=3/2$ では $s=1/2$ で半径が $1$ を横切るため
$$
u(3/2,0)=\int_0^{1/2}1\dd s+\int_{1/2}^1(3/2-s)\dd s=7/8.
$$
空間幅を無視して外力の時間積分だけを使うと、$t<2$ の違いが失われる。
:::

### 三次元の遅延ポテンシャル

三次元で零初期値、外力 $f$ とする。Kirchhoff の初速度作用素を Duhamel へ入れ、時間変数を $r=c(t-s)$ に変えると
$$
\begin{aligned}
u(t,x)&=\int_0^t\frac{1}{4\pi c^2(t-s)}\int_{|y-x|=c(t-s)} f(s,y)\dd S_y\dd s\\
&=\frac{1}{4\pi c^2}\int_{|y-x|<ct}\frac{f(t-|x-y|/c,y)}{|x-y|}\dd y .
\end{aligned}
$$
実際、第一式の面積は $4\pi r^2$、$\dd s=-\dd r/c$、$c^2(t-s)=cr$ である。従って球面を半径で積分すれば第二式になる。外力が観測時刻より前に働いていたとしても、同じ領域のすべての時刻が寄与するのではなく、各位置 $y$ について伝播に要する時間 $|x-y|/c$ だけ前の時刻に評価される。

:::example[球内の定数源が定常な値に近づく]
$c=1$、$f(s,y)=1$（$s>0,|y|<2$）、その他で零とする。滑らかな有界近似で遅延ポテンシャルの表示に優収束を適用できる。中心の応答は $\int_0^{\min(t,2)}r\dd r=\tfrac12\min(t,2)^2$ である。従って $t=1$ では $1/2$、$t\ge2$ では $2$ になる。半径 $2$ の位置からも効果が届いた後には、中心の値は $-\Delta$ の Newton ポテンシャルの値と一致する。ただし全空間では、まだ信号が来ていない場所が任意の有限時刻で残る。
:::

:::example[三角関数の外力に対する非共鳴応答]
区間 $(0,\pi)$ の固定端で $c=1$ とし、第2モードに外力 $3\cos t$、初期変位 $0$、初期速度 $0$ を与える。時間方程式は $a^{\prime\prime}+4a=3\cos t$ である。

時間変換の変数を $s$ とし、$\mathcal L[f](s)=\int_0^\infty e^{-st}f(t)\dd t$ とする。大きい $s$ の範囲で積分が収束する。時間微分を部分積分すると、$t=0$ の値が分子の項になる。まず外力の変換、微分から来る分母、初期データから来る分子を分ける。

$$
\begin{aligned}
f(t)&=3 \cos{\left(t \right)}\\
\widetilde f(s)&=\frac{3 s}{s^{2} + 1}\\
\text{微分の多項式}&=s^{2} + 4\\
\text{初期項}&=0\\
\widetilde a(s)&=\frac{3 s}{s^{4} + 5 s^{2} + 4}
\end{aligned}
$$

次に分母を因数分解する。一つの極が重なるときは単純な指数の係数に加え、時間のべきが掛かる項を求める。例えば $(s+\lambda)^{-2}$ は $te^{-\lambda t}$、分母の二次式 $(s+\lambda)^2+\omega^2$ は減衰する正弦または余弦に対応する。

$$
\begin{aligned}
\widetilde a(s)&=- \frac{s}{s^{2} + 4} + \frac{s}{s^{2} + 1}\\
a(t)&=\cos{\left(t \right)} - \cos{\left(2 t \right)}\\
a^{\prime}&=\frac{d}{d t} \left(\cos{\left(t \right)} - \cos{\left(2 t \right)}\right)=\left(4 \cos{\left(t \right)} - 1\right) \sin{\left(t \right)}\\
a^{\prime\prime}&=\frac{d^{2}}{d t^{2}} \left(\cos{\left(t \right)} - \cos{\left(2 t \right)}\right)=- \cos{\left(t \right)} + 4 \cos{\left(2 t \right)}
\end{aligned}
$$

最後に時刻を戻して初期条件を調べ、次に得た時間関数を元の左辺に代入する。変換領域の式を満たしても、逆変換で定数や初期項を落とすとこの二つの確認のどちらかに失敗する。

$$
\begin{aligned}
a(0)&=0,\qquad t=0\\
a^{\prime}(0)&=0,\qquad t=0\\
\text{時間方程式の左辺}&=3 \cos{\left(t \right)}\\
a(t_*)&=1,\qquad t=\frac{\pi}{2}\\
a^{\prime}(t_*)&=-1,\qquad t=\frac{\pi}{2}
\end{aligned}
$$

得られる $a(t)=\cos t-\cos2t$ は、外力と同じ振動数の特解に固有振動を足して初期変位を零にした形である。両振動数は異なるため $t$ の線形成長はなく、時刻 $\pi/2$ で値は $1$ となる。
:::

## 有界領域における固有関数展開

斉次 Dirichlet 条件を持つ $-\Delta$ の正規直交固有系を $\phi_j$、固有値を $\lambda_j>0$ とする。存在は [楕円型方程式の弱解](/math/pde/weak-solutions/) を参照する。$\omega_j=c\sqrt{\lambda_j}$ とすると $u_j''+\omega_j^2u_j=f_j$ であり、

$$
u_j(t)=g_j\cos(\omega_jt)+h_j\frac{\sin(\omega_jt)}{\omega_j}
+\int_0^t\frac{\sin(\omega_j(t-s))}{\omega_j}f_j(s)\,\dd s
\label{eq-pde-wave-eigen}
$$

である。各係数を微分すれば ODE と初期値が確認できる。斉次問題で各モードのエネルギーは一定であり、$\sum_j(\lambda_j|g_j|^2+|h_j|^2)<\infty$ なら部分和はエネルギー空間で時間について一様 Cauchy である。外力が $L^1_tL^2_x$ に入る場合はエネルギーの平方根の増加を $\norm f_2$ の積分で抑えると、外力の部分和も収束する。[双曲型方程式](/math/pde/hyperbolic/) ではこの極限を弱解として扱う。

:::example[矩形の重複固有値と Neumann の零固有値]
$(0,\pi)^2$ の固定端では $\sin(mx_1)\sin(nx_2)$ の固有値は $m^2+n^2$ である。$(m,n)=(1,2),(2,1)$ は $c=2$ ならともに角振動数 $2\sqrt5$ を持つ。他方、Neumann の定数モードでは $\lambda=0$ である。式の $\sin(\omega t)/\omega$ を極限 $t$ で読むと、一定な初期条件 $g=1$、$h=2$ から $u=1+2t$ が得られる。
:::

:::example[片側だけ自由な弦]
区間 $(0,L)$ で左を固定端、右を自由端とすると $X(0)=0,X'(L)=0$ である。固有関数は $\sin(kx)$ で $kL=(j+1/2)\pi$（$j\ge0$）となる。$L=\pi$、$c=2$ では最低角振動数は $1$ であり、両端固定の最低角振動数 $2$ より小さい。反射を使う場合、長さ $2L$ の平行移動で符号が反転し、同じ波数が得られる。
:::

### エネルギー空間での級数収束

対称な正の Dirichlet Laplace 作用素の固有値を $\lambda_j$、正規直交基底を $\phi_j$ とする。$g\in H_0^1$ と $h\in L^2$ なら $\sum_j(\lambda_j|g_j|^2+|h_j|^2)<\infty$ である。各振動子のエネルギー $|u_j'|^2+c^2\lambda_j|u_j|^2$ は時間によらない。それを $j>N$ について足すと $u_N\to u$ は $C([0,T];H_0^1)$、$u_N'\to u_t$ は $C([0,T];L^2)$ で収束する。弱式を固定した試験関数に対して有限和から極限へ移せる。対応する外力 $f\in L^1_tL^2_x$ を加える場合、有限和の Duhamel 部分のエネルギーノルムは $\int_0^t\norm{P_Nf(s)}_2\dd s$ 以下である。従って射影の差にも同じ評価を用いれば、優収束により全級数が同じ空間で収束する。

:::example[二階微分が二乗可積分でない初期変位]
$(0,\pi)$ で正規直交正弦基底の係数を $g_j=1/j^2$、$h_j=0$ とする。$\lambda_j=j^2$ なので勾配エネルギーは $\sum 1/j^2$、$\norm{\Delta g}_2^2$ に必要な和は $\sum 1$ である。従ってエネルギー解はあるが、原点時刻に二階空間微分を $L^2$ に要求できない。$c=1$、$t=2\pi$ では全ての $\cos(jt)$ が $1$ なので、初期変位の正則性のまま戻る。熱方程式のように正時刻で高周波を減衰させる仕組みはない。
:::

:::example[二つの速度モードが周期的に戻る]
$(0,\pi)$ で $c=2$、初期変位零、初期速度 $4\sin x+6\sin3x$ とする。

まず初期値を固有関数ごとに分け、対応する時間方程式を解く。外力がある場合も、同じ固有関数に属する外力だけを時間方程式に入れる。以下では得られた表示の微分を取って、式の各項の関係を確かめる。

$$
\begin{aligned}
u&=2 \sin{\left(2 t \right)} \sin{\left(x \right)} + \sin{\left(6 t \right)} \sin{\left(3 x \right)}\\
u_{tt}&=\frac{\partial^{2}}{\partial t^{2}} \left(2 \sin{\left(2 t \right)} \sin{\left(x \right)} + \sin{\left(6 t \right)} \sin{\left(3 x \right)}\right)=- 4 \left(2 \sin{\left(2 t \right)} \sin{\left(x \right)} + 9 \sin{\left(6 t \right)} \sin{\left(3 x \right)}\right)\\
u_{xx}&=\frac{\partial^{2}}{\partial x^{2}} \left(2 \sin{\left(2 t \right)} \sin{\left(x \right)} + \sin{\left(6 t \right)} \sin{\left(3 x \right)}\right)=- 2 \sin{\left(2 t \right)} \sin{\left(x \right)} - 9 \sin{\left(6 t \right)} \sin{\left(3 x \right)}\\
f&=0
\end{aligned}
$$

境界条件を各モードで確かめてから、時刻零のデータに戻る。端点値がすべて零でも時間の係数は任意に選べない。時刻を零に戻して初めて各係数の初期条件が確認できる。

$$
\begin{aligned}
u(0,x)&=0,\qquad t=0\\
u_t(0,x)&=4 \sin{\left(x \right)} + 6 \sin{\left(3 x \right)}=2 \left(2 \sin{\left(x \right)} + 3 \sin{\left(3 x \right)}\right),\qquad t=0\\
u(t,0)&=0,\qquad x=0\\
u(t,L)&=0,\qquad x=\pi
\end{aligned}
$$

エネルギーを空間積分で測る。異なるモードの積の積分は零であり、一つずつ積分して足し合わせられる。その時間微分を取り、解の式を使って積分する方法と比較する。

$$
\begin{aligned}
e(t,x)&=2 \left(2 \sin{\left(2 t \right)} \cos{\left(x \right)} + 3 \sin{\left(6 t \right)} \cos{\left(3 x \right)}\right)^{2} + \frac{\left(4 \sin{\left(x \right)} \cos{\left(2 t \right)} + 6 \sin{\left(3 x \right)} \cos{\left(6 t \right)}\right)^{2}}{2}\\
\mathcal E(t)&=13 \pi\\
\mathcal E^{\prime}&=\frac{d}{d t} 13 \pi=0\\
\text{入力と散逸の和}&=\int\limits_{0}^{\pi} 0\, dx=0
\end{aligned}
$$

最後に一つの観測点で変位または温度と勾配を計算する。点で値が相殺されることと、区間全体のエネルギーが零になることは区別する。

$$
\begin{aligned}
u&=3,\qquad t=\frac{\pi}{4},x=\frac{\pi}{2}\\
u_x&=0,\qquad t=\frac{\pi}{4},x=\frac{\pi}{2}\\
\mathcal E&=13 \pi,\qquad t=\frac{\pi}{4}
\end{aligned}
$$

時刻 $\pi/4$ では速度は両モードとも零だが変位は零ではない。従ってエネルギーは全て勾配の項に入る。時刻 $\pi$ では変位は零、速度は初期値と一致する。固有周波数は整数比なので共通の帰還時刻を持つ。
:::

:::example[物理変数と固有モードのエネルギーを比べる]
二つの同じ質量の振動子を考え、剛性行列を $\left(\begin{smallmatrix}5&-4\\-4&5\end{smallmatrix}\right)$ とする。これは二つの空間モードだけを取った波を表す模型でもある。

各行の質量係数と剛性係数を区別して、候補を連立方程式に代入する。$K$ の非対角成分は別の成分を結び付ける。従って個々の微分の値を求めるだけではなく、行列を掛けた結果も確認する。

$$
\begin{aligned}
M&=\left[\begin{matrix}1 & 0\\0 & 1\end{matrix}\right]\\
K&=\left[\begin{matrix}5 & -4\\-4 & 5\end{matrix}\right]\\
U&=\left[\begin{matrix}\cos{\left(t \right)} + \cos{\left(3 t \right)}\\\cos{\left(t \right)} - \cos{\left(3 t \right)}\end{matrix}\right]\\
U_{tt}&=\left[\begin{matrix}- (\cos{\left(t \right)} + 9 \cos{\left(3 t \right)})\\- \cos{\left(t \right)} + 9 \cos{\left(3 t \right)}\end{matrix}\right]\\
F&=\left[\begin{matrix}0\\0\end{matrix}\right]
\end{aligned}
$$

初期時刻へ戻して、先に指定したデータと比較する。次に正の時刻における成分の値を読む。成分が混ざっていても、初期データは方程式を積分する前の元の座標で指定される。

$$
\begin{aligned}
U(0)&=\left[\begin{matrix}2\\0\end{matrix}\right],\qquad t=0\\
U_t(0)&=\left[\begin{matrix}0\\0\end{matrix}\right],\qquad t=0\\
U&=\left[\begin{matrix}0\\\sqrt{2}\end{matrix}\right],\qquad t=\frac{\pi}{4}\\
U_t&=\left[\begin{matrix}- 2 \sqrt{2}\\\sqrt{2}\end{matrix}\right],\qquad t=\frac{\pi}{4}
\end{aligned}
$$

二次量を具体的に求め、微分した値と流束または外力の寄与を照合する。質量行列が単位行列でない場合は係数の普通の二乗和がエネルギーになるとは限らない。

$$
\begin{aligned}
e&=\tfrac12(U_t^\top M U_t+U^\top K U)\\
e_t&=U_t^\top(MU_{tt}+KU)=0\\
e(0)&=\tfrac12(2,0)K(2,0)^\top=10
\end{aligned}
$$

同相の固有方向は固有値 $1$、逆相の方向は固有値 $9$ を持つので、振動数は $1$ と $3$ である。例えば $t=\pi/4$ で第一成分の変位は零に打ち消される一方、速度と第二成分に運動は残る。従って一つの位置で零になってもエネルギーは初期値と同じである。
:::

:::example[固有周波数が零になる Neumann モード]
両端 Neumann の定数モードに外力 $f=6t$、初期変位 $1$、初期速度 $2$ を与える。時間方程式は $a^{\prime\prime}=6t$ である。

時間変換の変数を $s$ とし、$\mathcal L[f](s)=\int_0^\infty e^{-st}f(t)\dd t$ とする。大きい $s$ の範囲で積分が収束する。時間微分を部分積分すると、$t=0$ の値が分子の項になる。まず外力の変換、微分から来る分母、初期データから来る分子を分ける。

$$
\begin{aligned}
f(t)&=6 t\\
\widetilde f(s)&=\frac{6}{s^{2}}\\
\text{微分の多項式}&=s^{2}\\
\text{初期項}&=s + 2\\
\widetilde a(s)&=\frac{s^{3} + 2 s^{2} + 6}{s^{4}}
\end{aligned}
$$

次に分母を因数分解する。一つの極が重なるときは単純な指数の係数に加え、時間のべきが掛かる項を求める。例えば $(s+\lambda)^{-2}$ は $te^{-\lambda t}$、分母の二次式 $(s+\lambda)^2+\omega^2$ は減衰する正弦または余弦に対応する。

$$
\begin{aligned}
\widetilde a(s)&=\frac{1}{s} + \frac{2}{s^{2}} + \frac{6}{s^{4}}\\
a(t)&=t^{3} + 2 t + 1\\
a^{\prime}&=\frac{d}{d t} \left(t^{3} + 2 t + 1\right)=3 t^{2} + 2\\
a^{\prime\prime}&=\frac{d^{2}}{d t^{2}} \left(t^{3} + 2 t + 1\right)=6 t
\end{aligned}
$$

最後に時刻を戻して初期条件を調べ、次に得た時間関数を元の左辺に代入する。変換領域の式を満たしても、逆変換で定数や初期項を落とすとこの二つの確認のどちらかに失敗する。

$$
\begin{aligned}
a(0)&=1,\qquad t=0\\
a^{\prime}(0)&=2,\qquad t=0\\
\text{時間方程式の左辺}&=6 t\\
a(t_*)&=4,\qquad t=1\\
a^{\prime}(t_*)&=5,\qquad t=1
\end{aligned}
$$

零固有値では固有周波数で割る表示の極限が必要である。外力の Duhamel 核は $\sin(\omega (t-s))/\omega\to t-s$ となり、$\int_0^t(t-s)6s\dd s=t^3$ が得られる。初期速度は $2t$ を加える。勾配エネルギーは零でも運動エネルギーは残る。
:::

## 演習問題

::::exercise{level=basic}
一次元で $c=1$、$g=0$、$h(x)=x$ の解を求めよ。

:::solution
積分は $((x+t)^2-(x-t)^2)/4=xt$ である。$u_{tt}=u_{xx}=0$ で、$u_t(0,x)=x$ である。
:::
::::

::::exercise{level=basic}
$(0,\pi)$ の固定端で $c=2$、初期値 $(\sin3x,0)$ の解を求めよ。

:::solution
固有値は $9$ で $\omega=6$ だから $u=\cos6t\sin3x$ である。端点と初速度の条件を満たす。
:::
::::

::::exercise{level=basic}
三次元で $c=1$、$g=0$、$h=|x|^2$ の解を求めよ。

:::solution
球面平均が $|x|^2+t^2$ なので $u=t|x|^2+t^3$ である。$u_{tt}=\Delta u=6t$ が成立する。
:::
::::

::::exercise{level=basic}
一次元で初期台が $[-2,2]$、$c=3$ のとき $x=8$ に到達し得る最初の時刻はいつか。

:::solution
最短距離は $8-2=6$ である。円錐が初期台に接する最初の時刻は $6/3=2$ である。
:::
::::

::::exercise{level=advanced}
$c=1$、$(0,\pi)$ の固定端で $u=\sin t\sin x$ の全エネルギーを求めよ。

:::solution
$\int_0^\pi\sin^2x=\int_0^\pi\cos^2x=\pi/2$ なのでエネルギーは $\pi(\sin^2t+\cos^2t)/4=\pi/4$ である。
:::
::::

::::exercise{level=advanced}
一次元で外力 $f=t$ と零初期値から生じる波を求めよ。

:::solution
空間一定の初速度 $s$ は経過時間 $t-s$ で $(t-s)s$ の変位を与える。積分すると $u=t^3/6$ であり、$u_{tt}=t$ と零初期条件を満たす。
:::
::::

::::exercise{level=advanced}
円錐境界のエネルギー流束が $0$ になる条件を求めよ。

:::solution
:eqref[eq-pde-wave-cone-flux] の二つの平方が消える必要十分条件は $u_t=c\partial_\nu u$ と $D_{\rm tan}u=0$ である。一次元の右側境界では $u=A(x+ct)$ が例になる。
:::
::::

::::exercise{level=advanced}
二次元で $g=0$、$h=x_1$ の解を Poisson 公式で求めよ。

:::solution
$y=x+z$ とおくと $z_1$ の寄与は対称性で消える。定数初速度の規格化で残りは $tx_1$ である。
:::
::::

## 確認テスト

::::quiz{answer="6" tol="1e-6"}
$c=3$、Laplace 固有値 $4$ の角振動数を求めよ。

:::explanation
係数は $a''+c^2\lambda a=0$ を満たすので $3\sqrt4=6$ である。
:::
::::

::::quiz
三次元の自由波の零変位からの初速度項はどこを積分するか。

- [x] 半径 $ct$ の球面
- [ ] 半径 $ct$ の球の内部
- [ ] 全空間

:::explanation
Kirchhoff の初速度項は $tM_{ct}h$ であり、球面平均である。二次元へ降下させた場合は円盤の積分になる。
:::
::::

::::quiz{answer="2" tol="1e-6"}
一次元で $c=2$、初期台が $[-1,1]$ のとき $x=5$ に到達し得る最初の時刻を求めよ。

:::explanation
初期台からの距離は $4$ なので速度 $2$ で割って $2$ である。
:::
::::

::::quiz
固定端の自由波の全エネルギーはどう変化するか。

- [x] 保たれる
- [ ] 常に減少する
- [ ] 常に増加する

:::explanation
境界で $u_t=0$ なので境界流束 $c^2u_t\partial_\nu u$ が消える。外力がなければ内部に仕事項もなく、時間微分は $0$ である。
:::
::::

::::quiz{answer="0.875" tol="1e-6"}
$u_{tt}-u_{xx}=f$、初期変位・速度は零とする。$f(t,x)=1$ は $0<t<1$ かつ $|x|<1$ の範囲、他は零である。$u(3/2,0)$ を求めよ。

:::explanation
過去の時刻 $s$ における円錐の断面は半径 $3/2-s$ である。Duhamel の式は $0<s<1$ で、これと入力区間が重なる長さの半分を積分する。$s=1/2$ までは半幅 $1$、その後は $3/2-s$ なので、$1/2+\int_{1/2}^1(3/2-s)\dd s=1/2+3/8=7/8$ である。
:::
::::


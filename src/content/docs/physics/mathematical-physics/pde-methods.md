---
title: "偏微分方程式の解法"
description: "偏微分方程式の主部を分類し、領域とデータに合わせて分離、展開、変換、特性曲線を用いる。"
notation:
  - symbol: '\Delta'
    meaning: "ユークリッド空間の Laplacian"
  - symbol: '\kappa'
    meaning: "拡散係数"
  - symbol: 'c'
    meaning: "波動の伝播速度"
  - symbol: '\partial_n'
    meaning: "外向き法線微分"
---

## この章の目標

- 二階方程式の主部と与えるべきデータを区別する。
- 直交・円筒・球座標で空間モードを分離する。
- 時間または空間の変換で非斉次方程式を解く。
- 特性曲線と放射条件から解の選択を説明する。

## 2 階線形偏微分方程式の分類と境界条件

### 主部の符号と境界で指定する量

二変数の二階線形方程式の主部を $A u_{xx}+2B u_{xy}+C u_{yy}$ とする。点を固定した係数行列は実対称で、線形座標変換で実対角化できる。一般の滑らかな変数変換でも、二階項の係数は Jacobian による合同変換になる。連鎖律で Jacobian を微分した分は一階項に入るためである。従って慣性指数で分類でき、$B^2-AC<0$ は楕円型、$>0$ は双曲型、$=0$ で主部が非零なら放物型である。この最後の呼び方は主部の階数の分類であり、前進熱方程式の性質を得るには時間の一階項も確認する。例えば $\Delta u=0$、$u_{tt}-c^2u_{xx}=0$、$u_t-\kappa u_{xx}=0$ は順に三型の代表である。波速 $c$ の単位は $\mathrm{m\,s^{-1}}$、拡散係数 $\kappa$ は $\mathrm{m^2\,s^{-1}}$ である。

境界上で値を指定する条件を Dirichlet、外向き法線微分を指定する条件を Neumann、$\partial_n u+\alpha u=g$ を Robin と呼ぶ。例えば有界連結領域 $\Omega$ の境界が区分的に滑らかで、$\Delta u=f$、$\partial_n u=g$ を課すとする。境界まで十分滑らかな解には、発散定理から $\int_\Omega f\,\dd x=\int_{\partial\Omega}g\,\dd S$ が必要である。二解の差 $v$ は $\Delta v=0$、$\partial_n v=0$ だから
$$
\int_\Omega|\grad v|^2\,\dd x
=\int_{\partial\Omega}\bar v\partial_n v\,\dd S-\int_\Omega \bar v\Delta v\,\dd x=0
$$
となり、連結性から定数になる。Dirichlet 問題の差にも同じ等式を使え、境界値が零なのでその定数も零である。時間発展では波動方程式に初期変位と初期速度、熱方程式に初期温度を与える。この違いは、後で分離したモードの時間方程式の階数に現れる。

### 双曲型の特性方向と波のエネルギー

二階の主部が $A\partial_x^2+2B\partial_x\partial_y+C\partial_y^2$ のとき、新しい座標の一つを $\xi(x,y)$ とすると、$\partial_\xi^2$ の係数は $A\xi_x^2+2B\xi_x\xi_y+C\xi_y^2$ である。これが零となる等位線を特性曲線という。例えば $u_{tt}-c^2u_{xx}=0$ は $\xi=x-ct,\eta=x+ct$ で両方の純二階項が消え、$-4c^2 u_{\xi\eta}=0$ になる。積分して $u=F(\xi)+G(\eta)$ と表せる。これは積分変換で左右に進む解が現れることを局所的に説明する。

有限区間の両端を固定して波動方程式に $u_t$ を掛けると、$E(t)=\tfrac12\int(u_t^2+c^2u_x^2)\dd x$ の微分は $c^2[u_tu_x]_a^b=0$ である。固定端で $u_t=0$ になるためである。初期データがともに零なら $E(0)=0$ だから $u_x,u_t$ が零、さらに固定端から $u=0$ である。従ってこの条件のもとで初期値解は一意である。

### 主部の固有軸へ移る

実対称な主部行列を直交変換で対角化すれば混合項が消える。固有値の符号を読み、一方が零の場合は階数が落ちる方向を調べる。例えば $z=(\xi,\eta)$ の Jacobian を $J_{ai}=\partial_i z_a$ と書くと、新しい主部行列は $J\begin{pmatrix}A&B\\B&C\end{pmatrix}J^{\mathsf T}$ である。$\partial_i\partial_j z_a$ が掛かる項は一階項に入る。

双曲型が近傍全体で成り立ち、係数が滑らかなら、主部の二つの零方向から独立な特性座標 $\xi,\eta$ を作れる。それぞれの方向に接する曲線を局所的に解き、横断する曲線上の出発点を座標にすれば、各座標の勾配が対応する零方向の共変ベクトルになる。対角成分が両方零で行列式が非零なので、混合成分で割ると主部は $u_{\xi\eta}$ になる。定数係数の場合、対応する特性座標を線形に取ることができ、放物型の主部は $u_{\xi\xi}$ に、楕円型の主部は $u_{\xi\xi}+u_{\eta\eta}$ にできる。楕円型では二固有値の絶対値の平方根で各固有軸の座標を割り、必要なら式全体の符号を反転する。

:::example[型が変わる Tricomi 方程式]
$y u_{xx}+u_{yy}=0$ の判別式は $-y$ である。従って $y>0$ では楕円型、$y<0$ では双曲型であり、$y=0$ で主部の階数が落ちる。$y<0$ で
$$
\xi=x+\frac23(-y)^{3/2},\qquad
\eta=x-\frac23(-y)^{3/2}
$$
と置けば $\xi_y=-\sqrt{-y}$、$\eta_y=\sqrt{-y}$ である。純二階項は消えるが、変数変換後の式は
$$
4y u_{\xi\eta}+\frac{u_\xi-u_\eta}{2\sqrt{-y}}=0
$$
となる。座標の Jacobian の行列式 $2\sqrt{-y}$ が零に近づくため、この表示を型の境界まで延ばすことはできない。
:::

:::example[Robin 条件と放熱]
有界連結領域で $u_t=\kappa\Delta u$、$\partial_nu+\alpha u=0$、$\alpha\ge0$ とする。滑らかな実解に対して部分積分すると
$$
\frac{\dd}{\dd t}\frac12\int_\Omega u^2\,\dd x
=-\kappa\int_\Omega|\grad u|^2\,\dd x
-\kappa\int_{\partial\Omega}\alpha u^2\,\dd S
$$
を得る。二解の差にこの式を適用すれば、同じ初期値に対する一意性が従う。定常解は勾配が零の定数であり、$\alpha>0$ の境界部分が正の面積を持てば、その定数も零になる。
:::

:::example[二方向の勾配と主部（1）]
$u_{xx}+2m u_{xy}+u_{yy}$ の主部を $X=(x+y)/\sqrt2,Y=(x-y)/\sqrt2$ で対角化する。

$$
m=1
$$

対角化の前の行列を並べる。

$$
\begin{aligned}
A\big|_{m=1} &= \begin{pmatrix}1 & 1\\1 & 1\end{pmatrix}
\end{aligned}
$$

二つの固有軸の主係数を求める。

$$
\begin{aligned}
(\lambda_X,\lambda_Y)\big|_{m=1} &= \left( 2, \  0\right)
\end{aligned}
$$

判別式は混合係数の半分を二乗してから対角係数の積を引く。

$$
\begin{aligned}
B^2-AC\big|_{m=1} &= 0
\end{aligned}
$$
:::

:::example[二方向の勾配と主部（2）]
$u_{xx}+2m u_{xy}+u_{yy}$ の主部を $X=(x+y)/\sqrt2,Y=(x-y)/\sqrt2$ で対角化する。

$$
m=2
$$

対角化の前の行列を並べる。

$$
\begin{aligned}
A\big|_{m=2} &= \begin{pmatrix}1 & 2\\2 & 1\end{pmatrix}
\end{aligned}
$$

二つの固有軸の主係数を求める。

$$
\begin{aligned}
(\lambda_X,\lambda_Y)\big|_{m=2} &= \left( 3, \  -1\right)
\end{aligned}
$$

判別式は混合係数の半分を二乗してから対角係数の積を引く。

$$
\begin{aligned}
B^2-AC\big|_{m=2} &= 3
\end{aligned}
$$
:::

:::example[二方向の勾配と主部（3）]
$u_{xx}+2m u_{xy}+u_{yy}$ の主部を $X=(x+y)/\sqrt2,Y=(x-y)/\sqrt2$ で対角化する。

$$
m=3
$$

対角化の前の行列を並べる。

$$
\begin{aligned}
A\big|_{m=3} &= \begin{pmatrix}1 & 3\\3 & 1\end{pmatrix}
\end{aligned}
$$

二つの固有軸の主係数を求める。

$$
\begin{aligned}
(\lambda_X,\lambda_Y)\big|_{m=3} &= \left( 4, \  -2\right)
\end{aligned}
$$

判別式は混合係数の半分を二乗してから対角係数の積を引く。

$$
\begin{aligned}
B^2-AC\big|_{m=3} &= 8
\end{aligned}
$$
:::

:::example[二方向の勾配と主部（4）]
$u_{xx}+2m u_{xy}+u_{yy}$ の主部を $X=(x+y)/\sqrt2,Y=(x-y)/\sqrt2$ で対角化する。

$$
m=4
$$

対角化の前の行列を並べる。

$$
\begin{aligned}
A\big|_{m=4} &= \begin{pmatrix}1 & 4\\4 & 1\end{pmatrix}
\end{aligned}
$$

二つの固有軸の主係数を求める。

$$
\begin{aligned}
(\lambda_X,\lambda_Y)\big|_{m=4} &= \left( 5, \  -3\right)
\end{aligned}
$$

判別式は混合係数の半分を二乗してから対角係数の積を引く。

$$
\begin{aligned}
B^2-AC\big|_{m=4} &= 15
\end{aligned}
$$
:::

## 変数分離法（直交座標・円筒座標・球座標）

### 座標の境界に沿って分離する

直方体の Dirichlet 固有問題で $u=X(x)Y(y)Z(z)$ を $\Delta u=-\lambda u$ に入れ、積で割ると $X''/X+Y''/Y+Z''/Z=-\lambda$ となる。各項は他の変数を固定して変化させても一定だから分離定数である。例えば辺長 $(a,b,d)$ では $\sin(n\pi x/a)\sin(m\pi y/b)\sin(\ell\pi z/d)$ の固有値は $(n\pi/a)^2+(m\pi/b)^2+(\ell\pi/d)^2$ である。奇拡張の Fourier 級数の完全性を各変数に適用すると積基底が完全になる。

円筒では $\Delta$ の角と軸を $e^{im\phi}$ と $\sin(n\pi z/L)$ に分け、側壁 $r=a$ で零なら正則動径部 $J_{|m|}(j_{|m|,s}r/a)$ を取る。ここで $m\in\Z$、$n,s\ge1$ で、$j_{\nu,s}$ は $J_\nu$ の第 $s$ 正零点である。円筒座標の Laplacian は式 :eqref[eq-mp-orthogonal-operators] に尺度因子 $1,r,1$ を入れて得られ、動径方程式は
$$
(rR')'+\left[k_\perp^2r-\frac{m^2}{r}\right]R=0,\qquad k_\perp^2=\lambda-(n\pi/L)^2
$$
が従う。両辺に $\bar R$ を掛けて部分積分すると、$k_\perp^2\int_0^a r|R|^2$ は $\int_0^a(r|R'|^2+m^2|R|^2/r)$ に等しい。原点では正則性、側壁では Dirichlet 条件により端点項が消える。非自明解では右辺は正だから $k_\perp>0$ である。第二種の解は $r=0$ で発散し、第一種の解に側壁条件を課すと $k_\perp=j_{|m|,s}/a$ となる。従って固有値は $(j_{|m|,s}/a)^2+(n\pi/L)^2$ である。

球では :ref[def-mp-spherical-harmonic] で角部分を対角化する。$\Delta_{S^2}Y_{\ell m}=-\ell(\ell+1)Y_{\ell m}$ なので、$\lambda=k^2>0$ に対して
$$
(r^2R')'+[k^2r^2-\ell(\ell+1)]R=0
$$
を得る。原点で正則な $j_\ell(kr)$ を選び、外端で $j_\ell(ka)=0$ とする。例えば球の $\ell=0$ の最初のモードは $\sin(\pi r/a)/(\pi r/a)$ で、$r=0$ での極限は $1$、$r=a$ での値は零である。境界が座標面と一致しなければ、積表示を境界にそのまま入れて分離定数を決められない場合がある。

:::example[球面温度から内部の調和関数を求める]
$\Delta u=0$ の球内部では動径方程式を $k=0$ として解く。$R=r^\beta$ を代入すると $\beta(\beta+1)=\ell(\ell+1)$ であり、二つの解 $r^\ell,r^{-\ell-1}$ のうち原点で正則な前者を選ぶ。半径 $a$ の球面で $g=2Y_{00}-Y_{20}$ を指定すると、解は $2Y_{00}-(r/a)^2Y_{20}$ である。一般の有限和 $g=\sum g_{\ell m}Y_{\ell m}$ にも同じ計算が使える。外部問題で無限遠に零を課すなら動径因子は $(a/r)^{\ell+1}$ であり、同じデータは $2(a/r)Y_{00}-(a/r)^3Y_{20}$ を与える。
:::

### 規格化した分離解

各座標方向の固有関数は、境界条件を満たしたものを選ぶ。長方形で四辺 Dirichlet なら正弦の積を使える。直積領域でのノルムは各方向のノルムの積である。固有値は和になり、二つの添字の組が同じ和を与える場合は縮退になる。円筒では角方向の整数と Bessel の零点の組、球では角運動量の次数と球 Bessel の零点の組が同じ役割を持つ。

円板の動径因子の内積には面積要素に由来する重み $r$ が付く。$R_\rho(r)=J_\nu(\rho r)$ とし、異なる $\rho,\tau$ の方程式を相互に掛けて引くと
$$
(\rho^2-\tau^2)\int_0^a rR_\rho R_\tau\,\dd r
=[r(R_\rho R_\tau'-R_\rho'R_\tau)]_0^a
$$
を得る。$\rho,\tau$ がともに Dirichlet 波数なら右辺が消える。$\rho=j_{\nu,s}/a$ を固定し、$\tau$ を $\rho$ へ近づければ
$$
\int_0^a r J_\nu(\rho r)^2\,\dd r=\frac{a^2}{2}J_\nu'(j_{\nu,s})^2
$$
となる。右辺の端点差は $-a R_\rho'(a)R_\tau(a)$ であり、$R_\tau(a)=a(\tau-\rho)J_\nu'(j_{\nu,s})+o(\tau-\rho)$ を代入すればよい。球の場合は同様に重み $r^2$ を使う。従って境界の形だけでなく、積分測度も分離係数に影響する。

:::example[長方形で一方の量子数を変える（1）]
横の辺長 $2$、縦の辺長 $1$ とし、縦の第 1 モードに固定する。

$$
m=1
$$

横の Dirichlet 条件から正弦の波数が決まる。

$$
\begin{aligned}
k_x\big|_{m=1} &= \frac{\pi}{2}
\end{aligned}
$$

固有関数は二つの正弦の積であり、規格化に $\sqrt2$ を掛ける。

$$
\begin{aligned}
\phi_{m1}\big|_{m=1} &= \sqrt{2} \sin{\left(\frac{\pi x}{2} \right)} \sin{\left(\pi y \right)}
\end{aligned}
$$

$-\Delta$ の固有値は波数の二乗の和である。

$$
\begin{aligned}
\lambda_{m1}\big|_{m=1} &= \frac{5 \pi^{2}}{4}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
\phi_{m1}\big|_{x=\frac{1}{3}} &= \frac{\sqrt{2} \sin{\left(\pi y \right)}}{2} \\
\frac{\partial^1 (\phi_{m1})}{\partial x^1} &= \frac{\sqrt{2} \pi \sin{\left(\pi y \right)} \cos{\left(\frac{\pi x}{2} \right)}}{2} \\
\left.\frac{\partial^1 (\phi_{m1})}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{\sqrt{6} \pi \sin{\left(\pi y \right)}}{4} \\
\frac{\partial^2 (\phi_{m1})}{\partial x^2} &= - \frac{\sqrt{2} \pi^{2} \sin{\left(\frac{\pi x}{2} \right)} \sin{\left(\pi y \right)}}{4} \\
\left.\frac{\partial^2 (\phi_{m1})}{\partial x^2}\right|_{x=\frac{1}{3}} &= - \frac{\sqrt{2} \pi^{2} \sin{\left(\pi y \right)}}{8}
\end{aligned}
$$
:::

:::example[長方形で一方の量子数を変える（2）]
横の辺長 $2$、縦の辺長 $1$ とし、縦の第 1 モードに固定する。

$$
m=2
$$

横の Dirichlet 条件から正弦の波数が決まる。

$$
\begin{aligned}
k_x\big|_{m=2} &= \pi
\end{aligned}
$$

固有関数は二つの正弦の積であり、規格化に $\sqrt2$ を掛ける。

$$
\begin{aligned}
\phi_{m1}\big|_{m=2} &= \sqrt{2} \sin{\left(\pi x \right)} \sin{\left(\pi y \right)}
\end{aligned}
$$

$-\Delta$ の固有値は波数の二乗の和である。

$$
\begin{aligned}
\lambda_{m1}\big|_{m=2} &= 2 \pi^{2}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
\phi_{m1}\big|_{x=\frac{1}{3}} &= \frac{\sqrt{6} \sin{\left(\pi y \right)}}{2} \\
\frac{\partial^1 (\phi_{m1})}{\partial x^1} &= \sqrt{2} \pi \sin{\left(\pi y \right)} \cos{\left(\pi x \right)} \\
\left.\frac{\partial^1 (\phi_{m1})}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{\sqrt{2} \pi \sin{\left(\pi y \right)}}{2} \\
\frac{\partial^2 (\phi_{m1})}{\partial x^2} &= - \sqrt{2} \pi^{2} \sin{\left(\pi x \right)} \sin{\left(\pi y \right)} \\
\left.\frac{\partial^2 (\phi_{m1})}{\partial x^2}\right|_{x=\frac{1}{3}} &= - \frac{\sqrt{6} \pi^{2} \sin{\left(\pi y \right)}}{2}
\end{aligned}
$$
:::

:::example[長方形で一方の量子数を変える（3）]
横の辺長 $2$、縦の辺長 $1$ とし、縦の第 1 モードに固定する。

$$
m=3
$$

横の Dirichlet 条件から正弦の波数が決まる。

$$
\begin{aligned}
k_x\big|_{m=3} &= \frac{3 \pi}{2}
\end{aligned}
$$

固有関数は二つの正弦の積であり、規格化に $\sqrt2$ を掛ける。

$$
\begin{aligned}
\phi_{m1}\big|_{m=3} &= \sqrt{2} \sin{\left(\frac{3 \pi x}{2} \right)} \sin{\left(\pi y \right)}
\end{aligned}
$$

$-\Delta$ の固有値は波数の二乗の和である。

$$
\begin{aligned}
\lambda_{m1}\big|_{m=3} &= \frac{13 \pi^{2}}{4}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
\phi_{m1}\big|_{x=\frac{1}{3}} &= \sqrt{2} \sin{\left(\pi y \right)} \\
\frac{\partial^1 (\phi_{m1})}{\partial x^1} &= \frac{3 \sqrt{2} \pi \sin{\left(\pi y \right)} \cos{\left(\frac{3 \pi x}{2} \right)}}{2} \\
\left.\frac{\partial^1 (\phi_{m1})}{\partial x^1}\right|_{x=\frac{1}{3}} &= 0 \\
\frac{\partial^2 (\phi_{m1})}{\partial x^2} &= - \frac{9 \sqrt{2} \pi^{2} \sin{\left(\frac{3 \pi x}{2} \right)} \sin{\left(\pi y \right)}}{4} \\
\left.\frac{\partial^2 (\phi_{m1})}{\partial x^2}\right|_{x=\frac{1}{3}} &= - \frac{9 \sqrt{2} \pi^{2} \sin{\left(\pi y \right)}}{4}
\end{aligned}
$$
:::

:::example[長方形で一方の量子数を変える（4）]
横の辺長 $2$、縦の辺長 $1$ とし、縦の第 1 モードに固定する。

$$
m=4
$$

横の Dirichlet 条件から正弦の波数が決まる。

$$
\begin{aligned}
k_x\big|_{m=4} &= 2 \pi
\end{aligned}
$$

固有関数は二つの正弦の積であり、規格化に $\sqrt2$ を掛ける。

$$
\begin{aligned}
\phi_{m1}\big|_{m=4} &= \sqrt{2} \sin{\left(2 \pi x \right)} \sin{\left(\pi y \right)}
\end{aligned}
$$

$-\Delta$ の固有値は波数の二乗の和である。

$$
\begin{aligned}
\lambda_{m1}\big|_{m=4} &= 5 \pi^{2}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
\phi_{m1}\big|_{x=\frac{1}{3}} &= \frac{\sqrt{6} \sin{\left(\pi y \right)}}{2} \\
\frac{\partial^1 (\phi_{m1})}{\partial x^1} &= 2 \sqrt{2} \pi \sin{\left(\pi y \right)} \cos{\left(2 \pi x \right)} \\
\left.\frac{\partial^1 (\phi_{m1})}{\partial x^1}\right|_{x=\frac{1}{3}} &= - \sqrt{2} \pi \sin{\left(\pi y \right)} \\
\frac{\partial^2 (\phi_{m1})}{\partial x^2} &= - 4 \sqrt{2} \pi^{2} \sin{\left(2 \pi x \right)} \sin{\left(\pi y \right)} \\
\left.\frac{\partial^2 (\phi_{m1})}{\partial x^2}\right|_{x=\frac{1}{3}} &= - 2 \sqrt{6} \pi^{2} \sin{\left(\pi y \right)}
\end{aligned}
$$
:::

### 分離した境界データを係数にする

長方形 $0<x<L,0<y<H$ の三辺を零に固定し、上辺だけ $g(x)$ とする。下辺を零にするには各正弦モードの垂直因子を $\sinh(k_n y)$ に取る。上辺の値を合わせるには $\sinh(k_nH)$ で割る。したがって内部では高周波が $e^{-k_n(H-y)}$ の率で弱まる。$g$ の級数が境界で点ごとに収束するための条件と、任意のコンパクトな内部で微分できるための条件は異なる。後者には、係数が高々多項式的な大きさしか持たなければ、上の指数減衰を使える。前者を連続な全境界値として扱う場合は、上辺の端の値を側辺と一致させる。
:::example[長方形の上端で第 1 モードを指定する]
辺長をともに $1$ とし、Laplacian を作用させたときの二つの方向の寄与を分ける。

上端で $\sin(l\pi x)$ になり、下端で零になる鉛直方向の因子を選ぶ。

$$
\begin{aligned}
u &= \frac{\sin{\left(\pi x \right)} \sinh{\left(\pi y \right)}}{\sinh{\left(\pi \right)}}
\end{aligned}
$$

二つの二階微分は互いに反対符号であり、和が零である。

$$
\begin{aligned}
u_{xx} &= - \frac{\pi^{2} \sin{\left(\pi x \right)} \sinh{\left(\pi y \right)}}{\sinh{\left(\pi \right)}} \\
u_{yy} &= \frac{\pi^{2} \sin{\left(\pi x \right)} \sinh{\left(\pi y \right)}}{\sinh{\left(\pi \right)}}
\end{aligned}
$$

上端から半分の高さでは境界振幅に掛かる比が次の形になる。

$$
\begin{aligned}
R(1/2) &= \frac{1}{2 \cosh{\left(\frac{\pi}{2} \right)}}
\end{aligned}
$$
:::


:::example[長方形の上端で第 2 モードを指定する]
辺長をともに $1$ とし、Laplacian を作用させたときの二つの方向の寄与を分ける。

上端で $\sin(l\pi x)$ になり、下端で零になる鉛直方向の因子を選ぶ。

$$
\begin{aligned}
u &= \frac{\sin{\left(2 \pi x \right)} \sinh{\left(2 \pi y \right)}}{\sinh{\left(2 \pi \right)}}
\end{aligned}
$$

二つの二階微分は互いに反対符号であり、和が零である。

$$
\begin{aligned}
u_{xx} &= - \frac{4 \pi^{2} \sin{\left(2 \pi x \right)} \sinh{\left(2 \pi y \right)}}{\sinh{\left(2 \pi \right)}} \\
u_{yy} &= \frac{4 \pi^{2} \sin{\left(2 \pi x \right)} \sinh{\left(2 \pi y \right)}}{\sinh{\left(2 \pi \right)}}
\end{aligned}
$$

上端から半分の高さでは境界振幅に掛かる比が次の形になる。

$$
\begin{aligned}
R(1/2) &= \frac{1}{2 \cosh{\left(\pi \right)}}
\end{aligned}
$$
:::


:::example[長方形の上端で第 3 モードを指定する]
辺長をともに $1$ とし、Laplacian を作用させたときの二つの方向の寄与を分ける。

上端で $\sin(l\pi x)$ になり、下端で零になる鉛直方向の因子を選ぶ。

$$
\begin{aligned}
u &= \frac{\sin{\left(3 \pi x \right)} \sinh{\left(3 \pi y \right)}}{\sinh{\left(3 \pi \right)}}
\end{aligned}
$$

二つの二階微分は互いに反対符号であり、和が零である。

$$
\begin{aligned}
u_{xx} &= - \frac{9 \pi^{2} \sin{\left(3 \pi x \right)} \sinh{\left(3 \pi y \right)}}{\sinh{\left(3 \pi \right)}} \\
u_{yy} &= \frac{9 \pi^{2} \sin{\left(3 \pi x \right)} \sinh{\left(3 \pi y \right)}}{\sinh{\left(3 \pi \right)}}
\end{aligned}
$$

上端から半分の高さでは境界振幅に掛かる比が次の形になる。

$$
\begin{aligned}
R(1/2) &= \frac{1}{2 \cosh{\left(\frac{3 \pi}{2} \right)}}
\end{aligned}
$$
:::


:::example[長方形の上端で第 4 モードを指定する]
辺長をともに $1$ とし、Laplacian を作用させたときの二つの方向の寄与を分ける。

上端で $\sin(l\pi x)$ になり、下端で零になる鉛直方向の因子を選ぶ。

$$
\begin{aligned}
u &= \frac{\sin{\left(4 \pi x \right)} \sinh{\left(4 \pi y \right)}}{\sinh{\left(4 \pi \right)}}
\end{aligned}
$$

二つの二階微分は互いに反対符号であり、和が零である。

$$
\begin{aligned}
u_{xx} &= - \frac{16 \pi^{2} \sin{\left(4 \pi x \right)} \sinh{\left(4 \pi y \right)}}{\sinh{\left(4 \pi \right)}} \\
u_{yy} &= \frac{16 \pi^{2} \sin{\left(4 \pi x \right)} \sinh{\left(4 \pi y \right)}}{\sinh{\left(4 \pi \right)}}
\end{aligned}
$$

上端から半分の高さでは境界振幅に掛かる比が次の形になる。

$$
\begin{aligned}
R(1/2) &= \frac{1}{2 \cosh{\left(2 \pi \right)}}
\end{aligned}
$$
:::

## 固有関数展開による非斉次問題

### 非斉次項を同じ基底へ射影する

$-\Delta\phi_n=\lambda_n\phi_n$ の正規直交基底を持つ境界条件を選ぶ。$u_t-\kappa\Delta u=f$ に $u=\sum u_n(t)\phi_n$ を入れ、内積を取ると $u_n'+\kappa\lambda_n u_n=f_n$ である。積分因子で
$$
u_n(t)=e^{-\kappa\lambda_n t}u_n(0)+\int_0^t e^{-\kappa\lambda_n(t-s)}f_n(s)\dd s
\label{eq-mp-heat-modes}
$$
を得る。有限部分和では右辺も同じ範囲に射影した方程式が成り立つ。$u(0)\in L^2$、$f\in L^1([0,T];L^2)$、$\lambda_n\ge0$ とする。先頭 $N$ モードへの射影を $P_N$ とすると、残りのモードの $L^2$ ノルムは $\|(I-P_N)u(0)\|+\int_0^T\|(I-P_N)f(s)\|\dd s$ 以下である。完全性と優収束定理からこの上界は零に収束する。従って部分和は時間について一様に $L^2$ で収束し、コンパクトな台を持つ試験関数との積分で方程式を極限に移せる。古典解としての項別微分には追加の滑らかさを確認する。

波の場合は $u_n''+c^2\lambda_nu_n=f_n$ となる。$\omega_n=c\sqrt{\lambda_n}>0$ とし、初期変位と速度の係数を $a_n,b_n$ とすれば、
$$
u_n(t)=a_n\cos(\omega_nt)+b_n\frac{\sin(\omega_nt)}{\omega_n}
+\int_0^t\frac{\sin(\omega_n(t-s))}{\omega_n} f_n(s)\,\dd s
$$
である。積分項は一回微分すると上端が零、二回微分すると上端が $f_n(t)$ になる。零モードでは $\sin(\omega_nt)/\omega_n$ と積分核を、それぞれ $t$ と $t-s$ に置き換える。例えば $f_n(t)=F\sin(\omega_nt)$、初期値がともに零なら
$$
u_n(t)=\frac{F}{2\omega_n^2}\sin(\omega_nt)
-\frac{F}{2\omega_n}t\cos(\omega_nt)
$$
である。右辺を二回微分すると、$t\cos(\omega_nt)$ の微分から $F\sin(\omega_nt)$ が残り、二項の初期速度は相殺する。

:::example[固定温度と熱源]
$0<x<1$ で両端零、$u(0,x)=0$、$f=\sin\pi x$ なら一モードだけが励起され、$u=(1-e^{-\kappa\pi^2t})\sin\pi x/(\kappa\pi^2)$ となる。端の値が時間依存の場合は、その値に一致する滑らかな関数 $v(t,x)$ を先に引き、$f-v_t+\kappa\Delta v$ を展開する。
:::

### 持続する源と共鳴を区別する

熱方程式の各成分は $a_n^\prime+\kappa\lambda_n a_n=f_n$ である。一定の源なら、積分で $(1-e^{-\kappa\lambda_n t})/(\kappa\lambda_n)$ が出る。零固有値では極限を取れば $t$ であり、常に定常解へ近づくとは限らない。波動方程式では固有角周波数と同じ正弦入力で、畳み込みに線形の時間因子が出る。これらは時間の微分の階数と固有値が異なる例である。

時間を含まない非斉次問題 $(-\Delta-\mu)u=f$ では、固有関数の係数が $(\lambda_n-\mu)u_n=f_n$ を満たす。$-\Delta$ の自己共役実現が離散スペクトルと完全な固有基底を持つ場合、$\mu$ が固有値と異なれば $u_n=f_n/(\lambda_n-\mu)$ が一意解を与える。$f\in L^2$ に対して $\lambda_n/(\lambda_n-\mu)$ は有界なので $\sum_n|\lambda_nu_n|^2<\infty$ となり、解は作用素の定義域に属する。$\mu$ が固有値なら、該当する全ての $n$ で $f_n=0$ が必要であり、それらの $u_n$ を任意に取れる。

:::example[長方形の定常問題で共鳴する源]
単位正方形の四辺で $u=0$ とする。$f=\sin(\pi x)\sin(2\pi y)$ の固有値は $5\pi^2$ なので、$(-\Delta-5\pi^2)u=f$ は解を持たない。左辺を $f$ と内積すると自己共役性で零だが、右辺と $f$ の内積は $1/4$ だからである。代わりに右辺を $\sin(\pi x)\sin(\pi y)$ とすれば、$-1/(3\pi^2)$ 倍したものが特解である。$5\pi^2$ の固有空間は $(n,m)=(1,2),(2,1)$ に対応する二次元なので、斉次解としてその任意の組合せを加えられる。
:::

:::example[モードごとの定常値と過渡項（1）]
初期値零で、独立な減衰率 $m^2$ のモードに一定の強さ $m$ の熱源を加える。

$$
m=1
$$

積分因子で強制応答を求める。

$$
\begin{aligned}
a_m(t)\big|_{m=1} &= 1 - e^{- t}
\end{aligned}
$$

定常値との差は最初の符号を保って指数で減衰する。

$$
\begin{aligned}
a_m(t)-a_m(\infty)\big|_{m=1} &= - e^{- t}
\end{aligned}
$$

微分は時間の初期点で源の強さと一致する。

$$
\begin{aligned}
a_m^\prime(t)\big|_{m=1} &= e^{- t}
\end{aligned}
$$

観測点 $t=1$ での値と導関数を計算する。

$$
\begin{aligned}
a_m(t)\big|_{t=1} &= 1 - e^{-1} \\
\frac{\partial^1 (a_m(t))}{\partial t^1} &= e^{- t} \\
\left.\frac{\partial^1 (a_m(t))}{\partial t^1}\right|_{t=1} &= e^{-1} \\
\frac{\partial^2 (a_m(t))}{\partial t^2} &= - e^{- t} \\
\left.\frac{\partial^2 (a_m(t))}{\partial t^2}\right|_{t=1} &= - \frac{1}{e}
\end{aligned}
$$

観測点 $t=1$ での値と導関数を計算する。

$$
\begin{aligned}
a_m(t)-a_m(\infty)\big|_{t=1} &= - \frac{1}{e} \\
\frac{\partial^1 (a_m(t)-a_m(\infty))}{\partial t^1} &= e^{- t} \\
\left.\frac{\partial^1 (a_m(t)-a_m(\infty))}{\partial t^1}\right|_{t=1} &= e^{-1} \\
\frac{\partial^2 (a_m(t)-a_m(\infty))}{\partial t^2} &= - e^{- t} \\
\left.\frac{\partial^2 (a_m(t)-a_m(\infty))}{\partial t^2}\right|_{t=1} &= - \frac{1}{e}
\end{aligned}
$$

観測点 $t=1$ での値と導関数を計算する。

$$
\begin{aligned}
a_m^\prime(t)\big|_{t=1} &= e^{-1} \\
\frac{\partial^1 (a_m^\prime(t))}{\partial t^1} &= - e^{- t} \\
\left.\frac{\partial^1 (a_m^\prime(t))}{\partial t^1}\right|_{t=1} &= - \frac{1}{e} \\
\frac{\partial^2 (a_m^\prime(t))}{\partial t^2} &= e^{- t} \\
\left.\frac{\partial^2 (a_m^\prime(t))}{\partial t^2}\right|_{t=1} &= e^{-1}
\end{aligned}
$$
:::

:::example[モードごとの定常値と過渡項（2）]
初期値零で、独立な減衰率 $m^2$ のモードに一定の強さ $m$ の熱源を加える。

$$
m=2
$$

積分因子で強制応答を求める。

$$
\begin{aligned}
a_m(t)\big|_{m=2} &= \frac{1}{2} - \frac{e^{- 4 t}}{2}
\end{aligned}
$$

定常値との差は最初の符号を保って指数で減衰する。

$$
\begin{aligned}
a_m(t)-a_m(\infty)\big|_{m=2} &= - \frac{e^{- 4 t}}{2}
\end{aligned}
$$

微分は時間の初期点で源の強さと一致する。

$$
\begin{aligned}
a_m^\prime(t)\big|_{m=2} &= 2 e^{- 4 t}
\end{aligned}
$$

観測点 $t=1$ での値と導関数を計算する。

$$
\begin{aligned}
a_m(t)\big|_{t=1} &= - \frac{1 - e^{4}}{2 e^{4}} \\
\frac{\partial^1 (a_m(t))}{\partial t^1} &= 2 e^{- 4 t} \\
\left.\frac{\partial^1 (a_m(t))}{\partial t^1}\right|_{t=1} &= \frac{2}{e^{4}} \\
\frac{\partial^2 (a_m(t))}{\partial t^2} &= - 8 e^{- 4 t} \\
\left.\frac{\partial^2 (a_m(t))}{\partial t^2}\right|_{t=1} &= - \frac{8}{e^{4}}
\end{aligned}
$$

観測点 $t=1$ での値と導関数を計算する。

$$
\begin{aligned}
a_m(t)-a_m(\infty)\big|_{t=1} &= - \frac{1}{2 e^{4}} \\
\frac{\partial^1 (a_m(t)-a_m(\infty))}{\partial t^1} &= 2 e^{- 4 t} \\
\left.\frac{\partial^1 (a_m(t)-a_m(\infty))}{\partial t^1}\right|_{t=1} &= \frac{2}{e^{4}} \\
\frac{\partial^2 (a_m(t)-a_m(\infty))}{\partial t^2} &= - 8 e^{- 4 t} \\
\left.\frac{\partial^2 (a_m(t)-a_m(\infty))}{\partial t^2}\right|_{t=1} &= - \frac{8}{e^{4}}
\end{aligned}
$$

観測点 $t=1$ での値と導関数を計算する。

$$
\begin{aligned}
a_m^\prime(t)\big|_{t=1} &= \frac{2}{e^{4}} \\
\frac{\partial^1 (a_m^\prime(t))}{\partial t^1} &= - 8 e^{- 4 t} \\
\left.\frac{\partial^1 (a_m^\prime(t))}{\partial t^1}\right|_{t=1} &= - \frac{8}{e^{4}} \\
\frac{\partial^2 (a_m^\prime(t))}{\partial t^2} &= 32 e^{- 4 t} \\
\left.\frac{\partial^2 (a_m^\prime(t))}{\partial t^2}\right|_{t=1} &= \frac{32}{e^{4}}
\end{aligned}
$$
:::

:::example[モードごとの定常値と過渡項（3）]
初期値零で、独立な減衰率 $m^2$ のモードに一定の強さ $m$ の熱源を加える。

$$
m=3
$$

積分因子で強制応答を求める。

$$
\begin{aligned}
a_m(t)\big|_{m=3} &= \frac{1}{3} - \frac{e^{- 9 t}}{3}
\end{aligned}
$$

定常値との差は最初の符号を保って指数で減衰する。

$$
\begin{aligned}
a_m(t)-a_m(\infty)\big|_{m=3} &= - \frac{e^{- 9 t}}{3}
\end{aligned}
$$

微分は時間の初期点で源の強さと一致する。

$$
\begin{aligned}
a_m^\prime(t)\big|_{m=3} &= 3 e^{- 9 t}
\end{aligned}
$$

観測点 $t=1$ での値と導関数を計算する。

$$
\begin{aligned}
a_m(t)\big|_{t=1} &= - \frac{1 - e^{9}}{3 e^{9}} \\
\frac{\partial^1 (a_m(t))}{\partial t^1} &= 3 e^{- 9 t} \\
\left.\frac{\partial^1 (a_m(t))}{\partial t^1}\right|_{t=1} &= \frac{3}{e^{9}} \\
\frac{\partial^2 (a_m(t))}{\partial t^2} &= - 27 e^{- 9 t} \\
\left.\frac{\partial^2 (a_m(t))}{\partial t^2}\right|_{t=1} &= - \frac{27}{e^{9}}
\end{aligned}
$$

観測点 $t=1$ での値と導関数を計算する。

$$
\begin{aligned}
a_m(t)-a_m(\infty)\big|_{t=1} &= - \frac{1}{3 e^{9}} \\
\frac{\partial^1 (a_m(t)-a_m(\infty))}{\partial t^1} &= 3 e^{- 9 t} \\
\left.\frac{\partial^1 (a_m(t)-a_m(\infty))}{\partial t^1}\right|_{t=1} &= \frac{3}{e^{9}} \\
\frac{\partial^2 (a_m(t)-a_m(\infty))}{\partial t^2} &= - 27 e^{- 9 t} \\
\left.\frac{\partial^2 (a_m(t)-a_m(\infty))}{\partial t^2}\right|_{t=1} &= - \frac{27}{e^{9}}
\end{aligned}
$$

観測点 $t=1$ での値と導関数を計算する。

$$
\begin{aligned}
a_m^\prime(t)\big|_{t=1} &= \frac{3}{e^{9}} \\
\frac{\partial^1 (a_m^\prime(t))}{\partial t^1} &= - 27 e^{- 9 t} \\
\left.\frac{\partial^1 (a_m^\prime(t))}{\partial t^1}\right|_{t=1} &= - \frac{27}{e^{9}} \\
\frac{\partial^2 (a_m^\prime(t))}{\partial t^2} &= 243 e^{- 9 t} \\
\left.\frac{\partial^2 (a_m^\prime(t))}{\partial t^2}\right|_{t=1} &= \frac{243}{e^{9}}
\end{aligned}
$$
:::

:::example[モードごとの定常値と過渡項（4）]
初期値零で、独立な減衰率 $m^2$ のモードに一定の強さ $m$ の熱源を加える。

$$
m=4
$$

積分因子で強制応答を求める。

$$
\begin{aligned}
a_m(t)\big|_{m=4} &= \frac{1}{4} - \frac{e^{- 16 t}}{4}
\end{aligned}
$$

定常値との差は最初の符号を保って指数で減衰する。

$$
\begin{aligned}
a_m(t)-a_m(\infty)\big|_{m=4} &= - \frac{e^{- 16 t}}{4}
\end{aligned}
$$

微分は時間の初期点で源の強さと一致する。

$$
\begin{aligned}
a_m^\prime(t)\big|_{m=4} &= 4 e^{- 16 t}
\end{aligned}
$$

観測点 $t=1$ での値と導関数を計算する。

$$
\begin{aligned}
a_m(t)\big|_{t=1} &= - \frac{1 - e^{16}}{4 e^{16}} \\
\frac{\partial^1 (a_m(t))}{\partial t^1} &= 4 e^{- 16 t} \\
\left.\frac{\partial^1 (a_m(t))}{\partial t^1}\right|_{t=1} &= \frac{4}{e^{16}} \\
\frac{\partial^2 (a_m(t))}{\partial t^2} &= - 64 e^{- 16 t} \\
\left.\frac{\partial^2 (a_m(t))}{\partial t^2}\right|_{t=1} &= - \frac{64}{e^{16}}
\end{aligned}
$$

観測点 $t=1$ での値と導関数を計算する。

$$
\begin{aligned}
a_m(t)-a_m(\infty)\big|_{t=1} &= - \frac{1}{4 e^{16}} \\
\frac{\partial^1 (a_m(t)-a_m(\infty))}{\partial t^1} &= 4 e^{- 16 t} \\
\left.\frac{\partial^1 (a_m(t)-a_m(\infty))}{\partial t^1}\right|_{t=1} &= \frac{4}{e^{16}} \\
\frac{\partial^2 (a_m(t)-a_m(\infty))}{\partial t^2} &= - 64 e^{- 16 t} \\
\left.\frac{\partial^2 (a_m(t)-a_m(\infty))}{\partial t^2}\right|_{t=1} &= - \frac{64}{e^{16}}
\end{aligned}
$$

観測点 $t=1$ での値と導関数を計算する。

$$
\begin{aligned}
a_m^\prime(t)\big|_{t=1} &= \frac{4}{e^{16}} \\
\frac{\partial^1 (a_m^\prime(t))}{\partial t^1} &= - 64 e^{- 16 t} \\
\left.\frac{\partial^1 (a_m^\prime(t))}{\partial t^1}\right|_{t=1} &= - \frac{64}{e^{16}} \\
\frac{\partial^2 (a_m^\prime(t))}{\partial t^2} &= 1024 e^{- 16 t} \\
\left.\frac{\partial^2 (a_m^\prime(t))}{\partial t^2}\right|_{t=1} &= \frac{1024}{e^{16}}
\end{aligned}
$$
:::

## 積分変換による解法

### 無限領域を周波数で解く

実線で $u_t=\kappa u_{xx}$、$\kappa>0$ とし、初期値 $f$ を Schwartz 関数に取る。Fourier 変換は $\widehat f(k)=\int_\R f(x)e^{-ikx}\dd x$ とする。空間変換で $\partial_t\widehat u=-\kappa k^2\widehat u$ となり、逆変換すると $u=G_t*f$、$G_t=(4\pi\kappa t)^{-1/2}e^{-x^2/(4\kappa t)}$ である。核の積分は $1$、分散は $2\kappa t$ である。初期値への収束は近似単位の議論で得られる。例えば $f=e^{-a x^2}$、$a>0$ なら平方完成で $u=(1+4a\kappa t)^{-1/2}e^{-a x^2/(1+4a\kappa t)}$ となる。

実線の波動方程式では変換後の初期値解は $\cos(ckt)\widehat f+\sin(ckt)\widehat g/(ck)$ である。余弦は二つの平行移動の平均、後者は $(2c)^{-1}\mathbf1_{[-ct,ct]}$ の変換である。従って $u(t,x)=[f(x-ct)+f(x+ct)]/2+(2c)^{-1}\int_{x-ct}^{x+ct}g(y)\dd y$ になる。速度データの影響はこの区間に限られる。

:::example[周期境界での熱の平均]
長さ $L$ の周期区間で $f(x)=T_0+A\cos(2\pi x/L)$ を初期温度とする。各 Fourier モードは独立なので $u(t,x)=T_0+A e^{-4\pi^2\kappa t/L^2}\cos(2\pi x/L)$ となる。第 0 モードは時間変化しない。これは $u_t=\kappa u_{xx}$ を一周積分すると境界の熱流が周期性で相殺されることと一致する。例えば $L=2\pi\,\mathrm m$、$\kappa=1\,\mathrm{m^2\,s^{-1}}$ なら振幅は $1\,\mathrm s$ 後に $e^{-1}$ 倍になる。
:::

### 初期 Gaussian と熱の時刻

熱方程式の Fourier 変換は $\hat u(t,k)=e^{-\kappa k^2 t}\hat u(0,k)$ である。初期分布が分散 $\sigma_0^2$ の Gaussian なら、指数の二次係数が加算されて分散は $\sigma_0^2+2\kappa t$ になる。時間が増すと幅が広がり、質量を保つために中心の高さが下がる。単位は $\kappa$ が長さの二乗毎時、$k$ が長さの逆数であり、指数は無次元になる。

### 片側 Laplace 変換で境界からの信号を追う

半直線 $x>0$ で熱方程式を解き、初期温度を零、端の温度を $g(t)$ とする。$g$ は区分的に滑らかで指数位数を持つとし、変換が存在する右半平面で $U(x,s)=\int_0^\infty e^{-st}u(x,t)\dd t$ と置く。時間微分を部分積分すると $\kappa U_{xx}=sU$ となる。$\sqrt{s}$ は $\Re\sqrt{s}>0$ の枝を取り、$x\to\infty$ で有界な解を選ぶと $U(x,s)=\widetilde g(s)e^{-x\sqrt{s/\kappa}}$ である。

:::example[境界温度を一定値に切り替える]
$g(t)=1$、$t>0$ では $\widetilde g(s)=1/s$ である。逆変換の候補を
$$
u(t,x)=\operatorname{erfc}\!\left(\frac{x}{2\sqrt{\kappa t}}\right),
\qquad \operatorname{erfc}(a)=\frac{2}{\sqrt\pi}\int_a^\infty e^{-v^2}\,\dd v
$$
とする。$u_x=-(\pi\kappa t)^{-1/2}e^{-x^2/(4\kappa t)}$ を微分すれば $u_t=\kappa u_{xx}$ を確認できる。$x>0$ で $t\downarrow0$ に $u\to0$、$t>0$ で $u(t,0)=1$ である。この有界解の変換は上の常微分方程式を満たすので $e^{-x\sqrt{s/\kappa}}/s$ に一致する。原点 $(t,x)=(0,0)$ では初期値と境界値の極限が一致しない。
:::

:::example[上半平面の Poisson 核]
$y>0$ で $\Delta u=0$、$u(x,0)=g(x)$ とし、$g$ を Schwartz 関数とする。$x$ の Fourier 変換を取ると $\widehat u_{yy}-k^2\widehat u=0$ である。$y$ が大きい側で増大する指数を除くと $\widehat u(k,y)=e^{-|k|y}\widehat g(k)$ を得る。逆変換を $k\ge0$ と $k<0$ に分けると
$$
P_y(x)=\frac1\pi\int_0^\infty e^{-ky}\cos(kx)\,\dd k
=\frac1\pi\Re\frac1{y-ix}=\frac{y}{\pi(x^2+y^2)}
$$
となる。$\int_\R P_y=1$ であり、$|x|>\varepsilon$ 上の積分が $2y/(\pi\varepsilon)$ 以下なので $P_y*g$ は $y\downarrow0$ で $g$ に戻る。
:::

:::example[熱核で幅を増やした分布（1）]
初期分散を $m$、拡散係数を単位の $1$ とする。

$$
m=1
$$

変換に熱の乗数を掛ける。

$$
\begin{aligned}
\hat u(t,k)\big|_{m=1} &= e^{- \frac{k^{2} \left(2 t + 1\right)}{2}}
\end{aligned}
$$

反転して規格化した熱の分布を得る。

$$
\begin{aligned}
u(t,x)\big|_{m=1} &= \frac{\sqrt{2} e^{- \frac{x^{2}}{2 \left(2 t + 1\right)}}}{2 \sqrt{\pi} \sqrt{2 t + 1}}
\end{aligned}
$$

中心の高さを求める。

$$
\begin{aligned}
u(t,0)\big|_{m=1} &= \frac{\sqrt{2}}{2 \sqrt{\pi} \sqrt{2 t + 1}}
\end{aligned}
$$

観測点 $t=2$ での値と導関数を計算する。

$$
\begin{aligned}
\hat u(t,k)\big|_{t=2} &= e^{- \frac{5 k^{2}}{2}} \\
\frac{\partial^1 (\hat u(t,k))}{\partial t^1} &= - k^{2} e^{- \frac{k^{2} \left(2 t + 1\right)}{2}} \\
\left.\frac{\partial^1 (\hat u(t,k))}{\partial t^1}\right|_{t=2} &= - k^{2} e^{- \frac{5 k^{2}}{2}}
\end{aligned}
$$

観測点 $t=2$ での値と導関数を計算する。

$$
\begin{aligned}
u(t,x)\big|_{t=2} &= \frac{\sqrt{10} e^{- \frac{x^{2}}{10}}}{10 \sqrt{\pi}} \\
\frac{\partial^1 (u(t,x))}{\partial t^1} &= \frac{\sqrt{2} \left(- 2 t + x^{2} - 1\right) e^{- \frac{x^{2}}{2 \left(2 t + 1\right)}}}{2 \sqrt{\pi} \left(2 t + 1\right)^{\frac{5}{2}}} \\
\left.\frac{\partial^1 (u(t,x))}{\partial t^1}\right|_{t=2} &= \frac{\sqrt{10} \left(x^{2} - 5\right) e^{- \frac{x^{2}}{10}}}{250 \sqrt{\pi}}
\end{aligned}
$$

観測点 $t=2$ での値と導関数を計算する。

$$
\begin{aligned}
u(t,0)\big|_{t=2} &= \frac{\sqrt{10}}{10 \sqrt{\pi}} \\
\frac{\partial^1 (u(t,0))}{\partial t^1} &= - \frac{\sqrt{2}}{2 \sqrt{\pi} \left(2 t + 1\right)^{\frac{3}{2}}} \\
\left.\frac{\partial^1 (u(t,0))}{\partial t^1}\right|_{t=2} &= - \frac{\sqrt{10}}{50 \sqrt{\pi}}
\end{aligned}
$$
:::

:::example[熱核で幅を増やした分布（2）]
初期分散を $m$、拡散係数を単位の $1$ とする。

$$
m=2
$$

変換に熱の乗数を掛ける。

$$
\begin{aligned}
\hat u(t,k)\big|_{m=2} &= e^{- k^{2} \left(t + 1\right)}
\end{aligned}
$$

反転して規格化した熱の分布を得る。

$$
\begin{aligned}
u(t,x)\big|_{m=2} &= \frac{e^{- \frac{x^{2}}{4 \left(t + 1\right)}}}{2 \sqrt{\pi} \sqrt{t + 1}}
\end{aligned}
$$

中心の高さを求める。

$$
\begin{aligned}
u(t,0)\big|_{m=2} &= \frac{1}{2 \sqrt{\pi} \sqrt{t + 1}}
\end{aligned}
$$

観測点 $t=2$ での値と導関数を計算する。

$$
\begin{aligned}
\hat u(t,k)\big|_{t=2} &= e^{- 3 k^{2}} \\
\frac{\partial^1 (\hat u(t,k))}{\partial t^1} &= - k^{2} e^{- k^{2} \left(t + 1\right)} \\
\left.\frac{\partial^1 (\hat u(t,k))}{\partial t^1}\right|_{t=2} &= - k^{2} e^{- 3 k^{2}}
\end{aligned}
$$

観測点 $t=2$ での値と導関数を計算する。

$$
\begin{aligned}
u(t,x)\big|_{t=2} &= \frac{\sqrt{3} e^{- \frac{x^{2}}{12}}}{6 \sqrt{\pi}} \\
\frac{\partial^1 (u(t,x))}{\partial t^1} &= \frac{\left(- 2 t + x^{2} - 2\right) e^{- \frac{x^{2}}{4 \left(t + 1\right)}}}{8 \sqrt{\pi} \left(t + 1\right)^{\frac{5}{2}}} \\
\left.\frac{\partial^1 (u(t,x))}{\partial t^1}\right|_{t=2} &= \frac{\sqrt{3} \left(x^{2} - 6\right) e^{- \frac{x^{2}}{12}}}{216 \sqrt{\pi}}
\end{aligned}
$$

観測点 $t=2$ での値と導関数を計算する。

$$
\begin{aligned}
u(t,0)\big|_{t=2} &= \frac{\sqrt{3}}{6 \sqrt{\pi}} \\
\frac{\partial^1 (u(t,0))}{\partial t^1} &= - \frac{1}{4 \sqrt{\pi} \left(t + 1\right)^{\frac{3}{2}}} \\
\left.\frac{\partial^1 (u(t,0))}{\partial t^1}\right|_{t=2} &= - \frac{\sqrt{3}}{36 \sqrt{\pi}}
\end{aligned}
$$
:::

:::example[熱核で幅を増やした分布（3）]
初期分散を $m$、拡散係数を単位の $1$ とする。

$$
m=3
$$

変換に熱の乗数を掛ける。

$$
\begin{aligned}
\hat u(t,k)\big|_{m=3} &= e^{- \frac{k^{2} \left(2 t + 3\right)}{2}}
\end{aligned}
$$

反転して規格化した熱の分布を得る。

$$
\begin{aligned}
u(t,x)\big|_{m=3} &= \frac{\sqrt{2} e^{- \frac{x^{2}}{2 \left(2 t + 3\right)}}}{2 \sqrt{\pi} \sqrt{2 t + 3}}
\end{aligned}
$$

中心の高さを求める。

$$
\begin{aligned}
u(t,0)\big|_{m=3} &= \frac{\sqrt{2}}{2 \sqrt{\pi} \sqrt{2 t + 3}}
\end{aligned}
$$

観測点 $t=2$ での値と導関数を計算する。

$$
\begin{aligned}
\hat u(t,k)\big|_{t=2} &= e^{- \frac{7 k^{2}}{2}} \\
\frac{\partial^1 (\hat u(t,k))}{\partial t^1} &= - k^{2} e^{- \frac{k^{2} \left(2 t + 3\right)}{2}} \\
\left.\frac{\partial^1 (\hat u(t,k))}{\partial t^1}\right|_{t=2} &= - k^{2} e^{- \frac{7 k^{2}}{2}}
\end{aligned}
$$

観測点 $t=2$ での値と導関数を計算する。

$$
\begin{aligned}
u(t,x)\big|_{t=2} &= \frac{\sqrt{14} e^{- \frac{x^{2}}{14}}}{14 \sqrt{\pi}} \\
\frac{\partial^1 (u(t,x))}{\partial t^1} &= \frac{\sqrt{2} \left(- 2 t + x^{2} - 3\right) e^{- \frac{x^{2}}{2 \left(2 t + 3\right)}}}{2 \sqrt{\pi} \left(2 t + 3\right)^{\frac{5}{2}}} \\
\left.\frac{\partial^1 (u(t,x))}{\partial t^1}\right|_{t=2} &= \frac{\sqrt{14} \left(x^{2} - 7\right) e^{- \frac{x^{2}}{14}}}{686 \sqrt{\pi}}
\end{aligned}
$$

観測点 $t=2$ での値と導関数を計算する。

$$
\begin{aligned}
u(t,0)\big|_{t=2} &= \frac{\sqrt{14}}{14 \sqrt{\pi}} \\
\frac{\partial^1 (u(t,0))}{\partial t^1} &= - \frac{\sqrt{2}}{2 \sqrt{\pi} \left(2 t + 3\right)^{\frac{3}{2}}} \\
\left.\frac{\partial^1 (u(t,0))}{\partial t^1}\right|_{t=2} &= - \frac{\sqrt{14}}{98 \sqrt{\pi}}
\end{aligned}
$$
:::

:::example[熱核で幅を増やした分布（4）]
初期分散を $m$、拡散係数を単位の $1$ とする。

$$
m=4
$$

変換に熱の乗数を掛ける。

$$
\begin{aligned}
\hat u(t,k)\big|_{m=4} &= e^{- k^{2} \left(t + 2\right)}
\end{aligned}
$$

反転して規格化した熱の分布を得る。

$$
\begin{aligned}
u(t,x)\big|_{m=4} &= \frac{e^{- \frac{x^{2}}{4 \left(t + 2\right)}}}{2 \sqrt{\pi} \sqrt{t + 2}}
\end{aligned}
$$

中心の高さを求める。

$$
\begin{aligned}
u(t,0)\big|_{m=4} &= \frac{1}{2 \sqrt{\pi} \sqrt{t + 2}}
\end{aligned}
$$

観測点 $t=2$ での値と導関数を計算する。

$$
\begin{aligned}
\hat u(t,k)\big|_{t=2} &= e^{- 4 k^{2}} \\
\frac{\partial^1 (\hat u(t,k))}{\partial t^1} &= - k^{2} e^{- k^{2} \left(t + 2\right)} \\
\left.\frac{\partial^1 (\hat u(t,k))}{\partial t^1}\right|_{t=2} &= - k^{2} e^{- 4 k^{2}}
\end{aligned}
$$

観測点 $t=2$ での値と導関数を計算する。

$$
\begin{aligned}
u(t,x)\big|_{t=2} &= \frac{e^{- \frac{x^{2}}{16}}}{4 \sqrt{\pi}} \\
\frac{\partial^1 (u(t,x))}{\partial t^1} &= \frac{\left(- 2 t + x^{2} - 4\right) e^{- \frac{x^{2}}{4 \left(t + 2\right)}}}{8 \sqrt{\pi} \left(t + 2\right)^{\frac{5}{2}}} \\
\left.\frac{\partial^1 (u(t,x))}{\partial t^1}\right|_{t=2} &= \frac{\left(x^{2} - 8\right) e^{- \frac{x^{2}}{16}}}{256 \sqrt{\pi}}
\end{aligned}
$$

観測点 $t=2$ での値と導関数を計算する。

$$
\begin{aligned}
u(t,0)\big|_{t=2} &= \frac{1}{4 \sqrt{\pi}} \\
\frac{\partial^1 (u(t,0))}{\partial t^1} &= - \frac{1}{4 \sqrt{\pi} \left(t + 2\right)^{\frac{3}{2}}} \\
\left.\frac{\partial^1 (u(t,0))}{\partial t^1}\right|_{t=2} &= - \frac{1}{32 \sqrt{\pi}}
\end{aligned}
$$
:::

## 特性曲線法と 1 階偏微分方程式

### 曲線を選んで一変数の微分にする

$a(t,x,u)u_x+u_t=b(t,x,u)$ を考え、$a,b,f$ は $C^1$ とする。:term[特性曲線法]{en="method of characteristics" yomi="とくせいきょくせんほう"} では $\dot X=a(t,X,U)$、$\dot U=b(t,X,U)$、$X(0)=\xi$、$U(0)=f(\xi)$ を解く。初期線上で $\partial_\xi X=1$ だから、逆関数定理で局所的に $(t,\xi)$ を $(t,x)$ に戻せる。連鎖律で $u(t,X(t,\xi))=U(t,\xi)$ の導関数は $u_t+a u_x$ なので、得た関数は方程式を満たす。逆に滑らかな解があれば、その曲線に沿って同じ常微分方程式を解くため局所的に一意である。一般の初期曲線は、初期曲線の接線と平面上の特性ベクトルが一次独立な点で同じように扱える。

:::example[Hamilton–Jacobi 方程式の特性]
$S_t+H(t,x,S_x)=0$ とし、滑らかな解の $p=S_x$ を考える。$x$ 微分した式から $p_t+H_p p_x=-H_x$ が従う。よって滑らかな特性は $\dot X=H_p$、$\dot p=-H_x$、$\dot S=pH_p-H$ に従う。$H=p^2/2$、$S(0,x)=a x^2/2$ なら $p=a\xi$、$X=(1+at)\xi$、$S=(1+at)a\xi^2/2$ である。従って $1+at\ne0$ のうち初期線から連続に到達できる区間で $S=a x^2/[2(1+at)]$ となる。両偏微分を代入すると $-a^2x^2/[2(1+at)^2]+a^2x^2/[2(1+at)^2]=0$ である。
:::

:::example[Burgers 方程式の特性の交差]
$b=0,a=u$ では $U=f(\xi)$、$X=\xi+t f(\xi)$ である。例えば $f=-\sin \xi$ なら $\partial_\xi X=1-t\cos\xi$ で、正時刻で最初に零になるのは $t=1$、$\xi=2\pi n$（$n\in\Z$）である。そこへ近づくと $u_x=f'/(1+tf')$ が発散する。以後に弱解を選ぶには保存形 $u_t+(u^2/2)_x=0$ に加えてエントロピー条件が必要である。条件の導入は [保存則](/math/pde/conservation-laws/) で扱う。
:::

### 特性座標で積分因子を使う

線形輸送 $u_t+a u_x+b u=f(t,x)$ では、特性 $x(\tau)=x-a(t-\tau)$ に沿って $U^\prime+bU=f$ となる。積分因子 $e^{b\tau}$ で解き、初期値は $x-at$ で評価する。一方準線形方程式では $a$ が $u$ に依存する。初期点からの写像の導関数が零になると、特性座標の逆が取れず古典解の勾配が発散し得る。

:::example[減衰を伴う輸送波（1）]
速度 $m$、減衰率 $2$、外力零、初期値 $\cos x$ を選ぶ。

$$
m=1
$$

点 $(x,t)$ から初期線へ特性を戻す。

$$
\begin{aligned}
x_0\big|_{m=1} &= - t + x
\end{aligned}
$$

初期値を指数で減衰させる。

$$
\begin{aligned}
u(t,x)\big|_{m=1} &= e^{- 2 t} \cos{\left(t - x \right)}
\end{aligned}
$$

時間微分と移流項に含まれる位相微分は逆符号である。

$$
\begin{aligned}
u_t+m u_x\big|_{m=1} &= - 2 e^{- 2 t} \cos{\left(t - x \right)}
\end{aligned}
$$

観測点 $t=1$ での値と導関数を計算する。

$$
\begin{aligned}
x_0\big|_{t=1} &= x - 1 \\
\frac{\partial^1 (x_0)}{\partial t^1} &= -1 \\
\left.\frac{\partial^1 (x_0)}{\partial t^1}\right|_{t=1} &= -1 \\
\frac{\partial^2 (x_0)}{\partial t^2} &= 0 \\
\left.\frac{\partial^2 (x_0)}{\partial t^2}\right|_{t=1} &= 0
\end{aligned}
$$

観測点 $t=1$ での値と導関数を計算する。

$$
\begin{aligned}
u(t,x)\big|_{t=1} &= \frac{\cos{\left(x - 1 \right)}}{e^{2}} \\
\frac{\partial^1 (u(t,x))}{\partial t^1} &= - \left(\sin{\left(t - x \right)} + 2 \cos{\left(t - x \right)}\right) e^{- 2 t} \\
\left.\frac{\partial^1 (u(t,x))}{\partial t^1}\right|_{t=1} &= \frac{\sin{\left(x - 1 \right)} - 2 \cos{\left(x - 1 \right)}}{e^{2}} \\
\frac{\partial^2 (u(t,x))}{\partial t^2} &= \left(4 \sin{\left(t - x \right)} + 3 \cos{\left(t - x \right)}\right) e^{- 2 t} \\
\left.\frac{\partial^2 (u(t,x))}{\partial t^2}\right|_{t=1} &= \frac{- 4 \sin{\left(x - 1 \right)} + 3 \cos{\left(x - 1 \right)}}{e^{2}}
\end{aligned}
$$

観測点 $t=1$ での値と導関数を計算する。

$$
\begin{aligned}
u_t+m u_x\big|_{t=1} &= - \frac{2 \cos{\left(x - 1 \right)}}{e^{2}} \\
\frac{\partial^1 (u_t+m u_x)}{\partial t^1} &= 2 \left(\sin{\left(t - x \right)} + 2 \cos{\left(t - x \right)}\right) e^{- 2 t} \\
\left.\frac{\partial^1 (u_t+m u_x)}{\partial t^1}\right|_{t=1} &= \frac{2 \left(- \sin{\left(x - 1 \right)} + 2 \cos{\left(x - 1 \right)}\right)}{e^{2}} \\
\frac{\partial^2 (u_t+m u_x)}{\partial t^2} &= - \left(8 \sin{\left(t - x \right)} + 6 \cos{\left(t - x \right)}\right) e^{- 2 t} \\
\left.\frac{\partial^2 (u_t+m u_x)}{\partial t^2}\right|_{t=1} &= \frac{2 \left(4 \sin{\left(x - 1 \right)} - 3 \cos{\left(x - 1 \right)}\right)}{e^{2}}
\end{aligned}
$$
:::

:::example[減衰を伴う輸送波（2）]
速度 $m$、減衰率 $2$、外力零、初期値 $\cos x$ を選ぶ。

$$
m=2
$$

点 $(x,t)$ から初期線へ特性を戻す。

$$
\begin{aligned}
x_0\big|_{m=2} &= - 2 t + x
\end{aligned}
$$

初期値を指数で減衰させる。

$$
\begin{aligned}
u(t,x)\big|_{m=2} &= e^{- 2 t} \cos{\left(2 t - x \right)}
\end{aligned}
$$

時間微分と移流項に含まれる位相微分は逆符号である。

$$
\begin{aligned}
u_t+m u_x\big|_{m=2} &= - 2 e^{- 2 t} \cos{\left(2 t - x \right)}
\end{aligned}
$$

観測点 $t=1$ での値と導関数を計算する。

$$
\begin{aligned}
x_0\big|_{t=1} &= x - 2 \\
\frac{\partial^1 (x_0)}{\partial t^1} &= -2 \\
\left.\frac{\partial^1 (x_0)}{\partial t^1}\right|_{t=1} &= -2 \\
\frac{\partial^2 (x_0)}{\partial t^2} &= 0 \\
\left.\frac{\partial^2 (x_0)}{\partial t^2}\right|_{t=1} &= 0
\end{aligned}
$$

観測点 $t=1$ での値と導関数を計算する。

$$
\begin{aligned}
u(t,x)\big|_{t=1} &= \frac{\cos{\left(x - 2 \right)}}{e^{2}} \\
\frac{\partial^1 (u(t,x))}{\partial t^1} &= 2 \left(- \sin{\left(2 t - x \right)} - \cos{\left(2 t - x \right)}\right) e^{- 2 t} \\
\left.\frac{\partial^1 (u(t,x))}{\partial t^1}\right|_{t=1} &= - \frac{2 \sqrt{2} \cos{\left(x - 2 + \frac{\pi}{4} \right)}}{e^{2}} \\
\frac{\partial^2 (u(t,x))}{\partial t^2} &= 8 e^{- 2 t} \sin{\left(2 t - x \right)} \\
\left.\frac{\partial^2 (u(t,x))}{\partial t^2}\right|_{t=1} &= - \frac{8 \sin{\left(x - 2 \right)}}{e^{2}}
\end{aligned}
$$

観測点 $t=1$ での値と導関数を計算する。

$$
\begin{aligned}
u_t+m u_x\big|_{t=1} &= - \frac{2 \cos{\left(x - 2 \right)}}{e^{2}} \\
\frac{\partial^1 (u_t+m u_x)}{\partial t^1} &= 4 \left(\sin{\left(2 t - x \right)} + \cos{\left(2 t - x \right)}\right) e^{- 2 t} \\
\left.\frac{\partial^1 (u_t+m u_x)}{\partial t^1}\right|_{t=1} &= \frac{4 \sqrt{2} \cos{\left(x - 2 + \frac{\pi}{4} \right)}}{e^{2}} \\
\frac{\partial^2 (u_t+m u_x)}{\partial t^2} &= - 16 e^{- 2 t} \sin{\left(2 t - x \right)} \\
\left.\frac{\partial^2 (u_t+m u_x)}{\partial t^2}\right|_{t=1} &= \frac{16 \sin{\left(x - 2 \right)}}{e^{2}}
\end{aligned}
$$
:::

:::example[減衰を伴う輸送波（3）]
速度 $m$、減衰率 $2$、外力零、初期値 $\cos x$ を選ぶ。

$$
m=3
$$

点 $(x,t)$ から初期線へ特性を戻す。

$$
\begin{aligned}
x_0\big|_{m=3} &= - 3 t + x
\end{aligned}
$$

初期値を指数で減衰させる。

$$
\begin{aligned}
u(t,x)\big|_{m=3} &= e^{- 2 t} \cos{\left(3 t - x \right)}
\end{aligned}
$$

時間微分と移流項に含まれる位相微分は逆符号である。

$$
\begin{aligned}
u_t+m u_x\big|_{m=3} &= - 2 e^{- 2 t} \cos{\left(3 t - x \right)}
\end{aligned}
$$

観測点 $t=1$ での値と導関数を計算する。

$$
\begin{aligned}
x_0\big|_{t=1} &= x - 3 \\
\frac{\partial^1 (x_0)}{\partial t^1} &= -3 \\
\left.\frac{\partial^1 (x_0)}{\partial t^1}\right|_{t=1} &= -3 \\
\frac{\partial^2 (x_0)}{\partial t^2} &= 0 \\
\left.\frac{\partial^2 (x_0)}{\partial t^2}\right|_{t=1} &= 0
\end{aligned}
$$

観測点 $t=1$ での値と導関数を計算する。

$$
\begin{aligned}
u(t,x)\big|_{t=1} &= \frac{\cos{\left(x - 3 \right)}}{e^{2}} \\
\frac{\partial^1 (u(t,x))}{\partial t^1} &= - \left(3 \sin{\left(3 t - x \right)} + 2 \cos{\left(3 t - x \right)}\right) e^{- 2 t} \\
\left.\frac{\partial^1 (u(t,x))}{\partial t^1}\right|_{t=1} &= \frac{3 \sin{\left(x - 3 \right)} - 2 \cos{\left(x - 3 \right)}}{e^{2}} \\
\frac{\partial^2 (u(t,x))}{\partial t^2} &= \left(12 \sin{\left(3 t - x \right)} - 5 \cos{\left(3 t - x \right)}\right) e^{- 2 t} \\
\left.\frac{\partial^2 (u(t,x))}{\partial t^2}\right|_{t=1} &= - \frac{12 \sin{\left(x - 3 \right)} + 5 \cos{\left(x - 3 \right)}}{e^{2}}
\end{aligned}
$$

観測点 $t=1$ での値と導関数を計算する。

$$
\begin{aligned}
u_t+m u_x\big|_{t=1} &= - \frac{2 \cos{\left(x - 3 \right)}}{e^{2}} \\
\frac{\partial^1 (u_t+m u_x)}{\partial t^1} &= 2 \left(3 \sin{\left(3 t - x \right)} + 2 \cos{\left(3 t - x \right)}\right) e^{- 2 t} \\
\left.\frac{\partial^1 (u_t+m u_x)}{\partial t^1}\right|_{t=1} &= \frac{2 \left(- 3 \sin{\left(x - 3 \right)} + 2 \cos{\left(x - 3 \right)}\right)}{e^{2}} \\
\frac{\partial^2 (u_t+m u_x)}{\partial t^2} &= 2 \left(- 12 \sin{\left(3 t - x \right)} + 5 \cos{\left(3 t - x \right)}\right) e^{- 2 t} \\
\left.\frac{\partial^2 (u_t+m u_x)}{\partial t^2}\right|_{t=1} &= \frac{2 \left(12 \sin{\left(x - 3 \right)} + 5 \cos{\left(x - 3 \right)}\right)}{e^{2}}
\end{aligned}
$$
:::

:::example[減衰を伴う輸送波（4）]
速度 $m$、減衰率 $2$、外力零、初期値 $\cos x$ を選ぶ。

$$
m=4
$$

点 $(x,t)$ から初期線へ特性を戻す。

$$
\begin{aligned}
x_0\big|_{m=4} &= - 4 t + x
\end{aligned}
$$

初期値を指数で減衰させる。

$$
\begin{aligned}
u(t,x)\big|_{m=4} &= e^{- 2 t} \cos{\left(4 t - x \right)}
\end{aligned}
$$

時間微分と移流項に含まれる位相微分は逆符号である。

$$
\begin{aligned}
u_t+m u_x\big|_{m=4} &= - 2 e^{- 2 t} \cos{\left(4 t - x \right)}
\end{aligned}
$$

観測点 $t=1$ での値と導関数を計算する。

$$
\begin{aligned}
x_0\big|_{t=1} &= x - 4 \\
\frac{\partial^1 (x_0)}{\partial t^1} &= -4 \\
\left.\frac{\partial^1 (x_0)}{\partial t^1}\right|_{t=1} &= -4 \\
\frac{\partial^2 (x_0)}{\partial t^2} &= 0 \\
\left.\frac{\partial^2 (x_0)}{\partial t^2}\right|_{t=1} &= 0
\end{aligned}
$$

観測点 $t=1$ での値と導関数を計算する。

$$
\begin{aligned}
u(t,x)\big|_{t=1} &= \frac{\cos{\left(x - 4 \right)}}{e^{2}} \\
\frac{\partial^1 (u(t,x))}{\partial t^1} &= 2 \left(- 2 \sin{\left(4 t - x \right)} - \cos{\left(4 t - x \right)}\right) e^{- 2 t} \\
\left.\frac{\partial^1 (u(t,x))}{\partial t^1}\right|_{t=1} &= \frac{2 \left(2 \sin{\left(x - 4 \right)} - \cos{\left(x - 4 \right)}\right)}{e^{2}} \\
\frac{\partial^2 (u(t,x))}{\partial t^2} &= 4 \left(4 \sin{\left(4 t - x \right)} - 3 \cos{\left(4 t - x \right)}\right) e^{- 2 t} \\
\left.\frac{\partial^2 (u(t,x))}{\partial t^2}\right|_{t=1} &= - \frac{16 \sin{\left(x - 4 \right)} + 12 \cos{\left(x - 4 \right)}}{e^{2}}
\end{aligned}
$$

観測点 $t=1$ での値と導関数を計算する。

$$
\begin{aligned}
u_t+m u_x\big|_{t=1} &= - \frac{2 \cos{\left(x - 4 \right)}}{e^{2}} \\
\frac{\partial^1 (u_t+m u_x)}{\partial t^1} &= 4 \left(2 \sin{\left(4 t - x \right)} + \cos{\left(4 t - x \right)}\right) e^{- 2 t} \\
\left.\frac{\partial^1 (u_t+m u_x)}{\partial t^1}\right|_{t=1} &= \frac{4 \left(- 2 \sin{\left(x - 4 \right)} + \cos{\left(x - 4 \right)}\right)}{e^{2}} \\
\frac{\partial^2 (u_t+m u_x)}{\partial t^2} &= 8 \left(- 4 \sin{\left(4 t - x \right)} + 3 \cos{\left(4 t - x \right)}\right) e^{- 2 t} \\
\left.\frac{\partial^2 (u_t+m u_x)}{\partial t^2}\right|_{t=1} &= \frac{8 \left(4 \sin{\left(x - 4 \right)} + 3 \cos{\left(x - 4 \right)}\right)}{e^{2}}
\end{aligned}
$$
:::

## Helmholtz 方程式と散乱問題

### 外向きの波を選ぶ

時間因子 $e^{-i\omega t}$ を波動方程式に入れると $(\Delta+k^2)u=0$、$k=\omega/c$ である。外部領域では両方の球 Hankel 解が遠方で $1/r$ に減るため、単に零へ近づくという条件だけでは向きを選べない。三次元では $r(\partial_r-ik)u\to0$ を Sommerfeld の放射条件とする。球対称解 $Ae^{ikr}/r+Be^{-ikr}/r$ では第一項のこの量は $O(r^{-1})$、第二項は $-2ikBe^{-ikr}+O(r^{-1})$ なので $B=0$ が必要である。

散乱では入射解 $u_{\rm in}$ を先に指定し、散乱部分に放射条件を課す。球 $r=a$ で全場が零という境界を考え、入射部分の $\ell=0$ を $j_0(kr)$ とする。対応する散乱成分は $A h_0^{(1)}(kr)$、$h_0^{(1)}=j_0+iy_0=-ie^{ikr}/(kr)$ なので $A=-j_0(ka)/h_0^{(1)}(ka)=-i e^{-ika}\sin(ka)$ となる。線形性により各角モードを別に扱える。遠方の係数は角度別の振幅になり、$u_{\rm sc}\sim e^{ikr}f(\Omega)/r$ の $f$ は入射振幅を無次元とすれば長さの単位を持つ。

### 球散乱の位相差

球対称な散乱体では、外部で各 $\ell$ の動径解を $j_\ell\cos\delta_\ell-y_\ell\sin\delta_\ell$ と取れる。ここでは吸収を含まない実の境界条件を仮定する。遠方の $j_\ell,y_\ell$ の主項を入れると、相は $kr-\ell\pi/2$ から $\delta_\ell$ だけ進む。Dirichlet の球なら $r=a$ の値を零にして $\tan\delta_\ell=j_\ell(ka)/y_\ell(ka)$ を得る。比に極があるときも、条件 $j_\ell\cos\delta_\ell-y_\ell\sin\delta_\ell=0$ は使える。$\ell=0$ では $\tan\delta_0=-\tan(ka)$、従って $\delta_0=-ka$ を $\pi$ を法として取れる。前の散乱係数と同じことは、$j_\ell=(h^{(1)}_\ell+h^{(2)}_\ell)/2$、$y_\ell=(h^{(1)}_\ell-h^{(2)}_\ell)/(2i)$ で入射成分を合わせれば確認できる。入射の $j_\ell$ に足す外向き係数は $(e^{2i\delta_\ell}-1)/2$ であり、$\ell=0$ では $-i e^{-ika}\sin(ka)$ になる。

### 部分波の内向き成分と外向き成分

中心力散乱で $\ell=0$ の自由動径方程式は $v=ru$ とおくと $v^{\prime\prime}+k^2v=0$ になる。硬い球で $v(a)=0$ なら $v\propto\sin k(r-a)$ であり、位相差を $\delta_0=-ka$ と読める。位相差は $\pi$ を加えても同じ外向き・内向きの比を与える。一つの部分波の断面積は球の全面積そのものに固定されず、$4\pi\sin^2\delta_0/k^2$ である。

### 散乱振幅と全断面積

$k>0$ とし、回転対称な散乱体に平面波 $e^{ikz}$ を入射させる。:term[部分波]{en="partial wave" yomi="ぶぶんは"} は球面上の角運動量 $\ell$ ごとの成分である。角度の余弦を $\mu=\cos\theta$ とすれば、Legendre 直交性により入射波の展開係数は $(2\ell+1)i^\ell j_\ell(kr)$ である。実際 Rodrigues の公式を積分に入れ、端点で消える因子 $(\mu^2-1)^\ell$ を使って $\ell$ 回部分積分すると
$$
\frac12\int_{-1}^1 e^{iz\mu}P_\ell(\mu)\,\dd\mu
=\frac{(iz)^\ell}{2^{\ell+1}\ell!}\int_{-1}^1e^{iz\mu}(1-\mu^2)^\ell\,\dd\mu
=i^\ell j_\ell(z)
$$
となる。最後の等号は指数を展開し、奇数項を消し、偶数項を Beta 積分 $\int_{-1}^1\mu^{2n}(1-\mu^2)^\ell\dd\mu=B(n+1/2,\ell+1)$ に帰着すれば球 Bessel の級数と一致する。

$S_\ell=e^{2i\delta_\ell}$ と置けば、各入射 $j_\ell$ に加えるのは $(S_\ell-1)h_\ell^{(1)}/2$ である。$h_\ell^{(1)}(kr)\sim(-i)^{\ell+1}e^{ikr}/(kr)$ を用いて、遠方振幅の Legendre 係数を得る。遠方展開と角度積分が収束する散乱解では、直交性 $\int_{S^2}P_\ell P_m\,\dd\Omega=4\pi\delta_{\ell m}/(2\ell+1)$ から
$$
\begin{aligned}
f(\theta)&=\frac{1}{2ik}\sum_{\ell\ge0}(2\ell+1)(S_\ell-1)P_\ell(\cos\theta),\\
\frac{\dd\sigma}{\dd\Omega}&=|f(\theta)|^2,\qquad
\sigma=\frac{4\pi}{k^2}\sum_{\ell\ge0}(2\ell+1)\sin^2\delta_\ell
\end{aligned}
$$
となる。微分断面積は、遠方での散乱波の径方向の流束を、平面波の単位面積あたりの入射流束で割ったものである。例えば流束を $\Im(\bar u\grad u)$ と同じ規格化で計れば、両者に共通の $k$ が消える。さらに前方で級数が収束し、全ての位相差が実なら、$P_\ell(1)=1$ と $\Im[(e^{2i\delta_\ell}-1)/(2i)]=\sin^2\delta_\ell$ から、:term[光学定理]{en="optical theorem" yomi="こうがくていり"} $\sigma=(4\pi/k)\Im f(0)$ が従う。

:::example[低周波の Dirichlet 球]
$ka\ll1$ とする。$\ell=0$ だけの断面積は $4\pi a^2[\sin(ka)/(ka)]^2$ であり、$k\downarrow0$ で $4\pi a^2$ に近づく。固定した $\ell\ge1$ では $j_\ell(z)=O(z^\ell)$、$y_\ell(z)$ の主項は $-(2\ell-1)!!z^{-\ell-1}$ なので $\tan\delta_\ell=O((ka)^{2\ell+1})$ である。従って各高次部分波は固定した $\ell$ ごとに零次部分波より高次の抑制を受ける。
:::

:::example[硬い球の零次部分波（1）]
半径 $1$ で、$ka=m\pi/4$ の場合を比べる。

$$
m=1
$$

境界を零にする動径関数を選ぶ。

$$
\begin{aligned}
v(r)\big|_{m=1} &= \sin{\left(\pi \left(\frac{r}{4} - \frac{1}{4}\right) \right)}
\end{aligned}
$$

位相差を求める。

$$
\begin{aligned}
\delta_0\big|_{m=1} &= - \frac{\pi}{4}
\end{aligned}
$$

零次分だけの全断面積を計算する。

$$
\begin{aligned}
\sigma_0\big|_{m=1} &= \frac{32}{\pi}
\end{aligned}
$$

観測点 $r=2$ での値と導関数を計算する。

$$
\begin{aligned}
v(r)\big|_{r=2} &= \frac{\sqrt{2}}{2} \\
\frac{\partial^1 (v(r))}{\partial r^1} &= \frac{\pi \cos{\left(\pi \left(\frac{r}{4} - \frac{1}{4}\right) \right)}}{4} \\
\left.\frac{\partial^1 (v(r))}{\partial r^1}\right|_{r=2} &= \frac{\sqrt{2} \pi}{8}
\end{aligned}
$$
:::

:::example[硬い球の零次部分波（2）]
半径 $1$ で、$ka=m\pi/4$ の場合を比べる。

$$
m=2
$$

境界を零にする動径関数を選ぶ。

$$
\begin{aligned}
v(r)\big|_{m=2} &= \sin{\left(\pi \left(\frac{r}{2} - \frac{1}{2}\right) \right)}
\end{aligned}
$$

位相差を求める。

$$
\begin{aligned}
\delta_0\big|_{m=2} &= - \frac{\pi}{2}
\end{aligned}
$$

零次分だけの全断面積を計算する。

$$
\begin{aligned}
\sigma_0\big|_{m=2} &= \frac{16}{\pi}
\end{aligned}
$$

観測点 $r=2$ での値と導関数を計算する。

$$
\begin{aligned}
v(r)\big|_{r=2} &= 1 \\
\frac{\partial^1 (v(r))}{\partial r^1} &= \frac{\pi \sin{\left(\frac{\pi r}{2} \right)}}{2} \\
\left.\frac{\partial^1 (v(r))}{\partial r^1}\right|_{r=2} &= 0
\end{aligned}
$$
:::

:::example[硬い球の零次部分波（3）]
半径 $1$ で、$ka=m\pi/4$ の場合を比べる。

$$
m=3
$$

境界を零にする動径関数を選ぶ。

$$
\begin{aligned}
v(r)\big|_{m=3} &= \sin{\left(\pi \left(\frac{3 r}{4} - \frac{3}{4}\right) \right)}
\end{aligned}
$$

位相差を求める。

$$
\begin{aligned}
\delta_0\big|_{m=3} &= - \frac{3 \pi}{4}
\end{aligned}
$$

零次分だけの全断面積を計算する。

$$
\begin{aligned}
\sigma_0\big|_{m=3} &= \frac{32}{9 \pi}
\end{aligned}
$$

観測点 $r=2$ での値と導関数を計算する。

$$
\begin{aligned}
v(r)\big|_{r=2} &= \frac{\sqrt{2}}{2} \\
\frac{\partial^1 (v(r))}{\partial r^1} &= \frac{3 \pi \cos{\left(\pi \left(\frac{3 r}{4} - \frac{3}{4}\right) \right)}}{4} \\
\left.\frac{\partial^1 (v(r))}{\partial r^1}\right|_{r=2} &= - \frac{3 \sqrt{2} \pi}{8}
\end{aligned}
$$
:::

:::example[硬い球の零次部分波（4）]
半径 $1$ で、$ka=m\pi/4$ の場合を比べる。

$$
m=4
$$

境界を零にする動径関数を選ぶ。

$$
\begin{aligned}
v(r)\big|_{m=4} &= \sin{\left(\pi \left(r - 1\right) \right)}
\end{aligned}
$$

位相差を求める。

$$
\begin{aligned}
\delta_0\big|_{m=4} &= - \pi
\end{aligned}
$$

零次分だけの全断面積を計算する。

$$
\begin{aligned}
\sigma_0\big|_{m=4} &= 0
\end{aligned}
$$

観測点 $r=2$ での値と導関数を計算する。

$$
\begin{aligned}
v(r)\big|_{r=2} &= 0 \\
\frac{\partial^1 (v(r))}{\partial r^1} &= \pi \cos{\left(\pi \left(r - 1\right) \right)} \\
\left.\frac{\partial^1 (v(r))}{\partial r^1}\right|_{r=2} &= - \pi
\end{aligned}
$$
:::

## 演習問題

::::exercise{level=basic}
$u_{xx}+4u_{xy}+u_{yy}=0$ を分類せよ。
:::solution
$A=C=1,B=2$ なので $B^2-AC=3>0$ で双曲型である。係数行列の固有値は $3,-1$ で符号が異なる。
:::
::::

::::exercise{level=basic}
辺長 $2,1$ の長方形の第一 Dirichlet 固有値を求めよ。
:::solution
各方向の最低波数は $\pi/2,\pi$ であり、二乗和は $5\pi^2/4$ である。
:::
::::

::::exercise{level=basic}
熱方程式の一モード $u_n(0)=3$、$\lambda_n=2$、$\kappa=1$ の斉次解を求めよ。
:::solution
$u_n'+2u_n=0$ を解くと $u_n(t)=3e^{-2t}$ である。
:::
::::

::::exercise{level=basic}
$u_t+2u_x=0$、$u(0,x)=\cos x$ を解け。
:::solution
特性に沿って $x=\xi+2t$、$u$ は一定なので $u=\cos(x-2t)$ である。
:::
::::

::::exercise{level=advanced}
零 Neumann 境界の下で一様な正の源と Laplace の定常方程式が両立しないことを示せ。
:::solution
$\Delta u=f>0$ を体積積分すると右辺は正であるが、左辺は境界の流束であり零である。従って解は存在しない。
:::
::::

::::exercise{level=advanced}
波動方程式で初期変位零、初期速度 $g(x)=x$ の点ごとの解を求めよ。
:::solution
積分核により $(2c)^{-1}[(x+ct)^2-(x-ct)^2]/2=xt$ である。有限区間の式なので、$g$ が全実線で可積分である必要はない。直接代入しても両二階微分は零である。
:::
::::

::::exercise{level=advanced}
熱核の平均と分散を計算せよ。
:::solution
平均は偶関数に $x$ を掛けるため零である。$\partial_xG_t=-xG_t/(2\kappa t)$ だから $\int x^2G_t=-2\kappa t\int x\partial_xG_t=2\kappa t$ である。Gaussian 減衰で境界項は消える。
:::
::::

::::exercise{level=advanced}
入射 $\ell=0$ の球散乱で $ka=\pi$ のとき散乱係数を求めよ。
:::solution
$j_0(\pi)=0$ なので入射モードだけで境界値が零である。$h_0^{(1)}(\pi)=i/\pi\ne0$ だから散乱係数は零である。
:::
::::

## 確認テスト

::::quiz{answer="3" tol="1e-9"}
$A=C=1,B=2$ の $B^2-AC$ を答えよ。

:::explanation
$2^2-1\cdot1=3$ である。正なので主部は異符号の二方向を持つ。
:::
::::

::::quiz
熱の初期値問題で通常指定するものはどれか。

- [x] 初期温度
- [ ] 初期温度と独立の初期時間微分
- [ ] 初期時間微分だけ

:::explanation
各モードは時間一階で、初期温度を与えると方程式で初期時間微分も決まる。
:::
::::

::::quiz{answer="0" tol="1e-9"}
$u(0)=0$ で強制の無い零固有値の熱モードの値を答えよ。

:::explanation
式 :eqref[eq-mp-heat-modes] で $\lambda_n=0$、$f_n=0$ なら初期値から変わらず零である。
:::
::::

::::quiz
$e^{-i\omega t}$ のもとで外向きの球波はどれか。

- [ ] $e^{-ikr}/r$
- [x] $e^{ikr}/r$
- [ ] $\cos(kr)/r$

:::explanation
$kr-\omega t$ の位相が一定の点は外へ動く。余弦には両方の指数が含まれるので放射条件を満たさない。
:::
::::

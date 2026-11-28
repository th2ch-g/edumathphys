---
title: '測地線と指数写像'
description: '測地流と指数写像を構成し、Gauss の補題から局所的な最短性と強凸近傍を得る。'
notation:
  - symbol: '\gamma_{p,v}'
    meaning: '初期値 (p,v) の測地線'
  - symbol: '\exp_p'
    meaning: 'p における指数写像'
  - symbol: 'r_p'
    meaning: 'p からの距離'
---

## この章の目標

- 測地線方程式を接束上の常微分方程式に直す。
- 指数写像を使って正規座標を作る。
- Gauss の補題で距離と最短性を証明する。
- 始点も動かして強凸近傍を構成する。

## 測地線方程式と測地流

まっすぐ進むことは速度が進行方向に沿って変わらないことに対応する。曲線に沿う微分を使うと、座標に依存せずこの条件を書ける。

:::definition[測地線]{#def-rg-geodesic}
Levi-Civita 接続に関して $D_t\dot\gamma=0$ を満たす滑らかな曲線を :term[測地線]{en="geodesic" yomi="そくちせん"} という。定数曲線も含める。
:::

:ref[prop-rg-curve-derivative] より、$v^i=\dot x^i$ と置くと
$$
\dot x^i=v^i,\qquad \dot v^i=-\Gamma^i_{jk}(x)v^jv^k.
\label{eq-rg-geodesic-flow}
$$
右辺は接束上で滑らかであり、座標変換について整合することは共変微分の不変性から従う。ODE の局所存在、一意性、初期値への滑らかな依存性から最大解 $\gamma_{p,v}$ がある。接束上で $(p,v)\mapsto(\gamma_{p,v}(t),\dot\gamma_{p,v}(t))$ と送る局所流が :term[測地流]{en="geodesic flow" yomi="そくちりゅう"} である。計量の両立性から $\partial_t\norm{\dot\gamma}^2=2g(D_t\dot\gamma,\dot\gamma)=0$ である。従って測地流は各速度の球面束を保つ。

$\eta(s)=\gamma(t(s))$ では $D_s\eta'=(t')^2D_t\dot\gamma+t''\dot\gamma$ である。非定数の測地線の再パラメータ化が測地線となるのは $t''=0$ の場合に限る。特に一意性から $\gamma_{p,av}(t)=\gamma_{p,v}(at)$ が共通の定義域で成り立つ。

:::example[Euclid 平面の極座標]
$g=dr^2+r^2d\theta^2$ では $\ddot r-r\dot\theta^2=0$、$\ddot\theta+2\dot r\dot\theta/r=0$ である。直線 $(x,y)=(1,t)$ について $r=\sqrt{1+t^2}$、$\theta=\arctan t$ と置くと $\ddot r=(1+t^2)^{-3/2}=r\dot\theta^2$ である。さらに $r^2\dot\theta=1$ が二つ目の式を与える。
:::

:::example[回転対称計量の保存量]
$dr^2+f(r)^2d\theta^2$ で角成分の式を $f^2$ 倍すると $(f^2\dot\theta)'=0$ になる。したがって $c=f^2\dot\theta$ を固定すると $\dot r^2+c^2/f^2$ は速度の二乗で一定である。例えば $f=e^r$、単位速度、$c=2$ の測地線は $e^r\ge 2$ の領域に留まる。
:::

### 共役運動量と測地流

接ベクトルの成分 $v^i$ の代わりに $p_i=g_{ij}v^j$ と置く。これは計量による接束と余接束の同型である。余接束の関数 $H(x,p)=\frac12g^{ij}(x)p_i p_j$ を考える。逆行列の微分 $\partial_k g^{ij}=-g^{ia}(\partial_k g_{ab})g^{bj}$ から、式 :eqref[eq-rg-geodesic-flow] は
$$
\dot x^i=\frac{\partial H}{\partial p_i},\qquad
\dot p_i=-\frac{\partial H}{\partial x^i}
\label{eq-rg-geodesic-hamiltonian}
$$
となる。実際、$\dot p_i=\partial_k g_{ij}\,v^k v^j+g_{ij}\dot v^j$ に Christoffel 記号を代入すると、第 1 項は Koszul 公式の二つの対称項に相殺され、$\frac12\partial_i g_{jk}\,v^jv^k$ が残る。この式に $v^j=g^{ja}p_a$ を入れれば右辺である。

二組の方程式を合わせた余接束上の発散は、座標密度 $|dx^1\cdots dx^n\,dp_1\cdots dp_n|$ に関して零である。実際 $\partial_{x^i}\partial_{p_i}H-\partial_{p_i}\partial_{x^i}H=0$ である。余接座標を変えると、底の Jacobian と繊維の逆転置の Jacobian の行列式が打ち消し合うので、この密度は座標によらない。従って測地流はこの密度を保つ。これは [シンプレクティック多様体](/math/manifolds/symplectic/) で扱う Hamilton 流の体積保存の例である。また $H$ がある座標 $x^i$ に依存しない場合には $p_i$ が保存されるため、その方向の微分方程式を最初に積分できる。

単位接束上でも保存される密度がある。各余接空間で $p=G^{1/2}w$ とすると $H=|w|^2/2$ であり、繊維の座標変換の Jacobian は $\sqrt{\det G}$ である。$w=\rho u$ と極座標を取れば、余接束の密度は $\rho^{n-1}d\rho\,\mu_g\,d\sigma(u)$ と分解される。ここで $d\sigma$ は単位方向の球面密度である。底に依存する正規直交標構を替えても各球面の測度が変わらないため、$\mu_g\,d\sigma$ は単位接束で定義できる。測地流は $\rho$ を保つので、座標 $(\rho,x,u)$ での Jacobian は $\rho$ を固定した変換の Jacobian に等しい。全空間の密度を保つことから、特に $\rho=1$ で $\mu_g\,d\sigma$ を保つ。平坦トーラスでは、単位方向ごとに底の点が一定量だけ平行移動するので、この保存則を直接確かめられる。

:::example[上半平面で保存される水平運動量]
$H=y^2(p_x^2+p_y^2)/2$ なので $p_x$ は一定である。単位速度なら $2H=1$ であり、$p_x\ne0$ の測地線は $y\le1/|p_x|$ に留まる。例えば $p_x=1/2$ なら高さの上限は $2$ である。高さが上限に達する点では $p_y=0$、$\dot y=0$ で、水平に向く。$p_x=0$ なら $x$ が一定となり、垂直線上を $y=e^{\pm t}y_0$ に沿って進む。
:::

### エネルギーからの導出

座標内の経路の作用を $E(x)=\frac12\int g_{ij}\dot x^i\dot x^j\,dt$ とする。両端で零になる変分 $w$ を取り、積の微分を行うと
$$
\left.\partial_sE(x+sw)\right|_{s=0}
=\int\left(\tfrac12\partial_k g_{ij}\dot x^i\dot x^j w^k
+g_{kj}\dot x^j\dot w^k\right)dt.
$$
第 2 項を部分積分すると境界項が消える。任意の $w$ に対して微分が零となるには、その各成分に掛かる係数が零でなければならない。係数がある点で正なら、そこに台を持つ非負の $w^k$ を選んで矛盾するからである。従って
$$
\frac{d}{dt}(g_{kj}\dot x^j)-\tfrac12\partial_k g_{ij}\dot x^i\dot x^j=0.
$$
第 1 項の積を微分し、$\dot x^i\dot x^j$ の対称性で $\partial_i g_{kj}$ を二項の平均に置き換える。逆計量 $g^{\ell k}$ を掛けると座標の測地線方程式が得られる。逆に方程式が成立するなら積分は零となる。複数の座標域をまたぐ場合も、曲線の近くの変分を時間方向の分割で局所化できる。このとき座標内の条件は接続で書いた不変な方程式と同値なので、全区間で同じ条件になる。

:::example[双曲平面の半円を定速でたどる]
$g=(dx^2+dy^2)/y^2$ の上半平面で
$$
x(t)=b+a\tanh t,\qquad y(t)=a\,\operatorname{sech}t,\qquad a>0
$$
とする。$(x-b)^2+y^2=a^2$ であり、$\dot x=a\,\operatorname{sech}^2t$、$\dot y=-a\,\operatorname{sech}t\tanh t$ だから $\dot x/y^2=1/a$ で一定である。速度の二乗も $(\dot x^2+\dot y^2)/y^2=\operatorname{sech}^2t+\tanh^2t=1$ となる。$y$ の方程式は $y\ddot y+\dot x^2-\dot y^2=0$ である。ここで $\ddot y=y(\tanh^2t-\operatorname{sech}^2t)$ を代入すると三項が消える。例えば $a=2,b=1,t=\log2$ では $\tanh t=3/5$、$\operatorname{sech}t=4/5$ なので位置は $(11/5,8/5)$ である。
:::

### Lie 微分から得る保存量

場 $X$ が $\mathcal L_Xg=0$ を満たすとする。任意の $Y,Z$ について Lie 微分を展開し、$[X,Y]=\nabla_XY-\nabla_YX$ と両立式を使うと
$$
(\mathcal L_Xg)(Y,Z)=g(\nabla_YX,Z)+g(Y,\nabla_ZX)
$$
となる。測地線で $g(X,\dot\gamma)$ を微分すると $g(\nabla_{\dot\gamma}X,\dot\gamma)$ が残り、これは $(\mathcal L_Xg)(\dot\gamma,\dot\gamma)/2=0$ である。従って流れが等長的な方向場の成分は保存される。この議論では $X$ が特定の座標方向である必要はない。

:::example[球面の角運動量]
$S^2(a)$ 上で $X(q)=e_3\times q$ とする。流れは第 3 軸のまわりの回転なので等長である。$g(X,\dot q)=e_3\cdot(q\times\dot q)$ が一定となる。同じことを任意の軸に適用すると $q\times\dot q$ の成分がことごとく一定である。例えば $q(0)=(a,0,0)$、$\dot q(0)=(0,3/5,4/5)$ のとき、そのベクトルは $(0,-4a/5,3a/5)$ である。従って軌道は平面 $-4y+3z=0$ と球面の交わりに含まれる。加速度が球の法方向であることから直接微分しても同じ保存則を得る。
:::

:::example[円錐の頂点を除いた領域の角方向から出る測地線]
$g=dr^2+f^2d\theta^2$ に $f(r)=2 r$ を入れ、初期点 $(r,\theta)=(1,0)$ で角方向に単位速、動径速度零で出る。ここでは $r>0$ とする。角を $2\theta$ に替えると局所 Euclid 表示を得るが、角の周期は二倍になる。

角運動量は $c=f^2\dot\theta$、速度の二乗は $\dot r^2+f^2\dot\theta^2$ である。初期角速度は $1/f(r_0)$ だから $c=f(r_0)=2$ となる。後の時刻でも
$$
\begin{aligned}
\dot\theta&=\frac{2}{f(r)^2}\\
\dot r^2&=1-\frac{4}{f(r)^2}
\end{aligned}
$$
を満たす。従って到達可能な領域では $f(r)\ge f(r_0)$ が必要である。二乗した動径方程式だけでは、初期動径速度が零のときにどちらへ進むかは決まらない。二階の元の方程式 $\ddot r=ff'\dot\theta^2$ を使うと
$$
\ddot r(0)=f'(r_0)/f(r_0)=1,\qquad
\ddot\theta(0)=-2(f'(r_0)/f(r_0))\dot r(0)\dot\theta(0)=0
$$
となる。

$f'(r_0)\ne0$ なので円は測地線の像にならない。$r(\tau)=r_0+\frac12(1)\tau^2+O(\tau^3)$ であり、加速度は $f$ が増える側を向く。実際 $f(r(\tau))-f(r_0)=\frac{f'(r_0)^2}{2f(r_0)}\tau^2+O(\tau^3)$ で、保存量から得た到達可能条件を満たす。角方向に接する初期条件でも、測地線を円に拘束したことにはならない。
:::

:::example[カテノイドのくびれの外側の角方向から出る測地線]
$g=dr^2+f^2d\theta^2$ に $f(r)=\sqrt{r^{2} + 1}$ を入れ、初期点 $(r,\theta)=(1,0)$ で角方向に単位速、動径速度零で出る。半径は $r=0$ で最小になる。左右の同じ半径の円では $f'$ の符号が異なる。

角運動量は $c=f^2\dot\theta$、速度の二乗は $\dot r^2+f^2\dot\theta^2$ である。初期角速度は $1/f(r_0)$ だから $c=f(r_0)=\sqrt{2}$ となる。後の時刻でも
$$
\begin{aligned}
\dot\theta&=\frac{\sqrt{2}}{f(r)^2}\\
\dot r^2&=1-\frac{2}{f(r)^2}
\end{aligned}
$$
を満たす。従って到達可能な領域では $f(r)\ge f(r_0)$ が必要である。二乗した動径方程式だけでは、初期動径速度が零のときにどちらへ進むかは決まらない。二階の元の方程式 $\ddot r=ff'\dot\theta^2$ を使うと
$$
\ddot r(0)=f'(r_0)/f(r_0)=\frac{1}{2},\qquad
\ddot\theta(0)=-2(f'(r_0)/f(r_0))\dot r(0)\dot\theta(0)=0
$$
となる。

$f'(r_0)\ne0$ なので円は測地線の像にならない。$r(\tau)=r_0+\frac12(\frac{1}{2})\tau^2+O(\tau^3)$ であり、加速度は $f$ が増える側を向く。実際 $f(r(\tau))-f(r_0)=\frac{f'(r_0)^2}{2f(r_0)}\tau^2+O(\tau^3)$ で、保存量から得た到達可能条件を満たす。角方向に接する初期条件でも、測地線を円に拘束したことにはならない。
:::

:::example[双曲円柱の軸から離れた位置の角方向から出る測地線]
$g=dr^2+f^2d\theta^2$ に $f(r)=\cosh{\left(r \right)}$ を入れ、初期点 $(r,\theta)=(\log{\left(2 \right)},0)$ で角方向に単位速、動径速度零で出る。座標 $r$ は閉測地線の軸からの符号付き距離である。角方向を周期で割っても局所計算は同じである。

角運動量は $c=f^2\dot\theta$、速度の二乗は $\dot r^2+f^2\dot\theta^2$ である。初期角速度は $1/f(r_0)$ だから $c=f(r_0)=\frac{5}{4}$ となる。後の時刻でも
$$
\begin{aligned}
\dot\theta&=\frac{\frac{5}{4}}{f(r)^2}\\
\dot r^2&=1-\frac{\frac{25}{16}}{f(r)^2}
\end{aligned}
$$
を満たす。従って到達可能な領域では $f(r)\ge f(r_0)$ が必要である。二乗した動径方程式だけでは、初期動径速度が零のときにどちらへ進むかは決まらない。二階の元の方程式 $\ddot r=ff'\dot\theta^2$ を使うと
$$
\ddot r(0)=f'(r_0)/f(r_0)=\frac{3}{5},\qquad
\ddot\theta(0)=-2(f'(r_0)/f(r_0))\dot r(0)\dot\theta(0)=0
$$
となる。

$f'(r_0)\ne0$ なので円は測地線の像にならない。$r(\tau)=r_0+\frac12(\frac{3}{5})\tau^2+O(\tau^3)$ であり、加速度は $f$ が増える側を向く。実際 $f(r(\tau))-f(r_0)=\frac{f'(r_0)^2}{2f(r_0)}\tau^2+O(\tau^3)$ で、保存量から得た到達可能条件を満たす。角方向に接する初期条件でも、測地線を円に拘束したことにはならない。
:::

:::example[回転トーラスの内側の赤道の角方向から出る測地線]
$g=dr^2+f^2d\theta^2$ に $f(r)=\cos{\left(r \right)} + 2$ を入れ、初期点 $(r,\theta)=(\pi,0)$ で角方向に単位速、動径速度零で出る。内側の赤道では角方向の半径が最小になる。外側の赤道 $r=0$ と比べると二階微分の符号が逆である。

角運動量は $c=f^2\dot\theta$、速度の二乗は $\dot r^2+f^2\dot\theta^2$ である。初期角速度は $1/f(r_0)$ だから $c=f(r_0)=1$ となる。後の時刻でも
$$
\begin{aligned}
\dot\theta&=\frac{1}{f(r)^2}\\
\dot r^2&=1-\frac{1}{f(r)^2}
\end{aligned}
$$
を満たす。従って到達可能な領域では $f(r)\ge f(r_0)$ が必要である。二乗した動径方程式だけでは、初期動径速度が零のときにどちらへ進むかは決まらない。二階の元の方程式 $\ddot r=ff'\dot\theta^2$ を使うと
$$
\ddot r(0)=f'(r_0)/f(r_0)=0,\qquad
\ddot\theta(0)=-2(f'(r_0)/f(r_0))\dot r(0)\dot\theta(0)=0
$$
となる。

$f'(r_0)=0$ なので $r=r_0,\ \theta=\tau/f(r_0)$ は二階方程式を満たす。ODE の一意性からこの初期値では円を離れない。$f''(r_0)=1$ はその円の近くで角方向の半径がどのように変わるかを記録するが、円そのものの測地性の判定に使うのは一階の $f'$ である。
:::

:::example[指数的に縮む円柱の端の角方向から出る測地線]
$g=dr^2+f^2d\theta^2$ に $f(r)=e^{- \frac{r}{2}}$ を入れ、初期点 $(r,\theta)=(2 \log{\left(2 \right)},0)$ で角方向に単位速、動径速度零で出る。底の弧長を増やしたとき円周が一定の割合で縮む。座標ベクトルを使うと角方向の基底の長さが場所によって異なる。

角運動量は $c=f^2\dot\theta$、速度の二乗は $\dot r^2+f^2\dot\theta^2$ である。初期角速度は $1/f(r_0)$ だから $c=f(r_0)=\frac{1}{2}$ となる。後の時刻でも
$$
\begin{aligned}
\dot\theta&=\frac{\frac{1}{2}}{f(r)^2}\\
\dot r^2&=1-\frac{\frac{1}{4}}{f(r)^2}
\end{aligned}
$$
を満たす。従って到達可能な領域では $f(r)\ge f(r_0)$ が必要である。二乗した動径方程式だけでは、初期動径速度が零のときにどちらへ進むかは決まらない。二階の元の方程式 $\ddot r=ff'\dot\theta^2$ を使うと
$$
\ddot r(0)=f'(r_0)/f(r_0)=- \frac{1}{2},\qquad
\ddot\theta(0)=-2(f'(r_0)/f(r_0))\dot r(0)\dot\theta(0)=0
$$
となる。

$f'(r_0)\ne0$ なので円は測地線の像にならない。$r(\tau)=r_0+\frac12(- \frac{1}{2})\tau^2+O(\tau^3)$ であり、加速度は $f$ が増える側を向く。実際 $f(r(\tau))-f(r_0)=\frac{f'(r_0)^2}{2f(r_0)}\tau^2+O(\tau^3)$ で、保存量から得た到達可能条件を満たす。角方向に接する初期条件でも、測地線を円に拘束したことにはならない。
:::

## 指数写像

:::definition[指数写像]{#def-rg-exponential}
時刻 $1$ まで解が存在する $(p,v)$ の集合を $\mathcal D$ とする。:term[指数写像]{en="exponential map" yomi="しすうしゃぞう"} を $\exp_p(v)=\gamma_{p,v}(1)$ と定める。
:::

ODE の滑らかな依存性から $\mathcal D$ は零切断を含む開集合であり、$\exp$ は滑らかである。時刻のスケーリングにより各繊維の定義域は $0$ に関して星形であり、$\exp_p(tv)=\gamma_{p,v}(t)$ である。この式を $t=0$ で微分すると $(d\exp_p)_0=\id_{T_pM}$ となる。ここで $T_0(T_pM)$ をベクトル空間 $T_pM$ 自身と同一視した。

:::example[球面の指数写像]
球面 $S^n(a)$ 上で、接ベクトル $v$ が非零なら
$$
\exp_p(v)=\cos(\norm v/a)p+\frac{a}{\norm v}\sin(\norm v/a)v.
\label{eq-rg-sphere-exp}
$$
実際 $v$ を $tv$ に替えた曲線の Euclid 加速度は $-\norm v^2\gamma/a^2$ であり、接成分が零である。接成分への射影として得た接続は積の微分と捩率零を保つため、:ref[thm-rg-levi-civita] により Levi-Civita 接続である。$a=2$、$\norm v=\pi$ のときは像が $2v/\pi$ になる。
:::

:::example[開区間の定義域]
$(0,3)$ に $dx^2$ を入れる。$p=1$ では $\exp_p(v)=1+v$、定義域は $-1<v<2$ である。接空間全体は $\R$ だが、解が区間の外に出る速度では時刻 $1$ を使えない。
:::

### 積と定数倍

積計量の接続は因子ごとに作用する。各因子に依存する座標で交差 Christoffel 記号が零であることから、積の測地線方程式は二つの方程式に分かれる。従って初期値が $(p_1,p_2;v_1,v_2)$ のとき
$$
\exp_{(p_1,p_2)}(v_1,v_2)=(\exp_{p_1}v_1,\exp_{p_2}v_2).
$$
時刻 $1$ まで存在するには各成分で同時に存在することが必要十分であるから、繊維での定義域も直積である。一方、計量の定数倍 $\widetilde g=a^2g$ では、座標式で逆計量の $a^{-2}$ と偏微分の $a^2$ が消える。従って接ベクトルを固定した指数写像は同じだが、単位接ベクトルは元の $1/a$ 倍になる。

:::example[球面と開区間を同時に進む]
$S^2(2)\times(0,3)$ の基点 $((2,0,0),1)$ で初速度 $((0,\pi,0),3/2)$ を与える。球面の角速度は $\pi/2$ だから時刻 $1$ で $(0,2,0)$、区間の座標は $1+3/2=5/2$ となる。区間の制限から、この初期値の最大時刻区間は $(-2/3,4/3)$ である。積の長さは $1$ の間に一定速度 $\sqrt{\pi^2+9/4}$ を積分した値である。区間の初速度だけを $2$ に変えると時刻 $1$ には解が存在しない。
:::

:::example[双曲面の指数写像]
$B(p,p)=-a^2$ の双曲面で $B(v,p)=0$、$\ell^2=B(v,v)$ とする。$v=0$ では定数曲線である。$v\ne0$ のとき曲線
$$
\gamma(t)=\cosh(\ell t/a)p+\frac{a}{\ell}\sinh(\ell t/a)v
$$
は $B(\gamma,\gamma)=-a^2$、$\dot\gamma(0)=v$ を満たす。第 0 成分は途中で零にならないので同じ双曲面の成分に留まる。二階微分は $\ell^2\gamma/a^2$ で $B$ に直交する法方向だから、接加速度は零である。ここでも接射影による接続は対称性と両立式から Levi-Civita 接続である。$a=1$、$p=(1,0,0)$、$v=(0,\log2,0)$ なら $\exp_p(v)=(5/4,3/4,0)$ である。球の場合の三角関数と異なり、$\cosh$ は正の半径で狭義単調なので動径が重複しない。
:::

:::example[同じ垂直線上の二点の初期速度と指数写像]
計量 $(dx^2+dy^2)/y^2$ の上半平面で $p=\left( 0, \  1\right)$、$q=\left( 0, \  4\right)$ とする。同じ $x$ 座標なので測地線は垂直である。双曲面の計算が $\log y$ の弧長座標と一致するかを確認する。双曲面の Lorentz 内積を $B(V,W)=-V_0W_0+V_1W_1+V_2W_2$ とする。上半平面からの等長写像 :ref[prop-rg-half-space-isometry] を代入すると
$$
\begin{aligned}
P&=\begin{pmatrix}1\\0\\0\end{pmatrix}\\
Q&=\begin{pmatrix}\frac{17}{8}\\0\\\frac{15}{8}\end{pmatrix}\\
C=-B(P,Q)&=\frac{17}{8}
\end{aligned}
$$
になる。両点は $B(P,P)=B(Q,Q)=-1$ を満たし、第 0 成分が正なので同じ双曲面上にある。

$Q-CP$ は $B(P,Q-CP)=-C+C=0$ なので接ベクトルである。その二乗長は
$$
B(Q-CP,Q-CP)=-1+2C^2-C^2=C^2-1=\frac{225}{64}
$$
である。従って単位初期方向は $U=(Q-CP)/\sqrt{C^2-1}$ であり、この例では $\begin{pmatrix}0\\0\\1\end{pmatrix}$ となる。曲線
$$
\gamma(t)=\cosh t\,P+\sinh t\,U
$$
は $B(\gamma,\gamma)=-1$、$B(\dot\gamma,\dot\gamma)=1$ を満たす。さらに $\ddot\gamma=\gamma$ は周囲での法方向だから、接加速度は零である。$d=\log(C+\sqrt{C^2-1})$ とすると $C-\sqrt{C^2-1}=(C+\sqrt{C^2-1})^{-1}$ により $\cosh d=C$、$\sinh d=\sqrt{C^2-1}$ であり、$\gamma(d)=Q$ を得る。

この曲線は全体の最短弧である。実際、Lorentz 基底 $(P,U,\text{残りの単位接方向})$ を標準 Lorentz 基底に移すと、$P$ は頂点、$Q$ は動径距離 $d$ の点になる。極座標の速度は動径速度の絶対値以上であるから、どの経路も長さ $d$ 以上である。動径で等号が成立し、この例では $d=\log{\left(4 \right)}$ である。従って指定した初期ベクトルを指数写像に入れるなら $dU$、弧長 $t$ の初速度を指定するなら $U$ を用いる。

中点の計算には $P+Q=(1+\cosh d)P+\sinh d\,U$ を用いる。半角公式で右辺は $2\cosh(d/2)(\cosh(d/2)P+\sinh(d/2)U)$ になる。従って Lorentz 二乗長を $-1$ に正規化した
$$
M=\frac{P+Q}{\sqrt{2+2C}}=\begin{pmatrix}\frac{5}{4}\\0\\\frac{3}{4}\end{pmatrix}
$$
が中点である。上半平面へ戻す式 $y=1/(M_0-M_2)$、$x=M_1/(M_0-M_2)$ から中点は $\left( 0, \  2\right)$ となる。上半平面では垂直線 $x=0$ が軌跡である。弧長座標 $\log y$ での中点に対応するので、高さは二つの高さの幾何平均である。

同じ表示で一般の比率 $\lambda\in[0,1]$ の点は
$$
\gamma(\lambda d)
=\frac{\sinh((1-\lambda)d)}{\sinh d}P+
\frac{\sinh(\lambda d)}{\sinh d}Q
$$
である。$U=(Q-\cosh d\,P)/\sinh d$ を代入し、$\cosh(\lambda d)\sinh d-\sinh(\lambda d)\cosh d=\sinh((1-\lambda)d)$ を使った。二つの係数は通常の線形補間の $1-\lambda,\lambda$ とは異なる。これは $\gamma$ を双曲面上に保ちつつ弧長を比例させるためである。
:::

:::example[辺の一方が高さの異なる二点の初期速度と指数写像]
計量 $(dx^2+dy^2)/y^2$ の上半平面で $p=\left( 0, \  1\right)$、$q=\left( 2, \  2\right)$ とする。半円の中心を平面の中点で選ぶと、境界に直交する円にならない。二つの模型の計算を照合する。双曲面の Lorentz 内積を $B(V,W)=-V_0W_0+V_1W_1+V_2W_2$ とする。上半平面からの等長写像 :ref[prop-rg-half-space-isometry] を代入すると
$$
\begin{aligned}
P&=\begin{pmatrix}1\\0\\0\end{pmatrix}\\
Q&=\begin{pmatrix}\frac{9}{4}\\1\\\frac{7}{4}\end{pmatrix}\\
C=-B(P,Q)&=\frac{9}{4}
\end{aligned}
$$
になる。両点は $B(P,P)=B(Q,Q)=-1$ を満たし、第 0 成分が正なので同じ双曲面上にある。

$Q-CP$ は $B(P,Q-CP)=-C+C=0$ なので接ベクトルである。その二乗長は
$$
B(Q-CP,Q-CP)=-1+2C^2-C^2=C^2-1=\frac{65}{16}
$$
である。従って単位初期方向は $U=(Q-CP)/\sqrt{C^2-1}$ であり、この例では $\begin{pmatrix}0\\\frac{4 \sqrt{65}}{65}\\\frac{7 \sqrt{65}}{65}\end{pmatrix}$ となる。曲線
$$
\gamma(t)=\cosh t\,P+\sinh t\,U
$$
は $B(\gamma,\gamma)=-1$、$B(\dot\gamma,\dot\gamma)=1$ を満たす。さらに $\ddot\gamma=\gamma$ は周囲での法方向だから、接加速度は零である。$d=\log(C+\sqrt{C^2-1})$ とすると $C-\sqrt{C^2-1}=(C+\sqrt{C^2-1})^{-1}$ により $\cosh d=C$、$\sinh d=\sqrt{C^2-1}$ であり、$\gamma(d)=Q$ を得る。

この曲線は全体の最短弧である。実際、Lorentz 基底 $(P,U,\text{残りの単位接方向})$ を標準 Lorentz 基底に移すと、$P$ は頂点、$Q$ は動径距離 $d$ の点になる。極座標の速度は動径速度の絶対値以上であるから、どの経路も長さ $d$ 以上である。動径で等号が成立し、この例では $d=\log{\left(\frac{\sqrt{65}}{4} + \frac{9}{4} \right)}$ である。従って指定した初期ベクトルを指数写像に入れるなら $dU$、弧長 $t$ の初速度を指定するなら $U$ を用いる。

中点の計算には $P+Q=(1+\cosh d)P+\sinh d\,U$ を用いる。半角公式で右辺は $2\cosh(d/2)(\cosh(d/2)P+\sinh(d/2)U)$ になる。従って Lorentz 二乗長を $-1$ に正規化した
$$
M=\frac{P+Q}{\sqrt{2+2C}}=\begin{pmatrix}\frac{\sqrt{26}}{4}\\\frac{\sqrt{26}}{13}\\\frac{7 \sqrt{26}}{52}\end{pmatrix}
$$
が中点である。上半平面へ戻す式 $y=1/(M_0-M_2)$、$x=M_1/(M_0-M_2)$ から中点は $\left( \frac{2}{3}, \  \frac{\sqrt{26}}{3}\right)$ となる。上半平面では軌跡は中心 $c=\frac{7}{4}$、半径 $\frac{\sqrt{65}}{4}$ の半円である。中心は両端までの Euclid 二乗距離を等しいと置いた式から求めた。中点を代入しても $(x-c)^2+y^2=R^2$ となる。

同じ表示で一般の比率 $\lambda\in[0,1]$ の点は
$$
\gamma(\lambda d)
=\frac{\sinh((1-\lambda)d)}{\sinh d}P+
\frac{\sinh(\lambda d)}{\sinh d}Q
$$
である。$U=(Q-\cosh d\,P)/\sinh d$ を代入し、$\cosh(\lambda d)\sinh d-\sinh(\lambda d)\cosh d=\sinh((1-\lambda)d)$ を使った。二つの係数は通常の線形補間の $1-\lambda,\lambda$ とは異なる。これは $\gamma$ を双曲面上に保ちつつ弧長を比例させるためである。
:::

:::example[垂直線を途中の点からたどる初期速度と指数写像]
計量 $(dx^2+dy^2)/y^2$ の上半平面で $p=\left( 3, \  \frac{1}{3}\right)$、$q=\left( 3, \  3\right)$ とする。垂直線が原点でない位置にある場合でも、水平方向の平行移動は等長である。双曲面の Lorentz 内積を $B(V,W)=-V_0W_0+V_1W_1+V_2W_2$ とする。上半平面からの等長写像 :ref[prop-rg-half-space-isometry] を代入すると
$$
\begin{aligned}
P&=\begin{pmatrix}\frac{91}{6}\\9\\\frac{73}{6}\end{pmatrix}\\
Q&=\begin{pmatrix}\frac{19}{6}\\1\\\frac{17}{6}\end{pmatrix}\\
C=-B(P,Q)&=\frac{41}{9}
\end{aligned}
$$
になる。両点は $B(P,P)=B(Q,Q)=-1$ を満たし、第 0 成分が正なので同じ双曲面上にある。

$Q-CP$ は $B(P,Q-CP)=-C+C=0$ なので接ベクトルである。その二乗長は
$$
B(Q-CP,Q-CP)=-1+2C^2-C^2=C^2-1=\frac{1600}{81}
$$
である。従って単位初期方向は $U=(Q-CP)/\sqrt{C^2-1}$ であり、この例では $\begin{pmatrix}- \frac{89}{6}\\-9\\- \frac{71}{6}\end{pmatrix}$ となる。曲線
$$
\gamma(t)=\cosh t\,P+\sinh t\,U
$$
は $B(\gamma,\gamma)=-1$、$B(\dot\gamma,\dot\gamma)=1$ を満たす。さらに $\ddot\gamma=\gamma$ は周囲での法方向だから、接加速度は零である。$d=\log(C+\sqrt{C^2-1})$ とすると $C-\sqrt{C^2-1}=(C+\sqrt{C^2-1})^{-1}$ により $\cosh d=C$、$\sinh d=\sqrt{C^2-1}$ であり、$\gamma(d)=Q$ を得る。

この曲線は全体の最短弧である。実際、Lorentz 基底 $(P,U,\text{残りの単位接方向})$ を標準 Lorentz 基底に移すと、$P$ は頂点、$Q$ は動径距離 $d$ の点になる。極座標の速度は動径速度の絶対値以上であるから、どの経路も長さ $d$ 以上である。動径で等号が成立し、この例では $d=\log{\left(9 \right)}$ である。従って指定した初期ベクトルを指数写像に入れるなら $dU$、弧長 $t$ の初速度を指定するなら $U$ を用いる。

中点の計算には $P+Q=(1+\cosh d)P+\sinh d\,U$ を用いる。半角公式で右辺は $2\cosh(d/2)(\cosh(d/2)P+\sinh(d/2)U)$ になる。従って Lorentz 二乗長を $-1$ に正規化した
$$
M=\frac{P+Q}{\sqrt{2+2C}}=\begin{pmatrix}\frac{11}{2}\\3\\\frac{9}{2}\end{pmatrix}
$$
が中点である。上半平面へ戻す式 $y=1/(M_0-M_2)$、$x=M_1/(M_0-M_2)$ から中点は $\left( 3, \  1\right)$ となる。上半平面では垂直線 $x=3$ が軌跡である。弧長座標 $\log y$ での中点に対応するので、高さは二つの高さの幾何平均である。

同じ表示で一般の比率 $\lambda\in[0,1]$ の点は
$$
\gamma(\lambda d)
=\frac{\sinh((1-\lambda)d)}{\sinh d}P+
\frac{\sinh(\lambda d)}{\sinh d}Q
$$
である。$U=(Q-\cosh d\,P)/\sinh d$ を代入し、$\cosh(\lambda d)\sinh d-\sinh(\lambda d)\cosh d=\sinh((1-\lambda)d)$ を使った。二つの係数は通常の線形補間の $1-\lambda,\lambda$ とは異なる。これは $\gamma$ を双曲面上に保ちつつ弧長を比例させるためである。
:::

## Gauss の補題

指数写像で半径と向きを変えたときの微分を区別する。計量に関して両者が直交すれば、向きを変える曲線の長さは半径だけを変える場合以上になる。

:::lemma[Gauss の補題]{#lem-rg-gauss}
$v\in\mathcal D\cap T_pM$、$w\in T_pM$ に対し
$$
g\bigl((d\exp_p)_v v,(d\exp_p)_v w\bigr)=g_p(v,w).
\label{eq-rg-gauss-lemma}
$$
:::

:::proof
$\alpha(s,t)=\exp_p(t(v+sw))$、$T=\partial_t\alpha$、$S=\partial_s\alpha$ と置く。各 $s$ について $D_tT=0$、捩率零から $D_tS=D_sT$ である。それゆえ
$$
\partial_t g(S,T)=g(D_sT,T)=\tfrac12\partial_s g(T,T)=g_p(v+sw,w).
$$
最後は測地線の速度が一定であることを用いた。$S(s,0)=0$ なので $t=0$ から $1$ まで積分し、$s=0$ とすればよい。
:::

### 極座標の計量

指数写像が局所微分同相の領域で $v=ru$、$\norm u=1$ と書くと $\partial_r$ は単位で、球面に接する $w$ の像と直交する。従って計量は $dr^2+h_r$ と書ける。曲線の速さは $|\dot r|$ 以上であり、$r>0$ で等号なら角方向の速度は零である。

:::example[球面の角方向の縮小]
式 :eqref[eq-rg-sphere-exp]で固定した $r$ に対して単位方向 $u$ を変えると、微分の長さは $a|\sin(r/a)|\norm{\dot u}$ である。動径方向の長さは $1$ である。これに対し $n\ge2$ では、角方向の倍率は $r=\pi a/2$ で $a$、$r=\pi a$ で $0$ になる。従って後者では指数写像は局所微分同相でない。$S^1(a)$ には角方向がなく、対蹠点への二つの初期ベクトルでも微分は非特異である。
:::

:::example[双曲面で動径と角を分ける]
単位双曲面で基点 $(1,0)$ を取ると、正規極座標は $(\cosh r,\sinh r\,u)$ である。動径微分は $(\sinh r,\cosh r\,u)$、角微分は $(0,\sinh r\,w)$ で、$u\cdot w=0$ だから両者は直交する。$r=\log2$ では角方向の倍率は $\sinh(\log2)=3/4$ である。同じ $r$ の球面では $\sin r$ となり、原点近くでは両者とも一階で $r$ に一致するが、三階の符号が逆である。
:::

### 動径を保つ逆像の情報

正規座標の接空間を $x=(x^1,\ldots,x^n)$ とし、標準内積で成分を下げる。Gauss の補題に $v=x$ と任意の $w$ を入れると、引き戻した計量は $g_{ij}(x)x^j=x^i$ を満たす。この式から $F(x)=|x|^2/2$ の勾配は $x^i\partial_i$ になる。実際、そのベクトルを任意の $Y$ に内積で評価すると $x^iY^i=dF(Y)$ である。このようにして逆計量を成分で求めずに $F$ を微分できる。

動径からの角方向のずれを評価するには、正規極座標で $h_r$ を角方向の二次形式とする。任意の区分的滑らかな曲線で
$$
\norm{\dot\gamma}^2=\dot r^2+h_r(\dot u,\dot u)
$$
である。角方向の項は半径が零でない範囲で正定値であるため、その区間で $\dot u$ が一度でも非零なら、開時間区間で長さは半径の全変動より厳密に大きくなる。特に一定半径の球面上では、動径距離の変化が零でも曲線長は角計量で正に測られる。

:::example[球面上の円の長さと半径方向の尺度]
単位球面の正規極座標で $r=\pi/3$ の円を一周する。$h_r=\sin^2 r\,d\theta^2$ なので長さは $2\pi\sin(\pi/3)=\sqrt3\pi$ である。同じ半径の Euclid 円の長さは $2\pi^2/3$ である。半径方向の長さが一致しても、固定した角の差を持つ二点を円に沿って結ぶ長さは一致しない。角成分に $\sin r/r$ が現れ、$0<r<\pi$ ではこの比が $1$ 未満であるためである。
:::

:::example[単位半円の異なる側の初期速度と指数写像]
計量 $(dx^2+dy^2)/y^2$ の上半平面で $p=\left( - \frac{3}{5}, \  \frac{4}{5}\right)$、$q=\left( \frac{3}{5}, \  \frac{4}{5}\right)$ とする。二点は半径 $1$ の円上にあり、対称性から測地線の中点はその頂点になる。双曲面の Lorentz 内積を $B(V,W)=-V_0W_0+V_1W_1+V_2W_2$ とする。上半平面からの等長写像 :ref[prop-rg-half-space-isometry] を代入すると
$$
\begin{aligned}
P&=\begin{pmatrix}\frac{5}{4}\\- \frac{3}{4}\\0\end{pmatrix}\\
Q&=\begin{pmatrix}\frac{5}{4}\\\frac{3}{4}\\0\end{pmatrix}\\
C=-B(P,Q)&=\frac{17}{8}
\end{aligned}
$$
になる。両点は $B(P,P)=B(Q,Q)=-1$ を満たし、第 0 成分が正なので同じ双曲面上にある。

$Q-CP$ は $B(P,Q-CP)=-C+C=0$ なので接ベクトルである。その二乗長は
$$
B(Q-CP,Q-CP)=-1+2C^2-C^2=C^2-1=\frac{225}{64}
$$
である。従って単位初期方向は $U=(Q-CP)/\sqrt{C^2-1}$ であり、この例では $\begin{pmatrix}- \frac{3}{4}\\\frac{5}{4}\\0\end{pmatrix}$ となる。曲線
$$
\gamma(t)=\cosh t\,P+\sinh t\,U
$$
は $B(\gamma,\gamma)=-1$、$B(\dot\gamma,\dot\gamma)=1$ を満たす。さらに $\ddot\gamma=\gamma$ は周囲での法方向だから、接加速度は零である。$d=\log(C+\sqrt{C^2-1})$ とすると $C-\sqrt{C^2-1}=(C+\sqrt{C^2-1})^{-1}$ により $\cosh d=C$、$\sinh d=\sqrt{C^2-1}$ であり、$\gamma(d)=Q$ を得る。

この曲線は全体の最短弧である。実際、Lorentz 基底 $(P,U,\text{残りの単位接方向})$ を標準 Lorentz 基底に移すと、$P$ は頂点、$Q$ は動径距離 $d$ の点になる。極座標の速度は動径速度の絶対値以上であるから、どの経路も長さ $d$ 以上である。動径で等号が成立し、この例では $d=\log{\left(4 \right)}$ である。従って指定した初期ベクトルを指数写像に入れるなら $dU$、弧長 $t$ の初速度を指定するなら $U$ を用いる。

中点の計算には $P+Q=(1+\cosh d)P+\sinh d\,U$ を用いる。半角公式で右辺は $2\cosh(d/2)(\cosh(d/2)P+\sinh(d/2)U)$ になる。従って Lorentz 二乗長を $-1$ に正規化した
$$
M=\frac{P+Q}{\sqrt{2+2C}}=\begin{pmatrix}1\\0\\0\end{pmatrix}
$$
が中点である。上半平面へ戻す式 $y=1/(M_0-M_2)$、$x=M_1/(M_0-M_2)$ から中点は $\left( 0, \  1\right)$ となる。上半平面では軌跡は中心 $c=0$、半径 $1$ の半円である。中心は両端までの Euclid 二乗距離を等しいと置いた式から求めた。中点を代入しても $(x-c)^2+y^2=R^2$ となる。

同じ表示で一般の比率 $\lambda\in[0,1]$ の点は
$$
\gamma(\lambda d)
=\frac{\sinh((1-\lambda)d)}{\sinh d}P+
\frac{\sinh(\lambda d)}{\sinh d}Q
$$
である。$U=(Q-\cosh d\,P)/\sinh d$ を代入し、$\cosh(\lambda d)\sinh d-\sinh(\lambda d)\cosh d=\sinh((1-\lambda)d)$ を使った。二つの係数は通常の線形補間の $1-\lambda,\lambda$ とは異なる。これは $\gamma$ を双曲面上に保ちつつ弧長を比例させるためである。
:::

## 正規座標と正規近傍

:::definition[正規座標]{#def-rg-normal-coordinate}
$\exp_p$ が原点の星形近傍から像に微分同相となるとき、その像を :term[正規近傍]{en="normal neighborhood" yomi="せいききんぼう"} という。$T_pM$ の正規直交基底で $\exp_p^{-1}$ の成分を取る座標を :term[正規座標]{en="normal coordinates" yomi="せいきざひょう"} という。
:::

$(d\exp_p)_0=\id$ と逆関数定理により、十分小さい開球で正規座標が作れる。座標原点では $g_{ij}=\delta_{ij}$ である。任意方向の測地線が $x^i(t)=tv^i$ を満たすので :eqref[eq-rg-geodesic-flow] より $\Gamma^i_{jk}(0)v^jv^k=0$ である。下添字の対称性と偏極から $\Gamma^i_{jk}(0)=0$ であり、両立式から $\partial_k g_{ij}(0)=0$ を得る。逆に基点でこれらの等式が満たされるだけでは座標線が動径測地線とは限らない。例えば直線上で $y=x+x^3$ とすれば原点で $dy/dx=1$、$d^2y/dx^2=0$ だが $y$ は弧長ではない。

:::example[球面で失われる単射性]
式 :eqref[eq-rg-sphere-exp]では $0<r<\pi a$ の範囲に対し、$p$ との内積から $r$ が決まり、$v$ の方向は接成分から決まる。従って $\exp_p$ は半径 $\pi a$ の開球から対蹠点を除いた球面への微分同相である。境界 $r=\pi a$ はすべて対蹠点に写る。
:::

:::example[非線形座標での正規化]
直線に $g=e^{2x}dx^2$ を入れる。$x=0$ を基点とすれば $y=e^x-1$ は $dy=e^x dx$ を満たす弧長座標で、$dy/dx=1$ が原点で成り立つので正規座標である。指数写像は $v>-1$ に対し $x=\log(1+v)$ となる。これに対し $y=x+x^2/2$ は同じ一次と二次の微分を原点で持つが、$dy/dx=1+x$ が $e^x$ と一致しないので原点近傍の正規座標ではない。
:::

### 正規座標での球面の計量行列

単位球面の動径計量は $dr^2+\sin^2r\,g_{S^{n-1}}$ である。一方、接空間の Euclid 計量は $dr^2+r^2g_{S^{n-1}}$ である。$x=ru$ が零でないとき、接ベクトル $w$ を $w_\parallel=\inner{w}{u}u$ と $w_\perp=w-w_\parallel$ に分けると、正規座標の計量は
$$
g_x(w,w)=|w_\parallel|^2+\left(\frac{\sin r}{r}\right)^2|w_\perp|^2.
$$
動径成分の倍率は $1$、角方向では $du=w_\perp/r$ だからである。従って添字表示は
$$
g_{ij}(x)=\delta_{ij}-\tfrac13(r^2\delta_{ij}-x^i x^j)+O(r^4)
$$
となる。ここで $\sin r/r=1-r^2/6+O(r^4)$ を二乗した。$x=0$ での一階微分は全て零で、二階微分には球面の曲がり方が残る。単位双曲面では $\sin$ を $\sinh$ に変えるので、二次項の符号が逆になる。

:::example[変換で消せる一次項と残る二次項]
二次元球面で $x=(r,0)$ とする。$g_{11}=1$、$g_{12}=0$、$g_{22}=(\sin r/r)^2$ である。$r=\pi/2$ なら $g_{22}=4/\pi^2$ となり、同じ座標幅の二方向の長さが異なる。中心で $(\partial_1^2g_{22})(0)=-2/3$ であり、原点の基底だけ正規直交させた条件ではこの値を決められない。
:::

### 超曲面に垂直な測地座標

埋め込まれた超曲面 $\Sigma$ の一点の近くで単位法線 $\nu$ を選ぶ。写像 $\Phi(t,y)=\exp_y(t\nu_y)$ の $(0,y)$ での微分は、$\partial_t$ を $\nu_y$ へ、$\Sigma$ の接方向をそのままの接ベクトルへ送る。従って逆関数定理により、$\Sigma$ の局所座標と $t$ を組み合わせて周囲の座標が得られる。この表示を :term[半測地座標]{en="semigeodesic coordinates" yomi="はんそくちざひょう"} という。

$T=\Phi_*\partial_t$、$Y_a=\Phi_*\partial_{y^a}$ とすると $g(T,T)=1$、$\nabla_TT=0$ である。また座標場は可換で、捩率が零なので
$$
Tg(T,Y_a)=g(T,\nabla_TY_a)=g(T,\nabla_{Y_a}T)
=\tfrac12Y_ag(T,T)=0.
$$
$t=0$ では $g(T,Y_a)=0$ だから、この内積は座標域全体で零となる。従って計量は $dt^2+h_{ab}(t,y)dy^a dy^b$ と書け、$\grad t=T$ である。座標域内で $\Sigma$ から $(t,y)$ へ至る曲線の長さは $|t|$ 以上であり、$y$ を固定して進むと等号になる。これは座標域内での最短性である。超曲面が別のところから近づいている場合は、超曲面全体への距離をこの局所表示だけで決めることはできない。

:::example[双曲平面で水平線から垂直に進む]
上半平面で $\Sigma=\{y=1\}$ を取り、上向きの単位法線から出ると $(x,y)=(u,e^t)$ となる。この座標で $g=dt^2+e^{-2t}du^2$ である。$s=\log y$ は全上半平面で単位勾配を持つため、水平線のどの点から出ても長さは $|t|$ 以上となる。この例では座標域の外への近道を考える必要もなく、$\Sigma$ までの距離は全域で $|\log y|$ である。
:::

## 測地線の局所極小性

:::theorem[動径測地線の最短性]{#thm-rg-radial-minimal}
$\exp_p$ が半径 $R$ より大きい開球上で微分同相であり、$q=\exp_p(v)$、$\norm v<R$ とする。このとき $d_g(p,q)=\norm v$ で、長さを実現する区分的 $C^1$ 曲線の像は動径測地線の像に一致する。
:::

:::proof
半径 $R$ の像内では :ref[lem-rg-gauss] により長さは半径の全変動以上である。原点で極座標が使えない場合は、曲線が半径 $\epsilon$ の球を最後に出る時刻から積分し $\epsilon\to0$ とすれば下界が得られる。半径 $R$ の像から出る曲線の最初の境界到達までの長さは同じ議論で $R$ 以上である。動径は長さ $\norm v$ を実現する。

等号となる曲線では半径は単調非減少であり、正の半径の部分で角速度が零である。実際、等号のためには $|\dot r|$ の積分と半径差が等しく、また速さと $|\dot r|$ の差の積分が零でなければならない。差は各滑らかな区間で非負連続なので零である。従って動径の向きが一定である。
:::

測地線をある時刻から進めると、十分短い弧はこの定理で最短となる。始点をその時刻の前後に動かしても同じ下界を使えることは、次節の逆関数定理の議論で示す。この一様性を使うと、定速でパラメータ化した最短曲線は内部の各時刻をまたぐ小区間で動径測地線に一致する。実際、その区間が最短でなければ、より短い弧に取り替えて全体を短縮できる。従って最短曲線は滑らかな測地線であり、異なる速度の向きをつないだ折れ角は残らない。

:::example[長い弧と短い弧]
単位円の角 $0$ と $\pi/2$ を結ぶ二つの弧は、それぞれ長さ $\pi/2$ と $3\pi/2$ の測地線である。長い弧の任意の十分短い部分は最短であるが、端点全体では反対回りの弧が短い。
:::

:::example[長さと同じ像の再パラメータ化]
Euclid 直線の $\gamma(t)=t^2$、$0\le t\le1$ は端点を結ぶ曲線の最短長 $1$ を実現する。しかし共変加速度は $2$ で、与えられたパラメータで測地線ではない。動径像の一意性とパラメータの一意性を区別するには、速度を一定に固定する。パラメータ区間が $[0,1]$ なら、同じ向きの定速最短弧は $\gamma(t)=t$ に決まる。
:::

### 角のある二本の弧の短縮

二本の単位速測地線が $q$ でつながり、入射速度 $u$ と出射速度 $v$ が異なるとする。$q$ の正規座標で、$q$ の手前と先の距離 $\epsilon$ の点は $-\epsilon u,\epsilon v$ である。これらを座標線分で結ぶと、線分は座標半径 $\epsilon$ の球に入り、その Euclid 長は $\epsilon|u+v|$ である。$u\ne v$ なら $|u+v|^2=2+2g_q(u,v)<4$ である。連続性により球内では $\norm{w}_g\le(1+o(1))|w|$ とできるので、十分小さい $\epsilon$ で線分の Riemann 長は $2\epsilon$ 未満となる。$u=-v$ なら二端点は同じであり、定数曲線への置換が短縮になる。従って最短曲線に真の折れ角はない。この議論は特に切断点で異なる最短弧を取り替えるときに使う。

:::example[折れ角の大きさを長さに変える]
Euclid 平面で $q=0$、$u=e_1$、$v=e_2$ とする。折れ線は長さ $2\epsilon$、対角線は $\sqrt2\,\epsilon$ であり、短縮は $(2-\sqrt2)\epsilon$ である。Riemann 計量でも $q$ で正規直交する二方向なら同じ一次の差が得られる。この差は $\epsilon$ に比例するため、原点近くの計量成分の誤差から出る $o(\epsilon)$ より大きい。
:::

### エネルギーとしての一意性

時間幅 $T$ を固定すると、任意の曲線のエネルギーは $L^2/(2T)\ge d^2/(2T)$ である。動径測地線を定速で進めば、この両方の不等式で等号となる。エネルギー最小の曲線では長さの等号から動径像、Cauchy–Schwarz の等号から一定速度が得られるので、定数曲線を含めパラメータまで決まる。初期点から距離 $3$ の点が正規最短性の範囲にあるとき、時間幅 $2$ の最小エネルギーは $9/4$ である。幅 $1$ に変えると $9/2$ になるが、長さ自体は $3$ である。

:::example[ともに高い位置にある二点の初期速度と指数写像]
計量 $(dx^2+dy^2)/y^2$ の上半平面で $p=\left( 1, \  2\right)$、$q=\left( 3, \  4\right)$ とする。水平距離と高さの変化がともにある。Euclid の線分と測地線の中点を区別する。双曲面の Lorentz 内積を $B(V,W)=-V_0W_0+V_1W_1+V_2W_2$ とする。上半平面からの等長写像 :ref[prop-rg-half-space-isometry] を代入すると
$$
\begin{aligned}
P&=\begin{pmatrix}\frac{3}{2}\\\frac{1}{2}\\1\end{pmatrix}\\
Q&=\begin{pmatrix}\frac{13}{4}\\\frac{3}{4}\\3\end{pmatrix}\\
C=-B(P,Q)&=\frac{3}{2}
\end{aligned}
$$
になる。両点は $B(P,P)=B(Q,Q)=-1$ を満たし、第 0 成分が正なので同じ双曲面上にある。

$Q-CP$ は $B(P,Q-CP)=-C+C=0$ なので接ベクトルである。その二乗長は
$$
B(Q-CP,Q-CP)=-1+2C^2-C^2=C^2-1=\frac{5}{4}
$$
である。従って単位初期方向は $U=(Q-CP)/\sqrt{C^2-1}$ であり、この例では $\begin{pmatrix}\frac{2 \sqrt{5}}{5}\\0\\\frac{3 \sqrt{5}}{5}\end{pmatrix}$ となる。曲線
$$
\gamma(t)=\cosh t\,P+\sinh t\,U
$$
は $B(\gamma,\gamma)=-1$、$B(\dot\gamma,\dot\gamma)=1$ を満たす。さらに $\ddot\gamma=\gamma$ は周囲での法方向だから、接加速度は零である。$d=\log(C+\sqrt{C^2-1})$ とすると $C-\sqrt{C^2-1}=(C+\sqrt{C^2-1})^{-1}$ により $\cosh d=C$、$\sinh d=\sqrt{C^2-1}$ であり、$\gamma(d)=Q$ を得る。

この曲線は全体の最短弧である。実際、Lorentz 基底 $(P,U,\text{残りの単位接方向})$ を標準 Lorentz 基底に移すと、$P$ は頂点、$Q$ は動径距離 $d$ の点になる。極座標の速度は動径速度の絶対値以上であるから、どの経路も長さ $d$ 以上である。動径で等号が成立し、この例では $d=\log{\left(\frac{\sqrt{5}}{2} + \frac{3}{2} \right)}$ である。従って指定した初期ベクトルを指数写像に入れるなら $dU$、弧長 $t$ の初速度を指定するなら $U$ を用いる。

中点の計算には $P+Q=(1+\cosh d)P+\sinh d\,U$ を用いる。半角公式で右辺は $2\cosh(d/2)(\cosh(d/2)P+\sinh(d/2)U)$ になる。従って Lorentz 二乗長を $-1$ に正規化した
$$
M=\frac{P+Q}{\sqrt{2+2C}}=\begin{pmatrix}\frac{19 \sqrt{5}}{20}\\\frac{\sqrt{5}}{4}\\\frac{4 \sqrt{5}}{5}\end{pmatrix}
$$
が中点である。上半平面へ戻す式 $y=1/(M_0-M_2)$、$x=M_1/(M_0-M_2)$ から中点は $\left( \frac{5}{3}, \  \frac{4 \sqrt{5}}{3}\right)$ となる。上半平面では軌跡は中心 $c=5$、半径 $2 \sqrt{5}$ の半円である。中心は両端までの Euclid 二乗距離を等しいと置いた式から求めた。中点を代入しても $(x-c)^2+y^2=R^2$ となる。

同じ表示で一般の比率 $\lambda\in[0,1]$ の点は
$$
\gamma(\lambda d)
=\frac{\sinh((1-\lambda)d)}{\sinh d}P+
\frac{\sinh(\lambda d)}{\sinh d}Q
$$
である。$U=(Q-\cosh d\,P)/\sinh d$ を代入し、$\cosh(\lambda d)\sinh d-\sinh(\lambda d)\cosh d=\sinh((1-\lambda)d)$ を使った。二つの係数は通常の線形補間の $1-\lambda,\lambda$ とは異なる。これは $\gamma$ を双曲面上に保ちつつ弧長を比例させるためである。
:::

## 凸近傍の存在（Whitehead の定理）

:::definition[強凸近傍]{#def-rg-strongly-convex}
任意の二点を結ぶ最短測地線が多様体全体で一意に存在し、その像が近傍に含まれるとき、その開近傍を :term[強凸近傍]{en="strongly convex neighborhood" yomi="きょうとつきんぼう"} という。
:::

:::theorem[Whitehead の凸近傍定理]{#thm-rg-whitehead}
各点は任意に小さい強凸近傍を持つ。十分小さい距離球をそのような近傍にできる。
:::

:::proof
$\Psi(x,v)=(x,\exp_xv)$ を考える。$(p,0)$ での微分は適当な自明化で $(\xi,\eta)\mapsto(\xi,\xi+\eta)$ である。逆関数定理から、この点の近くで $\Psi$ は微分同相となる。自明化を縮めると、近くの全ての $x$ について半径 $R>0$ の閉球がこの定義域に入り、$(x,y)$ が $(p,p)$ に十分近いとき逆写像のベクトルは長さ $R/2$ 未満となる。従って :ref[thm-rg-radial-minimal] により近い $x,y$ を結ぶ最短弧は一意であり、端点に滑らかに依存する。

$p$ の正規座標で $F(x)=\norm{\exp_p^{-1}(x)}^2$ と置く。$\operatorname{Hess}_pF=2g_p$ なので連続性により近傍 $W$ で $\operatorname{Hess} F$ が正定値である。弧 $(x,y,t)\mapsto\exp_x(t\Psi^{-1}_2(x,y))$ は連続で、$(p,p,t)$ では常に $p$ である。区間 $0\le t\le1$ のコンパクト性により $x,y$ をさらに $p$ に近づければ弧全体が $W$ にある。弧に沿い $(F\circ\gamma)''=\operatorname{Hess}F(\dot\gamma,\dot\gamma)\ge0$ だから、途中の値は端点の値の最大以下である。従って十分小さい $F<\epsilon^2$ で弧が留まる。これは動径の最短性により距離球である。
:::

:::example[強凸球の半径]
単位球面で半径 $\epsilon<\pi/2$ の北極球は $\inner{p}{x}>\cos\epsilon$ である。二点 $x,y$ は対蹠点でなく、大円の短い弧は $(1-t)x+ty$ を正規化して得られる。分子の $p$ との内積は $\cos\epsilon$ より大きく、分母は $1$ 以下なので弧は球内に留まる。一方、半径 $\pi/2$ より大きい球は赤道上の対蹠点を含み最短弧が一意でない。
:::

:::example[円の強凸球の境界]
単位円で $p$ から距離 $\epsilon$ の開球を $(-\epsilon,\epsilon)$ と表す。$\epsilon\le\pi/2$ なら球内の任意の二角の差の絶対値は $\pi$ 未満だから、その区間内の弧が一意な最短弧である。$\epsilon>\pi/2$ なら $\pm\pi/2$ が球内にあり、互いに対蹠点なので最短弧は二つある。従って強凸な開球の半径はこの場合ちょうど $\pi/2$ までである。
:::

### 強凸領域上の補間と縮約

同じ強凸な距離球 $U$ を取り、その中の二点 $x,y$ を結ぶ一意な定速最短測地線を $\eta(x,y,t)$ と書く。Whitehead の定理の証明で選んだ球では、この写像は $x,y,t$ に滑らかである。$\eta(x,y,t)$ は $U$ の中にあるので、$p\in U$ を固定し $H(x,t)=\eta(p,x,1-t)$ と置けば、$U$ は $p$ へ滑らかに縮約される。$H(p,t)=p$ なので強い変形収縮である。有限個のこうした球の交わりが非空なら、二点を結ぶ大域的な最短弧は各球で同じ弧になり、交わりに留まる。共通点を固定して一つの球で定義した $\eta$ を使えば交わりも縮約できる。

この性質は、座標球の被覆で交わりが複数の成分に分かれる問題を避けるのに使える。ただし局所座標で Euclid 的に凸な開集合を選ぶだけでは、Riemann 距離に対する強凸性は保証されない。証明では距離の二乗の Hessian が正であることと、最短弧がその Hessian を制御できる場所に留まることの両方を用いた。

:::example[半球の交わり]
単位球面で $x_1>0$ と $x_3>0$ の開半球を取る。各半球内の二点は対蹠点でなく、$(1-t)x+ty$ の正規化で得る最短弧は正の座標成分を保つ。従って交わりでも最短弧が留まる。例えば $p=(1,0,1)/\sqrt2$ を固定すれば
$$
H(x,t)=\frac{(1-t)x+tp}{|(1-t)x+tp|}
$$
が交わりの縮約を与える。これは測地線の像を使うが、一般に $t$ は弧長に比例しない。それでも縮約には連続性と両端の値があれば十分である。
:::

## 距離関数と距離球

$r_p(q)=d_g(p,q)$ は三角不等式から $|r_p(q)-r_p(q')|\le d_g(q,q')$ を満たす。小さい正規球から $p$ を除けば、:ref[thm-rg-radial-minimal] により $r_p$ は正規座標の半径であり滑らかである。:ref[lem-rg-gauss] からその勾配は $\partial_r$ なので $\norm{\grad r_p}=1$ である。$r_p^2$ は正規座標の成分の二乗和なので $p$ でも滑らかである。$r_p$ 自身は非零方向の直線上で $|t|\norm v$ となるので、正の次元では $p$ で微分できない。

### 距離球面と境界

十分小さい $r>0$ で $r_p^{-1}(r)$ は正則値の定理により滑らかな超曲面である。外向き単位法線は $\grad r_p$ である。正規球の範囲を出ると、距離球面が複数の最短方向を接合することがある。[切断跡](/math/riemannian-geometry/index-cut-locus/) で微分可能性が失われる位置を述べる。

:::example[一次元の距離関数]
半径 $2$ の円で初期点との最小角差を $\theta\in[-\pi,\pi]$ とすれば $r_p=2|\theta|$ である。$\theta=\pi$ では同じ点の反対側で $\theta=-\pi$ に貼り合わさるため片側微分が異なる。
:::

:::example[平面の Hessian]
原点以外の Euclid 平面で $r=\sqrt{x^2+y^2}$ とする。$\partial_i r=x^i/r$、$\partial_j\partial_i r=\delta_{ij}/r-x^ix^j/r^3$ だから $\operatorname{Hess}r=(g-dr\otimes dr)/r$ である。$r=2$ では円に接する単位方向の値が $1/2$、動径方向は $0$ になる。
:::

### 二点の距離の微分

Whitehead の定理で選んだ近傍で $F(x,y)=d(x,y)^2/2$ とする。一意な初期ベクトル $v(x,y)=\exp_x^{-1}(y)$ は $(x,y)$ に滑らかに依存し、$F=\norm{v(x,y)}^2/2$ だから対角上でも滑らかである。端点の変分 $(X,Y)$ に沿って、定速最短弧 $\gamma$ を $[0,1]$ で変化させる。速度を $T$、変分場を $V$ とすると、捩率零により $D_sT=D_tV$ である。エネルギーの微分は $\int_0^1g(D_tV,T)\,dt$ となり、$D_tT=0$ なので境界項だけが残る。従って
$$
dF_{(x,y)}(X,Y)=-g_x(v(x,y),X)+g_y(\dot\gamma(1),Y).
$$
とくに、$y$ から $x$ への初速度は $-\dot\gamma(1)$ である。対角での一階微分は零であり、$y=x$ に固定したとき $x$ 方向の Hessian は $g_x$ となる。二点とも動く場合、対角での二次微分は $g_x(Y-X,Y-X)$ である。実際、逆関数定理で $dv_{(x,x)}(X,Y)=Y-X$、二乗ノルムの微分で残るのはこの線形項の二乗だけである。

:::example[二点を同時に動かす]
Euclid 平面で両点を $q=(0,0)$ から出し、$x(t)=(t,0)$、$y(t)=(0,2t)$ とする。$F(x(t),y(t))=5t^2/2$ なので二階微分は $5$ である。同じ初速度 $(1,0),(0,2)$ で一般の曲面上の二曲線を動かしても、基点で正規直交成分として与えれば二階微分は $5$ になる。曲率による違いを見分けるには、対角で二次より高い次数を見る必要がある。
:::

### Eikonal 方程式から動径を復元する

開集合で滑らかな関数 $r$ が $\norm{\grad r}=1$ を満たすとする。この式を :term[Eikonal 方程式]{en="eikonal equation" yomi="あいこなーるほうていしき"} という。任意の $X$ で $X\norm{\grad r}^2=2\operatorname{Hess}r(X,\grad r)=0$ である。Hessian の対称性より $g(\nabla_{\grad r}\grad r,X)=0$ がすべての $X$ で成立する。従って $\grad r$ の積分曲線は単位速測地線である。またその曲線に沿う $r$ の導関数は $\norm{\grad r}^2=1$ である。正規球から中心を除いた場合、この計算は距離が単位速動径に沿って増えることを別の方法で示す。ただし Eikonal 方程式だけから $r$ を一点からの距離と同定できない。例えば Euclid 平面の $r(x,y)=x$ は同じ式を満たし、水準集合は円でなく平行線である。

:::example[同じ高さの二点の初期速度と指数写像]
計量 $(dx^2+dy^2)/y^2$ の上半平面で $p=\left( -1, \  1\right)$、$q=\left( 1, \  1\right)$ とする。水平な線分は一般に最短でない。中点の高さを計算すると、最短弧がより高い位置を通ることが分かる。双曲面の Lorentz 内積を $B(V,W)=-V_0W_0+V_1W_1+V_2W_2$ とする。上半平面からの等長写像 :ref[prop-rg-half-space-isometry] を代入すると
$$
\begin{aligned}
P&=\begin{pmatrix}\frac{3}{2}\\-1\\\frac{1}{2}\end{pmatrix}\\
Q&=\begin{pmatrix}\frac{3}{2}\\1\\\frac{1}{2}\end{pmatrix}\\
C=-B(P,Q)&=3
\end{aligned}
$$
になる。両点は $B(P,P)=B(Q,Q)=-1$ を満たし、第 0 成分が正なので同じ双曲面上にある。

$Q-CP$ は $B(P,Q-CP)=-C+C=0$ なので接ベクトルである。その二乗長は
$$
B(Q-CP,Q-CP)=-1+2C^2-C^2=C^2-1=8
$$
である。従って単位初期方向は $U=(Q-CP)/\sqrt{C^2-1}$ であり、この例では $\begin{pmatrix}- \frac{3 \sqrt{2}}{4}\\\sqrt{2}\\- \frac{\sqrt{2}}{4}\end{pmatrix}$ となる。曲線
$$
\gamma(t)=\cosh t\,P+\sinh t\,U
$$
は $B(\gamma,\gamma)=-1$、$B(\dot\gamma,\dot\gamma)=1$ を満たす。さらに $\ddot\gamma=\gamma$ は周囲での法方向だから、接加速度は零である。$d=\log(C+\sqrt{C^2-1})$ とすると $C-\sqrt{C^2-1}=(C+\sqrt{C^2-1})^{-1}$ により $\cosh d=C$、$\sinh d=\sqrt{C^2-1}$ であり、$\gamma(d)=Q$ を得る。

この曲線は全体の最短弧である。実際、Lorentz 基底 $(P,U,\text{残りの単位接方向})$ を標準 Lorentz 基底に移すと、$P$ は頂点、$Q$ は動径距離 $d$ の点になる。極座標の速度は動径速度の絶対値以上であるから、どの経路も長さ $d$ 以上である。動径で等号が成立し、この例では $d=\log{\left(2 \sqrt{2} + 3 \right)}$ である。従って指定した初期ベクトルを指数写像に入れるなら $dU$、弧長 $t$ の初速度を指定するなら $U$ を用いる。

中点の計算には $P+Q=(1+\cosh d)P+\sinh d\,U$ を用いる。半角公式で右辺は $2\cosh(d/2)(\cosh(d/2)P+\sinh(d/2)U)$ になる。従って Lorentz 二乗長を $-1$ に正規化した
$$
M=\frac{P+Q}{\sqrt{2+2C}}=\begin{pmatrix}\frac{3 \sqrt{2}}{4}\\0\\\frac{\sqrt{2}}{4}\end{pmatrix}
$$
が中点である。上半平面へ戻す式 $y=1/(M_0-M_2)$、$x=M_1/(M_0-M_2)$ から中点は $\left( 0, \  \sqrt{2}\right)$ となる。上半平面では軌跡は中心 $c=0$、半径 $\sqrt{2}$ の半円である。中心は両端までの Euclid 二乗距離を等しいと置いた式から求めた。中点を代入しても $(x-c)^2+y^2=R^2$ となる。

同じ表示で一般の比率 $\lambda\in[0,1]$ の点は
$$
\gamma(\lambda d)
=\frac{\sinh((1-\lambda)d)}{\sinh d}P+
\frac{\sinh(\lambda d)}{\sinh d}Q
$$
である。$U=(Q-\cosh d\,P)/\sinh d$ を代入し、$\cosh(\lambda d)\sinh d-\sinh(\lambda d)\cosh d=\sinh((1-\lambda)d)$ を使った。二つの係数は通常の線形補間の $1-\lambda,\lambda$ とは異なる。これは $\gamma$ を双曲面上に保ちつつ弧長を比例させるためである。
:::

## 演習問題

::::exercise{level=basic}
$g=dx^2$ の直線で $\gamma(t)=2+3t^2$ は測地線か。

:::solution
$D_t\dot\gamma=\ddot\gamma\partial_x=6\partial_x$ であり、零でない。正則な部分の像は直線の区間だが、与えられたパラメータでは速度が一定でない。
:::
::::

::::exercise{level=basic}
半径 $2$ の球で $\norm v=\pi$ のとき $\exp_p(v)$ と $p$ の内積を求めよ。

:::solution
式 :eqref[eq-rg-sphere-exp] で $\norm v/a=\pi/2$ なので像は $2v/\pi$ となる。$v\perp p$ だから内積は $0$ である。
:::
::::

::::exercise{level=basic}
$g=dy^2/y^2$ の半直線で $\gamma(0)=1,\dot\gamma(0)=2$ の測地線を求めよ。

:::solution
$s=\log y$ で $g=ds^2$ となる。$s$ の初速度も $2$ なので $s(t)=2t$、すなわち $y(t)=e^{2t}$ である。
:::
::::

::::exercise{level=basic}
平面の原点からの距離の Hessian を $(3,4)$ において $(4/5,-3/5)$ に二度評価せよ。

:::solution
$r=5$ であり、指定したベクトルは単位かつ動径と直交する。「距離関数と距離球」の式により値は $1/r=1/5$ である。
:::
::::

::::exercise{level=advanced}
計量を正定数倍しても指数写像とその定義域が変わらないことを示せ。

:::solution
Koszul の公式 :eqref[eq-rg-koszul] では両辺が同じ定数倍になるので接続が変わらない。よって測地線の初期値問題が同じであり、最大存在区間と時刻 $1$ の値も一致する。ただし接空間の球の半径を数値で固定した場合は球自体が変わる。
:::
::::

::::exercise{level=advanced}
$dr^2+f(r)^2d\theta^2$ において $r=r_0$ の円が非定数測地線となる条件を求めよ。例として $f(r)=2+\cos r$ を調べよ。

:::solution
$\ddot\theta=0$ なので $\theta$ は非零の定速度で動く。動径方程式は $ff'\dot\theta^2=0$ なので $f'(r_0)=0$ が必要十分である。例では $-\sin r_0=0$、すなわち $r_0\in\pi\Z$ である。
:::
::::

::::exercise{level=advanced}
最短曲線が等速なら $[0,1]$ 上でエネルギー $\frac12\int_0^1\norm{\dot\gamma}^2$ も最小であることを示せ。

:::solution
任意の候補について Cauchy–Schwarz により $L^2\le\int_0^1\norm{\dot\gamma}^2$ である。等速の最短曲線では $L=d_g$ かつ Cauchy–Schwarz が等号なのでエネルギーは下界 $d_g^2/2$ を実現する。
:::
::::

::::exercise{level=advanced}
小さい正規球で $r_p^2/2$ の Hessian が中心で計量と一致することを証明せよ。

:::solution
正規座標で $r_p^2/2=\sum_i(x^i)^2/2$ である。中心では一次微分が零であり、Hessian の Christoffel 項も零となる。二階偏微分は $\delta_{ij}$、これは中心の計量成分に等しい。
:::
::::

## 確認テスト

::::quiz
測地線の定義を満たす条件はどれか。

- [x] $D_t\dot\gamma=0$
- [ ] どの座標でも $\ddot x^i=0$

:::explanation
共変加速度が零という条件である。成分には $\Gamma^i_{jk}\dot x^j\dot x^k$ が加わるので、極座標の直線のように $\ddot x^i$ が零でないこともある。
:::
::::

::::quiz{answer="3" tol="1e-6"}
Euclid 直線で $p=1,v=2$ のとき $\exp_p(v)$ を求めよ。

:::explanation
$\gamma(t)=1+2t$ が初期条件を満たすので時刻 $1$ の値は $3$ である。
:::
::::

::::quiz{answer="1" tol="1e-6"}
十分小さい正規距離球から中心を除いた部分における $\norm{\grad r_p}$ を求めよ。

:::explanation
Gauss の補題で $\partial_r$ は単位で角方向に直交する。ゆえに $\grad r_p=\partial_r$ であり長さは $1$ である。
:::
::::

::::quiz
どの主張が成り立つか。

- [x] 十分小さい距離球は強凸である。
- [ ] 任意の測地線弧は大域的に最短である。

:::explanation
:ref[thm-rg-whitehead] で両端点を動かしても十分小さい球内に一意の最短弧がある。円を長く回る測地線弧は大域的な最短性を失うので二番目は誤りである。
:::
::::

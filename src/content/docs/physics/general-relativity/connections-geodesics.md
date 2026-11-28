---
title: "共変微分と測地線"
description: "Levi-Civita 接続を構成し、平行移動、測地線、局所慣性座標と加速観測者を扱う。"
notation:
  - symbol: "\\nabla_a"
    meaning: "共変微分"
  - symbol: "\\Gamma^\\rho{}_{\\mu\\nu}"
    meaning: "座標基底での接続係数"
  - symbol: "u^a"
    meaning: "単位四元速度"
  - symbol: "a^a"
    meaning: "固有加速度"
  - symbol: "\\exp_p"
    meaning: "指数写像"
---

## この章の目標

- 座標での成分微分を共変微分に補正する。
- Levi-Civita 接続の一意性と測地線方程式を導く。
- 平行移動と固有時の変分の関係を説明する。
- 正規座標の有効範囲を述べ、加速観測者の枠を輸送する。

## 共変微分と接続係数

異なる点のベクトルは別々の接空間に属する。座標基底の変化を取り入れるため、微分の比較方法を指定する。[接続](/math/riemannian-geometry/connections/) の定義は計量の署名に依存しない。

:::definition[接続]{#def-gr-connection}
:term[接続]{en="connection" yomi="せつぞく"} $\nabla_XY$ は実双線形であり、関数 $f$ に対し $\nabla_{fX}Y=f\nabla_XY$、$\nabla_X(fY)=X(f)Y+f\nabla_XY$ を満たす。座標基底では $\nabla_{\partial_\mu}\partial_\nu=\Gamma^\rho{}_{\mu\nu}\partial_\rho$ と書く。
:::

余接ベクトルについては $\partial_\mu(\omega_\nu V^\nu)=(\nabla_\mu\omega)_\nu V^\nu+\omega_\nu\nabla_\mu V^\nu$ を要求する。任意の $V$ で比較すると

$$
\nabla_\mu V^\rho=\partial_\mu V^\rho+\Gamma^\rho{}_{\mu\sigma}V^\sigma,\qquad
\nabla_\mu\omega_\nu=\partial_\mu\omega_\nu-\Gamma^\sigma{}_{\mu\nu}\omega_\sigma .
$$

一般のテンソルにも積の微分則を課すと、上添字ごとに正、下添字ごとに負の接続項が加わる。新座標 $y$ で $\partial'_\alpha=K^\mu{}_\alpha\partial_\mu$ を代入して積を微分すると、

$$
\Gamma'^\gamma{}_{\alpha\beta}
=J^\gamma{}_\rho (K^\mu{}_\alpha K^\nu{}_\beta\Gamma^\rho{}_{\mu\nu}
+\partial'_\alpha K^\rho{}_\beta).
\label{eq-gr-connection-change}
$$

最後の項のため接続係数はテンソルの成分ではない。二つの接続の差ではこの項が消えるため、その差は $(1,2)$ 型テンソルになる。

:::example[一次元の座標]
$x>0$ の直線で旧接続を零とし、$y=x^2$ を用いる。$J=2x$、$K=1/(2x)$、$\partial_yK=-1/(4x^3)$ なので $\Gamma^y{}_{yy}=-1/(2y)$ である。旧座標で一定なベクトル $\partial_x$ の新成分 $2\sqrt y$ は、微分 $1/\sqrt y$ と接続項 $-1/\sqrt y$ が打ち消し合う。
:::

### 差テンソルと計量を使わない接続

$C(X,Y)=\widehat\nabla_XY-\nabla_XY$ は、接続の第二引数に関する積の法則を差し引くと $C(X,fY)=fC(X,Y)$ である。第一引数についても関数線形であるから、この差は各点の $X,Y$ だけで決まる。逆に任意の滑らかな $(1,2)$ 型テンソル $C$ に対し $\nabla+C$ は接続の公理を満たす。この構成に計量は要らない。例えば $\R^2$ で平坦接続に $C(X,Y)=\alpha(X)Y$、$\alpha=\dd x$ を加えると、$\widehat\nabla_{\partial_x}\partial_y=\partial_y$、$\widehat\nabla_{\partial_y}\partial_x=0$ になる。平行移動も元の接続の移動とは異なり、$x$ 方向に距離 $L$ 進んだときは成分が $e^{-L}$ 倍になる。実際、曲線 $x=s$ 上で $\widehat\nabla_{\partial_x}V=V'+V$ だからである。次節で課す二つの条件は、この差の自由を消す。

:::example[一形式の共変微分と Hessian]
$f$ の Hessian を $H(X,Y)=(\nabla_X\dd f)(Y)=X(Yf)-(\nabla_XY)f$ と定める。$H(X,Y)-H(Y,X)=-T(X,Y)f$ なので、ねじれが零なら対称である。平面の極座標で $f=r^2/2$ とすると $\partial_r f=r$、$\partial_\phi f=0$ である。従って $H_{rr}=1$、$H_{\phi\phi}=-\Gamma^r{}_{\phi\phi}r=r^2$ であり、$H_{r\phi}=0$ になる。この例では $H=g$ であり、直交座標で $f=(x^2+y^2)/2$ を直接微分した結果とも一致する。単なる二階偏微分では角成分を失う。
:::

## Levi-Civita 接続と Christoffel 記号

### 計量とねじれ

二つのベクトル場の Lie 括弧は、関数 $f$ への作用を $[X,Y](f)=X(Y(f))-Y(X(f))$ として定める。二階偏微分は差で消えるので、これはベクトル場である。テンソル $T(X,Y)=\nabla_XY-\nabla_YX-[X,Y]$ を :term[ねじれ]{en="torsion" yomi="ねじれ"} という。例えば $X\mapsto fX$ では $[fX,Y]=f[X,Y]-Y(f)X$ なので $T(fX,Y)=fT(X,Y)$ である。反対称性から第二引数でも同様である。座標基底では括弧が零なので、ねじれが零であることは $\Gamma^\rho{}_{\mu\nu}=\Gamma^\rho{}_{\nu\mu}$ と同値である。

:::theorem[Levi-Civita 接続]{#thm-gr-levi-civita}
非退化な滑らかな計量に対し、ねじれが零で $\nabla g=0$ を満たす :term[Levi-Civita 接続]{en="Levi-Civita connection" yomi="れゔぃちゔぃたせつぞく"} がただ一つ存在する。その :term[Christoffel 記号]{en="Christoffel symbols" yomi="くりすとっふぇるきごう"} は

$$
\Gamma^\rho{}_{\mu\nu}=\frac12g^{\rho\sigma}
(\partial_\mu g_{\sigma\nu}+\partial_\nu g_{\sigma\mu}-\partial_\sigma g_{\mu\nu}).
\label{eq-gr-christoffel}
$$

:::

:::proof
$\nabla_\mu g_{\nu\sigma}=0$ を $\partial_\mu g_{\nu\sigma}=g_{\rho\sigma}\Gamma^\rho{}_{\mu\nu}+g_{\nu\rho}\Gamma^\rho{}_{\mu\sigma}$ と書く。$\mu,\nu$ を交換した式を加え、$\sigma$ を微分添字とする式を引けば、下二添字の対称性により $2g_{\rho\sigma}\Gamma^\rho{}_{\mu\nu}$ が残る。これで一意性が分かる。逆に各座標内でこの式から接続を定めると、下二添字は対称で、右辺を計量微分に代入すると $\nabla g=0$ となる。重なりで座標変換した接続も同じ二つのテンソル条件を満たすため、一意性により一致する。
:::

:::example[平面の極座標]
$g=\dd r^2+r^2\dd\theta^2$ では $\partial_r g_{\theta\theta}=2r$ のみが非零である。よって $\Gamma^r{}_{\theta\theta}=-r$、$\Gamma^\theta{}_{r\theta}=\Gamma^\theta{}_{\theta r}=1/r$ である。平坦な計量でも係数は零とは限らない。
:::

:::example[静止時計の加速度]
$\dd s^2=-A(z)^2\dd t^2+\dd z^2$、$A>0$ では $\Gamma^z{}_{tt}=AA'$、$\Gamma^t{}_{zt}=A'/A$ である。静止観測者の四元速度は $u^t=1/A$ なので、$a^z=(\nabla_u u)^z=A'/A$ となる。
:::

:::example[球面上で計量整合性を確かめる]
単位二球面の $0<\theta<\pi$ で接続を作る。Lorentz 計量と正定値計量で同じ公式が使えるので、空間部分で符号を確かめられる。 座標の順序を $(\theta,\phi)$ と固定する。成分行列と逆行列は次の組になる。

$$
\begin{aligned}
(g_{ab})&=\left[\begin{matrix}1 & 0\\0 & \sin^{2}{\left(\theta \right)}\end{matrix}\right]\\
(g^{ab})&=\left[\begin{matrix}1 & 0\\0 & \frac{1}{\sin^{2}{\left(\theta \right)}}\end{matrix}\right]\\
\det g&=\sin^{2}{\left(\theta \right)}
\end{aligned}
$$

逆行列を掛けて単位行列を得ることを確かめてから、:eqref[eq-gr-christoffel] を使う。同じ下添字の順序を二度計算せず、$b\le c$ の組を選ぶ。中央の列には、座標に依存する計量成分の微分だけを残した。最後の列は微分を実行した値であり、表に出ない成分は添字交換で得られるものを除いて零である。

| 成分 | 微分への代入 | 結果 |
|---|---|---|
| $\Gamma^{\theta}_{\phi \phi}$ | $\frac{1}{2}\left(-\partial_{\theta}\left(\sin^{2}{\left(\theta \right)}\right)\right)$ | $- \frac{\sin{\left(2 \theta \right)}}{2}$ |
| $\Gamma^{\phi}_{\theta \phi}$ | $\frac{1}{2 \sin^{2}{\left(\theta \right)}}\left(\partial_{\theta}\left(\sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{1}{\tan{\left(\theta \right)}}$ |

アフィン径数の微分を点で表し、$w_a$ を座標 $a$ の速度成分とする。ここでは添字を下げる操作を行っていない。測地線の右辺 $-\Gamma^a{}_{bc}w_bw_c$ では $b\ne c$ の項を二度数える。

$$
\begin{aligned}
\ddot{\theta}&=\frac{w_{\phi}^{2} \sin{\left(2 \theta \right)}}{2}\\
\ddot{\phi}&=- \frac{2 w_{\phi} w_{\theta}}{\tan{\left(\theta \right)}}
\end{aligned}
$$

また計量ノルムは

$$
g(\dot x,\dot x)=w_{\phi}^{2} \sin^{2}{\left(\theta \right)} + w_{\theta}^{2}
$$
である。これを径数で微分して上の加速度を代入すると、係数の微分と速度の微分の和が零になる。従って初期のノルムが保存される。

$\nabla_\theta g_{\phi\phi}=\partial_\theta(\sin^2\theta)-2\Gamma^\phi{}_{\theta\phi}\sin^2\theta=0$ である。また $\nabla_\phi g_{\theta\phi}=-\Gamma^\phi{}_{\phi\theta}\sin^2\theta-\Gamma^\theta{}_{\phi\phi}=0$ である。後者は $g_{\theta\phi}=0$ でも二つの接続項がそれぞれ非零になる例である。
:::

:::example[円柱座標の接続と自由粒子]
Minkowski 計量の空間部で $X=r\cos\phi,Y=r\sin\phi$ とする。$\dd X^2+\dd Y^2=\dd r^2+r^2\dd\phi^2$ より、$r>0$ で以下の成分を得る。 座標の順序を $(t,r,\phi,z)$ と固定する。成分行列と逆行列は次の組になる。

$$
\begin{aligned}
(g_{ab})&=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & r^{2} & 0\\0 & 0 & 0 & 1\end{matrix}\right]\\
(g^{ab})&=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & \frac{1}{r^{2}} & 0\\0 & 0 & 0 & 1\end{matrix}\right]\\
\det g&=- r^{2}
\end{aligned}
$$

逆行列を掛けて単位行列を得ることを確かめてから、:eqref[eq-gr-christoffel] を使う。同じ下添字の順序を二度計算せず、$b\le c$ の組を選ぶ。中央の列には、座標に依存する計量成分の微分だけを残した。最後の列は微分を実行した値であり、表に出ない成分は添字交換で得られるものを除いて零である。

| 成分 | 微分への代入 | 結果 |
|---|---|---|
| $\Gamma^{r}_{\phi \phi}$ | $\frac{1}{2}\left(-\partial_{r}\left(r^{2}\right)\right)$ | $- r$ |
| $\Gamma^{\phi}_{r \phi}$ | $\frac{1}{2 r^{2}}\left(\partial_{r}\left(r^{2}\right)\right)$ | $\frac{1}{r}$ |

アフィン径数の微分を点で表し、$w_a$ を座標 $a$ の速度成分とする。ここでは添字を下げる操作を行っていない。測地線の右辺 $-\Gamma^a{}_{bc}w_bw_c$ では $b\ne c$ の項を二度数える。

$$
\begin{aligned}
\ddot{t}&=0\\
\ddot{r}&=r w_{\phi}^{2}\\
\ddot{\phi}&=- \frac{2 w_{\phi} w_{r}}{r}\\
\ddot{z}&=0
\end{aligned}
$$

また計量ノルムは

$$
g(\dot x,\dot x)=r^{2} w_{\phi}^{2} + w_{r}^{2} - w_{t}^{2} + w_{z}^{2}
$$
である。これを径数で微分して上の加速度を代入すると、係数の微分と速度の微分の和が零になる。従って初期のノルムが保存される。

さらに角方向の式に $r^2$ を掛けて積の微分を戻せば $(r^2\dot\phi)^{\cdot}=0$ になる。この計量は直交座標の平坦な計量から得たものである。角方向の接続係数が非零であっても、軌道と保存量を直交座標へ戻せば自由な直線運動を回復する。平坦性を調べる一般の方法は [曲率](/physics/general-relativity/curvature/) で扱う。

$r=2,\phi=0$ から $\dot r=0,\dot\phi=1$ で出る空間軌道は $X=2,Y=2s$ である。実際 $r(s)=2\sqrt{1+s^2}$、$\phi(s)=\arctan s$ を微分すると初期の $\ddot r$ は $2$ である。この非零の成分は直線運動を曲線座標で記録した結果である。
:::

### 座標基底を使わない Koszul の公式

$\nabla g=0$ なら $Xg(Y,Z)=g(\nabla_XY,Z)+g(Y,\nabla_XZ)$ である。これを $(X,Y,Z)$ で巡回させ、初めの二式を加え、$Zg(X,Y)$ の式を引く。ねじれが零なので、引数の順序を交換した接続の差を Lie 括弧に置き換えられる。残る式は
$$
\begin{aligned}
2g(\nabla_XY,Z)
={}&Xg(Y,Z)+Yg(Z,X)-Zg(X,Y)\\
&-g(X,[Y,Z])+g(Y,[Z,X])+g(Z,[X,Y]).
\end{aligned}
$$
これを :term[Koszul の公式]{en="Koszul formula" yomi="こしゅーるのこうしき"} と呼ぶ。座標枠のときは括弧が全て消え、:eqref[eq-gr-christoffel] に戻る。一方、正規直交枠では内積を微分する三項が消える。従って、どちらの種類の枠を選ぶかに応じて計算を短くできる。

:::example[球面の単位枠の接続を求める]
単位球で $e_1=\partial_\theta,e_2=\csc\theta\,\partial_\phi$ とする。$[e_1,e_2]=-\cot\theta\,e_2$ だから、Koszul の式に $X=Y=e_2,Z=e_1$ を入れると $2g(\nabla_{e_2}e_2,e_1)=-2\cot\theta$ となる。$e_2$ との内積は単位長の微分により零なので $\nabla_{e_2}e_2=-\cot\theta\,e_1$ である。直交性を微分すれば $\nabla_{e_2}e_1=\cot\theta\,e_2$ を得る。同じ公式を $X=e_1$ に適用すると $\nabla_{e_1}e_1=\nabla_{e_1}e_2=0$ となる。従って子午線方向に進む枠は平行だが、通常の緯線をたどる枠は回転する。
:::

### 共変発散と密度

接続の添字を縮約すると、発散を計量の行列式から計算できる。Christoffel 記号の式で $\Gamma^\mu{}_{\mu\nu}=\frac12g^{\mu\rho}\partial_\nu g_{\rho\mu}$ である。最初と最後の微分項は $\mu,\rho$ の交換により相殺する。行列式を一列ずつ微分すると、余因子による式 $\partial_\nu D=D g^{\mu\rho}\partial_\nu g_{\rho\mu}$ が出る。これを使うと $\Gamma^\mu{}_{\mu\nu}=\partial_\nu\log\sqrt{|D|}$ である。従って任意のベクトル場 $V$ に対し
$$
\nabla_\mu V^\mu=\frac{1}{\sqrt{|D|}}\partial_\mu(\sqrt{|D|}\,V^\mu)
$$
を得る。これは前章で体積形式から定めた発散と同じであり、計算には向き付けの選び方を使わない。

:::example[球面上の波動演算子の角部分]
半径 $r$ の球面で $f(\theta,\phi)$ の勾配は $r^{-2}f_\theta\partial_\theta+(r^2\sin^2\theta)^{-1}f_\phi\partial_\phi$ である。発散の式を適用すると
$$
\Delta_{S^2_r}f=\frac{1}{r^2}\left[
\frac{1}{\sin\theta}\partial_\theta(\sin\theta\,\partial_\theta f)
+\frac{1}{\sin^2\theta}\partial_\phi^2 f\right].
$$
$f=\cos\theta$ に対しては $-2\cos\theta/r^2$ を得る。$f=\cos\theta$ は球面の全域で滑らかなので、式の見かけ上の特異性は球面の極での座標の問題である。例えば北極に近づけると $\partial_\theta f=-\sin\theta$ が分母と相殺して有限な値に達する。
:::

## 平行移動

### 曲線に沿う比較

曲線 $\gamma(\lambda)$ に沿うベクトルを $V(\lambda)$ とし、
$(D_\lambda V)^\mu=\dot V^\mu+\Gamma^\mu{}_{\nu\rho}\dot x^\nu V^\rho$ とする。:term[平行移動]{en="parallel transport" yomi="へいこういどう"} は $D_\lambda V=0$ の解である。連続な係数をもつ線形常微分方程式なので、初期値を定めると曲線の全区間で一意に存在する。逆向きの曲線では方程式全体が再パラメータ化の微分因子を受けるだけである。そのため順方向と逆方向の移動は互いに逆となる。

並行な $V,W$ に対して $\dv{}{ \lambda}g(V,W)=g(D_\lambda V,W)+g(V,D_\lambda W)=0$ である。したがって移動は長さと内積を保つ。方程式で $D_\lambda V$ が零になるという条件は単調なパラメータ替えに依存しない。

:::example[極座標の円周]
平面の $r=R$ の円周で、規格化した成分を $(V^r,RV^\theta)=(A,B)$ とおく。接続係数から $A'=B$、$B'=-A$ となる。初期値 $(1,0)$ の解は $(\cos\theta,-\sin\theta)$ である。座標基底は回転するが、直交座標でベクトルは一定である。
:::

:::example[球面の緯線]
単位球の $g=\dd\theta^2+\sin^2\theta\,\dd\phi^2$ で $\Gamma^\theta{}_{\phi\phi}=-\sin\theta\cos\theta$、$\Gamma^\phi{}_{\phi\theta}=\cot\theta$ である。$\theta$ 一定の曲線上では $(V^\theta,\sin\theta\,V^\phi)$ が角度 $-\phi\cos\theta$ だけ回転する。$\theta=\pi/3$ なら一周後に角度が $-\pi$ だけ変わる。
:::

### 平行移動を行列の初期値問題として解く

曲線 $\gamma(s)$ に沿って $B^i{}_j(s)=\Gamma^i{}_{kj}\dot\gamma^k$ と置けば、平行移動は $V'=-BV$ である。初期行列 $U(s_0)=I$ に $U'=-BU$ を課す。各列が初期基底の平行移動を表し、線形性から $V(s)=U(s)V(s_0)$ である。行列 $W'=WB$、$W(s_0)=I$ も解くと $(WU)'=0$ なので $WU=I$ となり、移動は可逆である。係数が時刻に依存すると一般には $B(s_1)B(s_2)\ne B(s_2)B(s_1)$ である。そこで積分形 $U=I-\int B U$ に自分自身を代入して
$$
U=I-\int_{s_0}^s B(t_1)\dd t_1
+\int_{s_0\le t_2\le t_1\le s} B(t_1)B(t_2)\dd t_1\dd t_2-\cdots
$$
を得る。コンパクト区間で $\norm B\le C$ なら $n$ 次の項は $C^n(s-s_0)^n/n!$ 以下であり、一様収束する。積分形に戻すと方程式と初期条件を満たす。特にすべての $B(t)$ が互いに可換なら積分領域を対称化でき、通常の行列指数 $\exp[-\int B]$ が得られる。

:::example[球面の緯線に沿う一周]
極角 $\theta_0$ を一定とし、$e_1=\partial_\theta,e_2=(\sin\theta_0)^{-1}\partial_\phi$ という単位枠を取る。接続から $\nabla_\phi e_1=\cos\theta_0\,e_2$、$\nabla_\phi e_2=-\cos\theta_0\,e_1$ なので、枠成分を複素数 $w=V^1+iV^2$ と書くと $w'=-i\cos\theta_0\,w$ である。$w(2\pi)=e^{-2\pi i\cos\theta_0}w(0)$ となり、角度を $2\pi$ で同一視すると回転は $2\pi(1-\cos\theta_0)$ である。$\cos\theta_0=3/4$ では $\pi/2$ 回転し、最初に $e_1$ を選べば最後に $e_2$ が得られる。
:::

### 平行移動が経路に依存するかを調べる

点 $p$ から $q$ へ二つの道 $\gamma_1,\gamma_2$ を取り、その移動を $P_1,P_2$ とする。逆道の移動が逆写像なので、$P_1=P_2$ であることと、$\gamma_1$ の後に $\gamma_2$ を逆にたどる閉曲線で移動が恒等写像になることは同値である。各道の計量が同じというだけでは、端点で二つの移動の間に回転が残る可能性を除けない。

:::example[球面上の直角三角形を二つの道で通る]
単位球を $\R^3$ に置き、$A=(1,0,0)$、$B=(0,1,0)$、$C=(0,0,1)$ とする。それぞれを大円の四分の一弧で結ぶ。大円の単位接線は平行である。曲線の Euclid 加速度が半径方向なので、球面へ接射影した加速度が零になるからである。球面上の接続がこの射影で与えられることは、積の法則、計量整合性、ねじれ零を直接確認して一意性を使えばよい。大円を含む平面の固定単位法線も球面の接平面にあり、周囲で一定なので平行になる。

$A$ で初期ベクトルを $(0,1,0)$ とする。$A\to B$ では弧の接線であり、$B$ で $(-1,0,0)$ になる。$B\to C$ の弧の平面は $yz$ 平面なので、このベクトルは固定法線のままで $C$ に達する。一方、$A\to C$ を直接進むと初期ベクトルは $xz$ 平面の法線であり、$C$ でも $(0,1,0)$ になる。二つの結果の角度は $\pi/2$ である。球面では同じ初期ベクトルを同じ端点に運んでも、経路の違いが測れる。
:::

### 逆行列で共変微分を読み取る

曲線に沿う平行移動を $P_{s\to s_0}$ とすれば、曲線上のベクトル場 $V$ の微分は
$$
D_s V(s_0)=\lim_{h\to0}\frac{P_{s_0+h\to s_0}V(s_0+h)-V(s_0)}{h}
$$
で得られる。実際、行列 $U(s_0+h)=I-hB(s_0)+o(h)$ の逆は $I+hB(s_0)+o(h)$ なので、差商は $V'(s_0)+B(s_0)V(s_0)$ に収束する。これにより接続の係数を直接測る代わりに、近くへ運んだ二本を一つの接空間で比較するという意味が明確になる。

## 測地線方程式と変分原理

### エネルギー汎関数

:term[測地線]{en="geodesic" yomi="そくちせん"} は $D_\lambda\dot\gamma=0$ を満たす曲線である。この $\lambda$ をアフィンパラメータという。計量整合性により $g(\dot\gamma,\dot\gamma)$ は一定なので、定数曲線でない測地線の因果型は途中で変わらない。

:::proposition[測地線の変分]{#prop-gr-geodesic-variation}
端点とパラメータの区間を固定した $\mathcal E[\gamma]=\tfrac12\int g_{\mu\nu}\dot x^\mu\dot x^\nu\dd\lambda$ の停留曲線は、アフィンパラメータで表した測地線である。
:::

:::proof
変分を $\delta x^\rho=\xi^\rho$ とおく。$g$ の変化は $\partial_\rho g_{\mu\nu}\xi^\rho$、速度の変化は $\dot\xi^\mu$ である。後者の項 $g_{\rho\nu}\dot x^\nu\dot\xi^\rho$ を部分積分すると端点項は零であり、

$$
\delta\mathcal E=\int\left[\frac12\partial_\rho g_{\mu\nu}\dot x^\mu\dot x^\nu-
\dv{}{\lambda}(g_{\rho\nu}\dot x^\nu)\right]\xi^\rho\,\dd\lambda
=-\int g_{\rho\sigma}(\ddot x^\sigma+\Gamma^\sigma{}_{\mu\nu}\dot x^\mu\dot x^\nu)\xi^\rho\,\dd\lambda .
$$

任意の内部変分を取れるため係数は零であり、計量の逆を掛けると測地線方程式を得る。逆に測地線では積分が零である。
:::

$\lambda=as+b$ で $a$ が定数なら測地線方程式は $a^2$ 倍になる。一般の $\lambda(s)$ では $D_s x'=(\lambda''/\lambda')x'$ である。したがって任意の座標時刻をパラメータに使う場合は右辺に接線方向の項を残す。

:::example[角運動量の保存]
平面極座標のラグランジアンに $\theta$ が陽に現れないので、変分方程式は $(r^2\dot\theta)^{\cdot}=0$ を与える。直線 $(x,y)=(s,2)$ では $r^2=s^2+4$、$\dot\theta=-2/(s^2+4)$ だから保存量は $-2$ である。
:::

### 任意の径数で作用を書く

非光的な曲線について、長さの被積分関数 $\sqrt{|g(\dot x,\dot x)|}$ は速度の一次同次関数である。このため単調に径数を変えても積分の値は変わらない。他方 $\mathcal E$ の被積分関数は二次同次なので、端点の位置を保つ径数の変更でも一般には値が変わる。固有時の停留条件からアフィン測地線を得るには、変分を求めた後に単位速度を径数に使う。例えば時刻 $s$ で $\ell(s)=\sqrt{-g(\dot x,\dot x)}>0$ なら $u=\dot x/\ell$ である。$\nabla_u u=0$ を $s$ で書くと $\nabla_{\dot x}\dot x=(\ell'/\ell)\dot x$ となる。この接線項は、四元加速度が零の運動を非一様な径数で記録した結果である。

:::example[時刻の目盛りを変えたときの二つの作用]
Minkowski 平面で $t=s^2,x=0$、$1\le s\le2$ とする。この世界線の固有時は $\int_1^2 2s\,\dd s=3$ である。一方 $\mathcal E=-\tfrac12\int_1^2(2s)^2\dd s=-14/3$ になる。同じ径数区間を用い、$t=3s-2$ と直せば端点は同じで固有時も $3$ だが、$\mathcal E=-9/2$ になる。第一の表示の共変加速度は接線の $1/s$ 倍であり、第二は共変加速度が零の表示である。長さとエネルギーのどちらを変分しているかを指定すると、求めた方程式で径数をどう選ぶべきかが分かる。
:::

### 質量ゼロの粒子に使う補助変数

光的曲線では長さ汎関数の被積分関数が零なので、その平方根を微分して運動方程式を求める手続きは使えない。そこで正の補助変数 $e(s)$ を置き、
$$
S[x,e]=\frac12\int(e^{-1}g_{ab}\dot x^a\dot x^b-e\,m^2)\dd s
$$
を変分する。$e$ の変分で $g(\dot x,\dot x)=-e^2m^2$ を得る。$x$ の変分では端点で変分を零とし、$\dv{}{s}(e^{-1}g_{ab}\dot x^b)$ と $e^{-1}\partial_ag_{bc}\dot x^b\dot x^c/2$ を等置する。逆計量で縮約すると $\nabla_{\dot x}\dot x=(\dot e/e)\dot x$ である。$\dd\lambda=e\,\dd s$ によって $e$ を吸収すると $\nabla_{x'}x'=0$、$g(x',x')=-m^2$ となる。従って $m=0$ でもアフィン径数の光的測地線を導ける。一般の径数による余分な加速度は接線方向に限られる。

:::example[同じ零直線を二つの径数で書く]
Minkowski 平面で $t=x=e^s$ を取る。速度のノルムは零で、各座標の二階微分は速度と同じである。従って右辺の係数は $1$ であり、$e=e^s$ と置くと補助変数の式に一致する。$\lambda=e^s$ に替えると座標はともに $\lambda$、加速度はともに零になる。光的条件は最初から成立しているため、光的条件だけではアフィン性を判定できない。
:::

### Hamilton 形式で座標運動量を得る

正則な二次ラグランジアン $L=\frac12g_{\mu\nu}\dot x^\mu\dot x^\nu$ に対し、共役運動量は $p_\mu=\partial L/\partial \dot x^\mu=g_{\mu\nu}\dot x^\nu$ である。計量の非退化性によりこの関係は $\dot x^\mu=g^{\mu\nu}p_\nu$ と解ける。Legendre 変換 $H=p_\mu\dot x^\mu-L$ から
$$
H=\frac12g^{\mu\nu}p_\mu p_\nu,\qquad
\dot x^\mu=g^{\mu\nu}p_\nu,\qquad
\dot p_\mu=-\frac12(\partial_\mu g^{\nu\rho})p_\nu p_\rho
$$
を得る。最後の式が Euler–Lagrange 方程式と一致することは、逆行列の微分 $\partial_\mu g^{\nu\rho}=-g^{\nu a}g^{\rho b}\partial_\mu g_{ab}$ を代入すれば分かる。$H$ は位相空間の関数であり、$\dot H=\partial_xH\,\partial_p H-\partial_p H\,\partial_x H=0$ である。座標 $x^\mu$ が計量に陽に現れなければ $p_\mu$ も定数となる。光的測地線では $H=0$ を初期値に課す。

:::example[直線の曲座標表示を保存量から積分する]
Euclid 平面で単位速度とし、$p_\phi=b>0$ とする。$H=1/2$ から $\dot r^2=1-b^2/r^2$ が従う。最近点を $s=0$、$r=b$ として、外向きでは $s=\int_b^r r'\dd r'/\sqrt{r'^2-b^2}=\sqrt{r^2-b^2}$ である。内向きの脚も含めると $r^2=b^2+s^2$ である。$\dot\phi=b/(b^2+s^2)$ を積分して $\phi=\phi_0+\arctan(s/b)$ となる。$\phi_0=0$ なら直交座標は $(X,Y)=(b,s)$ である。転回点は測地線の停止ではなく、動径だけが停留する地点である。
:::

### 速度を変えて同じ像をたどる

非定数解 $\gamma$ の接線は零にならない。もし途中で零になれば、その位置に留まる定数解と初期値が一致し、局所一意性から同じ解になるからである。このため同じ有向の測地線を二つの滑らかなパラメータで表し、どちらもアフィンとすると、再パラメータ化の式から $\lambda''=0$ になる。時間的な場合には単位ノルムを課すと傾きは $1$ に固定される。光的な場合はノルム零なので、この規格化で傾きを決められない。光子の波数として使うときは、ある観測者が測る周波数を指定して傾きを定められる。

:::example[時間座標で書いた光路には接線項が残る]
$\dd s^2=-A(z)^2\dd t^2+\dd z^2$、$A>0$ で上向きの光路は $\dd z/\dd t=A(z)$ を満たす。$A=1+az$ なら $z(t)=(e^{at}-1)/a$ が $z(0)=0$ の解である。アフィンな測地線では $E=A^2\dot t$ が一定で、光的条件から $\dot z=E/A$ である。$t$ をパラメータにした式では $t$ 成分の加速度は $2A'(\dd z/\dd t)/A=2a$ となり、これは共変加速度が接線に比例する係数である。従って $z$ の二階座標微分が $aA$ でも自由な光路と両立する。$a\to0$ の極限では $z=t$ に戻り、$t$ 自身もアフィンに取れる。
:::

## 固有時の極値性と自由落下

### 固有時と自由粒子

$m>0$ の構造を無視できる自由粒子に作用 $-m\int\dd\tau$、$\dd\tau=\sqrt{-g(\dot x,\dot x)}\dd\lambda$ を用いる。四元速度 $u=\dd x/\dd\tau$ で $g(u,u)=-1$ である。$\xi$ による変分では、ねじれが零なので $D_\xi\dot x=D_\lambda\xi$ である。平方根を微分して固有時に変数を替えると、端点固定で

$$
\delta\tau=-\int g(u,D_\tau\xi)\dd\tau=\int g(\nabla_u u,\xi)\dd\tau .
$$

よって $u$ について $\nabla_u u=0$ を得る。光では固有時が零なのでこの作用をそのまま使えず、アフィンなエネルギー汎関数の停留条件に零ノルムを課す。

### 局所最大と長い区間

Minkowski 時空では端点間の静止系で、到着までの座標時間を $T$ とすれば $\tau=\int_0^T\sqrt{1-\abs{\vb v}^2}\dd t\le T$ である。等号は静止曲線でのみ成り立つ。一般の時空でも十分短い時間的測地線は近くの曲線の中で固有時を最大にする。次節の指数座標で基点 $p$ からの未来時間的な領域に $\rho=\sqrt{-g_p(v,v)}$ を取り、$v=\rho n$ と書く。放射測地線の単位接ベクトルを $U=\partial_\rho$ とする。方向を変える変分 $W$ について $[U,W]=0$ なので、

$$
\partial_\rho g(U,W)=g(U,\nabla_W U)=\tfrac12 W(g(U,U))=0.
$$

基点で $W=0$ だから $g(U,W)=0$ であり、$g=-\dd\rho^2+h_\rho$、$h_\rho$ は正定値となる。その結果 $\dd\tau\le\dd\rho$ で、等号には方向が一定であることが必要である。この比較は放射座標が定まる近傍内の時間的曲線に適用する。長い区間では曲線の族が再交差し、局所的な最大性も失われ得る。[因果構造と特異点定理](/physics/general-relativity/causal-structure/) でその判定を扱う。

:::example[往復の時計]
慣性系で座標時間 $10\,\mathrm s$ の前半を速度 $3c/5$、後半を逆向きに同じ速度で動く。折り返しの時間を無視すると固有時は $10\sqrt{1-9/25}=8\,\mathrm s$ である。直進する比較時計は $10\,\mathrm s$ 進む。
:::

### 時間的軌道の局所最大を不等式で示す

Minkowski 時空で同じ始点と終点を結ぶ未来時間的曲線を比較する。Lorentz 変換で両端の空間位置を同じにすれば、座標時刻の差を $T>0$ として
$$
\tau=\int_0^T\sqrt{1-\abs{\vb v(t)}^2}\dd t\le T
$$
である。等号は速度が恒等的に零のときである。曲がった時空でも十分小さい正規近傍では、始点からの時間的放射測地線の固有時を $\rho$ として同様の比較をする。$\gamma_\epsilon(s)=\exp_p(s(v+\epsilon w))$ の変分場を $J$ とする。測地線の速度 $U$ と $J$ は交換し、$\nabla_UJ=\nabla_JU$ なので
$$
\frac{\dd}{\dd s}g(U,J)=\tfrac12\partial_\epsilon g(U,U)=g_p(v,w).
$$
$J(0)=0$ より $g(U,J)=s g_p(v,w)$ となる。初期速度を単位双曲面上で動かせば $g_p(v,w)=0$ であり、放射方向は $\rho=\mathrm{const}$ 面に直交する。こうして計量は時間的円錐内で $-\dd\rho^2+h_\rho$ と書け、$h_\rho$ は各切片で正定値である。未来時間的曲線には $\dd\rho>0$ で、被積分関数 $\sqrt{\dot\rho^2-h_\rho(\dot y,\dot y)}$ は $\dot\rho$ 以下である。始点で $\rho\to0$ の極限を取れば、正規近傍内の放射測地線が最大になる。より長い区間では、放射座標が使える範囲を別に調べる必要がある。

:::example[折れ曲がった軌道で時計の差を計算する]
座標時刻の半分ずつを速度 $3/5$ と $-3/5$ で進むと、始点と終点の空間位置は等しい。$T=10\,\mathrm s$ なら各脚の固有時は $5\sqrt{1-9/25}=4\,\mathrm s$、和は $8\,\mathrm s$ で、静止軌道の $10\,\mathrm s$ より短い。折れ曲がりを滑らかにした場合も、速度が非零の区間があるため不等号は厳しい。
:::

### 正規近傍で時間的な円錐から出られないこと

先ほどの最大化で用いた放射座標の範囲を確かめる。基点の十分小さい正規近傍を取り、そこで正規座標の時間成分は未来時間的な曲線に沿って増えるようにする。これは基点で時間成分の微分が時間的であり、その性質が近傍でも保たれるためである。この近傍内で基準点 $p$ から見た $q=\exp_p(v)$ に対して $\sigma(q)=g_p(v,v)/2$ とする。Gauss の補題から $\grad\sigma$ はアフィン時刻 $1$ での放射測地線の接線である。実際、速度 $v$ を $v+\epsilon w$ に変えた測地線族の変分ベクトル $J$ について、終点で $\dd\sigma(J)=g_p(v,w)=g(U,J)$ である。未来光円錐上で $q\ne p$ なら $\grad\sigma$ は未来光的で、未来時間的な接線 $T$ に対し $\dd\sigma(T)=g(\grad\sigma,T)<0$ となる。$p$ から未来時間的に出る滑らかな曲線は、初速が未来時間的なので初めは正規座標の未来時間的円錐に入る。この後、円錐を出る最初の点では $\sigma$ が負から零に到達するため左からの微分が非負である。一方、その境界で未来時間的な接線に対して上の微分が負であり、矛盾する。従って曲線が近傍にある限り、放射固有時 $\rho=\sqrt{-2\sigma}$ を使える。区分的に滑らかな時間的曲線でも境界直前の左接線で同じ議論を使う。これにより、近傍の中の全ての時間的比較曲線へ固有時不等式を適用できる。

### 時間的な長さの逆三角不等式

同じ点の二つの未来時間的ベクトル $u,v$ を取る。$u$ の静止枠で $u=(\ell,0,0,0)$ と書けば $-g(u,v)=\ell v^0\ge \ell\sqrt{-g(v,v)}$ である。等号は $v$ の空間成分が零のときに限る。これを $(u+v)^2$ に使えば
$$
\sqrt{-g(u+v,u+v)}\ge \sqrt{-g(u,u)}+\sqrt{-g(v,v)}
$$
を得る。Minkowski 時空で折れ線の端点を直接結ぶ方が固有時が長いという命題は、変位ベクトルを足してこの式を繰り返せば従う。連続な時間的曲線についても、区分直線で近似して極限を取るか、端点の静止枠で被積分関数を直接比較できる。

:::example[空間が周期的な平坦時空で長い区間を比べる]
$\dd s^2=-\dd t^2+\dd x^2$ に $x\sim x+\ell$ を課す。$p=(0,0)$ と $q=(T,0)$ を結ぶ測地線を被覆空間に持ち上げると、終点が $(T,n\ell)$ の直線となる。$|n|\ell<T$ なら時間的であり、固有時は $\sqrt{T^2-n^2\ell^2}$ である。例えば $T=5,\ell=3$ では $n=0$ の固有時は $5$、$n=\pm1$ では $4$ となる。巻き数を固定すると各直線は、その巻き数での時間的曲線の中で最大である。しかし同じ端点だけを固定すると、巻き数の異なる曲線と比較する必要がある。局所的な測地線方程式の一意性は、離れた二端点を結ぶ測地線の一意性を意味しない。
:::

### 自由運動と積分した固有時の二次変分

Minkowski 時空で同じ二端点を結ぶ静止直線から、空間方向だけを $\boldsymbol x(t)=\epsilon\boldsymbol\xi(t)$ と変分する。$\boldsymbol\xi(0)=\boldsymbol\xi(T)=0$ とし、比較曲線は未来時間的となるほど $\epsilon$ を小さく取る。すると
$$
\tau(\epsilon)=T-\frac{\epsilon^2}{2}\int_0^T|\dot{\boldsymbol\xi}|^2\dd t+\order{\epsilon^4}.
$$
従って第一変分は零、第二変分は $\boldsymbol\xi$ が恒等的に零の場合を除いて負である。例えば $\boldsymbol\xi(t)=B\sin(\pi t/T)\boldsymbol e_1$ なら、固有時の減少は $\epsilon^2 B^2\pi^2/(4T)$ が先頭になる。一次では時間的測地線と近くの曲線を区別できないが、二次では変分の速度が寄与する。曲がった時空ではこの二次式に曲率が加わるため、長区間での符号には追加の条件を要する。[共役点と焦点](/physics/general-relativity/causal-structure/) でその条件を扱う。

:::example[力があると固有時の第一変分を検出できる]
世界線が単位速度 $u$ と加速度 $a=\nabla_u u$ を持つとする。端点を固定した変分場を $\xi=\chi(\tau)a$ とし、非負の滑らかな $\chi$ の台を $a$ が非零の小区間に取る。変分場は指数写像で曲線を少し移すことにより実現でき、変分が小さければ曲線は時間的なままである。第一変分は $\delta\tau=\int\chi\,g(a,a)\dd\tau>0$ になる。加速中の世界線の近くには、同じ端点を結び固有時がより長い曲線が存在することが分かる。この検査では座標の二階微分ではなく、四元加速度を比較する。
:::

## Riemann 正規座標

:term[指数写像]{en="exponential map" yomi="しすうしゃぞう"} $\exp_p(v)$ は初期値 $(p,v)$ の測地線をパラメータ $1$ まで進めた点である。速度を $sv$ に替えた解は元の解をパラメータ $s$ まで進めたものなので、$\frac{\dd}{\dd s}\exp_p(sv)|_{s=0}=v$ である。初期値への滑らかな依存性と逆関数定理から、原点近傍で $\exp_p$ は微分同相となる。$T_pM$ で正規直交基底を選び、$v$ の成分を $\exp_p(v)$ の座標に取る。これを :term[Riemann 正規座標]{en="Riemann normal coordinates" yomi="りーまんせいきざひょう"} という。

以下では、基点からの指数写像が原点を含む星形領域上の微分同相になる近傍を正規近傍と呼ぶ。この座標で放射測地線は $x^\mu=s v^\mu$ であり、$\Gamma^\mu{}_{\nu\rho}(sv)v^\nu v^\rho=0$ となる。原点ではすべての $v$ で成立するので、下添字の対称性と偏極により $\Gamma^\mu{}_{\nu\rho}(0)=0$ である。計量整合性から $\partial_\rho g_{\mu\nu}(0)=0$ となり、:ref[def-gr-local-inertial]{full} を幾何学的に構成できた。原点での二次の変化には次章の曲率が入る。

:::example[慣性座標と角度座標]
平面の $p=(2,0)$ から初速度 $v=(1,1)$ の測地線をアフィンパラメータ $s$ で表すと、直交座標では $(2+s,s)$ である。$s=1$ の点は極座標で $(\sqrt{10},\arctan(1/3))$ だが、$p$ を基点とする正規座標では $(1,1)$ である。
:::

### 放射測地線の三次までの展開

任意の滑らかな座標 $y^a$ の原点から初速度 $v$ で出る測地線について、$\dot y(0)=v$ とする。測地線方程式をそのまま用いると $\ddot y^a(0)=-\Gamma^a{}_{bc}(0)v^bv^c$ である。式を一度微分し、二つの速度の微分を接続の対称性でまとめれば
$$
y^{a(3)}(0)=
\bigl[-\partial_d\Gamma^a{}_{bc}
+2\Gamma^a{}_{bk}\Gamma^k{}_{cd}\bigr]_0 v^bv^cv^d .
$$
従って初速度を座標に採用する正規座標から元の座標への変換は、一階では同じ基底、二階では $\Gamma$ の補正、三階ではその微分と積を含む。二階の補正だけで一点の接続を消せるが、同じ座標で全近傍の接続が消えるかは曲率に依存する。

:::example[極座標で直線の初速度を座標にする]
点 $(X,Y)=(2,0)$ を基準に、Cartesian 成分 $(p,q)$ の初速度を取る。$X=2+sp,Y=sq$ だから $r(s)=\sqrt{(2+sp)^2+s^2q^2}$ である。展開すると $r=2+ps+q^2s^2/4-pq^2s^3/8+\order{s^4}$ となる。極座標での初期角速度は $q/2$ で、$\Gamma^r{}_{\phi\phi}=-2$ から $\ddot r(0)=q^2/2$ となり、展開の二次項と一致する。
:::

### Lorentz 署名にも通用する Gauss の補題

$v,w\in T_pM$ を取り、$\gamma_\epsilon(s)=\exp_p(s(v+\epsilon w))$ とする。$V=\partial_s\gamma$、$W=\partial_\epsilon\gamma$ と書く。ねじれ零より $D_sW=D_\epsilon V$ であり、各曲線で $D_sV=0$ なので
$$
\partial_s g(V,W)=\tfrac12\partial_\epsilon g(V,V)
=g_p(v+\epsilon w,w).
$$
$s=0$ では $W=0$ である。$\epsilon=0,s=1$ で積分した結果が
$$
g_{\exp_p(v)}\big((d\exp_p)_v(v),(d\exp_p)_v(w)\big)=g_p(v,w)
$$
となる。$\exp_p(sv)$ を $s$ で微分すると第一因子、$v$ を $w$ 方向へ変分すると第二因子が出る。この証明には正定値性を用いていないので、光的な $v$ も含む。例えば $g_p(v,v)=0$ のとき、放射接線は光円錐の接空間に属し、その全ての接線と直交する。固有時の局所最大性を示すときに用いた $\grad\sigma$ の同定も、この式の右辺 $g_p(v,w)$ から従う。

:::example[球面の正規座標が大域で破綻する]
単位球で北極から極角 $\theta=s$ の子午線を進む。Gauss の補題による放射計量は $\dd s^2+\sin^2s\dd\phi^2$ であり、同じ $s$ で初期方向を変えた速度の長さは $\sin s$ である。$s\to0$ で角方向の長さは $s$ に比例するため、初期速度から作る直交型の正規座標は原点まで滑らかである。他方 $s=\pi$ では全子午線が南極に着き、初期方向を変えた終点の微分が零になる。初期値問題の解はそこでも続くが、指数写像はその速度で局所微分同相ではない。
:::

### 正規座標同士を比較する

同じ点 $p$ で二つの正規直交基底 $e_A,e'_A$ を選ぶ。両者の正規座標は、同じ $v\in T_pM$ の二組の成分だから、重なりでは一定の Lorentz 変換で結ばれる。一方、別の基点で作ると、その指数写像の基点と初期値が変わる。座標変換の二階微分は両座標での接続の差を含み、一般には一定の Lorentz 変換にできない。各点での局所慣性枠が存在することと、時空全体に同じ慣性枠を選べることは区別する必要がある。

:::example[球面の二つの基点を使う]
単位球上で点 $p$ と、その子午線上で $p$ から距離 $\pi/4$ の点 $q$ を取る。子午線の接線に沿う初期方向と垂直な方向の変分を考える。$p$ の正規座標で、曲線の始点の方向を変えず半径 $s=\pi/4$ だけ進んだ終点が $q$ である。その半径を固定し初期方向を短い角度 $\dd\phi$ だけ変えると、$p$ の正規座標での横変位は $s\dd\phi$、実際の固有長は $\sin s\,\dd\phi$ である。$q$ の正規座標ではこの固有長が横成分の一次変化に等しい。従って二つの正規座標の変換の横方向の倍率は $\sin(\pi/4)/(\pi/4)$ である。正規座標が原点で正規直交になるという条件だけでは、他の点の係数を保てない。
:::

## Fermi–Walker 輸送と加速度観測者

加速観測者では $a=\nabla_u u$ が零でない。$u\cdot u=-1$ を微分すると $u\cdot a=0$ である。観測者とともに進み固有回転を加えない枠は次の式で定める。

:::definition[Fermi–Walker 輸送]{#def-gr-fermi-walker}
ベクトル $V$ の :term[Fermi–Walker 輸送]{en="Fermi–Walker transport" yomi="ふぇるみうぉーかーゆそう"} は

$$
\nabla_u V^a=(u^aa_b-a^au_b)V^b
\label{eq-gr-fermi-walker}
$$

を満たす輸送である。
:::

右辺の演算子は下添字に下ろすと反対称なので、二つの輸送ベクトルの内積は一定である。$V=u$ とおくと右辺は $a$ になり、$u$ 自身を輸送する。$V\cdot u=0$ なら右辺は $u(a\cdot V)$ であり、観測者の空間への射影は零である。したがって空間三脚には観測者が感じる回転がない。逆に空間への射影が零で $V\cdot u=0$ を維持するなら、$\nabla_u V=b u$ とおいて $0=-b+V\cdot a$ から同じ式が得られる。

:::example[一様固有加速度]
Minkowski 平面で $u=(\cosh(a_0\tau),\sinh(a_0\tau))$ とする。$e=(\sinh(a_0\tau),\cosh(a_0\tau))$ に対し $u'=a_0e$、$e'=a_0u$ である。$g(a,e)=a_0$ なので $e$ は Fermi–Walker 輸送される。$a_0\tau=\log 2$ では $e=(3/4,5/4)$ であり、$g(e,e)=1$ を保つ。
:::

### 加速枠に付随する座標

Minkowski 時空の観測者 $z^a(\tau)$ に沿って、$u=\dot z$ と :eqref[eq-gr-fermi-walker] で輸送する正規空間枠 $e_i$ を取る。近くの事象を $x^a=z^a(\tau)+X^ie_i{}^a(\tau)$ と書くと、$e_i'=a_i u$ により $\partial_\tau x=(1+a_iX^i)u$、$\partial_{X^i}x=e_i$ である。従って

$$
\dd s^2=-(1+a_i(\tau)X^i)^2\dd\tau^2+\delta_{ij}\dd X^i\dd X^j .
$$

これは観測者の瞬間静止面を貼り合わせた局所座標である。$1+a_iX^i=0$ に達すると Jacobian が退化し、加速枠を保ったままこの範囲を越えられない。一般の曲がった時空では $x=\exp_{z(\tau)}(X^ie_i)$ を使い、中心線上の同じ一次項に加えて $X$ の二次以上に曲率が入る。

:::example[同時面での時計比較]
$a_1=a_0$ が一定、他が零とする。$X^1=1/(2a_0)$ に静止する時計は中心の $\Delta\tau=2$ に対し固有時 $3$ を測る。この時計自身の固有加速度を求めるには $\Gamma^1{}_{\tau\tau}=(1+a_0 X^1)a_0$ に $(u^\tau)^2=(1+a_0 X^1)^{-2}$ を掛ける。結果は $a_0/(1+a_0X^1)=2a_0/3$ である。座標網で一定の間隔を保つには、すべての時計が同じ固有加速度を持つわけにはいかない。
:::

:::example[曲線に沿う再パラメータ化]
$D_s x'=F(s)x'$ を満たす曲線からアフィンパラメータ $\lambda$ を作るには $\lambda''/\lambda'=F(s)$ を解く。例えば $F(s)=2/s$、$s>0$ なら $\lambda'=Cs^2$、$\lambda=Cs^3/3+D$ である。$\lambda$ の二つの積分定数はアフィンパラメータの原点と尺度を表す。$\lambda'$ の符号で進行方向を選べる。
:::

### 回転する実験室と局所の慣性軸

Minkowski 時空の回転座標 $\phi=\phi_{\rm in}-\Omega t$ では
$$
\dd s^2=-(1-\Omega^2r^2)\dd t^2+2\Omega r^2\dd t\,\dd\phi+\dd r^2+r^2\dd\phi^2+\dd z^2 .
$$
一定 $r,\phi,z$ の観測者は $|\Omega|r<1$ で未来時間的で、$u^t=\gamma=(1-\Omega^2r^2)^{-1/2}$ である。慣性座標で $X=r\cos(\Omega t),Y=r\sin(\Omega t)$ を二回固有時微分すると空間加速度は内向きで、その大きさは $\gamma^2\Omega^2r$ になる。加速があると、平行輸送した軸は一般には観測者の静止空間から外れる。実際、平行な軸 $V$ に対して $g(u,V)$ の微分は $g(a,V)$ である。例えば $u\cdot V=0$ なら Fermi–Walker 輸送は $\nabla_u V=(a\cdot V)u$ であり、加速度の向きの変化を含めて直交を保つ。観測者の体に固定された枠には、さらに静止空間内での回転が加わり得る。

:::example[静止空間の二本の軸の内積]
Fermi–Walker 輸送する二本の軸を $V,W$ とすると、$V\cdot u=W\cdot u=0$ のとき $\dv{}{\tau}(V\cdot W)=(a\cdot V)u\cdot W+(a\cdot W)u\cdot V=0$ である。例えば最初の角度が $\pi/3$、双方のノルムが $1$ なら、円運動のどの位置でも内積は $1/2$ である。座標成分が変わっても、時計と共に運ぶ軸が示す内積は固定される。
:::

:::example[回転する座標で接続と回転項を計算する]
$\phi_{\rm in}=\phi+\Omega t$ である回転座標を使い、$r>0$ とする。ここで $\Omega$ は定数である。観測者が座標に静止できるのは $|\Omega|r<1$ であるが、座標変換と計量自体はこの条件を越えて定義できる。 座標の順序を $(t,r,\phi,z)$ と固定する。成分行列と逆行列は次の組になる。

$$
\begin{aligned}
(g_{ab})&=\left[\begin{matrix}\Omega^{2} r^{2} - 1 & 0 & \Omega r^{2} & 0\\0 & 1 & 0 & 0\\\Omega r^{2} & 0 & r^{2} & 0\\0 & 0 & 0 & 1\end{matrix}\right]\\
(g^{ab})&=\left[\begin{matrix}-1 & 0 & \Omega & 0\\0 & 1 & 0 & 0\\\Omega & 0 & - \frac{\left(\Omega r - 1\right) \left(\Omega r + 1\right)}{r^{2}} & 0\\0 & 0 & 0 & 1\end{matrix}\right]\\
\det g&=- r^{2}
\end{aligned}
$$

逆行列を掛けて単位行列を得ることを確かめてから、:eqref[eq-gr-christoffel] を使う。同じ下添字の順序を二度計算せず、$b\le c$ の組を選ぶ。中央の列には、座標に依存する計量成分の微分だけを残した。最後の列は微分を実行した値であり、表に出ない成分は添字交換で得られるものを除いて零である。

| 成分 | 微分への代入 | 結果 |
|---|---|---|
| $\Gamma^{r}_{t t}$ | $\frac{1}{2}\left(-\partial_{r}\left(\Omega^{2} r^{2} - 1\right)\right)$ | $- \Omega^{2} r$ |
| $\Gamma^{r}_{t \phi}$ | $\frac{1}{2}\left(-\partial_{r}\left(\Omega r^{2}\right)\right)$ | $- \Omega r$ |
| $\Gamma^{r}_{\phi \phi}$ | $\frac{1}{2}\left(-\partial_{r}\left(r^{2}\right)\right)$ | $- r$ |
| $\Gamma^{\phi}_{t r}$ | $\frac{\Omega}{2}\left(\partial_{r}\left(\Omega^{2} r^{2} - 1\right)\right)- \frac{\left(\Omega r - 1\right) \left(\Omega r + 1\right)}{2 r^{2}}\left(\partial_{r}\left(\Omega r^{2}\right)\right)$ | $\frac{\Omega}{r}$ |
| $\Gamma^{\phi}_{r \phi}$ | $\frac{\Omega}{2}\left(\partial_{r}\left(\Omega r^{2}\right)\right)- \frac{\left(\Omega r - 1\right) \left(\Omega r + 1\right)}{2 r^{2}}\left(\partial_{r}\left(r^{2}\right)\right)$ | $\frac{1}{r}$ |

アフィン径数の微分を点で表し、$w_a$ を座標 $a$ の速度成分とする。ここでは添字を下げる操作を行っていない。測地線の右辺 $-\Gamma^a{}_{bc}w_bw_c$ では $b\ne c$ の項を二度数える。

$$
\begin{aligned}
\ddot{t}&=0\\
\ddot{r}&=r \left(\Omega w_{t} + w_{\phi}\right)^{2}\\
\ddot{\phi}&=- \frac{2 w_{r} \left(\Omega w_{t} + w_{\phi}\right)}{r}\\
\ddot{z}&=0
\end{aligned}
$$

また計量ノルムは

$$
g(\dot x,\dot x)=w_{\phi} \left(\Omega r^{2} w_{t} + r^{2} w_{\phi}\right) + w_{r}^{2} + w_{t} \left(\Omega r^{2} w_{\phi} + w_{t} \left(\Omega^{2} r^{2} - 1\right)\right) + w_{z}^{2}
$$
である。これを径数で微分して上の加速度を代入すると、係数の微分と速度の微分の和が零になる。従って初期のノルムが保存される。

角方向の方程式を $\phi_{\rm in}=\phi+\Omega t$ で書けば $(r^2\dot\phi_{\rm in})^{\cdot}=0$ である。また動径の式は $\ddot r=r\dot\phi_{\rm in}^2$ となる。従って二つの方程式を合わせると、慣性円柱座標での自由な直線運動の式に戻る。たとえば $r=R$、$\dot r=0$ で瞬間的に回転座標に静止して放した粒子の動径の座標加速度は $R\Omega^2\dot t^2$ である。座標に固定されて回る観測者を保つには、この自由軌道からずらす加速度を加える必要がある。

$\Gamma^r{}_{tt}=-\Omega^2r$ は円運動を保つための向心成分を与える。自由粒子の動径式ではこれが右辺に移り、$t$ を径数にすると $\ddot r=r(\dot\phi+\Omega)^2$ となる。$\ddot\phi=-2\dot r(\dot\phi+\Omega)/r$ も同時に使うため、角速度を固定した運動と自由な運動を区別する。$\Omega$ を零にすれば回転を含まない極座標の式となる。
:::

### 円運動の体に固定した枠と Thomas 歳差

一定の半径 $R$ と慣性系の角速度 $\Omega>0$ の円運動を取る。$v=R\Omega<1$、$\gamma=(1-v^2)^{-1/2}$ である。慣性枠の $e_t$ と時刻とともに回る空間ベクトル $e_r,e_\phi$ を使い、運動者の単位時間ベクトルを $u=\gamma(e_t+v e_\phi)$ とする。静止枠の三軸は $E_1=e_r,E_2=\gamma(ve_t+e_\phi),E_3=e_z$ とできる。$D_\tau e_r=\gamma\Omega e_\phi$、$D_\tau e_\phi=-\gamma\Omega e_r$ であり、$e_\phi=\gamma(E_2-vu)$ を代入すると
$$
D_\tau E_1=\gamma^2\Omega E_2-\gamma^2v\Omega u,\qquad
D_\tau E_2=-\gamma^2\Omega E_1.
$$
固有加速度は $a=-\gamma^2v\Omega E_1$ である。Fermi–Walker 輸送する軸 $V=V^1E_1+V^2E_2$ では $D_\tau V=(a\cdot V)u$ だから、空間成分の式は $(V^1)'=\gamma^2\Omega V^2$、$(V^2)'=-\gamma^2\Omega V^1$ となる。初期に $V=E_1$ なら成分の角は $-\gamma^2\Omega\tau$ である。一周の固有時は $2\pi/(\gamma\Omega)$ なので、出発の枠へ戻ったとき残る角は $-2\pi\gamma$、$2\pi$ を法として $-2\pi(\gamma-1)$ である。$v\ll1$ では $-\pi v^2$ が先頭となる。座標が回るだけの場合の $-2\pi$ を除くと、相対論的な歳差が残る。

:::example[少し速い円運動での典型値]
$v=3/5$ なら $\gamma=5/4$ であり、一周の後に輸送軸は体の枠に対して $-\pi/2$ の角をなす。局所慣性軸の変化は世界線の全体に沿って輸送した結果であり、その瞬間の加速度の大きさだけからは決まらない。例えば同じ加速度の大きさを持つ直線的な加速では静止空間内の二軸の回転はこの円運動とは異なる。
:::

### 加速観測者の枠が回転するとき

Minkowski 時空の中心観測者に $u$ と正規空間三脚 $e_i$ を沿わせる。$a_i=g(e_i,a)$ とし、$e_i'=a_i u+W^j{}_i e_j$ と書く。直交性の微分により $W_{ji}=-W_{ij}$ なので、三次元では $W^j{}_i X^i$ を $\boldsymbol\omega\times\boldsymbol X$ と表せる。この回転する三脚を使い、$x=z(\tau)+X^ie_i(\tau)$ と置く。時間方向の微分は $(1+a_iX^i)u+(\boldsymbol\omega\times\boldsymbol X)^j e_j$、空間方向の微分は $e_i$ であるから
$$
\dd s^2=-(1+a_iX^i)^2\dd\tau^2
+\sum_i\big[\dd X^i+(\boldsymbol\omega\times\boldsymbol X)^i\dd\tau\big]^2.
$$
三脚を Fermi–Walker 輸送すると $\boldsymbol\omega=0$ になり、先の加速座標を得る。中心が慣性運動でも三脚だけを回転させれば混合成分が出る。従って回転による時計の効果は、中心の固有加速度の有無と別に計算する必要がある。

:::example[慣性中心のまわりで回る時計]
$a_i=0$、$\boldsymbol\omega=(0,0,\Omega)$ とする。$X,Y,Z$ を固定した時計の率は $\sqrt{1-\Omega^2(X^2+Y^2)}$ である。例えば $\Omega R=3/5$ の円周上の時計は中心の時間 $10$ に対して $8$ の固有時を測る。一方、同じ円周を回転枠に対して運動する光の接線では、$r=R$ を固定した案内路を仮定すると $0=-1+R^2(\dd\phi/\dd\tau+\Omega)^2$ となる。一周に要する中心時刻は二方向で $2\pi R/(1+\Omega R)$ と $2\pi R/(1-\Omega R)$ になる。これは空間内を自由に飛ぶ光の軌道とは異なり、鏡や導波路が光を円周に沿わせる場合の計算である。
:::

### 加速度計と枠を別々に指定する

光や粒子を測るとき、世界線が与える $u,a$ と、空間方向をどのように追跡するかは別のデータである。例えば中心を慣性運動させたまま三脚だけを回すと加速度計は零を読むが、外部の物標の座標は時間とともに変わる。逆に直線上の一様固有加速では加速度計は非零だが、Fermi–Walker 三脚を選べば静止空間内の回転は零である。上の計量の $g_{\tau i}$ は三脚の回転に、中心線上の $\partial_i g_{\tau\tau}=-2a_i$ は固有加速度に対応する。曲がった時空でも中心線上のこの関係は、$X$ 方向の放射測地線で座標を作って $X=0$ を評価すれば得られる。そこでは $\nabla_{\partial_i}\partial_j=0$、$\nabla_{\partial_i}\partial_\tau=\nabla_u e_i$ が成り立ち、後者はねじれ零と座標基底の可換性を使った結果である。

### 加速座標での同時性を慣性座標に戻す

一様固有加速度 $a>0$ の観測者に沿う枠で $x^a=z^a(\tau)+X e_1{}^a(\tau)$ と置くと、$T=(a^{-1}+X)\sinh(a\tau)$、$Z=(a^{-1}+X)\cosh(a\tau)$ である。従って $(T,Z)$ で被覆する範囲は $Z>|T|$ で、逆は $X=\sqrt{Z^2-T^2}-1/a$、$\tau=a^{-1}\operatorname{artanh}(T/Z)$ である。一定 $\tau$ は原点を通る空間的な直線の一部であり、中心観測者の速度に直交する。例えば $X=1/a$ の観測者の固有加速度は $a/2$ で、慣性座標では $Z^2-T^2=4/a^2$ をたどる。

一方、慣性座標での動径光 $Z=Z_0-T$、$Z_0>0$ は $(T,Z)=(Z_0/2,Z_0/2)$ でこの領域を出る。加速座標では $1+aX=a\sqrt{Z_0^2-2Z_0T}$ が零へ近づき、$\tau$ は無限大に発散する。光の慣性座標での接線は一定であり、アフィン径数が有限のところで座標の境界を通る。加速枠の時間の範囲と、その領域に突き当たる光の径数の有限性は矛盾しない。測地線の延長を調べるには、別の座標で曲線と計量の正則性を確認する。

:::example[同じ観測者に沿う二つの空間枠]
初期の速度 $u$ と同じ時間軸を持つ二つの Fermi–Walker 枠は、初期の空間回転を全区間で保つ。実際、一方の三脚 $E_i$ と他方の $F_j$ の内積を微分すると、それぞれの共変微分は $u$ に平行であるため零になる。初期に $F_1=(E_1+E_2)/\sqrt2$ なら、後の各時刻にも同じ関係が成り立つ。加速を変えたときに世界線は変わり得るが、同じ世界線の上で空間軸の初期の取り方を変えただけなら、この一定の回転で二つの記録を照合できる。
:::

## 演習問題

::::exercise{level=basic}
$\dd s^2=-\dd t^2+\dd r^2+r^2\dd\phi^2$ で $\Gamma^\phi{}_{r\phi}$ を $r=2$ で求めよ。

:::solution
$g^{\phi\phi}=r^{-2}$ と $\partial_r g_{\phi\phi}=2r$ を :eqref[eq-gr-christoffel] に代入すると $1/r$、したがって $1/2$ である。
:::
::::

::::exercise{level=basic}
$A(z)=1+2z$ の静止観測者の $z=1$ における固有加速度を求めよ。

:::solution
$a^z=A'/A=2/3$ であり、空間係数が $1$ なのでその大きさも $2/3$ である。
:::
::::

::::exercise{level=basic}
$\lambda=s^2$、$s>0$ で測地線を再表示したときの非アフィン項を求めよ。

:::solution
$\lambda''/\lambda'=2/(2s)=1/s$ である。したがって $D_s x'=x'/s$ であり、元の曲線の像は変わらない。
:::
::::

::::exercise{level=basic}
球面の赤道を一周する平行移動を求めよ。

:::solution
$\theta=\pi/2$ で $\cos\theta=0$ なので、規格化した両成分の微分方程式の右辺は零である。局所の枠は一周で元に戻り、移動は恒等写像である。
:::
::::

::::exercise{level=advanced}
アフィン径数で表した光的測地線の接線と、測地線に沿う任意の平行移動ベクトルとの内積が一定であることを証明せよ。

:::solution
光的な測地線にはアフィン径数を取り、その接線を $k$ とする。この $k$ は $\nabla_k k=0$、平行移動された $V$ は $\nabla_k V=0$ である。計量整合性で $k(g(k,V))=g(\nabla_k k,V)+g(k,\nabla_kV)=0$ となる。$k$ のノルムの符号はこの議論に影響しない。
:::
::::

::::exercise{level=advanced}
同じ座標に関する $g$ と $\widetilde g=e^{2f}g$ の接続係数の差を求めよ。

:::solution
$\widetilde g^{ab}=e^{-2f}g^{ab}$、$\partial_c\widetilde g_{ab}=e^{2f}(\partial_cg_{ab}+2f_cg_{ab})$ である。Christoffel の式の計量微分項は元と一致し、残りは $\delta^a_b f_c+\delta^a_c f_b-g_{bc}g^{ad}f_d$ である。
:::
::::

::::exercise{level=advanced}
前問の変換に対して零測地線の軌道が変わらないことを示せ。

:::solution
接続の差を光的接線 $k^b k^c$ で縮約すると $2k(f)k^a$ となる。したがって変換後の加速度は接線に平行である。新接線を $\widetilde k=e^{-2f}k$ とすれば、積の微分則で追加項が打ち消され $\widetilde\nabla_{\widetilde k}\widetilde k=0$ となる。
:::
::::

::::exercise{level=advanced}
Minkowski 時空で $a_0\tau=\log 3$ の一様加速枠の $u$ と $e$ を求め、直交性を確かめよ。

:::solution
$\cosh(\log3)=5/3$、$\sinh(\log3)=4/3$ なので $u=(5/3,4/3)$、$e=(4/3,5/3)$ である。内積は $-20/9+20/9=0$ である。
:::
::::

::::exercise{level=advanced}
面の計量を $\dd x^2+F(x)^2\dd y^2$、$F>0$ とする。$x=x_0$ を保ちながら $y$ が増える曲線に沿って、開始時に $\partial_x$ を向く単位ベクトルの平行移動を求めよ。またこの曲線が測地線になる条件を求めよ。

:::solution
$\Gamma^x{}_{yy}=-FF'$、$\Gamma^y{}_{xy}=F'/F$ なので、$e_x=\partial_x,e_y=F^{-1}\partial_y$ は $\nabla_{\partial_y}e_x=F'e_y$、$\nabla_{\partial_y}e_y=-F'e_x$ を満たす。$x$ を固定すると $F'(x_0)$ は定数である。従って移動したベクトルは
$$
V=\cos[F'(x_0)(y-y_0)]e_x-\sin[F'(x_0)(y-y_0)]e_y
$$
であり、微分して接続項を加えると零になる。$y$ を一定速度で進むと加速度は $-FF'\dot y^2\partial_x$ であり、接線と直交する。従って再パラメータ化で消せる接線方向の項ではなく、$\dot y\ne0$ で測地線になるのは $F'(x_0)=0$ の場合である。このとき枠の平行移動も恒等となる。
:::
::::

## 確認テスト

::::quiz
接続係数をテンソルの成分として変換すると失われる項を選べ。

- [ ] 逆計量そのもの
- [x] 座標変換の二階微分
- [ ] 縮約の恒等写像

:::explanation
:eqref[eq-gr-connection-change] には $\partial'_\alpha K^\rho{}_\beta=\partial^2x^\rho/(\partial y^\alpha\partial y^\beta)$ がある。これは新しい座標基底を微分する際に生じる。
:::
::::

::::quiz{answer="-3" tol="1e-9"}
平面極座標の $r=3$ における $\Gamma^r{}_{\theta\theta}$ を求めよ。

:::explanation
$\Gamma^r{}_{\theta\theta}=-\tfrac12\partial_r(r^2)=-r$ なので $-3$ である。
:::
::::

::::quiz
$\nabla g=0$ の接続で同じ曲線に沿う二本の平行移動ベクトルについて正しい命題を選べ。

- [x] 内積が一定である。
- [ ] すべての座標成分が一定である。
- [ ] 接線と必ず一致する。

:::explanation
内積を微分すると二本の共変微分の項に分かれ、両方が零である。座標に依存する成分は変わり得る。初期値は任意なので接線との一致も一般には要求されない。
:::
::::

::::quiz{answer="0" tol="1e-9"}
正規座標の原点における $\Gamma^0{}_{11}$ を求めよ。

:::explanation
原点を通る放射測地線は座標で直線である。初速度を第一空間軸に取ると $\Gamma^0{}_{11}(0)=0$ が測地線方程式から従う。
:::
::::

::::quiz{answer="1" tol="1e-7"}
緯線に沿う一周の平行移動で $\cos\theta_0=1/2$ とする。回転角 $2\pi(1-\cos\theta_0)$ を $\pi$ 単位で答えよ。

:::explanation
緯線に沿う枠成分の方程式は $w'=-i\cos\theta_0w$ である。移動後の角は $-2\pi\cos\theta_0=-\pi$ であり、$2\pi$ を加えた $\pi$ に等しい。
:::
::::

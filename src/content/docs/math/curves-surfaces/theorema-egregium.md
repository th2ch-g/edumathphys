---
title: Gauss の驚異の定理と曲面論の基本定理
description: 計量の微分から Gauss 曲率を取り出し、二つの基本形式が曲面を定めるための条件を証明する。
notation:
  - symbol: 'g_{ij},h_{ij}'
    meaning: 第一および第二基本形式の成分
  - symbol: '\Gamma^k_{ij}'
    meaning: Christoffel 記号
  - symbol: 'R^l_{kij}'
    meaning: 曲率作用素の座標成分
  - symbol: '\omega'
    meaning: 正規直交標構の接続形式
---

## この章の目標

- 等長写像と共形写像を、計量の変換により判定する。
- Christoffel 記号を計量から計算し、Gauss 方程式と Codazzi 方程式を導く。
- 曲面の Gauss 曲率が計量だけで決まる理由を示す。
- 二つの基本形式から曲面を再構成し、等温座標で曲率を求める。

## 等長写像と共形写像

円柱の一部分は引き伸ばさずに平面へ開ける。球面を立体射影で平面へ移すと長さは変わるが、曲線が交わる角は保たれる。何が保存されるかを計量の等式で表す。

:::definition[局所等長写像と共形写像]{#def-cs-isometric-conformal}
滑らかな局所微分同相 $f:S\to S'$ が $g'_{f(p)}(df_p w,df_p z)=g_p(w,z)$ を満たすとき、$f$ を :term[局所等長写像]{en="local isometry" yomi="きょくしょとうちょうしゃぞう"} という。右辺に滑らかな正関数 $\lambda(p)^2$ を掛けたものに等しいとき、$f$ を :term[共形写像]{en="conformal map" yomi="きょうけいしゃぞう"} という。局所等長写像が全域で微分同相なら等長写像である。
:::

:::proposition[長さと角の保存]{#prop-cs-isometry-measure}
局所等長写像は曲線の長さを保ち、共形写像は接ベクトル間の角を保つ。逆に曲線の長さをすべて保つ滑らかな局所微分同相は局所等長写像であり、角をすべて保つものは共形写像である。
:::

:::proof
前向きは曲線の速さと角の余弦の式に代入する。逆に長さを保つとき、任意の点を出る任意の接速度で短い曲線を取り、両辺の長さをパラメータ区間の長さで割って区間を縮める。連続性から出発点での速さが一致する。二次形式の偏極恒等式で内積も一致する。角を保つ場合は正規直交基底 $e_1,e_2$ の像が直交する。$e_1+e_2$ と $e_1-e_2$ も直交するため、像の内積から $\norm{df(e_1)}^2-\norm{df(e_2)}^2=0$ を得る。二方向の倍率は等しく、局所基底で計量の跡を取れば倍率の二乗が滑らかであることも分かる。
:::

:::example[円柱の局所展開]
半径 $2$ の円柱で $f(X(u,v))=(2u,v)$ とすれば計量はいずれも $4\,\dd u^2+\dd v^2$ である。$u$ は局所角度なのでこの写像は局所的である。角座標を一周して同じ点に戻ると、平面側の第一座標は $4\pi$ 変わる。
:::

:::example[立体射影の倍率]
[正則曲面と第一基本形式](/math/curves-surfaces/surfaces/) で求めた単位球面の立体射影の逆写像 $X=(2u,2v,u^2+v^2-1)/(1+u^2+v^2)$ に対し、計量は

$$
\mathrm I=\frac{4(\dd u^2+\dd v^2)}{(1+u^2+v^2)^2}.
$$

実際 $D=1+u^2+v^2$ とし、分子を $Y=(2u,2v,u^2+v^2-1)$ とすると $\inner{Y}{Y_i}=DD_i$ である。$\inner{X_i}{X_j}=(\inner{Y_i}{Y_j}-D_iD_j)/D^2$ を得て、分子は $4\delta_{ij}$ となる。原点では長さが $2$ 倍になるが、$u=1,v=0$ では倍率は $1$ である。
:::

### 等温座標の間の写像を複素微分で調べる

:ref[def-cs-isothermal] の等温座標では二本の座標軸の長さが等しく、かつ直交している。この場合、共形写像の微分が二本の座標速度にどう作用するかは、Cauchy–Riemann の方程式で読み取れる。

:::proposition[等温座標による共形性の判定]
定義域と値域に等温座標を取り、前者の計量を $\lambda^2(\dd u^2+\dd v^2)$、後者の計量を $\mu^2(\dd x^2+\dd y^2)$ とする。向きを保つ局所微分同相 $(u,v)\mapsto(x,y)$ は、共形であるときかつそのときに限り、$x_u=y_v$、$x_v=-y_u$ を満たす。その際の長さの倍率は $(\mu\circ f)\sqrt{x_u^2+y_u^2}/\lambda$ である。
:::

:::proof
二つの計量はいずれも正のスカラー倍の Euclid 計量なので、微分の列を座標上で比べればよい。共形性により列は直交し長さが等しい。行列式が正なので第二列は第一列を正の直角に回したもの、すなわち $(x_v,y_v)=(-y_u,x_u)$ となる。逆にこの成分の関係があれば二列は直交し長さも等しい。第一列の Euclid 長を値域の長さへ変換するには $\mu\circ f$ を掛け、元の座標速度の長さ $\lambda$ で割ると倍率になる。
:::

:::example[複素写像を球面計量で測る]
立体射影で球を表示した座標 $(u,v)$ に、平面の複素座標を $w=x+iy$ として $u+iv=w^2$ を代入する。$u=x^2-y^2$、$v=2xy$ なので Cauchy–Riemann の式を満たす。Jacobian は $4(x^2+y^2)$ であり、原点以外で局所微分同相になる。単位球側の等温倍率は $2/(1+|w|^4)$、微分の Euclid 長さの倍率は $2|w|$ だから、平面の Euclid 計量から球への長さの倍率は
$$
\frac{4|w|}{1+|w|^4}
$$
である。例えば $|w|=1$ では倍率が $2$ である。原点では写像は滑らかだが微分が零になり、局所微分同相としての共形写像の条件から外れる。共形性を調べるときは角の保存だけでなく、微分の階数も確認する。
:::

:::example[円柱の展開に高さのせん断を加える]
半径 $3$ の円柱の座標は $(u,v)=(\text{角度},\text{高さ})$ である。新しい平面の座標で曲面へ写す。

点の写像をまず座標で書き、入力の二つの方向を個別に微分する。微分の列はそれぞれ $\partial_x,\partial_y$ の像である。

$$
(u,v)=\left(\frac{x}{3},x + y\right),\qquad D=\begin{pmatrix}\frac{1}{3} & 0\\1 & 1\end{pmatrix},\qquad \det D=\frac{1}{3}
$$

入力と出力の座標値が違うので、出力側の計量は像の点で評価する。二本の像ベクトルの内積を取って入力側へ戻す。

$$
\begin{aligned}
&g_{\mathrm{target}}\circ F=\begin{pmatrix}9 & 0\\0 & 1\end{pmatrix},\\
&F^*g=D^\top(g_{\mathrm{target}}\circ F)D=\begin{pmatrix}2 & 1\\1 & 1\end{pmatrix}
\end{aligned}
$$

二つの座標速度 $(1,0)$ と $(1,1)$ を例に取ると、引き戻した行列からそのなす角を計算できる。微分の像の長さと内積を使う計算と一致する。

$$
\cos\angle((1,0),(1,1))=\frac{3 \sqrt{10}}{10},\qquad \det(F^*g)=1
$$

入力側に標準 Euclid 計量を入れると、この行列は恒等でもそのスカラー倍でもない。従ってこの写像には長さを変える方向と角を変える組がある。入力側の計量を引き戻しに選べば、同じ写像が局所等長になる。
:::

:::example[複素逆数が作る共形倍率]
原点を除いた平面で、複素数 $x+iy$ をその逆数へ送る写像を実変数で調べる。

点の写像をまず座標で書き、入力の二つの方向を個別に微分する。微分の列はそれぞれ $\partial_x,\partial_y$ の像である。

$$
\begin{aligned}
&(u,v)=\left(\frac{x}{x^{2} + y^{2}},- \frac{y}{x^{2} + y^{2}}\right),\\
&D=\begin{pmatrix}- \frac{\left(x - y\right) \left(x + y\right)}{\left(x^{2} + y^{2}\right)^{2}} & - \frac{2 x y}{\left(x^{2} + y^{2}\right)^{2}}\\\frac{2 x y}{\left(x^{2} + y^{2}\right)^{2}} & - \frac{\left(x - y\right) \left(x + y\right)}{\left(x^{2} + y^{2}\right)^{2}}\end{pmatrix},\\
&\det D=\frac{1}{\left(x^{2} + y^{2}\right)^{2}}
\end{aligned}
$$

入力と出力の座標値が違うので、出力側の計量は像の点で評価する。二本の像ベクトルの内積を取って入力側へ戻す。

$$
\begin{aligned}
&g_{\mathrm{target}}\circ F=\begin{pmatrix}1 & 0\\0 & 1\end{pmatrix},\\
&F^*g=D^\top(g_{\mathrm{target}}\circ F)D=\begin{pmatrix}\frac{1}{\left(x^{2} + y^{2}\right)^{2}} & 0\\0 & \frac{1}{\left(x^{2} + y^{2}\right)^{2}}\end{pmatrix}
\end{aligned}
$$

二つの座標速度 $(1,0)$ と $(1,1)$ を例に取ると、引き戻した行列からそのなす角を計算できる。微分の像の長さと内積を使う計算と一致する。

$$
\cos\angle((1,0),(1,1))=\frac{\sqrt{2}}{2},\qquad \det(F^*g)=\frac{1}{\left(x^{2} + y^{2}\right)^{4}}
$$

引き戻しは点ごとに同じ正の係数を二方向へ掛ける。従って角は保つが、原点からの距離が $1$ の点以外では接ベクトルの長さが変わる。
:::

## Christoffel 記号

曲面の二階微分を接平面と法線へ分ける。法線成分は前章で求めたが、接成分にも座標基底の変化という情報がある。多くの項が現れるものの、計量の一階微分からすべて求められる。

:::definition[Christoffel 記号]{#def-cs-christoffel}
$x^1=u,x^2=v$、$X_i=\partial_iX$ と書く。$X_{ij}=\Gamma^k_{ij}X_k+h_{ij}N$ と展開したときの接成分を :term[Christoffel 記号]{en="Christoffel symbols" yomi="くりすとっふぇるきごう"} という。この章では同じ項に上下で二回現れる添字を二方向にわたって足す。
:::

:::proposition[計量から記号を求める公式]{#prop-cs-christoffel-formula}
逆計量を $(g^{ij})$ とすると

$$
\Gamma^k_{ij}=\tfrac12 g^{kl}(\partial_i g_{jl}+\partial_jg_{il}-\partial_l g_{ij}).
\label{eq-cs-christoffel}
$$

:::

:::proof
$\partial_i g_{jl}=\inner{X_{ij}}{X_l}+\inner{X_j}{X_{il}}$ を使う。$\partial_i g_{jl}+\partial_jg_{il}-\partial_lg_{ij}$ では、偏微分の可換性により後二種の内積が相殺され、$2\inner{X_{ij}}{X_l}=2\Gamma^k_{ij}g_{kl}$ が残る。左から逆計量を掛けると式が出る。この導出から $\Gamma^k_{ij}=\Gamma^k_{ji}$ と $\partial_i g_{jl}=\Gamma^k_{ij}g_{kl}+\Gamma^k_{il}g_{jk}$ も従う。
:::

座標を変えると Christoffel 記号には二階微分から来る項が加わる。具体的に旧座標 $x^i=x^i(y)$ に対して

$$
\widetilde \Gamma^c_{ab}
=\frac{\partial y^c}{\partial x^k}
 \left(\frac{\partial^2x^k}{\partial y^a\partial y^b}
 +\Gamma^k_{ij}\frac{\partial x^i}{\partial y^a}\frac{\partial x^j}{\partial y^b}\right).
$$

これは $\partial_a\partial_b(X\circ x)$ を連鎖律で展開して接成分を新基底に直した式である。この二階の項があるため、$\Gamma$ の値そのものでは内在的な曲がりを測れない。

:::example[平面の極座標]
$g=\dd r^2+r^2\dd\theta^2$ では $\Gamma^r_{\theta\theta}=-r$、$\Gamma^\theta_{r\theta}=\Gamma^\theta_{\theta r}=1/r$ が非零である。例えば $\Gamma^r_{\theta\theta}=-\partial_r(r^2)/2=-r$ である。直交直線座標では全成分が零になる。
:::

:::example[回転面]
$g=\dd s^2+r(s)^2\dd\theta^2$ では同じ計算で $\Gamma^s_{\theta\theta}=-rr'$、$\Gamma^\theta_{s\theta}=r'/r$ を得る。単位球面の余緯度では $r=\sin s$ であり、$s=\pi/2$ では三つの非零候補がすべて零になる。
:::

### 一点で接続係数を消す座標を選ぶ

Christoffel 記号に座標の二階微分が入ることを逆に利用すると、任意の一点で全成分が零になる座標を作れる。これにより、その点で曲率を計算する式から接続係数の積が消える。ただしその近くでの係数の微分まで消えるとは限らない。

:::proposition[計量の一次微分が消える座標]
任意の点 $p$ において、$g_{ij}(p)=\delta_{ij}$ かつ $\partial_k g_{ij}(p)=0$ を満たす局所座標が存在する。
:::

:::proof
最初に旧座標 $x$ で $x(p)=0$、$g(p)=I$ になるよう線形変換する。次に新座標 $y$ を
$$
x^k=y^k-\frac12\sum_{a,b}\Gamma^k_{ab}(p)y^a y^b
$$
で指定する。一次微分は原点で恒等なので逆関数定理で $y$ は座標になる。二階微分は $\Gamma^k_{ab}=\Gamma^k_{ba}$ により $-\Gamma^k_{ab}(p)$ となる。座標変換の式に代入すると元の Christoffel 記号が相殺され、新しい記号は全て零である。計量の一次微分の式 $\partial_i g_{jl}=\Gamma^k_{ij}g_{kl}+\Gamma^k_{il}g_{jk}$ により新しい計量の一次微分は零である。
:::

:::example[平面の極座標で原点を取り直す]
極半径 $r=2$、角 $\theta=0$ の点で $x^1=r-2,\ x^2=2\theta$ とする。計量は $(\dd x^1)^2+(1+x^1/2)^2(\dd x^2)^2$ なので、原点で $\Gamma^1_{22}=-1/2,\ \Gamma^2_{12}=1/2$。上の構成で $x^1=y^1+(y^2)^2/4$、$x^2=y^2-y^1y^2/2$ と取る。例えば新しい第二方向は第一成分に $y^2/2$、第二成分に $1-y^1/2$ を持つ。後者と元の倍率 $1+x^1/2$ を掛けると一次の $y^1$ 項が相殺される。残る交差項も、第一方向側の第二成分 $-y^2/2$ が第二方向側の第一成分を打ち消すため、計量の一次微分は零になる。
:::

## Gauss 方程式と Codazzi–Mainardi 方程式

座標方向の順序を変えて三階微分しても、空間の位置ベクトルの結果は同じである。接成分と法線成分をそれぞれ比較すると、二つの基本形式の間に条件が生じる。

:::definition[曲率成分と第二基本形式の共変微分]{#def-cs-coordinate-curvature}
座標から作る量を

$$
\begin{aligned}
R^l_{kij}&=\partial_i\Gamma^l_{jk}-\partial_j\Gamma^l_{ik}
+\Gamma^m_{jk}\Gamma^l_{im}-\Gamma^m_{ik}\Gamma^l_{jm},\\
(\nabla_i h)_{jk}&=\partial_i h_{jk}-\Gamma^m_{ij}h_{mk}-\Gamma^m_{ik}h_{jm}
\end{aligned}
$$

と書く。この $R$ は、順序 $i,j$ の接方向微分の交換子に対応する。
:::

:::theorem[Gauss 方程式と Codazzi–Mainardi 方程式]{#thm-cs-gauss-codazzi}
正則曲面の二つの基本形式について

$$
\begin{aligned}
g_{lm}R^m_{kij}&=h_{jk}h_{il}-h_{ik}h_{jl},\\
(\nabla_i h)_{jk}&=(\nabla_j h)_{ik}
\end{aligned}
\label{eq-cs-gauss-codazzi}
$$

が成り立つ。
:::

:::proof
Weingarten の式は $N_i=-A_i^l X_l$、$A_i^l=g^{lm}h_{im}$ である。$\partial_i X_{jk}$ を $\partial_i(\Gamma^m_{jk}X_m+h_{jk}N)$ として展開すると、接成分の係数は $\partial_i\Gamma^l_{jk}+\Gamma^m_{jk}\Gamma^l_{im}-h_{jk}A_i^l$、法線成分は $\partial_i h_{jk}+\Gamma^m_{jk}h_{im}$ である。$i,j$ を反転した式との差は零。接成分を整理し計量を掛ければ第一式となる。法線成分を移項し、$\Gamma^m_{ij}=\Gamma^m_{ji}$ を用いれば第二式になる。
:::

:::example[平面計量と不可能な第二基本形式]
$g_{ij}=\delta_{ij}$ では $\Gamma=0$、$R=0$ である。$h=\dd u^2+\dd v^2$ なら $h_{22}h_{11}-h_{12}^2=1$ で Gauss 方程式に反する。従ってこの組を同時に満たす曲面はない。
:::

:::example[Codazzi の独立な制約]
同じ計量で $h=v\,\dd u^2$ なら $\det h=0$ なので Gauss 方程式には適合する。しかし $\partial_v h_{11}=1$ に対し $\partial_u h_{21}=0$ なので Codazzi 方程式に反する。
:::

### 曲率テンソルに残る一つの成分

$\nabla_XY$ は接方向に空間微分を射影したものであり、$[X,Y]$ は関数に対する方向微分の交換子である。これらを用いて $R(X,Y)Z=\nabla_X\nabla_YZ-\nabla_Y\nabla_XZ-\nabla_{[X,Y]}Z$ と置く。最後の項には二方向の微分を作用素として交換した際に生じる一次の項を除く役目がある。座標基底では $[\partial_u,\partial_v]=0$ なので、定義は :ref[def-cs-coordinate-curvature] に一致する。

:::proposition[曲面上の曲率作用素]{#prop-cs-curvature-operator}
曲面の正規直交基底 $(e_1,e_2)$ に対して $R(e_1,e_2)e_2=Ke_1$、$R(e_1,e_2)e_1=-Ke_2$ である。従って
$$
R(X,Y)Z=K\bigl(g(Y,Z)X-g(X,Z)Y\bigr)
$$
が成り立つ。
:::

:::proof
積の法則を $R(fX,Y)Z$ で展開すると、$\nabla_Y(f\nabla_X Z)$ から出る $Y(f)\nabla_XZ$ と $[fX,Y]=f[X,Y]-Y(f)X$ から出る項が相殺される。$Y$ については反対称性で同じである。$Z$ については $R(X,Y)(fZ)$ の二階微分が $X(Yf)-Y(Xf)-[X,Y]f=0$、一階微分に掛かる接続の項も逆符号で相殺される。従って点ごとの三重線形性がある。

計量の積の法則を二回使って差を取ると $g(R(X,Y)Z,W)+g(Z,R(X,Y)W)=0$。実際、中間に出る $g(\nabla_Y Z,\nabla_XW)$ などの交差項は順序を反転した式にも一度ずつ現れ、残るスカラーの交換子は $[X,Y]g(Z,W)$ である。よって $R(e_1,e_2)$ は反対称作用素であり、二次元では $g(R(e_1,e_2)e_2,e_1)$ の一成分だけで決まる。この成分は Gauss 方程式で $K$ になる。一般形は基底上で比較して三重線形性を使えばよい。
:::

:::example[球面の基本形式の適合条件]
$0<u<\pi$ で半径 $2$ の球面の座標を取る。外向き法線で $h=-g/2$ を候補とする。

第一基本形式から内積の行列を読み、逆行列を出す。逆行列を使うのは、内積で求めた下付き成分から接ベクトルの座標成分へ戻すためである。

$$
\begin{aligned}
&g=\begin{pmatrix}4 & 0\\0 & 4 \sin^{2}{\left(u \right)}\end{pmatrix},\\
&\det(g)=16 \sin^{2}{\left(u \right)},\\
&g^{-1}=\begin{pmatrix}\frac{1}{4} & 0\\0 & \frac{1}{4 \sin^{2}{\left(u \right)}}\end{pmatrix}
\end{aligned}
$$

逆行列と計量の一階微分を Christoffel の公式に代入する。零でない成分を次に示す。下の二添字を交換した成分は同じである。

$$
\Gamma^{u}_{vv}=- \frac{\sin{\left(2 u \right)}}{2},\quad \Gamma^{v}_{uv}=\frac{\cos u}{\sin u}
$$

さらに二方向の微分の交換子を計算する。$R(\partial_u,\partial_v)\partial_v$ の座標列に左から $g$ の第一行を掛ければ、$\partial_u$ との内積になる。

$$
\begin{aligned}
&[R(\partial_u,\partial_v)\partial_v]=\begin{pmatrix}\sin^2u\\0\end{pmatrix},\\
&K=\frac{g_{1m}R^m_{212}}{\det g}=\frac{1}{4}
\end{aligned}
$$

候補の第二基本形式は、まず行列式を比較して Gauss 方程式を確認する。さらに各 $j$ について右端の微分と接続の項を足し、Codazzi の差が零であることを確かめる。

$$
\begin{aligned}
&h=\begin{pmatrix}-2 & 0\\0 & - 2 \sin^{2}{\left(u \right)}\end{pmatrix},\\
&\det h=4 \sin^{2}{\left(u \right)},\\
&\partial_u h_{2j}-\partial_v h_{1j}+h_{1m}\Gamma^m_{2j}-h_{2m}\Gamma^m_{1j}=0
\end{aligned}
$$

Codazzi が成立するのは $h$ が定数倍の計量であり、その共変微分が消えるからである。座標成分が一定である必要はない。
:::

:::example[回転形計量の負曲率と第二基本形式]
$g=\dd u^2+e^{2u}\dd v^2$ に対して $h=2e^u \dd u\dd v$ を代入し、Gauss 方程式の符号を確認する。

第一基本形式から内積の行列を読み、逆行列を出す。逆行列を使うのは、内積で求めた下付き成分から接ベクトルの座標成分へ戻すためである。

$$
\begin{aligned}
&g=\begin{pmatrix}1 & 0\\0 & e^{2 u}\end{pmatrix},\\
&\det(g)=e^{2 u},\\
&g^{-1}=\begin{pmatrix}1 & 0\\0 & e^{- 2 u}\end{pmatrix}
\end{aligned}
$$

逆行列と計量の一階微分を Christoffel の公式に代入する。零でない成分を次に示す。下の二添字を交換した成分は同じである。

$$
\Gamma^{u}_{vv}=- e^{2 u},\quad \Gamma^{v}_{uv}=1
$$

さらに二方向の微分の交換子を計算する。$R(\partial_u,\partial_v)\partial_v$ の座標列に左から $g$ の第一行を掛ければ、$\partial_u$ との内積になる。

$$
[R(\partial_u,\partial_v)\partial_v]=\begin{pmatrix}- e^{2 u}\\0\end{pmatrix},\qquad K=\frac{g_{1m}R^m_{212}}{\det g}=-1
$$

行列式が合っても、この $h$ は Codazzi に適合しない。実際 $h_{12}=e^u$ の $u$ 微分は $e^u$ であり、交換子の $j=1$ の差は $2e^u$ となる。曲率だけから第二基本形式は決まらない。
:::

:::example[平面の斜め方向だけに曲率を持たせる]
Euclid 計量のまま $h=(1+(u+v)^2)(\dd u+\dd v)^2$ とする。零でない成分が増えても適合条件を確かめられる。

第一基本形式から内積の行列を読み、逆行列を出す。逆行列を使うのは、内積で求めた下付き成分から接ベクトルの座標成分へ戻すためである。

$$
g=\begin{pmatrix}1 & 0\\0 & 1\end{pmatrix},\qquad \det(g)=1,\qquad g^{-1}=\begin{pmatrix}1 & 0\\0 & 1\end{pmatrix}
$$

逆行列と計量の一階微分を Christoffel の公式に代入する。零でない成分を次に示す。下の二添字を交換した成分は同じである。

$$
\Gamma^k_{ij}=0\quad\text{for all }i,j,k
$$

さらに二方向の微分の交換子を計算する。$R(\partial_u,\partial_v)\partial_v$ の座標列に左から $g$ の第一行を掛ければ、$\partial_u$ との内積になる。

$$
[R(\partial_u,\partial_v)\partial_v]=\begin{pmatrix}0\\0\end{pmatrix},\qquad K=\frac{g_{1m}R^m_{212}}{\det g}=0
$$

候補の第二基本形式は、まず行列式を比較して Gauss 方程式を確認する。さらに各 $j$ について右端の微分と接続の項を足し、Codazzi の差が零であることを確かめる。

$$
\begin{aligned}
&h=\begin{pmatrix}u^{2} + 2 u v + v^{2} + 1 & u^{2} + 2 u v + v^{2} + 1\\u^{2} + 2 u v + v^{2} + 1 & u^{2} + 2 u v + v^{2} + 1\end{pmatrix},\\
&\det h=0,\\
&\partial_u h_{2j}-\partial_v h_{1j}+h_{1m}\Gamma^m_{2j}-h_{2m}\Gamma^m_{1j}=0
\end{aligned}
$$

行列は階数一で、全成分の $u$ 微分と $v$ 微分が一致する。従って対角だけを残す必要はない。$(u+v)/\sqrt2$ を横断面の弧長、$(u-v)/\sqrt2$ を母線方向とすれば、横断面だけに曲率を持つ円柱型の実現を作れる。
:::

## Gauss の驚異の定理

Gauss 曲率を主曲率の積として定義したときは、空間の中の法線を使った。しかし二階の差分を取ると、同じ量が曲面の上の測定だけで表せる。

:::theorem[Gauss の驚異の定理]{#thm-cs-egregium}
Gauss 曲率は第一基本形式とその二階までの座標微分から決まる。特に局所等長写像は対応する点で Gauss 曲率を保存する。
:::

:::proof
:ref[thm-cs-gauss-codazzi] に $(l,k,i,j)=(1,2,1,2)$ を代入すると

$$
K=\frac{g_{1m}R^m_{212}}{\det(g_{ij})}.
\label{eq-cs-intrinsic-k}
$$

右辺は :eqref[eq-cs-christoffel] から計量だけで求められる。局所等長写像があるときは片方の座標を写像で運び、同じ座標関数を使えば計量成分が等しい。従って両辺の値が等しい。
:::

座標域に正定値計量だけを与えた場合にも、右辺を曲率の定義にできる。:eqref[eq-cs-christoffel] から作った係数を用いて接ベクトル場の成分を微分すると、座標変換則の二階項が成分を変換してから微分する際の項と相殺する。従って接方向の微分 $\nabla$ は座標に依存しない。:ref[prop-cs-curvature-operator] の証明のうち、$R$ の三重線形性と反対称性には空間内の法線を使わない。よって正規直交基底での $\inner{R(e_1,e_2)e_2}{e_1}$ は基底によらず、:eqref[eq-cs-intrinsic-k] に一致する。この定義は、座標域に与えた計量が空間中の曲面から得られる場合には従来の Gauss 曲率と一致する。

:::example[平面と球面の相違]
半径 $2$ の球面の曲率は $1/4$、平面は $0$ である。従って球面のどれほど小さな開集合も平面へ局所等長には移せない。立体射影は角を保つが局所等長ではない。
:::

:::example[内在的な曲率と平均曲率の差]
円柱と平面は局所等長なのでいずれも $K=0$ となる。一方で半径 $2$ の円柱の平均曲率の絶対値は $1/4$、平面では零である。平均曲率は第一基本形式だけでは決まらない。
:::

:::example[斜交計量の接続項が零になる場合]
平面に線形な座標を取り、計量に交差項があっても微分がどうなるかを見る。

第一基本形式から内積の行列を読み、逆行列を出す。逆行列を使うのは、内積で求めた下付き成分から接ベクトルの座標成分へ戻すためである。

$$
\begin{aligned}
&g=\begin{pmatrix}5 & 1\\1 & 2\end{pmatrix},\\
&\det(g)=9,\\
&g^{-1}=\begin{pmatrix}\frac{2}{9} & - \frac{1}{9}\\- \frac{1}{9} & \frac{5}{9}\end{pmatrix}
\end{aligned}
$$

逆行列と計量の一階微分を Christoffel の公式に代入する。零でない成分を次に示す。下の二添字を交換した成分は同じである。

$$
\Gamma^k_{ij}=0\quad\text{for all }i,j,k
$$

さらに二方向の微分の交換子を計算する。$R(\partial_u,\partial_v)\partial_v$ の座標列に左から $g$ の第一行を掛ければ、$\partial_u$ との内積になる。

$$
[R(\partial_u,\partial_v)\partial_v]=\begin{pmatrix}0\\0\end{pmatrix},\qquad K=\frac{g_{1m}R^m_{212}}{\det g}=0
$$

座標線がなす角は直角ではないが、基底ベクトルは空間で一定である。従ってこの場合に一次微分の接続項は現れない。
:::

:::example[極座標の係数が半径で変わる場合]
$u>0$ で $u$ を半径、$v$ を角度に取る。接続係数の零でない項と曲率を区別する。

第一基本形式から内積の行列を読み、逆行列を出す。逆行列を使うのは、内積で求めた下付き成分から接ベクトルの座標成分へ戻すためである。

$$
\begin{aligned}
&g=\begin{pmatrix}1 & 0\\0 & u^{2}\end{pmatrix},\\
&\det(g)=u^{2},\\
&g^{-1}=\begin{pmatrix}1 & 0\\0 & \frac{1}{u^{2}}\end{pmatrix}
\end{aligned}
$$

逆行列と計量の一階微分を Christoffel の公式に代入する。零でない成分を次に示す。下の二添字を交換した成分は同じである。

$$
\Gamma^{u}_{vv}=- u,\quad \Gamma^{v}_{uv}=\frac{1}{u}
$$

さらに二方向の微分の交換子を計算する。$R(\partial_u,\partial_v)\partial_v$ の座標列に左から $g$ の第一行を掛ければ、$\partial_u$ との内積になる。

$$
[R(\partial_u,\partial_v)\partial_v]=\begin{pmatrix}0\\0\end{pmatrix},\qquad K=\frac{g_{1m}R^m_{212}}{\det g}=0
$$

$\Gamma^u_{vv}$ と $\Gamma^v_{uv}$ はいずれも零でない。それでも交換子では微分と積の寄与が相殺し、平面の曲率は零になる。
:::

:::example[原点で座標が一致しても二階項が残る]
$u,v$ を新しい座標とし、平面の直交座標を $(u,v+u^2)$ とする。接続係数を計量だけから復元する。

第一基本形式から内積の行列を読み、逆行列を出す。逆行列を使うのは、内積で求めた下付き成分から接ベクトルの座標成分へ戻すためである。

$$
\begin{aligned}
&g=\begin{pmatrix}4 u^{2} + 1 & 2 u\\2 u & 1\end{pmatrix},\\
&\det(g)=1,\\
&g^{-1}=\begin{pmatrix}1 & - 2 u\\- 2 u & 4 u^{2} + 1\end{pmatrix}
\end{aligned}
$$

逆行列と計量の一階微分を Christoffel の公式に代入する。零でない成分を次に示す。下の二添字を交換した成分は同じである。

$$
\Gamma^{v}_{uu}=2
$$

さらに二方向の微分の交換子を計算する。$R(\partial_u,\partial_v)\partial_v$ の座標列に左から $g$ の第一行を掛ければ、$\partial_u$ との内積になる。

$$
[R(\partial_u,\partial_v)\partial_v]=\begin{pmatrix}0\\0\end{pmatrix},\qquad K=\frac{g_{1m}R^m_{212}}{\det g}=0
$$

$u=0$ で計量が $I$ に戻っても $\Gamma^v_{uu}=2$ が残る。座標の基底をその点で正規直交にしたことと、近くでも平行にしたことを区別できる。
:::

### 等長性を否定するには曲率の範囲を調べる

二つの曲面のうち一方の曲率が一定なら、他方が候補となるためには対応する開集合で同じ定数を取る必要がある。点での曲率の値が一つ一致しても、近傍の等長性が従うわけではない。従って曲率の座標関数を一度求め、点から離れたときの変化を調べるとよい。

:::example[一点の曲率は同じでも近くで異なる]
$g_1=\dd u^2+\dd v^2$ と $g_2=\dd u^2+(2+u^3)^2\dd v^2$ を $|u|<1$ で比べる。前者は $K_1=0$ である。後者は $f=2+u^3$ として Christoffel の公式を使うと $\Gamma^u_{vv}=-ff'$、$\Gamma^v_{uv}=f'/f$ であり、曲率成分は $-ff''$。よって $K_2=-6u/(2+u^3)$ となる。$u=0$ の線では $K_1=K_2$ だが、その線のどの開近傍にも $K_2\ne0$ の点がある。従ってその近傍は平面のどの開集合とも局所等長にならない。
:::

:::example[同じ計量を持つ懸垂面と螺旋面]
共役な二つの極小面の等温計量 $\cosh^2u(\dd u^2+\dd v^2)$ について曲率を求める。

第一基本形式から内積の行列を読み、逆行列を出す。逆行列を使うのは、内積で求めた下付き成分から接ベクトルの座標成分へ戻すためである。

$$
\begin{aligned}
&g=\begin{pmatrix}\cosh^{2}{\left(u \right)} & 0\\0 & \cosh^{2}{\left(u \right)}\end{pmatrix},\\
&\det(g)=\cosh^{4}{\left(u \right)},\\
&g^{-1}=\begin{pmatrix}\frac{1}{\cosh^{2}{\left(u \right)}} & 0\\0 & \frac{1}{\cosh^{2}{\left(u \right)}}\end{pmatrix}
\end{aligned}
$$

逆行列と計量の一階微分を Christoffel の公式に代入する。零でない成分を次に示す。下の二添字を交換した成分は同じである。

$$
\begin{aligned}
&\Gamma^{u}_{uu}=\tanh{\left(u \right)},\\
&\Gamma^{u}_{vv}=- \tanh{\left(u \right)},\\
&\Gamma^{v}_{uv}=\tanh{\left(u \right)}
\end{aligned}
$$

さらに二方向の微分の交換子を計算する。$R(\partial_u,\partial_v)\partial_v$ の座標列に左から $g$ の第一行を掛ければ、$\partial_u$ との内積になる。

$$
\begin{aligned}
&[R(\partial_u,\partial_v)\partial_v]=\begin{pmatrix}- \frac{1}{\cosh^{2}{\left(u \right)}}\\0\end{pmatrix},\\
&K=\frac{g_{1m}R^m_{212}}{\det g}=- \frac{1}{\cosh^{4}{\left(u \right)}}
\end{aligned}
$$

どちらの面も同じ内在曲率を持つ。先に法線や二つの主方向を計算せずに、局所の等長性が曲率の一致を説明する。
:::

:::example[零曲率の非線形な目盛り]
計量が $e^{4u}(\dd u^2+\dd v^2)$ となる領域を調べる。

第一基本形式から内積の行列を読み、逆行列を出す。逆行列を使うのは、内積で求めた下付き成分から接ベクトルの座標成分へ戻すためである。

$$
\begin{aligned}
&g=\begin{pmatrix}e^{4 u} & 0\\0 & e^{4 u}\end{pmatrix},\\
&\det(g)=e^{8 u},\\
&g^{-1}=\begin{pmatrix}e^{- 4 u} & 0\\0 & e^{- 4 u}\end{pmatrix}
\end{aligned}
$$

逆行列と計量の一階微分を Christoffel の公式に代入する。零でない成分を次に示す。下の二添字を交換した成分は同じである。

$$
\Gamma^{u}_{uu}=2,\quad \Gamma^{u}_{vv}=-2,\quad \Gamma^{v}_{uv}=2
$$

さらに二方向の微分の交換子を計算する。$R(\partial_u,\partial_v)\partial_v$ の座標列に左から $g$ の第一行を掛ければ、$\partial_u$ との内積になる。

$$
[R(\partial_u,\partial_v)\partial_v]=\begin{pmatrix}0\\0\end{pmatrix},\qquad K=\frac{g_{1m}R^m_{212}}{\det g}=0
$$

複素写像 $u+iv\mapsto e^{2(u+iv)}/2$ の微分の絶対値は $e^{2u}$ である。従って十分小さな領域では、平面への実際の局所等長写像も作れる。
:::

:::example[同じ平坦計量の対数展開]
平坦な共形計量 $e^{2x}(\dd x^2+\dd y^2)$ の展開を作り、局所等長性の定義を確認する。

点の写像をまず座標で書き、入力の二つの方向を個別に微分する。微分の列はそれぞれ $\partial_x,\partial_y$ の像である。

$$
\begin{aligned}
&(u,v)=\left(e^{x} \cos{\left(y \right)},e^{x} \sin{\left(y \right)}\right),\\
&D=\begin{pmatrix}e^{x} \cos{\left(y \right)} & - e^{x} \sin{\left(y \right)}\\e^{x} \sin{\left(y \right)} & e^{x} \cos{\left(y \right)}\end{pmatrix},\\
&\det D=e^{2 x}
\end{aligned}
$$

入力と出力の座標値が違うので、出力側の計量は像の点で評価する。二本の像ベクトルの内積を取って入力側へ戻す。

$$
\begin{aligned}
&g_{\mathrm{target}}\circ F=\begin{pmatrix}1 & 0\\0 & 1\end{pmatrix},\\
&F^*g=D^\top(g_{\mathrm{target}}\circ F)D=\begin{pmatrix}e^{2 x} & 0\\0 & e^{2 x}\end{pmatrix}
\end{aligned}
$$

二つの座標速度 $(1,0)$ と $(1,1)$ を例に取ると、引き戻した行列からそのなす角を計算できる。微分の像の長さと内積を使う計算と一致する。

$$
\cos\angle((1,0),(1,1))=\frac{\sqrt{2}}{2},\qquad \det(F^*g)=e^{4 x}
$$

この写像は $y$ を $2\pi$ 増やすと同じ像を取る。したがって領域の角の幅を絞ると座標になる。曲率零という条件は局所の結論を与え、定義域全体で単射になるかは別に確かめる。
:::

## 曲面論の基本定理（Bonnet の定理）

曲線論では曲率と捩率から標構を積分した。曲面の場合は二方向に積分でき、結果が順序によらないことが必要となる。その条件がちょうど前節の二式である。

:::lemma[二方向の線形方程式の可積分条件]{#lem-cs-flat-system}
滑らかな正方行列 $M_1,M_2$ が

$$
\partial_1 M_2-\partial_2M_1+M_1M_2-M_2M_1=0
$$

を満たすとする。長方形の近傍で $\partial_i F=FM_i$ は初期値 $F(p)$ を与えると一意に解ける。
:::

:::proof
第一座標方向で初期点を通る直線上に $\partial_1F=FM_1$ を解き、次に各点から第二方向へ $\partial_2F=FM_2$ を解く。解は初期値とパラメータに滑らかに依存する。$D=\partial_1F-FM_1$ と置き、第二方向に微分する。順序交換と可積分条件により $\partial_2D=DM_2$ であり、初期直線で $D=0$。従って一意性により長方形全体で $D=0$ である。別の解も同じ二方向の初期値問題を満たすので一致する。
:::

:::theorem[曲面論の基本定理]{#thm-cs-bonnet-surface}
平面の単連結な領域 $U$ 上で滑らかな正定値対称行列 $g$ と滑らかな対称行列 $h$ が :eqref[eq-cs-gauss-codazzi] を満たすとする。すると第一基本形式 $g$、選んだ単位法線に対する第二基本形式 $h$ を持つはめ込み $X:U\to\R^3$ が存在する。初期点の位置と向き付き標構を指定すると一意である。ここでは座標に対する正法線を採用するので、残る自由度は向きを保つ剛体運動である。
:::

:::proof
$F=(X_1,X_2,N)$ と考えて未知の行列を取り、$\partial_iF=FM_i$ の係数を

$$
M_i=\begin{pmatrix}(\Gamma^l_{ij})_{l,j}&(-A_i^l)_l\\
                     (h_{ij})_j&0\end{pmatrix}
$$

と定める。$C=\partial_1M_2-\partial_2M_1+[M_1,M_2]$ をブロックに分ける。左上は $R^l_{j12}-A_1^l h_{2j}+A_2^l h_{1j}$、左下は $\partial_1h_{2j}-\partial_2h_{1j}+h_{1m}\Gamma^m_{2j}-h_{2m}\Gamma^m_{1j}$ だから、Gauss と Codazzi により零である。残りも零である。実際、$G_0=\diag(g,1)$ とすると :eqref[eq-cs-christoffel] と $g_{kl}A_i^l=h_{ik}$ から $\partial_i G_0=M_i^\top G_0+G_0 M_i$。その混合偏微分を引くと $C^\top G_0+G_0 C=0$。左二ブロックが零なら、この反対称性が右上と右下も零にする。

:ref[lem-cs-flat-system] で局所解を得る。任意の経路を有限個の座標小長方形に分けて解を延長できる。経路を端点を固定して連続変形するとき、変形のパラメータ正方形を細かく分割し、各小片の像が一つの長方形に入るようにする。各小片の周りで局所解が一意だから、延長結果は変形で変わらない。単連結性により初期点からの経路に依存せず、大域解 $F$ を得る。

初期値を $F^\top F=G_0$ かつ $\det F>0$ となるように取る。$F^\top F$ も $G_0$ と同じ線形方程式を満たすので、任意の経路に沿う一意性により全域で一致する。従って列は一次独立であり、第三列は最初の二列に対する正単位法線である。第一列を $P$、第二列を $Q$ と書くと、$\Gamma^k_{ij}$ と $h_{ij}$ の対称性により $\partial_2P=\partial_1Q$ である。よって $P\,\dd x^1+Q\,\dd x^2$ は閉じたベクトル値一形式となり、単連結性から積分して $X$ を得る。微分が $P,Q$ なので階数 $2$、計量は $g$、二階微分の法線成分は $h$ である。一意性は $F$ の一意性と、最後の積分の加法定数だけが自由であることによる。
:::

:::example[平面と円柱を作る二組のデータ]
$g=\dd u^2+\dd v^2$、$h=0$ は平面で実現する。$h=-\tfrac12\dd u^2$ でも $\det h=0$ と成分の一定性から両適合条件が成り立ち、$X=(2\cos(u/2),2\sin(u/2),v)$ がこれを実現する。領域が全平面なら後者ははめ込みだが単射でない。
:::

:::example[一方向だけ曲率が変わるデータ]
$g=\dd u^2+\dd v^2$、$h=(1+u)\dd u^2$ は両適合条件を満たす。実際 $\det h=0$、$\Gamma=0$ であり、唯一変化する成分は $u$ だけに依存するため Codazzi の左右が一致する。このデータを空間で構成するには、単位速さ平面曲線 $\alpha$ を符号付き曲率 $-(1+u)$ から :ref[thm-cs-plane-fundamental] で求め、$X=(\alpha(u),v)$ とする。$\alpha'$ は長さ $1$ で高さ方向と直交し、正法線は $\alpha'$ を時計回りに直角回転した方向であり、反時計回りの法線で測った符号付き曲率の負が $h_{11}=1+u$ になる。従って二つの形式が与えたものと一致する。$u=-1$ では形作用素が零になっても、計量が正定値なので表示は正則である。
:::

### 単連結性を外した場合の二つの周期

Bonnet の構成でまず積分するのは標構の行列であり、その後で位置を積分する。多重連結な領域では、閉路を一周した標構の戻りと、位置の積分の戻りを別に調べる必要がある。標構が一周で一致しても、位置に一定の並進が残る可能性がある。

:::example[平坦円柱をその計量のまま平面に置けない]
抽象領域 $(\R/2\pi\Z)\times\R$ に $g=\dd u^2+\dd v^2$ と $h=0$ を指定する。二つの適合条件は各点で成立する。被覆の平面で構成すると $X=(u,v,0)$ を得る。標構は定値なのに $u\mapsto u+2\pi$ で位置は $(2\pi,0,0)$ だけずれる。従ってこの $X$ は商へ降りない。また被覆での任意の実現はこの $X$ に向きを保つ剛体運動を施したものなので、非零な並進がなくなる実現も存在しない。局所では問題ないデータに、大域の閉路が新しい制約を課している。

同じ計量で $h=-\dd u^2$ に変えると $X=(\cos u,\sin u,v)$ が実現となり、$u$ 方向の一周で接平面も位置も戻る。これにより、基礎となる計量だけでなく、指定した第二基本形式も閉路の条件に関与することが分かる。
:::


:::example[回転が戻る条件も第二基本形式で変わる]
周期 $2\pi$ の座標 $u$ と実数値の座標 $v$ を持つ抽象円柱で $g=\dd u^2+\dd v^2$ とする。$h=-k\,\dd u^2$、$k>0$ とすると局所的な実現は
$$
X(u,v)=\left(\frac1k\cos(ku),\frac1k\sin(ku),v\right)
$$
である。第一基本形式は速度の長さがともに一で直交するため指定の式になり、円周方向の曲率は外向きで $-k$ となる。$u$ の一周で標構も位置も戻るためには $2\pi k\in2\pi\Z$ が必要十分である。従って $k=1/2$ は局所では実現できるが商全体へ降りず、$k=2$ は全体へ降りる。後者では $u$ を $\pi$ 増やすだけで同じ像を取るため、はめ込みにはなるが単射性は失われる。二つの基本形式の適合、大域での周期、実現の単射性はこの順に別の条件を加える。
:::

:::example[円柱の座標を変えて基本定理を確認する]
半径 $2$ の円柱を $X(u,v)=(2\cos u,2\sin u,v+u)$ で表す。外向き法線を取る。

第一基本形式から内積の行列を読み、逆行列を出す。逆行列を使うのは、内積で求めた下付き成分から接ベクトルの座標成分へ戻すためである。

$$
\begin{aligned}
&g=\begin{pmatrix}5 & 1\\1 & 1\end{pmatrix},\\
&\det(g)=4,\\
&g^{-1}=\begin{pmatrix}\frac{1}{4} & - \frac{1}{4}\\- \frac{1}{4} & \frac{5}{4}\end{pmatrix}
\end{aligned}
$$

逆行列と計量の一階微分を Christoffel の公式に代入する。零でない成分を次に示す。下の二添字を交換した成分は同じである。

$$
\Gamma^k_{ij}=0\quad\text{for all }i,j,k
$$

さらに二方向の微分の交換子を計算する。$R(\partial_u,\partial_v)\partial_v$ の座標列に左から $g$ の第一行を掛ければ、$\partial_u$ との内積になる。

$$
[R(\partial_u,\partial_v)\partial_v]=\begin{pmatrix}0\\0\end{pmatrix},\qquad K=\frac{g_{1m}R^m_{212}}{\det g}=0
$$

候補の第二基本形式は、まず行列式を比較して Gauss 方程式を確認する。さらに各 $j$ について右端の微分と接続の項を足し、Codazzi の差が零であることを確かめる。

$$
\begin{aligned}
&h=\begin{pmatrix}-2 & 0\\0 & 0\end{pmatrix},\\
&\det h=0,\\
&\partial_u h_{2j}-\partial_v h_{1j}+h_{1m}\Gamma^m_{2j}-h_{2m}\Gamma^m_{1j}=0
\end{aligned}
$$

非直交な座標でも $\det h=0$ で、成分が一定なので Codazzi も成り立つ。直接求めると $X_{uu}=(-2\cos u,-2\sin u,0)$ だから、法線成分は $-2$ で一致する。
:::

:::example[同じ球面を二つの縮尺で調べる]
半径 $3$ の球の立体射影から得る $g=36(1+u^2+v^2)^{-2}I$ に、$h=-g/3$ を組み合わせる。

第一基本形式から内積の行列を読み、逆行列を出す。逆行列を使うのは、内積で求めた下付き成分から接ベクトルの座標成分へ戻すためである。

$$
\begin{aligned}
&g=\begin{pmatrix}\frac{36}{\left(u^{2} + v^{2} + 1\right)^{2}} & 0\\0 & \frac{36}{\left(u^{2} + v^{2} + 1\right)^{2}}\end{pmatrix},\\
&\det(g)=\frac{1296}{\left(u^{2} + v^{2} + 1\right)^{4}},\\
&g^{-1}=\begin{pmatrix}\frac{\left(u^{2} + v^{2} + 1\right)^{2}}{36} & 0\\0 & \frac{\left(u^{2} + v^{2} + 1\right)^{2}}{36}\end{pmatrix}
\end{aligned}
$$

逆行列と計量の一階微分を Christoffel の公式に代入する。零でない成分を次に示す。下の二添字を交換した成分は同じである。

$$
\begin{aligned}
&\Gamma^{u}_{uu}=- \frac{2 u}{u^{2} + v^{2} + 1},\\
&\Gamma^{u}_{uv}=- \frac{2 v}{u^{2} + v^{2} + 1},\\
&\Gamma^{u}_{vv}=\frac{2 u}{u^{2} + v^{2} + 1},\\
&\Gamma^{v}_{uu}=\frac{2 v}{u^{2} + v^{2} + 1},\\
&\Gamma^{v}_{uv}=- \frac{2 u}{u^{2} + v^{2} + 1},\\
&\Gamma^{v}_{vv}=- \frac{2 v}{u^{2} + v^{2} + 1}
\end{aligned}
$$

さらに二方向の微分の交換子を計算する。$R(\partial_u,\partial_v)\partial_v$ の座標列に左から $g$ の第一行を掛ければ、$\partial_u$ との内積になる。

$$
\begin{aligned}
&[R(\partial_u,\partial_v)\partial_v]=\begin{pmatrix}\frac{4}{\left(u^{2} + v^{2} + 1\right)^{2}}\\0\end{pmatrix},\\
&K=\frac{g_{1m}R^m_{212}}{\det g}=\frac{1}{9}
\end{aligned}
$$

候補の第二基本形式は、まず行列式を比較して Gauss 方程式を確認する。さらに各 $j$ について右端の微分と接続の項を足し、Codazzi の差が零であることを確かめる。

$$
\begin{aligned}
&h=\begin{pmatrix}- \frac{12}{\left(u^{2} + v^{2} + 1\right)^{2}} & 0\\0 & - \frac{12}{\left(u^{2} + v^{2} + 1\right)^{2}}\end{pmatrix},\\
&\det h=\frac{144}{\left(u^{2} + v^{2} + 1\right)^{4}},\\
&\partial_u h_{2j}-\partial_v h_{1j}+h_{1m}\Gamma^m_{2j}-h_{2m}\Gamma^m_{1j}=0
\end{aligned}
$$

単連結な座標域で与えた二形式は、基本定理により初期値を合わせれば剛体運動を除き一意に実現される。第三列の法線が正になるよう、立体射影の座標の順序を選ぶ。
:::

## 等温座標と直交座標による Gauss 曲率の公式

平面の座標軸が曲面上で直交するなら、正規化して正規直交標構を作れる。この標構が接平面内でどれだけ回転するかを一形式で書くと、曲率の計算が短くなる。

:::proposition[直交座標と等温座標の曲率]{#prop-cs-orthogonal-k}
$g=a^2\,\dd u^2+b^2\,\dd v^2$、$a,b>0$ なら

$$
K=-\frac1{ab}\left(\partial_u(b_u/a)+\partial_v(a_v/b)\right).
\label{eq-cs-orthogonal-k}
$$

特に $g=e^{2\phi}(\dd u^2+\dd v^2)$ なら $K=-e^{-2\phi}(\phi_{uu}+\phi_{vv})$ である。
:::

:::proof
正規直交基底を $e_1=X_u/a$、$e_2=X_v/b$ とし、$\omega_i=\inner{\partial_i e_1}{e_2}$ と置く。計量を微分すると :eqref[eq-cs-christoffel] から $\omega_u=-a_v/b,\ \omega_v=b_u/a$ である。法線成分を $\eta_i=\inner{\partial_i e_1}N$、$\xi_i=\inner{\partial_i e_2}N$ とすると

$$
\partial_i e_1=\omega_i e_2+\eta_iN,\quad
\partial_i e_2=-\omega_i e_1+\xi_iN,\quad
\partial_iN=-\eta_i e_1-\xi_i e_2.
$$

$\partial_u\partial_v e_1-\partial_v\partial_u e_1=0$ の $e_2$ 成分から $\partial_u\omega_v-\partial_v\omega_u=\eta_v\xi_u-\eta_u\xi_v$。ここで $\eta_u=ah(e_1,e_1)$、$\xi_v=bh(e_2,e_2)$、$\eta_v=bh(e_1,e_2)$、$\xi_u=ah(e_1,e_2)$。従って右辺は $-abK$ であり、公式を得る。等温座標では $a=b=e^\phi$ なので $(b_u/a,a_v/b)=(\phi_u,\phi_v)$ となる。
:::

同じ関係は $\omega=\omega_u\dd u+\omega_v\dd v$ に対して $d\omega=-K\,\dd A$ と書ける。向きを保って基底を角 $\theta$ だけ回すと、新しい $\omega$ は $\omega+\dd\theta$ になる。この式は回転した基底を微分して第二成分を取ると出る。二度の混合偏微分が相殺するため $d\omega$ は基底の選び方に依存しない。

この外微分の式は、はめ込みを使わずにも導ける。正向き正規直交標構に対して、内積の微分から $\nabla_Y e_1=\omega(Y)e_2$、$\nabla_Y e_2=-\omega(Y)e_1$ である。曲率作用素を $e_1$ に作用させると、$e_1$ 方向の二つの積 $\omega(Y)\omega(Z)$ は相殺し、$e_2$ 方向には $Y(\omega(Z))-Z(\omega(Y))-\omega([Y,Z])=d\omega(Y,Z)$ が残る。$Y=e_1,Z=e_2$ とするとこの成分は $-K$ なので、$d\omega=-K\,\dd A$ が従う。従って正定値計量だけを指定した直交座標でも曲率の公式が使える。

:::example[零曲率のとき直交単位標構を積分する]
計量の曲率が小円板上で恒等的に零なら、正向き正規直交標構の接続形式は $d\omega=0$ を満たす。円板上で $d\theta=-\omega$ を積分し、標構を角 $\theta$ だけ回すと接続形式が零になる。回した標構の二つの双対一形式を $\eta^1,\eta^2$ とする。座標ベクトルの接微分は順序を交換できるので、$\partial_u(\eta^i(\partial_v))-\partial_v(\eta^i(\partial_u))=0$ である。実際、零の接続形式に対して標構の微分が消え、残る差は二つの座標ベクトルの混合微分の差だけである。従って $d\eta^i=0$ であり、これらを積分して $\eta^i=dy^i$ と書ける。二形式は点ごとに独立なので $(y^1,y^2)$ は局所座標となり、計量は $(\dd y^1)^2+(\dd y^2)^2$ になる。これで零曲率から平面への局所等長写像を、等温座標をあらかじめ求めずに構成できる。
:::

:::example[懸垂面と Enneper 曲面]
懸垂面で $\phi=\log a+\log\cosh v$ なら $\Delta\phi=\operatorname{sech}^2v$ なので $K=-1/(a^2\cosh^4v)$ となる。Enneper 曲面では $\phi=\log(1+u^2+v^2)$ から $\Delta\phi=4/(1+u^2+v^2)^2$、従って $K=-4/(1+u^2+v^2)^4$ である。後者の原点の値は第二基本形式による計算とも一致する。
:::

:::example[半平面の負曲率計量]
$v>0$ で $g=(\dd u^2+\dd v^2)/v^2$ とすると $\phi=-\log v$ であり、$\Delta\phi=1/v^2$。従ってこの計量の曲率は $-1$ である。ここでの公式は計量の微分だけを用いるため、空間へのはめ込みを前提にせず曲率を定められる。
:::

### 斜交座標での曲率計算

$F$ が零でない場合でも、正規直交な余標構を取れば曲率を一形式の外微分から計算できる。座標を実際に直交座標へ作り直す必要はない。

:::proposition[斜交した計量から求める曲率]{#prop-cs-oblique-k}
$g=E\,\dd u^2+2F\,\dd u\,\dd v+G\,\dd v^2$ とし

$$
a=\sqrt E,\quad b=F/\sqrt E,\quad c=\sqrt{EG-F^2}/\sqrt E
$$

と置く。$P=(b_u-a_v)/c$、$Q=(c_u+bP)/a$ に対し

$$
K=-\frac{\partial_u Q-\partial_v P}{ac}.
\label{eq-cs-oblique-k}
$$

:::

:::proof
$\eta^1=a\,\dd u+b\,\dd v$、$\eta^2=c\,\dd v$ と置けば $g=(\eta^1)^2+(\eta^2)^2$ であり、これは正向き正規直交余標構である。双対標構を $e_1,e_2$、$\omega(Y)=\inner{\nabla_Ye_1}{e_2}$ とする。二つの座標ベクトルに対する接微分を交換すると、座標ベクトルの $e_1$ 成分について $d\eta^1=\omega\wedge\eta^2$、$e_2$ 成分について $d\eta^2=-\omega\wedge\eta^1$ となる。例えば第一式は

$$
\partial_u\eta^1_v-\partial_v\eta^1_u-\omega_u\eta^2_v+\omega_v\eta^2_u=0
$$

であり、$X_{uv}$ の $e_1$ 成分を二通りに微分した差に等しい。$\omega=P\,\dd u+Q\,\dd v$ と書き、この二式に代入すると $Pc=b_u-a_v$ と $Qa-Pb=c_u$ を得る。従って $d\omega=(Q_u-P_v)\,\dd u\wedge\dd v$、一方 $\dd A=ac\,\dd u\wedge\dd v$ だから $d\omega=-K\,\dd A$ より結論を得る。
:::

:::example[斜交した球面座標]
球面の余緯度 $x$ と経度 $y$ を $x=u,\ y=u+v$ と書き換えると、$E=1+\sin^2u,F=\sin^2u,G=\sin^2u$ となる。この斜交計量には、正向き余標構 $\vartheta^1=\dd u$、$\vartheta^2=\sin u(\dd u+\dd v)$ も取れる。この余標構では $d\vartheta^1=0$、$d\vartheta^2=\cos u\,\dd u\wedge\dd v$ なので $\omega=\cos u(\dd u+\dd v)$。従って $d\omega=-\sin u\,\dd u\wedge\dd v$、面積密度は $\sin u$ だから $K=1$ となる。斜交座標の $F$ が非零でも、正規直交余標構は簡単な形で選べる場合がある。
:::

:::example[平坦な共形計量の局所展開]
共形計量 $e^{2\phi}(\dd u^2+\dd v^2)$ が平坦なら $\Delta\phi=0$。単連結な小領域で $\psi_u=-\phi_v$、$\psi_v=\phi_u$ を満たす調和共役 $\psi$ が取れる。右辺の一形式が閉じる条件は $\phi_{uu}+\phi_{vv}=0$ だからである。そこで $F(z)=\int e^{\phi+i\psi}\,\dd z$ と置くと $F'(z)\ne0$ で、Euclid 計量の引き戻しは $|F'|^2|\dd z|^2=e^{2\phi}|\dd z|^2$。例えば $\phi=2u$ なら $F(z)=e^{2z}/2$ を選べる。曲率が零であるという微分条件を、局所等長写像の具体的な構成に変えられる。
:::

### 曲率を共形変形で変える

:::proposition[共形変形前後の Gauss 曲率]
$\widehat g=e^{2\psi}g$ とすると、元の曲率 $K$ と変形後の曲率 $\widehat K$ には
$$
e^{2\psi}\widehat K=K-\Delta_g\psi
$$
の関係がある。
:::

:::proof
:ref[thm-cs-isothermal-existence] により局所で等温座標を選び、$g=e^{2\phi}(\dd u^2+\dd v^2)$ とする。同じ座標で $\widehat g=e^{2(\phi+\psi)}(\dd u^2+\dd v^2)$ である。等温座標の曲率公式を両者に使うと $e^{2\psi}\widehat K=-e^{-2\phi}(\partial_u^2+\partial_v^2)(\phi+\psi)$ となる。この第一項は $K$、第二項は等温座標の Laplace–Beltrami 作用素 $-\Delta_g\psi$ である。
:::

:::example[平坦計量の一部を正曲率にする]
Euclid 平面の計量を $e^{2\psi}$ 倍し、$\psi=\log2-\log(1+u^2+v^2)$ とする。平面の Laplace 作用素を取ると $\Delta\psi=-4/(1+u^2+v^2)^2$。従って $\widehat K=-e^{-2\psi}\Delta\psi=1$ である。局所的には球面の計量が復元されるが、原点から無限遠へ向かう径方向の長さは $\int_0^\infty 2/(1+u^2)\,\dd u=\pi$ である。従ってこの表示における無限遠は有限距離にあり、平面全体を球の全体と見なすことはできない。
:::


### 独立した一変数の係数で平坦性を確認する

:::example[両方の座標を伸び縮みさせる平坦計量]
$(u,v)$ 上の $g=e^{2u}\dd u^2+(1+v^2)^2\dd v^2$ では、$a=e^u$ は $u$ のみ、$b=1+v^2$ は $v$ のみに依存する。直交座標の公式に現れるのは $a_v,b_u$ なので両方が零となる。実際、変数
$$
x=e^u,\qquad y=v+v^3/3
$$
を取れば、$\dd x=e^u\dd u$、$\dd y=(1+v^2)\dd v$ であり、$g=\dd x^2+\dd y^2$ になる。従って曲率だけでなく展開先も計算できる。$x$ の値域は $x>0$、$y$ の値域は全実数である。$u$ が任意の実数を取っても、曲面はこの計量では平面全体でなく半平面と等長である。

例えば $(u,v)=(0,0)$ と $(\log 2,1)$ は、展開図では $(1,0),(2,4/3)$ となる。それらを結ぶ線分は $x>0$ に含まれるので距離は $\sqrt{1+(4/3)^2}=5/3$ である。ここでは座標の形式を直したことで、微分幾何の局所計算から距離の大域計算へ進める。
:::

:::example[一方向だけ曲率の符号が変わる計量]
$g=\dd u^2+(3+\cos u)^2\dd v^2$ では $a=1$、$b=3+\cos u$ である。$b_{uu}=-\cos u$ だから曲率は $\cos u/(3+\cos u)$。$u=0,\pi/2,\pi$ での値は $1/4,0,-1/2$ である。例えば $u=\pi/2$ の点では $K=0$ でも、小さな近傍全体を平面に展開できない。座標のうち一方向の係数だけが変化するときも、その変化の二階微分が曲率を生む。
:::

:::example[半径係数の二階微分で曲率を判断する]
$u>0$ で半径係数 $f(u)=1+u^2$ を選び、点で負曲率になることを計算する。

径方向の長さは座標差そのものである。一方、角方向に長さ一だけ進むには、その点の半径係数の逆数だけ座標を動かす。単位方向と非零になり得る接続係数を並べる。

$$
\begin{aligned}
&g=\dd u^2+\left(u^{2} + 1\right)^2\dd v^2,\\
&e_1=\partial_u,\\
&e_2=\frac{1}{u^{2} + 1}\partial_v,\\
&\Gamma^u_{vv}=- 2 u \left(u^{2} + 1\right),\\
&\Gamma^v_{uv}=\Gamma^v_{vu}=\frac{2 u}{u^{2} + 1}
\end{aligned}
$$

曲率を計算するには角方向の長さの変化をもう一度微分する。係数の一次微分の二乗が相殺することは、接続係数で計算した成分を使っても確かめられる。

$$
\begin{aligned}
&R^u{}_{vuv}=\partial_u\Gamma^u_{vv}-\Gamma^u_{vv}\Gamma^v_{uv}=- 2 \left(u^{2} + 1\right),\\
&K=-\frac{\partial_u^2 f}{f}=- \frac{2}{u^{2} + 1}
\end{aligned}
$$

この計量の $u$ 一定の線を $v$ の増加向きにたどると、加速度の単位径方向成分は $-f'/f$ である。この接成分が零になる場合の意味は [共変微分と測地線](/math/curves-surfaces/geodesics/) で扱う。

$$
\nabla_{e_2}e_2=- \frac{2 u}{u^{2} + 1}e_1,\qquad K\big|_{u=1}=-1
$$

係数の一次微分は $u>0$ で正だが、曲率の符号を決める式には二階微分が入る。計量が角座標に依存しないことだけから零曲率は従わない。
:::

:::example[二次の共形因子から曲率を求める]
$g=e^{2(u^2-v^2)}(\dd u^2+\dd v^2)$ を調べる。符号が異なる二つの二階微分から曲率を求める。

第一基本形式から内積の行列を読み、逆行列を出す。逆行列を使うのは、内積で求めた下付き成分から接ベクトルの座標成分へ戻すためである。

$$
\begin{aligned}
&g=\begin{pmatrix}e^{2 u^{2}} e^{- 2 v^{2}} & 0\\0 & e^{2 u^{2}} e^{- 2 v^{2}}\end{pmatrix},\\
&\det(g)=e^{4 u^{2}} e^{- 4 v^{2}},\\
&g^{-1}=\begin{pmatrix}e^{- 2 u^{2}} e^{2 v^{2}} & 0\\0 & e^{- 2 u^{2}} e^{2 v^{2}}\end{pmatrix}
\end{aligned}
$$

逆行列と計量の一階微分を Christoffel の公式に代入する。零でない成分を次に示す。下の二添字を交換した成分は同じである。

$$
\begin{aligned}
&\Gamma^{u}_{uu}=2 u,\\
&\Gamma^{u}_{uv}=- 2 v,\\
&\Gamma^{u}_{vv}=- 2 u,\\
&\Gamma^{v}_{uu}=2 v,\\
&\Gamma^{v}_{uv}=2 u,\\
&\Gamma^{v}_{vv}=- 2 v
\end{aligned}
$$

さらに二方向の微分の交換子を計算する。$R(\partial_u,\partial_v)\partial_v$ の座標列に左から $g$ の第一行を掛ければ、$\partial_u$ との内積になる。

$$
[R(\partial_u,\partial_v)\partial_v]=\begin{pmatrix}0\\0\end{pmatrix},\qquad K=\frac{g_{1m}R^m_{212}}{\det g}=0
$$

$u^2-v^2$ は調和なので、曲率の零は共形座標の式からも分かる。逆に $u^2+v^2$ に替えると Laplace 作用素の値が $4$ になり、曲率は負になる。
:::

:::example[半径と角を斜めに混ぜた円柱状の計量]
$g=\dd u^2+f(u)^2(\dd v+\dd u)^2$、$f=2+\cos u$ を計算する。

第一基本形式から内積の行列を読み、逆行列を出す。逆行列を使うのは、内積で求めた下付き成分から接ベクトルの座標成分へ戻すためである。

$$
\begin{aligned}
&g=\begin{pmatrix}\cos^{2}{\left(u \right)} + 4 \cos{\left(u \right)} + 5 & \left(\cos{\left(u \right)} + 2\right)^{2}\\\left(\cos{\left(u \right)} + 2\right)^{2} & \left(\cos{\left(u \right)} + 2\right)^{2}\end{pmatrix},\\
&\det(g)=\left(\cos{\left(u \right)} + 2\right)^{2},\\
&g^{-1}=\begin{pmatrix}1 & -1\\-1 & \frac{\cos^{2}{\left(u \right)} + 4 \cos{\left(u \right)} + 5}{\left(\cos{\left(u \right)} + 2\right)^{2}}\end{pmatrix}
\end{aligned}
$$

逆行列と計量の一階微分を Christoffel の公式に代入する。零でない成分を次に示す。下の二添字を交換した成分は同じである。

$$
\begin{aligned}
&\Gamma^{u}_{uu}=\left(\cos{\left(u \right)} + 2\right) \sin{\left(u \right)},\\
&\Gamma^{u}_{uv}=\left(\cos{\left(u \right)} + 2\right) \sin{\left(u \right)},\\
&\Gamma^{u}_{vv}=\left(\cos{\left(u \right)} + 2\right) \sin{\left(u \right)},\\
&\Gamma^{v}_{uu}=- \frac{\left(- \sin^{2}{\left(u \right)} + 4 \cos{\left(u \right)} + 7\right) \sin{\left(u \right)}}{\cos{\left(u \right)} + 2},\\
&\Gamma^{v}_{uv}=- \frac{\left(\cos^{2}{\left(u \right)} + 4 \cos{\left(u \right)} + 5\right) \sin{\left(u \right)}}{\cos{\left(u \right)} + 2},\\
&\Gamma^{v}_{vv}=- \left(\cos{\left(u \right)} + 2\right) \sin{\left(u \right)}
\end{aligned}
$$

さらに二方向の微分の交換子を計算する。$R(\partial_u,\partial_v)\partial_v$ の座標列に左から $g$ の第一行を掛ければ、$\partial_u$ との内積になる。

$$
\begin{aligned}
&[R(\partial_u,\partial_v)\partial_v]=\begin{pmatrix}\left(\cos{\left(u \right)} + 2\right) \cos{\left(u \right)}\\\sin^{2}{\left(u \right)} - 2 \cos{\left(u \right)} - 1\end{pmatrix},\\
&K=\frac{g_{1m}R^m_{212}}{\det g}=\frac{\cos{\left(u \right)}}{\cos{\left(u \right)} + 2}
\end{aligned}
$$

変数 $V=u+v$ に直せば直交計量となる。その式で $-f''/f$ を求めても同じ値になるため、交差項がある座標での曲率計算を独立に確かめられる。
:::

## 演習問題

::::exercise{level=basic}
$g=\dd u^2+e^{2u}\dd v^2$ の零でない Christoffel 記号を求めよ。

:::solution
$\partial_u g_{22}=2e^{2u}$ が唯一の非零の計量微分である。従って $\Gamma^1_{22}=-e^{2u}$、$\Gamma^2_{12}=\Gamma^2_{21}=1$。他の成分は :eqref[eq-cs-christoffel] に零を代入して得る。
:::
::::

::::exercise{level=basic}
$g=\dd u^2+\sin^2u\,\dd v^2$、$0<u<\pi$ の $K$ を求めよ。

:::solution
直交座標公式で $a=1$、$b=\sin u$。第二の微分項は零であり、第一は $b_{uu}=-\sin u$ なので $K=1$ となる。
:::
::::

::::exercise{level=basic}
$g=e^{2u}(\dd u^2+\dd v^2)$ の曲率を求め、Euclid 平面への局所等長写像を一つ与えよ。

:::solution
$\phi=u$ は調和なので $K=0$。$f(u,v)=(e^u\cos v,e^u\sin v)$ の二つの偏微分は直交し、ともに長さ $e^u$ である。従って平面計量をこの写像で引き戻すと指定の計量であり、階数も $2$ である。
:::
::::

::::exercise{level=basic}
$g=\dd u^2+\dd v^2$、$h=2\,\dd u^2+3\,\dd v^2$ は実現可能か。

:::solution
第一基本形式から $K=0$ となるが、第二基本形式の行列式は $6$ で第一の行列式は $1$。Gauss 方程式が同じ $K$ に $6$ を要求するので、この組は実現できない。
:::
::::

::::exercise{level=advanced}
$g$ が Euclid 計量で、$h$ が定数 $a$ と滑らかな関数 $q(u,v)$ を用いて $h=a\,\dd u^2+q\,\dd v^2$ と書かれるとする。$a\ne0$ のとき両適合条件を満たす $q$ を決定せよ。

:::solution
Gauss 方程式は $a q=0$ だから $q=0$ を要求する。このとき $h$ は定数行列なので Codazzi も成立する。逆に $q=0$ は例えば半径 $1/|a|$ の円柱で、法線の符号を合わせて実現できる。
:::
::::

::::exercise{level=advanced}
計量 $g=e^{2\phi}(\dd u^2+\dd v^2)$ を定数 $c^2>0$ 倍したとき $K$ がどう変わるか示せ。

:::solution
新しい共形指数は $\phi+\log c$ であり、Laplacian は変わらない。係数 $e^{-2\phi}$ だけが $c^{-2}$ 倍になるので曲率も $c^{-2}$ 倍になる。これは球の半径を $c$ 倍したときの $1/R^2$ の変化と一致する。
:::
::::

::::exercise{level=advanced}
$g=\dd u^2+(2+u)^2\dd v^2$、$u>-2$ で $K$ を求め、局所的な展開を与えよ。

:::solution
$a=1$、$b=2+u$ で $b_{uu}=a_v=0$ だから $K=0$。$(x,y)=((2+u)\cos v,(2+u)\sin v)$ とすれば極座標の計量が指定の式になる。半径が正なのでこれは局所微分同相であり、角変数が周期的であるため必要に応じて範囲を制限する。
:::
::::

::::exercise{level=advanced}
$g=\dd u^2+\dd v^2$、$h=0$ が与えられたとする。基本定理の構成を直接解け。

:::solution
全ての $\Gamma$ と $A$ が零なので $M_i=0$。従って $F$ は定数であり、最初の二列を直交単位ベクトル $P,Q$ とすれば $X(u,v)=p+uP+vQ$。第三列は $P\times Q$ で、この平面の法線も一定である。
:::
::::

::::exercise{level=advanced}
$|\!u\!|<1$ で計量 $g=\dd u^2+(2+u^2)^2(\dd u+\dd v)^2$ が与えられた。この計量の曲率を求め、原点の近くで回転面として実現できる座標を一組与えよ。
:::solution
$w=u+v$ と置くと $g=\dd u^2+f(u)^2\dd w^2$、$f=2+u^2$。直交座標の公式で $K=-f''/f=-2/(2+u^2)$ となる。実現する回転面を $X(u,w)=(f(u)\cos w,f(u)\sin w,z(u))$ の形で探す。$X_u$ と $X_w$ は直交し、$|X_w|=f$、$|X_u|^2=f'^2+z'^2$ だから、$z'(u)=\sqrt{1-4u^2}$ とすればよい。これは $|u|<1/2$ で正であり、
$$
z(u)=\tfrac12u\sqrt{1-4u^2}+\tfrac14\arcsin(2u)
$$
はその原始関数である。$X$ の外積は長さ $f>0$ を持つので正則である。元の座標には $w=u+v$ を戻せばよい。この実現方法が要求したのは $|f'|<1$ であり、曲率計算そのものはこれより広い領域で有効である。
:::
::::

::::exercise{level=advanced}
平面の領域 $v>0$ に $g=(\dd u^2+\dd v^2)/v^2$、$h=v^{-1}(\dd u\dd v+\dd v\dd u)$ を指定する。Gauss 方程式と Codazzi 方程式のどちらが成立するかを成分まで調べよ。ここで対称積の表示は $h_{12}=h_{21}=1/v$ を意味する。
:::solution
計量の曲率は $K=-1$、行列式は $v^{-4}$。候補の第二基本形式の行列式は $-v^{-2}$ なので、Gauss 方程式は $v=1$ の点で一致するものの、領域全体では成立しない。Codazzi 方程式も成分を計算する。非零な接続係数は $\Gamma^u_{uv}=\Gamma^u_{vu}=-1/v$、$\Gamma^v_{uu}=1/v$、$\Gamma^v_{vv}=-1/v$ である。差
$$
C_j=\partial_u h_{2j}-\partial_v h_{1j}
   + h_{1m}\Gamma^m_{2j}-h_{2m}\Gamma^m_{1j}
$$
で $j=1$ の微分項は零で、接続項も零である。$j=2$ の微分項は $v^{-2}$、二つの接続項は $-v^{-2}$ と $+v^{-2}$ で相殺する。従って $C_2=v^{-2}\ne0$ である。単に曲率が負である条件を満たしても、そのことから基本形式の組を実現できるわけではない。
:::
::::

## 確認テスト

::::quiz
局所等長写像が常に保つ量を選べ。

- [x] Gauss 曲率
- [ ] 平均曲率
- [ ] 空間の単位法線

:::explanation
Gauss 曲率は同じ座標で比較した計量とその微分から求まる。平面と円柱は局所等長だが法線の点ごとの変化が異なり、平均曲率も異なる。
:::
::::

::::quiz{answer="-1" tol="1e-8"}
$v>0$ 上の計量 $(\dd u^2+\dd v^2)/v^2$ の Gauss 曲率を求めよ。

:::explanation
$\phi=-\log v$ なので $\Delta\phi=1/v^2$。等温座標の式 $K=-e^{-2\phi}\Delta\phi$ に代入して $K=-v^2/v^2=-1$ を得る。
:::
::::

::::quiz
曲面論の基本定理で二つの基本形式の適合条件となるものを選べ。

- [x] Gauss 方程式
- [x] Codazzi–Mainardi 方程式
- [ ] 平均曲率が零であること

:::explanation
標構の二つの微分が可換であることを要求すると、接成分と法線成分から Gauss と Codazzi がそれぞれ出る。この条件が線形方程式を積分できる条件でもある。平均曲率は形作用素の跡の半分であり、零と仮定しない。
:::
::::

::::quiz{answer="1" tol="1e-8"}
$g=\dd r^2+r^2\dd\theta^2$ に対し $r=1$ で $\Gamma^\theta_{r\theta}$ を求めよ。

:::explanation
Christoffel 記号の公式の該当成分は $\tfrac12r^{-2}\partial_r(r^2)=1/r$ である。従って $r=1$ では値は $1$。
:::
::::

---
title: "曲線座標とテンソル"
description: "座標に合わせてベクトル場の微分を計算し、計量・共変微分・向きの変換則を導く。"
notation:
  - symbol: 'q^i'
    meaning: "第 i 座標"
  - symbol: 'h_i'
    meaning: "直交曲線座標の尺度因子"
  - symbol: 'g_{ij}'
    meaning: "計量の座標成分"
  - symbol: '\Gamma^i_{jk}'
    meaning: "接続の座標成分"
  - symbol: '\varepsilon_{ijk}'
    meaning: "向きを表す体積形式の成分"
---

## この章の目標

- 座標線の接ベクトルから尺度因子と体積要素を求める。
- 積分で表される流束と循環から曲線座標の微分公式を導く。
- 座標成分、双対基底および共変微分を使って座標に依存しない量を計算する。
- 反転での符号を調べ、交代記号・体積形式・密度を区別する。

## 直交曲線座標と尺度因子

### 座標写像と小さな変位

円筒の壁を調べるとき、位置を半径と角度で指定すると壁の方程式が簡単になる。その代わり、角度を同じだけ変えたときにも移動距離は半径によって変わる。微分演算にこの距離を取り入れるため、座標線の接ベクトルを使う。

:::definition[尺度因子]{#def-mp-scale-factors}
開集合 $U\subset\R^n$ からの滑らかな座標写像 $\vb r(q^1,\ldots,q^n)\in\R^n$ が局所微分同相であるとする。$\partial_i=\partial/\partial q^i$、$\vb e_i=\partial_i\vb r$ と書く。異なる $\vb e_i$ が直交するとき、$h_i=\abs{\vb e_i}>0$ を :term[尺度因子]{en="scale factor" yomi="しゃくどいんし"}、$\hat{\vb e}_i=\vb e_i/h_i$ を対応する単位ベクトルという。
:::

全微分から $\dd\vb r=\sum_i h_i\hat{\vb e}_i\dd q^i$ であり、内積を取ると $\dd s^2=\sum_i h_i^2(\dd q^i)^2$ となる。接ベクトルの内積を $g_{ij}=\vb e_i\cdot\vb e_j$ と書くと、直交座標では $g_{ij}$ は対角で、対角成分は $h_i^2$ である。以下の例ではこの行列を $g$ と書き、行列式と区別するため $\det g$ と明記する。微小直方体の平行六面体近似の体積は接ベクトルの行列式の絶対値である。互いに直交するため、Jacobian の絶対値は $\prod_i h_i$ となる。三次元では $q^i$ を固定した面の面素は $h_jh_k\,\dd q^j\dd q^k$ である。$q^i$ が角度なら $h_i$ の単位は $\mathrm m$、$q^i$ が長さなら無次元であり、$h_i\dd q^i$ はいずれも長さである。

:::example[円筒座標]
$\vb r=(\rho\cos\phi,\rho\sin\phi,z)$ を $\rho>0$ で使う。接ベクトルは $(\cos\phi,\sin\phi,0)$、$\rho(-\sin\phi,\cos\phi,0)$、$(0,0,1)$ なので $h=(1,\rho,1)$ である。たとえば $\rho=2\,\mathrm m$ で $\Delta\phi=0.1$ の円弧は $0.2\,\mathrm m$、$q^1$ を固定した壁の面素は $\rho\,\dd\phi\dd z$ である。
:::

:::example[球座標]
$\vb r=r(\sin\theta\cos\phi,\sin\theta\sin\phi,\cos\theta)$ とし、$\theta$ は正の $z$ 軸からの角度とする。$\partial_r\vb r=\vb r/r$、$\partial_\theta\vb r=r(\cos\theta\cos\phi,\cos\theta\sin\phi,-\sin\theta)$、$\partial_\phi\vb r=r\sin\theta(-\sin\phi,\cos\phi,0)$ だから $h=(1,r,r\sin\theta)$ である。$r=2\,\mathrm m$ の球面の北半分では面積が $4\int_0^{\pi/2}\sin\theta\,\dd\theta\int_0^{2\pi}\dd\phi=8\pi\,\mathrm{m^2}$ となる。
:::

:::example[放物線座標]
平面で $x=(u^2-v^2)/2$、$y=uv$ とすると、接ベクトル $(u,v)$ と $(-v,u)$ は直交する。両方の尺度因子は $\sqrt{u^2+v^2}$ である。$(u,v)=(1,2)$ における面積の倍率は $5$、長さの倍率は $\sqrt5$ である。
:::

円筒の軸や球座標の極では接ベクトルの一部が零になる。その点ではこの座標写像を微分同相として使えない。座標式に $1/\rho$ が出る場合、元の場が滑らかかどうかは直交座標へ戻して判定する。

### 円筒座標以外の分離に向いた座標

平面に二つの焦点があるときは楕円座標 $x=a\cosh u\cos v$、$y=a\sinh u\sin v$（$u>0$、$a>0$）が適している。例えば $v=0$ を固定して $u$ を変えると $x>a$ の半直線を得る。一方、$u$ を固定して $v$ を一周させると、長半径 $a\cosh u$、短半径 $a\sinh u$ の楕円を描く。接ベクトルの内積は
$$
a^2[\sinh u\cos v(-\cosh u\sin v)
+\cosh u\sin v\,\sinh u\cos v]=0
$$
である。二つの長さの二乗はともに $a^2(\sinh^2 u+\sin^2 v)$ になる。従って $u$ が一定の楕円の内側と外側で問題を分けるとき、境界は座標の一つを固定して与えられる。$u\to0$ で楕円は焦点間の線分に近づくので、その線分をまたいで場を定義するには別の座標または元の Cartesian 座標で確認する。

:::example[平面の極座標を立体へ上げる]
放物線座標に $z=w$ を足す。$(u,v)=(1,2)$ で、辺の座標幅がそれぞれ $\Delta u=\Delta v=\Delta w=0.01$ の箱の体積の主項は $(1^2+2^2)(0.01)^3=5\cdot10^{-6}$ である。$u,v$ に長さの平方根の次元を持たせた場合は $h_u,h_v$ にも同じ次元がある。そのため $h_uh_v\,\dd u\dd v\dd w$ は長さの三乗の次元になる。
:::

### 座標線と座標特異点

曲線座標で微小な長さを扱うには、各座標を単独で動かしたときの接線を比較する。単位正方形または単位立方体を写す線形近似は Jacobian である。独立な接線がそろうところでのみ、その像は面積または体積をもつ。例えば球座標の極では方位角をどれだけ変えても点は動かない。そこでは $h_\phi=0$ であり、極付近の滑らかな場を調べるには直交座標などに戻る。

直交座標では $g_{ij}=0$ が $i\ne j$ で成り立つ。したがって任意の二点間を座標ごとに離したときの一次近似は $\dd s^2=\sum_i h_i^2(\dd q^i)^2$ で与えられる。曲線全体の長さを求めるときは、その曲線の速度をこの式に代入してから積分する。座標値の差に一点での尺度因子を掛ける計算は、座標差が微小な場合に限る。

:::example[極座標の円周方向]
円周方向の角度差を弧長へ換算するため、平面の極座標を調べる。

座標の取り方を写像で指定し、一つずつ座標を動かして接ベクトルを得る。それらを列にした行列を $J$ とする。

$$
\begin{aligned}\vb X&=\begin{pmatrix}r \cos{\left(\phi \right)}\\r \sin{\left(\phi \right)}\end{pmatrix}\\ J&=\begin{pmatrix}\cos{\left(\phi \right)} & - r \sin{\left(\phi \right)}\\\sin{\left(\phi \right)} & r \cos{\left(\phi \right)}\end{pmatrix}\end{aligned}
$$

接ベクトルの内積を並べると、変位とその長さの関係が分かる。逆計量は双対基底から長さを読むときに使う。

$$
\begin{aligned}g=J^{\mathsf T}J&=\begin{pmatrix}1 & 0\\0 & r^{2}\end{pmatrix}\\ g^{-1}&=\begin{pmatrix}1 & 0\\0 & \frac{1}{r^{2}}\end{pmatrix}\end{aligned}
$$

体積係数は行列式の絶対値で与えられる。符号を残した値は正向きか逆向きかを区別する。

$$
\begin{aligned}\det J&=r \left(\sin^{2}{\left(\phi \right)} + \cos^{2}{\left(\phi \right)}\right)\\ \sqrt{\det g}&=r \left(\sin^{2}{\left(\phi \right)} + \cos^{2}{\left(\phi \right)}\right)\end{aligned}
$$

具体的に座標点 $\left( 2, \  \frac{\pi}{4}\right)$ で、位置と単位変位の計量を求める。

$$
\begin{aligned}\vb X&=\begin{pmatrix}\sqrt{2}\\\sqrt{2}\end{pmatrix}\\ g&=\begin{pmatrix}1 & 0\\0 & 4\end{pmatrix}\end{aligned}
$$

行列が対角なら、各成分の平方根が尺度因子である。非対角成分があれば、二方向へ同時に進むときに交差項が寄与する。
:::

:::example[二次元の放物線座標]
放物線座標を $(x,y)=((u^2-v^2)/2,uv)$、$u>0$ で取る。

座標の取り方を写像で指定し、一つずつ座標を動かして接ベクトルを得る。それらを列にした行列を $J$ とする。

$$
\begin{aligned}\vb X&=\begin{pmatrix}\frac{u^{2}}{2} - \frac{v^{2}}{2}\\u v\end{pmatrix}\\ J&=\begin{pmatrix}u & - v\\v & u\end{pmatrix}\end{aligned}
$$

接ベクトルの内積を並べると、変位とその長さの関係が分かる。逆計量は双対基底から長さを読むときに使う。

$$
\begin{aligned}g=J^{\mathsf T}J&=\begin{pmatrix}u^{2} + v^{2} & 0\\0 & u^{2} + v^{2}\end{pmatrix}\\ g^{-1}&=\begin{pmatrix}\frac{1}{u^{2} + v^{2}} & 0\\0 & \frac{1}{u^{2} + v^{2}}\end{pmatrix}\end{aligned}
$$

体積係数は行列式の絶対値で与えられる。符号を残した値は正向きか逆向きかを区別する。

$$
\begin{aligned}\det J&=u^{2} + v^{2}\\ \sqrt{\det g}&=u^{2} + v^{2}\end{aligned}
$$

具体的に座標点 $\left( 2, \  1\right)$ で、位置と単位変位の計量を求める。

$$
\begin{aligned}\vb X&=\begin{pmatrix}\frac{3}{2}\\2\end{pmatrix}\\ g&=\begin{pmatrix}5 & 0\\0 & 5\end{pmatrix}\end{aligned}
$$

行列が対角なら、各成分の平方根が尺度因子である。非対角成分があれば、二方向へ同時に進むときに交差項が寄与する。
:::

:::example[対数半径を使う座標]
対数半径は相対的な半径の増分を一定に刻むために使える。

座標の取り方を写像で指定し、一つずつ座標を動かして接ベクトルを得る。それらを列にした行列を $J$ とする。

$$
\begin{aligned}\vb X&=\begin{pmatrix}e^{u} \cos{\left(v \right)}\\e^{u} \sin{\left(v \right)}\end{pmatrix}\\ J&=\begin{pmatrix}e^{u} \cos{\left(v \right)} & - e^{u} \sin{\left(v \right)}\\e^{u} \sin{\left(v \right)} & e^{u} \cos{\left(v \right)}\end{pmatrix}\end{aligned}
$$

接ベクトルの内積を並べると、変位とその長さの関係が分かる。逆計量は双対基底から長さを読むときに使う。

$$
\begin{aligned}g=J^{\mathsf T}J&=\begin{pmatrix}e^{2 u} & 0\\0 & e^{2 u}\end{pmatrix}\\ g^{-1}&=\begin{pmatrix}e^{- 2 u} & 0\\0 & e^{- 2 u}\end{pmatrix}\end{aligned}
$$

体積係数は行列式の絶対値で与えられる。符号を残した値は正向きか逆向きかを区別する。

$$
\begin{aligned}\det J&=\left(\sin^{2}{\left(v \right)} + \cos^{2}{\left(v \right)}\right) e^{2 u}\\ \sqrt{\det g}&=\left(\sin^{2}{\left(v \right)} + \cos^{2}{\left(v \right)}\right) e^{2 u}\end{aligned}
$$

具体的に座標点 $\left( \log{\left(3 \right)}, \  \frac{\pi}{3}\right)$ で、位置と単位変位の計量を求める。

$$
\begin{aligned}\vb X&=\begin{pmatrix}\frac{3}{2}\\\frac{3 \sqrt{3}}{2}\end{pmatrix}\\ g&=\begin{pmatrix}9 & 0\\0 & 9\end{pmatrix}\end{aligned}
$$

行列が対角なら、各成分の平方根が尺度因子である。非対角成分があれば、二方向へ同時に進むときに交差項が寄与する。
:::

## 勾配・発散・回転・Laplacian の曲線座標表示

### 微小な箱と面を使った導出

$\vb A=\sum_i A_{\hat i}\hat{\vb e}_i$ のハット付き添字は単位基底に関する成分を表す。勾配、発散、回転は :ref[def-anl-vector-operators] の演算である。以下では $(q^1,q^2,q^3)$ が右手系となる領域に限る。

:::proposition[直交座標の微分演算]{#prop-mp-orthogonal-operators}
$f$ は $C^2$、$\vb A$ は $C^1$ とし、$J=h_1h_2h_3$ とする。$(i,j,k)$ が $(1,2,3)$ の巡回順列のとき、次が成り立つ。

$$
\begin{aligned}
(\grad f)_{\hat i}&=\frac{\partial_i f}{h_i},\\
\divg\vb A&=\frac1J\sum_i\partial_i(J A_{\hat i}/h_i),\\
(\rot\vb A)_{\hat i}&=\frac{\partial_j(h_k A_{\hat k})-\partial_k(h_j A_{\hat j})}{h_jh_k},\\
\Delta f&=\frac1J\sum_i\partial_i(Jh_i^{-2}\partial_i f).
\end{aligned}
\label{eq-mp-orthogonal-operators}
$$
:::

:::proof
$\dd f=\grad f\cdot\dd\vb r$ に $\dd\vb r=\sum h_i\hat{\vb e}_i\dd q^i$ を入れて係数を比較すると勾配を得る。座標で辺長が $\eta$ の箱をとる。$i$ に垂直な二面での外向き流束の和は、積分の差を微分で表すと $\eta^3\partial_i(h_jh_k A_{\hat i})+o(\eta^3)$ である。箱の体積は $J\eta^3+o(\eta^3)$ なので、:ref[thm-anl-divergence] を使い割り算して極限を取ると発散になる。

$q^i$ を固定した面を、境界で $+j,+k,-j,-k$ の順に回る。$+j$ と $-j$ の辺の寄与の和は $-\eta^2\partial_k(h_j A_{\hat j})+o(\eta^2)$ であり、残りの二辺から $+\eta^2\partial_j(h_k A_{\hat k})+o(\eta^2)$ が出る。:ref[thm-anl-stokes] によって、和を面積 $h_jh_k\eta^2+o(\eta^2)$ で割った極限は回転の第 $i$ 成分である。最後に $\Delta f=\divg\grad f$ に勾配と発散の式を順に使えばよい。
:::

:::example[円柱内の熱伝導]
角度と高さに依存しない温度 $T(\rho)$ では $\Delta T=(\rho T')'/\rho$ である。熱伝導率が一定で、熱源のない定常状態では $\rho T'=C$ なので、$0<a<\rho<b$ として
$T=T_a+(T_b-T_a)\log(\rho/a)/\log(b/a)$ を得る。熱伝導率 $\kappa$ の単位は $\mathrm{W\,m^{-1}K^{-1}}$ で、長さ $L$ の円柱を外向きに通る電力換算の熱流は $-2\pi\kappa L C$（$\mathrm W$）である。$b=2a$ なら温度が両端の平均となる位置は $\sqrt2 a$ である。
:::

:::example[回転する剛体]
角速度 $\Omega$ の剛体回転の速度 $\vb A=\Omega\rho\hat{\vb e}_\phi$ を考える。$\divg\vb A=0$、$\rot\vb A=\rho^{-1}\partial_\rho(\Omega\rho^2)\hat{\vb e}_z=2\Omega\hat{\vb e}_z$ となる。$\Omega=3\,\mathrm{s^{-1}}$ では渦度は $6\,\mathrm{s^{-1}}$ である。$\rho$ を一つ落として微分すると、この値を取り損なう。
:::

:::example[球面上の角度依存性]
球座標では $\Delta f=r^{-2}\partial_r(r^2\partial_r f)+(r^2\sin\theta)^{-1}\partial_\theta(\sin\theta\,\partial_\theta f)+(r\sin\theta)^{-2}\partial_\phi^2 f$ である。$f=r\cos\theta$ の動径部分は $2\cos\theta/r$、角度部分は $-2\cos\theta/r$ なので和は零である。実際、この関数は直交座標では $z$ である。
:::

### 回転を線積分から求める

面積素片を向き付きで取ると、回転の法線成分は面を囲む微小ループの循環を面積で割った値になる。例えば $(q^2,q^3)$ 面の正辺では線積分の被積分関数がそれぞれ $h_2A_{\hat2}$ と $h_3A_{\hat3}$ である。対向する二辺の差を整理すると、$\partial_2(h_3A_{\hat3})-\partial_3(h_2A_{\hat2})$ が得られる。分母はその面の面積係数 $h_2h_3$ である。この計算を $(3,1),(1,2)$ に巡回置換すれば三成分がそろう。

同じ座標で閉じた微小箱を取ると、発散は向かい合う面を通る流束の差から得られる。この場合は線素ではなく面積を掛けるため、例えば第 $1$ 成分の微分の内側に $h_2h_3$ が入る。面の外向き法線を負側で反転することが、二面の差を作る。両方の方法の積分形は :ref[thm-anl-stokes] と :ref[thm-anl-divergence] に対応する。

次の例では接ベクトル $\vb e_i$ に関する成分 $A^i=A_{\hat i}/h_i$ も使う。勾配の成分は同じ換算で $(\grad f)^i=h_i^{-2}\partial_i f$ となる。これを $\nabla^if$ と書く。また体積因子を $\sqrt g=\prod_i h_i$ と書く。ここでの $g$ は $g_{ij}$ の行列式である。平面の場合も、同じ式を第三座標 $z$ に依存しない場へ適用すれば二成分で計算できる。

:::example[円筒の軸依存性]
角速度が $2$ の回転場に、半径に比例する動径成分と高さの二乗に比例する軸方向成分を加える。

スカラーとベクトルを指定する。列の成分は座標基底に関する。

$$
\begin{aligned} g&=\begin{pmatrix}1 & 0 & 0\\0 & r^{2} & 0\\0 & 0 & 1\end{pmatrix}\\ f&=r^{2} z\\ (A^i)&=\begin{pmatrix}r\\2\\z^{2}\end{pmatrix}\end{aligned}
$$

まずスカラーの勾配を求める。偏導関数を並べた列は共変成分であるため、逆行列で添字を上げる。

$$
\begin{aligned}(\partial_i f)&=\begin{pmatrix}2 r z\\0\\r^{2}\end{pmatrix}\\ (\nabla^if)&=\begin{pmatrix}2 r z\\0\\r^{2}\end{pmatrix}\end{aligned}
$$

体積因子を各成分に掛け、各成分を同じ添字の座標で微分してから足す。

$$
\begin{aligned}\sqrt{g}&=r\\ (\sqrt{g}A^i)&=\begin{pmatrix}r^{2}\\2 r\\r z^{2}\end{pmatrix}\\ \operatorname{div}A&=\frac{2 r z + 2 r}{r}=2 z + 2\end{aligned}
$$

最後に、ベクトルをスカラーの勾配で置き換えると Laplacian が求まる。

$$
\begin{aligned}(\sqrt{g}\nabla^if)&=\begin{pmatrix}2 r^{2} z\\0\\r^{3}\end{pmatrix}\\ \Delta f&=4 z\end{aligned}
$$

発散はスカラーなので、別の座標で計算して同じ点に戻れば値が一致する。これに対して勾配の成分を比較するには基底変換がいる。
:::

:::example[極座標での非軸対称モード]
極座標で二方向の微分がともに必要になる例として、直交座標の $x^2-y^2$ を調べる。

スカラーとベクトルを指定する。列の成分は座標基底に関する。

$$
\begin{aligned} g&=\begin{pmatrix}1 & 0\\0 & r^{2}\end{pmatrix}\\ f&=r^{2} \cos{\left(2 \phi \right)}\\ (A^i)&=\begin{pmatrix}r \cos{\left(\phi \right)}\\\sin{\left(\phi \right)}\end{pmatrix}\end{aligned}
$$

まずスカラーの勾配を求める。偏導関数を並べた列は共変成分であるため、逆行列で添字を上げる。

$$
\begin{aligned}(\partial_i f)&=\begin{pmatrix}2 r \cos{\left(2 \phi \right)}\\- 2 r^{2} \sin{\left(2 \phi \right)}\end{pmatrix}\\ (\nabla^if)&=\begin{pmatrix}2 r \cos{\left(2 \phi \right)}\\- 2 \sin{\left(2 \phi \right)}\end{pmatrix}\end{aligned}
$$

体積因子を各成分に掛け、各成分を同じ添字の座標で微分してから足す。

$$
\begin{aligned}\sqrt{g}&=r\\ (\sqrt{g}A^i)&=\begin{pmatrix}r^{2} \cos{\left(\phi \right)}\\r \sin{\left(\phi \right)}\end{pmatrix}\\ \operatorname{div}A&=\frac{3 r \cos{\left(\phi \right)}}{r}=3 \cos{\left(\phi \right)}\end{aligned}
$$

最後に、ベクトルをスカラーの勾配で置き換えると Laplacian が求まる。

$$
\begin{aligned}(\sqrt{g}\nabla^if)&=\begin{pmatrix}2 r^{2} \cos{\left(2 \phi \right)}\\- 2 r \sin{\left(2 \phi \right)}\end{pmatrix}\\ \Delta f&=0\end{aligned}
$$

発散はスカラーなので、別の座標で計算して同じ点に戻れば値が一致する。これに対して勾配の成分を比較するには基底変換がいる。
:::

:::example[対数半径での発散]
距離の単位が各点で指数的に変化する座標を使う。

スカラーとベクトルを指定する。列の成分は座標基底に関する。

$$
\begin{aligned} g&=\begin{pmatrix}e^{2 u} & 0\\0 & e^{2 u}\end{pmatrix}\\ f&=u v\\ (A^i)&=\begin{pmatrix}u\\v\end{pmatrix}\end{aligned}
$$

まずスカラーの勾配を求める。偏導関数を並べた列は共変成分であるため、逆行列で添字を上げる。

$$
\begin{aligned}(\partial_i f)&=\begin{pmatrix}v\\u\end{pmatrix}\\ (\nabla^if)&=\begin{pmatrix}v e^{- 2 u}\\u e^{- 2 u}\end{pmatrix}\end{aligned}
$$

体積因子を各成分に掛け、各成分を同じ添字の座標で微分してから足す。

$$
\begin{aligned}\sqrt{g}&=e^{2 u}\\ (\sqrt{g}A^i)&=\begin{pmatrix}u e^{2 u}\\v e^{2 u}\end{pmatrix}\\ \operatorname{div}A&=\frac{2 u e^{2 u} + 2 e^{2 u}}{e^{2 u}}=2 u + 2\end{aligned}
$$

最後に、ベクトルをスカラーの勾配で置き換えると Laplacian が求まる。

$$
\begin{aligned}(\sqrt{g}\nabla^if)&=\begin{pmatrix}v\\u\end{pmatrix}\\ \Delta f&=0\end{aligned}
$$

発散はスカラーなので、別の座標で計算して同じ点に戻れば値が一致する。これに対して勾配の成分を比較するには基底変換がいる。
:::

## 一般曲線座標と計量テンソル

### 長さと角度を一つの行列にまとめる

座標線が斜めなら、線素には $\dd q^i\dd q^j$ の混合項が出る。斜交格子上で長さや流束を計算するには、接ベクトルの長さに加えて、その内積も保持する。

:::definition[計量]{#def-mp-metric-components}
座標基底の :term[計量テンソル]{en="metric tensor" yomi="けいりょうてんそる"} の成分を $g_{ij}=\vb e_i\cdot\vb e_j$ とする。逆行列を $(g^{ij})$、行列式を $g=\det(g_{ij})$ と書く。双対基底 $\dd q^i$ は $\dd q^i(\vb e_j)=\delta^i_j$ で定まる。
:::

この節以降、同じ項で一度上、一度下に出る添字は総和を取る。$\dd s^2=g_{ij}\dd q^i\dd q^j$ で、実ベクトル $\vb v=v^i\vb e_i$ の長さの二乗は $g_{ij}v^iv^j$ である。$B$ を接ベクトルを列に並べた行列とすると $g_{ij}=(B^\top B)_{ij}$ であり、$g=(\det B)^2>0$ となる。このため $n$ 次元の体積要素は $\sqrt g\,\dd^n q$ である。

:::proposition[勾配と Laplacian]{#prop-mp-metric-laplacian}
任意の滑らかな座標で、$C^2$ 級の関数 $f$ の勾配と Laplacian は次式で与えられる。

$$
(\grad f)^i=g^{ij}\partial_j f,\qquad
\Delta f=\frac1{\sqrt g}\partial_i(\sqrt g\,g^{ij}\partial_j f).
\label{eq-mp-metric-laplacian}
$$
:::

:::proof
内積 $\grad f\cdot\vb e_j=\partial_j f$ に逆行列を作用させると勾配を得る。台が座標領域内に収まる滑らかな関数 $\psi$ を掛けて、直交座標で部分積分すると $\int\psi\divg\vb A\,\dd V=-\int \vb A\cdot\grad\psi\,\dd V$ である。右辺を座標表示すると $-\int A^i\partial_i\psi\,\sqrt g\,\dd^n q$ であり、もう一度部分積分すると $\int\psi\partial_i(\sqrt g A^i)\,\dd^n q$ となる。差が連続で全ての $\psi$ と積分して零になるため、差が一箇所でも正または負ならその近くに非負の $\psi$ を選んで矛盾する。よって $\divg\vb A=g^{-1/2}\partial_i(\sqrt g A^i)$ であり、勾配を代入して結論を得る。
:::

:::example[せん断格子]
$x=u+2v$、$y=v$、$z=w$ とする。$u,v$ 部分の計量は $\begin{pmatrix}1&2\\2&5\end{pmatrix}$、逆行列は $\begin{pmatrix}5&-2\\-2&1\end{pmatrix}$ である。$\sqrt g=1$ なので $\Delta=5\partial_u^2-4\partial_u\partial_v+\partial_v^2+\partial_w^2$ となる。$u=x-2y$ から直接 $u^2$ を微分しても $\Delta u^2=10$ を得る。
:::

平面で計量が $\lambda(u,v)^2(\dd u^2+\dd v^2)$ なら $\sqrt g=\lambda^2$、$g^{ij}=\lambda^{-2}\delta^{ij}$ である。したがって $\Delta=\lambda^{-2}(\partial_u^2+\partial_v^2)$ となる。たとえば放物線座標では $u^2-v^2$ が調和関数になり、これは直交座標の $2x$ と一致する。

:::example[斜交座標の勾配と直交する等位線]
$x=u+2v,y=v$ で $f(u,v)=uv$ とする。逆計量を掛けると $\grad f=(5v-2u)\vb e_u+(u-2v)\vb e_v$ になる。$\vb e_u=(1,0)$、$\vb e_v=(2,1)$ なので Cartesian 成分は $(v,u-2v)$ である。直接 $f=(x-2y)y$ を微分して $(y,x-4y)$ を得ても同じである。例えば $(u,v)=(2,1)$ では勾配は $(1,0)$ となり、この点の等位線の接線は鉛直である。
:::

### 双対計量と面の測度

接ベクトルを列に並べた行列を $J$ とする。$\vb e_i$ を基底とするベクトルの第 $i$ 成分を内積で取り出すには、双対基底を計量でベクトルに直せばよい。$\vb e^i=g^{ij}\vb e_j$ とおけば $\vb e^i\cdot\vb e_k=\delta^i_k$ になる。この式から、曲面が $q^1=\mathrm{const.}$ で表されるとき法線は $\vb e^1$ に平行である。一般に $\vb e_1$ はその曲面に直交しない。

一方、面積係数には接線だけを使う。三次元の $q^1=\mathrm{const.}$ で接線が $\vb e_2,\vb e_3$ なら、平行四辺形の面積の二乗は $g_{22}g_{33}-g_{23}^2$ である。余因子公式からこの式は $g\,g^{11}$ に等しい。法線の単位化に $\sqrt{g^{11}}$ が現れるため、面積を掛けた流束の被積分関数は正側で $\sqrt{g}A^1$ になる。これを全座標面に適用して発散を求めれば、座標が斜交していても同じ表示を使える。

:::example[平面の斜交格子]
座標線の傾きが一致しない平行四辺形を作る。

座標の取り方を写像で指定し、一つずつ座標を動かして接ベクトルを得る。それらを列にした行列を $J$ とする。

$$
\begin{aligned}\vb X&=\begin{pmatrix}2 u + v\\3 v\end{pmatrix}\\ J&=\begin{pmatrix}2 & 1\\0 & 3\end{pmatrix}\end{aligned}
$$

接ベクトルの内積を並べると、変位とその長さの関係が分かる。逆計量は双対基底から長さを読むときに使う。

$$
\begin{aligned}g=J^{\mathsf T}J&=\begin{pmatrix}4 & 2\\2 & 10\end{pmatrix}\\ g^{-1}&=\begin{pmatrix}\frac{5}{18} & - \frac{1}{18}\\- \frac{1}{18} & \frac{1}{9}\end{pmatrix}\end{aligned}
$$

体積係数は行列式の絶対値で与えられる。符号を残した値は正向きか逆向きかを区別する。

$$
\begin{aligned}\det J&=6\\ \sqrt{\det g}&=6\end{aligned}
$$

具体的に座標点 $\left( 1, \  2\right)$ で、位置と単位変位の計量を求める。

$$
\begin{aligned}\vb X&=\begin{pmatrix}4\\6\end{pmatrix}\\ g&=\begin{pmatrix}4 & 2\\2 & 10\end{pmatrix}\end{aligned}
$$

行列が対角なら、各成分の平方根が尺度因子である。非対角成分があれば、二方向へ同時に進むときに交差項が寄与する。
:::

:::example[曲がる斜交格子]
せん断が位置によって変わる場合に逆計量を調べる。

座標の取り方を写像で指定し、一つずつ座標を動かして接ベクトルを得る。それらを列にした行列を $J$ とする。

$$
\begin{aligned}\vb X&=\begin{pmatrix}u\\u^{2} + v\end{pmatrix}\\ J&=\begin{pmatrix}1 & 0\\2 u & 1\end{pmatrix}\end{aligned}
$$

接ベクトルの内積を並べると、変位とその長さの関係が分かる。逆計量は双対基底から長さを読むときに使う。

$$
\begin{aligned}g=J^{\mathsf T}J&=\begin{pmatrix}4 u^{2} + 1 & 2 u\\2 u & 1\end{pmatrix}\\ g^{-1}&=\begin{pmatrix}1 & - 2 u\\- 2 u & 4 u^{2} + 1\end{pmatrix}\end{aligned}
$$

体積係数は行列式の絶対値で与えられる。符号を残した値は正向きか逆向きかを区別する。

$$
\begin{aligned}\det J&=1\\ \sqrt{\det g}&=1\end{aligned}
$$

具体的に座標点 $\left( 2, \  1\right)$ で、位置と単位変位の計量を求める。

$$
\begin{aligned}\vb X&=\begin{pmatrix}2\\5\end{pmatrix}\\ g&=\begin{pmatrix}17 & 4\\4 & 1\end{pmatrix}\end{aligned}
$$

行列が対角なら、各成分の平方根が尺度因子である。非対角成分があれば、二方向へ同時に進むときに交差項が寄与する。
:::

:::example[斜交格子の勾配と発散]
変位の長さを一定の斜交格子で測る。

スカラーとベクトルを指定する。列の成分は座標基底に関する。

$$
\begin{aligned} g&=\begin{pmatrix}1 & 1\\1 & 5\end{pmatrix}\\ f&=u v\\ (A^i)&=\begin{pmatrix}u + v\\v^{2}\end{pmatrix}\end{aligned}
$$

まずスカラーの勾配を求める。偏導関数を並べた列は共変成分であるため、逆行列で添字を上げる。

$$
\begin{aligned}(\partial_i f)&=\begin{pmatrix}v\\u\end{pmatrix}\\ (\nabla^if)&=\begin{pmatrix}- \frac{u}{4} + \frac{5 v}{4}\\\frac{u}{4} - \frac{v}{4}\end{pmatrix}\end{aligned}
$$

体積因子を各成分に掛け、各成分を同じ添字の座標で微分してから足す。

$$
\begin{aligned}\sqrt{g}&=2\\ (\sqrt{g}A^i)&=\begin{pmatrix}2 u + 2 v\\2 v^{2}\end{pmatrix}\\ \operatorname{div}A&=\frac{4 v + 2}{2}=2 v + 1\end{aligned}
$$

最後に、ベクトルをスカラーの勾配で置き換えると Laplacian が求まる。

$$
\begin{aligned}(\sqrt{g}\nabla^if)&=\begin{pmatrix}- \frac{u}{2} + \frac{5 v}{2}\\\frac{u}{2} - \frac{v}{2}\end{pmatrix}\\ \Delta f&=- \frac{1}{2}\end{aligned}
$$

発散はスカラーなので、別の座標で計算して同じ点に戻れば値が一致する。これに対して勾配の成分を比較するには基底変換がいる。
:::

## テンソルの変換則と添字記法

### 入力の種類と成分の変換

同じ変位を二つの座標で記述すると基底が変わる。変位そのものを保つには、成分を基底と逆向きに変換する必要がある。この考えを多重線形写像に適用すると、添字の数にかかわらず変換則を得られる。

:::definition[テンソルの型]{#def-mp-tensor-type}
接空間を $V$、双対空間を $V^*$ とする。$(r,s)$ 型 :term[テンソル]{en="tensor" yomi="てんそる"} は、$r$ 個の共変ベクトルと $s$ 個のベクトルを受け取って実数を返す多重線形写像である。この空間は :ref[thm-linalg-tensor-basis] により $V^{\otimes r}\otimes(V^*)^{\otimes s}$ と同一視できる。成分 $T^{i_1\cdots i_r}{}_{j_1\cdots j_s}$ は、入力に双対基底と座標基底をこの順で入れた値である。
:::

新座標を $Q^a$ とし、$A^a{}_i=\partial Q^a/\partial q^i$、$B^i{}_a=\partial q^i/\partial Q^a$ とする。連鎖律で $\dd Q^a=A^a{}_i\dd q^i$、$\vb E_a=B^i{}_a\vb e_i$ となる。多重線形性から
$$
T'^{a_1\cdots a_r}{}_{b_1\cdots b_s}
=A^{a_1}{}_{i_1}\cdots A^{a_r}{}_{i_r}
B^{j_1}{}_{b_1}\cdots B^{j_s}{}_{b_s}
T^{i_1\cdots i_r}{}_{j_1\cdots j_s}
\label{eq-mp-tensor-transform}
$$
を得る。上の添字を反変、下の添字を共変という。上下の一組を縮約すると $\sum_a B^j{}_aA^a{}_i=\delta^j_i$ が働くので、残った添字がテンソルの変換則に従う。たとえば $K^i{}_i$ は線形写像 $K$ のトレースである。計量を用いれば $v_i=g_{ij}v^j$、$v^i=g^{ij}v_j$ となる。$A_{\hat i}$ への換算とは区別する。

:::example[一次形式の読み取り]
$x=u+2v$、$y=v$ とし、$\alpha=3\dd x-\dd y$ とする。$\alpha=3\dd u+5\dd v$ である。変位 $\vb v=\partial_v$ は直交座標で $2\partial_x+\partial_y$ だから、$\alpha(\vb v)=5$ はどちらの表示でも同じである。
:::

:::example[線形写像の行列]
直交座標の二方向をそれぞれ $2$ 倍と $3$ 倍に伸ばす $K$ の行列は $\diag(2,3)$ である。上のせん断座標では $B^{-1}KB=\begin{pmatrix}2&-2\\0&3\end{pmatrix}$ となる。トレース $5$ と行列式 $6$ が変わらない一方、成分 $K^u{}_v$ は変化する。
:::

### 双対基底と変換の逆向き

接基底が $\vb E_a=B^i{}_a\vb e_i$ と変わるとき、双対基底は $\dd Q^a=(B^{-1})^a{}_i\dd q^i$ である。等式 $\dd Q^a(\vb E_b)=\delta^a_b$ を検算に使える。接ベクトルの成分が双対基底と同じ行列で変換されるのは、ベクトルの成分をその双対基底で評価するからである。例えば $q^1=2Q^1$ なら $\vb E_1=2\vb e_1$ なので、$\vb e_1$ 自体の $Q^1$ 成分は $1/2$ である。一方 $\dd q^1=2\dd Q^1$ なので、同じ一次形式の $Q^1$ 成分は $2$ になる。内積 $\dd q^1(\vb e_1)=1$ は積 $2\cdot(1/2)$ で保たれる。

テンソルの対称部分と交代部分は、同じ型の添字についてのみ座標に依存せず分けられる。例えば $C_{ij}$ の対称部分は $(C_{ij}+C_{ji})/2$、交代部分は差の半分であり、双線形写像の入力の順を交換して定義できる。混合成分 $T^i{}_j$ に同じ操作をすると、上と下で別の変換行列が掛かる。線形写像の対称性を扱うときは計量を使って $g_{ik}T^k{}_j$ に直し、その共変添字の対称性を調べる。

### 縮約とテンソル積

座標変換を $q\mapsto q'$ とし、その微分を $J$ とする。型 $(p,q)$ のテンソルを変換するとき、各上付きスロットに $J$、各下付きスロットに $J^{-1}$ を挿入する。上下の一組を縮約すると、その対について $\sum_a (J^{-1})^i{}_a J^a{}_j=\delta^i_j$ となる。したがって先に縮約してから変換する計算と、変換してから縮約する計算は一致する。これは座標表示での定義を別の座標系でも使えることを示す。

例えば $V\otimes\alpha$ の作用は、入力 $w$ に対して $w\mapsto V\,\alpha(w)$ である。像は $V$ の張る直線に入る。そのトレースは $\alpha(V)$ であり、$V$ と $\alpha$ がともに非零でもトレースだけ零になることがある。次の例では二つの成分を別々に変換してからこの縮約を確かめる。

:::example[せん断と直交する組]
せん断の前後で同じ写像を表す。

旧座標から新座標への Jacobian と、旧座標でのベクトルと共ベクトルを選ぶ。

$$
\begin{aligned}J&=\begin{pmatrix}1 & 2\\0 & 1\end{pmatrix}\\ V&=\begin{pmatrix}1\\2\end{pmatrix}\\ \alpha&=\begin{pmatrix}2 & -1\end{pmatrix}\end{aligned}
$$

ベクトルは列、共ベクトルは行で書く。行の方は逆行列を右から掛ける。

$$
\begin{aligned}J^{-1}&=\begin{pmatrix}1 & -2\\0 & 1\end{pmatrix}\\ V'&=\begin{pmatrix}5\\2\end{pmatrix}\\ \alpha'&=\begin{pmatrix}2 & -5\end{pmatrix}\end{aligned}
$$

外積として作った写像を二つの座標で計算する。写像を受け取るベクトルにも $J$ が掛かるため、同じ幾何学的操作である。

$$
\begin{aligned}T=V\alpha&=\begin{pmatrix}2 & -1\\4 & -2\end{pmatrix}\\ T'=J T J^{-1}&=\begin{pmatrix}10 & -25\\4 & -10\end{pmatrix}\end{aligned}
$$

最後に縮約を行う。旧成分の縮約と新成分の縮約が同じであることを確認する。

$$
\alpha'V'=\begin{pmatrix}0\end{pmatrix}=\alpha V
$$
:::

:::example[方向ごとに異なる伸縮]
第 1 座標は二倍、第 2 座標は三倍にして成分を比べる。

旧座標から新座標への Jacobian と、旧座標でのベクトルと共ベクトルを選ぶ。

$$
\begin{aligned}J&=\begin{pmatrix}2 & 0\\0 & 3\end{pmatrix}\\ V&=\begin{pmatrix}2\\1\end{pmatrix}\\ \alpha&=\begin{pmatrix}1 & 3\end{pmatrix}\end{aligned}
$$

ベクトルは列、共ベクトルは行で書く。行の方は逆行列を右から掛ける。

$$
\begin{aligned}J^{-1}&=\begin{pmatrix}\frac{1}{2} & 0\\0 & \frac{1}{3}\end{pmatrix}\\ V'&=\begin{pmatrix}4\\3\end{pmatrix}\\ \alpha'&=\begin{pmatrix}\frac{1}{2} & 1\end{pmatrix}\end{aligned}
$$

外積として作った写像を二つの座標で計算する。写像を受け取るベクトルにも $J$ が掛かるため、同じ幾何学的操作である。

$$
\begin{aligned}T=V\alpha&=\begin{pmatrix}2 & 6\\1 & 3\end{pmatrix}\\ T'=J T J^{-1}&=\begin{pmatrix}2 & 4\\\frac{3}{2} & 3\end{pmatrix}\end{aligned}
$$

最後に縮約を行う。旧成分の縮約と新成分の縮約が同じであることを確認する。

$$
\alpha'V'=\begin{pmatrix}5\end{pmatrix}=\alpha V
$$
:::

:::example[座標軸の交換]
二本の座標軸を交換する変換を使う。

旧座標から新座標への Jacobian と、旧座標でのベクトルと共ベクトルを選ぶ。

$$
\begin{aligned}J&=\begin{pmatrix}0 & 1\\1 & 0\end{pmatrix}\\ V&=\begin{pmatrix}-1\\3\end{pmatrix}\\ \alpha&=\begin{pmatrix}2 & 1\end{pmatrix}\end{aligned}
$$

ベクトルは列、共ベクトルは行で書く。行の方は逆行列を右から掛ける。

$$
\begin{aligned}J^{-1}&=\begin{pmatrix}0 & 1\\1 & 0\end{pmatrix}\\ V'&=\begin{pmatrix}3\\-1\end{pmatrix}\\ \alpha'&=\begin{pmatrix}1 & 2\end{pmatrix}\end{aligned}
$$

外積として作った写像を二つの座標で計算する。写像を受け取るベクトルにも $J$ が掛かるため、同じ幾何学的操作である。

$$
\begin{aligned}T=V\alpha&=\begin{pmatrix}-2 & -1\\6 & 3\end{pmatrix}\\ T'=J T J^{-1}&=\begin{pmatrix}3 & 6\\-1 & -2\end{pmatrix}\end{aligned}
$$

最後に縮約を行う。旧成分の縮約と新成分の縮約が同じであることを確認する。

$$
\alpha'V'=\begin{pmatrix}1\end{pmatrix}=\alpha V
$$
:::

## Christoffel 記号と共変微分

### 基底を微分してから成分を読む

一定のベクトルでも、円筒の角度を変えて観測すると成分が変わる。この変化をベクトルの微分に混ぜないよう、基底自身の変化を加算する。

:::definition[共変微分]{#def-mp-covariant-derivative}
Euclid 空間で $\partial_j\vb e_k=\Gamma^i_{jk}\vb e_i$ とし、$\Gamma^i_{jk}$ を :term[クリストッフェル記号]{en="Christoffel symbol" yomi="くりすとっふぇるきごう"} と呼ぶ。$\partial_j(v^k\vb e_k)$ の成分を :term[共変微分]{en="covariant derivative" yomi="きょうへんびぶん"} $\nabla_jv^i=\partial_j v^i+\Gamma^i_{jk}v^k$ とする。
:::

混合偏微分が交換できるので $\Gamma^i_{jk}=\Gamma^i_{kj}$ である。内積の微分により
$\partial_j g_{k\ell}=g_{i\ell}\Gamma^i_{jk}+g_{ki}\Gamma^i_{j\ell}$ である。この式を下添字の $j,k,\ell$ の順を変えて三つ書き、最初の二つの和から三つ目を引くと $2g_{i\ell}\Gamma^i_{jk}$ が残る。逆行列を掛けて
$$
\Gamma^i_{jk}=\tfrac12 g^{i\ell}(\partial_j g_{k\ell}+\partial_k g_{j\ell}-\partial_\ell g_{jk})
\label{eq-mp-christoffel}
$$
を得る。逆にこの式の右辺を使えば下添字に対称で、上の内積の微分式も満たす。このため、曲がった空間でも同じ式で捩れがなく計量を保つ接続を定められる。

共変ベクトル $\alpha$ では、$\partial_j(\alpha_i v^i)$ の積の微分を適用すると $(\nabla_j\alpha)_i=\partial_j\alpha_i-\Gamma^k_{ji}\alpha_k$ が必要になる。一般には上添字それぞれに足す項、下添字それぞれに引く項を一つずつ取り、積の微分で組み合わせる。とくに $\nabla_jg_{k\ell}=0$ であるから、添字の上げ下げは共変微分と交換できる。座標変換後の $\partial_b\vb E_c$ を微分すると、$B^i{}_c$ の微分から $A^a{}_i\partial_b\partial_c q^i$ が加わる。この項があるので $\Gamma$ 単体は式 :eqref[eq-mp-tensor-transform] に従わない。

:::example[円筒座標の加速度]
$\Gamma^\rho_{\phi\phi}=-\rho$、$\Gamma^\phi_{\rho\phi}=\Gamma^\phi_{\phi\rho}=1/\rho$ が非零である。軌道の二階微分を計算すると $\vb a=(\ddot q^i+\Gamma^i_{jk}\dot q^j\dot q^k)\vb e_i$ となる。特に $\rho=2\,\mathrm m$、$\dot\phi=3\,\mathrm{s^{-1}}$ の等速円運動では $\vb a=-18\hat{\vb e}_\rho\,\mathrm{m\,s^{-2}}$ である。回転成分が一定でも加速度が出る理由が基底の微分に現れている。
:::

:::example[一定ベクトルの微分]
$\partial_x=\cos\phi\,\vb e_\rho-\rho^{-1}\sin\phi\,\vb e_\phi$ である。$\nabla_\phi v^\rho=-\sin\phi-\rho(-\rho^{-1}\sin\phi)=0$、$\nabla_\phi v^\phi=-\rho^{-1}\cos\phi+\rho^{-1}\cos\phi=0$ となる。他の方向も定義から零となり、直交座標での結論と一致する。
:::

### 任意の道に沿った微分

軌道 $q(t)$ 上のベクトル場 $\vb V(t)=V^i(t)\vb e_i(q(t))$ の微分は、連鎖律によって $D_tV^i=\dot V^i+\Gamma^i_{jk}\dot q^jV^k$ となる。$D_tV=0$ を平行移動、$D_t\dot q=0$ を測地線の式と呼ぶ。Euclid 空間では前者は Cartesian 成分が一定、後者は加速度が零である。

:::example[極座標で記述した直線]
平面の直線 $(x,y)=(2,t)$ の極座標は $r=\sqrt{4+t^2}$、$\phi=\arctan(t/2)$ である。$\dot r=t/r$、$\ddot r=4/r^3$ と $\dot\phi=2/r^2$ より $\ddot r-r\dot\phi^2=0$ である。また $\ddot\phi=-4t/r^4$ なので $\ddot\phi+2\dot r\dot\phi/r=0$ である。座標の各成分が非線形でも、共変加速度は零である。
:::

### 計量を保つ接続

速度の長さの変化を基底によらず求めるには $\partial_k g_{ij}=\Gamma^\ell_{ki}g_{\ell j}+\Gamma^\ell_{kj}g_{i\ell}$ が必要である。これは $\partial_k(\vb e_i\cdot\vb e_j)$ に積の微分を適用した式である。同時に座標基底では $\partial_j\vb e_k=\partial_k\vb e_j$ なので、下二添字を交換しても $\Gamma$ は変わらない。計量の微分の三式を、$\partial_i g_{jk}+\partial_j g_{ik}-\partial_k g_{ij}$ と組み合わせると、右辺で目的の二項だけが残る。逆計量を掛ければ既出の Christoffel 記号の式になる。

また $\nabla$ は積に対して Leibniz 則に従う。共ベクトルとベクトルを縮約した $\alpha_i V^i$ を微分するとき、上付き添字の接続項と下付き添字の接続項が相殺する。したがって下付きの接続項の符号は上付きと反対になる。これを繰り返せば任意の型で共変微分を計算できる。

:::example[極座標での直線に沿う方向]
直交座標で一定の場 $\vb e_x$ を極座標で書く。

計量と対象のベクトルを座標基底で与える。

$$
\begin{aligned}g&=\begin{pmatrix}1 & 0\\0 & r^{2}\end{pmatrix}\\(V^i)&=\begin{pmatrix}\cos{\left(\phi \right)}\\- \frac{\sin{\left(\phi \right)}}{r}\end{pmatrix}\end{aligned}
$$

計量を微分して Christoffel 記号を作る。上付き添字を固定し、下付き添字を行と列にする。

$$
\begin{aligned}(\Gamma^{1}_{jk})&=\begin{pmatrix}0 & 0\\0 & - r\end{pmatrix}\\(\Gamma^{2}_{jk})&=\begin{pmatrix}0 & \frac{1}{r}\\\frac{1}{r} & 0\end{pmatrix}\end{aligned}
$$

列 $j$ は $q^j$ に沿う共変微分である。成分だけの微分と、基底の微分を並べると違いが見える。

$$
\begin{aligned}(\partial_jV^i)&=\begin{pmatrix}0 & - \sin{\left(\phi \right)}\\\frac{\sin{\left(\phi \right)}}{r^{2}} & - \frac{\cos{\left(\phi \right)}}{r}\end{pmatrix}\\(\nabla_jV^i)&=\begin{pmatrix}0 & 0\\0 & 0\end{pmatrix}\end{aligned}
$$

この行列に $V$ の列を掛けると流線に沿う変化になる。縮約すれば発散になる。

$$
\begin{aligned}\nabla_V V&=\begin{pmatrix}0\\0\end{pmatrix}\\ \nabla_iV^i&=0\end{aligned}
$$
:::

:::example[極座標での回転]
一定角速度 $2$ の剛体回転を座標成分で与える。

計量と対象のベクトルを座標基底で与える。

$$
\begin{aligned}g&=\begin{pmatrix}1 & 0\\0 & r^{2}\end{pmatrix}\\(V^i)&=\begin{pmatrix}0\\2\end{pmatrix}\end{aligned}
$$

計量を微分して Christoffel 記号を作る。上付き添字を固定し、下付き添字を行と列にする。

$$
\begin{aligned}(\Gamma^{1}_{jk})&=\begin{pmatrix}0 & 0\\0 & - r\end{pmatrix}\\(\Gamma^{2}_{jk})&=\begin{pmatrix}0 & \frac{1}{r}\\\frac{1}{r} & 0\end{pmatrix}\end{aligned}
$$

列 $j$ は $q^j$ に沿う共変微分である。成分だけの微分と、基底の微分を並べると違いが見える。

$$
\begin{aligned}(\partial_jV^i)&=\begin{pmatrix}0 & 0\\0 & 0\end{pmatrix}\\(\nabla_jV^i)&=\begin{pmatrix}0 & - 2 r\\\frac{2}{r} & 0\end{pmatrix}\end{aligned}
$$

この行列に $V$ の列を掛けると流線に沿う変化になる。縮約すれば発散になる。

$$
\begin{aligned}\nabla_V V&=\begin{pmatrix}- 4 r\\0\end{pmatrix}\\ \nabla_iV^i&=0\end{aligned}
$$
:::

:::example[非線形せん断での加速度]
写像 $(x,y)=(u,v+u^2)$ を用い、$u$ を一定速度で進める。

計量と対象のベクトルを座標基底で与える。

$$
\begin{aligned}g&=\begin{pmatrix}4 u^{2} + 1 & 2 u\\2 u & 1\end{pmatrix}\\(V^i)&=\begin{pmatrix}1\\0\end{pmatrix}\end{aligned}
$$

計量を微分して Christoffel 記号を作る。上付き添字を固定し、下付き添字を行と列にする。

$$
\begin{aligned}(\Gamma^{1}_{jk})&=\begin{pmatrix}0 & 0\\0 & 0\end{pmatrix}\\(\Gamma^{2}_{jk})&=\begin{pmatrix}2 & 0\\0 & 0\end{pmatrix}\end{aligned}
$$

列 $j$ は $q^j$ に沿う共変微分である。成分だけの微分と、基底の微分を並べると違いが見える。

$$
\begin{aligned}(\partial_jV^i)&=\begin{pmatrix}0 & 0\\0 & 0\end{pmatrix}\\(\nabla_jV^i)&=\begin{pmatrix}0 & 0\\2 & 0\end{pmatrix}\end{aligned}
$$

この行列に $V$ の列を掛けると流線に沿う変化になる。縮約すれば発散になる。

$$
\begin{aligned}\nabla_V V&=\begin{pmatrix}0\\2\end{pmatrix}\\ \nabla_iV^i&=0\end{aligned}
$$
:::

### スカラーの Hessian とベクトルの Laplacian

スカラーの一階微分 $\dd f$ は共変ベクトルである。そのため二階微分をテンソルにするには、下付き添字に対する接続項を含める。
$$
(\nabla\dd f)_{ij}=\partial_i\partial_j f-\Gamma^k_{ij}\partial_k f .
$$
これを :term[ヘッセ行列]{en="Hessian" yomi="へっせぎょうれつ"} の座標成分という。接続の下付き添字の対称性から二つの添字を交換できる。極座標で $f=x=r\cos\phi$ とすると $\partial_\phi^2 f=-r\cos\phi$ であるが、$-\Gamma^r_{\phi\phi}\partial_r f=r\cos\phi$ が相殺し、Hessian の $\phi\phi$ 成分は零になる。直交座標で $x$ が一次式であることと一致する。

体積因子の微分は $\partial_i\sqrt g=\Gamma^k_{ik}\sqrt g$ に従う。この式は、行列式の微分 $\partial_i g=g\,g^{jk}\partial_i g_{jk}$ に計量の積の微分を代入して得られる。そのため
$$
\nabla_iA^i=\partial_i A^i+\Gamma^i_{ik}A^k
=\frac1{\sqrt g}\partial_i(\sqrt g\,A^i).
$$
さらに $\nabla_i g^{jk}=0$ から、$A=\grad f$ の場合には $\nabla_iA^i=g^{ij}(\nabla\dd f)_{ij}$ となる。従って Hessian の計量による縮約は式 :eqref[eq-mp-metric-laplacian] の Laplacian に等しい。

Euclid 空間のベクトル Laplacian は、固定した Cartesian 基底で各成分にスカラーの $\Delta$ を作用させて定義できる。曲線座標では基底も微分される。このため、三次元では $\Delta\vb A=\grad\divg\vb A-\rot\rot\vb A$ を用いると、既出の微分公式だけで計算できる。この恒等式は次節の交代記号の縮約から従う。例えば $\vb A=\rho\,\hat{\vb e}_\rho$ は Cartesian 表示では $(x,y,0)$ なので $\Delta\vb A=0$ である。円筒座標でも発散が定数 $2$、回転が零であるため、同じ値を得る。対して $\rho$ 自体のスカラー Laplacian は $\rho^{-1}$ である。両者の差は $\hat{\vb e}_\rho$ の角度方向の微分に由来する。

### 基底を規格化した場合の接続係数

極座標で正規直交基底を $E_1=\partial_r$、$E_2=r^{-1}\partial_\phi$ と取る。新しい基底では $[E_1,E_2]=-r^{-1}E_2$ となる。共変微分の第二引数について積の微分を使うと
$$
\begin{aligned}
\nabla_{E_1}E_1&=0,& \nabla_{E_1}E_2&=0,\\
\nabla_{E_2}E_1&=r^{-1}E_2,& \nabla_{E_2}E_2&=-r^{-1}E_1
\end{aligned}
$$
となる。例えば $\nabla_{\partial_r}(r^{-1}\partial_\phi)$ では $-r^{-2}\partial_\phi$ と $r^{-1}\Gamma^\phi_{r\phi}\partial_\phi$ が消し合う。続いて $\nabla_{E_2}E_1$ では第二引数に尺度因子が掛からないので $r^{-1}\Gamma^\phi_{\phi r}\partial_\phi=r^{-1}E_2$ が残る。従って正規直交基底での下二つの添字に、座標基底の対称性をそのまま使えない。差 $\nabla_{E_1}E_2-\nabla_{E_2}E_1$ は実際に上で求めた括弧と一致する。

半径 $r=2\,\mathrm m$、接線速度 $3\,\mathrm{m\,s^{-1}}$ の円運動では、速度が $3E_2$ なので加速度は $9\nabla_{E_2}E_2=-4.5 E_1\,\mathrm{m\,s^{-2}}$ である。角速度を使う計算では $v_\phi/r=1.5\,\mathrm{s^{-1}}$ を座標成分に取り、同じ値が得られる。単位基底の接線速度と座標成分の角速度は単位が異なるため、最初にどちらを使っているか確かめる。

## Levi-Civita 記号と ε–δ 恒等式

### 交代記号の縮約

右手系で外積の成分を計算すると符号の並びが繰り返し現れる。これを符号表として定義しておくと、二重外積を行列式の恒等式から処理できる。

:::definition[交代記号]{#def-mp-alternating-symbol}
$[i_1\cdots i_n]$ を、$(1,\ldots,n)$ の偶置換で $1$、奇置換で $-1$、添字に重複があるとき $0$ とする。この :term[レヴィチヴィタ記号]{en="Levi-Civita symbol" yomi="れゔぃちゔぃたきごう"} を用いて行列式を展開できる。
:::

:::proposition[交代記号の積]{#prop-mp-epsilon-delta}
3 次元の正向き正規直交基底では $\sum_k[ijk][\ell mk]=\delta_{i\ell}\delta_{jm}-\delta_{im}\delta_{j\ell}$ が成り立つ。
:::

:::proof
$i=j$ または $\ell=m$ なら両辺は零である。その他の場合、左辺が非零となるには二組の添字が同じ集合である必要がある。$(i,j)=(\ell,m)$ なら並びの積は $1$、$(i,j)=(m,\ell)$ なら符号は $-1$ であり、右辺も同じ値を取る。より一般に二つの交代記号の積は $\det(\delta_{i_a j_b})_{a,b=1}^n$ に等しい。重複があれば行または列が一致し、重複がなければ置換行列の行列式になるためである。
:::

:::example[三重ベクトル積]
右手系で $(\vb a\times\vb b)_i=[ijk]a_jb_k$ と書く。$\sum_k[ijk][k\ell m]=\delta_{i\ell}\delta_{jm}-\delta_{im}\delta_{j\ell}$ を使うと $\vb a\times(\vb b\times\vb c)=\vb b(\vb a\cdot\vb c)-\vb c(\vb a\cdot\vb b)$ となる。$\vb a=(1,1,0)$、$\vb b=(0,1,0)$、$\vb c=(0,0,2)$ の場合は $(0,0,-2)$ となる。
:::

:::example[回転の回転]
同じ縮約を $\partial_j\partial_\ell A_m$ に使うと、$\rot\rot\vb A=\grad\divg\vb A-\Delta\vb A$ となる。ここで右辺の $\Delta\vb A$ は直交座標で成分ごとに計算する。たとえば $\vb A=(0,x^2,0)$ では両辺が $(0,-2,0)$ になる。
:::

### 行列式と交代性

交代記号の自由添字を両側に三つずつ残すと、積は Kronecker のデルタでできた $3\times3$ 行列式である。左の添字列に重複があれば行が重なり、右の添字列に重複があれば列が重なる。どちらも重複がないとき、その行列式は置換行列の符号に等しい。したがってこの式を一つ縮約すると二項の恒等式が得られ、さらに縮約すると $\sum_{j,k}\epsilon_{ijk}\epsilon_{\ell jk}=2\delta_{i\ell}$ を得る。$i=\ell$ にして和を取れば値は $6$ となる。

直交座標で $A,B$ が滑らかなとき、$(A\times B)_i=\epsilon_{ijk}A_jB_k$ を微分すると
$$
\partial_i(A\times B)_i
 =B_k\epsilon_{kij}\partial_i A_j-A_j\epsilon_{jik}\partial_i B_k
 = B\cdot(\rot A)-A\cdot(\rot B)
$$
となる。第 1 項では三添字の巡回置換、第 2 項では一回の交換を使った。曲線座標で同じ推論を行うときは、体積形式から作った交代テンソルと共変微分を使う。

:::example[斜めの平行四辺形]
二方向が直交しない場合に、外積のノルムから面積を求める。

成分で内積と外積を計算する。外積の第一成分は第二・第三成分で作った二次の小行列式である。

$$
\begin{aligned}A&=\begin{pmatrix}1\\2\\0\end{pmatrix}\\ B&=\begin{pmatrix}1\\0\\2\end{pmatrix}\\ A\cdot B&=1\\ A\times B&=\begin{pmatrix}4\\-2\\-2\end{pmatrix}\end{aligned}
$$

二つの外積を内積で結ぶ体積の恒等式を使う。右辺は $2\times2$ の Gram 行列式である。

$$
\begin{aligned}\norm{A\times B}^2&=24\\ (A\cdot A)(B\cdot B)-(A\cdot B)^2&=24\end{aligned}
$$

第三のベクトル $C$ を入れる。体積は列が $A,B,C$ の行列式で得られる。

$$
\begin{aligned} C&=\begin{pmatrix}1\\1\\1\end{pmatrix}\\ (A\times B)\cdot C&=0\\ \det[A\ B\ C]&=0\end{aligned}
$$

体積の符号には並べる順番が必要である。順番を一回交換すると符号が逆になる。
:::

:::example[体積が零の三本]
三本のうち一つが他の二本の和である場合を調べる。

成分で内積と外積を計算する。外積の第一成分は第二・第三成分で作った二次の小行列式である。

$$
\begin{aligned}A&=\begin{pmatrix}2\\1\\0\end{pmatrix}\\ B&=\begin{pmatrix}0\\1\\1\end{pmatrix}\\ A\cdot B&=1\\ A\times B&=\begin{pmatrix}1\\-2\\2\end{pmatrix}\end{aligned}
$$

二つの外積を内積で結ぶ体積の恒等式を使う。右辺は $2\times2$ の Gram 行列式である。

$$
\begin{aligned}\norm{A\times B}^2&=9\\ (A\cdot A)(B\cdot B)-(A\cdot B)^2&=9\end{aligned}
$$

第三のベクトル $C$ を入れる。体積は列が $A,B,C$ の行列式で得られる。

$$
\begin{aligned} C&=\begin{pmatrix}2\\2\\1\end{pmatrix}\\ (A\times B)\cdot C&=0\\ \det[A\ B\ C]&=0\end{aligned}
$$

体積の符号には並べる順番が必要である。順番を一回交換すると符号が逆になる。
:::

:::example[右手系の符号]
向きを入れ替えた順序で有向体積を測る。

成分で内積と外積を計算する。外積の第一成分は第二・第三成分で作った二次の小行列式である。

$$
\begin{aligned}A&=\begin{pmatrix}0\\1\\0\end{pmatrix}\\ B&=\begin{pmatrix}1\\0\\0\end{pmatrix}\\ A\cdot B&=0\\ A\times B&=\begin{pmatrix}0\\0\\-1\end{pmatrix}\end{aligned}
$$

二つの外積を内積で結ぶ体積の恒等式を使う。右辺は $2\times2$ の Gram 行列式である。

$$
\begin{aligned}\norm{A\times B}^2&=1\\ (A\cdot A)(B\cdot B)-(A\cdot B)^2&=1\end{aligned}
$$

第三のベクトル $C$ を入れる。体積は列が $A,B,C$ の行列式で得られる。

$$
\begin{aligned} C&=\begin{pmatrix}0\\0\\2\end{pmatrix}\\ (A\times B)\cdot C&=-2\\ \det[A\ B\ C]&=-2\end{aligned}
$$

体積の符号には並べる順番が必要である。順番を一回交換すると符号が逆になる。
:::

## 擬テンソルとテンソル密度

### 向きと面積の変換を分ける

積分の変数を交換しても体積は正である。しかし循環の向きまで交換すると符号が変わる。この違いは、Jacobian を絶対値で取るかどうかに現れる。

:::definition[密度と向き]{#def-mp-density}
式 :eqref[eq-mp-tensor-transform] の右辺にさらに $|\det B|^w$ が掛かる成分を、重み $w$ の :term[テンソル密度]{en="tensor density" yomi="てんそるみつど"} という。さらに $\sgn(\det B)$ を掛ける規約のものを擬テンソル密度と呼ぶ。重み $0$ の擬テンソル密度は :term[擬テンソル]{en="pseudotensor" yomi="ぎてんそる"} である。
:::

$g'=B^\top g B$ という行列式の関係から $\sqrt{g'}=|\det B|\sqrt g$ であり、$\sqrt g$ はスカラー密度の成分である。空間に固定した向きをもつ体積形式の成分は $\varepsilon_{ijk}=\sigma\sqrt g[ijk]$ である。ここで $\sigma=1$ は正向きの座標、$\sigma=-1$ は逆向きの座標である。実際、三個の $B$ で $[ijk]$ を変換すると $\det B$ が掛かり、$\sigma'=\sigma\sgn(\det B)$ と一致する。固定の符号表 $[ijk]$ だけをどの座標でも使うと、テンソルとして変換した値とは一致しない。

この規約では、下付きの固定符号表 $[ijk]$ は重み $-1$ の擬テンソル密度とみなせる。通常の下付き三添字の変換を施すと $\det B$ が掛かるが、擬テンソル密度の追加因子 $\sgn(\det B)|\det B|^{-1}$ を掛けると二つの因子の積は $1$ になる。したがって新しい座標でも同じ符号表が得られる。一方、$\varepsilon_{ijk}$ は計量と固定した向きから決まる通常のテンソルである。

:::example[反射する軸性ベクトル]
直交行列 $O$ に対して $(O\vb a)\times(O\vb b)=(\det O)O(\vb a\times\vb b)$ である。両辺と任意の $O\vb c$ の内積を取れば、三重積の行列式の変換から証明できる。$O=\diag(-1,1,1)$ とし、$\vb a=(1,0,0)$、$\vb b=(0,2,0)$ とすると外積は $(0,0,2)$ から $(0,0,-2)$ に変わる。角運動量 $\vb r\times\vb p$ はこの軸性ベクトルの例である。
:::

:::example[局所的な座標の伸縮]
$x=2u$、$y=3v$、$z=-w$ のとき $\det B=-6$ なので $\dd V$ の密度は $6$、$\dd x\wedge\dd y\wedge\dd z$ の成分は $-6$ となる。単位座標立方体の体積は $6$ だが、$(u,v,w)$ の順で指定した向きは元の向きと逆である。
:::

### 積分で変換則を確かめる

旧座標を $q$、新座標を $q'$ とし $J=\partial q'/\partial q$ とすると、$\dd^n q=|\det J|^{-1}\dd^n q'$ となる。座標体積当たりの質量 $\rho$ を積分する場合には、質量が同じになる条件から $\rho'(q')=|\det J|^{-1}\rho(q)$ を得る。非線形な座標変換では Jacobian が位置によって変わるので、積分では被積分関数の一部として点ごとに評価する。計量から得られる体積因子も $g'=J^{-\mathsf T}gJ^{-1}$ と $\sqrt{\det g'}=|\det J|^{-1}\sqrt{\det g}$ に従う。したがって物理的体積当たりの質量 $\rho/\sqrt{\det g}$ はスカラーである。

有向体積で比較する場合は絶対値を外す。古い向きに固定した体積形式を新座標で書くと、その係数は $\det J$ が負の場合に反転する。これと、空間自体を反射する操作の結果を区別するとよい。後者では、例えば二本の極性ベクトルを反射した後の外積が、反射行列をそのまま外積に掛けた結果と逆向きになる。

:::example[一方向の反射]
第 1 軸だけを反射して向きを調べる。

旧座標から新座標へ行列 $J$ を掛ける。変換の向きを決める符号と、体積の倍率を分けて読む。

$$
\begin{aligned}J&=\begin{pmatrix}-1 & 0 & 0\\0 & 1 & 0\\0 & 0 & 1\end{pmatrix}\\ \det J&=-1\\ |\det J|&=1\end{aligned}
$$

旧座標の密度を $6$ とする。新座標での密度と、旧座標の単位箱の像に含まれる積分値を求める。

$$
\begin{aligned}\rho^\prime&=\frac{6}{|\det J|}=6\\ \int_{\mathrm{image}}\rho^\prime\,\dd^3 q^\prime&=6\end{aligned}
$$

旧座標の計量を単位行列とする。新座標の計量で単位箱の像の物理的な体積を計算する。

$$
\begin{aligned}g^\prime&=\begin{pmatrix}1 & 0 & 0\\0 & 1 & 0\\0 & 0 & 1\end{pmatrix}\\ \sqrt{g^\prime}|\det J|&=1\end{aligned}
$$

有向体積を使うなら係数は $1/\det J$ になる。密度の積分には体積の倍率の絶対値を使う。
:::

:::example[二方向の倍率]
二つの軸の長さの単位を変える。

旧座標から新座標へ行列 $J$ を掛ける。変換の向きを決める符号と、体積の倍率を分けて読む。

$$
\begin{aligned}J&=\begin{pmatrix}2 & 0 & 0\\0 & 3 & 0\\0 & 0 & 1\end{pmatrix}\\ \det J&=6\\ |\det J|&=6\end{aligned}
$$

旧座標の密度を $6$ とする。新座標での密度と、旧座標の単位箱の像に含まれる積分値を求める。

$$
\begin{aligned}\rho^\prime&=\frac{6}{|\det J|}=1\\ \int_{\mathrm{image}}\rho^\prime\,\dd^3 q^\prime&=6\end{aligned}
$$

旧座標の計量を単位行列とする。新座標の計量で単位箱の像の物理的な体積を計算する。

$$
\begin{aligned}g^\prime&=\begin{pmatrix}\frac{1}{4} & 0 & 0\\0 & \frac{1}{9} & 0\\0 & 0 & 1\end{pmatrix}\\ \sqrt{g^\prime}|\det J|&=1\end{aligned}
$$

有向体積を使うなら係数は $1/\det J$ になる。密度の積分には体積の倍率の絶対値を使う。
:::

:::example[反射と拡大の合成]
符号と体積倍率がともに変わる場合を計算する。

旧座標から新座標へ行列 $J$ を掛ける。変換の向きを決める符号と、体積の倍率を分けて読む。

$$
\begin{aligned}J&=\begin{pmatrix}0 & 2 & 0\\1 & 0 & 0\\0 & 0 & 3\end{pmatrix}\\ \det J&=-6\\ |\det J|&=6\end{aligned}
$$

旧座標の密度を $6$ とする。新座標での密度と、旧座標の単位箱の像に含まれる積分値を求める。

$$
\begin{aligned}\rho^\prime&=\frac{6}{|\det J|}=1\\ \int_{\mathrm{image}}\rho^\prime\,\dd^3 q^\prime&=6\end{aligned}
$$

旧座標の計量を単位行列とする。新座標の計量で単位箱の像の物理的な体積を計算する。

$$
\begin{aligned}g^\prime&=\begin{pmatrix}\frac{1}{4} & 0 & 0\\0 & 1 & 0\\0 & 0 & \frac{1}{9}\end{pmatrix}\\ \sqrt{g^\prime}|\det J|&=1\end{aligned}
$$

有向体積を使うなら係数は $1/\det J$ になる。密度の積分には体積の倍率の絶対値を使う。
:::

## 演習問題

::::exercise{level=basic #exr-mp-cylindrical-volume}
$1<\rho<2$、$0<z<3$ の円筒殻の体積を求めよ。座標の長さは $\mathrm m$ とする。
:::solution
体積要素は $\rho\,\dd\rho\dd\phi\dd z$ である。積分は $(4-1)/2$、$2\pi$、$3$ の積なので $9\pi\,\mathrm{m^3}$ となる。
:::
::::

::::exercise{level=basic}
円筒座標で $\vb A=\rho^3\hat{\vb e}_\rho$ の発散を求めよ。
:::solution
式 :eqref[eq-mp-orthogonal-operators] で $J=\rho$ と置く。角度と高さの成分がないので、$\divg\vb A=\rho^{-1}\partial_\rho(\rho^4)=4\rho^2$ である。
:::
::::

::::exercise{level=basic}
三次元の球座標で $r>0$ とする。$\Delta r^\alpha$ を求め、$\alpha=-1$ の場合を確認せよ。
:::solution
角度の微分は零である。$\partial_r(r^2\partial_r r^\alpha)=\alpha(\alpha+1)r^\alpha$ を $r^2$ で割ると $\Delta r^\alpha=\alpha(\alpha+1)r^{\alpha-2}$ であり、$\alpha=-1$ で零になる。計算の領域は $r>0$ なので、原点での分布の寄与は判定していない。
:::
::::

::::exercise{level=basic}
$x=u+v$、$y=2v$ における $g_{ij}$ と $\sqrt g$ を求めよ。
:::solution
列 $(1,0)$ と $(1,2)$ の内積から $g_{ij}$ は $\begin{pmatrix}1&1\\1&5\end{pmatrix}$ となる。行列式は $5-1=4$ で、$\sqrt g=2$ である。これは接ベクトルの行列式の絶対値とも一致する。
:::
::::

::::exercise{level=advanced}
$\partial_i\sqrt g=\sqrt g\,\Gamma^k_{ik}$ を示せ。
:::solution
逆行列を掛けて $g^{k\ell}\partial_i g_{k\ell}$ を取ると、計量の積の微分式の二項は対称性によりそれぞれ $\Gamma^k_{ik}$ になる。行列式を一成分ずつ微分すると係数はその余因子であるから $\partial_i g=g\,g^{k\ell}\partial_i g_{k\ell}$ である。両辺を $2\sqrt g$ で割れば結論を得る。
:::
::::

::::exercise{level=advanced}
球座標で $\Gamma^r_{\phi\phi}$、$\Gamma^\theta_{\phi\phi}$ を求めよ。
:::solution
$g_{\phi\phi}=r^2\sin^2\theta$ と $g_{rr}=1$、$g_{\theta\theta}=r^2$ を用いる。式 :eqref[eq-mp-christoffel] では対角成分の微分だけが残り、$\Gamma^r_{\phi\phi}=-r\sin^2\theta$、$\Gamma^\theta_{\phi\phi}=-\sin\theta\cos\theta$ となる。半径と角度の両方が緯線に沿う基底の変化に寄与している。
:::
::::

::::exercise{level=advanced}
正向きの一般座標で $\varepsilon_{ijk}\varepsilon^{\ell mk}=\delta_i^\ell\delta_j^m-\delta_i^m\delta_j^\ell$ を示せ。
:::solution
両辺は上二つ下二つの添字を持つテンソルである。左辺は計量で上げた体積形式と元の体積形式の積を縮約したものである。右辺は恒等写像の成分の積である。点ごとに正向き正規直交基底を選べば左辺は :ref[prop-mp-epsilon-delta] の左辺になり、右辺と等しい。従って任意の基底でも等しい。
:::
::::

::::exercise{level=advanced}
半径 $a$ の円筒側面を通る流束を、$\vb A=(x,y,0)$ について高さ $L$ の範囲で求めよ。発散定理でも確かめよ。
:::solution
$\vb A=a\hat{\vb e}_\rho$、面積は $2\pi aL$ なので流束は $2\pi a^2L$ である。上下の面には法線成分がない。全表面を使って直交座標で計算すると発散は $2$ であるため、体積 $\pi a^2L$ を掛けて同じ答えを得る。
:::
::::

## 確認テスト

::::quiz{answer="3" tol="1e-9"}
球座標で $r=3$、$\theta=\pi/2$ のとき、$h_\phi$ はいくらか。

:::explanation
$\partial_\phi\vb r$ の長さは $r\sin\theta$ であり、$\sin(\pi/2)=1$ なので $3$ である。
:::
::::

::::quiz
一般の座標で $\divg\vb A$ を計算するとき、$A^i$ が座標成分ならどの式を使うか。

- [ ] $\partial_i A^i$
- [x] $g^{-1/2}\partial_i(\sqrt g\, A^i)$
- [ ] $\partial_i(\sqrt g\, A^i)$

:::explanation
座標体積 $\dd^3q$ に対する流束差を $\partial_i(\sqrt g A^i)$ で集め、物理体積の倍率 $\sqrt g$ で割る。$\partial_i A^i$ は $\sqrt g$ が一定のときに使える。
:::
::::

::::quiz{answer="4" tol="1e-9"}
円筒座標で $\Delta(\rho^2)$ はいくらか。

:::explanation
$\partial_\rho(\rho\,\partial_\rho\rho^2)=\partial_\rho(2\rho^2)=4\rho$ となるため、$\rho$ で割ると $4$ である。
:::
::::

::::quiz
座標変換後もテンソルの通常の成分変換則に従うものはどれか。

- [x] $g_{ij}$
- [ ] $\Gamma^i_{jk}$
- [ ] どの座標でも同じ数値とした $[ijk]$

:::explanation
計量は二つの基底ベクトルの内積で各入力に線形なので、二つの逆 Jacobian が掛かる。Christoffel 記号には座標変換の二階微分が加わる。固定した交代記号は行列式の倍率を補う変換を持つ。
:::
::::

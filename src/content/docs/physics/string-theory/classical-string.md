---
title: 相対論的古典弦
description: 世界面の面積から古典弦の運動を導き、境界条件、保存量、回転解と補助計量の関係を調べる。
notation:
  - symbol: \sigma^a=(\tau,\sigma)
    meaning: 世界面の座標
  - symbol: X^\mu
    meaning: 標的時空への写像
  - symbol: h_{ab}
    meaning: 誘導計量
  - symbol: \gamma_{ab}
    meaning: 独立な世界面計量
---

## この章の目標

- 誘導計量から面積の作用を作る。
- 開弦の境界項を計算し、端点の条件を選ぶ。
- 静的な弦の波と回転する弦の保存量を求める。
- 補助計量の変分から拘束を求める。

## Nambu–Goto 作用

点粒子の世界線を弦に拡張すると、弦の位置を指定する $\sigma$ が加わる。時空は平坦で計量は $(-,+,\ldots,+)$ とする。世界面の座標は、標的時空の時刻や距離をそのまま用いる例を除き、無次元とする。点粒子の補助場は :ref[prop-string-particle-einbein] で導入した。

:::definition[世界面の面積作用]{#def-string-nambu-goto}
写像 $X^\mu(\tau,\sigma)$ が描く面を :term[世界面]{en="worldsheet" yomi="せかいめん"} といい、$h_{ab}=\partial_aX\cdot\partial_bX$ を誘導計量とする。$\det h<0$ の領域では
$$
S_{\rm NG}=-T\int\dd^2\sigma\,\sqrt{-\det h},
\qquad T=(2\pi\alpha')^{-1}
\label{eq-string-ng}
$$
を Nambu–Goto 作用という。
:::

計量の座標変換 $h\mapsto J^{\mathsf T}hJ$ で平方根に $|\det J|$ が現れ、測度を相殺する。標的時空の Poincaré 変換は内積を保つ。面積密度は $\sqrt{(\dot X\cdot X')^2-\dot X^2X'^2}$ である。例えば $X=(\tau,\sigma,0)$ を長さ・時間の座標で書くと密度は 1 である。同じ直線が横方向に速さ $v$ で動けば密度は $\sqrt{1-v^2}$ となる。したがって小片の時間成分には点粒子と同じ因子が掛かる。

### 面積密度と時間方向の選択

$\sigma$ が一定の曲線に沿う速度 $\dot X$ を接線 $X'$ に平行な部分と直交部分に分ける。$X'^2>0$ の領域で $\dot X_\perp=\dot X-X'(\dot X\cdot X')/X'^2$ とすれば $-\det h=-X'^2\dot X_\perp^2$ である。従って世界面が時間的という条件は $\dot X_\perp^2<0$ と同値である。標的時空の時間 $t$ を使って $X=(t,\vb x(t,\sigma))$ と書けば、面積密度は

$$
\sqrt{-h}=|\vb x'|\sqrt{1-|\dot{\vb x}_\perp|^2}.
$$

この関係は $\dot{\vb x}\mapsto \dot{\vb x}+u\vb x'$ で面積が変わらない理由を示す。固定した $\sigma$ に付けた点を弦の上で滑らせても、同じ世界面を描くためである。この分解で Lagrange 密度を $\dot{\vb x}$ について微分すると、運動量密度は $T|\vb x'|\dot{\vb x}_\perp/\sqrt{1-|\dot{\vb x}_\perp|^2}$ である。$\dot{\vb x}$ との内積を取り Lagrange 密度を引けば、エネルギー密度は $T|\vb x'|/\sqrt{1-|\dot{\vb x}_\perp|^2}$ となる。

:::example[直線に沿う座標速度]
$X=(\tau,2\sigma+3\tau,\tau/2)$ を長さの単位を選んで書いたとする。$h_{\tau\tau}=-1+9+1/4=33/4$、$h_{\tau\sigma}=6$、$h_{\sigma\sigma}=4$ より $-\det h=3$ である。$\tau$ の座標線自体は空間的であるが、直線に平行な速度を引くと横速度は $1/2$ なので世界面は時間的である。斜交座標で $h_{\tau\tau}<0$ を必須にすると、この世界面を誤って排除する。
:::

:::example[平行移動する円弧]
半径 $a$ の円弧を $X=(t,a\cos\sigma,a\sin\sigma,vt)$、$0<\sigma<\theta_0$ とする。$h_{tt}=-(1-v^2)$、$h_{t\sigma}=0$、$h_{\sigma\sigma}=a^2$ である。時間間隔 $\Delta t$ で作用は $-Ta\theta_0\Delta t\sqrt{1-v^2}$ になる。円弧を固定する力が必要なので、この写像は運動方程式を満たすと仮定していない。例えば $v=3/5$ では世界面の面積は、同じ標的時空の時間の静止円弧の $4/5$ 倍である。
:::

### 誘導計量の階数と直接の変分

弦の長さを座標から読み取るには、まず写像の微分が独立な二方向を作るかを調べる。例えば $X=(t,x,0)$ を $t=u+v$、$x=u-v$ で書くと、$h_{uu}=h_{vv}=0$ となる。それでも $h_{uv}=-2$ なので、世界面には面積がある。一方 $X=(\tau+\sigma,\tau+\sigma,0)$ では二本の接ベクトルは同じ光的ベクトルであり、$h$ は零行列になる。後者では $h^{ab}$ を用いた Euler–Lagrange 方程式を立てられない。

時空を符号 $(-,+,\ldots,+)$ の曲がった Lorentz 多様体へ広げ、局所座標で $G_{\mu\nu}(X)$ と書く。このとき
$$
h_{ab}=G_{\mu\nu}(X)\partial_aX^\mu\partial_bX^\nu
$$
である。$\delta X^\rho=Y^\rho$ に対する計量の変化を、その点の標的時空における共変微分 $D_aY^\mu=\partial_aY^\mu+\Gamma^\mu{}_{\rho\nu}\partial_aX^\rho Y^\nu$ で書く。$\nabla G=0$ から
$$
\delta h_{ab}=G(D_aY,\partial_bX)+G(\partial_aX,D_bY)
$$
を得る。実際、右辺の二つの接続項を足せば $\partial_\rho G_{\mu\nu}Y^\rho\partial_aX^\mu\partial_bX^\nu$ に一致する。行列式の微分 $\delta\det h=\det h\,h^{ab}\delta h_{ab}$ を使い、$G$ に関する部分積分を行うと、内部の変分は
$$
\delta S_{\mathrm{int}}
=T\int\dd^2\sigma\,\sqrt{-h}\,G\!\left(
\frac1{\sqrt{-h}}D_a(\sqrt{-h}\,h^{ab}\partial_bX),Y\right)
$$
となる。ここで $D_a$ は密度の標的時空の添字だけに接続を加え、世界面の発散には普通の座標微分を用いる。

この式の括弧の中は、世界面の平均曲率ベクトルの 2 倍である。これを見るために $D_a\partial_bX$ を接線方向と法線方向に分ける。$\partial_c G(\partial_aX,\partial_bX)$ を三通り並べ、$D_a\partial_bX=D_b\partial_aX$ を用いて整理すると、
$$
G(D_a\partial_bX,\partial_cX)
=\frac12(\partial_a h_{bc}+\partial_b h_{ac}-\partial_c h_{ab})
$$
となる。従って接線方向は誘導計量の Christoffel 記号 $\Gamma(h)^c{}_{ab}\partial_cX$ に等しい。$B_{ab}=D_a\partial_bX-\Gamma(h)^c{}_{ab}\partial_cX$ は法線方向にあり、発散式の括弧は $h^{ab}B_{ab}$ となる。ここで二次元の跡を $2$ で割った $H=h^{ab}B_{ab}/2$ を :term[平均曲率ベクトル]{en="mean curvature vector" yomi="へいきんきょくりつべくとる"} という。従って時間的で非退化な世界面の内部の方程式は $H=0$ である。

接線方向の変分 $Y=\xi^a\partial_aX$ は座標の取り替えである。上の結果から内部の変分の係数は任意の写像でこの方向と直交しており、接線方向に独立な運動方程式は生じない。自由端や固定端を動かす場合は残った境界項も調べる。面積を停留させる条件であるため、Lorentz 世界面について作用が常に極小になるとは言えない。

:::example[静止した弦は空間の測地線を描く]
時空計量が $-\dd t^2+g_{ij}(x)\dd x^i\dd x^j$ で、弦が $X=(t,x(\sigma))$ と静止しているとする。誘導計量は $\diag(-1,|\partial_\sigma x|_g^2)$ である。時間間隔 $\Delta t$ の作用は $-T\Delta t$ と空間内の弧長の積になる。弧長を $\sigma$ と選ぶと $h_{ab}=\diag(-1,1)$ で、時間方向の $B_{tt}$ はゼロ、$B_{\sigma\sigma}=\nabla_\sigma x'$ である。従って固定された二端の間の静的な内部解は空間の測地線になる。

例えば空間を半径 $R$ の二次元球面とすれば、その赤道を弧長 $s$ で $x(s)=R(\cos(s/R),\sin(s/R),0)$ と書ける。埋め込んだ三次元 Euclid 空間の加速度 $x''=-x/R^2$ は球面に垂直であるため、球面内の共変加速度はゼロである。これに対し余緯度 $\vartheta$ を固定した緯線では、単位接ベクトル $\vb e_\phi$ の共変微分の $\vb e_\vartheta$ 成分は $-\cot\vartheta/R$ となる。実際、三次元での微分は水平方向の単位動径に $-1/(R\sin\vartheta)$ を掛けたもので、これと $\vb e_\vartheta$ の内積は $-\cos\vartheta/(R\sin\vartheta)$ である。赤道以外の緯線を保つには弦に外力が必要である。
:::

:::example[平面を小さく曲げたときの平均曲率]
Minkowski 空間のグラフ $X=(t,x,\delta f(t,x))$ を考え、$\delta$ を小さく取る。誘導計量は $\diag(-1,1)+O(\delta^2)$、法線の横成分は $1+O(\delta^2)$ である。従って平均曲率の横成分は $\delta(-f_{tt}+f_{xx})/2+O(\delta^3)$ となる。例えば $f=t^2+x^2$ では一次まで消えるが、$f=t^2$ では $-\delta$ が残る。個々の座標曲線に曲率があっても、時間と空間の寄与が反対符号になるため、世界面の曲率の跡は消えうる。
:::

## 運動方程式と境界条件（Neumann 条件・Dirichlet 条件）

### 変分と運動量流束

$\delta\sqrt{-h}=\sqrt{-h}h^{ab}\partial_a X\cdot\partial_b\delta X$ である。$P^a_\mu=-T\sqrt{-h}h^{ab}\partial_bX_\mu$ と置いて部分積分すると
$$
\delta S=\int\dd\tau [P^\sigma_\mu\delta X^\mu]_{0}^{\pi}
-\int\dd^2\sigma\,\partial_aP^a_\mu\,\delta X^\mu .
\label{eq-string-boundary-variation}
$$
始時と終時で変分を固定した。内部の任意性は $\partial_aP^a=0$ を与える。閉弦では周期条件で端の差が消える。開弦では自由な座標の流束をゼロとする :term[Neumann 条件]{en="Neumann boundary condition" yomi="のいまんじょうけん"}、または変分そのものをゼロとする :term[Dirichlet 条件]{en="Dirichlet boundary condition" yomi="でぃりくれじょうけん"} を方向ごとに課せる。後者では固定物体に運動量が渡りうる。

共形座標では $\sqrt{-h}h^{ab}=\diag(-1,1)$ となるから内部で $\ddot X-X''=0$ である。例えば Neumann 方向の $\cos(n\sigma)\cos(n\tau)$、Dirichlet 方向の $\sin(n\sigma)\cos(n\tau)$ はどちらも整数 $n$ で波動方程式と各境界条件を満たす。固定位置が両端で異なる場合は、Dirichlet の解に端点の座標を結ぶ一次関数を加える。ここで例示した個々の成分から弦の解を作るには計量の拘束も満たす必要がある。

### 固定端で反射する波と混合端

共形座標で内部解を $X=F(\tau+\sigma)+G(\tau-\sigma)$ と書く。$\sigma=0$ の Neumann 条件は $F'(u)=G'(u)$、Dirichlet 条件は $F'(u)=-G'(u)$ を与える。Dirichlet 端の位置は積分定数 $F+G$ で指定する。$\sigma=\pi$ に Neumann 条件を選べば $F'(\tau+\pi)=G'(\tau-\pi)$ である。従って両端 Neumann なら $F'$ は $2\pi$ 周期であるが、$\sigma=0$ が Dirichlet なら $F'$ は $2\pi$ の移動で符号が変わる。後者の Fourier 指数は $e^{-ir u}$、$r\in\Z+1/2$ である。

:::example[片端を固定した基音]
$Y=A\sin(\sigma/2)\cos(\tau/2)$ は $\sigma=0$ でゼロ、$\sigma=\pi$ で $Y'=A\cos(\pi/2)\cos(\tau/2)/2=0$ である。$\ddot Y=Y''=-Y/4$ なので、振動数は $1/2$ である。混合端の振動を $\sin(n\sigma)$ の整数列で表そうとすると、自由端で勾配がゼロにならない。
:::

:::example[二端間の張力]
$X=(L\tau/\pi,L\sigma/\pi,0)$ の端点は距離 $L$ だけ離れている。誘導計量は $(L/\pi)^2\diag(-1,1)$ であるから、$P^\sigma_1=-TL/\pi$ となる。この流束を標的時刻当たりに直すと $\dd\tau/\dd t=\pi/L$ より $-T$ である。従って弦は両端を内向きに大きさ $T$ の力で引く。
:::

### 境界の質量を変分する

弦の端が外部の壁で固定される場合とは別に、質量 $m_i$ の粒子を端に付けることができる。端の作用は $-m_i\int\sqrt{-\dot X_i^2}\dd\tau$ である。時間の両端を固定するとその変分は $-\int\dot p_{i\mu}\delta X_i^\mu\dd\tau$、$p_i=m_i\dot X_i/\sqrt{-\dot X_i^2}$ となる。:eqref[eq-string-boundary-variation] と足すと $\dot p_\pi=P^\sigma(\pi)$、$\dot p_0=-P^\sigma(0)$ である。世界面の保存式を積分した $\dot P_{\rm string}=-P^\sigma(\pi)+P^\sigma(0)$ と足せば、弦と端の粒子全体の運動量が保存する。

:::example[低速度の端が与える動的境界条件]
静的直線の長さを $L$ とし、片端 $x=0$ の横変位を固定し、もう一端には縦位置 $x=L$ を保ったまま横に動ける質量 $m$ を付ける。横変位の二次まで粒子の作用は $\int m\dot Y(L)^2/2\dd t$ である。従って $m\ddot Y(L)=-T Y'(L)$ である。$Y\propto e^{-ikt}\sin(kx)$ を入れると $Tk\cos(kL)=mk^2\sin(kL)$、すなわち $\cot(kL)=mk/T$ となる。$m=0$ では $\cos(kL)=0$ に戻り、$m$ が非常に大きいと正の非零の高い根は $\sin(kL)=0$ に近づく。最低の根には端の慣性が効き、$kL\ll1$ の $\cot(kL)\sim1/(kL)$ から $k^2\sim T/(mL)$ となる。
:::

### 同じ方程式を境界ごとの Fourier 基底で解く

波動方程式の一成分 $Y$ に対し、境界条件が時間によらないとする。Dirichlet 端に指定された値の差は $\sigma$ の一次関数で取り出し、残った斉次問題を解く。$\sigma=0,\pi$ の両方が Dirichlet なら $\sin(n\sigma)$、両方が Neumann なら $1,\cos(n\sigma)$ が基底である。左 Dirichlet、右 Neumann なら $\sin((n+1/2)\sigma)$、左 Neumann、右 Dirichlet なら $\cos((n+1/2)\sigma)$ になる。前二者で非零の整数モードは $n\ge1$、後二者は $n\ge0$ で数える。

これらが単なる特解の一覧に留まらないことは反射で確かめられる。左端で偶反射する場合は導関数がゼロ、奇反射する場合は関数がゼロである。右端でも同じように反射すると、二回の反射は距離 $2\pi$ の平行移動になる。二端が同じ種類なら符号は正で $2\pi$ 周期、異なるなら符号は負で $4\pi$ 周期になる。後者を $4\pi$ 周期で Fourier 展開すれば偶数モードは消え、半整数波数だけが残る。左端の反射の偶奇で余弦か正弦かが決まる。

初期値 $Y_0,Y_1$ をこれらの基底に射影し、$L^2$ 内積を $\int_0^\pi f(\sigma)g(\sigma)\dd\sigma$ とする。正の固有値 $\nu^2$ に対応する基底関数 $\varphi_\nu$ のノルム二乗は $\pi/2$ なので、係数は
$$
q_\nu(\tau)=\frac2\pi\int_0^\pi\dd\sigma\,\varphi_\nu(\sigma)
\left[Y_0(\sigma)\cos(\nu\tau)+Y_1(\sigma)\frac{\sin(\nu\tau)}{\nu}\right].
$$
Neumann の定数モードでは係数は $\pi^{-1}\int_0^\pi(Y_0+\tau Y_1)\dd\sigma$ である。滑らかな周期延長を持つ初期値なら、部分積分を繰り返して Fourier 係数の減衰を示せるので、微分まで項別に行える。例えば一定値の Dirichlet 境界を持つ滑らかな解を求めるには、端で $Y_0$ と $Y_1$ の斉次値に加え、方程式で決まる高階微分についての両立条件も必要になる。初期値の境界との両立を確認せずに Fourier 級数の滑らかさを仮定することはできない。

:::example[離れた固定端の零モードと第一振動]
時間と一つの空間方向を $X^0=L\tau/\pi$、$X^1=L\sigma/\pi$ とする。$X^1$ の差 $L$ は Dirichlet 条件で固定され、$X^1$ の定数速度のモードはない。もう一方向で Dirichlet 値を等しくし、$Y_0=A\sin\sigma$、$Y_1=0$ とすれば、斉次波動方程式の成分解は $Y=A\sin\sigma\cos\tau$ である。$Y'^2+\dot Y^2$ は一般にはゼロでないので、この $Y$ を元の $X^0,X^1$ にそのまま足した写像は共形拘束を二次で破る。小振幅なら横の一次解として使い、有限振幅なら $X^0$ や他の成分の補正も求める。
:::

:::example[基音の初速度を持つ混合端の成分]
左固定、右自由で、$Y_0=0$、$Y_1=B\sin(\sigma/2)$ とすると、解は $Y=2B\sin(\tau/2)\sin(\sigma/2)$ である。空間に一様な初速度を選ぶと左端で $Y_1(0)\ne0$ になり、滑らかな固定端の解を作れない。実際、左端の $Y(\tau,0)=0$ を時間で微分すると $Y_1(0)=0$ が必要である。初速度を与えるときには端点が実際に動ける方向を先に指定する。
:::

### 境界の背後にある対称性

境界の値を点ごとに固定する代わりに、端が一定の平坦な部分空間 $W$ の上を動けるようにしてもよい。$\delta X$ は $W$ に接し、境界項が任意の接ベクトルに対して消えるためには $P^\sigma$ が $W$ に直交しなければならない。$X^\perp$ の境界値は固定され、$X^\parallel$ の境界値は任意である。このような境界を持つ弦では、$W$ に平行な並進の運動量だけが弦の内部で保存する。二つの端が別々の面にある場合は、両面に接する方向の並進が保存則を与える。

共形座標で平坦な固定面の境界を考える。Neumann 成分は $X'^\parallel=0$、Dirichlet 成分は $\dot X^\perp=0$ なので、$\dot X\cdot X'=0$ が境界で自動的に成立する。これは世界面の空間方向へエネルギーが流れない条件である。Dirichlet 成分の $X'^\perp$ は一般に非零なので、標的時空の横運動量は面に渡る。世界面の混合応力と標的の並進カレントは異なる量であるため、両条件は両立する。

両端が自由で座標を $[0,\pi]$ に取ると、共形ゲージで各成分は
$$
X^\mu=x^\mu+2\alpha'p^\mu\tau
+i\sqrt{2\alpha'}\sum_{n\ne0}\frac{\alpha_n^\mu}{n}e^{-in\tau}\cos(n\sigma)
\label{eq-string-classical-open-modes}
$$
と書ける。実条件は $\alpha_{-n}^\mu=\bar\alpha_n^\mu$ である。これは $F$ と $G$ の微分を同じ周期関数にした解で、定数の差は $x^\mu$ に吸収した。$p^\mu=T\int_0^\pi\dot X^\mu\dd\sigma$ であり、$T\pi(2\alpha')=1$ が零モードの係数を決める。Dirichlet 方向では端を結ぶ一次項に正弦級数を加え、自由並進の零モードを外す。この展開で決まるのは運動方程式と境界条件の解であり、物理的な世界面になるには、後の節で求める二つの二次拘束を満たす係数を選ぶ。

## 静的ゲージと弦の運動

### 張った弦の横波

$X^0=t$、$X^1=x$ と選べる局所領域を考え、横変位を $Y(t,x)$ とする。$-\det h=1+Y'^2-\dot Y^2-\dot Y^2Y'^2+(\dot Y\cdot Y')^2$ である。勾配が小さいと仮定すると
$$
S=-T\int\dd t\dd x+\frac T2\int\dd t\dd x\,(\dot Y^2-Y'^2)+O((\partial Y)^4).
$$
Euler–Lagrange 方程式は $\ddot Y-Y''=0$ である。横方向の速さは 1 であり、基準の直線に沿う変位は座標の変更に吸収された。固定端間距離 $L$ なら $Y=A\sin(\pi x/L)\cos(\pi t/L)$ が基本振動である。その二次エネルギーは $TA^2\pi^2/(4L)$ となる。これは二乗積分 $\int_0^L\sin^2(\pi x/L)\dd x=L/2$ から得られる。

単一の横成分で $Y=f(t-x)$ なら平方根の中は厳密に 1 になる。ただし作用をその値に置き換えて変分する手順は使えない。変分後の式 $\partial_t(\dot Y/\sqrt{-h})-\partial_x(Y'/\sqrt{-h})=0$ に代入すれば、任意の滑らかな $f$ が内部の解になる。

### 有限振幅での Hamilton 密度

横成分一つの静的ゲージでは $\mathcal L=-T\sqrt{1+Y'^2-\dot Y^2}$ である。$P_Y=T\dot Y/\sqrt{1+Y'^2-\dot Y^2}$ と置いて二乗し、速度の符号を $P_Y$ と合わせると $\dot Y=P_Y\sqrt{1+Y'^2}/\sqrt{T^2+P_Y^2}$ となる。Legendre 変換の結果は

$$
\mathcal H=\sqrt{(T^2+P_Y^2)(1+Y'^2)}.
$$

このハミルトニアンの $\dot Y=\partial\mathcal H/\partial P_Y$ と $\dot P_Y=\partial_x(\partial\mathcal H/\partial Y')$ は、元の変分方程式と一致する。小振幅の規準は $|\partial Y|\ll1$ であり、特に波長が短ければ振幅の長さだけを見て近似を決められない。

:::example[傾いたまま動く直線]
$Y=3x/4+t/2$ で $1+Y'^2-\dot Y^2=21/16$ である。$P_Y=2T/\sqrt{21}$、エネルギー密度は $25T/(4\sqrt{21})$ となる。$\dot Y$ の標的時空の直線に垂直な成分の速さは $(1/2)/\sqrt{1+(3/4)^2}=2/5$ であり、曲線の長さの係数 $5/4$ と組み合わせると同じ密度となる。
:::

:::example[孤立した波束のエネルギー]
$Y=f(t-x)$ をコンパクトな台を持つ滑らかな波束とする。$\dot Y=f'$、$Y'=-f'$ より $\mathcal H=T(1+(f')^2)$ である。直線分のエネルギーを引くと $\Delta E=T\int(f'(t-x))^2\dd x$ となり、$u=t-x$ の置換で時刻に依存しない。またエネルギー流束 $-TY'\dot Y/\sqrt{1+Y'^2-\dot Y^2}$ は $T(f')^2$ であるため、余剰エネルギーは右向きに単位速度で進む。
:::

### 静的ゲージの完全な運動方程式

静的ゲージで複数の横方向を一度に残す。$u=\partial_t\vb Y$、$v=\partial_x\vb Y$ と書き、$\mathcal D=(1-u^2)(1+v^2)+(u\cdot v)^2$ と置く。作用密度 $-T\sqrt{\mathcal D}$ の微分から、横運動量と流束は
$$
\begin{aligned}
P_Y &= \frac{T}{\sqrt{\mathcal D}}\bigl[(1+v^2)u-(u\cdot v)v\bigr],\\
F_Y &=-\frac{T}{\sqrt{\mathcal D}}\bigl[(1-u^2)v+(u\cdot v)u\bigr].
\end{aligned}
$$
したがって方程式は $\partial_t P_Y+\partial_x F_Y=0$ である。これは標的時空の横成分の変分そのものである。世界面の再パラメータ化不変性により、全成分の変分の係数 $E_\mu$ は恒等的に $E_\mu\partial_aX^\mu=0$ を満たす。$E_i=0$ が横成分で成立すれば、$\partial_t X^0=1$、$\partial_x X^1=1$ より残りの $E_0,E_1$ もゼロになる。グラフ表示が成り立つ領域では、横方程式だけで内部の全方程式を満たす。

時空の時間の並進から得られる密度は $\mathcal H=P_Y\cdot u+T\sqrt{\mathcal D}=T(1+v^2)/\sqrt{\mathcal D}$ である。$P_Y^2+(P_Y\cdot v)^2$ を上の式で計算すると、$T^2(1+v^2)[(1+v^2)u^2-(u\cdot v)^2]/\mathcal D$ になる。$\mathcal D=1+v^2-[(1+v^2)u^2-(u\cdot v)^2]$ を用いると、完全な Hamilton 密度は
$$
\mathcal H=\sqrt{T^2(1+v^2)+P_Y^2+(P_Y\cdot v)^2}
$$
である。$P_Y$ に関する微分は $(P_Y+(P_Y\cdot v)v)/\mathcal H$ となり、先の $P_Y$ の式を逆にした $u$ を与える。$P_Y,v$ が小さければ $\mathcal H=T+(P_Y^2+T^2v^2)/(2T)+\cdots$ となる。

運動量の方程式を使えば、エネルギーにも局所保存則がある。$\partial_t\mathcal H=(\partial_tP_Y)\cdot u-(\partial_v\mathcal L)\cdot\partial_t v$ であり、$\partial_t v=\partial_xu$ だから $\partial_t\mathcal H+\partial_x(F_Y\cdot u)=0$ となる。従って固定端 $u=0$ でも、力を受けない横自由端 $F_Y=0$ でもエネルギーは保存する。端を外部装置で動かすと $F_Y\cdot u$ が境界で非零になり、装置と弦の間で仕事が行われる。

:::example[一方向の波は横偏極がいくつあっても厳密解である]
$\vb Y=\vb f(t-x)$ とすると $v=-u$ なので $\mathcal D=1$、$P_Y=Tu$、$F_Y=Tu$ である。従って $\partial_tP_Y+\partial_xF_Y=0$ が成り立つ。また $\mathcal H=T(1+u^2)$、$F_Y\cdot u=Tu^2$ である。波の余剰エネルギーは右向きの流束と等しい。例えば $\vb f(s)=A(\cos ks,\sin ks)$ なら $u^2=A^2k^2$ が一定であるから、単位長さ当たり $TA^2k^2$ の余剰エネルギーが移動する。
:::

:::example[逆向きの波を加えると非線形項が現れる]
横成分を一つにして $Y=f(t-x)+g(t+x)$ と試す。$u=f'+g'$、$v=-f'+g'$ のため $\mathcal D=1-4f'g'$ である。$a=f'$、$b=g'$ と書けば
$$
\partial_t(u/\sqrt{\mathcal D})-\partial_x(v/\sqrt{\mathcal D})
=\frac{4(a^2b'+b^2a')}{(1-4ab)^{3/2}} .
$$
この式は $\partial_t-\partial_x=2\partial_{t-x}$、$\partial_t+\partial_x=2\partial_{t+x}$ を使って得る。例えば $f(s)=g(s)=\varepsilon\sin s$ なら、一般の点で右辺は $O(\varepsilon^3)$ で非零である。共形座標では標的の各成分は左右進行波に分かれたが、静的ゲージでは $t,x$ の両方を標的座標に合わせたため、一般の二波をそのまま重ねることはできない。
:::

### 特性方向とゲージの破れ

横変位一成分では、世界面の光的方向の速度 $\lambda=\dd x/\dd t$ は
$$
(-1+u^2)+2uv\lambda+(1+v^2)\lambda^2=0,\qquad
\lambda_\pm=\frac{-uv\pm\sqrt{1+v^2-u^2}}{1+v^2}
$$
である。運動方程式の最高階だけを書いても $(1+v^2)Y_{tt}-2uvY_{tx}+(-1+u^2)Y_{xx}=0$ となり、同じ速度が出る。例えば $u=0$ の傾いた直線では $\lambda_\pm=\pm(1+v^2)^{-1/2}$ である。標的時空で線に沿って進む速度は $\sqrt{1+v^2}|\lambda_\pm|=1$ なので、座標 $x$ の速度が小さいのは線の投影による。

$\mathcal D>0$ では二つの特性が異なり、$\mathcal D\to0$ で一致する。一方、世界面が滑らかでもその $(t,x)$ への投影の階数が落ちれば、曲面全体を同じグラフで書けなくなる。回転する弦が $x$ 軸に垂直になる時刻がその例である。その場合は標的時刻と別の空間方向に合わせるか、弦全体に沿う世界面座標を用いる。このように座標が使えない状況と誘導計量が退化する状況は別に判定する。

## 世界面のカレントと保存量

### 並進と回転

並進のカレントは :eqref[eq-string-boundary-variation] に現れた $P^a_\mu$ である。Lorentz 変換には
$$
J^{a\mu\nu}=X^\mu P^{a\nu}-X^\nu P^{a\mu}
$$
が対応する。発散をとると $\partial_aP^a$ は運動方程式で消え、残りは $h^{ab}(\partial_aX^\mu\partial_bX^\nu-\partial_aX^\nu\partial_bX^\mu)$ なので対称性から消える。$P^\mu=\int P^{\tau\mu}\dd\sigma$ に対し $\partial_\tau P^\mu=-[P^{\sigma\mu}]_0^\pi$ である。従って閉弦では全方向で保存し、開弦では端点が保つ対称性の成分が保存する。たとえば端点を $X^1=\text{一定}$ に固定しても $X^0$ を自由にすれば、エネルギー流束はゼロとなる。固定壁に垂直な運動量には同じ議論を使えない。

### 時間局所化した変換から保存則を求める

標的時空の定数ベクトル $a^\mu$ に対する並進を $\delta X^\mu=\varepsilon(\tau,\sigma)a^\mu$ と局所化する。定数 $\varepsilon$ では変化ゼロだが、局所化すると $\delta S=\int P^a\cdot a\,\partial_a\varepsilon\,\dd^2\sigma$ となる。内部に台を持つ $\varepsilon$ で部分積分すれば運動方程式上 $\partial_a(P^a\cdot a)=0$ である。$\varepsilon$ が時間のみに依存する場合、端点からの角運動量の流出入も

$$
\frac{\dd J^{\mu\nu}}{\dd\tau}
=-[X^\mu P^{\sigma\nu}-X^\nu P^{\sigma\mu}]_0^\pi
$$

と求まる。固定したブレインの面を保つ変換だけを考えれば、端の変位はその面に接し、Neumann 条件が境界項を消す。

:::example[力は残るがトルクは消える配置]
弦を $x^1$ 軸上の二点に固定した静的解を選ぶ。各端で $P^\sigma$ の空間成分は同じ軸に沿うため、空間の反対称積 $X^iP^{\sigma j}-X^jP^{\sigma i}$ はゼロである。しかし各固定点には有限な力が必要である。角運動量の流束がゼロであることと、各点で運動量を交換しないことは別の条件である。
:::

:::example[同じ軸に沿って並進する端]
両端が平行なブレインにあり、$x^2$ が両方で Neumann とする。状態が時間とともに形を変えても、$P^{\sigma 2}(0)=P^{\sigma 2}(\pi)=0$ なので全運動量の第 2 成分は一定である。時間も Neumann ならエネルギーについても同じ式が成り立つ。ブレインが時間に依存して傾く場合はこの境界対称性を改めて確認する。
:::

### 保存電荷を生成子として確かめる

等時面で $P_\mu$ を $X^\mu$ の正準運動量とし、
$$
p^\mu=\int\dd\sigma\,P^\mu,\qquad
J^{\mu\nu}=\int\dd\sigma\,(X^\mu P^\nu-X^\nu P^\mu)
$$
と定める。Poisson 括弧を $\{X^\rho(\sigma),P_\lambda(\sigma')\}=\delta^\rho{}_\lambda\delta(\sigma-\sigma')$ とする。デルタ関数を積分すれば
$$
\begin{aligned}
\{X^\rho,p^\mu\}&=\eta^{\rho\mu},\\
\{X^\rho,J^{\mu\nu}\}&=X^\mu\eta^{\rho\nu}-X^\nu\eta^{\rho\mu},\\
\{P^\rho,J^{\mu\nu}\}&=P^\mu\eta^{\rho\nu}-P^\nu\eta^{\rho\mu}
\end{aligned}
$$
を得る。第 1 式は一定の変位、第 2 式と第 3 式はベクトルの Lorentz 変換を表す。そのため積の括弧を使うと
$$
\begin{aligned}
\{p^\rho,p^\mu\}&=0,\qquad
\{p^\rho,J^{\mu\nu}\}=p^\mu\eta^{\rho\nu}-p^\nu\eta^{\rho\mu},\\
\{J^{\rho\lambda},J^{\mu\nu}\}
&=\eta^{\rho\nu}J^{\mu\lambda}-\eta^{\rho\mu}J^{\nu\lambda}
+\eta^{\lambda\nu}J^{\rho\mu}-\eta^{\lambda\mu}J^{\rho\nu}
\end{aligned}
$$
となる。これらは境界が変換を許す成分の正準代数である。この事実を使えば、量子化後に同じ代数が残るかを、対応する演算子の交換子で検査できる。

角運動量の値は標的時空の原点を選ぶ。$X\mapsto X+a$ なら $J^{\mu\nu}\mapsto J^{\mu\nu}+a^\mu p^\nu-a^\nu p^\mu$ である。空間の全運動量がゼロなら空間の角運動量は原点によらない。一方 $J^{0i}=X^0p^i-\int X^iP^0\dd\sigma$ は一般には原点を移すと変わる。標的時刻 $X^0=t$ の等時面でエネルギーの重心を $x_{\rm E}^i=E^{-1}\int X^iP^0\dd\sigma$ とすると、$J^{0i}=tp^i-E x_{\rm E}^i$ から $\dd x_{\rm E}^i/\dd t=p^i/E$ である。ここではエネルギーと $J^{0i}$ を保つ閉弦または全方向自由な開弦を考えている。

:::example[並進させた回転解]
静止系で回る弦を定数 $\vb a$ だけずらす。$p^i=0$ なので $J^{12}$ は変わらないが、$J^{0i}$ は $-Ea^i$ だけ変わる。全エネルギーの重心だけが $\vb a$ へ移ったためである。同じ弦を Lorentz ブーストすると $p^i$ は一般に非零になり、空間角運動量にも原点の選択が入る。従って軌跡の比較には静止系の $J$ を使う。
:::

### 時間に依存する拘束と境界の仕事

時間に依存する外部の端を $X=(t,\vb x(t,\sigma))$ で記述する。標的の並進対称性が作用の内部にあるため、端が動いても内部の保存式 $\partial_aP^{a0}=0$ は成立する。時間に沿って積分すると
$$
\frac{\dd E}{\dd t}=-P^{\sigma0}(t,\pi)+P^{\sigma0}(t,0)
$$
である。共役運動量の定義より $P^\sigma_\mu\dot X^\mu=-T\sqrt{-h}h^{\sigma a}h_{at}=0$ である。$\dot X^0=1$ を入れると $P^{\sigma0}=\vb P^\sigma\cdot\dot{\vb x}$ となり、エネルギー変化は端に作用する外部の力の仕事率に一致する。静止した Dirichlet 端ではこの仕事率はゼロだが、空間の運動量流束は残る。

:::example[引き離す壁に張力がする仕事]
$X=(t,\sigma L(t)/\pi,0)$、$L(t)>0$ とし、左端を固定して右端を $L(t)$ まで動かす。$|\dot L|<1$ と仮定すれば端の世界線は時間的である。世界面は $(t,x)$ 平面の領域であり、面積密度は $L/\pi$、弦のエネルギーは $TL$ である。作用を $L$ で変分すると弦が右の壁に及ぼす力は $-T$ である。壁が $\dot L>0$ で移動するとき、外部装置は $T\dot L$ の仕事率で弦にエネルギーを供給する。点の速度のうち接線方向の部分は面積密度を変えないが、ここでは境界そのものの時空領域が変わる。
:::

:::example[閉弦における周期条件と流束]
弦が $\sigma=0$ で向きを変えて見える座標表示でも、$X$ とその微分が周期的なら $P^\sigma(0)=P^\sigma(2\pi)$ である。積分区間の切れ目を別の場所へ動かしても、流出した運動量は同じ点の反対側から入る。コンパクトな標的方向に巻き付く場合には座標表示の $X$ 自体に整数周期の差がありうるが、微分は周期的である。その方向の並進運動量は保存し、コンパクト化で破れる回転については保存を別に検討する。
:::

## 回転する開弦と Regge 軌跡

### 回転解とその物理量

開弦で $X^0=A\tau$、$X^1+iX^2=A e^{i\tau}\cos\sigma$ とする。$0<\sigma<\pi$ では $h=A^2\sin^2\sigma\,\diag(-1,1)$ である。$\partial_\sigma X=0$ が両端で成り立ち、運動方程式も満たす。端点の速さは 1 であり、端で縮退する誘導計量は内部からの極限で扱う。

$$
E=\pi TA,\qquad J_{12}=TA^2\int_0^\pi\cos^2\sigma\,\dd\sigma=\frac{\pi TA^2}{2},
\qquad J_{12}=\alpha'E^2 .
\label{eq-string-classical-regge}
$$

空間運動量は $\int_0^\pi\cos\sigma\dd\sigma=0$ で消えるので、ここで $E=M$ である。$A=2\ell_s$ なら $E=1/\ell_s$、$J_{12}=1$ である。$A$ を 2 倍するとエネルギーは 2 倍、角運動量は 4 倍となる。古典解での切片はゼロであり、量子化では真空の振動が切片を変える。

### 折れた閉弦と傾きの比

同じ局所解を $0\le\sigma<2\pi$ で周期化すると、$\sigma=0,\pi$ が折り返し点であり、$\cos\sigma$ は開弦の線分を往復する。折り返しの点で誘導計量は縮退するが、Polyakov 作用と共形拘束では $X$ とその微分は滑らかである。$P^{\tau 0}=TA$、$J^{\tau12}=TA^2\cos^2\sigma$ を積分すると、閉弦で $E=2\pi TA$、$J=\pi TA^2$ である。従って $J=(\alpha'/2)E^2$ となる。エネルギーと角運動量が同時に 2 倍になるため、同じエネルギーに直した傾きは半分である。

:::example[同じ質量における比較]
$E=2/\ell_s$ とする。開弦では $A=4\ell_s$ で $J=4$、閉弦では $A=2\ell_s$ で $J=2$ となる。閉弦の折れ曲がりを一点の境界と見なして Neumann 境界項を加える必要はなく、世界面は円周のままである。
:::

:::example[端点に近い場所のエネルギー]
開弦の回転平面で中心からの符号付き距離を $r=A\cos\sigma$ とすれば、$\dd E=TA|\dd\sigma|=T\dd r/\sqrt{1-r^2/A^2}$ である。内側の $|r|<A/2$ のエネルギーは $2TA\arcsin(1/2)=\pi TA/3$ であり、全体の $1/3$ である。座標 $r$ の密度は端で発散するが、発散は平方根型なので積分できる。
:::

### 端点質量のある剛体回転

同じ質量 $m$ を両端に取り付け、物理時刻 $t$ と半径の符号付き座標 $r$ によって $X=(t,r\cos\omega t,r\sin\omega t)$ とする。$-a\le r\le a$ で $v=\omega a<1$ とする。内部の面積作用は $-T\int\dd t\int_{-a}^{a}\sqrt{1-\omega^2 r^2}\dd r$ であり、両端の作用は $-2m\int\sqrt{1-\omega^2a^2}\dd t$ である。内部の空間向きの変分方程式は、径方向の単位ベクトルを $e_r=(\cos\omega t,\sin\omega t)$ とすると $\partial_t P^t=-T\omega^2r e_r/\sqrt{1-\omega^2r^2}$、$\partial_r P^r=+T\omega^2r e_r/\sqrt{1-\omega^2r^2}$ から満たす。時間方向は各カレントが時間に依存せず、$P^{r0}=0$ である。固定した $\omega$ で $a$ を変分すると端の条件は $T(1-v^2)=m\omega v$ となる。実際、各端からの運動量流束の物理時刻当たりの大きさは $T\sqrt{1-v^2}$、円運動する粒子の運動量の時間微分は $m\omega v/\sqrt{1-v^2}$ であり、同じ等式を得る。

弦のエネルギー密度は $T/\sqrt{1-\omega^2r^2}$、角運動量密度は $T\omega r^2/\sqrt{1-\omega^2r^2}$ である。$u=\arcsin(\omega r)$ を使って積分し、二端の粒子を足すと
$$
\begin{aligned}
E&=\frac{2T}{\omega}\arcsin v+\frac{2m}{\sqrt{1-v^2}},\\
J&=\frac{T}{\omega^2}\bigl(\arcsin v-v\sqrt{1-v^2}\bigr)+
\frac{2m v^2}{\omega\sqrt{1-v^2}} .
\end{aligned}
$$
$T,\omega$ を固定し $m\to0$ とすると $v\to1$、$m/\sqrt{1-v^2}=\sqrt{mT/(\omega v)}\to0$ である。従って $E\to\pi T/\omega$、$J\to\pi T/(2\omega^2)$ で、:eqref[eq-string-classical-regge] に戻る。

:::example[端の速度が半分の模型]
$m\omega/T=3/2$ なら $1-v^2=3v/2$ の正根は $v=1/2$ である。弦のエネルギーは $\pi T/(3\omega)$、二端のエネルギーは $2\sqrt3\,T/\omega$ になる。端の粒子を無視した開弦の Regge 軌跡ではこれらの和を表せない。
:::

### 古典軌跡の最大角運動量を調べる

剛体回転の解が作れるだけでは、与えられた質量に対して最も角運動量の大きい解かどうかは分からない。そこで正のエネルギー密度を持つ共形座標で時間を $X^0=A\tau$（$A>0$）とし、空間の運動量がゼロの系を考える。開弦の全方向が Neumann なら、:eqref[eq-string-classical-open-modes] から重心以外の空間成分は $\int_0^\pi \vb X_{\rm osc}\dd\sigma=0$ を満たす。また $\vb X'^2+\dot{\vb X}^2=A^2$ であるので
$$
\int_0^\pi(\vb X'^2+\dot{\vb X}^2)\dd\sigma=\pi A^2,\qquad E=\pi TA .
$$
開弦の $\sigma$ に対する正の最小波数は 1 であるから、余弦展開で $\int |\vb X_{\rm osc}|^2\dd\sigma\le\int|\vb X'|^2\dd\sigma$ である。これは各係数について $1\le n^2$ を足したものになる。

平面への射影を添字 $\parallel$ で表し、その面内の $\pi/2$ 回転を $R_{\pi/2}$ とする。$J_{12}=T\int (R_{\pi/2}\vb X_{{\rm osc},\parallel})\cdot\dot{\vb X}_\parallel\dd\sigma$ である。絶対値を取り Cauchy–Schwarz を使えば
$$
|J_{12}|\le
T\left(\int|\vb X'|^2\dd\sigma\int|\dot{\vb X}|^2\dd\sigma\right)^{1/2}
\le \frac{T}{2}\int(|\vb X'|^2+|\dot{\vb X}|^2)\dd\sigma
=\alpha' E^2 .
$$
最後の不等式は二つの正数の幾何平均と算術平均の関係である。等号には運動が選んだ面にあり、$\vb X_{\rm osc}$ が第一モードのみを含み、$\dot{\vb X}=\pm R_{\pi/2}\vb X_{\rm osc}$ でなければならない。従って、重心の位置と面内の位相を選べば剛体回転解を得る。$2\pi$ 周期の閉弦でも、定数を除いた最初の波数は 1 である。同じ計算で $\pi$ を $2\pi$ に替えると $|J_{12}|\le\alpha'E^2/2$ になる。折れた閉弦がこの上限に達する。

:::example[高調波だけで回すと傾きが下がる]
開弦で $X^0=A\tau$、$X^1+iX^2=(A/n)e^{in\tau}\cos(n\sigma)$ とし、整数 $n\ge1$ を取る。波動方程式と Neumann 条件を満たし、$\dot X\cdot X'=0$、$\dot X^2+X'^2=0$ である。$E=\pi TA$ は $n$ に依存しないが、$J_{12}=T(A^2/n)\int_0^\pi\cos^2(n\sigma)\dd\sigma=\alpha'E^2/n$ になる。$n=2$ ではエネルギーが等しい第一モードの半分の角運動量を持つ。この解は内部にも折り返しがあり、その点では共形計量の方程式で記述する。
:::

### 回転を慣性系の時刻で測る

剛体回転解の標的時刻は $t=A\tau$ なので、観測される角速度は $\omega=1/A$ である。半径 $r$ の小片の横速度は $\omega|r|$ であり、座標 $r$ で運動量を計算すると単位長さ当たりのエネルギーは $T/\sqrt{1-\omega^2r^2}$、角運動量は $Tr^2\omega/\sqrt{1-\omega^2r^2}$ になる。$r=A\sin u$ を代入すると
$$
E=TA\int_{-\pi/2}^{\pi/2}\dd u,\qquad
J_{12}=TA^2\int_{-\pi/2}^{\pi/2}\sin^2u\,\dd u .
$$
結果は :eqref[eq-string-classical-regge] と一致する。両端の速度を先に小さいと近似すると、端に近い部分のエネルギーの増え方を失う。回転の全体を比較する場合は光的な端も含む積分を使う。

:::example[半径の半分が担う角運動量]
$|r|<A/2$ の内側では $|u|<\pi/6$ である。そこが担う $J_{12}$ は $TA^2[\pi/6-\sqrt3/4]$ であり、全体に対する比は $1/3-\sqrt3/(2\pi)$ である。エネルギーの比は $1/3$ であった。角運動量には半径二乗の重みも掛かるため、中心に近い部分の寄与はさらに小さい。
:::

## Polyakov 作用と世界面計量

### 補助計量と共形ゲージ

:::definition[Polyakov 作用]{#def-string-polyakov}
独立な Lorentz 計量 $\gamma_{ab}$ を使う作用を
$$
S_{\rm P}=-\frac T2\int\dd^2\sigma\,\sqrt{-\gamma}\gamma^{ab}h_{ab}
\label{eq-string-polyakov}
$$
と定める。
:::

$\gamma^{ab}$ を変分すると $\delta\sqrt{-\gamma}=-\sqrt{-\gamma}\gamma_{ab}\delta\gamma^{ab}/2$ であるから、$h_{ab}-\gamma_{ab}\gamma^{cd}h_{cd}/2=0$ を得る。$\gamma$ と $h$ の時間方向が合い、$h$ が非退化なら $h_{ab}=\lambda\gamma_{ab}$、$\lambda>0$ である。代入すると $\sqrt{-h}=\lambda\sqrt{-\gamma}$ より :eqref[eq-string-ng] に戻る。

二次元では $\gamma\mapsto e^{2\omega}\gamma$ のもとで $\sqrt{-\gamma}\gamma^{ab}$ が不変である。これが :term[Weyl 変換]{en="Weyl transformation" yomi="わいるへんかん"} である。二つの光的方向に沿う座標を局所的にとり、Weyl 変換で交差成分を規格化すれば $\gamma=\diag(-1,1)$ にできる。この座標では変分の拘束は
$$
(\dot X+X')^2=(\dot X-X')^2=0 .
\label{eq-string-conformal-constraints}
$$
例えば直線 $X^0=a\tau$、$X^1=b\sigma$ では $a^2=b^2$ が必要である。波動方程式だけなら $a,b$ は独立だが、独立な計量を変分すると両者を結ぶ条件が得られる。

### 正準形式と古典拘束の括弧

面積作用の正準運動量を $P_\mu=P^\tau_\mu$ とする。$h^{\tau a}h_{a\sigma}=0$ より $P\cdot X'=0$ であり、$P^2=T^2(-h)h^{\tau\tau}=-T^2 X'^2$ である。そこで二つの拘束を $\mathcal D=P\cdot X'$、$\mathcal H=(P^2+T^2X'^2)/(2T)$ と書く。点粒子と同様に拘束を Lagrange 乗数で課す一次作用は $\int (P\cdot\dot X-N\mathcal H-N^\sigma\mathcal D)\dd^2\sigma$ である。$P$ の変分は $\dot X=N P/T+N^\sigma X'$ であり、$N=1,N^\sigma=0$ なら $P=T\dot X$ と光的な二条件を再現する。

閉弦で正準 Poisson 括弧 $\{X^\mu(\sigma),P_\nu(\sigma')\}=\delta^\mu_\nu\delta(\sigma-\sigma')$ を採用し、$\mathcal D[f]=\int f\mathcal D\dd\sigma$、$\mathcal H[f]=\int f\mathcal H\dd\sigma$ とする。周期的な $f$ を使えば部分積分で端は現れない。汎関数微分は

$$
\begin{aligned}
\fdv{\mathcal D[f]}{P}&=fX',&
\fdv{\mathcal D[f]}{X}&=-(fP)',\\
\fdv{\mathcal H[f]}{P}&=fP/T,&
\fdv{\mathcal H[f]}{X}&=-T(fX')'.
\end{aligned}
$$

これを $\{A,B\}=\int(\delta_X A\cdot\delta_P B-\delta_P A\cdot\delta_X B)\dd\sigma$ に代入する。$\{\mathcal D[f],\mathcal D[g]\}$ では $P'$ の二項が消え、$\mathcal D[fg'-f'g]$ が残る。$\{\mathcal H[f],\mathcal H[g]\}$ でも $X''$ の二項が消え、同じ $\mathcal D[fg'-f'g]$ を得る。混合括弧では
$-fg P'\cdot P/T+fgT X'\cdot X''-f'gP^2/T+fg'T X'^2$ を積分する。前半の二項を部分積分すると $\{\mathcal D[f],\mathcal H[g]\}=\mathcal H[fg'-f'g]$ になる。例えば $f=1,g=\sin\sigma$ ならこの混合括弧は $\mathcal H[\cos\sigma]$ である。拘束どうしの括弧は再び拘束であり、量子化ではこの代数に中心項が現れる。

### 左右の拘束と古典 Witt 代数

$\mathcal C_\pm=(\mathcal H\pm\mathcal D)/2=(P\pm TX')^2/(4T)$ とする。上の三つの括弧に代入すれば

$$
\{\mathcal C_\pm[f],\mathcal C_\pm[g]\}
=\pm\mathcal C_\pm[fg'-f'g],\qquad
\{\mathcal C_+[f],\mathcal C_-[g]\}=0.
$$

具体的に同符号では二つの $\mathcal D$ 項と、混合括弧の二つの $\pm\mathcal H$ 項が残る。異符号ではそれぞれが逆符号で打ち消す。$f=e^{im\sigma}$、$g=e^{in\sigma}$ なら $fg'-f'g=i(n-m)e^{i(m+n)\sigma}$ であるので、Poisson 括弧は Fourier モードの番号を加える。この代数を :term[Witt 代数]{en="Witt algebra" yomi="うぃっとだいすう"} と呼ぶ。

:::example[拘束が変える座標]
$\{X,\mathcal D[f]\}=fX'$ は $\sigma\to\sigma+f$ の無限小変換である。$\{X,\mathcal H[f]\}=fP/T$ は共形ゲージでは $f\dot X$ になる。$f=1$ なら、それぞれ円周と時間の一様な並進である。$\{P,\mathcal D[f]\}=(fP)'$ なので $P$ は $X$ と違い密度として変換する。
:::

:::example[周期的に振動する円形閉弦]
$X^0=A\tau$、$X^1+iX^2=A\cos\tau\,e^{i\sigma}$ とする。$\ddot X=X''$ であり、$\dot X\cdot X'=0$、$\dot X^2+X'^2=-A^2+A^2\sin^2\tau+A^2\cos^2\tau=0$ である。$\tau=\pi/2$ では半径がゼロになり、誘導計量も全周で縮退する。共形計量は固定されたままなので、Polyakov の波動方程式ではその時刻をまたいで解を続けられる。
:::

### 補助計量を正準作用の乗数に分ける

独立計量を
$$
\gamma_{ab}\dd\sigma^a\dd\sigma^b
=-n^2\dd\tau^2+\rho^2(\dd\sigma+s\dd\tau)^2,\qquad n,\rho>0
$$
と書く。ここで $n$ は時間方向の間隔、$s$ は隣り合う等時面の空間座標のずれを指定する。逆計量の成分は $\gamma^{\tau\tau}=-n^{-2}$、$\gamma^{\tau\sigma}=s/n^2$、$\gamma^{\sigma\sigma}=\rho^{-2}-s^2/n^2$、体積密度は $n\rho$ である。従って :eqref[eq-string-polyakov] の密度は
$$
\mathcal L_{\rm P}
=\frac{T}{2}\left[\frac{\rho}{n}(\dot X-sX')^2-\frac{n}{\rho}X'^2\right].
$$
ここから $P=T(\rho/n)(\dot X-sX')$ を求めて Legendre 変換すると、Hamilton 密度は $(n/\rho)\mathcal H+s\mathcal D$ になる。従って古典の一次形式の $N=n/\rho$、$N^\sigma=s$ は独立な計量の成分から得られる。Weyl 変換で $n,\rho$ は同じ因子を受け、$s$ は変わらない。したがって比 $n/\rho$ と $s$ を固定してから運動方程式を書く場合も、その二つの変分が与える拘束を別に課す。

:::example[斜交座標で同じ共形拘束を解く]
$N=1$、$s=s_0$ を定数とする。Hamilton 方程式より $\dot X-s_0 X'=P/T$ なので、$\mathcal C_\pm=0$ は $[\dot X+( \pm1-s_0)X']^2=0$ である。新しい空間変数 $\hat\sigma=\sigma+s_0\tau$ では $\partial_\tau|_{\hat\sigma}=\partial_\tau|_\sigma-s_0\partial_\sigma$ となるため、共形拘束は $\hat\tau=\tau$、$\hat\sigma$ における :eqref[eq-string-conformal-constraints] に戻る。この座標変換は閉弦で周期を保つ。開弦では座標の境界を動かすので、端の位置も同時に変換する。
:::

### 閉弦の進行波と周期性をそろえる

共形ゲージの閉弦を $X=F(u)+G(v)$、$u=\tau+\sigma$、$v=\tau-\sigma$ と書く。標的時空を非コンパクトとする。$X(\tau,\sigma+2\pi)=X(\tau,\sigma)$ であるため、$F(u+2\pi)-F(u)$ と $G(v)-G(v-2\pi)$ は共通の定数ベクトル $B$ に等しい。これは $u,v$ を独立に変えても両辺が等しくなければならないためである。$P^\mu=T\int_0^{2\pi}\dot X^\mu\dd\sigma=2TB^\mu$ より $B=\pi\alpha'p$ である。従って一般の周期解は
$$
X=x+\alpha'p\tau
+i\sqrt{\frac{\alpha'}{2}}\sum_{n\ne0}\frac{1}{n}\left(
\alpha_n e^{-inu}+\tilde\alpha_n e^{-inv}
\right)
\label{eq-string-classical-closed-modes}
$$
である。$F',G'$ は独立な周期関数だが、その平均は同じ $\alpha'p/2$ である。それぞれの Fourier 係数の実条件は $\alpha_{-n}=\bar\alpha_n$、$\tilde\alpha_{-n}=\bar{\tilde\alpha}_n$ である。

拘束は $F'(u)^2=G'(v)^2=0$ である。$X^0=A\tau$ のゲージを採れば $F'^0=G'^0=A/2$ なので、空間の二曲線を単位速さにそろえて
$$
\vb X(\tau,\sigma)=\frac{A}{2}[\vb a(u)+\vb b(v)],\qquad
|\vb a'|=|\vb b'|=1
$$
と書ける。静止系の閉弦では $\vb a',\vb b'$ の平均はともにゼロとなり、積分した $\vb a,\vb b$ をそれぞれ周期的に取れる。そこで
$$
\dot{\vb X}=\frac A2(\vb a'+\vb b'),\qquad
\vb X'=\frac A2(\vb a'-\vb b')
$$
である。二つの単位ベクトルが一致する点では $\vb X'=0$、$|\dot{\vb X}|=A$ となり、光的な速度に達して誘導計量が退化する。逆に $\vb a'=-\vb b'$ ではその点の速度がゼロとなり、$|\vb X'|=A$ である。

:::example[滑らかな周期解から縮退点を作る]
標的の空間三方向で $\vb a'(u)=(\cos u,\sin u,0)$、$\vb b'(v)=(\cos v,0,\sin v)$ とする。いずれも平均ゼロの単位ベクトルである。$u=v=0$ では両者が $(1,0,0)$ に一致し、その世界面点で $\vb X'=0$ となる。少し離れた点では一般に一致しないので、円形閉弦が同時に一周分縮む例とは異なり、縮退は孤立して生じる。写像の成分は滑らかで波動方程式を満たす。退化点で面積作用の逆計量を用いる代わりに、独立計量の方程式で解を扱う。
:::

:::example[同じ向きに動く零モードと周期性]
$\vb a',\vb b'$ に共通の平均 $\vb m$ があると、それぞれの積分は $2\pi$ 当たり $2\pi\vb m$ 増える。$\sigma$ の一周期では $u$ と $v$ が逆向きに動くため、その増分は相殺する。反対に二つの平均の差を残すと非コンパクトな標的で弦が閉じない。$X^0=A\tau$ では平均 $\vb m$ は全運動量と全エネルギーの比 $\vb p/E$ になる。
:::

### 共形座標の存在と使える範囲

二次元 Lorentz 計量で時空の向きと時間の向きを局所的に選ぶと、光的な線の方向は二つある。それぞれの滑らかな方向場に沿って積分曲線を引く。互いに横切る二つの族に対し、一方の葉で一定な座標を $u$、もう一方で一定な座標を $v$ と選ぶ。座標線が光的なので $\gamma_{uu}=\gamma_{vv}=0$ で、非退化性より $\gamma_{uv}\ne0$ である。平坦な Minkowski 計量に合わせて向きを選べば $\dd s^2=-e^{2\omega}\dd u\dd v$ となり、Weyl 変換で $e^{2\omega}$ を除ける。この構成は局所の主張であり、世界面の貼り合わせ方や周期については別に条件を課す。

$u\mapsto f(u)$、$v\mapsto g(v)$ なら光的座標の性質を保ち、$f',g'>0$ の場合は計量の変化を Weyl 因子で戻せる。開弦の二端を $\sigma=0,\pi$ に保つには一方の境界で $f=g$、もう一方で $f(u+2\pi)=f(u)+2\pi$ が必要である。閉弦で $\sigma$ の周期を $2\pi$ とした場合には $f$ と $g$ にそれぞれ同様の周期条件を課す。ゲージを選んだ後にもこれらの変換が残るため、横振動と縦振動を区別するには残余の変換まで調べる。具体的な時間方向の選択は [光円錐ゲージ量子化](/physics/string-theory/light-cone-quantization/) で扱う。

### 初期値から時間ゲージの解を作る

共形座標の初期面 $\tau=0$ において、$X^0=0$、$\dot X^0=A>0$、$\vb X=\vb f(\sigma)$、$\dot{\vb X}=\vb v(\sigma)$ を与える。閉弦なら $\vb f$ と $\vb v$ は $2\pi$ 周期である。空間成分の波動方程式の解は
$$
\vb X(\tau,\sigma)=
\frac{\vb f(\sigma+\tau)+\vb f(\sigma-\tau)}2+
\frac12\int_{\sigma-\tau}^{\sigma+\tau}\vb v(x)\dd x .
$$
この式を微分すれば、波動方程式と初期値がそのまま確かめられる。時間成分は $X^0=A\tau$ である。二つの拘束を初期面で課す条件は
$$
\vb v\cdot\vb f'=0,\qquad |\vb v|^2+|\vb f'|^2=A^2.
$$
実際、$\dot{\vb X}+\vb X'=\vb v(\sigma+\tau)+\vb f'(\sigma+\tau)$、$\dot{\vb X}-\vb X'=\vb v(\sigma-\tau)-\vb f'(\sigma-\tau)$ である。従って初期面の条件が成り立つと、各平方は特性に沿って $A^2$ のままになる。閉弦で一周期の始点を変えても式は同じである。全方向 Neumann の開弦では、初期値を反射で偶に延長して同じ公式を用いる。すると $\sigma=0,\pi$ で $X'=0$ が保たれる。

任意の滑らかな空間曲線と横速度からこの初期値を作ることもできる。曲線を媒介する変数を $s$、$|\vb v(s)|<1$ を満たす物理速度を曲線に直交するように与える。$|\partial_s\vb f|\ne0$ とし、
$$
\frac{\dd\sigma}{\dd s}
=\frac{|\partial_s\vb f|}{A\sqrt{1-|\vb v(s)|^2}}
$$
で新しい座標を定義する。右辺は正である。$\tau=t/A$ での初速度は $A\vb v$、座標 $\sigma$ での接線の長さは $A\sqrt{1-|\vb v|^2}$ になるので、上の二条件を満たす。閉弦で $\sigma$ の長さを $2\pi$ とするには $A=(2\pi)^{-1}\int |\partial_s\vb f|\dd s/\sqrt{1-|\vb v|^2}$ を選ぶ。これは $E=2\pi T A$ と、面積作用から得た静的時刻でのエネルギーが一致することも示す。自由端で速度が 1 になる解ではこの局所構成を端そのものに適用せず、内部から延長する。

:::example[楕円を静止させて離す]
平面上の楕円を初期形状に取り、各点を静止させる。周長を $L$ とすると $\sigma$ は初期の弧長を $2\pi/L$ 倍したもので、$A=L/(2\pi)$ である。解は $\vb X(\tau,\sigma)=[\vb f(\sigma+\tau)+\vb f(\sigma-\tau)]/2$ となる。中心を原点に置いた楕円の弧長表示では $\vb f(\sigma+\pi)=-\vb f(\sigma)$ であるから、$\tau=\pi/2$ では二項が打ち消し合い、全点が原点に来る。その標的時刻は $t=L/4$ である。円だけでなく、初期形状がこの半周期の反転対称性を持てば同じ議論が成立する。
:::

:::example[初期速度が直交しないときに起きること]
円を初期形状とし、その接線に沿って各点に同じ非零の速度を割り当てると、$\vb v\cdot\vb f'$ は非零になる。これを共形時間ゲージの初期速度としてそのまま用いることはできない。面積作用では、曲線の上でラベルを動かすことにより初期速度に接線方向を足せるので、まず接線速度を除く。この例の物理的な横速度はゼロであり、時間ゲージでの初期値は静止した円になる。
:::

### Weyl 対称性が働く次元

独立計量を導入する方法を $d_{\mathrm w}$ 次元の世界体積に拡張する。張力を $T_{\mathrm w}$ とし、
$$
S=-\frac{T_{\mathrm w}}2\int\dd^{d_{\mathrm w}}\sigma\,\sqrt{-\gamma}\,
\left[\gamma^{ab}h_{ab}-(d_{\mathrm w}-2)\right]
$$
を考える。$\gamma^{ab}$ の変分は
$$
h_{ab}-\frac12\gamma_{ab}\bigl[\gamma^{cd}h_{cd}-(d_{\mathrm w}-2)\bigr]=0
$$
を与える。$d_{\mathrm w}\ne2$ なら両辺のトレースを取って $\gamma^{ab}h_{ab}=d_{\mathrm w}$ が得られるため、$\gamma=h$ に固定される。従って非退化な場合には面積に相当する世界体積の作用に戻る。世界体積が二次元のときはトレースの式が恒等式になり、計量の全体の尺度が未定となる。これは運動項の $\sqrt{-\gamma}\gamma^{ab}$ が Weyl 変換で $e^{(d_{\mathrm w}-2)\omega}$ 倍になることと一致する。

:::example[世界線と三次元の世界体積]
$d_{\mathrm w}=1$ で $\gamma_{\lambda\lambda}=-e^2$ と書くと、補助作用は $(T_{\mathrm w}/2)\int(e^{-1}\dot X^2-e)\dd\lambda$ である。$T_{\mathrm w}=m$、$e_{\rm particle}=e/m$ とすれば :eqref[eq-string-particle-einbein] と同じ形になる。$d_{\mathrm w}=3$ では運動項と定数項がそれぞれ $e^\omega$ と $e^{3\omega}$ 倍になるため、二次元と同じ方法で全体の尺度を除くことはできない。
:::

### 再パラメータ化による拘束の伝播

共形ゲージで二つの拘束を $C_\pm=(\dot X\pm X')^2$ とする。波動方程式より $(\partial_\tau\mp\partial_\sigma)C_\pm=0$ であるから、初期条件として与えた $C_\pm=0$ は特性に沿って保たれる。開弦の端では $\dot X\cdot X'=0$ より $C_+=C_-$ であるため、片側の特性が端へ到達すると反射側に同じ拘束値が渡る。従って、内部の初期値と境界の両立を確かめれば、時間発展の途中で新しい独立拘束を足す必要はない。

### 自由端の光的運動を拘束から直接得る

全方向 Neumann の開弦を共形ゲージに置き、端で $X'^\mu=0$ とする。:eqref[eq-string-conformal-constraints] を端へ連続的に延長すると $\dot X^2=0$ になる。標的時刻が $X^0=A\tau$ と増える場合、端の物理速度は $|\dot{\vb X}|/A=1$ である。従って自由端の光的な運動は剛体回転に固有ではない。端に外力や質量があれば $P^\sigma$ の境界条件が変わるため、同じ議論では端の速度を固定できない。自由端の点で $h_{\tau\sigma}=h_{\sigma\sigma}=0$ となるので、補助計量を持たない式の逆計量を端で直接評価することはできない。内部では $\sqrt{-h}h^{ab}$ が共形座標で一定となり、その極限から $P^\sigma$ の境界条件を読む。

:::example[自由端に見える縦の座標速度]
時間ゲージで二つの特性が合流する端を取り、端の直近で $X'(\tau,\sigma)$ が一次でゼロに近づくとする。拘束を使うと $1-|\dot{\vb X}|^2/A^2=|\vb X'|^2/A^2$ なので、速度の 1 との差は端までの座標距離の二次である。一方 $P^{\tau0}=TA$ は曲線の長さ当たりでなく座標 $\sigma$ 当たりの密度であり、有限のままである。弧長の測度へ変えると $\dd E/\dd s=TA/|\vb X'|$ となり、弧長に対しては端で平方根型の発散が出る。端に近い弦のエネルギーは発散せず、どの測度で密度を表したかで見え方が変わる。
:::

## 演習問題

::::exercise[符号]{#exr-string-classical-string-1 level=basic}
$X^0=2\tau$、$X^1=3\sigma$ の面積密度を求めよ。

:::solution
誘導計量は $\diag(-4,9)$ なので $-\det h=36$、密度は 6 である。
:::
::::

::::exercise[固定端]{#exr-string-classical-string-2 level=basic}
両端を $Y(0)=1$、$Y(\pi)=3$ に固定した静的な波動方程式の解を求めよ。

:::solution
$Y''=0$ より $Y=a+b\sigma$ である。端点条件から $Y=1+2\sigma/\pi$ となる。
:::
::::

::::exercise[基本振動]{#exr-string-classical-string-3 level=basic}
間隔 $L=4$ の固定端で第 2 モードの角振動数を求めよ。

:::solution
固定端では波数 $k_n=n\pi/L$ であり、波動方程式から $\omega=k_n$ なので $\pi/2$ となる。
:::
::::

::::exercise[角運動量]{#exr-string-classical-string-4 level=basic}
回転解で $A=4\ell_s$ の $J_{12}$ を求めよ。

:::solution
:eqref[eq-string-classical-regge] を使う。$E=2/\ell_s$ なので $J_{12}=4$ である。
:::
::::

::::exercise[局所的な進行波]{#exr-string-classical-string-5 level=advanced}
単一の横変位 $Y=f(t+x)$ が静的ゲージで厳密解になることを示せ。

:::solution
$\dot Y=Y'=f'$ だから平方根の中は 1 である。変分方程式はこの解上で $\partial_t f'-\partial_x f'=f''-f''=0$ となる。
:::
::::

::::exercise[混合した境界条件]{#exr-string-classical-string-6 level=advanced}
波動方程式で $Y(0)=0$、$Y'(\pi)=0$ の最低の正の角振動数を求めよ。

:::solution
空間部分を $A\sin(k\sigma)+B\cos(k\sigma)$ とすると $B=0$ である。$\cos(k\pi)=0$ から $k=n+1/2$ であり、最低値は $1/2$ となる。$k=0$ の一次解は両条件からゼロである。
:::
::::

::::exercise[回転解の計量]{#exr-string-classical-string-7 level=advanced}
章の回転解で $h_{\tau\sigma}=0$ を確かめよ。

:::solution
二次元の空間ベクトルを複素数で書くと、$\partial_\tau X=iAe^{i\tau}\cos\sigma$、$\partial_\sigma X=-Ae^{i\tau}\sin\sigma$ である。内積は共役と積をとった実部なのでゼロである。時間成分も $X^{0\prime}=0$ より寄与しない。
:::
::::

::::exercise[Weyl 自由度]{#exr-string-classical-string-8 level=advanced}
$\gamma=e^{2\omega}h$ を Polyakov 作用に代入し、$\omega$ が消えることを示せ。

:::solution
二次元なので $\sqrt{-\gamma}=e^{2\omega}\sqrt{-h}$、$\gamma^{ab}h_{ab}=2e^{-2\omega}$ である。積が $2\sqrt{-h}$ となり、面積作用を得る。
:::
::::

::::exercise[初期形状から進行波を作る]{#exr-string-classical-string-9 level=advanced}
張った弦の線形近似で、$0\le x\le L$ の両端を固定する。初期値 $Y(0,x)=0$、$\dot Y(0,x)=v_0\sin(2\pi x/L)$ に対する解と二次の余剰エネルギーを求めよ。

:::solution
$k=2\pi/L$ とおくと、$\sin(kx)$ は端でゼロであり空間の演算子 $-\partial_x^2$ の固有値は $k^2$ である。従って $Y(t,x)=A\sin(kt)\sin(kx)+B\cos(kt)\sin(kx)$ の形で初期値を合わせられる。$B=0$、$kA=v_0$ より $Y=(v_0/k)\sin(kt)\sin(kx)$ となる。二つの解の差で初期値がゼロなら、同じ二次エネルギーの保存と正値性から解はゼロであり、この解は一意である。実際、エネルギーの時間微分は部分積分すると端の $T[\dot Y Y']$ だけであり、固定端では $\dot Y=0$ である。初期時刻に運動項だけが残るので $\Delta E=(T/2)v_0^2\int_0^L\sin^2(kx)\dd x=Tv_0^2L/4$ である。
:::
::::

## 確認テスト

::::quiz{answer="4" tol="1e-6"}
共形条件を満たす静的直線で $X^0=4\tau$ とする。$X^1=b\sigma$ の $b>0$ を求めよ。

:::explanation
拘束は $-4^2+b^2=0$ であり、正の根は 4 である。
:::
::::

::::quiz
自由端で消えるものを選べ。

- [ ] 並進カレントの時間成分
- [x] 運動量の法線方向への流束
- [ ] 張力

:::explanation
境界項は $P^\sigma\delta X$ である。$\delta X$ が任意なので係数の $P^\sigma$ がゼロになる。
:::
::::

::::quiz{answer="9" tol="1e-6"}
同じ張力で回転解のエネルギーを 3 倍にすると角運動量は何倍か。

:::explanation
$J=\alpha'E^2$ を使う。$\alpha'$ は張力が同じなので一定であり、比は $3^2=9$ である。
:::
::::

::::quiz
共形座標での拘束に含まれるものをすべて選べ。

- [x] $\dot X\cdot X'=0$
- [x] $\dot X^2+X'^2=0$
- [ ] $\dot X^2-X'^2=0$

:::explanation
二つの光的条件の差をとると交差項がゼロになり、和をとると二乗の和がゼロとなる。
:::
::::

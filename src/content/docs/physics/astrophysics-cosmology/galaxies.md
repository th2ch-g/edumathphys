---
title: 銀河の力学
description: 銀河のポテンシャルにおける星の運動を分布関数で記述し、回転曲線と渦状構造を力学から解釈する。
notation:
  - symbol: '\Phi'
    meaning: 重力ポテンシャル
  - symbol: 'f'
    meaning: 位相空間分布関数
  - symbol: '\kappa'
    meaning: エピサイクル角周波数
  - symbol: '\beta'
    meaning: 速度分散の異方性
  - symbol: 'v_c'
    meaning: 円速度
---

## この章の目標

- ポテンシャルと軌道の安定性を計算する。
- 無衝突の運動方程式から Jeans 方程式を導く。
- 緩和、暗黒物質、密度波の観測への影響を区別する。

## ポテンシャル理論と質量分布

体積密度 $\rho$ に対し $\nabla^2\Phi=4\pi G\rho$ と置く。孤立系では $\Phi(\vb x)=-G\int\rho(\vb y)/|\vb x-\vb y|\dd^3y$ である。核 $1/r$ は原点以外で調和で、小球上の法線微分の面積積分が $-4\pi$ なのでこの積分は Poisson 式を満たす。球対称なら原点で $r^2\Phi'\to0$ を用いて $r^2\Phi'=GM(<r)$ を得る。軸対称な薄円盤は球殻定理に従わない。厚さをまたいで Poisson 式を積分すると、面密度 $\Sigma$ に対して $\partial_z\Phi(0^+)-\partial_z\Phi(0^-)=4\pi G\Sigma$ となる。

:::example[一様球]
内部の $\Phi'=4\pi G\rho r/3$ により円速度は半径に比例する。半径二倍で円速度も二倍である。
:::


:::example[平面の力]
反射対称なら直上の加速度は $-2\pi G\Sigma$ である。$\Sigma$ が二倍の無限平面では加速度も二倍になる。
:::

### 有限質量に正規化できる球状模型

Plummer ポテンシャルを $\Phi=-GM/\sqrt{r^2+b^2}$、$b>0$ とする。球対称で $m=r^2\Phi'/G$ だから $m=Mr^3/(r^2+b^2)^{3/2}$ である。これを微分して
$$
\rho(r)=\frac{3Mb^2}{4\pi(r^2+b^2)^{5/2}}
$$
を得る。中心で有限、遠方で $r^{-5}$ となり、総質量は $M$ である。視線方向を $z$、空の半径を $R$ と置き、$z=\sqrt{R^2+b^2}\tan u$ を使うと $\int_{-\infty}^\infty(1+t^2)^{-5/2}\dd t=4/3$ である。従って面密度は $\Sigma(R)=Mb^2/[\pi(R^2+b^2)^2]$ となる。

:::example[射影半質量半径]
$M_{\rm proj}(<R)=2\pi\int_0^R\Sigma(R')R'\dd R'=MR^2/(R^2+b^2)$ である。$R=b$ は射影質量の半分を囲む一方、球半径 $b$ の質量は $M/(2\sqrt2)$ である。
:::

### 軸対称円盤の面密度から力を求める

円盤が無限に薄く、$\rho(R,z)=\Sigma(R)\delta(z)$ と表せる場合、$z\ne0$ では $\Phi$ が Laplace 方程式を満たす。円筒座標の軸対称 Laplace 演算子の動径部分を $L_R=R^{-1}\partial_R(R\partial_R)$ とする。Bessel 方程式により $L_R J_0(kR)=-k^2J_0(kR)$ である。両側で減衰する解と Hankel 変換を
$$
\widehat\Sigma(k)=\int_0^\infty R'\Sigma(R')J_0(kR')\dd R',
\qquad
\Sigma(R)=\int_0^\infty k\widehat\Sigma(k)J_0(kR)\dd k
$$
とする。これは平面上の二次元 Fourier 変換の角度を積分したものである。実際、角度平均 $(2\pi)^{-1}\int_0^{2\pi}e^{ikR\cos\varphi}\dd\varphi$ は $J_0(kR)$ である。従って Poisson 式の解は
$$
\Phi(R,z)=-2\pi G\int_0^\infty\widehat\Sigma(k)J_0(kR)e^{-k|z|}\dd k .
$$
地上と地下のそれぞれで Laplace 方程式を満たし、$\partial_z\Phi$ の飛びは $4\pi G\Sigma(R)$ になる。境界で零になる解を取れば、別の解との差は全空間で調和かつ無限遠で零となるため、最大値原理で差は零である。円盤面で $J_0'(x)=-J_1(x)$ を用いれば
$$
v_c^2(R)=2\pi G R\int_0^\infty k\widehat\Sigma(k)J_1(kR)\dd k
$$
となる。面密度が内側の円に含まれる部分だけで決まらないことは、変換の積分範囲が円盤全体に及ぶことから分かる。

:::example[Kuzmin 円盤]
$b>0$ とし、$\Phi(R,z)=-GM/\sqrt{R^2+(b+|z|)^2}$ とする。$z>0$ では点質量を $z=-b$ に置いた調和関数であり、$z<0$ では鏡映である。面での微分の飛びは $2GMb/(R^2+b^2)^{3/2}$ なので、$\Sigma=Mb/[2\pi(R^2+b^2)^{3/2}]$ である。面での円速度二乗は $GMR^2/(R^2+b^2)^{3/2}$ となり、$R=\sqrt2 b$ で最大となる。一方、面密度の累積は $M[1-b/\sqrt{R^2+b^2}]$ である。$R=b$ の円盤速度に $Rv_c^2/G$ を適用すると $M/(2\sqrt2)$ になるが、実際の円盤面内の質量は $M(1-1/\sqrt2)$ である。
:::

:::example[指数円盤の変換]
$\Sigma=\Sigma_0e^{-R/R_d}$ とする。$a=R_d^{-1}$ とおくと $\int_0^\infty e^{-aR}J_0(kR)\dd R=(a^2+k^2)^{-1/2}$ である。角度積分の形を用いると左辺は $(2\pi)^{-1}\int (a-ik\cos\varphi)^{-1}\dd\varphi$ であり、虚部は相殺する。実部を四分円へ縮め、$u=\tan\varphi$ とすると $(2/\pi)\int_0^\infty a/(a^2(1+u^2)+k^2)\dd u$ となり、右辺が得られる。$a$ で微分して $\widehat\Sigma(k)=\Sigma_0 R_d^2[1+(kR_d)^2]^{-3/2}$ を得る。$k=0$ の $2\pi\widehat\Sigma(0)$ は総円盤質量であり、$R=2R_d$ の内側の含有率は $1-3e^{-2}\simeq0.594$ である。
:::

### 球面調和関数で球状からのずれを扱う

ポテンシャルと密度を規格化した球面調和関数 $Y_{\ell m}$ によって $\Phi=\sum_{\ell m}\Phi_{\ell m}(r)Y_{\ell m}$、$\rho=\sum_{\ell m}\rho_{\ell m}(r)Y_{\ell m}$ と展開する。球面上の Laplace 演算の固有値は $-\ell(\ell+1)$ なので、径の方程式は
$$
\frac{\dd}{\dd r}(r^2\Phi_{\ell m}')-\ell(\ell+1)\Phi_{\ell m}
=4\pi G r^2\rho_{\ell m}.
$$
斉次解は $r^\ell$ と $r^{-\ell-1}$ である。中心で正則、外で減衰する解を作るため、$r$ の内側の源に外側の斉次解、外側の源に内側の斉次解を掛ける。$r=r'$ で解が連続し、$r^2\Phi'$ に強さ $4\pi G r'^2\rho_{\ell m}(r')\dd r'$ の差があるとすると
$$
\Phi_{\ell m}(r)=-\frac{4\pi G}{2\ell+1}
\left[r^{-\ell-1}\int_0^r\rho_{\ell m}(s)s^{\ell+2}\dd s
+r^\ell\int_r^\infty\rho_{\ell m}(s)s^{1-\ell}\dd s\right].
$$
係数の $2\ell+1$ は、接続する二本の斉次解の対数勾配が $\ell$ と $-\ell-1$ である差から出た。コンパクトな密度分布の外側では第二項がなくなり、最も遅く減る $\ell=0$ が総質量を表す。分布の中心を重心に置くと、$\ell=1$ は質量一次モーメントなので零である。従って定常な非球形を遠方から区別する先頭項は一般に $\ell=2$ である。

:::example[外側のポテンシャルが球状に近づく速さ]
外側の単極子は $r^{-1}$、四重極子は $r^{-3}$ である。半径を四倍にすると四重極と単極の比は十六分の一になる。よって内側の円盤速度に球状の換算を使う際の誤差が大きくても、全質量の十分外側では形の影響が減る。この議論は有限の外端を仮定し、外側にもハロー質量が多く分布する範囲にはそのまま使えない。
:::

### Abel の逆変換で観測面から立体分布へ戻る

球対称な数密度 $\nu(r)$ が距離 $R$ の位置に作る投影 $\Sigma(R)$ は、視線積分で
$$
\Sigma(R)=2\int_R^\infty\frac{\nu(r)r\,\dd r}{\sqrt{r^2-R^2}}
$$
である。減衰が速く、必要な微分と積分を交換できる滑らかな分布を考える。$u=R^2$、$s=r^2$ に替えると $S(u)=\Sigma(\sqrt u)=\int_u^\infty \nu(\sqrt s)\dd s/\sqrt{s-u}$ である。これに $(u-v)^{-1/2}$ を掛け $u\ge v$ で積分し、順を交換すると内積分 $\int_v^s[(u-v)(s-u)]^{-1/2}\dd u=\pi$ が出る。残りは $\pi\int_v^\infty\nu(\sqrt s)\dd s$ である。$u=v+w$ として左辺を $v$ で微分し、物理半径へ戻すと
$$
\nu(r)=-\frac1\pi\int_r^\infty\frac{\dd\Sigma}{\dd R}
\,\frac{\dd R}{\sqrt{R^2-r^2}} .
$$
投影に微分が掛かるので、観測点のノイズを直接差分すると振動しやすい。球対称が疑わしいときは、投影を同じ $\Sigma$ に変える三次元密度が一つに定まらない。

:::example[射影が Gaussian の場合]
$\Sigma(R)=\Sigma_0e^{-R^2/(2b^2)}$ とする。$\dd\Sigma/\dd R=-R\Sigma/b^2$ を入れ $z=\sqrt{R^2-r^2}$ を使うと
$$
\nu(r)=\frac{\Sigma_0 e^{-r^2/(2b^2)}}{\pi b^2}
\int_0^\infty e^{-z^2/(2b^2)}\dd z
=\frac{\Sigma_0}{\sqrt{2\pi}b}e^{-r^2/(2b^2)} .
$$
三次元 Gaussian を視線に沿って普通に積分しても元の面密度に戻る。中心の密度は同じ表面密度でも奥行き $b$ の逆数に依存する。
:::

## 銀河内の星の軌道

時間に依存しない球対称ポテンシャルでは、エネルギー $E=v^2/2+\Phi$ と比角運動量 $j$ が保存され、$\dot r^2/2+\Phi+j^2/(2r^2)=E$ となる。円軌道は有効ポテンシャルの停留点、安定性は二階微分で判定する。軸対称円盤の面内でも $j=R^2\Omega$ 一定の半径摂動から

$$
\kappa^2=\Phi''+3\Phi'/R=R\,\frac{\dd\Omega^2}{\dd R}+4\Omega^2,\quad \Omega^2=\Phi'/R . \label{eq-cosmo-epicycle}
$$

この小振動を :term[エピサイクル運動]{en="epicyclic motion" yomi="えぴさいくるうんどう"} という。$v_c\propto R^s$ の場合は $\kappa^2=2(s+1)\Omega^2$ であり、$s>-1$ で面内の小振動が安定である。

:::example[点質量と平坦回転]
Kepler 極限は $s=-1/2$ なので $\kappa/\Omega=1$、平坦なら $\sqrt2$ である。
:::


:::example[剛体回転]
$s=1$ で $\kappa=2\Omega$ となり、軌道一周に径方向の小振動二周期が入る。
:::

### エピサイクルによる歳差

ほぼ円軌道で径方向の振動周期は $2\pi/\kappa$、平均方位角速度は $\Omega$ である。一径方向周期の間の方位角増分は $2\pi\Omega/\kappa$ である。中心力の模型で局所的に円速度 $v_c\propto R^s$ とすれば、$\Omega^2\propto R^{2s-2}$ だから :eqref[eq-cosmo-epicycle] は $\kappa^2=2(s+1)\Omega^2$ となる。$s>-1$ の範囲では小さい径方向変位が振動する。

:::example[閉軌道と非閉軌道]
点質量で $s=-1/2$ なら $\kappa=\Omega$ で近地点は一周ごとに戻る。平坦な回転曲線では $s=0$、一径周期の方位角は $\sqrt2\pi\simeq4.44\,\mathrm{rad}$ であり、円に近い軌道の近地点は一周で閉じない。
:::

### 動径周期と方位角の進みを積分する

球対称の束縛軌道で $j>0$、二つの単根 $r_p<r_a$ を持つとする。$u(r)=\sqrt{2(E-\Phi(r))-j^2/r^2}$ と書くと、動径速度は行きで $+u$、戻りで $-u$ である。従って動径周期と一周期の方位差は
$$
T_r=2\int_{r_p}^{r_a}\frac{\dd r}{u},
\qquad
\Delta\phi=2\int_{r_p}^{r_a}\frac{j\dd r}{r^2u}.
$$
根の近傍で $u$ は根までの距離の平方根に比例するため、端点の特異性は積分できる。遠点と近点の両方で数値的に安定な変数として $r=(r_a+r_p)/2-(r_a-r_p)\cos\theta/2$ を使えばよい。分子の $\sin\theta$ が端点での平方根の零を相殺する。

$J_r=\pi^{-1}\int_{r_p}^{r_a}u\,\dd r$ を :term[動径作用]{en="radial action" yomi="どうけいさよう"} と呼ぶ。積分端で $u=0$ なので $E,j$ で微分したとき端点の項は消える。$\partial_Eu=1/u$、$\partial_ju=-j/(r^2u)$ を使って
$$
\pdv{J_r}{E}=\frac{T_r}{2\pi},\qquad
\pdv{J_r}{j}=-\frac{\Delta\phi}{2\pi}
$$
となる。従って $(J_r,j)$ を固定点の近傍で独立変数にすれば、$\partial E/\partial J_r$ は動径周波数、$J_r$ 一定の $\partial E/\partial j$ は平均方位角速度である。二つの角周波数の比が有理数なら、軌道は有限回の径振動の後に閉じる。

:::example[等方調和ポテンシャル]
$\Phi=\omega^2r^2/2$ なら直交座標の各成分が同じ $\omega$ で振動する。軌道面で時刻と座標軸の原点を選べば $x=A\cos\omega t,y=B\sin\omega t$ と書ける。$A\ne B$、$A,B>0$ のとき $r^2$ は $2\omega$ で振動し、一径周期で方位が $\pi$ 進む。二径周期で閉じるので、$\kappa=2\Omega$ の円軌道近似と一致する。この模型では近円軌道の外でも同じ比が成立する。
:::

### 鉛直運動と軸対称性

定常な軸対称ポテンシャルで Lagrangian は $(\dot R^2+R^2\dot\phi^2+\dot z^2)/2-\Phi(R,z)$ である。$\phi$ は循環座標なので $L_z=R^2\dot\phi$ が一定となり、$(R,z)$ は有効ポテンシャル $\Phi+L_z^2/(2R^2)$ に従う。対称面 $z=0$ で $\Phi(R,z)=\Phi(R,-z)$ なら、$\Phi_{Rz}(R,0)=0$ である。従って円軌道では、一次で径振動と鉛直振動が分離する。鉛直周波数の二乗は $\nu_z^2=\Phi_{zz}(R,0)$ である。薄い面密度がある位置には力の跳びがあり、その点の鉛直運動を滑らかな調和振動で近似できない。滑らかな対称面では Poisson 式と円速度の定義から
$$
\nu_z^2=4\pi G\rho(R,0)-\frac1R\dv{v_c^2}{R}.
$$
この計算は鉛直周波数から面の体積密度を推定できる理由を示す。回転曲線が平坦なら第二項は消える。

:::example[平坦な球と同じ回転曲線の円盤]
球ポテンシャルなら $\partial_z\Phi=\Phi'(r)z/r$ なので、$z=0$ で $\nu_z^2=\Phi'(R)/R=\Omega^2$ である。球を平たくした模型が同じ $v_c$ を持つ場合でも、$\rho(R,0)$ の違いで $\nu_z$ が異なり得る。従って鉛直速度のデータは、回転曲線だけでは分けにくい質量分布を制約する。
:::

### 潮汐が星団の束縛を弱める尺度

質量 $M_s$ の小さな衛星が、滑らかな球ポテンシャル $\Phi_h(R)$ の半径 $R_0$ で円運動するとする。衛星の質量は背景に比べ小さいとし、$\Omega^2=\Phi_h'(R_0)/R_0$ を保つ。衛星とともに回る座標の径方向に $x$ だけ離れる点で、背景の重力と遠心力の差を展開する。零次は円運動の条件で消え、一次は $[\Omega^2-\Phi_h''(R_0)]x$ である。点質量としての衛星の力 $-GM_s x/|x|^3$ とこの項が釣り合う位置は
$$
r_t=\left[\frac{G M_s}{\Omega^2-\Phi_h''(R_0)}\right]^{1/3}
=\left[\frac{GM_s}{4\Omega^2-\kappa^2}\right]^{1/3}
$$
である。分母が正で $r_t\ll R_0$ の場合の近似であり、その等ポテンシャル面が力学的な境界の尺度を与える。正確な脱出時間には星の位相、衛星の形、軌道の離心率も効く。コアに近い調和ポテンシャルでは分母が零となり、この一次の式で有限の境界を得られない。

:::example[主銀河を点質量で近似する]
$\Phi_h=-GM_h/R$ なら $\Phi_h''=-2\Omega^2$ であり、$r_t=R_0[M_s/(3M_h)]^{1/3}$ である。$M_s/M_h=3\times10^{-6}$ なら $r_t/R_0=0.01$ である。同じ衛星を平坦な回転曲線の位置に置き、局所の円速度と囲む質量を揃えた場合、分母は $2\Omega^2$ となって半径は $(3/2)^{1/3}$ 倍になる。
:::

### 軸対称でも運動が三つの積分を必要とする理由

軸対称な定常場ではエネルギー $E$ と角運動量 $L_z$ は必ず保存する。その二つを指定した後の運動は有効ポテンシャル中の $(R,z)$ の二自由度である。$E,L_z$ に依存する分布は定常であるが、同じ二量を持つ軌道の鉛直運動や径運動は必ずしも同じでない。分離できる例として $\Phi(R,z)=\Phi_R(R)+\Phi_z(z)$ なら、$E_z=\dot z^2/2+\Phi_z(z)$ が独立に保存する。実際、$\dot E_z=\dot z(\ddot z+\Phi_z')=0$ である。従って $f(E,L_z,E_z)$ とすれば、同じ全エネルギーでも薄い軌道と厚い軌道の重みを変えられる。完全に分離できない場合でも、安定な円軌道の近くでは $J_z\simeq E_z/\nu_z$ と $J_r\simeq E_r/\kappa$ が先頭で独立な作用になる。

:::example[低い鉛直周波数で厚く見える軌道]
調和鉛直運動なら最大高度 $z_{\rm max}$ は $E_z=\nu_z^2 z_{\rm max}^2/2$ から決まる。固定の $E_z$ で $\nu_z$ が半分なら高度は二倍、平均する一周期も二倍である。固定する量が断熱的な $J_z$ の場合は $z_{\rm max}=\sqrt{2J_z/\nu_z}$ なので、高度は $\sqrt2$ 倍になる。外側へゆっくり軌道が移る場合と異なる初期条件の軌道を比べる場合で、どの量を固定するかが違う。
:::

### ほぼ円の軌道の面内の形を求める

基準円を $R=R_g$、$L_z=R_g^2\Omega_g$ とする。小さい径変位を $x=A\cos(\kappa t+\theta_0)$ として、保存則から方位角速度は $\dot\phi=L_z/(R_g+x)^2\simeq\Omega_g(1-2x/R_g)$ である。円とともに回る系の接線距離 $y=R_g(\phi-\Omega_gt)$ は、位相の原点を揃えると
$$
y=-\frac{2\Omega_g}{\kappa}A\sin(\kappa t+\theta_0)
$$
である。従って回転系の小軌道は楕円で、接線方向の半径は径方向の $2\Omega_g/\kappa$ 倍である。局所の円速度から測った方位速度のずれは、回転系の $\dot y$ と違う。物理方位速度は $v_\phi=L_z/(R_g+x)\simeq R_g\Omega_g-\Omega_gx$ であり、その点の円速度 $v_c(R_g+x)$ を引くと
$$
v_\phi-v_c(R_g+x)\simeq-[2\Omega_g+R_g\Omega_g']x
=-\frac{\kappa^2}{2\Omega_g}x .
$$
位相が一様で振幅の分布も狭い局所近似では $\sigma_\phi^2/\sigma_R^2=\kappa^2/(4\Omega_g^2)$ となる。使ったのは速度が小さい星であり、異なるガイド中心を持つ広い分布の平均を正確に決めた式ではない。

:::example[平坦な曲線でのエピサイクルの偏り]
$\kappa=\sqrt2\Omega_g$ なので回転系の軌道の長軸は接線方向で径方向の $\sqrt2$ 倍である。一方、局所の円速度に対する速度分散は $\sigma_\phi^2/\sigma_R^2=1/2$ となる。座標の振幅が大きい方向と、円運動からの速度のずれが大きい方向が一致する必要はない。
:::

## 無衝突 Boltzmann 方程式と Jeans 方程式

星が一つずつの遭遇でほとんど曲がらない場合、任意相空間体積で個数を数え、速度場 $(\vb v,-\nabla\Phi)$ の発散が零であることを使う。この保存を分布 $f$ について表したものが :term[無衝突 Boltzmann 方程式]{en="collisionless Boltzmann equation" yomi="むしょうとつぼるつまんほうていしき"} である。

$$
\partial_t f+\vb v\cdot\nabla_x f-\nabla\Phi\cdot\nabla_v f=0 \label{eq-cosmo-collisionless}
$$

これは :eqref[eq-neq-liouville] を一粒子の平均場へ適用した式である。$\nu=\int f\dd^3v$ とし、$1$ と $v_i$ を掛けて積分する。速度端で境界項が消えれば $\partial_t\nu+\partial_j(\nu\bar v_j)=0$、$\partial_t(\nu\bar v_i)+\partial_j(\nu\overline{v_iv_j})=-\nu\partial_i\Phi$ となる。球対称、定常、平均流なしで応力の発散の径成分を取ると

$$
\dv{(\nu\sigma_r^2)}r+\frac{2\beta}r\nu\sigma_r^2=-\nu\Phi',\qquad \beta=1-\frac{\sigma_\theta^2+\sigma_\phi^2}{2\sigma_r^2}. \label{eq-cosmo-jeans-spherical}
$$

一定の $\beta$ なら積分因子 $r^{2\beta}$ で解ける。観測が主に視線速度なので、異方性と質量を同時に一意に決められないことがある。

:::example[等温トレーサー]
平坦な円速度 $v_c$、$\nu\propto r^{-3}$、$\beta=0$ の一定分散解は $\sigma_r^2=v_c^2/3$ である。
:::


:::example[位置に依らない保存量]
$f=F(E)$ は軌道上で $E$ が一定なので定常で、速度に関して球対称である。
:::

### 視線速度に投影して観測と比べる

半径 $r$ の球対称で定常な星系を、投影半径 $R$ に沿って観測する。視線方向を $\ell$ とすれば $r^2=R^2+\ell^2$ で、前後を足して $\Sigma_*(R)=2\int_R^\infty\nu(r)r\dd r/\sqrt{r^2-R^2}$ となる。球対称の分散テンソルに対し、視線と半径の間の角度の二乗余弦は $1-R^2/r^2$ である。接線方向の分散は $\sigma_\theta^2=\sigma_\phi^2=(1-\beta)\sigma_r^2$ なので、局所での視線分散は $\sigma_r^2(1-\beta R^2/r^2)$ になる。従って
$$
\Sigma_*(R)\sigma_{\rm los}^2(R)=2\int_R^\infty
\frac{\nu(r)\sigma_r^2(r)(1-\beta R^2/r^2)r}{\sqrt{r^2-R^2}}\dd r .
$$
同じ質量分布でも $\beta$ を変えると :eqref[eq-cosmo-jeans-spherical] の径分散と投影因子の双方が変わる。測定した $\sigma_{\rm los}$ から質量を求めるには別の運動成分やトレーサーを合わせる理由がここにある。

:::example[平坦曲線で等方な等速分散]
$\nu\propto r^{-3}$、$v_c$ 一定、$\beta=0$ の初等解を用いる。投影式の被積分で $\sigma_r^2=v_c^2/3$ を外へ出せるので、$\sigma_{\rm los}=v_c/\sqrt3$ は投影半径によらない。
:::

### 密度から等方的な分布を復元する

孤立した球が静止し、脱出運動の星を含まないとする。相対ポテンシャル $\psi(r)=-\Phi(r)+\Phi(\infty)\ge0$ と相対結合エネルギー $\mathcal E=\psi-v^2/2$ を置く。等方な分布 $f(\mathcal E)$ は $\mathcal E\le0$ で零とする。$\nu=4\pi\int_0^{\sqrt{2\psi}}v^2 f(\psi-v^2/2)\dd v$ を $\mathcal E$ へ変えると
$$
\nu(\psi)=4\pi\sqrt2\int_0^\psi f(\mathcal E)\sqrt{\psi-\mathcal E}\dd\mathcal E .
$$
微分し、$\dd\nu/\dd\psi$ を $1/\sqrt{\mathcal E-\psi}$ で重み付けして $0<\psi<\mathcal E$ で積分する。積分順序を交換したとき、内側の $\int_u^\mathcal E \dd\psi/\sqrt{(\mathcal E-\psi)(\psi-u)}$ は $\pi$ となる。$\psi=u+(\mathcal E-u)\sin^2\theta$ で確かめられる。よって一回微分すれば
$$
f(\mathcal E)=\frac1{\sqrt8\pi^2}\frac{\dd}{\dd\mathcal E}
\int_0^\mathcal E\frac{\dd\nu/\dd\psi}{\sqrt{\mathcal E-\psi}}\dd\psi .
$$
これが :term[Eddington の逆変換]{en="Eddington inversion" yomi="えでぃんとんのぎゃくへんかん"} である。滑らかな単調 $\psi(r)$ から $\nu(\psi)$ を作る。計算した $f$ が負なら、与えた密度とポテンシャルはこの等方な平衡の仮定を満たさない。

:::example[べき密度]
$\nu(\psi)=A\psi^2$ なら積分は $2A\int_0^\mathcal E\psi(\mathcal E-\psi)^{-1/2}\dd\psi=(8A/3)\mathcal E^{3/2}$ となる。よって $f(\mathcal E)=\sqrt2 A\sqrt{\mathcal E}/\pi^2$ である。
:::

### 積分因子で Jeans 方程式を解く

異方性を $\beta(r)$ とし、$J(r)=\exp[\int_{r_*}^r 2\beta(s)\dd s/s]$ と置く。:eqref[eq-cosmo-jeans-spherical] に掛けると $\partial_r[J\nu\sigma_r^2]=-J\nu\Phi'$ である。外端 $r_b$ での値を $Q_b=(J\nu\sigma_r^2)(r_b)$ とすれば
$$
\sigma_r^2(r)=\frac{Q_b+\int_r^{r_b}J(s)\nu(s)\Phi'(s)\dd s}{J(r)\nu(r)} .
$$
無限遠に零の境界を課すには積分が収束し、$J\nu\sigma_r^2$ が零へ行く必要がある。分散が正でも、それを持つ位相空間分布 $f\ge0$ の存在は別に確かめる。

:::example[分散と異方性の縮退]
平坦回転 $v_c$、$\nu\propto r^{-\gamma}$、一定 $\beta$ で $\gamma>2\beta$ なら上の無限遠解は $\sigma_r^2=v_c^2/(\gamma-2\beta)$ となる。$\gamma=3$ で $\beta=0$ と $\beta=1/2$ を比べると、同じ $v_c$ に対し分散は $v_c^2/3$ と $v_c^2/2$ になる。
:::

### Jeans 方程式と分散の輸送

速度分散テンソルを $C_{ij}=\overline{v_iv_j}-\bar v_i\bar v_j$ と定める。連続方程式を用いると、一次モーメントは
$$
(\partial_t+\bar v_j\partial_j)\bar v_i
=-\partial_i\Phi-\nu^{-1}\partial_j(\nu C_{ij})
$$
となる。これを得るには二次モーメントを $\bar v_i\bar v_j+C_{ij}$ に分け、積の微分中の $\bar v_i[\partial_t\nu+\partial_j(\nu\bar v_j)]$ を消す。二次の重み $v_iv_j$ を無衝突方程式に掛けて積分すると、その式には $\partial_k(\nu\overline{v_iv_jv_k})$ が出る。一般の衝突しない分布では三次モーメントが二次だけでは決まらず、流体の局所平衡の状態式を自動的に使えない。定常な分布を保証する一つの方法は、独立な運動の積分 $I_1,\ldots,I_s$ の関数 $f=F(I_1,\ldots,I_s)$ を取ることである。連鎖律で軌道に沿う $\dot f$ が零になる。さらに自己重力で閉じるなら、この $f$ の速度積分から定まる密度も Poisson 式と一致させる。

球対称で平均流がない場合、$Q=\nu C$ は球の直交基底で $Q_{rr}\hat{\vb r}\hat{\vb r}+Q_t(I-\hat{\vb r}\hat{\vb r})$ である。$\nabla\cdot\hat{\vb r}=2/r$、$(\hat{\vb r}\cdot\nabla)\hat{\vb r}=0$ より $(\nabla\cdot Q)_r=Q_{rr}'+2(Q_{rr}-Q_t)/r$ である。$Q_t=(1-\beta)Q_{rr}$ を入れれば :eqref[eq-cosmo-jeans-spherical] を再現する。

:::example[Plummer 球の等方な速度分散]
自己重力 Plummer 模型で $\nu$ は質量密度に比例し、$\Phi'=GMr/(r^2+b^2)^{3/2}$ である。外側で $\nu\sigma_r^2=0$ を境界条件とすると、$\nu\Phi'$ は $r(r^2+b^2)^{-4}$ に比例し、その積分は $(r^2+b^2)^{-3}/6$ である。従って $\sigma_r^2=GM/(6\sqrt{r^2+b^2})$ となる。同じ球で $\nu\propto\psi^5$ だから Eddington 式の積分は $\mathcal E^{9/2}$ に比例し、$f\propto\mathcal E^{7/2}\ge0$ である。従って分散だけでなく、正の分布関数を持つ平衡模型でもある。
:::

### 全体の運動エネルギーに対するビリアル関係

粒子の質量を $m_i$ とし $I=\sum_i m_i|\vb x_i|^2$、$T=\sum_i m_i|\vb v_i|^2/2$ とする。二回微分すると $\ddot I/2=2T+\sum_i\vb x_i\cdot\vb F_i$ である。孤立した Newton 重力の対ごとの項は
$$
\vb x_i\cdot\vb F_{ij}+\vb x_j\cdot\vb F_{ji}
=-\frac{Gm_im_j}{|\vb x_i-\vb x_j|}.
$$
従って力の総和は重力のポテンシャルエネルギー $W$ である。有限な $I$ を持つ定常系、あるいは長時間平均で端点の $\dot I$ の差を時間幅で割ると零になる系では $2\langle T\rangle+\langle W\rangle=0$ を得る。有限の観測領域だけを取れば星が境界を通過し、応力の表面項が残り得る。固定された外場中のトレーサーなら、力の項は $-\int \rho_* \vb r\cdot\nabla\Phi\,\dd^3x$ として計算する。

:::example[一様球の見積もり]
一定密度の球で $M(<r)=M(r/R)^3$ とすれば、殻を組み上げる仕事の総和は $W=-\int_0^R GM(<r)\dd M/r=-3GM^2/(5R)$ である。圧力や粒子運動の境界項が零で、全体の三方向の速度二乗平均が等しければ $2T=3M\sigma^2$ だから $\sigma^2=GM/(5R)$ となる。ここで $\sigma^2$ は全体平均であり、どの半径でも同じ分散を仮定する必要はない。
:::

### 軸対称の Jeans 式と恒星の回転補正

軸対称な定常分布で $\bar v_R=\bar v_z=0$ とする。円筒座標の速度成分を使うと、一粒子の動径加速度は $\dot v_R=v_\phi^2/R-\partial_R\Phi$、$\dot v_\phi=-v_Rv_\phi/R$ である。分布の運動式を $v_R$ で重み付けして速度にわたり積分する。速度空間の部分積分から $\nu(\overline{v_R^2}-\overline{v_\phi^2})/R$ を得て、
$$
\partial_R(\nu\sigma_R^2)+\partial_z(\nu\overline{v_Rv_z})
+\frac{\nu(\sigma_R^2-\overline{v_\phi^2})}{R}
=-\nu\partial_R\Phi.
$$
ここで $v_c^2=R\partial_R\Phi$、$\overline{v_\phi^2}=\bar v_\phi^2+\sigma_\phi^2$ とおくと
$$
v_c^2-\bar v_\phi^2=\sigma_R^2\!\left[-\partial_{\log R}\log(\nu\sigma_R^2)
-1+\frac{\sigma_\phi^2}{\sigma_R^2}\right]
-\frac R\nu\partial_z(\nu\overline{v_Rv_z})
$$
である。:term[非対称ドリフト]{en="asymmetric drift" yomi="ひたいしょうどりふと"} は、円速度と平均方位速度の差である。速度分散の径勾配だけでなく、速度楕円の傾きの鉛直変化もこの差へ入る。対称面で $\overline{v_Rv_z}=0$ となっても、鉛直微分が消えるとは限らない。

:::example[指数分布の遅い回転]
$\nu\propto e^{-R/R_d}$、$\sigma_R$ 一定、鉛直の混合項を零、$\sigma_\phi^2/\sigma_R^2=1/2$ と仮定する。$R=3R_d$ では $v_c^2-\bar v_\phi^2=(5/2)\sigma_R^2$ である。$v_c=220\,\mathrm{km\,s^{-1}}$、$\sigma_R=40\,\mathrm{km\,s^{-1}}$ なら $\bar v_\phi\simeq211\,\mathrm{km\,s^{-1}}$ である。分散の大きい恒星と冷たいガスの回転速度を同じ円速度へ変換するには、それぞれの支持力を入れる。
:::

### 速度積分を持つ平衡分布を実際に正規化する

Plummer 球の質量密度を $\rho=A\psi^5$、$A=3b^2/(4\pi G^5M^4)$、$\psi=GM/\sqrt{r^2+b^2}$ と書く。$f(\mathcal E)=F_0\mathcal E^{7/2}$ を仮定し、速度積分式に入れて $\mathcal E=u\psi$ とする。すると
$$
\rho=4\pi\sqrt2 F_0\psi^5
\int_0^1 u^{7/2}(1-u)^{1/2}\dd u
=\frac{7\sqrt2\pi^2}{64}F_0\psi^5 .
$$
最後の積分は Beta 積分 $B(9/2,3/2)=7\pi/256$ である。従って $F_0=24\sqrt2 b^2/(7\pi^3G^5M^4)$ になる。$f\ge0$ で質量密度を再現し、運動に沿う $\mathcal E$ が一定なので、この式は自己重力の定常解を与える。

:::example[中心で許される速さ]
中心の相対ポテンシャルは $GM/b$ であるから、速度密度は $v^2(GM/b-v^2/2)^{7/2}$ に比例し、$v\ge\sqrt{2GM/b}$ では零である。$x=v^2/(2GM/b)$ で書くと関数は $x(1-x)^{7/2}$ なので極大は $x=2/9$ にある。この極大は速さ $v$ 当たりの分布であり、$x$ 当たりの確率密度に変数変換すると Jacobian のため位置が変わる。
:::

### Jeans 方程式で観測分散から質量を読む限界

球で :eqref[eq-cosmo-jeans-spherical] を $GM(<r)/r^2=-\nu^{-1}\partial_r(\nu\sigma_r^2)-2\beta\sigma_r^2/r$ と解く。$\nu,\sigma_r^2>0$ の区間で対数微分を使うと
$$
M(<r)=-\frac{r\sigma_r^2}{G}\left[
\frac{\dd\log\nu}{\dd\log r}+
\frac{\dd\log\sigma_r^2}{\dd\log r}+2\beta\right].
$$
右辺が正であるには括弧が負である必要がある。また質量が外へ向けて減るなら、分布の非負性や測定・仮定に疑いがある。一方、射影の速度だけから径方向の分散を得るには異方性が要る。等方なら積分式で $\Sigma_*\sigma_{\rm los}^2$ が $\nu\sigma_r^2$ の Abel 投影なので、測光と速度二つの投影を独立に逆変換してから代入できる。異方性が既知でない場合は重み $1-\beta R^2/r^2$ がさらに未知となる。

:::example[べき則での正質量条件]
$\nu\propto r^{-\gamma}$、$\sigma_r^2\propto r^{-\alpha_s}$、$\beta$ 一定なら $M(<r)=(\gamma+\alpha_s-2\beta)r\sigma_r^2/G$ である。$\gamma+\alpha_s>2\beta$ で正とし、$\alpha_s>1$ の広い区間を持つ仮定を置くと、推定質量は $r^{1-\alpha_s}$ で減少する。独立に当てはめた分散と密度が、それだけで物理的な平衡模型になるとは限らない具体例である。
:::

### 定常な分布を軌道ごとに定める

定常な平均場で :eqref[eq-cosmo-collisionless] は、$f$ が Hamilton 軌道に沿って一定という条件である。従って定常解は、位相空間内で一つの軌道に属する点へ同じ値を割り当てる。特に、正則な作用・角変数 $(\vb J,\boldsymbol\theta)$ を持つ可積分領域を考え、$\dot{\vb J}=0$、$\dot{\boldsymbol\theta}=\boldsymbol\Omega(\vb J)$ とする。連続な定常解を角度で Fourier 展開したとき、$\boldsymbol\Omega\cdot\vb m\ne0$ のモードの係数は $\boldsymbol\Omega\cdot\nabla_\theta f=0$ により零となる。よって全ての零でない整数ベクトル $\vb m$ に対しこの不等式が成り立つトーラスでは $f$ は角度に依存せず、作用だけで表せる。この結論には可積分性と非共鳴の仮定がある。閉軌道や共鳴したトーラスでは、軌道上で一定でも一定の角度の結合に依存する余地が残る。さらに自己無撞着な銀河模型なら、得られた分布の速度積分を重力の Poisson 方程式へ入れたときに同じ平均場を再現する必要がある。

:::example[角変数の比が有理数の場合]
二角の周波数が $\Omega_2=2\Omega_1$ なら、$2\theta_1-\theta_2$ は軌道上で一定である。共鳴トーラスに限れば $\cos(2\theta_1-\theta_2)$ を含む分布も定常条件を満たす。周波数が少しずれてこの組合せが零でなくなると、その分布は同じ軌道上に沿って変化し続ける。従って整った近円軌道があることだけから、一般の定常分布に全ての角度を自由に入れることはできない。
:::

## 緩和時間と力学的摩擦

星一個の質量を $m$、相対速度を $v$ とする。衝突径数 $b$ で直線近似すると横速度キックは $2Gm/(bv)$ である。遭遇率 $2\pi b\,\dd b\,nv$ にその二乗を掛けると $\dd\langle\Delta v^2\rangle/\dd t\simeq8\pi G^2m^2n\log\Lambda/v$、$\Lambda=b_{\max}/b_{\min}$ である。速度が典型値分だけ拡散する時間は $t_{\rm rel}\sim v^3/(8\pi G^2m^2n\log\Lambda)$ となる。自己重力系の尺度 $v^2\sim GNm/R$ と $n\sim N/R^3$ を代入すると $t_{\rm rel}/t_{\rm cross}\sim N/(8\pi\log\Lambda)$ である。大質量天体は背景に後方過密を作り、進行方向と逆の重力を受ける。等方一様な Maxwell 背景での :term[Chandrasekhar の近似]{en="Chandrasekhar dynamical friction" yomi="ちゃんどらせかーるのきんじ"} は $\dot{\vb v}=-4\pi G^2(M+m)\rho\log\Lambda\,\mathcal F(<v)\vb v/v^3$ である。$\mathcal F$ は速さ $v$ 未満の質量分率で、その速度方向を含む積分は後で導出する。

:::example[遭遇の数]
$\log\Lambda$ 固定で密度が二倍なら緩和時間は半分になる。
:::


:::example[摩擦の質量依存]
$M\gg m$ で同じ軌道速度と背景なら、衛星質量を十倍にすると摩擦加速度も十倍になる。
:::

### 小角散乱から速度空間の積分へ

試験天体 $M$ と背景粒子 $m$ の相対速度を $\vb w=\vb v-\vb u$ とする。Newton の双曲線散乱で $\tan(\chi/2)=G(M+m)/(bw^2)$ である。これは比エネルギー $w^2/2$ と比角運動量 $bw$ から離心率 $e^2=1+b^2w^4/[G^2(M+m)^2]$ が出て、漸近線から $\sin(\chi/2)=1/e$ となるためである。重心系で試験体の速度は $m\vb w/(M+m)$ なので、方位角平均した一回の速度変化は $-m(1-\cos\chi)\vb w/(M+m)$ である。衝突頻度 $2\pi b\dd b\,w$ を掛ける。$b_0=G(M+m)/w^2$ とすれば $1-\cos\chi=2b_0^2/(b^2+b_0^2)$ であり、$b$ の積分は $\frac12\log[(b_{\max}^2+b_0^2)/(b_{\min}^2+b_0^2)]$ を与える。

速度に依存するこの対数を $\log\Lambda$ で近似すると、残る方向積分は $\int f(\vb u)(\vb v-\vb u)|\vb v-\vb u|^{-3}\dd^3u$ である。等方 $f$ なら速度空間の球殻定理により $|\vb u|>v$ の殻が消え、遅い粒子の数 $N(<v)$ を用いて $N(<v)\vb v/v^3$ になる。質量密度を $\rho$ と表せば先の Chandrasekhar 式を得る。この導出では背景の軌道・自己重力の集団応答を無視した。

:::example[遅い天体の減速]
滑らかな等方 $f$ で $v$ が小さければ $N(<v)=4\pi f(0)v^3/3+\order{v^5}$ である。従って加速度は先頭で $-\vb v$ に比例する。大きい速度で得た $v^{-2}$ の加速度の大きさを低速まで延長すると、逆の極限になる。
:::

### Maxwell 背景で摩擦の速度依存を求める

一成分速度分散 $\sigma^2$ の等方 Maxwell 背景では、速さ $v$ 未満の質量分率は
$$
\mathcal F(<v)=\frac4{\sqrt\pi}\int_0^X u^2e^{-u^2}\dd u
=\operatorname{erf}X-\frac{2X}{\sqrt\pi}e^{-X^2},\qquad
X=\frac{v}{\sqrt2\sigma}.
$$
$2u^2e^{-u^2}=e^{-u^2}-\partial_u(ue^{-u^2})$ を使った。$X\ll1$ では $\mathcal F=4X^3/(3\sqrt\pi)+\order{X^5}$ なので加速度の大きさは $v$ に比例する。$X\gg1$ では $\mathcal F$ は一へ近づき、加速度の大きさは $v^{-2}$ に近づく。背景密度と Coulomb 対数が固定という仮定の下での両極限である。

:::example[等温球内の円軌道の減衰]
球ハローが $v_c$ 一定で等方かつ自己重力の特異等温球なら、Jeans 式から $\sigma^2=v_c^2/2$ である。円軌道上を動く一定質量 $M\gg m$ の衛星では $X=1$、$\mathcal F(1)\simeq0.428$ である。$\rho=v_c^2/(4\pi Gr^2)$ なので単位衛星質量当たりのトルクは $-GM\log\Lambda\,\mathcal F(1)/r$ になる。ゆっくり縮む円軌道で $j=rv_c$ を時間微分すると $r\dot r=-GM\log\Lambda\,\mathcal F(1)/v_c$ である。従って $r_i$ から $r_f<r_i$ までの時間は $v_c(r_i^2-r_f^2)/(2GM\log\Lambda\,\mathcal F(1))$ となる。核付近までこの式を適用するには、背景の密度勾配や衛星の潮汐損失が小さいかを別に確かめる。
:::

### 十分に強い散乱と有限の衝突径数

相対速度 $w$ の二体問題で $b_{90}=G(M+m)/w^2$ とすると $\tan(\chi/2)=b_{90}/b$ である。$b\gg b_{90}$ では $\chi\simeq 2b_{90}/b$ から小角の速度キックの積分が対数になる。一方、$b\ll b_{90}$ では $\chi$ は $\pi$ に近づくので、一回の散乱が無限大のキックを与えることはない。質点として数える範囲で正確な衝突積分を行えば、$b_{\min}=0$ でも力学的摩擦に現れる対数は $\frac12\log[1+(b_{\max}/b_{90})^2]$ になる。大きな距離をどこで切るかは、有限な系の局所近似をどこまで使うかで決まる。引力による近接二体を長時間追う場合、多数の弱い遭遇を独立に扱って背景を一様とする仮定も点検する。

:::example[緩和時間と公転の数]
同じ質量の恒星を $N=10^{10}$ 個とし $\log\Lambda=20$ を仮定する。本文の概算で $t_{\rm rel}/t_{\rm cross}\simeq2.0\times10^7$ である。$t_{\rm cross}=5\times10^7\,\mathrm{yr}$ なら $t_{\rm rel}\sim10^{15}\,\mathrm{yr}$ となり、二体緩和は宇宙年齢より長い。一方、$N=10^5$、$\log\Lambda=10$、交差が $10^6\,\mathrm{yr}$ なら同じ見積もりで約 $4\times10^8\,\mathrm{yr}$ となる。この比と時間の規格化は分布に依存するが、銀河と星団で二体緩和の重要度が大きく違う理由を示す。
:::

## 回転曲線と暗黒物質の証拠

ガスの Doppler 速度や星の運動から $v_c(R)$ を求める。薄円盤の視線速度は円運動だけなら $v_{\rm los}-v_{\rm sys}=v_c\sin i\cos\varphi$ である。$\varphi$ は円盤面で主軸から測る方位角である。光がほとんど増えない外側でも平坦なら、球形の近似で

$$
M(<r)=rv_c^2/G,\qquad \rho(r)=v_c^2/(4\pi Gr^2) \label{eq-cosmo-rotation-mass}
$$

となる。二式目は $v_c$ 一定の範囲で質量を微分して得る。可視成分の質量が有限なら十分外で $v_c\propto r^{-1/2}$ が予想されるため、バリオンの分布と合わせて追加質量を推定できる。円盤が非軸対称、傾斜が小さい、圧力支持が効く場合は速度から円速度への変換を点検する。

:::example[外縁比較]
一定の $v_c$ で半径二倍なら内包質量二倍、局所密度は四分の一になる。
:::


:::example[傾斜角誤差]
同じ観測振幅で $i=30^\circ$ と $90^\circ$ を仮定すると前者の円速度は二倍、球質量は四倍になる。
:::

### 円盤の速度を球として逆算する場合

軸対称な銀河では各成分のポテンシャルが線形に足し合わされるため、平面の円速度は $v_c^2=\sum_iR\partial_R\Phi_i=\sum_i v_{c,i}^2$ である。円盤自体が球の重力法則を満たすわけではないが、光から推定した円盤、バルジ、ガスのポテンシャルを先に計算し、余った速度二乗を球ハロー模型へ割り当てることはできる。この分解は恒星の質量光度比に依存する。独立した鉛直運動やレンズの測定は分解の制約を増やす。

:::example[速度を直接差し引く場合との比較]
ある半径で $v_c=200\,\mathrm{km\,s^{-1}}$、バリオン模型の $v_b=120\,\mathrm{km\,s^{-1}}$ とする。残りの球ハローの寄与は $\sqrt{200^2-120^2}=160\,\mathrm{km\,s^{-1}}$ である。$80\,\mathrm{km\,s^{-1}}$ を代入したときの球内質量は正しい残差の四分の一になる。
:::

### 表面輝度を質量へ移す

近赤外の表面輝度を $I(R)$、恒星の帯域別質量光度比を $\Upsilon_*$ とすると、単純な恒星円盤模型は $\Sigma_*=\Upsilon_* I$ である。Poisson 式は密度に線形なので、固定した $I$ のもとで円盤の円速度二乗への寄与は $\Upsilon_*$ に比例する。ガス円盤は原子・分子の放射線から面密度を見積もるが、元素組成や光学的厚さによって線強度から質量への換算が変わる。さらに軸対称な定常ガスの径方向の Euler 式は
$$
v_\phi^2=v_c^2+\frac{R}{\rho}\partial_R P .
$$
円速度と平均回転速度の差には圧力勾配の符号が入る。爆発や流入で径速度が大きい場合、この定常釣り合いも検査が必要である。恒星の円速度へはこれに対応する Jeans 式の速度分散補正を行う。

:::example[圧力補正]
一定の速度分散で $P=c_s^2\rho$、$\rho\propto e^{-R/R_d}$ なら $v_c^2-v_\phi^2=c_s^2R/R_d$ である。$R=3R_d$、$c_s=10\,\mathrm{km\,s^{-1}}$、$v_\phi=30\,\mathrm{km\,s^{-1}}$ では $v_c=\sqrt{1200}\simeq34.6\,\mathrm{km\,s^{-1}}$ となる。球模型で質量を換算する場合、補正なしの値の $4/3$ 倍である。
:::

### 球状ハローの形を回転曲線へ写す

一つの模型として :term[Navarro–Frenk–White 分布]{en="Navarro–Frenk–White profile" yomi="なばろふれんくほわいとぶんぷ"} $\rho(r)=\rho_s/[x(1+x)^2]$、$x=r/r_s$ を使う。$\rho_s,r_s>0$ とする。積分すると
$$
M(<r)=4\pi\rho_sr_s^3m(x),\qquad m(x)=\log(1+x)-\frac{x}{1+x}
$$
である。$m'=x/(1+x)^2$、$m(0)=0$ を確認すれば、質量微分が $4\pi r^2\rho$ に戻る。従って $v_c^2=4\pi G\rho_sr_s^2m(x)/x$ である。$x\ll1$ では $m=x^2/2+\order{x^3}$ なので $v_c\propto r^{1/2}$ である。$x\gg1$ では $m=\log x-1+\order{x^{-1}}$ であり、総質量を有限とするには外側の定義または切断を別に指定する。ある外半径 $r_\Delta$ までの平均密度を $\Delta\rho_{\rm ref}$ とするなら、集中度 $c_\Delta=r_\Delta/r_s$ について $\rho_s=\Delta\rho_{\rm ref}c_\Delta^3/[3m(c_\Delta)]$ である。$\rho_{\rm ref}$ に臨界密度と平均物質密度のどちらを用いたかを記録する必要がある。

:::example[丸い中心を持つ模型との違い]
中心の密度が有限で零でない球では $M(<r)\simeq(4\pi/3)\rho_0r^3$ なので $v_c\propto r$ となる。密度が $r^{-1}$ の中心では $M(<r)\propto r^2$、従って $v_c\propto\sqrt r$ である。最も内側のデータが傾斜角、有限の分解能、棒構造の非円運動によって偏ると、この違いも弱く見える。
:::

### 傾斜したガス円盤の速度場を当てはめる

薄い円盤の見かけの長軸に沿う座標を $x$、短軸に沿う座標を $y$ とする。一定の傾斜 $i$、$0<i<\pi/2$ なら $R^2=x^2+(y/\cos i)^2$、$\cos\varphi=x/R$ である。面内の平均速度を $v_R\hat{\vb e}_R+v_\phi\hat{\vb e}_\phi$ と分解すると、視線速度の符号を一方の主軸で正に選んだとき
$$
v_{\rm los}=v_{\rm sys}+\sin i(v_\phi\cos\varphi+v_R\sin\varphi)
$$
と書ける。従って傾斜角と長軸方向を固定できる場合、主軸の速度は回転、短軸の速度は径向き運動を調べるのに使える。円盤の傾きや主軸が半径ごとに変わると、短軸の信号を径向き速度だけに割り当てられない。$i\to0$ では投影が小さいため、観測角を精密に測るほど重要になる。

:::example[低傾斜での微小誤差]
主軸に沿う測定を固定し $v_\phi\propto1/\sin i$ を対数微分すると $\dd v_\phi/v_\phi=-\cot i\,\dd i$ である。$i=15^\circ$ で $2^\circ$ の誤差は一次で約 $13\%$ の速度誤差、球状の質量換算では約 $26\%$ に相当する。同じ角度幅でも $i=60^\circ$ なら速度誤差は約 $2\%$ である。
:::

### 有限なコアを持つハローと回転の上限

球状で $\rho(r)=\rho_0/[1+(r/r_c)^2]$ とする。$x=r/r_c$ を使い $x^2/(1+x^2)=1-1/(1+x^2)$ を積分すれば
$$
M(<r)=4\pi\rho_0r_c^3(x-\arctan x),\qquad
v_c^2=4\pi G\rho_0r_c^2\left(1-\frac{\arctan x}x\right).
$$
$x\ll1$ で $\arctan x=x-x^3/3+\order{x^5}$ なら内部で剛体回転型へ、$x\gg1$ なら一定速度へ近づく。密度の外側が $r^{-2}$ なので無限遠までの総質量は発散する。観測範囲で平坦な曲線を記述しても、全半径での質量を同じべきで決めたことにはならない。

:::example[速度が一定値へどの速さで近づくか]
$x=1$ では $v_c/v_\infty=\sqrt{1-\pi/4}\simeq0.463$、$x=10$ では約 $0.924$ である。観測から一つの平坦値だけを取ると $4\pi G\rho_0r_c^2$ の積が決まり、$\rho_0$ と $r_c$ が別々に定まらない。上りの区間まで測れば、形の関数がそれらを分ける情報を与える。
:::

## 渦状構造と密度波理論

円盤上の密度の山は星が同じ腕に固定されていなくても伝わる。薄い気体円盤で摂動 $\propto e^{i(kR+m\phi-\omega t)}$、$|kR|\gg m$ を仮定する。共回転周波数 $s=\omega-m\Omega$ に対し、連続式 $s\delta\Sigma=k\Sigma u_R$ と、Coriolis 項を含む速度式 $-isu_R-2\Omega u_\phi=-ik(c_s^2\delta\Sigma/\Sigma+\delta\Phi)$、$-isu_\phi+\kappa^2u_R/(2\Omega)=0$ を使う。薄面の Poisson 式を面外で指数減衰する解に接続し $\delta\Phi=-2\pi G\delta\Sigma/|k|$ を得る。消去すれば

$$
s^2=\kappa^2+c_s^2 k^2-2\pi G\Sigma|k|. \label{eq-cosmo-disk-wave}
$$

軸対称摂動の最小は $|k|=\pi G\Sigma/c_s^2$ にあり、指数成長する局所軸対称モードがない条件は $Q_g=c_s\kappa/(\pi G\Sigma)\ge1$ となる。等号ではある波数で $s^2=0$ となり、限界安定である。星だけの円盤では星の無衝突軌道を積分する必要があり、ここで導いた気体の式を直接適用できない。反応核と非軸対称な増幅は :cite[binney-tremaine] を参照する。

:::example[局所不安定性]
$Q_g=1/2$ のとき最小の $s^2$ は $-3\kappa^2$ である。
:::


:::example[共回転]
パターン速度 $\Omega_p=\omega/m$ に対して $\Omega=\Omega_p$ では $s=0$ となる。腕を横切る角速度差が零になる。
:::

### パターン速度と共鳴

角度依存が $\cos[m(\phi-\Omega_p t)]$ の弱い摂動を、基準の円軌道に沿って観察する。強制振動数は $m|\Omega-\Omega_p|$ である。径変位の線形方程式 $\ddot\xi+\kappa^2\xi=f_0\cos[m(\Omega-\Omega_p)t]$ の特解には分母 $\kappa^2-m^2(\Omega-\Omega_p)^2$ が現れる。$\Omega-\Omega_p=\pm\kappa/m$ を Lindblad 共鳴、$\Omega=\Omega_p$ を共回転という。共回転ではこの径強制振動の分母は発散しないが、パターンが軌道に対して止まるため持続的なトルクが生じうる。

:::example[平坦曲線での二本腕]
$v_c$ 一定、$m=2$ とすると、共回転半径 $R_c=v_c/\Omega_p$ に対し Lindblad 半径は $R_c(1\mp1/\sqrt2)$ となる。内側は約 $0.293R_c$、外側は約 $1.707R_c$ である。この位置は背景回転曲線を指定して初めて決まる。
:::

### パターンの回転座標で保存する量

回転する非軸対称ポテンシャルを $\Phi(R,\phi-\Omega_p t,z)$ とし、一定パターン速度 $\Omega_p$ を仮定する。粒子の比機械エネルギー $E=|\vb v|^2/2+\Phi$ の時間微分は固定位置での $\partial_t\Phi$ である。軸まわりの比角運動量 $L_z$ の変化は $-\partial_\phi\Phi$ なので、$\dot E=\Omega_p\dot L_z$ となる。従って :term[Jacobi 積分]{en="Jacobi integral" yomi="やこびせきぶん"} $E_J=E-\Omega_pL_z$ が保存する。局所で回転速度を差し引いた速度を $\vb v_{\rm rot}$ として平方を展開すれば、$E_J=|\vb v_{\rm rot}|^2/2+\Phi-\Omega_p^2R^2/2$ である。従ってある $E_J$ の軌道は、有効ポテンシャルが $E_J$ 以下の領域だけを動ける。複数のパターン速度が同時に作用すれば、単一の $\Omega_p$ に対する保存性は一般に失われる。

:::example[エネルギーと角運動量の交換]
$\Omega_p=25\,\mathrm{km\,s^{-1}\,kpc^{-1}}$ のパターンから、星が $\Delta L_z=40\,\mathrm{kpc\,km\,s^{-1}}$ の比角運動量を得るとする。保存式から $\Delta E=\Omega_p\Delta L_z=1000\,\mathrm{(km\,s^{-1})^2}$ である。軌道にトルクが働くと慣性系でのエネルギーも変わるが、その二つの変化は同じパターン速度で結ばれる。
:::

### 局所的な腕の伝播と巻き上がり

各星と一緒に動く模様の初期位相を $\phi_0(R)$ とすると、時刻 $t$ の位相は $\phi(R,t)=\phi_0(R)+\Omega(R)t$ である。$t$ とともに $\partial_R\phi=\phi_0'(R)+t\Omega'(R)$ が変化する。腕の接線が方位方向となす角を $i_p$ とすれば $|\tan i_p|=|1/(R\,\partial_R\phi)|$ である。従って差動回転の下では物質とともに動く模様は巻き上がる。密度波のパターンは位相が座標を進むため、同じ粒子が位相の山に留まる仮定が要らない。

気体の局所分散関係で実周波数の枝を $\omega=m\Omega\pm \sqrt{\kappa^2+c_s^2 k^2-2\pi G\Sigma|k|}$ と書けば、$k\ne0$ で波束の径方向群速度は
$$
v_g=\pdv{\omega}{k}=\frac{c_s^2k-\pi G\Sigma\,\operatorname{sgn}k}{\omega-m\Omega}.
$$
群速度は位相の角速度 $\Omega_p$ と別に、波束の包絡が運ぶ擾乱の方向を与える。分母が零に近くなる場合や波長が背景の変化尺度に近づく場合、この局所近似で大域の反射条件を求められない。

:::example[平坦回転曲線の巻き上がり]
初期に腕が径方向で $v_c$ が一定なら、$\Omega'= -v_c/R^2$ だから $|\tan i_p|=1/(\Omega t)$ となる。基準半径の一周後には $|\tan i_p|=1/(2\pi)$ で $i_p$ は約 $9^\circ$ である。観測の腕の形を一つの物質線として長時間保つ場合に問題となる時間尺度を、この計算が与える。
:::

### 衝突しない円盤の自己重力と速度分散

気体の圧力項に恒星の分散を代入すると、軌道の位相に沿う平均が消える。星の円盤について、薄い円盤、近円軌道、局所軸対称波、Gaussian 型の径速度分布を仮定する。$s_*=(\omega-m\Omega)/\kappa$、$\chi=k^2\sigma_R^2/\kappa^2$ とする。無衝突方程式を非摂動のエピサイクルに沿って積分すると、自己重力の分散関係は
$$
s_*^2=1-\frac{2\pi G\Sigma|k|}{\kappa^2}\mathcal F(s_*,\chi),\qquad
\mathcal F(s_*,\chi)=
\frac{2(1-s_*^2)e^{-\chi}}{\chi}
\sum_{n=1}^\infty\frac{I_n(\chi)}{1-s_*^2/n^2}.
$$
$I_n$ は変形 Bessel 関数であり、整数共鳴以外での形式である。この応答核は後で軸対称の場合について導く。共鳴を含む大域的な応答の処方は :cite[binney-tremaine] を参照する。ここでは中立点の条件を求める。母関数 $e^{\chi\cos\theta}=I_0(\chi)+2\sum_{n\ge1}I_n(\chi)\cos(n\theta)$ を $\theta=0$ で評価すれば
$$
\mathcal F(0,\chi)=\frac{1-e^{-\chi}I_0(\chi)}{\chi}.
$$
$u=|k|\sigma_R/\kappa$ とおくと、波数ごとの中立条件は
$$
\frac{\kappa\sigma_R}{G\Sigma}=
2\pi\frac{1-e^{-u^2}I_0(u^2)}{u}.
$$
右辺は $u\to0$ で $2\pi u$、$u\to\infty$ で高々 $2\pi/u$ となる。最大値は約 $3.36$ である。成長率 $\gamma>0$ の場合は $s_*=i\gamma/\kappa$ とおき、分散関係を $1-s_*^2$ で割る。右辺の級数の各正項は $\gamma$ の増加とともに減り、無限大で零へ行く。従って $\gamma=0$ の右辺が一を超える波数に指数成長解が一つある。この局所軸対称の判定に使う :term[Toomre 安定度パラメータ]{en="Toomre stability parameter" yomi="とぅーむれあんていどぱらめーた"} は $Q_*=\kappa\sigma_R/(3.36G\Sigma)$ である。非軸対称の過渡増幅や棒形不安定性については、この条件だけで結論できない。

:::example[二枚の円盤で異なる分散]
同じ $\Sigma,\kappa$ の星円盤で $\sigma_R$ を二倍にすると $Q_*$ は二倍になる。$Q_*=0.8$ から $1.6$ へ移れば局所軸対称成長の有無が切り替わる。ただし気体と星が同時にある円盤では、両成分の密度摂動が同じポテンシャルに入り、別々の $Q$ だけでは混合系を判定できない。
:::

### 二つの気体成分の局所安定条件

平らな円盤で二つの流体が同じ $\Omega,\kappa$ の円運動をする理想化を考える。面密度と音速を $\Sigma_i,c_i$ とする。半径の局所波数が $k$ の軸対称モードでは、速度を消去した線形式から $(\omega^2-\kappa^2-c_i^2k^2)\delta\Sigma_i=k^2\Sigma_i\delta\Phi$ である。Poisson 式によって $\delta\Phi=-2\pi G(\delta\Sigma_1+\delta\Sigma_2)/|k|$ である。全密度摂動が零でないモードについて
$$
1=2\pi G|k|\sum_i\frac{\Sigma_i}{\kappa^2+c_i^2k^2-\omega^2}
$$
となる。純成長率を $\omega=i\gamma$ と置くと右辺は $\gamma>0$ で狭義に減少し、無限大で零へ行く。従って中立点 $\gamma=0$ で右辺が一を超える波数に成長がある。密度摂動の和が零で重力を作らないモードでは、各成分は正の音速と $\kappa^2>0$ の復元力に従い、純成長はない。両方の成分が単独で安定でも、同じ摂動ポテンシャルを作るため和が一を超え得る。

:::example[各成分単独で安定でも混合系で不安定]
$c_1=c_2=c$、$\Sigma_1=\Sigma_2$ とすると、両方を同じ位相で動かすモードは面密度 $\Sigma_1+\Sigma_2$ の一流体に等しい。各成分の $Q_g=3/2$ なら、同位相のモードは合計面密度の $Q=3/4$ となって不安定である。逆位相で面密度摂動を相殺するモードにはこの自重力の強化がない。
:::

### 恒星の局所応答核を作用変数から導く

軸対称の波に限り、円運動のガイド中心を $X$、径作用を $J_R$、位相を $\theta_R$ とする。非摂動 Hamiltonian の径方向部分は $\kappa J_R$、径変位は $\sqrt{2J_R/\kappa}\cos\theta_R$ である。背景の尺度に比べ波長が短い局所近似で、ガイド中心の密度と $\kappa$ を一定に扱い、$J_R$ の確率密度を
$$
P(J_R)=\frac{\kappa}{\sigma_R^2}e^{-\kappa J_R/\sigma_R^2}
$$
とする。$R-X$ と径速度が調和振動子の位相空間で Gaussian となる分布である。摂動のポテンシャルを $\Phi_k e^{ikR-i\omega t}$ とする。円運動に沿って読むと径方向因子は、$A=\sqrt{2J_R/\kappa}$ を使って $e^{ikA\cos\theta_R}$ である。整数の通常 Bessel 関数 $J_n$ による展開 $e^{ikA\cos\theta_R}=\sum_n i^nJ_n(kA)e^{in\theta_R}$ を使う。係数は一周の角度積分で定義できる。

線形化した無衝突方程式は
$$
(\partial_t+\kappa\partial_{\theta_R}) f_1
=(\partial_{J_R} f_0)\partial_{\theta_R}\Phi_1
$$
である。左辺に非摂動 Hamiltonian との Poisson 括弧を入れ、右辺は摂動 Hamiltonian と $f_0(J_R)$ の括弧を移項した。従って $n$ 次成分の $f_1$ は $\partial_{J_R} f_0$ とポテンシャルの同じ成分の積に $n/(n\kappa-\omega)$ を掛けた値になる。応答を再び固定した $R$ の密度へ戻すと、位相積分で $e^{-ikA\cos\theta_R}$ と掛け合わせるため $J_n(kA)^2$ が残る。比 $s_*=\omega/\kappa$ とし、$n$ と $-n$ を合わせると
$$
\frac{\delta\Sigma_k}{\Phi_k}
=-\frac{2\Sigma}{\sigma_R^2}
\sum_{n\ge1}\frac{\langle J_n(kA)^2\rangle_{J_R}}{1-s_*^2/n^2}
$$
となる。ここでは $\partial_{J_R}f_0=-(\kappa/\sigma_R^2)f_0$ を使った。成長モードも含め $\omega$ の実部が共鳴に重なるときは、実軸の極をどう通るかに初期条件が入る。この節では $\omega$ が純虚数か零の近傍にある安定性を主に扱う。

作用平均は角度の差 $\delta$ を導入して求められる。位相を一周平均した $e^{ikA[\cos(\theta_R+\delta)-\cos\theta_R]}$ の Fourier 係数は $J_n(kA)^2$ である。さらに $P(J_R)$ を使って作用平均すると、$A\cos\theta_R,A\sin\theta_R$ がともに分散 $\sigma_R^2/\kappa^2$ の独立な Gaussian であることから、期待値は
$$
\exp[-k^2\sigma_R^2(1-\cos\delta)/\kappa^2]
=e^{-\chi}e^{\chi\cos\delta}
$$
である。右辺の $n$ 次係数が $e^{-\chi}I_n(\chi)$ なので、これが $\langle J_n(kA)^2\rangle$ である。自己重力の $\Phi_k=-2\pi G\delta\Sigma_k/|k|$ を代入すると
$$
1=\frac{4\pi G\Sigma |k|}{\kappa^2\chi}e^{-\chi}
\sum_{n\ge1}\frac{I_n(\chi)}{1-s_*^2/n^2}
$$
となり、上の $\mathcal F$ を使った式に戻る。冷たい極限 $\chi\to0$ では $I_1=\chi/2+\order{\chi^3}$、$I_{n\ge2}=\order{\chi^2}$ なので、$\omega^2=\kappa^2-2\pi G\Sigma|k|$ が戻る。有限の分散では各星が波を横切って異なる位相の力を受けるため、応答の強さが変わる。

:::example[相対的に小さい波長]
中立点で $\chi\gg1$ なら $e^{-\chi}I_0(\chi)$ は一より小さく、$\mathcal F(0,\chi)\le1/\chi$ である。従って $2\pi G\Sigma |k|\mathcal F/\kappa^2$ は高波数で高々 $|k|^{-1}$ で減る。冷たい円盤では同じ量が $|k|$ に比例して増える。ここでは速い星が有限の距離を一周で移動することが、短波長の自己重力応答の抑制へつながる。
:::

## 演習問題

::::exercise{#exr-cosmo-galaxies-1 level=basic}
一様球内部の半径三倍で重力加速度は何倍か。
:::solution
内包質量が二十七倍、面積因子が九倍なので三倍となる。
:::
::::

::::exercise{#exr-cosmo-galaxies-2 level=advanced}
剛体回転の $\kappa$ を求めよ。
:::solution
$\Omega$ の径微分が零なので :eqref[eq-cosmo-epicycle] から $\kappa=2\Omega$ である。
:::
::::

::::exercise{#exr-cosmo-galaxies-3 level=basic}
平坦回転で $v_c=200\,\mathrm{km\,s^{-1}}$、$r=10\,\mathrm{kpc}$ の質量を太陽質量で概算せよ。
:::solution
$G=4.30\times10^{-6}\,\mathrm{kpc\, (km\,s^{-1})^2}/M_\odot$ として $rv_c^2/G\simeq9.3\times10^{10}M_\odot$ である。
:::
::::

::::exercise{#exr-cosmo-galaxies-4 level=advanced}
一定異方性で $\nu\propto r^{-4}$、平坦曲線の一定分散解を求めよ。
:::solution
:eqref[eq-cosmo-jeans-spherical] の微分を入れれば $(4-2\beta)\sigma_r^2=v_c^2$ となる。従って $\sigma_r^2=v_c^2/(4-2\beta)$ である。
:::
::::

::::exercise{#exr-cosmo-galaxies-5 level=advanced}
一定速度近似で遭遇の寄与が対数になることを示せ。
:::solution
キックの二乗には $b^{-2}$、遭遇率には $b\,\dd b$ があるので積分は $\log(b_{\max}/b_{\min})$ となる。
:::
::::

::::exercise{#exr-cosmo-galaxies-6 level=basic}
他を固定して速度二倍で緩和時間は何倍か。
:::solution
キックが小さくなるうえ到達すべき分散が増え、$t_{\rm rel}\propto v^3$ から八倍になる。
:::
::::

::::exercise{#exr-cosmo-galaxies-7 level=advanced}
一定 $v_c$ の球で二半径間のポテンシャル差を求めよ。
:::solution
$\Phi'=v_c^2/r$ を積分し $v_c^2\log(r_2/r_1)$ となる。差だけに有限な基準を置ける。
:::
::::

::::exercise{#exr-cosmo-galaxies-8 level=advanced}
式 :eqref[eq-cosmo-disk-wave] で最も不安定な軸対称波数を求めよ。
:::solution
$x=\abs{k}\ge0$ で右辺の二階微分は $2c_s^2>0$ である。微分 $2c_s^2x-2\pi G\Sigma$ を零と置けば $x=\pi G\Sigma/c_s^2$ を得る。
:::
::::

::::exercise{#exr-cosmo-galaxies-9 level=advanced}
軸対称で滑らかな銀河の対称面で $v_c\propto R^s$ とする。$\Omega=v_c/R$、$\nu_z$ を鉛直角周波数として、Poisson 式から面の体積密度を求めよ。$s=-1/2$ の点質量の外部で結果を検算せよ。

:::solution
$R^{-1}\dd(v_c^2)/\dd R=2s v_c^2/R^2=2s\Omega^2$ であり、鉛直項は $\Phi_{zz}=\nu_z^2$ である。従って $4\pi G\rho=\nu_z^2+2s\Omega^2$ である。球対称の点質量の外部では $\nu_z^2=\Omega^2$ で、$s=-1/2$ なので右辺が零になり真空と一致する。従って $\nu_z>0$ だけでその場所に質量密度が正だとは結論できない。外側からの重力の曲率も鉛直振動へ入る。
:::
::::

## 確認テスト

::::quiz
球対称な系で円速度が一定の半径区間では、内包質量と局所密度にどの関係が成り立つか。

- [x] 内包質量は半径に比例する。
- [ ] 内包質量は一定である。
- [ ] 局所密度は一定である。

:::explanation
円速度一定なら $M(<r)=rv_c^2/G$ が半径に比例し、外側にも質量が必要となる。
:::
::::

::::quiz{answer="1" tol="1e-6"}
Kepler 円盤の $\kappa/\Omega$ はいくらか。
:::explanation
$s=-1/2$ なので $2(s+1)=1$ である。
:::
::::

::::quiz{answer="0.5" tol="1e-6"}
回転曲線が平坦な球で、$\nu\propto r^{-2}$ の等方定常トレーサーを考える。一定分散の Jeans 解について $\sigma_r^2/v_c^2$ を求めよ。
:::explanation
一定分散の仮定では $\dd\log\nu/\dd\log r=-2$ を Jeans 方程式に入れ、比 $1/2$ を得る。
:::
::::

::::quiz
同じ音速と $\kappa$ の薄い気体円盤で、軸対称不安定が起こりやすいのはどれか。

- [x] 面密度が大きい。
- [ ] 面密度が小さい。
- [ ] 面密度は関係しない。

:::explanation
面密度が増すと自己重力が強まり $Q_g$ が下がる。
:::
::::

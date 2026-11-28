---
title: 構造形成
description: 密度摂動の成長からハローの質量関数までを導き、銀河分布の音響尺度を計算に結び付ける。
notation:
  - symbol: '\delta'
    meaning: 物質密度の相対摂動
  - symbol: 'D(a)'
    meaning: 線形成長因子
  - symbol: 'T(k)'
    meaning: 物質の転送関数
  - symbol: '\sigma(M)'
    meaning: 質量尺度で平滑化した線形密度の標準偏差
  - symbol: '\delta_c'
    meaning: 球対称崩壊の線形閾値
---

## この章の目標

- 膨張下の不安定性と成長解を導く。
- 初期曲率から物質パワーへの変換を説明する。
- 球対称崩壊と確率分布からハローの存在量を計算する。

## Jeans 不安定性と線形成長

自然単位、符号 $(-,+,+,+)$ を用いる。物理波長が $H^{-1}$ より十分短く、流れが非相対論的なら共動位置 $\vb x$ と固有特異速度 $\vb u=a\dot{\vb x}$ を使う。一次連続式 $\dot\delta=-a^{-1}\nabla_x\cdot\vb u$ と Euler 式 $\dot{\vb u}+H\vb u=-a^{-1}\nabla_x(\Phi+c_s^2\delta)$、Poisson 式 $\nabla_x^2\Phi=4\pi Ga^2\bar\rho_m\delta$ を消去すると

$$
\ddot\delta+2H\dot\delta+\left(\frac{c_s^2 k^2}{a^2}-4\pi G\bar\rho_m\right)\delta=0 . \label{eq-cosmo-growth-fluid}
$$

:eqref[eq-cosmo-jeans-relation] の静止背景版に膨張による減衰項を加え、圧力項の波数を物理波数 $k/a$ に置き換えた式である。一成分で $c_s>0$ なら重力項と圧力項が同じ大きさになる波数は $k_J=a\sqrt{4\pi G\bar\rho_m}/c_s$ である。冷たい物質、平坦物質優勢で $H=2/(3t)$、$4\pi G\bar\rho_m=2/(3t^2)$ を入れ $\delta\propto t^p$ と置くと根は $p=2/3,-1$ である。従って成長モードは $a$、減衰モードは $a^{-3/2}$ に比例する。

:::example[重力と圧力の差]
$k=k_J/2$ では圧力は重力項の四分の一になる。
:::


:::example[線形領域の終わり]
平坦な物質宇宙で $\delta_i=0.01$ の成長モードが $\delta=1$ に達する線形予測では $a/a_i=100$ である。
:::

### 背景に寄与するが集まらない成分

物質優勢の平坦背景で、一部 $f$ の物質が考えている波数上では密度摂動をほとんど作らないとする。圧力なしの残りの成分 $\delta_c$ に対する重力源を $(1-f)\bar\rho_m\delta_c$ と置く。$x=\log a$ で $H\propto a^{-3/2}$ を使うと
$$
\partial_x^2\delta_c+\tfrac12\partial_x\delta_c-\tfrac32(1-f)\delta_c=0.
$$
$e^{px}$ を代入すると成長根は $p=(\sqrt{25-24f}-1)/4=1-3f/5+\order{f^2}$ である。この模型は $f$ が時間一定で、滑らかな成分の応答を零とみなせる尺度に限る。自由流するニュートリノのように波数によって集まり方が違う場合は、その境界で $p$ も変わる。

:::example[数パーセントの滑らかな成分]
$f=0.05$ なら $p\simeq0.9696$ である。同じ初期振幅で $a$ が $100$ 倍になると、$p=1$ での振幅に対して $100^{p-1}\simeq0.869$ 倍になる。線形パワーの比にはこの二乗が入る。
:::

### 膨張する背景から Newton 方程式を作る

固有位置 $\vb r=a\vb x$ を用い、全速度を $H\vb r+\vb u$ に分ける。流体に沿った微分は固有位置で $\partial_t|_{\vb r}+(H\vb r+\vb u)\cdot\nabla_r$ であり、共動位置では $\partial_t|_{\vb x}+(\vb u/a)\cdot\nabla_x$ になる。密度を $\rho=\bar\rho(1+\delta)$ とし、質量保存から $\dot{\bar\rho}+3H\bar\rho=0$ を引くと
$$
\dot\delta+\frac1a\nabla_x\cdot[(1+\delta)\vb u]=0 .
$$
流体要素の加速度 $\dd^2(a\vb x)/\dd t^2$ から背景の $(\ddot a/a)\vb r$ を除くと、残りは $\dot{\vb u}+H\vb u+(\vb u\cdot\nabla_x)\vb u/a$ である。従って摂動 Euler 式は $-\nabla_x\Psi/a-\nabla_x p/(a\rho)$ を右辺に持つ。圧力を $\delta p=c_s^2\bar\rho\delta$ と線形化し、遅い運動での Poisson 式を使えば :eqref[eq-cosmo-growth-fluid] を得る。両計量ポテンシャルが等しい Newton 極限を想定しているため、本節の $\Phi$ は運動を駆動する $\Psi$ に等しい。

$k_J$ は瞬間的な釣合いを示す。$k/k_J$ が時間とともに変化すれば、同じ Fourier モードが圧力優勢と重力優勢を移る。例として断熱的な非相対論単原子気体だけが平坦背景を支配し、温度が $a^{-2}$ に減るとする。$c_s\propto a^{-1}$、$\bar\rho\propto a^{-3}$ より $k_J\propto a^{1/2}$ である。$a_i$ で $k=2k_J(a_i)$ だったモードは $a=4a_i$ で釣り合う。瞬間的な閾値を越えた時刻の密度と速度を、線形成長式で連続につなぐ必要がある。

### 二つの圧力なし成分を組み合わせる

同じ Newton ポテンシャルを受ける二成分を取り、背景の質量分率を $f_1,f_2$ とする。$f_i$ が一定なら全密度 $\delta_m=f_1\delta_1+f_2\delta_2$ は通常の成長式に従う。一方、二式の差から重力が消え、相対差 $s=\delta_1-\delta_2$ は $\ddot s+2H\dot s=0$ を満たす。従って $s=s_i+a_i^2\dot s_i\int_{t_i}^t\dd t'/a^2(t')$ となる。物質期 $a\propto t^{2/3}$ で時間依存部は $a^{-1/2}$ である。初期に全密度が零でも各成分の密度が正負にずれた補償モードは残り得る。片方だけに圧力がある間は、相対式の右辺に圧力差が出て全密度の式とも再び結合する。

### 自由流が Jeans の近似と異なる点

圧力支持の議論では、力学を単一の密度と速度で閉じた。一方、無衝突の粒子は同じ空間の点に異なる速度で到達できる。局所的に $a$ を定数とみなし、自己重力を止めた短い区間で $\vb x(\eta)=\vb x_i+\vb v(\eta-\eta_i)$ とする。全ての速度ビンに同じ初期の分数密度 $\delta_i(\vb k)$ があるなら、後の密度は
$$
\delta(\vb k,\eta)=\delta_i(\vb k)
\int f_0(\vb v)e^{-i\vb k\cdot\vb v(\eta-\eta_i)}\dd^3v .
$$
各速度の位相がずれるため、平均速度が零でも密度は減り得る。各方向の速度分散が $\sigma_v^2$ の等方 Gaussian なら積分は $e^{-k^2\sigma_v^2(\eta-\eta_i)^2/2}$ である。この減少は各粒子のエネルギー損失によるものではなく、速度平均で起きる。膨張が効く場合は速度ごとの共動変位を積分に入れ、重力を戻す場合には各速度ビンを駆動するポテンシャルも同時に求める。従って自由流の境界を $k_J$ と同じ定数の音速で精密に表すことはできない。

:::example[異なる波数の位相混合]
$\sigma_v\Delta\eta=1\,\mathrm{Mpc}$ とする。$k=0.1\,\mathrm{Mpc^{-1}}$ で振幅は $e^{-0.005}\simeq0.995$ 倍、$k=2\,\mathrm{Mpc^{-1}}$ では $e^{-2}\simeq0.135$ 倍になる。重力を止めた局所模型なので、この数値は宇宙全体の転送関数の予測ではない。
:::

### 膨張する気体の Jeans 解

平坦物質模型で $\mathcal H=2/\eta$、断熱非相対論的な気体の音速を $c_s=c_i(\eta_i/\eta)^2$ とする。圧力の背景への寄与を無視できるが、小尺度の圧力勾配は残ると仮定する。成長方程式を共形時刻に直すと $\delta''+\mathcal H\delta'+(c_s^2 k^2-4\pi Ga^2\bar\rho_m)\delta=0$ である。$b=kc_i\eta_i^2$、$x=b/\eta$ とおけば
$$
\delta'=-\frac{x}{\eta}\delta_x,\qquad
\delta''=\frac{x^2}{\eta^2}\delta_{xx}
+\frac{2x}{\eta^2}\delta_x .
$$
摩擦項で右端が相殺し、$4\pi Ga^2\bar\rho_m=6/\eta^2$ なので $\delta_{xx}+(1-6/x^2)\delta=0$ となる。この解は $xj_2(x)$ と $xy_2(x)$ の線形結合である。ここで $j_2,y_2$ は球 Bessel 関数と第二種の球 Bessel 関数であり、$v=\delta/x$ と置けば $v''+2v'/x+(1-6/x^2)v=0$ を得る。早い時刻の $x\gg1$ では振動し、遅い $x\ll1$ では $xj_2\sim x^3/15$、$xy_2\sim-3/x^2$ なので $\eta^{-3}$ の減衰と $\eta^2$ の成長が表れる。同じモードでどの枝が残るかは、音の位相と速度の初期値で決まる。

:::example[圧力と重力の境界]
上の模型で圧力の係数は $b^2/\eta^4$、重力の係数は $6/\eta^2$ なので境界は $x=\sqrt6$ である。初期に $x=2\sqrt6$ なら、時刻 $\eta$ が二倍、尺度因子が四倍で釣り合う。それ以前の音響振動から、遅い成長の係数へ連続につなぐ必要がある。
:::

## 物質パワースペクトルと転送関数

:term[転送関数]{en="transfer function" yomi="てんそうかんすう"} $T(k)$ は初期条件からの波数依存の進化を記録する。断熱長波長から物質期へ入り、圧力を無視すると :eqref[eq-cosmo-einstein-scalar] と $\Phi=-3\zeta/5$ から $\delta_m=(2/5)k^2a\zeta/(\Omega_{m0}H_0^2)$ となる。ここで $\delta_m$ は共動密度であり、Newton ゲージの密度と同一視するには亜地平線近似を使う。$D$ を早い物質期で $D\simeq a$ に規格化して一般化すると

$$
P_m(k,a)=\left[\frac{2k^2 T(k)D(a)}{5\Omega_{m0}H_0^2}\right]^2\frac{2\pi^2}{k^3}\mathcal P_\zeta(k). \label{eq-cosmo-matter-power}
$$

等密度時より後に入るモードはこの規格化で $T\to1$ となる。先に入るモードでは放射が重力の背景を支え、冷たい物質の自己重力が弱い。放射のポテンシャル振動が減った後に $\ddot\delta+2H\dot\delta\simeq0$、$a\propto t^{1/2}$ から $\delta\sim A+B\log t$ が得られ、小尺度が抑制される。大きい $k$ で冷たい物質の漸近形は概ね $\log(k/k_{\rm eq})/(k/k_{\rm eq})^2$ である。係数と対数の引数、バリオンの振動、質量を持つニュートリノの自由流は連立 Boltzmann 解で決まる。遅い時期にも成長が波数依存なら上式の $D(a)$ への分離は近似である。

:::example[低波数の傾き]
$\mathcal P_\zeta\propto k^{n_s-1}$、$T=1$ なら $P_m\propto k^{n_s}$ となる。
:::


:::example[小尺度の傾き]
$n_s=1$ として上の漸近形を入れると $P_m\propto k^{-3}\log^2(k/k_{\rm eq})$ となる。
:::

### 放射と物質を含む背景での成長

放射が膨張に寄与し、その密度摂動を無視できる亜地平線の近似を取る。$y=a/a_{\rm eq}$ とすると $\Omega_m(a)=y/(1+y)$ で、$H^2\propto y^{-4}(1+y)$ である。$\dd/\dd\log a=y\,\dd/\dd y$ を :eqref[eq-cosmo-growth-fluid] の圧力なしの場合へ使う。
$$
\delta_{yy}+\frac{2+3y}{2y(1+y)}\delta_y-\frac{3}{2y(1+y)}\delta=0.
$$
この方程式を :term[Mészáros 方程式]{en="Meszaros equation" yomi="めーさーろしほうていしき"} と呼ぶ。$u=1+3y/2$ は代入すると解である。$\delta=uv$ と置き、$v'$ の一次式を積分すると $v'\propto[y\sqrt{1+y}\,u^2]^{-1}$ である。従って独立解を
$$
\delta_-(y)=\left(1+\frac32 y\right)\log\frac{\sqrt{1+y}+1}{\sqrt{1+y}-1}-3\sqrt{1+y}
$$
と選べる。微分して前の式に代入するか、$v'$ が同じ積分因子を持つことを確かめればよい。$y\ll1$ では $\delta_-=\log(4/y)-3+\order{y\log y}$ で、放射期にある定数・対数の基底を成す。$y\gg1$ では $u$ が $a$ に比例し、他方は $y^{-3/2}$ に減衰する。放射期で $aH\propto a^{-1}$ のため、地平線進入は $a_{\rm ent}\propto k^{-1}$ である。進入時と等密度時の尺度比を対数へ入れると、小波長側の $\log(k/k_{\rm eq})$ が現れる。進入時の厳密な係数には放射摂動を復元する必要がある。

:::example[等密度での成長解]
$y=1$ で成長解は $u=5/2$ である。この規格化では $y\to0$ の値は $1$ で、放射期の先頭では成長解単独が定数に近い。
:::

### 球窓が波数を選ぶ仕組み

半径 $R$ の球内で一定、外で零の窓を体積で正規化する。Fourier 変換は角度積分 $2\sin(kr)/(kr)$ の後
$$
W(kR)=\frac3{R^3}\int_0^R r^2\frac{\sin(kr)}{kr}\dd r
=\frac{3[\sin(kR)-kR\cos(kR)]}{(kR)^3}
$$
である。分散は $\sigma_R^2=(2\pi^2)^{-1}\int k^2P(k)W(kR)^2\dd k$ となる。$P=Ak^n$ の場合、変数 $x=kR$ で $\sigma_R^2\propto R^{-(3+n)}$ を得る。積分が存在するには、小さい $x$ で $W\to1$ より $n>-3$、大きい $x$ で $W=\order{x^{-2}}$ より $n<1$ が必要である。

:::example[平滑化半径を大きくする]
$n=-2$ では $R$ を二倍にすると分散は半分、標準偏差は $1/\sqrt2$ 倍である。質量は同じ背景なら $R^3$ に比例して八倍になるため、$\sigma(M)\propto M^{-1/6}$ となる。
:::

### 等密度波数と転送関数の規格化

質量変化が効く成分を無視し、$a_{\rm eq}=\Omega_{r0}/\Omega_{m0}$ とする。平坦な物質・放射背景で $k_{\rm eq}=a_{\rm eq}H(a_{\rm eq})$ と定義すれば
$$
k_{\rm eq}=H_0\sqrt{2\Omega_{m0}/a_{\rm eq}}
           =H_0\sqrt2\,\Omega_{m0}/\sqrt{\Omega_{r0}} .
$$
等密度では物質も放射も $H^2$ に同じ量を寄与するため $\sqrt2$ が付く。$k>k_{\rm eq}$ のモードは放射期に進入する。一方で、十分小さい $k$ は放射の圧力による抑制を受ける前に物質が優勢となる。物質優勢に入った時点で低波数の振幅を規格化すれば $T(k\to0)=1$ となる。遅い時刻に観測と比べるなら、同じ $D$ の規格化を維持するか、今日 $D(1)=1$ に置き換えた分をパワーの振幅に繰り込む。

:::example[波数の単位を復元する]
$H_0=70\,\mathrm{km\,s^{-1}\,Mpc^{-1}}$、$\Omega_{m0}=0.3$、$a_{\rm eq}=1/3400$ とすれば $k_{\rm eq}\simeq 0.0105\,\mathrm{Mpc^{-1}}$ である。$H_0/c$ を使うことで時間の逆数を長さの逆数へ変換した。$h\,\mathrm{Mpc^{-1}}$ で表す数値はこれを $h=0.7$ で割る。
:::

### 実空間相関と無次元パワー

統計的に一様・等方な密度場で $\langle\delta_{\vb k}\delta_{\vb k'}^*\rangle=(2\pi)^3\delta^3(\vb k-\vb k')P_m(k)$ とする。二点関数の角度積分を実行すると
$$
\xi_m(r)=\int_0^\infty \frac{k^2P_m(k)}{2\pi^2}\frac{\sin(kr)}{kr}\dd k .
$$
$\Delta_m^2=k^3P_m/(2\pi^2)$ は対数波数当たりの密度分散であり、$P_m$ は長さの三乗の次元を持つ。特に線形解が $P_m\propto k^{-3}\log^2(k/k_{\rm eq})$ なら $\Delta_m^2$ は対数で緩やかに増える。そのため次元付きパワーが小尺度で下がっても、相対密度の非線形性が小さくなるとは限らない。物理的な自由流や温度による小尺度の切断がなければ、平滑化前の点ごとの分散には注意が必要である。

### 放射の重力が作る対数成長の係数

平坦な放射宇宙に少量の冷たい物質を試験成分として加える。放射は完全流体、異方性応力は零とし、$x=k\eta/\sqrt3$ で $\Phi=3\Phi_i(\sin x-x\cos x)/x^3$ とする。物質の Euler 式 $(\eta\theta_c)'=k^2\eta\Phi$ を正則初期値から積分すると
$$
\theta_c=\frac{9\Phi_i}{\eta}(1-\sin x/x).
$$
実際、$k^2\int^\eta \eta'\Phi(\eta')\dd\eta'=3\int_0^x u\Phi(u)\dd u$ であり、$(\sin u-u\cos u)/u^2=-(\sin u/u)'$ を使った。$\theta_c$ は地平線の外で $k^2\eta\Phi_i/2$ に近づく。これは初期の $\Psi$ が一定とした正則条件と一致する。

連続式を積分して $\delta_c=-3\Phi_i/2+3(\Phi-\Phi_i)-9\Phi_i I(x)$ とおくと
$$
I(x)=\int_0^x\frac{1-\sin u/u}{u}\dd u
=\log x+\gamma_{\rm E}-\operatorname{Ci}(x)+\frac{\sin x}{x}-1 .
$$
ここで $\gamma_{\rm E}$ は Euler 定数、$\operatorname{Ci}'(x)=\cos x/x$、$\operatorname{Ci}(x)\simeq\gamma_{\rm E}+\log x$ である。積分は $\int_0^x(1-\cos u)\dd u/u$ と $(\sin u/u)'$ の和に分けると得られる。$x\gg1$ では $\Phi$ と $\operatorname{Ci}$ の振動部分は小さく、$\delta_c\simeq-9\Phi_i(\log x+\gamma_{\rm E}-1/2)$ となる。物質の自己重力を無視した式から対数の存在だけを求めると、その係数は任意の積分定数として残る。ここでは放射の駆動を保持したため、初期曲率で係数を定められた。

等密度付近では物質の自己重力を戻す必要がある。放射期の対数成長と物質期の $D\propto a$ を結ぶと、等密度までの $x_{\rm eq}\propto k/k_{\rm eq}$ により高波数で $\delta_c/\zeta$ に $\log(k/k_{\rm eq})$ が現れる。低波数の規格化は $\delta_c/\zeta\propto k^2$ であるから、転送関数には $\log(k/k_{\rm eq})/k^2$ の比が入る。ただし中間の係数を求めるには放射と物質がともに重力源となる時期を解く。

### 異なる窓で分散を比較する

実空間の球窓は半径で質量を指定しやすいが、Fourier 空間では振動する尾を持つ。代わりに Gaussian 窓 $W_G(kR)=e^{-(kR)^2/2}$ を使うと、べき型スペクトル $P_m(k)=Ak^n$ の分散は
$$
\sigma_G^2=\frac{A}{2\pi^2}\int_0^\infty k^{n+2}e^{-(kR)^2}\dd k
=\frac{A}{4\pi^2R^{n+3}}\Gamma\!\left(\frac{n+3}{2}\right),\qquad n>-3
$$
である。$u=(kR)^2$ に置き換えれば Gamma 関数が現れる。ここでは高波数は指数で切れるため、球窓の $n<1$ の制限はない。$R$ のべきが両窓で同じでも係数は違い、$R$ を同じ質量と解釈するには窓ごとに定義を決める必要がある。鋭い波数窓 $W=\mathbf 1_{k<K}$ なら $\sigma^2=A K^{3+n}/[2\pi^2(3+n)]$ である。この窓で独立な波数殻を足せることは、後の初通過計算に重要であるが、実空間では単位密度の球に対応しない。

## 成長因子

成長因子 $D$ は圧力のない線形方程式の成長解である。$\log a$ 微分を添字 $x$ とすると

$$
D_{xx}+\left(2+\dv{\log H}{\log a}\right)D_x-\tfrac32\Omega_m(a)D=0 . \label{eq-cosmo-growth-factor}
$$

平坦で物質と $\Lambda$ だけなら $\dot H=-4\pi G\bar\rho_m$ なので $\ddot H=-3H\dot H$、よって $D_-=H$ は元の時間式を満たす。$D=Hv$ として代入すると $\dd(a^2H^2\dot v)/\dd t=0$ である。$D\simeq a$ の規格化を得るよう初期積分の係数を合わせると
$$
D(a)=\frac{5\Omega_{m0}H_0^2}{2}H(a)\int_0^a\frac{\dd a'}{(a'H(a'))^3}.
$$
:term[成長率]{en="growth rate" yomi="せいちょうりつ"} $f=\dd\log D/\dd\log a$ は線形の速度発散 $\nabla\cdot\vb u=-aH f\delta$ を与える。異なる重力模型や成分で $H$ を同じにしても、源項が違えば成長は異なり得る。

:::example[物質期の検算]
$H=H_0a^{-3/2}$、$\Omega_{m0}=1$ を代入すると積分は $2a^{5/2}/(5H_0^3)$、$D=a$ となる。
:::


:::example[真空優勢な将来]
$H$ が定数へ近づくと、積分の後半は $\int^\infty a^{-3}\dd a$ 型で収束し、成長は有限の値へ近づく。
:::

### 成長率の指数近似を得る

$f=\dd\log D/\dd\log a$ とおく。:eqref[eq-cosmo-growth-factor] を $D$ で割ると $f_x+f^2+(2+H_x/H)f=3\Omega_m/2$ である。平坦で物質と一様な一定 $w$ の成分だけを含むなら $\Omega_{m,x}=3w\Omega_m(1-\Omega_m)$、$H_x/H=-(3/2)[1+w(1-\Omega_m)]$ である。$f=\Omega_m^\gamma\simeq1-\gamma(1-\Omega_m)$ と展開すると、一次の係数は
$$
(3w-5/2)\gamma-\tfrac32w=-\tfrac32.
$$
従って $\gamma=3(1-w)/(5-6w)$ となる。これは $\Omega_m\simeq1$ から得た近似であり、現在まで厳密な定数であるとは限らない。暗黒成分が集まる場合や重力法則が変わる場合は出発の成長方程式も変わる。

:::example[宇宙定数と成長率]
$w=-1$ では $\gamma=6/11$ である。$\Omega_m=0.3$ に便宜上外挿すると $f\simeq0.3^{6/11}\simeq0.518$ である。観測の精度を比較する場合は、この外挿を成長方程式の数値積分と照合する。
:::

### 成長因子の積分と赤方偏移空間での振幅

平坦 $\Lambda$CDM で $\lambda=\Omega_{\Lambda0}/\Omega_{m0}$ とすれば、$H=H_0\sqrt{\Omega_{m0}}a^{-3/2}\sqrt{1+\lambda a^3}$ である。成長の積分解の中身は $a'^{3/2}(1+\lambda a'^3)^{-3/2}/(H_0^3\Omega_{m0}^{3/2})$ である。$\lambda a^3\ll1$ で二項展開して積分し、外の $H$ も展開すると
$$
D(a)=a\left[1-\frac{2}{11}\lambda a^3+\order{(\lambda a^3)^2}\right],
\qquad f=1-\frac6{11}\lambda a^3+\order{(\lambda a^3)^2}.
$$
実際、内側の相対一次補正は $-(3/2)(5/11)\lambda a^3$、外の $\sqrt{1+\lambda a^3}$ は $+\lambda a^3/2$ である。$\Omega_m=1-\lambda a^3+\cdots$ と比較すると $f\simeq\Omega_m^{6/11}$ の先頭と一致する。この展開は $\lambda a^3$ が小さい時の検算用であり、今日の値は元の積分で求める。

線形パワーが $D^2$ と分離すると、球窓で測った標準偏差も $D$ に比例する。$\sigma_8(a)$ を半径 $8h^{-1}\,\mathrm{Mpc}$ で定義し、$\delta_m$ のパワーを $\sigma_8^2$ で割った形を固定する。そのとき銀河の速度と密度から見える振幅は $b_g\sigma_8$ と $f\sigma_8$ である。例えば重力と初期振幅を変えて $f$ と $\sigma_8$ を逆向きに動かすと、速度の振幅だけは保てる。幾何と成長を両方測るときは、測定から物理距離への変換も同時に指定する。

### 膨張が同じでも成長の源が変わる場合

物質が膨張に寄与する分率 $\Omega_m$ のうち、ある波数で重力を作る比を $F(k,a)$ とする簡単な模型を考える。圧力のない集まる成分に対し、Newton の Poisson 式に $F\bar\rho_m\delta$ を入れると、:eqref[eq-cosmo-growth-factor] の源は $3F\Omega_mD/2$ となる。従って $\dot H$ が同じでも成長は同じとは限らない。この一係数は各成分が重力源として受動的に応答するとみなす近似である。現実の質量を持つニュートリノでは過去の速度分布が先に決まり、他成分に対する比 $\delta_\nu/\delta_c$ は波数と時刻で変わる。圧力のない成分をまとめた転送関数だけで総物質の成長を表せないことがある。

膨張率の異なる模型を比較する際は、成長の初期値をそろえる必要がある。$D(1)=1$ に各模型を別々に規格化すると、現在の成長の差は規格化へ吸収される。原始曲率の振幅を同じにするときは早期の成長規格化と転送を固定し、今日の $\sigma_8$ を計算する。逆に二つの模型の今日の $\sigma_8$ をそろえて比べれば、より早い時刻の振幅や速度の違いを調べる比較になる。

### 成長方程式を積分解と照合する

物質と $\Lambda$ の平坦模型について、$s=\log a$ を時間変数にして $D,U=\partial_sD$ を同時に積分する。EdS の十分早い時点 $a_i$ で $D=U=a_i$ と置く。ここでは放射のない背景を厳密に仮定したので $a_i$ を小さくできる。現実に放射が効く時刻まで初期化を早めると、同じ初期値は使えない。比較する積分解の規格化は $D/a\to1$ でそろえる。

~~~python
from math import exp, log, sqrt

def integrate_growth(omega_m, initial_scale, steps):
    omega_lambda = 1 - omega_m
    start = log(initial_scale)
    step = -start / steps

    def rate(time, state):
        growth, slope = state
        matter = omega_m / (omega_m + omega_lambda * exp(3 * time))
        return slope, 1.5 * matter * growth - (2 - 1.5 * matter) * slope

    def shift(state, increment, scale):
        return tuple(a + scale * b for a, b in zip(state, increment))

    state = (initial_scale, initial_scale)
    for index in range(steps):
        time = start + index * step
        first = rate(time, state)
        second = rate(time + step / 2, shift(state, first, step / 2))
        third = rate(time + step / 2, shift(state, second, step / 2))
        fourth = rate(time + step, shift(state, third, step))
        state = tuple(
            value + step * (a + 2 * b + 2 * c + d) / 6
            for value, a, b, c, d in zip(state, first, second, third, fourth)
        )
    growth, slope = state
    return growth, slope / growth

def growth_quadrature(omega_m, intervals=10000):
    def integrand(a):
        return a ** 1.5 / (omega_m + (1 - omega_m) * a ** 3) ** 1.5

    odd = sum(integrand(i / intervals) for i in range(1, intervals, 2))
    even = sum(integrand(i / intervals) for i in range(2, intervals, 2))
    return 2.5 * omega_m * (4 * odd + 2 * even + integrand(1)) / (3 * intervals)

for omega_m in (1.0, 0.3):
    growth, rate = integrate_growth(omega_m, initial_scale=1e-3, steps=1000)
    print(omega_m, f"{growth:.7f}", f"{growth_quadrature(omega_m):.7f}", f"{rate:.7f}")
~~~

物質だけなら $D(1)=1$、$f(1)=1$ である。$\Omega_{m0}=0.3$ なら $D(1)\simeq0.779$、$f(1)\simeq0.513$ となる。成長率に対する $\Omega_m^{6/11}$ の外挿は約 $0.518$ なので、ここでは相対約 $1\%$ の差がある。積分解は $D$ を直接与え、微分方程式は同じ時刻の傾きも与える。物質以外の滑らかな成分が宇宙定数から逸脱する場合は、積分解の $H$ だけを新しい背景に入れ替えることは一般にできない。

## 球対称崩壊と非線形構造

均一な球状過密領域を圧力なしの独立な球殻として近似する。$\Lambda=0$ で $\ddot R=-GM/R^2$ であり、束縛解 $R=A(1-\cos\vartheta)$、$t=B(\vartheta-\sin\vartheta)$、$A^3=GMB^2$ は微分すれば運動式を満たす。背景 Einstein–de Sitter の $\bar\rho=1/(6\pi Gt^2)$ と球内密度を比較すると $1+\delta_{\rm nl}=(9/2)(\vartheta-\sin\vartheta)^2/(1-\cos\vartheta)^3$ となる。小さい $\vartheta$ で $\delta_{\rm nl}\simeq3\vartheta^2/20$ である。$t^{2/3}$ に比例する線形成長を外挿すると $\delta_{\rm lin}=(3/5)[3(\vartheta-\sin\vartheta)/4]^{2/3}$ であり、形式崩壊 $\vartheta=2\pi$ で $\delta_c=(3/5)(3\pi/2)^{2/3}\simeq1.686$ となる。

最大半径 $R_{\rm ta}=2A$ では全エネルギーがポテンシャルだけである。エネルギー保存とビリアル関係 $2K+\Omega=0$ により、同形の平衡球の $\Omega_{\rm vir}=2\Omega_{\rm ta}$、$R_{\rm vir}=R_{\rm ta}/2$ となる。背景密度を形式崩壊時刻で評価すると、転回の四分の一なので密度比 $1+\delta_{\rm vir}=\rho_{\rm vir}/\bar\rho(t_{\rm col})$ は $18\pi^2$ となる。殻の交差や非球対称性がある実際のハローでは、半径を測る密度閾値を明記する。

:::example[転回での密度]
$\vartheta=\pi$ の密度比は $9\pi^2/16\simeq5.55$ である。
:::


:::example[ビリアル化した球]
半径が転回の半分なら密度八倍、背景密度四分の一なので比は転回時の三十二倍である。
:::

### ビリアル半径を求める幾何学

EdS 背景で一様な球の殻交差を無視する。最大膨張時の半径 $R_{\rm ta}$ で運動エネルギーは零である。保存総エネルギーは $\Omega_{\rm ta}=-3GM^2/(5R_{\rm ta})$ である。最終的な一様ビリアル球で $2K+\Omega=0$ を仮定すると $E=\Omega_{\rm vir}/2$ なので $R_{\rm vir}=R_{\rm ta}/2$ となる。形式上の崩壊時刻は転回時刻の二倍で、背景密度が四分の一になる。この間に球内密度は八倍になるため、崩壊時の背景比は転回での背景比の三十二倍である。転回で $(3\pi/4)^2$ を使うと $18\pi^2$ を得る。軌道混合を扱わない球模型のため、最後の半径はビリアル条件による推定である。

:::example[密度基準を変える]
EdS の背景比 $18\pi^2\simeq178$ に対し、同じ質量を $\Delta=200$ で定義する球の半径は $(18\pi^2/200)^{1/3}\simeq0.961$ 倍である。背景密度と臨界密度が異なる模型で定義を比較する際には、分母の密度も確認する。
:::

### 非線形の密度式と二次の成長

球内の局所膨張率を $H_{\rm loc}=\dot R/R$ とする。球の質量保存から $\dot\delta=-3(1+\delta)(H_{\rm loc}-H)$ であり、$\dot H_{\rm loc}+H_{\rm loc}^2=-4\pi G\bar\rho_m(1+\delta)/3$ である。これは圧力なし、宇宙定数零で中心から同じ倍率で全ての殻が動く解である。背景の加速度式との差を取り、$H_{\rm loc}-H=-\dot\delta/[3(1+\delta)]$ を使うと
$$
\ddot\delta+2H\dot\delta-\frac43\frac{\dot\delta^2}{1+\delta}
=4\pi G\bar\rho_m\delta(1+\delta) .
$$
$\delta$ が小さいと非線形項を落とした通常の成長式になる。EdS の時間変数 $x=\log a$ で成長する線形解を $\delta_1\propto a$、二次解を $b\,\delta_1^2$ とする。左辺の線形作用素を二次解に掛けた係数は $4+1-3/2=7/2$ であり、右辺へ移した非線形項の係数は $4/3+3/2=17/6$ である。従って $\delta=\delta_1+(17/21)\delta_1^2+\cdots$ となる。過密と不足を同じ絶対値の $\delta_1$ から始めても、二次項は両方で正なので両者は対称にならない。

### 宇宙定数がある場合の球のエネルギー

質量 $M$ の圧力なし相似球で、宇宙定数は一様とする。殻の運動は $\ddot R=-GM/R^2+\Lambda R/3$ なので、単位質量のエネルギーは $\dot R^2/2-GM/R-\Lambda R^2/6$ である。球の自己重力のエネルギーは $U_G=-3GM^2/(5R)$、一様な $\Lambda$ から受ける外部エネルギーは $U_\Lambda=-\Lambda MR^2/10$ である。後者は $-\Lambda r^2/6$ を一様球全体で積分し $\langle r^2\rangle=3R^2/5$ を用いた。ビリアル条件 $2K=\sum_i\vb r_i\cdot\nabla_i U$ で、位置を一定倍率だけ変えたとき $U_G$ は次数 $-1$、$U_\Lambda$ は $2$ なので $2K+U_G-2U_\Lambda=0$ となる。

転回時に $K=0$、ビリアル時に $E=U_G/2+2U_\Lambda$ と置く。$x=R_{\rm vir}/R_{\rm ta}$、$q=U_\Lambda(R_{\rm ta})/U_G(R_{\rm ta})=\Lambda R_{\rm ta}^3/(6GM)$ とすれば
$$
1+q=\frac1{2x}+2qx^2 .
$$
$q=0$ と連続な根は $x=1/2-q/4+\order{q^2}$ である。これが EdS の半径を任意の背景へそのまま使えない理由である。また $w(a)\ne-1$ の一様成分では外部ポテンシャルが時間に依存するので、同じエネルギー保存を仮定できない。

### 粒子流の写像と殻交差

圧力のない流体要素に初期の位置ラベル $\vb q$ を付け、移動先を $\vb x=\vb q+D(t)\vb s(\vb q)$ と近似する。初期の一様な質量要素 $a^3\bar\rho_m\dd^3q$ を移すと、各流れが一対一である間は
$$
1+\delta(\vb x)=\left[\det(\delta^i_j+D\,\partial s^i/\partial q^j)\right]^{-1}.
$$
線形では $\delta=-D\nabla_q\cdot\vb s$ となる。初期が渦なしなら $\vb s$ はポテンシャル勾配であり、一次の運動方程式は時間依存を成長因子 $D$ にまとめられる。これが :term[Zeldovich 近似]{en="Zeldovich approximation" yomi="ぜるどびっちきんじ"} である。局所の変位勾配の固有値が $\lambda_i$ なら、体積比は $\prod_i(1+D\lambda_i)$ である。一つが零になると密度が発散し、殻交差が始まる。その後は複数の流れが同じ位置を占め、単一の流体速度による Euler 方程式では閉じない。

:::example[一軸だけが先に潰れる領域]
$\lambda_1=-2,\lambda_2=-0.5,\lambda_3=0$ なら $D=1/2$ で第一軸が最初に潰れる。$D=1/4$ での密度比は $[(1-1/2)(1-1/8)]^{-1}=16/7$ である。球模型のように三軸が同時に小さくなるとは限らない。
:::

### Euler 表示の二次摂動

EdS で圧力なし、渦度なし、単一流のまま成長する時期を考える。時刻変数を $s=\log a$、速度の無次元発散を $\vartheta=-\nabla\cdot\vb u/(aH)$ とする。線形では $\delta=\vartheta\propto a$ である。非線形の連続式の右辺には $-\nabla\cdot(\delta\vb u)/(aH)$ があり、Euler の発散には $-\nabla\cdot[(\vb u\cdot\nabla)\vb u]/a$ がある。$\vb u_{\vb k}=i(aH)\vb k\,\vartheta_{\vb k}/k^2$ を代入する。積の Fourier 変換を、$\vb k_1+\vb k_2=\vb k$ に対する畳み込み $\int_{\vb k_1+\vb k_2=\vb k}$ と略記すると
$$
\begin{aligned}
\partial_s\delta_{\vb k}-\vartheta_{\vb k}
 &=\int_{\vb k_1+\vb k_2=\vb k}
 \alpha(\vb k_1,\vb k_2)\vartheta_{\vb k_1}\delta_{\vb k_2},\\
\partial_s\vartheta_{\vb k}+\tfrac12\vartheta_{\vb k}-\tfrac32\delta_{\vb k}
 &=\int_{\vb k_1+\vb k_2=\vb k}
 \beta(\vb k_1,\vb k_2)\vartheta_{\vb k_1}\vartheta_{\vb k_2},\\
\alpha(\vb k_1,\vb k_2)&=\frac{(\vb k_1+\vb k_2)\cdot\vb k_1}{k_1^2},\qquad
\beta(\vb k_1,\vb k_2)=\frac{|\vb k_1+\vb k_2|^2(\vb k_1\cdot\vb k_2)}{2k_1^2 k_2^2}.
\end{aligned}
$$
$\beta$ は二つの速度を交換して平均した核である。$1/2$ の膨張係数は $\dd\log(aH)/\dd s=-1/2$ と元の速度の膨張減衰から出る。一次解を初期の $\delta_0$ と因子 $a$ の積とし、二次の密度・発散をそれぞれ $a^2 F_2$、$a^2 G_2$ の核で $\delta_0\delta_0$ に作用させる。すると $2F_2-G_2=\alpha_s$、$5G_2/2-3F_2/2=\beta$ となる。$\alpha_s=[\alpha(1,2)+\alpha(2,1)]/2$ である。連立式を解けば
$$
F_2=\frac57\alpha_s+\frac27\beta
=\frac57+\frac{\mu_{12}}2\left(\frac{k_1}{k_2}+\frac{k_2}{k_1}\right)
+\frac27\mu_{12}^2,\qquad
G_2=\frac37\alpha_s+\frac47\beta .
$$
$\mu_{12}=\hat{\vb k}_1\cdot\hat{\vb k}_2$ である。角度平均では $\mu_{12}$ が消え、$\langle\mu_{12}^2\rangle=1/3$ により $\langle F_2\rangle=17/21$ となる。これは球対称の二次係数と一致する。角度依存を残すと、二つのモードの向きが作る潮汐効果を含む。

:::example[反対向きの波数と平均密度]
$k_1=k_2$、$\mu_{12}=-1$ では $F_2=5/7-1+2/7=0$ である。従って一次に平均零の密度場へ、反対向きの二波から零波数の密度を作らない。これは質量保存の検算である。一方、等波数で垂直な二波では $F_2=5/7$ となり、その和の波数の密度を二次で作れる。
:::

初期の $\delta_0$ が Gaussian なら、線形の三点は零である。しかし二次場に含まれる二本と残り二本の一次場を Wick 収縮すると $2F_2(\vb k_1,\vb k_2)P_L(k_1)P_L(k_2)$ が生じる。他の二通りの頂点を合わせれば、線形成長の後だけでも密度の三点は非零となる。原始の三点を観測の三点から推定するには、この後の重力の寄与も予測する。

### 粒子で非線形領域を進める

殻交差後も粒子の位相空間の流れを追うには、共動位置 $\vb x$ と単位質量当たりの共役運動量 $\vb p=a^2\dot{\vb x}$ を用いる。Hamiltonian は $\vb p^2/(2a^2)+\Psi(\vb x,t)$ なので $\dot{\vb x}=\vb p/a^2$、$\dot{\vb p}=-\nabla_x\Psi$ となる。時刻が依存する系でも、粒子の移動と運動量更新をそれぞれこの式で積分できる。例えば前半の運動量更新、位置移動、後半の運動量更新からなる二次精度の対称積分では、各更新を $\int\dd t$ または $\int\dd t/a^2$ で重み付けする。重力はその段階の粒子密度から Poisson 式で求める。

周期箱では平均密度を引いた Fourier 成分 $\delta_{\vb k}$ に $-4\pi Ga^2\bar\rho_m/k^2$ を掛けてポテンシャルを求める。零波数は膨張背景に含め、逆 Laplace を適用しない。計算箱の長さ $L$ により最低波数は $2\pi/L$、格子間隔 $\Delta x$ により一方向の Nyquist 波数は $\pi/\Delta x$ である。粒子数、力の解像度、箱の大きさ、時間刻みを別々に変え、対象とする質量や波数の統計量が収束するか確認する。線形の成長と初期の位相、球崩壊の閾値のどれに一致したかを記録すれば、数値誤差の出所を切り分けられる。

### 転回の存在と数値的な崩壊時刻

正の宇宙定数を持つ球対称の殻について $F(R)=2E+2GM/R+\Lambda R^2/3$ と置く。速度は $\dot R=\pm\sqrt{F(R)}$ である。小半径と大半径で $F>0$ になり、正の半径の極小点は $R_s=(3GM/\Lambda)^{1/3}$ である。$F'(R)=-2GM/R^2+2\Lambda R/3$ よりそこでの値は $2E+3GM/R_s$ である。初期に小半径を膨張する殻が有限時間に転回するには、この極小値が負である必要がある。負なら二根の小さい方が転回半径であり、正なら途中で速度は零にならない。極小値が零なら二重根へ近づく時間は対数で発散する。この三つの違いは $F$ の根の多重度から決まる。

小さい方の根を $R_{\rm ta}$ とすると、共通の Big Bang 時刻を選んだ膨張枝の転回時刻は $t_{\rm ta}=\int_0^{R_{\rm ta}}\dd R/\sqrt{F(R)}$ である。$F$ が時間に依存しないので落下にも同じ時間がかかる。単純根の積分端では $F(R)\simeq |F'(R_{\rm ta})|(R_{\rm ta}-R)$ だから $R=R_{\rm ta}(1-u^2)$ で平方根の特異性を取り除ける。崩壊時刻を固定したい場合、初期の小振幅 $\delta_i$ と成長解の速度をそろえて $E$ を計算し、$2t_{\rm ta}$ が目標になるまで $\delta_i$ を調整する。その初期振幅を線形成長で外挿したものが模型依存の $\delta_c$ である。転回できる最小の初期振幅と、指定時刻までに崩壊する振幅は同じではない。

## Press–Schechter 理論と質量関数

質量 $M=(4\pi/3)\bar\rho_{m0}R^3$ の共動球で線形密度を平均する。規格化した窓の Fourier 変換は $W(x)=3(\sin x-x\cos x)/x^3$、$x=kR$ である。これは球の角度平均 $\sin(kr)/(kr)$ を径方向に積分して得る。分散は $\sigma^2(M)=\int k^2P_m(k)W^2(kR)\dd k/(2\pi^2)$ である。Gaussian な密度が閾値 $\delta_c$ を超える確率は $\frac12\operatorname{erfc}(\delta_c/(\sqrt2\sigma))$ である。:term[Press–Schechter 理論]{en="Press–Schechter theory" yomi="ぷれすしぇくたーりろん"} はこれを二倍し、ハローへの質量分率とする。鋭い $k$ 窓で独立な Gaussian 波数殻を加える拡散過程では、定数障壁の反射原理により初通過確率が終点越え確率の二倍になる。実空間球窓では歩みの相関があるため、二倍をそのまま厳密な初通過とは解釈できない。質量分率を微分すると

$$
\frac{\dd n}{\dd M}=\sqrt{\frac2\pi}\frac{\bar\rho_{m0}}{M^2}\frac{\delta_c}{\sigma}\left|\dv{\log\sigma}{\log M}\right|e^{-\delta_c^2/(2\sigma^2)} . \label{eq-cosmo-press-schechter}
$$

実際の高精度質量関数は N 体計算とハロー定義に依存する。相関した窓や崩壊障壁が変化する場合の初通過問題は :cite[peebles-cosmology] の構造形成の枠組みを補う確率過程を要する。

:::example[希少ハロー]
$\delta_c/\sigma$ が $2$ から $4$ になると指数因子の比は $e^{-6}$ である。
:::


:::example[平滑化の尺度]
球の半径を二倍にすると対応付ける質量は八倍になる。
:::

### 初通過で質量分率を得る

Gaussian 場を波数の鋭い上限で平滑化し、新しく加わる独立な Fourier 殻の分散を時間変数 $S$ とする。密度の変化は増分分散 $\dd S$ の独立な Gaussian なので、障壁 $\delta_c$ に触れていない確率密度 $p(\delta,S)$ は $\partial_S p=\tfrac12\partial_\delta^2p$、$p(\delta_c,S)=0$ に従う。初期に $\delta=0$ とし、熱核と鏡像を差し引くと
$$
p(\delta,S)=\frac{e^{-\delta^2/(2S)}-e^{-(2\delta_c-\delta)^2/(2S)}}{\sqrt{2\pi S}},
\qquad \delta<\delta_c.
$$
各項が拡散式を満たし、差が障壁で零となる。初期の鏡像は定義域の外にあるので初期密度は原点に集中する。生存確率を積分すると $\operatorname{erf}(\delta_c/\sqrt{2S})$ である。微分の負は $\delta_c e^{-\delta_c^2/(2S)}/(\sqrt{2\pi}S^{3/2})$ となり、障壁に初めて達する密度を与える。実空間の球窓では増分が相関するため、この像の解をそのまま使えない。

:::example[初通過と一点での超過の差]
$S=\delta_c^2$ までに初通過した分率は $\operatorname{erfc}(1/\sqrt2)\simeq0.3173$ である。一点の Gaussian 分布がその時点で閾値を越える確率は約 $0.1587$ であり、過去に一度越えて戻る道を数えない分だけ小さい。
:::

### 長い密度モードによるハロー数の応答

ハロー形成の小尺度分散 $\sigma^2$ のほかに、調べる領域でほぼ一定の線形な長波長密度 $\delta_L$ があるとする。小尺度摂動が必要とする閾値は $\delta_c-\delta_L$ と置く。分散は変化しない近似では、:eqref[eq-cosmo-press-schechter] で閾値に依存する部分は $(\delta_c-\delta_L)\exp[-(\delta_c-\delta_L)^2/(2\sigma^2)]$ である。$\delta_L=0$ で対数微分すると
$$
b_L=\left.\pdv{\log(\dd n/\dd M)}{\delta_L}\right|_0
=\frac{\nu^2-1}{\delta_c},\qquad \nu=\frac{\delta_c}{\sigma}.
$$
$b_L$ は初期共動領域の単位体積で数えた応答である。物質の流れに沿いハロー数を保存すると、最終の単位体積での密度には因子 $1+\delta_L$ が加わるため、線形では $b_E=1+b_L$ である。これは一定閾値と Gaussian な Press–Schechter 模型の応答であり、現実のハローの選択や形成史への依存を含めるには補正が要る。

:::example[希少なハローほど長波長へ強く応答する]
$\delta_c=1.686$、$\nu=3$ なら $b_L=8/1.686\simeq4.74$、$b_E\simeq5.74$ となる。$\nu=1$ では $b_L=0$ でも $b_E=1$ である。初期体積と観測時の体積で数えた応答を比べる際は、どちらを用いたか明示する。
:::

### バリアを越える質量と質量関数の換算

$S=\sigma^2$ とし、固定した時刻での障壁の高さを $\delta_c$ とする。単調に減少する $\sigma(M)$ について $F(>M)=\operatorname{erfc}[\delta_c/(\sqrt2\sigma(M))]$ であり、区間 $(M,M+\dd M)$ の質量分率は $-\partial_MF\,\dd M$ である。一個のハローが質量 $M$ を担うので、数密度に直すには $\bar\rho_{m0}/M$ を掛ける。$\partial_M(\delta_c/\sigma)=-(\delta_c/\sigma)\partial_M\log\sigma$ と $\operatorname{erfc}'(x)=-2e^{-x^2}/\sqrt\pi$ を用いれば :eqref[eq-cosmo-press-schechter] を再現する。$\nu=\delta_c/\sigma$ を変数にすれば質量分率は $\sqrt{2/\pi}e^{-\nu^2/2}\dd\nu$ である。$\nu$ が $0$ から $\infty$ に及べば積分は $1$ になる。小尺度のパワーに切断があり $\sigma$ が有限で止まる場合、$\nu$ の下限が正であり、同じ模型で数えたハロー質量分率は $1$ 未満になる。

二つの滑らかさ $S_1<S_2$ で同じ鋭い波数窓を使い、$\delta_1<\delta_c$ から初通過を始める。独立増分のため、それ以前に障壁へ触れていないという条件の下では、$S_2-S_1$ の過程に障壁高さ $\delta_c-\delta_1$ を入れればよい。像法の熱核で $\delta_c$ と $S$ をこの二量に置き換えると条件付き初通過密度が得られる。この性質が長波長密度に応じた障壁の移動を正当化する。実空間窓では、一般に長い側の履歴が短い側の増分にも関係する。

### 線形に動く障壁での初通過

鋭い波数窓の Gaussian 過程 $\delta(S)$ と、障壁 $\delta_c+\beta S$ を考える。$X(S)=\delta(S)-\beta S$ は一定障壁 $\delta_c$ に対し、平均増分 $-\beta\dd S$ と分散増分 $\dd S$ を持つ。吸収前の確率密度 $p(X,S)$ は $\partial_S p=\beta\partial_Xp+\partial_X^2p/2$ に従う。$p=e^{-\beta X-\beta^2S/2}p_0$ を代入すれば、$p_0$ は通常の拡散式と同じゼロ境界条件、原点のデルタ初期値を持つ。従って $p_0$ に定数障壁の鏡像の式を使える。吸収の確率流束は $-\partial_X p/2$ を障壁で評価する。そこで $p_0=0$ なので、定数障壁の初通過密度に $e^{-\beta\delta_c-\beta^2S/2}$ を掛けた
$$
f(S)=\frac{\delta_c}{\sqrt{2\pi}S^{3/2}}
       \exp[-(\delta_c+\beta S)^2/(2S)]
$$
となる。前の係数は $\delta_c+\beta S$ ではなく $\delta_c$ である。実空間の球窓や非線形の崩壊条件では独立増分や線形障壁の仮定を改める必要があるが、変動する閾値が質量関数を変える仕組みをこの解で具体的に見られる。

### 質量で平均したハローバイアス

定数障壁模型で全質量がハローに属するとし、$\nu$ の範囲を $[0,\infty)$ とする。質量で重み付けした Euler バイアスの平均は
$$
\int_0^\infty\sqrt{2/\pi}e^{-\nu^2/2}
 \left[1+\frac{\nu^2-1}{\delta_c}\right]\dd\nu=1 .
$$
半正規分布の二乗平均が $1$ なので第二項は零である。全ハローの質量を足せば物質自身になることと整合する。数で重み付けする場合は $\bar\rho_m/M$ の換算が残るので、一般には同じ平均にならない。観測銀河の選択ではハロー質量ごとの銀河数や明るさが入るため、観測バイアスをこの質量平均に置き換えない。

### べき型のスペクトルで質量関数の形を求める

線形の $P_m=Ak^n$ を球窓で平滑化し、$-3<n<1$ とする。質量 $M\propto R^3$ より $\sigma(M)\propto M^{-s}$、$s=(n+3)/6>0$ である。$\sigma(M_*)=\delta_c$ と定義すれば $\nu=(M/M_*)^s$ となる。質量関数に代入すると
$$
\frac{\dd n}{\dd M}=\sqrt{\frac2\pi}\frac{\bar\rho_{m0}}{M^2}
s\left(\frac M{M_*}\right)^s
\exp[-(M/M_*)^{2s}/2].
$$
$M\ll M_*$ では微分数密度のべきは $-2+s$ であり、$M\gg M_*$ では指数抑制が効く。初期のパワーが同じなら $\sigma$ は $D$ に比例し、$M_*(a)\propto D(a)^{1/s}$ である。従って同じ平均物質密度の宇宙でも、成長すると特徴質量が大きい側へ動く。小質量で数密度が発散しても、質量分率の積分は $s>0$ により収束する。現実のパワーでは傾きが一定でなく、自由流などで切断があればこの小質量の外挿も変わる。

:::example[質量と成長の非線形な換算]
$n=-2$ なら $s=1/6$ である。成長因子が二倍になると定数障壁を越える特徴質量は $2^6=64$ 倍になる。これは単純なべき型の分散と障壁固定の算例であり、各ハローの質量が同じ倍率で増える意味ではない。
:::

### 障壁が上がり続けるときの未崩壊分率

前の線形障壁で $\beta>0$ とする。定数障壁に対し平均変位が $-\beta\dd S$ の過程が、位置 $X$ からいつか障壁 $\delta_c$ に達する確率を $u(X)$ とする。短い $\dd S$ だけ進んだ後の到達確率を平均すると、一次までで $u(X)=u(X)-\beta u'(X)\dd S+u''(X)\dd S/2$ である。従って $u''/2-\beta u'=0$ となる。障壁上で $u=1$、遠く下からは $u\to0$ なので $u(X)=e^{2\beta(X-\delta_c)}$ である。初期 $X=0$ なら総初通過分率は $e^{-2\beta\delta_c}$ であり、$S$ を無限に増やしても初通過しない部分が残る。これは最初の分布を二倍して常に全質量に規格化する処方が、動く障壁には無条件で成り立たないことを示す。この計算は分散を任意に大きくできる鋭い波数窓の理想化である。

## バリオン音響振動と銀河分布

光子とバリオンの音波の波面は、バリオンが光子の抵抗から離れる時刻 $\eta_d$ までに $r_d=\int_0^{\eta_d}c_s\dd\eta$ を進む。後の銀河二点相関にはこの尺度の近くに弱い過剰が残る。これを :term[バリオン音響振動]{en="baryon acoustic oscillations" yomi="ばりおんおんきょうしんどう"} と呼ぶ。小角と狭い赤方偏移幅で共動横距離差は $D_M\Delta\theta$、径方向差は $c\Delta z/H$ なので、角度は $r_d/D_M$、赤方偏移幅は $r_dH/c$ を測る。$D_M=(1+z)D_A$ である。観測される銀河の密度は物質へのバイアスを持ち、線形で $\delta_g=b_g\delta$ とする。速度による座標ずれ $s_\parallel=x_\parallel+u_\parallel/(aH)$ を線形化し個数を保存すると $\delta_g^{\,s}=(b_g+f\mu^2)\delta$ になる。$\mu$ は視線と波数の余弦である。

:::example[横方向の尺度]
$r_d=150\,\mathrm{Mpc}$、$D_M=1500\,\mathrm{Mpc}$ なら角度は $0.1\,\mathrm{rad}$ である。
:::


:::example[速度効果]
$b_g=2,f=1$ で径方向と横方向の線形パワー比は $(3/2)^2=2.25$ である。
:::

### 薄い球殻と振動するスペクトル

音波が距離 $r_d$ を伝わってできる薄い共動球殻を $\delta_{\rm D}(r-r_d)/(4\pi r_d^2)$ と理想化する。観測中心からの角度平均の Fourier 因子は $\frac12\int_{-1}^1e^{ikr\mu}\dd\mu=\sin(kr)/(kr)$ だから、球殻の変換は $\sin(kr_d)/(kr_d)$ である。座標空間で標準偏差 $\Sigma$ の等方三次元 Gaussian と畳み込むと、Fourier で $e^{-k^2\Sigma^2/2}$ が掛かる。この因子は標準偏差一つの変位模型であり、現実の非線形成長の全補正ではない。

:::example[波数で見た間隔と減衰]
$r_d=150\,\mathrm{Mpc}$ なら最初の正の零点は $k=\pi/r_d\simeq0.0209\,\mathrm{Mpc^{-1}}$ である。$\Sigma=5\,\mathrm{Mpc}$ とすると $k=0.2\,\mathrm{Mpc^{-1}}$ の振幅を $e^{-1/2}\simeq0.607$ 倍にする。半径と幅は位相と振幅に異なる効果を持つ。
:::

### 抵抗から離れる時刻と角度・径の尺度

バリオン Euler 式の散乱項の率は $\dot\kappa_T/R_b$ である。光子の温度場が外部から与えられる近似で、バリオン速度が初期値を忘れる因子は $\exp[-\int(\dot\kappa_T/R_b)\dd\eta]$ である。従って後の時刻までの抵抗の光学的厚さ $\tau_d(\eta)=\int_\eta^{\eta_0}\dot\kappa_T/R_b\,\dd\eta'$ が $1$ 前後になる時刻を $\eta_d$ の目安に取る。光子の最終散乱の深さは重み $1/R_b$ を持たず、同じ時刻である必要はない。BAO の $r_d$ と CMB の $r_s(\eta_*)$ はそれぞれの時刻を用いて計算する。

観測で試行した距離と真の距離が違えば、同じ角度と赤方偏移から推定した横・径の共動距離は違う。基準の量を添字 $\mathrm{fid}$ で表す。横方向の音響尺度に対応する伸び率は $\alpha_\perp=(D_M/r_d)/(D_M/r_d)_{\rm fid}$、径方向は $\alpha_\parallel=[c/(Hr_d)]/[c/(Hr_d)]_{\rm fid}$ である。等方な一変数に圧縮するとき、共動体積要素は $D_M^2(c/H)\dd\Omega\dd z$ を持つから $\alpha_\perp^2\alpha_\parallel$ の三乗根が自然な体積平均となる。実際の角度依存を捨てれば、二つの伸びを個別には決められない。

### 赤方偏移空間の多重極

遠方観測者近似で各銀河への視線方向を固定する。座標移動の Jacobian は $1+\partial_\parallel u_\parallel/(aH)$ なので、個数保存は $\delta_g^s=\delta_g-\partial_\parallel u_\parallel/(aH)$ を与える。渦のない線形速度は $\vb u=i aHf\,\vb k\,\delta/k^2$ であり、$\partial_\parallel u_\parallel/(aH)=-f\mu^2\delta$ である。従って $P_g^s=(b_g+f\mu^2)^2P_m$ となる。この線形の異方性は :term[Kaiser 効果]{en="Kaiser effect" yomi="かいざーこうか"} と呼ばれる。$\mu^2=1/3+(2/3)P_2(\mu)$、$\mu^4=1/5+(4/7)P_2(\mu)+(8/35)P_4(\mu)$ を代入すると
$$
P_0^s=\left(b_g^2+\frac23 b_gf+\frac15 f^2\right)P_m,\quad
P_2^s=\left(\frac43 b_gf+\frac47f^2\right)P_m,\quad
P_4^s=\frac8{35}f^2P_m .
$$
線形の角度依存のためこれより上の多重極は零である。銀河団内のランダム速度、選択効果、広い観測角ではこの性質が崩れ得る。長波長での速度バイアスを無視し、銀河密度が線形かつ決定論的に物質へ比例する場合の結果である。

### 音響ピークの幅に効く相対変位

初期に二点が共動距離 $\vb r$ だけ離れ、それぞれの重力変位が $\vb s_1,\vb s_2$ なら、後の分離は $\vb r+\Delta\vb s$、$\Delta\vb s=\vb s_2-\vb s_1$ である。位置相関の狭いピークを変位でならすとき、Fourier 因子は $\langle e^{i\vb k\cdot\Delta\vb s}\rangle$ となる。等方で平均零、共分散 $A_{ij}=\langle\Delta s_i\Delta s_j\rangle$ の Gaussian なら、この因子は $e^{-k_i k_j A_{ij}/2}$ である。二点を同じだけ移した一様な変位は $\Delta\vb s$ に入らない。大きな流れの絶対速度を二乗して幅へ加えると、その共通変位まで数えてしまう。

線形の初期変位を $\vb s_{\vb q}=i\vb q\,\delta_{\vb q}/q^2$ とすれば
$$
A_{ij}(\vb r)=2\int\frac{\dd^3q}{(2\pi)^3}
\frac{q_i q_j}{q^4}P_m(q)[1-\cos(\vb q\cdot\vb r)]
$$
である。大きな波長では角括弧に $(q r)^2$ の抑制が掛かる。実際の減衰は向きと距離で変わるので、前の一定幅の球殻ではこの行列の主要な効果だけを一つの $\Sigma$ にまとめた。銀河分布から長波長の変位を推定して位置を戻す方法は、推定が正しければ音響ピークの幅を小さくできる。銀河バイアスや雑音を含む推定では、補正後のピークの誤差も同時に評価する。

### 銀河相関の統計誤差と BAO の距離

大きな体積 $V$ の周期箱で実 Gaussian な線形密度を考え、ある狭い波数殻に正負の波数を合わせて $N$ 個のモードがあるとする。一方の独立な複素係数の二乗は実・虚二つの Gaussian の二乗の和で、分散は平均二乗に等しい。実場の共役を独立と数えられないので、殻の平均パワーの分散は $2P^2/N$ である。銀河を独立に Poisson 抽出する近似では点の自己対により $1/\bar n_g$ がパワーへ足される。平均からこのショット雑音を引いても、推定量の先頭分散には $2(P+1/\bar n_g)^2/N$ が入る。Gaussian でないハローの分布や窓関数を持つ観測では、追加の共分散を入れる必要がある。

一定共動体積の例で、殻の幅 $\Delta k$ が $k$ より十分小さければ $N\simeq V4\pi k^2\Delta k/(2\pi)^3$ である。$V$ を二倍にすると、同じ波数範囲の標準偏差は理想条件で $1/\sqrt2$ 倍になる。銀河数密度を二倍にしても $P\bar n_g\gg1$ の範囲ではすでに宇宙分散が支配するため、同じ改善は得られない。BAO の縦横距離を測る場合は、音響のある波数と視線角の範囲を保ちつつ、観測窓が混ぜるモードを扱う。

### 音響尺度を相関とパワーの両方で表す

実空間の二点関数 $\xi_g(r)$ は、平均からの超過確率 $\dd P=\bar n_g^2[1+\xi_g(r)]\dd V_1\dd V_2$ で定義できる。等方で境界のない観測なら、狭い距離の殻にある組の数は、殻の体積 $4\pi r^2\dd r$ と $1+\xi_g(r)$ に比例する。従って音響の弱い肩を測るには、幾何が作る $r^2$ の増大から背景を分ける必要がある。実際の観測では不均一な選択関数と空のマスクがあるので、同じ幾何・選択で独立な点のカタログを作り、そのペアと比較する。線形の銀河バイアスが一定なら $\xi_g=b_g^2\xi_m$ であり、振幅だけが変わる。小尺度の形成史が波数依存する場合は、それが相関の広い形を変えるため、ピークの位置だけを固定倍の模型で合わせない。

短い距離で出発した密度の中心と、半径 $r_d$ に広がった薄い球殻を線形の同じ初期振幅が作るとする。転送の模型を $T(k)=T_{\rm sm}(k)+\epsilon_b j_0(kr_d)$ とすれば、パワーは
$$
|T(k)|^2=|T_{\rm sm}(k)|^2+
2\epsilon_b\operatorname{Re}T_{\rm sm}(k)j_0(kr_d)
+\epsilon_b^2 j_0(kr_d)^2
$$
である。バリオンの効果が小さい場合、振動の主項は滑らかな成分との交差である。殻の寄与だけを二乗して足すと、この一次の振動を落とす。実際のプラズマでは音速が時間で変わり、拡散、重力、ニュートリノの位相応答も入るため、薄い一殻で転送の全てを与えられない。この模型は周期と干渉の符号を理解するために使う。

:::example[パワーでの波数間隔]
$r_d=150\,\mathrm{Mpc}$ の薄い殻では $j_0$ の隣り合う零点の間隔は $\pi/r_d\simeq0.0209\,\mathrm{Mpc^{-1}}$ である。$j_0$ の符号は零点ごとに反転するため、同符号の山どうしの代表間隔はその約二倍である。パワーの交差項が主なら零点と振幅の符号を保つが、二乗項だけにすると符号が消え、位相の読み方が変わる。
:::

## 演習問題

::::exercise{#exr-cosmo-structure-formation-1 level=basic}
平坦物質宇宙で尺度因子三倍の成長モード比はいくらか。
:::solution
成長解は $D\propto a$ なので三倍になる。
:::
::::

::::exercise{#exr-cosmo-structure-formation-2 level=advanced}
物質宇宙の成長の時間指数を二つ求めよ。
:::solution
$\delta=t^p$ で $p(p-1)+(4/3)p-2/3=0$ となり、$(p-2/3)(p+1)=0$ から $2/3$ と $-1$ である。
:::
::::

::::exercise{#exr-cosmo-structure-formation-3 level=basic}
$\mathcal P_\zeta$ 一定、$T=1$ で波数三倍の $P_m$ は何倍か。
:::solution
:eqref[eq-cosmo-matter-power] の $k^4/k^3=k$ から三倍である。
:::
::::

::::exercise{#exr-cosmo-structure-formation-4 level=advanced}
放射背景で自己重力を無視した物質の常微分方程式を解け。
:::solution
$H=1/(2t)$ により $(t\dot\delta)'=0$ となるので、解は $A+B\log t$ である。
:::
::::

::::exercise{#exr-cosmo-structure-formation-5 level=basic}
球の転回時刻と崩壊時刻の比を求めよ。
:::solution
$t=B(\vartheta-\sin\vartheta)$ に $\pi$ と $2\pi$ を入れ、$1/2$ を得る。
:::
::::

::::exercise{#exr-cosmo-structure-formation-6 level=advanced}
球窓 $W(y)$ の小さい $y$ の展開を二次まで求めよ。
:::solution
$\sin y-y\cos y=y^3/3-y^5/30+\order{y^7}$ なので $W=1-y^2/10+\order{y^4}$ となる。
:::
::::

::::exercise{#exr-cosmo-structure-formation-7 level=basic}
$\delta_c/\sigma=3$ での質量関数の指数因子を求めよ。
:::solution
指数は $-9/2$、因子は $e^{-9/2}$ である。
:::
::::

::::exercise{#exr-cosmo-structure-formation-8 level=advanced}
$H=100\,\mathrm{km\,s^{-1}\,Mpc^{-1}}$、$r_d=150\,\mathrm{Mpc}$ の径方向 BAO の赤方偏移幅を求めよ。
:::solution
$c\simeq3\times10^5\,\mathrm{km\,s^{-1}}$ として $\Delta z=Hr_d/c=0.05$ となる。
:::
::::

::::exercise{#exr-cosmo-structure-formation-9 level=advanced}
二次の重力成長の $F_2$ を、等しい波数の二本が $60^\circ$ で交わる場合に求めよ。正反対の場合と比較し、零波数への意味を述べよ。
:::solution
$\mu_{12}=1/2$、比 $k_1/k_2=1$ なので $F_2=5/7+1/2+(2/7)(1/4)=9/7$ である。正反対では $\mu_{12}=-1$ から $5/7-1+2/7=0$ になる。正反対かつ等波数なら合成波数は零であり、平均の密度コントラストをこの対が作らないという質量保存を確認できる。
:::
::::

::::exercise{#exr-cosmo-structure-formation-10 level=advanced}
線形赤方偏移空間で $b_g=2,f=1$ とする。$P_2^s/P_0^s$ と $P_4^s/P_0^s$ を求めよ。
:::solution
$P_0^s/P_m=4+4/3+1/5=83/15$、$P_2^s/P_m=8/3+4/7=68/21$、$P_4^s/P_m=8/35$ である。従って $P_2^s/P_0^s=340/581\simeq0.585$、$P_4^s/P_0^s=24/581\simeq0.0413$ となる。径と横の二方向だけのパワー比と異なり、多重極は全方向を重み付けして平均した量である。
:::
::::

::::exercise{#exr-cosmo-structure-formation-11 level=advanced}
EdS で初期変位の主軸が二方向のみにあり、変位勾配の固有値が $-\lambda,-\lambda,0$、$\lambda>0$ とする。Zeldovich 近似の密度を $D\lambda$ の二次まで展開し、殻交差に達する成長因子も求めよ。
:::solution
体積比は $(1-D\lambda)^2$ なので密度比は $(1-D\lambda)^{-2}$ である。$D\lambda\ll1$ で $\delta=2D\lambda+3(D\lambda)^2+\order{(D\lambda)^3}$ となる。線形の $\delta_1=2D\lambda$ で書けば二次は $3\delta_1^2/4$ である。第一・第二の主軸が同時に零になる $D=1/\lambda$ が殻交差である。この変位の幾何と近似は球対称の力学から得た $17/21$ の場合と異なるため、二次係数を混用できない。
:::
::::

## 確認テスト

::::quiz
圧力を無視した冷たい物質の成長式に入る減衰項はどれか。

- [x] 膨張減衰
- [ ] 冷たい成分の輻射圧
- [ ] 対流の Schwarzschild 判定条件

:::explanation
圧力を持つ場合は $c_s^2k^2/a^2$ が重力と対抗するが、冷たい物質では消え、膨張減衰 $2H\dot\delta$ が残る。
:::
::::

::::quiz{answer="1.68647" tol="1e-4"}
球対称 EdS 崩壊の線形閾値を求めよ。
:::explanation
摂動の線形成長を $\vartheta=2\pi$ まで外挿し、$(3/5)(3\pi/2)^{2/3}$ を得る。
:::
::::

::::quiz{answer="2.25" tol="1e-6"}
$b_g=1,f=0.5$ で径方向と横方向のパワー比はいくらか。
:::explanation
$(b_g+f)^2/b_g^2=2.25$ である。
:::
::::

::::quiz
BAO が測る径方向の比を選べ。

- [x] $Hr_d/c$
- [ ] $r_d/D_M$
- [ ] $D_MH/r_d$

:::explanation
共動径差が $c\Delta z/H$ なので音響長に対応する $\Delta z$ は $Hr_d/c$ である。
:::
::::

---
title: Brown 運動と Langevin 方程式
description: 熱浴中の粒子の運動を解き、摩擦、熱雑音、拡散と移動度の関係を導く。
notation:
  - symbol: 'm'
    meaning: 粒子の質量
  - symbol: '\gamma'
    meaning: 速度に比例する摩擦係数
  - symbol: '\tau=m/\gamma'
    meaning: 速度の緩和時間
  - symbol: 'D'
    meaning: 一方向の拡散係数
  - symbol: '\mu'
    meaning: 静的移動度
---

## この章の目標

- 粒子流のつり合いから Einstein の関係を導く。
- Langevin 方程式を解いて速度相関と変位分散を求める。
- 観測時間に応じて慣性の効果を判定する。
- 定常外力と振動外力に対する移動度を求める。

## Brown 運動の現象論と Einstein の関係

### 観測と粗視化

液体中の微粒子を同じ位置から何度も放すと、終点は試行ごとに異なる。衝突の詳細を捨て、終点の確率密度または粒子の数密度 $n$ を予測するのが :term[Brown 運動]{en="Brownian motion" yomi="ぶらうんうんどう"} の記述である。以下では一方向を扱う。座標を $x$、時刻を $t$ とし、温度 $T$ は空間的に一様とする。

拡散流を $-D\partial_x n$、弱い力 $F$ による流れを $\mu F n$ とおく。粒子一個に対する :term[移動度]{en="mobility" yomi="いどうど"} $\mu$ は、定常平均速度の力に対する一次の係数である。:term[拡散係数]{en="diffusion coefficient" yomi="かくさんけいすう"} $D$ は長時間の変位分散を $2t$ で割った極限である。粒子の蓄積を $\partial_t n=-\partial_x J$ と表し、両者を加えて $J=\mu Fn-D\partial_x n$ を得る。係数が一定で $F=0$ なら $\partial_t n=D\partial_x^2 n$ である。

### 平衡と輸送係数

:::proposition[Einstein の関係]{#prop-neq-einstein}
希薄な非相互作用粒子を一様な熱浴に入れ、恒等的には零でない保存力 $F=-U'$ の下で位置によらない係数 $D,\mu$ による局所的な流束則が成り立つとき、無流束の平衡条件は $D=\mu k_{\mathrm B}T$ を与える。
:::

:::proof
:ref[def-stm-canonical] により、平衡密度は $n_{\mathrm{eq}}\propto e^{-\beta U}$、$\beta=(k_{\mathrm B}T)^{-1}$ である。このため $\partial_x n_{\mathrm{eq}}=\beta F n_{\mathrm{eq}}$ であり、$J=(\mu-\beta D)Fn_{\mathrm{eq}}$ となる。力が零でない領域で $J=0$ を課すと $\mu=\beta D$ である。係数を一様と仮定したので、その関係は全領域で成り立つ。
:::

:::example[重力による沈降]
浮力を含めた質量を $m_{\mathrm b}>0$、上向き座標を $x$ とすれば $F=-m_{\mathrm b}g$ である。無流束の壁の上で密度は $n(0)e^{-x/\ell}$ となり、$\ell=D/(\mu m_{\mathrm b}g)=k_{\mathrm B}T/(m_{\mathrm b}g)$ である。$m_{\mathrm b}$ は沈降に使う質量であり、慣性を記述する $m$ と区別する。
:::

:::example[水中の球]
連続体・粘着境界・低 Reynolds 数の Stokes 抵抗 $6\pi\eta a$ を用いる。半径 $a=0.5\,\mu\mathrm m$、水の粘性率を $10^{-3}\,\mathrm{Pa\,s}$、$T=300\,\mathrm K$ とすると $\gamma\simeq9.42\times10^{-9}\,\mathrm{kg\,s^{-1}}$、$D\simeq4.4\times10^{-13}\,\mathrm{m^2\,s^{-1}}$ である。三方向の拡散が独立なら $1\,\mathrm s$ での二乗平均平方根は $\sqrt{6Dt}\simeq1.6\,\mu\mathrm m$ である。抵抗則の流体力学は [粘性流体](/physics/continuum-mechanics/viscous-fluids/) で扱う。
:::

### ランダムウォークから連続の拡散へ

観測の刻みを $h$、一歩の長さを $\ell$ とする。各刻みで独立に $\pm\ell$ をとると、$N$ 歩後の変位 $Y_N$ の平均は零、分散は $N\ell^2$ である。$Nh=t$、$\ell^2=2Dh$ を保つ極限では、分散が $2Dt$ に保たれる。分布については、任意の波数 $k$ で特性関数が
$$
\mathbb E e^{ikY_N}=[\cos(k\sqrt{2Dh})]^N
=\exp\bigl[N(-Dk^2h+O(h^2))\bigr]\to e^{-Dk^2t}
$$
となる。最後の関数は平均零、分散 $2Dt$ の正規分布の特性関数である。互いに交わらない時刻区間では刻みの集合を分けられるため、有限次元分布の極限も独立増分をもつ。連続な経路の構成は [確率過程の基礎](/physics/nonequilibrium/stochastic-processes/) で扱う。

各刻みに平均 $\mu Fh$ を加えるには、一定力 $F$ のとき右へ進む確率を $1/2+\mu Fh/(2\ell)$ とすればよい。左へ進む確率との和は $1$ であり、$h$ が十分小さければ両者とも正になる。一歩の分散は $2Dh-(\mu Fh)^2$ なので、有限時刻の分散への補正は $O(h)$ で消える。一歩の特性関数は $\cos(k\ell)+i(\mu Fh/\ell)\sin(k\ell)=1+(ik\mu F-Dk^2)h+O(h^2)$ である。従って $t/h$ 歩の特性関数は $e^{ik\mu Ft-Dk^2t}$ へ収束し、極限は平均 $\mu Ft$、分散 $2Dt$ の正規分布になる。

:::example[確率の偏りと流れの向き]
$D=2\,\mu\mathrm m^2/\mathrm s$、$\mu F=1\,\mu\mathrm m/\mathrm s$、$h=0.01\,\mathrm s$ とする。一歩の長さは $0.2\,\mu\mathrm m$、右向きの確率は $0.525$ である。$100$ 歩後の平均は $1\,\mu\mathrm m$ であり、分散は $100(0.2^2-0.01^2)=3.99\,\mu\mathrm m^2$ である。連続近似の分散 $4\,\mu\mathrm m^2$ との差は、一歩の平均の二乗を引いた効果である。
:::

### 化学ポテンシャルを用いた流束則

理想希薄系の数密度に対する自由エネルギーは、基準数密度 $n_*>0$ を用いて
$$
\mathcal F[n]=\int\{k_{\mathrm B}T\, n[\log(n/n_*)-1]+Un\}\,\dd x
$$
と書ける。$n>0$ とし、任意の変分 $\delta n$ について $n\log n$ を微分すれば $\delta\mathcal F/\delta n=k_{\mathrm B}T\log(n/n_*)+U$ である。化学ポテンシャルの勾配に比例して粒子が流れ、その比例係数が単粒子の移動度と密度の積であると仮定する。すると
$$
J=-\mu n\partial_x\frac{\delta\mathcal F}{\delta n}
=-\mu k_{\mathrm B}T\,\partial_x n-\mu nU'
$$
となる。拡散の項は密度による化学ポテンシャルの変化から出る。壁で $J=0$ とすれば、部分積分から
$$
\dv{\mathcal F}{t}=\int J \partial_x\frac{\delta\mathcal F}{\delta n}\,\dd x
=-\int \frac{J^2}{\mu n}\,\dd x
$$
である。固定した総粒子数のもとで定常となるには化学ポテンシャルが一定になる必要があり、可積分な領域では Boltzmann 密度に戻る。

:::example[同じ拡散係数で異なる流束]
ある点で $n=2n_*/3$、$\partial_xn=n_*/(3L)$ とする。そこで $U'=k_{\mathrm B}T/(2L)$ なら、拡散と外力はそれぞれ $-Dn_*/(3L)$ を寄与する。$U'$ を $-k_{\mathrm B}T/(2L)$ に変えれば、二つの流れが相殺してその点の $J$ は零となる。定常性を判定するには一点での相殺に加え、境界条件と全領域での流束も調べる必要がある。
:::

### 連続初期分布と拡散核

点源の拡散を観測すると、密度の最大値は時間とともに低下する。ピークの低下と幅の増加を一つの式で結び付けるため、実直線上で時刻零に $x_0$ から出発する粒子を考える。密度を $p(t,x)$ とし、$\widehat p(t,k)=\int_{\R} e^{-ikx}p(t,x)\,\dd x$ と定義する。$\partial_t p=D\partial_x^2p$ を二回部分積分すると $\partial_t\widehat p=-Dk^2\widehat p$ となり、初期値は $e^{-ikx_0}$ である。従って $\widehat p(t,k)=e^{-ikx_0-Dk^2t}$ である。逆変換のため、$I(y)=\int_\R e^{-Dtk^2+iky}\dd k$ とおく。$y$ 微分と $k$ の部分積分により $I'=-yI/(2Dt)$ であり、$I(0)=\sqrt{\pi/(Dt)}$ だから
$$
G_t(y)=\frac{I(y)}{2\pi}=\frac{1}{\sqrt{4\pi Dt}}e^{-y^2/(4Dt)}
\label{eq-neq-brownian-heat-kernel}
$$
である。$G_t$ を :term[熱核]{en="heat kernel" yomi="ねつかく"} という。規格化は Gaussian 積分で従い、平均は対称性で零である。分散は $y G_t=-2Dt\,G_t'$ を使い、$\int y^2G_t=-2Dt\int yG_t'=2Dt$ と計算できる。任意の初期確率密度 $p_0$ からは $p(t,x)=\int_\R G_t(x-y)p_0(y)\dd y$ が得られる。実際、Fourier 変換は初期変換に $e^{-Dk^2t}$ を掛けるので、これは全波数で同じ初期値問題を解いている。

初期分布が幅を持つ実験では、$X_t=X_0+Y_t$ と分けられる。独立な $Y_t$ の分散は $2Dt$ であるから、初期分布の分散を測り忘れると初期幅まで拡散で広がったと解釈してしまう。一様に流れる媒質では、速度 $u$ のドリフトが加わり、初期分布に同じ核を畳み込んだ後で全体を $ut$ だけ移動させる。これは移動座標 $y=x-ut$ で $\partial_t p+u\partial_xp=D\partial_x^2p$ を熱方程式に戻すことからも確かめられる。

:::example[初期幅と流速を同時に推定する]
$\mu\mathrm m$ と秒を単位に用い、初期位置を平均 $1$、分散 $4$ の正規分布から選ぶ。$D=3$、$u=-2$ なら、$t=2$ の分布は平均 $-3$、分散 $16$ である。正規分布の和は正規分布であり、密度のピークは $1/\sqrt{32\pi}$ になる。逆に同じ初期分布を用いて時刻 $2$ に平均 $-3$、分散 $16$ が得られたなら、$u=(-3-1)/2=-2$、$D=(16-4)/4=3$ と復元する。幅の変化から流速を求めることはできず、平均も測る必要がある。
:::

### 壁を鏡像で表す

$0\le x<\infty$ で $x_0>0$ から始める場合、壁が粒子を反射するか吸収するかで境界条件が変わる。反射なら $J(0)=-D\partial_xp(0)=0$、吸収なら $p(0)=0$ を課す。:eqref[eq-neq-brownian-heat-kernel] の偶対称性を用い、二つの解を
$$
p_{\rm refl}(t,x)=G_t(x-x_0)+G_t(x+x_0),\qquad
p_{\rm abs}(t,x)=G_t(x-x_0)-G_t(x+x_0)
$$
と構成する。いずれの項も熱方程式を満たす。和の微分は $x=0$ で零になり、差は $x=0$ で零になる。初期時刻では第二項の点源が領域の外にあるので、どちらも $x_0$ での点源に収束する。この方法は :term[鏡像法]{en="method of images" yomi="きょうぞうほう"} である。十分な減衰を持つ二つの滑らかな解が同じ初期状態から出発したなら、その差 $q$ の二乗積分は $\dv{}{t}\int q^2=-2D\int(q')^2$ となるので、境界条件をそろえれば一意である。点源の場合は滑らかな初期分布を点へ集中させた極限で解を読む。

:::example[吸収壁で失った総数]
$D=1$、$x_0=2$ とし、時刻 $t=2$ を考える。吸収されず残った確率を求めるには密度の差を半直線で積分する。一つ目に $y=x-x_0$、二つ目に $y=x+x_0$ を使うと差は $\int_{-x_0}^{x_0}G_t(y)\dd y$ になる。この例の標準偏差は $\sqrt{2Dt}=2$ なので残存確率は $2\Phi(1)-1\simeq0.683$ である。$\Phi$ は標準正規分布の累積分布関数である。反射なら、第二項が第一項の壁の外に出た分を補うので、残存確率は常に $1$ である。
:::

:::example[吸収壁に達する時刻の密度]
吸収解を微分すると $\partial_x p_{\rm abs}(t,0)=x_0G_t(x_0)/(Dt)$ である。左向きの流出量は $D\partial_xp_{\rm abs}(t,0)=x_0G_t(x_0)/t$ となる。$t$ に関する対数微分は $-3/(2t)+x_0^2/(4Dt^2)$ であり、最大は $t=x_0^2/(6D)$ で生じる。例えば $D=1$、$x_0=2$ なら最頻時刻は $2/3$ である。一方、密度の尾は $t^{-3/2}$ であり、分布全体の平均時刻は発散する。最も多く観測される時刻と平均時刻は異なる。
:::

## Langevin 方程式とその解

### 熱浴を表す確率方程式

液体の力を、平均的な抵抗と平均零の残差に分ける。熱浴の相関時間が粒子の速度緩和時間より短い領域を考え、粒子に働く残差の積分を Wiener 過程 $W_t$ で表す。$W_0=0$ とし、$W_t-W_s$ は平均零、分散 $t-s$ の正規変数である。互いに交わらない時間区間の増分は独立であり、初期状態も将来の増分と独立に選ぶ。次式を :term[Langevin 方程式]{en="Langevin equation" yomi="らんじゅばんほうていしき"} と呼ぶ。

$$
\dd X_t=V_t\,\dd t,\qquad m\,\dd V_t=(-\gamma V_t+F(t))\,\dd t+\sqrt{2A}\,\dd W_t. \label{eq-neq-langevin}
$$

ここで $\gamma>0$、$A\ge0$、$F(t)$ は決定論的な力である。雑音の力を形式的に $\xi=\sqrt{2A}\dot W$ と書くと $\mathbb E[\xi(t)\xi(s)]=2A\delta(t-s)$ である。$\delta$ は時間の逆数の次元をもつので $A$ の単位は $\mathrm{N^2\,s}$ である。

### 積分因子と変位

$\tau=m/\gamma$ とする。$e^{t/\tau}$ は滑らかな決定論的関数なので、積の微分で雑音との二次変分は零である。よって :eqref[eq-neq-langevin] に積分因子をかけると

$$
V_t=e^{-t/\tau}V_0+\frac1m\int_0^t e^{-(t-s)/\tau}F(s)\,\dd s+\frac{\sqrt{2A}}m\int_0^t e^{-(t-s)/\tau}\,\dd W_s. \label{eq-neq-brownian-velocity}
$$

二つの解の差は初期差を $e^{-t/\tau}$ 倍した値に等しい。したがって同じ雑音と初期値を用いた解は一意である。確定した被積分関数 $g$ に対し正規増分の和の分散は $\sum g(s_i)^2\Delta s_i$ であり、極限で $\mathbb E[(\int g\,\dd W)^2]=\int g^2\,\dd s$ を得る。このため固定した $V_0$ から出発した速度は正規分布で、分散は $A(1-e^{-2t/\tau})/(m\gamma)$ である。

位置を求めるには、もとの速度方程式を積分して $X_t-X_0=\int_0^t V_s\,\dd s$ を消去する。すると

$$
X_t-X_0=\tau(V_0-V_t)+\gamma^{-1}\int_0^t F(s)\,\dd s+\frac{\sqrt{2A}}{\gamma}W_t. \label{eq-neq-brownian-displacement}
$$

自由粒子で $V_0$ を固定すれば、変位の雑音部は $\sqrt{2A}\gamma^{-1}\int_0^t(1-e^{-(t-s)/\tau})\,\dd W_s$ となる。変位と速度は同じ増分を含むので、一般に独立でない。

:::example[静止から放す試行]
$V_0=0$ として $F=0$ のとき、短時間で変位分散は $2A t^3/(3m^2)+O(t^4)$ となる。実際、積分核の先頭は $(t-s)/\tau$ であり、その二乗を積分すると $t^3/(3\tau^2)$ になる。次節の弾道的な $t^2$ は初速度も熱平衡から選んだ場合の結果である。
:::

:::example[力積による初速の変化]
時刻零で力積 $I$ を加えると $V_0$ に $I/m$ が加わる。応答の差は平均速度で $I e^{-t/\tau}/m$、平均位置で $I(1-e^{-t/\tau})/\gamma$ である。後者の極限は力積を摩擦係数で割った長さである。
:::

### 位置と速度を同時に標本化する

時間刻み $h>0$ で、左端の状態 $(X,V)=(x,v)$ から自由粒子を進める。ここで $D=A/\gamma^2$ と書く。この値が長時間の拡散係数になることは次節で示す。$r=e^{-h/\tau}$ とすると条件付き平均は $(x+\tau(1-r)v,rv)$ である。雑音部を $(\eta_x,\eta_v)$ と書くと、両者は同じ Wiener 増分の積分である。核の積を積分すれば、離散時刻の正確な標本生成に使う共分散は

$$
\begin{aligned}
\operatorname{Var}\eta_v&=(D/\tau)(1-r^2),\\
\operatorname{Cov}(\eta_x,\eta_v)&=D(1-r)^2,\\
\operatorname{Var}\eta_x&=2D[h-2\tau(1-r)+\tau(1-r^2)/2]
\end{aligned}
$$

となる。交差項については $(2D/\tau)\int_0^h(1-e^{-s/\tau})e^{-s/\tau}\dd s=D[2(1-r)-(1-r^2)]$ である。$\eta_x,\eta_v$ を独立に生成すると同じ時刻の分散は合うが、次の時刻へ持ち越す相関が変わる。$D=0$ なら雑音部は両方とも零である。$D>0$ なら標準正規変数 $Z_1,Z_2$ を独立に選び、$\eta_v=s_vZ_1$、$\eta_x=c_{xv}Z_1/s_v+\sqrt{s_x^2-c_{xv}^2/s_v^2}\,Z_2$ とすればよい。平方根の中の非負性は積分核に対する Cauchy–Schwarz の不等式から従う。

:::example[緩和時間より短い標本]
$h/\tau\ll1$ では $\operatorname{Var}\eta_v\simeq 2Dh/\tau^2$、$\operatorname{Var}\eta_x\simeq 2Dh^3/(3\tau^2)$、共分散は $Dh^2/\tau^2$ である。従って雑音部の相関係数は $\sqrt3/2$ に近づく。短い刻みほど相関を無視できるとは限らない。
:::

### 線形系の遷移行列と外力の畳み込み

調和力 $-\kappa X$ と駆動 $f(t)$ を同時に加える。$Z=(X,V)^\mathsf T$ とし、
$$
M=\begin{pmatrix}0&1\\-\kappa/m&-\gamma/m\end{pmatrix},
\qquad B=\binom{0}{\sqrt{2A}/m},\qquad b(t)=\binom{0}{f(t)/m}
$$
とおけば、$\dd Z=(MZ+b)\dd t+B\dd W$ である。$e^{-Mt}$ を掛けて積分すると、平均 $e^{Mt}Z_0+\int_0^te^{M(t-s)}b(s)\dd s$、条件付き共分散 $\int_0^te^{Ms}BB^\mathsf T e^{M^\mathsf Ts}\dd s$ をもつ正規遷移が得られる。ここでは $Z_0$ を固定した。初期値も確率変数で、未来の雑音と独立なら、共分散には初期共分散 $\Sigma_0$ から伝わった $e^{Mt}\Sigma_0 e^{M^\mathsf Tt}$ を加える。初期値の分布が非正規なら一般に遷移後の周辺分布も非正規である。

$C(t)$ を共分散とする。積分の上端を微分し、$\dv{}{s}(e^{Ms}BB^\mathsf Te^{M^\mathsf Ts})=M e^{Ms}BB^\mathsf Te^{M^\mathsf Ts}+e^{Ms}BB^\mathsf Te^{M^\mathsf Ts}M^\mathsf T$ を用いると
$$
\dot C=MC+CM^\mathsf T+BB^\mathsf T
$$
となる。$\kappa,\gamma>0$ なら $M$ の根の和は負、積は正なので、実根でも複素根でも実部が負となる。従って初期共分散の項は消え、積分は無限大まで収束する。この積分を $C_\infty$ とすれば、有限刻みのノイズ共分散は $C_\infty-e^{Mh}C_\infty e^{M^\mathsf Th}$ である。

:::example[初期値のばらつきと浴のばらつき]
自由粒子で $m=\gamma=k_{\mathrm B}T=1$ とし、$V_0=\pm2$ を等確率でとる。$t=\log 2$ では条件付き速度分布が平均 $\pm1$、分散 $3/4$ の正規分布である。二つを混ぜた速度の分散は $3/4+1=7/4$ である。四次モーメントは $1+6(3/4)+3(3/4)^2=115/16$ となる。分散 $7/4$ の単一の正規分布なら $147/16$ になるため、分散が分かっても過渡的な分布は決まらない。
:::

:::example[力の印加と遮断]
自由粒子を平均速度零から出発させ、$0<t<t_1$ に一定力 $F$ を加える。$t\le t_1$ の平均速度は $(F/\gamma)(1-e^{-t/\tau})$ であり、$t\ge t_1$ では遮断直前の値に $e^{-(t-t_1)/\tau}$ を掛ける。遮断後の変位を積分すると $(Ft_1/\gamma)-\mathbb E[X_{t_1}-X_0]=(F\tau/\gamma)(1-e^{-t_1/\tau})$ である。例えば $t_1=\tau\log2$ では、全変位 $F\tau\log2/\gamma$ のうち $F\tau/(2\gamma)$ が遮断後に生じる。
:::

### Euler 法の定常分散の偏り

時間発展を数値的に追う場合、力と摩擦を各刻みの左端だけで評価する方法がある。自由粒子の速度に対する Euler–Maruyama 更新は $V_{n+1}=(1-q)V_n+\sqrt{2Ah}\,Z_n/m$、$q=h/\tau$ である。$Z_n$ は独立な標準正規変数である。平均偏差は $r_E=1-q$ を刻みごとに掛けた値だから、$0<q<2$ が緩和の条件である。$A>0$ かつこの条件の下で分散 $C_n$ は $C_{n+1}=r_E^2 C_n+2Ah/m^2$ に従う。反復して幾何級数を足すと、定常値は
$$
C_E=\frac{2Ah/m^2}{1-(1-q)^2}=\frac{A}{m\gamma}\frac1{1-q/2}
$$
である。従って分散は真の値より大きい。$q$ が $1$ を超えると、真の速度相関は常に正であるにもかかわらず、交互に符号が変わる。正確な遷移では係数が $r=e^{-q}>0$ であり、雑音の分散も $A(1-r^2)/(m\gamma)$ に合わせるので偏りがない。

刻みを固定して位置を $X_{n+1}=X_n+hV_n$ と進め、定常な Euler 速度の相関を用いる。拡散係数の離散版は
$$
D_E=\frac h2 C_E\left(1+2\sum_{j\ge1}r_E^j\right)
=\frac h2 C_E \frac{1+r_E}{1-r_E}=\frac A{\gamma^2}.
$$
第一式は $h\sum_{j=0}^{N-1}V_j$ の分散を $2Nh$ で割り、時差ごとの個数を数えて得た。定常分散の偏りと相関のずれが長時間の傾きでは相殺するため、拡散だけの比較では不適切な刻みを見抜けない。

:::example[拡散は合うが速度は合わない刻み]
$m=\gamma=A=1$ とする。$h=3/2$ なら $r_E=-1/2$、$C_E=4$ であり、一刻み遅れの共分散は $-2$ である。真の共分散は $e^{-3/2}>0$ で、定常分散も $1$ である。両計算とも $D=1$ を与えるが、短時間の二乗変位や速度分布を評価するには正確な更新か小さな刻みが必要である。
:::

### 線形安定な拘束での正確な遷移

線形系の遷移行列には、固有値が一致する場合も含めて使える表式がある。$M$ を調和拘束で導入した行列とし、$c=\gamma/(2m)$、$\omega_0^2=\kappa/m$、$N=M+cI$ とおく。行列の積を計算すると $N^2=(c^2-\omega_0^2)I$ である。$\delta^2=c^2-\omega_0^2$ として指数関数の級数を偶数項と奇数項に分けると
$$
e^{Mh}=e^{-ch}\left[I\cosh(\delta h)+N\frac{\sinh(\delta h)}{\delta}\right]
$$
となる。$\delta=0$ では $\sinh(\delta h)/\delta\to h$ と読む。$\delta$ が純虚数なら、$\cosh$ と $\sinh/\delta$ は実の $\cos$ と $\sin$ の組になる。固有ベクトルを求めずに計算できるので、二つの根が接近したときにも同じ式を用いられる。

平衡浴なら $C_\infty=\operatorname{diag}(k_{\mathrm B}T/\kappa,k_{\mathrm B}T/m)$ であることは、$MC_\infty+C_\infty M^\mathsf T+BB^\mathsf T=0$ に代入すれば確認できる。正確な刻みでは、発生させる雑音共分散 $C_h$ が $C_\infty-e^{Mh}C_\infty e^{M^\mathsf Th}$ である。この差が数値的に小さいときは積分表示に戻る。特に位置方向の分散は短い刻みで $O(h^3)$、速度方向は $O(h)$ なので、桁落ちの影響が両者で違う。

:::example[臨界減衰での共分散]
$m=\kappa=k_{\mathrm B}T=1$、$\gamma=2$ とすると $C_\infty=I$、$c=1$、$N=\begin{pmatrix}1&1\\-1&-1\end{pmatrix}$ である。$h=1$ なら $e^{Mh}=e^{-1}\begin{pmatrix}2&1\\-1&0\end{pmatrix}$ であり、$C_h=I-e^{-2}\begin{pmatrix}5&-2\\-2&1\end{pmatrix}$ になる。従って交差共分散は $2e^{-2}$ である。条件付き平均で予測した運動にこの相関を持つ雑音を足すと、定常共分散から始めた場合にも、次の刻みで定常共分散が正確に再現される。
:::

## 速度相関関数と拡散係数

### 相関を二回積分する

速度を定常な平均零の分布から始める。自由粒子について :eqref[eq-neq-brownian-velocity] の未来の雑音は過去と独立なので、$s\ge0$ に対して $C_V(s)=\mathbb E[V_{t+s}V_t]=C_V(0)e^{-s/\tau}$ となる。分散を時刻によらず一定とおけば $C_V(0)=A/(m\gamma)$ を得る。実確率変数の積は可換なので定常相関は偶関数である。

:::proposition[変位分散と速度相関]{#prop-neq-brownian-msd}
定常な平均零の速度が二次モーメントをもち、相関 $C_V$ が正の半直線で絶対可積分なら

$$
\mathbb E[(X_t-X_0)^2]=2\int_0^t (t-s)C_V(s)\,\dd s,\qquad
D=\int_0^\infty C_V(s)\,\dd s. \label{eq-neq-brownian-green-kubo}
$$

:::

:::proof
$\mathbb E[(\int_0^t V_s\,\dd s)^2]$ を時刻の正方形上で積分する。Cauchy–Schwarz により被積分関数の絶対期待値は $C_V(0)$ 以下なので積分順を交換できる。正方形を対角線で分け、上半分で時差を $s$ とすると、同じ時差に対応する区間の長さは $t-s$ である。これで第一式が得られる。$2t$ で割ると重み $(1-s/t)\mathbf1_{\{s<t\}}$ が点ごとに $1$ に収束する。$|C_V|$ による優収束で第二式となる。
:::

指数相関を代入し、$\int_0^t(t-s)e^{-s/\tau}\,\dd s=\tau[t-\tau(1-e^{-t/\tau})]$ を用いると $D=A/\gamma^2$ である。変位分散の正確な式は $2D[t-\tau(1-e^{-t/\tau})]$ となる。$t\ll\tau$ で $C_V(0)t^2$、$t\gg\tau$ で $2Dt-2D\tau+o(1)$ である。

:::example[有限時間の傾き]
二乗変位を $2$ で割ったものの瞬間の傾きは $D(1-e^{-t/\tau})$ である。$t=\tau\log 2$ では傾きが $D/2$ になる。実験で有限時刻の傾きから $D$ を読むときには、速度がまだ相関している時間帯を避ける。
:::

:::example[方向ごとの違い]
摩擦係数が主軸方向で $\gamma_1,\gamma_2$、雑音が独立で温度が等しければ、各方向の $D_i$ は異なる。二次元の全二乗変位の傾きは $2(D_1+D_2)$ である。等方的な場合に限り、これを $4D$ と書ける。
:::

### 変位の相関と位置測定誤差

平衡自由粒子では観測刻み $\Delta$ で集めた変位 $Y_j=X_{(j+1)\Delta}-X_{j\Delta}$ は定常である。ただし、$Y_0$ と $Y_n$ の独立性は速度の相関が残る限り成立しない。$n\ge1$ に対して二つの区間の順序が定まり、
$$
\begin{aligned}
\mathbb E[Y_0Y_n]
&=\frac D\tau\int_0^\Delta\dd s\int_{n\Delta}^{(n+1)\Delta}\dd u\, e^{-(u-s)/\tau}\\
&=D\tau(1-e^{-\Delta/\tau})^2e^{-(n-1)\Delta/\tau}.
\end{aligned}
$$
従って長い軌道の変位を独立な標本とみなして平均や分散の誤差を評価すると、誤差を過小評価しうる。和の分散では対角成分が $N$ 個、時差 $n$ の非対角成分が往復で $2(N-n)$ 個あるので、両方を足す。$N\to\infty$ で $N$ 当たりの分散は $\mathbb E[Y_0^2]+2\sum_{n\ge1}\mathbb E[Y_0Y_n]$ となる。幾何級数を足すと、第一項が含む $-2D\tau(1-e^{-\Delta/\tau})$ と第二項が相殺し、結果は $2D\Delta$ である。

:::example[位置測定誤差による反相関]
独立な平均零、分散 $s^2$ の測定誤差 $\epsilon_j$ が位置に入ると、観測変位は $Y_j+\epsilon_{j+1}-\epsilon_j$ となる。隣接する誤差差分の共分散は $-s^2$ であり、二つ以上離れた差分では共有項がない。$D=\tau=1$、$\Delta=\log2$、$s^2=1/2$ とすると、真の隣接共分散は $1/4$ であり、観測した共分散は $-1/4$ である。この符号だけで拘束や反発力を推定すると、測定誤差と力の効果を混同する。
:::

:::example[長時間尾が積分できない場合]
定常速度の相関が $C_V(s)=c(1+s/t_0)^{-1}$ だとする。:eqref[eq-neq-brownian-green-kubo] の有限時間の式に入れ、$u=1+s/t_0$ とおくと分散は $2ct_0[(t+t_0)\log(1+t/t_0)-t]$ となる。$t\gg t_0$ では $t\log t$ の増加であり、有限の拡散係数は得られない。一方、尾が $c(1+s/t_0)^{-2}$ であれば積分は $ct_0$ となる。指数減衰かどうかよりも、時間積分の収束が長時間の傾きを決める。
:::

### 位置観測の時間平均と独立な標本の数

実験では速度を瞬時に測るのが難しく、変位の時間平均から速度を推定する。定常な自由粒子の区間速度 $\bar V_\Delta=(X_{t+\Delta}-X_t)/\Delta$ は正規であり、分散は
$$
\operatorname{Var}\bar V_\Delta
=\frac{2D}{\Delta}\left[1-\frac{\tau}{\Delta}(1-e^{-\Delta/\tau})\right]
$$
となる。$\Delta\to0$ では $D/\tau=k_{\mathrm B}T/m$ に戻り、長い区間では $2D/\Delta$ に近づく。長い区間では速度を平均するため、瞬時の熱平衡の分散から大きく減る。

ラグ $\Delta$ を固定し、非重複な $N$ 区間の二乗変位の平均 $\hat M=N^{-1}\sum_{j=0}^{N-1}Y_j^2$ を考える。分布が中心化した同時正規であれば、$\operatorname{Cov}(Y_i^2,Y_j^2)=2\operatorname{Cov}(Y_i,Y_j)^2$ である。これを示すため、分散 $M>0$ を共通とし、$Y_j=cY_i/M+Z$ と置く。共分散 $c=\operatorname{Cov}(Y_i,Y_j)$ としたので、$Z$ と $Y_i$ の相関は零であり、正規性により独立である。$\mathbb EY_i^4=3M^2$ と $\mathbb E Z^2=M-c^2/M$ を代入すると $\mathbb E[Y_i^2Y_j^2]=M^2+2c^2$ になる。よって
$$
\operatorname{Var}\hat M=\frac{2M^2}{N}+\frac4{N^2}\sum_{j=1}^{N-1}(N-j)c_j^2,
\qquad c_j=\operatorname{Cov}(Y_0,Y_j).
$$
二乗変位の誤差には共分散の二乗が入る。従って測定誤差が作った負の共分散も、二乗変位の平均の不確かさを大きくする。

:::example[二区間を一つずつ平均する]
$D=\tau=1$、$\Delta=\log2$ とする。$M=2\log2-1$、$c_1=1/4$ が得られる。二区間の二乗変位の平均では、独立なら分散は $M^2$ であるが、実際は $M^2+1/16$ である。比は約 $1.419$ である。例えば独立な試行ごとに一つの変位を採る設計なら、この追加分は出ない。
:::

### 指数相関と複数の緩和過程

速度相関を $C(t)=\sum_{j=1}^N a_j e^{-|t|/\tau_j}$、$a_j,\tau_j>0$ とする。これは独立な定常 Ornstein–Uhlenbeck 速度成分を足す模型から得られる。:eqref[eq-neq-brownian-green-kubo] に項ごとに代入すると、
$$
M(t)=2\sum_j a_j\tau_j[t-\tau_j(1-e^{-t/\tau_j})],\qquad D=\sum_j a_j\tau_j
$$
である。弾道の係数は $C(0)=\sum_j a_j$、長時間の切片は $-2\sum_j a_j\tau_j^2$ である。従って短時間の二乗変位と長時間の傾きだけから単一時定数を合わせると、切片が合うとは限らない。実際、単一時定数の推定値は $\bar\tau=\sum_j a_j\tau_j/\sum_j a_j$ である。Cauchy–Schwarz から $\sum_j a_j\tau_j^2\ge(\sum_j a_j\tau_j)^2/\sum_j a_j$ なので、切片の絶対値は単一指数で置いた値以上になる。等号はすべての時定数が同一のときに限る。

:::example[始めと終わりを合わせても途中が変わる]
$a_1=a_2=1$、$\tau_1=1$、$\tau_2=3$ とする。$C(0)=2$、$D=4$ なので単一指数なら $\bar\tau=2$ である。しかし長時間の切片は二時定数なら $-2(1+9)=-20$、単一指数なら $-2\cdot2\cdot2^2=-16$ となる。時刻の原点や初期位置を誤差の原因にする前に、複数の緩和時間が候補になる。
:::

## 揺動散逸関係

### 散逸と雑音強度を照合する

摩擦だけなら速度は零に向かう。温度 $T$ の熱浴と平衡にある速度は、運動エネルギー $m v^2/2$ によって決まる Maxwell 分布をもつ。正規積分を微分すれば $\mathbb E[mV^2]=k_{\mathrm B}T$ を得る。前節の定常分散を等置すると

$$
A=\gamma k_{\mathrm B}T,\qquad
\mathbb E[\xi(t)\xi(s)]=2\gamma k_{\mathrm B}T\,\delta(t-s). \label{eq-neq-brownian-fdt}
$$

この関係は :term[揺動散逸関係]{en="fluctuation–dissipation relation" yomi="ようどうさんいつかんけい"} の一例である。雑音は摩擦によって粒子が失ったエネルギーを平均として補う。

検算として、二乗の微分を増分で計算する。$\dd(V^2)=2V\,\dd V+(\dd V)^2$ とし、$(\dd W)^2=\dd t$ を用いると、自由粒子のエネルギーは

$$
\dv{}{t}\mathbb E\!\left[\frac{mV^2}{2}\right]=-\gamma\mathbb E[V^2]+\frac{A}{m}
$$

を満たす。二次変分の項を落とすと雑音からのエネルギー流入が消えてしまう。雑音の強さを測定値から自由に選べる場合もあるが、平衡浴を表すなら :eqref[eq-neq-brownian-fdt] を満たす必要がある。

:::example[雑音の強度と平衡温度]
$\gamma$ と $m$ を固定し、$A$ だけを二倍にする。定常速度分散は二倍、標準偏差は $\sqrt2$ 倍である。単一の自由粒子の定常分布から推定する温度は二倍になる。
:::

:::example[保存力を含めた平衡]
$F=-U'(x)$ とする。位相空間密度 $\rho\propto e^{-\beta(mv^2/2+U)}$ について、保存力による微分 $-v\partial_x\rho+(U'/m)\partial_v\rho$ は相殺する。速度の散逸流 $-(\gamma/m)v\rho-(A/m^2)\partial_v\rho$ も $A=\gamma/\beta$ で零になる。位置の流れ $v\rho$ が位相空間で零になる必要はなく、平衡での時間反転は $v$ の符号も変える。
:::

### 摩擦が異方的な浴と反対称な力

$d$ 次元で $m$ をスカラー、摩擦を実対称正定値行列 $G$ とする。速度には保存的な速度依存力 $Hv$ も働き、$H^\mathsf T=-H$ とする。この力の仕事率は $v^\mathsf THv=0$ である。$\dd V=m^{-1}(-G+H)V\dd t+m^{-1}B\dd W$ の共分散 $C$ を Itô の積の公式で微分すると、
$$
\dot C=m^{-1}[(-G+H)C+C(-G-H)]+m^{-2}BB^\mathsf T
$$
となる。熱平衡の $C=(k_{\mathrm B}T/m)I$ を代入すると、雑音の共分散は $BB^\mathsf T=2k_{\mathrm B}TG$ と決まる。反対称力の寄与は二項で打ち消し合う。任意の複素固有ベクトル $z$ に対し、固有値の実部は $-z^\dagger G z/(m z^\dagger z)<0$ なので、この共分散への緩和が起きる。

:::example[横方向の応答と拡散]
二次元で $G=\gamma I$、$H=b\begin{pmatrix}0&-1\\1&0\end{pmatrix}$ とする。一定力に対する定常移動度は $(G-H)^{-1}=(G+H)/(\gamma^2+b^2)$ である。位置の長時間共分散は、速度相関を積分した行列とその転置の和に時間を掛けたものである。実際、正の時差の相関は $e^{(-G+H)t/m}C$ であり、積分は $k_{\mathrm B}T(G-H)^{-1}$ となる。従って拡散行列は $[k_{\mathrm B}T\gamma/(\gamma^2+b^2)]I$ となる。$b=\gamma$ なら縦移動度は $1/(2\gamma)$、横移動度の大きさも $1/(2\gamma)$ であるが、横移動度は位置の共分散に寄与しない。
:::

### 複数の熱浴によるエネルギー流

同じ速度に浴 $i$ が独立な雑音と摩擦を加え、各浴で $A_i=\gamma_i k_{\mathrm B}T_i$ とする。全速度分散の定常値は $(\sum_i A_i)/(m\sum_i\gamma_i)$ である。Itô の式における浴ごとの寄与は、平均で $-\gamma_i\mathbb E[V^2]+A_i/m$ となる。従って速度だけの定常正規分布に温度を割り当てられても、浴ごとの平均エネルギー流が零になるとは限らない。

:::example[同じ摩擦を持つ二浴]
$m=\gamma_1=\gamma_2=k_{\mathrm B}=1$、$T_1=3$、$T_2=1$ とする。定常分散は $2$ なので、高温浴から粒子に流れ込むエネルギーは単位時間に $3-2=1$、低温浴からの流入は $1-2=-1$ である。粒子の平均エネルギーは一定だが、浴の合計エントロピーは単位時間に $-1/3+1=2/3$ だけ増える。
:::

### 位相空間で可逆な流れと不可逆な流れを分ける

速度を含む Langevin 過程では、平衡でも $x$ 方向の流れが $v p$ である。この流れをそのまま非平衡の散逸と読むと、運動を続ける平衡粒子まで散逸していることになる。有界な領域での鏡面反射または十分速い減衰を仮定し、平衡密度を $\pi(x,v)=Z^{-1}e^{-\beta(mv^2/2+U(x))}$ とする。前進方程式を
$$
\partial_t p=L_H^*p-\partial_v J_{\rm irr},\qquad
L_H^*p=-v\partial_xp+(U'/m)\partial_vp,\qquad
J_{\rm irr}=-\frac{\gamma}{m}vp-\frac{A}{m^2}\partial_vp
$$
と書く。$L_H^*$ は発散零の Hamiltonian 流であり、$\pi$ を保つ。$A=\gamma/\beta$ で $J_{\rm irr}=-(A/m^2)p\,\partial_v\log(p/\pi)$ である。相対エントロピー $\mathcal D(p\Vert\pi)=\int p\log(p/\pi)$ を微分し、前進式を代入する。Hamiltonian の項は、流に沿う微分 $L_H$ が $\log\pi$ を消し、$\int L_H p=0$ になるため零である。散逸流の項には部分積分で $J_{\rm irr}\partial_v\log(p/\pi)$ が残り、
$$
\dv{}{t}\mathcal D(p\Vert\pi)
=-\frac A{m^2}\int p\left|\partial_v\log(p/\pi)\right|^2\dd x\,\dd v
$$
となる。従って相対エントロピーを直接減らすのは速度方向の流れである。式の右辺がある一時刻に零でも $p=\pi$ とは限らない。$p/\pi$ が $x$ の関数なら速度微分は零だが、Hamiltonian 流が位置方向の差を速度の差へ移す。

:::example[位置だけずらして用意した状態]
調和井戸で速度は平衡分布、位置は同じ平衡分散を持つが平均が $q_0\ne0$ の正規分布とする。力学は線形なので、以後も共分散は平衡値であり平均だけが $q(t),\dot q(t)$ へ移る。二つの同分散 Gaussian の相対エントロピーを計算すると $\beta[\kappa q(t)^2+m\dot q(t)^2]/2$ である。微分は $\beta\dot q(\kappa q+m\ddot q)=-\beta\gamma\dot q^2$ になる。例えば $\dot q(0)=0$ なら初期微分は零だが、初期加速度 $-\kappa q_0/m$ が非零なので運動の後に散逸が起きる。
:::

### 時間反転の符号を生成作用素から読む

平衡測度の下で実の観測量の積を平均し、内積 $\langle f,g\rangle_\pi=\int fg\pi$ とする。後退生成作用素を $L=L_H+L_D$ に分けると、
$$
L_H f=v\partial_xf-(U'/m)\partial_vf,\qquad
L_D f=(A/m^2)\pi^{-1}\partial_v(\pi\partial_vf)
$$
である。第一項は位相流の発散が零で $\pi$ を保つから反対称、第二項は部分積分して対称である。速度反転 $R f(x,v)=f(x,-v)$ は $RL_HR=-L_H$、$RL_DR=L_D$ を満たす。従って重み付き随伴は $L^\dagger=RLR$ である。十分滑らかな試験関数と、これらの演算が保たれる拡散半群の領域で、$\langle f,e^{tL}g\rangle_\pi=\langle Rg,e^{tL}Rf\rangle_\pi$ が従う。これは平衡でも速度を反転する必要がある理由を示す。反射壁がある場合は、壁の反射則も速度反転と可換であることを要する。

:::example[位置と速度の遅れの符号]
調和拘束中の定常粒子で、$f=x$ は偶、$g=v$ は奇である。従って $\langle X_0 V_t\rangle=-\langle V_0 X_t\rangle$ である。$m=\kappa=1$、$\gamma=2$ の臨界減衰なら $C_X(t)=k_{\mathrm B}T(1+t)e^{-t}$ である。時刻 $t=1$ で $\langle X_0V_1\rangle=C_X'(1)=-k_{\mathrm B}T/e$、逆の組は $+k_{\mathrm B}T/e$ である。位置同士の相関は偶でも、位置と速度の交差相関は同じ交換に対して符号が変わる。
:::

### 調和振動子の浴から摩擦と雑音を導く

粒子の座標と運動量を $x,P$、浴の振動子の座標と運動量を $q_j,p_j$ とする。質量 $m_j>0$、角振動数 $\omega_j>0$ の振動子を粒子に線形に結合し、
$$
H=\frac{P^2}{2m}+U(x)+\sum_j
\left[\frac{p_j^2}{2m_j}+\frac{m_j\omega_j^2}{2}(q_j-\alpha_j x)^2\right]
$$
とする。平方の形を採ったので、$x$ を固定して浴の Gibbs 重みを積分しても、その積分は $x$ によらず、粒子の平衡ポテンシャルは $U$ である。正準方程式は $m\ddot x=-U'(x)+\sum_j m_j\omega_j^2\alpha_j(q_j-\alpha_j x)$ と $\ddot q_j+\omega_j^2 q_j=\omega_j^2\alpha_j x$ である。後者を強制振動として解き、$\omega_j\sin(\omega_j(t-s))=\partial_s\cos(\omega_j(t-s))$ で部分積分すれば
$$
q_j(t)-\alpha_j x(t)
=(q_j(0)-\alpha_j x(0))\cos\omega_jt
+\frac{p_j(0)}{m_j\omega_j}\sin\omega_j t
-\alpha_j\int_0^t\cos(\omega_j(t-s))\dot x(s)\dd s .
$$
従って粒子方程式の遅延摩擦は $-\int_0^t K(t-s)\dot x(s)\dd s$ であり、$K(t)=\sum_j m_j\omega_j^2\alpha_j^2\cos\omega_j t$ となる。残った二つの初期値の項が力雑音 $\xi(t)$ になる。

浴を粒子の初期位置に対する条件付き Gibbs 分布から選ぶと、$q_j(0)-\alpha_j x(0)$ と $p_j(0)/(m_j\omega_j)$ は独立で、分散はどちらも $k_{\mathrm B}T/(m_j\omega_j^2)$ である。異なる $j$ の相関は零だから、$\cos\omega_jt\cos\omega_js+\sin\omega_jt\sin\omega_js=\cos(\omega_j(t-s))$ を用いて $\langle\xi(t)\xi(s)\rangle=k_{\mathrm B}T K(t-s)$ を得る。有限個の振動子だけでは $K$ は一般に減衰しない。白色近似には、浴の周波数が広く分布し、粒子の速度が変わる時間より短い時間に相関が実効的に減衰する条件が要る。

### 記憶のある熱浴を補助変数で表す

摩擦が遅れて働く模型を、速度と浴の変数を合わせて解く。自由粒子に対し $\gamma,b,m>0$ とし、
$$
m\,\dd V=(-Y+F(t))\dd t,\qquad
\dd Y=(\gamma bV-bY)\dd t+b\sqrt{2\gamma k_{\mathrm B}T}\,\dd W
$$
とする。無外力で平均零、共分散 $\diag(k_{\mathrm B}T/m,\gamma b k_{\mathrm B}T)$ の同時正規分布から出発する。積の微分は $\dv{}{t}\langle V^2\rangle=-2\langle VY\rangle/m$、$\dv{}{t}\langle VY\rangle=\gamma b\langle V^2\rangle-\langle Y^2\rangle/m-b\langle VY\rangle$、$\dv{}{t}\langle Y^2\rangle=2\gamma b\langle VY\rangle-2b\langle Y^2\rangle+2\gamma b^2k_{\mathrm B}T$ を与える。指定した共分散では三つとも零であり、線形方程式で正規性も保たれるため、この初期分布は定常である。

$Y$ を積分因子で解いて速度の式に代入すると、摩擦は $-\int_0^t K(t-s)V(s)\dd s$、$K(t)=\gamma b e^{-bt}$ である。残りの力は
$$
\xi(t)=-e^{-bt}Y_0-b\sqrt{2\gamma k_{\mathrm B}T}\int_0^t e^{-b(t-s)}\dd W_s .
$$
$\xi$ の等時刻分散は $\gamma b k_{\mathrm B}T$ である。$t\ge s$ で $\xi(t)=e^{-b(t-s)}\xi(s)$ に独立な平均零の積分が加わるので、$\langle\xi(t)\xi(s)\rangle=k_{\mathrm B}T K(|t-s|)$ となる。$Y_0$ の分散が過去の雑音に相当するため、$Y_0=0$ に固定するとこの相関は非定常になる。

無外力の速度方程式に独立な $V_0$ を掛けて平均すると、$C(t)=\langle V_tV_0\rangle$ は $m C'(t)=-\int_0^t K(t-s)C(s)\dd s$、$C(0)=k_{\mathrm B}T/m$ を満たす。片側 Laplace 変換では $\widetilde C(z)=k_{\mathrm B}T/(mz+\widetilde K(z))$ である。分母を $z+b$ 倍した多項式 $mz^2+mbz+\gamma b$ の根は左半平面にあるため、相関は可積分である。従って拡散係数は $\widetilde C(0)=k_{\mathrm B}T/\gamma$ となる。

:::example[記憶が速度相関の符号を変える]
$m=b=\gamma=k_{\mathrm B}T=1$ とする。$\widetilde C(z)=(z+1)/(z^2+z+1)$ の分母を $(z+1/2)^2+3/4$ と書くと、$C(t)=e^{-t/2}[\cos(\sqrt3 t/2)+\sin(\sqrt3 t/2)/\sqrt3]$（$t\ge0$）を得る。$t=2\pi/\sqrt3$ では相関は $-e^{-\pi/\sqrt3}$ になる。瞬時摩擦の相関は $e^{-t}$ であるが、どちらも相関の積分は $1$ である。
:::

## 過減衰極限

### 慣性を消去する手順

速度が緩和した後の位置の力学を得るため、$\gamma,T>0$ を固定して $m\to0$ とする。この操作を :term[過減衰極限]{en="overdamped limit" yomi="かげんすいきょくげん"} という。自由粒子または決定論的な時間依存力を考え、初速度の二次モーメントが $O(1/m)$ 以下とする。:eqref[eq-neq-brownian-displacement] の端点項の二乗期待値は $O(m)$ になる。実際、$\tau^2=m^2/\gamma^2$ に対し、:eqref[eq-neq-brownian-velocity] で速度分散の最大の項は $A/(m\gamma)$ である。力が有限区間で有界ならその速度への寄与は $\sup |F|/\gamma$ 以下である。このため各固定時刻の位置は平均二乗で

$$
X_t=X_0+\mu\int_0^t F(s)\,\dd s+\sqrt{2D}\,W_t,\qquad \mu=\gamma^{-1},\quad D=k_{\mathrm B}T/\gamma
$$

に収束する。位置の極限は存在するが、定常速度分散は $1/m$ に比例して発散する。そのため極限の位置を時間微分して通常の瞬間速度に戻すことはできない。

位置依存力 $F$ が大域 Lipschitz 連続で、初期位置の二次モーメントが有限であるとする。速度の解を時間積分すると、位置で力に掛かる核は $\mu(1-e^{-(t-s)/\tau})$、雑音に掛かる核は $\sqrt{2D}(1-e^{-(t-s)/\tau})$ になる。従って任意の有限な $t_*$ と $t\le t_*$ について
$$
\|X_t\|_2\le\|X_0\|_2+\tau\|V_0\|_2+\sqrt{2Dt_*}
+\mu |F(0)|t_*+\mu L\int_0^t\|X_s\|_2\dd s .
$$
ここで $L$ は Lipschitz 定数である。この積分不等式から $\sup_{s\le t_*}\|X_s\|_2$ と $\sup_{s\le t_*}\|F(X_s)\|_2$ は質量によらず有界である。速度の積分表示に戻ると、力の項のノルムは $\gamma^{-1}\sup_{s\le t_*}\|F(X_s)\|_2$ 以下、雑音の項は $\sqrt{A/(m\gamma)}$ 以下だから、$\sup_{s\le t_*}\tau\|V_s\|_2=O(\sqrt m)$ となる。

同じ初期位置と雑音を使う極限過程を $\dd\bar X=\mu F(\bar X)\dd t+\sqrt{2D}\dd W$ とする。積分した運動方程式を引くと $\|X_t-\bar X_t\|_2\le\tau\|V_0-V_t\|_2+\mu L\int_0^t\|X_s-\bar X_s\|_2\dd s$ となる。従って再び積分不等式から $\sup_{t\le t_*}\|X_t-\bar X_t\|_2=O(\sqrt m)$ を得る。

:::example[ばねのある過減衰模型]
$F=-\kappa X$ なら平均は $\mathbb E[X_t]=X_0e^{-\kappa t/\gamma}$ である。位置の分散は $(k_{\mathrm B}T/\kappa)(1-e^{-2\kappa t/\gamma})$ となる。慣性を含む平均の特性方程式 $mr^2+\gamma r+\kappa=0$ では、$m\kappa/\gamma^2\ll1$ のとき遅い根が $-\kappa/\gamma$ に近づく。
:::

:::example[観測刻みの選び方]
半径 $0.5\,\mu\mathrm m$ の水と同密度の球の質量は約 $5.2\times10^{-16}\,\mathrm{kg}$ である。先の抵抗係数を用いると $\tau\simeq5.6\times10^{-8}\,\mathrm s$ となる。$\mathrm{ms}$ 刻みの位置観測ではこの速度緩和を分解できない。短時間の流体の記憶も無視した見積もりである。
:::

### 位置依存の移動度と熱ドリフト

狭い通路や壁の近くでは、摩擦は位置によって変わる。$T$ を一定とし、無流束で Boltzmann 密度に戻る局所流束 $J=\mu Fn-k_{\mathrm B}T\mu\partial_xn$ を使う。ここで $\mu(x)>0$ は滑らかとする。Itô 方程式 $\dd X=a(X)\dd t+\sqrt{2D(X)}\,\dd W$ で試験関数の平均の微分は Itô 公式から $\int (af'+Df'')n\,\dd x$ である。部分積分すると連続式の流束は $J=an-\partial_x(Dn)$ となる。$D=k_{\mathrm B}T\mu$ とおいて二つを照合すると、
$$
a(x)=\mu(x)F(x)+D'(x)
$$
となる。$D'$ も速度の次元をもち、これは :term[熱ドリフト]{en="thermal drift" yomi="ねつどりふと"} と呼ばれる。[確率過程の基礎](/physics/nonequilibrium/stochastic-processes/) で導く Stratonovich 表記の補正は $\sqrt{2D}(\sqrt{2D})'/2=D'/2$ なので、式に書くドリフトは $\mu F+D'/2$ となる。同じ局所流束を表す二つの式を比較しているため、異なる解釈でドリフトの値をそのまま使い回すと、定常密度も変わる。

:::example[通路の半分を数える]
反射壁のある $0\le x\le L$ で $F=0$、$D(x)=D_0(1+x/L)$ とする。熱ドリフトは $D_0/L$ であり、対応する定常密度は一様なので右半分にいる確率は $1/2$ である。Itô ドリフトを零にすると、無流束条件は $(Dn)'=0$ となる。規格化した密度は $[L\log2(1+x/L)]^{-1}$ であり、右半分にいる確率は $\log(4/3)/\log2$ へ変わる。異なる定常分布として実験で区別できる。
:::

:::example[急な外力と速度の境界層]
初速度 $V_0=v_0$ の粒子に一定力 $F$ を加える。過減衰解からの平均位置の差は $\tau(v_0-F/\gamma)(1-e^{-t/\tau})$ である。$v_0=0$ では $t=\tau$ の変位が $(F\tau/\gamma)e^{-1}$ となり、慣性を無視した値 $F\tau/\gamma$ との差が同程度である。それでも固定した $t>0$ では、$m\to0$ に従い差は $\tau$ に比例して消える。切り替え直後の相対誤差と、有限時間の絶対誤差の両方を区別する必要がある。
:::

### 速度緩和を消去する前後の保存則

慣性を残した密度を $p(x,v,t)$ とする。速度について $n=\int p\dd v$、$J=\int vp\dd v$、$S=\int v^2p\dd v$ と定める。位置の移流項と速度の摩擦・力の項を積分し、無限遠で境界項が消えるとすると
$$
\partial_t n=-\partial_x J,\qquad
m\partial_t J=-m\partial_x S+Fn-\gamma J.
$$
第二式の速度拡散項が零になるのは、$v$ を重みにした積分で二回微分の一回目を部分積分し、残る $\partial_vp$ の積分も零になるためである。緩和時間が観測時間より十分短く、$p$ の速度部分が局所 Maxwell 分布に近いとき、$mS\simeq k_{\mathrm B}Tn$ と $m\partial_tJ\simeq0$ が成り立つ。これを第二式に入れて $J\simeq(Fn-k_{\mathrm B}T\partial_xn)/\gamma$ を得る。位置の確率保存を保ったまま速度を消去できる。局所 Maxwell 分布から $J$ を直接求めると零になるが、その式は高々 $S$ の先頭を計算するために使った。$J$ の先頭の非零項は力と圧力勾配のずれで決まる。

:::example[壁の近くで二つの項が釣り合う]
上向きに高さ $x$ をとり、$F=-f$、$f>0$ とする。温度一定で $n=n_0e^{-x/a}$ を代入すると $J=n(-f+k_{\mathrm B}T/a)/\gamma$ である。$a=2k_{\mathrm B}T/f$ の場合、下向きの外力の半分だけが拡散で相殺し、残りの $-fn/(2\gamma)$ が流れる。$a=k_{\mathrm B}T/f$ の場合だけ無流束になる。
:::

### 小質量の長時間極限と短時間極限

熱平衡から出発した自由粒子について、過減衰分散 $2Dt$ からの正確な相対差は
$$
\frac{2Dt-\operatorname{Var}(X_t-X_0)}{2Dt}
=\frac{1-e^{-q}}q,\qquad q=t/\tau .
$$
$q\to\infty$ では $1/q$ に減り、$q\downarrow0$ では $1-q/2+O(q^2)$ になる。有限の $t$ で $m\to0$ とした場合は $q\to\infty$ である。しかし、時刻を $t=c\,m/\gamma$ と質量に合わせて短くすると $q=c$ で一定となり、相対差は消えない。平均二乗で位置が収束することと、短時間で規格化した変位が収束することとは条件が違う。

有限差分速度は過減衰式なら $\operatorname{Var}[(X_{t+h}-X_t)/h]=2D/h$ である。慣性込みなら、この分散は $h\downarrow0$ で $D/\tau$ に留まる。$h\to0$ を先にとり $m$ を固定すれば、差分速度は平均二乗で瞬時速度へ収束する。実際、差分と $V_t$ の差は $h^{-1}\int_0^h(V_{t+s}-V_t)\dd s$ であり、Cauchy–Schwarz で分散は $h^{-1}\int_0^h 2C_V(0)(1-e^{-s/\tau})\dd s$ 以下なので零へ行く。$m\to0$ を先にとった位置過程にはこの有限な速度がない。

:::example[誤差を一割に抑える観測時間]
相対差は $(1-e^{-q})/q<1/q$ なので $q\ge10$ なら一割未満である。$q=2$ なら $0.432$ 程度、$q=10$ ではほぼ $0.1$ である。例えば $\tau=10^{-6}\,\mathrm s$ の粒子を $2\,\mu\mathrm s$ 幅で観測する場合、過減衰の分散をそのまま当てはめると拡散係数を約 $43\%$ 過小評価する。同じ粒子に対して $10\,\mu\mathrm s$ の幅を使うと過小評価は一割未満になる。
:::

### 初期層の力積と位置の初期値

$m\to0$ のとき、平衡初速度の分散は $1/m$ 程度であるが、位置への初期速度の効果は $\tau V_0$ 程度なので平均二乗で消える。これに対し、初期の運動量 $P_0=mV_0$ を質量と独立に固定する極限では、自由粒子の変位に $P_0(1-e^{-t/\tau})/\gamma$ が残る。これは :eqref[eq-neq-brownian-velocity] の初期値項を時間積分したものである。$P_0$ が有限で $t>0$ を固定すれば、この項は $P_0/\gamma$ へ収束する。従って過減衰の解は $t=0$ で与えた位置から連続に始まる解に収束せず、初期層の後の位置を初期値として選ぶ必要がある。

:::example[質量とともに初速度を上げる]
$\gamma=2\,\mathrm{kg\,s^{-1}}$、$P_0=1\,\mathrm{kg\,m\,s^{-1}}$ とする。摩擦で運動量が失われるまでに平均位置は最終的に $0.5\,\mathrm m$ 進む。粒子の質量を半分にすれば速度緩和の時間は半分になるが、全変位は同じである。$t=\tau\log2$ までにその半分の $0.25\,\mathrm m$ が生じる。これは固定した温度の平衡速度から始める実験とは異なる初期条件である。
:::

## 外力下の Brown 粒子と移動度

### 定常応答と複素移動度

一定力 $F$ に対し $\mathbb E[V_t]=F/\gamma+(\mathbb E[V_0]-F/\gamma)e^{-t/\tau}$ なので、静的移動度は $\mu=1/\gamma$ である。さらに $D=\gamma k_{\mathrm B}T/\gamma^2$ より :ref[prop-neq-einstein] を再現する。平均の直線ドリフトを引いた変位分散は、力が時刻だけの関数なら自由粒子と同じである。

力を $\Re(F_\omega e^{-i\omega t})$ とすると、減衰する過渡項を除いた速度は $\Re(\mu(\omega)F_\omega e^{-i\omega t})$ となる。微分を $-i\omega$ に置けば

$$
\mu(\omega)=\frac1{\gamma-i\omega m},\qquad
\overline{F V}=\frac{|F_\omega|^2}{2}\frac{\gamma}{\gamma^2+\omega^2m^2}.
$$

上線は周期平均である。積のうち周波数 $2\omega$ の項は平均で消え、残る係数が $\Re\mu(\omega)$ である。$\omega m/\gamma\ll1$ なら応答はほぼ同位相となり、逆の領域では応答振幅が $1/(\omega m)$ に従って減る。

:::example[半分の吸収]
$\omega=\tau^{-1}$ では $\Re\mu=1/(2\gamma)$、$\Im\mu=1/(2\gamma)$ となる。速度は力に対して $\pi/4$ 遅れる。符号の読み方は $e^{-i\omega t}$ の約束に依存し、周期平均の散逸は変わらない。
:::

:::example[共分散から摩擦を推定する]
定常な自由粒子を測り、相関が初期値の $1/2$ になる時刻を $t_{1/2}$ とすれば $\gamma=m\log 2/t_{1/2}$ である。独立に測った長時間拡散係数が $k_{\mathrm B}T/\gamma$ と一致すれば、単一の時間尺度と平衡浴を仮定した模型を検査できる。
:::

### 平衡での調和拘束を慣性込みで解く

調和力 $F=-\kappa X$（$\kappa>0$）を加えたとき、平衡位置の相関 $C_X(t)$ を求める。$t>0$ の雑音は初期位置と独立なので、粒子の運動方程式を $X_0$ と掛けて平均すると $m C_X''+\gamma C_X'+\kappa C_X=0$ である。等分配から $C_X(0)=k_{\mathrm B}T/\kappa$、位置と速度が独立で $C_X'(0)=0$ である。根 $r_\pm=(-\gamma\pm\sqrt{\gamma^2-4m\kappa})/(2m)$ が異なれば

$$
\frac{C_X(t)}{C_X(0)}=\frac{r_+e^{r_-t}-r_-e^{r_+t}}{r_+-r_-}.
$$

重根 $r=-\gamma/(2m)$ なら極限は $(1-rt)e^{rt}$ である。弱減衰で根が複素共役の場合は実の減衰振動になる。自由粒子と異なり、定常変位の分散は $2[C_X(0)-C_X(t)]$ なので $2k_{\mathrm B}T/\kappa$ で飽和する。そのため、長時間の変位分散を時間で割った比は零に収束する。

:::example[臨界減衰での観測]
$m=\kappa=1$、$\gamma=2$ とした一貫した単位では $r=-1$ である。$t=1$ で規格化相関は $2/e$ であり、同じ減衰率の単一指数 $1/e$ より大きい。初期傾きが零という条件が差を作る。
:::

### 周期ポテンシャルにおける移動度と拡散

周期 $L$ の滑らかなポテンシャル $U(x)$ に一定駆動力 $F$ を加える。過減衰で $D=\mu/\beta$ を一定とすると、円上の定常密度 $p$ は $\int_0^L p=1$、$p(x+L)=p(x)$ を満たす。傾いたポテンシャルを $V=U-Fx$ と書くと、定常流は $J=-De^{-\beta V}(e^{\beta V}p)'$ である。これを $x$ から $x+L$ まで積分し、$V(x+L)=V(x)-FL$ を使えば、
$$
p(x)=\frac{J e^{-\beta V(x)}}{D(1-e^{-\beta FL})}\int_x^{x+L}e^{\beta V(y)}\,\dd y
$$
となる。$F=0$ は極限で読む。規格化と平均速度 $v=LJ$ から、
$$
\left.\dv{v}{F}\right|_{F=0}
=\frac{\mu L^2}{\int_0^L e^{\beta U(x)}\,\dd x\int_0^L e^{-\beta U(x)}\,\dd x}.
$$
分母は Cauchy–Schwarz により $L^2$ 以上であり、一定の $U$ のときに限り等号になる。

この駆動実験とは別に、無駆動で直線上の長時間の分散を求める。補正関数 $\chi(x)$ を周期的にとり、$Y=X+\chi(X)$ のドリフトを零にする。Itô 公式により、$g=1+\chi'$ は $Dg'-D\beta U'g=0$ を満たす。周期性の必要条件 $\int_0^L g=L$ から $g=Le^{\beta U}/\int_0^Le^{\beta U}$ と決まる。$\chi'$ の周期積分は零なので、この $g$ から有界な周期関数 $\chi$ を構成できる。円上の初期位相を定常密度 $e^{-\beta U}/\int_0^Le^{-\beta U}$ から選ぶと、$Y$ の増分は二次変分の期待値 $2Dt\int_0^Lg^2p\,\dd x$ を分散にもつ。有界な補正の分散は $O(1)$、共分散は Cauchy–Schwarz で $O(\sqrt t)$ なので、$X$ と $Y$ の長時間分散の傾きは一致する。従って
$$
D_{\rm eff}=\frac{DL^2}{\int_0^L e^{\beta U}\,\dd x\int_0^L e^{-\beta U}\,\dd x}
=k_{\mathrm B}T\left.\dv{v}{F}\right|_{F=0}.
$$

:::example[弱い周期力での低下率]
$U(x)=u\cos(2\pi x/L)$ とする。$u$ が小さいとき、二つの積分を $L$ で割った値はともに $1+(\beta u)^2/4+O((\beta u)^4)$ である。一次と三次は周期平均で消え、二次の係数は $\cos^2$ の平均 $1/2$ から出る。よって $D_{\rm eff}/D=1-(\beta u)^2/2+O((\beta u)^4)$ となる。先頭の低下率が $0.02$ である場合、摂動式から $\beta|u|\simeq0.2$ と読める。
:::

### 調和井戸の周波数応答を測る

質量 $m$ と摩擦 $\gamma$ を持つ粒子を固定したばねの中に置き、$f(t)=\Re(f_\omega e^{-i\omega t})$ を加える。定常平均位置を $\Re(x_\omega e^{-i\omega t})$ として運動方程式へ入れると
$$
\chi_x(\omega)=\frac{x_\omega}{f_\omega}
=\frac1{\kappa-m\omega^2-i\gamma\omega},\qquad
\mu_{\rm trap}(\omega)=-i\omega\chi_x(\omega)
$$
である。電力に相当する平均の仕事率は力と速度の積なので $|f_\omega|^2\Re\mu_{\rm trap}/2$ である。$\omega>0$ なら
$$
\Re\mu_{\rm trap}=\frac{\gamma\omega^2}{(\kappa-m\omega^2)^2+\gamma^2\omega^2}\ge0.
$$
$\omega\to0$ では $\chi_x\to1/\kappa$ だが、$\mu_{\rm trap}\to0$ である。一定力は粒子の平均位置を移しても、井戸内で一定の平均速度を保てない。自由粒子と比較する場合、$\kappa\to0$ と $\omega\to0$ の順によって速度応答の極限が異なる。

吸収係数は $\Re\mu_{\rm trap}$ である。前式の分母と分子を $\omega^2$ で割ると
$$
\Re\mu_{\rm trap}=\frac{\gamma}{(\kappa/\omega-m\omega)^2+\gamma^2}
$$
である。従って $\sqrt{\kappa/m}$ で最大の $1/\gamma$ をとる。一方、位置の振幅 $|\chi_x|$ は分母の二乗を $\omega^2$ で微分して最小化するので、正の周波数で最大になるには $\kappa/m-\gamma^2/(2m^2)>0$ が必要である。吸収の最大と位置振幅の最大は同じ測定ではない。

:::example[吸収の最大と変位の最大]
$m=\kappa=\gamma=1$ の一貫した単位を使う。吸収を最大にする角周波数は $1$、位置振幅を最大にする角周波数は $1/\sqrt2$ である。$\omega=1$ では $\chi_x=i$、$\mu_{\rm trap}=1$ なので速度は力と同位相になる。$\omega=1/\sqrt2$ では吸収係数 $\Re\mu_{\rm trap}=2/3$ である。最大の変位を生む周波数で駆動しても、同じ力振幅での平均仕事率は吸収最大の三分の二である。
:::

### 位置スペクトルから拘束と摩擦を求める

調和拘束の平衡では熱雑音の両側スペクトルが $2\gamma k_{\mathrm B}T$ である。安定な線形系では時刻 $-\infty$ からの畳み込みで平衡位置過程が得られ、周波数ごとに $\chi_x$ が掛かる。従って
$$
S_X(\omega)=\frac{2\gamma k_{\mathrm B}T}{(\kappa-m\omega^2)^2+\gamma^2\omega^2}.
$$
この式において、零周波数値は $2\gamma k_{\mathrm B}T/\kappa^2$、高周波では $\omega^4S_X\to2\gamma k_{\mathrm B}T/m^2$ となる。等時刻分散は同じ式の全周波数積分を $2\pi$ で割ったものである。独立に平衡の分散が $k_{\mathrm B}T/\kappa$ と分かるため、零周波数と分散を使えば $\gamma/\kappa=S_X(0)/(2\operatorname{Var} X)$ となる。高周波側に測定の帯域制限があれば、先頭の $\omega^{-4}$ が見えない場合もある。

:::example[ばねの位置記録の二つの較正]
$k_{\mathrm B}T$ をエネルギーの単位として $1$ とし、$\operatorname{Var}X=1/4$、$S_X(0)=1/2$ とする。$\kappa=4$ であり、$S_X(0)/(2\operatorname{Var}X)=1$ から $\gamma=4$ となる。さらに高周波側の $\omega^4S_X$ が $2$ へ近づくなら、$2\gamma/m^2=2$ より $m=2$ である。静的な分散だけ測った場合、$\gamma$ と $m$ はどちらも定まらない。
:::

### 遅れて働く摩擦の下での静的移動度

速度を消去せず、過去の速度への因果的な応答を $m\dot V(t)=-\int_0^t K(t-s)V(s)\dd s+F(t)+\xi(t)$ とおく。ここでは $K$ は絶対可積分で、雑音の平均は $\mathbb E\xi(t)=0$ とする。平均応答の片側 Laplace 変換をとると
$$
[mz+\widetilde K(z)]\widetilde{\mathbb E V}(z)
=m\mathbb E V_0+\widetilde F(z).
$$
平均速度が一定力に対して定常極限へ収束すると仮定する。そのとき $z\widetilde{\mathbb E V}(z)$ は $z\downarrow0$ でその値に収束する。実際、定常値を引いた後の積分は十分遅い時刻で任意に小さくでき、残る有界な時間での積分には前に $z$ が掛かるから消える。従って $\Gamma=\int_0^\infty K(t)\dd t>0$ なら定常移動度は $1/\Gamma$ である。周波数ごとの摩擦は $\Gamma$ だけで決まらないため、この関係は時間依存力の追随を予測しない。

:::example[指数記憶の定常応答と過渡応答]
$K(t)=\gamma b e^{-bt}$ とし、$m=\gamma=b=1$ とする。定数力に対する移動度は $1$ である。一方 Laplace 変換の応答は $(z+1)/(z^2+z+1)$ であり、極は $-1/2\pm i\sqrt3/2$ である。そのため段差力に対する偏差は減衰振動の形をとる。局所摩擦 $\gamma=1$ の偏差は $e^{-t}$ に従うので、定常応答が一致しても過渡の形は異なる。
:::

## 演習問題

::::exercise{level=basic}
$D=2\,\mathrm{mm^2\,s^{-1}}$ の等方拡散について、長時間式で $3\,\mathrm s$ 後の二次元二乗変位を求めよ。
:::solution
各成分が $2Dt=12\,\mathrm{mm^2}$ を寄与するので、和は $24\,\mathrm{mm^2}$ である。
:::
::::
::::exercise{level=basic}
$m=2\,\mathrm{kg}$、$\gamma=4\,\mathrm{kg\,s^{-1}}$、$F=12\,\mathrm N$ とする。静止から出発して平均速度が定常速度の半分に達する時刻を求めよ。
:::solution
$\tau=1/2\,\mathrm s$、定常速度は $3\,\mathrm{m\,s^{-1}}$ である。$V$ の平均は $3(1-e^{-2t})$ であり、半分になるのは $(\log2)/2\,\mathrm s$ である。
:::
::::
::::exercise{level=basic}
初速度を固定した自由粒子の変位分散を全時刻で求めよ。
:::solution
変位の積分核は $1-e^{-(t-s)/\tau}$ なので、積分を展開して分散 $2D[t-2\tau(1-e^{-t/\tau})+\tau(1-e^{-2t/\tau})/2]$ を得る。初速度が平均だけに寄与する条件での式である。
:::
::::
::::exercise{level=basic}
Maxwell 分布から選ぶ初速度の分散が $k_{\mathrm B}T/m$ となることを示せ。
:::solution
$a=\beta m/2$ と置く。$\int e^{-av^2}\dd v=\sqrt{\pi/a}$ を $a$ で微分すると、二乗の重み付き積分は $\sqrt\pi/(2a^{3/2})$ となる。正規化積分で割れば $1/(2a)=k_{\mathrm B}T/m$ である。
:::
::::
::::exercise{level=advanced}
定常自由粒子の $\operatorname{Cov}(X_t-X_0,V_t)$ を求めよ。
:::solution
$\operatorname{Cov}=\int_0^t C_V(t-s)\dd s=C_V(0)\tau(1-e^{-t/\tau})=D(1-e^{-t/\tau})$ である。単位は $\mathrm{m^2\,s^{-1}}$ である。
:::
::::
::::exercise{level=advanced}
力の履歴 $F(s)=F_0e^{-s/\tau}$（$F_0>0$）と $V_0=0$ に対する平均速度とその最大時刻を求めよ。
:::solution
:eqref[eq-neq-brownian-velocity] の積分内の指数の積は $e^{-t/\tau}$ なので平均は $(F_0/m)t e^{-t/\tau}$ である。$F_0>0$ とすると微分の符号が $1-t/\tau$ に従い、最大は $t=\tau$ となる。
:::
::::
::::exercise{level=advanced}
位置を長さ $L$ の円周上で観測する。過減衰の一定力 $F$ で一様密度が定常になることと、一周あたりの平均通過率を示せ。
:::solution
密度 $\rho=1/L$ のとき拡散流は零で、全流束 $J=\mu F/L$ は位置によらず $\partial_t\rho=-\partial_x J=0$ である。任意の切断面を横切る正味の確率流が一周あたりの平均通過率であり、$\mu F/L$ となる。
:::
::::
::::exercise{level=advanced}
同じ雑音で二つの調和井戸の過減衰運動を解く。井戸の係数と浴が同じで初期位置だけが違う場合、位置の差を求めよ。
:::solution
差 $\Delta_t$ では雑音が相殺して $\dd\Delta_t=-(\kappa/\gamma)\Delta_t\dd t$ になる。積分因子をかけると $\Delta_t=e^{-\kappa t/\gamma}\Delta_0$ であり、初期差を失う速さが分かる。
:::
::::

## 確認テスト

::::quiz
平衡初速度を用いた自由粒子で、変位分散の短時間の主要項はどれか。

- [ ] $t$ に比例する。
- [x] $t^2$ に比例する。
- [ ] $t^3$ に比例する。

:::explanation
$C_V(s)$ は原点で連続なので、:eqref[eq-neq-brownian-green-kubo] の短時間近似で $2C_V(0)\int_0^t(t-s)\dd s=C_V(0)t^2$ となる。
:::
::::
::::quiz{answer="3" tol="1e-6"}
$\gamma=2\,\mathrm{kg\,s^{-1}}$ で一定力 $6\,\mathrm N$ を受ける粒子の定常平均速度を $\mathrm{m\,s^{-1}}$ で答えよ。

:::explanation
平均加速度が零になるため $\gamma\mathbb E[V]=F$ であり、$6/2=3$ である。
:::
::::
::::quiz
定常速度分散を平衡値にする雑音係数 $A$ を選べ。力の雑音を $\sqrt{2A}\,\dd W$ とする。

- [ ] $A=m k_{\mathrm B}T$
- [x] $A=\gamma k_{\mathrm B}T$
- [ ] $A=2\gamma k_{\mathrm B}T$

:::explanation
定常分散 $A/(m\gamma)$ を $k_{\mathrm B}T/m$ に等置して $A=\gamma k_{\mathrm B}T$ を得る。
:::
::::
::::quiz{answer="2" tol="1e-6"}
$D=1\,\mathrm{m^2\,s^{-1}}$ とする。一方向の拡散の長時間での分散の傾きを $\mathrm{m^2\,s^{-1}}$ で答えよ。

:::explanation
$\operatorname{Var}(X_t-X_0)=2Dt+O(1)$ であり、傾きの極限は $2D=2$ である。
:::
::::

---
title: 連続体の運動学
description: 物質点と空間点を区別し、流れの幾何、有限・微小ひずみ、速度勾配と質量保存を導く。
notation:
  - symbol: '\vb X,\vb x'
    meaning: 物質点の参照位置と現在位置
  - symbol: '\chi_t'
    meaning: 時刻 t の運動写像
  - symbol: 'F,J'
    meaning: 変形勾配とその行列式
  - symbol: 'E,\eps'
    meaning: Green–Lagrange ひずみと微小ひずみ
  - symbol: 'L,D,W'
    meaning: 速度勾配とその対称・反対称部分
  - symbol: 'D_t'
    meaning: 物質微分
---

## この章の目標

- 物質点を追う記述と固定点で測る記述を相互に変換する。
- 線素、面素、体積の変化を変形勾配から計算する。
- 有限変形のひずみと微小変形の一次近似を使い分ける。
- 移動する領域の輸送式から質量保存の局所形を導く。

## Lagrange 記述と Euler 記述

### 連続体と配置

水の速度を測るとき、分子の速度は衝突ごとに変わる。多数の分子を含む小領域で運動量を平均すると、測定位置と時刻に対応する速度を定められる。平均領域の長さは分子間隔より十分大きく、速度が変わる巨視的な長さより十分小さくとる。これら二つの尺度の間隔が広いとき、平均値を位置の滑らかな関数で表す近似が使える。希薄気体では平均自由行程も比較すべき微視的尺度である。

:::definition[運動と配置]{#def-cont-motion}
参照領域 $B_0\subset\R^3$ の点 $\vb X$ を :term[物質点]{en="material point" yomi="ぶっしつてん"} のラベルとする。時刻 $t$ にその点を $\vb x=\chi_t(\vb X)$ に置く。$\chi_t$ を運動写像、その像 $B_t$ を現在配置という。本章では $\chi$ は十分滑らかで、各 $\chi_t$ は向きを保つ微分同相とする。参照配置を初期配置に選んだときは $\chi_0=\id$ である。
:::

位置は $\mathrm m$、時刻は $\mathrm s$、密度は $\mathrm{kg\,m^{-3}}$ で表す。デカルト座標で同じ項に二度現れる添字には和をとる。大文字添字 $A,B$ は参照配置、小文字添字 $i,j$ は現在配置の成分を表す。両者は同じ次元でも異なる変数である。

### 場を引き戻す

温度計が水槽の一点に固定されていれば、その読みは空間点を固定した関数である。水とともに流される温度計には、同じ温度場を粒子の位置に代入する。

:::definition[二つの記述]{#def-cont-descriptions}
場 $f(\vb x,t)$ を現在位置で指定する方法を :term[Euler 記述]{en="Eulerian description" yomi="おいらーきじゅつ"}、$\widehat f(\vb X,t)=f(\chi_t(\vb X),t)$ で同じ値を指定する方法を :term[Lagrange 記述]{en="Lagrangian description" yomi="らぐらんじゅきじゅつ"} という。ここで両記述のベクトル成分は同じ固定基底で測る。
:::

Euler 速度 $\vb v$ は $\partial_t\chi_t(\vb X)=\vb v(\chi_t(\vb X),t)$ で定める。逆に時刻 $0$ のラベルを実位置にとり、この常微分方程式を $\chi_0(\vb X)=\vb X$ から解けば流れ写像を得る。空間に関して局所 Lipschitz な速度場では、解が存在する範囲で初期値に対して解が一意である。終時刻から逆に解くと逆写像を得るので、異なる物質点が同じ空間点に来ることはない。有限時間で領域の外へ達する場合は、写像の定義域をその時刻まで留まる粒子に制限する。

:::example[平行移動に乗る温度計]
$\chi_t(\vb X)=\vb X+(2\,\mathrm{m\,s^{-1}})t\vb e_1$ とし、$\widehat T=300\,\mathrm K+(3\,\mathrm{K\,m^{-1}})X_1$ とする。逆写像から $X_1=x_1-(2\,\mathrm{m\,s^{-1}})t$ なので、固定点では温度が毎秒 $6\,\mathrm K$ 下がる。物質点に固定した温度計では $X_1$ が一定なので温度は変わらない。
:::

:::example[膨張と逆写像]
$\chi_t=(e^{at}X_1,X_2,X_3)$、$a=2\,\mathrm{s^{-1}}$ とする。$\partial_t\chi_1=ae^{at}X_1$ に $X_1=e^{-at}x_1$ を代入すると $v_1=ax_1$ となる。$x_1=0.5\,\mathrm m$ を通過する粒子は通過時に $1\,\mathrm{m\,s^{-1}}$ で動く。個々の粒子は $x_1$ が増えるにつれて加速する。
:::

### 粒子を識別できる時間範囲

逆写像の存在には、$J\ne 0$ による局所条件と、異なる粒子が同じ点を占めないという大域条件が要る。例えば座標を無次元化し、写像を $\chi(X_1,X_2,X_3)=(e^{X_1}\cos X_2,e^{X_1}\sin X_2,X_3)$ とする。Jacobian の行列式は常に $e^{2X_1}>0$ であるが、$X_2$ と $X_2+2\pi$ は同じ点に写る。$X_2$ を幅 $2\pi$ 以下の開区間に制限すれば、この例では偏角を一価に選べる。参照配置を取り替える場合にも一対一性が必要である。

運動の構成を速度場から始めるときは、初期条件 $\chi_{t_0}(X)=X$ を与え、$\partial_t\chi_t(X)=\vb v(\chi_t(X),t)$ を解く。空間について局所 Lipschitz な速度場では、解が同じ点に合流したと仮定して時間を逆にたどると、同じ終点を持つ初期値問題の一意性から同じ初期点に戻る。したがって軌道が存在する間は単射である。逆写像は時刻 $t$ の点を初期時刻までたどって与えられる。領域から流出する場合は、終点をそこまでたどれる粒子に限定する。

### 移流と拡大を別に測る

一次元の運動を $x=b t+(1+\alpha t)X$ とする。$b$ は基準点の速度、$\alpha$ は倍率の時間変化率であり、$1+\alpha t>0$ の範囲を扱う。$X$ で微分して隣接粒子の距離の倍率を求め、$t$ で微分して個々の粒子速度を得る。次に $X=(x-bt)/(1+\alpha t)$ を速度式に代入すると、固定点で観測する速度が分かる。測定時の座標を粒子座標と混同すると、たとえ速度の式が一次式でも加速度を誤る。

Euler 記述では $v=(b+\alpha x)/(1+\alpha t)$ である。このため $\partial_t v=-\alpha (b+\alpha x)/(1+\alpha t)^2$ と $v\partial_x v=\alpha(b+\alpha x)/(1+\alpha t)^2$ が相殺する。すなわち各粒子は等速で進むが、固定点を通過する粒子の速度は変化する。同時刻に二点で測れば速度勾配 $L=\alpha/(1+\alpha t)$ が分かり、粒子を二つ追えば倍率 $F=1+\alpha t$ が分かる。

:::example[伸長中の固定点を通過する粒子]
$b=2 \mathrm{m\,s^{-1}}$、$\alpha=\frac{1}{2} \mathrm{s^{-1}}$、$t=2 \mathrm{s}$、$x=10 \mathrm{m}$ とする。以下では$F$ が伸長倍率、$X$ が観測点に来た粒子の初期位置、$v$ が観測点の速度、$L$ が速度勾配である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
F&=\alpha t + 1=2\\
X&=\frac{- b t + x}{\alpha t + 1}=3\quad[\mathrm{m}]\\
v&=\frac{\alpha x + b}{\alpha t + 1}=\frac{7}{2}\quad[\mathrm{m\,s^{-1}}]\\
L&=\frac{\alpha}{\alpha t + 1}=\frac{1}{4}\quad[\mathrm{s^{-1}}]
\end{aligned}
$$

初期距離 $1\,\mathrm m$ の二粒子はこの時刻に $2\,\mathrm m$ 離れ、同じ場所を後から通る粒子ほど遅い。後者は固定点での時間微分が負であることに対応する。
:::

:::example[圧縮中に初期位置をたどる]
$b=1 \mathrm{m\,s^{-1}}$、$\alpha=- \frac{1}{4} \mathrm{s^{-1}}$、$t=2 \mathrm{s}$、$x=0 \mathrm{m}$ とする。以下では$F$ が伸長倍率、$X$ が観測点に来た粒子の初期位置、$v$ が観測点の速度、$L$ が速度勾配である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
F&=\alpha t + 1=\frac{1}{2}\\
X&=\frac{- b t + x}{\alpha t + 1}=-4\quad[\mathrm{m}]\\
v&=\frac{\alpha x + b}{\alpha t + 1}=2\quad[\mathrm{m\,s^{-1}}]\\
L&=\frac{\alpha}{\alpha t + 1}=- \frac{1}{2}\quad[\mathrm{s^{-1}}]
\end{aligned}
$$

観測点は原点だが、そこへ来た粒子の初期位置は負である。$1+\alpha t=0$ では全粒子が一点に集まるため、この式を使える時間は $4\,\mathrm s$ より前に限られる。
:::

:::example[速度勾配から変形率を推定する]
$b=2 \mathrm{m\,s^{-1}}$、$t=2 \mathrm{s}$、$x=10 \mathrm{m}$ に固定し、倍率の時間変化率 $\alpha$ を $\alpha=\frac{1}{2} \mathrm{s^{-1}}$ の近くで変える。速度勾配 $L$ の局所的な変化を求める。固定した値を先に代入すると、残る一変数の関数は次のようになる。以下の微分では数値を先の単位系で表す。

$$
\begin{aligned}
L(\alpha)&=\frac{\alpha}{2 \alpha + 1},\\
\frac{\mathrm dL}{\mathrm d\alpha}&=\frac{1}{\left(2 \alpha + 1\right)^{2}}.
\end{aligned}
$$

基準値で $L=\frac{1}{4}$、微分は $\frac{1}{4}$ である。したがって一次の展開は、基準値からの相対変化を $\alpha\mapsto \alpha(1+\delta)$ と書くと

$$
\frac{\Delta L}{L}
=\left(\frac{\alpha}{L}
\frac{\mathrm dL}{\mathrm d\alpha}\right)\delta+\order{\delta^2}
=\frac{1}{2}\,\delta+\order{\delta^2}
$$

となる。例えば入力を $1\%$ 増やすとき、一次式で得る出力の変化率は $\frac{1}{2}\%$ である。元の式に直接代入した変化率は $\frac{100}{201}\%$ なので、この比較で一次化の誤差を点検できる。

逆に、出力の測定値が基準値から $1\%$ 上がったとき、適用した模型が正しく、ほかの入力が固定されていれば、入力の変化率は一次で $2\%$ と推定される。同時刻に倍率も観測すれば $\alpha=FL$ を使える。その測定を省く場合は、経過時間とモデルの仮定に推定が依存する。
:::

### 基準にする時刻を変えて相対変形を追う

ある時刻 $t_0$ に付けたラベルで $\chi_t$ を書いているとする。後の時刻 $s$ の配置を新しい基準に取ると、$s$ の位置 $\vb Y$ から $t$ の位置を返す写像は $\chi_{t,s}=\chi_t\circ\chi_s^{-1}$ である。$r$ の位置から $s$ を経由して $t$ に進める合成は、$\chi_s^{-1}\circ\chi_s$ が消えるので $\chi_{t,s}\circ\chi_{s,r}=\chi_{t,r}$ となる。連鎖律から新しい変形勾配は $\grad_Y\chi_{t,s}(\vb Y)=F_t(\vb X)F_s(\vb X)^{-1}$、$\vb X=\chi_s^{-1}(\vb Y)$ である。ここで $F_t=\grad_X\chi_t$ と $J_t=\det F_t$ を用いた。相対体積比 $J_t/J_s$ を求める際も、同じ粒子のラベル $\vb X$ で両者を評価する。

:::example[途中時刻を基準にとったせん断]
$F_t=I+\gamma t\,\vb e_1\otimes\vb e_2$、$\gamma=2\,\mathrm{s^{-1}}$ とする。$N=\vb e_1\otimes\vb e_2$ は $N^2=0$ を満たすので $F_s^{-1}=I-\gamma sN$。$s=1\,\mathrm s$ から $t=3\,\mathrm s$ の増分は $I+4N$ である。各時点で $J=1$ なので、相対変形でも体積を保つ。同じ操作でも、一般に $F_t F_s^{-1}$ の順序を逆にすることはできない。例えば $F_s=\diag(2,1,1)$、$F_t=I+N$ では $F_tF_s^{-1}=\diag(1/2,1,1)+N$ だが、逆順では最後の係数が $1/2$ になる。
:::

## 物質微分と流線・流跡線・流脈線

### 粒子に沿って測る変化率

:::proposition[物質微分]{#prop-cont-material-derivative}
滑らかなスカラー場または固定基底で表したベクトル場 $f$ について、粒子に沿う :term[物質微分]{en="material derivative" yomi="ぶっしつびぶん"} は

$$
D_tf=\partial_t f+v_j\partial_j f
\label{eq-cont-material-derivative}
$$

である。特に粒子加速度は $D_t\vb v$ である。
:::

:::proof
$\widehat f(\vb X,t)=f(\chi_t(\vb X),t)$ を $\vb X$ 固定で微分する。合成関数の微分により $\partial_t\widehat f=(\partial_tf+\partial_j f\,\partial_t\chi_j)\circ\chi_t$ となる。$\partial_t\chi_j=v_j\circ\chi_t$ を代入する。加速度については $\partial_t^2\chi_t=\partial_t(\vb v\circ\chi_t)$ を用いる。
:::

:::example[定常加速]
速度が $v_1=ax_1$、$a=2\,\mathrm{s^{-1}}$ のとき、$\partial_t v_1=0$ である。一方 $D_tv_1=a^2x_1$ なので $x_1=0.5\,\mathrm m$ での加速度は $2\,\mathrm{m\,s^{-2}}$ である。固定点で一定の値を観測しても、通過する粒子の加速は残る。
:::

### 三つの曲線を作る操作

:::definition[流れを示す曲線]{#def-cont-flow-curves}
時刻 $t_*$ を固定した $\mathrm d\vb r/\mathrm ds=\vb v(\vb r,t_*)$ の積分曲線を :term[流線]{en="streamline" yomi="りゅうせん"} という。特定の粒子 $\chi_t(\vb X)$ の時間軌道を :term[流跡線]{en="pathline" yomi="りゅうせきせん"} という。固定位置 $\vb x_s$ を放出時刻 $\tau$ に通った粒子を、観察時刻 $t_*$ に集めて得る集合を :term[流脈線]{en="streakline" yomi="りゅうみゃくせん"} という。その点は $\chi_{t_*}(\chi_\tau^{-1}(\vb x_s))$ で表せる。
:::

流線は一瞬の矢印をつなぎ、流跡線は一粒の軌跡をつなぐ。これに対して流脈線は連続注入した色素の形である。定常流では、放出位置から時刻差 $t_*-\tau$ だけ同じ自律方程式を解くため、三者は同じ積分曲線の対応する区間になる。よどみ点では定数軌道を含み、非零速度の点における向き付き曲線の議論とは区別する。

:::example[三曲線を描き分ける]
単位を $\mathrm m$ と $\mathrm s$ に固定した数値表示で $\vb v=(1,t,0)$ とする。$(0,0)$ から $t=0$ に出た粒子の軌道は $x=t,y=t^2/2$ である。$t_*=2$ に原点を通る流線は $y=2x$ である。同じ時刻に原点へ連続注入した色素では、$x=2-\tau$、$y=(4-\tau^2)/2$ なので $y=2x-x^2/2$、$0\le x\le2$ である。二つの放物線は注入口 $(0,0)$ と、最初に放出した粒子の観察時の位置 $(2,2)$ を共有する。流脈線の注入口での傾きは $2$ で、その瞬間の流線と一致する。
:::

### Euler の加速度の別表示

Levi-Civita 記号 $\epsilon_{ijk}$ は $\epsilon_{123}=1$ の完全反対称記号とする。$\vb\omega=\rot\vb v$ とおくと、$\epsilon_{ijk}\epsilon_{klm}=\delta_{il}\delta_{jm}-\delta_{im}\delta_{jl}$ より

$$
(\vb v\times\vb\omega)_i=v_j\partial_i v_j-v_j\partial_j v_i.
$$

従って $D_t\vb v=\partial_t\vb v+\grad(|\vb v|^2/2)-\vb v\times\vb\omega$ である。後に [完全流体](/physics/continuum-mechanics/ideal-fluids/) で、回転がある場合の Bernoulli 関数の空間変化に使う。

### 速度計を置く位置と移流項

時不変の場でも物質微分は消えない。定常場 $v(x)=a x$ の流跡線は、$x>0$ で両辺を $x$ で割って積分すると $x=Xe^{at}$ となる。$X=0$ は別に静止解であり、負の $X$ にも同じ指数式が成り立つ。加速度はこの解の二階微分 $a^2Xe^{at}$ である。一方、固定点の速度計で $\partial_t v$ を測ればゼロであり、粒子の加速度との差はすべて $vv_x$ に入る。

### 共に運ばれるスカラーの空間勾配

$D_t c=0$ のとき $c(\chi_t(X),t)=c_0(X)$ である。この等式を $X_A$ で微分すると $c_{,i}F_{iA}=c_{0,A}$ を得る。したがって $\grad_x c=F^{-T}\grad_X c_0$。たとえば線形流に乗る色の濃さを $c_0(X)=g_0 X$ とすると、固定点では $c=g_0xe^{-at}$ である。伸長で勾配の絶対値が小さくなるのは、同じ濃度差をもつ二粒子の間隔が増すためである。濃度を単位体積あたりの質量として扱う場合には、体積変化も数える必要がある。ここでの $c$ は粒子のもつ無次元の標識である。

:::example[粒子の加速と標識のなだらかさ]
$a=\frac{1}{2} \mathrm{s^{-1}}$、$t=2 \mathrm{s}$、$X=1 \mathrm m$、$g_{0}=\frac{1}{4} \mathrm{m^{-1}}$ とする。以下では$x$ が粒子の現在位置、$v$ が粒子の速度、$A$ が粒子の加速度、$G$ が現在の標識勾配である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
x&=X e^{a t}=e\simeq 2.7183\quad[\mathrm m]\\
v&=X a e^{a t}=\frac{e}{2}\simeq 1.3591\quad[\mathrm{m\,s^{-1}}]\\
A&=X a^{2} e^{a t}=\frac{e}{4}\simeq 0.67957\quad[\mathrm{m\,s^{-2}}]\\
G&=g_{0} e^{- a t}=\frac{1}{4 e}\simeq 0.09197\quad[\mathrm{m^{-1}}]
\end{aligned}
$$

同じ粒子の標識値は $Gx=g_0X=1/4$ である。標識は粒子とともに保存されるが、隣り合う粒子の距離は $e$ 倍になる。
:::

:::example[よどみ点へ近づく粒子]
$a=- \frac{1}{4} \mathrm{s^{-1}}$、$t=4 \mathrm{s}$、$X=2 \mathrm m$、$g_{0}=\frac{1}{2} \mathrm{m^{-1}}$ とする。以下では$x$ が粒子の現在位置、$v$ が粒子の速度、$A$ が粒子の加速度、$G$ が現在の標識勾配である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
x&=X e^{a t}=\frac{2}{e}\simeq 0.73576\quad[\mathrm m]\\
v&=X a e^{a t}=- \frac{1}{2 e}\simeq -0.18394\quad[\mathrm{m\,s^{-1}}]\\
A&=X a^{2} e^{a t}=\frac{1}{8 e}\simeq 0.045985\quad[\mathrm{m\,s^{-2}}]\\
G&=g_{0} e^{- a t}=\frac{e}{2}\simeq 1.3591\quad[\mathrm{m^{-1}}]
\end{aligned}
$$

粒子は左向きに進むが、加速度は右向きである。近づくにつれて速度の絶対値が指数的に小さくなるため、有限時間内によどみ点を通過しない。
:::

:::example[標識勾配を用いた速度勾配の推定]
$t=2 \mathrm{s}$、$X=1 \mathrm m$、$g_{0}=\frac{1}{4} \mathrm{m^{-1}}$ に固定し、定常速度の勾配 $a$ を $a=\frac{1}{2} \mathrm{s^{-1}}$ の近くで変える。現在の標識勾配 $G$ の局所的な変化を求める。固定した値を先に代入すると、残る一変数の関数は次のようになる。以下の微分では数値を先の単位系で表す。

$$
\begin{aligned}
G(a)&=\frac{e^{- 2 a}}{4},\\
\frac{\mathrm dG}{\mathrm da}&=- \frac{e^{- 2 a}}{2}.
\end{aligned}
$$

基準値で $G=\frac{1}{4 e}\simeq 0.09197$、微分は $- \frac{1}{2 e}\simeq -0.18394$ である。したがって一次の展開は、基準値からの相対変化を $a\mapsto a(1+\delta)$ と書くと

$$
\frac{\Delta G}{G}
=\left(\frac{a}{G}
\frac{\mathrm dG}{\mathrm da}\right)\delta+\order{\delta^2}
=-1\,\delta+\order{\delta^2}
$$

となる。例えば入力を $1\%$ 増やすとき、一次式で得る出力の変化率は $-1\%$ である。元の式に直接代入した変化率は $-100 + \frac{100}{e^{\frac{1}{100}}}\simeq -0.99502\%$ なので、この比較で一次化の誤差を点検できる。

逆に、出力の測定値が基準値から $1\%$ 上がったとき、適用した模型が正しく、ほかの入力が固定されていれば、入力の変化率は一次で $-1\%$ と推定される。負の相対感度は、同じ時間だけ正の伸長率を増すと標識勾配がさらに弱まることを表している。
:::

## 変形勾配と有限ひずみテンソル

### 線素と体積

:::definition[変形勾配]{#def-cont-deformation-gradient}
:term[変形勾配]{en="deformation gradient" yomi="へんけいこうばい"} $F$ を $F_{iA}=\partial\chi_i/\partial X_A$ と定め、$J=\det F>0$ とおく。$F,J$ は無次元である。
:::

Taylor 展開により、参照線素 $\mathrm d\vb X$ は現在線素 $\mathrm d\vb x=F\mathrm d\vb X$ へ移る。三本の線素が張る向き付き体積は行列式であるから、体積要素は $\mathrm dv=J\mathrm dV$ となる。互いに独立な二本 $\vb A,\vb B$ の張る面素には変換

$$
\vb n\,\mathrm da=JF^{-T}\vb N\,\mathrm dA
\label{eq-cont-nanson}
$$

が成り立つ。実際、任意の $\vb C$ について $(F\vb C)\cdot((F\vb A)\times(F\vb B))=J\vb C\cdot(\vb A\times\vb B)$ が行列式の積から従う。$F\vb C$ が任意のベクトルを動くので $(F\vb A)\times(F\vb B)=JF^{-T}(\vb A\times\vb B)$ である。この面素の式を Nanson の公式という。

:::example[面積と体積の違い]
$F=\diag(2,3,1/6)$ では $J=1$ である。参照法線が $\vb e_3$ の面で $JF^{-T}\vb e_3=6\vb e_3$ なので面積は 6 倍になる。この面に垂直な厚さが $1/6$ になり、体積を保つ。
:::

### 計量からひずみを求める

:::definition[有限ひずみ]{#def-cont-finite-strain}
右 Cauchy–Green テンソルを $C=F^TF$、左 Cauchy–Green テンソルを $B=FF^T$ とする。:term[Green–Lagrange ひずみ]{en="Green–Lagrange strain" yomi="ぐりーんらぐらんじゅひずみ"} は $E=(C-I)/2$ である。Euler–Almansi ひずみは $e=(I-B^{-1})/2$ とする。
:::

長さの二乗の差は $\mathrm d\vb x\cdot\mathrm d\vb x-\mathrm d\vb X\cdot\mathrm d\vb X=2\mathrm d\vb X\cdot E\mathrm d\vb X$ である。代わりに $\mathrm d\vb X=F^{-1}\mathrm d\vb x$ を使えば同じ差は $2\mathrm d\vb x\cdot e\mathrm d\vb x$ と書ける。二本の線素の内積にも同じ計算を適用できる。現在配置に剛体回転 $Q$ を重ねると $F$ が $QF$ に変わるが、$C$ は $F^TQ^TQF=C$ で保たれる。$\vb x$ に加える一定の平行移動は最初から $F$ に現れない。

:::theorem[右極分解]{#thm-cont-polar}
$\det F>0$ なら一意に $F=RU$ と書ける。ここで $R$ は向きを保つ直交行列、$U$ は実対称正定値行列である。
:::

:::proof
$C$ は $\vb a^TC\vb a=|F\vb a|^2>0$ を満たす。:ref[thm-linalg-real-spectral] を適用し、正の固有値の正平方根で $U=C^{1/2}$ を構成する。$R=FU^{-1}$ とすると $R^TR=U^{-1}CU^{-1}=I$ であり、$\det R>0$ である。別の分解の正定値因子 $\widetilde U$ も $\widetilde U^2=C$ を満たす。$\widetilde U$ は $C$ と可換なので $C$ の各固有空間を保つ。そこで $\widetilde U$ の正の固有値は当該固有値の正平方根に限られるから $\widetilde U=U$ である。これで $R$ も定まる。
:::

$U$ の固有値 $\lambda_i$ を主伸長という。対応する単位方向では $|\mathrm d\vb x|/|\mathrm d\vb X|=\lambda_i$ であり、$E$ の固有値は $(\lambda_i^2-1)/2$ である。また $F=(RUR^T)R$ とも書ける。$J=\lambda_1\lambda_2\lambda_3$ なので、例えば有限の等容変形で保存される伸長の関数は $\sum_i\log\lambda_i$ であり、一般に $\tr E$ は保存されない。

:::example[単純せん断]
$\chi_1=X_1+\gamma X_2$、$\chi_2=X_2$、$\chi_3=X_3$ とすると

$$
E=\frac12\begin{pmatrix}0&\gamma&0\\\gamma&\gamma^2&0\\0&0&0\end{pmatrix},
\qquad J=1.
$$

$\gamma=1$ では第 2 軸方向の線素は $\sqrt2$ 倍となり、初め直交していた第 1 軸との角度は $45^\circ$ になる。面が保たれても、長さと角度は変わる。
:::

### 単純せん断の回転と伸長を分離する

平面の単純せん断 $F=\left(\begin{smallmatrix}1&\gamma\\0&1\end{smallmatrix}\right)$ を考え、面外はそのままにする。右極分解では $R=\left(\begin{smallmatrix}c&-s\\s&c\end{smallmatrix}\right)$ を用い、$U=R^TF$ の非対角成分を一致させる。第 $(1,2)$ 成分は $c\gamma+s$、第 $(2,1)$ 成分は $-s$ なので $2s=-c\gamma$。正定値性から $c>0$ の枝を取り、

$$
R=\frac1{\sqrt{4+\gamma^2}}\begin{pmatrix}2&\gamma\\-\gamma&2\end{pmatrix},
\qquad U=\frac1{\sqrt{4+\gamma^2}}\begin{pmatrix}2&\gamma\\\gamma&2+\gamma^2\end{pmatrix}.
$$

$U$ の左上成分は正、行列式は 1 なので正定値である。この結果は単純せん断のまま格子を眺めるだけでは分かりにくい。$\gamma=2$ なら極分解の回転は時計回りに $45^\circ$ であり、$U$ の主伸長は $\sqrt2+1,\sqrt2-1$ となる。実際、跡 $2\sqrt2$、行列式 1 を用いた二次方程式の根である。従って大きい主伸長の方向が圧縮される主方向と直交し、体積は変わらない。

### 伸長を加算するひずみ

一方向だけを相次いで $\lambda_1,\lambda_2$ 倍に伸ばすと、伸長比は $\lambda_1\lambda_2$ になる。対数ひずみ $\log\lambda$ はこのように同軸の伸長を合成すると加法的になる。$U$ が対称正定値なら、各固有値の実対数をとった $\log U$ を定義できる。固有分解で $\tr(\log U)=\sum_i\log\lambda_i=\log J$ であり、任意の有限変形で体積変化を正確に表す。互いに異なる軸を持つ二つの伸長では積の中に回転も含まれるので、そのまま対数ひずみを足すことはできない。

例えば引張で $\lambda=1.5$、縮める段階で同じ方向に $\lambda=2/3$ を掛けると合成は恒等である。二段階の対数ひずみは相殺するが、各段を基準として計算した Green–Lagrange ひずみを足すと $5/8-5/18=25/72$ でゼロにならない。二つの有限ひずみを足すときは、参照配置が同じかを確かめる必要がある。

### 有限変形で面積と角度を復元する

短い辺を $\vb A,\vb B$ とする長方形を考える。線素は $F$ で写るから、現在の内積は $(F\vb A)\cdot(F\vb B)=\vb A^T C\vb B$。また、任意の $\vb d$ について三重積の行列式を展開すると

$$
(F\vb d)\cdot[(F\vb A)\times(F\vb B)]
=J\,\vb d\cdot(\vb A\times\vb B).
$$

これはあらゆる $\vb d$ について成り立つため、$F^T[(F\vb A)\times(F\vb B)]=J(\vb A\times\vb B)$ である。両辺に $F^{-T}$ を掛けると :eqref[eq-cont-nanson] を得る。$J$ が全体の体積倍率でも、どの面で測るかによって面積倍率は変わる。

次に $F=\diag(a,b,c)$ で $a,b,c>0$ を考える。座標軸は主伸長方向である。参照方向 $\vb N=(1,1,0)/\sqrt2$ に沿う線素の伸長は $\sqrt{(a^2+b^2)/2}$。同じ $\vb N$ を法線とする単位面素は $JF^{-T}\vb N=(bc,ac,0)/\sqrt2$ に移る。面積倍率はこのベクトルの長さで求まり、単位法線はベクトルを面積倍率で割ったものである。例えば参照接線 $(1,-1,0)/\sqrt2$ の像と $F\vb N$ の内積は $(a^2-b^2)/2$ なので、$a\ne b$ なら $F\vb N$ は変形後の面の法線にならない。

:::example[体積保存変形の斜面]
$a=2$、$b=1$、$c=\frac{1}{2}$ とする。以下では$J$ が体積倍率、$E_{NN}$ が参照方向に沿う Green–Lagrange ひずみ、$\ell$ が線素の伸長比、$S$ が指定面の面積倍率である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
J&=a b c=1\\
E_{NN}&=\frac{a^{2}}{4} + \frac{b^{2}}{4} - \frac{1}{2}=\frac{3}{4}\\
\ell&=\sqrt{\frac{a^{2}}{2} + \frac{b^{2}}{2}}=\frac{\sqrt{10}}{2}\simeq 1.5811\\
S&=c \sqrt{\frac{a^{2}}{2} + \frac{b^{2}}{2}}=\frac{\sqrt{10}}{4}\simeq 0.79057
\end{aligned}
$$

$J=1$ でも斜めの面は縮む。線素は伸びるので、線素倍率とその線素に垂直な面の倍率が互いに逆数になるとは限らない。
:::

:::example[厚さで体積変化を補償する変形]
$a=\frac{3}{2}$、$b=\frac{1}{2}$、$c=\frac{4}{3}$ とする。以下では$J$ が体積倍率、$E_{NN}$ が参照方向に沿う Green–Lagrange ひずみ、$\ell$ が線素の伸長比、$S$ が指定面の面積倍率である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
J&=a b c=1\\
E_{NN}&=\frac{a^{2}}{4} + \frac{b^{2}}{4} - \frac{1}{2}=\frac{1}{8}\\
\ell&=\sqrt{\frac{a^{2}}{2} + \frac{b^{2}}{2}}=\frac{\sqrt{5}}{2}\simeq 1.118\\
S&=c \sqrt{\frac{a^{2}}{2} + \frac{b^{2}}{2}}=\frac{2 \sqrt{5}}{3}\simeq 1.4907
\end{aligned}
$$

同じ体積保存でも面積は増える。第３軸も面の接線方向に含まれ、$c$ が $S$ に直接掛かるためである。
:::

:::example[主伸長比と面積の測定]
$b=1$、$c=\frac{1}{2}$ に固定し、第１軸の伸長比 $a$ を $a=2$ の近くで変える。指定面の面積倍率 $S$ の局所的な変化を求める。固定した値を先に代入すると、残る一変数の関数は次のようになる。以下の微分では数値を先の単位系で表す。

$$
\begin{aligned}
S(a)&=\frac{\sqrt{\frac{a^{2}}{2} + \frac{1}{2}}}{2},\\
\frac{\mathrm dS}{\mathrm da}&=\frac{\sqrt{2} a}{4 \sqrt{a^{2} + 1}}.
\end{aligned}
$$

基準値で $S=\frac{\sqrt{10}}{4}\simeq 0.79057$、微分は $\frac{\sqrt{10}}{10}\simeq 0.31623$ である。したがって一次の展開は、基準値からの相対変化を $a\mapsto a(1+\delta)$ と書くと

$$
\frac{\Delta S}{S}
=\left(\frac{a}{S}
\frac{\mathrm dS}{\mathrm da}\right)\delta+\order{\delta^2}
=\frac{4}{5}\,\delta+\order{\delta^2}
$$

となる。例えば入力を $1\%$ 増やすとき、一次式で得る出力の変化率は $\frac{4}{5}\%$ である。元の式に直接代入した変化率は $-100 + \frac{2 \sqrt{63505}}{5}\simeq 0.80079\%$ なので、この比較で一次化の誤差を点検できる。

逆に、出力の測定値が基準値から $1\%$ 上がったとき、適用した模型が正しく、ほかの入力が固定されていれば、入力の変化率は一次で $\frac{5}{4}\%$ と推定される。感度は $a^2/(a^2+b^2)$ であり、二つの主伸長のうち大きい方の相対変化が斜面積に強く現れる。
:::

### 変形勾配が変形から生じる条件

各点で正の行列式を持つ行列場 $F_{iA}(\vb X)$ を測ったとする。これがある運動写像の勾配となるには、各成分の混合偏微分が一致する必要がある。これはひずみだけを与える場合よりも直接的な条件である。等価性には、領域の穴の有無も関わる。

:::proposition[変形勾配の局所的な適合条件]{#prop-cont-gradient-compatibility}
単連結な開領域上の $C^1$ 行列場 $F$ が $\partial_BF_{iA}=\partial_AF_{iB}$ を満たすとする。このとき $F=\grad_X\chi$ となる $C^2$ 写像 $\chi$ があり、連結領域では並進を除いて一意である。さらに $\det F>0$ ならば $\chi$ は局所的に向きを保つ微分同相である。
:::

:::proof
各 $i$ について微分形式 $F_{iA}\,\dd X_A$ は仮定により閉じている。:ref[prop-anl-potential] の経路独立性を使い、基準点からの線積分で $\chi_i$ を定める。終点を座標方向へ短距離だけ動かすと短い線分上の積分が増すため、微分は $F_{iA}$ である。二つの解の差は全成分の勾配がゼロなので連結領域で一定である。行列式の条件を加えれば、逆関数定理によって各点の近傍で局所可逆である。
:::

:::example[行列式が一定でも適合しない場]
座標を基準長さで無次元化し、$F=\diag(1,1+X_1,(1+X_1)^{-1})$ を $X_1>-1$ で与える。$J=1$ だが $\partial_1F_{22}=1$、$\partial_2F_{21}=0$ である。第二成分の線積分を $X_1,X_2$ 平面の単位正方形に沿って反時計回りに取ると、鉛直辺の寄与は $2$ と $-1$ で和は $1$ になる。閉路上の勾配の積分ならゼロになるので、この行列場は勾配にならない。
:::

:::example[穴の周囲で積分が戻らない場]
$X_1^2+X_2^2>0$ で $r^2=X_1^2+X_2^2$ とし、$F=I+\vb e_3\otimes(-X_2/r^2,X_1/r^2,0)$ とする。$J=1$ であり、第一・第二行は一定、第三行の回転も各点でゼロである。それでも $X_1=r\cos\theta,X_2=r\sin\theta$ の閉路に沿う第三行の積分は $\int_0^{2\pi}\dd\theta=2\pi$。第三成分に偏角を加える写像は切断面を入れれば定義できるが、全領域を連続な一価の写像で覆うことはできない。
:::

### 面積変換と発散

$F$ が滑らかな運動写像から来る場合、余因子行列 $G=JF^{-T}$ の各行は参照座標で発散ゼロになる。成分表示

$$
G_{iA}=\tfrac12\epsilon_{ijk}\epsilon_{ABC}F_{jB}F_{kC}
$$

を $X_A$ で微分する。第一因子の微分 $\partial_A F_{jB}$ は $A,B$ に関して対称であり、$\epsilon_{ABC}$ との縮約が消える。第二因子の微分も同様に $A,C$ に関して対称なので $\partial_A G_{iA}=0$ である。そこで現在配置のベクトル場 $\vb a$ に対し積の微分と連鎖律を使うと

$$
\partial_A\bigl(G_{iA}(a_i\circ\chi_t)\bigr)
=G_{iA}F_{jA}(\partial_j a_i)\circ\chi_t
=J(\divg\vb a)\circ\chi_t .
$$

この :term[Piola 変換]{en="Piola transformation" yomi="ぴおらへんかん"} は、:eqref[eq-cont-nanson] で面流束を参照配置へ移す操作に対応する。例えば $\chi_t=(X_1+f(X_2),X_2,X_3)$ なら $J=1$ で、$G=I-f'(X_2)\vb e_2\otimes\vb e_1$ となる。非定数成分は $X_1$ 微分に入るので、参照発散が実際にゼロになる。

## 微小ひずみテンソルと局所回転

### 何を微小とするか

$\vb u(\vb X,t)=\chi_t(\vb X)-\vb X$ とし、$H_{iA}=\partial_A u_i$ とする。$F=I+H$ を :ref[def-cont-finite-strain] に入れると

$$
E=\frac{H+H^T}{2}+\frac{H^TH}{2}.
\label{eq-cont-strain-expansion}
$$

したがって $|H|\ll1$ の範囲で二次項を落とせる。変位の大きさそのものは基準にならない。全体を遠くへ並進しても $H$ はゼロである。

:::definition[微小ひずみと回転]{#def-cont-small-strain}
:term[微小ひずみ]{en="infinitesimal strain" yomi="びしょうひずみ"} は $\eps=(H+H^T)/2$ であり、微小回転テンソルは $\Omega=(H-H^T)/2$ である。$\Omega\vb a=\vb\theta\times\vb a$ と書くと $\vb\theta=(\rot_{\!X}\vb u)/2$ となる。
:::

最後の関係は $\Omega_{ij}a_j=(\partial_j u_i-\partial_i u_j)a_j/2$ に二つの Levi-Civita 記号の縮約を適用すると得る。$\Omega$ は $\vb a\cdot\Omega\vb a=0$ なので一次で長さを変えない。単位方向 $\vb m$ の相対伸びは $\vb m^T\eps\vb m$ を一次項とする。直交方向 $\vb m,\vb n$ の変形後の内積は $2\vb m^T\eps\vb n+\order{|H|^2}$ である。従って元の直角からの減少量は一次で $2\eps_{mn}$ である。工学せん断ひずみを用いるときの係数 2 はここから来る。

:::example[回転の打ち切り誤差]
$xy$ 平面を角度 $\theta$ だけ剛体回転すると $H=Q-I$ であり、$\eps_{11}=\eps_{22}=\cos\theta-1$ となる。$\theta=\pi/3$ では $\eps_{11}=-1/2$ を得る一方、$C=I$ なので真の有限ひずみはゼロである。二次項は $\frac12(Q-I)^T(Q-I)$ であり、ちょうど $1-\cos\theta$ を各平面対角成分に与える。
:::

:::example[傾斜した方向の伸び]
$\eps_{11}=10^{-3}$、$\eps_{22}=-10^{-3}$、$\eps_{12}=2\times10^{-3}$ とする。$(\vb e_1+\vb e_2)/\sqrt2$ 方向の相対伸びは $(\eps_{11}+\eps_{22})/2+\eps_{12}=2\times10^{-3}$ である。せん断成分も斜め方向の伸びに寄与する。
:::

### ひずみがゼロの変位

連結な開領域上の $C^2$ 変位について、微分を入れ替えると

$$
\partial_j\partial_k u_i=\partial_j\eps_{ik}+\partial_k\eps_{ij}-\partial_i\eps_{jk}
\label{eq-cont-displacement-hessian}
$$

となる。右辺にひずみの定義を代入すると $\partial_j\partial_i u_k$ と $\partial_i\partial_j u_k$ が消え、$u_j$ の項も消えて $u_i$ の二項が加わる。$\eps=0$ なら二階微分が消え、連結性から $H$ は一定の反対称行列となる。したがって $\vb u=\vb c+\Omega\vb X$ である。静力学の変位の非一意性はこの並進と微小回転に限られることを [弾性体の基礎方程式](/physics/continuum-mechanics/elasticity/) で使う。

### 曲げとせん断を変位から区別する

細長い物体の平面変位 $u_x=-\kappa xy$、$u_y=\kappa x^2/2$ を考える。$\kappa$ は中立線 $y=0$ の微小曲率である。勾配は

$$
H=\begin{pmatrix}-\kappa y&-\kappa x\\\kappa x&0\end{pmatrix}.
$$

したがって $\eps_{xx}=-\kappa y$、$\eps_{yy}=\eps_{xy}=0$ で、スピンに対応する傾きは $\kappa x$。斜めの成分が勾配にあっても、対称化によって相殺する。曲げの外側と内側で伸縮の符号が異なる。

微小近似の誤差は $E-\eps=H^TH/2$ である。この例では $E_{xx}=-\kappa y+\kappa^2(x^2+y^2)/2$、$E_{yy}=\kappa^2x^2/2$、$E_{xy}=\kappa^2xy/2$ である。$y=0$ では線形理論の伸びがゼロなので相対誤差は定義できず、絶対誤差を調べる。中立線を有限回転で長さ一定に曲げたいなら、変位を二次以上の項まで定め直す必要がある。現在の $u$ は小変位の曲げを表す一つの場である。

:::example[表側の板で曲げの誤差を調べる]
$\kappa=\frac{1}{10} \mathrm{m^{-1}}$、$x=\frac{1}{5} \mathrm m$、$y=\frac{1}{50} \mathrm m$ とする。以下では$e$ が軸方向の微小ひずみ、$\theta$ が微小回転角、$D$ が有限ひずみとの差、$E_{xx}$ が軸方向の有限ひずみである。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
e&=- \kappa y=- \frac{1}{500}\\
\theta&=\kappa x=\frac{1}{50}\\
D&=\frac{\kappa^{2} \left(x^{2} + y^{2}\right)}{2}=\frac{101}{500000}\\
E_{xx}&=\frac{\kappa^{2} \left(x^{2} + y^{2}\right)}{2} - \kappa y=- \frac{899}{500000}
\end{aligned}
$$

この点では有限ひずみの補正は正で、圧縮の絶対値を小さくする。小さい回転勾配でも、線形の圧縮ひずみと比べると補正は約１割ある。
:::

:::example[裏側で回転角が増えた場合]
$\kappa=\frac{1}{5} \mathrm{m^{-1}}$、$x=\frac{1}{2} \mathrm m$、$y=- \frac{1}{50} \mathrm m$ とする。以下では$e$ が軸方向の微小ひずみ、$\theta$ が微小回転角、$D$ が有限ひずみとの差、$E_{xx}$ が軸方向の有限ひずみである。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
e&=- \kappa y=\frac{1}{250}\\
\theta&=\kappa x=\frac{1}{10}\\
D&=\frac{\kappa^{2} \left(x^{2} + y^{2}\right)}{2}=\frac{313}{62500}\\
E_{xx}&=\frac{\kappa^{2} \left(x^{2} + y^{2}\right)}{2} - \kappa y=\frac{563}{62500}
\end{aligned}
$$

線形ひずみと有限補正はともに正である。この位置では補正が線形ひずみを上回るため、線形モデルによる伸長の推定に注意を要する。
:::

:::example[有限ひずみを曲率の指標にする]
$x=\frac{1}{5} \mathrm m$、$y=\frac{1}{50} \mathrm m$ に固定し、曲率 $\kappa$ を $\kappa=\frac{1}{10} \mathrm{m^{-1}}$ の近くで変える。軸方向の有限ひずみ $E_{xx}$ の局所的な変化を求める。固定した値を先に代入すると、残る一変数の関数は次のようになる。以下の微分では数値を先の単位系で表す。

$$
\begin{aligned}
E_{xx}(\kappa)&=\frac{101 \kappa^{2}}{5000} - \frac{\kappa}{50},\\
\frac{\mathrm dE_{xx}}{\mathrm d\kappa}&=\frac{101 \kappa - 50}{2500}.
\end{aligned}
$$

基準値で $E_{xx}=- \frac{899}{500000}$、微分は $- \frac{399}{25000}$ である。したがって一次の展開は、基準値からの相対変化を $\kappa\mapsto \kappa(1+\delta)$ と書くと

$$
\frac{\Delta E_{xx}}{E_{xx}}
=\left(\frac{\kappa}{E_{xx}}
\frac{\mathrm dE_{xx}}{\mathrm d\kappa}\right)\delta+\order{\delta^2}
=\frac{798}{899}\,\delta+\order{\delta^2}
$$

となる。例えば入力を $1\%$ 増やすとき、一次式で得る出力の変化率は $\frac{798}{899}\%$ である。元の式に直接代入した変化率は $\frac{79699}{89900}\%$ なので、この比較で一次化の誤差を点検できる。

逆に、出力の測定値が基準値から $1\%$ 上がったとき、適用した模型が正しく、ほかの入力が固定されていれば、入力の変化率は一次で $\frac{899}{798}\%$ と推定される。一次項と二次項が逆符号を持つ領域では、測点によって感度がゼロに近づく。曲率を求めるときは、反対面で測ったひずみも併用できる。
:::

## 速度勾配テンソル（変形速度と渦度）

### 瞬間の分解

:::definition[変形速度とスピン]{#def-cont-velocity-gradient}
$L_{ij}=\partial_jv_i$ とし、$D=(L+L^T)/2$ を :term[変形速度テンソル]{en="rate of deformation tensor" yomi="へんけいそくどてんそる"}、$W=(L-L^T)/2$ をスピンテンソルと呼ぶ。:term[渦度]{en="vorticity" yomi="うずど"} $\vb\omega=\rot\vb v$ について $W\vb a=\vb\omega\times\vb a/2$ である。これらの単位は $\mathrm{s^{-1}}$ である。
:::

:::proposition[物質線と体積の時間発展]{#prop-cont-gradient-evolution}
参照位置固定で $\dot F=LF$、$\dot J=J\tr L$ および $\dot E=F^TDF$ が成り立つ。右辺の空間場は現在位置 $\chi_t(\vb X)$ で評価する。
:::

:::proof
$\partial_A\partial_t\chi_i=(\partial_jv_i)F_{jA}$ より最初の式を得る。次に $F(t+h)=(I+hL+\order{h^2})F(t)$ と書く。行列式を展開すると一次項は $\det(I+hL)=1+h\tr L+\order{h^2}$ なので二番目の式を得る。$F^TF$ を積の微分で計算すれば $\dot E=(F^TL^TF+F^TLF)/2$ である。
:::

物質線素 $\vb\ell$ は $\dot{\vb\ell}=L\vb\ell$ を満たす。長さの変化率は $\mathrm d\log|\vb\ell|/\mathrm dt=\vb m^TD\vb m$、方向 $\vb m=\vb\ell/|\vb\ell|$ は $\dot{\vb m}=(I-\vb m\vb m^T)D\vb m+W\vb m$ を満たす。よって方向変化に対称部分も寄与する。渦度の半分は速度勾配の反対称部分を表すため、任意の線素の角速度と同一視するとせん断による方向変化を取り逃がす。

現在配置のテンソル $B=FF^T$ を同じ粒子に沿って微分すると、:ref[prop-cont-gradient-evolution] より $D_t B=LB+BL^T$ となる。ここで $B$ は $F$ の積を $\chi_t^{-1}$ で現在配置へ移した場である。逆行列の微分を使えば $D_t B^{-1}=-B^{-1}L-L^T B^{-1}$。従って Euler–Almansi ひずみは $D_t e=D-eL-L^Te$ を満たす。既に変形した物体では、$D$ だけを時間積分しても一般に現在配置の有限ひずみは求まらない。

:::example[せん断を伸長と回転に分ける]
$\vb v=(ay,0,0)$、$a=4\,\mathrm{s^{-1}}$ では $\omega_3=-4\,\mathrm{s^{-1}}$、$D_{12}=2\,\mathrm{s^{-1}}$ である。第 1 軸に沿う線素は $L\vb e_1=0$ なので回らない。第 2 軸に沿う線素は初め $L\vb e_2=4\vb e_1\,\mathrm{s^{-1}}$ で傾く。二方向の瞬間角速度の平均は $-2\,\mathrm{s^{-1}}$ となる。
:::

:::example[体積変化と形の変化]
$\vb v=(ax,-ay,0)$ は $\divg\vb v=0$、$\vb\omega=0$ を満たすが、$D=\diag(a,-a,0)$ である。$a=1\,\mathrm{s^{-1}}$ として $t=\log 2\,\mathrm s$ まで待つと、第 1 軸の長さは 2 倍、第 2 軸の長さは $1/2$ になる。体積と回転の情報だけでは形の変化を決められない。
:::

### 速度場を局所形から読み取る

測定点 $\vb x_0$ の近傍では $\vb v(\vb x_0+\vb r)=\vb v(\vb x_0)+L\vb r+\order{|\vb r|^2}$ と展開できる。共に移動する座標では定数項を除ける。さらに

$$
L=\tfrac13(\tr L)I+D_{\rm d}+W,\qquad
D_{\rm d}=D-\tfrac13(\tr D)I
$$

と分ければ、球状の微小流体要素に生じる等方膨張、純粋な形の変化、回転が順に分かる。$D_{\rm d}$ は回転を伴わない一次の変形を表す一方、主軸の向きは時間や位置で変化できるため、時間積分した変形を求めるには $F$ の方程式を解く必要がある。

:::example[二段階のせん断を順に加える]
第 1 段階で $x\mapsto x+2y$、第 2 段階で $y\mapsto y+x$ とする。平面の変形勾配を $F_1,F_2$ とすると、合成の勾配は $F_2F_1=\left(\begin{smallmatrix}1&2\\1&3\end{smallmatrix}\right)$ である。順序を逆にすれば $F_1F_2=\left(\begin{smallmatrix}3&2\\1&1\end{smallmatrix}\right)$ となる。両者の行列式は 1 だが、第 1 軸の伸長比は前者で $\sqrt2$、後者で $\sqrt{10}$ である。観測される最終ひずみには変形の履歴が入る。
:::

### 参照ラベルを取り替える

物質ラベルを $\vb Y=\Psi(\vb X)$ に取り替え、$\Psi$ を時刻によらず向きを保つ微分同相とする。同じ運動の表示は $\widetilde\chi_t=\chi_t\circ\Psi^{-1}$ であり、$\widetilde F=F\,\grad_{\!Y}\Psi^{-1}$、$\widetilde J=J\det\grad_{\!Y}\Psi^{-1}$ となる。このとき $\widetilde J\,\mathrm d^3Y=J\,\mathrm d^3X$ なので現在体積はラベルに依存しない。向きを反転するラベルを使う場合は体積要素に行列式の絶対値を使う。例えば $\vb Y$ が長さの単位を変えた座標であれば、$F$ の値は変わるが物質点の軌道は同じである。Green–Lagrange ひずみを新しい座標へ移すときは、現在配置から引き戻した計量と参照配置の計量を同じ座標系で比較する。任意のラベルに対して参照計量を単位行列とすると、変形とは無関係なラベル変更をひずみとして数えてしまう。

### 伸長する方向も時間とともに変わる

平面速度場 $v_x=\gamma y$、$v_y=0$ を考える。変形速度 $D$ の非対角成分は $\gamma/2$ で、反対称成分 $W_{12}=\gamma/2$ は局所回転を与える。時刻ゼロに二点が $(0,\ell_0)$ だけ離れていれば、$y$ が変わらないから現在の差は $\vb\ell=\ell_0(\gamma t,1)$ となる。長さの対数を直接微分すると

$$
\frac{\mathrm d}{\mathrm dt}\log|\vb\ell|
=\frac{\gamma^2 t}{1+\gamma^2t^2}.
$$

$D$ の最大固有値は $\abs{\gamma}/2$ だが、この線素が常にその速さで伸びるわけではない。瞬間の伸び率は、現在の単位方向 $\vb n=\vb\ell/|\vb\ell|$ を用いて $\vb n^TD\vb n$ と計算する。実際に代入すると上式と一致する。$\gamma>0$ なら $\gamma t=1$ で線素が最大主伸長速度方向を通り、長時間では流れの向きに近づく。

方向自体の速度も、$\dot{\vb\ell}=L\vb\ell$ と長さの式から $\dot{\vb n}=L\vb n-(\vb n^TD\vb n)\vb n$ である。$W$ は方向を回すが、$D$ にも $\vb n$ に垂直な成分があれば方向は変わる。場の渦度だけでは、個々の線素の角度変化を決められない。

:::example[せん断の初期に観察した線素]
$\gamma=2 \mathrm{s^{-1}}$、$t=\frac{1}{4} \mathrm s$、$\ell_{0}=\frac{1}{10} \mathrm m$ とする。以下では$L$ が現在の線素の長さ、$A$ が向きの角速度、$G$ が現在の相対伸び率である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
L&=\ell_{0} \sqrt{\gamma^{2} t^{2} + 1}=\frac{\sqrt{5}}{20}\simeq 0.1118\quad[\mathrm m]\\
A&=- \frac{\gamma}{\gamma^{2} t^{2} + 1}=- \frac{8}{5}\quad[\mathrm{s^{-1}}]\\
G&=\frac{\gamma^{2} t}{\gamma^{2} t^{2} + 1}=\frac{4}{5}\quad[\mathrm{s^{-1}}]
\end{aligned}
$$

現在の線素はまだ鉛直に近い。鉛直から右に傾くため $x$ 軸から測る角度は減る。
:::

:::example[主伸長速度方向を通り過ぎた線素]
$\gamma=2 \mathrm{s^{-1}}$、$t=1 \mathrm s$、$\ell_{0}=\frac{1}{10} \mathrm m$ とする。以下では$L$ が現在の線素の長さ、$A$ が向きの角速度、$G$ が現在の相対伸び率である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
L&=\ell_{0} \sqrt{\gamma^{2} t^{2} + 1}=\frac{\sqrt{5}}{10}\simeq 0.22361\quad[\mathrm m]\\
A&=- \frac{\gamma}{\gamma^{2} t^{2} + 1}=- \frac{2}{5}\quad[\mathrm{s^{-1}}]\\
G&=\frac{\gamma^{2} t}{\gamma^{2} t^{2} + 1}=\frac{4}{5}\quad[\mathrm{s^{-1}}]
\end{aligned}
$$

先の例と相対伸び率は同じであるが、線素はより長く、回転も遅い。伸び率が等しくても現在の形や向きは一意に定まらない。
:::

:::example[観察時刻の違いが伸び率に与える影響]
$\gamma=2 \mathrm{s^{-1}}$、$\ell_{0}=\frac{1}{10} \mathrm m$ に固定し、経過時間 $t$ を $t=\frac{1}{4} \mathrm s$ の近くで変える。現在の相対伸び率 $G$ の局所的な変化を求める。固定した値を先に代入すると、残る一変数の関数は次のようになる。以下の微分では数値を先の単位系で表す。

$$
\begin{aligned}
G(t)&=\frac{4 t}{4 t^{2} + 1},\\
\frac{\mathrm dG}{\mathrm dt}&=- \frac{4 \left(2 t - 1\right) \left(2 t + 1\right)}{\left(4 t^{2} + 1\right)^{2}}.
\end{aligned}
$$

基準値で $G=\frac{4}{5}$、微分は $\frac{48}{25}$ である。したがって一次の展開は、基準値からの相対変化を $t\mapsto t(1+\delta)$ と書くと

$$
\frac{\Delta G}{G}
=\left(\frac{t}{G}
\frac{\mathrm dG}{\mathrm dt}\right)\delta+\order{\delta^2}
=\frac{3}{5}\,\delta+\order{\delta^2}
$$

となる。例えば入力を $1\%$ 増やすとき、一次式で得る出力の変化率は $\frac{3}{5}\%$ である。元の式に直接代入した変化率は $\frac{29900}{50201}\%$ なので、この比較で一次化の誤差を点検できる。

逆に、出力の測定値が基準値から $1\%$ 上がったとき、適用した模型が正しく、ほかの入力が固定されていれば、入力の変化率は一次で $\frac{5}{3}\%$ と推定される。この一次推定は最大値の手前で行っている。最大値をまたぐと時間に対する増減が反転するため、全時間で成り立つ逆関数はない。
:::

## Reynolds の輸送定理と連続の方程式

### 物質領域の積分

密度の時間変化を測るには、領域の中身と領域の大きさの両方を追う必要がある。動く領域の積分を参照配置に戻すと、固定された積分域で微分できる。

:::theorem[Reynolds の輸送定理]{#thm-cont-transport}
$B_t$ を滑らかな速度 $\vb w$ で動く有界領域とし、$f$ はその近傍で $C^1$ とする。$B_t$ の外向き単位法線を $\vb n$ とすると

$$
\frac{\mathrm d}{\mathrm dt}\int_{B_t} f\,\mathrm dv
=\int_{B_t}\partial_t f\,\mathrm dv+
\int_{\partial B_t}f\vb w\cdot\vb n\,\mathrm da.
\label{eq-cont-transport}
$$

:::

:::proof
境界の運動を領域内へ滑らかに延長し、その速度場の写像を $\Psi_t$ とする。変数変換で積分は $\int_{B_0}(f\circ\Psi_t)\det\nabla\Psi_t\,\mathrm dV$ となる。:ref[prop-cont-gradient-evolution] を速度 $\vb w$ に適用して微分すると、現在領域の被積分関数は $\partial_t f+\vb w\cdot\grad f+f\divg\vb w$ となる。最後の二項は $\divg(f\vb w)$ であるから、:ref[thm-anl-divergence] で境界積分に変換できる。境界における法線速度だけが残る。
:::

### 保存則の三つの形

:::law[質量保存]{#law-cont-mass}
物質領域では総質量が一定であり、滑らかな場では次の三条件は同値である。

$$
\begin{aligned}
\rho(\chi_t(\vb X),t)J(\vb X,t)&=\rho_0(\vb X),\\
D_t\rho+\rho\divg\vb v&=0,\\
\partial_t\rho+\divg(\rho\vb v)&=0.
\end{aligned}
\label{eq-cont-mass}
$$

:::

任意の物質部分体の質量を参照配置へ戻すと $\int_{B_0}\rho(\chi_t,t)J\,\mathrm dV=\int_{B_0}\rho_0\,\mathrm dV$ となる。その差がどこかで非零なら連続性により同じ符号を保つ小球がとれ、積分が等しいことに矛盾する。よって一行目が従う。一行目の微分と $\dot J=J\divg\vb v$ が二行目を与え、その積の微分が三行目を与える。二行目を $\rho J$ について積分すれば一行目に戻る。

一般の移動領域では :eqref[eq-cont-transport] と三行目から

$$
\frac{\mathrm d}{\mathrm dt}\int_{B_t}\rho\,\mathrm dv
=-\int_{\partial B_t}\rho(\vb v-\vb w)\cdot\vb n\,\mathrm da.
$$

右辺は境界を横切って外に出る質量流率の符号を逆にしたものである。物質境界では相対法線速度がゼロである。密度が正なら $D_t\rho=0$ と $\divg\vb v=0$ は同値である。粒子ごとに異なる初期密度を与えれば、体積保存の流れでも空間密度に勾配が残る。

:::example[固定検査面の収支]
体積 $0.1\,\mathrm{m^3}$ の槽に、密度 $1000\,\mathrm{kg\,m^{-3}}$ の水が毎秒 $2\times10^{-3}\,\mathrm{m^3}$ 入り、毎秒 $10^{-3}\,\mathrm{m^3}$ 出るなら、槽内水の質量は毎秒 $1\,\mathrm{kg}$ 増える。体積一定の満水槽ならこの流量の組は持続しない。空きのある槽なら自由表面が上昇し、同じ収支が成り立つ。
:::

:::example[移流する密度分布]
$\vb v=(2,0,0)\,\mathrm{m\,s^{-1}}$ とし、$\rho_0(X_1)=1000+10\sin(X_1/(1\,\mathrm m))$ を $\mathrm{kg\,m^{-3}}$ で表す。$J=1$ なので $\rho=1000+10\sin((x_1-v_1t)/(1\,\mathrm m))$ である。粒子の密度は一定であるが、固定点では時間とともに変わる。
:::

### 面を横切る質量と不連続

連続の方程式は衝撃波や異なる液体の界面を含む問題でも、積分形で意味を持つ。固定の検査体積内を滑らかな界面 $S_t$ が動き、両側の密度と速度は界面まで片側から連続に延びるとする。界面の単位法線 $\vb n$ は側 $-$ から側 $+$ を向き、法線速度を $s$ とする。側 $-$ の部分の外向き法線は $\vb n$、側 $+$ では $-\vb n$ である。各部分体の質量微分に :eqref[eq-cont-transport] を用いると、界面の寄与の和は

$$
\int_{S_t}\{\rho_+(\vb v_+\cdot\vb n-s)
-\rho_-(\vb v_-\cdot\vb n-s)\}\,\mathrm da
$$

となる。外側境界からの流入だけで全質量が変わると仮定し、界面に面密度を持つ質量を置かないなら、検査体積を界面近くで任意に選べるので $\rho_-(\vb v_-\cdot\vb n-s)=\rho_+(\vb v_+\cdot\vb n-s)$ である。これは移動する面に対する法線質量流束が同じであることを表す。界面で $\rho$ が不連続なときは、連続の方程式の非保存形で $\rho\divg\vb v$ を個別に評価する代わりにこの流束を使う。

:::example[動く面を貫く圧縮流]
一次元で左から密度 $\rho_-=1\,\mathrm{kg\,m^{-3}}$、速度 $v_-=4\,\mathrm{m\,s^{-1}}$ の流れが来る。界面が右向きに $s=1\,\mathrm{m\,s^{-1}}$ で動き、右側で $\rho_+=3\,\mathrm{kg\,m^{-3}}$ とする。面に固定した系での左側の流束は $1(4-1)=3\,\mathrm{kg\,m^{-2}\,s^{-1}}$。右側でも同じ流束を運ぶには $3(v_+-1)=3$ より $v_+=2\,\mathrm{m\,s^{-1}}$ である。これは質量だけの条件であり、圧力を決めるには運動量も課す必要がある。
:::

物質界面では両側で $\vb v_\pm\cdot\vb n=s$ である。従って密度が跳んでも質量収支は自動で満たされる。また、質量保存は接線成分を制約しないため、物質面を挟んで速度の接線成分が異なる模型を立てられる。完全流体の渦層はこの例である。

### 管の断面が変わる場合の物質収支

長さ方向位置 $x$ で断面積 $A(x)$ をもつ固定管を考える。側壁は不透過とし、各断面で密度 $\rho$ と軸方向速度 $u$ をほぼ一様とみなす。短い管区間 $[x,x+\Delta x]$ に含まれる質量は $\int_x^{x+\Delta x}\rho A\,\dd x'$ である。その時間変化率を左右断面の質量流率 $\rho Au$ の差と等置し、$\Delta x$ で割って極限をとると

$$
\partial_t(\rho A)+\partial_x(\rho Au)=0.
$$

断面内の分布が無視できない場合、単位長さ当たりの質量は $\int_A\rho\,\dd A$、質量流率は $\int_A\rho u\,\dd A$ になる。密度と速度の両方に分布があるとき、各々の断面平均を掛けて断面積を乗じても、一般には正しい質量流率を得られない。

定常で発散がゼロなら体積流量 $Q=Au$ は一定である。$A=A_0(1+\alpha x)$ なら $u=Q/[A_0(1+\alpha x)]$。粒子の通過時間は速度の逆数を位置で積分するので、長さ $L$ までに $t_f=A_0(L+\alpha L^2/2)/Q$ かかる。これに $Q$ を掛けると、ちょうど同じ区間の管内容積 $\int_0^L A\,\dd x$ になる。拡大・縮小が緩やかで断面速度が一様に近いときの一次元モデルである。

:::example[拡大管での通過時間]
$A_{0}=\frac{1}{100} \mathrm{m^2}$、$Q=\frac{1}{50} \mathrm{m^3\,s^{-1}}$、$\alpha=\frac{1}{10} \mathrm{m^{-1}}$、$L=2 \mathrm m$ とする。以下では$A_{L}$ が出口断面積、$u_{L}$ が出口速度、$V$ が管内容積、$t_{f}$ が通過時間である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
A_{L}&=A_{0} \left(L \alpha + 1\right)=\frac{3}{250}\quad[\mathrm{m^2}]\\
u_{L}&=\frac{Q}{A_{0} \left(L \alpha + 1\right)}=\frac{5}{3}\quad[\mathrm{m\,s^{-1}}]\\
V&=A_{0} \left(\frac{L^{2} \alpha}{2} + L\right)=\frac{11}{500}\quad[\mathrm{m^3}]\\
t_{f}&=\frac{A_{0} \left(\frac{L^{2} \alpha}{2} + L\right)}{Q}=\frac{11}{10}\quad[\mathrm s]
\end{aligned}
$$

入口速度は $2\,\mathrm{m\,s^{-1}}$ である。出口に向かって速度が下がるので、通過時間は入口速度で全長を割った $1\,\mathrm s$ より長い。
:::

:::example[縮小管での通過時間]
$A_{0}=\frac{1}{100} \mathrm{m^2}$、$Q=\frac{1}{50} \mathrm{m^3\,s^{-1}}$、$\alpha=- \frac{1}{20} \mathrm{m^{-1}}$、$L=4 \mathrm m$ とする。以下では$A_{L}$ が出口断面積、$u_{L}$ が出口速度、$V$ が管内容積、$t_{f}$ が通過時間である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
A_{L}&=A_{0} \left(L \alpha + 1\right)=\frac{1}{125}\quad[\mathrm{m^2}]\\
u_{L}&=\frac{Q}{A_{0} \left(L \alpha + 1\right)}=\frac{5}{2}\quad[\mathrm{m\,s^{-1}}]\\
V&=A_{0} \left(\frac{L^{2} \alpha}{2} + L\right)=\frac{9}{250}\quad[\mathrm{m^3}]\\
t_{f}&=\frac{A_{0} \left(\frac{L^{2} \alpha}{2} + L\right)}{Q}=\frac{9}{5}\quad[\mathrm s]
\end{aligned}
$$

出口速度は入口より速い。そのため通過時間は、管長を入口速度で割った $2\,\mathrm s$ より短い。
:::

:::example[管区間を延ばすときの通過時間]
$A_{0}=\frac{1}{100} \mathrm{m^2}$、$Q=\frac{1}{50} \mathrm{m^3\,s^{-1}}$、$\alpha=\frac{1}{10} \mathrm{m^{-1}}$ に固定し、管区間の長さ $L$ を $L=2 \mathrm m$ の近くで変える。通過時間 $t_{f}$ の局所的な変化を求める。固定した値を先に代入すると、残る一変数の関数は次のようになる。以下の微分では数値を先の単位系で表す。

$$
\begin{aligned}
t_{f}(L)&=\frac{L^{2}}{40} + \frac{L}{2},\\
\frac{\mathrm dt_{f}}{\mathrm dL}&=\frac{L + 10}{20}.
\end{aligned}
$$

基準値で $t_{f}=\frac{11}{10}$、微分は $\frac{3}{5}$ である。したがって一次の展開は、基準値からの相対変化を $L\mapsto L(1+\delta)$ と書くと

$$
\frac{\Delta t_{f}}{t_{f}}
=\left(\frac{L}{t_{f}}
\frac{\mathrm dt_{f}}{\mathrm dL}\right)\delta+\order{\delta^2}
=\frac{12}{11}\,\delta+\order{\delta^2}
$$

となる。例えば入力を $1\%$ 増やすとき、一次式で得る出力の変化率は $\frac{12}{11}\%$ である。元の式に直接代入した変化率は $\frac{1201}{1100}\%$ なので、この比較で一次化の誤差を点検できる。

逆に、出力の測定値が基準値から $1\%$ 上がったとき、適用した模型が正しく、ほかの入力が固定されていれば、入力の変化率は一次で $\frac{11}{12}\%$ と推定される。絶対的な増分の係数は出口速度の逆数に等しい。時間変化を解釈するときは、管長を変えた際に入口流量も変化していないかを確認する。
:::

## 演習問題

::::exercise{level=basic #exr-cont-inverse}
$\chi_t=(X_1+atX_2,X_2,X_3)$ の逆写像と Euler 速度を求めよ。$a$ は一定で単位は $\mathrm{s^{-1}}$ である。

:::solution
第 2、第 3 成分から $X_2=x_2,X_3=x_3$、第 1 成分から $X_1=x_1-atx_2$ である。$\partial_t\chi_t=(aX_2,0,0)$ に代入して $\vb v=(ax_2,0,0)$ を得る。
:::
::::

::::exercise{level=basic #exr-cont-material-temperature}
SI 数値で $\vb v=(2,0,0)$、$T=300+4x_1-3t$ とする。物質微分を求めよ。

:::solution
$\partial_tT=-3$、$\partial_1T=4$ である。:eqref[eq-cont-material-derivative] より $D_tT=-3+2\cdot4=5\,\mathrm{K\,s^{-1}}$ である。
:::
::::

::::exercise{level=basic #exr-cont-streak}
$\vb v=(2,t,0)$ の下で、原点に $t=0$ から連続注入した色素の $t=3$ における形を求めよ。座標の数値は SI で表す。

:::solution
放出時刻 $\tau$ の粒子は $x=2(3-\tau)$、$y=(9-\tau^2)/2$ にある。$\tau=3-x/2$ を消去すると $y=3x/2-x^2/8$ となる。$\tau$ は $0\le\tau\le3$ なので $0\le x\le6$ である。
:::
::::

::::exercise{level=basic #exr-cont-strain-diagonal}
$F=\diag(3/2,2/3,1)$ の $J,E,e$ を求めよ。

:::solution
積より $J=1$ である。$F^TF=\diag(9/4,4/9,1)$ なので $E=\diag(5/8,-5/18,0)$ である。$B^{-1}=\diag(4/9,9/4,1)$ を代入して $e=\diag(5/18,-5/8,0)$ と求まる。
:::
::::

::::exercise{level=basic #exr-cont-rate}
$\vb v=(-3y,3x,2z)$ を SI 数値で表す。$D,\vb\omega$ と瞬間の相対体積増加率を求めよ。

:::solution
$L_{12}=-3,L_{21}=3,L_{33}=2$ で、それ以外はゼロである。$D=\diag(0,0,2)\,\mathrm{s^{-1}}$、$\vb\omega=(0,0,6)\,\mathrm{s^{-1}}$ である。:ref[prop-cont-gradient-evolution] より $\dot J/J=\tr D=2\,\mathrm{s^{-1}}$ となる。
:::
::::

::::exercise{level=advanced #exr-cont-face}
$F=\left(\begin{smallmatrix}1&2&0\\0&1&0\\0&0&1\end{smallmatrix}\right)$ で、$\vb N=\vb e_1$ の面の面積比と変形後の単位法線を求めよ。

:::solution
$J=1$、$F^{-T}\vb e_1=(1,-2,0)^T$ である。:eqref[eq-cont-nanson] より面積比は $\sqrt5$、単位法線は $(1,-2,0)^T/\sqrt5$ となる。実際、接線 $(2,1,0)^T,(0,0,1)^T$ との内積はともにゼロである。
:::
::::

::::exercise{level=advanced #exr-cont-density-affine}
$\vb v=(ax,by,0)$、$\rho_0=\rho_*(1+X_1^2/\ell^2)$ とする。密度を求め、連続の方程式を直接確かめよ。

:::solution
流れ写像は $\chi_t=(e^{at}X_1,e^{bt}X_2,X_3)$、$J=e^{(a+b)t}$ である。従って $\rho=\rho_*e^{-(a+b)t}(1+e^{-2at}x^2/\ell^2)$。括弧の時刻微分は $-2ae^{-2at}x^2/\ell^2$、移流微分はその逆符号なので $D_t\rho=-(a+b)\rho$。発散が $a+b$ のため :eqref[eq-cont-mass] を満たす。
:::
::::

::::exercise{level=advanced #exr-cont-objective-length}
二つの変形勾配 $F_1,F_2$ が正の行列式を持ち、すべての線素の変形後の長さが等しいとする。$F_2=QF_1$ となる固有回転 $Q$ を求めよ。

:::solution
差 $S=F_2^TF_2-F_1^TF_1$ は対称であり、仮定は $\vb a^TS\vb a=0$ を意味する。$\vb a=\vb e_i$ と $\vb e_i+\vb e_j$ を選び、対角成分と非対角成分がすべてゼロと分かる。$Q=F_2F_1^{-1}$ なら $Q^TQ=F_1^{-T}F_2^TF_2F_1^{-1}=I$、$\det Q>0$ である。$F_1$ が可逆なので $Q$ は一意である。
:::
::::

## 確認テスト

::::quiz
固定点で測る速度が時間に依存しないとき、必ず成り立つ式を選べ。

- [x] $\partial_t\vb v=0$
- [ ] $D_t\vb v=0$
- [ ] $L=0$

:::explanation
定常性は固定点での時間微分が消えることを指す。例えば $\vb v=(ax,0,0)$ は定常で $L_{11}=a$、$D_t v_1=a^2x$ であるから、後二者は一般には成り立たない。
:::
::::

::::quiz{answer="6" tol="1e-6"}
$F=\diag(2,3,1)$ の $J$ を求めよ。

:::explanation
対角行列の行列式は対角成分の積なので $J=2\cdot3\cdot1=6$ である。
:::
::::

::::quiz{answer="-5" tol="1e-6"}
$\vb v=(5y,0,0)$ を SI 数値で表す。渦度の第 3 成分を $\mathrm{s^{-1}}$ で答えよ。

:::explanation
$\omega_3=\partial_xv_2-\partial_yv_1=0-5=-5\,\mathrm{s^{-1}}$ である。この負号は $xy$ 平面の時計回りに対応する。
:::
::::

::::quiz
体積保存の速度場で、質量保存の結果として成り立つものを選べ。

- [x] 物質点に沿う密度が一定である。
- [ ] 異なるすべての物質点で密度が等しい。
- [ ] 固定点での密度が必ず時間一定である。

:::explanation
:ref[law-cont-mass] で $\divg\vb v=0$ とすると $D_t\rho=0$ となる。初期に異なる密度の粒子があれば、その違いは運ばれる。したがって異なる粒子の通過する固定点では時間変化がありうる。
:::
::::

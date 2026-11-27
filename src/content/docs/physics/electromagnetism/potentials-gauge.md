---
title: ポテンシャルとゲージ変換
description: 電場と磁束密度をポテンシャルで表し、ゲージの自由度、縦・横成分、波動方程式の因果的な解を導く。
notation:
  - symbol: '\phi'
    meaning: スカラーポテンシャル
  - symbol: '\vb{A}'
    meaning: ベクトルポテンシャル
  - symbol: '\chi'
    meaning: ゲージ関数
  - symbol: 't_r'
    meaning: 遅延時刻
  - symbol: 'G_{\rm ret}'
    meaning: 遅延 Green 関数
---

## この章の目標

- 電磁場を二つのポテンシャルから復元する。
- 場を変えずにゲージ条件を課す方程式を解く。
- Coulomb ゲージでの瞬時項と因果性を両立させる。
- 遅延 Green 関数を構成し、源の時刻微分まで含む場を求める。

## スカラーポテンシャルとベクトルポテンシャル

### 二つの拘束式を自動的に満たす

磁束には始点も終点もなく、変動する磁束は周囲に回転する電場をつくる。
二つの関係を別々に解く代わりに、まず $\nabla\cdot\vb{B}=0$ を使う。
星形領域で滑らかな発散零の場なら $\vb{B}=\nabla\times\vb{A}$ と書ける。

原点が領域の星の中心であれば、具体的に

$$
\vb{A}(\vb{r},t)=-\vb{r}\times\int_0^1s\vb{B}(s\vb{r},t)\,\dd s
$$

ととれる。
実際

$$
\nabla\times[\vb{B}(s\vb{r})\times\vb{r}]=2\vb{B}(s\vb{r})+s(\vb{r}\cdot\nabla')\vb{B}(s\vb{r})
$$

であり、$\nabla'$ は場の引数に関する微分である。

これに $s$ を掛けると $\partial_s[s^2\vb{B}(s\vb{r})]$ になる。
積分すれば元の場が戻る。

次に Faraday 則より

$$
\nabla\times(\vb{E}+\partial_t\vb{A})=0
$$

である。
単連結な領域では閉路の周回積分が零になるので、基準点からの積分で一価のスカラー $\phi$ を定められる。
従って

$$
\vb{B}=\nabla\times\vb{A},\qquad \vb{E}=-\nabla\phi-\partial_t\vb{A}.
\label{eq-em-potentials}
$$

一般領域では、閉曲面を通る磁束がすべて零であることや、閉路に沿う積分が消えることを大域的にも確認する。
局所的な微分条件だけでは、穴を回る経路や取り除いた点を囲む面を検出できない。

:::example[軸方向に一様な変動磁場]
$\vb{B}=B(t)\hat{\vb{z}}$ に対して $\vb{A}=B(t)s\hat{\vb{\varphi}}/2$、$\phi=0$ と置く。
円柱座標で $B_z=s^{-1}\partial_s(sA_\varphi)=B(t)$、$E_\varphi=-s\dot B/2$ である。

円周積分は $-\pi s^2\dot B$ となり Faraday 則と一致する。
:::

:::example[静電場を時間依存する関数で表す]
$\phi=0,\ \vb{A}=-tE_0\hat{\vb{x}}$ とすれば $\vb{B}=0$、$\vb{E}=E_0\hat{\vb{x}}$ である。
同じ場には $\phi=-E_0x,\ \vb{A}=0$ も使える。
ポテンシャルが時間とともに変わることだけでは、観測される場の時間変化は決まらない。
:::

### 境界条件を加えて縦成分と横成分を決める

滑らかで十分速く減衰する空間場 $\vb{F}$ を考える。
発散と回転を測ってから場を再構成するには、自由空間の基本解を各成分に畳み込む。
$\vb{W}=(4\pi)^{-1}\int\vb{F}(\vb{r}')/R\,\dd^3r'$ とおけば $-\nabla^2\vb{W}=\vb{F}$ であり、二重回転の恒等式から

$$
\vb{F}=-\nabla V+\nabla\times\vb{U},\qquad V=\frac1{4\pi}\int\frac{\nabla'\cdot\vb{F}}{R}\dd^3r',\qquad \vb{U}=\frac1{4\pi}\int\frac{\nabla'\times\vb{F}}{R}\dd^3r'.
\label{eq-em-helmholtz-decomposition}
$$

ここで $V=\nabla\cdot\vb{W},\ \vb{U}=\nabla\times\vb{W}$ の微分を部分積分で源へ移した。
第一項は回転が零、第二項は発散が零であり、この分解を :term[Helmholtz 分解]{en="Helmholtz decomposition" yomi="へるむほるつぶんかい"} という。
無限遠で消えない一様場や有界領域では境界項を別に残す必要がある。

:::example[二種類の空間変化をもつ Fourier 成分]
波数 $\vb{k}\ne0$ の一成分では $\nabla\to i\vb{k}$ なので $\vb{F}_L=\vb{k}(\vb{k}\cdot\vb{F})/k^2$ と $\vb{F}_T=\vb{F}-\vb{F}_L$ に分かれる。

例えば $\vb{k}=(k,0,k)$、振幅 $\vb{F}=(F,0,0)$ なら $\vb{F}_L=(F/2,0,F/2)$、$\vb{F}_T=(F/2,0,-F/2)$ である。
前者の回転と後者の発散はそれぞれ零で、二成分の内積も零になる。
:::

### 閉じた曲面の磁束と大域ポテンシャル

原点を除いた三次元領域を考える。
そこでは、任意の閉曲線を一点へ滑らかに縮められる。
それでも発散が局所的に零の場すべてに一価で滑らかなベクトルポテンシャルがあるわけではない。
数学的な反例として $\vb{B}=b\hat{\vb{r}}/r^2$ を与える。

原点の外では発散が零だが、原点を囲む球面の磁束は $4\pi b$ である。

球面全体の近傍で滑らかな一つの $\vb{A}$ が存在すると仮定し、赤道で南北の半球に分ける。
半球それぞれに Stokes の定理を使うと、全磁束は赤道を逆向きに回る二つの線積分の和で零になる。
$b\ne0$ ならこれは矛盾である。
スカラーポテンシャルには閉路での循環、ベクトルポテンシャルには閉曲面での磁束という異なる大域条件が対応する。

この反例は物理的な磁気単極子の存在を仮定するものではなく、領域についての条件を試すための場である。

### 電気 Hertz ベクトルから保存電流の場を作る

積分を何度もやり直さずにアンテナの向きや時間波形を変えたい。そのために、まず任意の既知のベクトル密度 $\vb P_s$ から保存源 $\rho=-\nabla\cdot\vb P_s,\ \vb J=\partial_t\vb P_s$ を作る。この構成は恒等的に連続の式を満たす。$\vb P_s$ は媒質の構成式で決まる分極に限らず、局在した電流源を記述するためにも使える。電気 :term[Hertz ベクトル]{en="Hertz vector" yomi="へるつべくとる"} を

$$
(\nabla^2-c^{-2}\partial_t^2)\boldsymbol\Pi=-\vb P_s/\varepsilon_0,\qquad
\phi=-\nabla\cdot\boldsymbol\Pi,\qquad \vb A=c^{-2}\partial_t\boldsymbol\Pi
\label{eq-em-electric-hertz}
$$

で定める。実際に発散と時間微分を交換すれば $\nabla\cdot\vb A+c^{-2}\partial_t\phi=0$ である。波動演算子と同じ微分を交換すれば、$\phi$ の方程式の右辺は $\nabla\cdot\vb P_s/\varepsilon_0=-\rho/\varepsilon_0$、$\vb A$ では $-\partial_t\vb P_s/(\varepsilon_0c^2)=-\mu_0\vb J$ となる。従って過去に自由波がなければ一度だけ遅延積分を行い、

$$
\boldsymbol\Pi(\vb r,t)=\frac{1}{4\pi\varepsilon_0}
\int\frac{\vb P_s(\vb r',t-R/c)}{R}\dd^3r'
$$

から全成分が得られる。$\vb E=\nabla(\nabla\cdot\boldsymbol\Pi)-c^{-2}\partial_t^2\boldsymbol\Pi$、$\vb B=c^{-2}\nabla\times\partial_t\boldsymbol\Pi$ である。真空の源の外では最初の式は $\nabla\times\nabla\times\boldsymbol\Pi$ と等しいが、源の内部では $-\vb P_s/\varepsilon_0$ の局所項が加わる。源に重なる場を計算するときには、この区別が必要となる。


:::example[有限長ソレノイドの中心付近で一様磁場を重ねる]
局所的に一様な $B_z(t)$ に対する $A_\varphi=sB_z(t)/2$ は、回転をとると $s^{-1}\partial_s(sA_\varphi)=B_z$ となる。半径 $s$ の円周に生じる起電力は $-\pi s^2\dot B_z$ である。回転対称な配置で単位長さあたりの仕事を得るには円周長で割る。

$x=s/a$、$a=0.100\,\mathrm m$、$\dot B_z=2.00\,\mathrm{T/s}$ とし、$Y=E_\varphi$、$Y_*=a\dot B_z/2=0.100\,\mathrm{V/m}$ とする。

$$
f(x)=- x,\qquad Y=Y_* f(x).
$$

変数を少し動かしたときの読みの変化は、接線の傾きから計算できる。基準値 $x_0$ の近くで

$$
\begin{aligned}
\frac{\Delta Y}{Y_*}= f\!\left(x_0+\Delta x\right)-f(x_0)= f'(x_0)\Delta x+O((\Delta x)^2),\\
f'(x)=-1.
\end{aligned}
$$

同じ尺度で値と傾きを並べる。

| $x$ | $f(x)$ | $Y$ $[\mathrm{V/m}]$ | $\dd Y/\dd x$ $[\mathrm{V/m}]$ |
|---:|---:|---:|---:|
| $0$ | $0$ | $0$ | $-0.1$ |
| $\frac{1}{4}$ | $-0.25$ | $-0.025$ | $-0.1$ |
| $\frac{1}{2}$ | $-0.5$ | $-0.05$ | $-0.1$ |
| $1$ | $-1$ | $-0.1$ | $-0.1$ |

時間依存するベクトルポテンシャルによる場は、場の増加を打ち消す向きを向く。同じ磁場が局在するソレノイドの外では囲む磁束は一定となるため、同じ一次式を遠くまで延長できない。
:::

## ゲージ変換

### 任意関数による書き換え

:::definition[ゲージ変換]{#def-em-gauge}
空間と時間について十分滑らかな一価関数 $\chi$ を使う :term[ゲージ変換]{en="gauge transformation" yomi="げーじへんかん"} は

$$
\vb{A}'=\vb{A}+\nabla\chi,\qquad \phi'=\phi-\partial_t\chi
\label{eq-em-gauge}
$$

である。
$\chi$ の単位は $\mathrm{V\,s}$ である。
:::

$\nabla\times\nabla\chi=0$ なので $\vb{B}'=\vb{B}$ である。
また $\partial_t\nabla\chi=\nabla\partial_t\chi$ なので $\vb{E}'=\vb{E}$ となる。

逆に単連結域で同じ $\vb{B}$ をもつ二組を比べれば $\vb{A}'-\vb{A}=\nabla\chi$ と書ける。
同じ $\vb{E}$ を使うと

$$
\nabla(\phi'-\phi+\partial_t\chi)=0
$$

となるので、括弧は時間だけの関数である。

$\chi$ に時間だけの関数を足してこれを吸収すれば、二組は :eqref[eq-em-gauge] で結ばれる。

:::example[一様場のゲージ]
$\vb{A}_1=B(-y,x,0)/2$ に $\chi=Bxy/2$ を使うと $\vb{A}_2=(0,Bx,0)$ となる。
$B$ が時間一定なら両者の電位は同じでよい。
$B=B(t)$ なら同時に $\phi_2=\phi_1-\dot Bxy/2$ とする必要がある。
:::

### 粒子の作用に現れる自由度

質量 $m$、電荷 $q$ の低速粒子に $L=m v^2/2+q\vb{v}\cdot\vb{A}-q\phi$ を用いる。
運動方程式を成分ごとに計算すると $\partial L/\partial v_i=mv_i+qA_i$、

$$
\partial L/\partial x_i=qv_j\partial_iA_j-q\partial_i\phi
$$

である。

時間微分の $qv_j\partial_jA_i$ を引くと

$$
m\dot v_i=q[E_i+v_j(\partial_iA_j-\partial_jA_i)]
$$

となり Lorentz 力を再現する。
変換による差は $L'-L=q(\partial_t\chi+\vb{v}\cdot\nabla\chi)=q\,\dd\chi/\dd t$ なので、端点を固定した作用の変分は変化しない。

一方、正準運動量 $\vb{p}_{\rm can}=m\vb{v}+q\vb{A}$ の値はゲージに依存する。

:::example[孔のある領域]
ソレノイド外の $A_\varphi=\Phi/(2\pi s)$ は局所的に $\nabla(\Phi\varphi/(2\pi))$ である。
しかし角度が一周すると候補の関数は $\Phi$ だけ変わる。
一価な $\chi$ による変換では $\oint\nabla\chi\cdot\dd\vb{l}=0$ なので、この周回積分を消すことはできない。

量子力学における位相の観測には波動関数の記述が必要となるため、この古典的計算の範囲を越える。
:cite[griffiths-em]
:::

### 時間成分を消すゲージ

初期時刻 $t_1$ を選び、$\chi(\vb{r},t)=\int_{t_1}^t\phi(\vb{r},s)\dd s$ とすれば $\phi'=0$ を満たす。
これは時間ゲージと呼ばれる。
このゲージでは $\vb{E}=-\partial_t\vb{A}'$ となり、Gauss 則は $-\partial_t\nabla\cdot\vb{A}'=\rho/\varepsilon_0$ という拘束式に姿を変える。

電荷がある一般の問題で時間ゲージと Coulomb ゲージを同時に課すと、その拘束と矛盾する。

:::example[時間ゲージ中の点電荷]
静止した点電荷の Coulomb 電位 $\phi=k_eq/r$、$\vb{A}=0$ から $\chi=(t-t_1)k_eq/r$ を選ぶ。
得られるポテンシャルは

$$
\vb{A}'=-(t-t_1)k_eq\hat{\vb{r}}/r^2
$$

である。

$\vb{E}$ はその時間微分から復元でき、$\vb{B}$ は回転が零なので消える。
時間依存するゲージを許すと静電場を扱うためにもベクトルポテンシャルが必要になりうる。
:::

### ゲージを固定してから使える初期条件

ポテンシャルに四つの成分があっても、ゲージを選ぶと四成分の初期値を独立には与えられない。
Lorenz 条件の左辺を

$$
g=\nabla\cdot\vb{A}+c^{-2}\dot\phi
$$

とする。
初期に $g=0$ は $\dot\phi=-c^2\nabla\cdot\vb{A}$ を要求する。

さらに $\dot g=0$ について、スカラー波動方程式 $c^{-2}\ddot\phi=\nabla^2\phi+\rho/\varepsilon_0$ を使うと

$$
\dot g=\nabla\cdot\dot{\vb{A}}+\nabla^2\phi+\rho/\varepsilon_0=-\nabla\cdot\vb{E}+\rho/\varepsilon_0.
$$

従って初期電場への Gauss 則が、Lorenz 条件を持続させるための第二の条件になっている。
例えば $\rho=0$ のとき、任意の初期 $\phi$ と $\vb{A}$ を選ぶことはできるが、時間微分も無条件に零と指定すると、$\nabla\cdot\vb{A}=0,\ \nabla^2\phi=0$ を満たさない限り物理的な解にはならない。

### ゲージ関数が一価であるという条件

無限に長いソレノイドの外部で、二つのベクトルポテンシャルの差が $k\hat{\vb e}_\varphi/s$ であるとする。この差の回転は $s>0$ で零であり、局所的には $\chi=k\varphi$ の勾配である。しかしソレノイドを一周したときの $\chi$ の変化は $2\pi k$ である。一価の滑らかな関数の勾配なら周回積分が零でなければならないため、$k\ne0$ なら全外部で許された一価ゲージ変換では結べない。軸を切り抜いた領域での磁束情報は、点ごとの回転だけからは復元できない。

時間依存するソレノイドを準静的に考えると、この周回量が外部で測定できる。外部の円周で $\oint \vb A\cdot\dd\vb l=\Phi(t)$ を選び、一価スカラー電位の勾配の周回が零であることを使うと $\oint\vb E\cdot\dd\vb l=-\dot\Phi(t)$ である。$\Phi$ が変わるとき $\vb A$ を局所的に消した表示では、同じ起電力は切断面の両側でのスカラー電位の差に現れる。全周を一つの滑らかな座標パッチとみなしてその寄与を捨てると、Faraday の法則と整合しなくなる。

:::example[開いた弧の電位差と閉路の起電力]
一様な角度方向の誘導場をもつ円周を、角幅 $\alpha$ の弧に分ける。$E_\varphi=-\dot\Phi/(2\pi s)$ なら、弧に沿って単位電荷に働く仕事は $-\alpha\dot\Phi/(2\pi)$ である。磁束が毎秒 $2.00\,\mathrm{mWb}$ 増す場合、半円で $-1.00\,\mathrm{mV}$、一周で $-2.00\,\mathrm{mV}$ となる。同じ二点を逆の半円からつなぐ場合との仕事の差を取ると一周分が得られる。このため、時間変化する磁束を囲む二本の測定リードでは、接続点が同じでも測定経路を指定する必要がある。
:::


:::example[時間ゲージで読む一様電場の有限パルス]
空間に一様な電場 $E_x=E_0\sin(\pi t/T)$ を $0<t<T$ の間だけ加える。パルスの初めに $A_x=0$、時間ゲージで $\phi=0$ とすれば、$A_x=-\int_0^t E_x(s)\dd s$ である。

$x=t/T$、$E_0=100\,\mathrm{V/m}$、$T=1.00\,\mathrm{ns}$、$Y=A_x$、$Y_*=E_0T/\pi=31.83\,\mathrm{nT\,m}$ とする。適用範囲は電場がほぼ一様とみなせる実験領域である。

$$
f(x)=\cos{\left(\pi x \right)} - 1,\qquad Y=Y_* f(x).
$$

二点の読みの差を求めるには、基準を同じに保ち、振幅の換算を最後に一度だけ行う。最初と最後の設定で

$$
\Delta Y/Y_* = f\!\left(1\right)-f\!\left(0\right)=-2.
$$

第四列は初期設定からの変化である。

| $x$ | $f(x)$ | $Y$ $[\mathrm{nT\,m}]$ | $Y-Y(x_0)$ $[\mathrm{nT\,m}]$ |
|---:|---:|---:|---:|
| $0$ | $0$ | $0$ | $0$ |
| $\frac{1}{4}$ | $-0.2929$ | $-9.323$ | $-9.323$ |
| $\frac{1}{2}$ | $-1$ | $-31.83$ | $-31.83$ |
| $\frac{3}{4}$ | $-1.707$ | $-54.34$ | $-54.34$ |
| $1$ | $-2$ | $-63.66$ | $-63.66$ |

終了後に $A_x$ の定数が残るが、電場は零である。試験電荷の力積は $q\int E_x\dd t=-q\Delta A_x$ であり、この時間ゲージと空間一様性の下では積分結果を両端の値から求められる。
:::

## Coulomb ゲージと Lorenz ゲージ

### 方程式からゲージを選ぶ

真空で :eqref[eq-em-potentials] を Gauss 則と Ampère–Maxwell 則に代入し、

$$
\mathcal L=\nabla^2-c^{-2}\partial_t^2
$$

と置くと

$$
\nabla^2\phi+\partial_t(\nabla\cdot\vb{A})=-\rho/\varepsilon_0,\qquad \mathcal L\vb{A}=-\mu_0\vb{J}+\nabla(\nabla\cdot\vb{A}+c^{-2}\partial_t\phi).
$$

どの成分を先に解くかでゲージの選択が変わる。

:term[Coulomb ゲージ]{en="Coulomb gauge" yomi="くーろんげーじ"} は $\nabla\cdot\vb{A}_C=0$ と定める。
任意の $\vb{A}$ からは $\nabla^2\chi=-\nabla\cdot\vb{A}$ を各時刻で解けばよい。
無限遠で消える $\phi_C$ はその時刻の電荷による Coulomb 積分となる。

一方、

$$
\mathcal L\vb{A}_C=-\mu_0\vb{J}_T,\qquad \vb{J}_T=\vb{J}-\varepsilon_0\nabla\partial_t\phi_C
$$

である。
$\nabla^2\phi_C=-\rho/\varepsilon_0$ と連続の式を用いれば $\nabla\cdot\vb{J}_T=0$ となる。
従って横成分の電流が $\vb{A}_C$ の源となる。
残った自由度は各時刻で調和関数の $\chi$ であり、正則性と無限遠での条件を加えると空間勾配は消える。

:term[Lorenz ゲージ]{en="Lorenz gauge" yomi="ろーれんつげーじ"} は $\nabla\cdot\vb{A}_L+c^{-2}\partial_t\phi_L=0$ と定める。
この条件の左辺は変換で $\mathcal L\chi$ だけ変わるので、元の左辺の負を源として波動方程式を解けば条件を課せる。

結果は

$$
\mathcal L\phi_L=-\rho/\varepsilon_0,\qquad \mathcal L\vb{A}_L=-\mu_0\vb{J}.
\label{eq-em-lorenz-wave}
$$

$\mathcal L\chi=0$ を満たす波は残余自由度である。
連続の式から $G=\nabla\cdot\vb{A}_L+c^{-2}\partial_t\phi_L$ も $\mathcal LG=0$ を満たす。
初期に $G=\partial_tG=0$ なら条件は保たれる。

:::example[源のない平面波]
複素振幅を $(\phi_0,\vb{A}_0)e^{i(\vb{k}\cdot\vb{r}-\omega t)}$ とする。
Lorenz 条件は

$$
\vb{k}\cdot\vb{A}_0=\omega\phi_0/c^2
$$

である。
$\omega=c|\vb{k}|\ne0$ なら同じ位相の $\chi_0=-\phi_0/(i\omega)$ を選ぶと $\phi'_0=0$ となり、$\vb{k}\cdot\vb{A}'_0=0$ も従う。

物理的には独立な横偏光二成分が残る。
:::

:::example[Coulomb 電位の瞬時項]
Lorenz ゲージの遅延解から Coulomb ゲージへの変換では $\partial_t\chi=\phi_L-\phi_C$ である。
すると

$$
-\partial_t\vb{A}_C=-\partial_t\vb{A}_L-\nabla\phi_L+\nabla\phi_C
$$

となり、$\vb{E}$ に加える $-\nabla\phi_C$ が相殺される。

瞬時に決まる関数を使っても、測る電場は遅延解と同じである。
:::

### Coulomb ゲージへの具体的な変換関数

自由空間に局在する滑らかな源について、遅延解から Coulomb ゲージへ移す関数を実際に構成できる。
源が過去に消えていると仮定して

$$
\chi_C(\vb{r},t)=-\frac1{4\pi\varepsilon_0}\int\frac{\dd^3r'}R\int_{t-R/c}^{t}\rho(\vb{r}',s)\,\dd s
$$

と置く。
上端と下端を微分すると

$$
\partial_t\chi_C=\phi_L-\phi_C
$$

である。
従って変換後の電位は所望の $\phi_C$ になる。
変換後も Gauss 則が成り立つので、$\nabla^2\phi_C=-\rho/\varepsilon_0$ を代入すると

$$
\partial_t\nabla\cdot\vb{A}_C=0
$$

である。

過去の場と変換関数が消える条件を使えば全時刻で $\nabla\cdot\vb{A}_C=0$ とわかる。
電位の遅延を消した分だけ、$\vb{A}_C$ は遅延時刻から現在に至る源の履歴を含むことになる。

局在した源が永久に存在する場合は、過去に零とした証明をそのまま適用できない。
その場合にも $\partial_t\chi_C$ から電位の変換は確かめられるが、残った時間独立な発散には初期の空間ゲージを指定する必要がある。

### 場の応答にゲージの自由度が残らないことを波数空間で確かめる

$e^{i\vb{k}\cdot\vb{r}-i\omega t}$ の一成分を取り、$k\ne0$、$d_k=k^2-\omega^2/c^2\ne0$ とする。
Lorenz ゲージでは $\phi_L=\rho/(\varepsilon_0d_k)$、$\vb{A}_L=\mu_0\vb{J}/d_k$ である。

電場は二項の差になるが、縦成分では連続の式が $\vb{J}_L=\omega\rho\vb{k}/k^2$ を与えるため

$$
\vb{E}_L=-\frac{i\rho\vb{k}}{\varepsilon_0d_k}\left(1-\frac{\omega^2}{c^2k^2}\right)=-\frac{i\rho\vb{k}}{\varepsilon_0 k^2}.
$$

この式の添字 $L$ は縦成分を表す。
右辺は Coulomb ゲージで $-i\vb{k}\phi_C$ として得られる場と一致する。
一方、横成分は

$$
\vb{E}_T=i\omega\mu_0\vb{J}_T/d_k
$$

であり、伝播波の分母が残る。

縦・横に分ける射影自体が全空間の情報を使うので、一方の成分が瞬時の電荷と結びつくことから局所信号速度を読み取ることはできない。

### 強制振動する電流がつくる電場の三成分を一度に求める

源の時間因子を $e^{-i\omega t}$、波数を $k_0=\omega/c\ne0$ とする。
自由空間で $G=e^{ik_0R}/(4\pi R)$ は $(\nabla^2+k_0^2)G=-\delta^3$ を満たす。
Lorenz ゲージの遅延解に $\nabla'\cdot\vb{J}=i\omega\rho$ を代入し、電位の積分で部分積分すると

$$
\vb{E}=i\omega\mu_0\left(\vb{I}+\frac{\nabla\nabla}{k_0^2}\right)\cdot\int G(\vb{r}-\vb{r}')\vb{J}(\vb{r}')\,\dd^3r' .
$$

ここで $\vb{I}$ は単位テンソル、$\nabla\nabla\cdot\vb{f}$ は $\nabla(\nabla\cdot\vb{f})$ である。
実際

$$
\int G\nabla'\cdot\vb{J}=-\int(\nabla'G)\cdot\vb{J}=\nabla\cdot\int G\vb{J}
$$

であり、$\vb{E}=i\omega\vb{A}-\nabla\phi$ の二項がこの形にまとまる。

従ってテンソル核

$$
\mathcal G_{ij}=\delta_{ij}G+\frac1{k_0^2}\partial_i\partial_jG
$$

が電流の向きと測定する電場成分を結ぶ。

検算には元のポテンシャルを使わず、波動演算子を作用させる。
$\vb{f}=\int G\vb{J}$ とすれば

$$
(\nabla^2+k_0^2)\vb{f}=-\vb{J}
$$

で、勾配の回転は零である。
従って

$$
(\nabla\times\nabla\times-k_0^2)\left[\vb{f}+k_0^{-2}\nabla(\nabla\cdot\vb{f})\right]=-(\nabla^2+k_0^2)\vb{f}=\vb{J} .
$$

ここで偏微分は観測座標に関してとる。
源と観測点が重なる場合、$\partial_i\partial_j G$ を普通の関数として積分すると原点の寄与を失うため、分布微分として扱う。
遠方では $\partial_i\partial_jG\simeq-k_0^2n_i n_jG$ なので $\mathcal G_{ij}\simeq(\delta_{ij}-n_i n_j)G$ となり、横成分だけが放射されることが直接わかる。

### 瞬時のスカラー電位が遅延電場を与える仕組み

局在した保存源の Lorenz ゲージを $\phi_L,\vb A_L$ とし、瞬時の Coulomb 電位を $\phi_C=k_e\int\rho(\vb r',t)\dd^3r'/R$ とする。$\phi_C=\phi_L-\partial_t\chi$ から時間を積分して、十分昔に源がなかった場合に

$$
\chi(\vb r,t)=-k_e\int\frac{\dd^3r'}{R}\int_{t-R/c}^t\rho(\vb r',s)\dd s.
$$

上限と下限を微分すれば $\partial_t\chi=\phi_L-\phi_C$ が直ちに確かめられる。$R$ による空間微分では下端も変わる。$H=\int_{t-R/c}^t\rho(s)\dd s$ と置くと $\partial_RH=\rho(t-R/c)/c$ であり、

$$
\nabla\chi=k_e\int\hat{\vb R}
\left[\frac H{R^2}-\frac{\rho(t-R/c)}{cR}\right]\dd^3r'.
$$

従って $\vb A_C=\vb A_L+\nabla\chi$ は一般に現在と遅延時刻の間にわたる電荷の履歴を含む。この補正の時間微分と $-\nabla\phi_C$ の瞬時項は、電場を作るときに厳密に相殺する。さらに $\nabla^2\partial_t\chi=\nabla^2(\phi_L-\phi_C)=c^{-2}\partial_t^2\phi_L$ なので、過去の積分定数が零なら $\nabla^2\chi=c^{-2}\partial_t\phi_L=-\nabla\cdot\vb A_L$ となる。従ってこの積分が作るポテンシャルは実際に Coulomb 条件を満たす。

## 遅延ポテンシャル

### 波面上にある点源の応答

無限の三次元真空で、過去からの自由な入射波を零とする。
$\mathcal L G=-\delta^3(\vb{R})\delta(t-t')$ の因果的な基本解を構成する。
$\tau=t-t'$ とし、空間 Fourier 変換をとると、波数 $k$ ごとに

$$
(\partial_\tau^2+c^2k^2)\widetilde G=c^2\delta(\tau)
$$

である。

$\tau<0$ で零、$\tau=0$ で関数が連続、時間微分の跳びが $c^2$ であることから

$$
\widetilde G=c\sin(ck\tau)\Theta(\tau)/k
$$

を得る。
逆変換の角度積分は $4\pi\sin(kR)/(kR)$ なので動径積分に $\int_0^\infty\sin(kR)\sin(kc\tau)\dd k=\pi[\delta(R-c\tau)-\delta(R+c\tau)]/2$ が現れる。

これは余弦の Fourier 表現を差し引いた分布等式である。
従って

$$
G_{\rm ret}(\vb{R},\tau)=\frac{\delta(\tau-R/c)}{4\pi R}.
\label{eq-em-retarded-green}
$$

$\delta(\tau+R/c)/(4\pi R)$ も方程式を満たすが、将来の源から現在の場を定める先進解である。
選択は時間境界条件による。

式 :eqref[eq-em-retarded-green] の核を :eqref[eq-em-lorenz-wave] の源と時空で畳み込み、$t'$ を積分すると

$$
\phi_L(\vb{r},t)=\frac1{4\pi\varepsilon_0}\int\frac{\rho(\vb{r}',t_r)}{R}\,\dd^3r',\qquad \vb{A}_L(\vb{r},t)=\frac{\mu_0}{4\pi}\int\frac{\vb{J}(\vb{r}',t_r)}{R}\,\dd^3r',\quad t_r=t-R/c.
\label{eq-em-retarded}
$$

これを :term[遅延ポテンシャル]{en="retarded potential" yomi="ちえんぽてんしゃる"} と呼ぶ。
初期に外部波があれば対応する同次解を加える。
境界がある場合には自由空間核に境界条件を満たす補正が必要である。

### 遅延時刻も微分する

角括弧 $[f]$ を $t_r$ で評価した量、$\vb{n}=(\vb{r}-\vb{r}')/R$ とする。
場の点を動かすと $\nabla t_r=-\vb{n}/c$ なので、$1/R$ の微分だけでは足りない。
積の勾配と回転をとると

$$
\begin{aligned}
\vb{E}&=\frac1{4\pi\varepsilon_0}\int\left(\frac{[\rho]\vb{n}}{R^2}+\frac{[\dot\rho]\vb{n}}{cR}-\frac{[\dot{\vb{J}}]}{c^2R}\right)\dd^3r',\\
\vb{B}&=\frac{\mu_0}{4\pi}\int\left(\frac{[\vb{J}]}{R^2}+\frac{[\dot{\vb{J}}]}{cR}\right)\times\vb{n}\,\dd^3r'.
\end{aligned}
$$

これらを Jefimenko の式という。
源の時刻変化が零なら Coulomb 則と Biot–Savart 則に戻る。
時間尺度が $\tau$ のとき、$[f]\simeq f(t)-(R/c)\dot f(t)$ という展開は $R/(c\tau)\ll1$ で使える。

源の広がりが小さくても、遠方の観測時刻から伝播遅延全体を捨てることはできない。

:::example[電流パルスが届く時刻]
長さ $0.30\,\mathrm m$ 程度の局在源から $300\,\mathrm m$ 離れた点では、中心からの遅延が $1.00\,\mathrm{\mu s}$、源の両端からの到着差は最大で約 $1.00\,\mathrm{ns}$ である。

$10\,\mathrm{ns}$ の変動なら両端の遅延差を小さな補正にできるが、中心からの遅延は波形の到着時刻を決める。
:::

:::example[因果的な応答と正弦定常応答の対応]
局在源を過去から断熱的に立ち上げたとき、周波数 $\omega$ の応答は遅延核の時間 Fourier 変換で得られる。
実際に

$$
G_\omega(R)=\int_{-\infty}^{\infty}\frac{\delta(\tau-R/c)}{4\pi R}e^{i\omega\tau}\dd\tau=\frac{e^{ikR}}{4\pi R},\qquad k=\omega/c.
$$

時間因子を $e^{-i\omega t}$ とすると位相一定面は $R$ が大きい側へ移動する。
$\omega\to0$ とすれば静電場の基本解に戻る。
例えば $f=100\,\mathrm{MHz}$ では距離を約 $1.50\,\mathrm m$ 増やすと位相は $\pi$ 進み、同一時刻の振幅は逆向きとなる。
:::

### 誘導場と放射場の大きさを比べる

遅延のある電流から磁場を計算する際、源の外にいても $R^{-2}$ の項と $R^{-1}$ の項が同時に寄与する。
局在した電流パルスを遠くから観測するとき、変化時間を $\tau$ として両項の代表比は $R/(c\tau)$ である。
例えば $R=30\,\mathrm m,\ \tau=1.0\,\mathrm{\mu s}$ なら約 $0.10$、$R=3.0\,\mathrm{km}$ なら約 $10$ となる。

ただし源の異なる部分で電流の向きが打ち消し合えば、この尺度だけで放射項の係数は決まらない。
次章以降で使う多重極展開は、この打ち消しを源の空間積分に残す方法である。

### Green 関数に境界で反射する波を加える

半空間 $z>0$ の波動方程式のスカラー解 $u$ に $u(x,y,0,t)=0$ を課す。
$\bar{\vb{r}}'=(x',y',-z')$ として自由空間遅延核の差

$$
G_D(\vb{r},\vb{r}',\tau)=\frac{\delta(\tau-|\vb{r}-\vb{r}'|/c)}{4\pi|\vb{r}-\vb{r}'|}-\frac{\delta(\tau-|\vb{r}-\bar{\vb{r}}'|/c)}{4\pi|\vb{r}-\bar{\vb{r}}'|}
$$

を作る。
実源と反射源から板上の各点までの距離が同じなので、二つのパルスは板上で同時に打ち消し合う。
半空間内の特異点は元の源だけで、過去からの入射もない。
従って所望の因果的な核である。
法線微分を壁で零にしたければ、差を和に変えると鏡映による法線微分の反転で条件が満たされる。

:::example[二本の到着経路を区別する]
源を $z'=1.00\,\mathrm m$、観測点を同じ垂線の $z=3.00\,\mathrm m$ に置く。
短いパルスは直接の距離 $2.00\,\mathrm m$ と、板で反射した距離 $4.00\,\mathrm m$ を通る。

最初の到着は $6.67\,\mathrm{ns}$、次は $13.34\,\mathrm{ns}$ で、Dirichlet 条件なら後者の符号は逆、振幅の絶対値は半分である。
電磁場をこの方法で扱う場合には、ベクトルの各成分に同じ符号を機械的に選ばず、導体面の接線電場と法線磁束の条件を確認する。
:::

### 初期波の球面平均による表示

三次元の無源波で初期変位 $f$、初期速度 $g$ があるとする。
前節の Fourier 核を逆変換すると、$g$ に掛かる $\sin(ckt)/(ck)$ の空間核は $G_{\rm ret}/c^2$ である。
半径方向のデルタを積分すれば、$g$ の寄与は $t M_g(\vb{r},ct)$ となる。

ここで $M_g(\vb{r},a)=(4\pi)^{-1}\int g(\vb{r}+a\vb{n})\dd\Omega$ は球面平均である。
初期変位に掛かる $\cos(ckt)$ は同じ核の時間微分なので、解は

$$
u(\vb{r},t)=\partial_t[t M_f(\vb{r},ct)]+t M_g(\vb{r},ct).
$$

$t\downarrow0$ では $M_f(\vb{r},ct)=f(\vb{r})+O(t^2)$ なので所定の二初期値を再現する。
観測点の内部にある球全体ではなく、過去光円錐との交線である球面上の初期値だけが現れるのは、無限の三次元真空の波動方程式の性質である。

### 遅延解を用いた双極子層の境界表示

静電場と同じように、時間調和するスカラー場でも境界値とその法線微分から内部解を再構成できる。
滑らかな有界領域 $V$ で $(\nabla^2+k^2)u=-s$ とする。
$V$ の内部点に置いた自由空間 Helmholtz 核 $G$ と $u$ に Green の第二恒等式を適用すると

$$
u(\vb{r})=\int_V Gs\,\dd^3r'+\oint_{\partial V}(G\partial_{n'}u-u\partial_{n'}G)\,\dd S'.
$$

$k^2 G u$ の項が打ち消し合い、核のデルタから観測点の値が残る。
電位に相当する $u$ も $\partial_{n'}u$ も必要な表示であって、両方を任意に指定できるという意味ではない。
与えられた境界データの一方だけで内部解を表すには、その境界条件に適合する核を選ぶか、もう一方の境界値を積分方程式から求める。

:::example[空洞の固有周波数で境界値だけから解が決まらなくなる例]
一次元 $0<x<a$ の無源 Helmholtz 方程式に $u(0)=0,\ u(a)=U$ を与える。
$\sin ka\ne0$ なら $u=U\sin(kx)/\sin(ka)$ である。
$ka=n\pi$、$n\ge1$ では右端は振幅にかかわらず零となり、$U=0$ なら任意振幅の固有関数、$U\ne0$ なら有限振幅の正弦定常解が存在しない。

静電場の一意性をそのまま任意の実周波数に延長できない理由である。
時間領域で初期値を与えた問題はそれでも解け、厳密な共鳴駆動では振幅が時間とともに増える。
:::

### 正弦波電流から電場を直接求める Green 関数

一点で電流の向きを変えると、観測点での電場の向きも変わる。スカラー波の Green 関数だけを各電流成分に掛けると電荷の寄与を落とすので、一般の電流から電場を求める応答には二つの添字が必要になる。真空中の局在源を考え、時間因子を $e^{-i\omega t}$、$k=\omega/c>0$ とする。外向き条件に従うスカラー核は $g(R)=e^{ikR}/(4\pi R)$ であり、$(\nabla^2+k^2)g=-\delta^3$ を満たす。

電流密度の複素振幅を $\vb{J}$ とすると、連続の式より $\rho=\nabla\cdot\vb{J}/(i\omega)$ である。Lorenz ゲージのポテンシャルを $\vb{A}=\mu_0\int g\vb{J}\,\dd^3r'$、$\phi=\varepsilon_0^{-1}\int g\rho\,\dd^3r'$ と書く。スカラー電位の積分を部分積分すると、$\nabla' g=-\nabla g$ の符号と部分積分の符号が打ち消し合う。従って $\vb{E}=i\omega\vb{A}-\nabla\phi$ は

$$
E_i(\vb{r})=i\omega\mu_0\int G_{ij}(\vb{r}-\vb{r}')J_j(\vb{r}')\,\dd^3r',
\qquad G_{ij}=\delta_{ij}g+\frac{1}{k^2}\partial_i\partial_j g
$$

となる。添字 $j$ については和を取る。応答を検算するため、任意の一定ベクトル $\vb{a}$ に対する $\vb{G}\vb{a}$ に $\nabla\times\nabla\times-k^2$ を作用させる。勾配項に回転を二度掛けたものは零であり、残った $\nabla(\vb{a}\cdot\nabla g)$ は $-k^2$ を勾配項に掛けた寄与と相殺する。その結果は $-\vb{a}(\nabla^2+k^2)g=\vb{a}\delta^3$ である。

源と観測点が異なるときは、$x=kR$、$\vb{n}=(\vb r-\vb r')/R$ を使い、動径関数の Hessian を $ (g''-g'/R)n_i n_j+(g'/R)\delta_{ij}$ と書ける。$g'/g=ik-1/R$ を代入すると

$$
G_{ij}=g\left[
\left(1+\frac{i}{x}-\frac{1}{x^2}\right)\delta_{ij}
+\left(-1-\frac{3i}{x}+\frac{3}{x^2}\right)n_i n_j
\right]\qquad(R>0).
$$

遠方の $1/R$ 項では電流の横成分だけが残る。一方、連続した電流分布の内部にもこの式を積分するときは、球を除く主値積分に $-\delta_{ij}\delta^3/(3k^2)$ を加える。これは $g$ の特異な部分が $1/(4\pi R)$ であるためで、その二階微分のデルタ項は静電双極子の場合と同じである。従って、源の外側で成立する通常関数の式と、源をまたいで積分できる式には区別が必要である。

:::example[短い電流素片を軸上と赤道面から観測する]
$kR=1$ の二つの観測点を取り、一方では電流を $\vb{n}$ に平行、他方では垂直に見る。共通の $g$ と電流振幅を除くと、平行成分の係数は二つの括弧の和 $2-2i$、垂直成分は $i$ である。この距離では軸上の電場も大きい。これに対して $kR\gg1$ では平行成分の係数は $O(1/(kR))$ に落ち、垂直成分は $1$ に近づく。放射域で横波になることと、近くで縦の電場を測れることが一つの核から得られる。
:::

### 平面上の短い電流パルスと因果的な跳び

一様な面電流 $K(t)\hat{\vb x}$ を真空の $z=0$ に流す。面内の発散が零なので $\rho=0$、$\phi=0$ とできる。遅延解は左右へ対称に出ていき、$A_x$ は $t-|z|/c$ だけに依存する。波動方程式を面の直前直後で積分すると $[\partial_z A_x]=-\mu_0 K$ である。$A_x=F(t-|z|/c)$ なら左辺は $-2F'/c$ なので

$$
A_x=\frac{Z_0}{2}\int_{-\infty}^{t-|z|/c}K(s)\dd s,\qquad
E_x=-\frac{Z_0}{2}K(t-|z|/c),\quad B_y=\frac{\sgn z}{c}E_x.
$$

$Z_0=\mu_0 c$ である。パルスが去った後、電流の時間積分が非零なら $A_x$ は一定値にとどまるが、場は零に戻る。この一定ベクトルは一様な空間勾配として除ける。従って波形の到着をポテンシャルの絶対値だけで判断することはできない。電流を維持する外部源は、自己電場に逆らって仕事をし、単位面積あたり $-K E_x=Z_0 K^2/2$ を供給する。これは左右それぞれへ流れる $E_x^2/Z_0$ の二倍に等しい。

### Fourier 成分で比較する二つのゲージ

空間と時間に $e^{i(kz-\omega t)}$ を用い、$k>0$ とする。$x=\omega/(ck)$、$d=1-x^2$、$\rho=\rho_* r$、$\vb J=c\rho_* (j_t,0,xr)$ と置く。最後の成分は連続の式 $kJ_z=\omega\rho$ で決まる。$V_*=\rho_*/(\varepsilon_0 k^2)$ で $\phi$ を、$V_*/c$ で $\vb A$ を規格化する。$d\ne0$ の定常な特解を選ぶと Lorenz ゲージでは $(v,a_t,a_z)=(r,j_t,xr)/d$ である。Coulomb ゲージで同じ三成分は $(r,j_t/d,0)$ になる。ゲージ関数の振幅を $V_*/(ck)$ の単位で $h$ と置けば $v_C=v_L+ixh$、$a_{z,C}=a_{z,L}+ih$ であり、$h=ia_{z,L}$ を選べばよい。

観測する電場と磁場の規格化は $kV_*$、$kV_*/c$ である。$e_t=ix a_t$、$e_z=-iv+ix a_z$、$b_y=ia_t$ を計算し、特に縦電場に $1/d$ の見かけの極が残るかを調べる。これは外部から指定された一つの Fourier 源の特解であり、初期条件から定まる自由波を加える余地がある。

:::example[静止した電荷格子と横電流]
空間的に正弦波状の電荷と、それに直交する定常電流を重ねる。電極が作る電場と閉じた帰路を持つ電流の磁場を別々に静的に求めた結果に戻るかを調べる。具体的には $x=0$、$r=1$、$j_t=2$ とする。$k=10.0\,\mathrm{m^{-1}}$、$V_*=1.00\,\mathrm V$ とすれば $\rho_*=\varepsilon_0 k^2V_*=0.885\,\mathrm{nC/m^3}$ である。まず縦電流を保存条件から求める。次に $d$ で割って各ポテンシャルの振幅を得る。

$$
\begin{aligned}
j_z&=x r=0,\quad d=1-x^2=1,\\
(v,a_t,a_z)_L&=\left(\begin{matrix}1 & 2 & 0\end{matrix}\right).
\end{aligned}
$$

空間と時間に同じ位相規約を使ってゲージ関数を作用させると、二つの補正は

$$
\begin{aligned}
h&=ia_{z,L}=0,\quad (ixh,0,ih)=\left(\begin{matrix}0 & 0 & 0\end{matrix}\right),\\
(v,a_t,a_z)_C&=\left(\begin{matrix}1 & 2 & 0\end{matrix}\right).
\end{aligned}
$$

観測量のうち縦電場には電位の勾配とベクトルポテンシャルの時間微分がともに現れる。それぞれの値を並べると

$$
\begin{aligned}
e_{z,L}&=\left(- i\right)+\left(0\right)=- i,\\
e_{z,C}&=- i,\quad (e_t,b_y)=\left(0,2 i\right).
\end{aligned}
$$

実数の場として同じ時空点で比較する。位相 $\theta=kz-\omega t$ の各行は、振幅に $e^{i\theta}$ を掛けてから実部を取った値である。

| $\theta$ | $E_x$ $[\mathrm{V/m}]$ | $E_z$ $[\mathrm{V/m}]$ | $B_y$ $[\mathrm{nT}]$ |
|---:|---:|---:|---:|
| $0^\circ$ | $0$ | $0$ | $0$ |
| $45^\circ$ | $0$ | $7.071$ | $-47.17$ |
| $90^\circ$ | $0$ | $10$ | $-66.71$ |
| $180^\circ$ | $0$ | $0$ | $0$ |

独立な検算として Gauss の法則は $ie_z=r$、Faraday の法則は $e_t=xb_y$ を要求する。ここで得た振幅はどちらも満たす。時間変化が消えると両ゲージは同じ特解になる。電位と縦電場、ベクトルポテンシャルと横磁場がそれぞれ別の源で決まる。
:::

:::example[ゆっくり移動する縦の電荷模様]
電荷の山と谷が光より遅く移る。連続の式に必要な縦電流だけを流し、横電流は零とする。具体的には $x=\frac{1}{5}$、$r=1$、$j_t=0$ とする。$k=10.0\,\mathrm{m^{-1}}$、$V_*=1.00\,\mathrm V$ とすれば $\rho_*=\varepsilon_0 k^2V_*=0.885\,\mathrm{nC/m^3}$ である。まず縦電流を保存条件から求める。次に $d$ で割って各ポテンシャルの振幅を得る。

$$
\begin{aligned}
j_z&=x r=\frac{1}{5},\quad d=1-x^2=\frac{24}{25},\\
(v,a_t,a_z)_L&=\left(\begin{matrix}\frac{25}{24} & 0 & \frac{5}{24}\end{matrix}\right).
\end{aligned}
$$

空間と時間に同じ位相規約を使ってゲージ関数を作用させると、二つの補正は

$$
\begin{aligned}
h&=ia_{z,L}=\frac{5 i}{24},\quad (ixh,0,ih)=\left(\begin{matrix}- \frac{1}{24} & 0 & - \frac{5}{24}\end{matrix}\right),\\
(v,a_t,a_z)_C&=\left(\begin{matrix}1 & 0 & 0\end{matrix}\right).
\end{aligned}
$$

観測量のうち縦電場には電位の勾配とベクトルポテンシャルの時間微分がともに現れる。それぞれの値を並べると

$$
\begin{aligned}
e_{z,L}&=\left(- \frac{25 i}{24}\right)+\left(\frac{i}{24}\right)=- i,\\
e_{z,C}&=- i,\quad (e_t,b_y)=\left(0,0\right).
\end{aligned}
$$

位相 $\theta=kz-\omega t$ に対応する実場を示す。

| $\theta$ | $E_x$ $[\mathrm{V/m}]$ | $E_z$ $[\mathrm{V/m}]$ | $B_y$ $[\mathrm{nT}]$ |
|---:|---:|---:|---:|
| $0^\circ$ | $0$ | $0$ | $0$ |
| $45^\circ$ | $0$ | $7.071$ | $0$ |
| $90^\circ$ | $0$ | $10$ | $0$ |
| $180^\circ$ | $0$ | $0$ | $0$ |

この源では電流が非零でも磁場は零である。縦電流と縦変位電流が相殺するため、磁場の回転を必要としない。
:::

:::example[横電流だけで励起する波]
電荷の偏りを作らず横方向の電流を周期的に流す。面全体で電流を同位相に駆動すれば、空間的な位相差は軸方向だけにある。具体的には $x=\frac{3}{5}$、$r=0$、$j_t=1$ とする。$k=10.0\,\mathrm{m^{-1}}$、$V_*=1.00\,\mathrm V$ とすれば $\rho_*=\varepsilon_0 k^2V_*=0.885\,\mathrm{nC/m^3}$ である。まず縦電流を保存条件から求める。次に $d$ で割って各ポテンシャルの振幅を得る。

$$
\begin{aligned}
j_z&=x r=0,\quad d=1-x^2=\frac{16}{25},\\
(v,a_t,a_z)_L&=\left(\begin{matrix}0 & \frac{25}{16} & 0\end{matrix}\right).
\end{aligned}
$$

空間と時間に同じ位相規約を使ってゲージ関数を作用させると、二つの補正は

$$
\begin{aligned}
h&=ia_{z,L}=0,\quad (ixh,0,ih)=\left(\begin{matrix}0 & 0 & 0\end{matrix}\right),\\
(v,a_t,a_z)_C&=\left(\begin{matrix}0 & \frac{25}{16} & 0\end{matrix}\right).
\end{aligned}
$$

観測量のうち縦電場には電位の勾配とベクトルポテンシャルの時間微分がともに現れる。それぞれの値を並べると

$$
\begin{aligned}
e_{z,L}&=\left(0\right)+\left(0\right)=0,\\
e_{z,C}&=0,\quad (e_t,b_y)=\left(\frac{15 i}{16},\frac{25 i}{16}\right).
\end{aligned}
$$

位相 $\theta=kz-\omega t$ に対応する実場を示す。

| $\theta$ | $E_x$ $[\mathrm{V/m}]$ | $E_z$ $[\mathrm{V/m}]$ | $B_y$ $[\mathrm{nT}]$ |
|---:|---:|---:|---:|
| $0^\circ$ | $0$ | $0$ | $0$ |
| $45^\circ$ | $-6.629$ | $0$ | $-36.85$ |
| $90^\circ$ | $-9.375$ | $0$ | $-52.12$ |
| $180^\circ$ | $0$ | $0$ | $0$ |

縦のゲージ自由度で横ポテンシャルは消せない。電荷のない源領域でも電場はあり、ここでは時間変化する横電流から生まれる。
:::

:::example[縦源と横源を四分の一周期ずらす]
同じ波数において、横電流の駆動時刻を電荷密度の最大より四分の一周期ずらす。複素係数の虚部を保ったまま足して最後に実部を取る。具体的には $x=\frac{3}{4}$、$r=1$、$j_t=i$ とする。$k=10.0\,\mathrm{m^{-1}}$、$V_*=1.00\,\mathrm V$ とすれば $\rho_*=\varepsilon_0 k^2V_*=0.885\,\mathrm{nC/m^3}$ である。まず縦電流を保存条件から求める。次に $d$ で割って各ポテンシャルの振幅を得る。

$$
\begin{aligned}
j_z&=x r=\frac{3}{4},\quad d=1-x^2=\frac{7}{16},\\
(v,a_t,a_z)_L&=\left(\begin{matrix}\frac{16}{7} & \frac{16 i}{7} & \frac{12}{7}\end{matrix}\right).
\end{aligned}
$$

空間と時間に同じ位相規約を使ってゲージ関数を作用させると、二つの補正は

$$
\begin{aligned}
h&=ia_{z,L}=\frac{12 i}{7},\quad (ixh,0,ih)=\left(\begin{matrix}- \frac{9}{7} & 0 & - \frac{12}{7}\end{matrix}\right),\\
(v,a_t,a_z)_C&=\left(\begin{matrix}1 & \frac{16 i}{7} & 0\end{matrix}\right).
\end{aligned}
$$

観測量のうち縦電場には電位の勾配とベクトルポテンシャルの時間微分がともに現れる。それぞれの値を並べると

$$
\begin{aligned}
e_{z,L}&=\left(- \frac{16 i}{7}\right)+\left(\frac{9 i}{7}\right)=- i,\\
e_{z,C}&=- i,\quad (e_t,b_y)=\left(- \frac{12}{7},- \frac{16}{7}\right).
\end{aligned}
$$

位相 $\theta=kz-\omega t$ に対応する実場を示す。

| $\theta$ | $E_x$ $[\mathrm{V/m}]$ | $E_z$ $[\mathrm{V/m}]$ | $B_y$ $[\mathrm{nT}]$ |
|---:|---:|---:|---:|
| $0^\circ$ | $-17.14$ | $0$ | $-76.24$ |
| $45^\circ$ | $-12.12$ | $7.071$ | $-53.91$ |
| $90^\circ$ | $0$ | $10$ | $0$ |
| $180^\circ$ | $17.14$ | $0$ | $76.24$ |

相対位相は観測可能であり、一つのゲージで虚数に見えた項を一律に捨てる操作では保存できない。
:::

:::example[光円錐の近くで縦電場を取り出す]
光速に非常に近い位相速度を持つ強制源を指定する。弱い横駆動も分母の小ささによって増幅される一方、縦ポテンシャルの大きな二項は電場を作る際に相殺する。具体的には $x=\frac{99}{100}$、$r=1$、$j_t=\frac{1}{100}$ とする。$k=10.0\,\mathrm{m^{-1}}$、$V_*=1.00\,\mathrm V$ とすれば $\rho_*=\varepsilon_0 k^2V_*=0.885\,\mathrm{nC/m^3}$ である。まず縦電流を保存条件から求める。次に $d$ で割って各ポテンシャルの振幅を得る。

$$
\begin{aligned}
j_z&=x r=\frac{99}{100},\quad d=1-x^2=\frac{199}{10000},\\
(v,a_t,a_z)_L&=\left(\begin{matrix}\frac{10000}{199} & \frac{100}{199} & \frac{9900}{199}\end{matrix}\right).
\end{aligned}
$$

空間と時間に同じ位相規約を使ってゲージ関数を作用させると、二つの補正は

$$
\begin{aligned}
h&=ia_{z,L}=\frac{9900 i}{199},\quad (ixh,0,ih)=\left(\begin{matrix}- \frac{9801}{199} & 0 & - \frac{9900}{199}\end{matrix}\right),\\
(v,a_t,a_z)_C&=\left(\begin{matrix}1 & \frac{100}{199} & 0\end{matrix}\right).
\end{aligned}
$$

観測量のうち縦電場には電位の勾配とベクトルポテンシャルの時間微分がともに現れる。それぞれの値を並べると

$$
\begin{aligned}
e_{z,L}&=\left(- \frac{10000 i}{199}\right)+\left(\frac{9801 i}{199}\right)=- i,\\
e_{z,C}&=- i,\quad (e_t,b_y)=\left(\frac{99 i}{199},\frac{100 i}{199}\right).
\end{aligned}
$$

位相 $\theta=kz-\omega t$ に対応する実場を示す。

| $\theta$ | $E_x$ $[\mathrm{V/m}]$ | $E_z$ $[\mathrm{V/m}]$ | $B_y$ $[\mathrm{nT}]$ |
|---:|---:|---:|---:|
| $0^\circ$ | $0$ | $0$ | $0$ |
| $45^\circ$ | $-3.518$ | $7.071$ | $-11.85$ |
| $90^\circ$ | $-4.975$ | $10$ | $-16.76$ |
| $180^\circ$ | $0$ | $0$ | $0$ |

有限時間の駆動では極に近い系が定常に達するまで時間がかかる。厳密に光円錐上でこの有理式の分母を零と置くことはできない。
:::

:::example[位相速度が光速を超える駆動模様]
空間の各場所の独立した源に指定した位相差を与え、波数を小さくする。位相模様の速度は個々の担体の移動速度を意味しない。具体的には $x=\frac{3}{2}$、$r=1$、$j_t=1$ とする。$k=10.0\,\mathrm{m^{-1}}$、$V_*=1.00\,\mathrm V$ とすれば $\rho_*=\varepsilon_0 k^2V_*=0.885\,\mathrm{nC/m^3}$ である。まず縦電流を保存条件から求める。次に $d$ で割って各ポテンシャルの振幅を得る。

$$
\begin{aligned}
j_z&=x r=\frac{3}{2},\quad d=1-x^2=- \frac{5}{4},\\
(v,a_t,a_z)_L&=\left(\begin{matrix}- \frac{4}{5} & - \frac{4}{5} & - \frac{6}{5}\end{matrix}\right).
\end{aligned}
$$

空間と時間に同じ位相規約を使ってゲージ関数を作用させると、二つの補正は

$$
\begin{aligned}
h&=ia_{z,L}=- \frac{6 i}{5},\quad (ixh,0,ih)=\left(\begin{matrix}\frac{9}{5} & 0 & \frac{6}{5}\end{matrix}\right),\\
(v,a_t,a_z)_C&=\left(\begin{matrix}1 & - \frac{4}{5} & 0\end{matrix}\right).
\end{aligned}
$$

観測量のうち縦電場には電位の勾配とベクトルポテンシャルの時間微分がともに現れる。それぞれの値を並べると

$$
\begin{aligned}
e_{z,L}&=\left(\frac{4 i}{5}\right)+\left(- \frac{9 i}{5}\right)=- i,\\
e_{z,C}&=- i,\quad (e_t,b_y)=\left(- \frac{6 i}{5},- \frac{4 i}{5}\right).
\end{aligned}
$$

位相 $\theta=kz-\omega t$ に対応する実場を示す。

| $\theta$ | $E_x$ $[\mathrm{V/m}]$ | $E_z$ $[\mathrm{V/m}]$ | $B_y$ $[\mathrm{nT}]$ |
|---:|---:|---:|---:|
| $0^\circ$ | $0$ | $0$ | $0$ |
| $45^\circ$ | $8.485$ | $7.071$ | $18.87$ |
| $90^\circ$ | $12$ | $10$ | $26.69$ |
| $180^\circ$ | $0$ | $0$ | $0$ |

横応答の分母は符号を変えるが、縦の Gauss 拘束にはこの符号変化がない。情報の到達速度を調べるには単色特解でなく駆動開始時刻を持つ遅延解を見る。
:::


:::example[双極子の二点の遅延を同一時刻に丸める誤差]
遠方の観測者へ向かう線上に間隔 $d$ で二つの同位相の単色素片を置く。振幅の $1/R$ 因子は共通とするが遅延の差は保つ。それぞれの位相は中心の位相から $\pm kd/2$ ずれ、和は中心の振幅の二倍に $\cos(kd/2)$ を掛けたものとなる。ここでは両素片が同じ非零の横放射振幅を観測方向に作るとする。

$x=kd/2$、$Y$ は受信電場の符号付き同相振幅、同位置に置いたときの基準値を $Y_*=2.00\,\mathrm{mV/m}$ とする。

$$
f(x)=\cos{\left(x \right)},\qquad Y=Y_* f(x).
$$

$x=0$ のまわりで先頭の二項を残した近似は

$$
f_{\rm a}(x)=1 - \frac{x^{2}}{2}.
$$

次の表の差は $Y_*[f(x)-f_{\rm a}(x)]$ である。

| $x$ | $f(x)$ | $Y$ $[\mathrm{mV/m}]$ | 厳密値と近似値の差 $[\mathrm{mV/m}]$ |
|---:|---:|---:|---:|
| $\frac{1}{10}$ | $0.995$ | $1.99$ | $8.331\times10^{-6}$ |
| $\frac{1}{2}$ | $0.8776$ | $1.755$ | $0.005165$ |
| $1$ | $0.5403$ | $1.081$ | $0.0806$ |
| $\frac{3}{2}$ | $0.07074$ | $0.1415$ | $0.3915$ |

距離が電流源の大きさより十分大きくても、電流源の両端の位相差が無視できるとは限らない。遠方近似と長波長近似で必要とする長さの比較が異なることを、振幅の減少から直接判定できる。
:::

:::example[近接する鏡像からの反射パルスを分離する]
平面に Dirichlet 条件を課したスカラー波の点源と受信器を同じ高さ $h$ に置き、水平距離を $l$ とする。直接波の距離は $l$、鏡像からは $\sqrt{l^2+4h^2}$ なので到着差は両者の差を $c$ で割った時間である。

$x=h/l$、$l=3.00\,\mathrm m$、$Y$ は到着差、$Y_*=l/c=10.01\,\mathrm{ns}$ とする。鏡像の符号は負だが到着時刻を変えない。

$$
f(x)=\sqrt{4 x^{2} + 1} - 1,\qquad Y=Y_* f(x).
$$

変数を少し動かしたときの読みの変化は、接線の傾きから計算できる。基準値 $x_0$ の近くで

$$
\begin{aligned}
\frac{\Delta Y}{Y_*}= f\!\left(x_0+\Delta x\right)-f(x_0)= f'(x_0)\Delta x+O((\Delta x)^2),\\
f'(x)=\frac{4 x}{\sqrt{4 x^{2} + 1}}.
\end{aligned}
$$


| $x$ | $f(x)$ | $Y$ $[\mathrm{ns}]$ | $\dd Y/\dd x$ $[\mathrm{ns}]$ |
|---:|---:|---:|---:|
| $0$ | $0$ | $0$ | $0$ |
| $\frac{1}{10}$ | $0.0198$ | $0.1982$ | $3.925$ |
| $\frac{1}{2}$ | $0.4142$ | $4.145$ | $14.15$ |
| $1$ | $1.236$ | $12.37$ | $17.9$ |

パルス幅が到着差より大きければ二つは重なる。高さが低いほど二つの距離が近づくので、反射の振幅だけでなく時間分解能を考える必要がある。電磁場へ応用するときは各成分の接線・法線の境界条件を先に確かめる。
:::

### 円筒波で遅延核が波面の後にも残る理由

三次元では点源の瞬間的な変化は鋭い球面上に届く。一方、無限に長い線源を横から見ると、同時刻に変化した線の各点からの距離が違うため、最初の到着後にも信号が続く。三次元の核を軸方向に積分すれば、この違いを正確に求められる。軸からの距離を $s>0$、経過時間を $t$ とし、波動演算子を $\nabla^2-c^{-2}\partial_t^2$ とする。点源と同じ規格化で二次元の遅延核は

$$
G_2(s,t)=\int_{-\infty}^\infty
\frac{\delta(t-\sqrt{s^2+z^2}/c)}{4\pi\sqrt{s^2+z^2}}\dd z .
$$

$t<s/c$ には根がない。$t>s/c$ には $z=\pm\sqrt{c^2t^2-s^2}$ の二根があり、デルタの変数変換で各点に $c\sqrt{s^2+z^2}/|z|$ が掛かる。従って

$$
G_2(s,t)=\frac{c\,\Theta(t-s/c)}{2\pi\sqrt{c^2t^2-s^2}}.
$$

波頭では $(t-s/c)^{-1/2}$ に比例するが、有限時間幅で積分すれば有限である。時間が長いと $1/(2\pi t)$ に近づく。さらに一方向を積分し、残した直線座標を $x$ とすれば、$ct>|x|$ では $\int_{-b}^{b}\dd y/\sqrt{b^2-y^2}=\pi$、$b^2=c^2t^2-x^2$ より $G_1(x,t)=c\Theta(t-|x|/c)/2$ となる。したがって平面上の電流パルスに対するベクトルポテンシャルに、電流の時間積分が現れたことも説明できる。

この構成は無限に広い同時駆動源の理想化であって、情報が線上を瞬間に伝わったわけではない。有限長 $2b_0$ の一様線分を中央の垂線上で見ると、端からの最後の到着は $t=\sqrt{s^2+b_0^2}/c$ であり、デルタ状の源に対する尾はそこで切れる。

```tikz caption="同じ信号時刻を与える線源上の二点。最初の到着後は、より遠い線素が順に寄与する。"
\begin{tikzpicture}
\draw[thick] (-2.3,0)--(2.3,0) node[right] {$z$};
\fill (0,1.45) circle (2pt) node[above] {$P$};
\draw[dashed] (0,0)--(0,1.45) node[midway,right] {$s$};
\draw (-1.45,0)--(0,1.45)--(1.45,0) node[midway,above right] {$ct$};
\fill (-1.45,0) circle (2pt);\fill (1.45,0) circle (2pt);
\end{tikzpicture}
```

## 演習問題

::::exercise{level=basic}
$\chi=\alpha xt$ のゲージ変換で $\phi,\vb{A}$ の変化を求めよ。

:::solution
時間微分は $\alpha x$、勾配は $\alpha t\hat{\vb{x}}$ である。
従って $\Delta\phi=-\alpha x$、

$$
\Delta\vb{A}=\alpha t\hat{\vb{x}}
$$

となり、二項が電場で打ち消し合う。
:::
::::

::::exercise{level=basic}
$\vb{A}=bx\hat{\vb{x}}$ を Coulomb ゲージに移す時間独立な $\chi$ を一つ求めよ。

:::solution
発散は $b$ なので $\nabla^2\chi=-b$ が条件である。
$\chi=-bx^2/2$ なら $\vb{A}'=0$ で条件を満たす。
:::
::::

::::exercise{level=basic}
光が真空中の $1.00\,\mathrm{km}$ を進む時間を求めよ。

:::solution
遅延時刻との差は $R/c$ である。

$$
1000/(299792458)=3.34\times10^{-6}\,\mathrm s
$$

となる。
:::
::::

::::exercise{level=basic}
$\vb{A}=A_0\cos(kz-\omega t)\hat{\vb{x}}$、$\phi=0$ の場を求めよ。

:::solution
$\vb{E}=-\partial_t\vb{A}=-\omega A_0\sin(kz-\omega t)\hat{\vb{x}}$、

$$
\vb{B}=\partial_zA_x\hat{\vb{y}}=-kA_0\sin(kz-\omega t)\hat{\vb{y}}
$$

となる。

真空進行波なら $\omega=ck$ で振幅比は $E/B=c$ である。
:::
::::

::::exercise{level=advanced}
Lorenz ゲージを保ちつつポテンシャルを変える非定数の関数を構成せよ。

:::solution
$\chi=f(z-ct)$ を選ぶと $\partial_z^2\chi=f''$、$c^{-2}\partial_t^2\chi=f''$ である。
従って $\mathcal L\chi=0$ となり、勾配や時間微分が非零でも Lorenz 条件を保つ。
:::
::::

::::exercise{level=advanced}
空間的に局在する保存電流の遅延解が Lorenz 条件を満たすことを、畳み込みを用いて証明せよ。

:::solution
時空並進に依存する核では、観測点に関する微分は源座標の微分の負である。
部分積分で微分を源へ移し、空間無限遠と過去境界の項を零とすれば

$$
\nabla\cdot\vb{A}+c^{-2}\partial_t\phi=\mu_0G_{\rm ret}*(\nabla'\cdot\vb{J}+\partial_{t'}\rho)=0
$$

となる。
:::
::::

::::exercise{level=advanced}
$\rho=0$ の領域で、任意の真空波形 $f(z-ct)$ に対して $\vb{A}=f\hat{\vb{z}}$、$\phi=cf$ が作る場を求めよ。

:::solution
縦向きポテンシャルは $z,t$ だけに依存するので回転は零である。
電場の二項は $-cf'\hat{\vb{z}}$ と $cf'\hat{\vb{z}}$ で相殺する。
実際に $\chi=\int^{z-ct} f(s)\dd s$ を零ポテンシャルに施すとこの組が得られる。
:::
::::

::::exercise{level=advanced}
単色源 $\rho(\vb{r}',t)=\Re[\rho_\omega(\vb{r}')e^{-i\omega t}]$ の電位振幅の核を求め、遠方での境界条件を述べよ。

:::solution
時刻を $t-R/c$ にすれば位相因子に $e^{ikR}$ が掛かる。
従って核は $e^{ikR}/(4\pi R)$、$k=\omega/c$ である。
局在源に対して無限遠では $\phi_\omega\sim F(\Omega)e^{ikr}/r$ となり、$r(\partial_r-ik)\phi_\omega\to0$ を満たす。

これは外向き波の Sommerfeld 条件である。
:::
::::

::::exercise{level=advanced}
有限な周期箱で空間 Fourier 成分を用い、Coulomb ゲージのポテンシャルを、波数 $\vb{k}\ne0$、周波数 $\omega$ の電荷と電流の振幅から求めよ。
時間因子は $e^{-i\omega t}$ とし、まず $\omega^2\ne c^2 k^2$ の定常解を求めよ。さらに共鳴時に横電流が自由波に及ぼす効果を述べよ。

:::solution
$\nabla\cdot\vb{A}_C=0$ なので $\vb{k}\cdot\vb{A}_C=0$ である。
Poisson 方程式は

$$
-k^2\phi_C=-\rho/\varepsilon_0
$$

を与える。
一方、電流の横成分は

$$
\vb{J}_T=\vb{J}-\vb{k}(\vb{k}\cdot\vb{J})/k^2
$$

である。

従って

$$
\phi_C=\frac{\rho}{\varepsilon_0k^2},\qquad \vb{A}_C=\frac{\mu_0\vb{J}_T}{k^2-\omega^2/c^2}.
$$

$$
\vb{E}=-i\vb{k}\phi_C+i\omega\vb{A}_C
$$

より縦成分に対して $i\vb{k}\cdot\vb{E}=\rho/\varepsilon_0$ が再現される。
連続の式 $\vb{k}\cdot\vb{J}=\omega\rho$ によって縦電流と電荷は独立ではない。

共鳴時には時間領域の式 $\ddot{\vb A}_T+c^2 k^2\vb A_T=\mu_0 c^2\vb J_T e^{-i\omega t}$ に戻る。$\omega^2=c^2 k^2$ では $t e^{-i\omega t}$ に左辺の微分演算子を掛けると $-2i\omega e^{-i\omega t}$ となる。従って $\vb J_T\ne0$ なら特解は $-\mu_0 c^2\vb J_T t e^{-i\omega t}/(2i\omega)$ であり、有限振幅の単色定常解はない。$\vb J_T=0$ なら共鳴モードは自由波である。いずれの場合も周期箱での斉次解は初期条件で決まる。なお、除外した $\vb k=0$ に対しては、Poisson 方程式を箱全体で積分すると全電荷が零という整合条件が得られる。
:::
::::

::::exercise{level=advanced}
非斉次波動方程式 $\mathcal L f=-s$ について、零の初期値から始める場合と、同じ源に非零の初期値を加えた場合を比較せよ。
正則で減衰する場を考える。

:::solution
最初の時刻 $t_1$ より前で源を零に延長し、自由空間の $G_{\rm ret}$ と畳み込むと、因果性から $t<t_1$ では $f=0$ となる。
有限で滑らかな源の立ち上がりなら $f(t_1)=\dot f(t_1)=0$ を満たす。

別の解との差 $h$ は $\mathcal Lh=0$ に従う。
各空間 Fourier 成分では $\ddot h_{\vb{k}}+c^2k^2 h_{\vb{k}}=0$ なので、

$$
h_{\vb{k}}(t)=h_{\vb{k}}(t_1)\cos(ck\Delta t)+\dot h_{\vb{k}}(t_1)\frac{\sin(ck\Delta t)}{ck}.
$$

$\Delta t=t-t_1$ であり、$k=0$ では第二因子を極限の $\Delta t$ とする。
遅延積分にこの斉次解を足すと任意の初期値に合わせられる。
電磁場ではさらに Gauss 則とゲージ条件への適合を確かめる必要がある。
:::
::::

::::exercise{level=advanced}
真空中で電気双極子の点源を $\rho(\vb{r},t)=-\vb{p}(t)\cdot\nabla\delta^3(\vb{r})$、$\vb{J}(\vb{r},t)=\dot{\vb{p}}(t)\delta^3(\vb{r})$ と表す。

電荷保存を確かめ、Lorenz ゲージで $r>0$ の電位を求めよ。

:::solution
$\partial_t\rho=-\dot{\vb{p}}\cdot\nabla\delta^3$ と $\nabla\cdot\vb{J}=\dot{\vb{p}}\cdot\nabla\delta^3$ が消し合う。
空間・時間を独立な積分変数として、スカラーポテンシャルの空間デルタの微分を部分積分で核へ移す。

$\vb{n}=\vb{r}/r$、$t_0=t-r/c$ とおくと、源時刻積分後に

$$
\phi=\frac1{4\pi\varepsilon_0}\left[\frac{\vb{n}\cdot\vb{p}(t_0)}{r^2}+\frac{\vb{n}\cdot\dot{\vb{p}}(t_0)}{cr}\right],\qquad \vb{A}=\frac{\mu_0}{4\pi r}\dot{\vb{p}}(t_0)
$$

を得る。
第一項は $1/R$ の空間微分、第二項は遅延核の時間デルタを微分した寄与である。
静止した双極子なら第二項と $\vb{A}$ が消える。
:::
::::

::::exercise{level=advanced}
長さ $6.00\,\mathrm m$ の線分全体に瞬間的なスカラー源を置く。中央から垂直距離 $4.00\,\mathrm m$ の点で応答が続く時間を求めよ。

:::solution
最も短い経路は中央からの $4.00\,\mathrm m$、最長は端からの $\sqrt{4^2+3^2}=5.00\,\mathrm m$ である。従って時間幅は $(5.00-4.00)\,\mathrm m/c=3.34\,\mathrm{ns}$ である。この幅は源そのものの持続時間を零にしても残る。線分を一点に置き換える近似は、時刻をこの幅より粗くしか測れない場合に限って妥当である。
:::
::::

## 確認テスト

::::quiz{answer="3.34" tol="0.01"}
真空で $1.00\,\mathrm{km}$ 先に届く遅延を $\mathrm{\mu s}$ で答えよ。

:::explanation
距離を光速で割ると

$$
1000/c=3.3356\times10^{-6}\,\mathrm s
$$

であり、$3.34\,\mathrm{\mu s}$ となる。
:::
::::

::::quiz
Coulomb ゲージ条件を選べ。

- [ ] $\phi=0$
- [x] $\nabla\cdot\vb{A}=0$
- [ ] $\nabla\times\vb{A}=0$

:::explanation
Coulomb ゲージはベクトルポテンシャルの発散を零にする。
その結果スカラーポテンシャルが各時刻の Poisson 方程式を満たす。
:::
::::

::::quiz{answer="2" tol="1e-8"}
真空を伝わる非零周波数の電磁平面波の独立な偏光の数はいくつか。

:::explanation
残余ゲージで $\phi_0=0$ とすると $\vb{k}\cdot\vb{A}_0=0$ である。
三つの空間成分のうち波数に垂直な二成分が残る。
:::
::::

::::quiz
遅延核の時間引数を選べ。
$\tau=t-t'$ とする。

- [ ] $\delta(\tau+R/c)$
- [x] $\delta(\tau-R/c)$
- [ ] $\delta(\tau)$

:::explanation
過去の時刻 $t'=t-R/c$ にある源から影響が届くので $\tau=R/c$ を選ぶデルタ分布になる。
:::
::::

::::quiz
$e^{-i\omega t}$ の時間規約で、三次元自由空間の外向き球面波を表す動径因子を選べ。

- [ ] $e^{-ikr}/r$
- [ ] $e^{ikr}r$
- [x] $e^{ikr}/r$

:::explanation
位相 $kr-\omega t$ を一定にすると時間の増加につれて半径が増える。
振幅 $1/r$ は球面積の増加と合わせて全電力を一定にする。
:::
::::

::::quiz{answer="0.5" tol="1e-8"}
点状のスカラー源が半空間の接地面で反射する。
観測点まで直接は $2\,\mathrm m$、鏡映位置からは $4\,\mathrm m$ である。
反射パルスの振幅の絶対値を直接波で割れ。

:::explanation
自由空間遅延核の振幅は $1/(4\pi R)$ である。
像は符号が逆で絶対値は同じなので、振幅比は $(1/4)/(1/2)=1/2$ となる。
:::
::::

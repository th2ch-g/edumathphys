---
title: Maxwell 方程式
description: 変位電流を導入して場の方程式を閉じ、電荷・エネルギー・運動量の局所保存を具体的な収支から読み取る。
notation:
  - symbol: '\vb{J}_d'
    meaning: 変位電流密度
  - symbol: 'u'
    meaning: 場のエネルギー密度
  - symbol: '\vb{S}'
    meaning: エネルギー流密度
  - symbol: '\vb{g}'
    meaning: 運動量密度
  - symbol: 'T_{ij}'
    meaning: Maxwell 応力テンソル
---

## この章の目標

- 充電中の電極に張る二つの面で周回積分を一致させる。
- Maxwell 方程式の積分形と時間発展の整合性を示す。
- 場から電荷に渡る仕事と運動量を積分して求める。
- 応力テンソルで静電圧力と放射圧を計算する。

## 変位電流

### 充電電流が通る二つの面

コンデンサーのリード線を囲む閉路に面を張る。
一つの面はリード線を切り、もう一つは線を避けて極板間を通る。
前者には電流 $I$ が流れ、後者には真空なら伝導電流がない。
同じ周回積分に二つの答えが出るため、定常場の周回積分則 :eqref[eq-em-ampere] の適用条件を調べ直す必要がある。

電荷保存と $\nabla\cdot\vb{D}=\rho_f$ より

$$
\nabla\cdot(\vb{J}_f+\partial_t\vb{D})=0
$$

となる。
従って補正項 :term[変位電流密度]{en="displacement current density" yomi="へんいでんりゅうみつど"} $\vb{J}_d=\partial_t\vb{D}$ を加えれば、どちらの面でも同じ総電流を与える。

実験に合う場の法則は

$$
\nabla\times\vb{H}=\vb{J}_f+\partial_t\vb{D}.
\label{eq-em-ampere-maxwell}
$$

電荷保存が決めるのは補正項の発散である。
この形そのものには電磁場についての経験的内容がある。
真空中では $\vb{D}=\varepsilon_0\vb{E}$、$\vb{H}=\vb{B}/\mu_0$ である。
物質中では分極電流 $\partial_t\vb{P}$ も $\partial_t\vb{D}$ に含まれる。

:::example[円形極板を流れる変位電流]
半径 $a$ の平行平板で端を無視し、$I(t)=\dot Q$ とする。
間隙の $D_z=Q/(\pi a^2)$ は一様である。
半径 $s<a$ の円なら $I_d=I s^2/a^2$、従って $B_\varphi=\mu_0Is/(2\pi a^2)$ となる。

$s=a=5.00\,\mathrm{cm}$、$I=0.500\,\mathrm A$ なら

$$
B_\varphi=2.00\,\mathrm{\mu T}
$$

である。
変化が速く板間の場を一様とみなせない場合には、この近似解は使えない。
:::

:::example[伝導と変位の大きさ]
線形で Ohm 則に従う媒質に角周波数 $\omega$ の場を加えると、電流振幅の比は

$$
J_d/J_c=\omega\varepsilon/\sigma_c
$$

である。
$\varepsilon_r=4$、$\sigma_c=10^{-6}\,\mathrm{S/m}$ と仮定すれば両者が等しい周波数は

$$
\sigma_c/(2\pi\varepsilon)=4.49\,\mathrm{kHz}
$$

である。

二つの電流には位相差もある。
:::

### 電極の表面で電流がどうつながるか

固定界面に自由面電荷 $\sigma_f$ と面電流 $\vb{K}_f$ を許す。
界面上の小片を厚さの消える柱で囲んで電荷保存を適用すると

$$
\partial_t\sigma_f+\nabla_s\cdot\vb{K}_f+\vb{n}\cdot(\vb{J}_{f2}-\vb{J}_{f1})=0.
$$

$\nabla_s\cdot$ は面内での発散である。
$\vb{n}$ を $1$ から $2$ へ向けたとき、面電荷は下側から流入した量で増え、上側への流出と小片の縁を横切る面電流で減る。
Gauss 則の境界条件を時間微分すれば

$$
\partial_t\sigma_f=\vb{n}\cdot\partial_t(\vb{D}_2-\vb{D}_1)
$$

である。

従って横へ流れる電流がなければ、伝導電流が電極で途切れても、変位電流を加えた法線電流密度は連続する。

:::example[漏れのある平板コンデンサー]
一定の $\varepsilon,\sigma_c$ をもつ媒質が極板間を満たすと、周波数 $\omega$ の電流は $e^{i\omega t}$ の規約で $I=(G+i\omega C)V$、

$$
G=\sigma_c A/d,\ C=\varepsilon A/d
$$

である。

実効電圧 $V_{\rm rms}$ で発熱する平均電力は $G V_{\rm rms}^2$ である。
理想容量に入る電流は位相が四分の一周期進むため、一周期の正味の仕事は零になる。
:::

### 一本の導線を流れる電荷とギャップの磁束を同じ時間軸で追う

抵抗 $R$ を介し、初め無電荷の理想容量 $C$ を直流電源 $V_0$ につなぐ。
回路が小さくインダクタンスが無視できれば $R\dot Q+Q/C=V_0$ である。
従って

$$
Q(t)=CV_0(1-e^{-t/(RC)}),\qquad I_c(t)=(V_0/R)e^{-t/(RC)}.
$$

真空のギャップで面積分した電気変位は $Q(t)$ なので $I_d(t)=\dot Q(t)=I_c(t)$ である。
電池から出たすべての電子が真空を横切る必要はない。
電荷の移動は極板で止まり、板間の電場が変わる。

例えば $R=1.00\,\mathrm{k\Omega},\ C=1.00\,\mathrm{nF},\ V_0=5.00\,\mathrm V$ では、$t=1.00\,\mathrm{\mu s}$ に伝導電流と変位電流はともに $1.84\,\mathrm{mA}$、電荷は $3.16\,\mathrm{nC}$ である。

電源仕事は最終的に $CV_0^2$ となり、$CV_0^2/2$ がコンデンサーに蓄えられる。
残りは

$$
\int_0^\infty R I_c^2\dd t=CV_0^2/2
$$

として抵抗内で熱になる。

### 異なる緩和時間を持つ二つの層を充電する

平行板の間を、導電率 $\sigma_i$、誘電率 $\epsilon_i$、厚さ $d_i$ の均一な二層で満たす。界面の法線を一層目から二層目へ向け、電流と電場を同じ法線に沿う符号付き量で表す。側面から漏れがなければ、各層で共通なのは伝導電流でなく $j=\sigma_i E_i+\epsilon_i\dot E_i$ である。一方、電圧 $V$ は $d_1 E_1+d_2 E_2$ である。最初は場が零で、時刻零に電圧をステップ状に $V_0$ に変える。実際には試料を横切る伝播時間より遅く、抵抗緩和より速い立上りを理想化している。

直後には伝導電流の有限値を微小時間だけ積分しても有限な界面電荷を運べず、$\epsilon_2 E_2(0^+)=\epsilon_1 E_1(0^+)$ である。十分時間が経つと変位電流は消え、$\sigma_2 E_2(\infty)=\sigma_1 E_1(\infty)$ になる。電圧拘束から $E_2$ を消すと、

$$
(d_2\epsilon_1+d_1\epsilon_2)\dot E_1+
(d_2\sigma_1+d_1\sigma_2)E_1=\sigma_2 V_0,
$$

従って $E_1=E_{1,\infty}+(E_{1,0}-E_{1,\infty})e^{-t/\tau}$、$\tau=(d_2\epsilon_1+d_1\epsilon_2)/(d_2\sigma_1+d_1\sigma_2)$ である。界面自由電荷 $q_s=\epsilon_2E_2-\epsilon_1E_1$ の微分は $\sigma_1E_1-\sigma_2E_2$ に等しい。各層の $\epsilon/\sigma$ が同じなら直後の分配が定常分配にもなり、界面に自由電荷は蓄積しない。

:::example[同じ誘電率でも界面が帯電する]
$d_1=d_2=1.00\,\mathrm{mm}$、$\epsilon_1=\epsilon_2=2\varepsilon_0$、$\sigma_1=1.00\,\mathrm{nS/m}$、$\sigma_2=3.00\,\mathrm{nS/m}$、$V_0=4.00\,\mathrm V$ とする。直後の両層の場は $2.00\,\mathrm{kV/m}$、定常値は順に $3.00,1.00\,\mathrm{kV/m}$ である。従って $q_s$ は零から $-35.4\,\mathrm{nC/m^2}$ に時定数 $8.85\,\mathrm{ms}$ で近づく。負号は、初期に二層目へ出ていく伝導電流が一層目からやって来る電流より大きいことと整合する。
:::


:::example[漏れ電流と変位電流を電流計で区別する]
平板ギャップに導電率 $\sigma$、誘電率 $\epsilon$ の一様な物質を入れる。交流電場の実効値を $E$ とすれば、伝導成分 $\sigma E$ と変位成分 $\omega\epsilon E$ は直交位相である。従って実効電流密度は二成分の二乗和の平方根となる。

$x=\omega\epsilon/\sigma$、$\sigma=1.00\,\mathrm{\mu S/m}$、$E=1.00\,\mathrm{kV/m}$ とし、$Y_*=\sigma E=1.00\,\mathrm{mA/m^2}$、$Y$ を測定する実効値とする。

$$
f(x)=\sqrt{x^{2} + 1},\qquad Y=Y_* f(x).
$$

二点の読みの差を求めるには、基準を同じに保ち、振幅の換算を最後に一度だけ行う。最初と最後の設定で

$$
\Delta Y/Y_* = f\!\left(2\right)-f\!\left(0\right)=1.236.
$$

第四列は初期設定からの変化である。

| $x$ | $f(x)$ | $Y$ $[\mathrm{mA/m^2}]$ | $Y-Y(x_0)$ $[\mathrm{mA/m^2}]$ |
|---:|---:|---:|---:|
| $0$ | $1$ | $1$ | $0$ |
| $\frac{1}{2}$ | $1.118$ | $1.118$ | $0.118$ |
| $1$ | $1.414$ | $1.414$ | $0.4142$ |
| $2$ | $2.236$ | $2.236$ | $1.236$ |

$x=1$ で両成分の実効値が等しくなるが、電流計は二倍でなく $\sqrt2$ 倍を指す。熱として消費される平均電力は $\sigma E^2$ であり、実効電流の全量に $E$ を掛けたものとは異なる。
:::

## Maxwell 方程式

### 四本の式と力の法則

:::law[Maxwell 方程式]{#law-em-maxwell}
真空を基準に全電荷・全電流を $\rho,\vb{J}$ で表すと、局所方程式は

$$
\begin{aligned}
\nabla\cdot\vb{E}&=\rho/\varepsilon_0,&\nabla\cdot\vb{B}&=0,\\
\nabla\times\vb{E}&=-\partial_t\vb{B},&\nabla\times\vb{B}&=\mu_0\vb{J}+c^{-2}\partial_t\vb{E},\qquad c^{-2}=\varepsilon_0\mu_0.
\end{aligned}
\label{eq-em-maxwell}
$$

電荷の運動は Lorentz 力 $q(\vb{E}+\vb{v}\times\vb{B})$ に従う。
:::

巨視的な方程式では二つの電気的源を自由電荷に限って、$\nabla\cdot\vb{D}=\rho_f$、:eqref[eq-em-ampere-maxwell] とする。
他の二式は同じである。
実際に $\rho=\rho_f-\nabla\cdot\vb{P}$、$\vb{J}=\vb{J}_f+\partial_t\vb{P}+\nabla\times\vb{M}$ を代入すればよい。

電荷・電流を自由成分と全成分のどちらで数えているかで、使う物質量が決まる。

固定体積 $V$ では発散定理により $\oint_{\partial V}\vb{D}\cdot\dd\vb{S}=Q_f$、

$$
\oint_{\partial V}\vb{B}\cdot\dd\vb{S}=0
$$

である。
固定した向き付き面 $S$ では Stokes の定理により

$$
\oint_{\partial S}\vb{E}\cdot\dd\vb{l}=-\frac{\dd}{\dd t}\int_S\vb{B}\cdot\dd\vb{S},\qquad \oint_{\partial S}\vb{H}\cdot\dd\vb{l}=I_f(S)+\frac{\dd}{\dd t}\int_S\vb{D}\cdot\dd\vb{S}.
$$

この二つの積分定理の仮定と向き付けは [解析学: 線積分・面積分とベクトル解析](/math/analysis/vector-calculus/) に従う。
特異な面電荷・面電流がある場合は、場が滑らかな領域ごとに積分して面の寄与を足す。

### 初期値の制約と波

回転の発散は恒等的に零であるから、Faraday 則の発散をとると $\partial_t\nabla\cdot\vb{B}=0$ である。
同様に電荷保存を使えば

$$
\partial_t(\nabla\cdot\vb{E}-\rho/\varepsilon_0)=0
$$

となる。

従って二つの発散条件は初期に満たせば時間発展で保たれる。
逆に与えた源が連続の式を満たさなければ解はない。

真空の無源領域で Faraday 則の回転をとると、左辺は $\nabla(\nabla\cdot\vb{E})-\nabla^2\vb{E}$、右辺は $-c^{-2}\partial_t^2\vb{E}$ である。
従って

$$
\nabla^2\vb{E}-c^{-2}\partial_t^2\vb{E}=0
$$

となる。

磁束密度も同じ方程式に従う。
各成分を独立な波として選ぶ自由はなく、発散条件と二つの回転の式を同時に満たす必要がある。

:::example[波形を直接代入する]
$\vb{E}=f(z-ct)\hat{\vb{x}}$ とし、$\vb{B}=f(z-ct)\hat{\vb{y}}/c$ を組にする。
発散はともに零、

$$
\nabla\times\vb{E}=f'\hat{\vb{y}}=-\partial_t\vb{B}
$$

である。

残る式も

$$
\nabla\times\vb{B}=-f'\hat{\vb{x}}/c=c^{-2}\partial_t\vb{E}
$$

となる。
正弦波だけでなく任意の滑らかな波形が伝播できる。
:::

:::example[初期場が両方の発散条件を満たしていても静止するとは限らない]
真空で $t=0$ に $\vb{E}=E_0\cos(kz)\hat{\vb{x}}$、$\vb{B}=0$ を与える。
電荷も電流も零なので発散条件は満たす。
Faraday 則は $\partial_t B_y(0)=kE_0\sin(kz)$、Ampère–Maxwell 則は $\partial_t E_x(0)=0$ を指定する。

これらを満たす解は

$$
E_x=E_0\cos(kz)\cos(ckt),\qquad B_y=(E_0/c)\sin(kz)\sin(ckt).
$$

磁場が一度も変わらないと仮定すれば Faraday 則を破る。
初期値だけから最初の時間微分まで与えられることが、二階の波動方程式だけを解く場合との違いである。
:::

### 準静的近似の二つの極限

典型的な長さ $a$ と時間 $\tau$ をとる。
Faraday 則で生じる回転電場の尺度は $aB/\tau$ なので、これが主要な電場より小さければ、電場を各瞬間のスカラーポテンシャルから求められる。
一方、導線で変位電流の尺度 $\varepsilon E/\tau$ が伝導電流に比べて小さければ、磁場はその時の電流から求められる。

前者は小さい容量素子、後者は良導体の渦電流の計算などに対応する。
どちらも $a/(c\tau)$ が小さい条件で使われるが、コンデンサーの間隙には伝導電流がないため、回路全体から変位電流を一律に落とすことはできない。

### 動く不連続面での Maxwell 方程式

面が静止すると仮定して得た境界条件を、速度 $v_n$ で法線方向へ動く面にそのまま用いることはできない。
表面を $s(\vb{r},t)=0$ と表し、面で

$$
\nabla s=\vb{n},\ \partial_t s=-v_n
$$

とする。

$\vb{n}$ は媒質 $1$ から $2$ へ向く。
跳びを $[F]=F_2-F_1$ と書き、$F=F_1+[F]\Theta(s)$ とおくと、特異部分は $\nabla F\supset[F]\vb{n}\delta(s)$、$\partial_t F\supset-v_n[F]\delta(s)$ である。

正則な場に加え、自由電荷と電流に面分布 $\sigma_f\delta(s),\ \vb{K}_{\rm tot}\delta(s)$ があるとする。
各式のデルタ係数を比較すれば

$$
\begin{aligned}
\vb{n}\cdot[\vb{B}]&=0,&\vb{n}\cdot[\vb{D}]&=\sigma_f,\\
\vb{n}\times[\vb{E}]&=v_n[\vb{B}],&\vb{n}\times[\vb{H}]&=\vb{K}_{\rm tot}-v_n[\vb{D}].
\end{aligned}
$$

ここで $\vb{K}_{\rm tot}$ は移動する面電荷が運ぶ法線電流 $\sigma_f v_n\vb{n}$ も含む。
最後の式の右辺の法線成分は確かに相殺する。
$v_n=0$ とすると静止界面の四条件に戻る。
運動する媒質の構成式は別途必要である。

:::example[波面そのものを不連続面とみなす]
真空で電荷・電流がなく、静止した背景に微小な階段波面が進む場合を考える。
法線条件から二つの場の跳びは横向きである。
回転の二条件は $\vb{n}\times[\vb{E}]=v_n[\vb{B}]$、

$$
\vb{n}\times[\vb{B}]=-\mu_0\varepsilon_0v_n[\vb{E}]
$$

である。

第一式にもう一度 $\vb{n}\times$ を施すと

$$
-[\vb{E}]=-\mu_0\varepsilon_0 v_n^2[\vb{E}]
$$

である。
非零の跳びを運べる速さは $|v_n|=c$ に限られる。
:::

### 初期値を二つの伝播方向に分ける

真空で、すべての量が $z,t$ だけの関数で $E_x,B_y$ だけが非零とする。回転方程式は $\partial_t E_x=-c^2\partial_z B_y$、$\partial_t B_y=-\partial_z E_x$ となる。$f=(E_x+cB_y)/2$、$h=(E_x-cB_y)/2$ を作ると、$(\partial_t+c\partial_z)f=0$、$(\partial_t-c\partial_z)h=0$ が別々に成り立つ。従って初期値全体を $E_x(z,0)=e(z)$、$B_y(z,0)=b(z)$ とすれば

$$
E_x(z,t)=\tfrac12[e(z-ct)+cb(z-ct)+e(z+ct)-cb(z+ct)],
$$

磁場は同じ二つの波形の差を $c$ で割ったものとなる。たとえ初期の電場の形がわかっていても、磁場を指定しなければどちらへどれだけのエネルギーが流れるかは決まらない。これは波動方程式だけで時間二階の初期値を指定するとき、電場の初期時間微分に磁場の回転という制約があることの別表現である。

$p=f(z-ct),q=h(z+ct)$ と略記すれば、$u=\varepsilon_0(p^2+q^2)$、$S_z=\varepsilon_0 c(p^2-q^2)$ となる。電場同士には $2pq$ の干渉項があるが、磁場同士の交差項は $-2pq$ であり、全エネルギーで相殺する。電場だけを測る器具は定在波の腹と節を区別する一方、両方向のエネルギーを足した量には同じ節が現れるとは限らない。


:::example[初期磁場を変えて波束の進行方向を選ぶ]
初期に $cB_y=xe(z)$、$E_x=e(z)$ と準備する。右向きと左向きの電場は $(1+x)e/2$、$(1-x)e/2$ である。各波のエネルギーを積分して右向きの割合をとると波形の積分が相殺される。

$Y$ はパーセント表示した右向きエネルギー、$Y_*=100$ とする。例えば初期の総エネルギーを各設定で $1.00\,\mathrm{mJ}$ にそろえれば百分率の読みから右向きのエネルギーを求められる。

$$
f(x)=\frac{\left(x + 1\right)^{2}}{2 \left(x^{2} + 1\right)},\qquad Y=Y_* f(x).
$$

変数を少し動かしたときの読みの変化は、接線の傾きから計算できる。基準値 $x_0$ の近くで

$$
\begin{aligned}
\frac{\Delta Y}{Y_*}= f\!\left(x_0+\Delta x\right)-f(x_0)= f'(x_0)\Delta x+O((\Delta x)^2),\\
f'(x)=- \frac{\left(x - 1\right) \left(x + 1\right)}{\left(x^{2} + 1\right)^{2}}.
\end{aligned}
$$

同じ尺度で値と傾きを並べる。

| $x$ | $f(x)$ | $Y$ $[\%]$ | $\dd Y/\dd x$ $[\%]$ |
|---:|---:|---:|---:|
| $-1$ | $0$ | $0$ | $0$ |
| $- \frac{1}{2}$ | $0.1$ | $10$ | $48$ |
| $0$ | $0.5$ | $50$ | $100$ |
| $\frac{1}{2}$ | $0.9$ | $90$ | $48$ |
| $1$ | $1$ | $100$ | $0$ |

初期電場だけを指定し磁場を零にすればエネルギーは等分される。$x$ の符号を反転すると二つの進行方向が交換される。純進行波に近い場合、逆方向の小さな振幅はエネルギーには二次でしか現れない。
:::

## 電磁場のエネルギーと Poynting ベクトル

### 局所的な仕事をたどる

Lorentz 力の仕事率は一粒子で $q\vb{v}\cdot\vb{E}$、全粒子について体積あたり $\vb{J}\cdot\vb{E}$ である。
この供給源を Maxwell 方程式から取り出す。
Ampère–Maxwell 則と $\vb{E}$ の内積から、Faraday 則と $\vb{B}/\mu_0$ の内積を組み合わせる。

発散公式を使えば

$$
\partial_tu+\nabla\cdot\vb{S}=-\vb{J}\cdot\vb{E},\qquad u=\tfrac12(\varepsilon_0 E^2+B^2/\mu_0),\qquad \vb{S}=\vb{E}\times\vb{B}/\mu_0.
\label{eq-em-poynting}
$$

実際に $\nabla\cdot(\vb{E}\times\vb{B})=\vb{B}\cdot(\nabla\times\vb{E})-\vb{E}\cdot(\nabla\times\vb{B})$ の右辺を代入すると、時間微分二項と電流の仕事が現れる。

$\vb{S}$ は :term[Poynting ベクトル]{en="Poynting vector" yomi="ぽいんてぃんぐべくとる"} で単位は $\mathrm{W/m^2}$ である。
固定体積で積分したとき、場のエネルギーの増加と外へ出たエネルギーの和が、場が電荷にした仕事の負に等しい。

物質中では同じ計算を $\vb{H},\vb{D}$ で行うと

$$
\nabla\cdot(\vb{E}\times\vb{H})=-\vb{E}\cdot\vb{J}_f-\vb{E}\cdot\partial_t\vb{D}-\vb{H}\cdot\partial_t\vb{B}
$$

となる。

固定・非分散・無損失で対称な線形媒質なら後二項は $-\partial_t(\vb{E}\cdot\vb{D}+\vb{B}\cdot\vb{H})/2$ である。
分散があると物質の振動エネルギーも別に追跡する必要がある。

:::example[抵抗に横から入る電力]
円柱抵抗に $+z$ 向きの電流を流す。
側面では $E_z=V/\ell$、$H_\varphi=I/(2\pi a)$ なので $\vb{E}\times\vb{H}$ は半径の内向きである。
側面積を掛けた流入は

$$
(V/\ell)[I/(2\pi a)]2\pi a\ell=VI
$$

であり、内部の Joule 発熱に一致する。
:::

:::example[充電中の真空コンデンサー]
半径 $a$、間隔 $d\ll a$ の平行な円形板を、時間尺度 $\tau\gg a/c$ で充電する。縁の場の乱れを無視すると、準静的近似の最低次で板の縁において

$$
S_r=-EB_\varphi/\mu_0=-EI/(2\pi a)
$$

である。
間隔 $d$ とすれば縁からの流入は $2\pi ad|S_r|=EdI=VI$ となる。
一定容量では $I=C\dot V$ なので、この次数で流入は $\dd(CV^2/2)/\dd t$ に等しい。近似の精度を調べるには、極板間の $B_\varphi(s)=\mu_0 I s/(2\pi a^2)$ を使えばよい。$B_\varphi^2/(2\mu_0)$ を体積 $2\pi s d\,\dd s$ で積分すると $U_B=\mu_0 d I^2/(16\pi)$ である。$Q\ne0$ では $U_B/U_E=(aI/cQ)^2/8$ となる。電流が時間変化すると、この磁場も変化する。従って同じ次数の補正まで収支を取るには、Faraday の法則が要求する電場の空間変化も含める必要がある。
:::

### 正弦定常場で流れ込む複素電力

時間因子を $e^{-i\omega t}$ として振幅に添字 $0$ を付ける。
真空では $\nabla\times\vb{E}_0=i\omega\mu_0\vb{H}_0$、

$$
\nabla\times\vb{H}_0=\vb{J}_0-i\omega\varepsilon_0\vb{E}_0
$$

である。

複素共役をとる箇所を保ったまま発散を計算すると

$$
\nabla\cdot\left(\tfrac12\vb{E}_0\times\vb{H}_0^*\right)=-\tfrac12\vb{E}_0\cdot\vb{J}_0^*+2i\omega(\overline w_m-\overline w_e).
$$

ここで $\overline w_m=\mu_0|\vb{H}_0|^2/4$、

$$
\overline w_e=\varepsilon_0|\vb{E}_0|^2/4
$$

である。
実部は時間平均した仕事と流出の収支である。
虚部は磁気・電気の蓄積量の差を反映し、熱や遠方への正味の輸送を表す量ではない。

電場が大きい場所でも磁場と位相が四分の一周期ずれていると、平均電力は流れない。

:::example[空洞での電気・磁気エネルギー]
無源で損失のない完全導体空洞の固有振動を考える。
壁に接する電場が零なので、複素電力の面積分は虚部も含めて零である。
従って正の固有周波数では

$$
\int\overline w_e\,\dd^3 r=\int\overline w_m\,\dd^3r
$$

である。

これは体積全体を平均した結果であり、定在波の電場の節では局所的な二つの値は一般に異なる。
:::

### 回路を電力が伝わる有限の速さ

理想同軸線の短い区間 $\dd z$ に分布容量 $C'\dd z$ と分布インダクタンス $L'\dd z$ を割り当てる。
漏れと抵抗を無視する。
隣り合う断面での電位差の減少は直列インダクタの電圧降下、電流の減少は区間の充電に等しい。
従って

$$
\partial_z V=-L'\partial_t I,\qquad \partial_z I=-C'\partial_t V.
$$

第一式に $I$、第二式に $V$ を掛けて足すと

$$
\partial_z(VI)+\partial_t\left(\tfrac12L'I^2+\tfrac12C'V^2\right)=0.
$$

この一次元の保存則は、真空の Poynting の定理を線路の断面で積分した結果と同じである。
二式を組み合わせると $V,I$ のいずれも速さ $1/\sqrt{L'C'}$ の波に従う。
伝播方向が $+z$ の単独の波では $V/I=\sqrt{L'/C'}$ となる。

真空同軸線で電力が $c$ で伝わることは、回路を流れる電子の漂移速度が小さいことと両立する。

### 局在した自由電磁場のエネルギー中心

電磁場が速く減衰し、無限遠でエネルギー流に位置を掛けた境界項も消える場合を考える。
真空中で全エネルギー $U=\int u\,\dd^3r>0$、エネルギー中心 $\vb{R}_U=U^{-1}\int\vb{r} u\,\dd^3r$ を定める。

物質との交換がないなら $U$ は一定である。
局所保存則に座標 $x_i$ を掛けると

$$
U\dot R_{U,i}=-\int x_i\partial_j S_j\dd^3r=\int S_i\dd^3r.
$$

第二の等号は部分積分による。
運動量の節で導く真空の運動量密度 $\vb{g}=\vb{S}/c^2$ を用いれば

$$
U\dot{\vb{R}}_U=c^2\vb{P}_{\rm field}
$$

となる。
さらに電磁力を受ける物質がなく、外部への運動量流出もなければ $\vb{P}_{\rm field}$ は一定であり、エネルギー中心は等速で進む。

:::example[二方向へ分かれる光パルス]
互いに重ならず反対方向へ進む二つのほぼ平面状のパルスに $U_+,U_-$ のエネルギーがあるとする。
縦運動量は $(U_+-U_-)/c$ なので、エネルギー中心の速度は $c(U_+-U_-)/(U_++U_-)$ である。
$U_+=3U_-$ なら $c/2$ となる。

二つの波の各先端は真空で $c$ で進むが、エネルギーを重みとして平均した位置は異なる速さで移動する。
:::

### 導電シートが吸収できる電力

真空の同一の二半空間の間に、厚さを無視できる薄膜を置く。各周波数で接線電流を $K=G_s E_t$ と書けるとし、単位は $G_s$ が $\mathrm S$、$K$ が $\mathrm{A/m}$ である。左から垂直入射する平面波の電場振幅を一とし、反射・透過振幅を $r,t$ とする。右進行波では $B_y=E_x/c$ なので、$E_x/H_y=\mu_0c=Z_0$ が真空の波動インピーダンスである。左進行波ではこの比の符号が反転する。従って接線電場の連続から $1+r=t$、磁場の跳びから $(1-r-t)/Z_0=G_s t$ を得る。$g=Z_0G_s$ と置けば

$$
r=-\frac{g}{2+g},\qquad t=\frac{2}{2+g},\qquad
A=1-|r|^2-|t|^2=\frac{4\Re g}{|2+g|^2}.
\label{eq-em-sheet-absorption}
$$

これは体積の Poynting の定理をシートの厚さまで積分した結果とも比較できる。膜に直接渡る平均電力密度は $\Re G_s |tE_i|^2/2$ であり、入射強度 $|E_i|^2/(2Z_0)$ で割れば上式に一致する。ここで振幅は実効値でなくピーク値である。実数 $g>0$ で吸収が最大になるのは微分から $g=2$ のときで、その最大値は $1/2$ である。高い導電率で電場がほぼ零になると、電流が流れやすいのに吸収は減り、反射へエネルギーが移る。両側から位相をそろえて照射する場合には、一方向だけから得たこの上限を適用できない。


:::example[導電膜を通る電力と発熱の整合]
単一の入射波に対する実数のシートコンダクタンスを $G_s$ とし、$x=Z_0G_s/2$ と置く。接線電場と面電流を整合すると一方向からの吸収率は $2x/(1+x)^2$ である。

入射強度を $Y_*=100\,\mathrm{W/m^2}$ とする。$Y$ は単位面積あたりに発生する熱であり、$x=1$ に対応するシート抵抗は $Z_0/2\simeq188.4\,\Omega$ である。

$$
f(x)=\frac{2 x}{\left(x + 1\right)^{2}},\qquad Y=Y_* f(x).
$$

変数を少し動かしたときの読みの変化は、接線の傾きから計算できる。基準値 $x_0$ の近くで

$$
\begin{aligned}
\frac{\Delta Y}{Y_*}= f\!\left(x_0+\Delta x\right)-f(x_0)= f'(x_0)\Delta x+O((\Delta x)^2),\\
f'(x)=- \frac{2 \left(x - 1\right)}{\left(x + 1\right)^{3}}.
\end{aligned}
$$


| $x$ | $f(x)$ | $Y$ $[\mathrm{W/m^2}]$ | $\dd Y/\dd x$ $[\mathrm{W/m^2}]$ |
|---:|---:|---:|---:|
| $\frac{1}{10}$ | $0.1653$ | $16.53$ | $135.2$ |
| $\frac{1}{2}$ | $0.4444$ | $44.44$ | $29.63$ |
| $1$ | $0.5$ | $50$ | $0$ |
| $2$ | $0.4444$ | $44.44$ | $-7.407$ |
| $10$ | $0.1653$ | $16.53$ | $-1.352$ |

大きい抵抗と小さい抵抗で同じ発熱が起こりうる。二つの膜を区別するには、吸収電力とともに反射振幅を測ればよい。最大の近くではコンダクタンスが少しずれても吸収の変化は二次である。
:::

## 電磁場の運動量と Maxwell 応力テンソル

### 面を通じて働く力

電荷が受ける単位体積あたりの力は

$$
\vb{f}=\rho\vb{E}+\vb{J}\times\vb{B}
$$

である。
源を Maxwell 方程式から消去すると電場の項は $\varepsilon_0\vb{E}(\nabla\cdot\vb{E})$、磁場の項は $(\nabla\times\vb{B})\times\vb{B}/\mu_0-\varepsilon_0\partial_t\vb{E}\times\vb{B}$ になる。

ここへ $\partial_t(\varepsilon_0\vb{E}\times\vb{B})$ を足し、$\partial_t\vb{B}=-\nabla\times\vb{E}$ を代入する。
各場で $(\nabla\times\vb{X})\times\vb{X}=(\vb{X}\cdot\nabla)\vb{X}-\nabla X^2/2$ を使えば

$$
\begin{aligned}
f_i+\partial_tg_i&=\sum_j\partial_j T_{ij},\qquad \vb{g}=\varepsilon_0\vb{E}\times\vb{B}=\vb{S}/c^2,\\
T_{ij}&=\varepsilon_0(E_iE_j-\tfrac12\delta_{ij}E^2)+\mu_0^{-1}(B_iB_j-\tfrac12\delta_{ij}B^2).
\end{aligned}
\label{eq-em-stress}
$$

$T_{ij}$ を :term[Maxwell 応力テンソル]{en="Maxwell stress tensor" yomi="まくすうぇるおうりょくてんそる"} と呼ぶ。
成分は圧力の単位 $\mathrm{N/m^2}$ をもつ。
固定体積内の物質に働く力は $\oint\vb{T}\vb{n}\,\dd S-\dd\int\vb{g}\,\dd^3r/\dd t$ であり、外向き運動量流束はこの符号規約で $-T_{ij}$ となる。

物質を囲む面を真空中に選ぶと、物質内の運動量の分け方に依存せず全力を求められる。
対称な $T_{ij}$ では、位置ベクトルとの外積をとっても追加の体積トルクが残らず、角運動量密度は $\vb{r}\times\vb{g}$ である。

:::example[コンデンサー極板の引力]
導体面近くで電場が法線方向で磁場が零なら、真空側の $T_{nn}=\varepsilon_0E^2/2$ である。
導体内部の応力は零なので、表面は外向きにこの大きさで引かれる。
向かい合う二極板は互いに近づく向きに力を受ける。
:::

:::example[平面波が押す面]
真空の進行波は横波で $B=E/c$、$S=cu$ なので、運動量密度は $g=u/c$ である。
吸収面は入射エネルギー $W$ に対して $W/c$ を受け取り、静止した完全反射面は $2W/c$ を受け取る。
垂直入射の平均強度が $1.00\,\mathrm{kW/m^2}$ なら圧力はそれぞれ $3.34$ と $6.67\,\mathrm{\mu Pa}$ である。
:::

### 応力の積分で二電荷間の力を再現する

等しい正電荷 $q$ を $z=\pm a$ に置く。
上半空間を囲む面を赤道面と遠方の半球から作る。
遠方では $T_{ij}=O(r^{-4})$ なので半球の寄与は消える。
赤道面上で $E_z=0$、$E_s=2k_eqs/(s^2+a^2)^{3/2}$ であり、上半空間の外向き法線は $-\hat{\vb{z}}$ である。

従って面から受ける力は上向きで

$$
F_z=\int_0^\infty\frac{\varepsilon_0E_s^2}{2}\,2\pi s\,\dd s=4\pi\varepsilon_0 k_e^2q^2\int_0^\infty\frac{s^3\dd s}{(s^2+a^2)^3}.
$$

$u=s^2$ と置くと残る積分は $1/(4a^2)$ であり、$F_z=k_eq^2/(2a)^2$ を得る。
点電荷の位置における自己場を掛ける操作は一度も使っていない。
点電荷を半径の小さい球と考え、その球を含む半空間への力を求めてから半径を零にすれば、この表面積分の意味が明確になる。

### 運動量の流束と角運動量の流束

場の角運動量密度を $l_i=\epsilon_{ijk}x_jg_k$ とする。
運動量保存則の $k$ 成分に $\epsilon_{ijk}x_j$ を掛ける。
積の微分を展開すると

$$
\epsilon_{ijk}x_j\partial_m T_{km}=\partial_m(\epsilon_{ijk}x_jT_{km})-\epsilon_{imk}T_{km}.
$$

Maxwell 応力は対称なので最後の項は零になる。
従って物質に対するトルクの体積密度を

$$
\boldsymbol{\tau}=\vb{r}\times\vb{f}
$$

とすると

$$
\int_V\boldsymbol{\tau}\,\dd^3r=\oint_{\partial V}\vb{r}\times(\vb{T}\vb{n})\,\dd S-\frac{\dd}{\dd t}\int_V\vb{r}\times\vb{g}\,\dd^3r.
$$

二つの荷電物体への力を別々に計算しなくても、両方を含む面での応力からその全トルクを求められる。
原点を移すと角運動量もトルクも変わるが、力と運動量の保存則に対応する項が両側に加わるため、収支式の形は同じである。

### 斜めに入る光が鏡に及ぼす圧力

鏡の法線から角 $\theta$ で平面波が入射する。
入射強度 $I$ は光の進行方向に直角な単位面積を通る電力とする。
鏡の実面積 $A$ が受ける電力は射影により $IA\cos\theta$ である。
吸収なら単位エネルギーあたりの法線運動量は $\cos\theta/c$、鏡面反射ならこの成分が反転する。

従って静止した完全反射鏡では

$$
F_n/A=\frac{2I}{c}\cos^2\theta,\qquad F_t/A=0.
$$

接線力が零なのは、反射しても接線運動量が変わらないためである。
完全吸収面では $F_n/A=I\cos^2\theta/c$ に加え、入射方向に沿った接線成分 $F_t/A=I\sin\theta\cos\theta/c$ が残る。
吸収後の熱がどの方向へ再放射されるかは、ここで計算した入射光の寄与とは別に指定する必要がある。

$45^\circ$ の鏡面反射では同じ強度の垂直入射に比べて圧力は半分になる。
射影 $A\cos\theta$ を忘れて運動量の向きだけを変えると、誤って $\cos\theta$ 一個の式になる。

### 静止した電荷と磁場が持つ角運動量を回収する

電場も磁場も時間に依存しないときでも、両者が同じ場所に存在すれば電磁場に運動量がありうる。円筒形の装置では、対称性で全運動量は打ち消し合うものの、軸まわりの角運動量は残る。その値が機械的に測れることを、磁場を消す操作で確かめる。

無限に長い同軸の二つの薄い絶縁円筒に、単位長さ当たり $+\lambda$、$-\lambda$ を一様に固定する。半径は $a,b$ の順とする。その間に半径 $R$ の長い理想ソレノイドを配置し、$a<R<b$ とする。初めの磁束密度は $s<R$ で $B_0\hat{\vb{z}}$、外では零である。電場は $a<s<b$ でのみ $\lambda\hat{\vb{s}}/(2\pi\varepsilon_0s)$ となる。従って二つの場が重なる $a<s<R$ で、運動量密度は

$$
\vb{g}=\varepsilon_0\vb{E}\times\vb{B}
=-\frac{\lambda B_0}{2\pi s}\hat{\boldsymbol{\varphi}}
$$

である。円周一周の $\vb{g}$ の積分は零であるのに対し、$\vb{r}\times\vb{g}$ の軸方向成分は一周のどこでも同符号である。単位長さ当たりの角運動量は

$$
L_{z,\mathrm{field}}'=\int_a^R s g_\varphi\,2\pi s\,\dd s
=-\frac{\lambda B_0}{2}(R^2-a^2).
$$

次にソレノイド電流をゆっくり零にする。円筒を固定する軸受に伝わるトルクを測れば、円筒が動いて新しい磁場を作る効果も避けられる。任意の半径 $s$ での Faraday 則を時間積分すると

$$
\int E_\varphi(s,t)\,\dd t
=\frac{\Phi_B(s,\mathrm{initial})-\Phi_B(s,\mathrm{final})}{2\pi s}
=\frac{B_0}{2s}\min(s^2,R^2).
$$

従って、場が電荷に加えた角力積は、内筒では単位長さ当たり $+\lambda B_0 a^2/2$、外筒では $-\lambda B_0R^2/2$ である。合計は初めに電磁場に蓄えられていた $L_{z,\mathrm{field}}'$ に等しい。内筒と外筒はどちらも力を受けるが、外筒の位置で初めの磁場は零であった。外筒に働く力を生じたのは、その円周が囲む磁束の減少である。

:::example[角力積の符号を測る]
$a=1.00\,\mathrm{cm}$、$R=3.00\,\mathrm{cm}$、$\lambda=1.00\,\mu\mathrm{C/m}$、$B_0=0.200\,\mathrm T$ とする。内筒と外筒の単位長さ当たりの角力積は、それぞれ $+1.00\times10^{-11}$、$-9.00\times10^{-11}\,\mathrm{kg\,m/s}$ である。これを足した $-8.00\times10^{-11}\,\mathrm{kg\,m/s}$ は、場の角運動量の体積積分と一致する。$b$ を増やすと誘導電場は $1/b$ に比例して減るが、腕の長さは $b$ に比例して増すので合計は変わらない。
:::

### 薄膜の力を運動量流束から求める

一方向からの照射でも反射と吸収は別々に圧力に寄与する。入射強度を $I_i$ と置くと、反射波は逆向きに戻るため、単位面積・単位時間あたりに膜へ渡る運動量は $(1+|r|^2-|t|^2)I_i/c$ である。従って吸収係数 $A$ と反射率 $R$ を使えば $p=(A+2R)I_i/c$ と書ける。熱として吸収する割合が同じ二枚の膜でも、反射率が違えば支持台に加わる荷重が異なる。

この式を場だけから確かめるには、面に平行な $E_x,B_y$ に対して $T_{zz}=-(\varepsilon_0 E_x^2+B_y^2/\mu_0)/2$ であることを使う。左右に厚さの小さい箱をとり、上と下の面の法線の符号に注意して平均 $T_{zz}$ の差を計算する。入射と反射の交差項は電気と磁気で消え、圧力の式が得られる。シート内部に蓄えられる反応性エネルギーが時間で振動しても、周期平均では保存則の時間微分が消える。短パルスの途中で同じ差を用いる場合には、箱の中で変化する場の運動量を差し引かなければならない。

### 真空中の二つの波の局所応力

横電場を $\vb E=E_*\vb e$、$c\vb B=E_*\vb b$ とし、同じ周波数の複素ピーク振幅を用いる。実場は $\Re(\vb E e^{-i\omega t})$ とする。平均に用いる積は $\langle E_i E_j\rangle=E_*^2\Re(e_i e_j^*)/2$ である。$U_*=\varepsilon_0 E_*^2/2$ と書けば、

$$
\frac{\langle u\rangle}{U_*}=\frac{|\vb e|^2+|\vb b|^2}{2},\quad
\frac{\langle\vb S\rangle}{cU_*}=\Re(\vb e\times\vb b^*),\quad
\frac{\langle T_{ij}\rangle}{U_*}=\Re(e_i e_j^*+b_i b_j^*)-
\delta_{ij}\frac{\langle u\rangle}{U_*}.
$$

同じ観測点で重なる二本の波について、束と応力を具体的に比較する。各進行方向 $\vb n$ は単位ベクトルで、個々の波は $\vb b=\vb n\times\vb e$ を満たす。振幅を足してから積をとるので、異なる向きの光を独立な粒子流として足した値とは一般に異なる。以下の数値には $E_*=100\,\mathrm{V/m}$ を用いる。

:::example[同じ偏光の反射波が作る腹と節]
完全に反対向きの進行波を等振幅で重ねる。位相が合う位置では電場の腹、反対の位置では磁場の腹ができる。二つの進行方向、電場、およびそれぞれの波の磁場を列ベクトルで示す。

$$
\begin{aligned}
(\vb n_1,\vb n_2)&=\left(\begin{matrix}0 & 0\\0 & 0\\1 & -1\end{matrix}\right),\qquad (\vb e_1,\vb e_2)=\left(\begin{matrix}1 & 1\\0 & 0\\0 & 0\end{matrix}\right),\\
(\vb b_1,\vb b_2)&=\left(\begin{matrix}0 & 0\\1 & -1\\0 & 0\end{matrix}\right).
\end{aligned}
$$

列はこの順で対応する。第二波に対する相対位相を $\delta$、$C=\cos\delta$ と置く。電気と磁気の各エネルギーは、それぞれの二本の振幅の内積を交差項に含む。

$$
\begin{aligned}
\langle u_E\rangle/U_*&=C + 1,\quad \langle u_B\rangle/U_*=1 - C,\\
\langle\vb S\rangle/(cU_*)&=\left(\begin{matrix}0\\0\\0\end{matrix}\right).
\end{aligned}
$$

次に応力の対角成分からは全エネルギーを引き、非対角成分では向きの異なる成分の積を残す。法線 $\hat{\vb n}$ を指定した片面で測ると

$$
\begin{aligned}
\langle T\rangle/U_*&=\left(\begin{matrix}2 C & 0 & 0\\0 & - 2 C & 0\\0 & 0 & -2\end{matrix}\right),\\
\hat{\vb n}&=\left(\begin{matrix}0\\0\\1\end{matrix}\right),\qquad \langle T\rangle\hat{\vb n}/U_*=\left(\begin{matrix}0\\0\\-2\end{matrix}\right).
\end{aligned}
$$

この符号は外向き法線に対する応力である。物体に加わる全力は、この値を一つの開いた面で測るだけでは決まらず、物体を囲む他の面も含めて積分する。位相の異なる位置で法線への射影を物理単位で比べると次の通りである。

| 相対位相 | 全エネルギー密度 $[\mathrm{nJ/m^3}]$ | 法線電力流 $[\mathrm{W/m^2}]$ | 法線応力 $[\mathrm{nPa}]$ |
|---:|---:|---:|---:|
| $0^\circ$ | $88.54$ | $0$ | $-88.54$ |
| $60^\circ$ | $88.54$ | $0$ | $-88.54$ |
| $90^\circ$ | $88.54$ | $0$ | $-88.54$ |
| $120^\circ$ | $88.54$ | $0$ | $-88.54$ |
| $180^\circ$ | $88.54$ | $0$ | $-88.54$ |

閉じた共振器で平均電力が流れていないことは、壁に圧力が働かないことを意味しない。縦方向応力が相対位相によらない点がそれを示す。
:::

:::example[不完全な反射をもつ進行波]
右向きの波に振幅が半分の戻り波を加える。電場の包絡線には極小が現れるが、完全な節まで下がらない。二つの進行方向、電場、およびそれぞれの波の磁場を列ベクトルで示す。

$$
\begin{aligned}
(\vb n_1,\vb n_2)&=\left(\begin{matrix}0 & 0\\0 & 0\\1 & -1\end{matrix}\right),\qquad (\vb e_1,\vb e_2)=\left(\begin{matrix}1 & \frac{1}{2}\\0 & 0\\0 & 0\end{matrix}\right),\\
(\vb b_1,\vb b_2)&=\left(\begin{matrix}0 & 0\\1 & - \frac{1}{2}\\0 & 0\end{matrix}\right).
\end{aligned}
$$

列はこの順で対応する。第二波に対する相対位相を $\delta$、$C=\cos\delta$ と置く。電気と磁気の各エネルギーは、それぞれの二本の振幅の内積を交差項に含む。

$$
\begin{aligned}
\langle u_E\rangle/U_*&=\frac{4 C + 5}{8},\quad \langle u_B\rangle/U_*=- \frac{4 C - 5}{8},\\
\langle\vb S\rangle/(cU_*)&=\left(\begin{matrix}0\\0\\\frac{3}{4}\end{matrix}\right).
\end{aligned}
$$

次に応力の対角成分からは全エネルギーを引き、非対角成分では向きの異なる成分の積を残す。法線 $\hat{\vb n}$ を指定した片面で測ると

$$
\begin{aligned}
\langle T\rangle/U_*&=\left(\begin{matrix}C & 0 & 0\\0 & - C & 0\\0 & 0 & - \frac{5}{4}\end{matrix}\right),\\
\hat{\vb n}&=\left(\begin{matrix}0\\0\\1\end{matrix}\right),\qquad \langle T\rangle\hat{\vb n}/U_*=\left(\begin{matrix}0\\0\\- \frac{5}{4}\end{matrix}\right).
\end{aligned}
$$

面に垂直な流れと応力を、位相ごとに比較する。

| 相対位相 | 全エネルギー密度 $[\mathrm{nJ/m^3}]$ | 法線電力流 $[\mathrm{W/m^2}]$ | 法線応力 $[\mathrm{nPa}]$ |
|---:|---:|---:|---:|
| $0^\circ$ | $55.34$ | $9.954$ | $-55.34$ |
| $60^\circ$ | $55.34$ | $9.954$ | $-55.34$ |
| $90^\circ$ | $55.34$ | $9.954$ | $-55.34$ |
| $120^\circ$ | $55.34$ | $9.954$ | $-55.34$ |
| $180^\circ$ | $55.34$ | $9.954$ | $-55.34$ |

進む電力の差と縦応力を決める電力の和を区別できる。反射振幅が半分のとき戻る電力は四分の一である。
:::

:::example[直交する方向から同じ偏光で照射する]
直角に交わる二本の波の電場を、ともに入射面に垂直にする。二等分面を通るエネルギーと運動量は干渉によって場所ごとに変化する。二つの進行方向、電場、およびそれぞれの波の磁場を列ベクトルで示す。

$$
\begin{aligned}
(\vb n_1,\vb n_2)&=\left(\begin{matrix}0 & 1\\0 & 0\\1 & 0\end{matrix}\right),\qquad (\vb e_1,\vb e_2)=\left(\begin{matrix}0 & 0\\1 & 1\\0 & 0\end{matrix}\right),\\
(\vb b_1,\vb b_2)&=\left(\begin{matrix}-1 & 0\\0 & 0\\0 & 1\end{matrix}\right).
\end{aligned}
$$

列はこの順で対応する。第二波に対する相対位相を $\delta$、$C=\cos\delta$ と置く。電気と磁気の各エネルギーは、それぞれの二本の振幅の内積を交差項に含む。

$$
\begin{aligned}
\langle u_E\rangle/U_*&=C + 1,\quad \langle u_B\rangle/U_*=1,\\
\langle\vb S\rangle/(cU_*)&=\left(\begin{matrix}C + 1\\0\\C + 1\end{matrix}\right).
\end{aligned}
$$

次に応力の対角成分からは全エネルギーを引き、非対角成分では向きの異なる成分の積を残す。法線 $\hat{\vb n}$ を指定した片面で測ると

$$
\begin{aligned}
\langle T\rangle/U_*&=\left(\begin{matrix}- C - 1 & 0 & - C\\0 & C & 0\\- C & 0 & - C - 1\end{matrix}\right),\\
\hat{\vb n}&=\left(\begin{matrix}\frac{\sqrt{2}}{2}\\0\\\frac{\sqrt{2}}{2}\end{matrix}\right),\qquad \langle T\rangle\hat{\vb n}/U_*=\left(\begin{matrix}- \frac{\sqrt{2} \left(2 C + 1\right)}{2}\\0\\- \frac{\sqrt{2} \left(2 C + 1\right)}{2}\end{matrix}\right).
\end{aligned}
$$

面に垂直な流れと応力を、位相ごとに比較する。

| 相対位相 | 全エネルギー密度 $[\mathrm{nJ/m^3}]$ | 法線電力流 $[\mathrm{W/m^2}]$ | 法線応力 $[\mathrm{nPa}]$ |
|---:|---:|---:|---:|
| $0^\circ$ | $132.8$ | $37.54$ | $-132.8$ |
| $60^\circ$ | $110.7$ | $28.15$ | $-88.54$ |
| $90^\circ$ | $88.54$ | $18.77$ | $-44.27$ |
| $120^\circ$ | $66.41$ | $9.385$ | $0$ |
| $180^\circ$ | $44.27$ | $0$ | $44.27$ |

逆位相で電場が消える場所では磁場のエネルギーが残る。局所のエネルギー密度をそのまま光速倍して平均エネルギー流束に替えることはできない。
:::

:::example[入射面内に偏光する交差波]
今度は電場を二つの伝播方向を含む平面内に取り、磁場の向きをそろえる。斜めの検査面の法線を二つの伝播方向の差に沿って選ぶ。二つの進行方向、電場、およびそれぞれの波の磁場を列ベクトルで示す。

$$
\begin{aligned}
(\vb n_1,\vb n_2)&=\left(\begin{matrix}0 & 1\\0 & 0\\1 & 0\end{matrix}\right),\qquad (\vb e_1,\vb e_2)=\left(\begin{matrix}1 & 0\\0 & 0\\0 & -1\end{matrix}\right),\\
(\vb b_1,\vb b_2)&=\left(\begin{matrix}0 & 0\\1 & 1\\0 & 0\end{matrix}\right).
\end{aligned}
$$

列はこの順で対応する。第二波に対する相対位相を $\delta$、$C=\cos\delta$ と置く。電気と磁気の各エネルギーは、それぞれの二本の振幅の内積を交差項に含む。

$$
\begin{aligned}
\langle u_E\rangle/U_*&=1,\quad \langle u_B\rangle/U_*=C + 1,\\
\langle\vb S\rangle/(cU_*)&=\left(\begin{matrix}C + 1\\0\\C + 1\end{matrix}\right).
\end{aligned}
$$

次に応力の対角成分からは全エネルギーを引き、非対角成分では向きの異なる成分の積を残す。法線 $\hat{\vb n}$ を指定した片面で測ると

$$
\begin{aligned}
\langle T\rangle/U_*&=\left(\begin{matrix}- C - 1 & 0 & - C\\0 & C & 0\\- C & 0 & - C - 1\end{matrix}\right),\\
\hat{\vb n}&=\left(\begin{matrix}\frac{\sqrt{2}}{2}\\0\\- \frac{\sqrt{2}}{2}\end{matrix}\right),\qquad \langle T\rangle\hat{\vb n}/U_*=\left(\begin{matrix}- \frac{\sqrt{2}}{2}\\0\\\frac{\sqrt{2}}{2}\end{matrix}\right).
\end{aligned}
$$

面に垂直な流れと応力を、位相ごとに比較する。

| 相対位相 | 全エネルギー密度 $[\mathrm{nJ/m^3}]$ | 法線電力流 $[\mathrm{W/m^2}]$ | 法線応力 $[\mathrm{nPa}]$ |
|---:|---:|---:|---:|
| $0^\circ$ | $132.8$ | $0$ | $-44.27$ |
| $60^\circ$ | $110.7$ | $0$ | $-44.27$ |
| $90^\circ$ | $88.54$ | $0$ | $-44.27$ |
| $120^\circ$ | $66.41$ | $0$ | $-44.27$ |
| $180^\circ$ | $44.27$ | $0$ | $-44.27$ |

対称面に垂直なエネルギー輸送が零でも、そこでの応力ベクトルは一般に零ではない。向きの異なる運動量を運ぶ二つの流れが寄与する。
:::

:::example[互いに直交する偏光とせん断応力]
逆向きの波の電場を直交させる。電場の強度には干渉縞が現れないが、応力の非対角成分は偏光間の相関を検出する。二つの進行方向、電場、およびそれぞれの波の磁場を列ベクトルで示す。

$$
\begin{aligned}
(\vb n_1,\vb n_2)&=\left(\begin{matrix}0 & 0\\0 & 0\\1 & -1\end{matrix}\right),\qquad (\vb e_1,\vb e_2)=\left(\begin{matrix}1 & 0\\0 & 1\\0 & 0\end{matrix}\right),\\
(\vb b_1,\vb b_2)&=\left(\begin{matrix}0 & 1\\1 & 0\\0 & 0\end{matrix}\right).
\end{aligned}
$$

列はこの順で対応する。第二波に対する相対位相を $\delta$、$C=\cos\delta$ と置く。電気と磁気の各エネルギーは、それぞれの二本の振幅の内積を交差項に含む。

$$
\begin{aligned}
\langle u_E\rangle/U_*&=1,\quad \langle u_B\rangle/U_*=1,\\
\langle\vb S\rangle/(cU_*)&=\left(\begin{matrix}0\\0\\0\end{matrix}\right).
\end{aligned}
$$

次に応力の対角成分からは全エネルギーを引き、非対角成分では向きの異なる成分の積を残す。法線 $\hat{\vb n}$ を指定した片面で測ると

$$
\begin{aligned}
\langle T\rangle/U_*&=\left(\begin{matrix}0 & 2 C & 0\\2 C & 0 & 0\\0 & 0 & -2\end{matrix}\right),\\
\hat{\vb n}&=\left(\begin{matrix}1\\0\\0\end{matrix}\right),\qquad \langle T\rangle\hat{\vb n}/U_*=\left(\begin{matrix}0\\2 C\\0\end{matrix}\right).
\end{aligned}
$$

面に垂直な流れと応力を、位相ごとに比較する。

| 相対位相 | 全エネルギー密度 $[\mathrm{nJ/m^3}]$ | 法線電力流 $[\mathrm{W/m^2}]$ | 法線応力 $[\mathrm{nPa}]$ |
|---:|---:|---:|---:|
| $0^\circ$ | $88.54$ | $0$ | $0$ |
| $60^\circ$ | $88.54$ | $0$ | $0$ |
| $90^\circ$ | $88.54$ | $0$ | $0$ |
| $120^\circ$ | $88.54$ | $0$ | $0$ |
| $180^\circ$ | $88.54$ | $0$ | $0$ |

面の法線に直交する成分はせん断方向の運動量輸送を表す。向きを測らない強度検出器だけではこの項を区別できない。
:::

:::example[斜め入射する二光束と対称面]
同じ縦方向速度成分を持つ二光束を左右から重ねる。対称面では横向きの入射運動量が相殺するが、光束間の横方向の位相差が局所強度を変える。二つの進行方向、電場、およびそれぞれの波の磁場を列ベクトルで示す。

$$
\begin{aligned}
(\vb n_1,\vb n_2)&=\left(\begin{matrix}\frac{3}{5} & - \frac{3}{5}\\0 & 0\\\frac{4}{5} & \frac{4}{5}\end{matrix}\right),\qquad (\vb e_1,\vb e_2)=\left(\begin{matrix}0 & 0\\1 & 1\\0 & 0\end{matrix}\right),\\
(\vb b_1,\vb b_2)&=\left(\begin{matrix}- \frac{4}{5} & - \frac{4}{5}\\0 & 0\\\frac{3}{5} & - \frac{3}{5}\end{matrix}\right).
\end{aligned}
$$

列はこの順で対応する。第二波に対する相対位相を $\delta$、$C=\cos\delta$ と置く。電気と磁気の各エネルギーは、それぞれの二本の振幅の内積を交差項に含む。

$$
\begin{aligned}
\langle u_E\rangle/U_*&=C + 1,\quad \langle u_B\rangle/U_*=\frac{7 C + 25}{25},\\
\langle\vb S\rangle/(cU_*)&=\left(\begin{matrix}0\\0\\\frac{8 \left(C + 1\right)}{5}\end{matrix}\right).
\end{aligned}
$$

次に応力の対角成分からは全エネルギーを引き、非対角成分では向きの異なる成分の積を残す。法線 $\hat{\vb n}$ を指定した片面で測ると

$$
\begin{aligned}
\langle T\rangle/U_*&=\left(\begin{matrix}- \frac{18}{25} & 0 & 0\\0 & \frac{18 C}{25} & 0\\0 & 0 & - \frac{2 \left(25 C + 16\right)}{25}\end{matrix}\right),\\
\hat{\vb n}&=\left(\begin{matrix}0\\0\\1\end{matrix}\right),\qquad \langle T\rangle\hat{\vb n}/U_*=\left(\begin{matrix}0\\0\\- \frac{2 \left(25 C + 16\right)}{25}\end{matrix}\right).
\end{aligned}
$$

面に垂直な流れと応力を、位相ごとに比較する。

| 相対位相 | 全エネルギー密度 $[\mathrm{nJ/m^3}]$ | 法線電力流 $[\mathrm{W/m^2}]$ | 法線応力 $[\mathrm{nPa}]$ |
|---:|---:|---:|---:|
| $0^\circ$ | $145.2$ | $42.47$ | $-145.2$ |
| $60^\circ$ | $116.9$ | $31.85$ | $-100.9$ |
| $90^\circ$ | $88.54$ | $21.24$ | $-56.67$ |
| $120^\circ$ | $60.21$ | $10.62$ | $-12.4$ |
| $180^\circ$ | $31.88$ | $0$ | $31.88$ |

広い検出器で多くの干渉周期を平均すると位相に比例する項が落ちる。一点での局所応力と面全体で積分した荷重は、ここで初めて結び付く。
:::


:::example[反射損失を含む光圧]
透過のない面に真空から垂直に光を照射し、エネルギー反射率を $x$ とする。吸収光は $I(1-x)/c$、反射光は $2Ix/c$ だけ単位面積あたりの力に寄与する。面からの熱放射が後で運び去る運動量はこの瞬間の計算には含めない。

入射強度 $I=1.00\,\mathrm{kW/m^2}$、$Y$ は光圧、$Y_*=I/c=3.336\,\mathrm{\mu Pa}$ とする。

$$
f(x)=x + 1,\qquad Y=Y_* f(x).
$$

二点の読みの差を求めるには、基準を同じに保ち、振幅の換算を最後に一度だけ行う。最初と最後の設定で

$$
\Delta Y/Y_* = f\!\left(1\right)-f\!\left(0\right)=1.
$$

第四列は初期設定からの変化である。

| $x$ | $f(x)$ | $Y$ $[\mathrm{\mu Pa}]$ | $Y-Y(x_0)$ $[\mathrm{\mu Pa}]$ |
|---:|---:|---:|---:|
| $0$ | $1$ | $3.336$ | $0$ |
| $\frac{1}{4}$ | $1.25$ | $4.17$ | $0.8339$ |
| $\frac{1}{2}$ | $1.5$ | $5.003$ | $1.668$ |
| $\frac{9}{10}$ | $1.9$ | $6.338$ | $3.002$ |
| $1$ | $2$ | $6.671$ | $3.336$ |

反射できない成分も面で止まるため、反射率を零にしても荷重は零にならない。同じ実験を斜入射に替えると、法線単位面積に当たる光量と法線運動量への射影の両方に余弦が入る。
:::

### 光パルスが箱を通過する間の力積

真空中を $+z$ 向きに進む断面積 $A$ の平面波パルスを、$0<z<L$ の固定した箱で追う。横方向の端の効果を無視し、入射強度を時刻の任意関数 $I(t)$ とする。右端での強度は $I(t-L/c)$ である。物質を入れていなくても、入射開始直後は左端の応力だけが非零である。これを内部物質の力と解釈しないために、三つの項を個別に求める。

局所の運動量密度は $g_z=I(t-z/c)/c^2$ であり、箱に蓄えられた全運動量は

$$
P_{\rm field}(t)=\frac A{c^2}\int_0^L I(t-z/c)\dd z
=\frac A c\int_{t-L/c}^{t}I(s)\dd s.
$$

時間微分すると $\dot P_{\rm field}=(A/c)[I(t)-I(t-L/c)]$ となる。応力は $T_{zz}=-I/c$、左端の法線が $-\hat{\vb z}$、右端が $+\hat{\vb z}$ なので、応力の全表面積分も同じ値である。従って物体への力は任意の波形の任意の時刻で厳密に零となる。パルスが箱より十分短く、全体が一時箱に収まるなら、その瞬間の運動量はパルスエネルギーを $c$ で割った値である。

右端に完全吸収体を置けば、物質へ渡る力が $AI(t-L/c)/c$ だけ現れる。反射がなく、熱放射が観測時間内に無視できる場合には、時間積分した力積は $U_{\rm pulse}/c$ である。負荷に届く前には左端から運動量が入っていても力はまだ生じない。定常波の圧力を時間依存する場合に適用するときは、照射面での遅延した強度を用いる必要がある。

```tikz caption="固定した検査体積をパルスが横切る。左右の面の外向き法線は逆向きである。"
\begin{tikzpicture}
\draw[dashed] (0,-.7) rectangle (4,1.0);
\draw[->] (-.1,.76)--(-.7,.76) node[left] {$-\hat z$};
\draw[->] (4.1,.76)--(4.7,.76) node[right] {$+\hat z$};
\draw[thick] (-1.3,-.35)--(.5,-.35)--(1.1,.5)--(1.7,-.35)--(4.6,-.35);
\draw[->] (1.8,.3)--(2.8,.3) node[midway,above] {$c$};
\node[below] at (0,-.7) {$0$};
\node[below] at (4,-.7) {$L$};
\end{tikzpicture}
```

## 演習問題

::::exercise{level=basic}
容量 $100\,\mathrm{pF}$ の素子に $2.0\times10^6\,\mathrm{V/s}$ で昇圧する。
変位電流を求めよ。

:::solution
電極間を通る変位電流は $\dot Q=C\dot V$ であり、

$$
100\times10^{-12}\times2.0\times10^6=0.200\,\mathrm{mA}
$$

となる。
:::
::::

::::exercise{level=basic}
一様な $300\,\mathrm{V/m}$ の静電場の真空エネルギー密度を求めよ。

:::solution
磁場がないので

$$
u=\varepsilon_0E^2/2=3.98\times10^{-7}\,\mathrm{J/m^3}
$$

である。
:::
::::

::::exercise{level=basic}
真空の平面波が $E_{\rm rms}=10.0\,\mathrm{V/m}$ のとき平均強度を求めよ。

:::solution
$B=E/c$ より

$$
\langle S\rangle=\varepsilon_0c\langle E^2\rangle
$$

である。
実効値は二乗平均の平方根なので $0.265\,\mathrm{W/m^2}$ となる。
:::
::::

::::exercise{level=basic}
$10.0\,\mathrm W$ の光が垂直入射して完全吸収される。
力を求めよ。

:::solution
毎秒吸収されるエネルギーは $10.0\,\mathrm J$ である。
進行波の運動量はエネルギーを $c$ で割った値なので、力は

$$
10.0/c=3.34\times10^{-8}\,\mathrm N
$$

である。
:::
::::

::::exercise{level=advanced}
真空同軸線の内筒半径 $a$、外筒半径 $b$ に電位差 $V$ と電流 $I$ がある。
断面を流れる電力を求めよ。

:::solution
$E_s=V/[s\ln(b/a)]$、$H_\varphi=I/(2\pi s)$ である。
向きを合わせれば $S_z=VI/[2\pi s^2\ln(b/a)]$、環状断面の積分は $\int_a^b S_z2\pi s\dd s=VI$ となる。
:::
::::

::::exercise{level=advanced}
内部に電荷・電流がなく完全導体に囲まれた空洞で、同じ初期値からの二つの Maxwell 解が一致することを示せ。

:::solution
差の場にも真空の保存則が成立する。
壁では電場の接線成分が零なので、差の Poynting 流も法線成分が零である。
従って差の場の全エネルギーは初期値の零を保つ。
非負な電磁場の二乗積分が零なので、正則な差の場は至る所で零である。
:::
::::

::::exercise{level=advanced}
真空の静磁場 $\vb{B}=B\hat{\vb{z}}$ の応力の主値を求めよ。

:::solution
磁気項の外積行列の対角は $(0,0,B^2)$ で、等方項 $B^2/2$ を各対角から引く。
従って固有値は $(-B^2,-B^2,B^2)/(2\mu_0)$ である。
軸方向に張力、垂直に圧力が働く符号である。
:::
::::

::::exercise{level=advanced}
真空の任意の場について $|\vb{S}|\le cu$ を証明し、等号の条件を述べよ。

:::solution
$|\vb{E}\times\vb{B}|\le EB$ と $2(\sqrt{\varepsilon_0}E)(B/\sqrt{\mu_0})\le\varepsilon_0E^2+B^2/\mu_0$ を掛け合わせると主張を得る。
非零場での等号には場が直交し、$E=cB$ であることが必要十分である。
:::
::::

::::exercise{level=advanced}
半径 $a,b$、長さ $\ell$ の長い同軸円筒をコンデンサーとする。
内筒に $+Q$、外筒に $-Q$ があり、隙間に一様な静磁場 $B_0\hat{\vb{z}}$ がある。
端効果を無視してこの隙間の電磁場の軸方向角運動量を求めよ。

:::solution
電場は

$$
E_s=Q/(2\pi\varepsilon_0\ell s)
$$

である。
運動量密度は

$$
\vb{g}=\varepsilon_0\vb{E}\times\vb{B}=-QB_0\hat{\vb{\varphi}}/(2\pi\ell s)
$$

となる。

従って角運動量密度の軸成分 $s g_\varphi$ は一定で、

$$
L_z=\int_a^b(sg_\varphi)\,2\pi s\ell\,\dd s=-\frac{QB_0}{2}(b^2-a^2).
$$

磁場をゆっくり零へ下げる場合を考える。
両円筒を囲む磁束の差は $\pi(b^2-a^2)B$ である。
Faraday 則から両表面へのトルクの和は

$$
Q(aE_\varphi(a)-bE_\varphi(b))=Q(b^2-a^2)\dot B/2
$$

である。

時間積分するとちょうど初期の場の角運動量が物質へ移る。
:::
::::

::::exercise{level=advanced}
真空中で、エネルギー密度 $u$ と運動量密度 $\vb{g}$ について $u^2-c^2g^2$ を $\vb{E}^2,\vb{B}^2,\vb{E}\cdot\vb{B}$ だけで表せ。
どの場で零になるか。

:::solution
$\vb{g}=\varepsilon_0\vb{E}\times\vb{B}$ と $|\vb{E}\times\vb{B}|^2=E^2B^2-(\vb{E}\cdot\vb{B})^2$ を用いて展開すると

$$
u^2-c^2g^2=\frac14\left(\varepsilon_0E^2-\frac{B^2}{\mu_0}\right)^2+\frac{\varepsilon_0}{\mu_0}(\vb{E}\cdot\vb{B})^2.
$$

右辺は二つの平方和である。
非零場では $\vb{E},\vb{B}$ が直交して同時に $E=cB$ のときに限り零である。
一方向の真空平面波はこれを満たす。
逆向きに進む二波の合成はそれぞれが平面波でも、一般にはこの条件を満たさない。
:::
::::

::::exercise{level=advanced}
電場のピーク振幅 $E_0$ の真空平面波が完全導体面に垂直入射する。
入射・反射を重ねた壁直前の瞬時応力から平均圧力を求めよ。

:::solution
壁は接線電場を零にするので、壁直前の入射波と反射波の電場は打ち消し合う。
磁場は進行方向の反転も受けるため、壁での接線磁場は $B_t=2(E_0/c)\cos\omega t$ と書ける。
法線成分は零なので $T_{nn}=-B_t^2/(2\mu_0)$ である。

導体へ向かう圧力の大きさは

$$
\overline p=\frac1{2\mu_0}\left\langle\frac{4E_0^2}{c^2}\cos^2\omega t\right\rangle=\varepsilon_0 E_0^2=\frac{2I_{\rm in}}c.
$$

壁では平均 Poynting 流が零でも、応力の平均は零にならない。
反射波が運動量の向きを変えるためである。
:::
::::

::::exercise{level=advanced}
幅の広いレーザーパルスが真空中の吸収板に垂直に入射する。全エネルギーは $3.00\,\mathrm{mJ}$ で、板はエネルギーの $20.0\%$ を正反射し、残りを吸収する。板に渡る力積を求めよ。

:::solution
入射時の法線運動量は $U/c$、外へ戻る光は $-0.200U/c$ を持つ。透過はなく、反射終了後には箱内の場も消える。従って板への力積は差の $1.20 U/c=1.20\times10^{-11}\,\mathrm{N\,s}$ であり、板の内側へ向かう。吸収された部分と反射部分で分けても、$(0.800+2\times0.200)U/c$ と一致する。
:::
::::

## 確認テスト

::::quiz{answer="2" tol="0.01"}
$0.50\,\mathrm A$ で充電中の半径 $5.0\,\mathrm{cm}$ の理想円形極板の縁の磁場を $\mathrm{\mu T}$ で答えよ。

:::explanation
周長 $2\pi a$ に沿う磁場を積分すると $\mu_0I$ になるため、

$$
B=\mu_0I/(2\pi a)=2.00\,\mathrm{\mu T}
$$

である。
:::
::::

::::quiz
真空で場から電荷へ渡る単位体積・単位時間あたりのエネルギーを選べ。

- [x] $\vb{J}\cdot\vb{E}$
- [ ] $\vb{J}\cdot\vb{B}$
- [ ] $\nabla\cdot\vb{E}$

:::explanation
粒子の仕事率 $q\vb{v}\cdot(\vb{E}+\vb{v}\times\vb{B})$ の磁気項は零である。
体積あたりで足すと $\vb{J}\cdot\vb{E}$ になる。
:::
::::

::::quiz{answer="2" tol="1e-8"}
同じ強度の光が垂直に入射するとき、完全反射面の圧力は完全吸収面の何倍か。

:::explanation
反射で法線運動量が正から負へ反転するので運動量変化は吸収の場合の二倍である。
:::
::::

::::quiz
静電場と静磁場が直交して共存するときに成立するものを選べ。

- [ ] 静的なので運動量密度は零である。
- [x] 運動量密度は $\varepsilon_0\vb{E}\times\vb{B}$ である。
- [ ] 応力テンソルは反対称である。

:::explanation
運動量密度には場の時間微分が含まれない。
直交した非零の電場と磁場は、時間一定でも非零の運動量密度をもつ。
:::
::::

::::quiz
真空中の閉曲面で Maxwell 応力を積分した量が、そのまま内部物質への全力になる条件を選べ。

- [x] 曲面内部の場の全運動量が時間的に一定である。
- [ ] 閉曲面内に自由電荷があるだけで十分である。
- [ ] 電場と磁場が直交すれば常に成立する。

:::explanation
物質への力は $\oint\vb{T}\vb{n}\,\dd S-\partial_t\int\vb{g}\,\dd^3r$ である。
第二項が零であれば応力の積分だけが残る。
静的な場はその十分条件の一つである。
:::
::::

::::quiz{answer="3.34" tol="0.01"}
$1.00\,\mathrm{kW/m^2}$ の真空進行波を垂直に完全吸収する面の圧力を $\mathrm{\mu Pa}$ で答えよ。

:::explanation
エネルギーとともに運ばれる法線運動量はエネルギーを $c$ で割った量である。
従って圧力は

$$
1000/c=3.34\times10^{-6}\,\mathrm{Pa}
$$

である。
:::
::::

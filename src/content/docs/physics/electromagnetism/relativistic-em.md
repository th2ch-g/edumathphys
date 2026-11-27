---
title: 電磁気学と特殊相対論
description: 時空の計量と四元ベクトルを定義し、電場と磁場の変換、電磁場の不変量、相対論的な荷電粒子の軌道を導く。
notation:
  - symbol: '\eta_{\mu\nu}'
    meaning: Minkowski 計量
  - symbol: 'F^{\mu\nu}'
    meaning: 電磁場テンソル
  - symbol: 'J^\mu'
    meaning: 四元電流
  - symbol: 'U^\mu'
    meaning: 四元速度
  - symbol: 'T^{\mu\nu}'
    meaning: エネルギー運動量テンソル
---

## この章の目標

- 時間成分と空間成分の単位をそろえ、添字を上げ下げする。
- 共変な二組の方程式から Maxwell の四式を取り出す。
- 電磁場の変換と不変量を用いて場を分類する。
- 四元 Lorentz 力から一様な電場・磁場中の運動を解く。

## 電磁場テンソル

### 事象と慣性系

電荷とともに動く観測者は純粋な静電場を見る。
同じ電荷を通り過ぎる観測者は電流も測るので、磁場を含めて同じ現象を記述する必要がある。
そのために、空間座標だけでなく測定時刻も変換する。
四元ベクトルの変換則には :ref[def-sr-fourvector] の規約を用いる。以下では原点を共通とする二つの慣性系を考える。

事象の四座標は $x^\mu=(ct,x,y,z)$、計量は

$$
\eta_{\mu\nu}=\operatorname{diag}(-1,+1,+1,+1)
$$

とする。
ギリシャ文字の添字は $0,1,2,3$、ラテン文字の添字は $1,2,3$ をとる。

上下に一度ずつ現れる添字は和をとる。
四元ベクトルの下付き成分は $a_\mu=\eta_{\mu\nu}a^\nu=(-a^0,\vb{a})$、内積は

$$
a_\mu b^\mu=-a^0b^0+\vb{a}\cdot\vb{b}
$$

である。

微分は $\partial_\mu=(c^{-1}\partial_t,\nabla)$、

$$
\partial^\mu=(-c^{-1}\partial_t,\nabla)
$$

となる。

慣性系 $S'$ が $S$ から見て $+x$ に速度 $v$ で動く受動変換を採る。
$\beta=v/c$、$\gamma=(1-\beta^2)^{-1/2}$ とすれば

$$
x'^0=\gamma(x^0-\beta x^1),\qquad x'^1=\gamma(x^1-\beta x^0),\qquad x'^2=x^2,\quad x'^3=x^3.
\label{eq-em-boost}
$$

事象の差に代入すると

$$
-(\Delta x'^0)^2+(\Delta x'^1)^2=-(\Delta x^0)^2+(\Delta x^1)^2
$$

であり、時空間隔が不変である。
この計量を保存する行列を $\Lambda$ と書き、四元ベクトルは $a'^\mu=\Lambda^\mu{}_\nu a^\nu$ で変換する。

質量のある粒子に沿う固有時間は

$$
\dd\tau=\sqrt{-\eta_{\mu\nu}\dd x^\mu\dd x^\nu}/c=\dd t/\gamma_u
$$

で、$\gamma_u$ は粒子の速さ $u$ に対応する。

### ポテンシャルから反対称テンソルを作る

四元ポテンシャルを $A^\mu=(\phi/c,\vb{A})$ とする。
時間成分の単位は他の三成分と同じ $\mathrm{T\,m}$ である。
Lorenz 条件は $\partial_\mu A^\mu=0$ となる。
二階テンソルでは各上付き添字が独立に行列 $\Lambda$ で変換する。

この変換に従う :term[電磁場テンソル]{en="electromagnetic field tensor" yomi="でんじばてんそる"} を

$$
F^{\mu\nu}=\partial^\mu A^\nu-\partial^\nu A^\mu,\qquad F^{0i}=E_i/c,\qquad F^{ij}=\epsilon_{ijk} B_k
\label{eq-em-field-tensor}
$$

と定める。
ここで $\epsilon_{123}=+1$ であり、空間添字では上下の数値は同じである。
例えば $F^{01}=-(\partial_t A_x+\partial_x\phi)/c$、

$$
F^{12}=\partial_xA_y-\partial_yA_x=B_z
$$

である。

行列にすれば

$$
F^{\mu\nu}=\begin{pmatrix}
0&E_x/c&E_y/c&E_z/c\\
-E_x/c&0&B_z&-B_y\\
-E_y/c&-B_z&0&B_x\\
-E_z/c&B_y&-B_x&0
\end{pmatrix}.
$$

全成分は $\mathrm T$ の次元をもつ。
両添字を下げると時間添字を一個含む成分の符号だけが変わる。
$A^\mu\to A^\mu+\partial^\mu\chi$ では二階偏微分が交換するため $F^{\mu\nu}$ は不変である。

:::example[一様電場を二つの四元ポテンシャルで表す]
$\phi=-E_0x,\ \vb{A}=0$ と $\phi'=0,\ \vb{A}'=-E_0t\hat{\vb{x}}$ は $\chi=-E_0xt$ で移り合う。
どちらでも $F^{01}=E_0/c$、$F^{10}=-E_0/c$ で、他の成分は零である。

$F_{01}$ は $-E_0/c$ となる。
:::

### なぜ力を表す行列が反対称なのか

粒子の四元速度に線形な局所的な力 $K^\mu=qG^{\mu\nu}U_\nu$ を考える。
静止質量が一定なら $\dd(U_\mu U^\mu)/\dd\tau=0$ であるから、$U_\mu G^{\mu\nu}U_\nu=0$ がすべての時間的速度について成立しなければならない。

任意の時間的ベクトルは四元速度の定数倍なので、この二次式は時間的ベクトルの開集合で零である。
二次多項式が開集合で恒等的に零ならすべての係数が零であるため、$G$ の対称部分は消える。
そこで残る六つの独立成分に、静止電荷の力の三成分と速度に直交する磁気力の三成分を対応させると :eqref[eq-em-field-tensor] の成分が定まる。

:::example[観測者が測る電場を四元量から取り出す]
四元速度 $W^\mu$ の観測者に対して $e^\mu=F^{\mu\nu}W_\nu$ と置く。
反対称性から $e^\mu W_\mu=0$ であり、その観測者の静止系では $W_\nu=(-c,0,0,0)$ なので $e^\mu=(0,\vb{E})$ となる。

向きを $\epsilon^{0123}=+1$ とし、$\displaystyle {}^*F^{\mu\nu}=\epsilon^{\mu\nu\alpha\beta}F_{\alpha\beta}/2$ と定義する。計量の時間成分が負なので $\epsilon_{0123}=-1$ である。この双対テンソルを用いて $b^\mu={}^*F^{\mu\nu}W_\nu/c$ とすると静止系で $(0,\vb{B})$ となる。
後者では空間の向きから磁束密度の符号も決まる。
この計算では、誰の時計・定規で測る場かが $W^\mu$ によって明示される。
:::

### 異なる観測者が場を分解しても同じテンソルに戻ること

観測者の四元速度を $U^\mu U_\mu=-c^2$ とし、$e^\mu=F^{\mu\nu}U_\nu$、$b^\mu={}^*F^{\mu\nu}U_\nu/c$ を作る。本章の符号規約では、静止している観測者について $e^\mu=(0,\vb E)$、$b^\mu=(0,\vb B)$ である。反対称性により $e^\mu U_\mu=b^\mu U_\mu=0$ が恒等的に成り立つ。従ってこの二本に含まれる独立な量はそれぞれ三成分で、合わせて元の六成分と同じである。逆に

$$
F^{\mu\nu}=\frac{U^\mu e^\nu-U^\nu e^\mu}{c^2}
-\frac1c\epsilon^{\mu\nu\rho\sigma}U_\rho b_\sigma
\label{eq-em-observer-reconstruction}
$$

と再構成できる。観測者の静止系では第一項から $F^{0i}=E_i/c$ が出る。第二項では $U_0=-c$ のみが残り $F^{ij}=\epsilon^{ij0k}B_k=\epsilon_{ijk}B_k$ となる。両辺がテンソルで一つの慣性系で全成分が一致するので、任意の慣性系で同じである。この証明は、異なる慣性系の $\vb E$ と $\vb B$ から一本ずつ選んでテンソルを作れない理由も示す。二本とも同じ観測者の空間へ射影した量でなければならない。


:::example[同じ固有電荷を数える測定器]
一様な電荷雲の静止系密度を $\rho_0$ とする。体積素片と共に動く系では電流は零であり、四元電流を変換すると実験室の密度は $\gamma\rho_0$、電流密度は $\gamma\rho_0\vb v$ になる。

$x=v/c$、$Y=\rho$、$Y_*=\rho_0=1.00\,\mathrm{nC/m^3}$ とする。有限の切片の両端は実験室系で同時に測るものとする。

$$
f(x)=\frac{1}{\sqrt{1 - x} \sqrt{x + 1}},\qquad Y=Y_* f(x).
$$

変数を少し動かしたときの読みの変化は、接線の傾きから計算できる。基準値 $x_0$ の近くで

$$
\begin{aligned}
\frac{\Delta Y}{Y_*}= f\!\left(x_0+\Delta x\right)-f(x_0)= f'(x_0)\Delta x+O((\Delta x)^2),\\
f'(x)=\frac{x}{\left(1 - x\right)^{\frac{3}{2}} \left(x + 1\right)^{\frac{3}{2}}}.
\end{aligned}
$$

同じ尺度で値と傾きを並べる。

| $x$ | $f(x)$ | $Y$ $[\mathrm{nC/m^3}]$ | $\dd Y/\dd x$ $[\mathrm{nC/m^3}]$ |
|---:|---:|---:|---:|
| $0$ | $1$ | $1$ | $0$ |
| $\frac{3}{10}$ | $1.048$ | $1.048$ | $0.3456$ |
| $\frac{3}{5}$ | $1.25$ | $1.25$ | $1.172$ |
| $\frac{4}{5}$ | $1.667$ | $1.667$ | $3.704$ |

密度の増加とともに進行方向の長さが逆比で縮むため、同じ粒子群を包む総電荷は変わらない。等しい座標体積を両系で形式的に掛けると、異なる粒子群を数えてしまう。
:::

## Maxwell 方程式の共変形

### 四元電流と二組の方程式

電荷 $q$ の世界線 $z^\mu(\tau)$ の接線を :term[四元速度]{en="four-velocity" yomi="よんげんそくど"} と呼び、

$$
U^\mu=\dd z^\mu/\dd\tau=\gamma_u(c,\vb{u})
$$

とする。
:term[四元電流]{en="four-current" yomi="よんげんでんりゅう"} を $J^\mu(x)=qc\int U^\mu\delta^4(x-z(\tau))\,\dd\tau$ と書くと、時間のデルタを積分して $J^\mu=(c\rho,\vb{J})$ が得られる。

$\delta^4$ は Lorentz 変換の行列式の絶対値が $1$ なので不変であり、従って $J^\mu$ は四元ベクトルである。
個々の粒子の電荷が一定なら $\partial_\mu J^\mu=0$、すなわち :eqref[eq-em-continuity] である。

:::law[Maxwell 方程式の共変形]{#law-em-covariant}
上で定めた $F^{\mu\nu}$ では、真空の全電荷・全電流を源として

$$
\partial_\mu F^{\nu\mu}=\mu_0J^\nu,\qquad \partial_\lambda F_{\mu\nu}+\partial_\mu F_{\nu\lambda}+\partial_\nu F_{\lambda\mu}=0.
\label{eq-em-covariant}
$$

:::

第一式で $\nu=0$ なら $\nabla\cdot\vb{E}/c=\mu_0 c\rho$、$\nu=i$ なら

$$
-\partial_t E_i/c^2+(\nabla\times\vb{B})_i=\mu_0J_i
$$

となる。
第二式で $(\lambda,\mu,\nu)=(1,2,3)$ を選ぶと $\nabla\cdot\vb{B}=0$、$(0,1,2)$ を選ぶと

$$
\partial_tB_z+\partial_xE_y-\partial_yE_x=0
$$

が得られる。

他の二組の時間・空間添字は巡回置換で得る。
第二式は $F$ をポテンシャルで表せば二階微分が対消滅する恒等式である。
第一式の発散をとると微分は対称、$F$ は反対称なので $\partial_\nu J^\nu=0$ が必要となる。

四次元の完全反対称記号を $\epsilon^{0123}=+1$ とし、双対を

$$
\displaystyle {}^*F^{\mu\nu}=\epsilon^{\mu\nu\alpha\beta}F_{\alpha\beta}/2
$$

と置く。

すると

$$
\displaystyle {}^*F^{0i}=B_i,\ {}^*F^{ij}=-\epsilon_{ijk}E_k/c
$$

である。
第二式は $\partial_\mu{}^*F^{\nu\mu}=0$ に等しい。
Lorenz ゲージで第一式を $A^\mu$ で書けば $\partial_\alpha\partial^\alpha A^\nu=-\mu_0J^\nu$、

$$
\partial_\alpha\partial^\alpha=\nabla^2-c^{-2}\partial_t^2
$$

となり、:eqref[eq-em-lorenz-wave] の時間・空間成分をまとめられる。

### エネルギーと運動量を合わせる

共変な量

$$
T^{\mu\nu}=\frac1{\mu_0}\left(F^{\mu\alpha}F^\nu{}_{\alpha}-\tfrac14\eta^{\mu\nu}F_{\alpha\beta}F^{\alpha\beta}\right)
\label{eq-em-energy-momentum}
$$

は対称である。
成分を入れると $T^{00}=u$、$T^{0i}=T^{i0}=S_i/c$、$T^{ij}$ は :eqref[eq-em-stress] の Maxwell 応力の負となる。
$f^\nu=F^{\nu\lambda}J_\lambda$ は $(\vb{J}\cdot\vb{E}/c,\rho\vb{E}+\vb{J}\times\vb{B})$ の成分をもつ。

従って Poynting の定理と運動量収支を並べれば

$$
\partial_\mu T^{\mu\nu}=-f^\nu
$$

が成り立つ。
$\eta_{\mu\nu}T^{\mu\nu}=0$ は二項の縮約が等しいためで、電磁場のエネルギー密度が零という意味ではない。

:::example[単独の平面波]
$+z$ へ進む波では $E=cB$、$S_z=cu$ である。
非零成分は $T^{00}=T^{03}=T^{30}=T^{33}=u$ の四つとなる。
波のエネルギーが $c$ で運ばれ、法線運動量流束がエネルギー密度と等しいことを一つの行列から読み取れる。
:::

### 作用から場の源の符号を確認する

時空体積を $\dd^4x=c\,\dd t\,\dd^3r$ とする。
外から与えた保存電流に対する場の作用を

$$
\mathcal S[A]=\frac1c\int\left[-\frac1{4\mu_0}F_{\mu\nu}F^{\mu\nu}+J^\mu A_\mu\right]\dd^4x
$$

と書く。
後項は空間と時間で書くと $\vb{J}\cdot\vb{A}-\rho\phi$ である。
境界で消える $\delta A_\nu$ をとると、反対称性により場の項の変分は $-F^{\mu\nu}\partial_\mu\delta A_\nu/\mu_0$ になる。

部分積分すると変分の係数は $\partial_\mu F^{\mu\nu}/\mu_0+J^\nu$ なので、作用が停留する条件は :eqref[eq-em-covariant] と同じになる。
ゲージ変換による源の作用の変化は $\int J^\mu\partial_\mu\chi\,\dd^4x/c$ である。

境界項を除けば

$$
-\int\chi\partial_\mu J^\mu\,\dd^4x/c=0
$$

となり、源との結合のゲージ不変性と電荷保存が結びつく。

:::example[電流密度が増えても全電荷は変わらない]
体積 $V_0=A L_0$ の短い直方体が、静止系で一様な電荷密度 $\rho_0$ をもつ。
長辺方向へ速さ $v$ の系から同時刻に両端を測ると長さは $L_0/\gamma$、断面積は $A$ である。
一方、$J^\mu=(c\rho_0,\vb 0)$ の時間成分は $c\gamma\rho_0$ になる。

従って

$$
\rho' V'=\gamma\rho_0V_0/\gamma=\rho_0V_0
$$

であり、測定する同時面が変わっても全電荷は同じである。
:::

### 点電荷のポテンシャルを共変量だけで書く

観測事象 $x^\mu$ と遅延した電荷の位置 $z^\mu(\tau_*)$ の差を $R^\mu=x^\mu-z^\mu(\tau_*)$ とする。
過去光円錐上の条件は $R_\mu R^\mu=0,\ R^0>0$ である。
空間距離を $R$ とすれば $R^\mu=(R,R\vb{n})$ なので、正のスカラー

$$
D=-R_\mu U^\mu/c=\gamma R(1-\vb{n}\cdot\vb{\beta})
$$

が得られる。
加速している粒子にも適用できることを、遅延 Green 関数から確かめる。$x^0=ct$ とし、$R^\mu=x^\mu-z^\mu(\tau)$ と書く。$\Box=\nabla^2-\partial_0^2$ に対して $\Box G=-\delta^{(4)}$ を満たす遅延核は $G=\Theta(R^0)\delta(R_\mu R^\mu)/(2\pi)$ である。実際、$R=|\vb R|>0$ とすると、デルタ関数の二つの根のうち未来向きの根だけが残り、$G=\delta(R^0-R)/(4\pi R)$ になる。これは三次元の遅延核を長さ座標 $x^0$ で書いたものである。

点電荷の四元電流は $J^\mu(x')=qc\int U^\mu(\tau)\delta^{(4)}(x'-z(\tau))\,\dd\tau$ である。これを $A^\mu=\mu_0\int GJ^\mu\,\dd^4x'$ に入れる。遅延点では $\dd(R_\nu R^\nu)/\dd\tau=-2R_\nu U^\nu=2cD$ なので、残る固有時間積分の係数は $\mu_0qc/(2\pi\,2cD)$ となる。
従って四元的な形は

$$
A^\mu=\frac{\mu_0 q}{4\pi}\frac{U^\mu}{D}.
$$

実験室で $U^\mu=\gamma(c,\vb{v})$ を代入すれば $\gamma$ が約分され、時間成分と空間成分は :eqref[eq-em-lw] を再現する。
ここで $U^\mu$ と $D$ はいずれも同じ遅延点で評価する。

分母に現れる $1-\vb{n}\cdot\vb{\beta}$ は、時刻の変数変換からも Lorentz スカラーからも得られる。

### 電磁場のエネルギー運動量の発散を成分に依存せず計算する

場が荷電物質へ渡す四元力密度は $f^\nu=F^{\nu\lambda}J_\lambda$ である。成分の保存則を別々に確かめる代わりに、$\partial_\mu T^{\mu\nu}$ を展開する。一方の $F$ だけを微分すると Maxwell 方程式より $\partial_\mu F^{\mu\alpha}=-\mu_0J^\alpha$ が現れ、この項は $-f^\nu$ になる。他方を微分した残りは Bianchi の恒等式から消える。具体的には

$$
F^{\mu\alpha}(\partial_\mu F_{\lambda\alpha}
+\partial_\lambda F_{\alpha\mu}+\partial_\alpha F_{\mu\lambda})=0.
$$

第一項と第三項は、第三項で $\mu$ と $\alpha$ を交換し、二つの反対称性を使うと等しい。第二項は $-\partial_\lambda(F_{\mu\alpha}F^{\mu\alpha})/2$ である。従って $F^{\mu\alpha}\partial_\mu F_{\lambda\alpha}=\partial_\lambda(F^{\mu\alpha}F_{\mu\alpha})/4$ となり、$T$ の跡を引く項の微分に等しい。結果として $\partial_\mu T^{\mu\nu}=-f^\nu$ を得た。$\nu=0$ で $cf^0=\vb J\cdot\vb E$、空間成分では $\vb f=\rho\vb E+\vb J\times\vb B$ である。電場による仕事と磁場を含む運動量の受け渡しが一組の局所方程式にまとまる。

### 同時刻面を変えても保存電荷が一致する理由

電荷の密度はスカラーではない。それにもかかわらず総電荷が同じなのは、異なる観測者が異なる同時刻面を積分するからである。未来向き単位法線を $n^\mu$、その面上の固有体積を $\dd\Sigma$ と書けば

$$
Q(\Sigma)=-\frac1c\int_\Sigma J^\mu n_\mu\,\dd\Sigma.
$$

$ct=\text{一定}$ の超平面では $n_\mu=(-1,0,0,0)$ であり、通常の $\int\rho\dd^3r$ に戻る。局在した系で二つの同時刻面と遠方の側面に囲まれた四次元領域を考える。発散定理を $\partial_\mu J^\mu=0$ に適用すると、側面を流れる電荷が零である限り上面と下面の積分が等しい。面が交差するときは交点で二領域に分ければ同じ結論となる。無限導線の一部分のように側面で電荷が出入りする系では、部分区間だけの電荷が不変であるという結論は得られない。


:::example[向かい合う二つの光束を同じ強度に見る速度]
正の軸を進む光束のエネルギー密度を $U_+$、負の軸を進む非干渉光束の密度を $U_-$ とする。観測者を正の方向へ動かすと各密度には $(1-x)/(1+x)$ とその逆数が掛かる。新しい系の正味電力流はその差の光速倍である。

$x=v/c$、$U_+=4U_-$、$Y$ は符号付き電力流、$Y_*=cU_-=1.00\,\mathrm{W/m^2}$ とする。

$$
f(x)=- \frac{\left(x - 3\right) \left(3 x - 1\right)}{\left(x - 1\right) \left(x + 1\right)},\qquad Y=Y_* f(x).
$$

変数を少し動かしたときの読みの変化は、接線の傾きから計算できる。基準値 $x_0$ の近くで

$$
\begin{aligned}
\frac{\Delta Y}{Y_*}= f\!\left(x_0+\Delta x\right)-f(x_0)= f'(x_0)\Delta x+O((\Delta x)^2),\\
f'(x)=- \frac{2 \left(5 x^{2} - 6 x + 5\right)}{\left(x - 1\right)^{2} \left(x + 1\right)^{2}}.
\end{aligned}
$$


| $x$ | $f(x)$ | $Y$ $[\mathrm{W/m^2}]$ | $\dd Y/\dd x$ $[\mathrm{W/m^2}]$ |
|---:|---:|---:|---:|
| $0$ | $3$ | $3$ | $-10$ |
| $\frac{1}{5}$ | $1.167$ | $1.167$ | $-8.681$ |
| $\frac{1}{3}$ | $0$ | $0$ | $-9$ |
| $\frac{3}{5}$ | $-3$ | $-3$ | $-15.62$ |

零点は $2(1-x)=1+x$ より $x=1/3$ である。単一の光束を追い抜けなくても、逆向きの光束もあれば総運動量が消える慣性系を選べる。光束の差や和を変換する前に双方の進行方向を区別する必要がある。
:::

## 場の Lorentz 変換

### 観測者を替える

行列の形では $F'=\Lambda F\Lambda^\top$ である。
例えば $x$ ブーストでは $F'^{02}=\gamma(F^{02}-\beta F^{12})=\gamma(E_y/c-\beta B_z)$、

$$
F'^{12}=\gamma(F^{12}-\beta F^{02})=\gamma(B_z-\beta E_y/c)
$$

となる。

$F'^{01}=F^{01}$ は対応する二次小行列の行列式が $1$ のためである。
平行・垂直成分に分け、回転して任意の $\vb{v}$ に直せば

$$
\begin{aligned}
\vb{E}'_\parallel&=\vb{E}_\parallel,&\vb{E}'_\perp&=\gamma(\vb{E}_\perp+\vb{v}\times\vb{B}),\\
\vb{B}'_\parallel&=\vb{B}_\parallel,&\vb{B}'_\perp&=\gamma(\vb{B}_\perp-\vb{v}\times\vb{E}/c^2).
\end{aligned}
\label{eq-em-field-boost}
$$

ポテンシャルも四元ベクトルなので $\phi'=\gamma(\phi-\vb{v}\cdot\vb{A})$、

$$
A'_\parallel=\gamma(A_\parallel-v\phi/c^2)
$$

となる。
異なるゲージで比べるときは :ref[def-em-gauge] の変換も含める。

:::example[中性の電流線]
実験室で中性の長い直線に $+x$ 向きの電流 $I$ がある。
四元電流の断面積分を変換すると移動系の線電荷は

$$
\lambda'=\gamma(\lambda-vI/c^2)=-\gamma vI/c^2
$$

である。
$y>0$ の点では $B_z=\mu_0I/(2\pi y)$ なので $E'_y=-\gamma vB_z$ となり、$\lambda'/(2\pi\varepsilon_0 y)$ と一致する。

移動系で静止する試験電荷が受ける電気力を、この系の電荷密度からも求められる。
電流担体を静止から加速する過程でどの線電荷が残るかは、回路全体の電荷の供給条件を別に指定する問題である。
:cite[purcell-em]
:::

### 二つの不変量

テンソルを完全に縮約すると座標系によらない。
非零成分を順に足すと

$$
\tfrac12F_{\mu\nu}F^{\mu\nu}=B^2-E^2/c^2,\qquad \tfrac14F_{\mu\nu}{}^*F^{\mu\nu}=-\vb{E}\cdot\vb{B}/c.
\label{eq-em-invariants}
$$

第二量は空間反転で符号を変える擬スカラーで、通常のブーストや回転では不変である。
$\vb{E}\cdot\vb{B}\ne0$ なら片方だけを消す系は存在しない。
直交場で $E<cB$ なら $\vb{v}=\vb{E}\times\vb{B}/B^2$ を選べる。

$|\vb{v}|=E/B<c$、$\vb{v}\times\vb{B}=-\vb{E}$ なので $\vb{E}'=0$ である。
逆に $E>cB$ なら $\vb{v}=c^2\vb{E}\times\vb{B}/E^2$ で $\vb{B}'=0$ となる。

$E=cB$ で直交する非零場は両不変量が零の場である。
片方が零なら他方も零となり、可逆な変換で非零場を完全に消せないため、いずれの成分も消去できない。

:::example[波を追う観測者]
$+x$ 方向へ進む平面波で $E_y=cB_z$ とする。
$+x$ に速度 $v$ の観測者には $E'_y=\gamma(1-\beta)E_y$ となる。
位相 $kx-\omega t$ が不変なので $\omega'=\gamma(\omega-vk)=\gamma(1-\beta)\omega$ でもある。

$\beta=0.60$ では振幅と周波数は半分、強度は四分の一である。
:::

### 一般の場を平行にする系

直交しない $\vb{E},\vb{B}$ では一方を消せなくても両者を平行にできる。
局所的な場を $\vb{e}=\vb{E}/c,\ \vb{b}=\vb{B}$ と同じ単位にし、$\vb{e}\times\vb{b}$ が非零としてその向きの単位ベクトルを $\vb{n}$ とする。

$$
s=|\vb{e}\times\vb{b}|,\ w=e^2+b^2
$$

とする。
$\vb{v}=c\beta\vb{n}$ のブーストでは

$$
\vb{e}'=\gamma(\vb{e}+\beta\vb{n}\times\vb{b}),\qquad \vb{b}'=\gamma(\vb{b}-\beta\vb{n}\times\vb{e}).
$$

両者の外積を展開すると

$$
\vb{e}'\times\vb{b}'=\gamma^2[(1+\beta^2)s-\beta w]\vb{n}
$$

である。
小さい方の根は

$$
\beta=\frac{2s}{w+\sqrt{(b^2-e^2)^2+4(\vb{e}\cdot\vb{b})^2}}.
$$

電磁場が非零の平面波と同じ不変量をもつ場合を除けば $2s<w$ なので、この根は $1$ より小さい。
この系ではエネルギー流 $\vb{E}'\times\vb{B}'/\mu_0$ が零になる。
初めから二つの場が平行な場合はブーストを行う必要がない。

:::example[斜めに交わる電磁場]
同じ単位で $\vb{e}=(b_0,0,0)$、$\vb{b}=(b_0,b_0,0)$ とする。
$s=b_0^2,\ w=3b_0^2$ なので $\beta=(3-\sqrt5)/2\simeq0.382$ の $+z$ ブーストで平行になる。

$\vb{e}\cdot\vb{b}=b_0^2$ は非零なので、どの慣性系でも両方の場が測定される。
:::

### 一様運動する電荷の電場を静止系から求める

電荷の静止系を $S'$、電荷が $+x$ 向きに速さ $v$ で進む実験室を $S$ とする。
実験室の同一時刻における電荷から観測点までの差を $\vb{s}=(x-vt,y,z)$ とすると、静止系での座標は

$$
x'=\gamma s_x,\ y'=s_y,\ z'=s_z
$$

である。

静止系の $\vb{E}'=k_eq\vb{r}'/r'^3$ と $\vb{B}'=0$ を逆変換すれば

$$
\vb{E}=\frac{k_eq(1-\beta^2)\vb{s}}{[s_x^2+(1-\beta^2)(s_y^2+s_z^2)]^{3/2}},\qquad \vb{B}=\frac{\vb{v}\times\vb{E}}{c^2}.
$$

分母を遅延時刻を用いた等速場と比較すれば、同じ場の別表現であることがわかる。
電場線は実験室で瞬時の電荷から放射状に伸びるが、角度ごとの密度が違う。
電荷が加速すれば、加速の影響は光円錐上を伝わるのでこの静止系からの一定ブーストによる計算は使えない。

例えば $\beta=0.80$ のとき、同じ瞬時距離で横方向と縦方向の場の比は $\gamma^3=(5/3)^3=4.63$ である。
半径 $s$ の球面を通る総電束は変わらない。
実際 $u=\cos\theta$ と置くと角度積分に $(1-\beta^2)\int_{-1}^1[1-\beta^2+\beta^2u^2]^{-3/2}\dd u=2$ が現れる。

荷電粒子の速さにかかわらず Gauss の法則が成立する。

### 応力も含めて動く測定器のエネルギー密度を求める

観測者が測るエネルギー密度は $u(U)=T^{\mu\nu}U_\mu U_\nu/c^2$ である。$x$ 軸へ速度 $v$ で動く観測者について $U_\mu/c=(-\gamma,\gamma\beta,0,0)$ を代入すれば、

$$
u'=\gamma^2(T^{00}-2\beta T^{01}+\beta^2T^{11}).
$$

エネルギー密度と運動量密度だけでは変換が閉じない。積分面が傾くだけでなく、同じ時空の面を横切る運動量の流れもエネルギーの読みへ寄与する。二つの波が逆向きに伝わって元の総運動量が零でも、$T^{11}$ が正ならこの項は残る。

以下の例では $E_* =2.99792458\times10^5\,\mathrm{V/m}$、$B_*=E_*/c=1.00\,\mathrm{mT}$ とし、$\vb e=\vb E/E_*$、$\vb b=\vb B/B_*$ を用いる。場のテンソルの全成分を $B_*$ で割ったものを $F$、$T^{\mu\nu}$ を $u_*=B_*^2/\mu_0=0.796\,\mathrm{J/m^3}$ で割ったものを $T$ と略記する。実場のある一事象での瞬間値を扱い、光の周期平均に特有な $1/2$ は付け加えない。受動変換行列を $\Lambda$ とすると、上付きの二階テンソルは $F'=\Lambda F\Lambda^{\mathsf T}, T'=\Lambda T\Lambda^{\mathsf T}$ である。$\eta$ を間に挟むのは縮約の場合であり、二つの上付き添字の変換自体には挟まない。

:::example[場に平行に動く観測者]
平行板の中央に近い一様電場を考える。観測者を電場に沿って走らせると電場の平行成分は変わらない。エネルギー密度を単なる移動粒子の集まりとして変換した値と比較する。入力を順に $\vb e=\left(\begin{matrix}2 & 0 & 0\end{matrix}\right)$、$\vb b=\left(\begin{matrix}0 & 0 & 0\end{matrix}\right)$、$\boldsymbol\beta=\left(\begin{matrix}\frac{3}{5} & 0 & 0\end{matrix}\right)$ とする。$\gamma=\frac{5}{4}$ を用いる受動変換とその結果は

$$
\begin{aligned}
\Lambda&=\left(\begin{matrix}\frac{5}{4} & - \frac{3}{4} & 0 & 0\\- \frac{3}{4} & \frac{5}{4} & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right),\\
F'^{\mu\nu}&=\left(\begin{matrix}0 & 2 & 0 & 0\\-2 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{matrix}\right).
\end{aligned}
$$

第一行の空間部分から電場を読み、空間成分は $(23,31,12)$ の順で磁場を読む。その後に二つの不変量を計算すると

$$
\begin{aligned}
\vb e'&=\left(\begin{matrix}2 & 0 & 0\end{matrix}\right),\qquad \vb b'=\left(\begin{matrix}0 & 0 & 0\end{matrix}\right),\\
|\vb b'|^2-|\vb e'|^2&=-4,\qquad \vb e'\cdot\vb b'=0.
\end{aligned}
$$

同じブーストをエネルギー運動量テンソルに適用する。上付きの空間成分が Maxwell 応力の負になる規約で、

$$
\begin{aligned}
T&=\left(\begin{matrix}2 & 0 & 0 & 0\\0 & -2 & 0 & 0\\0 & 0 & 2 & 0\\0 & 0 & 0 & 2\end{matrix}\right),\\
T'&=\left(\begin{matrix}2 & 0 & 0 & 0\\0 & -2 & 0 & 0\\0 & 0 & 2 & 0\\0 & 0 & 0 & 2\end{matrix}\right).
\end{aligned}
$$

第一行を SI 単位へ戻す。$T^{0i}$ は $S_i/c$ であり、電力密度への換算には $cu_*$、運動量密度への換算には $u_*/c$ を使う。ここでは後者の三成分とエネルギー密度を示す。

| 観測量 | 変換前 | 変換後 |
|---|---:|---:|
| エネルギー密度 $[\mathrm{J/m^3}]$ | $1.592$ | $1.592$ |
| $x$ 運動量密度 $[\mathrm{nN\,s/m^3}]$ | $0$ | $0$ |
| $y$ 運動量密度 $[\mathrm{nN\,s/m^3}]$ | $0$ | $0$ |
| $z$ 運動量密度 $[\mathrm{nN\,s/m^3}]$ | $0$ | $0$ |

二つの行列の Minkowski 跡はともに零である。変換後の三次元場を直接 $u$ と $\vb E\times\vb B/\mu_0$ に代入しても同じ第一行が得られる。張力に対応する負の運動量流束が、観測者の運動による二次の項をちょうど補償する。これは場の向きに沿うブーストに固有の性質である。
:::

:::example[磁場を横切る観測者の電力流]
静止系で導体の外に純磁場しかない事象を選ぶ。横に動く測定器は電場も測るため、新しい系での Poynting ベクトルの向きを求める。入力を順に $\vb e=\left(\begin{matrix}0 & 0 & 0\end{matrix}\right)$、$\vb b=\left(\begin{matrix}0 & 0 & 2\end{matrix}\right)$、$\boldsymbol\beta=\left(\begin{matrix}\frac{3}{5} & 0 & 0\end{matrix}\right)$ とする。$\gamma=\frac{5}{4}$ を用いる受動変換とその結果は

$$
\begin{aligned}
\Lambda&=\left(\begin{matrix}\frac{5}{4} & - \frac{3}{4} & 0 & 0\\- \frac{3}{4} & \frac{5}{4} & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right),\\
F'^{\mu\nu}&=\left(\begin{matrix}0 & 0 & - \frac{3}{2} & 0\\0 & 0 & \frac{5}{2} & 0\\\frac{3}{2} & - \frac{5}{2} & 0 & 0\\0 & 0 & 0 & 0\end{matrix}\right).
\end{aligned}
$$

第一行の空間部分から電場を読み、空間成分は $(23,31,12)$ の順で磁場を読む。その後に二つの不変量を計算すると

$$
\begin{aligned}
\vb e'&=\left(\begin{matrix}0 & - \frac{3}{2} & 0\end{matrix}\right),\qquad \vb b'=\left(\begin{matrix}0 & 0 & \frac{5}{2}\end{matrix}\right),\\
|\vb b'|^2-|\vb e'|^2&=4,\qquad \vb e'\cdot\vb b'=0.
\end{aligned}
$$

同じブーストをエネルギー運動量テンソルに適用する。上付きの空間成分が Maxwell 応力の負になる規約で、

$$
\begin{aligned}
T&=\left(\begin{matrix}2 & 0 & 0 & 0\\0 & 2 & 0 & 0\\0 & 0 & 2 & 0\\0 & 0 & 0 & -2\end{matrix}\right),\\
T'&=\left(\begin{matrix}\frac{17}{4} & - \frac{15}{4} & 0 & 0\\- \frac{15}{4} & \frac{17}{4} & 0 & 0\\0 & 0 & 2 & 0\\0 & 0 & 0 & -2\end{matrix}\right).
\end{aligned}
$$

第一行から得るエネルギー密度と運動量密度は次の通りである。

| 観測量 | 変換前 | 変換後 |
|---|---:|---:|
| エネルギー密度 $[\mathrm{J/m^3}]$ | $1.592$ | $3.382$ |
| $x$ 運動量密度 $[\mathrm{nN\,s/m^3}]$ | $0$ | $-9.954$ |
| $y$ 運動量密度 $[\mathrm{nN\,s/m^3}]$ | $0$ | $0$ |
| $z$ 運動量密度 $[\mathrm{nN\,s/m^3}]$ | $0$ | $0$ |

元の系にエネルギー輸送が見えなくても、観測者の交代で場の運動量密度が現れる。速度の符号は受動変換の定義から決まる。
:::

:::example[光を追うときの場と強度]
正の軸方向へ進む平面波を後ろから追う観測者を考える。電場・磁場はいずれも小さくなるが、両方の不変量が零のため、有限速度で場を完全に消すことはできない。入力を順に $\vb e=\left(\begin{matrix}0 & 1 & 0\end{matrix}\right)$、$\vb b=\left(\begin{matrix}0 & 0 & 1\end{matrix}\right)$、$\boldsymbol\beta=\left(\begin{matrix}\frac{4}{5} & 0 & 0\end{matrix}\right)$ とする。$\gamma=\frac{5}{3}$ を用いる受動変換とその結果は

$$
\begin{aligned}
\Lambda&=\left(\begin{matrix}\frac{5}{3} & - \frac{4}{3} & 0 & 0\\- \frac{4}{3} & \frac{5}{3} & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right),\\
F'^{\mu\nu}&=\left(\begin{matrix}0 & 0 & \frac{1}{3} & 0\\0 & 0 & \frac{1}{3} & 0\\- \frac{1}{3} & - \frac{1}{3} & 0 & 0\\0 & 0 & 0 & 0\end{matrix}\right).
\end{aligned}
$$

第一行の空間部分から電場を読み、空間成分は $(23,31,12)$ の順で磁場を読む。その後に二つの不変量を計算すると

$$
\begin{aligned}
\vb e'&=\left(\begin{matrix}0 & \frac{1}{3} & 0\end{matrix}\right),\qquad \vb b'=\left(\begin{matrix}0 & 0 & \frac{1}{3}\end{matrix}\right),\\
|\vb b'|^2-|\vb e'|^2&=0,\qquad \vb e'\cdot\vb b'=0.
\end{aligned}
$$

同じブーストをエネルギー運動量テンソルに適用する。上付きの空間成分が Maxwell 応力の負になる規約で、

$$
\begin{aligned}
T&=\left(\begin{matrix}1 & 1 & 0 & 0\\1 & 1 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{matrix}\right),\\
T'&=\left(\begin{matrix}\frac{1}{9} & \frac{1}{9} & 0 & 0\\\frac{1}{9} & \frac{1}{9} & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{matrix}\right).
\end{aligned}
$$

第一行から得るエネルギー密度と運動量密度は次の通りである。

| 観測量 | 変換前 | 変換後 |
|---|---:|---:|
| エネルギー密度 $[\mathrm{J/m^3}]$ | $0.7958$ | $0.08842$ |
| $x$ 運動量密度 $[\mathrm{nN\,s/m^3}]$ | $2.654$ | $0.2949$ |
| $y$ 運動量密度 $[\mathrm{nN\,s/m^3}]$ | $0$ | $0$ |
| $z$ 運動量密度 $[\mathrm{nN\,s/m^3}]$ | $0$ | $0$ |

エネルギー密度の倍率は電場振幅の倍率の二乗である。単一光子のエネルギーの Doppler 因子だけを掛けても密度にはならず、到来頻度の変化も含まれる。
:::

:::example[電場を消すドリフト系]
直交する二場で磁気成分が優勢な配置を選ぶ。電場による力と観測者の運動に伴う磁気力が相殺する速度を選ぶと、荷電粒子をその系で静止させられる。入力を順に $\vb e=\left(\begin{matrix}0 & 3 & 0\end{matrix}\right)$、$\vb b=\left(\begin{matrix}0 & 0 & 5\end{matrix}\right)$、$\boldsymbol\beta=\left(\begin{matrix}\frac{3}{5} & 0 & 0\end{matrix}\right)$ とする。$\gamma=\frac{5}{4}$ を用いる受動変換とその結果は

$$
\begin{aligned}
\Lambda&=\left(\begin{matrix}\frac{5}{4} & - \frac{3}{4} & 0 & 0\\- \frac{3}{4} & \frac{5}{4} & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right),\\
F'^{\mu\nu}&=\left(\begin{matrix}0 & 0 & 0 & 0\\0 & 0 & 4 & 0\\0 & -4 & 0 & 0\\0 & 0 & 0 & 0\end{matrix}\right).
\end{aligned}
$$

第一行の空間部分から電場を読み、空間成分は $(23,31,12)$ の順で磁場を読む。その後に二つの不変量を計算すると

$$
\begin{aligned}
\vb e'&=\left(\begin{matrix}0 & 0 & 0\end{matrix}\right),\qquad \vb b'=\left(\begin{matrix}0 & 0 & 4\end{matrix}\right),\\
|\vb b'|^2-|\vb e'|^2&=16,\qquad \vb e'\cdot\vb b'=0.
\end{aligned}
$$

同じブーストをエネルギー運動量テンソルに適用する。上付きの空間成分が Maxwell 応力の負になる規約で、

$$
\begin{aligned}
T&=\left(\begin{matrix}17 & 15 & 0 & 0\\15 & 17 & 0 & 0\\0 & 0 & 8 & 0\\0 & 0 & 0 & -8\end{matrix}\right),\\
T'&=\left(\begin{matrix}8 & 0 & 0 & 0\\0 & 8 & 0 & 0\\0 & 0 & 8 & 0\\0 & 0 & 0 & -8\end{matrix}\right).
\end{aligned}
$$

第一行から得るエネルギー密度と運動量密度は次の通りである。

| 観測量 | 変換前 | 変換後 |
|---|---:|---:|
| エネルギー密度 $[\mathrm{J/m^3}]$ | $13.53$ | $6.366$ |
| $x$ 運動量密度 $[\mathrm{nN\,s/m^3}]$ | $39.82$ | $0$ |
| $y$ 運動量密度 $[\mathrm{nN\,s/m^3}]$ | $0$ | $0$ |
| $z$ 運動量密度 $[\mathrm{nN\,s/m^3}]$ | $0$ | $0$ |

新系で磁場の大きさが元の値より小さくなることは第一不変量からもわかる。電場が消える速度の横成分を求めた後、残る磁場に平行な追加ブーストも許される。
:::

:::example[電場と磁場が直交しない場合]
測定器を基準系の第三軸に沿って動かす。ブースト方向の成分がある場合、すべての横成分だけを同じ符号で増幅する手続きではテンソルの変換を再現できない。入力を順に $\vb e=\left(\begin{matrix}1 & 1 & 0\end{matrix}\right)$、$\vb b=\left(\begin{matrix}0 & 1 & 1\end{matrix}\right)$、$\boldsymbol\beta=\left(\begin{matrix}0 & 0 & \frac{3}{5}\end{matrix}\right)$ とする。$\gamma=\frac{5}{4}$ を用いる受動変換とその結果は

$$
\begin{aligned}
\Lambda&=\left(\begin{matrix}\frac{5}{4} & 0 & 0 & - \frac{3}{4}\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\- \frac{3}{4} & 0 & 0 & \frac{5}{4}\end{matrix}\right),\\
F'^{\mu\nu}&=\left(\begin{matrix}0 & \frac{1}{2} & \frac{5}{4} & 0\\- \frac{1}{2} & 0 & 1 & - \frac{1}{2}\\- \frac{5}{4} & -1 & 0 & \frac{3}{4}\\0 & \frac{1}{2} & - \frac{3}{4} & 0\end{matrix}\right).
\end{aligned}
$$

第一行の空間部分から電場を読み、空間成分は $(23,31,12)$ の順で磁場を読む。その後に二つの不変量を計算すると

$$
\begin{aligned}
\vb e'&=\left(\begin{matrix}\frac{1}{2} & \frac{5}{4} & 0\end{matrix}\right),\qquad \vb b'=\left(\begin{matrix}\frac{3}{4} & \frac{1}{2} & 1\end{matrix}\right),\\
|\vb b'|^2-|\vb e'|^2&=0,\qquad \vb e'\cdot\vb b'=1.
\end{aligned}
$$

同じブーストをエネルギー運動量テンソルに適用する。上付きの空間成分が Maxwell 応力の負になる規約で、

$$
\begin{aligned}
T&=\left(\begin{matrix}2 & 1 & -1 & 1\\1 & 1 & -1 & 0\\-1 & -1 & 0 & -1\\1 & 0 & -1 & 1\end{matrix}\right),\\
T'&=\left(\begin{matrix}\frac{29}{16} & \frac{5}{4} & - \frac{1}{2} & - \frac{11}{16}\\\frac{5}{4} & 1 & -1 & - \frac{3}{4}\\- \frac{1}{2} & -1 & 0 & - \frac{1}{2}\\- \frac{11}{16} & - \frac{3}{4} & - \frac{1}{2} & \frac{13}{16}\end{matrix}\right).
\end{aligned}
$$

第一行から得るエネルギー密度と運動量密度は次の通りである。

| 観測量 | 変換前 | 変換後 |
|---|---:|---:|
| エネルギー密度 $[\mathrm{J/m^3}]$ | $1.592$ | $1.442$ |
| $x$ 運動量密度 $[\mathrm{nN\,s/m^3}]$ | $2.654$ | $3.318$ |
| $y$ 運動量密度 $[\mathrm{nN\,s/m^3}]$ | $-2.654$ | $-1.327$ |
| $z$ 運動量密度 $[\mathrm{nN\,s/m^3}]$ | $2.654$ | $-1.825$ |

第二不変量が非零なので、どの慣性系でも電場と磁場の双方が残る。両場が存在することと、ある系で電力の流れが存在することは別の条件である。
:::

:::example[伝播方向と斜交する観測者]
斜めに走る観測者に対し、光が来る方向と偏光方向を同時に求める。運動を各軸の独立な一次元ブーストとみなすと、有限速度では回転まで混入する。入力を順に $\vb e=\left(\begin{matrix}0 & 1 & 0\end{matrix}\right)$、$\vb b=\left(\begin{matrix}0 & 0 & 1\end{matrix}\right)$、$\boldsymbol\beta=\left(\begin{matrix}\frac{9}{25} & \frac{12}{25} & 0\end{matrix}\right)$ とする。$\gamma=\frac{5}{4}$ を用いる受動変換とその結果は

$$
\begin{aligned}
\Lambda&=\left(\begin{matrix}\frac{5}{4} & - \frac{9}{20} & - \frac{3}{5} & 0\\- \frac{9}{20} & \frac{109}{100} & \frac{3}{25} & 0\\- \frac{3}{5} & \frac{3}{25} & \frac{29}{25} & 0\\0 & 0 & 0 & 1\end{matrix}\right),\\
F'^{\mu\nu}&=\left(\begin{matrix}0 & \frac{12}{25} & \frac{16}{25} & 0\\- \frac{12}{25} & 0 & \frac{4}{5} & 0\\- \frac{16}{25} & - \frac{4}{5} & 0 & 0\\0 & 0 & 0 & 0\end{matrix}\right).
\end{aligned}
$$

第一行の空間部分から電場を読み、空間成分は $(23,31,12)$ の順で磁場を読む。その後に二つの不変量を計算すると

$$
\begin{aligned}
\vb e'&=\left(\begin{matrix}\frac{12}{25} & \frac{16}{25} & 0\end{matrix}\right),\qquad \vb b'=\left(\begin{matrix}0 & 0 & \frac{4}{5}\end{matrix}\right),\\
|\vb b'|^2-|\vb e'|^2&=0,\qquad \vb e'\cdot\vb b'=0.
\end{aligned}
$$

同じブーストをエネルギー運動量テンソルに適用する。上付きの空間成分が Maxwell 応力の負になる規約で、

$$
\begin{aligned}
T&=\left(\begin{matrix}1 & 1 & 0 & 0\\1 & 1 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{matrix}\right),\\
T'&=\left(\begin{matrix}\frac{16}{25} & \frac{64}{125} & - \frac{48}{125} & 0\\\frac{64}{125} & \frac{256}{625} & - \frac{192}{625} & 0\\- \frac{48}{125} & - \frac{192}{625} & \frac{144}{625} & 0\\0 & 0 & 0 & 0\end{matrix}\right).
\end{aligned}
$$

第一行から得るエネルギー密度と運動量密度は次の通りである。

| 観測量 | 変換前 | 変換後 |
|---|---:|---:|
| エネルギー密度 $[\mathrm{J/m^3}]$ | $0.7958$ | $0.5093$ |
| $x$ 運動量密度 $[\mathrm{nN\,s/m^3}]$ | $2.654$ | $1.359$ |
| $y$ 運動量密度 $[\mathrm{nN\,s/m^3}]$ | $0$ | $-1.019$ |
| $z$ 運動量密度 $[\mathrm{nN\,s/m^3}]$ | $0$ | $0$ |

変換後も電力流の速さは光速であるが向きは変わる。到来方向の変化と振幅の変化は同じ場の変換から同時に得られる。
:::

## 荷電粒子の相対論的運動

### 四元力から仕事と運動方程式を得る

静止質量 $m$、四元運動量 $p^\mu=mU^\mu=(\mathcal E/c,\vb{p})$ について $\mathcal E=\gamma_umc^2$、$\vb{p}=\gamma_um\vb{u}$、$\mathcal E^2=m^2c^4+p^2c^2$ である。

荷電粒子の運動方程式を

$$
\frac{\dd p^\mu}{\dd\tau}=qF^{\mu\nu}U_\nu
\label{eq-em-four-force}
$$

と書く。
空間成分は $U_0=-\gamma_uc$、$U_i=\gamma_u u_i$ より $\gamma_u\dd\vb{p}/\dd t=\gamma_u q(\vb{E}+\vb{u}\times\vb{B})$、時間成分は

$$
\dd\mathcal E/\dd t=q\vb{E}\cdot\vb{u}
$$

である。

反対称性から $U_\mu F^{\mu\nu}U_\nu=0$ であり、$p_\mu p^\mu=-m^2c^2$ が保存する。
速度に平行な力は $F_\parallel=\gamma_u^3ma_\parallel$、垂直なら $F_\perp=\gamma_u ma_\perp$ である。

これは $\dd\gamma_u/\dd t=\gamma_u^3\vb{u}\cdot\vb{a}/c^2$ を $\dd(\gamma_um\vb{u})/\dd t$ に代入して得られる。

### 一様場中の軌道

静止から $x$ 方向の一様場で加速する場合を考える。
$qE>0$、$a_0=qE/m$ とすれば $p_x=qEt$ である。
エネルギーと速度の関係 $u=c^2p_x/\mathcal E$ を使って積分すると

$$
u(t)=\frac{a_0t}{\sqrt{1+(a_0t/c)^2}},\qquad x(t)=\frac{c^2}{a_0}\left[\sqrt{1+(a_0t/c)^2}-1\right].
\label{eq-em-hyperbola}
$$

従って

$$
(x+c^2/a_0)^2-(ct)^2=(c^2/a_0)^2
$$

であり、$t=(c/a_0)\sinh(a_0\tau/c)$、

$$
x=(c^2/a_0)[\cosh(a_0\tau/c)-1]
$$

と書ける。
固有加速度は一定でも実験室での加速度は $a_0/\gamma_u^3$ に減少する。

一様磁場 $\vb{B}=B\hat{\vb{z}}$ のみなら仕事率が零で $\gamma_u$ は一定である。
運動方程式は $\dot p_x=qBp_y/(\gamma_um)$、$\dot p_y=-qBp_x/(\gamma_um)$、$\dot p_z=0$ である。

$p_x+ip_y$ は $e^{-iqBt/(\gamma_um)}$ に比例するため軌道は螺旋となり、半径は $p_\perp/(|q|B)$、実験室の角周波数の大きさは $|q|B/(\gamma_um)$ である。
放射による損失を無視できる時間に限る。

:::example[電子の軌道半径]
運動エネルギー $K=1.602\times10^{-13}\,\mathrm J$ の電子が $0.100\,\mathrm T$ の場に垂直に入る。

$$
m_ec^2=8.187\times10^{-14}\,\mathrm J
$$

より $pc=\sqrt{K(K+2m_ec^2)}=2.278\times10^{-13}\,\mathrm J$、$r=p/(eB)=4.74\,\mathrm{cm}$ となる。

非相対論式 $p=\sqrt{2m_eK}$ では半径を小さく見積もる。
:::

:::example[直交場中の案内中心]
$E<cB$ の一様直交場では $\vb{v}_d=\vb{E}\times\vb{B}/B^2$ の系で純磁場になる。
その系で磁場に垂直な円を一周する二事象は同じ空間位置で起きる。
実験室に戻すと二事象の変位は $\vb{v}_d$ と時間差の積である。

従って実験室での一周あたりの平均速度は $\vb{v}_d$ となる。
電荷の符号によって回転方向は逆になるが案内中心の漂移方向は同じである。
:::

### 相対論的な正準運動量と保存エネルギー

粒子の作用を

$$
\mathcal S_p=-mc^2\int\dd\tau+q\int A_\mu\dd x^\mu
$$

と置く。
世界線に沿って $A_\mu\dd x^\mu=-\phi\,\dd t+\vb{A}\cdot\dd\vb{r}$ なので、実験室時刻を変数とすれば

$$
L=-mc^2\sqrt{1-u^2/c^2}+q\vb{u}\cdot\vb{A}-q\phi.
$$

速度に関する微分は

$$
\vb{P}=\gamma_um\vb{u}+q\vb{A}
$$

となる。
時間微分から $q(\partial_t+\vb{u}\cdot\nabla)\vb{A}$ を取り出し、位置に関する偏微分との差をとると、空間成分の Lorentz 力が得られる。

$\mathcal H=\vb{P}\cdot\vb{u}-L$ を $\vb{P}$ の関数に直せば

$$
\mathcal H(\vb{r},\vb{P},t)=\sqrt{m^2c^4+c^2|\vb{P}-q\vb{A}|^2}+q\phi.
$$

時間に陽に依存しないポテンシャルではこの $\mathcal H$ が保存する。
静電加速で $\Delta K=-q\Delta\phi$ となるのは、非相対論領域だけでなく、このエネルギー保存が成立する全速度域で正しい。
正準運動量と運動学的な $\gamma_um\vb{u}$ は、ゲージ変換したときの振る舞いが異なる。

:::example[磁場を横切る軌道の保存量]
一様な $B_z=B$ に対して $A_y=Bx$、他成分と $\phi$ は零とする。
$L$ は $y$ に陽に依存しないので $P_y=p_y+qBx$ が保存する。
従って $x_c=x+p_y/(qB)$ は一定である。
同じ運動方程式から $y_c=y-p_x/(qB)$ も一定で、

$$
(x-x_c)^2+(y-y_c)^2=p_\perp^2/(q^2B^2)
$$

となる。

ゲージを変えて正準運動量の表示が変わっても、軌道中心と半径は同じである。
:::

### 放射電力と固有加速度

四元加速度 $a^\mu=\dd U^\mu/\dd\tau$ の時間成分は $\gamma_u^4(\vb{u}\cdot\dot{\vb{u}})/c$、空間成分は $\gamma_u^2\dot{\vb{u}}+\gamma_u^4(\vb{u}\cdot\dot{\vb{u}})\vb{u}/c^2$ である。

運動方向を $x$ とする瞬間には、

$$
a^0=\gamma_u^4\beta a_\parallel,\ a^1=\gamma_u^4a_\parallel,\ \vb{a}_\perp^{(4)}=\gamma_u^2\vb{a}_\perp
$$

となる。

従って

$$
a_\mu a^\mu=\gamma_u^6a_\parallel^2+\gamma_u^4a_\perp^2.
$$

これは瞬間静止系の三次元加速度の二乗と等しい。
[電磁放射](/physics/electromagnetism/radiation/) で得た放出電力は

$$
P=q^2a_\mu a^\mu/(6\pi\varepsilon_0c^3)
$$

と書ける。

瞬間静止系では微小時間の放射の全運動量が零なので、外向き放射の四元運動量の増加率は

$$
\dd p^\mu_{\rm rad}/\dd\tau=(P/c^2)U^\mu
$$

である。
これは遠方へ運び去られる運動量の収支であり、粒子自身への四元力とは区別する必要がある。

そのまま負にした量は $U^\mu$ に平行で、質量一定の粒子に働く四元力が満たす直交性と両立しないためである。

### 三次元の力は四元ベクトルの空間成分か

同じ粒子に働く力を二つの観測者が測るとき、四元力と三次元の力を混同すると余分な $\gamma$ が現れる。
粒子速度を $\vb{u}$、三次元の力を $\vb{f}=\dd\vb{p}/\dd t$ とする。
四元力は

$$
K^\mu=\gamma_u(\vb{f}\cdot\vb{u}/c,\vb{f})
$$

である。
$x$ 方向に速度 $v$ の $S'$ に対し、時間の変換より

$$
\gamma_{u'}=\gamma_v\gamma_u(1-vu_x/c^2)
$$

を得る。
$K'^\mu=\Lambda^\mu{}_\nu K^\nu$ の空間成分を $\gamma_{u'}$ で割ると

$$
f'_x=\frac{f_x-v(\vb{f}\cdot\vb{u})/c^2}{1-vu_x/c^2},\qquad \vb{f}'_\perp=\frac{\vb{f}_\perp}{\gamma_v(1-vu_x/c^2)}.
$$

例えば粒子が $S$ で瞬間的に静止しているなら平行力は変わらず、垂直力は $1/\gamma_v$ 倍になる。
この場合には $S'$ で電場の垂直成分が $\gamma_v$ 倍になる一方、同じ粒子に新しい磁気力も働く。
その両方を合計すると上の力になる。

:::example[純電場による横向きの力を二つの方法で求める]
$S$ で $\vb{E}=E\hat{\vb{y}},\ \vb{B}=0$、粒子は静止していて $\vb{f}=qE\hat{\vb{y}}$ とする。
$S'$ では粒子は $-v\hat{\vb{x}}$ で進み、$E'_y=\gamma_v E,\ B'_z=-\gamma_v vE/c^2$ なので

$$
f'_y=q(E'_y+vB'_z)=q\gamma_vE(1-v^2/c^2)=qE/\gamma_v.
$$

$\beta=0.80$ なら力は $0.60$ 倍である。
粒子の静止系以外で電場成分だけから力を計算すると、速度に依存する第二項が抜けてしまう。
:::

### 平面波パルス中の粒子の軌道を厳密に求める

平面波の磁場は粒子に仕事をしないが、横向きに加速された粒子を波の進行方向へ曲げる。そのため、電場だけを残した運動方程式では粒子の縦方向の変位を求められない。真空中を $+z$ へ進む平面波を、$\eta=t-z/c$、$\phi=0$、$\vb{A}=A(\eta)\hat{\vb{x}}$ と書く。電場は $E_x=-A'(\eta)$、磁場は $B_y=E_x/c$ である。波が到達する前は $A=0$、粒子の質量は $m$、電荷は $q$、初期速度は零とし、放射による運動の変化を無視する。

$x$ 方向の運動方程式は $\dd p_x/\dd t=q(E_x-v_z B_y)=q E_x(1-v_z/c)$ である。一方、軌道上で $\dd\eta/\dd t=1-v_z/c$ なので、直ちに $p_x+qA=0$ を得る。$y$ 方向の運動量も零のままである。残る二つの運動方程式を別々に積分する必要はない。エネルギー変化 $\dd(\gamma mc^2)/\dd t=qv_x E_x$ と縦運動量変化 $\dd p_z/\dd t=qv_x B_y$ が同じ場を含むため、

$$
\frac{\dd}{\dd t}\left(\gamma-\frac{p_z}{mc}\right)=0,\qquad
\gamma-\frac{p_z}{mc}=1
$$

が従う。四元的には、波の四元波数 $k^\mu=(\omega/c,0,0,\omega/c)$ に対する $k_\mu p^\mu$ の保存である。ここでの $\omega$ は位相変数を無次元にするための任意の基準角周波数でよく、パルスを単色とは仮定していない。

$u_i=p_i/(mc)$、$\xi(\eta)=-qA(\eta)/(mc)$ と置く。質量殻条件 $\gamma^2=1+u_x^2+u_z^2$ と $\gamma-u_z=1$ を連立すれば、解は

$$
u_x=\xi,\qquad u_z=\xi^2/2,\qquad \gamma=1+\xi^2/2
$$

である。さらに $\dd\eta/\dd\tau=\gamma-u_z=1$ なので、粒子に沿った固有時間と $\eta$ の増分は同じになる。従って軌道そのものは $\dd x/\dd\eta=c\xi$、$\dd z/\dd\eta=c\xi^2/2$、$\dd t/\dd\eta=1+\xi^2/2$ を一度積分すれば得られる。この結果は弱い場で横運動が場の一次、縦運動が二次から始まることも示している。

:::example[波を通過した後に残るもの]
$\xi=0.200\sin^2(\pi\eta/T)$ を $0<\eta<T$ に限った理想パルスとし、前後は $\xi=0$ とする。中央で電子の運動エネルギーは $(\gamma-1)m_ec^2=0.0200\,m_ec^2=10.2\,\mathrm{keV}$ に達する。通過後は $\xi=0$ に戻るため、横・縦の両運動量も零に戻る。しかしその間の縦変位は正である。$\int_0^T\sin^4(\pi\eta/T)\,\dd\eta=3T/8$ なので、$T=1.00\,\mathrm{ps}$ なら $\Delta z=3c(0.200)^2 T/16=2.25\,\mu\mathrm m$ である。従って運動の前後のエネルギーだけを測ると、この変位は分からない。集束して横方向にも変化する実際の波や、入射前後でベクトルポテンシャルが異なる場合には、同じ終状態を推定できない。
:::

### 平行で一定な電場と磁場中の固有時間解

$\vb E=E\hat{\vb z}$、$\vb B=B\hat{\vb z}$ とすれば、混合テンソル $F^\mu{}_\nu$ は $(0,3)$ の二成分を結ぶ対称ブロックと、$(1,2)$ の二成分を回す反対称ブロックに分かれる。固有時間で $\dd U^\mu/\dd\tau=(q/m)F^\mu{}_\nu U^\nu$ を書き、$\alpha=qE/(mc),\ \omega=qB/m$ と置くと

$$
\begin{aligned}
\dot U^0&=\alpha U^3,&\dot U^3&=\alpha U^0,\\
\dot U^1&=\omega U^2,&\dot U^2&=-\omega U^1 .
\end{aligned}
$$

上段の和と差はそれぞれ $\pm\alpha$ 倍の微分方程式、下段の $U^1+iU^2$ は $-i\omega$ 倍の微分方程式である。従って初期値に添字 $i$ を付ければ

$$
U^0=U_i^0\cosh\alpha\tau+U_i^3\sinh\alpha\tau,\quad
U^3=U_i^3\cosh\alpha\tau+U_i^0\sinh\alpha\tau,\quad
U^1+iU^2=(U_i^1+iU_i^2)e^{-i\omega\tau}.
$$

横運動の振幅は一定でも、縦方向の電場が正の仕事をすると $\gamma=U^0/c$ が増える。位相 $\arg(U^1+iU^2)$ を実験室時間で微分すると、$\dd t/\dd\tau=\gamma$ なので、符号を含む角速度は $-\omega/\gamma$、その大きさは $|\omega|/\gamma$ となる。四元速度の規格化は、縦ブロックで $-(U^0)^2+(U^3)^2$、横ブロックで $(U^1)^2+(U^2)^2$ がそれぞれ一定であることから確認できる。軌道を得るにはこれらを $\tau$ で積分すればよい。例えば $U_i^3=0$、$U_i^1=u_\perp$、$U_i^2=0$ なら、初期点を原点として

$$
ct=U_i^0\frac{\sinh\alpha\tau}{\alpha},\quad
z=U_i^0\frac{\cosh\alpha\tau-1}{\alpha},\quad
x=u_\perp\frac{\sin\omega\tau}{\omega},\quad
y=u_\perp\frac{\cos\omega\tau-1}{\omega}.
$$

ここで $u_\perp$ は三速度ではなく四元速度の成分で、$U_i^0=\sqrt{c^2+u_\perp^2}$ である。電荷の符号を変えると横回転と縦加速はともに反転するが、同じ初期条件に対する $t(\tau)$ は変わらない。


:::example[一定電場による運動を非相対論近似と比較する]
静止粒子に一様な電場を加えると、実験室時間では $p=qEt$ である。運動エネルギーは全エネルギーと静止エネルギーの差 $K=mc^2[\sqrt{1+(p/(mc))^2}-1]$ となる。

$x=|qE|t/(mc)$、$Y=K$、電子なら $Y_*=mc^2=511.0\,\mathrm{keV}$ とする。場の放射損失は無視できると仮定する。

$$
f(x)=\sqrt{x^{2} + 1} - 1,\qquad Y=Y_* f(x).
$$

$x=0$ のまわりで先頭の二項を残した近似は

$$
f_{\rm a}(x)=- \frac{x^{4}}{8} + \frac{x^{2}}{2}.
$$

次の表の差は $Y_*[f(x)-f_{\rm a}(x)]$ である。

| $x$ | $f(x)$ | $Y$ $[\mathrm{keV}]$ | 厳密値と近似値の差 $[\mathrm{keV}]$ |
|---:|---:|---:|---:|
| $\frac{1}{20}$ | $0.001249$ | $0.6384$ | $4.982\times10^{-7}$ |
| $\frac{1}{4}$ | $0.03078$ | $15.73$ | $0.007505$ |
| $\frac{1}{2}$ | $0.118$ | $60.32$ | $0.4326$ |
| $1$ | $0.4142$ | $211.7$ | $20.04$ |

小さい運動量では $p^2/(2m)$ が最初の項である。高次項は一定運動量での相対論的な運動エネルギーを低下させる。一定速度で古典式と比べたときの補正の符号とは区別する必要がある。
:::

:::example[平面波パルス中の縦運動量とエネルギー]
初め静止していた電荷が一方向の平面波に入ると、$p^0-p^z=mc$ が保存される。横正準運動量も保存され、初期のポテンシャルを零に選ぶと横運動量の大きさは $|qA_\perp|$ となる。質量殻条件 $-(p^0)^2+p_z^2+p_\perp^2=-m^2c^2$ から $p_z/(mc)=x^2/2$ である。

$x=|q A_\perp|/(mc)$、$Y$ をその瞬間の運動エネルギー、$Y_*=mc^2=511.0\,\mathrm{keV}$（電子）とする。ここで $p^0=mc+p^z$ を用いる。

$$
f(x)=\frac{x^{2}}{2},\qquad Y=Y_* f(x).
$$

二点の読みの差を求めるには、基準を同じに保ち、振幅の換算を最後に一度だけ行う。最初と最後の設定で

$$
\Delta Y/Y_* = f\!\left(2\right)-f\!\left(0\right)=2.
$$

第四列は初期設定からの変化である。

| $x$ | $f(x)$ | $Y$ $[\mathrm{keV}]$ | $Y-Y(x_0)$ $[\mathrm{keV}]$ |
|---:|---:|---:|---:|
| $0$ | $0$ | $0$ | $0$ |
| $\frac{1}{4}$ | $0.03125$ | $15.97$ | $15.97$ |
| $\frac{1}{2}$ | $0.125$ | $63.88$ | $63.88$ |
| $1$ | $0.5$ | $255.5$ | $255.5$ |
| $2$ | $2$ | $1022$ | $1022$ |

波の中では粒子に縦運動量が生まれるが、パルスの前後で横ポテンシャルが同じなら、自由な平面波モデルの最終運動量も初期値に戻る。有限焦点や他の外力を含む系にこの最終値を外挿することはできない。
:::

### 零不変量の一定場を固有時間で積分する

$\vb E=E\hat{\vb x}$、$\vb B=(E/c)\hat{\vb y}$ の一定交差場を考える。両不変量が零で、有限速度で純電場にも純磁場にも移れない場合である。それでも運動方程式は厳密に解ける。$a=qE/(mc)$ とし、固有時間を独立変数にとると、混合テンソルによる運動方程式は

$$
\dot U^0=aU^x,\quad \dot U^x=a(U^0-U^z),\quad \dot U^y=0,\quad \dot U^z=aU^x .
$$

ここでは電場が $+x$、Poynting ベクトルが $+z$ を向くことが符号の照合になる。第一式と第四式の差が零だから $s=U^0-U^z$ は一定である。初めに粒子が静止していれば $s=c$ であり、$U^x=ac\tau$、$U^z=a^2c\tau^2/2$、$U^0=c+U^z$、$U^y=0$ が順に得られる。直接二乗すれば $-(U^0)^2+(U^x)^2+(U^z)^2=-c^2$ である。さらに一度積分すれば

$$
t=\tau+a^2\tau^3/6,\quad x=ac\tau^2/2,\quad z=a^2c\tau^3/6 .
$$

$x$ 方向の運動は電荷の符号を変えると反転するが、$z$ 方向は変わらない。$x$ 速度と電荷の両方の符号が変わるため $q\vb v\times\vb B$ の縦成分が保たれるからである。大きな $|a\tau|$ で三速度は縦に光速へ近づくが、四元速度の規格化からどの有限の固有時間でも光速には達しない。有限長の装置に適用する場合は、上の軌道が一様場の領域に収まる時間までに限る。

```tikz caption="交差場中で正の電荷が静止から進む軌道。初めは電場の向きへ、その後は電力流の向きにも進む。"
\begin{tikzpicture}
\draw[->] (0,0)--(3.4,0) node[right] {$z$};
\draw[->] (0,0)--(0,2.5) node[above] {$x$};
\draw[->,thick] (0,0) .. controls (.04,.65) and (.7,1.5) .. (2.8,2.2);
\draw[->] (2.7,.3)--(2.7,1.1) node[right] {$\mathbf E$};
\node at (1.1,.2) {$\mathbf S$};
\draw[->] (.6,.4)--(1.5,.4);
\end{tikzpicture}
```

## 演習問題

::::exercise{level=basic}
$v=0.60c$ で移動する粒子の実験室時間 $10.0\,\mathrm{ns}$ に対する固有時間を求めよ。

:::solution
$\gamma=1/\sqrt{1-0.60^2}=1.25$ である。
従って

$$
\Delta\tau=\Delta t/\gamma=8.00\,\mathrm{ns}
$$

である。
:::
::::

::::exercise{level=basic}
$B_z=0.20\,\mathrm T$、$\vb{E}=0$ の場で $F^{12},F_{12},F^{21}$ を求めよ。

:::solution
$\epsilon_{123}=+1$ なので $F^{12}=B_z=0.20\,\mathrm T$ である。
空間添字を下げても符号は変わらず $F_{12}=0.20\,\mathrm T$、反対称性から $F^{21}=-0.20\,\mathrm T$ となる。
:::
::::

::::exercise{level=basic}
$\vb{E}=0,\ B_z=0.100\,\mathrm T$ の系に対して $+x$ 方向へ $0.60c$ で動く観測者の場を求めよ。

:::solution
$\gamma=1.25$、

$$
\hat{\vb{x}}\times\hat{\vb{z}}=-\hat{\vb{y}}
$$

である。
従って $E'_y=-\gamma v B=-2.25\times10^7\,\mathrm{V/m}$、

$$
B'_z=\gamma B=0.125\,\mathrm T
$$

で、他の成分は零である。
:::
::::

::::exercise{level=basic}
$\vb{E}=3.00\times10^5\hat{\vb{x}}\,\mathrm{V/m}$、$\vb{B}=0.100\hat{\vb{z}}\,\mathrm T$ で電場が消える系の速度を求めよ。

:::solution
直交して $E<cB$ なので

$$
\vb{v}=\vb{E}\times\vb{B}/B^2=-3.00\times10^6\hat{\vb{y}}\,\mathrm{m/s}
$$

である。
$\vb{v}\times\vb{B}=-\vb{E}$ となることも確認できる。
:::
::::

::::exercise{level=advanced}
共変な Maxwell 方程式から電荷保存が必要であることを添字で証明せよ。

:::solution
第一式を $\partial_\nu$ で微分すると

$$
\mu_0\partial_\nu J^\nu=\partial_\nu\partial_\mu F^{\nu\mu}
$$

である。
ダミー添字を交換し偏微分を入れ替えると右辺は元の負になるので零である。
:::
::::

::::exercise{level=advanced}
$E>cB$ の直交場で、純電場の系における電場の大きさを求めよ。

:::solution
$\vb{v}=c^2\vb{E}\times\vb{B}/E^2$ なら速度は $c$ より小さく、変換で $\vb{B}'=0$ となる。
不変量 $B^2-E^2/c^2=-E'^2/c^2$ から $E'=\sqrt{E^2-c^2B^2}$ である。
:::
::::

::::exercise{level=advanced}
一定電場によって静止から $p=mc$ に達する時刻、速さ、移動距離を求めよ。
$qE>0$ とする。

:::solution
$p=qEt$ なので $t=mc/(qE)$ である。
$\mathcal E=\sqrt2mc^2$ より $u=c/\sqrt2$、仕事 $qEx=\mathcal E-mc^2$ より $x=(\sqrt2-1)mc^2/(qE)$ となる。
:::
::::

::::exercise{level=advanced}
磁場中で $p_\perp,p_\parallel$ をもつ粒子の螺旋の一周の軸方向変位を求め、エネルギーへの依存を説明せよ。

:::solution
周期は $2\pi\gamma_um/(|q|B)$、軸方向速度は $p_\parallel/(\gamma_um)$ である。
積は $2\pi p_\parallel/(|q|B)$ となる。
$p_\parallel$ を固定して横運動量を増すと周期は延びるが、軸速度が同じ割合で下がるため変位は変わらない。
:::
::::

::::exercise{level=advanced}
四元電流が時間的な場合と空間的な場合の違いを、$\vb{J}$ が零になる系または $\rho$ が零になる系の存在で説明せよ。
$J^\mu$ 自身は非零とする。

:::solution
不変量は $J_\mu J^\mu=J^2-c^2\rho^2$ である。
$\vb{v}$ を電流軸に沿って選ぶと

$$
\vb{J}'_\parallel=\gamma(\vb{J}-\rho\vb{v})
$$

である。
時間的なら $\rho\ne0$ で $|\vb{J}/\rho|<c$ なので $\vb{v}=\vb{J}/\rho$ とすれば電流は零になる。

空間的なら $J\ne0$ かつ $c|\rho|<J$ である。
$\rho'=\gamma(\rho-\vb{v}\cdot\vb{J}/c^2)$ に対し $\vb{v}=c^2\rho\vb{J}/J^2$ とすれば亜光速で密度を零にできる。

零ノルムの場合には必要速度の大きさが $c$ になる。
非零の四元電流を Lorentz 変換で零にすることもできないため、どちらの系も存在しない。
:::
::::

::::exercise{level=advanced}
光の波数四元ベクトルと観測者の固有時間から、任意の方向に進む波の Doppler 効果と光行差を求めよ。
観測者の移動方向を $+x$ とし、実験室で波と $+x$ のなす角を $\theta$ とする。

:::solution
位相を $k_\mu x^\mu$ と書くと $k^\mu=(\omega/c,\vb{k})$、$k_\mu k^\mu=0$ である。
四元速度 $W^\mu$ の観測者が自身の時計で測る角周波数は位相の減少率 $-k_\mu W^\mu$ なので、

$$
\omega'=\gamma\omega(1-\beta\cos\theta).
$$

空間成分を変換すると $k'_x=\gamma(\omega/c)(\cos\theta-\beta)$、$k'_y=(\omega/c)\sin\theta$ である。
$|\vb{k}'|=\omega'/c$ を用いると

$$
\cos\theta'=\frac{\cos\theta-\beta}{1-\beta\cos\theta},\qquad \sin\theta'=\frac{\sin\theta}{\gamma(1-\beta\cos\theta)}.
$$

実験室で横から届く $\theta=\pi/2$ の光なら $\omega'=\gamma\omega$ である。
どの系で角度を定義したかを変えると、「横から」の条件自体が変わる。
:::
::::

::::exercise{level=advanced}
相対論的な一様電場の軌道を、四元運動方程式から固有時間を変数として直接解け。
最初に静止とし、$\vb{E}=E\hat{\vb{x}}$、$qE>0$ とする。

:::solution
$U_0=-U^0,\ U_1=U^1$ と $F^{01}=E/c$ から、$\alpha=qE/(mc)$ とすれば

$$
\frac{\dd U^0}{\dd\tau}=\alpha U^1,\qquad \frac{\dd U^1}{\dd\tau}=\alpha U^0.
$$

$U^0(0)=c,\ U^1(0)=0$ より $U^0=c\cosh(\alpha\tau)$、$U^1=c\sinh(\alpha\tau)$ である。
$\dd t/\dd\tau=U^0/c,\ \dd x/\dd\tau=U^1$ を積分すると

$$
t=\frac{\sinh(\alpha\tau)}{\alpha},\qquad x=\frac c\alpha[\cosh(\alpha\tau)-1].
$$

$\cosh^2-\sinh^2=1$ より双曲線軌道が得られる。
固有時間発展は $(U^0,U^1)$ 平面内を Lorentz ブーストで進む運動であり、速さは $c\tanh(\alpha\tau)$ のため $c$ を超えない。
:::
::::

::::exercise{level=advanced}
質量 $m$、電荷 $q>0$ の粒子が、一定な交差場 $\vb E=E\hat{\vb x}$、$\vb B=(E/c)\hat{\vb y}$ の中で静止から動き出す。$a=qE/(mc)$ とし、固有時間が $a\tau=2$ に達した瞬間の Lorentz 因子と三速度の二成分を求めよ。仕事によって得た運動エネルギーとも照合せよ。

:::solution
$U^0/c=1+(a\tau)^2/2=3$、$U^x/c=a\tau=2$、$U^z/c=(a\tau)^2/2=2$ である。三速度は空間成分を $\gamma=3$ で割るので $v_x=v_z=2c/3$ となり、速さの二乗は $8c^2/9$ となる。電場の仕事は $qEx=qEac\tau^2/2=mc^2(a\tau)^2/2=2mc^2$ であり、$(\gamma-1)mc^2$ に一致する。磁場は軌道を曲げても仕事をしない。
:::
::::

## 確認テスト

::::quiz{answer="1.25" tol="1e-8"}
速さ $0.60c$ の Lorentz 因子を答えよ。

:::explanation

$$
\gamma=(1-0.60^2)^{-1/2}=1/0.80=1.25
$$

である。
:::
::::

::::quiz
本章の $F^{0i}=E_i/c$ と計量で $F_{0i}$ を選べ。

- [ ] $E_i/c$
- [x] $-E_i/c$
- [ ] $-cE_i$

:::explanation
時間の添字を下げると $\eta_{00}=-1$ が一度掛かり、空間添字からの符号変化はない。
:::
::::

::::quiz{answer="0.5" tol="1e-8"}
$+x$ へ進む光を $+x$ に $0.60c$ で追う観測者にとって、周波数は元の何倍か。

:::explanation
波の四元波数を変換すると

$$
\omega'/\omega=\gamma(1-\beta)=1.25\times0.40=0.50
$$

である。
:::
::::

::::quiz
どの慣性系でも磁場と電場がともに非零となる場を選べ。

- [x] 非零で直交し $E=cB$ の場
- [ ] 非零の純静電場
- [ ] 非零の純静磁場

:::explanation
二つの不変量が零であるため、片方を消せば他方も零となる。
Lorentz 変換は可逆であり、非零テンソルを零へ変換できない。
:::
::::

::::quiz
自由電荷が静止している慣性系と、その電荷だけが動いて見える慣性系で異なりうる量を選べ。

- [ ] 個々の粒子の電荷
- [x] 一定の実験室体積で測る電荷密度
- [ ] 四元電流の自己縮約

:::explanation
同時に測る体積は観測者によって異なり、密度 $J^0/c$ は空間成分と混ざる。
個々の電荷はスカラーであり、$J_\mu J^\mu$ は完全縮約なので慣性系で不変である。
:::
::::

::::quiz{answer="4.63" tol="0.01"}
速さ $0.80c$ の等速点電荷について、同じ瞬時距離での横方向の電場を縦方向の電場で割れ。

:::explanation
静止場の大きさに対し横方向は $\gamma$ 倍、縦方向は $\gamma^{-2}$ 倍なので比は $\gamma^3$ である。
$\gamma=5/3$ より $125/27=4.63$ となる。
:::
::::

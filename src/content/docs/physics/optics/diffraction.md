---
title: 回折理論
description: スカラー波の境界積分から回折を導き、近距離と遠距離の回折像および分光分解能を求める。
notation:
  - symbol: 'U'
    meaning: スカラー場の複素振幅
  - symbol: 'G'
    meaning: 外向き Helmholtz Green 関数
  - symbol: 'k,\lambda'
    meaning: 観測媒質内の波数と波長
  - symbol: '\boldsymbol\rho,\boldsymbol\rho'''
    meaning: 観測面と開口面の横座標
  - symbol: '\operatorname{sinc}u'
    meaning: 原点で 1 とした sin u / u
---

## この章の目標

- Green の恒等式を用いて回折積分を導く。
- Fresnel 近似と Fraunhofer 近似の位相誤差を見積もる。
- 開口の形と回折像の零点を関係づける。
- 格子の分散と有限開口による線幅を比較する。

## Huygens–Fresnel の原理

### 波面上の振幅から次の面へ
一様な媒質で $k=2\pi/\lambda>0$、時間因子を $e^{-i\omega t}$ とする。偏光の変化を無視できる場合、特定の横電場成分の複素振幅をスカラー $U$ として $\nabla^2U+k^2U=0$ を解く。開口が波長程度で金属の縁の接線場が重要なときには、この置き換えに誤差が出る。:term[Huygens–Fresnel の原理]{en="Huygens–Fresnel principle" yomi="ほいへんすふれねるのげんり"} では、通過した波面の各点から出る二次波を複素振幅で足す。斜方向への放射には、境界条件に応じた傾斜因子が掛かる。

距離 $R=|\vb r-\vb r'|$ の球面因子 $G(\vb r,\vb r')=e^{ikR}/(4\pi R)$ は、$\vb r\ne\vb r'$ で Helmholtz 方程式を満たす（:eqref[eq-optics-spherical]）。中心の小球面で $\partial_R G\simeq-1/(4\pi R^2)$ の流出を積分すると $-1$ である。従って $\vb r'$ の変数に関して $(\nabla'^2+k^2)G=-\delta(\vb r'-\vb r)$。$G$ が二次波に対応し、境界での値と勾配をどう与えるかが次節の問題となる。

:::example[二点での打ち消し]
同じ距離と同じ傾斜因子で届く等振幅の二次波が重なる点では、出発位相が $\pi$ 異なれば複素和はゼロである。光強度を先に足す計算ではこの暗点を得られない。
:::
:::example[遠方の小開口]
開口の端から観測点への距離差が $\lambda/100$ 以下なら位相差は $2\pi/100$ 以下である。照明が一様で傾斜因子と距離の逆数が開口全体でほぼ一定なら、面積を二倍にすると観測点で振幅がほぼ二倍になる。
:::

### スカラー近似の中で使う二種類の源
式 :eqref[eq-optics-spherical] の外向き解は、点から等方的に出る Green 関数に対応する。法線で微分した $\partial_{n'}G$ は、この点を法線に沿って少しずらした二つの点源を、反対の符号で足した極限である。観測領域の外向き法線を $-\hat{\vb z}$ とし、$\theta_o$ を $+\hat{\vb z}$ から測る。$kR\gg1$ なら微分は球面因子に $ik\cos\theta_o$ を掛けたものに近づく。したがって傾斜因子は、開口場を値だけで与えるか勾配も与えるかで変わる。後続の積分に二次波の波面と振幅を両方入れる理由である。

:::example[点の位相差で暗点を探す]
二つの等振幅点が同相で、観測点までの距離差が $0.25\,\mathrm{\mu m}$ とする。波長 $0.50\,\mathrm{\mu m}$ なら位相差は $\pi$ なので、振幅因子の差を無視できる遠方では暗点になる。波長を $1.0\,\mathrm{\mu m}$ に変えると差は $\pi/2$ で、二つの照度の和を観測する。
:::

### 複素振幅の加法で暗点の深さを決める

二つの二次源が同じ初期位相で、観測点までの距離がそれぞれ $R$ と $R+\Delta$ とする。傾斜因子を等しいとみなせば球面因子の加法から、第一の場で規格化した合成振幅は $1+q e^{ik\Delta}$、$q=R/(R+\Delta)$ である。従って $\Delta$ が半波長の奇数倍でも強度は $(1-q)^2$ だけ残る。逆に、二源の強さを補って $q=1$ にすれば暗点を深くできるが、光路誤差が $\delta$ なら暗点の近くで $2[1-\cos(k\delta)]\simeq(k\delta)^2$ が残る。振幅不均衡と光路誤差を独立に測る意味がここにある。

:::example[遠方の半周期]
等しい源の半波長の距離差で打ち消しを調べる。

入力は第一の距離 $R=\frac{1}{10}\,\mathrm m$、距離差 $\Delta\simeq 2.5\times 10^{-7}\,\mathrm m$、波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$ である。

逆距離による振幅低下を含める。

$$
\begin{aligned}
q&=\frac{R}{R + \Delta},\\
q&\simeq \frac{\left(\frac{1}{10}\right)}{\left(\frac{1}{10}\right) + \left(2.5\times 10^{-7}\right)}\simeq 1.
\end{aligned}
$$

光路の差だけが位相差を作る。

$$
\begin{aligned}
\phi&=\frac{2 \pi \Delta}{\lambda},\\
\phi/(\mathrm{rad})&\simeq \frac{2 \cdot \pi \cdot \left(2.5\times 10^{-7}\right)}{\left(5\times 10^{-7}\right)}\simeq 3.1416.
\end{aligned}
$$

振幅を足してから二乗する。

$$
\begin{aligned}
I/I_1&=q^{2} + 2 q \cos{\left(\phi \right)} + 1,\\
I/I_1&\simeq \left(0.9999975\right)^{2} + 2 \cdot \left(0.9999975\right) \cdot \cos{\left(\left(3.141592654\right) \right)} + 1\simeq 6.25\times 10^{-12}.
\end{aligned}
$$

同じ振幅で位相だけを調整したときの上限を求める。

$$
\begin{aligned}
I_{\max}/I_1&=\left(q + 1\right)^{2},\\
I_{\max}/I_1&\simeq \left(\left(0.9999975\right) + 1\right)^{2}\simeq 4.
\end{aligned}
$$

二源の暗い位置を明るい位置と同じ規格で比較する。

$$
\begin{aligned}
I/I_{\max}&=\frac{I/I_1}{I_{\max}/I_1},\\
I/I_{\max}&\simeq \frac{\left(6.24996875\times 10^{-12}\right)}{\left(3.99999\right)}\simeq 1.5625\times 10^{-12}.
\end{aligned}
$$

距離の差による振幅不均衡だけが残る。
:::

:::example[同じ位相差で距離を縮める]
源から数波長の位置で、同じ距離差を持つ素片を比較する。

入力は第一の距離 $R\simeq 5\times 10^{-6}\,\mathrm m$、距離差 $\Delta\simeq 2.5\times 10^{-7}\,\mathrm m$、波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$ である。

逆距離による振幅低下を含める。

$$
\begin{aligned}
q&=\frac{R}{R + \Delta},\\
q&\simeq \frac{\left(5\times 10^{-6}\right)}{\left(5\times 10^{-6}\right) + \left(2.5\times 10^{-7}\right)}=\frac{20}{21}.
\end{aligned}
$$

光路の差だけが位相差を作る。

$$
\begin{aligned}
\phi&=\frac{2 \pi \Delta}{\lambda},\\
\phi/(\mathrm{rad})&\simeq \frac{2 \cdot \pi \cdot \left(2.5\times 10^{-7}\right)}{\left(5\times 10^{-7}\right)}\simeq 3.1416.
\end{aligned}
$$

振幅を足してから二乗する。

$$
\begin{aligned}
I/I_1&=q^{2} + 2 q \cos{\left(\phi \right)} + 1,\\
I/I_1&\simeq \left(\frac{20}{21}\right)^{2} + 2 \cdot \left(\frac{20}{21}\right) \cdot \cos{\left(\left(3.141592654\right) \right)} + 1=\frac{1}{441}.
\end{aligned}
$$

同じ振幅で位相だけを調整したときの上限を求める。

$$
\begin{aligned}
I_{\max}/I_1&=\left(q + 1\right)^{2},\\
I_{\max}/I_1&= \left(\left(\frac{20}{21}\right) + 1\right)^{2}=\frac{1681}{441}.
\end{aligned}
$$

二源の暗い位置を明るい位置と同じ規格で比較する。

$$
\begin{aligned}
I/I_{\max}&=\frac{I/I_1}{I_{\max}/I_1},\\
I/I_{\max}&= \frac{\left(\frac{1}{441}\right)}{\left(\frac{1681}{441}\right)}=\frac{1}{1681}.
\end{aligned}
$$

理想化した球面因子のモデルでも、距離の比を一と置くと暗点を深く見積もる。
:::

:::example[光路差をわずかに増やす]
遠方へ戻して、距離差を半波長からずらす。

入力は第一の距離 $R=\frac{1}{10}\,\mathrm m$、距離差 $\Delta\simeq 2.6\times 10^{-7}\,\mathrm m$、波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$ である。

逆距離による振幅低下を含める。

$$
\begin{aligned}
q&=\frac{R}{R + \Delta},\\
q&\simeq \frac{\left(\frac{1}{10}\right)}{\left(\frac{1}{10}\right) + \left(2.6\times 10^{-7}\right)}\simeq 1.
\end{aligned}
$$

光路の差だけが位相差を作る。

$$
\begin{aligned}
\phi&=\frac{2 \pi \Delta}{\lambda},\\
\phi/(\mathrm{rad})&\simeq \frac{2 \cdot \pi \cdot \left(2.6\times 10^{-7}\right)}{\left(5\times 10^{-7}\right)}\simeq 3.2673.
\end{aligned}
$$

振幅を足してから二乗する。

$$
\begin{aligned}
I/I_1&=q^{2} + 2 q \cos{\left(\phi \right)} + 1,\\
I/I_1&\simeq \left(0.9999974\right)^{2} + 2 \cdot \left(0.9999974\right) \cdot \cos{\left(\left(3.26725636\right) \right)} + 1\simeq 0.015771.
\end{aligned}
$$

同じ振幅で位相だけを調整したときの上限を求める。

$$
\begin{aligned}
I_{\max}/I_1&=\left(q + 1\right)^{2},\\
I_{\max}/I_1&\simeq \left(\left(0.9999974\right) + 1\right)^{2}\simeq 4.
\end{aligned}
$$

二源の暗い位置を明るい位置と同じ規格で比較する。

$$
\begin{aligned}
I/I_{\max}&=\frac{I/I_1}{I_{\max}/I_1},\\
I/I_{\max}&\simeq \frac{\left(0.01577055637\right)}{\left(3.9999896\right)}\simeq 0.0039426.
\end{aligned}
$$

この条件では振幅不均衡より位相の誤差が大きく効く。
:::

## Kirchhoff の回折積分

### Green の定理と開口での近似
体積 $V$ の内点に観測点を置く。$\divg(U\grad'G-G\grad'U)=U\nabla'^2G-G\nabla'^2U$ を積分し、外向き法線 $\hat{\vb n}'$ を用いると
$$
U(\vb r)=\int_{\partial V}(G\partial_{n'}U-U\partial_{n'}G)\dd S'.
\label{eq-optics-green-surface}
$$
上半空間 $z>0$ と大半球を境界にする。外向き場の遠方条件 $\partial_R U-ikU=o(R^{-1})$、$U=O(R^{-1})$ の下で、球面上の二項の先頭次数は相殺されるので、半径を広げると寄与はゼロになる。

平面 $z'=0$ の法線を $-\hat{\vb z}$ とする。:term[Kirchhoff の回折積分]{en="Kirchhoff diffraction integral" yomi="きるひほっふのかいせきせきぶん"} では、開口上で $U$ と $\partial_{n'}U$ をスクリーンがない照明場で置き換え、遮蔽部で両方ゼロとする。平面照明で $\partial_{n'}U\simeq-ik\cos\theta_i U$、観測方向の法線に対する角を $\theta_o$ として $\partial_{n'}G=(ik-1/R)G\cos\theta_o$ である。$kR\gg1$ の主項は
$$
U(\vb r)\simeq\int_A U(\vb r')\frac{e^{ikR}}{i\lambda R}
\frac{\cos\theta_i+\cos\theta_o}{2}\dd^2\rho' .
$$
有限開口で値と法線微分を独立に指定すると一般には両立しない。特に縁の近傍では両量が同時に急変するという近似を信用できない。

:::example[正面と斜方向]
正面照明、正面観測なら係数は 1。観測が $\theta_o=60^\circ$ なら係数は $3/4$ で、振幅の比較と強度の比較を混同してはならない。
:::
:::example[面の向き]
$G$ の法線微分に用いる法線は観測領域から外向きである。$+\hat{\vb z}$ を使う場合は式 :eqref[eq-optics-green-surface] の平面積分全体の符号が反転する。
:::

### 点照明での入射角
開口の点 $\vb r'$ から軸上の点光源までの距離を $R_i$ とする。入射波は $U_i=A_i e^{ikR_i}/R_i$ なので、観測領域から外向きの法線方向への微分は $-(ik-1/R_i)U_i\cos\theta_i$。式 :eqref[eq-optics-green-surface] から、近似する前の積分の括弧には $-[(ik-1/R_i)\cos\theta_i+(ik-1/R)\cos\theta_o]G U_i$ が入る。両距離が波長に比べて大きければ $1/R_i,1/R$ を落として同じ傾斜因子を得るが、積分の位相には $k(R_i+R)$ を残す。

:::example[傾斜因子の両側を評価する]
点源からの入射角が $60^\circ$、観測角が $0^\circ$ の開口素片では、Kirchhoff の遠方の傾斜因子は $(1/2+1)/2=3/4$。観測角も $60^\circ$ にすれば $1/2$ になる。同じ素片の比較でも入射角だけでは係数は決まらない。
:::

### 点照明の近接項を二つに分ける

点照明の一素片の Kirchhoff 核では、二つの法線微分から $ik(\cos i+\cos o)$ と $-(\cos i/R_i+\cos o/R_o)$ が現れる。実数の角と距離なら二つは直交位相である。$\cos i+\cos o>0$ とすれば、主項に対する近接項の比の大きさは $(\cos i/R_i+\cos o/R_o)/(k(\cos i+\cos o))$。従って各点の核での相対振幅の補正はこの量の二乗から始まる一方、位相補正は一次から始まる。積分全体の相対誤差はこれより大きくなることがある。主項を足し合わせたときにほぼ打ち消し合う位置では、小さい補正項でも強度を決めるためである。

:::example[遠い点源と遠い観測点]
波長に比べて両側の距離が長い素片を評価する。

入力は波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、点源から素片までの距離 $R_i=\frac{1}{10}\,\mathrm m$、素片から観測点までの距離 $R_o=\frac{1}{5}\,\mathrm m$、入射角 $\theta_i\simeq 1.0472\,\mathrm{rad}$、観測角 $\theta_o\simeq 0.5236\,\mathrm{rad}$ である。

入射側と観測側の余弦を足す。

$$
\begin{aligned}
C&=\frac{\cos{\left(\theta_i \right)}}{2} + \frac{\cos{\left(\theta_o \right)}}{2},\\
C&\simeq \frac{\cos{\left(\left(1.047197551\right) \right)}}{2} + \frac{\cos{\left(\left(0.5235987756\right) \right)}}{2}\simeq 0.68301.
\end{aligned}
$$

両方の距離の逆数を法線の投影で重み付けする。

$$
\begin{aligned}
r&=\frac{\lambda \left(\frac{\cos{\left(\theta_o \right)}}{R_o} + \frac{\cos{\left(\theta_i \right)}}{R_i}\right)}{4 \pi C},\\
r&\simeq \frac{\left(5\times 10^{-7}\right) \cdot \left(\frac{\cos{\left(\left(0.5235987756\right) \right)}}{\left(\frac{1}{5}\right)} + \frac{\cos{\left(\left(1.047197551\right) \right)}}{\left(\frac{1}{10}\right)}\right)}{4 \cdot \pi \cdot \left(0.6830127019\right)}\simeq 5.4352\times 10^{-7}.
\end{aligned}
$$

元の括弧は主項に $1+ir$ を掛けたものなので、符号は正である。

$$
\begin{aligned}
\delta\phi&=\operatorname{atan}{\left(r \right)},\\
\delta\phi/(\mathrm{rad})&\simeq \operatorname{atan}{\left(\left(5.435242385\times 10^{-7}\right) \right)}\simeq 5.4352\times 10^{-7}.
\end{aligned}
$$

近接項は位相が四分の一周期異なるため二乗して加える。

$$
\begin{aligned}
m&=\sqrt{r^{2} + 1},\\
m&\simeq \sqrt{\left(\left(5.435242385\times 10^{-7}\right)\right)^{2} + 1}\simeq 1.
\end{aligned}
$$

近接項は小さく、傾斜因子だけを残した核を使える。
:::

:::example[点源が開口に近い]
観測距離は長いが点照明が近い条件を取る。

入力は波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、点源から素片までの距離 $R_i\simeq 2\times 10^{-6}\,\mathrm m$、素片から観測点までの距離 $R_o=\frac{1}{5}\,\mathrm m$、入射角 $\theta_i\simeq 1.0472\,\mathrm{rad}$、観測角 $\theta_o\simeq 0.5236\,\mathrm{rad}$ である。

入射側と観測側の余弦を足す。

$$
\begin{aligned}
C&=\frac{\cos{\left(\theta_i \right)}}{2} + \frac{\cos{\left(\theta_o \right)}}{2},\\
C&\simeq \frac{\cos{\left(\left(1.047197551\right) \right)}}{2} + \frac{\cos{\left(\left(0.5235987756\right) \right)}}{2}\simeq 0.68301.
\end{aligned}
$$

両方の距離の逆数を法線の投影で重み付けする。

$$
\begin{aligned}
r&=\frac{\lambda \left(\frac{\cos{\left(\theta_o \right)}}{R_o} + \frac{\cos{\left(\theta_i \right)}}{R_i}\right)}{4 \pi C},\\
r&\simeq \frac{\left(5\times 10^{-7}\right) \cdot \left(\frac{\cos{\left(\left(0.5235987756\right) \right)}}{\left(\frac{1}{5}\right)} + \frac{\cos{\left(\left(1.047197551\right) \right)}}{\left(2\times 10^{-6}\right)}\right)}{4 \cdot \pi \cdot \left(0.6830127019\right)}\simeq 0.014564.
\end{aligned}
$$

元の括弧は主項に $1+ir$ を掛けたものなので、符号は正である。

$$
\begin{aligned}
\delta\phi&=\operatorname{atan}{\left(r \right)},\\
\delta\phi/(\mathrm{rad})&\simeq \operatorname{atan}{\left(\left(0.01456394033\right) \right)}\simeq 0.014563.
\end{aligned}
$$

近接項は位相が四分の一周期異なるため二乗して加える。

$$
\begin{aligned}
m&=\sqrt{r^{2} + 1},\\
m&\simeq \sqrt{\left(0.01456394033\right)^{2} + 1}\simeq 1.0001.
\end{aligned}
$$

観測側だけを遠くしても入射側の近接項は消えない。
:::

:::example[入射が斜めで観測が近い]
同じ波長で距離の短い側を観測側に移す。

入力は波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、点源から素片までの距離 $R_i=\frac{1}{10}\,\mathrm m$、素片から観測点までの距離 $R_o\simeq 2\times 10^{-6}\,\mathrm m$、入射角 $\theta_i\simeq 1.0472\,\mathrm{rad}$、観測角 $\theta_o\simeq 0.5236\,\mathrm{rad}$ である。

入射側と観測側の余弦を足す。

$$
\begin{aligned}
C&=\frac{\cos{\left(\theta_i \right)}}{2} + \frac{\cos{\left(\theta_o \right)}}{2},\\
C&\simeq \frac{\cos{\left(\left(1.047197551\right) \right)}}{2} + \frac{\cos{\left(\left(0.5235987756\right) \right)}}{2}\simeq 0.68301.
\end{aligned}
$$

両方の距離の逆数を法線の投影で重み付けする。

$$
\begin{aligned}
r&=\frac{\lambda \left(\frac{\cos{\left(\theta_o \right)}}{R_o} + \frac{\cos{\left(\theta_i \right)}}{R_i}\right)}{4 \pi C},\\
r&\simeq \frac{\left(5\times 10^{-7}\right) \cdot \left(\frac{\cos{\left(\left(0.5235987756\right) \right)}}{\left(2\times 10^{-6}\right)} + \frac{\cos{\left(\left(1.047197551\right) \right)}}{\left(\frac{1}{10}\right)}\right)}{4 \cdot \pi \cdot \left(0.6830127019\right)}\simeq 0.025225.
\end{aligned}
$$

元の括弧は主項に $1+ir$ を掛けたものなので、符号は正である。

$$
\begin{aligned}
\delta\phi&=\operatorname{atan}{\left(r \right)},\\
\delta\phi/(\mathrm{rad})&\simeq \operatorname{atan}{\left(\left(0.02522533897\right) \right)}\simeq 0.02522.
\end{aligned}
$$

近接項は位相が四分の一周期異なるため二乗して加える。

$$
\begin{aligned}
m&=\sqrt{r^{2} + 1},\\
m&\simeq \sqrt{\left(0.02522533897\right)^{2} + 1}\simeq 1.0003.
\end{aligned}
$$

近い側の余弦に重みが移る。近接項の大きさは最短距離だけでは決まらない。
:::

## Rayleigh–Sommerfeld の回折理論

境界で場の値だけを与える半空間問題を考える。実際の遮蔽板で境界値を求める問題と、与えた境界値を伝播する問題を分ける。後者では像点を $z<0$ に置き、$G_D=G(\vb r,\vb r')-G(\vb r_{\mathrm{im}},\vb r')$ とすると平面上で $G_D=0$ である。式 :eqref[eq-optics-green-surface] に入れると平面上で $\partial_{n'}G_D=2(ik-1/R)G(z/R)$ であり、
$$
U(\boldsymbol\rho,z)=\int_{\R^2} U(\boldsymbol\rho',0)
\frac{e^{ikR}}{i\lambda R}\frac{z}{R}
\left(1-\frac{1}{ikR}\right)\dd^2\rho',\qquad
R=\sqrt{z^2+|\boldsymbol\rho-\boldsymbol\rho'|^2}.
\label{eq-optics-rayleigh-sommerfeld}
$$
これが第一の :term[Rayleigh–Sommerfeld 回折積分]{en="Rayleigh–Sommerfeld diffraction integral" yomi="れいりーぞんまーふぇるとかいせきせきぶん"} である。法線微分 $g=\partial_{n'}U$ を与える場合は、像を足した Green 関数 $G_N$ の法線微分が平面で消える。従って第二の積分は $U(\vb r)=2\int_{\R^2}G(\vb r,\vb r')g(\boldsymbol\rho')\dd^2\rho'$ となる。開口で値と勾配に同じ Kirchhoff の近似を入れると、第一の積分と第二の積分の算術平均が Kirchhoff の式になる。第一の積分には値の項が 2 倍で入り、第二の積分には勾配の項が 2 倍で入るためである。Dirichlet の式と Neumann の式の入力データは異なる。開口上のデータを入射波で代用するなら、その段階で近似が入る。

:::example[正面に限った一致]
正面照明と近軸観測を選び、$kR\gg1$ とすると傾斜因子 $\cos\theta_o$ は 1 に近づく。したがって Kirchhoff の積分と同じ近軸核になる。
:::
:::example[平面波の境界データ]
境界値が一定なら並進対称な外向き解は $U(z)=U(0)e^{ikz}$ である。反対の符号を選んだ解は半空間の奥から波が来ることを意味する。
:::

:::example[最も近い点での補正の尺度]
波長 $500\,\mathrm{nm}$、開口点までの距離 $R=5\,\mathrm{\mu m}$ なら、核で落とす $1/(ikR)$ の絶対値は $1/(20\pi)\simeq0.016$。ただし積分の主項が干渉で消える位置では、補正の相対効果はこの数より大きくなり得る。
:::

### 円形ディスク状の境界値で正面を厳密に積分する
半径 $a$ の内側で $U(\boldsymbol\rho',0)=U_0$、外側でゼロとする Dirichlet データを考える。観測を軸上に置き、$R=\sqrt{z^2+r^2}$ で $r\dd r=R\dd R$。式 :eqref[eq-optics-rayleigh-sommerfeld] は
$$
U(0,z)=-U_0 z\int_z^{\sqrt{z^2+a^2}}
\frac{\dd}{\dd R}\left(\frac{e^{ikR}}R\right)\dd R
=U_0\left[e^{ikz}-\frac{z}{\sqrt{z^2+a^2}}e^{ik\sqrt{z^2+a^2}}\right].
$$
この式はスカラー境界値問題内での結果であり、実在のスクリーンがこのデータを正確に作ることまでは仮定しない。$a/z\ll1$ として分数を 1 とし、二次位相を残すと円孔の Fresnel 式が再現される。固定の $z>0$ で $a\to\infty$ とすると第 2 項は $z/a$ の振幅で消え、平面波に収束する。これに対し Fresnel 式のまま $a\to\infty$ とすると、仮定 $a/z\ll1$ を超えてしまう。

### 角スペクトルで境界値の解の一意性をみる

横方向に Fourier 変換すると、横波数 $K$ ごとの縦振幅は $\partial_z^2 U+(k^2-K^2)U=0$ を満たす。$K<k$ なら外向きの $e^{i\sqrt{k^2-K^2}z}$、$K>k$ なら減衰する $e^{-\sqrt{K^2-k^2}z}$ だけを残す。$K=k$ は有界な定数解を取る。各振幅は $z=0$ の値で決まるので、同じ境界値を持つ二つの解の差はゼロである。これは平面で Fourier 変換できる場のクラスでの一意性である。:eqref[eq-optics-rayleigh-sommerfeld] は同じ方程式と境界に対する積分解であり、外向きの枝は一致する。

場強度の比と面を横切る流れの比を区別するには、スカラー流 $j_z=\Im(U^*\partial_z U)$ を使える。Helmholtz 方程式とその共役の差から $\divg \Im(U^*\grad U)=0$ が従う。伝播成分では $j_z=k_z|U|^2$、単独の指数減衰成分では $\partial_z U$ が $U$ の負の実数倍なので $j_z=0$ である。境界周期がちょうど波長の成分でも $k_z=0$ なので、場強度が一定のまま法線方向の流れはゼロとなる。

:::example[境界のゆっくりした変調]
波長の二倍の横周期を与える。

入力は波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、境界の正弦波の周期 $\Lambda\simeq 1\times 10^{-6}\,\mathrm m$、伝播距離 $z\simeq 1\times 10^{-6}\,\mathrm m$ である。

媒質で許される波数の大きさを求める。

$$
\begin{aligned}
k&=\frac{2 \pi}{\lambda},\\
k/(\mathrm{m^{-1}})&\simeq \frac{2 \cdot \pi}{\left(5\times 10^{-7}\right)}\simeq 1.2566\times 10^{7}.
\end{aligned}
$$

データの横位相周期を波数へ換算する。

$$
\begin{aligned}
K&=\frac{2 \pi}{\Lambda},\\
K/(\mathrm{m^{-1}})&\simeq \frac{2 \cdot \pi}{\left(1\times 10^{-6}\right)}\simeq 6.2832\times 10^{6}.
\end{aligned}
$$

正の実部または正の虚部の枝を取る。

$$
\begin{aligned}
k_z&=\sqrt{- K^{2} + k^{2}},\\
k_z/(\mathrm{m^{-1}})&\simeq \sqrt{- \left(\left(6.283185307\times 10^{6}\right)\right)^{2} + \left(\left(1.256637061\times 10^{7}\right)\right)^{2}}\simeq 1.0883\times 10^{7}.
\end{aligned}
$$

独立な一つの横成分を進める。

$$
\begin{aligned}
T&=e^{i k_z z},\\
T&\simeq e^{i \cdot \left(1.088279619\times 10^{7}\right) \cdot \left(1\times 10^{-6}\right)}\simeq -0.11254-0.99365i.
\end{aligned}
$$

伝達因子の絶対値を二乗すると、境界に対する場強度の比が得られる。

$$
\begin{aligned}
Q&=\left|{T}\right|^{2},\\
Q&\simeq \left|{\left(-0.1125391852-0.9936472874i\right)}\right|^{2}=1.
\end{aligned}
$$

縦方向への位相移動が生じ、二乗の比は一である。
:::

:::example[ちょうど波長の境界周期]
横方向で波数の全てを使う。

入力は波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、境界の正弦波の周期 $\Lambda\simeq 5\times 10^{-7}\,\mathrm m$、伝播距離 $z\simeq 1\times 10^{-6}\,\mathrm m$ である。

媒質で許される波数の大きさを求める。

$$
\begin{aligned}
k&=\frac{2 \pi}{\lambda},\\
k/(\mathrm{m^{-1}})&\simeq \frac{2 \cdot \pi}{\left(5\times 10^{-7}\right)}\simeq 1.2566\times 10^{7}.
\end{aligned}
$$

データの横位相周期を波数へ換算する。

$$
\begin{aligned}
K&=\frac{2 \pi}{\Lambda},\\
K/(\mathrm{m^{-1}})&\simeq \frac{2 \cdot \pi}{\left(5\times 10^{-7}\right)}\simeq 1.2566\times 10^{7}.
\end{aligned}
$$

正の実部または正の虚部の枝を取る。

$$
\begin{aligned}
k_z&=\sqrt{- K^{2} + k^{2}},\\
k_z/(\mathrm{m^{-1}})&\simeq \sqrt{- \left(\left(1.256637061\times 10^{7}\right)\right)^{2} + \left(\left(1.256637061\times 10^{7}\right)\right)^{2}}=0.
\end{aligned}
$$

独立な一つの横成分を進める。

$$
\begin{aligned}
T&=e^{i k_z z},\\
T&\simeq e^{i \cdot \left(0\right) \cdot \left(1\times 10^{-6}\right)}=1.
\end{aligned}
$$

伝達因子の絶対値を二乗すると、境界に対する場強度の比が得られる。

$$
\begin{aligned}
Q&=\left|{T}\right|^{2},\\
Q&= \left|{\left(1\right)}\right|^{2}=1.
\end{aligned}
$$

有界解は縦座標に依存しない。
:::

:::example[波長以下の変調]
横周期を半波長に縮める。

入力は波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、境界の正弦波の周期 $\Lambda\simeq 2.5\times 10^{-7}\,\mathrm m$、伝播距離 $z\simeq 1\times 10^{-7}\,\mathrm m$ である。

媒質で許される波数の大きさを求める。

$$
\begin{aligned}
k&=\frac{2 \pi}{\lambda},\\
k/(\mathrm{m^{-1}})&\simeq \frac{2 \cdot \pi}{\left(5\times 10^{-7}\right)}\simeq 1.2566\times 10^{7}.
\end{aligned}
$$

データの横位相周期を波数へ換算する。

$$
\begin{aligned}
K&=\frac{2 \pi}{\Lambda},\\
K/(\mathrm{m^{-1}})&\simeq \frac{2 \cdot \pi}{\left(2.5\times 10^{-7}\right)}\simeq 2.5133\times 10^{7}.
\end{aligned}
$$

正の実部または正の虚部の枝を取る。

$$
\begin{aligned}
k_z&=\sqrt{- K^{2} + k^{2}},\\
k_z/(\mathrm{m^{-1}})&\simeq \sqrt{- \left(\left(2.513274123\times 10^{7}\right)\right)^{2} + \left(\left(1.256637061\times 10^{7}\right)\right)^{2}}\simeq 2.1766\times 10^{7}i.
\end{aligned}
$$

独立な一つの横成分を進める。

$$
\begin{aligned}
T&=e^{i k_z z},\\
T&\simeq e^{i \cdot \left(2.176559237\times 10^{7}i\right) \cdot \left(1\times 10^{-7}\right)}\simeq 0.11343.
\end{aligned}
$$

伝達因子の絶対値を二乗すると、境界に対する場強度の比が得られる。

$$
\begin{aligned}
Q&=\left|{T}\right|^{2},\\
Q&\simeq \left|{\left(0.1134311497\right)}\right|^{2}\simeq 0.012867.
\end{aligned}
$$

横の周期は保たれるが、振幅は距離とともに急速に減る。
:::

### 回折場を開口の縁で積分する

平面波が垂直に入射し、境界値を開口 $A$ で $U_0$、外部でゼロとする。開口面で $\boldsymbol q=\boldsymbol\rho'-\boldsymbol\rho$、$R=\sqrt{z^2+|\boldsymbol q|^2}$ とおく。$\boldsymbol q\ne0$ では、面内のベクトル場 $\boldsymbol F=(\boldsymbol q/|\boldsymbol q|^2)e^{ikR}/R$ の発散は $R^{-1}(\dd/\dd R)(e^{ikR}/R)$。実際、面内極座標の半径を $q$ とすると発散は $q^{-1}\dd(qF_q)/\dd q$ で、$\dd R/\dd q=q/R$ を使えばよい。したがって式 :eqref[eq-optics-rayleigh-sommerfeld] の被積分関数は $-U_0z/(2\pi)$ をこの発散に掛けたものになる。

観測点の射影 $\boldsymbol\rho$ が開口内にあれば、その点を囲む小円を除いて発散定理を適用する。穴の外向き法線は $-\boldsymbol q/q$ なので、小円からは極限で $-2\pi e^{ikz}/z$ が得られる。射影が開口外にあれば小円は要らない。縁の外向き単位法線を $\hat{\boldsymbol n}_A$ とすると、射影が縁上にない点で
$$
\frac{U}{U_0}=\boldsymbol 1_A(\boldsymbol\rho)e^{ikz}
-\frac z{2\pi}\oint_{\partial A}\frac{e^{ikR}}{R}
 \frac{\boldsymbol q\cdot\hat{\boldsymbol n}_A}{|\boldsymbol q|^2}\dd s .
$$
第一項は光線が開口を通るときの幾何光学場、第二項は縁の寄与である。この形は第一 Rayleigh–Sommerfeld 積分と垂直平面波に対する :term[境界回折波]{en="boundary diffraction wave" yomi="きょうかいかいせつは"} の表示である。円孔の軸上では縁に沿って $R$ が一定であり、$\dd s=a\dd\varphi$ と $\boldsymbol q\cdot\hat{\boldsymbol n}_A=a$ を入れると、縁の項は $-(z/R)e^{ikR}$ になる。

## Fresnel 回折と Fresnel 帯

### 二次位相の伝播
$q=|\boldsymbol\rho-\boldsymbol\rho'|$ とし、$q/z\ll1$ のとき $R=z+q^2/(2z)-q^4/(8z^3)+\cdots$。振幅側を $1/z$ に、位相を二次まで近似するには $kq^4/(8z^3)\ll1$ も要る。$kz\gg1$ なら式 :eqref[eq-optics-rayleigh-sommerfeld] は
$$
U(\boldsymbol\rho,z)\simeq\frac{e^{ikz}}{i\lambda z}
\int_A U(\boldsymbol\rho',0)
e^{i k|\boldsymbol\rho-\boldsymbol\rho'|^2/(2z)}\dd^2\rho' .
\label{eq-optics-fresnel-diffraction}
$$
これが :term[Fresnel 回折]{en="Fresnel diffraction" yomi="ふれねるかいせき"} である。半径 $a$ に対して $N_F=a^2/(\lambda z)$ を Fresnel 数と呼ぶ。軸上では開口の端と中心の二次位相差は $\pi N_F$ である。幾何学的な小角条件だけでは、この位相を落とせない。

### 半周期ずつ分ける
軸上観測で $k r^2/(2z)=m\pi$ となる半径 $r_m=\sqrt{m\lambda z}$ を境界にする。同じ幅 $\lambda z$ の $r^2$ 区間なので各帯の面積は $\pi\lambda z$ である。平面波の照明で傾斜因子の変化を捨てれば、隣接帯の位相は $\pi$ 違い、積分は交互符号になる。例えば円孔の軸上では $2\pi r\dd r$ を積分して
$$
U(0,z)\simeq U_0 e^{ikz}[1-e^{i\pi N_F}] .
$$
これは有限開口の式であり、同じ近似で半径を無限大に送って通常の極限を取ることはできない。平面波の積分には収束因子を入れてから外すなどの扱いが必要である。

:::example[一帯だけを開く]
$N_F=1$ なら軸上の振幅は遮りのない平面波の約 2 倍、強度は 4 倍である。面全体のエネルギーが増えるのではなく、軸の近くへ再分配される。
:::
:::example[二帯を開く]
$N_F=2$ なら二帯の寄与が打ち消し合い、同近似の軸上強度はゼロとなる。遮蔽を減らしても固定した位置の強度は単調に増えない。
:::

### 隔帯を通すゾーンプレート
入射が平面波で、軸上の設計距離を $f$ とする。半径 $\sqrt{m\lambda f}$ で区切り、奇数番目の帯だけを通すと、隣り合う開帯間の位相周期は $2\pi$。設計点に向かう軸上振幅は各帯の同相の和になる。一つの開帯の積分は遮りのない平面波の 2 倍に等しいので、$M$ 個開けば理想 Fresnel 計算で振幅が $2M$ 倍になる。

透過関数を $t=r^2/(\lambda f)$ の周期 2 の矩形波とみなし、$e^{-i\pi j t}$ の係数を求めると $c_j=\frac12\int_0^1 e^{i\pi j t}\dd t$。$j$ が正の奇数なら $|c_j|=1/(\pi j)$ で、対応する因子 $e^{-i\pi j t}$ の二次位相は焦点距離 $f/j$ の薄レンズと等しい。従って多数の帯があり各次数の焦点を分離できる近似では、帯板全域への入射パワーのうち各次数が受け取る割合は $|c_j|^2$ である。一方、$f/j$ での軸上では一帯の位相幅が $j\pi$ なので、奇数 $j$ で軸上振幅の比はなお $2M$ となる。高次で集光パワーが小さくなることと、軸上照度そのものが小さくなることは同値ではない。短い焦点距離では同じ開口の焦点スポットが狭くなる。

:::example[第一帯の半径]
$\lambda=500\,\mathrm{nm},f=0.20\,\mathrm m$ の設計で第一帯境界は $\sqrt{\lambda f}\simeq316\,\mathrm{\mu m}$。第 4 境界はその二倍となる。
:::

### 点照明で Fresnel 帯の幅を変える
軸上の点光源から開口面までを $z_i>0$、開口から軸上の観測点までを $z_o>0$ とする。開口半径方向 $r$ で位相は $k(\sqrt{z_i^2+r^2}+\sqrt{z_o^2+r^2})$。$r$ の二次までで位相差は $kr^2/(2z_{\rm eff})$、$z_{\rm eff}^{-1}=z_i^{-1}+z_o^{-1}$。従って第 $m$ 帯の境界は $r_m=\sqrt{m\lambda z_{\rm eff}}$。ゾーンプレートの固定周期に対して $z_{\rm eff}$ を固定すれば、物距離を変えたときの像距離が薄レンズの結像式と同じ関係で動く。光源が開口に近いときは照明振幅もずれるので、帯の面積が一定というだけでは各帯の寄与は一致しない。

### 円孔の軸外の Fresnel 場

半径 $a$ の一様な円孔を平面波で照らす。観測半径を $\rho$ とし、開口の極角 $\varphi$ を $\boldsymbol\rho$ の方向から測ると、Fresnel 核の交差項は $e^{-ik\rho r\cos\varphi/z}$ になる。その角平均を $J_0(k\rho r/z)$ と書けば、軸外の積分は
$$
\frac{U(\rho,z)}{U_0}=
\frac{2\pi e^{ikz+ik\rho^2/(2z)}}{i\lambda z}
\int_0^a r e^{ikr^2/(2z)}J_0(k\rho r/z)\dd r .
$$
$\rho=0$ では $J_0(0)=1$ なので、$r^2$ を変数にすれば先の軸上式が得られる。逆に $ka^2/(2z)\ll1$ のときは被積分関数の二次位相を 1 に近似でき、そのとき残る積分が遠方の Airy パターンを与える。二次位相が無視できない場合は、この一次元積分に残して計算する。

### Fresnel 帯板の高次焦点と効率

近軸の帯板の座標を $t=r^2/(\lambda f)$ とする。半周期で交互に開く帯板は $t$ に周期 $2$ の透過率を持つ。$0<t<1$ を開く周期関数の Fourier 係数は $c_0=1/2$、$c_m=(1-e^{-im\pi})/(2im\pi)$ である。これは係数の積分を $[0,1]$ 上で行った結果である。Fresnel 伝播の軸上核は $e^{i\pi(f/z)t}$ を掛け、面積測度は $\dd t$ に比例する。従って $z=f/m$ で指数を打ち消す $c_{-m}$ が正の高次焦点を作る。多数の帯があるとき隣の次数が分離し、奇数次の理想パワー効率は $1/(m\pi)^2$。実際の縁の有限幅や帯数は焦点の広がりを生む。

:::example[主焦点の効率]
多数の帯を持つ振幅帯板で一次焦点を選ぶ。

入力は設計波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、一次焦点距離 $f=\frac{1}{5}\,\mathrm m$、最外側の帯番号 $N=100$、調べる正の次数 $m=1$ である。

設計した半周期帯の条件から半径を取る。

$$
\begin{aligned}
a&=\sqrt{N f \lambda},\\
a/(\mathrm m)&\simeq \sqrt{\left(100\right) \cdot \left(\frac{1}{5}\right) \cdot \left(5\times 10^{-7}\right)}\simeq 0.0031623.
\end{aligned}
$$

横二次位相を打ち消す距離を求める。

$$
\begin{aligned}
z_m&=\frac{f}{m},\\
z_m/(\mathrm m)&= \frac{\left(\frac{1}{5}\right)}{\left(1\right)}=\frac{1}{5}.
\end{aligned}
$$

Fourier 係数の二乗を用いる。

$$
\begin{aligned}
\eta_m&=\frac{\sin^{2}{\left(\frac{\pi m}{2} \right)}}{\pi^{2} m^{2}},\\
\eta_m&= \frac{\sin^{2}{\left(\frac{\pi \cdot \left(1\right)}{2} \right)}}{\pi^{2} \cdot \left(1\right)^{2}}\simeq 0.10132.
\end{aligned}
$$

開口外端の角を評価して、近軸の次数公式が使えるか調べる。

$$
\begin{aligned}
u&=\frac{a}{z_m},\\
u&\simeq \frac{\left(0.00316227766\right)}{\left(\frac{1}{5}\right)}\simeq 0.015811.
\end{aligned}
$$

一次へ向かうパワーは入射の約一割である。
:::

:::example[偶数次の位置]
同じ板の半分の焦点距離を調べる。

入力は設計波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、一次焦点距離 $f=\frac{1}{5}\,\mathrm m$、最外側の帯番号 $N=100$、調べる正の次数 $m=2$ である。

設計した半周期帯の条件から半径を取る。

$$
\begin{aligned}
a&=\sqrt{N f \lambda},\\
a/(\mathrm m)&\simeq \sqrt{\left(100\right) \cdot \left(\frac{1}{5}\right) \cdot \left(5\times 10^{-7}\right)}\simeq 0.0031623.
\end{aligned}
$$

横二次位相を打ち消す距離を求める。

$$
\begin{aligned}
z_m&=\frac{f}{m},\\
z_m/(\mathrm m)&= \frac{\left(\frac{1}{5}\right)}{\left(2\right)}=\frac{1}{10}.
\end{aligned}
$$

Fourier 係数の二乗を用いる。

$$
\begin{aligned}
\eta_m&=\frac{\sin^{2}{\left(\frac{\pi m}{2} \right)}}{\pi^{2} m^{2}},\\
\eta_m&= \frac{\sin^{2}{\left(\frac{\pi \cdot \left(2\right)}{2} \right)}}{\pi^{2} \cdot \left(2\right)^{2}}=0.
\end{aligned}
$$

開口外端の角を評価して、近軸の次数公式が使えるか調べる。

$$
\begin{aligned}
u&=\frac{a}{z_m},\\
u&\simeq \frac{\left(0.00316227766\right)}{\left(\frac{1}{10}\right)}\simeq 0.031623.
\end{aligned}
$$

理想の周期関数には二次の係数がない。
:::

:::example[三次の短い焦点]
同じ板の三次焦点を使う。

入力は設計波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、一次焦点距離 $f=\frac{1}{5}\,\mathrm m$、最外側の帯番号 $N=100$、調べる正の次数 $m=3$ である。

設計した半周期帯の条件から半径を取る。

$$
\begin{aligned}
a&=\sqrt{N f \lambda},\\
a/(\mathrm m)&\simeq \sqrt{\left(100\right) \cdot \left(\frac{1}{5}\right) \cdot \left(5\times 10^{-7}\right)}\simeq 0.0031623.
\end{aligned}
$$

横二次位相を打ち消す距離を求める。

$$
\begin{aligned}
z_m&=\frac{f}{m},\\
z_m/(\mathrm m)&= \frac{\left(\frac{1}{5}\right)}{\left(3\right)}=\frac{1}{15}.
\end{aligned}
$$

Fourier 係数の二乗を用いる。

$$
\begin{aligned}
\eta_m&=\frac{\sin^{2}{\left(\frac{\pi m}{2} \right)}}{\pi^{2} m^{2}},\\
\eta_m&= \frac{\sin^{2}{\left(\frac{\pi \cdot \left(3\right)}{2} \right)}}{\pi^{2} \cdot \left(3\right)^{2}}\simeq 0.011258.
\end{aligned}
$$

開口外端の角を評価して、近軸の次数公式が使えるか調べる。

$$
\begin{aligned}
u&=\frac{a}{z_m},\\
u&\simeq \frac{\left(0.00316227766\right)}{\left(\frac{1}{15}\right)}\simeq 0.047434.
\end{aligned}
$$

角は一次の三倍で、近軸模型でのパワー効率は一次の九分の一になる。ただし球面波の位相を二次で打ち切ると外端の誤差は $ka^4/(8z_3^3)\simeq0.53\,\mathrm{rad}$。従ってこの条件で実際の焦点を精密に求めるには、角だけで判定せず高次の位相を含めて計算する。
:::

## Fresnel 積分と Cornu 螺旋

### 直線の縁での積分
一次元の開口 $a<x'<b$ で、$y$ 方向に一様な平面波を二次位相の式で伝播する。変数 $v=\sqrt{2/(\lambda z)}(x'-x)$ を使うと位相は $\pi v^2/2$。:term[Fresnel 積分]{en="Fresnel integrals" yomi="ふれねるせきぶん"} は
$$
C(v)=\int_0^v\cos(\pi s^2/2)\dd s,\qquad S(v)=\int_0^v\sin(\pi s^2/2)\dd s .
$$
全実軸の $\int e^{i\pi v^2/2}\dd v=1+i$ は $\int e^{-(\epsilon-i\pi/2)v^2}\dd v$ を評価し $\epsilon\downarrow0$ として得る。したがって遮りのない場に対する振幅比は $[(C+iS)(v_b)-(C+iS)(v_a)]/(1+i)$ である。:term[Cornu 螺旋]{en="Cornu spiral" yomi="こるにゅらせん"} は $(C(v),S(v))$ の軌跡で、$v_a,v_b$ を結ぶ有向弦が分子を表す。弧長を $v$ で測れるのは、導関数 $(\cos,\sin)$ の長さが 1 だからである。

:::example[半平面の境界]
縁が $x'=0$ で正の半面が開いているとき、幾何光学的な影の境界 $x=0$ では積分区間は $[0,\infty)$。被積分関数は偶関数なので振幅は $1/2$、照度は $1/4$ である。
:::
:::example[対称なスリット]
$x=0$、$a=-b$ なら $C,S$ は奇関数なので弦は $2(C(v_b),S(v_b))$ である。瞳を狭くすると $C(v_b)\simeq v_b$、$S(v_b)=O(v_b^3)$ となる。
:::

:::example[影の深さを無次元化する]
波長 $500\,\mathrm{nm}$、距離 $0.50\,\mathrm m$ では Fresnel 積分の引数が 1 になる横距離は $\sqrt{\lambda z/2}\simeq0.354\,\mathrm{mm}$ である。同じ直線の縁でも遠くで見ると、この横の尺度は $\sqrt z$ に従って拡大する。
:::

### 二つの縁を重ねてスリットを考える
幅 $b-a$ のスリットは、$x'>a$ の半平面と $x'>b$ の半平面の差である。各端で位相を保った回折場を引くと、Cornu 螺旋の弦が端点間の差になる。強度を直接引いても一般に同じ答えにならない。逆に細い不透明帯では、無遮蔽の場からこの弦に対応する場を引く。Fresnel 帯と同程度の幅では二端のパターンが重なるので、独立な二つの片縁強度として足さない。

### 半平面の縁を基準に像の尺度を出す

開いている半平面を $x'>0$、観測座標を $x$ とする。Fresnel 積分で $u=\sqrt{2/(\lambda z)}(x'-x)$ と変数変換し、偶関数の被積分関数を用いると、遮蔽がないときの軸上振幅に対する比は $[(1/2+C(v))+i(1/2+S(v))]/(1+i)$、$v=\sqrt{2/(\lambda z)}x$ である。分母は全直線の積分で、Gauss 積分を微小の正の実減衰を入れて評価すると $1+i$。従って $v=0$ で強度は四分の一、開いた側で $v\to\infty$ なら一へ近づく。影の側を正に取る定義では $C,S$ の符号が逆になる。

:::example[幾何学的な境界上]
半平面の縁の真後ろで照度を測る。

入力は波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、伝播距離 $z=1\,\mathrm m$、縁の投影からの開口側座標 $x=0\,\mathrm m$ である。

縁がぼける幅の自然な尺度を求める。

$$
\begin{aligned}
d&=\frac{\sqrt{2} \sqrt{\lambda z}}{2},\\
d/(\mathrm m)&\simeq \frac{\sqrt{2} \cdot \sqrt{\left(5\times 10^{-7}\right) \cdot \left(1\right)}}{2}=\frac{1}{2000}.
\end{aligned}
$$

縁の投影を原点とする。

$$
\begin{aligned}
v&=\frac{x}{d},\\
v&= \frac{\left(0\right)}{\left(\frac{1}{2000}\right)}=0.
\end{aligned}
$$

実部は余弦積分の変化を含む。

$$
\begin{aligned}
A&=C\left(v\right) + \frac{1}{2},\\
A&= C\left(\left(0\right)\right) + \frac{1}{2}=\frac{1}{2}.
\end{aligned}
$$

虚部は正弦積分の変化を含む。

$$
\begin{aligned}
B&=S\left(v\right) + \frac{1}{2},\\
B&= S\left(\left(0\right)\right) + \frac{1}{2}=\frac{1}{2}.
\end{aligned}
$$

全直線の振幅の二乗で規格化する。

$$
\begin{aligned}
I/I_0&=\frac{A^{2}}{2} + \frac{B^{2}}{2},\\
I/I_0&= \frac{\left(\frac{1}{2}\right)^{2}}{2} + \frac{\left(\frac{1}{2}\right)^{2}}{2}=\frac{1}{4}.
\end{aligned}
$$

振幅の半分という結果を照度へ換えると四分の一になる。
:::

:::example[影側へ一尺度]
縁から遮蔽側へ移る。

入力は波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、伝播距離 $z=1\,\mathrm m$、縁の投影からの開口側座標 $x=- \frac{1}{2000}\,\mathrm m$ である。

縁がぼける幅の自然な尺度を求める。

$$
\begin{aligned}
d&=\frac{\sqrt{2} \sqrt{\lambda z}}{2},\\
d/(\mathrm m)&\simeq \frac{\sqrt{2} \cdot \sqrt{\left(5\times 10^{-7}\right) \cdot \left(1\right)}}{2}=\frac{1}{2000}.
\end{aligned}
$$

縁の投影を原点とする。

$$
\begin{aligned}
v&=\frac{x}{d},\\
v&= \frac{\left(- \frac{1}{2000}\right)}{\left(\frac{1}{2000}\right)}=-1.
\end{aligned}
$$

実部は余弦積分の変化を含む。

$$
\begin{aligned}
A&=C\left(v\right) + \frac{1}{2},\\
A&= C\left(\left(-1\right)\right) + \frac{1}{2}\simeq -0.27989.
\end{aligned}
$$

虚部は正弦積分の変化を含む。

$$
\begin{aligned}
B&=S\left(v\right) + \frac{1}{2},\\
B&= S\left(\left(-1\right)\right) + \frac{1}{2}\simeq 0.061741.
\end{aligned}
$$

全直線の振幅の二乗で規格化する。

$$
\begin{aligned}
I/I_0&=\frac{A^{2}}{2} + \frac{B^{2}}{2},\\
I/I_0&\simeq \frac{\left(-0.2798934004\right)^{2}}{2} + \frac{\left(0.06174085261\right)^{2}}{2}\simeq 0.041076.
\end{aligned}
$$

照度はゼロにはならず、影の側へ回り込む。
:::

:::example[開口側へ一尺度]
縁から照明側へ同じ距離だけ移る。

入力は波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、伝播距離 $z=1\,\mathrm m$、縁の投影からの開口側座標 $x=\frac{1}{2000}\,\mathrm m$ である。

縁がぼける幅の自然な尺度を求める。

$$
\begin{aligned}
d&=\frac{\sqrt{2} \sqrt{\lambda z}}{2},\\
d/(\mathrm m)&\simeq \frac{\sqrt{2} \cdot \sqrt{\left(5\times 10^{-7}\right) \cdot \left(1\right)}}{2}=\frac{1}{2000}.
\end{aligned}
$$

縁の投影を原点とする。

$$
\begin{aligned}
v&=\frac{x}{d},\\
v&= \frac{\left(\frac{1}{2000}\right)}{\left(\frac{1}{2000}\right)}=1.
\end{aligned}
$$

実部は余弦積分の変化を含む。

$$
\begin{aligned}
A&=C\left(v\right) + \frac{1}{2},\\
A&= C\left(\left(1\right)\right) + \frac{1}{2}\simeq 1.2799.
\end{aligned}
$$

虚部は正弦積分の変化を含む。

$$
\begin{aligned}
B&=S\left(v\right) + \frac{1}{2},\\
B&= S\left(\left(1\right)\right) + \frac{1}{2}\simeq 0.93826.
\end{aligned}
$$

全直線の振幅の二乗で規格化する。

$$
\begin{aligned}
I/I_0&=\frac{A^{2}}{2} + \frac{B^{2}}{2},\\
I/I_0&\simeq \frac{\left(1.2798934\right)^{2}}{2} + \frac{\left(0.9382591474\right)^{2}}{2}\simeq 1.2592.
\end{aligned}
$$

照度は一を越え得る。幾何学的な開閉の割合ではこの明るさを説明できない。
:::

## Fraunhofer 回折（矩形開口・円形開口と Airy パターン）

### 二次位相を捨てられる距離
:eqref[eq-optics-fresnel-diffraction] の指数を $|\boldsymbol\rho|^2-2\boldsymbol\rho\cdot\boldsymbol\rho'+|\boldsymbol\rho'|^2$ に分ける。開口内で $k|\boldsymbol\rho'|^2/(2z)\ll1$ なら最後の項を無視できる。これが :term[Fraunhofer 回折]{en="Fraunhofer diffraction" yomi="ふらうんほーふぁーかいせき"} であり、正面近くでは開口場の Fourier 変換を $\boldsymbol q=k\boldsymbol\rho/z$ で評価する。一般の遠方角では、開口の中心からの観測距離を $R_0=|\vb r|$、観測方向を $\hat{\vb r}=\vb r/R_0$ として $|\vb r-\vb r'|=R_0-\hat{\vb r}\cdot\vb r'+O(|\vb r'|^2/R_0)$ と展開する。従って $k|\vb r'|^2/R_0\ll1$ では $q_x=k\hat r_x,\ q_y=k\hat r_y$ を使える。この場合は核の傾斜因子も角度に応じて残る。以下の強度比は、近軸でその変化を無視できる場合、またはその共通因子を除いた回折因子について述べる。

開口を一様振幅、同位相で照明する。幅 $a,b$ の中心付き矩形の積分は独立な一次元積分に分かれ、
$$
\frac{U(q_x,q_y)}{U(0,0)}=\operatorname{sinc}(q_x a/2)\operatorname{sinc}(q_y b/2),\qquad
\operatorname{sinc}u=\sin(u)/u,\quad \operatorname{sinc}0=1 .
$$
ここでは観測面での共通の曲率位相を除いた。円孔の半径を $a$ として極変数を取ると、角度積分 $2\pi J_0(qr)=\int_0^{2\pi}e^{-iqr\cos\varphi}\dd\varphi$ が現れる。$J_0$ のこの式を指数展開すると $J_0(u)=\sum_{j\ge0}(-1)^j(u/2)^{2j}/(j!)^2$ になる。項別積分で $\int_0^a J_0(qr)r\dd r=aJ_1(qa)/q$ を得る。従って振幅は中心を 1 として $2J_1(qa)/(qa)$、二乗が :term[Airy パターン]{en="Airy pattern" yomi="えありーぱたーん"} である。$J_1$ の第一正根 $3.8317$ から、最初の暗環の角は小角で $1.22\lambda/D$、$D=2a$ となる。

:::example[矩形開口]
$a=0.10\,\mathrm{mm}$、$\lambda=500\,\mathrm{nm}$ では $x$ 方向第一暗線の角は約 $0.005\,\mathrm{rad}$。幅を半分にすると暗線の角は二倍となる。
:::
:::example[円孔と像の半径]
$D=5\,\mathrm{mm}$、$\lambda=500\,\mathrm{nm}$ の円孔の暗環角は $1.22\times10^{-4}\,\mathrm{rad}$。焦点距離 $100\,\mathrm{mm}$ の系なら半径は約 $12.2\,\mathrm{\mu m}$ である。
:::

### 円形点像の中心に入るエネルギー
規格化された半径 $x=qa$ で Airy 強度を $[2J_1(x)/x]^2$ とする。半径 $X$ までのエネルギーの割合は $1-J_0(X)^2-J_1(X)^2$ である。これを示すには $J_0'=-J_1$、$J_1'=J_0-J_1/x$ を用いる。前の級数から項別微分で両式が得られ、$[1-J_0^2-J_1^2]'=2J_1^2/x$ となる。これは強度を面積測度 $x\dd x$ で積分したものに比例する。$J_0(0)=1,J_1(0)=0$、十分遠くで両関数がゼロへ収束することから、上の積分の全域での値が 1 になる。遠方の収束は角積分を、位相微分が消える方向の小近傍とそれ以外に分ければ示せる。近傍の寄与は幅を小さくして抑え、それ以外は部分積分で $O(X^{-1})$ となる。

:::example[第一暗環の内側]
$X=3.8317$ で $J_1(X)\simeq0$、$J_0(X)\simeq-0.4028$。従って中心の明領域には全エネルギーの約 $83.8$ パーセントが入る。
:::

### 二点の分離と第一暗環
独立な同強度の二点を撮れば二点像強度の和になる。Rayleigh 基準は片方の中心と他方の第一暗環を重ねる間隔である。中点は両中心から同じ距離なので、$x$ を Airy の規格化半径として中点強度は $2[2J_1(x_1/2)/(x_1/2)]^2$、$x_1\simeq3.8317$。片方の像中心では他方の寄与がゼロなので、強度は単独ピークに等しい。中点はその約 $0.735$ であり、二つの像の間にくぼみが残る。完全なコヒーレント結像では二点が同位相か逆位相かで中点の交差項が変わるため、この割合をそのまま使えない。

### 矩形瞳を斜めに見るときの二方向の幅

矩形の幅を $a,b$ とする。Fraunhofer 方向の単位ベクトルの横成分を $s_x,s_y$ とすれば、積分は二つの一次元積分の積で、中心で規格化した振幅は $\operatorname{sinc}(\pi a s_x/\lambda)\operatorname{sinc}(\pi b s_y/\lambda)$。第一の零線は $s_x=\pm\lambda/a$ および $s_y=\pm\lambda/b$ である。スクリーンで大きな角を考えるなら、$x/z=s_x/\sqrt{1-s_x^2-s_y^2}$ として幾何学で換算する。観測面が傾くと素片の投影や核の傾斜因子も変わるので、ここでは開口の Fourier 因子だけを比較する。

:::example[横の暗線に乗る点]
横幅の最初の暗線と縦方向の斜め位置を組み合わせる。

入力は横幅 $a\simeq 1\times 10^{-4}\,\mathrm m$、縦幅 $b=\frac{1}{5000}\,\mathrm m$、波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、横方向余弦 $s_x=\frac{1}{200}$、縦方向余弦 $s_y=\frac{1}{1000}$ である。

横幅の両端の位相差の半分を取る。

$$
\begin{aligned}
u&=\frac{\pi a s_x}{\lambda},\\
u&\simeq \frac{\pi \cdot \left(1\times 10^{-4}\right) \cdot \left(\frac{1}{200}\right)}{\left(5\times 10^{-7}\right)}\simeq 3.1416.
\end{aligned}
$$

同じ操作を縦幅で行う。

$$
\begin{aligned}
v&=\frac{\pi b s_y}{\lambda},\\
v&\simeq \frac{\pi \cdot \left(\frac{1}{5000}\right) \cdot \left(\frac{1}{1000}\right)}{\left(5\times 10^{-7}\right)}\simeq 1.2566.
\end{aligned}
$$

独立な二方向の積分を掛ける。

$$
\begin{aligned}
A&=\operatorname{sinc}{\left(u \right)} \operatorname{sinc}{\left(v \right)},\\
A&\simeq \operatorname{sinc}{\left(\left(3.141592654\right) \right)} \cdot \operatorname{sinc}{\left(\left(1.256637061\right) \right)}=0.
\end{aligned}
$$

複素振幅を作ってから二乗する。

$$
\begin{aligned}
I/I_0&=A^{2},\\
I/I_0&= \left(0\right)^{2}=0.
\end{aligned}
$$

開口が細い方向で零線がどれだけ離れるか求める。

$$
\begin{aligned}
s_{y,1}&=\frac{\lambda}{b},\\
s_{y,1}&\simeq \frac{\left(5\times 10^{-7}\right)}{\left(\frac{1}{5000}\right)}=\frac{1}{400}.
\end{aligned}
$$

縦成分がゼロでなくても横の積分だけで完全に消える。
:::

:::example[二軸とも明るい領域]
同じ開口で横も縦も第一暗線より内側へ取る。

入力は横幅 $a\simeq 1\times 10^{-4}\,\mathrm m$、縦幅 $b=\frac{1}{5000}\,\mathrm m$、波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、横方向余弦 $s_x=\frac{1}{1000}$、縦方向余弦 $s_y=\frac{1}{1000}$ である。

横幅の両端の位相差の半分を取る。

$$
\begin{aligned}
u&=\frac{\pi a s_x}{\lambda},\\
u&\simeq \frac{\pi \cdot \left(1\times 10^{-4}\right) \cdot \left(\frac{1}{1000}\right)}{\left(5\times 10^{-7}\right)}\simeq 0.62832.
\end{aligned}
$$

同じ操作を縦幅で行う。

$$
\begin{aligned}
v&=\frac{\pi b s_y}{\lambda},\\
v&\simeq \frac{\pi \cdot \left(\frac{1}{5000}\right) \cdot \left(\frac{1}{1000}\right)}{\left(5\times 10^{-7}\right)}\simeq 1.2566.
\end{aligned}
$$

独立な二方向の積分を掛ける。

$$
\begin{aligned}
A&=\operatorname{sinc}{\left(u \right)} \operatorname{sinc}{\left(v \right)},\\
A&\simeq \operatorname{sinc}{\left(\left(0.6283185307\right) \right)} \cdot \operatorname{sinc}{\left(\left(1.256637061\right) \right)}\simeq 0.708.
\end{aligned}
$$

複素振幅を作ってから二乗する。

$$
\begin{aligned}
I/I_0&=A^{2},\\
I/I_0&\simeq \left(0.7080032943\right)^{2}\simeq 0.50127.
\end{aligned}
$$

開口が細い方向で零線がどれだけ離れるか求める。

$$
\begin{aligned}
s_{y,1}&=\frac{\lambda}{b},\\
s_{y,1}&\simeq \frac{\left(5\times 10^{-7}\right)}{\left(\frac{1}{5000}\right)}=\frac{1}{400}.
\end{aligned}
$$

両軸の減少を掛けるため、どちらか一軸だけの値より暗い。
:::

:::example[縦を細くした開口]
横幅は変えずに縦幅だけを縮める。

入力は横幅 $a\simeq 1\times 10^{-4}\,\mathrm m$、縦幅 $b\simeq 5\times 10^{-5}\,\mathrm m$、波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、横方向余弦 $s_x=\frac{1}{1000}$、縦方向余弦 $s_y=\frac{1}{1000}$ である。

横幅の両端の位相差の半分を取る。

$$
\begin{aligned}
u&=\frac{\pi a s_x}{\lambda},\\
u&\simeq \frac{\pi \cdot \left(1\times 10^{-4}\right) \cdot \left(\frac{1}{1000}\right)}{\left(5\times 10^{-7}\right)}\simeq 0.62832.
\end{aligned}
$$

同じ操作を縦幅で行う。

$$
\begin{aligned}
v&=\frac{\pi b s_y}{\lambda},\\
v&\simeq \frac{\pi \cdot \left(5\times 10^{-5}\right) \cdot \left(\frac{1}{1000}\right)}{\left(5\times 10^{-7}\right)}\simeq 0.31416.
\end{aligned}
$$

独立な二方向の積分を掛ける。

$$
\begin{aligned}
A&=\operatorname{sinc}{\left(u \right)} \operatorname{sinc}{\left(v \right)},\\
A&\simeq \operatorname{sinc}{\left(\left(0.6283185307\right) \right)} \cdot \operatorname{sinc}{\left(\left(0.3141592654\right) \right)}\simeq 0.92018.
\end{aligned}
$$

複素振幅を作ってから二乗する。

$$
\begin{aligned}
I/I_0&=A^{2},\\
I/I_0&\simeq \left(0.9201768613\right)^{2}\simeq 0.84673.
\end{aligned}
$$

開口が細い方向で零線がどれだけ離れるか求める。

$$
\begin{aligned}
s_{y,1}&=\frac{\lambda}{b},\\
s_{y,1}&\simeq \frac{\left(5\times 10^{-7}\right)}{\left(5\times 10^{-5}\right)}=\frac{1}{100}.
\end{aligned}
$$

縦方向での回折幅が広がり、同じ角では中心に近い強度になる。
:::

## 回折格子と分光器の分解能

### 周期と有限本数の効果
間隔 $d$ の同じ開口が $N\ge2$ 本あり、正面から同じ位相で照明するとする。角 $\theta$ の遠方で隣接する位相差を $\delta=kd\sin\theta$ とすると、単一開口の振幅に $\sum_{j=0}^{N-1}e^{-ij\delta}$ が掛かる。幾何級数を整理すれば強度因子は $\sin^2(N\delta/2)/\sin^2(\delta/2)$。この有限列因子の峰の条件は
$$
d\sin\theta_m=m\lambda ,
\label{eq-optics-grating}
$$
その周囲での第一零点は $|\delta-2\pi m|=2\pi/N$ である。$m\ne0$ で、単一開口の包絡線が峰の周りで緩やかに変わり、$\cos\theta_m\ne0$ とする。$d\cos\theta_m\dd\theta=m\dd\lambda$ だから、二波長の峰の間隔が片方の第一零点までの角幅に等しい条件は $\Delta\lambda\simeq\lambda/(|m|N)$。この :term[分解能]{en="resolving power" yomi="ぶんかいのう"} $\lambda/\Delta\lambda\simeq|m|N$ は均一に照明した理想格子の Rayleigh 基準である。実機では入射スリットの幅、結像収差や検出画素も波長幅に換算する。

:::example[スペクトル線の幅]
照明された本数 $N=1000$、$m=2$ で $600\,\mathrm{nm}$ の近傍を測る。理想条件で $\Delta\lambda=0.30\,\mathrm{nm}$ である。
:::
:::example[次数の重なり]
$600\,\mathrm{nm}$ の 1 次と $300\,\mathrm{nm}$ の 2 次の $m\lambda$ は同じなので角も同じになる。回折角単独から二波長を識別できない。
:::

### 斜め入射と開口率
接線方向の進行を同じ座標で測り、入射角と観測角を符号付きで $\theta_i,\theta_o$ とする。隣接開口で照明の位相が $kd\sin\theta_i$ だけ進むため、遠方では $d(\sin\theta_o-\sin\theta_i)=m\lambda$ になる。矩形の各スリットの幅を $a$ とすると、正面照明で各次数 $m$ の振幅には $\operatorname{sinc}(\pi m a/d)$ が掛かる。周期で許される次数でも、この因子の零点に重なると消える。

:::example[欠級]
$a/d=1/2$ なら $m=2,4,\dots$ の非零偶数次数で $\sin(\pi m/2)=0$。これらは周期の条件を満たすが、各スリット内の積分が打ち消し合って現れない。
:::

### 分光器の位置分散
レンズ焦点面で角位置を $x=F\tan\theta_m$ と表すモデルでは $\dd x/\dd\lambda=Fm/(d\cos^3\theta_m)$ になる。これは :eqref[eq-optics-grating] の角微分に $\dd x/\dd\theta=F/\cos^2\theta$ を掛けた結果である。小角では $Fm/d$ へ近づく。画素幅を $p$ とすると、一画素に相当する波長幅は $p/|\dd x/\dd\lambda|$。格子の第一零点がこれより細かい場合でも、そのまま検出側で分離できるとは限らない。線像をいくつの画素で標本化するかを別に評価する。

### 格子で線幅と次数間隔を別に数える

周期 $d$、照らされた溝数 $N$ の格子では隣の溝の位相差を $\psi=2\pi d\sin\theta/\lambda$ とする。和の因子 $\sum_{j=0}^{N-1}e^{ij\psi}$ は等比級数である。極大 $2\pi m$ の隣の零点は分子が初めてゼロになる $\psi=2\pi m\pm2\pi/N$。従って角の片側の幅は小幅で $\lambda/(Nd\cos\theta_m)$、一方格子式を微分した角分散の大きさは $|m|/(d\cos\theta_m)$。非零次数で比を取ると分離できる最小波長差の目安 $\lambda/(|m|N)$ を得る。違う次数どうしが同じ角に出る問題はこの線幅では解けず、次数を選ぶフィルターが要る。

:::example[一次スペクトル]
一次の細線を分ける目安を数える。

入力は溝の周期 $d\simeq 2\times 10^{-6}\,\mathrm m$、照明される溝数 $N=5000$、波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、正の次数 $m=1$ である。

垂直入射の格子条件を使う。

$$
\begin{aligned}
\sin\theta_m&=\frac{\lambda m}{d},\\
\sin\theta_m&\simeq \frac{\left(5\times 10^{-7}\right) \cdot \left(1\right)}{\left(2\times 10^{-6}\right)}=\frac{1}{4}.
\end{aligned}
$$

波長微分を角微分に換える。

$$
\begin{aligned}
\dd\theta_m/\dd\lambda&=\frac{m}{d \sqrt{1-(\sin\theta_m)^2}},\\
\dd\theta_m/\dd\lambda/(\mathrm{rad/m})&\simeq \frac{\left(1\right)}{\left(2\times 10^{-6}\right) \cdot \sqrt{1 - \left(\frac{1}{4}\right)^{2}}}\simeq 5.164\times 10^{5}.
\end{aligned}
$$

隣の溝の位相差を $2\pi/N$ 動かす。

$$
\begin{aligned}
\Delta\theta&=\frac{\lambda}{N d \sqrt{1-(\sin\theta_m)^2}},\\
\Delta\theta/(\mathrm{rad})&\simeq \frac{\left(5\times 10^{-7}\right)}{\left(5000\right) \cdot \left(2\times 10^{-6}\right) \cdot \sqrt{1 - \left(\frac{1}{4}\right)^{2}}}\simeq 5.164\times 10^{-5}.
\end{aligned}
$$

同じ角幅を作る波長変化を求める。

$$
\begin{aligned}
\Delta\lambda&=\frac{\Delta\theta}{\dd\theta_m/\dd\lambda},\\
\Delta\lambda/(\mathrm m)&\simeq \frac{\left(5.163977795\times 10^{-5}\right)}{\left(5.163977795\times 10^{5}\right)}\simeq 1\times 10^{-10}.
\end{aligned}
$$

角度の同じ別次数は $m\lambda$ を維持する。

$$
\begin{aligned}
\lambda'&=\frac{\lambda m}{m + 1},\\
\lambda'/(\mathrm m)&\simeq \frac{\left(5\times 10^{-7}\right) \cdot \left(1\right)}{\left(1\right) + 1}\simeq 2.5\times 10^{-7}.
\end{aligned}
$$

線幅と、別次数の波長が重なる問題は別々に確認する。
:::

:::example[同じ格子の二次]
照明される溝数を変えずに次数を上げる。

入力は溝の周期 $d\simeq 2\times 10^{-6}\,\mathrm m$、照明される溝数 $N=5000$、波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、正の次数 $m=2$ である。

垂直入射の格子条件を使う。

$$
\begin{aligned}
\sin\theta_m&=\frac{\lambda m}{d},\\
\sin\theta_m&\simeq \frac{\left(5\times 10^{-7}\right) \cdot \left(2\right)}{\left(2\times 10^{-6}\right)}=\frac{1}{2}.
\end{aligned}
$$

波長微分を角微分に換える。

$$
\begin{aligned}
\dd\theta_m/\dd\lambda&=\frac{m}{d \sqrt{1-(\sin\theta_m)^2}},\\
\dd\theta_m/\dd\lambda/(\mathrm{rad/m})&\simeq \frac{\left(2\right)}{\left(2\times 10^{-6}\right) \cdot \sqrt{1 - \left(\frac{1}{2}\right)^{2}}}\simeq 1.1547\times 10^{6}.
\end{aligned}
$$

隣の溝の位相差を $2\pi/N$ 動かす。

$$
\begin{aligned}
\Delta\theta&=\frac{\lambda}{N d \sqrt{1-(\sin\theta_m)^2}},\\
\Delta\theta/(\mathrm{rad})&\simeq \frac{\left(5\times 10^{-7}\right)}{\left(5000\right) \cdot \left(2\times 10^{-6}\right) \cdot \sqrt{1 - \left(\frac{1}{2}\right)^{2}}}\simeq 5.7735\times 10^{-5}.
\end{aligned}
$$

同じ角幅を作る波長変化を求める。

$$
\begin{aligned}
\Delta\lambda&=\frac{\Delta\theta}{\dd\theta_m/\dd\lambda},\\
\Delta\lambda/(\mathrm m)&\simeq \frac{\left(5.773502692\times 10^{-5}\right)}{\left(1.154700538\times 10^{6}\right)}\simeq 5\times 10^{-11}.
\end{aligned}
$$

角度の同じ別次数は $m\lambda$ を維持する。

$$
\begin{aligned}
\lambda'&=\frac{\lambda m}{m + 1},\\
\lambda'/(\mathrm m)&\simeq \frac{\left(5\times 10^{-7}\right) \cdot \left(2\right)}{\left(2\right) + 1}\simeq 3.3333\times 10^{-7}.
\end{aligned}
$$

波長分解幅は半分になる。重なる次次数の波長は元の波長へ近づく。
:::

:::example[照明幅が狭い二次]
二次を使うまま照明する溝数を半分にする。

入力は溝の周期 $d\simeq 2\times 10^{-6}\,\mathrm m$、照明される溝数 $N=2500$、波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、正の次数 $m=2$ である。

垂直入射の格子条件を使う。

$$
\begin{aligned}
\sin\theta_m&=\frac{\lambda m}{d},\\
\sin\theta_m&\simeq \frac{\left(5\times 10^{-7}\right) \cdot \left(2\right)}{\left(2\times 10^{-6}\right)}=\frac{1}{2}.
\end{aligned}
$$

波長微分を角微分に換える。

$$
\begin{aligned}
\dd\theta_m/\dd\lambda&=\frac{m}{d \sqrt{1-(\sin\theta_m)^2}},\\
\dd\theta_m/\dd\lambda/(\mathrm{rad/m})&\simeq \frac{\left(2\right)}{\left(2\times 10^{-6}\right) \cdot \sqrt{1 - \left(\frac{1}{2}\right)^{2}}}\simeq 1.1547\times 10^{6}.
\end{aligned}
$$

隣の溝の位相差を $2\pi/N$ 動かす。

$$
\begin{aligned}
\Delta\theta&=\frac{\lambda}{N d \sqrt{1-(\sin\theta_m)^2}},\\
\Delta\theta/(\mathrm{rad})&\simeq \frac{\left(5\times 10^{-7}\right)}{\left(2500\right) \cdot \left(2\times 10^{-6}\right) \cdot \sqrt{1 - \left(\frac{1}{2}\right)^{2}}}\simeq 1.1547\times 10^{-4}.
\end{aligned}
$$

同じ角幅を作る波長変化を求める。

$$
\begin{aligned}
\Delta\lambda&=\frac{\Delta\theta}{\dd\theta_m/\dd\lambda},\\
\Delta\lambda/(\mathrm m)&\simeq \frac{\left(1.154700538\times 10^{-4}\right)}{\left(1.154700538\times 10^{6}\right)}\simeq 1\times 10^{-10}.
\end{aligned}
$$

角度の同じ別次数は $m\lambda$ を維持する。

$$
\begin{aligned}
\lambda'&=\frac{\lambda m}{m + 1},\\
\lambda'/(\mathrm m)&\simeq \frac{\left(5\times 10^{-7}\right) \cdot \left(2\right)}{\left(2\right) + 1}\simeq 3.3333\times 10^{-7}.
\end{aligned}
$$

分散自体は変わらないが、回折線幅は広がり、一次の例と同程度の波長分解幅になる。
:::

## Babinet の原理

### 相補開口と遮りのない波
スカラー伝播演算子を $\mathcal P$ とする。開口の透過関数 $A$ と相補開口 $1-A$ で同じ照明 $U_0$ を用いると、線形性により
$$
\mathcal P[AU_0]+\mathcal P[(1-A)U_0]=\mathcal P[U_0].
$$
これが :term[Babinet の原理]{en="Babinet's principle" yomi="ばびねのげんり"} のスカラー形である。平面波の遠方では右辺が直進方向に集中するので、直進成分を外した位置では相補開口の場は符号が逆で、強度は等しい。ベクトル電磁波と実在の導体にこの式を使う場合、偏光を含む境界問題を別に検討する必要がある。

:::example[細い線とスリット]
直径に相当する幅が $a$ のまっすぐな不透明線を平面波で照らす。直進成分を除けば幅 $a$ のスリットと同じ $\operatorname{sinc}^2$ の強度分布を持つ。最初の暗線から $a$ を測れる。
:::
:::example[円板の影の中央]
Fresnel 近似で同半径の円孔の軸上場が $U_0e^{ikz}(1-e^{i\pi N_F})$ なら、円板だけを入れた場は $U_0e^{ikz}e^{i\pi N_F}$ である。照度は軸上で遮りがない場合と同じになり、幾何学的な影の中心に明点が出る。
:::

:::example[開口と遮蔽の強度は一般には相補でない]
無遮蔽振幅を 1 とし、ある点で開口の振幅が $1+i$ とする。相補な遮蔽物の振幅は $-i$。強度の和は $2+1=3$ となり、無遮蔽強度 1 と一致しない。振幅を足してから二乗すると、交差項 $-2$ も含めて 1 が得られる。
:::

### 円板の軸上の影を複素振幅で評価する

円板と円孔に複素線形な伝播を適用すると、軸上で $U_{\rm disk}+U_{\rm hole}=U_0e^{ikz}$。円孔の厳密なスカラー Dirichlet 解 :eqref[eq-optics-rayleigh-sommerfeld] の軸上積分を引けば、円板側は $U_0(z/R)e^{ikR}$、$R=\sqrt{z^2+a^2}$ となる。従って影の中心強度は遮蔽前の $z^2/(z^2+a^2)$。二次近似で $a/z\ll1$ なら明るさはもととほぼ同じだが、中心で測る位相は開口を通した場合と別である。ベクトルの電磁境界を厳密に解く円板問題では、ここで指定したスカラー境界データをそのまま置けることを仮定しない。

:::example[遠い影の中央]
半径が距離に比べて十分小さい円板を用いる。

入力は円板の半径 $a=\frac{1}{1000}\,\mathrm m$、円板と観測点の距離 $z=1\,\mathrm m$、波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$ である。

全方位の縁が同じ距離を持つ。

$$
\begin{aligned}
R&=\sqrt{a^{2} + z^{2}},\\
R/(\mathrm m)&= \sqrt{\left(\frac{1}{1000}\right)^{2} + \left(1\right)^{2}}\simeq 1.
\end{aligned}
$$

近い二距離の差を有理化して求める。

$$
\begin{aligned}
\Delta R&=\frac{a^{2}}{R + z},\\
\Delta R/(\mathrm m)&\simeq \frac{\left(\frac{1}{1000}\right)^{2}}{\left(1.0000005\right) + \left(1\right)}\simeq 5\times 10^{-7}.
\end{aligned}
$$

遮蔽前の平面波と比較した位相を求める。

$$
\begin{aligned}
\Delta\phi&=\frac{2 \pi \Delta R}{\lambda},\\
\Delta\phi/(\mathrm{rad})&\simeq \frac{2 \cdot \pi \cdot \left(4.99999875\times 10^{-7}\right)}{\left(5\times 10^{-7}\right)}\simeq 6.2832.
\end{aligned}
$$

相補開口の振幅を引いて残った振幅の二乗を取る。

$$
\begin{aligned}
I/I_0&=\frac{z^{2}}{R^{2}},\\
I/I_0&\simeq \frac{\left(1\right)^{2}}{\left(1.0000005\right)^{2}}\simeq 1.
\end{aligned}
$$

強度がほぼ一でも位相は変化する。
:::

:::example[同じ板の近くの中心]
円板を固定し、観測距離だけを短くする。

入力は円板の半径 $a=\frac{1}{1000}\,\mathrm m$、円板と観測点の距離 $z=\frac{1}{250}\,\mathrm m$、波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$ である。

全方位の縁が同じ距離を持つ。

$$
\begin{aligned}
R&=\sqrt{a^{2} + z^{2}},\\
R/(\mathrm m)&= \sqrt{\left(\frac{1}{1000}\right)^{2} + \left(\frac{1}{250}\right)^{2}}\simeq 0.0041231.
\end{aligned}
$$

近い二距離の差を有理化して求める。

$$
\begin{aligned}
\Delta R&=\frac{a^{2}}{R + z},\\
\Delta R/(\mathrm m)&\simeq \frac{\left(\frac{1}{1000}\right)^{2}}{\left(0.004123105626\right) + \left(\frac{1}{250}\right)}\simeq 1.2311\times 10^{-4}.
\end{aligned}
$$

遮蔽前の平面波と比較した位相を求める。

$$
\begin{aligned}
\Delta\phi&=\frac{2 \pi \Delta R}{\lambda},\\
\Delta\phi/(\mathrm{rad})&\simeq \frac{2 \cdot \pi \cdot \left(1.231056256\times 10^{-4}\right)}{\left(5\times 10^{-7}\right)}\simeq 1547.
\end{aligned}
$$

相補開口の振幅を引いて残った振幅の二乗を取る。

$$
\begin{aligned}
I/I_0&=\frac{z^{2}}{R^{2}},\\
I/I_0&\simeq \frac{\left(\frac{1}{250}\right)^{2}}{\left(0.004123105626\right)^{2}}=\frac{16}{17}.
\end{aligned}
$$

近軸の中心強度一から少し外れる。
:::

:::example[半径と距離が同程度]
指定されたスカラー境界値の結果を比較する。

入力は円板の半径 $a=\frac{1}{1000}\,\mathrm m$、円板と観測点の距離 $z=\frac{1}{1000}\,\mathrm m$、波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$ である。

全方位の縁が同じ距離を持つ。

$$
\begin{aligned}
R&=\sqrt{a^{2} + z^{2}},\\
R/(\mathrm m)&= \sqrt{\left(\frac{1}{1000}\right)^{2} + \left(\frac{1}{1000}\right)^{2}}\simeq 0.0014142.
\end{aligned}
$$

近い二距離の差を有理化して求める。

$$
\begin{aligned}
\Delta R&=\frac{a^{2}}{R + z},\\
\Delta R/(\mathrm m)&\simeq \frac{\left(\frac{1}{1000}\right)^{2}}{\left(0.001414213562\right) + \left(\frac{1}{1000}\right)}\simeq 4.1421\times 10^{-4}.
\end{aligned}
$$

遮蔽前の平面波と比較した位相を求める。

$$
\begin{aligned}
\Delta\phi&=\frac{2 \pi \Delta R}{\lambda},\\
\Delta\phi/(\mathrm{rad})&\simeq \frac{2 \cdot \pi \cdot \left(4.142135624\times 10^{-4}\right)}{\left(5\times 10^{-7}\right)}\simeq 5205.2.
\end{aligned}
$$

相補開口の振幅を引いて残った振幅の二乗を取る。

$$
\begin{aligned}
I/I_0&=\frac{z^{2}}{R^{2}},\\
I/I_0&\simeq \frac{\left(\frac{1}{1000}\right)^{2}}{\left(0.001414213562\right)^{2}}=\frac{1}{2}.
\end{aligned}
$$

厳密核では二分の一である。Fresnel 近似はこの幾何では使えない。
:::

## 演習問題

::::exercise{level=basic}
$G=e^{ikR}/(4\pi R)$ の小球の法線勾配の流出を求めよ。

:::solution
勾配は $(ik-1/R)G$ である。球面積を掛けると $(ikR-1)e^{ikR}$ になり、$R\to0$ の極限は $-1$ である。
:::
::::

::::exercise{level=basic}
$\lambda=500\,\mathrm{nm},z=1\,\mathrm m,a=1\,\mathrm{mm}$ の Fresnel 数を求めよ。

:::solution
$N_F=a^2/(\lambda z)=10^{-6}/(5\times10^{-7})=2$ である。
:::
::::

::::exercise{level=basic}
一次元半平面の影の境界で、無遮蔽時に対する照度比を求めよ。

:::solution
Fresnel 核は観測点から測った横距離の偶関数なので、半分の積分は遮りのない場の半分になる。照度比は $1/4$ である。
:::
::::

::::exercise{level=basic}
$d=2\,\mathrm{\mu m},\lambda=500\,\mathrm{nm}$ の格子で $m=1$ の $\sin\theta$ を求めよ。

:::solution
式 :eqref[eq-optics-grating] より $\sin\theta=\lambda/d=1/4$ である。
:::
::::

::::exercise{level=advanced}
二次回折核を作る際に落とす最初の位相項を求めよ。

:::solution
$R=z\sqrt{1+q^2/z^2}=z+q^2/(2z)-q^4/(8z^3)+O(q^6/z^5)$。これに $k$ を掛けるので落とす最初の位相は $-kq^4/(8z^3)$ である。
:::
::::

::::exercise{level=advanced}
同じ格子を使い、周期と観測次数を固定したまま照明本数が半分になると Rayleigh 基準の波長差はどうなるか。

:::solution
近接した峰間の角分散は格子間隔と次数で決まり、本数にはよらない。第一零点までの角幅は $1/N$ に比例するので、分離に必要な波長差は二倍になる。
:::
::::

::::exercise{level=advanced}
$D=10\,\mathrm{mm},\lambda=500\,\mathrm{nm}$ の円孔で第一暗環角を求めよ。

:::solution
小角で $1.22\lambda/D=6.10\times10^{-5}\,\mathrm{rad}$ となる。
:::
::::

::::exercise{level=advanced}
Fresnel 近似の円孔の軸上振幅を使い、$N_F=3$ の場合の円板の照度比を求めよ。

:::solution
円孔の振幅比は $1-e^{i3\pi}=2$。相補な円板の振幅比は $1-2=-1$ であり、その二乗絶対値は 1 である。
:::
::::

::::exercise{level=advanced}
入射光の位相を半周期帯ごとに $0,\pi$ と切り替える二値位相板を作る。半径二乗に比例する無次元座標 $t$ について、一周期を $[0,2)$ とする。一次焦点の理想パワー効率を、交互に遮る振幅帯板と比べよ。

:::solution
二値位相板の透過振幅は $[0,1)$ で $1$、$[1,2)$ で $-1$。一次の係数は $\frac12\int_0^2 a(t)e^{-i\pi t}\dd t$ である。後半を $t-1$ で置換すると、透過率の負号と指数の負号が打ち消され、前半と同じ積分が加わる。従って係数は $2/(i\pi)$、理想効率は $4/\pi^2$。後半を遮る振幅帯板の係数はその半分なので効率は $1/\pi^2$ であり、位相板は四倍になる。
:::
::::

## 確認テスト

::::quiz
第一暗線が小角にある Fraunhofer 領域で、矩形開口の幅を半分にすると、この暗線の角はどうなるか。

- [x] 約二倍になる。
- [ ] 半分になる。
- [ ] 変わらない。

:::explanation
角が小さいと $\sin\theta\simeq\theta$。第一零点の条件は $a\sin\theta=\lambda$ なので、$a$ を半分にすると角は約二倍になる。
:::
::::

::::quiz{answer="4" tol="1e-6"}
Fresnel 近似で、$N_F=1$ の円孔の軸上照度を基準平面波の倍数で答えよ。

:::explanation
振幅比は $1-e^{i\pi}=2$ であり、照度比は絶対値の二乗で 4 である。
:::
::::

::::quiz
スカラー Babinet の原理で足す量はどれか。

- [x] 複素振幅
- [ ] 強度
- [ ] 強度の平方根の絶対値

:::explanation
伝播は場に対し線形である。したがって照明と二透過関数から得る場の位相を保ったまま足す。
:::
::::

::::quiz{answer="1000" tol="1e-6"}
$m=2,N=500$ の理想格子の分解能を答えよ。

:::explanation
第一零点までの幅と角分散を比べて $\abs{m}N$ を得るので $2(500)=1000$ である。
:::
::::

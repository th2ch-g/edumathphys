---
title: Fourier 光学
description: 横方向のスペクトルを用いて自由伝播、結像と空間フィルターを計算し、参照光による複素振幅の記録を導く。
notation:
  - symbol: '\boldsymbol f'
    meaning: 空間周波数
  - symbol: 'H_z'
    meaning: 距離 z の自由空間伝達関数
  - symbol: 'P'
    meaning: 瞳関数
  - symbol: 'h'
    meaning: 点入力の振幅応答
  - symbol: 'H_{\mathrm{inc}}'
    meaning: 強度結像の光学伝達関数
---

## この章の目標

- 角スペクトルから一様媒質の自由伝播を求める。
- レンズの前後焦点面を Fourier 変換で結ぶ。
- 振幅結像と強度結像の帯域を比較する。
- 位相差顕微鏡とホログラフィーで位相が測れる理由を説明する。

## 角スペクトル表示と自由空間の伝達関数

### 横座標の微分を周波数に換える
スカラー場を $U(\boldsymbol\rho,z)$ とし、横座標について $\widehat U(\boldsymbol f,z)=\int U(\boldsymbol\rho,z)e^{-2\pi i\boldsymbol f\cdot\boldsymbol\rho}\dd^2\rho$ と変換する。:term[空間周波数]{en="spatial frequency" yomi="くうかんしゅうはすう"} $\boldsymbol f$ は単位長さ当たりの周期数であり、規約は :ref[def-four-euclid-transform] と同じである。逆変換では符号を反転する。部分積分を二回行うと Helmholtz 方程式は $\partial_z^2\widehat U+(k^2-4\pi^2|\boldsymbol f|^2)\widehat U=0$ となる。半空間 $z>0$ で外向きまたは減衰する解を残すと、
$$
\widehat U(\boldsymbol f,z)=H_z(\boldsymbol f)\widehat U(\boldsymbol f,0),\qquad
H_z=e^{iz\sqrt{k^2-4\pi^2|\boldsymbol f|^2}} .
\label{eq-optics-angular-spectrum}
$$
平方根には、引数が非負なら非負の実根、引数が負なら虚部が正の根を取る。逆変換は傾斜の異なる平面波と減衰波の和であり、これが :term[角スペクトル]{en="angular spectrum" yomi="かくすぺくとる"} である。$k=2\pi/\lambda$、$\lambda|\boldsymbol f|\ll1$ なら、平方根を展開して $H_z\simeq e^{ikz-i\pi\lambda z|\boldsymbol f|^2}$。Gauss 積分で二次指数を逆変換すると :eqref[eq-optics-fresnel-diffraction] の畳み込み核を得る。

:::example[伝播の境界]
$\lambda=500\,\mathrm{nm}$ なら伝播成分の上限は $2\times10^6\,\mathrm{m^{-1}}$。周期 $1\,\mathrm{\mu m}$ の基本波は伝播できるが、周期 $0.25\,\mathrm{\mu m}$ では減衰する。
:::
:::example[逆伝播の不安定性]
$f=2/\lambda$ で $z=1/k$ なら前進因子は $e^{-\sqrt3}$、形式的な逆伝播は $e^{\sqrt3}$ 倍である。さらに高周波の測定雑音ほど強く増幅される。
:::

### 放射パワーの重み
均一媒質で Maxwell の横波条件を満たす、伝播波だけのベクトル角スペクトル $\widehat{\vb E}(\boldsymbol f)$ を考える。各成分は :eqref[eq-optics-plane] に従う。横面全体での平均流を積分すると、異なる横周波数間の交差項は Fourier の直交性で消える。従って法線パワーは $\int |\widehat{\vb E}|^2(k_z/k)\dd^2 f/(2\eta)$ である。明示的な $k_z/k$ は斜めの成分の法線への射影である。スカラーの $L^2$ ノルムをパワーと同一視する近軸近似は、この重みを 1 に置き換えている。

### Fresnel 核と伝達関数が同じ近似であること
近軸因子 $e^{-i\pi\lambda z f^2}$ の一次元逆変換を計算する。小さい正の $\epsilon$ を加えて $\int e^{-(\epsilon+i\pi\lambda z)f^2+2\pi ifx}\dd f$ とすると、平方完成した結果は $\sqrt{\pi/(\epsilon+i\pi\lambda z)}\exp[-\pi^2x^2/(\epsilon+i\pi\lambda z)]$。$\epsilon\to0^+$ で平方根の枝は正実部からの連続として固定する。二横方向の積は $(i\lambda z)^{-1}e^{ik\rho^2/(2z)}$。これに搬送位相 $e^{ikz}$ を掛ければ式 :eqref[eq-optics-fresnel-diffraction] の核が得られる。伝達関数の積が距離を加える性質から、途中の面を挿入して畳み込んでも同じ近似の距離を一度伝播した結果になる。

:::example[伝播因子の位相を求める]
$\lambda f=3/5$ の成分は $k_z/k=4/5$。一波長だけ軸方向へ伝播すると、その成分の位相は $8\pi/5$ 進む。一様な平面成分の $2\pi$ との差 $-2\pi/5$ が、二成分を重ねた像の変化を作る。
:::

### 近軸伝達関数の位相誤差を厳密な平方根で測る

横周期 $\Lambda$ の場に対し $\xi=\lambda/\Lambda<1$ とおく。軸位相を除いた厳密な位相は $kz(\sqrt{1-\xi^2}-1)$、Fresnel 伝達関数は $-kz\xi^2/2$ である。差の絶対値を有理化すると $kz\xi^4/[2(1+\sqrt{1-\xi^2})^2]$。従って近軸角を支配する $\xi$ の小ささに加えて、許容する誤差 $\varepsilon$ に対して $z$ が十分短いかを調べる。角分布に支持の上限 $\xi_{\max}<1$ があれば、この差は $\xi$ の増加関数なので最大の横周波数を入れればよい。伝達因子の差は $|e^{ia}-e^{ib}|\le |a-b|$ で抑えられるため、Plancherel の等式（:ref[thm-four-plancherel]）を適用すると、この最大の位相誤差が場の相対誤差 $\|U_{\rm exact}-U_{\rm F}\|_2/\|U(\cdot,0)\|_2$ の上限になる。

:::example[小さい横波数を近距離へ送る]
長い周期の場で位相誤差を測る。

入力は波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、横周期 $\Lambda\simeq 5\times 10^{-5}\,\mathrm m$、伝播距離 $z=\frac{1}{10}\,\mathrm m$ である。

分散面を波数 $k$ で規格化する。

$$
\begin{aligned}
\xi&=\frac{\lambda}{\Lambda},\\
\xi&\simeq \frac{\left(5\times 10^{-7}\right)}{\left(5\times 10^{-5}\right)}=\frac{1}{100}.
\end{aligned}
$$

位相誤差を実際の長さに換算する。

$$
\begin{aligned}
k&=\frac{2 \pi}{\lambda},\\
k/(\mathrm{m^{-1}})&\simeq \frac{2 \cdot \pi}{\left(5\times 10^{-7}\right)}\simeq 1.2566\times 10^{7}.
\end{aligned}
$$

縦波数を横波数の二次まで展開する。

$$
\begin{aligned}
\phi_{\rm F}&=- \frac{k \xi^{2} z}{2},\\
\phi_{\rm F}/(\mathrm{rad})&\simeq - \frac{\left(1.256637061\times 10^{7}\right) \cdot \left(\frac{1}{100}\right)^{2} \cdot \left(\frac{1}{10}\right)}{2}\simeq -62.832.
\end{aligned}
$$

差を有理化して、近い二数の引き算を避ける。

$$
\begin{aligned}
\delta_\phi&=\frac{k \xi^{4} z}{2 \left(\sqrt{1 - \xi^{2}} + 1\right)^{2}},\\
\delta_\phi/(\mathrm{rad})&\simeq \frac{\left(1.256637061\times 10^{7}\right) \cdot \left(\frac{1}{100}\right)^{4} \cdot \left(\frac{1}{10}\right)}{2 \cdot \left(\sqrt{1 - \left(\frac{1}{100}\right)^{2}} + 1\right)^{2}}\simeq 0.0015709.
\end{aligned}
$$

誤差が伝播距離に比例することを用いる。

$$
\begin{aligned}
z_{1}&=\frac{z}{\delta_\phi},\\
z_{1}/(\mathrm m)&\simeq \frac{\left(\frac{1}{10}\right)}{\left(0.001570874872\right)}\simeq 63.659.
\end{aligned}
$$

この距離では二次の位相までで精度を保てる。
:::

:::example[同じ角成分を遠くへ送る]
横周期は固定し、伝播距離だけを増やす。

入力は波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、横周期 $\Lambda\simeq 5\times 10^{-5}\,\mathrm m$、伝播距離 $z=100\,\mathrm m$ である。

分散面を波数 $k$ で規格化する。

$$
\begin{aligned}
\xi&=\frac{\lambda}{\Lambda},\\
\xi&\simeq \frac{\left(5\times 10^{-7}\right)}{\left(5\times 10^{-5}\right)}=\frac{1}{100}.
\end{aligned}
$$

位相誤差を実際の長さに換算する。

$$
\begin{aligned}
k&=\frac{2 \pi}{\lambda},\\
k/(\mathrm{m^{-1}})&\simeq \frac{2 \cdot \pi}{\left(5\times 10^{-7}\right)}\simeq 1.2566\times 10^{7}.
\end{aligned}
$$

縦波数を横波数の二次まで展開する。

$$
\begin{aligned}
\phi_{\rm F}&=- \frac{k \xi^{2} z}{2},\\
\phi_{\rm F}/(\mathrm{rad})&\simeq - \frac{\left(1.256637061\times 10^{7}\right) \cdot \left(\frac{1}{100}\right)^{2} \cdot \left(100\right)}{2}\simeq -62832.
\end{aligned}
$$

差を有理化して、近い二数の引き算を避ける。

$$
\begin{aligned}
\delta_\phi&=\frac{k \xi^{4} z}{2 \left(\sqrt{1 - \xi^{2}} + 1\right)^{2}},\\
\delta_\phi/(\mathrm{rad})&\simeq \frac{\left(1.256637061\times 10^{7}\right) \cdot \left(\frac{1}{100}\right)^{4} \cdot \left(100\right)}{2 \cdot \left(\sqrt{1 - \left(\frac{1}{100}\right)^{2}} + 1\right)^{2}}\simeq 1.5709.
\end{aligned}
$$

誤差が伝播距離に比例することを用いる。

$$
\begin{aligned}
z_{1}&=\frac{z}{\delta_\phi},\\
z_{1}/(\mathrm m)&\simeq \frac{\left(100\right)}{\left(1.570874872\right)}\simeq 63.659.
\end{aligned}
$$

小角でも伝播で位相誤差が蓄積する。
:::

:::example[短い周期を近距離へ送る]
横成分が大きい変調を使い、近軸核を選べるか判定する。

入力は波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、横周期 $\Lambda\simeq 5\times 10^{-6}\,\mathrm m$、伝播距離 $z=\frac{1}{10}\,\mathrm m$ である。

分散面を波数 $k$ で規格化する。

$$
\begin{aligned}
\xi&=\frac{\lambda}{\Lambda},\\
\xi&\simeq \frac{\left(5\times 10^{-7}\right)}{\left(5\times 10^{-6}\right)}=\frac{1}{10}.
\end{aligned}
$$

位相誤差を実際の長さに換算する。

$$
\begin{aligned}
k&=\frac{2 \pi}{\lambda},\\
k/(\mathrm{m^{-1}})&\simeq \frac{2 \cdot \pi}{\left(5\times 10^{-7}\right)}\simeq 1.2566\times 10^{7}.
\end{aligned}
$$

縦波数を横波数の二次まで展開する。

$$
\begin{aligned}
\phi_{\rm F}&=- \frac{k \xi^{2} z}{2},\\
\phi_{\rm F}/(\mathrm{rad})&\simeq - \frac{\left(1.256637061\times 10^{7}\right) \cdot \left(\frac{1}{10}\right)^{2} \cdot \left(\frac{1}{10}\right)}{2}\simeq -6283.2.
\end{aligned}
$$

差を有理化して、近い二数の引き算を避ける。

$$
\begin{aligned}
\delta_\phi&=\frac{k \xi^{4} z}{2 \left(\sqrt{1 - \xi^{2}} + 1\right)^{2}},\\
\delta_\phi/(\mathrm{rad})&\simeq \frac{\left(1.256637061\times 10^{7}\right) \cdot \left(\frac{1}{10}\right)^{4} \cdot \left(\frac{1}{10}\right)}{2 \cdot \left(\sqrt{1 - \left(\frac{1}{10}\right)^{2}} + 1\right)^{2}}\simeq 15.787.
\end{aligned}
$$

誤差が伝播距離に比例することを用いる。

$$
\begin{aligned}
z_{1}&=\frac{z}{\delta_\phi},\\
z_{1}/(\mathrm m)&\simeq \frac{\left(\frac{1}{10}\right)}{\left(15.78699742\right)}\simeq 0.0063343.
\end{aligned}
$$

絶対位相を保つ必要があるなら厳密な角スペクトル表示を使う。
:::

## Fraunhofer 回折と Fourier 変換

Fresnel 核の指数を $\rho^2-2\boldsymbol\rho\cdot\boldsymbol\rho'+\rho'^2$ に分け、開口内で $k\rho'^2/(2z)\ll1$ とする。最後の項を除くと
$$
U(\boldsymbol\rho,z)\simeq \frac{e^{ikz+ik\rho^2/(2z)}}{i\lambda z}\widehat U\!\left(\frac{\boldsymbol\rho}{\lambda z},0\right).
\label{eq-optics-far-transform}
$$
これは :term[遠方場]{en="far field" yomi="えんぽうば"} における振幅の対応である。強度を求めるなら外側の位相は消える。場を再び伝播させる場合には観測面の二次位相も必要である。平行移動で変換結果に一次位相だけが掛かるのは、照明を含む開口場全体を平行移動した場合である。開口だけを不均一な照明の下で動かすと、変換に入る場も変わる。逆に開口の振幅に一次位相を掛けると、分布そのものが遠方で移動する。

:::example[三つの回折次数]
有限な窓 $W$ 内の振幅を $W[1+m\cos(2\pi f_g x)]$ とする。指数表示を展開すると中心の $\widehat W$ と、$\pm f_g$ に移動した $m\widehat W/2$ の和になる。
:::
:::example[位相勾配の位置換算]
$\lambda=500\,\mathrm{nm}$、傾き $f_0=10^4\,\mathrm{m^{-1}}$ なら距離 $1\,\mathrm m$ の遠方で分布は $\lambda z f_0=5\,\mathrm{mm}$ 移動する。
:::

:::example[位相の情報が消える例]
実振幅の狭い開口を $x=x_0$ と $x=-x_0$ に別々に置くと、各遠方強度は同じである。しかし遠方の複素振幅は互いに $e^{-4\pi i f_x x_0}$ の因子だけ異なる。遠方強度だけを一枚記録しても、一般に入力位置を一意に求められない。
:::

### 回折強度から分かる開口の自己相関

一次元で有限支持の開口 $t(x)$ の Fourier 変換を $\hat t(f)=\int t(x)e^{-2\pi ifx}\dd x$ とする。強度 $|\hat t|^2$ を二重積分で書き、差 $u=x-y$ を変数に取れば、その逆変換は $A(u)=\int t(x+u)t^*(x)\dd x$ である。よって異なる開口が同じ自己相関を持てば強度だけでは区別できない。平行移動や複素共役を伴う反転は、この曖昧さを生む。横幅 $a$ の一様スリット二本、中心距離 $d>a$ なら、自己相関は原点の高さ $2a$ の三角形と、$\pm d$ に中心を持つ高さ $a$ の三角形を加えたもの。それぞれの三角形はスリットどうしの重なり長さから得る。

:::example[中心の自己相関]
差の座標をゼロにして全開口の二乗積分を確認する。

入力はスリット幅 $a\simeq 1\times 10^{-4}\,\mathrm m$、二スリットの距離 $d=\frac{1}{2500}\,\mathrm m$、自己相関のずれ $u=0\,\mathrm m$ である。

同じスリット上の二つの区間の共通部分を求める。

$$
\begin{aligned}
A&=\max\left(0, a - \left|{u}\right|\right),\\
A/(\mathrm m)&\simeq \max\left(0, \left(1\times 10^{-4}\right) - \left|{\left(0\right)}\right|\right)\simeq 1\times 10^{-4}.
\end{aligned}
$$

別々のスリットの位置差を含めて重なりを求める。

$$
\begin{aligned}
B&=\max\left(0, a - \left|{d - u}\right|\right),\\
B/(\mathrm m)&\simeq \max\left(0, \left(1\times 10^{-4}\right) - \left|{\left(\frac{1}{2500}\right) - \left(0\right)}\right|\right)=0.
\end{aligned}
$$

逆の向きに組み合わせた区間も同じ重なりの形になる。

$$
\begin{aligned}
C&=\max\left(0, a - \left|{d + u}\right|\right),\\
C/(\mathrm m)&\simeq \max\left(0, \left(1\times 10^{-4}\right) - \left|{\left(\frac{1}{2500}\right) + \left(0\right)}\right|\right)=0.
\end{aligned}
$$

各スリットの自己項を二回数え、二つの交差項を足す。

$$
\begin{aligned}
A_{\rm norm}&=\frac{2 A + B + C}{2 a},\\
A_{\rm norm}&\simeq \frac{2 \cdot \left(1\times 10^{-4}\right) + \left(0\right) + \left(0\right)}{2 \cdot \left(1\times 10^{-4}\right)}=1.
\end{aligned}
$$

中心の値はスリット二本の面積に比例する。
:::

:::example[開口間隔に相当する差]
別々のスリットが重なる差で相関を評価する。

入力はスリット幅 $a\simeq 1\times 10^{-4}\,\mathrm m$、二スリットの距離 $d=\frac{1}{2500}\,\mathrm m$、自己相関のずれ $u=\frac{1}{2500}\,\mathrm m$ である。

同じスリット上の二つの区間の共通部分を求める。

$$
\begin{aligned}
A&=\max\left(0, a - \left|{u}\right|\right),\\
A/(\mathrm m)&\simeq \max\left(0, \left(1\times 10^{-4}\right) - \left|{\left(\frac{1}{2500}\right)}\right|\right)=0.
\end{aligned}
$$

別々のスリットの位置差を含めて重なりを求める。

$$
\begin{aligned}
B&=\max\left(0, a - \left|{d - u}\right|\right),\\
B/(\mathrm m)&\simeq \max\left(0, \left(1\times 10^{-4}\right) - \left|{\left(\frac{1}{2500}\right) - \left(\frac{1}{2500}\right)}\right|\right)\simeq 1\times 10^{-4}.
\end{aligned}
$$

逆の向きに組み合わせた区間も同じ重なりの形になる。

$$
\begin{aligned}
C&=\max\left(0, a - \left|{d + u}\right|\right),\\
C/(\mathrm m)&\simeq \max\left(0, \left(1\times 10^{-4}\right) - \left|{\left(\frac{1}{2500}\right) + \left(\frac{1}{2500}\right)}\right|\right)=0.
\end{aligned}
$$

各スリットの自己項を二回数え、二つの交差項を足す。

$$
\begin{aligned}
A_{\rm norm}&=\frac{2 A + B + C}{2 a},\\
A_{\rm norm}&\simeq \frac{2 \cdot \left(0\right) + \left(1\times 10^{-4}\right) + \left(0\right)}{2 \cdot \left(1\times 10^{-4}\right)}=\frac{1}{2}.
\end{aligned}
$$

一本分だけが重なるため、中心の半分になる。
:::

:::example[開口の間隔より少しずらす]
交差相関の三角形の肩で値を求める。

入力はスリット幅 $a\simeq 1\times 10^{-4}\,\mathrm m$、二スリットの距離 $d=\frac{1}{2500}\,\mathrm m$、自己相関のずれ $u\simeq 4.5\times 10^{-4}\,\mathrm m$ である。

同じスリット上の二つの区間の共通部分を求める。

$$
\begin{aligned}
A&=\max\left(0, a - \left|{u}\right|\right),\\
A/(\mathrm m)&\simeq \max\left(0, \left(1\times 10^{-4}\right) - \left|{\left(4.5\times 10^{-4}\right)}\right|\right)=0.
\end{aligned}
$$

別々のスリットの位置差を含めて重なりを求める。

$$
\begin{aligned}
B&=\max\left(0, a - \left|{d - u}\right|\right),\\
B/(\mathrm m)&\simeq \max\left(0, \left(1\times 10^{-4}\right) - \left|{\left(\frac{1}{2500}\right) - \left(4.5\times 10^{-4}\right)}\right|\right)\simeq 5\times 10^{-5}.
\end{aligned}
$$

逆の向きに組み合わせた区間も同じ重なりの形になる。

$$
\begin{aligned}
C&=\max\left(0, a - \left|{d + u}\right|\right),\\
C/(\mathrm m)&\simeq \max\left(0, \left(1\times 10^{-4}\right) - \left|{\left(\frac{1}{2500}\right) + \left(4.5\times 10^{-4}\right)}\right|\right)=0.
\end{aligned}
$$

各スリットの自己項を二回数え、二つの交差項を足す。

$$
\begin{aligned}
A_{\rm norm}&=\frac{2 A + B + C}{2 a},\\
A_{\rm norm}&\simeq \frac{2 \cdot \left(0\right) + \left(5\times 10^{-5}\right) + \left(0\right)}{2 \cdot \left(1\times 10^{-4}\right)}=\frac{1}{4}.
\end{aligned}
$$

重なりが半分になり、中心と比較すると四分の一である。
:::

### 周期物体が近距離で自己像を作る

周期 $d$ の場 $U(x,0)=\sum_m a_m e^{2\pi i m x/d}$ を Fresnel 伝播させる。搬送位相を除くと第 $m$ 成分に $e^{-i\pi\lambda z m^2/d^2}$ が掛かる。距離 $z_{\rm T}=2d^2/\lambda$ では、全成分の指数が $-2\pi i$ の整数倍となり、元の複素振幅が戻る。半分の距離では $e^{-i\pi m^2}=(-1)^m$ なので、場は横に $d/2$ ずれたものになる。これは :term[Talbot 効果]{en="Talbot effect" yomi="たるぼこうか"} である。近軸の位相誤差は伝播する距離と帯域で決まるため、開口で格子を切った場合や高次成分を多く含む場合には、自己像の周辺や細部にずれが出る。

:::example[半周期だけ位置がずれる格子像]
$d=50\,\mathrm{\mu m}$、$\lambda=500\,\mathrm{nm}$ なら $z_{\rm T}=10\,\mathrm{mm}$ である。その半分の $5\,\mathrm{mm}$ で現れる像は $25\,\mathrm{\mu m}$ ずれる。
:::

## レンズの Fourier 変換作用

### 前焦点面と後焦点面
焦点距離 $F$ の薄レンズはレンズ上の横座標 $\boldsymbol v$ に位相 $-k|\boldsymbol v|^2/(2F)$ を加える。これは軸上焦点までの距離差 $|\boldsymbol v|^2/(2F)$ を打ち消す補正である。レンズ直前から後焦点面までの伝播では Fresnel 核の入力の二次項がこのレンズ位相で消える。したがって、距離を大きくしなくても開口を含む入力の変換が得られる。

前焦点面の点 $\boldsymbol u$、後焦点面の点 $\boldsymbol x$ を結ぶ場合、指数は $ik/(2F)$ に
$$
|\boldsymbol v-\boldsymbol u|^2-|\boldsymbol v|^2+|\boldsymbol x-\boldsymbol v|^2
=|\boldsymbol v-\boldsymbol u-\boldsymbol x|^2-2\boldsymbol u\cdot\boldsymbol x
$$
を掛けたものである。無限レンズの近似で $\boldsymbol v$ の積分は $i\lambda F$ となるので、残りは $e^{2ikF}\widehat U(\boldsymbol x/(\lambda F),0)/(i\lambda F)$。この配置では出力二次位相も残らない。有限レンズでは $\boldsymbol v$ 積分に瞳が入り、この関係には誤差が出る。

:::example[格子周期と焦点面]
$\lambda=500\,\mathrm{nm}$、$F=0.20\,\mathrm m$、格子周期 $100\,\mathrm{\mu m}$ なら一次成分は $\pm1.0\,\mathrm{mm}$ に来る。
:::
:::example[前焦点面の点]
理想系で軸外の点は傾斜した平面波になる。点の位置が $\boldsymbol u$ なら後焦点面の位相勾配は $-k\boldsymbol u/F$ である。
:::

### 変換を連続して適用する
二つのレンズの焦点距離を $F_1,F_2$ とし、第一の後焦点面と第二の前焦点面を共有させる。両瞳を無限とした近軸式では、 入力座標 $\boldsymbol u$ と共有面の座標 $\boldsymbol v$ の位相は $-k\boldsymbol u\cdot\boldsymbol v/F_1$、次は $-k\boldsymbol v\cdot\boldsymbol x/F_2$。$\boldsymbol v$ の無限面積分は $\boldsymbol u/F_1+\boldsymbol x/F_2=0$ を支持とするデルタ関数である。従って像は $\boldsymbol x=-(F_2/F_1)\boldsymbol u$ に現れる。共通位相を除いた振幅の倍率の絶対値は $F_1/F_2$ なので、面積倍率 $(F_2/F_1)^2$ と相殺し全パワーが保たれる。

:::example[二枚の焦点面の実際の間隔]
$F_1=50\,\mathrm{mm},F_2=100\,\mathrm{mm}$ のレンズ二枚を焦点面共有で使う。レンズ間隔は $150\,\mathrm{mm}$、入力は第一レンズの前 $50\,\mathrm{mm}$、出力は第二レンズの後 $100\,\mathrm{mm}$ である。倍率は $-2$ なので、入力の $+0.20\,\mathrm{mm}$ の点は出力の $-0.40\,\mathrm{mm}$ に写る。
:::

### Fourier 面で二つの次数を選んで像を作る

周期 $d$ の薄い格子を前焦点面に置く。格子の Fourier 級数の第 $m$ 項は小角で後焦点面の $X_m=\lambda f_1 m/d$ に集まる。第二レンズを共通焦点面から $f_2$ 離して置き、さらにその後焦点面で観測すれば、二回の Fourier 変換により物体座標は $M=-f_2/f_1$ 倍される。従って調整した焦点面で二つの次数 $m,n$ だけを通すと、像の振幅は二つの複素指数の和になる。強度の交差項の周波数は $|m-n|/(|M|d)$。同じ絶対値の正負一次だけを通した場合、強度周期は物体周期の $|M|/2$ 倍であり、透過する次数を変えると像の周期も変わり得る。

:::example[零次と一次を通す]
一方に基準となる直流成分を残す。

入力は波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、格子周期 $d\simeq 2\times 10^{-5}\,\mathrm m$、第一レンズの焦点距離 $f_{1}=\frac{1}{10}\,\mathrm m$、第二レンズの焦点距離 $f_{2}=\frac{1}{10}\,\mathrm m$、選ぶ次数 $m=0$、もう一つ選ぶ次数 $n=1$ である。

前焦点面の格子の横周波数を後焦点面の位置に換算する。

$$
\begin{aligned}
X_m&=\frac{f_{1} \lambda m}{d},\\
X_m/(\mathrm m)&\simeq \frac{\left(\frac{1}{10}\right) \cdot \left(5\times 10^{-7}\right) \cdot \left(0\right)}{\left(2\times 10^{-5}\right)}=0.
\end{aligned}
$$

もう一つの次数も同じ尺度で位置を決める。

$$
\begin{aligned}
X_n&=\frac{f_{1} \lambda n}{d},\\
X_n/(\mathrm m)&\simeq \frac{\left(\frac{1}{10}\right) \cdot \left(5\times 10^{-7}\right) \cdot \left(1\right)}{\left(2\times 10^{-5}\right)}=\frac{1}{400}.
\end{aligned}
$$

共通焦点面の開口位置を決める。

$$
\begin{aligned}
D&=\left|{X_m - X_n}\right|,\\
D/(\mathrm m)&= \left|{\left(0\right) - \left(\frac{1}{400}\right)}\right|=\frac{1}{400}.
\end{aligned}
$$

二回の Fourier 変換で反転と拡縮が起こる。

$$
\begin{aligned}
M&=- \frac{f_{2}}{f_{1}},\\
M&= - \frac{\left(\frac{1}{10}\right)}{\left(\frac{1}{10}\right)}=-1.
\end{aligned}
$$

像面で二振幅を重ねると差の周波数になる。

$$
\begin{aligned}
\Lambda_I&=\frac{d \left|{M}\right|}{\left|{m - n}\right|},\\
\Lambda_I/(\mathrm m)&\simeq \frac{\left(2\times 10^{-5}\right) \cdot \left|{\left(-1\right)}\right|}{\left|{\left(0\right) - \left(1\right)}\right|}\simeq 2\times 10^{-5}.
\end{aligned}
$$

像の周期は物体と等しい。
:::

:::example[正負一次だけを通す]
零次を遮り、対称な二次数を選ぶ。

入力は波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、格子周期 $d\simeq 2\times 10^{-5}\,\mathrm m$、第一レンズの焦点距離 $f_{1}=\frac{1}{10}\,\mathrm m$、第二レンズの焦点距離 $f_{2}=\frac{1}{10}\,\mathrm m$、選ぶ次数 $m=1$、もう一つ選ぶ次数 $n=-1$ である。

前焦点面の格子の横周波数を後焦点面の位置に換算する。

$$
\begin{aligned}
X_m&=\frac{f_{1} \lambda m}{d},\\
X_m/(\mathrm m)&\simeq \frac{\left(\frac{1}{10}\right) \cdot \left(5\times 10^{-7}\right) \cdot \left(1\right)}{\left(2\times 10^{-5}\right)}=\frac{1}{400}.
\end{aligned}
$$

もう一つの次数も同じ尺度で位置を決める。

$$
\begin{aligned}
X_n&=\frac{f_{1} \lambda n}{d},\\
X_n/(\mathrm m)&\simeq \frac{\left(\frac{1}{10}\right) \cdot \left(5\times 10^{-7}\right) \cdot \left(-1\right)}{\left(2\times 10^{-5}\right)}=- \frac{1}{400}.
\end{aligned}
$$

共通焦点面の開口位置を決める。

$$
\begin{aligned}
D&=\left|{X_m - X_n}\right|,\\
D/(\mathrm m)&= \left|{\left(\frac{1}{400}\right) - \left(- \frac{1}{400}\right)}\right|=\frac{1}{200}.
\end{aligned}
$$

二回の Fourier 変換で反転と拡縮が起こる。

$$
\begin{aligned}
M&=- \frac{f_{2}}{f_{1}},\\
M&= - \frac{\left(\frac{1}{10}\right)}{\left(\frac{1}{10}\right)}=-1.
\end{aligned}
$$

像面で二振幅を重ねると差の周波数になる。

$$
\begin{aligned}
\Lambda_I&=\frac{d \left|{M}\right|}{\left|{m - n}\right|},\\
\Lambda_I/(\mathrm m)&\simeq \frac{\left(2\times 10^{-5}\right) \cdot \left|{\left(-1\right)}\right|}{\left|{\left(1\right) - \left(-1\right)}\right|}\simeq 1\times 10^{-5}.
\end{aligned}
$$

同じ光学倍率で強度周期が半分になる。
:::

:::example[再結像で倍率を変える]
選ぶ次数は維持し、第二レンズの焦点距離を半分にする。

入力は波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、格子周期 $d\simeq 2\times 10^{-5}\,\mathrm m$、第一レンズの焦点距離 $f_{1}=\frac{1}{10}\,\mathrm m$、第二レンズの焦点距離 $f_{2}=\frac{1}{20}\,\mathrm m$、選ぶ次数 $m=1$、もう一つ選ぶ次数 $n=-1$ である。

前焦点面の格子の横周波数を後焦点面の位置に換算する。

$$
\begin{aligned}
X_m&=\frac{f_{1} \lambda m}{d},\\
X_m/(\mathrm m)&\simeq \frac{\left(\frac{1}{10}\right) \cdot \left(5\times 10^{-7}\right) \cdot \left(1\right)}{\left(2\times 10^{-5}\right)}=\frac{1}{400}.
\end{aligned}
$$

もう一つの次数も同じ尺度で位置を決める。

$$
\begin{aligned}
X_n&=\frac{f_{1} \lambda n}{d},\\
X_n/(\mathrm m)&\simeq \frac{\left(\frac{1}{10}\right) \cdot \left(5\times 10^{-7}\right) \cdot \left(-1\right)}{\left(2\times 10^{-5}\right)}=- \frac{1}{400}.
\end{aligned}
$$

共通焦点面の開口位置を決める。

$$
\begin{aligned}
D&=\left|{X_m - X_n}\right|,\\
D/(\mathrm m)&= \left|{\left(\frac{1}{400}\right) - \left(- \frac{1}{400}\right)}\right|=\frac{1}{200}.
\end{aligned}
$$

二回の Fourier 変換で反転と拡縮が起こる。

$$
\begin{aligned}
M&=- \frac{f_{2}}{f_{1}},\\
M&= - \frac{\left(\frac{1}{20}\right)}{\left(\frac{1}{10}\right)}=- \frac{1}{2}.
\end{aligned}
$$

像面で二振幅を重ねると差の周波数になる。

$$
\begin{aligned}
\Lambda_I&=\frac{d \left|{M}\right|}{\left|{m - n}\right|},\\
\Lambda_I/(\mathrm m)&\simeq \frac{\left(2\times 10^{-5}\right) \cdot \left|{\left(- \frac{1}{2}\right)}\right|}{\left|{\left(1\right) - \left(-1\right)}\right|}\simeq 5\times 10^{-6}.
\end{aligned}
$$

共通焦点面で選ぶ位置は同じだが、像の強度周期は縮む。
:::

## コヒーレント結像と瞳関数

二枚のレンズで変換を二回行うと像が反転する。像座標を反転し倍率を物体に換算した上で、途中の Fourier 面の複素振幅透過率を $P$ とする。これをこの配置での :term[瞳関数]{en="pupil function" yomi="ひとみかんすう"} と呼ぶ。:term[コヒーレント結像]{en="coherent imaging" yomi="こひーれんとけつぞう"} では $\widehat U_i=P\widehat U_o$ である。畳み込み定理（:ref[thm-four-euclid-convolution]）により
$$
U_i=h*U_o,\qquad h=\mathcal F^{-1}P,\qquad
(h*U_o)(\boldsymbol\rho)=\int h(\boldsymbol\rho-\boldsymbol\rho')U_o(\boldsymbol\rho')\dd^2\rho' .
\label{eq-optics-coherent-imaging}
$$
瞳の角度範囲を $\alpha$ とすると、物体側の周波数上限は $\sin\alpha/\lambda=\mathrm{NA}/\lambda_0$。波面収差 $W$ は $P$ に位相 $e^{ik_0W}$ として入る。瞳を移動するだけなら像の複素振幅応答には傾きが生じる。

:::example[二つの物点]
対称な点応答で等振幅逆相の二物点を写すと、像の中点で振幅が消える。各点の応答が中点で非零なら、像強度の和はそこで正になる。
:::
:::example[遮断周波数]
$\mathrm{NA}=0.10,\lambda_0=500\,\mathrm{nm}$ なら上限は $2.0\times10^5\,\mathrm{m^{-1}}$。周期 $4\,\mathrm{\mu m}$ の成分は理想の円瞳を通らない。
:::

:::example[一次成分を失った格子]
物体の振幅を $1+\epsilon\cos(2\pi f_0x)$ とする。円瞳が零次のみを通すと像の振幅は一様になる。三成分とも通る理想系の強度は $1+2\epsilon\cos(2\pi f_0x)+\epsilon^2\cos^2(2\pi f_0x)$。ここに現れる $2f_0$ は像の検出で作られた高調波であり、物体振幅の $2f_0$ 成分が瞳を通ったことを意味しない。
:::

### 斜め照明で物体スペクトルをずらす
薄い物体の透過率を $t(\boldsymbol\rho)$、入射平面波を $e^{2\pi i\boldsymbol s\cdot\boldsymbol\rho}$ とする。物体直後の変換は $\widehat t(\boldsymbol f-\boldsymbol s)$。従って通る物体周波数を $\boldsymbol q=\boldsymbol f-\boldsymbol s$ で表すと $P(\boldsymbol q+\boldsymbol s)$ が掛かる。軸上照明では瞳外だった高周波が、適切な $\boldsymbol s$ の符号で内部へ入る。反対側の周波数は同じ露光で失われ得る。複数の照明方向で測った場を組み合わせるには、振幅の位相と、照明から生じる座標ずれの校正が必要である。

この考え方で、格子の零次と一次をともに通す条件を求められる。一次元で瞳が $|f|\le f_c$ を通し、物体の周期が $d$ なら、照明 $s$ に対する零次と正一次は $s$ と $s+1/d$ にある。両者を収める $s$ は $1/d\le2f_c$ の場合に存在し、間隔がちょうど $2f_c$ なら両者を瞳の対向する端に置く。一方、正面照明では $1/d\le f_c$ が必要である。これが同じ対物瞳で斜め照明により細かい周期を観察できる理由である。実際には照明側が供給できる $|s|$ の範囲も課し、少なくとも二成分が通ったかを調べる。

### 瞳の二値位相が像の対称性を変える

横振幅を表す線形なコヒーレント像で、瞳を $0<f<F$ で $H(f)=1$、$-F<f<0$ で $H(f)=e^{i\phi}$ とする。点像振幅は正負の積分をそれぞれ実行して $h(x)=F\operatorname{sinc}(\pi Fx)[e^{i\pi Fx}+e^{i\phi-i\pi Fx}]$ となる。従って瞳の大きさを変えなくても、中心の強度は $\cos^2(\phi/2)$ 倍になり、$\phi=\pi$ では消える。ここで $|h(x)|^2$ は単一物点を測った関数である。複数の同相物点の像は $h$ を足すため、この強度関数との畳み込みでは表せない。

:::example[分割しない一様瞳]
両半分の位相を合わせ、像の中心を選ぶ。

入力は片側の遮断周波数 $F=100000\,\mathrm{m^{-1}}$、負周波数側の位相 $\phi=0\,\mathrm{rad}$、像の横位置 $x=0\,\mathrm m$ である。

半瞳を積分した sinc の引数を求める。

$$
\begin{aligned}
v&=\pi F x,\\
v&= \pi \cdot \left(100000\right) \cdot \left(0\right)=0.
\end{aligned}
$$

全開口の無収差中心振幅 $2F$ を基準にする。

$$
\begin{aligned}
h_+&=\frac{e^{i v} \operatorname{sinc}{\left(v \right)}}{2},\\
h_+&= \frac{e^{i \cdot \left(0\right)} \cdot \operatorname{sinc}{\left(\left(0\right) \right)}}{2}=\frac{1}{2}.
\end{aligned}
$$

位置による位相は正側と反対で、瞳素子の位相をさらに加える。

$$
\begin{aligned}
h_-&=\frac{e^{i \left(\phi - v\right)} \operatorname{sinc}{\left(v \right)}}{2},\\
h_-&= \frac{e^{i \cdot \left(\left(0\right) - \left(0\right)\right)} \cdot \operatorname{sinc}{\left(\left(0\right) \right)}}{2}=\frac{1}{2}.
\end{aligned}
$$

二つの半瞳の振幅を重ねてから絶対値の二乗を取る。

$$
\begin{aligned}
I/I_{\rm ref}&=\left|{h_+ + h_-}\right|^{2},\\
I/I_{\rm ref}&= \left|{\left(\frac{1}{2}\right) + \left(\frac{1}{2}\right)}\right|^{2}=1.
\end{aligned}
$$

無収差瞳の振幅に一致する。
:::

:::example[位相段差で中心を消す]
片側だけ半周期遅らせる。

入力は片側の遮断周波数 $F=100000\,\mathrm{m^{-1}}$、負周波数側の位相 $\phi\simeq 3.1416\,\mathrm{rad}$、像の横位置 $x=0\,\mathrm m$ である。

半瞳を積分した sinc の引数を求める。

$$
\begin{aligned}
v&=\pi F x,\\
v&= \pi \cdot \left(100000\right) \cdot \left(0\right)=0.
\end{aligned}
$$

全開口の無収差中心振幅 $2F$ を基準にする。

$$
\begin{aligned}
h_+&=\frac{e^{i v} \operatorname{sinc}{\left(v \right)}}{2},\\
h_+&= \frac{e^{i \cdot \left(0\right)} \cdot \operatorname{sinc}{\left(\left(0\right) \right)}}{2}=\frac{1}{2}.
\end{aligned}
$$

位置による位相は正側と反対で、瞳素子の位相をさらに加える。

$$
\begin{aligned}
h_-&=\frac{e^{i \left(\phi - v\right)} \operatorname{sinc}{\left(v \right)}}{2},\\
h_-&\simeq \frac{e^{i \cdot \left(\left(3.141592654\right) - \left(0\right)\right)} \cdot \operatorname{sinc}{\left(\left(0\right) \right)}}{2}=- \frac{1}{2}.
\end{aligned}
$$

二つの半瞳の振幅を重ねてから絶対値の二乗を取る。

$$
\begin{aligned}
I/I_{\rm ref}&=\left|{h_+ + h_-}\right|^{2},\\
I/I_{\rm ref}&= \left|{\left(\frac{1}{2}\right) + \left(- \frac{1}{2}\right)}\right|^{2}=0.
\end{aligned}
$$

両半分は等しい強度を通しても、中心で逆位相になる。
:::

:::example[中心の外で二成分を合わせる]
片側の遅相を維持して中心から移る。

入力は片側の遮断周波数 $F=100000\,\mathrm{m^{-1}}$、負周波数側の位相 $\phi\simeq 3.1416\,\mathrm{rad}$、像の横位置 $x\simeq 5\times 10^{-6}\,\mathrm m$ である。

半瞳を積分した sinc の引数を求める。

$$
\begin{aligned}
v&=\pi F x,\\
v&\simeq \pi \cdot \left(100000\right) \cdot \left(5\times 10^{-6}\right)\simeq 1.5708.
\end{aligned}
$$

全開口の無収差中心振幅 $2F$ を基準にする。

$$
\begin{aligned}
h_+&=\frac{e^{i v} \operatorname{sinc}{\left(v \right)}}{2},\\
h_+&\simeq \frac{e^{i \cdot \left(1.570796327\right)} \cdot \operatorname{sinc}{\left(\left(1.570796327\right) \right)}}{2}\simeq 0.31831i.
\end{aligned}
$$

位置による位相は正側と反対で、瞳素子の位相をさらに加える。

$$
\begin{aligned}
h_-&=\frac{e^{i \left(\phi - v\right)} \operatorname{sinc}{\left(v \right)}}{2},\\
h_-&\simeq \frac{e^{i \cdot \left(\left(3.141592654\right) - \left(1.570796327\right)\right)} \cdot \operatorname{sinc}{\left(\left(1.570796327\right) \right)}}{2}\simeq 0.31831i.
\end{aligned}
$$

二つの半瞳の振幅を重ねてから絶対値の二乗を取る。

$$
\begin{aligned}
I/I_{\rm ref}&=\left|{h_+ + h_-}\right|^{2},\\
I/I_{\rm ref}&\simeq \left|{\left(0.3183098862i\right) + \left(0.3183098862i\right)}\right|^{2}\simeq 0.40528.
\end{aligned}
$$

位置による位相差が段差を打ち消す。エネルギーは瞳で失われず、点像の形が変わる。
:::

## インコヒーレント結像と光学伝達関数

### 点像強度と瞳の自己相関
光源点の振幅が互いに無相関な場合、光源点を離散的に置くと像強度の和で異なる点の交差項が消える。連続模型では物体面の場をデルタ相関を持つ一般化確率場とし、この条件を $\langle U_o(\boldsymbol\rho_1)U_o^*(\boldsymbol\rho_2)\rangle=I_o(\boldsymbol\rho_1)\delta^{(2)}(\boldsymbol\rho_1-\boldsymbol\rho_2)$ と置く。$I_o$ はデルタ関数の係数であり、点ごとの有限な場の分散ではない。式 :eqref[eq-optics-coherent-imaging] を絶対値二乗した二重積分に代入すると、デルタ関数により積分が一つになって $I_i=|h|^2*I_o$ を得る。これが :term[インコヒーレント結像]{en="incoherent imaging" yomi="いんこひーれんとけつぞう"} である。$h$ と共役の逆変換を掛けると、出力周波数は瞳の二周波数の差になる。零周波数で規格化した :term[光学伝達関数]{en="optical transfer function" yomi="こうがくでんたつかんすう"} は
$$
H_{\mathrm{inc}}(\boldsymbol f)=
\frac{\int P(\boldsymbol s+\boldsymbol f/2)P^*(\boldsymbol s-\boldsymbol f/2)\dd^2s}{\int |P(\boldsymbol s)|^2\dd^2s} .
\label{eq-optics-otf}
$$
その絶対値が :term[変調伝達関数]{en="modulation transfer function" yomi="へんちょうでんたつかんすう"}（MTF）で、位相は強度変調の位置のずれを表す。半径 $f_c$ の一様円瞳では分子は二円の重なり面積である。$u=|\boldsymbol f|/(2f_c)$ に対し、片円の弓形は扇形 $f_c^2\arccos u$ から三角形 $f_c^2 u\sqrt{1-u^2}$ を除いたもの。従って $0\le u\le1$ で $H_{\mathrm{inc}}=2(\arccos u-u\sqrt{1-u^2})/\pi$、$u>1$ でゼロとなる。

:::example[中間周波数]
$u=1/2$ の伝達は約 $0.391$。物体の正弦強度変調が 10 パーセントなら、像では 3.91 パーセントとなる。
:::
:::example[位相反転]
収差で $H_{\mathrm{inc}}$ が負の実数となれば、実正弦の明暗は反転する。絶対値だけではこの違いを表示できない。
:::

### 二次位相による伝達の低下
一次元の矩形瞳を $P(s)=1$（$|s|\le f_c$）の内側に位相 $e^{i a s^2}$ を掛けたものとする。二つの瞳の重なりは $|s|\le f_c-|f|/2$ であり、その内部で位相差は $a[(s+f/2)^2-(s-f/2)^2]=2asf$ となる。重なりの半幅を $w$ とすると :eqref[eq-optics-otf] の分子は $2w\operatorname{sinc}(2afw)$、従って伝達は $(w/f_c)\operatorname{sinc}(2afw)$ である。幾何学的な帯域は変わらないが、内部での位相干渉が伝達を下げる。特定の周波数で $\operatorname{sinc}$ の符号も変わり得る。

### 照明方向を平均した透過物体の像
薄物体を方向 $\boldsymbol s$ の互いに独立な平面波で照明し、方向強度を $J(\boldsymbol s)$ とする。各方向での像振幅は、絶対値に効かない位相 $e^{2\pi i\boldsymbol s\cdot\boldsymbol \rho}$ を除いて $\int P(\boldsymbol q+\boldsymbol s)\widehat t(\boldsymbol q)e^{2\pi i\boldsymbol q\cdot\boldsymbol\rho}\dd^2q$。二乗絶対値を方向で足すと
$$
\begin{aligned}
I(\boldsymbol\rho)&=\iint K(\boldsymbol q_1,\boldsymbol q_2)
\widehat t(\boldsymbol q_1)\widehat t^*(\boldsymbol q_2)
e^{2\pi i(\boldsymbol q_1-\boldsymbol q_2)\cdot\boldsymbol\rho}\dd^2q_1\dd^2q_2,\\
K(\boldsymbol q_1,\boldsymbol q_2)&=\int J(\boldsymbol s)
P(\boldsymbol q_1+\boldsymbol s)P^*(\boldsymbol q_2+\boldsymbol s)\dd^2s .
\end{aligned}
$$
一方向しかなければ $K$ は因数分解でき、コヒーレント結像を再現する。方向分布が、関係する二つの移動瞳の重なりの全域で一定なら、変数をずらして $K$ は周波数差だけの関数になる。このとき $\widehat t(\boldsymbol q_1)\widehat t^*(\boldsymbol q_2)$ を差周波数一定で積分したものは $|t|^2$ の Fourier 変換なので、強度結像の式が現れる。有限な照明角でその中間にあるとき、単一の振幅伝達または強度伝達だけでは一般に記述できない。

:::example[一次元の二つのカットオフ]
実で平坦な矩形瞳の半幅を $f_c=40\,\mathrm{mm^{-1}}$ とする。$f=60\,\mathrm{mm^{-1}}$ の物体強度の正弦は、瞳の重なりが $2f_c-f=20\,\mathrm{mm^{-1}}$ なので変調が $1/4$ 残る。同周波数の物体振幅をコヒーレントに写す場合は振幅瞳の外側となる。
:::

### デフォーカスのある矩形瞳の OTF

一次元で $P(f)=\boldsymbol 1_{[-F,F]}(f)e^{ib f^2}$ とする。周波数差 $s\ge0$ の規格化自己相関は、$s<2F$ なら重なり区間 $[-F+s/2,F-s/2]$ 上で位相差 $b[(f+s/2)^2-(f-s/2)^2]=2bsf$ を積分し、$2F$ で割った値。従って OTF は $(1-s/(2F))\operatorname{sinc}(bs(2F-s))$ である。$s=0$ では一、$s\ge2F$ ではゼロで、デフォーカスで支持が広がることはない。区間内部で sinc が負になると正弦状の強度模様の明暗が逆転する。MTF はこの符号を除いた大きさである。

:::example[焦点が合った場合]
瞳の遮断周波数に等しい入力変調を通す。

入力は瞳の片側遮断周波数 $F=200000\,\mathrm{m^{-1}}$、入力強度の空間周波数 $s=200000\,\mathrm{m^{-1}}$、瞳縁での二次位相 $\phi_{\rm edge}=0\,\mathrm{rad}$ である。

縁での位相を横周波数の二乗で割る。

$$
\begin{aligned}
b&=\frac{\phi_{\rm edge}}{F^{2}},\\
b/(\mathrm{m^2})&= \frac{\left(0\right)}{\left(200000\right)^{2}}=0.
\end{aligned}
$$

瞳が重なる幅を全幅で割る。

$$
\begin{aligned}
a&=1 - \frac{s}{2 F},\\
a&= 1 - \frac{\left(200000\right)}{2 \cdot \left(200000\right)}=\frac{1}{2}.
\end{aligned}
$$

自己相関の位相差を重なりの外端で評価する。

$$
\begin{aligned}
v&=b s \left(2 F - s\right),\\
v/(\mathrm{rad})&= \left(0\right) \cdot \left(200000\right) \cdot \left(2 \cdot \left(200000\right) - \left(200000\right)\right)=0.
\end{aligned}
$$

重なりの幅と区間内での振動の平均を掛ける。

$$
\begin{aligned}
H&=a \operatorname{sinc}{\left(v \right)},\\
H&= \left(\frac{1}{2}\right) \cdot \operatorname{sinc}{\left(\left(0\right) \right)}=\frac{1}{2}.
\end{aligned}
$$

像と物体の強度変調の振幅比には絶対値を取る。

$$
\begin{aligned}
M&=\left|{H}\right|,\\
M&= \left|{\left(\frac{1}{2}\right)}\right|=\frac{1}{2}.
\end{aligned}
$$

OTF は瞳の重なりだけで決まる。
:::

:::example[位相差で相殺する場合]
同じ周波数でデフォーカスを増し、区間内で位相が一回転するようにする。

入力は瞳の片側遮断周波数 $F=200000\,\mathrm{m^{-1}}$、入力強度の空間周波数 $s=200000\,\mathrm{m^{-1}}$、瞳縁での二次位相 $\phi_{\rm edge}\simeq 3.1416\,\mathrm{rad}$ である。

縁での位相を横周波数の二乗で割る。

$$
\begin{aligned}
b&=\frac{\phi_{\rm edge}}{F^{2}},\\
b/(\mathrm{m^2})&\simeq \frac{\left(3.141592654\right)}{\left(200000\right)^{2}}\simeq 7.854\times 10^{-11}.
\end{aligned}
$$

瞳が重なる幅を全幅で割る。

$$
\begin{aligned}
a&=1 - \frac{s}{2 F},\\
a&= 1 - \frac{\left(200000\right)}{2 \cdot \left(200000\right)}=\frac{1}{2}.
\end{aligned}
$$

自己相関の位相差を重なりの外端で評価する。

$$
\begin{aligned}
v&=b s \left(2 F - s\right),\\
v/(\mathrm{rad})&\simeq \left(7.853981634\times 10^{-11}\right) \cdot \left(200000\right) \cdot \left(2 \cdot \left(200000\right) - \left(200000\right)\right)\simeq 3.1416.
\end{aligned}
$$

重なりの幅と区間内での振動の平均を掛ける。

$$
\begin{aligned}
H&=a \operatorname{sinc}{\left(v \right)},\\
H&\simeq \left(\frac{1}{2}\right) \cdot \operatorname{sinc}{\left(\left(3.141592654\right) \right)}=0.
\end{aligned}
$$

像と物体の強度変調の振幅比には絶対値を取る。

$$
\begin{aligned}
M&=\left|{H}\right|,\\
M&= \left|{\left(0\right)}\right|=0.
\end{aligned}
$$

重なりはあるが変調は伝わらない。
:::

:::example[伝達が負になる場合]
デフォーカスをさらに増す。

入力は瞳の片側遮断周波数 $F=200000\,\mathrm{m^{-1}}$、入力強度の空間周波数 $s=200000\,\mathrm{m^{-1}}$、瞳縁での二次位相 $\phi_{\rm edge}\simeq 4.7124\,\mathrm{rad}$ である。

縁での位相を横周波数の二乗で割る。

$$
\begin{aligned}
b&=\frac{\phi_{\rm edge}}{F^{2}},\\
b/(\mathrm{m^2})&\simeq \frac{\left(4.71238898\right)}{\left(200000\right)^{2}}\simeq 1.1781\times 10^{-10}.
\end{aligned}
$$

瞳が重なる幅を全幅で割る。

$$
\begin{aligned}
a&=1 - \frac{s}{2 F},\\
a&= 1 - \frac{\left(200000\right)}{2 \cdot \left(200000\right)}=\frac{1}{2}.
\end{aligned}
$$

自己相関の位相差を重なりの外端で評価する。

$$
\begin{aligned}
v&=b s \left(2 F - s\right),\\
v/(\mathrm{rad})&\simeq \left(1.178097245\times 10^{-10}\right) \cdot \left(200000\right) \cdot \left(2 \cdot \left(200000\right) - \left(200000\right)\right)\simeq 4.7124.
\end{aligned}
$$

重なりの幅と区間内での振動の平均を掛ける。

$$
\begin{aligned}
H&=a \operatorname{sinc}{\left(v \right)},\\
H&\simeq \left(\frac{1}{2}\right) \cdot \operatorname{sinc}{\left(\left(4.71238898\right) \right)}\simeq -0.1061.
\end{aligned}
$$

像と物体の強度変調の振幅比には絶対値を取る。

$$
\begin{aligned}
M&=\left|{H}\right|,\\
M&\simeq \left|{\left(-0.1061032954\right)}\right|\simeq 0.1061.
\end{aligned}
$$

MTF は正である一方、OTF の負符号が明暗の反転を示す。
:::

## 空間フィルタリングと位相差顕微鏡

Fourier 面で中心付近を遮れば一様な照明を止めて散乱成分を見られる。これを :term[暗視野]{en="dark field" yomi="あんしや"} という。逆に中心だけを通すと高周波が消え、細部が見えなくなる。

弱い実位相 $\phi$ を加える試料の透過は $e^{i\phi}\simeq 1+i\phi$ である。位相の空間平均をゼロとする。合焦した無収差の実偶対称な瞳では、伝達された位相分布も実数なので、像強度に位相の一次項は出ない。試料の帯域をすべて通す理想系なら、この関係は $\Re(i\phi)=0$ で表せる。:term[位相差顕微鏡]{en="phase contrast microscope" yomi="いそうさけんびきょう"} では零次へ $+\pi/2$ を加えて、理想系での像振幅を $i+i\phi$ にする。振幅の絶対値を二乗すると $1+2\phi+O(\phi^2)$ になり、一次変化が見える。位相板が有限なら零次の周りに入る低周波成分も影響を受ける。

:::example[位相差の測定]
$\phi=0.03\,\mathrm{rad}$ の部分は理想の $+\pi/2$ の位相板で、強度が一次で $0.06$ 増加する。
:::
:::example[非対称なフィルター]
正弦格子の正の一次だけを零次と通すと、振幅には $e^{2\pi i f_g x}$ が残る。負の一次も通す場合と違い、変調成分は実関数にならない。
:::

### 照明を抑えた位相差
一様成分の透過振幅を $\eta e^{i\theta}$（$0<\eta\le1$）とし、散乱成分 $i\phi$ の透過は 1 とする。一次強度は $\eta^2+2\eta\phi\sin\theta$。強度変化と背景の比は $2\phi\sin\theta/\eta$ である。零次を弱くすると相対コントラストは増えるが、捨てた $\phi^2$ が $\eta|\phi|$ に対して小さくなければ一次式は使えない。暗視野は $\eta=0$ の場合である。散乱成分をすべて通すなら、強度の先頭は $\phi^2$ となるため位相の符号は区別できない。有限瞳では、通過した散乱場の絶対値二乗を用いる。

:::example[零次を弱めた場合の測定値]
$\phi=0.020,\eta=0.5,\theta=\pi/2$ の透過物体では、背景は $0.25$、一次の局所増加は $2(0.5)(0.020)=0.020$ である。背景からの割合は約 8 パーセント。捨てた $\phi^2=0.0004$ は一次増加の 2 パーセントなので、この比で近似の誤差を見積もれる。
:::

### 位相物体を二つの対称な露光で測る

弱い位相物体の場を $1+i\varphi(x)+O(\varphi^2)$ とする。零次だけに位相 $\pm\pi/2$ と実振幅透過 $\eta$ を与え、残りの成分を通す。平均の位相をゼロに選べば像場は $\pm i\eta+i\varphi$ である。強度差はこの一次場モデルでは $4\eta\varphi$、和は $2(\eta^2+\varphi^2)$。実際の透過場には二次補正があるため、強度の和を二次まで求めるときは通過した零次振幅も補正する必要がある。位相以外の変調がなく残りの成分をすべて通す場合、$D$ を $e^{i\varphi}$ の空間平均、$s=e^{i\varphi}-D$ とすると、厳密な強度差は $|i\eta D+s|^2-|-i\eta D+s|^2=4\eta\Im(D^*s)$ である。平均ゼロの小さい $\varphi$ なら $\Im(D^*s)=\varphi+O(\|\varphi\|_\infty^3)$ となり、二次の実摂動は対称な差から消える。$\eta$ を小さくすれば背景に対する変調は大きいが、絶対の強度差は小さくなる。$|\varphi|/\eta$ が大きい場合は、強度単独を一次に近似すると負の強度を予測し得る。

:::example[強い参照で位相を読む]
小さい位相差を十分な零次光と重ねる。

入力は零次の振幅透過率 $\eta=1$、ある点の小さい位相 $\varphi=\frac{1}{50}\,\mathrm{rad}$ である。

一次場を正符号の参照と足す。

$$
\begin{aligned}
I_+&=\left(\eta + \varphi\right)^{2},\\
I_+&= \left(\left(1\right) + \left(\frac{1}{50}\right)\right)^{2}=\frac{2601}{2500}.
\end{aligned}
$$

参照の位相を反転して同じ点を測る。

$$
\begin{aligned}
I_-&=\left(- \eta + \varphi\right)^{2},\\
I_-&= \left(- \left(1\right) + \left(\frac{1}{50}\right)\right)^{2}=\frac{2401}{2500}.
\end{aligned}
$$

共通の背景と一次場の二乗の共通部分が消える。

$$
\begin{aligned}
D&=I_+ - I_-,\\
D&= \left(\frac{2601}{2500}\right) - \left(\frac{2401}{2500}\right)=\frac{2}{25}.
\end{aligned}
$$

既知の参照振幅を用いて位相へ戻す。

$$
\begin{aligned}
\hat\varphi&=\frac{D}{4 \eta},\\
\hat\varphi/(\mathrm{rad})&= \frac{\left(\frac{2}{25}\right)}{4 \cdot \left(1\right)}=\frac{1}{50}.
\end{aligned}
$$

一枚の像を一次強度だけで近似できるか調べる。

$$
\begin{aligned}
R&=\left|{\frac{\varphi}{\eta}}\right|,\\
R&= \left|{\frac{\left(\frac{1}{50}\right)}{\left(1\right)}}\right|=\frac{1}{50}.
\end{aligned}
$$

背景に対する変調は小さいが、差で符号を取得できる。
:::

:::example[参照を減らす]
物体を保って零次を減光する。

入力は零次の振幅透過率 $\eta=\frac{1}{5}$、ある点の小さい位相 $\varphi=\frac{1}{50}\,\mathrm{rad}$ である。

一次場を正符号の参照と足す。

$$
\begin{aligned}
I_+&=\left(\eta + \varphi\right)^{2},\\
I_+&= \left(\left(\frac{1}{5}\right) + \left(\frac{1}{50}\right)\right)^{2}=\frac{121}{2500}.
\end{aligned}
$$

参照の位相を反転して同じ点を測る。

$$
\begin{aligned}
I_-&=\left(- \eta + \varphi\right)^{2},\\
I_-&= \left(- \left(\frac{1}{5}\right) + \left(\frac{1}{50}\right)\right)^{2}=\frac{81}{2500}.
\end{aligned}
$$

共通の背景と一次場の二乗の共通部分が消える。

$$
\begin{aligned}
D&=I_+ - I_-,\\
D&= \left(\frac{121}{2500}\right) - \left(\frac{81}{2500}\right)=\frac{2}{125}.
\end{aligned}
$$

既知の参照振幅を用いて位相へ戻す。

$$
\begin{aligned}
\hat\varphi&=\frac{D}{4 \eta},\\
\hat\varphi/(\mathrm{rad})&= \frac{\left(\frac{2}{125}\right)}{4 \cdot \left(\frac{1}{5}\right)}=\frac{1}{50}.
\end{aligned}
$$

一枚の像を一次強度だけで近似できるか調べる。

$$
\begin{aligned}
R&=\left|{\frac{\varphi}{\eta}}\right|,\\
R&= \left|{\frac{\left(\frac{1}{50}\right)}{\left(\frac{1}{5}\right)}}\right|=\frac{1}{10}.
\end{aligned}
$$

相対的な縞は強まるが、二露光の差の絶対値は減る。
:::

:::example[参照に近い位相振幅]
零次をさらに弱める。

入力は零次の振幅透過率 $\eta=\frac{3}{100}$、ある点の小さい位相 $\varphi=\frac{1}{50}\,\mathrm{rad}$ である。

一次場を正符号の参照と足す。

$$
\begin{aligned}
I_+&=\left(\eta + \varphi\right)^{2},\\
I_+&= \left(\left(\frac{3}{100}\right) + \left(\frac{1}{50}\right)\right)^{2}=\frac{1}{400}.
\end{aligned}
$$

参照の位相を反転して同じ点を測る。

$$
\begin{aligned}
I_-&=\left(- \eta + \varphi\right)^{2},\\
I_-&= \left(- \left(\frac{3}{100}\right) + \left(\frac{1}{50}\right)\right)^{2}\simeq 1\times 10^{-4}.
\end{aligned}
$$

共通の背景と一次場の二乗の共通部分が消える。

$$
\begin{aligned}
D&=I_+ - I_-,\\
D&\simeq \left(\frac{1}{400}\right) - \left(1\times 10^{-4}\right)=\frac{3}{1250}.
\end{aligned}
$$

既知の参照振幅を用いて位相へ戻す。

$$
\begin{aligned}
\hat\varphi&=\frac{D}{4 \eta},\\
\hat\varphi/(\mathrm{rad})&= \frac{\left(\frac{3}{1250}\right)}{4 \cdot \left(\frac{3}{100}\right)}=\frac{1}{50}.
\end{aligned}
$$

一枚の像を一次強度だけで近似できるか調べる。

$$
\begin{aligned}
R&=\left|{\frac{\varphi}{\eta}}\right|,\\
R&= \left|{\frac{\left(\frac{1}{50}\right)}{\left(\frac{3}{100}\right)}}\right|=\frac{2}{3}.
\end{aligned}
$$

一次場モデルの強度差は使えるが、一露光の強度の二次項は無視しにくい。
:::

### 有限な位相板での位相伝達

物体の透過率を $1-a(\boldsymbol\rho)+i\varphi(\boldsymbol\rho)$ と近似する。$a,\varphi$ は実で平均ゼロの微小変調とする。フィルター込みの振幅伝達を $H(\boldsymbol f)$、参照の振幅を $H_0=H(\boldsymbol 0)$ と置けば、像強度の一次変化は $H_0^*\delta U+H_0\delta U^*$ である。第一項を Fourier 変換すると $H_0^* H(\boldsymbol f)(-\widehat a+i\widehat\varphi)$、第二項には $H^*(-\boldsymbol f)$ が入る。従って
$$
\begin{aligned}
\widehat{\delta I}(\boldsymbol f)
&=-[H_0^*H(\boldsymbol f)+H_0H^*(-\boldsymbol f)]\widehat a(\boldsymbol f)\\
&\quad+i[H_0^*H(\boldsymbol f)-H_0H^*(-\boldsymbol f)]\widehat\varphi(\boldsymbol f).
\end{aligned}
$$
実偶対称な無収差瞳では第二行が消える。一方、参照だけに $+\pi/2$ の位相を加え、対象の側帯を振幅 1 で通す位相板では、第一行が消えて位相項は $2\widehat\varphi$ になる。有限な位相板では、参照と同じ位相が掛かった中心近傍に位相伝達の低い帯域ができる。

## ホログラフィー

### 参照波との交差項
物体波 $O$ と参照波 $R$ を重ね、強度 $I=|O|^2+|R|^2+OR^*+O^*R$ を感光面に記録する。:term[ホログラフィー]{en="holography" yomi="ほろぐらふぃー"} では交差項が位相情報を持つ。記録に比例する透過を $R$ で照明すると $(|O|^2+|R|^2)R+|R|^2O+R^2O^*$。第二項が元の物体波、第三項が共役の像を作る。$R$ に線形位相を加えれば交差項を Fourier 面の別々の位置に移せる。帯域の重なりがなければ片方を取り出し、既知の参照因子で割って複素振幅を求められる。

:::example[二つの曲率]
物体波が収束球面なら、その共役は発散球面である。二項を分離しなければ再生した像と双像が別の深さに現れる。
:::
:::example[標本化]
画素間隔 $5\,\mathrm{\mu m}$ の Nyquist 周波数は $100\,\mathrm{mm^{-1}}$。物体波の帯域半幅が $10\,\mathrm{mm^{-1}}$ で参照の搬送が $60\,\mathrm{mm^{-1}}$ なら、各交差項は中心から 50 から $70\,\mathrm{mm^{-1}}$ の範囲に入り、標本化帯域内に収まる。
:::

### 参照の強度と縞の大きさ
ある画素で物体と参照が同じ偏光の定常場なら、位相差を動かしたとき縞の可視度は $2|O||R|/(|O|^2+|R|^2)$。これは :eqref[eq-optics-two-beam] の固定振幅の場合である。参照強度を物体強度の 9 倍にすると可視度は $3/5$ であるが、交差項の絶対振幅は参照を強くすると大きくなる。露出はその交差項をセンサーが分解でき、かつ直流分が飽和しない範囲で選ぶ。

### 参照位相を順に変えて複素数を復元する
物体を動かさず $R e^{i\phi}$ の位相を $0,\pi/2,\pi,3\pi/2$ と変え、四枚の強度 $I_\phi$ を測る。$C=OR^*$ と置けば $I_\phi=|O|^2+|R|^2+2\Re(Ce^{-i\phi})$ である。従って
$$
C=\tfrac14[(I_0-I_\pi)+i(I_{\pi/2}-I_{3\pi/2})].
$$
定数項は差し引きで消え、二つの位相差が実部と虚部を与える。$R\ne0$ の画素で $C$ を $R^*$ で割れば $O$ が得られる。各露光で物体と参照の振幅が同じことが必要で、振動で位相が指定値からずれると共役項が残る。

:::example[四回露光で位相を求める]
同じ規格化で $I_0=9,I_\pi=1,I_{\pi/2}=I_{3\pi/2}=5$ とする。交差積 $OR^*=(9-1)/4=2$ で実正である。参照の振幅が実正で 1 なら物体振幅は 2。二強度の平均は 5 であり、$1+2^2$ と整合する。
:::

### オフ軸ホログラムの標本化と次数の分離

物体場の一次元スペクトルの支持を $[-B,B]$ とし、参照を $R e^{2\pi i f_c x}$ とする。記録強度の自己項 $|U|^2$ のスペクトルは $[-2B,2B]$ にある。これは積の Fourier 変換が畳み込みになるためである。交差項 $U R^*e^{-2\pi if_c x}$ とその共役の帯域はそれぞれ中心 $\mp f_c$、半幅 $B$。零次と交差項を分離する十分条件は $f_c>3B$。画素間隔 $p$ で折り返さず記録するには $f_c+B<1/(2p)$ も必要である。従って単純な周波数窓で復元する場合には、参照角を増やしすぎても小さすぎても失敗する。

:::example[側帯が分離する設定]
物体帯域を既知として参照角と画素間隔を選ぶ。

入力は物体場の帯域半幅 $B=20000\,\mathrm{m^{-1}}$、参照の搬送周波数 $f_c=80000\,\mathrm{m^{-1}}$、画素間隔 $p\simeq 4\times 10^{-6}\,\mathrm m$ である。

標本化周波数の半分を求める。

$$
\begin{aligned}
f_{\rm N}&=\frac{1}{2 p},\\
f_{\rm N}/(\mathrm{m^{-1}})&\simeq \frac{1}{2 \cdot \left(4\times 10^{-6}\right)}=125000.
\end{aligned}
$$

側帯の下端と零次帯域の上端を比較する。

$$
\begin{aligned}
\Delta_{\rm DC}&=- 3 B + f_c,\\
\Delta_{\rm DC}/(\mathrm{m^{-1}})&= - 3 \cdot \left(20000\right) + \left(80000\right)=20000.
\end{aligned}
$$

信号の最高周波数が折り返し端を越えないか調べる。

$$
\begin{aligned}
\Delta_{\rm N}&=- B + f_{\rm N} - f_c,\\
\Delta_{\rm N}/(\mathrm{m^{-1}})&= - \left(20000\right) + \left(125000\right) - \left(80000\right)=25000.
\end{aligned}
$$

参照と中心物体波が作る縞の幅を得る。

$$
\begin{aligned}
\Lambda_c&=\frac{1}{f_c},\\
\Lambda_c/(\mathrm m)&= \frac{1}{\left(80000\right)}\simeq 1.25\times 10^{-5}.
\end{aligned}
$$

搬送縞だけの標本数も比較する。ただし側帯の端には別の制約がある。

$$
\begin{aligned}
M&=\frac{\Lambda_c}{p},\\
M&\simeq \frac{\left(1.25\times 10^{-5}\right)}{\left(4\times 10^{-6}\right)}=\frac{25}{8}.
\end{aligned}
$$

二つの余裕が正なので一方の側帯を窓で選べる。
:::

:::example[参照角が小さい設定]
画素を変えず搬送周波数を下げる。

入力は物体場の帯域半幅 $B=20000\,\mathrm{m^{-1}}$、参照の搬送周波数 $f_c=50000\,\mathrm{m^{-1}}$、画素間隔 $p\simeq 4\times 10^{-6}\,\mathrm m$ である。

標本化周波数の半分を求める。

$$
\begin{aligned}
f_{\rm N}&=\frac{1}{2 p},\\
f_{\rm N}/(\mathrm{m^{-1}})&\simeq \frac{1}{2 \cdot \left(4\times 10^{-6}\right)}=125000.
\end{aligned}
$$

側帯の下端と零次帯域の上端を比較する。

$$
\begin{aligned}
\Delta_{\rm DC}&=- 3 B + f_c,\\
\Delta_{\rm DC}/(\mathrm{m^{-1}})&= - 3 \cdot \left(20000\right) + \left(50000\right)=-10000.
\end{aligned}
$$

信号の最高周波数が折り返し端を越えないか調べる。

$$
\begin{aligned}
\Delta_{\rm N}&=- B + f_{\rm N} - f_c,\\
\Delta_{\rm N}/(\mathrm{m^{-1}})&= - \left(20000\right) + \left(125000\right) - \left(50000\right)=55000.
\end{aligned}
$$

参照と中心物体波が作る縞の幅を得る。

$$
\begin{aligned}
\Lambda_c&=\frac{1}{f_c},\\
\Lambda_c/(\mathrm m)&= \frac{1}{\left(50000\right)}\simeq 2\times 10^{-5}.
\end{aligned}
$$

搬送縞だけの標本数も比較する。ただし側帯の端には別の制約がある。

$$
\begin{aligned}
M&=\frac{\Lambda_c}{p},\\
M&\simeq \frac{\left(2\times 10^{-5}\right)}{\left(4\times 10^{-6}\right)}=5.
\end{aligned}
$$

標本化は足りるが、無条件に側帯を切り離すことはできない。
:::

:::example[参照角が大きい設定]
逆に搬送周波数を大きくする。

入力は物体場の帯域半幅 $B=20000\,\mathrm{m^{-1}}$、参照の搬送周波数 $f_c=115000\,\mathrm{m^{-1}}$、画素間隔 $p\simeq 4\times 10^{-6}\,\mathrm m$ である。

標本化周波数の半分を求める。

$$
\begin{aligned}
f_{\rm N}&=\frac{1}{2 p},\\
f_{\rm N}/(\mathrm{m^{-1}})&\simeq \frac{1}{2 \cdot \left(4\times 10^{-6}\right)}=125000.
\end{aligned}
$$

側帯の下端と零次帯域の上端を比較する。

$$
\begin{aligned}
\Delta_{\rm DC}&=- 3 B + f_c,\\
\Delta_{\rm DC}/(\mathrm{m^{-1}})&= - 3 \cdot \left(20000\right) + \left(115000\right)=55000.
\end{aligned}
$$

信号の最高周波数が折り返し端を越えないか調べる。

$$
\begin{aligned}
\Delta_{\rm N}&=- B + f_{\rm N} - f_c,\\
\Delta_{\rm N}/(\mathrm{m^{-1}})&= - \left(20000\right) + \left(125000\right) - \left(115000\right)=-10000.
\end{aligned}
$$

参照と中心物体波が作る縞の幅を得る。

$$
\begin{aligned}
\Lambda_c&=\frac{1}{f_c},\\
\Lambda_c/(\mathrm m)&= \frac{1}{\left(115000\right)}\simeq 8.6957\times 10^{-6}.
\end{aligned}
$$

搬送縞だけの標本数も比較する。ただし側帯の端には別の制約がある。

$$
\begin{aligned}
M&=\frac{\Lambda_c}{p},\\
M&\simeq \frac{\left(8.695652174\times 10^{-6}\right)}{\left(4\times 10^{-6}\right)}=\frac{50}{23}.
\end{aligned}
$$

搬送縞自体を記録できても、側帯の右端は折り返す。
:::

## 演習問題

::::exercise{level=basic}
$\lambda=500\,\mathrm{nm}$ の伝播成分の最大横周波数を求めよ。

:::solution
式 :eqref[eq-optics-angular-spectrum] で平方根の引数が非負となる上限は $1/\lambda=2\times10^6\,\mathrm{m^{-1}}$。
:::
::::

::::exercise{level=basic}
$\lambda=500\,\mathrm{nm},F=10\,\mathrm{cm}$ で $f=2\times10^4\,\mathrm{m^{-1}}$ を焦点面座標に変換せよ。

:::solution
比例式 $x=\lambda Ff$ に入れ $1.0\,\mathrm{mm}$ を得る。
:::
::::

::::exercise{level=basic}
一様な平面波照明の下で、開口の横移動によって Fraunhofer 強度は変わるか。

:::solution
変数変換で複素振幅は絶対値 1 の位相因子だけを受ける。二乗で消えるので強度分布は変わらない。
:::
::::

::::exercise{level=basic}
弱位相試料で位相板を $-\pi/2$ にすると一次の像強度はどうなるか。

:::solution
振幅 $-i+i\phi=-i(1-\phi)$ の二乗絶対値で、一次まで $1-2\phi$ である。
:::
::::

::::exercise{level=advanced}
$0<\int|P|^2\dd^2s<\infty$ を満たす瞳に対し $H_{\mathrm{inc}}(0)=1$、$\abs{H_{\mathrm{inc}}}\le1$ を示せ。

:::solution
零では分子と分母が等しい。一般には瞳の二つの平行移動に Cauchy–Schwarz を適用する。両ノルムが元の瞳と一致するので、分子の絶対値は分母以下となる。
:::
::::

::::exercise{level=advanced}
一様で無収差な円瞳の半径が $50\,\mathrm{mm^{-1}}$ のとき強度伝達が初めてゼロとなる正の周波数を求めよ。

:::solution
重なり面積が消える距離は直径なので $100\,\mathrm{mm^{-1}}$ となる。
:::
::::

::::exercise{level=advanced}
物体波の周波数帯が $[-B,B]$ で参照の搬送周波数が $f_R>0$ のとき、支持の上限だけを使って直流帯域と正側交差項の分離を保証する条件を求めよ。

:::solution
物体強度の帯域は差周波数より $[-2B,2B]$。交差項は $[f_R-B,f_R+B]$ なので $f_R-B>2B$、従って $f_R>3B$ が十分である。
:::
::::

::::exercise{level=advanced}
$f=2/\lambda$ の成分が $z=\lambda/(2\pi)$ 進む振幅比を求めよ。

:::solution
$k_z=ik\sqrt3$ で $kz=1$ なので、振幅比は $e^{-\sqrt3}$ となる。
:::
::::

::::exercise{level=advanced}
オフ軸ホログラムで物体場の周波数支持が $[-B,B]$、画素間隔が $5\,\mathrm{\mu m}$ とする。零次と正側の側帯を重なりなく分け、折返しも避ける搬送周波数が存在するための $B$ の条件を求めよ。

:::solution
必要な二つの余裕は $f_c-3B>0$ と $f_{\rm N}-(f_c+B)>0$。従って $3B<f_c<f_{\rm N}-B$ が空でない条件は $4B<f_{\rm N}$ である。画素間隔から $f_{\rm N}=1/(2p)=10^5\,\mathrm{m^{-1}}$。従って $B<2.5\times10^4\,\mathrm{m^{-1}}$ なら二つの条件を同時に満たす搬送周波数を選べる。等号では帯域同士が接するので、正の余裕は取れない。
:::
::::

## 確認テスト

::::quiz
コヒーレント結像で線形に写す量を選べ。

- [x] 複素振幅
- [ ] 強度
- [ ] 位相だけ

:::explanation
式 :eqref[eq-optics-coherent-imaging] の畳み込みの入力は場である。検出強度は、その出力の絶対値を二乗した結果である。
:::
::::

::::quiz{answer="40" tol="1e-6"}
$\lambda=500\,\mathrm{nm},F=0.10\,\mathrm m$ で $x=2\,\mathrm{mm}$ に対応する周波数を $\mathrm{mm^{-1}}$ で答えよ。

:::explanation
$x/(\lambda F)=4.0\times10^4\,\mathrm{m^{-1}}$ であり、mm に換算して 40 である。
:::
::::

::::quiz
強度結像の光学伝達関数に現れる相関を選べ。

- [x] 複素瞳の自己相関
- [ ] 物体強度の自己相関
- [ ] 時間の電場自己相関

:::explanation
点像強度が $\abs{h}^2$ である。これを Fourier 変換すると、二つの瞳の周波数差を固定した積分となる。
:::
::::

::::quiz{answer="2" tol="1e-6"}
一様円瞳の強度結像の遮断周波数は振幅結像の何倍か。

:::explanation
移動した二円の重なりが消えるのは移動量が半径の二倍のときで、振幅の上限は半径だからである。
:::
::::

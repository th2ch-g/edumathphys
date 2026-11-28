---
title: コヒーレンス
description: 光の相関関数を時間と空間で定義し、スペクトル、天体の角径、強度相関との関係を導く。
notation:
  - symbol: '\Gamma_{12}(\tau)'
    meaning: 二点の相互コヒーレンス関数
  - symbol: '\gamma_{12}(\tau)'
    meaning: 規格化した相互コヒーレンス
  - symbol: 'S(\nu)'
    meaning: 正周波数側のスペクトル密度
  - symbol: '\boldsymbol B'
    meaning: 受光点間の基線
  - symbol: 'g^{(2)}'
    meaning: 規格化強度相関
---

## この章の目標

- 実電場から解析信号と包絡線を定義する。
- 時間相関から線幅と縞の見え方を計算する。
- 遠方の光源分布を基線上の複素相関と結ぶ。
- 熱光とコヒーレント光の強度相関を区別する。

## 準単色光と解析信号

### 実信号から正周波数を選ぶ
実電場 $E(t)$ の Fourier 変換を $\widehat E(\nu)=\int E(t)e^{2\pi i\nu t}\dd t$ とする。実数値なので $\widehat E(-\nu)=\widehat E(\nu)^*$ である。正周波数を二倍して逆変換した $u(t)=2\int_0^\infty\widehat E(\nu)e^{-2\pi i\nu t}\dd\nu$ を、この時間規約の :term[解析信号]{en="analytic signal" yomi="かいせきしんごう"} と呼ぶ。直流成分は別に扱う。負周波数の積分は正周波数の共役なので $E=\Re u$ となる。Fourier 変換の反転については :ref[thm-four-euclid-inversion] を参照する。定常ランダム場には全時間の Fourier 積分をそのまま使えない場合がある。その場合は後述するように、有限時間のスペクトルの二次モーメントからスペクトルを求める。

正周波数の帯域が $\nu_0$ の近くで相対的に狭い光を :term[準単色光]{en="quasi-monochromatic light" yomi="じゅんたんしょくこう"} という。$u(t)=A(t)e^{-2\pi i\nu_0t}$ と置くと、$A$ のスペクトルは離調だけを含むので光周期に比べてゆっくり変わる。実電場を二乗すると $|u|^2/2$ と $\pm2\nu_0$ 付近の振動項に分かれ、光周期より長い平均で振動項は消える。以後は媒質による一定の換算係数を振幅へ吸収し、強度を $I(t)=|u(t)|^2$ と書く。

:::example[単一周波数]
$E(t)=E_0\cos(2\pi\nu_0t+\phi)$ の解析信号は $u=E_0e^{-i(2\pi\nu_0t+\phi)}$。負周波数も別の自由度として足すと振幅を二重に数える。
:::
:::example[二周波の包絡線]
$\nu_0\pm\Delta\nu/2$ の等振幅波を足すと $A=2a\cos(\pi\Delta\nu t)$。強度の周期は $1/\Delta\nu$ で、搬送周期 $1/\nu_0$ より長い。
:::

### 解析信号の虚部を実信号から計算する
減衰の十分速い実信号について、Hilbert 変換を $\mathcal HE=\frac1\pi\operatorname{PV}\int E(t')/(t-t')\dd t'$ と定義する。この変換は $e^{-i\omega t}$ に対して $i\operatorname{sgn}\omega$ を掛ける。これは奇核の Fourier 積分を正負半軸に分け、$\int_0^\infty e^{-\epsilon x}\sin(\omega x)/x\dd x=\arctan(\omega/\epsilon)$ を使えば得られる。最後の積分は $\omega$ で微分すると $\epsilon/(\epsilon^2+\omega^2)$ となり、$\omega=0$ でゼロである。従って $u=E-i\mathcal HE$ にすると、正周波数の係数が二倍、負周波数の係数がゼロになる。この符号は用いる指数の規約に依存する。

:::example[一つの光周期と包絡線の周期]
中心 $5\times10^{14}\,\mathrm{Hz}$ に対して間隔 $10^9\,\mathrm{Hz}$ の二色を重ねる。搬送周期は $2\,\mathrm{fs}$、強度のうなり周期は $1\,\mathrm{ns}$ である。その間に $5\times10^5$ 周期の搬送波が入るので、光周期を平均しても強度の変化が残る時間窓を設定できる。
:::

### 二つの鋭い線を搬送波から分離する

実電場が振幅 $a,b$ の二つの余弦の和で、正周波数を $\nu_0\pm\delta$ とする。各余弦から正周波数だけを取って係数を二倍すると、ここでは搬送波を除いた包絡線を $u(t)$ と書くと、$u(t)=a e^{-2\pi i\delta t}+b e^{2\pi i\delta t}$ になる。従って $|u|^2=a^2+b^2+2ab\cos(4\pi\delta t)$。包絡線の位相は二成分の干渉で時間に依存し、包絡線がゼロになった時刻では位相も不定になる。包絡線を低周波とみなすには $\delta\ll\nu_0$ が要るが、この条件で包絡線が定数になるわけではない。

:::example[等しい二線の明点]
二線が同相になる時刻を見る。

入力は上側の振幅 $a=1$、下側の振幅 $b=1$、搬送波からの線のずれ $\delta=1000000\,\mathrm{Hz}$、観測時刻 $t=0\,\mathrm s$ である。

$\nu_0$ の位相を取り除いて二線の周波数差を見る。

$$
\begin{aligned}
f_{\rm beat}&=2 \delta,\\
f_{\rm beat}/(\mathrm{Hz})&= 2 \cdot \left(1000000\right)=2000000.
\end{aligned}
$$

周波数の逆数で包絡線強度の繰返しを得る。

$$
\begin{aligned}
T_{\rm beat}&=\frac{1}{f_{\rm beat}},\\
T_{\rm beat}/(\mathrm s)&= \frac{1}{\left(2000000\right)}\simeq 5\times 10^{-7}.
\end{aligned}
$$

差周波数を観測時刻まで積分する。

$$
\begin{aligned}
\phi&=2 \pi f_{\rm beat} t,\\
\phi/(\mathrm{rad})&= 2 \cdot \pi \cdot \left(2000000\right) \cdot \left(0\right)=0.
\end{aligned}
$$

正周波数を抽出した後も二成分の干渉は残る。

$$
\begin{aligned}
\lvert u\rvert^2&=a^{2} + 2 a b \cos{\left(\phi \right)} + b^{2},\\
\lvert u\rvert^2&= \left(1\right)^{2} + 2 \cdot \left(1\right) \cdot \left(1\right) \cdot \cos{\left(\left(0\right) \right)} + \left(1\right)^{2}=4.
\end{aligned}
$$

交差項の振幅を平均の強度で割る。

$$
\begin{aligned}
V&=\frac{2 a b}{a^{2} + b^{2}},\\
V&= \frac{2 \cdot \left(1\right) \cdot \left(1\right)}{\left(1\right)^{2} + \left(1\right)^{2}}=1.
\end{aligned}
$$

強度が単独二線の強度和の二倍になる。
:::

:::example[等しい二線の暗点]
同じ光を半うなり周期だけ後で見る。

入力は上側の振幅 $a=1$、下側の振幅 $b=1$、搬送波からの線のずれ $\delta=1000000\,\mathrm{Hz}$、観測時刻 $t\simeq 2.5\times 10^{-7}\,\mathrm s$ である。

$\nu_0$ の位相を取り除いて二線の周波数差を見る。

$$
\begin{aligned}
f_{\rm beat}&=2 \delta,\\
f_{\rm beat}/(\mathrm{Hz})&= 2 \cdot \left(1000000\right)=2000000.
\end{aligned}
$$

周波数の逆数で包絡線強度の繰返しを得る。

$$
\begin{aligned}
T_{\rm beat}&=\frac{1}{f_{\rm beat}},\\
T_{\rm beat}/(\mathrm s)&= \frac{1}{\left(2000000\right)}\simeq 5\times 10^{-7}.
\end{aligned}
$$

差周波数を観測時刻まで積分する。

$$
\begin{aligned}
\phi&=2 \pi f_{\rm beat} t,\\
\phi/(\mathrm{rad})&\simeq 2 \cdot \pi \cdot \left(2000000\right) \cdot \left(2.5\times 10^{-7}\right)\simeq 3.1416.
\end{aligned}
$$

正周波数を抽出した後も二成分の干渉は残る。

$$
\begin{aligned}
\lvert u\rvert^2&=a^{2} + 2 a b \cos{\left(\phi \right)} + b^{2},\\
\lvert u\rvert^2&\simeq \left(1\right)^{2} + 2 \cdot \left(1\right) \cdot \left(1\right) \cdot \cos{\left(\left(3.141592654\right) \right)} + \left(1\right)^{2}=0.
\end{aligned}
$$

交差項の振幅を平均の強度で割る。

$$
\begin{aligned}
V&=\frac{2 a b}{a^{2} + b^{2}},\\
V&= \frac{2 \cdot \left(1\right) \cdot \left(1\right)}{\left(1\right)^{2} + \left(1\right)^{2}}=1.
\end{aligned}
$$

解析信号はそこでゼロになる。搬送波の位相から独立な包絡線の零点である。
:::

:::example[振幅が釣り合わない二線]
相対位相は逆相のまま一方の振幅を下げる。

入力は上側の振幅 $a=1$、下側の振幅 $b=\frac{1}{2}$、搬送波からの線のずれ $\delta=1000000\,\mathrm{Hz}$、観測時刻 $t\simeq 2.5\times 10^{-7}\,\mathrm s$ である。

$\nu_0$ の位相を取り除いて二線の周波数差を見る。

$$
\begin{aligned}
f_{\rm beat}&=2 \delta,\\
f_{\rm beat}/(\mathrm{Hz})&= 2 \cdot \left(1000000\right)=2000000.
\end{aligned}
$$

周波数の逆数で包絡線強度の繰返しを得る。

$$
\begin{aligned}
T_{\rm beat}&=\frac{1}{f_{\rm beat}},\\
T_{\rm beat}/(\mathrm s)&= \frac{1}{\left(2000000\right)}\simeq 5\times 10^{-7}.
\end{aligned}
$$

差周波数を観測時刻まで積分する。

$$
\begin{aligned}
\phi&=2 \pi f_{\rm beat} t,\\
\phi/(\mathrm{rad})&\simeq 2 \cdot \pi \cdot \left(2000000\right) \cdot \left(2.5\times 10^{-7}\right)\simeq 3.1416.
\end{aligned}
$$

正周波数を抽出した後も二成分の干渉は残る。

$$
\begin{aligned}
\lvert u\rvert^2&=a^{2} + 2 a b \cos{\left(\phi \right)} + b^{2},\\
\lvert u\rvert^2&\simeq \left(1\right)^{2} + 2 \cdot \left(1\right) \cdot \left(\frac{1}{2}\right) \cdot \cos{\left(\left(3.141592654\right) \right)} + \left(\frac{1}{2}\right)^{2}=\frac{1}{4}.
\end{aligned}
$$

交差項の振幅を平均の強度で割る。

$$
\begin{aligned}
V&=\frac{2 a b}{a^{2} + b^{2}},\\
V&= \frac{2 \cdot \left(1\right) \cdot \left(\frac{1}{2}\right)}{\left(1\right)^{2} + \left(\frac{1}{2}\right)^{2}}=\frac{4}{5}.
\end{aligned}
$$

暗い点にも有限の強度が残る。
:::

## 時間的コヒーレンスと自己相関関数

### 遅延した同じ光を重ねる
平均値が一定で二時刻の二次相関が時刻差だけで決まる二乗可積分な確率場を :term[広義定常]{en="wide-sense stationarity" yomi="こうぎていじょう"} と呼ぶ。以下の $\langle\cdot\rangle$ は集団平均を表す。時間平均一回でこれを置き換えるときは、対象とするモーメントについてエルゴード性などを確認する。:term[時間的コヒーレンス]{en="temporal coherence" yomi="じかんてきこひーれんす"} は、同一点で時間差 $\tau$ の場について
$$
\Gamma(\tau)=\langle u^*(t)u(t+\tau)\rangle,\qquad
g^{(1)}(\tau)=\Gamma(\tau)/\Gamma(0)
\label{eq-optics-time-correlation}
$$
で測る。順序を交換して定常性を適用すると $\Gamma(-\tau)=\Gamma(\tau)^*$、Cauchy–Schwarz から $|g^{(1)}|\le1$。干渉計で同じ場を遅延させた二波の強度を釣り合わせると、式 :eqref[eq-optics-two-beam] で可視度が $|g^{(1)}(\tau)|$ になる。

:term[コヒーレンス時間]{en="coherence time" yomi="こひーれんすじかん"} の定義には半値幅や積分幅がある。本章では特に指定がなければ $\tau_c=\int_{-\infty}^\infty |g^{(1)}(\tau)|^2\dd\tau$ とする。真空中のコヒーレンス長は $c\tau_c$ である。装置の腕長と比較する際は、Michelson では片腕の変位が往復で二倍の光路差を作る。

:::example[完全な単色波]
$\Gamma(\tau)=I_0 e^{-2\pi i\nu_0\tau}$ なら $\abs{g^{(1)}}=1$。定義上 $\tau_c$ は無限大であり、鏡を動かすと位相だけが変わる。
:::
:::example[指数相関]
$g^{(1)}(\tau)=e^{-a|\tau|}e^{-2\pi i\nu_0\tau}$ とすると $\tau_c=2\int_0^\infty e^{-2a\tau}\dd\tau=1/a$。振幅の片側 $1/e$ 幅もこの例では同じ値である。
:::

:::example[光路差を遅延へ換算する]
真空で遅延を $10\,\mathrm{ps}$ 変えるには光路差を約 $3.0\,\mathrm{mm}$ 変える。Michelson の一枚の鏡の位置としてはその半分の $1.5\,\mathrm{mm}$ である。可視度の包絡線を測るときは、この遅い尺度と光周期の縞を区別する。
:::

### コヒーレンス時間の複数の定義を換算する
相関の包絡線を $|g^{(1)}(\tau)|=e^{-a|\tau|}$ とすると、片側の振幅半値点は $\ln2/a$、片側の強度半値点は $\ln2/(2a)$、両側二乗積分幅は $1/a$ である。同じ光源を測っても、可視度の半値幅で報告するか、その二乗の積分で報告するかで波形に依存する定数が変わる。スペクトル側の半値幅 $\Delta\nu=a/\pi$ に対して本章の $\tau_c\Delta\nu=1/\pi$ になる。一方、幅 $\Delta\nu$ の矩形スペクトルでは Parseval の等式で $\tau_c=1/\Delta\nu$ となる。時間と帯域の積を比較する際は、スペクトルの形と幅の定義を一緒に記す。

### 時間コヒーレンスと腕の強度比を分離する

二腕で共通の光源を観測し、片腕のパワーを一定比 $q$ にしたとする。光路差を固定して位相だけを掃引すると、強度の平均は $(1+q)I_0$、変動の振幅は $2\sqrt q I_0|g^{(1)}(\tau)|$。従って可視度から時間相関を取り出すには $2\sqrt q/(1+q)$ を割り戻す。Michelson では鏡の変位 $x$ と時間差は $2x/c$ であり、鏡の移動と光路差を区別する。狭帯域 Lorentz スペクトルの半値全幅が $w$ なら、包絡線は $e^{-\pi w|\tau|}$。これは $(w/2)^2+(\nu-\nu_0)^2$ の逆数を Fourier 積分した結果である。

:::example[腕の強度が釣り合う条件]
遅れによって可視度が下がる割合を測る。

入力は腕の強度比 $q=1$、片鏡の変位 $x=\frac{1}{100}\,\mathrm m$、スペクトル半値全幅 $w=1000000000\,\mathrm{Hz}$ である。

鏡の変位を往復の光路差へ換算する。

$$
\begin{aligned}
\tau&=\frac{2x}{c},\\
\tau/(\mathrm s)&= \frac{\left(\frac{1}{100}\right)}{149896229}\simeq 6.6713\times 10^{-11}.
\end{aligned}
$$

Lorentz 線形状の Fourier 変換を使う。

$$
\begin{aligned}
G&=e^{- \pi w \left|{\tau}\right|},\\
G&\simeq e^{- \pi \cdot \left(1000000000\right) \cdot \left|{\left(6.671281904\times 10^{-11}\right)}\right|}\simeq 0.81092.
\end{aligned}
$$

時間遅れがゼロでも強度が釣り合わなければ一にはならない。

$$
\begin{aligned}
\mathcal V_0&=\frac{2 \sqrt{q}}{q + 1},\\
\mathcal V_0&= \frac{2 \cdot \sqrt{\left(1\right)}}{\left(1\right) + 1}=1.
\end{aligned}
$$

二つの原因を掛けて測定値と比較する。

$$
\begin{aligned}
\mathcal V&=G \mathcal V_0,\\
\mathcal V&\simeq \left(0.8109211119\right) \cdot \left(1\right)\simeq 0.81092.
\end{aligned}
$$

二乗の指数の積分を往復変位でなく光路長へ換算する。

$$
\begin{aligned}
L_c&=\frac{c}{\pi w},\\
L_c/(\mathrm m)&= \frac{299792458}{\pi \cdot \left(1000000000\right)}\simeq 0.095427.
\end{aligned}
$$

可視度をそのまま時間相関に使える。
:::

:::example[腕に損失がある条件]
同じ遅れと線幅で一腕だけ減光する。

入力は腕の強度比 $q=\frac{1}{4}$、片鏡の変位 $x=\frac{1}{100}\,\mathrm m$、スペクトル半値全幅 $w=1000000000\,\mathrm{Hz}$ である。

鏡の変位を往復の光路差へ換算する。

$$
\begin{aligned}
\tau&=\frac{2x}{c},\\
\tau/(\mathrm s)&= \frac{\left(\frac{1}{100}\right)}{149896229}\simeq 6.6713\times 10^{-11}.
\end{aligned}
$$

Lorentz 線形状の Fourier 変換を使う。

$$
\begin{aligned}
G&=e^{- \pi w \left|{\tau}\right|},\\
G&\simeq e^{- \pi \cdot \left(1000000000\right) \cdot \left|{\left(6.671281904\times 10^{-11}\right)}\right|}\simeq 0.81092.
\end{aligned}
$$

時間遅れがゼロでも強度が釣り合わなければ一にはならない。

$$
\begin{aligned}
\mathcal V_0&=\frac{2 \sqrt{q}}{q + 1},\\
\mathcal V_0&= \frac{2 \cdot \sqrt{\left(\frac{1}{4}\right)}}{\left(\frac{1}{4}\right) + 1}=\frac{4}{5}.
\end{aligned}
$$

二つの原因を掛けて測定値と比較する。

$$
\begin{aligned}
\mathcal V&=G \mathcal V_0,\\
\mathcal V&\simeq \left(0.8109211119\right) \cdot \left(\frac{4}{5}\right)\simeq 0.64874.
\end{aligned}
$$

二乗の指数の積分を往復変位でなく光路長へ換算する。

$$
\begin{aligned}
L_c&=\frac{c}{\pi w},\\
L_c/(\mathrm m)&= \frac{299792458}{\pi \cdot \left(1000000000\right)}\simeq 0.095427.
\end{aligned}
$$

時間コヒーレンスは前例と同じで、目に見える縞は弱くなる。
:::

:::example[より広い線幅を測る]
腕を再び釣り合わせて光源だけを変える。

入力は腕の強度比 $q=1$、片鏡の変位 $x=\frac{1}{100}\,\mathrm m$、スペクトル半値全幅 $w\simeq 1\times 10^{10}\,\mathrm{Hz}$ である。

鏡の変位を往復の光路差へ換算する。

$$
\begin{aligned}
\tau&=\frac{2x}{c},\\
\tau/(\mathrm s)&= \frac{\left(\frac{1}{100}\right)}{149896229}\simeq 6.6713\times 10^{-11}.
\end{aligned}
$$

Lorentz 線形状の Fourier 変換を使う。

$$
\begin{aligned}
G&=e^{- \pi w \left|{\tau}\right|},\\
G&\simeq e^{- \pi \cdot \left(1\times 10^{10}\right) \cdot \left|{\left(6.671281904\times 10^{-11}\right)}\right|}\simeq 0.12297.
\end{aligned}
$$

時間遅れがゼロでも強度が釣り合わなければ一にはならない。

$$
\begin{aligned}
\mathcal V_0&=\frac{2 \sqrt{q}}{q + 1},\\
\mathcal V_0&= \frac{2 \cdot \sqrt{\left(1\right)}}{\left(1\right) + 1}=1.
\end{aligned}
$$

二つの原因を掛けて測定値と比較する。

$$
\begin{aligned}
\mathcal V&=G \mathcal V_0,\\
\mathcal V&\simeq \left(0.1229662905\right) \cdot \left(1\right)\simeq 0.12297.
\end{aligned}
$$

二乗の指数の積分を往復変位でなく光路長へ換算する。

$$
\begin{aligned}
L_c&=\frac{c}{\pi w},\\
L_c/(\mathrm m)&\simeq \frac{299792458}{\pi \cdot \left(1\times 10^{10}\right)}\simeq 0.0095427.
\end{aligned}
$$

同じ鏡の変位で時間相関が小さくなる。
:::

## Wiener–Khinchin の定理とスペクトル

### 有限区間のスペクトルからの導出
区間 $[-T/2,T/2]$ 上で $\widehat u_T(\nu)=\int u(t)e^{2\pi i\nu t}\dd t$ とする。絶対値を二乗して平均し、二時刻の差を $\tau$ と置く。二時刻が区間に収まる長さは $T-|\tau|$ なので
$$
\frac{\langle|\widehat u_T(\nu)|^2\rangle}{T}=
\int_{-T}^T(1-|\tau|/T)\Gamma(\tau)e^{2\pi i\nu\tau}\dd\tau .
$$
$\Gamma$ が可積分なら支配収束定理により極限は $S(\nu)=\int\Gamma(\tau)e^{2\pi i\nu\tau}\dd\tau$。左辺は非負なので $S\ge0$ である。$S$ も可積分で $\Gamma$ が連続なら Fourier 反転を使って
$$
\Gamma(\tau)=\int S(\nu)e^{-2\pi i\nu\tau}\dd\nu .
\label{eq-optics-wiener-khinchin}
$$
これが :term[Wiener–Khinchin の定理]{en="Wiener–Khinchin theorem" yomi="うぃーなーひんちんのていり"} の可積分な場合である。準単色の解析信号では正周波数側に分布する。理想単色波のように $\Gamma$ が可積分でない場合は $S$ をスペクトル測度として扱い、単色波は強度を重みに持つ点質量となる。離調を全実軸に延ばす狭帯域模型では負の絶対周波数への尾を無視する。

:::example[矩形帯域]
幅 $\Delta\nu$ の一様スペクトルを中心 $\nu_0$ に置く。積分で $g^{(1)}=e^{-2\pi i\nu_0\tau}\operatorname{sinc}(\pi\Delta\nu\tau)$。可視度の最初の零点は $|\tau|=1/\Delta\nu$ にある。
:::
:::example[指数相関の線形状]
前節の相関では、正負の $\tau$ に分けて積分すると $S/I_0=2a/[a^2+4\pi^2(\nu-\nu_0)^2]$。半値の離調は $a/(2\pi)$、周波数の全半値幅は $a/\pi$ である。
:::

:::example[Gaussian スペクトル]
離調の標準偏差を $\sigma_\nu$ とする規格化 Gaussian スペクトルは、二次指数の積分で $g^{(1)}(\tau)=e^{-2\pi i\nu_0\tau}e^{-2\pi^2\sigma_\nu^2\tau^2}$ を与える。両側の二乗積分幅は $1/(2\sqrt\pi\,\sigma_\nu)$。周波数の標準偏差を二倍にするとこの時間は半分になり、矩形スペクトルと違って相関に有限の零点はない。
:::

### 相関が満たす正値条件
任意の複素係数 $c_1,\ldots,c_M$ と時刻 $t_1,\ldots,t_M$ について、$\sum_{jk} c_j^* \Gamma(t_k-t_j)c_k=\langle|\sum_k c_k u(t_k)|^2\rangle\ge0$ である。従って測定した $\Gamma$ が Hermite 対称でも、すべての時刻の組に対してこの条件を満たさなければ定常場の相関にはならない。反対に非負で可積分な $S$ から :eqref[eq-optics-wiener-khinchin] を作ると、二次形式は $\int S(\nu)|\sum_k c_k e^{-2\pi i\nu t_k}|^2\dd\nu\ge0$。これにより実験で仮定する線形状と時間相関が矛盾しないかを確認できる。

### 二つのピークを持つスペクトルの相関

離調について同じ標準偏差 $\sigma$ の Gaussian 線二つを $\pm\Delta/2$ に置き、それぞれの積分強度を半分にする。スペクトルの並進が変換後の指数になるので、規格化相関から搬送位相を除いた値は $\exp(-2\pi^2\sigma^2\tau^2)\cos(\pi\Delta\tau)$。零点は二本の線の間隔で決まり、長い遅延における振幅の上限は個々の線幅で決まる。二乗積分では $\cos^2 x=(1+\cos2x)/2$ を使い、二項の Gaussian 積分から $\tau_c=[1+e^{-\Delta^2/(4\sigma^2)}]/(4\sqrt\pi\sigma)$ を得る。線幅が同じでも、二本の混合はコヒーレンス時間を変える。

:::example[線が重なる場合]
二つの中心が一致したとする。

入力は各線の標準偏差 $\sigma=100000000\,\mathrm{Hz}$、二線の間隔 $\Delta=0\,\mathrm{Hz}$、時間差 $\tau\simeq 2.5\times 10^{-10}\,\mathrm s$ である。

一本の線幅で決まる減衰を入れる。

$$
\begin{aligned}
A&=e^{- 2 \pi^{2} \sigma^{2} \tau^{2}},\\
A&\simeq e^{- 2 \cdot \pi^{2} \cdot \left(100000000\right)^{2} \cdot \left(\left(2.5\times 10^{-10}\right)\right)^{2}}\simeq 0.98774.
\end{aligned}
$$

二つの中心周波数の位相因子を等しく足す。

$$
\begin{aligned}
B&=\cos{\left(\pi \Delta \tau \right)},\\
B&\simeq \cos{\left(\pi \cdot \left(0\right) \cdot \left(2.5\times 10^{-10}\right) \right)}=1.
\end{aligned}
$$

搬送波の位相を除いた相関の絶対値を取る。

$$
\begin{aligned}
\mathcal V&=\left|{A B}\right|,\\
\mathcal V&\simeq \left|{\left(0.9877387834\right) \cdot \left(1\right)}\right|\simeq 0.98774.
\end{aligned}
$$

二つの線の重なりによる補正を含める。

$$
\begin{aligned}
\tau_c&=\frac{1 + e^{- \frac{\Delta^{2}}{4 \sigma^{2}}}}{4 \sqrt{\pi} \sigma},\\
\tau_c/(\mathrm s)&= \frac{1 + e^{- \frac{\left(0\right)^{2}}{4 \cdot \left(100000000\right)^{2}}}}{4 \cdot \sqrt{\pi} \cdot \left(100000000\right)}\simeq 2.8209\times 10^{-9}.
\end{aligned}
$$

単一の Gaussian 線と同じ積分幅を得る。
:::

:::example[分離した線の最初の零点]
二線の間隔を線幅より大きくし、包絡線の零点へ遅延を調整する。

入力は各線の標準偏差 $\sigma=100000000\,\mathrm{Hz}$、二線の間隔 $\Delta=2000000000\,\mathrm{Hz}$、時間差 $\tau\simeq 2.5\times 10^{-10}\,\mathrm s$ である。

一本の線幅で決まる減衰を入れる。

$$
\begin{aligned}
A&=e^{- 2 \pi^{2} \sigma^{2} \tau^{2}},\\
A&\simeq e^{- 2 \cdot \pi^{2} \cdot \left(100000000\right)^{2} \cdot \left(\left(2.5\times 10^{-10}\right)\right)^{2}}\simeq 0.98774.
\end{aligned}
$$

二つの中心周波数の位相因子を等しく足す。

$$
\begin{aligned}
B&=\cos{\left(\pi \Delta \tau \right)},\\
B&\simeq \cos{\left(\pi \cdot \left(2000000000\right) \cdot \left(2.5\times 10^{-10}\right) \right)}=0.
\end{aligned}
$$

搬送波の位相を除いた相関の絶対値を取る。

$$
\begin{aligned}
\mathcal V&=\left|{A B}\right|,\\
\mathcal V&\simeq \left|{\left(0.9877387834\right) \cdot \left(0\right)}\right|=0.
\end{aligned}
$$

二つの線の重なりによる補正を含める。

$$
\begin{aligned}
\tau_c&=\frac{1 + e^{- \frac{\Delta^{2}}{4 \sigma^{2}}}}{4 \sqrt{\pi} \sigma},\\
\tau_c/(\mathrm s)&= \frac{1 + e^{- \frac{\left(2000000000\right)^{2}}{4 \cdot \left(100000000\right)^{2}}}}{4 \cdot \sqrt{\pi} \cdot \left(100000000\right)}\simeq 1.4105\times 10^{-9}.
\end{aligned}
$$

一本ずつの相関は大きくても、二線の寄与が逆位相で消える。
:::

:::example[零点より後で相関が戻る]
同じスペクトルで遅延を二倍にする。

入力は各線の標準偏差 $\sigma=100000000\,\mathrm{Hz}$、二線の間隔 $\Delta=2000000000\,\mathrm{Hz}$、時間差 $\tau\simeq 5\times 10^{-10}\,\mathrm s$ である。

一本の線幅で決まる減衰を入れる。

$$
\begin{aligned}
A&=e^{- 2 \pi^{2} \sigma^{2} \tau^{2}},\\
A&\simeq e^{- 2 \cdot \pi^{2} \cdot \left(100000000\right)^{2} \cdot \left(\left(5\times 10^{-10}\right)\right)^{2}}\simeq 0.95185.
\end{aligned}
$$

二つの中心周波数の位相因子を等しく足す。

$$
\begin{aligned}
B&=\cos{\left(\pi \Delta \tau \right)},\\
B&\simeq \cos{\left(\pi \cdot \left(2000000000\right) \cdot \left(5\times 10^{-10}\right) \right)}=-1.
\end{aligned}
$$

搬送波の位相を除いた相関の絶対値を取る。

$$
\begin{aligned}
\mathcal V&=\left|{A B}\right|,\\
\mathcal V&\simeq \left|{\left(0.9518498074\right) \cdot \left(-1\right)}\right|\simeq 0.95185.
\end{aligned}
$$

二つの線の重なりによる補正を含める。

$$
\begin{aligned}
\tau_c&=\frac{1 + e^{- \frac{\Delta^{2}}{4 \sigma^{2}}}}{4 \sqrt{\pi} \sigma},\\
\tau_c/(\mathrm s)&= \frac{1 + e^{- \frac{\left(2000000000\right)^{2}}{4 \cdot \left(100000000\right)^{2}}}}{4 \cdot \sqrt{\pi} \cdot \left(100000000\right)}\simeq 1.4105\times 10^{-9}.
\end{aligned}
$$

相関は符号を反転して回復する。可視度だけでは負号は見えない。
:::

### 連続スペクトルと線スペクトルを同じ式で扱う

$\Gamma$ が連続で正値型なら、有限な非負測度 $\mu$ を用いて $\Gamma(\tau)=\int e^{-2\pi i\nu\tau}\dd\mu(\nu)$ と表せる。これは Bochner の定理である。非零の場合に $\Gamma(0)$ で割れば、確率測度を与える実直線上の Bochner の定理を適用できる。その証明は [実直線上の Bochner の定理](/math/probability/characteristic-functions/) を参照せよ。光の相関とスペクトルの関係は :cite[born-wolf] でも扱われる。全質量は $\Gamma(0)$、密度のある場合は $\dd\mu=S(\nu)\dd\nu$ となる。

この表示を有限時間の測定式に入れると、$\langle|\widehat u_T|^2\rangle/T$ は $\mu$ と $K_T(\nu)=T\operatorname{sinc}^2(\pi T\nu)$ の畳み込みである。$K_T$ は長さ $T$ の窓の Fourier 変換を二乗して $T$ で割ったものなので、Plancherel の等式から積分が 1 になる。さらに $|\nu|\ge\delta>0$ の質量は $2/(\pi^2T\delta)$ 以下である。従って連続な密度では窓を広げると連続なスペクトルを得る。点質量では一定面積の狭いピークとなって測度の意味で収束する。

## 空間的コヒーレンスと相互コヒーレンス関数

### 二点での相関
:term[相互コヒーレンス関数]{en="mutual coherence function" yomi="そうごこひーれんすかんすう"} を $\Gamma_{12}(\tau)=\langle u^*(\vb r_1,t)u(\vb r_2,t+\tau)\rangle$ とする。$\tau=0$ の値を相互強度と呼ぶ。$I_j=\Gamma_{jj}(0)>0$ なら $\gamma_{12}=\Gamma_{12}/\sqrt{I_1I_2}$ と規格化する。Cauchy–Schwarz より大きさは 1 以下で、二点の場が一定の複素係数で結ばれているときに等号を取る。同時刻の多点相関を行列に並べると半正定値になる。実際、任意の重み $c_j$ に対する二次形式は $\langle|\sum_j c_j u_j|^2\rangle$ である。

例えば単一のランダムな時間信号 $a(t)$ を各位置で $u_j=b_j a(t)$ と受ければ空間的な規格化相関の絶対値は 1 である。時間相関が短い光でも、幾何学的な遅延が揃えば空間相関が高くなり得る。

:::example[二スリットと相関]
等強度の二受光点で $\gamma_{12}(0)=i/2$ なら、位相を掃引したとき可視度は $1/2$、山を作る制御位相は $-\pi/2$ である。
:::
:::example[相関と伝播位相]
平面波を同時刻で二点測ると $\gamma_{12}=e^{i\vb k\cdot(\vb r_2-\vb r_1)}$。二点が離れても絶対値は 1 であり、位相が波面の向きを表す。
:::

### コヒーレンスの伝播
入力面と出力面の複素振幅を $u_o(\boldsymbol x)=\int h(\boldsymbol x,\boldsymbol s)u_i(\boldsymbol s)\dd^2s$ で結ぶ。二点の積を作って平均すると
$$
\Gamma_o(\boldsymbol x_1,\boldsymbol x_2)=
\iint h^*(\boldsymbol x_1,\boldsymbol s_1) h(\boldsymbol x_2,\boldsymbol s_2)
\Gamma_i(\boldsymbol s_1,\boldsymbol s_2)\dd^2s_1\dd^2s_2 .
$$
ここでは単色または共通遅延の近似で時刻差を省いた。有限エネルギーで期待値と積分を交換できると仮定した式である。入力が完全な空間コヒーレンスを持ち $\Gamma_i=f^*(\boldsymbol s_1)f(\boldsymbol s_2)$ なら、出力でも積は二つの独立した積分に分かれる。一方、非コヒーレントなデルタ相関では二変数が同じ位置に固定され、出力相関は一つの積分になる。この形から次節の van Cittert–Zernike の式を導く。

:::example[光量差から相関度を逆算する]
二腕の強度を $4,1$、可視度を 0.4 と測った場合、強度の不釣合いの係数は $2\sqrt4/5=0.8$。従って規格化相関の絶対値は $0.4/0.8=0.5$ となる。可視度が低い原因を空間相関だけに帰す前に、片腕ずつ遮って強度を測る。
:::

### 周波数ごとに相関を伝播させる

一般の定常場では周波数ごとに伝播核が異なる。時間差の Fourier 変換 $W_{12}(\nu)=\int\Gamma_{12}(\tau)e^{2\pi i\nu\tau}\dd\tau$ を :term[相互スペクトル密度]{en="cross-spectral density" yomi="そうごすぺくとるみつど"} と呼ぶ。二点に窓付き変換を取り、その積 $\widehat u_{1,T}(\nu)^*\widehat u_{2,T}(\nu)/T$ を平均すると、自己相関の場合と同じ三角の時間窓が掛かる。各 $\Gamma_{jk}$ が可積分なら、その窓を広げた極限が各 $W_{jk}$ である。有限窓の相関行列が半正定値なので、極限も半正定値であり、$|W_{12}|^2\le W_{11}W_{22}$ が従う。$W_{jj}$ が正の周波数では $W_{12}/\sqrt{W_{11}W_{22}}$ がその周波数でのコヒーレンス度になる。

時間に依存しない線形な光学系では異なる周波数は混ざらない。各 $\nu$ で先の二重積分の $\Gamma_i$ を $W_i(\nu)$ に、伝播核を $h_\nu$ に置き換える。周波数積分で $\Gamma_o(\tau)$ が戻るため、帯域で $h_\nu$ の位相が大きく変わる光学系でも遅延とコヒーレンスを同時に追える。

### 二つの独立モードの重なりとして相互相関を測る

二つの独立で平均ゼロのモード振幅 $a_j$ を分散一に規格化し、二点の場を $u_1=a_1\cos\theta+a_2\sin\theta$、$u_2=a_1\cos\theta+a_2e^{i\phi}\sin\theta$ とする。独立性により交差平均 $\langle a_1^*a_2\rangle=0$ なので、各点の照度は一、相互相関は $\cos^2\theta+e^{i\phi}\sin^2\theta$ になる。各モード内では二点が同一振幅を共有するのに、相関を足すと二モードの位相差で消え得る。一般形は、各モードの分散を重みにした二点の振幅ベクトルの内積である。絶対値が一になる条件は Cauchy–Schwarz の等号条件から、二点で見たモード振幅のベクトルが比例することである。

:::example[共通位相の二点]
二点で同じ線形結合を測る。

入力はモード分配の角 $\theta\simeq 0.7854\,\mathrm{rad}$、第二モードに付く位相 $\phi=0\,\mathrm{rad}$ である。

二点に共通に入る第一モードの重みを取る。

$$
\begin{aligned}
A&=\cos^{2}{\left(\theta \right)},\\
A&\simeq \cos^{2}{\left(\left(0.7853981634\right) \right)}=\frac{1}{2}.
\end{aligned}
$$

独立な第二モードの寄与を取る。

$$
\begin{aligned}
B&=\sin^{2}{\left(\theta \right)},\\
B&\simeq \sin^{2}{\left(\left(0.7853981634\right) \right)}=\frac{1}{2}.
\end{aligned}
$$

独立成分では電場を相関させた後で寄与を足す。

$$
\begin{aligned}
\gamma_{12}&=A + B e^{i \phi},\\
\gamma_{12}&= \left(\frac{1}{2}\right) + \left(\frac{1}{2}\right) \cdot e^{i \cdot \left(0\right)}=1.
\end{aligned}
$$

位相掃引で現れる可視度を求める。

$$
\begin{aligned}
V&=\left|{\gamma_{12}}\right|,\\
V&= \left|{\left(1\right)}\right|=1.
\end{aligned}
$$

各モードは独立でも二点での場は完全に対応する。
:::

:::example[相関の打ち消し]
等しい重みで第二モードだけに逆位相を入れる。

入力はモード分配の角 $\theta\simeq 0.7854\,\mathrm{rad}$、第二モードに付く位相 $\phi\simeq 3.1416\,\mathrm{rad}$ である。

二点に共通に入る第一モードの重みを取る。

$$
\begin{aligned}
A&=\cos^{2}{\left(\theta \right)},\\
A&\simeq \cos^{2}{\left(\left(0.7853981634\right) \right)}=\frac{1}{2}.
\end{aligned}
$$

独立な第二モードの寄与を取る。

$$
\begin{aligned}
B&=\sin^{2}{\left(\theta \right)},\\
B&\simeq \sin^{2}{\left(\left(0.7853981634\right) \right)}=\frac{1}{2}.
\end{aligned}
$$

独立成分では電場を相関させた後で寄与を足す。

$$
\begin{aligned}
\gamma_{12}&=A + B e^{i \phi},\\
\gamma_{12}&\simeq \left(\frac{1}{2}\right) + \left(\frac{1}{2}\right) \cdot e^{i \cdot \left(3.141592654\right)}=0.
\end{aligned}
$$

位相掃引で現れる可視度を求める。

$$
\begin{aligned}
V&=\left|{\gamma_{12}}\right|,\\
V&= \left|{\left(0\right)}\right|=0.
\end{aligned}
$$

各点の照度を変えずに複素相関をゼロにできる。
:::

:::example[片方のモードを強める]
位相は前例のまま第一モードの寄与を多くする。

入力はモード分配の角 $\theta\simeq 0.5236\,\mathrm{rad}$、第二モードに付く位相 $\phi\simeq 3.1416\,\mathrm{rad}$ である。

二点に共通に入る第一モードの重みを取る。

$$
\begin{aligned}
A&=\cos^{2}{\left(\theta \right)},\\
A&\simeq \cos^{2}{\left(\left(0.5235987756\right) \right)}=\frac{3}{4}.
\end{aligned}
$$

独立な第二モードの寄与を取る。

$$
\begin{aligned}
B&=\sin^{2}{\left(\theta \right)},\\
B&\simeq \sin^{2}{\left(\left(0.5235987756\right) \right)}=\frac{1}{4}.
\end{aligned}
$$

独立成分では電場を相関させた後で寄与を足す。

$$
\begin{aligned}
\gamma_{12}&=A + B e^{i \phi},\\
\gamma_{12}&\simeq \left(\frac{3}{4}\right) + \left(\frac{1}{4}\right) \cdot e^{i \cdot \left(3.141592654\right)}=\frac{1}{2}.
\end{aligned}
$$

位相掃引で現れる可視度を求める。

$$
\begin{aligned}
V&=\left|{\gamma_{12}}\right|,\\
V&= \left|{\left(\frac{1}{2}\right)}\right|=\frac{1}{2}.
\end{aligned}
$$

打ち消しが不完全で、可視度が戻る。
:::

## van Cittert–Zernike の定理

空間的に独立な微小光源を平面に並べ、面から距離 $L$ の受光面を考える。角度広がりを小さくし、準単色で各光源から受光点への伝播振幅の変化を無視する。光源の横座標を $\boldsymbol s$、二受光点を $\boldsymbol\rho_j$ とすると、光路長 $R_j\simeq L+|\boldsymbol\rho_j-\boldsymbol s|^2/(2L)$ である。異なる光源間で平均交差項が消えるので、同じ光源からの積を強度密度 $B(\boldsymbol s)$ で重み付けして積分する。$\boldsymbol B=\boldsymbol\rho_2-\boldsymbol\rho_1$ と置いて
$$
\gamma_{12}(0)=e^{ik(\rho_2^2-\rho_1^2)/(2L)}
\frac{\int B(\boldsymbol s)e^{-ik\boldsymbol B\cdot\boldsymbol s/L}\dd^2s}
{\int B(\boldsymbol s)\dd^2s} .
\label{eq-optics-vcz}
$$
これが :term[van Cittert–Zernike の定理]{en="van Cittert–Zernike theorem" yomi="ばんしったーとぜるにけのていり"} の近軸形である。曲率位相を除くと、光源の角座標 $\boldsymbol s/L$ での輝度を空間周波数 $\boldsymbol B/\lambda$ で変換した値になる。相関を落とさないため、帯域による光路差の平均化が無視できるかも確認する。

:::example[二点源]
同強度の二点が角 $\pm\alpha$ にある一次元模型で、曲率補正後の $\gamma$ は $\cos(2\pi B\alpha/\lambda)$。$B=\lambda/(4\alpha)$ で初めてゼロになる。
:::
:::example[角度の中心をずらす]
輝度を角 $\alpha_0$ だけ平行移動すれば相関に $e^{-2\pi i B\alpha_0/\lambda}$ が掛かる。絶対値だけでは移動を測れない。
:::

:::example[矩形の非コヒーレント光源]
一次元に一様な幅 $a$ の光源が距離 $L$ にある。式 :eqref[eq-optics-vcz] で曲率位相を補正し $s\in[-a/2,a/2]$ を積分すると $\gamma(B)=\operatorname{sinc}(\pi Ba/(\lambda L))$。最初の零点は $\lambda L/a$ なので、光源の角幅を二倍にすると受光面の相関範囲は半分になる。
:::

### 帯域と幾何学遅延を同時に残す
各光源要素のスペクトル密度を $B_\nu(\boldsymbol s)$ とする。二受光点への光路差を $\Delta R(\boldsymbol s)=R_2-R_1$ とし、伝播振幅を共通とみなす。スペクトルごとの相関を足せば、未規格化の時間相関は
$$
\Gamma_{12}(\tau)\propto\int\!\!\int B_\nu(\boldsymbol s)
e^{-2\pi i\nu[\tau-\Delta R(\boldsymbol s)/c]}\dd\nu\,\dd^2s .
$$
$\Gamma$ のピークを読むなら、中央方向の幾何学遅延に $\tau$ を合わせる。残る遅延は角方向ごとに異なり、狭い帯域でその差がコヒーレンス時間より小さいとき、式 :eqref[eq-optics-vcz] の単色近似を使える。広帯域のまま単一の有効波長で相関を予測すると、離れた方向からの縞の低下を説明できない。

### 一様な線光源の角幅を基線から測る

遠方の小さい角幅を持つ一次元光源を角 $\vartheta$ で表す。長さ $B$ の基線に沿う射影の位相因子は $e^{-2\pi iB\vartheta/\lambda}$ である。光源が中心角 $\vartheta_0$、全角幅 $\Theta$ の一様区間なら、区間の積分から相関は $e^{-2\pi iB\vartheta_0/\lambda}\operatorname{sinc}(\pi B\Theta/\lambda)$。基線の向きと角の正方向は式 :eqref[eq-optics-vcz] と同じである。基線を反転すると共役になる。相関の絶対値を測ると中心の平行移動は分からず、位相も校正すると中心位置を求められる。零点基線 $\lambda/\Theta$ は円盤光源の零点の係数とは異なる。

:::example[中心を測定軸に置く]
基線が最初の零点より短い一様光源を観測する。

入力は基線 $B=5\,\mathrm m$、波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、光源の全角幅 $\Theta\simeq 5\times 10^{-8}\,\mathrm{rad}$、中心角 $\vartheta_0=0\,\mathrm{rad}$ である。

光源の端の間の相対位相の半分を計算する。

$$
\begin{aligned}
v&=\frac{\pi B \Theta}{\lambda},\\
v&\simeq \frac{\pi \cdot \left(5\right) \cdot \left(5\times 10^{-8}\right)}{\left(5\times 10^{-7}\right)}\simeq 1.5708.
\end{aligned}
$$

角幅と独立に中心の位相を求める。

$$
\begin{aligned}
\phi_c&=-\frac{2 \pi B \vartheta_0}{\lambda},\\
\phi_c/(\mathrm{rad})&\simeq \frac{2 \cdot \pi \cdot \left(5\right) \cdot \left(0\right)}{\left(5\times 10^{-7}\right)}=0.
\end{aligned}
$$

中心位相と一様光源の Fourier 因子を掛ける。

$$
\begin{aligned}
\gamma&=e^{i \phi_c} \operatorname{sinc}{\left(v \right)},\\
\gamma&\simeq e^{i \cdot \left(0\right)} \cdot \operatorname{sinc}{\left(\left(1.570796327\right) \right)}\simeq 0.63662.
\end{aligned}
$$

中心位相は可視度に寄与しない。

$$
\begin{aligned}
\mathcal V&=\left|{\gamma}\right|,\\
\mathcal V&\simeq \left|{\left(0.6366197724\right)}\right|\simeq 0.63662.
\end{aligned}
$$

最初の sinc の零点を基線へ換算する。

$$
\begin{aligned}
B_1&=\frac{\lambda}{\Theta},\\
B_1/(\mathrm m)&\simeq \frac{\left(5\times 10^{-7}\right)}{\left(5\times 10^{-8}\right)}=10.
\end{aligned}
$$

実数の正の相関が得られる。
:::

:::example[中心を平行移動する]
角幅と波長を保ったまま見かけの中心をずらす。

入力は基線 $B=5\,\mathrm m$、波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、光源の全角幅 $\Theta\simeq 5\times 10^{-8}\,\mathrm{rad}$、中心角 $\vartheta_0\simeq 2.5\times 10^{-8}\,\mathrm{rad}$ である。

光源の端の間の相対位相の半分を計算する。

$$
\begin{aligned}
v&=\frac{\pi B \Theta}{\lambda},\\
v&\simeq \frac{\pi \cdot \left(5\right) \cdot \left(5\times 10^{-8}\right)}{\left(5\times 10^{-7}\right)}\simeq 1.5708.
\end{aligned}
$$

角幅と独立に中心の位相を求める。

$$
\begin{aligned}
\phi_c&=-\frac{2 \pi B \vartheta_0}{\lambda},\\
\phi_c/(\mathrm{rad})&\simeq -\frac{2 \cdot \pi \cdot \left(5\right) \cdot \left(2.5\times 10^{-8}\right)}{\left(5\times 10^{-7}\right)}\simeq -1.5708.
\end{aligned}
$$

中心位相と一様光源の Fourier 因子を掛ける。

$$
\begin{aligned}
\gamma&=e^{i \phi_c} \operatorname{sinc}{\left(v \right)},\\
\gamma&\simeq e^{i \cdot \left(-1.570796327\right)} \cdot \operatorname{sinc}{\left(\left(1.570796327\right) \right)}\simeq -0.63662i.
\end{aligned}
$$

中心位相は可視度に寄与しない。

$$
\begin{aligned}
\mathcal V&=\left|{\gamma}\right|,\\
\mathcal V&\simeq \left|{\left(-0.6366197724i\right)}\right|\simeq 0.63662.
\end{aligned}
$$

最初の sinc の零点を基線へ換算する。

$$
\begin{aligned}
B_1&=\frac{\lambda}{\Theta},\\
B_1/(\mathrm m)&\simeq \frac{\left(5\times 10^{-7}\right)}{\left(5\times 10^{-8}\right)}=10.
\end{aligned}
$$

可視度は変わらず、複素位相だけが変わる。
:::

:::example[基線を伸ばして零点を測る]
基線を二倍にして光源幅の条件を読む。

入力は基線 $B=10\,\mathrm m$、波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、光源の全角幅 $\Theta\simeq 5\times 10^{-8}\,\mathrm{rad}$、中心角 $\vartheta_0\simeq 2.5\times 10^{-8}\,\mathrm{rad}$ である。

光源の端の間の相対位相の半分を計算する。

$$
\begin{aligned}
v&=\frac{\pi B \Theta}{\lambda},\\
v&\simeq \frac{\pi \cdot \left(10\right) \cdot \left(5\times 10^{-8}\right)}{\left(5\times 10^{-7}\right)}\simeq 3.1416.
\end{aligned}
$$

角幅と独立に中心の位相を求める。

$$
\begin{aligned}
\phi_c&=-\frac{2 \pi B \vartheta_0}{\lambda},\\
\phi_c/(\mathrm{rad})&\simeq -\frac{2 \cdot \pi \cdot \left(10\right) \cdot \left(2.5\times 10^{-8}\right)}{\left(5\times 10^{-7}\right)}\simeq -3.1416.
\end{aligned}
$$

中心位相と一様光源の Fourier 因子を掛ける。

$$
\begin{aligned}
\gamma&=e^{i \phi_c} \operatorname{sinc}{\left(v \right)},\\
\gamma&\simeq e^{i \cdot \left(-3.141592654\right)} \cdot \operatorname{sinc}{\left(\left(3.141592654\right) \right)}=0.
\end{aligned}
$$

中心位相は可視度に寄与しない。

$$
\begin{aligned}
\mathcal V&=\left|{\gamma}\right|,\\
\mathcal V&= \left|{\left(0\right)}\right|=0.
\end{aligned}
$$

最初の sinc の零点を基線へ換算する。

$$
\begin{aligned}
B_1&=\frac{\lambda}{\Theta},\\
B_1/(\mathrm m)&\simeq \frac{\left(5\times 10^{-7}\right)}{\left(5\times 10^{-8}\right)}=10.
\end{aligned}
$$

絶対値がゼロなので、この基線では中心位置を位相から読めない。
:::

## 恒星干渉計

### 円盤の角直径を測る
離れて置いた二望遠鏡で同じ天体を集め、腕の遅延を補償して干渉させる装置を :term[恒星干渉計]{en="stellar interferometer" yomi="こうせいかんしょうけい"} と呼ぶ。基線は天体方向に垂直な面への射影で測る。一様円盤の角直径を $\Theta$ とすると、式 :eqref[eq-optics-vcz] の規格化変換は回折の円孔と同じ積分なので
$$
\gamma(B)=\frac{2J_1(\pi B\Theta/\lambda)}{\pi B\Theta/\lambda}.
$$
第一零点は $B\simeq1.22\lambda/\Theta$ である。縁の減光がある天体ではこの円盤モデルがずれるので、複数基線の測定を輝度モデルに当てはめる。同様に連星では基線の向きを変えると、相関の周期から射影分離を測れる。大気で二腕の位相が変動すると、符号や複素位相は強度の可視度より測りにくい。

:::example[角直径の換算]
$\lambda=500\,\mathrm{nm}$、第一零点が $B=100\,\mathrm m$ なら $\Theta\simeq6.1\times10^{-9}\,\mathrm{rad}$、約 $1.26$ ミリ秒角である。
:::
:::example[同じ物体で波長を変える]
角直径を固定して波長を倍にすると、第一零点に達する基線長も倍になる。物体の角直径が波長で変わる場合は、この単純な換算は使えない。
:::

### 周縁減光がある円盤の可視度

一様円盤の代わりに、角半径を 1 に規格化した半径 $r$ で、輝度を $I(r)=I_0[(1-u)+u\sqrt{1-r^2}]$ とする。$0\le u\le1$ は :term[周縁減光]{en="limb darkening" yomi="しゅうえんげんこう"} の強さを表す。$x=\pi B\Theta/\lambda$ とすれば、角積分により規格化相関は
$$
\gamma(x)=\frac{(1-u)J_1(x)/x+u(\sin x-x\cos x)/x^3}{(1-u)/2+u/3}.
$$
第二項の積分 $\int_0^1 r\sqrt{1-r^2}J_0(xr)\dd r$ は、単位球に一様な密度を置き、その視線方向の積分を Fourier 変換したものに比例する。球の断面積 $\pi(1-t^2)$ を別の軸 $t$ に沿って積分すれば、この項は $\tfrac14\int_{-1}^1(1-t^2)e^{-ixt}\dd t=(\sin x-x\cos x)/x^3$ となる。両項の $x\to0$ の極限は $1/2,1/3$ である。$u=0$ では一様円盤に戻り、$u=1$ では第一正根が $\tan x=x$ の約 $4.4934$ に移る。同じ角直径でも、この二つの模型では第一零点を観測する基線が異なる。

### 三つの基線で位相を組み合わせる
受光点 $j$ ごとに大気が未知の位相 $\phi_j$ を足すと、定義した順序では観測される $\Gamma_{jk}$ に $e^{i(\phi_k-\phi_j)}$ が掛かる。閉じた三角形の積 $\Gamma_{12}\Gamma_{23}\Gamma_{31}$ では位相の差を足してゼロになる。その積の偏角は各望遠鏡に固有な位相の影響を受けない。この閉包位相は三つの基線の位相和を測る量であり、個々の基線の位相を一つずつ取り出すものではない。

:::example[基線の射影を換算する]
物理的な基線が $100\,\mathrm m$ で天体の方向となす角が $30^\circ$ なら、波面に平行な有効成分は $50\,\mathrm m$、幾何学的遅延を作る成分は $100\cos30^\circ\,\mathrm m$ である。円盤径の式には前者を、遅延線には後者を使う。
:::

### 等光度の連星から周期と射影分離を求める

二つの点光源の中心を角座標の原点とし、基線に射影した角距離を $\Theta$ とする。等光度なら規格化された二つの位相因子の和が $\cos(\pi B\Theta/\lambda)$。最小の正の零点は $B=\lambda/(2\Theta)$ である。光度が違う場合、第一光源の割合を $f$ とすれば相関の二乗は $1-4f(1-f)\sin^2(\pi B\Theta/\lambda)$ になる。これは相関を逆方向の位相因子の加重和とし、絶対値の二乗を取ったもの。最小値は $(2f-1)^2$ なので、不等光度では零点がなくても極小の基線は分離を与える。基線に直交する距離成分は、向きの異なる基線を使わないと決まらない。

:::example[等光度の連星]
基線を最初の極小に選ぶ。

入力は基線長 $B=25\,\mathrm m$、波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、基線方向の角距離 $\Theta\simeq 1\times 10^{-8}\,\mathrm{rad}$、第一の星の光度割合 $f=\frac{1}{2}$ である。

二つの位置の基線への射影を位相に換算する。

$$
\begin{aligned}
v&=\frac{\pi B \Theta}{\lambda},\\
v&\simeq \frac{\pi \cdot \left(25\right) \cdot \left(1\times 10^{-8}\right)}{\left(5\times 10^{-7}\right)}\simeq 1.5708.
\end{aligned}
$$

相関の二乗の振れは光度の積に比例する。

$$
\begin{aligned}
g&=4 f \left(1 - f\right),\\
g&= 4 \cdot \left(\frac{1}{2}\right) \cdot \left(1 - \left(\frac{1}{2}\right)\right)=1.
\end{aligned}
$$

二つの点からの複素相関を加えて絶対値を取る。

$$
\begin{aligned}
V^2&=- g \sin^{2}{\left(v \right)} + 1,\\
V^2&\simeq - \left(1\right) \cdot \sin^{2}{\left(\left(1.570796327\right) \right)} + 1=0.
\end{aligned}
$$

振幅干渉計の可視度を求める。

$$
\begin{aligned}
V&=\sqrt{V^2},\\
V&= \sqrt{\left(0\right)}=0.
\end{aligned}
$$

同じ射影角距離で逆位相になる基線を求める。

$$
\begin{aligned}
B_{\min}&=\frac{\lambda}{2 \Theta},\\
B_{\min}/(\mathrm m)&\simeq \frac{\left(5\times 10^{-7}\right)}{2 \cdot \left(1\times 10^{-8}\right)}=25.
\end{aligned}
$$

二つの寄与は同じ大きさで完全に相殺する。
:::

:::example[不等光度の連星]
同じ分離で一方の光度を低くする。

入力は基線長 $B=25\,\mathrm m$、波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、基線方向の角距離 $\Theta\simeq 1\times 10^{-8}\,\mathrm{rad}$、第一の星の光度割合 $f=\frac{4}{5}$ である。

二つの位置の基線への射影を位相に換算する。

$$
\begin{aligned}
v&=\frac{\pi B \Theta}{\lambda},\\
v&\simeq \frac{\pi \cdot \left(25\right) \cdot \left(1\times 10^{-8}\right)}{\left(5\times 10^{-7}\right)}\simeq 1.5708.
\end{aligned}
$$

相関の二乗の振れは光度の積に比例する。

$$
\begin{aligned}
g&=4 f \left(1 - f\right),\\
g&= 4 \cdot \left(\frac{4}{5}\right) \cdot \left(1 - \left(\frac{4}{5}\right)\right)=\frac{16}{25}.
\end{aligned}
$$

二つの点からの複素相関を加えて絶対値を取る。

$$
\begin{aligned}
V^2&=- g \sin^{2}{\left(v \right)} + 1,\\
V^2&\simeq - \left(\frac{16}{25}\right) \cdot \sin^{2}{\left(\left(1.570796327\right) \right)} + 1=\frac{9}{25}.
\end{aligned}
$$

振幅干渉計の可視度を求める。

$$
\begin{aligned}
V&=\sqrt{V^2},\\
V&= \sqrt{\left(\frac{9}{25}\right)}=\frac{3}{5}.
\end{aligned}
$$

同じ射影角距離で逆位相になる基線を求める。

$$
\begin{aligned}
B_{\min}&=\frac{\lambda}{2 \Theta},\\
B_{\min}/(\mathrm m)&\simeq \frac{\left(5\times 10^{-7}\right)}{2 \cdot \left(1\times 10^{-8}\right)}=25.
\end{aligned}
$$

極小の位置は同じだが、相関が消え残る。
:::

:::example[別の角度から連星を見る]
基線を傾けて射影された分離を小さくする。

入力は基線長 $B=25\,\mathrm m$、波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、基線方向の角距離 $\Theta\simeq 5\times 10^{-9}\,\mathrm{rad}$、第一の星の光度割合 $f=\frac{4}{5}$ である。

二つの位置の基線への射影を位相に換算する。

$$
\begin{aligned}
v&=\frac{\pi B \Theta}{\lambda},\\
v&\simeq \frac{\pi \cdot \left(25\right) \cdot \left(5\times 10^{-9}\right)}{\left(5\times 10^{-7}\right)}\simeq 0.7854.
\end{aligned}
$$

相関の二乗の振れは光度の積に比例する。

$$
\begin{aligned}
g&=4 f \left(1 - f\right),\\
g&= 4 \cdot \left(\frac{4}{5}\right) \cdot \left(1 - \left(\frac{4}{5}\right)\right)=\frac{16}{25}.
\end{aligned}
$$

二つの点からの複素相関を加えて絶対値を取る。

$$
\begin{aligned}
V^2&=- g \sin^{2}{\left(v \right)} + 1,\\
V^2&\simeq - \left(\frac{16}{25}\right) \cdot \sin^{2}{\left(\left(0.7853981634\right) \right)} + 1=\frac{17}{25}.
\end{aligned}
$$

振幅干渉計の可視度を求める。

$$
\begin{aligned}
V&=\sqrt{V^2},\\
V&= \sqrt{\left(\frac{17}{25}\right)}\simeq 0.82462.
\end{aligned}
$$

同じ射影角距離で逆位相になる基線を求める。

$$
\begin{aligned}
B_{\min}&=\frac{\lambda}{2 \Theta},\\
B_{\min}/(\mathrm m)&\simeq \frac{\left(5\times 10^{-7}\right)}{2 \cdot \left(5\times 10^{-9}\right)}=50.
\end{aligned}
$$

同じ基線長では極小に達しない。射影角を変えると別の曲線になる。
:::

## 強度干渉と Hanbury Brown–Twiss 効果

### 強度相関と Gaussian 場
二検出器で $I_1(t),I_2(t+\tau)$ を測り、:term[強度相関]{en="intensity correlation" yomi="きょうどそうかん"} を
$$
g^{(2)}_{12}(\tau)=\frac{\langle I_1(t)I_2(t+\tau)\rangle}{\langle I_1\rangle\langle I_2\rangle}
$$
と定義する。二時刻の複素振幅 $u_1,u_2$ が中心化された円対称複素 Gaussian とする。円対称性から $\langle u_i u_j\rangle=0$。四次モーメントは二次モーメントの対の積の和である。実 Gaussian の母関数 $\exp(\vb t^\mathsf T C\vb t/2)$ を四回微分すると三通りの対分けが現れ、それを実部と虚部の線形結合に適用すればよい。従って $\langle u_1^*u_1 u_2^*u_2\rangle=\langle I_1\rangle\langle I_2\rangle+|\Gamma_{12}|^2$ となり、
$$
g^{(2)}_{12}(\tau)=1+|\gamma_{12}(\tau)|^2
\label{eq-optics-siegert}
$$
を得る。これが Siegert の関係である。熱光の :term[Hanbury Brown–Twiss 効果]{en="Hanbury Brown–Twiss effect" yomi="はんべりーぶらうんついすこうか"} では、この余分な同時計数から空間コヒーレンスの大きさを推定する。相互相関の複素位相はこの二点の測定からは復元できない。また検出時間分解能が相関幅より悪いと、観測される相関の山は時間窓で平均され、高さが低下する。二検出器で光子を数える場合、各検出器固有のショット雑音と、二検出器間の相関を区別する。

:::example[同一点の熱光と一定強度光]
単一偏光単一モード熱光の同時刻では $\gamma_{11}(0)=1$、従って $g^{(2)}(0)=2$。一定強度の場では直接定義から $g^{(2)}=1$ であり、式 :eqref[eq-optics-siegert] の Gaussian 仮定を満たさない。
:::
:::example[独立な二偏光]
同強度で独立な熱光二偏光を識別せずに測る。各強度の分散はその平均の二乗であり、和の分散は二倍だが、平均の二乗は四倍になる。従って $g^{(2)}(0)=1+2/4=3/2$ である。
:::

### 有限検出時間の計算
二検出器が同じ幅 $\Delta$ の時間窓で強度を平均する。測定の零遅延の余分な相関は、式 :eqref[eq-optics-siegert] の余分な部分を $(\Delta-|\tau|)/\Delta^2$（$|\tau|\le\Delta$）で重み付けした積分になる。二つの積分変数の差を固定すると残りの区間長が $\Delta-|\tau|$ だからである。相関幅が $\Delta$ より十分短いなら、重みは相関が大きい部分でほぼ $1/\Delta$。従って観測される高さの余分は $\Delta^{-1}\int|\gamma_{12}(\tau)|^2\dd\tau$ に近づく。同一点の測定なら、この積分は本章で定義した $\tau_c$ である。瞬時の規格化相関が 2 の熱光でも、長い露光で得られる値は 1 に近い。

:::example[積分時間の尺度]
例えば単一熱モードの $\tau_c=1\,\mathrm{ps}$ に対して $\Delta=1\,\mathrm{ns}$ なら、余分な相関はおよそ $10^{-3}$ に薄まる。空間相関の絶対値二乗がさらに $1/4$ であれば、余分は約 $2.5\times10^{-4}$ となる。
:::

### 光子計数での規格化の違い
同じ検出器の短い窓の計数を $N$ とすると、異なる二事象の組の数は $N(N-1)$ である。$N^2$ では同じ事象を二回数える $N$ が加わる。従って同一モードの零遅延の量子光の規格化相関は $\langle N(N-1)\rangle/\langle N\rangle^2$ を用いる。Poisson 分布で平均 $\mu$ なら階乗モーメントは $\mu^2$、値は 1。各試行で必ず一光子なら $N(N-1)=0$、値は 0 である。さらに、非負の古典強度で瞬時の同一点の比を定義すれば、分散の非負性により $\langle I^2\rangle/\langle I\rangle^2\ge1$。従って 1 を下回る量子相関は、古典のランダム強度だけの混合では表せない。

### 有限露光による強度相関の平均

熱光の正規化強度相関が $g^{(2)}(t)-1=e^{-2|t|/\tau}$ とする。同じ長さ $T$ の窓で二検出器の強度を平均すると、ゼロ遅延の規格化共分散は $T^{-2}\int_0^T\int_0^T e^{-2|t-s|/\tau}\dd s\dd t$。差時刻 $u=|t-s|$ を固定した測度は $2(T-u)\dd u$ なので、結果は $2(x-1+e^{-x})/x^2$、$x=2T/\tau$ である。$T\ll\tau$ なら一に近く、$T\gg\tau$ なら $\tau/T$ 程度まで薄まる。独立な検出器の共分散を用いれば同じ検出イベントの自己相関を避けられるが、時間平均による低下は残る。

:::example[十分短い計数窓]
熱光の相関の時間に比べて短い窓を選ぶ。

入力は窓の長さ $T\simeq 1\times 10^{-11}\,\mathrm s$、場の指数相関の時間 $\tau\simeq 1\times 10^{-9}\,\mathrm s$ である。

強度相関の指数率で窓の長さを換算する。

$$
\begin{aligned}
x&=\frac{2 T}{\tau},\\
x&\simeq \frac{2 \cdot \left(1\times 10^{-11}\right)}{\left(1\times 10^{-9}\right)}=\frac{1}{50}.
\end{aligned}
$$

最大の時刻差での相関を計算する。

$$
\begin{aligned}
q&=e^{- x},\\
q&= e^{- \left(\frac{1}{50}\right)}\simeq 0.9802.
\end{aligned}
$$

窓内での差時刻の三角形の測度を入れる。

$$
\begin{aligned}
g_T^{(2)}(0)-1&=\frac{2 q + 2 x - 2}{x^{2}},\\
g_T^{(2)}(0)-1&\simeq \frac{2 \cdot \left(0.9801986733\right) + 2 \cdot \left(\frac{1}{50}\right) - 2}{\left(\frac{1}{50}\right)^{2}}\simeq 0.99337.
\end{aligned}
$$

独立な背景の積を表す一を足す。

$$
\begin{aligned}
g_T^{(2)}(0)&=g_T^{(2)}(0)-1 + 1,\\
g_T^{(2)}(0)&\simeq \left(0.9933665338\right) + 1\simeq 1.9934.
\end{aligned}
$$

ゼロ遅延のバンチング値二に近い。
:::

:::example[相関時間と同じ窓]
相関の変化を一窓の中で平均する。

入力は窓の長さ $T\simeq 1\times 10^{-9}\,\mathrm s$、場の指数相関の時間 $\tau\simeq 1\times 10^{-9}\,\mathrm s$ である。

強度相関の指数率で窓の長さを換算する。

$$
\begin{aligned}
x&=\frac{2 T}{\tau},\\
x&\simeq \frac{2 \cdot \left(1\times 10^{-9}\right)}{\left(1\times 10^{-9}\right)}=2.
\end{aligned}
$$

最大の時刻差での相関を計算する。

$$
\begin{aligned}
q&=e^{- x},\\
q&= e^{- \left(2\right)}\simeq 0.13534.
\end{aligned}
$$

窓内での差時刻の三角形の測度を入れる。

$$
\begin{aligned}
g_T^{(2)}(0)-1&=\frac{2 q + 2 x - 2}{x^{2}},\\
g_T^{(2)}(0)-1&\simeq \frac{2 \cdot \left(0.1353352832\right) + 2 \cdot \left(2\right) - 2}{\left(2\right)^{2}}\simeq 0.56767.
\end{aligned}
$$

独立な背景の積を表す一を足す。

$$
\begin{aligned}
g_T^{(2)}(0)&=g_T^{(2)}(0)-1 + 1,\\
g_T^{(2)}(0)&\simeq \left(0.5676676416\right) + 1\simeq 1.5677.
\end{aligned}
$$

測定値は二より小さくなる。熱光であることと測定値が二になることは同じ条件ではない。
:::

:::example[多数の独立な時間成分を含む窓]
露光をさらに長くする。

入力は窓の長さ $T\simeq 1\times 10^{-8}\,\mathrm s$、場の指数相関の時間 $\tau\simeq 1\times 10^{-9}\,\mathrm s$ である。

強度相関の指数率で窓の長さを換算する。

$$
\begin{aligned}
x&=\frac{2 T}{\tau},\\
x&\simeq \frac{2 \cdot \left(1\times 10^{-8}\right)}{\left(1\times 10^{-9}\right)}=20.
\end{aligned}
$$

最大の時刻差での相関を計算する。

$$
\begin{aligned}
q&=e^{- x},\\
q&= e^{- \left(20\right)}\simeq 2.0612\times 10^{-9}.
\end{aligned}
$$

窓内での差時刻の三角形の測度を入れる。

$$
\begin{aligned}
g_T^{(2)}(0)-1&=\frac{2 q + 2 x - 2}{x^{2}},\\
g_T^{(2)}(0)-1&\simeq \frac{2 \cdot \left(2.061153622\times 10^{-9}\right) + 2 \cdot \left(20\right) - 2}{\left(20\right)^{2}}\simeq 0.095.
\end{aligned}
$$

独立な背景の積を表す一を足す。

$$
\begin{aligned}
g_T^{(2)}(0)&=g_T^{(2)}(0)-1 + 1,\\
g_T^{(2)}(0)&\simeq \left(0.09500000001\right) + 1\simeq 1.095.
\end{aligned}
$$

規格化した余剰の大きさは相関時間と窓の長さの比で見積もれる。
:::

### 不等強度の独立モードと有効なモード数

独立な熱モードの瞬時強度を $I_j$ とし、その平均を $\mu_j>0$ とする。複素円対称 Gaussian 場の実部と虚部が独立で同分散なら、各強度は二個の標準正規変数の二乗和に比例する。極座標で二次元の正規密度を積分すると、$I_j$ の密度は $\mu_j^{-1}e^{-I_j/\mu_j}$ である。この密度から $\langle I_j^2\rangle=2\mu_j^2$、従って分散は $\mu_j^2$ となる。異なる偏光や識別しない直交モードを検出器が強度として加える場合は、独立性により分散を足せる。この条件で、同一点・零遅延の値は
$$
g^{(2)}(0)=1+\frac{\sum_j\mu_j^2}{(\sum_j\mu_j)^2}
=1+\frac{1}{M_{\rm eff}},\qquad
M_{\rm eff}=\frac{(\sum_j\mu_j)^2}{\sum_j\mu_j^2}.
\label{eq-optics-effective-thermal-modes}
$$
Cauchy–Schwarz で $(\sum_j\mu_j)^2\le M\sum_j\mu_j^2$ なので $M_{\rm eff}\le M$。強度が均一なときに限り等号となる。強度が一モードに集中する極限では $M_{\rm eff}\to1$。有限露光でも余剰相関が減るため、時間平均を補正せず観測値だけからモード数を逆算すると混同が起こる。

:::example[偏光を区別しない二モード]
直交する二偏光の熱光が独立で、平均強度の比が $9:1$ とする。式 :eqref[eq-optics-effective-thermal-modes] の余剰は $(81+1)/100=0.82$、従って $g^{(2)}(0)=1.82$。同強度の二偏光なら値は $1.5$ である。同じ二モードでも強度比が大きいと、相関の読み出しは一モードの値に近い。
:::

## 演習問題

::::exercise{level=basic}
同じ強度の二光束で規格化相関の絶対値が $0.4$ のとき可視度を求めよ。

:::solution
可視度は相関の絶対値に $2\sqrt{I_1I_2}/(I_1+I_2)$ を掛ける。この係数が 1 なので $0.4$ である。
:::
::::

::::exercise{level=basic}
矩形のスペクトル幅が $10^9\,\mathrm{Hz}$ のとき可視度の第一零点の遅延を求めよ。

:::solution
$\operatorname{sinc}(\pi\Delta\nu\tau)$ の第一零点で $\Delta\nu\tau=1$。よって $1\,\mathrm{ns}$ である。
:::
::::

::::exercise{level=basic}
定常場の自己相関について $\Gamma(-\tau)=\Gamma(\tau)^*$ を示せ。

:::solution
$\Gamma(-\tau)=\langle u^*(t)u(t-\tau)\rangle=\langle u^*(t+\tau)u(t)\rangle$。最後は時刻の原点を移した式で、これは $\Gamma(\tau)^*$ である。
:::
::::

::::exercise{level=basic}
熱光で基線上の $\abs{\gamma}=1/2$ のとき $g^{(2)}$ を求めよ。

:::solution
式 :eqref[eq-optics-siegert] の交差項は $(1/2)^2$ なので $g^{(2)}=5/4$ である。
:::
::::

::::exercise{level=advanced}
指数相関の係数が $a=2\times10^6\,\mathrm{s^{-1}}$ のとき、定義したコヒーレンス時間を求めよ。

:::solution
絶対値の二乗は $e^{-2a\abs{\tau}}$。両側を積分すると $1/a=0.5\,\mathrm{\mu s}$ となる。
:::
::::

::::exercise{level=advanced}
$\gamma_{12}$ の絶対値が 1 のとき場の関係を証明せよ。

:::solution
固定した遅延 $\tau$ において $u_1=u(\vb r_1,t)$、$u_2=u(\vb r_2,t+\tau)$ とする。$a=\Gamma_{12}(\tau)/I_1$ と置く。二乗誤差 $\langle\abs{u_2-au_1}^2\rangle=I_2-\abs{\Gamma_{12}}^2/I_1$ が仮定でゼロなので、二次モーメントの意味で $u_2=au_1$ となる。
:::
::::

::::exercise{level=advanced}
相互に独立な等強度の熱光を $M$ モード識別せずに足す。正規化強度相関の零遅延値を求めよ。

:::solution
各モードの平均を $\mu$ とすれば分散は $\mu^2$。和の平均は $M\mu$、分散は $M\mu^2$。二次モーメントは平均の二乗と分散の和なので、割ると $1+1/M$ となる。
:::
::::

::::exercise{level=advanced}
離れて配置した同強度二点源の角の差が $\theta$ で、同方向の基線を走査する。第一零点の基線を求めよ。

:::solution
角 $\pm\theta/2$ と置ける。式 :eqref[eq-optics-vcz] で曲率位相を除くと $\cos(\pi B\theta/\lambda)$。最小の正根から $B=\lambda/(2\theta)$ である。
:::
::::

::::exercise{level=advanced}
モードを分離しない検出器に、平均強度の比が $2:1:1$ の独立な熱モードを入れる。時間平均による低下がないとして $M_{\rm eff}$ と零遅延の $g^{(2)}$ を求めよ。

:::solution
式 :eqref[eq-optics-effective-thermal-modes] において平均強度の和に比例する値は $2+1+1=4$、二乗和に比例する値は $4+1+1=6$。従って $M_{\rm eff}=16/6=8/3$ であり、$g^{(2)}(0)=1+3/8=11/8$ となる。三つの強度が等しいときの $M_{\rm eff}=3$ より小さい理由は、強いモードの分散が平均の二乗で効くためである。
:::
::::

## 確認テスト

::::quiz
時間相関とスペクトルの関係を選べ。

- [x] Fourier 変換の対である。
- [ ] 常に同じ数値になる。
- [ ] 線幅が広いほど相関時間が長い。

:::explanation
式 :eqref[eq-optics-wiener-khinchin] の積分でスペクトルの各成分は異なる位相で足される。帯域を広げると、同じ形を相似に拡大した場合は時間相関が短くなる。
:::
::::

::::quiz{answer="0.5" tol="1e-6"}
$\Delta\nu=2\,\mathrm{GHz}$ の矩形帯域の第一零点の遅延を ns 単位で答えよ。

:::explanation
$\pi\Delta\nu\tau=\pi$ が最初の零点なので $\tau=1/(2\times10^9)\,\mathrm s=0.5\,\mathrm{ns}$。
:::
::::

::::quiz
単一モードの熱光の強度相関から一般に復元できないものを選べ。

- [x] 複素相互相関の位相
- [ ] 複素相互相関の絶対値
- [ ] 無相関時の基線の高さ

:::explanation
式 :eqref[eq-optics-siegert] は $\abs{\gamma}^2$ だけを含む。位相が異なっても同じ絶対値なら同じ強度相関になる。
:::
::::

::::quiz{answer="2" tol="1e-6"}
単一モードの熱光で同一点・零遅延の $g^{(2)}$ を答えよ。

:::explanation
自己相関の規格化で $\gamma(0)=1$ である。従って $g^{(2)}(0)=1+1=2$。
:::
::::

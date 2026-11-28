---
title: 偏光と結晶光学
description: 偏光の複素振幅と相関を区別し、波長板、結晶の複屈折および磁場による偏光回転を導く。
notation:
  - symbol: 'J'
    meaning: Jones ベクトル
  - symbol: 'S_0,S_1,S_2,S_3'
    meaning: Stokes パラメータ
  - symbol: '\eps_r'
    meaning: 比誘電率テンソル
  - symbol: 'n_o,n_e'
    meaning: 常光と主異常光の屈折率
---

## この章の目標

- 偏光楕円を二成分の振幅と位相から求める。
- Jones 行列と Mueller 行列を用いて素子を計算する。
- 波面法線を与えて結晶の固有波を求める。
- 偏光回転の往復を説明する。

## 偏光状態と Jones ベクトル

### 電場の軌跡
単色平面波は :eqref[eq-optics-plane] の横波条件を満たす。等方な媒質を $+z$ に進む波の固定位置で、全体の位相を時間原点へ移すと、実電場は $x=a\cos\omega t,\ y=b\cos(\omega t-\delta)$ となる。複素成分を並べた $J=(a,be^{i\delta})^\mathsf T$ を :term[Jones ベクトル]{en="Jones vector" yomi="じょーんずべくとる"} という。$a,b>0$ に対し $X=x/a,Y=y/b$ とすると $Y-X\cos\delta=\sin\omega t\sin\delta$ である。二乗して
$$
X^2+Y^2-2XY\cos\delta=\sin^2\delta
$$
を得る。$\sin\delta=0$ では直線となり、$a=b,\ \cos\delta=0$ なら円となる。一般の楕円の軸は、実対称行列 $\begin{pmatrix}a^2&ab\cos\delta\\ab\cos\delta&b^2\end{pmatrix}$ の固有ベクトルである。この行列は $(\cos\omega t,\sin\omega t)$ から電場を作る実行列とその転置の積なので、固有値が半軸長の二乗になる。

:::example[回転の向き]
$J=(1,i)^\mathsf T/\sqrt2$ では時間とともに電場が $+x$ から $+y$ へ回る。本章ではこれを正の円偏光とする。円偏光の呼称は観察する向きによって入れ替わるため、成分を明示する。
:::
:::example[楕円]
$J=(2,i)^\mathsf T$ なら軌跡は $x^2/4+y^2=1$ で、長軸は $x$ 方向である。
:::

:::example[楕円軸を対角化する]
$a=b=1$、$\delta=\pi/3$ なら二乗半軸の行列の対角は $1,1$、非対角は $1/2$ である。固有ベクトルは $(1,\pm1)/\sqrt2$、半軸長は $\sqrt{3/2},\sqrt{1/2}$。元の $x,y$ 振幅が等しくても円とは限らない。
:::

### 楕円を時間波形から再構成する

電場を $C(\cos\omega t,\sin\omega t)^{\mathsf T}$ と書くと、$C$ の特異値が楕円の半軸になる。これは $CC^{\mathsf T}$ の直交固有基底で、$C$ の二つの行が直交し、長さがそれぞれ固有値の平方根になるためである。$CC^{\mathsf T}$ の跡は $a^2+b^2$、行列式は $a^2b^2\sin^2\delta$ だから、固有方程式を解けば半軸を得る。半軸長が一致して方位が不定になる条件と、短半軸が消えて回転の向きが不定になる条件は異なる。長短半軸を $a^\prime\ge b^\prime$ とする。長半軸が正であれば、符号付き楕円率角 $\chi$ を $\tan\chi=\operatorname{sgn}(\sin\delta)b^\prime/a^\prime$、$-\pi/4\le\chi\le\pi/4$ で定める。

:::example[両振幅が違う直交位相]
二成分の位相差を四分の一周期とする。

入力は横電場の振幅 $a=3\,\mathrm{V/m}$、縦電場の振幅 $b=1\,\mathrm{V/m}$、相対位相 $\delta\simeq 1.5708\,\mathrm{rad}$ である。

行列の跡を取る。

$$
\begin{aligned}
T&=a^{2} + b^{2},\\
T/(\mathrm{V^2/m^2})&= \left(3\right)^{2} + \left(1\right)^{2}=10.
\end{aligned}
$$

特性多項式の判別式の平方根である。

$$
\begin{aligned}
\Delta&=\sqrt{4 a^{2} b^{2} \cos^{2}{\left(\delta \right)} + \left(a^{2} - b^{2}\right)^{2}},\\
\Delta/(\mathrm{V^2/m^2})&\simeq \sqrt{4 \cdot \left(3\right)^{2} \cdot \left(1\right)^{2} \cdot \cos^{2}{\left(\left(1.570796327\right) \right)} + \left(\left(3\right)^{2} - \left(1\right)^{2}\right)^{2}}=8.
\end{aligned}
$$

大きい固有値の平方根を取る。

$$
\begin{aligned}
a'&=\sqrt{\frac{\Delta}{2} + \frac{T}{2}},\\
a'/(\mathrm{V/m})&= \sqrt{\frac{\left(8\right)}{2} + \frac{\left(10\right)}{2}}=3.
\end{aligned}
$$

特異値の積が行列式の絶対値に等しいことを使う。

$$
\begin{aligned}
b'&=\frac{\left|{a b \sin{\left(\delta \right)}}\right|}{a'},\\
b'/(\mathrm{V/m})&\simeq \frac{\left|{\left(3\right) \cdot \left(1\right) \cdot \sin{\left(\left(1.570796327\right) \right)}}\right|}{\left(3\right)}=1.
\end{aligned}
$$

短半軸の比に回転方向の符号を与える。

$$
\begin{aligned}
\tan\chi&=\frac{a b \sin{\left(\delta \right)}}{(a')^{2}},\\
\tan\chi&\simeq \frac{\left(3\right) \cdot \left(1\right) \cdot \sin{\left(\left(1.570796327\right) \right)}}{\left(3\right)^{2}}=\frac{1}{3}.
\end{aligned}
$$

楕円の軸は測定軸と一致する。
:::

:::example[測定軸より傾く楕円]
二成分の振幅が等しく位相差が小さい状態を選ぶ。

入力は横電場の振幅 $a=2\,\mathrm{V/m}$、縦電場の振幅 $b=2\,\mathrm{V/m}$、相対位相 $\delta\simeq 0.5236\,\mathrm{rad}$ である。

行列の跡を取る。

$$
\begin{aligned}
T&=a^{2} + b^{2},\\
T/(\mathrm{V^2/m^2})&= \left(2\right)^{2} + \left(2\right)^{2}=8.
\end{aligned}
$$

特性多項式の判別式の平方根である。

$$
\begin{aligned}
\Delta&=\sqrt{4 a^{2} b^{2} \cos^{2}{\left(\delta \right)} + \left(a^{2} - b^{2}\right)^{2}},\\
\Delta/(\mathrm{V^2/m^2})&\simeq \sqrt{4 \cdot \left(2\right)^{2} \cdot \left(2\right)^{2} \cdot \cos^{2}{\left(\left(0.5235987756\right) \right)} + \left(\left(2\right)^{2} - \left(2\right)^{2}\right)^{2}}\simeq 6.9282.
\end{aligned}
$$

大きい固有値の平方根を取る。

$$
\begin{aligned}
a'&=\sqrt{\frac{\Delta}{2} + \frac{T}{2}},\\
a'/(\mathrm{V/m})&\simeq \sqrt{\frac{\left(6.92820323\right)}{2} + \frac{\left(8\right)}{2}}\simeq 2.7321.
\end{aligned}
$$

特異値の積が行列式の絶対値に等しいことを使う。

$$
\begin{aligned}
b'&=\frac{\left|{a b \sin{\left(\delta \right)}}\right|}{a'},\\
b'/(\mathrm{V/m})&\simeq \frac{\left|{\left(2\right) \cdot \left(2\right) \cdot \sin{\left(\left(0.5235987756\right) \right)}}\right|}{\left(2.732050808\right)}\simeq 0.73205.
\end{aligned}
$$

短半軸の比に回転方向の符号を与える。

$$
\begin{aligned}
\tan\chi&=\frac{a b \sin{\left(\delta \right)}}{(a')^{2}},\\
\tan\chi&\simeq \frac{\left(2\right) \cdot \left(2\right) \cdot \sin{\left(\left(0.5235987756\right) \right)}}{\left(2.732050808\right)^{2}}\simeq 0.26795.
\end{aligned}
$$

長軸は $45^\circ$ に傾く。成分振幅が等しくても長短は一致しない。
:::

:::example[回転方向を逆にする]
前例で電場の一方の位相を反転する。

入力は横電場の振幅 $a=2\,\mathrm{V/m}$、縦電場の振幅 $b=2\,\mathrm{V/m}$、相対位相 $\delta\simeq -0.5236\,\mathrm{rad}$ である。

行列の跡を取る。

$$
\begin{aligned}
T&=a^{2} + b^{2},\\
T/(\mathrm{V^2/m^2})&= \left(2\right)^{2} + \left(2\right)^{2}=8.
\end{aligned}
$$

特性多項式の判別式の平方根である。

$$
\begin{aligned}
\Delta&=\sqrt{4 a^{2} b^{2} \cos^{2}{\left(\delta \right)} + \left(a^{2} - b^{2}\right)^{2}},\\
\Delta/(\mathrm{V^2/m^2})&\simeq \sqrt{4 \cdot \left(2\right)^{2} \cdot \left(2\right)^{2} \cdot \cos^{2}{\left(\left(-0.5235987756\right) \right)} + \left(\left(2\right)^{2} - \left(2\right)^{2}\right)^{2}}\simeq 6.9282.
\end{aligned}
$$

大きい固有値の平方根を取る。

$$
\begin{aligned}
a'&=\sqrt{\frac{\Delta}{2} + \frac{T}{2}},\\
a'/(\mathrm{V/m})&\simeq \sqrt{\frac{\left(6.92820323\right)}{2} + \frac{\left(8\right)}{2}}\simeq 2.7321.
\end{aligned}
$$

特異値の積が行列式の絶対値に等しいことを使う。

$$
\begin{aligned}
b'&=\frac{\left|{a b \sin{\left(\delta \right)}}\right|}{a'},\\
b'/(\mathrm{V/m})&\simeq \frac{\left|{\left(2\right) \cdot \left(2\right) \cdot \sin{\left(\left(-0.5235987756\right) \right)}}\right|}{\left(2.732050808\right)}\simeq 0.73205.
\end{aligned}
$$

短半軸の比に回転方向の符号を与える。

$$
\begin{aligned}
\tan\chi&=\frac{a b \sin{\left(\delta \right)}}{(a')^{2}},\\
\tan\chi&\simeq \frac{\left(2\right) \cdot \left(2\right) \cdot \sin{\left(\left(-0.5235987756\right) \right)}}{\left(2.732050808\right)^{2}}\simeq -0.26795.
\end{aligned}
$$

軌跡の形は同じで、時間とともにたどる向きだけが逆になる。
:::

## 偏光素子と Jones 行列

線形素子が入出力の横波成分だけを結ぶとき、出力を $MJ$ と書く。この $M$ が :term[Jones 行列]{en="Jones matrix" yomi="じょーんずぎょうれつ"} である。軸 $\vb u=(\cos\theta,\sin\theta)^\mathsf T$ の理想偏光子は $M=\vb u\vb u^\mathsf T$ である。直線偏光と軸の差を $\alpha$ とすると振幅の射影が $\cos\alpha$、強度比が $\cos^2\alpha$ になる。これが Malus の法則である。

反対に、二つの直交成分を消さず、位相差だけを与える素子を考える。厚さ $d$ の無反射の板の両軸で屈折率が $n_u,n_v$ なら、相対位相は $\delta=k_0d(n_v-n_u)$ である。軸が角 $\theta$ をなす場合の行列は $R(\theta)\diag(1,e^{i\delta})R(-\theta)$ であり、$R$ は実二次元回転である。ここでは最初に通る素子を積の右に置く。

:::example[四分の一波長板]
$(1,1)^\mathsf T/\sqrt2$ に $\diag(1,i)$ を作用させると正の円偏光になる。入射が軸と平行なら相対位相は観測できない。
:::
:::example[半波長板]
$\diag(1,-1)$ は軸座標で直線の角を $\alpha$ から $-\alpha$ に移す。実験室での板の軸角が $\theta$、入射が $\psi$ なら、出射角は $2\theta-\psi$ になる。
:::

### 素子の交換で出力が変わる理由
軸 $\theta$ の偏光子は階数 1 なので、直前の素子で偏光が回ると透過強度も変わる。波長板はユニタリであり、出力強度を変えないが成分の相対位相を変える。二素子の可換性は、どちらも同じ直交固有軸を持つ場合には成り立つ。しかし軸が違う場合、一般には固有軸の一方を他方が混ぜるので可換でない。

:::example[板の前後の偏光子]
$x$ 入力に、軸角 $45^\circ$ の半波長板を作用させると $y$ 出力になる。後段の $y$ 偏光子なら全透過する。一方、この偏光子を先に置けば $x$ 入力を全部遮り、板には光が届かない。
:::

### 固有軸が異なる遅相器を組み合わせる
すべての無損失 Jones 行列は共通位相を除けば行列式 1 のユニタリ行列になる。二固有値は単位円上で互いに逆数なので、ある直交する二固有偏光の上で $e^{\pm i\delta/2}$ と書ける。単一の波長板では固有偏光は直線であるが、軸の異なる板を組み合わせると楕円の固有偏光も得られる。例えば $\diag(1,i)$ を $Q$、$45^\circ$ の板を $Q_{45}$ と置く。$Q_{45}=R(\pi/4)Q R(-\pi/4)$ では非対角が $(1-i)/2$ でゼロではない。従って $Q Q_{45}$ と $Q_{45}Q$ の非対角は $i$ の因子が付く位置が異なり、二つの配置は一般に別の素子になる。

### 平行検光子と直交検光子で遅相を測る

$x$ 直線を、軸角 $\theta$ の無損失遅相板へ入れる。固有軸を $\vb u=(c,s)$、$\vb v=(-s,c)$ とすれば、出力は $c\vb u-e^{i\delta}s\vb v$。$x$ 検光子には $c^2+e^{i\delta}s^2$、$y$ 検光子には $cs(1-e^{i\delta})$ が入る。従って直交配置での透過率は $\sin^2(2\theta)\sin^2(\delta/2)$。二測定の和はユニタリ性により一である。直交強度だけでは $\delta$ と $-\delta$ が区別できないため、遅相の符号を求めるなら円検光も加える。

:::example[半波長板を測る]
光軸を入射の直線から $45^\circ$ 傾ける。

入力は板の軸角 $\theta\simeq 0.7854\,\mathrm{rad}$、相対遅相 $\delta\simeq 3.1416\,\mathrm{rad}$、入射パワー $P=\frac{1}{250}\,\mathrm W$ である。

二軸に等しく入ったときに最も遅相が効く。

$$
\begin{aligned}
g&=\sin^{2}{\left(2 \theta \right)},\\
g&\simeq \sin^{2}{\left(2 \cdot \left(0.7853981634\right) \right)}=1.
\end{aligned}
$$

固有成分を直交検光したときの位相差項である。

$$
\begin{aligned}
h&=\sin^{2}{\left(\frac{\delta}{2} \right)},\\
h&\simeq \sin^{2}{\left(\frac{\left(3.141592654\right)}{2} \right)}=1.
\end{aligned}
$$

軸角と位相の二つの因子を掛ける。

$$
\begin{aligned}
T_\perp&=g h,\\
T_\perp&= \left(1\right) \cdot \left(1\right)=1.
\end{aligned}
$$

板から出た光を $y$ へ射影する。

$$
\begin{aligned}
P_y&=P T_\perp,\\
P_y/(\mathrm W)&= \left(\frac{1}{250}\right) \cdot \left(1\right)=\frac{1}{250}.
\end{aligned}
$$

直交する二出力の合計パワーはもとの入射に等しい。

$$
\begin{aligned}
P_x&=P - P_y,\\
P_x/(\mathrm W)&= \left(\frac{1}{250}\right) - \left(\frac{1}{250}\right)=0.
\end{aligned}
$$

出力の直線は $y$ 方向である。
:::

:::example[四分の一波長板を測る]
軸角を保って遅相を半分にする。

入力は板の軸角 $\theta\simeq 0.7854\,\mathrm{rad}$、相対遅相 $\delta\simeq 1.5708\,\mathrm{rad}$、入射パワー $P=\frac{1}{250}\,\mathrm W$ である。

二軸に等しく入ったときに最も遅相が効く。

$$
\begin{aligned}
g&=\sin^{2}{\left(2 \theta \right)},\\
g&\simeq \sin^{2}{\left(2 \cdot \left(0.7853981634\right) \right)}=1.
\end{aligned}
$$

固有成分を直交検光したときの位相差項である。

$$
\begin{aligned}
h&=\sin^{2}{\left(\frac{\delta}{2} \right)},\\
h&\simeq \sin^{2}{\left(\frac{\left(1.570796327\right)}{2} \right)}=\frac{1}{2}.
\end{aligned}
$$

軸角と位相の二つの因子を掛ける。

$$
\begin{aligned}
T_\perp&=g h,\\
T_\perp&= \left(1\right) \cdot \left(\frac{1}{2}\right)=\frac{1}{2}.
\end{aligned}
$$

板から出た光を $y$ へ射影する。

$$
\begin{aligned}
P_y&=P T_\perp,\\
P_y/(\mathrm W)&= \left(\frac{1}{250}\right) \cdot \left(\frac{1}{2}\right)=\frac{1}{500}.
\end{aligned}
$$

直交する二出力の合計パワーはもとの入射に等しい。

$$
\begin{aligned}
P_x&=P - P_y,\\
P_x/(\mathrm W)&= \left(\frac{1}{250}\right) - \left(\frac{1}{500}\right)=\frac{1}{500}.
\end{aligned}
$$

両出力は等しい。円であることの確認には円検光を用いる。
:::

:::example[軸を固定してはいけない例]
半波長の遅相で軸を入射方向に合わせる。

入力は板の軸角 $\theta=0\,\mathrm{rad}$、相対遅相 $\delta\simeq 3.1416\,\mathrm{rad}$、入射パワー $P=\frac{1}{250}\,\mathrm W$ である。

二軸に等しく入ったときに最も遅相が効く。

$$
\begin{aligned}
g&=\sin^{2}{\left(2 \theta \right)},\\
g&= \sin^{2}{\left(2 \cdot \left(0\right) \right)}=0.
\end{aligned}
$$

固有成分を直交検光したときの位相差項である。

$$
\begin{aligned}
h&=\sin^{2}{\left(\frac{\delta}{2} \right)},\\
h&\simeq \sin^{2}{\left(\frac{\left(3.141592654\right)}{2} \right)}=1.
\end{aligned}
$$

軸角と位相の二つの因子を掛ける。

$$
\begin{aligned}
T_\perp&=g h,\\
T_\perp&= \left(0\right) \cdot \left(1\right)=0.
\end{aligned}
$$

板から出た光を $y$ へ射影する。

$$
\begin{aligned}
P_y&=P T_\perp,\\
P_y/(\mathrm W)&= \left(\frac{1}{250}\right) \cdot \left(0\right)=0.
\end{aligned}
$$

直交する二出力の合計パワーはもとの入射に等しい。

$$
\begin{aligned}
P_x&=P - P_y,\\
P_x/(\mathrm W)&= \left(\frac{1}{250}\right) - \left(0\right)=\frac{1}{250}.
\end{aligned}
$$

直交出力は消える。遅相が無い場合との区別には軸角を変える。
:::

## Stokes パラメータと Poincaré 球

準単色の光を考え、$J$ を横電場の複素包絡線とする。光周期より長い時間平均または集合平均を用い、$C=\langle JJ^\dagger\rangle$ を相関行列とする。$J$ が電場の単位を持つとき、以下の強度は共通の換算係数を除いた電場二乗の尺度である。放射照度へ換算するときは、同じ媒質を同じ方向に進む成分に $1/(2\eta)$ を掛ける。総強度 $S_0$ と、直交する三対の検光子で得られる強度差 $S_1,S_2,S_3$ が :term[Stokes パラメータ]{en="Stokes parameters" yomi="すとーくすぱらめーた"} である。対は $x/y$、$45^\circ/135^\circ$、正負の円偏光とする。例えば円偏光の強度差は、$J$ に $(1,\pm i)/\sqrt2$ の共役転置を掛けて絶対値を二乗すると $2\Im C_{yx}$ となる。従って
$$
C=\frac12\begin{pmatrix}S_0+S_1&S_2-iS_3\\ S_2+iS_3&S_0-S_1\end{pmatrix}.
\label{eq-optics-stokes}
$$
$C$ は半正定値で、その固有値は $(S_0\pm|\vb S|)/2$、$\vb S=(S_1,S_2,S_3)$ である。実際、跡は $S_0$ で行列式は $(S_0^2-|\vb S|^2)/4$ となる。完全偏光では階数 1 なので $|\vb S|=S_0$ である。$S_0>0$ での $\vb S/S_0$ の単位球を :term[Poincaré 球]{en="Poincaré sphere" yomi="ぽあんかれきゅう"} という。

:::example[赤道と極]
直線 $(\cos\psi,\sin\psi)$ は $(\cos2\psi,\sin2\psi,0)$、正円偏光は $(0,0,1)$ に対応する。直線偏光では $\psi$ と $\psi+\pi$ が同一状態であることが二倍角に現れる。
:::
:::example[偏光角の推定]
規格化した三成分が $(0,1,0)$ なら $2\psi=\pi/2$ で、直線は $45^\circ$ を向く。
:::

### 楕円の方位と楕円率
完全偏光の長短半軸を $a'\ge b'\ge0$ とし、長軸方位を $\psi$、回転の符号を $\operatorname{sgn}S_3$ とする。長軸座標の Jones ベクトルは共通位相を除き $(a',\pm i b')$ に取れる。そこで $\tan\chi=\pm b'/a'$ と置くと、$S_1/S_0=\cos2\chi$、$S_3/S_0=\sin2\chi$、$S_2=0$。実座標を $\psi$ 回転して相関行列へ代入すると $(S_1,S_2,S_3)/S_0=(\cos2\chi\cos2\psi,\cos2\chi\sin2\psi,\sin2\chi)$ となる。球の極では方位を決められない。円の長軸が定まらないことに対応する。

### Stokes 空間の回転
$C'=MCM^\dagger$ で $M$ がユニタリなら、$\tr C$ と $\det C$ が不変である。式 :eqref[eq-optics-stokes] より $S_0$ と $\abs{\vb S}$ も不変。作用は $\vb S$ に線形で、任意の $\vb S$ が十分大きい $S_0$ と組み合わせられるため、三次元空間で直交変換となる。ユニタリ行列は固有位相を連続にゼロへ動かせるので、その直交変換の行列式は $+1$。よって Poincaré 球上の回転である。例えば軸の合った波長板 $M=\diag(1,e^{i\delta})$ では $S_2+iS_3$ が $e^{i\delta}$ 倍となり、$S_1$ 軸回りに角 $\delta$ だけ回転する。

:::example[三組の検光子で相関を復元する]
検光子出力の共通規格で $I_x=3,I_y=1,I_{45^\circ}=2.5,I_{135^\circ}=1.5$、円偏光の二出力が $2,2$ とする。差から $(S_0,S_1,S_2,S_3)=(4,2,1,0)$、式 :eqref[eq-optics-stokes] から $C_{xy}=1/2$ を得る。偏光度は $\sqrt5/4$。方位は $2\psi=\arctan(1/2)$ であり、完全偏光部分の直線は $x$ から約 $13.3^\circ$ 傾く。
:::

### 二つの偏光の球面上の距離と検光

規格化した純粋 Jones 状態 $J,K$ の相関行列は階数一の射影である。ここで用いる Pauli 行列の順序は $\sigma_1=\diag(1,-1)$、$\sigma_2=\begin{pmatrix}0&1\\1&0\end{pmatrix}$、$\sigma_3=\begin{pmatrix}0&-i\\i&0\end{pmatrix}$ とし、$\tr(\sigma_j)=0$、$\tr(\sigma_j\sigma_k)=2\delta_{jk}$ は行列の積から確かめられる。:eqref[eq-optics-stokes] をそれぞれ $(I+\vb s\cdot\vb\sigma)/2$ と $(I+\vb t\cdot\vb\sigma)/2$ に展開する。積の跡は一方で $|J^\dagger K|^2$、他方で Pauli 行列の跡の直交性から $(1+\vb s\cdot\vb t)/2$。球の中心角を $\alpha$ とすれば透過率は $\cos^2(\alpha/2)$ である。従って球面上の対蹠点は直交した電場を表す。一方、Stokes ベクトルの直交は電場の直交に対応せず、透過率が二分の一になる。

:::example[直線検光子に直交する入力]
横直線の入力を縦直線の検光子で測る。

入力は入射の楕円方位 $\psi=0\,\mathrm{rad}$、入射の符号付き楕円率角 $\chi=0\,\mathrm{rad}$、検光子の方位 $\psi_{\rm a}\simeq 1.5708\,\mathrm{rad}$、検光子の楕円率角 $\chi_{\rm a}=0\,\mathrm{rad}$ である。

二つの赤道射影の方位角の差は実空間の方位角の差の二倍である。

$$
\begin{aligned}
h&=\cos{\left(2 \chi \right)} \cos{\left(2 \chi_{\rm a} \right)} \cos{\left(2 \psi - 2 \psi_{\rm a} \right)},\\
h&\simeq \cos{\left(2 \cdot \left(0\right) \right)} \cdot \cos{\left(2 \cdot \left(0\right) \right)} \cdot \cos{\left(2 \cdot \left(0\right) - 2 \cdot \left(1.570796327\right) \right)}=-1.
\end{aligned}
$$

円成分の強度差を表す縦成分の積を求める。

$$
\begin{aligned}
g&=\sin{\left(2 \chi \right)} \sin{\left(2 \chi_{\rm a} \right)},\\
g&= \sin{\left(2 \cdot \left(0\right) \right)} \cdot \sin{\left(2 \cdot \left(0\right) \right)}=0.
\end{aligned}
$$

球の水平成分と縦成分を足す。

$$
\begin{aligned}
\vb s\cdot\vb t&=g + h,\\
\vb s\cdot\vb t&= \left(0\right) + \left(-1\right)=-1.
\end{aligned}
$$

二状態が球面上でどれだけ離れているかを角で表す。

$$
\begin{aligned}
\alpha&=\operatorname{acos}{\left(\vb s\cdot\vb t \right)},\\
\alpha/(\mathrm{rad})&= \operatorname{acos}{\left(\left(-1\right) \right)}\simeq 3.1416.
\end{aligned}
$$

Jones 射影の二乗を球上の内積から求める。

$$
\begin{aligned}
T&=\frac{\vb s\cdot\vb t}{2} + \frac{1}{2},\\
T&= \frac{\left(-1\right)}{2} + \frac{1}{2}=0.
\end{aligned}
$$

赤道の対蹠点なので透過しない。
:::

:::example[円入力を直線で検光]
円偏光を斜めの直線検光子に入れる。

入力は入射の楕円方位 $\psi=0\,\mathrm{rad}$、入射の符号付き楕円率角 $\chi=\pi/4\,\mathrm{rad}$、検光子の方位 $\psi_{\rm a}\simeq 0.4488\,\mathrm{rad}$、検光子の楕円率角 $\chi_{\rm a}=0\,\mathrm{rad}$ である。

二つの赤道射影の方位角の差は実空間の方位角の差の二倍である。

$$
\begin{aligned}
h&=\cos{\left(2 \chi \right)} \cos{\left(2 \chi_{\rm a} \right)} \cos{\left(2 \psi - 2 \psi_{\rm a} \right)},\\
h&\simeq \cos{\left(2 \cdot \left(\pi/4\right) \right)} \cdot \cos{\left(2 \cdot \left(0\right) \right)} \cdot \cos{\left(2 \cdot \left(0\right) - 2 \cdot \left(0.4487989505\right) \right)}=0.
\end{aligned}
$$

円成分の強度差を表す縦成分の積を求める。

$$
\begin{aligned}
g&=\sin{\left(2 \chi \right)} \sin{\left(2 \chi_{\rm a} \right)},\\
g&\simeq \sin{\left(2 \cdot \left(\pi/4\right) \right)} \cdot \sin{\left(2 \cdot \left(0\right) \right)}=0.
\end{aligned}
$$

球の水平成分と縦成分を足す。

$$
\begin{aligned}
\vb s\cdot\vb t&=g + h,\\
\vb s\cdot\vb t&= \left(0\right) + \left(0\right)=0.
\end{aligned}
$$

二状態が球面上でどれだけ離れているかを角で表す。

$$
\begin{aligned}
\alpha&=\operatorname{acos}{\left(\vb s\cdot\vb t \right)},\\
\alpha/(\mathrm{rad})&= \operatorname{acos}{\left(\left(0\right) \right)}\simeq 1.5708.
\end{aligned}
$$

Jones 射影の二乗を球上の内積から求める。

$$
\begin{aligned}
T&=\frac{\vb s\cdot\vb t}{2} + \frac{1}{2},\\
T&= \frac{\left(0\right)}{2} + \frac{1}{2}=\frac{1}{2}.
\end{aligned}
$$

方位によらず半分が出る。極と赤道の中心角は直角である。
:::

:::example[同じ回転向きの楕円を読む]
円検光子で正の楕円率を持つ入力を測る。

入力は入射の楕円方位 $\psi=0\,\mathrm{rad}$、入射の符号付き楕円率角 $\chi\simeq 0.2618\,\mathrm{rad}$、検光子の方位 $\psi_{\rm a}=0\,\mathrm{rad}$、検光子の楕円率角 $\chi_{\rm a}=\pi/4\,\mathrm{rad}$ である。

二つの赤道射影の方位角の差は実空間の方位角の差の二倍である。

$$
\begin{aligned}
h&=\cos{\left(2 \chi \right)} \cos{\left(2 \chi_{\rm a} \right)} \cos{\left(2 \psi - 2 \psi_{\rm a} \right)},\\
h&\simeq \cos{\left(2 \cdot \left(0.2617993878\right) \right)} \cdot \cos{\left(2 \cdot \left(\pi/4\right) \right)} \cdot \cos{\left(2 \cdot \left(0\right) - 2 \cdot \left(0\right) \right)}=0.
\end{aligned}
$$

円成分の強度差を表す縦成分の積を求める。

$$
\begin{aligned}
g&=\sin{\left(2 \chi \right)} \sin{\left(2 \chi_{\rm a} \right)},\\
g&\simeq \sin{\left(2 \cdot \left(0.2617993878\right) \right)} \cdot \sin{\left(2 \cdot \left(\pi/4\right) \right)}=\frac{1}{2}.
\end{aligned}
$$

球の水平成分と縦成分を足す。

$$
\begin{aligned}
\vb s\cdot\vb t&=g + h,\\
\vb s\cdot\vb t&= \left(\frac{1}{2}\right) + \left(0\right)=\frac{1}{2}.
\end{aligned}
$$

二状態が球面上でどれだけ離れているかを角で表す。

$$
\begin{aligned}
\alpha&=\operatorname{acos}{\left(\vb s\cdot\vb t \right)},\\
\alpha/(\mathrm{rad})&= \operatorname{acos}{\left(\left(\frac{1}{2}\right) \right)}\simeq 1.0472.
\end{aligned}
$$

Jones 射影の二乗を球上の内積から求める。

$$
\begin{aligned}
T&=\frac{\vb s\cdot\vb t}{2} + \frac{1}{2},\\
T&= \frac{\left(\frac{1}{2}\right)}{2} + \frac{1}{2}=\frac{3}{4}.
\end{aligned}
$$

正の円検光子へ半分より多くが出る。
:::

## 部分偏光と Mueller 行列

$C$ の小さい固有値を $\lambda_-$ とすると、$C=\lambda_-I+(C-\lambda_-I)$ である。後項は半正定値で階数 1 以下なので、非偏光成分と完全偏光成分に分けられる。完全偏光部分の強度の割合が偏光度 $|\vb S|/S_0$ である。この分解の一意性は、非偏光とした行列を引いた残りにゼロ固有値が必要なことによる。

:eqref[eq-optics-stokes] を $C=\frac12\sum_\nu S_\nu\sigma_\nu$ と書く。$\sigma_0=I$、$\sigma_1=\diag(1,-1)$、$\sigma_2$ は実非対角成分が 1、$\sigma_3$ は右上 $-i$、左下 $i$ の行列である。直接計算から $\tr(\sigma_\mu\sigma_\nu)=2\delta_{\mu\nu}$ である。出力 $C'=MCM^\dagger$ に $\sigma_\mu$ を掛けて跡を取れば
$$
S'_\mu=\sum_\nu \mathcal M_{\mu\nu}S_\nu,\qquad
\mathcal M_{\mu\nu}=\tfrac12\tr(\sigma_\mu M\sigma_\nu M^\dagger)
$$
となる。これが :term[Mueller 行列]{en="Mueller matrix" yomi="みゅらーぎょうれつ"} である。素子が時間で変動するとき、入力と素子の変動が独立なら Mueller 行列の平均を取れる。Jones 行列を振幅のまま平均する計算とは異なる。

:::example[非偏光から直線偏光へ]
$C=I/2$ を $x$ 偏光子に入れると $C'=\diag(1/2,0)$ である。総強度は半分で偏光度は 1 になる。
:::
:::example[等強度の独立成分]
互いに独立な $x,y$ 直線に同じ強度を与えると $C$ は単位行列に比例する。従って検光子を回しても透過強度が変わらない。
:::

### 位相差がランダムな素子
$M(\delta)=\diag(1,e^{i\delta})$ に対して $C'_{yx}=e^{i\delta}C_{yx}$、対角は不変である。$\delta$ を独立に平均すれば、対角を保ったまま非対角だけに $\langle e^{i\delta}\rangle$ が掛かる。これは一つの実現では強度も偏光度も変えない板が、平均では偏光度を下げ得る機構である。

:::example[偏光解消の具体例]
$\delta=0,\pi$ が等確率なら、入力 $J=(1,1)/\sqrt2$ の平均出力は $I/2$ となる。単純に Jones 行列を平均した $\diag(1,0)$ を作用させると偏光した半強度の出力となるため、統計平均の計算では誤る。
:::

### 不完全な検光子の校正
直交する二軸の振幅透過率を正実数 $t_x,t_y$ とすると、$C'_{xx}=t_x^2C_{xx}$、$C'_{yy}=t_y^2C_{yy}$、$C'_{yx}=t_xt_yC_{yx}$。$S_0,S_1$ は行列 $\frac12\begin{pmatrix}t_x^2+t_y^2&t_x^2-t_y^2\\t_x^2-t_y^2&t_x^2+t_y^2\end{pmatrix}$ で変換され、$S_2,S_3$ は $t_xt_y$ 倍である。両軸の透過が異なると非偏光からも $S_1$ が生じる。

:::example[有限の漏れ]
軸ごとのパワー透過が $1,0.01$ の検光子に非偏光を入れると規格化強度は $0.505$、偏光度は $0.99/1.01\simeq0.980$。漏れが小さくても理想偏光子の値 1 とは異なる。
:::

### 二色性素子を通る部分偏光

直線固有軸でのパワー透過率を $T_x,T_y$ とし、Jones 行列を $\diag(\sqrt{T_x},\sqrt{T_y}e^{i\delta})$ とする。対角成分の和と差に :eqref[eq-optics-stokes] を入れると $S'_0=aS_0+bS_1$、$S'_1=bS_0+aS_1$、$a=(T_x+T_y)/2,b=(T_x-T_y)/2$。非対角成分から $S_2'+iS_3'=\sqrt{T_xT_y}e^{i\delta}(S_2+iS_3)$ を得る。従って $(S_0')^2-|\vb S'|^2=T_xT_y(S_0^2-|\vb S|^2)$ である。以下の例で $s_j=S_j/S_0$、$s_j'=S_j'/S_0$ とし、出力の成分も入射の総強度で規格化する。$s_0^\prime$ が正であれば、三成分のノルムを $s_0^\prime$ で割って出力の偏光度を得る。損失の偏光差は偏光度を変えるが、固定 Jones 行列が生き残った純粋入力の純粋性を失わせることはない。

:::example[無偏光に差が生まれる]
二軸で透過が違う素子へ無偏光を入れる。

入力は横透過率 $T_x=\frac{9}{10}$、縦透過率 $T_y=\frac{1}{10}$、入力の $S_1/S_0$ $s_1=0$、入力の $S_2/S_0$ $s_2=0$ である。

無偏光入力には平均の透過率が効く。

$$
\begin{aligned}
a&=\frac{T_x}{2} + \frac{T_y}{2},\\
a&= \frac{\left(\frac{9}{10}\right)}{2} + \frac{\left(\frac{1}{10}\right)}{2}=\frac{1}{2}.
\end{aligned}
$$

透過率の非対称な寄与を取る。

$$
\begin{aligned}
b&=\frac{T_x}{2} - \frac{T_y}{2},\\
b&= \frac{\left(\frac{9}{10}\right)}{2} - \frac{\left(\frac{1}{10}\right)}{2}=\frac{2}{5}.
\end{aligned}
$$

強く透過する軸の入力に応じて増減する。

$$
\begin{aligned}
s'_0&=a + b s_1,\\
s'_0&= \left(\frac{1}{2}\right) + \left(\frac{2}{5}\right) \cdot \left(0\right)=\frac{1}{2}.
\end{aligned}
$$

入力が無偏光でも $b$ に応じて差が生まれる。

$$
\begin{aligned}
s'_1&=a s_1 + b,\\
s'_1&= \left(\frac{1}{2}\right) \cdot \left(0\right) + \left(\frac{2}{5}\right)=\frac{2}{5}.
\end{aligned}
$$

以下では遅相と入力の円成分をゼロとする。

$$
\begin{aligned}
s'_2&=s_2 \sqrt{T_x T_y},\\
s'_2&= \left(0\right) \cdot \sqrt{\left(\frac{9}{10}\right) \cdot \left(\frac{1}{10}\right)}=0.
\end{aligned}
$$

球内の半径を、出力の総強度で規格化して求める。

$$
\begin{aligned}
\mathcal P'&=\frac{\sqrt{(s'_1)^{2} + (s'_2)^{2}}}{s'_0},\\
\mathcal P'&= \frac{\sqrt{\left(\frac{2}{5}\right)^{2} + \left(0\right)^{2}}}{\left(\frac{1}{2}\right)}=\frac{4}{5}.
\end{aligned}
$$

入射では相関がなくても出射の偏光度は正になる。
:::

:::example[部分偏光と選択軸を合わせる]
入射の偏光が透過の大きい軸へ傾くとする。

入力は横透過率 $T_x=\frac{9}{10}$、縦透過率 $T_y=\frac{1}{10}$、入力の $S_1/S_0$ $s_1=\frac{1}{2}$、入力の $S_2/S_0$ $s_2=\frac{3}{10}$ である。

無偏光入力には平均の透過率が効く。

$$
\begin{aligned}
a&=\frac{T_x}{2} + \frac{T_y}{2},\\
a&= \frac{\left(\frac{9}{10}\right)}{2} + \frac{\left(\frac{1}{10}\right)}{2}=\frac{1}{2}.
\end{aligned}
$$

透過率の非対称な寄与を取る。

$$
\begin{aligned}
b&=\frac{T_x}{2} - \frac{T_y}{2},\\
b&= \frac{\left(\frac{9}{10}\right)}{2} - \frac{\left(\frac{1}{10}\right)}{2}=\frac{2}{5}.
\end{aligned}
$$

強く透過する軸の入力に応じて増減する。

$$
\begin{aligned}
s'_0&=a + b s_1,\\
s'_0&= \left(\frac{1}{2}\right) + \left(\frac{2}{5}\right) \cdot \left(\frac{1}{2}\right)=\frac{7}{10}.
\end{aligned}
$$

入力が無偏光でも $b$ に応じて差が生まれる。

$$
\begin{aligned}
s'_1&=a s_1 + b,\\
s'_1&= \left(\frac{1}{2}\right) \cdot \left(\frac{1}{2}\right) + \left(\frac{2}{5}\right)=\frac{13}{20}.
\end{aligned}
$$

以下では遅相と入力の円成分をゼロとする。

$$
\begin{aligned}
s'_2&=s_2 \sqrt{T_x T_y},\\
s'_2&= \left(\frac{3}{10}\right) \cdot \sqrt{\left(\frac{9}{10}\right) \cdot \left(\frac{1}{10}\right)}=\frac{9}{100}.
\end{aligned}
$$

球内の半径を、出力の総強度で規格化して求める。

$$
\begin{aligned}
\mathcal P'&=\frac{\sqrt{(s'_1)^{2} + (s'_2)^{2}}}{s'_0},\\
\mathcal P'&= \frac{\sqrt{\left(\frac{13}{20}\right)^{2} + \left(\frac{9}{100}\right)^{2}}}{\left(\frac{7}{10}\right)}\simeq 0.93743.
\end{aligned}
$$

出力の横の差成分が強まり、第二成分も残る。
:::

:::example[純粋状態を減衰させる]
同じ素子で直線偏光を減衰させる。

入力は横透過率 $T_x=\frac{9}{10}$、縦透過率 $T_y=\frac{1}{10}$、入力の $S_1/S_0$ $s_1=\frac{3}{5}$、入力の $S_2/S_0$ $s_2=\frac{4}{5}$ である。

無偏光入力には平均の透過率が効く。

$$
\begin{aligned}
a&=\frac{T_x}{2} + \frac{T_y}{2},\\
a&= \frac{\left(\frac{9}{10}\right)}{2} + \frac{\left(\frac{1}{10}\right)}{2}=\frac{1}{2}.
\end{aligned}
$$

透過率の非対称な寄与を取る。

$$
\begin{aligned}
b&=\frac{T_x}{2} - \frac{T_y}{2},\\
b&= \frac{\left(\frac{9}{10}\right)}{2} - \frac{\left(\frac{1}{10}\right)}{2}=\frac{2}{5}.
\end{aligned}
$$

強く透過する軸の入力に応じて増減する。

$$
\begin{aligned}
s'_0&=a + b s_1,\\
s'_0&= \left(\frac{1}{2}\right) + \left(\frac{2}{5}\right) \cdot \left(\frac{3}{5}\right)=\frac{37}{50}.
\end{aligned}
$$

入力が無偏光でも $b$ に応じて差が生まれる。

$$
\begin{aligned}
s'_1&=a s_1 + b,\\
s'_1&= \left(\frac{1}{2}\right) \cdot \left(\frac{3}{5}\right) + \left(\frac{2}{5}\right)=\frac{7}{10}.
\end{aligned}
$$

以下では遅相と入力の円成分をゼロとする。

$$
\begin{aligned}
s'_2&=s_2 \sqrt{T_x T_y},\\
s'_2&= \left(\frac{4}{5}\right) \cdot \sqrt{\left(\frac{9}{10}\right) \cdot \left(\frac{1}{10}\right)}=\frac{6}{25}.
\end{aligned}
$$

球内の半径を、出力の総強度で規格化して求める。

$$
\begin{aligned}
\mathcal P'&=\frac{\sqrt{(s'_1)^{2} + (s'_2)^{2}}}{s'_0},\\
\mathcal P'&= \frac{\sqrt{\left(\frac{7}{10}\right)^{2} + \left(\frac{6}{25}\right)^{2}}}{\left(\frac{37}{50}\right)}=1.
\end{aligned}
$$

横と縦の釣合いは変わるが、偏光度は一のままである。
:::

## 異方性媒質と屈折率楕円体

正定値な実対称誘電率 $\eps_r$ を持つ無損失非磁性結晶を考える。波数方向 $\hat{\vb s}$ に垂直な射影を $P_\perp=I-\hat{\vb s}\hat{\vb s}^\mathsf T$ とする。Maxwell 方程式の回転を平面波に適用すると $n^2P_\perp\vb E=\vb D/\eps_0$、発散式から $\vb D\perp\hat{\vb s}$ となる。$\vb E=\eps_r^{-1}\vb D/\eps_0$ を入れると
$$
P_\perp\eps_r^{-1}P_\perp\vb D=n^{-2}\vb D .
\label{eq-optics-indicatrix}
$$
この横平面での固有値問題を、:term[屈折率楕円体]{en="index ellipsoid" yomi="くっせつりつだえんたい"} $\vb x^\mathsf T\eps_r^{-1}\vb x=1$ の中心断面として解釈できる。固有基底で断面を表すと $x_1^2/n_1^2+x_2^2/n_2^2=1$ なので、半軸長が屈折率、方向が電束密度の向きである。電場は一般に $\hat{\vb s}$ と垂直にならない。

:::example[主軸]
$\eps_r=\diag(4,9,16)$ で $z$ へ進む波の横断面は $x^2/4+y^2/9=1$。固有屈折率は $2,3$ である。
:::
:::example[等方媒質]
$\eps_r=n^2I$ なら断面はどの向きでも円である。二偏光の位相速度は一致する。
:::

### 二軸性で円断面を探す
主屈折率が $n_x<n_y<n_z$ とする。$xz$ 面内の伝播方向 $\hat{\vb s}=(\sin\theta,0,\cos\theta)$ の横断面は二直交基底 $\hat{\vb y},(\cos\theta,0,-\sin\theta)$ で対角になる。円断面の条件は $n_y^{-2}=n_x^{-2}\cos^2\theta+n_z^{-2}\sin^2\theta$。右辺は $\sin^2\theta$ の一次式なので、$n_y^{-2}$ が二端の間にあると鏡映に対応する二つの伝播軸が得られる。

ほかの方向に円断面がないことも確かめる。$\eps_r^{-1}$ の対角要素を $a>b>c$ とする。任意の二次元断面は $xy$ 面、$yz$ 面とそれぞれ少なくとも一本の直線で交わる。従って断面上の Rayleigh 商の最大値は $b$ 以上、最小値は $b$ 以下である。円断面で両値が一致するなら値は $b$。その断面上の任意の点は $(a-b)x^2+(c-b)z^2=0$ を満たす。これは $(\sqrt{a-b}x-\sqrt{b-c}z)(\sqrt{a-b}x+\sqrt{b-c}z)=0$ と因数分解できる。二つの平面の和集合に含まれる平面はどちらかに一致する。実際、どちらにも含まれなければ、一方にだけ属するベクトルと他方にだけ属するベクトルを足したとき、両方の条件が破れる。従って円断面の法線は上で求めた二方向に限られる。

:::example[二つの軸の方向]
$n_x=\sqrt2,n_y=\sqrt{8/3},n_z=2$ とすると逆二乗は $1/2,3/8,1/4$ である。円断面には $\sin^2\theta=1/2$ が必要なので、軸は $z$ から $x$ 方向へ $\pm45^\circ$ にある。
:::

### 任意方向の楕円断面を跡と行列式で解く

$A=\eps_r^{-1}=\diag(a,b,c)$ とし、単位波面法線を $\vb s=(l,m,n)$ とする。横平面の実直交基底を行列 $Q$ の二列に並べると、逆二乗屈折率は $Q^{\mathsf T}A Q$ の固有値である。跡は $a+b+c-\vb s^{\mathsf T}A\vb s$。行列式は $bc\,l^2+ac\,m^2+ab\,n^2$ である。後者は $Q^{\mathsf T}A Q$ の行列式を二次の小行列で展開すると、$Q$ の第 $i$、$j$ 行の小行列式の二乗に対応する対角積が掛かり、小行列式は二基底の外積の成分になるためである。従って断面の二屈折率を求めるだけなら、具体的な横基底を構成する必要がない。

:::example[主軸方向で検算]
第三主軸へ進む波の横断面を解く。

入力は第一主軸の逆二乗屈折率 $a=\frac{1}{4}$、第二主軸の逆二乗屈折率 $b=\frac{1}{9}$、第三主軸の逆二乗屈折率 $c=\frac{1}{16}$、波面法線の第一成分 $l=0$、波面法線の第二成分 $m=0$、波面法線の第三成分 $n=1$ である。

三次元の跡から法線成分への寄与を引く。

$$
\begin{aligned}
T&=- a l^{2} + a - b m^{2} + b - c n^{2} + c,\\
T&= - \left(\frac{1}{4}\right) \cdot \left(0\right)^{2} + \left(\frac{1}{4}\right) - \left(\frac{1}{9}\right) \cdot \left(0\right)^{2} + \left(\frac{1}{9}\right) - \left(\frac{1}{16}\right) \cdot \left(1\right)^{2} + \left(\frac{1}{16}\right)=\frac{13}{36}.
\end{aligned}
$$

横基底の外積の各成分を二乗して重みを取る。

$$
\begin{aligned}
D&=a b n^{2} + a c m^{2} + b c l^{2},\\
D&= \left(\frac{1}{4}\right) \cdot \left(\frac{1}{9}\right) \cdot \left(1\right)^{2} + \left(\frac{1}{4}\right) \cdot \left(\frac{1}{16}\right) \cdot \left(0\right)^{2} + \left(\frac{1}{9}\right) \cdot \left(\frac{1}{16}\right) \cdot \left(0\right)^{2}=\frac{1}{36}.
\end{aligned}
$$

二次の特性多項式の判別式を求める。

$$
\begin{aligned}
g&=\sqrt{- 4 D + T^{2}},\\
g&= \sqrt{- 4 \cdot \left(\frac{1}{36}\right) + \left(\frac{13}{36}\right)^{2}}=\frac{5}{36}.
\end{aligned}
$$

固有値が大きい方で、逆平方根は小さくなる。

$$
\begin{aligned}
n_{\min}&=\sqrt{2} \sqrt{\frac{1}{T + g}},\\
n_{\min}&= \sqrt{2} \cdot \sqrt{\frac{1}{\left(\frac{13}{36}\right) + \left(\frac{5}{36}\right)}}=2.
\end{aligned}
$$

もう一方の固有偏光の屈折率を求める。

$$
\begin{aligned}
n_{\max}&=\sqrt{2} \sqrt{\frac{1}{T - g}},\\
n_{\max}&= \sqrt{2} \cdot \sqrt{\frac{1}{\left(\frac{13}{36}\right) - \left(\frac{5}{36}\right)}}=3.
\end{aligned}
$$

横の二主軸の屈折率がそのまま出る。
:::

:::example[全成分を持つ方向]
三主軸に等しい角をなす方向を選ぶ。

入力は第一主軸の逆二乗屈折率 $a=\frac{1}{4}$、第二主軸の逆二乗屈折率 $b=\frac{1}{9}$、第三主軸の逆二乗屈折率 $c=\frac{1}{16}$、波面法線の第一成分 $l\simeq 0.57735$、波面法線の第二成分 $m\simeq 0.57735$、波面法線の第三成分 $n\simeq 0.57735$ である。

三次元の跡から法線成分への寄与を引く。

$$
\begin{aligned}
T&=- a l^{2} + a - b m^{2} + b - c n^{2} + c,\\
T&\simeq - \left(\frac{1}{4}\right) \cdot \left(0.5773502692\right)^{2} + \left(\frac{1}{4}\right) - \left(\frac{1}{9}\right) \cdot \left(0.5773502692\right)^{2} + \left(\frac{1}{9}\right) - \left(\frac{1}{16}\right) \cdot \left(0.5773502692\right)^{2} + \left(\frac{1}{16}\right)=\frac{61}{216}.
\end{aligned}
$$

横基底の外積の各成分を二乗して重みを取る。

$$
\begin{aligned}
D&=a b n^{2} + a c m^{2} + b c l^{2},\\
D&\simeq \left(\frac{1}{4}\right) \cdot \left(\frac{1}{9}\right) \cdot \left(0.5773502692\right)^{2} + \left(\frac{1}{4}\right) \cdot \left(\frac{1}{16}\right) \cdot \left(0.5773502692\right)^{2} + \left(\frac{1}{9}\right) \cdot \left(\frac{1}{16}\right) \cdot \left(0.5773502692\right)^{2}=\frac{29}{1728}.
\end{aligned}
$$

二次の特性多項式の判別式を求める。

$$
\begin{aligned}
g&=\sqrt{- 4 D + T^{2}},\\
g&= \sqrt{- 4 \cdot \left(\frac{29}{1728}\right) + \left(\frac{61}{216}\right)^{2}}\simeq 0.11236.
\end{aligned}
$$

固有値が大きい方で、逆平方根は小さくなる。

$$
\begin{aligned}
n_{\min}&=\sqrt{2} \sqrt{\frac{1}{T + g}},\\
n_{\min}&\simeq \sqrt{2} \cdot \sqrt{\frac{1}{\left(\frac{61}{216}\right) + \left(0.1123579731\right)}}\simeq 2.2508.
\end{aligned}
$$

もう一方の固有偏光の屈折率を求める。

$$
\begin{aligned}
n_{\max}&=\sqrt{2} \sqrt{\frac{1}{T - g}},\\
n_{\max}&\simeq \sqrt{2} \cdot \sqrt{\frac{1}{\left(\frac{61}{216}\right) - \left(0.1123579731\right)}}\simeq 3.4295.
\end{aligned}
$$

どちらの屈折率も元の楕円体の最大軸と最小軸の間にある。
:::

:::example[円断面を検算]
光学軸になる方向を指定し、二根が一致することを確かめる。

入力は第一主軸の逆二乗屈折率 $a=\frac{1}{2}$、第二主軸の逆二乗屈折率 $b=\frac{3}{8}$、第三主軸の逆二乗屈折率 $c=\frac{1}{4}$、波面法線の第一成分 $l=1/\sqrt{2}$、波面法線の第二成分 $m=0$、波面法線の第三成分 $n=1/\sqrt{2}$ である。

三次元の跡から法線成分への寄与を引く。

$$
\begin{aligned}
T&=- a l^{2} + a - b m^{2} + b - c n^{2} + c,\\
T&\simeq - \left(\frac{1}{2}\right) \cdot \left(1/\sqrt{2}\right)^{2} + \left(\frac{1}{2}\right) - \left(\frac{3}{8}\right) \cdot \left(0\right)^{2} + \left(\frac{3}{8}\right) - \left(\frac{1}{4}\right) \cdot \left(1/\sqrt{2}\right)^{2} + \left(\frac{1}{4}\right)=\frac{3}{4}.
\end{aligned}
$$

横基底の外積の各成分を二乗して重みを取る。

$$
\begin{aligned}
D&=a b n^{2} + a c m^{2} + b c l^{2},\\
D&\simeq \left(\frac{1}{2}\right) \cdot \left(\frac{3}{8}\right) \cdot \left(1/\sqrt{2}\right)^{2} + \left(\frac{1}{2}\right) \cdot \left(\frac{1}{4}\right) \cdot \left(0\right)^{2} + \left(\frac{3}{8}\right) \cdot \left(\frac{1}{4}\right) \cdot \left(1/\sqrt{2}\right)^{2}=\frac{9}{64}.
\end{aligned}
$$

二次の特性多項式の判別式を求める。

$$
\begin{aligned}
g&=\sqrt{- 4 D + T^{2}},\\
g&= \sqrt{- 4 \cdot \left(\frac{9}{64}\right) + \left(\frac{3}{4}\right)^{2}}=0.
\end{aligned}
$$

固有値が大きい方で、逆平方根は小さくなる。

$$
\begin{aligned}
n_{\min}&=\sqrt{2} \sqrt{\frac{1}{T + g}},\\
n_{\min}&= \sqrt{2} \cdot \sqrt{\frac{1}{\left(\frac{3}{4}\right) + \left(0\right)}}\simeq 1.633.
\end{aligned}
$$

もう一方の固有偏光の屈折率を求める。

$$
\begin{aligned}
n_{\max}&=\sqrt{2} \sqrt{\frac{1}{T - g}},\\
n_{\max}&= \sqrt{2} \cdot \sqrt{\frac{1}{\left(\frac{3}{4}\right) - \left(0\right)}}\simeq 1.633.
\end{aligned}
$$

判別式がゼロになった。ここでは屈折率の値から偏光の固有方向を区別できない。
:::

## 一軸性結晶と複屈折

光学軸を $z$ とし、$\eps_r=\diag(n_o^2,n_o^2,n_e^2)$ とする。方向 $(\sin\theta,0,\cos\theta)$ に垂直な基底 $\hat{\vb y},(\cos\theta,0,-\sin\theta)$ に式 :eqref[eq-optics-indicatrix] を制限すると行列は対角であり、
$$
n_{\mathrm{ord}}=n_o,\qquad n_{\mathrm{ext}}^{-2}=\cos^2\theta/n_o^2+\sin^2\theta/n_e^2
$$
を得る。この二偏光で屈折率が異なる性質が :term[複屈折]{en="birefringence" yomi="ふくくっせつ"} である。異常光の分散面は $k_x^2/n_e^2+k_z^2/n_o^2=k_0^2$ である。屈折率自体の周波数依存を無視して微分すると、群速度は $(k_x/n_e^2,0,k_z/n_o^2)$ に平行になる。従って位相法線とエネルギーの進む方向にはずれが生じ得る。

:::example[光学軸上]
$\theta=0$ では $n_{\mathrm{ext}}=n_o$ で、二偏光の波数は同じである。
:::
:::example[遅相板の厚さ]
軸に垂直な入射で $\lambda_0=600\,\mathrm{nm}$、$n_e-n_o=0.01$ なら、四分の一波長板の最小厚は $\lambda_0/[4(n_e-n_o)]=15\,\mathrm{\mu m}$ である。
:::

### 位相速度とエネルギー流の角度
異常波の群速度の方向を $z$ から $\theta_g$ とすると、分散面を微分した結果から $\tan\theta_g=(n_o^2/n_e^2)\tan\theta$ である。二つの角の差がウォークオフ角であり、板の内部でビーム中心を横にずらす。位相整合を満たしても有限幅のビームの重なりはこのずれで減少し得る。

:::example[角度のずれ]
$n_o=1.5,n_e=1.6,\theta=45^\circ$ なら $\theta_g=\arctan[(1.5/1.6)^2]\simeq41.3^\circ$。異常波のエネルギーの方向は位相法線から約 $3.7^\circ$ 光学軸側にずれる。
:::

### 異常光のエネルギー方向を場から確認する
光学軸のある $xz$ 面に波数を取ると、異常光では $\vb H=H_y\hat{\vb y}$ である。実際、$\vb k=(k_x,0,k_z)$ に対して $\vb D=-\vb k\times\vb H/\omega=(k_z,0,-k_x)H_y/\omega$ が :eqref[eq-optics-indicatrix] の異常光の方向と一致する。$\vb E=\eps_r^{-1}\vb D/\eps_0$ を :eqref[eq-optics-irradiance] の Poynting ベクトルに代入すれば
$$
\langle\vb S\rangle=
\frac{|H_y|^2}{2\omega\eps_0}
\left(\frac{k_x}{n_e^2},0,\frac{k_z}{n_o^2}\right)
$$
を得る。これは屈折率を一定として求めた群速度の方向と一致する。界面では :eqref[eq-optics-snell] の位相整合と同様に、接線方向の波数を保存する。次に常光の球または異常光の楕円から法線波数を求め、エネルギー流が媒質内部へ向く根を取る。

:::example[入射面に光学軸がある結晶]
空気から斜入射し、結晶面の法線と光学軸がともに $z$ を向くとする。$n_o=1.5$、$n_e=2$、入射角 $30^\circ$ では $k_x/k_0=1/2$。常光と異常光の法線波数はそれぞれ $k_0\sqrt{2}$ と $k_0(3\sqrt{15}/8)$ である。異常光の流れの傾きは $\langle S_x\rangle/\langle S_z\rangle=(n_o/n_e)^2 k_x/k_z=3/(4\sqrt{15})$ となる。ここでは位相法線とエネルギー方向を別々に求めた。
:::

### 位相遅れと横ずれを同時に評価する

波面の法線が板の法線に一致したとする。位相の差は $2\pi(n_{\rm ext}-n_o)d/\lambda_0$ だが、強度包絡線の中心は群速度方向に移動する。異常光のずれは板の厚さ $d$ に対して $d\tan(\theta_g-\theta)$、常光ではゼロである。これを使って遅相の大きさとビームの重なりを同時に見積もる。$n_e>n_o$ なら $0<\theta<\pi/2$ で $\theta_g<\theta$、負一軸性なら不等号が逆になる。結晶と外空間での反射と回折を無視できる場合にこの幾何学的なずれを使う。

:::example[正一軸性の斜めカット]
板の法線を結晶軸から傾けて入射する。

入力は常光の主屈折率 $n_o=\frac{3}{2}$、異常光の主屈折率 $n_e=\frac{8}{5}$、光学軸から見た波面法線角 $\theta\simeq 0.7854\,\mathrm{rad}$、厚さ $d=\frac{1}{1000}\,\mathrm m$、真空波長 $\lambda_0\simeq 6\times 10^{-7}\,\mathrm m$ である。

指定方向で楕円を切った長さである。

$$
\begin{aligned}
n_{\rm ext}&=\frac{1}{\sqrt{\frac{\cos^{2}{\left(\theta \right)}}{n_o^{2}} + \frac{\sin^{2}{\left(\theta \right)}}{n_e^{2}}}},\\
n_{\rm ext}&\simeq \frac{1}{\sqrt{\frac{\cos^{2}{\left(\left(0.7853981634\right) \right)}}{\left(\frac{3}{2}\right)^{2}} + \frac{\sin^{2}{\left(\left(0.7853981634\right) \right)}}{\left(\frac{8}{5}\right)^{2}}}}\simeq 1.5476.
\end{aligned}
$$

軸に対する異常光の群速度方向を求める。

$$
\begin{aligned}
\theta_g&=\operatorname{atan}{\left(\frac{n_o^{2} \tan{\left(\theta \right)}}{n_e^{2}} \right)},\\
\theta_g/(\mathrm{rad})&\simeq \operatorname{atan}{\left(\frac{\left(\frac{3}{2}\right)^{2} \cdot \tan{\left(\left(0.7853981634\right) \right)}}{\left(\frac{8}{5}\right)^{2}} \right)}\simeq 0.72104.
\end{aligned}
$$

波面法線からの符号付き変位角を求める。

$$
\begin{aligned}
\rho&=\theta_g - \theta,\\
\rho/(\mathrm{rad})&\simeq \left(0.7210381109\right) - \left(0.7853981634\right)\simeq -0.06436.
\end{aligned}
$$

板を厚さ方向に進む間の横変位へ換算する。

$$
\begin{aligned}
\Delta x&=d \tan{\left(\rho \right)},\\
\Delta x/(\mathrm m)&\simeq \left(\frac{1}{1000}\right) \cdot \tan{\left(\left(-0.06436005254\right) \right)}\simeq -6.4449\times 10^{-5}.
\end{aligned}
$$

二偏光の位相差は法線波数の差で評価する。

$$
\begin{aligned}
\Delta\phi/(2\pi)&=\frac{d \left(n_{\rm ext} - n_o\right)}{\lambda_0},\\
\Delta\phi/(2\pi)&\simeq \frac{\left(\frac{1}{1000}\right) \cdot \left(\left(1.547582113\right) - \left(\frac{3}{2}\right)\right)}{\left(6\times 10^{-7}\right)}\simeq 79.304.
\end{aligned}
$$

エネルギーは法線から軸寄りへ外れ、位相は常光より進む。
:::

:::example[負一軸性へ変える]
厚さと方向を固定して異常光の主屈折率を下げる。

入力は常光の主屈折率 $n_o=\frac{3}{2}$、異常光の主屈折率 $n_e=\frac{7}{5}$、光学軸から見た波面法線角 $\theta\simeq 0.7854\,\mathrm{rad}$、厚さ $d=\frac{1}{1000}\,\mathrm m$、真空波長 $\lambda_0\simeq 6\times 10^{-7}\,\mathrm m$ である。

指定方向で楕円を切った長さである。

$$
\begin{aligned}
n_{\rm ext}&=\frac{1}{\sqrt{\frac{\cos^{2}{\left(\theta \right)}}{n_o^{2}} + \frac{\sin^{2}{\left(\theta \right)}}{n_e^{2}}}},\\
n_{\rm ext}&\simeq \frac{1}{\sqrt{\frac{\cos^{2}{\left(\left(0.7853981634\right) \right)}}{\left(\frac{3}{2}\right)^{2}} + \frac{\sin^{2}{\left(\left(0.7853981634\right) \right)}}{\left(\frac{7}{5}\right)^{2}}}}\simeq 1.4474.
\end{aligned}
$$

軸に対する異常光の群速度方向を求める。

$$
\begin{aligned}
\theta_g&=\operatorname{atan}{\left(\frac{n_o^{2} \tan{\left(\theta \right)}}{n_e^{2}} \right)},\\
\theta_g/(\mathrm{rad})&\simeq \operatorname{atan}{\left(\frac{\left(\frac{3}{2}\right)^{2} \cdot \tan{\left(\left(0.7853981634\right) \right)}}{\left(\frac{7}{5}\right)^{2}} \right)}\simeq 0.85417.
\end{aligned}
$$

波面法線からの符号付き変位角を求める。

$$
\begin{aligned}
\rho&=\theta_g - \theta,\\
\rho/(\mathrm{rad})&\simeq \left(0.8541731332\right) - \left(0.7853981634\right)\simeq 0.068775.
\end{aligned}
$$

板を厚さ方向に進む間の横変位へ換算する。

$$
\begin{aligned}
\Delta x&=d \tan{\left(\rho \right)},\\
\Delta x/(\mathrm m)&\simeq \left(\frac{1}{1000}\right) \cdot \tan{\left(\left(0.06877496978\right) \right)}\simeq 6.8884\times 10^{-5}.
\end{aligned}
$$

二偏光の位相差は法線波数の差で評価する。

$$
\begin{aligned}
\Delta\phi/(2\pi)&=\frac{d \left(n_{\rm ext} - n_o\right)}{\lambda_0},\\
\Delta\phi/(2\pi)&\simeq \frac{\left(\frac{1}{1000}\right) \cdot \left(\left(1.447415585\right) - \left(\frac{3}{2}\right)\right)}{\left(6\times 10^{-7}\right)}\simeq -87.641.
\end{aligned}
$$

ずれと位相差の符号がともに反対になる。
:::

:::example[光学軸に垂直に近い板]
法線が軸にほぼ垂直な向きを調べる。

入力は常光の主屈折率 $n_o=\frac{3}{2}$、異常光の主屈折率 $n_e=\frac{8}{5}$、光学軸から見た波面法線角 $\theta\simeq 1.5533\,\mathrm{rad}$、厚さ $d=\frac{1}{1000}\,\mathrm m$、真空波長 $\lambda_0\simeq 6\times 10^{-7}\,\mathrm m$ である。

指定方向で楕円を切った長さである。

$$
\begin{aligned}
n_{\rm ext}&=\frac{1}{\sqrt{\frac{\cos^{2}{\left(\theta \right)}}{n_o^{2}} + \frac{\sin^{2}{\left(\theta \right)}}{n_e^{2}}}},\\
n_{\rm ext}&\simeq \frac{1}{\sqrt{\frac{\cos^{2}{\left(\left(1.553343034\right) \right)}}{\left(\frac{3}{2}\right)^{2}} + \frac{\sin^{2}{\left(\left(1.553343034\right) \right)}}{\left(\frac{8}{5}\right)^{2}}}}\simeq 1.6.
\end{aligned}
$$

軸に対する異常光の群速度方向を求める。

$$
\begin{aligned}
\theta_g&=\operatorname{atan}{\left(\frac{n_o^{2} \tan{\left(\theta \right)}}{n_e^{2}} \right)},\\
\theta_g/(\mathrm{rad})&\simeq \operatorname{atan}{\left(\frac{\left(\frac{3}{2}\right)^{2} \cdot \tan{\left(\left(1.553343034\right) \right)}}{\left(\frac{8}{5}\right)^{2}} \right)}\simeq 1.5509.
\end{aligned}
$$

波面法線からの符号付き変位角を求める。

$$
\begin{aligned}
\rho&=\theta_g - \theta,\\
\rho/(\mathrm{rad})&\simeq \left(1.550938952\right) - \left(1.553343034\right)\simeq -0.0024041.
\end{aligned}
$$

板を厚さ方向に進む間の横変位へ換算する。

$$
\begin{aligned}
\Delta x&=d \tan{\left(\rho \right)},\\
\Delta x/(\mathrm m)&\simeq \left(\frac{1}{1000}\right) \cdot \tan{\left(\left(-0.002404082031\right) \right)}\simeq -2.4041\times 10^{-6}.
\end{aligned}
$$

二偏光の位相差は法線波数の差で評価する。

$$
\begin{aligned}
\Delta\phi/(2\pi)&=\frac{d \left(n_{\rm ext} - n_o\right)}{\lambda_0},\\
\Delta\phi/(2\pi)&\simeq \frac{\left(\frac{1}{1000}\right) \cdot \left(\left(1.599966429\right) - \left(\frac{3}{2}\right)\right)}{\left(6\times 10^{-7}\right)}\simeq 166.61.
\end{aligned}
$$

位相差は大きい一方、横ずれは小さい。遅相と分離の大きさは比例しない。
:::

### 二つの偏光を角度で分ける Wollaston プリズム
同じ一軸結晶を二つのくさびにし、外面を平行にして接合する。法線入射の波数を $+z$、接合面の法線の傾きを $\alpha$ とし、くさびの面内を $xz$ 平面とする。前側の光学軸を $x$、後側を $y$ に取ると、垂直伝播では前側で $x$ 偏光が異常光、$y$ 偏光が常光となり、後側ではその役割が交換する。誘電率テンソルは両側で $xz$ 面に関する鏡映対称性を持つ。接合面の法線も $xz$ 面にあるため、電場がこの面内にある成分と $y$ 方向にある成分は境界で混ざらない。

くさび角が小さいとき、異常光の屈折率は法線に垂直な光学軸での $n_e$ から二次以上しか変わらない。従って一次では、接合面で接線波数を保存して $n_{\rm before}(-\alpha)=n_{\rm after}(\theta-\alpha)$ と書ける。出口の平面でも同様に接線波数を保存する。外側を空気とすると、出射角は $x$ 偏光で $(n_o-n_e)\alpha$、$y$ 偏光で $(n_e-n_o)\alpha$ となる。これが :term[Wollaston プリズム]{en="Wollaston prism" yomi="うぉらすとんぷりずむ"} の偏光分離である。$n_o-n_e=0.10$、$\alpha=5^\circ$ なら、近軸の二出射角のなす角は約 $1.0^\circ$ となる。

## 旋光性と Faraday 効果

横誘電率が $\begin{pmatrix}a&-ig\\ig&a\end{pmatrix}$、$a>|g|$ のとき、二円偏光の固有値は $a\pm g$ である。各成分が位相 $\phi_\pm$ を持つようになった後、$x$ 直線の入射を再構成する。等振幅の和の $x$ 成分は $(e^{i\phi_+}+e^{i\phi_-})/2$、$y$ 成分は $i(e^{i\phi_+}-e^{i\phi_-})/2$ なので、出射は共通位相を除いて $(\cos\Delta,-\sin\Delta)$、$\Delta=(\phi_+-\phi_-)/2$ となる。

相反的な :term[旋光性]{en="optical activity" yomi="せんこうせい"} の最低次模型では $g$ が波数に比例し、逆進すると二円成分の屈折率が交換する。これに対して :term[Faraday 効果]{en="Faraday effect" yomi="ふぁらでーこうか"} では $g$ が固定磁場に比例する。従って同じ実験室座標で、逆進しても固有値の割当ては保たれる。小さい外場で回転角を $VBL$ と書いたとき、$V$ は $\mathrm{rad/(T\,m)}$ 単位の Verdet 定数である。

:::example[往復の測定]
鏡の Jones 行列が固定横座標で定数倍の単位行列なら、片道 $20^\circ$ の Faraday 回転は往復で $40^\circ$ となる。相反的旋光では回転行列が互いに逆で、往復角はゼロになる。
:::
:::example[磁場の向き]
$V=10\,\mathrm{rad/(T\,m)}$、$B=-0.5\,\mathrm T,L=0.02\,\mathrm m$ なら角度は $-0.10\,\mathrm{rad}$ である。
:::

:::example[偏光選別と非相反性]
入口の理想偏光子を $x$、出口を $45^\circ$ にして、間に順行で $+45^\circ$ の Faraday 素子を入れる。順行の偏光は出口の軸に一致する。逆から出口の偏光子を通った $45^\circ$ の状態には、実験室座標でさらに $+45^\circ$ が加わる。入口では $y$ となるので遮断される。実際のアイソレータの隔離度には偏光子の漏れも入る。
:::

### 磁場で円偏光の共鳴が分かれる
軸方向の静磁場 $B\hat{\vb z}$ と、束縛電荷の Lorentz 振動子を考える。変位を $\vb X e^{-i\omega t}$ と書くと磁気力は $q(-i\omega\vb X)\times B\hat{\vb z}=i\omega qB\hat{\vb z}\times\vb X$。従って $D(\omega)=\omega_0^2-\omega^2-i\gamma\omega$、$\Omega_B=qB/m$ と置けば
$$
[D(\omega)I-i\omega\Omega_B(\hat{\vb z}\times)]\vb X=(q/m)\vb E .
$$
二円偏光で $\hat{\vb z}\times$ の固有値は $\mp i$ である。数密度 $N$ なら感受率は $\chi_\pm=F/(D\mp\omega\Omega_B)$、$F=Nq^2/(\eps_0m)$。磁場の一次では差が $\chi_+-\chi_-\simeq 2F\omega\Omega_B/D^2$ となる。透明域では屈折率差が偏光を回転させ、吸収域では虚部の差が二円成分の振幅も変えて出力を楕円にする。この模型で磁場を反転すると共鳴の割当てが交換するが、実験室座標で進行方向だけを反転しても係数に $k_z$ がないので交換しない。

### Faraday 回転子と偏光子で往復を比較する

偏光子を角 $0,\pi/4$ に置き、その間の回転子が順方向で $\theta$ 回すとする。入力を最初の偏光子に平行とすれば、順方向の出力は $\cos^2(\theta-\pi/4)$。逆方向の光が第二偏光子に平行なら、固定座標で同じ方向へ回す Faraday 回転子では第一偏光子の透過率が $\cos^2(\theta+\pi/4)$。これは二つの単位直線ベクトルの内積の二乗から求まる。相反な旋光子では逆進の回転が $-\theta$ となるため、逆方向の透過率は順方向と一致する。有限の偏光子の消光や吸収は別の損失として扱う。

:::example[理想的な四十五度回転]
磁場と長さを逆方向が遮られる値に合わせる。

入力は符号込みの Verdet 定数 $V=100\,\mathrm{rad/(T\,m)}$、軸に沿った磁束密度 $B\simeq 0.7854\,\mathrm T$、回転子の長さ $L=\frac{1}{100}\,\mathrm m$ である。

弱磁場で回転が磁場と光路に比例する近似を使う。

$$
\begin{aligned}
\theta&=B L V,\\
\theta/(\mathrm{rad})&\simeq \left(0.7853981634\right) \cdot \left(\frac{1}{100}\right) \cdot \left(100\right)\simeq 0.7854.
\end{aligned}
$$

出力軸に対するずれを求める。

$$
\begin{aligned}
\delta_{\to}&=\theta - \frac{\pi}{4},\\
\delta_{\to}/(\mathrm{rad})&\simeq \left(0.7853981634\right) - \frac{\pi}{4}=0.
\end{aligned}
$$

二つの偏光の軸への射影を取る。

$$
\begin{aligned}
T_{\to}&=\cos^{2}{\left(\delta_{\to} \right)},\\
T_{\to}&= \cos^{2}{\left(\left(0\right) \right)}=1.
\end{aligned}
$$

逆進でも固定座標で同じ向きに回る。

$$
\begin{aligned}
T_{\gets}&=\cos^{2}{\left(\theta + \frac{\pi}{4} \right)},\\
T_{\gets}&\simeq \cos^{2}{\left(\left(0.7853981634\right) + \frac{\pi}{4} \right)}=0.
\end{aligned}
$$

理想偏光子で逆進波を遮る割合を求める。

$$
\begin{aligned}
A_{\gets}&=1 - T_{\gets},\\
A_{\gets}&= 1 - \left(0\right)=1.
\end{aligned}
$$

順方向は全透過し、逆方向は遮られる。
:::

:::example[磁場の大きさを下げる]
極性を保ち、磁場だけを減らす。

入力は符号込みの Verdet 定数 $V=100\,\mathrm{rad/(T\,m)}$、軸に沿った磁束密度 $B=\frac{7}{10}\,\mathrm T$、回転子の長さ $L=\frac{1}{100}\,\mathrm m$ である。

弱磁場で回転が磁場と光路に比例する近似を使う。

$$
\begin{aligned}
\theta&=B L V,\\
\theta/(\mathrm{rad})&= \left(\frac{7}{10}\right) \cdot \left(\frac{1}{100}\right) \cdot \left(100\right)=\frac{7}{10}.
\end{aligned}
$$

出力軸に対するずれを求める。

$$
\begin{aligned}
\delta_{\to}&=\theta - \frac{\pi}{4},\\
\delta_{\to}/(\mathrm{rad})&= \left(\frac{7}{10}\right) - \frac{\pi}{4}\simeq -0.085398.
\end{aligned}
$$

二つの偏光の軸への射影を取る。

$$
\begin{aligned}
T_{\to}&=\cos^{2}{\left(\delta_{\to} \right)},\\
T_{\to}&\simeq \cos^{2}{\left(\left(-0.0853981634\right) \right)}\simeq 0.99272.
\end{aligned}
$$

逆進でも固定座標で同じ向きに回る。

$$
\begin{aligned}
T_{\gets}&=\cos^{2}{\left(\theta + \frac{\pi}{4} \right)},\\
T_{\gets}&= \cos^{2}{\left(\left(\frac{7}{10}\right) + \frac{\pi}{4} \right)}\simeq 0.0072751.
\end{aligned}
$$

理想偏光子で逆進波を遮る割合を求める。

$$
\begin{aligned}
A_{\gets}&=1 - T_{\gets},\\
A_{\gets}&\simeq 1 - \left(0.007275135006\right)\simeq 0.99272.
\end{aligned}
$$

軸からの誤差により順方向の損失と逆方向の漏れが現れる。
:::

:::example[磁場を反転する]
四十五度回転の設定から磁場の符号を反転する。

入力は符号込みの Verdet 定数 $V=100\,\mathrm{rad/(T\,m)}$、軸に沿った磁束密度 $B\simeq -0.7854\,\mathrm T$、回転子の長さ $L=\frac{1}{100}\,\mathrm m$ である。

弱磁場で回転が磁場と光路に比例する近似を使う。

$$
\begin{aligned}
\theta&=B L V,\\
\theta/(\mathrm{rad})&\simeq \left(-0.7853981634\right) \cdot \left(\frac{1}{100}\right) \cdot \left(100\right)\simeq -0.7854.
\end{aligned}
$$

出力軸に対するずれを求める。

$$
\begin{aligned}
\delta_{\to}&=\theta - \frac{\pi}{4},\\
\delta_{\to}/(\mathrm{rad})&\simeq \left(-0.7853981634\right) - \frac{\pi}{4}\simeq -1.5708.
\end{aligned}
$$

二つの偏光の軸への射影を取る。

$$
\begin{aligned}
T_{\to}&=\cos^{2}{\left(\delta_{\to} \right)},\\
T_{\to}&\simeq \cos^{2}{\left(\left(-1.570796327\right) \right)}=0.
\end{aligned}
$$

逆進でも固定座標で同じ向きに回る。

$$
\begin{aligned}
T_{\gets}&=\cos^{2}{\left(\theta + \frac{\pi}{4} \right)},\\
T_{\gets}&\simeq \cos^{2}{\left(\left(-0.7853981634\right) + \frac{\pi}{4} \right)}=1.
\end{aligned}
$$

理想偏光子で逆進波を遮る割合を求める。

$$
\begin{aligned}
A_{\gets}&=1 - T_{\gets},\\
A_{\gets}&= 1 - \left(1\right)=0.
\end{aligned}
$$

透過しやすい向きと遮られる向きが入れ替わる。
:::

## 演習問題

::::exercise{level=basic}
$(1,-i)^\mathsf T/\sqrt2$ の $S_3$ を求めよ。

:::solution
$C_{yx}=-i/2$ なので $S_3=2\Im C_{yx}=-1$ である。
:::
::::

::::exercise{level=basic}
軸から $60^\circ$ ずれた直線の偏光子透過率を求めよ。

:::solution
射影は $\cos60^\circ=1/2$ であり、強度は $1/4$ 倍になる。
:::
::::

::::exercise{level=basic}
同じ四分の一波長板を二枚、軸を合わせて置く。行列を求めよ。

:::solution
$\diag(1,i)$ を二乗すると $\diag(1,-1)$ である。
:::
::::

::::exercise{level=basic}
$x,y$ への非干渉強度が $3,1$ のとき偏光度を求めよ。

:::solution
$C=\diag(3,1)$ なので、差を和で割って $(3-1)/4=1/2$ である。
:::
::::

::::exercise{level=advanced}
$C=\begin{pmatrix}2&1\\1&2\end{pmatrix}$ の完全偏光成分を求めよ。

:::solution
固有値は $3,1$。小さい方の単位行列倍を引くと全成分が 1 になる。これは $45^\circ$ の偏光を強度 2 だけ持つ相関行列である。
:::
::::

::::exercise{level=advanced}
一軸結晶を $\hat{\vb x}$ に進む二固有波を求めよ。

:::solution
横平面は $yz$ であり、式 :eqref[eq-optics-indicatrix] の行列は $\diag(n_o^{-2},n_e^{-2})$。従って $y$ が常光、$z$ が異常光である。
:::
::::

::::exercise{level=advanced}
$V=10\,\mathrm{rad/(T\,m)}, B=0.5\,\mathrm T,L=0.02\,\mathrm m$ の回転角を求めよ。

:::solution
積を取って $10(0.5)(0.02)=0.10\,\mathrm{rad}$ である。
:::
::::

::::exercise{level=advanced}
両対角が正の対角 $C$ が Jones ベクトル一つでは表せないことを示せ。

:::solution
Jones ベクトル一つの外積は階数 1 で行列式ゼロである。与えた行列の行列式は両対角の積で正なので矛盾する。
:::
::::

::::exercise{level=advanced}
純粋状態と検光子の規格化 Stokes の三成分がそれぞれ $(0,3/5,4/5)$、$(0,0,-1)$ である。透過率を求めよ。次に検光子を対蹠点の状態に換え、二透過率の和が一になることを示せ。

:::solution
二つの状態の内積は $-4/5$。Jones 射影の強度は Stokes の内積を $d$ として $(1+d)/2$ なので、透過率は $1/10$ である。検光子を対蹠点へ移すと内積が $+4/5$ になり、透過率は $9/10$。一般に $d$ と $-d$ の二透過率は $(1+d)/2+(1-d)/2=1$ を満たす。
:::
::::

## 確認テスト

::::quiz
本章の正円偏光について、規格化した Stokes の三成分 $(S_1,S_2,S_3)/S_0$ を選べ。

- [x] $(0,0,1)$
- [ ] $(1,0,0)$
- [ ] $(0,1,0)$

:::explanation
$J=(1,i)^\mathsf T/\sqrt2$ を入れると $C_{yx}=i/2$ で両対角が $1/2$ だから、規格化後は $(0,0,1)$ になる。
:::
::::

::::quiz{answer="0.5" tol="1e-6"}
直線偏光を $45^\circ$ ずれた偏光子に入れる。強度比を答えよ。

:::explanation
振幅が $\cos45^\circ$ 倍となるので強度はその二乗で $1/2$ 倍である。
:::
::::

::::quiz
一軸性結晶で光学軸に沿う固有屈折率はどうなるか。

- [x] 二偏光で一致する。
- [ ] 差が最大となる。
- [ ] 一方は伝播しない。

:::explanation
横断面には $n_o$ の二主軸が入り、断面が円となるため二偏光の屈折率がともに $n_o$ になる。
:::
::::

::::quiz{answer="30" tol="1e-6"}
片道 $15^\circ$ の Faraday 素子の往復角を度で答えよ。

:::explanation
二円偏光の相対位相差が同じ符号で二度加わるので、偏光方向は実験室座標で $2(15^\circ)=30^\circ$ 回る。
:::
::::

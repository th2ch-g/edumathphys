---
title: 調和振動子と応答
description: 減衰振動と強制振動から因果的な応答と散逸を導き、係数変調および弱い非線形性の効果を調べる。
notation:
  - symbol: '\omega_0,\omega_d'
    meaning: 無減衰の固有角振動数と不足減衰時の角振動数
  - symbol: '\gamma'
    meaning: 振幅の減衰率
  - symbol: 'G,\chi'
    meaning: 遅延インパルス応答と複素応答率
  - symbol: 'Q'
    meaning: 品質係数
---

## この章の目標

- 初期変位と初速度を一つの複素振幅で表し、実際の運動を復元する。
- 応答率の極、共鳴の線幅、吸収電力を同じ運動方程式から導く。
- 遅延応答の畳み込みから、切り替えを含む外力に対する運動を計算する。
- 外力による共鳴、係数変調による不安定化、非線形振動の振幅選択を区別する。

## 単振動の複素表示と位相

### 振動の時刻原点と初期条件

質量 $m>0$ の質点が、平衡点からの変位 $x$ に比例する復元力 $-kx$ を受けるとする。:ref[def-mech-harmonic]{full} に従い、$k>0$ として $\omega_0=\sqrt{k/m}$ と置く。変位は $\mathrm{m}$、質量は $\mathrm{kg}$、ばね定数は $\mathrm{N\,m^{-1}}$、角振動数は $\mathrm{s^{-1}}$ で数える。角度はラジアンである。変位の符号と速度の符号を同時に追うため、以下では時間因子を $e^{-i\omega t}$ に統一する。

:::definition[複素振幅と初期位相]{#def-wave-phasor}
実数値の :term[単振動]{en="simple harmonic motion" yomi="たんしんどう"} を $x(t)=\Re(Xe^{-i\omega_0t})$ と表したとき、$X$ を :term[複素振幅]{en="complex amplitude" yomi="ふくそしんぷく"} と呼ぶ。$X\ne0$ なら、$X=Ae^{i\phi}$ と書いて $A=|X|$ を振幅、$\phi$ を初期位相と呼ぶ。位相は $2\pi$ の整数倍を法として定まる。
:::

初期条件を $x(0)=x_0$、$\dot x(0)=v_0$ とすると、$\Re X=x_0$、$\Re(-i\omega_0X)=v_0$ である。したがって

$$
X=x_0+i\frac{v_0}{\omega_0},\qquad
x(t)=x_0\cos\omega_0t+\frac{v_0}{\omega_0}\sin\omega_0t.
$$

この式を二回微分すると $\ddot x+\omega_0^2x=0$ に戻る。一意性は二つの解の差 $y$ について $(\dot y^2+\omega_0^2y^2)/2$ を微分すれば分かる。初期値がともに零の差は、この非負量が全時刻で零であるため恒等的に零である。

:::example[位置から出発する振動と速度から出発する振動]
$\omega_0=4\,\mathrm{s^{-1}}$、$x_0=3\,\mathrm{mm}$、$v_0=16\,\mathrm{mm\,s^{-1}}$ なら $X=(3+4i)\,\mathrm{mm}$ である。振幅は $5\,\mathrm{mm}$、初期位相は $\arctan(4/3)$ である。最初の正の最大変位は $t=\arctan(4/3)/4\simeq0.232\,\mathrm{s}$ で生じる。
:::

:::example[測定の時刻原点をずらす]
$x(t)=2\,\mathrm{mm}\cos(10t)$ を $t=\pi/20\,\mathrm{s}$ から数え直す。新しい時刻を $\tau=t-\pi/20$ と置くと $X'=2e^{-i\pi/2}\,\mathrm{mm}=-2i\,\mathrm{mm}$ である。新原点で $x=0$、$\dot x=-20\,\mathrm{mm\,s^{-1}}$ となる。振幅は変わらず、位相だけが回転する。
:::

### 位相平面とエネルギー

複素量 $z=x+i\dot x/\omega_0$ は $\dot z=-i\omega_0 z$ を満たす。したがって $(x,\dot x/\omega_0)$ 平面では状態が円周上を時計回りに進む。力学的エネルギーは $E=m\omega_0^2|z|^2/2$ であり、円の半径が一定であることとエネルギー保存は同じ内容である。速度を変位と同じ単位に直すための $\omega_0$ を省くと、位相軌道は円でなく楕円になる。

測定器が実変位を表示する場合、$|Xe^{-i\omega t}|$ を取ると一定値になり、振動そのものを消してしまう。変位の復元には実部、最大変位の計算には絶対値を使う。二つの同周波数成分の合成では実部を取る前に $X_1+X_2$ を作る。例えば $3\,\mathrm{mm}$ と $4i\,\mathrm{mm}$ の合成振幅は $7\,\mathrm{mm}$ でなく $5\,\mathrm{mm}$ である。異なる周波数の成分の和には一般に単一の位相を割り当てられない。

### 直交する振動を重ねる

変位を二つの独立な方向で測るときは、複素振幅を足す代わりに、それぞれの成分から軌道を作る。例えば $x=A\cos\omega t$、$y=B\cos(\omega t-\delta)$ とし、$A,B>0$ とする。第 2 式を加法定理で展開すると $y/B-(x/A)\cos\delta=\sin\omega t\sin\delta$ である。この式を二乗して $\sin^2\omega t=1-x^2/A^2$ を入れると

$$
\frac{x^2}{A^2}+\frac{y^2}{B^2}
-2\frac{xy}{AB}\cos\delta=\sin^2\delta
$$

となる。$\sin\delta\ne0$ なら円周の可逆な線形変換による楕円で、$\delta=0,\pi$ なら線分である。例えば $A=B$、$\delta=\pi/2$ では、$t=0$ に正の $x$ 軸上から正の $y$ 方向へ動く円となる。$\delta=-\pi/2$ に変えると同じ円を逆向きに回るので、軌道の形だけでは回転方向は決まらない。平行と直交の重ね合わせについては :cite[french-vibrations] も参照できる。

### 二つの測定から初期値を復元する

初速度を直接測れないときは、既知の $\omega_0$ に対し二時刻の変位を使う。$t=0$ と $t=\tau$ での値を $x_0,x_\tau$ とすれば、$\sin\omega_0\tau\ne0$ のとき

$$
v_0=\omega_0\frac{x_\tau-x_0\cos\omega_0\tau}{\sin\omega_0\tau}
$$

である。分子の第 2 項は、初期変位だけから予測される $\tau$ での寄与である。測定誤差をそれぞれ $\delta_0,\delta_\tau$ とすると、復元速度の誤差は $\omega_0(\delta_\tau-\delta_0\cos\omega_0\tau)/\sin\omega_0\tau$ である。半周期の整数倍では二つの変位が同じ成分しか含まないため、速度の情報を失う。四分の一周期では $\cos\omega_0\tau=0$ で、二つの測定値がそのまま複素振幅の実部と虚部になる。

:::example[二時刻の読みを複素平面に置く]
$\omega_0=8\,\mathrm{s^{-1}}$、$x(0)=1\,\mathrm{mm}$、$x(\pi/16)= -2\,\mathrm{mm}$ とする。$\omega_0\tau=\pi/2$ なので $X=(1-2i)\,\mathrm{mm}$、$v_0=-16\,\mathrm{mm\,s^{-1}}$ である。例えば $m=0.25\,\mathrm{kg}$ ならエネルギーは $m\omega_0^2|X|^2/2=40\,\mu\mathrm J$ である。$t=\pi/8$ の測定値は $-1\,\mathrm{mm}$ となるが、これを $t=\pi/16$ の代わりに使っても虚部は決まらない。
:::

### 位相平面で力積の作用を追う

摩擦のない振動子について、状態を $(x,v/\omega_0)$ で表す。各成分は長さの単位を持つ。時間が $\tau$ だけ進むと $X=x+iv/\omega_0$ は $Xe^{-i\omega_0\tau}$ へ移る。実部と虚部を取れば、この写像は

$$
\begin{pmatrix}x(t+\tau)\\v(t+\tau)/\omega_0\end{pmatrix}
=
\begin{pmatrix}
\cos\omega_0\tau&\sin\omega_0\tau\\
-\sin\omega_0\tau&\cos\omega_0\tau
\end{pmatrix}
\begin{pmatrix}x(t)\\v(t)/\omega_0\end{pmatrix}
$$

となる。この行列の列は長さ $1$ で互いに直交する。したがって、軌道の半径と $E=m\omega_0^2\abs{X}^2/2$ が保存される。力積 $J$ を瞬間的に与えた場合は、運動方程式をその時刻の前後で積分すると $v$ が $J/m$ 跳ぶ。有限な変位によるばねの力積は積分区間を縮めると消え、変位は連続である。よってこの操作は複素平面で $X\mapsto X+iJ/(m\omega_0)$ となる。

:::example[速度を打ち消す向き]
$m=1\,\mathrm{kg}$、$\omega_0=2\,\mathrm{s^{-1}}$、$X=(3+4i)\,\mathrm{mm}$ の瞬間に $J=-8\,\mathrm{mN\,s}$ を与えると、虚部は $4-8/2=0$ となる。振幅は $5\,\mathrm{mm}$ から $3\,\mathrm{mm}$ に減り、エネルギーの比は $9/25$ である。速度と逆向きの力積を使うと振動を弱められる。
:::

:::example[二つの力積の位相差]
静止から同じ正の力積 $J$ を二度加える。間隔を $T/4$（$T=2\pi/\omega_0$）に取ると、二度目の直後は $X=iJ(-i+1)/(m\omega_0)$ である。振幅は一度だけの場合の $\sqrt2$ 倍である。これに対し間隔 $T/2$ では最初の虚部が反転し、二度目の力積で静止に戻る。
:::

## 減衰振動と強制振動の複素振幅

### 自由振動の三つの場合

速度に比例する抵抗 $-b\dot x$ を仮定し、$\gamma=b/(2m)\ge0$ とする。このモデルは粘性抵抗が線形の範囲に適用する。乾性摩擦や、速さの二乗に比例する抗力には別の方程式が必要である。外力を $F(t)$ として運動方程式は

$$
m(\ddot x+2\gamma\dot x+\omega_0^2x)=F(t).
$$

自由振動では $x=e^{-\gamma t}y$ を代入すると $y''+(\omega_0^2-\gamma^2)y=0$ になる。$s_0=v_0+\gamma x_0$ と置けば、初期条件 $y(0)=x_0$、$y'(0)=s_0$ から次の三式が出る。

| 条件 | $y(t)$ |
| --- | --- |
| $\gamma<\omega_0$（不足減衰） | $x_0\cos\omega_dt+(s_0/\omega_d)\sin\omega_dt$、$\omega_d=\sqrt{\omega_0^2-\gamma^2}$ |
| $\gamma=\omega_0$（臨界減衰） | $x_0+s_0t$ |
| $\gamma>\omega_0$（過減衰） | $x_0\cosh dt+(s_0/d)\sinh dt$、$d=\sqrt{\gamma^2-\omega_0^2}$ |

いずれの行も初期値および変換後の方程式を満たす。臨界値には両側から同じ極限がある。過減衰の遅い指数は $-\gamma+d=-\omega_0^2/(\gamma+d)$ である。$\gamma$ を極端に大きくするとその絶対値が小さくなるので、強い抵抗で静止位置への戻りが遅くなることがある。二つの実指数の係数の一方が初期条件によって零になる場合は、残った指数だけで減衰する。

### 力に対する変位の比

$F(t)=\Re(F_0e^{-i\omega t})$ を代入し、線形演算をしてから実部を取る。微分は $-i\omega$ の乗算に置き換わるので、定常解は

$$
X(\omega)=\chi(\omega)F_0,\qquad
\chi(\omega)=\frac{1}{m(\omega_0^2-\omega^2-2i\gamma\omega)}
\label{eq-wave-osc-susceptibility}
$$

となる。$\chi$ を変位の :term[応答率]{en="susceptibility" yomi="おうとうりつ"} と呼ぶ。単位は $\mathrm{m\,N^{-1}}$ である。$\gamma>0$ なら時間が経つと同次解が減衰するので、定常解は初期条件によらない。$\gamma>0$、$F_0>0$、$\omega>0$ とすると、振幅と遅れは

$$
A=\frac{F_0/m}{\sqrt{(\omega_0^2-\omega^2)^2+4\gamma^2\omega^2}},\qquad
\delta=\arg(\omega_0^2-\omega^2+2i\gamma\omega),\quad 0<\delta<\pi.
$$

この符号規約では $x=A\cos(\omega t-\delta)$ なので、変位の山は力の山から $\delta/\omega$ 遅れる。分母の実部が負になる高周波側では $\tan\delta$ だけで象限を判定できない。高周波では $X\simeq-F_0/(m\omega^2)$、低周波では $X\simeq F_0/k$ である。

:::example[臨界減衰の時刻歴]
$\gamma=\omega_0=5\,\mathrm{s^{-1}}$、$x_0=2\,\mathrm{mm}$、$v_0=0$ なら $x=(2\,\mathrm{mm})(1+5t)e^{-5t}$ である。$t=0.2\,\mathrm{s}$ では $4/e\,\mathrm{mm}$ となる。指数因子だけを残すと、初速度が零という条件を失う。
:::

:::example[変位と速度の位相]
$m=2\,\mathrm{kg}$、$\omega_0=3\,\mathrm{s^{-1}}$、$\gamma=1\,\mathrm{s^{-1}}$ とし、角振動数 $3\,\mathrm{s^{-1}}$ の力の振幅を $12\,\mathrm{N}$ とする。$X=i\,\mathrm m$ だから変位は力から $\pi/2$ 遅れる。一方、速度振幅 $-i\omega X$ は $3\,\mathrm{m\,s^{-1}}$ の正実数で、力と同位相である。この数値は線形モデルの例であり、実際のばねに $1\,\mathrm{m}$ の伸びを許すとは限らない。
:::

### 無減衰の極限と支点加振

$\gamma=0$、$\omega=\omega_0$ では代数式に有限値の解がない。$t\sin\omega_0t$ に $\partial_t^2+\omega_0^2$ を作用させると $2\omega_0\cos\omega_0t$ になるので、静止状態からの解は $x=F_0t\sin\omega_0t/(2m\omega_0)$ である。$\omega\ne\omega_0$ なら同じ初期状態の解は $F_0(\cos\omega t-\cos\omega_0t)/[m(\omega_0^2-\omega^2)]$ であり、$\omega\to\omega_0$ を取ると上の解に収束する。減衰を零に近づける操作と、過渡項を待って消す操作を交換できないのは、待ち時間 $1/\gamma$ が発散するためである。

ばねとダンパーの共通支点 $y$ を動かすと、右辺の力は $ky+b\dot y$ である。実験室系での変位比は $X/Y=(k-ib\omega)\chi$ である一方、相対変位比は $(X-Y)/Y=m\omega^2\chi$ である。例えば $m=1\,\mathrm{kg}$、$k=4\,\mathrm{N\,m^{-1}}$、$b=1\,\mathrm{kg\,s^{-1}}$、$\omega=4\,\mathrm{s^{-1}}$ なら $|X/Y|^2=32/160=1/5$ である。分子に入力の仕方が入るため、同じ機械でも加えた力と支点変位では共鳴曲線が異なる。

### 一定力への立ち上がり

時刻 $0$ に正の一定力 $F_s$ を加え、初期値が零であるとする。力を加えた後の静的平衡は $x_s=F_s/k$ である。$y=x-x_s$ は自由方程式を満たし、$y(0)=-x_s$、$\dot y(0)=0$ である。不足減衰なら

$$
\begin{aligned}
x(t)&=x_s\left[1-e^{-\gamma t}\left(\cos\omega_dt+\frac{\gamma}{\omega_d}\sin\omega_dt\right)\right],\\
\dot x(t)&=x_s\frac{\omega_0^2}{\omega_d}e^{-\gamma t}\sin\omega_dt .
\end{aligned}
$$

したがって最初の最大は $t=\pi/\omega_d$ で、平衡値を超えた割合は $e^{-\pi\gamma/\omega_d}$ である。臨界減衰では $\dot x=x_s\omega_0^2te^{-\omega_0t}\ge0$、過減衰では $\dot x=x_s(\omega_0^2/d)e^{-\gamma t}\sinh dt\ge0$ だから行き過ぎない。調和力での共鳴条件は定常変位の周波数微分で決まる。一方、段階入力での行き過ぎは $x-x_s$ が符号を変えるかで決まる。

:::example[ばねを急に押したときの最大変位]
$m=0.20\,\mathrm{kg}$、$\omega_0=10\,\mathrm{s^{-1}}$、$\gamma=6\,\mathrm{s^{-1}}$、$F_s=0.060\,\mathrm N$ とする。$k=20\,\mathrm{N\,m^{-1}}$、$x_s=3.0\,\mathrm{mm}$、$\omega_d=8\,\mathrm{s^{-1}}$ である。最初の最大変位は $(3.0\,\mathrm{mm})(1+e^{-3\pi/4})\simeq3.28\,\mathrm{mm}$ で、時刻は約 $0.393\,\mathrm s$ である。同じ力をゆっくり加えて過渡振動が十分小さければ、変位はほぼ $3.0\,\mathrm{mm}$ に留まる。
:::

### 一定角振動数で力を入れ始める

任意の駆動周波数での過渡応答は、定常解から初期値を引くと整理できる。$F(t)=F_0\cos\omega t$ を時刻 $0$ から印加し、$X=F_0\chi(\omega)$、$x_s(t)=\Re(Xe^{-i\omega t})$ とする。初期に静止なら自由項の開始値が $-\Re X$、初速度が $-\omega\Im X$ である。共鳴中心 $\omega=\omega_0$ では $X=iA$、$A=F_0/(2m\gamma\omega_0)$ なので不足減衰の全解は

$$
x(t)=A\left[\sin\omega_0t-\frac{\omega_0}{\omega_d}e^{-\gamma t}\sin\omega_dt\right].
$$

式を時間微分すると $t=0$ で二項が打ち消し合う。弱減衰で $\gamma t=O(1)$ なら $\omega_0-\omega_d=O(\gamma^2/\omega_0)$ より二つの位相差は $O(\gamma/\omega_0)$ であり、包絡は $A(1-e^{-\gamma t})$ に近い。共鳴振幅が $\gamma^{-1}$ に比例して増えると同時に、振幅が定常値へ立ち上がる時間も $\gamma^{-1}$ に比例して延びる。

### 振動絶縁で変位比を設計する

質量の絶対変位を $x$、支点の変位を $y$ とする。ばねとダンパがともにこの二点を結ぶなら、$m\ddot x=-k(x-y)-b(\dot x-\dot y)$ である。したがって支点の複素振幅 $Y$ から質量の複素振幅 $X$ への比は

$$
\frac{X}{Y}=\frac{k-i b\omega}{k-m\omega^2-i b\omega}.
$$

同じ問題を相対変位 $z=x-y$ で解くと $m\ddot z+b\dot z+kz=-m\ddot y$ なので、$Z/Y=m\omega^2/(k-m\omega^2-i b\omega)$ となる。この比に $1$ を足すと絶対変位比を回復する。低周波では $x$ が台に追随する一方で、高周波では振幅比 $\abs{X/Y}$ が小さくなり、相対変位はほぼ $-y$ になる。

分子と分母の絶対値の二乗の差は $k^2-(k-m\omega^2)^2=m\omega^2(2k-m\omega^2)$ である。従って $\omega>0$ で $\abs{X/Y}<1$ となる条件は $\omega>\sqrt2\omega_0$ である。また $b>0$ を固定して $\omega\to\infty$ とすると $\abs{X/Y}\sim b/(m\omega)$ である。$b=0$ のときには $k/(m\omega^2)$ になる。ダンパを強めると共鳴域の増幅は抑えられるが、高周波では台の速度がダンパを介して質量に力を及ぼす。

:::example[絶縁が始まる周波数]
$m=2\,\mathrm{kg}$、$k=200\,\mathrm{N/m}$、$b=10\,\mathrm{kg/s}$ では $\omega_0=10\,\mathrm{s^{-1}}$ である。$\omega=10\,\mathrm{s^{-1}}$ での比の大きさは $\sqrt{200^2+100^2}/100=\sqrt5$ であり、$\omega=20\,\mathrm{s^{-1}}$ では $\sqrt{200^2+200^2}/\sqrt{600^2+200^2}=1/\sqrt5$ になる。境界の角振動数は $10\sqrt2\,\mathrm{s^{-1}}$ である。
:::

## 共鳴曲線と Q 値・線幅

### 変位と吸収で異なるピーク

ここでは $\gamma>0$ とし、非零の外力振幅を一定にして周波数を掃引する。

:term[共鳴曲線]{en="resonance curve" yomi="きょうめいきょくせん"} は、何を描くかによって形が変わる。変位振幅を求めるときは、:eqref[eq-wave-osc-susceptibility] の分母の絶対値二乗を $m^2$ で割った $D=(\omega_0^2-\omega^2)^2+4\gamma^2\omega^2$ を最小にする。$D'=4\omega(\omega^2-\omega_0^2+2\gamma^2)$ だから、$\gamma<\omega_0/\sqrt2$ なら $\omega_r=\sqrt{\omega_0^2-2\gamma^2}$ が正の変位ピークである。$\gamma\ge\omega_0/\sqrt2$ では $D$ が正周波数で単調増加し、最大変位は静的な値になる。不足減衰の存在条件 $\gamma<\omega_0$ とは異なる。

吸収電力は後で求めるように $\omega^2/D$ に比例する。$\omega>0$ で $D/\omega^2=(\omega_0^2/\omega-\omega)^2+4\gamma^2$ となるので、最大は厳密に $\omega_0$ である。半値の方程式は $|\omega_0^2-\omega^2|=2\gamma\omega$ である。低周波側では $\omega^2+2\gamma\omega-\omega_0^2=0$、高周波側では $\omega^2-2\gamma\omega-\omega_0^2=0$ となり、正の根は

$$
\omega_\pm=\sqrt{\omega_0^2+\gamma^2}\pm\gamma,\qquad
\Delta\omega=\omega_+-\omega_-=2\gamma,\qquad
Q=\frac{\omega_0}{2\gamma}.
\label{eq-wave-osc-width}
$$

$\Delta\omega$ が吸収の :term[半値全幅]{en="full width at half maximum" yomi="はんちぜんはば"}、$Q$ が :term[品質係数]{en="quality factor" yomi="ひんしつけいすう"} である。半値の中点は一般に $\omega_0$ より高く、幾何平均は $\sqrt{\omega_+\omega_-}=\omega_0$ になる。

### 弱減衰の近似と自由減衰からの測定

$\omega=\omega_0+\eta$ として $|\eta|,\gamma\ll\omega_0$ の範囲を拡大すると、$\omega_0^2-\omega^2=-2\omega_0\eta-\eta^2$ である。吸収をピークで規格化した曲線は $\gamma^2/(\eta^2+\gamma^2)$ に近づく。曲線が半値になる点の間隔は $2\gamma$、振幅の二乗のピークも同じ近似で扱える。変位振幅そのものを半値にすると、規格化した二乗が $1/4$ であり、幅はこの近似で $2\sqrt3\gamma$ となる。

不足減衰の連続する同符号の極大の間隔は $T_d=2\pi/\omega_d$ で、極大の絶対値は毎回 $e^{-\gamma T_d}$ 倍になる。これは波形を一周期ずらすと正弦部分が同じであることから従う。:term[対数減衰率]{en="logarithmic decrement" yomi="たいすうげんすいりつ"} $\Lambda=\log(A_j/A_{j+1})$ とすると $\gamma=\Lambda/T_d$、$\omega_0^2=(2\pi/T_d)^2+\gamma^2$ である。

:::example[線幅から時定数を求める]
吸収を測ると $f_+-f_-=4.0\,\mathrm{Hz}$、ピーク $f_0=200\,\mathrm{Hz}$ であった。式 :eqref[eq-wave-osc-width] から $Q=50$、$\gamma=\pi(4.0)\,\mathrm{s^{-1}}$ である。自由振動の振幅が $1/e$ になる時間は $1/\gamma\simeq0.080\,\mathrm{s}$ となる。
:::

:::example[目に見える自由振動と平坦な応答]
$\omega_0=5\,\mathrm{s^{-1}}$、$\gamma=4\,\mathrm{s^{-1}}$ では $\omega_d=3\,\mathrm{s^{-1}}$ である。これに対し $\omega_0^2-2\gamma^2=-7\,\mathrm{s^{-2}}$ なので変位曲線に正周波数のピークはない。吸収の半値は $\sqrt{41}\pm4\,\mathrm{s^{-1}}$ で、幅は $8\,\mathrm{s^{-1}}$ である。
:::

### 変位振幅の二乗で定めた場合の線幅

吸収の線幅は式 :eqref[eq-wave-osc-width] である。変位を使う測定では、別に分母の最小値を求める必要がある。$\gamma<\omega_0/\sqrt2$ で $D=(\omega^2-\omega_r^2)^2+4\gamma^2(\omega_0^2-\gamma^2)$ と平方完成できる。従って $A^2$ がピークの半値になる二つの二乗周波数の候補は

$$
\omega^2=\omega_0^2-2\gamma^2\pm 2\gamma\sqrt{\omega_0^2-\gamma^2}.
$$

下側が正である条件は $\gamma^2/\omega_0^2<(2-\sqrt2)/4$ である。実際 $\omega_0^2>2\gamma^2$ のもとで不等式を二乗し、$8(\gamma^2/\omega_0^2)^2-8(\gamma^2/\omega_0^2)+1>0$ として小さい根を選べばよい。共鳴ピークが正の周波数にあっても、低い周波数へ下げたとき静的な変位二乗が大きければ半値に届かない。吸収では $\omega\to0$ で電力が零になるため、二つの半値が常にある。

:::example[変位と吸収の幅を比較する]
$\omega_0=10\,\mathrm{s^{-1}}$、$\gamma=3\,\mathrm{s^{-1}}$ なら、変位二乗の半値の二乗周波数は $82\pm6\sqrt{91}$ である。幅は $\sqrt{82+6\sqrt{91}}-\sqrt{82-6\sqrt{91}}\simeq6.82\,\mathrm{s^{-1}}$ となる。吸収の幅は $6.00\,\mathrm{s^{-1}}$ である。さらに $\gamma=4\,\mathrm{s^{-1}}$ まで上げると $68-8\sqrt{84}<0$ となり、変位の正周波数ピークはまだあるのに、その低周波側に半値がなくなる。
:::

### 線幅の逆問題と低周波側の半値

吸収の半値を二点測れば、ピークを直接走査しなくても $\omega_0=\sqrt{\omega_+\omega_-}$ と $\gamma=(\omega_+-\omega_-)/2$ が決まる。式 :eqref[eq-wave-osc-width] の積と差を取ればよい。この推定には外力振幅を周波数にわたり一定にするか、その二乗で吸収を割る必要がある。駆動装置の周波数応答を含む電力曲線をそのまま使うと、半値の位置が移る。

変位二乗に対しては低い方の半値根が実際に存在する条件を求められる。$r=\gamma^2/\omega_0^2$ とすると、ピークが正周波数にある範囲 $0<r<1/2$ でこの条件は $1-2r\ge2\sqrt{r(1-r)}$ である。両辺が非負なので二乗でき、$8r^2-8r+1\ge0$ となる。許される解は $r\le (1-1/\sqrt2)/2$ である。等号では下の半値が $\omega=0$ となり、それより減衰が大きいと正周波数側の交点が一つになる。

:::example[ピークがあっても下の半値がない曲線]
$\omega_0=5\,\mathrm{s^{-1}}$、$\gamma=2\,\mathrm{s^{-1}}$ では $r=4/25$ である。変位二乗のピークは $\omega_r^2=17\,\mathrm{s^{-2}}$ にあり、下側の候補は $\omega^2=17-4\sqrt{21}<0$ である。一方、吸収には $\omega_\pm=\sqrt{29}\pm2\,\mathrm{s^{-1}}$ の二点がある。観測した縦軸の定義を確かめれば、同じ装置でこの違いが生じることを説明できる。
:::

## 応答関数と因果性

### 力積から任意の力へ

系を時刻 $0$ より前に静止させ、時刻 $0$ に力積 $J$ を与える。運動方程式を時刻 $0$ の前後で積分すると $m[\dot x]_{0-}^{0+}=J$ となる。有限の力積で変位は跳ばない。もし変位に跳びがあれば加速度にデルタ関数の微分が生じ、それを右辺のデルタ関数だけで打ち消せないからである。単位力積に対する変位 $G$ を不足減衰の場合に書くと

$$
G(t)=H(t)\frac{e^{-\gamma t}\sin\omega_dt}{m\omega_d},
\qquad
x(t)=x_{\mathrm{free}}(t)+\int_0^t G(t-s)F(s)\dd s.
\label{eq-wave-osc-duhamel}
$$

$H$ は $t<0$ で零、$t>0$ で $1$ の階段関数である。$G$ は :term[遅延 Green 関数]{en="retarded Green function" yomi="ちえんぐりーんかんすう"} と呼ばれる。臨界減衰では $\sin\omega_dt/\omega_d$ を $t$、過減衰では $\sinh dt/d$ に置き換える。$G$ の単位は $\mathrm{m/(N\,s)}$ であり、畳み込みで変位に戻る。

畳み込みの微分で $G(0+)=0$ を使うと、一階微分は $\int_0^tG'(t-s)F(s)\dd s$ になる。二階微分の上端項は $G'(0+)F(t)=F(t)/m$ である。積分内では $G$ が同次方程式を満たすため、作用素を掛けた結果はちょうど $F(t)$ になる。初期値が零の特解がこれで得られ、一般の初期値には自由振動を足す。

### 解析性と実部・虚部の関係

$\gamma>0$ とする。時間を $e^{-i\omega t}$ で展開する規約では

$$
\int_0^\infty G(t)e^{izt}\dd t
=\frac{1}{m(\omega_0^2-z^2-2i\gamma z)}
$$

である。不足減衰では正弦を二つの指数に分けて積分でき、三つの場合に共通する導出には二回の部分積分を使える。すなわち $\int_0^\infty G''e^{izt}\dd t=-1/m-z^2\int_0^\infty Ge^{izt}\dd t$、$\int_0^\infty G'e^{izt}\dd t=-iz\int_0^\infty Ge^{izt}\dd t$ を同次方程式へ入れる。$z$ が上半平面にあると $e^{izt}$ が指数減衰するため、コンパクト部分集合では $G(t)e^{izt}$ とその $z$ 微分を可積分関数で支配できる。従って $\chi(z)$ は上半平面で正則である。極は不足減衰で $z=\pm\omega_d-i\gamma$、過減衰で $z=-i(\gamma\mp d)$ にあり、いずれも下半平面にある。

実軸上の応答率を $\chi'+i\chi''$ と書く。$\chi(-\omega)=\conj{\chi(\omega)}$ は $G$ が実数値であることから従い、$\chi'$ は偶、$\chi''$ は奇である。さらにこの有理関数は無限遠で $O(|z|^{-2})$ なので、:ref[thm-cpx-cauchy-formula]{full} を用いて $z=\omega+i\eta$ に対し $\chi(z)=(2\pi i)^{-1}\int_{\R}\chi(u)/(u-z)\dd u$ と書ける。$\eta\downarrow0$ で核を $(u-\omega)/((u-\omega)^2+\eta^2)+i\eta/((u-\omega)^2+\eta^2)$ と分ける。前者は $\omega$ の近くで定数項を除くと主値に収束する。後者は幅 $\eta$ の近くに積分値 $i\pi$ を集中させるので $i\pi\chi(\omega)$ を与える。これを左辺に移すと

$$
\chi'(\omega)=\frac1\pi\,\mathrm{PV}\int_{\R}\frac{\chi''(u)}{u-\omega}\dd u,\qquad
\chi''(\omega)=-\frac1\pi\,\mathrm{PV}\int_{\R}\frac{\chi'(u)}{u-\omega}\dd u.
\label{eq-wave-osc-kk}
$$

ここで $\mathrm{PV}$ は特異点の両側から等しい長さの区間を除いて極限を取ることを意味する。これらを :term[Kramers–Kronig 関係]{en="Kramers–Kronig relations" yomi="くらまーすくろーにひかんけい"} という。積分が成立するためには上半平面での正則性と無限遠の減衰を使っている。別の系で瞬時応答 $\chi_\infty$ が残るなら、その定数を引いた応答に対して同じ導出をする。

:::example[一周期だけ力を掛ける]
無減衰の振動子を静止させ、$0<s<D$ で一定力 $F$ を加える。$t>D$ では :eqref[eq-wave-osc-duhamel] の積分範囲が $0$ から $D$ となり、結果は $(F/k)[\cos\omega_0(t-D)-\cos\omega_0t]$ である。$D=2\pi/\omega_0$ なら二項が打ち消し合い、切った後は静止する。$D=\pi/\omega_0$ なら残る振幅は $2|F|/k$ である。
:::

:::example[低周波応答の符号]
式 :eqref[eq-wave-osc-kk] の $\omega=0$ を見ると、$\chi'(0)=(2/\pi)\int_0^\infty\chi''(u)/u\,\dd u$ である。粘性減衰の振動子では $\chi''(u)>0$ であり、これは静的に押した向きへ動く $\chi'(0)=1/k>0$ と整合する。
:::

### 無減衰応答と吸収の極限

$\gamma\downarrow0$ では、実数の周波数が固有振動数から離れていれば $\chi''$ は零へ近づく。しかし共鳴に沿う領域を積分すると、一般に零にならない。正周波数の共鳴近くで $\omega=\omega_0+\gamma s$ とおくと

$$
\gamma\chi''(\omega_0+\gamma s)\longrightarrow
\frac{1}{2m\omega_0(1+s^2)}
$$

である。固定した小区間 $|\omega-\omega_0|<\omega_0/2$ では、$|\omega_0^2-\omega^2|=|\omega_0+\omega|\,|\omega-\omega_0|$ と $\omega\ge\omega_0/2$ を使い、左辺を $C/(1+s^2)$ で抑えられる。区間外の有界領域では分母が零から離れているので積分は零へ近づき、十分大きい周波数でも $O(\gamma/\omega^3)$ で抑えられる。試験関数を掛けて積分した極限は、正周波数で $\pi/(2m\omega_0)$ を重みとするデルタ関数である。負周波数にある共鳴の重みは $\chi''$ の奇性で符号が逆になる。無減衰の式に通常の関数として虚部がないことと、弱く減衰した極限の積分が残ることを両立させるには、どの極限を取るかを指定する必要がある。

### 運動方程式を Laplace 変換する

$t\ge0$ の Laplace 変換を $\widetilde x(s)=\int_0^\infty e^{-st}x(t)\dd t$ とする。収束する右半平面で部分積分すると

$$
m(s^2+2\gamma s+\omega_0^2)\widetilde x
=\widetilde F+m[(s+2\gamma)x_0+v_0]
$$

である。分母が作る極は自由振動の指数である。入力と初期値のどちらから生じた項であるかは分子で区別できる。$s=-iz$ とおけば $1/[m(s^2+2\gamma s+\omega_0^2)]$ が式 :eqref[eq-wave-osc-susceptibility] へ移る。Laplace 領域の右半平面と Fourier 領域の上半平面が同じ因果的な積分に対応する。

:::example[二つの短い押し方を区別する]
同じ無減衰振動子へ、時刻 $0$ に力積 $J$、時刻 $D$ に反対向きの力積 $-J$ を与える。$t>D$ の解は $J[\sin\omega_0t-\sin\omega_0(t-D)]/(m\omega_0)$ で、残留振幅は $2|J\sin(\omega_0D/2)|/(m\omega_0)$ である。半周期で二つの力積を与えると、総外力積は零でも残留振幅は単発の二倍となる。最初の入力による速度が半周期後に逆向きに回っており、二つ目の力積がこの速度を増すためである。
:::

### 線形応答の周波数積分と短時間応答

式 :eqref[eq-wave-osc-kk] は実軸上の情報から静的剛性を求める一例を与える。さらに $z\chi(z)$ は上半円の上で $-1/(mz)+O(\abs{z}^{-2})$ となる。上半平面で極を持たないので、半径 $R$ の上半円弧と $[-R,R]$ に沿う積分の和は零である。弧を $z=Re^{i\theta}$、$0\le\theta\le\pi$ と置くと、その積分は $-i\pi/m+O(1/R)$ になる。実軸の積分で虚部を取り、$\omega\chi''(\omega)$ が偶関数で $O(\omega^{-2})$ であることを使うと

$$
\int_0^\infty \omega\chi''(\omega)\dd\omega=\frac{\pi}{2m}.
$$

この関係を :term[総和則]{en="sum rule" yomi="そうわそく"} と呼ぶ。無減衰に近いと吸収は狭い帯域に集中するが、同じ質量ならこの積分は変化しない。逆に、単一の粘性抵抗モデルが測定に適合するなら、積分値を質量の見積もりに使える。周波数を $\mathrm{Hz}$ で積分するときは $\omega=2\pi f$ を微分にも適用する。

:::example[臨界減衰で積分を直接確かめる]
$\gamma=\omega_0=a>0$ では $\chi''(\omega)=2a\omega/[m(a^2+\omega^2)^2]$ である。$\omega=a\tan\theta$ とすると $\omega\chi''(\omega)\dd\omega=(2/m)\sin^2\theta\dd\theta$ となる。$0\le\theta<\pi/2$ で積分すると $\pi/(2m)$ である。例えば $m=2\,\mathrm{kg}$ なら積分は $\pi/4\,\mathrm{kg^{-1}}$ になる。
:::

## エネルギー吸収と散逸

### 瞬時収支と周期平均

運動方程式に $\dot x$ を掛けると、運動とばねのエネルギー $E=m\dot x^2/2+kx^2/2$ は

$$
\dot E=F\dot x-b\dot x^2
\label{eq-wave-osc-energy}
$$

を満たす。右辺の第 1 項は外力がした仕事率、第 2 項の絶対値は抵抗による散逸率であり、どちらも $\mathrm{W}$ である。定常調和解では周期平均 $\langle\dot E\rangle=0$ なので、平均入力と平均散逸が等しい。$\Re(Ae^{-i\omega t})\Re(Be^{-i\omega t})$ を一周期積分すると、周波数 $2\omega$ の項が消え、平均は $\Re(A\conj B)/2$ である。$\dot x$ の複素振幅を $V=-i\omega X$ とすると

$$
\langle P\rangle
=\frac12\Re(F_0\conj V)
=\frac{\omega|F_0|^2}{2}\chi''(\omega)
=\frac{b\omega^2|X|^2}{2}.
\label{eq-wave-osc-power}
$$

この符号規約では受動的な系で $\omega\chi''(\omega)\ge0$ となる。速度は変位より $\pi/2$ 進むので、力と変位が $\pi/2$ ずれる $\omega=\omega_0$ では力と速度が同相となる。動作点の複素位相を変えると、外力と速度の積に負の瞬間が現れる。それでも定常一周期で正味の仕事は非負である。

### 速度インピーダンスと一周期の仕事

非零の調和速度に対する力の比を $Z_{\mathrm m}=F_0/V$ と定義すると

$$
Z_{\mathrm m}=b+i\left(\frac{k}{\omega}-m\omega\right).
$$

実部 $b$ は散逸、虚部はばねと慣性の間の可逆なエネルギー交換に対応する。定常振動の $F$ を $(kx+m\ddot x)+b\dot x$ と分けて一周期積分すると、前者の仕事は端点のエネルギー差で零であり、後者の仕事は $b\int\dot x^2\dd t$ である。従って $F$ 対 $x$ の閉軌道の向き付き面積は $\oint F\dd x=\pi b\omega |X|^2$ となる。

共鳴時には速度の最大時と変位の最大時のエネルギーが等しく、$E=m\omega_0^2|X|^2/2$ は周期内で一定である。$W_{\mathrm{cyc}}=\langle P\rangle(2\pi/\omega_0)$ とすれば $2\pi E/W_{\mathrm{cyc}}=m\omega_0/b=Q$ になる。自由減衰で一周期あたりの損失を使う定義は、弱減衰でこれと一致する。実際、周期の端点で全エネルギーが $e^{-2\gamma T_d}$ 倍となり、損失率 $1-e^{-2\gamma T_d}$ は $2\gamma T_d$ に近い。

### 直列回路で電流の共鳴を測る

インダクタンス $L_{\mathrm e}$、抵抗 $R_{\mathrm e}$、静電容量 $C_{\mathrm e}$ を直列につなぐ。コンデンサーの電荷を $q$、電流を $I=\dot q$ とすると、各素子の電圧の和は $L_{\mathrm e}\ddot q+R_{\mathrm e}\dot q+q/C_{\mathrm e}$ であり、電源電圧に等しい。従って機械系の $(m,b,k,x,F)$ はそれぞれ $(L_{\mathrm e},R_{\mathrm e},1/C_{\mathrm e},q,V)$ に対応する。電荷が変位に、電流が速度に対応するので、電流の振幅を測ると吸収のピーク位置が得られる。電圧と電流の複素振幅の比は $R_{\mathrm e}+i(1/(\omega C_{\mathrm e})-\omega L_{\mathrm e})$ である。また電源を零電圧にすると、エネルギー $L_{\mathrm e}I^2/2+q^2/(2C_{\mathrm e})$ は $R_{\mathrm e} I^2$ の割合で減る。この対応は :cite[crawford-waves] の電気振動にも現れる。

:::example[抵抗を変えて帯域を調べる]
$L_{\mathrm e}=10\,\mathrm{mH}$、$C_{\mathrm e}=1.0\,\mu\mathrm F$、$R_{\mathrm e}=10\,\Omega$ なら、電流の共鳴は $\omega_0=10^4\,\mathrm{s^{-1}}$、電力の半値全幅は $R_{\mathrm e}/L_{\mathrm e}=10^3\,\mathrm{s^{-1}}$ で、$Q=10$ である。電源の電圧振幅が $2.0\,\mathrm V$ なら、共鳴時に流れる電流の振幅は $0.20\,\mathrm A$、抵抗に入る平均電力は $0.20\,\mathrm W$ になる。同じ回路でもコンデンサー両端の電圧は電荷に比例するので、そのピークには変位の共鳴条件を使う。
:::

:::example[平均電力と一周期の仕事]
$m=1\,\mathrm{kg}$、$b=2\,\mathrm{kg\,s^{-1}}$、$\omega_0=4\,\mathrm{s^{-1}}$ で、共鳴時に $F_0=4\,\mathrm N$ を加える。速度振幅は $F_0/b=2\,\mathrm{m\,s^{-1}}$ だから平均電力は $4\,\mathrm W$ である。一周期の仕事は $2\pi\,\mathrm J$ であり、$E=2\,\mathrm J$ から求める $Q$ は $2$ である。
:::

:::example[変位を測って吸収を推定する]
$b=0.4\,\mathrm{kg\,s^{-1}}$、$\omega=50\,\mathrm{s^{-1}}$、$A=3\,\mathrm{mm}$ なら速度振幅は $0.15\,\mathrm{m\,s^{-1}}$ となる。定常吸収は $b(0.15)^2/2=4.5\,\mathrm{mW}$ である。実効速度を使う測定なら $v_{\mathrm{rms}}=0.15/\sqrt2\,\mathrm{m\,s^{-1}}$ なので $\langle P\rangle=bv_{\mathrm{rms}}^2$ と書ける。
:::

### 過渡外力が残すエネルギー

実数値の力 $F$ を有限の時刻区間 $[0,D]$ に加え、減衰のない振動子を静止から動かす。力が消えた後、式 :eqref[eq-wave-osc-duhamel] の $\sin\omega_0(t-s)$ を二つの三角関数の積へ分けると

$$
E(t>D)=\frac1{2m}\left|\int_0^D F(s)e^{i\omega_0s}\dd s\right|^2 .
$$

実際、$\sin\omega_0t$ の係数が $\int F\cos\omega_0s/(m\omega_0)$、$\cos\omega_0t$ の係数が $-\int F\sin\omega_0s/(m\omega_0)$ である。エネルギーの二乗和にすると $m\omega_0^2$ が打ち消される。力積は $\omega=0$ における同じ積分であり、エネルギーを決めるのは固有振動数での値である。このため、時間波形を設計して任意の有界な正負の力を使うときは、力の総和だけで残留振動を予測できない。

:::example[半周期だけの一定力]
$m=0.50\,\mathrm{kg}$、$\omega_0=20\,\mathrm{s^{-1}}$、$F_s=0.20\,\mathrm N$ として $D=\pi/\omega_0$ の間に力を加える。積分の絶対値は $2F_s/\omega_0=0.020\,\mathrm{N\,s}$ なので、残るエネルギーは $4.0\times10^{-4}\,\mathrm J$ である。力を同じ大きさで一周期加える場合はこの積分が零になる。前半と後半で複素位相が逆向きになり、寄与を打ち消すためである。
:::

### 変位と力をグラフにしたときの仕事

過渡を含む測定の $F(t)$ と $x(t)$ からは、時間によらず式 :eqref[eq-wave-osc-energy] を積分できる。特に $F_0$ と $\omega$ が固定された周期応答 $x=A\cos(\omega t-\delta)$ で、横軸を $x$、縦軸を $F$ とした曲線を得たとする。$\theta=\omega t-\delta$ と置くと

$$
x/A=\cos\theta,\qquad
F/F_0=\cos\theta\cos\delta-\sin\theta\sin\delta.
$$

$\sin\delta\ne0$ なら、第二式から第一式の $\cos\delta$ 倍を除いて二乗すると楕円の式になる。仕事は $F\dd x=-AF_0(\cos\theta\cos\delta-\sin\theta\sin\delta)\sin\theta\dd\theta$ を一周積分した $\pi A F_0\sin\delta$ である。力と変位が同相または逆相なら楕円は線分に縮み、囲む面積は零になる。粘性抵抗では $F_0\sin\delta=b\omega A$ なので、どの周波数でも $\pi b\omega A^2$ に一致する。例えば $F_0=2.0\,\mathrm N$、$A=1.0\,\mathrm{mm}$、$\delta=\pi/6$ なら一周期の仕事は $\pi\times10^{-3}\,\mathrm J$ であり、$100\,\mathrm{Hz}$ での平均電力は約 $0.314\,\mathrm W$ である。この測り方では平均電力の符号と位相遅れの符号を同時に確認できる。

$F_0$ と $A$ の単位を SI で掛けた結果は $\mathrm{N\,m}=\mathrm J$ であり、周波数を掛けると $\mathrm W$ になる。仕事を楕円の面積として測る場合も、各軸の目盛をこの単位に換算する。

### 複数周波数が作る散逸

線形かつ時間に依存しない系では、入力を角振動数ごとに分けて応答を求められる。ここでは $\omega_j=n_j\Omega$、相異なる正整数 $n_j$ の有限和 $F(t)=\sum_j\Re(F_j e^{-in_j\Omega t})$ を用いる。過渡が消えた後の速度は $v(t)=\sum_j\Re(V_j e^{-in_j\Omega t})$、$V_j=-in_j\Omega\chi(n_j\Omega)F_j$ である。$v^2$ を展開すると、各積は周波数 $(n_j+n_\ell)\Omega$ または $(n_j-n_\ell)\Omega$ をもつ。共通周期 $2\pi/\Omega$ で積分すると、$n_j=n_\ell$ の差周波数だけが残る。したがって

$$
\langle Fv\rangle=b\langle v^2\rangle
=\frac{b}{2}\sum_j\abs{V_j}^2 .
$$

同じ角振動数の二つの力は、先に複素振幅を足す。そうすると、相対位相によって仕事が増減することも同じ式に含まれる。過渡を含む時間区間では $\int Fv\dd t=\Delta E+b\int v^2\dd t$ を使う。

:::example[共鳴に別の調波を加える]
$m=1\,\mathrm{kg}$、$k=4\,\mathrm{N/m}$、$b=1\,\mathrm{kg/s}$ の系に $F=2\cos(2t)+3\cos(4t)\,\mathrm N$ を加える。$\omega=2\,\mathrm{s^{-1}}$ の速度振幅は $2\,\mathrm{m/s}$ である。$\omega=4\,\mathrm{s^{-1}}$ では $\abs{V}^2=12^2/(12^2+4^2)=9/10\,\mathrm{m^2/s^2}$ になる。従って平均散逸は $2+9/20=2.45\,\mathrm W$ である。二つの位相をずらしても各成分の二乗は変わらない。
:::

## パラメトリック共振

### 剛性を周期的に変える

外力がなくても、剛性が時間によって変わると振動へ仕事を供給できる。$k(t)=m\omega_0^2[1+h\cos(2\nu t)]$ と置いた

$$
\ddot x+2\gamma\dot x+\omega_0^2[1+h\cos(2\nu t)]x=0
$$

は減衰を含む :term[Mathieu 方程式]{en="Mathieu equation" yomi="まちゅーほうていしき"} である。$h$ は無次元であり、$|h|\ll1$ を仮定する。剛性が時間に依存する状態で力学的エネルギーを微分すると $\dot E=\dot k(t)x^2/2-b\dot x^2$ となる。外力で直接引く場合と異なり、$x=0$、$\dot x=0$ はいつでも厳密な解である。不安定性とは、そこからの微小なずれが増幅されることを指す。

### 主不安定帯の導出

$\gamma/\nu$、$|h|$、$|\omega_0-\nu|/\nu$ を同程度の微小量とする。実変位を $x=a(t)e^{i\nu t}+\conj{a(t)}e^{-i\nu t}$ と表す。$a$ は時間 $1/\nu$ でほとんど変わらないとして $\ddot a$ と $\gamma\dot a$ を一次近似で落とす。$h\cos2\nu t$ の積では、$\conj a e^{-i\nu t}$ に $e^{2i\nu t}$ を掛けた項が $e^{i\nu t}$ になり、同じ位相で蓄積する。$e^{\pm3i\nu t}$ の項は一周期平均で消え、一次のゆっくりした振幅方程式に寄与しない。そのため

$$
\begin{aligned}
2i\nu\dot a+(2i\gamma\nu+\omega_0^2-\nu^2)a
+\frac{h\omega_0^2}{2}\conj a&=0,\\
\dot a&=(-\gamma+i\Delta)a+ip\conj a,\qquad
\Delta=\frac{\omega_0^2-\nu^2}{2\nu},\quad p=\frac{h\omega_0^2}{4\nu}.
\end{aligned}
$$

共役の方程式も書くと固有指数の特性式は $(\lambda+\gamma)^2+\Delta^2-p^2=0$ である。成長率が正になる一次近似の条件は $p^2>\Delta^2+\gamma^2$ となる。中心 $\nu=\omega_0$ では $|h|>4\gamma/\omega_0$ が閾値であり、無減衰ではどんな小さな非零変調でも中心付近で成長する。これは近似式の条件である。高い次数や別の帯を調べるには、変調一周期 $T=\pi/\nu$ の基本解行列を求める。一周期の行列が $S$ なら周期性から $n$ 周期で $S^n$ になるため、固有値の絶対値が $1$ を超えるかどうかで指数的な不安定性を判定できる。絶対値 $1$ の固有値では、Jordan ブロックがあると多項式成長が起こりうる。

:::example[閾値と離調]
$\omega_0=100\,\mathrm{s^{-1}}$、$\gamma=1\,\mathrm{s^{-1}}$ なら中心で閾値は $0.04$ である。$h=0.08$ なら中心で $p=2\,\mathrm{s^{-1}}$ となり、成長率は $1\,\mathrm{s^{-1}}$ である。一次の離調幅は $|\Delta|<\sqrt3\,\mathrm{s^{-1}}$ である。振幅が大きくなると、線形の剛性変調モデルだけでは運動を定められない場合がある。
:::

:::example[固有の成長方向と減衰方向]
中心で $p>0$ とし $a=u+iv$ と分けると $\dot u=-\gamma u+pv$、$\dot v=pu-\gamma v$ である。$u=v$ の成分は $p-\gamma$、$u=-v$ の成分は $-p-\gamma$ で変化する。同じ変調を受けても、外乱が最初にどの位相を持つかで過渡的な増減が異なる。
:::

### 一周期に供給されるエネルギーと変調位相

弱い変調のもとで $x\simeq A\cos(\nu t+\phi)$ と近似する。$\dot k=-2\nu m\omega_0^2h\sin2\nu t$、$x^2=A^2[1+\cos(2\nu t+2\phi)]/2$ であり、時刻平均は

$$
\left\langle\frac{\dot k\,x^2}{2}\right\rangle
=\frac{m\omega_0^2h\nu A^2}{4}\sin2\phi,\qquad
\langle b\dot x^2\rangle=\frac{b\nu^2 A^2}{2}
$$

となる。$\sin2\nu t\cos(2\nu t+2\phi)$ の周期平均が $-\sin2\phi/2$ であることを使った。中心で $h>0$ の場合、$\phi=\pi/4$ でエネルギーの供給が最大、$\phi=-\pi/4$ で供給が負になる。この位相は主不安定帯の固有方向に一致する。中心で最大のとき $E\simeq m\omega_0^2A^2/2$ と置けば、$A'/A=h\omega_0/4-\gamma$ となり、複素振幅からの計算を再現する。

### 周期写像の行列式

変調を含む減衰の方程式で、二つの実解を $x_1,x_2$ とする。$W=x_1\dot x_2-x_2\dot x_1$ を微分すると、復元力の項が打ち消し合って $\dot W=-2\gamma W$ である。従って一周期の状態行列の行列式は $e^{-2\gamma T}$ である。無減衰なら二つの乗数は相反数であり、減衰が正なら乗数の積は $1$ 未満となる。実係数の二次の方程式で複素共役の乗数がある場合、その絶対値は $e^{-\gamma T}$ となる。片方が指数成長するには、それが実数の乗数で絶対値 $1$ を超える必要がある。変調周波数が $2\omega_0$ 近くでは、時間因子 $e^{i\nu T}=-1$ を含むため、成長率が正の場合にも主帯の乗数は負になりうる。乗数の符号は一周期後の位相、絶対値は振幅の倍率を表す。

### 零解の安定性を一周期で判定する

周期 $T_m=\pi/\nu$ の線形方程式で二つの独立な初期条件 $(1,0)$ と $(0,1)$ を一周期進め、その列を並べた実行列 $M$ を作る。第 $j$ 周期の方程式は第 $0$ 周期と同じなので、$n$ 周期後の状態は $M^n$ を初期状態に掛けたものになる。この $M$ の固有値 $\lambda$ を :term[Floquet 乗数]{en="Floquet multiplier" yomi="ふろけじょうすう"} という。固有方向の振幅は整数周期ごとに $\abs{\lambda}$ 倍になる。両固有値が単位円内なら、重根で対角化できない場合にも Jordan 形の項 $n\lambda^{n-1}$ が零に収束するので全解が減衰する。$\abs{\lambda}=1$ 上の非対角化可能な重根は多項式成長を許す。

$\gamma>0$ の場合は $d=\det M=e^{-2\gamma T_m}\in(0,1)$ である。$\tau=\tr M$ と書くと特性多項式は $P(\lambda)=\lambda^2-\tau\lambda+d$ である。根が非実なら共役で $\abs{\lambda}^2=d$、従って減衰する。実根なら積が $d>0$ なので同符号であり、二根がともに $1$ より大きい場合や、ともに $-1$ より小さい場合は積が $1$ を超えてしまう。従って実根がともに $(-1,1)$ に入るための必要十分条件は $P(1)>0$ かつ $P(-1)>0$ である。非実根にもこの条件は成り立つので、判定は $\abs{\tau}<1+d$ にまとめられる。等号では $\lambda=\pm1$ ともう一つの絶対値 $d$ の根が現れる。

:::example[周期ごとの符号反転]
ある変調で $d=1/4$、$\tau=-3/2$ とする。$P(-1)=1+\tau+d=-1/4$ なので $-1$ より小さい根がある。二根は $(-3\pm\sqrt5)/4$ であり、成長する方向では周期ごとに符号も反転する。この検査は弱変調の平均化とは独立に使える。
:::

## 非線形振動子（Duffing 振動子と van der Pol 振動子）

### Duffing 振動子の周波数移動

復元力の奇対称性がある系では、線形の次の項として $-m\alpha x^3$ が現れうる。:term[Duffing 振動子]{en="Duffing oscillator" yomi="だふぃんぐしんどうし"} の式を

$$
\ddot x+2\gamma\dot x+\omega_0^2x+\alpha x^3=f\cos\omega t
$$

とする。$f=F_0/m$ は加速度、$\alpha$ の単位は $\mathrm{m^{-2}s^{-2}}$ である。$\gamma=f=0$ なら $m(\dot x^2+\omega_0^2x^2)/2+m\alpha x^4/4$ が保存される。$\alpha<0$ のモデルを大変位まで使うと四次ポテンシャルが下に発散するため、原点付近の振動に範囲を限定する。

振幅 $A$ で $\varepsilon=|\alpha|A^2/\omega_0^2\ll1$ と仮定し、$\theta=\Omega t$、$x=A\cos\theta+y(\theta)$、$y=O(\varepsilon A)$ とする。$x^3$ の一次の寄与は $A^3(3\cos\theta+\cos3\theta)/4$ である。残差の $\cos\theta$ 成分を零にすると、共鳴した補正 $\theta\sin\theta$ が現れず、

$$
\Omega^2=\omega_0^2+\frac34\alpha A^2+O(\varepsilon^2\omega_0^2),\qquad
\Omega=\omega_0+\frac{3\alpha A^2}{8\omega_0}+O(\varepsilon^2\omega_0)
$$

を得る。残る第三高調波について $(\partial_\theta^2+1)\cos3\theta=-8\cos3\theta$ なので、補正は $\alpha A^3(\cos3\theta-\cos\theta)/(32\omega_0^2)$ と取れる。第 2 項は $x(0)=A$ に合わせるためである。正の $\alpha$ は硬化、負の $\alpha$ は軟化を示す。

### 外力がある場合の振幅方程式

弱減衰で駆動が基本周波数に近い場合、$x\simeq A\cos(\omega t-\delta)$ を代入する。三次項の基本波だけを残して余弦と正弦の係数を比較すると

$$
f\cos\delta=A(\omega_0^2-\omega^2+3\alpha A^2/4),\qquad
f\sin\delta=2\gamma\omega A.
$$

従って $B=A^2$ とすると、$B[(\omega_0^2-\omega^2+3\alpha B/4)^2+4\gamma^2\omega^2]=f^2$ という三次方程式になる。たとえば $d=\omega_0^2-\omega^2$、$a=3\alpha/4$、$r=2\gamma\omega$ とおく。左辺の $B$ 微分は $3a^2B^2+4adB+d^2+r^2$ である。停止点が二つあるには判別式 $4a^2(d^2-3r^2)>0$ と、二根がともに正であることが必要である。後者は $ad<0$ と同値である。三つの正根がある場合も、それぞれの安定性は振幅方程式の根を数えただけでは決まらない。過渡振動を含む時間発展や、周期解に沿う変分方程式が要る。

### van der Pol 振動子の振幅選択

ばねの非線形性と別に、振幅に応じて抵抗の符号が変わる系を考える。無次元の変位 $u$ と時刻 $\tau=\omega_0t$ を採用した :term[van der Pol 振動子]{en="van der Pol oscillator" yomi="ふぁんでるぽーるしんどうし"} は

$$
u''-\mu(1-u^2)u'+u=0
$$

である。$\mu>0$ では $|u|<1$ でエネルギーを受け取り、$|u|>1$ では失う。実際 $E=(u'^2+u^2)/2$ と置けば $E'=\mu(1-u^2)u'^2$ である。$0<\mu\ll1$ で $u\simeq R(\tau)\cos(\tau+\phi)$ を使うと、一周期平均は $\langle u'^2\rangle=R^2/2$、$\langle u^2u'^2\rangle=R^4/8$ である。後者は $\cos^2\theta\sin^2\theta=(1-\cos4\theta)/8$ から分かる。$\langle E'\rangle\simeq R R'$ として

$$
R'=\frac{\mu}{2}R\left(1-\frac{R^2}{4}\right)
$$

を得る。$0<R<2$ で振幅が増え、$R>2$ で減る。$s=R^2$ とおけば $s'=\mu s(1-s/4)$ なので、$s(0)>0$ では $s=4/[1+(4/s(0)-1)e^{-\mu\tau}]$ となる。$R\simeq2$ は弱非線形近似の結果である。元の方程式では $\mu>0$ に対して安定な周期軌道が一つ存在する。その証明は :ref[thm-ode-lienard]{full} による。ここでは弱非線形での振幅と、強非線形での周期の先頭近似をそれぞれ求める。

:::example[硬化型の振動数変化]
$\omega_0=10\,\mathrm{s^{-1}}$、$\alpha=10^4\,\mathrm{m^{-2}s^{-2}}$、$A=0.020\,\mathrm{m}$ では $\alpha A^2/\omega_0^2=0.04$ である。一次の角振動数移動は $3\alpha A^2/(8\omega_0)=0.15\,\mathrm{s^{-1}}$ である。振幅を半分にすると移動は $0.0375\,\mathrm{s^{-1}}$ になる。
:::

:::example[自励振動の初期成長]
$\mu=0.1$、$R=1$ での平均変化率は $R'=0.0375$ である。同じ $\mu$ で $R=3$ なら $R'=-0.1875$ である。どちらから出発しても平均化された振幅は $2$ に近づくが、原点に厳密に置いた初期状態は原点に留まる。
:::

### Duffing 振動子の有限振幅の周期積分

外力と粘性抵抗が零で、原点を囲む周期軌道を考える。最大変位を $A>0$ とする。無次元の $\varepsilon=\alpha A^2/\omega_0^2$ に対し、軌道が原点の谷の中にある十分条件は $\varepsilon>-1$ である。エネルギー保存式より、$0\le x\le A$ での速さの二乗は

$$
\dot x^2=(A^2-x^2)\left[\omega_0^2+\frac{\alpha}{2}(A^2+x^2)\right].
$$

周期の四分の一は原点から $A$ までの移動時間であり、$x=A\sin\theta$ とすると

$$
T=\frac4{\omega_0}\int_0^{\pi/2}
\left[1+\frac{\varepsilon}{2}(1+\sin^2\theta)\right]^{-1/2}\dd\theta .
$$

加振応答の正弦近似から独立して、この積分で自由振動の周期を調べられる。$|\varepsilon|\ll1$ では括弧に二項展開を一様に使える。$\int_0^{\pi/2}(1+\sin^2\theta)\dd\theta=3\pi/4$ なので、$T=(2\pi/\omega_0)[1-3\varepsilon/8+O(\varepsilon^2)]$ となる。これを逆数にすると前に求めた振動数移動を再現する。$\alpha>0$ なら $\varepsilon$ とともに被積分関数が全域で小さくなるので、周期は振幅とともに厳密に短くなる。軟化型で $\varepsilon\downarrow-1$ では頂点の近くの速度がさらに小さくなり、周期積分の端で発散する。

:::example[振幅を変えたときの周期]
$\omega_0=10\,\mathrm{s^{-1}}$、$\varepsilon=0.04$ なら、一次近似で周期は $(2\pi/10)(1-0.015)\simeq0.619\,\mathrm s$ である。周期を二桁の測定で比べるなら、無振幅へ外挿した $0.628\,\mathrm s$ との差は約 $9.4\,\mathrm{ms}$ になる。周期の差を比較するときには、角振動数の差をそのまま時間へ換算するのでなく逆数を取る。
:::

### van der Pol 振動子の強非線形極限

無次元の式で $F(u)=u^3/3-u$、$y=F(u)+u'/\mu$ とおく。すると $u'=\mu[y-F(u)]$、$y'=-u/\mu$ である。遅い時刻 $s=\tau/\mu$ で書くと

$$
\frac1{\mu^2}\frac{\dd u}{\dd s}=y-F(u),\qquad
\frac{\dd y}{\dd s}=-u .
$$

$\mu\gg1$ で $y$ を一時的に固定すると、$y=F(u)$ の根の近くの速い変化は $-F'(u)$ に比例する。$|u|>1$ なら $F'=u^2-1>0$ で根へ戻り、$|u|<1$ なら離れる。そのため遅い区間では外側の枝を動き、$u=\pm1$ に達する近くで速い区間へ移る。正の枝で $u\to1$ になると $y\to-2/3$ であり、$F(-2)=-2/3$ だから速い区間では $u\simeq-2$ へ移る。負の枝からは対称に $-1$ の近くから $2$ へ移る。これは極限の軌道の構成であり、分岐点の近くで $y=F(u)$ という近似をそのまま微分することはできない。

正の遅い区間では $F'(u)\dd u/\dd s=-u$ なので、$2$ から $1$ へ移る遅い時刻の増加は $\int_1^2(u-1/u)\dd u=3/2-\log2$ である。もう一方の遅い区間も同じで、速い切り替えにかかる時間を先頭項から省くと、周期は無次元時刻 $\tau$ で $\mu(3-2\log2)$ に近い。弱非線形の平均化で振幅が約 $2$ に近づいた一方、ここではほぼ同じ大きさでも波形と周期が大きく異なる。

:::example[二つの近似が使える範囲]
$\mu=0.05$ では初期の小振幅の成長率は無次元で約 $0.025$ であり、一周期 $2\pi$ で指数は約 $0.157$ だけ変わる。$\mu=20$ で同じ平均化を用いると、一周期で大変化するという予測になり、ゆっくりした振幅の仮定を失う。大きい $\mu$ の周期の先頭項は $20(3-2\log2)\simeq32.3$ であり、外側の枝をたどる時間が周期の主要部分になる。
:::

### Duffing 応答の折れ曲がり

一調波の近似で $x=A\cos(\omega t-\delta)$、$s=A^2$ とする。$\beta=3\alpha/4$、$d=\omega_0^2-\omega^2$、$c=2\gamma\omega$、$f=F_0/m$ とおくと

$$
f^2=g(s)=s[(d+\beta s)^2+c^2],\qquad s\ge0.
$$

硬化型 $\beta>0$ を考える。近似曲線で入力を横軸にすると、垂直に折れ曲がる候補は $g'(s)=3\beta^2s^2+4d\beta s+d^2+c^2=0$ である。その二根は

$$
s_\pm=\frac{-2d\pm\sqrt{d^2-3c^2}}{3\beta}.
$$

従って正の相異なる二点があるのは $d<-\sqrt3\,c$ の場合である。$g$ は $s<s_-$ と $s>s_+$ で増加し、その間は減少する。従って $g(s_+)<f^2<g(s_-)$ なら一調波の振幅候補が三つある。これは応答曲線の近似から数えた候補である。各周期解への摂動が増えるかどうかを判定するには、時間依存の振幅方程式または元の方程式を周期解のまわりで線形化する。

:::example[折れ曲がりの入力を挟む]
無次元化した調和振幅式で $\beta=1$、$d=-2$、$c=1$ とする。このとき $s_-=1$、$s_+=5/3$ であり、$g(1)=2$、$g(5/3)=50/27$ になる。$f^2=19/10$ はこの区間にあるため、三つの正の候補を得る。$f^2$ を $2$ より大きくすると、二つの小振幅側の候補が失われる。
:::

### 遅い複素振幅で折れ曲がりの安定性を判定する

駆動に近い振動を $x=\Re[Z(t)e^{-i\omega t}]$ と置く。$|\omega_0-\omega|/\omega$、$\gamma/\omega$、$|\beta||Z|^2/\omega^2$、$|f|/(\omega^2|Z|)$ が小さく、$Z$ が振動周期より長い時間で変わる範囲を考える。$\ddot Z$ と $2\gamma\dot Z$ を高次として落とし、三次項の基本波 $\beta |Z|^2Z$ を残すと

$$
\dot Z=-\gamma Z-\frac{i}{2\omega}(d+\beta|Z|^2)Z+\frac{if}{2\omega}
$$

となる。固定点の式を絶対値二乗にすると、上の $g(s)=f^2$ が得られる。固定点の位相を座標の回転で除き、そこで $Z_*=\sqrt{s}$ として微小変化を $u+iv$ と書く。$|Z|^2Z$ の一次変化は $3su+isv$ だから、線形化した行列は

$$
\begin{pmatrix}
-\gamma & (d+\beta s)/(2\omega)\\
-(d+3\beta s)/(2\omega) & -\gamma
\end{pmatrix}
$$

である。トレースは $-2\gamma$、行列式は $\gamma^2+(d+\beta s)(d+3\beta s)/(4\omega^2)=g'(s)/(4\omega^2)$ である。従って $\gamma>0$ なら $g'(s)>0$ の固定点は線形化で減衰し、$g'(s)<0$ の固定点には増加方向が一つある。三つの候補がある場合、この近似では小振幅と大振幅の枝が安定で、中間の枝は不安定である。$g'(s)=0$ では固有値が零なので一次の安定性判定は使えない。元の運動方程式の周期解について、近似誤差が無視できない範囲で判定するときは、その解に沿って一周期写像を調べる。

:::example[van der Pol 振動子の周期平均を測る]
無次元方程式 $u''+u=\mu(1-u^2)u'$ に対し、$u$ を一周期で $R\cos\tau$ に近似する。$E=(u^2+u'^2)/2$ の平均変化率は $\mu(R^2/2-R^4/8)$ である。実際、$\langle\sin^2\tau\rangle=1/2$、$\langle\sin^2\tau\cos^2\tau\rangle=1/8$ である。$R=1$ なら $3\mu/8$、$R=3$ なら $-45\mu/8$ となる。この符号の変化が弱い非線形性での定常振幅 $R=2$ に対応する。
:::

## 演習問題

::::exercise{level=basic}
$\omega_0=5\,\mathrm{s^{-1}}$、初期変位 $-3\,\mathrm{mm}$、初速度 $20\,\mathrm{mm\,s^{-1}}$ の単振動の複素振幅と実変位を求めよ。
:::solution
複素振幅の実部を初期変位、虚部を $v_0/\omega_0$ に合わせて $X=(-3+4i)\,\mathrm{mm}$ を得る。よって $x(t)=[-3\cos(5t)+4\sin(5t)]\,\mathrm{mm}$ であり、最大変位は $5\,\mathrm{mm}$ である。
:::
::::

::::exercise{level=basic}
$\gamma=5\,\mathrm{s^{-1}}$、$\omega_0=4\,\mathrm{s^{-1}}$ で自由振動させる。初期変位 $1\,\mathrm{mm}$、初速度零の場合を求めよ。
:::solution
特性根は $-5\pm3$ なので $x=(Ae^{-2t}+Be^{-8t})\,\mathrm{mm}$ と置ける。$A+B=1$、$-2A-8B=0$ を解き $A=4/3$、$B=-1/3$ を得る。原点から十分時間が経つと $e^{-2t}$ の項が支配的になる。
:::
::::

::::exercise{level=basic}
$m=1\,\mathrm{kg}$、$k=9\,\mathrm{N\,m^{-1}}$、$b=2\,\mathrm{kg\,s^{-1}}$ の振動子を $\omega=3\,\mathrm{s^{-1}}$、$F_0=6\,\mathrm{N}$ で駆動する。定常変位と平均電力を求めよ。
:::solution
実際の剛性項と慣性項が相殺するので $X=6/(-6i)=i\,\mathrm m$ となる。従って $x=\sin(3t)\,\mathrm m$、速度振幅は $3\,\mathrm{m\,s^{-1}}$ であり、$b(3)^2/2=9\,\mathrm W$ を散逸する。
:::
::::

::::exercise{level=basic}
吸収の半値周波数が $8\,\mathrm{s^{-1}}$ と $18\,\mathrm{s^{-1}}$ である。$\omega_0$、$\gamma$、$Q$ を求めよ。
:::solution
線幅は $10\,\mathrm{s^{-1}}$ だから $\gamma=5\,\mathrm{s^{-1}}$ である。また $\omega_+\omega_-=\omega_0^2=144\,\mathrm{s^{-2}}$ より $\omega_0=12\,\mathrm{s^{-1}}$、$Q=6/5$ である。
:::
::::

::::exercise{level=advanced}
不足減衰の振動子が平衡位置に静止している。$0$ と $\pi/\omega_d$ に同じ向きの力積を与えて二つ目の直後に静止させたい。必要な比 $J_2/J_1$ を求めよ。
:::solution
$T=\pi/\omega_d$ と置く。第 1 力積で $x(t)=J_1e^{-\gamma t}\sin(\omega_dt)/(m\omega_d)$ となるので、$T-$ で変位は零、速度は $-J_1e^{-\gamma T}/m$ である。第 2 力積は速度を $J_2/m$ だけ増すため $J_2/J_1=e^{-\gamma T}$ とすれば零になる。それ以後の力は零であり、初期値が両方とも零なので静止が続く。
:::
::::

::::exercise{level=advanced}
粘性係数 $b>0$ の強制振動で、非零の外力振幅を一定に保つ。吸収半値において $Z_{\mathrm m}$ の虚部の絶対値が実部に等しいことを示し、低周波側と高周波側の速度位相を求めよ。
:::solution
定常電力は $b|F_0|^2/[2(b^2+X_{\mathrm m}^2)]$、$X_{\mathrm m}=k/\omega-m\omega$ である。最大の半分になるには分母が $2b^2$ であり、$X_{\mathrm m}=\pm b$ が必要十分である。低周波側は $X_{\mathrm m}>0$ なので $V/F_0=1/[b(1+i)]$ の偏角は $-\pi/4$ であり、$e^{-i\omega t}$ の規約では速度が $\pi/4$ 進む。高周波側では偏角が $+\pi/4$ となり、速度が同じ角度だけ遅れる。
:::
::::

::::exercise{level=advanced}
Mathieu 方程式の中心で $\omega_0=80\,\mathrm{s^{-1}}$、$\gamma=0.8\,\mathrm{s^{-1}}$、$h=0.05$ とする。一次近似の成長・減衰率と、対応する複素振幅の位相を求めよ。
:::solution
$\nu=\omega_0$ なので $p=h\omega_0/4=1\,\mathrm{s^{-1}}$ である。実部と虚部が等しい $a$、例えば $\arg a=\pi/4$ の方向は $p-\gamma=0.2\,\mathrm{s^{-1}}$ で成長する。$u=-v$、例えば $\arg a=-\pi/4$ の方向は $-p-\gamma=-1.8\,\mathrm{s^{-1}}$ で減衰する。それぞれの位相には $\pi$ を加えても同じ固有方向を表す。
:::
::::

::::exercise{level=advanced}
Duffing の自由振動で、最大変位 $A>0$ とし、$\omega_0^2+\alpha A^2>0$ を仮定する。最大速度をエネルギーから求め、硬化型での最大速度の一次補正を周波数の一次補正と比べよ。
:::solution
原点に向かう区間でポテンシャルが単調に下がる仮定のもと、最大速度は $x=0$ で生じる。エネルギー保存より $v_{\max}^2=\omega_0^2A^2+\alpha A^4/2$、従って $v_{\max}/A=\omega_0+\alpha A^2/(4\omega_0)+O(\varepsilon^2\omega_0)$ である。角振動数の一次補正は $3\alpha A^2/(8\omega_0)$ であり、係数が異なる。非線形の軌道は正弦波だけでは表せないためである。
:::
::::

::::exercise{level=advanced}
平均化された van der Pol 振動子を $R(0)=1$ から出発させる。$R=\sqrt2$ になる無次元時刻を求めよ。
:::solution
$s=4/(1+3e^{-\mu\tau})$ に $s=2$ を代入すると $3e^{-\mu\tau}=1$ である。従って $\tau=(\log3)/\mu$ であり、実時間は $t=\tau/\omega_0$ である。
:::
::::

::::exercise{level=advanced}
調和入力の下で $\omega_0=5\,\mathrm{s^{-1}}$ とする。変位振幅の二乗が零周波数でちょうどピークの半分になる $\gamma>0$ を求め、吸収の半値全幅と比べよ。
:::solution
ピークが正周波数にあるとし、$D(0)=\omega_0^4$ と $D_{\min}=4\gamma^2(\omega_0^2-\gamma^2)$ を使う。半値条件は $D(0)=2D_{\min}$ なので、$r=\gamma^2/\omega_0^2$ は $8r(1-r)=1$ を満たす。$\omega_r^2>0$ の条件 $r<1/2$ から $\gamma=5\sqrt{(1-1/\sqrt2)/2}\,\mathrm{s^{-1}}$ を選ぶ。このとき変位二乗の下の半値は零である。一方、吸収の半値全幅 $2\gamma$ は二つの正周波数の差であり、吸収の二つの半値周波数はいずれも正である。
:::
::::

## 確認テスト

::::quiz
$e^{-i\omega t}$ の規約で、純粋に正の虚数の変位応答 $X/F_0$ はどの位相を表すか。

- [ ] 変位は力に対して $\pi/2$ 進む。
- [x] 変位は力に対して $\pi/2$ 遅れる。
- [ ] 変位と力は同位相である。

:::explanation
正の実数 $B$ に対し $\Re(iBe^{-i\omega t})=B\sin\omega t=B\cos(\omega t-\pi/2)$ である。力を $\cos\omega t$ とすれば、変位が最大となるのは力の最大より $\pi/(2\omega)$ 後である。
:::
::::

::::quiz{answer="20" tol="1e-6"}
吸収ピークが $40\,\mathrm{Hz}$、半値全幅が $2\,\mathrm{Hz}$ のとき $Q$ はいくつか。
:::explanation
$\omega_0=2\pi(40)\,\mathrm{s^{-1}}$、$\Delta\omega=2\pi(2)\,\mathrm{s^{-1}}$ より、比 $Q=\omega_0/\Delta\omega=20$ を得る。
:::
::::

::::quiz
$\gamma>0$ の振動子の遅延応答の性質として正しいものをすべて選べ。

- [x] 力積を加える前の応答は零である。
- [ ] その周波数応答は上半平面に極をもつ。
- [x] 実周波数軸上で $\chi''$ は奇関数である。

:::explanation
$G(t)=0$ が $t<0$ で成り立つ。$\chi(z)$ の積分核には $e^{izt}$ があり、上半平面で収束して正則になる。さらに $G$ が実数値なので $\chi(-\omega)=\conj{\chi(\omega)}$ であり、虚部は奇である。
:::
::::

::::quiz{answer="0.08" tol="1e-6"}
$\omega_0=50\,\mathrm{s^{-1}}$、$\gamma=1\,\mathrm{s^{-1}}$ で、主パラメトリック共振の中心における $|h|$ の一次近似の閾値を求めよ。
:::explanation
中心では $|\Delta|=0$ であり、$|p|=|h|\omega_0/4$ が $\gamma$ を超えると成長する。境界は $|h|=4\gamma/\omega_0=0.08$ である。
:::
::::

::::quiz
弱い非線形性の効果について正しいものを選べ。

- [x] 正の三次剛性により Duffing 振動子の固有振動数は振幅とともに増える。
- [ ] $\mu>0$ の van der Pol 振動子は原点近くで全時刻にわたってエネルギーを失う。
- [ ] 無減衰の Duffing 振動子では、復元力が非線形であるだけで毎周期エネルギーが失われる。

:::explanation
周波数の一次補正は $3\alpha A^2/(8\omega_0)$ である。これに対して van der Pol 振動子のエネルギー変化は $\mu(1-u^2)u'^2$ で、原点付近では非負である。Duffing 振動子の復元力はポテンシャルの微分なので、外力と摩擦がなければエネルギーを保存する。
:::
::::

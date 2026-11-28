---
title: "Schwarzschild 解"
description: "真空球対称計量を求め、保存量で粒子と光の軌道を解析する。"
notation:
  - symbol: "m"
    meaning: "質量 M に対応する長さ GM"
  - symbol: "f(r)"
    meaning: "1−2m/r"
  - symbol: "E,L"
    meaning: "単位質量当たりまたは光のアフィン接線の保存量"
  - symbol: "V_\\kappa"
    meaning: "動径有効ポテンシャル"
---

## この章の目標

- 球対称な計量を面積半径で表す。
- 真空方程式を積分し、Birkhoff の定理の局所的主張を示す。
- 質量をもつ粒子の円軌道とその安定性を求める。
- 光子球と捕獲の臨界衝突径数を求める。

## 静的球対称計量の一般形

球対称性とは、空間回転群が等長変換として作用し、一般の軌道が二次元球面になることである。球面上の一点を固定する回転が接平面の全回転を含むので、そこから球面に垂直な方向への混合成分は零である。また球面内の対称形式は標準球面計量の倍となる。同じ球面上で係数が一定なので、面積を $4\pi r^2$ と書くと $r$ が座標不変に定まる。これを :term[面積半径]{en="areal radius" yomi="めんせきはんけい"} という。

静的で $\nabla r$ が空間的な領域では直交面を用い、時間座標に依存しない

$$
\dd s^2=-e^{2\Phi(r)}\dd t^2+e^{2\lambda(r)}\dd r^2+r^2(\dd\theta^2+\sin^2\theta\,\dd\phi^2)
$$

を取れる。半径の差は一般に固有距離でなく、同じ時刻での尺は $e^\lambda\dd r$ である。例えば $e^{2\lambda}=4$ なら、半径が $0.1\,\mathrm m$ だけ異なる近接球面の固有距離は約 $0.2\,\mathrm m$ である。

### 面積半径が局所座標になる条件

球対称な領域では回転群の二次元軌道を一枚の球として扱う。各球の誘導計量は回転で不変だから、その一点の接平面では等方性群の回転で不変な二次形式である。従って単位球の計量に、球ごとに一定な正の係数を掛けたものになる。面積が $4\pi r^2$ となるよう $r>0$ を選べば、その係数は $r^2$ である。また球の接平面に不変な非零の一次形式はないため、軌道商の座標と球面座標の混合成分は消える。$dr\ne0$ の点では逆関数定理で $r$ を商の座標にできる。ただし $g^{-1}(\dd r,\dd r)>0$ が必要になるのは、これと直交する座標を時間方向に選ぶ段階である。

:::example[固有半径と面積半径]
空間計量を $\dd\ell^2+L^2\sin^2(\ell/L)\dd\Omega^2$ とする。球面積は $4\pi L^2\sin^2(\ell/L)$ なので $r=L\sin(\ell/L)$ である。$\ell=\pi L/2$ では $dr/d\ell=0$ となり、$\ell$ を $r$ に置き換えると係数 $\dd\ell^2=\dd r^2/(1-r^2/L^2)$ が発散する。元の $\ell$ 座標はここで滑らかである。さらに $\ell$ と $\pi L-\ell$ では同じ $r$ になるので、一つの $r$ が空間全体を識別できない。
:::

### 光的な切片に適合させる

球対称で $dr\ne0$ の領域で、光的な時刻関数 $v$ を選ぶ。$v$ 一定の曲線の接線 $\partial_r$ は光的になるので、軌道商の計量を
$$
q=-e^{2\psi(v,r)}f(v,r)\dd v^2+2e^{\psi(v,r)}\dd v\dd r
$$
の形にできる。交差項は Lorentz 非退化性から非零で、ここでは向きを選んで正とした。逆計量は $q^{vv}=0,q^{vr}=e^{-\psi},q^{rr}=f$ なので $f=(\nabla r)^2$ である。$\Gamma^r{}_{rr}=\partial_r\psi$ だから、真空で商の曲率の $rr$ 成分を使うと $R_{rr}=(2/r)\partial_r\psi=0$ である。商の二次元 Ricci 自体のこの成分は $q_{rr}=0$ により消える。従って $\psi$ は $v$ の関数であり、$\dd\bar v=e^{\psi(v)}\dd v$ により零に戻せる。この座標は $f=0$ の点でも非退化なので、真空解を地平線の外だけで扱う必要がない。

### 静的 Killing 場と回転の相性

光的な切片を使う座標は、軌道商上の二本の光的方向のうち、$\dd r$ が消えない方を選ぶと構成できる。その方向の積分曲線を局所に引き、曲線ごとに一定となるラベルを $v$ に取れば、$(v,r)$ の微分は独立になる。静的な領域で対角座標を取る場合は、時間的 Killing 場が球対称性に適合し、球の面積半径を保存する場を用いる。球面に沿う成分を持たないその場は商の一方向を張る。$\nabla r$ が空間的なら直交する時間方向を選べるため、混合成分を消しても $r$ を座標として保てる。

:::example[交差成分がある静的な球対称座標]
二次元部分を $-A(r)\dd s^2+2B(r)\dd s\dd r+C(r)\dd r^2$ とする。$A>0$ かつ $AC+B^2>0$ なら Lorentz 署名で、$\partial_s$ は時間的である。$t=s-\int B/A\,\dd r$ と置くと
$$
q=-A\dd t^2+(C+B^2/A)\dd r^2.
$$
従って $\dd r$ のノルムは $A/(AC+B^2)$ になる。例えば $A=1$、$B=1$、$C=0$ なら $t=s-r$ で Minkowski の球対称表示を得る。一方 $A\to0$ で $B$ が零でない場合、元の二次元計量の行列式は非零のままだが、対角化に使った積分は発散し得る。光的な切片に適合させる座標はこの割り算を必要としない。
:::

### 面積半径と等方座標での半径を使い分ける

空間の球対称な正定値計量 $A(r)^2\dd r^2+r^2\dd\Omega^2$ では、$\dd\varrho/\varrho=A(r)\dd r/r$ を積分すると、新しい半径 $\varrho$ で空間計量を $(r/\varrho)^2(\dd\varrho^2+\varrho^2\dd\Omega^2)$ にできる。$A>0,r>0$ の範囲では $\varrho$ を単調に取れる。これは空間の各方向を同じ倍率で伸縮する形であり、$\varrho$ で計算した座標球の面積 $4\pi\varrho^2$ はそのまま物理的な面積にはならない。例えば平坦な場合は $A=1$ なので $\varrho$ と $r$ は定数倍であり、無限遠で両者の尺度をそろえれば一致する。湾曲した静的切片では以下の例のように二つの半径の差が重力の次数で現れる。

## Schwarzschild 解の導出

### 真空方程式の積分

$f=e^{-2\lambda}$ と置く。:eqref[eq-gr-christoffel] より $\Gamma^t{}_{tr}=\Phi'$、$\Gamma^r{}_{tt}=f e^{2\Phi}\Phi'$、$\Gamma^r{}_{rr}=-f'/2f$、$\Gamma^r{}_{\theta\theta}=-rf$、$\Gamma^\theta{}_{r\theta}=1/r$ であり、$\phi$ 成分には球面の係数が入る。Ricci の角成分は $R_{\theta\theta}=1-f-rf\Phi'-rf'/2$ である。時間・動径の縮約には

$$
R^t{}_t=-f(\Phi''+\Phi'^2+2\Phi'/r)-f'\Phi'/2,\qquad
R^r{}_r=-f(\Phi''+\Phi'^2)-f'\Phi'/2-f'/r
$$

を得る。従って Einstein テンソルは

$$
G^t{}_t=\frac{f-1+rf'}{r^2},\qquad
G^r{}_r=\frac{f-1+2rf\Phi'}{r^2}.
\label{eq-gr-spherical-field}
$$

本節では $\Lambda=0$ かつ $T_{ab}=0$ とする。$G^t{}_t=0$ から $(rf)'=1$、従って $f=1+C/r$ である。$G^r{}_r=0$ との比較から $2f\Phi'=f'$ なので $e^{2\Phi}=A f$、定数 $A>0$ を時刻の尺度に吸収する。角方程式は $G^\theta{}_\theta=f''/2+f'/r=0$ で、自動的に成立する。遠方で $g_{tt}=-(1+2\Phi_{\rm N})$、$\Phi_{\rm N}=-GM/r$ と照合すると $C=-2GM$ である。よって :term[Schwarzschild 計量]{en="Schwarzschild metric" yomi="しゅわるつしるとけいりょう"} は

$$
\dd s^2=-f\dd t^2+f^{-1}\dd r^2+r^2\dd\Omega^2,\qquad f=1-\frac{2m}{r},\quad m=GM.
\label{eq-gr-schwarzschild}
$$

SI では $m=GM/c^2$ で、太陽質量に対して約 $1.48\,\mathrm{km}$ である。

:::example[一関数の球対称計量で真空方程式を縮約する]
$f>0,r>0$ の範囲で $g_{tt}=-f,g_{rr}=1/f$ とする。真空方程式の一部を使って時間の規格化を選んだ後、残る関数を決める手順である。 座標の順序を $(t,r,\theta,\phi)$ と固定する。成分行列と逆行列は次の組になる。

$$
\begin{aligned}
(g_{ab})&=\left[\begin{matrix}- f{\left(r \right)} & 0 & 0 & 0\\0 & \frac{1}{f{\left(r \right)}} & 0 & 0\\0 & 0 & r^{2} & 0\\0 & 0 & 0 & r^{2} \sin^{2}{\left(\theta \right)}\end{matrix}\right]\\
(g^{ab})&=\left[\begin{matrix}- \frac{1}{f{\left(r \right)}} & 0 & 0 & 0\\0 & f{\left(r \right)} & 0 & 0\\0 & 0 & \frac{1}{r^{2}} & 0\\0 & 0 & 0 & \frac{1}{r^{2} \sin^{2}{\left(\theta \right)}}\end{matrix}\right]\\
\det g&=- r^{4} \sin^{2}{\left(\theta \right)}
\end{aligned}
$$

逆行列を掛けて単位行列を得ることを確かめてから、:eqref[eq-gr-christoffel] を使う。同じ下添字の順序を二度計算せず、$b\le c$ の組を選ぶ。中央の列には、座標に依存する計量成分の微分だけを残した。最後の列は微分を実行した値であり、表に出ない成分は添字交換で得られるものを除いて零である。

| 成分 | 微分への代入 | 結果 |
|---|---|---|
| $\Gamma^{t}_{t r}$ | $- \frac{1}{2 f{\left(r \right)}}\left(\partial_{r}\left(- f{\left(r \right)}\right)\right)$ | $\frac{f'{\left(r \right)}}{2 f{\left(r \right)}}$ |
| $\Gamma^{r}_{t t}$ | $\frac{f{\left(r \right)}}{2}\left(-\partial_{r}\left(- f{\left(r \right)}\right)\right)$ | $\frac{f{\left(r \right)} f'{\left(r \right)}}{2}$ |
| $\Gamma^{r}_{r r}$ | $\frac{f{\left(r \right)}}{2}\left(\partial_{r}\left(\frac{1}{f{\left(r \right)}}\right)+\partial_{r}\left(\frac{1}{f{\left(r \right)}}\right)-\partial_{r}\left(\frac{1}{f{\left(r \right)}}\right)\right)$ | $- \frac{f'{\left(r \right)}}{2 f{\left(r \right)}}$ |
| $\Gamma^{r}_{\theta \theta}$ | $\frac{f{\left(r \right)}}{2}\left(-\partial_{r}\left(r^{2}\right)\right)$ | $- r f{\left(r \right)}$ |
| $\Gamma^{r}_{\phi \phi}$ | $\frac{f{\left(r \right)}}{2}\left(-\partial_{r}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $- r f{\left(r \right)} \sin^{2}{\left(\theta \right)}$ |
| $\Gamma^{\theta}_{r \theta}$ | $\frac{1}{2 r^{2}}\left(\partial_{r}\left(r^{2}\right)\right)$ | $\frac{1}{r}$ |
| $\Gamma^{\theta}_{\phi \phi}$ | $\frac{1}{2 r^{2}}\left(-\partial_{\theta}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $- \frac{\sin{\left(2 \theta \right)}}{2}$ |
| $\Gamma^{\phi}_{r \phi}$ | $\frac{1}{2 r^{2} \sin^{2}{\left(\theta \right)}}\left(\partial_{r}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{1}{r}$ |
| $\Gamma^{\phi}_{\theta \phi}$ | $\frac{1}{2 r^{2} \sin^{2}{\left(\theta \right)}}\left(\partial_{\theta}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{1}{\tan{\left(\theta \right)}}$ |

曲率の計算では :eqref[eq-gr-curvature-components] の微分差と積差を分ける。第二列は $\partial_c\Gamma^a{}_{db}$、第三列は $-\partial_d\Gamma^a{}_{cb}$ であり、第四列は接続の積を符号付きで足したものである。三列を足して最後の列を得る。

| 成分 | 第一の微分 | 第二の微分 | 接続の積 | 結果 |
|---|---|---|---|---|
| $R^{t}_{r t r}$ | $0$ | $- \frac{f''{\left(r \right)}}{2 f{\left(r \right)}} + \frac{\left(f'{\left(r \right)}\right)^{2}}{2 f^{2}{\left(r \right)}}$ | $- \frac{\left(f'{\left(r \right)}\right)^{2}}{4 f^{2}{\left(r \right)}} - \frac{\left(f'{\left(r \right)}\right)^{2}}{4 f^{2}{\left(r \right)}}$ | $- \frac{f''{\left(r \right)}}{2 f{\left(r \right)}}$ |
| $R^{t}_{\theta t \theta}$ | $0$ | $0$ | $- \frac{r f'{\left(r \right)}}{2}$ | $- \frac{r f'{\left(r \right)}}{2}$ |
| $R^{t}_{\phi t \phi}$ | $0$ | $0$ | $- \frac{r \sin^{2}{\left(\theta \right)} f'{\left(r \right)}}{2}$ | $- \frac{r \sin^{2}{\left(\theta \right)} f'{\left(r \right)}}{2}$ |
| $R^{r}_{\theta r \theta}$ | $- r f'{\left(r \right)} - f{\left(r \right)}$ | $0$ | $\frac{r f'{\left(r \right)}}{2} + f{\left(r \right)}$ | $- \frac{r f'{\left(r \right)}}{2}$ |
| $R^{r}_{\phi r \phi}$ | $- r \sin^{2}{\left(\theta \right)} f'{\left(r \right)} - f{\left(r \right)} \sin^{2}{\left(\theta \right)}$ | $0$ | $\frac{r \sin^{2}{\left(\theta \right)} f'{\left(r \right)}}{2} + f{\left(r \right)} \sin^{2}{\left(\theta \right)}$ | $- \frac{r \sin^{2}{\left(\theta \right)} f'{\left(r \right)}}{2}$ |
| $R^{\theta}_{\phi \theta \phi}$ | $- \cos{\left(2 \theta \right)}$ | $0$ | $- f{\left(r \right)} \sin^{2}{\left(\theta \right)} + \frac{\sin{\left(2 \theta \right)}}{2 \tan{\left(\theta \right)}}$ | $- \left(f{\left(r \right)} - 1\right) \sin^{2}{\left(\theta \right)}$ |

第一添字を計量で下ろすと $R_{abcd}=g_{ak}R^k{}_{bcd}$ となる。対の反対称性と交換対称性を用いて、非零の独立な組を下表に示す。

| 成分 | 結果 |
|---|---|
| $R_{t r t r}$ | $\frac{f''{\left(r \right)}}{2}$ |
| $R_{t \theta t \theta}$ | $\frac{r f{\left(r \right)} f'{\left(r \right)}}{2}$ |
| $R_{t \phi t \phi}$ | $\frac{r f{\left(r \right)} \sin^{2}{\left(\theta \right)} f'{\left(r \right)}}{2}$ |
| $R_{r \theta r \theta}$ | $- \frac{r f'{\left(r \right)}}{2 f{\left(r \right)}}$ |
| $R_{r \phi r \phi}$ | $- \frac{r \sin^{2}{\left(\theta \right)} f'{\left(r \right)}}{2 f{\left(r \right)}}$ |
| $R_{\theta \phi \theta \phi}$ | $- r^{2} \left(f{\left(r \right)} - 1\right) \sin^{2}{\left(\theta \right)}$ |

Ricci テンソルは第一・第三添字を縮約して得る。曲率が残る場合も相殺が起こるので、各項を落とさずに加える。

| 成分 | 縮約の項 | 代入 | 結果 |
|---|---|---|---|
| $R_{t t}$ | $R^{r}_{t r t}+R^{\theta}_{t \theta t}+R^{\phi}_{t \phi t}$ | $\frac{f{\left(r \right)} f''{\left(r \right)}}{2} + \frac{f{\left(r \right)} f'{\left(r \right)}}{2 r} + \frac{f{\left(r \right)} f'{\left(r \right)}}{2 r}$ | $\frac{\left(r f''{\left(r \right)} + 2 f'{\left(r \right)}\right) f{\left(r \right)}}{2 r}$ |
| $R_{r r}$ | $R^{t}_{r t r}+R^{\theta}_{r \theta r}+R^{\phi}_{r \phi r}$ | $- \frac{f''{\left(r \right)}}{2 f{\left(r \right)}} - \frac{f'{\left(r \right)}}{2 r f{\left(r \right)}} - \frac{f'{\left(r \right)}}{2 r f{\left(r \right)}}$ | $- \frac{r f''{\left(r \right)} + 2 f'{\left(r \right)}}{2 r f{\left(r \right)}}$ |
| $R_{\theta \theta}$ | $R^{t}_{\theta t \theta}+R^{r}_{\theta r \theta}+R^{\phi}_{\theta \phi \theta}$ | $- \frac{r f'{\left(r \right)}}{2} - \frac{r f'{\left(r \right)}}{2} - \left(f{\left(r \right)} - 1\right)$ | $- r f'{\left(r \right)} - f{\left(r \right)} + 1$ |
| $R_{\phi \phi}$ | $R^{t}_{\phi t \phi}+R^{r}_{\phi r \phi}+R^{\theta}_{\phi \theta \phi}$ | $- \frac{r \sin^{2}{\left(\theta \right)} f'{\left(r \right)}}{2} - \frac{r \sin^{2}{\left(\theta \right)} f'{\left(r \right)}}{2} - \left(f{\left(r \right)} - 1\right) \sin^{2}{\left(\theta \right)}$ | $- \left(r f'{\left(r \right)} + f{\left(r \right)} - 1\right) \sin^{2}{\left(\theta \right)}$ |

逆計量と縮約すると Ricci スカラーとなる。さらに $G_{ab}=R_{ab}-g_{ab}R/2$ を成分ごとに計算する。

$$
\begin{aligned}
R&=\sum_{a,b}g^{ab}R_{ab}\\
&=- \frac{r f''{\left(r \right)} + 2 f'{\left(r \right)}}{2 r} - \frac{r f''{\left(r \right)} + 2 f'{\left(r \right)}}{2 r} + \frac{- r f'{\left(r \right)} - f{\left(r \right)} + 1}{r^{2}} - \frac{r f'{\left(r \right)} + f{\left(r \right)} - 1}{r^{2}}\\
&=- \frac{r^{2} f''{\left(r \right)} + 4 r f'{\left(r \right)} + 2 f{\left(r \right)} - 2}{r^{2}}
\end{aligned}
$$

非零の Einstein テンソル成分は次のとおりである。

| 成分 | Ricci 項 | 跡の項 | 差 |
|---|---|---|---|
| $G_{t t}$ | $\frac{\left(r f''{\left(r \right)} + 2 f'{\left(r \right)}\right) f{\left(r \right)}}{2 r}$ | $\frac{\left(r^{2} f''{\left(r \right)} + 4 r f'{\left(r \right)} + 2 f{\left(r \right)} - 2\right) f{\left(r \right)}}{2 r^{2}}$ | $- \frac{\left(r f'{\left(r \right)} + f{\left(r \right)} - 1\right) f{\left(r \right)}}{r^{2}}$ |
| $G_{r r}$ | $- \frac{r f''{\left(r \right)} + 2 f'{\left(r \right)}}{2 r f{\left(r \right)}}$ | $- \frac{r^{2} f''{\left(r \right)} + 4 r f'{\left(r \right)} + 2 f{\left(r \right)} - 2}{2 r^{2} f{\left(r \right)}}$ | $\frac{r f'{\left(r \right)} + f{\left(r \right)} - 1}{r^{2} f{\left(r \right)}}$ |
| $G_{\theta \theta}$ | $- r f'{\left(r \right)} - f{\left(r \right)} + 1$ | $- \frac{r^{2} f''{\left(r \right)}}{2} - 2 r f'{\left(r \right)} - f{\left(r \right)} + 1$ | $\frac{r \left(r f''{\left(r \right)} + 2 f'{\left(r \right)}\right)}{2}$ |
| $G_{\phi \phi}$ | $- \left(r f'{\left(r \right)} + f{\left(r \right)} - 1\right) \sin^{2}{\left(\theta \right)}$ | $- \frac{\left(r^{2} f''{\left(r \right)} + 4 r f'{\left(r \right)} + 2 f{\left(r \right)} - 2\right) \sin^{2}{\left(\theta \right)}}{2}$ | $\frac{r \left(r f''{\left(r \right)} + 2 f'{\left(r \right)}\right) \sin^{2}{\left(\theta \right)}}{2}$ |

$G^t{}_t=(rf'+f-1)/r^2$ だから、零にする式は $(rf)'=1$ である。従って $f=1+C/r$ となる。角成分は $f''/2+f'/r=0$ であり、この解で実際に消える。無限遠で時間を規格化し、弱い場で $\Phi\simeq C/(2r)$ を点質量のポテンシャルと比べれば $C=-2GM$ である。
:::

### 弱い場の規格化と質量の測り方

球対称真空の積分定数 $m$ は長さである。遠方で時刻の正規化を $g_{tt}\to-1$ として固定すると、静止時計について $\sqrt f=1-m/r+\order{(m/r)^2}$ となる。従って二つの遠い球 $r_1,r_2$ に置いた静止時計の比較では $\nu_2/\nu_1=1-m(1/r_1-1/r_2)+\order{(m/r_1)^2+(m/r_2)^2}$ となる。他方、時間的な円軌道の角速度からも $m=r^3\Omega^2$ を得る。例えば $m/r=10^{-6}$ の球を共通に用い、半径 $r$ と $2r$ の静止時計を比較すれば周波数比のずれは約 $5\times10^{-7}$ である。このように $m$ の物理的な同定には時間の尺度をどこで固定したかが重要である。時刻を定数倍すると同じ世界線でも $\Omega$ の値が変わるため、無限遠での時計の規格化が必要になる。

### 正規直交枠の曲率を不変量へまとめる

$f>0$ の外部で $\theta^{\hat0}=\sqrt f\,\dd t$、$\theta^{\hat1}=\dd r/\sqrt f$、$\theta^{\hat2}=r\dd\theta$、$\theta^{\hat3}=r\sin\theta\,\dd\phi$ とする。先の座標成分に枠の係数を掛けると
$$
\begin{aligned}
R_{\hat0\hat1\hat0\hat1}&=f''/2,&
R_{\hat0\hat I\hat0\hat I}&=f'/(2r),\\
R_{\hat1\hat I\hat1\hat I}&=-f'/(2r),&
R_{\hat2\hat3\hat2\hat3}&=(1-f)/r^2,\qquad I=2,3 .
\end{aligned}
$$
時間添字は各非零成分に零個または二個ずつ現れるため、曲率二乗の縮約では各二乗に正号が付く。各成分の二つの添字対の順序を四通り数えると
$$
R_{abcd}R^{abcd}=(f'')^2+4(f'/r)^2+4[(1-f)/r^2]^2 .
$$
$f=1-2m/r$ では厳密に $48m^2/r^6$ である。$m\ne0$ なら $r\to0$ でこの値は発散するが、$m>0$ の $r\to2m$ では $3/(4m^4)$ に近づく。後者の有限性だけから正則性を証明することはできない。先に導入した光的な切片の座標に計量を移して、行列式が非零のまま滑らかに延長できることを併せて確認する。

:::example[時計の加速度と潮汐の違いを測る]
静止枠の時間軸に対する電気型曲率は $\diag(-2,1,1)m/r^3$ である。ある事象で静止枠と同じ速度を持つように粒子を自由に放せば、その事象での隣接測地線の相対加速度はこの行列の負になる。同じ事象で支持された静止観測者の加速度計は $m/(r^2\sqrt{f})$ を読む。例えば $r=4m$ では潮汐の動径係数は $1/(32m^2)$、静止に要する固有加速度は $\sqrt2/(16m)$ である。相対加速度の係数は長さの逆二乗であり、基線の長さを掛けて初めて加速度になる。静止加速度とは異なる測定を行っている。
:::

### 半径と時計の規格化を別々に読む

外部の静止切片で $r=r_1$ から $r_2$ までの動径固有距離は $\int_{r_1}^{r_2}\dd r/\sqrt f$ である。$m>0$ で $r=2m\cosh^2u$ と置くと $\sqrt f=\tanh u$、$\dd r=4m\cosh u\sinh u\,\dd u$ である。従って下限 $2m$ からの積分は $2mu+m\sinh(2u)$ で記述できる。元の半径で書けば
$$
\ell(r)=\sqrt{r(r-2m)}+
2m\log\!\left(\frac{\sqrt r+\sqrt{r-2m}}{\sqrt{2m}}\right).
$$
$r=2m+\delta$、$0<\delta\ll m$ なら $\ell\simeq2\sqrt{2m\delta}$ であり、$\sqrt f\simeq \ell/(4m)$ となる。時刻の率が零へ近づいても、静止切片上の固有距離の極限は有限である。この計算は外部の静止切片を選んだ結果であり、移動する観測者の固有時はその世界線に沿って別に積分する。

### Schwarzschild 外部を等方半径で表す

外部の計量では $A=(1-2m/r)^{-1/2}$ である。$r=2m\cosh^2u$ と置くと $A\dd r/r=2\dd u$ だから $\varrho=(m/2)e^{2u}$ と選べる。この定数は $\varrho/r\to1$ で固定した。$\cosh^2u$ を指数で書き直すと
$$
r=\varrho\left(1+\frac m{2\varrho}\right)^2,\qquad
\dd s^2=-\left(\frac{1-m/(2\varrho)}{1+m/(2\varrho)}\right)^2\dd t^2
+\left(1+\frac m{2\varrho}\right)^4(\dd\varrho^2+\varrho^2\dd\Omega^2).
$$
$r$ の $\varrho$ 微分は $1-m^2/(4\varrho^2)$ で、外部を表すには $\varrho>m/2$ とする。$\varrho=m/2$ で時間係数が零になり、面積半径も座標として独立でなくなる。この形を $0<\varrho<m/2$ に代入しても $r>2m$ になり、同じ外部の面積半径を再び記録する。内部 $r<2m$ をこの式の値域から得ることはできない。

弱い場 $\varrho\gg m$ では時間成分は $-(1-2m/\varrho)$、空間の共形因子は $1+2m/\varrho$ が一次である。従って Newton のポテンシャルが $-m/\varrho$ であることと同時に、空間の一次の補正も読み取れる。例えば $\varrho=2m$ では $r=25m/8$、静止時計の率は $3/5$ である。これを面積半径の表示に戻しても $1-2m/r=9/25$ となり、時計の測定値は一致する。変換をするとき、半径の数値をそのまま両方に代入して比較しないようにする。

## Birkhoff の定理

### 時間依存性の制限

球対称真空で $\nabla r$ が非零かつ非光的なら、二次元の軌道商でそれに直交する座標を選べる。$f=(\nabla r)^2$、$\eps=\sgn(f)$ とし、時刻ラベルを $s$ として商の計量を $-\eps e^{2\Phi(s,r)}\dd s^2+f(s,r)^{-1}\dd r^2$ と書く。$f<0$ の領域ではこの $s$ は空間的であり、時間的な方向は $r$ の変化に沿う。二次元計量を $q$ とすると、球面方向の接続は $\Gamma^i{}_{Aj}=r^{-1}\partial_A r\delta^i_j$、$\Gamma^A{}_{ij}=-r\nabla^A r\,\gamma_{ij}$ である。曲率を縮約すると $R_{AB}[g]=R_{AB}[q]-2r^{-1}\nabla_A\nabla_B r$、$R_{ij}[g]=(1-f-r\Box_q r)\gamma_{ij}$ が得られる。前者の追加項は、球面添字を縮約した $\partial_B(2\partial_A r/r)$ と接続の積を合わせたものである。二次元の Ricci は $q$ に比例するので、$R_{sr}=2\Gamma^r{}_{sr}/r=-\partial_s f/(rf)$ となる。真空で $\partial_s f=0$ であり、残る対角方程式は :eqref[eq-gr-spherical-field] と同じ形になる。したがって $f=1-2m/r$、$e^{2\Phi}=A(s)\abs f$ となり、$\dd t=\sqrt{A(s)}\dd s$ で時間依存を消せる。

これが :term[Birkhoff の定理]{en="Birkhoff theorem" yomi="びるくほふのていり"} の非光的な面積半径を持つ領域での主張である。外部で $f>0$ なら局所的に静的であり、内部で $f<0$ なら球対称性から得た同じ Killing 場は空間的になる。$\nabla r$ が光的な面の扱いには正則な座標が必要である。真空の開領域で $f$ が恒等的に零、$\dd r\ne0$ なら $\nabla_{\nabla r}\nabla r=0$ となり、二次元でこの零ベクトルの方向を組に取れば $\Box_qr=0$ となる。すると $R_{ij}=\gamma_{ij}$ で矛盾する。$r$ が定数の開領域でも同じ矛盾が得られる。従って真空の球対称な脈動は、外部の質量が一定である限り外部に時変な球対称重力波を作らない。

### 動径勾配が光的になる点

二次元の軌道商の Ricci は $q_{AB}$ に比例するので、真空の $AB$ 方程式の跡なし部分から $\nabla_A\nabla_Br=(\Box_qr)q_{AB}/2$ を得る。角方程式は $\Box_qr=(1-f)/r$ である。従って $\nabla_A f=2\nabla^B r\,\nabla_A\nabla_B r=(1-f)\partial_A r/r$ となり、$\nabla_A[r(1-f)/2]=0$ である。これにより商の全域で $m=r(1-f)/2$ が局所的な定数となる。この計算は $f$ で割らないので $f=0$ でも成り立つ。

二次元の計量体積形式 $\epsilon_{AB}$ を用い、$K_A=\epsilon_{AB}\nabla^B r$ とする。体積形式は接続で平行であり、上で得た Hessian は $q$ の倍なので $\nabla_{(A}K_{B)}=0$ である。$K(r)=0$ なので球面の $r^2\dd\Omega^2$ に対しても Killing 場になる。また $q(K,K)=-f$ である。$K\ne0$ の零面では $r$ が座標だから、$K=\partial_v$ とする座標で商計量を $-f\dd v^2+2A\dd v\dd r+B\dd r^2$ と書ける。$K_A=\epsilon_{AB}\nabla^B r$ によりこのブロックの行列式は $-1$ であり、$A^2+fB=1$ である。$A\simeq1$ の向きを選び $\dd v=\dd v'+h(r)\dd r$、$h=-B/(1+A)$ とすると $B+2Ah-fh^2=0$、新しい混合係数は $A-fh=1$ となる。こうして零面上でも $-f\dd v'^2+2\dd v'\dd r$ の正則な形を得る。$\dd r=0$ の二分岐球面の近くでは、二つの零座標を用いる構成を [Schwarzschild ブラックホール](/physics/general-relativity/black-holes/) で扱う。

### 商の面積形式から時間対称性を作る

商の二次元 Lorentz 計量を $q$、面積形式を $\epsilon_{AB}$ とする。真空方程式の $AB$ 成分の跡を除くと、$\nabla_A\nabla_B r=(\Box_q r)q_{AB}/2$ である。ここで $\nabla$ は商の接続である。$\epsilon$ は計量接続で平行なので、$K_A=\epsilon_A{}^B\nabla_B r$ に対し $\nabla_AK_B$ は $\epsilon_{BA}$ の倍になり、対称化すると零になる。さらに $K(r)=0$ なので、この場で球面計量 $r^2\gamma$ も変わらない。従って $K$ は四次元の Killing 場である。$K^2=-(\nabla r)^2$ は Lorentz 直交枠で $\epsilon$ が時間・空間成分を交換することから従う。零地平線ではこの場は光的になるが、構成には $1/f$ を使っていない。

:::example[正則な零座標で Killing 場を確認する]
$q=-f(r)\dd v^2+2\dd v\dd r$ では $\det q=-1$ である。$\epsilon_{rv}=1$ を取ると $K^A=q^{AC}\epsilon_C{}^r$ の成分は $(1,0)$ であり、$K=\partial_v$ を得る。計量の係数に $v$ がなく、ノルムは $-f$ である。従って $f=0$ で時間的でなくなっても Killing 方程式自体には問題がない。
:::

### 質量関数を真空方程式で保存させる

軌道商で $F=q^{AB}\partial_A r\partial_B r$、$m(r)=r(1-F)/2$ とする。真空で $\nabla_A\nabla_Br=(\Box r)q_{AB}/2$ なら $\partial_A F=2(\nabla_A\nabla_Br)\nabla^Br=(\Box r)\partial_A r$ である。また角方向の方程式は $1-F=r\Box r$ である。両者を微分
$$
\partial_A m=\tfrac12[(1-F)\partial_A r-r\partial_A F]
$$
へ代入すると零になる。この証明は $F$ の符号を使わず、$F=0$ の点も含む。例えば非真空な模型で $F=1-2m(v)/r$ とすると、$\partial_v m\ne0$ があり得る。ここで真空の $AB$ 方程式が破れていることを調べれば、球対称だけから時間依存を排除できない理由が分かる。球の面積から定義した $r$ を使うため、$m$ の値は角座標の付け方に依存しない。

### 軌道商の曲率を縮約する手順

軌道商の添字を $A,B$、球面の添字を $i,j$ とする。$g=q+r^2\gamma$ で接続の混合成分は本文に示した二種類だけである。曲率の $AB$ 成分の縮約で、球面方向の追加項は
$$
-\partial_B(2\partial_A r/r)
+\Gamma^C{}_{AB}(q)\frac{2\partial_C r}{r}
-\frac{2\partial_A r\partial_B r}{r^2}
=-\frac{2}{r}\nabla_A\nabla_Br
$$
である。角方向の式では、外側の商添字を持つ曲率が $-r\,\nabla_A\nabla_Br\,\gamma_{ij}$ となる。これは $\partial_A(-r\nabla^B r)$ を微分し、$\Gamma^i{}_{Aj}=(\partial_A r/r)\delta^i_j$ との積で一階微分の積を消す結果である。一方、球面内の曲率には単位球の成分と $\Gamma^A{}_{ij}\Gamma^k{}_{Al}$ の積がある。この差は
$$
R^i{}_{jkl}=(1-f)(\delta^i_k\gamma_{jl}-\delta^i_l\gamma_{jk})
$$
になる。球面の二次元でこれを縮約すると $(1-f)\gamma_{ij}$、先ほどの商の曲率を縮約すると $-r\Box_qr\,\gamma_{ij}$ である。従って $R_{ij}=(1-f-r\Box_qr)\gamma_{ij}$ となり、混合した Ricci 成分は球対称性から零である。この形を使えば、座標の選択に依存せず $f$ と質量の式を導ける。

### 面積半径の微分が零になる点の近傍

$r>0$ の球軌道で $\dd r=0$ の点を扱う。この点で $f=0$ なので質量関数から $r=2m$、従って $m>0$ である。二次元 Lorentz 計量には、異なる二本の光的方向の積分曲線を用いて局所の二重零座標 $u,v$ を選べる。原点を考えている点とし、
$$
q=-2e^{2\sigma(u,v)}\dd u\,\dd v
$$
と書く。$uu,vv$ の Hessian 条件は $r_{uu}-2\sigma_u r_u=0$、$r_{vv}-2\sigma_v r_v=0$ である。従って $v=0$ では $e^{-2\sigma}r_u$、$u=0$ では $e^{-2\sigma}r_v$ が軸上で一定であり、原点では両者が零である。よって両方の軸上で $r=2m$ を得る。さらに $u,v$ をそれぞれ単独で再パラメータ化し、軸上で $\sigma=\sigma_0$ とできる。例えば $U'(u)=e^{2(\sigma(u,0)-\sigma_0)}$、$V'(v)=e^{2(\sigma(0,v)-\sigma_0)}$ と取れば原点の微分が $1$ で、この条件を両軸で満たす。

残る方程式は
$$
r_{uv}=-\frac{e^{2\sigma}+2r_u r_v}{2r},\qquad
\sigma_{uv}=-\frac{r_{uv}}{r}.
$$
第一式は $1-f=r\Box_qr$ に $f=-2e^{-2\sigma}r_ur_v$ と $\Box_qr=-2e^{-2\sigma}r_{uv}$ を入れたものである。第二式は $R_{uv}[q]=-2\sigma_{uv}$ を $R_{uv}[g]=0$ に代入した。両軸上の値を与えると、この二つの式の滑らかな解は小さい箱内で一意になる。実際、第一式の右辺を $F(r,r_u,r_v,\sigma)$ とすれば、$r_u$ は $v$ に関する $F$ の積分、$r_v$ は $u$ に関する積分であり、$r$ と $\sigma$ はそれぞれ $F$ と $-F/r$ の二重積分で求まる。二解の差の最大値を取ると、$r>0$ の小箱では右辺が Lipschitz 連続なので、差は箱の幅に比例する定数倍の自分自身以下になる。幅を十分小さくすると定数が $1$ 未満になり、差は零である。

比較する解は $U,V$ を用いて
$$
UV=(1-r/(2m))e^{r/(2m)},\qquad
q=-\frac{32m^3}{r}e^{-r/(2m)}\dd U\,\dd V
$$
と置いたものである。右辺の $r$ に関する一階微分は $-r e^{r/(2m)}/(4m^2)$ で、$r=2m$ でも非零であるため $r(U,V)$ は滑らかである。例えば外部の一つの象限で $r_*=r+2m\log(r/(2m)-1)$ とし、$U=-e^{-(t-r_*)/(4m)}$、$V=e^{(t+r_*)/(4m)}$ と取ると :eqref[eq-gr-schwarzschild] からこの式を得る。他の象限でも、$UV\ne0$ で $t=2m\log|V/U|$ と $r$ を使えば同じ Schwarzschild 形に戻る。それゆえ真空方程式は各象限で成立し、その滑らかさから軸上でも成立する。両軸で $r=2m$、面積係数は定数である。$U,V$ の定数倍で $\sigma_0$ を合わせれば、上の一意性から、任意の球対称真空解はこの点の近傍でも Schwarzschild の正則な局所表示になる。これで $\dd r=0$ を含む $r>0$ の球軌道上で Birkhoff の局所主張を示せた。大域的な延長と各領域の因果関係は [Schwarzschild ブラックホール](/physics/general-relativity/black-holes/) で扱う。

### 真空と球対称性を分けて検査する

$q=-[1-2m(v)/r]\dd v^2+2\dd v\,\dd r$ を考える。球面の面積半径は同じだが、$\partial_vf=-2m'(v)/r$ であり、混合する時刻依存性を検査できる。逆計量から $\Gamma^r{}_{vv}=(ff_r-f_v)/2$、$\Gamma^r{}_{vr}=-f_r/2$ である。従って $AB$ 部分の Ricci は $R_{AB}=R_{AB}[q]+2\Gamma^r{}_{AB}/r$ を満たす。$R_{AB}[q]=-f_{rr}q_{AB}/2$ であり、$f=1-2m(v)/r$ を代入すると $R_{vr}=0$、$R_{rr}=0$、$R_{vv}=2m'(v)/r^2$ が残る。角成分も零なので $R=0$ である。よって $T_{vv}=m'(v)/(4\pi Gr^2)$、他の共変成分は零である。$m'(v)>0$ なら未来向きの入射光流を表す。半径の変化から質量の変化を読み取れるが、真空を仮定したときのみ $m$ が定数になる。

### 宇宙項がある場合と一定半径の例外

球対称で $R_{ab}=\Lambda g_{ab}$ とする。軌道商の $AB$ 成分の跡なし部分は同じだから $\nabla_A\nabla_Br=(\Box_qr)q_{AB}/2$ である。角成分は $1-f-r\Box_qr=\Lambda r^2$ になり、$\dd f=(1-f-\Lambda r^2)\dd r/r$ を得る。これを使うと、質量関数
$$
m_\Lambda=\frac r2\left(1-f-\frac{\Lambda r^2}{3}\right)
$$
の微分は零で、$f=1-2m_\Lambda/r-\Lambda r^2/3$ である。$\dd r\ne0$ なら前と同様の Killing 場 $K$ が零でないので、局所の計量をこの $f$ で書ける。一方、$\dd r=0$ が開領域全体で成立する場合は $\Lambda r^2=1$ となり、$\Lambda>0$ で可能である。球半径は $1/\sqrt\Lambda$、商の計量は $R_{AB}[q]=\Lambda q_{AB}$ である。この直積解は、$r$ を座標にできる枝とは別である。$r$ が一定の四次元計量を無理に面積半径座標に入れるとこの解を落とす。

:::example[一定半径を持つ球対称真空を構成する]
$L=1/\sqrt\Lambda$、$|x|<L$ とし、
$$
q=-(1-x^2/L^2)\dd t^2+(1-x^2/L^2)^{-1}\dd x^2
$$
を選ぶ。この二次元計量のスカラー曲率は $-(1-x^2/L^2)''=2/L^2$ である。二次元では Ricci はその半分を計量に掛けたものだから $R_{AB}[q]=q_{AB}/L^2$ である。球面の半径を $L$ に固定すると、その Ricci も計量の $1/L^2$ 倍になる。直積は混合成分を持たないので四次元の $R_{ab}=\Lambda g_{ab}$ を満たす。Schwarzschild の $\Lambda=0$ の証明で半径一定の開領域を排除できたのは、そのとき角方向の正曲率を打ち消す項がないためである。
:::

### 正則な中心を含む真空

球軌道が縮んで中心に達する場合、中心自体には上で用いた二次元軌道の仮定を適用しない。中心へ近づける 中心を除いた真空近傍で得た質量 $m$ は一定であり、$m\ne0$ なら曲率二乗 $48m^2/r^6$ が発散する。滑らかな中心で曲率は有限でなければならないため、この枝では $m=0$ が必要である。$m=0$ の近傍では $g$ は Minkowski の球座標表示であり、Cartesian 座標にすればその中心は滑らかに含まれる。この結論は真空の中心が対象である。物質を含む星の中心と外の真空で同じ質量関数が一定だと仮定することはできない。

## Schwarzschild 時空の測地線と有効ポテンシャル

赤道に初期位置と速度を取る。反射 $\theta\mapsto\pi-\theta$ が等長変換なので、測地線の一意性により赤道内に留まる。回転を選べば一般の非動径運動もこの形で表せる。:eqref[eq-gr-killing-charge] から $E=f\dot t$、$L=r^2\dot\phi$ を保存する。規格化 $g(\dot x,\dot x)=-\kappa$ で $\kappa=1$（時間的）または $0$（光的）とすると

$$
\dot r^2+V_\kappa(r)=E^2,\qquad
V_\kappa=f(\kappa+L^2/r^2).
\label{eq-gr-schwarz-potential}
$$

動径方向で $\ddot r=-V_\kappa'/2$ となることは、動径の測地線式に $E,L$ を代入して得る。従って折返し点の前後でもこの加速度式は成立する。例えば $L=0,E=1$ なら $\dot r=-\sqrt{2m/r}$ を選べ、$r_1$ から $r_2$ までの固有時は $2(r_1^{3/2}-r_2^{3/2})/(3\sqrt{2m})$ である。一方、$r>2m$ に静止した加速度計は $\sqrt{g_{rr}}\Gamma^r{}_{tt}(u^t)^2=m/(r^2\sqrt f)$ を測る。

### 有限離心率の束縛軌道

円軌道の微小振動以外にも、二つの転回半径から保存量を決められる。半通径の無次元値を $p$、離心率を $0<e<1$ とし、$r_{\rm p}=mp/(1+e)$、$r_{\rm a}=mp/(1-e)$ とする。軌道面の向きで $L>0$ を選び、$x=m/r$、$\ell=L/m$ と置けば $V_1=(1-2x)(1+\ell^2 x^2)$ である。二つの転回点 $x_\pm=(1\pm e)/p$ で $V_1$ が等しいので、その差を $x_+-x_-$ で割る。$x_++x_-=2/p$ と $x_+^2+x_+x_-+x_-^2=(3+e^2)/p^2$ により
$$
\ell^2=\frac{p^2}{p-3-e^2},\qquad E^2=\frac{(p-2)^2-4e^2}{p(p-3-e^2)}
$$
を得る。後式は前式を一方の転回点の $V_1$ に代入したものである。$E^2-V_1$ の三次項と二次項から三根の和は $1/2$ なので、残りの根は $x_3=1/2-2/p$ となる。動径運動には $E^2-V_1=2\ell^2(x-x_-)(x-x_+)(x-x_3)\ge0$ が必要である。二つの転回点の間を往復するには $x_3>x_+$、すなわち $p>6+2e$ が必要であり、等号では近点側の根が重なって軌道が近点に漸近する。

近点からの位相を $\chi$ として $x=(1+e\cos\chi)/p$ と置く。:eqref[eq-gr-schwarz-potential] と $\dot\phi=L/r^2$ により $(\dd x/\dd\phi)^2=(E^2-V_1)/\ell^2$ なので、三根表示との比較から
$$
\dv{\phi}{\chi}=\sqrt{\frac{p}{p-6-2e\cos\chi}}
$$
を得る。$\chi$ を $0$ から $2\pi$ まで積分し $2\pi$ を引けば、有限離心率での一周の近点移動になる。$p\gg1$ なら積分の一次項は $\int_0^{2\pi}(3+e\cos\chi)\dd\chi/p=6\pi/p$ である。例えば $e=1/2$ では、往復軌道の境界は $p=7$、円の極限 $e\to0$ では $p=6$ に戻る。

### 転回点の型を一階式と二階式で判定する

$E^2=V_\kappa(r_0)$ かつ $V_\kappa'(r_0)\ne0$ なら、アフィン径数 $s_0$ で $r=r_0,\dot r=0$ として
$$
r(s)=r_0-\frac{V_\kappa'(r_0)}4(s-s_0)^2+\order{(s-s_0)^3}
$$
である。$V'$ が負なら最小半径、正なら最大半径になる。これに対し二重根では一階式の右辺を平方根にすると、近くの許容領域では $|\dot r|$ が距離 $|r-r_0|$ に比例する。従って転回点までの径数積分は対数的に発散し、一般に有限径数で折り返す解にはならない。

:::example[無限遠で静止する粒子の動径条件]
$E=1$ とし、$\ell=L/m$、$x=m/r$ を用いると $1-V_1=2x-\ell^2x^2+2\ell^2x^3$ である。$x>0$ の根は $2\ell^2x^2-\ell^2x+2=0$ で、判別式は $\ell^2(\ell^2-16)$ になる。従って $\abs\ell<4$ では折返しがなく、$\abs\ell>4$ で外から入る粒子は小さい $x$ の根で戻る。$\abs\ell=4$ では二重根 $x=1/4$ となり、$r=4m$ の不安定円軌道に漸近する。
:::

### 座標時刻と落下時計を別々に積分する

$E=1,L=0$ の時間的粒子の内向きの脚で、$s=\sqrt{r/(2m)}>1$ と置く。固有時に関して $\dot r=-1/s$、$\dot t=1/f=s^2/(s^2-1)$ である。従って $\dd\tau/\dd s=-4m s^2$、$\dd t/\dd s=-4m s^4/(s^2-1)$ となる。後者の分数を $s^2+1+1/(s^2-1)$ に分けると
$$
t(s)=-4m\!\left[\frac{s^3}{3}+s+\frac12\log\frac{s-1}{s+1}\right]+C
$$
である。$s\to1^+$ で対数が発散する一方、固有時の差は $4m(s_1^3-s^3)/3$ で有限である。進入光の時刻 $v=t+r_*$ で見れば
$$
\dot v=\frac{1+\dot r}{f}=\frac{s}{s+1}
$$
となり、極限は $1/2$ である。従って同じ落下で発散するのは、外部の静止時刻を使った座標の値である。

:::example[光と粒子の径数を比較する]
動径入射光では正則な $v=t+r_*$ を一定に取れ、$E=-k\cdot\partial_v=-\dot r>0$ なので $r(\lambda)=r_1-E(\lambda-\lambda_1)$ である。同じ半径から出る $E=1$ の時間的粒子では $\dot r=-\sqrt{2m/r}$ である。動径光のアフィン径数の尺度は任意だが、時間的粒子では固有時の規格化がある。例えば $r_1=8m$ から $2m$ まで、光の径数差を $E=1$ で $6m$ と表せる一方、粒子が測る固有時は $28m/3$ である。二つの数をそのまま同じ種類の時計の読みとして比較することはできない。
:::

### 捕獲と束縛を有効ポテンシャルで分類する

$m>0$、$L\ne0$ で $V_1'=(2/r^4)(mr^2-L^2r+3mL^2)$ である。$L^2<12m^2$ なら二次式は正で、外部の $V_1$ は単調増加する。$L^2>12m^2$ なら極値の半径は
$$
r_\mp=\frac{L^2}{2m}\left(1\mp\sqrt{1-\frac{12m^2}{L^2}}\right).
$$
$V_1'$ は外側と内側で正、その間で負なので $r_-$ は極大、$r_+$ は極小である。$V_1(2m)=0$、遠方で $V_1\to1$ である。例えば $V_1(r_+)<E^2<\min(1,V_1(r_-))$ なら、$r_+$ を挟んで二つの転回点を持つ外側の許容領域がある。内側にも地平線につながる別の許容領域があり、保存量だけが同じでも初期位置で運動が変わる。$E^2>1$ の粒子が無限遠から入るときは、極大の高さが $E^2$ より高ければ散乱し、低ければ捕獲される。

:::example[同じ角運動量で安定と不安定の二つの円を作る]
$L=4m$ とすると $r_-=4m$、$r_+=12m$ である。$r=4m$ では $V_1=1$、$r=12m$ では $V_1=25/27$ である。従って例えば $E^2=26/27$ なら外側に束縛軌道があり、無限遠からは同じエネルギーの粒子は入ってこられない。円軌道と転回点の存在を分けて調べたため、エネルギー値だけから運動を決める誤りを避けられる。
:::

### 無限遠から来る質量粒子の捕獲境界

遠方で速さ $v_\infty>0$ を持つ入射粒子では $E=(1-v_\infty^2)^{-1/2}>1$ であり、角運動量の大きさを $L=b\sqrt{E^2-1}$ と書ける。実際、遠方の直線近似で運動量の単位質量当たりの大きさは $E v_\infty=\sqrt{E^2-1}$、直線と中心との距離が $b$ である。捕獲の境界では外側から来る粒子が有効ポテンシャルの極大に漸近する。$x=r/m$ とおくと円軌道の式は $E^2=(x-2)^2/[x(x-3)]$ であり、解は
$$
x_c=\frac{8}{4-3E^2+E\sqrt{9E^2-8}},\qquad
b_c=\frac{m x_c}{\sqrt{(x_c-3)(E^2-1)}}.
$$
ここで $3<x_c<4$ であることは、二次式 $(E^2-1)x^2+(4-3E^2)x-4$ の正根が $3$ と $4$ で符号を変えることで確かめられる。従って不安定円軌道である。固定した $E$ で $L$ を増すと $f>0$ の各 $r$ で $V_1$ が増えるので、$b>b_c$ では外から来た粒子は折り返し、$b<b_c$ では地平線まで進む。小さい $v_\infty$ で $L_c\to4m$ のため $b_c\simeq4m/v_\infty$ となる。一方、$E\to\infty$ では $(x_c-3)E^2=(x_c-2)^2/x_c\to1/3$、従って $b_c\to3\sqrt3m$ となる。

:::example[不安定円から入射条件を逆算する]
$r_c=7m/2$ では $E^2=9/7$、$L_c^2=49m^2/2$、$b_c^2=343m^2/4$ である。この入射粒子の遠方速度は $v_\infty=\sqrt{2}/3$ である。円での $L_c$ をそのまま衝突径数にすると、遠方の速さの因子を落とす。同じ $E$ で $L$ がわずかに大きい場合は、ポテンシャルの極大が $E^2$ を超えて外側に転回点が生じる。小さい場合は極大付近を通過して落下する。
:::

### 近点と遠点を周期で区別する

有限離心率の束縛軌道で $\chi$ を使うと、$\dot\phi=L/r^2$ より
$$
\frac{\dd\tau}{\dd\chi}
=\frac{m p^{3/2}\sqrt{p-3-e^2}}
{(1+e\cos\chi)^2\sqrt{p-6-2e\cos\chi}} .
$$
座標時刻の周期はこれに $E/f(r(\chi))$ を掛けて $\chi=0$ から $2\pi$ まで積分したものである。状況に応じて、方位角の増分、固有時の周期、遠方の時刻の周期を別々に計算する。例えば $e\to0$ における座標周期は、円軌道の $r_0=mp$ を使って $2\pi/\Omega_r$ へ近づく。回転角を $2\pi$ 固定した時間ではなく、動径が一往復する時間を表している。

離心率 $0<e<1$ を固定し、境界 $p_s=6+2e$ に対して $p=p_s+\delta$ とする。$\dd\phi/\dd\chi$ の分母の二乗は $\delta+2e(1-\cos\chi)$ である。$\chi$ が小さい部分では $\delta+e\chi^2$ が先頭であり、積分には $e^{-1/2}\operatorname{arsinh}(\chi\sqrt{e/\delta})$ が現れる。区間の反対端 $2\pi$ でも同じ寄与がある。両端の寄与を合わせると全角度は $-\sqrt{p_s/e}\log\delta+\order{1}$ で発散する。端から離れた区間で分母は正の下限を持ち、その寄与は有限である。端で $\cos\chi$ を二次に置き換えた誤差の積分も、次の項が $\chi^4$ で分母が三乗にしかならないため有界である。従って外側の往復が残っていても、近点付近では不安定円軌道を何周もする場合がある。この極限と $e=0$ の小振動の極限は、分母の発散する部分が異なる。

## 最内安定円軌道

以下、軌道を論じるときは $m>0$ とし、まず外部領域 $r>2m$ を用いる。時間的円軌道では $V_1'=0$ と $E^2=V_1$ の両方が必要であり、

$$
L^2=\frac{mr^2}{r-3m},\qquad E^2=\frac{(1-2m/r)^2}{1-3m/r},\qquad
\left(\frac{\dd\phi}{\dd t}\right)^2=\frac{m}{r^3}.
\label{eq-gr-circular-orbit}
$$

第一式は $V_1'=2m/r^2-2L^2/r^3+6mL^2/r^4=0$ を解いたもの、第三式は $Lf/(Er^2)$ を二乗したものである。$r>3m$ が必要で、$L$ を固定して二度微分してから値を代入すると $V_1''=2m(r-6m)/(r^3(r-3m))$ である。線形の小振動が安定な領域は $r>6m$ であり、その境界を :term[最内安定円軌道]{en="innermost stable circular orbit" yomi="さいないあんていえんきどう"} と呼ぶ。$r=6m$ では $L=2\sqrt3m$、$E=\sqrt{8/9}$、静止時の質量エネルギーに対する結合の割合は $1-\sqrt{8/9}\simeq0.0572$ となる。$r=4m$ では $E=1$ だが、$V_1''<0$ なので不安定である。

### 動径振動と方位角の振動数

安定円軌道 $r_0>6m$ の近くで、$L$ を円軌道の値に固定して $\delta r$ を一次まで残すと $\delta\ddot r=-V_1''(r_0)\delta r/2$ である。固有時の振動数は $\sqrt{V_1''/2}$ で、座標時刻の振動数はこれに円軌道上の $\dd\tau/\dd t=f/E$ を掛ける。時刻換算の一次補正と $\delta r$ の積は二次なので捨てられる。:eqref[eq-gr-circular-orbit] より
$$
\Omega_r^2=\frac{m}{r_0^3}\left(1-\frac{6m}{r_0}\right),\qquad
\Omega_\phi^2=\frac{m}{r_0^3}.
$$
近点から次の近点までの方位角は $2\pi\Omega_\phi/\Omega_r$ である。$m/r_0\ll1$ で $2\pi$ を引くと $6\pi m/r_0$ が先頭になり、有限離心率の式の円極限と一致する。

:::example[同じ周波数を使う Newton 近似が破れる半径]
$r_0=8m$ では $1-6m/r_0=1/4$ なので $\Omega_r=\Omega_\phi/2$ である。微小離心率で動径が一周期する間に方位角は $4\pi$ 進み、近点の進みは $2\pi$ になる。これは弱い場の一次展開をそのまま代入した角度と異なるため、比較には背景半径での正確な振動数を使う。
:::

### 境界の円軌道では三次項が残る

$L^2=12m^2$ のとき $V_1'$ の分子は $2m(r-6m)^2$ である。従って $r=6m$ は極小でも極大でもなく停留変曲点である。$r_0=6m$、$y=r-r_0$ として積分すると、
$$
V_1(r_0+y)=\frac89+\frac{2m}{3r_0^4}y^3+\order{y^4/m^4}
$$
である。線形の動径振動数は零であり、最初の加速度は $-m y^2/r_0^4$ で内向きになる。例えば $y<0$ にわずかにずらし、速度を零にして $L$ を固定すると、二階式から粒子は内向きに動き始める。従って $r=6m$ を最内安定円軌道と呼ぶときは、安定円軌道の族の下限であるという意味である。

### エネルギーと角運動量の関係を円軌道に沿って微分する

未来向きで $L>0$ の円軌道の族に沿うと $E^2=V_1(r,L(r))$ である。偏微分 $(\partial_r V_1)_L=0$ だから、
$$
2E E'=(\partial_L V_1)_r L'=\frac{2fL}{r^2}L',
\qquad E'=\Omega_\phi L'.
$$
これは $r=6m$ の前後で $E'$ と $L'$ の符号が同時に変わる理由でもある。安定な円軌道の列で角運動量を少し失うと、別の円へ移るためには対応するエネルギーも減る。実際の降着流の放射効率を $1-E(6m)$ と同一視するには、遠方で粒子がほとんど静止し、円軌道に沿ってエネルギーを失い、最内安定円軌道より内側で失うエネルギーを無視する模型を使う。円軌道の比エネルギー自体はこれらの流体の仮定なしに求められる。

## 光子球と光線の軌道

光では $V_0=L^2(r-2m)/r^3$ とする。$V_0'=2L^2(3m-r)/r^4$ だから $m>0$ かつ $L\ne0$ なら $r=3m$ に極大がある。ここにある不安定な円軌道の球面を :term[光子球]{en="photon sphere" yomi="こうしきゅう"} と呼ぶ。無限遠から来る光の $b=L/E$ を衝突径数とすると、極大値と $E^2$ が等しい臨界では $|b|=3\sqrt3m$ である。$b$ の絶対値がこれより大きければ外側に折返しがあり、小さければ折り返さない。軌道方程式は $u=1/r$ を使うと $(\dd u/\dd\phi)^2=b^{-2}-u^2+2mu^3$ である。動径の二階式 $\ddot r=-V_0'/2$ と $\dd u/\dd\phi=-\dot r/L$ を使うと $u''+u=3mu^2$ が得られる。この導出は転回点と円軌道でも有効である。$m=0$ で $u=\sin\phi/b$ は直線である。$m>0$ での摂動解を次章で光の湾曲に用いる。

### 臨界衝突径数の外側で光が巻き付く理由

$b>0$ の場合、$b_c=3\sqrt3m$ とし、$b=b_c(1+\delta)$、$0<\delta\ll1$ とする。$V_0/E^2=b^2 f/r^2$ の最大付近で $r=3m+y$ と置くと、
$$
1-b^2f/r^2=-2\delta+\frac{y^2}{3m^2}
+\order{\delta^2}+\order{\delta y^2/m^2}+\order{y^3/m^3}
$$
となる。従って外側転回点は先頭で $y_0=m\sqrt{6\delta}$ である。$\dd\phi/\dd r=b/(r^2\sqrt{1-b^2f/r^2})$ にこの展開を入れると、近傍での角度の先頭は $\int_{y_0}^{y_1}\dd y/\sqrt{y^2-y_0^2}=\operatorname{arcosh}(y_1/y_0)$ となる。対数項の係数を確かめるには、$h=f/r^2$ とし、外側近傍で $z=\sqrt{1-b_c^2h}$ を積分変数に取る。上の二次項から $z=y/(\sqrt3m)+\order{y^2/m^2}$ である。$z$ は $y=0$ でも片側から滑らかで、そこで微分が非零になる。実際、平方根の中から $y^2$ をくくれば、残る滑らかな因子は正である。すると
$$
1-b^2 h=(1+\delta)^2(z^2-z_0^2),\qquad
z_0^2=1-(1+\delta)^{-2}
$$
となり、厳密な転回点は $z=z_0$ である。角度の積分に残る係数 $B(z)=(b_c/r(z)^2)\dv{r}{z}$ は $B(0)=1$ なので、$1+\order{z}$ と書ける。固定した $z_1>0$ まで積分すると、$1$ の項は $\operatorname{arcosh}(z_1/z_0)=-\log\delta/2+\order{1}$、$\order{z}$ の寄与は $\int_{z_0}^{z_1}z\,\dd z/\sqrt{z^2-z_0^2}\le z_1$ により有界である。近傍外では分母は零にならず、遠方で被積分関数は $r^{-2}$ に比例するので残る角度も有界である。外向きの脚も同じ角度を持つため、全角度は $-\log\delta+\order{1}$ である。従って小さい $\delta$ に対して一次の弱重力偏向では表せない多数回の巻き付きが現れる。

:::example[観測地での影の角半径]
$r_o>3m$ の静止枠で、中心方向からの角度を $\alpha$ とする。光の局所エネルギーは $E/\sqrt f$、横運動量の大きさは $|L|/r_o$ なので $\sin\alpha=|b|\sqrt f/r_o$ である。臨界を代入し $r_o=6m$ とすると $\sin\alpha=1/\sqrt2$、従って内向きの半球での影の角半径は $\pi/4$ になる。遠方の小角近似 $b_c/r_o$ には静止枠への変換の $\sqrt f$ が含まれないため、この距離ではそのまま使えない。
:::

### 散乱角を一つの積分で表す

無限遠から来る $b>b_c$ の光の外側転回半径を $r_{\min}$ とする。$\dot r$ の符号が異なる二脚は角度の大きさが等しいので、偏向角は
$$
\Delta\phi=2\int_{r_{\min}}^\infty\frac{b\,\dd r}{r^2\sqrt{1-b^2 f/r^2}}-\pi .
$$
このとき $b=r_{\min}/\sqrt{1-2m/r_{\min}}$ である。$z=r_{\min}/r$、$\epsilon=m/r_{\min}$ と置くと半分の角度は
$$
I=\int_0^1\frac{\dd z}{\sqrt{(1-z^2)-2\epsilon(1-z^3)}} .
$$
転回点を固定した $z=1$ に移したため、$\epsilon\ll1$ なら $\epsilon(1-z^3)/(1-z^2)$ は積分全体で一様に小さい。一次項は $\int_0^1(1-z^3)/(1-z^2)^{3/2}\dd z$ である。$z/\sqrt{1-z^2}$ の微分が $(1-z^2)^{-3/2}$ なので、部分積分して $3\int_0^1 z^3/\sqrt{1-z^2}\dd z=2$ となる。従って $\Delta\phi=4m/r_{\min}+\order{(m/r_{\min})^2}=4m/b+\order{(m/b)^2}$ である。$b$ を固定して動く転回点をそのまま展開すると端点で発散する項が現れ得るため、先に転回点を固定する。観測との比較は [一般相対論の検証](/physics/general-relativity/solar-system-tests/) で扱う。

### 光的な円軌道の動径摂動

$b=b_c$ のとき $V_0''(3m)=-2L^2/(81m^4)$ である。動径の二階式を $r=3m+y$ で一次化すると $y''_\lambda=L^2 y/(81m^4)$ となる。よってアフィン径数当たりの指数率は $|L|/(9m^2)$ である。円軌道で $E=|L|/(3\sqrt3 m)$、$\dd t/\dd\lambda=3E$ だから、座標時刻での指数率は $1/(3\sqrt3 m)$ になる。同じ軌道の $\dd\phi/\dd t=L/(27m^2E)$ の絶対値も同じである。円軌道からずれた光が近傍に長く滞在する場合、方位角が一周進むごとに線形なずれの増加モードは $e^{2\pi}$ 倍になる。この比較は指数率を線形化した近傍での結果である。

:::example[静止観測者から見た脱出円錐]
ブラックホールの外部 $2m<r_o<3m$ で、外向きの動径方向から角度 $\beta$ を持つ光を放つ。$b=r_o\sin\beta/\sqrt{f(r_o)}$ の絶対値が $b_c$ より小さく、かつ外向きなら光は光子球の極大を越えて脱出する。外向きでも $b>b_c$ なら内側の転回点で戻る。従って脱出する方向の立体角は $2\pi(1-\cos\beta_c)$、$\sin\beta_c=b_c\sqrt f/r_o$ である。例えば $r_o=9m/4$ なら $f=1/9$、$\sin\beta_c=4\sqrt3/9$、$\cos\beta_c=\sqrt{33}/9$ となる。半径を $2m$ に近づけると静止枠での脱出円錐は狭くなる。ただし、その位置を保つ静止観測者には前に求めた固有加速度が必要である。
:::

### 光路を空間の計量で表す

静的な外部で $f>0$ とすると、光的条件は赤道で
$$
\dd t=\sqrt{f^{-2}\dd r^2+(r^2/f)\dd\phi^2}
$$
である。右辺の正定値計量を光学計量と呼ぶ。光的アフィン接線では $E=f\dot t$ が保存するので、光学計量での速さは座標時刻 $t$ に関して $1$ であり、角の運動量は $(r^2/f)(\dd\phi/\dd t)=L/E$ である。光学半径を用いればこの幾何は、動径の固有長を $r_*$ とする回転対称な空間の測地線の問題になる。逆に光学計量の角運動量と単位速度条件から
$$
\left(\frac{\dd r}{\dd t}\right)^2=f^2\left(1-\frac{b^2 f}{r^2}\right)
$$
を得る。これは光の動径一階式に $\dot t=E/f$ を入れたものと等しい。さらに光学計量のラグランジアンを $r$ で変分した式は、角運動量と単位速度条件を代入すると動径の二階の式を与える。転回点では一階式だけを解いて枝を任意に張り替えず、この二階の式で連続な運動を選ぶ。

:::example[光学計量の円周の極小と光子球]
光学計量で動径の長さを $\ell_*=r_*$ とすると、球面の光学半径は $R_{\mathrm{opt}}=r/\sqrt f$ である。ニュートンの粒子の有効ポテンシャルにおける遠心項と同様に、光の動径式には $b^2/R_{\mathrm{opt}}^2$ が入る。$\dd R_{\mathrm{opt}}/\dd r=f^{-3/2}(1-3m/r)$ だから $r=3m$ で $R_{\mathrm{opt}}$ が極小となり、値は $3\sqrt3m$ である。円周の半径が極小でも、固定された $b$ に対する項 $b^2/R_{\mathrm{opt}}^2$ は極大なので光の円軌道は不安定である。
:::

### 有限半径での動径光の周波数

外部の静止枠を $n=\partial_t/\sqrt f$、外向き単位軸を $e_r=\sqrt f\partial_r$ とする。入射光は $k=(E_\gamma/\sqrt f)(n-e_r)$ と分解できる。$E_{\rm fall}=1,L=0$ の落下粒子の枠での動径速度は $v=-\sqrt{2m/r}$、Lorentz 因子は $\gamma=1/\sqrt f$ である。従って落下者の観測エネルギーは
$$
-g(k,u)=\frac{E_\gamma}{f}(1-\sqrt{2m/r})
=\frac{E_\gamma}{1+\sqrt{2m/r}}.
$$
静止枠の周波数が $f^{-1/2}$ で発散しても、同方向に落ちる観測者では $r\to2m$ において $E_\gamma/2$ に近づく。この結果は $v=t+r_*$ の表示でも得られる。入射光の共変運動量は $k_a=-E_\gamma(\dd v)_a$ なので、測定値は $E_\gamma\dot v$ であり、その半径に対応する $\dot v=s/(s+1)$ を代入すればよい。例えば $r=8m$ では $2E_\gamma/3$ である。同じ測定を二つの枠から計算することで、光と落下者のそれぞれの規格化を照合できる。

### 臨界値での吸収と散乱の境目を積分する

$b=3\sqrt3m$ とし、$x=r/m$ を用いる。動径の平方根の中身は $(x^3-27x+54)/x^3=(x-3)^2(x+6)/x^3$ となる。$r>3m$ から内向きに進めば $\dot r$ は $-(x-3)$ に正の滑らかな係数を掛けたものなので、光は有限なアフィン径数で $3m$ に達しない。例えば $x$ を $3+\varepsilon$ から $3+\varepsilon/2$ まで下げるには、$\varepsilon\to0^+$ で一定の正の径数差が必要であり、そのような区間を無限に通る。円軌道と同じ $b$ でも、最初から円上にある場合と、無限遠側から近づく場合では世界線が異なる。衝突径数の等号を、散乱か捕獲のどちらかに含めるとこの第三の場合を見落とす。

## 演習問題

::::exercise{level=basic}
$m=3\,\mathrm{km}$ の地平線候補と光子球の面積半径を求めよ。

:::solution
$f=0$ は $r=2m=6\,\mathrm{km}$、$V_0'=0$ は $r=3m=9\,\mathrm{km}$ である。
:::
::::

::::exercise{level=basic}
$r=8m$ の時間的円軌道の $L^2$ を求めよ。

:::solution
:eqref[eq-gr-circular-orbit] から $L^2=m(8m)^2/(5m)=64m^2/5$ である。$8m>6m$ なので安定である。
:::
::::

::::exercise{level=basic}
$r=4m$ の静止時計の率を無限遠に対して求めよ。

:::solution
$f=1/2$ なので $\dd\tau/\dd t=\sqrt f=1/\sqrt2$ である。時刻 $t$ は無限遠で静止時計の固有時に一致する。
:::
::::

::::exercise{level=basic}
質量に対応する長さが $m>0$ のブラックホールに $b=6m$ で光を入射すると捕獲されるか。

:::solution
$6m>3\sqrt3m$ である。無限遠から入る光はポテンシャルの極大に達する前に外側で折り返す。
:::
::::

::::exercise{level=advanced}
円軌道 $r_0>6m$ の固有時での動径小振動数を求めよ。

:::solution
$\ddot r=-V_1'/2$ を $r=r_0+\delta r$ で一次化すると $\delta\ddot r+V_1''(r_0)\delta r/2=0$ である。従って角振動数の二乗は $m(r_0-6m)/(r_0^3(r_0-3m))$ である。
:::
::::

::::exercise{level=advanced}
$E=1,L=0$ の落下で $r=8m$ から $r=2m$ までの固有時を求めよ。

:::solution
$\dd\tau=-\sqrt{r/(2m)}\dd r$ を積分し、$2((8m)^{3/2}-(2m)^{3/2})/(3\sqrt{2m})=28m/3$ を得る。
:::
::::

::::exercise{level=advanced}
静的真空方程式の $f$ に $1-2m/r+q/r^2$ を代入し、$q\ne0$ の場合の $G^t{}_t$ を求めよ。

:::solution
$f'=2m/r^2-2q/r^3$ なので $f-1+rf'=-q/r^2$、従って $G^t{}_t=-q/r^4$ である。
:::
::::

::::exercise{level=advanced}
$E=1$ の円軌道の位置と安定性を求めよ。

:::solution
$(1-2m/r)^2=1-3m/r$ から有限半径では $4m^2/r^2=m/r$、よって $r=4m$ である。$r-6m<0$ なので不安定である。
:::
::::

::::exercise{level=advanced}
円軌道の比エネルギー $E=(r-2m)/\sqrt{r(r-3m)}$ を $r>3m$ で微分し、結合エネルギー $1-E$ が極大になる半径を求めよ。その円軌道の局所速度を静止枠に対して求めよ。

:::solution
対数微分より
$$
\frac{E'}E=\frac1{r-2m}-\frac1{2r}-\frac1{2(r-3m)}
=\frac{m(r-6m)}{2r(r-2m)(r-3m)} .
$$
従って $E$ は $3m<r<6m$ で減少、$r>6m$ で増加し、$1-E$ は $6m$ で極大である。静止枠での速度は各角度の固有長 $r\dd\phi$ を静止時計の $\sqrt f\dd t$ で割って $v=r\Omega/\sqrt f=\sqrt{m/(r-2m)}$ となる。$r=6m$ では $v=1/2$ である。この速度は座標の $r\Omega$ と異なり、時計の率の補正を含む。
:::
::::

## 確認テスト

::::quiz
Birkhoff の定理を適用できる領域を選べ。

- [x] 脈動する球対称星の外側の球対称真空領域
- [ ] 一般の時間依存する非球対称な真空領域
- [ ] 圧力が非零の球対称流体内部

:::explanation
定理は球対称性と真空の両方を用いる。脈動していても星の外がこれらを満たせば局所的に Schwarzschild 形となる。
:::
::::

::::quiz{answer="6" tol="1e-9"}
最内安定円軌道の $r/m$ を求めよ。

:::explanation
安定性が変わるのは $V_1''=2m(r-6m)/(r^3(r-3m))$ が零になる $r=6m$ である。
:::
::::

::::quiz
時間的円軌道の安定性の判定に用いるものを選べ。

- [x] $L$ を固定した $V_1''$ の符号
- [ ] $f$ の符号だけ
- [ ] $L(r)$ を先に代入した $V_1''$

:::explanation
動径摂動の前後で保存量 $L$ は同じである。従って固定した $L$ に対する $V_1$ の極小は安定、極大は不安定である。
:::
::::

::::quiz{answer="27" tol="1e-9"}
臨界衝突径数の二乗 $b^2/m^2$ を求めよ。

:::explanation
$r=3m$ で $V_0=L^2/(27m^2)$ が $E^2$ と等しいので $b^2/m^2=27$ である。
:::
::::

::::quiz{answer="0.5" tol="1e-7"}
$r_0=8m$ の時間的円軌道において、微小振幅の動径振動数と方位角振動数の比 $\Omega_r/\Omega_\phi$ を答えよ。

:::explanation
$\Omega_r^2/\Omega_\phi^2=1-6m/r_0$ に $r_0=8m$ を代入すると $1/4$ となる。正の振動数どうしの比はその平方根 $1/2$ である。
:::
::::

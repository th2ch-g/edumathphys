---
title: "Schwarzschild ブラックホール"
description: "正則な零座標で Schwarzschild 解を拡張し、地平線と特異点を因果関係と潮汐測定から区別する。"
notation:
  - symbol: "m"
    meaning: "GM"
  - symbol: "r_*"
    meaning: "亀座標"
  - symbol: "u,v"
    meaning: "外向きと内向きの零座標"
  - symbol: "U,V"
    meaning: "Kruskal–Szekeres 座標"
---

## この章の目標

- 曲率不変量と座標の退化を区別する。
- Eddington–Finkelstein 座標で未来地平線を横切る。
- Kruskal–Szekeres 拡張と Penrose 図の各領域を説明する。
- 橋が通過できない理由と自由落下中の潮汐を計算する。

## Schwarzschild 半径と座標特異点

:eqref[eq-gr-schwarzschild] で $f=1-2m/r$ とすると $g_{rr}=f^{-1}$ は $r=2m$ で発散する。この面積半径 $2m$ を :term[Schwarzschild 半径]{en="Schwarzschild radius" yomi="しゅわるつしるとはんけい"} と呼ぶ。ここでは $m>0$ とする。計量成分の発散だけでは、曲率そのものが発散するか判断できない。外部で $\theta^0=\sqrt f\dd t$、$\theta^1=\dd r/\sqrt f$、$\theta^2=r\dd\theta$、$\theta^3=r\sin\theta\dd\phi$ の枠を取る。:eqref[eq-gr-cartan] から $\omega^0{}_1=(f'/2)\dd t$、$\omega^2{}_1=\sqrt f\dd\theta$、$\omega^3{}_1=\sqrt f\sin\theta\dd\phi$、$\omega^3{}_2=\cos\theta\dd\phi$ である。微分と外積を取ると対空間で非零の対角成分は $R_{0101}=f''/2$、$R_{0i0i}=f'/2r$、$R_{1i1i}=-f'/2r$（$i=2,3$）、$R_{2323}=(1-f)/r^2$ になる。例えば $\dd\omega^0{}_1=(f''/2)\dd r\wedge\dd t$ が最初の成分を与え、$\omega^0{}_1\wedge\omega^1{}_2$ が時間・角成分を与える。二つの反対称対による四倍を数えると :term[Kretschmann スカラー]{en="Kretschmann scalar" yomi="くれっちまんすからー"} は

$$
R_{abcd}R^{abcd}=f''{}^2+4(f'/r)^2+4((1-f)/r^2)^2=\frac{48m^2}{r^6}.
\label{eq-gr-schwarz-kretschmann}
$$

$r=2m$ で $3/(4m^4)$ と有限であり、$r=0$ では無限大になる。前者での正則性は次の変換で示す。後者の発散は滑らかな非退化計量での延長を妨げる。

### 固有距離で見た外部の端

静的な切片の $r=2m$ から $r=R>2m$ までの距離は $\ell=\int_{2m}^R \sqrt{r/(r-2m)}\dd r$ である。$r=2m\cosh^2 s$ と置くと $\dd\ell=4m\cosh^2s\,\dd s$ なので、
$$
\ell=\sqrt{R(R-2m)}+2m\operatorname{arcosh}\sqrt{R/(2m)} .
$$
$R=2m+\epsilon$ では $\ell\simeq2\sqrt{2m\epsilon}$ である。従って $g_{rr}$ の発散は固有距離を無限にはしない。また $\sqrt{1-2m/R}\simeq\ell/(4m)$ なので、同じ近傍で静止時計の固有時率は固有距離に比例する。$g_{rr}$ が無限、lapse が零、曲率スカラーが有限という三つの性質を同じ点の極限で確認できる。

### 発散の場所を計量の成分によらず比べる

:::example[地平線と内部で曲率スカラーを比べる]
:eqref[eq-gr-schwarz-kretschmann] を $r=2m$ に入れると $3/(4m^4)$、$r=m$ では $48/m^4$ である。従って地平線を越えて半径が半分になったところで、同じ不変量は $64$ 倍になる。$r=2m$ の $g_{rr}$ の発散はこの値に反映しないが、$r$ をさらに $r=0$ に近づけるとこのスカラー自体が発散する。正則な計量と二階まで連続な微分を持つ点には、その発散する系列をそのまま収束させられない。
:::

### 静止を保つための力と座標の限界

外部で半径を一定に保つ観測者の四元速度は $u=f^{-1/2}\partial_t$ である。$\Gamma^r{}_{tt}=ff'/2$ から四元加速度の動径成分は $a^r=f'/2$、固有加速度の大きさは $m/(r^2\sqrt f)$ である。この加速度は $r\downarrow2m$ で発散する。一方、静止時計と無限遠の時計の率の比 $\sqrt f$ を掛ければ $m/r^2\to1/(4m)$ で有限となる。後に地平線の Killing 生成子から得るのと同じ値が、静止を保つ力を遠方の時刻で規格化した極限にも出る。

喉までの固有距離 $\ell$ を用いて、$r-2m=\ell^2/(8m)+\order{\ell^4/m^3}$ と展開する。すると動径計量は
$$
\dd s^2=-\frac{\ell^2}{16m^2}[1+\order{(\ell/m)^2}]\dd t^2+\dd\ell^2
$$
となる。角方向も小領域に限れば、平坦時空の加速座標を近似として使える。動径部の先頭で $X=\ell\cosh(t/(4m))$、$T=\ell\sinh(t/(4m))$ とすると、$-\dd T^2+\dd X^2=\dd\ell^2-\ell^2\dd t^2/(16m^2)$ である。従って固有時率が零に近づくことだけでは、そこで時空の曲率が発散するとは言えない。曲率を求めるには捨てた高次項や球面の計量も要る。

:::example[地平線近傍の二つの静止時計]
$\ell$ と $2\ell$ の位置に静止時計を置き、$\ell\ll m$ とする。無限遠での同じ座標時間の間に進む固有時は、それぞれ $\ell/(4m)$ と $2\ell/(4m)$ に比例する。外側の時計は約二倍進む。一方、静止に必要な固有加速度は先頭で $\ell^{-1}$ と $(2\ell)^{-1}$ である。地平線の接近時に「何が発散し何が零になるか」は、測定者の運動と時間の規格化まで指定して比較する。
:::

## Eddington–Finkelstein 座標と事象の地平線

$\dd r_*=\dd r/f$ を積分し、$r_*=r+2m\log\abs{r/(2m)-1}$ とする。これは :term[亀座標]{en="tortoise coordinate" yomi="かめざひょう"} である。入射零座標 $v=t+r_*$ に替えると

$$
\dd s^2=-f\dd v^2+2\dd v\dd r+r^2\dd\Omega^2 .
\label{eq-gr-ingoing-ef}
$$

:term[Eddington–Finkelstein 座標]{en="Eddington–Finkelstein coordinates" yomi="えでぃんとんふぃんけるすたいんざひょう"} の二次元ブロックの行列式は $-1$ であり、$r=2m$ でも非退化である。内向き光は $v$ 一定、外向き光は $\dd r/\dd v=f/2$ である。未来向き因果曲線で $\dd v>0$ なら、角方向の正の項を捨てても $\dd r/\dd v\le f/2$ となる。$v$ 一定の未来向き光でも $r$ は下がる。従って未来側の $r<2m$ から $r>2m$ に出られず、外部からは外向き光が無限遠に到達する。未来零無限遠を $\mathscr I^+$ と書き、そこへ信号を伝えられる領域の境界を :term[事象の地平線]{en="event horizon" yomi="じしょうのちへいせん"} と呼ぶ。この拡張では $r=2m$ が未来の境界である。例えば $r=m$ で外向き光にも $\dd r/\dd v=-1/2$ が成り立つ。

### 自由落下時刻の座標

$\dd T=\dd t+\sqrt{2m/r}\,\dd r/f$ と置く。$\dd t$ を :eqref[eq-gr-schwarzschild] に入れると、$\dd r^2$ の係数は $(1-2m/r)/f=1$ となる。動径部を平方完成すれば

$$
\dd s^2=-\dd T^2+(\dd r+\sqrt{2m/r}\,\dd T)^2+r^2\dd\Omega^2
$$

である。これが :term[Painlevé–Gullstrand 座標]{en="Painlevé–Gullstrand coordinates" yomi="ぱんるべぐるすとらんざひょう"} である。$T$ 一定面の誘導計量は Euclid の球座標形で、時空の曲率は交差項や時間の埋め込みに現れる。動径の零曲線は $\dd r/\dd T=-\sqrt{2m/r}\pm1$ である。例えば $r=m/2$ では外向きの枝でも傾きは $-1$ である。$E=1$ の落下では $\dot T=(1-2m/r)/f=1$ となるので、$T$ は各落下時計の固有時と定数だけ異なる。

:::example[地平線を横切る座標で接続を求める]
$\theta\in(0,\pi),r>0$ とし、:eqref[eq-gr-ingoing-ef] の $f$ を一旦任意の滑らかな関数として扱う。 座標の順序を $(v,r,\theta,\phi)$ と固定する。成分行列と逆行列は次の組になる。

$$
\begin{aligned}
(g_{ab})&=\left[\begin{matrix}- f{\left(r \right)} & 1 & 0 & 0\\1 & 0 & 0 & 0\\0 & 0 & r^{2} & 0\\0 & 0 & 0 & r^{2} \sin^{2}{\left(\theta \right)}\end{matrix}\right]\\
(g^{ab})&=\left[\begin{matrix}0 & 1 & 0 & 0\\1 & f{\left(r \right)} & 0 & 0\\0 & 0 & \frac{1}{r^{2}} & 0\\0 & 0 & 0 & \frac{1}{r^{2} \sin^{2}{\left(\theta \right)}}\end{matrix}\right]\\
\det g&=- r^{4} \sin^{2}{\left(\theta \right)}
\end{aligned}
$$

逆行列を掛けて単位行列を得ることを確かめてから、:eqref[eq-gr-christoffel] を使う。同じ下添字の順序を二度計算せず、$b\le c$ の組を選ぶ。中央の列には、座標に依存する計量成分の微分だけを残した。最後の列は微分を実行した値であり、表に出ない成分は添字交換で得られるものを除いて零である。

| 成分 | 微分への代入 | 結果 |
|---|---|---|
| $\Gamma^{v}_{v v}$ | $\frac{1}{2}\left(-\partial_{r}\left(- f{\left(r \right)}\right)\right)$ | $\frac{f'{\left(r \right)}}{2}$ |
| $\Gamma^{v}_{\theta \theta}$ | $\frac{1}{2}\left(-\partial_{r}\left(r^{2}\right)\right)$ | $- r$ |
| $\Gamma^{v}_{\phi \phi}$ | $\frac{1}{2}\left(-\partial_{r}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $- r \sin^{2}{\left(\theta \right)}$ |
| $\Gamma^{r}_{v v}$ | $\frac{f{\left(r \right)}}{2}\left(-\partial_{r}\left(- f{\left(r \right)}\right)\right)$ | $\frac{f{\left(r \right)} f'{\left(r \right)}}{2}$ |
| $\Gamma^{r}_{v r}$ | $\frac{1}{2}\left(\partial_{r}\left(- f{\left(r \right)}\right)\right)$ | $- \frac{f'{\left(r \right)}}{2}$ |
| $\Gamma^{r}_{\theta \theta}$ | $\frac{f{\left(r \right)}}{2}\left(-\partial_{r}\left(r^{2}\right)\right)$ | $- r f{\left(r \right)}$ |
| $\Gamma^{r}_{\phi \phi}$ | $\frac{f{\left(r \right)}}{2}\left(-\partial_{r}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $- r f{\left(r \right)} \sin^{2}{\left(\theta \right)}$ |
| $\Gamma^{\theta}_{r \theta}$ | $\frac{1}{2 r^{2}}\left(\partial_{r}\left(r^{2}\right)\right)$ | $\frac{1}{r}$ |
| $\Gamma^{\theta}_{\phi \phi}$ | $\frac{1}{2 r^{2}}\left(-\partial_{\theta}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $- \frac{\sin{\left(2 \theta \right)}}{2}$ |
| $\Gamma^{\phi}_{r \phi}$ | $\frac{1}{2 r^{2} \sin^{2}{\left(\theta \right)}}\left(\partial_{r}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{1}{r}$ |
| $\Gamma^{\phi}_{\theta \phi}$ | $\frac{1}{2 r^{2} \sin^{2}{\left(\theta \right)}}\left(\partial_{\theta}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{1}{\tan{\left(\theta \right)}}$ |

曲率の計算では :eqref[eq-gr-curvature-components] の微分差と積差を分ける。第二列は $\partial_c\Gamma^a{}_{db}$、第三列は $-\partial_d\Gamma^a{}_{cb}$ であり、第四列は接続の積を符号付きで足したものである。三列を足して最後の列を得る。

| 成分 | 第一の微分 | 第二の微分 | 接続の積 | 結果 |
|---|---|---|---|---|
| $R^{v}_{\theta v \theta}$ | $0$ | $0$ | $- \frac{r f'{\left(r \right)}}{2}$ | $- \frac{r f'{\left(r \right)}}{2}$ |
| $R^{v}_{\phi v \phi}$ | $0$ | $0$ | $- \frac{r \sin^{2}{\left(\theta \right)} f'{\left(r \right)}}{2}$ | $- \frac{r \sin^{2}{\left(\theta \right)} f'{\left(r \right)}}{2}$ |
| $R^{r}_{r v r}$ | $0$ | $\frac{f''{\left(r \right)}}{2}$ | $0$ | $\frac{f''{\left(r \right)}}{2}$ |
| $R^{r}_{\theta r \theta}$ | $- r f'{\left(r \right)} - f{\left(r \right)}$ | $0$ | $\frac{r f'{\left(r \right)}}{2} + f{\left(r \right)}$ | $- \frac{r f'{\left(r \right)}}{2}$ |
| $R^{r}_{\phi r \phi}$ | $- r \sin^{2}{\left(\theta \right)} f'{\left(r \right)} - f{\left(r \right)} \sin^{2}{\left(\theta \right)}$ | $0$ | $\frac{r \sin^{2}{\left(\theta \right)} f'{\left(r \right)}}{2} + f{\left(r \right)} \sin^{2}{\left(\theta \right)}$ | $- \frac{r \sin^{2}{\left(\theta \right)} f'{\left(r \right)}}{2}$ |
| $R^{\theta}_{\phi \theta \phi}$ | $- \cos{\left(2 \theta \right)}$ | $0$ | $- f{\left(r \right)} \sin^{2}{\left(\theta \right)} + \frac{\sin{\left(2 \theta \right)}}{2 \tan{\left(\theta \right)}}$ | $- \left(f{\left(r \right)} - 1\right) \sin^{2}{\left(\theta \right)}$ |

第一添字を計量で下ろすと $R_{abcd}=g_{ak}R^k{}_{bcd}$ となる。対の反対称性と交換対称性を用いて、非零の独立な組を下表に示す。

| 成分 | 結果 |
|---|---|
| $R_{v r v r}$ | $\frac{f''{\left(r \right)}}{2}$ |
| $R_{v \theta v \theta}$ | $\frac{r f{\left(r \right)} f'{\left(r \right)}}{2}$ |
| $R_{v \theta r \theta}$ | $- \frac{r f'{\left(r \right)}}{2}$ |
| $R_{v \phi v \phi}$ | $\frac{r f{\left(r \right)} \sin^{2}{\left(\theta \right)} f'{\left(r \right)}}{2}$ |
| $R_{v \phi r \phi}$ | $- \frac{r \sin^{2}{\left(\theta \right)} f'{\left(r \right)}}{2}$ |
| $R_{\theta \phi \theta \phi}$ | $- r^{2} \left(f{\left(r \right)} - 1\right) \sin^{2}{\left(\theta \right)}$ |

Ricci テンソルは第一・第三添字を縮約して得る。曲率が残る場合も相殺が起こるので、各項を落とさずに加える。

| 成分 | 縮約の項 | 代入 | 結果 |
|---|---|---|---|
| $R_{v v}$ | $R^{r}_{v r v}+R^{\theta}_{v \theta v}+R^{\phi}_{v \phi v}$ | $\frac{f{\left(r \right)} f''{\left(r \right)}}{2} + \frac{f{\left(r \right)} f'{\left(r \right)}}{2 r} + \frac{f{\left(r \right)} f'{\left(r \right)}}{2 r}$ | $\frac{\left(r f''{\left(r \right)} + 2 f'{\left(r \right)}\right) f{\left(r \right)}}{2 r}$ |
| $R_{v r}$ | $R^{v}_{v v r}+R^{\theta}_{v \theta r}+R^{\phi}_{v \phi r}$ | $- \frac{f''{\left(r \right)}}{2} - \frac{f'{\left(r \right)}}{2 r} - \frac{f'{\left(r \right)}}{2 r}$ | $- \frac{r f''{\left(r \right)} + 2 f'{\left(r \right)}}{2 r}$ |
| $R_{\theta \theta}$ | $R^{v}_{\theta v \theta}+R^{r}_{\theta r \theta}+R^{\phi}_{\theta \phi \theta}$ | $- \frac{r f'{\left(r \right)}}{2} - \frac{r f'{\left(r \right)}}{2} - \left(f{\left(r \right)} - 1\right)$ | $- r f'{\left(r \right)} - f{\left(r \right)} + 1$ |
| $R_{\phi \phi}$ | $R^{v}_{\phi v \phi}+R^{r}_{\phi r \phi}+R^{\theta}_{\phi \theta \phi}$ | $- \frac{r \sin^{2}{\left(\theta \right)} f'{\left(r \right)}}{2} - \frac{r \sin^{2}{\left(\theta \right)} f'{\left(r \right)}}{2} - \left(f{\left(r \right)} - 1\right) \sin^{2}{\left(\theta \right)}$ | $- \left(r f'{\left(r \right)} + f{\left(r \right)} - 1\right) \sin^{2}{\left(\theta \right)}$ |

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
| $G_{v v}$ | $\frac{\left(r f''{\left(r \right)} + 2 f'{\left(r \right)}\right) f{\left(r \right)}}{2 r}$ | $\frac{\left(r^{2} f''{\left(r \right)} + 4 r f'{\left(r \right)} + 2 f{\left(r \right)} - 2\right) f{\left(r \right)}}{2 r^{2}}$ | $- \frac{\left(r f'{\left(r \right)} + f{\left(r \right)} - 1\right) f{\left(r \right)}}{r^{2}}$ |
| $G_{v r}$ | $- \frac{r f''{\left(r \right)} + 2 f'{\left(r \right)}}{2 r}$ | $- \frac{r^{2} f''{\left(r \right)} + 4 r f'{\left(r \right)} + 2 f{\left(r \right)} - 2}{2 r^{2}}$ | $\frac{r f'{\left(r \right)} + f{\left(r \right)} - 1}{r^{2}}$ |
| $G_{\theta \theta}$ | $- r f'{\left(r \right)} - f{\left(r \right)} + 1$ | $- \frac{r^{2} f''{\left(r \right)}}{2} - 2 r f'{\left(r \right)} - f{\left(r \right)} + 1$ | $\frac{r \left(r f''{\left(r \right)} + 2 f'{\left(r \right)}\right)}{2}$ |
| $G_{\phi \phi}$ | $- \left(r f'{\left(r \right)} + f{\left(r \right)} - 1\right) \sin^{2}{\left(\theta \right)}$ | $- \frac{\left(r^{2} f''{\left(r \right)} + 4 r f'{\left(r \right)} + 2 f{\left(r \right)} - 2\right) \sin^{2}{\left(\theta \right)}}{2}$ | $\frac{r \left(r f''{\left(r \right)} + 2 f'{\left(r \right)}\right) \sin^{2}{\left(\theta \right)}}{2}$ |

アフィン径数の微分を点で表し、$w_a$ を座標 $a$ の速度成分とする。ここでは添字を下げる操作を行っていない。測地線の右辺 $-\Gamma^a{}_{bc}w_bw_c$ では $b\ne c$ の項を二度数える。

$$
\begin{aligned}
\ddot{v}&=\frac{2 r w_{\phi}^{2} \sin^{2}{\left(\theta \right)} + 2 r w_{\theta}^{2} - w_{v}^{2} f'{\left(r \right)}}{2}\\
\ddot{r}&=\frac{2 r w_{\phi}^{2} f{\left(r \right)} \sin^{2}{\left(\theta \right)} + 2 r w_{\theta}^{2} f{\left(r \right)} + 2 w_{r} w_{v} f'{\left(r \right)} - w_{v}^{2} f{\left(r \right)} f'{\left(r \right)}}{2}\\
\ddot{\theta}&=\frac{r w_{\phi}^{2} \sin{\left(2 \theta \right)} - 4 w_{\theta} w_{r}}{2 r}\\
\ddot{\phi}&=- \frac{2 w_{\phi} \left(r w_{\theta} + w_{r} \tan{\left(\theta \right)}\right)}{r \tan{\left(\theta \right)}}
\end{aligned}
$$

また計量ノルムは

$$
g(\dot x,\dot x)=r^{2} w_{\phi}^{2} \sin^{2}{\left(\theta \right)} + r^{2} w_{\theta}^{2} + w_{r} w_{v} + w_{v} \left(w_{r} - w_{v} f{\left(r \right)}\right)
$$
である。これを径数で微分して上の加速度を代入すると、係数の微分と速度の微分の和が零になる。従って初期のノルムが保存される。

逆計量に $1/f$ は現れない。従って $f$ が零であっても計量と接続は正則であり、測地線の初期値問題をその点から解ける。角方向の速度が零なら $\dot v=0$ の枝は $\ddot r=0$ である。未来向きは $\dot r<0$ で、有限の径数で地平線を横切る。
:::

### Eddington–Finkelstein 座標で外向き光の解を積分する

外部で $dv/dr=2/f$ の枝を取ると、二半径 $r_1,r_2>2m$ の間の入射時刻の差は
$$
v_2-v_1=2(r_2-r_1)+4m\log\frac{r_2-2m}{r_1-2m}.
$$
これは座標 $v$ での外向き光の式で、外向き光に沿う遅れ時刻 $u=v-2r_*$ は一定である。$r_1\downarrow2m$ を固定した $v_1$ で取れば、外側の定数半径へ到達する $v_2$ は対数的に遅くなる。一方、入射光は $v$ 一定であり、地平線も横切れる。どちらの式も $v$ が正則な座標であることと整合する。

:::example[地平線近くの赤方偏移と遅れ]
明るさの比較を除き、二本の外向き光を同じ入射零面 $v=v_1$ で $r_1=2m+\epsilon$ と $2m+2\epsilon$ から発射する。固定した受信半径 $r_2$ での到着時刻の差は $2\epsilon+4m\log2$ であり、小さい半径から発した光が遅れる。$\epsilon$ を小さくしてもこの差は $4m\log2$ に近づく。一本ごとの遅れが $-\log\epsilon$ を含むので、半径の差より比に敏感になる。
:::

### 二つの零膨張率から球面を分類する

入射座標で球面に垂直な二つの未来零ベクトルを、$\ell=\partial_v+(f/2)\partial_r$、$n=-\partial_r$ と取る。$\ell\cdot n=-1$ であり、各方向への面積 $A=4\pi r^2$ の相対微分を膨張率と呼ぶ。計算すると
$$
\theta_{(\ell)}=\frac{\ell(A)}A=\frac fr,\qquad
\theta_{(n)}=\frac{n(A)}A=-\frac2r .
$$
$r>2m$ では外向きに面積が増え、内向きには減る。$r<2m$ では両方で減るため、未来 :term[捕捉面]{en="trapped surface" yomi="ほそくめん"} である。地平線の球では外向きの膨張率が零になる。この一致は静的な真空解での結果である。[相対論的星](/physics/general-relativity/relativistic-stars/) の崩壊時空では、膨張率が示す局所の境界と、外へ届く信号から定義する境界を別々に追える。

:::example[同じ半径を過去側の座標で調べる]
外向きの零座標 $u$ では $\dd s^2=-f\dd u^2-2\dd u\dd r+r^2\dd\Omega^2$ である。右外部から延長した未来零方向を $\partial_r$ と $\partial_u-(f/2)\partial_r$ と選べば、膨張率は $2/r$ と $-f/r$ である。従って $r<2m$ で両方が正になる。入射座標の内側では両方が負だったため、同じ $r$ でも異なる時空領域を覆うことが分かる。この違いが Kruskal–Szekeres 拡張での未来と過去の内側の区別に対応する。
:::

### 落下時計を時刻にして動径の交差項を消す

Painlevé–Gullstrand の時刻 $T$ を使い、$E=1$ で落ちる動径測地線にラベル $\varrho$ を付ける。ラベルを
$$
r^{3/2}=\frac32\sqrt{2m}(\varrho-T)
$$
で定めれば、$\varrho$ 一定で $\partial_T r=-\sqrt{2m/r}$ になる。$\dd r=\sqrt{2m/r}(\dd\varrho-\dd T)$ を時刻の計量に入れると
$$
\dd s^2=-\dd T^2+\frac{2m}{r(T,\varrho)}\dd\varrho^2+r(T,\varrho)^2\dd\Omega^2 .
$$
これが :term[Lemaître 座標]{en="Lemaître coordinates" yomi="るめーとるざひょう"} である。$\varrho$ を固定した時計が測る $T$ は固有時であり、$r=2m$ での動径空間係数は $1$ である。地平線は $\varrho-T=4m/3$ と表される。光の動径速度 $\dd\varrho/\dd T=\pm\sqrt{r/(2m)}$ の正の枝は、この直線に沿って進む。

:::example[同じ落下ラベルで固有時を確認する]
$r=2m$ で $\varrho-T_H=4m/3$、$r=m/2$ で $\varrho-T_1=m/6$ である。従って $T_1-T_H=7m/6$ であり、同じ世界線に沿った固有時を与える。$\varrho$ が異なる時計との同時性はこの落下系の $T$ で定めている。外部の静止系の同時面とは一致しない。
:::

## Kruskal–Szekeres 座標と最大拡張

外部で $u=t-r_*$、$v=t+r_*$ とすると動径計量は $-f\dd u\dd v$ である。$U=-e^{-u/(4m)}$、$V=e^{v/(4m)}$ と置くと

$$
-UV=(r/(2m)-1)e^{r/(2m)},\qquad
\dd s^2=-\frac{32m^3}{r}e^{-r/(2m)}\dd U\dd V+r^2\dd\Omega^2.
\label{eq-gr-kruskal}
$$

右辺の $r$ の関数は $r>0$ で導関数 $r e^{r/(2m)}/(4m^2)>0$、値域 $(-1,\infty)$ である。したがって $UV<1$ に一意に滑らかな $r$ を定める。これが :term[Kruskal–Szekeres 座標]{en="Kruskal–Szekeres coordinates" yomi="くるすかるせけれしゅざひょう"} である。未来を $U,V$ が増える向きに取る。領域 I は $U<0<V$、領域 III は $V<0<U$ でともに外部、領域 II は $U,V>0$ で未来ブラックホール、領域 IV は $U,V<0$ で過去ホワイトホールである。$U=V=0$ は面積半径 $2m$ の二分岐球面である。$UV=1$ の二枝のうち未来側のものは未来因果曲線の終端となる。$r>0$ の全域を含み、特異境界を曲率が発散しない点として付け加えることはできない。星の崩壊でできる時空には、一般に領域 III と IV を備えた永遠解の構造は現れない。

### 指数変換の係数と四領域

外部で $U=-e^{-u/(4m)},V=e^{v/(4m)}$ と取ると、$-UV=e^{r_*/(2m)}=(r/(2m)-1)e^{r/(2m)}$ になる。$y=r/(2m)>0$ とすると右辺の導関数は $ye^y>0$ であるため、左辺の値から正の $r$ が一意に決まる。値域は $(-1,\infty)$ であり、従って Lorentz 領域は $UV<1$、$r=0$ は $UV=1$ である。地平線は $UV=0$ だが、$U=0,V>0$ と $V=0,U<0$ は未来と過去で異なる。未来向き動径因果曲線では $U,V$ とも非減少であり、例えば $U>0,V>0$ から $U<0,V>0$ の右外部へは到達できない。一方、$U<0,V<0$ の白色穴領域からはどちらの外部にも未来に到達できる。

:::example[特異境界に到達する座標方向]
$U=1/2$ の未来向き動径光線では $V$ が増える。$V=0$ で左側の地平線を横切ると $r=2m$ で、$V=2$ で $UV=1$ に近づくと $r\to0$ である。途中で $r$ が小さくなることは $y\mapsto(y-1)e^y$ の単調性から分かる。有限の座標値 $V=2$ であることだけでは測地線完備性は決まらず、アフィン径数は別に積分する。
:::

### 地平線の Killing 接線をアフィン化する

未来右地平線で $U=0,V>0$ とする。無限遠で規格化した時間並進は $\xi=(-U\partial_U+V\partial_V)/(4m)$ なので、地平線上では $\xi=(V/(4m))\partial_V$ である。動径計量の係数は $UV$ にだけ依存するため、$U=0$ の生成子上では $V$ で微分しても変化しない。角速度を零としたとき $\Gamma^V{}_{VV}$ はその係数の対数微分であり、地平線上で零になる。従って $k=\partial_V$ はアフィン接線である。一方、$\nabla_\xi\xi=(1/(4m))\xi$ であり、Killing 時刻はアフィンではない。例えば $V$ を二倍にする区間は Killing 時刻で $4m\log2$ と常に同じだが、アフィン長は初めの $V$ に比例する。

### 光的方向でのアフィン径数と境界

動径計量を $-\mathcal F(r)\dd U\dd V$ と書く。$U$ 一定の零線に沿って $\Gamma^V{}_{VV}=\partial_V\log\mathcal F$ なので、アフィン径数 $\lambda$ には $\mathcal F\dd V/\dd\lambda=\mathrm{const.}$ が成り立つ。$U\ne0$ なら関係 $-UV=(r/(2m)-1)e^{r/(2m)}$ を微分して
$$
\dd V=-\frac{r\,e^{r/(2m)}}{4m^2U}\dd r,\qquad
\mathcal F\dd V=-\frac{8m}{U}\dd r .
$$
従って $r$ 自身がアフィン径数の一次関数である。$U<0$ を固定して $V\to\infty$ に向かう外向き光のアフィン長は無限になる。これに対し $U>0$ を固定して未来内部で $V\to1/U$ へ向かう光は有限のアフィン径数で終わる。有限な座標の範囲へ図示することと、有限のアフィン径数で境界に近づくことは別々に判定できる。$U=0$ の地平線では先の結果のとおり $V$ 自身がアフィン径数であり、未来へ無限に進める。$V$ が一定の零線でも二座標を交換した式が成り立つ。

### 準備した拡張の外にどのような端があり得るか

座標パッチを拡張できるか確かめるには、曲がった座標の端に達した因果曲線の振る舞いを調べる。球面方向の座標特異点には複数の角座標を用いる。$U,V$ が有限で $UV<1$ なら、既に計量が滑らかであり、$U=0$ や $V=0$ でも測地線の局所解を延長できる。$UV\to1$ では曲率が発散する。無限遠について、時間的または零の測地線に対し $\varepsilon=1,0$ と置けば
$$
\dot r^2=E^2-f(\varepsilon+L^2/r^2).
$$
$r$ が大きい所では右辺は有界であるため、$r\to\infty$ に進むには無限の固有時またはアフィン径数が要る。$r$ が正の有限値で保たれるのに $|U|$ または $|V|$ が発散する場合には、$t\to\pm\infty$ に対応するか、$r\to2m$ で座標の指数の片側だけが発散する。前者では $f$ が零から離れていれば $\dot t=E/f$ が有界である。後者では $E\ne0$ の横断測地線について、横断の向きに応じて $\dot v=(E+\dot r)/f$ または $\dot u=(E-\dot r)/f$ が有界になり、一方の正則な座標で有限の点に到達する。$E=0$ の時間的測地線が内側から $r=2m$ に近づく場合は $t$ が一定で、$U,V$ はともに零に近づく。零測地線で地平線に留まる場合には既に求めたアフィン径数を使う。

従ってこの延長には、有限アフィン径数で達する正則な因果測地線の端が残っていない。滑らかな時空をさらに延長できるなら、延長領域の境界を横切る有限径数の時間的測地線が局所正規近傍で取れるので、そのような端を少なくとも一つ必要とする。因果測地線の端が $r=0$ にある場合は曲率不変量の発散が連続な曲率を妨げる。この意味で Kruskal–Szekeres 時空は滑らかな計量に関して最大である。曲率を定義できない低い正則性での延長問題には、この論証をそのまま使えない。

### 拡張の全体を横切る初期面

$T=(U+V)/2,X=(V-U)/2$ では $UV=T^2-X^2<1$ である。従って $T=0$ 面はすべての $X\in\R$ を含む。因果的な曲線では $T$ が狭義に増加し、角方向の速度を含めた場合にも $|\dd X/\dd T|\le1$ である。$T$ の有限な区間では $X$ に有限の極限がある。もし $r$ の極限も正なら、球面上の速度は因果性から $\sqrt{\mathcal F}/r$ 以下に抑えられるため、角方向にも極限を持つ。このとき曲線は正則な計量の点から先へ延ばせる。

$T>0$ の領域で未来向きのまま $T=0$ から遠ざかる場合と、過去向きに $T=0$ に向かう場合を区別する。後者が $T=0$ に達するより前に $r=0$ に近づくには、未来内部で $UV\to1$ となる必要がある。しかし未来内部を過去向きに進むと $U,V$ とも非増加なので、両方が正の間は積が減る。従ってその端へは近づけない。$T<0$ の曲線を未来向きに延ばす場合も、過去内部 $U,V<0$ で積が減るので同じである。よって任意の延長不能な因果曲線は $T=0$ に一度交わる。二つの外部と喉を含むこの空間切片は、この意味で :term[Cauchy 面]{en="Cauchy surface" yomi="こーしーめん"} である。この面のデータを一方の外部だけのデータに置き換えると、他方の外部や過去側から来る因果曲線を取り逃す。

:::example[一方の外部の静止時刻の切片だけを使う場合]
右外部で $t=0$ を取ると $U=-V<0$ の半直線だけが得られる。左外部の一定半径の世界線では $U>0,V<0$ がずっと成り立つため、この半直線に一度も交わらない。半直線の端に二分岐球面を足しても、左外部の世界線を捕まえることはできない。初期面を拡張全体に取る場合には $T=0$ のもう半分も必要である。
:::

## 共形図（Penrose 図）

:term[Penrose 図]{en="Penrose diagram" yomi="ぺんろーずず"} は球面方向を省略し、零方向を保ったまま無限遠までを有限の図に描く。$\widehat U=\arctan U$、$\widehat V=\arctan V$ とすると動径計量は正の関数と $-\dd\widehat U\dd\widehat V$ の積である。横軸を $(\widehat V-\widehat U)/2$、縦軸を $(\widehat V+\widehat U)/2$ とすると、光は傾き $\pm1$ を持つ。未来内部で $UV=1$ なら $\widehat U+\widehat V=\pi/2$ なので特異点は水平線となる。過去内部の対応線は $-\pi/2$ である。右外部の $\mathcal I^+$ は $\widehat V=\pi/2$、$-\pi/2<\widehat U<0$、$\mathcal I^-$ は $\widehat U=-\pi/2$、$0<\widehat V<\pi/2$ となる。右外部の空間無限遠 $i^0$ は $(\widehat U,\widehat V)=(-\pi/2,\pi/2)$ である。このコンパクト化は因果方向を表すが、図の線分の長さから固有時や面積を読めない。

### 共形図で傾きを保持する写像

$U,V$ に対して $\bar U=\arctan U,\bar V=\arctan V$ と置き、$\bar T=(\bar V+\bar U)/2,\bar X=(\bar V-\bar U)/2$ とする。動径計量は正の係数を掛けた $-\dd\bar T^2+\dd\bar X^2$ となる。従って $U$ または $V$ が一定の光線は傾き $\pm1$ の直線に写る。未来特異境界では $U,V>0,UV=1$ なので $\bar U+\bar V=\pi/2$、すなわち $\bar T=\pi/4$ である。右の未来光的無限遠は $V\to\infty$ かつ $U<0$ を固定する境界であり、$r\to\infty$ となる。同じ $V$ の極限でも $U$ が正であれば $UV<1$ を満たせないので、無限遠として加える辺の範囲は座標の値域だけからは決められない。

### 四つの端と二つの外部を一枚に描く

$\widehat T=(\widehat V+\widehat U)/2$、$\widehat X=(\widehat V-\widehat U)/2$ を軸に取る。右外部の未来時間的無限遠 $i^+$ は $(\widehat U,\widehat V)=(0,\pi/2)$ で、過去側の $i^-$ は $(-\pi/2,0)$ である。左外部では $\widehat X$ を反転して読める。未来と過去の特異境界の線分は、それぞれ左右の $i^+$、左右の $i^-$ を結ぶ。これは図上での接続を表すもので、端点を滑らかな物理時空の点として含めるものではない。

```tikz caption="最大拡張した Schwarzschild 時空の共形図。太い水平線は曲率特異境界、破線は地平線を表す。"
\begin{tikzpicture}[scale=1.25]
  \draw[very thick] (-1,1)--(1,1) (-1,-1)--(1,-1);
  \draw (-1,1)--(-2,0)--(-1,-1) (1,1)--(2,0)--(1,-1);
  \draw[dashed] (-1,-1)--(1,1) (-1,1)--(1,-1);
  \node at (1.05,0) {I};
  \node at (0,0.62) {II};
  \node at (-1.05,0) {III};
  \node at (0,-0.62) {IV};
  \node[above] at (0,1) {$r=0$};
  \node[below] at (0,-1) {$r=0$};
  \node[above right] at (1,1) {$i^+$};
  \node[below right] at (1,-1) {$i^-$};
  \node[right] at (2,0) {$i^0$};
  \node[right] at (1.52,0.55) {$\mathcal I^+$};
  \node[right] at (1.52,-0.55) {$\mathcal I^-$};
  \node[above left] at (-1,1) {$i^+$};
  \node[below left] at (-1,-1) {$i^-$};
  \node[left] at (-2,0) {$i^0$};
\end{tikzpicture}
```

図では各内点が球面を表す。従って図上の未来向き因果曲線の傾きは $\pm1$ の間でも、上に省略した角速度があれば図上ではその分だけ時間的になる。右の $\mathscr I^+$ を選んでその因果的過去を取ると $U<0$ の領域が含まれる。左の $\mathscr I^+$ まで含めれば、$V<0$ も含まれる。二端の和を使うと、未来ブラックホール領域は $U,V\ge0$ の側にあり、過去ホワイトホール領域からは両端のいずれへも出られる。

:::example[二つの無限遠を選んで光を追う]
$U=-1,V=-1/2$ の IV の点を取る。$U$ を固定して $V$ を増やすと $V=0$ で右外部に入り、さらに右の $\mathscr I^+$ に達する。$V$ を固定して $U$ を増やすと左外部を経由して左の $\mathscr I^+$ に達する。一方、I の $U=-1,V=1/2$ から左の $\mathscr I^+$ を目指すと、途中で $V$ を負にする必要があるので不可能である。
:::

## Einstein–Rosen 橋とワームホール

Kruskal の時間 $T=(V+U)/2$ を零とする面では $V=-U$ なので $UV=-U^2$ であり、$r$ は二分岐球面で最小になる。この空間面の動径計量は、両側で $\dd r^2/f$ と書ける。喉の近くでは $\ell=\pm\int_{2m}^r\dd s/\sqrt{f(s)}\simeq\pm2\sqrt{2m(r-2m)}$ であり、$r=2m+\ell^2/(8m)+\order{\ell^4}$ と滑らかである。二つの外部を接続するこの空間幾何を :term[Einstein–Rosen 橋]{en="Einstein–Rosen bridge" yomi="あいんしゅたいんろーぜんきょう"} と呼ぶ。ただし領域 I から III に移るには $V>0$ から $V<0$ へ移動する必要がある。未来因果曲線では $\dot U,\dot V\ge0$ なので不可能である。空間切片が連結であることと、時空を光や粒子が横断できることは別の条件である。

### 面積の最小と因果的な行き来

時刻反転で固定される空間切片で、喉を $\ell=0$ とする。外部の切片における $r(\ell)$ は $(dr/d\ell)^2=1-2m/r$ を満たす。喉の近くで $r=2m+\ell^2/(8m)+\order{\ell^4}$ なので、球の面積は $\ell=0$ で局所最小となる。切片の二つの端で $r\to\infty$ となるが、この切片上の道は空間的である。因果的な行き来を調べるには、右外部から左外部を結ぶ曲線で $U,V$ が非減少であるかを確認すればよい。右外部は $V>0$、左外部は $V<0$ なので右から左の未来因果曲線はない。逆向きには $U$ が減少しなければならないため、左から右も同じ理由で結べない。

### 喉の面積の変化と横断条件

:::example[二つの端の中間に現れる最小球]
空間切片の固有距離 $\ell$ を喉で零とし、$r=2m+\ell^2/(8m)+\cdots$ を面積 $A=4\pi r^2$ に代入する。$A=16\pi m^2+2\pi\ell^2+\cdots$ なので $A'(0)=0$、$A''(0)=4\pi>0$ である。半径 $r$ を座標にすると $A'=8\pi r$ は零にならないが、この座標は喉で左右を識別できない。従って面積の極値を検査する際は、極値の候補をまたいで滑らかな座標 $\ell$ を使う。
:::

:::example[横断できる喉を持つ計量と必要な応力]
Einstein–Rosen 橋との比較に、$\dd s^2=-\dd t^2+\dd\ell^2+(\ell^2+L^2)\dd\Omega^2$、$-\infty<\ell<\infty,L>0$ を用いる。最小の球は $\ell=0$ で面積が $4\pi L^2$ である。時間項が $-1$ なので、動径の零曲線は $\dd\ell/\dd t=\pm1$ であり、有限な $t$ で喉を横切る。次にこの幾何を要求する源を場の方程式から読む。 座標の順序を $(t,\ell,\theta,\phi)$ と固定する。成分行列と逆行列は次の組になる。

$$
\begin{aligned}
(g_{ab})&=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & L^{2} + \ell^{2} & 0\\0 & 0 & 0 & \left(L^{2} + \ell^{2}\right) \sin^{2}{\left(\theta \right)}\end{matrix}\right]\\
(g^{ab})&=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & \frac{1}{L^{2} + \ell^{2}} & 0\\0 & 0 & 0 & \frac{1}{\left(L^{2} + \ell^{2}\right) \sin^{2}{\left(\theta \right)}}\end{matrix}\right]\\
\det g&=- \left(L^{2} + \ell^{2}\right)^{2} \sin^{2}{\left(\theta \right)}
\end{aligned}
$$

逆行列を掛けて単位行列を得ることを確かめてから、:eqref[eq-gr-christoffel] を使う。同じ下添字の順序を二度計算せず、$b\le c$ の組を選ぶ。中央の列には、座標に依存する計量成分の微分だけを残した。最後の列は微分を実行した値であり、表に出ない成分は添字交換で得られるものを除いて零である。

| 成分 | 微分への代入 | 結果 |
|---|---|---|
| $\Gamma^{\ell}_{\theta \theta}$ | $\frac{1}{2}\left(-\partial_{\ell}\left(L^{2} + \ell^{2}\right)\right)$ | $- \ell$ |
| $\Gamma^{\ell}_{\phi \phi}$ | $\frac{1}{2}\left(-\partial_{\ell}\left(\left(L^{2} + \ell^{2}\right) \sin^{2}{\left(\theta \right)}\right)\right)$ | $- \ell \sin^{2}{\left(\theta \right)}$ |
| $\Gamma^{\theta}_{\ell \theta}$ | $\frac{1}{2 \left(L^{2} + \ell^{2}\right)}\left(\partial_{\ell}\left(L^{2} + \ell^{2}\right)\right)$ | $\frac{\ell}{L^{2} + \ell^{2}}$ |
| $\Gamma^{\theta}_{\phi \phi}$ | $\frac{1}{2 \left(L^{2} + \ell^{2}\right)}\left(-\partial_{\theta}\left(\left(L^{2} + \ell^{2}\right) \sin^{2}{\left(\theta \right)}\right)\right)$ | $- \frac{\sin{\left(2 \theta \right)}}{2}$ |
| $\Gamma^{\phi}_{\ell \phi}$ | $\frac{1}{2 \left(L^{2} + \ell^{2}\right) \sin^{2}{\left(\theta \right)}}\left(\partial_{\ell}\left(\left(L^{2} + \ell^{2}\right) \sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{\ell}{L^{2} + \ell^{2}}$ |
| $\Gamma^{\phi}_{\theta \phi}$ | $\frac{1}{2 \left(L^{2} + \ell^{2}\right) \sin^{2}{\left(\theta \right)}}\left(\partial_{\theta}\left(\left(L^{2} + \ell^{2}\right) \sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{1}{\tan{\left(\theta \right)}}$ |

曲率の計算では :eqref[eq-gr-curvature-components] の微分差と積差を分ける。第二列は $\partial_c\Gamma^a{}_{db}$、第三列は $-\partial_d\Gamma^a{}_{cb}$ であり、第四列は接続の積を符号付きで足したものである。三列を足して最後の列を得る。

| 成分 | 第一の微分 | 第二の微分 | 接続の積 | 結果 |
|---|---|---|---|---|
| $R^{\ell}_{\theta \ell \theta}$ | $-1$ | $0$ | $\frac{\ell^{2}}{L^{2} + \ell^{2}}$ | $- \frac{L^{2}}{L^{2} + \ell^{2}}$ |
| $R^{\ell}_{\phi \ell \phi}$ | $- \sin^{2}{\left(\theta \right)}$ | $0$ | $\frac{\ell^{2} \sin^{2}{\left(\theta \right)}}{L^{2} + \ell^{2}}$ | $- \frac{L^{2} \sin^{2}{\left(\theta \right)}}{L^{2} + \ell^{2}}$ |
| $R^{\theta}_{\phi \theta \phi}$ | $- \cos{\left(2 \theta \right)}$ | $0$ | $- \frac{\ell^{2} \sin^{2}{\left(\theta \right)}}{L^{2} + \ell^{2}} + \frac{\sin{\left(2 \theta \right)}}{2 \tan{\left(\theta \right)}}$ | $\frac{L^{2} \sin^{2}{\left(\theta \right)}}{L^{2} + \ell^{2}}$ |

第一添字を計量で下ろすと $R_{abcd}=g_{ak}R^k{}_{bcd}$ となる。対の反対称性と交換対称性を用いて、非零の独立な組を下表に示す。

| 成分 | 結果 |
|---|---|
| $R_{\ell \theta \ell \theta}$ | $- \frac{L^{2}}{L^{2} + \ell^{2}}$ |
| $R_{\ell \phi \ell \phi}$ | $- \frac{L^{2} \sin^{2}{\left(\theta \right)}}{L^{2} + \ell^{2}}$ |
| $R_{\theta \phi \theta \phi}$ | $L^{2} \sin^{2}{\left(\theta \right)}$ |

Ricci テンソルは第一・第三添字を縮約して得る。曲率が残る場合も相殺が起こるので、各項を落とさずに加える。

| 成分 | 縮約の項 | 代入 | 結果 |
|---|---|---|---|
| $R_{\ell \ell}$ | $R^{\theta}_{\ell \theta \ell}+R^{\phi}_{\ell \phi \ell}$ | $- \frac{L^{2}}{\left(L^{2} + \ell^{2}\right)^{2}} - \frac{L^{2}}{\left(L^{2} + \ell^{2}\right)^{2}}$ | $- \frac{2 L^{2}}{\left(L^{2} + \ell^{2}\right)^{2}}$ |
| $R_{\theta \theta}$ | $R^{\ell}_{\theta \ell \theta}+R^{\phi}_{\theta \phi \theta}$ | $- \frac{L^{2}}{L^{2} + \ell^{2}} + \frac{L^{2}}{L^{2} + \ell^{2}}$ | $0$ |
| $R_{\phi \phi}$ | $R^{\ell}_{\phi \ell \phi}+R^{\theta}_{\phi \theta \phi}$ | $- \frac{L^{2} \sin^{2}{\left(\theta \right)}}{L^{2} + \ell^{2}} + \frac{L^{2} \sin^{2}{\left(\theta \right)}}{L^{2} + \ell^{2}}$ | $0$ |

逆計量と縮約すると Ricci スカラーとなる。さらに $G_{ab}=R_{ab}-g_{ab}R/2$ を成分ごとに計算する。

$$
\begin{aligned}
R&=\sum_{a,b}g^{ab}R_{ab}\\
&=- \frac{2 L^{2}}{\left(L^{2} + \ell^{2}\right)^{2}}\\
&=- \frac{2 L^{2}}{\left(L^{2} + \ell^{2}\right)^{2}}
\end{aligned}
$$

非零の Einstein テンソル成分は次のとおりである。

| 成分 | Ricci 項 | 跡の項 | 差 |
|---|---|---|---|
| $G_{t t}$ | $0$ | $\frac{L^{2}}{\left(L^{2} + \ell^{2}\right)^{2}}$ | $- \frac{L^{2}}{\left(L^{2} + \ell^{2}\right)^{2}}$ |
| $G_{\ell \ell}$ | $- \frac{2 L^{2}}{\left(L^{2} + \ell^{2}\right)^{2}}$ | $- \frac{L^{2}}{\left(L^{2} + \ell^{2}\right)^{2}}$ | $- \frac{L^{2}}{\left(L^{2} + \ell^{2}\right)^{2}}$ |
| $G_{\theta \theta}$ | $0$ | $- \frac{L^{2}}{L^{2} + \ell^{2}}$ | $\frac{L^{2}}{L^{2} + \ell^{2}}$ |
| $G_{\phi \phi}$ | $0$ | $- \frac{L^{2} \sin^{2}{\left(\theta \right)}}{L^{2} + \ell^{2}}$ | $\frac{L^{2} \sin^{2}{\left(\theta \right)}}{L^{2} + \ell^{2}}$ |

光的ベクトル $k=\partial_t+\partial_\ell$ を取れば $g(k,k)=0$ である。表から $R_{ab}k^ak^b=-2L^2/(\ell^2+L^2)^2<0$ なので、:eqref[eq-gr-einstein-field] の光的縮約による $8\pi G T_{ab}k^ak^b$ も負である。よってこの計量を保つ物質は零エネルギー条件に反する。喉で $(\dd r/\dd\ell)=0$ かつ $(\dd^2r/\dd\ell^2)>0$ であるため、喉を越える光束は面積の減少から増加に変わる。計算した負の光的収束はこの変化と一致する。
:::

### 喉をもつ空間切片の断面曲率

$\dd\ell^2+r(\ell)^2\dd\Omega^2$ の三次元計量をとる。動径と角方向の断面曲率は $K_{\ell\theta}=-r''/r$、角方向だけの断面曲率は $K_{\theta\phi}=(1-r'^2)/r^2$ である。実際、$\Gamma^\ell{}_{\theta\theta}=-rr'$ と $\Gamma^\theta{}_{\ell\theta}=r'/r$ から $R_{\ell\theta\ell\theta}=-rr''$ となり、球面固有の曲率には動径接続の寄与 $-r'^2/r^2$ が加わる。橋では $r'^2=1-2m/r$ を微分して $r''=m/r^2$ となる。従って動径を含む断面の曲率は $-m/r^3$、角断面は $2m/r^3$ である。三次元のスカラー曲率は各断面曲率の和の二倍なので零になる。橋の空間はこの縮約が零でも各断面には曲率がある。

赤道だけを取り出し Euclid の回転面 $\dd r^2+\dd Z^2+r^2\dd\phi^2$ に等長に埋め込む。内在する計量の動径係数を一致させるには $(\dd Z/\dd r)^2=(1/f)-1=2m/(r-2m)$ が必要である。積分すると
$$
Z=\pm2\sqrt{2m(r-2m)},\qquad r=2m+\frac{Z^2}{8m}.
$$
この回転面は喉の最小半径と二つの側を描く。しかしここで選んだのは時間反転で固定される一枚の空間切片であり、$Z$ が時空の時間座標になるわけではない。

:::example[喉と少し離れた球面の断面曲率]
喉 $r=2m$ では $K_{\ell\theta}=-1/(8m^2)$、$K_{\theta\phi}=1/(4m^2)$ である。$r=4m$ の球では各々が $1/8$ になる。埋め込みの高さはこの半径で $|Z|=4m$ だが、喉からの内在的距離は $\sqrt8m+2m\operatorname{arcosh}\sqrt2$ で、約 $4.59m$ である。高さと内在的距離は異なるため、距離の計算には誘導計量を使う。
:::

## 地平線への落下と潮汐効果

無限遠で静止から落ちる動径粒子は $E=1,L=0$ を持つので $\dot r=-\sqrt{2m/r}$ である。地平線から特異点までの固有時は $2(2m)^{3/2}/(3\sqrt{2m})=4m/3$、SI では $4GM/(3c^3)$ である。静的座標の $\dot t=1/f$ は発散するが、入射座標では $\dot v=(1-\sqrt{2m/r})/f=1/(1+\sqrt{2m/r})$ であり、地平線で $1/2$ となる。動径自由落下枠は外部静止枠を $01$ 面で boost したものである。$R_{0101}$ はその平面の面積二形式が不変なので変わらず、角方向では $R_{0i0i}=-R_{1i1i}$ が boost の $\cosh^2-\sinh^2=1$ を残す。従って測地線偏差の固有値は $(2,-1,-1)m/r^3$ である。地平線で半径方向の係数は $1/(4m^2)$ で、質量が大きいほど小さくなる。例えば $M=10M_\odot$ の地平線では $1\,\mathrm m$ 離れた自由落下点の半径方向相対加速度が約 $10^8\,\mathrm{m\,s^{-2}}$ であるが、$10^6M_\odot$ なら約 $0.01\,\mathrm{m\,s^{-2}}$ である。

### 地平線の内側にいられる固有時

未来内部では $r$ が時間的で、全ての未来時間的曲線に沿って減る。$x=2m/r-1>0$ とすると動径・時間部は $x\dd t^2-\dd r^2/x$ である。未来曲線の $\dd\tau^2$ は $\dd r^2/x$ から空間方向の寄与を引いたものなので、地平線から特異点までの固有時は
$$
\tau\le \int_0^{2m}\sqrt{\frac{r}{2m-r}}\dd r=\pi m .
$$
最後の積分では $r=2m\sin^2\chi$ と置けば被積分が $4m\sin^2\chi\dd\chi$ になる。等号は内部で $t,\theta,\phi$ が一定となる曲線であり、永遠解では二分岐球面から出る動径測地線に対応する。通常の外部から未来地平線を横切る落下では、この値は上界である。例えば無限遠から静止で落下した $E=1$ の場合は $4m/3$ であり、最大値に対する比は $4/(3\pi)\simeq0.424$ となる。半径の減少は時間的方向に進むため、推進を加えても内部に任意の固有時間とどまれない。

### 平行枠での潮汐を静止枠から得る

自由落下観測者の速度は局所的に静止枠の $0,1$ 平面の boost で表せる。$e'_0=\gamma(e_0+ve_1)$、$e'_1=\gamma(ve_0+e_1)$ とすると $e'_0\wedge e'_1=e_0\wedge e_1$ なので $R_{1'0'1'0'}=R_{1010}$ である。横方向 $i=2,3$ では $R_{i0i1}=0$、$R_{i1i1}=-R_{i0i0}$ を使い $R_{i0'i0'}=\gamma^2(1-v^2)R_{i0i0}=R_{i0i0}$ である。従って潮汐の固有値は動径落下の瞬間の速度に依存しない。静止枠は $r=2m$ までしか使えないが、落下枠での値は $m/r^3$ の倍として有限なので、正則な座標で連続に延長できる。

:::example[二つの質量で地平線の潮汐を比べる]
地平線の動径偏差係数は $2m/(2m)^3=1/(4m^2)$ である。SI の時間を秒で取れば $c^2/(4m^2)$ となる。$m$ が $10^6$ 倍の穴では、この係数は $10^{-12}$ 倍なので、同じ長さの装置が横切るときの相対加速度も $10^{-12}$ 倍となる。地平線を越えた後に $r$ が減れば係数は $r^{-3}$ で増加する。
:::

:::example[落下の時刻と潮汐の増え方]
$E=1,L=0$ の動径落下では $\dd r/\dd\tau=-\sqrt{2m/r}$ なので、$r$ から $0$ までの残りの固有時は $2r^{3/2}/(3\sqrt{2m})$ である。動径の潮汐係数は $2m/r^3$ である。$x=r/(2m)$ と置き、時刻と曲率をそれぞれ $m$ と $m^{-2}$ で無次元化する。

ここでは$x$ は面積半径を地平線半径で割った値であるとする。計算に用いる式を先に並べる。

$$
\begin{aligned}
\Delta\tau/m&=\frac{4 x^{\frac{3}{2}}}{3}\\
m^2 \abs{R_{1010}}&=\frac{1}{4 x^{3}}\\
m^2 \abs{R_{2020}}&=\frac{1}{8 x^{3}}
\end{aligned}
$$

次の条件を比較する。入力を区別してから計算し、定数の規格化をそろえた値を表に示す。

| $x$ | $\Delta\tau/m$ | $m^2 \abs{R_{1010}}$ | $m^2 \abs{R_{2020}}$ |
|---|---|---|---|
| $4$ | $\frac{32}{3}$ | $\frac{1}{256}$ | $\frac{1}{512}$ |
| $1$ | $\frac{4}{3}$ | $\frac{1}{4}$ | $\frac{1}{8}$ |
| $\frac{1}{4}$ | $\frac{1}{6}$ | $16$ | $8$ |
| $\frac{1}{16}$ | $\frac{1}{48}$ | $1024$ | $512$ |

地平線から $r=0$ までの固有時は $4m/3$ である。例えば $r$ が $2m$ の四分の一になると残りの固有時は八分の一、潮汐は六十四倍になる。二つの式から半径を消すと、動径潮汐係数は $4/[9(\Delta\tau)^2]$ となる。
:::

### 落下中に見える光の振動数

光の無限遠での保存エネルギーを $\mathcal E=-k_t>0$ とする。外部の動径光では $k^r=\sigma\mathcal E$、$\sigma=+1$ を外向き、$\sigma=-1$ を内向きとする。外部の動径落下者は $u^t=E/f$、$u^r=-\sqrt{E^2-f}$ である。局所角振動数は $-k\cdot u$ なので
$$
\frac{\omega_{\rm ff}}{\mathcal E}
=\frac{E+\sigma\sqrt{E^2-f}}f .
$$
内向きの光では $1/(E+\sqrt{E^2-f})$ に等しく、地平線で $1/(2E)$ に近づく。無限遠から静止で落ちる人なら、外から入る動径光を有限の振動数で測る。静止観測者の光の振動数 $\mathcal E/\sqrt f$ の発散は、同じ光でも落下観測者へ移る boost と相殺する。一方、外向き光を固定した $\mathcal E$ で比較すると発散する。その外向き光を地平線直前の落下装置が有限の固有振動数で放ったとすれば、$\mathcal E$ が零に近づくため、遠方の測定では大きく赤方偏移する。

:::example[入射光と出射光の比]
$E=1$ の落下者に対して $r=8m$ では $\sqrt{1-f}=1/2$ である。同じ $\mathcal E$ の動径光は内向きで $2\mathcal E/3$、外向きで $2\mathcal E$ になる。従って二方向の周波数比は $3$ である。地平線へ近づいた極限で比は大きくなるが、ここでは同じ保存エネルギーの二方向の光を比較した。
:::

### 落下枠で分離ベクトルを終端まで積分する

$E=1$ の動径粒子で $\delta=\tau_{\rm sing}-\tau$ とおくと $r^3=(9/2)m\delta^2$ である。従って半径方向の Jacobi 方程式は $\dd^2\xi_r/\dd\delta^2=4\xi_r/(9\delta^2)$、角方向は $\dd^2\xi_\perp/\dd\delta^2=-2\xi_\perp/(9\delta^2)$ となる。$\xi=\delta^\alpha$ を入れると $\alpha(\alpha-1)$ が右辺の係数に等しい。ゆえに
$$
\xi_r=A_r\delta^{4/3}+B_r\delta^{-1/3},\qquad
\xi_\perp=A_\perp\delta^{2/3}+B_\perp\delta^{1/3}.
$$
角方向はどちらの独立解も縮む。動径では一般に伸びるモードの係数 $B_r$ が非零である。ただし測地線偏差が扱うのは無限小分離なので、有限サイズの物体の変形を終端までこの一次式で読み取れない。

:::example[初めに相対速度が零の分離]
$\delta=\delta_0$ で $\xi_r=\xi_0$、$\dd\xi_r/\dd\tau=0$ とすると、$A_r\delta_0^{4/3}=\xi_0/5$、$B_r\delta_0^{-1/3}=4\xi_0/5$ である。例えば $\delta=\delta_0/8$ では $\xi_r/\xi_0=(1/5)(1/16)+(4/5)2=129/80$ である。同じ初期条件の横方向では $A_\perp\delta_0^{2/3}=-\xi_0$、$B_\perp\delta_0^{1/3}=2\xi_0$ となり、$\xi_\perp/\xi_0=-1/4+1=3/4$ になる。二つの符号が異なる潮汐の効果を、時刻を有限に進めた比で確かめられる。
:::

### 光子球と落下観測者の脱出円錐

外部の静止枠で光のエネルギーを $\omega_{\rm stat}$、外向き動径方向との角を $\alpha$ とする。保存量は $\mathcal E=\sqrt f\,\omega_{\rm stat}$、$|L|=r\,\omega_{\rm stat}\sin\alpha$ なので、衝突径数は $b=r\sin\alpha/\sqrt f$ である。零測地線の動径式は $(\dot r/\mathcal E)^2=1-b^2 f/r^2$ である。$f/r^2$ の極大は $r=3m$ にあり、値は $1/(27m^2)$ である。従って $2m<r<3m$ から外向きに放つ光は $b<3\sqrt3m$ の場合に無限遠へ出る。等号で外向きに放つと、極大に漸近して有限のアフィン径数ではその先へ進まない。内向きの光には内側に転回点がないため落ちる。終端での符号を調べるため、指定した半径がポテンシャル極大の内側か外側かが要る。$r>3m$ で外向きに出た光はすべて脱出する。$r>3m$ から内向きに出た光は $b>3\sqrt3m$ なら極大までの間で転回し、それより小さければ内側へ進む。

内側から放射する静止観測者の脱出円錐の境界を $\alpha_c$ とすると $\sin\alpha_c=3\sqrt3m\sqrt f/r$ である。この角は地平線近くで零になる。次に同じ点で $E=1$ の動径落下者を考える。静止枠に対する落下速度を $-v$、$v=\sqrt{2m/r}$ とすると、光の動径成分の Lorentz 変換から
$$
\cos\alpha_{\rm ff}=\frac{\cos\alpha+v}{1+v\cos\alpha},\qquad
1-\cos\alpha_{\rm ff,c}
=\frac{(1-v)(1-\cos\alpha_c)}{1+v\cos\alpha_c}
$$
となる。落下観測者の脱出円錐は静止観測者よりも狭い。地平線近くでは $1-v\simeq f/2$、$1-\cos\alpha_c\simeq27f/8$ であるため、$1-\cos\alpha_{\rm ff,c}\simeq27f^2/32$ となる。球面の正規直交枠で等方に放った光のうち、この円錐に入る割合は $(1-\cos\alpha_{\rm ff,c})/2$ である。これは光子の脱出割合を幾何から数えたものであり、遠方でのエネルギー流束には一光子当たりの赤方偏移も入る。

:::example[光子球の内側で衝突径数を指定する]
$r=5m/2$ から外向きに光を放つ。$f=1/5$ なので境界の $\sin\alpha_c=6\sqrt{15}/25\simeq0.930$、$\cos\alpha_c=\sqrt{85}/25$ である。例えば半径方向から $30^\circ$ の光は $\sin\alpha=1/2$ なので脱出する。$80^\circ$ では $\sin\alpha\simeq0.985$ で境界より大きいため、外向きであっても途中で転回する。後者の衝突径数は $b=5\sqrt5m\sin(80^\circ)/2$ で、$3\sqrt3m$ を越える。地平線の外にいることは、外へ出られる方向が少なくとも一つあることを意味する。
:::

### 収縮の大きさを落下する局所枠で測る

動径から外れた時間的測地線を、正則な入射座標の定数 $E,L$ で追う。赤道で $\dot r<0$ とすれば $f\dot v-\dot r=E$、$\dot\phi=L/r^2$、$\dot r^2=E^2-f(1+L^2/r^2)$ である。地平線で $E>0$ を固定すると $\dot r\to-E$ であり、
$$
\dot v=\frac{1+L^2/r^2}{E+\sqrt{E^2-f(1+L^2/r^2)}}
\ \longrightarrow\ \frac{1+L^2/(4m^2)}{2E}.
$$
すべて有限だが、球面方向を含む落下枠では動径方向だけの boost を使った潮汐の議論をそのまま適用できない。外部で定めた衝突径数や単位質量当たりの角運動量は地平線を越えても保存されるが、局所で測る速度の成分は観測者の枠に依存する。

:::example[角運動量を持つ粒子の地平線通過]
$E=1,L=2m$ とし、外から落ちる枝を選ぶ。このとき $\dot r^2=2m(r^2-2mr+4m^2)/r^3>0$ であり、正の半径に転回点がない。$r=2m$ で $\dot r=-1$、$\dot v=1$、$\dot\phi=1/(2m)$ である。入射計量に代入すると $2\dot v\dot r+r^2\dot\phi^2=-2+1=-1$ となる。時間的な正規化を保ったまま横断することを、特異な $t$ を経由せず確認できる。
:::

:::example[縦向きの落下時刻を SI へ戻す]
$M=10M_\odot$ で $GM_\odot/c^3=4.93\,\mathrm{\mu s}$ を使うと、$E=1$ の軌道で地平線から特異境界までの固有時は約 $65.7\,\mathrm{\mu s}$ である。$10^6M_\odot$ なら約 $6.57\,\mathrm s$ になる。同じ半径比で軌道を比べたとき、固有時は質量に比例するが、地平線の潮汐係数は質量の二乗に反比例する。静止座標の時刻が発散することから、これらの固有時の大きさは得られない。
:::

## 演習問題

::::exercise{level=basic}
$r=4m$ の外向き光の $\dd r/\dd v$ を求めよ。

:::solution
$f=1/2$ なので :eqref[eq-gr-ingoing-ef] の零条件から $\dd r/\dd v=f/2=1/4$ である。
:::
::::

::::exercise{level=basic}
$r=2m$ で $m^4R_{abcd}R^{abcd}$ を求めよ。

:::solution
:eqref[eq-gr-schwarz-kretschmann] から $48/2^6=3/4$ である。
:::
::::

::::exercise{level=basic}
$U=1/2,V=-1$ の Kruskal 点の領域を答えよ。

:::solution
$UV=-1/2<0$ で外部であり、$U>0>V$ なので左外部 III である。
:::
::::

::::exercise{level=basic}
$E=1,L=0$ の粒子の $r=8m$ での $\dot r$ を求めよ。

:::solution
$\dot r=-\sqrt{2m/(8m)}=-1/2$ である。
:::
::::

::::exercise{level=advanced}
未来地平線の法線が光的であることを示せ。

:::solution
入射座標の動径ブロックの逆は $\begin{pmatrix}0&1\\1&f\end{pmatrix}$ なので、面 $r=2m$ の法線 $\dd r$ のノルムは $g^{rr}=f=0$ である。
:::
::::

::::exercise{level=advanced}
外向きの零座標 $u=t-r_*$ の計量と正則性を求めよ。

:::solution
$\dd t=\dd u+\dd r/f$ を代入すると $\dd s^2=-f\dd u^2-2\dd u\dd r+r^2\dd\Omega^2$ となる。動径ブロックの行列式は $-1$ であり、過去側の地平線を滑らかに覆う。
:::
::::

::::exercise{level=advanced}
Kruskal の時間平行移動場を求め、ノルムが外部の値と一致することを確かめよ。

:::solution
$\partial_t U=-U/(4m)$、$\partial_t V=V/(4m)$ なので $\partial_t=(-U\partial_U+V\partial_V)/(4m)$ である。:eqref[eq-gr-kruskal] に入れるとノルムは $(2m/r)e^{-r/(2m)}UV=-(1-2m/r)$ となる。
:::
::::

::::exercise{level=advanced}
$E=1$ の動径落下で $r=2m$ から $r=m/2$ までの固有時を求めよ。

:::solution
$2((2m)^{3/2}-(m/2)^{3/2})/(3\sqrt{2m})=(4/3-1/6)m=7m/6$ である。
:::
::::

::::exercise{level=advanced}
$r_0>2m$ で外部の静止点から放す時間的な動径測地線について、未来地平線上での $\dot v$ を $r_0,m$ で表せ。無限遠から放す極限と地平線に近い位置から放す極限を比較せよ。

:::solution
静止初期条件で $E=\sqrt{1-2m/r_0}>0$、内向きに $\dot r=-\sqrt{E^2-f}$ を選ぶ。入射時刻の微分は $\dot v=(E+\dot r)/f$ であり、分子と分母に $E-\dot r$ を掛ければ
$$
\dot v=\frac1{E+\sqrt{E^2-f}}
$$
となる。この式は $r=2m$ でも値が有限で、$\dot v_H=1/(2E)=[2\sqrt{1-2m/r_0}]^{-1}$ である。$r_0\to\infty$ で $1/2$ へ行くが、$r_0\downarrow2m$ ではこの値は発散する。各固定した $r_0>2m$ で座標が正則であることと、軌道の族の規格化した接線に一様な上界がないことは両立する。
:::
::::

## 確認テスト

::::quiz
地平線で成り立つものをすべて選べ。

- [x] 面 $r=2m$ は光的である。
- [x] 入射座標で計量は正則である。
- [ ] Kretschmann スカラーが発散する。

:::explanation
入射座標で $g^{rr}=f=0$ だがブロックの行列式は $-1$ である。したがって法線は光的で、曲率不変量 $48m^2/r^6$ も有限である。
:::
::::

::::quiz{answer="0.5" tol="1e-9"}
$E=1,L=0$ の入射粒子の地平線での $\dd v/\dd\tau$ を求めよ。

:::explanation
$\dot v=1/(1+\sqrt{2m/r})$ に $r=2m$ を代入すると $1/2$ である。
:::
::::

::::quiz
右外部から左外部へ未来向き信号が進めない理由を選べ。

- [x] その移動は $V$ の符号を正から負に変える。
- [ ] 空間面に喉が存在しない。
- [ ] 光は動径平面の零曲線を進めない。

:::explanation
未来向きの因果性から $V$ は非減少である。空間面自体は喉で接続しているが、右外部の $V>0$ を左外部の $V<0$ に変えられない。
:::
::::

::::quiz{answer="0.25" tol="1e-9"}
地平線での半径方向潮汐係数を $1/m^2$ 単位で答えよ。

:::explanation
半径方向の相対加速度係数は $2m/r^3$ で、$r=2m$ では $1/(4m^2)$ である。
:::
::::

::::quiz{answer="1.33333333333" tol="1e-7"}
$E=1,L=0$ の落下で $r=2m$ から $r=0$ の端までの固有時を $m$ 単位で求めよ。

:::explanation
$-\dot r=\sqrt{2m/r}$ なので $\int_0^{2m}\sqrt{r/(2m)}\dd r=2(2m)^{3/2}/(3\sqrt{2m})=4m/3$ である。
:::
::::

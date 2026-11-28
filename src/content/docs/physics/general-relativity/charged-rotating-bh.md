---
title: "荷電ブラックホールと回転ブラックホール"
description: "電荷と角運動量を持つ定常解の地平線を調べ、エネルギー抽出と Kerr 測地線の変数分離を導く。"
notation:
  - symbol: "m"
    meaning: "GM"
  - symbol: "a"
    meaning: "単位質量当たりの角運動量 J/M"
  - symbol: "q"
    meaning: "幾何単位の電荷長"
  - symbol: "\\Sigma,\\Delta"
    meaning: "回転解の計量関数"
  - symbol: "\\mathcal Q"
    meaning: "Carter 定数"
---

## この章の目標

- 電磁場の応力から球対称荷電解を求める。
- Kerr 解の地平線、エルゴ領域、環特異点を区別する。
- エネルギー抽出の符号と測地線の分離を計算する。
- 一意性定理の適用条件を説明する。

## Reissner–Nordström 解と内部地平線

静的球対称な電磁真空で、電荷以外の電流や磁気単極子を置かない。真空 Maxwell 方程式 $\partial_r(\sqrt{-g}F^{rt})=0$ に $g_{tt}=-e^{2\Phi}$、$g_{rr}=f^{-1}$ を入れると、$r^2e^{-\Phi}\sqrt f F_{tr}$ が定数となる。これを $-Q_0$ とすれば、正規直交枠で $F_{01}=-Q_0/r^2$ である。$T^t{}_t=T^r{}_r=-Q_0^2/(2r^4)$ なので :eqref[eq-gr-spherical-field] の差から $g_{tt}=-f$、$g_{rr}=1/f$ とできる。$q^2=4\pi G Q_0^2$ とし、$(rf)'=1-q^2/r^2$ を積分すると

$$
f=1-2m/r+q^2/r^2,\qquad r_\pm=m\pm\sqrt{m^2-q^2}.
\label{eq-gr-rn}
$$

:term[Reissner–Nordström 解]{en="Reissner–Nordström solution" yomi="らいすなーのるどしゅとれーむかい"} はこれを Schwarzschild 形に入れたものである。実際 $G^\theta{}_\theta=f''/2+f'/r=q^2/r^4$ であり、角圧も満たす。$0<q<m$ では内外二面が光的になる。ここで $q=\sqrt{4\pi G}\abs{Q_0}\ge0$ であり、電荷の符号は $Q_0$ が担う。入射座標の光の傾きは $\dd r/\dd v=f/2$ で、外では正、二面の間では負、内側では正である。理想的な解析拡張では $r=r_-$ が外側の初期面に関する Cauchy 地平線となる。$f'(r_-)<0$ なので、その近傍の外向き光は $r-r_-\propto e^{f'(r_-)v/2}$ に従い、同じ有限の $v$ 間隔が正則な座標では指数的に短い間隔へ写る。周波数としての青方偏移は、後で正則な座標を作って確認する。この青方偏移は内部摂動を小さいままと見なす仮定を脅かす。例えば $q=4m/5$ では $r_+=8m/5,r_-=2m/5$ である。

:::example[電場に対応する Ricci テンソルを成分で確認する]
$f=1-2m/r+q^2/r^2$ を静的な外部で用いる。$q^2$ は計量に現れる長さの二乗である。質量の項と電場の項が曲率の縮約で異なる働きをすることを確かめる。 座標の順序を $(t,r,\theta,\phi)$ と固定する。成分行列と逆行列は次の組になる。

$$
\begin{aligned}
(g_{ab})&=\left[\begin{matrix}\frac{2 m}{r} - \frac{q^{2}}{r^{2}} - 1 & 0 & 0 & 0\\0 & \frac{1}{- \frac{2 m}{r} + \frac{q^{2}}{r^{2}} + 1} & 0 & 0\\0 & 0 & r^{2} & 0\\0 & 0 & 0 & r^{2} \sin^{2}{\left(\theta \right)}\end{matrix}\right]\\
(g^{ab})&=\left[\begin{matrix}- \frac{r^{2}}{- 2 m r + q^{2} + r^{2}} & 0 & 0 & 0\\0 & \frac{- 2 m r + q^{2} + r^{2}}{r^{2}} & 0 & 0\\0 & 0 & \frac{1}{r^{2}} & 0\\0 & 0 & 0 & \frac{1}{r^{2} \sin^{2}{\left(\theta \right)}}\end{matrix}\right]\\
\det g&=- r^{4} \sin^{2}{\left(\theta \right)}
\end{aligned}
$$

逆行列を掛けて単位行列を得ることを確かめてから、:eqref[eq-gr-christoffel] を使う。同じ下添字の順序を二度計算せず、$b\le c$ の組を選ぶ。中央の列には、座標に依存する計量成分の微分だけを残した。最後の列は微分を実行した値であり、表に出ない成分は添字交換で得られるものを除いて零である。

| 成分 | 微分への代入 | 結果 |
|---|---|---|
| $\Gamma^{t}_{t r}$ | $- \frac{r^{2}}{2 \left(- 2 m r + q^{2} + r^{2}\right)}\left(\partial_{r}\left(\frac{2 m}{r} - \frac{q^{2}}{r^{2}} - 1\right)\right)$ | $- \frac{- m r + q^{2}}{r \left(- 2 m r + q^{2} + r^{2}\right)}$ |
| $\Gamma^{r}_{t t}$ | $\frac{- 2 m r + q^{2} + r^{2}}{2 r^{2}}\left(-\partial_{r}\left(\frac{2 m}{r} - \frac{q^{2}}{r^{2}} - 1\right)\right)$ | $- \frac{\left(- m r + q^{2}\right) \left(- 2 m r + q^{2} + r^{2}\right)}{r^{5}}$ |
| $\Gamma^{r}_{r r}$ | $\frac{- 2 m r + q^{2} + r^{2}}{2 r^{2}}\left(\partial_{r}\left(\frac{1}{- \frac{2 m}{r} + \frac{q^{2}}{r^{2}} + 1}\right)+\partial_{r}\left(\frac{1}{- \frac{2 m}{r} + \frac{q^{2}}{r^{2}} + 1}\right)-\partial_{r}\left(\frac{1}{- \frac{2 m}{r} + \frac{q^{2}}{r^{2}} + 1}\right)\right)$ | $\frac{- m r + q^{2}}{r \left(- 2 m r + q^{2} + r^{2}\right)}$ |
| $\Gamma^{r}_{\theta \theta}$ | $\frac{- 2 m r + q^{2} + r^{2}}{2 r^{2}}\left(-\partial_{r}\left(r^{2}\right)\right)$ | $- \frac{- 2 m r + q^{2} + r^{2}}{r}$ |
| $\Gamma^{r}_{\phi \phi}$ | $\frac{- 2 m r + q^{2} + r^{2}}{2 r^{2}}\left(-\partial_{r}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $- \frac{\left(- 2 m r + q^{2} + r^{2}\right) \sin^{2}{\left(\theta \right)}}{r}$ |
| $\Gamma^{\theta}_{r \theta}$ | $\frac{1}{2 r^{2}}\left(\partial_{r}\left(r^{2}\right)\right)$ | $\frac{1}{r}$ |
| $\Gamma^{\theta}_{\phi \phi}$ | $\frac{1}{2 r^{2}}\left(-\partial_{\theta}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $- \frac{\sin{\left(2 \theta \right)}}{2}$ |
| $\Gamma^{\phi}_{r \phi}$ | $\frac{1}{2 r^{2} \sin^{2}{\left(\theta \right)}}\left(\partial_{r}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{1}{r}$ |
| $\Gamma^{\phi}_{\theta \phi}$ | $\frac{1}{2 r^{2} \sin^{2}{\left(\theta \right)}}\left(\partial_{\theta}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{1}{\tan{\left(\theta \right)}}$ |

曲率の計算では :eqref[eq-gr-curvature-components] の微分差と積差を分ける。第二列は $\partial_c\Gamma^a{}_{db}$、第三列は $-\partial_d\Gamma^a{}_{cb}$ であり、第四列は接続の積を符号付きで足したものである。三列を足して最後の列を得る。

| 成分 | 第一の微分 | 第二の微分 | 接続の積 | 結果 |
|---|---|---|---|---|
| $R^{t}_{r t r}$ | $0$ | $- \frac{m}{r \left(- 2 m r + q^{2} + r^{2}\right)} + \frac{\left(2 m - 2 r\right) \left(- m r + q^{2}\right)}{r \left(- 2 m r + q^{2} + r^{2}\right)^{2}} - \frac{- m r + q^{2}}{r^{2} \left(- 2 m r + q^{2} + r^{2}\right)}$ | $- \frac{\left(- m r + q^{2}\right)^{2}}{r^{2} \left(- 2 m r + q^{2} + r^{2}\right)^{2}} - \frac{\left(- m r + q^{2}\right)^{2}}{r^{2} \left(- 2 m r + q^{2} + r^{2}\right)^{2}}$ | $- \frac{- 2 m r + 3 q^{2}}{r^{2} \left(- 2 m r + q^{2} + r^{2}\right)}$ |
| $R^{t}_{\theta t \theta}$ | $0$ | $0$ | $\frac{- m r + q^{2}}{r^{2}}$ | $\frac{- m r + q^{2}}{r^{2}}$ |
| $R^{t}_{\phi t \phi}$ | $0$ | $0$ | $\frac{\left(- m r + q^{2}\right) \sin^{2}{\left(\theta \right)}}{r^{2}}$ | $\frac{\left(- m r + q^{2}\right) \sin^{2}{\left(\theta \right)}}{r^{2}}$ |
| $R^{r}_{\theta r \theta}$ | $- \frac{- 2 m + 2 r}{r} + \frac{- 2 m r + q^{2} + r^{2}}{r^{2}}$ | $0$ | $- \frac{- m r + q^{2}}{r^{2}} + \frac{- 2 m r + q^{2} + r^{2}}{r^{2}}$ | $\frac{- m r + q^{2}}{r^{2}}$ |
| $R^{r}_{\phi r \phi}$ | $- \frac{\left(- 2 m + 2 r\right) \sin^{2}{\left(\theta \right)}}{r} + \frac{\left(- 2 m r + q^{2} + r^{2}\right) \sin^{2}{\left(\theta \right)}}{r^{2}}$ | $0$ | $- \frac{\left(- m r + q^{2}\right) \sin^{2}{\left(\theta \right)}}{r^{2}} + \frac{\left(- 2 m r + q^{2} + r^{2}\right) \sin^{2}{\left(\theta \right)}}{r^{2}}$ | $\frac{\left(- m r + q^{2}\right) \sin^{2}{\left(\theta \right)}}{r^{2}}$ |
| $R^{\theta}_{\phi \theta \phi}$ | $- \cos{\left(2 \theta \right)}$ | $0$ | $\frac{\sin{\left(2 \theta \right)}}{2 \tan{\left(\theta \right)}} - \frac{\left(- 2 m r + q^{2} + r^{2}\right) \sin^{2}{\left(\theta \right)}}{r^{2}}$ | $- \frac{\left(- 2 m r + q^{2}\right) \sin^{2}{\left(\theta \right)}}{r^{2}}$ |

第一添字を計量で下ろすと $R_{abcd}=g_{ak}R^k{}_{bcd}$ となる。対の反対称性と交換対称性を用いて、非零の独立な組を下表に示す。

| 成分 | 結果 |
|---|---|
| $R_{t r t r}$ | $\frac{- 2 m r + 3 q^{2}}{r^{4}}$ |
| $R_{t \theta t \theta}$ | $- \frac{\left(- m r + q^{2}\right) \left(- 2 m r + q^{2} + r^{2}\right)}{r^{4}}$ |
| $R_{t \phi t \phi}$ | $- \frac{\left(- m r + q^{2}\right) \left(- 2 m r + q^{2} + r^{2}\right) \sin^{2}{\left(\theta \right)}}{r^{4}}$ |
| $R_{r \theta r \theta}$ | $\frac{- m r + q^{2}}{- 2 m r + q^{2} + r^{2}}$ |
| $R_{r \phi r \phi}$ | $\frac{\left(- m r + q^{2}\right) \sin^{2}{\left(\theta \right)}}{- 2 m r + q^{2} + r^{2}}$ |
| $R_{\theta \phi \theta \phi}$ | $- \left(- 2 m r + q^{2}\right) \sin^{2}{\left(\theta \right)}$ |

Ricci テンソルは第一・第三添字を縮約して得る。曲率が残る場合も相殺が起こるので、各項を落とさずに加える。

| 成分 | 縮約の項 | 代入 | 結果 |
|---|---|---|---|
| $R_{t t}$ | $R^{r}_{t r t}+R^{\theta}_{t \theta t}+R^{\phi}_{t \phi t}$ | $\frac{\left(- 2 m r + 3 q^{2}\right) \left(- 2 m r + q^{2} + r^{2}\right)}{r^{6}} - \frac{\left(- m r + q^{2}\right) \left(- 2 m r + q^{2} + r^{2}\right)}{r^{6}} - \frac{\left(- m r + q^{2}\right) \left(- 2 m r + q^{2} + r^{2}\right)}{r^{6}}$ | $\frac{q^{2} \left(- 2 m r + q^{2} + r^{2}\right)}{r^{6}}$ |
| $R_{r r}$ | $R^{t}_{r t r}+R^{\theta}_{r \theta r}+R^{\phi}_{r \phi r}$ | $- \frac{- 2 m r + 3 q^{2}}{r^{2} \left(- 2 m r + q^{2} + r^{2}\right)} + \frac{- m r + q^{2}}{r^{2} \left(- 2 m r + q^{2} + r^{2}\right)} + \frac{- m r + q^{2}}{r^{2} \left(- 2 m r + q^{2} + r^{2}\right)}$ | $- \frac{q^{2}}{r^{2} \left(- 2 m r + q^{2} + r^{2}\right)}$ |
| $R_{\theta \theta}$ | $R^{t}_{\theta t \theta}+R^{r}_{\theta r \theta}+R^{\phi}_{\theta \phi \theta}$ | $- \frac{- 2 m r + q^{2}}{r^{2}} + \frac{- m r + q^{2}}{r^{2}} + \frac{- m r + q^{2}}{r^{2}}$ | $\frac{q^{2}}{r^{2}}$ |
| $R_{\phi \phi}$ | $R^{t}_{\phi t \phi}+R^{r}_{\phi r \phi}+R^{\theta}_{\phi \theta \phi}$ | $- \frac{\left(- 2 m r + q^{2}\right) \sin^{2}{\left(\theta \right)}}{r^{2}} + \frac{\left(- m r + q^{2}\right) \sin^{2}{\left(\theta \right)}}{r^{2}} + \frac{\left(- m r + q^{2}\right) \sin^{2}{\left(\theta \right)}}{r^{2}}$ | $\frac{q^{2} \sin^{2}{\left(\theta \right)}}{r^{2}}$ |

逆計量と縮約すると Ricci スカラーとなる。さらに $G_{ab}=R_{ab}-g_{ab}R/2$ を成分ごとに計算する。

$$
\begin{aligned}
R&=\sum_{a,b}g^{ab}R_{ab}\\
&=- \frac{q^{2}}{r^{4}} - \frac{q^{2}}{r^{4}} + \frac{q^{2}}{r^{4}} + \frac{q^{2}}{r^{4}}\\
&=0
\end{aligned}
$$

非零の Einstein テンソル成分は次のとおりである。

| 成分 | Ricci 項 | 跡の項 | 差 |
|---|---|---|---|
| $G_{t t}$ | $\frac{q^{2} \left(- 2 m r + q^{2} + r^{2}\right)}{r^{6}}$ | $0$ | $\frac{q^{2} \left(- 2 m r + q^{2} + r^{2}\right)}{r^{6}}$ |
| $G_{r r}$ | $- \frac{q^{2}}{r^{2} \left(- 2 m r + q^{2} + r^{2}\right)}$ | $0$ | $- \frac{q^{2}}{r^{2} \left(- 2 m r + q^{2} + r^{2}\right)}$ |
| $G_{\theta \theta}$ | $\frac{q^{2}}{r^{2}}$ | $0$ | $\frac{q^{2}}{r^{2}}$ |
| $G_{\phi \phi}$ | $\frac{q^{2} \sin^{2}{\left(\theta \right)}}{r^{2}}$ | $0$ | $\frac{q^{2} \sin^{2}{\left(\theta \right)}}{r^{2}}$ |

混合成分は $G^t{}_t=G^r{}_r=-q^2/r^4$、$G^\theta{}_\theta=G^\phi{}_\phi=q^2/r^4$ となる。従って $8\pi G$ で割れば静止枠で、電場方向の圧力が密度の負、横の圧力が密度に等しい。混合成分の跡は零なので $R=0$ となるが、電場が非零なら Ricci は零でない。
:::

### 二つの地平線での亀座標の対数

$0<q<m$ の非退化荷電解で $dr_*/dr=f^{-1}=r^2/[(r-r_+)(r-r_-)]$ とする。部分分数分解では、定数部分 $1$ を除いた二項の係数は、分母の各根で $r^2/(r-r_\mp)$ を評価して $r_+^2/(r_+-r_-)$ と $-r_-^2/(r_+-r_-)$ である。従って
$$
r_*=r+\frac{r_+^2}{r_+-r_-}\log|r-r_+|
-\frac{r_-^2}{r_+-r_-}\log|r-r_-|+\mathrm{const}.
$$
対数の引数は共通の長さで割り、変わる定数は最後の積分定数に含める。両根で対数項の符号が逆になるのは $f'(r_\pm)=(r_\pm-r_\mp)/r_\pm^2$ の符号が異なるためである。$q=m$ の極限解では分母が二重根で、直接積分すると $r_*=r+2m\log|r-m|-m^2/(r-m)+\mathrm{const}$ になる。非退化の地平線の指数拡張を考えるときは、対数の係数が発散する極限をそのまま使わず、二重根の式から出発する。

### 内部地平線を横切る零座標と周波数

非退化な二地平線をもつ場合に $\kappa_-=(r_+-r_-)/(2r_-^2)>0$ と置く。二根の間で $u=t-r_*$、$v=t+r_*$ とし、$u$ を有限に保って $v\to+\infty$ に向かう枝を考える。部分分数の式から $r_*=-(2\kappa_-)^{-1}\log(r-r_-)+\order{1}$ である。従って $r-r_-\propto e^{-\kappa_-(v-u)}$ となる。新しい座標を $V_-=-e^{-\kappa_-v}$ とすれば、$\dd v=-\dd V_-/(\kappa_- V_-)$ である。$-f\dd u\dd v$ の係数 $f/(\kappa_-V_-)$ は $V_-\to0$ で非零の有限値を持つ。従ってこの枝は理想解では正則であり、$v$ が大きくなるたびに同じ $\dd v$ の間隔が小さな $\dd V_-$ の間隔へ写る。

:::example[内側境界へ向かう観測者が入射光を測る]
入射する動径の波の位相を $-\omega_v v$、$\omega_v>0$ とする。$V_-$ が零へ近づくときに $\dd V_-/\dd\tau$ が有限かつ非零の観測者では、測定角振動数は $\omega_{\rm obs}=\omega_v\dd v/\dd\tau$ である。従って $\omega_{\rm obs}$ は $e^{\kappa_-v}$ に比例する。例えば $q=4m/5$ なら $r_+=8m/5,r_-=2m/5$ より $\kappa_-=15/(4m)$ である。観測者を同じ枝でさらに $\Delta v=(4m/15)\log2$ だけ進めれば、先頭の青方偏移の因子は二倍になる。外部の入射振幅が $v$ の冪で減るだけなら、この指数効果との積は小さくならない。完全な電磁真空解の座標延長が可能であることと、摂動を入れた時空でも同じ正則性が保たれることは別の問題である。
:::

### 極限解の喉と電荷の規格化

電磁真空で $q=m>0$ の場合は $f=(1-m/r)^2$ である。$r=m+\epsilon y$、$t=T/\epsilon$ と書き、$\epsilon$ を無次元とする。$y>0$ を固定して $\epsilon\to0$ を取ると
$$
\dd s^2\longrightarrow-\frac{y^2}{m^2}\dd T^2
+\frac{m^2}{y^2}\dd y^2+m^2\dd\Omega^2 .
$$
初めの二項で $z=m^2/y$ とすれば $(m^2/z^2)(-\dd T^2+\dd z^2)$ となる。二次元の計量から $R=-2/m^2$、球面から $R=+2/m^2$ なので、四次元で両者の和は零である。この打ち消しは電磁場の跡が零であることと一致する。極限の球面半径は一定で、静的な切片で $r=m$ から一定の半径まで測る距離は $\int\dd r/\sqrt f$ の対数発散を持つ。非退化な根の近傍では $\int\dd r/\sqrt{r-r_+}$ が収束するため、地平線付近の静的な切片は異なる。

SI で電荷を $Q_{\rm SI}$ と書く場合、電場のエネルギー密度は $\epsilon_0 E_{\rm SI}^2/2$、$E_{\rm SI}=Q_{\rm SI}/(4\pi\epsilon_0r^2)$ である。質量関数の式へ入れるには密度を $c^2$ で割り、計量へさらに $G/c^2$ を掛ける。従って $q^2=GQ_{\rm SI}^2/(4\pi\epsilon_0c^4)$ である。$m=GM/c^2$ と比べれば、極限電荷の大きさは $\sqrt{4\pi\epsilon_0 G}\,M$ になる。単位系を変えても比較するのは長さの比 $q/m$ である。

### 電荷が内側の曲率を変える

Reissner–Nordström の $f=1-2m/r+q^2/r^2$ を球対称の曲率の二乗縮約 $f''{}^2+4(f'/r)^2+4[(1-f)/r^2]^2$ に入れる。$f'=2m/r^2-2q^2/r^3$、$f''=-4m/r^3+6q^2/r^4$ を使えば
$$
R_{abcd}R^{abcd}=\frac{48m^2}{r^6}-\frac{96mq^2}{r^7}+\frac{56q^4}{r^8}
$$
となる。外部で得た式は、根をまたいでも入射座標の滑らかな曲率の式と一致する。$r_\pm>0$ では有限だが、$q\ne0$ なら中心に向かう主項は $r^{-8}$ である。電磁場の応力の跡が零で $R=0$ であっても、Riemann テンソルの全体は発散する。

さらに $r<r_-$ では $f>0$ なので、曲面 $r=\mathrm{const.}$ の法線 $\dd r$ は空間的である。従って $r=0$ に近づく等半径面は時間的になる。Schwarzschild の未来内部では $f<0$ で等半径面が空間的だった。二つの内側を同じ因果図で描けないのは、内側の根が法線の符号をもう一度変えるためである。

:::example[内側に反射点を持つ中性粒子]
$E=1,L=0$ の時間的粒子の動径式は $\dot r^2=1-f=2m/r-q^2/r^2$ である。外から落ちる枝には $r_{\rm turn}=q^2/(2m)$ で転回点がある。$r_+r_-=q^2$ と $r_+<2m$ から $r_{\rm turn}<r_-$ である。例えば $q=4m/5$ なら $r_{\rm turn}=8m/25$ で、$r_-=2m/5$ より小さい。転回点で $\ddot r=-f'/2$ は正であり、理想解の粒子は曲率特異性に達する前に外向きとなる。これは特定の測地線の結果であって、全ての因果曲線が中心を避けることは意味しない。例えば保存エネルギーが非零の動径光は $\dot r^2=\mathcal E^2$ であり、この転回点を持たない。
:::

## Kerr 解とその構造

回転解には楕円体に適合した Boyer–Lindquist 座標を用いる。$\Sigma=r^2+a^2\cos^2\theta$、$\Delta=r^2-2mr+a^2$、$s=\sin\theta$ と置く。:term[Kerr 解]{en="Kerr solution" yomi="かーかい"} の計量は

$$
\dd s^2=-\frac{\Delta}{\Sigma}(\dd t-as^2\dd\phi)^2
+\frac{s^2}{\Sigma}((r^2+a^2)\dd\phi-a\dd t)^2
+\Sigma\left(\frac{\dd r^2}{\Delta}+\dd\theta^2\right).
\label{eq-gr-kerr}
$$

$a=J/M$ で、遠方で $g_{t\phi}=-2GJ\sin^2\theta/r+\order{r^{-2}}$ となる。$a=0$ では Schwarzschild に戻る。Einstein 方程式の厳密な積分による解の構成は定常軸対称真空の境界値問題を要するので :cite[wald-gr] を参照する。この章では計量から地平線と軌道を導く。$g^{rr}=\Delta/\Sigma$ より、$0<\abs a<m$ では二つの正の根 $r_\pm=m\pm\sqrt{m^2-a^2}$ がある。$a=0$ での根 $r_-=0$ には曲率特異点があり、地平線として扱えない。変換 $\dd v=\dd t+(r^2+a^2)\dd r/\Delta$、$\dd\tilde\phi=\dd\phi+a\dd r/\Delta$ を代入すると、一つ目の括弧は $\dd v-as^2\dd\tilde\phi-\Sigma\dd r/\Delta$、二つ目は $((r^2+a^2)\dd\tilde\phi-a\dd v)$ となる。$\dd r^2/\Delta$ が打ち消され、$\Sigma>0$ で地平線を越える。$\Sigma=0$ の $r=0,\theta=\pi/2$ は環特異点である。例えば赤道上で曲率二乗縮約は $48m^2/r^6$ となる。$t,\phi$ ブロックの行列式は $-\Delta s^2$、全体は $-\Sigma^2s^2$ であり、極の $s=0$ は角座標の退化である。

### 回転の測定と二種類の円運動

外側で角運動量が零の観測者を $\partial_t+\omega\partial_\phi$ の方向に取り、$u_\phi=0$ とすると $\omega=-g_{t\phi}/g_{\phi\phi}$ である。遠方で $g_{\phi\phi}\simeq r^2\sin^2\theta$、$g_{t\phi}\simeq-2ma\sin^2\theta/r$ なので $\omega\simeq2ma/r^3$ である。地平線に近づけると $\Delta$ を含む項が消え $\omega\to a/(r_+^2+a^2)$ となる。回転とともに動くこの観測者と、角座標が固定された観測者の四元速度は異なる。

赤道内の円**測地線**の角速度 $\Omega=\dd\phi/\dd t$ は、一般にこの $\omega$ と一致しない。一定の $r,\theta$ に対し変分の動径式は $g_{tt,r}+2\Omega g_{t\phi,r}+\Omega^2 g_{\phi\phi,r}=0$ である。Kerr の赤道では $g_{tt}=-1+2m/r$、$g_{t\phi}=-2ma/r$、$g_{\phi\phi}=r^2+a^2+2ma^2/r$ であり、代入すると $(r^3-ma^2)\Omega^2+2ma\Omega-m=0$ になる。従って

$$
\Omega_\pm=\frac{\pm\sqrt m}{r^{3/2}\pm a\sqrt m}.
$$

円運動が時間的かどうかはこの候補を $g(\partial_t+\Omega\partial_\phi,\partial_t+\Omega\partial_\phi)$ に入れて判定する。例えば $r=4m,a=m/2$ では $\Omega_+=2/(17m)$、$\Omega_-=-2/(15m)$ となる。$a=0$ では大きさが一致し、:eqref[eq-gr-circular-orbit] に戻る。これらは無限遠で規格化した時刻での角速度である。

### 地平線の共回転速度

Kerr の計量の二つの平方を使うと、$r=r_+$ では $\Delta$ を係数に持つ平方が消える。$\chi=\partial_t+\Omega\partial_\phi$ の二乗ノルムは残る平方の係数に $[(r_+^2+a^2)\Omega-a]^2$ を掛けたものになる。従って軸以外で $g(\chi,\chi)=0$ とするには $\Omega_H=a/(r_+^2+a^2)$ が必要である。正則な入射座標では $\chi$ の双対が地平線の $\dd r$ に平行であるため接線であるとともに法線でもある。Killing 場 $\partial_t$ 単独のノルムが零になる面とは一般に異なる。

:::example[外部地平線と静止限界の間隔]
$a=3m/5$ なら $r_+=9m/5$、$r_-=m/5$ である。赤道の静止限界は $2m$ なので間隔は $m/5$、一方で軸の静止限界は $r_+$ に一致する。$\Omega_H=a/(2mr_+)=1/(6m)$ である。$r_+^2+a^2=2mr_+$ は根の方程式を使ったため、電荷が加わると同じ書換えには電荷項が入る。
:::

### 外部で時刻の勾配と回転する観測者を調べる

$g_{\phi\phi}=\mathcal A\sin^2\theta/\Sigma$、$\mathcal A=(r^2+a^2)^2-\Delta a^2\sin^2\theta$ と置く。計量の $t,\phi$ 部分を平方完成すると
$$
-\alpha^2\dd t^2+g_{\phi\phi}(\dd\phi-\omega\dd t)^2,\qquad
\alpha^2=\frac{\Delta\Sigma}{\mathcal A},\qquad
\omega=\frac{2mar}{\mathcal A}.
$$
実際、交差項から $\omega=-g_{t\phi}/g_{\phi\phi}$、行列式から $\alpha^2 g_{\phi\phi}=\Delta\sin^2\theta$ である。外部で $\Delta>0$、$\mathcal A>0$ なら $g^{tt}=-1/\alpha^2<0$ である。従って時刻一形式 $\dd t$ はエルゴ領域でも時間的であり、一般の未来因果曲線に沿って $t$ が増える。一方、座標ベクトル $\partial_t$ のノルムは $-\alpha^2+g_{\phi\phi}\omega^2$ であり、正にもなる。時刻の勾配と座標ベクトルのノルムが違うのは交差項があるためである。

:::example[零角運動量の観測者による光の向きの測定]
観測者を $n=\alpha^{-1}(\partial_t+\omega\partial_\phi)$ とし、外部の局所直交枠の方位成分を $p_{\hat\phi}=L/\sqrt{g_{\phi\phi}}$ とする。この観測者の測るエネルギーは $\mathcal E=(E-\omega L)/\alpha$ である。例えば局所でちょうど方位方向に進む光なら $p_{\hat\phi}=\pm\mathcal E$ なので $E=\mathcal E(\alpha\pm\omega\sqrt{g_{\phi\phi}})$ である。$a>0$ のエルゴ領域では $\omega\sqrt{g_{\phi\phi}}>\alpha$ なので、逆向きの光は $E<0$ となる。局所エネルギー $\mathcal E$ はどちらも正である。地平線付近で角座標の速度は $\omega$ に近づくため、局所で逆向きに進む光でも $\dd\phi/\dd t$ が正になり得る。
:::

### Kerr の地平線上の非アフィン性

入射座標で $\chi=\partial_v+\Omega_H\partial_{\tilde\phi}$ を使う。計量から地平線上で $\chi_b=C(\theta)(\dd r)_b$、$C=1-a\Omega_H\sin^2\theta=\Sigma_+/(r_+^2+a^2)$ である。Killing 方程式により $\chi^a\nabla_a\chi_b=-\partial_b(\chi^2)/2$ となる。$\chi^2$ の二つの平方のうち、$\Delta$ が掛からない方は $r=r_+$ で二乗の中身が零なので、そこでの一次微分には寄与しない。従って $\partial_r(\chi^2)=-\Delta'(r_+) C^2/\Sigma_+$ であり、
$$
\chi^a\nabla_a\chi_b=\kappa_H\chi_b,\qquad
\kappa_H=\frac{\Delta'(r_+)}{2(r_+^2+a^2)}.
$$
非退化なら $\Delta'(r_+)=r_+-r_->0$ である。緯度依存が相殺するため、すべての生成子で同じ係数になる。非退化でのアフィン化には Killing 流のパラメータに沿う指数因子が必要である。極限で $\kappa_H=0$ になると同じ指数による変換をそのまま座標として用いられない。

## エルゴ領域と Penrose 過程

外側で $g_{tt}=-(\Delta-a^2s^2)/\Sigma$ が零となる面は $r_{\rm e}(\theta)=m+\sqrt{m^2-a^2\cos^2\theta}$ である。$r_+<r<r_{\rm e}$ を :term[エルゴ領域]{en="ergoregion" yomi="えるごりょういき"} という。静止方向 $\partial_t$ が空間的なので、未来時間的な運動量 $p$ でも $E=-p_t$ を負にできる。実際、空間的ベクトルと内積が正になる未来時間的ベクトルは、その空間方向に十分大きな速度を選ぶと得られる。粒子がエルゴ領域で分裂し負の $E$ の片方が落ち、もう片方が脱出すれば、保存 $E_{\rm in}=E_{\rm fall}+E_{\rm esc}$ から $E_{\rm esc}>E_{\rm in}$ となる。これが :term[Penrose 過程]{en="Penrose process" yomi="ぺんろーずかてい"} の原理である。分裂が任意の初期値で実現するとは限らず、運動量保存と両軌道の許容領域も必要である。

地平線生成子 $\chi=\partial_t+\Omega_H\partial_\phi$ では、:eqref[eq-gr-kerr] の二つ目の括弧が零になるよう $\Omega_H=a/(r_+^2+a^2)$ を選ぶ。未来因果的粒子には $-p\cdot\chi=E-\Omega_H L\ge0$ が必要である。$a>0$ で負エネルギーを持ち込む粒子の $L$ は負で、回転も減少する。例えば $a=4m/5$ なら $r_+=8m/5$、$\Omega_H=1/(4m)$ である。赤道の外側静止限界は $2m$ なので、二面の差は $2m/5$ になる。

### 地平線を横切る粒子の不等式

未来地平線の生成子 $\chi=\partial_t+\Omega_H\partial_\phi$ と未来因果的な中性粒子の運動量 $p$ に対し $g(p,\chi)\le0$ である。局所直交枠で二つの未来因果ベクトルの空間成分の積は時間成分の積以下なので、この符号は場の方程式を用いずに従う。保存量 $E=-p_t,L=p_\phi$ により条件は $E-\Omega_HL\ge0$ となる。従って $a>0$ で $E<0$ の粒子が落ちるなら $L<0$ が必要である。分裂の前後で粒子の総運動量が保存すれば、逃げる側の $E$ は入射側の $E$ より大きくなり得る。ただしこの不等式は横切るための必要条件であり、分裂点からの軌道が実際に地平線または無限遠へ達するかは動径ポテンシャルで判定する。

### 静止枠を持てなくなるが未来因果的方向は残る

:::example[エルゴ領域で角速度の符号を判定する]
赤道で $a=4m/5,r=9m/5$ とする。$r_+=8m/5<r<2m$ なので $g_{tt}=-1+2m/r=1/9>0$ である。さらに $g_{\phi\phi}=r^2+a^2+2ma^2/r>0$、$g_{t\phi}=-2ma/r<0$ である。$\partial_t+\Omega\partial_\phi$ が時間的になる範囲の両端は二次式の根である。根の積は $g_{tt}/g_{\phi\phi}>0$、和は $-2g_{t\phi}/g_{\phi\phi}>0$ なので、両端の角速度は正である。$\Delta>0$ により判別式も正であり、二根の間に未来向きに選べる時間的方向が残る。従ってこの点で時間的に座標半径を保つには、回転と同じ向きに $\phi$ を変える必要がある。
:::

### 保存エネルギーの分配で局所的な上限を求める

質量 $\mu$ の親粒子を二つの光へ分裂させる模型を考える。分裂点で親の四元速度を $u$、$u$ に直交する単位空間方向を $s$ とすれば、両方の光の運動量は $p_\pm=(\mu/2)(u\pm s)$ である。局所運動量保存を満たし、一方の光だけの局所エネルギーを独立に変える自由度はない。$\xi=\partial_t$、$E_0=-\mu u\cdot\xi$ とすると
$$
E_\pm=\frac{E_0}{2}\mp\frac{\mu}{2}s\cdot\xi,\qquad
|\xi+(u\cdot\xi)u|^2=\xi^2+(E_0/\mu)^2 .
$$
従って、$s$ を空間射影された $\xi$ に沿わせると $|s\cdot\xi|$ は最大になり、一方の光の保存エネルギーは
$$
E_{\rm max}=\frac12\left[E_0+\sqrt{E_0^2+\mu^2g_{tt}}\right]
$$
を越えない。未束縛の親粒子で $E_0\ge\mu$ とすれば、エルゴ領域の中でも利得の割合は有限である。

:::example[赤道の地平線に近づけた局所上界]
赤道では $g_{tt}=2m/r-1$ であり、外部で可能な最大値は極限 Kerr の $r\downarrow m$ で $1$ へ近づく。遠方で静止から来た親粒子の $E_0=\mu$ を使うと、二光子への分裂模型での局所上界は $E_{\rm max}/E_0\le(1+\sqrt2)/2$ である。余分に得られるエネルギーの比は約 $0.207$ 以下である。ただしこの上界は一点での運動量保存だけを使った。外へ出る光が分離式の許容領域に沿って無限遠へ届くか、負の側が落ちるかは別の条件である。この局所計算は任意の軌道に対して利得を達成できることを保証しない。
:::

### 面積が減少しない向きに回転を変える

ここでは Kerr 族のパラメータだけを微分する。$j=GJ=ma$ と書くと $\Delta(r_+)=r_+^2-2mr_++j^2/m^2=0$ であり、$A=4\pi(r_+^2+a^2)$ である。$m,a$ を変えて根の式を微分すると $(r_+-m)\delta r_+=r_+\delta m-a\delta a$ である。また $\delta j=a\delta m+m\delta a$ なので
$$
\frac{\kappa_H}{8\pi}\delta A=\delta m-\Omega_H\delta j .
$$
この等式には非退化な外側根を使った。捕獲された中性の微小試験粒子の保存量を $\delta m=GE$、$\delta j=GL$ と同一視すれば、地平線での $E-\Omega_H L\ge0$ から、吸収前後の定常解の面積差が一次で非負になる。これは定常解族と試験粒子についての微分計算である。任意の時間依存過程の面積に関する定理は [ブラックホール熱力学](/physics/general-relativity/bh-thermodynamics/) で扱う。

$M_{\rm irr}=\sqrt{A/(16\pi G^2)}$ と定義すると、$A/(16\pi)=mr_+/2$ であり、$r_+^2-2mr_++a^2=0$ から
$$
M^2=M_{\rm irr}^2+\frac{J^2}{4G^2M_{\rm irr}^2}
$$
が従う。例えば極限で $r_+=m$ なら $M_{\rm irr}=M/\sqrt2$ である。面積を保ったまま自転を取り除ければ、最終の静止穴の質量は $M/\sqrt2$ になる。失える質量の相対差 $1-1/\sqrt2\simeq0.293$ は、先の一回の分裂で得られる割合と異なる。前者は Kerr 族で回転を消す全体の変更、後者は一事象での粒子エネルギーの分配に由来する。

## Kerr 時空の測地線と Carter 定数

:term[Hamilton–Jacobi 方程式]{en="Hamilton–Jacobi equation" yomi="はみるとんやこびほうていしき"} $g^{ab}\partial_a S\partial_b S=-\mu^2$ を $S=-Et+L\phi+S_r(r)+S_\theta(\theta)$ で解く。計量を双対枠で反転し、$\Sigma$ を掛けると

$$
\Delta(S_r')^2-\frac{P(r)^2}{\Delta}+
(S_\theta')^2+(L/\sin\theta-aE\sin\theta)^2+\mu^2\Sigma=0,\qquad
P=E(r^2+a^2)-aL .
$$

$r,\theta$ の項を分け、角側の定数から $(L-aE)^2$ を引いたものを $\mathcal Q$ とする。$p_r=\Sigma\dot r/\Delta$、$p_\theta=\Sigma\dot\theta$ なので

$$
\begin{aligned}
(\Sigma\dot r)^2 &=P^2-\Delta[\mu^2r^2+(L-aE)^2+\mathcal Q],\\
(\Sigma\dot\theta)^2 &=\mathcal Q-\cos^2\theta[a^2(\mu^2-E^2)+L^2/\sin^2\theta].
\end{aligned}
\label{eq-gr-carter}
$$

任意の初期運動量で分離定数を選べ、Hamilton–Jacobi の特性曲線では一定である。直接にも Hamiltonian を $\Sigma H=(A_r+A_\theta)/2$ と分けて、$\mathcal Q=A_\theta-2Ha^2\cos^2\theta-(L-aE)^2$ と置けば、Hamilton の式で $\dot A_\theta=2H\frac{\dd}{\dd\lambda}(a^2\cos^2\theta)$ が従う。従って :term[Carter 定数]{en="Carter constant" yomi="かーたーていすう"} は運動量の二次形式であり、[対称性と Killing ベクトル](/physics/general-relativity/symmetries-killing/) の Killing テンソルの判定法を使える。$a=0$ では $\mathcal Q=L_{\rm total}^2-L^2$ となる。赤道面内では $\mathcal Q=0$、逆に赤道を横切るときは $\mathcal Q=p_\theta^2$ である。

### 時刻と方位角を積分する

Hamiltonian $H=g^{ab}p_ap_b/2$ を使うと $\dot t=\partial H/\partial p_t=-\partial H/\partial E$、$\dot\phi=\partial H/\partial L$ である。分離式のうち $E,L$ を含む部分は $2\Sigma H=-P^2/\Delta+(L/\sin\theta-aE\sin\theta)^2+\cdots$ である。従って
$$
\Sigma\dot t=\frac{r^2+a^2}{\Delta}P+a(L-aE\sin^2\theta),\qquad
\Sigma\dot\phi=\frac{a}{\Delta}P+\frac{L}{\sin^2\theta}-aE.
$$
必要なら $\dd\lambda=\Sigma\,\dd\gamma$ と新しい径数 $\gamma$ を取る。このとき $r(\gamma)$ と $\theta(\gamma)$ は :eqref[eq-gr-carter] の各平方根で独立に進み、残る二式の右辺はそれぞれ $r$ の関数と $\theta$ の関数の和になる。転回点では平方根の符号を切り替え、位置と運動量を連続にする。例えば $a=0$ で赤道に制限すると $\Sigma=r^2,\Delta=r^2f$ なので $\dot t=E/f$、$\dot\phi=L/r^2$ が復元する。地平線に近い点ではこれらの座標微分は発散し得るため、地平線の横断には先に作った入射座標を用いる。

### 赤道面から外れる初期速度

Carter の分離定数を $Q=p_\theta^2+\cos^2\theta[a^2(\mu^2-E^2)+L^2/\sin^2\theta]$ とする。赤道 $\theta=\pi/2$ では $Q=p_\theta^2$ であるため、赤道で $p_\theta\ne0$ から始める軌道は $Q>0$ が必要となる。逆に $Q=0$ かつ初期点が赤道なら $p_\theta=0$ であり、反射対称性と測地線の一意性によって赤道に留まる。赤道からの小さな角変位 $\delta$ について $\cos^2\theta=\delta^2+\order{\delta^4}$、$\sin^{-2}\theta=1+\delta^2+\order{\delta^4}$ なので、横方向ポテンシャルの二次係数は $a^2(\mu^2-E^2)+L^2$ である。例えば $E^2<\mu^2$ の束縛軌道で $a,L$ の少なくとも一方が非零なら、この係数は正になり、小さな $Q$ は有限の緯度振幅に対応する。

### 円軌道の半径を保存量の二条件から求める

赤道の中性時間的粒子を単位静止質量で規格化する。円軌道では $\dot\phi=\Omega\dot t$ とし、時刻・角成分の規格化から
$$
u^t=(-g_{tt}-2\Omega g_{t\phi}-\Omega^2g_{\phi\phi})^{-1/2},\qquad
E=-(g_{tt}+\Omega g_{t\phi})u^t,\quad L=(g_{t\phi}+\Omega g_{\phi\phi})u^t .
$$
動径の Euler–Lagrange 式では $\partial_r g_{tt}+2\Omega\partial_r g_{t\phi}+\Omega^2\partial_r g_{\phi\phi}=0$ である。$g_{tt}=-1+2m/r,g_{t\phi}=-2ma/r,g_{\phi\phi}=r^2+a^2+2ma^2/r$ を代入すると $r^3\Omega^2=m(1-a\Omega)^2$ である。符号 $s=\pm1$ の二枝は $\Omega=s\sqrt m/(r^{3/2}+s a\sqrt m)$ となる。これを規格化に戻すと、$x=\sqrt{r/m}$、$j=a/m$ として
$$
E=\frac{x^3-2x+sj}{x^{3/2}\sqrt{x^3-3x+2sj}},
\qquad
L=\frac{s m(x^4-2sjx+j^2)}{x^{3/2}\sqrt{x^3-3x+2sj}} .
$$
分母の根号が正である範囲を使い、安定性は固定した $E,L$ での動径関数 $\mathcal R''$ で判定する。円軌道条件を先に $\mathcal R$ に代入してから $r$ で微分すると、別の円軌道へ移ったときの $E,L$ の微分まで含む。その結果は固定した保存量での摂動を表さない。

### 保存量を用いた緯度の転回点

:::example[非回転極限で軌道面の傾きを求める]
$a=0$ で $\mathcal Q=L_z^2$ とする。角方程式は $p_\theta^2=\mathcal Q-L_z^2\cot^2\theta$ なので、転回では $\cot^2\theta=1$ である。従って $\theta=\pi/4,3\pi/4$ の間を運動する。全角運動量の二乗は $L_{\rm tot}^2=\mathcal Q+L_z^2=2L_z^2$ なので、軌道面の法線と対称軸の角は $L_z>0$ のとき $\arccos(L_z/L_{\rm tot})=\pi/4$ である。球対称性によって面を回転すれば赤道運動に戻るが、元の軸を固定して記述した軌道では $\theta$ は一定でない。
:::

### 安定な円軌道の端を固定保存量の微分から選ぶ

赤道面の中性時間的軌道で $\mu=1$ とし、$B=L-aE$ を置く。$P=Er^2-aB$ なので :eqref[eq-gr-carter] の動径多項式は
$$
\mathcal R=(E^2-1)r^4+2mr^3-(B^2+2aEB+a^2)r^2+2mB^2r .
$$
正の半径では $\mathcal R/r$ を使って $\mathcal R=\mathcal R'=0$ と同値な二式を作れる。円軌道において $(\mathcal R/r)''=6(E^2-1)r+4m$ であり、$\mathcal R''=r(\mathcal R/r)''$ である。微小な動径変位 $\delta r$ には $\delta\ddot r=\mathcal R''\delta r/(2r^4)$ となるので、安定な円運動は $\mathcal R''<0$ を要する。等号が動径安定性の端を表す。ここで点は粒子の固有時であり、円軌道から離れたときの $E,L$ の変化は微分に入れない。

前に求めた $E$ を入れる。$x=\sqrt{r/m}$、$j=a/m$、$s=\pm1$ とすると、$\mathcal R''=0$ は
$$
x^4-6x^2+8sjx-3j^2=0
$$
と同値になる。途中で割る $x^3-3x+2sj$ は、時間的円軌道が存在する範囲で正である。無限遠の安定枝から半径を下げて最初に安定性を失う根を選ぶ。例えば $j=0$ なら正根は $x=\sqrt6$ である。$j=1$ の順行枝の多項式は $(x-1)^3(x+3)$ だが、$x=1$ では同時に分母が零になる。この根を通常の外部の時間的円軌道に直接代入できない。非極限の順行枝を連続に $j\uparrow1$ として調べる。一方、$j=1,s=-1$ の式は $(x+1)^3(x-3)$ であり、逆行枝の端は $r=9m$ となる。

:::example[半径 $4m$ の順行軌道が限界になるスピン]
$x=2,s=+1$ を条件に入れると $-8+16j-3j^2=0$ である。$0\le j\le1$ にある根は $j=(8-2\sqrt{10})/3\simeq0.5585$ である。このとき $x^3-3x+2j=2+2j>0$ なので時間的軌道のエネルギーが定義できる。逆行側では同じ半径で多項式は $-8-16j-3j^2<0$ なので、そこは安定枝の端ではない。円軌道としての存在と、動径摂動に対する安定性を個別に確かめる。
:::

### 光子の球面軌道と遠方の影の境界

Kerr では光が一定半径を保つ場合にも $\theta$ は変化できる。$\mu=0$、$E>0$ とし、$\xi=L/E$、$\eta=\mathcal Q/E^2$ と置く。動径式で $\mathcal R=\mathcal R'=0$ とすれば
$$
P_0^2=\Delta K,\qquad 4rP_0=\Delta'K,\qquad
P_0=r^2+a^2-a\xi,\qquad K=(\xi-a)^2+\eta .
$$
ここでは外側で $0<|a|<m$、$\Delta>0$ とする。候補軌道には角方向の許容条件も必要なので、$K$ は $[p_\theta^2+(L/\sin\theta-aE\sin\theta)^2]/E^2\ge0$ である。$K=0$ なら両項が零で $\xi=a\sin^2\theta$ となり、$P_0=r^2+a^2\cos^2\theta>0$ に矛盾する。従って $P_0\ne0$ として割り算ができる。$\Delta'=2(r-m)$ により
$$
\xi(r)=\frac{r^2+a^2-4r\Delta/\Delta'}{a},
\qquad
\eta(r)=\frac{16r^2\Delta}{(\Delta')^2}-(\xi(r)-a)^2
$$
となる。$\xi,\eta$ を求めても、$\eta+a^2\cos^2\theta-\xi^2\cot^2\theta\ge0$ の角範囲があるかを確かめる。$\mathcal Q$ はこの定義では常に非負と仮定できない。

十分遠方で傾き $\theta_{\rm o}$ にいる観測者の視野に、軸と直交する横向きの座標 $\alpha$ と、もう一つの横向き座標 $\beta$ を取る。像の横座標を角度と観測半径の積で規格化すれば、$\alpha=-\xi/\sin\theta_{\rm o}$ である。これは静止枠の方位運動量が先頭で $p_\phi/(r_{\rm o}\sin\theta_{\rm o})$、エネルギーが $E$ になるためである。同様に $\beta=\pm p_\theta/E$ だから
$$
\beta=\pm\sqrt{\eta+a^2\cos^2\theta_{\rm o}-\xi^2\cot^2\theta_{\rm o}}.
$$
不安定な球面軌道に漸近する光は、遠方から逆向きに追った場合に、地平線へ落ちる光と外へ戻る光の境界を与える。係数の計算では、球面軌道がその観測角から到達できることも確かめる。例えば $\eta<0$ の候補は赤道で $\beta^2<0$ になるので、赤道の観測者が見るこの境界には含めない。

:::example[赤道の視野に出る二つの端]
赤道を通る球面軌道で視野の高さが零なら $\eta=0$ である。これを上式へ代入すると $r(r-3m)^2=4ma^2$ が従う。$x=\sqrt{r/m}$ で分ければ、赤道の順行・逆行の枝は $x^3-3x+2s a/m=0$ の外部根である。$a=0$ では両枝とも $r=3m$ に退化する。この極限で $\xi$ の式に $a=0$ を直接入れると不定になるので、Schwarzschild の動径式に戻る。その極限の関係は $\xi^2+\eta=27m^2$ であり、視野の半径も $3\sqrt3m$ となる。回転時には一定の $r$ の軌道でも傾きを変えられるため、視野の境界を一つの円軌道だけでは与えられない。
:::

### 定常解の対称性から経路の数値積分を組み立てる

径数 $\gamma$ を $\dd\lambda=\Sigma\dd\gamma$ で選ぶ。$R=(\dd r/\dd\gamma)^2$ と $\Theta=(\dd\theta/\dd\gamma)^2$ を二つの分離ポテンシャルで表すと、途中で $\dd r/\dd\gamma=0$ を通る数値積分では平方根の向きを決める必要がある。単純根をまたぐには $(\dd^2 r/\dd\gamma^2)=R'(r)/2$ と二階式で進めるのが一つの方法である。導出は $\dd r/\dd\gamma\ne0$ の区間で微分し、単純根まで連続延長する。$\theta$ も同様である。構成した二階式では $[\dd r/\dd\gamma]^2-R(r)$ の微分が恒等的に零なので、初期データが一次式を満たせば厳密解でも満たす。数値解で残差が大きくなると、刻みや初期値を見直せる。

$(r,\dd r/\dd\gamma)$ の初期値が多項式の二重根 $(r_0,0)$ なら、二階式の一意性により $r=r_0$ に留まる。微小変位に対する式は $\delta r''=R''(r_0)\delta r/2$ であり、$R''>0$ で指数的に広がる。平方根式で単純根と二重根を同じアルゴリズムで扱うと、二重根で静止すべき解を誤って枝変更させる場合がある。また位置が地平線に近い場合、$t,\phi$ には座標による発散が現れる。通過を含めるには入射座標の微分式に変換してから積分する。

## Kerr–Newman 解

電荷と回転が共存する :term[Kerr–Newman 解]{en="Kerr–Newman solution" yomi="かーにゅーまんかい"} では :eqref[eq-gr-kerr] の $\Delta$ を $r^2-2mr+a^2+q^2$ にし、電磁ポテンシャルを $A=-(Q_0r/\Sigma)(\dd t-a\sin^2\theta\,\dd\phi)$ とする。$F=\dd A$ と置くと $a=0$ で $F_{tr}=-Q_0/r^2$ になり、$q=0$ で Maxwell 場が消える。地平線が実根を持つには $m^2\ge a^2+q^2$ が必要である。例えば $a=3m/5,q=4m/5$ では二根が $m$ に重なる。計量の行列式と入射座標での相殺は同じである。面積は地平線上の角計量から $\sqrt{g_{\theta\theta}g_{\phi\phi}}=(r_+^2+a^2)\sin\theta$ を積分して $4\pi(r_+^2+a^2)$ となる。中性粒子の分離は $\Delta$ の置換だけで保たれる。電荷 $e$ の粒子の場合、力学的運動量は $\partial_aS-eA_a$ なので $P$ が $P-eQ_0r$ に替わる。角側では $A_\phi+a\sin^2\theta A_t=0$ なので追加項は相殺する。

### 電磁真空の方程式を枠で確認する

地平線の外で $\sigma_0=\dd t-a\sin^2\theta\dd\phi$、$\sigma_3=(r^2+a^2)\dd\phi-a\dd t$ とし、$\rho_\Sigma=\sqrt{\Sigma}$、$D=\sqrt{\Delta}$ と置く。この計量の正規直交枠は $e^0=D\sigma_0/\rho_\Sigma$、$e^1=\rho_\Sigma\dd r/D$、$e^2=\rho_\Sigma\dd\theta$、$e^3=\sin\theta\,\sigma_3/\rho_\Sigma$ である。$A=-Q_0r\sigma_0/\Sigma$ を微分する。$\Sigma\dd\phi=\sigma_3+a\sigma_0$ を使うと、角微分の $\dd\theta\wedge\sigma_0$ が消えるので

$$
F=-\frac{Q_0(r^2-a^2\cos^2\theta)}{\Sigma^2}e^0\wedge e^1+
\frac{2Q_0ra\cos\theta}{\Sigma^2}e^2\wedge e^3 .
$$

二係数の二乗和は $Q_0^2/\Sigma^2$ である。従って枠の応力は $T_{AB}=Q_0^2\diag(1,-1,1,1)/(2\Sigma^2)$ である。Maxwell 方程式の確認には、逆計量で両添字を上げて $\partial_r(\Sigma\sin\theta F^{r b})+\partial_\theta(\Sigma\sin\theta F^{\theta b})=0$ を計算する。例えば $b=t$ では二つの被微分項は

$$
\Sigma\sin\theta F^{rt}=-\frac{Q_0(r^2+a^2)(r^2-a^2\cos^2\theta)\sin\theta}{\Sigma^2},\qquad
\Sigma\sin\theta F^{\theta t}=\frac{2Q_0ra^2\cos\theta\sin^2\theta}{\Sigma^2}.
$$

$r,\theta$ の微分を取ると符号が反対である。$b=\phi$ では第一項の $r^2+a^2$ が $a$ に、第二項の $a^2\sin^2\theta$ が $a$ に替わって同じ打ち消しが起こる。残る成分は $F^{r\theta}=0$ および定常・軸対称性から零である。

重力側を再現するため、$\sin\theta=s,\cos\theta=c_\theta$ と略し、次の六つの関数を置く。

$$
\begin{aligned}
A_1&=(\Delta'/(2\Delta)-r/\Sigma)D/\rho_\Sigma,&
A_2&=a^2s c_\theta/\rho_\Sigma^3,&
A_3&=2a c_\theta D/\rho_\Sigma^3,\\
A_4&=rD/\rho_\Sigma^3,&
A_5&=2ras/\rho_\Sigma^3,&
A_6&=\cot\theta/\rho_\Sigma+A_2 .
\end{aligned}
$$

$\dd r=(D/\rho_\Sigma)e^1$ と $\dd\theta=e^2/\rho_\Sigma$ を使って各枠を微分すると

$$
\begin{aligned}
\dd e^0&=A_1 e^1\wedge e^0-A_2e^2\wedge e^0-A_3e^2\wedge e^3,\\
\dd e^1&=-A_2e^2\wedge e^1,\qquad
\dd e^2=A_4e^1\wedge e^2,\\
\dd e^3&=A_4e^1\wedge e^3+A_6e^2\wedge e^3+A_5e^1\wedge e^0.
\end{aligned}
$$

双対ベクトルを $E_A$ とし、$C^A{}_{BC}=-\dd e^A(E_B,E_C)$ と置く。:eqref[eq-gr-cartan] の第一式と計量整合性を解くと $\omega_{AB}(E_C)=(C_{ACB}-C_{CBA}+C_{BAC})/2$ である。次に

$$
R^A{}_{BCD}=E_C\omega^A{}_B(E_D)-E_D\omega^A{}_B(E_C)
-C^J{}_{CD}\omega^A{}_B(E_J)
+\omega^A{}_J(E_C)\omega^J{}_B(E_D)-\omega^A{}_J(E_D)\omega^J{}_B(E_C)
$$

に代入する。係数は $r,\theta$ の関数だから微分は $E_1=(D/\rho_\Sigma)\partial_r$、$E_2=\rho_\Sigma^{-1}\partial_\theta$ のみでよい。例えば Ricci の $00$ 成分に残る組は

$$
R_{00}=E_1A_1-E_2A_2+A_1^2+2A_1A_4+2A_2^2-A_2A_6+(A_3^2-A_5^2)/2 .
$$

$\Delta'=2(r-m)$、$\Delta''=2$ を使って整理すると $R_{00}=q^2/\Sigma^2$ であり、同じ縮約で $R_{11}=-R_{00}$、$R_{22}=R_{33}=R_{00}$、異なる添字の成分は零になる。従って $R=0$、$R_{AB}=8\pi GT_{AB}$ である。$Q_0=0$ なら同じ枠の計算で Ricci がすべて零になり、Kerr の真空方程式も確認できる。

### 地平線上で電位を評価する

定常な規格化でポテンシャルを $A=-Q_0 r(\dd t-a\sin^2\theta\dd\phi)/\Sigma$ とする。無限遠で $A_t\to0$ であり、地平線で $\chi=\partial_t+\Omega_H\partial_\phi$ と縮約すると
$$
-A(\chi)=\frac{Q_0r_+}{\Sigma_+}\left(1-\frac{a^2\sin^2\theta}{r_+^2+a^2}\right)
=\frac{Q_0r_+}{r_+^2+a^2}.
$$
二つ目の等式の分子に $\Sigma_+=r_+^2+a^2\cos^2\theta$ があるため、緯度依存は消える。$a=0$ では $Q_0/r_+$ に戻り、$Q_0=0$ なら零である。電位はゲージの定数変更で移り得るため、第一法則で使うときは無限遠との差を採用する。

:::example[電荷と回転の効果を同じ質量で比較する]
$j=a/m$、$e=q/m$ とし、$j^2+e^2\le1$ とする。外側根は $r_+/m=1+\sqrt{1-j^2-e^2}$ である。面積は $4\pi(r_+^2+a^2)$、角速度は $a/(r_+^2+a^2)$ である。$\kappa$ は地平線の生成子を無限遠の時刻で規格化して $(r_+-r_-)/(2(r_+^2+a^2))$ を使う。

ここでは$j$ は無次元スピン、$e$ は無次元電荷であるとする。計算に用いる式を先に並べる。

$$
\begin{aligned}
r_+/m&=\sqrt{- e^{2} - j^{2} + 1} + 1\\
A_H/(4\pi m^2)&=j^{2} + \left(\sqrt{- e^{2} - j^{2} + 1} + 1\right)^{2}\\
m\Omega_H&=\frac{j}{j^{2} + \left(\sqrt{- e^{2} - j^{2} + 1} + 1\right)^{2}}\\
m\kappa&=\frac{2 \sqrt{- e^{2} - j^{2} + 1}}{2 j^{2} + 2 \left(\sqrt{- e^{2} - j^{2} + 1} + 1\right)^{2}}
\end{aligned}
$$

次の条件を比較する。入力を区別してから計算し、定数の規格化をそろえた値を表に示す。

| $j$ | $e$ | $r_+/m$ | $A_H/(4\pi m^2)$ | $m\Omega_H$ | $m\kappa$ |
|---|---|---|---|---|---|
| $0$ | $0$ | $2$ | $4$ | $0$ | $\frac{1}{4}$ |
| $\frac{3}{5}$ | $0$ | $\frac{9}{5}$ | $\frac{18}{5}$ | $\frac{1}{6}$ | $\frac{2}{9}$ |
| $0$ | $\frac{3}{5}$ | $\frac{9}{5}$ | $\frac{81}{25}$ | $0$ | $\frac{20}{81}$ |
| $\frac{3}{5}$ | $\frac{4}{5}$ | $1$ | $\frac{34}{25}$ | $\frac{15}{34}$ | $0$ |
| $1$ | $0$ | $1$ | $2$ | $\frac{1}{2}$ | $0$ |

回転だけ、電荷だけの例で $r_+$ が一致しても、面積は回転の $a^2$ 項によって異なる。極限では $r_+=r_-=m$ で $\kappa=0$ だが、面積は有限で、$a$ が非零なら角速度も有限である。
:::

### 電荷のある粒子を地平線へ運ぶ条件

力学的な運動量を $\pi_a=p_a-eA_a$、正準運動量を $p_a=\partial_aS$ とする。背景の定常性と軸対称性に適合するゲージで $E=-p_t,L=p_\phi$ を保存する。地平線で $\chi=\partial_t+\Omega_H\partial_\phi$ と $\Phi_H=-A(\chi)$ を用いると、未来因果性が要求するのは
$$
-\pi_a\chi^a=E-\Omega_H L-e\Phi_H\ge0
$$
である。正準運動量をそのまま因果条件に使うと、最後の項が落ちる。実際、先の動径の分離式に現れる $P-eQ_0r$ は地平線で $(r_+^2+a^2)(E-\Omega_H L-e\Phi_H)$ である。

:::example[静電的反発と動径到達条件]
非回転で $Q_0>0$、$e>0$ とする。動径の式は $\dot r^2=(E-eQ_0/r)^2-\mu^2 f$ となる。$E$ の符号だけでは地平線上の未来向きを判定できない。例えば $E=eQ_0/(2r_+)>0$ は $E-eQ_0/r_+<0$ なので、平方した動径式が地平線で正でも未来向きに地平線を横切れない。逆に $E-eQ_0/r_+>0$ でも途中に転回点があれば、外の初期点からそこまで進めない。
:::

### Kerr の面積の微分を電磁真空へ拡張する

$q^2=4\pi GQ_0^2$ の規格化で、根の式を微分すると $(r_+-m)\delta r_+=r_+\delta m-a\delta a-q\delta q$ である。$A=4\pi(r_+^2+a^2)$ と $j=ma$ を用いれば
$$
\delta m=\frac{\kappa_H}{8\pi}\delta A+
\Omega_H\delta j+\frac{r_+}{r_+^2+a^2}\,q\delta q.
$$
$q\delta q=4\pi GQ_0\delta Q_0$ なので、$F=\dd A$ の単位系で電束から得る電荷 $Q=4\pi Q_0$ を用いると、最後の項は $G\Phi_H\delta Q$ となる。$G$ を割って物理的な質量と角運動量に直した式には、$\kappa_H\delta A/(8\pi G)$、$\Omega_H\delta J$、$\Phi_H\delta Q$ が入る。試験粒子を吸収して $\delta M=E,\delta J=L,\delta Q=e$ とした一次の面積変化は、先の未来因果性条件と同じ組合せに従う。極限解を直接微分する場合は、二重根が単純根へ変わるときの微分可能性を別に調べる。

## 無毛定理

:term[無毛定理]{en="black hole uniqueness theorem" yomi="むもうていり"} は、定常なブラックホールの外部解を一定の仮定のもとで分類する。定常流が完備で、外部の通信領域が大域的に双曲的とする。さらにその領域に、漸近端を含む連結な空間的かつ相互に因果的でない切片があり、その閉包の境界は未来地平線の各生成子と一度だけ交わるとする。閉包はコンパクトな領域に有限個の漸近端を付けた境界付き多様体とする。この一組の仮定を $I^+$ 正則性と呼ぶ。四次元の漸近平坦で解析的な定常電磁真空で、この条件と連結な非退化の事象地平線を仮定する形では、外部の通信領域は Kerr–Newman の対応する領域に等長である。真空の非退化の場合は Kerr となる。磁気電荷を許せば双対回転で扱う。定常性から必要な軸対称性を導く剛性定理と、二次元に縮約した非線形境界値問題の一意性が必要になる。この章では定理の結論を用いる。定常な事象の地平線を Killing 地平線として扱う剛性の議論と、軸対称な外部解の一意性の考え方は :cite[wald-gr] を参照する。この主張は一般の時間依存する外部に重力波や物質が存在しないという意味ではない。例えば二つのブラックホールが接近する過程の時空には、定常・連結という仮定を適用できない。

### 一意性の仮定を別々に確認する

同じ局所方程式でも、その解に許す漸近条件によって一意性の結論は変わる。例えば断熱的に平衡へ近づく過程では、時刻依存の放射が残っている間は定常性を仮定できない。また二つの極限荷電穴が静止する模型では連結性と非退化性の仮定が破れる。Kerr–Newman 族内では $m^2>a^2+q^2$ のとき二つの単純根を持ち、等号で重根となる。従って同じ族で連続な極限を取れることから、非退化を仮定した一意性定理をその極限の他の時空へ拡張することはできない。ある解を同定する手順は、電荷を含む方程式、外部の正則性、漸近端、地平線の性質を順に検査してから定理の対応する形を選ぶことである。

### 電磁真空以外へ結論を広げる前に調べること

:::example[宇宙項がある球対称解で漸近条件を確かめる]
局所で $f=1-2m/r-\Lambda r^2/3$ を持つ球対称な計量を考える。曲率の縮約は $R_{ab}=\Lambda g_{ab}$ なので、$\Lambda\ne0$ では真空の方程式の意味も変わる。例えば $\Lambda<0$ なら大きな $r$ で $f\simeq|\Lambda|r^2/3$ となり、無限遠で $g_{tt}\to-1$ にする定数倍の時刻変換は存在しない。また $R=4\Lambda$ が遠方でも残る。従って静的・球対称・局所的に真空であっても、漸近平坦という仮定をこの計量に適用できない。電荷・回転以外の仮定の確認が必要になる具体例である。
:::

### 無毛定理の連結性を外した定常解の例

四次元の電磁真空で、平坦な三座標 $\vb x$ 上の正の調和関数 $H$ により
$$
\dd s^2=-H^{-2}\dd t^2+H^2\dd\vb x^2,\qquad
A_t=\frac{H^{-1}-1}{\sqrt{4\pi G}}
$$
とする。$A_i=0$ であり、空間微分を $\partial_i$ とする。$\sqrt{-g}=H^2$ で、$H^2 F^{it}=\partial_iH/\sqrt{4\pi G}$ だから Maxwell 方程式は $\Delta H=0$ に等しい。Einstein 側を確かめるには三計量 $\gamma_{ij}=H^2\delta_{ij}$ と静的 lapse $N=H^{-1}$ に分ける。時間を含む接続は $\Gamma^t{}_{ti}=N^{-1}\partial_iN$ と $\Gamma^i{}_{tt}=N\gamma^{ij}\partial_jN$ であり、縮約すると $R_{tt}=N D^iD_iN$、$R_{ij}=R_{ij}[\gamma]-N^{-1}D_iD_jN$ となる。

$b_i=\partial_i\log H$ とおく。空間接続は $\Gamma^i{}_{jk}=\delta^i_jb_k+\delta^i_k b_j-\delta_{jk}b^i$ である。これを微分して積を取ると
$$
R_{ij}[\gamma]=-\partial_i b_j-\delta_{ij}\partial_k b_k+b_i b_j-\delta_{ij}b_kb_k,\qquad
N^{-1}D_iD_j N=-\partial_i b_j+3b_i b_j-\delta_{ij}b_kb_k.
$$
$\Delta H=0$ により $\partial_i b_i=-b_i b_i$ である。従って $R_{ij}=-2b_i b_j+\delta_{ij}b_kb_k$ となる。もう一方の式では $N^{-1}D^iD_iN=H^{-2}b_kb_k$ なので $R_{tt}=H^{-4}b_kb_k$ である。正規直交枠の電場は符号を除いて $E_i=H^{-1}b_i/\sqrt{4\pi G}$ であり、$8\pi GT_{\hat0\hat0}=H^{-2}b_kb_k$、$8\pi GT_{\hat i\hat j}=H^{-2}(\delta_{ij}b_kb_k-2b_i b_j)$ となる。従って $R=0$ と合わせて電磁真空の Einstein 方程式を満たす。

:::example[一つの極と二つ以上の極]
$H=1+m_1/|\vb x|$、$m_1>0$ なら面積半径は $r=|\vb x|+m_1$ であり、$g_{tt}=-(1-m_1/r)^2$ となる。従って極限 Reissner–Nordström の外部に一致する。調和関数を $H=1+\sum_i m_i/|\vb x-\vb x_i|$ とすれば、各極を除いた領域で同じ検算が成り立つ。各極に近づくと $H=m_i/|\vb x-\vb x_i|+\order{1}$ であり、固有距離は対数的に発散し、面積は $4\pi m_i^2$ に近づく。lapse の二乗も二次で零になるため、各極は非退化な単一球面の近傍とは異なる。こうした多中心解の地平線の滑らかな延長を確認するには角依存の高次項まで扱う必要があり、局所の方程式だけでは判定できない。ここで計算したのは、連結かつ非退化という一意性の仮定から離れた定常な電磁真空の外部である。
:::

## 演習問題

::::exercise{level=basic}
$q=3m/5$ の球対称荷電解の地平線を求めよ。

:::solution
$\sqrt{m^2-q^2}=4m/5$ であり、$r_+=9m/5,r_-=m/5$ である。
:::
::::

::::exercise{level=basic}
$a=0$ の Kerr の外部地平線面積を求めよ。

:::solution
$r_+=2m$ だから $4\pi(r_+^2+a^2)=16\pi m^2$ となる。
:::
::::

::::exercise{level=basic}
$a=4m/5$ の Kerr の赤道のエルゴ領域の半径範囲を求めよ。

:::solution
$r_+=8m/5$ であり、赤道では $\cos\theta=0$ より $r_{\rm e}=2m$ である。範囲は $8m/5<r<2m$ である。
:::
::::

::::exercise{level=basic}
赤道を $p_\theta=3$ で横切る測地線の $\mathcal Q$ を求めよ。

:::solution
角の式で $\cos\theta=0$ なので $\mathcal Q=p_\theta^2=9$ である。
:::
::::

::::exercise{level=advanced}
球対称荷電解の非退化な内外地平線での $f'$ を求めよ。

:::solution
$f=(r-r_+)(r-r_-)/r^2$ なので $f'(r_\pm)=(r_\pm-r_\mp)/r_\pm^2$ である。外では正、内では負となる。
:::
::::

::::exercise{level=advanced}
Kerr で $a>0$ とする。$L=0$ の未来時間的粒子が負の $E$ のまま地平線を越えられないことを示せ。

:::solution
地平線上では非零の未来光的 $\chi$ と未来時間的粒子の内積は厳密に負である。従って $-p\cdot\chi=E-\Omega_H L=E>0$ が必要である。
:::
::::

::::exercise{level=advanced}
赤道で定半径の Kerr 光軌道に必要な動径多項式の条件を説明せよ。

:::solution
:eqref[eq-gr-carter] で $\mu=\mathcal Q=0$ とし、右辺を $\mathcal R(r)$ とする。$\dot r=0$ から $\mathcal R=0$、動径 Hamilton 方程式で $\dot p_r=0$ から $\mathcal R'=0$ を要する。半径微分では $E,L$ を固定する。
:::
::::

::::exercise{level=advanced}
$a=0$ で角方程式から通常の角運動量二乗を復元せよ。

:::solution
$\mathcal Q=p_\theta^2+L^2\cot^2\theta$ なので、$\mathcal Q+L^2=p_\theta^2+L^2/\sin^2\theta=L_{\rm total}^2$ となる。
:::
::::

::::exercise{level=advanced}
Reissner–Nordström の中性光の円軌道半径を $m,q$ で求めよ。$0<q<m$ とし、大きい方の根が外側地平線より大きいことを示せ。

:::solution
光のポテンシャル $L^2f/r^2$ を微分すると $-2L^2(r^2-3mr+2q^2)/r^5$ となる。従って根は $(3m\pm\sqrt{9m^2-8q^2})/2$ であり、外側根を $r_{\rm ph}$ とする。地平線 $r_+$ で分子の二次式は、$q^2=2mr_+-r_+^2$ を使い $r_+(m-r_+)<0$ となる。上に開く二次式が負であるから $r_+$ は二根の間にあり、$r_{\rm ph}>r_+$ である。動径ポテンシャルの微分は大きい根で正から負に変わるため不安定である。$q\to0$ の外側根は $3m$、$q\to m$ で $2m$ に近づく。内側の根を同じ外部の観測者が周回する軌道と解釈することはできない。
:::
::::

## 確認テスト

::::quiz
Kerr の赤道のエルゴ領域で成立する命題を選べ。

- [x] $\partial_t$ は空間的である。
- [ ] 未来因果曲線はすべて外へ出られない。
- [ ] $E=-p_t$ は未来粒子の局所エネルギーと必ず一致する。

:::explanation
領域では $g_{tt}>0$ である。事象の地平線はその内側にあり、局所エネルギーには時間的な観測者との縮約を使うため $-p_t$ と同一視できない。
:::
::::

::::quiz{answer="1.6" tol="1e-9"}
$a/m=4/5$ の Kerr の外側地平線で $r_+/m$ を答えよ。

:::explanation
$\Delta=0$ の大きい根は $r_+/m=1+\sqrt{1-(4/5)^2}=1+3/5=1.6$ である。
:::
::::

::::quiz
分離定数 $\mathcal Q$ を求めるのに用いる量を選べ。

- [x] $\theta$ と $p_\theta$、および $E,L,\mu$
- [ ] $r$ だけ
- [ ] 事象の地平線面積だけ

:::explanation
角運動の式を解くと $\mathcal Q=p_\theta^2+\cos^2\theta[a^2(\mu^2-E^2)+L^2/\sin^2\theta]$ である。測地線の運動量を指定する必要がある。
:::
::::

::::quiz{answer="0.25" tol="1e-9"}
$a/m=4/5$ の Kerr の $m\Omega_H$ を答えよ。

:::explanation
$r_+/m=8/5$ なので $m\Omega_H=(4/5)/[(8/5)^2+(4/5)^2]=1/4$ である。
:::
::::

::::quiz{answer="1.36" tol="1e-7"}
$a=3m/5,q=4m/5$ の Kerr–Newman で $A_H/(4\pi m^2)$ を求めよ。

:::explanation
$a^2+q^2=m^2$ なので外側根は $r_+=m$ である。$A_H=4\pi(r_+^2+a^2)$ であり、規格化した値は $1+9/25=34/25$ である。
:::
::::

---
title: "相対論的星"
description: "静水圧平衡式から相対論的星の内部を求め、密度・縮退圧・崩壊が質量や半径を制限する仕組みを調べる。"
notation:
  - symbol: "M(r)"
    meaning: "半径 r の内側の重力質量"
  - symbol: "\\rho,p"
    meaning: "全エネルギー密度と圧力"
  - symbol: "R"
    meaning: "圧力が零となる表面"
  - symbol: "n"
    meaning: "粒子の数密度"
---

## この章の目標

- 完全流体と真空の境界で計量を接続する。
- TOV 方程式を導き、内部解を中心から積分する。
- 一様密度解から中心圧を求め、Buchdahl の不等式を証明する。
- 縮退星の質量尺度と塵の崩壊の違いを説明する。

## 静的球対称な完全流体

静止流体の四元速度を $u=e^{-\Phi}\partial_t$、面積半径を $r$ とする。:eqref[eq-gr-perfect-fluid] により混合成分は $\diag(-\rho,p,p,p)$ であり、$T^t{}_r=0$ で熱流はない。計量を $g_{tt}=-e^{2\Phi}$、$g_{rr}=(1-2GM(r)/r)^{-1}$ と置く。$M(r)$ は面積半径から決めた質量関数であり、密度を固有体積で積分した量とは区別する。例えば静止の粒子数は $4\pi\int n r^2(1-2GM/r)^{-1/2}\dd r$ であり、$M'$ の式にはこの根号がない。中心の滑らかさには $M(0)=0$、有限の $\rho,p$ が必要である。表面 $r=R$ では $p=0$ とし、真空 Schwarzschild 解に接ぐ。時刻を外部で規格化すると $e^{2\Phi(R)}=1-2GM(R)/R$ となる。薄い表面層を置かない場合は表面の外的曲率も連続であり、以下の方程式では $M$ と $p$ の境界値でその条件が満たされる。表面の密度は有限の跳びを許す。

### 表面の外的曲率と圧力

表面 $r=R$ の外向き単位法線は $n=\sqrt f\,\partial_r$ である。法線に沿う計量の Lie 微分の半分を接方向に制限すると、表面の正規直交枠で $K_{\hat t\hat t}=-\sqrt f\,\Phi'$、$K_{\hat\theta\hat\theta}=K_{\hat\phi\hat\phi}=\sqrt f/R$ となる。角成分が連続なら $f$ は連続であり、外部の質量は $M(R)$ である。一方、時間成分が連続なら $\Phi'$ が連続である。動径の Einstein 方程式から $\Phi'$ の内外の差は $4\pi GRp(R)/f(R)$ となるため、真空と滑らかにつなぐ条件は $p(R)=0$ である。$M'$ はこの二つの外的曲率成分に含まれないので、有限の表面密度はこれらの条件に矛盾しない。

:::example[静止時の局所エネルギーと無限遠のエネルギー]
表面の静止時計の lapse は $e^{\Phi(R)}=\sqrt{1-2GM/R}$ である。例えば $2GM/R=5/9$ なら表面で振動数 $\nu_R$ の光は無限遠で $2\nu_R/3$ となり、赤方偏移は $z=1/2$ である。搭載した時計を表面に静止させるために必要な加速度の大きさは $\sqrt f\,\Phi'=GM/(R^2\sqrt f)$ で、$GM/R^2$ の $3/2$ 倍となる。同じ lapse の因子が時間の規格化と加速度の両方に現れる。
:::

:::example[二関数の球対称計量から内部の条件を作る]
星の静的領域で $P=e^\Phi>0$、$B=1-2GM/r>0$ とする。外部の真空では二つの関数に関係があるが、内部でその関係を先に課すと独立な圧力方程式を失うため、ここでは別の関数とする。 座標の順序を $(t,r,\theta,\phi)$ と固定する。成分行列と逆行列は次の組になる。

$$
\begin{aligned}
(g_{ab})&=\left[\begin{matrix}- P^{2}{\left(r \right)} & 0 & 0 & 0\\0 & \frac{1}{B{\left(r \right)}} & 0 & 0\\0 & 0 & r^{2} & 0\\0 & 0 & 0 & r^{2} \sin^{2}{\left(\theta \right)}\end{matrix}\right]\\
(g^{ab})&=\left[\begin{matrix}- \frac{1}{P^{2}{\left(r \right)}} & 0 & 0 & 0\\0 & B{\left(r \right)} & 0 & 0\\0 & 0 & \frac{1}{r^{2}} & 0\\0 & 0 & 0 & \frac{1}{r^{2} \sin^{2}{\left(\theta \right)}}\end{matrix}\right]\\
\det g&=- \frac{r^{4} P^{2}{\left(r \right)} \sin^{2}{\left(\theta \right)}}{B{\left(r \right)}}
\end{aligned}
$$

逆行列を掛けて単位行列を得ることを確かめてから、:eqref[eq-gr-christoffel] を使う。同じ下添字の順序を二度計算せず、$b\le c$ の組を選ぶ。中央の列には、座標に依存する計量成分の微分だけを残した。最後の列は微分を実行した値であり、表に出ない成分は添字交換で得られるものを除いて零である。

| 成分 | 微分への代入 | 結果 |
|---|---|---|
| $\Gamma^{t}_{t r}$ | $- \frac{1}{2 P^{2}{\left(r \right)}}\left(\partial_{r}\left(- P^{2}{\left(r \right)}\right)\right)$ | $\frac{P'{\left(r \right)}}{P{\left(r \right)}}$ |
| $\Gamma^{r}_{t t}$ | $\frac{B{\left(r \right)}}{2}\left(-\partial_{r}\left(- P^{2}{\left(r \right)}\right)\right)$ | $B{\left(r \right)} P{\left(r \right)} P'{\left(r \right)}$ |
| $\Gamma^{r}_{r r}$ | $\frac{B{\left(r \right)}}{2}\left(\partial_{r}\left(\frac{1}{B{\left(r \right)}}\right)+\partial_{r}\left(\frac{1}{B{\left(r \right)}}\right)-\partial_{r}\left(\frac{1}{B{\left(r \right)}}\right)\right)$ | $- \frac{B'{\left(r \right)}}{2 B{\left(r \right)}}$ |
| $\Gamma^{r}_{\theta \theta}$ | $\frac{B{\left(r \right)}}{2}\left(-\partial_{r}\left(r^{2}\right)\right)$ | $- r B{\left(r \right)}$ |
| $\Gamma^{r}_{\phi \phi}$ | $\frac{B{\left(r \right)}}{2}\left(-\partial_{r}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $- r B{\left(r \right)} \sin^{2}{\left(\theta \right)}$ |
| $\Gamma^{\theta}_{r \theta}$ | $\frac{1}{2 r^{2}}\left(\partial_{r}\left(r^{2}\right)\right)$ | $\frac{1}{r}$ |
| $\Gamma^{\theta}_{\phi \phi}$ | $\frac{1}{2 r^{2}}\left(-\partial_{\theta}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $- \frac{\sin{\left(2 \theta \right)}}{2}$ |
| $\Gamma^{\phi}_{r \phi}$ | $\frac{1}{2 r^{2} \sin^{2}{\left(\theta \right)}}\left(\partial_{r}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{1}{r}$ |
| $\Gamma^{\phi}_{\theta \phi}$ | $\frac{1}{2 r^{2} \sin^{2}{\left(\theta \right)}}\left(\partial_{\theta}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{1}{\tan{\left(\theta \right)}}$ |

曲率の計算では :eqref[eq-gr-curvature-components] の微分差と積差を分ける。第二列は $\partial_c\Gamma^a{}_{db}$、第三列は $-\partial_d\Gamma^a{}_{cb}$ であり、第四列は接続の積を符号付きで足したものである。三列を足して最後の列を得る。

| 成分 | 第一の微分 | 第二の微分 | 接続の積 | 結果 |
|---|---|---|---|---|
| $R^{t}_{r t r}$ | $0$ | $- \frac{P''{\left(r \right)}}{P{\left(r \right)}} + \frac{\left(P'{\left(r \right)}\right)^{2}}{P^{2}{\left(r \right)}}$ | $- \frac{\left(P'{\left(r \right)}\right)^{2}}{P^{2}{\left(r \right)}} - \frac{B'{\left(r \right)} P'{\left(r \right)}}{2 B{\left(r \right)} P{\left(r \right)}}$ | $- \frac{2 B{\left(r \right)} P''{\left(r \right)} + B'{\left(r \right)} P'{\left(r \right)}}{2 B{\left(r \right)} P{\left(r \right)}}$ |
| $R^{t}_{\theta t \theta}$ | $0$ | $0$ | $- \frac{r B{\left(r \right)} P'{\left(r \right)}}{P{\left(r \right)}}$ | $- \frac{r B{\left(r \right)} P'{\left(r \right)}}{P{\left(r \right)}}$ |
| $R^{t}_{\phi t \phi}$ | $0$ | $0$ | $- \frac{r B{\left(r \right)} \sin^{2}{\left(\theta \right)} P'{\left(r \right)}}{P{\left(r \right)}}$ | $- \frac{r B{\left(r \right)} \sin^{2}{\left(\theta \right)} P'{\left(r \right)}}{P{\left(r \right)}}$ |
| $R^{r}_{\theta r \theta}$ | $- r B'{\left(r \right)} - B{\left(r \right)}$ | $0$ | $\frac{r B'{\left(r \right)}}{2} + B{\left(r \right)}$ | $- \frac{r B'{\left(r \right)}}{2}$ |
| $R^{r}_{\phi r \phi}$ | $- r \sin^{2}{\left(\theta \right)} B'{\left(r \right)} - B{\left(r \right)} \sin^{2}{\left(\theta \right)}$ | $0$ | $\frac{r \sin^{2}{\left(\theta \right)} B'{\left(r \right)}}{2} + B{\left(r \right)} \sin^{2}{\left(\theta \right)}$ | $- \frac{r \sin^{2}{\left(\theta \right)} B'{\left(r \right)}}{2}$ |
| $R^{\theta}_{\phi \theta \phi}$ | $- \cos{\left(2 \theta \right)}$ | $0$ | $- B{\left(r \right)} \sin^{2}{\left(\theta \right)} + \frac{\sin{\left(2 \theta \right)}}{2 \tan{\left(\theta \right)}}$ | $- \left(B{\left(r \right)} - 1\right) \sin^{2}{\left(\theta \right)}$ |

第一添字を計量で下ろすと $R_{abcd}=g_{ak}R^k{}_{bcd}$ となる。対の反対称性と交換対称性を用いて、非零の独立な組を下表に示す。

| 成分 | 結果 |
|---|---|
| $R_{t r t r}$ | $\frac{\left(2 B{\left(r \right)} P''{\left(r \right)} + B'{\left(r \right)} P'{\left(r \right)}\right) P{\left(r \right)}}{2 B{\left(r \right)}}$ |
| $R_{t \theta t \theta}$ | $r B{\left(r \right)} P{\left(r \right)} P'{\left(r \right)}$ |
| $R_{t \phi t \phi}$ | $r B{\left(r \right)} P{\left(r \right)} \sin^{2}{\left(\theta \right)} P'{\left(r \right)}$ |
| $R_{r \theta r \theta}$ | $- \frac{r B'{\left(r \right)}}{2 B{\left(r \right)}}$ |
| $R_{r \phi r \phi}$ | $- \frac{r \sin^{2}{\left(\theta \right)} B'{\left(r \right)}}{2 B{\left(r \right)}}$ |
| $R_{\theta \phi \theta \phi}$ | $- r^{2} \left(B{\left(r \right)} - 1\right) \sin^{2}{\left(\theta \right)}$ |

Ricci テンソルは第一・第三添字を縮約して得る。曲率が残る場合も相殺が起こるので、各項を落とさずに加える。

| 成分 | 縮約の項 | 代入 | 結果 |
|---|---|---|---|
| $R_{t t}$ | $R^{r}_{t r t}+R^{\theta}_{t \theta t}+R^{\phi}_{t \phi t}$ | $\frac{\left(2 B{\left(r \right)} P''{\left(r \right)} + B'{\left(r \right)} P'{\left(r \right)}\right) P{\left(r \right)}}{2} + \frac{B{\left(r \right)} P{\left(r \right)} P'{\left(r \right)}}{r} + \frac{B{\left(r \right)} P{\left(r \right)} P'{\left(r \right)}}{r}$ | $\frac{\left(2 r B{\left(r \right)} P''{\left(r \right)} + r B'{\left(r \right)} P'{\left(r \right)} + 4 B{\left(r \right)} P'{\left(r \right)}\right) P{\left(r \right)}}{2 r}$ |
| $R_{r r}$ | $R^{t}_{r t r}+R^{\theta}_{r \theta r}+R^{\phi}_{r \phi r}$ | $- \frac{2 B{\left(r \right)} P''{\left(r \right)} + B'{\left(r \right)} P'{\left(r \right)}}{2 B{\left(r \right)} P{\left(r \right)}} - \frac{B'{\left(r \right)}}{2 r B{\left(r \right)}} - \frac{B'{\left(r \right)}}{2 r B{\left(r \right)}}$ | $- \frac{2 r B{\left(r \right)} P''{\left(r \right)} + r B'{\left(r \right)} P'{\left(r \right)} + 2 P{\left(r \right)} B'{\left(r \right)}}{2 r B{\left(r \right)} P{\left(r \right)}}$ |
| $R_{\theta \theta}$ | $R^{t}_{\theta t \theta}+R^{r}_{\theta r \theta}+R^{\phi}_{\theta \phi \theta}$ | $- \frac{r B{\left(r \right)} P'{\left(r \right)}}{P{\left(r \right)}} - \frac{r B'{\left(r \right)}}{2} - \left(B{\left(r \right)} - 1\right)$ | $- \frac{2 r B{\left(r \right)} P'{\left(r \right)} + r P{\left(r \right)} B'{\left(r \right)} + 2 B{\left(r \right)} P{\left(r \right)} - 2 P{\left(r \right)}}{2 P{\left(r \right)}}$ |
| $R_{\phi \phi}$ | $R^{t}_{\phi t \phi}+R^{r}_{\phi r \phi}+R^{\theta}_{\phi \theta \phi}$ | $- \frac{r B{\left(r \right)} \sin^{2}{\left(\theta \right)} P'{\left(r \right)}}{P{\left(r \right)}} - \frac{r \sin^{2}{\left(\theta \right)} B'{\left(r \right)}}{2} - \left(B{\left(r \right)} - 1\right) \sin^{2}{\left(\theta \right)}$ | $- \frac{\left(2 r B{\left(r \right)} P'{\left(r \right)} + r P{\left(r \right)} B'{\left(r \right)} + 2 B{\left(r \right)} P{\left(r \right)} - 2 P{\left(r \right)}\right) \sin^{2}{\left(\theta \right)}}{2 P{\left(r \right)}}$ |

逆計量と縮約すると Ricci スカラーとなる。さらに $G_{ab}=R_{ab}-g_{ab}R/2$ を成分ごとに計算する。

$$
\begin{aligned}
R&=\sum_{a,b}g^{ab}R_{ab}\\
&=- \frac{2 r B{\left(r \right)} P''{\left(r \right)} + r B'{\left(r \right)} P'{\left(r \right)} + 4 B{\left(r \right)} P'{\left(r \right)}}{2 r P{\left(r \right)}} - \frac{2 r B{\left(r \right)} P''{\left(r \right)} + r B'{\left(r \right)} P'{\left(r \right)} + 2 P{\left(r \right)} B'{\left(r \right)}}{2 r P{\left(r \right)}} - \frac{2 r B{\left(r \right)} P'{\left(r \right)} + r P{\left(r \right)} B'{\left(r \right)} + 2 B{\left(r \right)} P{\left(r \right)} - 2 P{\left(r \right)}}{2 r^{2} P{\left(r \right)}} - \frac{2 r B{\left(r \right)} P'{\left(r \right)} + r P{\left(r \right)} B'{\left(r \right)} + 2 B{\left(r \right)} P{\left(r \right)} - 2 P{\left(r \right)}}{2 r^{2} P{\left(r \right)}}\\
&=- \frac{2 r^{2} B{\left(r \right)} P''{\left(r \right)} + r^{2} B'{\left(r \right)} P'{\left(r \right)} + 4 r B{\left(r \right)} P'{\left(r \right)} + 2 r P{\left(r \right)} B'{\left(r \right)} + 2 B{\left(r \right)} P{\left(r \right)} - 2 P{\left(r \right)}}{r^{2} P{\left(r \right)}}
\end{aligned}
$$

非零の Einstein テンソル成分は次のとおりである。

| 成分 | Ricci 項 | 跡の項 | 差 |
|---|---|---|---|
| $G_{t t}$ | $\frac{\left(2 r B{\left(r \right)} P''{\left(r \right)} + r B'{\left(r \right)} P'{\left(r \right)} + 4 B{\left(r \right)} P'{\left(r \right)}\right) P{\left(r \right)}}{2 r}$ | $\frac{\left(2 r^{2} B{\left(r \right)} P''{\left(r \right)} + r^{2} B'{\left(r \right)} P'{\left(r \right)} + 4 r B{\left(r \right)} P'{\left(r \right)} + 2 r P{\left(r \right)} B'{\left(r \right)} + 2 B{\left(r \right)} P{\left(r \right)} - 2 P{\left(r \right)}\right) P{\left(r \right)}}{2 r^{2}}$ | $- \frac{\left(r B'{\left(r \right)} + B{\left(r \right)} - 1\right) P^{2}{\left(r \right)}}{r^{2}}$ |
| $G_{r r}$ | $- \frac{2 r B{\left(r \right)} P''{\left(r \right)} + r B'{\left(r \right)} P'{\left(r \right)} + 2 P{\left(r \right)} B'{\left(r \right)}}{2 r B{\left(r \right)} P{\left(r \right)}}$ | $- \frac{2 r^{2} B{\left(r \right)} P''{\left(r \right)} + r^{2} B'{\left(r \right)} P'{\left(r \right)} + 4 r B{\left(r \right)} P'{\left(r \right)} + 2 r P{\left(r \right)} B'{\left(r \right)} + 2 B{\left(r \right)} P{\left(r \right)} - 2 P{\left(r \right)}}{2 r^{2} B{\left(r \right)} P{\left(r \right)}}$ | $\frac{2 r B{\left(r \right)} P'{\left(r \right)} + B{\left(r \right)} P{\left(r \right)} - P{\left(r \right)}}{r^{2} B{\left(r \right)} P{\left(r \right)}}$ |
| $G_{\theta \theta}$ | $- \frac{2 r B{\left(r \right)} P'{\left(r \right)} + r P{\left(r \right)} B'{\left(r \right)} + 2 B{\left(r \right)} P{\left(r \right)} - 2 P{\left(r \right)}}{2 P{\left(r \right)}}$ | $- \frac{2 r^{2} B{\left(r \right)} P''{\left(r \right)} + r^{2} B'{\left(r \right)} P'{\left(r \right)} + 4 r B{\left(r \right)} P'{\left(r \right)} + 2 r P{\left(r \right)} B'{\left(r \right)} + 2 B{\left(r \right)} P{\left(r \right)} - 2 P{\left(r \right)}}{2 P{\left(r \right)}}$ | $\frac{r \left(2 r B{\left(r \right)} P''{\left(r \right)} + r B'{\left(r \right)} P'{\left(r \right)} + 2 B{\left(r \right)} P'{\left(r \right)} + P{\left(r \right)} B'{\left(r \right)}\right)}{2 P{\left(r \right)}}$ |
| $G_{\phi \phi}$ | $- \frac{\left(2 r B{\left(r \right)} P'{\left(r \right)} + r P{\left(r \right)} B'{\left(r \right)} + 2 B{\left(r \right)} P{\left(r \right)} - 2 P{\left(r \right)}\right) \sin^{2}{\left(\theta \right)}}{2 P{\left(r \right)}}$ | $- \frac{\left(2 r^{2} B{\left(r \right)} P''{\left(r \right)} + r^{2} B'{\left(r \right)} P'{\left(r \right)} + 4 r B{\left(r \right)} P'{\left(r \right)} + 2 r P{\left(r \right)} B'{\left(r \right)} + 2 B{\left(r \right)} P{\left(r \right)} - 2 P{\left(r \right)}\right) \sin^{2}{\left(\theta \right)}}{2 P{\left(r \right)}}$ | $\frac{r \left(2 r B{\left(r \right)} P''{\left(r \right)} + r B'{\left(r \right)} P'{\left(r \right)} + 2 B{\left(r \right)} P'{\left(r \right)} + P{\left(r \right)} B'{\left(r \right)}\right) \sin^{2}{\left(\theta \right)}}{2 P{\left(r \right)}}$ |

静止枠に直して $G^t{}_t=-8\pi G\rho$、$G^r{}_r=8\pi Gp$、$G^\theta{}_\theta=G^\phi{}_\phi=8\pi Gp$ と読む。最初の式は $B$ と $B'$ だけを含む。第二式には $P'/P$ が現れる。角の式にはさらに $P''$ があるため、密度と動径圧で決めた解に対して角の式と整合するかを、保存則と等方圧を使って確認できる。
:::

### 固有体積で量った質量と重力質量

球殻ごとの固有エネルギーを同じ静的切片で積分した量を $M_{\rm prop}=\int_0^R4\pi r^2\rho/\sqrt f\,\dd r$ とする。$\rho\ge0$ なら $0<f\le1$ であるから $M_{\rm prop}\ge M(R)$ である。重力質量との差には空間の曲率が入る。一方、粒子数 $N$ に粒子の静止質量 $m_b$ を掛けた $m_bN$ には、熱や粒子間の相互作用によるエネルギーを含めない。従って $M_{\rm prop}-M(R)$ と $m_bN-M(R)$ の間には内部エネルギー分の差がある。

弱い場で $1/\sqrt f=1+GM(r)/r+\order{(GM(r)/r)^2}$ と展開すると
$$
M_{\rm prop}-M(R)=\int_0^R\frac{GM(r)}r\,\dd M(r)+\cdots
$$
になる。Newton の自己重力エネルギーは半径 $r$ の殻を外から積むときの仕事 $-GM(r)\dd M(r)/r$ を足したものであり、この極限の差はその絶対値に一致する。例えば一様密度では $M(r)=M(R)(r/R)^3$ なので $3GM(R)^2/(5R)$ である。

:::example[一様密度星の固有エネルギーを全次数で計算する]
$C=2GM(R)/R$、$x=r/R$ とすれば $M_{\rm prop}/M(R)=3\int_0^1x^2/\sqrt{1-Cx^2}\,\dd x$ である。$y=\arcsin(\sqrt C\,x)$ を使うと、積分は $3C^{-3/2}\int_0^{\arcsin\sqrt C}\sin^2y\,\dd y$ になり、
$$
\frac{M_{\rm prop}}{M(R)}=\frac{3}{2C^{3/2}}
\left(\arcsin\sqrt C-\sqrt{C(1-C)}\right)
$$
を得る。$C\ll1$ では $1+3C/10+9C^2/56+\order{C^3}$ であり、一次の差は先の自己重力エネルギーと一致する。例えば $C=3/4$ なら $M_{\rm prop}/M(R)=4\pi/(3\sqrt3)-1$ で、約 $1.42$ である。
:::

## Tolman–Oppenheimer–Volkoff 方程式

:eqref[eq-gr-spherical-field] に $f=1-2GM/r$ を代入すると、時間成分から $M'=4\pi r^2\rho$、動径成分から $\Phi'=G(M+4\pi r^3p)/(r^2f)$ である。静止流体の加速度の共変動径成分は $\Phi'$ なので、:eqref[eq-gr-fluid-conservation] が $p'=-(\rho+p)\Phi'$ を与える。まとめると :term[Tolman–Oppenheimer–Volkoff 方程式]{en="Tolman–Oppenheimer–Volkoff equation" yomi="とるまんおっぺんはいまーヴぉるこふほうていしき"} は

$$
\dv{M}{r}=4\pi r^2\rho,\qquad
\dv{p}{r}=-\frac{G(\rho+p)(M+4\pi r^3p)}{r^2(1-2GM/r)}.
\label{eq-gr-tov}
$$

圧力は慣性側 $\rho+p$ と場の源側 $M+4\pi r^3p$ の両方に入る。弱重力 $2GM/r\ll1$ と低圧 $p\ll\rho$、$4\pi r^3p\ll M$ のもとで $p'\simeq-G\rho M/r^2$ となる。密度が外向きに非増加であれば、後者の圧力条件は $p\ll\rho$ から従う。状態方程式 $\rho(p)$ と中心圧 $p_c$ を指定すれば外向きに積分できる。数値計算では原点の分母を直接計算せず、$M=(4\pi/3)\rho_c r^3+\order{r^5}$、$p=p_c-2\pi G(\rho_c+p_c)(\rho_c/3+p_c)r^2+\order{r^4}$ を小半径で使う。例えば $p_c=\rho_c/3$ では二次項の係数が $-16\pi G\rho_c^2/9$ となる。

### 中心展開で初期条件をそろえる

中心で $\rho(p)$ が二回微分可能なら、$b=(\dd\rho/\dd p)_c$ とし $p=p_c+p_2r^2+p_4r^4+\cdots$ と置ける。密度展開 $\rho_c+bp_2r^2+\cdots$ を質量の式に入れると $M=M_3r^3+M_5r^5+\cdots$、$M_3=4\pi\rho_c/3,M_5=4\pi bp_2/5$ である。圧力式で $1/f=1+2GM_3r^2+\cdots$ とし、係数を比較すると
$$
\begin{aligned}
2p_2&=-G(\rho_c+p_c)(M_3+4\pi p_c),\\
4p_4&=-G\{(1+b)p_2(M_3+4\pi p_c)
+(\rho_c+p_c)[M_5+4\pi p_2+2GM_3(M_3+4\pi p_c)]\}.
\end{aligned}
$$
従って状態方程式の局所微分が圧力に現れるのは四次からである。この展開で初期点の半径 $\epsilon$ を小さくしても、積分の刻みを粗くすれば表面位置の誤差は残る。逆に刻みを細かくしても中心展開を粗く切れば初期値の誤差が残るので、二つの尺度を別に変えて結果の収束を確認する。

:::example[表面での圧力の傾き]
自己束縛模型で $\rho(p=0)=\rho_s>0$ とする。表面から内側へ小距離 $\delta$ で $p(R-\delta)=G\rho_s M(R)\delta/[R^2(1-2GM(R)/R)]+\order{\delta^2}$ である。例えば $2GM(R)/R=1/2$ なら $p(R-\delta)/\rho_s=\delta/(2R)+\order{(\delta/R)^2}$ となる。これに対して表面で密度も零になる模型では同じ代入の線形項が消え、表面の位置を補間する際に傾きが小さくなり得る。
:::

### 角成分の方程式が満たされることを確かめる

計量から作った $S^a{}_b=G^a{}_b-8\pi G T^a{}_b$ に対し、Bianchi 恒等式と保存則から $\nabla_aS^a{}_b=0$ である。球対称・静的で $S^a{}_b=\diag(S_t,S_r,S_\Omega,S_\Omega)$ とすると、
$$
0=\nabla_a S^a{}_r=S'_r+\Phi'(S_r-S_t)+\frac2r(S_r-S_\Omega).
$$
この式は、斜めの縮約で $\Gamma^a{}_{ar}S_r-\Gamma^a{}_{ar}S_a$ が入ることから得られる。$tt$ と $rr$ の Einstein 方程式を解いて $S_t=S_r=0$ とすれば、$r>0$ で $S_\Omega=0$ が従う。中心は正則解の連続性で扱う。従って保存則から導いた圧力の式と、独立な二つの Einstein 方程式で、角成分まで満たす内部解になる。

### 中心の正則な初期値を積分方程式で与える

$M(r)=4\pi\int_0^r s^2\rho(p(s))\,\dd s$ を圧力方程式に代入し、$p(r)=p_c+\int_0^r p'(s)\dd s$ と書く。状態方程式 $\rho(p)$ が $p_c$ の近傍で有界かつ Lipschitz 連続であれば、この形から正則な局所解が得られる。実際、$M(r)/r^3=4\pi\int_0^1u^2\rho(p(ru))\,\dd u$ は有界で、二つの圧力候補の差に対して一様ノルムで Lipschitz 連続である。圧力方程式は $p'=-rG(\rho+p)(M/r^3+4\pi p)/(1-2Gr^2 M/r^3)$ であるから、$r\le\epsilon$ を十分小さく取ると分母を $1/2$ 以上にできる。

この範囲では圧力の積分作用素は、$p_c$ の近くの閉球を自分自身に写す。理由は $|p(r)-p_c|\le K\epsilon^2/2$ と見積もれるためである。また二つの候補への作用の差は $L\epsilon^2/2$ と候補の一様ノルム差の積以下である。$\epsilon$ をさらに小さくすれば縮小写像になり、一意な固定点が存在する。こうして得た圧力は積分方程式から微分可能になる。中心から離れ、$p$ が状態方程式の範囲にあり、$f>0$ の区間では通常の常微分方程式として延長する。この局所議論だけでは圧力が有限半径で零になるとは言えない。例えば $p=w\rho$ のように $\rho+p$ が $p$ に比例する模型では、正則な有限半径の点まで係数が有界なら $\log p$ の微分も有界であり、そこで突然 $p$ は零にならない。

## 一様密度星と Buchdahl の限界

### 一定密度の解

$\rho=\rho_0$、$m=GM(R)$、$y(r)=\sqrt{1-2mr^2/R^3}$、$y_R=y(R)$ とする。$\Phi'$ と $p'$ から $(\rho_0+p)e^\Phi=\rho_0 y_R$ である。$w=p/\rho_0$ とすると :eqref[eq-gr-tov] が $w'/[(1+w)(1+3w)]=y'/(2y)$ になる。積分で $\log((1+3w)/(1+w))=\log(y/y_R)$ を得る。従って

$$
p(r)=\rho_0\frac{y(r)-y_R}{3y_R-y(r)},\qquad
e^{\Phi(r)}=\tfrac12(3y_R-y(r)).
\label{eq-gr-uniform-star}
$$

中心の有限な正圧には $3y_R>1$、すなわち $2m/R<8/9$ が必要である。$2m/R=3/4$ なら $y_R=1/2$ で $p_c=\rho_0$ となり、Newton 近似の低圧条件が破れる。圧力が変わっても密度が一定の模型では圧縮率が零で、音速が無限大となるため、因果的な物質の状態方程式としては使わない。

### 単調な密度の上限

一般の正則な等方圧の星で $\rho\ge0$ が外向き非増加、$f>0$ と仮定する。$m(r)=GM(r)$、$y=\sqrt f$、$\zeta=e^\Phi$ と置く。角成分は接続を代入して $G^\theta{}_\theta=f\zeta''/\zeta+(f'/2+f/r)\zeta'/\zeta+f'/(2r)$、動径成分は $G^r{}_r=(f-1)/r^2+2f\zeta'/(r\zeta)$ となる。両者を等置して $\zeta'$ の項をまとめると

$$
y\dv{}{r}\left(\frac{y\zeta'}r\right)=\zeta\dv{}{r}\left(\frac{m(r)}{r^3}\right).
$$

平均密度は局所密度以上なので右辺は非正であり、$y\zeta'/r$ は外向き非増加である。表面値は $m(R)/R^3$ だから $\zeta_R-\zeta(0)\ge(m(R)/R^3)\int_0^R r/y(r)\,\dd r$ となる。さらに $m(r)/r^3\ge m(R)/R^3$ により積分は一様密度の値以上で、右辺は $(1-y_R)/2$ 以上となる。従って $\zeta(0)\le(3y_R-1)/2$ である。正則な中心は $\zeta(0)>0$ なので $2GM(R)/R<8/9$ を得る。これを :term[Buchdahl の限界]{en="Buchdahl bound" yomi="ぶふだーるのげんかい"} という。圧力に異方性を許すと角成分と動径成分の等置が成り立たない。

### 一様密度解の固有体積

表面のコンパクトネスを $C=2GM/R$ とする。一様密度星の固有エネルギー積分は
$$
M_{\rm prop}=4\pi\rho_0 R^3\int_0^1\frac{x^2\dd x}{\sqrt{1-Cx^2}}
=\frac{2\pi\rho_0 R^3}{C^{3/2}}\left[\arcsin\sqrt C-\sqrt{C(1-C)}\right]
$$
である。置換 $x=\sin u/\sqrt C$ により積分は $C^{-3/2}\int_0^{\arcsin\sqrt C}\sin^2u\,\dd u$ となる。$M=4\pi\rho_0R^3/3$ との比を小さな $C$ で展開すると $M_{\rm prop}/M=1+3C/10+\order{C^2}$ である。従って先頭の差 $M_{\rm prop}-M$ は $3GM^2/(5R)$ であり、Newton の球の結合エネルギーの絶対値と一致する。ここで $M_{\rm prop}$ は全エネルギー密度で作った量であり、流体の静止粒子質量だけを足した場合には内部エネルギーも区別する。

:::example[中心圧の発散を測る]
$C=8/9-\epsilon$ とすると、$3\sqrt{1-C}-1=(9/2)\epsilon+\order{\epsilon^2}$ である。中心圧の分子 $1-\sqrt{1-C}$ は $2/3+\order\epsilon$ なので、$p_c/\rho_0$ は $4/(27\epsilon)$ のように発散する。表面の redshift は $z=(1-C)^{-1/2}-1$ で、同じ極限でも $2$ に留まる。発散する量と有限な量を別の観測で区別できる。
:::

:::example[一様密度星で弱い場から高いコンパクトネスまで]
$C=2GM/R$ とする。表面 lapse は $y_R=\sqrt{1-C}$、中心 lapse は $\zeta_c=(3y_R-1)/2$ である。表面圧は零なので、内部の静的な流体で $(\rho_0+p)\zeta=\rho_0 y_R$ から中央圧を復元する。表面から無限遠への光の赤方偏移は $y_R^{-1}-1$ である。

ここでは$C$ は表面のコンパクトネスであるとする。計算に用いる式を先に並べる。

$$
\begin{aligned}
y_R&=\sqrt{1 - C}\\
\zeta_c&=\frac{3 \sqrt{1 - C}}{2} - \frac{1}{2}\\
p_c/\rho_0&=\frac{1 - \sqrt{1 - C}}{3 \sqrt{1 - C} - 1}\\
z_R&=-1 + \frac{1}{\sqrt{1 - C}}
\end{aligned}
$$

次の条件を比較する。入力を区別してから計算し、定数の規格化をそろえた値を表に示す。

| $C$ | $y_R$ | $\zeta_c$ | $p_c/\rho_0$ | $z_R$ |
|---|---|---|---|---|
| $\frac{19}{100}$ | $\frac{9}{10}$ | $\frac{17}{20}$ | $\frac{1}{17}$ | $\frac{1}{9}$ |
| $\frac{5}{9}$ | $\frac{2}{3}$ | $\frac{1}{2}$ | $\frac{1}{3}$ | $\frac{1}{2}$ |
| $\frac{3}{4}$ | $\frac{1}{2}$ | $\frac{1}{4}$ | $1$ | $1$ |
| $\frac{21}{25}$ | $\frac{2}{5}$ | $\frac{1}{10}$ | $3$ | $\frac{3}{2}$ |
| $\frac{551}{625}$ | $\frac{\sqrt{74}}{25}$ | $- \frac{1}{2} + \frac{3 \sqrt{74}}{50}$ | $\frac{403}{41} + \frac{50 \sqrt{74}}{41}$ | $-1 + \frac{25 \sqrt{74}}{74}$ |

最後は $y_R= \sqrt{74}/25$ で $3y_R$ が $1$ に近く、中心圧は密度を大きく上回る。これに対して表面の $z_R$ は有限なままである。$C\ge8/9$ を同じ式に入れると、中心から表面まで正の lapse と有限圧力を保つという前提を満たさなくなる。
:::

### 単調性の鎖でどこに等号が入るか

表面の質量を $m=GM(R)$ とする。非増加の密度なら $m(r)/r^3\ge m/R^3$ であるから、$y(r)\le\sqrt{1-2mr^2/R^3}$ である。また $y\zeta'/r\ge m/R^3$ である。二つを組み合わせると各 $r>0$ で
$$
\zeta'(r)\ge \frac{m}{R^3}\frac{r}{\sqrt{1-2mr^2/R^3}} .
$$
従って $\zeta(0)\le(3\sqrt{1-2m/R}-1)/2$ である。この不等式の右辺は同じ表面質量と半径の一様密度星の中心時刻率である。もし $\rho$ が非増加で一様でなければ、ある半径区間で平均密度が表面までの平均より厳密に大きい。そこで $y(r)<\sqrt{1-2mr^2/R^3}$ となるため、中心の不等式も厳密になる。曲率や圧力を考慮した中心時計の比較で、一様密度が特別な位置を占める理由である。

:::example[表面の赤方偏移から中心時計の率を絞る]
表面赤方偏移が $z_R=1$ なら $\sqrt{1-2m/R}=1/2$ である。単調密度の等方圧星では $\zeta(0)\le1/4$、従って中心から無限遠への赤方偏移は $z_c\ge3$ となる。一様密度なら量は両方で等号を取る。局所の密度分布を指定しなくても、仮定を満たす星の中心について表面の測定から下界を得られる。
:::

球対称でも動径圧 $p_r$ と横圧 $p_\perp$ が異なると、この証明の途中が変わる。保存則の動径成分では $p_r'=-(\rho+p_r)\Phi'+2(p_\perp-p_r)/r$ となる。これは角方向の接続 $\Gamma^\theta{}_{r\theta}=\Gamma^\phi{}_{r\phi}=1/r$ が応力の差に掛かるためである。場の角・動径成分の差を $8\pi G(p_\perp-p_r)$ として先ほどの式を整理すると
$$
y\dv{}{r}\left(\frac{y\zeta'}r\right)
=\zeta\left[\dv{}{r}\left(\frac{m(r)}{r^3}\right)+\frac{8\pi G}{r}(p_\perp-p_r)\right].
$$
従って横圧が動径圧より大きいと、密度が単調でも右辺が正になる場合がある。等方性は単調性の証明で使った条件である。例えば一定密度で有限の正の $p_\perp-p_r$ がある区間では、右辺の符号は正になる。この指摘は等方圧の限界がそのまま従わない理由を示すもので、異方的な物質でどこまで圧縮できるかを特定するには応力の追加条件が必要である。

## 白色矮星と Chandrasekhar 質量

[Fermi 気体](/physics/statistical-mechanics/fermi-gas/) の零温度極限で電子を記述し、電子数密度を $n_e=Y_e\rho/m_u$ とする。この節で $\rho$ はイオンが支配する質量密度、$m_u$ は原子質量単位、$Y_e$ は重粒子一個当たりの電子数である。Pauli 原理の状態密度から $n_e=p_F^3/(3\pi^2\hbar^3)$ である。超相対論的電子の速度を $c$ とし、方向平均した運動量流を積分すると $p=(c/3\pi^2\hbar^3)\int_0^{p_F}q^3\dd q=cp_F^4/(12\pi^2\hbar^3)$ となる。従って $p=K\rho^{4/3}$、$K=(3\pi^2)^{1/3}\hbar c(Y_e/m_u)^{4/3}/4$ である。

Newton 平衡で $p'=-G\rho M/r^2$ を使い、$\rho=\rho_c\theta^3$、$r=a\xi$、$a^2=K\rho_c^{-2/3}/(\pi G)$ と取ると $(\xi^2\theta')'=-\xi^2\theta^3$ となる。$\theta(0)=1,\theta'(0)=0$ で積分して最初の零 $\xi_1$ を表面とすれば、質量は

$$
M=4\pi(K/(\pi G))^{3/2}[-\xi_1^2\theta'(\xi_1)] .
$$

中心密度が式から消える理由は、超相対論的縮退圧が $\rho^{4/3}$ に比例するためである。数値積分では $\theta=1-\xi^2/6+\xi^4/40+\order{\xi^6}$ から始める。$\xi_1\simeq6.897$、$-\xi_1^2\theta'(\xi_1)\simeq2.018$ で、$Y_e=1/2$ なら $M\simeq1.46M_\odot$ である。これが理想模型の :term[Chandrasekhar 質量]{en="Chandrasekhar mass" yomi="ちゃんどらせかーるしつりょう"} である。有限電子質量では低密度で $p\propto\rho^{5/3}$ に変わる。電子捕獲、温度、回転は別の補正である。

### 低密度の白色矮星と数値積分

$p_F\ll m_e c$ では粒子の速さを $q/m_e$ として運動量流を積分する。同じ状態密度を使うと $p=p_F^5/(15\pi^2\hbar^3m_e)=K_{\rm nr}\rho^{5/3}$、$K_{\rm nr}=(3\pi^2)^{2/3}\hbar^2(Y_e/m_u)^{5/3}/(5m_e)$ である。一般に $p=K\rho^{1+1/n}$ では $\rho=\rho_c\theta^n$、$a^2=(n+1)K\rho_c^{1/n-1}/(4\pi G)$ と取ると、圧力の導関数の $\theta^n$ が密度と相殺し $(\xi^2\theta')'=-\xi^2\theta^n$ になる。$n=3/2$ では $a\propto\rho_c^{-1/6}$、$M\propto\rho_c a^3\propto\rho_c^{1/2}$ である。したがって同じ組成で中心密度を四倍にすると、近似の範囲では質量が二倍、半径が $2^{-1/3}$ 倍になる。$M$ が増えるほど半径が小さくなり、電子が超相対論的になりやすい。

中心を含む数値積分では $\theta=1-\xi^2/6+n\xi^4/120+\order{\xi^6}$ を使う。この係数は $\theta^n=1+n(\theta-1)+\order{\xi^4}$ を Lane–Emden 方程式の右辺に入れて比較すると得られる。$\xi$ の小さい正の値で $\theta,\theta'$ を初期化し、系 $(y_1',y_2')=(y_2,-2y_2/\xi-y_1^n)$ を積分する。$n=3$ の例では刻みを $0.002$ にして四次の Runge–Kutta 法を用い、符号が変わる区間内で $\theta$ の一次補間により最初の零点を得る。そこで $-\xi^2\theta'$ を評価し、先の質量式に入れる。$n$ が整数でないときは、未知関数が負になった後に $y_1^n$ を計算せず、表面までで積分を止める。

### 非相対論と超相対論の間の電子圧

零温度の電子の速度は $v(q)=cq/\sqrt{q^2+m_e^2c^2}$ である。方向平均で $qv/3$ を一状態当たりの圧力寄与として積分し、$x=p_F/(m_ec)$ とすると
$$
p=\frac{m_e^4c^5}{3\pi^2\hbar^3}\int_0^x\frac{u^4\dd u}{\sqrt{1+u^2}}
=\frac{m_e^4c^5}{24\pi^2\hbar^3}
[x(2x^2-3)\sqrt{1+x^2}+3\operatorname{arsinh}x]
$$
である。$u=\sinh w$ と置けば積分は $\sinh^4 w$ になり、$\sinh^4w=(\cosh4w-4\cosh2w+3)/8$ を積分して戻せる。$x\ll1$ では被積分関数は $u^4$ が先頭で $p\propto x^5$、$x\gg1$ では $u^3$ が先頭で $p\propto x^4$ である。密度が $x^3$ に比例するので、二つの指数 $5/3,4/3$ を同じ式から得られる。

:::example[相対論的になる密度尺度]
$x=1$ に対応する質量密度は $\rho_*=(m_u/Y_e)(m_ec/\hbar)^3/(3\pi^2)$ である。$Y_e=1/2$ とし、$\hbar/(m_ec)=3.86\times10^{-13}\,\mathrm m$、$m_u=1.66\times10^{-27}\,\mathrm{kg}$ とすれば $\rho_*\simeq1.95\times10^9\,\mathrm{kg\,m^{-3}}$ である。$\rho=8\rho_*$ なら $x=2$ なので、厳密な圧力式を使えば二つの漸近式の間にある圧力を評価できる。
:::

### 重力と縮退圧を同時に拡大縮小する

白色矮星の Newton 模型で、総質量を固定して密度を $\rho_\lambda(\vb x)=\lambda^3\rho(\lambda\vb x)$ と変える。半径を縮める向きは $\lambda>1$ である。圧力が $p=K\rho^\gamma$ なら静止質量を除く内部エネルギーは $U=K\int\rho^\gamma\dd^3x/(\gamma-1)$ であり、これは $\lambda^{3(\gamma-1)}U$ になる。Newton の重力エネルギーは
$$
W=-\frac G2\int\frac{\rho(\vb x)\rho(\vb y)}{|\vb x-\vb y|}
\,\dd^3x\,\dd^3y
$$
なので、分母と体積を変換して $\lambda W$ となる。従って $E(\lambda)=\lambda^{3(\gamma-1)}U+\lambda W$ である。平衡で $E'(1)=0$ は $3(\gamma-1)U+W=0$ を与え、$E''(1)=3(\gamma-1)(3\gamma-4)U$ である。

$\gamma=5/3$ でこの変形に対する二次のエネルギーが正になるのに対し、$\gamma=4/3$ では二次が零になる。その場合、両エネルギーが縮小率に同じように比例するので、平衡を保つ総質量が固定され、中心密度だけを変えた相似解ができる。$\gamma<4/3$ なら、この一つの許容される変形だけで二次エネルギーを負にできる。逆にこの変形で正になっても、他の全ての摂動が安定とは直ちにいえない。

:::example[冷たい電子圧の有効指数]
有限質量の厳密式を $p=C_e I(x)$、$I(x)=\int_0^x u^4/\sqrt{1+u^2}\dd u$ と書く。$\rho\propto x^3$ から、同じ組成での指数は
$$
\Gamma_{\rm eff}=\dv{\log p}{\log\rho}=\frac{x^5}{3I(x)\sqrt{1+x^2}} .
$$
例えば $x=1$ では $I(1)=(-\sqrt2+3\operatorname{arsinh}1)/8$ であり、$\Gamma_{\rm eff}\simeq1.534$ である。指数は非相対論的な $5/3$ と超相対論的な $4/3$ の中間にある。星の内部では $x$ が半径とともに変わるため、この一箇所の指数を全体のポリトロープ指数に置き換えるときは近似になる。
:::

### 低密度の Newton 極限で動径振動を計算する

圧力 $p(r)$ の平衡に対し、球殻を $r+\xi(t,r)$ にずらす。殻の質量座標 $M(r)$ を固定して $\Delta$ を殻に追随する変化とすると、体積の相対変化は $\xi'+2\xi/r$ である。粒子数保存と断熱性より
$$
\Delta\rho=-\rho(\xi'+2\xi/r),\qquad
\Delta p=-\Gamma_1 p(\xi'+2\xi/r),\qquad
\Gamma_1=(\partial\log p/\partial\log\rho)_s .
$$
質量座標の運動方程式 $\ddot r=-4\pi r^2\partial_M p-GM/r^2$ を一次化する。平衡の $\partial_Mp=-GM/(4\pi r^4)$ と $\partial_M\Delta p=(\Delta p)'/(4\pi r^2\rho)$ を使うと
$$
\rho\ddot\xi=-(\Delta p)'-\frac{4p'}r\,\xi
$$
である。ここで $\xi=r\,u(r)e^{i\omega t}$ と置き、$\Delta p=-\Gamma_1p(ru'+3u)e^{i\omega t}$ を代入する。$r^3$ を掛けて整理すると、中心で $u$ が有限で、表面で $\Delta p=0$ の動径振動方程式
$$
(\Gamma_1 p r^4 u')'+r^3[(3\Gamma_1-4)p]'u+\rho r^4\omega^2u=0
$$
を得る。$\Gamma_1=4/3$ が一定なら $u=\mathrm{const.}$ は $\omega^2=0$ の解である。静的な相似変形の自由度と振動方程式の零モードが同じ内容を表す。

:::example[一様密度の Newton 星の相似振動]
$\rho=\rho_0$、$p'=-(4\pi G/3)\rho_0^2r$ とする。摂動時の $\Gamma_1$ を定数として $u$ を一定に取ると、方程式は $\omega^2=(4\pi G\rho_0/3)(3\Gamma_1-4)$ を与える。$\Gamma_1=5/3$ なら周期は $2\pi/\sqrt{4\pi G\rho_0/3}$ である。一様密度という平衡の指定と、動いたときの断熱応答 $\Gamma_1$ は別の仮定である。この例では摂動のときも密度が常に一定という条件を課していない。
:::

## 中性子星と質量の上限

中性子星では重力も強く、相互作用のある核物質の全エネルギー密度 $\rho$ と圧力を :eqref[eq-gr-tov] に入れる。冷たい一成分流体なら第一法則 $\dd(\rho/n)=-p\,\dd(1/n)$ から $\dd\rho=(\rho+p)\dd n/n$ を得る。局所の音速の二乗は $c_s^2=\dd p/\dd\rho$ であり、安定性と因果性から $0\le c_s^2\le1$ を課すモデルを考える。例えば $p=w(\rho-\rho_s)$、$0<w\le1$ で無次元変数を $r=(G\rho_s)^{-1/2}\hat r$、$GM=(G\rho_s)^{-1/2}\hat m$ と置くと方程式から $\rho_s$ が消える。よって同じ $w$ の質量・半径は密度尺度の逆平方根に比例する。状態方程式を固定して $p_c$ を動かすと質量系列が得られ、安定枝は動径の固有振動数の二乗が非負の部分である。最初の極大質量と動径モードの安定性変化を関係付けるには、保存粒子数を固定した摂動固有値問題が必要である。以下で動径モードの固有値問題を導く。極大での安定性変化を平衡系列全体から判定する定理には、摂動の保存量と零モードの対応が必要であり、その大域的な扱いは :cite[mtw-gravitation] を参照する。Buchdahl の限界は質量と半径の比を制約するもので、状態方程式を指定せず質量の一意な数値を与えるものではない。

### 粒子数を圧力から復元する

冷たい線形状態方程式 $p=w(\rho-\rho_s)$ に $\dd\rho/(\rho+p)=\dd n/n$ を用いる。$p=0$ で $n=n_s$ と規格化し積分すると
$$
n=n_s\left(\frac{(1+w)\rho-w\rho_s}{\rho_s}\right)^{1/(1+w)},\qquad
\mu=\frac{\rho+p}{n}.
$$
$\dd\log\mu=\dd p/(\rho+p)$ はこの式の微分から求まる。TOV で $\dd p/(\rho+p)=-\dd\Phi$ なので、静的平衡中で $\mu e^\Phi$ は一定である。例えば $w=1/3,\rho=4\rho_s$ の点では $(\rho+p)/\rho_s=5$ だから $n/n_s=5^{3/4}$、$\mu/\mu_s=5^{1/4}$ である。従って中心に近いこの点で lapse は表面の $5^{-1/4}$ 倍でなければならない。

### 低密度と高密度で質量半径系列が変わる理由

冷たい一成分ポリトロープを $p=K n^\gamma$、$\rho=m_b n+p/(\gamma-1)$ と取ると、$\dd\rho/\dd n=m_b+\gamma K n^{\gamma-1}/(\gamma-1)$ であり、$(\rho+p)/n$ と一致する。音速は
$$
c_s^2=\frac{\gamma p}{\rho+p}.
$$
低密度では $c_s^2\simeq\gamma K n^{\gamma-1}/m_b$ で、圧力は小さい。一方、高密度では $c_s^2\to\gamma-1$ となる。従って $\gamma>2$ の単純な模型は十分高密度で因果性条件を破る。そのまま TOV に入れれば形式的な平衡系列は得られるが、物質模型として使える圧力範囲を別に検査する。例えば $\gamma=3$ で $p/(\rho+p)>1/3$ なら超光速の音速である。$p=K n^3,\rho=m_b n+p/2$ なので、境界は $p=2m_bn/3$ となる。

:::example[表面を既知の独立変数にして TOV を数値積分する]
圧力が零になる点を探索する代わりに、比エンタルピーの対数に相当する量
$$
h(p)=\int_0^p\frac{\dd p'}{\rho(p')+p'}
$$
を独立変数とする。ここでは $p=w(\rho-\rho_s)$、$w=1/3$、$\rho_s>0$ の模型を取る。$\rho=\rho_s+p/w$ を入れて積分すると $h=w\log[1+(1+w)p/(w\rho_s)]/(1+w)$ で、従って $p/\rho_s=w(e^{(1+w)h/w}-1)/(1+w)$ である。表面は $h=0$ であり、中心値 $h_c>0$ を一つずつ選ぶ。

長さの単位を $(G\rho_s)^{-1/2}$、質量に $G$ を掛けて同じ長さの単位で測り、ここでは無次元値にも $r,M$ を用いる。$\rho,p$ は $\rho_s$ で割った数とする。:eqref[eq-gr-tov] を $h'(r)=p'(r)/(\rho+p)$ で割ると
$$
\begin{aligned}
\dv r h&=-\frac{r(r-2M)}{M+4\pi r^3p},\\
\dv M h&=-\frac{4\pi r^3(r-2M)\rho}{M+4\pi r^3p}.
\end{aligned}
$$
$h$ を中心から表面へ減少させると半径・質量が増える。中心は $0/0$ なので $h_i=h_c-\epsilon$ から始める。$h_c-h=2\pi(\rho_c/3+p_c)r^2+\order{r^4}$、$M=(4\pi/3)\rho_c r^3+\order{r^5}$ を使って初期値を作る。

積分では中心近くで微分が発散しないように $X=r^2$、$Y=M/r$ を使う。上式を微分し直すと
$$
\dv Xh=-\frac{2X(1-2Y)}{Y+4\pi Xp},\qquad
\dv Yh=-\frac{(1-2Y)(4\pi X\rho-Y)}{Y+4\pi Xp}.
$$
$\epsilon=10^{-8}$ とし、上で得た小半径の初期値を $X,Y$ に変換し、$h_i$ から $0$ を等分する。状態の列を $Z=(X,Y)$、右辺を $F(h,Z)$、負の刻み $d$ について $k_1=F(h,Z)$、$k_2=F(h+d/2,Z+dk_1/2)$、$k_3=F(h+d/2,Z+dk_2/2)$、$k_4=F(h+d,Z+dk_3)$ を求める。次の状態は $Z+d(k_1+2k_2+2k_3+k_4)/6$ である。これを繰り返し、$r>2M$ が保たれているかも途中で確認する。表の $R,M_R$ は $8000$ 等分の結果である。最後の二列は $2000$ 等分の値との差の絶対値で、初期値の $\epsilon$ による誤差とは別である。

| $h_c$ | $p_c$ | $\rho_c$ | $R$ | $M_R$ | $2M_R/R$ | $\abs{\Delta R}$ | $\abs{\Delta M_R}$ |
|---|---|---|---|---|---|---|---|
| 0.05 | 0.0553507 | 1.16605 | 0.135445 | 0.0110632 | 0.16336 | 2.63123e-14 | 9.03028e-14 |
| 0.1 | 0.122956 | 1.36887 | 0.170173 | 0.0233614 | 0.27456 | 1.11638e-12 | 4.52416e-14 |
| 0.2 | 0.306385 | 1.91916 | 0.195245 | 0.0401234 | 0.411005 | 5.5928e-12 | 1.83058e-12 |
| 0.4 | 0.988258 | 3.96477 | 0.194593 | 0.0513721 | 0.527996 | 1.38659e-11 | 8.05148e-12 |
| 0.6 | 2.50579 | 8.51738 | 0.178289 | 0.049916 | 0.559945 | 1.13459e-11 | 1.15476e-11 |
| 0.8 | 5.88313 | 18.6494 | 0.161979 | 0.0449039 | 0.554442 | 7.25373e-12 | 8.38768e-12 |

軽い枝では $h_c$ を上げると質量が増えるが、さらに高い中心値まで比べると質量が減る組がある。圧力を増しても圧力自身が慣性と場の源に入るため、同じ $\rho_s$ で質量がどこまでも増加するわけではない。ただしこの表は静水圧平衡だけを解いた結果であり、各解の動径安定性は振動方程式を解いて調べる。次の表は $h_c=0.4$ の一つの解を中心から表面まで追った結果である。

| $h$ | $r$ | $M(r)$ | $p$ |
|---|---|---|---|
| 0.3 | 0.0856141 | 0.0083471 | 0.580029 |
| 0.2 | 0.12553 | 0.0211698 | 0.306385 |
| 0.1 | 0.160373 | 0.035614 | 0.122956 |
| 4.99004e-17 | 0.194593 | 0.0513721 | 5.55112e-17 |

表面で圧力が零になってもこの状態方程式では密度は $\rho_s$ に留まる。外部は真空であるから表面で密度の跳びがあり、$M'$ の内外も異なる。この模型の表面への接続では $M$ と圧力の値が要るため、密度の跳びだけでは薄い殻を追加しない。物理単位へ戻す際には、例えば $\rho_s$ を質量密度で指定した SI の場合、長さは $c/\sqrt{G\rho_s}$、質量は $c^3/(G^{3/2}\sqrt{\rho_s})$ を掛ける。
:::

### 保存粒子数のもとで相対論的な動径モードを調べる

静的背景を $\dd s^2=-e^{2\Phi}\dd t^2+e^{2\lambda}\dd r^2+r^2\dd\Omega^2$ とし、$e^{-2\lambda}=f$ とする。以下では冷たい一成分流体が背景と摂動で同じ状態方程式を満たすとする。球対称の摂動を面積半径の座標で取り、流体要素の変位を $\xi(t,r)$ とする。固定した座標での変化は $\delta$、要素に追随する変化は $\Delta=\delta+\xi\partial_r$ と区別する。要素の粒子数は一次で $nr^2e^\lambda\dd r$ に比例するため
$$
\frac{\Delta n}n=-\xi'-\left(\frac2r+\lambda'\right)\xi-\delta\lambda .
$$
時間・動径の Einstein 方程式を一次化すると $\delta\dot\lambda=-4\pi Gr e^{2\lambda}(\rho+p)\dot\xi$ となる。実際、動径流の混合成分が一次で $T^r{}_t=-(\rho+p)\dot\xi$、幾何側が $G^r{}_t=2e^{-2\lambda}\dot\lambda/r$ だからである。ここでは時間依存モードを取り、時間微分を取り除く。背景式から $\Phi'+\lambda'=4\pi Gr e^{2\lambda}(\rho+p)$ なので $\delta\lambda=-(\Phi'+\lambda')\xi$ となる。

$\zeta=r^2e^{-\Phi}\xi$ とおけば、粒子数と第一法則から
$$
\frac{\Delta n}n=-\frac{e^\Phi}{r^2}\zeta',\qquad
\Delta p=-\Gamma_1p\frac{e^\Phi}{r^2}\zeta',\qquad
\Delta\rho=-(\rho+p)\frac{e^\Phi}{r^2}\zeta',
\quad
\Gamma_1p=(\rho+p)c_s^2 .
$$
同じ測地線距離で測った殻の変形に比べて $\Phi'$ の項が入るのは、重力場の摂動 $\delta\lambda$ も同時に解いたためである。

動径の場の方程式と Euler 方程式を残す。$G^r{}_r=(f-1)/r^2+2f\Phi'/r$ の変分は
$$
\delta\Phi'=4\pi G r e^{2\lambda}\delta p+(\tfrac1r+2\Phi')\delta\lambda .
$$
一方、Euler 方程式は
$$
(\rho+p)e^{2\lambda-2\Phi}\ddot\xi
=-\delta p'-(\delta\rho+\delta p)\Phi'-(\rho+p)\delta\Phi'.
$$
後者は動径加速度の一次変化が $e^{2\lambda-2\Phi}\ddot\xi+\delta\Phi'$ であることから得られる。静的な圧力勾配は $p'=-(\rho+p)\Phi'$ で消えている。

$\xi=e^\Phi\zeta/r^2$、$\delta p=\Delta p-\xi p'$、$\delta\rho=\Delta\rho-\xi\rho'$ を右辺に入れる。$\zeta'$ に比例する項は $\zeta'$ を含む $\delta p'$ の微分と $(\lambda'+2\Phi')\Gamma_1p e^\Phi\zeta'/r^2$ の和になる。$\zeta$ に比例する項をまとめるには $p''=-(\rho'+p')\Phi'-(\rho+p)\Phi''$ を使い、さらに角成分の背景式で $\Phi''$ を消す。全体に $e^{\lambda+2\Phi}$ を掛け、$\ddot\zeta=-\omega^2\zeta$ とすれば
$$
(\Pi\zeta')'+(Q+\omega^2 W)\zeta=0,\quad
\begin{aligned}
\Pi&=\frac{e^{\lambda+3\Phi}}{r^2}\Gamma_1p,\\
Q&=\frac{e^{\lambda+3\Phi}}{r^2}\left[
\frac{(p')^2}{\rho+p}-\frac{4p'}r-8\pi G e^{2\lambda}p(\rho+p)\right],\\
W&=\frac{e^{3\lambda+\Phi}}{r^2}(\rho+p).
\end{aligned}
$$
中心では滑らかな動径ベクトルに対して $\xi=\order{r}$、従って $\zeta=\order{r^3}$ である。表面の流体要素は移動後も真空と接するから $\Delta p=0$ を課す。内部で $\rho+p>0$、$c_s^2>0$ なら $\Pi,W>0$ である。$\zeta$ を方程式に掛けて積分し、境界項を消すと
$$
\omega^2=\frac{\displaystyle\int_0^R(\Pi|\zeta'|^2-Q|\zeta|^2)\,\dd r}
{\displaystyle\int_0^R W|\zeta|^2\,\dd r}.
$$
分子の第一項は圧縮に伴う復元力を含む。外部で規格化した時刻 $t$ に関して $\omega^2<0$ なら指数的に成長する解がある。単一の試行関数で分子を負にできれば、同じ境界条件の自己共役な動径問題の最低固有値も負になる。

:::example[動径振動式の弱い場への極限]
$e^\Phi,e^\lambda\simeq1$、$p/\rho\ll1$、$GM/r\ll1$ とすると $\Pi\simeq\Gamma_1p/r^2$、$Q\simeq-4p'/r^3$、$W\simeq\rho/r^2$ である。$Q$ の残りの二項は、それぞれ圧力勾配の二乗と $G p\rho$ に比例し、Newton の静水圧の次数より高い。$\zeta\simeq r^3u$ を代入して $r^3$ を掛ければ、前節の Newton 方程式になる。従って $\Gamma_1=4/3$ の Newton 模型の中立モードには相対論的補正が効く。最大質量を調べるには平衡系列と固有値の両方で、同じ状態方程式と表面条件を用いる。
:::

## 重力崩壊と Oppenheimer–Snyder 解

圧力を零にした一様球の崩壊では、内部を閉じた同質等方計量 $-\dd\tau^2+a(\tau)^2(\dd\chi^2+\sin^2\chi\,\dd\Omega^2)$ とする。縮退圧のような平衡を支える力がないので、境界 $\chi=\chi_b$ も粒子に沿って動く。密度保存は $\rho a^3=\rho_0 a_0^3$、時間成分の Einstein 方程式は $\dot a^2+1=(8\pi G/3)\rho a^2$ である。この章で用いる閉じた塵の内部方程式は節末で直接導く。一般の流体と空間曲率を含む形は [相対論的宇宙論](/physics/general-relativity/cosmology/) で扱う。$\dot a=0$ の $a=a_0$ から始めると $\dot a^2=a_0/a-1$ である。$a=a_0(1+\cos\eta)/2$、$\dd\tau=a\dd\eta$ を代入して崩壊側を選ぶと $\tau=a_0(\eta+\sin\eta)/2$、$\tau_{\rm sing}=\pi a_0/2$ となる。

境界の面積半径は $R=a\sin\chi_b$ である。内部の $g^{ab}\partial_a r\partial_b r=\cos^2\chi-\dot a^2\sin^2\chi$ を外部の $1-2m/R$ と比較すると $m=a_0\sin^3\chi_b/2$、すなわち $M=4\pi\rho R^3/3$ が一定となる。$0<\chi_b<\pi/2$ を取る。外部の境界軌道は $\dot R^2=2m/R-\sin^2\chi_b$ であり、$E=\cos\chi_b$ の動径測地線である。内外で誘導計量は $-\dd\tau^2+R^2\dd\Omega^2$、法線方向の角外的曲率は $n(r)/R=\cos\chi_b/R$、時間成分は測地運動のため零である。従って薄い殻を要しない。これが :term[Oppenheimer–Snyder 崩壊]{en="Oppenheimer–Snyder collapse" yomi="おっぺんはいまーすないだーほうかい"} である。表面が $R=2m$ に達する条件は $a/a_0=\sin^2\chi_b$ で、$a=0$ より前である。

### 崩壊内部で地平線の生成子を追う

表面が外部地平線と交わるのは $(1+\cos\eta_H)/2=\sin^2\chi_b$ より $\eta_H=\pi-2\chi_b$ である。内部の動径計量は $a^2(-\dd\eta^2+\dd\chi^2)$ なので外向き光は $\chi-\eta$ が一定である。表面の交点 $(\eta_H,\chi_b)$ を通る光を過去に追うと、中心での時刻は $\eta_H-\chi_b=\pi-3\chi_b$ である。$\chi_b<\pi/3$ なら崩壊の初期面より後に中心を出て、表面と $R=2m$ で出会う。これより早く中心を出た光は表面を $R>2m$ で横切れる。遅ければ表面を横切るときにはすでに外部の未来地平線の内側である。従って事象の地平線の決定には表面の時空全体での運動が必要である。

:::example[外からの半径と内部の捕捉面]
$\chi_b=\pi/6$ なら表面との交点は $\eta_H=2\pi/3$、中心から出るのは $\eta=\pi/2$ である。一方、内部の球の両光的膨張率の符号は $k_\pm(r)=\dot a\sin\chi\pm\cos\chi$ で決まる。崩壊では $\dot a<0$ であり、外向きの膨張が零になる条件は $\dot a=-\cot\chi$ である。$\dot a^2=a_0/a-1$ から $a/a_0=\sin^2\chi$ を得る。この局所条件を満たす面と、中心を先に通る事象地平線の生成子を区別できる。
:::

### 密度と残りの崩壊時間

Oppenheimer–Snyder の崩壊で $s=\pi-\eta$ を小さく取る。$a=a_0(1-\cos s)/2$ は $a_0s^2/4$ が先頭で、$\tau_{\rm sing}-\tau=a_0(s-\sin s)/2$ は $a_0s^3/12$ が先頭である。従って $a\propto(\tau_{\rm sing}-\tau)^{2/3}$ であり、$\rho\propto a^{-3}$ から密度は残り固有時の逆二乗で発散する。係数まで求めるには $(8\pi G/3)\rho a^3=a_0$ を使う。$a^3\simeq(9a_0/4)(\tau_{\rm sing}-\tau)^2$ なので $\rho\simeq1/[6\pi G(\tau_{\rm sing}-\tau)^2]$ となる。境界の角度はこの局所の係数には入らないが、境界が地平線に入る時刻には入る。

### 崩壊の境界に沿った時刻と赤方偏移

境界の固有時 $\tau$ に沿って、外部の入射時刻 $v=t+r_*$ と遅延時刻 $u=t-r_*$ を評価する。外部の動径運動の定数 $E=\cos\chi_b$、$\dot R=-\sqrt{E^2-f(R)}$ を使うと
$$
\dot v=\frac{E+\dot R}{f}=\frac1{E-\dot R},\qquad
\dot u=\frac{E-\dot R}{f}=\frac1{E+\dot R}.
$$
表面が地平線を通る時刻を $\tau_H$ とすれば、$\dot R(\tau_H)=-E$ である。従って $\dot v$ は $1/(2E)$ に近づくが、$\dot u$ は発散する。外向き光では $u$ が保たれるため、無限遠で規格化した座標時刻で受け取るパルス間隔は $\dot u$ 倍になる。実際、$R-2m\simeq E(\tau_H-\tau)$、$f\simeq(R-2m)/(2m)$ から $\dot u\simeq4m/(\tau_H-\tau)$ である。有限半径の静止観測者なら、さらに受信側の固有時率 $\sqrt{f(r_{\rm o})}$ を掛ける。

:::example[表面から地平線直前に出た二つの信号]
$\tau_H-\tau_1=2\epsilon$、$\tau_H-\tau_2=\epsilon$ とする。外部の遅延時刻の差は先頭で $4m\int_\epsilon^{2\epsilon}\dd s/s=4m\log2$ であり、放射側の固有時差 $\epsilon$ を小さくしても有限の値に近づく。これに対して同じ半径に到着した二つの波頂の瞬間の周波数の比には、固有時率の逆数が入る。放射側の固有振動数を一定に保っても、遠方での受信周波数は境界が地平線へ近づくにつれて零に近づく。
:::

### 初期半径が崩壊時間を決める

初期の半径は $R_0=a_0\sin\chi_b$ なので $m=a_0\sin^3\chi_b/2$ を使うと $\sin^2\chi_b=2m/R_0$ である。従って初期の静止時刻から中心の終端までの固有時は
$$
\tau_{\rm sing}=\frac{\pi}{2}\sqrt{\frac{R_0^3}{2m}} .
$$
この形は動径自由落下の Newton の時間尺度に一致する。内部の一様な塵では殻それぞれの初期質量と半径の比 $R_0^3/m$ が同じなので、全ての殻で $a=0$ が同じ共動時刻に来る。一様性を外すと殻ごとの崩壊時間も変わるため、ここから一般の崩壊にも同時性を仮定することはできない。

:::example[初期のコンパクトネスが四分の一の球]
$R_0=8m$ なら $\sin\chi_b=1/2$、$\chi_b=\pi/6$、$a_0=16m$ である。従って中心の特異境界までの固有時は $8\pi m$ である。中心の地平線の生成子は $\eta=\pi/2$ から出るので、その時刻は $8m(\pi/2+1)$ である。表面が地平線に入る時刻は $8m(2\pi/3+\sqrt3/2)$ であり、中心と表面での地平線の時刻は同じでない。ここで比べた時刻は、共動する各時計を一様球の初期面でそろえた値である。
:::

### 閉じた塵の内部方程式を直接求める

内部の計量を $\dd s^2=-\dd\tau^2+a^2(\dd\chi^2+\sin^2\chi\,\dd\Omega^2)$ とする。空間単位球面の計量を $\gamma_{ij}$ とすると、時間を含む非零の接続は $\Gamma^0{}_{ij}=a\dot a\gamma_{ij}$ と $\Gamma^i{}_{0j}=(\dot a/a)\delta^i_j$ である。球面の断面曲率は $+1$ であり、これに時間接続の積が加わって空間曲率項は $1+\dot a^2$ になる。具体的に $R_{00}=-3\ddot a/a$、$R_{ij}=[a\ddot a+2(1+\dot a^2)]\gamma_{ij}$ である。従って $G_{00}=3(1+\dot a^2)/a^2$、$G^i{}_j=-[2a\ddot a+1+\dot a^2]\delta^i_j/a^2$ である。塵の $T^i{}_j=0$ から $\ddot a=-(1+\dot a^2)/(2a)$ となり、静止から始めても $\ddot a<0$ なので崩壊側へ動く。時間成分は $\dot a^2+1=(8\pi G/3)\rho a^2$ を与える。圧力が零なので $\rho a^3$ は保存され、パラメータ解は空間成分も満たす。

:::example[境界から見える内側の質量を球面ごとに確かめる]
半径 $r=a\sin\chi$ の任意の共動する球で、$1-(\nabla r)^2=(1+\dot a^2)\sin^2\chi$ である。この量に $r/2$ を掛ければ幾何学的質量長は $G(4\pi/3)\rho r^3$ となり、固定した $\chi$ で一定である。例えば二殻の初期面積半径比が $1:2$ なら、内包質量の比は一貫して $1:8$ である。各殻で $2G M/r$ は $a^{-1}$ に比例して増えるため、外側の殻が先に $2GM/r=1$ を満たす。事象地平線の外向き生成子はこの局所条件で決めた軌跡とは異なる。
:::

## 演習問題

::::exercise{level=basic}
一様密度 $\rho_0$ の星の $M(r)$ を求めよ。

:::solution
$M(0)=0$ と $M'=4\pi r^2\rho_0$ から $M=4\pi\rho_0 r^3/3$ である。
:::
::::

::::exercise{level=basic}
一様密度星の $2GM/R=3/4$ の中心圧を密度の倍数で答えよ。

:::solution
:eqref[eq-gr-uniform-star] に $y(0)=1,y_R=1/2$ を代入すると $p_c/\rho_0=(1-1/2)/(3/2-1)=1$ である。
:::
::::

::::exercise{level=basic}
$p=w(\rho-\rho_s)$、$w=1/4$ の音速を光速の倍数で求めよ。

:::solution
$c_s^2=w=1/4$ なので $c_s=1/2$ である。
:::
::::

::::exercise{level=basic}
超相対論的縮退電子の $Y_e$ が半分になると質量尺度は何倍か。

:::solution
$K\propto Y_e^{4/3}$、$M\propto K^{3/2}$ なので $M\propto Y_e^2$ である。従って $1/4$ 倍になる。
:::
::::

::::exercise{level=advanced}
一様密度星の小さな $GM/R$ での中心圧を一次まで求めよ。

:::solution
$y_R=1-GM/R+\order{(GM/R)^2}$ なので $p_c/\rho_0=(GM/R)/(2-3GM/R)+\order{(GM/R)^2}=GM/(2R)+\order{(GM/R)^2}$ である。これは $p'=-4\pi G\rho_0^2r/3$ の積分と一致する。
:::
::::

::::exercise{level=advanced}
密度が非増加なら $M(r)/r^3$ も非増加であることを示せ。

:::solution
$3M(r)/(4\pi r^3)=3r^{-3}\int_0^r s^2\rho(s)\dd s\ge \rho(r)$ である。よって $(M/r^3)'=4\pi\rho(r)/r-3M/r^4\le0$ となる。
:::
::::

::::exercise{level=advanced}
塵の崩壊で $\chi_b=\pi/6$ の表面が $R=2m$ となる $\eta$ を求めよ。

:::solution
$(1+\cos\eta)/2=\sin^2\chi_b=1/4$ より $\cos\eta=-1/2$ である。崩壊区間 $0<\eta<\pi$ では $\eta=2\pi/3$ となる。
:::
::::

::::exercise{level=advanced}
TOV の中心圧が $p_c=\rho_c$ であるとき、$p$ の二次係数を求めよ。

:::solution
$M=(4\pi/3)\rho_c r^3+\order{r^5}$ により $p'=-4\pi G(2\rho_c)(4\rho_c/3)r+\order{r^3}$ である。積分すると $p=p_c-(16\pi G/3)\rho_c^2r^2+\order{r^4}$ となる。
:::
::::

::::exercise{level=advanced}
コンパクトネス $C$ の一様密度星で $p(r)=\rho_0/3$ となる点の面積半径を求めよ。$C=3/4$ の例でも計算し、どの範囲の $C$ でこの点が存在するかを判定せよ。

:::solution
:eqref[eq-gr-uniform-star] の式で、$p/\rho_0=1/3$ は $3(y-y_R)=3y_R-y$、従って $y=3y_R/2$ に等しい。$y=\sqrt{1-Cr^2/R^2}$、$y_R=\sqrt{1-C}$ を代入すると
$$
\frac{r^2}{R^2}=\frac{9C-5}{4C}
$$
である。$C=3/4$ では $r/R=\sqrt{7/12}$ となる。$C<5/9$ では中心圧も $\rho_0/3$ に届かず、$C=5/9$ で初めて中心に現れる。$5/9<C<8/9$ では $0<r<R$ である。$C\ge8/9$ では元の正則な一様密度星の仮定を満たさない。
:::
::::

## 確認テスト

::::quiz
密度が外向きに非増加の星について、TOV の Newton 極限に必要な小さい無次元量を選べ。

- [x] $p/\rho$
- [x] $2GM/r$
- [ ] 質量密度そのもの

:::explanation
$p/\rho$ が小さいと圧力の慣性と源の寄与を低次で落とせ、$2GM/r$ が小さいと分母が $1$ に近づく。密度には単位があり、密度の数値単独では近似の精度を定められない。
:::
::::

::::quiz{answer="0.888888889" tol="1e-7"}
Buchdahl の限界での $2GM/R$ の上限を答えよ。

:::explanation
中心で $\zeta(0)>0$ なら $3\sqrt{1-2GM/R}>1$ である。従って $2GM/R<8/9$ で、上限は約 $0.888889$ である。
:::
::::

::::quiz
TOV に必要な中心の正則条件を選べ。

- [x] $M(0)=0$
- [ ] $p(0)=0$ が常に必要である。
- [ ] 中心の密度が必ず無限大である。

:::explanation
有限密度なら球の質量は $r^3$ に比例するので $M(0)=0$ である。非零の有限中心圧から積分し、表面で $p=0$ とする。
:::
::::

::::quiz{answer="0.25" tol="1e-9"}
$Y_e=1/2$ を $1/4$ に変えたときの理想 Chandrasekhar 質量の比を求めよ。

:::explanation
$M\propto K^{3/2}\propto Y_e^2$ なので、比は $(1/2)^2=1/4$ である。
:::
::::

::::quiz{answer="0.5" tol="1e-7"}
一様密度で $2GM/R=5/9$ の星の表面から無限遠への赤方偏移 $z_R$ を求めよ。

:::explanation
表面 lapse は $\sqrt{1-5/9}=2/3$ である。従って $1+z_R$ はその逆数 $3/2$ で、$z_R=1/2$ である。
:::
::::

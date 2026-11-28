---
title: "作用原理と Hamilton 形式"
description: "境界を含めた計量変分を行い、空間切片から ADM の正準運動量と拘束を導く。"
notation:
  - symbol: "h_{ij}"
    meaning: "空間切片上の誘導計量"
  - symbol: "K_{ij}"
    meaning: "未来法線に関する外的曲率"
  - symbol: "N,\\beta^i"
    meaning: "lapse と shift"
  - symbol: "\\pi^{ij}"
    meaning: "誘導計量の正準運動量"
---

## この章の目標

- 誘導計量と外的曲率を定義し、曲率の射影を計算する。
- 境界項を含む Einstein–Hilbert 作用を変分する。
- Palatini 変分が Levi-Civita 接続を与える条件を確認する。
- ADM 作用から拘束と正準運動量を求める。

## 超曲面と誘導計量・外的曲率

空間的な超曲面の未来向き単位法線を $n^a$ とする。誘導計量は $h_{ab}=g_{ab}+n_an_b$ を切片に制限したものであり、$h^a{}_b$ は接空間への射影である。:term[外的曲率]{en="extrinsic curvature" yomi="がいてききょくりつ"} を $K_{ab}=h_a{}^ch_b{}^d\nabla_c n_d$ とする。接ベクトル $X,Y$ では $K(X,Y)-K(Y,X)=-n\cdot[X,Y]=0$ である。$\mathcal L_n g$ を二方向で射影すると $2K$ となる。また $\nabla_XY$ の法線との内積は $-K(X,Y)$ だから、$\nabla_XY=D_XY+K(X,Y)n$ である。接成分 $D$ はねじれ零で $h$ を保つので :ref[thm-gr-levi-civita] により切片の接続である。例えば Minkowski 時空の $t$ 一定面で $K=0$、Milne の双曲面では $K_{ij}=h_{ij}/t$ となる。

### 法線を表すレベル面の関数

超曲面を $F=\mathrm{const}$、$\nabla F$ を過去向き時間的とする。未来法線を $n_a=-N_F\partial_aF$、$N_F=(-g^{ab}\partial_aF\partial_bF)^{-1/2}$ と取る。射影した微分では $\partial_aN_F$ の項が $\partial_bF$ を含んで消え、
$$
K_{ab}=-N_F h_a{}^c h_b{}^d\nabla_c\nabla_dF
$$
である。Hessian が対称なので $K$ も対称である。同じ面を $F\mapsto f(F)$、$f'>0$ で表し直すと、射影後の Hessian は $f'$ 倍、$N_F$ は $1/f'$ 倍になる。従って外的曲率は面のラベルの尺度に依存しない。

:::example[Milne の切片]
Minkowski の未来円錐で $\tau=(T^2-X^2-Y^2-Z^2)^{1/2}$ を使う。$\partial_a\tau=-x_a/\tau$、従って $n_a=x_a/\tau$ である。$\partial_b n_a=\eta_{ab}/\tau+x_ax_b/\tau^3$ となるから $K_{ab}=h_{ab}/\tau$ である。時空の曲率は零だが、埋め込まれた切片の外的曲率は零ではない。
:::

### 時間的境界と光的境界

一般の非光的面の単位法線を $s$、$\epsilon=g(s,s)$ とする。$h_{ab}=g_{ab}-\epsilon s_as_b$、$h^a{}_b=\delta^a_b-\epsilon s^as_b$ とすれば、$h^a{}_bs^b=0$ かつ $h^a{}_ch^c{}_b=h^a{}_b$ である。時間的な面では $\epsilon=+1$ で、誘導計量は Lorentz 計量になる。以下、面に接する二つの場 $X,Y$ に対し $K_s(X,Y)=g(\nabla_X s,Y)$ とする。直交条件を微分すると $g(s,\nabla_XY)=-K_s(X,Y)$ であるから
$$
\nabla_XY=D_XY-\epsilon K_s(X,Y)s.
$$
従って法線の向きを反転すると $K_s$ の符号が反転する一方、法線成分 $-\epsilon K_s(X,Y)s$ は変わらない。曲面上だけで $s$ を定義しても、$X$ が接方向なので $\nabla_Xs$ を評価できる。

光的な面では例えば Minkowski 時空の $u=t-x=0$ を取る。法線 $\ell^a=-g^{ab}\partial_bu$ は面に接し、誘導二次形式の核を張る。このため $h^a{}_b$ に相当する一意な直交射影を法線一本からは作れない。もう一本の光的方向 $k$ を $g(\ell,k)=-1$ と取れば $q^a{}_b=\delta^a_b+\ell^ak_b+k^a\ell_b$ が両者と直交する二方向への射影となる。実際、$q(\ell)=q(k)=0$ であり、両者に直交する $X$ では $q(X)=X$ である。$q_a{}^c q_b{}^d\nabla_c\ell_d$ の跡は光束の断面積の対数微分である。ここでは非光的境界の変分を導き、光束の計算は [因果構造と特異点定理](/physics/general-relativity/causal-structure/) で扱う。

:::example[傾けた平面と曲げた時刻面]
Minkowski 時空で $t=vx$、$|v|<1$ とすると法線一形式は $(-\dd t+v\dd x)/\sqrt{1-v^2}$ で一定である。切片の計量は $(1-v^2)\dd x^2+\dd y^2+\dd z^2$、外的曲率は零である。$v\to1$ では誘導計量の $xx$ 成分が零となり、単位法線はこの極限では使えない。一方、$t=F(x)$、$|F'|<1$ なら、接方向 $e_x=\partial_x+F'\partial_t$ で $\nabla_{e_x}e_x=F''\partial_t$ である。未来法線との縮約から $K_{xx}=F''/\sqrt{1-(F')^2}$、他の成分は零である。固有長で測った外的曲率の唯一の非零候補は $K^x{}_x=F''/(1-(F')^2)^{3/2}$ となる。誘導計量の係数は $X(x)=\int_0^x\sqrt{1-(F'(s))^2}\dd s$ で定数に直せるので、切片の内的曲率は零のままである。
:::

### 適合座標から断面の変化を読む

非光的面から出る単位法線測地線の固有距離または固有時 $y$ を使うと、法線指数写像が局所微分同相である近傍では計量を $\epsilon\dd y^2+h_{ij}(y,x)\dd x^i\dd x^j$ と書ける。$\partial_y$ のノルムは測地線に沿って一定である。また $[\partial_y,\partial_i]=0$ なので $\partial_y g(\partial_y,\partial_i)=g(\partial_y,\nabla_i\partial_y)=\partial_i g(\partial_y,\partial_y)/2=0$ であり、初期面の直交性が保存される。この座標を :term[Gaussian 法線座標]{en="Gaussian normal coordinates" yomi="がうすほうせんざひょう"} という。

計量整合性と交換関係から $K_{ij}=\partial_yh_{ij}/2$、また $\partial_y\log\sqrt{|\det h|}=h^{ij}\partial_yh_{ij}/2=K$ となる。従って $K$ は法線方向へ移った面の体積密度の対数変化率である。ただし、$K=0$ でも各方向の長さが不変とは限らない。例えば局所計量 $-\dd y^2+e^{2by}\dd x^2+e^{-2by}\dd z^2+\dd w^2$ では $K^i{}_j$ の固有値は $(b,-b,0)$ である。面の体積は一定だが、一方向の辺の長さは増え、もう一方向では減る。

## Gauss–Codazzi 方程式

:term[Gauss–Codazzi 方程式]{en="Gauss–Codazzi equations" yomi="がうすこだっつぃほうていしき"} は時空曲率の射影を記述する。$\nabla_YZ=D_YZ+K(Y,Z)n$ を $\nabla_X$ で微分すると、接成分は $D_XD_YZ+K(Y,Z)\nabla_Xn$、法線係数は $K(X,D_YZ)+X(K(Y,Z))$ である。$X,Y$ を交換して引き、括弧の項を引くと

$$
R_{ijkl}^{\rm proj}=R_{ijkl}^{(3)}+K_{ik}K_{jl}-K_{il}K_{jk},
\qquad n_aR^a{}_{jkl}=-(D_kK_{lj}-D_lK_{kj}).
$$

第一式を二度縮約した空間成分を $S$ とすると、逆計量 $g^{ab}=h^{ab}-n^an^b$ から $R=S-2R_{nn}$ なので $G_{nn}=S/2$ である。第二式の対を交換し、空間添字を縮約すると

$$
2G_{nn}=R^{(3)}+K^2-K_{ij}K^{ij},\qquad
G_{ni}=D_jK^j{}_i-D_iK .
\label{eq-gr-gauss-constraints}
$$

例えば Milne 面では $R^{(3)}=-6/t^2$ と $K^2-K_{ij}K^{ij}=6/t^2$ が相殺し、真空の拘束を満たす。

### 一様だが異方的な切片

$g=-\dd t^2+\sum_i a_i(t)^2(\dd x^i)^2$ とする。$H_i=\dot a_i/a_i$、$K^i{}_j=H_i\delta^i_j$ で、各切片の内的曲率は零である。Gauss の式から $G_{nn}=\tfrac12[(\sum H_i)^2-\sum H_i^2]=\sum_{i<j}H_iH_j$ となる。運動量拘束の空間微分も全て零である。真空で $a_i=t^{p_i}$ とすると $R^i{}_i=p_i(\sum_jp_j-1)/t^2$、$R_{tt}=(\sum_i p_i-\sum_i p_i^2)/t^2$ である。従って $\sum p_i=\sum p_i^2=1$ の Kasner 指数は真空発展を与える。この式は計量の時間接続 $\Gamma^t{}_{ii}=a_i\dot a_i$、$\Gamma^i{}_{ti}=H_i$ を曲率に入れ、$\dot H_i$ と $H_i\sum H_j$ を足して得られる。例えば $(-1/3,2/3,2/3)$ はこの二条件を満たし、一軸が縮み二軸が伸びる。

:::example[異方的な計量で Gauss の二重縮約を検算する]
$a(t),b(t),c(t)>0$ の Bianchi I 計量とする。切片の座標は直交し、各切片の内的曲率は零であるが、外的曲率の固有値はそれぞれ $a'/a,b'/b,c'/c$ である。時空から直接 $G_{tt}$ を作って照合する。 座標の順序を $(t,x,y,z)$ と固定する。成分行列と逆行列は次の組になる。

$$
\begin{aligned}
(g_{ab})&=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & a^{2}{\left(t \right)} & 0 & 0\\0 & 0 & b^{2}{\left(t \right)} & 0\\0 & 0 & 0 & c^{2}{\left(t \right)}\end{matrix}\right]\\
(g^{ab})&=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & \frac{1}{a^{2}{\left(t \right)}} & 0 & 0\\0 & 0 & \frac{1}{b^{2}{\left(t \right)}} & 0\\0 & 0 & 0 & \frac{1}{c^{2}{\left(t \right)}}\end{matrix}\right]\\
\det g&=- a^{2}{\left(t \right)} b^{2}{\left(t \right)} c^{2}{\left(t \right)}
\end{aligned}
$$

逆行列を掛けて単位行列を得ることを確かめてから、:eqref[eq-gr-christoffel] を使う。同じ下添字の順序を二度計算せず、$b\le c$ の組を選ぶ。中央の列には、座標に依存する計量成分の微分だけを残した。最後の列は微分を実行した値であり、表に出ない成分は添字交換で得られるものを除いて零である。

| 成分 | 微分への代入 | 結果 |
|---|---|---|
| $\Gamma^{t}_{x x}$ | $- \frac{1}{2}\left(-\partial_{t}\left(a^{2}{\left(t \right)}\right)\right)$ | $a{\left(t \right)} a'{\left(t \right)}$ |
| $\Gamma^{t}_{y y}$ | $- \frac{1}{2}\left(-\partial_{t}\left(b^{2}{\left(t \right)}\right)\right)$ | $b{\left(t \right)} b'{\left(t \right)}$ |
| $\Gamma^{t}_{z z}$ | $- \frac{1}{2}\left(-\partial_{t}\left(c^{2}{\left(t \right)}\right)\right)$ | $c{\left(t \right)} c'{\left(t \right)}$ |
| $\Gamma^{x}_{t x}$ | $\frac{1}{2 a^{2}{\left(t \right)}}\left(\partial_{t}\left(a^{2}{\left(t \right)}\right)\right)$ | $\frac{a'{\left(t \right)}}{a{\left(t \right)}}$ |
| $\Gamma^{y}_{t y}$ | $\frac{1}{2 b^{2}{\left(t \right)}}\left(\partial_{t}\left(b^{2}{\left(t \right)}\right)\right)$ | $\frac{b'{\left(t \right)}}{b{\left(t \right)}}$ |
| $\Gamma^{z}_{t z}$ | $\frac{1}{2 c^{2}{\left(t \right)}}\left(\partial_{t}\left(c^{2}{\left(t \right)}\right)\right)$ | $\frac{c'{\left(t \right)}}{c{\left(t \right)}}$ |

曲率の計算では :eqref[eq-gr-curvature-components] の微分差と積差を分ける。第二列は $\partial_c\Gamma^a{}_{db}$、第三列は $-\partial_d\Gamma^a{}_{cb}$ であり、第四列は接続の積を符号付きで足したものである。三列を足して最後の列を得る。

| 成分 | 第一の微分 | 第二の微分 | 接続の積 | 結果 |
|---|---|---|---|---|
| $R^{t}_{x t x}$ | $a{\left(t \right)} a''{\left(t \right)} + \left(a'{\left(t \right)}\right)^{2}$ | $0$ | $- \left(a'{\left(t \right)}\right)^{2}$ | $a{\left(t \right)} a''{\left(t \right)}$ |
| $R^{t}_{y t y}$ | $b{\left(t \right)} b''{\left(t \right)} + \left(b'{\left(t \right)}\right)^{2}$ | $0$ | $- \left(b'{\left(t \right)}\right)^{2}$ | $b{\left(t \right)} b''{\left(t \right)}$ |
| $R^{t}_{z t z}$ | $c{\left(t \right)} c''{\left(t \right)} + \left(c'{\left(t \right)}\right)^{2}$ | $0$ | $- \left(c'{\left(t \right)}\right)^{2}$ | $c{\left(t \right)} c''{\left(t \right)}$ |
| $R^{x}_{y x y}$ | $0$ | $0$ | $\frac{b{\left(t \right)} a'{\left(t \right)} b'{\left(t \right)}}{a{\left(t \right)}}$ | $\frac{b{\left(t \right)} a'{\left(t \right)} b'{\left(t \right)}}{a{\left(t \right)}}$ |
| $R^{x}_{z x z}$ | $0$ | $0$ | $\frac{c{\left(t \right)} a'{\left(t \right)} c'{\left(t \right)}}{a{\left(t \right)}}$ | $\frac{c{\left(t \right)} a'{\left(t \right)} c'{\left(t \right)}}{a{\left(t \right)}}$ |
| $R^{y}_{z y z}$ | $0$ | $0$ | $\frac{c{\left(t \right)} b'{\left(t \right)} c'{\left(t \right)}}{b{\left(t \right)}}$ | $\frac{c{\left(t \right)} b'{\left(t \right)} c'{\left(t \right)}}{b{\left(t \right)}}$ |

第一添字を計量で下ろすと $R_{abcd}=g_{ak}R^k{}_{bcd}$ となる。対の反対称性と交換対称性を用いて、非零の独立な組を下表に示す。

| 成分 | 結果 |
|---|---|
| $R_{t x t x}$ | $- a{\left(t \right)} a''{\left(t \right)}$ |
| $R_{t y t y}$ | $- b{\left(t \right)} b''{\left(t \right)}$ |
| $R_{t z t z}$ | $- c{\left(t \right)} c''{\left(t \right)}$ |
| $R_{x y x y}$ | $a{\left(t \right)} b{\left(t \right)} a'{\left(t \right)} b'{\left(t \right)}$ |
| $R_{x z x z}$ | $a{\left(t \right)} c{\left(t \right)} a'{\left(t \right)} c'{\left(t \right)}$ |
| $R_{y z y z}$ | $b{\left(t \right)} c{\left(t \right)} b'{\left(t \right)} c'{\left(t \right)}$ |

Ricci テンソルは第一・第三添字を縮約して得る。曲率が残る場合も相殺が起こるので、各項を落とさずに加える。

| 成分 | 縮約の項 | 代入 | 結果 |
|---|---|---|---|
| $R_{t t}$ | $R^{x}_{t x t}+R^{y}_{t y t}+R^{z}_{t z t}$ | $- \frac{c''{\left(t \right)}}{c{\left(t \right)}} - \frac{b''{\left(t \right)}}{b{\left(t \right)}} - \frac{a''{\left(t \right)}}{a{\left(t \right)}}$ | $- \frac{a{\left(t \right)} b{\left(t \right)} c''{\left(t \right)} + a{\left(t \right)} c{\left(t \right)} b''{\left(t \right)} + b{\left(t \right)} c{\left(t \right)} a''{\left(t \right)}}{a{\left(t \right)} b{\left(t \right)} c{\left(t \right)}}$ |
| $R_{x x}$ | $R^{t}_{x t x}+R^{y}_{x y x}+R^{z}_{x z x}$ | $a{\left(t \right)} a''{\left(t \right)} + \frac{a{\left(t \right)} a'{\left(t \right)} c'{\left(t \right)}}{c{\left(t \right)}} + \frac{a{\left(t \right)} a'{\left(t \right)} b'{\left(t \right)}}{b{\left(t \right)}}$ | $\frac{\left(b{\left(t \right)} c{\left(t \right)} a''{\left(t \right)} + b{\left(t \right)} a'{\left(t \right)} c'{\left(t \right)} + c{\left(t \right)} a'{\left(t \right)} b'{\left(t \right)}\right) a{\left(t \right)}}{b{\left(t \right)} c{\left(t \right)}}$ |
| $R_{y y}$ | $R^{t}_{y t y}+R^{x}_{y x y}+R^{z}_{y z y}$ | $b{\left(t \right)} b''{\left(t \right)} + \frac{b{\left(t \right)} b'{\left(t \right)} c'{\left(t \right)}}{c{\left(t \right)}} + \frac{b{\left(t \right)} a'{\left(t \right)} b'{\left(t \right)}}{a{\left(t \right)}}$ | $\frac{\left(a{\left(t \right)} c{\left(t \right)} b''{\left(t \right)} + a{\left(t \right)} b'{\left(t \right)} c'{\left(t \right)} + c{\left(t \right)} a'{\left(t \right)} b'{\left(t \right)}\right) b{\left(t \right)}}{a{\left(t \right)} c{\left(t \right)}}$ |
| $R_{z z}$ | $R^{t}_{z t z}+R^{x}_{z x z}+R^{y}_{z y z}$ | $c{\left(t \right)} c''{\left(t \right)} + \frac{c{\left(t \right)} b'{\left(t \right)} c'{\left(t \right)}}{b{\left(t \right)}} + \frac{c{\left(t \right)} a'{\left(t \right)} c'{\left(t \right)}}{a{\left(t \right)}}$ | $\frac{\left(a{\left(t \right)} b{\left(t \right)} c''{\left(t \right)} + a{\left(t \right)} b'{\left(t \right)} c'{\left(t \right)} + b{\left(t \right)} a'{\left(t \right)} c'{\left(t \right)}\right) c{\left(t \right)}}{a{\left(t \right)} b{\left(t \right)}}$ |

逆計量と縮約すると Ricci スカラーとなる。さらに $G_{ab}=R_{ab}-g_{ab}R/2$ を成分ごとに計算する。

$$
\begin{aligned}
R&=\sum_{a,b}g^{ab}R_{ab}\\
&=\frac{a{\left(t \right)} b{\left(t \right)} c''{\left(t \right)} + a{\left(t \right)} c{\left(t \right)} b''{\left(t \right)} + b{\left(t \right)} c{\left(t \right)} a''{\left(t \right)}}{a{\left(t \right)} b{\left(t \right)} c{\left(t \right)}} + \frac{a{\left(t \right)} b{\left(t \right)} c''{\left(t \right)} + a{\left(t \right)} b'{\left(t \right)} c'{\left(t \right)} + b{\left(t \right)} a'{\left(t \right)} c'{\left(t \right)}}{a{\left(t \right)} b{\left(t \right)} c{\left(t \right)}} + \frac{a{\left(t \right)} c{\left(t \right)} b''{\left(t \right)} + a{\left(t \right)} b'{\left(t \right)} c'{\left(t \right)} + c{\left(t \right)} a'{\left(t \right)} b'{\left(t \right)}}{a{\left(t \right)} b{\left(t \right)} c{\left(t \right)}} + \frac{b{\left(t \right)} c{\left(t \right)} a''{\left(t \right)} + b{\left(t \right)} a'{\left(t \right)} c'{\left(t \right)} + c{\left(t \right)} a'{\left(t \right)} b'{\left(t \right)}}{a{\left(t \right)} b{\left(t \right)} c{\left(t \right)}}\\
&=\frac{2 \left(a{\left(t \right)} b{\left(t \right)} c''{\left(t \right)} + a{\left(t \right)} c{\left(t \right)} b''{\left(t \right)} + a{\left(t \right)} b'{\left(t \right)} c'{\left(t \right)} + b{\left(t \right)} c{\left(t \right)} a''{\left(t \right)} + b{\left(t \right)} a'{\left(t \right)} c'{\left(t \right)} + c{\left(t \right)} a'{\left(t \right)} b'{\left(t \right)}\right)}{a{\left(t \right)} b{\left(t \right)} c{\left(t \right)}}
\end{aligned}
$$

非零の Einstein テンソル成分は次のとおりである。

| 成分 | Ricci 項 | 跡の項 | 差 |
|---|---|---|---|
| $G_{t t}$ | $- \frac{a{\left(t \right)} b{\left(t \right)} c''{\left(t \right)} + a{\left(t \right)} c{\left(t \right)} b''{\left(t \right)} + b{\left(t \right)} c{\left(t \right)} a''{\left(t \right)}}{a{\left(t \right)} b{\left(t \right)} c{\left(t \right)}}$ | $- \frac{a{\left(t \right)} b{\left(t \right)} c''{\left(t \right)} + a{\left(t \right)} c{\left(t \right)} b''{\left(t \right)} + a{\left(t \right)} b'{\left(t \right)} c'{\left(t \right)} + b{\left(t \right)} c{\left(t \right)} a''{\left(t \right)} + b{\left(t \right)} a'{\left(t \right)} c'{\left(t \right)} + c{\left(t \right)} a'{\left(t \right)} b'{\left(t \right)}}{a{\left(t \right)} b{\left(t \right)} c{\left(t \right)}}$ | $\frac{a{\left(t \right)} b'{\left(t \right)} c'{\left(t \right)} + b{\left(t \right)} a'{\left(t \right)} c'{\left(t \right)} + c{\left(t \right)} a'{\left(t \right)} b'{\left(t \right)}}{a{\left(t \right)} b{\left(t \right)} c{\left(t \right)}}$ |
| $G_{x x}$ | $\frac{\left(b{\left(t \right)} c{\left(t \right)} a''{\left(t \right)} + b{\left(t \right)} a'{\left(t \right)} c'{\left(t \right)} + c{\left(t \right)} a'{\left(t \right)} b'{\left(t \right)}\right) a{\left(t \right)}}{b{\left(t \right)} c{\left(t \right)}}$ | $\frac{\left(a{\left(t \right)} b{\left(t \right)} c''{\left(t \right)} + a{\left(t \right)} c{\left(t \right)} b''{\left(t \right)} + a{\left(t \right)} b'{\left(t \right)} c'{\left(t \right)} + b{\left(t \right)} c{\left(t \right)} a''{\left(t \right)} + b{\left(t \right)} a'{\left(t \right)} c'{\left(t \right)} + c{\left(t \right)} a'{\left(t \right)} b'{\left(t \right)}\right) a{\left(t \right)}}{b{\left(t \right)} c{\left(t \right)}}$ | $- \frac{\left(b{\left(t \right)} c''{\left(t \right)} + c{\left(t \right)} b''{\left(t \right)} + b'{\left(t \right)} c'{\left(t \right)}\right) a^{2}{\left(t \right)}}{b{\left(t \right)} c{\left(t \right)}}$ |
| $G_{y y}$ | $\frac{\left(a{\left(t \right)} c{\left(t \right)} b''{\left(t \right)} + a{\left(t \right)} b'{\left(t \right)} c'{\left(t \right)} + c{\left(t \right)} a'{\left(t \right)} b'{\left(t \right)}\right) b{\left(t \right)}}{a{\left(t \right)} c{\left(t \right)}}$ | $\frac{\left(a{\left(t \right)} b{\left(t \right)} c''{\left(t \right)} + a{\left(t \right)} c{\left(t \right)} b''{\left(t \right)} + a{\left(t \right)} b'{\left(t \right)} c'{\left(t \right)} + b{\left(t \right)} c{\left(t \right)} a''{\left(t \right)} + b{\left(t \right)} a'{\left(t \right)} c'{\left(t \right)} + c{\left(t \right)} a'{\left(t \right)} b'{\left(t \right)}\right) b{\left(t \right)}}{a{\left(t \right)} c{\left(t \right)}}$ | $- \frac{\left(a{\left(t \right)} c''{\left(t \right)} + c{\left(t \right)} a''{\left(t \right)} + a'{\left(t \right)} c'{\left(t \right)}\right) b^{2}{\left(t \right)}}{a{\left(t \right)} c{\left(t \right)}}$ |
| $G_{z z}$ | $\frac{\left(a{\left(t \right)} b{\left(t \right)} c''{\left(t \right)} + a{\left(t \right)} b'{\left(t \right)} c'{\left(t \right)} + b{\left(t \right)} a'{\left(t \right)} c'{\left(t \right)}\right) c{\left(t \right)}}{a{\left(t \right)} b{\left(t \right)}}$ | $\frac{\left(a{\left(t \right)} b{\left(t \right)} c''{\left(t \right)} + a{\left(t \right)} c{\left(t \right)} b''{\left(t \right)} + a{\left(t \right)} b'{\left(t \right)} c'{\left(t \right)} + b{\left(t \right)} c{\left(t \right)} a''{\left(t \right)} + b{\left(t \right)} a'{\left(t \right)} c'{\left(t \right)} + c{\left(t \right)} a'{\left(t \right)} b'{\left(t \right)}\right) c{\left(t \right)}}{a{\left(t \right)} b{\left(t \right)}}$ | $- \frac{\left(a{\left(t \right)} b''{\left(t \right)} + b{\left(t \right)} a''{\left(t \right)} + a'{\left(t \right)} b'{\left(t \right)}\right) c^{2}{\left(t \right)}}{a{\left(t \right)} b{\left(t \right)}}$ |

アフィン径数の微分を点で表し、$w_a$ を座標 $a$ の速度成分とする。ここでは添字を下げる操作を行っていない。測地線の右辺 $-\Gamma^a{}_{bc}w_bw_c$ では $b\ne c$ の項を二度数える。

$$
\begin{aligned}
\ddot{t}&=- w_{x}^{2} a{\left(t \right)} a'{\left(t \right)} - w_{y}^{2} b{\left(t \right)} b'{\left(t \right)} - w_{z}^{2} c{\left(t \right)} c'{\left(t \right)}\\
\ddot{x}&=- \frac{2 w_{t} w_{x} a'{\left(t \right)}}{a{\left(t \right)}}\\
\ddot{y}&=- \frac{2 w_{t} w_{y} b'{\left(t \right)}}{b{\left(t \right)}}\\
\ddot{z}&=- \frac{2 w_{t} w_{z} c'{\left(t \right)}}{c{\left(t \right)}}
\end{aligned}
$$

また計量ノルムは

$$
g(\dot x,\dot x)=- w_{t}^{2} + w_{x}^{2} a^{2}{\left(t \right)} + w_{y}^{2} b^{2}{\left(t \right)} + w_{z}^{2} c^{2}{\left(t \right)}
$$
である。これを径数で微分して上の加速度を代入すると、係数の微分と速度の微分の和が零になる。従って初期のノルムが保存される。

表の $G_{tt}$ は三組 $(a'/a)(b'/b),(a'/a)(c'/c),(b'/b)(c'/c)$ の和であり、Gauss の二重縮約の半分に等しい。例えば $a=b=c$ では同じ積が三つあり $3H^2$ になる。逆に $a=t,b=c=1$ では一つの外的曲率だけが非零で、外的項は零である。後者は $(T,X)=(t\cosh x,t\sinh x)$ の座標変換をしただけなので、直接計算の全曲率も零である。
:::

### 法線を二本含む曲率と測地線の束

Gaussian 時刻 $t$ で $g=-\dd t^2+h_{ij}\dd x^i\dd x^j$ とする。$R^0{}_{i0j}$ の微分項は $\dot K_{ij}$、積項は $-K_{jk}K^k{}_i$ である。第一添字を下ろして
$$
R_{0i0j}=-\dot K_{ij}+K_{ik}K^k{}_j
$$
となる。添字を上げた $K^i{}_j$ の時間変化には $\dot h^{ij}=-2K^{ij}$ が寄与するため、同じ式は $R^i{}_{0j0}=-\partial_tK^i{}_j-K^i{}_kK^k{}_j$ である。この二式を比べれば、時間微分が添字を上げる操作と可換でないことを確認できる。

例えば $h_{ij}=a(t)^2\delta_{ij}$ なら $R_{0i0j}=-(\ddot a/a)h_{ij}$ である。$a=e^{Ht}$ の場合は $-\ddot a/a=-H^2$ となる一方、Gauss の空間成分は $H^2(h_{ik}h_{jl}-h_{il}h_{jk})$ となる。両者が同じ定曲率 $H^2$ を表すことは、計量の時間成分が $-1$ であることを使うと分かる。

### 拘束と通常の空間曲率を区別する

時間的な超曲面も含めて接方向の交換子を計算すると、$R^{\mathrm{proj}}_{ijkl}=R^{(h)}_{ijkl}-\epsilon(K_{ik}K_{jl}-K_{il}K_{jk})$ である。符号は $\nabla_XY$ に含まれる $-\epsilon K(X,Y)s$ から来る。特に、平坦な時空に埋め込まれた時間的円筒 $r=R$ では $K^i{}_j$ の固有値は $(0,1/R,1/R)$ である。従って円筒の時間・角方向の混合曲率は零、二つの直交角方向の断面曲率は $1/R^2$ となる。誘導計量 $-\dd t^2+R^2\dd\Omega_2^2$ の直積構造から計算しても同じ結果である。

空間的な面について $\epsilon=-1$ とし、$K=0$ を仮定する。真空で $\Lambda=0$ なら縮約した Gauss 式が課す条件は $R^{(3)}=0$ である。これは $R^{(3)}_{ij}=0$ の六成分をまとめた条件ではない。例えば三次元で Ricci 固有値が局所的に $(\lambda,-\lambda,0)$ であれば、スカラー曲率が零でも断面曲率は残る。直交枠で対角成分を $R^{(3)}_{11}=k_{12}+k_{13}$、$R^{(3)}_{22}=k_{12}+k_{23}$、$R^{(3)}_{33}=k_{13}+k_{23}$ と書けば、$k_{12}=0,k_{13}=\lambda,k_{23}=-\lambda$ を得る。従って拘束から空間の各成分を決めるには追加の条件が要る。

### 運動量拘束の独立な計算

空間依存する $K$ を持つ簡単な例として $g=-\dd t^2+a(t,x)^2(\dd x^2+\dd y^2+\dd z^2)$ を取る。$K^j{}_i=H(t,x)\delta^j_i$、$H=\partial_t\log a$ である。$D_jK^j{}_i=\partial_i H$ と $D_iK=3\partial_i H$ から $G_{0i}=-2\partial_i H$ を得る。時空の接続で直接確かめるには、$R_{0i}$ の第一の微分項 $\partial_j\Gamma^j{}_{i0}$ が $\partial_i H$、第二の微分項 $-\partial_i\Gamma^j{}_{j0}$ が $-3\partial_i H$ となることを使う。積項では $H\Gamma^j{}_{ji}$ が二つ現れて相殺する。$g_{0i}=0$ なので跡を引いても値は同じである。

従って例えば $a=e^{btx}$ では $G_{0x}=-2b$、この法線観測者に対する $j_x=b/(4\pi G)$ となる。$b\ne0$ なら、各点で外的曲率が等方的であることだけから物質の運動量密度を零とはできない。

## Einstein–Hilbert 作用と変分原理

:term[Einstein–Hilbert 作用]{en="Einstein–Hilbert action" yomi="あいんしゅたいんひるべるとさよう"} を $S_g=(16\pi G)^{-1}\int\sqrt{-g}(R-2\Lambda)\dd^4x$ とする。逆計量の変分に対し $\delta g_{ab}=-g_{ac}g_{bd}\delta g^{cd}$、行列式の微分から $\delta\sqrt{-g}=-\sqrt{-g}g_{ab}\delta g^{ab}/2$ である。:eqref[eq-gr-curvature-components] を微分して、接続差の共変微分にまとめると

$$
\delta R_{ab}=\nabla_c\delta\Gamma^c{}_{ab}-\nabla_b\delta\Gamma^c{}_{ac},\qquad
g^{ab}\delta R_{ab}=\nabla_c V^c,\quad V^c=g^{ab}\delta\Gamma^c{}_{ab}-g^{ac}\delta\Gamma^b{}_{ab}.
$$

従って内部の係数は $(G_{ab}+\Lambda g_{ab})/(16\pi G)$ である。物質作用の変分を $\delta S_m=-\tfrac12\int\sqrt{-g}T_{ab}\delta g^{ab}$ と定義すれば、コンパクトな台の変分から :eqref[eq-gr-einstein-field] を得る。例えばスカラー作用の $\delta g^{ab}$ は運動項の縮約と体積密度に作用し、$T_{ab}=\partial_a\phi\partial_b\phi-g_{ab}((\nabla\phi)^2+m^2\phi^2)/2$ になる。物質の運動方程式が成り立つとき、微小微分同相での $\delta g^{ab}=-2\nabla^{(a}\xi^{b)}$ を代入し $\xi$ の微分を部分積分すると $\nabla^aT_{ab}=0$ である。これは :ref[thm-am-field-noether]{full} に対応する保存則である。

### 共形変分で跡方程式を取り出す

$\delta g_{ab}=2f g_{ab}$、従って $\delta g^{ab}=-2f g^{ab}$ とする。四次元では $\delta\sqrt{-g}=4f\sqrt{-g}$ である。接続の変分は $\delta\Gamma^c{}_{ab}=\delta^c_a\partial_bf+\delta^c_b\partial_af-g_{ab}\nabla^cf$ なので、縮約した $V^c$ は $-6\nabla^cf$ となる。従って $\delta R=-2fR-6\Box_g f$ である。$f$ の台が内部にあるなら $\Box_g f$ の積分は零で、重力部分の変分は $(16\pi G)^{-1}\int \sqrt{-g}f(2R-8\Lambda)$ になる。物質部分は $\int \sqrt{-g} fT$ なので、この変分だけで $R=4\Lambda-8\pi G T$ が従う。他の成分を得るには計量の跡を変えない変分も必要である。

### 接続変分を計量だけで表す

同じ座標で二つの Levi-Civita 接続を比べると、接続係数の差はテンソルである。$v_{ab}=\delta g_{ab}$、$v=g^{ab}v_{ab}$ とする。計量整合性 $\nabla_a g_{bc}=0$ を微分すると $\nabla_a v_{bc}=\delta\Gamma^d{}_{ab}g_{dc}+\delta\Gamma^d{}_{ac}g_{bd}$ である。この式で $a,b,c$ を巡回させ、最初の二式を加えて第三式を引くと、ねじれ零の対称性を使って
$$
\delta\Gamma^c{}_{ab}=\tfrac12 g^{cd}(\nabla_a v_{bd}+\nabla_b v_{ad}-\nabla_d v_{ab})
$$
を得る。$g^{ab}\delta\Gamma^c{}_{ab}=\nabla_a v^{ac}-\nabla^cv/2$、$\delta\Gamma^b{}_{ab}=\nabla_av/2$ であるから、発散のベクトルは $V^c=\nabla_a v^{ac}-\nabla^cv$ である。従って $v_{ab}$ を地域の境界で零にしても $V^c$ の法線成分は残り得る。逆に $v_{ab}$ が内部に台を持つなら $\int\nabla_cV^c$ は零であり、基本補題を対称成分ごとに適用できる。

### 保存則を作用の対称性から得る

コンパクトな台を持つ $\xi$ による微分同相を考える。重力作用は密度の積分であり、その Lie 微分は全発散となる。$E^g_{ab}=G_{ab}+\Lambda g_{ab}$ と書けば、逆計量変分 $\delta g^{ab}=-2\nabla^{(a}\xi^{b)}$ に対して
$$
0=-\frac{1}{8\pi G}\int\sqrt{-g}\,E^g_{ab}\nabla^a\xi^b\,\dd^4x
=\frac{1}{8\pi G}\int\sqrt{-g}\,\xi^b\nabla^a E^g_{ab}\,\dd^4x .
$$
任意の $\xi$ に対して成り立つので、$\nabla^a E^g_{ab}=0$ は運動方程式を用いずに成立する。例えば最小結合スカラー場では、物質の変分が運動方程式を満たさないときでも直接微分して
$$
\nabla^aT_{ab}=(\Box_g\phi-m^2\phi)\partial_b\phi
$$
を得る。実際、$\nabla^a\phi\nabla_a\nabla_b\phi$ は運動項の跡の微分と相殺する。従って場の方程式を使うと保存則が従い、重力側の恒等式と整合する。逆に $\nabla T=0$ からスカラー方程式を読む場合、$\dd\phi=0$ の開集合上では積を割れないので、物質変分自体が必要である。

:::example[定数スカラー場の有効な宇宙項]
物質ラグランジアンを $-\tfrac12(\nabla\phi)^2-U(\phi)$ とし、$\phi=\phi_0$ が $U'(\phi_0)=0$ を満たすとする。運動項は零、体積密度の変分だけが寄与して $T_{ab}=-U(\phi_0)g_{ab}$ である。従って重力方程式の左辺に移せば $\Lambda_{\rm eff}=\Lambda+8\pi G U(\phi_0)$ となる。例えば平坦切片の指数膨張 $a(t)=e^{Ht}$ では $G_{ab}=-3H^2g_{ab}$ であり、$3H^2=\Lambda_{\rm eff}$ のとき解になる。$U'(\phi_0)\ne0$ の定数場でも応力だけは発散零であるが、独立な物質変分はその値を許さない。
:::

### コンパクトな台の変分で境界を避ける

:::example[コンパクトな台の共形変分で跡を照合する]
平坦時空の計量を $e^{2\epsilon f}\eta_{ab}$ で変分し、$f$ は滑らかでコンパクトな台を持つとする。一次の $R$ は $-6\epsilon\Box f$ で、作用の一次変分はこの全発散の積分になり零である。一方、宇宙項を入れると $\sqrt{-g}=1+4\epsilon f+\cdots$ なので $\delta S_\Lambda=-8\Lambda\int f/(16\pi G)$ となる。$\Lambda\ne0$ なら任意の $f$ について零にならない。従って物質がない平坦時空は、宇宙項が非零の作用の停留点にならない。
:::

## Gibbons–Hawking–York 境界項

滑らかな非光的境界の外向き単位法線を $s$、$\epsilon=s\cdot s=\pm1$ とする。境界を固定した法線座標 $g=\epsilon\dd y^2+h_{ij}\dd x^i\dd x^j$ で $s=\partial_y$ とする。誘導計量を固定する変分は $\delta h_{ij}=0$ だが、$\partial_y\delta h_{ij}$ は任意である。法線座標条件を保つ代表を取ると、$\delta\Gamma^y{}_{ij}=-\epsilon\partial_y\delta h_{ij}/2$、$\delta\Gamma^i{}_{yi}=h^{ij}\partial_y\delta h_{ij}/2$ である。従って $V^y=-2\epsilon\delta K$、外向き流束は $\epsilon s_aV^a=V^y$ となる。法線座標条件は境界で恒等の微分同相を使って変分した計量ごとに整えられる。よって

$$
S_{\rm GHY}=\frac1{8\pi G}\int_{\partial M}\epsilon\sqrt{\abs{\det h}}\,K\,\dd^3x
$$

を加えると法線微分が相殺する。これを :term[Gibbons–Hawking–York 境界項]{en="Gibbons–Hawking–York term" yomi="ぎぼんずほーきんぐよーくきょうかいこう"} という。例えば Minkowski 時空で固定半径 $r=R$ の円筒に外向き法線を取ると角方向二つが $1/R$ ずつ、時間方向が零なので $K=2/R$ である。角や零境界がある場合はこの法線座標の議論を直接適用できない。

### 平坦な切片の作用から端の加速度を除く

有限な共動体積 $V_0$ を持つ平坦な空間トーラスで $g=-\dd t^2+a(t)^2\delta_{ij}\dd x^i\dd x^j$ とする。$R=6(\ddot a/a+\dot a^2/a^2)$ なので $a^3R=6[\dd(a^2\dot a)/\dd t-a\dot a^2]$ である。未来端の法線は $+\partial_t$、過去端は $-\partial_t$ で、$\epsilon=-1$、$K=\pm3\dot a/a$ を GHY 項に入れると端で
$$
S_{\rm GHY}=-\frac{6V_0}{16\pi G}[a^2\dot a]_{t_1}^{t_2}
$$
を得る。Einstein–Hilbert 項の全微分と相殺するため、端の変分には $\delta\dot a$ が残らない。残った作用は $-6V_0\int a\dot a^2\dd t/(16\pi G)$ である。lapse を先に $1$ と固定したので、この作用を $a$ だけで変分しても lapse の拘束は復元できない。

### 円筒の法線を内外で変更する

:::example[平坦な球殻領域で境界項の符号を比べる]
Minkowski で $R_1<r<R_2$、$0<t<T$ の領域を取る。動径境界の法線は外壁で $+\partial_r$、内壁で $-\partial_r$ だから、$K$ はそれぞれ $+2/R_2,-2/R_1$ である。両方で $\epsilon=+1$、面の体積は $4\pi R_i^2 T$ なので、動径の二面の GHY 項の和は $T(R_2-R_1)/G$ である。初めと終わりの平坦な時刻面では $K=0$ である。この計算は滑らかな各面からの寄与を求めたものであり、一般の境界の接合では二面の法線の相対角度も境界変分に入る。
:::

### 境界に残る有限の変分

一般の $\delta h_{ij}$ を許し、先と同じ Gaussian 代表で $v_{yy}=v_{yi}=0$ を取る。$K_{ij}=\partial_yh_{ij}/2$、$\delta h^{ij}=-v^{ij}$ だから
$$
\delta K=\tfrac12h^{ij}\partial_yv_{ij}-K^{ij}v_{ij},\qquad
V^y=\epsilon(K^{ij}v_{ij}-h^{ij}\partial_yv_{ij}).
$$
第二式は $\nabla_a v^{ay}=-\epsilon K_{ij}v^{ij}$ と $\nabla^yv=\epsilon(h^{ij}\partial_yv_{ij}-2K^{ij}v_{ij})$ の差である。さらに $\delta\sqrt{|\det h|}=\sqrt{|\det h|}h^{ij}v_{ij}/2$ を使うと、境界上の全変分は
$$
\delta(S_g+S_{\rm GHY})\big|_{\partial M}
=\frac{1}{16\pi G}\int_{\partial M}\epsilon\sqrt{|\det h|}
(Kh^{ij}-K^{ij})v_{ij}\,\dd^3x.
$$
特に未来の時刻面では $\epsilon=-1$ なので、係数は $\sqrt h(K^{ij}-Kh^{ij})/(16\pi G)$ である。これは次節以降で速度の偏微分から得る正準運動量と一致する。過去面では外向き法線の反転によって逆符号になるため、二つの端の寄与は通常の力学と同様に終点と始点の差になる。

### 接合面における境界項と薄い殻

時間的な面 $y=0$ に沿って二つの時空を接着し、$h_{ij}$ は連続だが $\partial_y h_{ij}$ に有限な跳びがあるとする。法線を両側で同じ $+\partial_y$ とし、$[K_{ij}]=K_{ij}(0^+)-K_{ij}(0^-)$ と定義する。両側の作用に外向き法線による GHY 項を付けると、境界変分の和では左の法線が $+\partial_y$、右の法線が $-\partial_y$ なので、係数は $([K^{ij}]-[K]h^{ij})/(16\pi G)$ になる。面の物質作用が $\delta S_{\rm shell}=\tfrac12\int\sqrt{|h|}S^{ij}v_{ij}$ と変分される場合、面での停留条件は
$$
8\pi G S_{ij}=-[K_{ij}]+[K]h_{ij}.
$$
これは :term[Israel 接合条件]{en="Israel junction conditions" yomi="いすらえるせつごうじょうけん"} である。連続な計量を分布として微分しても、$R_{ij}$ の特異部 $-[K_{ij}]\delta(y)$、$R_{yy}$ の特異部 $-[K]\delta(y)$ から同じ $G_{ij}$ を得る。実際 $\Gamma^y{}_{ij}=-K_{ij}$、$\Gamma^i{}_{yj}=K^i{}_j$ の $y$ 微分だけがデルタ関数を含み、$R$ の特異部は $-2[K]\delta(y)$ である。接合が滑らかで $[K_{ij}]=0$ なら、二つの境界変分は消え、薄い殻の物質を置く必要はない。

:::example[静的な球殻の面密度]
面積半径 $R>2m$ の静止殻の内側を Minkowski、外側を Schwarzschild とする。両側の時刻を面の固有時に合わせると、誘導計量は $-\dd\tau^2+R^2\dd\Omega_2^2$ で連続である。外向きに共通の向きの法線を取り、外側で $s^r=\sqrt{1-2m/R}$、内側で $s^r=1$ である。角方向は $K^\theta{}_\theta=K^\phi{}_\phi=s^r/R$ であるから、$\sigma=-S^\tau{}_\tau=-(1/(4\pi G))[K^\theta{}_\theta]$ は
$$
\sigma=\frac{1-\sqrt{1-2m/R}}{4\pi G R}
$$
となる。殻の固有質量 $M_0=4\pi R^2\sigma$ を使って両辺を二乗すると $m/G=M_0-GM_0^2/(2R)$ である。外から測る質量との差は負の結合エネルギーである。この構成では静止を保つ面圧も必要であり、時間成分の跳びを含めると $S^\theta{}_\theta=([K^\tau{}_\tau]+[K^\theta{}_\theta])/(8\pi G)$ となる。外の $K^\tau{}_\tau=f'(R)/(2\sqrt f)$ と内の零を代入すれば面圧も一意に決まる。
:::

### 滑らかでない境界における接合の寄与

非光的な二つの境界面が角で交わる場合、法線座標は角を横切って一つには取れない。例えば時間的な角断面を持つ二つの時間的壁では、角に直交する二平面が正定値であり、壁を小さい円弧で丸めることができる。円弧の半径を $\delta$、外向き法線の符号付き回転角を $\eta$ とする。丸めた部分で円弧方向の外的曲率が符号付きで $1/\delta$、固有長要素が $\delta\,\dd\eta$ なので、その積分の極限は $\eta$ である。角に沿った曲率が有界なら、その寄与は $\delta\to0$ で消える。従って $\epsilon=+1$ のこの場合、外向き法線の回転角による接合項は $(8\pi G)^{-1}\int\sqrt{|\det q|}\,\eta\,\dd^2x$ となる。ここで $q$ は角に誘導される計量である。例えば直角の凸な壁なら $\eta=\pi/2$、滑らかにつながる面なら $\eta=0$ である。角に直交する平面が Lorentz 署名の場合は回転角の代わりに boost の急速度が必要になり、光的接合には法線の規格化も関わる。これらを同じユークリッド円弧で代用しない。非光的面を組み合わせた作用の取り扱いは :cite[poisson-toolkit] を参照する。

## Palatini 形式

:term[Palatini 形式]{en="Palatini formulation" yomi="ぱらてぃーにけいしき"} では $g$ とねじれ零の $\Gamma$ を独立に変分する。物質は接続に結合しないとし、$P^{ab}=\sqrt{-g}g^{ab}$ とする。上の $\delta R$ の式から部分積分すると、$\delta\Gamma^c{}_{ab}$ の下二添字が対称なので

$$
\nabla_cP^{ab}-\delta_c^{(a}\nabla_dP^{b)d}=0
$$

を得る。この積分には重み一のベクトル密度の共変発散が偏発散に等しいことを使った。$c=a$ の縮約は $(1-n)\nabla_dP^{bd}/2=0$ を与える。$n=4$ では $\nabla_cP^{ab}=0$ となる。$g_{ab}$ でさらに縮約すると、行列式に対する密度の微分から $(1-n/2)g_{ab}\nabla_cg^{ab}=0$ が従う。従って $n=4$ では $\nabla g=0$ であり、ねじれ零の接続は Levi-Civita 接続になる。計量変分は独立接続の Ricci の対称部分を与え、接続の解を代入すると通常の方程式に戻る。

### 二次元への縮約の例外

$P^{ab}=\sqrt{-g}g^{ab}$ の重み一の共変微分は
$$
\nabla_cP^{ab}=\sqrt{-g}\left(\nabla_cg^{ab}-\tfrac12g^{ab}g_{ij}\nabla_cg^{ij}\right)
$$
である。計量の行列式に Jacobi の公式を使い、密度に対する接続の跡を引くとこの式になる。従って $n=2$ では $\nabla_cg^{ab}=W_c g^{ab}$ であっても $\nabla_cP^{ab}=0$ である。具体的に Levi-Civita 接続に $D^a{}_{bc}=\delta^a_b A_c+\delta^a_c A_b-g_{bc}A^a$ を足すと $\nabla_c g^{ab}=2A_cg^{ab}$ となる。従って $n=2$ では任意の $A$ がこの条件を満たし、四次元での一意性をそのまま主張できない。

### 同じ計量を二つの方法で微分する

:::example[重み一の逆計量を保つ接続]
二次元平坦時空に $A=\lambda\,\dd x$ を取り、Levi-Civita 接続に $D$ を足す。もとの計量は座標で定数であるが、$D$ は $\nabla_x g^{tt}=2\lambda g^{tt}$ を作る。同時に接続の跡は $D^i{}_{xi}=2\lambda$ なので、この密度重みの補正が $\nabla_x P^{tt}$ で相殺する。四次元で同じ形の $D$ を採用すると接続の跡は $4\lambda$ で、$P$ には $-2\lambda P^{tt}$ が残る。このため二次元の非一意性は高次元の Palatini 接続には持ち越されない。
:::

### 接続を自由にした場合の射影対称性

先の証明は接続の下二添字の対称性を使った。ここでは $\nabla_c v^a=\partial_cv^a+\Gamma^a{}_{cb}v^b$ の順で添字を書き、ねじれを許す。任意の一形式 $A_c$ による $\Gamma^a{}_{cb}\mapsto\Gamma^a{}_{cb}+\delta^a_b A_c$ を考える。行列 $(\Gamma_c)^a{}_b$ で見ると $A_c I$ の追加である。曲率 $\partial_c\Gamma_d-\partial_d\Gamma_c+[\Gamma_c,\Gamma_d]$ の交換子は変化せず、微分項は $I(\partial_cA_d-\partial_dA_c)$ だけ増える。このため Ricci の変化は $\partial_b A_d-\partial_d A_b$ という反対称部分であり、$g^{bd}R_{bd}$ は変わらない。従って接続に結合しない物質と計量による作用だけでこの自由度を固定することはできない。これを :term[射影対称性]{en="projective symmetry" yomi="しゃえいたいしょうせい"} という。

例えば Minkowski 時空で元の接続を零とし、$A=y\,\dd x$ を追加すると曲率二形式は $\dd y\wedge\dd x$ の恒等行列倍になる。Ricci の反対称部分は非零だが、スカラー曲率は零のままである。自己平行線の方程式には $A(\dot x)\dot x^a$ が足される。増えた加速度は速度と平行なので、アフィン径数を変えることで同じ曲線の像を得る。これが対称性の名称の由来である。この例の接続では一般にねじれが非零であるため、先の対称接続を独立変分した結論と両立する。

### 曲率を非線形に含むと何が変わるか

計量と対称接続を独立としたまま、作用の $R$ を $F(\mathcal R)$ に置き換える。ただし $\mathcal R=g^{ab}R_{ab}(\Gamma)$ とし、$F'(\mathcal R)>0$ の枝だけを考える。計量変分は
$$
F'(\mathcal R)R_{(ab)}-\tfrac12F(\mathcal R)g_{ab}=8\pi G T_{ab},
\qquad F'(\mathcal R)\mathcal R-2F(\mathcal R)=8\pi G T .
$$
接続を変分する部分積分では $P^{ab}=\sqrt{-g}F' g^{ab}$ が現れる。四次元で $\hat g_{ab}=F' g_{ab}$ とすれば $\sqrt{-\hat g}\,\hat g^{ab}=\sqrt{-g}F' g^{ab}$ となり、先の縮約から $\Gamma$ は $\hat g$ の Levi-Civita 接続である。例えば $F(\mathcal R)=\mathcal R+b\mathcal R^2$ なら跡は $-\mathcal R=8\pi G T$ となり、$T$ が空間で変わると $\hat g$ の共形因子も変わる。一方、真空では $\mathcal R=0$、$F'=1$ であり、通常の真空解に戻る。接続を始めから $g$ の Levi-Civita 接続とした曲率二乗作用では、部分積分で $\nabla_a\nabla_bF'$ が残り、同じ代数的跡方程式にはならない。

## 3+1 分解と ADM 形式

:term[ADM 分解]{en="ADM decomposition" yomi="えーでぃーえむぶんかい"} は $\partial_t=Nn+\beta^i\partial_i$ と分ける。$N>0$ を lapse、$\beta^i$ を shift と呼ぶ。前者は法線の固有時間、後者は座標格子の接方向変位を表す。

$$
\dd s^2=-N^2\dd t^2+h_{ij}(\dd x^i+\beta^i\dd t)(\dd x^j+\beta^j\dd t),\qquad
K_{ij}=(\dot h_{ij}-D_i\beta_j-D_j\beta_i)/(2N).
$$

枠の変換が三角なので $\sqrt{-g}=N\sqrt h$ である。法線加速度 $a=\nabla_n n$ を用いる。$\nabla_a\nabla_b n^a-\nabla_b\nabla_a n^a=R_{db}n^d$ を $n^b$ と縮約し、積の微分則を使うと $R_{nn}=-n(K)-K_{ij}K^{ij}+\nabla_a a^a$ である。Gauss の式と $\nabla_a(Kn^a)=n(K)+K^2$ により

$$
R=R^{(3)}+K_{ij}K^{ij}-K^2+2\nabla_a(Kn^a-a^a).
$$

時間切片の発散項を境界項で相殺し、空間境界を別に扱えば、バルクの作用密度は $N\sqrt h(R^{(3)}+K_{ij}K^{ij}-K^2-2\Lambda)/(16\pi G)$ となる。$\dot h$ で微分すると $\pi^{ij}=\sqrt h(K^{ij}-Kh^{ij})/(16\pi G)$、従って $\pi=-2K\sqrt h/(16\pi G)$ である。$\pi^{ij}\dot h_{ij}$ に $\dot h=2NK+\mathcal L_\beta h$ を入れて shift の微分を部分積分すると、Hamilton 密度は $N\mathcal H_\perp+\beta^i\mathcal H_i$、ただし

$$
\mathcal H_\perp=\frac{16\pi G}{\sqrt h}(\pi_{ij}\pi^{ij}-\tfrac12\pi^2)
-\frac{\sqrt h}{16\pi G}(R^{(3)}-2\Lambda),\qquad \mathcal H_i=-2D_j\pi^j{}_i
$$

である。従って $N,\beta$ の変分は真空の拘束を課す。

### lapse を残した一変数の模型

空間が平坦で $h_{ij}=a^2\delta_{ij}$、shift を零とし、$N(t)$ を未定にする。$K^i{}_j=(\dot a/(Na))\delta^i_j$ であり、宇宙項だけがあるバルク作用は
$$
S=\frac{V_0}{16\pi G}\int\left(-\frac{6a\dot a^2}{N}-2\Lambda N a^3\right)\dd t .
$$
$N$ で変分して $6a\dot a^2/N^2-2\Lambda a^3=0$ を得る。従って $3(\dot a/(Na))^2=\Lambda$ であり、$N=1$ を選ぶのは拘束を得た後でよい。運動量は $p_a=-12(V_0/16\pi G)a\dot a/N$ である。$B=V_0/(16\pi G)$ とすると Legendre 変換は $N[-p_a^2/(24Ba)+2B\Lambda a^3]$ となり、$N$ は速度が入らない乗数である。例えば一定 $N=N_0$ なら膨張解は $a=a_*e^{N_0H(t-t_*)}$ で、固有時に換算すると $H=\sqrt{\Lambda/3}$ である。

### 逆計量と法線時間の求め方

ADM 計量で行列の時間・空間ブロックを展開すると
$$
g_{tt}=-N^2+h_{ij}\beta^i\beta^j,\qquad
g_{ti}=h_{ij}\beta^j,\qquad g_{ij}=h_{ij}.
$$
$h^{ij}$ を三次元の逆行列とする。候補 $g^{tt}=-N^{-2},g^{ti}=\beta^i/N^2,g^{ij}=h^{ij}-\beta^i\beta^j/N^2$ を掛けると、時間行×時間列で $(-N^2+\beta^2)(-N^{-2})+\beta^2/N^2=1$ となる。時間と空間の混合は相殺し、空間ブロックには $h_{ik}h^{kj}=\delta_i{}^j$ が残る。従って $n_a=-N(\dd t)_a$ の規格化は $g^{tt}N^2=-1$ であり、$n^a=N^{-1}(\partial_t-\beta^i\partial_i)^a$ になる。

:::example[座標の静止と切片の法線を区別する]
平坦な二次元で $\dd s^2=-\dd t^2+(\dd x+B\dd t)^2$ とし、$B$ を定数とする。法線は $\partial_t-B\partial_x$ でノルム $-1$ である。座標を $X=x+Bt$ に戻すと法線は慣性時刻の $\partial_t|_X$ に一致する。一方、$x=\mathrm{const}$ の曲線は $X$ 方向の速度 $B$ を持ち、$|B|\ge1$ では時間的でない。lapse が常に $1$ でも静止した座標が時計を置ける軌道を与えるとは限らない。
:::

### 物質を含む Legendre 変換

最小結合スカラー場を例に、計量の拘束と物質の速度の関係を確かめる。$w=n^a\partial_a\phi=N^{-1}(\dot\phi-\beta^i\partial_i\phi)$ とするとラグランジアン密度は
$$
\mathcal L_m=N\sqrt h\bigl(\tfrac12w^2-\tfrac12h^{ij}\partial_i\phi\partial_j\phi-U(\phi)\bigr).
$$
従って $p_\phi=\partial\mathcal L_m/\partial\dot\phi=\sqrt h\,w$ である。$p_\phi\dot\phi$ に $\dot\phi=Np_\phi/\sqrt h+\beta^i\partial_i\phi$ を入れれば
$$
\begin{aligned}
\mathcal H_{\perp,m}
&=\frac{p_\phi^2}{2\sqrt h}+\sqrt h\left(\tfrac12 h^{ij}\partial_i\phi\partial_j\phi+U(\phi)\right),\\
\mathcal H_{i,m}&=p_\phi\partial_i\phi .
\end{aligned}
$$
ここで $\rho=(w^2+|D\phi|^2)/2+U$、$j_i=-w\partial_i\phi$ なので、$\mathcal H_{\perp,m}=\sqrt h\,\rho$、$\mathcal H_{i,m}=-\sqrt h\,j_i$ となる。重力と物質の Hamilton 密度を足して $\beta^i$ で変分すれば $-2D_j\pi^j{}_i-\sqrt h\,j_i=0$ である。$\pi^j{}_i/\sqrt h=(K^j{}_i-K\delta^j_i)/(16\pi G)$ と $D_j\sqrt h=0$ から、Gauss–Codazzi の運動量拘束と符号まで一致する。

例えば空間トーラスで $\phi=\phi(t)$、$U=0$、shift 零なら、$\partial_t(a^3\dot\phi/N)=0$ である。これは作用が定数の $\phi$ のずれに不変であることからも従う。$N=1$ で $a^3\dot\phi=P$ と置くと $\rho=P^2/(2a^6)$ で、$\Lambda=0$ の拘束が $H^2=4\pi G P^2/(3a^6)$ を与える。膨張枝で $(a^3)^{\cdot}=\sqrt{12\pi G}|P|$ なので $a\propto(t-t_0)^{1/3}$ である。この解は $p=\rho$ の完全流体と同じ膨張率を持つ。

### 拘束は空間微分同相と法線変形を生成する

境界のない切片、または内部に台を持つ試験関数に限る。Poisson 括弧を
$$
\{F,G\}=\int_\Sigma\left(
\frac{\delta F}{\delta h_{ij}}\frac{\delta G}{\delta\pi^{ij}}
-\frac{\delta G}{\delta h_{ij}}\frac{\delta F}{\delta\pi^{ij}}\right)\dd^3x
$$
と定義し、$D[\beta]=\int\beta^i\mathcal H_i$、$H[N]=\int N\mathcal H_\perp$ と書く。$D[\beta]=\int\pi^{ij}\mathcal L_\beta h_{ij}$ は空間部分積分によって得られる。従って $\{h,D[\beta]\}=\mathcal L_\beta h$ であり、$\int \pi^{ij}\mathcal L_\beta(\delta h_{ij})=-\int(\mathcal L_\beta\pi^{ij})\delta h_{ij}$ から $\{\pi,D[\beta]\}=\mathcal L_\beta\pi$ も従う。後者では $\pi$ が重み一の密度であることを使う。

$\{D[\beta],D[\gamma]\}$ は、$\gamma$ を固定し $D[\gamma]$ の作用に沿って $D[\beta]$ を変分すれば計算できる。$\pi$ と $h$ の変分を使い、全 Lie 微分の積分を除くと $\int\pi^{ij}(\mathcal L_\beta\mathcal L_\gamma-\mathcal L_\gamma\mathcal L_\beta)h_{ij}$ が残る。従って $\{D[\beta],D[\gamma]\}=D[[\beta,\gamma]]$ である。同様に $\mathcal H_\perp$ はスカラー密度なので $\{H[N],D[\beta]\}=\int N\mathcal L_\beta\mathcal H_\perp=-H[\beta^i\partial_iN]$ となる。

最後に $\{H[N],H[M]\}$ を求める。速度に対応する変分は $\delta H[M]/\delta\pi^{ij}=2MK_{ij}$ である。$\delta H[N]/\delta h_{ij}$ のうち $N$ の導関数を持つ項は、空間曲率の変分を二度部分積分して
$$
-\frac{\sqrt h}{16\pi G}(D^iD^j N-h^{ij}\Delta_h N)
$$
となる。$N$ を微分しない項は $N,M$ を交換して引くと消える。$\pi^{ij}=\sqrt h(K^{ij}-Kh^{ij})/(16\pi G)$ を使うと
$$
\begin{aligned}
\{H[N],H[M]\}
&=2\int\pi^{ij}(N D_iD_jM-MD_iD_jN)\dd^3x\\
&=-2\int (D_i\pi^{ij})(N\partial_jM-M\partial_jN)\dd^3x\\
&=D[h^{ij}(N\partial_jM-M\partial_jN)] .
\end{aligned}
$$
二行目で $D_iN D_jM-D_iM D_jN$ が $\pi^{ij}$ と縮約して消える。右辺の係数に動的な $h^{ij}$ が入るのは、隣り合う二回の法線変形を比べるとき、面に沿った勾配の向きを誘導計量で決めるからである。

:::example[法線変形の交換を一次で計算する]
平坦な切片上で、小さな領域だけを考えて $N=1$、$M=x$ とする。規定した括弧の右辺は $D[\partial_x]$ である。$t=\delta s\,x$ へ一次だけ変形した面には、法線方向 $\partial_t+\delta s\,\partial_x$ が立つ。従って、先に面を傾けてから法線方向へ動かす場合には、法線を保ったまま動かす順序とは違う接方向の移動が生じる。切片の外まで試験関数を使うときは、これらの値を保つ領域の外で滑らかに零へ落とせば部分積分の仮定に従う。
:::

### 一次拘束と物理自由度の数

$N,\beta^i$ も座標として位相空間に含めると、その正準運動量 $p_N,p_{\beta,i}$ は零である。ラグランジアンに $\dot N,\dot\beta^i$ がないためである。これらが :term[一次拘束]{en="primary constraint" yomi="いちじこうそく"} である。全 Hamiltonian を $H[N]+D[\beta]+\int(\lambda p_N+\lambda^ip_{\beta,i})\dd^3x$ とすると、$\dot p_N=-\mathcal H_\perp$、$\dot p_{\beta,i}=-\mathcal H_i$ である。従って一次拘束を保存するために二次拘束 $\mathcal H_\perp=\mathcal H_i=0$ が出る。上で求めた括弧は二次拘束の線形結合になるため、さらに独立な拘束は生じない。

正則な点で局所自由度を数えるときは、$N,\beta^i$ を乗数として除いた十二次元の正準変数 $h_{ij},\pi^{ij}$ だけを考えてよい。各独立な第一級拘束は、制約面を一成分狭めるとともに、その Hamilton ベクトル場が張るゲージ方向を一成分与える。四つの拘束で $12-2\cdot4=4$ 個の位相空間成分が残り、配置自由度としては二つである。例えば平坦背景の線形重力波では、二つの偏光がこれに対応する。この数え方では局所的に独立な拘束を数えているため、境界の電荷や対称性のために作用が退化する特別なデータの大域自由度は別に扱う。

## 拘束方程式と初期値問題

$\rho=T_{nn}$、$j_i=-T_{ni}$ とすると :eqref[eq-gr-gauss-constraints] は

$$
R^{(3)}+K^2-K_{ij}K^{ij}=16\pi G\rho+2\Lambda,\qquad D_jK^j{}_i-D_iK=-8\pi G j_i .
\label{eq-gr-adm-constraints}
$$

接方向の方程式が $K_{ij}$ の時間発展を与える。$E_{ab}=G_{ab}+\Lambda g_{ab}-8\pi G T_{ab}$ とすれば、Bianchi と物質の保存から $\nabla_a E^{ab}=0$ である。例えば $E_{ij}=0$ を発展で課すと、Gaussian 座標での主部は $-\partial_tE_{00}+h^{ij}\partial_iE_{j0}=0$、$-\partial_tE_{0i}=0$ である。従って拘束の成分が線形同次の式を満たす。真空の初期値定理は、滑らかな制約付きデータ $(h,K)$ に対して最大大域双曲発展が存在し、初期面の埋め込みを保つ等長写像を除いて一意であると主張する。調和座標で準線形波動系に還元して局所解を構成し、幾何学的一意性を用いて発展を貼り合わせる。非線形波動系のエネルギー評価と最大発展の構成は :cite[wald-gr] を参照する。例えば $\Lambda=0$ の時間対称真空データでは $K=0$ なので $R^{(3)}=0$ が条件であり、空間の全曲率の消失までは必要ない。

### 外的曲率の発展式

まず $\beta=0$ の座標では $\Gamma^0{}_{ij}=K_{ij}/N$、$\Gamma^i{}_{0j}=NK^i{}_j$、$\Gamma^0{}_{0i}=\partial_i\log N$、$\Gamma^i{}_{00}=ND^iN$ となる。$\Gamma^0{}_{00}=\dot N/N$ と接続の空間成分も含めて $R_{ij}$ を縮約する。時間微分の $\partial_t(K_{ij}/N)$ に $\Gamma^0{}_{00}\Gamma^0{}_{ij}$ が加わり lapse の時間微分が相殺する。$\Gamma^0{}_{0i}$ の空間微分と積は $-D_iD_jN/N$、残りの時間添字を含む積は $KK_{ij}-2K_i{}^kK_{kj}$ となる。shift を戻すと時間微分が $\partial_t-\mathcal L_\beta$ に替わる。物質の接射影を $S_{ij}$、その跡を $S$ とすると

$$
\begin{aligned}
(\partial_t-\mathcal L_\beta)K_{ij}
&=D_iD_jN-N(R_{ij}^{(3)}+KK_{ij}-2K_i{}^kK_{kj})\\
&\quad+N\left[8\pi G(S_{ij}-\tfrac12h_{ij}(S-\rho))+\Lambda h_{ij}\right].
\end{aligned}
$$

最後の項は :eqref[eq-gr-trace-einstein] の接射影である。従ってここで課した残差は $E_{ij}-h_{ij}g^{ab}E_{ab}/2=0$ であり、拘束をまだ仮定しない段階では $E_{ij}=0$ と同じ条件ではない。例えば平坦な共動切片で $N=1,\beta=0,K_{ij}=Hh_{ij}$ とすれば、左辺は $(\dot H+2H^2)h_{ij}$、真空で右辺は $(-H^2+\Lambda)h_{ij}$ となる。拘束の $3H^2=\Lambda$ を使えば $\dot H=0$ であり、正の $\Lambda$ の指数膨張と一致する。

### 時間対称な共形初期データ

時間対称 $K=0$、$\Lambda=0$ のデータで $h_{ij}=\psi^4\delta_{ij}$ とする。三次元の共形変換で $f=2\log\psi$ を用いると $R^{(3)}=e^{-2f}(-4\Delta f-2|\grad f|^2)$ である。$\Delta f=2\Delta\psi/\psi-2|\grad\psi|^2/\psi^2$、$|\grad f|^2=4|\grad\psi|^2/\psi^2$ を代入すると、勾配二乗が消えて
$$
R^{(3)}=-8\psi^{-5}\Delta\psi.
$$
従って真空拘束は $\Delta\psi=0$ になる。例えば $\R^3\setminus\{0\}$ 上の $\psi=1+m/(2r)$ は正で調和的である。$\bar r=m^2/(4r)$ と置くと $r\to0$ の端でも計量が同じ形に戻る。従って原点を除いたことは有限距離に穴を開けたこととは異なり、二つ目の漸近端を持つデータである。

### 共形自由度を残した真空拘束

三次元で $\tilde h_{ij}=\psi^4h_{ij}$ とし、外的曲率を跡 $\tilde K$ と無跡の $\tilde A^{ij}$ に分ける。$h_{ij}A^{ij}=0$、$\tilde A^{ij}=\psi^{-10}A^{ij}$ とする。接続差を $C^i{}_{jk}=2(\delta^i_j\partial_k\log\psi+\delta^i_k\partial_j\log\psi-h_{jk}D^i\log\psi)$ と書けば、$\tilde D_j\tilde A^{ij}$ の微分から $-10\psi^{-10}A^{ij}\partial_j\log\psi$ が出る。二つの接続差はそれぞれ係数 $4$ と $6$ の同じ項を与え、跡項は零である。従って $\tilde D_j\tilde A^{ij}=\psi^{-10}D_j A^{ij}$ となる。$\tilde K$ が定数なら、真空の運動量拘束には $D_j A^{ij}=0$ を選べばよい。Hamilton 拘束は $\tilde A_{ij}\tilde A^{ij}=\psi^{-12}A_{ij}A^{ij}$ と $\tilde R=\psi^{-5}(R\psi-8\Delta_h\psi)$ を代入して
$$
-8\Delta_h\psi+R\psi-A_{ij}A^{ij}\psi^{-7}+
(2\tilde K^2/3-2\Lambda)\psi^5=0.
$$
従って計量の拘束の一部を、正の共形因子に対する非線形楕円方程式として扱える。例えば平坦で $A=0,\tilde K^2=3\Lambda$ なら定数 $\psi>0$ が解になる。任意の自由データで正の解があるとは限らないため、一般の場合は境界条件と可解性を調べる。

### 拘束の時間微分を Gaussian 座標で求める

$g=-\dd t^2+h_{ij}(t,x)\dd x^i\dd x^j$ とし、$E_{ab}=G_{ab}+\Lambda g_{ab}-8\pi G T_{ab}$ が対称かつ発散零とする。空間成分 $E_{ij}=0$ を発展方程式で課す。すると $C=E_{00},F_i=E_{0i}$ が残る。非零の混合時間接続は $\Gamma^0{}_{ij}=K_{ij},\Gamma^i{}_{0j}=K^i{}_j$ である。共変発散 $\nabla^a E_{ab}=0$ の $b=0$ 成分は
$$
-\dot C+h^{ij}(\partial_iF_j-\Gamma^k{}_{ij}F_k-K_{ij}C)=0 .
$$
従って $\dot C=D_iF^i-K C$ である。$b=k$ では時間側が $-\dot F_k+K^i{}_kF_i$、空間側が $-K F_k-K^i{}_kF_i$ だから、混合項は相殺して $\dot F_k=-K F_k$ である。よって初期に $F=0$ なら各座標点での線形常微分方程式からそのまま零、さらに初期に $C=0$ なら $C$ も零である。この計算は物質方程式が保存則を満たすことと、真に $E_{ij}=0$ を発展として課したことの両方を使う。

### 初期データの端を見分ける

:::example[共形因子を半径反転で検算する]
$h=(1+m/(2r))^4(\dd r^2+r^2\dd\Omega^2)$、$m>0$ とする。$R=m^2/(4r)$ では $1+m/(2r)=(2R/m)(1+m/(2R))$ であり、$\dd r^2+r^2\dd\Omega^2=(m/(2R))^4(\dd R^2+R^2\dd\Omega^2)$ である。掛けると $R$ についても同じ形になる。従って $r\to0$ は $R\to\infty$ の漸近端である。また $r=m/2$ は反転で固定され、ここで球の面積半径 $r(1+m/(2r))^2$ は $2m$ となる。座標の原点を外した以外に内部の境界条件を課さず、二つの端を持つ時間対称なデータが得られる。
:::

### 跡反転した発展式でも拘束が保存する

同じ Gaussian 座標で $C=E_{00},F_i=E_{0i}$ とする。接方向で $E_{ij}=h_{ij}(-C+h^{kl}E_{kl})/2$ を課して縮約すると $h^{ij}E_{ij}=3C$、従って $E_{ij}=h_{ij}C$ となる。前の発散計算で $E_{ij}$ の寄与を残すと
$$
\partial_t C=D_iF^i-2KC,\qquad
\partial_t F_i=D_i C-KF_i.
$$
第一式の追加項は $-K^{ij}E_{ij}=-KC$ であり、第二式では $D^jE_{ji}=D_iC$ が加わる。従って先の三角型の常微分方程式の議論は、この発展式には直接使えない。

例えば空間がコンパクトで境界がない領域では、$\mathcal E(t)=\tfrac12\int(C^2+h^{ij}F_iF_j)\sqrt h\,\dd^3x$ を使う。空間微分の寄与 $C D_iF^i+F^iD_iC=D_i(CF^i)$ の積分は零である。$\partial_t\sqrt h=K\sqrt h$、$\partial_t h^{ij}=-2K^{ij}$ による項と、発展式中の $K$ の項は、有限の滑らかな時間区間で $|\mathcal E'(t)|\le c\,\mathcal E(t)$ と評価できる。$\mathcal E(t_0)=0$ なら $\partial_t(e^{-ct}\mathcal E(t))\le0$ から未来にも零である。局所的には時間的な境界からの自由な流入を許すと結論は成り立たない。初期領域の依存領域で同じ評価を行う場合は、外向きの光的境界で $|C F^i s_i|\le(C^2+|F|^2)/2$ を使うと、未知の流入を要しない。

### 波動座標とゲージ拘束の伝播

座標関数の波動演算子 $\Box_g x^a=-g^{bc}\Gamma^a{}_{bc}$ を零とする条件を :term[調和座標条件]{en="harmonic coordinate condition" yomi="ちょうわざひょうじょうけん"} という。局所座標で $C^a=g^{bc}\Gamma^a{}_{bc}$ と置く。座標を固定した還元方程式を考える際は、この成分から作る場 $C$ をその座標近傍でベクトルとして扱う。Ricci の第二微分を集めると
$$
R_{ab}=-\tfrac12g^{cd}\partial_c\partial_dg_{ab}+\nabla_{(a}C_{b)}+Q_{ab}(g,\partial g)
$$
となる。実際、$\partial_c\Gamma^c{}_{ab}-\partial_b\Gamma^c{}_{ac}$ の第二微分部分は $g^{cd}(\partial_a\partial_cg_{bd}+\partial_b\partial_cg_{ad}-\partial_c\partial_dg_{ab}-\partial_a\partial_bg_{cd})/2$ で、第一、第二、第四項が $\partial_{(a}C_{b)}$ の主部である。従って $R_{ab}-\nabla_{(a}C_{b)}=\Lambda g_{ab}$ は真空の準線形波動系になる。

還元系の解では $E_{ab}=\nabla_{(a}C_{b)}-\tfrac12g_{ab}\nabla_cC^c$ である。Bianchi 恒等式でこの発散を取り、$\nabla^a\nabla_b C_a-\nabla_b\nabla^a C_a=R_{ba}C^a$ を使うと
$$
\Box_g C_b+R_{ba}C^a=0
$$
を得る。初期面で座標条件によって $C=0$ とすれば、その接方向微分は零である。正規直交枠で $E_{ni}=0$ は $\nabla_n C_i/2=0$、$E_{nn}=0$ は $\nabla_nC_n/2=0$ を与える。従って還元系の解は、拘束を満たす初期値から出れば $C$ の波動方程式の一意性によって元の Einstein 方程式を満たす。例えば平坦背景の一次変分 $v_{ab}$ では $C_b=\partial^a v_{ab}-\partial_bv/2$ であり、$\Box C_b=0$ となる。初期の $C_b$ と $\partial_t C_b$ がともに零でなければ、還元した計量方程式だけから調和条件を保証できない。

### 共形法で一定平均曲率のデータを組み立てる

平均曲率が一定なら、共形法では無跡で発散零の $A^{ij}$ を選び、残ったスカラー方程式を解ける。与えた無跡テンソル $B^{ij}$ から発散零のものを作る場合、ベクトル $W$ の変形
$$
(LW)^{ij}=D^iW^j+D^jW^i-\tfrac23 h^{ij}D_kW^k
$$
を引いて $A=B-LW$ と置く。$D_j(LW)^{ij}=D_jB^{ij}$ を解けばよい。交換子を縮約すると作用素は $\Delta_hW^i+\tfrac13D^iD_jW^j+R^i{}_jW^j$ となる。主シンボルは、非零の共変ベクトル $\zeta$ に直交する成分で $-|\zeta|^2$、平行な成分で $-4|\zeta|^2/3$ であり、楕円型である。コンパクトな境界なしの切片では部分積分して
$$
\int W_iD_j(LW)^{ij}\,\dd V_h
=-\tfrac12\int (LW)_{ij}(LW)^{ij}\,\dd V_h
$$
なので、核は共形 Killing 場である。源 $D_jB^{ij}$ はこの核と直交する。実際、共形 Killing 場 $Z$ に対して積分を移すと $-\int B^{ij}(LZ)_{ij}/2=0$ である。

:::example[トーラスで縦の成分を取り除く]
各辺 $2\pi$ の平坦なトーラスで $B^{ij}=\operatorname{diag}(2,-1,-1)\cos x$ とする。$W^x=(3/2)\sin x$、$W^y=W^z=0$ なら $D_iW^i=(3/2)\cos x$、従って $(LW)^{xx}=2\cos x$、$(LW)^{yy}=(LW)^{zz}=-\cos x$ である。ゆえにこの $B$ の全てが縦の成分であり、$A=0$ となる。対照として $B^{ij}=\operatorname{diag}(0,1,-1)\cos x$ なら、ただ一つの微分方向 $x$ の列が零なので初めから発散零である。後者を自由データとする場合も、Hamilton 拘束の可解性は別途確かめる必要がある。例えば $\Lambda=\tilde K=0$ では $-8\Delta\psi=A_{ij}A^{ij}\psi^{-7}$ を全トーラスで積分すると左辺が零で右辺が正になり、正の滑らかな解は存在しない。
:::

## 演習問題

::::exercise{level=basic}
$K_{ij}=Hh_{ij}$ の三次元面で $K$ と $K_{ij}K^{ij}$ を求めよ。

:::solution
$h_{ij}h^{ij}=3$ なので $K=3H$、$K_{ij}K^{ij}=3H^2$ である。
:::
::::

::::exercise{level=basic}
$\beta=0,N=2$ の法線観測者の $\Delta t=3$ に対する固有時を求めよ。

:::solution
$\dd\tau=N\dd t$ なので $\Delta\tau=6$ である。
:::
::::

::::exercise{level=basic}
時間対称真空面で $\Lambda=3/L^2$ の空間スカラー曲率を求めよ。

:::solution
:eqref[eq-gr-adm-constraints] に $K_{ij}=\rho=0$ を入れると $R^{(3)}=2\Lambda=6/L^2$ である。
:::
::::

::::exercise{level=basic}
Minkowski 直交座標の $t$ 一定面で $\pi^{ij}$ を求めよ。

:::solution
$\dot h=0,N=1,\beta=0$ なので $K_{ij}=0$ であり、$\pi^{ij}=0$ となる。
:::
::::

::::exercise{level=advanced}
$K_{ij}=H(t)h_{ij}$ の運動量拘束の左辺が零であることを示せ。

:::solution
$D_iH=0$、$Dh=0$ により $D_jK^j{}_i=0$ である。また $K=3H$ は空間で一定なので $D_iK=0$ となる。
:::
::::

::::exercise{level=advanced}
平坦な切片と $K_{ij}=Hh_{ij}$ の宇宙で $H^2$ を拘束から求めよ。

:::solution
$K^2-K_{ij}K^{ij}=6H^2$ なので $6H^2=16\pi G\rho+2\Lambda$、従って $H^2=8\pi G\rho/3+\Lambda/3$ である。
:::
::::

::::exercise{level=advanced}
$\pi^{ij}$ から $K^{ij}$ を求めよ。

:::solution
$\pi=-2K\sqrt h/(16\pi G)$ を元の定義へ代入し、$K^{ij}=(16\pi G/\sqrt h)(\pi^{ij}-h^{ij}\pi/2)$ となる。
:::
::::

::::exercise{level=advanced}
Palatini 変分で $n=2$ の場合に計量整合性が従わない理由を述べよ。

:::solution
$\nabla_cP^{ab}=0$ の縮約は $(1-n/2)g_{ab}\nabla_cg^{ab}=0$ である。$n=2$ では係数が零なので非計量性の跡を零と結論できない。
:::
::::

::::exercise{level=advanced}
空間計量 $h_{ij}=\delta_{ij}$、一定 $\Lambda=0$ で $K^i{}_j=\diag(H,H,J)$、$H\ne0$ を空間定数として指定する。真空の Hamilton 拘束と運動量拘束から $J$ を決め、その初期値を含む Kasner 枝を一つ求めよ。

:::solution
運動量拘束はすべて空間微分なので零である。Hamilton 拘束の外的項は $(2H+J)^2-(2H^2+J^2)=2H^2+4HJ$ である。従って $J=-H/2$ となる。$H>0$ なら $p=(2/3,2/3,-1/3)$、$t_0=2/(3H)$ として $a_i(t)=(t/t_0)^{p_i}$ を取る。切片 $t=t_0$ では $h_{ij}=\delta_{ij}$、$K^i{}_i=p_i/t_0$ が指定値になる。$\sum p=1,\sum p^2=1$ なので真空の発展である。$H<0$ の場合には正の $|H|$ に対するこの解の時間方向を反転すると外的曲率の符号が反転する。
:::
::::

## 確認テスト

::::quiz
Hamilton 拘束を課す変分を選べ。

- [x] $N$ の変分
- [ ] $\pi^{ij}$ の変分
- [ ] 空間方向の $h_{ij}$ の微分

:::explanation
正準作用に $-N\mathcal H_\perp$ があり、$N$ の時間微分はない。従って真空で $N$ を変分すると $\mathcal H_\perp=0$ である。
:::
::::

::::quiz{answer="6" tol="1e-9"}
$K_{ij}=h_{ij}$ の空間切片で $K^2-K_{ij}K^{ij}$ を求めよ。

:::explanation
$K=3$、$K_{ij}K^{ij}=3$ なので $3^2-3=6$ である。
:::
::::

::::quiz
GHY 項を加える目的を選べ。

- [x] 固定した誘導計量の法線微分に関する変分を相殺する。
- [ ] 物質の運動方程式を変える。
- [ ] 内部の Ricci を零と仮定する。

:::explanation
Einstein–Hilbert 作用の発散は境界上で $\delta K$ を含む。GHY 項が相殺し、内部の変分係数は変わらない。
:::
::::

::::quiz{answer="4" tol="1e-9"}
三次元の初期面で Einstein 方程式の拘束の成分数を答えよ。

:::explanation
法線二本の射影が一つ、法線一本と接方向の射影が三つで、計 $4$ である。
:::
::::

::::quiz{answer="1" tol="1e-7"}
Kasner 指数が $(-1/3,2/3,2/3)$ のとき $\sum_i p_i^2$ を求めよ。

:::explanation
二乗は順に $1/9,4/9,4/9$ で合計は $1$ である。$\sum_i p_i=1$ も同時に成り立つため、この三つは真空の Kasner 条件を満たす。
:::
::::

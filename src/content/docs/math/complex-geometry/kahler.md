---
title: Kähler 多様体
description: 閉じた基本形式、局所ポテンシャル、Lefschetz 作用素と Kähler 恒等式を結び付ける
notation:
  - symbol: '\omega'
    meaning: Hermite 計量の基本形式
  - symbol: 'L,\Lambda'
    meaning: ω との外積およびその随伴
  - symbol: '\rho'
    meaning: Ricci 形式
---

## この章の目標

- Hermite 計量が Kähler になる条件を局所座標で確かめる。
- Fubini–Study 形式を局所ポテンシャルで計算する。
- Lefschetz 分解を外積代数から導く。
- Kähler 恒等式を座標で証明し、曲率と Ricci 形式を求める。

## Kähler 計量と Kähler 形式

:::definition[Kähler 計量]{#def-cg-kahler-metric}
複素構造 $J$ を保つ Riemann 計量 $g$ の基本形式を $\omega(V,W)=g(JV,W)$ とする。$\dd\omega=0$ のとき $g$ を :term[Kähler 計量]{en="Kähler metric" yomi="けーらーけいりょう"}、$\omega$ を :term[Kähler 形式]{en="Kähler form" yomi="けーらーけいしき"} と呼ぶ。局所座標で
$$
g=2\Re\sum_{j,k}g_{j\bar k}\dd z_j\otimes\dd\bar z_k,\qquad
\omega=i\sum_{j,k}g_{j\bar k}\dd z_j\wedge\dd\bar z_k .
$$
:::

行列 $(g_{j\bar k})$ が Hermite 正定値であることは、$V\ne0$ で $\omega(V,JV)=g(V,V)>0$ と同値である。逆に正定値な実 $(1,1)$ 形式から $g(V,W)=\omega(V,JW)$ で元の計量を復元できる。Kähler 形式は型が $(1,1)$ なので、閉じる条件は $\partial\omega=\bar\partial\omega=0$ と同値である。係数の条件は
$$
\partial_\ell g_{j\bar k}=\partial_j g_{\ell\bar k}
$$
とその共役である。これは $\partial\omega$ の $\dd z_\ell\wedge\dd z_j\wedge\dd\bar z_k$ の係数を読めば従う。

:::example[曲線では自動的に閉じる]
複素次元 $1$ では $\omega=i h(z)\dd z\wedge\dd\bar z$、$h>0$ は実最高次形式なので $\dd\omega=0$ である。複素次元 $2$ では $\omega=i(1+|z_2|^2)\dd z_1\wedge\dd\bar z_1+i\dd z_2\wedge\dd\bar z_2$ は正定値だが、$\partial_2g_{1\bar1}=\bar z_2$、$\partial_1g_{2\bar1}=0$ となるため Kähler ではない。
:::

:::example[平坦トーラス]
$\C^n$ 上の定数係数 $\omega=i\sum_j h_j\dd z_j\wedge\dd\bar z_j$、$h_j>0$ は平行移動で不変なので、格子商に降りる。係数が定数のため閉じており、計量の接続係数もこの座標で零である。
:::

### 共形倍と Kähler の閉性

$\omega$ が Kähler 形式で $f$ が正の実関数なら、$\mathrm d(f\omega)=\mathrm df\wedge\omega$ である。複素次元 $n\ge2$ では、これが零なら $f$ は各連結成分で定数となる。実際、$\omega^{n-2}$ を掛ければ $\mathrm df\wedge\omega^{n-1}=0$ である。この写像は実 1 形式上で単射である。$\omega$ を $\sum_j e_j\wedge f_j$ とする実基底では、1 形式の $e_j$ 成分を掛けて残る項が $\{e_k,f_k:k\ne j\}$ の積と $e_j$ であり、$f_j$ を欠く最高次の一歩手前の基底になる。他の各成分も異なる基底を生むので単射である。

:::example[同じ双曲型の座標でも次元で変わる]
$\mathbb H$ 上では $y^{-2}\mathrm dx\wedge\mathrm dy$ は最高次形式なので閉じる。これを積 $\mathbb H\times\C$ の全ての方向に共形倍として用いると、$\mathrm df=-2y^{-3}\mathrm dy$ が第二因子の面積形式との積に残り、閉じない。第一因子だけを $y^{-2}$ 倍した積計量なら両因子の面積形式の和なので閉じる。
:::

### 局所の閉性と正値性

:::example[ポテンシャルから体積と Ricci を追う]
$\varphi=\bar z_1 z_{1} \left(\bar z_1 + \bar z_2\right) \left(z_{1} + z_{2}\right) + \bar z_1 z_{1} + \bar z_2 z_{2}$ とし、点 $p=(0,1)$ で複素 Hessian と Ricci 形式の係数を求める。$\partial_{z_j}$ は $\bar z_k$ を固定して微分するため、計算中は両者を独立な変数として扱う。

$$
\begin{aligned}
G=(\varphi_{j\bar k})&=\left[\begin{matrix}\bar z_1 z_{1} + \bar z_1 \left(z_{1} + z_{2}\right) + z_{1} \left(\bar z_1 + \bar z_2\right) + \left(\bar z_1 + \bar z_2\right) \left(z_{1} + z_{2}\right) + 1 & \bar z_1 \left(2 z_{1} + z_{2}\right)\\z_{1} \left(2 \bar z_1 + \bar z_2\right) & \bar z_1 z_{1} + 1\end{matrix}\right] \\
\det G&=5 \bar z_1 z_{1} + 2 \bar z_1 z_{2} + 2 \bar z_2 z_{1} + \bar z_2 z_{2} + 1 \\
G(p)&=\left[\begin{matrix}2 & 0\\0 & 1\end{matrix}\right],\qquad\det G(p)=2
\end{aligned}
$$

点での先頭小行列式 $G_{1\bar1}(p)=2$ と行列式がともに正なので近傍で Kähler 計量になる。体積密度は $4\det G$ である。Ricci 形式の係数 $R$ は $-\partial_j\partial_{\bar k}\log\det G$ だから、商の微分で次の式を得る。

$$
\begin{aligned}
R_{j\bar k}&=-\frac{d_{j\bar k}}d+\frac{d_jd_{\bar k}}{d^2}\qquad(d=\det G) \\
R(p)&=\left[\begin{matrix}- \frac{3}{2} & - \frac{1}{2}\\- \frac{1}{2} & - \frac{1}{4}\end{matrix}\right] \\
\Lambda_\omega\rho(p)&=\tr(G(p)^{-1}R(p))=-1
\end{aligned}
$$

対数行列式の二階微分は、点での行列式を先に数へ置き換えると失われる。そこで行列式を関数として微分してから、その点の値を代入した。正の複素 Hessian を取ったので、閉性は二階偏微分の交換から従う。
:::

:::example[ポテンシャルから体積と Ricci を追う]
$\varphi=\bar z_1 \bar z_2 z_{1} z_{2} + 2 \bar z_1 z_{1} + \bar z_2 z_{2}$ とし、点 $p=(1,2)$ で複素 Hessian と Ricci 形式の係数を求める。$\partial_{z_j}$ は $\bar z_k$ を固定して微分するため、計算中は両者を独立な変数として扱う。

$$
\begin{aligned}
G=(\varphi_{j\bar k})&=\left[\begin{matrix}\bar z_2 z_{2} + 2 & \bar z_1 z_{2}\\\bar z_2 z_{1} & \bar z_1 z_{1} + 1\end{matrix}\right] \\
\det G&=2 \bar z_1 z_{1} + \bar z_2 z_{2} + 2 \\
G(p)&=\left[\begin{matrix}6 & 2\\2 & 2\end{matrix}\right],\qquad\det G(p)=8
\end{aligned}
$$

点での先頭小行列式 $G_{1\bar1}(p)=6$ と行列式がともに正なので近傍で Kähler 計量になる。体積密度は $4\det G$ である。Ricci 形式の係数 $R$ は $-\partial_j\partial_{\bar k}\log\det G$ だから、商の微分で次の式を得る。

$$
\begin{aligned}
R_{j\bar k}&=-\frac{d_{j\bar k}}d+\frac{d_jd_{\bar k}}{d^2}\qquad(d=\det G) \\
R(p)&=\left[\begin{matrix}- \frac{3}{16} & \frac{1}{16}\\\frac{1}{16} & - \frac{1}{16}\end{matrix}\right] \\
\Lambda_\omega\rho(p)&=\tr(G(p)^{-1}R(p))=- \frac{1}{8}
\end{aligned}
$$

対数行列式の二階微分は、点での行列式を先に数へ置き換えると失われる。そこで行列式を関数として微分してから、その点の値を代入した。混合成分を持つので、正値性を対角成分だけで判断しない。
:::

:::example[対角形式の外積作用素を計算する]
複素次元 $2$ の一点で正則余接枠 $\eta_j$ を取り、$a_j=i\eta_j\wedge\bar\eta_j$ とする。$a_j^2=0$ である。最高次形式 $\beta$ の $a_1\wedge\cdots\wedge a_n$ に関する係数を $[\beta]_{\mathrm{top}}$ と書く。基本形式には次の正係数の線形結合を使う。

$$
\begin{aligned}
\omega&=2 a_{1} + 3 a_{2}
\end{aligned}
$$

外積作用素 $L=\omega\wedge$ の $n-2$ 乗は、$a_j$ の張る実空間を $n-1$ 個の添字の積が張る空間に写す。上位側の基底を昇順の添字の積とし、$a_i$ に $\omega^{n-2}$ を掛けると、同じ添字に二度出会う項が消える。残りの各添字集合の係数には並べ替えの数 $(n-2)!$ が掛かる。

$$
\begin{aligned}
[L^{0}]&=\left[\begin{matrix}1 & 0\\0 & 1\end{matrix}\right] \\
\det[L^{0}]&=1
\end{aligned}
$$

最高次の係数 $Q(a,b)=[a\wedge b\wedge\omega^{n-2}]_{\mathrm{top}}$ は、行と列が異なるときに残りの添字の重みの積から求まる。$Q(a,\omega)=0$ となる空間を $P$ とし、基底を列に並べた行列を $N$ とする。

$$
\begin{aligned}
[Q]&=\left[\begin{matrix}0 & 1\\1 & 0\end{matrix}\right],\qquad Q(\omega,\omega)=12 \\
N&=\left[\begin{matrix}- \frac{2}{3}\\1\end{matrix}\right] \\
Q|_P&=\left[\begin{matrix}- \frac{4}{3}\end{matrix}\right]
\end{aligned}
$$

$Q(a,\omega)=0$ は $a\wedge\omega^{n-1}=0$ と同値である。実際、$\omega=\sum_j w_j a_j$ とすると $w_j a_j$ がユニタリ枠に対応する基本形式で、$\Lambda a_j=1/w_j$ である。一方 $Q(a_j,\omega)=(n-1)!\prod_{k\ne j}w_k$ であり、この値は $1/w_j$ の一定倍である。従って $P$ は原始成分に一致する。行列 $Q|_P=N^\top QN$ は、原始成分の外積の符号を与える。
:::

## 例：Fubini–Study 計量

$\mathbb P^n$ の $U_i$ に $\varphi_i=\log(\sum_j|Z_j/Z_i|^2)$ を置く。重なりでは $\varphi_j-\varphi_i=\log|Z_i/Z_j|^2$ である。正則な非消失関数の対数絶対値二乗は局所正則関数の実部なので、$\partial\bar\partial$ が零である。従って
$$
\omega_{\mathrm{FS}}=i\partial\bar\partial\varphi_i
$$
は大域形式になる。$\partial^2=\bar\partial^2=0$ で閉じており、:ref[def-cg-positive-line] に続く Hessian の計算から正定値である。これを :term[Fubini–Study 計量]{en="Fubini–Study metric" yomi="ふびにすたでぃけいりょう"} と呼ぶ。ここでは $[\omega_{\mathrm{FS}}/(2\pi)]=c_1(\mathcal O(1))$ となる規格化を使う。

:::example[原点での成分と面積]
$U_0$ の原点で $g_{j\bar k}=\delta_{jk}$ である。$\mathbb P^1$ に制限すると $\omega_{\mathrm{FS}}=2r\dd r\wedge\dd\theta/(1+r^2)^2$ なので全体の面積は $2\pi$ である。これは第一 Chern 類の積分 $1$ に $2\pi$ を掛けた値である。
:::

:::example[部分多様体に制限する]
滑らかな射影部分多様体 $Y\hookrightarrow\mathbb P^n$ へ引き戻した形式は閉じており、接空間への制限も正定値である。例えば二次 Veronese 写像で $\mathbb P^1$ に引き戻すとコホモロジー類は $2[\omega_{\mathrm{FS}}]$ なので面積は $4\pi$ である。これは引き戻し束が $\mathcal O(2)$ になることから従う。
:::

### ユニタリ回転と球面積分

$S=1+\sum_j|z_j|^2$ とすれば
$$
g_{j\bar k}=S^{-1}\delta_{jk}-S^{-2}\bar z_jz_k.
$$
計量二次形式では第二項が $|\sum_j\bar z_jv_j|^2/S^2$ を引く。Cauchy–Schwarz から $g(v,\bar v)\ge\|v\|^2/S^2$ である。$v$ が $z$ に比例すると等号が成立し、$z$ に直交すると第二項は零になる。従って半径方向とそれ以外の固有値は別々である。このことから $\det g=S^{-n-1}$ が得られる。$\mathbb P^n$ の体積を $\omega_{\rm FS}^n/n!$ で測ると、コホモロジーの積分より $(2\pi)^n/n!$ である。アフィンチャートで $\det g$ を積分しても一致する。

:::example[射影平面を極座標で積分する]
$\C^2$ において $\omega_{\rm FS}^2/2=4(1+r^2)^{-3}\mathrm dV_{\R^4}$ であり、単位 $3$ 球面の体積は $2\pi^2$ である。従って体積は
$$
8\pi^2\int_0^\infty\frac{r^3}{(1+r^2)^3}\mathrm dr
=8\pi^2\cdot\frac14=2\pi^2.
$$
無限遠の超平面は実余次元二なので積分の値を変えない。
:::

### 射影積の体積を因子ごとに取る

:::example[対角形式の外積作用素を計算する]
複素次元 $3$ の一点で正則余接枠 $\eta_j$ を取り、$a_j=i\eta_j\wedge\bar\eta_j$ とする。$a_j^2=0$ である。最高次形式 $\beta$ の $a_1\wedge\cdots\wedge a_n$ に関する係数を $[\beta]_{\mathrm{top}}$ と書く。基本形式には次の正係数の線形結合を使う。

$$
\begin{aligned}
\omega&=a_{1} + 2 a_{2} + 3 a_{3}
\end{aligned}
$$

外積作用素 $L=\omega\wedge$ の $n-2$ 乗は、$a_j$ の張る実空間を $n-1$ 個の添字の積が張る空間に写す。上位側の基底を昇順の添字の積とし、$a_i$ に $\omega^{n-2}$ を掛けると、同じ添字に二度出会う項が消える。残りの各添字集合の係数には並べ替えの数 $(n-2)!$ が掛かる。

$$
\begin{aligned}
[L^{1}]&=\left[\begin{matrix}2 & 1 & 0\\3 & 0 & 1\\0 & 3 & 2\end{matrix}\right] \\
\det[L^{1}]&=-12
\end{aligned}
$$

最高次の係数 $Q(a,b)=[a\wedge b\wedge\omega^{n-2}]_{\mathrm{top}}$ は、行と列が異なるときに残りの添字の重みの積から求まる。$Q(a,\omega)=0$ となる空間を $P$ とし、基底を列に並べた行列を $N$ とする。

$$
\begin{aligned}
[Q]&=\left[\begin{matrix}0 & 3 & 2\\3 & 0 & 1\\2 & 1 & 0\end{matrix}\right],\qquad Q(\omega,\omega)=36 \\
N&=\left[\begin{matrix}- \frac{1}{2} & - \frac{1}{3}\\1 & 0\\0 & 1\end{matrix}\right] \\
Q|_P&=\left[\begin{matrix}-3 & -1\\-1 & - \frac{4}{3}\end{matrix}\right]
\end{aligned}
$$

$Q(a,\omega)=0$ は $a\wedge\omega^{n-1}=0$ と同値である。実際、$\omega=\sum_j w_j a_j$ とすると $w_j a_j$ がユニタリ枠に対応する基本形式で、$\Lambda a_j=1/w_j$ である。一方 $Q(a_j,\omega)=(n-1)!\prod_{k\ne j}w_k$ であり、この値は $1/w_j$ の一定倍である。従って $P$ は原始成分に一致する。行列 $Q|_P=N^\top QN$ は、原始成分の外積の符号を与える。
:::

:::example[対角形式の外積作用素を計算する]
複素次元 $4$ の一点で正則余接枠 $\eta_j$ を取り、$a_j=i\eta_j\wedge\bar\eta_j$ とする。$a_j^2=0$ である。最高次形式 $\beta$ の $a_1\wedge\cdots\wedge a_n$ に関する係数を $[\beta]_{\mathrm{top}}$ と書く。基本形式には次の正係数の線形結合を使う。

$$
\begin{aligned}
\omega&=2 a_{1} + 3 a_{2} + a_{3} + 4 a_{4}
\end{aligned}
$$

外積作用素 $L=\omega\wedge$ の $n-2$ 乗は、$a_j$ の張る実空間を $n-1$ 個の添字の積が張る空間に写す。上位側の基底を昇順の添字の積とし、$a_i$ に $\omega^{n-2}$ を掛けると、同じ添字に二度出会う項が消える。残りの各添字集合の係数には並べ替えの数 $(n-2)!$ が掛かる。

$$
\begin{aligned}
[L^{2}]&=\left[\begin{matrix}6 & 4 & 12 & 0\\24 & 16 & 0 & 12\\8 & 0 & 16 & 4\\0 & 8 & 24 & 6\end{matrix}\right] \\
\det[L^{2}]&=-27648
\end{aligned}
$$

最高次の係数 $Q(a,b)=[a\wedge b\wedge\omega^{n-2}]_{\mathrm{top}}$ は、行と列が異なるときに残りの添字の重みの積から求まる。$Q(a,\omega)=0$ となる空間を $P$ とし、基底を列に並べた行列を $N$ とする。

$$
\begin{aligned}
[Q]&=\left[\begin{matrix}0 & 8 & 24 & 6\\8 & 0 & 16 & 4\\24 & 16 & 0 & 12\\6 & 4 & 12 & 0\end{matrix}\right],\qquad Q(\omega,\omega)=576 \\
N&=\left[\begin{matrix}- \frac{2}{3} & -2 & - \frac{1}{2}\\1 & 0 & 0\\0 & 1 & 0\\0 & 0 & 1\end{matrix}\right] \\
Q|_P&=\left[\begin{matrix}- \frac{32}{3} & -16 & -4\\-16 & -96 & -12\\-4 & -12 & -6\end{matrix}\right]
\end{aligned}
$$

$Q(a,\omega)=0$ は $a\wedge\omega^{n-1}=0$ と同値である。実際、$\omega=\sum_j w_j a_j$ とすると $w_j a_j$ がユニタリ枠に対応する基本形式で、$\Lambda a_j=1/w_j$ である。一方 $Q(a_j,\omega)=(n-1)!\prod_{k\ne j}w_k$ であり、この値は $1/w_j$ の一定倍である。従って $P$ は原始成分に一致する。行列 $Q|_P=N^\top QN$ は、原始成分の外積の符号を与える。
:::

## Kähler ポテンシャルと局所的特徴付け

:::proposition[局所 Kähler ポテンシャル]{#prop-cg-kahler-potential}
実 $(1,1)$ 形式が閉じているなら、各点の近くで $\omega=i\partial\bar\partial\varphi$ と書ける。$\varphi$ は滑らかな実関数であり、二つの選択の差は局所的に正則関数の実部である。
:::

:::proof
局所 Poincaré の補題で実 $1$ 形式 $\eta$ を $\omega=\dd\eta$ となるように選ぶ。型 $(0,2)$ が零なので $\bar\partial\eta^{0,1}=0$ であり、:ref[thm-cg-local-dolbeault] から $\eta^{0,1}=\bar\partial f$ と書ける。$\eta$ の実性から $\eta^{1,0}=\partial\bar f$ である。従って $\dd\eta=\partial\bar\partial(f-\bar f)=i\partial\bar\partial(2\Im f)$ となる。差 $u$ で $\partial\bar\partial u=0$ なら $2\partial u$ は閉正則 $1$ 形式で、局所原始関数 $F$ を持つ。$u$ の実性により $\dd(\Re F)=\dd u$ となり、定数を調整すれば $u=\Re F$ である。
:::

定数を加える自由度だけでなく、例えば $\varphi$ と $\varphi+\Re(z_1^2+3z_2)$ は同じ形式を与える。一方、コンパクト正次元 Kähler 多様体には、大域ポテンシャルは存在しない。あれば $\omega=\dd(i\bar\partial\varphi)$ であり、$\omega^n$ を完全形式として積分すると零になるが、正の体積を与えるため矛盾する。

:::proposition[Kähler 正規座標]{#prop-cg-kahler-normal}
任意の点 $p$ で、$g_{j\bar k}(p)=\delta_{jk}$ かつ全ての一次導関数が $p$ で零になる正則座標を取れる。Kähler 計量の Levi-Civita 接続は $J$ を保ち、正則接束上で Chern 接続に一致する。
:::

:::proof
まず線形変更で $g(p)$ を単位行列にする。$p=0$ とし $C^k_{ij}=\partial_i g_{j\bar k}(0)$ と置くと $C^k_{ij}=C^k_{ji}$ である。新座標 $w_k=z_k+\frac12\sum_{ij}C^k_{ij}z_i z_j$ を選ぶと、その逆の二次係数は符号が逆になる。計量の変換式を一回微分すると、逆変換の二次係数から $-C^k_{ij}$ が加わり、元の $\partial_i g_{j\bar k}(0)$ を消す。共役の導関数も零となる。

実座標での計量の一次導関数が零なので、Levi-Civita 接続の係数もその点で零となる。正則座標で $J$ は定数行列であり、従って $\nabla J$ はその点で零である。各点で同じことが言える。正則接束の Chern 接続の行列もその点で零なので、両接続の差は各点で零となり一致する。逆に Levi-Civita 接続が $J$ を保てば $\nabla\omega=0$ であり、Levi-Civita の捩率がないので $\dd\omega$ は $\nabla\omega$ の交代化となって消える。
:::

### 実計量と複素 Hessian の規格化を揃える

平坦な $g_{j\bar k}=\delta_{jk}$ の座標では $g=2\sum_j(\dd x_j^2+\dd y_j^2)$、$\omega=2\sum_j\dd x_j\wedge\dd y_j$ である。従って $\omega^n/n!$ は $2^n\prod_j(\dd x_j\wedge\dd y_j)$ となり、実計量の体積と一致する。一般の点でも係数をユニタリ枠にすれば同じ式になる。この規格化では別の文献で $\omega=(i/2)\sum_j\dd z_j\wedge\dd\bar z_j$ とする平坦計量と、実計量の倍率が二倍異なる。後の Ricci 形式や曲率の数値を比較するときにこの倍率を確認する。

:::example[原点で同じ計量でも曲率は異なる]
$\varphi_1=|z|^2$ と $\varphi_2=\log(1+|z|^2)$ はともに原点で $g_{z\bar z}=1$ となる。しかし後者では $g_{z\bar z}=(1+|z|^2)^{-2}$ で、原点の近くの二次係数が $-2$ になる。前者は係数が定数なので、原点での計量が等しくても曲率は一致しない。
:::

### ポテンシャルを合成したときの平均値

実関数 $\varphi$ に対して、正則写像 $F$ は $\partial,\bar\partial$ と可換なので、引き戻した形式の局所ポテンシャルは $\varphi\circ F$ である。特に $\varphi=\sum |w_j|^2$ なら、複素 Hessian は Jacobian の列の内積になり、階数は Jacobian の階数に等しい。ポテンシャルが実の凸関数である必要はなく、複素 Hessian の正定値性を判定する。

:::example[グラフに誘導されるポテンシャル]
$F(z)=(z,2z^2)$ に平坦なポテンシャルを引き戻すと $\varphi=|z|^2+4|z|^4$ である。複素 Hessian は $1+16|z|^2$ であり、$|F'_1|^2+|F'_2|^2$ と一致する。第一成分を除いた $F(z)=2z^2$ だけを使うと原点で Hessian が零になり、近くで非負でも Kähler 計量にはならない。
:::

### ポテンシャルの選び方と微分

:::example[ポテンシャルから体積と Ricci を追う]
$\varphi=4 \bar z_1 z_{1} + 2 \bar z_2 z_{2} + \left(\bar z_1 - 2 \bar z_2\right)^{2} \left(z_{1} - 2 z_{2}\right)^{2}$ とし、点 $p=(0,1)$ で複素 Hessian と Ricci 形式の係数を求める。$\partial_{z_j}$ は $\bar z_k$ を固定して微分するため、計算中は両者を独立な変数として扱う。

$$
\begin{aligned}
G=(\varphi_{j\bar k})&=\left[\begin{matrix}4 \left(\left(\bar z_1 - 2 \bar z_2\right) \left(z_{1} - 2 z_{2}\right) + 1\right) & - 8 \left(\bar z_1 - 2 \bar z_2\right) \left(z_{1} - 2 z_{2}\right)\\- 8 \left(\bar z_1 - 2 \bar z_2\right) \left(z_{1} - 2 z_{2}\right) & 2 \left(8 \left(\bar z_1 - 2 \bar z_2\right) \left(z_{1} - 2 z_{2}\right) + 1\right)\end{matrix}\right] \\
\det G&=8 \left(9 \bar z_1 z_{1} - 18 \bar z_1 z_{2} - 18 \bar z_2 z_{1} + 36 \bar z_2 z_{2} + 1\right) \\
G(p)&=\left[\begin{matrix}20 & -32\\-32 & 66\end{matrix}\right],\qquad\det G(p)=296
\end{aligned}
$$

点での先頭小行列式 $G_{1\bar1}(p)=20$ と行列式がともに正なので近傍で Kähler 計量になる。体積密度は $4\det G$ である。Ricci 形式の係数 $R$ は $-\partial_j\partial_{\bar k}\log\det G$ だから、商の微分で次の式を得る。

$$
\begin{aligned}
R_{j\bar k}&=-\frac{d_{j\bar k}}d+\frac{d_jd_{\bar k}}{d^2}\qquad(d=\det G) \\
R(p)&=\left[\begin{matrix}- \frac{9}{1369} & \frac{18}{1369}\\\frac{18}{1369} & - \frac{36}{1369}\end{matrix}\right] \\
\Lambda_\omega\rho(p)&=\tr(G(p)^{-1}R(p))=- \frac{81}{202612}
\end{aligned}
$$

対数行列式の二階微分は、点での行列式を先に数へ置き換えると失われる。そこで行列式を関数として微分してから、その点の値を代入した。ポテンシャルに任意の正則多項式の実部を足しても、表の Hessian は変わらない。
:::

:::example[ポテンシャルから体積と Ricci を追う]
$\varphi=\bar z_1^{2} z_{1}^{2} + \bar z_1 z_{1} + \bar z_2^{2} z_{2}^{2} + 2 \bar z_2 z_{2}$ とし、点 $p=(1,1)$ で複素 Hessian と Ricci 形式の係数を求める。$\partial_{z_j}$ は $\bar z_k$ を固定して微分するため、計算中は両者を独立な変数として扱う。

$$
\begin{aligned}
G=(\varphi_{j\bar k})&=\left[\begin{matrix}4 \bar z_1 z_{1} + 1 & 0\\0 & 2 \left(2 \bar z_2 z_{2} + 1\right)\end{matrix}\right] \\
\det G&=2 \left(4 \bar z_1 z_{1} + 1\right) \left(2 \bar z_2 z_{2} + 1\right) \\
G(p)&=\left[\begin{matrix}5 & 0\\0 & 6\end{matrix}\right],\qquad\det G(p)=30
\end{aligned}
$$

点での先頭小行列式 $G_{1\bar1}(p)=5$ と行列式がともに正なので近傍で Kähler 計量になる。体積密度は $4\det G$ である。Ricci 形式の係数 $R$ は $-\partial_j\partial_{\bar k}\log\det G$ だから、商の微分で次の式を得る。

$$
\begin{aligned}
R_{j\bar k}&=-\frac{d_{j\bar k}}d+\frac{d_jd_{\bar k}}{d^2}\qquad(d=\det G) \\
R(p)&=\left[\begin{matrix}- \frac{4}{25} & 0\\0 & - \frac{2}{9}\end{matrix}\right] \\
\Lambda_\omega\rho(p)&=\tr(G(p)^{-1}R(p))=- \frac{233}{3375}
\end{aligned}
$$

対数行列式の二階微分は、点での行列式を先に数へ置き換えると失われる。そこで行列式を関数として微分してから、その点の値を代入した。このポテンシャルは因子の和なので混合係数が零になる。
:::

:::example[平坦トーラスの体積指定を積の形で解く]
複素次元 $2$ の正方トーラスで、座標 $z_j=x_j+iy_j$ の両実周期を $1$ とする。初期の計量行列を $I$ とし、体積比として次の積を指定する。

$$
\begin{aligned}
e^f&=\left(1 - \frac{\cos{\left(4 \pi x_{2} \right)}}{3}\right) \left(\frac{\cos{\left(2 \pi x_{1} \right)}}{2} + 1\right)
\end{aligned}
$$

各因子は正で、周期積分では余弦の平均が零なので因子ごとの平均は $1$ である。従って積の全平均も $1$ となる。ポテンシャルは各一変数の関数を足した形に取る。混合偏微分が消えるため、行列式が指定の積になるには対角成分だけを調整すればよい。

$$
\begin{aligned}
\varphi&=- \frac{\cos{\left(2 \pi x_{1} \right)}}{2 \pi^{2}} + \frac{\cos{\left(4 \pi x_{2} \right)}}{12 \pi^{2}} \\
g_\varphi&=\left[\begin{matrix}\frac{\cos{\left(2 \pi x_{1} \right)}}{2} + 1 & 0\\0 & 1 - \frac{\cos{\left(4 \pi x_{2} \right)}}{3}\end{matrix}\right] \\
\det(g_\varphi)&=\left(1 - \frac{\cos{\left(4 \pi x_{2} \right)}}{3}\right) \left(\frac{\cos{\left(2 \pi x_{1} \right)}}{2} + 1\right)
\end{aligned}
$$

各方向について $u_{z_j\bar z_j}=u_{x_jx_j}/4$ を計算すれば対角成分を再現する。従って新しい複素 Hessian は正定値で、同じ Kähler 類にある。さらに対数行列式の二階微分から、新しい Ricci 形式の係数は次の対角行列になる。

$$
\begin{aligned}
R_{1\overline{1}}&=\frac{\pi^{2} \left(\sin^{2}{\left(2 \pi x_{1} \right)} + \cos^{2}{\left(2 \pi x_{1} \right)} + 2 \cos{\left(2 \pi x_{1} \right)}\right)}{\left(\cos{\left(2 \pi x_{1} \right)} + 2\right)^{2}} \\
R_{2\overline{2}}&=\frac{4 \pi^{2} \left(\sin^{2}{\left(4 \pi x_{2} \right)} + \cos^{2}{\left(4 \pi x_{2} \right)} - 3 \cos{\left(4 \pi x_{2} \right)}\right)}{\left(\cos{\left(4 \pi x_{2} \right)} - 3\right)^{2}}
\end{aligned}
$$

新しい体積が旧体積と同じ平均を持つことと、Ricci が点ごとに零であることは異なる。各 $\lambda_j$ の正値性から計量は定義されるが、一般の点で Ricci は零にならない。大域に定義されるのは初期形式との差のポテンシャルである。初期形式そのものの大域ポテンシャルを仮定していない。
:::

:::example[平坦トーラスの体積指定を積の形で解く]
複素次元 $2$ の正方トーラスで、座標 $z_j=x_j+iy_j$ の両実周期を $1$ とする。初期の計量行列を $I$ とし、体積比として次の積を指定する。

$$
\begin{aligned}
e^f&=\left(1 - \frac{\cos{\left(4 \pi x_{1} \right)}}{2}\right) \left(\frac{\cos{\left(2 \pi x_{2} \right)}}{4} + 1\right)
\end{aligned}
$$

各因子は正で、周期積分では余弦の平均が零なので因子ごとの平均は $1$ である。従って積の全平均も $1$ となる。ポテンシャルは各一変数の関数を足した形に取る。混合偏微分が消えるため、行列式が指定の積になるには対角成分だけを調整すればよい。

$$
\begin{aligned}
\varphi&=\frac{\cos{\left(4 \pi x_{1} \right)}}{8 \pi^{2}} - \frac{\cos{\left(2 \pi x_{2} \right)}}{4 \pi^{2}} \\
g_\varphi&=\left[\begin{matrix}1 - \frac{\cos{\left(4 \pi x_{1} \right)}}{2} & 0\\0 & \frac{\cos{\left(2 \pi x_{2} \right)}}{4} + 1\end{matrix}\right] \\
\det(g_\varphi)&=\left(1 - \frac{\cos{\left(4 \pi x_{1} \right)}}{2}\right) \left(\frac{\cos{\left(2 \pi x_{2} \right)}}{4} + 1\right)
\end{aligned}
$$

各方向について $u_{z_j\bar z_j}=u_{x_jx_j}/4$ を計算すれば対角成分を再現する。従って新しい複素 Hessian は正定値で、同じ Kähler 類にある。さらに対数行列式の二階微分から、新しい Ricci 形式の係数は次の対角行列になる。

$$
\begin{aligned}
R_{1\overline{1}}&=\frac{4 \pi^{2} \left(\sin^{2}{\left(4 \pi x_{1} \right)} + \cos^{2}{\left(4 \pi x_{1} \right)} - 2 \cos{\left(4 \pi x_{1} \right)}\right)}{\left(\cos{\left(4 \pi x_{1} \right)} - 2\right)^{2}} \\
R_{2\overline{2}}&=\frac{\pi^{2} \left(\sin^{2}{\left(2 \pi x_{2} \right)} + \cos^{2}{\left(2 \pi x_{2} \right)} + 4 \cos{\left(2 \pi x_{2} \right)}\right)}{\left(\cos{\left(2 \pi x_{2} \right)} + 4\right)^{2}}
\end{aligned}
$$

新しい体積が旧体積と同じ平均を持つことと、Ricci が点ごとに零であることは異なる。各 $\lambda_j$ の正値性から計量は定義されるが、一般の点で Ricci は零にならない。各座標に周期を入れたため、局所座標の選び直しでも差のポテンシャルはそのまま貼り合う。
:::

### 正則部分多様体へ正規座標を制限する

正則はめ込み $F:Y\to X$ に対し、Kähler 形式の引き戻しは閉である。非零接ベクトル $v$ の像が非零であり $\mathrm dF(Jv)=J\mathrm dF(v)$ なので、$F^*\omega(v,Jv)>0$ である。従って誘導計量も Kähler である。ただし環境の正規座標の一部を $Y$ 上に制限したものが既に $Y$ の正規座標になるとは限らない。計量成分が一定でも、$F$ の Jacobian が変動し、誘導計量の導関数を生むからである。

:::example[はめ込みの二次項が計量に現れる]
$F(z)=(z,z+z^2)$ を平坦な $\C^2$ に入れると、誘導計量の複素成分は $1+|1+2z|^2=2+2(z+\bar z)+4|z|^2$ である。$z=0$ の $z$ 微分は $2$ であり、環境の計量が定数でも零ではない。まず $u=\sqrt2z$ として中心の計量を $1$ にし、次に $w=u+u^2/(2\sqrt2)$ とおけば、逆変換の二次係数が一次の微分を打ち消す。
:::

## Hermite 線形代数と Lefschetz 作用素

一点でユニタリな $(1,0)$ 余接枠 $\theta_j$ を選び、$\omega=i\sum_j\theta_j\wedge\bar\theta_j$ とする。$e_j,\bar e_j$ をそれぞれ $\theta_j,\bar\theta_j$ の外積、$\iota_j,\bar\iota_j$ をその Hermite 随伴とする。外積の基底に作用させると $\iota_j e_k+e_k\iota_j=\delta_{jk}$、$\bar\iota_j\bar e_k+\bar e_k\bar\iota_j=\delta_{jk}$、混合する外積と収縮の反交換子は零となる。添字を挿入、除去する二順序を比較すれば、異なる添字では符号が逆で、同じ添字では片方の順序だけが恒等となる。

:::definition[Lefschetz 作用素と原始形式]{#def-cg-lefschetz-operators}
$L=\omega\wedge=i\sum_j e_j\bar e_j$ を :term[Lefschetz 作用素]{en="Lefschetz operator" yomi="れふしぇっつさようそ"}、その随伴を $\Lambda=-i\sum_j\bar\iota_j\iota_j$ とする。$\Lambda\alpha=0$ となる形式 $\alpha$ を :term[原始形式]{en="primitive form" yomi="げんしけいしき"} という。
:::

:::proposition[Lefschetz 分解の線形代数]{#prop-cg-lefschetz-linear}
全次数 $k$ では $[\Lambda,L]=(n-k)\id$ である。次数 $k\le n$ の原始形式 $\alpha$ について
$$
\Lambda L^j\alpha=j(n-k-j+1)L^{j-1}\alpha.
\label{eq-cg-lefschetz-string}
$$
外積代数全体は、原始形式から始まる列 $\alpha,L\alpha,\ldots,L^{n-k}\alpha$ の直交直和に分解する。また $0\le k\le n$ で $L^{n-k}:\wedge^k\to\wedge^{2n-k}$ は同型である。
:::

:::proof
各添字の対 $(\theta_j,\bar\theta_j)$ を含む数が $0,1,2$ のとき、その対の交換子は基底上で $1,0,-1$ となる。添字ごとに足すと $n-k$ である。原始 $\alpha$ に対し $\Lambda L^j$ に交換子を順に挿入すると $\sum_{a=0}^{j-1}(n-k-2a)$ が現れ、:eqref[eq-cg-lefschetz-string] になる。内積で $\|L^j\alpha\|^2=j(n-k-j+1)\|L^{j-1}\alpha\|^2$ を得る。よって非零原始形式では $k\le n$ であり、列は $j=n-k+1$ で零になる。列は $L,\Lambda$ の両方で不変なので、その直交補も両方で不変である。残った空間の最低次数の非零ベクトルは $\Lambda$ で消える。有限次元なのでこの操作で直和が尽くされる。各列では $L^{n-k}$ は相補的な次数を非零係数で結ぶため同型である。
:::

:::example[複素次元二の原始成分]
$\omega_j=i\theta_j\wedge\bar\theta_j$ とすると $\Lambda\omega_j=1$ である。従って $3\omega_1+\omega_2=2\omega+(\omega_1-\omega_2)$ であり、後ろの項が原始成分である。その外積 $\omega\wedge(\omega_1-\omega_2)$ は零になる。
:::

### Lefschetz の列のノルムと原始射影

次数 $k$ の原始形式 $\alpha$ に対し、前の漸化式を $j$ 回掛けると
$$
\|L^j\alpha\|^2=j!\frac{(n-k)!}{(n-k-j)!}\|\alpha\|^2,\qquad 0\le j\le n-k.
$$
原始形式を互いに直交する正規化した基底から選ぶと、同じ次数の別の列の内積は $\Lambda$ を移して繰り返し下げれば零になる。異なる次数の形式はもともと直交する。次数二では $\Lambda\omega=n$ であるから、原始射影は $\alpha-(\Lambda\alpha/n)\omega$ と明示できる。

:::example[複素三次元の対角形式]
$\alpha=4\omega_1+\omega_2-2\omega_3$ の縮約は $3$ なので、原始部分は $3\omega_1-3\omega_3$ である。これを $\omega^2$ と外積すると、$\omega_1\omega_2\omega_3$ の係数は $2(3-3)=0$ になる。一方 $\alpha$ 自身の $\omega^2$ との外積の係数は $2(4+1-2)=6$ である。
:::

### Lefschetz の作用を基底で追う

:::example[対角形式の外積作用素を計算する]
複素次元 $4$ の一点で正則余接枠 $\eta_j$ を取り、$a_j=i\eta_j\wedge\bar\eta_j$ とする。$a_j^2=0$ である。最高次形式 $\beta$ の $a_1\wedge\cdots\wedge a_n$ に関する係数を $[\beta]_{\mathrm{top}}$ と書く。基本形式には次の正係数の線形結合を使う。

$$
\begin{aligned}
\omega&=a_{1} + a_{2} + a_{3} + a_{4}
\end{aligned}
$$

外積作用素 $L=\omega\wedge$ の $n-2$ 乗は、$a_j$ の張る実空間を $n-1$ 個の添字の積が張る空間に写す。上位側の基底を昇順の添字の積とし、$a_i$ に $\omega^{n-2}$ を掛けると、同じ添字に二度出会う項が消える。残りの各添字集合の係数には並べ替えの数 $(n-2)!$ が掛かる。

$$
\begin{aligned}
[L^{2}]&=\left[\begin{matrix}2 & 2 & 2 & 0\\2 & 2 & 0 & 2\\2 & 0 & 2 & 2\\0 & 2 & 2 & 2\end{matrix}\right] \\
\det[L^{2}]&=-48
\end{aligned}
$$

最高次の係数 $Q(a,b)=[a\wedge b\wedge\omega^{n-2}]_{\mathrm{top}}$ は、行と列が異なるときに残りの添字の重みの積から求まる。$Q(a,\omega)=0$ となる空間を $P$ とし、基底を列に並べた行列を $N$ とする。

$$
\begin{aligned}
[Q]&=\left[\begin{matrix}0 & 2 & 2 & 2\\2 & 0 & 2 & 2\\2 & 2 & 0 & 2\\2 & 2 & 2 & 0\end{matrix}\right],\qquad Q(\omega,\omega)=24 \\
N&=\left[\begin{matrix}-1 & -1 & -1\\1 & 0 & 0\\0 & 1 & 0\\0 & 0 & 1\end{matrix}\right] \\
Q|_P&=\left[\begin{matrix}-4 & -2 & -2\\-2 & -4 & -2\\-2 & -2 & -4\end{matrix}\right]
\end{aligned}
$$

$Q(a,\omega)=0$ は $a\wedge\omega^{n-1}=0$ と同値である。実際、$\omega=\sum_j w_j a_j$ とすると $w_j a_j$ がユニタリ枠に対応する基本形式で、$\Lambda a_j=1/w_j$ である。一方 $Q(a_j,\omega)=(n-1)!\prod_{k\ne j}w_k$ であり、この値は $1/w_j$ の一定倍である。従って $P$ は原始成分に一致する。行列 $Q|_P=N^\top QN$ は、原始成分の外積の符号を与える。
:::

:::example[対角形式の外積作用素を計算する]
複素次元 $3$ の一点で正則余接枠 $\eta_j$ を取り、$a_j=i\eta_j\wedge\bar\eta_j$ とする。$a_j^2=0$ である。最高次形式 $\beta$ の $a_1\wedge\cdots\wedge a_n$ に関する係数を $[\beta]_{\mathrm{top}}$ と書く。基本形式には次の正係数の線形結合を使う。

$$
\begin{aligned}
\omega&=a_{1} + 2 a_{2} + a_{3}
\end{aligned}
$$

外積作用素 $L=\omega\wedge$ の $n-2$ 乗は、$a_j$ の張る実空間を $n-1$ 個の添字の積が張る空間に写す。上位側の基底を昇順の添字の積とし、$a_i$ に $\omega^{n-2}$ を掛けると、同じ添字に二度出会う項が消える。残りの各添字集合の係数には並べ替えの数 $(n-2)!$ が掛かる。

$$
\begin{aligned}
[L^{1}]&=\left[\begin{matrix}2 & 1 & 0\\1 & 0 & 1\\0 & 1 & 2\end{matrix}\right] \\
\det[L^{1}]&=-4
\end{aligned}
$$

最高次の係数 $Q(a,b)=[a\wedge b\wedge\omega^{n-2}]_{\mathrm{top}}$ は、行と列が異なるときに残りの添字の重みの積から求まる。$Q(a,\omega)=0$ となる空間を $P$ とし、基底を列に並べた行列を $N$ とする。

$$
\begin{aligned}
[Q]&=\left[\begin{matrix}0 & 1 & 2\\1 & 0 & 1\\2 & 1 & 0\end{matrix}\right],\qquad Q(\omega,\omega)=12 \\
N&=\left[\begin{matrix}- \frac{1}{2} & -1\\1 & 0\\0 & 1\end{matrix}\right] \\
Q|_P&=\left[\begin{matrix}-1 & -1\\-1 & -4\end{matrix}\right]
\end{aligned}
$$

$Q(a,\omega)=0$ は $a\wedge\omega^{n-1}=0$ と同値である。実際、$\omega=\sum_j w_j a_j$ とすると $w_j a_j$ がユニタリ枠に対応する基本形式で、$\Lambda a_j=1/w_j$ である。一方 $Q(a_j,\omega)=(n-1)!\prod_{k\ne j}w_k$ であり、この値は $1/w_j$ の一定倍である。従って $P$ は原始成分に一致する。行列 $Q|_P=N^\top QN$ は、原始成分の外積の符号を与える。
:::

:::example[曲率交換子の正・零・負の成分を分ける]
正規直交枠で直線束の曲率の固有値を $(\lambda_1,\ldots,\lambda_n)=(1,1,1)$ とする。$(p,q)$ 型形式の基底元を、正則添字 $I$ と反正則添字 $J$ で示す。

$iF\wedge=\sum_j\lambda_j L_j$、$\Lambda=\sum_j\Lambda_j$ とおくと、異なる添字同士の交換子は零である。$L_j,\Lambda_j$ が作用する部分に添字が二つあるとき、$\Lambda_j$ を先に作用させれば元へ戻り、$L_j$ を先に作用させれば零となる。添字がないときはその二順序が逆転し、片方だけのときには両順序が零である。

$$
\begin{aligned}
\lambda_{\mathrm{op}}(I,J)&=\sum_{j\in I\cap J}\lambda_j-\sum_{j\notin I\cup J}\lambda_j \\
I=\{1,2,3\},\ J=\{1,2\}&\colon\quad 2-0=2 \\
I=\{1\},\ J=\{2\}&\colon\quad 0-1=-1 \\
I=\{1\},\ J=\{1\}&\colon\quad 1-2=-1 \\
I=\{\},\ J=\{\}&\colon\quad 0-3=-3
\end{aligned}
$$

従って $I$ が全添字なら負の和が消え、固有値は $J$ の和になる。このとき $q$ が固定されているなら、その型全体の下限は固有値の小さい方から $q$ 個の和である。任意の $I,J$ を許すと負の寄与も出現する。ここで $iF=\omega$ と見れば交換子は $[L,\Lambda]$ である。全次数 $k$ によって固有値が $k-n$ になる。
:::

## Kähler 恒等式

滑らかな形式の $L^2$ 内積について $\partial^*,\bar\partial^*$ を形式随伴とする。ここでの等式は局所微分作用素の等式であり、コンパクト性は不要である。

:::theorem[Kähler 恒等式]{#thm-cg-kahler-identities}
Kähler 多様体では
$$
[\Lambda,\partial]=i\bar\partial^*,\qquad
[\Lambda,\bar\partial]=-i\partial^*
\label{eq-cg-kahler-commutators}
$$
が成り立つ。また $\Delta_\dd=2\Delta_\partial=2\Delta_{\bar\partial}$ である。ここで $\Delta_D=DD^*+D^*D$ とする。
:::

:::proof
任意の点に :ref[prop-cg-kahler-normal] の座標を取る。計量の一次微分が零なので、その点では $\partial=\sum_j e_j\partial_j$、$\bar\partial=\sum_j\bar e_j\partial_{\bar j}$、$\partial^*=-\sum_j\iota_j\partial_{\bar j}$、$\bar\partial^*=-\sum_j\bar\iota_j\partial_j$ と計算できる。随伴の式は、複素微分を積分で移すと共役微分に負号がつくことによる。$\Lambda$ の係数の一次微分も零である。外積と収縮の反交換により $[\Lambda,e_j]=-i\bar\iota_j$、$[\Lambda,\bar e_j]=i\iota_j$ なので最初の二式が従う。

$\{A,B\}=AB+BA$ とする。$\partial^2=\bar\partial^2=0$ により $\{\partial,[\Lambda,\partial]\}=\{\bar\partial,[\Lambda,\bar\partial]\}=0$ である。従って $\{\partial,\bar\partial^*\}=\{\bar\partial,\partial^*\}=0$ となり、$\Delta_\dd=\Delta_\partial+\Delta_{\bar\partial}$ である。また
$$
\{\partial,[\Lambda,\bar\partial]\}+\{\bar\partial,[\Lambda,\partial]\}
=[\Lambda,\{\partial,\bar\partial\}]=0
$$
なので $-i\Delta_\partial+i\Delta_{\bar\partial}=0$ となる。
:::

$[L,\partial]=[L,\bar\partial]=0$ は $\omega$ が各微分で閉じていることから従う。さらに :eqref[eq-cg-kahler-commutators] の随伴を取ると $[L,\bar\partial^*]=-i\partial$ となる。よって $[L,\Delta_{\bar\partial}]=-i\{\bar\partial,\partial\}=0$ であり、随伴で $[\Lambda,\Delta_{\bar\partial}]=0$ も成り立つ。次章の調和形式にはこの両作用が降りる。

### 閉性がなくなると交換子に項が加わる

一般の Hermite 計量では $\omega$ は正定値な $(1,1)$ 形式であり、外積作用素 $L$ とその随伴は定義できる。それでも $[L,\partial]\alpha=-(\partial\omega)\wedge\alpha$ である。これは $\partial(\omega\wedge\alpha)=(\partial\omega)\wedge\alpha+\omega\wedge\partial\alpha$ から従う。例えば前節の $\omega=i(1+|z_2|^2)\dd z_1\wedge\dd\bar z_1+i\dd z_2\wedge\dd\bar z_2$ では $[L,\partial]1=-i\bar z_2\dd z_2\wedge\dd z_1\wedge\dd\bar z_1$ であり、一般点で非零である。従って Lefschetz の線形代数だけから $L$ が調和形式を保つことは従わず、微分作用素の交換には Kähler 条件が要る。

### 交換子を二次形式に作用させる

平坦座標で $\partial^*$ は $-\sum_j\iota_j\partial_{\bar j}$ である。この式では微分と収縮を両方行う。例えば関数だけに $\partial^*$ を作用させても次数が負になるため零である。一方、関数 $f$ に $\partial$ を掛けた後なら
$$
\Delta_\partial f=-\sum_j \partial_{\bar j}\partial_j f.
$$
同様に $\Delta_{\bar\partial} f$ も同じである。実計量が $2$ 倍の Euclid 計量であるため、$\Delta_{\mathrm d}$ は正の規約で $-\frac12\sum_j(\partial_{x_j}^2+\partial_{y_j}^2)$ となり、二倍の恒等式と一致する。

:::example[平坦トーラスでの固有値]
$z=x+iy$、周期 $1$ の両実座標を持つ平坦トーラスで、$f=e^{2\pi i(2x-y)}$ とする。$\Delta_\partial f=5\pi^2 f$ である。実 Laplacian では $\partial_x^2f=-16\pi^2f$、$\partial_y^2f=-4\pi^2f$ だから $\Delta_{\mathrm d}f=10\pi^2 f$ となる。
:::

### 平坦な Kähler 計量を使った微分

:::example[平均と非零 Fourier モードを分離する]
$X=\C/(\Z+2i\Z)$ 上で $\alpha=a\,\mathrm d\bar z$ とし、周期座標を $z=x+iy$ とする。係数と微分の規約は次である。

$$
\begin{aligned}
a(x,y)&=2 e^{2 i \pi \left(x + \frac{y}{2}\right)} + e^{2 i \pi \left(2 x - \frac{y}{2}\right)} + 3 e^{- i \pi y} \\
\partial_{\bar z}&=\tfrac12(\partial_x+i\partial_y)
\end{aligned}
$$

周期積分では全微分の平均が零になる。したがって定数モードは $\bar\partial$ の方程式の障害である。他方、指数関数 $e^{2\pi i(mx+ny/L)}$ の二つの微分で $2\pi im,2\pi in/L$ が掛かるので、固有値は $\pi(im-n/L)$ になる。

$$
\begin{aligned}
(1,1) &\colon\quad \lambda=\pi \left(- \frac{1}{2} + i\right),\quad a=2,\quad a/\lambda=\frac{4 \left(-1 - 2 i\right)}{5 \pi} \\
(2,-1) &\colon\quad \lambda=\pi \left(\frac{1}{2} + 2 i\right),\quad a=1,\quad a/\lambda=\frac{2 \left(1 - 4 i\right)}{17 \pi} \\
(0,-1) &\colon\quad \lambda=\frac{\pi}{2},\quad a=3,\quad a/\lambda=\frac{6}{\pi} \\
a_{00}&=0
\end{aligned}
$$

定数モードを除いて各係数を固有値で割る。これで得られる関数は平均零で、微分を取り直すと次を満たす。

$$
\begin{aligned}
u&=\frac{8 \left(- \frac{1}{2} - i\right) e^{2 i \pi \left(x + \frac{y}{2}\right)}}{5 \pi} + \frac{4 \left(\frac{1}{2} - 2 i\right) e^{2 i \pi \left(2 x - \frac{y}{2}\right)}}{17 \pi} + \frac{6 e^{- i \pi y}}{\pi} \\
u_{\bar z}&=a-0
\end{aligned}
$$

全てのモードの周期が元のトーラスと一致するため $u$ は大域的である。Kähler 恒等式で関数の二種類の Laplacian の固有値は $2$ 倍だけ異なる。
:::

:::example[平均と非零 Fourier モードを分離する]
$X=\C/(\Z+3i\Z)$ 上で $\alpha=a\,\mathrm d\bar z$ とし、周期座標を $z=x+iy$ とする。係数と微分の規約は次である。

$$
\begin{aligned}
a(x,y)&=- 2 e^{2 i \pi x} + 2 e^{\frac{2 i \pi y}{3}} + 2 e^{2 i \pi \left(x - \frac{2 y}{3}\right)} \\
\partial_{\bar z}&=\tfrac12(\partial_x+i\partial_y)
\end{aligned}
$$

周期積分では全微分の平均が零になる。したがって定数モードは $\bar\partial$ の方程式の障害である。他方、指数関数 $e^{2\pi i(mx+ny/L)}$ の二つの微分で $2\pi im,2\pi in/L$ が掛かるので、固有値は $\pi(im-n/L)$ になる。

$$
\begin{aligned}
(0,1) &\colon\quad \lambda=- \frac{\pi}{3},\quad a=2,\quad a/\lambda=- \frac{6}{\pi} \\
(1,-2) &\colon\quad \lambda=\pi \left(\frac{2}{3} + i\right),\quad a=2,\quad a/\lambda=\frac{6 \left(2 - 3 i\right)}{13 \pi} \\
(1,0) &\colon\quad \lambda=i \pi,\quad a=-2,\quad a/\lambda=\frac{2 i}{\pi} \\
a_{00}&=0
\end{aligned}
$$

定数モードを除いて各係数を固有値で割る。これで得られる関数は平均零で、微分を取り直すと次を満たす。

$$
\begin{aligned}
u&=\frac{2 i e^{2 i \pi x}}{\pi} - \frac{6 e^{\frac{2 i \pi y}{3}}}{\pi} + \frac{18 \left(\frac{2}{3} - i\right) e^{2 i \pi \left(x - \frac{2 y}{3}\right)}}{13 \pi} \\
u_{\bar z}&=a-0
\end{aligned}
$$

全てのモードの周期が元のトーラスと一致するため $u$ は大域的である。複素微分を一つ解くことと、実 Laplace 方程式を解くことは同じ操作ではない。
:::

### 外積と収縮の二つの順序

:::example[曲率交換子の正・零・負の成分を分ける]
正規直交枠で直線束の曲率の固有値を $(\lambda_1,\ldots,\lambda_n)=(1,1,1,1)$ とする。$(p,q)$ 型形式の基底元を、正則添字 $I$ と反正則添字 $J$ で示す。

$iF\wedge=\sum_j\lambda_j L_j$、$\Lambda=\sum_j\Lambda_j$ とおくと、異なる添字同士の交換子は零である。$L_j,\Lambda_j$ が作用する部分に添字が二つあるとき、$\Lambda_j$ を先に作用させれば元へ戻り、$L_j$ を先に作用させれば零となる。添字がないときはその二順序が逆転し、片方だけのときには両順序が零である。

$$
\begin{aligned}
\lambda_{\mathrm{op}}(I,J)&=\sum_{j\in I\cap J}\lambda_j-\sum_{j\notin I\cup J}\lambda_j \\
I=\{1\},\ J=\{2,3\}&\colon\quad 0-1=-1 \\
I=\{1,4\},\ J=\{1,2\}&\colon\quad 1-1=0 \\
I=\{1,2,3\},\ J=\{3,4\}&\colon\quad 1-0=1 \\
I=\{1,2,3,4\},\ J=\{1,4\}&\colon\quad 2-0=2
\end{aligned}
$$

従って $I$ が全添字なら負の和が消え、固有値は $J$ の和になる。このとき $q$ が固定されているなら、その型全体の下限は固有値の小さい方から $q$ 個の和である。任意の $I,J$ を許すと負の寄与も出現する。これは零階の交換子である。微分作用素を混ぜる Kähler 恒等式には、さらに計量の一次微分を消す正規座標が必要になる。
:::

## Kähler 多様体の曲率

Kähler 正規座標の計量で、Chern 接続の公式を接束に適用すると曲率の成分は
$$
R_{i\bar j k\bar\ell}
=-\partial_i\partial_{\bar j}g_{k\bar\ell}
+\sum_{a,b}g^{a\bar b}
  (\partial_i g_{k\bar b})(\partial_{\bar j}g_{a\bar\ell})
$$
となる。ここで逆行列の添字は縮約が行列の逆になるように置いている。実際、:ref[thm-cg-chern-connection] の $H^{-1}\partial H$ を $\bar\partial$ で微分し、$\dd\bar z_j\wedge\dd z_i$ を逆順にする。逆行列の微分で積の項が生じる。Kähler 条件で $i,k$、共役条件で $j,\ell$ を入れ替えても式は変わらない。

曲率の跡に $i$ を掛けたものが :term[Ricci 形式]{en="Ricci form" yomi="りっちけいしき"} $\rho$ である。$\tr(H^{-1}\partial H)=\partial\log\det H$ なので
$$
\rho=-i\partial\bar\partial\log\det(g_{j\bar k}),\qquad
[\rho/(2\pi)]=c_1(T_X).
\label{eq-cg-ricci-form}
$$
従って Ricci 形式が正であることと、個々の全曲率成分が正であることは異なる条件である。

:::example[Fubini–Study の Ricci 形式]
複素 Hessian の固有値は $z$ の方向で $(1+|z|^2)^{-2}$、直交方向で $(1+|z|^2)^{-1}$ である。よって行列式は $(1+|z|^2)^{-(n+1)}$、従って $\rho=(n+1)\omega_{\mathrm{FS}}$ である。一方平坦トーラスでは $\det g$ は定数であり $\rho=0$ である。
:::

### Fubini–Study 計量の正則断面曲率

複素直線方向の実二次元平面における断面曲率を正則断面曲率という。$\omega_{\mathrm{FS}}$ は $\C^{n+1}$ の標準 Hermite 内積からできるのでユニタリ変換で不変である。点とそこでの複素直線方向は、原点 $[1:0:\cdots:0]$ と $z_1$ 軸方向へ移せる。射影直線 $\{Z_2=\cdots=Z_n=0\}$ は残りの座標を一斉に符号反転する等長変換の固定集合の一成分である。この直線に接して出発する測地線を等長変換で送ると同じ初期条件を持つので、測地線の一意性により固定される。従って射影直線は全測地的であり、この方向の断面曲率は誘導計量の Gauss 曲率と等しい。

その直線の実計量は $2(1+|z|^2)^{-2}(\dd x^2+\dd y^2)$ である。$u=(\log2)/2-\log(1+x^2+y^2)$ とすれば $\Delta u=-4/(1+x^2+y^2)^2$ であり、曲率 $-e^{-2u}\Delta u$ は $2$ である。従ってこの規格化では全ての正則断面曲率が $2$ になる。計量を $c>0$ 倍すると Levi-Civita 接続は同じで、断面曲率は $1/c$ 倍になるため、例えば面積を $\pi$ に規格化した射影直線では値は $4$ になる。

### 部分多様体の曲率と第二基本形式

正則部分多様体 $Y\subset X$ に誘導計量を入れる。:ref[prop-cg-subbundle-curvature] を $T_Y\subset T_X|_Y$ に適用すると、Chern 接続と Levi-Civita 接続の一致により Kähler 版の Gauss の式を得る。$B$ を接方向から法方向への第二基本形式とし、$u,v\in T^{1,0}Y$ を取る。$(B^*\wedge B)(u,\bar u)=-B(u)^*B(u)$ なので、
$$
R_Y(u,\bar u,v,\bar v)=R_X(u,\bar u,v,\bar v)-\|B(u)v\|^2
$$
となる。ここで $R$ は本節の成分規約を使う。$v=u$ とすると分子が正則断面曲率の分子になる。環境と部分多様体で $u$ のノルムは同じなので、誘導した正則断面曲率は環境の同じ方向の値以下である。等号になるかどうかは、ここでは $B(u)u$ が零かどうかで決まる。

:::example[平坦空間内の放物線]
$F(z)=(z,az^2)$ を標準の平坦計量を持つ $\C^2$ に入れる。原点では $F'(0)=(1,0)$ は単位の正則接ベクトルで、$F''(0)=(0,2a)$ は法方向である。従って $B(\partial_z)\partial_z=(0,2a)$ で、Gauss の式は原点の曲率を $-4|a|^2$ と与える。同じ値は $g_{z\bar z}=1+4|a|^2|z|^2$ を直接微分しても得られる。
:::

### 積計量の Ricci 形式と Einstein 条件

各因子の正則座標を合わせると、積計量の行列式は各因子の行列式の積である。従って $-\partial\bar\partial$ をその対数に作用させると、Ricci 形式は各因子の引き戻しの和になる。定数による計量の倍は対数行列式に定数を加えるだけなので、Ricci 形式は変わらない。従って二つの Kähler–Einstein 因子の積を Einstein 計量にするには、両因子の Einstein 定数が同じになるよう倍率を選ぶ必要がある。

:::example[二つの射影空間で倍率を合わせる]
$\mathbb P^1\times\mathbb P^2$ で $\omega=a\omega_1+b\omega_2$ とする。$\rho=2\omega_1+3\omega_2$ なので $\rho=\lambda\omega$ となる条件は $\lambda a=2,\lambda b=3$ である。例えば $(a,b)=(2,3)$ なら $\lambda=1$、$(a,b)=(1,1)$ なら一つの $\lambda$ は存在しない。曲率が正の因子の積でも Einstein になるかは倍率によって異なる。
:::

### 曲率を持つ計量の比較

:::example[ポテンシャルから体積と Ricci を追う]
$\varphi=2 \bar z_1^{2} z_{1}^{2} + \bar z_1 z_{1} + 3 \bar z_2^{2} z_{2}^{2} + \bar z_2 z_{2}$ とし、点 $p=(0,0)$ で複素 Hessian と Ricci 形式の係数を求める。$\partial_{z_j}$ は $\bar z_k$ を固定して微分するため、計算中は両者を独立な変数として扱う。

$$
\begin{aligned}
G=(\varphi_{j\bar k})&=\left[\begin{matrix}8 \bar z_1 z_{1} + 1 & 0\\0 & 12 \bar z_2 z_{2} + 1\end{matrix}\right] \\
\det G&=\left(8 \bar z_1 z_{1} + 1\right) \left(12 \bar z_2 z_{2} + 1\right) \\
G(p)&=\left[\begin{matrix}1 & 0\\0 & 1\end{matrix}\right],\qquad\det G(p)=1
\end{aligned}
$$

点での先頭小行列式 $G_{1\bar1}(p)=1$ と行列式がともに正なので近傍で Kähler 計量になる。体積密度は $4\det G$ である。Ricci 形式の係数 $R$ は $-\partial_j\partial_{\bar k}\log\det G$ だから、商の微分で次の式を得る。

$$
\begin{aligned}
R_{j\bar k}&=-\frac{d_{j\bar k}}d+\frac{d_jd_{\bar k}}{d^2}\qquad(d=\det G) \\
R(p)&=\left[\begin{matrix}-8 & 0\\0 & -12\end{matrix}\right] \\
\Lambda_\omega\rho(p)&=\tr(G(p)^{-1}R(p))=-20
\end{aligned}
$$

対数行列式の二階微分は、点での行列式を先に数へ置き換えると失われる。そこで行列式を関数として微分してから、その点の値を代入した。原点の Hessian は単位行列だが、Ricci は平坦計量と一致しない。
:::

:::example[ポテンシャルから体積と Ricci を追う]
$\varphi=\bar z_1 z_{1} + \bar z_2 z_{2} + \left(\bar z_1 + \bar z_2\right)^{2} \left(z_{1} + z_{2}\right)^{2}$ とし、点 $p=(0,0)$ で複素 Hessian と Ricci 形式の係数を求める。$\partial_{z_j}$ は $\bar z_k$ を固定して微分するため、計算中は両者を独立な変数として扱う。

$$
\begin{aligned}
G=(\varphi_{j\bar k})&=\left[\begin{matrix}4 \left(\bar z_1 + \bar z_2\right) \left(z_{1} + z_{2}\right) + 1 & 4 \left(\bar z_1 + \bar z_2\right) \left(z_{1} + z_{2}\right)\\4 \left(\bar z_1 + \bar z_2\right) \left(z_{1} + z_{2}\right) & 4 \left(\bar z_1 + \bar z_2\right) \left(z_{1} + z_{2}\right) + 1\end{matrix}\right] \\
\det G&=8 \bar z_1 z_{1} + 8 \bar z_1 z_{2} + 8 \bar z_2 z_{1} + 8 \bar z_2 z_{2} + 1 \\
G(p)&=\left[\begin{matrix}1 & 0\\0 & 1\end{matrix}\right],\qquad\det G(p)=1
\end{aligned}
$$

点での先頭小行列式 $G_{1\bar1}(p)=1$ と行列式がともに正なので近傍で Kähler 計量になる。体積密度は $4\det G$ である。Ricci 形式の係数 $R$ は $-\partial_j\partial_{\bar k}\log\det G$ だから、商の微分で次の式を得る。

$$
\begin{aligned}
R_{j\bar k}&=-\frac{d_{j\bar k}}d+\frac{d_jd_{\bar k}}{d^2}\qquad(d=\det G) \\
R(p)&=\left[\begin{matrix}-8 & -8\\-8 & -8\end{matrix}\right] \\
\Lambda_\omega\rho(p)&=\tr(G(p)^{-1}R(p))=-16
\end{aligned}
$$

対数行列式の二階微分は、点での行列式を先に数へ置き換えると失われる。そこで行列式を関数として微分してから、その点の値を代入した。同じ原点の計量に対して、非対角の Ricci 形式が生じる場合である。
:::

:::example[ポテンシャルから体積と Ricci を追う]
$\varphi=- \log{\left(- \bar z_1 z_{1} + 1 \right)} + \log{\left(\bar z_2 z_{2} + 1 \right)}$ とし、点 $p=(0,1)$ で複素 Hessian と Ricci 形式の係数を求める。$\partial_{z_j}$ は $\bar z_k$ を固定して微分するため、計算中は両者を独立な変数として扱う。

$$
\begin{aligned}
G=(\varphi_{j\bar k})&=\left[\begin{matrix}\frac{\frac{\bar z_1 z_{1}}{\bar z_1 z_{1} - 1} - 1}{\bar z_1 z_{1} - 1} & 0\\0 & \frac{- \frac{\bar z_2 z_{2}}{\bar z_2 z_{2} + 1} + 1}{\bar z_2 z_{2} + 1}\end{matrix}\right] \\
\det G&=\frac{1}{\left(\bar z_1 z_{1} - 1\right)^{2} \left(\bar z_2 z_{2} + 1\right)^{2}} \\
G(p)&=\left[\begin{matrix}1 & 0\\0 & \frac{1}{4}\end{matrix}\right],\qquad\det G(p)=\frac{1}{4}
\end{aligned}
$$

点での先頭小行列式 $G_{1\bar1}(p)=1$ と行列式がともに正なので近傍で Kähler 計量になる。体積密度は $4\det G$ である。Ricci 形式の係数 $R$ は $-\partial_j\partial_{\bar k}\log\det G$ だから、商の微分で次の式を得る。

$$
\begin{aligned}
R_{j\bar k}&=-\frac{d_{j\bar k}}d+\frac{d_jd_{\bar k}}{d^2}\qquad(d=\det G) \\
R(p)&=\left[\begin{matrix}-2 & 0\\0 & \frac{1}{2}\end{matrix}\right] \\
\Lambda_\omega\rho(p)&=\tr(G(p)^{-1}R(p))=0
\end{aligned}
$$

対数行列式の二階微分は、点での行列式を先に数へ置き換えると失われる。そこで行列式を関数として微分してから、その点の値を代入した。第一因子では単位円板の範囲を取る。二因子の Ricci の符号が異なるため、同じ定数による Einstein 条件を満たさない。
:::

:::example[平坦トーラスの体積指定を積の形で解く]
複素次元 $2$ の正方トーラスで、座標 $z_j=x_j+iy_j$ の両実周期を $1$ とする。初期の計量行列を $I$ とし、体積比として次の積を指定する。

$$
\begin{aligned}
e^f&=\left(\frac{\cos{\left(2 \pi x_{1} \right)}}{4} + 1\right) \left(\frac{\cos{\left(2 \pi x_{2} \right)}}{3} + 1\right)
\end{aligned}
$$

各因子は正で、周期積分では余弦の平均が零なので因子ごとの平均は $1$ である。従って積の全平均も $1$ となる。ポテンシャルは各一変数の関数を足した形に取る。混合偏微分が消えるため、行列式が指定の積になるには対角成分だけを調整すればよい。

$$
\begin{aligned}
\varphi&=- \frac{\cos{\left(2 \pi x_{1} \right)}}{4 \pi^{2}} - \frac{\cos{\left(2 \pi x_{2} \right)}}{3 \pi^{2}} \\
g_\varphi&=\left[\begin{matrix}\frac{\cos{\left(2 \pi x_{1} \right)}}{4} + 1 & 0\\0 & \frac{\cos{\left(2 \pi x_{2} \right)}}{3} + 1\end{matrix}\right] \\
\det(g_\varphi)&=\left(\frac{\cos{\left(2 \pi x_{1} \right)}}{4} + 1\right) \left(\frac{\cos{\left(2 \pi x_{2} \right)}}{3} + 1\right)
\end{aligned}
$$

各方向について $u_{z_j\bar z_j}=u_{x_jx_j}/4$ を計算すれば対角成分を再現する。従って新しい複素 Hessian は正定値で、同じ Kähler 類にある。さらに対数行列式の二階微分から、新しい Ricci 形式の係数は次の対角行列になる。

$$
\begin{aligned}
R_{1\overline{1}}&=\frac{\pi^{2} \left(\sin^{2}{\left(2 \pi x_{1} \right)} + \cos^{2}{\left(2 \pi x_{1} \right)} + 4 \cos{\left(2 \pi x_{1} \right)}\right)}{\left(\cos{\left(2 \pi x_{1} \right)} + 4\right)^{2}} \\
R_{2\overline{2}}&=\frac{\pi^{2} \left(\sin^{2}{\left(2 \pi x_{2} \right)} + \cos^{2}{\left(2 \pi x_{2} \right)} + 3 \cos{\left(2 \pi x_{2} \right)}\right)}{\left(\cos{\left(2 \pi x_{2} \right)} + 3\right)^{2}}
\end{aligned}
$$

新しい体積が旧体積と同じ平均を持つことと、Ricci が点ごとに零であることは異なる。各 $\lambda_j$ の正値性から計量は定義されるが、一般の点で Ricci は零にならない。Ricci 形式は完全で総第一類が零でも、局所の値は零でない。
:::

## 演習問題

::::exercise{level=basic}
$\varphi=2|z|^2+3|w|^2+\Re(zw)$ の複素 Hessian を求めよ。
:::solution
$\Re(zw)$ は正則関数の実部なので混合微分が零である。二つのノルム二乗の混合微分は対角成分だけを与え、Hessian は $\diag(2,3)$ である。従ってこのポテンシャルは Kähler 形式を与える。
:::
::::

::::exercise{level=basic}
$n=3$ で次数 $2$ の形式に対する $[\Lambda,L]$ の値を求めよ。
:::solution
交換子は $(n-k)\id$ なので $n=3,k=2$ を代入して恒等作用素となる。
:::
::::

::::exercise{level=basic}
$\mathbb P^2$ の Fubini–Study 計量の Ricci 形式を $\omega_{\mathrm{FS}}$ で表せ。
:::solution
$n=2$ の Hessian の行列式は $(1+|z|^2)^{-3}$ である。対数を取り $-i\partial\bar\partial$ を作用させると $3i\partial\bar\partial\log(1+|z|^2)=3\omega_{\mathrm{FS}}$ となる。
:::
::::

::::exercise{level=basic}
$\omega=i(e^{|w|^2}\dd z\wedge\dd\bar z+\dd w\wedge\dd\bar w)$ が Kähler であるか判定せよ。
:::solution
$g_{z\bar z}=e^{|w|^2}$ なので $\partial_wg_{z\bar z}=\bar w e^{|w|^2}$ である。一方 $g_{w\bar z}=0$ で、その $z$ 微分は零である。両者は開領域上で一致しないので閉形式ではなく、計量は Kähler ではない。
:::
::::

::::exercise{level=advanced}
コンパクトな Kähler 多様体の $[\omega]^k$ が $0\le k\le n$ で非零であることを示せ。
:::solution
$k=0$ は定数 $1$ の類で非零である。もし $\omega^k=\dd\eta$、$1\le k\le n$ なら $\omega^{n-k}$ は閉じているので $\omega^n=\dd(\eta\wedge\omega^{n-k})$ となる。Stokes の定理で積分は零だが、$\omega^n/n!$ は正の体積形式なので矛盾する。
:::
::::

::::exercise{level=advanced}
複素次元 $3$ で原始 $1$ 形式 $\alpha$ に対して $\Lambda L^2\alpha$ と $L^3\alpha$ を求めよ。
:::solution
:eqref[eq-cg-lefschetz-string] に $k=1,j=2$ を入れると係数は $2(3-1-2+1)=2$ なので $\Lambda L^2\alpha=2L\alpha$ である。列の最後は $L^{n-k}\alpha=L^2\alpha$ なので $L^3\alpha=0$ である。
:::
::::

::::exercise{level=advanced}
任意の Kähler 計量について $\Delta_\dd$ が $(p,q)$ 型を保つことを示せ。
:::solution
$\bar\partial$ は第二次数を上げ、$\bar\partial^*$ は下げる。よって $\Delta_{\bar\partial}$ の二項はいずれも各型を保つ。:ref[thm-cg-kahler-identities] で $\Delta_\dd=2\Delta_{\bar\partial}$ だから実外微分の Laplacian も各型を保つ。
:::
::::

::::exercise{level=advanced}
$\omega_{\mathrm{FS}}$ を $c>0$ 倍した計量の Ricci 形式と、Ricci 形式の新しい Kähler 形式に対する比例定数を求めよ。
:::solution
行列式は $c^n$ 倍になり、対数の差 $n\log c$ は定数なので Ricci 形式はもとの $(n+1)\omega_{\mathrm{FS}}$ のままである。新しい形式 $c\omega_{\mathrm{FS}}$ に対する比例定数は $(n+1)/c$ である。
:::
::::

::::exercise{level=advanced}
平坦な $\C^2$ へ $F(z)=(z,cz^2)$、$c\ne0$ で埋め込んだ曲線の誘導 Kähler 計量と原点での Gauss 曲率を求めよ。ここで環境の Kähler 形式は $i\sum_j\mathrm dz_j\wedge\mathrm d\bar z_j$ とする。

:::solution
各成分の引き戻しを加えると複素係数は
$$
h=|1|^2+|2cz|^2=1+4|c|^2|z|^2.
$$
従って実計量は $2h(\mathrm dx^2+\mathrm dy^2)$ である。$u=\log(2h)/2$ とおくと $\Delta u=2\partial_z\partial_{\bar z}\log h$ であり、Gauss 曲率は $-e^{-2u}\Delta u=-h^{-1}\partial_z\partial_{\bar z}\log h$ となる。$A=4|c|^2$ と書けば
$$
\partial_z\partial_{\bar z}\log(1+A|z|^2)=\frac A{(1+A|z|^2)^2}
$$
なので $K=-A/(1+A|z|^2)^3$ である。原点では $-4|c|^2$ となる。環境は平坦でも、埋め込みの二次の成分が誘導曲率を生む。
:::
::::

::::exercise{level=basic}
複素次元 $4$ で $\omega=\omega_1+\cdots+\omega_4$ とする。$\alpha=\omega_1+\omega_2-2\omega_3$ が原始であることを確かめ、$\Lambda L\alpha$ を求めよ。

:::solution
各 $\Lambda\omega_j=1$ なので $\Lambda\alpha=1+1-2=0$ である。全次数は $2$ であり、交換子 $[\Lambda,L]=(4-2)\id$ を使うと $\Lambda L\alpha=2\alpha+L\Lambda\alpha=2\alpha$ となる。
:::
::::

## 確認テスト

::::quiz{answer="3" tol="1e-6"}
複素次元 $3$ で $\Lambda\omega$ の値を求めよ。
:::explanation
$\Lambda1=0$ なので $\Lambda L1=[\Lambda,L]1=3$ となる。各ユニタリ座標の対から $1$ ずつ寄与する。
:::
::::

::::quiz
Hermite 計量が Kähler であるための条件を選べ。

- [x] 基本形式 $\omega$ が閉じている。
- [ ] 局所ポテンシャルが定数である。
- [ ] $\omega$ が全域で完全である。

:::explanation
閉性 $\dd\omega=0$ が定義の条件である。局所的には $i\partial\bar\partial\varphi$ と書けるが、コンパクトな場合に全域で完全なら体積積分が零になり矛盾する。
:::
::::

::::quiz{answer="2" tol="1e-6"}
Kähler 多様体で $\Delta_\dd=c\Delta_{\bar\partial}$ となる定数 $c$ を求めよ。
:::explanation
混合する反交換子が Kähler 恒等式で消え、$\Delta_\dd=\Delta_\partial+\Delta_{\bar\partial}$ となる。さらに二つの Laplacian が等しいので $c=2$ である。
:::
::::

::::quiz
$\omega=i\partial\bar\partial\varphi$ のポテンシャルに局所的に加えても $\omega$ を変えない関数を選べ。

- [ ] $|z_1|^2$
- [x] $\Re(z_1^2)$
- [ ] $\Re(z_1\bar z_1)$

:::explanation
$\Re(z_1^2)$ は正則関数の実部で、$\partial\bar\partial$ が消える。残りはいずれも $|z_1|^2$ であり、複素 Hessian の対角成分を $1$ だけ変える。
:::
::::

::::quiz{answer="6" tol="1e-6"}
複素三次元で $\alpha=4\omega_1+\omega_2-2\omega_3$、$\omega=\sum_j\omega_j$ とする。$\alpha\wedge\omega^2$ の $\omega_1\wedge\omega_2\wedge\omega_3$ の係数を求めよ。

:::explanation
$\omega_j^2=0$ であり、異なる二因子の積は順序を二通り選べる。各係数 $4,1,-2$ に $2$ が掛かるため、合計は $2(4+1-2)=6$ である。
:::
::::

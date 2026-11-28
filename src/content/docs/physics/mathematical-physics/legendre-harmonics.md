---
title: "Legendre 関数と球面調和関数"
description: "Legendre 多項式から球面調和関数を構成し、球対称な境界値問題と多重極展開を解く。"
notation:
  - symbol: 'P_\ell'
    meaning: "Legendre 多項式"
  - symbol: '\mathsf P_\ell^m'
    meaning: "区間内の Ferrers 関数"
  - symbol: 'Y_{\ell m}'
    meaning: "正規化した球面調和関数"
  - symbol: 'Q_\ell'
    meaning: "区間外の第 2 種 Legendre 関数"
---

## この章の目標

- Legendre 方程式の多項式解とその規格化を導く。
- 陪関数から球面の正規直交基底を構成する。
- 加法定理から向きに依存する Coulomb 核を展開する。
- 球面上の境界値から内外の調和関数を求める。

## Legendre 方程式と Legendre 多項式

### 極で有限な軸対称解

球面上で軸対称な場は $x=\cos\theta$ だけに依存する。角方向の Laplacian は $\partial_x[(1-x^2)\partial_x]$ なので、$-\Delta_{S^2}$ の固有値を $\ell(\ell+1)$ とした方程式を :term[ルジャンドル方程式]{en="Legendre equation" yomi="るじゃんどるほうていしき"} という。$\ell=0,1,\ldots$ について
$$
P_\ell(x)=\frac1{2^\ell\ell!}\dvn{(x^2-1)^\ell}{x}{\ell}
\label{eq-mp-legendre-rodrigues}
$$
と置く。$F=(x^2-1)^\ell$ は $(x^2-1)F'=2\ell xF$ を満たす。この恒等式を $\ell+1$ 回微分すると、$y=F^{(\ell)}$ に対して $(1-x^2)y''-2xy'+\ell(\ell+1)y=0$ になる。$(x-1)^\ell(x+1)^\ell$ を $x=1$ で微分すると $P_\ell(1)=1$ である。例えば $P_0=1,P_1=x,P_2=(3x^2-1)/2$ であり、偶奇性は $P_\ell(-x)=(-1)^\ell P_\ell(x)$ である。

$n<\ell$ の多項式と $P_\ell$ の積分では、$\ell$ 回部分積分すると他方の微分が零になる。端では $F$ の $\ell-1$ 階以下の微分が全部消える。先頭係数 $c_\ell=(2\ell)!/[2^\ell(\ell!)^2]$ を使い、$P_\ell$ 自身を掛けて同じ計算を行うと、ノルムの二乗は $c_\ell2^{-\ell}\int_{-1}^1(1-x^2)^\ell\dd x=2/(2\ell+1)$ になる。積分は $x=2u-1$ と置くと $2^{2\ell+1}B(\ell+1,\ell+1)$ である。

### 端点の正則性と固有値の選択

軸対称の角方程式で、端点付近の第二解は一般に対数を含む。正則な多項式を選ぶには $((1-x^2)u^\prime)^\prime+\lambda u=0$ の係数を比べる。最高次の項の係数は $\lambda-m(m+1)$ なので、次数 $m$ の非零多項式を持つなら $\lambda=m(m+1)$ である。別の次数の解との積分では、部分積分の境界項に $1-x^2$ が掛かり、正則な解では消える。

:::example[Legendre 方程式の各項（1）]
多項式の二階項と一階項を分けて計算し、最高次の項だけでなく全係数が相殺することを確かめる。

$$
m=1
$$

使用する規格化は $P_m(1)=1$ である。

$$
\begin{aligned}
P_m(x)\big|_{m=1} &= x
\end{aligned}
$$

演算子の二階微分の部分を求める。

$$
\begin{aligned}
(1-x^2)P_m^{\prime\prime}\big|_{m=1} &= 0
\end{aligned}
$$

一階微分の部分と零階項を合わせる。二階部分の反対になる。

$$
\begin{aligned}
-2xP_m^\prime+m(m+1)P_m\big|_{m=1} &= 0
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
P_m(x)\big|_{x=\frac{1}{2}} &= \frac{1}{2} \\
\frac{\partial^1 (P_m(x))}{\partial x^1} &= 1 \\
\left.\frac{\partial^1 (P_m(x))}{\partial x^1}\right|_{x=\frac{1}{2}} &= 1
\end{aligned}
$$
:::

:::example[Legendre 方程式の各項（2）]
多項式の二階項と一階項を分けて計算し、最高次の項だけでなく全係数が相殺することを確かめる。次数を一つ上げる。

$$
m=2
$$

使用する規格化は $P_m(1)=1$ である。

$$
\begin{aligned}
P_m(x)\big|_{m=2} &= \frac{3 x^{2}}{2} - \frac{1}{2}
\end{aligned}
$$

演算子の二階微分の部分を求める。

$$
\begin{aligned}
(1-x^2)P_m^{\prime\prime}\big|_{m=2} &= 3 - 3 x^{2}
\end{aligned}
$$

一階微分の部分と零階項を合わせる。二階部分の反対になる。

$$
\begin{aligned}
-2xP_m^\prime+m(m+1)P_m\big|_{m=2} &= 3 x^{2} - 3
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
P_m(x)\big|_{x=\frac{1}{2}} &= - \frac{1}{8} \\
\frac{\partial^1 (P_m(x))}{\partial x^1} &= 3 x \\
\left.\frac{\partial^1 (P_m(x))}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{3}{2}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
(1-x^2)P_m^{\prime\prime}\big|_{x=\frac{1}{2}} &= \frac{9}{4} \\
\frac{\partial^1 ((1-x^2)P_m^{\prime\prime})}{\partial x^1} &= - 6 x \\
\left.\frac{\partial^1 ((1-x^2)P_m^{\prime\prime})}{\partial x^1}\right|_{x=\frac{1}{2}} &= -3
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
-2xP_m^\prime+m(m+1)P_m\big|_{x=\frac{1}{2}} &= - \frac{9}{4} \\
\frac{\partial^1 (-2xP_m^\prime+m(m+1)P_m)}{\partial x^1} &= 6 x \\
\left.\frac{\partial^1 (-2xP_m^\prime+m(m+1)P_m)}{\partial x^1}\right|_{x=\frac{1}{2}} &= 3
\end{aligned}
$$
:::

:::example[Legendre 方程式の各項（3）]
多項式の二階項と一階項を分けて計算し、最高次の項だけでなく全係数が相殺することを確かめる。さらに次数を上げて、低次の結果と照合する。

$$
m=3
$$

使用する規格化は $P_m(1)=1$ である。

$$
\begin{aligned}
P_m(x)\big|_{m=3} &= \frac{x \left(5 x^{2} - 3\right)}{2}
\end{aligned}
$$

演算子の二階微分の部分を求める。

$$
\begin{aligned}
(1-x^2)P_m^{\prime\prime}\big|_{m=3} &= 15 x \left(1 - x^{2}\right)
\end{aligned}
$$

一階微分の部分と零階項を合わせる。二階部分の反対になる。

$$
\begin{aligned}
-2xP_m^\prime+m(m+1)P_m\big|_{m=3} &= 15 x \left(x^{2} - 1\right)
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
P_m(x)\big|_{x=\frac{1}{2}} &= - \frac{7}{16} \\
\frac{\partial^1 (P_m(x))}{\partial x^1} &= \frac{15 x^{2}}{2} - \frac{3}{2} \\
\left.\frac{\partial^1 (P_m(x))}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{3}{8}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
(1-x^2)P_m^{\prime\prime}\big|_{x=\frac{1}{2}} &= \frac{45}{8} \\
\frac{\partial^1 ((1-x^2)P_m^{\prime\prime})}{\partial x^1} &= 15 - 45 x^{2} \\
\left.\frac{\partial^1 ((1-x^2)P_m^{\prime\prime})}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{15}{4}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
-2xP_m^\prime+m(m+1)P_m\big|_{x=\frac{1}{2}} &= - \frac{45}{8} \\
\frac{\partial^1 (-2xP_m^\prime+m(m+1)P_m)}{\partial x^1} &= 45 x^{2} - 15 \\
\left.\frac{\partial^1 (-2xP_m^\prime+m(m+1)P_m)}{\partial x^1}\right|_{x=\frac{1}{2}} &= - \frac{15}{4}
\end{aligned}
$$
:::

:::example[Legendre 方程式の各項（4）]
多項式の二階項と一階項を分けて計算し、最高次の項だけでなく全係数が相殺することを確かめる。四次の結果まで具体化する。

$$
m=4
$$

使用する規格化は $P_m(1)=1$ である。

$$
\begin{aligned}
P_m(x)\big|_{m=4} &= \frac{35 x^{4}}{8} - \frac{15 x^{2}}{4} + \frac{3}{8}
\end{aligned}
$$

演算子の二階微分の部分を求める。

$$
\begin{aligned}
(1-x^2)P_m^{\prime\prime}\big|_{m=4} &= - \frac{105 x^{4}}{2} + 60 x^{2} - \frac{15}{2}
\end{aligned}
$$

一階微分の部分と零階項を合わせる。二階部分の反対になる。

$$
\begin{aligned}
-2xP_m^\prime+m(m+1)P_m\big|_{m=4} &= \frac{105 x^{4}}{2} - 60 x^{2} + \frac{15}{2}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
P_m(x)\big|_{x=\frac{1}{2}} &= - \frac{37}{128} \\
\frac{\partial^1 (P_m(x))}{\partial x^1} &= \frac{5 x \left(7 x^{2} - 3\right)}{2} \\
\left.\frac{\partial^1 (P_m(x))}{\partial x^1}\right|_{x=\frac{1}{2}} &= - \frac{25}{16}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
(1-x^2)P_m^{\prime\prime}\big|_{x=\frac{1}{2}} &= \frac{135}{32} \\
\frac{\partial^1 ((1-x^2)P_m^{\prime\prime})}{\partial x^1} &= - 210 x^{3} + 120 x \\
\left.\frac{\partial^1 ((1-x^2)P_m^{\prime\prime})}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{135}{4}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
-2xP_m^\prime+m(m+1)P_m\big|_{x=\frac{1}{2}} &= - \frac{135}{32} \\
\frac{\partial^1 (-2xP_m^\prime+m(m+1)P_m)}{\partial x^1} &= 210 x^{3} - 120 x \\
\left.\frac{\partial^1 (-2xP_m^\prime+m(m+1)P_m)}{\partial x^1}\right|_{x=\frac{1}{2}} &= - \frac{135}{4}
\end{aligned}
$$
:::

### 両端で有界な解はなぜ多項式になるか

$((1-x^2)y')'+\lambda y=0$ の解が $x=\pm1$ で有界だとする。各端は指数 $0$ の重根を持つ確定特異点である。階数低下の式 :eqref[eq-mp-reduction-order] により第 2 解は非零係数の対数を含むので、有界解は各端で正則な枝に属する。従って $(1-x^2)y'$ は両端で零である。方程式に $P_n$ を掛けて積分し、二回部分積分すると
$$
(\lambda-n(n+1))\int_{-1}^1 y(x)P_n(x)\dd x=0
$$
となる。多項式は連続関数の一様近似を通じて $L^2(-1,1)$ に稠密であり、$P_n$ は全多項式を張る。よって $\lambda$ がどの $n(n+1)$ とも違えば $y=0$ となる。$\lambda=n(n+1)$ のときは $P_n$ に対応する係数だけ残るので、非零解は $P_n$ の定数倍である。これで端点条件から固有値を選ぶ向きも示された。

## 母関数と多重極展開

### 逆距離の母関数

:term[母関数]{en="generating function" yomi="ぼかんすう"} は
$$
(1-2xt+t^2)^{-1/2}=\sum_{\ell=0}^\infty P_\ell(x)t^\ell,\qquad |t|<1,\quad |x|\le1.
\label{eq-mp-legendre-generating}
$$
左辺は $t=0$ で $1$ となる枝を取り、複素 $t$ の特異点は単位円上にある。微分して $(1-2xt+t^2)\partial_t G=(x-t)G$ の係数を比較すると $(\ell+1)P_{\ell+1}=(2\ell+1)xP_\ell-\ell P_{\ell-1}$ を得る。左辺の級数が確かに $P_\ell$ となるのは、直交多項式の式 :eqref[eq-mp-orthogonal-recurrence] に偶奇性、先頭係数、ノルムを代入すると同じ漸化式を得るためである。

$\gamma$ を二つの位置ベクトルの間の角とし、$r_>$ と $r_<$ を大きい方と小さい方の半径とする。逆距離は $r_>^{-1}G(\cos\gamma,r_</r_>)$ である。例えば $z=\pm a$ の同符号の点源では、奇数項が打ち消され、外部で核の和は $2/r+2a^2P_2(\cos\theta)/r^3+O(a^4/r^5)$ となる。$a$ と $r$ をメートルで測れば各項の単位は $\mathrm{m^{-1}}$ で同じである。精度評価は $r>a$ を固定し $a/r\to0$ とする。

### 異なる極性の源を重ねる

軸上の源を $z=\pm a$ に置くと、二つの逆距離は角余弦を $x$ と $-x$ にして表せる。$P_m(-x)=(-1)^mP_m(x)$ なので同符号の和では偶数次数だけ、異符号の差では奇数次数だけが残る。この性質は境界の対称性から多重極の最初の非零項を予測するのに使える。展開では $a/r<1$ に固定する。$r=a$ の近くでは幾何級数の残差が小さくなるのに多くの項が要る。

:::example[軸上二点源の成分（1）]
外部の半径を $r$、源の位置を $\pm1$ として係数を比べる。

$$
m=1
$$

正側の源の次数 $m$ の成分を求める。

$$
\begin{aligned}
V_{+,m}\big|_{m=1} &= \frac{x}{r^{2}}
\end{aligned}
$$

負側の源を引いた場合の角度因子には偶奇性が入る。

$$
\begin{aligned}
V_{+,m}-V_{-,m}\big|_{m=1} &= \frac{2 x}{r^{2}}
\end{aligned}
$$

足した場合には逆の偶奇性が残る。

$$
\begin{aligned}
V_{+,m}+V_{-,m}\big|_{m=1} &= 0
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
V_{+,m}\big|_{x=\frac{1}{2}} &= \frac{1}{2 r^{2}} \\
\frac{\partial^1 (V_{+,m})}{\partial x^1} &= \frac{1}{r^{2}} \\
\left.\frac{\partial^1 (V_{+,m})}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{1}{r^{2}}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
V_{+,m}-V_{-,m}\big|_{x=\frac{1}{2}} &= \frac{1}{r^{2}} \\
\frac{\partial^1 (V_{+,m}-V_{-,m})}{\partial x^1} &= \frac{2}{r^{2}} \\
\left.\frac{\partial^1 (V_{+,m}-V_{-,m})}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{2}{r^{2}}
\end{aligned}
$$
:::

:::example[軸上二点源の成分（2）]
外部の半径を $r$、源の位置を $\pm1$ として係数を比べる。次数を一つ上げる。

$$
m=2
$$

正側の源の次数 $m$ の成分を求める。

$$
\begin{aligned}
V_{+,m}\big|_{m=2} &= \frac{3 x^{2} - 1}{2 r^{3}}
\end{aligned}
$$

負側の源を引いた場合の角度因子には偶奇性が入る。

$$
\begin{aligned}
V_{+,m}-V_{-,m}\big|_{m=2} &= 0
\end{aligned}
$$

足した場合には逆の偶奇性が残る。

$$
\begin{aligned}
V_{+,m}+V_{-,m}\big|_{m=2} &= \frac{3 x^{2} - 1}{r^{3}}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
V_{+,m}\big|_{x=\frac{1}{2}} &= - \frac{1}{8 r^{3}} \\
\frac{\partial^1 (V_{+,m})}{\partial x^1} &= \frac{3 x}{r^{3}} \\
\left.\frac{\partial^1 (V_{+,m})}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{3}{2 r^{3}}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
V_{+,m}+V_{-,m}\big|_{x=\frac{1}{2}} &= - \frac{1}{4 r^{3}} \\
\frac{\partial^1 (V_{+,m}+V_{-,m})}{\partial x^1} &= \frac{6 x}{r^{3}} \\
\left.\frac{\partial^1 (V_{+,m}+V_{-,m})}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{3}{r^{3}}
\end{aligned}
$$
:::

:::example[軸上二点源の成分（3）]
外部の半径を $r$、源の位置を $\pm1$ として係数を比べる。さらに次数を上げて、低次の結果と照合する。

$$
m=3
$$

正側の源の次数 $m$ の成分を求める。

$$
\begin{aligned}
V_{+,m}\big|_{m=3} &= \frac{x \left(5 x^{2} - 3\right)}{2 r^{4}}
\end{aligned}
$$

負側の源を引いた場合の角度因子には偶奇性が入る。

$$
\begin{aligned}
V_{+,m}-V_{-,m}\big|_{m=3} &= \frac{x \left(5 x^{2} - 3\right)}{r^{4}}
\end{aligned}
$$

足した場合には逆の偶奇性が残る。

$$
\begin{aligned}
V_{+,m}+V_{-,m}\big|_{m=3} &= 0
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
V_{+,m}\big|_{x=\frac{1}{2}} &= - \frac{7}{16 r^{4}} \\
\frac{\partial^1 (V_{+,m})}{\partial x^1} &= \frac{3 \left(5 x^{2} - 1\right)}{2 r^{4}} \\
\left.\frac{\partial^1 (V_{+,m})}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{3}{8 r^{4}}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
V_{+,m}-V_{-,m}\big|_{x=\frac{1}{2}} &= - \frac{7}{8 r^{4}} \\
\frac{\partial^1 (V_{+,m}-V_{-,m})}{\partial x^1} &= \frac{3 \left(5 x^{2} - 1\right)}{r^{4}} \\
\left.\frac{\partial^1 (V_{+,m}-V_{-,m})}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{3}{4 r^{4}}
\end{aligned}
$$
:::

:::example[軸上二点源の成分（4）]
外部の半径を $r$、源の位置を $\pm1$ として係数を比べる。四次の結果まで具体化する。

$$
m=4
$$

正側の源の次数 $m$ の成分を求める。

$$
\begin{aligned}
V_{+,m}\big|_{m=4} &= \frac{35 x^{4} - 30 x^{2} + 3}{8 r^{5}}
\end{aligned}
$$

負側の源を引いた場合の角度因子には偶奇性が入る。

$$
\begin{aligned}
V_{+,m}-V_{-,m}\big|_{m=4} &= 0
\end{aligned}
$$

足した場合には逆の偶奇性が残る。

$$
\begin{aligned}
V_{+,m}+V_{-,m}\big|_{m=4} &= \frac{35 x^{4} - 30 x^{2} + 3}{4 r^{5}}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
V_{+,m}\big|_{x=\frac{1}{2}} &= - \frac{37}{128 r^{5}} \\
\frac{\partial^1 (V_{+,m})}{\partial x^1} &= \frac{35 x^{3} - 15 x}{2 r^{5}} \\
\left.\frac{\partial^1 (V_{+,m})}{\partial x^1}\right|_{x=\frac{1}{2}} &= - \frac{25}{16 r^{5}}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
V_{+,m}+V_{-,m}\big|_{x=\frac{1}{2}} &= - \frac{37}{64 r^{5}} \\
\frac{\partial^1 (V_{+,m}+V_{-,m})}{\partial x^1} &= \frac{35 x^{3} - 15 x}{r^{5}} \\
\left.\frac{\partial^1 (V_{+,m}+V_{-,m})}{\partial x^1}\right|_{x=\frac{1}{2}} &= - \frac{25}{8 r^{5}}
\end{aligned}
$$
:::

:::example[奇数の多重極だけを選ぶ]
所定の重みで直交射影し、与えた関数のうち低次のモードがどこまで回収されるか調べる。

投影する関数 $f=x^{5}$ を取り、同じ区間と重みで各成分の積分を計算する。正規化済みでない基底も使用できる。

$$
c_j=\frac{\int_{-1}^{1} f\phi_j\,w\,\dd x}{\int_{-1}^{1}\phi_j^2w\,\dd x},\qquad w=1
$$

第 1 の方向に対し、分子と規格化のための分母を別々に計算する。

$$
\begin{aligned}
\phi_{0} &= 1 \\
\int f\phi_{0}w\,\dd x &= 0 \\
\int \phi_{0}^2w\,\dd x &= 2 \\
c_{0} &= 0
\end{aligned}
$$

第 2 の方向に対し、分子と規格化のための分母を別々に計算する。

$$
\begin{aligned}
\phi_{1} &= x \\
\int f\phi_{1}w\,\dd x &= \frac{2}{7} \\
\int \phi_{1}^2w\,\dd x &= \frac{2}{3} \\
c_{1} &= \frac{3}{7}
\end{aligned}
$$

第 3 の方向に対し、分子と規格化のための分母を別々に計算する。

$$
\begin{aligned}
\phi_{2} &= \frac{3 x^{2}}{2} - \frac{1}{2} \\
\int f\phi_{2}w\,\dd x &= 0 \\
\int \phi_{2}^2w\,\dd x &= \frac{2}{5} \\
c_{2} &= 0
\end{aligned}
$$

第 4 の方向に対し、分子と規格化のための分母を別々に計算する。

$$
\begin{aligned}
\phi_{3} &= \frac{5 x^{3}}{2} - \frac{3 x}{2} \\
\int f\phi_{3}w\,\dd x &= \frac{8}{63} \\
\int \phi_{3}^2w\,\dd x &= \frac{2}{7} \\
c_{3} &= \frac{4}{9}
\end{aligned}
$$

第 5 の方向に対し、分子と規格化のための分母を別々に計算する。

$$
\begin{aligned}
\phi_{4} &= \frac{35 x^{4}}{8} - \frac{15 x^{2}}{4} + \frac{3}{8} \\
\int f\phi_{4}w\,\dd x &= 0 \\
\int \phi_{4}^2w\,\dd x &= \frac{2}{9} \\
c_{4} &= 0
\end{aligned}
$$

第 6 の方向に対し、分子と規格化のための分母を別々に計算する。

$$
\begin{aligned}
\phi_{5} &= \frac{63 x^{5}}{8} - \frac{35 x^{3}}{4} + \frac{15 x}{8} \\
\int f\phi_{5}w\,\dd x &= \frac{16}{693} \\
\int \phi_{5}^2w\,\dd x &= \frac{2}{11} \\
c_{5} &= \frac{8}{63}
\end{aligned}
$$

成分を足すと、指定された部分空間への射影が得られる。残差の内積は、使った各基底方向で零になる。

$$
\Pi f=x^{5},\qquad f-\Pi f=0
$$
:::

## 陪 Legendre 関数

### 方位角を持つ解

整数 $0\le m\le\ell$ に対し、区間 $-1<x<1$ の Ferrers 関数を $\mathsf P_\ell^m=(-1)^m(1-x^2)^{m/2}P_\ell^{(m)}$ とする。$(-1)^m$ は Condon–Shortley 位相である。$u=P_\ell^{(m)}$ に対して Legendre 方程式を $m$ 回微分すると $(1-x^2)u''-2(m+1)xu'+[\ell(\ell+1)-m(m+1)]u=0$ となる。$y=(-1)^m(1-x^2)^{m/2}u$ に戻すと
$$
[(1-x^2)y']'+\left[\ell(\ell+1)-\frac{m^2}{1-x^2}\right]y=0
$$
となる。負の次数は $\mathsf P_\ell^{-m}=(-1)^m(\ell-m)!/(\ell+m)!\,\mathsf P_\ell^m$ とする。

ノルムは微分した方程式を $[(1-x^2)^{m+1}u']'+[\ell(\ell+1)-m(m+1)](1-x^2)^mu=0$ と書いて求められる。$u$ を掛けて積分すると、$m+1$ 段のノルムは $m$ 段の $(\ell-m)(\ell+m+1)$ 倍になる。境界項は $(1-x^2)^{m+1}$ で消える。従って
$$
\int_{-1}^1(\mathsf P_\ell^m)^2\dd x=\frac{2(\ell+m)!}{(2\ell+1)(\ell-m)!}.
$$
例えば $\mathsf P_1^1=-\sqrt{1-x^2}$、$\mathsf P_2^1=-3x\sqrt{1-x^2}$ である。符号を変えた流儀では $Y_{\ell m}$ の定義も一緒に調整する必要がある。

:::example[陪関数の方程式を直接確かめる]
$\ell=1,m=1$ の $y=-\sqrt{1-x^2}$ は $y'=x/\sqrt{1-x^2}$ である。従って $[(1-x^2)y']'=(1-2x^2)/\sqrt{1-x^2}$ となる。ポテンシャル項 $[2-1/(1-x^2)]y$ と足せば零である。$y^2=1-x^2$ の積分は $4/3$ で、規格化の式の $\ell=m=1$ の値と一致する。極に近いところで $y'$ は発散するが、$\theta$ で書いた $y=-\sin\theta$ と球面調和の組は滑らかである。
:::

### 陪関数の端点と位相

整数の $m\ge0$ について Ferrers 関数を $(-1)^m(1-x^2)^{m/2}P_\ell^{(m)}(x)$ で作る。$\theta$ が小さいとき $\sqrt{1-x^2}=\sin\theta$ が小さくなる。従って角座標では $\sin^m\theta$ が極の方位角依存を抑える。負の $m$ では規格化に階乗比が入り、次数を同じにしても $m$ と $-m$ の関数は単純な符号だけでは結ばれない。

:::example[二階の陪関数（1）]
$\ell=m+1$ として次数と位数を同時に変え、導関数の次数を比較する。

$$
m=1
$$

微分する前の Legendre 多項式を求める。

$$
\begin{aligned}
P_{m+1}\big|_{m=1} &= \frac{3 x^{2}}{2} - \frac{1}{2}
\end{aligned}
$$

最高次とその次からの寄与を比べる。$m$ 回微分では後者の次数が不足する。

$$
\begin{aligned}
\partial_x^m P_{m+1}\big|_{m=1} &= 3 x
\end{aligned}
$$

Condon–Shortley 位相を掛ける。

$$
\begin{aligned}
\mathsf P_{m+1}^{m}\big|_{m=1} &= - 3 x \sqrt{1 - x^{2}}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
P_{m+1}\big|_{x=\frac{1}{2}} &= - \frac{1}{8} \\
\frac{\partial^1 (P_{m+1})}{\partial x^1} &= 3 x \\
\left.\frac{\partial^1 (P_{m+1})}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{3}{2}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
\partial_x^m P_{m+1}\big|_{x=\frac{1}{2}} &= \frac{3}{2} \\
\frac{\partial^1 (\partial_x^m P_{m+1})}{\partial x^1} &= 3 \\
\left.\frac{\partial^1 (\partial_x^m P_{m+1})}{\partial x^1}\right|_{x=\frac{1}{2}} &= 3
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
\mathsf P_{m+1}^{m}\big|_{x=\frac{1}{2}} &= - \frac{3 \sqrt{3}}{4} \\
\frac{\partial^1 (\mathsf P_{m+1}^{m})}{\partial x^1} &= \frac{3 \left(2 x^{2} - 1\right)}{\sqrt{1 - x^{2}}} \\
\left.\frac{\partial^1 (\mathsf P_{m+1}^{m})}{\partial x^1}\right|_{x=\frac{1}{2}} &= - \sqrt{3}
\end{aligned}
$$
:::

:::example[二階の陪関数（2）]
$\ell=m+1$ として次数と位数を同時に変え、導関数の次数を比較する。次数を一つ上げる。

$$
m=2
$$

微分する前の Legendre 多項式を求める。

$$
\begin{aligned}
P_{m+1}\big|_{m=2} &= \frac{x \left(5 x^{2} - 3\right)}{2}
\end{aligned}
$$

最高次とその次からの寄与を比べる。$m$ 回微分では後者の次数が不足する。

$$
\begin{aligned}
\partial_x^m P_{m+1}\big|_{m=2} &= 15 x
\end{aligned}
$$

Condon–Shortley 位相を掛ける。

$$
\begin{aligned}
\mathsf P_{m+1}^{m}\big|_{m=2} &= 15 x \left(1 - x^{2}\right)
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
P_{m+1}\big|_{x=\frac{1}{2}} &= - \frac{7}{16} \\
\frac{\partial^1 (P_{m+1})}{\partial x^1} &= \frac{15 x^{2}}{2} - \frac{3}{2} \\
\left.\frac{\partial^1 (P_{m+1})}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{3}{8}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
\partial_x^m P_{m+1}\big|_{x=\frac{1}{2}} &= \frac{15}{2} \\
\frac{\partial^1 (\partial_x^m P_{m+1})}{\partial x^1} &= 15 \\
\left.\frac{\partial^1 (\partial_x^m P_{m+1})}{\partial x^1}\right|_{x=\frac{1}{2}} &= 15
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
\mathsf P_{m+1}^{m}\big|_{x=\frac{1}{2}} &= \frac{45}{8} \\
\frac{\partial^1 (\mathsf P_{m+1}^{m})}{\partial x^1} &= 15 - 45 x^{2} \\
\left.\frac{\partial^1 (\mathsf P_{m+1}^{m})}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{15}{4}
\end{aligned}
$$
:::

:::example[二階の陪関数（3）]
$\ell=m+1$ として次数と位数を同時に変え、導関数の次数を比較する。さらに次数を上げて、低次の結果と照合する。

$$
m=3
$$

微分する前の Legendre 多項式を求める。

$$
\begin{aligned}
P_{m+1}\big|_{m=3} &= \frac{35 x^{4}}{8} - \frac{15 x^{2}}{4} + \frac{3}{8}
\end{aligned}
$$

最高次とその次からの寄与を比べる。$m$ 回微分では後者の次数が不足する。

$$
\begin{aligned}
\partial_x^m P_{m+1}\big|_{m=3} &= 105 x
\end{aligned}
$$

Condon–Shortley 位相を掛ける。

$$
\begin{aligned}
\mathsf P_{m+1}^{m}\big|_{m=3} &= - 105 x \left(1 - x^{2}\right)^{\frac{3}{2}}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
P_{m+1}\big|_{x=\frac{1}{2}} &= - \frac{37}{128} \\
\frac{\partial^1 (P_{m+1})}{\partial x^1} &= \frac{5 x \left(7 x^{2} - 3\right)}{2} \\
\left.\frac{\partial^1 (P_{m+1})}{\partial x^1}\right|_{x=\frac{1}{2}} &= - \frac{25}{16}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
\partial_x^m P_{m+1}\big|_{x=\frac{1}{2}} &= \frac{105}{2} \\
\frac{\partial^1 (\partial_x^m P_{m+1})}{\partial x^1} &= 105 \\
\left.\frac{\partial^1 (\partial_x^m P_{m+1})}{\partial x^1}\right|_{x=\frac{1}{2}} &= 105
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
\mathsf P_{m+1}^{m}\big|_{x=\frac{1}{2}} &= - \frac{315 \sqrt{3}}{16} \\
\frac{\partial^1 (\mathsf P_{m+1}^{m})}{\partial x^1} &= \sqrt{1 - x^{2}} \left(420 x^{2} - 105\right) \\
\left.\frac{\partial^1 (\mathsf P_{m+1}^{m})}{\partial x^1}\right|_{x=\frac{1}{2}} &= 0
\end{aligned}
$$
:::

:::example[二階の陪関数（4）]
$\ell=m+1$ として次数と位数を同時に変え、導関数の次数を比較する。四次の結果まで具体化する。

$$
m=4
$$

微分する前の Legendre 多項式を求める。

$$
\begin{aligned}
P_{m+1}\big|_{m=4} &= \frac{x \left(63 x^{4} - 70 x^{2} + 15\right)}{8}
\end{aligned}
$$

最高次とその次からの寄与を比べる。$m$ 回微分では後者の次数が不足する。

$$
\begin{aligned}
\partial_x^m P_{m+1}\big|_{m=4} &= 945 x
\end{aligned}
$$

Condon–Shortley 位相を掛ける。

$$
\begin{aligned}
\mathsf P_{m+1}^{m}\big|_{m=4} &= 945 x \left(x^{2} - 1\right)^{2}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
P_{m+1}\big|_{x=\frac{1}{2}} &= \frac{23}{256} \\
\frac{\partial^1 (P_{m+1})}{\partial x^1} &= \frac{315 x^{4}}{8} - \frac{105 x^{2}}{4} + \frac{15}{8} \\
\left.\frac{\partial^1 (P_{m+1})}{\partial x^1}\right|_{x=\frac{1}{2}} &= - \frac{285}{128}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
\partial_x^m P_{m+1}\big|_{x=\frac{1}{2}} &= \frac{945}{2} \\
\frac{\partial^1 (\partial_x^m P_{m+1})}{\partial x^1} &= 945 \\
\left.\frac{\partial^1 (\partial_x^m P_{m+1})}{\partial x^1}\right|_{x=\frac{1}{2}} &= 945
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
\mathsf P_{m+1}^{m}\big|_{x=\frac{1}{2}} &= \frac{8505}{32} \\
\frac{\partial^1 (\mathsf P_{m+1}^{m})}{\partial x^1} &= 4725 x^{4} - 5670 x^{2} + 945 \\
\left.\frac{\partial^1 (\mathsf P_{m+1}^{m})}{\partial x^1}\right|_{x=\frac{1}{2}} &= - \frac{2835}{16}
\end{aligned}
$$
:::

## 球面調和関数とその完全性

### 球面上の正規直交基底

:::definition[球面調和関数]{#def-mp-spherical-harmonic}
$0\le m\le\ell$ で
$$
Y_{\ell m}(\theta,\phi)
=\sqrt{\frac{2\ell+1}{4\pi}\frac{(\ell-m)!}{(\ell+m)!}}\,
\mathsf P_\ell^m(\cos\theta)e^{im\phi}
$$
を :term[球面調和関数]{en="spherical harmonic" yomi="きゅうめんちょうわかんすう"} とし、$Y_{\ell,-m}=(-1)^m\overline{Y_{\ell m}}$ とする。面積測度は $\dd\Omega=\sin\theta\dd\theta\dd\phi$ である。
:::

Laplacian の角部分に代入して $\Delta_{S^2}Y_{\ell m}=-\ell(\ell+1)Y_{\ell m}$ を得る。異なる $m$ は $\phi$ の積分で直交する。同じ $m$ で異なる $\ell$ は、両方の陪方程式に互いの解を掛けて差を取ると直交する。境界項 $(1-x^2)(y'z-yz')$ は、$y=(1-x^2)^{|m|/2}a$、$z=(1-x^2)^{|m|/2}b$ と書くと $(1-x^2)^{|m|+1}(a'b-ab')$ になる。共通因子の微分が打ち消し合い、$a,b$ は多項式なので両端で零となる。ノルムは前節の積分と方位角の $2\pi$ により $1$ である。

完全性を示すには、すべてに直交する $f\in L^2(S^2)$ の $\phi$ に関する Fourier 係数 $f_m(x)$ を取る。$P_\ell^{(|m|)}$ は $\ell=|m|,|m|+1,\ldots$ で全ての次数を一つずつ持ち、先頭係数が非零である。よって $f_m(x)(1-x^2)^{|m|/2}$ は全多項式に直交する。区間上の多項式は一様近似定理と連続関数の $L^2$ 稠密性から $L^2$ 稠密なので、この積は零である。重みは内部で正だから $f_m=0$ で、$\phi$ の Fourier 級数の完全性から $f=0$ となる。
例えば $Y_{00}=1/\sqrt{4\pi}$、$Y_{10}=\sqrt{3/(4\pi)}\cos\theta$ である。前者は球面平均を、後者は $z$ 方向の双極子成分を抽出する。

### 軸対称の射影と角の平均

球の面積測度は $\sin\theta\,\dd\theta\,\dd\phi$ である。軸対称な関数では $\phi$ を積分して $2\pi$ が出る。従って係数は $[-1,1]$ の単位重みの積分から求まる。相異なる方位角の指数因子の積分は零なので、軸対称関数に $m\ne0$ の成分を混ぜる必要はない。一般の関数では各方位角の Fourier 成分ごとに、この角変数の射影を行う。

:::example[単独の角固有モードの射影（1）]
規格化した軸対称モードと $P_\ell$ の間を行き来する。

$$
m=1
$$

角余弦を $x$ とし、多項式を決める。

$$
\begin{aligned}
P_m(x)\big|_{m=1} &= x
\end{aligned}
$$

球面での正規化には $2\pi$ と多項式のノルムの両方を使う。

$$
\begin{aligned}
Y_{m0}\big|_{m=1} &= \frac{\sqrt{3} x}{2 \sqrt{\pi}}
\end{aligned}
$$

軸対称度分布が $P_m$ のとき、その $Y_{m0}$ 成分の係数を求める。

$$
\begin{aligned}
\int_{S^2}P_m\overline{Y_{m0}}\,\dd\Omega\big|_{m=1} &= \frac{2 \sqrt{3} \sqrt{\pi}}{3}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
P_m(x)\big|_{x=\frac{1}{2}} &= \frac{1}{2} \\
\frac{\partial^1 (P_m(x))}{\partial x^1} &= 1 \\
\left.\frac{\partial^1 (P_m(x))}{\partial x^1}\right|_{x=\frac{1}{2}} &= 1
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
Y_{m0}\big|_{x=\frac{1}{2}} &= \frac{\sqrt{3}}{4 \sqrt{\pi}} \\
\frac{\partial^1 (Y_{m0})}{\partial x^1} &= \frac{\sqrt{3}}{2 \sqrt{\pi}} \\
\left.\frac{\partial^1 (Y_{m0})}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{\sqrt{3}}{2 \sqrt{\pi}}
\end{aligned}
$$
:::

:::example[単独の角固有モードの射影（2）]
規格化した軸対称モードと $P_\ell$ の間を行き来する。次数を一つ上げる。

$$
m=2
$$

角余弦を $x$ とし、多項式を決める。

$$
\begin{aligned}
P_m(x)\big|_{m=2} &= \frac{3 x^{2}}{2} - \frac{1}{2}
\end{aligned}
$$

球面での正規化には $2\pi$ と多項式のノルムの両方を使う。

$$
\begin{aligned}
Y_{m0}\big|_{m=2} &= \frac{\sqrt{5} \left(3 x^{2} - 1\right)}{4 \sqrt{\pi}}
\end{aligned}
$$

軸対称度分布が $P_m$ のとき、その $Y_{m0}$ 成分の係数を求める。

$$
\begin{aligned}
\int_{S^2}P_m\overline{Y_{m0}}\,\dd\Omega\big|_{m=2} &= \frac{2 \sqrt{5} \sqrt{\pi}}{5}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
P_m(x)\big|_{x=\frac{1}{2}} &= - \frac{1}{8} \\
\frac{\partial^1 (P_m(x))}{\partial x^1} &= 3 x \\
\left.\frac{\partial^1 (P_m(x))}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{3}{2}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
Y_{m0}\big|_{x=\frac{1}{2}} &= - \frac{\sqrt{5}}{16 \sqrt{\pi}} \\
\frac{\partial^1 (Y_{m0})}{\partial x^1} &= \frac{3 \sqrt{5} x}{2 \sqrt{\pi}} \\
\left.\frac{\partial^1 (Y_{m0})}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{3 \sqrt{5}}{4 \sqrt{\pi}}
\end{aligned}
$$
:::

:::example[単独の角固有モードの射影（3）]
規格化した軸対称モードと $P_\ell$ の間を行き来する。さらに次数を上げて、低次の結果と照合する。

$$
m=3
$$

角余弦を $x$ とし、多項式を決める。

$$
\begin{aligned}
P_m(x)\big|_{m=3} &= \frac{x \left(5 x^{2} - 3\right)}{2}
\end{aligned}
$$

球面での正規化には $2\pi$ と多項式のノルムの両方を使う。

$$
\begin{aligned}
Y_{m0}\big|_{m=3} &= \frac{\sqrt{7} x \left(5 x^{2} - 3\right)}{4 \sqrt{\pi}}
\end{aligned}
$$

軸対称度分布が $P_m$ のとき、その $Y_{m0}$ 成分の係数を求める。

$$
\begin{aligned}
\int_{S^2}P_m\overline{Y_{m0}}\,\dd\Omega\big|_{m=3} &= \frac{2 \sqrt{7} \sqrt{\pi}}{7}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
P_m(x)\big|_{x=\frac{1}{2}} &= - \frac{7}{16} \\
\frac{\partial^1 (P_m(x))}{\partial x^1} &= \frac{15 x^{2}}{2} - \frac{3}{2} \\
\left.\frac{\partial^1 (P_m(x))}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{3}{8}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
Y_{m0}\big|_{x=\frac{1}{2}} &= - \frac{7 \sqrt{7}}{32 \sqrt{\pi}} \\
\frac{\partial^1 (Y_{m0})}{\partial x^1} &= \frac{3 \sqrt{7} \left(5 x^{2} - 1\right)}{4 \sqrt{\pi}} \\
\left.\frac{\partial^1 (Y_{m0})}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{3 \sqrt{7}}{16 \sqrt{\pi}}
\end{aligned}
$$
:::

:::example[単独の角固有モードの射影（4）]
規格化した軸対称モードと $P_\ell$ の間を行き来する。四次の結果まで具体化する。

$$
m=4
$$

角余弦を $x$ とし、多項式を決める。

$$
\begin{aligned}
P_m(x)\big|_{m=4} &= \frac{35 x^{4}}{8} - \frac{15 x^{2}}{4} + \frac{3}{8}
\end{aligned}
$$

球面での正規化には $2\pi$ と多項式のノルムの両方を使う。

$$
\begin{aligned}
Y_{m0}\big|_{m=4} &= \frac{3 \left(35 x^{4} - 30 x^{2} + 3\right)}{16 \sqrt{\pi}}
\end{aligned}
$$

軸対称度分布が $P_m$ のとき、その $Y_{m0}$ 成分の係数を求める。

$$
\begin{aligned}
\int_{S^2}P_m\overline{Y_{m0}}\,\dd\Omega\big|_{m=4} &= \frac{2 \sqrt{\pi}}{3}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
P_m(x)\big|_{x=\frac{1}{2}} &= - \frac{37}{128} \\
\frac{\partial^1 (P_m(x))}{\partial x^1} &= \frac{5 x \left(7 x^{2} - 3\right)}{2} \\
\left.\frac{\partial^1 (P_m(x))}{\partial x^1}\right|_{x=\frac{1}{2}} &= - \frac{25}{16}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
Y_{m0}\big|_{x=\frac{1}{2}} &= - \frac{111}{256 \sqrt{\pi}} \\
\frac{\partial^1 (Y_{m0})}{\partial x^1} &= \frac{105 x^{3} - 45 x}{4 \sqrt{\pi}} \\
\left.\frac{\partial^1 (Y_{m0})}{\partial x^1}\right|_{x=\frac{1}{2}} &= - \frac{75}{32 \sqrt{\pi}}
\end{aligned}
$$
:::

:::example[四次の軸対称境界値の展開]
所定の重みで直交射影し、与えた関数のうち低次のモードがどこまで回収されるか調べる。

投影する関数 $f=x^{4}$ を取り、同じ区間と重みで各成分の積分を計算する。正規化済みでない基底も使用できる。

$$
c_j=\frac{\int_{-1}^{1} f\phi_j\,w\,\dd x}{\int_{-1}^{1}\phi_j^2w\,\dd x},\qquad w=1
$$

第 1 の方向に対し、分子と規格化のための分母を別々に計算する。

$$
\begin{aligned}
\phi_{0} &= 1 \\
\int f\phi_{0}w\,\dd x &= \frac{2}{5} \\
\int \phi_{0}^2w\,\dd x &= 2 \\
c_{0} &= \frac{1}{5}
\end{aligned}
$$

第 2 の方向に対し、分子と規格化のための分母を別々に計算する。

$$
\begin{aligned}
\phi_{1} &= x \\
\int f\phi_{1}w\,\dd x &= 0 \\
\int \phi_{1}^2w\,\dd x &= \frac{2}{3} \\
c_{1} &= 0
\end{aligned}
$$

第 3 の方向に対し、分子と規格化のための分母を別々に計算する。

$$
\begin{aligned}
\phi_{2} &= \frac{3 x^{2}}{2} - \frac{1}{2} \\
\int f\phi_{2}w\,\dd x &= \frac{8}{35} \\
\int \phi_{2}^2w\,\dd x &= \frac{2}{5} \\
c_{2} &= \frac{4}{7}
\end{aligned}
$$

第 4 の方向に対し、分子と規格化のための分母を別々に計算する。

$$
\begin{aligned}
\phi_{3} &= \frac{5 x^{3}}{2} - \frac{3 x}{2} \\
\int f\phi_{3}w\,\dd x &= 0 \\
\int \phi_{3}^2w\,\dd x &= \frac{2}{7} \\
c_{3} &= 0
\end{aligned}
$$

第 5 の方向に対し、分子と規格化のための分母を別々に計算する。

$$
\begin{aligned}
\phi_{4} &= \frac{35 x^{4}}{8} - \frac{15 x^{2}}{4} + \frac{3}{8} \\
\int f\phi_{4}w\,\dd x &= \frac{16}{315} \\
\int \phi_{4}^2w\,\dd x &= \frac{2}{9} \\
c_{4} &= \frac{8}{35}
\end{aligned}
$$

成分を足すと、指定された部分空間への射影が得られる。残差の内積は、使った各基底方向で零になる。

$$
\Pi f=x^{4},\qquad f-\Pi f=0
$$
:::

### 方位角の次数を上げる

整数 $0\le m<\ell$ について、Ferrers 関数の定義を $x$ で微分する。$x=\cos\theta$ を代入すると
$$
(\partial_\theta-m\cot\theta)\mathsf P_\ell^m(\cos\theta)
 =\mathsf P_\ell^{m+1}(\cos\theta)
$$
となる。実際、$\partial_\theta$ が $\sin^m\theta$ に作用した項を $m\cot\theta$ が打ち消し、$P_\ell$ の導関数に作用した項が $-\sin\theta$ をもう一つ生む。この式を $Y_{\ell m}$ に用い、規格化の比を取ると
$$
 e^{i\phi}(\partial_\theta+i\cot\theta\,\partial_\phi)Y_{\ell m}
 =\sqrt{(\ell-m)(\ell+m+1)}\,Y_{\ell,m+1}
$$
を得る。$m=\ell$ では $P_\ell^{(\ell+1)}=0$ なので左辺が零になる。内積を第 1 引数で線形とし、球面測度で部分積分すると、この作用素の形式随伴は $e^{-i\phi}(-\partial_\theta+i\cot\theta\,\partial_\phi)$ である。$\partial_\theta$ が測度に作用して生じる項は、$\partial_\phi$ が指数因子に作用して生じる項で消える。この随伴と $Y_{\ell,-m}=(-1)^m\overline{Y_{\ell m}}$ を使えば負の $m$ についても同じ上昇係数になる。例えば $Y_{10}$ を一段上げると $-\sqrt{3/(4\pi)}\sin\theta\,e^{i\phi}=\sqrt2\,Y_{11}$ となる。

## 球面調和関数の加法定理

### 頂点間の角だけに依存する核

:::theorem[加法定理]{#thm-mp-addition-theorem}
単位ベクトル $\Omega,\Omega'$ の間の角を $\gamma$ とすると
$$
\sum_{m=-\ell}^\ell Y_{\ell m}(\Omega)\overline{Y_{\ell m}(\Omega')}
=\frac{2\ell+1}{4\pi}P_\ell(\cos\gamma).
\label{eq-mp-spherical-addition}
$$
:::

:::proof
球面の計量を保つ回転は Laplacian と交換する。実際、滑らかな固有関数 $f$ を任意の $Y_{jm}$ と内積して球面上で部分積分すると、$(j(j+1)-\ell(\ell+1))\inner{f}{Y_{jm}}=0$ となる。完全性からその固有空間は次数 $\ell$ の $2\ell+1$ 次元である。回転するとその基底がユニタリ変換されるため、左辺の和は同時回転で不変である。$\Omega'$ を北極へ移すと $|m|>0$ の項は $(1-x^2)^{|m|/2}$ によって零になる。残りは $Y_{\ell0}(\Omega)\overline{Y_{\ell0}(\mathrm{north})}$ であり $P_\ell(1)=1$ を使うと示す式になる。
:::

$\Omega=\Omega'$ で $\sum_m|Y_{\ell m}|^2=(2\ell+1)/(4\pi)$ が分かる。従って Cauchy–Schwarz により $|P_\ell(x)|\le1$（$|x|\le1$）である。母関数に加法定理を入れると各次数が二点の角度を別々に含む形になり、任意の球面電荷分布から多重極係数を積分で取り出せる。

### 射影核と対角値

加法核を角関数の正規直交基底の積和で定める。任意の基底成分と積分するとその成分だけが戻るため、この核は次数 $\ell$ の射影核である。核の対角を球で積分すると、各正規化成分が $1$ を寄与して $2\ell+1$ になる。回転不変なので対角の値は定数であり、球面積 $4\pi$ で割って決まる。非対角の値は二方向間の角度によって決まる。

:::example[方向をずらしたときの加法核（1）]
一方の方向を固定し、二方向の内積を $x$ として計算する。

$$
m=1
$$

角度の依存性は多項式で与えられる。

$$
\begin{aligned}
K_m(x)\big|_{m=1} &= \frac{3 x}{4 \pi}
\end{aligned}
$$

同じ方向を重ねると多項式が $1$ になる。

$$
\begin{aligned}
K_m(1)\big|_{m=1} &= \frac{3}{4 \pi}
\end{aligned}
$$

反対方向では多項式の偶奇性の符号が加わる。

$$
\begin{aligned}
K_m(-1)\big|_{m=1} &= - \frac{3}{4 \pi}
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
K_m(x)\big|_{x=0} &= 0 \\
\frac{\partial^1 (K_m(x))}{\partial x^1} &= \frac{3}{4 \pi} \\
\left.\frac{\partial^1 (K_m(x))}{\partial x^1}\right|_{x=0} &= \frac{3}{4 \pi}
\end{aligned}
$$
:::

:::example[方向をずらしたときの加法核（2）]
一方の方向を固定し、二方向の内積を $x$ として計算する。次数を一つ上げる。

$$
m=2
$$

角度の依存性は多項式で与えられる。

$$
\begin{aligned}
K_m(x)\big|_{m=2} &= \frac{5 \left(3 x^{2} - 1\right)}{8 \pi}
\end{aligned}
$$

同じ方向を重ねると多項式が $1$ になる。

$$
\begin{aligned}
K_m(1)\big|_{m=2} &= \frac{5}{4 \pi}
\end{aligned}
$$

反対方向では多項式の偶奇性の符号が加わる。

$$
\begin{aligned}
K_m(-1)\big|_{m=2} &= \frac{5}{4 \pi}
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
K_m(x)\big|_{x=0} &= - \frac{5}{8 \pi} \\
\frac{\partial^1 (K_m(x))}{\partial x^1} &= \frac{15 x}{4 \pi} \\
\left.\frac{\partial^1 (K_m(x))}{\partial x^1}\right|_{x=0} &= 0
\end{aligned}
$$
:::

:::example[方向をずらしたときの加法核（3）]
一方の方向を固定し、二方向の内積を $x$ として計算する。さらに次数を上げて、低次の結果と照合する。

$$
m=3
$$

角度の依存性は多項式で与えられる。

$$
\begin{aligned}
K_m(x)\big|_{m=3} &= \frac{35 x^{3} - 21 x}{8 \pi}
\end{aligned}
$$

同じ方向を重ねると多項式が $1$ になる。

$$
\begin{aligned}
K_m(1)\big|_{m=3} &= \frac{7}{4 \pi}
\end{aligned}
$$

反対方向では多項式の偶奇性の符号が加わる。

$$
\begin{aligned}
K_m(-1)\big|_{m=3} &= - \frac{7}{4 \pi}
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
K_m(x)\big|_{x=0} &= 0 \\
\frac{\partial^1 (K_m(x))}{\partial x^1} &= \frac{21 \left(5 x^{2} - 1\right)}{8 \pi} \\
\left.\frac{\partial^1 (K_m(x))}{\partial x^1}\right|_{x=0} &= - \frac{21}{8 \pi}
\end{aligned}
$$
:::

:::example[方向をずらしたときの加法核（4）]
一方の方向を固定し、二方向の内積を $x$ として計算する。四次の結果まで具体化する。

$$
m=4
$$

角度の依存性は多項式で与えられる。

$$
\begin{aligned}
K_m(x)\big|_{m=4} &= \frac{9 \left(35 x^{4} - 30 x^{2} + 3\right)}{32 \pi}
\end{aligned}
$$

同じ方向を重ねると多項式が $1$ になる。

$$
\begin{aligned}
K_m(1)\big|_{m=4} &= \frac{9}{4 \pi}
\end{aligned}
$$

反対方向では多項式の偶奇性の符号が加わる。

$$
\begin{aligned}
K_m(-1)\big|_{m=4} &= \frac{9}{4 \pi}
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
K_m(x)\big|_{x=0} &= \frac{27}{32 \pi} \\
\frac{\partial^1 (K_m(x))}{\partial x^1} &= \frac{315 x^{3} - 135 x}{8 \pi} \\
\left.\frac{\partial^1 (K_m(x))}{\partial x^1}\right|_{x=0} &= 0
\end{aligned}
$$
:::

## 第 2 種 Legendre 関数

### 第 2 解と端点の対数

$x>1$ の第 2 種を $Q_\ell(x)=\frac12\int_{-1}^1 P_\ell(t)/(x-t)\dd t$ と定める。この核について $\partial_x[(1-x^2)\partial_x](x-t)^{-1}=\partial_t[(1-t^2)\partial_t](x-t)^{-1}$ である。右辺を積分で二回部分積分すると端点で $1-t^2$ が消え、Legendre 方程式から $Q_\ell$ も同じ方程式を満たす。$x\to+\infty$ では $1/(x-t)$ を幾何級数に展開し、直交性により初項は $x^{-\ell-1}\int P_\ell t^\ell/2$ で非零である。多項式解とは独立である。$-1<x<1$ の実の第 2 種 $\mathsf Q_\ell$ は同じ積分の主値で得られる。境界値の実部分と読めるのは式 :eqref[eq-mp-plemelj] による。
例えば $Q_0(x)=\frac12\log((x+1)/(x-1))$、$Q_1=xQ_0-1$ である。区間内の主値では $\mathsf Q_0(x)=\frac12\log((1+x)/(1-x))$ になる。この端点の対数を含む解は、球面全体で滑らかな固有関数からは除かれる。

### 外部領域での積分表示

実数 $z>1$ に対し $Q_\ell(z)=\frac12\int_{-1}^1 P_\ell(t)/(z-t)\dd t$ と表せる。$\ell=0$ では積分が対数を与える。三項漸化式を積分すると $t/(z-t)=z/(z-t)-1$ により $\ell=0$ 以外では定数との直交性から余分な項が消える。$\ell=0$ の余分な項が $Q_1=zQ_0-1$ を与える。従って漸化式によりすべての整数次数で同じ規格化になる。

:::example[第二種関数の先頭の遠方項（1）]
$1/(z-t)$ を $z^{-1}\sum_{j\ge0}t^j z^{-j}$ に展開する。

$$
m=1
$$

最初の非零モーメントの相手は次数が同じ冪である。

$$
\begin{aligned}
\frac12\int_{-1}^1P_m(t)t^m\dd t\big|_{m=1} &= \frac{1}{3}
\end{aligned}
$$

分母の展開でこのモーメントに $z^{-m-1}$ が掛かる。

$$
\begin{aligned}
Q_m(z)\ \text{の主項}\big|_{m=1} &= \frac{1}{3 u^{2}}
\end{aligned}
$$

次の奇偶が逆のモーメントは消えるため、指標は一つ飛ぶ。

$$
\begin{aligned}
\frac12\int_{-1}^1P_m(t)t^{m+1}\dd t &= 0
\end{aligned}
$$

$u=3$ で値と導関数を評価する。

$$
\begin{aligned}
Q_m(z)\ \text{の主項}\big|_{u=3} &= \frac{1}{27} \\
\frac{\partial^1 (Q_m(z)\ \text{の主項})}{\partial u^1} &= - \frac{2}{3 u^{3}} \\
\left.\frac{\partial^1 (Q_m(z)\ \text{の主項})}{\partial u^1}\right|_{u=3} &= - \frac{2}{81}
\end{aligned}
$$
:::

:::example[第二種関数の先頭の遠方項（2）]
$1/(z-t)$ を $z^{-1}\sum_{j\ge0}t^j z^{-j}$ に展開する。次数を一つ上げる。

$$
m=2
$$

最初の非零モーメントの相手は次数が同じ冪である。

$$
\begin{aligned}
\frac12\int_{-1}^1P_m(t)t^m\dd t\big|_{m=2} &= \frac{2}{15}
\end{aligned}
$$

分母の展開でこのモーメントに $z^{-m-1}$ が掛かる。

$$
\begin{aligned}
Q_m(z)\ \text{の主項}\big|_{m=2} &= \frac{2}{15 u^{3}}
\end{aligned}
$$

次の奇偶が逆のモーメントは消えるため、指標は一つ飛ぶ。

$$
\begin{aligned}
\frac12\int_{-1}^1P_m(t)t^{m+1}\dd t &= 0
\end{aligned}
$$

$u=3$ で値と導関数を評価する。

$$
\begin{aligned}
Q_m(z)\ \text{の主項}\big|_{u=3} &= \frac{2}{405} \\
\frac{\partial^1 (Q_m(z)\ \text{の主項})}{\partial u^1} &= - \frac{2}{5 u^{4}} \\
\left.\frac{\partial^1 (Q_m(z)\ \text{の主項})}{\partial u^1}\right|_{u=3} &= - \frac{2}{405}
\end{aligned}
$$
:::

:::example[第二種関数の先頭の遠方項（3）]
$1/(z-t)$ を $z^{-1}\sum_{j\ge0}t^j z^{-j}$ に展開する。さらに次数を上げて、低次の結果と照合する。

$$
m=3
$$

最初の非零モーメントの相手は次数が同じ冪である。

$$
\begin{aligned}
\frac12\int_{-1}^1P_m(t)t^m\dd t\big|_{m=3} &= \frac{2}{35}
\end{aligned}
$$

分母の展開でこのモーメントに $z^{-m-1}$ が掛かる。

$$
\begin{aligned}
Q_m(z)\ \text{の主項}\big|_{m=3} &= \frac{2}{35 u^{4}}
\end{aligned}
$$

次の奇偶が逆のモーメントは消えるため、指標は一つ飛ぶ。

$$
\begin{aligned}
\frac12\int_{-1}^1P_m(t)t^{m+1}\dd t &= 0
\end{aligned}
$$

$u=3$ で値と導関数を評価する。

$$
\begin{aligned}
Q_m(z)\ \text{の主項}\big|_{u=3} &= \frac{2}{2835} \\
\frac{\partial^1 (Q_m(z)\ \text{の主項})}{\partial u^1} &= - \frac{8}{35 u^{5}} \\
\left.\frac{\partial^1 (Q_m(z)\ \text{の主項})}{\partial u^1}\right|_{u=3} &= - \frac{8}{8505}
\end{aligned}
$$
:::

:::example[第二種関数の先頭の遠方項（4）]
$1/(z-t)$ を $z^{-1}\sum_{j\ge0}t^j z^{-j}$ に展開する。四次の結果まで具体化する。

$$
m=4
$$

最初の非零モーメントの相手は次数が同じ冪である。

$$
\begin{aligned}
\frac12\int_{-1}^1P_m(t)t^m\dd t\big|_{m=4} &= \frac{8}{315}
\end{aligned}
$$

分母の展開でこのモーメントに $z^{-m-1}$ が掛かる。

$$
\begin{aligned}
Q_m(z)\ \text{の主項}\big|_{m=4} &= \frac{8}{315 u^{5}}
\end{aligned}
$$

次の奇偶が逆のモーメントは消えるため、指標は一つ飛ぶ。

$$
\begin{aligned}
\frac12\int_{-1}^1P_m(t)t^{m+1}\dd t &= 0
\end{aligned}
$$

$u=3$ で値と導関数を評価する。

$$
\begin{aligned}
Q_m(z)\ \text{の主項}\big|_{u=3} &= \frac{8}{76545} \\
\frac{\partial^1 (Q_m(z)\ \text{の主項})}{\partial u^1} &= - \frac{8}{63 u^{6}} \\
\left.\frac{\partial^1 (Q_m(z)\ \text{の主項})}{\partial u^1}\right|_{u=3} &= - \frac{8}{45927}
\end{aligned}
$$
:::

### 第 2 種の規格化と境界値

積分表示を複素 $z\in\C\setminus[-1,1]$ に延ばす。$-1<x<1$ の上側と下側の極限は、式 :eqref[eq-mp-plemelj] を分母 $x-t\pm i0$ に適用して
$$
Q_\ell(x\pm i0)=\mathsf Q_\ell(x)\mp\frac{i\pi}{2}P_\ell(x)
$$
である。従って区間内の実関数 $\mathsf Q_\ell$ は境界値の平均である。例えば $Q_1$ なら平均は $x\mathsf Q_0(x)-1$ となる。ここでは整数次数・位数 0 の通常の $Q_\ell$ を用いている。

また $P_\ell(t)=P_\ell(x)+(t-x)S_{\ell-1}(t,x)$ と割り算すると
$$
\mathsf Q_\ell(x)=P_\ell(x)\mathsf Q_0(x)
 -\frac12\int_{-1}^1 S_{\ell-1}(t,x)\dd t
$$
であり、右端の積分は $x$ の多項式である。これにより $x\uparrow1$ で $P_\ell(1)=1$ が対数の係数を決める。さらに Abel の式 :eqref[eq-mp-abel-wronskian] により、区間内で $P_\ell\mathsf Q_\ell'-P_\ell'\mathsf Q_\ell=C/(1-x^2)$ となる。右端に近づくと多項式項は $(1-x^2)$ を掛けた積から消え、$(1-x^2)\mathsf Q_0'=1$ だけが残るので $C=1$ である。

## 球座標での Laplace 方程式の解

### 角度の固有値が動径の冪を決める

$r>0$ の球座標で $u=R(r)Y_{\ell m}$ とおくと $(r^2R')'=\ell(\ell+1)R$ である。$R=r^\alpha$ を代入して $(\alpha-\ell)(\alpha+\ell+1)=0$ より二解 $r^\ell,r^{-\ell-1}$ を得る。球の内側で原点に正則なら前者、無限遠で零となる外部解なら後者を選ぶ。境界 $r=a$ の係数が $f_{\ell m}=\int f\overline{Y_{\ell m}}\dd\Omega$ なら、内側は $\sum f_{\ell m}(r/a)^\ell Y_{\ell m}$、外側は $\sum f_{\ell m}(a/r)^{\ell+1}Y_{\ell m}$ である。境界から離れると、幾何減衰と式 :eqref[eq-mp-spherical-addition] の対角値から級数は局所一様収束する。調和な部分和の局所一様極限は平均値性によって調和である。境界では $L^2$ において各係数の倍率が $1$ に近づくため $f$ に収束する。連続な境界値については最大値原理で連続な解の一意性を示せる。外部での一意性には無限遠で一様に零へ近づく条件を用い、有限の外半径で囲んだ球殻に最大値原理を適用してから外半径を無限大にする。

:::example[四重極型の境界電位]
半径 $a=2\,\mathrm m$ の球面に $V_0P_2(\cos\theta)$ を与える。内側では $V_0(r/a)^2P_2$、外側では $V_0(a/r)^3P_2$ である。$V_0=8\,\mathrm V$ なら外側の軸上 $r=4\,\mathrm m$ で $1\,\mathrm V$、赤道上の同じ半径で $-1/2\,\mathrm V$ になる。
:::

### 球の Poisson 核

内側の解の各次数を式 :eqref[eq-mp-spherical-addition] で境界積分へ戻す。$t=r/a<1$ に対し、角距離の余弦を $\mu=\Omega\cdot\Omega'$ とすると核は
$$
\frac1{4\pi}\sum_{\ell\ge0}(2\ell+1)t^\ell P_\ell(\mu)
=\frac{1-t^2}{4\pi(1-2t\mu+t^2)^{3/2}}
$$
である。母関数に $1+2t\partial_t$ を作用させれば右辺を得る。積分すると次数零だけ残るので全質量は $1$、核は非負である。$t\uparrow1$ に $\Omega$ の任意の角近傍の外で分母は零から離れているから、外の質量は零へ近づく。従って連続な境界値に対し一様収束で境界まで延長できる。例えば $r=0$ では核は $1/(4\pi)$、内側の中心電位は境界の球面平均である。

### 内側と外側での境界値

外部解は内側の解の半径を反転して $u_{\mathrm{out}}(r,\Omega)=(a/r)\,u_{\mathrm{in}}(a^2/r,\Omega)$ とも書ける。実際、各モードの $(r/a)^\ell$ に代入すると $(a/r)^{\ell+1}$ になるので調和性が保たれる。従って内部の Poisson 核で示した境界への一様収束は、外部の連続解にも移る。この変換は球面を固定する :term[ケルビン変換]{en="Kelvin transform" yomi="けるびんへんかん"} である。球の境界で各モードの値を $g_{\ell m}$ とする。内側で原点に滑らかな解は $(r/a)^\ell$ が掛かる。外側で無限遠で零になる解では $(a/r)^{\ell+1}$ が掛かる。両者の半径方向の微分の差は球面上で $-(2\ell+1)g_{\ell m}/a$ となる。そのため同じ角モードで球表面に単位面積当たりの電荷を載せると、電位と表面電荷の係数を直接結び付けられる。

:::example[球面に載せた電荷の角モード（1）]
半径 $2$ の球で境界電位を角モード $P_m$ 一つにした場合を考える。

$$
m=1
$$

球の内側の動径因子を境界で規格化する。

$$
\begin{aligned}
V_{\rm in}\big|_{m=1} &= \frac{r x}{2}
\end{aligned}
$$

外側は動径の負の冪を選ぶ。

$$
\begin{aligned}
V_{\rm out}\big|_{m=1} &= \frac{4 x}{r^{2}}
\end{aligned}
$$

境界で半径方向の微分の差を求める。外向きに測った電場の跳びはこの符号を反転したものである。

$$
\begin{aligned}
\partial_rV_{\rm out}-\partial_rV_{\rm in}\big|_{m=1} &= - \frac{3 x}{2}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
V_{\rm in}\big|_{x=\frac{1}{2}} &= \frac{r}{4} \\
\frac{\partial^1 (V_{\rm in})}{\partial x^1} &= \frac{r}{2} \\
\left.\frac{\partial^1 (V_{\rm in})}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{r}{2}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
V_{\rm out}\big|_{x=\frac{1}{2}} &= \frac{2}{r^{2}} \\
\frac{\partial^1 (V_{\rm out})}{\partial x^1} &= \frac{4}{r^{2}} \\
\left.\frac{\partial^1 (V_{\rm out})}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{4}{r^{2}}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
\partial_rV_{\rm out}-\partial_rV_{\rm in}\big|_{x=\frac{1}{2}} &= - \frac{3}{4} \\
\frac{\partial^1 (\partial_rV_{\rm out}-\partial_rV_{\rm in})}{\partial x^1} &= - \frac{3}{2} \\
\left.\frac{\partial^1 (\partial_rV_{\rm out}-\partial_rV_{\rm in})}{\partial x^1}\right|_{x=\frac{1}{2}} &= - \frac{3}{2}
\end{aligned}
$$
:::

:::example[球面に載せた電荷の角モード（2）]
半径 $2$ の球で境界電位を角モード $P_m$ 一つにした場合を考える。次数を一つ上げる。

$$
m=2
$$

球の内側の動径因子を境界で規格化する。

$$
\begin{aligned}
V_{\rm in}\big|_{m=2} &= \frac{r^{2} \left(3 x^{2} - 1\right)}{8}
\end{aligned}
$$

外側は動径の負の冪を選ぶ。

$$
\begin{aligned}
V_{\rm out}\big|_{m=2} &= \frac{4 \left(3 x^{2} - 1\right)}{r^{3}}
\end{aligned}
$$

境界で半径方向の微分の差を求める。外向きに測った電場の跳びはこの符号を反転したものである。

$$
\begin{aligned}
\partial_rV_{\rm out}-\partial_rV_{\rm in}\big|_{m=2} &= \frac{5}{4} - \frac{15 x^{2}}{4}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
V_{\rm in}\big|_{x=\frac{1}{2}} &= - \frac{r^{2}}{32} \\
\frac{\partial^1 (V_{\rm in})}{\partial x^1} &= \frac{3 r^{2} x}{4} \\
\left.\frac{\partial^1 (V_{\rm in})}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{3 r^{2}}{8}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
V_{\rm out}\big|_{x=\frac{1}{2}} &= - \frac{1}{r^{3}} \\
\frac{\partial^1 (V_{\rm out})}{\partial x^1} &= \frac{24 x}{r^{3}} \\
\left.\frac{\partial^1 (V_{\rm out})}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{12}{r^{3}}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
\partial_rV_{\rm out}-\partial_rV_{\rm in}\big|_{x=\frac{1}{2}} &= \frac{5}{16} \\
\frac{\partial^1 (\partial_rV_{\rm out}-\partial_rV_{\rm in})}{\partial x^1} &= - \frac{15 x}{2} \\
\left.\frac{\partial^1 (\partial_rV_{\rm out}-\partial_rV_{\rm in})}{\partial x^1}\right|_{x=\frac{1}{2}} &= - \frac{15}{4}
\end{aligned}
$$
:::

:::example[球面に載せた電荷の角モード（3）]
半径 $2$ の球で境界電位を角モード $P_m$ 一つにした場合を考える。さらに次数を上げて、低次の結果と照合する。

$$
m=3
$$

球の内側の動径因子を境界で規格化する。

$$
\begin{aligned}
V_{\rm in}\big|_{m=3} &= \frac{r^{3} x \left(5 x^{2} - 3\right)}{16}
\end{aligned}
$$

外側は動径の負の冪を選ぶ。

$$
\begin{aligned}
V_{\rm out}\big|_{m=3} &= \frac{40 x^{3} - 24 x}{r^{4}}
\end{aligned}
$$

境界で半径方向の微分の差を求める。外向きに測った電場の跳びはこの符号を反転したものである。

$$
\begin{aligned}
\partial_rV_{\rm out}-\partial_rV_{\rm in}\big|_{m=3} &= \frac{7 x \left(3 - 5 x^{2}\right)}{4}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
V_{\rm in}\big|_{x=\frac{1}{2}} &= - \frac{7 r^{3}}{128} \\
\frac{\partial^1 (V_{\rm in})}{\partial x^1} &= \frac{3 r^{3} \left(5 x^{2} - 1\right)}{16} \\
\left.\frac{\partial^1 (V_{\rm in})}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{3 r^{3}}{64}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
V_{\rm out}\big|_{x=\frac{1}{2}} &= - \frac{7}{r^{4}} \\
\frac{\partial^1 (V_{\rm out})}{\partial x^1} &= \frac{24 \left(5 x^{2} - 1\right)}{r^{4}} \\
\left.\frac{\partial^1 (V_{\rm out})}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{6}{r^{4}}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
\partial_rV_{\rm out}-\partial_rV_{\rm in}\big|_{x=\frac{1}{2}} &= \frac{49}{32} \\
\frac{\partial^1 (\partial_rV_{\rm out}-\partial_rV_{\rm in})}{\partial x^1} &= \frac{21}{4} - \frac{105 x^{2}}{4} \\
\left.\frac{\partial^1 (\partial_rV_{\rm out}-\partial_rV_{\rm in})}{\partial x^1}\right|_{x=\frac{1}{2}} &= - \frac{21}{16}
\end{aligned}
$$
:::

:::example[球面に載せた電荷の角モード（4）]
半径 $2$ の球で境界電位を角モード $P_m$ 一つにした場合を考える。四次の結果まで具体化する。

$$
m=4
$$

球の内側の動径因子を境界で規格化する。

$$
\begin{aligned}
V_{\rm in}\big|_{m=4} &= \frac{r^{4} \left(35 x^{4} - 30 x^{2} + 3\right)}{128}
\end{aligned}
$$

外側は動径の負の冪を選ぶ。

$$
\begin{aligned}
V_{\rm out}\big|_{m=4} &= \frac{4 \left(35 x^{4} - 30 x^{2} + 3\right)}{r^{5}}
\end{aligned}
$$

境界で半径方向の微分の差を求める。外向きに測った電場の跳びはこの符号を反転したものである。

$$
\begin{aligned}
\partial_rV_{\rm out}-\partial_rV_{\rm in}\big|_{m=4} &= - \frac{315 x^{4}}{16} + \frac{135 x^{2}}{8} - \frac{27}{16}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
V_{\rm in}\big|_{x=\frac{1}{2}} &= - \frac{37 r^{4}}{2048} \\
\frac{\partial^1 (V_{\rm in})}{\partial x^1} &= \frac{5 r^{4} x \left(7 x^{2} - 3\right)}{32} \\
\left.\frac{\partial^1 (V_{\rm in})}{\partial x^1}\right|_{x=\frac{1}{2}} &= - \frac{25 r^{4}}{256}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
V_{\rm out}\big|_{x=\frac{1}{2}} &= - \frac{37}{4 r^{5}} \\
\frac{\partial^1 (V_{\rm out})}{\partial x^1} &= \frac{560 x^{3} - 240 x}{r^{5}} \\
\left.\frac{\partial^1 (V_{\rm out})}{\partial x^1}\right|_{x=\frac{1}{2}} &= - \frac{50}{r^{5}}
\end{aligned}
$$

$x=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
\partial_rV_{\rm out}-\partial_rV_{\rm in}\big|_{x=\frac{1}{2}} &= \frac{333}{256} \\
\frac{\partial^1 (\partial_rV_{\rm out}-\partial_rV_{\rm in})}{\partial x^1} &= \frac{45 x \left(3 - 7 x^{2}\right)}{4} \\
\left.\frac{\partial^1 (\partial_rV_{\rm out}-\partial_rV_{\rm in})}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{225}{32}
\end{aligned}
$$
:::

## 演習問題

::::exercise{level=basic}
$P_3$ を漸化式から求めよ。
:::solution
$3P_3=5xP_2-2P_1$ より $(5x^3-3x)/2$ である。
:::
::::

::::exercise{level=basic}
$\mathsf P_2^2(x)$ を求めよ。
:::solution
$P_2''=3$ で位相は正、従って $3(1-x^2)$ である。
:::
::::

::::exercise{level=basic}
$\int_{-1}^1[P_2(x)]^2\dd x$ を直接計算せよ。
:::solution
$\tfrac14(9x^4-6x^2+1)$ を積分して $\tfrac14(18/5-4+2)=2/5$ である。
:::
::::

::::exercise{level=basic}
$Y_{11}$ の $\phi$ 依存性と北極での値を求めよ。
:::solution
$\mathsf P_1^1(\cos\theta)=-\sin\theta$ なので $Y_{11}=-\sqrt{3/(8\pi)}\sin\theta\,e^{i\phi}$ である。北極で零になる。
:::
::::

::::exercise{level=advanced}
$Q_0,Q_1$ の先頭の無限遠漸近形を求めよ。
:::solution
$\tfrac12\log((1+1/x)/(1-1/x))=x^{-1}+\tfrac13 x^{-3}+O(x^{-5})$ である。従って $Q_0\sim x^{-1}$、$Q_1=xQ_0-1\sim(3x^2)^{-1}$ である。
:::
::::

::::exercise{level=advanced}
$r=a$ で $V=3\cos\theta+2$ のとき原点に有限な内側解を求めよ。
:::solution
$\cos\theta=P_1$、定数は $P_0$ であり、それぞれの動径倍率は $r/a,1$ である。よって $V=2+3(r/a)\cos\theta=2+3z/a$ となり直交座標でも Laplacian が零である。
:::
::::

::::exercise{level=advanced}
直交する二方向で次数 $\ell=2$ の加法核の値を求めよ。
:::solution
$\cos\gamma=0$ で $P_2(0)=-1/2$ だから :eqref[eq-mp-spherical-addition] より $-5/(8\pi)$ である。
:::
::::

::::exercise{level=advanced}
真空中で $z=a$ に $+q$、$z=-a$ に $-q$ の点電荷を置く。$r>a$ における電位の $a/r\to0$ での主項を求めよ。
:::solution
逆距離の差では $P_\ell(-x)=(-1)^\ell P_\ell(x)$ により偶数項が消える。従って $\ell=1$ が最初で、電荷を含む核は $2qa\cos\theta/r^2$ である。次は $\ell=3$ で $O(qa^3/r^4)$ である。真空の電位にはさらに $1/(4\pi\epsilon_0)$ が掛かる。
:::
::::

## 確認テスト

::::quiz{answer="7" tol="1e-9"}
$\ell=3$ の角固有空間の次元を求めよ。

:::explanation
$m=-3,-2,-1,0,1,2,3$ が正規直交基底を与えるため $7$ である。
:::
::::

::::quiz
$\mathsf P_1^1(x)$ の符号は Condon–Shortley 位相を含めるとどれか。

- [x] $-\sqrt{1-x^2}$
- [ ] $+\sqrt{1-x^2}$
- [ ] $x$

:::explanation
$P_1'=1$ に $(1-x^2)^{1/2}$ と $(-1)^1$ を掛ける。
:::
::::

::::quiz{answer="-0.5" tol="1e-9"}
$P_2(0)$ を答えよ。

:::explanation
$P_2(x)=(3x^2-1)/2$ への代入から $-1/2$ である。
:::
::::

::::quiz
$\ell=2$ の外部で無限遠に減衰する動径因子はどれか。

- [ ] $r^2$
- [ ] $r^{-2}$
- [x] $r^{-3}$

:::explanation
$\alpha(\alpha+1)=6$ の根は $2,-3$ である。負の根だけが無限遠で減衰する。
:::
::::

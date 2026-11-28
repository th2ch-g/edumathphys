---
title: "緩増加超関数と Fourier 変換"
description: "Schwartz 空間の双対に変換を拡張し、主値、Sobolev 空間、構造定理を扱う。"
notation:
  - symbol: '\mathcal S'''
    meaning: "緩増加超関数の空間"
  - symbol: '\langle\xi\rangle'
    meaning: "(1+4π²|ξ|²) の平方根"
  - symbol: 'H^s'
    meaning: "Fourier 変換で定義する Sobolev 空間"
  - symbol: '\operatorname{pv}(1/x)'
    meaning: "1/x の対称主値"
---

## この章の目標

- 緩増加性を Schwartz セミノルムから判定する。
- 双対性で Fourier 変換を定義し、デルタ、主値、斉次核を変換する。
- $H^s$ の双対、微分、埋め込みを周波数側で証明する。
- 緩増加超関数を連続関数の導関数で表す。

## 緩増加超関数

試験関数をコンパクト台から急減少関数へ広げると、無限遠における増大も制御しなければならない。局所的には超関数でも、急減少試験関数への作用が連続に延長できない例がある。

:::definition[緩増加超関数]{#def-four-tempered}
:ref[def-four-schwartz] の空間 $\mathcal S(\R^n)$ の連続な複素線形双対を $\mathcal S'(\R^n)$ と書き、その元を :term[緩増加超関数]{en="tempered distribution" yomi="かんぞうかちょうかんすう"} という。$\langle T,\phi\rangle$ に共役は含めない。$T\in\mathcal S'$ であるための連続性は

$$
|\langle T,\phi\rangle|\le C\max_{|\alpha|\le m}\sup_x(1+|x|)^N|\partial^\alpha\phi(x)|
\label{eq-four-tempered-bound}
$$

をある整数 $m,N\ge0$ と $C$ で満たすことと同値である。
:::

右辺の半ノルムは有限個の Schwartz セミノルムの和で抑えられ、逆に各セミノルムは右辺の形で抑えられる。線形汎関数の連続性は零近傍での有界性を意味するので、斉次性により評価と同値となる。$\mathcal D$ は $\mathcal S$ に稠密である。原点近くで $1$ の固定切断を $\chi(x/R)$ に拡大して掛ければ、差の各導関数は $|x|\gtrsim R$ で少なくとも一つ余分な減衰冪で抑えられる。その結果、制限 $\mathcal S'\to\mathcal D'$ は単射である。

:::example[多項式増大と指数増大]
$|f(x)|\le A(1+|x|)^q$ を満たす局所可積分関数の作用 $\int f\phi$ は、$N>q+n$ と取れば :eqref[eq-four-tempered-bound] で抑えられる。$\int (1+|x|)^{-N}d|\mu|<\infty$ となる Radon 測度も同様である。対して $e^x$ は $\R$ 上の緩増加超関数にならない。固定された非負で非零の $\psi\in\mathcal D$ を $R$ だけ右へ移すと、各 Schwartz セミノルムは多項式増大にとどまるが、作用は $e^R\int e^x\psi(x)dx$ である。
:::

:::example[無限個の点質量]
$\sum_{k\in\Z}\delta_k$ の作用は $\sum_k\phi(k)$ である。$|\phi(k)|\le C(1+|k|)^{-2}$ で和が収束し、その $C$ はセミノルムで取れる。任意のコンパクト台超関数も :ref[prop-four-compact-distribution] によって $\mathcal S$ 上に連続に作用する。
:::

### 測度の質量と緩増加性

点ごとに多項式で抑えられる関数は緩増加超関数を与える。しかし判定に直接必要なのは、試験関数が与える重み付きの積分である。特に正測度には相殺がないため、半径を大きくしたときの質量の評価から必要十分条件を得られる。

:::proposition[正測度の緩増加性]{#prop-four-tempered-positive-measure}
局所有限な正 Radon 測度 $\mu$ について、次は同値である。

1. $\mathcal D$ 上の測度の作用が $\mathcal S'$ に延長される。
2. ある $C,q\ge0$ について $\mu(B(0,R))\le CR^q$ が $R\ge1$ で成り立つ。
3. ある $N\ge0$ について $\int(1+|x|)^{-N}d\mu(x)<\infty$ である。
:::

:::proof
非負の試験関数 $\chi$ で単位球上で $1$、半径 $2$ の球の外でゼロになるものを取る。条件 1 の評価 :eqref[eq-four-tempered-bound] を $\chi_R(x)=\chi(x/R)$ に適用する。$R\ge1$ で $\partial^\alpha\chi_R$ は $R^{-|\alpha|}$ 倍の有界な導関数であり、台上の重みは $(1+2R)^N$ 以下である。$\mu$ が正であるため、$\mu(B(0,R))\le\langle\mu,\chi_R\rangle\le C'R^N$ である。これが条件 2 を与える。

条件 2 からは、$2^j\le|x|<2^{j+1}$ における重み付き積分を $C2^{(j+1)q-jN}$ で抑えられる。$N>q$ を選べばこの幾何級数が収束する。単位球上では局所有限性を使うので条件 3 となる。最後に条件 3 の積分を定数に取り、試験関数を $(1+|x|)^N|\phi(x)|$ の上限で評価すれば条件 1 である。
:::

:::example[高さは指数増大でも質量は有限]
$\R$ 上で互いに素な区間 $I_j=[j,j+2^{-2j}]$（$j\ge1$）を取り、$f=\sum_{j\ge1}2^j1_{I_j}$ とする。$f$ は各区間で急増する高さを持つが、積分は $\sum_j2^j2^{-2j}=\sum_j2^{-j}=1$ である。従って $f$ は $L^1$ を経由して緩増加超関数となる。各 $I_j$ は正の測度を持ち、そこで $2^j$ を $j$ の固定の冪で上から抑えることはできない。このため点ごとの多項式増大は必要条件ではない。
:::

:::example[符号の相殺による緩増加性]
$f(x)=e^x\cos(e^x)$ とする。通常微分で $f=(\sin(e^x))'$ であり、$\sin(e^x)$ は有界である。従って試験関数に対して
$$
\langle f,\phi\rangle=-\int \sin(e^x)\phi'(x)dx
$$
と作用を延長できる。右辺は $N>1$ の Schwartz セミノルムで抑えられる。一方 $|f|\,dx$ は緩増加でない。実際、$R>0$ に対し $\int_0^R |f(x)|dx=\int_1^{e^R}|\cos t|dt$ であり、周期 $\pi$ ごとに積分が $2$ なので $e^R$ と同程度に増える。正測度なら必要だった質量条件が、符号付きの作用では相殺のため不要になる。
:::

### Schwartz 空間での弱収束と局所収束

緩増加超関数の列 $T_j$ が弱く $T$ に収束するとは、各 $\phi\in\mathcal S$ に対して $\langle T_j,\phi\rangle\to\langle T,\phi\rangle$ となることをいう。通常の超関数としての収束は、ここで使う試験関数を $\mathcal D$ に制限したものである。両者の違いは、無限遠へ移った台を Schwartz 関数が検出できることにある。

:::proposition[共通の緩増加評価による弱収束]{#prop-four-tempered-uniform-bound}
$T_j\in\mathcal S'$ が全て共通の $C,m,N$ で :eqref[eq-four-tempered-bound] を満たすとする。各 $\phi\in\mathcal D$ に対し $\langle T_j,\phi\rangle$ が収束すれば、その極限は $\mathcal S'$ の元に一意に延長され、$T_j$ は $\mathcal S'$ の弱い意味でそれに収束する。
:::

:::proof
$\phi\in\mathcal S$ を固定し、$\psi_R=\chi(x/R)\phi$ と置き $\chi$ を原点近くで $1$ の切断とする。$\psi_R$ は $\mathcal D$ に属し、$\phi-\psi_R$ は全 Schwartz セミノルムでゼロへ収束する。共通評価により $\langle T_j,\phi-\psi_R\rangle$ の上限も $j$ に関して一様にゼロへ近づく。よって $R$ を先に固定してから $\langle T_j,\psi_R\rangle$ が収束することを用い、$\langle T_j,\phi\rangle$ が Cauchy 列になると分かる。その極限へ共通の評価を移せるので、極限汎関数は連続である。一意性は $\mathcal D$ の稠密性から従う。
:::

:::example[無限遠へ逃げる点で二つの収束を比べる]
$T_j=j^2\delta_j$ は $\mathcal D'$ でゼロへ収束する。さらに $(1+|x|)^2|\phi(x)|$ の上限で作用を共通に評価できるので、$\mathcal S'$ でもゼロへ収束する。対して $U_j=e^{j^2}\delta_j$ は $\mathcal D'$ では同じ極限を持つが、$\phi(x)=e^{-x^2}$ に対して作用が常に $1$ である。もし $\mathcal S'$ で極限が存在すれば、その $\mathcal D$ への制限はゼロなので $\mathcal S'$ の元としてもゼロになるはずであり、矛盾する。個々の $U_j$ は緩増加だが、$j$ にわたる共通評価は持たない。
:::

## 緩増加超関数の Fourier 変換

:::definition[双対による Fourier 変換]{#def-four-tempered-transform}
$T\in\mathcal S'$ に対し

$$
\langle\widehat T,\phi\rangle=\langle T,\widehat\phi\rangle,\qquad
\langle\check T,\phi\rangle=\langle T,\check\phi\rangle
\label{eq-four-tempered-transform}
$$

で定める。試験関数側も同じ符号の変換を使う。
:::

:::proposition[通常の変換との一致と微分則]{#prop-four-tempered-rules}
この変換は $\mathcal S'$ の可逆線形写像であり、$L^1$、$L^2$ の変換と一致する。さらに

$$
\mathcal F(\partial^\alpha T)=(2\pi i\xi)^\alpha\widehat T,\qquad
\mathcal F(x^\alpha T)=(-2\pi i)^{-|\alpha|}\partial_\xi^\alpha\widehat T .
$$
:::

:::proof
:ref[thm-four-euclid-inversion] の位相同型を双対に移すので、作用は連続かつ逆変換を持つ。$T=T_f$、$f\in L^1$ では $\iint|f(x)\phi(\xi)|\,dx\,d\xi=\|f\|_1\|\phi\|_1<\infty$ により Fubini を用いれば積分変換と一致する。$L^2$ では Schwartz 近似と Cauchy–Schwarz で作用が両側で収束する。微分則は $\langle\widehat{\partial_j T},\phi\rangle=-\langle T,\partial_j\widehat\phi\rangle=\langle T,\widehat{2\pi i\xi_j\phi}\rangle$ である。座標倍では $x_j\widehat\phi=\widehat{\partial_j\phi}/(2\pi i)$ と双対微分の符号を用いる。反復で多重指数の場合が従う。
:::

:::example[定数、点、周期的な点列]
$\widehat{\delta_a}(\xi)=e^{-2\pi ia\cdot\xi}$ であり、特に $\widehat{\delta_0}=1$ である。また $\widehat 1=\delta_0$ は $\int\widehat\phi=\phi(0)$ から従う。一次元で $\widehat{x}=(i/(2\pi))\delta_0'$ である。:ref[thm-four-poisson-summation] により $\widehat{\sum_k\delta_k}=\sum_k\delta_k$ である。実際、作用は $\sum_k\widehat\phi(k)=\sum_k\phi(k)$ になる。
:::

:::example[多項式が反転で点に集中する]
$\widehat1=\delta_0$ に座標倍の公式を三回適用すると $\widehat{x^3}=-i\delta_0'''/(8\pi^3)$ となる。従って $\xi^3\widehat{x^3}$ は $(-i/(8\pi^3))(-6)\delta_0$ である。実際、$\xi^3\delta_0'''$ の作用は $-(\xi^3\phi)'''|_{\xi=0}=-6\phi(0)$ で、係数の符号を直接追える。
:::

### コンパクト台を変換の増大から読む

コンパクト台の超関数は複素指数にも作用できる。作用が試験関数の有限階の微分で抑えられるため、変換には多項式因子が加わる。その一方で、虚軸方向の指数の増大率は台の半径を保つ。

:::theorem[Paley–Wiener–Schwartz の球版]{#thm-four-paley-wiener-schwartz}
$R\ge0$ とする。$\operatorname{supp} T\subset\{x:|x|\le R\}$ の超関数 $T$ の変換は $\C^n$ 上で整関数となり、ある $C,m$ について
$$
|F(z)|\le C(1+|z|)^m e^{2\pi R|\Im z|}
$$
を満たす。逆にこの評価を満たす整関数は、その閉球に台を持つ超関数の変換である。
:::

:::proof
:ref[prop-four-compact-distribution] によって $F(z)=\langle T,e^{-2\pi ix\cdot z}\rangle$ を定められる。固定切断の台上では指数の冪級数が全ての $x$ 導関数を含め局所一様に収束するため、作用と交換して $F$ が整関数であることが分かる。実軸上の値が $\widehat T$ を表すことは、Schwartz 試験関数の積分で指数を重ね、有限階の $x$ 導関数を含めて積分を作用の外へ出せば従う。

評価のため $\delta=(2+|z|)^{-1}$ とし、半径 $R+\delta/2$ の球の近くで $1$、$R+2\delta$ 以内に台を持つ切断を使う。この切断は半径 $R+\delta$ の球の指示関数を半径 $\delta/4$ の滑らかな核で畳み込めばよく、各 $k$ 階の導関数は $C_k\delta^{-k}$ 以下である。固定した半径 $R+2$ の球の局所評価を切断付き指数に適用する。$m$ 階以下の微分は $(1+|z|)^m e^{2\pi(R+2\delta)|\Im z|}$ 以下の定数倍となる。$\delta|\Im z|\le1$ なので余分な指数因子は定数で吸収できる。

逆に実軸上で $F$ は多項式増大なので、$T=\check F$ を緩増加超関数として定める。閉球の外の任意の点の近傍では、ある単位ベクトル $v$ と $\delta>0$ で $x\cdot v>R+\delta$ となる。そこに台を持つ試験関数 $\phi$ を取る。$z=\xi+iYv$、$Y\ge0$ で
$$
|\check\phi(z)|\le C_l(1+Y)^{2l}(1+|\xi|^2)^{-l}e^{-2\pi(R+\delta)Y}
$$
である。実部の指数に $1-\Delta_x$ を $l$ 回使い、括弧内の試験関数と $e^{-2\pi Y x\cdot v}$ に微分を移せばこの評価になる。$F$ を掛け、$2l>n+m$ とすると $\xi$ について積分できる。$v$ を第一座標に移し、その他の座標を実数のまま保って第一座標を上へ移す。Cauchy の積分定理により、各固定 $Y$ で縦辺の寄与は $|\xi_1|\to\infty$ とともにゼロとなる。上の可積分な評価で他座標の積分も交換できる。従って
$$
\langle T,\phi\rangle=\int_{\R^n} F(\xi+iYv)\check\phi(\xi+iYv)\,d\xi .
$$
右辺は $C(1+Y)^{m+2l}e^{-2\pi\delta Y}$ 以下となる。$Y\to\infty$ でゼロになり、台は閉球に入る。
:::

:::example[球の端にある点質量]
一次元で $a>0$ とし、$T=\delta_a'$ とすれば $F(z)=2\pi iz\,e^{-2\pi iaz}$ である。従って上半平面の虚軸での大きさは $2\pi Y e^{2\pi aY}$ となる。前置の一次多項式はデルタの一階微分から来ており、指数率は台が原点から距離 $a$ にあることに対応する。
:::

### 平行移動、変調および離散格子

超関数の平行移動 $T(\cdot-a)$ は $\langle T,\phi(\cdot+a)\rangle$ で定める。試験関数側では移動の方向が反転する。乗数 $e^{2\pi i b\cdot x}$ は全ての導関数が有界なので、試験関数に掛けると連続に Schwartz 空間に残る。従って変調も $\mathcal S'$ 上で定義できる。双対変換の定義に Schwartz 関数の移動・変調の式を代入すると
$$
\mathcal F\bigl(T(\cdot-a)\bigr)=e^{-2\pi i a\cdot\xi}\widehat T,\qquad
\mathcal F(e^{2\pi i b\cdot x}T)=\widehat T(\cdot-b)
$$
を得る。例えば第一式では $\widehat\phi(x+a)=\mathcal F(e^{-2\pi i a\cdot\xi}\phi)(x)$ を用いる。第二式も $\mathcal F(\phi(\cdot+b))(x)=e^{2\pi i b\cdot x}\widehat\phi(x)$ から得られる。

:::proposition[格子デルタの変換]{#prop-four-tempered-lattice-comb}
可逆行列 $A$ と $v\in\R^n$ に対し、$\Delta_{A,v}=\sum_{k\in\Z^n}\delta_{v+Ak}$ と置く。この和は緩増加超関数として収束し、
$$
\widehat{\Delta_{A,v}}
=|\det A|^{-1}\sum_{m\in\Z^n}
 e^{-2\pi i v\cdot A^{-\mathsf T}m}\delta_{A^{-\mathsf T}m}
$$
である。
:::

:::proof
格子点は $|k|$ が大きいとき $|v+Ak|\ge c|k|$ を満たすので、$N>n$ の減衰重みを用いれば作用の和は絶対収束し連続である。$v=0$ の場合、任意の $\phi\in\mathcal S$ に対して $\langle\widehat\Delta_{A,0},\phi\rangle=\sum_k\widehat\phi(Ak)$ である。:eqref[eq-four-lattice-poisson] を $\widehat\phi$ に適用すると、右辺は $|\det A|^{-1}\sum_m\phi(-A^{-\mathsf T}m)$ である。符号を添字の置換で消せるので主張の $v=0$ を得る。一般の $v$ は格子の移動であり、上の公式で位相因子が付く。
:::

:::example[半格子にずらした点列]
$\Delta_{2,1}=\sum_k\delta_{1+2k}$ なら $\widehat{\Delta_{2,1}}=\tfrac12\sum_m(-1)^m\delta_{m/2}$ である。原点に点質量がない列でも、変換の零周波数には $1/2$ が現れる。これは長さ $2$ のセルあたり質量 $1$ の密度である。二つの隣り合う剰余類の列 $\Delta_{2,0}+\Delta_{2,1}$ では、奇数 $m$ の項が左右で打ち消し合い、偶数 $m$ の係数は $1$ になるため、整数格子の Fourier 不変性に戻る。
:::

### 周期超関数の係数

周期を持つ超関数は、点質量のように通常関数として書けない場合も指数で展開できる。その係数は急減少でなく、多項式程度に増えてよい。

:::proposition[周期超関数の Fourier 展開]{#prop-four-periodic-distribution-transform}
$L>0$ とする。$\R$ 上の周期 $L$ の超関数 $T$ には一意な多項式増大の係数列 $(c_k)_{k\in\Z}$ があり、
$$
T=\sum_{k\in\Z} c_k e^{2\pi i kx/L},\qquad
\widehat T=\sum_{k\in\Z} c_k\delta_{k/L}
$$
が $\mathcal S'$ における弱収束として成り立つ。逆に多項式増大の任意の列から周期超関数を得る。
:::

:::proof
コンパクト台の非負な滑らかな関数 $b$ で $[-L,L]$ 上で正のものを取り、$\chi=b/\sum_{j\in\Z} b(\cdot-jL)$ とする。分母は周期的で正であり、各コンパクト上の和は有限である。従って $\chi\in\mathcal D$ で $\sum_j\chi(x-jL)=1$ となる。コンパクト台の $\phi$ に対しては、周期性とこの分割を使うと
$$
\langle T,\phi\rangle
=\left\langle T,\chi(x)\sum_j\phi(x+jL)\right\rangle
$$
である。例えば右辺の第 $j$ 項を $-jL$ だけ移動すると、$\langle T,\chi(x-jL)\phi(x)\rangle$ になる。有限和で元の作用を得る。

$c_k=L^{-1}\langle T,\chi(x)e^{-2\pi i kx/L}\rangle$ と置く。固定した台上での評価と Leibniz 則により $|c_k|\le C(1+|k|)^m$ である。$\phi$ の周期化の Fourier 係数は $L^{-1}\widehat\phi(k/L)$ である。$\phi$ が滑らかなコンパクト台を持つため、周期化の級数は全ての階の導関数まで一様に収束する。その展開を右辺に入れると
$$
\langle T,\phi\rangle=\sum_{k\in\Z}c_k\widehat\phi(-k/L)
$$
となる。右辺は $\phi\in\mathcal S$ にも Schwartz セミノルムによって連続に延長できる。従って $T$ は緩増加であり、級数の形から変換すれば点質量の列を得る。実際 $\widehat{\widehat\phi}(-k/L)=\phi(k/L)$ である。

逆に多項式増大の列から同じ絶対収束する式で作用を定義すると、各指数が周期 $L$ を持つため極限も周期 $L$ である。変換側で一つの格子点の近くにだけ台を持ち、その点で $1$ の試験関数を使えば対応する係数を読めるから、係数は一意である。
:::

:::example[周期的なデルタの微分]
$T=\sum_{j\in\Z}\delta_{2j}'$ では :ref[prop-four-tempered-lattice-comb] と微分公式から $c_k=(2\pi i(k/2))/2=\pi i k/2$ である。例えば $c_0=0$、$c_2=\pi i$ になる。係数は一次に増えるが、試験関数の変換の急減少により展開は絶対収束する。一つのセル内で点の位置を $1$ へ移せば、各係数がさらに $(-1)^k$ 倍になる。
:::

### コンパクト台の因子を持つ畳み込み

:::proposition[緩増加超関数の畳み込みと積]{#prop-four-tempered-compact-convolution}
$S$ をコンパクト台の超関数、$T\in\mathcal S'$ とする。$S*T$ は緩増加超関数であり、$\widehat{S*T}=\widehat S\,\widehat T$ である。右辺では $\widehat S$ を実軸上の滑らかな関数として掛ける。
:::

:::proof
$\phi\in\mathcal S$ に対して
$$
\Psi(y)=\langle S_x,\phi(x+y)\rangle
$$
と置く。$x$ は固定したコンパクト集合内で局所評価されるため、$y$ に関する差商を作用と交換できる。また $x$ がその範囲にあると、$(1+|y|)^N$ は定数倍の $(1+|x+y|)^N$ 以下である。従って $\Psi$ の各 Schwartz セミノルムは $\phi$ の有限階のセミノルムで抑えられ、$\langle T,\Psi\rangle$ は緩増加超関数を定める。$\phi$ がコンパクト台なら前章の畳み込みの定義と一致する。

$\widehat S(\xi)=\langle S_x,e^{-2\pi i x\cdot\xi}\rangle$ であり、$\xi$ 微分は $x$ の多項式を掛ける操作である。その各微分も固定台と有限階の評価によって多項式増大する。よって $\widehat S\phi\in\mathcal S$ である。$\widehat\phi(x+y)$ を指数の積分で展開すると、
$$
\langle S_x,\widehat\phi(x+y)\rangle
=\int_{\R^n} \widehat S(\xi)\phi(\xi)e^{-2\pi i y\cdot\xi}d\xi
=\mathcal F(\widehat S\phi)(y).
$$
積分と作用の交換は、有限階の $x$ 導関数を指数に取っても $\phi$ の急減少が可積分性を与えるためである。従って $\langle\widehat {S*T},\phi\rangle=\langle T,\mathcal F(\widehat S\phi)\rangle=\langle\widehat T,\widehat S\phi\rangle$ となる。
:::

:::example[二点の平均で離散周波数を消す]
$S=(\delta_0+\delta_1)/2$、$T=\Delta_{2,0}$ とする。空間側では $S*T$ は整数格子のデルタの半分である。周波数側では $\widehat S(\xi)=(1+e^{-2\pi i\xi})/2$ であり、$T$ の第 $m$ 点 $m/2$ での値は $(1+(-1)^m)/2$ である。奇数 $m$ の半整数点は消え、偶数 $m$ の整数点に元の係数 $1/2$ が残る。
:::

### Schwartz 関数との畳み込みで平滑化する

平滑化に Gauss 関数を使う場合、その台はコンパクトではない。緩増加超関数には急減少関数を試験関数として使えるため、台の代わりに Schwartz セミノルムを使って畳み込みを定義できる。

:::proposition[Schwartz 関数による平滑化]{#prop-four-tempered-schwartz-convolution}
$T\in\mathcal S'(\R^n)$、$g\in\mathcal S(\R^n)$ に対して $u(x)=\langle T_y,g(x-y)\rangle$ と置く。$u$ は滑らかで、各導関数は多項式増大を持つ。また $\widehat u=\widehat g\,\widehat T$ である。
:::

:::proof
$x$ をコンパクト集合に限ると、Taylor の剰余を各 $y$ 導関数に適用して、$x$ に関する差商が対応する導関数へ全 Schwartz セミノルムで収束すると分かる。重みは $(1+|y|)^N\le (1+|x|)^N(1+|x-y|)^N$ で扱える。従って $\partial^\beta u(x)=\langle T_y,\partial^\beta g(x-y)\rangle$ である。この重みの不等式を :eqref[eq-four-tempered-bound] に使えば、各導関数は $C_\beta(1+|x|)^N$ 以下となる。

$\phi\in\mathcal S$ に対し、$\langle\widehat u,\phi\rangle=\int \widehat\phi(x)\langle T_y,g(x-y)\rangle\,dx$ を考える。積分は Schwartz セミノルムごとに収束する。実際、$y$ に関するセミノルムでの移動の増大は $(1+|x|)^N$ 以下で、$\widehat\phi$ は急減少する。従って有界領域での Riemann 和から始め、次に積分の範囲を広げれば、積分を作用の内側へ移せる。その内側の積分に $\widehat\phi$ の表示を入れ、$v=x-y$ と置くと
$$
\int g(x-y)\widehat\phi(x)dx
=\int e^{-2\pi i y\cdot\xi}\widehat g(\xi)\phi(\xi)d\xi
=\mathcal F(\widehat g\phi)(y).
$$
ここでは $g,\phi\in L^1$ により Fubini を使った。従って $\langle\widehat u,\phi\rangle=\langle\widehat T,\widehat g\phi\rangle$ となる。
:::

例えば $T=\delta_0'$ と $g(x)=e^{-\pi x^2}$ なら $u=g'=-2\pi x e^{-\pi x^2}$ であり、変換は $2\pi i\xi e^{-\pi \xi^2}$ である。

## 例：主値・Heaviside 関数・斉次超関数

主値では原点の両側を同じ幅だけ除く。変換の符号と係数を追うため、まず原点を平滑化した有理関数を変換する。

:::proposition[主値と階段関数の変換]{#prop-four-tempered-pv-heaviside}
一次元で

$$
\mathcal F\left(\operatorname{pv}\frac1x\right)=-i\pi\sgn\xi,\qquad
\widehat H=\frac12\delta_0+\frac1{2\pi i}\operatorname{pv}\frac1{\xi}.
$$
:::

:::proof
$x/(x^2+a^2)$ は $a\downarrow0$ で $\operatorname{pv}(1/x)$ に作用の意味で収束する。正負の半軸でまとめると被積分関数は $x(\phi(x)-\phi(-x))/(x^2+a^2)$ である。$0<x<1$ では $2\sup|\phi'|$、$x\ge1$ では Schwartz 減衰で支配できる。この評価は極限が $\mathcal S'$ に入ることも示す。

しかし変換側では、:ref[thm-four-euclid-inversion] の Poisson 核の例より $(x^2+a^2)^{-1}$ の変換が $(\pi/a)e^{-2\pi a|\xi|}$ である。:ref[prop-four-tempered-rules] でこれを微分すると $-i\pi\sgn(\xi)e^{-2\pi a|\xi|}$ になる。折れ点の値は連続なので一階微分にデルタ項はない。優収束で右辺は $-i\pi\sgn\xi$ になる。

:ref[exm-four-step-derivative] により $2\pi i\xi\widehat H=1$ である。$\xi\operatorname{pv}(1/\xi)=1$ なので、差は $\xi$ 倍でゼロとなる。:ref[thm-four-point-supported] と $\xi\delta^{(k)}=-k\delta^{(k-1)}$ を使うと、差は $c\delta_0$ に限られる。$H(x)+H(-x)=1$ を変換し、主値の奇性を見ると $2c=1$ である。
:::

:::definition[斉次超関数]{#def-four-homogeneous-distribution}
$\lambda>0$ に対して $\langle T(\lambda\cdot),\phi\rangle=\lambda^{-n}\langle T,\phi(\cdot/\lambda)\rangle$ とする。$T(\lambda\cdot)=\lambda^dT$ のとき :term[斉次超関数]{en="homogeneous distribution" yomi="せいじちょうかんすう"} で次数 $d$ という。通常の関数なら変数変換から点ごとの斉次性に一致する。
:::

:eqref[eq-four-tempered-transform] と尺度変換により $\mathcal F(T(\lambda\cdot))=\lambda^{-n}\widehat T(\cdot/\lambda)$ である。従って斉次次数は $d\mapsto -n-d$ に移る。例えば $\delta_0$ の次数は $-n$、定数 $1$ の次数は $0$ である。$\operatorname{pv}(1/x)$ は対称切断を $\lambda$ 倍すれば次数 $-1$ になり、変換された符号関数は次数 $0$ である。

:::proposition[放射状の冪核]{#prop-four-riesz-kernel}
$0<\alpha<n$ なら

$$
\mathcal F(|x|^{-\alpha})(\xi)=
\pi^{\alpha-n/2}\frac{\Gamma((n-\alpha)/2)}{\Gamma(\alpha/2)}
|\xi|^{\alpha-n}.
$$
:::

:::proof
$x\ne0$ で

$$
|x|^{-\alpha}=\frac{\pi^{\alpha/2}}{\Gamma(\alpha/2)}
\int_0^\infty t^{\alpha/2-1}e^{-\pi t|x|^2}dt
$$

である。これは $u=\pi t|x|^2$ の置換による。$\widehat\phi$ と積分するとき、$t$ を先に積分した絶対値の上界は $C\int |x|^{-\alpha}|\widehat\phi(x)|dx$ で有限である。Gauss 変換で内側を $t^{-n/2}e^{-\pi|\xi|^2/t}$ に変える。$u=\pi|\xi|^2/t$ と置いて $t$ を先に積分すると定数を掛けた $|\xi|^{\alpha-n}$ になる。これは原点で局所可積分であり、$\phi$ は急減少するので二回目の Fubini も正当化される。
:::

:::example[三次元の Coulomb 核]
$n=3$、$\alpha=1$ では定数は $1/\pi$ であり、$\mathcal F((4\pi|x|)^{-1})=(4\pi^2|\xi|^2)^{-1}$ である。これに $4\pi^2|\xi|^2$ を掛けると $1$ になるので、三次元の $-\Delta$ の基本解を周波数側でも確認できる。
:::

:::example[原点の発散を正則化で追う]
$x/(x^2+a^2)$ の変換は原点の左右で符号が異なる。$a=1/(2\pi)$ なら $\xi=1$ で $-i\pi/e$、$\xi=-1$ で $i\pi/e$ である。$a\downarrow0$ ではこの絶対値が $\pi$ に増える。原点での関数値には超関数の作用として意味がなく、左右の差だけが残る。
:::

### 複素平面から境界へ近づける

原点の極を実軸の上側と下側へ少し動かすと、それぞれに対して普通の関数として積分できる。境界へ近づけると、極の対称的な部分と点に集中する部分が分かれる。

:::proposition[一階極の境界超関数]{#prop-four-sokhotski-plemelj}
$\eps>0$ から $\eps\downarrow0$ とすると、$\mathcal S'(\R)$ で
$$
\frac1{x\pm i\eps}\longrightarrow
\operatorname{pv}\frac1x\mp i\pi\delta_0
$$
である。従って二つの境界値の変換は、上の符号で $-2\pi iH(\xi)$、下の符号で $2\pi iH(-\xi)$ である。
:::

:::proof
$1/(x\pm i\eps)=x/(x^2+\eps^2)\mp i\eps/(x^2+\eps^2)$ と分ける。実部は :ref[prop-four-tempered-pv-heaviside] の証明と同じく主値に近づく。虚部の核に $x=\eps t$ を入れると
$$
\int_\R\frac{\eps}{x^2+\eps^2}\phi(x)dx
=\int_\R \frac{\phi(\eps t)}{1+t^2}dt
\longrightarrow \pi\phi(0).
$$
支配関数は $\|\phi\|_\infty/(1+t^2)$ である。主値とデルタの変換を代入すれば、上の符号では $-i\pi(\sgn\xi+1)$、下の符号では $-i\pi(\sgn\xi-1)$ となる。これを半直線の指示関数に書き直すと主張の式である。
:::

:::example[二階の極の符号]
$\partial_x(1/(x+i\eps))=-(x+i\eps)^{-2}$ であるから、弱極限を微分して
$$
(x+i0)^{-2}=-\partial_x\operatorname{pv}(1/x)+i\pi\delta_0'
$$
を得る。右辺の変換は、まず上側の境界値の変換 $-2\pi iH(\xi)$ を $-2\pi i\xi$ 倍するため $-4\pi^2\xi H(\xi)$ になる。実部に相当する $-\partial_x\operatorname{pv}(1/x)$ の変換は $-2\pi^2|\xi|$、デルタの導関数の変換は $-2\pi^2\xi$ であり、足し合わせると同じ片側の式になる。
:::

### Euler の式と臨界次数

通常の斉次関数の式を尺度について微分すると、Euler の式が得られる。超関数でも試験関数に対する作用として微分すれば同じ式になる。それにより、原点の特異性を取り除いたときに斉次性を保てたかを調べられる。

:::proposition[斉次性の微分による判定]{#prop-four-euler-homogeneity}
$T\in\mathcal S'$ が次数 $d$ で斉次であることと
$$
\sum_{j=1}^n x_j\partial_j T=dT
$$
は同値である。
:::

:::proof
$\lambda=e^t$ として $\langle T(e^t\cdot),\phi\rangle=e^{-nt}\langle T,\phi(e^{-t}\cdot)\rangle$ を $t$ で微分する。$t$ が有界区間を動くとき、合成とその差商は各 Schwartz セミノルムで連続であり、微分の公式は通常の連鎖律と Taylor の剰余評価から成り立つ。$t=0$ での微分は $\langle T,-n\phi-x\cdot\nabla\phi\rangle$ である。これは積の微分の定義により $\langle\sum_jx_j\partial_jT,\phi\rangle$ に等しい。斉次なら同じ微分は $d\langle T,\phi\rangle$ である。逆に Euler の式を仮定すると、尺度変換の群法則から任意の $t$ における微分も $d\langle T(e^t\cdot),\phi\rangle$ である。この一次常微分方程式を解くと $e^{dt}\langle T,\phi\rangle$ になる。
:::

:::example[ちょうど非可積分になる冪を延長する]
$\sigma_{n-1}$ を単位球面の面積とする。$|x|^{-n}$ の有限部分を
$$
\langle U,\phi\rangle=
\int_{|x|<1}|x|^{-n}(\phi(x)-\phi(0))dx+
\int_{|x|\ge1}|x|^{-n}\phi(x)dx
$$
と定める。最初の差は $|x|\sup|\nabla\phi|$ 以下なので半径積分が収束し、第二項は Schwartz 減衰で収束する。この式は $\eps$ 球を除いた積分から $\sigma_{n-1}\log(1/\eps)\phi(0)$ を引いて得られる極限である。従って尺度変換では
$$
U(\lambda\cdot)=\lambda^{-n}U+
 \lambda^{-n}\sigma_{n-1}\log\lambda\,\delta_0
$$
となる。実際、試験関数内の $x/\lambda$ を $y$ と置けば、除外半径は $\eps/\lambda$ になり、その差が $\log\lambda$ の項を残す。別の超関数への延長との差は原点に台を持つので、:ref[thm-four-point-supported] によってデルタの導関数の有限和である。この差に原点の近くで $1$ の試験関数を入れると、正の階数の導関数は尺度変換の前後とも消え、デルタ自身は次数 $-n$ のため変換の差に寄与しない。従っていずれの延長も対数項を消せず、次数 $-n$ の斉次延長は存在しない。
:::

平均ゼロの角度成分では、この対数項が消える。$\Omega$ を球面上の有界関数とし $\int_{S^{n-1}}\Omega=0$ とする。すると $K(x)=|x|^{-n}\Omega(x/|x|)$ を半径 $\eps$ の球の外で積分した作用は、$\eps\downarrow0$ で超関数に収束する。実際 $\eps<|x|<1$ の領域で $\phi(0)$ を引いても球面平均がゼロなので積分は変わらず、上と同じ収束評価が成り立つ。尺度の変更で新たに加わる環上の定数項も角度積分がゼロであるため、作用は次数 $-n$ の斉次性を持つ。例えば一次元で $\Omega(1)=1$、$\Omega(-1)=-1$ とすれば対称主値 $1/x$ に一致する。二次元で $\Omega(\omega)=\omega_1$ とすれば反射 $x_1\mapsto-x_1$ に関して偶な試験関数で作用はゼロになる。

## Fourier 変換による Sobolev 空間 Hˢ

実数 $s$ で高周波側を何冪まで制御するかを定める。整数の微分を組み合わせるよりも、この方法では負や非整数の指数も同じ式で扱える。

:::definition[二乗可積分型 Sobolev 空間]{#def-four-sobolev}
$\langle\xi\rangle=(1+4\pi^2|\xi|^2)^{1/2}$ とし、$s\in\R$ に対して

$$
H^s=\{T\in\mathcal S':\widehat T\text{ が関数で }\langle\xi\rangle^s\widehat T\in L^2\},\qquad
\|T\|_{H^s}=\|\langle\xi\rangle^s\widehat T\|_2 .
$$

これを :term[Sobolev 空間]{en="Sobolev space" yomi="そぼれふくうかん"} という。
:::

:::proposition[Hilbert 構造と整数階、双対]{#prop-four-sobolev-structure}
重み付き変換 $J_s:T\mapsto\langle\xi\rangle^s\widehat T$ は $H^s\to L^2$ の等長同型である。非負整数 $m$ では $H^m$ は $|\alpha|\le m$ の超関数微分が全て $L^2$ にある空間と一致し、ノルムは同値である。$H^s$ の連続線形双対は $H^{-s}$ と同定できる。
:::

:::proof
$w\in L^2$ に対して $\langle\xi\rangle^{-s}w$ は Schwartz 関数との積分で緩増加超関数を定める。Cauchy–Schwarz で $\|\langle\xi\rangle^{-s}\phi\|_2$ によって作用を評価できるからである。逆変換を取れば全射性を得る。完全性と Hilbert 内積は $L^2$ から移る。

整数 $m$ では $(1+4\pi^2|\xi|^2)^m$ を多項展開すると、$|(2\pi\xi)^\alpha|^2$（$|\alpha|\le m$）の正の係数の和になる。各微分の変換に :ref[prop-four-tempered-rules] と :ref[thm-four-plancherel] を使えば、ノルム同値性と空間の一致が従う。

$T\in H^{-s}$ は $\phi\mapsto\int \widehat T(\xi)\widehat\phi(-\xi)d\xi$ で $H^s$ 上の連続線形汎関数を定める。二つの逆の重みを入れて Cauchy–Schwarz を使う。逆に、任意の連続線形汎関数 $L$ を $J_s$ により $L^2$ へ移すと、:ref[thm-fa-riesz-representation] によって $L(\phi)=\int\langle\eta\rangle^s\widehat\phi(\eta)\overline{w(\eta)}\,d\eta$ と書ける。ただし $w\in L^2$ である。ここで $\widehat T(\xi)=\langle\xi\rangle^s\overline{w(-\xi)}$ と置けば $T\in H^{-s}$ であり、変数を反転するとちょうど上の組み合わせになる。Schwartz 関数との作用の一致は :eqref[eq-four-tempered-transform] と反転に従う。
:::

:::example[デルタが属する指数]
$\widehat{\delta_0}=1$ なので $\delta_0\in H^s$ は $2s<-n$ と同値である。球座標で半径 $r\ge1$ の積分が $\int r^{n-1+2s}dr$ と比較できるからである。$s=-n/2$ は対数発散する。一次元では $\|\delta_0\|_{H^{-1}}^2=\int(1+4\pi^2\xi^2)^{-1}d\xi=1/2$ である。
:::

:::example[一階ノルムの具体値]
$f(x)=e^{-\pi x^2}$ なら $\|f\|_{H^1}^2=\int(1+4\pi^2\xi^2)e^{-2\pi\xi^2}d\xi$ である。Gauss 積分とその係数の微分で $\int e^{-2\pi\xi^2}=1/\sqrt2$、$\int\xi^2e^{-2\pi\xi^2}=1/(4\sqrt2\pi)$ なので、結果は $(1+\pi)/\sqrt2$ である。空間側の $\|f\|_2^2+\|f'\|_2^2$ による計算とも一致する。
:::

### 指数の間の包含と近似

重み $\langle\xi\rangle$ が $1$ 以上なので、$s_1\ge s_0$ なら $\|u\|_{H^{s_0}}\le\|u\|_{H^{s_1}}$ である。中間の指数 $s=(1-\theta)s_0+\theta s_1$（$0<\theta<1$）では
$$
\|u\|_{H^s}\le\|u\|_{H^{s_0}}^{1-\theta}\|u\|_{H^{s_1}}^\theta
$$
である。二乗ノルムの被積分関数が、両端の二乗ノルムの被積分関数の $1-\theta$ 乗と $\theta$ 乗の積になるため、Hölder の不等式で従う。この式で指数が高くなると強い仮定になる理由は、座標微分についての正の整数の場合だけでなく負の指数でも同じである。

周波数側の滑らかなコンパクト台関数の逆変換からなる集合は、各 $H^s$ に稠密である。実際 $J_s$ によって $L^2$ に移し、$L^2$ 関数を $C_c^\infty$ で近似してから $\langle\xi\rangle^{-s}$ を掛けても、滑らかなコンパクト台関数のままである。従って Schwartz 関数は稠密である。さらに Schwartz 関数は前節の切断で $\mathcal D$ の関数によって $\mathcal S$ で近似できる。Schwartz 収束は Fourier 側でも成り立ち、十分高い減衰冪を取れば重み付き $L^2$ 収束を含む。ゆえに $C_c^\infty(\R^n)$ も各 $H^s$ に稠密である。

:::example[周波数切断で誤差率を得る]
$u\in H^{s+r}$、$r>0$ とし、$\widehat{u_R}=1_{\{|\xi|\le R\}}\widehat u$ とする。誤差には $|\xi|>R$ の成分しかないため $\|u-u_R\|_{H^s}\le(1+4\pi^2R^2)^{-r/2}\|u\|_{H^{s+r}}$ である。$s=0,r=2$ の場合、$H^2$ のノルムが分かっていれば、$R$ を大きくするにつれ $L^2$ 誤差を $O(R^{-2})$ に抑えられる。$u_R$ の値を局所的に切るのではなく周波数を除いているため、$u_R$ は一般にはコンパクト台を持たない。
:::

### 滑らかな関数を掛けて局所化する

空間の中である領域に話を限るには、切断関数が必要である。Fourier 側では掛け算が畳み込みになるため、もとの重みがこの操作で保たれるかを評価する。

:::proposition[滑らかな局所化の有界性]{#prop-four-sobolev-localization}
$\psi\in\mathcal S(\R^n)$、$s\in\R$ とする。$u\mapsto\psi u$ は $H^s$ 上の有界作用素であり、
$$
\|\psi u\|_{H^s}\le
2^{|s|/2}\|\langle\cdot\rangle^{|s|}\widehat\psi\|_1\|u\|_{H^s}.
$$
:::

:::proof
$\langle a+b\rangle\le\sqrt2\langle a\rangle\langle b\rangle$ は $1+4\pi^2|a+b|^2\le2(1+4\pi^2|a|^2)(1+4\pi^2|b|^2)$ から従う。$s\ge0$ ではこれを $a=\xi-\eta,b=\eta$ に使う。$s<0$ では $\langle\eta\rangle\le\sqrt2\langle\xi-\eta\rangle\langle\xi\rangle$ を使う。いずれの場合も
$$
\langle\xi\rangle^s\le
2^{|s|/2}\langle\xi-\eta\rangle^{|s|}\langle\eta\rangle^s
$$
となる。まず $u\in\mathcal S$ とし、$\widehat{\psi u}=\widehat\psi*\widehat u$ にこの評価を入れる。$L^1$ と $L^2$ の畳み込みを $L^2$ で評価して主張を得る。一般の $u$ へは稠密性で拡張する。$H^s$ 収束は試験関数の作用での収束を含み、$\phi\mapsto\psi\phi$ は $\mathcal S$ 上連続なので、拡張した作用は超関数としての積と一致する。
:::

:::example[負の指数で点の位置を局所化する]
一次元で $u=\delta_0+2\delta_1$ は $s<-1/2$ なら $H^s$ に属する。$\psi$ を原点の近くで $1$、点 $1$ の近くでゼロとすれば $\psi u=\delta_0$ である。同じ切断を $\delta_0'$ に掛けても $\psi'(0)=0$ なので $\delta_0'$ のままである。一方原点で $\psi(0)=0,\psi'(0)=1$ の乗数なら $\psi\delta_0'=-\delta_0$ となる。作用素評価の仮定は $\psi$ が滑らかであることにあり、$u$ が通常の関数である必要はない。
:::

### Bessel 作用素による指数の移動

$\mathcal J^t u=\mathcal F^{-1}(\langle\xi\rangle^t\widehat u)$ と定める。$\langle\xi\rangle^t$ とその導関数が全て多項式増大以下なので、この乗算は $\mathcal S'$ 上で定義できる。$u\in H^s$ なら $\|\mathcal J^t u\|_{H^{s-t}}=\|u\|_{H^s}$ であり、逆写像は $\mathcal J^{-t}$ である。特に $(1-\Delta)=\mathcal J^2$ は $H^{s+2}$ から $H^s$ への等長同型である。微分の次数と同じく高周波を二冪だけ増やすが、記号の下限が $1$ なので零周波数でも逆を持つ。

:::example[解いた方程式から二階分の正則性を得る]
一次元で $(1-\partial_x^2)u=\delta_0$ とする。周波数側で $u$ は $(1+4\pi^2\xi^2)^{-1}$ を変換に持つから、$H^s$ への所属条件は $\int_{|\xi|\ge1}|\xi|^{2s-4}d\xi<\infty$、すなわち $s<3/2$ である。右辺のデルタの条件が $s<-1/2$ であったことと比べると、二階分だけ上限が大きい。空間側の解は $e^{-|x|}/2$ であり、連続だが一階導関数は原点で跳ぶので、整数階では $H^1$ までに属する。周波数側の計算ではその二つの整数の間の正確な閾値も読める。
:::

## L² 型の Sobolev の埋め込み定理

高周波の重みが十分強ければ、Cauchy–Schwarz により変換が可積分になる。そのとき逆変換が連続代表元を与え、微分した場合にも必要な重みを計算できる。

:::theorem[連続性と Hölder 連続性]{#thm-four-sobolev-embedding}
整数 $k\ge0$ で $s>k+n/2$ なら、$H^s$ は導関数が $k$ 階まで連続で無限遠で消える関数の空間に連続に埋め込まれる。さらに $0<\gamma\le1$、$s>k+\gamma+n/2$ なら $k$ 階導関数は $\gamma$-Hölder 連続である。
:::

:::proof
$|\alpha|\le k$ について

$$
\|(2\pi i\xi)^\alpha\widehat T\|_1
\le\|T\|_{H^s}\left(\int |(2\pi\xi)^\alpha|^2\langle\xi\rangle^{-2s}\,d\xi\right)^{1/2}
$$

である。積分は原点付近で有限であり、無限遠では半径の指数が $n-1+2|\alpha|-2s<-1$ となる。このため逆積分の微分が許され、:ref[thm-four-euclid-rl] により各導関数は無限遠で消える。逆積分が元の超関数を表すことは試験関数との Fubini と反転から従う。差には $|e^{2\pi ih\cdot\xi}-1|\le C_\gamma|h|^\gamma|\xi|^\gamma$ を入れ、同じ評価の $|\alpha|$ を $|\alpha|+\gamma$ に置き換える。
:::

:::example[厳密な閾値]
$R>2$ で $\widehat u_R=1_{\{1\le|\xi|\le R\}}\langle\xi\rangle^{-2s}/A_R^{1/2}$、$A_R=\int_{1\le|\xi|\le R}\langle\xi\rangle^{-2s}d\xi$ とする。$\|u_R\|_{H^s}=1$ だが、$u_R(0)=A_R^{1/2}$ である。$s\le n/2$ なら $A_R\to\infty$ なので、$H^s$ ノルムで連続代表元の上限ノルムを一様に抑えることはできない。
:::

### 高い正則性で積を取る

:::corollary[Sobolev 空間の積の評価]{#cor-four-sobolev-algebra}
$s>n/2$ なら $u,v\in H^s$ の連続代表元の積も $H^s$ に属し、$\|uv\|_{H^s}\le C_{n,s}\|u\|_{H^s}\|v\|_{H^s}$ が成り立つ。
:::

:::proof
$\langle a+b\rangle\le\langle a\rangle+\langle b\rangle$ は、$(1,2\pi(a+b))=(1,2\pi a)+(0,2\pi b)$ と分けて三角不等式を使えば従う。従って $s>0$ に対して $\langle a+b\rangle^s\le C_s(\langle a\rangle^s+\langle b\rangle^s)$ である。Schwartz 関数では積の変換は畳み込みなので、重みをいずれか一方の変換に掛け、$L^2$ と $L^1$ の畳み込みで評価すると
$$
\|uv\|_{H^s}\le C_s\bigl(\|u\|_{H^s}\|\widehat v\|_1+
                             \|\widehat u\|_1\|v\|_{H^s}\bigr)
$$
となる。$s>n/2$ なら :ref[thm-four-sobolev-embedding] の証明中の Cauchy–Schwarz 評価で $\|\widehat u\|_1\le C_{n,s}\|u\|_{H^s}$ が成り立つ。$u,v$ をそれぞれ Schwartz 関数で近似すると、得られた双線形評価により積は $H^s$ で Cauchy 列となる。同時に、埋め込みによって両因子はそれぞれの連続代表元へ一様収束するので、$H^s$ の極限はこの積に一致する。
:::

### 二点間の差で分数階を測る

非整数の $s$ では個々の通常微分を数えることはできない。しかし関数を少し移した差を距離に応じて積分すれば、同じ周波数の重みを得る。

:::proposition[分数階の差のノルム]{#prop-four-sobolev-differences}
$0<s<1$、$f\in L^2(\R^n)$ とすると
$$
\int_{\R^n}\int_{\R^n}\frac{|f(x+y)-f(x)|^2}{|y|^{n+2s}}\,dx\,dy
=C_{n,s}\int_{\R^n}|\xi|^{2s}|\widehat f(\xi)|^2\,d\xi ,
$$
ただし $C_{n,s}=\int_{\R^n}|e^{2\pi iy_1}-1|^2|y|^{-n-2s}dy\in(0,\infty)$ である。従って左辺と $\|f\|_2^2$ の和の平方根は $H^s$ ノルムと同値である。
:::

:::proof
固定した $y$ で差の変換は $(e^{2\pi iy\cdot\xi}-1)\widehat f(\xi)$ である。Plancherel と Tonelli によって $\xi,y$ の積分を交換できる。$\xi\ne0$ では直交変換で方向を $e_1$ に合わせ、$y$ を $|\xi|$ 倍して $|\xi|^{2s}$ を外に出す。原点近くの内側の積分では分子が $O(|y|^2)$、遠方では $4$ 以下である。球座標の半径積分が前者で $s<1$、後者で $s>0$ によって収束する。積分は正測度集合で正なので定数は正である。最後に $1+|\xi|^{2s}$ と $\langle\xi\rangle^{2s}$ は定数倍の範囲で等しい。
:::

:::example[区間指示関数の閾値]
$f=1_{[0,1]}$ では $\int |f(x+y)-f(x)|^2dx=2\min(1,|y|)$ である。差は二つの区間の対称差で絶対値 $1$ になるからである。$0<s<1$ で上の二重積分は
$$
4\left(\int_0^1 y^{-2s}dy+\int_1^\infty y^{-1-2s}dy\right)
$$
であり、有限となるのは $s<1/2$ のときに限る。その場合の値は $4((1-2s)^{-1}+(2s)^{-1})$ である。跳躍がこの分数階での境界を決める。
:::

### 有限の指数への Sobolev 埋め込み

:::proposition[非臨界の $L^q$ 埋め込み]{#prop-four-sobolev-lq}
$2<q<\infty$、$s>n(1/2-1/q)$ なら $\|f\|_q\le C_{n,s,q}\|f\|_{H^s}$ である。
:::

:::proof
$1/q'=1/2+1/r$ と置くと $1/r=1/2-1/q$ である。仮定より $\langle\xi\rangle^{-s}\in L^r$ なので Hölder から $\|\widehat f\|_{q'}\le\|\langle\xi\rangle^{-s}\|_r\|f\|_{H^s}$ を得る。$q'<2$ に対して :ref[thm-four-hausdorff-young] を逆変換に適用すれば結論になる。まず Schwartz 関数上で行い、:ref[prop-four-sobolev-structure] の等長写像で周波数側を滑らかなコンパクト台関数で近似して一般の場合に移す。
:::

:::example[二次元での四乗可積分性]
$n=2,q=4$ なら条件は $s>2(1/2-1/4)=1/2$ である。例えば $H^{3/4}$ の関数は $L^4$ に属する。連続代表元を保証する前節の条件 $s>1$ より必要な指数が小さいのは、点の上限の代わりに積分ノルムを測っているためである。
:::

### 臨界指数での埋め込み

前の Hölder による証明では、ちょうど $s=n(1/2-1/q)$ にすると重みの積分が対数発散する。この指数では、空間側の冪核を近い部分と遠い部分に分け、その境界を入力に応じて選ぶ。

:::lemma[分数積分の評価]{#lem-four-fractional-integral}
$0<s<n/2$、$g\in\mathcal S(\R^n)$ とする。$I_sg=c_{n,s}|\,\cdot\,|^{s-n}*g$ とし、$c_{n,s}>0$ は :ref[prop-four-riesz-kernel] により変換が $|2\pi\xi|^{-s}$ になるように選ぶ。$q=2n/(n-2s)$ で $\|I_sg\|_q\le C_{n,s}\|g\|_2$ である。
:::

:::proof
$|y|\le R$ の積分を環 $2^{-j-1}R<|y|\le 2^{-j}R$ に分ける。各項は $C(2^{-j}R)^s Mg(x)$ 以下で、$M$ は :ref[def-meas-maximal] の極大関数である。従って内側は $CR^s Mg(x)$ 以下である。外側は Cauchy–Schwarz で $C R^{s-n/2}\|g\|_2$ 以下である。実際、核の二乗の半径積分は $\int_R^\infty r^{2s-n-1}dr$ である。$R^{n/2}=\|g\|_2/Mg(x)$ と選べば
$$
|I_sg(x)|\le C(Mg(x))^{1-2s/n}\|g\|_2^{2s/n}
$$
となる。$Mg(x)=0$ なら全ての球で $|g|$ の積分がゼロなので $g=0$ である。最後に $q(1-2s/n)=2$ を使う。$M$ の $L^2$ 有界性は :ref[thm-meas-maximal-weak] と自明な $L^\infty$ 評価を :ref[thm-meas-marcinkiewicz] で補間して得る。よって $q$ 乗して積分すれば主張が従う。
:::

:::theorem[臨界の $L^q$ 埋め込み]{#thm-four-sobolev-critical}
$0<s<n/2$、$q=2n/(n-2s)$ なら $\|f\|_q\le C\|f\|_{H^s}$ である。
:::

:::proof
$\widehat f$ が原点を避けた滑らかなコンパクト台関数の場合、$\widehat g=|2\pi\xi|^s\widehat f$ と置けば $g\in\mathcal S$ である。:ref[prop-four-riesz-kernel] によって $I_sg=f$ となる。変換の積公式は、Schwartz 関数と冪核の畳み込みを試験関数へ作用させ、積分順序を交換すると得られる。絶対可積分性は原点での冪核の可積分性と、他の二因子の急減少から従う。従って :ref[lem-four-fractional-integral] と Plancherel で $\|f\|_q\le C\|g\|_2\le C\|f\|_{H^s}$ となる。原点を避けた滑らかなコンパクト台関数は重み付き $L^2$ で稠密である。原点近くと遠方を切り、その間で近似すればよい。$H^s$ と $L^q$ の両極限は試験関数への作用で一致するので任意の $f$ に移る。
:::

:::example[尺度で指数を確かめる]
$n=3,s=1$ なら $q=6$ である。$f_\lambda(x)=f(\lambda x)$ とすると $\|f_\lambda\|_6=\lambda^{-1/2}\|f\|_6$、$\|\nabla f_\lambda\|_2=\lambda^{1-3/2}\|\nabla f\|_2=\lambda^{-1/2}\|\nabla f\|_2$ である。両辺の尺度が一致する。$q>6$ なら $\lambda\to\infty$ で $\|f_\lambda\|_q$ の減衰の方が遅く、全ての $f$ に対する同じ評価は成り立たない。
:::

### 超平面へ制限するための閾値

$H^s$ の元に対して、座標 $x_n=0$ で値を取る操作を調べる。空間全体の点で値を取るには $s>n/2$ が必要であったが、一つの法線方向だけを固定する場合は、周波数積分を一変数だけ実行する。この違いが $1/2$ の閾値を与える。

:::theorem[超平面の跡と延長]{#thm-four-sobolev-trace-plane}
$n\ge2$ とする。$s>1/2$ なら Schwartz 関数の制限 $u(x',0)$ は連続な全射
$$
\operatorname{Tr}:H^s(\R^n)\longrightarrow H^{s-1/2}(\R^{n-1})
$$
に一意に延長され、連続な線形右逆写像を持つ。$s\le1/2$ の場合、Schwartz 関数の制限を保つ $H^s(\R^n)\to\mathcal D'(\R^{n-1})$ の連続写像は存在しない。ここでの連続性は少なくとも全ての固定試験関数に対する作用の連続性を意味する。
:::

:::proof
$\xi=(\xi',t)$、$a=(1+4\pi^2|\xi'|^2)^{1/2}$ とする。Schwartz 関数で逆変換を先に法線方向へ積分すると $\mathcal F_{x'}(\operatorname{Tr}u)(\xi')=\int_\R\widehat u(\xi',t)dt$ である。また $s>1/2$ なら
$$
A(\xi')=\int_\R (a^2+4\pi^2t^2)^{-s}dt
=C_s a^{1-2s},\qquad
C_s=\frac1{2\pi}\int_\R(1+z^2)^{-s}dz .
$$
Cauchy–Schwarz により $\left|\int\widehat u\,dt\right|^2\le A(\xi')\int(a^2+4\pi^2t^2)^s|\widehat u|^2dt$ である。これに $a^{2s-1}$ を掛けて $\xi'$ を積分すると、連続性が従う。Schwartz 関数は稠密であるため、一意な連続延長となる。この式は一般の $u$ に対しても成立する。実際、法線周波数の積分は Cauchy–Schwarz によってほとんど全ての $\xi'$ で絶対収束し、その結果を $H^{s-1/2}$ の変換とみなす操作自体が上と同じ評価を持つ。

任意の $g\in H^{s-1/2}(\R^{n-1})$ に対して $u$ を
$$
\widehat u(\xi',t)=\frac{\widehat g(\xi')}{A(\xi')}(a^2+4\pi^2t^2)^{-s}
$$
で定める。法線方向を積分すると $\widehat g$ に戻り、$H^s$ ノルムの二乗は $\int |\widehat g|^2/A=C_s^{-1}\|g\|_{H^{s-1/2}}^2$ である。従ってこれが連続な右逆となる。さらにこの同じ式は Cauchy–Schwarz の等号を達成するので、同じ跡を持つ延長のうち $H^s$ ノルムが最小である。

$s\le1/2$ では非零の $g$ で $\widehat g\in C_c^\infty$ を取り、$0\le\chi_R\le1$ が $[2,R]$ で $1$、$[1,2R]$ の外でゼロになる滑らかな切断とする。$B_R=\int\chi_R(t)(1+4\pi^2t^2)^{-s}dt$ は無限大に向かう。$\widehat u_R(\xi',t)=\widehat g(\xi')\chi_R(t)(1+4\pi^2t^2)^{-s}/\sqrt{B_R}$ は滑らかなコンパクト台関数である。$\widehat g$ の台上で $a$ は有界なので、$(a^2+4\pi^2t^2)^s$ と $(1+4\pi^2t^2)^s$ は定数倍で比較できる。そのため $\chi_R^2\le\chi_R$ を用いて $\|u_R\|_{H^s}$ が一様に有界となる。一方 $\operatorname{Tr}u_R=\sqrt{B_R}g$ である。$v_R=u_R/\sqrt{B_R}$ と置けば $H^s$ でゼロへ収束しながら跡は常に $g$ なので、作用の連続性を持つ延長は存在しない。
:::

:::example[一階ノルムを最小にする延長]
$s=1$ では $C_1=1/2$ なので $A=1/(2a)$ である。従って最小ノルム延長の法線変数の逆変換は $\mathcal F_{x'}u(\xi',x_n)=\widehat g(\xi')e^{-a|x_n|}$ になる。実際 $e^{-a|t|}$ の変換は $2a/(a^2+4\pi^2\xi_n^2)$ である。$\|u\|_{H^1}^2=2\|g\|_{H^{1/2}}^2$ であり、各半空間で $x_n$ の符号を固定して二回微分すれば $(1-\Delta)u=0$ となる。境界から離れると法線方向に指数減衰し、接方向の周波数が高いほど減衰は速い。
:::

### 臨界指数近くでの連続度

Hölder の評価を得る基本証明では指数関数の差を全周波数で一つの冪により評価した。低周波でその差が小さく、高周波で有界であることを別々に使えば、$s-k-n/2$ が $0$ と $1$ の間にあるとき、境界の指数そのものを得られる。

:::proposition[周波数の分割による連続度]{#prop-four-sobolev-holder-sharp}
整数 $k\ge0$、$\gamma=s-k-n/2$ とし、$0<\gamma\le1$ とする。$u\in H^s$、$|\alpha|=k$、$0<|h|\le1/2$ のとき
$$
|\partial^\alpha u(x+h)-\partial^\alpha u(x)|
\le C\|u\|_{H^s}
\begin{cases}
|h|^\gamma &(0<\gamma<1),\\
|h|\sqrt{1+\log(1/|h|)} &(\gamma=1)
\end{cases}
$$
が全ての $x$ で成り立つ。
:::

:::proof
連続代表元は :ref[thm-four-sobolev-embedding] によって存在する。逆積分の差に Cauchy–Schwarz を使うと、$\|u\|_{H^s}$ に掛ける因子の二乗は
$$
C\int_0^\infty \min(|h|^2r^2,1)\,r^{n-1+2k}(1+r^2)^{-s}dr
$$
以下である。$r<1$ では $C|h|^2$ 以下、$1\le r\le |h|^{-1}$ では $C|h|^2\int_1^{|h|^{-1}}r^{1-2\gamma}dr$ 以下、$r>|h|^{-1}$ では $C\int_{|h|^{-1}}^\infty r^{-1-2\gamma}dr$ 以下である。中間の積分は $0<\gamma<1$ なら $C|h|^{2\gamma}$、$\gamma=1$ なら $C|h|^2\log(1/|h|)$ を与える。最後の積分は $C|h|^{2\gamma}$ である。各評価を足して平方根を取ると表示式が得られる。
:::

:::example[一次元で半階の Hölder 連続性を得る]
$n=1,k=0,s=1$ では $\gamma=1/2$ である。従って $H^1(\R)$ の連続代表元は $1/2$-Hölder 連続であり、同じ結論は $u'\in L^2$ から $|\int_x^{x+h}u'|\le |h|^{1/2}\|u'\|_2$ としても得られる。後者はまず滑らかな $u$ で示せ、$H^1$ 近似と連続代表元の一様収束で一般の元に移せる。二次元で $s=2,k=0$ なら $\gamma=1$ であり、周波数側の中間範囲の積分が対数になるので、同じ計算から直接 Lipschitz 評価を結論できない。
:::

## 緩増加超関数の構造定理

超関数の作用には試験関数の有限階の微分だけが必要である。これを関数の側へ移して、連続関数の有限個の微分として表す。無限遠の増大を先に重みで除くと、Sobolev 空間の双対が使える。

:::theorem[緩増加超関数の構造]{#thm-four-tempered-structure}
$T\in\mathcal S'(\R^n)$ は有限和 $\sum_{|\alpha|\le m}\partial^\alpha F_\alpha$ で表せる。各 $F_\alpha$ は連続で、多項式増大を持つ。逆に、この形の和は緩増加超関数である。
:::

:::proof
:eqref[eq-four-tempered-bound] の $m,N$ を取り、$W(x)=(1+|x|^2)^M$ を $2M\ge N$ となるように選ぶ。$U=W^{-1}T$ と置く。$W^{-1}$ の各導関数は $C_\alpha(1+|x|)^{-2M}$ 以下である。実際、微分するたび分母の冪が増え、分子に現れる座標の冪はそれを打ち消すほど増えない。Leibniz 則から

$$
|\langle U,\phi\rangle|\le C'\max_{|\alpha|\le m}\|\partial^\alpha\phi\|_\infty\le C''\|\phi\|_{H^s}
$$

となる。ただし $s>m+n/2$ とし、最後は :ref[thm-four-sobolev-embedding] を用いた。$J_s$ で滑らかなコンパクト台の周波数関数は Schwartz 関数の像となり $L^2$ で稠密なので、$U$ の作用は $H^s$ に連続に延長できる。前節の双対の同定より $U\in H^{-s}$ である。

整数 $l$ を $2l-s>n/2$ となるよう取り、$\widehat V=\langle\xi\rangle^{-2l}\widehat U$ とする。$V\in H^{2l-s}$ なので有界連続である。微分則で $U=(1-\Delta)^l V$ であり、$T=W(1-\Delta)^l V$ となる。積の微分を逆に使い $W\partial_j V=\partial_j(WV)-(\partial_j W)V$ と書く。この操作を反復すれば連続な多項式倍 $P V$ の有限階微分の和になる。逆方向は各 $F_\alpha$ の多項式増大と試験関数の急減少から :eqref[eq-four-tempered-bound] を直接得る。
:::

:::example[古典関数の微分で点を表す]
$\delta_0=\partial^2(|x|/2)$ は増大が一次の連続関数の二階微分である。同様に $\operatorname{pv}(1/x)=(\log|x|)'$ が成り立つ。後者では原点の両側で部分積分を行うと境界項が $\log\eps(\phi(\eps)-\phi(-\eps))$ となりゼロに近づく。連続な原始関数 $x\log|x|-x$ にさらに一階微分を移せば、構造定理の形になる。
:::

:::example[階段関数の連続な原始関数]
一次元で $x_+=\max(x,0)$ とすると、$\int x_+\phi'=-\int_0^\infty\phi$ である。従って $H=\partial_x x_+$ は多項式増大を持つ連続関数の一次微分である。さらに $\sgn x$ は $|x|$ の微分である。実際、正負の半軸で部分積分すれば境界項はどちらも $0$ であり、符号の違いだけが残る。
:::

### Bessel 核で一つの連続関数にまとめる

構造定理の表示では項ごとに異なる連続関数を使った。周波数の各成分を十分高い階数だけ減衰させれば、全体を一つの連続関数に移し、一定の微分作用素で戻すこともできる。その逆作用素は正の熱核の重ね合わせで表せる。

:::proposition[Bessel 核と一つの関数による表示]{#prop-four-tempered-bessel-structure}
整数 $l\ge1$ とし、熱核 $p_t(x)=(4\pi t)^{-n/2}e^{-|x|^2/(4t)}$ により
$$
G_l(x)=\frac1{\Gamma(l)}\int_0^\infty e^{-t}t^{l-1}p_t(x)dt
$$
と定める。これは非負の可積分関数で、$\widehat G_l(\xi)=\langle\xi\rangle^{-2l}$ である。また $T\in\mathcal S'$ なら、十分大きい整数 $l$ に対し $V=\mathcal J^{-2l}T$ は多項式増大の連続関数であり、$T=(1-\Delta)^lV$ である。
:::

:::proof
$\int p_t=1$ なので、Tonelli により $\int G_l=1$ である。また $\widehat p_t=e^{-4\pi^2t|\xi|^2}$ だから、Fubini 後に $u=t(1+4\pi^2|\xi|^2)$ と置くと $\widehat G_l=\langle\xi\rangle^{-2l}$ となる。

$|\alpha|\le m$ とし $l>m/2$ とする。このとき $\partial^\alpha G_l$ は積分内の $p_t$ を微分することで得られる関数であり、全ての $N\ge0$ について $\int(1+|x|)^N|\partial^\alpha G_l(x)|dx<\infty$ である。実際、$x=\sqrt t\,y$ と置くと
$$
\int(1+|x|)^N|\partial^\alpha p_t(x)|dx
\le C_{\alpha,N} t^{-|\alpha|/2}(1+t^{N/2}).
$$
$p_1$ の導関数は多項式と Gauss 関数の積なので、変数を移した後の積分定数は有限である。$t\le1$ では $l>|\alpha|/2$、$t\ge1$ では指数 $e^{-t}$ によって積分できる。これにより試験関数に対する積分で微分を移せるので、導関数の式も超関数として成り立つ。

:ref[thm-four-tempered-structure] によって $T=\sum_{|\alpha|\le m}\partial^\alpha F_\alpha$ とし、$V=\sum_{|\alpha|\le m}(\partial^\alpha G_l)*F_\alpha$ と置く。各積分は絶対収束する。$|F_\alpha(x-y)|\le C(1+|x|)^N(1+|y|)^N$ と評価して、上の重み付き可積分性を使えば $V$ は多項式増大する。また $x$ を有界集合に制限した場合は同じ可積分な支配関数が使え、$F_\alpha$ が連続なので優収束により $V$ は連続である。

最後に、全ての重み付き $L^1$ ノルムが有限な核 $K$ について、$\phi\in\mathcal S$ なら $\int K(y)\phi(\cdot+y)dy$ は Schwartz 空間に属する。各導関数に $(1+|x|)^N\le(1+|x+y|)^N(1+|y|)^N$ を使えばよい。この操作と Fourier 変換を試験関数上で交換すると $\widehat K$ の乗算になる。例えば変換済みの試験関数を $\widehat\phi(x+y)=\int\phi(\xi)e^{-2\pi i(x+y)\cdot\xi}d\xi$ と展開し、絶対収束で $y,\xi$ 積分を交換できる。従って多項式増大の連続関数 $F$ との畳み込みについて $\widehat{K*F}=\widehat K\,\widehat F$ が成り立つ。これを $K=\partial^\alpha G_l$ に適用し、全項を足すと $\widehat V=\langle\xi\rangle^{-2l}\widehat T$ である。ゆえに $V=\mathcal J^{-2l}T$ であり、$(1-\Delta)^lV=T$ を得る。
:::

:::example[一次元の二回の平滑化を計算する]
一次元で $G_1(x)=e^{-|x|}/2$ であり、$G_2=G_1*G_1$ である。後者は変換が $(1+4\pi^2\xi^2)^{-2}$ となることと一意性から分かる。$x\ge0$ で積分を $y<0$、$0\le y\le x$、$y>x$ に分けると、係数 $1/4$ を掛ける前の三つの値は $e^{-x}/2,xe^{-x},e^{-x}/2$ である。従って $G_2(x)=(1+|x|)e^{-|x|}/4$ となる。例えば $T=\delta_0'$ では $V=G_2'=-xe^{-|x|}/4$ が連続であり、$(1-\partial_x^2)^2 V=\delta_0'$ である。デルタの導関数を一つの連続関数と一定の四階作用素で表せた。
:::

## 演習問題

::::exercise{level=basic}
$\mathcal F(\delta_a')$ を一次元で求めよ。
:::solution
$\widehat{\delta_a}=e^{-2\pi ia\xi}$ なので微分則により $2\pi i\xi e^{-2\pi ia\xi}$ である。
:::
::::

::::exercise{level=basic}
$\widehat{x^2}$ を求めよ。
:::solution
$\widehat1=\delta_0$ に座標倍の式を二回使う。係数は $(-2\pi i)^{-2}=-1/(4\pi^2)$ なので $-\delta_0''/(4\pi^2)$ である。
:::
::::

::::exercise{level=basic}
$H(x)-H(-x)$ の変換を求めよ。
:::solution
これは a.e. で $\sgn x$ である。反射で $\delta_0$ の項は同じ、主値の項は符号が反転する。差は $(\pi i)^{-1}\operatorname{pv}(1/\xi)$ である。
:::
::::

::::exercise{level=basic}
一次元で $\|e^{-\pi x^2}\|_{H^0}^2$ を求めよ。
:::solution
$H^0=L^2$ は Plancherel から従う。$\int e^{-2\pi x^2}dx=1/\sqrt2$ である。
:::
::::

::::exercise{level=advanced}
$\partial^\alpha:H^s\to H^{s-|\alpha|}$ の作用素ノルムが $1$ 以下であることを示せ。
:::solution
$|(2\pi\xi)^\alpha|\le\langle\xi\rangle^{|\alpha|}$ である。これを Fourier 側のノルム式に代入すると重みの積は $\langle\xi\rangle^s$ 以下となり、主張が従う。
:::
::::

::::exercise{level=advanced}
一次元で $(1-\partial^2)u=\delta_0$ を満たす緩増加超関数を求めよ。
:::solution
変換すると $(1+4\pi^2\xi^2)\widehat u=1$ である。乗数の逆数とその導関数は多項式増大以下なので $\mathcal S'$ の積として割れる。$\widehat{e^{-|x|}}=2/(1+4\pi^2\xi^2)$ だから $u=e^{-|x|}/2$ が一意解である。
:::
::::

::::exercise{level=advanced}
$\delta_0^{(k)}\in H^s(\R)$ の条件を求めよ。
:::solution
変換の絶対値は $(2\pi|\xi|)^k$ であり、二乗ノルムの被積分関数は無限遠で $|\xi|^{2k+2s}$ と同程度である。原点の積分は有限なので、条件は $2k+2s<-1$、すなわち $s<-k-1/2$ である。
:::
::::

::::exercise{level=advanced}
$T=\sum_{j\in\Z} j\delta_j$ は緩増加超関数であることを示し、その変換を求めよ。
:::solution
$\sum_j |j|(1+|j|)^{-3}<\infty$ なので作用は Schwartz セミノルムで抑えられる。$T=x\sum_j\delta_j$ と書け、格子デルタの変換は自身である。座標倍の式により $\widehat T=(i/(2\pi))\sum_j\delta_j'$ である。微分は連続な試験関数操作なので和に移せる。
:::
::::

::::exercise{level=basic}
一次元で $\sum_{k\in\Z}\delta_{3k+1}$ の Fourier 変換を求めよ。
:::solution
:ref[prop-four-tempered-lattice-comb] で $A=3,v=1$ なので、答えは $3^{-1}\sum_{m\in\Z}e^{-2\pi im/3}\delta_{m/3}$ である。双対格子の間隔は $1/3$、一セルあたりの密度も $1/3$ である。例えば $m=0,3$ の位相は $1$ だが、$m=1$ では $e^{-2\pi i/3}$ である。
:::
::::

::::exercise{level=advanced}
一次元で $G_2*1$ と $(\partial_x G_2)*1$ を求めよ。畳み込みは通常の絶対収束積分で定める。
:::solution
$G_2$ は質量 $1$ なので $G_2*1=1$ であり、$\partial_xG_2=-xe^{-|x|}/4$ は可積分な奇関数なので $(\partial_xG_2)*1=0$ である。周波数側でも $\widehat 1=\delta_0$ に $G_2$ の変換を掛けると原点の値 $1$ が残り、微分した核では原点で $2\pi i\xi$ がゼロになることに対応する。
:::
::::

## 確認テスト

::::quiz{answer="0.5" tol="1e-6"}
一次元で $\|\delta_0\|_{H^{-1}}^2$ を求めよ。
:::explanation
$\widehat{\delta_0}=1$ なので $u=2\pi\xi$ と置くと $(2\pi)^{-1}\int(1+u^2)^{-1}du=\pi/(2\pi)=1/2$ となる。
:::
::::

::::quiz
$\widehat H$ に含まれる点質量の係数はどれか。

- [ ] $0$
- [x] $1/2$
- [ ] $1$

:::explanation
$H(x)+H(-x)=1$ である。変換の反射を足すと主値部分が消え、点質量は同じ係数で二つ足される。$\widehat 1=\delta_0$ だから各係数は $1/2$ である。
:::
::::

::::quiz{answer="1" tol="1e-6"}
$n=2$ で :ref[thm-four-sobolev-embedding] が $k=0$ の連続性を与える $s$ の下限を答えよ。
:::explanation
条件は $s>n/2$ である。$n=2$ なら $s>1$ なので、許容範囲の下限は $1$ である。
:::
::::

::::quiz
$\R^n$ で次数 $d$ の斉次超関数の変換の次数はどれか。

- [ ] $-d$
- [x] $-n-d$
- [ ] $d+n$

:::explanation
$\mathcal F(T(\lambda\cdot))=\lambda^{-n}\widehat T(\cdot/\lambda)$ の左辺は $\lambda^d\widehat T$ である。従って $\widehat T(\cdot/\lambda)=\lambda^{n+d}\widehat T$ となり、$1/\lambda$ を尺度と読めば次数は $-n-d$ である。
:::
::::

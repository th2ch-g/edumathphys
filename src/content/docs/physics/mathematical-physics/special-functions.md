---
title: "ガンマ関数と直交多項式"
description: "Gamma と Beta の積分を導き、直交多項式の漸化式・零点・求積を共通の仕組みで計算する。"
notation:
  - symbol: '\Gamma(z)'
    meaning: "Gamma 関数"
  - symbol: '\psi(z)'
    meaning: "Gamma 関数の対数微分"
  - symbol: '\zeta(s)'
    meaning: "Riemann ゼータ関数"
  - symbol: 'H_n,L_n^{(\alpha)},T_n'
    meaning: "Hermite・一般化 Laguerre・第 1 種 Chebyshev 多項式"
---

## この章の目標

- Gamma、Beta およびその対数微分を積分と級数から計算する。
- 三項漸化式と零点を直交性から導く。
- Hermite、Laguerre、Chebyshev の規格化を選んで展開する。
- 多項式の零点から Gauss 求積の重みを構成する。

## ガンマ関数とベータ関数

### 階乗を連続な引数へ延ばす

$\Re z>0$ の :term[ガンマ関数]{en="gamma function" yomi="がんまかんすう"} は $\Gamma(z)=\int_0^\infty t^{z-1}e^{-t}\dd t$ である。$t>0$ の実対数を用いる。端点を支配する積分を使えば局所一様収束と正則性が従う。部分積分で $\Gamma(z+1)=z\Gamma(z)$、直接 $\Gamma(1)=1$ であるから $\Gamma(n+1)=n!$ となる。漸化式で左へ延長すると $-n$（$n\ge0$）に単純極を持ち、$\Gamma(z+n+1)/(z(z+1)\cdots(z+n))$ をそこで評価して留数 $(-1)^n/n!$ を得る。

:term[ベータ関数]{en="beta function" yomi="べーたかんすう"} を $B(x,y)=\int_0^1 u^{x-1}(1-u)^{y-1}\dd u$ とする（$\Re x,\Re y>0$）。$\Gamma(x)\Gamma(y)$ の二重積分で $(s,t)=(vu,v(1-u))$ と置くと Jacobian は $v$、指数は $e^{-v}$ になる。絶対収束するので、先に $u$、次に $v$ で積分して $\Gamma(x)\Gamma(y)=B(x,y)\Gamma(x+y)$ を得る。

$x=z$、$y=1-z$ では右辺の Gamma は $\Gamma(1)=1$ である。$B(z,1-z)$ に $u=t/(1+t)$ を入れ、式 :eqref[eq-mp-keyhole-beta] を使うと $\Gamma(z)\Gamma(1-z)=\pi/\sin\pi z$ を得る。この反射公式はまず $0<\Re z<1$ で成り立ち、その後に有理型関数の恒等式として接続できる。右辺は非整数点で有限かつ非零で、左辺のどちらの因子にもそこで極はない。従って Gamma 関数は非整数点で零にならない。正整数での値も階乗で非零なので、先ほどの積の式を割って
$$
 B(x,y)=\frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)}
 \label{eq-mp-beta-gamma}
$$
を得る。例えば反射公式に $z=1/2$ を入れ、積分から正値を選ぶと $\Gamma(1/2)=\sqrt\pi$ である。漸化式により $\Gamma(5/2)=3\sqrt\pi/4$ になる。さらに $B(z,z)$ に $u=(1+v)/2$、続いて $v^2=w$ を使うと $B(z,z)=2^{1-2z}B(1/2,z)$ となる。Gamma との関係に戻せば倍角公式 $\Gamma(z)\Gamma(z+1/2)=2^{1-2z}\sqrt\pi\,\Gamma(2z)$ が得られる。

### 正値と対数凸性による特徴付け

実変数 $x>0$ では $\Gamma(x)>0$ である。積分表示と Hölder の不等式から、$0<t<1$ について $\Gamma((1-t)x+ty)\le\Gamma(x)^{1-t}\Gamma(y)^t$ となる。従って $\log\Gamma$ は凸である。逆に $f>0$、$f(1)=1$、$f(x+1)=xf(x)$、$\log f$ が凸という条件を課すと $f=\Gamma$ である。実際、$0<x\le1$、整数 $N\ge2$ について、$\log f$ の区間 $[N,N+x]$ の割線の傾きは $[N-1,N]$ と $[N,N+1]$ の傾きで挟まれる。従って
$$
 (N-1)^x(N-1)!\le f(N+x)\le N^x(N-1)! .
$$
両端を $x(x+1)\cdots(x+N-1)$ で割れば、真ん中は $f(x)$ になる。同じ上下界が $\Gamma(x)$ にも成り立ち、上下界の比が $(N/(N-1))^x\to1$ となるので $f(x)=\Gamma(x)$ である。漸化式によって全ての正の実数へ延長できる。これを :term[ボーア・モレルップの定理]{en="Bohr–Mollerup theorem" yomi="ぼーあもれるっぷのていり"} という。例えば $\Gamma(x)$ に正の周期 1 の因子を掛けると漸化式は保たれるが、規格化と対数凸性を同時に保つ因子は定数 $1$ に限られる。

### 引数を半分ずらす計算

Euler の積分から $\Gamma(z+1)=z\Gamma(z)$ が従う。正の半整数の初期値は $t=u^2$ と置けば $2\int_0^\infty e^{-u^2}\dd u=\sqrt{\pi}$ になる。この値を一段ずつ上げると整数の場合とは違う係数列を得る。逆に負の半整数では再帰で割り算を使う。分母が零になる整数点では、この手順を点の値として行えない。そこで局所的に再帰を使うと、$\Gamma$ の $-n$ における留数は $(-1)^n/n!$ になる。

:::example[半整数の Gamma と Beta（1）]
半整数の積分は、平方根の因子を除いた有理数を計算すればよい。

$$
m=1
$$

半整数の値を再帰で求める。最後は $\Gamma(1/2)=\sqrt{\pi}$ である。

$$
\begin{aligned}
\Gamma(m+1/2)\big|_{m=1} &= \frac{\sqrt{\pi}}{2}
\end{aligned}
$$

次に二つの引数を持つ Beta を Gamma の比で表す。分母の引数は二つの和である。

$$
\begin{aligned}
B(m+1/2,1/2)\big|_{m=1} &= \frac{\pi}{2}
\end{aligned}
$$

この値は $t=\sin^2\theta$ で円弧上のモーメントになる。Jacobian から $2$ が出る。

$$
\begin{aligned}
\int_0^{\pi/2}\sin^{2m}\theta\,\dd\theta\big|_{m=1} &= \frac{\pi}{4}
\end{aligned}
$$
:::

:::example[半整数の Gamma と Beta（2）]
半整数の積分は、平方根の因子を除いた有理数を計算すればよい。次数を一つ上げる。

$$
m=2
$$

半整数の値を再帰で求める。最後は $\Gamma(1/2)=\sqrt{\pi}$ である。

$$
\begin{aligned}
\Gamma(m+1/2)\big|_{m=2} &= \frac{3 \sqrt{\pi}}{4}
\end{aligned}
$$

次に二つの引数を持つ Beta を Gamma の比で表す。分母の引数は二つの和である。

$$
\begin{aligned}
B(m+1/2,1/2)\big|_{m=2} &= \frac{3 \pi}{8}
\end{aligned}
$$

この値は $t=\sin^2\theta$ で円弧上のモーメントになる。Jacobian から $2$ が出る。

$$
\begin{aligned}
\int_0^{\pi/2}\sin^{2m}\theta\,\dd\theta\big|_{m=2} &= \frac{3 \pi}{16}
\end{aligned}
$$
:::

:::example[半整数の Gamma と Beta（3）]
半整数の積分は、平方根の因子を除いた有理数を計算すればよい。さらに次数を上げて、低次の結果と照合する。

$$
m=3
$$

半整数の値を再帰で求める。最後は $\Gamma(1/2)=\sqrt{\pi}$ である。

$$
\begin{aligned}
\Gamma(m+1/2)\big|_{m=3} &= \frac{15 \sqrt{\pi}}{8}
\end{aligned}
$$

次に二つの引数を持つ Beta を Gamma の比で表す。分母の引数は二つの和である。

$$
\begin{aligned}
B(m+1/2,1/2)\big|_{m=3} &= \frac{5 \pi}{16}
\end{aligned}
$$

この値は $t=\sin^2\theta$ で円弧上のモーメントになる。Jacobian から $2$ が出る。

$$
\begin{aligned}
\int_0^{\pi/2}\sin^{2m}\theta\,\dd\theta\big|_{m=3} &= \frac{5 \pi}{32}
\end{aligned}
$$
:::

:::example[半整数の Gamma と Beta（4）]
半整数の積分は、平方根の因子を除いた有理数を計算すればよい。四次の結果まで具体化する。

$$
m=4
$$

半整数の値を再帰で求める。最後は $\Gamma(1/2)=\sqrt{\pi}$ である。

$$
\begin{aligned}
\Gamma(m+1/2)\big|_{m=4} &= \frac{105 \sqrt{\pi}}{16}
\end{aligned}
$$

次に二つの引数を持つ Beta を Gamma の比で表す。分母の引数は二つの和である。

$$
\begin{aligned}
B(m+1/2,1/2)\big|_{m=4} &= \frac{35 \pi}{128}
\end{aligned}
$$

この値は $t=\sin^2\theta$ で円弧上のモーメントになる。Jacobian から $2$ が出る。

$$
\begin{aligned}
\int_0^{\pi/2}\sin^{2m}\theta\,\dd\theta\big|_{m=4} &= \frac{35 \pi}{256}
\end{aligned}
$$
:::

## ディガンマ関数と Riemann ゼータ関数

### 積の対数微分

$\psi(z)=\Gamma'(z)/\Gamma(z)$ （$z\notin\{0,-1,-2,\ldots\}$）を :term[ディガンマ関数]{en="digamma function" yomi="でぃがんまかんすう"} という。Euler の極限式は $N^z B(z,N+1)\to\Gamma(z)$ から得られる。実際、$u=v/N$ と置くと $(1-v/N)^N\mathbf1_{0<v<N}\le e^{-v}$ で、積分と必要な $z$ 微分に優収束を使える。従って $\Re z>0$ で
$$
\Gamma(z)=\lim_{N\to\infty}\frac{N^zN!}{z(z+1)\cdots(z+N)},\qquad
\psi(z)=-\gamma+\sum_{n=0}^\infty\left(\frac1{n+1}-\frac1{n+z}\right).
$$
右辺は極を避けるコンパクト集合で局所一様に収束するので、ディガンマの式はその全域に延びる。ここで $\gamma=\lim_N(\sum_{n=1}^N1/n-\log N)$ である。この極限は積分との比較で単調有界なことから存在する。例えば $\psi(3)=3/2-\gamma$ である。項別微分すると $\psi^{(m)}(1)=(-1)^{m+1}m!\sum_{n\ge1}n^{-m-1}$（$m\ge1$）となる。

:term[リーマンゼータ関数]{en="Riemann zeta function" yomi="りーまんぜーたかんすう"} を $\Re s>1$ で $\zeta(s)=\sum_{n\ge1}n^{-s}$ と定める。この範囲では級数は絶対収束し、先ほどの導関数の値は $\psi^{(m)}(1)=(-1)^{m+1}m!\zeta(m+1)$ と書ける。そこで $\Gamma(s)\zeta(s)=\int_0^\infty t^{s-1}/(e^t-1)\dd t$ である。$1/(e^t-1)=\sum_{n\ge1}e^{-nt}$ と展開して積分し、$v=nt$ を使うと示せる。$0<t<1$ で $1/(e^t-1)=t^{-1}-1/2+O(t)$ を引き、引いた項を積分すれば $\Gamma(s)\zeta(s)$ の表示に $1/(s-1)-1/(2s)$ が現れる。残りは $\Re s>-1$ で正則なので $s=1$ で留数 $1$、$s=0$ で $\zeta(0)=-1/2$ を得る。さらに Taylor 項を引くと任意に左へ延長でき、$s=1$ 以外の分子の極は $\Gamma$ の極で相殺される。

### 逆 Gamma の積表示

Euler の極限式を使うと
$$
\frac1{\Gamma(z)}=z e^{\gamma z}\prod_{n=1}^\infty
\left(1+\frac zn\right)e^{-z/n}
$$
を得る。有限積の両辺を比較すると余分な指数は $z(\sum_{n=1}^N1/n-\log N)$ で、その極限が $z\gamma$ である。積の第 $n$ 因子の対数は、固定したコンパクト集合で十分大きい $n$ に対し $O(n^{-2})$ である。従って末尾の積は局所一様収束し零にならない。有限個の因子だけが零点を与えるため、$1/\Gamma(z)$ は $0,-1,-2,\ldots$ にだけ単純零点を持つ。正の引数では積の各因子が正、負の非整数では有限個の負因子が符号を決める。例えば $\Gamma(-1/2)=-2\sqrt\pi$ は漸化式でも積の符号でも確認できる。

### Riemann ゼータ関数の関数等式

Poisson の和公式 :eqref[eq-mp-poisson-sum] を Gaussian に適用する。$\vartheta(t)=\sum_{n\in\Z}e^{-\pi n^2t}$ と書くと、その式に $T=\sqrt t$ を入れて $\vartheta(t)=t^{-1/2}\vartheta(1/t)$ となる。$\Re s>1$ では項別積分が絶対収束し、
$$
\Lambda(s)=\pi^{-s/2}\Gamma(s/2)\zeta(s)
=\frac12\int_0^\infty(\vartheta(t)-1)t^{s/2-1}\dd t
$$
となる。$t<1$ で変換公式を入れ、$t=1/u$ と置く。定数項を積分すると、全体は
$$
\Lambda(s)=\frac1{s-1}-\frac1s+
\frac12\int_1^\infty(\vartheta(t)-1)
 \bigl(t^{s/2-1}+t^{(1-s)/2-1}\bigr)\dd t
$$
である。積分は指数減衰により $s$ の整関数である。右辺が $s\mapsto1-s$ で不変なので、この式が有理型接続と関数等式 $\Lambda(s)=\Lambda(1-s)$ を同時に与える。$s=-2,-4,\ldots$ で $\Gamma(s/2)$ が単純極、$\Lambda(s)=\Lambda(1-s)$ が有限で非零だから、$\zeta$ はこれらの点で単純零点を持つ。例えば $\zeta(-2)=0$ である。$s=0$ は $\Lambda$ 自体も極であるため同じ議論の対象外であり、上で得た $-1/2$ になる。

### ゼータ関数の Euler 積とポリガンマ関数

$\Re s>1$ において
$$
\zeta(s)=\prod_{p\ {\rm prime}}(1-p^{-s})^{-1}
$$
である。有限個の素数について幾何級数を掛ければ、その素数だけを素因数に持つ整数にわたる和になる。素因数分解は一意で、$n^{-\Re s}$ の和が有限なので、素数をすべて含む極限で $\zeta(s)$ を得る。また $\sum_p|p^{-s}|<\infty$ であるため、有限個を除いた積の対数は絶対収束する。各因子も非零なので、$\Re s>1$ で $\zeta(s)$ に零点はない。

$m\ge1$ に対し $\psi^{(m)}$ を :term[ポリガンマ関数]{en="polygamma function" yomi="ぽりがんまかんすう"} という。ディガンマ級数を $m$ 回微分すると
$$
\psi^{(m)}(z)=(-1)^{m+1}m!\sum_{n\ge0}(n+z)^{-m-1}
$$
となる。例えば $\psi^{(1)}(1/2)=4\sum_{n\ge0}(2n+1)^{-2}=3\zeta(2)$ である。さらに反射公式を対数微分すると $\psi(1-z)-\psi(z)=\pi\cot\pi z$ となる。$z=1/4$ では $\psi(3/4)-\psi(1/4)=\pi$、倍角公式を微分すれば和は $2\psi(1/2)-2\log 2$ である。従って両者を消去すると $\psi(1/4)=-\gamma-3\log 2-\pi/2$ になる。

### 積を微分して調和和を取り出す

再帰の対数微分から $\psi(z+1)=\psi(z)+1/z$、もう一度微分すると $\psi^\prime(z+1)=\psi^\prime(z)-1/z^2$ を得る。初期値 $\psi(1)=-\gamma$ は Weierstrass 積で $z=1$ を代入して求まる。$\psi^\prime(1)=\zeta(2)$ では対数微分の部分和が単調に収束し、項別微分で二乗の逆数和になる。従って整数での評価は有限の調和和と原点側の規格化に分けられる。

:::example[整数での対数微分（1）]
初期値から再帰して、一次と二次の逆数和を比べる。

$$
m=1
$$

一次の調和数を足す。下の添字は $m$ であるので引数は $m+1$ とする。

$$
\begin{aligned}
\psi(m+1)+\gamma\big|_{m=1} &= 1
\end{aligned}
$$

二階対数微分の差は、逆二乗和の符号を反転したものである。

$$
\begin{aligned}
\psi^\prime(m+1)-\zeta(2)\big|_{m=1} &= -1
\end{aligned}
$$

無限和からこの有限和を引くと、残る値は正である。

$$
\begin{aligned}
\psi^\prime(m+1)\big|_{m=1} &= -1 + \frac{\pi^{2}}{6}
\end{aligned}
$$
:::

:::example[整数での対数微分（2）]
初期値から再帰して、一次と二次の逆数和を比べる。次数を一つ上げる。

$$
m=2
$$

一次の調和数を足す。下の添字は $m$ であるので引数は $m+1$ とする。

$$
\begin{aligned}
\psi(m+1)+\gamma\big|_{m=2} &= \frac{3}{2}
\end{aligned}
$$

二階対数微分の差は、逆二乗和の符号を反転したものである。

$$
\begin{aligned}
\psi^\prime(m+1)-\zeta(2)\big|_{m=2} &= - \frac{5}{4}
\end{aligned}
$$

無限和からこの有限和を引くと、残る値は正である。

$$
\begin{aligned}
\psi^\prime(m+1)\big|_{m=2} &= - \frac{5}{4} + \frac{\pi^{2}}{6}
\end{aligned}
$$
:::

:::example[整数での対数微分（3）]
初期値から再帰して、一次と二次の逆数和を比べる。さらに次数を上げて、低次の結果と照合する。

$$
m=3
$$

一次の調和数を足す。下の添字は $m$ であるので引数は $m+1$ とする。

$$
\begin{aligned}
\psi(m+1)+\gamma\big|_{m=3} &= \frac{11}{6}
\end{aligned}
$$

二階対数微分の差は、逆二乗和の符号を反転したものである。

$$
\begin{aligned}
\psi^\prime(m+1)-\zeta(2)\big|_{m=3} &= - \frac{49}{36}
\end{aligned}
$$

無限和からこの有限和を引くと、残る値は正である。

$$
\begin{aligned}
\psi^\prime(m+1)\big|_{m=3} &= - \frac{49}{36} + \frac{\pi^{2}}{6}
\end{aligned}
$$
:::

:::example[整数での対数微分（4）]
初期値から再帰して、一次と二次の逆数和を比べる。四次の結果まで具体化する。

$$
m=4
$$

一次の調和数を足す。下の添字は $m$ であるので引数は $m+1$ とする。

$$
\begin{aligned}
\psi(m+1)+\gamma\big|_{m=4} &= \frac{25}{12}
\end{aligned}
$$

二階対数微分の差は、逆二乗和の符号を反転したものである。

$$
\begin{aligned}
\psi^\prime(m+1)-\zeta(2)\big|_{m=4} &= - \frac{205}{144}
\end{aligned}
$$

無限和からこの有限和を引くと、残る値は正である。

$$
\begin{aligned}
\psi^\prime(m+1)\big|_{m=4} &= - \frac{205}{144} + \frac{\pi^{2}}{6}
\end{aligned}
$$
:::

## 直交多項式の一般論（三項漸化式と零点）

### 直交性を漸化式に変える

区間に正の重み $w$ を置き、全モーメントが有限で $w>0$ がほとんど至るところで成り立つとする。内積 $\int f\overline g w$ に対して $1,x,\ldots$ を直交化し、先頭係数 $1$ の $P_n$ を取る。$h_n=\int P_n^2w>0$ とする。$xP_n$ の $P_j$ 係数は $j<n-1$ なら零である。$\int xP_nP_jw=\int P_n(xP_j)w$ であり $xP_j$ の次数が $n$ 未満だからである。$j=n-1$ なら $xP_{n-1}$ の先頭成分は $P_n$ だから係数は $h_n/h_{n-1}$ である。従って
$$
xP_n=P_{n+1}+\alpha_nP_n+\frac{h_n}{h_{n-1}}P_{n-1},\qquad
\alpha_n=h_n^{-1}\int xP_n^2w
\label{eq-mp-orthogonal-recurrence}
$$
を得る。$n=0$ で最後の項は省く。

$P_n$ が区間の内部で符号を変える相異なる点を全て集め、その点を零点とするモニック多項式 $S$ を作る。すると $P_nS$ は区間内で一定符号で、ほとんど至るところで非零である。もし $\deg S<n$ なら直交性で $\int P_nSw=0$ となり矛盾する。従って零点は全て実で内部にあり、単純である。さらに漸化式で $P_n$ の根における $P_{n+1}$ の符号は $-P_{n-1}$ の符号である。帰納的に $P_{n-1}$ の根が $P_n$ の隙間に一つずつあれば、その符号は隣接する根で交代し、両端にも先頭係数の符号と比較して一つずつ $P_{n+1}$ の根が入る。$P_0=1$ から始めれば根の交錯が示せる。

:::example[直交化の最初の三項]
区間 $[-1,1]$ で $w=1$ とする。偶奇で直交するから $P_0=1,P_1=x$ である。$x^2$ の定数への射影は $(2/3)/2=1/3$ なので $P_2=x^2-1/3$ である。ここで $\alpha_1=0$、$h_1/h_0=(2/3)/2=1/3$ であり、式 :eqref[eq-mp-orthogonal-recurrence] の $xP_1=P_2+P_0/3$ と一致する。$P_2$ の二根 $\pm1/\sqrt3$ の間に $P_1$ の根 $0$ が入る。
:::

### 三重対角行列と零点

$e_j=P_j/\sqrt{h_j}$ とおく。次数 $n-1$ 以下の空間で $x$ を掛けてから同じ空間に直交射影する作用素は、基底 $e_0,\ldots,e_{n-1}$ で実対称三重対角行列 $J_n$ になる。対角成分は $\alpha_j$、隣の成分は $\sqrt{h_{j+1}/h_j}$ である。行列式 $\det(xI-J_n)$ を最後の行と列で展開すると、$\det(xI-J_{n-1})$ に $x-\alpha_{n-1}$ を掛けた項から、$\det(xI-J_{n-2})$ に $h_{n-1}/h_{n-2}$ を掛けた項を引くことになる。これは $P_0=1$、$P_1=x-\alpha_0$ から始まる式 :eqref[eq-mp-orthogonal-recurrence] と同じなので、$\det(xI-J_n)=P_n(x)$ である。従って零点を求める問題は、対称行列の固有値問題として計算できる。先の交錯性は、$J_{n-1}$ が $J_n$ の主部分行列であることにも対応する。

単位重みの $[-1,1]$ では、この節のモニック列 $P_m$ とは別に $R_m(1)=1$ の規格化を取れる。実際、$R_m=(2^m m!)^{-1}\dvn{(x^2-1)^m}{x}{m}$ は $m$ 回の部分積分により低次と直交し、$x=1$ で微分すれば値が $1$ になる。[Legendre 多項式](/physics/mathematical-physics/legendre-harmonics/) では、この規格化を使う。ここで $R_m$ の先頭係数を $a_m=(2m)!/[2^m(m!)^2]$ とおくと、$P_m=R_m/a_m$ である。同じ部分積分を $\int R_m^2$ に用いると $2^{-m}a_m\int_{-1}^1 (1-x^2)^m\dd x$ となり、$x=2v-1$ と Beta 積分で評価して $2/(2m+1)$ を得る。従って以下の例では $h_m=2/[(2m+1)a_m^2]$ を使う。

:::example[Legendre 系のノルム比（1）]
重み $1$ の区間 $[-1,1]$ でモニックに直したとき、漸化式の下側の係数を計算する。

$$
m=1
$$

最高次係数を取り出し、これで割る。

$$
\begin{aligned}
a_m\big|_{m=1} &= 1
\end{aligned}
$$

規格化していないノルムからモニックのノルムに直す。

$$
\begin{aligned}
h_m\big|_{m=1} &= \frac{2}{3}
\end{aligned}
$$

前の次数との比を取る。

$$
\begin{aligned}
h_m/h_{m-1}\big|_{m=1} &= \frac{1}{3}
\end{aligned}
$$
:::

:::example[Legendre 系のノルム比（2）]
重み $1$ の区間 $[-1,1]$ でモニックに直したとき、漸化式の下側の係数を計算する。次数を一つ上げる。

$$
m=2
$$

最高次係数を取り出し、これで割る。

$$
\begin{aligned}
a_m\big|_{m=2} &= \frac{3}{2}
\end{aligned}
$$

規格化していないノルムからモニックのノルムに直す。

$$
\begin{aligned}
h_m\big|_{m=2} &= \frac{8}{45}
\end{aligned}
$$

前の次数との比を取る。

$$
\begin{aligned}
h_m/h_{m-1}\big|_{m=2} &= \frac{4}{15}
\end{aligned}
$$
:::

:::example[Legendre 系のノルム比（3）]
重み $1$ の区間 $[-1,1]$ でモニックに直したとき、漸化式の下側の係数を計算する。さらに次数を上げて、低次の結果と照合する。

$$
m=3
$$

最高次係数を取り出し、これで割る。

$$
\begin{aligned}
a_m\big|_{m=3} &= \frac{5}{2}
\end{aligned}
$$

規格化していないノルムからモニックのノルムに直す。

$$
\begin{aligned}
h_m\big|_{m=3} &= \frac{8}{175}
\end{aligned}
$$

前の次数との比を取る。

$$
\begin{aligned}
h_m/h_{m-1}\big|_{m=3} &= \frac{9}{35}
\end{aligned}
$$
:::

:::example[Legendre 系のノルム比（4）]
重み $1$ の区間 $[-1,1]$ でモニックに直したとき、漸化式の下側の係数を計算する。四次の結果まで具体化する。

$$
m=4
$$

最高次係数を取り出し、これで割る。

$$
\begin{aligned}
a_m\big|_{m=4} &= \frac{35}{8}
\end{aligned}
$$

規格化していないノルムからモニックのノルムに直す。

$$
\begin{aligned}
h_m\big|_{m=4} &= \frac{128}{11025}
\end{aligned}
$$

前の次数との比を取る。

$$
\begin{aligned}
h_m/h_{m-1}\big|_{m=4} &= \frac{16}{63}
\end{aligned}
$$
:::

:::example[Legendre 系への非多項式型の射影]
所定の重みで直交射影し、与えた関数のうち低次のモードがどこまで回収されるか調べる。

射影する関数を $f=|x|$ とし、基底には $R_j$ を使う。同じ区間と重みで各成分の積分を計算する。正規化済みでない基底も使用できる。

$$
c_j=\frac{\int_{-1}^{1} f\phi_j\,w\,\dd x}{\int_{-1}^{1}\phi_j^2w\,\dd x},\qquad w=1
$$

第 1 の方向に対し、分子と規格化のための分母を別々に計算する。

$$
\begin{aligned}
\phi_{0} &= 1 \\
\int f\phi_{0}w\,\dd x &= 1 \\
\int \phi_{0}^2w\,\dd x &= 2 \\
c_{0} &= \frac{1}{2}
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
\int f\phi_{2}w\,\dd x &= \frac{1}{4} \\
\int \phi_{2}^2w\,\dd x &= \frac{2}{5} \\
c_{2} &= \frac{5}{8}
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
\int f\phi_{4}w\,\dd x &= - \frac{1}{24} \\
\int \phi_{4}^2w\,\dd x &= \frac{2}{9} \\
c_{4} &= - \frac{3}{16}
\end{aligned}
$$

成分を足すと、指定された部分空間への射影が得られる。残差の内積は、使った各基底方向で零になる。

$$
\Pi f=- \frac{105 x^{4}}{128} + \frac{105 x^{2}}{64} + \frac{15}{128},\qquad f-\Pi f=\frac{105 x^{4}}{128} - \frac{105 x^{2}}{64} + \left|{x}\right| - \frac{15}{128}
$$
:::

### Christoffel–Darboux の公式

重み付き空間で次数 $0$ から $n$ までの多項式へ直交射影する積分核は
$$
K_n(x,y)=\sum_{j=0}^{n}\frac{P_j(x)P_j(y)}{h_j}
 =\frac{P_{n+1}(x)P_n(y)-P_n(x)P_{n+1}(y)}{h_n(x-y)}
 \quad (x\ne y)
$$
である。実際、和の各項に $x-y$ を掛けて式 :eqref[eq-mp-orthogonal-recurrence] を両変数に使うと、$\alpha_j$ の項は打ち消し合う。$h_j/h_{j-1}$ の項は一つ前の $P_jP_{j-1}$ の組と逆符号で現れるため、$j=n$ の上端だけが残る。$x\to y$ の極限では
$$
K_n(y,y)=\frac{P_{n+1}'(y)P_n(y)-P_n'(y)P_{n+1}(y)}{h_n}
$$
となる。射影という主張は、最初の和を $P_j(y)$ と掛けて積分すれば確かめられる。特に $n$ 点 Gauss 求積の相異なる二節点 $x_i,x_j$ では $K_{n-1}(x_i,x_j)=0$ である。このため $K_{n-1}(x,x_j)/K_{n-1}(x_j,x_j)$ は Lagrange 基底 $\ell_j$ に等しい。積分すると定数を再現する性質から求積重み $w_j=1/K_{n-1}(x_j,x_j)$ を得る。

## 母関数と Rodrigues の公式

### 微分で直交性を作る

重み $w$ と多項式 $\sigma$ に対して $w^{-1}\dvn{(\sigma^n w)}x n$ が次数 $n$ の多項式になる場合がある。これを定数倍して用いる表示を :term[ロドリゲスの公式]{en="Rodrigues formula" yomi="ろどりげすのこうしき"} という。境界項が零になる条件のもとで $n$ 回部分積分すれば次数 $n$ 未満との直交性が従う。母関数 $\sum c_nP_n(x)t^n$ は、微分した式の係数を比較して漸化式をまとめて得るのに使える。例えば $e^{2xt-t^2}=e^{x^2}e^{-(x-t)^2}$ を $t=0$ で展開すると、次節の Hermite 多項式が $c_n=1/n!$ で現れる。母関数に出る $n!$ は規格化によって変わるため、係数を比較する前に確認する。

### 指数母関数の係数抽出

母関数を使うと、多項式の微分を引数 $x$ 側と補助変数側に分けられる。例えば $e^{2xt-t^2}$ を $x$ で微分すると $2t$ 倍される。このため係数に $n!$ を含める規約では $H_n^\prime=2nH_{n-1}$ が得られる。一方 $t$ で微分すると $(2x-2t)$ 倍されて三項漸化式になる。従って母関数から Rodrigues の公式と同じ列が得られることを、初期条件と漸化式だけでも確かめられる。

:::example[Hermite 母関数から取り出す多項式（1）]
係数 $[t^m]$ を $m!$ 倍して Hermite 多項式にする。

$$
m=1
$$

補助変数の係数を求める。

$$
\begin{aligned}
[t^m]e^{2xt-t^2}\big|_{m=1} &= 2 x
\end{aligned}
$$

次数 $m$ の多項式を得る。最高次係数は $2^m$ である。

$$
\begin{aligned}
H_m(x)\big|_{m=1} &= 2 x
\end{aligned}
$$

微分の因子 $2t$ は添字を一つ下げ、階乗比 $m$ を掛ける。

$$
\begin{aligned}
2mH_{m-1}(x)\big|_{m=1} &= 2
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
[t^m]e^{2xt-t^2}\big|_{x=0} &= 0 \\
\frac{\partial^1 ([t^m]e^{2xt-t^2})}{\partial x^1} &= 2 \\
\left.\frac{\partial^1 ([t^m]e^{2xt-t^2})}{\partial x^1}\right|_{x=0} &= 2
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
H_m(x)\big|_{x=0} &= 0 \\
\frac{\partial^1 (H_m(x))}{\partial x^1} &= 2 \\
\left.\frac{\partial^1 (H_m(x))}{\partial x^1}\right|_{x=0} &= 2
\end{aligned}
$$
:::

:::example[Hermite 母関数から取り出す多項式（2）]
係数 $[t^m]$ を $m!$ 倍して Hermite 多項式にする。次数を一つ上げる。

$$
m=2
$$

補助変数の係数を求める。

$$
\begin{aligned}
[t^m]e^{2xt-t^2}\big|_{m=2} &= 2 x^{2} - 1
\end{aligned}
$$

次数 $m$ の多項式を得る。最高次係数は $2^m$ である。

$$
\begin{aligned}
H_m(x)\big|_{m=2} &= 4 x^{2} - 2
\end{aligned}
$$

微分の因子 $2t$ は添字を一つ下げ、階乗比 $m$ を掛ける。

$$
\begin{aligned}
2mH_{m-1}(x)\big|_{m=2} &= 8 x
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
[t^m]e^{2xt-t^2}\big|_{x=0} &= -1 \\
\frac{\partial^1 ([t^m]e^{2xt-t^2})}{\partial x^1} &= 4 x \\
\left.\frac{\partial^1 ([t^m]e^{2xt-t^2})}{\partial x^1}\right|_{x=0} &= 0
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
H_m(x)\big|_{x=0} &= -2 \\
\frac{\partial^1 (H_m(x))}{\partial x^1} &= 8 x \\
\left.\frac{\partial^1 (H_m(x))}{\partial x^1}\right|_{x=0} &= 0
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
2mH_{m-1}(x)\big|_{x=0} &= 0 \\
\frac{\partial^1 (2mH_{m-1}(x))}{\partial x^1} &= 8 \\
\left.\frac{\partial^1 (2mH_{m-1}(x))}{\partial x^1}\right|_{x=0} &= 8
\end{aligned}
$$
:::

:::example[Hermite 母関数から取り出す多項式（3）]
係数 $[t^m]$ を $m!$ 倍して Hermite 多項式にする。さらに次数を上げて、低次の結果と照合する。

$$
m=3
$$

補助変数の係数を求める。

$$
\begin{aligned}
[t^m]e^{2xt-t^2}\big|_{m=3} &= \frac{4 x^{3}}{3} - 2 x
\end{aligned}
$$

次数 $m$ の多項式を得る。最高次係数は $2^m$ である。

$$
\begin{aligned}
H_m(x)\big|_{m=3} &= 8 x^{3} - 12 x
\end{aligned}
$$

微分の因子 $2t$ は添字を一つ下げ、階乗比 $m$ を掛ける。

$$
\begin{aligned}
2mH_{m-1}(x)\big|_{m=3} &= 24 x^{2} - 12
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
[t^m]e^{2xt-t^2}\big|_{x=0} &= 0 \\
\frac{\partial^1 ([t^m]e^{2xt-t^2})}{\partial x^1} &= 4 x^{2} - 2 \\
\left.\frac{\partial^1 ([t^m]e^{2xt-t^2})}{\partial x^1}\right|_{x=0} &= -2
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
H_m(x)\big|_{x=0} &= 0 \\
\frac{\partial^1 (H_m(x))}{\partial x^1} &= 24 x^{2} - 12 \\
\left.\frac{\partial^1 (H_m(x))}{\partial x^1}\right|_{x=0} &= -12
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
2mH_{m-1}(x)\big|_{x=0} &= -12 \\
\frac{\partial^1 (2mH_{m-1}(x))}{\partial x^1} &= 48 x \\
\left.\frac{\partial^1 (2mH_{m-1}(x))}{\partial x^1}\right|_{x=0} &= 0
\end{aligned}
$$
:::

:::example[Hermite 母関数から取り出す多項式（4）]
係数 $[t^m]$ を $m!$ 倍して Hermite 多項式にする。四次の結果まで具体化する。

$$
m=4
$$

補助変数の係数を求める。

$$
\begin{aligned}
[t^m]e^{2xt-t^2}\big|_{m=4} &= \frac{2 x^{4}}{3} - 2 x^{2} + \frac{1}{2}
\end{aligned}
$$

次数 $m$ の多項式を得る。最高次係数は $2^m$ である。

$$
\begin{aligned}
H_m(x)\big|_{m=4} &= 16 x^{4} - 48 x^{2} + 12
\end{aligned}
$$

微分の因子 $2t$ は添字を一つ下げ、階乗比 $m$ を掛ける。

$$
\begin{aligned}
2mH_{m-1}(x)\big|_{m=4} &= 64 x^{3} - 96 x
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
[t^m]e^{2xt-t^2}\big|_{x=0} &= \frac{1}{2} \\
\frac{\partial^1 ([t^m]e^{2xt-t^2})}{\partial x^1} &= \frac{8 x^{3}}{3} - 4 x \\
\left.\frac{\partial^1 ([t^m]e^{2xt-t^2})}{\partial x^1}\right|_{x=0} &= 0
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
H_m(x)\big|_{x=0} &= 12 \\
\frac{\partial^1 (H_m(x))}{\partial x^1} &= 64 x^{3} - 96 x \\
\left.\frac{\partial^1 (H_m(x))}{\partial x^1}\right|_{x=0} &= 0
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
2mH_{m-1}(x)\big|_{x=0} &= 0 \\
\frac{\partial^1 (2mH_{m-1}(x))}{\partial x^1} &= 192 x^{2} - 96 \\
\left.\frac{\partial^1 (2mH_{m-1}(x))}{\partial x^1}\right|_{x=0} &= -96
\end{aligned}
$$
:::

## Hermite 多項式

### Gaussian 重みと Hermite 多項式

物理で使う :term[エルミート多項式]{en="Hermite polynomial" yomi="えるみーとたこうしき"} は $H_n=(-1)^ne^{x^2}\dvn{e^{-x^2}}x n$ とする。微分から $H_{n+1}=2xH_n-H_n'$、母関数の $x$ 微分から $H_n'=2nH_{n-1}$ が得られる。従って $H_0=1,H_1=2x,H_2=4x^2-2$ である。二式を合わせると $H_n''-2xH_n'+2nH_n=0$ になる。全ての境界項は Gaussian と多項式の積なので消える。$H_n$ と掛けて Rodrigues 表示を $n$ 回部分積分すると、$H_n^{(n)}=2^nn!$ より $\int e^{-x^2}H_n^2=\sqrt\pi\,2^nn!$ である。

:::example[無次元の調和振動子]
$u=e^{-x^2/2}H_n$ とおくと $-u''+x^2u=e^{-x^2/2}[-H_n''+2xH_n'+H_n]=(2n+1)u$ である。物理座標 $X$ から $x=\sqrt{m\omega/\hbar}X$ とすれば、次元をもつエネルギーは $\hbar\omega(n+1/2)$ になる。$n=1$ の無次元波動関数は $\sqrt2\,\pi^{-1/4}x e^{-x^2/2}$ であり、奇関数で原点に節を一つ持つ。
:::

### Gaussian を掛けた形の微分方程式

無次元の振動子 $-\psi^{\prime\prime}+x^2\psi=\epsilon\psi$ に $\psi=e^{-x^2/2}u$ を代入する。二階微分に生じる $x^2 u$ がポテンシャルの項と相殺され、$u^{\prime\prime}-2xu^\prime+(\epsilon-1)u=0$ になる。多項式の最高次 $m$ から $\epsilon=2m+1$ が決まる。多項式解を Gaussian と掛ければ二乗可積分である。固有値の完全性は、等長変換 $\psi\leftrightarrow u$ と、Gaussian 重み付き空間で多項式が稠密であることから従う。稠密性は、すべてのモーメントが零なら $u(x)e^{-x^2}$ の Fourier 変換が全関数で Taylor 係数がすべて零になることを使って示せる。

:::example[Hermite モードの零点と規格化（1）]
Gaussian を掛けた固有関数の原点での値を調べる。

$$
m=1
$$

規格化する前の固有関数を作る。

$$
\begin{aligned}
\psi_m\big|_{m=1} &= 2 x e^{- \frac{x^{2}}{2}}
\end{aligned}
$$

元のエネルギー演算子を作用させる。係数は無次元化した固有値になる。

$$
\begin{aligned}
(-\partial_x^2+x^2)\psi_m\big|_{m=1} &= 6 x e^{- \frac{x^{2}}{2}}
\end{aligned}
$$

ノルムの二乗は同じ多項式を Gaussian 重みで積分した値になる。

$$
\begin{aligned}
\int_\R |\psi_m|^2\dd x\big|_{m=1} &= 2 \sqrt{\pi}
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
\psi_m\big|_{x=0} &= 0 \\
\frac{\partial^1 (\psi_m)}{\partial x^1} &= 2 \left(1 - x^{2}\right) e^{- \frac{x^{2}}{2}} \\
\left.\frac{\partial^1 (\psi_m)}{\partial x^1}\right|_{x=0} &= 2
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
(-\partial_x^2+x^2)\psi_m\big|_{x=0} &= 0 \\
\frac{\partial^1 ((-\partial_x^2+x^2)\psi_m)}{\partial x^1} &= 6 \left(1 - x^{2}\right) e^{- \frac{x^{2}}{2}} \\
\left.\frac{\partial^1 ((-\partial_x^2+x^2)\psi_m)}{\partial x^1}\right|_{x=0} &= 6
\end{aligned}
$$
:::

:::example[Hermite モードの零点と規格化（2）]
Gaussian を掛けた固有関数の原点での値を調べる。次数を一つ上げる。

$$
m=2
$$

規格化する前の固有関数を作る。

$$
\begin{aligned}
\psi_m\big|_{m=2} &= 2 \left(2 x^{2} - 1\right) e^{- \frac{x^{2}}{2}}
\end{aligned}
$$

元のエネルギー演算子を作用させる。係数は無次元化した固有値になる。

$$
\begin{aligned}
(-\partial_x^2+x^2)\psi_m\big|_{m=2} &= 10 \left(2 x^{2} - 1\right) e^{- \frac{x^{2}}{2}}
\end{aligned}
$$

ノルムの二乗は同じ多項式を Gaussian 重みで積分した値になる。

$$
\begin{aligned}
\int_\R |\psi_m|^2\dd x\big|_{m=2} &= 8 \sqrt{\pi}
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
\psi_m\big|_{x=0} &= -2 \\
\frac{\partial^1 (\psi_m)}{\partial x^1} &= \left(- 4 x^{3} + 10 x\right) e^{- \frac{x^{2}}{2}} \\
\left.\frac{\partial^1 (\psi_m)}{\partial x^1}\right|_{x=0} &= 0
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
(-\partial_x^2+x^2)\psi_m\big|_{x=0} &= -10 \\
\frac{\partial^1 ((-\partial_x^2+x^2)\psi_m)}{\partial x^1} &= \left(- 20 x^{3} + 50 x\right) e^{- \frac{x^{2}}{2}} \\
\left.\frac{\partial^1 ((-\partial_x^2+x^2)\psi_m)}{\partial x^1}\right|_{x=0} &= 0
\end{aligned}
$$
:::

:::example[Hermite モードの零点と規格化（3）]
Gaussian を掛けた固有関数の原点での値を調べる。さらに次数を上げて、低次の結果と照合する。

$$
m=3
$$

規格化する前の固有関数を作る。

$$
\begin{aligned}
\psi_m\big|_{m=3} &= \left(8 x^{3} - 12 x\right) e^{- \frac{x^{2}}{2}}
\end{aligned}
$$

元のエネルギー演算子を作用させる。係数は無次元化した固有値になる。

$$
\begin{aligned}
(-\partial_x^2+x^2)\psi_m\big|_{m=3} &= \left(56 x^{3} - 84 x\right) e^{- \frac{x^{2}}{2}}
\end{aligned}
$$

ノルムの二乗は同じ多項式を Gaussian 重みで積分した値になる。

$$
\begin{aligned}
\int_\R |\psi_m|^2\dd x\big|_{m=3} &= 48 \sqrt{\pi}
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
\psi_m\big|_{x=0} &= 0 \\
\frac{\partial^1 (\psi_m)}{\partial x^1} &= 4 \left(- 2 x^{4} + 9 x^{2} - 3\right) e^{- \frac{x^{2}}{2}} \\
\left.\frac{\partial^1 (\psi_m)}{\partial x^1}\right|_{x=0} &= -12
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
(-\partial_x^2+x^2)\psi_m\big|_{x=0} &= 0 \\
\frac{\partial^1 ((-\partial_x^2+x^2)\psi_m)}{\partial x^1} &= 28 \left(- 2 x^{4} + 9 x^{2} - 3\right) e^{- \frac{x^{2}}{2}} \\
\left.\frac{\partial^1 ((-\partial_x^2+x^2)\psi_m)}{\partial x^1}\right|_{x=0} &= -84
\end{aligned}
$$
:::

:::example[Hermite モードの零点と規格化（4）]
Gaussian を掛けた固有関数の原点での値を調べる。四次の結果まで具体化する。

$$
m=4
$$

規格化する前の固有関数を作る。

$$
\begin{aligned}
\psi_m\big|_{m=4} &= 4 \left(4 x^{4} - 12 x^{2} + 3\right) e^{- \frac{x^{2}}{2}}
\end{aligned}
$$

元のエネルギー演算子を作用させる。係数は無次元化した固有値になる。

$$
\begin{aligned}
(-\partial_x^2+x^2)\psi_m\big|_{m=4} &= 36 \left(4 x^{4} - 12 x^{2} + 3\right) e^{- \frac{x^{2}}{2}}
\end{aligned}
$$

ノルムの二乗は同じ多項式を Gaussian 重みで積分した値になる。

$$
\begin{aligned}
\int_\R |\psi_m|^2\dd x\big|_{m=4} &= 384 \sqrt{\pi}
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
\psi_m\big|_{x=0} &= 12 \\
\frac{\partial^1 (\psi_m)}{\partial x^1} &= \left(- 16 x^{5} + 112 x^{3} - 108 x\right) e^{- \frac{x^{2}}{2}} \\
\left.\frac{\partial^1 (\psi_m)}{\partial x^1}\right|_{x=0} &= 0
\end{aligned}
$$

$x=0$ で値と導関数を評価する。

$$
\begin{aligned}
(-\partial_x^2+x^2)\psi_m\big|_{x=0} &= 108 \\
\frac{\partial^1 ((-\partial_x^2+x^2)\psi_m)}{\partial x^1} &= \left(- 144 x^{5} + 1008 x^{3} - 972 x\right) e^{- \frac{x^{2}}{2}} \\
\left.\frac{\partial^1 ((-\partial_x^2+x^2)\psi_m)}{\partial x^1}\right|_{x=0} &= 0
\end{aligned}
$$
:::

:::example[Gaussian 重みでの次数の整理]
所定の重みで直交射影し、与えた関数のうち低次のモードがどこまで回収されるか調べる。

投影する関数 $f=x^{6}$ を取り、同じ区間と重みで各成分の積分を計算する。正規化済みでない基底も使用できる。

$$
c_j=\frac{\int_{-\infty}^{\infty} f\phi_j\,w\,\dd x}{\int_{-\infty}^{\infty}\phi_j^2w\,\dd x},\qquad w=e^{- x^{2}}
$$

第 1 の方向に対し、分子と規格化のための分母を別々に計算する。

$$
\begin{aligned}
\phi_{0} &= 1 \\
\int f\phi_{0}w\,\dd x &= \frac{15 \sqrt{\pi}}{8} \\
\int \phi_{0}^2w\,\dd x &= \sqrt{\pi} \\
c_{0} &= \frac{15}{8}
\end{aligned}
$$

第 2 の方向に対し、分子と規格化のための分母を別々に計算する。

$$
\begin{aligned}
\phi_{1} &= 2 x \\
\int f\phi_{1}w\,\dd x &= 0 \\
\int \phi_{1}^2w\,\dd x &= 2 \sqrt{\pi} \\
c_{1} &= 0
\end{aligned}
$$

第 3 の方向に対し、分子と規格化のための分母を別々に計算する。

$$
\begin{aligned}
\phi_{2} &= 4 x^{2} - 2 \\
\int f\phi_{2}w\,\dd x &= \frac{45 \sqrt{\pi}}{2} \\
\int \phi_{2}^2w\,\dd x &= 8 \sqrt{\pi} \\
c_{2} &= \frac{45}{16}
\end{aligned}
$$

第 4 の方向に対し、分子と規格化のための分母を別々に計算する。

$$
\begin{aligned}
\phi_{3} &= 8 x^{3} - 12 x \\
\int f\phi_{3}w\,\dd x &= 0 \\
\int \phi_{3}^2w\,\dd x &= 48 \sqrt{\pi} \\
c_{3} &= 0
\end{aligned}
$$

第 5 の方向に対し、分子と規格化のための分母を別々に計算する。

$$
\begin{aligned}
\phi_{4} &= 16 x^{4} - 48 x^{2} + 12 \\
\int f\phi_{4}w\,\dd x &= 180 \sqrt{\pi} \\
\int \phi_{4}^2w\,\dd x &= 384 \sqrt{\pi} \\
c_{4} &= \frac{15}{32}
\end{aligned}
$$

第 6 の方向に対し、分子と規格化のための分母を別々に計算する。

$$
\begin{aligned}
\phi_{5} &= 32 x^{5} - 160 x^{3} + 120 x \\
\int f\phi_{5}w\,\dd x &= 0 \\
\int \phi_{5}^2w\,\dd x &= 3840 \sqrt{\pi} \\
c_{5} &= 0
\end{aligned}
$$

第 7 の方向に対し、分子と規格化のための分母を別々に計算する。

$$
\begin{aligned}
\phi_{6} &= 64 x^{6} - 480 x^{4} + 720 x^{2} - 120 \\
\int f\phi_{6}w\,\dd x &= 720 \sqrt{\pi} \\
\int \phi_{6}^2w\,\dd x &= 46080 \sqrt{\pi} \\
c_{6} &= \frac{1}{64}
\end{aligned}
$$

成分を足すと、指定された部分空間への射影が得られる。残差の内積は、使った各基底方向で零になる。

$$
\Pi f=x^{6},\qquad f-\Pi f=0
$$
:::

## Laguerre 多項式と陪 Laguerre 多項式

### 指数減衰重みと Laguerre 多項式

$\alpha>-1$、$x>0$ で $L_n^{(\alpha)}=x^{-\alpha}e^x\dvn{(e^{-x}x^{n+\alpha})}x n/n!$ とする。$L_n^{(\alpha)}$ を :term[一般化ラゲール多項式]{en="generalized Laguerre polynomial" yomi="いっぱんからげーるたこうしき"} といい、$L_n=L_n^{(0)}$ である。積の微分を展開すると
$$
L_n^{(\alpha)}(x)=\sum_{j=0}^n \frac{(-1)^j}{j!}\binom{n+\alpha}{n-j}x^j.
$$
二項係数の係数比から $xy''+(\alpha+1-x)y'+ny=0$ に代入して零になる。境界では $n-1$ 階以下の微分でも $x^{\alpha+1}$ 以上の冪が残るため部分積分の境界項が消える。先頭係数は $(-1)^n/n!$ なので $n$ 回部分積分により $\int_0^\infty x^\alpha e^{-x}(L_n^{(\alpha)})^2=\Gamma(n+\alpha+1)/n!$ である。例えば $L_1^{(2)}=3-x$ で、ノルムの二乗は $6$ である。係数を $n$ で足し、二項級数を使うと母関数 $(1-t)^{-\alpha-1}\exp[-xt/(1-t)]$ が得られる。整数 $m\ge0$ では係数を $m$ 回ずらすことで $(-1)^m \dvn{L_{n+m}}x m=L_n^{(m)}$ になる。

### 漸化式と重み付き空間での完全性

母関数を $G(x,t)$ とすると
$$
(1-t)^2\partial_tG=[(\alpha+1)(1-t)-x]G
$$
である。係数比較により
$$
 (n+1)L_{n+1}^{(\alpha)}=(2n+\alpha+1-x)L_n^{(\alpha)}
 -(n+\alpha)L_{n-1}^{(\alpha)}
$$
を得る。$n=0$ では最後の項を省く。一方、$x$ 微分は $\partial_xG=-tG/(1-t)$ である。別の添字 $\alpha+1$ の母関数は $G/(1-t)$ なので、$(L_n^{(\alpha)})'=-L_{n-1}^{(\alpha+1)}$ となる。

$w(x)=x^\alpha e^{-x}$、$\alpha>-1$ とし、$f\in L^2(w\dd x)$ が全ての多項式と直交するとする。$A(z)=\int_0^\infty e^{zx}f(x)w(x)\dd x$ は $\Re z<1/2$ で正則である。実際、Cauchy–Schwarz により、$\Re z\le1/2-\delta$ で任意の $x^j e^{zx}$ が $L^2(w\dd x)$ で一様に有界になる。$A^{(j)}(0)=0$ が全次数で成り立つため恒等定理から $A=0$ である。虚軸上の $A$ は $L^1$ 関数 $fw$ を負の半直線で零に延ばしたものの Fourier 変換に一致する。その一意性から $fw=0$、従って $f=0$ である。この議論により $L_n^{(\alpha)}$ を規格化した列は重み付き $L^2$ の正規直交基底になる。Hermite の場合も同じ論法で、Gaussian のため $A(z)$ は半平面に限らず整関数になる。

### 一般化 Laguerre の直交性と端での挙動

Rodrigues 表示の $n$ 階微分を多項式側に移すと、次数が $n$ 未満なら零となる。原点の境界項では $x^{\alpha+n}$ の微分を高々 $n-1$ 回行うため、少なくとも $x^{\alpha+1}$ が残る。よって $\alpha>-1$ なら境界項は消える。無限遠の境界項は指数関数が消す。相手を $L_n^{(\alpha)}$ にすると、$n$ 階微分は最高次係数から $(-1)^n$ になる。従ってノルムは $\Gamma(n+\alpha+1)/n!$ である。

:::example[重みを変えた Laguerre 多項式（1）]
次数を $2$ に固定し、重みの指数を変える。

$$
m=1
$$

原点の各係数を規格化から決める。

$$
\begin{aligned}
L_2^{(m)}(x)\big|_{m=1} &= \frac{x^{2}}{2} - 3 x + 3
\end{aligned}
$$

重みを掛けた二乗を積分する。展開すれば各項は Gamma 積分である。

$$
\begin{aligned}
\int_0^\infty (L_2^{(m)})^2 x^m e^{-x}\dd x\big|_{m=1} &= 3
\end{aligned}
$$

一次の多項式と内積を取ると、二度の部分積分で消える。

$$
\begin{aligned}
\langle L_2^{(m)},L_1^{(m)}\rangle &= 0
\end{aligned}
$$

$x=1$ で値と導関数を評価する。

$$
\begin{aligned}
L_2^{(m)}(x)\big|_{x=1} &= \frac{1}{2} \\
\frac{\partial^1 (L_2^{(m)}(x))}{\partial x^1} &= x - 3 \\
\left.\frac{\partial^1 (L_2^{(m)}(x))}{\partial x^1}\right|_{x=1} &= -2 \\
\frac{\partial^2 (L_2^{(m)}(x))}{\partial x^2} &= 1 \\
\left.\frac{\partial^2 (L_2^{(m)}(x))}{\partial x^2}\right|_{x=1} &= 1
\end{aligned}
$$
:::

:::example[重みを変えた Laguerre 多項式（2）]
次数を $2$ に固定し、重みの指数を変える。次数を一つ上げる。

$$
m=2
$$

原点の各係数を規格化から決める。

$$
\begin{aligned}
L_2^{(m)}(x)\big|_{m=2} &= \frac{x^{2}}{2} - 4 x + 6
\end{aligned}
$$

重みを掛けた二乗を積分する。展開すれば各項は Gamma 積分である。

$$
\begin{aligned}
\int_0^\infty (L_2^{(m)})^2 x^m e^{-x}\dd x\big|_{m=2} &= 12
\end{aligned}
$$

一次の多項式と内積を取ると、二度の部分積分で消える。

$$
\begin{aligned}
\langle L_2^{(m)},L_1^{(m)}\rangle &= 0
\end{aligned}
$$

$x=1$ で値と導関数を評価する。

$$
\begin{aligned}
L_2^{(m)}(x)\big|_{x=1} &= \frac{5}{2} \\
\frac{\partial^1 (L_2^{(m)}(x))}{\partial x^1} &= x - 4 \\
\left.\frac{\partial^1 (L_2^{(m)}(x))}{\partial x^1}\right|_{x=1} &= -3 \\
\frac{\partial^2 (L_2^{(m)}(x))}{\partial x^2} &= 1 \\
\left.\frac{\partial^2 (L_2^{(m)}(x))}{\partial x^2}\right|_{x=1} &= 1
\end{aligned}
$$
:::

:::example[重みを変えた Laguerre 多項式（3）]
次数を $2$ に固定し、重みの指数を変える。さらに次数を上げて、低次の結果と照合する。

$$
m=3
$$

原点の各係数を規格化から決める。

$$
\begin{aligned}
L_2^{(m)}(x)\big|_{m=3} &= \frac{x^{2}}{2} - 5 x + 10
\end{aligned}
$$

重みを掛けた二乗を積分する。展開すれば各項は Gamma 積分である。

$$
\begin{aligned}
\int_0^\infty (L_2^{(m)})^2 x^m e^{-x}\dd x\big|_{m=3} &= 60
\end{aligned}
$$

一次の多項式と内積を取ると、二度の部分積分で消える。

$$
\begin{aligned}
\langle L_2^{(m)},L_1^{(m)}\rangle &= 0
\end{aligned}
$$

$x=1$ で値と導関数を評価する。

$$
\begin{aligned}
L_2^{(m)}(x)\big|_{x=1} &= \frac{11}{2} \\
\frac{\partial^1 (L_2^{(m)}(x))}{\partial x^1} &= x - 5 \\
\left.\frac{\partial^1 (L_2^{(m)}(x))}{\partial x^1}\right|_{x=1} &= -4 \\
\frac{\partial^2 (L_2^{(m)}(x))}{\partial x^2} &= 1 \\
\left.\frac{\partial^2 (L_2^{(m)}(x))}{\partial x^2}\right|_{x=1} &= 1
\end{aligned}
$$
:::

:::example[重みを変えた Laguerre 多項式（4）]
次数を $2$ に固定し、重みの指数を変える。四次の結果まで具体化する。

$$
m=4
$$

原点の各係数を規格化から決める。

$$
\begin{aligned}
L_2^{(m)}(x)\big|_{m=4} &= \frac{x^{2}}{2} - 6 x + 15
\end{aligned}
$$

重みを掛けた二乗を積分する。展開すれば各項は Gamma 積分である。

$$
\begin{aligned}
\int_0^\infty (L_2^{(m)})^2 x^m e^{-x}\dd x\big|_{m=4} &= 360
\end{aligned}
$$

一次の多項式と内積を取ると、二度の部分積分で消える。

$$
\begin{aligned}
\langle L_2^{(m)},L_1^{(m)}\rangle &= 0
\end{aligned}
$$

$x=1$ で値と導関数を評価する。

$$
\begin{aligned}
L_2^{(m)}(x)\big|_{x=1} &= \frac{19}{2} \\
\frac{\partial^1 (L_2^{(m)}(x))}{\partial x^1} &= x - 6 \\
\left.\frac{\partial^1 (L_2^{(m)}(x))}{\partial x^1}\right|_{x=1} &= -5 \\
\frac{\partial^2 (L_2^{(m)}(x))}{\partial x^2} &= 1 \\
\left.\frac{\partial^2 (L_2^{(m)}(x))}{\partial x^2}\right|_{x=1} &= 1
\end{aligned}
$$
:::

:::example[指数重みでの多項式の分解]
所定の重みで直交射影し、与えた関数のうち低次のモードがどこまで回収されるか調べる。

投影する関数 $f=x^{5}$ を取り、同じ区間と重みで各成分の積分を計算する。正規化済みでない基底も使用できる。

$$
c_j=\frac{\int_{0}^{\infty} f\phi_j\,w\,\dd x}{\int_{0}^{\infty}\phi_j^2w\,\dd x},\qquad w=e^{- x}
$$

第 1 の方向に対し、分子と規格化のための分母を別々に計算する。

$$
\begin{aligned}
\phi_{0} &= 1 \\
\int f\phi_{0}w\,\dd x &= 120 \\
\int \phi_{0}^2w\,\dd x &= 1 \\
c_{0} &= 120
\end{aligned}
$$

第 2 の方向に対し、分子と規格化のための分母を別々に計算する。

$$
\begin{aligned}
\phi_{1} &= 1 - x \\
\int f\phi_{1}w\,\dd x &= -600 \\
\int \phi_{1}^2w\,\dd x &= 1 \\
c_{1} &= -600
\end{aligned}
$$

第 3 の方向に対し、分子と規格化のための分母を別々に計算する。

$$
\begin{aligned}
\phi_{2} &= \frac{x^{2}}{2} - 2 x + 1 \\
\int f\phi_{2}w\,\dd x &= 1200 \\
\int \phi_{2}^2w\,\dd x &= 1 \\
c_{2} &= 1200
\end{aligned}
$$

第 4 の方向に対し、分子と規格化のための分母を別々に計算する。

$$
\begin{aligned}
\phi_{3} &= - \frac{x^{3}}{6} + \frac{3 x^{2}}{2} - 3 x + 1 \\
\int f\phi_{3}w\,\dd x &= -1200 \\
\int \phi_{3}^2w\,\dd x &= 1 \\
c_{3} &= -1200
\end{aligned}
$$

第 5 の方向に対し、分子と規格化のための分母を別々に計算する。

$$
\begin{aligned}
\phi_{4} &= \frac{x^{4}}{24} - \frac{2 x^{3}}{3} + 3 x^{2} - 4 x + 1 \\
\int f\phi_{4}w\,\dd x &= 600 \\
\int \phi_{4}^2w\,\dd x &= 1 \\
c_{4} &= 600
\end{aligned}
$$

第 6 の方向に対し、分子と規格化のための分母を別々に計算する。

$$
\begin{aligned}
\phi_{5} &= - \frac{x^{5}}{120} + \frac{5 x^{4}}{24} - \frac{5 x^{3}}{3} + 5 x^{2} - 5 x + 1 \\
\int f\phi_{5}w\,\dd x &= -120 \\
\int \phi_{5}^2w\,\dd x &= 1 \\
c_{5} &= -120
\end{aligned}
$$

成分を足すと、指定された部分空間への射影が得られる。残差の内積は、使った各基底方向で零になる。

$$
\Pi f=x^{5},\qquad f-\Pi f=0
$$
:::

## Chebyshev 多項式と Gauss 求積

### 余弦と節点を使った積分

$x=\cos\theta$ と置き $T_n(x)=\cos(n\theta)$ とする。余弦の加法定理より $T_{n+1}=2xT_n-T_{n-1}$、$T_0=1,T_1=x$ なので多項式である。$w=(1-x^2)^{-1/2}$ での内積は $[0,\pi]$ の余弦の直交性に帰着し、ノルムの二乗は $n=0$ で $\pi$、それ以外で $\pi/2$ である。第 2 種 $U_n(\cos\theta)=\sin((n+1)\theta)/\sin\theta$ は同じ漸化式で $U_0=1,U_1=2x$ から始まり、重み $\sqrt{1-x^2}$ のもとでノルムの二乗は $\pi/2$ である。

:::proposition[Gauss 求積]{#prop-mp-gauss-quadrature}
$n$ 次直交多項式の根を $x_j$、その根に対する Lagrange 基底を $\ell_j$ とする。重み $w_j=\int w\ell_j$ とする求積は次数 $2n-1$ 以下で正確で、$w_j>0$ である。
:::

:::proof
多項式 $f=P_nq+r$ と除算すると $\deg q,\deg r<n$ で、直交性によって $\int wf=\int wr$ になる。根上では $f=r$ であり、Lagrange 補間で $\int wr=\sum_jw_jf(x_j)$ である。$\ell_j^2$ にも正確性が適用できて $w_j=\int w\ell_j^2>0$ を得る。
:::

Chebyshev の根は $x_j=\cos((2j-1)\pi/(2n))$ である。この場合 $w_j=\pi/n$ である。$1\le m<2n$ で $\sum_j\cos(m(2j-1)\pi/(2n))=0$ となることは、$m$ が奇数なら対称性、偶数なら等比和で分かる。従って等重みは $T_0,\ldots,T_{2n-1}$ を正確に積分する。例えば $n=2$ では節点 $\pm1/\sqrt2$ を使い、$x^2$ の積分は $\pi/2$ となる。

:::example[重みが一定の二点公式]
同じ区間で前の例の $P_2$ に対応する二点公式を作る。反射対称性により二つの重みは同じで、定数 $1$ の積分は $2$ だから各 $1$ である。従って $f=x^2$ では節点値の和が $2/3$ で厳密値と一致する。$x^4$ では節点値の和は $2/9$、厳密値は $2/5$ である。$x^4$ は二点公式の保証次数 $3$ を超える。この重みは Chebyshev の $w$ と異なり、節点も $\pm1/\sqrt2$ から移る。
:::

### 求積の誤差の符号

$f\in C^{2n}([a,b])$ とする。Gauss 節点 $x_1,\ldots,x_n$ で $f$ と $f'$ に一致する次数 $2n-1$ 以下の Hermite 補間多項式 $H$ は一意である。同じデータを持つ二つの差には $n$ 個の二重零点があり、次数の上限から零となるためである。係数とデータの次元がともに $2n$ なので存在も従う。任意の節点以外の $x$ で $f-H-A P_n^2$ が $x$ でも零となるよう $A$ を選ぶ。重複度込みで $2n+1$ 個の零があるため、Rolle の定理を $2n$ 回使うとある $\xi_x$ で $A=f^{(2n)}(\xi_x)/(2n)!$ となる。ここで $P_n$ はモニックである。よって $f^{(2n)}\ge0$ なら $f-H\ge0$ であり、Gauss 値は積分値以下となる。誤差の絶対値は $\sup|f^{(2n)}|\,h_n/(2n)!$ 以下である。例えば単位重みの $[-1,1]$ で $n=2,f(x)=x^4$ とすると、モニック多項式は $P_2=x^2-1/3$ である。二点求積は $2/9$、積分値は $2/5$ で、差は $8/45$ になる。$f-H=P_2^2$ なので、この差は $h_2=\int_{-1}^1(x^2-1/3)^2\dd x$ に等しい。

### Chebyshev 展開と微分方程式

$x=\cos\theta$ を使うと $\partial_\theta^2$ は $(1-x^2)\partial_x^2-x\partial_x$ になる。従って $T_n$ は $(1-x^2)T_n''-xT_n'+n^2T_n=0$ を満たす。また $T_n'(x)=n\sin(n\theta)/\sin\theta=nU_{n-1}(x)$ である。$1-x^2$ で割って式を正規形にする場合には、端点ではなく $-1<x<1$ を考える。

重み $(1-x^2)^{-1/2}$ における $f$ の直交展開の係数は $c_0=\pi^{-1}\int_0^\pi f(\cos\theta)\dd\theta$、$n\ge1$ で $c_n=2\pi^{-1}\int_0^\pi f(\cos\theta)\cos(n\theta)\dd\theta$ である。偶拡張に対する Fourier 級数の完全性から、和 $c_0+\sum_{n\ge1}c_nT_n$ は重み付き $L^2$ で $f$ に収束する。例えば $x^4=(3T_0+4T_2+T_4)/8$ は余弦の倍角公式から得られ、重み付き積分値は定数成分だけから $3\pi/8$ になる。

### 極座標の角を使った求積

Chebyshev 第 1 種の重みでは $x=\cos\theta$ により $w(x)\dd x$ が $\dd\theta$ となる。このとき $T_j(x)$ は $\cos(j\theta)$ になる。節点を $\theta=(2i-1)\pi/(2m)$ に取ると、$1\le j<2m$ の余弦の離散和は零になる。$j$ が奇数なら角を $\pi-\theta$ と組にでき、偶数なら等比級数で零になる。$j=0$ の和は $m$ なので重みを $\pi/m$ にすると定積分が合う。$T_0,\ldots,T_{2m-1}$ はその次数までの全多項式を張る。

:::example[Chebyshev 求積の次数と最初の誤差（1）]
$T_{2m}$ の連続積分と離散和を比べる。

$$
m=1
$$

節点を決める角を $j=1,\ldots,m$ で添字付けする。

$$
\begin{aligned}
\theta_j\big|_{m=1} &= \pi \left(n - \frac{1}{2}\right)
\end{aligned}
$$

ここで $T_{2m}(\cos\theta_j)=\cos((2j-1)\pi)=-1$ になる。

$$
\begin{aligned}
Q_m(T_{2m}) &= - \pi
\end{aligned}
$$

連続積分では $T_{2m}$ は定数と直交するため零である。

$$
\begin{aligned}
\int_{-1}^1T_{2m}(x)(1-x^2)^{-1/2}\dd x &= 0
\end{aligned}
$$
:::

:::example[Chebyshev 求積の次数と最初の誤差（2）]
$T_{2m}$ の連続積分と離散和を比べる。次数を一つ上げる。

$$
m=2
$$

節点を決める角を $j=1,\ldots,m$ で添字付けする。

$$
\begin{aligned}
\theta_j\big|_{m=2} &= \frac{\pi \left(2 n - 1\right)}{4}
\end{aligned}
$$

ここで $T_{2m}(\cos\theta_j)=\cos((2j-1)\pi)=-1$ になる。

$$
\begin{aligned}
Q_m(T_{2m}) &= - \pi
\end{aligned}
$$

連続積分では $T_{2m}$ は定数と直交するため零である。

$$
\begin{aligned}
\int_{-1}^1T_{2m}(x)(1-x^2)^{-1/2}\dd x &= 0
\end{aligned}
$$
:::

:::example[Chebyshev 求積の次数と最初の誤差（3）]
$T_{2m}$ の連続積分と離散和を比べる。さらに次数を上げて、低次の結果と照合する。

$$
m=3
$$

節点を決める角を $j=1,\ldots,m$ で添字付けする。

$$
\begin{aligned}
\theta_j\big|_{m=3} &= \frac{\pi \left(2 n - 1\right)}{6}
\end{aligned}
$$

ここで $T_{2m}(\cos\theta_j)=\cos((2j-1)\pi)=-1$ になる。

$$
\begin{aligned}
Q_m(T_{2m}) &= - \pi
\end{aligned}
$$

連続積分では $T_{2m}$ は定数と直交するため零である。

$$
\begin{aligned}
\int_{-1}^1T_{2m}(x)(1-x^2)^{-1/2}\dd x &= 0
\end{aligned}
$$
:::

:::example[Chebyshev 求積の次数と最初の誤差（4）]
$T_{2m}$ の連続積分と離散和を比べる。四次の結果まで具体化する。

$$
m=4
$$

節点を決める角を $j=1,\ldots,m$ で添字付けする。

$$
\begin{aligned}
\theta_j\big|_{m=4} &= \frac{\pi \left(2 n - 1\right)}{8}
\end{aligned}
$$

ここで $T_{2m}(\cos\theta_j)=\cos((2j-1)\pi)=-1$ になる。

$$
\begin{aligned}
Q_m(T_{2m}) &= - \pi
\end{aligned}
$$

連続積分では $T_{2m}$ は定数と直交するため零である。

$$
\begin{aligned}
\int_{-1}^1T_{2m}(x)(1-x^2)^{-1/2}\dd x &= 0
\end{aligned}
$$
:::

### 区間の平行移動と節点の重み

Gauss–Legendre 求積を $[a,b]$ へ移す場合は、$x=(a+b)/2+(b-a)t/2$ と変数変換する。単位重みはそのままで、Jacobian として係数に $(b-a)/2$ が掛かる。二点公式の節点は $t=\pm1/\sqrt3$、各重みは $1$ である。この対称配置で奇数次は相殺し、定数および二次との積分を合わせれば三次まで正確になる。次の例では、和の値が正確な次数と、最初に誤差が出る次数を同じ区間で比べる。
:::example[区間 0 から 1 の二点公式]
二節点を区間に移し、一次から四次までの項を積分する。

節点と両方の重みを列挙する。

$$
\begin{aligned}
x_- &= \frac{1}{2} - \frac{\sqrt{3}}{6} \\
x_+ &= \frac{\sqrt{3}}{6} + \frac{1}{2} \\
w &= \frac{1}{2}
\end{aligned}
$$

1 次について、二点で値を計算して重み付き和を取り、原始関数から得る積分と比較する。

$$
\begin{aligned}
Q(x^{1}) &= \frac{1}{2} \\
\int_a^b x^1\dd x &= \frac{1}{2}
\end{aligned}
$$

2 次について、二点で値を計算して重み付き和を取り、原始関数から得る積分と比較する。

$$
\begin{aligned}
Q(x^{2}) &= \frac{1}{3} \\
\int_a^b x^2\dd x &= \frac{1}{3}
\end{aligned}
$$

3 次について、二点で値を計算して重み付き和を取り、原始関数から得る積分と比較する。

$$
\begin{aligned}
Q(x^{3}) &= \frac{1}{4} \\
\int_a^b x^3\dd x &= \frac{1}{4}
\end{aligned}
$$

4 次について、二点で値を計算して重み付き和を取り、原始関数から得る積分と比較する。

$$
\begin{aligned}
Q(x^{4}) &= \frac{7}{36} \\
\int_a^b x^4\dd x &= \frac{1}{5}
\end{aligned}
$$
:::


:::example[区間 0 から 2 の二点公式]
二節点を区間に移し、一次から四次までの項を積分する。

節点と両方の重みを列挙する。

$$
\begin{aligned}
x_- &= 1 - \frac{\sqrt{3}}{3} \\
x_+ &= \frac{\sqrt{3}}{3} + 1 \\
w &= 1
\end{aligned}
$$

1 次について、二点で値を計算して重み付き和を取り、原始関数から得る積分と比較する。

$$
\begin{aligned}
Q(x^{1}) &= 2 \\
\int_a^b x^1\dd x &= 2
\end{aligned}
$$

2 次について、二点で値を計算して重み付き和を取り、原始関数から得る積分と比較する。

$$
\begin{aligned}
Q(x^{2}) &= \frac{8}{3} \\
\int_a^b x^2\dd x &= \frac{8}{3}
\end{aligned}
$$

3 次について、二点で値を計算して重み付き和を取り、原始関数から得る積分と比較する。

$$
\begin{aligned}
Q(x^{3}) &= 4 \\
\int_a^b x^3\dd x &= 4
\end{aligned}
$$

4 次について、二点で値を計算して重み付き和を取り、原始関数から得る積分と比較する。

$$
\begin{aligned}
Q(x^{4}) &= \frac{56}{9} \\
\int_a^b x^4\dd x &= \frac{32}{5}
\end{aligned}
$$
:::


:::example[区間 -1 から 2 の二点公式]
二節点を区間に移し、一次から四次までの項を積分する。

節点と両方の重みを列挙する。

$$
\begin{aligned}
x_- &= \frac{1}{2} - \frac{\sqrt{3}}{2} \\
x_+ &= \frac{1}{2} + \frac{\sqrt{3}}{2} \\
w &= \frac{3}{2}
\end{aligned}
$$

1 次について、二点で値を計算して重み付き和を取り、原始関数から得る積分と比較する。

$$
\begin{aligned}
Q(x^{1}) &= \frac{3}{2} \\
\int_a^b x^1\dd x &= \frac{3}{2}
\end{aligned}
$$

2 次について、二点で値を計算して重み付き和を取り、原始関数から得る積分と比較する。

$$
\begin{aligned}
Q(x^{2}) &= 3 \\
\int_a^b x^2\dd x &= 3
\end{aligned}
$$

3 次について、二点で値を計算して重み付き和を取り、原始関数から得る積分と比較する。

$$
\begin{aligned}
Q(x^{3}) &= \frac{15}{4} \\
\int_a^b x^3\dd x &= \frac{15}{4}
\end{aligned}
$$

4 次について、二点で値を計算して重み付き和を取り、原始関数から得る積分と比較する。

$$
\begin{aligned}
Q(x^{4}) &= \frac{21}{4} \\
\int_a^b x^4\dd x &= \frac{33}{5}
\end{aligned}
$$
:::

## 演習問題

::::exercise{level=basic}
$\Gamma(7/2)$ を求めよ。
:::solution
$\Gamma(1/2)=\sqrt\pi$ に漸化式で $1/2,3/2,5/2$ を掛けると $15\sqrt\pi/8$ になる。
:::
::::

::::exercise{level=basic}
$B(3,2)$ を求めよ。
:::solution
式 :eqref[eq-mp-beta-gamma] より $2!1!/4!=1/12$ である。$\int_0^1 x^2(1-x)\dd x=1/3-1/4$ とも一致する。
:::
::::

::::exercise{level=basic}
$H_3$ を求めよ。
:::solution
$H_3=2xH_2-4H_1=2x(4x^2-2)-8x=8x^3-12x$ である。
:::
::::

::::exercise{level=basic}
$T_3$ の根を求めよ。
:::solution
$3\theta=\pi/2,3\pi/2,5\pi/2$ より $\theta=\pi/6,\pi/2,5\pi/6$ なので根は $\sqrt3/2,0,-\sqrt3/2$ である。
:::
::::

::::exercise{level=advanced}
$\psi(1/2)=-\gamma-2\log 2$ を示せ。
:::solution
倍角公式を対数微分すると $\psi(z)+\psi(z+1/2)=-2\log2+2\psi(2z)$ である。$z=1/2$ で $\psi(1)=-\gamma$ を代入して整理する。
:::
::::

::::exercise{level=advanced}
$L_2^{(1)}$ の式とノルムの二乗を求めよ。
:::solution
係数式から $3-3x+x^2/2$ である。ノルムは $\Gamma(4)/2!=3$ である。展開して積分しても、各単項式の積分が階乗になり同じ値を返す。
:::
::::

::::exercise{level=advanced}
正重みについて $n$ 点 Gauss 求積が次数 $2n$ に対しては一般に正確でないことを示せ。
:::solution
$f=P_n^2$ を選ぶと全節点で零なので求積値は零である。一方 $\int wP_n^2=h_n>0$ である。この反例はどの正重みでも使える。
:::
::::

::::exercise{level=advanced}
$\zeta(2)$ を用いて $\psi'(1)$ を評価せよ。
:::solution
ディガンマ級数の導関数の第 $n$ 項は $1/(n+z)^2$ である。$z=1$ においてその和は $\zeta(2)$ であり、式 :eqref[eq-mp-cotangent-expansion] の比較から $\pi^2/6$ となる。
:::
::::

## 確認テスト

::::quiz{answer="6" tol="1e-9"}
$\Gamma(4)$ を求めよ。

:::explanation
$\Gamma(1)=1$ に $1,2,3$ を掛けて $6$ である。
:::
::::

::::quiz
Gauss 求積の節点の選び方はどれか。

- [x] 対応する重みに対する直交多項式の根
- [ ] 常に等間隔
- [ ] 直交多項式の極大点だけ

:::explanation
節点で $P_n$ を零にすれば、剰余 $r$ の補間に置き換えられる。積分でも $P_nq$ は直交性で消える。
:::
::::

::::quiz{answer="-2" tol="1e-9"}
$H_2(0)$ を求めよ。

:::explanation
$H_2=2x(2x)-2=4x^2-2$ の原点値は $-2$ である。
:::
::::

::::quiz
Chebyshev 第 1 種の重みはどれか。

- [ ] $e^{-x^2}$
- [ ] $1$
- [x] $(1-x^2)^{-1/2}$

:::explanation
$x=\cos\theta$ では $|\dd x|=\sin\theta\,\dd\theta$ なので、$\sin\theta=\sqrt{1-x^2}$ で割ると一定の測度 $\dd\theta$ を得る。
:::
::::

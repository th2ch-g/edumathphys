---
title: "Bessel 関数"
description: "円筒・球対称の固有方程式を解き、Bessel 関数の規格化、漸近形と変換を導く。"
notation:
  - symbol: 'J_\nu,Y_\nu'
    meaning: "第 1 種・第 2 種 Bessel 関数"
  - symbol: 'H_\nu^{(1)},H_\nu^{(2)}'
    meaning: "外向き・内向き波に使う Hankel 関数"
  - symbol: 'I_\nu,K_\nu'
    meaning: "変形 Bessel 関数"
  - symbol: 'j_\ell,y_\ell'
    meaning: "球 Bessel 関数"
  - symbol: 'j_{\nu,n}'
    meaning: "Jν の第 n 正零点"
---

## この章の目標

- Bessel 方程式を Frobenius 法で解く。
- 母関数から漸化式と角度積分を得る。
- 放射する波と指数減衰する解を区別する。
- 有限円板の展開と全平面の Hankel 変換を構成する。

## Bessel 方程式と第 1 種・第 2 種 Bessel 関数

### 円筒の半径方向を分離する

$\Delta u+k^2u=0$ を平面の極座標で考え、角度因子を $e^{im\phi}$ とすると $R''+r^{-1}R'+(k^2-m^2/r^2)R=0$ となる。$x=kr$ は無次元であり、$m$ を一般の次数 $\nu$ に替えた式 $x^2y''+xy'+(x^2-\nu^2)y=0$ を :term[ベッセル方程式]{en="Bessel equation" yomi="べっせるほうていしき"} という。

:::definition[円筒 Bessel 関数]{#def-mp-bessel}
$x>0$、複素次数 $\nu$ について第 1 種を
$$
J_\nu(x)=\sum_{j=0}^\infty\frac{(-1)^j(x/2)^{2j+\nu}}{j!\Gamma(j+\nu+1)}
\label{eq-mp-bessel-series}
$$
とする。第 2 種は非整数の $\nu$ で $Y_\nu=(\cos\pi\nu\,J_\nu-J_{-\nu})/\sin\pi\nu$、整数では次数での極限で定義する。複素引数を用いる場合は $-\pi<\arg x<\pi$ の主枝を取る。
:::

$x^{\nu+2j}$ の係数を $A_j$ とすると、方程式の同じ次数の係数が消える条件は $4j(j+\nu)A_j+A_{j-1}=0$ である。式 :eqref[eq-mp-bessel-series] の係数は逆 Gamma の恒等式 $z/\Gamma(z+1)=1/\Gamma(z)$ により、$j+\nu=0$ の場合もこの関係を満たす。分母が零にならない範囲では係数比は $-1/[4j(j+\nu)]$ で、Frobenius 法の式 :eqref[eq-mp-frobenius-recurrence] に一致する。十分大きい $j$ ではこの比が次数のコンパクト集合で一様に $O(j^{-2})$ なので、$J_\nu$ は次数について整関数であり、引数の微分と次数の極限は交換できる。負整数次数では逆 Gamma の零により初めの $n$ 項が消え、残りの添字をずらすと $J_{-n}=(-1)^nJ_n$ となる。従って $Y_\nu$ の定義式の分子も整数で零になり、その点での極限は存在する。非整数の $\nu$ に対し二解 $J_\nu,J_{-\nu}$ の Wronskian の先頭を取ると $-2\nu/[x\Gamma(1+\nu)\Gamma(1-\nu)]=-2\sin\pi\nu/(\pi x)$ である。Abel の公式 :eqref[eq-mp-abel-wronskian] によって $x$ を掛けた値は一定である。従って $W(J_\nu,Y_\nu)=2/(\pi x)$ は次数の極限でも成立し、両者は独立である。例えば $J_0=1-x^2/4+O(x^4)$、$J_1=x/2-x^3/16+O(x^5)$ となる。

### 整数次の二解

整数次数で $J_{-m}$ が $J_m$ に比例するのは、逆 Gamma の零が負の指数の先頭係数を消すためである。第 2 種を得るには次数を整数でない値から近づける。次数微分では $(x/2)^\nu$ から $\log(x/2)$ が現れる。従って $Y_0$ の対数は重根の階数低下に現れる対数に対応する。

:::example[Bessel 級数の打ち切り（1）]
原点側で最初の補正を求める。

$$
m=1
$$

先頭項は Gamma 関数で規格化する。

$$
\begin{aligned}
J_m^{[0]}(X)\big|_{m=1} &= \frac{X}{2}
\end{aligned}
$$

続く二項は二乗ごとに出る。

$$
\begin{aligned}
J_m^{[2]}(X)\big|_{m=1} &= \frac{X \left(X^{4} - 24 X^{2} + 192\right)}{384}
\end{aligned}
$$

最初に捨てる項の係数を取り出す。

$$
\begin{aligned}
[X^{m+6}]J_m(X)\big|_{m=1} &= - \frac{1}{18432}
\end{aligned}
$$

$X=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
J_m^{[0]}(X)\big|_{X=\frac{1}{2}} &= \frac{1}{4} \\
\frac{\partial^1 (J_m^{[0]}(X))}{\partial X^1} &= \frac{1}{2} \\
\left.\frac{\partial^1 (J_m^{[0]}(X))}{\partial X^1}\right|_{X=\frac{1}{2}} &= \frac{1}{2} \\
\frac{\partial^2 (J_m^{[0]}(X))}{\partial X^2} &= 0 \\
\left.\frac{\partial^2 (J_m^{[0]}(X))}{\partial X^2}\right|_{X=\frac{1}{2}} &= 0
\end{aligned}
$$

$X=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
J_m^{[2]}(X)\big|_{X=\frac{1}{2}} &= \frac{2977}{12288} \\
\frac{\partial^1 (J_m^{[2]}(X))}{\partial X^1} &= \frac{5 X^{4}}{384} - \frac{3 X^{2}}{16} + \frac{1}{2} \\
\left.\frac{\partial^1 (J_m^{[2]}(X))}{\partial X^1}\right|_{X=\frac{1}{2}} &= \frac{2789}{6144} \\
\frac{\partial^2 (J_m^{[2]}(X))}{\partial X^2} &= \frac{X \left(5 X^{2} - 36\right)}{96} \\
\left.\frac{\partial^2 (J_m^{[2]}(X))}{\partial X^2}\right|_{X=\frac{1}{2}} &= - \frac{139}{768}
\end{aligned}
$$
:::

:::example[Bessel 級数の打ち切り（2）]
原点側で最初の補正を求める。次数を一つ上げる。

$$
m=2
$$

先頭項は Gamma 関数で規格化する。

$$
\begin{aligned}
J_m^{[0]}(X)\big|_{m=2} &= \frac{X^{2}}{8}
\end{aligned}
$$

続く二項は二乗ごとに出る。

$$
\begin{aligned}
J_m^{[2]}(X)\big|_{m=2} &= \frac{X^{2} \left(X^{4} - 32 X^{2} + 384\right)}{3072}
\end{aligned}
$$

最初に捨てる項の係数を取り出す。

$$
\begin{aligned}
[X^{m+6}]J_m(X)\big|_{m=2} &= - \frac{1}{184320}
\end{aligned}
$$

$X=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
J_m^{[0]}(X)\big|_{X=\frac{1}{2}} &= \frac{1}{32} \\
\frac{\partial^1 (J_m^{[0]}(X))}{\partial X^1} &= \frac{X}{4} \\
\left.\frac{\partial^1 (J_m^{[0]}(X))}{\partial X^1}\right|_{X=\frac{1}{2}} &= \frac{1}{8} \\
\frac{\partial^2 (J_m^{[0]}(X))}{\partial X^2} &= \frac{1}{4} \\
\left.\frac{\partial^2 (J_m^{[0]}(X))}{\partial X^2}\right|_{X=\frac{1}{2}} &= \frac{1}{4}
\end{aligned}
$$

$X=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
J_m^{[2]}(X)\big|_{X=\frac{1}{2}} &= \frac{6017}{196608} \\
\frac{\partial^1 (J_m^{[2]}(X))}{\partial X^1} &= \frac{X^{5}}{512} - \frac{X^{3}}{24} + \frac{X}{4} \\
\left.\frac{\partial^1 (J_m^{[2]}(X))}{\partial X^1}\right|_{X=\frac{1}{2}} &= \frac{5891}{49152} \\
\frac{\partial^2 (J_m^{[2]}(X))}{\partial X^2} &= \frac{5 X^{4}}{512} - \frac{X^{2}}{8} + \frac{1}{4} \\
\left.\frac{\partial^2 (J_m^{[2]}(X))}{\partial X^2}\right|_{X=\frac{1}{2}} &= \frac{1797}{8192}
\end{aligned}
$$
:::

:::example[Bessel 級数の打ち切り（3）]
原点側で最初の補正を求める。さらに次数を上げて、低次の結果と照合する。

$$
m=3
$$

先頭項は Gamma 関数で規格化する。

$$
\begin{aligned}
J_m^{[0]}(X)\big|_{m=3} &= \frac{X^{3}}{48}
\end{aligned}
$$

続く二項は二乗ごとに出る。

$$
\begin{aligned}
J_m^{[2]}(X)\big|_{m=3} &= \frac{X^{3} \left(X^{4} - 40 X^{2} + 640\right)}{30720}
\end{aligned}
$$

最初に捨てる項の係数を取り出す。

$$
\begin{aligned}
[X^{m+6}]J_m(X)\big|_{m=3} &= - \frac{1}{2211840}
\end{aligned}
$$

$X=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
J_m^{[0]}(X)\big|_{X=\frac{1}{2}} &= \frac{1}{384} \\
\frac{\partial^1 (J_m^{[0]}(X))}{\partial X^1} &= \frac{X^{2}}{16} \\
\left.\frac{\partial^1 (J_m^{[0]}(X))}{\partial X^1}\right|_{X=\frac{1}{2}} &= \frac{1}{64} \\
\frac{\partial^2 (J_m^{[0]}(X))}{\partial X^2} &= \frac{X}{8} \\
\left.\frac{\partial^2 (J_m^{[0]}(X))}{\partial X^2}\right|_{X=\frac{1}{2}} &= \frac{1}{16}
\end{aligned}
$$

$X=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
J_m^{[2]}(X)\big|_{X=\frac{1}{2}} &= \frac{10081}{3932160} \\
\frac{\partial^1 (J_m^{[2]}(X))}{\partial X^1} &= \frac{X^{2} \left(7 X^{4} - 200 X^{2} + 1920\right)}{30720} \\
\left.\frac{\partial^1 (J_m^{[2]}(X))}{\partial X^1}\right|_{X=\frac{1}{2}} &= \frac{29927}{1966080} \\
\frac{\partial^2 (J_m^{[2]}(X))}{\partial X^2} &= \frac{X \left(21 X^{4} - 400 X^{2} + 1920\right)}{15360} \\
\left.\frac{\partial^2 (J_m^{[2]}(X))}{\partial X^2}\right|_{X=\frac{1}{2}} &= \frac{29141}{491520}
\end{aligned}
$$
:::

:::example[Bessel 級数の打ち切り（4）]
原点側で最初の補正を求める。四次の結果まで具体化する。

$$
m=4
$$

先頭項は Gamma 関数で規格化する。

$$
\begin{aligned}
J_m^{[0]}(X)\big|_{m=4} &= \frac{X^{4}}{384}
\end{aligned}
$$

続く二項は二乗ごとに出る。

$$
\begin{aligned}
J_m^{[2]}(X)\big|_{m=4} &= \frac{X^{4} \left(X^{4} - 48 X^{2} + 960\right)}{368640}
\end{aligned}
$$

最初に捨てる項の係数を取り出す。

$$
\begin{aligned}
[X^{m+6}]J_m(X)\big|_{m=4} &= - \frac{1}{30965760}
\end{aligned}
$$

$X=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
J_m^{[0]}(X)\big|_{X=\frac{1}{2}} &= \frac{1}{6144} \\
\frac{\partial^1 (J_m^{[0]}(X))}{\partial X^1} &= \frac{X^{3}}{96} \\
\left.\frac{\partial^1 (J_m^{[0]}(X))}{\partial X^1}\right|_{X=\frac{1}{2}} &= \frac{1}{768} \\
\frac{\partial^2 (J_m^{[0]}(X))}{\partial X^2} &= \frac{X^{2}}{32} \\
\left.\frac{\partial^2 (J_m^{[0]}(X))}{\partial X^2}\right|_{X=\frac{1}{2}} &= \frac{1}{128}
\end{aligned}
$$

$X=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
J_m^{[2]}(X)\big|_{X=\frac{1}{2}} &= \frac{15169}{94371840} \\
\frac{\partial^1 (J_m^{[2]}(X))}{\partial X^1} &= \frac{X^{3} \left(X^{4} - 36 X^{2} + 480\right)}{46080} \\
\left.\frac{\partial^1 (J_m^{[2]}(X))}{\partial X^1}\right|_{X=\frac{1}{2}} &= \frac{7537}{5898240} \\
\frac{\partial^2 (J_m^{[2]}(X))}{\partial X^2} &= \frac{X^{2} \left(7 X^{4} - 180 X^{2} + 1440\right)}{46080} \\
\left.\frac{\partial^2 (J_m^{[2]}(X))}{\partial X^2}\right|_{X=\frac{1}{2}} &= \frac{22327}{2949120}
\end{aligned}
$$
:::

### 整数次数での特異項と規格化

$n=0$ では $Y_\nu$ の分子を次数で微分すると、$J_\nu$ と $-J_{-\nu}$ の両方から $\partial_\nu J_\nu|_{\nu=0}$ が一つずつ現れる。従って $Y_0=(2/\pi)\partial_\nu J_\nu|_{\nu=0}$ である。$H_j=\sum_{r=1}^j 1/r$、$H_0=0$ とするとディガンマ関数は $\psi(j+1)=H_j-\gamma$ なので、
$$
Y_0(x)=\frac2\pi\bigl(\log(x/2)+\gamma\bigr)J_0(x)
 -\frac2\pi\sum_{j=1}^{\infty}\frac{(-1)^j H_j(x/2)^{2j}}{(j!)^2}
$$
となる。級数の最初の項を取れば $Y_0=(2/\pi)(\log(x/2)+\gamma)+O(x^2|\log x|)$ である。

$n\ge1$ では $J_{-\nu}$ の $j=0$ の項を整数の極限で調べる。$1-\nu=(1-n)-(\nu-n)$ であり、$\Gamma$ の $1-n$ における留数は $(-1)^{n-1}/(n-1)!$ である。従って
$$
\lim_{\nu\to n}\frac{1}{\Gamma(1-\nu)\sin\pi\nu}
 =\frac{(n-1)!}{\pi}.
$$
$j=1,\ldots,n-1$ からは $x^{2j-n}$ の項が、残る添字と $J_\nu$ からは $O(x^n(1+|\log x|))$ の項が出る。従って $Y_n(x)\sim-(n-1)!(x/2)^{-n}/\pi$ である。例えば $Y_1(x)\sim-2/(\pi x)$ なので、円筒の軸で有界性を課すと $J_1$ 側が残る。

## 母関数・漸化式と積分表示

### 二つの指数を掛けて係数を集める

$e^{xt/2}e^{-x/(2t)}$ の二重級数は $t\ne0$ で絶対収束し、$t^m$ の係数は式 :eqref[eq-mp-bessel-series] である。従って
$$
e^{x(t-t^{-1})/2}=\sum_{m\in\Z}J_m(x)t^m .
\label{eq-mp-bessel-generating}
$$
$t=e^{i\theta}$ とすれば $J_m(x)=(2\pi)^{-1}\int_{-\pi}^\pi e^{i(x\sin\theta-m\theta)}\dd\theta$ となる。$x$ 微分と $t\partial_t$ をそれぞれ係数比較すれば $2J_m'=J_{m-1}-J_{m+1}$、$2mJ_m/x=J_{m-1}+J_{m+1}$ である。一般次数でも級数と Gamma の漸化式で同じ関係が示せる。例えば $J_0'=-J_1$、$(xJ_1)'=xJ_0$ となり、動径積分 $\int_0^a rJ_0(kr)\dd r=aJ_1(ka)/k$ が得られる。

### 非整数次数に使える Poisson 積分

$\Re\nu>-1/2$、$x>0$ で
$$
J_\nu(x)=\frac{(x/2)^\nu}{\sqrt\pi\Gamma(\nu+1/2)}
\int_{-1}^1 e^{ixt}(1-t^2)^{\nu-1/2}\dd t
\label{eq-mp-bessel-poisson}
$$
となる。指数を展開すると奇数項は零、$2j$ 次の積分は $t^2=u$ とすると $B(j+1/2,\nu+1/2)$ である。Gamma の倍角公式で $\Gamma(j+1/2)/(2j)!=\sqrt\pi/(4^jj!)$ だから、積分は式 :eqref[eq-mp-bessel-series] を再現する。例えば $\nu=1/2$ では重みが $1$ になり、積分は $2\sin x/x$ である。これに前の係数を掛けると半整数の明示式を得る。

### 指数の Fourier 係数

母関数に $t=e^{i\theta}$ を入れると $e^{ix\sin\theta}$ になる。指数の Taylor 展開を行い、$\sin^j\theta=((e^{i\theta}-e^{-i\theta})/(2i))^j$ と書くと、第 $m$ モードには $j=|m|,|m|+2,\ldots$ が寄与する。$\theta$ の閉区間で一様収束するので項別積分が許される。

:::example[母関数の最初の係数（1）]
正負のモードと次の冪を比べる。

$$
m=1
$$

第一次の寄与には正の指数だけを選ぶ。

$$
\begin{aligned}
[x^m]J_m(x)\big|_{m=1} &= \frac{1}{2}
\end{aligned}
$$

負のモードでは負の指数を選ぶたびに符号が変わる。

$$
\begin{aligned}
[x^m]J_{-m}(x)\big|_{m=1} &= - \frac{1}{2}
\end{aligned}
$$

次の寄与では $m+2$ 個のうち一つを負の指数にする。

$$
\begin{aligned}
[x^{m+2}]J_m(x)\big|_{m=1} &= - \frac{1}{16}
\end{aligned}
$$
:::

:::example[母関数の最初の係数（2）]
正負のモードと次の冪を比べる。次数を一つ上げる。

$$
m=2
$$

第一次の寄与には正の指数だけを選ぶ。

$$
\begin{aligned}
[x^m]J_m(x)\big|_{m=2} &= \frac{1}{8}
\end{aligned}
$$

負のモードでは負の指数を選ぶたびに符号が変わる。

$$
\begin{aligned}
[x^m]J_{-m}(x)\big|_{m=2} &= \frac{1}{8}
\end{aligned}
$$

次の寄与では $m+2$ 個のうち一つを負の指数にする。

$$
\begin{aligned}
[x^{m+2}]J_m(x)\big|_{m=2} &= - \frac{1}{96}
\end{aligned}
$$
:::

:::example[母関数の最初の係数（3）]
正負のモードと次の冪を比べる。さらに次数を上げて、低次の結果と照合する。

$$
m=3
$$

第一次の寄与には正の指数だけを選ぶ。

$$
\begin{aligned}
[x^m]J_m(x)\big|_{m=3} &= \frac{1}{48}
\end{aligned}
$$

負のモードでは負の指数を選ぶたびに符号が変わる。

$$
\begin{aligned}
[x^m]J_{-m}(x)\big|_{m=3} &= - \frac{1}{48}
\end{aligned}
$$

次の寄与では $m+2$ 個のうち一つを負の指数にする。

$$
\begin{aligned}
[x^{m+2}]J_m(x)\big|_{m=3} &= - \frac{1}{768}
\end{aligned}
$$
:::

:::example[母関数の最初の係数（4）]
正負のモードと次の冪を比べる。四次の結果まで具体化する。

$$
m=4
$$

第一次の寄与には正の指数だけを選ぶ。

$$
\begin{aligned}
[x^m]J_m(x)\big|_{m=4} &= \frac{1}{384}
\end{aligned}
$$

負のモードでは負の指数を選ぶたびに符号が変わる。

$$
\begin{aligned}
[x^m]J_{-m}(x)\big|_{m=4} &= \frac{1}{384}
\end{aligned}
$$

次の寄与では $m+2$ 個のうち一つを負の指数にする。

$$
\begin{aligned}
[x^{m+2}]J_m(x)\big|_{m=4} &= - \frac{1}{7680}
\end{aligned}
$$
:::

## Hankel 関数と変形 Bessel 関数

### 進行波と減衰する波

:term[ハンケル関数]{en="Hankel function" yomi="はんけるかんすう"} は $H_\nu^{(1,2)}=J_\nu\pm iY_\nu$ とする。時間因子 $e^{-i\omega t}$ では大きな半径で $e^{ikr}$ を含む $H^{(1)}$ が外へ進む。一方、変形方程式 $x^2y''+xy'-(x^2+\nu^2)y=0$ は引数を $ix$ に替えて得られる。$I_\nu=e^{-i\pi\nu/2}J_\nu(ix)$ では式 :eqref[eq-mp-bessel-series] の因子 $(-1)^j$ と $i^{2j}$ が打ち消し合い、$K_\nu=\pi(I_{-\nu}-I_\nu)/(2\sin\pi\nu)$ は整数で次数の極限を取る。

減衰側の積分表示 $K_\nu(x)=\int_0^\infty e^{-x\cosh t}\cosh(\nu t)\dd t$ は、被積分関数の二階 $t$ 微分を二回部分積分すると変形方程式を満たす。$x\to+\infty$ の Laplace 評価で $\sqrt{\pi/(2x)}e^{-x}$ になる。この積分を $u=e^t$ で両側に延ばすと $\tfrac12\int_0^\infty e^{-x(u+1/u)/2}u^{\nu-1}\dd u$ である。$0<\nu<1$ では $x\downarrow0$ に $x^\nu K\to2^{\nu-1}\Gamma(\nu)$ となる。$v=xu/2$ と置き、$e^{-x^2/(4v)}$ と $1$ の差を $v=x^2/(4w)$ で積分して次の係数を取ると $2^{-\nu-1}\Gamma(-\nu)x^\nu$ を得る。実際、残る積分は $\int_0^\infty w^{-\nu-1}(e^{-w}-1)\dd w=\Gamma(-\nu)$ で、部分積分と優収束が使える。反射公式によりこの二係数は上の $I_{-\nu},I_\nu$ の組合せと一致する。次数の解析接続で他の実次数にも成立する。積分の $\cosh(\nu t)$ は次数のコンパクト集合で $C e^{M t}$ と抑えられ、$e^{-x\cosh t}$ がこれより速く減衰するので積分は次数について整関数である。

### 変形 Bessel 関数の漸化式と Wronskian

$I_\nu$ の級数を微分して、逆 Gamma の漸化式で係数をそろえると
$$
 I_\nu'-\frac{\nu}{x}I_\nu=I_{\nu+1},\qquad
 I_\nu'+\frac{\nu}{x}I_\nu=I_{\nu-1}
$$
となる。第 1 式では $2j$ が出て添字 $j$ が一つ下がり、第 2 式では $2(j+\nu)$ が出て Gamma の引数が一つ下がる。一方、$K_\nu$ の積分を $x$ で微分すると $\cosh t$ が掛かるため、積の加法公式により $K_\nu'=-(K_{\nu+1}+K_{\nu-1})/2$ である。また
$$
 K_{\nu+1}-K_{\nu-1}
 =2\int_0^\infty e^{-x\cosh t}\sinh(\nu t)\sinh t\dd t
 =\frac{2\nu}{x}K_\nu .
$$
最後の等号には $-x\sinh t\,e^{-x\cosh t}=\partial_t e^{-x\cosh t}$ と部分積分を使った。原点では $\sinh(\nu t)=0$、無限遠では指数減衰するので境界項は消える。例えば $K_{3/2}=K_{-1/2}+x^{-1}K_{1/2}$ であり、$K_{-\nu}=K_\nu$ と半整数の初期値から $K_{3/2}=\sqrt{\pi/(2x)}e^{-x}(1+1/x)$ である。繰り返せば全ての正半整数について指数と $1/x$ の有限多項式を得る。

非整数の $\nu$ に対し、$I_\nu,I_{-\nu}$ の先頭の冪から $x W(I_\nu,I_{-\nu})=-2\sin(\pi\nu)/\pi$ を得る。変形方程式の一階係数も $1/x$ なので、Abel の公式によりこれは $x$ によらない。定義の線形結合に戻すと $W(I_\nu,K_\nu)=-1/x$ となり、次数で接続すれば整数でも成立する。この符号は、増大側 $I_\nu$ を $K_\nu$ より先に取った順序に対応する。

### 変形方程式の正の解

引数を $ix$ にすると $x^2$ 項の符号が変わる。次数の冪から $i^\nu$ を除いた第 1 解が $I_\nu$ である。$\nu\ge0$、$x>0$ なら級数の各項は正なので $I_\nu>0$ である。正の実軸で指数的に減衰する方は $K_\nu$ である。半整数では $1/x$ の有限多項式に指数と平方根が掛かる。

:::example[変形 Bessel の原点側（1）]
正の整数次数で $I_m$ の増加が各補正にどう現れるか調べる。

$$
m=1
$$

第 1 種 $J$ と同じ先頭を持つ。

$$
\begin{aligned}
I_m^{[0]}(X)\big|_{m=1} &= \frac{X}{2}
\end{aligned}
$$

二乗の補正を二つ残した式を作る。

$$
\begin{aligned}
I_m^{[2]}(X)\big|_{m=1} &= \frac{X \left(X^{4} + 24 X^{2} + 192\right)}{384}
\end{aligned}
$$

先頭で割った比を求める。

$$
\begin{aligned}
I_m^{[2]}/I_m^{[0]}\big|_{m=1} &= \frac{X^{4}}{192} + \frac{X^{2}}{8} + 1
\end{aligned}
$$

$X=2$ で値と導関数を評価する。

$$
\begin{aligned}
I_m^{[0]}(X)\big|_{X=2} &= 1 \\
\frac{\partial^1 (I_m^{[0]}(X))}{\partial X^1} &= \frac{1}{2} \\
\left.\frac{\partial^1 (I_m^{[0]}(X))}{\partial X^1}\right|_{X=2} &= \frac{1}{2}
\end{aligned}
$$

$X=2$ で値と導関数を評価する。

$$
\begin{aligned}
I_m^{[2]}(X)\big|_{X=2} &= \frac{19}{12} \\
\frac{\partial^1 (I_m^{[2]}(X))}{\partial X^1} &= \frac{5 X^{4}}{384} + \frac{3 X^{2}}{16} + \frac{1}{2} \\
\left.\frac{\partial^1 (I_m^{[2]}(X))}{\partial X^1}\right|_{X=2} &= \frac{35}{24}
\end{aligned}
$$

$X=2$ で値と導関数を評価する。

$$
\begin{aligned}
I_m^{[2]}/I_m^{[0]}\big|_{X=2} &= \frac{19}{12} \\
\frac{\partial^1 (I_m^{[2]}/I_m^{[0]})}{\partial X^1} &= \frac{X \left(X^{2} + 12\right)}{48} \\
\left.\frac{\partial^1 (I_m^{[2]}/I_m^{[0]})}{\partial X^1}\right|_{X=2} &= \frac{2}{3}
\end{aligned}
$$
:::

:::example[変形 Bessel の原点側（2）]
正の整数次数で $I_m$ の増加が各補正にどう現れるか調べる。次数を一つ上げる。

$$
m=2
$$

第 1 種 $J$ と同じ先頭を持つ。

$$
\begin{aligned}
I_m^{[0]}(X)\big|_{m=2} &= \frac{X^{2}}{8}
\end{aligned}
$$

二乗の補正を二つ残した式を作る。

$$
\begin{aligned}
I_m^{[2]}(X)\big|_{m=2} &= \frac{X^{2} \left(X^{4} + 32 X^{2} + 384\right)}{3072}
\end{aligned}
$$

先頭で割った比を求める。

$$
\begin{aligned}
I_m^{[2]}/I_m^{[0]}\big|_{m=2} &= \frac{X^{4}}{384} + \frac{X^{2}}{12} + 1
\end{aligned}
$$

$X=2$ で値と導関数を評価する。

$$
\begin{aligned}
I_m^{[0]}(X)\big|_{X=2} &= \frac{1}{2} \\
\frac{\partial^1 (I_m^{[0]}(X))}{\partial X^1} &= \frac{X}{4} \\
\left.\frac{\partial^1 (I_m^{[0]}(X))}{\partial X^1}\right|_{X=2} &= \frac{1}{2}
\end{aligned}
$$

$X=2$ で値と導関数を評価する。

$$
\begin{aligned}
I_m^{[2]}(X)\big|_{X=2} &= \frac{11}{16} \\
\frac{\partial^1 (I_m^{[2]}(X))}{\partial X^1} &= \frac{X^{5}}{512} + \frac{X^{3}}{24} + \frac{X}{4} \\
\left.\frac{\partial^1 (I_m^{[2]}(X))}{\partial X^1}\right|_{X=2} &= \frac{43}{48}
\end{aligned}
$$

$X=2$ で値と導関数を評価する。

$$
\begin{aligned}
I_m^{[2]}/I_m^{[0]}\big|_{X=2} &= \frac{11}{8} \\
\frac{\partial^1 (I_m^{[2]}/I_m^{[0]})}{\partial X^1} &= \frac{X \left(X^{2} + 16\right)}{96} \\
\left.\frac{\partial^1 (I_m^{[2]}/I_m^{[0]})}{\partial X^1}\right|_{X=2} &= \frac{5}{12}
\end{aligned}
$$
:::

:::example[変形 Bessel の原点側（3）]
正の整数次数で $I_m$ の増加が各補正にどう現れるか調べる。さらに次数を上げて、低次の結果と照合する。

$$
m=3
$$

第 1 種 $J$ と同じ先頭を持つ。

$$
\begin{aligned}
I_m^{[0]}(X)\big|_{m=3} &= \frac{X^{3}}{48}
\end{aligned}
$$

二乗の補正を二つ残した式を作る。

$$
\begin{aligned}
I_m^{[2]}(X)\big|_{m=3} &= \frac{X^{3} \left(X^{4} + 40 X^{2} + 640\right)}{30720}
\end{aligned}
$$

先頭で割った比を求める。

$$
\begin{aligned}
I_m^{[2]}/I_m^{[0]}\big|_{m=3} &= \frac{X^{4}}{640} + \frac{X^{2}}{16} + 1
\end{aligned}
$$

$X=2$ で値と導関数を評価する。

$$
\begin{aligned}
I_m^{[0]}(X)\big|_{X=2} &= \frac{1}{6} \\
\frac{\partial^1 (I_m^{[0]}(X))}{\partial X^1} &= \frac{X^{2}}{16} \\
\left.\frac{\partial^1 (I_m^{[0]}(X))}{\partial X^1}\right|_{X=2} &= \frac{1}{4}
\end{aligned}
$$

$X=2$ で値と導関数を評価する。

$$
\begin{aligned}
I_m^{[2]}(X)\big|_{X=2} &= \frac{17}{80} \\
\frac{\partial^1 (I_m^{[2]}(X))}{\partial X^1} &= \frac{X^{2} \left(7 X^{4} + 200 X^{2} + 1920\right)}{30720} \\
\left.\frac{\partial^1 (I_m^{[2]}(X))}{\partial X^1}\right|_{X=2} &= \frac{59}{160}
\end{aligned}
$$

$X=2$ で値と導関数を評価する。

$$
\begin{aligned}
I_m^{[2]}/I_m^{[0]}\big|_{X=2} &= \frac{51}{40} \\
\frac{\partial^1 (I_m^{[2]}/I_m^{[0]})}{\partial X^1} &= \frac{X \left(X^{2} + 20\right)}{160} \\
\left.\frac{\partial^1 (I_m^{[2]}/I_m^{[0]})}{\partial X^1}\right|_{X=2} &= \frac{3}{10}
\end{aligned}
$$
:::

:::example[変形 Bessel の原点側（4）]
正の整数次数で $I_m$ の増加が各補正にどう現れるか調べる。四次の結果まで具体化する。

$$
m=4
$$

第 1 種 $J$ と同じ先頭を持つ。

$$
\begin{aligned}
I_m^{[0]}(X)\big|_{m=4} &= \frac{X^{4}}{384}
\end{aligned}
$$

二乗の補正を二つ残した式を作る。

$$
\begin{aligned}
I_m^{[2]}(X)\big|_{m=4} &= \frac{X^{4} \left(X^{4} + 48 X^{2} + 960\right)}{368640}
\end{aligned}
$$

先頭で割った比を求める。

$$
\begin{aligned}
I_m^{[2]}/I_m^{[0]}\big|_{m=4} &= \frac{X^{4}}{960} + \frac{X^{2}}{20} + 1
\end{aligned}
$$

$X=2$ で値と導関数を評価する。

$$
\begin{aligned}
I_m^{[0]}(X)\big|_{X=2} &= \frac{1}{24} \\
\frac{\partial^1 (I_m^{[0]}(X))}{\partial X^1} &= \frac{X^{3}}{96} \\
\left.\frac{\partial^1 (I_m^{[0]}(X))}{\partial X^1}\right|_{X=2} &= \frac{1}{12}
\end{aligned}
$$

$X=2$ で値と導関数を評価する。

$$
\begin{aligned}
I_m^{[2]}(X)\big|_{X=2} &= \frac{73}{1440} \\
\frac{\partial^1 (I_m^{[2]}(X))}{\partial X^1} &= \frac{X^{3} \left(X^{4} + 36 X^{2} + 480\right)}{46080} \\
\left.\frac{\partial^1 (I_m^{[2]}(X))}{\partial X^1}\right|_{X=2} &= \frac{1}{9}
\end{aligned}
$$

$X=2$ で値と導関数を評価する。

$$
\begin{aligned}
I_m^{[2]}/I_m^{[0]}\big|_{X=2} &= \frac{73}{60} \\
\frac{\partial^1 (I_m^{[2]}/I_m^{[0]})}{\partial X^1} &= \frac{X \left(X^{2} + 24\right)}{240} \\
\left.\frac{\partial^1 (I_m^{[2]}/I_m^{[0]})}{\partial X^1}\right|_{X=2} &= \frac{7}{30}
\end{aligned}
$$
:::

:::example[半整数次数 1-1/2 の方程式]
この次数では負の指数に掛かる有限和を作り、変形 Bessel の方程式に戻す。係数を読み、試行関数を微分して方程式へ戻す。

この例の係数を左から二階、一階、零階の順に書く。右辺がある場合は三つの項の和と比較する。

$$
\left( X^{2}, \  X, \  - X^{2} - \frac{1}{4}\right)
$$

関数とその二つの導関数を並べる。特に一つ下の階数の導関数には、一階項の係数を掛ける。

$$
\begin{aligned}u&=\frac{\sqrt{2} \sqrt{\pi} e^{- X}}{2 \sqrt{X}}\\u^{\prime}&=- \frac{\sqrt{2} \sqrt{\pi} e^{- X}}{2 \sqrt{X}} - \frac{\sqrt{2} \sqrt{\pi} e^{- X}}{4 X^{\frac{3}{2}}}\\u^{\prime\prime}&=\frac{\sqrt{2} \sqrt{\pi} \left(4 + \frac{4}{X} + \frac{3}{X^{2}}\right) e^{- X}}{8 \sqrt{X}}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=\frac{\sqrt{2} \sqrt{\pi} X^{\frac{3}{2}} \left(4 + \frac{4}{X} + \frac{3}{X^{2}}\right) e^{- X}}{8}+\left(X \left(- \frac{\sqrt{2} \sqrt{\pi} e^{- X}}{2 \sqrt{X}} - \frac{\sqrt{2} \sqrt{\pi} e^{- X}}{4 X^{\frac{3}{2}}}\right)\right)+\left(\frac{\sqrt{2} \sqrt{\pi} \left(- X^{2} - \frac{1}{4}\right) e^{- X}}{2 \sqrt{X}}\right)\\&=0\end{aligned}
$$

$X$ が $2$ のところで値を評価する。変換した方程式を使う場合も、初期値や境界値は元の変数で照合する。

$$
\begin{aligned}u&=\frac{\sqrt{\pi}}{2 e^{2}}\\ u^\prime&=- \frac{5 \sqrt{\pi}}{8 e^{2}}\\ Lu&=0\end{aligned}
$$
:::
:::example[半整数次数 2-1/2 の方程式]
この次数では負の指数に掛かる有限和を作り、変形 Bessel の方程式に戻す。係数を読み、試行関数を微分して方程式へ戻す。

この例の係数を左から二階、一階、零階の順に書く。右辺がある場合は三つの項の和と比較する。

$$
\left( X^{2}, \  X, \  - X^{2} - \frac{9}{4}\right)
$$

関数とその二つの導関数を並べる。特に一つ下の階数の導関数には、一階項の係数を掛ける。

$$
\begin{aligned}u&=\frac{\sqrt{2} \sqrt{\pi} e^{- X}}{2 \sqrt{X}} + \frac{\sqrt{2} \sqrt{\pi} e^{- X}}{2 X^{\frac{3}{2}}}\\u^{\prime}&=- \frac{\sqrt{2} \sqrt{\pi} e^{- X}}{2 \sqrt{X}} - \frac{3 \sqrt{2} \sqrt{\pi} e^{- X}}{4 X^{\frac{3}{2}}} - \frac{3 \sqrt{2} \sqrt{\pi} e^{- X}}{4 X^{\frac{5}{2}}}\\u^{\prime\prime}&=\frac{\sqrt{2} \sqrt{\pi} \left(\frac{1}{2} + \frac{1}{X} + \frac{15}{8 X^{2}} + \frac{15}{8 X^{3}}\right) e^{- X}}{\sqrt{X}}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=\sqrt{2} \sqrt{\pi} X^{\frac{3}{2}} \left(\frac{1}{2} + \frac{1}{X} + \frac{15}{8 X^{2}} + \frac{15}{8 X^{3}}\right) e^{- X}+\left(X \left(- \frac{\sqrt{2} \sqrt{\pi} e^{- X}}{2 \sqrt{X}} - \frac{3 \sqrt{2} \sqrt{\pi} e^{- X}}{4 X^{\frac{3}{2}}} - \frac{3 \sqrt{2} \sqrt{\pi} e^{- X}}{4 X^{\frac{5}{2}}}\right)\right)+\left(\left(- X^{2} - \frac{9}{4}\right) \left(\frac{\sqrt{2} \sqrt{\pi} e^{- X}}{2 \sqrt{X}} + \frac{\sqrt{2} \sqrt{\pi} e^{- X}}{2 X^{\frac{3}{2}}}\right)\right)\\&=0\end{aligned}
$$

$X$ が $2$ のところで値を評価する。変換した方程式を使う場合も、初期値や境界値は元の変数で照合する。

$$
\begin{aligned}u&=\frac{3 \sqrt{\pi}}{4 e^{2}}\\ u^\prime&=- \frac{17 \sqrt{\pi}}{16 e^{2}}\\ Lu&=0\end{aligned}
$$
:::
:::example[半整数次数 3-1/2 の方程式]
この次数では負の指数に掛かる有限和を作り、変形 Bessel の方程式に戻す。係数を読み、試行関数を微分して方程式へ戻す。

この例の係数を左から二階、一階、零階の順に書く。右辺がある場合は三つの項の和と比較する。

$$
\left( X^{2}, \  X, \  - X^{2} - \frac{25}{4}\right)
$$

関数とその二つの導関数を並べる。特に一つ下の階数の導関数には、一階項の係数を掛ける。

$$
\begin{aligned}u&=\frac{\sqrt{2} \sqrt{\pi} e^{- X}}{2 \sqrt{X}} + \frac{3 \sqrt{2} \sqrt{\pi} e^{- X}}{2 X^{\frac{3}{2}}} + \frac{3 \sqrt{2} \sqrt{\pi} e^{- X}}{2 X^{\frac{5}{2}}}\\u^{\prime}&=- \frac{\sqrt{2} \sqrt{\pi} e^{- X}}{2 \sqrt{X}} - \frac{7 \sqrt{2} \sqrt{\pi} e^{- X}}{4 X^{\frac{3}{2}}} - \frac{15 \sqrt{2} \sqrt{\pi} e^{- X}}{4 X^{\frac{5}{2}}} - \frac{15 \sqrt{2} \sqrt{\pi} e^{- X}}{4 X^{\frac{7}{2}}}\\u^{\prime\prime}&=\frac{\sqrt{2} \sqrt{\pi} \left(\frac{1}{2} + \frac{2}{X} + \frac{51}{8 X^{2}} + \frac{105}{8 X^{3}} + \frac{105}{8 X^{4}}\right) e^{- X}}{\sqrt{X}}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=\sqrt{2} \sqrt{\pi} X^{\frac{3}{2}} \left(\frac{1}{2} + \frac{2}{X} + \frac{51}{8 X^{2}} + \frac{105}{8 X^{3}} + \frac{105}{8 X^{4}}\right) e^{- X}+\left(X \left(- \frac{\sqrt{2} \sqrt{\pi} e^{- X}}{2 \sqrt{X}} - \frac{7 \sqrt{2} \sqrt{\pi} e^{- X}}{4 X^{\frac{3}{2}}} - \frac{15 \sqrt{2} \sqrt{\pi} e^{- X}}{4 X^{\frac{5}{2}}} - \frac{15 \sqrt{2} \sqrt{\pi} e^{- X}}{4 X^{\frac{7}{2}}}\right)\right)+\left(\left(- X^{2} - \frac{25}{4}\right) \left(\frac{\sqrt{2} \sqrt{\pi} e^{- X}}{2 \sqrt{X}} + \frac{3 \sqrt{2} \sqrt{\pi} e^{- X}}{2 X^{\frac{3}{2}}} + \frac{3 \sqrt{2} \sqrt{\pi} e^{- X}}{2 X^{\frac{5}{2}}}\right)\right)\\&=0\end{aligned}
$$

$X$ が $2$ のところで値を評価する。変換した方程式を使う場合も、初期値や境界値は元の変数で照合する。

$$
\begin{aligned}u&=\frac{13 \sqrt{\pi}}{8 e^{2}}\\ u^\prime&=- \frac{89 \sqrt{\pi}}{32 e^{2}}\\ Lu&=0\end{aligned}
$$
:::
:::example[半整数次数 4-1/2 の方程式]
この次数では負の指数に掛かる有限和を作り、変形 Bessel の方程式に戻す。係数を読み、試行関数を微分して方程式へ戻す。

この例の係数を左から二階、一階、零階の順に書く。右辺がある場合は三つの項の和と比較する。

$$
\left( X^{2}, \  X, \  - X^{2} - \frac{49}{4}\right)
$$

関数とその二つの導関数を並べる。特に一つ下の階数の導関数には、一階項の係数を掛ける。

$$
\begin{aligned}u&=\frac{\sqrt{2} \sqrt{\pi} e^{- X}}{2 \sqrt{X}} + \frac{3 \sqrt{2} \sqrt{\pi} e^{- X}}{X^{\frac{3}{2}}} + \frac{15 \sqrt{2} \sqrt{\pi} e^{- X}}{2 X^{\frac{5}{2}}} + \frac{15 \sqrt{2} \sqrt{\pi} e^{- X}}{2 X^{\frac{7}{2}}}\\u^{\prime}&=- \frac{\sqrt{2} \sqrt{\pi} e^{- X}}{2 \sqrt{X}} - \frac{13 \sqrt{2} \sqrt{\pi} e^{- X}}{4 X^{\frac{3}{2}}} - \frac{12 \sqrt{2} \sqrt{\pi} e^{- X}}{X^{\frac{5}{2}}} - \frac{105 \sqrt{2} \sqrt{\pi} e^{- X}}{4 X^{\frac{7}{2}}} - \frac{105 \sqrt{2} \sqrt{\pi} e^{- X}}{4 X^{\frac{9}{2}}}\\u^{\prime\prime}&=\frac{\sqrt{2} \sqrt{\pi} \left(4 + \frac{28}{X} + \frac{135}{X^{2}} + \frac{450}{X^{3}} + \frac{945}{X^{4}} + \frac{945}{X^{5}}\right) e^{- X}}{8 \sqrt{X}}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=\frac{\sqrt{2} \sqrt{\pi} X^{\frac{3}{2}} \left(4 + \frac{28}{X} + \frac{135}{X^{2}} + \frac{450}{X^{3}} + \frac{945}{X^{4}} + \frac{945}{X^{5}}\right) e^{- X}}{8}+\left(X \left(- \frac{\sqrt{2} \sqrt{\pi} e^{- X}}{2 \sqrt{X}} - \frac{13 \sqrt{2} \sqrt{\pi} e^{- X}}{4 X^{\frac{3}{2}}} - \frac{12 \sqrt{2} \sqrt{\pi} e^{- X}}{X^{\frac{5}{2}}} - \frac{105 \sqrt{2} \sqrt{\pi} e^{- X}}{4 X^{\frac{7}{2}}} - \frac{105 \sqrt{2} \sqrt{\pi} e^{- X}}{4 X^{\frac{9}{2}}}\right)\right)+\left(\left(- X^{2} - \frac{49}{4}\right) \left(\frac{\sqrt{2} \sqrt{\pi} e^{- X}}{2 \sqrt{X}} + \frac{3 \sqrt{2} \sqrt{\pi} e^{- X}}{X^{\frac{3}{2}}} + \frac{15 \sqrt{2} \sqrt{\pi} e^{- X}}{2 X^{\frac{5}{2}}} + \frac{15 \sqrt{2} \sqrt{\pi} e^{- X}}{2 X^{\frac{7}{2}}}\right)\right)\\&=0\end{aligned}
$$

$X$ が $2$ のところで値を評価する。変換した方程式を使う場合も、初期値や境界値は元の変数で照合する。

$$
\begin{aligned}u&=\frac{77 \sqrt{\pi}}{16 e^{2}}\\ u^\prime&=- \frac{643 \sqrt{\pi}}{64 e^{2}}\\ Lu&=0\end{aligned}
$$
:::
:::example[半整数次数 5-1/2 の方程式]
この次数では負の指数に掛かる有限和を作り、変形 Bessel の方程式に戻す。係数を読み、試行関数を微分して方程式へ戻す。

この例の係数を左から二階、一階、零階の順に書く。右辺がある場合は三つの項の和と比較する。

$$
\left( X^{2}, \  X, \  - X^{2} - \frac{81}{4}\right)
$$

関数とその二つの導関数を並べる。特に一つ下の階数の導関数には、一階項の係数を掛ける。

$$
\begin{aligned}u&=\frac{\sqrt{2} \sqrt{\pi} e^{- X}}{2 \sqrt{X}} + \frac{5 \sqrt{2} \sqrt{\pi} e^{- X}}{X^{\frac{3}{2}}} + \frac{45 \sqrt{2} \sqrt{\pi} e^{- X}}{2 X^{\frac{5}{2}}} + \frac{105 \sqrt{2} \sqrt{\pi} e^{- X}}{2 X^{\frac{7}{2}}} + \frac{105 \sqrt{2} \sqrt{\pi} e^{- X}}{2 X^{\frac{9}{2}}}\\u^{\prime}&=- \frac{\sqrt{2} \sqrt{\pi} e^{- X}}{2 \sqrt{X}} - \frac{21 \sqrt{2} \sqrt{\pi} e^{- X}}{4 X^{\frac{3}{2}}} - \frac{30 \sqrt{2} \sqrt{\pi} e^{- X}}{X^{\frac{5}{2}}} - \frac{435 \sqrt{2} \sqrt{\pi} e^{- X}}{4 X^{\frac{7}{2}}} - \frac{945 \sqrt{2} \sqrt{\pi} e^{- X}}{4 X^{\frac{9}{2}}} - \frac{945 \sqrt{2} \sqrt{\pi} e^{- X}}{4 X^{\frac{11}{2}}}\\u^{\prime\prime}&=\frac{\sqrt{2} \sqrt{\pi} \left(4 + \frac{44}{X} + \frac{303}{X^{2}} + \frac{1470}{X^{3}} + \frac{4935}{X^{4}} + \frac{10395}{X^{5}} + \frac{10395}{X^{6}}\right) e^{- X}}{8 \sqrt{X}}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=\frac{\sqrt{2} \sqrt{\pi} X^{\frac{3}{2}} \left(4 + \frac{44}{X} + \frac{303}{X^{2}} + \frac{1470}{X^{3}} + \frac{4935}{X^{4}} + \frac{10395}{X^{5}} + \frac{10395}{X^{6}}\right) e^{- X}}{8}+\left(X \left(- \frac{\sqrt{2} \sqrt{\pi} e^{- X}}{2 \sqrt{X}} - \frac{21 \sqrt{2} \sqrt{\pi} e^{- X}}{4 X^{\frac{3}{2}}} - \frac{30 \sqrt{2} \sqrt{\pi} e^{- X}}{X^{\frac{5}{2}}} - \frac{435 \sqrt{2} \sqrt{\pi} e^{- X}}{4 X^{\frac{7}{2}}} - \frac{945 \sqrt{2} \sqrt{\pi} e^{- X}}{4 X^{\frac{9}{2}}} - \frac{945 \sqrt{2} \sqrt{\pi} e^{- X}}{4 X^{\frac{11}{2}}}\right)\right)+\left(\left(- X^{2} - \frac{81}{4}\right) \left(\frac{\sqrt{2} \sqrt{\pi} e^{- X}}{2 \sqrt{X}} + \frac{5 \sqrt{2} \sqrt{\pi} e^{- X}}{X^{\frac{3}{2}}} + \frac{45 \sqrt{2} \sqrt{\pi} e^{- X}}{2 X^{\frac{5}{2}}} + \frac{105 \sqrt{2} \sqrt{\pi} e^{- X}}{2 X^{\frac{7}{2}}} + \frac{105 \sqrt{2} \sqrt{\pi} e^{- X}}{2 X^{\frac{9}{2}}}\right)\right)\\&=0\end{aligned}
$$

$X$ が $2$ のところで値を評価する。変換した方程式を使う場合も、初期値や境界値は元の変数で照合する。

$$
\begin{aligned}u&=\frac{591 \sqrt{\pi}}{32 e^{2}}\\ u^\prime&=- \frac{5935 \sqrt{\pi}}{128 e^{2}}\\ Lu&=0\end{aligned}
$$
:::

## 球 Bessel 関数

### 球面の角固有値を入れる

三次元の $\Delta+k^2$ の動径方程式は $R''+2R'/r+[k^2-\ell(\ell+1)/r^2]R=0$ である。$x=kr$、$R=x^{-1/2}u$ とすると次数 $\ell+1/2$ の円筒方程式になる。そこで $j_\ell=\sqrt{\pi/(2x)}J_{\ell+1/2}$、$y_\ell=\sqrt{\pi/(2x)}Y_{\ell+1/2}$ とする。Gamma の半整数の漸化式から $J_{1/2}=\sqrt{2/(\pi x)}\sin x$ と $J_{-1/2}=\sqrt{2/(\pi x)}\cos x$ が従う。よって $j_0=\sin x/x$、$y_0=-\cos x/x$、さらに $j_1=-j_0'=\sin x/x^2-\cos x/x$ である。$R$ が中心で滑らかであるためには $j_\ell$ 側を選ぶ。$y_\ell$ の先頭は $x^{-\ell-1}$ であるため除く。

### 球面波の原点規格化

三次元では動径関数に $x$ を掛けると一階微分が消え、遠心項が $\ell(\ell+1)/x^2$ になる。$x$ を掛ける前の原点指数は $\ell,-\ell-1$ である。滑らかな解を半整数の $J_{\ell+1/2}$ と同一視すると先頭係数は $1/(2\ell+1)!!$ になる。

:::example[球 Bessel 関数の原点付近（1）]
$x^\ell$ を除いたあとの補正を調べる。

$$
m=1
$$

最初の項を求める。

$$
\begin{aligned}
j_m^{[0]}(X)\big|_{m=1} &= \frac{X}{3}
\end{aligned}
$$

二乗の補正を加える。

$$
\begin{aligned}
j_m^{[1]}(X)\big|_{m=1} &= \frac{X \left(10 - X^{2}\right)}{30}
\end{aligned}
$$

先頭に対する補正の比を取る。

$$
\begin{aligned}
(j_m^{[1]}-j_m^{[0]})/j_m^{[0]}\big|_{m=1} &= - \frac{X^{2}}{10}
\end{aligned}
$$

$X=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
j_m^{[0]}(X)\big|_{X=\frac{1}{2}} &= \frac{1}{6} \\
\frac{\partial^1 (j_m^{[0]}(X))}{\partial X^1} &= \frac{1}{3} \\
\left.\frac{\partial^1 (j_m^{[0]}(X))}{\partial X^1}\right|_{X=\frac{1}{2}} &= \frac{1}{3} \\
\frac{\partial^2 (j_m^{[0]}(X))}{\partial X^2} &= 0 \\
\left.\frac{\partial^2 (j_m^{[0]}(X))}{\partial X^2}\right|_{X=\frac{1}{2}} &= 0
\end{aligned}
$$

$X=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
j_m^{[1]}(X)\big|_{X=\frac{1}{2}} &= \frac{13}{80} \\
\frac{\partial^1 (j_m^{[1]}(X))}{\partial X^1} &= \frac{1}{3} - \frac{X^{2}}{10} \\
\left.\frac{\partial^1 (j_m^{[1]}(X))}{\partial X^1}\right|_{X=\frac{1}{2}} &= \frac{37}{120} \\
\frac{\partial^2 (j_m^{[1]}(X))}{\partial X^2} &= - \frac{X}{5} \\
\left.\frac{\partial^2 (j_m^{[1]}(X))}{\partial X^2}\right|_{X=\frac{1}{2}} &= - \frac{1}{10}
\end{aligned}
$$

$X=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
(j_m^{[1]}-j_m^{[0]})/j_m^{[0]}\big|_{X=\frac{1}{2}} &= - \frac{1}{40} \\
\frac{\partial^1 ((j_m^{[1]}-j_m^{[0]})/j_m^{[0]})}{\partial X^1} &= - \frac{X}{5} \\
\left.\frac{\partial^1 ((j_m^{[1]}-j_m^{[0]})/j_m^{[0]})}{\partial X^1}\right|_{X=\frac{1}{2}} &= - \frac{1}{10} \\
\frac{\partial^2 ((j_m^{[1]}-j_m^{[0]})/j_m^{[0]})}{\partial X^2} &= - \frac{1}{5} \\
\left.\frac{\partial^2 ((j_m^{[1]}-j_m^{[0]})/j_m^{[0]})}{\partial X^2}\right|_{X=\frac{1}{2}} &= - \frac{1}{5}
\end{aligned}
$$
:::

:::example[球 Bessel 関数の原点付近（2）]
$x^\ell$ を除いたあとの補正を調べる。次数を一つ上げる。

$$
m=2
$$

最初の項を求める。

$$
\begin{aligned}
j_m^{[0]}(X)\big|_{m=2} &= \frac{X^{2}}{15}
\end{aligned}
$$

二乗の補正を加える。

$$
\begin{aligned}
j_m^{[1]}(X)\big|_{m=2} &= \frac{X^{2} \left(14 - X^{2}\right)}{210}
\end{aligned}
$$

先頭に対する補正の比を取る。

$$
\begin{aligned}
(j_m^{[1]}-j_m^{[0]})/j_m^{[0]}\big|_{m=2} &= - \frac{X^{2}}{14}
\end{aligned}
$$

$X=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
j_m^{[0]}(X)\big|_{X=\frac{1}{2}} &= \frac{1}{60} \\
\frac{\partial^1 (j_m^{[0]}(X))}{\partial X^1} &= \frac{2 X}{15} \\
\left.\frac{\partial^1 (j_m^{[0]}(X))}{\partial X^1}\right|_{X=\frac{1}{2}} &= \frac{1}{15} \\
\frac{\partial^2 (j_m^{[0]}(X))}{\partial X^2} &= \frac{2}{15} \\
\left.\frac{\partial^2 (j_m^{[0]}(X))}{\partial X^2}\right|_{X=\frac{1}{2}} &= \frac{2}{15}
\end{aligned}
$$

$X=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
j_m^{[1]}(X)\big|_{X=\frac{1}{2}} &= \frac{11}{672} \\
\frac{\partial^1 (j_m^{[1]}(X))}{\partial X^1} &= \frac{2 X \left(7 - X^{2}\right)}{105} \\
\left.\frac{\partial^1 (j_m^{[1]}(X))}{\partial X^1}\right|_{X=\frac{1}{2}} &= \frac{9}{140} \\
\frac{\partial^2 (j_m^{[1]}(X))}{\partial X^2} &= \frac{2}{15} - \frac{2 X^{2}}{35} \\
\left.\frac{\partial^2 (j_m^{[1]}(X))}{\partial X^2}\right|_{X=\frac{1}{2}} &= \frac{5}{42}
\end{aligned}
$$

$X=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
(j_m^{[1]}-j_m^{[0]})/j_m^{[0]}\big|_{X=\frac{1}{2}} &= - \frac{1}{56} \\
\frac{\partial^1 ((j_m^{[1]}-j_m^{[0]})/j_m^{[0]})}{\partial X^1} &= - \frac{X}{7} \\
\left.\frac{\partial^1 ((j_m^{[1]}-j_m^{[0]})/j_m^{[0]})}{\partial X^1}\right|_{X=\frac{1}{2}} &= - \frac{1}{14} \\
\frac{\partial^2 ((j_m^{[1]}-j_m^{[0]})/j_m^{[0]})}{\partial X^2} &= - \frac{1}{7} \\
\left.\frac{\partial^2 ((j_m^{[1]}-j_m^{[0]})/j_m^{[0]})}{\partial X^2}\right|_{X=\frac{1}{2}} &= - \frac{1}{7}
\end{aligned}
$$
:::

:::example[球 Bessel 関数の原点付近（3）]
$x^\ell$ を除いたあとの補正を調べる。さらに次数を上げて、低次の結果と照合する。

$$
m=3
$$

最初の項を求める。

$$
\begin{aligned}
j_m^{[0]}(X)\big|_{m=3} &= \frac{X^{3}}{105}
\end{aligned}
$$

二乗の補正を加える。

$$
\begin{aligned}
j_m^{[1]}(X)\big|_{m=3} &= \frac{X^{3} \left(18 - X^{2}\right)}{1890}
\end{aligned}
$$

先頭に対する補正の比を取る。

$$
\begin{aligned}
(j_m^{[1]}-j_m^{[0]})/j_m^{[0]}\big|_{m=3} &= - \frac{X^{2}}{18}
\end{aligned}
$$

$X=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
j_m^{[0]}(X)\big|_{X=\frac{1}{2}} &= \frac{1}{840} \\
\frac{\partial^1 (j_m^{[0]}(X))}{\partial X^1} &= \frac{X^{2}}{35} \\
\left.\frac{\partial^1 (j_m^{[0]}(X))}{\partial X^1}\right|_{X=\frac{1}{2}} &= \frac{1}{140} \\
\frac{\partial^2 (j_m^{[0]}(X))}{\partial X^2} &= \frac{2 X}{35} \\
\left.\frac{\partial^2 (j_m^{[0]}(X))}{\partial X^2}\right|_{X=\frac{1}{2}} &= \frac{1}{35}
\end{aligned}
$$

$X=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
j_m^{[1]}(X)\big|_{X=\frac{1}{2}} &= \frac{71}{60480} \\
\frac{\partial^1 (j_m^{[1]}(X))}{\partial X^1} &= - \frac{X^{4}}{378} + \frac{X^{2}}{35} \\
\left.\frac{\partial^1 (j_m^{[1]}(X))}{\partial X^1}\right|_{X=\frac{1}{2}} &= \frac{211}{30240} \\
\frac{\partial^2 (j_m^{[1]}(X))}{\partial X^2} &= \frac{2 X \left(27 - 5 X^{2}\right)}{945} \\
\left.\frac{\partial^2 (j_m^{[1]}(X))}{\partial X^2}\right|_{X=\frac{1}{2}} &= \frac{103}{3780}
\end{aligned}
$$

$X=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
(j_m^{[1]}-j_m^{[0]})/j_m^{[0]}\big|_{X=\frac{1}{2}} &= - \frac{1}{72} \\
\frac{\partial^1 ((j_m^{[1]}-j_m^{[0]})/j_m^{[0]})}{\partial X^1} &= - \frac{X}{9} \\
\left.\frac{\partial^1 ((j_m^{[1]}-j_m^{[0]})/j_m^{[0]})}{\partial X^1}\right|_{X=\frac{1}{2}} &= - \frac{1}{18} \\
\frac{\partial^2 ((j_m^{[1]}-j_m^{[0]})/j_m^{[0]})}{\partial X^2} &= - \frac{1}{9} \\
\left.\frac{\partial^2 ((j_m^{[1]}-j_m^{[0]})/j_m^{[0]})}{\partial X^2}\right|_{X=\frac{1}{2}} &= - \frac{1}{9}
\end{aligned}
$$
:::

:::example[球 Bessel 関数の原点付近（4）]
$x^\ell$ を除いたあとの補正を調べる。四次の結果まで具体化する。

$$
m=4
$$

最初の項を求める。

$$
\begin{aligned}
j_m^{[0]}(X)\big|_{m=4} &= \frac{X^{4}}{945}
\end{aligned}
$$

二乗の補正を加える。

$$
\begin{aligned}
j_m^{[1]}(X)\big|_{m=4} &= \frac{X^{4} \left(22 - X^{2}\right)}{20790}
\end{aligned}
$$

先頭に対する補正の比を取る。

$$
\begin{aligned}
(j_m^{[1]}-j_m^{[0]})/j_m^{[0]}\big|_{m=4} &= - \frac{X^{2}}{22}
\end{aligned}
$$

$X=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
j_m^{[0]}(X)\big|_{X=\frac{1}{2}} &= \frac{1}{15120} \\
\frac{\partial^1 (j_m^{[0]}(X))}{\partial X^1} &= \frac{4 X^{3}}{945} \\
\left.\frac{\partial^1 (j_m^{[0]}(X))}{\partial X^1}\right|_{X=\frac{1}{2}} &= \frac{1}{1890} \\
\frac{\partial^2 (j_m^{[0]}(X))}{\partial X^2} &= \frac{4 X^{2}}{315} \\
\left.\frac{\partial^2 (j_m^{[0]}(X))}{\partial X^2}\right|_{X=\frac{1}{2}} &= \frac{1}{315}
\end{aligned}
$$

$X=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
j_m^{[1]}(X)\big|_{X=\frac{1}{2}} &= \frac{29}{443520} \\
\frac{\partial^1 (j_m^{[1]}(X))}{\partial X^1} &= \frac{X^{3} \left(44 - 3 X^{2}\right)}{10395} \\
\left.\frac{\partial^1 (j_m^{[1]}(X))}{\partial X^1}\right|_{X=\frac{1}{2}} &= \frac{173}{332640} \\
\frac{\partial^2 (j_m^{[1]}(X))}{\partial X^2} &= \frac{X^{2} \left(44 - 5 X^{2}\right)}{3465} \\
\left.\frac{\partial^2 (j_m^{[1]}(X))}{\partial X^2}\right|_{X=\frac{1}{2}} &= \frac{19}{6160}
\end{aligned}
$$

$X=\frac{1}{2}$ で値と導関数を評価する。

$$
\begin{aligned}
(j_m^{[1]}-j_m^{[0]})/j_m^{[0]}\big|_{X=\frac{1}{2}} &= - \frac{1}{88} \\
\frac{\partial^1 ((j_m^{[1]}-j_m^{[0]})/j_m^{[0]})}{\partial X^1} &= - \frac{X}{11} \\
\left.\frac{\partial^1 ((j_m^{[1]}-j_m^{[0]})/j_m^{[0]})}{\partial X^1}\right|_{X=\frac{1}{2}} &= - \frac{1}{22} \\
\frac{\partial^2 ((j_m^{[1]}-j_m^{[0]})/j_m^{[0]})}{\partial X^2} &= - \frac{1}{11} \\
\left.\frac{\partial^2 ((j_m^{[1]}-j_m^{[0]})/j_m^{[0]})}{\partial X^2}\right|_{X=\frac{1}{2}} &= - \frac{1}{11}
\end{aligned}
$$
:::

## 漸近形と零点

### 遠方の位相と零点の間隔

整数 $m$ の積分表示では $\sin\theta$ の停留点 $\pm\pi/2$ を使う。二点の位相は $\pm(x-m\pi/2-\pi/4)$、振幅はそれぞれ $1/\sqrt{2\pi x}$ なので
$$
J_m(x)=\sqrt{\frac2{\pi x}}\cos(x-m\pi/2-\pi/4)+O(x^{-1}).
$$
固定実次数 $\nu$ の標準的な遠方規格化では $J_\nu,Y_\nu$ の主項はそれぞれ $\sqrt{2/(\pi x)}\cos\chi,\sqrt{2/(\pi x)}\sin\chi$ で、$\chi=x-\nu\pi/2-\pi/4$ である。一般次数の位相と誤差はこの節の後半で導く。上式により Hankel の主項は $\sqrt{2/(\pi x)}e^{\pm i\chi}$ である。原点から数えた正零点を $j_{\nu,n}$ と書く。正零点が単純なのは、値と一階微分が両方零なら正則点での初期値問題の一意性に反するためである。例えば $j_0$ の正零点は厳密に $n\pi$ であり、$J_{1/2}$ も同じ零点を持つ。

### 二つの端点から遠方位相を決める

式 :eqref[eq-mp-bessel-poisson] で $\nu>-1/2$ とし、区間を高さ $\eta>0$ の水平辺と左右の縦辺へ移す。重みの枝は実区間で正、内部では $\log(1-t^2)$ をその枝から延長する。端点を半径 $\epsilon$ の弧で除くと弧の積分は $O(\epsilon^{\nu+1/2})$ で消える。上辺は $O(e^{-\eta x})$ である。右端の縦辺では元の積分への寄与は $-i e^{ix}\int_0^\eta e^{-xv}(v^2-2iv)^{\nu-1/2}\dd v$ である。$v^2-2iv=(-2i)v(1+iv/2)$ と分け、:ref[thm-mp-watson]{text="Watson の補題"} を使うと積分の主項は $(-2i)^{\nu-1/2}\Gamma(\nu+1/2)x^{-\nu-1/2}$ である。従って $J_\nu$ への寄与は $e^{i(x-\nu\pi/2-\pi/4)}/\sqrt{2\pi x}$ になる。左端は実次数ならその共役であり、合わせて余弦項を与える。次の項との比は各端で $O(x^{-1})$ である。

任意の固定次数への延長は $J_{\nu-1}=2\nu J_\nu/x-J_{\nu+1}$ で行える。主項では右辺の第 1 項が $O(x^{-3/2})$ なので、第 2 項の符号反転だけで位相が一つ下の次数へ移る。非整数の $Y_\nu$ には $J_\nu,J_{-\nu}$ の余弦を定義に戻すと正弦項が残る。整数での極限に必要な誤差の制御は、積分の Taylor 剰余が $\Re\nu>-1/2$ のコンパクト集合で一様であることによる。各端の変数を $v=u/x$ とすれば $x^\nu x^{-\nu-1/2}=x^{-1/2}$ で、剰余はその $1/x$ 倍で抑えられる。漸化式で動かした後も誤差は同じく局所一様である。$Y_\nu$ の分子では整数で主項も誤差も消えるので、Cauchy の微分評価により $\sin\pi\nu$ で割った極限も $O(x^{-3/2})$ である。

### 正零点を原点から数える

$\nu\ge0$ の固定した値と $1/2$ を結ぶコンパクト区間を取る。そこでは $J_\nu(x)=(x/2)^\nu[1+O(x^2)]/\Gamma(\nu+1)$ が $x\downarrow0$ で一様で、十分小さい共通の $x$ の範囲で正である。正零点が単純であるため、その位置は次数を少しずらしても陰関数定理で連続に動く。固定した十分大きな整数 $m$ に対し $b_m(\nu)=(m+\nu/2+1/4)\pi$ では $\chi=m\pi$ である。遠方展開の一様な誤差 $O(x^{-3/2})$ により $J_\nu(b_m(\nu))\ne0$ である。よって $(0,b_m(\nu))$ の零点数は次数の区間で一定である。実際、原点付近も動く端点付近も零点を通さず、内部の各単純根は局所的に一本の枝をなすためである。$\nu=1/2$ では $b_m=(m+1/2)\pi$、正根は $n\pi$ なので数は $m$ 個である。

この議論を $m$ と $m-1$ に使うと、十分大きい $m$ の第 $m$ 正零点は $b_{m-1}(\nu)$ と $b_m(\nu)$ の間にある。その点では主項の余弦の絶対値が $O(1/x)$ でなければならない。二つの境界の間で余弦の零はただ一つで、そこでは導関数の絶対値が $1$ なので
$$
j_{\nu,m}=(m+\nu/2-1/4)\pi+O(m^{-1})
$$
を得る。遠方展開を次数の一点だけで見ると根の間隔は分かるが、原点からの個数までは決まらない。ここでは既知の半整数の根と次数の連続性が添字のずれを固定している。

### 隣接する次数の零点の交錯

$\nu\ge0$ とする。二つの漸化式から $(x^{-\nu}J_\nu)'=-x^{-\nu}J_{\nu+1}$ である。$F=x^{-\nu}J_\nu$ と書くと、Bessel 方程式は $(x^{2\nu+1}F')'=-x^{2\nu+1}F$ になる。$F$ の隣接零点の間では右辺の符号が一定かつ非零なので、$x^{2\nu+1}F'$ は狭義単調である。Rolle の定理により $F'$ は少なくとも一つ零になり、狭義単調性によりちょうど一つ零になる。原点から $F$ の最初の正零点まででは $F>0$ であり、$x^{2\nu+1}F'\to0$ だから $F'<0$ である。従ってこの区間に $J_{\nu+1}$ の正零点はない。さらに $J_\nu$ の正零点では単純性により $J_{\nu+1}=-J_\nu'\ne0$ である。以上を合わせると
$$
 j_{\nu,n}<j_{\nu+1,n}<j_{\nu,n+1}\qquad(n\ge1)
$$
となる。例えば $J_{1/2}$ の零点 $n\pi$ を使うと、$j_1(x)$ の各正零点が $(n\pi,(n+1)\pi)$ に一つずつあると分かる。

### 遠方での補正

解を $X^{-1/2}e^{iX}v(X)$ とすると $v^{\prime\prime}+2iv^\prime+(1/4-\nu^2)v/X^2=0$ になる。$v=\sum c_j X^{-j}$ と置くと、$c_{j+1}=[j(j+1)+1/4-\nu^2]c_j/[2i(j+1)]$ である。これは形式的な漸化式である。$J_\nu$ の実軸上の展開については、上の縦辺積分で $(1+iv/2)^{\nu-1/2}$ を任意の有限次数で打ち切り、Watson の補題の剰余評価を適用すれば同じ係数が得られる。端点の Gamma 因子の比は $\Gamma(\nu+1/2+j)/\Gamma(\nu+1/2)=(\nu+1/2)_j$ である。従って $c_0=1$ としたとき、第 1 補正は $c_1=i(4\nu^2-1)/8$ になる。右端と左端を足す際に $i e^{i\chi}$ の実部が $-\sin\chi$ になるので
$$
J_\nu(x)=\sqrt{\frac2{\pi x}}\left(\cos\chi-\frac{4\nu^2-1}{8x}\sin\chi+O(x^{-2})\right)
$$
である。任意の固定実次数への延長には前と同じ漸化式を使える。例えば $\nu=1/2$ では補正係数が零で、既知の半整数の厳密解と一致する。

:::example[遠方での形式係数（1）]
最初の二つの補正の位相を比べる。

$$
m=1
$$

最初の係数は純虚数である。

$$
\begin{aligned}
c_1\big|_{m=1} &= \frac{3 i}{8}
\end{aligned}
$$

次の係数は実数である。

$$
\begin{aligned}
c_2\big|_{m=1} &= \frac{15}{128}
\end{aligned}
$$

定数の規格化を $1$ とした有限和を作る。

$$
\begin{aligned}
v^{[2]}(X)\big|_{m=1} &= 1 + \frac{3 i}{8 X} + \frac{15}{128 X^{2}}
\end{aligned}
$$

$X=20$ で値と導関数を評価する。

$$
\begin{aligned}
v^{[2]}(X)\big|_{X=20} &= \frac{10243}{10240} + \frac{3 i}{160} \\
\frac{\partial^1 (v^{[2]}(X))}{\partial X^1} &= \frac{3 \left(- 8 i X - 5\right)}{64 X^{3}} \\
\left.\frac{\partial^1 (v^{[2]}(X))}{\partial X^1}\right|_{X=20} &= - \frac{3}{102400} - \frac{3 i}{3200}
\end{aligned}
$$
:::

:::example[遠方での形式係数（2）]
最初の二つの補正の位相を比べる。次数を一つ上げる。

$$
m=2
$$

最初の係数は純虚数である。

$$
\begin{aligned}
c_1\big|_{m=2} &= \frac{15 i}{8}
\end{aligned}
$$

次の係数は実数である。

$$
\begin{aligned}
c_2\big|_{m=2} &= - \frac{105}{128}
\end{aligned}
$$

定数の規格化を $1$ とした有限和を作る。

$$
\begin{aligned}
v^{[2]}(X)\big|_{m=2} &= 1 + \frac{15 i}{8 X} - \frac{105}{128 X^{2}}
\end{aligned}
$$

$X=20$ で値と導関数を評価する。

$$
\begin{aligned}
v^{[2]}(X)\big|_{X=20} &= \frac{10219}{10240} + \frac{3 i}{32} \\
\frac{\partial^1 (v^{[2]}(X))}{\partial X^1} &= \frac{15 \left(- 8 i X + 7\right)}{64 X^{3}} \\
\left.\frac{\partial^1 (v^{[2]}(X))}{\partial X^1}\right|_{X=20} &= \frac{21}{102400} - \frac{3 i}{640}
\end{aligned}
$$
:::

:::example[遠方での形式係数（3）]
最初の二つの補正の位相を比べる。さらに次数を上げて、低次の結果と照合する。

$$
m=3
$$

最初の係数は純虚数である。

$$
\begin{aligned}
c_1\big|_{m=3} &= \frac{35 i}{8}
\end{aligned}
$$

次の係数は実数である。

$$
\begin{aligned}
c_2\big|_{m=3} &= - \frac{945}{128}
\end{aligned}
$$

定数の規格化を $1$ とした有限和を作る。

$$
\begin{aligned}
v^{[2]}(X)\big|_{m=3} &= 1 + \frac{35 i}{8 X} - \frac{945}{128 X^{2}}
\end{aligned}
$$

$X=20$ で値と導関数を評価する。

$$
\begin{aligned}
v^{[2]}(X)\big|_{X=20} &= \frac{10051}{10240} + \frac{7 i}{32} \\
\frac{\partial^1 (v^{[2]}(X))}{\partial X^1} &= \frac{35 \left(- 8 i X + 27\right)}{64 X^{3}} \\
\left.\frac{\partial^1 (v^{[2]}(X))}{\partial X^1}\right|_{X=20} &= \frac{189}{102400} - \frac{7 i}{640}
\end{aligned}
$$
:::

:::example[遠方での形式係数（4）]
最初の二つの補正の位相を比べる。四次の結果まで具体化する。

$$
m=4
$$

最初の係数は純虚数である。

$$
\begin{aligned}
c_1\big|_{m=4} &= \frac{63 i}{8}
\end{aligned}
$$

次の係数は実数である。

$$
\begin{aligned}
c_2\big|_{m=4} &= - \frac{3465}{128}
\end{aligned}
$$

定数の規格化を $1$ とした有限和を作る。

$$
\begin{aligned}
v^{[2]}(X)\big|_{m=4} &= 1 + \frac{63 i}{8 X} - \frac{3465}{128 X^{2}}
\end{aligned}
$$

$X=20$ で値と導関数を評価する。

$$
\begin{aligned}
v^{[2]}(X)\big|_{X=20} &= \frac{9547}{10240} + \frac{63 i}{160} \\
\frac{\partial^1 (v^{[2]}(X))}{\partial X^1} &= \frac{63 \left(- 8 i X + 55\right)}{64 X^{3}} \\
\left.\frac{\partial^1 (v^{[2]}(X))}{\partial X^1}\right|_{X=20} &= \frac{693}{102400} - \frac{63 i}{3200}
\end{aligned}
$$
:::

## Fourier–Bessel 展開と Hankel 変換

### 円板の離散化と全平面の連続化

$\nu\ge0$ とし、$k_n=j_{\nu,n}/a$、$R_n(r)=J_\nu(k_n r)$ を $0<r<a$ に取る。原点で有界な側を選んでおり、$r=a$ で零である。$\nu$ が非整数なら原点での $r^\nu$ は一般に解析的ではない。方程式を $(rR_n')'+(k_n^2r-\nu^2/r)R_n=0$ と書き、異なる二解を交差して引き算すると、境界項が消えて重み $r$ で直交する。ノルムを得るには $u(r,k)=J_\nu(kr)$ の方程式を $k$ で微分し、元の式を掛けて引く。積分すると $2k\int r u^2\dd r=[r(\partial_k u)u'-ru\partial_k u']_0^a$ である。$u(a,k_n)=0$ で $u'(a,k_n)=k_nJ_\nu'(j_{\nu,n})$、$\partial_k u(a,k_n)=aJ_\nu'(j_{\nu,n})$ なのでノルムは $a^2J_{\nu+1}(j_{\nu,n})^2/2$ となる。従って展開係数は $\int_0^a r f(r)R_n(r)\dd r$ をこの値で割る。完全性は以下で逆作用素を構成して示す。[特異 Sturm–Liouville 問題](/physics/mathematical-physics/sturm-liouville/) では、この原点条件を端点の分類と関連付ける。

角依存 $f(r)e^{im\theta}$ の二次元 Fourier 変換は式 :eqref[eq-mp-bessel-generating] により $2\pi(-i)^m e^{im\phi}\int_0^\infty f(r)J_m(kr)r\dd r$ である。角積分で位相 $(-i)^m$ が出て、逆変換の角積分が逆の位相を打ち消す。従って整数次数 $m\ge0$ の :term[ハンケル変換]{en="Hankel transform" yomi="はんけるへんかん"} の逆核も $J_m(kr)$ で、逆の測度は $k\dd k$ である。Plancherel は重み $r\dd r$ と $k\dd k$ の間の等長性を与える。例えば円板の箱 $f=\mathbf1_{r<a}$ の零次変換は上の動径積分 $a J_1(ka)/k$ で、$k=0$ で $a^2/2$ になる。

### 有限半径での完全性

$Lu=r^{-1}[-(ru')'+(\nu^2/r)u]$ とおく。$Lu=0$ の原点で有界な解は $r^\nu$（$\nu=0$ では $1$）である。二つの零エネルギー解を用いて逆作用素の核を
$$
G(r,s)=(rs)^\nu\int_{\max(r,s)}^a t^{-1-2\nu}\dd t
$$
と取る。各 $s$ を固定すると左右で同次方程式を満たし、$r=a$ で零である。$s$ で連続であり、$r\partial_rG$ は右へ $-1$ だけ跳ぶ。従って $Kf(r)=\int_0^a G(r,s)f(s)s\dd s$ は $L Kf=f$ を満たし、$K$ は単射である。核は実対称で、$\nu>0$ では $0\le G\le1/(2\nu)$、$\nu=0$ では $G=\log(a/\max(r,s))$ である。後者の二乗も $\int_0^a s^3\log^2(a/s)\dd s<\infty$ により可積分である。従って $K$ は自己共役であり、核を有限和の積関数で $L^2$ 近似すれば有限階数近似によりコンパクトである。

Fubini の定理で積分の順序を変えると
$$
\inner{Kf}{f}_{r\dd r}=\int_0^a t^{-1-2\nu}
  \left|\int_0^t f(s)s^{\nu+1}\dd s\right|^2\dd t\ge0 .
$$
この交換は $|f|$ を入れた二重積分が Cauchy–Schwarz と核の二乗可積分性で有限になることによる。また $|Kf(r)|\le\|G(r,\cdot)\|_{L^2(s\dd s)}\|f\|_{L^2(s\dd s)}$ で、右辺は $\nu>0$ の核の上界または $G(r,s)\le\log(a/s)$ により $r$ に一様である。従って $K$ の非零固有値 $\mu$ に属する固有ベクトルは $Lu=\mu^{-1}u$ の原点で有界な解である。$\mu>0$ なのでこれは $J_\nu(r/\sqrt\mu)$ の倍数であり、$u(a)=0$ が正零点を選ぶ。逆にこれらの Bessel モードは、$KLu-u$ が原点で有界かつ外端で零の同次解なので、$K$ の固有ベクトルである。:ref[thm-fa-compact-spectral]{full} と $K$ の単射性から、これらのモードは $L^2((0,a),r\dd r)$ に完全である。

### 実次数の Hankel 変換と Laguerre 基底

$\nu\ge0$ とする。関数 $\Phi_n(r)=r^\nu e^{-r^2/2}L_n^{(\nu)}(r^2)$ を取る。$x=r^2$ への変換により、内積は重み $x^\nu e^{-x}$ の Laguerre 内積の半分になる。従って [Laguerre 多項式の完全性](/physics/mathematical-physics/special-functions/) から、$\Phi_n$ を規格化した列は $L^2((0,\infty),r\dd r)$ の正規直交基底である。

実部が正の $c$ について、式 :eqref[eq-mp-bessel-series] を Gaussian と掛けて項別積分すると
$$
\int_0^\infty r^{\nu+1}e^{-cr^2}J_\nu(kr)\dd r
  =\frac{k^\nu}{(2c)^{\nu+1}}e^{-k^2/(4c)}
$$
となる。実際、$j$ 番目の積分は $\Gamma(\nu+j+1)/(2c^{\nu+j+1})$ を与え、絶対値の和も $\Re c$ を用いた同じ指数級数で有限になる。Laguerre の母関数を掛けると $c=(1+t)/(2(1-t))$ になり、右辺の $k$ 側に残る母関数は $t$ を $-t$ に替えたものである。$|t|$ を固定して $1$ より小さく取れば Gaussian の減衰が一様なので、係数比較により
$$
\int_0^\infty \Phi_n(r)J_\nu(kr)r\dd r=(-1)^n\Phi_n(k)
$$
を得る。従って積分変換は基底上で符号を掛ける作用であり、全空間へ等長に延長できる。二度作用させると全成分の符号が戻るため、逆変換の核も $J_\nu$ である。例えば $r^\nu e^{-r^2/2}=\Phi_0(r)$ は変換しても同じ形を保つ。$f\in C_c^\infty(0,\infty)$ では実際の積分はこの等長拡張と一致する。原点側で積分は $O(k^\nu)$ である。遠方では Bessel 方程式と部分積分で $k^2$ を積分から取り出して $\int (Lf)(r)J_\nu(kr)r\dd r$ に移し、$J_\nu$ の漸近形から $O(k^{-5/2})$ を得る。従って積分は $L^2(k\dd k)$ に属する。その $\Phi_n$ 成分は Fubini の定理により $(-1)^n\inner{f}{\Phi_n}_{r\dd r}$ であり、基底上で定義した拡張の成分と同じである。ここでの核 $J_\nu(kr)$ と測度 $r\dd r$ は、$f(r)$ と変換の結果にそれぞれ $\sqrt r$ と $\sqrt k$ を掛けると、核 $\sqrt{kr}J_\nu(kr)$ と測度 $\dd r$ の形になる。一般の $L^2$ データに対しては、この稠密集合からの $L^2$ 極限で読む。

### 球面のモードごとの Hankel 変換

平面波 $e^{-ikr\cos\gamma}$ の角部分を Legendre 多項式に展開する。式 :eqref[eq-mp-legendre-rodrigues] を $-\!1<x<1$ で用いて $\ell$ 回部分積分すると、
$$
\int_{-1}^1 e^{-ikrx}P_\ell(x)\dd x
=\frac{(-ikr)^\ell}{2^\ell\ell!}\int_{-1}^1 e^{-ikrx}(1-x^2)^\ell\dd x
=2(-i)^\ell j_\ell(kr)
$$
となる。境界項は $(x^2-1)^\ell$ の零で消え、最後は式 :eqref[eq-mp-bessel-poisson] に $\nu=\ell+1/2$ を入れた。加法定理で角積分すると、$f(r)Y_{\ell m}(\Omega)$ の三次元 Fourier 変換は
$$
4\pi(-i)^\ell Y_{\ell m}(\hat k)\int_0^\infty r^2 f(r)j_\ell(kr)\dd r
$$
である。逆では $(2\pi)^{-3}(4\pi)^2=2/\pi$ が前に残り、$k^2\dd k$ の積分となる。例えば $\ell=0$ なら $j_0(kr)=\sin(kr)/(kr)$ なので、$r f(r)$ に対する正弦変換へ帰着する。円筒の場合の $r\dd r$ を球にそのまま用いるとこの規格化は得られない。

### 動径モードの測度

動径の式を $(rR^\prime)^\prime+(k^2r-m^2/r)R=0$ と書くと、固有値の係数が $r$ である。二解に相手の式を掛けて引くと $\int R_jR_\ell\,r\dd r$ が出る。表面の Dirichlet 境界と原点で境界項が消えるため、二つの波数が異なれば直交する。球の場合は三次元の体積測度により重みが $r^2$ に変わる。

:::example[球 Bessel の規格化（1）]
半径 $2$ の球の軸対称の次数零には $j_0(kr)=\sin(kr)/(kr)$ を使う。

$$
m=1
$$

境界で正弦が零になるようにする。

$$
\begin{aligned}
k_m\big|_{m=1} &= \frac{\pi}{2}
\end{aligned}
$$

球の動径重みを掛けると正弦の二乗になる。

$$
\begin{aligned}
\int_0^2 r^2j_0(k_mr)^2\dd r\big|_{m=1} &= \frac{4}{\pi^{2}}
\end{aligned}
$$

動径関数をノルムで割る。

$$
\begin{aligned}
\widetilde R_m(r)\big|_{m=1} &= \frac{\sin{\left(\frac{\pi r}{2} \right)}}{r}
\end{aligned}
$$

$r=1$ で値と導関数を評価する。

$$
\begin{aligned}
\widetilde R_m(r)\big|_{r=1} &= 1 \\
\frac{\partial^1 (\widetilde R_m(r))}{\partial r^1} &= \frac{\frac{\pi r \cos{\left(\frac{\pi r}{2} \right)}}{2} - \sin{\left(\frac{\pi r}{2} \right)}}{r^{2}} \\
\left.\frac{\partial^1 (\widetilde R_m(r))}{\partial r^1}\right|_{r=1} &= -1
\end{aligned}
$$
:::

:::example[球 Bessel の規格化（2）]
半径 $2$ の球の軸対称の次数零には $j_0(kr)=\sin(kr)/(kr)$ を使う。次数を一つ上げる。

$$
m=2
$$

境界で正弦が零になるようにする。

$$
\begin{aligned}
k_m\big|_{m=2} &= \pi
\end{aligned}
$$

球の動径重みを掛けると正弦の二乗になる。

$$
\begin{aligned}
\int_0^2 r^2j_0(k_mr)^2\dd r\big|_{m=2} &= \frac{1}{\pi^{2}}
\end{aligned}
$$

動径関数をノルムで割る。

$$
\begin{aligned}
\widetilde R_m(r)\big|_{m=2} &= \frac{\sin{\left(\pi r \right)}}{r}
\end{aligned}
$$

$r=1$ で値と導関数を評価する。

$$
\begin{aligned}
\widetilde R_m(r)\big|_{r=1} &= 0 \\
\frac{\partial^1 (\widetilde R_m(r))}{\partial r^1} &= \frac{\pi r \cos{\left(\pi r \right)} - \sin{\left(\pi r \right)}}{r^{2}} \\
\left.\frac{\partial^1 (\widetilde R_m(r))}{\partial r^1}\right|_{r=1} &= - \pi
\end{aligned}
$$
:::

:::example[球 Bessel の規格化（3）]
半径 $2$ の球の軸対称の次数零には $j_0(kr)=\sin(kr)/(kr)$ を使う。さらに次数を上げて、低次の結果と照合する。

$$
m=3
$$

境界で正弦が零になるようにする。

$$
\begin{aligned}
k_m\big|_{m=3} &= \frac{3 \pi}{2}
\end{aligned}
$$

球の動径重みを掛けると正弦の二乗になる。

$$
\begin{aligned}
\int_0^2 r^2j_0(k_mr)^2\dd r\big|_{m=3} &= \frac{4}{9 \pi^{2}}
\end{aligned}
$$

動径関数をノルムで割る。

$$
\begin{aligned}
\widetilde R_m(r)\big|_{m=3} &= \frac{\sin{\left(\frac{3 \pi r}{2} \right)}}{r}
\end{aligned}
$$

$r=1$ で値と導関数を評価する。

$$
\begin{aligned}
\widetilde R_m(r)\big|_{r=1} &= -1 \\
\frac{\partial^1 (\widetilde R_m(r))}{\partial r^1} &= \frac{\frac{3 \pi r \cos{\left(\frac{3 \pi r}{2} \right)}}{2} - \sin{\left(\frac{3 \pi r}{2} \right)}}{r^{2}} \\
\left.\frac{\partial^1 (\widetilde R_m(r))}{\partial r^1}\right|_{r=1} &= 1
\end{aligned}
$$
:::

:::example[球 Bessel の規格化（4）]
半径 $2$ の球の軸対称の次数零には $j_0(kr)=\sin(kr)/(kr)$ を使う。四次の結果まで具体化する。

$$
m=4
$$

境界で正弦が零になるようにする。

$$
\begin{aligned}
k_m\big|_{m=4} &= 2 \pi
\end{aligned}
$$

球の動径重みを掛けると正弦の二乗になる。

$$
\begin{aligned}
\int_0^2 r^2j_0(k_mr)^2\dd r\big|_{m=4} &= \frac{1}{4 \pi^{2}}
\end{aligned}
$$

動径関数をノルムで割る。

$$
\begin{aligned}
\widetilde R_m(r)\big|_{m=4} &= \frac{\sin{\left(2 \pi r \right)}}{r}
\end{aligned}
$$

$r=1$ で値と導関数を評価する。

$$
\begin{aligned}
\widetilde R_m(r)\big|_{r=1} &= 0 \\
\frac{\partial^1 (\widetilde R_m(r))}{\partial r^1} &= \frac{2 \pi r \cos{\left(2 \pi r \right)} - \sin{\left(2 \pi r \right)}}{r^{2}} \\
\left.\frac{\partial^1 (\widetilde R_m(r))}{\partial r^1}\right|_{r=1} &= 2 \pi
\end{aligned}
$$
:::

## 演習問題

::::exercise{level=basic}
$J_0$ の $x^4$ 係数を求めよ。
:::solution
式 :eqref[eq-mp-bessel-series] に $j=2,\nu=0$ を代入して $1/(2^4(2!)^2)=1/64$ である。
:::
::::

::::exercise{level=basic}
$J_{-2}$ と $J_2$ の関係を示せ。
:::solution
$j<2$ は逆 Gamma の零で消える。$j=n+2$ とすると符号は $(-1)^{n+2}=(-1)^n$、分母は $(n+2)!n!$ で $J_2$ の級数になる。
:::
::::

::::exercise{level=basic}
$j_1(x)$ の小 $x$ の主項を求めよ。
:::solution
$\sin x/x^2=x^{-1}-x/6+O(x^3)$、$\cos x/x=x^{-1}-x/2+O(x^3)$ なので差は $x/3+O(x^3)$ である。
:::
::::

::::exercise{level=basic}
$H_0^{(1)}$ と $H_0^{(2)}$ の和を求めよ。
:::solution
定義で虚部が反対符号なので和は $2J_0$ である。実 $x$ では二つは複素共役である。
:::
::::

::::exercise{level=advanced}
$J_0$ が最初の正零点 $j$ で減少して横切ることを示せ。
:::solution
$J_0(0)=1$ と零点の定義から直前は正なので、$J_0'(j)\le0$ である。単純性で等号は除かれ、$J_0'(j)<0$ となる。
:::
::::

::::exercise{level=advanced}
$\nu=1/2$ で $K_{1/2}$ を指数で表せ。
:::solution
半整数の級数の符号を正に替えると $I_{\pm1/2}$ は共通係数 $\sqrt{2/(\pi x)}$ に、正次数で $\sinh x$、負次数で $\cosh x$ を掛けたものになる。$\sin(\pi/2)=1$ なので $K_{1/2}=\sqrt{\pi/(2x)}e^{-x}$ である。
:::
::::

::::exercise{level=advanced}
半径 $a$ の円板で定数 $1$ の零次 Fourier–Bessel 係数を求めよ。
:::solution
分子は $\int_0^a rJ_0(j_{0,n}r/a)\dd r=a^2J_1(j_{0,n})/j_{0,n}$ である。ノルムで割って $2/[j_{0,n}J_1(j_{0,n})]$ となる。境界の一点は重み付き $L^2$ の展開に影響しない。
:::
::::

::::exercise{level=advanced}
$j_0$ の方程式の独立解として $-\cos x/x$ を確認せよ。
:::solution
$y=v/x$ と置けば $y''+2y'/x+y=(v''+v)/x$ となる。$v=\sin x$ と $v=-\cos x$ は $v^{\prime\prime}+v=0$ を満たし、この順で Wronskian は $1$ なので独立であるから、元の解も独立である。
:::
::::

## 確認テスト

::::quiz{answer="1" tol="1e-9"}
$J_0(0)$ の連続延長の値を求めよ。

:::explanation
定数項は $1$ で、$j\ge1$ の全項は $x\to0$ で消える。
:::
::::

::::quiz
時間因子 $e^{-i\omega t}$ に対し外向きの円筒波を選べ。

- [x] $H_\nu^{(1)}(kr)$
- [ ] $H_\nu^{(2)}(kr)$
- [ ] $K_\nu(kr)$

:::explanation
$H^{(1)}$ の位相は遠方で $kr-\omega t$ なので等位相点が時間とともに外へ動く。第 2 種では $-kr-\omega t$、$K_\nu$ は指数減衰である。
:::
::::

::::quiz{answer="2" tol="1e-9"}
$j_0$ の第 2 正零点を $\pi$ の倍数で答えよ。

:::explanation
正零点は $x>0$ で $\sin x=0$ となる $n\pi$ で、第 2 は $2\pi$ である。
:::
::::

::::quiz
円板の動径直交性で使う測度はどれか。

- [ ] $\dd r$
- [x] $r\dd r$
- [ ] $r^2\dd r$

:::explanation
極座標の面積は $r\dd r\dd\phi$ である。固有方程式を発散形で書いたときにも固有値に $r$ が掛かる。
:::
::::

---
title: "積分変換"
description: "Fourier・Laplace・Mellin・Hankel 変換を導き、周期信号、初期値問題、標本化を扱う。"
notation:
  - symbol: '\widehat f(k)'
    meaning: "負の指数を用いた Fourier 変換"
  - symbol: '\mathcal L f(s)'
    meaning: "片側 Laplace 変換"
  - symbol: '\mathcal M f(s)'
    meaning: "Mellin 変換"
  - symbol: '\mathcal H_0 f(k)'
    meaning: "零次 Hankel 変換"
---

## この章の目標

- 周期信号の係数から不連続点近くの挙動を求める。
- Fourier 変換によって微分と畳み込みを乗算に移す。
- Laplace 変換で初期条件を含む微分方程式を解く。
- 対数座標、動径座標および標本化に適した変換を選ぶ。

本章では :ref[def-four-euclid-transform]{text="Fourier 変換"} を初期値問題や対称性を持つ積分に使う。周波数 $k$ は角波数であり、$x$ が長さなら単位は $\mathrm{m^{-1}}$ である。

## Fourier 級数と Gibbs 現象

### 周期の内側積分で係数を選ぶ

$2\pi$ 周期関数の :term[フーリエ係数]{en="Fourier coefficient" yomi="ふーりえけいすう"} を $c_n=(2\pi)^{-1}\int_{-\pi}^{\pi}f(x)e^{-inx}\dd x$ とする。$\int_{-\pi}^{\pi}e^{i(n-m)x}\dd x=2\pi\delta_{nm}$ なので、これは有限三角多項式で $e^{inx}$ の係数を取り出す。有限部分和は
$$
S_N(x)=\frac1{2\pi}\int_{-\pi}^\pi D_N(t)f(x-t)\dd t,\qquad
D_N(t)=\frac{\sin((N+1/2)t)}{\sin(t/2)} .
$$
ここで $D_N=\sum_{-N}^N e^{int}$ を有限等比和で計算した。$D_N$ は偶関数で、積分は $2\pi$ である。$f$ が区分 $C^1$ なら、左右極限の平均を引いた積分を $0<t<\pi$ にまとめる。振動因子に掛かる振幅は $[f(x+t)+f(x-t)-f(x+)-f(x-)]/\sin(t/2)$ であり、原点で有界な可積分関数である。振動積分が零に近づくことは、この振幅を $L^1$ で階段関数に近似すれば示せる。階段ごとの積分は周波数の逆数以下であり、近似誤差は周波数によらず抑えられる。従って部分和の点ごとの極限は左右の平均である。

:::example[Gibbs 現象の大きさ]
$f(x)=\sgn(\sin x)$ の正弦係数は奇数 $n$ で $4/(\pi n)$、偶数で零である。実際、奇関数なので余弦項は零、正弦係数は $(2/\pi)\int_0^\pi\sin(nx)\dd x$ である。最初の $N$ 個の奇数項の和 $T_N$ は $T_N(0)=0$、$T_N'=(2/\pi)\sin(2Nx)/\sin x$ を満たす。後者は $\cos((2j-1)x)\sin x=[\sin(2jx)-\sin((2j-2)x)]/2$ を足して得る。従って最初の正の極大は $x=\pi/(2N)$ にあり、
$$
T_N\left(\frac{\pi}{2N}\right)\longrightarrow
\frac2\pi\int_0^\pi\frac{\sin u}{u}\dd u
\approx1.17898.
$$
これは $u=2Nx$ で積分変数を替え、$(2N)\sin(u/(2N))$ を $u$ に置き換える極限である。最大点は跳びに近づくが、超過分は消えない。超過分は跳び幅 $2$ の約 $8.949\%$ である。これを :term[ギブス現象]{en="Gibbs phenomenon" yomi="ぎぶすげんしょう"} という。
:::

### Fejér 平均で跳びの外側の収束を改善する

部分和 $S_0,\ldots,S_{N-1}$ の平均は畳み込み核 $F_N(t)=N^{-1}|\sum_{j=0}^{N-1}e^{ijt}|^2$ を使う。二重和の中で添字差が $m$ となる対は $N-|m|$ 個なので、Fourier 係数は $1-|m|/N$（$|m|<N$）である。従って平均と核の式が一致する。$F_N\ge0$、積分 $2\pi$ であるため、$-1\le f\le1$ なら平均も同じ範囲に入り、矩形波の超過は生じない。一方、跳びの点では核が偶で、両側の質量はそれぞれ $\pi$ である。原点から離れると $F_N(t)\le 1/[N\sin^2(t/2)]$ だから、近傍の外の質量は零に近づく。これを左右の片側極限と組み合わせると、左右平均に収束する。例えば $f=3+2\cos x$ の $N=2$ の平均は $3+\cos x$ になる。各モードの重みが $1$ へ近づくので、連続関数なら近似単位の議論で一様に元へ戻る。

### 連続性と Fourier 係数の減衰

周期関数を区分的に積分すると、周期端の値が合うときは一度の部分積分で境界項が消える。さらに導関数が区分的に滑らかなら、その跳びが次の部分積分の境界項になる。従って関数自体の跳びは $1/n$、導関数の跳びは条件の下で $1/n^2$ の減衰に対応する。絶対和が有限なら Weierstrass の判定により Fourier 級数は一様に収束し、跳びの直後に有限割合だけ行き過ぎる現象は起きない。

:::example[折れ線の周期延長（1）]
区間 $[-\pi,\pi]$ の $|x|$ を周期延長する。第 $m$ 余弦係数を調べる。基準の値で計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

偶関数なので正弦係数は零である。余弦係数には $2/\pi$ を掛ける。

$$
\begin{aligned}
a_m\big|_{m=1} &= - \frac{4}{\pi}
\end{aligned}
$$

一階導関数の正側は $1$ であり、最初の部分積分に残る積分を求める。

$$
\begin{aligned}
\int_0^\pi\sin(mx)\dd x\big|_{m=1} &= 2
\end{aligned}
$$

定数項を除いた第 $m$ モードを元の変数に戻す。

$$
\begin{aligned}
a_m\cos(mx)\big|_{m=1} &= - \frac{4 \cos{\left(x \right)}}{\pi}
\end{aligned}
$$

観測点 $x$=$\frac{\pi}{3}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
a_m\cos(mx)\big|_{x=\frac{\pi}{3}} &= - \frac{2}{\pi} \\
\frac{\partial^1 (a_m\cos(mx))}{\partial x^1} &= \frac{4 \sin{\left(x \right)}}{\pi} \\
\left.\frac{\partial^1 (a_m\cos(mx))}{\partial x^1}\right|_{x=\frac{\pi}{3}} &= \frac{2 \sqrt{3}}{\pi}
\end{aligned}
$$
:::

:::example[折れ線の周期延長（2）]
区間 $[-\pi,\pi]$ の $|x|$ を周期延長する。第 $m$ 余弦係数を調べる。値を変えて依存性を調べる。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

偶関数なので正弦係数は零である。余弦係数には $2/\pi$ を掛ける。

$$
\begin{aligned}
a_m\big|_{m=2} &= 0
\end{aligned}
$$

一階導関数の正側は $1$ であり、最初の部分積分に残る積分を求める。

$$
\begin{aligned}
\int_0^\pi\sin(mx)\dd x\big|_{m=2} &= 0
\end{aligned}
$$

定数項を除いた第 $m$ モードを元の変数に戻す。

$$
\begin{aligned}
a_m\cos(mx)\big|_{m=2} &= 0
\end{aligned}
$$
:::

:::example[折れ線の周期延長（3）]
区間 $[-\pi,\pi]$ の $|x|$ を周期延長する。第 $m$ 余弦係数を調べる。さらに値を変えて比較する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

偶関数なので正弦係数は零である。余弦係数には $2/\pi$ を掛ける。

$$
\begin{aligned}
a_m\big|_{m=3} &= - \frac{4}{9 \pi}
\end{aligned}
$$

一階導関数の正側は $1$ であり、最初の部分積分に残る積分を求める。

$$
\begin{aligned}
\int_0^\pi\sin(mx)\dd x\big|_{m=3} &= \frac{2}{3}
\end{aligned}
$$

定数項を除いた第 $m$ モードを元の変数に戻す。

$$
\begin{aligned}
a_m\cos(mx)\big|_{m=3} &= - \frac{4 \cos{\left(3 x \right)}}{9 \pi}
\end{aligned}
$$

観測点 $x$=$\frac{\pi}{3}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
a_m\cos(mx)\big|_{x=\frac{\pi}{3}} &= \frac{4}{9 \pi} \\
\frac{\partial^1 (a_m\cos(mx))}{\partial x^1} &= \frac{4 \sin{\left(3 x \right)}}{3 \pi} \\
\left.\frac{\partial^1 (a_m\cos(mx))}{\partial x^1}\right|_{x=\frac{\pi}{3}} &= 0
\end{aligned}
$$
:::

:::example[折れ線の周期延長（4）]
区間 $[-\pi,\pi]$ の $|x|$ を周期延長する。第 $m$ 余弦係数を調べる。最後にこの値でも計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=4
$$

偶関数なので正弦係数は零である。余弦係数には $2/\pi$ を掛ける。

$$
\begin{aligned}
a_m\big|_{m=4} &= 0
\end{aligned}
$$

一階導関数の正側は $1$ であり、最初の部分積分に残る積分を求める。

$$
\begin{aligned}
\int_0^\pi\sin(mx)\dd x\big|_{m=4} &= 0
\end{aligned}
$$

定数項を除いた第 $m$ モードを元の変数に戻す。

$$
\begin{aligned}
a_m\cos(mx)\big|_{m=4} &= 0
\end{aligned}
$$
:::

## Fourier 変換の性質と畳み込み

### 平行移動と畳み込み

:::definition[Fourier 変換]{#def-mp-fourier}
$f\in L^1(\R)$ に対して :term[フーリエ変換]{en="Fourier transform" yomi="ふーりえへんかん"} を $\widehat f(k)=\int_\R f(x)e^{-ikx}\dd x$ で定める。さらに $\widehat f\in L^1(\R)$ なら、$f$ の連続点での逆変換は次式である。

$$
f(x)=\frac1{2\pi}\int_\R\widehat f(k)e^{ikx}\dd k.
\label{eq-mp-fourier-convention}
$$
:::

変数変換から $f(x-a)\mapsto e^{-ika}\widehat f(k)$、$e^{ibx}f(x)\mapsto\widehat f(k-b)$、$f(ax)\mapsto|a|^{-1}\widehat f(k/a)$（実数 $a\ne0$）を得る。急減少する滑らかな関数では部分積分によって $f'\mapsto ik\widehat f$、積分内微分によって $xf\mapsto i\partial_k\widehat f$ である。

畳み込み $(f*g)(x)=\int f(x-y)g(y)\dd y$ の変換は $\widehat f\widehat g$ となる。$f,g\in L^1$ なら積分の絶対値の二重積分は $\|f\|_1\|g\|_1$ であり、Fubini を使って $u=x-y$ と置けば証明できる。$f,g\in L^1$ で、$g$ が連続かつ $\widehat g\in L^1$ なら、$g$ の逆変換を $\int f(x)g(x)e^{-ikx}\dd x$ に代入できる。絶対値の二重積分は $\|f\|_1\|\widehat g\|_1/(2\pi)$ で抑えられ、順序を交換すると $\widehat{fg}=(2\pi)^{-1}\widehat f*\widehat g$ を得る。

:::example[箱の平均化]
$f=\mathbf1_{[-1,1]}$ では $\widehat f(k)=2\sin k/k$ である。$f*f$ は $[-1,1]$ と $[x-1,x+1]$ の共通長なので $(2-|x|)_+$ となる。変換すると $4\sin^2 k/k^2$ であり、$k=0$ では三角形の面積 $4$ を返す。
:::

### 逆変換を Gaussian で正則化する

$\widehat f$ に $e^{-\epsilon k^2}$ を掛けた逆積分は Fubini によって $f*G_\epsilon$ となる。ここで $G_\epsilon(x)=(4\pi\epsilon)^{-1/2}e^{-x^2/(4\epsilon)}$ である。この Fourier 対は Gaussian の積分を複素平行移動するか、微分と部分積分で変換側の一階方程式を解いて得る。$G_\epsilon$ は積分が $1$ で任意の原点近傍の外の質量が零に近づく。$f\in L^1$ の連続点で畳み込みを近傍と外に分けると、近傍は連続性、外は $\|f\|_1$ に $G_\epsilon$ の上限を掛けて $f*G_\epsilon\to f$ を得る。$\widehat f\in L^1$ なら逆積分側にも優収束を使え、式 :eqref[eq-mp-fourier-convention] が成立する。

### 広がりと畳み込みの面積

可積分関数の畳み込みでは Fubini の定理を使える。実際、絶対値を積分すると各関数の $L^1$ ノルムの積で抑えられる。二つの関数を正規化した密度として使えば、畳み込みは独立な変数の和の密度である。このため平均は足し算になり、中心化した二乗を積分すると交差項が零なので分散も足し算になる。Fourier 側では零周波数で積分を、原点での二階微分で二次モーメントを読める。

:::example[正規化された Gaussian の幅（1）]
分散 $m$ の Gaussian に分散 $1$ の Gaussian を畳み込む。基準の値で計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

各入力の変換は平方完成で求まる。

$$
\begin{aligned}
\widehat f(k)\big|_{m=1} &= e^{- \frac{k^{2}}{2}}
\end{aligned}
$$

もう一方の変換を掛けると、二次の指数の係数が加わる。

$$
\begin{aligned}
\widehat {f*g}(k)\big|_{m=1} &= e^{- k^{2}}
\end{aligned}
$$

同じ変換対を逆に使う。面積は $1$ のままで、幅は分散の平方根になる。

$$
\begin{aligned}
(f*g)(x)\big|_{m=1} &= \frac{e^{- \frac{x^{2}}{4}}}{2 \sqrt{\pi}}
\end{aligned}
$$

観測点 $k$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\widehat f(k)\big|_{k=0} &= 1 \\
\frac{\partial^1 (\widehat f(k))}{\partial k^1} &= - k e^{- \frac{k^{2}}{2}} \\
\left.\frac{\partial^1 (\widehat f(k))}{\partial k^1}\right|_{k=0} &= 0 \\
\frac{\partial^2 (\widehat f(k))}{\partial k^2} &= \left(k^{2} - 1\right) e^{- \frac{k^{2}}{2}} \\
\left.\frac{\partial^2 (\widehat f(k))}{\partial k^2}\right|_{k=0} &= -1
\end{aligned}
$$

観測点 $k$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\widehat {f*g}(k)\big|_{k=0} &= 1 \\
\frac{\partial^1 (\widehat {f*g}(k))}{\partial k^1} &= - 2 k e^{- k^{2}} \\
\left.\frac{\partial^1 (\widehat {f*g}(k))}{\partial k^1}\right|_{k=0} &= 0 \\
\frac{\partial^2 (\widehat {f*g}(k))}{\partial k^2} &= 2 \left(2 k^{2} - 1\right) e^{- k^{2}} \\
\left.\frac{\partial^2 (\widehat {f*g}(k))}{\partial k^2}\right|_{k=0} &= -2
\end{aligned}
$$
:::

:::example[正規化された Gaussian の幅（2）]
分散 $m$ の Gaussian に分散 $1$ の Gaussian を畳み込む。値を変えて依存性を調べる。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

各入力の変換は平方完成で求まる。

$$
\begin{aligned}
\widehat f(k)\big|_{m=2} &= e^{- k^{2}}
\end{aligned}
$$

もう一方の変換を掛けると、二次の指数の係数が加わる。

$$
\begin{aligned}
\widehat {f*g}(k)\big|_{m=2} &= e^{- \frac{3 k^{2}}{2}}
\end{aligned}
$$

同じ変換対を逆に使う。面積は $1$ のままで、幅は分散の平方根になる。

$$
\begin{aligned}
(f*g)(x)\big|_{m=2} &= \frac{\sqrt{6} e^{- \frac{x^{2}}{6}}}{6 \sqrt{\pi}}
\end{aligned}
$$

観測点 $k$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\widehat f(k)\big|_{k=0} &= 1 \\
\frac{\partial^1 (\widehat f(k))}{\partial k^1} &= - 2 k e^{- k^{2}} \\
\left.\frac{\partial^1 (\widehat f(k))}{\partial k^1}\right|_{k=0} &= 0 \\
\frac{\partial^2 (\widehat f(k))}{\partial k^2} &= 2 \left(2 k^{2} - 1\right) e^{- k^{2}} \\
\left.\frac{\partial^2 (\widehat f(k))}{\partial k^2}\right|_{k=0} &= -2
\end{aligned}
$$

観測点 $k$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\widehat {f*g}(k)\big|_{k=0} &= 1 \\
\frac{\partial^1 (\widehat {f*g}(k))}{\partial k^1} &= - 3 k e^{- \frac{3 k^{2}}{2}} \\
\left.\frac{\partial^1 (\widehat {f*g}(k))}{\partial k^1}\right|_{k=0} &= 0 \\
\frac{\partial^2 (\widehat {f*g}(k))}{\partial k^2} &= 3 \left(3 k^{2} - 1\right) e^{- \frac{3 k^{2}}{2}} \\
\left.\frac{\partial^2 (\widehat {f*g}(k))}{\partial k^2}\right|_{k=0} &= -3
\end{aligned}
$$
:::

:::example[正規化された Gaussian の幅（3）]
分散 $m$ の Gaussian に分散 $1$ の Gaussian を畳み込む。さらに値を変えて比較する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

各入力の変換は平方完成で求まる。

$$
\begin{aligned}
\widehat f(k)\big|_{m=3} &= e^{- \frac{3 k^{2}}{2}}
\end{aligned}
$$

もう一方の変換を掛けると、二次の指数の係数が加わる。

$$
\begin{aligned}
\widehat {f*g}(k)\big|_{m=3} &= e^{- 2 k^{2}}
\end{aligned}
$$

同じ変換対を逆に使う。面積は $1$ のままで、幅は分散の平方根になる。

$$
\begin{aligned}
(f*g)(x)\big|_{m=3} &= \frac{\sqrt{2} e^{- \frac{x^{2}}{8}}}{4 \sqrt{\pi}}
\end{aligned}
$$

観測点 $k$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\widehat f(k)\big|_{k=0} &= 1 \\
\frac{\partial^1 (\widehat f(k))}{\partial k^1} &= - 3 k e^{- \frac{3 k^{2}}{2}} \\
\left.\frac{\partial^1 (\widehat f(k))}{\partial k^1}\right|_{k=0} &= 0 \\
\frac{\partial^2 (\widehat f(k))}{\partial k^2} &= 3 \left(3 k^{2} - 1\right) e^{- \frac{3 k^{2}}{2}} \\
\left.\frac{\partial^2 (\widehat f(k))}{\partial k^2}\right|_{k=0} &= -3
\end{aligned}
$$

観測点 $k$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\widehat {f*g}(k)\big|_{k=0} &= 1 \\
\frac{\partial^1 (\widehat {f*g}(k))}{\partial k^1} &= - 4 k e^{- 2 k^{2}} \\
\left.\frac{\partial^1 (\widehat {f*g}(k))}{\partial k^1}\right|_{k=0} &= 0 \\
\frac{\partial^2 (\widehat {f*g}(k))}{\partial k^2} &= 4 \left(4 k^{2} - 1\right) e^{- 2 k^{2}} \\
\left.\frac{\partial^2 (\widehat {f*g}(k))}{\partial k^2}\right|_{k=0} &= -4
\end{aligned}
$$
:::

:::example[正規化された Gaussian の幅（4）]
分散 $m$ の Gaussian に分散 $1$ の Gaussian を畳み込む。最後にこの値でも計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=4
$$

各入力の変換は平方完成で求まる。

$$
\begin{aligned}
\widehat f(k)\big|_{m=4} &= e^{- 2 k^{2}}
\end{aligned}
$$

もう一方の変換を掛けると、二次の指数の係数が加わる。

$$
\begin{aligned}
\widehat {f*g}(k)\big|_{m=4} &= e^{- \frac{5 k^{2}}{2}}
\end{aligned}
$$

同じ変換対を逆に使う。面積は $1$ のままで、幅は分散の平方根になる。

$$
\begin{aligned}
(f*g)(x)\big|_{m=4} &= \frac{\sqrt{10} e^{- \frac{x^{2}}{10}}}{10 \sqrt{\pi}}
\end{aligned}
$$

観測点 $k$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\widehat f(k)\big|_{k=0} &= 1 \\
\frac{\partial^1 (\widehat f(k))}{\partial k^1} &= - 4 k e^{- 2 k^{2}} \\
\left.\frac{\partial^1 (\widehat f(k))}{\partial k^1}\right|_{k=0} &= 0 \\
\frac{\partial^2 (\widehat f(k))}{\partial k^2} &= 4 \left(4 k^{2} - 1\right) e^{- 2 k^{2}} \\
\left.\frac{\partial^2 (\widehat f(k))}{\partial k^2}\right|_{k=0} &= -4
\end{aligned}
$$

観測点 $k$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\widehat {f*g}(k)\big|_{k=0} &= 1 \\
\frac{\partial^1 (\widehat {f*g}(k))}{\partial k^1} &= - 5 k e^{- \frac{5 k^{2}}{2}} \\
\left.\frac{\partial^1 (\widehat {f*g}(k))}{\partial k^1}\right|_{k=0} &= 0 \\
\frac{\partial^2 (\widehat {f*g}(k))}{\partial k^2} &= 5 \left(5 k^{2} - 1\right) e^{- \frac{5 k^{2}}{2}} \\
\left.\frac{\partial^2 (\widehat {f*g}(k))}{\partial k^2}\right|_{k=0} &= -5
\end{aligned}
$$
:::

## Parseval–Plancherel の等式

### エネルギーと周波数の規格化

:::theorem[Parseval–Plancherel の等式]{#thm-mp-plancherel}
$\inner f g=\int f\overline g$ の規約で、急減少する滑らかな $f,g$ に対し $\inner f g=(2\pi)^{-1}\inner{\widehat f}{\widehat g}$ である。従って $(2\pi)^{-1/2}$ 倍した変換は $L^2(\R)$ 上のユニタリ作用素に延びる。
:::

:::proof
逆変換した $g$ を複素共役して内積に代入する。$f$ と $\widehat g$ は可積分なので順序を交換でき、内側の積分が $\widehat f$ となる。急減少する滑らかな関数は $L^2$ 稠密だから、等長写像を完備化して一意に延長できる。部分積分と積分内微分によってこの関数空間は変換と逆変換で保たれる。従って範囲は稠密かつ閉であり、全射となる。周期の場合は有限部分和の直交性から $\|S_N\|_2^2/(2\pi)=\sum_{|n|\le N}|c_n|^2$ である。三角多項式は連続周期関数を一様近似し、連続周期関数は $L^2$ で稠密である。直交射影 $S_N$ は最良近似なので $\|f-S_N\|_2\to0$ となり、周期版の等式 $\|f\|_2^2/(2\pi)=\sum_{n\in\Z}|c_n|^2$ を得る。指数関数の張る空間の稠密性は、非負な Fejér 核 $N^{-1}|\sum_{j=0}^{N-1}e^{ijx}|^2$ で連続周期関数を一様近似できることから従う。この核は平均 $1$ で、原点から離れると等比和の評価により $O(1/N)$ となる。
:::

:::example[二乗可積分な箱]
前節の箱のノルムの二乗は $2$ である。Plancherel に入れると $(2\pi)^{-1}\int_\R 4\sin^2 k/k^2\dd k=2$ より $\int_\R(\sin k/k)^2\dd k=\pi$ となる。
:::

### 規格化を二つの空間で合わせる

Fourier 変換が内積を保つためには、指定の規約では周波数側に $1/(2\pi)$ を掛ける。局在した関数の尺度を変えると、一方の $L^2$ ノルムが増え、他方にも同じ倍率が掛かる。空間ごとの比較では、まず積分する変数と測度を合わせる。例えば $k$ を通常周波数 $\nu$ に直して $k=2\pi\nu$ とすれば、周波数側の測度 $\dd k/(2\pi)$ は $\dd\nu$ となる。

:::example[指数で局在する関数の二乗積分（1）]
関数 $f(x)=e^{-m|x|}$ に対し、両側で全二乗を計算する。基準の値で計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

原点の両側に分けて、$f$ の Fourier 変換を足す。

$$
\begin{aligned}
\widehat f(k)\big|_{m=1} &= \frac{2}{k^{2} + 1}
\end{aligned}
$$

空間側の全二乗は偶関数を積分すれば求まる。

$$
\begin{aligned}
\int_\R |f(x)|^2\dd x\big|_{m=1} &= 1
\end{aligned}
$$

周波数側では二乗分母になる。$k=m v$ と置くと $\int_\R (1+v^2)^{-2}\dd v=\pi/2$ を使える。

$$
\begin{aligned}
\frac1{2\pi}\int_\R |\widehat f(k)|^2\dd k\big|_{m=1} &= 1
\end{aligned}
$$

観測点 $k$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\widehat f(k)\big|_{k=0} &= 2 \\
\frac{\partial^1 (\widehat f(k))}{\partial k^1} &= - \frac{4 k}{\left(k^{2} + 1\right)^{2}} \\
\left.\frac{\partial^1 (\widehat f(k))}{\partial k^1}\right|_{k=0} &= 0 \\
\frac{\partial^2 (\widehat f(k))}{\partial k^2} &= \frac{4 \left(3 k^{2} - 1\right)}{\left(k^{2} + 1\right)^{3}} \\
\left.\frac{\partial^2 (\widehat f(k))}{\partial k^2}\right|_{k=0} &= -4
\end{aligned}
$$
:::

:::example[指数で局在する関数の二乗積分（2）]
関数 $f(x)=e^{-m|x|}$ に対し、両側で全二乗を計算する。値を変えて依存性を調べる。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

原点の両側に分けて、$f$ の Fourier 変換を足す。

$$
\begin{aligned}
\widehat f(k)\big|_{m=2} &= \frac{4}{k^{2} + 4}
\end{aligned}
$$

空間側の全二乗は偶関数を積分すれば求まる。

$$
\begin{aligned}
\int_\R |f(x)|^2\dd x\big|_{m=2} &= \frac{1}{2}
\end{aligned}
$$

周波数側では二乗分母になる。$k=m v$ と置くと $\int_\R (1+v^2)^{-2}\dd v=\pi/2$ を使える。

$$
\begin{aligned}
\frac1{2\pi}\int_\R |\widehat f(k)|^2\dd k\big|_{m=2} &= \frac{1}{2}
\end{aligned}
$$

観測点 $k$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\widehat f(k)\big|_{k=0} &= 1 \\
\frac{\partial^1 (\widehat f(k))}{\partial k^1} &= - \frac{8 k}{\left(k^{2} + 4\right)^{2}} \\
\left.\frac{\partial^1 (\widehat f(k))}{\partial k^1}\right|_{k=0} &= 0 \\
\frac{\partial^2 (\widehat f(k))}{\partial k^2} &= \frac{8 \left(3 k^{2} - 4\right)}{\left(k^{2} + 4\right)^{3}} \\
\left.\frac{\partial^2 (\widehat f(k))}{\partial k^2}\right|_{k=0} &= - \frac{1}{2}
\end{aligned}
$$
:::

:::example[指数で局在する関数の二乗積分（3）]
関数 $f(x)=e^{-m|x|}$ に対し、両側で全二乗を計算する。さらに値を変えて比較する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

原点の両側に分けて、$f$ の Fourier 変換を足す。

$$
\begin{aligned}
\widehat f(k)\big|_{m=3} &= \frac{6}{k^{2} + 9}
\end{aligned}
$$

空間側の全二乗は偶関数を積分すれば求まる。

$$
\begin{aligned}
\int_\R |f(x)|^2\dd x\big|_{m=3} &= \frac{1}{3}
\end{aligned}
$$

周波数側では二乗分母になる。$k=m v$ と置くと $\int_\R (1+v^2)^{-2}\dd v=\pi/2$ を使える。

$$
\begin{aligned}
\frac1{2\pi}\int_\R |\widehat f(k)|^2\dd k\big|_{m=3} &= \frac{1}{3}
\end{aligned}
$$

観測点 $k$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\widehat f(k)\big|_{k=0} &= \frac{2}{3} \\
\frac{\partial^1 (\widehat f(k))}{\partial k^1} &= - \frac{12 k}{\left(k^{2} + 9\right)^{2}} \\
\left.\frac{\partial^1 (\widehat f(k))}{\partial k^1}\right|_{k=0} &= 0 \\
\frac{\partial^2 (\widehat f(k))}{\partial k^2} &= \frac{36 \left(k^{2} - 3\right)}{\left(k^{2} + 9\right)^{3}} \\
\left.\frac{\partial^2 (\widehat f(k))}{\partial k^2}\right|_{k=0} &= - \frac{4}{27}
\end{aligned}
$$
:::

:::example[指数で局在する関数の二乗積分（4）]
関数 $f(x)=e^{-m|x|}$ に対し、両側で全二乗を計算する。最後にこの値でも計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=4
$$

原点の両側に分けて、$f$ の Fourier 変換を足す。

$$
\begin{aligned}
\widehat f(k)\big|_{m=4} &= \frac{8}{k^{2} + 16}
\end{aligned}
$$

空間側の全二乗は偶関数を積分すれば求まる。

$$
\begin{aligned}
\int_\R |f(x)|^2\dd x\big|_{m=4} &= \frac{1}{4}
\end{aligned}
$$

周波数側では二乗分母になる。$k=m v$ と置くと $\int_\R (1+v^2)^{-2}\dd v=\pi/2$ を使える。

$$
\begin{aligned}
\frac1{2\pi}\int_\R |\widehat f(k)|^2\dd k\big|_{m=4} &= \frac{1}{4}
\end{aligned}
$$

観測点 $k$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\widehat f(k)\big|_{k=0} &= \frac{1}{2} \\
\frac{\partial^1 (\widehat f(k))}{\partial k^1} &= - \frac{16 k}{\left(k^{2} + 16\right)^{2}} \\
\left.\frac{\partial^1 (\widehat f(k))}{\partial k^1}\right|_{k=0} &= 0 \\
\frac{\partial^2 (\widehat f(k))}{\partial k^2} &= \frac{16 \left(3 k^{2} - 16\right)}{\left(k^{2} + 16\right)^{3}} \\
\left.\frac{\partial^2 (\widehat f(k))}{\partial k^2}\right|_{k=0} &= - \frac{1}{16}
\end{aligned}
$$
:::

## Laplace 変換と Bromwich 積分

### 収束を保証する縦線

$\mathcal L f(s)=\int_0^\infty e^{-st}f(t)\dd t$ を :term[ラプラス変換]{en="Laplace transform" yomi="らぷらすへんかん"} とする。$f$ が $[0,\infty)$ 上で局所可積分かつ $|f(t)|\le Ce^{at}$ を十分大きな $t$ で満たすと、$\Re s>a$ で収束する。$t^n e^{-(\Re s-a)t}$ が内部コンパクト集合で可積分だから、積分内微分で正則である。$\gamma>a$ を固定すると $\mathcal L f(\gamma+ik)$ は $f(t)e^{-\gamma t}\mathbf1_{t\ge0}$ の Fourier 変換である。その逆変換により、例えば変換が縦線で可積分で $f$ が $t>0$ で連続なら
$$
f(t)=\frac1{2\pi i}\int_{\gamma-i\infty}^{\gamma+i\infty}e^{st}\mathcal L f(s)\dd s.
\label{eq-mp-bromwich}
$$
これを Bromwich 積分という。対称に打ち切った積分は、縦線上での絶対可積分性がなくても収束する場合がある。例えば $g(t)=f(t)e^{-\gamma t}\mathbf1_{t\ge0}$ が $t$ の左右で区分 $C^1$ 級で、一階導関数が $t$ の近くで片側ごとに有界とする。打ち切った逆 Fourier 積分は Fubini により
$$
\frac1\pi\int_\R g(t-u)\frac{\sin(Ru)}u\dd u
$$
となる。$u=0$ の近くで片側極限を引けば、残りの $\sin(Ru)$ に掛かる係数は局所可積分である。$u=0$ を避けた部分も $g\in L^1$ により可積分である。従って階段関数での近似を使うと、残りの積分は $R\to\infty$ で零になる。定数部分には :ref[exm-mp-dirichlet-integral] の $\int_0^\infty\sin v\,\dd v/v=\pi/2$ を用いる。従って $t>0$ での対称打ち切りは $[f(t-)+f(t+)]/2$ を返す。$t=0$ に同じ議論を用いると、零延長の値との平均 $f(0+)/2$ を返す。

:::example[二位極の反転]
$\mathcal L f(s)=(s+3)^{-2}$ とする。$t>0$ で縦線の左を閉じると、左半円は $|e^{st}|\le e^{\gamma t}$ と分母の二次減衰により消える。$-3$ での留数は $\partial_s e^{st}|_{s=-3}=te^{-3t}$ で、これが $f(t)$ である。$t<0$ では右を閉じ、極がなく零になる。
:::

### 縦線をずらすと寄与の現れ方が分かる

例えば $F(s)=1/[(s+1)(s+2)]$ を $\gamma=0$ の線で反転する。$\gamma$ を $-3/2$ へずらすと、長方形の水平辺は二次の減衰で消え、通過した $s=-1$ の留数 $e^{-t}$ が加わる。さらに線を $-3$ にずらすと $s=-2$ の留数 $-e^{-2t}$ が加わる。この線は両方の極より左にあるため、$t>0$ で左に閉じると残りは零である。従って元の片側の解は $e^{-t}-e^{-2t}$ である。右の収束線から出発するのは、単に特異点を避けるだけでなく、$t<0$ に応答を持たない逆を指定するためでもある。

### 縦線と因果性の選択

片側 Laplace 変換では、$f(t)$ を負時刻で零に延ばすと $e^{-\sigma t}f(t)$ の Fourier 変換として反転できる。$\sigma$ はその延長が可積分になるように選ぶ。正時刻では $e^{st}$ が左半平面で減衰し、負時刻では右半平面で減衰する。従って同じ有理関数でも、縦線をどの極の右に取るかが時間側の台を決める。

:::example[多重極の逆変換（1）]
原点から左にずれた極の次数を変える。縦線は極より右に取る。基準の値で計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

周波数の記号を $s$ とする。変換には $-2$ にある極が一つだけある。

$$
\begin{aligned}
F(s)\big|_{m=1} &= \frac{1}{s + 2}
\end{aligned}
$$

正時刻では $e^{st}$ の Taylor 展開のうち極の位数より一つ下の係数を取る。

$$
\begin{aligned}
\operatorname{Res}_{s=-2}(e^{st}F(s))\big|_{m=1} &= e^{- 2 t}
\end{aligned}
$$

積分を戻して確認する。正時刻のモーメント積分は factorial を打ち消す。

$$
\begin{aligned}
\int_0^\infty e^{-s t}f(t)\dd t\big|_{m=1} &= \frac{1}{s + 2}
\end{aligned}
$$

観測点 $t$=$1$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\operatorname{Res}_{s=-2}(e^{st}F(s))\big|_{t=1} &= e^{-2} \\
\frac{\partial^1 (\operatorname{Res}_{s=-2}(e^{st}F(s)))}{\partial t^1} &= - 2 e^{- 2 t} \\
\left.\frac{\partial^1 (\operatorname{Res}_{s=-2}(e^{st}F(s)))}{\partial t^1}\right|_{t=1} &= - \frac{2}{e^{2}} \\
\frac{\partial^2 (\operatorname{Res}_{s=-2}(e^{st}F(s)))}{\partial t^2} &= 4 e^{- 2 t} \\
\left.\frac{\partial^2 (\operatorname{Res}_{s=-2}(e^{st}F(s)))}{\partial t^2}\right|_{t=1} &= \frac{4}{e^{2}}
\end{aligned}
$$
:::

:::example[多重極の逆変換（2）]
原点から左にずれた極の次数を変える。縦線は極より右に取る。値を変えて依存性を調べる。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

周波数の記号を $s$ とする。変換には $-2$ にある極が一つだけある。

$$
\begin{aligned}
F(s)\big|_{m=2} &= \frac{1}{\left(s + 2\right)^{2}}
\end{aligned}
$$

正時刻では $e^{st}$ の Taylor 展開のうち極の位数より一つ下の係数を取る。

$$
\begin{aligned}
\operatorname{Res}_{s=-2}(e^{st}F(s))\big|_{m=2} &= t e^{- 2 t}
\end{aligned}
$$

積分を戻して確認する。正時刻のモーメント積分は factorial を打ち消す。

$$
\begin{aligned}
\int_0^\infty e^{-s t}f(t)\dd t\big|_{m=2} &= \frac{1}{\left(s + 2\right)^{2}}
\end{aligned}
$$

観測点 $t$=$1$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\operatorname{Res}_{s=-2}(e^{st}F(s))\big|_{t=1} &= e^{-2} \\
\frac{\partial^1 (\operatorname{Res}_{s=-2}(e^{st}F(s)))}{\partial t^1} &= \left(1 - 2 t\right) e^{- 2 t} \\
\left.\frac{\partial^1 (\operatorname{Res}_{s=-2}(e^{st}F(s)))}{\partial t^1}\right|_{t=1} &= - \frac{1}{e^{2}} \\
\frac{\partial^2 (\operatorname{Res}_{s=-2}(e^{st}F(s)))}{\partial t^2} &= 4 \left(t - 1\right) e^{- 2 t} \\
\left.\frac{\partial^2 (\operatorname{Res}_{s=-2}(e^{st}F(s)))}{\partial t^2}\right|_{t=1} &= 0
\end{aligned}
$$
:::

:::example[多重極の逆変換（3）]
原点から左にずれた極の次数を変える。縦線は極より右に取る。さらに値を変えて比較する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

周波数の記号を $s$ とする。変換には $-2$ にある極が一つだけある。

$$
\begin{aligned}
F(s)\big|_{m=3} &= \frac{1}{\left(s + 2\right)^{3}}
\end{aligned}
$$

正時刻では $e^{st}$ の Taylor 展開のうち極の位数より一つ下の係数を取る。

$$
\begin{aligned}
\operatorname{Res}_{s=-2}(e^{st}F(s))\big|_{m=3} &= \frac{t^{2} e^{- 2 t}}{2}
\end{aligned}
$$

積分を戻して確認する。正時刻のモーメント積分は factorial を打ち消す。

$$
\begin{aligned}
\int_0^\infty e^{-s t}f(t)\dd t\big|_{m=3} &= \frac{1}{\left(s + 2\right)^{3}}
\end{aligned}
$$

観測点 $t$=$1$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\operatorname{Res}_{s=-2}(e^{st}F(s))\big|_{t=1} &= \frac{1}{2 e^{2}} \\
\frac{\partial^1 (\operatorname{Res}_{s=-2}(e^{st}F(s)))}{\partial t^1} &= t \left(1 - t\right) e^{- 2 t} \\
\left.\frac{\partial^1 (\operatorname{Res}_{s=-2}(e^{st}F(s)))}{\partial t^1}\right|_{t=1} &= 0 \\
\frac{\partial^2 (\operatorname{Res}_{s=-2}(e^{st}F(s)))}{\partial t^2} &= \left(2 t^{2} - 4 t + 1\right) e^{- 2 t} \\
\left.\frac{\partial^2 (\operatorname{Res}_{s=-2}(e^{st}F(s)))}{\partial t^2}\right|_{t=1} &= - \frac{1}{e^{2}}
\end{aligned}
$$
:::

:::example[多重極の逆変換（4）]
原点から左にずれた極の次数を変える。縦線は極より右に取る。最後にこの値でも計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=4
$$

周波数の記号を $s$ とする。変換には $-2$ にある極が一つだけある。

$$
\begin{aligned}
F(s)\big|_{m=4} &= \frac{1}{\left(s + 2\right)^{4}}
\end{aligned}
$$

正時刻では $e^{st}$ の Taylor 展開のうち極の位数より一つ下の係数を取る。

$$
\begin{aligned}
\operatorname{Res}_{s=-2}(e^{st}F(s))\big|_{m=4} &= \frac{t^{3} e^{- 2 t}}{6}
\end{aligned}
$$

積分を戻して確認する。正時刻のモーメント積分は factorial を打ち消す。

$$
\begin{aligned}
\int_0^\infty e^{-s t}f(t)\dd t\big|_{m=4} &= \frac{1}{\left(s + 2\right)^{4}}
\end{aligned}
$$

観測点 $t$=$1$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\operatorname{Res}_{s=-2}(e^{st}F(s))\big|_{t=1} &= \frac{1}{6 e^{2}} \\
\frac{\partial^1 (\operatorname{Res}_{s=-2}(e^{st}F(s)))}{\partial t^1} &= \frac{t^{2} \left(3 - 2 t\right) e^{- 2 t}}{6} \\
\left.\frac{\partial^1 (\operatorname{Res}_{s=-2}(e^{st}F(s)))}{\partial t^1}\right|_{t=1} &= \frac{1}{6 e^{2}} \\
\frac{\partial^2 (\operatorname{Res}_{s=-2}(e^{st}F(s)))}{\partial t^2} &= \frac{t \left(2 t^{2} - 6 t + 3\right) e^{- 2 t}}{3} \\
\left.\frac{\partial^2 (\operatorname{Res}_{s=-2}(e^{st}F(s)))}{\partial t^2}\right|_{t=1} &= - \frac{1}{3 e^{2}}
\end{aligned}
$$
:::

## Laplace 変換による常微分方程式の解法

### 初期値も変換に含める

部分積分から $\mathcal L f'=s\mathcal L f-f(0+)$、$\mathcal L f''=s^2\mathcal L f-sf(0+)-f'(0+)$ である。上端の境界項が消えるよう、関数と必要な導関数に対して $\Re s$ を十分大きく取る。片側畳み込み $\int_0^t f(t-v)g(v)\dd v$ の変換は、変数を $v,t-v$ として正の象限に移せば積である。時間遅延 $\mathbf1_{t\ge a}f(t-a)$ は $e^{-as}\mathcal L f(s)$ に対応する。

:::example[二つの時定数]
$y''+3y'+2y=1$、$y(0)=y'(0)=0$ の時間を無次元とする。$(s+1)(s+2)Y=1/s$ より $Y=\tfrac1{2s}-1/(s+1)+1/[2(s+2)]$ となる。従って $y=\tfrac12-e^{-t}+\tfrac12e^{-2t}$ である。$t=0$ の値と導関数はともに零であり、定常値は $1/2$ である。時間を秒に戻すときは、係数 $3,2$ にそれぞれ $\mathrm{s^{-1}}$ と $\mathrm{s^{-2}}$ が対応する。
:::

### 同じ極を持つ入力と応答

線形常微分方程式で零初期値とすると、解の変換は外力の変換と伝達関数の積になる。外力が自由解と同じ指数を持つと、二つの極は同じ点に重なる。したがって反転すると時間の冪が一つ増える。時間の側では、二つの指数を畳み込むと積分変数への依存性が打ち消し合い、その区間の長さ $t$ が残る。

### 極と長時間の極限

変位 $y(t)$ が有限の初期値に右から連続で指数増大以下なら、$\lim_{s\to+\infty}sY(s)=y(0+)$ である。$t$ を小区間と残りに分けると、小区間では $s e^{-st}$ が質量を原点に集め、残りは指数増大の上界を使って零にできる。同様に $y(t)$ が $t\to\infty$ で有限の $L$ に近づくなら、$\lim_{s\downarrow0}sY(s)=L$ となる。$y-L$ の積分を有限区間とその外に分ければ、有限区間は先頭の $s$ によって零に近づき、外側は $|y-L|$ の上限で抑えられる。

逆向きには極の位置を調べる。例えば厳密に真分数の有理関数 $Y$ が、$s=0$ の単純極と $\Re s<0$ の極だけを持つとする。部分分数分解で原点の項は $L/s$、他の項は $(s-p)^{-j}$ の定数倍になる。逆変換では後者が $e^{pt}t^{j-1}/(j-1)!$ になり、全て零へ減衰する。従ってこの条件で最終値は $\lim_{s\downarrow0}sY(s)$ である。例えば $Y=(s^2+1)^{-1}$ ではこの代数的な極限が零でも、逆変換の $\sin t$ は振動を続ける。虚軸上の二つの極が残るためである。

:::example[指数外力が自由解と一致する場合（1）]
初期値 $y(0)=0$ の問題 $y^\prime+my=e^{-mt}$ を解く。基準の値で計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

外力の変換と一次の演算子から応答の変換を作る。

$$
\begin{aligned}
Y(s)\big|_{m=1} &= \frac{1}{\left(s + 1\right)^{2}}
\end{aligned}
$$

二重極を反転すると、自由解と同じ指数に $t$ が掛かる。

$$
\begin{aligned}
y(t)\big|_{m=1} &= t e^{- t}
\end{aligned}
$$

微分と零階項を足して元の外力に戻ることを確認する。

$$
\begin{aligned}
y^\prime+my\big|_{m=1} &= e^{- t}
\end{aligned}
$$

観測点 $t$=$1$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
y(t)\big|_{t=1} &= e^{-1} \\
\frac{\partial^1 (y(t))}{\partial t^1} &= \left(1 - t\right) e^{- t} \\
\left.\frac{\partial^1 (y(t))}{\partial t^1}\right|_{t=1} &= 0 \\
\frac{\partial^2 (y(t))}{\partial t^2} &= \left(t - 2\right) e^{- t} \\
\left.\frac{\partial^2 (y(t))}{\partial t^2}\right|_{t=1} &= - \frac{1}{e}
\end{aligned}
$$

観測点 $t$=$1$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
y^\prime+my\big|_{t=1} &= e^{-1} \\
\frac{\partial^1 (y^\prime+my)}{\partial t^1} &= - e^{- t} \\
\left.\frac{\partial^1 (y^\prime+my)}{\partial t^1}\right|_{t=1} &= - \frac{1}{e} \\
\frac{\partial^2 (y^\prime+my)}{\partial t^2} &= e^{- t} \\
\left.\frac{\partial^2 (y^\prime+my)}{\partial t^2}\right|_{t=1} &= e^{-1}
\end{aligned}
$$
:::

:::example[指数外力が自由解と一致する場合（2）]
初期値 $y(0)=0$ の問題 $y^\prime+my=e^{-mt}$ を解く。値を変えて依存性を調べる。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

外力の変換と一次の演算子から応答の変換を作る。

$$
\begin{aligned}
Y(s)\big|_{m=2} &= \frac{1}{\left(s + 2\right)^{2}}
\end{aligned}
$$

二重極を反転すると、自由解と同じ指数に $t$ が掛かる。

$$
\begin{aligned}
y(t)\big|_{m=2} &= t e^{- 2 t}
\end{aligned}
$$

微分と零階項を足して元の外力に戻ることを確認する。

$$
\begin{aligned}
y^\prime+my\big|_{m=2} &= e^{- 2 t}
\end{aligned}
$$

観測点 $t$=$1$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
y(t)\big|_{t=1} &= e^{-2} \\
\frac{\partial^1 (y(t))}{\partial t^1} &= \left(1 - 2 t\right) e^{- 2 t} \\
\left.\frac{\partial^1 (y(t))}{\partial t^1}\right|_{t=1} &= - \frac{1}{e^{2}} \\
\frac{\partial^2 (y(t))}{\partial t^2} &= 4 \left(t - 1\right) e^{- 2 t} \\
\left.\frac{\partial^2 (y(t))}{\partial t^2}\right|_{t=1} &= 0
\end{aligned}
$$

観測点 $t$=$1$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
y^\prime+my\big|_{t=1} &= e^{-2} \\
\frac{\partial^1 (y^\prime+my)}{\partial t^1} &= - 2 e^{- 2 t} \\
\left.\frac{\partial^1 (y^\prime+my)}{\partial t^1}\right|_{t=1} &= - \frac{2}{e^{2}} \\
\frac{\partial^2 (y^\prime+my)}{\partial t^2} &= 4 e^{- 2 t} \\
\left.\frac{\partial^2 (y^\prime+my)}{\partial t^2}\right|_{t=1} &= \frac{4}{e^{2}}
\end{aligned}
$$
:::

:::example[指数外力が自由解と一致する場合（3）]
初期値 $y(0)=0$ の問題 $y^\prime+my=e^{-mt}$ を解く。さらに値を変えて比較する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

外力の変換と一次の演算子から応答の変換を作る。

$$
\begin{aligned}
Y(s)\big|_{m=3} &= \frac{1}{\left(s + 3\right)^{2}}
\end{aligned}
$$

二重極を反転すると、自由解と同じ指数に $t$ が掛かる。

$$
\begin{aligned}
y(t)\big|_{m=3} &= t e^{- 3 t}
\end{aligned}
$$

微分と零階項を足して元の外力に戻ることを確認する。

$$
\begin{aligned}
y^\prime+my\big|_{m=3} &= e^{- 3 t}
\end{aligned}
$$

観測点 $t$=$1$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
y(t)\big|_{t=1} &= e^{-3} \\
\frac{\partial^1 (y(t))}{\partial t^1} &= \left(1 - 3 t\right) e^{- 3 t} \\
\left.\frac{\partial^1 (y(t))}{\partial t^1}\right|_{t=1} &= - \frac{2}{e^{3}} \\
\frac{\partial^2 (y(t))}{\partial t^2} &= 3 \left(3 t - 2\right) e^{- 3 t} \\
\left.\frac{\partial^2 (y(t))}{\partial t^2}\right|_{t=1} &= \frac{3}{e^{3}}
\end{aligned}
$$

観測点 $t$=$1$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
y^\prime+my\big|_{t=1} &= e^{-3} \\
\frac{\partial^1 (y^\prime+my)}{\partial t^1} &= - 3 e^{- 3 t} \\
\left.\frac{\partial^1 (y^\prime+my)}{\partial t^1}\right|_{t=1} &= - \frac{3}{e^{3}} \\
\frac{\partial^2 (y^\prime+my)}{\partial t^2} &= 9 e^{- 3 t} \\
\left.\frac{\partial^2 (y^\prime+my)}{\partial t^2}\right|_{t=1} &= \frac{9}{e^{3}}
\end{aligned}
$$
:::

:::example[指数外力が自由解と一致する場合（4）]
初期値 $y(0)=0$ の問題 $y^\prime+my=e^{-mt}$ を解く。最後にこの値でも計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=4
$$

外力の変換と一次の演算子から応答の変換を作る。

$$
\begin{aligned}
Y(s)\big|_{m=4} &= \frac{1}{\left(s + 4\right)^{2}}
\end{aligned}
$$

二重極を反転すると、自由解と同じ指数に $t$ が掛かる。

$$
\begin{aligned}
y(t)\big|_{m=4} &= t e^{- 4 t}
\end{aligned}
$$

微分と零階項を足して元の外力に戻ることを確認する。

$$
\begin{aligned}
y^\prime+my\big|_{m=4} &= e^{- 4 t}
\end{aligned}
$$

観測点 $t$=$1$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
y(t)\big|_{t=1} &= e^{-4} \\
\frac{\partial^1 (y(t))}{\partial t^1} &= \left(1 - 4 t\right) e^{- 4 t} \\
\left.\frac{\partial^1 (y(t))}{\partial t^1}\right|_{t=1} &= - \frac{3}{e^{4}} \\
\frac{\partial^2 (y(t))}{\partial t^2} &= 8 \left(2 t - 1\right) e^{- 4 t} \\
\left.\frac{\partial^2 (y(t))}{\partial t^2}\right|_{t=1} &= \frac{8}{e^{4}}
\end{aligned}
$$

観測点 $t$=$1$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
y^\prime+my\big|_{t=1} &= e^{-4} \\
\frac{\partial^1 (y^\prime+my)}{\partial t^1} &= - 4 e^{- 4 t} \\
\left.\frac{\partial^1 (y^\prime+my)}{\partial t^1}\right|_{t=1} &= - \frac{4}{e^{4}} \\
\frac{\partial^2 (y^\prime+my)}{\partial t^2} &= 16 e^{- 4 t} \\
\left.\frac{\partial^2 (y^\prime+my)}{\partial t^2}\right|_{t=1} &= \frac{16}{e^{4}}
\end{aligned}
$$
:::

## Mellin 変換と Hankel 変換

### 対数と動径に適した基底

:term[メリン変換]{en="Mellin transform" yomi="めりんへんかん"} を $\mathcal M f(s)=\int_0^\infty f(x)x^{s-1}\dd x$ とする。$x=e^u$ とすれば $\mathcal M f(c+ik)=\int e^{cu} f(e^u)e^{iku}\dd u$ である。収束する $c$ を固定して Fourier 逆変換を使うと $f(x)=(2\pi i)^{-1}\int_{c-i\infty}^{c+i\infty}\mathcal M f(s)x^{-s}\dd s$ を得る。成立条件は対数座標での Fourier 逆変換の条件である。例えば $\Re s>0$ で $\mathcal M(e^{-x})(s)=\Gamma(s)=\int_0^\infty x^{s-1}e^{-x}\dd x$ であり、$\mathcal M(\mathbf1_{(0,1)})(s)=1/s$ である。拡大 $f(ax)$（$a>0$）は $a^{-s}$ を掛け、$xf'$ は境界項が消える条件のもとで $-s$ を掛ける。

二次元平面で急減少する滑らかな動径関数 $f(r)$ を考える。このとき、角積分 $(2\pi)^{-1}\int_0^{2\pi}e^{-iz\cos\theta}\dd\theta=J_0(z)$ を使う。これは [Bessel 関数](/physics/mathematical-physics/bessel-functions/) の零次関数である。角積分を冪級数にすれば、奇数項が零、偶数項の平均が $\binom{2m}m 2^{-2m}$ であるため $J_0(z)=\sum(-z^2/4)^m/(m!)^2$ となる。二次元 Fourier 変換は $2\pi\int_0^\infty f(r)J_0(kr)r\dd r$ である。そこでこの積分の $2\pi$ を除いたものを零次 :term[ハンケル変換]{en="Hankel transform" yomi="はんけるへんかん"} と呼ぶ。逆 Fourier 変換でも角積分して $f(r)=\int_0^\infty \mathcal H_0f(k) J_0(kr)k\dd k$ となる。

:::example[動径 Gaussian]
$f(r)=e^{-ar^2}$（$a>0$）は直交座標で二つの Gaussian の積なので Fourier 変換は $\pi a^{-1}e^{-k^2/(4a)}$ である。従って Hankel 変換は $(2a)^{-1}e^{-k^2/(4a)}$、$a=2,k=0$ で $1/4$ である。
:::

### 角度モードごとの Hankel 変換

二次元の場が $f(r)e^{im\theta}$（$m\ge0$ は整数）の場合、周波数側の偏角を $\psi$ として $\theta-\psi=\alpha$ と置く。角積分の公式
$$
\frac1{2\pi}\int_0^{2\pi}e^{-iz\cos\alpha}e^{im\alpha}\dd\alpha
=(-i)^m J_m(z),\qquad
J_m(z)=\sum_{j=0}^\infty\frac{(-1)^j(z/2)^{2j+m}}{j!(j+m)!}
$$
は、$\exp[-iz e^{i\alpha}/2]\exp[-iz e^{-i\alpha}/2]$ の二つの冪級数を掛けて、指数が相殺する項だけを選べば得られる。$\alpha$ の積分区間がコンパクトなので、固定した $z$ に対する二重級数は絶対一様収束する。

従って二次元 Fourier 変換は $2\pi(-i)^m e^{im\psi}\mathcal H_m f(k)$ となる。ただし $\mathcal H_m f(k)=\int_0^\infty f(r)J_m(kr)r\dd r$ である。元の場が平面上で滑らかかつ急減少するなら、逆変換で角積分を再び行うと $i^m$ が出て位相が消える。その結果
$$
f(r)=\int_0^\infty \mathcal H_m f(k)J_m(kr)k\dd k .
$$
さらに二次元の Parseval の式で角度を先に積分すると $\int_0^\infty|f(r)|^2r\dd r=\int_0^\infty|\mathcal H_m f(k)|^2k\dd k$ である。この二次元の等式は一変数の等式を二方向に順に適用して得る。

平面の Laplacian をモードへ作用させると、動径側に $f''+r^{-1}f'-m^2r^{-2}f$ が残る。二次元 Fourier 変換で Laplacian は $-k^2$ の乗算になるため、この組合せを Hankel 変換しても $-k^2\mathcal H_m f$ となる。例えば $(x+iy)e^{-ar^2}$ の Fourier 変換は、Gaussian の変換に $i(\partial_{k_x}+i\partial_{k_y})$ を作用させて求められる。角度因子を除けば $\mathcal H_1(re^{-ar^2})(k)=k e^{-k^2/(4a)}/(4a^2)$ である。$k=0$ で値が零なのは、元の場の角方向の平均が零だからである。

### 拡大の畳み込みを Mellin 変換で乗算にする

$x>0$ で $(f*_Mg)(x)=\int_0^\infty f(x/y)g(y)\dd y/y$ と定める。ある実数 $c$ で $\int |f(x)|x^{c-1}\dd x$ と $g$ の同じ積分が有限なら、$\Re s=c$ の二重積分は絶対収束する。内側で $x=uy$ と置くと $x^{s-1}\dd x=u^{s-1}y^s\dd u$ である。従って $\mathcal M(f*_Mg)(s)=\mathcal M f(s)\mathcal M g(s)$ となる。これは $x=e^u$ で拡大が加法に移ることの積分表示である。

例えば $f=g=\mathbf1_{(0,1)}$ なら二つが同時に非零なのは $x<y<1$ であり、$0<x<1$ で畳み込みは $-\log x$、$x\ge1$ で零である。Mellin 変換は $\Re s>0$ で $1/s^2$ となり、$1/s$ の積と一致する。二重極を反転すると対数因子が出る。

### Mellin 変数を微分する

重み $x^{s-1}$ の $s$ を微分すると $\log x$ が掛かる。積分の微分を交換するには両端の評価が必要である。$e^{-mx}$ なら無限遠で指数減衰する。原点で実部 $\Re s$ の下限を正とすれば、$x^{\Re s-1}|\log x|^j$ は任意の有限の $j$ で可積分である。このため同じ収束帯の内部で対数を掛けた積分を計算できる。

次の例では、実数 $s>0$ での $\Gamma$ の対数導関数を $\operatorname{polygamma}(0,s)=\Gamma'(s)/\Gamma(s)$ と書く。これは :term[ディガンマ関数]{en="digamma function" yomi="でぃがんまかんすう"} である。$s=2$ における値は次の積分から求められる。
$$
\int_0^N x^{s-1}(1-x/N)^N\dd x
=\frac{N^sN!}{s(s+1)\cdots(s+N)}.
$$
$x=Ny$ とした後、$y^{s-1}$ の原始関数 $y^s/s$ を使って部分積分を繰り返すと右辺を得る。$[0,N]$ 上で $(1-x/N)^N\le e^{-x}$ だから、左辺とその $s$ 微分はそれぞれ $\Gamma(s)$ と $\Gamma'(s)$ に収束する。従って対数導関数は $\lim_{N\to\infty}(\log N-\sum_{j=0}^N1/(s+j))$ となる。Euler の定数を $\gamma=\lim_{N\to\infty}(\sum_{j=1}^N1/j-\log N)$ とすれば、$s=2$ でこの値は $1-\gamma$ である。$\gamma$ の極限は積分 $\int_1^N\dd x/x$ と調和和の差を比較すると分かる。具体的に差は下に有界であり、次の項を足したときの変化 $1/(N+1)-\int_N^{N+1}\dd x/x$ は負である。

:::example[対数を一つ含む Mellin 積分（1）]
指数関数を拡大し、Gamma 関数の微分から対数重みの積分を得る。基準の値で計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

まず $v=mx$ と置いて変換する。

$$
\begin{aligned}
\mathcal M f(s)\big|_{m=1} &= \Gamma\left(s\right)
\end{aligned}
$$

この式を $s$ で微分する。関数の微分と拡大の微分が両方現れる。

$$
\begin{aligned}
\mathcal M(f\log x)(s)\big|_{m=1} &= \Gamma\left(s\right) \operatorname{polygamma}{\left(0,s \right)}
\end{aligned}
$$

拡大しなかった場合との差は、$\log(mx)=\log m+\log x$ に対応する。

$$
\begin{aligned}
\mathcal M(f\log(mx))(s)\big|_{m=1} &= \Gamma\left(s\right) \operatorname{polygamma}{\left(0,s \right)}
\end{aligned}
$$

観測点 $s$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\mathcal M f(s)\big|_{s=2} &= 1
\end{aligned}
$$

観測点 $s$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\mathcal M(f\log x)(s)\big|_{s=2} &= 1 - \gamma
\end{aligned}
$$

観測点 $s$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\mathcal M(f\log(mx))(s)\big|_{s=2} &= 1 - \gamma
\end{aligned}
$$
:::

:::example[対数を一つ含む Mellin 積分（2）]
指数関数を拡大し、Gamma 関数の微分から対数重みの積分を得る。値を変えて依存性を調べる。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

まず $v=mx$ と置いて変換する。

$$
\begin{aligned}
\mathcal M f(s)\big|_{m=2} &= 2^{- s} \Gamma\left(s\right)
\end{aligned}
$$

この式を $s$ で微分する。関数の微分と拡大の微分が両方現れる。

$$
\begin{aligned}
\mathcal M(f\log x)(s)\big|_{m=2} &= 2^{- s} \left(\operatorname{polygamma}{\left(0,s \right)} - \log{\left(2 \right)}\right) \Gamma\left(s\right)
\end{aligned}
$$

拡大しなかった場合との差は、$\log(mx)=\log m+\log x$ に対応する。

$$
\begin{aligned}
\mathcal M(f\log(mx))(s)\big|_{m=2} &= 2^{- s} \Gamma\left(s\right) \operatorname{polygamma}{\left(0,s \right)}
\end{aligned}
$$

観測点 $s$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\mathcal M f(s)\big|_{s=2} &= \frac{1}{4}
\end{aligned}
$$

観測点 $s$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\mathcal M(f\log x)(s)\big|_{s=2} &= - \frac{\log{\left(2 \right)}}{4} - \frac{\gamma}{4} + \frac{1}{4}
\end{aligned}
$$

観測点 $s$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\mathcal M(f\log(mx))(s)\big|_{s=2} &= \frac{1}{4} - \frac{\gamma}{4}
\end{aligned}
$$
:::

:::example[対数を一つ含む Mellin 積分（3）]
指数関数を拡大し、Gamma 関数の微分から対数重みの積分を得る。さらに値を変えて比較する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

まず $v=mx$ と置いて変換する。

$$
\begin{aligned}
\mathcal M f(s)\big|_{m=3} &= 3^{- s} \Gamma\left(s\right)
\end{aligned}
$$

この式を $s$ で微分する。関数の微分と拡大の微分が両方現れる。

$$
\begin{aligned}
\mathcal M(f\log x)(s)\big|_{m=3} &= 3^{- s} \left(\operatorname{polygamma}{\left(0,s \right)} - \log{\left(3 \right)}\right) \Gamma\left(s\right)
\end{aligned}
$$

拡大しなかった場合との差は、$\log(mx)=\log m+\log x$ に対応する。

$$
\begin{aligned}
\mathcal M(f\log(mx))(s)\big|_{m=3} &= 3^{- s} \Gamma\left(s\right) \operatorname{polygamma}{\left(0,s \right)}
\end{aligned}
$$

観測点 $s$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\mathcal M f(s)\big|_{s=2} &= \frac{1}{9}
\end{aligned}
$$

観測点 $s$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\mathcal M(f\log x)(s)\big|_{s=2} &= - \frac{\log{\left(3 \right)}}{9} - \frac{\gamma}{9} + \frac{1}{9}
\end{aligned}
$$

観測点 $s$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\mathcal M(f\log(mx))(s)\big|_{s=2} &= \frac{1}{9} - \frac{\gamma}{9}
\end{aligned}
$$
:::

:::example[対数を一つ含む Mellin 積分（4）]
指数関数を拡大し、Gamma 関数の微分から対数重みの積分を得る。最後にこの値でも計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=4
$$

まず $v=mx$ と置いて変換する。

$$
\begin{aligned}
\mathcal M f(s)\big|_{m=4} &= 4^{- s} \Gamma\left(s\right)
\end{aligned}
$$

この式を $s$ で微分する。関数の微分と拡大の微分が両方現れる。

$$
\begin{aligned}
\mathcal M(f\log x)(s)\big|_{m=4} &= 4^{- s} \left(\operatorname{polygamma}{\left(0,s \right)} - \log{\left(4 \right)}\right) \Gamma\left(s\right)
\end{aligned}
$$

拡大しなかった場合との差は、$\log(mx)=\log m+\log x$ に対応する。

$$
\begin{aligned}
\mathcal M(f\log(mx))(s)\big|_{m=4} &= 4^{- s} \Gamma\left(s\right) \operatorname{polygamma}{\left(0,s \right)}
\end{aligned}
$$

観測点 $s$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\mathcal M f(s)\big|_{s=2} &= \frac{1}{16}
\end{aligned}
$$

観測点 $s$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\mathcal M(f\log x)(s)\big|_{s=2} &= - \frac{\log{\left(4 \right)}}{16} - \frac{\gamma}{16} + \frac{1}{16}
\end{aligned}
$$

観測点 $s$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\mathcal M(f\log(mx))(s)\big|_{s=2} &= \frac{1}{16} - \frac{\gamma}{16}
\end{aligned}
$$
:::

## 離散 Fourier 変換と標本化定理

### 有限個の標本と周波数の重なり

長さ $N$ の列 $x_j$ の :term[離散フーリエ変換]{en="discrete Fourier transform" yomi="りさんふーりえへんかん"} を $X_m=\sum_{j=0}^{N-1}x_j e^{-2\pi ijm/N}$ とする。$\sum_{m=0}^{N-1}e^{2\pi i(j-\ell)m/N}=N\delta_{j\ell}$ は有限等比和なので、逆は $x_j=N^{-1}\sum_m X_m e^{2\pi ijm/N}$ である。各添字を $N$ 周期で読むと、同じ計算により巡回畳み込みの変換が積、$\sum|x_j|^2=N^{-1}\sum|X_m|^2$ が成り立つ。

:::proposition[帯域制限信号の標本化]{#prop-mp-sampling}
$T>0$ とし、$\widehat f\in L^2(\R)$ が $[-\pi/T,\pi/T]$ に台を持ち、$f$ は逆変換で定めた連続代表元とする。すると
$$
f(t)=\sum_{n\in\Z} f(nT)\frac{\sin[\pi(t/T-n)]}{\pi(t/T-n)}
$$
となる。分母零の値は連続延長で $1$ である。
:::

:::proof
$\widehat f$ をこの区間で $e^{-iknT}$ の Fourier 級数に展開する。その係数は $T/(2\pi)\int\widehat f(k)e^{iknT}\dd k=T f(nT)$ である。逆積分して各指数項を計算すると表示式を得る。周波数の区間の長さが有限なので $L^2$ 収束は $L^1$ 収束を含意し、その逆積分は全 $t$ に一様収束する。
:::

標本で $e^{iknT}$ と $e^{i(k+2\pi/T)nT}$ は一致する。従って帯域を制限しなければ異なる周波数を区別できない。この折り返しを :term[エイリアシング]{en="aliasing" yomi="えいりあしんぐ"} という。例えば $T=1\,\mathrm{ms}$ では通常周波数 $200\,\mathrm{Hz}$ と $1200\,\mathrm{Hz}$ の複素指数は同じ標本を作る。

帯域端の角波数を $K$ とすると、上の命題で再構成できる標本間隔は $T\le\pi/K$ である。角波数 $k$ を通常周波数 $\nu=k/(2\pi)$ で書けば、標本化周波数 $1/T$ は最大周波数 $K/(2\pi)$ の二倍以上である。等号のとき、命題では $L^2$ のスペクトルを仮定するため端点だけに集中する成分は含まれない。例えば $f(t)=\sin(\pi t/T)$ は全標本が零であり、二つの帯域端の線スペクトルを許す場合はこの条件が必要になる。

### 周期化と Poisson の和公式

急減少する滑らかな $f$ を $P(x)=\sum_{n\in\Z}f(x+nT)$ で周期化する。和は導関数も含めて局所一様収束し、周期の積分を和の内側に移すと係数は $T^{-1}\widehat f(2\pi m/T)$ になる。変換も急減少するので級数は絶対収束する。その $x=0$ の値は
$$
\sum_{n\in\Z}f(nT)=\frac1T\sum_{m\in\Z}\widehat f(2\pi m/T)\label{eq-mp-poisson-sum}
$$
である。例えば $f(x)=e^{-\pi x^2}$ では $\widehat f(k)=e^{-k^2/(4\pi)}$ なので $\sum e^{-\pi n^2T^2}=T^{-1}\sum e^{-\pi m^2/T^2}$ となる。$T$ が小さいと右辺は $T^{-1}$ にきわめて近く、左の多数の項を短い和で計算できる。

偶数長 $N$ の DFT では偶数の標本と奇数の標本を分け、長さ $N/2$ の変換をそれぞれ $E_m,O_m$ とする。定義を代入して $X_m=E_m+e^{-2\pi im/N}O_m$、$X_{m+N/2}=E_m-e^{-2\pi im/N}O_m$ を得る。$N$ が二の冪なら同じ操作を反復でき、各段で $O(N)$ の演算、段数 $\log_2N$ である。例えば $(1,2,3,4)$ では $E=(4,-2)$、$O=(6,-2)$ なので $X=(10,-2+2i,-2,-2-2i)$ である。

### 余弦の二つのビンと折り返し

$N$ 個の標本を取ると、指数 $e^{2\pi i j n/N}$ は $j$ と $j+N$ で同じ値になる。余弦では正・負周波数の指数を半分ずつ足すため、二つの非零係数が出る。ただし $j=0$ と、偶数の $N$ での $j=N/2$ では二つのビンが重なる。同じ信号が異なる連続周波数を持ち得るのは、標本化で区別できない指数が存在するためである。

:::example[八点での余弦の DFT（1）]
番号 $n$ を $0$ から $7$ まで取り、入力 $v_n=\cos(2\pi m n/8)$ の係数を求める。基準の値で計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

標本値を番号順に並べる。

$$
\begin{aligned}
(v_n)_{n=0}^7\big|_{m=1} &= \left( 1, \  \frac{\sqrt{2}}{2}, \  0, \  - \frac{\sqrt{2}}{2}, \  -1, \  - \frac{\sqrt{2}}{2}, \  0, \  \frac{\sqrt{2}}{2}\right)
\end{aligned}
$$

二つの等比級数を足すと、添字 $m$ と $8-m$ のビンに $4$ ずつ入る。添字が重なれば足す。

$$
\begin{aligned}
(\widehat v_j)_{j=0}^7\big|_{m=1} &= \left( 0, \  4, \  0, \  0, \  0, \  0, \  0, \  4\right)
\end{aligned}
$$

有限の Parseval の等式では右辺の周波数二乗和を $8$ で割る。

$$
\begin{aligned}
\sum_{n=0}^7|v_n|^2\big|_{m=1} &= 4
\end{aligned}
$$
:::

:::example[八点での余弦の DFT（2）]
番号 $n$ を $0$ から $7$ まで取り、入力 $v_n=\cos(2\pi m n/8)$ の係数を求める。値を変えて依存性を調べる。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

標本値を番号順に並べる。

$$
\begin{aligned}
(v_n)_{n=0}^7\big|_{m=2} &= \left( 1, \  0, \  -1, \  0, \  1, \  0, \  -1, \  0\right)
\end{aligned}
$$

二つの等比級数を足すと、添字 $m$ と $8-m$ のビンに $4$ ずつ入る。添字が重なれば足す。

$$
\begin{aligned}
(\widehat v_j)_{j=0}^7\big|_{m=2} &= \left( 0, \  0, \  4, \  0, \  0, \  0, \  4, \  0\right)
\end{aligned}
$$

有限の Parseval の等式では右辺の周波数二乗和を $8$ で割る。

$$
\begin{aligned}
\sum_{n=0}^7|v_n|^2\big|_{m=2} &= 4
\end{aligned}
$$
:::

:::example[八点での余弦の DFT（3）]
番号 $n$ を $0$ から $7$ まで取り、入力 $v_n=\cos(2\pi m n/8)$ の係数を求める。さらに値を変えて比較する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

標本値を番号順に並べる。

$$
\begin{aligned}
(v_n)_{n=0}^7\big|_{m=3} &= \left( 1, \  - \frac{\sqrt{2}}{2}, \  0, \  \frac{\sqrt{2}}{2}, \  -1, \  \frac{\sqrt{2}}{2}, \  0, \  - \frac{\sqrt{2}}{2}\right)
\end{aligned}
$$

二つの等比級数を足すと、添字 $m$ と $8-m$ のビンに $4$ ずつ入る。添字が重なれば足す。

$$
\begin{aligned}
(\widehat v_j)_{j=0}^7\big|_{m=3} &= \left( 0, \  0, \  0, \  4, \  0, \  4, \  0, \  0\right)
\end{aligned}
$$

有限の Parseval の等式では右辺の周波数二乗和を $8$ で割る。

$$
\begin{aligned}
\sum_{n=0}^7|v_n|^2\big|_{m=3} &= 4
\end{aligned}
$$
:::

:::example[八点での余弦の DFT（4）]
番号 $n$ を $0$ から $7$ まで取り、入力 $v_n=\cos(2\pi m n/8)$ の係数を求める。最後にこの値でも計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=4
$$

標本値を番号順に並べる。

$$
\begin{aligned}
(v_n)_{n=0}^7\big|_{m=4} &= \left( 1, \  -1, \  1, \  -1, \  1, \  -1, \  1, \  -1\right)
\end{aligned}
$$

二つの等比級数を足すと、添字 $m$ と $8-m$ のビンに $4$ ずつ入る。添字が重なれば足す。

$$
\begin{aligned}
(\widehat v_j)_{j=0}^7\big|_{m=4} &= \left( 0, \  0, \  0, \  0, \  8, \  0, \  0, \  0\right)
\end{aligned}
$$

有限の Parseval の等式では右辺の周波数二乗和を $8$ で割る。

$$
\begin{aligned}
\sum_{n=0}^7|v_n|^2\big|_{m=4} &= 8
\end{aligned}
$$
:::

## 演習問題

::::exercise{level=basic}
$f(x)=e^{-2|x|}$ の Fourier 変換を求めよ。
:::solution
正側は $\int_0^\infty e^{-(2+ik)x}\dd x=(2+ik)^{-1}$、負側は $(2-ik)^{-1}$ である。和は $4/(4+k^2)$ であり、零周波数では $1$ となる。
:::
::::

::::exercise{level=basic}
$2\pi$ 周期関数 $f(x)=3+2\cos x$ の係数と平均二乗を求めよ。
:::solution
指数表示から $c_0=3,c_{\pm1}=1$ で他は零である。直交性から平均二乗は $9+1+1=11$ である。
:::
::::

::::exercise{level=basic}
$\mathcal L(t e^{-2t})(s)$ の収束半平面と値を求めよ。
:::solution
$\Re(s+2)>0$ で部分積分し、$\int_0^\infty t e^{-(s+2)t}\dd t=1/(s+2)^2$ を得る。境界 $\Re(s+2)=0$ では原始関数に増加する項があり収束しない。
:::
::::

::::exercise{level=basic}
列 $(1,0,-1,0)$ の四点 DFT を求めよ。
:::solution
和の非零項は $j=0,2$ だけなので $X_m=1-e^{-\pi im}=1-(-1)^m$ である。従って $(0,2,0,2)$ となる。
:::
::::

::::exercise{level=advanced}
$y'+2y=\mathbf1_{t\ge1}$、$y(0)=0$ を Laplace 変換で解け。
:::solution
$(s+2)Y=e^{-s}/s$ である。$1/[s(s+2)]=\tfrac12(1/s-1/(s+2))$ と分解し、逆変換して $y(t)=\frac12\mathbf1_{t\ge1}(1-e^{-2(t-1)})$ を得る。$t=1$ で解は連続で、左右の微分の跳びは強制力の跳び $1$ と一致する。
:::
::::

::::exercise{level=advanced}
$\mathcal M(e^{-3x})(s)$ を求め、$\Re s>0$ で $\mathcal M(xf')(s)=-s\mathcal M f(s)$ を確かめよ。
:::solution
$u=3x$ より $\mathcal M f=3^{-s}\Gamma(s)$ である。$xf'=-3xe^{-3x}$ なのでその変換は $-3^{-s}\Gamma(s+1)$ となる。端点で $u^s e^{-u}$ が消えるため部分積分で $\Gamma(s+1)=s\Gamma(s)$ であり、両辺が一致する。
:::
::::

::::exercise{level=advanced}
$f(x)=e^{-a x^2}$ と $g(x)=e^{-b x^2}$（$a,b>0$）の畳み込みを直接求めよ。
:::solution
$a(x-y)^2+by^2=(a+b)(y-ax/(a+b))^2+abx^2/(a+b)$ である。平行移動して Gaussian 積分により $(f*g)(x)=\sqrt{\pi/(a+b)}e^{-abx^2/(a+b)}$ を得る。$a=b=1$ のとき $x=0$ の値は $\sqrt{\pi/2}$ である。
:::
::::

::::exercise{level=advanced}
帯域を $[-\pi,\pi]$ として標本 $f(n)=\delta_{n,2}$ に対応する関数と、その $t=5/2$ での値を求めよ。
:::solution
:ref[prop-mp-sampling] で $T=1$ とすれば $f(t)=\sin(\pi(t-2))/(\pi(t-2))$ である。$t=5/2$ では分子が $1$、分母が $\pi/2$ であるから $2/\pi$ になる。
:::
::::

## 確認テスト

::::quiz{answer="2" tol="1e-9"}
$f=\mathbf1_{[-1,1]}$ について $\widehat f(0)$ を求めよ。

:::explanation
指数が $1$ になるので区間の長さを積分し、$2$ である。
:::
::::

::::quiz
畳み込みを Fourier 変換するとどうなるか。

- [ ] 二つの変換の和
- [x] 二つの変換の積
- [ ] 二つの変換の積の $1/(2\pi)$ 倍

:::explanation
$u=x-y$ と置けば指数が $e^{-iku}e^{-iky}$ に分かれ、二重積分が $\widehat f(k)\widehat g(k)$ になる。$1/(2\pi)$ はこの規約の逆変換に現れる。
:::
::::

::::quiz{answer="0.5" tol="1e-9"}
$\mathcal L(e^{-2t})(0)$ を求めよ。

:::explanation
周波数因子を代入して $\int_0^\infty e^{-2t}\dd t=1/2$ である。
:::
::::

::::quiz
Bromwich 積分の縦線をどこに取るか。

- [x] 変換の絶対収束半平面の内部
- [ ] 常に原点を通る虚軸
- [ ] 変換の極を通る縦線

:::explanation
減衰させた時間関数に Fourier 逆変換を適用するので、その減衰が時間関数を可積分にする $\gamma$ を選ぶ。極を通過すると積分の値が留数分だけ変わり得る。
:::
::::

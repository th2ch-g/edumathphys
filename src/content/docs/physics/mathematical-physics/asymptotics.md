---
title: "漸近展開と級数"
description: "漸近級数の誤差を定義し、積分の局所化、和の端点補正、発散級数の再総和を導く。"
notation:
  - symbol: 'f(x)\sim\sum a_n\phi_n(x)'
    meaning: "指定極限での漸近展開"
  - symbol: 'B_m'
    meaning: "Bernoulli 数"
  - symbol: '\mathcal B f'
    meaning: "形式級数の Borel 変換"
---

## この章の目標

- 収束と漸近性を剰余から判定する。
- 相の最小点、端点、停留点の近くで積分を近似する。
- Bernoulli 数から和と積分の差を評価する。
- 発散する展開から Borel 和を計算する。

## 漸近級数の定義と性質

### 有限個の項で誤差を制御する

大きなパラメータで小さい項だけを残すときは、無限和が収束するかより、打ち切った剰余がどれだけ小さいかが重要である。

:::definition[漸近展開]{#def-mp-asymptotic-scale}
$x\to x_0$ の指定した経路で、零でない関数列 $\phi_{n+1}=o(\phi_n)$ を漸近尺度という。任意の固定した $N$ に対し
$$
f(x)-\sum_{n=0}^N a_n\phi_n(x)=o(\phi_N(x))
$$
なら $f\sim\sum_{n\ge0}a_n\phi_n$ と書き、:term[漸近展開]{en="asymptotic expansion" yomi="ぜんきんてんかい"} と呼ぶ。ここで $h=o(\phi)$ は $h/\phi\to0$ を意味する。
:::

$a_0$ は $f/\phi_0$ の極限であり、既知の項を順に引いて $\phi_n$ で割れば $a_n$ が一意に決まる。二つの関数が同じ展開を持っても差は零とは限らない。例えば $e^{-x}=o(x^{-N})$ がすべての $N$ について成り立つので、$1$ と $1+e^{-x}$ は $x\to+\infty$ で同じ逆冪展開を持つ。逆冪の展開では剰余を $O(x^{-N-1})$ としても、次の係数まで展開することで定義と同値になる。

逆冪の展開は加法と積に関して有限打ち切りで計算できる。各因子を $N$ 次まで置き換えると、誤差と有界な因子の積が $O(x^{-N-1})$ のためである。定数項が零でないときは逆数も、積の係数が $1,0,\ldots,0$ となるように三角的に決められる。$f$ は大きい $x$ で零から離れているため、$fP_N-1=O(x^{-N-1})$ を $f$ で割れば逆数の誤差評価になる。一方、微分には別の評価が要る。$x^{-1}$ の全ての冪より小さい $e^{-x}\sin(e^{2x})$ の導関数は、$e^{2x}=2\pi n$ となる列に沿って $2e^x$ になる。

:::example[積分の有限打ち切り]
$I(x)=\int_0^\infty e^{-t}/(1+t/x)\dd t$ では有限等比和の剰余が $(-t/x)^{N+1}/(1+t/x)$ なので
$I(x)=\sum_{n=0}^N(-1)^n n!x^{-n}+R_N$、$|R_N|\le (N+1)!x^{-N-1}$ となる。階乗は $\int_0^\infty t^ne^{-t}\dd t$ の逐次部分積分から出る。$x=10$ で $N=2$ とすると $0.92$ に誤差 $0.006$ 以下である。
:::

### 幾何展開と一様性

小さい変数を $t$ とすると、幾何恒等式 $1/(1+mt)=\sum_{j=0}^{N-1}(-mt)^j+(-mt)^N/(1+mt)$ は有限和として厳密である。従って $t\to0^+$ で $m$ を有界に保てば、剰余は一様に $O(t^N)$ である。ただし $m$ も変わる問題で $mt$ が一定なら、同じ次数評価は使えない。実際、分母を展開する小さい量は $t$ 単独ではなく $mt$ である。

:::example[幾何級数の打ち切りと誤差（1）]
正の $t$ で有限多項式を作る。余りの符号から真値を挟む。基準の値で計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

展開する関数を決める。

$$
\begin{aligned}
f(t)\big|_{m=1} &= \frac{1}{t + 1}
\end{aligned}
$$

二次までの部分和を取る。

$$
\begin{aligned}
S_2(t)\big|_{m=1} &= t^{2} - t + 1
\end{aligned}
$$

差を一つの分数に戻す。

$$
\begin{aligned}
f-S_2\big|_{m=1} &= - \frac{t^{3}}{t + 1}
\end{aligned}
$$

差を $t^3$ で割ると、$t=0$ における最初の捨てた係数まで得られる。

$$
\begin{aligned}
(f-S_2)/t^3\big|_{m=1} &= - \frac{1}{t + 1}
\end{aligned}
$$

観測点 $t$=$\frac{1}{10}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
f(t)\big|_{t=\frac{1}{10}} &= \frac{10}{11} \\
\frac{\partial^1 (f(t))}{\partial t^1} &= - \frac{1}{\left(t + 1\right)^{2}} \\
\left.\frac{\partial^1 (f(t))}{\partial t^1}\right|_{t=\frac{1}{10}} &= - \frac{100}{121}
\end{aligned}
$$

観測点 $t$=$\frac{1}{10}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
S_2(t)\big|_{t=\frac{1}{10}} &= \frac{91}{100} \\
\frac{\partial^1 (S_2(t))}{\partial t^1} &= 2 t - 1 \\
\left.\frac{\partial^1 (S_2(t))}{\partial t^1}\right|_{t=\frac{1}{10}} &= - \frac{4}{5}
\end{aligned}
$$

観測点 $t$=$\frac{1}{10}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
f-S_2\big|_{t=\frac{1}{10}} &= - \frac{1}{1100} \\
\frac{\partial^1 (f-S_2)}{\partial t^1} &= \frac{t^{2} \left(- 2 t - 3\right)}{\left(t + 1\right)^{2}} \\
\left.\frac{\partial^1 (f-S_2)}{\partial t^1}\right|_{t=\frac{1}{10}} &= - \frac{16}{605}
\end{aligned}
$$

観測点 $t$=$\frac{1}{10}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
(f-S_2)/t^3\big|_{t=\frac{1}{10}} &= - \frac{10}{11} \\
\frac{\partial^1 ((f-S_2)/t^3)}{\partial t^1} &= \frac{1}{\left(t + 1\right)^{2}} \\
\left.\frac{\partial^1 ((f-S_2)/t^3)}{\partial t^1}\right|_{t=\frac{1}{10}} &= \frac{100}{121}
\end{aligned}
$$
:::

:::example[幾何級数の打ち切りと誤差（2）]
正の $t$ で有限多項式を作る。余りの符号から真値を挟む。値を変えて依存性を調べる。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

展開する関数を決める。

$$
\begin{aligned}
f(t)\big|_{m=2} &= \frac{1}{2 t + 1}
\end{aligned}
$$

二次までの部分和を取る。

$$
\begin{aligned}
S_2(t)\big|_{m=2} &= 4 t^{2} - 2 t + 1
\end{aligned}
$$

差を一つの分数に戻す。

$$
\begin{aligned}
f-S_2\big|_{m=2} &= - \frac{8 t^{3}}{2 t + 1}
\end{aligned}
$$

差を $t^3$ で割ると、$t=0$ における最初の捨てた係数まで得られる。

$$
\begin{aligned}
(f-S_2)/t^3\big|_{m=2} &= - \frac{8}{2 t + 1}
\end{aligned}
$$

観測点 $t$=$\frac{1}{10}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
f(t)\big|_{t=\frac{1}{10}} &= \frac{5}{6} \\
\frac{\partial^1 (f(t))}{\partial t^1} &= - \frac{2}{\left(2 t + 1\right)^{2}} \\
\left.\frac{\partial^1 (f(t))}{\partial t^1}\right|_{t=\frac{1}{10}} &= - \frac{25}{18}
\end{aligned}
$$

観測点 $t$=$\frac{1}{10}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
S_2(t)\big|_{t=\frac{1}{10}} &= \frac{21}{25} \\
\frac{\partial^1 (S_2(t))}{\partial t^1} &= 8 t - 2 \\
\left.\frac{\partial^1 (S_2(t))}{\partial t^1}\right|_{t=\frac{1}{10}} &= - \frac{6}{5}
\end{aligned}
$$

観測点 $t$=$\frac{1}{10}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
f-S_2\big|_{t=\frac{1}{10}} &= - \frac{1}{150} \\
\frac{\partial^1 (f-S_2)}{\partial t^1} &= \frac{t^{2} \left(- 32 t - 24\right)}{\left(2 t + 1\right)^{2}} \\
\left.\frac{\partial^1 (f-S_2)}{\partial t^1}\right|_{t=\frac{1}{10}} &= - \frac{17}{90}
\end{aligned}
$$

観測点 $t$=$\frac{1}{10}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
(f-S_2)/t^3\big|_{t=\frac{1}{10}} &= - \frac{20}{3} \\
\frac{\partial^1 ((f-S_2)/t^3)}{\partial t^1} &= \frac{16}{\left(2 t + 1\right)^{2}} \\
\left.\frac{\partial^1 ((f-S_2)/t^3)}{\partial t^1}\right|_{t=\frac{1}{10}} &= \frac{100}{9}
\end{aligned}
$$
:::

:::example[幾何級数の打ち切りと誤差（3）]
正の $t$ で有限多項式を作る。余りの符号から真値を挟む。さらに値を変えて比較する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

展開する関数を決める。

$$
\begin{aligned}
f(t)\big|_{m=3} &= \frac{1}{3 t + 1}
\end{aligned}
$$

二次までの部分和を取る。

$$
\begin{aligned}
S_2(t)\big|_{m=3} &= 9 t^{2} - 3 t + 1
\end{aligned}
$$

差を一つの分数に戻す。

$$
\begin{aligned}
f-S_2\big|_{m=3} &= - \frac{27 t^{3}}{3 t + 1}
\end{aligned}
$$

差を $t^3$ で割ると、$t=0$ における最初の捨てた係数まで得られる。

$$
\begin{aligned}
(f-S_2)/t^3\big|_{m=3} &= - \frac{27}{3 t + 1}
\end{aligned}
$$

観測点 $t$=$\frac{1}{10}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
f(t)\big|_{t=\frac{1}{10}} &= \frac{10}{13} \\
\frac{\partial^1 (f(t))}{\partial t^1} &= - \frac{3}{\left(3 t + 1\right)^{2}} \\
\left.\frac{\partial^1 (f(t))}{\partial t^1}\right|_{t=\frac{1}{10}} &= - \frac{300}{169}
\end{aligned}
$$

観測点 $t$=$\frac{1}{10}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
S_2(t)\big|_{t=\frac{1}{10}} &= \frac{79}{100} \\
\frac{\partial^1 (S_2(t))}{\partial t^1} &= 18 t - 3 \\
\left.\frac{\partial^1 (S_2(t))}{\partial t^1}\right|_{t=\frac{1}{10}} &= - \frac{6}{5}
\end{aligned}
$$

観測点 $t$=$\frac{1}{10}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
f-S_2\big|_{t=\frac{1}{10}} &= - \frac{27}{1300} \\
\frac{\partial^1 (f-S_2)}{\partial t^1} &= \frac{t^{2} \left(- 162 t - 81\right)}{\left(3 t + 1\right)^{2}} \\
\left.\frac{\partial^1 (f-S_2)}{\partial t^1}\right|_{t=\frac{1}{10}} &= - \frac{486}{845}
\end{aligned}
$$

観測点 $t$=$\frac{1}{10}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
(f-S_2)/t^3\big|_{t=\frac{1}{10}} &= - \frac{270}{13} \\
\frac{\partial^1 ((f-S_2)/t^3)}{\partial t^1} &= \frac{81}{\left(3 t + 1\right)^{2}} \\
\left.\frac{\partial^1 ((f-S_2)/t^3)}{\partial t^1}\right|_{t=\frac{1}{10}} &= \frac{8100}{169}
\end{aligned}
$$
:::

:::example[幾何級数の打ち切りと誤差（4）]
正の $t$ で有限多項式を作る。余りの符号から真値を挟む。最後にこの値でも計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=4
$$

展開する関数を決める。

$$
\begin{aligned}
f(t)\big|_{m=4} &= \frac{1}{4 t + 1}
\end{aligned}
$$

二次までの部分和を取る。

$$
\begin{aligned}
S_2(t)\big|_{m=4} &= 16 t^{2} - 4 t + 1
\end{aligned}
$$

差を一つの分数に戻す。

$$
\begin{aligned}
f-S_2\big|_{m=4} &= - \frac{64 t^{3}}{4 t + 1}
\end{aligned}
$$

差を $t^3$ で割ると、$t=0$ における最初の捨てた係数まで得られる。

$$
\begin{aligned}
(f-S_2)/t^3\big|_{m=4} &= - \frac{64}{4 t + 1}
\end{aligned}
$$

観測点 $t$=$\frac{1}{10}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
f(t)\big|_{t=\frac{1}{10}} &= \frac{5}{7} \\
\frac{\partial^1 (f(t))}{\partial t^1} &= - \frac{4}{\left(4 t + 1\right)^{2}} \\
\left.\frac{\partial^1 (f(t))}{\partial t^1}\right|_{t=\frac{1}{10}} &= - \frac{100}{49}
\end{aligned}
$$

観測点 $t$=$\frac{1}{10}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
S_2(t)\big|_{t=\frac{1}{10}} &= \frac{19}{25} \\
\frac{\partial^1 (S_2(t))}{\partial t^1} &= 32 t - 4 \\
\left.\frac{\partial^1 (S_2(t))}{\partial t^1}\right|_{t=\frac{1}{10}} &= - \frac{4}{5}
\end{aligned}
$$

観測点 $t$=$\frac{1}{10}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
f-S_2\big|_{t=\frac{1}{10}} &= - \frac{8}{175} \\
\frac{\partial^1 (f-S_2)}{\partial t^1} &= \frac{t^{2} \left(- 512 t - 192\right)}{\left(4 t + 1\right)^{2}} \\
\left.\frac{\partial^1 (f-S_2)}{\partial t^1}\right|_{t=\frac{1}{10}} &= - \frac{304}{245}
\end{aligned}
$$

観測点 $t$=$\frac{1}{10}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
(f-S_2)/t^3\big|_{t=\frac{1}{10}} &= - \frac{320}{7} \\
\frac{\partial^1 ((f-S_2)/t^3)}{\partial t^1} &= \frac{256}{\left(4 t + 1\right)^{2}} \\
\left.\frac{\partial^1 ((f-S_2)/t^3)}{\partial t^1}\right|_{t=\frac{1}{10}} &= \frac{6400}{49}
\end{aligned}
$$
:::

## Laplace の方法と Stirling の公式

### 内部の極小点を拡大する

:::theorem[Laplace の方法]{#thm-mp-laplace-method}
$[a,b]$ 上で滑らかな実関数 $\phi$ が唯一の最小点 $t_0\in(a,b)$ を持ち、$\phi''(t_0)>0$ とする。滑らかな $f$ に対して $I(\lambda)=\int_a^b f(t)e^{-\lambda\phi(t)}\dd t$ は $\lambda\to+\infty$ で
$$
I(\lambda)=e^{-\lambda\phi(t_0)}
\sqrt{\frac{2\pi}{\lambda\phi''(t_0)}}\bigl(f(t_0)+O(\lambda^{-1})\bigr)
\label{eq-mp-laplace-leading}
$$
となる。
:::

:::proof
最小点を含む小区間の外では $\phi\ge\phi(t_0)+\eta$ となる。中では $u=\sgn(t-t_0)\sqrt{2(\phi(t)-\phi(t_0))}$ を取る。$\phi(t)-\phi(t_0)=(t-t_0)^2 b(t)$、$b(t_0)=\phi''(t_0)/2>0$ なので滑らかな座標であり、$t'(0)=1/\sqrt{\phi''(t_0)}$ である。$u$ で対称な区間を取り、$F(u)=f(t(u))t'(u)$ を Taylor 展開する。奇数冪は対称積分で消え、偶数冪の積分を実線全体に延ばすと誤差は指数的に小さい。$\int e^{-\lambda u^2/2}\dd u=\sqrt{2\pi/\lambda}$ は二乗して平面で極座標を使えば得られる。二回ごとのモーメントは部分積分で $\int u^{2m}e^{-\lambda u^2/2}\dd u=(2m-1)!!\sqrt{2\pi}\lambda^{-m-1/2}$ となる。Taylor 剰余を絶対値で積分すれば、任意の固定次数まで誤差を評価できる。
:::

:::example[Stirling の公式]
$n!=\int_0^\infty t^ne^{-t}\dd t$ で $t=ns$ と置くと $n^{n+1}e^{-n}\int_0^\infty e^{-n(s-1-\log s)}\dd s$ となる。$\phi=s-1-\log s$ は $s=1$ で最小、$\phi''(1)=1$ である。中心近くを除く部分では $\phi\ge\eta$ かつ $\int_0^\infty e^{-\phi}\dd s<\infty$ なので同じ証明が使える。$h=s-1$ と置けば $2\phi=h^2-\tfrac23 h^3+\tfrac12 h^4+O(h^5)$ で、$u=\sgn(h)\sqrt{2\phi}$ を逆に解くと $h=u+u^2/3+u^3/36+O(u^4)$ である。従って $s'(u)=1+2u/3+u^2/12+O(u^3)$ となる。少なくとも四次まで展開して奇数項を落とせば
$$
n!=\sqrt{2\pi n}(n/e)^n\left(1+\frac1{12n}+O(n^{-2})\right)
$$
を得る。$n=10$ では相対補正 $1/120$ は約 $0.0083$ である。
:::

### 多変数の極小点と退化の影響

相 $\phi$ が $d$ 次元の積分領域の内部で唯一の厳密な最小値を取り、そこで Hessian $A$ が正定値であるとする。振幅 $f$ の台がコンパクトで、最小点 $x_0$ を含むとする。主項は $e^{-\lambda\phi(x_0)}(2\pi/\lambda)^{d/2}f(x_0)/\sqrt{\det A}$ である。実際、近傍の外は指数的に小さく、近傍で $\phi-\phi(x_0)\ge c|x-x_0|^2$ が成り立つ。$x=x_0+\lambda^{-1/2}u$ と置けば点ごとに相の差が $u^\top Au/2$ に収束し、$e^{-c|u|^2}$ が優関数になる。対称行列 $A$ を回転で対角化して Gaussian 積分を各方向に分けると行列式が出る。

例えば $\phi=2x^2+2xy+2y^2$ なら Hessian は $\begin{pmatrix}4&2\\2&4\end{pmatrix}$、行列式は $12$ であり、全平面の振幅 $1$ の積分は厳密に $\pi/(\sqrt3\lambda)$ である。二次係数の一つが零のときはこの式は使えない。例えば $\phi=x^4$ の実線積分は $u=\lambda^{1/4}x$ で $\lambda^{-1/4}\Gamma(1/4)/2$ になり、幅の冪が変わる。これは二次近似の行列式を零に近づけた式からは得られない。

:::example[極小点が積分の端へ近づく場合]
$I(\lambda,a)=\int_0^\infty e^{-\lambda(t-a)^2/2}\dd t$ を考える。$a>0$ が固定なら極小点は領域の内部にあるが、$a=0$ では端点に重なる。$u=\sqrt{\lambda}(t-a)$ と置くと、どちらの場合にも使える厳密な式
$$
I(\lambda,a)=\lambda^{-1/2}\int_{-a\sqrt\lambda}^\infty e^{-u^2/2}\dd u
$$
を得る。固定した $a>0$ では下限が $-\infty$ に近づいて $I\sim\sqrt{2\pi/\lambda}$ となり、$a=0$ では対称性でその半分となる。$a=A/\sqrt\lambda$ なら下限は $-A$ のままであり、$I$ の係数は $A$ に応じて連続的に変化する。従って $a$ が零に近づく範囲を一度に扱うには、下限を残した式を使う。

一方、固定した $a<0$ では相の最小値が端点にあり、その点で一階導関数が $|a|$ である。$t=v/(\lambda|a|)$ と置いた積分では、$e^{-\lambda a^2/2}/(\lambda|a|)$ をくくると被積分関数が $e^{-v}e^{-v^2/(2\lambda a^2)}$ になる。これを $e^{-v}$ で抑えて優収束を使えば、残りの積分は $1$ へ収束する。端点からの距離が $\lambda^{-1/2}$ 程度になると、この端点近似と内部極小の近似の両方が適用範囲を外れることが分かる。
:::

### Gaussian のモーメントと補正項

内部の二次極小で $s=\sqrt{\lambda}(x-x_0)$ と置くと、振幅の展開の奇数項は対称積分で零になる。$\int e^{-ms^2}s^{2j}\dd s$ は、$m$ で微分するか部分積分して前の次数を掛け直せば求まる。切り落とした領域の寄与を別に評価すれば、偶数項ごとの補正が $\lambda^{-1}$ 刻みで出る。次の模型では相が厳密に二次なので、振幅が多項式なら有限個の項で終わる。

:::example[Gaussian 積分の高次振幅（1）]
幅が $\lambda^{-1/2}$ になる積分で、振幅の零点と次数が主項にどう影響するか計算する。基準の値で計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

振幅を $t^{2m}$ とし、相を固定する。

$$
\begin{aligned}
f(t)\big|_{m=1} &= t^{2} e^{- t^{2} u}
\end{aligned}
$$

正側に移して $v=ut^2$ と置くと、Jacobian の因子も含めて Gamma 積分になる。

$$
\begin{aligned}
M_m(u)\big|_{m=1} &= \frac{\sqrt{\pi}}{2 u^{\frac{3}{2}}}
\end{aligned}
$$

前の偶数モーメントと比を取ると、振幅に $t^2$ を加えるたびに一段減衰が速くなる。

$$
\begin{aligned}
M_m/M_{m-1}\big|_{m=1} &= \frac{1}{2 u}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
f(t)\big|_{u=2} &= t^{2} e^{- 2 t^{2}} \\
\frac{\partial^1 (f(t))}{\partial u^1} &= - t^{4} e^{- t^{2} u} \\
\left.\frac{\partial^1 (f(t))}{\partial u^1}\right|_{u=2} &= - t^{4} e^{- 2 t^{2}}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
M_m(u)\big|_{u=2} &= \frac{\sqrt{2} \sqrt{\pi}}{8} \\
\frac{\partial^1 (M_m(u))}{\partial u^1} &= - \frac{3 \sqrt{\pi}}{4 u^{\frac{5}{2}}} \\
\left.\frac{\partial^1 (M_m(u))}{\partial u^1}\right|_{u=2} &= - \frac{3 \sqrt{2} \sqrt{\pi}}{32}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
M_m/M_{m-1}\big|_{u=2} &= \frac{1}{4} \\
\frac{\partial^1 (M_m/M_{m-1})}{\partial u^1} &= - \frac{1}{2 u^{2}} \\
\left.\frac{\partial^1 (M_m/M_{m-1})}{\partial u^1}\right|_{u=2} &= - \frac{1}{8}
\end{aligned}
$$
:::

:::example[Gaussian 積分の高次振幅（2）]
幅が $\lambda^{-1/2}$ になる積分で、振幅の零点と次数が主項にどう影響するか計算する。値を変えて依存性を調べる。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

振幅を $t^{2m}$ とし、相を固定する。

$$
\begin{aligned}
f(t)\big|_{m=2} &= t^{4} e^{- t^{2} u}
\end{aligned}
$$

正側に移して $v=ut^2$ と置くと、Jacobian の因子も含めて Gamma 積分になる。

$$
\begin{aligned}
M_m(u)\big|_{m=2} &= \frac{3 \sqrt{\pi}}{4 u^{\frac{5}{2}}}
\end{aligned}
$$

前の偶数モーメントと比を取ると、振幅に $t^2$ を加えるたびに一段減衰が速くなる。

$$
\begin{aligned}
M_m/M_{m-1}\big|_{m=2} &= \frac{3}{2 u}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
f(t)\big|_{u=2} &= t^{4} e^{- 2 t^{2}} \\
\frac{\partial^1 (f(t))}{\partial u^1} &= - t^{6} e^{- t^{2} u} \\
\left.\frac{\partial^1 (f(t))}{\partial u^1}\right|_{u=2} &= - t^{6} e^{- 2 t^{2}}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
M_m(u)\big|_{u=2} &= \frac{3 \sqrt{2} \sqrt{\pi}}{32} \\
\frac{\partial^1 (M_m(u))}{\partial u^1} &= - \frac{15 \sqrt{\pi}}{8 u^{\frac{7}{2}}} \\
\left.\frac{\partial^1 (M_m(u))}{\partial u^1}\right|_{u=2} &= - \frac{15 \sqrt{2} \sqrt{\pi}}{128}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
M_m/M_{m-1}\big|_{u=2} &= \frac{3}{4} \\
\frac{\partial^1 (M_m/M_{m-1})}{\partial u^1} &= - \frac{3}{2 u^{2}} \\
\left.\frac{\partial^1 (M_m/M_{m-1})}{\partial u^1}\right|_{u=2} &= - \frac{3}{8}
\end{aligned}
$$
:::

:::example[Gaussian 積分の高次振幅（3）]
幅が $\lambda^{-1/2}$ になる積分で、振幅の零点と次数が主項にどう影響するか計算する。さらに値を変えて比較する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

振幅を $t^{2m}$ とし、相を固定する。

$$
\begin{aligned}
f(t)\big|_{m=3} &= t^{6} e^{- t^{2} u}
\end{aligned}
$$

正側に移して $v=ut^2$ と置くと、Jacobian の因子も含めて Gamma 積分になる。

$$
\begin{aligned}
M_m(u)\big|_{m=3} &= \frac{15 \sqrt{\pi}}{8 u^{\frac{7}{2}}}
\end{aligned}
$$

前の偶数モーメントと比を取ると、振幅に $t^2$ を加えるたびに一段減衰が速くなる。

$$
\begin{aligned}
M_m/M_{m-1}\big|_{m=3} &= \frac{5}{2 u}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
f(t)\big|_{u=2} &= t^{6} e^{- 2 t^{2}} \\
\frac{\partial^1 (f(t))}{\partial u^1} &= - t^{8} e^{- t^{2} u} \\
\left.\frac{\partial^1 (f(t))}{\partial u^1}\right|_{u=2} &= - t^{8} e^{- 2 t^{2}}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
M_m(u)\big|_{u=2} &= \frac{15 \sqrt{2} \sqrt{\pi}}{128} \\
\frac{\partial^1 (M_m(u))}{\partial u^1} &= - \frac{105 \sqrt{\pi}}{16 u^{\frac{9}{2}}} \\
\left.\frac{\partial^1 (M_m(u))}{\partial u^1}\right|_{u=2} &= - \frac{105 \sqrt{2} \sqrt{\pi}}{512}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
M_m/M_{m-1}\big|_{u=2} &= \frac{5}{4} \\
\frac{\partial^1 (M_m/M_{m-1})}{\partial u^1} &= - \frac{5}{2 u^{2}} \\
\left.\frac{\partial^1 (M_m/M_{m-1})}{\partial u^1}\right|_{u=2} &= - \frac{5}{8}
\end{aligned}
$$
:::

:::example[Gaussian 積分の高次振幅（4）]
幅が $\lambda^{-1/2}$ になる積分で、振幅の零点と次数が主項にどう影響するか計算する。最後にこの値でも計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=4
$$

振幅を $t^{2m}$ とし、相を固定する。

$$
\begin{aligned}
f(t)\big|_{m=4} &= t^{8} e^{- t^{2} u}
\end{aligned}
$$

正側に移して $v=ut^2$ と置くと、Jacobian の因子も含めて Gamma 積分になる。

$$
\begin{aligned}
M_m(u)\big|_{m=4} &= \frac{105 \sqrt{\pi}}{16 u^{\frac{9}{2}}}
\end{aligned}
$$

前の偶数モーメントと比を取ると、振幅に $t^2$ を加えるたびに一段減衰が速くなる。

$$
\begin{aligned}
M_m/M_{m-1}\big|_{m=4} &= \frac{7}{2 u}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
f(t)\big|_{u=2} &= t^{8} e^{- 2 t^{2}} \\
\frac{\partial^1 (f(t))}{\partial u^1} &= - t^{10} e^{- t^{2} u} \\
\left.\frac{\partial^1 (f(t))}{\partial u^1}\right|_{u=2} &= - t^{10} e^{- 2 t^{2}}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
M_m(u)\big|_{u=2} &= \frac{105 \sqrt{2} \sqrt{\pi}}{512} \\
\frac{\partial^1 (M_m(u))}{\partial u^1} &= - \frac{945 \sqrt{\pi}}{32 u^{\frac{11}{2}}} \\
\left.\frac{\partial^1 (M_m(u))}{\partial u^1}\right|_{u=2} &= - \frac{945 \sqrt{2} \sqrt{\pi}}{2048}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
M_m/M_{m-1}\big|_{u=2} &= \frac{7}{4} \\
\frac{\partial^1 (M_m/M_{m-1})}{\partial u^1} &= - \frac{7}{2 u^{2}} \\
\left.\frac{\partial^1 (M_m/M_{m-1})}{\partial u^1}\right|_{u=2} &= - \frac{7}{8}
\end{aligned}
$$
:::

## Watson の補題

### 端点での冪から減衰の冪を決める

:::theorem[Watson の補題]{#thm-mp-watson}
$\alpha>0$、$N\ge1$ とし、$f(t)=t^{\alpha-1}\sum_{n=0}^{N-1}a_nt^n+O(t^{\alpha+N-1})$ が $t\downarrow0$ で成り立つとする。$f$ は $(0,\infty)$ で局所可積分とし、ある実数 $c$ で $\int_1^\infty |f(t)|e^{-ct}\dd t<\infty$ とする。このとき $x\to+\infty$ で
$$
\int_0^\infty e^{-xt}f(t)\dd t=\sum_{n=0}^{N-1}a_n\Gamma(\alpha+n)x^{-\alpha-n}+O(x^{-\alpha-N})
$$
である。ここで $\Gamma(p)=\int_0^\infty e^{-t}t^{p-1}\dd t$ とする。
:::

:::proof
$0<\epsilon<1$ で積分を分ける。$[\epsilon,\infty)$ では $e^{-xt}\le e^{-(x-c)\epsilon}e^{-ct}$ を用いる。$[0,\epsilon]$ では各単項式の積分を無限遠まで延ばして $u=xt$ と置く。切り落とした尾は指数的に小さい。剰余の絶対値は $C\int_0^\infty e^{-xt}t^{\alpha+N-1}\dd t=C\Gamma(\alpha+N)x^{-\alpha-N}$ で抑えられる。
:::

:::example[端点に平方根の特異性がある積分]
$f(t)=t^{-1/2}e^{-2t}$ では $a_0=1,a_1=-2,a_2=2$ である。$\Gamma(1/2)=\sqrt\pi$ は $t=u^2$ と Gaussian 積分から分かる。部分積分の漸化式から積分は $\sqrt\pi x^{-1/2}(1-x^{-1}+\tfrac32 x^{-2}+O(x^{-3}))$ で、厳密値 $\sqrt\pi/\sqrt{x+2}$ の二項展開と一致する。
:::

:::example[正則な端点と指数変数]
$I(x)=\int_0^\infty e^{-x(e^t-1)}\dd t$ で $u=e^t-1$ と置くと、振幅は $(1+u)^{-1}$ になる。端点で $1-u+u^2$ と展開すれば $I=x^{-1}-x^{-2}+2x^{-3}+O(x^{-4})$ である。元の変数の指数を直接高次まで展開する場合と違い、この変数選択では誤差は $(1+u)^{-1}$ の幾何級数の剰余から $6x^{-4}$ 以下とただちに分かる。
:::

### 端点の変数を変えて因子を分離する

$\phi(0)=0$、$\phi^\prime(0)>0$ の端点なら、$t=\phi(x)$ の局所逆関数が存在する。$x^\alpha g(x)\dd x$ を移すと $x(t)^\alpha g(x(t))x^\prime(t)\dd t$ であり、$x(t)/t$ が正に延びるので再び $t^\alpha$ と滑らかな振幅に分けられる。この議論で現れる導関数は、積分の向きを変えたときは符号も変わる。$\phi^\prime(0)=0$ なら逆関数の冪が変わるため、Watson の補題へ入れる指数も変化する。

端点で $f(t)=\sum_{j=0}^{J-1}a_jt^{\beta_j-1}+O(t^{\beta_J-1})$、$0<\beta_0<\cdots<\beta_J$ と仮定しても、各冪を同じ方法で積分できる。複素パラメータ $z$ を使うときは、固定した $\eta>0$ に対し $|\arg z|\le\pi/2-\eta$ の扇形で $|z|\to\infty$ とする。主値の対数で $z^{-\beta_j}$ を定めれば、$\int_0^\infty e^{-zt}t^{\beta_j-1}\dd t=\Gamma(\beta_j)z^{-\beta_j}$ である。両辺は $\Re z>0$ で正則で、正実軸上では変数変換で一致するため、一致の定理が使える。剰余には $|e^{-zt}|\le e^{-|z|t\sin\eta}$ を適用する。近傍の剰余が $O(|z|^{-\beta_J})$ となり、近傍の外の積分は指数的に小さい。この絶対値による評価は $\eta$ を固定して使う。

:::example[退化した端点（1）]
半直線で $\phi(x)=x^m$ とし、$u=x^m$ を入れる。基準の値で計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

逆変換を記し、Jacobian のべきが変わることを確認する。

$$
\begin{aligned}
x(u)\big|_{m=1} &= u
\end{aligned}
$$

振幅が $x$ のときは、逆変換の値と導関数の積になる。

$$
\begin{aligned}
x(u)x^\prime(u)\big|_{m=1} &= u
\end{aligned}
$$

この冪に Laplace の指数重みを掛ける。$m$ が大きいほど減衰が遅くなる。

$$
\begin{aligned}
\int_0^\infty e^{-s x^m}x\,\dd x\big|_{m=1} &= \frac{1}{s^{2}}
\end{aligned}
$$

観測点 $s$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\int_0^\infty e^{-s x^m}x\,\dd x\big|_{s=2} &= \frac{1}{4} \\
\frac{\partial^1 (\int_0^\infty e^{-s x^m}x\,\dd x)}{\partial s^1} &= - \frac{2}{s^{3}} \\
\left.\frac{\partial^1 (\int_0^\infty e^{-s x^m}x\,\dd x)}{\partial s^1}\right|_{s=2} &= - \frac{1}{4}
\end{aligned}
$$
:::

:::example[退化した端点（2）]
半直線で $\phi(x)=x^m$ とし、$u=x^m$ を入れる。値を変えて依存性を調べる。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

逆変換を記し、Jacobian のべきが変わることを確認する。

$$
\begin{aligned}
x(u)\big|_{m=2} &= \sqrt{u}
\end{aligned}
$$

振幅が $x$ のときは、逆変換の値と導関数の積になる。

$$
\begin{aligned}
x(u)x^\prime(u)\big|_{m=2} &= \frac{1}{2}
\end{aligned}
$$

この冪に Laplace の指数重みを掛ける。$m$ が大きいほど減衰が遅くなる。

$$
\begin{aligned}
\int_0^\infty e^{-s x^m}x\,\dd x\big|_{m=2} &= \frac{1}{2 s}
\end{aligned}
$$

観測点 $s$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\int_0^\infty e^{-s x^m}x\,\dd x\big|_{s=2} &= \frac{1}{4} \\
\frac{\partial^1 (\int_0^\infty e^{-s x^m}x\,\dd x)}{\partial s^1} &= - \frac{1}{2 s^{2}} \\
\left.\frac{\partial^1 (\int_0^\infty e^{-s x^m}x\,\dd x)}{\partial s^1}\right|_{s=2} &= - \frac{1}{8}
\end{aligned}
$$
:::

:::example[退化した端点（3）]
半直線で $\phi(x)=x^m$ とし、$u=x^m$ を入れる。さらに値を変えて比較する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

逆変換を記し、Jacobian のべきが変わることを確認する。

$$
\begin{aligned}
x(u)\big|_{m=3} &= \sqrt[3]{u}
\end{aligned}
$$

振幅が $x$ のときは、逆変換の値と導関数の積になる。

$$
\begin{aligned}
x(u)x^\prime(u)\big|_{m=3} &= \frac{1}{3 \sqrt[3]{u}}
\end{aligned}
$$

この冪に Laplace の指数重みを掛ける。$m$ が大きいほど減衰が遅くなる。

$$
\begin{aligned}
\int_0^\infty e^{-s x^m}x\,\dd x\big|_{m=3} &= \frac{\Gamma\left(\frac{2}{3}\right)}{3 s^{\frac{2}{3}}}
\end{aligned}
$$

観測点 $s$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\int_0^\infty e^{-s x^m}x\,\dd x\big|_{s=2} &= \frac{\sqrt[3]{2} \Gamma\left(\frac{2}{3}\right)}{6} \\
\frac{\partial^1 (\int_0^\infty e^{-s x^m}x\,\dd x)}{\partial s^1} &= - \frac{2 \Gamma\left(\frac{2}{3}\right)}{9 s^{\frac{5}{3}}} \\
\left.\frac{\partial^1 (\int_0^\infty e^{-s x^m}x\,\dd x)}{\partial s^1}\right|_{s=2} &= - \frac{\sqrt[3]{2} \Gamma\left(\frac{2}{3}\right)}{18}
\end{aligned}
$$
:::

:::example[退化した端点（4）]
半直線で $\phi(x)=x^m$ とし、$u=x^m$ を入れる。最後にこの値でも計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=4
$$

逆変換を記し、Jacobian のべきが変わることを確認する。

$$
\begin{aligned}
x(u)\big|_{m=4} &= \sqrt[4]{u}
\end{aligned}
$$

振幅が $x$ のときは、逆変換の値と導関数の積になる。

$$
\begin{aligned}
x(u)x^\prime(u)\big|_{m=4} &= \frac{1}{4 \sqrt{u}}
\end{aligned}
$$

この冪に Laplace の指数重みを掛ける。$m$ が大きいほど減衰が遅くなる。

$$
\begin{aligned}
\int_0^\infty e^{-s x^m}x\,\dd x\big|_{m=4} &= \frac{\sqrt{\pi}}{4 \sqrt{s}}
\end{aligned}
$$

観測点 $s$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\int_0^\infty e^{-s x^m}x\,\dd x\big|_{s=2} &= \frac{\sqrt{2} \sqrt{\pi}}{8} \\
\frac{\partial^1 (\int_0^\infty e^{-s x^m}x\,\dd x)}{\partial s^1} &= - \frac{\sqrt{\pi}}{8 s^{\frac{3}{2}}} \\
\left.\frac{\partial^1 (\int_0^\infty e^{-s x^m}x\,\dd x)}{\partial s^1}\right|_{s=2} &= - \frac{\sqrt{2} \sqrt{\pi}}{32}
\end{aligned}
$$
:::

## 停留位相法

### 相の変化が遅い点

複素振動では指数の絶対値は $1$ であっても、相が速く変わる領域の積分は打ち消し合う。$f$ をコンパクト台とし、その上で $\phi'\ne0$ なら
$\int fe^{i\lambda\phi}=-(i\lambda)^{-1}\int (f/\phi')'e^{i\lambda\phi}$ となる。台の端で境界項が零なので、この操作を任意回繰り返せる。

:::proposition[非退化停留点]{#prop-mp-stationary-phase}
$f\in C_c^\infty(\R)$、$\phi$ は実滑らかで、その台に停留点 $t_0$ を一つ持ち、$\phi''(t_0)\ne0$ とする。$\sigma=\sgn\phi''(t_0)$ とすれば
$$
\int f(t)e^{i\lambda\phi(t)}\dd t
=e^{i\lambda\phi(t_0)+i\sigma\pi/4}
\sqrt{\frac{2\pi}{\lambda|\phi''(t_0)|}} f(t_0)+O(\lambda^{-1}).
$$
:::

:::proof
停留点近くへ滑らかに切り落とし、$\phi=\phi(t_0)+\sigma u^2/2$ となる変数に替える。振幅 $F$ は $F(0)\eta(u)+uB(u)$ と書ける。$\eta=1$ を原点近くで満たす滑らかなカットオフとし、$B$ は滑らかで台を持つ。$uB$ の積分は $u e^{i\lambda\sigma u^2/2}=(i\lambda\sigma)^{-1}\partial_u e^{i\lambda\sigma u^2/2}$ により $O(\lambda^{-1})$ である。$\eta$ を $1$ に置き換えた振動積分との差も $O(\lambda^{-1})$ であり、$1-\eta$ を $u$ で割って部分積分するとその導関数が可積分になる。Gaussian 積分を $\epsilon-i\lambda\sigma$（$\epsilon>0$）に解析接続し、$\epsilon\downarrow0$ とすれば全実線の Fresnel 積分は $\sqrt{2\pi/\lambda}e^{i\sigma\pi/4}$ となる。この Abel 極限は無限遠の尾が部分積分で一様に $O((\lambda R)^{-1})$ となるので通常の振動積分とも一致する。最後に $F(0)=f(t_0)/\sqrt{|\phi''(t_0)|}$ を用いる。
:::

### 停留点の高次補正

証明と同じ二次の変数 $u$ を使う。切り取った振幅を $F_0(u)=F(u)$、原点の近くで $1$ のカットオフを $\eta(u)$ とし、
$$
F_{j+1}(u)=i\sigma\frac{\dd}{\dd u}\left(\frac{F_j(u)-F_j(0)\eta(u)}{u}\right)
$$
と定義する。括弧内は原点でも滑らかで台はコンパクトである。部分積分により $F_j$ の積分は $F_j(0)\eta$ の積分と $\lambda^{-1}F_{j+1}$ の積分の和になる。$\eta$ を $1$ に置き換えた差は、$(1-\eta)/u$ の微分による部分積分を繰り返すと任意の逆冪で抑えられる。最後に残る $F_N$ の積分は先の命題で $O(\lambda^{-1/2})$ である。従って
$$
\begin{aligned}
\int f(t)e^{i\lambda\phi(t)}\dd t
&=e^{i\lambda\phi(t_0)+i\sigma\pi/4}\sqrt{\frac{2\pi}{\lambda}}
\sum_{j=0}^{N-1}\frac{(i\sigma)^j F^{(2j)}(0)}{2^j j!\lambda^j}
+O(\lambda^{-N-1/2}).
\end{aligned}
$$
係数は、原点の近くで $F_{j+1}^{(\ell)}(0)=i\sigma F_j^{(\ell+2)}(0)/(\ell+2)$ となることから得られる。元の $t$ で相が厳密に二次の場合、$F$ の導関数は振幅の導関数を一定の倍率で変換すればよい。一般の相では変数変換の Jacobian も $F$ の導関数に入る。

:::example[二つの停留点]
円周で相 $\cos\theta$、振幅 $1$ とする。$\theta=0,\pi$ で相は $1,-1$、二階微分は $-1,1$ である。円周を分割して二点の結果を足すと $\int_0^{2\pi}e^{i\lambda\cos\theta}\dd\theta=2\sqrt{2\pi/\lambda}\cos(\lambda-\pi/4)+O(\lambda^{-1})$ となる。停留点が別々に寄与するため、振幅の主項が干渉で零になる $\lambda$ もある。
:::

### 停留点が無い端点でも境界項が出る

有限区間で振幅が端点を避けないときは、一度の部分積分が
$$
\int_a^b f e^{i\lambda\phi}\dd t=
\left[\frac{f e^{i\lambda\phi}}{i\lambda\phi'}\right]_a^b
-\frac{1}{i\lambda}\int_a^b (f/\phi')' e^{i\lambda\phi}\dd t
$$
となる。$\phi'$ が全区間で零でなく、係数が滑らかなら残りをもう一度部分積分して $O(\lambda^{-2})$ と評価できる。例えば $\int_0^1 t e^{i\lambda t}\dd t$ は $e^{i\lambda}/(i\lambda)-(e^{i\lambda}-1)/(i\lambda)^2$ である。台が内部にある場合は境界で全ての微分が消えるが、有限区間に切っただけの振幅では端から寄与が出る。

### 端点項の繰り返しと符号

相 $\phi$ の微分が区間で零にならなければ、$e^{i\lambda\phi}=(i\lambda\phi^\prime)^{-1}\partial_x e^{i\lambda\phi}$ とする。一度の部分積分で端点項が得られ、内側の振幅は $-\partial_x(f/\phi^\prime)$ になる。その振幅に同じ操作を行うと、二次以降の端点項が得られる。台が区間内部に収まればすべての端点項が消え、十分滑らかな場合は任意の負の冪で評価できる。

:::example[線形相の端点項（1）]
区間 $[0,1]$ の線形相 $ut$ で、$t^m$ の積分に部分積分を繰り返す。基準の値で計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

右端から来る一次の項は振幅の右端の値に比例する。

$$
\begin{aligned}
R_1(u) &= - \frac{i e^{i u}}{u}
\end{aligned}
$$

二次の項には振幅の微分が入り、部分積分の回数に応じて符号が反転する。

$$
\begin{aligned}
R_2(u)\big|_{m=1} &= \frac{e^{i u}}{u^{2}}
\end{aligned}
$$

多項式振幅では微分を続けると積分が消え、残った境界項の和が厳密解になる。

$$
\begin{aligned}
I_m(u)\big|_{m=1} &= - \frac{i e^{i u}}{u} + \frac{e^{i u}}{u^{2}} - \frac{1}{u^{2}}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
R_1(u)\big|_{u=2} &= - \frac{i e^{2 i}}{2}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
R_2(u)\big|_{u=2} &= \frac{e^{2 i}}{4}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
I_m(u)\big|_{u=2} &= - \frac{1}{4} + \frac{\left(1 - 2 i\right) e^{2 i}}{4}
\end{aligned}
$$
:::

:::example[線形相の端点項（2）]
区間 $[0,1]$ の線形相 $ut$ で、$t^m$ の積分に部分積分を繰り返す。値を変えて依存性を調べる。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

右端から来る一次の項は振幅の右端の値に比例する。

$$
\begin{aligned}
R_1(u) &= - \frac{i e^{i u}}{u}
\end{aligned}
$$

二次の項には振幅の微分が入り、部分積分の回数に応じて符号が反転する。

$$
\begin{aligned}
R_2(u)\big|_{m=2} &= \frac{2 e^{i u}}{u^{2}}
\end{aligned}
$$

多項式振幅では微分を続けると積分が消え、残った境界項の和が厳密解になる。

$$
\begin{aligned}
I_m(u)\big|_{m=2} &= - \frac{i e^{i u}}{u} + \frac{2 e^{i u}}{u^{2}} + \frac{2 i e^{i u}}{u^{3}} - \frac{2 i}{u^{3}}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
R_1(u)\big|_{u=2} &= - \frac{i e^{2 i}}{2}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
R_2(u)\big|_{u=2} &= \frac{e^{2 i}}{2}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
I_m(u)\big|_{u=2} &= - \frac{i}{4} + \frac{\left(2 - i\right) e^{2 i}}{4}
\end{aligned}
$$
:::

:::example[線形相の端点項（3）]
区間 $[0,1]$ の線形相 $ut$ で、$t^m$ の積分に部分積分を繰り返す。さらに値を変えて比較する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

右端から来る一次の項は振幅の右端の値に比例する。

$$
\begin{aligned}
R_1(u) &= - \frac{i e^{i u}}{u}
\end{aligned}
$$

二次の項には振幅の微分が入り、部分積分の回数に応じて符号が反転する。

$$
\begin{aligned}
R_2(u)\big|_{m=3} &= \frac{3 e^{i u}}{u^{2}}
\end{aligned}
$$

多項式振幅では微分を続けると積分が消え、残った境界項の和が厳密解になる。

$$
\begin{aligned}
I_m(u)\big|_{m=3} &= - \frac{i e^{i u}}{u} + \frac{3 e^{i u}}{u^{2}} + \frac{6 i e^{i u}}{u^{3}} - \frac{6 e^{i u}}{u^{4}} + \frac{6}{u^{4}}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
R_1(u)\big|_{u=2} &= - \frac{i e^{2 i}}{2}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
R_2(u)\big|_{u=2} &= \frac{3 e^{2 i}}{4}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
I_m(u)\big|_{u=2} &= \frac{3}{8} + \frac{\left(3 + 2 i\right) e^{2 i}}{8}
\end{aligned}
$$
:::

:::example[線形相の端点項（4）]
区間 $[0,1]$ の線形相 $ut$ で、$t^m$ の積分に部分積分を繰り返す。最後にこの値でも計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=4
$$

右端から来る一次の項は振幅の右端の値に比例する。

$$
\begin{aligned}
R_1(u) &= - \frac{i e^{i u}}{u}
\end{aligned}
$$

二次の項には振幅の微分が入り、部分積分の回数に応じて符号が反転する。

$$
\begin{aligned}
R_2(u)\big|_{m=4} &= \frac{4 e^{i u}}{u^{2}}
\end{aligned}
$$

多項式振幅では微分を続けると積分が消え、残った境界項の和が厳密解になる。

$$
\begin{aligned}
I_m(u)\big|_{m=4} &= - \frac{i e^{i u}}{u} + \frac{4 e^{i u}}{u^{2}} + \frac{12 i e^{i u}}{u^{3}} - \frac{24 e^{i u}}{u^{4}} - \frac{24 i e^{i u}}{u^{5}} + \frac{24 i}{u^{5}}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
R_1(u)\big|_{u=2} &= - \frac{i e^{2 i}}{2}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
R_2(u)\big|_{u=2} &= e^{2 i}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
I_m(u)\big|_{u=2} &= \frac{\left(-2 + i\right) e^{2 i}}{4} + \frac{3 i}{4}
\end{aligned}
$$
:::

## 鞍点法（最急降下法）

### 複素の相に沿って経路を選ぶ

解析的な相 $\Phi$ の非退化な臨界点 $z_0$ について、$\Phi(z)-\Phi(z_0)=-u(z)^2/2$ と書ける。実際、二次の零を除いた正則関数は近傍で零を持たず、局所的な対数から平方根を取れる。そこで $u$ が実数となる向きに経路を通せば $\Im\Phi$ は一定、$\Re\Phi$ は二次で減少する。局所の寄与は Laplace の方法と同じ Gaussian 積分になる。振幅は元のものに $z'(u)$ を掛ける。どの平方根を選ぶかは通過方向で決まり、$z'(0)$ の符号を経路の向きと揃える必要がある。大域的な変形には、極や切断を横切らないことと、経路の残りの積分が小さいことを別に確認する。

具体的に、変形した経路が鞍点付近で $u\in[-b,b]$ を増加方向に通り、残りで $\Re\Phi\le\Re\Phi(z_0)-\eta$ とする。残りの経路で $|a(z)\dd z|$ が可積分なら、そこの寄与は $O(e^{\lambda(\Re\Phi(z_0)-\eta)})$ である。局所の積分では $a(z(u))z'(u)$ を展開し、奇数項を対称性で消す。そこで
$$
\int a(z)e^{\lambda\Phi(z)}\dd z
= e^{\lambda\Phi(z_0)}\sqrt{\frac{2\pi}{\lambda}}
\left(a(z_0)z'(0)+O(\lambda^{-1})\right),
\qquad z'(0)^2=-\frac{1}{\Phi''(z_0)}
$$
を得る。経路の向きを逆にすると $z'(0)$ の符号も逆になる。複数の鞍点を通る場合は各近傍に分割し、残りに同じ指数の差があることを確認して寄与を足す。鞍点が合流して $\Phi''$ が零になる極限では二次の変数が局所座標にならないため、相の低次項をまとめて残す必要がある。

:::example[ずれた Gaussian]
$\lambda>0$、$a\in\R$ に対して $\Phi(z)=-z^2/2+iaz$ の鞍点は $ia$ である。実線を $\R+ia$ へ移せる。長方形の縦辺では $z=\pm R+iy$ で $\Re\Phi=-R^2/2+y^2/2-ay$ であり、$y$ の範囲が有界なので寄与は消える。$z=x+ia$ にすると $\Phi=-x^2/2-a^2/2$ である。従って積分は $\sqrt{2\pi/\lambda}e^{-\lambda a^2/2}$ となる。ここでは二次より高次の位相項がないため主項が厳密値である。
:::

:::example[係数抽出で現れる非二次の鞍点]
$e^{nz}$ の $z^n$ の係数 $c_n$ を Cauchy の積分で取り出す。半径 $R$ の円の積分では、指数の実部が最大となる正の実軸上で $R-\log R$ が現れる。この値を最小にする半径は $R=1$ である。そこで $z=e^{i\theta}$ と置くと
$$
c_n=\frac1{2\pi}\int_{-\pi}^{\pi}e^{n(e^{i\theta}-i\theta)}\dd\theta .
$$
相の $\theta$ 微分は $\theta=0$ で零であり、そこでの相は $1$、二階微分は $-1$ である。半径を選んだことで、円上の指数の実部が最大となる点で虚部の一階変化も消えた。$|\theta|\ge\epsilon>0$ では実部 $\cos\theta\le\cos\epsilon<1$ により寄与が指数的に小さい。

近傍では $\theta=u/\sqrt n$ と置いて $n(e^{i\theta}-i\theta-1)\to-u^2/2$ を得る。また $\cos\theta-1\le-c\theta^2$ となる $c>0$ があるので、拡大後の被積分関数の絶対値は $e^{-cu^2}$ で抑えられる。優収束で積分すると $c_n\sim e^n/\sqrt{2\pi n}$ となる。Taylor 係数からの厳密値 $n^n/n!$ に Stirling の公式を代入しても同じ主項を得る。この方法は係数の初等的な閉じた形が手に入らない場合にも、円の半径を調整して同様の鞍点を探す手順を与える。
:::

### 複素の平行移動が許される条件

二次相では平方完成して鞍点を出せる。直線の経路を平行移動するときには、有限の長方形の側辺が消えることを示す。例えば $e^{-mx^2+i\lambda x}$ の側辺を $x=\pm R+iy$ と置くと、有限範囲の $y$ で絶対値は $e^{-mR^2+my^2-\lambda y}$ となる。従って $R\to\infty$ で側辺は消える。鞍点を通る直線上では振動が指数の定数因子に移り、実 Gaussian が残る。

:::example[ずれた鞍点の比較（1）]
積分 $\int_\R e^{-mx^2+iu x}\dd x$ を考える。基準の値で計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

相を微分して鞍点の座標を求める。

$$
\begin{aligned}
z_s\big|_{m=1} &= \frac{i u}{2}
\end{aligned}
$$

平方完成した後、鞍点での指数の値を分離する。

$$
\begin{aligned}
-m z_s^2+iu z_s\big|_{m=1} &= - \frac{u^{2}}{4}
\end{aligned}
$$

鞍点に合わせて経路を移し、残った実 Gaussian を積分する。

$$
\begin{aligned}
I(u)\big|_{m=1} &= \sqrt{\pi} e^{- \frac{u^{2}}{4}}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
z_s\big|_{u=2} &= i \\
\frac{\partial^1 (z_s)}{\partial u^1} &= \frac{i}{2} \\
\left.\frac{\partial^1 (z_s)}{\partial u^1}\right|_{u=2} &= \frac{i}{2}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
-m z_s^2+iu z_s\big|_{u=2} &= -1 \\
\frac{\partial^1 (-m z_s^2+iu z_s)}{\partial u^1} &= - \frac{u}{2} \\
\left.\frac{\partial^1 (-m z_s^2+iu z_s)}{\partial u^1}\right|_{u=2} &= -1
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
I(u)\big|_{u=2} &= \frac{\sqrt{\pi}}{e} \\
\frac{\partial^1 (I(u))}{\partial u^1} &= - \frac{\sqrt{\pi} u e^{- \frac{u^{2}}{4}}}{2} \\
\left.\frac{\partial^1 (I(u))}{\partial u^1}\right|_{u=2} &= - \frac{\sqrt{\pi}}{e}
\end{aligned}
$$
:::

:::example[ずれた鞍点の比較（2）]
積分 $\int_\R e^{-mx^2+iu x}\dd x$ を考える。値を変えて依存性を調べる。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

相を微分して鞍点の座標を求める。

$$
\begin{aligned}
z_s\big|_{m=2} &= \frac{i u}{4}
\end{aligned}
$$

平方完成した後、鞍点での指数の値を分離する。

$$
\begin{aligned}
-m z_s^2+iu z_s\big|_{m=2} &= - \frac{u^{2}}{8}
\end{aligned}
$$

鞍点に合わせて経路を移し、残った実 Gaussian を積分する。

$$
\begin{aligned}
I(u)\big|_{m=2} &= \frac{\sqrt{2} \sqrt{\pi} e^{- \frac{u^{2}}{8}}}{2}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
z_s\big|_{u=2} &= \frac{i}{2} \\
\frac{\partial^1 (z_s)}{\partial u^1} &= \frac{i}{4} \\
\left.\frac{\partial^1 (z_s)}{\partial u^1}\right|_{u=2} &= \frac{i}{4}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
-m z_s^2+iu z_s\big|_{u=2} &= - \frac{1}{2} \\
\frac{\partial^1 (-m z_s^2+iu z_s)}{\partial u^1} &= - \frac{u}{4} \\
\left.\frac{\partial^1 (-m z_s^2+iu z_s)}{\partial u^1}\right|_{u=2} &= - \frac{1}{2}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
I(u)\big|_{u=2} &= \frac{\sqrt{2} \sqrt{\pi}}{2 e^{\frac{1}{2}}} \\
\frac{\partial^1 (I(u))}{\partial u^1} &= - \frac{\sqrt{2} \sqrt{\pi} u e^{- \frac{u^{2}}{8}}}{8} \\
\left.\frac{\partial^1 (I(u))}{\partial u^1}\right|_{u=2} &= - \frac{\sqrt{2} \sqrt{\pi}}{4 e^{\frac{1}{2}}}
\end{aligned}
$$
:::

:::example[ずれた鞍点の比較（3）]
積分 $\int_\R e^{-mx^2+iu x}\dd x$ を考える。さらに値を変えて比較する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

相を微分して鞍点の座標を求める。

$$
\begin{aligned}
z_s\big|_{m=3} &= \frac{i u}{6}
\end{aligned}
$$

平方完成した後、鞍点での指数の値を分離する。

$$
\begin{aligned}
-m z_s^2+iu z_s\big|_{m=3} &= - \frac{u^{2}}{12}
\end{aligned}
$$

鞍点に合わせて経路を移し、残った実 Gaussian を積分する。

$$
\begin{aligned}
I(u)\big|_{m=3} &= \frac{\sqrt{3} \sqrt{\pi} e^{- \frac{u^{2}}{12}}}{3}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
z_s\big|_{u=2} &= \frac{i}{3} \\
\frac{\partial^1 (z_s)}{\partial u^1} &= \frac{i}{6} \\
\left.\frac{\partial^1 (z_s)}{\partial u^1}\right|_{u=2} &= \frac{i}{6}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
-m z_s^2+iu z_s\big|_{u=2} &= - \frac{1}{3} \\
\frac{\partial^1 (-m z_s^2+iu z_s)}{\partial u^1} &= - \frac{u}{6} \\
\left.\frac{\partial^1 (-m z_s^2+iu z_s)}{\partial u^1}\right|_{u=2} &= - \frac{1}{3}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
I(u)\big|_{u=2} &= \frac{\sqrt{3} \sqrt{\pi}}{3 e^{\frac{1}{3}}} \\
\frac{\partial^1 (I(u))}{\partial u^1} &= - \frac{\sqrt{3} \sqrt{\pi} u e^{- \frac{u^{2}}{12}}}{18} \\
\left.\frac{\partial^1 (I(u))}{\partial u^1}\right|_{u=2} &= - \frac{\sqrt{3} \sqrt{\pi}}{9 e^{\frac{1}{3}}}
\end{aligned}
$$
:::

:::example[ずれた鞍点の比較（4）]
積分 $\int_\R e^{-mx^2+iu x}\dd x$ を考える。最後にこの値でも計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=4
$$

相を微分して鞍点の座標を求める。

$$
\begin{aligned}
z_s\big|_{m=4} &= \frac{i u}{8}
\end{aligned}
$$

平方完成した後、鞍点での指数の値を分離する。

$$
\begin{aligned}
-m z_s^2+iu z_s\big|_{m=4} &= - \frac{u^{2}}{16}
\end{aligned}
$$

鞍点に合わせて経路を移し、残った実 Gaussian を積分する。

$$
\begin{aligned}
I(u)\big|_{m=4} &= \frac{\sqrt{\pi} e^{- \frac{u^{2}}{16}}}{2}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
z_s\big|_{u=2} &= \frac{i}{4} \\
\frac{\partial^1 (z_s)}{\partial u^1} &= \frac{i}{8} \\
\left.\frac{\partial^1 (z_s)}{\partial u^1}\right|_{u=2} &= \frac{i}{8}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
-m z_s^2+iu z_s\big|_{u=2} &= - \frac{1}{4} \\
\frac{\partial^1 (-m z_s^2+iu z_s)}{\partial u^1} &= - \frac{u}{8} \\
\left.\frac{\partial^1 (-m z_s^2+iu z_s)}{\partial u^1}\right|_{u=2} &= - \frac{1}{4}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
I(u)\big|_{u=2} &= \frac{\sqrt{\pi}}{2 e^{\frac{1}{4}}} \\
\frac{\partial^1 (I(u))}{\partial u^1} &= - \frac{\sqrt{\pi} u e^{- \frac{u^{2}}{16}}}{16} \\
\left.\frac{\partial^1 (I(u))}{\partial u^1}\right|_{u=2} &= - \frac{\sqrt{\pi}}{8 e^{\frac{1}{4}}}
\end{aligned}
$$
:::

## Euler–Maclaurin の公式と Bernoulli 数

### 和の端点補正

Bernoulli 多項式 $B_n(t)$ を $s e^{ts}/(e^s-1)=\sum_{n\ge0}B_n(t)s^n/n!$ で定義し、$B_n=B_n(0)$ とする。母関数を $t$ で微分すれば $B_n'=nB_{n-1}$ となる。$t=1,0$ での差は $s$ なので、$n\ne1$ で端点値は一致する。$s/(e^s-1)+s/2$ は偶関数なので、$B_1=-1/2$、奇数添字 $n>1$ は零である。分母の級数を掛けて係数を比較すると $B_2=1/6$、$B_4=-1/30$ を得る。

:::proposition[Euler–Maclaurin の公式]{#prop-mp-euler-maclaurin}
整数 $a<b$、整数 $m\ge1$ と $f\in C^{2m}([a,b])$ に対して
$$
\begin{aligned}
\sum_{k=a}^b f(k)
&=\int_a^b f(t)\dd t+\tfrac12(f(a)+f(b))\\
&\quad+\sum_{j=1}^m\frac{B_{2j}}{(2j)!}
   (f^{(2j-1)}(b)-f^{(2j-1)}(a))
-\frac1{(2m)!}\int_a^b B_{2m}(\{t\})f^{(2m)}(t)\dd t .
\end{aligned}
$$
:::

:::proof
各単位区間で $(t-k-1/2)f'(t)$ を部分積分すると端点値の平均から積分を引いた値になる。区間を足すと誤差が $\int B_1(\{t\}) f'(t)\dd t$ である。$B_n'=n B_{n-1}$ により、各区間内で部分積分を順に行う。$n\ge2$ では内部端点の寄与が連続性で相殺され、全区間の両端だけが残る。奇数の Bernoulli 数は零なので、最後の偶数次数まで実行すると表示式になる。
:::

:::example[平方和]
$f(t)=t^2$、$a=0$、$b=N$ とする。$m=2$ なら最終剰余は零で、$\int f=N^3/3$、端点平均は $N^2/2$、最初の補正は $N/6$ となる。これらを足せば $N(N+1)(2N+1)/6$ であり、例えば $N=4$ では $30$ となる。
:::

### 無限和の尾を評価する

$f(t)=t^{-2}$、$a=N$、$b\to\infty$ を式に入れると $\sum_{k>N}k^{-2}=N^{-1}-(2N^2)^{-1}+(6N^3)^{-1}+O(N^{-5})$ である。$f^{(4)}=120t^{-6}$ なので $m=2$ の剰余は $B_4(\{t\})$ の上限を使って $O(N^{-5})$ になり、$B_4$ 端点項も同じ次数である。例えば $N=10$ では三項の和は $0.0951667$ 程度である。$k=N$ の項を尾に含めるなら $N^{-2}$ を足すため、二次の符号が変わる。

### Stirling の展開を対数で高次まで出す

$f(t)=\log t$ に $[1,N]$ で Euler–Maclaurin を適用する。積分と上側の端点補正の和は $(N+1/2)\log N-N+1$ である。上端で $f^{(2j-1)}(N)=(2j-2)!/N^{2j-1}$ なので、その Bernoulli 項は $B_{2j}/[2j(2j-1)N^{2j-1}]$ となる。下端から来た数はすべて $N$ に依存しない。$m+1$ 次まで取った剰余の積分は $f^{(2m+2)}=-(2m+1)!t^{-2m-2}$ のため $N\to\infty$ で定数に収束し、その尾は $O(N^{-2m-1})$ である。従って
$$
\log(N!)=(N+1/2)\log N-N+\tfrac12\log(2\pi)
+\sum_{j=1}^m\frac{B_{2j}}{2j(2j-1)N^{2j-1}}
+O(N^{-2m-1})
$$
となる。定数 $\tfrac12\log(2\pi)$ は本章の Laplace 法による主項と比較して決めた。例えば $m=2$ では補正は $1/(12N)-1/(360N^3)$ である。元の階乗に戻すときは指数関数を展開するので、$1/N^2$ の係数 $1/288$ は対数に二次項が無くても生じる。

### 和と積分の差の境界補正

整数 $k$ で区切り、$B_1(t-\lfloor t\rfloor)$ を使って各区間で部分積分する。この周期関数は区間の内部では導関数が $1$、左右端では値が $-1/2$ と $1/2$ になる。端で出る半分ずつの寄与を合成すると台形和が現れる。高次の周期 Bernoulli 多項式は端の値が一致する。そのため次の部分積分では内部節点が相殺され、外側二端の導関数だけが残る。

:::example[冪和の各補正（1）]
$k=0,\ldots,N$ の $k^m$ の和を、積分と端点項に分ける。基準の値で計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

連続積分を計算する。

$$
\begin{aligned}
I_N\big|_{m=1} &= \frac{N^{2}}{2}
\end{aligned}
$$

端点の平均と、$B_2/2!$ から生じる一階導関数の差を求める。

$$
\begin{aligned}
C_N\big|_{m=1} &= \frac{N}{2}
\end{aligned}
$$

有限和を評価して補正後の差を確かめる。次数が高い場合は四階以上の Bernoulli 数の項が必要になる。

$$
\begin{aligned}
S_N\big|_{m=1} &= \sum_{n=0}^{N} n = \frac{N \left(N + 1\right)}{2}
\end{aligned}
$$

観測点 $N$=$5$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
I_N\big|_{N=5} &= \frac{25}{2}
\end{aligned}
$$

観測点 $N$=$5$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
C_N\big|_{N=5} &= \frac{5}{2}
\end{aligned}
$$

観測点 $N$=$5$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
S_N\big|_{N=5} &= \sum_{n=0}^{5} n
\end{aligned}
$$
:::

:::example[冪和の各補正（2）]
$k=0,\ldots,N$ の $k^m$ の和を、積分と端点項に分ける。値を変えて依存性を調べる。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

連続積分を計算する。

$$
\begin{aligned}
I_N\big|_{m=2} &= \frac{N^{3}}{3}
\end{aligned}
$$

端点の平均と、$B_2/2!$ から生じる一階導関数の差を求める。

$$
\begin{aligned}
C_N\big|_{m=2} &= \frac{N \left(3 N + 1\right)}{6}
\end{aligned}
$$

有限和を評価して補正後の差を確かめる。次数が高い場合は四階以上の Bernoulli 数の項が必要になる。

$$
\begin{aligned}
S_N\big|_{m=2} &= \sum_{n=0}^{N} n^{2} = \frac{N \left(N + 1\right) \left(2 N + 1\right)}{6}
\end{aligned}
$$

観測点 $N$=$5$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
I_N\big|_{N=5} &= \frac{125}{3}
\end{aligned}
$$

観測点 $N$=$5$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
C_N\big|_{N=5} &= \frac{40}{3}
\end{aligned}
$$

観測点 $N$=$5$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
S_N\big|_{N=5} &= \sum_{n=0}^{5} n^{2}
\end{aligned}
$$
:::

:::example[冪和の各補正（3）]
$k=0,\ldots,N$ の $k^m$ の和を、積分と端点項に分ける。さらに値を変えて比較する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

連続積分を計算する。

$$
\begin{aligned}
I_N\big|_{m=3} &= \frac{N^{4}}{4}
\end{aligned}
$$

端点の平均と、$B_2/2!$ から生じる一階導関数の差を求める。

$$
\begin{aligned}
C_N\big|_{m=3} &= \frac{N^{2} \left(2 N + 1\right)}{4}
\end{aligned}
$$

有限和を評価して補正後の差を確かめる。次数が高い場合は四階以上の Bernoulli 数の項が必要になる。

$$
\begin{aligned}
S_N\big|_{m=3} &= \sum_{n=0}^{N} n^{3} = \frac{N^{2} \left(N + 1\right)^{2}}{4}
\end{aligned}
$$

観測点 $N$=$5$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
I_N\big|_{N=5} &= \frac{625}{4}
\end{aligned}
$$

観測点 $N$=$5$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
C_N\big|_{N=5} &= \frac{275}{4}
\end{aligned}
$$

観測点 $N$=$5$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
S_N\big|_{N=5} &= \sum_{n=0}^{5} n^{3}
\end{aligned}
$$
:::

:::example[冪和の各補正（4）]
$k=0,\ldots,N$ の $k^m$ の和を、積分と端点項に分ける。最後にこの値でも計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=4
$$

連続積分を計算する。

$$
\begin{aligned}
I_N\big|_{m=4} &= \frac{N^{5}}{5}
\end{aligned}
$$

端点の平均と、$B_2/2!$ から生じる一階導関数の差を求める。

$$
\begin{aligned}
C_N\big|_{m=4} &= \frac{N^{3} \left(3 N + 2\right)}{6}
\end{aligned}
$$

有限和を評価して補正後の差を確かめる。次数が高い場合は四階以上の Bernoulli 数の項が必要になる。

$$
\begin{aligned}
S_N\big|_{m=4} &= \sum_{n=0}^{N} n^{4} = \frac{N \left(N + 1\right) \left(2 N + 1\right) \left(3 N^{2} + 3 N - 1\right)}{30}
\end{aligned}
$$

観測点 $N$=$5$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
I_N\big|_{N=5} &= 625
\end{aligned}
$$

観測点 $N$=$5$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
C_N\big|_{N=5} &= \frac{2125}{6}
\end{aligned}
$$

観測点 $N$=$5$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
S_N\big|_{N=5} &= \sum_{n=0}^{5} n^{4} = 979
\end{aligned}
$$

この場合、$I_N+C_N$ は $S_N$ より $N/30$ 大きい。三階導関数が $24t$ なので、$B_4/4!$ の補正は $-24N/(30\cdot4!)=-N/30$ である。$m=3$ までの冪ではこの差が零だった。ここでは六階導関数が零なので、Euler–Maclaurin の公式を $B_6$ の項まで取れば剰余が厳密に零になる。$N=5$ なら $625+2125/6-1/6=979$ である。
:::

## 発散級数の和と Borel 和

### 階乗成長を割ってから積分する

:::definition[Borel 和]{#def-mp-borel}
形式級数 $\sum a_n x^n$ に対して $\widehat f(s)=\sum a_n s^n/n!$ を :term[ボレル変換]{en="Borel transform" yomi="ぼれるへんかん"} とする。原点近くで収束し、正実軸に沿って解析接続でき、$|\widehat f(s)|\le C e^{cs}$ なら、$x>0$ で $1/x>c$ に対して $x^{-1}\int_0^\infty e^{-s/x}\widehat f(s)\dd s$ を :term[ボレル和]{en="Borel sum" yomi="ぼれるわ"} とする。
:::

Watson の補題を $\alpha=1$ で適用すると、Borel 和の $x\downarrow0$ の展開は元の係数を再現する。$\Gamma(n+1)=n!$ が分母の階乗を打ち消す。Borel 変換に正実軸上の極がある場合は定義通りの積分ができない。上下から迂回する積分を使うなら、それぞれを指定する必要がある。

通常の冪級数が正の収束半径 $R$ を持つときは、$0<x<R$ で Borel 和もその和に等しい。$x<r<R$ を取り、係数を $|a_n|\le M_r r^{-n}$ で抑える。すると Borel 変換は全平面で収束し、$|\widehat f(s)|\le M_r e^{|s|/r}$ である。正実軸上では $e^{-s/x}$ を掛けた各項の絶対値の積分が合計 $\sum|a_n|x^n<\infty$ となる。従って積分と和を交換でき、各項の積分が $a_nx^n$ を返す。この議論は元の級数の収束範囲内で二つの和が一致することを保証する。

:::example[交代する階乗級数と極]
$\sum(-1)^n n! x^n$ の Borel 変換は $(1+s)^{-1}$ である。和は $\int_0^\infty e^{-t}/(1+xt)\dd t$ になる。係数が非交代の $\sum n!x^n$ では変換は $(1-s)^{-1}$ となる。正実軸上の $s=1$ の留数は $-e^{-1/x}/x$ なので、上下の迂回で得る値の差の大きさは $2\pi e^{-1/x}/x$ である。これは任意の $x^N$ より小さく、冪の係数から差を検出できない。
:::

### 最小項付近と積分の剰余

係数が $m^j j!$ の場合、隣り合う項の大きさの比は $m(j+1)t$ である。$t>0$ を固定して項数を増やすと、この比はやがて $1$ を超える。ところが $t\to0$ で項数を固定すれば、積分の剰余は次の項と同じ次数に抑えられる。これは幾何恒等式の剰余を $\int_0^\infty e^{-s}(\cdots)\dd s$ に入れ、分母 $1+mts\ge1$ を使った結果である。

:::example[Borel 積分での符号と誤差（1）]
正の $t$ で $\int_0^\infty e^{-s}/(1+mts)\dd s$ を使う。基準の値で計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

Borel 変換は幾何級数を和にしたものになる。

$$
\begin{aligned}
\mathcal B(u)\big|_{m=1} &= \frac{1}{u + 1}
\end{aligned}
$$

三次までの項を戻して有限和を作る。

$$
\begin{aligned}
S_3(t)\big|_{m=1} &= - 6 t^{3} + 2 t^{2} - t + 1
\end{aligned}
$$

次の項の大きさが上側の誤差境界になる。剰余の符号は正の場合である。

$$
\begin{aligned}
R_3(t)\text{ の上界}\big|_{m=1} &= 24 t^{4}
\end{aligned}
$$

観測点 $t$=$\frac{1}{20}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
S_3(t)\big|_{t=\frac{1}{20}} &= \frac{3817}{4000}
\end{aligned}
$$

観測点 $t$=$\frac{1}{20}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
R_3(t)\text{ の上界}\big|_{t=\frac{1}{20}} &= \frac{3}{20000}
\end{aligned}
$$
:::

:::example[Borel 積分での符号と誤差（2）]
正の $t$ で $\int_0^\infty e^{-s}/(1+mts)\dd s$ を使う。値を変えて依存性を調べる。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

Borel 変換は幾何級数を和にしたものになる。

$$
\begin{aligned}
\mathcal B(u)\big|_{m=2} &= \frac{1}{2 u + 1}
\end{aligned}
$$

三次までの項を戻して有限和を作る。

$$
\begin{aligned}
S_3(t)\big|_{m=2} &= - 48 t^{3} + 8 t^{2} - 2 t + 1
\end{aligned}
$$

次の項の大きさが上側の誤差境界になる。剰余の符号は正の場合である。

$$
\begin{aligned}
R_3(t)\text{ の上界}\big|_{m=2} &= 384 t^{4}
\end{aligned}
$$

観測点 $t$=$\frac{1}{20}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
S_3(t)\big|_{t=\frac{1}{20}} &= \frac{457}{500}
\end{aligned}
$$

観測点 $t$=$\frac{1}{20}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
R_3(t)\text{ の上界}\big|_{t=\frac{1}{20}} &= \frac{3}{1250}
\end{aligned}
$$
:::

:::example[Borel 積分での符号と誤差（3）]
正の $t$ で $\int_0^\infty e^{-s}/(1+mts)\dd s$ を使う。さらに値を変えて比較する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

Borel 変換は幾何級数を和にしたものになる。

$$
\begin{aligned}
\mathcal B(u)\big|_{m=3} &= \frac{1}{3 u + 1}
\end{aligned}
$$

三次までの項を戻して有限和を作る。

$$
\begin{aligned}
S_3(t)\big|_{m=3} &= - 162 t^{3} + 18 t^{2} - 3 t + 1
\end{aligned}
$$

次の項の大きさが上側の誤差境界になる。剰余の符号は正の場合である。

$$
\begin{aligned}
R_3(t)\text{ の上界}\big|_{m=3} &= 1944 t^{4}
\end{aligned}
$$

観測点 $t$=$\frac{1}{20}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
S_3(t)\big|_{t=\frac{1}{20}} &= \frac{3499}{4000}
\end{aligned}
$$

観測点 $t$=$\frac{1}{20}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
R_3(t)\text{ の上界}\big|_{t=\frac{1}{20}} &= \frac{243}{20000}
\end{aligned}
$$
:::

:::example[Borel 積分での符号と誤差（4）]
正の $t$ で $\int_0^\infty e^{-s}/(1+mts)\dd s$ を使う。最後にこの値でも計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=4
$$

Borel 変換は幾何級数を和にしたものになる。

$$
\begin{aligned}
\mathcal B(u)\big|_{m=4} &= \frac{1}{4 u + 1}
\end{aligned}
$$

三次までの項を戻して有限和を作る。

$$
\begin{aligned}
S_3(t)\big|_{m=4} &= - 384 t^{3} + 32 t^{2} - 4 t + 1
\end{aligned}
$$

次の項の大きさが上側の誤差境界になる。剰余の符号は正の場合である。

$$
\begin{aligned}
R_3(t)\text{ の上界}\big|_{m=4} &= 6144 t^{4}
\end{aligned}
$$

観測点 $t$=$\frac{1}{20}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
S_3(t)\big|_{t=\frac{1}{20}} &= \frac{104}{125}
\end{aligned}
$$

観測点 $t$=$\frac{1}{20}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
R_3(t)\text{ の上界}\big|_{t=\frac{1}{20}} &= \frac{24}{625}
\end{aligned}
$$
:::

## 演習問題

::::exercise{level=basic}
$(x+2)^{-1}$ を $x\to+\infty$ で三項まで展開せよ。
:::solution
$u=2/x$ として $(1+u)^{-1}=1-u+u^2-u^3/(1+u)$ を使う。全体に $x^{-1}$ を掛けると $x^{-1}-2x^{-2}+4x^{-3}+O(x^{-4})$ となる。
:::
::::

::::exercise{level=basic}
$\int_{-1}^1 e^{-2\lambda t^2}\dd t$ の主項を求めよ。
:::solution
唯一の極小点は $0$、$\phi''(0)=4$ なので式 :eqref[eq-mp-laplace-leading] で主項は $\sqrt{\pi/(2\lambda)}$ である。実線まで積分を延ばした誤差は $2\int_1^\infty t e^{-2\lambda t^2}\dd t=e^{-2\lambda}/(2\lambda)$ 以下である。
:::
::::

::::exercise{level=basic}
$\int_0^\infty e^{-xt}t(1+3t)\dd t$ を求め、Watson の補題と照合せよ。
:::solution
$u=xt$ と置くと $x^{-2}\int ue^{-u}\dd u+3x^{-3}\int u^2e^{-u}\dd u=x^{-2}+6x^{-3}$ となる。Watson の係数は $\alpha=2$ より $\Gamma(2)=1$ と $3\Gamma(3)=6$ で一致する。
:::
::::

::::exercise{level=basic}
Bernoulli 数 $B_1,B_2$ を母関数から求めよ。
:::solution
$(1+s/2+s^2/6+\cdots)(1+B_1s+B_2s^2/2+\cdots)=1$ である。一次で $B_1=-1/2$、二次で $B_2/2-1/4+1/6=0$ となり $B_2=1/6$ を得る。
:::
::::

::::exercise{level=advanced}
$\int_0^1 e^{-\lambda t}(1+t^2)\dd t$ の打ち切り誤差を求めよ。
:::solution
無限遠まで延ばした値は $\lambda^{-1}+2\lambda^{-3}$ である。尾では $t=1+u$ と置くと被積分式は $e^{-\lambda}e^{-\lambda u}(2+2u+u^2)$ になる。従って元の値はこの二項から $e^{-\lambda}(2/\lambda+2/\lambda^2+2/\lambda^3)$ を引いた値である。
:::
::::

::::exercise{level=advanced}
滑らかな振幅 $f$ が $0$ で零になり、台が有界とする。$\int f(t)e^{i\lambda t^2/2}\dd t=O(\lambda^{-1})$ を部分積分で示せ。
:::solution
$f(t)=tB(t)$ と書ける。$t=0$ では $B(0)=f'(0)$ とし、Taylor の積分剰余で滑らかさが分かる。$f$ の台の外は $B=0$ である。指数の導関数を使って積分を $-(i\lambda)^{-1}\int B'(t)e^{i\lambda t^2/2}\dd t$ と書け、$B'$ は可積分なので結論を得る。
:::
::::

::::exercise{level=advanced}
$\sum_{k=0}^N k^3$ を Euler–Maclaurin の公式で求めよ。
:::solution
積分は $N^4/4$、端点補正は $N^3/2$、$B_2$ の項は $3N^2/12=N^2/4$ である。$f'''$ は一定でその端点差は零、$f^{(4)}=0$ なので剰余も零である。和は $N^2(N+1)^2/4$ となる。
:::
::::

::::exercise{level=advanced}
$\sum_{n\ge0}(-2)^n n! x^n$ の Borel 和を求め、二次打ち切りの誤差上限を与えよ。
:::solution
Borel 変換は $(1+2s)^{-1}$ で、和は $\int_0^\infty e^{-t}/(1+2xt)\dd t$ である。有限等比和で $1-2xt+4x^2t^2$ を引いた残りは $-(2xt)^3/(1+2xt)$ なので、$1-2x+8x^2$ からの誤差は $8x^3\int e^{-t}t^3\dd t=48x^3$ 以下である。
:::
::::

## 確認テスト

::::quiz{answer="6" tol="1e-9"}
Watson の補題で $f(t)=t^3$ のとき、$x^{-4}$ の係数を求めよ。

:::explanation
$t=u/x$ と変換すると積分は $x^{-4}\int_0^\infty u^3 e^{-u}\dd u$ で、三回の部分積分で積分は $3!=6$ である。
:::
::::

::::quiz
非退化な一つの内部極小点が寄与するとき、振幅が零でなければ指数因子以外にどの冪が現れるか。

- [ ] $\lambda^{-1}$
- [x] $\lambda^{-1/2}$
- [ ] $\lambda^{-2}$

:::explanation
二次相 $\lambda\phi''(t_0)(t-t_0)^2/2$ が $1$ 程度になる積分幅は $\lambda^{-1/2}$ に比例する。変数をその幅で拡大すると積分測度からこの冪が出る。
:::
::::

::::quiz{answer="30" tol="1e-9"}
Euler–Maclaurin の平方和の式を $N=4$ に適用せよ。

:::explanation
$N^3/3+N^2/2+N/6$ に $4$ を代入し、$64/3+8+4/6=30$ を得る。
:::
::::

::::quiz
漸近展開の一意性から何が分かるか。

- [x] 関数と漸近尺度を固定したとき係数は一意である
- [ ] 係数を指定すると関数は一意である
- [ ] 係数の無限和は必ず収束する

:::explanation
順に残りを尺度で割って極限を取れば係数を決定できる。ただし全尺度より小さい項を関数に足せる場合があり、逆に関数は決まらない。階乗係数の例のように、有限打ち切りの評価は無限和の収束を仮定しない。
:::
::::

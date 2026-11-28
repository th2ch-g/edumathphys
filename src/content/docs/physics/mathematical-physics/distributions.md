---
title: "デルタ関数と超関数"
description: "デルタ関数を超関数として定義し、微分、Fourier 変換、変数変換と点源を計算する。"
notation:
  - symbol: '\langle T,\varphi\rangle'
    meaning: "超関数 T の試験関数への作用"
  - symbol: '\delta_a'
    meaning: "点 a のデルタ超関数"
  - symbol: '\operatorname{Fp}'
    meaning: "Hadamard の有限部分"
---

## この章の目標

- デルタ関数を近似の列と試験関数の作用で記述する。
- 跳びのある関数を微分して点に集中する項を計算する。
- 主値と有限部分の定義を使い分ける。
- 座標の体積要素を含めて点源を規格化する。

試験関数空間の一般論は :ref[def-four-distribution]{text="超関数"} と :ref[def-four-tempered]{text="緩増加超関数"} を参照する。本章では試験関数への作用を使い、跳びと点源が微分方程式に与える寄与を求める。

## デルタ関数の物理的導入と近似関数列

### 短い時間に加わる力

一定の力積 $I$（$\mathrm{N\,s}$）を幅 $\epsilon$ の時間内に与える。力が $I/\epsilon$ なら、運動量は時間幅によらず $I$ だけ変わる。時間幅を零へ縮めるには、力の高さでなく、時間に依存する係数を掛けた積分の極限を取る。

:::definition[デルタ超関数]{#def-mp-dirac}
滑らかな試験関数 $\varphi$ に対し $\langle\delta_a,\varphi\rangle=\varphi(a)$ と定める。これを :term[デルタ超関数]{en="Dirac delta distribution" yomi="でるたちょうかんすう"} と呼ぶ。
:::

$\rho\in L^1(\R)$、$\int\rho=1$ なら $\rho_\epsilon(x)=\epsilon^{-1}\rho((x-a)/\epsilon)$ を積分に入れると $\int\rho(u)\varphi(a+\epsilon u)\dd u$ となる。$|\rho(u)|\|\varphi\|_\infty$ が可積分なので極限は $\varphi(a)$ である。この意味で一様箱 $\mathbf1_{[-1/2,1/2]}$、Gaussian $\pi^{-1/2}e^{-u^2}$、Cauchy 核 $1/[\pi(1+u^2)]$ は同じ点源を与える。後二つの質量はそれぞれ Gaussian 積分、逆正接で確かめられる。$a$ が時刻なら $\delta_a$ の単位は $\mathrm{s^{-1}}$ で、$I\delta_a$ は力の単位を持つ。

### 近似単位の幅とモーメント

近似関数列を作るとき、全積分を $1$ にしたうえで、任意に固定した原点の外側に残る質量が零に近づくことを示す。非負の近似単位を使えば、試験関数の原点での値との差を、近くでの連続性と遠くでの質量に分けて抑えられる。対称な近似単位では一次モーメントが零になる。試験関数の二階微分が有界なら、Taylor の剰余から収束誤差を二次モーメントで抑えられる。

:::example[同じ幅で箱の偶数モーメントを比較（1）]
幅 $2u$、高さ $1/(2u)$ の対称な箱を近似単位に使う。基準の値で計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

二乗以降の偶数モーメントは正側を二倍して計算できる。

$$
\begin{aligned}
\int x^{2m}\rho_u(x)\dd x\big|_{m=1} &= \frac{u^{2}}{3}
\end{aligned}
$$

奇数の冪は正側と負側で正負が逆なので消える。

$$
\begin{aligned}
\int x^{2m-1}\rho_u(x)\dd x &= 0
\end{aligned}
$$

試験関数の Taylor 展開に入れると、偶数微分に掛かる重みは階乗でさらに割る。

$$
\begin{aligned}
c_{2m}\big|_{m=1} &= \frac{u^{2}}{6}
\end{aligned}
$$

観測点 $u$=$\frac{1}{2}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\int x^{2m}\rho_u(x)\dd x\big|_{u=\frac{1}{2}} &= \frac{1}{12} \\
\frac{\partial^1 (\int x^{2m}\rho_u(x)\dd x)}{\partial u^1} &= \frac{2 u}{3} \\
\left.\frac{\partial^1 (\int x^{2m}\rho_u(x)\dd x)}{\partial u^1}\right|_{u=\frac{1}{2}} &= \frac{1}{3} \\
\frac{\partial^2 (\int x^{2m}\rho_u(x)\dd x)}{\partial u^2} &= \frac{2}{3} \\
\left.\frac{\partial^2 (\int x^{2m}\rho_u(x)\dd x)}{\partial u^2}\right|_{u=\frac{1}{2}} &= \frac{2}{3}
\end{aligned}
$$

観測点 $u$=$\frac{1}{2}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
c_{2m}\big|_{u=\frac{1}{2}} &= \frac{1}{24} \\
\frac{\partial^1 (c_{2m})}{\partial u^1} &= \frac{u}{3} \\
\left.\frac{\partial^1 (c_{2m})}{\partial u^1}\right|_{u=\frac{1}{2}} &= \frac{1}{6} \\
\frac{\partial^2 (c_{2m})}{\partial u^2} &= \frac{1}{3} \\
\left.\frac{\partial^2 (c_{2m})}{\partial u^2}\right|_{u=\frac{1}{2}} &= \frac{1}{3}
\end{aligned}
$$
:::

:::example[同じ幅で箱の偶数モーメントを比較（2）]
幅 $2u$、高さ $1/(2u)$ の対称な箱を近似単位に使う。値を変えて依存性を調べる。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

二乗以降の偶数モーメントは正側を二倍して計算できる。

$$
\begin{aligned}
\int x^{2m}\rho_u(x)\dd x\big|_{m=2} &= \frac{u^{4}}{5}
\end{aligned}
$$

奇数の冪は正側と負側で正負が逆なので消える。

$$
\begin{aligned}
\int x^{2m-1}\rho_u(x)\dd x &= 0
\end{aligned}
$$

試験関数の Taylor 展開に入れると、偶数微分に掛かる重みは階乗でさらに割る。

$$
\begin{aligned}
c_{2m}\big|_{m=2} &= \frac{u^{4}}{120}
\end{aligned}
$$

観測点 $u$=$\frac{1}{2}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\int x^{2m}\rho_u(x)\dd x\big|_{u=\frac{1}{2}} &= \frac{1}{80} \\
\frac{\partial^1 (\int x^{2m}\rho_u(x)\dd x)}{\partial u^1} &= \frac{4 u^{3}}{5} \\
\left.\frac{\partial^1 (\int x^{2m}\rho_u(x)\dd x)}{\partial u^1}\right|_{u=\frac{1}{2}} &= \frac{1}{10} \\
\frac{\partial^2 (\int x^{2m}\rho_u(x)\dd x)}{\partial u^2} &= \frac{12 u^{2}}{5} \\
\left.\frac{\partial^2 (\int x^{2m}\rho_u(x)\dd x)}{\partial u^2}\right|_{u=\frac{1}{2}} &= \frac{3}{5}
\end{aligned}
$$

観測点 $u$=$\frac{1}{2}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
c_{2m}\big|_{u=\frac{1}{2}} &= \frac{1}{1920} \\
\frac{\partial^1 (c_{2m})}{\partial u^1} &= \frac{u^{3}}{30} \\
\left.\frac{\partial^1 (c_{2m})}{\partial u^1}\right|_{u=\frac{1}{2}} &= \frac{1}{240} \\
\frac{\partial^2 (c_{2m})}{\partial u^2} &= \frac{u^{2}}{10} \\
\left.\frac{\partial^2 (c_{2m})}{\partial u^2}\right|_{u=\frac{1}{2}} &= \frac{1}{40}
\end{aligned}
$$
:::

:::example[同じ幅で箱の偶数モーメントを比較（3）]
幅 $2u$、高さ $1/(2u)$ の対称な箱を近似単位に使う。さらに値を変えて比較する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

二乗以降の偶数モーメントは正側を二倍して計算できる。

$$
\begin{aligned}
\int x^{2m}\rho_u(x)\dd x\big|_{m=3} &= \frac{u^{6}}{7}
\end{aligned}
$$

奇数の冪は正側と負側で正負が逆なので消える。

$$
\begin{aligned}
\int x^{2m-1}\rho_u(x)\dd x &= 0
\end{aligned}
$$

試験関数の Taylor 展開に入れると、偶数微分に掛かる重みは階乗でさらに割る。

$$
\begin{aligned}
c_{2m}\big|_{m=3} &= \frac{u^{6}}{5040}
\end{aligned}
$$

観測点 $u$=$\frac{1}{2}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\int x^{2m}\rho_u(x)\dd x\big|_{u=\frac{1}{2}} &= \frac{1}{448} \\
\frac{\partial^1 (\int x^{2m}\rho_u(x)\dd x)}{\partial u^1} &= \frac{6 u^{5}}{7} \\
\left.\frac{\partial^1 (\int x^{2m}\rho_u(x)\dd x)}{\partial u^1}\right|_{u=\frac{1}{2}} &= \frac{3}{112} \\
\frac{\partial^2 (\int x^{2m}\rho_u(x)\dd x)}{\partial u^2} &= \frac{30 u^{4}}{7} \\
\left.\frac{\partial^2 (\int x^{2m}\rho_u(x)\dd x)}{\partial u^2}\right|_{u=\frac{1}{2}} &= \frac{15}{56}
\end{aligned}
$$

観測点 $u$=$\frac{1}{2}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
c_{2m}\big|_{u=\frac{1}{2}} &= \frac{1}{322560} \\
\frac{\partial^1 (c_{2m})}{\partial u^1} &= \frac{u^{5}}{840} \\
\left.\frac{\partial^1 (c_{2m})}{\partial u^1}\right|_{u=\frac{1}{2}} &= \frac{1}{26880} \\
\frac{\partial^2 (c_{2m})}{\partial u^2} &= \frac{u^{4}}{168} \\
\left.\frac{\partial^2 (c_{2m})}{\partial u^2}\right|_{u=\frac{1}{2}} &= \frac{1}{2688}
\end{aligned}
$$
:::

:::example[同じ幅で箱の偶数モーメントを比較（4）]
幅 $2u$、高さ $1/(2u)$ の対称な箱を近似単位に使う。最後にこの値でも計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=4
$$

二乗以降の偶数モーメントは正側を二倍して計算できる。

$$
\begin{aligned}
\int x^{2m}\rho_u(x)\dd x\big|_{m=4} &= \frac{u^{8}}{9}
\end{aligned}
$$

奇数の冪は正側と負側で正負が逆なので消える。

$$
\begin{aligned}
\int x^{2m-1}\rho_u(x)\dd x &= 0
\end{aligned}
$$

試験関数の Taylor 展開に入れると、偶数微分に掛かる重みは階乗でさらに割る。

$$
\begin{aligned}
c_{2m}\big|_{m=4} &= \frac{u^{8}}{362880}
\end{aligned}
$$

観測点 $u$=$\frac{1}{2}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\int x^{2m}\rho_u(x)\dd x\big|_{u=\frac{1}{2}} &= \frac{1}{2304} \\
\frac{\partial^1 (\int x^{2m}\rho_u(x)\dd x)}{\partial u^1} &= \frac{8 u^{7}}{9} \\
\left.\frac{\partial^1 (\int x^{2m}\rho_u(x)\dd x)}{\partial u^1}\right|_{u=\frac{1}{2}} &= \frac{1}{144} \\
\frac{\partial^2 (\int x^{2m}\rho_u(x)\dd x)}{\partial u^2} &= \frac{56 u^{6}}{9} \\
\left.\frac{\partial^2 (\int x^{2m}\rho_u(x)\dd x)}{\partial u^2}\right|_{u=\frac{1}{2}} &= \frac{7}{72}
\end{aligned}
$$

観測点 $u$=$\frac{1}{2}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
c_{2m}\big|_{u=\frac{1}{2}} &= \frac{1}{92897280} \\
\frac{\partial^1 (c_{2m})}{\partial u^1} &= \frac{u^{7}}{45360} \\
\left.\frac{\partial^1 (c_{2m})}{\partial u^1}\right|_{u=\frac{1}{2}} &= \frac{1}{5806080} \\
\frac{\partial^2 (c_{2m})}{\partial u^2} &= \frac{u^{6}}{6480} \\
\left.\frac{\partial^2 (c_{2m})}{\partial u^2}\right|_{u=\frac{1}{2}} &= \frac{1}{414720}
\end{aligned}
$$
:::

### 熱核の近似と弱極限

Gauss 核 $\rho_\epsilon=(4\pi\epsilon)^{-1/2}e^{-x^2/(4\epsilon)}$ は $\int\rho_\epsilon=1$、$\int x^2\rho_\epsilon=2\epsilon$ である。対称性により、十分滑らかで二階導関数が有界な $\varphi$ に対して
$$
\left|\int_\R \rho_\epsilon(x)\varphi(x)\dd x-\varphi(0)\right|
\le \epsilon\sup_{x\in\R}|\varphi''(x)|
$$
となる。実際、$\varphi(x)-\varphi(0)-x\varphi'(0)$ の絶対値は $\sup|\varphi''|\,x^2/2$ であり、一次項の積分は零である。この評価は試験関数に対する誤差であり、$\rho_\epsilon$ の関数値が一様収束することを意味しない。例えば $\rho_\epsilon(0)\to\infty$ であり、原点から離れた点では逆に零へ収束する。
:::example[Gaussian 近似単位の規格化と誤差尺度]
分散 $2$ の核で偶数モーメントと奇数モーメントを比べる。

求める積分を $m_j$ と置く。重みと冪を掛けてから、その変数で積分する。積の形で対称性が分かる場合は先に利用できる。

$$
m_j=\int_{-\infty}^{\infty} x^j\left(\frac{e^{- \frac{x^{2}}{4}}}{2 \sqrt{\pi}}\right)\,\dd x
$$

次数 $0$ の積を整理してから定積分を取る。

$$
m_{0}=\int_{-\infty}^{\infty}\left(\frac{e^{- \frac{x^{2}}{4}}}{2 \sqrt{\pi}}\right)\dd x=1
$$

次数 $1$ では次の積分になる。

$$
m_{1}=\int_{-\infty}^{\infty}\left(\frac{x e^{- \frac{x^{2}}{4}}}{2 \sqrt{\pi}}\right)\dd x=0
$$

次数 $2$ では次の積分になる。

$$
m_{2}=\int_{-\infty}^{\infty}\left(\frac{x^{2} e^{- \frac{x^{2}}{4}}}{2 \sqrt{\pi}}\right)\dd x=2
$$

次数 $3$ では次の積分になる。

$$
m_{3}=\int_{-\infty}^{\infty}\left(\frac{x^{3} e^{- \frac{x^{2}}{4}}}{2 \sqrt{\pi}}\right)\dd x=0
$$

次数 $4$ では次の積分になる。

$$
m_{4}=\int_{-\infty}^{\infty}\left(\frac{x^{4} e^{- \frac{x^{2}}{4}}}{2 \sqrt{\pi}}\right)\dd x=12
$$

次数 $5$ では次の積分になる。

$$
m_{5}=\int_{-\infty}^{\infty}\left(\frac{x^{5} e^{- \frac{x^{2}}{4}}}{2 \sqrt{\pi}}\right)\dd x=0
$$

次数 $6$ では次の積分になる。

$$
m_{6}=\int_{-\infty}^{\infty}\left(\frac{x^{6} e^{- \frac{x^{2}}{4}}}{2 \sqrt{\pi}}\right)\dd x=120
$$
:::

## 試験関数と超関数

### 何に作用するかを指定する

:::definition[超関数]{#def-mp-distribution}
開集合 $U\subset\R^n$ で $\mathcal D(U)=C_c^\infty(U)$ を :term[試験関数]{en="test function" yomi="しけんかんすう"} の空間とする。$\mathcal D(U)$ 上の複素線形汎関数 $T$ を考え、作用 $\langle T,\varphi\rangle$ では共役を取らない。各コンパクト集合 $K\subset U$ に対し、$K$ に台を持つ全ての $\varphi$ について
$$
|\langle T,\varphi\rangle|\le C_K\max_{|\alpha|\le m_K}\sup_K|\partial^\alpha\varphi|
$$
となる定数 $C_K\ge0$ と非負整数 $m_K$ が存在するとき、$T$ を :term[超関数]{en="distribution" yomi="ちょうかんすう"} という。
:::

局所可積分な $f$ は $\langle T_f,\varphi\rangle=\int f\varphi$ によって超関数になり、上式で $m_K=0$ が使える。$\delta_a$ も $m_K=0$ であるが局所可積分な関数からは得られない。実際、$a$ で $1$、小球外で $0$ の $0\le\varphi\le1$ を選ぶと、可積分な関数の場合の値は積分の絶対連続性で零に近づく。滑らかな $g$ との積は $\langle gT,\varphi\rangle=\langle T,g\varphi\rangle$ であり、積の微分の公式で連続性の評価を得る。たとえば $x\delta_2=2\delta_2$ である。二つの任意の超関数の積はこの定義からは作れない。

### 台で点の寄与を特定する

超関数の台とは、超関数が零に制限される開集合をすべて合わせた部分の補集合である。例えば $\delta'_a$ の台は $\{a\}$ である。$a$ を避ける試験関数では値は零だが、$a$ の任意に小さい近くで $\varphi'(a)\ne0$ のものを取れるからである。$\delta_a$ と $\delta'_a$ を同じ集中度と考えると作用の違いを見落とす。例えば $\varphi(x)=(x-a)\eta(x)$、$\eta(a)=1$ とすれば前者は零、後者は $-1$ である。

さらに次数 $m$ の評価を満たし台が $\{a\}$ に含まれる超関数は $\delta_a,\ldots,\delta^{(m)}_a$ の線形結合である。これを示すには全ての $m$ 階までの Taylor 係数が零の $\varphi$ を取る。$a$ の近くで $1$、幅 $O(\epsilon)$ の外で零となる $\chi_\epsilon$ を掛けても値は変わらない。積の微分で $j\le m$ 階の各項は $\chi_\epsilon$ の微分の $\epsilon$ の負冪と、Taylor の剰余 $O(\epsilon^{m+1-(j-\ell)})$ の積になり、全体で $O(\epsilon^{m+1-j})$ である。従って値は零となる。残る作用は $\varphi^{(j)}(a)$ のみで決まり、各評価は $(-1)^j\delta_a^{(j)}$ である。

### 台と特異台

超関数の :term[特異台]{en="singular support" yomi="とくいだい"} は、滑らかな関数で表せる最大の開集合の補集合である。例えば階段関数 $H=\mathbf1_{(0,\infty)}$ の台は $[0,\infty)$ で、特異台は $\{0\}$ である。原点を含む区間に滑らかな代表元があるとすれば、区間の負側では恒等的に零、正側では恒等的に $1$ となる。滑らかな二つの関数が同じ超関数を表すときは連続性により点ごとに一致するためである。ところが左右の値が異なるので原点で連続にならない。$\delta_0$ は原点の外で零であり、原点の近くで局所可積分な代表元を持たないことから、台と特異台がともに $\{0\}$ になる。

:::example[デルタ微分と乗算（1）]
係数 $x^2$ を $\delta_1^{(m)}$ に掛ける。積の微分を使い、高い階数から順に係数を求める。基準の値で計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

最高次には $x^2$ の点での値が掛かる。

$$
\begin{aligned}
[\delta_1^{(m)}] &= 1
\end{aligned}
$$

一段低い項では一度係数を微分し、符号が反転する。

$$
\begin{aligned}
[\delta_1^{(m-1)}]\big|_{m=1} &= -2
\end{aligned}
$$

$x^2$ を一度だけ微分するので、これで $x^2\delta'_1=\delta'_1-2\delta_1$ の全項がそろう。
:::

:::example[デルタ微分と乗算（2）]
係数 $x^2$ を $\delta_1^{(m)}$ に掛ける。積の微分を使い、高い階数から順に係数を求める。値を変えて依存性を調べる。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

最高次には $x^2$ の点での値が掛かる。

$$
\begin{aligned}
[\delta_1^{(m)}] &= 1
\end{aligned}
$$

一段低い項では一度係数を微分し、符号が反転する。

$$
\begin{aligned}
[\delta_1^{(m-1)}]\big|_{m=2} &= -4
\end{aligned}
$$

二段低い項には $m(m-1)$ が掛かる。$m=1$ ではその項自体を省く。

$$
\begin{aligned}
[\delta_1^{(m-2)}]\big|_{m=2} &= 2
\end{aligned}
$$
:::

:::example[デルタ微分と乗算（3）]
係数 $x^2$ を $\delta_1^{(m)}$ に掛ける。積の微分を使い、高い階数から順に係数を求める。さらに値を変えて比較する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

最高次には $x^2$ の点での値が掛かる。

$$
\begin{aligned}
[\delta_1^{(m)}] &= 1
\end{aligned}
$$

一段低い項では一度係数を微分し、符号が反転する。

$$
\begin{aligned}
[\delta_1^{(m-1)}]\big|_{m=3} &= -6
\end{aligned}
$$

二段低い項には $m(m-1)$ が掛かる。$m=1$ ではその項自体を省く。

$$
\begin{aligned}
[\delta_1^{(m-2)}]\big|_{m=3} &= 6
\end{aligned}
$$
:::

:::example[デルタ微分と乗算（4）]
係数 $x^2$ を $\delta_1^{(m)}$ に掛ける。積の微分を使い、高い階数から順に係数を求める。最後にこの値でも計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=4
$$

最高次には $x^2$ の点での値が掛かる。

$$
\begin{aligned}
[\delta_1^{(m)}] &= 1
\end{aligned}
$$

一段低い項では一度係数を微分し、符号が反転する。

$$
\begin{aligned}
[\delta_1^{(m-1)}]\big|_{m=4} &= -8
\end{aligned}
$$

二段低い項には $m(m-1)$ が掛かる。$m=1$ ではその項自体を省く。

$$
\begin{aligned}
[\delta_1^{(m-2)}]\big|_{m=4} &= 12
\end{aligned}
$$
:::

## 超関数の微分と階段関数

### 部分積分を定義にする

$\langle T',\varphi\rangle=-\langle T,\varphi'\rangle$ とする。符号は通常の関数の部分積分と一致する。$H(x)=\mathbf1_{x>0}$ では $-\int_0^\infty\varphi'=\varphi(0)$ なので $H'=\delta_0$ となる。有限個の分割点を除いて $C^1$ 級で、各分割点まで両側から $C^1$ 級に延びる $f$ に対し、各区間で部分積分すると
$$
T_f'=T_{f'_{\mathrm{reg}}}+\sum_j[f(a_j+)-f(a_j-)]\delta_{a_j}
\label{eq-mp-distribution-jump}
$$
となる。値 $f(a_j)$ 自体は積分に影響しない。例えば $(|x|)''=2\delta_0$ である。さらに $\langle \delta'_0,x\varphi\rangle=-\varphi(0)$ から $x\delta'_0=-\delta_0$ である。一般の滑らかな $g$ に対しては、試験関数上の積の微分を展開すると $(gT)'=g'T+gT'$ が従う。

### 境界での跳びが表面源になる

滑らかな境界 $\Sigma$ を挟み、法線 $n$ を内側 $\Omega_-$ から外側 $\Omega_+$ に向ける。ベクトル場が各側で境界まで $C^1$ 級に延び、その境界値が $A_\pm$ なら、超関数の発散は区分的な発散に $(A_+-A_-)\cdot n$ を面測度に掛けたものを加える。定義 $-\int\vb A\cdot\grad\varphi$ に対して各領域で部分積分すると、内側の法線は $n$、外側の内境界の法線は $-n$ であり、境界項の和は $\int_\Sigma(A_+-A_-)\cdot n\,\varphi\dd S$ になる。例えば $r<a$ で $\vb A=\vb r$、外で零なら、内部の発散は $3$、表面の密度は $-a$ である。全空間積分は $3(4\pi a^3/3)-a(4\pi a^2)=0$ となり、遠方に流束がないことと整合する。

### 複数回微分したときの跳び

区分的に滑らかな $f$ を点 $0$ の両側で部分積分すると、超関数の微分は区分微分と $[f]_0\delta_0$ の和になる。この式をもう一度使うと、二階微分には $[f^\prime]_0\delta_0+[f]_0\delta_0^\prime$ が加わる。偶奇性から特異項を推測したあとで、試験関数に作用させると符号を確かめられる。

:::example[正部分の冪の微分（1）]
正の半直線で $x^m$、負側で零の関数を順に微分する。基準の値で計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

正側の $m$ 階微分は階乗であり、負側は零である。

$$
\begin{aligned}
m\text{ 階微分の }H(x)\text{ の係数}\big|_{m=1} &= 1
\end{aligned}
$$

ここまでに境界項が出ないのは、階数 $m-1$ までの原点での値が両側とも零のためである。

$$
\begin{aligned}
[D^{m-1}x_+^m]_0 &= 0
\end{aligned}
$$

さらに一度微分すると階段関数がデルタになる。

$$
\begin{aligned}
[\delta_0]\partial^{m+1}x_+^m\big|_{m=1} &= 1
\end{aligned}
$$
:::

:::example[正部分の冪の微分（2）]
正の半直線で $x^m$、負側で零の関数を順に微分する。値を変えて依存性を調べる。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

正側の $m$ 階微分は階乗であり、負側は零である。

$$
\begin{aligned}
m\text{ 階微分の }H(x)\text{ の係数}\big|_{m=2} &= 2
\end{aligned}
$$

ここまでに境界項が出ないのは、階数 $m-1$ までの原点での値が両側とも零のためである。

$$
\begin{aligned}
[D^{m-1}x_+^m]_0 &= 0
\end{aligned}
$$

さらに一度微分すると階段関数がデルタになる。

$$
\begin{aligned}
[\delta_0]\partial^{m+1}x_+^m\big|_{m=2} &= 2
\end{aligned}
$$
:::

:::example[正部分の冪の微分（3）]
正の半直線で $x^m$、負側で零の関数を順に微分する。さらに値を変えて比較する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

正側の $m$ 階微分は階乗であり、負側は零である。

$$
\begin{aligned}
m\text{ 階微分の }H(x)\text{ の係数}\big|_{m=3} &= 6
\end{aligned}
$$

ここまでに境界項が出ないのは、階数 $m-1$ までの原点での値が両側とも零のためである。

$$
\begin{aligned}
[D^{m-1}x_+^m]_0 &= 0
\end{aligned}
$$

さらに一度微分すると階段関数がデルタになる。

$$
\begin{aligned}
[\delta_0]\partial^{m+1}x_+^m\big|_{m=3} &= 6
\end{aligned}
$$
:::

:::example[正部分の冪の微分（4）]
正の半直線で $x^m$、負側で零の関数を順に微分する。最後にこの値でも計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=4
$$

正側の $m$ 階微分は階乗であり、負側は零である。

$$
\begin{aligned}
m\text{ 階微分の }H(x)\text{ の係数}\big|_{m=4} &= 24
\end{aligned}
$$

ここまでに境界項が出ないのは、階数 $m-1$ までの原点での値が両側とも零のためである。

$$
\begin{aligned}
[D^{m-1}x_+^m]_0 &= 0
\end{aligned}
$$

さらに一度微分すると階段関数がデルタになる。

$$
\begin{aligned}
[\delta_0]\partial^{m+1}x_+^m\big|_{m=4} &= 24
\end{aligned}
$$
:::

### 跳びの位置を動かした場合

折れ線 $f(x)=|x-a|-|x-b|$、$a<b$ を考える。$x<a$ と $x>b$ では傾きが零であり、内側では傾きが $2$ となる。関数は両端で連続なので一階微分に点質量はない。一方、傾きは $a$ で $2$ 増え、$b$ で $2$ 減る。従って $f''=2\delta_a-2\delta_b$ である。例えば $a=1,b=3$ なら試験関数への作用は $2\varphi(1)-2\varphi(3)$ になる。この二つの係数が相殺するため、全実軸で一階微分の両端の値が同じになることとも整合する。

## 緩増加超関数の Fourier 変換

### 急減少する関数を試験する

:term[シュワルツ空間]{en="Schwartz space" yomi="しゅわるつくうかん"} $\mathcal S(\R)$ は、任意の $m,n\ge0$ に対し $\sup|x^m\varphi^{(n)}|$ が有限な滑らかな関数の空間である。これらの有限個の半ノルムの和で絶対値が抑えられる線形汎関数を緩増加超関数という。多項式程度の増大の関数および $\delta_a$ はその例である。Fourier 変換は :eqref[eq-mp-fourier-convention] の規約で、$\langle\widehat T,\varphi\rangle=\langle T,\widehat\varphi\rangle$ と定める。複素共役を掛けない双対性では積分の指数が対称であることに対応する。$\mathcal S$ が保たれることは、積分内微分で出る冪と部分積分で出る冪を可積分関数のノルムで抑えれば分かる。

定義から $\widehat\delta_a=e^{-ika}$、逆変換の公式から $\widehat1=2\pi\delta_0$ である。また $T'$ の変換は $ik\widehat T$ である。例えば $\widehat{\delta'_2}=ik e^{-2ik}$ となる。

### 周期的な点源

$T>0$ とし、$\mathrm{III}_T=\sum_{n\in\Z}\delta_{nT}$ を周期的なデルタ列とする。Schwartz 関数の値は急減少するので作用の和は絶対収束する。Poisson の和公式 :eqref[eq-mp-poisson-sum] を $\widehat\varphi$ に適用し、$\widehat{\widehat\varphi}(x)=2\pi\varphi(-x)$ を使うと
$$
\widehat{\mathrm{III}_T}=\frac{2\pi}{T}\sum_{m\in\Z}\delta_{2\pi m/T}
$$
となる。例えば $T=2$ の列は周波数側で間隔 $\pi$、重み $\pi$ となる。間隔を半分にすると周波数側の点の間隔と重みはともに二倍になる。

### 超関数と平行移動の位相

Schwartz 試験関数に対する Fourier 変換は空間を保つ。定義を双対で延長すると、点質量は位相因子に移り、点での微分は $ik$ の冪を掛ける。点を実軸上で複数並べると位相因子が干渉する。同じ規約で逆変換すれば $1/(2\pi)$ が必ず周波数側に入る。

:::example[反対符号の点源を二つ置く（1）]
源を $\delta_m-\delta_{-m}$ とする。基準の値で計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

各点の位相を差し引く。

$$
\begin{aligned}
\widehat{T}(k)\big|_{m=1} &= - 2 i \sin{\left(k \right)}
\end{aligned}
$$

原点の値が零になるのは全質量が零であるためである。

$$
\begin{aligned}
\widehat T(0) &= 0
\end{aligned}
$$

空間微分を行うと両方の点源の微分が現れ、変換側では $ik$ が掛かる。

$$
\begin{aligned}
\widehat{\partial_x T}(k)\big|_{m=1} &= 2 k \sin{\left(k \right)}
\end{aligned}
$$

観測点 $k$=$\frac{\pi}{4}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\widehat{T}(k)\big|_{k=\frac{\pi}{4}} &= - \sqrt{2} i \\
\frac{\partial^1 (\widehat{T}(k))}{\partial k^1} &= - 2 i \cos{\left(k \right)} \\
\left.\frac{\partial^1 (\widehat{T}(k))}{\partial k^1}\right|_{k=\frac{\pi}{4}} &= - \sqrt{2} i \\
\frac{\partial^2 (\widehat{T}(k))}{\partial k^2} &= 2 i \sin{\left(k \right)} \\
\left.\frac{\partial^2 (\widehat{T}(k))}{\partial k^2}\right|_{k=\frac{\pi}{4}} &= \sqrt{2} i
\end{aligned}
$$

観測点 $k$=$\frac{\pi}{4}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\widehat{\partial_x T}(k)\big|_{k=\frac{\pi}{4}} &= \frac{\sqrt{2} \pi}{4} \\
\frac{\partial^1 (\widehat{\partial_x T}(k))}{\partial k^1} &= 2 k \cos{\left(k \right)} + 2 \sin{\left(k \right)} \\
\left.\frac{\partial^1 (\widehat{\partial_x T}(k))}{\partial k^1}\right|_{k=\frac{\pi}{4}} &= \frac{\sqrt{2} \left(\pi + 4\right)}{4} \\
\frac{\partial^2 (\widehat{\partial_x T}(k))}{\partial k^2} &= - 2 k \sin{\left(k \right)} + 4 \cos{\left(k \right)} \\
\left.\frac{\partial^2 (\widehat{\partial_x T}(k))}{\partial k^2}\right|_{k=\frac{\pi}{4}} &= \frac{\sqrt{2} \left(8 - \pi\right)}{4}
\end{aligned}
$$
:::

:::example[反対符号の点源を二つ置く（2）]
源を $\delta_m-\delta_{-m}$ とする。値を変えて依存性を調べる。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

各点の位相を差し引く。

$$
\begin{aligned}
\widehat{T}(k)\big|_{m=2} &= - 2 i \sin{\left(2 k \right)}
\end{aligned}
$$

原点の値が零になるのは全質量が零であるためである。

$$
\begin{aligned}
\widehat T(0) &= 0
\end{aligned}
$$

空間微分を行うと両方の点源の微分が現れ、変換側では $ik$ が掛かる。

$$
\begin{aligned}
\widehat{\partial_x T}(k)\big|_{m=2} &= 2 k \sin{\left(2 k \right)}
\end{aligned}
$$

観測点 $k$=$\frac{\pi}{4}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\widehat{T}(k)\big|_{k=\frac{\pi}{4}} &= - 2 i \\
\frac{\partial^1 (\widehat{T}(k))}{\partial k^1} &= - 4 i \cos{\left(2 k \right)} \\
\left.\frac{\partial^1 (\widehat{T}(k))}{\partial k^1}\right|_{k=\frac{\pi}{4}} &= 0 \\
\frac{\partial^2 (\widehat{T}(k))}{\partial k^2} &= 8 i \sin{\left(2 k \right)} \\
\left.\frac{\partial^2 (\widehat{T}(k))}{\partial k^2}\right|_{k=\frac{\pi}{4}} &= 8 i
\end{aligned}
$$

観測点 $k$=$\frac{\pi}{4}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\widehat{\partial_x T}(k)\big|_{k=\frac{\pi}{4}} &= \frac{\pi}{2} \\
\frac{\partial^1 (\widehat{\partial_x T}(k))}{\partial k^1} &= 4 k \cos{\left(2 k \right)} + 2 \sin{\left(2 k \right)} \\
\left.\frac{\partial^1 (\widehat{\partial_x T}(k))}{\partial k^1}\right|_{k=\frac{\pi}{4}} &= 2 \\
\frac{\partial^2 (\widehat{\partial_x T}(k))}{\partial k^2} &= - 8 k \sin{\left(2 k \right)} + 8 \cos{\left(2 k \right)} \\
\left.\frac{\partial^2 (\widehat{\partial_x T}(k))}{\partial k^2}\right|_{k=\frac{\pi}{4}} &= - 2 \pi
\end{aligned}
$$
:::

:::example[反対符号の点源を二つ置く（3）]
源を $\delta_m-\delta_{-m}$ とする。さらに値を変えて比較する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

各点の位相を差し引く。

$$
\begin{aligned}
\widehat{T}(k)\big|_{m=3} &= - 2 i \sin{\left(3 k \right)}
\end{aligned}
$$

原点の値が零になるのは全質量が零であるためである。

$$
\begin{aligned}
\widehat T(0) &= 0
\end{aligned}
$$

空間微分を行うと両方の点源の微分が現れ、変換側では $ik$ が掛かる。

$$
\begin{aligned}
\widehat{\partial_x T}(k)\big|_{m=3} &= 2 k \sin{\left(3 k \right)}
\end{aligned}
$$

観測点 $k$=$\frac{\pi}{4}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\widehat{T}(k)\big|_{k=\frac{\pi}{4}} &= - \sqrt{2} i \\
\frac{\partial^1 (\widehat{T}(k))}{\partial k^1} &= - 6 i \cos{\left(3 k \right)} \\
\left.\frac{\partial^1 (\widehat{T}(k))}{\partial k^1}\right|_{k=\frac{\pi}{4}} &= 3 \sqrt{2} i \\
\frac{\partial^2 (\widehat{T}(k))}{\partial k^2} &= 18 i \sin{\left(3 k \right)} \\
\left.\frac{\partial^2 (\widehat{T}(k))}{\partial k^2}\right|_{k=\frac{\pi}{4}} &= 9 \sqrt{2} i
\end{aligned}
$$

観測点 $k$=$\frac{\pi}{4}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\widehat{\partial_x T}(k)\big|_{k=\frac{\pi}{4}} &= \frac{\sqrt{2} \pi}{4} \\
\frac{\partial^1 (\widehat{\partial_x T}(k))}{\partial k^1} &= 6 k \cos{\left(3 k \right)} + 2 \sin{\left(3 k \right)} \\
\left.\frac{\partial^1 (\widehat{\partial_x T}(k))}{\partial k^1}\right|_{k=\frac{\pi}{4}} &= \frac{\sqrt{2} \left(4 - 3 \pi\right)}{4} \\
\frac{\partial^2 (\widehat{\partial_x T}(k))}{\partial k^2} &= - 18 k \sin{\left(3 k \right)} + 12 \cos{\left(3 k \right)} \\
\left.\frac{\partial^2 (\widehat{\partial_x T}(k))}{\partial k^2}\right|_{k=\frac{\pi}{4}} &= \sqrt{2} \left(- \frac{9 \pi}{4} - 6\right)
\end{aligned}
$$
:::

:::example[反対符号の点源を二つ置く（4）]
源を $\delta_m-\delta_{-m}$ とする。最後にこの値でも計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=4
$$

各点の位相を差し引く。

$$
\begin{aligned}
\widehat{T}(k)\big|_{m=4} &= - 2 i \sin{\left(4 k \right)}
\end{aligned}
$$

原点の値が零になるのは全質量が零であるためである。

$$
\begin{aligned}
\widehat T(0) &= 0
\end{aligned}
$$

空間微分を行うと両方の点源の微分が現れ、変換側では $ik$ が掛かる。

$$
\begin{aligned}
\widehat{\partial_x T}(k)\big|_{m=4} &= 2 k \sin{\left(4 k \right)}
\end{aligned}
$$

観測点 $k$=$\frac{\pi}{4}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\widehat{T}(k)\big|_{k=\frac{\pi}{4}} &= 0 \\
\frac{\partial^1 (\widehat{T}(k))}{\partial k^1} &= - 8 i \cos{\left(4 k \right)} \\
\left.\frac{\partial^1 (\widehat{T}(k))}{\partial k^1}\right|_{k=\frac{\pi}{4}} &= 8 i \\
\frac{\partial^2 (\widehat{T}(k))}{\partial k^2} &= 32 i \sin{\left(4 k \right)} \\
\left.\frac{\partial^2 (\widehat{T}(k))}{\partial k^2}\right|_{k=\frac{\pi}{4}} &= 0
\end{aligned}
$$

観測点 $k$=$\frac{\pi}{4}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\widehat{\partial_x T}(k)\big|_{k=\frac{\pi}{4}} &= 0 \\
\frac{\partial^1 (\widehat{\partial_x T}(k))}{\partial k^1} &= 8 k \cos{\left(4 k \right)} + 2 \sin{\left(4 k \right)} \\
\left.\frac{\partial^1 (\widehat{\partial_x T}(k))}{\partial k^1}\right|_{k=\frac{\pi}{4}} &= - 2 \pi \\
\frac{\partial^2 (\widehat{\partial_x T}(k))}{\partial k^2} &= - 32 k \sin{\left(4 k \right)} + 16 \cos{\left(4 k \right)} \\
\left.\frac{\partial^2 (\widehat{\partial_x T}(k))}{\partial k^2}\right|_{k=\frac{\pi}{4}} &= -16
\end{aligned}
$$
:::

:::example[変換の原点微分からモーメントを読む]
二点の符号付き質量を 1 と 2 に置き、重みをそれぞれ $1$ と $-1$ にする。

Fourier 変換は平行移動の位相を引いた式になる。

$$
\begin{aligned}
\widehat T(k) &= e^{- i k} - e^{- 2 i k}
\end{aligned}
$$

各導関数を零周波数に戻す。質量と、一次・二次の符号付きモーメントが別々に現れる。

$$
\begin{aligned}
\widehat T(0) &= 0 \\
\widehat T^\prime(0) &= i \\
\widehat T^{\prime\prime}(0) &= 3
\end{aligned}
$$

緩増加超関数としては、空間の二階微分で $-k^2$ が掛かる。

$$
\begin{aligned}
\widehat{D^2T}(k) &= - k^{2} \left(e^{- i k} - e^{- 2 i k}\right)
\end{aligned}
$$
:::


:::example[変換の原点微分からモーメントを読む]
二点の符号付き質量を -2 と 1 に置き、重みをそれぞれ $1$ と $-1$ にする。

Fourier 変換は平行移動の位相を引いた式になる。

$$
\begin{aligned}
\widehat T(k) &= e^{2 i k} - e^{- i k}
\end{aligned}
$$

各導関数を零周波数に戻す。質量と、一次・二次の符号付きモーメントが別々に現れる。

$$
\begin{aligned}
\widehat T(0) &= 0 \\
\widehat T^\prime(0) &= 3 i \\
\widehat T^{\prime\prime}(0) &= -3
\end{aligned}
$$

緩増加超関数としては、空間の二階微分で $-k^2$ が掛かる。

$$
\begin{aligned}
\widehat{D^2T}(k) &= - k^{2} \left(e^{2 i k} - e^{- i k}\right)
\end{aligned}
$$
:::


:::example[変換の原点微分からモーメントを読む]
二点の符号付き質量を -2 と 2 に置き、重みをそれぞれ $1$ と $-1$ にする。

Fourier 変換は平行移動の位相を引いた式になる。

$$
\begin{aligned}
\widehat T(k) &= e^{2 i k} - e^{- 2 i k}
\end{aligned}
$$

各導関数を零周波数に戻す。質量と、一次・二次の符号付きモーメントが別々に現れる。

$$
\begin{aligned}
\widehat T(0) &= 0 \\
\widehat T^\prime(0) &= 4 i \\
\widehat T^{\prime\prime}(0) &= 0
\end{aligned}
$$

緩増加超関数としては、空間の二階微分で $-k^2$ が掛かる。

$$
\begin{aligned}
\widehat{D^2T}(k) &= - k^{2} \left(e^{2 i k} - e^{- 2 i k}\right)
\end{aligned}
$$
:::

## 主値と有限部分

### 特異な積分から発散分を差し引く

$\operatorname{PV}(1/x)$ は $\int_{|x|>\epsilon}\varphi(x)/x\,\dd x$ の極限である。原点の近くで偶なカットオフを掛けた $\varphi(0)$ を引くと、残りの比は有界である。二乗の特異性では偶部分が残るため、:term[有限部分]{en="finite part" yomi="ゆうげんぶぶん"} を
$$
\left\langle\operatorname{Fp}\frac1{x^2},\varphi\right\rangle
=\lim_{\epsilon\downarrow0}\left[\int_{|x|>\epsilon}\frac{\varphi(x)}{x^2}\dd x-\frac{2\varphi(0)}{\epsilon}\right]
$$
と定める。部分積分すると $-(\varphi(\epsilon)+\varphi(-\epsilon))/\epsilon$ が $\langle-(\operatorname{PV}1/x)',\varphi\rangle$ の境界項になる。二階までの Taylor 展開で差が $O(\epsilon)$ なので、有限部分は $-(\operatorname{PV}1/x)'$ に等しい。

式 :eqref[eq-mp-plemelj] は $(x\pm i0)^{-1}=\operatorname{PV}(1/x)\mp i\pi\delta_0$ と読める。$H(x)e^{-\epsilon x}$ の Fourier 変換 $(\epsilon+ik)^{-1}$ にこの式を適用すると $\widehat H=\pi\delta_0-i\operatorname{PV}(1/k)$ である。指数を除く極限は試験関数の優収束で得る。これにより $ik\widehat H=1$ となり $H'=\delta_0$ を確かめられる。

### 主値の Fourier 変換と片側の周波数

$\sgn x=2H(x)-1$ の変換は $-2i\operatorname{PV}(1/k)$ である。全変換を二回適用すると $\widehat{\widehat T}(x)=2\pi T(-x)$ だから、奇性を使って
$$
\widehat{\operatorname{PV}(1/x)}=-i\pi\sgn k
$$
を得る。従って有限部分 $\operatorname{Fp}(1/x^2)$ の変換は $-ik(-i\pi\sgn k)=-\pi|k|$ である。式 :eqref[eq-mp-plemelj] に戻すと $(x+i0)^{-1}$ の変換は $-i\pi(\sgn k+1)$ で、正の周波数側だけに現れる。符号の検算には $\epsilon>0$ の $x+i\epsilon$ が下半平面に極を持ち、負指数 Fourier 変換で正の $k$ のときだけ下側を閉じることが使える。

### 有限部分の高階化

$\operatorname{Fp}(x^{-m})$ は $\operatorname{PV}(1/x)$ を $m-1$ 回微分して $(-1)^{m-1}/(m-1)!$ 倍したものと定められる。この規約では乗算で $m\ge2$ なら $x\operatorname{Fp}(x^{-m})=\operatorname{Fp}(x^{1-m})$ が成り立つ。$m=1$ の出発点は $x\operatorname{PV}(1/x)=1$ である。証明は $xT^\prime=(xT)^\prime-T$ を繰り返すことで得られる。特異点を含む通常の積分をそのまま分けると各項が発散する。先に同じ試験関数に作用させ、差し引きと境界の極限を合わせる。

同じ規約では高階の境界値も
$$
(x\pm i0)^{-m}
=\operatorname{Fp}(x^{-m})\mp
\frac{i\pi(-1)^{m-1}}{(m-1)!}\delta_0^{(m-1)}
$$
となる。$\epsilon>0$ で $(x\pm i\epsilon)^{-1}$ を $m-1$ 回微分すると $(-1)^{m-1}(m-1)!(x\pm i\epsilon)^{-m}$ になる。超関数として微分した極限は、試験関数をその導関数に置き換えて一次の境界値公式を使えばよい。例えば二次の境界値は $\operatorname{Fp}(x^{-2})\pm i\pi\delta'_0$ である。そのため試験関数の値だけでなく傾きも、極を避ける側によって効いてくる。

:::example[有限部分の Fourier 乗数（1）]
一階の主値の変換から高階に移す。正の $k$ 側を調べる。基準の値で計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

主値の変換は $-i\pi\operatorname{sgn}k$ である。正側では符号が $+1$ になる。

$$
\begin{aligned}
\widehat{\operatorname{PV}(1/x)}(k) &= - i \pi
\end{aligned}
$$

高階に移す規格化を求める。

$$
\begin{aligned}
(-1)^{m-1}/(m-1)!\big|_{m=1} &= 1
\end{aligned}
$$

微分の乗数を掛けて整理する。

$$
\begin{aligned}
\widehat{\operatorname{Fp}(x^{-m})}(k)\big|_{m=1} &= - i \pi
\end{aligned}
$$
:::

:::example[有限部分の Fourier 乗数（2）]
一階の主値の変換から高階に移す。正の $k$ 側を調べる。値を変えて依存性を調べる。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

主値の変換は $-i\pi\operatorname{sgn}k$ である。正側では符号が $+1$ になる。

$$
\begin{aligned}
\widehat{\operatorname{PV}(1/x)}(k) &= - i \pi
\end{aligned}
$$

高階に移す規格化を求める。

$$
\begin{aligned}
(-1)^{m-1}/(m-1)!\big|_{m=2} &= -1
\end{aligned}
$$

微分の乗数を掛けて整理する。

$$
\begin{aligned}
\widehat{\operatorname{Fp}(x^{-m})}(k)\big|_{m=2} &= - \pi k
\end{aligned}
$$

観測点 $k$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\widehat{\operatorname{Fp}(x^{-m})}(k)\big|_{k=2} &= - 2 \pi \\
\frac{\partial^1 (\widehat{\operatorname{Fp}(x^{-m})}(k))}{\partial k^1} &= - \pi \\
\left.\frac{\partial^1 (\widehat{\operatorname{Fp}(x^{-m})}(k))}{\partial k^1}\right|_{k=2} &= - \pi
\end{aligned}
$$
:::

:::example[有限部分の Fourier 乗数（3）]
一階の主値の変換から高階に移す。正の $k$ 側を調べる。さらに値を変えて比較する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

主値の変換は $-i\pi\operatorname{sgn}k$ である。正側では符号が $+1$ になる。

$$
\begin{aligned}
\widehat{\operatorname{PV}(1/x)}(k) &= - i \pi
\end{aligned}
$$

高階に移す規格化を求める。

$$
\begin{aligned}
(-1)^{m-1}/(m-1)!\big|_{m=3} &= \frac{1}{2}
\end{aligned}
$$

微分の乗数を掛けて整理する。

$$
\begin{aligned}
\widehat{\operatorname{Fp}(x^{-m})}(k)\big|_{m=3} &= \frac{i \pi k^{2}}{2}
\end{aligned}
$$

観測点 $k$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\widehat{\operatorname{Fp}(x^{-m})}(k)\big|_{k=2} &= 2 i \pi \\
\frac{\partial^1 (\widehat{\operatorname{Fp}(x^{-m})}(k))}{\partial k^1} &= i \pi k \\
\left.\frac{\partial^1 (\widehat{\operatorname{Fp}(x^{-m})}(k))}{\partial k^1}\right|_{k=2} &= 2 i \pi
\end{aligned}
$$
:::

:::example[有限部分の Fourier 乗数（4）]
一階の主値の変換から高階に移す。正の $k$ 側を調べる。最後にこの値でも計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=4
$$

主値の変換は $-i\pi\operatorname{sgn}k$ である。正側では符号が $+1$ になる。

$$
\begin{aligned}
\widehat{\operatorname{PV}(1/x)}(k) &= - i \pi
\end{aligned}
$$

高階に移す規格化を求める。

$$
\begin{aligned}
(-1)^{m-1}/(m-1)!\big|_{m=4} &= - \frac{1}{6}
\end{aligned}
$$

微分の乗数を掛けて整理する。

$$
\begin{aligned}
\widehat{\operatorname{Fp}(x^{-m})}(k)\big|_{m=4} &= \frac{\pi k^{3}}{6}
\end{aligned}
$$

観測点 $k$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\widehat{\operatorname{Fp}(x^{-m})}(k)\big|_{k=2} &= \frac{4 \pi}{3} \\
\frac{\partial^1 (\widehat{\operatorname{Fp}(x^{-m})}(k))}{\partial k^1} &= \frac{\pi k^{2}}{2} \\
\left.\frac{\partial^1 (\widehat{\operatorname{Fp}(x^{-m})}(k))}{\partial k^1}\right|_{k=2} &= 2 \pi
\end{aligned}
$$
:::

### 積の規約を固定する

例えば $\delta_0$ に連続な関数を掛けるなら、点の値を用いて $f(0)\delta_0$ と定義できる。しかし 関数 $f$ が $0$ で微分可能なら、テイラー係数から $f\delta'_0=f(0)\delta'_0-f'(0)\delta_0$ と定められる。さらに $H\delta_0$ のように双方が原点で特異なら、どちらをどのように近似したかで極限が変わる。$H_\epsilon(x)=\int_{-\infty}^x\rho_\epsilon(t)\dd t$ とし、$\rho_\epsilon$ を非負で積分が $1$ の近似単位とする。すると $H_\epsilon\rho_\epsilon=(H_\epsilon^2)'/2$ であり、$0\le H_\epsilon\le1$ と原点外での収束から、試験関数に対して $\delta_0/2$ へ収束する。一方、$\rho_\epsilon$ の台が $[-\epsilon,\epsilon]$ に収まるよう選び、階段を右へ $3\epsilon$ だけ移せば、その台上で $H_\epsilon(x-3\epsilon)=0$ である。左へ同じだけ移せば $1$ となる。どちらの階段も $H$ に収束するが、積の極限はそれぞれ零と $\delta_0$ になる。従ってこの正則化の指定なしに積を公式として使うことはできない。

## 多次元のデルタ関数と座標変換

### 点源を曲線座標に移す

$n$ 次元で $\langle\delta^{(n)}_a,\varphi\rangle=\varphi(a)$ とする。この節の上付き $(n)$ は空間の次元を表す。微分同相 $x=F(q)$、$q_0=F^{-1}(a)$ では
$$
\delta^{(n)}(F(q)-a)=\frac{\delta^{(n)}(q-q_0)}{|\det DF(q_0)|}
\label{eq-mp-delta-jacobian}
$$
である。可積分な近似を左辺に入れて $x$ に変数変換すると、試験関数は $\varphi(F^{-1}(x))/|\det DF(F^{-1}(x))|$ に変わるので証明できる。一次元で $g$ の零点が単純かつ有限個なら、零点近傍で同じ議論、残りで台が零点を避けることから $\delta(g)=\sum_{g(a)=0}\delta_a/|g'(a)|$ となる。重根ではこの証明の逆写像が存在しない。例えば $\delta(x^2-4)=(\delta_2+\delta_{-2})/4$ である。球座標の原点には座標の微分同相が無いので、体積要素の逆数を形式的に原点へ代入することはできない。

### 面に分布する源

滑らかな実関数 $F$ の零集合 $F=0$ で $\grad F\ne0$ なら
$$
\int \delta(F(x))\varphi(x)\dd^nx
=\int_{F=0}\frac{\varphi(x)}{|\grad F(x)|}\dd S
$$
である。局所的に $\partial_{x_n}F\ne0$ として、$x'=(x_1,\ldots,x_{n-1})$ を保ち $s=F(x)$ に変数変換する。$s$ のデルタが $0$ だけを取り、測度の倍率は $1/|\partial_{x_n}F|$ である。零面をグラフ $x_n=g(x')$ とすると、$\partial_i g=-\partial_iF/\partial_{x_n}F$ より $\dd S=|\grad F|\,\dd x'/|\partial_{x_n}F|$ である。試験関数の台を有限の座標近傍に分割して足すと成立する。例えば三次元で $F=r^2-a^2$ なら半径 $a$ の球面上で倍率は $1/(2a)$ であり、全質量は $2\pi a$ となる。半径そのものを $F=r-a$ とした場合は全質量が $4\pi a^2$ になる。

### 面に載せた源と半径の重み

写像 $g:\R^d\to\R$ の零集合で $\grad g\ne0$ とする。零集合の局所座標と $g$ を座標に取ると、面積測度の法線方向の速度は $1/|\grad g|$ になる。そこで $g$ の積分を点質量で実行すると面上の積分が残る。$\delta(g)$ を面の面積密度にしたい場合は、$|\grad g|$ を掛けて規格化する。原点のように $\grad g$ が消える値は別に扱う。

:::example[球殻のデルタ源（1）]
三次元で $g=r^2-m^2$ とする。球殻と内側の点を区別する。基準の値で計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

零集合の半径での法線勾配の大きさを求める。

$$
\begin{aligned}
|\nabla g|_{r=m}\big|_{m=1} &= 2
\end{aligned}
$$

動径積分には $r^2$ が掛かるので、単位試験関数に対する質量が得られる。

$$
\begin{aligned}
\int_{\R^3}\delta(g)\dd^3x\big|_{m=1} &= 2 \pi
\end{aligned}
$$

例えば試験関数が半径の二乗なら、球上で定数になりさらに $m^2$ が掛かる。

$$
\begin{aligned}
\int_{\R^3}r^2\delta(g)\dd^3x\big|_{m=1} &= 2 \pi
\end{aligned}
$$
:::

:::example[球殻のデルタ源（2）]
三次元で $g=r^2-m^2$ とする。球殻と内側の点を区別する。値を変えて依存性を調べる。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

零集合の半径での法線勾配の大きさを求める。

$$
\begin{aligned}
|\nabla g|_{r=m}\big|_{m=2} &= 4
\end{aligned}
$$

動径積分には $r^2$ が掛かるので、単位試験関数に対する質量が得られる。

$$
\begin{aligned}
\int_{\R^3}\delta(g)\dd^3x\big|_{m=2} &= 4 \pi
\end{aligned}
$$

例えば試験関数が半径の二乗なら、球上で定数になりさらに $m^2$ が掛かる。

$$
\begin{aligned}
\int_{\R^3}r^2\delta(g)\dd^3x\big|_{m=2} &= 16 \pi
\end{aligned}
$$
:::

:::example[球殻のデルタ源（3）]
三次元で $g=r^2-m^2$ とする。球殻と内側の点を区別する。さらに値を変えて比較する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

零集合の半径での法線勾配の大きさを求める。

$$
\begin{aligned}
|\nabla g|_{r=m}\big|_{m=3} &= 6
\end{aligned}
$$

動径積分には $r^2$ が掛かるので、単位試験関数に対する質量が得られる。

$$
\begin{aligned}
\int_{\R^3}\delta(g)\dd^3x\big|_{m=3} &= 6 \pi
\end{aligned}
$$

例えば試験関数が半径の二乗なら、球上で定数になりさらに $m^2$ が掛かる。

$$
\begin{aligned}
\int_{\R^3}r^2\delta(g)\dd^3x\big|_{m=3} &= 54 \pi
\end{aligned}
$$
:::

:::example[球殻のデルタ源（4）]
三次元で $g=r^2-m^2$ とする。球殻と内側の点を区別する。最後にこの値でも計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=4
$$

零集合の半径での法線勾配の大きさを求める。

$$
\begin{aligned}
|\nabla g|_{r=m}\big|_{m=4} &= 8
\end{aligned}
$$

動径積分には $r^2$ が掛かるので、単位試験関数に対する質量が得られる。

$$
\begin{aligned}
\int_{\R^3}\delta(g)\dd^3x\big|_{m=4} &= 8 \pi
\end{aligned}
$$

例えば試験関数が半径の二乗なら、球上で定数になりさらに $m^2$ が掛かる。

$$
\begin{aligned}
\int_{\R^3}r^2\delta(g)\dd^3x\big|_{m=4} &= 128 \pi
\end{aligned}
$$
:::

## 超関数としての Green 関数

### 微分作用素に対する点応答

作用素 $L$ の基本解とは $LG=\delta_0$ を満たす超関数である。$L$ が定数係数なら、台がコンパクトな滑らかな $f$ に対し $L(G*f)=(LG)*f=f$ となる。畳み込みを $\langle G_y,f(x-y)\rangle$ と書くと、$x$ 微分は $f$ の $y$ 微分の負に等しく、超関数の微分の定義で $G$ へ移るからである。境界や因果律があれば、それを満たす核を選ぶ必要がある。

:::example[一次元と三次元の点源]{#exm-mp-fundamental-laplace}
一次元で $G=|x|/2$ は $G''=\delta_0$ を満たす。三次元では $G=1/(4\pi r)$ として $-\Delta G=\delta^{(3)}_0$ である。後者を示すには、半径 $\epsilon$ の球を除いて Green の恒等式を使う。領域の内側境界の外向き法線は $-\hat{\vb e}_r$ なので $\partial_nG=1/(4\pi\epsilon^2)$ である。$\langle G,-\Delta\varphi\rangle$ は $\int_{r=\epsilon}(\varphi\partial_n G-G\partial_n\varphi)\dd S$ に収束し、第一項は $\varphi(0)$、第二項は $O(\epsilon)$ で消える。除いた球での $G\Delta\varphi$ の積分も $O(\epsilon^2)$ で消える。
:::

二次元で $-\Delta$ の基本解は $G(r)=-\log(r/r_0)/(2\pi)$ である。$r_0>0$ は対数の引数を無次元にする基準長さである。原点を除けば $(rG')'=0$ であり、原点近くでは $r|\log(r/r_0)|$ が可積分なので、$G$ は局所可積分関数として超関数を定める。円を除いて Green の恒等式を使うと、内側の法線微分が $1/(2\pi\epsilon)$ になり、$\int_{r=\epsilon}\varphi\partial_nG\,\dd S$ は円周平均となる。これは $\varphi(0)$ へ収束し、もう一方の境界項は $O(\epsilon|\log\epsilon|)$ で消える。従って $-\Delta G=\delta^{(2)}_0$ である。$r_0$ の変更は定数を足す操作になり、点源の強さには影響しない。

### 因果条件のある点応答

時間の調和振動子 $L=\partial_t^2+\omega_0^2$（$\omega_0>0$）で、時刻零の衝撃に対する遅延基本解を考える。$G(t)=H(t)\sin(\omega_0t)/\omega_0$ は $t<0$ で零であり、$t>0$ では同次方程式を満たす。原点で連続なので $\delta'_0$ の項は現れず、導関数の跳びが $1$ なので $LG=\delta_0$ である。

台が $[0,\infty)$ に含まれる基本解が二つあるとし、その差を $u$ とする。$L=(\partial_t+i\omega_0)(\partial_t-i\omega_0)$ と因数分解し、$v=(\partial_t-i\omega_0)u$ と置くと $(e^{i\omega_0t}v)'=0$ となる。微分が零の超関数は定数である。実際、積分が零の試験関数はコンパクト台を持つ原始関数の微分なので、その超関数の作用が零になる。従って作用は試験関数の積分だけに比例する。$v$ は負の時刻で零であるため定数も零になる。さらに $(e^{-i\omega_0t}u)'=0$ に同じ議論を使うと $u=0$ となり、遅延基本解の一意性が分かる。

単位質量の振動子に時刻 $a$ で速度を $I$ だけ変える衝撃が加わると、静止状態からの解は $I G(t-a)$ になる。$\omega_0=2\,\mathrm{s^{-1}}$、$I=3\,\mathrm{m\,s^{-1}}$ なら、衝撃の後の振幅は $1.5\,\mathrm m$ である。因果条件で選んだこの核は $t$ に関する減衰を仮定しておらず、振動は持続する。

### 二階作用素の跳びと減衰条件

一次元の演算子 $L=-\partial_x^2+m^2$（$m>0$）の点源解を作る。左右の同次解のうち無限遠で減衰する枝を選び、原点で同じ値にする。導関数の跳びは $[G^\prime]_0=-1$ でなければならない。すると原点の外では $LG=0$、原点ではデルタの係数が $1$ になる。Fourier 変換で $(k^2+m^2)\widehat G=1$ としても同じ規格化に達する。

:::example[減衰核の積分と幅（1）]
正側の指数の係数を導関数の跳びから求める。基準の値で計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

正側の半直線に減衰する解を置く。

$$
\begin{aligned}
G_+(x)\big|_{m=1} &= \frac{e^{- x}}{2}
\end{aligned}
$$

両側で積分すると、定数右辺に対する応答が得られる。

$$
\begin{aligned}
\int_\R G(x)\dd x\big|_{m=1} &= 1
\end{aligned}
$$

Fourier 側での逆作用素の乗数を求める。

$$
\begin{aligned}
\widehat G(k)\big|_{m=1} &= \frac{1}{k^{2} + 1}
\end{aligned}
$$

観測点 $x$=$1$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
G_+(x)\big|_{x=1} &= \frac{1}{2 e} \\
\frac{\partial^1 (G_+(x))}{\partial x^1} &= - \frac{e^{- x}}{2} \\
\left.\frac{\partial^1 (G_+(x))}{\partial x^1}\right|_{x=1} &= - \frac{1}{2 e} \\
\frac{\partial^2 (G_+(x))}{\partial x^2} &= \frac{e^{- x}}{2} \\
\left.\frac{\partial^2 (G_+(x))}{\partial x^2}\right|_{x=1} &= \frac{1}{2 e}
\end{aligned}
$$
:::

:::example[減衰核の積分と幅（2）]
正側の指数の係数を導関数の跳びから求める。値を変えて依存性を調べる。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

正側の半直線に減衰する解を置く。

$$
\begin{aligned}
G_+(x)\big|_{m=2} &= \frac{e^{- 2 x}}{4}
\end{aligned}
$$

両側で積分すると、定数右辺に対する応答が得られる。

$$
\begin{aligned}
\int_\R G(x)\dd x\big|_{m=2} &= \frac{1}{4}
\end{aligned}
$$

Fourier 側での逆作用素の乗数を求める。

$$
\begin{aligned}
\widehat G(k)\big|_{m=2} &= \frac{1}{k^{2} + 4}
\end{aligned}
$$

観測点 $x$=$1$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
G_+(x)\big|_{x=1} &= \frac{1}{4 e^{2}} \\
\frac{\partial^1 (G_+(x))}{\partial x^1} &= - \frac{e^{- 2 x}}{2} \\
\left.\frac{\partial^1 (G_+(x))}{\partial x^1}\right|_{x=1} &= - \frac{1}{2 e^{2}} \\
\frac{\partial^2 (G_+(x))}{\partial x^2} &= e^{- 2 x} \\
\left.\frac{\partial^2 (G_+(x))}{\partial x^2}\right|_{x=1} &= e^{-2}
\end{aligned}
$$
:::

:::example[減衰核の積分と幅（3）]
正側の指数の係数を導関数の跳びから求める。さらに値を変えて比較する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

正側の半直線に減衰する解を置く。

$$
\begin{aligned}
G_+(x)\big|_{m=3} &= \frac{e^{- 3 x}}{6}
\end{aligned}
$$

両側で積分すると、定数右辺に対する応答が得られる。

$$
\begin{aligned}
\int_\R G(x)\dd x\big|_{m=3} &= \frac{1}{9}
\end{aligned}
$$

Fourier 側での逆作用素の乗数を求める。

$$
\begin{aligned}
\widehat G(k)\big|_{m=3} &= \frac{1}{k^{2} + 9}
\end{aligned}
$$

観測点 $x$=$1$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
G_+(x)\big|_{x=1} &= \frac{1}{6 e^{3}} \\
\frac{\partial^1 (G_+(x))}{\partial x^1} &= - \frac{e^{- 3 x}}{2} \\
\left.\frac{\partial^1 (G_+(x))}{\partial x^1}\right|_{x=1} &= - \frac{1}{2 e^{3}} \\
\frac{\partial^2 (G_+(x))}{\partial x^2} &= \frac{3 e^{- 3 x}}{2} \\
\left.\frac{\partial^2 (G_+(x))}{\partial x^2}\right|_{x=1} &= \frac{3}{2 e^{3}}
\end{aligned}
$$
:::

:::example[減衰核の積分と幅（4）]
正側の指数の係数を導関数の跳びから求める。最後にこの値でも計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=4
$$

正側の半直線に減衰する解を置く。

$$
\begin{aligned}
G_+(x)\big|_{m=4} &= \frac{e^{- 4 x}}{8}
\end{aligned}
$$

両側で積分すると、定数右辺に対する応答が得られる。

$$
\begin{aligned}
\int_\R G(x)\dd x\big|_{m=4} &= \frac{1}{16}
\end{aligned}
$$

Fourier 側での逆作用素の乗数を求める。

$$
\begin{aligned}
\widehat G(k)\big|_{m=4} &= \frac{1}{k^{2} + 16}
\end{aligned}
$$

観測点 $x$=$1$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
G_+(x)\big|_{x=1} &= \frac{1}{8 e^{4}} \\
\frac{\partial^1 (G_+(x))}{\partial x^1} &= - \frac{e^{- 4 x}}{2} \\
\left.\frac{\partial^1 (G_+(x))}{\partial x^1}\right|_{x=1} &= - \frac{1}{2 e^{4}} \\
\frac{\partial^2 (G_+(x))}{\partial x^2} &= 2 e^{- 4 x} \\
\left.\frac{\partial^2 (G_+(x))}{\partial x^2}\right|_{x=1} &= \frac{2}{e^{4}}
\end{aligned}
$$
:::

### 外力を畳み込みで求める

$L=-\partial_x^2+\mu^2$（$\mu>0$）に対し、減衰条件を選んだ基本解を $G=e^{-\mu|x|}/(2\mu)$ とする。コンパクトな台を持つ滑らかな $f$ については $L(G*f)=(LG)*f=f$ となる。これは試験関数に対する二階微分を積分の内側で $G$ へ移し、デルタの作用で $f$ を選ぶためである。入力も点の有限和なら $G$ の平行移動の和になる。例えば $f=\delta_a-\delta_{-a}$ では、遠くで指数減衰し、かつ奇対称になる。原点で解が連続で零であっても、源のある二点では導関数が跳ぶ。源の符号が異なるため、その跳びの符号も逆である。非減衰な同次解を足すと同じ分布方程式を満たしても無限遠条件を失う。
:::example[二つの点源の間での核の重ね合わせ]
減衰率 $2$、点源を $+1,-1$ に置く。区間 $-1<x<1$ に限り、符号を決めた二指数の差を調べる。係数を読み、試行関数を微分して方程式へ戻す。

この例の係数を左から二階、一階、零階の順に書く。右辺がある場合は三つの項の和と比較する。

$$
\left( -1, \  0, \  4\right)
$$

関数とその二つの導関数を並べる。特に一つ下の階数の導関数には、一階項の係数を掛ける。

$$
\begin{aligned}u&=- \frac{e^{- 2 x - 2}}{4} + \frac{e^{2 x - 2}}{4}\\u^{\prime}&=\frac{e^{- 2 x - 2}}{2} + \frac{e^{2 x - 2}}{2}\\u^{\prime\prime}&=- e^{- 2 x - 2} + e^{2 \left(x - 1\right)}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=e^{- 2 x - 2} - e^{2 x - 2}+\left(0\right)+\left(- e^{- 2 x - 2} + e^{2 x - 2}\right)\\&=0\end{aligned}
$$

$x$ が $0$ のところで値を評価する。変換した方程式を使う場合も、初期値や境界値は元の変数で照合する。

$$
\begin{aligned}u&=0\\ u^\prime&=e^{-2}\\ Lu&=0\end{aligned}
$$
:::

## 演習問題

::::exercise{level=basic}
箱 $\mathbf1_{[-2,2]}$ の微分を求めよ。
:::solution
左端の跳びは $1$、右端は $-1$ であり、式 :eqref[eq-mp-distribution-jump] で $\delta_{-2}-\delta_2$ となる。
:::
::::

::::exercise{level=basic}
$x^2\delta'_3$ を簡約せよ。
:::solution
試験関数への作用は $-(x^2\varphi)'(3)=-9\varphi'(3)-6\varphi(3)$ なので $9\delta'_3-6\delta_3$ である。
:::
::::

::::exercise{level=basic}
$\delta(3x-6)$ を簡約せよ。
:::solution
零点は $2$、導関数は $3$ なので $\delta_2/3$ である。
:::
::::

::::exercise{level=basic}
$\widehat{\delta_{-1}-\delta_1}$ を求めよ。
:::solution
平行移動の指数を差し引くと $e^{ik}-e^{-ik}=2i\sin k$ になる。
:::
::::

::::exercise{level=advanced}
$\delta(x^3-x)$ の三つの重みを求めよ。
:::solution
零点は $-1,0,1$、導関数 $3x^2-1$ の絶対値は順に $2,1,2$ である。従って $(\delta_{-1}+\delta_1)/2+\delta_0$ になる。
:::
::::

::::exercise{level=advanced}
$x^2\operatorname{Fp}(1/x^2)=1$ を示せ。
:::solution
有限部分に $x^2\varphi$ を入れると原点値が零なので、差し引きは不要になる。残りは $\int_{|x|>\epsilon}\varphi$ で、極限は $\int\varphi$ である。
:::
::::

::::exercise{level=advanced}
$G=e^{-2|x|}/4$ が $(-\partial_x^2+4)G=\delta_0$ を満たすことを示せ。
:::solution
原点外では $G''=4G$ である。$G$ は連続だが、導関数は左で $1/2$、右で $-1/2$ なので跳びは $-1$ である。従って $G''=4G-\delta_0$ となり結論を得る。
:::
::::

::::exercise{level=advanced}
$\varphi=e^{-x^2}$ に対する $\operatorname{Fp}(1/x^2)$ の値を求めよ。
:::solution
微分との関係から $\langle\operatorname{PV}(1/x),\varphi'\rangle=-2\int_\R e^{-x^2}\dd x=-2\sqrt\pi$ である。Gaussian は急減少なので無限遠の境界項は消える。
:::
::::

## 確認テスト

::::quiz{answer="2" tol="1e-9"}
$|x|''=a\delta_0$ の $a$ はいくらか。

:::explanation
一階導関数は左右で $-1,+1$ なので、跳びは $1-(-1)=2$ である。
:::
::::

::::quiz
$\widehat1$ はどれか。

- [ ] $1$
- [x] $2\pi\delta_0$
- [ ] $\delta_0/(2\pi)$

:::explanation
逆変換を原点で評価すると $\int\widehat\varphi(x)\dd x=2\pi\varphi(0)$ である。
:::
::::

::::quiz{answer="0.5" tol="1e-9"}
$\delta(2x-4)=a\delta_2$ の重み $a$ を求めよ。

:::explanation
絶対 Jacobian は $2$ なので $a=1/2$ である。
:::
::::

::::quiz
$\delta_a$ の微分の試験関数への作用はどれか。

- [x] $-\varphi'(a)$
- [ ] $\varphi'(a)$
- [ ] $\varphi(a)$

:::explanation
部分積分に合わせた定義から $\langle\delta'_a,\varphi\rangle=-\langle\delta_a,\varphi'\rangle=-\varphi'(a)$ である。
:::
::::

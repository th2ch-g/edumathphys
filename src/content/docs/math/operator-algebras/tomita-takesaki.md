---
title: Tomita–Takesaki 理論
description: 巡回分離ベクトルからモジュラー作用素を作り、自己同型群、KMS 条件、型の不変量へ進む。
notation:
  - symbol: 'S'
    meaning: Tomita 作用素
  - symbol: '\Delta'
    meaning: モジュラー作用素
  - symbol: 'J'
    meaning: モジュラー共役
  - symbol: '\sigma_t^\varphi'
    meaning: 状態 φ のモジュラー自己同型
---

## この章の目標

- 巡回分離ベクトルから閉じた反線形作用素を作る。
- 有限次元で $\Delta$ と $J$ を具体的に計算する。
- モジュラー群と KMS 条件の符号を合わせる。
- 状態の変更と III 型の不変量の関係を理解する。

## 巡回分離ベクトルと標準形

複素 Hilbert 空間 $H$ 上の von Neumann 環 $M$ に巡回分離単位ベクトル $\xi$ があるとする。各 $x$ の作用 $x\xi$ を知れば分離性で $x$ が決まり、巡回性でこの種のベクトルは稠密である。従って代数の対合を稠密な線形領域上の写像として運べる。内積は第 1 引数について線形である。

忠実正規状態 $\varphi$ から出発しても、その GNS 表現を使えばこの状況を作れる。$\norm{x\xi}^2=\varphi(x^*x)$ だから忠実性により $\xi$ は分離的であり、表現も忠実である。巡回性は GNS 構成による。巡回軌道上の行列係数 $\varphi(b^*xa)$ は $x$ に関して σ 弱連続であり、表現の単位球では一様有界性により任意の行列係数も σ 弱連続である。従って表現の単位球は WOT コンパクトであり、:ref[thm-oa-kaplansky] により像は von Neumann 環である。

:::definition[標準形]{#def-oa-standard-form}
von Neumann 環の :term[標準形]{en="standard form" yomi="ひょうじゅんけい"} $(M,H,J,\mathcal P)$ では、$J$ は反ユニタリ対合、$\mathcal P$ は自己双対な閉錐で、$JMJ=M'$、$JzJ=z^*$（$z\in Z(M)$）、$J$ が $\mathcal P$ を点ごとに固定し、$xJxJ\mathcal P\subset\mathcal P$（$x\in M$）が成り立つ。自己双対とは $\mathcal P=\{\eta\in H:\inner{\eta}{\zeta}\in[0,\infty)\text{ for all }\zeta\in\mathcal P\}$ が成り立つことである。
:::

:::example[行列環の標準形]
$H=M_n$ に $\inner a b=\tr(b^*a)$ を入れ、$M$ を左乗法、$J(a)=a^*$、$\mathcal P=M_n^+$ とする。正錐の自己双対性は階数一射影で試せば従う。$J$ は錐を固定し、$JL_xJ=R_{x^*}$、$L_xJL_xJ(a)=xax^*$ である。左作用の可換子の $T$ には $T(a)=aT(1)$ が成り立つので右作用だけである。中心の左スカラーは $J$ によって共役スカラーに移る。
:::

一般の巡回分離表現では、後で得る $\Delta$ を用いて :term[自然正錐]{en="natural positive cone" yomi="しぜんせいすい"} $\mathcal P=\overline{\{\Delta^{1/4}x\xi:x\in M_+\}}$ を作ると標準形が得られる。各 $x\xi$ は $\Delta^{1/2}$ の定義域に入るので $\Delta^{1/4}$ を作用させられる。全ての von Neumann 環に標準形が存在し、標準形では任意の正規正値汎関数が $\mathcal P$ に属する一意なベクトルのベクトル汎関数として表される。また正規 *-同型は錐を保つ一意なユニタリで実装される。これらの一般形には自然正錐の自己双対性と、忠実正規状態を持たない環の構成が必要なので、:cite[takesaki-structure] の標準双加群を参照する。

### 標準形のベクトルと状態

行列環の標準形では、正錐の単位ベクトル $D^{1/2}$ が状態 $\varphi(x)=\tr(Dx)$ を与える。逆に正錐のベクトル $a$ が同じ状態なら $\tr(a^2x)=\tr(Dx)$ が全行列で成り立つため、$a^2=D$ となる。正平方根の一意性から $a=D^{1/2}$ である。正錐を指定することで、同じ状態を表す複数のベクトルから一つを選べる。

:::example[状態に忠実性がないと分離性を失う]
$H=M_2$、$\xi=E_{11}$ とすると $M\xi$ は第 2 列が零の行列だけからなる。この空間は二次元で、$H$ 全体は四次元なので巡回的でない。また $E_{22}\xi=0$ なので分離的でもない。これに対し $\eta=(E_{11}+E_{22})/\sqrt2$ なら任意の $a\in H$ は $\sqrt2 a\eta$ と書ける。さらに $a\eta=0$ なら $a=0$ なので、巡回分離的である。
:::

:::example[二つの密度に異なる重みを置く]
$M_d$ を Hilbert–Schmidt 内積の行列空間に左から作用させる。状態の密度 $D$ に対する巡回分離ベクトルと、行列 $a$ をそのベクトルに作用させる方法を計算する。

密度は正の対角行列でトレースは $1$ である。標準形の正錐に入る $\xi=D^{1/2}$ を取ると、$\norm{\xi}^2=\tr D=1$ である。

$$
\begin{aligned}
D&=\begin{pmatrix}\frac{1}{4} & 0\\0 & \frac{3}{4}\end{pmatrix}\\
\xi&=\begin{pmatrix}\frac{1}{2} & 0\\0 & \frac{\sqrt{3}}{2}\end{pmatrix}\\
a&=\begin{pmatrix}1 & 1\\1 & 1\end{pmatrix}
\end{aligned}
$$

$a\xi$ では $a$ の両側からではなく、右から平方根を掛けるため、列がそれぞれ異なる倍率を受ける。状態の値 $\varphi(a)$ は $\inner{a\xi}{\xi}=\tr(Da)$ になる。

$$
\begin{aligned}
a\xi&=\begin{pmatrix}\frac{1}{2} & \frac{\sqrt{3}}{2}\\\frac{1}{2} & \frac{\sqrt{3}}{2}\end{pmatrix}\\
\varphi(a)&=1\\
\norm{a\xi}^2&=2
\end{aligned}
$$

$a$ は $a\xi$ に $\xi^{-1}$ を右から掛けて復元できる。この逆行列があるため $a\xi=0$ なら $a=0$ であり、任意の Hilbert–Schmidt ベクトルも $a\xi$ と書ける。

$$
\begin{aligned}
\xi^{-1}&=\begin{pmatrix}2 & 0\\0 & \frac{2 \sqrt{3}}{3}\end{pmatrix}\\
(a\xi)\xi^{-1}&=\begin{pmatrix}1 & 1\\1 & 1\end{pmatrix}
\end{aligned}
$$

一般の環ではこの復元写像が有界になるとは限らない。有限模型での逆密度に相当する部分がモジュラー作用素の非有界性へつながる。
:::

:::example[逆密度で第二区画を伸ばす]
$M_d$ を Hilbert–Schmidt 内積の行列空間に左から作用させる。状態の密度 $D$ に対する巡回分離ベクトルと、行列 $a$ をそのベクトルに作用させる方法を計算する。

密度は正の対角行列でトレースは $1$ である。標準形の正錐に入る $\xi=D^{1/2}$ を取ると、$\norm{\xi}^2=\tr D=1$ である。

$$
\begin{aligned}
D&=\begin{pmatrix}\frac{4}{5} & 0\\0 & \frac{1}{5}\end{pmatrix}\\
\xi&=\begin{pmatrix}\frac{2 \sqrt{5}}{5} & 0\\0 & \frac{\sqrt{5}}{5}\end{pmatrix}\\
a&=\begin{pmatrix}1 & 2\\2 & 4\end{pmatrix}
\end{aligned}
$$

$a\xi$ では $a$ の両側からではなく、右から平方根を掛けるため、列がそれぞれ異なる倍率を受ける。状態の値 $\varphi(a)$ は $\inner{a\xi}{\xi}=\tr(Da)$ になる。

$$
\begin{aligned}
a\xi&=\begin{pmatrix}\frac{2 \sqrt{5}}{5} & \frac{2 \sqrt{5}}{5}\\\frac{4 \sqrt{5}}{5} & \frac{4 \sqrt{5}}{5}\end{pmatrix}\\
\varphi(a)&=\frac{8}{5}\\
\norm{a\xi}^2&=8
\end{aligned}
$$

$a$ は $a\xi$ に $\xi^{-1}$ を右から掛けて復元できる。この逆行列があるため $a\xi=0$ なら $a=0$ であり、任意の Hilbert–Schmidt ベクトルも $a\xi$ と書ける。

$$
\begin{aligned}
\xi^{-1}&=\begin{pmatrix}\frac{\sqrt{5}}{2} & 0\\0 & \sqrt{5}\end{pmatrix}\\
(a\xi)\xi^{-1}&=\begin{pmatrix}1 & 2\\2 & 4\end{pmatrix}
\end{aligned}
$$

一般の環ではこの復元写像が有界になるとは限らない。有限模型での逆密度に相当する部分がモジュラー作用素の非有界性へつながる。
:::

:::example[三つのスペクトル成分を持つ状態]
$M_d$ を Hilbert–Schmidt 内積の行列空間に左から作用させる。状態の密度 $D$ に対する巡回分離ベクトルと、行列 $a$ をそのベクトルに作用させる方法を計算する。

密度は正の対角行列でトレースは $1$ である。標準形の正錐に入る $\xi=D^{1/2}$ を取ると、$\norm{\xi}^2=\tr D=1$ である。

$$
\begin{aligned}
D&=\begin{pmatrix}\frac{2}{5} & 0 & 0\\0 & \frac{2}{5} & 0\\0 & 0 & \frac{1}{5}\end{pmatrix}\\
\xi&=\begin{pmatrix}\frac{\sqrt{10}}{5} & 0 & 0\\0 & \frac{\sqrt{10}}{5} & 0\\0 & 0 & \frac{\sqrt{5}}{5}\end{pmatrix}\\
a&=\begin{pmatrix}1 & 2 & 1\\2 & 4 & 2\\1 & 2 & 1\end{pmatrix}
\end{aligned}
$$

$a\xi$ では $a$ の両側からではなく、右から平方根を掛けるため、列がそれぞれ異なる倍率を受ける。状態の値 $\varphi(a)$ は $\inner{a\xi}{\xi}=\tr(Da)$ になる。

$$
\begin{aligned}
a\xi&=\begin{pmatrix}\frac{\sqrt{10}}{5} & \frac{2 \sqrt{10}}{5} & \frac{\sqrt{5}}{5}\\\frac{2 \sqrt{10}}{5} & \frac{4 \sqrt{10}}{5} & \frac{2 \sqrt{5}}{5}\\\frac{\sqrt{10}}{5} & \frac{2 \sqrt{10}}{5} & \frac{\sqrt{5}}{5}\end{pmatrix}\\
\varphi(a)&=\frac{11}{5}\\
\norm{a\xi}^2&=\frac{66}{5}
\end{aligned}
$$

$a$ は $a\xi$ に $\xi^{-1}$ を右から掛けて復元できる。この逆行列があるため $a\xi=0$ なら $a=0$ であり、任意の Hilbert–Schmidt ベクトルも $a\xi$ と書ける。

$$
\begin{aligned}
\xi^{-1}&=\begin{pmatrix}\frac{\sqrt{10}}{2} & 0 & 0\\0 & \frac{\sqrt{10}}{2} & 0\\0 & 0 & \sqrt{5}\end{pmatrix}\\
(a\xi)\xi^{-1}&=\begin{pmatrix}1 & 2 & 1\\2 & 4 & 2\\1 & 2 & 1\end{pmatrix}
\end{aligned}
$$

一般の環ではこの復元写像が有界になるとは限らない。有限模型での逆密度に相当する部分がモジュラー作用素の非有界性へつながる。
:::

:::example[規格化トレースでの標準形]
$M_d$ を Hilbert–Schmidt 内積の行列空間に左から作用させる。状態の密度 $D$ に対する巡回分離ベクトルと、行列 $a$ をそのベクトルに作用させる方法を計算する。

密度は正の対角行列でトレースは $1$ である。標準形の正錐に入る $\xi=D^{1/2}$ を取ると、$\norm{\xi}^2=\tr D=1$ である。

$$
\begin{aligned}
D&=\begin{pmatrix}\frac{1}{3} & 0 & 0\\0 & \frac{1}{3} & 0\\0 & 0 & \frac{1}{3}\end{pmatrix}\\
\xi&=\begin{pmatrix}\frac{\sqrt{3}}{3} & 0 & 0\\0 & \frac{\sqrt{3}}{3} & 0\\0 & 0 & \frac{\sqrt{3}}{3}\end{pmatrix}\\
a&=\begin{pmatrix}1 & 1 & 1\\1 & 1 & 1\\1 & 1 & 1\end{pmatrix}
\end{aligned}
$$

$a\xi$ では $a$ の両側からではなく、右から平方根を掛けるため、列がそれぞれ異なる倍率を受ける。状態の値 $\varphi(a)$ は $\inner{a\xi}{\xi}=\tr(Da)$ になる。

$$
\begin{aligned}
a\xi&=\begin{pmatrix}\frac{\sqrt{3}}{3} & \frac{\sqrt{3}}{3} & \frac{\sqrt{3}}{3}\\\frac{\sqrt{3}}{3} & \frac{\sqrt{3}}{3} & \frac{\sqrt{3}}{3}\\\frac{\sqrt{3}}{3} & \frac{\sqrt{3}}{3} & \frac{\sqrt{3}}{3}\end{pmatrix}\\
\varphi(a)&=1\\
\norm{a\xi}^2&=3
\end{aligned}
$$

$a$ は $a\xi$ に $\xi^{-1}$ を右から掛けて復元できる。この逆行列があるため $a\xi=0$ なら $a=0$ であり、任意の Hilbert–Schmidt ベクトルも $a\xi$ と書ける。

$$
\begin{aligned}
\xi^{-1}&=\begin{pmatrix}\sqrt{3} & 0 & 0\\0 & \sqrt{3} & 0\\0 & 0 & \sqrt{3}\end{pmatrix}\\
(a\xi)\xi^{-1}&=\begin{pmatrix}1 & 1 & 1\\1 & 1 & 1\\1 & 1 & 1\end{pmatrix}
\end{aligned}
$$

一般の環ではこの復元写像が有界になるとは限らない。有限模型での逆密度に相当する部分がモジュラー作用素の非有界性へつながる。
:::

:::example[右作用と左作用の区別]
$M_d$ を Hilbert–Schmidt 内積の行列空間に左から作用させる。状態の密度 $D$ に対する巡回分離ベクトルと、行列 $a$ をそのベクトルに作用させる方法を計算する。

密度は正の対角行列でトレースは $1$ である。標準形の正錐に入る $\xi=D^{1/2}$ を取ると、$\norm{\xi}^2=\tr D=1$ である。

$$
\begin{aligned}
D&=\begin{pmatrix}\frac{1}{7} & 0 & 0\\0 & \frac{2}{7} & 0\\0 & 0 & \frac{4}{7}\end{pmatrix}\\
\xi&=\begin{pmatrix}\frac{\sqrt{7}}{7} & 0 & 0\\0 & \frac{\sqrt{14}}{7} & 0\\0 & 0 & \frac{2 \sqrt{7}}{7}\end{pmatrix}\\
a&=\begin{pmatrix}1 & 2 & 0\\2 & 4 & 0\\0 & 0 & 0\end{pmatrix}
\end{aligned}
$$

$a\xi$ では $a$ の両側からではなく、右から平方根を掛けるため、列がそれぞれ異なる倍率を受ける。状態の値 $\varphi(a)$ は $\inner{a\xi}{\xi}=\tr(Da)$ になる。

$$
\begin{aligned}
a\xi&=\begin{pmatrix}\frac{\sqrt{7}}{7} & \frac{2 \sqrt{14}}{7} & 0\\\frac{2 \sqrt{7}}{7} & \frac{4 \sqrt{14}}{7} & 0\\0 & 0 & 0\end{pmatrix}\\
\varphi(a)&=\frac{9}{7}\\
\norm{a\xi}^2&=\frac{45}{7}
\end{aligned}
$$

$a$ は $a\xi$ に $\xi^{-1}$ を右から掛けて復元できる。この逆行列があるため $a\xi=0$ なら $a=0$ であり、任意の Hilbert–Schmidt ベクトルも $a\xi$ と書ける。

$$
\begin{aligned}
\xi^{-1}&=\begin{pmatrix}\sqrt{7} & 0 & 0\\0 & \frac{\sqrt{14}}{2} & 0\\0 & 0 & \frac{\sqrt{7}}{2}\end{pmatrix}\\
(a\xi)\xi^{-1}&=\begin{pmatrix}1 & 2 & 0\\2 & 4 & 0\\0 & 0 & 0\end{pmatrix}
\end{aligned}
$$

一般の環ではこの復元写像が有界になるとは限らない。有限模型での逆密度に相当する部分がモジュラー作用素の非有界性へつながる。
:::

## モジュラー作用素とモジュラー共役

:::proposition[対合の閉包と極分解]{#prop-oa-tomita-operator}
$S_0(x\xi)=x^*\xi$ は閉包可能である。その閉包は $S=J\Delta^{1/2}$ と極分解でき、$J^2=1$、$J\Delta J=\Delta^{-1}$ である。ここで $\Delta=S^*S$ は正自己共役作用素で、核が零である。第 1 引数を線形とする内積では、反線形作用素の随伴は $\inner{Su}{v}=\inner{S^*v}{u}$ で定める。$J$ と $\Delta$ をそれぞれ :term[モジュラー共役]{en="modular conjugation" yomi="もじゅらーきょうやく"} と :term[モジュラー作用素]{en="modular operator" yomi="もじゅらーさようそ"} という。
:::
:::proof
:ref[prop-oa-cyclic-separating] により $M'\xi$ は稠密である。$y\in M'$ に対し $\inner{S_0x\xi}{y\xi}=\inner{y^*\xi}{x\xi}$ なので、$x_n\xi\to0$、$S_0x_n\xi\to\eta$ なら $\inner{\eta}{y\xi}=0$ が全 $y$ で成立して $\eta=0$ である。これが閉包可能性である。グラフは座標交換で不変な $S_0$ のグラフの閉包であり、従って $S^{-1}=S$ である。特に $S$ は単射で値域は稠密である。

反線形の場合の極分解は、任意の反ユニタリ対合 $C_0$ で $T_0=C_0S$ を線形にしてから得られる。ここで閉線形作用素の正部分が存在する理由を確認する。$T_0$ の閉グラフへの $(\eta,0)$ の直交射影を $(x,T_0x)$ とする。グラフとの直交性から $T_0x$ は $T_0^*$ の定義域に入り、$(1+T_0^*T_0)x=\eta$ となる。従って $1+T_0^*T_0$ は全射である。$T_0^*T_0$ は正で、定義域に直交する $\eta$ にこの全射性を使うと $0=\inner{\eta}{x}=\norm x^2+\norm{T_0x}^2$ から $\eta=0$ となる。よって稠密定義である。$\Delta=T_0^*T_0$ の対称性と $\Img(1+\Delta)=H$ により自己共役性が出る。実際、$z\in\operatorname{Dom}\Delta^*$ に対し $(1+\Delta)y=(1+\Delta^*)z$ を解けば $z-y$ は $\Ker(1+\Delta^*)$ に入り、この核は全射性により零である。

$\operatorname{Dom}\Delta$ 上では $\norm{T_0x}^2=\inner{\Delta x}{x}=\norm{\Delta^{1/2}x}^2$ である。この領域は $T_0$ のグラフノルムで稠密である。実際、グラフ内でその部分に直交する $(v,T_0v)$ に対して、任意の $x\in\operatorname{Dom}\Delta$ で $\inner v{(1+\Delta)x}=0$ となるため、全射性から $v=0$ である。同じ領域はスペクトルを有界区間へ切る近似により、$\Delta^{1/2}$ のグラフノルムでも稠密である。よって $T_0$ と $\Delta^{1/2}$ の定義域とグラフノルムが一致する。$\Delta^{1/2}x\mapsto T_0x$ は等長に延び、両方の値域が稠密なのでユニタリとなる。左に $C_0$ を掛けた反ユニタリを $J$ とすれば $S=J\Delta^{1/2}$ である。この等長な対応は稠密集合上で決まるので一意である。

$S^{-1}=\Delta^{-1/2}J^{-1}=J^{-1}(J\Delta^{-1/2}J^{-1})$ を極分解の形で比較すると $J=J^{-1}$ と $\Delta^{1/2}=J\Delta^{-1/2}J^{-1}$ が従う。ここで $T_0=C_0S$ のため $T_0^*T_0=S^*S$ である。
:::

:::example[忠実な密度行列]
$D=\diag(d_1,\dots,d_n)>0$、$\tr D=1$ とし、前節の模型で $\xi=D^{1/2}$ とする。可逆性で巡回分離である。$S(aD^{1/2})=a^*D^{1/2}$ から $S(b)=D^{-1/2}b^*D^{1/2}$ となる。$\Delta(E_{ij})=(d_i/d_j)E_{ij}$、$J(b)=b^*$ とすると $J\Delta^{1/2}=S$ である。$D=\diag(1/4,3/4)$ では固有値は $1,1,1/3,3$ となる。
:::

### 二つの密度を同じ Hilbert 空間で比較する

有限次元の場合、$\xi=D^{1/2}$ とし $a=xD^{1/2}$ から $x=aD^{-1/2}$ を復元できる。従って
$$
S_D(a)=D^{-1/2}a^*D^{1/2},\quad
\Delta_D(a)=DaD^{-1},\quad J(a)=a^*
$$
である。最初の式は $x^*D^{1/2}$ の展開であり、第 2 式の正平方根が $a\mapsto D^{1/2}aD^{-1/2}$ であることは $D$ の固有基底で確かめられる。モジュラー作用素は密度と共に変化するが、$J$ は標準形の上では同じである。

:::example[対数から加法的な周波数を得る]
$D=\diag(1/7,2/7,4/7)$ とする。$\log\Delta_D(E_{ij})=(\log d_i-\log d_j)E_{ij}$ であり、$E_{13}$ の固有値は $-\log4$、$E_{31}$ の固有値は $\log4$ である。$J$ で二つの行列単位が入れ替わるので、$\log\Delta_D$ の符号が反転する。対角成分は $\log\Delta_D$ で零になり、密度の異なる固有空間を結ぶ方向にだけ非零の周波数が現れる。
:::

:::example[最小固有値を結ぶ方向]
$H=M_d$ の正規直交基底 $(E_{ij})$ に対する Tomita 作用素 $S$ と正部分 $\Delta^{1/2}$ を求める。密度は対角なので各行列単位の作用が独立に計算できる。

$a=E_{ij}$ に対し $D^{-1/2}a^*D^{1/2}=\sqrt{d_i/d_j}E_{ji}$ となる。従って $S$ は随伴で行と列を交換し、交換する前の添字の比でスケールする。

$$
\begin{aligned}
D&=\begin{pmatrix}\frac{1}{5} & 0\\0 & \frac{4}{5}\end{pmatrix}\\
(d_i/d_j)_{ij}&=\begin{pmatrix}1 & \frac{1}{4}\\4 & 1\end{pmatrix}
\end{aligned}
$$

$\Delta(E_{ij})=(d_i/d_j)E_{ij}$ なので、この表の全成分を集めたものがスペクトルである。指定した一つの方向を選び、$S$ と $\Delta^{1/2}$ の違いを確認する。

$$
\begin{aligned}
(i,j)&=\left( 1, \  2\right)\\
\Delta^{1/2} E_{ij}&=\begin{pmatrix}0 & \frac{1}{2}\\0 & 0\end{pmatrix}\\
S(E_{ij})&=\begin{pmatrix}0 & 0\\\frac{1}{2} & 0\end{pmatrix}
\end{aligned}
$$

$J$ は随伴なので $J\Delta^{1/2}$ はこの $S$ に一致する。$S$ を二度使うと二つの倍率の積は $\sqrt{d_i/d_j}\sqrt{d_j/d_i}=1$ になる。

$$
\begin{aligned}
\sigma(\Delta)&=\left\{\frac{1}{4}, 1, 4\right\}\\
\norm{\Delta}&=4
\end{aligned}
$$

対角の行列単位はすべて固定される。密度が同じ値の座標が二つあれば、その間の非対角の行列単位も固有値一を持つ。
:::

:::example[モジュラー倍率が一より大きい方向]
$H=M_d$ の正規直交基底 $(E_{ij})$ に対する Tomita 作用素 $S$ と正部分 $\Delta^{1/2}$ を求める。密度は対角なので各行列単位の作用が独立に計算できる。

$a=E_{ij}$ に対し $D^{-1/2}a^*D^{1/2}=\sqrt{d_i/d_j}E_{ji}$ となる。従って $S$ は随伴で行と列を交換し、交換する前の添字の比でスケールする。

$$
\begin{aligned}
D&=\begin{pmatrix}\frac{3}{4} & 0\\0 & \frac{1}{4}\end{pmatrix}\\
(d_i/d_j)_{ij}&=\begin{pmatrix}1 & 3\\\frac{1}{3} & 1\end{pmatrix}
\end{aligned}
$$

$\Delta(E_{ij})=(d_i/d_j)E_{ij}$ なので、この表の全成分を集めたものがスペクトルである。指定した一つの方向を選び、$S$ と $\Delta^{1/2}$ の違いを確認する。

$$
\begin{aligned}
(i,j)&=\left( 1, \  2\right)\\
\Delta^{1/2} E_{ij}&=\begin{pmatrix}0 & \sqrt{3}\\0 & 0\end{pmatrix}\\
S(E_{ij})&=\begin{pmatrix}0 & 0\\\sqrt{3} & 0\end{pmatrix}
\end{aligned}
$$

$J$ は随伴なので $J\Delta^{1/2}$ はこの $S$ に一致する。$S$ を二度使うと二つの倍率の積は $\sqrt{d_i/d_j}\sqrt{d_j/d_i}=1$ になる。

$$
\begin{aligned}
\sigma(\Delta)&=\left\{\frac{1}{3}, 1, 3\right\}\\
\norm{\Delta}&=3
\end{aligned}
$$

対角の行列単位はすべて固定される。密度が同じ値の座標が二つあれば、その間の非対角の行列単位も固有値一を持つ。
:::

:::example[三つの倍数からなるスペクトル]
$H=M_d$ の正規直交基底 $(E_{ij})$ に対する Tomita 作用素 $S$ と正部分 $\Delta^{1/2}$ を求める。密度は対角なので各行列単位の作用が独立に計算できる。

$a=E_{ij}$ に対し $D^{-1/2}a^*D^{1/2}=\sqrt{d_i/d_j}E_{ji}$ となる。従って $S$ は随伴で行と列を交換し、交換する前の添字の比でスケールする。

$$
\begin{aligned}
D&=\begin{pmatrix}\frac{1}{7} & 0 & 0\\0 & \frac{2}{7} & 0\\0 & 0 & \frac{4}{7}\end{pmatrix}\\
(d_i/d_j)_{ij}&=\begin{pmatrix}1 & \frac{1}{2} & \frac{1}{4}\\2 & 1 & \frac{1}{2}\\4 & 2 & 1\end{pmatrix}
\end{aligned}
$$

$\Delta(E_{ij})=(d_i/d_j)E_{ij}$ なので、この表の全成分を集めたものがスペクトルである。指定した一つの方向を選び、$S$ と $\Delta^{1/2}$ の違いを確認する。

$$
\begin{aligned}
(i,j)&=\left( 2, \  3\right)\\
\Delta^{1/2} E_{ij}&=\begin{pmatrix}0 & 0 & 0\\0 & 0 & \frac{\sqrt{2}}{2}\\0 & 0 & 0\end{pmatrix}\\
S(E_{ij})&=\begin{pmatrix}0 & 0 & 0\\0 & 0 & 0\\0 & \frac{\sqrt{2}}{2} & 0\end{pmatrix}
\end{aligned}
$$

$J$ は随伴なので $J\Delta^{1/2}$ はこの $S$ に一致する。$S$ を二度使うと二つの倍率の積は $\sqrt{d_i/d_j}\sqrt{d_j/d_i}=1$ になる。

$$
\begin{aligned}
\sigma(\Delta)&=\left\{\frac{1}{4}, \frac{1}{2}, 1, 2, 4\right\}\\
\norm{\Delta}&=4
\end{aligned}
$$

対角の行列単位はすべて固定される。密度が同じ値の座標が二つあれば、その間の非対角の行列単位も固有値一を持つ。
:::

:::example[重複した密度の間では倍率が一]
$H=M_d$ の正規直交基底 $(E_{ij})$ に対する Tomita 作用素 $S$ と正部分 $\Delta^{1/2}$ を求める。密度は対角なので各行列単位の作用が独立に計算できる。

$a=E_{ij}$ に対し $D^{-1/2}a^*D^{1/2}=\sqrt{d_i/d_j}E_{ji}$ となる。従って $S$ は随伴で行と列を交換し、交換する前の添字の比でスケールする。

$$
\begin{aligned}
D&=\begin{pmatrix}\frac{1}{5} & 0 & 0\\0 & \frac{1}{5} & 0\\0 & 0 & \frac{3}{5}\end{pmatrix}\\
(d_i/d_j)_{ij}&=\begin{pmatrix}1 & 1 & \frac{1}{3}\\1 & 1 & \frac{1}{3}\\3 & 3 & 1\end{pmatrix}
\end{aligned}
$$

$\Delta(E_{ij})=(d_i/d_j)E_{ij}$ なので、この表の全成分を集めたものがスペクトルである。指定した一つの方向を選び、$S$ と $\Delta^{1/2}$ の違いを確認する。

$$
\begin{aligned}
(i,j)&=\left( 1, \  2\right)\\
\Delta^{1/2} E_{ij}&=\begin{pmatrix}0 & 1 & 0\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix}\\
S(E_{ij})&=\begin{pmatrix}0 & 0 & 0\\1 & 0 & 0\\0 & 0 & 0\end{pmatrix}
\end{aligned}
$$

$J$ は随伴なので $J\Delta^{1/2}$ はこの $S$ に一致する。$S$ を二度使うと二つの倍率の積は $\sqrt{d_i/d_j}\sqrt{d_j/d_i}=1$ になる。

$$
\begin{aligned}
\sigma(\Delta)&=\left\{\frac{1}{3}, 1, 3\right\}\\
\norm{\Delta}&=3
\end{aligned}
$$

対角の行列単位はすべて固定される。密度が同じ値の座標が二つあれば、その間の非対角の行列単位も固有値一を持つ。
:::

:::example[密度の端から端に移る作用]
$H=M_d$ の正規直交基底 $(E_{ij})$ に対する Tomita 作用素 $S$ と正部分 $\Delta^{1/2}$ を求める。密度は対角なので各行列単位の作用が独立に計算できる。

$a=E_{ij}$ に対し $D^{-1/2}a^*D^{1/2}=\sqrt{d_i/d_j}E_{ji}$ となる。従って $S$ は随伴で行と列を交換し、交換する前の添字の比でスケールする。

$$
\begin{aligned}
D&=\begin{pmatrix}\frac{2}{11} & 0 & 0\\0 & \frac{3}{11} & 0\\0 & 0 & \frac{6}{11}\end{pmatrix}\\
(d_i/d_j)_{ij}&=\begin{pmatrix}1 & \frac{2}{3} & \frac{1}{3}\\\frac{3}{2} & 1 & \frac{1}{2}\\3 & 2 & 1\end{pmatrix}
\end{aligned}
$$

$\Delta(E_{ij})=(d_i/d_j)E_{ij}$ なので、この表の全成分を集めたものがスペクトルである。指定した一つの方向を選び、$S$ と $\Delta^{1/2}$ の違いを確認する。

$$
\begin{aligned}
(i,j)&=\left( 1, \  3\right)\\
\Delta^{1/2} E_{ij}&=\begin{pmatrix}0 & 0 & \frac{\sqrt{3}}{3}\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix}\\
S(E_{ij})&=\begin{pmatrix}0 & 0 & 0\\0 & 0 & 0\\\frac{\sqrt{3}}{3} & 0 & 0\end{pmatrix}
\end{aligned}
$$

$J$ は随伴なので $J\Delta^{1/2}$ はこの $S$ に一致する。$S$ を二度使うと二つの倍率の積は $\sqrt{d_i/d_j}\sqrt{d_j/d_i}=1$ になる。

$$
\begin{aligned}
\sigma(\Delta)&=\left\{\frac{1}{3}, \frac{1}{2}, \frac{2}{3}, 1, \frac{3}{2}, 2, 3\right\}\\
\norm{\Delta}&=3
\end{aligned}
$$

対角の行列単位はすべて固定される。密度が同じ値の座標が二つあれば、その間の非対角の行列単位も固有値一を持つ。
:::

:::example[非有界なモジュラー作用素]
$K=\ell^2(\N)$ とし、基底を $j\ge1$ で番号付けする。$D=\diag(2^{-j})$ は核が零でトレースが $1$ である。$B(K)$ を Hilbert–Schmidt 作用素の空間に左乗法で作用させ、$\xi=D^{1/2}$ とする。$xD^{1/2}=0$ なら $D^{1/2}$ の稠密な値域を $x$ が消すため $x=0$ である。また各行列単位 $E_{ij}$ は $E_{ij}D^{1/2}$ の定数倍なので $\xi$ は巡回的である。

$b=x\xi$ の行列成分を $b_{ij}=2^{-j/2}x_{ij}$ と書くと、$S_0 b$ の $(j,i)$ 成分は $2^{(j-i)/2}\overline{b_{ij}}$ である。この式は
$$
\operatorname{Dom}S=\left\{b:\sum_{i,j}2^{j-i}|b_{ij}|^2<\infty\right\}
$$
を定義域とする閉作用素を与える。実際、通常の二乗和とこの重み付き二乗和を合わせたものがグラフノルムの二乗であり、行列の有限切断は両方のノルムで近似できる。有限行列は $B(K)\xi$ に属し、$S_0$ は全定義域で同じ成分式を満たすため、その閉包はこの作用素に一致する。すると $Jb=b^*$、$\Delta E_{ij}=2^{j-i}E_{ij}$ である。$j-i$ が大きいと固有値が無限に増え、逆に $i-j$ が大きいと固有値は零に近づく。従って $\Delta$ は核が零でありながら、$\Delta$ と $\Delta^{-1}$ の両方が非有界である。
:::

## Tomita の基本定理

閉作用素 $S$ をそのまま元の両側から掛けると、積の定義域を失うことがある。そこで Rieffel と Van Daele の方法に従い、$\Delta$ の有界な関数
$$
R=(1+\Delta)^{-1},\qquad T=\Delta^{1/2}(1+\Delta)^{-1},\qquad C=TJ
$$
を用いる。$T$ は単射で稠密な値域を持つ。関数 $\sqrt t/(1+t)$ は $t\mapsto t^{-1}$ で不変なので $T$ は $J$ と可換である。従って反線形作用素 $C$ は $\inner{C u}{v}=\inner{C v}{u}$ を満たす。一方、独立に定義した $M'$ の対合は $S^*$ の制限であり、$S^*(y'\xi)=y'^*\xi$（$y'\in M'$）となる。これは閉包可能性を示した内積の等式で、反線形随伴の定義を読んだものである。$S=J\Delta^{1/2}$、$S^*=\Delta^{1/2}J$ から
$$
C y\xi=(1-R)y^*\xi,\qquad C y'\xi=R y'^*\xi\quad(y\in M,\ y'\in M')
$$
を得る。この式には $S^*$ の定義域全体を $M'\xi$ で近似できるという仮定は要らない。

:::lemma[レゾルベントを挟んだ作用素方程式]{#lem-oa-modular-resolvent}
$x'\in M'_+$ および $\Re\lambda>0$ に対し、ある $x=x^*\in M$ が
$$
C x'C=\lambda(1-R)xR+\bar\lambda R x(1-R)
$$
を満たす。
:::

:::proof
まず $y=y^*\in M$ に対する実線形汎関数 $\inner{y\xi}{x'\xi}$ を考える。自己共役元のノルム球 $\norm x\le\norm{x'}/(2\Re\lambda)$ を取り、各 $y$ に値 $2\Re(\lambda\inner{y\xi}{x\xi})$ を対応させる。この像は点ごとの位相でコンパクト凸である。実際、球は WOT コンパクトで、それぞれの値は WOT 連続である。もし目的の汎関数が像になければ、有限個の座標による実分離から、ある $y=y^*$ でその値が像の上限を越える。ところが
$$
\inner{y\xi}{x'\xi}\le\inner{\abs y\,\xi}{x'\xi}\le\norm{x'}\inner{\abs y\,\xi}{\xi}
$$
である。$x'$ は $M$ と可換なので、最初の不等式は $x'^{1/2}$ を両側から掛けた正性、後の不等式は $\abs y^{1/2}$ を両側から掛けた正性による。上限は $x=\norm{x'}\operatorname{sgn}(y)/(2\Re\lambda)$ を選べばこの値以上である。符号関数は Borel カルキュラスにより $M$ に属し、$y\operatorname{sgn}(y)=\abs y$ となる。従って分離は不可能で、目的の $x$ が存在する。

$y=y^*$ では対になった二項を共役で書ける。任意の $y$ に複素線形に延長して
$$
\inner{y\xi}{x'\xi}
=\lambda\inner{y\xi}{x\xi}+\bar\lambda\inner{x\xi}{y^*\xi}
$$
を得る。特に $\lambda=1$ なら、この式は $x\xi$ が $\operatorname{Dom} S^*$ に属し $S^*x\xi=x'\xi-x\xi$ であることを意味する。$x=x^*$ から $Sx\xi=x\xi$ であるため、$(1+\Delta)x\xi=x'\xi$、すなわち $x\xi=R x'\xi=C x'\xi$ である。正元の線形包を使えば、任意の $u'\in M'$ に対し $y\in M$ で $y\xi=C u'\xi$、$y^*\xi=C u'^*\xi$ となるものを取れる。自己共役元では正負分解を、一般元では自己共役成分の分解と $C$ の反線形性を使う。

最初の任意の $\lambda$ に戻り、等式に $z^*y$ を代入する。$x'$ と $M$ の可換性から
$$
\inner{y\xi}{x'z\xi}
=\lambda\inner{y\xi}{zx\xi}+\bar\lambda\inner{yx\xi}{z\xi}.
$$
ここで $y\xi=C u'\xi$、$z\xi=C v'\xi$ となる $y,z\in M$ を取る。左辺は $\inner{v'\xi}{Cx'C u'\xi}$ である。第 1 項の内積は $\inner{C(zx\xi)}{u'\xi}$ となり、$C(zx\xi)=(1-R)x z^*\xi=(1-R)xR v'\xi$ である。第 2 項の内積は $\inner{v'\xi}{C(yx\xi)}$ となり、同様に $C(yx\xi)=(1-R)xR u'\xi$ である。係数を第 2 引数へ移すと共役が付くため、両辺の行列係数は主張の式になる。$M'\xi$ は稠密だから作用素が一致する。
:::

:::lemma[レゾルベント方程式の積分表示]{#lem-oa-modular-fourier}
:ref[lem-oa-modular-resolvent] で $\lambda=e^{i\theta}$、$\abs\theta<\pi/2$ とすると、その解 $x$ は
$$
x=\int_\R \frac{e^{-2\theta t}}{2\cosh(\pi t)}
\,\Delta^{it}Jx'J\Delta^{-it}\,\dd t
$$
である。積分は有界作用素の行列係数で取る。
:::

:::proof
スカラー恒等式
$$
\int_\R\frac{e^{itu}}{2\cosh(\pi t)}\,\dd t=\frac1{2\cosh(u/2)}
$$
を先に示す。$u$ が実数のとき、積分路を高さ $1$ の水平線に移す。$t=i/2$ にある単純極の留数は $e^{-u/2}/(2\pi i)$ である。高さ $1$ の被積分関数は実軸のものの $-e^{-u}$ 倍になり、長方形の鉛直辺は横幅を広げると指数減衰で消える。よって実軸積分に $1+e^{-u}$ を掛けた値が $e^{-u/2}$ となる。積分は $\abs{\Im u}<\pi$ で正則なので、恒等定理でその帯まで等式が延びる。$u$ を $u+2i\theta$ に替えると、積分の核に重み $e^{-2\theta t}$ が付く。

スペクトル射影 $P_m=1_{[1/m,m]}(\Delta)$ で圧縮する。その空間では $T$ が可逆なので、前の式を両側から $T^{-1}$ で割れる。$P=P_m$、$D=\Delta|_{PH}$、$X=PxP$、$A=PJx'JP$ とすると
$$
A=e^{i\theta}D^{1/2}XD^{-1/2}+e^{-i\theta}D^{-1/2}XD^{1/2}
$$
となる。まず $\log D$ が有限スペクトルだとする。固有値 $h_i,h_j$ の間のブロックでは、$A$ は $X$ の対応ブロックに $2\cosh((h_i-h_j)/2+i\theta)$ を掛けたものである。一方、積分で $D^{it}$ と $D^{-it}$ を作用させると $e^{it(h_i-h_j)}$ が掛かる。スカラー恒等式から積分による乗数は上の乗数の逆数であり、積分の値は $X$ になる。

一般の $\log D$ はスペクトルを短い区間に分割した有限スペクトル作用素でノルム近似できる。各近似については $X$ を固定し、右辺から近似した $A$ を定める。その $A$ はノルムで収束し、虚数冪も各 $t$ に対してノルム収束する。核は積分可能で共役はノルムを保つため、優収束で $PH$ 上の式が成立する。最後に $P_m$ は強収束で $1$ に近づく。共役は $P_m$ と可換であり、両側の行列係数の極限を取れば全空間で成立する。
:::

:::theorem[Tomita の基本定理]{#thm-oa-tomita}
巡回分離ベクトルに対応する $J,\Delta$ に対し $JMJ=M'$ および $\Delta^{it}M\Delta^{-it}=M$（$t\in\R$）である。
:::
:::proof
$x'\in M'_+$、$b'\in M'$ を固定する。$\Delta^{it}Jx'J\Delta^{-it}$ と $b'$ の交換子の任意の行列係数を $g(t)$ とする。これは有界連続である。:ref[lem-oa-modular-fourier] の積分の値は $M$ に属するので、$\int_\R e^{-2\theta t}g(t)/(2\cosh(\pi t))\,\dd t=0$ が全 $\abs\theta<\pi/2$ で成立する。指数の係数 $2\theta$ を複素数 $z$ に替えると積分は $\abs{\Re z}<\pi$ で正則である。実区間で消えるため恒等定理により虚軸でも消える。これは $g(t)/(2\cosh(\pi t))$ の Fourier 変換が零という意味であり、Fourier 変換の単射性と $g$ の連続性から $g=0$ である。全行列係数で成立し、$M'$ の正元は線形に $M'$ を張る。従って全 $x'\in M'$ で $\Delta^{it}Jx'J\Delta^{-it}\in M$ である。$t=0$ で $JM'J\subset M$ を得る。

逆向きの可換性には、まず $\Delta\xi=\xi$、$J\xi=\xi$ を使う。$S\xi=S^*\xi=\xi$ は各対合の定義で分かり、極分解に代入すれば両式が出る。$x,y\in M$ が自己共役なら $Jx\xi=\Delta^{1/2}x\xi$ と $Jy\xi=\Delta^{1/2}y\xi$ である。$\Delta^{1/2}$ の自己共役性で $\inner{Jx\xi}{y\xi}=\inner{x\xi}{Jy\xi}$ であり、
$$
\inner{y(JxJ)\xi}{\xi}=\inner{\xi}{x(JyJ)\xi}
$$
と書ける。両辺は $x$ に共役線形、$y$ に線形なので一般の $x,y$ にも成り立つ。$b'\in M'$ として $y$ を $y(Jb'J)$ に替えてよい。元の $x,y$ を再び自己共役とすると、置換後の左辺は
$$
\inner{J b' x\xi}{y\xi}=\inner{Jy\xi}{b'x\xi}
=\inner{x(JyJ)\xi}{b'\xi}
$$
であり、右辺は $\inner{(JyJ)x\xi}{b'\xi}$ になる。$M'\xi$ は稠密だから $x(JyJ)\xi=(JyJ)x\xi$ であり、自己共役成分に分ければ全 $x,y$ で成り立つ。任意の $w\in M$ を使い、$xw$ での等式と $w$ での等式を比べると $x$ と $JyJ$ は $w\xi$ 上で可換である。巡回性で $JyJ\in M'$ となり、$JMJ\subset M'$ を得る。両包含と $J^2=1$ によって等号である。先に得た $\Delta^{it}JM'J\Delta^{-it}\subset M$ に適用し、$-t$ にも同じ包含を用いれば $\Delta^{it}M\Delta^{-it}=M$ となる。
:::

### 行列模型で二つの結論を直接検査する

$M=L(M_d)$ とする。$JL_xJ(a)=(xa^*)^*=ax^*$ なので $JMJ=R(M_d)$ である。逆に左作用と可換な $T$ は $T(a)=T(L_a1)=aT(1)$ となるので、この右作用に尽きる。さらに
$$
\Delta^{it}L_x\Delta^{-it}(a)=D^{it}xD^{-it}a
$$
である。ここで第 2 因子の $D^{it}$ と $D^{-it}$ は打ち消し合う。従って $\Delta^{it}$ 自身の作用は左右からの積だが、左作用への共役は左側の密度だけを残す。一般の証明は無限次元の作用素でこの打ち消しを保証する。特に $\Delta$ が非有界な場合、$M\xi$ における形式計算の延長には上で用いた稠密性とレゾルベントの議論が要る。

:::example[非対角の位相を半回転させる]
密度を $D=c^{-1}\diag(b^{n_j})$ とする。時刻 $t=\pi r/\log b$ での $\Delta^{it}L_a\Delta^{-it}$ を一つの行列単位上で直接計算し、左からの作用になるか調べる。

$D^{it}$ には全座標に共通の $c^{-it}$ が付く。共役の積ではこれを打ち消せるため、位相だけの対角 $U=\diag(e^{i\pi r n_j})$ を使う。

$$
\begin{aligned}
(b,r)&=\left( 2, \  1\right)\\
D&=\begin{pmatrix}\frac{1}{3} & 0\\0 & \frac{2}{3}\end{pmatrix}\\
U&=\begin{pmatrix}1 & 0\\0 & -1\end{pmatrix}\\
a&=\begin{pmatrix}1 & 2\\3 & 4\end{pmatrix}
\end{aligned}
$$

$\Delta^{-it}(v)=U^*vU$ を先に計算する。その後 $a$ を左から掛け、最後にもう一度共役する。この順序を行列行列積として示す。

$$
\begin{aligned}
v&=\begin{pmatrix}0 & 1\\0 & 0\end{pmatrix}\\
\Delta^{-it}(v)&=\begin{pmatrix}0 & -1\\0 & 0\end{pmatrix}\\
L_a\Delta^{-it}(v)&=\begin{pmatrix}0 & -1\\0 & -3\end{pmatrix}
\end{aligned}
$$

残る右側の二ユニタリは互いに打ち消し合う。従って同じ結果は $\sigma_t(a)=UaU^*$ を $v$ の左に掛けても得られる。

$$
\begin{aligned}
\sigma_t(a)&=\begin{pmatrix}1 & -2\\-3 & 4\end{pmatrix}\\
\sigma_t(a)v&=\begin{pmatrix}0 & 1\\0 & -3\end{pmatrix}
\end{aligned}
$$

この計算で任意の $v$ に対する式も見える。一般の Tomita 定理では、モジュラー作用素による共役が同じ環へ戻ることを解析的に証明する。
:::

:::example[四分の一回転での左作用]
密度を $D=c^{-1}\diag(b^{n_j})$ とする。時刻 $t=\pi r/\log b$ での $\Delta^{it}L_a\Delta^{-it}$ を一つの行列単位上で直接計算し、左からの作用になるか調べる。

$D^{it}$ には全座標に共通の $c^{-it}$ が付く。共役の積ではこれを打ち消せるため、位相だけの対角 $U=\diag(e^{i\pi r n_j})$ を使う。

$$
\begin{aligned}
(b,r)&=\left( 3, \  \frac{1}{2}\right)\\
D&=\begin{pmatrix}\frac{1}{4} & 0\\0 & \frac{3}{4}\end{pmatrix}\\
U&=\begin{pmatrix}1 & 0\\0 & i\end{pmatrix}\\
a&=\begin{pmatrix}2 & 1\\1 & 3\end{pmatrix}
\end{aligned}
$$

$\Delta^{-it}(v)=U^*vU$ を先に計算する。その後 $a$ を左から掛け、最後にもう一度共役する。この順序を行列行列積として示す。

$$
\begin{aligned}
v&=\begin{pmatrix}0 & 1\\0 & 0\end{pmatrix}\\
\Delta^{-it}(v)&=\begin{pmatrix}0 & i\\0 & 0\end{pmatrix}\\
L_a\Delta^{-it}(v)&=\begin{pmatrix}0 & 2 i\\0 & i\end{pmatrix}
\end{aligned}
$$

残る右側の二ユニタリは互いに打ち消し合う。従って同じ結果は $\sigma_t(a)=UaU^*$ を $v$ の左に掛けても得られる。

$$
\begin{aligned}
\sigma_t(a)&=\begin{pmatrix}2 & - i\\i & 3\end{pmatrix}\\
\sigma_t(a)v&=\begin{pmatrix}0 & 2\\0 & i\end{pmatrix}
\end{aligned}
$$

この計算で任意の $v$ に対する式も見える。一般の Tomita 定理では、モジュラー作用素による共役が同じ環へ戻ることを解析的に証明する。
:::

:::example[二段離れた位相が戻る]
密度を $D=c^{-1}\diag(b^{n_j})$ とする。時刻 $t=\pi r/\log b$ での $\Delta^{it}L_a\Delta^{-it}$ を一つの行列単位上で直接計算し、左からの作用になるか調べる。

$D^{it}$ には全座標に共通の $c^{-it}$ が付く。共役の積ではこれを打ち消せるため、位相だけの対角 $U=\diag(e^{i\pi r n_j})$ を使う。

$$
\begin{aligned}
(b,r)&=\left( 2, \  1\right)\\
D&=\begin{pmatrix}\frac{1}{7} & 0 & 0\\0 & \frac{2}{7} & 0\\0 & 0 & \frac{4}{7}\end{pmatrix}\\
U&=\begin{pmatrix}1 & 0 & 0\\0 & -1 & 0\\0 & 0 & 1\end{pmatrix}\\
a&=\begin{pmatrix}1 & 1 & 2\\1 & 2 & 3\\2 & 3 & 4\end{pmatrix}
\end{aligned}
$$

$\Delta^{-it}(v)=U^*vU$ を先に計算する。その後 $a$ を左から掛け、最後にもう一度共役する。この順序を行列行列積として示す。

$$
\begin{aligned}
v&=\begin{pmatrix}0 & 0 & 1\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix}\\
\Delta^{-it}(v)&=\begin{pmatrix}0 & 0 & 1\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix}\\
L_a\Delta^{-it}(v)&=\begin{pmatrix}0 & 0 & 1\\0 & 0 & 1\\0 & 0 & 2\end{pmatrix}
\end{aligned}
$$

残る右側の二ユニタリは互いに打ち消し合う。従って同じ結果は $\sigma_t(a)=UaU^*$ を $v$ の左に掛けても得られる。

$$
\begin{aligned}
\sigma_t(a)&=\begin{pmatrix}1 & -1 & 2\\-1 & 2 & -3\\2 & -3 & 4\end{pmatrix}\\
\sigma_t(a)v&=\begin{pmatrix}0 & 0 & 1\\0 & 0 & -1\\0 & 0 & 2\end{pmatrix}
\end{aligned}
$$

この計算で任意の $v$ に対する式も見える。一般の Tomita 定理では、モジュラー作用素による共役が同じ環へ戻ることを解析的に証明する。
:::

:::example[密度が同じブロック内での作用]
密度を $D=c^{-1}\diag(b^{n_j})$ とする。時刻 $t=\pi r/\log b$ での $\Delta^{it}L_a\Delta^{-it}$ を一つの行列単位上で直接計算し、左からの作用になるか調べる。

$D^{it}$ には全座標に共通の $c^{-it}$ が付く。共役の積ではこれを打ち消せるため、位相だけの対角 $U=\diag(e^{i\pi r n_j})$ を使う。

$$
\begin{aligned}
(b,r)&=\left( 2, \  1\right)\\
D&=\begin{pmatrix}\frac{1}{4} & 0 & 0\\0 & \frac{1}{4} & 0\\0 & 0 & \frac{1}{2}\end{pmatrix}\\
U&=\begin{pmatrix}1 & 0 & 0\\0 & 1 & 0\\0 & 0 & -1\end{pmatrix}\\
a&=\begin{pmatrix}0 & 1 & 2\\1 & 0 & 3\\2 & 3 & 4\end{pmatrix}
\end{aligned}
$$

$\Delta^{-it}(v)=U^*vU$ を先に計算する。その後 $a$ を左から掛け、最後にもう一度共役する。この順序を行列行列積として示す。

$$
\begin{aligned}
v&=\begin{pmatrix}0 & 0 & 1\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix}\\
\Delta^{-it}(v)&=\begin{pmatrix}0 & 0 & -1\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix}\\
L_a\Delta^{-it}(v)&=\begin{pmatrix}0 & 0 & 0\\0 & 0 & -1\\0 & 0 & -2\end{pmatrix}
\end{aligned}
$$

残る右側の二ユニタリは互いに打ち消し合う。従って同じ結果は $\sigma_t(a)=UaU^*$ を $v$ の左に掛けても得られる。

$$
\begin{aligned}
\sigma_t(a)&=\begin{pmatrix}0 & 1 & -2\\1 & 0 & -3\\-2 & -3 & 4\end{pmatrix}\\
\sigma_t(a)v&=\begin{pmatrix}0 & 0 & 0\\0 & 0 & 1\\0 & 0 & -2\end{pmatrix}
\end{aligned}
$$

この計算で任意の $v$ に対する式も見える。一般の Tomita 定理では、モジュラー作用素による共役が同じ環へ戻ることを解析的に証明する。
:::

:::example[二段目と複素成分]
密度を $D=c^{-1}\diag(b^{n_j})$ とする。時刻 $t=\pi r/\log b$ での $\Delta^{it}L_a\Delta^{-it}$ を一つの行列単位上で直接計算し、左からの作用になるか調べる。

$D^{it}$ には全座標に共通の $c^{-it}$ が付く。共役の積ではこれを打ち消せるため、位相だけの対角 $U=\diag(e^{i\pi r n_j})$ を使う。

$$
\begin{aligned}
(b,r)&=\left( 3, \  \frac{1}{2}\right)\\
D&=\begin{pmatrix}\frac{1}{10} & 0\\0 & \frac{9}{10}\end{pmatrix}\\
U&=\begin{pmatrix}1 & 0\\0 & -1\end{pmatrix}\\
a&=\begin{pmatrix}2 & 1 + i\\1 - i & 3\end{pmatrix}
\end{aligned}
$$

$\Delta^{-it}(v)=U^*vU$ を先に計算する。その後 $a$ を左から掛け、最後にもう一度共役する。この順序を行列行列積として示す。

$$
\begin{aligned}
v&=\begin{pmatrix}0 & 1\\0 & 0\end{pmatrix}\\
\Delta^{-it}(v)&=\begin{pmatrix}0 & -1\\0 & 0\end{pmatrix}\\
L_a\Delta^{-it}(v)&=\begin{pmatrix}0 & -2\\0 & -1 + i\end{pmatrix}
\end{aligned}
$$

残る右側の二ユニタリは互いに打ち消し合う。従って同じ結果は $\sigma_t(a)=UaU^*$ を $v$ の左に掛けても得られる。

$$
\begin{aligned}
\sigma_t(a)&=\begin{pmatrix}2 & -1 - i\\-1 + i & 3\end{pmatrix}\\
\sigma_t(a)v&=\begin{pmatrix}0 & 2\\0 & -1 + i\end{pmatrix}
\end{aligned}
$$

この計算で任意の $v$ に対する式も見える。一般の Tomita 定理では、モジュラー作用素による共役が同じ環へ戻ることを解析的に証明する。
:::

## モジュラー自己同型群

$\varphi(x)=\inner{x\xi}{\xi}$ に付随する :term[モジュラー自己同型群]{en="modular automorphism group" yomi="もじゅらーじこどうけいぐん"} は $\sigma_t^\varphi(x)=\Delta^{it}x\Delta^{-it}$ である。:ref[thm-oa-tomita] で各写像が $M$ の自己同型となり、スペクトル積分から $\Delta^{it}\Delta^{is}=\Delta^{i(t+s)}$ である。$\Delta$ の核は零なので $\lambda^{it}$ がスペクトル上ほとんど至る所で定まり、絶対値が $1$ である。優収束で $t\mapsto\Delta^{it}$ は強連続で、共役作用も各 $x$ で強連続である。$S\xi=S^*\xi=\xi$ より $\Delta\xi=\xi$ であり、$\varphi$ は群で不変である。ここで $S^*\xi=\xi$ は閉包可能性の証明の式で $y=1$ として得られる。

:::example[行列単位の周波数]
行列模型で $\sigma_t^\varphi(a)=D^{it}aD^{-it}$ である。$D=\diag(1/4,3/4)$ の $E_{12}$ には $e^{-it\log3}$ が掛かり、対角行列は固定される。トレース状態では $D$ はスカラーなので全て固定される。一般に巡回分離ベクトルの状態がトレースなら $\norm{x\xi}=\norm{x^*\xi}$ であるため $S$ は反ユニタリ、$\Delta=1$ となる。
:::

:::example[有限因子環でのトレースの変更]
有限因子環の忠実正規トレースを $\tau$ とし、$\tau(D)=1$ で正かつ可逆な $D\in M$ を取る。トレースの GNS 空間において $\xi=D^{1/2}$ を使う。$xD^{1/2}$ は $M$ を動くため巡回的である。$\norm{xD^{1/2}}_2=0$ なら $\tau(D^{1/2}x^*xD^{1/2})=0$ で、忠実性と可逆性により $x=0$ となるため分離的でもある。右乗法 $R_y$ はトレース性から $\norm{xy}_2\le\norm y\norm x_2$ を満たし、$J_0(x)=x^*$ は反ユニタリである。

$S(x)=D^{-1/2}x^*D^{1/2}$ と計算でき、これは $J_0 L_{D^{1/2}}R_{D^{-1/2}}$ である。正部分は可換な正作用素の積だから $J=J_0$、$\Delta=L_D R_{D^{-1}}$ となる。よって有限次元と同じ $\sigma_t^\varphi(a)=D^{it}aD^{-it}$ を得る。例えば $\tau(p)=1/2$ の射影 $p$ を取り $D=(3/2)p+(1/2)(1-p)$ とすると $\tau(D)=1$ である。$pM(1-p)$ では位相因子 $3^{it}$、$(1-p)Mp$ では $3^{-it}$ が掛かる。従ってこの状態の固定点は二つの対角の角の直和である。
:::

### 中心化環

状態のモジュラー自己同型群で固定される部分環を :term[中心化環]{en="centralizer" yomi="ちゅうしんかかん"} という。行列模型では固定性は $xD=Dx$ に等しい。$D$ の固有空間間の行列成分に対する倍率 $(d_i/d_j)^{it}$ がすべての $t$ で $1$ となるのは $d_i=d_j$ の場合だけである。また全 $y$ に対する $\varphi(xy)=\varphi(yx)$ は $\tr((Dx-xD)y)=0$ と同値であり、行列単位で検査すればやはり $Dx=xD$ が出る。

:::example[固有値の重複による固定点の増加]
$D=\diag(1/4,1/4,1/2)$ なら固定点は $M_2\oplus\C$ である。例えば $E_{12}$ は固定されるが、$E_{13}$ は $2^{-it}E_{13}$ に変わる。その時間平均の係数は
$$
\frac1{2T}\int_{-T}^T2^{-it}\dd t
=\frac{\sin(T\log2)}{T\log2}
$$
である。$T\to\infty$ でこの成分は消える。各行列単位で同じ計算を行えば、時間平均が固有空間ごとのブロック期待値になる。
:::

:::example[三つの固有空間の周期平均]
行列の密度の固有値を $c^{-1}b^{n_j}$ とする。モジュラー自己同型群をちょうど一周期平均すると、どの行列成分が残るか求める。

$\sigma_t(E_{ij})=e^{it(n_i-n_j)\log b}E_{ij}$ である。周波数を $\log b$ で割った整数の表にすると、平均すべき成分が見える。

$$
\begin{aligned}
(b,n_1,n_2,n_3)&=\left( 2, \  0, \  1, \  2\right)\\
((n_i-n_j))_{ij}&=\begin{pmatrix}0 & -1 & -2\\1 & 0 & -1\\2 & 1 & 0\end{pmatrix}
\end{aligned}
$$

$2\pi/\log b$ の区間で積分すると、整数周波数が非零なら始点と終点の位相が同じで積分が零になり、零なら積分を区間長で割った値は $1$ になる。

$$
\begin{aligned}
a&=\begin{pmatrix}1 & 2 & 1\\2 & 2 & 3\\1 & 3 & 3\end{pmatrix}\\
\frac{\log b}{2\pi}\int_0^{2\pi/\log b}\sigma_t(a)\,\dd t&=\begin{pmatrix}1 & 0 & 0\\0 & 2 & 0\\0 & 0 & 3\end{pmatrix}
\end{aligned}
$$

平均は密度の固有空間ごとのブロックをそのまま残す。各ブロックでの行と列を数えれば、固定点環の線形次元は多重度の二乗和になる。

$$
\begin{aligned}
\dim M^{\sigma^\varphi}&=3
\end{aligned}
$$

時間平均が全対角への圧縮になるのは固有値が重複しない場合である。重複があると非対角でも周波数が零になり、平均で消えない。
:::

:::example[二次元ブロックを残す平均]
行列の密度の固有値を $c^{-1}b^{n_j}$ とする。モジュラー自己同型群をちょうど一周期平均すると、どの行列成分が残るか求める。

$\sigma_t(E_{ij})=e^{it(n_i-n_j)\log b}E_{ij}$ である。周波数を $\log b$ で割った整数の表にすると、平均すべき成分が見える。

$$
\begin{aligned}
(b,n_1,n_2,n_3)&=\left( 2, \  0, \  0, \  1\right)\\
((n_i-n_j))_{ij}&=\begin{pmatrix}0 & 0 & -1\\0 & 0 & -1\\1 & 1 & 0\end{pmatrix}
\end{aligned}
$$

$2\pi/\log b$ の区間で積分すると、整数周波数が非零なら始点と終点の位相が同じで積分が零になり、零なら積分を区間長で割った値は $1$ になる。

$$
\begin{aligned}
a&=\begin{pmatrix}1 & 2 & 1\\2 & 2 & 3\\1 & 3 & 3\end{pmatrix}\\
\frac{\log b}{2\pi}\int_0^{2\pi/\log b}\sigma_t(a)\,\dd t&=\begin{pmatrix}1 & 2 & 0\\2 & 2 & 0\\0 & 0 & 3\end{pmatrix}
\end{aligned}
$$

平均は密度の固有空間ごとのブロックをそのまま残す。各ブロックでの行と列を数えれば、固定点環の線形次元は多重度の二乗和になる。

$$
\begin{aligned}
\dim M^{\sigma^\varphi}&=5
\end{aligned}
$$

時間平均が全対角への圧縮になるのは固有値が重複しない場合である。重複があると非対角でも周波数が零になり、平均で消えない。
:::

:::example[トレースでは時間平均が恒等]
行列の密度の固有値を $c^{-1}b^{n_j}$ とする。モジュラー自己同型群をちょうど一周期平均すると、どの行列成分が残るか求める。

$\sigma_t(E_{ij})=e^{it(n_i-n_j)\log b}E_{ij}$ である。周波数を $\log b$ で割った整数の表にすると、平均すべき成分が見える。

$$
\begin{aligned}
(b,n_1,n_2,n_3)&=\left( 3, \  1, \  1, \  1\right)\\
((n_i-n_j))_{ij}&=\begin{pmatrix}0 & 0 & 0\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix}
\end{aligned}
$$

$2\pi/\log b$ の区間で積分すると、整数周波数が非零なら始点と終点の位相が同じで積分が零になり、零なら積分を区間長で割った値は $1$ になる。

$$
\begin{aligned}
a&=\begin{pmatrix}1 & 2 & 1\\2 & 2 & 3\\1 & 3 & 3\end{pmatrix}\\
\frac{\log b}{2\pi}\int_0^{2\pi/\log b}\sigma_t(a)\,\dd t&=\begin{pmatrix}1 & 2 & 1\\2 & 2 & 3\\1 & 3 & 3\end{pmatrix}
\end{aligned}
$$

平均は密度の固有空間ごとのブロックをそのまま残す。各ブロックでの行と列を数えれば、固定点環の線形次元は多重度の二乗和になる。

$$
\begin{aligned}
\dim M^{\sigma^\varphi}&=9
\end{aligned}
$$

時間平均が全対角への圧縮になるのは固有値が重複しない場合である。重複があると非対角でも周波数が零になり、平均で消えない。
:::

:::example[離れた座標をまとめる固定点]
行列の密度の固有値を $c^{-1}b^{n_j}$ とする。モジュラー自己同型群をちょうど一周期平均すると、どの行列成分が残るか求める。

$\sigma_t(E_{ij})=e^{it(n_i-n_j)\log b}E_{ij}$ である。周波数を $\log b$ で割った整数の表にすると、平均すべき成分が見える。

$$
\begin{aligned}
(b,n_1,n_2,n_3)&=\left( 3, \  0, \  2, \  0\right)\\
((n_i-n_j))_{ij}&=\begin{pmatrix}0 & -2 & 0\\2 & 0 & 2\\0 & -2 & 0\end{pmatrix}
\end{aligned}
$$

$2\pi/\log b$ の区間で積分すると、整数周波数が非零なら始点と終点の位相が同じで積分が零になり、零なら積分を区間長で割った値は $1$ になる。

$$
\begin{aligned}
a&=\begin{pmatrix}1 & 2 & 1\\2 & 2 & 3\\1 & 3 & 3\end{pmatrix}\\
\frac{\log b}{2\pi}\int_0^{2\pi/\log b}\sigma_t(a)\,\dd t&=\begin{pmatrix}1 & 0 & 1\\0 & 2 & 0\\1 & 0 & 3\end{pmatrix}
\end{aligned}
$$

平均は密度の固有空間ごとのブロックをそのまま残す。各ブロックでの行と列を数えれば、固定点環の線形次元は多重度の二乗和になる。

$$
\begin{aligned}
\dim M^{\sigma^\varphi}&=5
\end{aligned}
$$

時間平均が全対角への圧縮になるのは固有値が重複しない場合である。重複があると非対角でも周波数が零になり、平均で消えない。
:::

:::example[最大周波数が最小の三倍になる]
行列の密度の固有値を $c^{-1}b^{n_j}$ とする。モジュラー自己同型群をちょうど一周期平均すると、どの行列成分が残るか求める。

$\sigma_t(E_{ij})=e^{it(n_i-n_j)\log b}E_{ij}$ である。周波数を $\log b$ で割った整数の表にすると、平均すべき成分が見える。

$$
\begin{aligned}
(b,n_1,n_2,n_3)&=\left( 2, \  0, \  1, \  3\right)\\
((n_i-n_j))_{ij}&=\begin{pmatrix}0 & -1 & -3\\1 & 0 & -2\\3 & 2 & 0\end{pmatrix}
\end{aligned}
$$

$2\pi/\log b$ の区間で積分すると、整数周波数が非零なら始点と終点の位相が同じで積分が零になり、零なら積分を区間長で割った値は $1$ になる。

$$
\begin{aligned}
a&=\begin{pmatrix}1 & 2 & 1\\2 & 2 & 3\\1 & 3 & 3\end{pmatrix}\\
\frac{\log b}{2\pi}\int_0^{2\pi/\log b}\sigma_t(a)\,\dd t&=\begin{pmatrix}1 & 0 & 0\\0 & 2 & 0\\0 & 0 & 3\end{pmatrix}
\end{aligned}
$$

平均は密度の固有空間ごとのブロックをそのまま残す。各ブロックでの行と列を数えれば、固定点環の線形次元は多重度の二乗和になる。

$$
\begin{aligned}
\dim M^{\sigma^\varphi}&=3
\end{aligned}
$$

時間平均が全対角への圧縮になるのは固有値が重複しない場合である。重複があると非対角でも周波数が零になり、平均で消えない。
:::

## KMS 条件

:::definition[上半分の帯での KMS 条件]{#def-oa-kms}
一径数自己同型群 $\alpha$ に対し、状態 $\varphi$ が逆温度 $\beta>0$ の :term[KMS 条件]{en="KMS condition" yomi="けーえむえすじょうけん"} を満たすとは、各 $a,b$ に対し $0\le\Im z\le\beta$ で有界連続、内部で正則な $F$ があり、
$$
F(t)=\varphi(a\alpha_t(b)),\qquad F(t+i\beta)=\varphi(\alpha_t(b)a)
$$
となることである。
:::

ここでは von Neumann 環上の群を σ 弱連続とし、状態を正規とする。C*-環で同じ定義を使う場合は、通常は各元の軌道にノルム連続性を仮定する。有限次元で Hamiltonian を $H$、Gibbs 密度を $D=e^{-\beta H}/\tr(e^{-\beta H})$ とすると、$D^{-it/\beta}aD^{it/\beta}=e^{itH}ae^{-itH}$ である。物理時間を用いる場合は $t$ を $t/\hbar$ に置き換えるため、帯の幅も $\beta$ から $\hbar\beta$ になる。

:::proposition[モジュラー状態の KMS 性]{#prop-oa-modular-kms}
この上側の帯と順序を用いると、巡回分離単位ベクトル $\xi$ の状態 $\varphi$ は $\alpha_t=\sigma_{-t}^\varphi$ に関し $\beta=1$ の KMS 状態である。
:::

:::proof
:ref[thm-fa-unbounded-spectral] で得る $\Delta$ のスペクトル測度を $E$ とし、複素測度 $\nu(B)=\inner{E(B)b\xi}{a^*\xi}$ を取る。$b\xi,a^*\xi$ はともに $\Delta^{1/2}$ の定義域に属する。分割による Cauchy–Schwarz を使うと、$\int\lambda\,\dd\abs\nu\le\norm{\Delta^{1/2}b\xi}\norm{\Delta^{1/2}a^*\xi}<\infty$、同様に $\int\dd\abs\nu<\infty$ である。従って $F(z)=\int\lambda^{-iz}\,\dd\nu(\lambda)$ は $0\le\Im z\le1$ で有界連続である。$\lambda=0$ に質量はない。内部では $\lambda^{\Im z}\abs{\log\lambda}$ が $1+\lambda$ の定数倍で抑えられるため、コンパクトな部分帯上で微分でき、正則である。

下の境界で $F(t)=\inner{\Delta^{-it}b\xi}{a^*\xi}=\varphi(a\alpha_t(b))$ となる。上の境界では
$$
F(t+i)=\inner{\Delta^{-it}\Delta^{1/2}b\xi}{\Delta^{1/2}a^*\xi}.
$$
$S=J\Delta^{1/2}$ から $\Delta^{1/2}b\xi=Jb^*\xi$ と $\Delta^{1/2}a^*\xi=Ja\xi$ である。また $J$ は反線形なので $J\Delta J=\Delta^{-1}$ は $J\Delta^{it}J=\Delta^{it}$ を導く。よって上式は $\inner{a\xi}{\Delta^{-it}b^*\xi}=\varphi(\alpha_t(b)a)$ となる。
:::

行列環では同じ値を直接計算できる。$\alpha_{t+i}(b)=D\alpha_t(b)D^{-1}$ なので、$\tr(DaD\alpha_t(b)D^{-1})=\tr(D\alpha_t(b)a)$ となる。

符号を検算するには $a=E_{21},b=E_{12}$ と置く。$F(0)=d_2$ に対し $F(i)=(d_1/d_2)d_2=d_1$ となる。規約を替えて $\sigma_t$ を直接使うときは下側の帯を用いる。

### KMS 条件が時間の尺度を決める

:::proposition[忠実正規 KMS 状態が定める自己同型群]{#prop-oa-kms-uniqueness}
忠実正規状態 $\varphi$ が、σ 弱連続な自己同型群 $\alpha$ に関して上側の帯の幅 $\beta>0$ で KMS 条件を満たすとする。このとき $\alpha_t=\sigma^\varphi_{-t/\beta}$ である。
:::

:::proof
まず $a=1$ とすると上下の境界値はともに $\varphi(\alpha_t(b))$ である。この帯を虚方向に整数回貼り合わせれば、境界での連続性と Morera の定理によって有界な整関数になる。Liouville の定理により定数で、$\varphi\circ\alpha_t=\varphi$ である。GNS 空間で $U_tx\xi=\alpha_t(x)\xi$ と置くと等長作用素で、逆は $U_{-t}$ である。各 $x\xi$ 上の強連続性は $\norm{U_tx\xi-x\xi}^2=2\varphi(x^*x)-2\Re\varphi(x^*\alpha_t(x))$ から従う。稠密性と一様有界性で全空間にも広がる。:ref[thm-fa-stone]{text="Stone の定理"} で $U_t=e^{itL}$ と書く。

$r>0$ に対して $b_r=(4\pi r)^{-1/2}\int_{\R}e^{-t^2/(4r)}\alpha_t(b)\,\dd t$ を σ 弱積分で取る。核を $e^{-(t-z)^2/(4r)}$ に替えると、コンパクトな $z$ の集合上で積分可能な共通上界を持つので $\alpha_z(b_r)$ の整関数延長を与える。閉帯では核の絶対値の積分が $e^{(\Im z)^2/(4r)}$ 以下なので、延長は帯全体で有界である。スペクトル定理と Gauss 核の Fourier 変換から $b_r\xi=e^{-rL^2}b\xi$、$\alpha_{i\beta}(b_r)\xi=e^{-\beta L}b_r\xi$ である。帯での解析関数は下境界での値から一意である。差が下境界で零なら Schwarz の反射で延長でき、恒等定理が適用できるからである。従って KMS 条件を上端で使うと
$$
\inner{e^{-\beta L}b_r\xi}{a^*\xi}
=\varphi(b_ra)
=\inner{\Delta^{1/2}b_r\xi}{\Delta^{1/2}a^*\xi}
$$
となる。最後の等式では $Sx\xi=x^*\xi$ と $S=J\Delta^{1/2}$ を使った。$M\xi$ は $\Delta^{1/2}$ のグラフノルムで稠密な定義域なので、随伴の定義から $b_r\xi\in\operatorname{Dom}\Delta$ であり $\Delta b_r\xi=e^{-\beta L}b_r\xi$ となる。

最後に $e^{-rL^2}$ と $e^{-\beta L}e^{-rL^2}$ は有界なので、各 $r$ について稠密な $M\xi$ を全空間へ置き換えても、作られるベクトルは $e^{-\beta L}$ のグラフノルムで近似される。$r\downarrow0$ で $e^{-rL^2}$ は $e^{-\beta L}$ の定義域上でグラフノルム収束する。よって $\Delta$ は $e^{-\beta L}$ のグラフ全体を含む。両方が自己共役なので等しい。従って $U_t=\Delta^{-it/\beta}$ となり、共役により自己同型群の式が得られる。
:::

:::example[上側で小さい密度へ切り替わる]
忠実な対角密度 $D$ の状態を取り、$\alpha_t=\sigma^\varphi_{-t}$ の上側の帯で KMS 相関を計算する。

$a=E_{ji}$、$b=E_{ij}$ とする。$ab=E_{jj}$、$ba=E_{ii}$ なので、二つの順序は異なる準位の密度を読む。

$$
\begin{aligned}
D&=\begin{pmatrix}\frac{1}{4} & 0\\0 & \frac{3}{4}\end{pmatrix}\\
(i,j)&=\left( 1, \  2\right)\\
\varphi(ab)&=\frac{3}{4}\\
\varphi(ba)&=\frac{1}{4}
\end{aligned}
$$

$D^{-iz}bD^{iz}=(d_i/d_j)^{-iz}b$ である。 そこで下境界では $F(t)$ は $d_j$ にこの位相を掛けたものになる。上端 $i1$ では位相の指数が $1$ になり、比 $d_i/d_j$ が掛かる。

$$
\begin{aligned}
d_i/d_j&=\frac{1}{3}\\
F(0)&=\frac{3}{4}\\
F(i1)&=\frac{1}{4}
\end{aligned}
$$

時間 $t$ を加えても双方に同じ実時間の位相が掛かるので、$F(t+i1)=\varphi(\alpha_t(b)a)$ になる。有限次元の行列単位は線形基底なので、一般の $a,b$ に対しては各項を足せばよい。

下側と上側の値を比較するには、解析延長の位相だけでなく、境界で演算子を掛ける順序の変化も追う必要がある。
:::

:::example[逆の行列単位で上側の値が増える]
忠実な対角密度 $D$ の状態を取り、$\alpha_t=\sigma^\varphi_{-t}$ の上側の帯で KMS 相関を計算する。

$a=E_{ji}$、$b=E_{ij}$ とする。$ab=E_{jj}$、$ba=E_{ii}$ なので、二つの順序は異なる準位の密度を読む。

$$
\begin{aligned}
D&=\begin{pmatrix}\frac{3}{4} & 0\\0 & \frac{1}{4}\end{pmatrix}\\
(i,j)&=\left( 1, \  2\right)\\
\varphi(ab)&=\frac{1}{4}\\
\varphi(ba)&=\frac{3}{4}
\end{aligned}
$$

$D^{-iz}bD^{iz}=(d_i/d_j)^{-iz}b$ である。 そこで下境界では $F(t)$ は $d_j$ にこの位相を掛けたものになる。上端 $i1$ では位相の指数が $1$ になり、比 $d_i/d_j$ が掛かる。

$$
\begin{aligned}
d_i/d_j&=3\\
F(0)&=\frac{1}{4}\\
F(i1)&=\frac{3}{4}
\end{aligned}
$$

時間 $t$ を加えても双方に同じ実時間の位相が掛かるので、$F(t+i1)=\varphi(\alpha_t(b)a)$ になる。有限次元の行列単位は線形基底なので、一般の $a,b$ に対しては各項を足せばよい。

下側と上側の値を比較するには、解析延長の位相だけでなく、境界で演算子を掛ける順序の変化も追う必要がある。
:::

:::example[三準位間の遷移の相関]
忠実な対角密度 $D$ の状態を取り、$\alpha_t=\sigma^\varphi_{-t}$ の上側の帯で KMS 相関を計算する。

$a=E_{ji}$、$b=E_{ij}$ とする。$ab=E_{jj}$、$ba=E_{ii}$ なので、二つの順序は異なる準位の密度を読む。

$$
\begin{aligned}
D&=\begin{pmatrix}\frac{1}{7} & 0 & 0\\0 & \frac{2}{7} & 0\\0 & 0 & \frac{4}{7}\end{pmatrix}\\
(i,j)&=\left( 1, \  3\right)\\
\varphi(ab)&=\frac{4}{7}\\
\varphi(ba)&=\frac{1}{7}
\end{aligned}
$$

$D^{-iz}bD^{iz}=(d_i/d_j)^{-iz}b$ である。 そこで下境界では $F(t)$ は $d_j$ にこの位相を掛けたものになる。上端 $i1$ では位相の指数が $1$ になり、比 $d_i/d_j$ が掛かる。

$$
\begin{aligned}
d_i/d_j&=\frac{1}{4}\\
F(0)&=\frac{4}{7}\\
F(i1)&=\frac{1}{7}
\end{aligned}
$$

時間 $t$ を加えても双方に同じ実時間の位相が掛かるので、$F(t+i1)=\varphi(\alpha_t(b)a)$ になる。有限次元の行列単位は線形基底なので、一般の $a,b$ に対しては各項を足せばよい。

下側と上側の値を比較するには、解析延長の位相だけでなく、境界で演算子を掛ける順序の変化も追う必要がある。
:::

:::example[同じ重みを持つ二つの準位]
忠実な対角密度 $D$ の状態を取り、$\alpha_t=\sigma^\varphi_{-t}$ の上側の帯で KMS 相関を計算する。

$a=E_{ji}$、$b=E_{ij}$ とする。$ab=E_{jj}$、$ba=E_{ii}$ なので、二つの順序は異なる準位の密度を読む。

$$
\begin{aligned}
D&=\begin{pmatrix}\frac{3}{8} & 0 & 0\\0 & \frac{3}{8} & 0\\0 & 0 & \frac{1}{4}\end{pmatrix}\\
(i,j)&=\left( 1, \  2\right)\\
\varphi(ab)&=\frac{3}{8}\\
\varphi(ba)&=\frac{3}{8}
\end{aligned}
$$

$D^{-iz}bD^{iz}=(d_i/d_j)^{-iz}b$ である。 そこで下境界では $F(t)$ は $d_j$ にこの位相を掛けたものになる。上端 $i1$ では位相の指数が $1$ になり、比 $d_i/d_j$ が掛かる。

$$
\begin{aligned}
d_i/d_j&=1\\
F(0)&=\frac{3}{8}\\
F(i1)&=\frac{3}{8}
\end{aligned}
$$

時間 $t$ を加えても双方に同じ実時間の位相が掛かるので、$F(t+i1)=\varphi(\alpha_t(b)a)$ になる。有限次元の行列単位は線形基底なので、一般の $a,b$ に対しては各項を足せばよい。

下側と上側の値を比較するには、解析延長の位相だけでなく、境界で演算子を掛ける順序の変化も追う必要がある。
:::

:::example[順序を入れ替えた相関]
忠実な対角密度 $D$ の状態を取り、$\alpha_t=\sigma^\varphi_{-t}$ の上側の帯で KMS 相関を計算する。

$a=E_{ji}$、$b=E_{ij}$ とする。$ab=E_{jj}$、$ba=E_{ii}$ なので、二つの順序は異なる準位の密度を読む。

$$
\begin{aligned}
D&=\begin{pmatrix}\frac{1}{3} & 0 & 0\\0 & \frac{1}{6} & 0\\0 & 0 & \frac{1}{2}\end{pmatrix}\\
(i,j)&=\left( 3, \  2\right)\\
\varphi(ab)&=\frac{1}{6}\\
\varphi(ba)&=\frac{1}{2}
\end{aligned}
$$

$D^{-iz}bD^{iz}=(d_i/d_j)^{-iz}b$ である。 そこで下境界では $F(t)$ は $d_j$ にこの位相を掛けたものになる。上端 $i1$ では位相の指数が $1$ になり、比 $d_i/d_j$ が掛かる。

$$
\begin{aligned}
d_i/d_j&=3\\
F(0)&=\frac{1}{6}\\
F(i1)&=\frac{1}{2}
\end{aligned}
$$

時間 $t$ を加えても双方に同じ実時間の位相が掛かるので、$F(t+i1)=\varphi(\alpha_t(b)a)$ になる。有限次元の行列単位は線形基底なので、一般の $a,b$ に対しては各項を足せばよい。

下側と上側の値を比較するには、解析延長の位相だけでなく、境界で演算子を掛ける順序の変化も追う必要がある。
:::

### 中心化環と状態を保存する圧縮

:::proposition[中心化環のトレース条件]{#prop-oa-centralizer-trace}
巡回分離ベクトルの状態 $\varphi$ に対し、$\sigma_t^\varphi(x)=x$ が全 $t$ で成り立つことと、$\varphi(xy)=\varphi(yx)$ が全 $y\in M$ で成り立つことは同値である。
:::

:::proof
$\alpha_t=\sigma_{-t}^\varphi$ とする。$x$ が固定されるなら $a=y,b=x$ に対する KMS 関数の下端は一定である。下端からの一意性により上端も同じ定数なので $\varphi(yx)=\varphi(xy)$ である。逆に $x$ がトレース条件を満たせば、共役を取ることで $x^*$ も満たす。$a=x^*,b=x$ に対する KMS 関数は、上端 $\varphi(\alpha_t(x)x^*)$ と下端 $\varphi(x^*\alpha_t(x))$ が一致する。帯を虚方向に貼り合わせると有界な整関数になり、従って一定である。よって $\varphi(x^*\alpha_t(x))=\varphi(x^*x)$ である。状態の不変性を用いて $\norm{(\alpha_t(x)-x)\xi}^2$ を展開すると零になり、分離性により $\alpha_t(x)=x$ である。
:::

この結果により $\varphi$ を中心化環に制限したものは忠実正規トレースである。一般の部分環 $N$ についても、モジュラー群での不変性が状態を保つ期待値の存在を決める。

:::theorem[Takesaki の条件付き期待値定理の状態版]{#thm-oa-takesaki-expectation}
$M$ に巡回分離単位ベクトル $\xi$ があり、その状態を $\varphi$ とする。$M$ と単位を共有する von Neumann 部分環 $N$ への $\varphi$ を保存する条件付き期待値が存在することは、$\sigma_t^\varphi(N)=N$ が全 $t$ で成り立つことと同値である。存在するとき期待値は一意で、忠実かつ正規である。
:::

:::proof
$H_N=\overline{N\xi}$ への直交射影を $P$ とする。$\sigma_t^\varphi$ が $N$ を保つなら、$H_N$ は $\Delta^{it}$ で還元される。$N$ の $H_N$ 上の表現は忠実正規であり、$\xi$ は巡回分離的である。制限状態は制限された群に KMS 条件を満たすので、:ref[prop-oa-kms-uniqueness] によりそのモジュラー群は $\sigma_t^\varphi|_N$ となる。両群の実装ユニタリは $N\xi$ で一致するため、$\Delta_N^{it}=\Delta^{it}|_{H_N}$ である。これから $\Delta_N^{1/2}$ も $\Delta^{1/2}$ の制限であり、$N\xi$ 上の対合が一致することから $J_N\Delta^{1/2}n\xi=J\Delta^{1/2}n\xi$ を得る。$\Delta_N^{1/2}N\xi$ は稠密である。実際、対合は $N\xi$ をそれ自身へ写し、$J_N$ は反ユニタリである。従って $J_N=J|_{H_N}$ である。

任意の $x\in M$ の圧縮 $P x|_{H_N}$ は $J_NNJ_N$ と可換である。実際、$JNJ\subset M'$ は $x$ と可換であり、$H_N$ と直交補を保つため $P$ とも可換である。$N$ に対する Tomita の基本定理から $P x|_{H_N}$ は $N$ の像に属する。その元を $E(x)$ とすると、圧縮なので正規かつ正であり、$N$ を固定し両側加群性を満たす。また $\xi\in H_N$ なので $\varphi(E(x))=\varphi(x)$ であり、状態の忠実性から $E$ も忠実である。

逆に $\varphi$ を保存する期待値があれば、$b\in N$ に対して $\varphi(b^*E(x))=\varphi(b^*x)$ なので $E(x)\xi=Px\xi$ である。これと分離性から一意性も従う。対合を保つ $E$ は $M\xi$ 上で $P S_0=S_0 P$ を満たす。グラフを閉じると $P$ と $1-P$ は $S$ を還元するため $S^*$ も還元し、$\Delta$ と可換である。従って $\sigma_t^\varphi(n)\xi=\Delta^{it}n\xi$ は $H_N$ に入る。これに $P$ を掛けて分離性を使うと $E(\sigma_t^\varphi(n))=\sigma_t^\varphi(n)$、すなわち $\sigma_t^\varphi(n)\in N$ である。$-t$ にも同じことを適用して等号を得る。
:::

## Connes のコサイクル定理と III 型因子環の分類（概観）

:::lemma[状態と可換な射影による角]{#lem-oa-modular-corner}
忠実正規状態 $\theta$ と射影 $p\in M$ が $\theta(px)=\theta(xp)$（全 $x\in M$）を満たすとする。$\sigma_t^\theta$ は $p$ を固定し、$pMp$ への制限は $\theta|_{pMp}/\theta(p)$ のモジュラー群となる。ただし $p\ne0$ とする。
:::

:::proof
忠実正規状態 $\theta$ の GNS 表現で議論する。冒頭で示したように、GNS ベクトル $\xi$ は巡回分離的であり、この表現での像は von Neumann 環である。

任意の $x$ について $\theta(px(1-p))=\theta((1-p)xp)=0$ である。従って $\norm{xp\xi}^2=\theta(px^*xp)\le\theta(x^*x)$ なので、$R_p(x\xi)=xp\xi$ は有界に延びる。$xp\xi$ と $x(1-p)\xi$ が互いに直交するため $R_p$ は直交射影であり、左射影 $L_p$ と可換である。対合を定義した稠密な領域 $M\xi$ 上で $S L_p=R_p S$、$S R_p=L_p S$ となる。グラフの閉包を取ると両式が $S$ の定義域へ延びる。随伴の定義で $S^*R_p=L_p S^*$、$S^*L_p=R_p S^*$ もその定義域で成立する。よって $L_p,R_p$ は $\Delta=S^*S$ を還元する。これで $\sigma_t^\theta(p)=p$ である。

積射影 $L_pR_p$ の像は $\overline{pMp\,\xi}$ である。その空間は $S$ とその随伴を還元する。実際、積射影は $S$ と可換で、補射影にも同じことが言える。$M\xi$ のグラフ近似を積射影で圧縮すれば、角上の $S$ は $pMp\,\xi$ で定義した対合の閉包である。ゆえに正部分も角に制限したものであり、その虚数冪が角のモジュラー群を実装する。ベクトル $p\xi$ を規格化しても対合の動作は変わらない。
:::

:::theorem[状態変更のコサイクル]{#thm-oa-connes-cocycle}
同じ von Neumann 環上の忠実正規状態 $\psi,\varphi$ に対し、強連続ユニタリ族 $u_t=(D\psi:D\varphi)_t\in M$ があり、
$$
u_{s+t}=u_s\sigma_s^\varphi(u_t),\qquad \sigma_t^\psi(x)=u_t\sigma_t^\varphi(x)u_t^*
$$
を満たす。
:::
:::proof
$M_2(M)$ 上に $\theta([a_{ij}])=(\psi(a_{11})+\varphi(a_{22}))/2$ を取る。これは忠実正規状態である。正行列を平方根の積で書けば、対角成分が両方零のとき全ての列が零になるので忠実性が分かる。二つの対角射影と $\theta$ は :ref[lem-oa-modular-corner] の仮定を満たす。従って $\sigma_t^\theta$ はこれらを固定し、二つの対角の角ではそれぞれ $\sigma_t^\psi,\sigma_t^\varphi$ である。

$e_{12}$ の像は同じ初期射影と終射影を持つ部分等長作用素なので、右上成分だけを持ち、その成分 $u_t$ は $M$ のユニタリである。$e_{12}$ に下側の角の $a$ を右から掛けたものは、上側の角の $a$ を $e_{12}$ に左から掛けたものと等しい。この等式に $\sigma_t^\theta$ を適用すれば $u_t\sigma_t^\varphi(a)=\sigma_t^\psi(a)u_t$ となる。さらに $\sigma_{s+t}^\theta(e_{12})=\sigma_s^\theta(\sigma_t^\theta(e_{12}))$ の右上成分を取ると $u_{s+t}=u_s\sigma_s^\varphi(u_t)$ である。モジュラー群の強連続性から $u_t$ は強連続である。強連続性は忠実正規表現を替えても有界集合上で保たれる。これは $\norm{(u_t-u_s)\eta}^2$ を表す正規汎関数を使えば確認できる。

行列環で密度が $D_\psi,D_\varphi$ なら、二つを対角に置いた密度行列の共役で $e_{12}$ の像を計算でき、$u_t=D_\psi^{it}D_\varphi^{-it}$ となる。忠実正規半有限重みまで含める一般化や、標準形の相対モジュラー作用素による構成は :cite[takesaki-structure] を参照する。

三つの忠実正規状態 $\chi,\psi,\varphi$ に対しては、三つの状態の平均を $M_3(M)$ の対角に入れる。各二つの対角射影の和は :ref[lem-oa-modular-corner] の条件を満たし、その角で二状態の場合の式が使える。モジュラー自己同型は $e_{12}e_{23}=e_{13}$ を保つため、各像の右上成分を掛け合わせると
$$
(D\chi:D\psi)_t\,(D\psi:D\varphi)_t=(D\chi:D\varphi)_t
$$
となる。これは状態を二段階で変更したときの連鎖則である。同じ状態を二つ使う場合、スカラー行列角の交換ユニタリと状態が可換なので、対角成分と合わせて $e_{12}$ は固定される。従って $(D\varphi:D\varphi)_t=1$ であり、状態の順序を逆にしたコサイクルは同じ時刻で随伴になる。
:::

因子環に対する Connes の $S$ 不変量は、忠実正規半有限重みのモジュラー作用素のスペクトルの共通部分 $S(M)=\bigcap_\omega\sigma(\Delta_\omega)$ である。重みは状態と異なり正元上で無限大を取ってもよい。III 型では $S(M)=\{0,1\}$ が $\mathrm{III}_0$、$S(M)=\{0\}\cup\lambda^\Z$（$0<\lambda<1$）が $\mathrm{III}_\lambda$、$S(M)=[0,\infty)$ が $\mathrm{III}_1$ に対応する。この分類は型を区別するものであり、一般の因子環の同型を一意に定めるものではない。この対応と取り得る $S(M)$ の証明は半有限重みとモジュラー接合積の構造を用いるため :cite[takesaki-structure] を参照する。

コサイクルの式から、モジュラー自己同型を内部自己同型で割った類 $[\sigma_t^\varphi]$ は状態に依存しない。実際、状態を変更する式の二つの自己同型はユニタリ $u_t$ による内部自己同型の分だけ異なる。従って $T(M)=\{t:[\sigma_t^\varphi]=1\}$ も環だけで決まる。例えば行列因子では $\sigma_t^\varphi=\operatorname{Ad}(D^{it})$ なので $T(M)=\R$ である。忠実正規な半有限トレースを持つ因子では、そのトレースに対応するモジュラー作用素が恒等になるため $S(M)=\{1\}$ である。ここで $1$ が各モジュラー作用素のスペクトルに入ることは一般の重みに対するモジュラー理論の結果である。

さらに $\sigma^\varphi$ による [von Neumann 環の接合積](/math/operator-algebras/crossed-products/) $M\rtimes_{\sigma^\varphi}\R$ の中心への双対作用を :term[重みの流れ]{en="flow of weights" yomi="おもみのながれ"} と呼ぶ。接合積では $M$ と、時刻 $t$ での自己同型を実装するユニタリ $\lambda(t)$ を使う。双対作用は $M$ を固定し、$\lambda(t)$ に $e^{-ist}$ を掛ける。異なる状態の接合積は、$\sigma_t^\psi$ の実装ユニタリを $u_t\lambda(t)$ へ写して同型になる。コサイクルの式がユニタリ群の積を保ち、状態変更の共役式が $M$ への作用を保つためである。この同型は双対作用も保つので、中心に制限した流れも状態に依存しない。特に $\mathrm{III}_0$ 型の内部を比較するとき、$S(M)$ の値は共通となるため、重みの流れを調べる意義がある。

### コサイクルは群準同型とは限らない

行列密度が $D,E$ なら $u_t=E^{it}D^{-it}$ である。これは
$$
u_t\sigma_t^\varphi(u_s)
=E^{it}D^{-it}D^{it}E^{is}D^{-is}D^{-it}
=E^{i(t+s)}D^{-i(t+s)}=u_{t+s}
$$
を満たす。中間の $\sigma_t^\varphi$ がない積では $D^{-it}$ と $E^{is}$ の順序を交換できない。従って密度が可換なら $u_t$ はユニタリ群になるが、一般のコサイクルには上の修正項が要る。この等式は二つの状態のモジュラー群を比べる際、順序を追うために使われる。

:::example[二準位の固有軸を交換する]
二つの密度を $D=\diag(1,b)/(1+b)$、$E=RDR^*$ とする。実直交行列 $R$ は固有軸を回しているので、$D$ と $E$ を同時には対角化できない場合がある。

$t=\pi r/\log b$ として共通の位相を除く。すると $D^{it}$ の共役作用は $U=\diag(1,e^{i\pi r})$ によるもので、$E^{it}$ の位相部分は $V=RUR^*$ である。

$$
\begin{aligned}
(b,r)&=\left( 2, \  1\right)\\
R&=\begin{pmatrix}\frac{3}{5} & - \frac{4}{5}\\\frac{4}{5} & \frac{3}{5}\end{pmatrix}\\
U&=\begin{pmatrix}1 & 0\\0 & -1\end{pmatrix}\\
V&=\begin{pmatrix}- \frac{7}{25} & \frac{24}{25}\\\frac{24}{25} & \frac{7}{25}\end{pmatrix}
\end{aligned}
$$

状態変更のコサイクルでは $u_t=VU^*$ である。$t$ を二倍すると $u_{2t}=V^2(U^*)^2$ になる。一般の二乗 $(VU^*)^2$ も比較する。

$$
\begin{aligned}
u_t&=\begin{pmatrix}- \frac{7}{25} & - \frac{24}{25}\\\frac{24}{25} & - \frac{7}{25}\end{pmatrix}\\
u_{2t}&=\begin{pmatrix}1 & 0\\0 & 1\end{pmatrix}\\
(u_t)^2&=\begin{pmatrix}- \frac{527}{625} & \frac{336}{625}\\- \frac{336}{625} & - \frac{527}{625}\end{pmatrix}
\end{aligned}
$$

モジュラー群による修正を入れれば $u_t\sigma_t^\varphi(u_t)=(VU^*)U(VU^*)U^*=V^2(U^*)^2$ であり、正しい二倍の値になる。

$$
\begin{aligned}
u_t\sigma_t^\varphi(u_t)&=\begin{pmatrix}1 & 0\\0 & 1\end{pmatrix}
\end{aligned}
$$

相対位相の共通スカラーが消えるように、ここでは二つの密度に同じ固有値を用いた。固有軸の違いだけでもコサイクルの非可換な積が現れる。
:::

:::example[四分の一周期でコサイクルを求める]
二つの密度を $D=\diag(1,b)/(1+b)$、$E=RDR^*$ とする。実直交行列 $R$ は固有軸を回しているので、$D$ と $E$ を同時には対角化できない場合がある。

$t=\pi r/\log b$ として共通の位相を除く。すると $D^{it}$ の共役作用は $U=\diag(1,e^{i\pi r})$ によるもので、$E^{it}$ の位相部分は $V=RUR^*$ である。

$$
\begin{aligned}
(b,r)&=\left( 2, \  \frac{1}{2}\right)\\
R&=\begin{pmatrix}\frac{3}{5} & - \frac{4}{5}\\\frac{4}{5} & \frac{3}{5}\end{pmatrix}\\
U&=\begin{pmatrix}1 & 0\\0 & i\end{pmatrix}\\
V&=\begin{pmatrix}\frac{9}{25} + \frac{16 i}{25} & \frac{12}{25} - \frac{12 i}{25}\\\frac{12}{25} - \frac{12 i}{25} & \frac{16}{25} + \frac{9 i}{25}\end{pmatrix}
\end{aligned}
$$

状態変更のコサイクルでは $u_t=VU^*$ である。$t$ を二倍すると $u_{2t}=V^2(U^*)^2$ になる。一般の二乗 $(VU^*)^2$ も比較する。

$$
\begin{aligned}
u_t&=\begin{pmatrix}\frac{9}{25} + \frac{16 i}{25} & \frac{12 i \left(-1 + i\right)}{25}\\\frac{12}{25} - \frac{12 i}{25} & \frac{9}{25} - \frac{16 i}{25}\end{pmatrix}\\
u_{2t}&=\begin{pmatrix}- \frac{7}{25} & - \frac{24}{25}\\\frac{24}{25} & - \frac{7}{25}\end{pmatrix}\\
(u_t)^2&=\begin{pmatrix}- \frac{463}{625} + \frac{288 i}{625} & - \frac{216}{625} - \frac{216 i}{625}\\\frac{216}{625} - \frac{216 i}{625} & - \frac{463}{625} - \frac{288 i}{625}\end{pmatrix}
\end{aligned}
$$

モジュラー群による修正を入れれば $u_t\sigma_t^\varphi(u_t)=(VU^*)U(VU^*)U^*=V^2(U^*)^2$ であり、正しい二倍の値になる。

$$
\begin{aligned}
u_t\sigma_t^\varphi(u_t)&=\begin{pmatrix}- \frac{7}{25} & - \frac{24}{25}\\\frac{24}{25} & - \frac{7}{25}\end{pmatrix}
\end{aligned}
$$

相対位相の共通スカラーが消えるように、ここでは二つの密度に同じ固有値を用いた。固有軸の違いだけでもコサイクルの非可換な積が現れる。
:::

:::example[異なる密度の固有軸が傾く]
二つの密度を $D=\diag(1,b)/(1+b)$、$E=RDR^*$ とする。実直交行列 $R$ は固有軸を回しているので、$D$ と $E$ を同時には対角化できない場合がある。

$t=\pi r/\log b$ として共通の位相を除く。すると $D^{it}$ の共役作用は $U=\diag(1,e^{i\pi r})$ によるもので、$E^{it}$ の位相部分は $V=RUR^*$ である。

$$
\begin{aligned}
(b,r)&=\left( 3, \  1\right)\\
R&=\begin{pmatrix}\frac{5}{13} & - \frac{12}{13}\\\frac{12}{13} & \frac{5}{13}\end{pmatrix}\\
U&=\begin{pmatrix}1 & 0\\0 & -1\end{pmatrix}\\
V&=\begin{pmatrix}- \frac{119}{169} & \frac{120}{169}\\\frac{120}{169} & \frac{119}{169}\end{pmatrix}
\end{aligned}
$$

状態変更のコサイクルでは $u_t=VU^*$ である。$t$ を二倍すると $u_{2t}=V^2(U^*)^2$ になる。一般の二乗 $(VU^*)^2$ も比較する。

$$
\begin{aligned}
u_t&=\begin{pmatrix}- \frac{119}{169} & - \frac{120}{169}\\\frac{120}{169} & - \frac{119}{169}\end{pmatrix}\\
u_{2t}&=\begin{pmatrix}1 & 0\\0 & 1\end{pmatrix}\\
(u_t)^2&=\begin{pmatrix}- \frac{239}{28561} & \frac{28560}{28561}\\- \frac{28560}{28561} & - \frac{239}{28561}\end{pmatrix}
\end{aligned}
$$

モジュラー群による修正を入れれば $u_t\sigma_t^\varphi(u_t)=(VU^*)U(VU^*)U^*=V^2(U^*)^2$ であり、正しい二倍の値になる。

$$
\begin{aligned}
u_t\sigma_t^\varphi(u_t)&=\begin{pmatrix}1 & 0\\0 & 1\end{pmatrix}
\end{aligned}
$$

相対位相の共通スカラーが消えるように、ここでは二つの密度に同じ固有値を用いた。固有軸の違いだけでもコサイクルの非可換な積が現れる。
:::

:::example[逆順の積が異なる模型]
二つの密度を $D=\diag(1,b)/(1+b)$、$E=RDR^*$ とする。実直交行列 $R$ は固有軸を回しているので、$D$ と $E$ を同時には対角化できない場合がある。

$t=\pi r/\log b$ として共通の位相を除く。すると $D^{it}$ の共役作用は $U=\diag(1,e^{i\pi r})$ によるもので、$E^{it}$ の位相部分は $V=RUR^*$ である。

$$
\begin{aligned}
(b,r)&=\left( 3, \  \frac{1}{2}\right)\\
R&=\begin{pmatrix}\frac{3}{5} & - \frac{4}{5}\\\frac{4}{5} & \frac{3}{5}\end{pmatrix}\\
U&=\begin{pmatrix}1 & 0\\0 & i\end{pmatrix}\\
V&=\begin{pmatrix}\frac{9}{25} + \frac{16 i}{25} & \frac{12}{25} - \frac{12 i}{25}\\\frac{12}{25} - \frac{12 i}{25} & \frac{16}{25} + \frac{9 i}{25}\end{pmatrix}
\end{aligned}
$$

状態変更のコサイクルでは $u_t=VU^*$ である。$t$ を二倍すると $u_{2t}=V^2(U^*)^2$ になる。一般の二乗 $(VU^*)^2$ も比較する。

$$
\begin{aligned}
u_t&=\begin{pmatrix}\frac{9}{25} + \frac{16 i}{25} & \frac{12 i \left(-1 + i\right)}{25}\\\frac{12}{25} - \frac{12 i}{25} & \frac{9}{25} - \frac{16 i}{25}\end{pmatrix}\\
u_{2t}&=\begin{pmatrix}- \frac{7}{25} & - \frac{24}{25}\\\frac{24}{25} & - \frac{7}{25}\end{pmatrix}\\
(u_t)^2&=\begin{pmatrix}- \frac{463}{625} + \frac{288 i}{625} & - \frac{216}{625} - \frac{216 i}{625}\\\frac{216}{625} - \frac{216 i}{625} & - \frac{463}{625} - \frac{288 i}{625}\end{pmatrix}
\end{aligned}
$$

モジュラー群による修正を入れれば $u_t\sigma_t^\varphi(u_t)=(VU^*)U(VU^*)U^*=V^2(U^*)^2$ であり、正しい二倍の値になる。

$$
\begin{aligned}
u_t\sigma_t^\varphi(u_t)&=\begin{pmatrix}- \frac{7}{25} & - \frac{24}{25}\\\frac{24}{25} & - \frac{7}{25}\end{pmatrix}
\end{aligned}
$$

相対位相の共通スカラーが消えるように、ここでは二つの密度に同じ固有値を用いた。固有軸の違いだけでもコサイクルの非可換な積が現れる。

ここで $(VU^*)_{12}=(-12-12i)/25$ だが、$(U^*V)_{12}=(12-12i)/25$ である。右から $U^*$ を掛ける場合は第 2 列に $-i$ が掛かり、左から掛ける場合は第 1 行がそのまま残るためである。従って密度の順序は、共通の固有値だけを見ても決められない。
:::

## 演習問題

::::exercise{level=basic}
$D=\diag(2/5,3/5)$ の $\Delta(E_{12})$ を求めよ。
:::solution
$\Delta=L_D R_{D^{-1}}$ で左から第 1 固有値を掛け、右から第 2 固有値で割る。従って $(2/3)E_{12}$ である。
:::
::::
::::exercise{level=basic}
$H=M_2$ の標準形で $JL_{E_{12}}J$ を求めよ。
:::solution
$b\mapsto(E_{12}b^*)^*=bE_{21}$ なので $R_{E_{21}}$ である。
:::
::::
::::exercise{level=basic}
$\sigma_t^\varphi(E_{21})$ を $D=\diag(1/3,2/3)$ について求めよ。
:::solution
左右の行列積から $(d_2/d_1)^{it}E_{21}=2^{it}E_{21}$ である。
:::
::::
::::exercise{level=basic}
$D=\diag(1/3,2/3)$、$\alpha_t=\sigma_{-t}^\varphi$ とする。$a=E_{21},b=E_{12}$ について KMS 関数の $F(0),F(i)$ を求めよ。
:::solution
$a b=E_{22}$、$b a=E_{11}$ なので境界値は順に $2/3$、$1/3$ である。実際 $F(z)=\frac23\,2^{iz}$ は $z=i$ で半分になる。
:::
::::
::::exercise{level=advanced}
忠実な行列状態のモジュラー群の全固定点を記述せよ。
:::solution
$E_{ij}$ の係数には $(d_i/d_j)^{it}$ が掛かる。全 $t$ で変わらないのは $d_i=d_j$ の場合であり、固定点は固有値ごとのブロック対角環である。これは $D$ の可換子でもある。
:::
::::
::::exercise{level=advanced}
密度行列が非可逆なとき $\xi=D^{1/2}$ が左作用に分離的でないことを示せ。
:::solution
$\Ker D$ への非零射影 $p$ に対し $pD^{1/2}=0$ である。$L_p$ は非零なのに $\xi$ を消すため、分離性がない。
:::
::::
::::exercise{level=advanced}
$D_\psi,D_\varphi$ が可換な行列密度の場合、コサイクルが通常の一径数ユニタリ群であることを示せ。
:::solution
二つの対数は同じ同時対角化で可換である。そのため $u_t=e^{it(\log D_\psi-\log D_\varphi)}$ であり、指数法則で $u_{s+t}=u_su_t$ となる。
:::
::::
::::exercise{level=advanced}
$\varphi$ の状態がトレースである場合に $\Delta=1$ となることをグラフから示せ。
:::solution
$\norm{S_0x\xi}^2=\varphi(xx^*)=\varphi(x^*x)=\norm{x\xi}^2$ である。稠密な領域の反線形等長写像は全空間に延び、値域は $M\xi$ を含むため全射である。その極分解の正部分は恒等なので $\Delta=1$ である。
:::
::::

## 確認テスト

::::quiz{answer="4" tol="1e-6"}
$D=\diag(1/5,4/5)$ から作る $\Delta$ の最大固有値を答えよ。
:::explanation
固有値は $d_i/d_j$ であり、最大は $(4/5)/(1/5)=4$ である。
:::
::::
::::quiz
Tomita の基本定理の等式はどれか。

- [x] $JMJ=M'$
- [x] $\Delta^{it}M\Delta^{-it}=M$
- [ ] $\Delta=1$ が常に成り立つ。

:::explanation
$J$ は環を可換子に送り、$\Delta$ の虚数冪は環の自己同型を実装する。密度に異なる固有値がある行列状態では $\Delta$ に $1$ 以外の固有値がある。
:::
::::
::::quiz{answer="0.25" tol="1e-6"}
$D=\diag(1/4,3/4)$ の $\varphi$ で $a=E_{21},b=E_{12}$ とする。$\alpha_t=\sigma_{-t}^\varphi$ の帯の境界値 $F(i)$ を求めよ。
:::explanation
KMS 上側の境界は $\varphi(ba)$ で、$ba=E_{11}$ である。従って値は $d_1=1/4$ である。
:::
::::
::::quiz
行列標準形で $J$ が $b$ に作用した結果はどれか。

- [x] $b^*$
- [ ] $DbD^{-1}$
- [ ] $b$ が全ての $b$ で成り立つ。

:::explanation
$J$ は反線形な随伴であり $J^2=1$ である。$b\mapsto DbD^{-1}$ は正線形作用素 $\Delta$ である。例えば $b=E_{12}$ の随伴は $E_{21}$ であり元とは異なる。
:::
::::

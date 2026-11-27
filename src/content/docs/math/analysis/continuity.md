---
title: 関数の極限と連続性
description: ε–δ 条件と点列による判定から出発し、無限遠での極限と Cauchy の判定条件、不連続点の分類、中間値の定理、最大値・最小値の定理、一様連続性と連続拡張、単調関数の不連続点と逆関数の連続性を扱う。
notation:
  - symbol: '\lim_{x\to a}f(x)'
    meaning: 定義域の点が a に近づくときの f の極限
  - symbol: 'f(a-),\ f(a+)'
    meaning: a における左極限と右極限
  - symbol: 'f(D)'
    meaning: 定義域 D の像
  - symbol: 'f^{-1}'
    meaning: 単射 f の像の上で定義される逆関数
  - symbol: 'U_r(\alpha)'
    meaning: α ∈ ℝ ∪ {±∞} の r-近傍
---

## この章の目標

- 関数の極限と連続性を ε–δ 条件で定義し、点列を使った判定、極限の演算法則、Cauchy の判定条件を証明する。無限遠での極限と無限大の極限を同じ枠組みで扱う。
- 不連続点を分類し、のこぎり波・Dirichlet 関数・Thomae 関数などの例で連続性の局所的な性質を理解する。
- 区間上の連続関数について、中間値の定理と最大値・最小値の定理を証明し、方程式の解の存在に応用する。
- 一様連続性を誤差の選び方から理解し、Heine–Cantor の定理と有界区間からの連続拡張を証明する。
- 単調関数の不連続点が可算個の跳びに限ることを示し、連続な単射が狭義単調であること、逆関数が連続であることを証明する。

数列では「ある添字より後のすべての項」を調べた。関数では変数の値に、近づく順番があらかじめ付いていない。
そこで「ある距離より近いすべての点」で値の誤差を調べる。この定義を使えば、近づき方を一つの数列として選ぶ方法も使えるようになる。

## 関数の極限（ε-δ）

### 定義域から近づく

$D\subset\R$、$f:D\to\R$ とする。各 $r>0$ に対し $0<\abs{x-a}<r$ となる $x\in D$ があるとき、
$a$ を $D$ の :term[集積点]{en="accumulation point" yomi="しゅうせきてん"} という。
$a$ 自身は $D$ の元でも、そうでなくてもよい。例えば $0$ は $(0,1]$ の集積点であり、$\N$ は集積点をもたない。

:::definition[関数の極限]{#def-anl-function-limit}
$a$ を $D$ の集積点とする。実数 $L$ について

$$
\forall\eps>0\quad\exists\delta>0\quad\forall x\in D,\qquad
0<\abs{x-a}<\delta\ \Longrightarrow\ \abs{f(x)-L}<\eps
\label{eq-anl-function-limit}
$$

が成り立つとき、$x$ が $a$ に近づくときの $f(x)$ の極限は $L$ であるといい、$\lim_{x\to a}f(x)=L$ または $x\to a$ のとき $f(x)\to L$ と書く。
:::

式 :eqref[eq-anl-function-limit] では、要求された出力の誤差 $\eps$ に応じて入力の幅 $\delta$ を選ぶ。
ある $\delta$ で条件が成り立てば、それより小さい正の幅も使える。このことは複数の条件を同時に課すときに役立つ。
$0<\abs{x-a}$ としたので、$f(a)$ の値を変えてもこの極限は変化しない。また $a$ が集積点であることを仮定したのは、条件を満たす $x$ が一つもないと任意の $L$ が極限になってしまうからである。

:::example[二乗の極限]
$\lim_{x\to2}x^2=4$ を定義から証明する。$\abs{x^2-4}=\abs{x-2}\abs{x+2}$ なので、まず $\abs{x-2}<1$ と制限すれば $\abs{x+2}<5$ となる。
そこで指定された $\eps>0$ に対して

$$
\delta=\min\{1,\eps/5\}
$$

ととる。$0<\abs{x-2}<\delta$ ならば $\abs{x^2-4}<5\delta\le\eps$ となる。
例えば出力の誤差を $1/100$ 未満にするには、$\delta=1/500$ で十分である。
:::

### 片側からの極限

$x$ を $D\cap(-\infty,a)$ に限った極限を :term[左極限]{en="left-hand limit" yomi="ひだりきょくげん"} $f(a-)$、
$D\cap(a,\infty)$ に限った極限を :term[右極限]{en="right-hand limit" yomi="みぎきょくげん"} $f(a+)$ という。
いずれも、その側で $a$ が集積点である場合に用いる記号であり、$\lim_{x\to a-0}f(x)$、$\lim_{x\uparrow a}f(x)$ などとも書く。

両側から近づけるなら、極限が $L$ であることと両方の片側極限が $L$ であることは同値である。
実際、両側の条件から片側への制限は直ちに従う。逆に同じ $\eps$ に対して左と右で得た幅の小さい方をとれば、どちら側の $x$ でも誤差が $\eps$ 未満となる。

:::example[片側で異なる値に近づく関数]
$f(x)=x/\abs{x}$ を $D=\R\setminus\{0\}$ で定義する。$x<0$ なら値は $-1$、$x>0$ なら値は $1$ である。
従って $f(0-)=-1,\ f(0+)=1$ である。両方に共通する極限はない。
:::

### 点列で極限を判定する

:::theorem[極限の点列による特徴付け]{#thm-anl-sequential-limit}
$a$ を $D$ の集積点とする。$\lim_{x\to a}f(x)=L$ であるための必要十分条件は、
$D\setminus\{a\}$ の点からなる任意の列 $x_n\to a$ について $f(x_n)\to L$ となることである。
:::

:::proof
極限が $L$ とする。任意の $\eps>0$ に対して :ref[def-anl-function-limit] の $\delta$ をとる。
$x_n\to a$ なので $n\ge N$ なら $\abs{x_n-a}<\delta$ となる $N$ がある。全項で $x_n\ne a$ だから、その範囲で $\abs{f(x_n)-L}<\eps$ である。

逆に関数の極限の条件が破れると仮定する。否定を量化記号に沿って書くと、ある $\eps_0>0$ に対し、**どの** $\delta>0$ にも

$$
0<\abs{x-a}<\delta,\qquad \abs{f(x)-L}\ge\eps_0
$$

を満たす $x\in D$ がある。各 $n$ で $\delta=1/n$ としてそのような点 $x_n$ を一つ選ぶ。
すると $0<\abs{x_n-a}<1/n$ より $x_n\to a$ であるが、$f(x_n)$ は $L$ に収束しない。これは点列に関する仮定に反する。
:::

この定理から極限の一意性も従う。集積点の定義により $0<\abs{x_n-a}<1/n$ となる $x_n\in D$ を選べる。
二つの値がともに極限の定義を満たせば、同じ数列 $f(x_n)$ が両方の値に収束するので、:ref[prop-anl-limit-basics] によって両者は等しい。

:::proposition[関数の極限の演算と順序]{#prop-anl-function-limit-laws}
$a$ を $D$ の集積点とし、$f,g,h:D\to\R$ について $\lim_{x\to a}f(x)=L$、$\lim_{x\to a}g(x)=M$ とする。

1. $f\pm g\to L\pm M$、$fg\to LM$ である。
2. $L>0$ ならば、ある $\delta>0$ があり、$x\in D$、$0<\abs{x-a}<\delta$ で $f(x)>L/2$ である。
3. $M\ne0$ ならば、$a$ の近くで $g(x)\ne0$ であり、$f/g\to L/M$ である。
4. ある $\delta_0>0$ について $x\in D$、$0<\abs{x-a}<\delta_0$ で $f(x)\le g(x)$ ならば $L\le M$ である。
5. 同じ範囲で $f(x)\le h(x)\le g(x)$ かつ $L=M$ ならば $h(x)\to L$ である。
:::

:::proof
2 は定義で $\eps=L/2$ とすればよい。$\abs g\to\abs M$（三角不等式による）に 2 を適用すると、$a$ の近くで $\abs{g(x)}>\abs M/2>0$ となり、3 の前半を得る。
その範囲の点からなる集合 $D'$ についても $a$ は集積点である。

残りは点列で示す。$D\setminus\{a\}$（3 では $D'\setminus\{a\}$）の点列 $x_n\to a$ を任意にとると、:ref[thm-anl-sequential-limit] により $f(x_n)\to L$、$g(x_n)\to M$ である。
:ref[thm-anl-limit-laws] により $f(x_n)\pm g(x_n)\to L\pm M$ などが成り立ち、再び :ref[thm-anl-sequential-limit] を用いて 1 と 3 を得る。
4 と 5 では、十分後の $n$ で $\abs{x_n-a}<\delta_0$ となるので :ref[prop-anl-squeeze] を適用すればよい。
:::

:::example[穴の開いた有理関数]
$x\ne1$ で定まる $h(x)=(x^2-1)/(x-1)$ は、定義域のどの点でも $x+1$ に等しい。
従って $\abs{h(x)-2}=\abs{x-1}$ であり、$\delta=\eps$ と選べば $\lim_{x\to1}h(x)=2$ を得る。
この計算は $h(1)$ の値を使っていない。次節では、点 $1$ を定義域に加える方法を考える。
:::

### 無限遠での極限と無限大の極限

$x\to+\infty$ のときの極限や、値が $+\infty$ に発散する場合も、近傍の言葉を使えば同じ形で定義できる。

:::definition[拡張実数での近傍と極限]{#def-anl-extended-limit}
$r>0$ に対し、$\alpha\in\R$ なら $U_r(\alpha)=(\alpha-r,\alpha+r)$、$U_r(+\infty)=(1/r,+\infty)$、$U_r(-\infty)=(-\infty,-1/r)$ とおき、$\alpha$ の $r$-近傍という。

1. 実数列 $(x_n)$ が $\alpha\in\overline{\R}$ に向かうとは、任意の $r>0$ に対し、十分後のすべての $n$ で $x_n\in U_r(\alpha)$ となることである。
2. $\alpha\in\overline{\R}$ が $D$ の集積点であるとは、任意の $r>0$ について $U_r(\alpha)$ が $D\setminus\{\alpha\}$ の点を含むことである。$+\infty$ が集積点であることは $D$ が上に有界でないことと同値である。
3. $\alpha$ を $D$ の集積点、$\beta\in\overline{\R}$ とする。任意の $\eps>0$ に対してある $r>0$ があり、$x\in(D\setminus\{\alpha\})\cap U_r(\alpha)$ ならば $f(x)\in U_\eps(\beta)$ となるとき、$\lim_{x\to\alpha}f(x)=\beta$ と書く。
:::

$\alpha,\beta\in\R$ なら 3 は :ref[def-anl-function-limit] と一致し、1 は数列の収束（:ref[def-anl-sequence-limit]）、$\alpha=\pm\infty$ なら :ref[def-anl-diverge-infinity] と一致する（$M\le0$ の場合は $M=1$ の条件から従う）。
例えば $\lim_{x\to+\infty}f(x)=L\in\R$ は「任意の $\eps>0$ に対してある $R$ があり、$x\in D$、$x>R$ ならば $\abs{f(x)-L}<\eps$」を意味し、$\lim_{x\to a}f(x)=+\infty$ は「任意の $M$ に対してある $\delta>0$ があり、$x\in D$、$0<\abs{x-a}<\delta$ ならば $f(x)>M$」を意味する。
$r\le r'$ なら $U_r(\alpha)\subset U_{r'}(\alpha)$ であることに注意する。

:::proposition[拡張された点列判定]{#prop-anl-sequential-extended}
$\alpha\in\overline{\R}$ を $D$ の集積点、$\beta\in\overline{\R}$ とする。$\lim_{x\to\alpha}f(x)=\beta$ であるための必要十分条件は、$D\setminus\{\alpha\}$ の点からなり $\alpha$ に向かう任意の列 $(x_n)$ について、$f(x_n)$ が $\beta$ に向かうことである。
:::

:::proof
極限が $\beta$ とし、$\eps>0$ に対する $r$ をとる。$x_n\to\alpha$ なので十分後で $x_n\in U_r(\alpha)$ となり、そこで $f(x_n)\in U_\eps(\beta)$ である。

逆に極限が $\beta$ でないとすると、ある $\eps_0>0$ について、どの $r>0$ にも $f(x)\notin U_{\eps_0}(\beta)$ となる $x\in(D\setminus\{\alpha\})\cap U_r(\alpha)$ がある。
$r=1/n$ としてそのような $x_n$ を選ぶ。任意の $r>0$ について $1/n\le r$ となる $n$ 以降で $x_n\in U_{1/n}(\alpha)\subset U_r(\alpha)$ なので $x_n\to\alpha$ であるが、$f(x_n)$ はどの $n$ でも $U_{\eps_0}(\beta)$ に入らず、$\beta$ に向かわない。
:::

この命題により、数列についての :ref[thm-anl-limit-laws] と :ref[prop-anl-infinite-limits] は、そのまま関数の極限に移る。

:::example[無限遠での極限と無限大の極限]
(1) $f(x)=(2x^2+1)/(x^2-3)$（$x>2$）について $\lim_{x\to+\infty}f(x)=2$ である。$x_n\to+\infty$ なら $1/x_n^2\to0$（:ref[prop-anl-infinite-limits] の 4 と積の法則）であり

$$
f(x_n)=\frac{2+1/x_n^2}{1-3/x_n^2}\to\frac{2+0}{1-0}=2
$$

となるので、:ref[prop-anl-sequential-extended] から従う。

(2) $\lim_{x\to0}1/x^2=+\infty$ である。$M>0$ に対し $\delta=1/\sqrt M$ とすれば、$0<\abs x<\delta$ で $1/x^2>M$ となる。

(3) $1/x$ は $x\to0$ で極限をもたない。右からは $+\infty$、左からは $-\infty$ に向かうからである。
:::

### Cauchy の判定条件

極限値の候補がわからないときにも、極限の存在を判定できる。

:::theorem[関数の極限に関する Cauchy の判定条件]{#thm-anl-function-cauchy}
$\alpha\in\overline{\R}$ を $D$ の集積点とする。有限な極限 $\lim_{x\to\alpha}f(x)$ が存在するための必要十分条件は、任意の $\eps>0$ に対してある $r>0$ があり

$$
x,y\in(D\setminus\{\alpha\})\cap U_r(\alpha)\quad\Longrightarrow\quad\abs{f(x)-f(y)}<\eps
$$

となることである。
:::

:::proof
極限 $L$ があれば、誤差 $\eps/2$ に対する $r$ をとると $\abs{f(x)-f(y)}\le\abs{f(x)-L}+\abs{L-f(y)}<\eps$ である。

逆に条件を仮定する。$D\setminus\{\alpha\}$ の点列 $x_n\to\alpha$ をとる。$\eps>0$ に対する $r$ をとれば、十分後の $m,n$ で $x_m,x_n\in U_r(\alpha)$ なので $\abs{f(x_m)-f(x_n)}<\eps$ である。
従って $f(x_n)$ は Cauchy 列であり、:ref[thm-anl-cauchy-complete] により収束する。

極限が点列によらないことを示す。$\alpha$ に向かう二つの列 $(x_n),(y_n)$ をとり、交互に並べた列 $z=(x_1,y_1,x_2,y_2,\dots)$ を作る。$z$ も $\alpha$ に向かうので $f(z_n)$ は収束し、その部分列 $f(x_n)$ と $f(y_n)$ は同じ極限をもつ。
この共通の極限を $L$ とすれば、:ref[prop-anl-sequential-extended] により $\lim_{x\to\alpha}f(x)=L$ である。
:::

## 連続関数

### 一点での連続性

:::definition[連続性]{#def-anl-continuity}
$a\in D$ とする。各 $\eps>0$ に対してある $\delta>0$ が存在し、

$$
x\in D,\quad\abs{x-a}<\delta
\quad\Longrightarrow\quad \abs{f(x)-f(a)}<\eps
$$

となるとき、$f$ は $a$ で :term[連続]{en="continuous" yomi="れんぞく"} であるという。
$D$ のすべての点で連続な $f$ を $D$ 上の連続関数という。
$D\cap[a,\infty)$ に制限した関数が $a$ で連続なとき右連続、$D\cap(-\infty,a]$ に制限した関数が連続なとき左連続という。
:::

$a\in D$ が集積点なら、連続性とは極限が存在して $f(a)$ と一致することである。
一方、$a$ 以外の $D$ の点がある幅の中に一つもないなら、その幅を $\delta$ に選ぶだけで連続性は成り立つ。
このような $a$ を孤立点という。閉区間 $[a,b]$ の端では区間の内側だけを調べればよい。
両側に定義域の点がある場合、連続であることと、左連続かつ右連続であることは同値である。

:::corollary[連続性の点列判定]{#cor-anl-sequential-continuity}
$f:D\to\R$ が $a\in D$ で連続であることと、$D$ 内の任意の $x_n\to a$ について $f(x_n)\to f(a)$ となることは同値である。
:::

:::proof
連続ならば $\eps$ に対応する $\delta$ を選び、列が $a$ の $\delta$ 以内に入った後を見ればよい。
逆に連続でないなら、ある $\eps_0>0$ と、各 $n$ に対する $x_n\in D$ があり、
$\abs{x_n-a}<1/n$ かつ $\abs{f(x_n)-f(a)}\ge\eps_0$ である。この列は点列の条件に反する。
:::

従って、先の $h(x)$ に $h(1)=2$ と補うと $1$ でも連続になる。別の値を補うと、$x\ne1$ からの極限 $2$ と値が一致しない。

:::proposition[連続関数の演算と合成]{#prop-anl-continuous-algebra}
$f,g:D\to\R$ が $a\in D$ で連続なら、和、差、積も $a$ で連続である。
$g(a)\ne0$ なら商は $g\ne0$ である定義域で $a$ において連続である。
さらに $f(D)\subset E$、$h:E\to\R$ が $f(a)$ で連続なら、合成 $h\circ f$ は $a$ で連続である。
:::

:::proof
任意の $x_n\to a$ に対し、:ref[cor-anl-sequential-continuity] より $f(x_n)\to f(a),\ g(x_n)\to g(a)$ である。
:ref[thm-anl-limit-laws] を適用し、再び同じ点列判定を用いれば和、差、積、商の結論を得る。
なお $g(a)\ne0$ なら、$a$ の十分近くで $\abs{g(x)-g(a)}<\abs{g(a)}/2$ となり、そこでは分母は消えない。

合成では幅を選ぶ順番を明示できる。まず $\eps>0$ に対し、$\abs{y-f(a)}<\eta$ なら $\abs{h(y)-h(f(a))}<\eps$ となる $\eta>0$ を選ぶ。
次に $f$ の連続性により $\abs{x-a}<\delta$ なら $\abs{f(x)-f(a)}<\eta$ となる $\delta>0$ を選ぶ。
このとき $y=f(x)$ は最初の条件を満たすので、結論を得る。
:::

定数関数は出力の差が $0$ なので連続であり、恒等関数 $x\mapsto x$ は $\delta=\eps$ ととれるので連続である。
従って上の命題を有限回用いれば、多項式は実数全体で、有理関数は分母が消えない各点で連続である。
$\bigl|\abs x-\abs a\bigr|\le\abs{x-a}$ より、絶対値 $x\mapsto\abs x$ も連続である。
連続関数 $f$ が $f(a)>0$ を満たせば、:ref[prop-anl-function-limit-laws] の 2 と同じ議論により $a$ の近くで $f(x)>f(a)/2>0$ となる。この「符号の保存」は以下で繰り返し使う。

### 合成関数の極限

:::proposition[合成関数の極限]{#prop-anl-composite-limit}
$f:D\to E$、$g:E\to\R$ とし、$\alpha\in\overline{\R}$ を $D$ の集積点、$\lim_{x\to\alpha}f(x)=b\in\R$ とする。次のいずれかの場合に $\lim_{x\to\alpha}g(f(x))=c$ が成り立つ。

1. $b\in E$ で $g$ が $b$ において連続であり、$c=g(b)$ である。
2. $\lim_{y\to b}g(y)=c$ であり、ある $r_0>0$ について $x\in(D\setminus\{\alpha\})\cap U_{r_0}(\alpha)$ ならば $f(x)\ne b$ である。
:::

:::proof
$D\setminus\{\alpha\}$ の点列 $(x_n)$ で $\alpha$ に向かうものをとると、:ref[prop-anl-sequential-extended] により $f(x_n)\to b$ である。
1 の場合は :ref[cor-anl-sequential-continuity] により $g(f(x_n))\to g(b)$ である。
2 の場合は、十分後の $n$ で $x_n\in U_{r_0}(\alpha)$ なので $f(x_n)\in E\setminus\{b\}$ であり、この部分の列に :ref[thm-anl-sequential-limit] を $g$ について適用して $g(f(x_n))\to c$ を得る。
いずれの場合も、再び :ref[prop-anl-sequential-extended] により結論が従う。
:::

条件 2 の $f(x)\ne b$ は外せない。$f\equiv0$、$g(y)=1$（$y\ne0$）、$g(0)=0$ とすると $\lim_{y\to0}g(y)=1$ であるが、$g(f(x))=0$ はすべての $x$ で $0$ である。

### 不連続点の分類

:::definition[不連続点の分類]{#def-anl-discontinuity}
$I$ を区間、$c$ を $I$ の内点（$c$ の両側に $I$ の点がある）、$f:I\to\R$ が $c$ で不連続であるとする。

1. $f(c-)$ と $f(c+)$ がともに（有限な値として）存在するとき、$c$ を第一種の不連続点という。
特に $f(c-)=f(c+)$ のとき :term[除去可能]{en="removable" yomi="じょきょかのう"} であるといい、$f(c-)\ne f(c+)$ のとき $c$ で :term[跳び]{en="jump" yomi="とび"} があるという。
2. 第一種でない不連続点を第二種の不連続点という。
:::

除去可能な不連続点では、$f(c)$ を共通の極限値に定め直すと連続になる。跳びのある点では、値をどう定め直しても連続にはならない。

:::example[のこぎり波と振動する関数]{#exm-anl-sawtooth}
$\varphi(t)=\min\{t-\lfloor t\rfloor,\ \lfloor t\rfloor+1-t\}$ は $t$ から最も近い整数までの距離であり、$0\le\varphi(t)\le1/2$、$\varphi(t+1)=\varphi(t)$ を満たす。
実際、整数 $k\le\lfloor t\rfloor-1$ なら $t-k>t-\lfloor t\rfloor$、$k\ge\lfloor t\rfloor+2$ なら $k-t>\lfloor t\rfloor+1-t$ なので、$\varphi(t)=\min_{k\in\Z}\abs{t-k}$ である。
任意の $s,t$ と整数 $k$ について $\varphi(s)\le\abs{s-k}\le\abs{s-t}+\abs{t-k}$ であり、$k$ を $t$ に最も近い整数にとれば $\varphi(s)\le\abs{s-t}+\varphi(t)$ となる。$s,t$ を入れ替えて

$$
\abs{\varphi(s)-\varphi(t)}\le\abs{s-t}
$$

を得る。従って $\varphi$ は連続である。

(1) $g(x)=\varphi(1/x)$（$x\ne0$）は $x\to0$ で極限をもたない。$x_n=1/n$ では $g(x_n)=\varphi(n)=0$、$y_n=1/(n+1/2)$ では $g(y_n)=\varphi(n+1/2)=1/2$ だからである。
右からも左からも同じことが起こるので、$g(0)$ をどう定めても $0$ は第二種の不連続点である。

(2) $f(x)=x\varphi(1/x)$（$x\ne0$）、$f(0)=0$ とすると $\abs{f(x)-f(0)}\le\abs x/2$ なので、$\delta=2\eps$ とすれば $f$ は $0$ で連続である。
$1/x$ が $0$ の近くで激しく変化しても、係数 $x$ が関数値の振れ幅を抑える。
[初等関数の厳密な構成](/math/analysis/elementary-functions/) で正弦関数を定義した後は、$\sin(1/x)$ と $x\sin(1/x)$ が同じ振る舞いをすることが同じ議論でわかる。
:::

### 稠密な二つの集合で値を変える

:::example[Dirichlet 関数]
有理数で $1$、無理数で $0$ をとる関数を Dirichlet 関数という。
:ref[thm-anl-rational-density] と :ref[cor-anl-irrational-density] により、各 $n$ について $(a,a+1/n)$ に有理数 $r_n$ と無理数 $s_n$ を一つずつ選べる。
どちらの列も $a$ に近づくが、関数値の列はそれぞれ定数 $1$ と定数 $0$ である。
従ってどの $a\in\R$ でも極限が存在せず、片側極限も存在しない。すべての点が第二種の不連続点である。
:::

:::example[Thomae 関数]{#exm-anl-thomae}
有理数 $x$ に対し、$nx\in\Z$ となる最小の自然数 $n$ を $q(x)$ と書く（:ref[thm-anl-well-ordering] により存在する）。これは $x$ を既約分数 $p/q$（$q\in\N$）で表したときの分母である。実数全体で

$$
T(x)=
\begin{cases}
1/q(x) & (x\in\Q),\\
0 & (x\notin\Q)
\end{cases}
$$

と定める。この Thomae 関数は、無理数のすべてで連続であり、有理数のすべてで不連続である。

実際、点 $a$ と誤差 $\eps>0$ に対して、$Q=\lfloor1/\eps\rfloor+1$ と選ぶ。
$[a-1,a+1]$ 内の有理数で $q(x)\le Q$ となるものは有限個である。$q(x)=q$ を固定すると $x=p/q$、$p=qx$ は有界区間 $[q(a-1),q(a+1)]$ 内の整数に限られ、そのような整数は有限個だからである。
その有限集合から $a$ 自身を除いたものを $F$ とし、

$$
\delta=\min\bigl(\{1\}\cup\{\abs{r-a}\mid r\in F\}\bigr)>0
$$

とおく。$0<\abs{x-a}<\delta$ ならば $x\notin F$ であるから、$x$ が有理数の場合には $q(x)>Q$、すなわち $T(x)<1/Q<\eps$ である。無理数なら値は $0$ であるから、どちらも $0\le T(x)<\eps$ を満たす。

従ってすべての点で $\lim_{x\to a}T(x)=0$ である。この値が $T(a)$ と一致するのは、ちょうど $a$ が無理数のときである。
有理数 $a$ では $T(a)=1/q(a)>0$ と極限 $0$ が異なるので、除去可能な不連続点である。
:::

下図では、分母が $6$ 以下の点を抜き出した。高さが $1/4$ 以上の点をすべて避けるには、分母が $4$ 以下の有限個の点だけを調べればよい。
証明で有限集合 $F$ と点 $a$ の距離をとったのは、この操作を任意の高さで行うためである。

```tikz caption="区間 [0,1] における Thomae 関数の値のうち、既約分母が 6 以下の点。破線の高さは 1/4 である。"
\begin{tikzpicture}[x=6cm,y=3cm]
\draw[->] (-0.02,0) -- (1.09,0) node[right] {$x$};
\draw[->] (0,-0.02) -- (0,1.12) node[above] {$T(x)$};
\draw[dashed] (0,0.25) -- (1,0.25);
\foreach \x/\v in {0/0,0.5/{1/2},1/1}
  \node[below] at (\x,0) {$\v$};
\foreach \y/\v in {0.25/{1/4},0.5/{1/2},1/1}
  \node[left] at (0,\y) {$\v$};
\foreach \p/\q in {0/1,1/1,1/2,1/3,2/3,1/4,3/4,1/5,2/5,3/5,4/5,1/6,5/6}
  \fill ({\p/\q},{1/\q}) circle (1.6pt);
\end{tikzpicture}
```

## 中間値の定理

### 区間の特徴付け

[実数の連続性](/math/analysis/real-numbers/) では区間を端点で表した九種類の集合として定めた。連続関数の像を調べるには、次の特徴付けが便利である。

:::definition[区間]{#def-anl-interval}
$I\subset\R$ が、任意の $x,y\in I$ と $x<z<y$ を満たす任意の $z$ について $z\in I$ を満たすとき、$I$ は凸であるという。
:::

:::proposition[区間の特徴付け]{#prop-anl-interval-characterization}
空でない $I\subset\R$ が凸であるための必要十分条件は、$I$ が $\R$ の区間（一点からなる $[a,a]$ を含む）であることである。
:::

:::proof
区間が凸であることは定義から明らかである。逆に $I$ を空でない凸集合とし、拡張実数の意味で $a=\inf I$、$b=\sup I$ とおく（:ref[prop-anl-extended-sup]）。$a\le b$ である。
$a<z<b$ とすると、$z$ は下界でないので $x<z$ となる $x\in I$ があり、上界でないので $z<y$ となる $y\in I$ がある。凸性から $z\in I$ である。
従って $(a,b)\subset I\subset[a,b]$（$a,b$ が $\pm\infty$ のときは $[a,b]$ の該当する端点を除く）であり、$I$ は $(a,b)$ に端点 $a,b$ のうち $I$ に属するものを付け加えた集合である。
これは区間の九種類の形のいずれかである（$a=b$ なら $I=\{a\}$）。
:::

以後、区間といえば凸集合の意味で用いる。端点を除いた $(a,b)$ の点を $I$ の内点という。

### 値を飛び越えないこと

:::theorem[中間値の定理]{#thm-anl-intermediate-value}
$a<b$、$f:[a,b]\to\R$ は連続とする。$y$ が $f(a)$ と $f(b)$ の間にあれば、$f(c)=y$ を満たす $c\in[a,b]$ が存在する。
$y$ が両端の値の間に厳密にあれば、$c$ は $(a,b)$ にとれる。
:::

:::proof
端点で等号が成り立てばその端点をとる。$f(a)>f(b)$ の場合は $-f,-y$ を考えればよいので、以下 $f(a)<y<f(b)$ とする。
集合 $A=\{x\in[a,b]\mid f(x)\le y\}$ は $a$ を含み、$b$ を上界とする。従って $c=\sup A$ があり $a\le c\le b$ である。

各 $n$ で $c-1/n<x_n\le c$ となる $x_n\in A$ を選ぶと $x_n\to c$ である。
連続性により $f(x_n)\to f(c)$ なので、$f(x_n)\le y$ と :ref[prop-anl-squeeze] から $f(c)\le y$ となる。特に $f(b)>y$ なので $c<b$ である。

もし $f(c)<y$ なら、$\eps=y-f(c)>0$ として連続性から幅 $\delta>0$ が得られる。
そこで $x=c+\tfrac12\min\{\delta,b-c\}$ とすれば $c<x<b$ かつ $f(x)<f(c)+\eps=y$ である。
つまり $x\in A$ なのに $c$ より大きく、上界の定義に反する。よって $f(c)=y$ であり、$f(a)<y$ より $c>a$ でもある。
:::

この証明では上限の公理が本質的である。$\Q$ 上の連続関数 $x^2-2$ は $0$ で負、$2$ で正であるが、$\Q$ の中に零点をもたない。

:::corollary[区間の連続像]{#cor-anl-continuous-image-interval}
区間 $I$ 上の連続関数 $f$ の像 $f(I)$ は区間である。
:::

:::proof
:ref[prop-anl-interval-characterization] により、$f(I)$ が凸であることを示せばよい。$u=f(x)$、$v=f(y)$ を $f(I)$ の二点とし、$u<w<v$ とする。
$x,y$ を端点とする閉区間は $I$ に含まれ（凸性）、その上で $f$ は連続であるから、:ref[thm-anl-intermediate-value] により $f(c)=w$ となる $c$ がその閉区間にある。従って $w\in f(I)$ である。
:::

### 零点を挟んで探す

:::example[三次方程式の零点]
$f(x)=x^3+x-1$ は連続で、$f(0)=-1,\ f(1)=1$ である。従って $(0,1)$ のどこかに零点がある。
さらに、$[0,1]$ 上で $x^3$ は増加し $x$ は狭義増加するので、$f$ は狭義増加する。このため零点は一つである。

区間を二分すると $f(1/2)=-3/8$ なので零点は $(1/2,1)$ にある。次に $f(3/4)=11/64$ なので $(1/2,3/4)$ にある。
さらに $f(5/8)=-67/512$ なので $(5/8,3/4)$ にある。
この操作を $n$ 回行えば、途中で零点そのものが得られるか、幅 $2^{-n}$ の区間が得られる。その中点を近似値にすれば誤差は $2^{-n-1}$ 以下である。
:::

ここで端点の符号は零点の**存在**を示し、狭義増加性は**一意性**を示した。
例えば連続関数 $x\mapsto x(x-1)(x+1)$ は複数の零点をもつので、中間値の定理だけから一意性は従わない。

:::proposition[奇数次多項式の実根]{#prop-anl-odd-degree-root}
$n$ が奇数のとき、実係数の多項式 $p(x)=x^n+c_{n-1}x^{n-1}+\cdots+c_1x+c_0$ は実数の零点をもつ。
:::

:::proof
$R=1+\sum_{k=0}^{n-1}\abs{c_k}$ とおく。$\abs x\ge R$ なら $\abs x\ge1$ より $\abs x^k\le\abs x^{n-1}$（$k\le n-1$）なので

$$
\left|\sum_{k=0}^{n-1}c_kx^k\right|\le\sum_{k=0}^{n-1}\abs{c_k}\abs x^{n-1}<R\abs x^{n-1}\le\abs x^n
$$

である。従って $p(R)\ge R^n-\bigl|\sum c_kR^k\bigr|>0$ である。$n$ は奇数なので $(-R)^n=-R^n$ であり、同様に $p(-R)<0$ となる。
:ref[thm-anl-intermediate-value] により $(-R,R)$ に零点がある。
:::

最高次の係数が $1$ でない場合も、その係数で割れば同じ結論を得る。偶数次では $x^2+1$ のように実根をもたないことがある。
また $t>0$ の $n$ 乗根の存在（:ref[thm-anl-nth-root]）も、$x^n-t$ に $[0,1+t]$ で中間値の定理を適用すれば別証明が得られる。

## 最大値・最小値の定理

### 有界になることと値をとること

:::theorem[最大値・最小値の定理]{#thm-anl-extreme-value}
$a\le b$ とし、$f:[a,b]\to\R$ は連続とする。ある $p,q\in[a,b]$ が存在して、すべての $x\in[a,b]$ について

$$
f(p)\le f(x)\le f(q)
$$

が成り立つ。
:::

:::proof
まず $f$ は上に有界である。もし上に有界でなければ、各 $n$ について $f(x_n)>n$ を満たす $x_n\in[a,b]$ を選べる。
有界列 $(x_n)$ には :ref[thm-anl-bolzano-weierstrass] によって収束部分列 $x_{n_k}\to c$ がある。
$a\le x_{n_k}\le b$ を極限に移すと $a\le c\le b$ なので、連続性が使える。
すると $f(x_{n_k})\to f(c)$ であるが、$f(x_{n_k})>n_k\ge k$ は有界でなく、:ref[prop-anl-limit-basics] に反する。

次に $M=\sup f([a,b])$ とおく。$M-1/n<f(y_n)\le M$ を満たす $y_n\in[a,b]$ を選ぶと、再び収束部分列 $y_{n_k}\to q\in[a,b]$ がある。
挟み撃ちから $f(y_{n_k})\to M$、連続性から $f(y_{n_k})\to f(q)$ である。
極限の一意性より $f(q)=M$ となる。連続関数 $-f$ の最大値について同じ結論を適用すると、$f$ の最小値をとる点 $p$ も得られる。
:::

この証明で定義域の有界性は部分列を選ぶために、閉区間の端点を含むことは極限 $q$ を定義域内に保つために用いた。
最後に連続性を使って、点の極限での関数値を決定した。
「有界閉集合の点列は、その集合内に極限をもつ部分列を含む」という性質は点列コンパクト性と呼ばれ、[ℝⁿ の位相と多変数関数の連続性](/math/analysis/euclidean-space/) で一般の次元に拡張される。

:::example[仮定を一つずつ外す]
$(0,1)$ 上の $f(x)=x$ は有界で連続だが、上限 $1$ も下限 $0$ も値としてはとらない。
定義域を $[0,\infty)$ にすれば同じ関数は最小値 $0$ をとるが、今度は上に有界でない。
$(0,1]$ 上の $1/x$ は連続であるが上に有界でない。

閉区間 $[0,1]$ でも、$g(x)=x$（$x<1$）、$g(1)=0$ と定めれば、像は $[0,1)$ となり最大値がない。
実際 $x_n=1-1/n$ では $g(x_n)\to1\ne g(1)$ となる。ここでは点列の極限を関数値の極限に移す段階で連続性が失われている。
:::

:::corollary[閉区間の連続像]{#cor-anl-interval-image}
空でない有界閉区間 $I$ と連続関数 $f:I\to\R$ について、$f(I)=[m,M]$ である。ここで $m$ は最小値、$M$ は最大値である。
:::

:::proof
:ref[thm-anl-extreme-value] により $m,M$ が存在して $f(I)\subset[m,M]$ となり、両端の値は達せられる。
それらの値を与える二点の間の閉区間に :ref[thm-anl-intermediate-value] を適用すれば、すべての中間の値も像に入る。
$m=M$ ならば像はその一点である。
:::

特に $[a,b]$ 上の連続関数 $f$ が各点で $f(x)>0$ ならば、最小値 $m=f(p)>0$ をもつので $f(x)\ge m>0$ が一様に成り立ち、$1/f$ は有界である。
各点での正値性が一様な下界に強まることは、閉区間上の連続関数の典型的な性質である。

## 一様連続性

### 定義域全体で一つの幅を使う

:::definition[一様連続性]{#def-anl-uniform-continuity}
$f:D\to\R$ が $D$ 上で :term[一様連続]{en="uniformly continuous" yomi="いちようれんぞく"} であるとは、

$$
\forall\eps>0\quad\exists\delta>0\quad\forall x,y\in D,\qquad
\abs{x-y}<\delta\ \Longrightarrow\ \abs{f(x)-f(y)}<\eps
\label{eq-anl-uniform}
$$

が成り立つことである。
:::

各点での連続性では基準点 $a$ を固定した後で $\delta$ を選ぶ。
式 :eqref[eq-anl-uniform] では、$x,y$ のどちらを決めるよりも前に $\delta$ を選ぶ。
従って一様連続な関数は、$y=a$ を固定すれば各点でも連続である。

:::definition[Lipschitz 連続性]{#def-anl-lipschitz}
定数 $K\ge0$ があり、すべての $x,y\in D$ について $\abs{f(x)-f(y)}\le K\abs{x-y}$ となるとき、$f$ は :term[Lipschitz 連続]{en="Lipschitz continuous" yomi="りぷしっつれんぞく"} であるという。
:::

Lipschitz 連続な関数は一様連続である。$K>0$ として $\delta=\eps/K$ とすればよい（$K=0$ なら定数関数である）。
例えばのこぎり波 $\varphi$（:ref[exm-anl-sawtooth]）や $x\mapsto\abs x$ は $K=1$ で Lipschitz 連続である。

:::example[区間によって幅が変わる二乗関数]
$R>0$ を固定する。$x,y\in[-R,R]$ ならば

$$
\abs{x^2-y^2}=\abs{x-y}\abs{x+y}\le2R\abs{x-y}
$$

であるから、$\delta=\eps/(2R)$ が区間全体で使える。
一方、$D=\R$ では $x_n=n,\ y_n=n+1/n$ とおくと $\abs{x_n-y_n}=1/n\to0$ にもかかわらず

$$
\abs{y_n^2-x_n^2}=2+1/n^2>2.
$$

従って $\eps=1$ に使える共通の幅はない。任意の $\delta>0$ に対して $n>1/\delta$ を選べば、その幅に収まる反例になるからである。
:::

:::example[Lipschitz 連続でない一様連続関数]
$[0,\infty)$ 上の $\sqrt x$ を考える。$x\ge y\ge0$ ならば $\sqrt x-\sqrt y\le\sqrt x+\sqrt y$ なので

$$
(\sqrt x-\sqrt y)^2\le(\sqrt x-\sqrt y)(\sqrt x+\sqrt y)=x-y
$$

であり、$\abs{\sqrt x-\sqrt y}\le\sqrt{\abs{x-y}}$ が成り立つ。従って $\delta=\eps^2$ がすべての点で使え、$\sqrt x$ は一様連続である。
しかし $\abs{\sqrt x-\sqrt0}/\abs{x-0}=1/\sqrt x$ は $x\to0$ で有界でないので、$[0,1]$ 上でさえ Lipschitz 連続ではない。
:::

:::proposition[一様連続性の点列判定]{#prop-anl-uniform-sequences}
$f:D\to\R$ が一様連続であるための必要十分条件は、$D$ の任意の二列 $(x_n),(y_n)$ について、
$\abs{x_n-y_n}\to0$ ならば $\abs{f(x_n)-f(y_n)}\to0$ となることである。
:::

:::proof
一様連続ならば、$\eps>0$ に対応する $\delta$ を選ぶ。十分後で $\abs{x_n-y_n}<\delta$ となるから、関数値の差は $\eps$ 未満となる。
逆に一様連続でなければ、ある $\eps_0>0$ に対してすべての $\delta>0$ で反例の二点が存在する。
$\delta=1/n$ とするたびにその二点を $x_n,y_n$ に選ぶと、$\abs{x_n-y_n}<1/n$ でありながら関数値の差は常に $\eps_0$ 以上となり、点列の条件に反する。
:::

:::proposition[一様連続関数と Cauchy 列]{#prop-anl-uniform-cauchy-seq}
$f:D\to\R$ が一様連続ならば、$D$ の任意の Cauchy 列 $(x_n)$ について $(f(x_n))$ も Cauchy 列である。
:::

:::proof
$\eps>0$ に対して一様連続性の $\delta$ をとり、$m,n\ge N$ で $\abs{x_m-x_n}<\delta$ となる $N$ をとる。そのとき $\abs{f(x_m)-f(x_n)}<\eps$ である。
:::

各点で連続なだけではこの性質は成り立たない。$(0,1]$ 上の $1/x$ は連続であるが、Cauchy 列 $1/n$ を Cauchy 列でない $n$ に写す。従って $1/x$ は $(0,1]$ で一様連続でない。

### 有界閉区間の場合

:::theorem[Heine–Cantor の定理]{#thm-anl-heine-cantor}
有界閉区間 $[a,b]$ 上で連続な実数値関数は、この区間で一様連続である。
:::

:::proof
一様連続でないとすると、ある $\eps_0>0$ と $[a,b]$ 内の二列があり、

$$
\abs{x_n-y_n}<1/n,\qquad \abs{f(x_n)-f(y_n)}\ge\eps_0
$$

を満たす。:ref[thm-anl-bolzano-weierstrass] により $x_{n_k}\to c$ となる部分列を選べる。閉区間なので $c\in[a,b]$ である。
同じ添字で $y_{n_k}$ をとると、$\abs{y_{n_k}-c}\le1/n_k+\abs{x_{n_k}-c}\to0$ である。

$c$ における連続性により $f(x_{n_k}),f(y_{n_k})$ はともに $f(c)$ に収束する。
従って十分後で両方の値が $f(c)$ の $\eps_0/2$ 以内に入り、両者の差は $\eps_0$ 未満となる。これは上の選び方に反する。
:::

例えば区間を細かく分け、各小区間の中で関数値がどれだけ変わるかを抑えたいとする。
この定理によれば、連続関数に対しては、どの小区間にも共通の幅を誤差から選べる。
この事実は [Riemann 積分](/math/analysis/riemann-integral/) で連続関数の可積分性を示す際に用いる。

:::theorem[有界区間からの連続拡張]{#thm-anl-uniform-extension}
$a<b$ とし、$f:(a,b)\to\R$ は連続とする。次の三条件は同値である。

1. $f$ は $(a,b)$ で一様連続である。
2. 片側極限 $f(a+)$ と $f(b-)$ がともに有限な値として存在する。
3. $[a,b]$ 上の連続関数 $F$ で、$(a,b)$ 上で $f$ に一致するものが存在する。

このとき $F$ は一意であり、$F(a)=f(a+)$、$F(b)=f(b-)$ である。
:::

:::proof
**1 ⇒ 2。** $\eps>0$ に対して一様連続性の $\delta$ をとる。$x,y\in(a,a+\delta/2)$ ならば $\abs{x-y}<\delta$ なので $\abs{f(x)-f(y)}<\eps$ である。
:ref[thm-anl-function-cauchy] を $D=(a,b)$、$\alpha=a$ に適用すれば $f(a+)$ が存在する。$f(b-)$ も同様である。

**2 ⇒ 3。** $F(a)=f(a+)$、$F(b)=f(b-)$、内部では $F=f$ と定める。$a$ での $F$ の連続性は右極限が $F(a)$ に等しいことそのものであり、$b$ でも同様である。内部の点では $f$ の連続性による。

**3 ⇒ 1。** :ref[thm-anl-heine-cantor] により $F$ は $[a,b]$ で一様連続であり、その制限 $f$ も一様連続である。

一意性：$[a,b]$ 上の連続な拡張 $G$ について、$a$ での連続性から $G(a)=\lim_{x\to a+}G(x)=f(a+)$ である。$b$ でも同様なので $G=F$ である。
:::

例えば $(0,1]$ 上の $1/x$ は $f(0+)$ が存在しないので、$(0,1)$ でも一様連続でない。逆に $x\varphi(1/x)$（:ref[exm-anl-sawtooth]）は $(0,1)$ で一様連続である。

## 単調関数と逆関数

### 片側極限と跳び

区間 $I$ で、$x<y$ ならば $f(x)\le f(y)$ を満たす関数を :term[単調増加関数]{en="nondecreasing function" yomi="たんちょうぞうかかんすう"} という。
常に $f(x)<f(y)$ なら狭義増加という。向きを逆にして単調減少、狭義減少を定める。

:::proposition[単調関数の片側極限]{#prop-anl-monotone-limits}
$f:I\to\R$ が単調増加で、$c$ が区間 $I$ の内点ならば、両方の片側極限が存在し

$$
f(c-)=\sup_{\substack{x\in I\\x<c}} f(x)
\ \le f(c)\le\
\inf_{\substack{x\in I\\x>c}} f(x)=f(c+)
$$

である。特に $c$ で連続であるための必要十分条件は $f(c-)=f(c+)$ である。
:::

:::proof
左側の値の集合は空でなく $f(c)$ で上に有界である。その上限を $l$ とし、$\eps>0$ に対して $l-\eps<f(u)$ となる $u<c$ を選ぶ。
$\delta=c-u$ とすれば、$x\in I,\ c-\delta<x<c$ で $l-\eps<f(u)\le f(x)\le l$ である。従って左極限は $l$ である。
右側の値の集合は $f(c)$ で下に有界である。その下限を $r$ として $f(v)<r+\eps$ となる $v>c$ をとると、$c<x<v$ で $r\le f(x)\le f(v)<r+\eps$ となる。
ゆえに右極限は $r$ である。$l\le f(c)\le r$ なので、両極限が等しければそれらは $f(c)$ に等しく、左右から連続である。
逆に連続なら左右からの極限はいずれも $f(c)$ となる。
:::

端点 $c=\min I$ でも同じ議論で $f(c+)=\inf_{x>c}f(x)\ge f(c)$ が存在し、$c=\max I$ では $f(c-)\le f(c)$ が存在する。
従って単調関数の不連続点はすべて第一種であり、内点では跳びである。

:::example[原点での跳び]
$f(x)=x$（$x<0$）、$f(x)=x+1$（$x\ge0$）とおくと、$f$ は狭義増加である。
ところが $f(0-)=0,\ f(0)=f(0+)=1$ であり、像からは $[0,1)$ が抜ける。
従って単調であっても値が跳ぶことがあり、連続性を調べるには左右の極限が一致するかを確かめる必要がある。
:::

:::theorem[単調関数の不連続点]{#thm-anl-monotone-discontinuities}
区間 $I$ 上の単調関数の不連続点全体は高々可算である。すなわち、空集合であるか、$\N$ からの全射の像として表される。
:::

:::proof
$f$ を単調増加とする（減少なら $-f$ を考える）。$I$ の内点で $f$ が不連続な点の集合を $S$ とする。
$c\in S$ なら :ref[prop-anl-monotone-limits] により $f(c-)<f(c+)$ であり、:ref[thm-anl-rational-density] によって $f(c-)<r(c)<f(c+)$ となる有理数 $r(c)$ を一つ選ぶ。

$c<c'$ を $S$ の二点とし、$c<m<c'$ をとると

$$
f(c+)=\inf_{x>c}f(x)\le f(m)\le\sup_{x<c'}f(x)=f(c'-)
$$

であるから、開区間 $(f(c-),f(c+))$ と $(f(c'-),f(c'+))$ は共通部分をもたず、$r(c)\ne r(c')$ である。
従って $r:S\to\Q$ は単射である。$\Q$ は可算であるから（[可算集合](/math/sets-logic/cardinality/)）、$S$ は高々可算である。
不連続点は $S$ に $I$ の端点（高々二点）を加えたものに含まれるので、やはり高々可算である。
:::

この定理は、単調関数の不連続点が「少ない」ことを述べている。例えば Dirichlet 関数のようにすべての点で不連続な関数は、どの区間（一点でないもの）に制限しても単調にならない。そのような区間は :ref[thm-anl-uncountable] と同じ議論により可算でないからである。

### 像が区間であることと連続性

:::proposition[単調関数の連続性の判定]{#prop-anl-monotone-interval-image}
区間 $I$ 上の単調関数 $f$ について、像 $f(I)$ が区間ならば $f$ は連続である。
:::

:::proof
$f$ を単調増加とし、内点 $c$ で不連続と仮定する。:ref[prop-anl-monotone-limits] により $f(c-)<f(c)$ または $f(c)<f(c+)$ である。
前者の場合、$f(c-)<w<f(c)$ となる $w$ をとる。$x_0<c$ となる $x_0\in I$ について $f(x_0)\le f(c-)<w<f(c)$ であり、$f(x_0),f(c)\in f(I)$ と凸性から $w\in f(I)$ である。
しかし $x<c$ なら $f(x)\le f(c-)<w$、$x\ge c$ なら $f(x)\ge f(c)>w$ なので、$f(x)=w$ となる $x$ はない。これは矛盾である。
後者の場合は $f(c)<w<f(c+)$ となる $w$ と $x_1>c$ となる $x_1\in I$ をとる。$f(c)<w<f(c+)\le f(x_1)$ より $w\in f(I)$ であるが、$x\le c$ なら $f(x)\le f(c)<w$、$x>c$ なら $f(x)\ge f(c+)>w$ となり、やはり矛盾する。

端点 $c=\min I$ で不連続なら $f(c)<f(c+)$ であり、$f(c)<w<f(c+)$ となる $w$ について同じ矛盾が生じる。$c=\max I$ でも同様である。
:::

:::theorem[連続な単射は狭義単調]{#thm-anl-injective-monotone}
区間 $I$ 上の連続な単射 $f$ は狭義単調である。
:::

:::proof
$x<y$ と $x'<y'$ を $I$ の二組の点とし、$t\in[0,1]$ に対して

$$
x(t)=(1-t)x+tx',\qquad y(t)=(1-t)y+ty',\qquad h(t)=f(y(t))-f(x(t))
$$

とおく。$x(t),y(t)$ は $I$ の二点を結ぶ線分上にあるので凸性により $I$ に属し、$x(t)<y(t)$ である（$(1-t)(y-x)+t(y'-x')>0$）。
$f$ は単射なので $h(t)\ne0$ であり、$h$ は連続関数の合成として $[0,1]$ 上で連続である。
もし $h(0)$ と $h(1)$ の符号が異なれば、:ref[thm-anl-intermediate-value] により $h(t)=0$ となる $t$ が存在して矛盾する。
従って $f(y)-f(x)$ と $f(y')-f(x')$ は同符号である。すなわち、ある一組で $f$ が増加していればすべての組で増加しており、$f$ は狭義増加または狭義減少である。
:::

連続性は外せない。$[0,2]$ 上で $f(x)=x$（$0\le x\le1$）、$f(x)=4-x$（$1<x\le2$）とすると、像はそれぞれ $[0,1]$ と $[2,3)$ で交わらないので $f$ は単射であるが、単調でない。

### 逆関数が連続になる理由

$f:I\to\R$ が狭義単調なら、異なる二点は異なる値をもつ。従って像 $J=f(I)$ の各 $y$ にただ一つの原像 $x$ が対応する。
この対応を :term[逆関数]{en="inverse function" yomi="ぎゃくかんすう"} $f^{-1}:J\to I$ という。$f^{-1}$ も $f$ と同じ向きに狭義単調である。

:::theorem[単調関数の逆関数]{#thm-anl-continuous-inverse}
区間 $I$ 上の狭義単調関数 $f$ の逆関数は、像 $J=f(I)$ 上で連続である。
特に $f$ が連続でもあれば、$J$ は区間であり、$f$ は $I$ と $J$ の間の両方向に連続な全単射を与える。
:::

:::proof
$f$ が狭義増加の場合を考える。$y_0=f(x_0)$、$\eps>0$ を固定する。
$x_0$ より左に $I$ の点があれば、区間の凸性から $x_0-\eps<u<x_0$ となる $u\in I$ をとれる。
右に点があれば同様に $x_0<v<x_0+\eps$ となる $v\in I$ をとる。

$1$ と、存在する側についての差 $y_0-f(u)$、$f(v)-y_0$ の最小を $\delta$ とする。狭義増加性により $\delta>0$ である。
$y\in J$ が $\abs{y-y_0}<\delta$ を満たすとし、$x=f^{-1}(y)$ とおく。左側がある場合は $f(u)<y$ より $u<x$（$x\le u$ なら $y=f(x)\le f(u)$ となるため）、従って $x_0-\eps<x$ である。
左側がない場合も $x\ge x_0$ なので同じ不等式が成り立つ。右側についても、存在すれば $y<f(v)$ より $x<v<x_0+\eps$、存在しなければ $x\le x_0$ である。
以上より $\abs{x-x_0}<\eps$ となる。端点や一点だけの区間もこの議論に含まれる。

狭義減少の場合には $g=-f$ が狭義増加である。逆関数は $f^{-1}(y)=g^{-1}(-y)$ と表せるので連続である。
最後に $f$ が連続なら、:ref[cor-anl-continuous-image-interval] より像 $J$ は区間である。
:::

先の跳びをもつ例では、逆関数の定義域は $(-\infty,0)\cup[1,\infty)$ である。値 $1$ にこの定義域の中で近づくとき、十分近い点はすべて $1$ 以上である。
従って逆関数は $1$ でも連続である。逆関数をどの集合で考えているかがここで効いている。
:ref[thm-anl-injective-monotone] と合わせると、区間上の連続な単射は、像の区間への連続な逆関数をもつ。
区間上の連続性と単調性についての別の記述は :cite[sugiura-analysis] にもある。

例えば $x\mapsto x^n$ は $[0,\infty)$ で連続かつ狭義増加であり（:ref[thm-anl-nth-root] の段階 1）、像は $[0,\infty)$ である（$n$ 乗根の存在による）。
従って $n$ 乗根 $t\mapsto\sqrt[n]t$ は $[0,\infty)$ 上の連続関数である。

:::example[立方根の誤差を直接抑える]
$f(x)=x^3$ の逆関数を具体的に調べる。次の恒等式が、狭義増加性と誤差の評価を同時に与える。

$$
x^3-y^3=(x-y)(x^2+xy+y^2),\qquad
x^2+xy+y^2=\tfrac14(x-y)^2+\tfrac34(x+y)^2 .
$$

$x>y$ なら右の因子は正なので $x^3>y^3$ である。また任意の $u\in\R$ に対し $M=1+\abs{u}$ とすると $-M^3<u<M^3$ である。
連続な $f$ に中間値の定理を用いれば、$x^3=u$ を満たす実数があり、狭義増加性によりただ一つである。その $x$ を $\sqrt[3]{u}$ と書く。

上の恒等式から $\abs{x^3-y^3}\ge\abs{x-y}^3/4$ が従う。そこで $x=\sqrt[3]{u},\ y=\sqrt[3]{v}$ とおくと、
$\abs{u-v}<\eps^3/4$ ならば $\abs{x-y}^3\le4\abs{u-v}<\eps^3$、従って $\abs{x-y}<\eps$ となる。
従って $\delta=\eps^3/4$ はすべての $u,v$ に使える。立方根は原点を含む実数全体で一様連続である。
:::

## 演習問題

::::exercise{level=basic}
$f(x)=1/x$ の $a=2$ における連続性を、$\eps$ に対する具体的な $\delta$ を与えて示せ。
:::solution
$\abs{x-2}<1$ に限れば $x>1$ である。従って

$$
\abs{\frac1x-\frac12}=\frac{\abs{x-2}}{2\abs{x}}\le\frac{\abs{x-2}}{2}.
$$

そこで $\delta=\min\{1,2\eps\}$ と選ぶ。$\abs{x-2}<\delta$ なら右辺は $\delta/2\le\eps$ 未満なので、定義が満たされる。
:::
::::

::::exercise{level=basic}
有理数で $F(x)=x$、無理数で $F(x)=-x$ と定める。連続な点をすべて求めよ。
:::solution
$\abs{F(x)-F(0)}=\abs{x}$ なので、$a=0$ では $\delta=\eps$ と選べば連続である。
$a\ne0$ とする。:ref[thm-anl-rational-density] と :ref[cor-anl-irrational-density] により $a$ に収束する有理数列 $r_n$ と無理数列 $s_n$ をとると、$F(r_n)=r_n\to a$、$F(s_n)=-s_n\to-a$ である。
$a\ne-a$ なので、:ref[cor-anl-sequential-continuity] により $F$ は $a$ で連続でない。連続な点は $0$ のみである。
:::
::::

::::exercise{level=basic}
$f(x)=x^2-4x+1$ の $[0,3]$ での像を求めよ。最大値と最小値をどこでとるかも答えよ。
:::solution
平方完成すると $f(x)=(x-2)^2-3$ である。$x-2$ は $[-2,1]$ を動くので、その絶対値は $x=2$ で最小値 $0$、$x=0$ で最大値 $2$ をとる。
従って $f$ の最小値は $f(2)=-3$、最大値は $f(0)=1$ である。
:ref[cor-anl-interval-image] により像は $[-3,1]$ である。
:::
::::

::::exercise{level=basic}
床関数 $x\mapsto\lfloor x\rfloor$ と $g(x)=x-\lfloor x\rfloor$ の不連続点を求め、その種類を答えよ。
:::solution
整数でない $a$ では、$\delta=\min\{a-\lfloor a\rfloor,\lfloor a\rfloor+1-a\}>0$ とすれば $\abs{x-a}<\delta$ で $\lfloor x\rfloor=\lfloor a\rfloor$ である。従って床関数は $a$ の近くで定数であり連続である。$g$ も連続関数の差なので連続である。

整数 $k$ では、$k-1<x<k$ で $\lfloor x\rfloor=k-1$、$k\le x<k+1$ で $\lfloor x\rfloor=k$ なので、左極限は $k-1$、右極限と値は $k$ である。
従って各整数は跳びのある第一種の不連続点であり、床関数はそこで右連続である。
$g$ については左極限が $k-(k-1)=1$、右極限と値が $0$ であり、やはり各整数で跳び（大きさ $-1$）をもち右連続である。
:::
::::

::::exercise{level=basic}
$\lim_{x\to+\infty}\bigl(\sqrt{x^2+x}-x\bigr)$ を求めよ。
:::solution
$x>0$ で分子を有理化すると

$$
\sqrt{x^2+x}-x=\frac{x}{\sqrt{x^2+x}+x}=\frac1{\sqrt{1+1/x}+1}
$$

である（$\sqrt{x^2+x}=x\sqrt{1+1/x}$ は、両辺が非負で平方が等しいことによる）。
$t\ge0$ について $\abs{\sqrt{1+t}-1}=t/(\sqrt{1+t}+1)\le t$ であるから、$0\le\sqrt{1+1/x}-1\le1/x$ となり、$x\to+\infty$ で $\sqrt{1+1/x}\to1$ である。
これは $1+1/x\to1$ と平方根の連続性から :ref[prop-anl-composite-limit] の 1 によっても得られる。
従って極限の演算法則（:ref[prop-anl-sequential-extended] により :ref[thm-anl-limit-laws] が $x\to+\infty$ の極限に移る）から、極限は $1/(1+1)=1/2$ である。
:::
::::

::::exercise{level=advanced}
連続関数 $f:[0,1]\to[0,1]$ について $f(c)=c$ となる点が存在することを示せ。
さらに定義域と値域をどちらも $(0,1)$ にすると、この主張が成り立たない例を与えよ。
:::solution
$g(x)=f(x)-x$ は連続であり、値域の仮定から $g(0)=f(0)\ge0,\ g(1)=f(1)-1\le0$ となる。
端で等号があればそこで $f(c)=c$ であり、いずれも等号でなければ中間値の定理で $g(c)=0$ となる内点が得られる。

一方、$(0,1)$ 上で $f(x)=(x+1)/2$ とすると $x<f(x)<1$ である。よってこの連続関数は $(0,1)$ を $(0,1)$ に写し、どの点でも $f(x)\ne x$ である。
:::
::::

::::exercise{level=advanced}
$f(x)=1/x$ は $(0,1]$ では一様連続でなく、任意の $c>0$ に対して $[c,\infty)$ では一様連続であることを示せ。
:::solution
前者では $x_n=1/n,\ y_n=1/(n+1)$ ととる。入力の差は $1/(n(n+1))\to0$ であるが、出力の差の絶対値は $1$ なので :ref[prop-anl-uniform-sequences] の条件に反する。
後者では $x,y\ge c$ より $xy\ge c^2$ なので

$$
\abs{\frac1x-\frac1y}=\frac{\abs{x-y}}{xy}\le\frac{\abs{x-y}}{c^2}
$$

であり、$1/x$ は $[c,\infty)$ で Lipschitz 連続である。従ってどの二点にも共通に $\delta=c^2\eps$ が使える。
:::
::::

::::exercise{level=advanced}
連続関数 $f:\R\to\R$ がすべての $x,y$ で $f(x+y)=f(x)+f(y)$ を満たすならば、$c=f(1)$ として $f(x)=cx$ であることを示せ。
:::solution
$x=y=0$ とすると $f(0)=2f(0)$ より $f(0)=0$ であり、$y=-x$ とすると $f(-x)=-f(x)$ である。
帰納法により自然数 $n$ について $f(nx)=nf(x)$ であり、奇関数であることから整数 $n$ についても成り立つ。
$m\in\Z$、$n\in\N$ について $nf(m/n)=f(m)=mf(1)$ なので $f(m/n)=c\,m/n$、すなわち有理数 $r$ で $f(r)=cr$ である。

任意の実数 $x$ に対し、:ref[thm-anl-rational-density] により $x$ に収束する有理数列 $r_n$ がとれる（$x-1/n<r_n<x$）。
連続性と :ref[cor-anl-sequential-continuity] により $f(x)=\lim f(r_n)=\lim cr_n=cx$ である。
:::
::::

::::exercise{level=advanced}
$D\subset\R$ が有界で、$f:D\to\R$ が一様連続ならば、$f$ は有界であることを示せ。また $D$ が有界でないと成り立たない例を挙げよ。
:::solution
$D\ne\varnothing$ としてよい。$\eps=1$ に対する一様連続性の $\delta$ をとる。$D\subset[\alpha,\beta]$ とし、Archimedes の性質により $K>(\beta-\alpha)/\delta$ となる自然数 $K$ を選んで、$[\alpha,\beta]$ を長さ $(\beta-\alpha)/K<\delta$ の $K$ 個の閉区間 $J_1,\dots,J_K$ に分ける。
$D\cap J_i$ が空でない各 $i$ について点 $x_i\in D\cap J_i$ を一つ選ぶ。
任意の $x\in D$ はある $J_i$ に属し、$\abs{x-x_i}<\delta$ なので $\abs{f(x)}<\abs{f(x_i)}+1$ である。
従って $\abs{f(x)}\le1+\max_i\abs{f(x_i)}$（有限個の最大値）となり、$f$ は有界である。

$D=\R$、$f(x)=x$ は Lipschitz 連続なので一様連続であるが、有界でない。
:::
::::

::::exercise{level=advanced}
$f:\R\to\R$ は連続で $\lim_{x\to+\infty}f(x)=\lim_{x\to-\infty}f(x)=0$ を満たし、ある $x_0$ で $f(x_0)>0$ とする。$f$ は $\R$ 上で最大値をとることを示せ。
:::solution
$\eps=f(x_0)/2>0$ とする。極限の仮定から、ある $R>0$ があって $\abs x>R$ なら $\abs{f(x)}<\eps$ となる。$R$ は $\abs{x_0}$ より大きくとってよい。
:ref[thm-anl-extreme-value] により $f$ は $[-R,R]$ 上で最大値 $f(q)$ をとり、$x_0\in[-R,R]$ なので $f(q)\ge f(x_0)=2\eps$ である。
$\abs x>R$ なら $f(x)<\eps<f(q)$ であるから、$f(q)$ は $\R$ 全体での最大値である。
:::
::::

## 確認テスト

::::quiz
$a$ は $D$ の集積点である。$\lim_{x\to a}f(x)=L$ と同値なものはどれか。

- [ ] $a$ に近づく一つの列で関数値が $L$ に近づく。
- [x] $D\setminus\{a\}$ で $a$ に近づくすべての列で関数値が $L$ に近づく。
- [ ] $f(a)$ が定義され、その値が $L$ である。

:::explanation
(B) は :ref[thm-anl-sequential-limit] そのものである。(A) は不十分で、例えば $x/\abs x$ は列 $1/n$ に沿って $1$ に近づくが、列 $-1/n$ に沿っては $-1$ に近づき、$0$ での極限をもたない。
(C) は誤りで、極限は $0<\abs{x-a}$ の点だけで決まり、$a$ での値には依存しない（定義されている必要もない）。
:::
::::

::::quiz{answer="0.2" tol="1e-9"}
$f(x)=5x-1$ に対して、$\abs{x-y}<\delta$ ならば $\abs{f(x)-f(y)}<1$ がすべての実数 $x,y$ で成り立つ。使える $\delta$ の最大値を求めよ。

:::explanation
出力の差は入力の差の $5$ 倍、すなわち $\abs{f(x)-f(y)}=5\abs{x-y}$ なので、$\delta=1/5$ なら $\abs{x-y}<1/5$ から $\abs{f(x)-f(y)}<1$ が従う。
これより大きい幅 $\delta>1/5$ では、距離がちょうど $1/5$ の二点が条件 $\abs{x-y}<\delta$ に入り、その出力の差が $1$ となって条件に反する。従って最大値は $0.2$ である。
:::
::::

::::quiz{answer="0.25" tol="1e-9"}
:ref[exm-anl-thomae] の関数で $T(6/8)$ を求めよ。

:::explanation
$6/8=3/4$ であり、$n\cdot\tfrac34\in\Z$ となる最小の自然数は $n=4$ である（$n=1,2,3$ では $3/4,3/2,9/4$ で整数でない）。従って $q(6/8)=4$、$T(6/8)=1/4$ である。
:::
::::

::::quiz
連続な $f:[a,b]\to\R$（$a<b$）について必ず成り立つものをすべて選べ。

- [x] 最大値をとる。
- [ ] 最大値をとる点は一つである。
- [x] 一様連続である。
- [x] 像は閉区間である。

:::explanation
最大値の存在は :ref[thm-anl-extreme-value]、一様連続性は :ref[thm-anl-heine-cantor]、像については :ref[cor-anl-interval-image] による（定数関数では一点 $[m,m]$ になる）。
定数関数はすべての点で最大値をとるので、最大点の一意性は従わない。
:::
::::

::::quiz{answer="0.01" tol="1e-9"}
$[0,1]$ 上の $f(x)=\sqrt x$ について、「$\abs{x-y}<\delta$ ならば $\abs{\sqrt x-\sqrt y}<1/10$」がすべての $x,y\in[0,1]$ で成り立つような $\delta$ の最大値を求めよ。

:::explanation
本文の例により $\abs{\sqrt x-\sqrt y}\le\sqrt{\abs{x-y}}$ なので、$\delta=1/100$ なら $\abs{\sqrt x-\sqrt y}<\sqrt{1/100}=1/10$ である。
$\delta>1/100$ なら $y=0$、$x=1/100$ は $\abs{x-y}<\delta$ を満たすが $\sqrt x-\sqrt y=1/10$ となり条件に反する。従って最大値は $0.01$ である。
:::
::::

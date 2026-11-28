---
title: Riemann 面上の微分形式と Dolbeault の補題
description: 複素微分形式の型、周期と留数、Cauchy 核による局所可解性を扱う
notation:
  - symbol: '\mathcal A^{p,q}_X'
    meaning: 滑らかな型 (p,q) の複素微分形式の層
  - symbol: '\partial,\bar\partial'
    meaning: 外微分の型 (1,0) と (0,1) 方向の成分
  - symbol: '\operatorname{Res}_p\eta'
    meaning: 有理型微分 η の点 p での留数
  - symbol: '\Delta'
    meaning: 局所座標では ∂²/∂x² + ∂²/∂y²
---

## この章の目標

- 正則座標の変換に対して微分形式の型が保たれることを示す。
- 微分の位数、留数、周期を局所展開と線積分から計算する。
- Cauchy 核を用いて $\bar\partial$ 方程式を局所的に解く。
- 弱い意味の調和性から滑らかさを導く。

## 複素構造による微分形式の型分解

$z=x+iy$ を正則座標とすると、$\dd z=\dd x+i\dd y$、$\dd\bar z=\dd x-i\dd y$ は複素化した余接空間の基底である。正則座標 $w$ に取り替えると $\dd w=w'(z)\dd z$、$\dd\bar w=\overline{w'(z)}\dd\bar z$ となる。従って二つの直線が混ざらない。

:::definition[型と Wirtinger 微分]{#def-cg-form-type}
$\dd z$ で張られる形式を :term[型]{en="type" yomi="かた"} $(1,0)$、$\dd\bar z$ で張られる形式を型 $(0,1)$ と呼ぶ。関数は型 $(0,0)$、$\dd z\wedge\dd\bar z$ の倍数は型 $(1,1)$ である。型 $(p,q)$ の滑らかな形式全体を開集合ごとに集めた層を $\mathcal A_X^{p,q}$ と書く。外微分を型で分け、$p$ を増やす成分を $\partial$、$q$ を増やす成分を $\bar\partial$ とする。
:::

基底を反転すると
$$
\begin{aligned}
u_z&=\tfrac12(u_x-iu_y), & u_{\bar z}&=\tfrac12(u_x+iu_y),\\
\dd u&=u_z\dd z+u_{\bar z}\dd\bar z, &
\dd z\wedge\dd\bar z&=-2i\,\dd x\wedge\dd y.
\end{aligned}
\label{eq-cg-wirtinger}
$$
正則座標の実 Jacobian は $|w'|^2>0$ だから、$\dd x\wedge\dd y$ の正の向きは座標によらない。一般の $1$ 形式 $\alpha=a\dd z+b\dd\bar z$ について積の微分を適用すれば、
$$
\dd\alpha=(b_z-a_{\bar z})\,\dd z\wedge\dd\bar z.
$$
$\dd^2=0$ を型で比べれば $\partial^2=\bar\partial^2=0$ と $\partial\bar\partial+\bar\partial\partial=0$ を得る。ここで $u$ の混合偏微分は等しいが、外積の順序の交換により作用素に負号が付く。

:::example[型を分けると相殺が見える]
$\alpha=\bar z\dd z+z\dd\bar z=\dd|z|^2$ は閉じている。実際、$\bar\partial(\bar z\dd z)=-\dd z\wedge\dd\bar z$、$\partial(z\dd\bar z)=\dd z\wedge\dd\bar z$ になる。閉じた形式の二つの型成分が個別に閉じているとは限らない。
:::

:::example[面積形式と座標の倍率]
$w=2z$ では $(i/2)\dd w\wedge\dd\bar w=4\dd x\wedge\dd y$ である。正向きが座標で一致しても、面積要素には座標の倍率が入る。
:::

### 共役と Hodge 星作用素

$1$ 形式の複素共役は型を交換する。従って実 $1$ 形式は $a\dd z+\bar a\dd\bar z$ と一意に書ける。複素多様体の正則座標で共形な計量 $g=\lambda^2(\dd x^2+\dd y^2)$ を取る。Hodge 星作用素を $\alpha\wedge *\beta=g(\alpha,\beta)\dd V_g$ で定め、複素線形に延長する。$g(\dd x,\dd x)=\lambda^{-2}$、$\dd V_g=\lambda^2\dd x\wedge\dd y$ だから $*\dd x=\dd y$、$*\dd y=-\dd x$ となる。$1$ 形式上では $\lambda$ が相殺し、$*\dd z=-i\dd z$、$*\dd\bar z=i\dd\bar z$ となる。このため型分解は、向きと共形構造だけからも取り出せる。

:::example[実形式から正則形式を取り出す]
実 $1$ 形式 $\alpha$ について $\alpha+i*\alpha$ は型 $(1,0)$ である。$\alpha$ と $*\alpha$ が閉じていればこの形式も閉じる。型 $(1,0)$ では外微分の $(1,1)$ 成分が $\bar\partial$ なので、これは正則形式である。逆に正則形式 $\eta$ は実部 $\alpha$ を持ち、$*\alpha=\Im\eta$ である。例えば $\eta=z\dd z$ なら $\alpha=x\dd x-y\dd y=\dd((x^2-y^2)/2)$、$*\alpha=y\dd x+x\dd y=\dd(xy)$ となる。
:::

### 反変換と実形式

座標を $w=h(z)$ とすると、$\mathrm dw=h'(z)\mathrm dz$、$\mathrm d\bar w=\overline{h'(z)}\mathrm d\bar z$ である。例えば $(1,0)$ 形式 $a_z\mathrm dz=a_w\mathrm dw$ の係数には $a_w=a_z/h'$ が必要となる。$(1,1)$ 形式の係数は $|h'|^{-2}$ 倍されるので、正の面積密度という条件は座標によらない。$\mathrm dz\wedge\mathrm d\bar z=-2i\,\mathrm dx\wedge\mathrm dy$ であるため、実形式は $i\,b\,\mathrm dz\wedge\mathrm d\bar z$ の形で $b$ が実数であるものに限る。

:::example[平面の向きを確かめる]
$w=iz$ なら $\mathrm dw=i\mathrm dz$、$\mathrm d\bar w=-i\mathrm d\bar z$ である。$\mathrm d(\Re w)\wedge\mathrm d(\Im w)=(-\mathrm dy)\wedge\mathrm dx=\mathrm dx\wedge\mathrm dy$ となり、座標回転は向きを保つ。$w=\bar z$ では $\mathrm dw=\mathrm d\bar z$ と型が入れ替わり、実面積の符号も反転する。後者は複素座標の正則な遷移に使えない。
:::

:::example[反転で係数に極が現れる]
$w=1/z$ のとき $\mathrm dz=-w^{-2}\mathrm dw$ である。平面で定数係数の $\mathrm dz$ にも球面では無限遠で二位の極がある。実面積 $\frac i2\mathrm dz\wedge\mathrm d\bar z$ は $\frac i{2|w|^4}\mathrm dw\wedge\mathrm d\bar w$ となる。一方、これに $(1+|z|^2)^{-2}$ を掛けた形式の $w$ 表示は $\frac i{2(1+|w|^2)^2}\mathrm dw\wedge\mathrm d\bar w$ なので滑らかに延びる。
:::

### 1 形式の閉性を成分で判定する

$\alpha=a\mathrm dz+b\mathrm d\bar z$ と置くと
$$
\mathrm d\alpha=(b_z-a_{\bar z})\,\mathrm dz\wedge\mathrm d\bar z.
$$
従って閉性は $a_{\bar z}=b_z$ であり、それぞれが零であることを要求しない。一方、型が $(1,0)$ で閉なら $a_{\bar z}=0$ となるので正則である。$\mathrm d^2=0$ は関数に適用すれば $f_{z\bar z}=f_{\bar zz}$ に帰着し、これを型ごとに書くと $\partial^2=\bar\partial^2=0$ と $\partial\bar\partial+\bar\partial\partial=0$ が得られる。

:::example[閉形式だが各成分は閉じない]
$u=z^2\bar z+\bar z^2z$ と置く。$a=2z\bar z+\bar z^2$、$b=z^2+2z\bar z$ に対し $\mathrm du=a\mathrm dz+b\mathrm d\bar z$ である。$a_{\bar z}=b_z=2(z+\bar z)$ なので全体は閉じる。各成分の外微分は互いに逆であり、一般の点では零でない。
:::

### 局所の座標変換と余接束の貼り合わせ

二枚のチャート $U,V$ の座標 $z,w$ が $w=f(z)$ でつながるとする。正則余接束の枠を $\mathrm dz,\mathrm dw$ と取れば、$V$ から $U$ への係数変換は $f'$ を掛ける操作である。例えば $g(w)\mathrm dw$ の $U$ 表示は $(g\circ f)f'\mathrm dz$ になる。三枚の重なりでは連鎖律 $(h\circ f)'=(h'\circ f)f'$ が貼り合わせの条件を与える。接ベクトルでは逆の倍率になる。実際、$V$ における接ベクトル $\partial_w$ に対し $\mathrm dw(\partial_w)=1$ であるから、$U$ 表示は $(f')^{-1}\partial_z$ である。

:::example[ベクトルと微分の縮約]
$w=(z-1)/(z+1)$ を $z=2$ で使うと $w=1/3$、$f'=2/9$ である。$z^2\partial_z$ はこの点で $4\cdot(2/9)\partial_w=(8/9)\partial_w$ となる。一方、$\eta=z^{-1}\mathrm dz$ は同じ点で $(1/2)(9/2)\mathrm dw=(9/4)\mathrm dw$ となる。縮約は $z$ 側で $4/2=2$、$w$ 側で $(8/9)(9/4)=2$ で一致する。
:::

:::example[共役の係数の変換]
$w=z^2$ を $z=1+i$ の近傍の座標変換として使う。この点では $w=2i$、$|f'|^2=8$ である。$\eta=\bar z\,\mathrm d\bar z$ の $w$ での係数は $\bar z/(2\bar z)=1/2$ である。実二次形式 $i(1+|z|^2)\mathrm dz\wedge\mathrm d\bar z$ の $w$ での係数は $(1+2)/8=3/8$ となる。第一の形式は反正則の Jacobian だけを使い、第二の形式は正則・反正則の双方の Jacobian を使う。
:::

## 正則 1 形式と有理型微分

:::definition[正則微分と有理型微分]{#def-cg-holomorphic-differential}
型 $(1,0)$ の形式 $\eta=a(z)\dd z$ が :term[正則微分]{en="holomorphic differential" yomi="せいそくびぶん"} であるとは、係数 $a$ が正則であることをいう。係数が有理型なら :term[有理型微分]{en="meromorphic differential" yomi="ゆうりけいびぶん"} という。正則微分の層は $\Omega_X^1$ または $K_X$ と書く。零でない微分の位数 $\operatorname{ord}_p\eta$ は、$a$ の位数である。
:::

係数は座標の逆変換の導関数だけを余分に掛けられ、その導関数は零にならない。よって定義は座標によらない。$\bar\partial\eta=0$ は係数の Cauchy–Riemann 方程式である。$\partial\eta=0$ は型 $(2,0)$ が存在しないことから従う。したがって正則微分は閉じている。円板では係数の Taylor 級数を積分して原始関数が得られる。

有理型関数 $f=z^k u(z)$、$u(0)\ne0$ に対し、$k\ne0$ なら $\dd f$ の位数は $k-1$ である。$k=0$ では $u'$ を調べる必要がある。また Laurent 展開を項別微分しても $\dd z/z$ の係数は生じない。対数微分については $\dd f/f=k\,\dd z/z+\dd u/u$ となる。

:::example[トーラスに降りる微分]
$X=\C/(\Z+2i\Z)$ では $\dd z$ が平行移動で不変なので降りる。逆に正則微分の係数を $\C$ へ持ち上げると二重周期の整関数になり、基本平行四辺形で有界である。Liouville の定理で定数だから正則微分の空間は $\C\dd z$ である。
:::

:::example[二枚被覆での微分]
$y^2=x(x-1)(x-2)$ の有限分岐点では $y$ が座標で、$2y\dd y=P'(x)\dd x$ から $\dd x/y=2\dd y/P'(x)$ は正則かつ非零である。無限遠では $x$ が二位、$y$ が三位の極なので $\dd x/y$ の位数は $-3-(-3)=0$ である。したがってこの微分は全域で零をもたない。
:::

### 超楕円曲線で正則な微分を見つける

平方自由な多項式 $P(x)$ の次数を $2g+1$ または $2g+2$ とし、$y^2=P(x)$ の滑らかな完成面を考える。有限の非分岐点では $\dd x/y$ は $x$ 座標で非零正則である。有限分岐点でも $2\dd y/P'(x)$ という式により同じ性質を持つ。違いが出るのは無限遠である。

奇数次数 $2g+1$ のときは $x=t^{-2}$、$y=t^{-(2g+1)}v(t)$、$v(0)\ne0$ と置ける。従って $\dd x/y=-2t^{2g-2}\dd t/v(t)$ である。偶数次数 $2g+2$ のときは $s=1/x$ を無限遠の各点の座標にでき、$y=s^{-(g+1)}v(s)$ と書ける。このとき $\dd x/y=-s^{g-1}\dd s/v(s)$ である。従って両場合とも $x^j\dd x/y$ は $0\le j\le g-1$ で正則となる。

これらがすべてを張ることも局所計算から示せる。反転 $\iota:(x,y)\mapsto(x,-y)$ で不変な正則微分は $x$ 球面へ降りる。分岐の近くで $\iota(t)=-t$ とすれば、不変条件は係数 $h(t)$ が奇関数となることである。従って $h(t)\dd t=tH(t^2)\dd t=\frac12H(s)\dd s$ と書け、底でも正則である。球面上の正則微分は零である。実際、有限座標で $f(z)\dd z$ と書くと $f$ は整関数で、無限遠での正則性から $O(|z|^{-2})$ であり、Liouville の定理を適用できる。このため任意の正則微分は $\iota$ で符号を変える。$\dd x/y$ で割った比は不変な有理型関数で、分岐座標では Laurent 展開が偶数べきだけになるため底に降りる。有限点で極がないので比は多項式である。無限遠の位数の条件から次数は $g-1$ 以下となる。

:::example[奇数次の種数二曲線]
$y^2=x^5-1$ では $x=t^{-2},y=t^{-5}(1-t^{10})^{1/2}$ である。$\dd x/y$ の無限遠での位数は $2$、$x\dd x/y$ では $0$ である。$x^2\dd x/y$ になると位数は $-2$ なので正則微分にはならない。
:::

### 分岐座標での引き戻しと微分の因子

$f:X\to Y$ を非定数正則写像とし、$p$ で局所次数 $e$ とする。局所標準形で値域座標 $w=z^e$ を使い、微分を $w^mu(w)\mathrm dw$、$u(0)\ne0$ と書く。引き戻しは $e z^{em+e-1}u(z^e)\mathrm dz$ なので位数は $e(m+1)-1$ である。この式を全点に適用すると
$$
(f^*\eta)=f^*(\eta)+R_f,\qquad R_f=\sum_p(e_p(f)-1)[p]
$$
となる。左辺は微分の因子、右辺の第一項は因子の引き戻しである。$m=-1$ の場合は分岐しても単純極のままである。$m\ge0$ の場合、微分が値域で非零でも分岐点では零になりうる。

:::example[超楕円微分が分岐で消えない理由]
$y^2=x(x-1)(x-3)$ の $x=0$ の上で $y=t$ を座標に取る。右辺の $x$ 微分は原点で $3$ だから $x=t^2/3+O(t^4)$ である。従って $\mathrm dx/y=(2/3+O(t^2))\mathrm dt$ は正則かつ非零である。$\mathrm dx$ 自体には位数 $1$ の零があり、それを $y$ の零でちょうど割っている。
:::

:::example[偶数次数の超楕円曲線の無限遠]
$y^2=x^8+1$ の無限遠では $s=1/x$、$y=s^{-4}v$、$v^2=1+s^8$ である。$v(0)=1,-1$ の二枝で $\mathrm dx/y=-s^2\mathrm ds/v$ となる。従って $\mathrm dx/y,x\mathrm dx/y,x^2\mathrm dx/y$ の位数は各枝で $2,1,0$ である。有限の分岐点では先の局所計算を使えるため、三つは大域正則微分である。
:::

:::example[奇数次数では無限遠が分岐する]
$y^2=x^7-1$ では $x=t^{-2}$、$y=t^{-7}v$、$v^2=1-t^{14}$ と置ける。$v(0)=1$ を取れば $t$ と $-t$ が二根を動かし、無限遠上には一点しかない。$x^j\mathrm dx/y=-2t^{4-2j}\mathrm dt/v$ だから $j=0,1,2$ のとき位数は $4,2,0$ となる。
:::

### 微分の位数を局所次数から計算する

:::example[3 枚の根の被覆を完成する]
$y^3=x \left(x - 3\right) \left(x - 2\right) \left(x - 1\right)$ の根の被覆を取り、完成で加わる点と各点での局所次数を調べる。

根に $e^{2\pi i/3}$ を掛けて番号を付ける。指数 $m$ の点では、非消失因子の $3$ 乗根を選んで除く。一周は番号を $m$ だけ進め、$d=\gcd(3,m)$ 個の長さ $e=3/d$ の軌道に分ける。無限遠では $s=1/x$ として指数の和の負を使う。

$$
\begin{aligned}
0\colon\quad m=1,\quad d=\gcd(3,1)=1,\quad e=3,\quad R=2 \\
1\colon\quad m=1,\quad d=\gcd(3,1)=1,\quad e=3,\quad R=2 \\
2\colon\quad m=1,\quad d=\gcd(3,1)=1,\quad e=3,\quad R=2 \\
3\colon\quad m=1,\quad d=\gcd(3,1)=1,\quad e=3,\quad R=2 \\
\infty\colon\quad m=-4,\quad d=\gcd(3,-4)=1,\quad e=3,\quad R=2
\end{aligned}
$$

完成では各軌道に一点を加え、底の中心座標を $t^e$ とする。全体は $1$ 成分であり、番号を一つ進める置換は各成分を入れ替える。従って各成分の次数は $3$、分岐寄与は全寄与を $1$ で割った $10$ である。分岐値ごとの補正は $d$ 個の逆像を合わせた $d(e-1)$ である。

$$
\begin{aligned}
2g-2&=-2\cdot 3+10=4 \\
g&=3
\end{aligned}
$$

有限点では $x-a=t^e$、無限遠では $x=t^{-e}$ を微分する。$y$ の位数を差し引くと次の表となる。

$$
\begin{aligned}
0\colon&\quad \operatorname{ord}\mathrm dx=2,\quad \operatorname{ord}y=1,\quad\operatorname{ord}(\mathrm dx/y)=1. \\
1\colon&\quad \operatorname{ord}\mathrm dx=2,\quad \operatorname{ord}y=1,\quad\operatorname{ord}(\mathrm dx/y)=1. \\
2\colon&\quad \operatorname{ord}\mathrm dx=2,\quad \operatorname{ord}y=1,\quad\operatorname{ord}(\mathrm dx/y)=1. \\
3\colon&\quad \operatorname{ord}\mathrm dx=2,\quad \operatorname{ord}y=1,\quad\operatorname{ord}(\mathrm dx/y)=1. \\
\infty\colon&\quad \operatorname{ord}\mathrm dx=-4,\quad \operatorname{ord}y=-4,\quad\operatorname{ord}(\mathrm dx/y)=0.
\end{aligned}
$$

表にない値では $x$ が局所座標で $y$ は非消失なので、$\mathrm dx/y$ は正則かつ非零である。従ってこの微分の因子は表から完全に決まる。
:::

:::example[4 枚の根の被覆を完成する]
$y^4=x \left(x - 2\right)^{3} \left(x - 1\right)^{2}$ の根の被覆を取り、完成で加わる点と各点での局所次数を調べる。

根に $e^{2\pi i/4}$ を掛けて番号を付ける。指数 $m$ の点では、非消失因子の $4$ 乗根を選んで除く。一周は番号を $m$ だけ進め、$d=\gcd(4,m)$ 個の長さ $e=4/d$ の軌道に分ける。無限遠では $s=1/x$ として指数の和の負を使う。

$$
\begin{aligned}
0\colon\quad m=1,\quad d=\gcd(4,1)=1,\quad e=4,\quad R=3 \\
1\colon\quad m=2,\quad d=\gcd(4,2)=2,\quad e=2,\quad R=2 \\
2\colon\quad m=3,\quad d=\gcd(4,3)=1,\quad e=4,\quad R=3 \\
\infty\colon\quad m=-6,\quad d=\gcd(4,-6)=2,\quad e=2,\quad R=2
\end{aligned}
$$

完成では各軌道に一点を加え、底の中心座標を $t^e$ とする。全体は $1$ 成分であり、番号を一つ進める置換は各成分を入れ替える。従って各成分の次数は $4$、分岐寄与は全寄与を $1$ で割った $10$ である。分岐値ごとの補正は $d$ 個の逆像を合わせた $d(e-1)$ である。

$$
\begin{aligned}
2g-2&=-2\cdot 4+10=2 \\
g&=2
\end{aligned}
$$

有限点では $x-a=t^e$、無限遠では $x=t^{-e}$ を微分する。$y$ の位数を差し引くと次の表となる。

$$
\begin{aligned}
0\colon&\quad \operatorname{ord}\mathrm dx=3,\quad \operatorname{ord}y=1,\quad\operatorname{ord}(\mathrm dx/y)=2. \\
1\colon&\quad \operatorname{ord}\mathrm dx=1,\quad \operatorname{ord}y=1,\quad\operatorname{ord}(\mathrm dx/y)=0. \\
2\colon&\quad \operatorname{ord}\mathrm dx=3,\quad \operatorname{ord}y=3,\quad\operatorname{ord}(\mathrm dx/y)=0. \\
\infty\colon&\quad \operatorname{ord}\mathrm dx=-3,\quad \operatorname{ord}y=-3,\quad\operatorname{ord}(\mathrm dx/y)=0.
\end{aligned}
$$

表にない値では $x$ が局所座標で $y$ は非消失なので、$\mathrm dx/y$ は正則かつ非零である。従ってこの微分の因子は表から完全に決まる。
:::

:::example[5 枚の根の被覆を完成する]
$y^5=x^{2} \left(x - 1\right) \left(x + 1\right)^{2}$ の根の被覆を取り、完成で加わる点と各点での局所次数を調べる。

根に $e^{2\pi i/5}$ を掛けて番号を付ける。指数 $m$ の点では、非消失因子の $5$ 乗根を選んで除く。一周は番号を $m$ だけ進め、$d=\gcd(5,m)$ 個の長さ $e=5/d$ の軌道に分ける。無限遠では $s=1/x$ として指数の和の負を使う。

$$
\begin{aligned}
0\colon\quad m=2,\quad d=\gcd(5,2)=1,\quad e=5,\quad R=4 \\
1\colon\quad m=1,\quad d=\gcd(5,1)=1,\quad e=5,\quad R=4 \\
-1\colon\quad m=2,\quad d=\gcd(5,2)=1,\quad e=5,\quad R=4 \\
\infty\colon\quad m=-5,\quad d=\gcd(5,-5)=5,\quad e=1,\quad R=0
\end{aligned}
$$

完成では各軌道に一点を加え、底の中心座標を $t^e$ とする。全体は $1$ 成分であり、番号を一つ進める置換は各成分を入れ替える。従って各成分の次数は $5$、分岐寄与は全寄与を $1$ で割った $12$ である。分岐値ごとの補正は $d$ 個の逆像を合わせた $d(e-1)$ である。

$$
\begin{aligned}
2g-2&=-2\cdot 5+12=2 \\
g&=2
\end{aligned}
$$

有限点では $x-a=t^e$、無限遠では $x=t^{-e}$ を微分する。$y$ の位数を差し引くと次の表となる。

$$
\begin{aligned}
0\colon&\quad \operatorname{ord}\mathrm dx=4,\quad \operatorname{ord}y=2,\quad\operatorname{ord}(\mathrm dx/y)=2. \\
1\colon&\quad \operatorname{ord}\mathrm dx=4,\quad \operatorname{ord}y=1,\quad\operatorname{ord}(\mathrm dx/y)=3. \\
-1\colon&\quad \operatorname{ord}\mathrm dx=4,\quad \operatorname{ord}y=2,\quad\operatorname{ord}(\mathrm dx/y)=2. \\
\infty\colon&\quad \operatorname{ord}\mathrm dx=-2,\quad \operatorname{ord}y=-1,\quad\operatorname{ord}(\mathrm dx/y)=-1.
\end{aligned}
$$

表にない値では $x$ が局所座標で $y$ は非消失なので、$\mathrm dx/y$ は正則かつ非零である。従ってこの微分の因子は表から完全に決まる。
:::

### 有限点と無限遠の位数を比較する

:::example[3 枚の根の被覆を完成する]
$y^3=x \left(x - 3\right) \left(x - 1\right)^{2}$ の根の被覆を取り、完成で加わる点と各点での局所次数を調べる。

根に $e^{2\pi i/3}$ を掛けて番号を付ける。指数 $m$ の点では、非消失因子の $3$ 乗根を選んで除く。一周は番号を $m$ だけ進め、$d=\gcd(3,m)$ 個の長さ $e=3/d$ の軌道に分ける。無限遠では $s=1/x$ として指数の和の負を使う。

$$
\begin{aligned}
0\colon\quad m=1,\quad d=\gcd(3,1)=1,\quad e=3,\quad R=2 \\
1\colon\quad m=2,\quad d=\gcd(3,2)=1,\quad e=3,\quad R=2 \\
3\colon\quad m=1,\quad d=\gcd(3,1)=1,\quad e=3,\quad R=2 \\
\infty\colon\quad m=-4,\quad d=\gcd(3,-4)=1,\quad e=3,\quad R=2
\end{aligned}
$$

完成では各軌道に一点を加え、底の中心座標を $t^e$ とする。全体は $1$ 成分であり、番号を一つ進める置換は各成分を入れ替える。従って各成分の次数は $3$、分岐寄与は全寄与を $1$ で割った $8$ である。分岐値ごとの補正は $d$ 個の逆像を合わせた $d(e-1)$ である。

$$
\begin{aligned}
2g-2&=-2\cdot 3+8=2 \\
g&=2
\end{aligned}
$$

有限点では $x-a=t^e$、無限遠では $x=t^{-e}$ を微分する。$y$ の位数を差し引くと次の表となる。

$$
\begin{aligned}
0\colon&\quad \operatorname{ord}\mathrm dx=2,\quad \operatorname{ord}y=1,\quad\operatorname{ord}(\mathrm dx/y)=1. \\
1\colon&\quad \operatorname{ord}\mathrm dx=2,\quad \operatorname{ord}y=2,\quad\operatorname{ord}(\mathrm dx/y)=0. \\
3\colon&\quad \operatorname{ord}\mathrm dx=2,\quad \operatorname{ord}y=1,\quad\operatorname{ord}(\mathrm dx/y)=1. \\
\infty\colon&\quad \operatorname{ord}\mathrm dx=-4,\quad \operatorname{ord}y=-4,\quad\operatorname{ord}(\mathrm dx/y)=0.
\end{aligned}
$$

表にない値では $x$ が局所座標で $y$ は非消失なので、$\mathrm dx/y$ は正則かつ非零である。従ってこの微分の因子は表から完全に決まる。
:::

:::example[有限部と無限遠の Laurent データを照合する]
$\eta=R(z)\,\mathrm dz$、$R=- \frac{2}{z + 1} + \frac{2}{z - 1} + \frac{1}{\left(z - 1\right)^{3}}$ とする。有限点は分母の根から得られる。無限遠では $w=1/z$ の微分の因子 $-w^{-2}$ も含める。

$$
\begin{aligned}
R(z)&=- \frac{2}{z + 1} + \frac{2}{z - 1} + \frac{1}{\left(z - 1\right)^{3}} \\
\eta&=\left(\frac{5 w^{2} - 7 w + 4}{w^{4} - 2 w^{3} + 2 w - 1}\right)\,\mathrm dw\qquad(z=1/w)
\end{aligned}
$$

$$
\begin{aligned}
z=-1\colon&\quad m=1,\quad a_{-m}=-2,\quad \operatorname{Res}\eta=-2 \\
z=1\colon&\quad m=3,\quad a_{-m}=1,\quad \operatorname{Res}\eta=2 \\
\operatorname{Res}_\infty\eta&=0
\end{aligned}
$$

表の $m$ は極の位数、$a_{-m}$ は最上位の主要係数である。中心 $1$、半径 $1$ の円の内側で留数を合計すると $2$ となる。正向きに一周する積分は $4 i \pi$ である。

積分の途中で枝を一つ選べば係数の原始関数は次で与えられる。対数を含む項がある場合は、これを一つの大域有理型関数と解釈できない。

$$
\begin{aligned}
\int R(z)\,\mathrm dz&=2 \log{\left(z - 1 \right)} - 2 \log{\left(z + 1 \right)} - \frac{1}{2 z^{2} - 4 z + 2}+C \\
\frac{\mathrm d}{\mathrm dz}\int R(z)\,\mathrm dz&=- \frac{2}{z + 1} + \frac{2}{z - 1} + \frac{1}{\left(z - 1\right)^{3}}
\end{aligned}
$$
:::

## 線積分と留数定理

区分的滑らかな道 $\gamma$ に沿う $1$ 形式の積分は $\int \gamma^*\eta$ で定める。座標を変えても、$\gamma^*$ の連鎖律で値は同じである。閉じた形式の閉曲線積分を :term[周期]{en="period" yomi="しゅうき"} と呼ぶ。$1$ 形式が完全、すなわち $\eta=\dd f$ なら周期は零である。逆に全周期が零なら、基点からの積分が道によらず、局所原始関数と定数だけ違うため大域原始関数になる。境界をもつ滑らかな曲面片について Stokes の定理を用いると、正則形式の周期はホモロジー類だけで決まる。

:::definition[留数]{#def-cg-residue}
$p$ で有理型の微分 $\eta$ の :term[留数]{en="residue" yomi="りゅうすう"} $\operatorname{Res}_p\eta$ は、中心座標での Laurent 展開の $z^{-1}\dd z$ の係数である。これは $p$ を一回正向きに囲む小円の積分を $2\pi i$ で割った値に等しい。
:::

Laurent 級数を小円上で一様に積分すればこの等式が得られる。異なる座標で取った小円は穴あき近傍内で同じ向きに一周するため、間に Stokes の定理を使うと積分が一致する。したがって留数は座標によらない。

:::theorem[コンパクト面上の留数定理]{#thm-cg-global-residue}
コンパクト Riemann 面 $X$ 上の有理型微分 $\eta$ に対し $\sum_{p\in X}\operatorname{Res}_p\eta=0$ である。
:::

:::proof
極は閉離散集合なので有限である。各極の小座標円板を互いに重ならないように選んで除く。残った面では $\dd\eta=0$ なので、Stokes の定理により境界積分の和は零である。境界は各小円の負向きなので $0=-2\pi i\sum_p\operatorname{Res}_p\eta$ となる。
:::

:::example[球面の無限遠留数]
$\eta=(2/z+3/(z-1)^2)\dd z$ は $0$ で留数 $2$、$1$ では $0$ である。$s=1/z$ とすると $\eta=(-2/s-3/(1-s)^2)\dd s$ となり無限遠の留数は $-2$ である。
:::

:::example[引き戻しと留数の倍率]
$\eta=(3/w+2/w^2)\dd w$ を $w=z^4$ で引き戻すと $(12/z+8/z^5)\dd z$ である。一般に局所次数 $e$ の引き戻しで留数は $e$ 倍になる。標準形 $w=z^e$ を代入すると、$w^k\dd w$ は $e z^{e(k+1)-1}\dd z$ となって $k=-1$ だけが留数に寄与する。
:::

対数微分に :ref[thm-cg-global-residue] を適用すると、零点の位数の和と極の位数の和が一致する。周期と留数は異なる大域データを測る。例えばトーラスの $\dd z$ は留数をもたず、格子生成元に沿う周期は非零である。

### 有限写像に沿う微分の跡

固有な非定数正則写像 $f:X\to Y$ に対し、非分岐値の近くで逆枝 $\sigma_j$ に沿う微分を全て足して $\operatorname{Tr}_f\eta=\sum_j\sigma_j^*\eta$ とする。逆枝を取り替えても項の順序が変わるだけなので、この式は非分岐部分全体で定まる。

分岐点の近くでは、局所標準形 $w=z^e$ を取り、$\eta=\sum_k a_k z^k\dd z$ と書く。穴あき円板上の逆枝を $z=\zeta^j w^{1/e}$、$\zeta=e^{2\pi i/e}$ とすると、$z^k\dd z$ の引き戻しの和は $e\mid k+1$ の場合に $w^{(k+1)/e-1}\dd w$、それ以外では零である。実際、各項は $\zeta^{j(k+1)}w^{(k+1)/e-1}\dd w/e$ となり、単位根の和で選別される。Laurent 級数は小円環上で局所一様に収束するので、項ごとの計算で跡が有理型に延長する。微分が正則なら残る指数は全て非負であり、跡も正則である。留数については $k=-1$ だけが寄与するから、
$$
\operatorname{Res}_q(\operatorname{Tr}_f\eta)
=\sum_{p\in f^{-1}(q)}\operatorname{Res}_p\eta
$$
となる。特に、跡では逆像の点の留数をそのまま足す。引き戻しの場合には各点で局所次数が掛かる。

:::example[二枚の枝の和が零になる微分]
$w=z^2$ では $\operatorname{Tr}_f(\dd z)=0$、$\operatorname{Tr}_f(z\dd z)=\dd w$、$\operatorname{Tr}_f(\dd z/z)=\dd w/w$ である。最後の式を引き戻すと $2\dd z/z$ になる。これは一般の等式 $\operatorname{Tr}_f(f^*\beta)=(\deg f)\beta$ とも一致する。
:::

### 原始関数の特異性は留数で分かる

穴あき円板の有理型微分 $\eta=(\sum_{k\ge-m}a_k z^k)\dd z$ は、$k=-1$ 以外の項を項別積分できる。従って局所的には
$$
\eta=\dd\left(\sum_{k\ne-1}\frac{a_k}{k+1}z^{k+1}\right)
    +a_{-1}\frac{\dd z}{z}
$$
である。第一項の原始関数は穴を越えて有理型に延長する。第二項は $a_{-1}\log z$ の微分で、一周すると原始関数が $2\pi ia_{-1}$ 変わる。従って穴あき円板上の単価原始関数が存在すること、そこで有理型に延長する原始関数が存在すること、留数が零であることは同値である。

:::example[留数を持たない二位の極]
$\eta=(3/z^2+5z^4)\dd z$ は極を持つが、$-3/z+z^5$ の外微分である。$\eta+\dd z/z$ では同じ二位の極が残る一方、周期は $2\pi i$ となる。従って極の最高位数だけでは原始関数の有無を判定できない。
:::

### 積分の値を変えない変形

極を避ける面で正則微分 $\eta$ は閉じる。二つの閉曲線の差が向き付けた 2 チェインの境界なら、Stokes の定理で積分が等しい。極を横切る変形では、その点を囲む小円をチェインから除いて Stokes を使う。この小円が差に現れる回数と向きが留数の係数になる。すなわち、局所で書いた Laurent 級数の $z^{-1}$ 項だけが、変形後の積分に残る。

:::example[二つの極を持つ微分の留数表]
$\eta=(2z+1)/(z^2(z-2))\,\mathrm dz$ を球面で考える。部分分数分解は
$$
\frac{2z+1}{z^2(z-2)}=-\frac1{2z^2}-\frac5{4z}+\frac5{4(z-2)}.
$$
従って半径 $1$ の正向きの円の積分は $-5\pi i/2$、半径 $3$ では $0$ である。無限遠では係数が $O(z^{-2})$ なので微分は正則である。
:::

:::example[重極だが積分は零]
$\eta=\mathrm dz/(z-1)^3$ には原始関数 $-1/(2(z-1)^2)$ がある。極の周りの閉曲線でも積分は零である。三位の極の存在だけでは周期が決まらず、Laurent 級数の $-1$ 次係数を調べる必要がある。
:::

### 三種類の有理型微分

有理型微分のうち、極を持たないものを第一種、全留数が零で極を持つものを第二種、高々単純極で極を持つものを第三種と呼ぶ。任意の有理型微分がこの三種類のどれかに直接属するわけではない。例えば $\eta=\mathrm dz/z^2+\mathrm dz/z$ は第二種と第三種の和である。留数総和が零になるので第三種の微分は極を少なくとも二点に持つ。第二種の微分を局所的に積分すると有理型な原始関数が得られるが、大域では周期を持つことがある。

:::example[三つの単純極を指定する]
$\eta=(\mathrm dz/z+2\mathrm dz/(z-1))$ は $0,1,\infty$ でそれぞれ $1,2,-3$ の留数を持つ。$w=1/z$ で書くと
$$
\eta=-\left(\frac1w+\frac{2}{w(1-w)}\right)\mathrm dw
$$
となるので、無限遠の留数も直接読める。$\exp\int\eta=z(z-1)^2$ を非零定数倍した関数が得られ、三点での零または極の位数が留数と一致する。
:::

### 主要部と周期の照合

:::example[有限部と無限遠の Laurent データを照合する]
$\eta=R(z)\,\mathrm dz$、$R=\frac{2 z + 5}{z^{3} \left(z - 2\right)}$ とする。有限点は分母の根から得られる。無限遠では $w=1/z$ の微分の因子 $-w^{-2}$ も含める。

$$
\begin{aligned}
R(z)&=\frac{9}{8 \left(z - 2\right)} - \frac{9}{8 z} - \frac{9}{4 z^{2}} - \frac{5}{2 z^{3}} \\
\eta&=\left(\frac{5 w^{2} + 2 w}{2 w - 1}\right)\,\mathrm dw\qquad(z=1/w)
\end{aligned}
$$

$$
\begin{aligned}
z=2\colon&\quad m=1,\quad a_{-m}=\frac{9}{8},\quad \operatorname{Res}\eta=\frac{9}{8} \\
z=0\colon&\quad m=3,\quad a_{-m}=- \frac{5}{2},\quad \operatorname{Res}\eta=- \frac{9}{8} \\
\operatorname{Res}_\infty\eta&=0
\end{aligned}
$$

表の $m$ は極の位数、$a_{-m}$ は最上位の主要係数である。中心 $0$、半径 $1$ の円の内側で留数を合計すると $- \frac{9}{8}$ となる。正向きに一周する積分は $- \frac{9 i \pi}{4}$ である。

積分の途中で枝を一つ選べば係数の原始関数は次で与えられる。対数を含む項がある場合は、これを一つの大域有理型関数と解釈できない。

$$
\begin{aligned}
\int R(z)\,\mathrm dz&=- \frac{9 \log{\left(z \right)}}{8} + \frac{9 \log{\left(z - 2 \right)}}{8} + \frac{9 z + 5}{4 z^{2}}+C \\
\frac{\mathrm d}{\mathrm dz}\int R(z)\,\mathrm dz&=\frac{2 z + 5}{z^{3} \left(z - 2\right)}
\end{aligned}
$$
:::

:::example[有限部と無限遠の Laurent データを照合する]
$\eta=R(z)\,\mathrm dz$、$R=\frac{z^{2} + 1}{\left(z - 1\right)^{2}}$ とする。有限点は分母の根から得られる。無限遠では $w=1/z$ の微分の因子 $-w^{-2}$ も含める。

$$
\begin{aligned}
R(z)&=1 + \frac{2}{z - 1} + \frac{2}{\left(z - 1\right)^{2}} \\
\eta&=\left(\frac{- w^{2} - 1}{w^{4} - 2 w^{3} + w^{2}}\right)\,\mathrm dw\qquad(z=1/w)
\end{aligned}
$$

$$
\begin{aligned}
z=1\colon&\quad m=2,\quad a_{-m}=2,\quad \operatorname{Res}\eta=2 \\
\operatorname{Res}_\infty\eta&=-2
\end{aligned}
$$

表の $m$ は極の位数、$a_{-m}$ は最上位の主要係数である。中心 $0$、半径 $2$ の円の内側で留数を合計すると $2$ となる。正向きに一周する積分は $4 i \pi$ である。

積分の途中で枝を一つ選べば係数の原始関数は次で与えられる。対数を含む項がある場合は、これを一つの大域有理型関数と解釈できない。

$$
\begin{aligned}
\int R(z)\,\mathrm dz&=z + 2 \log{\left(z - 1 \right)} - \frac{2}{z - 1}+C \\
\frac{\mathrm d}{\mathrm dz}\int R(z)\,\mathrm dz&=\frac{z^{2} + 1}{\left(z - 1\right)^{2}}
\end{aligned}
$$
:::

:::example[有限部と無限遠の Laurent データを照合する]
$\eta=R(z)\,\mathrm dz$、$R=- \frac{1}{\left(z - 2\right)^{2}} + \frac{3}{z^{2}}$ とする。有限点は分母の根から得られる。無限遠では $w=1/z$ の微分の因子 $-w^{-2}$ も含める。

$$
\begin{aligned}
R(z)&=- \frac{1}{\left(z - 2\right)^{2}} + \frac{3}{z^{2}} \\
\eta&=\left(\frac{- 12 w^{2} + 12 w - 2}{4 w^{2} - 4 w + 1}\right)\,\mathrm dw\qquad(z=1/w)
\end{aligned}
$$

$$
\begin{aligned}
z=2\colon&\quad m=2,\quad a_{-m}=-1,\quad \operatorname{Res}\eta=0 \\
z=0\colon&\quad m=2,\quad a_{-m}=3,\quad \operatorname{Res}\eta=0 \\
\operatorname{Res}_\infty\eta&=0
\end{aligned}
$$

表の $m$ は極の位数、$a_{-m}$ は最上位の主要係数である。中心 $0$、半径 $1$ の円の内側で留数を合計すると $0$ となる。正向きに一周する積分は $0$ である。

積分の途中で枝を一つ選べば係数の原始関数は次で与えられる。対数を含む項がある場合は、これを一つの大域有理型関数と解釈できない。

$$
\begin{aligned}
\int R(z)\,\mathrm dz&=\frac{6 - 2 z}{z^{2} - 2 z}+C \\
\frac{\mathrm d}{\mathrm dz}\int R(z)\,\mathrm dz&=- \frac{1}{\left(z - 2\right)^{2}} + \frac{3}{z^{2}}
\end{aligned}
$$
:::

:::example[有限部と無限遠の Laurent データを照合する]
$\eta=R(z)\,\mathrm dz$、$R=\frac{3 z + 1}{z \left(z - 1\right) \left(z + 1\right)}$ とする。有限点は分母の根から得られる。無限遠では $w=1/z$ の微分の因子 $-w^{-2}$ も含める。

$$
\begin{aligned}
R(z)&=- \frac{1}{z + 1} + \frac{2}{z - 1} - \frac{1}{z} \\
\eta&=\left(\frac{w + 3}{w^{2} - 1}\right)\,\mathrm dw\qquad(z=1/w)
\end{aligned}
$$

$$
\begin{aligned}
z=-1\colon&\quad m=1,\quad a_{-m}=-1,\quad \operatorname{Res}\eta=-1 \\
z=1\colon&\quad m=1,\quad a_{-m}=2,\quad \operatorname{Res}\eta=2 \\
z=0\colon&\quad m=1,\quad a_{-m}=-1,\quad \operatorname{Res}\eta=-1 \\
\operatorname{Res}_\infty\eta&=0
\end{aligned}
$$

表の $m$ は極の位数、$a_{-m}$ は最上位の主要係数である。中心 $1$、半径 $1/2$ の円の内側で留数を合計すると $2$ となる。正向きに一周する積分は $4 i \pi$ である。

積分の途中で枝を一つ選べば係数の原始関数は次で与えられる。対数を含む項がある場合は、これを一つの大域有理型関数と解釈できない。

$$
\begin{aligned}
\int R(z)\,\mathrm dz&=2 \log{\left(z - 1 \right)} - \log{\left(z^{2} + z \right)}+C \\
\frac{\mathrm d}{\mathrm dz}\int R(z)\,\mathrm dz&=\frac{3 z + 1}{z \left(z - 1\right) \left(z + 1\right)}
\end{aligned}
$$
:::

## Dolbeault の補題

Cauchy の積分公式は正則関数を境界値から復元する。非正則関数に対しては領域内の $\bar\partial$ が加わり、その項の核が $\bar\partial$ の右逆を与える。

:::lemma[Cauchy 核の基本解]{#lem-cg-cauchy-kernel}
$\C$ 上の面積を $\dd A=\dd x\dd y$ とする。局所可積分関数 $k(z)=1/(\pi z)$ は分布として $\partial k/\partial\bar z=\delta_0$ を満たす。
:::

:::proof
試験関数 $\phi$ の台を含む円板から半径 $\epsilon$ の円板を除く。$z\ne0$ で $\dd((\phi/z)\dd z)=2i(\phi_{\bar z}/z)\dd A$ なので、Stokes の定理により
$$
-\frac1\pi\int_{|z|>\epsilon}\frac{\phi_{\bar z}}{z}\dd A
=\frac{1}{2\pi i}\int_{|z|=\epsilon}\frac{\phi(z)}z\dd z.
$$
右辺は円上の $\phi$ の平均であり $\phi(0)$ に収束する。左辺の絶対可積分性は $\int_{|z|<\epsilon}|z|^{-1}\dd A=2\pi\epsilon$ で従う。
:::

:::theorem[Dolbeault の補題]{#thm-cg-dolbeault-curve}
任意の滑らかな $(0,1)$ 形式は各点の近傍で $\bar\partial u$ と書ける。任意の滑らかな $(1,1)$ 形式は各点の近傍で $\bar\partial\eta$ と書け、$\eta$ は滑らかな $(1,0)$ 形式に取れる。
:::

:::proof
円板内で係数 $a$ に、対象の点の近くで $1$ となる滑らかな台付き関数を掛け、平面上の滑らかな台付き関数 $b$ を得る。畳み込み
$$
u(z)=(k*b)(z)=\frac1\pi\int_{\C}\frac{b(\zeta)}{z-\zeta}\dd A(\zeta)
$$
を取る。変数を $\xi=z-\zeta$ に替えると、$z$ の微分はすべて $b(z-\xi)$ に移る。コンパクトな $z$ 集合上では積分領域を固定コンパクト集合に収められ、$k$ はそこで可積分だから項別微分できる。:ref[lem-cg-cauchy-kernel] より $\partial_{\bar z}u=b$ であり、$b=a$ の近傍で $\bar\partial u=a\dd\bar z$ となる。$(1,1)$ 形式を $a\dd z\wedge\dd\bar z$ と書くなら $\eta=-u\dd z$ とすればよい。
:::

:::example[多項式の右逆]
$\bar\partial u=(2z\bar z+\bar z^2)\dd\bar z$ には $u=z\bar z^2+\bar z^3/3$ が一つの解である。さらに二つの解の差は $\bar\partial$ で消えるので正則である。従って全解は、この解に任意の正則関数を足したものである。
:::

:::example[トーラスでの大域的な障害]
$X=\C/(\Z+i\Z)$ 上の $\dd\bar z$ は局所的には $\bar\partial\bar z$ である。大域関数 $u$ で $\bar\partial u=\dd\bar z$ と書けるなら $\dd(u\dd z)=-\dd z\wedge\dd\bar z$ となる。右辺の積分は $2i$ で非零だが左辺は Stokes の定理で零である。よって局所解があっても貼り合わせられるとは限らない。
:::

### 積分表示から境界と内部の寄与を分ける

閉円板 $\bar D$ の近くで滑らかな $u$ に対し、点 $z\in D$ を中心とする小円を除いて Stokes の定理を適用する。$\zeta\ne z$ では
$\dd((u(\zeta)/(\zeta-z))\dd\zeta)=2i u_{\bar\zeta}(\zeta)\dd A(\zeta)/(\zeta-z)$
である。内側の円を縮めるとその正向き積分は $2\pi i u(z)$ になる。従って
$$
u(z)=\frac1{2\pi i}\int_{\partial D}\frac{u(\zeta)}{\zeta-z}\dd\zeta
+\frac1\pi\int_D\frac{u_{\bar\zeta}(\zeta)}{z-\zeta}\dd A(\zeta).
$$
これは :ref[prop-cpx-cauchy-pompeiu] の式である。第一項は $z$ について正則であり、第二項だけで $\bar\partial$ の値を再現する。正則関数の場合、第二項が零になって Cauchy の公式へ戻る。

:::example[境界積分が零の例]
単位円板で $u(z)=\bar z$ とする。境界では $\bar\zeta=1/\zeta$ である。$z\ne0$ なら $1/(\zeta(\zeta-z))$ の極 $0,z$ の留数は $-1/z,1/z$ で相殺し、$z=0$ でも $\zeta^{-2}$ の留数は零である。従って境界項は零である。内部項は $u_{\bar\zeta}=1$ を使って $\bar z$ を再現する。この例では解の正則部分を境界項が含まない。
:::

### 円板で境界の Fourier 成分を読む

Cauchy–Pompeiu の境界項を単位円で計算してみる。$u(z)=z^a\bar z^b$、$a,b$ は非負整数とする。境界上で $u(\zeta)=\zeta^{a-b}$、核 $1/(\zeta-z)$ を $\sum_{k\ge0}z^k\zeta^{-k-1}$ に展開できる。固定した $|z|<1$ に対してこれは境界上で一様収束するため、積分後に $k=a-b$ の項だけが残る。従って境界項は $a\ge b$ なら $z^{a-b}$、$a<b$ なら $0$ である。残りが内部項となり、$\bar\partial u=b z^a\bar z^{b-1}\mathrm d\bar z$ を再現する。

:::example[同じ境界値を持つ二種類の延長]
$u=z^3\bar z^2$ の境界値は $z$ に一致する。境界積分は $z$、内部積分は $z^3\bar z^2-z$ である。内部の解は境界で零になるが、$\bar\partial$ の解がいつも境界零を満たせるわけではない。$u=\bar z^2$ の場合、境界の Fourier 指数は $-2$ なので境界項は零となり、内部項の境界値は $\bar z^2$ のままである。
:::

:::example[係数を積分して型の符号を調べる]
$\alpha=(6z^2\bar z+4\bar z^3)\mathrm dz\wedge\mathrm d\bar z$ とする。$u=3z^2\bar z^2+\bar z^4$ について $u_{\bar z}=6z^2\bar z+4\bar z^3$ である。$\bar\partial(u\mathrm dz)=u_{\bar z}\mathrm d\bar z\wedge\mathrm dz=-\alpha$ なので、求める原始形式は $-u\mathrm dz$ である。関数の $\bar\partial$ の式をそのまま係数に適用するとこの符号を落としやすい。
:::

### 大域解が存在する場合の一意性

連結コンパクトな面上で $\bar\partial u=\alpha$ の大域解が二つあれば、差は正則であり最大値原理から定数である。従って平均値零の条件を加えれば解は一意である。非コンパクト面では差が非定数になれる。例えば $\C$ 上で $\alpha=\mathrm d\bar z$ なら $\bar z$ と $\bar z+z^2$ は異なる解である。大域の構成で一意性の条件を用いるときは、空間のコンパクト性と境界条件を確認する必要がある。

### 周期条件が局所解の貼り合わせを左右する

:::example[平均と非零 Fourier モードを分離する]
$X=\C/(\Z+3i\Z)$ 上で $\alpha=a\,\mathrm d\bar z$ とし、周期座標を $z=x+iy$ とする。係数と微分の規約は次である。

$$
\begin{aligned}
a(x,y)&=3 e^{2 i \pi \left(x + \frac{2 y}{3}\right)} + 2 + e^{- \frac{2 i \pi y}{3}} \\
\partial_{\bar z}&=\tfrac12(\partial_x+i\partial_y)
\end{aligned}
$$

周期積分では全微分の平均が零になる。したがって定数モードは $\bar\partial$ の方程式の障害である。他方、指数関数 $e^{2\pi i(mx+ny/L)}$ の二つの微分で $2\pi im,2\pi in/L$ が掛かるので、固有値は $\pi(im-n/L)$ になる。

$$
\begin{aligned}
(1,2) &\colon\quad \lambda=\pi \left(- \frac{2}{3} + i\right),\quad a=3,\quad a/\lambda=\frac{9 \left(-2 - 3 i\right)}{13 \pi} \\
(0,-1) &\colon\quad \lambda=\frac{\pi}{3},\quad a=1,\quad a/\lambda=\frac{3}{\pi} \\
a_{00}&=2
\end{aligned}
$$

定数モードを除いて各係数を固有値で割る。これで得られる関数は平均零で、微分を取り直すと次を満たす。

$$
\begin{aligned}
u&=\frac{27 \left(- \frac{2}{3} - i\right) e^{2 i \pi \left(x + \frac{2 y}{3}\right)}}{13 \pi} + \frac{3 e^{- \frac{2 i \pi y}{3}}}{\pi} \\
u_{\bar z}&=a-2
\end{aligned}
$$

全てのモードの周期が元のトーラスと一致するため $u$ は大域的である。平均が二であるから、局所原始関数をそのまま全体に貼ることはできない。
:::

:::example[平均と非零 Fourier モードを分離する]
$X=\C/(\Z+2i\Z)$ 上で $\alpha=a\,\mathrm d\bar z$ とし、周期座標を $z=x+iy$ とする。係数と微分の規約は次である。

$$
\begin{aligned}
a(x,y)&=- 2 e^{2 i \pi y} + e^{2 i \pi \left(- x - \frac{y}{2}\right)} + 2 e^{2 i \pi \left(x - \frac{y}{2}\right)} \\
\partial_{\bar z}&=\tfrac12(\partial_x+i\partial_y)
\end{aligned}
$$

周期積分では全微分の平均が零になる。したがって定数モードは $\bar\partial$ の方程式の障害である。他方、指数関数 $e^{2\pi i(mx+ny/L)}$ の二つの微分で $2\pi im,2\pi in/L$ が掛かるので、固有値は $\pi(im-n/L)$ になる。

$$
\begin{aligned}
(1,-1) &\colon\quad \lambda=\pi \left(\frac{1}{2} + i\right),\quad a=2,\quad a/\lambda=\frac{4 \left(1 - 2 i\right)}{5 \pi} \\
(-1,-1) &\colon\quad \lambda=\pi \left(\frac{1}{2} - i\right),\quad a=1,\quad a/\lambda=\frac{2 \left(1 + 2 i\right)}{5 \pi} \\
(0,2) &\colon\quad \lambda=- \pi,\quad a=-2,\quad a/\lambda=\frac{2}{\pi} \\
a_{00}&=0
\end{aligned}
$$

定数モードを除いて各係数を固有値で割る。これで得られる関数は平均零で、微分を取り直すと次を満たす。

$$
\begin{aligned}
u&=\frac{2 e^{2 i \pi y}}{\pi} + \frac{4 \left(\frac{1}{2} + i\right) e^{2 i \pi \left(- x - \frac{y}{2}\right)}}{5 \pi} + \frac{8 \left(\frac{1}{2} - i\right) e^{2 i \pi \left(x - \frac{y}{2}\right)}}{5 \pi} \\
u_{\bar z}&=a-0
\end{aligned}
$$

全てのモードの周期が元のトーラスと一致するため $u$ は大域的である。平均が零の場合には、微分を取り直した係数は元の右辺全体である。
:::

### 多項式上の局所右逆の違いを見つける

$u(z,\bar z)=\sum_{a,b}u_{ab}z^a\bar z^b$ を多項式と見て、$\bar z$ についての微分を取る。$\frac{\partial}{\partial\bar z}(z^a\bar z^{b+1}/(b+1))=z^a\bar z^b$ なので、多項式入力には積分核を使わずに原始関数を求められる。核は $\bar z$ の正次数を持たない多項式である。円板上の滑らかな解を許すと、この核は正則関数の全体に広がる。

:::example[次数と一点の規格化を比べる]
$\alpha=(4z^2+3z\bar z^2+2\bar z^3)\mathrm d\bar z$ なら、一つの原始関数は $u=4z^2\bar z+z\bar z^3+\bar z^4/2$ である。$u(0)=0$ が成り立つが、この規格化だけでは $u+z$ も解になる。単位円周上の非負 Fourier 成分を除く規格化を採れば、$4z^2\bar z$ の境界の指数だけが $1$ で、他の二項の指数は負である。従って $u-4z$ となる。この差は確かに正則である。
:::

## 調和関数と Weyl の補題

$u$ が :term[調和関数]{en="harmonic function" yomi="ちょうわかんすう"} であるとは、座標で $\Delta u=u_{xx}+u_{yy}=0$ となることをいう。$\Delta=4\partial_z\partial_{\bar z}$ である。正則座標の変更 $w=w(z)$ に対して $\Delta_z(u\circ w)=|w'|^2(\Delta_wu)\circ w$ となる。これは連鎖律で、$w$ の $\bar z$ 微分が消えることから従う。実数値の調和関数では $2\partial u$ が正則微分で、その局所原始関数の実部が $u$ と定数だけ異なる。

:::example[調和共役の周期]
$\C^\times$ で $u=\log|z|$ は調和で $2\partial u=\dd z/z$ となる。単位円の周期は $2\pi i$ なので、$u$ を実部に持つ大域単価な正則関数は存在しない。
:::

:::theorem[Weyl の補題]{#thm-cg-weyl}
平面の開集合 $U$ で $u\in L^1_{\mathrm{loc}}(U)$ が弱い意味で $\Delta u=0$、すなわち全試験関数 $\phi$ に対し $\int u\Delta\phi\,\dd A=0$ なら、$u$ は滑らかな調和関数とほとんど至る所で等しい。さらに $\partial_{\bar z}u=0$ が分布の意味で成り立つなら、その代表は正則である。
:::

:::proof
非負で回転対称な滑らかな単位積分の関数 $\rho$ を単位円板に台を持つように選び、$\rho_\epsilon(z)=\epsilon^{-2}\rho(z/\epsilon)$ とする。$u_\epsilon=u*\rho_\epsilon$ は境界から $\epsilon$ 以上離れた部分で滑らかであり、弱い調和性を試験関数 $\rho_\epsilon(z-\cdot)$ に使うと $\Delta u_\epsilon=0$ となる。また $u_\epsilon\to u$ は局所 $L^1$ 収束する。これは平行移動の $L^1$ 連続性と $\rho_\epsilon$ の台の縮小から従う。

任意の比較的小さい閉円板 $B$ を取り、その近傍の固定半径 $r$ の円が $U$ に入るように選ぶ。調和関数の円周平均を重みにより積分すれば $u_\epsilon=u_\epsilon*\rho_r$ が $B$ 上で成立する。右辺は局所 $L^1$ 収束と滑らかな核の全導関数の有界性により、$B$ 上で全導関数について $u*\rho_r$ に一様収束する。よって $u*\rho_r$ は滑らかで調和であり、$L^1$ 極限を比較して $u$ とほとんど至る所で一致する。円板を替えても滑らかな代表は連続性で一致する。最後に $\partial_{\bar z}u=0$ なら $\Delta u=4\partial_z\partial_{\bar z}u=0$ なので同じ代表が得られ、その Cauchy–Riemann 方程式を使えば正則である。
:::

:::example[点に特異性が集中する場合]
$\log|z|$ は局所可積分で穴あき円板では調和である。しかし半径 $r$ の円で外向き微分の積分は $(1/r)2\pi r=2\pi$ である。試験関数に対する Green の公式を小円を除いて適用すると $\Delta\log|z|=2\pi\delta_0$ となる。小円上の $\log r$ と試験関数の微分の積分は $O(r|\log r|)$ で消え、$1/r$ を掛けた試験関数の平均だけが残る。従って原点を含む円板では弱い調和性が成り立たない。
:::

次章では、局所解が大域的につながらない差を $1$ 次の層コホモロジーで記録する。

### 可除特異点と弱形式を区別する

穴あき円板上の実調和関数 $u$ が原点で有界なら、穴を越えて調和に延びる。まず半径 $R$ の円周上で $u$ と一致する Poisson 解 $h$ を取る。$u-h$ は外周で零で、内半径 $\epsilon$ の円上では一様定数 $M$ で抑えられる。比較原理で
$$
|u(z)-h(z)|\le M\frac{\log(|z|/R)}{\log(\epsilon/R)}
\qquad(\epsilon<|z|<R)
$$
を得る。固定した $z\ne0$ で $\epsilon\to0$ とすると右辺は零なので $u=h$ となる。有界性を失うと $\log|z|$ のような対数項が残り得る。局所可積分性だけでは $\Delta u=0$ を点を含む弱形式に延長できないことは、原点の質量 $2\pi$ の計算から分かる。

:::example[ノルムの差が現れる例]
半径 $1$ の円板で $\log|z|$ は $L^2$ に属する。半径積分は $2\pi\int_0^1r(\log r)^2\dd r=\pi/2$ で有限である。一方、勾配の絶対値は $1/r$ なのでその二乗積分は $2\pi\int_0^1\dd r/r=\infty$ になる。従って関数の可積分性と弱微分の正則性を分けて確認する。
:::

### Dirichlet エネルギーの座標不変性

実数値の滑らかな関数 $u$ に対し、座標円板上のエネルギー密度は $(u_x^2+u_y^2)\mathrm dx\wedge\mathrm dy$ である。別の正則座標に移ると勾配のノルムの二乗に $|w'|^2$ が掛かり、面積要素の変換で同じ因子が現れるので密度は座標によらない。$u+tv$ を代入して $t=0$ で微分すると $2\int\nabla u\cdot\nabla v\,\mathrm dA=-2\int (\Delta u)v\,\mathrm dA$ となる。境界に達しない台を持つ $v$ のすべてに対してこれが零なら弱い意味で調和である。Weyl の補題を適用できる局所可積分性のもとでは、古典的に調和な代表が得られる。

:::example[調和多項式のエネルギー]
単位円板で $u=\Re(z^3)=x^3-3xy^2$ を取る。$u_x=3(x^2-y^2)$、$u_y=-6xy$ なので $|\nabla u|^2=9(x^2+y^2)^2$ となる。極座標で積分すると $2\pi\cdot9\int_0^1r^5\mathrm dr=3\pi$ である。この値は境界積分 $\int_{\partial D}u\,\partial_\nu u\,\mathrm ds$ からも、$3\int_0^{2\pi}\cos^2(3\theta)\mathrm d\theta=3\pi$ と得られる。
:::

:::example[対数のエネルギーが発散する速さ]
$u=\log|z|$ に対して、$r_0<|z|<1$ のエネルギーは $2\pi\int_{r_0}^1r^{-1}\mathrm dr=-2\pi\log r_0$ である。$r_0\downarrow0$ で発散する。一方 $u$ 自体の絶対値の局所積分は $\int_0^1r|\log r|\mathrm dr=1/4$ なので有限である。Weyl の仮定が失敗する箇所は局所可積分性ではなく、試験関数に対する Laplace 方程式である。
:::

### 調和関数の単価性を周期で調べる

実調和関数 $u$ について $2\partial u$ の実部は $\mathrm du$ である。従ってこの微分の閉曲線上の周期は純虚数である。全部が零なら大域の正則関数 $F$ がその原始関数になり、定数を調整して $u=\Re F$ とできる。周期が $2\pi i\Z$ に入る場合は $e^{\int 2\partial u}$ が大域正則関数となるが、積分そのものは単価である必要がない。

:::example[円環の回転対称な調和関数]
$u(z)=U(r)$ が $a<r<b$ で調和なら、$\Delta u=U''(r)+U'(r)/r=0$ なので $(rU')'=0$ である。従って $u=A\log r+B$ となる。例えば $1<r<e^2$ で境界値 $0,3$ を指定すると $A=3/2,B=0$ である。対応する周期は $3\pi i$ なので大域の調和共役は存在しないが、$e^{2\int 2\partial u}$ は定数倍を除いて $z^3$ である。
:::

### 領域の形と Dirichlet エネルギー

:::example[円環の角モード 1 を解く]
内半径 $1$、外半径 $2$ の円環で、内周の値を $0\cos{\left(\theta \right)}$、外周の値を $2\cos{\left(\theta \right)}$ とする。角方向のモードを固定して $u(r,\theta)=v(r)\cos{\left(\theta \right)}$ と置く。円環には内外の境界があるため、正と負の半径のべきを許す。

$$
\begin{aligned}
0&=v''+r^{-1}v'-1^2r^{-2}v \\
v(r)&=c_1r^1+c_2r^{-1} \\
\left[\begin{matrix}1 & 1\\2 & \frac{1}{2}\end{matrix}\right]\begin{pmatrix}c_1\\c_2\end{pmatrix}&=\left[\begin{matrix}0\\2\end{matrix}\right] \\
\begin{pmatrix}c_1\\c_2\end{pmatrix}&=\left[\begin{matrix}\frac{4}{3}\\- \frac{4}{3}\end{matrix}\right] \\
v(r)&=\frac{4 r}{3} - \frac{4}{3 r}
\end{aligned}
$$

内外の半径を代入すると $v(1)=0$、$v(2)=2$ となる。角方向の二階微分が $-1^2$ を掛けるため、半径方向の方程式との和で Laplacian が零になる。二解の差は零境界値の調和関数なので、最大値原理で零になる。

二つの半径の幾何平均 $\sqrt{2}$ での振幅は $\frac{2 \sqrt{2}}{3}$ である。Dirichlet エネルギーを境界積分で求める。内周の外向き微分は $-\partial_r$ であり、角モードの二乗の積分は $\pi$ である。

$$
\begin{aligned}
\int |\nabla u|^2\,\mathrm dA&=\pi\bigl[Rv(R)v'(R)-av(a)v'(a)\bigr] \\
&=\frac{20 \pi}{3}
\end{aligned}
$$

同じ量は正の積分 $\pi\int_a^R(r(v')^2+k^2v^2/r)\,\mathrm dr$ でも得られる。境界と内部で値を照合すると、法線の符号が正しいことを確かめられる。
:::

:::example[円環の角モード 2 を解く]
内半径 $1$、外半径 $3$ の円環で、内周の値を $1\sin{\left(2 \theta \right)}$、外周の値を $0\sin{\left(2 \theta \right)}$ とする。角方向のモードを固定して $u(r,\theta)=v(r)\sin{\left(2 \theta \right)}$ と置く。円環には内外の境界があるため、正と負の半径のべきを許す。

$$
\begin{aligned}
0&=v''+r^{-1}v'-2^2r^{-2}v \\
v(r)&=c_1r^2+c_2r^{-2} \\
\left[\begin{matrix}1 & 1\\9 & \frac{1}{9}\end{matrix}\right]\begin{pmatrix}c_1\\c_2\end{pmatrix}&=\left[\begin{matrix}1\\0\end{matrix}\right] \\
\begin{pmatrix}c_1\\c_2\end{pmatrix}&=\left[\begin{matrix}- \frac{1}{80}\\\frac{81}{80}\end{matrix}\right] \\
v(r)&=- \frac{r^{2}}{80} + \frac{81}{80 r^{2}}
\end{aligned}
$$

内外の半径を代入すると $v(1)=1$、$v(3)=0$ となる。角方向の二階微分が $-2^2$ を掛けるため、半径方向の方程式との和で Laplacian が零になる。二解の差は零境界値の調和関数なので、最大値原理で零になる。

二つの半径の幾何平均 $\sqrt{3}$ での振幅は $\frac{3}{10}$ である。Dirichlet エネルギーを境界積分で求める。内周の外向き微分は $-\partial_r$ であり、角モードの二乗の積分は $\pi$ である。

$$
\begin{aligned}
\int |\nabla u|^2\,\mathrm dA&=\pi\bigl[Rv(R)v'(R)-av(a)v'(a)\bigr] \\
&=\frac{41 \pi}{20}
\end{aligned}
$$

同じ量は正の積分 $\pi\int_a^R(r(v')^2+k^2v^2/r)\,\mathrm dr$ でも得られる。境界と内部で値を照合すると、法線の符号が正しいことを確かめられる。
:::

:::example[円環の角モード 3 を解く]
内半径 $2$、外半径 $3$ の円環で、内周の値を $1\cos{\left(3 \theta \right)}$、外周の値を $1\cos{\left(3 \theta \right)}$ とする。角方向のモードを固定して $u(r,\theta)=v(r)\cos{\left(3 \theta \right)}$ と置く。円環には内外の境界があるため、正と負の半径のべきを許す。

$$
\begin{aligned}
0&=v''+r^{-1}v'-3^2r^{-2}v \\
v(r)&=c_1r^3+c_2r^{-3} \\
\left[\begin{matrix}8 & \frac{1}{8}\\27 & \frac{1}{27}\end{matrix}\right]\begin{pmatrix}c_1\\c_2\end{pmatrix}&=\left[\begin{matrix}1\\1\end{matrix}\right] \\
\begin{pmatrix}c_1\\c_2\end{pmatrix}&=\left[\begin{matrix}\frac{1}{35}\\\frac{216}{35}\end{matrix}\right] \\
v(r)&=\frac{r^{3}}{35} + \frac{216}{35 r^{3}}
\end{aligned}
$$

内外の半径を代入すると $v(2)=1$、$v(3)=1$ となる。角方向の二階微分が $-3^2$ を掛けるため、半径方向の方程式との和で Laplacian が零になる。二解の差は零境界値の調和関数なので、最大値原理で零になる。

二つの半径の幾何平均 $\sqrt{6}$ での振幅は $\frac{12 \sqrt{6}}{35}$ である。Dirichlet エネルギーを境界積分で求める。内周の外向き微分は $-\partial_r$ であり、角モードの二乗の積分は $\pi$ である。

$$
\begin{aligned}
\int |\nabla u|^2\,\mathrm dA&=\pi\bigl[Rv(R)v'(R)-av(a)v'(a)\bigr] \\
&=\frac{114 \pi}{35}
\end{aligned}
$$

同じ量は正の積分 $\pi\int_a^R(r(v')^2+k^2v^2/r)\,\mathrm dr$ でも得られる。境界と内部で値を照合すると、法線の符号が正しいことを確かめられる。
:::

## 演習問題

::::exercise{level=basic}
$\alpha=z^2\bar z\dd z+z\bar z^2\dd\bar z$ の外微分を求めよ。
:::solution
$\dd(a\dd z+b\dd\bar z)=(b_z-a_{\bar z})\dd z\wedge\dd\bar z$ を使う。$b_z=\bar z^2$、$a_{\bar z}=z^2$ なので $(\bar z^2-z^2)\dd z\wedge\dd\bar z$ である。
:::
::::

::::exercise{level=basic}
球面で $\eta=(1/z+2/(z-2))\dd z$ の留数をすべて求めよ。
:::solution
$z=0,2$ で他方の項が正則なので留数は順に $1,2$ である。$s=1/z$ では $\eta=-(1+2/(1-2s))\dd s/s$ だから $s=0$ の留数は $-3$ である。
:::
::::

::::exercise{level=basic}
円板で $\bar\partial u=z^2\bar z^3\dd\bar z$ を解け。
:::solution
$u=z^2\bar z^4/4$ を $\bar z$ で微分すると右辺の係数になる。全解はこれと任意の正則関数の和であり、差が正則であることは Cauchy–Riemann 方程式で従う。
:::
::::

::::exercise{level=basic}
$u=x^3-3xy^2$ が調和であることを確認し、$u$ の調和共役を求めよ。
:::solution
$u_{xx}=6x,u_{yy}=-6x$ なので調和である。$(x+iy)^3$ の実部が $u$、虚部が $3x^2y-y^3$ だから、調和共役は虚部に任意の実定数を加えたものである。
:::
::::

::::exercise{level=advanced}
トーラス $X=\C/(\Z+2i\Z)$ 上の閉形式 $\dd z+c\dd\bar z$ が完全になる $c\in\C$ はあるか。
:::solution
格子ベクトル $1,2i$ に沿う周期はそれぞれ $1+c,2i(1-c)$ である。完全なら両方が零だが、前者は $c=-1$、後者は $c=1$ を要求する。従ってそのような $c$ はない。
:::
::::

::::exercise{level=advanced}
コンパクト Riemann 面 $X$ の有理型関数 $f$ が零点と極をもたないなら、定数であることを二通り示せ。
:::solution
極がないから $f:X\to\C$ は正則であり、コンパクト性と最大値原理から定数となる。別法として $u=\log|f|$ を取る。零点と極がないので局所対数が存在し $u$ は大域調和関数になる。コンパクト性で最大点があり、平均値性質で最大値を取る集合は開となる。連結性で $u$ は定数となり、開写像定理から $f$ も定数である。
:::
::::

::::exercise{level=advanced}
$\C^\times$ 上の $a z^{-2}\dd z+b z^{-1}\dd z$ が正則原始関数を持つ条件を求めよ。
:::solution
単位円積分は $2\pi i b$ なので $b=0$ が必要である。$b=0$ なら $\dd(-a/z)=az^{-2}\dd z$ となるので十分である。
:::
::::

::::exercise{level=advanced}
$u\in L^1_{\mathrm{loc}}(U)$ が弱い意味で $\partial_{\bar z}u=v$ を満たし、$v$ が滑らかであるとする。$u$ は滑らかな代表をもつことを示せ。
:::solution
各点の近くで Dolbeault の補題を使い、滑らかな $w$ を $\partial_{\bar z}w=v$ となるように選ぶ。$u-w$ の $\bar z$ 微分は分布として零なので :ref[thm-cg-weyl] により正則な代表をもつ。$w$ を戻すと $u$ の滑らかな代表を得る。重なる領域上の二つの滑らかな代表はほとんど至る所で一致し、連続性から全点で一致する。
:::
::::

## 確認テスト

::::quiz{answer="-5" tol="1e-6"}
$\eta=(3/z+2/(z-1))\dd z$ の無限遠での留数を求めよ。
:::explanation
有限点の留数は $3,2$ であり、球面での留数の総和は零である。よって無限遠では $-5$ である。
:::
::::

::::quiz
$u$ を滑らかな関数とする。型分解の式として正しいものを選べ。

- [ ] $\bar\partial u=u_z\dd z$。
- [x] $\bar\partial u=u_{\bar z}\dd\bar z$。
- [ ] $\bar\partial u=u_{\bar z}\dd z$。

:::explanation
連鎖律で $\dd u=u_z\dd z+u_{\bar z}\dd\bar z$ となる。$\bar\partial$ は二番目の型の成分を取るので、$u_{\bar z}\dd\bar z$ である。
:::
::::

::::quiz{answer="4" tol="1e-6"}
$w=z^2$ による $2\dd w/w$ の引き戻しの原点での留数を求めよ。
:::explanation
$\dd w=2z\dd z$ を代入すると $4\dd z/z$ であり、留数は $4$ である。
:::
::::

::::quiz
Dolbeault の補題から直ちに得られる結論を選べ。

- [x] 各点の近くで滑らかな $(0,1)$ 形式は $\bar\partial u$ と書ける。
- [ ] コンパクト面上の全ての $(0,1)$ 形式は大域関数の $\bar\partial$ である。
- [ ] 任意の局所解は一意である。

:::explanation
Cauchy 核で係数に合わせた局所解を作れる。解に正則関数を加えても $\bar\partial$ は同じである。平方トーラス上の $\dd\bar z$ の例では局所解が大域関数を与えない。
:::
::::

::::quiz{answer="6" tol="1e-6"}
$z=t^3$ による $2\,\mathrm dz/z$ の引き戻しの $t=0$ での留数を求めよ。

:::explanation
$\mathrm dz=3t^2\mathrm dt$ を代入すると、$2\cdot3t^2\mathrm dt/t^3=6\mathrm dt/t$ である。$t^{-1}$ の係数 $6$ が留数である。
:::
::::

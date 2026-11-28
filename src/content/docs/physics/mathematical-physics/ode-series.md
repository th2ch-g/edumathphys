---
title: "常微分方程式の級数解法"
description: "二階線形方程式の局所解を、正則点、確定特異点、不確定特異点に分けて構成する。"
notation:
  - symbol: 'P,Q'
    meaning: "正規形の一階項と零階項の係数"
  - symbol: 'I(r)'
    meaning: "確定特異点における決定多項式"
  - symbol: 'W(y_1,y_2)'
    meaning: "二解の Wronskian"
---

## この章の目標

- 正則点で初期値に対応する冪級数を構成する。
- Frobenius 解の指数と共鳴時の対数項を求める。
- 無限遠を含む特異点を分類する。
- Wronskian を使って第 2 解を作る。

[級数解法](/math/ode/series-solutions/) をもとに、係数や独立変数を複素数にしても使える方法を導く。実数の例では、分数冪や対数は特に断らない限り $x>0$ で扱う。

## 正則点と冪級数解

### 初期値から係数を決める

正規形 $y''+P(z)y'+Q(z)y=0$ で、係数 $P,Q$ が $z_0$ の近くで正則なら $z_0$ を :term[正則点]{en="ordinary point" yomi="せいそくてん"} という。$x=z-z_0$、$P=\sum p_nx^n$、$Q=\sum q_nx^n$、$y=\sum a_nx^n$ と置く。係数を比較すると
$$
(n+2)(n+1)a_{n+2}=-\sum_{j=0}^n[(j+1)a_{j+1}p_{n-j}+a_jq_{n-j}].
$$
従って $a_0=y(z_0),a_1=y'(z_0)$ が決まれば全係数が定まる。

収束は一階系 $Y'=\mathcal A(x)Y$、$Y=(y,y')^\top$ と直して示せる。係数が正則な任意の閉円板 $|x|\le R$ で $\|\mathcal A\|\le M$ とする。積分作用素 $Tf(x)=\int_0^x\mathcal A(\zeta)f(\zeta)\dd\zeta$ を直線分で定める。$T$ は正則関数を正則関数に移し、帰納的に $\|T^jY_0(x)\|\le(M|x|)^j\|Y_0\|/j!$ である。従って $\sum_{j\ge0}T^jY_0$ は一様収束し、積分方程式から微分方程式を満たす。二解の差に同じ反復評価を適用すると一意性を得る。Taylor 級数は係数から得た級数と一致する。

:::example[空間に比例する係数]
$y''-xy=0$ では $a_2=0$、$(n+2)(n+1)a_{n+2}=a_{n-1}$（$n\ge1$）である。$(a_0,a_1)=(1,0)$ なら $1+x^3/6+x^6/180+\cdots$、$(0,1)$ なら $x+x^4/12+x^7/504+\cdots$ となる。係数は整関数なので両級数は全平面で収束する。
:::

### 局所収束の範囲を見積もる

中心を $x=0$ に移し、正規形を $y''+P(x)y'+Q(x)y=R(x)$ とする。係数の冪級数が半径 $d$ の円板で正則なら、その円板内で解を作れる。これを見るため $Y=(y,y')^{\mathsf T}$ とし、一次系 $Y'=B(x)Y+g(x)$ にする。半径 $\rho<d$ の閉円板では $B,g$ が有界である。放射状積分による逐次近似を使うと、差の第 $j$ 項は $\norm{Y_1-Y_0}_\infty(M\rho)^{j-1}/(j-1)!$ 以下になる。これは $j-1$ 重積分の領域が単体であるためである。したがって級数は一様に収束する。各項は正則であり、積分方程式を満たす極限も正則である。差の反復に同じ階乗評価を適用すると、同じ初期値を持つ二解の差は零になる。Taylor 係数の逐次決定は、この正則解の係数を計算している。

係数の最も近い特異点までの距離は、解の収束半径の下限である。個別の解では級数がその距離より遠くまで収束することがある。例えば後述の逆正接の方程式では、同じ正規形から生じる定数解と逆正接で収束半径が異なる。

:::example[Airy 型の第一系列]
$y''=xy$ の $y(0)=1,y'(0)=0$ の解では $a_0=1$、$a_1=a_2=0$ である。漸化式 $a_{j+3}=a_j/[(j+3)(j+2)]$ から、3 の倍数の次数にだけ非零係数が残る。最初の三つは $a_0=1$、$a_3=1/6$、$a_6=1/180$ である。次の非零項は 9 次なので
$$
y(x)=1+\frac{x^3}{6}+\frac{x^6}{180}+O(x^9)
$$
となる。以下の例では、この剰余と打ち切った多項式を方程式へ代入したときの残差を区別する。
:::

:::example[Airy 型の残差]
この有限和を二階方程式へ戻すと、相殺後の残差の次数が分かる。

係数を二階、一階、零階の順に並べる。

$$
\left( 1, \  0, \  - x\right)
$$

$$
\begin{aligned}u&=\frac{x^{6}}{180} + \frac{x^{3}}{6} + 1\\u^{\prime}&=\frac{x^{5}}{30} + \frac{x^{2}}{2}\\u^{\prime\prime}&=x \left(\frac{x^{3}}{6} + 1\right)\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=x \left(\frac{x^{3}}{6} + 1\right)+\left(0\right)+\left(- x \left(\frac{x^{6}}{180} + \frac{x^{3}}{6} + 1\right)\right)\\&=- \frac{x^{7}}{180}\end{aligned}
$$

$x$ が $0$ のところで値を評価する。

$$
\begin{aligned}u&=1\\ u^\prime&=0\\ Lu&=0\end{aligned}
$$
:::

:::example[初期傾斜を選んだ解]
同じ方程式で $y(0)=0,y^\prime(0)=1$ とし、非零係数の並びを変える。

係数を二階、一階、零階の順に並べる。

$$
\left( 1, \  0, \  - x\right)
$$

$$
\begin{aligned}u&=\frac{x^{7}}{504} + \frac{x^{4}}{12} + x\\u^{\prime}&=\frac{x^{6}}{72} + \frac{x^{3}}{3} + 1\\u^{\prime\prime}&=x^{2} \left(\frac{x^{3}}{12} + 1\right)\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=x^{2} \left(\frac{x^{3}}{12} + 1\right)+\left(0\right)+\left(- x \left(\frac{x^{7}}{504} + \frac{x^{4}}{12} + x\right)\right)\\&=- \frac{x^{8}}{504}\end{aligned}
$$

$x$ が $0$ のところで値を評価する。

$$
\begin{aligned}u&=0\\ u^\prime&=1\\ Lu&=0\end{aligned}
$$
:::

:::example[多項式で止まる場合]
Legendre 型の係数で終端する最高次が異なる例を取る。

係数を二階、一階、零階の順に並べる。

$$
\left( 1 - x^{2}, \  - 2 x, \  6\right)
$$

$$
\begin{aligned}u&=\frac{3 x^{2}}{2} - \frac{1}{2}\\u^{\prime}&=3 x\\u^{\prime\prime}&=3\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=3 - 3 x^{2}+\left(- 6 x^{2}\right)+\left(9 x^{2} - 3\right)\\&=0\end{aligned}
$$

$x$ が $0$ のところで値を評価する。

$$
\begin{aligned}u&=- \frac{1}{2}\\ u^\prime&=0\\ Lu&=0\end{aligned}
$$
:::

### 偶奇の保存と係数列

係数が $P(-x)=-P(x)$、$Q(-x)=Q(x)$ を満たす斉次方程式では、$v(x)=u(-x)$ も解になる。初期傾斜が零の場合には $v(0)=u(0)$ と $v'(0)=u'(0)=0$ がそろうので、初期値問題の一意性から $v=u$ となる。初期値が零のときは $-v$ と $u$ の二条件がそろうため、解は奇関数である。したがって係数の計算を半分に減らせる。ただし中心を移動すると $P,Q$ の偶奇は変わり得るので、移した中心での対称条件を使う。

例えば $(1+x^2)u''+2xu'=0$ では $(1+x^2)u'$ が一定である。初期値 $(0,1)$ で $u=\arctan x$ を得る。級数を用いる場合は $a_1=1$ とし、$(n+2)(n+1)a_{n+2}=-n(n+1)a_n$ である。従って偶数係数はすべて零、奇数係数は $a_{2j+1}=(-1)^j/(2j+1)$ となる。$|x|<1$ でこの級数を一度微分すると幾何級数 $1-x^2+x^4-\cdots=(1+x^2)^{-1}$ になる。$x=i,-i$ で導関数が極を持つため、原点まわりの収束半径は正確に $1$ である。他方、初期値 $(1,0)$ なら定数解になり、同じ係数の方程式でも解自体の収束半径は無限である。
:::example[実軸上の初期傾斜と複素特異点]
正則点の中心を原点に取り、上で分離した奇関数の級数を具体化する。

まず展開する関数と中心を定める。中心で微分してから階乗で割ると各係数が得られる。

$$
f(x)=\operatorname{atan}{\left(x \right)},\qquad x_0=0
$$

計算した係数を、定数項から順に並べる。

$$
\begin{aligned}
\frac{f^{(0)}(0)}{0!} &= 0 \\
\frac{f^{(1)}(0)}{1!} &= 1 \\
\frac{f^{(2)}(0)}{2!} &= 0 \\
\frac{f^{(3)}(0)}{3!} &= - \frac{1}{3} \\
\frac{f^{(4)}(0)}{4!} &= 0 \\
\frac{f^{(5)}(0)}{5!} &= \frac{1}{5} \\
\frac{f^{(6)}(0)}{6!} &= 0 \\
\frac{f^{(7)}(0)}{7!} &= - \frac{1}{7}
\end{aligned}
$$

これらは次の有限展開を与える。偶数係数はすべて零なので、次の非零項は 9 次である。

$$
\arctan x = x-\frac{x^3}{3}+\frac{x^5}{5}-\frac{x^7}{7}+O(x^9)
$$
:::
:::example[Arctan の打ち切りを代入する]
七次まで残した多項式の残差から、係数がどの次数まで消えたかを調べる。

係数を二階、一階、零階の順に並べる。

$$
\left( x^{2} + 1, \  2 x, \  0\right)
$$

$$
\begin{aligned}u&=- \frac{x^{7}}{7} + \frac{x^{5}}{5} - \frac{x^{3}}{3} + x\\u^{\prime}&=- x^{6} + x^{4} - x^{2} + 1\\u^{\prime\prime}&=2 x \left(- 3 x^{4} + 2 x^{2} - 1\right)\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=2 x \left(x^{2} + 1\right) \left(- 3 x^{4} + 2 x^{2} - 1\right)+\left(2 x \left(- x^{6} + x^{4} - x^{2} + 1\right)\right)+\left(0\right)\\&=- 8 x^{7}\end{aligned}
$$

$x$ が $\frac{1}{2}$ のところで値を評価する。

$$
\begin{aligned}u&=\frac{6229}{13440}\\ u^\prime&=\frac{51}{64}\\ Lu&=- \frac{1}{16}\end{aligned}
$$
:::

## 確定特異点と Frobenius の方法

### 零点の近くで先頭の指数を選ぶ

原点で $p=xP$ と $q=x^2Q$ が正則なら、原点が正則点でない場合に :term[確定特異点]{en="regular singular point" yomi="かくていとくいてん"} という。$y=x^r\sum_{n\ge0}a_nx^n$（$a_0\ne0$）を代入する。複素平面では原点から切断を入れ、$x^r=e^{r\log x}$ の枝を固定する。最低次から :term[決定方程式]{en="indicial equation" yomi="けっていほうていしき"} $I(r)=r(r-1)+p_0r+q_0=0$ を得る。残りは
$$
I(r+n)a_n=-\sum_{j=1}^n[p_j(r+n-j)+q_j]a_{n-j}.
\label{eq-mp-frobenius-recurrence}
$$
$I(r+n)\ne0$ が全ての $n\ge1$ で成り立つ根なら解を作れる。

収束には任意の係数の収束円板内で $R>0$ を選び、$|p_j|+|q_j|\le MR^{-j}$ を使う。$0<\rho<R$ とし $b_n=|a_n|\rho^n$ と置く。$I(r+n)$ は大きな $n$ で絶対値が $c n^2$ 以上であるため、漸化式から $b_n\le Cn^{-1}\sum_{j=1}^n(\rho/R)^j b_{n-j}$ である。十分大きい $n$ では係数の総和が $1$ 以下になる。初項から有限個の最大値を取って帰納すれば全 $b_n$ は有界で、$|x|<\rho$ で絶対収束する。$\rho$ は $R$ に任意に近づけられる。

:::example[Euler 型と一般の違い]
$x^2y''-2y=0$ では $I(r)=(r-2)(r+1)$、$p_j=q_j=0$（$j\ge1$）なので $x^2$、$x^{-1}$ が独立解である。整数差だけから対数の有無は決まらない。
:::

:::example[原点の指数と Bessel 型の係数]
$x^2y''+xy'+(x^2-1/4)y=0$ では $I(r)=r^2-1/4$ で、根は $\pm1/2$ である。大きい根を選ぶと $I(1/2+n)=n(n+1)$ なので、$a_0=1$ に対し $a_2=-1/6$、$a_4=1/120$ となる。小さい根では差が $1$ であるため $n=1$ の分母は零だが、ここでは分子も零である。小さい根からは $a_1=0$ として $a_2=-1/2$、$a_4=1/24$ を得る。級数はそれぞれ $\sin x/\sqrt x$、$\cos x/\sqrt x$ に一致する。整数差でも、式が偶数だけを結ぶと対数が現れない場合を具体的に検算できる。
:::

### 根が共鳴しないときの収束

$xP=\sum p_jx^j$、$x^2Q=\sum q_jx^j$ とする。指数 $r$ を決定方程式 $I(r)=0$ で選び、各 $n>0$ で $I(r+n)\ne0$ なら係数が一意に決まる。収束を見るため $0<\rho<d$ を取り、$|p_j|+|q_j|\le C\rho^{-j}$ とする。十分大きい $n$ では $|I(r+n)|\ge C_1 n^2$ なので、漸化式を絶対値で取ると
$$
|a_n|\rho^n\le \frac{C_2}{n}\sum_{k=0}^{n-1}|a_k|\rho^k
$$
となる。$B_n=\sum_{k=0}^n|a_k|\rho^k$ とおくと $B_n\le(1+C_2/n)B_{n-1}$ である。$\log(1+v)\le v$ と調和和の対数評価を用いると $B_n\le C_3 n^{C_2}$ を得る。したがって $|x|<\rho$ では係数の増大を幾何級数が抑える。$\rho$ を係数の半径に近づければ所要の収束範囲を得る。

因子 $x^r$ については、実軸なら $x>0$、複素平面なら対数を一価に選んだ扇形を使う。二指数の差が整数でなければ二つの係数列が作れる。各先頭係数を $1$ とすれば、Wronskian の先頭項は $(r_2-r_1)x^{r_1+r_2-1}$ となる。

:::example[非整数の正指数]
正の半直線で $r=3/2$ の単純な指数解を作る。

係数を二階、一階、零階の順に並べる。

$$
\left( x^{2}, \  x, \  - \frac{9}{4}\right)
$$

$$
\begin{aligned}u&=x^{\frac{3}{2}}\\u^{\prime}&=\frac{3 \sqrt{x}}{2}\\u^{\prime\prime}&=\frac{3}{4 \sqrt{x}}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=\frac{3 x^{\frac{3}{2}}}{4}+\left(\frac{3 x^{\frac{3}{2}}}{2}\right)+\left(- \frac{9 x^{\frac{3}{2}}}{4}\right)\\&=0\end{aligned}
$$

$x$ が $4$ のところで値を評価する。

$$
\begin{aligned}u&=8\\ u^\prime&=3\\ Lu&=0\end{aligned}
$$
:::

:::example[小さい方の非整数指数]
同じ演算子を負の指数の解に適用し、原点での発散と途中の相殺を区別する。

係数を二階、一階、零階の順に並べる。

$$
\left( x^{2}, \  x, \  - \frac{9}{4}\right)
$$

$$
\begin{aligned}u&=\frac{1}{x^{\frac{3}{2}}}\\u^{\prime}&=- \frac{3}{2 x^{\frac{5}{2}}}\\u^{\prime\prime}&=\frac{15}{4 x^{\frac{7}{2}}}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=\frac{15}{4 x^{\frac{3}{2}}}+\left(- \frac{3}{2 x^{\frac{3}{2}}}\right)+\left(- \frac{9}{4 x^{\frac{3}{2}}}\right)\\&=0\end{aligned}
$$

$x$ が $4$ のところで値を評価する。

$$
\begin{aligned}u&=\frac{1}{8}\\ u^\prime&=- \frac{3}{64}\\ Lu&=0\end{aligned}
$$
:::

:::example[収束する級数の最初の項]
Bessel 型で一階下の項がなく、原点の指数から二つおきに補正が出る。

係数を二階、一階、零階の順に並べる。

$$
\left( x^{2}, \  x, \  x^{2} - \frac{1}{4}\right)
$$

$$
\begin{aligned}u&=\sqrt{x} \left(\frac{x^{4}}{120} - \frac{x^{2}}{6} + 1\right)\\u^{\prime}&=\sqrt{x} \left(\frac{x^{3}}{30} - \frac{x}{3}\right) + \frac{\frac{x^{4}}{120} - \frac{x^{2}}{6} + 1}{2 \sqrt{x}}\\u^{\prime\prime}&=\frac{2 \sqrt{x} \left(x^{2} - 10\right) + 2 \sqrt{x} \left(3 x^{2} - 10\right) - \frac{x^{4} - 20 x^{2} + 120}{8 x^{\frac{3}{2}}}}{60}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=\frac{x^{2} \left(2 \sqrt{x} \left(x^{2} - 10\right) + 2 \sqrt{x} \left(3 x^{2} - 10\right) - \frac{x^{4} - 20 x^{2} + 120}{8 x^{\frac{3}{2}}}\right)}{60}+\left(x \left(\sqrt{x} \left(\frac{x^{3}}{30} - \frac{x}{3}\right) + \frac{\frac{x^{4}}{120} - \frac{x^{2}}{6} + 1}{2 \sqrt{x}}\right)\right)+\left(\sqrt{x} \left(x^{2} - \frac{1}{4}\right) \left(\frac{x^{4}}{120} - \frac{x^{2}}{6} + 1\right)\right)\\&=\frac{x^{\frac{13}{2}}}{120}\end{aligned}
$$

$x$ が $1$ のところで値を評価する。

$$
\begin{aligned}u&=\frac{101}{120}\\ u^\prime&=\frac{29}{240}\\ Lu&=\frac{1}{120}\end{aligned}
$$
:::

:::example[指数が複素の場合]
実係数の Euler 型方程式で、対数半径に周期を持つ実解を作る。

係数を二階、一階、零階の順に並べる。

$$
\left( x^{2}, \  x, \  4\right)
$$

$$
\begin{aligned}u&=\cos{\left(2 \log{\left(x \right)} \right)}\\u^{\prime}&=- \frac{2 \sin{\left(2 \log{\left(x \right)} \right)}}{x}\\u^{\prime\prime}&=\frac{2 \left(\sin{\left(2 \log{\left(x \right)} \right)} - 2 \cos{\left(2 \log{\left(x \right)} \right)}\right)}{x^{2}}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=2 \sin{\left(2 \log{\left(x \right)} \right)} - 4 \cos{\left(2 \log{\left(x \right)} \right)}+\left(- 2 \sin{\left(2 \log{\left(x \right)} \right)}\right)+\left(4 \cos{\left(2 \log{\left(x \right)} \right)}\right)\\&=0\end{aligned}
$$

$x$ が $1$ のところで値を評価する。

$$
\begin{aligned}u&=1\\ u^\prime&=0\\ Lu&=0\end{aligned}
$$
:::

### 異なる側での指数の意味

確定特異点 $0$ の近くに複素引数を取る。収束級数 $A(x)$ は原点周りの円板内で一価なので、$u(x)=x^\rho A(x)$ を原点の周りに正の向きに一周接続すると、対数が $2\pi i$ 増え、$u$ は $e^{2\pi i\rho}u$ になる。例えば指数が $1/2,-1/2$ なら両方の乗数が $-1$ になる。整数の指数を持つ場合でも第 2 解に対数があれば、接続は単なる恒等にならない。この違いを調べるには解空間上に作用する行列を求める。

実係数の場合、正の半直線では $x^\rho=e^{\rho\log x}$ を使える。負の半直線に同じ式を持ち込む場合には、$\arg x=\pi$ と $-\pi$ の選び方で $e^{\pm i\pi\rho}$ が異なる。実軸での実解だけが欲しいときは $x=-s$ と置いて $s>0$ の問題を立て直すこともできる。その場合は左右からの延長の関係を指定しない限り、二側の積分定数は独立である。
:::example[実部と虚部で二つの実解を作る]
正の半直線で指数が $-1+i,-1-i$ の場合の実成分を調べる。

係数を二階、一階、零階の順に並べる。

$$
\left( u^{2}, \  3 u, \  2\right)
$$

$$
\begin{aligned}u&=\frac{\cos{\left(\log{\left(u \right)} \right)}}{u}\\u^{\prime}&=- \frac{\sin{\left(\log{\left(u \right)} \right)}}{u^{2}} - \frac{\cos{\left(\log{\left(u \right)} \right)}}{u^{2}}\\u^{\prime\prime}&=\frac{3 \sin{\left(\log{\left(u \right)} \right)} + \cos{\left(\log{\left(u \right)} \right)}}{u^{3}}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=\frac{3 \sin{\left(\log{\left(u \right)} \right)} + \cos{\left(\log{\left(u \right)} \right)}}{u}+\left(3 u \left(- \frac{\sin{\left(\log{\left(u \right)} \right)}}{u^{2}} - \frac{\cos{\left(\log{\left(u \right)} \right)}}{u^{2}}\right)\right)+\left(\frac{2 \cos{\left(\log{\left(u \right)} \right)}}{u}\right)\\&=0\end{aligned}
$$

$u$ が $2$ のところで値を評価する。

$$
\begin{aligned}u&=\frac{\cos{\left(\log{\left(2 \right)} \right)}}{2}\\ u^\prime&=- \frac{\cos{\left(\log{\left(2 \right)} \right)}}{4} - \frac{\sin{\left(\log{\left(2 \right)} \right)}}{4}\\ Lu&=0\end{aligned}
$$
:::
:::example[複素指数の虚部]
同じ指数から得た虚成分を別の初期データで確認する。

係数を二階、一階、零階の順に並べる。

$$
\left( u^{2}, \  3 u, \  2\right)
$$

$$
\begin{aligned}u&=\frac{\sin{\left(\log{\left(u \right)} \right)}}{u}\\u^{\prime}&=- \frac{\sin{\left(\log{\left(u \right)} \right)}}{u^{2}} + \frac{\cos{\left(\log{\left(u \right)} \right)}}{u^{2}}\\u^{\prime\prime}&=\frac{\sin{\left(\log{\left(u \right)} \right)} - 3 \cos{\left(\log{\left(u \right)} \right)}}{u^{3}}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=\frac{\sin{\left(\log{\left(u \right)} \right)} - 3 \cos{\left(\log{\left(u \right)} \right)}}{u}+\left(3 u \left(- \frac{\sin{\left(\log{\left(u \right)} \right)}}{u^{2}} + \frac{\cos{\left(\log{\left(u \right)} \right)}}{u^{2}}\right)\right)+\left(\frac{2 \sin{\left(\log{\left(u \right)} \right)}}{u}\right)\\&=0\end{aligned}
$$

$u$ が $2$ のところで値を評価する。

$$
\begin{aligned}u&=\frac{\sin{\left(\log{\left(2 \right)} \right)}}{2}\\ u^\prime&=- \frac{\sin{\left(\log{\left(2 \right)} \right)}}{4} + \frac{\cos{\left(\log{\left(2 \right)} \right)}}{4}\\ Lu&=0\end{aligned}
$$
:::

## 決定方程式の根が整数差の場合と第 2 解

### 漸化式が止まる場合

差 $r_1-r_2=m$ が非負整数のとき、$r_1$ の解 $y_1=x^{r_1}A(x)$ は上の方法で作れる。$A(0)\ne0$ とし、零点のない小領域で $y_2=y_1v$ を入れると $v''+(2y_1'/y_1+P)v'=0$ となる。よって
$$
y_2=y_1\int^x\frac{e^{-\int^t P(s)\dd s}}{y_1(t)^2}\dd t .
\label{eq-mp-reduction-order}
$$
積分の中身は $t^{-p_0-2r_1}B(t)=t^{-m-1}B(t)$、$B$ は正則、$B(0)\ne0$ である。後者の指数には根の和 $r_1+r_2=1-p_0$ を使った。$B$ の第 $m$ 係数だけが $\log x$ を与える。$m=0$ では必ず対数があり、$m>0$ ではその係数が零なら両方とも冪級数型である。積分定数は $y_1$ の倍を足すだけである。

:::example[重根と異なる整数根]
$x^2y''+xy'=0$ は $I=r^2$ を持ち、$y_1=1$ から積分は $\int\dd t/t=\log x$ になる。$xy''-y'+xy=0$ では $I=r(r-2)$ で、$r=0$ の漸化式の $n=2$ に $0=-a_0$ が現れる。この場合は小さい根から純粋な級数は作れず、第 2 解に対数が必要である。
:::

### 階数低下で対数項を見つける

既知の解を $u=x^{r_1}A(x)$ とし、$A(0)=1$ とする。第 2 解を $v=uH$ とおくと、$H''+(2u'/u+P)H'=0$ である。$P=p_0/x+P_1(x)$ の場合に積分すると
$$
H'=x^{-p_0-2r_1}B(x),\qquad B(x)=A(x)^{-2}\exp\left(-\int_0^x P_1(s)\dd s\right).
$$
根の和は $r_1+r_2=1-p_0$ なので、先頭の冪は $r_2-r_1-1$ になる。$N=r_1-r_2\in\N\cup\{0\}$ とし $B=\sum b_jx^j$ と書くと、$H'$ の中で指数が $-1$ になるのは $j=N$ である。したがって対数の係数は $b_N$ であり、重根の場合は必ず $b_0=1$ が残る。$N>0$ のときは $b_N=0$ なら二つとも対数を含まない。積分定数は $u$ の倍数に吸収できる。

:::example[重根で対数を持つ解]
Euler 型の重根に対し、第 2 解の対数を微分で確かめる。

係数を二階、一階、零階の順に並べる。

$$
\left( x^{2}, \  x, \  0\right)
$$

$$
\begin{aligned}u&=\log{\left(x \right)}\\u^{\prime}&=\frac{1}{x}\\u^{\prime\prime}&=- \frac{1}{x^{2}}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=-1+\left(1\right)+\left(0\right)\\&=0\end{aligned}
$$

$x$ が $2$ のところで値を評価する。

$$
\begin{aligned}u&=\log{\left(2 \right)}\\ u^\prime&=\frac{1}{2}\\ Lu&=0\end{aligned}
$$
:::

:::example[整数差で対数のない解]
整数差が二でも、負の指数の冪だけで方程式を満たす例がある。

係数を二階、一階、零階の順に並べる。

$$
\left( x^{2}, \  x, \  -1\right)
$$

$$
\begin{aligned}u&=\frac{1}{x}\\u^{\prime}&=- \frac{1}{x^{2}}\\u^{\prime\prime}&=\frac{2}{x^{3}}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=\frac{2}{x}+\left(- \frac{1}{x}\right)+\left(- \frac{1}{x}\right)\\&=0\end{aligned}
$$

$x$ が $2$ のところで値を評価する。

$$
\begin{aligned}u&=\frac{1}{2}\\ u^\prime&=- \frac{1}{4}\\ Lu&=0\end{aligned}
$$
:::

:::example[非定数係数が対数を作る例]
既知解 $u=1$ とすると $H^\prime=x^{-2}(1+x)$ となる方程式を取る。

係数を二階、一階、零階の順に並べる。

$$
\left( x^{2} \left(x + 1\right), \  x \left(x + 2\right), \  0\right)
$$

$$
\begin{aligned}u&=\log{\left(x \right)} - \frac{1}{x}\\u^{\prime}&=\frac{1}{x} + \frac{1}{x^{2}}\\u^{\prime\prime}&=- \frac{1 + \frac{2}{x}}{x^{2}}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=- \left(1 + \frac{2}{x}\right) \left(x + 1\right)+\left(x \left(\frac{1}{x} + \frac{1}{x^{2}}\right) \left(x + 2\right)\right)+\left(0\right)\\&=0\end{aligned}
$$

$x$ が $1$ のところで値を評価する。

$$
\begin{aligned}u&=-1\\ u^\prime&=2\\ Lu&=0\end{aligned}
$$
:::

:::example[共鳴係数の抽出]
上の第 2 解では $B(x)=1+x$ の一次係数が対数の係数になる。別の例として $B(x)=(1+x)^2 e^{-x}$ を取る。共鳴段 $N$ が分かれば、その次数の係数を読むだけで対数項が求まる。

まず展開する関数と中心を定める。中心で微分してから階乗で割ると各係数が得られる。

$$
f(x)=\left(x + 1\right)^{2} e^{- x},\qquad x_0=0
$$

計算した係数を、定数項から順に並べる。

$$
\begin{aligned}
\frac{f^{(0)}(0)}{0!} &= 1 \\
\frac{f^{(1)}(0)}{1!} &= 1 \\
\frac{f^{(2)}(0)}{2!} &= - \frac{1}{2} \\
\frac{f^{(3)}(0)}{3!} &= - \frac{1}{6} \\
\frac{f^{(4)}(0)}{4!} &= \frac{5}{24} \\
\frac{f^{(5)}(0)}{5!} &= - \frac{11}{120}
\end{aligned}
$$

これらは次の有限展開を与える。剰余に含まれる冪は、最後に残した項の一つ次の次数である。

$$
\left(x + 1\right)^{2} e^{- x} = - \frac{11 x^{5}}{120} + \frac{5 x^{4}}{24} - \frac{x^{3}}{6} - \frac{x^{2}}{2} + x + 1+O\bigl((x)^{6}\bigr)
$$
:::

### 指数を微分して第 2 解を作る

決定多項式を $I$、方程式を $Ly=x^2y''+xp(x)y'+q(x)y=0$ と書く。式 :eqref[eq-mp-frobenius-recurrence] に従い、指数 $s$ を動かしながら $F(s,x)=x^s\sum_{j\ge0}a_j(s)x^j$、$a_0(s)=1$ を作る。$n\ge1$ の係数だけを消したので $LF(s,x)=I(s)x^s$ が残る。収束証明の下限 $|I(s+n)|\ge c n^2$ は $s$ のコンパクト集合で取れるため、有限個の分母の零を避ければ級数は局所一様収束し、Cauchy の積分公式により $s$ で微分できる。

$r$ が重根なら $I(r)=I'(r)=0$ で、$n\ge1$ の分母は $r$ の近くで零にならない。従って $\partial_s F(r,x)$ は同次解である。$x^s$ を微分した項から $F(r,x)\log x$ が現れる。

$r_1-r_2=N>0$ が整数なら、$r_2$ の近くで零になり得る分母は $I(s+N)$ だけで、その零は単純である。従って $G(s,x)=(s-r_2)F(s,x)$ の全係数は $r_2$ の近くで正則である。最初の $N+1$ 個の係数を調べてから同じ収束評価を用いれば、$G$ に対しても $s$ 微分と $L$ を交換できる。$LG=(s-r_2)I(s)x^s=O((s-r_2)^2)$ より $G(r_2,x)$ と $\partial_sG(r_2,x)$ はともに同次解である。前者は $r_1$ の解の定数倍で、共鳴段の分子が零なら恒等的に零である。後者は $x^{r_2}$ の係数が $1$ で、対数項は $G(r_2,x)\log x$ である。この式によって、共鳴段の分子と対数項の有無が直接対応する。

例えば $xy''-y'+xy=0$ を $x$ 倍すると $I(s)=s(s-2)$ である。$a_2(s)=-1/[s(s+2)]$ なので $G(0,x)$ は $-x^2/2$ から始まる。$x^2$ で始まる解を $y_1$ と規格化すれば、第 2 解の対数項は $-\frac12 y_1\log x$ になる。

### 対数解を一周接続する

階数低下で $v=x^{\rho_2}B(x)+c u(x)\log x$ を得たとする。指数差が整数なので、$u$ と前半に掛かる一周の位相は同じ $\lambda=e^{2\pi i\rho_1}$ である。従って $v$ の像は $\lambda v+2\pi i c\lambda u$ になる。基底を $(u,v)$ として、列に各基底の像を並べれば、接続行列は対角が $\lambda$、右上に $2\pi i c\lambda$ がある Jordan 型になる。対数の係数 $c$ が零なら対角行列になり、零でなければ基底を変えても対角化できない。実際、固有ベクトルの次元が前者で二、後者で一だからである。
:::example[整数の重根での一周]
基底を $1,\log z$ とした場合、第二列に対数の増分が入る。

各列は解析接続した基底ベクトルの座標である。

$$
\begin{aligned}
A &= \begin{pmatrix}1 & 2 i \pi\\0 & 1\end{pmatrix}
\end{aligned}
$$

固有値を得るには $tI-A$ の行列式を求める。

$$
\begin{aligned}
\det(tI-A) &= (t-1)^2
\end{aligned}
$$

固有値は $1$ で、代数的重複度は $2$ である。

$$
\begin{aligned}
(\lambda_j)_j &= \left( 1,\right)
\end{aligned}
$$
:::
:::example[半整数の重根での一周]
基底を $z^{1/2},z^{1/2}\log z$ とすると、対数だけでなく冪も一周で符号を変える。

各列は解析接続した基底ベクトルの座標である。

$$
\begin{aligned}
A &= \begin{pmatrix}-1 & - 2 i \pi\\0 & -1\end{pmatrix}
\end{aligned}
$$

固有値を得るには $tI-A$ の行列式を求める。

$$
\begin{aligned}
\det(tI-A) &= (t+1)^2
\end{aligned}
$$

固有値は $-1$ で、代数的重複度は $2$ である。

$$
\begin{aligned}
(\lambda_j)_j &= \left( -1,\right)
\end{aligned}
$$
:::

## Fuchs 型方程式と Riemann の P 記号

### 無限遠も局所座標で調べる

有限点では上の条件を使い、無限遠では $w=1/z$ とおく。$y_z=-w^2y_w$、$y_{zz}=w^4y_{ww}+2w^3y_w$ だから $w$ での係数は $2/w-P(1/w)/w^2$、$Q(1/w)/w^4$ となる。従って無限遠が正則点または確定特異点なら $P=O(z^{-1})$、$Q=O(z^{-2})$ である。Riemann 球面上で有理係数がこの条件と全ての有限極での確定性を満たす方程式を :term[フックス型方程式]{en="Fuchsian equation" yomi="ふっくすがたほうていしき"} という。

特異点が $0,1,\infty$ だけで、指数を順に $(\alpha,\alpha'),(\beta,\beta'),(\gamma,\gamma')$ とした解空間を $P\{0,1,\infty;\alpha,\beta,\gamma;\alpha',\beta',\gamma';z\}$ と表す。規約では無限遠の指数 $\gamma$ は $z^{-\gamma}$ である。部分分数で $P(z)=A/z+B/(z-1)$、$Q=C/z^2+D/(z-1)^2+E/[z(z-1)]$ と書ける。各有限点の根の和は $1-A,1-B$、無限遠での根の和は $A+B-1$ で、全指数の和は $1$ になる。各有限点の根の積は $C,D$、無限遠の根の積は $C+D+E$ なので、この三点の場合は指数から係数が定まる。例えば [Gauss の超幾何方程式](/physics/mathematical-physics/hypergeometric/) は指数が $(0,1-c),(0,c-a-b),(a,b)$ である。

### 特異点が四つあると指数以外の量が残る

点 $0,1,a,\infty$（$a\ne0,1$）を考え、$P=0$、$Q=\kappa/[x(x-1)(x-a)]$ とする。各有限点で $Q$ はせいぜい単純極なので決定多項式は $r(r-1)$ である。無限遠でも $Q=O(x^{-3})$ だから指数は $r(r+1)=0$ の根である。従って $\kappa\ne0$ の範囲で全指数が同じでも、方程式は $\kappa$ によって異なる。三特異点で指数だけから決まるという性質は、この四点の例には延びない。例えば $a=2$ なら、原点近くの小さい根 $r=0$ の漸化式の $n=1$ では右辺に $-\kappa a_0/2$ が残る。従って $\kappa\ne0$ では対数が必要であり、その係数は値に依存する。

### 三つの特異点から係数を組み立てる

$0,1,\infty$ だけが確定特異点の候補で、解を $z$ の多価関数として扱う。正規形の一階項は $A/z+B/(z-1)$、零階項は $C/z^2+D/(z-1)^2+E/(z(z-1))$ である。これは有限極で主部を引くと整関数が残り、無限遠で $P=O(1/z),Q=O(1/z^2)$ によってその整関数が零になるからである。有限点の指数の和と積から $A,B,C,D$ が求まり、無限遠での指数の積から $C+D+E$ が求まる。最後に無限遠での指数和が $A+B-1$ となる。この値と有限点の和 $(1-A)+(1-B)$ を足すと $1$ を得る。

$y=gv$、$g=z^\alpha(1-z)^\beta$ と置くと、$v$ の係数に $g'/g=\alpha/z+\beta/(z-1)$ が現れる。$\alpha,\beta$ に有限点の指数を一つずつ選ぶと、新しい零階項 $Q+Pg'/g+g''/g$ の二位の極の係数は各点の決定多項式の値になり、零となる。従って $v$ の指数は $0$ で $(0,\alpha'-\alpha)$、$1$ で $(0,\beta'-\beta)$、無限遠で $(\gamma+\alpha+\beta,\gamma'+\alpha+\beta)$ になる。三点の係数の一意性から $v$ は Gauss の方程式を満たし、そのパラメータは $c=1+\alpha-\alpha'$、$a=\gamma+\alpha+\beta$、$b=\gamma'+\alpha+\beta$ である。三点からの再構成は :ref[thm-ode-riemann-p]、異なる点の解を結ぶ係数は [超幾何関数とその他の特殊関数](/physics/mathematical-physics/hypergeometric/) を参照せよ。

### 有限特異点が任意個ある場合の Fuchs の関係式

相異なる有限点 $z_1,\ldots,z_m$ で $P$ の留数を $A_i$ とする。$P=\sum_i A_i/(z-z_i)$ であるから、有限点の指数和をすべて加えると $m-\sum_i A_i$ となる。無限遠の局所座標では一階項の留数が $2-\sum_i A_i$ なので、無限遠の指数和は $\sum_i A_i-1$ である。従って $m+1$ 点での指数をすべて加えると $m-1$ である。点が実は正則な場合も、指数の対 $(0,1)$ を含めて数えれば同じ計算が成り立つ。

有限の三点 $a,b,c$ を $0,1,\infty$ へ移すには $\xi=(z-a)(b-c)/[(z-c)(b-a)]$ と取る。各点の局所座標は非零の一階係数を持って変換されるため、Frobenius 解の先頭指数は保存される。対数の引数は、変換前後の局所座標の比の対数だけずれる。その比は点の近傍で正則かつ非零なので、このずれは正則因子に吸収できる。従って逆変換によって任意の位置にある三つの特異点にも上の構成を適用できる。

### 指数から方程式を復元する例

:::example[有限点の分数冪と無限遠の規約]
指数を $0$ で $(1/2,0)$、$1$ で $(1/3,0)$、無限遠で $(1/6,0)$ とする。指数和は $1/2+1/3+1/6=1$ であり、有限点の積と無限遠の積は全て零である。従って $A=1/2$、$B=2/3$、$C=D=E=0$ であり、正規形は
$$
 y''+\left(\frac1{2z}+\frac2{3(z-1)}\right)y'=0
$$
となる。一つの解は定数 $1$ である。もう一つは $0<z<1$ で $y_2(z)=\int_0^z t^{-1/2}(1-t)^{-2/3}\dd t$ と取れる。実際、一階の積分因子を微分すると $y_2'$ が方程式を満たす。原点では $y_2(z)=2z^{1/2}(1+O(z))$ であり、$z\uparrow1$ では別の定数 $d=\lim_{z\uparrow1}y_2(z)$ を引くと $y_2(z)-d=-3(1-z)^{1/3}(1+O(1-z))$ である。同じ二次元空間でも、局所指数に対応する基底を選び直す点に注意せよ。

無限遠で枝を固定すると $y_2'$ は非零定数倍の $z^{-7/6}(1+O(1/z))$ である。従って定数を引いた $y_2$ の先頭は $z^{-1/6}$ となる。これは無限遠の局所座標 $1/z$ では正の指数 $1/6$ に対応し、有限点の指数と同じ符号の規約で読める。
:::

:::example[特異点に見えるが両解が正則な場合]
$y''-y'/z=0$ の正規形では $z=0$ に極がある。決定多項式は $r(r-2)$ で、独立解は $1,z^2$ だから両方が原点でも正則である。逆に $z\ne0$ の二解を初期値から作ると、この二解の線形結合になるので、どの解でも原点の一周で変化しない。このように方程式の係数には特異点があっても、全局所解が極のない一価関数として延びる場合がある。整数差という情報だけでは、この例と対数を含む例を区別できない。
:::

:::example[Gauss 方程式の終端例]
指数から作った方程式が一次の多項式を持つ例を取る。

係数を二階、一階、零階の順に並べる。

$$
\left( x \left(1 - x\right), \  2 - 2 x, \  2\right)
$$

$$
\begin{aligned}u&=1 - x\\u^{\prime}&=-1\\u^{\prime\prime}&=0\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=0+\left(2 x - 2\right)+\left(2 - 2 x\right)\\&=0\end{aligned}
$$

$x$ が $0$ のところで値を評価する。

$$
\begin{aligned}u&=1\\ u^\prime&=-1\\ Lu&=0\end{aligned}
$$
:::

:::example[有限点に平方根を持つ解]
有限点 $0$ の指数 $1/2$ と $1$ の指数 $0$ を選ぶ。

係数を二階、一階、零階の順に並べる。

$$
\left( x \left(1 - x\right), \  \frac{1}{2} - \frac{3 x}{2}, \  \frac{1}{2}\right)
$$

$$
\begin{aligned}u&=\sqrt{x}\\u^{\prime}&=\frac{1}{2 \sqrt{x}}\\u^{\prime\prime}&=- \frac{1}{4 x^{\frac{3}{2}}}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=- \frac{1 - x}{4 \sqrt{x}}+\left(\frac{\frac{1}{2} - \frac{3 x}{2}}{2 \sqrt{x}}\right)+\left(\frac{\sqrt{x}}{2}\right)\\&=0\end{aligned}
$$

$x$ が $\frac{1}{4}$ のところで値を評価する。

$$
\begin{aligned}u&=\frac{1}{2}\\ u^\prime&=1\\ Lu&=0\end{aligned}
$$
:::

:::example[無限遠での指数と多項式]
次数二の多項式の無限遠指数は $-2$ である。この規格化は $x=0$ で値 $1$ を取る。

係数を二階、一階、零階の順に並べる。

$$
\left( x \left(1 - x\right), \  1 - 2 x, \  6\right)
$$

$$
\begin{aligned}u&=6 x^{2} - 6 x + 1\\u^{\prime}&=12 x - 6\\u^{\prime\prime}&=12\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=12 x \left(1 - x\right)+\left(\left(1 - 2 x\right) \left(12 x - 6\right)\right)+\left(36 x^{2} - 36 x + 6\right)\\&=0\end{aligned}
$$

$x$ が $0$ のところで値を評価する。

$$
\begin{aligned}u&=1\\ u^\prime&=-6\\ Lu&=0\end{aligned}
$$
:::

:::example[一点での残差だけでは解を判定できない]
$u=x(1-x)$ を試行関数に取る。両端の零点を与えただけでは途中で方程式を満たすとは限らないので、残差を区間全体で調べる。

係数を二階、一階、零階の順に並べる。

$$
\left( x \left(1 - x\right), \  1 - 2 x, \  2\right)
$$

$$
\begin{aligned}u&=x \left(1 - x\right)\\u^{\prime}&=1 - 2 x\\u^{\prime\prime}&=-2\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=- 2 x \left(1 - x\right)+\left(\left(1 - 2 x\right)^{2}\right)+\left(2 x \left(1 - x\right)\right)\\&=4 x^{2} - 4 x + 1\end{aligned}
$$

$x$ が $\frac{1}{2}$ のところで値を評価する。

$$
\begin{aligned}u&=\frac{1}{4}\\ u^\prime&=0\\ Lu&=0\end{aligned}
$$
従って $x=1/2$ で値と傾斜を与え、そこで方程式まで満たしたように見えても、$Lu=(2x-1)^2$ は恒等的に零ではない。正則点での初期値問題の一意性を使って解を特定するには、近傍全体で方程式が成り立つことを確認する必要がある。
:::

### 変数の逆数で無限遠を調べる

$z=1/x$、$y(x)=v(z)$ とする。$\partial_x=-z^2\partial_z$ なので、$y''=z^4v''+2z^3v'$ である。従って $v$ の一階係数は $2/z-P(1/z)/z^2$、零階係数は $Q(1/z)/z^4$ になる。有限点の判定をこの式に使うと、無限遠での確定性は $P(x)=O(x^{-1}),Q(x)=O(x^{-2})$ に対応する。ただしここでは $P,Q$ が無限遠で有理型であることを仮定する。$xP(x)\to A$、$x^2Q(x)\to B$ なら、無限遠の決定式は $\rho(\rho-1)+(2-A)\rho+B=0$ になる。$x$ の正の冪 $x^m$ は $z^{-m}$ なので、無限遠での指数は $-m$ である。
:::example[無限遠座標で冪が入れ替わる]
$P(x)=3/x$、$Q(x)=0$ を逆数座標に移し、元の減衰解 $x^{-2}$ を調べる。

係数を二階、一階、零階の順に並べる。

$$
\left( u^{2}, \  - u, \  0\right)
$$

$$
\begin{aligned}u&=u^{2}\\u^{\prime}&=2 u\\u^{\prime\prime}&=2\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=2 u^{2}+\left(- 2 u^{2}\right)+\left(0\right)\\&=0\end{aligned}
$$

$u$ が $1$ のところで値を評価する。

$$
\begin{aligned}u&=1\\ u^\prime&=2\\ Lu&=0\end{aligned}
$$
:::
:::example[Legendre 型の無限遠側]
$(1-x^2)y^{\prime\prime}-2xy^\prime+2y=0$ を逆数座標へ移し、$y=x$ に対応する枝を調べる。

係数を二階、一階、零階の順に並べる。

$$
\left( u^{2} \left(u^{2} - 1\right), \  2 u^{3}, \  2\right)
$$

$$
\begin{aligned}u&=\frac{1}{u}\\u^{\prime}&=- \frac{1}{u^{2}}\\u^{\prime\prime}&=\frac{2}{u^{3}}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=\frac{2 \left(u^{2} - 1\right)}{u}+\left(- 2 u\right)+\left(\frac{2}{u}\right)\\&=0\end{aligned}
$$

$u$ が $\frac{1}{2}$ のところで値を評価する。

$$
\begin{aligned}u&=2\\ u^\prime&=-4\\ Lu&=0\end{aligned}
$$
:::

## 不確定特異点と漸近解

### 級数が収束しない点での計算

$xP$ または $x^2Q$ が正則に延びない孤立特異点を :term[不確定特異点]{en="irregular singular point" yomi="ふかくていとくいてん"} という。指数因子を先に分けても残りが収束するとは限らない。例えば $x^2y'+y=x$ の同次解は $e^{1/x}$ であり、形式解 $y=\sum_{n\ge1}a_nx^n$ は $a_1=1,a_n=-(n-1)a_{n-1}$ なので $a_n=(-1)^{n-1}(n-1)!$ となる。これは任意の $x\ne0$ で発散する。$x>0$ では $y=x\int_0^\infty e^{-t}/(1+xt)\dd t$ が実際の解である。$u=t+1/x$ とすると $y=e^{1/x}\int_{1/x}^\infty e^{-u}\dd u/u$ なので、微分すると確かめられる。有限等比和の剰余により形式級数がこの解の漸近展開になる。方向を変えると $e^{1/x}$ の大小が変わるため、形式級数だけでは全方向の境界条件は決まらない。次にこの方程式で積分路を変え、方向による差を求める。

### 積分路が極を越えると生じる Stokes の跳び

$u=xt$ と変数を替えれば、正実軸での解は $\int_0^\infty e^{-u/x}(1+u)^{-1}\dd u$ である。$u$ の積分路を偏角 $\alpha$ の半直線 $\mathcal C_\alpha$ へ移し、
$$
 y_\alpha(x)=\int_{\mathcal C_\alpha}\frac{e^{-u/x}}{1+u}\dd u
$$
とおく。$\alpha\not\equiv\pi\pmod{2\pi}$ で、同じ枝で $|\alpha-\arg x|<\pi/2$ なら積分は絶対収束する。積分内で $x^2\partial_x+1$ を作用させると、分母を打ち消して $\int_{\mathcal C_\alpha}e^{-u/x}\dd u=x$ となる。したがってどの方向の積分も同じ非斉次方程式を満たす。

有限等比和によって $1/(1+u)=\sum_{j=0}^{N-1}(-u)^j+(-u)^N/(1+u)$ と分ける。各単項式の Laplace 積分は、指数の減衰域で積分路を回せば $j!x^{j+1}$ になる。半直線に沿って $|1+u|$ は正の下限を持つ。さらに $|\alpha-\arg x|\le\pi/2-\delta$（$\delta>0$）では指数が $e^{-c|u|/|x|}$ で抑えられるため、剰余は $O(|x|^{N+1})$ である。したがってこの範囲で、どの積分にも同じ形式級数が対応する。

負の実軸を上から避ける積分を $y_+$、下から避ける積分を $y_-$ とする。$x$ が負実軸の近くにあれば両者の指数は遠方で減衰する。上側を原点から左へ進み、下側を逆に戻る閉曲線は極 $u=-1$ を正向きに囲む。従って
$$
y_+(x)-y_-(x)=2\pi i\,e^{1/x}.
$$
右辺は同次解であり、この方向では任意の代数次数より速く零へ近づく。従って各冪の漸近係数は同じでも、積分路によって実際の解が変わる。このようなセクターごとの解の接続を :term[ストークス現象]{en="Stokes phenomenon" yomi="すとーくすげんしょう"} という。一般の系のセクター解の存在は、本章の積分法が使えない場合を含むため、:cite[coddington-levinson] を参照せよ。

### 支配的な指数を先に求める

$y''=x y$ の $x\to+\infty$ では $y=e^{S}$ と置けば $(S')^2+S''=x$ である。最初に $(S')^2=x$ とすれば $S'=\pm\sqrt x$ で、$S''$ は $O(x^{-1/2})$ となり主項 $x$ より小さい。$S'=\pm\sqrt x+g$ を戻すと、次の次数では $\pm2\sqrt x g\pm(2\sqrt x)^{-1}=0$ となる。従って $g=-1/(4x)$、候補は $x^{-1/4}e^{\pm2x^{3/2}/3}$ である。この代入だけで真の解との差の評価は出ない。減衰側については [Airy 関数の輪郭積分](/physics/mathematical-physics/hypergeometric/) で、Gaussian の局所化から実際にこの冪と指数を得る。

### 指数因子を求めてから剰余を分ける

二階方程式を $y=e^{S}v$ とおくと、$y'/e^S=v'+S'v$、$y''/e^S=v''+2S'v'+(S''+(S')^2)v$ となる。このため、零階項で不確定特異点を生む強い極を $(S')^2$ と相殺させることができる。例えば $y''=x^{-4}y$ で支配項を釣り合わせると $S'=\pm x^{-2}$ となる。$S''$ の分を振幅で補正すると $v=x$ が得られ、二つの解が厳密になる。正の半直線で一方は増大し、他方は任意の代数次数より速く減衰する。

一般には形式計算の残差が小さくても、誤差が小さいと結論するには境界条件が要る。一次方程式 $v'+A v=r$ に積分因子を掛けると、解の誤差は残差の積分と同次解の和になる。このうち同次解の係数は形式級数から決まらない。係数を決めるには、どの方向のどの端で成長を抑えるかを指定する。

:::example[正の半直線の減衰解]
釣合いで得た指数のうち、$x\to0^+$ で減衰する方を検証する。

係数を二階、一階、零階の順に並べる。

$$
\left( x^{4}, \  0, \  -1\right)
$$

$$
\begin{aligned}u&=x e^{- \frac{1}{x}}\\u^{\prime}&=e^{- \frac{1}{x}} + \frac{e^{- \frac{1}{x}}}{x}\\u^{\prime\prime}&=\frac{e^{- \frac{1}{x}}}{x^{3}}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=x e^{- \frac{1}{x}}+\left(0\right)+\left(- x e^{- \frac{1}{x}}\right)\\&=0\end{aligned}
$$

$x$ が $1$ のところで値を評価する。

$$
\begin{aligned}u&=e^{-1}\\ u^\prime&=\frac{2}{e}\\ Lu&=0\end{aligned}
$$
:::

:::example[もう一方の指数解]
同じ微分方程式でも初期データによって増大する成分が入る。

係数を二階、一階、零階の順に並べる。

$$
\left( x^{4}, \  0, \  -1\right)
$$

$$
\begin{aligned}u&=x e^{\frac{1}{x}}\\u^{\prime}&=e^{\frac{1}{x}} - \frac{e^{\frac{1}{x}}}{x}\\u^{\prime\prime}&=\frac{e^{\frac{1}{x}}}{x^{3}}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=x e^{\frac{1}{x}}+\left(0\right)+\left(- x e^{\frac{1}{x}}\right)\\&=0\end{aligned}
$$

$x$ が $1$ のところで値を評価する。

$$
\begin{aligned}u&=e\\ u^\prime&=0\\ Lu&=0\end{aligned}
$$
:::

:::example[階乗係数に対応する残差]
一次方程式 $x^2y^\prime+y=x$ に四項を戻す。零階項と微分項の和を $x$ と比べる。

係数を二階、一階、零階の順に並べる。

$$
\left( 0, \  x^{2}, \  1\right)
$$

$$
\begin{aligned}u&=- 6 x^{4} + 2 x^{3} - x^{2} + x\\u^{\prime}&=- 24 x^{3} + 6 x^{2} - 2 x + 1\\u^{\prime\prime}&=2 \left(- 36 x^{2} + 6 x - 1\right)\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=0+\left(x^{2} \left(- 24 x^{3} + 6 x^{2} - 2 x + 1\right)\right)+\left(- 6 x^{4} + 2 x^{3} - x^{2} + x\right)\\&=- 24 x^{5} + x\end{aligned}
$$

$x$ が $1$ のところで値を評価する。

$$
\begin{aligned}u&=-4\\ u^\prime&=-19\\ Lu&=-23\end{aligned}
$$
:::

:::example[振幅を省いたときの残差]
同じ指数を使っても $v=1$ では $S^{\prime\prime}$ を打ち消せないことを確かめる。

係数を二階、一階、零階の順に並べる。

$$
\left( x^{4}, \  0, \  -1\right)
$$

$$
\begin{aligned}u&=e^{- \frac{1}{x}}\\u^{\prime}&=\frac{e^{- \frac{1}{x}}}{x^{2}}\\u^{\prime\prime}&=\frac{\left(-2 + \frac{1}{x}\right) e^{- \frac{1}{x}}}{x^{3}}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=x \left(-2 + \frac{1}{x}\right) e^{- \frac{1}{x}}+\left(0\right)+\left(- e^{- \frac{1}{x}}\right)\\&=- 2 x e^{- \frac{1}{x}}\end{aligned}
$$

$x$ が $1$ のところで値を評価する。

$$
\begin{aligned}u&=e^{-1}\\ u^\prime&=e^{-1}\\ Lu&=- \frac{2}{e}\end{aligned}
$$
:::

### 正規形で極の強さを測る

$y''+Py'+Qy=0$ で $y=e^{-\frac12\int P}v$ と置けば、$v''=Vv$、$V=P'/2+P^2/4-Q$ である。これは局所的に積分の枝を固定して使う変換である。$V(x)=c x^{-m}(1+O(x))$、$m>2$、$c\ne0$ の場合、$v=e^S$ の式 $S''+(S')^2=V$ では、$S'\sim\pm \sqrt c\,x^{-m/2}$ と仮定すると $S''/(S')^2=O(x^{m/2-1})$ となる。従って先頭の指数関数の指数は
$$
 S(x)\sim\pm\frac{\sqrt c}{1-m/2}\,x^{1-m/2}
$$
と予想される。指数が分数冪を含むときには、$x$ の対数の枝をここでも固定する。$m=2$ では微分項と二乗項が同じ次数なので、片方だけで釣り合わせず決定方程式に戻る。

例えば $x^3y''=y$ では $m=3$、$c=1$ で、指数の候補は $\pm2/\sqrt x$ である。$x^{-m/2}$ の二乗と導関数が相殺しないため、候補そのままの指数関数は厳密解にならない。次の変換を適用するときは $Q=x^{-3}$ を使うので振幅は $x^{3/4}$ になる。さらに $\xi=2/\sqrt x$ と置き、$x\downarrow0$ に向けて $\xi$ が増えるよう座標の符号を取れば、変換後の補正は $\psi=3/(4\xi^2)$ である。その積分は無限遠で有限なので、減衰側は $x^{3/4}e^{-2/\sqrt x}(1+O(\sqrt x))$ となる。

### 正の係数に対する Liouville–Green 近似

$y''=Q(x)y$、$Q>0$、$Q\in C^2$ とし、$\xi=\int_{x_0}^x\sqrt{Q(t)}\dd t$、$y=Q^{-1/4}v(\xi)$ と置く。微分すると $v_\xi$ の係数は $2(Q^{-1/4})'\sqrt Q+Q^{-1/4}(\sqrt Q)'=0$ となる。残る方程式は
$$
v_{\xi\xi}=(1+\psi)v,\qquad
\psi=\frac{Q''}{4Q^2}-\frac{5(Q')^2}{16Q^3}
$$
である。$\xi\to\infty$ で $\int^\infty|\psi(\eta)|\dd\eta<\infty$ と仮定する。減衰解を $v=e^{-\xi}h$ と書き、$h\to1$、$h'\to0$ を目標に $h''-2h'=\psi h$ を積分すると
$$
h(\xi)=1+\frac12\int_\xi^\infty(1-e^{-2(\eta-\xi)})\psi(\eta)h(\eta)\dd\eta
$$
を得る。この式から解を作るには、右辺の積分作用素を $T$ とし、$\sum_{n\ge0}T^n1$ を取る。各積分領域が $\xi\le\eta_1\le\cdots\le\eta_n$ であるため、$|T^n1|\le (\int_\xi^\infty|\psi|/2)^n/n!$ である。従って和は一様収束し、式を微分すれば所定の方程式を満たす。特に $h-1=O(\int_\xi^\infty|\psi|)$ である。十分遠方では $v$ は非零なので、第 2 解を $2v(\xi)\int_{\xi_0}^\xi v(\eta)^{-2}\dd\eta$ と取れる。積分と $e^{2\xi}/2$ の比は、導関数の比 $h(\xi)^{-2}\to1$ により $1$ に近づく。従って第 2 解は $e^\xi(1+o(1))$ である。

$Q=x$ なら定数を調整して $\xi=2x^{3/2}/3$、$\psi=-5/(16x^3)=-5/(36\xi^2)$ と書ける。よって前に得た減衰側の候補は、相対誤差 $O(x^{-3/2})$ で真の解に対応する。転回点 $Q=0$ では変数変換の微分と振幅因子の仮定が破れるため、この評価をそこまで延長することはできない。

### 指数の向きを実部で比較する

$x=re^{i\theta}$ として $e^{1/x}$ の絶対値を取ると $\exp(\cos\theta/r)$ になる。従って $-\pi/2<\theta<\pi/2$ では増大し、$\pi/2<\theta<3\pi/2$ では減衰する。境界の二方向では絶対値が $1$ なので、同じ代数次数の剰余評価を閉扇形全体で使う場合には境界を余計に含められない。例えば $x e^{-1/x}$ は正の軸に沿って任意の $N$ に対し $O(x^N)$ である。この評価は $t=1/x$ を入れ、$t^{N-1}e^{-t}$ が有界であることから従う。ところが負の軸上ではその絶対値が発散する。形式的な零級数を持つ解と恒等的に零の解は、減衰方向の有限次数の係数だけでは区別できない。
:::example[実指数の振幅の違い]
二つの指数を足した場合にも線形性で厳密解になり、片方の成分だけで誤差を見積もるときには方向の情報が要る。

係数を二階、一階、零階の順に並べる。

$$
\left( u^{4}, \  0, \  -1\right)
$$

$$
\begin{aligned}u&=u \left(2 e^{\frac{1}{u}} + e^{- \frac{1}{u}}\right)\\u^{\prime}&=u \left(- \frac{2 e^{\frac{1}{u}}}{u^{2}} + \frac{e^{- \frac{1}{u}}}{u^{2}}\right) + 2 e^{\frac{1}{u}} + e^{- \frac{1}{u}}\\u^{\prime\prime}&=\frac{\frac{2 e^{\frac{1}{u}}}{u} + \frac{e^{- \frac{1}{u}}}{u}}{u^{2}}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=u^{2} \left(\frac{2 e^{\frac{1}{u}}}{u} + \frac{e^{- \frac{1}{u}}}{u}\right)+\left(0\right)+\left(- u \left(2 e^{\frac{1}{u}} + e^{- \frac{1}{u}}\right)\right)\\&=0\end{aligned}
$$

$u$ が $1$ のところで値を評価する。

$$
\begin{aligned}u&=e^{-1} + 2 e\\ u^\prime&=\frac{2}{e}\\ Lu&=0\end{aligned}
$$
:::

## Wronskian と Abel の公式

### 二解の独立性を調べる

:term[ロンスキアン]{en="Wronskian" yomi="ろんすきあん"} を $W=y_1y_2'-y_1'y_2$ とする。微分すると $W'=y_1y_2''-y_1''y_2=-PW$ である。積分すれば
$$
W(x)=W(x_0)\exp\left[-\int_{x_0}^xP(t)\dd t\right]
\label{eq-mp-abel-wronskian}
$$
を得る。これが Abel の公式である。実変数なら係数が連続な区間で、複素変数なら係数が正則で積分路を固定した領域で用いる。領域内の一点で $W$ が零なら、初期値ベクトルが一次従属であり、一意性から二解は一次従属である。零でなければ全点で零でない。式 :eqref[eq-mp-reduction-order] は $W=y_1^2v'$ と読んでも導ける。例えば $P=1/x$ の $1,\log x$ の Wronskian は $1/x$ である。位置に依存する Wronskian の倍率は同じ基底で積分を行えば相殺される。

### Wronskian から初期値を合わせる

$y=c_1u+c_2v$ とし、点 $x_0$ で $(y,y')=(a,b)$ を指定する。二成分を並べた行列の行列式は $W(x_0)$ である。非零なら Cramer の公式により $c_1=(av'-bv)/W$、$c_2=(bu-au')/W$ となる。右辺の全関数を同じ点 $x_0$ で評価する。さらに右辺 $f(x)$ を持つ問題で係数を関数とし、$c_1'u+c_2'v=0$ と置けば、第二行は $c_1'u'+c_2'v'=f$ になる。したがって $c_1'=-v f/W$、$c_2'=u f/W$ である。積分の下限を $x_0$ にすると、非斉次応答の初期値が両方とも零になる。先頭係数が $1$ でない場合は先に右辺をその係数で割る。

:::example[実指数の対]
異なる指数率が混ざると、Wronskian 自体も指数関数になる。

先に関数と微分を並べる。微分する順と関数の順を固定しておくと符号が決まる。

$$
\begin{aligned}
u &= e^{x} \\
v &= e^{2 x} \\
u^\prime &= e^{x} \\
v^\prime &= 2 e^{2 x}
\end{aligned}
$$

積 $uv^\prime$ から $u^\prime v$ を引いて簡約する。

$$
\begin{aligned}
uv^\prime &= 2 e^{3 x} \\
u^\prime v &= e^{3 x} \\
W &= e^{3 x}
\end{aligned}
$$

Abel の式を逆に使うと、共通の方程式の一階項の係数が求まる。$W$ が零にならない区間で計算する。

$$
-W^\prime/W=-3
$$
:::

:::example[重根の場合]
指数因子を共有するときは $x$ を掛けた解と組み合わせる。

先に関数と微分を並べる。微分する順と関数の順を固定しておくと符号が決まる。

$$
\begin{aligned}
u &= e^{- 2 x} \\
v &= x e^{- 2 x} \\
u^\prime &= - 2 e^{- 2 x} \\
v^\prime &= - 2 x e^{- 2 x} + e^{- 2 x}
\end{aligned}
$$

積 $uv^\prime$ から $u^\prime v$ を引いて簡約する。

$$
\begin{aligned}
uv^\prime &= \left(- 2 x e^{- 2 x} + e^{- 2 x}\right) e^{- 2 x} \\
u^\prime v &= - 2 x e^{- 4 x} \\
W &= e^{- 4 x}
\end{aligned}
$$

Abel の式を逆に使うと、共通の方程式の一階項の係数が求まる。$W$ が零にならない区間で計算する。

$$
-W^\prime/W=4
$$
:::

:::example[確定特異点での二つの冪]
点 $0$ を除いた区間で計算する。Wronskian が発散しても区間内での独立性は保たれる。

先に関数と微分を並べる。微分する順と関数の順を固定しておくと符号が決まる。

$$
\begin{aligned}
u &= x \\
v &= \frac{1}{x} \\
u^\prime &= 1 \\
v^\prime &= - \frac{1}{x^{2}}
\end{aligned}
$$

積 $uv^\prime$ から $u^\prime v$ を引いて簡約する。

$$
\begin{aligned}
uv^\prime &= - \frac{1}{x} \\
u^\prime v &= \frac{1}{x} \\
W &= - \frac{2}{x}
\end{aligned}
$$

Abel の式を逆に使うと、共通の方程式の一階項の係数が求まる。$W$ が零にならない区間で計算する。

$$
-W^\prime/W=\frac{1}{x}
$$
:::

:::example[対数を含む対]
重根の階数低下で得た定数と対数を組にする。

先に関数と微分を並べる。微分する順と関数の順を固定しておくと符号が決まる。

$$
\begin{aligned}
u &= 1 \\
v &= \log{\left(x \right)} \\
u^\prime &= 0 \\
v^\prime &= \frac{1}{x}
\end{aligned}
$$

積 $uv^\prime$ から $u^\prime v$ を引いて簡約する。

$$
\begin{aligned}
uv^\prime &= \frac{1}{x} \\
u^\prime v &= 0 \\
W &= \frac{1}{x}
\end{aligned}
$$

Abel の式を逆に使うと、共通の方程式の一階項の係数が求まる。$W$ が零にならない区間で計算する。

$$
-W^\prime/W=\frac{1}{x}
$$
:::

### 非斉次の零初期応答

$u,v$ を独立な同次解とすると、$x\ge s$ で $K(x,s)=[u(s)v(x)-v(s)u(x)]/W(s)$ を作ることができる。$K(s,s)=0$、$\partial_x K(s,s)=1$ であるため、$y(x)=\int_{x_0}^xK(x,s)f(s)\dd s$ を微分すると、一回目には端項がなく、二回目に $f(x)$ が一つ出る。積分部は各 $s$ について同次方程式を満たす。従って $y$ は正規化された二階式の零初期応答である。初期点をずらすと、積分下限の違いは二つの同次解の結合として現れる。
:::example[定数の源を加える]
零初期条件で $y^{\prime\prime}+4y=1$ を解く。積分核の正弦をまとめた式をもとの方程式と比較する。

係数を二階、一階、零階の順に並べる。

$$
\left( 1, \  0, \  4\right)
$$

$$
\begin{aligned}u&=\frac{1}{4} - \frac{\cos{\left(2 x \right)}}{4}\\u^{\prime}&=\frac{\sin{\left(2 x \right)}}{2}\\u^{\prime\prime}&=\cos{\left(2 x \right)}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=\cos{\left(2 x \right)}+\left(0\right)+\left(1 - \cos{\left(2 x \right)}\right)\\&=1\end{aligned}
$$

$x$ が $0$ のところで値を評価する。

$$
\begin{aligned}u&=0\\ u^\prime&=0\\ Lu&=1\end{aligned}
$$
:::
:::example[自由解と同じ周波数の源]
同じ演算子で $f(x)=\cos(2x)$ とすると、一次の振幅因子が現れる。

係数を二階、一階、零階の順に並べる。

$$
\left( 1, \  0, \  4\right)
$$

$$
\begin{aligned}u&=\frac{x \sin{\left(2 x \right)}}{4}\\u^{\prime}&=\frac{x \cos{\left(2 x \right)}}{2} + \frac{\sin{\left(2 x \right)}}{4}\\u^{\prime\prime}&=- x \sin{\left(2 x \right)} + \cos{\left(2 x \right)}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=- x \sin{\left(2 x \right)} + \cos{\left(2 x \right)}+\left(0\right)+\left(x \sin{\left(2 x \right)}\right)\\&=\cos{\left(2 x \right)}\end{aligned}
$$

$x$ が $0$ のところで値を評価する。

$$
\begin{aligned}u&=0\\ u^\prime&=0\\ Lu&=1\end{aligned}
$$
:::

### 一次系の行列式

一般の $d$ 成分一次系 $Y'=A(x)Y$ で、$d$ 本の解を列にした行列を $\Phi$ とする。$\det\Phi$ の微分は、各行を一行ずつ微分した行列式の和になる。第 $i$ 行の微分は $\sum_j A_{ij}\Phi_{j,:}$ であり、$j\ne i$ の項は同じ行が二つ現れて消える。このため、$\Phi$ が可逆かどうかにかかわらず $(\det\Phi)'=(\tr A)\det\Phi$ である。二階正規形では行列の trace が $-P$ なので Abel の式が戻る。この証明では $\Phi^{-1}$ を使わないため、従属な解を並べた場合にも成り立つ。

## 演習問題

::::exercise{level=basic}
$y''+4y=0$、$y(0)=1,y'(0)=0$ の最初の三つの偶数係数を求めよ。
:::solution
$a_{n+2}=-4a_n/[(n+1)(n+2)]$ である。$a_0=1,a_2=-2,a_4=2/3$ となり、奇数係数は全て零である。
:::
::::

::::exercise{level=basic}
$x^2y''+3xy'-3y=0$ の指数と解を求めよ。
:::solution
$I=r(r-1)+3r-3=(r-1)(r+3)$ である。係数に高次項が無いので独立解は $x,x^{-3}$ である。
:::
::::

::::exercise{level=basic}
$y''+x^{-1}y'=0$ の無限遠での指数を求めよ。
:::solution
$w=1/x$ にすると係数は $2/w-(w/w^2)=1/w$ で、同じ型になる。従って決定多項式は $r^2$ で指数は重根 $0$ である。
:::
::::

::::exercise{level=basic}
$(\cos 3x,\sin 3x)$ の Wronskian を求めよ。
:::solution
$W=3\cos^2 3x+3\sin^2 3x=3$ である。$P=0$ の Abel の公式とも一致する。
:::
::::

::::exercise{level=advanced}
$x^2y''+xy'-y=0$ で $y_1=x$ として階数低下を実行せよ。
:::solution
$e^{-\int P}=1/x$、分母 $y_1^2=x^2$ なので $v' = x^{-3}$ と取れる。$v=-1/(2x^2)$ より $y_2=-1/(2x)$ となり、$W=x^{-1}$ で独立である。
:::
::::

::::exercise{level=advanced}
$x^2y''+x(1+x)y'=0$ の重根で第 2 解の対数と一次項を求めよ。
:::solution
$y_1=1$、$I=r^2$ である。$v'=e^{-x}/x=x^{-1}-1+x/2+\cdots$ を項別積分すると $y_2=\log x-x+O(x^2)$ となる。
:::
::::

::::exercise{level=advanced}
$P=A/x+B/(x-1)$、$Q=0$ で $A=1/2,B=1/3$ の全指数を求めよ。
:::solution
有限点の根は $(0,1/2),(0,2/3)$ である。無限遠の根は積が零、和は $A+B-1=-1/6$ なので $(0,-1/6)$ である。和は $1/2+2/3-1/6=1$ となる。
:::
::::

::::exercise{level=advanced}
$x^2y'+y=x$ の上の実解に対し、$x-x^2+2x^3$ までの誤差を評価せよ。
:::solution
$x/(1+xt)=x-x^2t+x^3t^2-x^4t^3/(1+xt)$ なので、積分の剰余は $x>0$ で負、絶対値は $x^4\int_0^\infty e^{-t}t^3\dd t=6x^4$ 以下である。
:::
::::

## 確認テスト

::::quiz{answer="3" tol="1e-9"}
$y''=0$ の解 $1,3x$ の Wronskian を求めよ。

:::explanation
第一の導関数は零、第二は $3$ なので $W=3$ である。
:::
::::

::::quiz
二つの指数の差が正整数であるとき何が言えるか。

- [ ] 必ず対数が現れる
- [x] 小さい根の漸化式で分母が消える段がある
- [ ] 大きい根の級数は作れない

:::explanation
小さい根にその整数を足すと大きい根になるため、式 :eqref[eq-mp-frobenius-recurrence] の分母がその段で零になる。分子も零なら対数は不要な場合がある。
:::
::::

::::quiz{answer="-2" tol="1e-9"}
$x^2y''-6y=0$ の小さい指数を求めよ。

:::explanation
$I=r(r-1)-6=(r-3)(r+2)$ なので $-2$ である。
:::
::::

::::quiz
通常の正則点の Taylor 解を指定するデータはどれか。

- [x] $y$ と $y'$ のその点での値
- [ ] $y$ と $y''$ の任意の値
- [ ] 値を指定できない

:::explanation
級数の定数と一次係数は自由で、残りは二階方程式から順に決まる。二次係数は既にこの二つから制約される。
:::
::::

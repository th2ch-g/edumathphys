---
title: 留数解析と定積分の計算
description: 実積分を複素平面上の閉積分に移し、極、分枝、経路の向きと極限操作を照合して値を決める。
notation:
  - symbol: '\operatorname{PV}'
    meaning: Cauchy の主値
  - symbol: 'C_R'
    meaning: 半径 R の補助円弧
  - symbol: '\operatorname{Res}'
    meaning: 留数
---

## この章の目標

- 実軸、単位円、扇形、鍵穴型積分路を被積分関数に応じて選ぶ。
- 補助円弧が消える条件を長さと指数減衰から導く。
- 通常の広義積分と主値を区別し、半留数の符号を計算する。
- 積分核の極を整数に配置して無限級数を評価する。

## 三角関数の有理式の積分

$z=e^{it}$ とすると実変数の周期 $2\pi$ は単位円の正の一周になる。$\dd t=\dd z/(iz)$、$\cos t=(z+z^{-1})/2$、$\sin t=(z-z^{-1})/(2i)$ なので三角関数の有理式は $z$ の有理式に変わる。単位円上の分母の零点は、約分後も残るかを先に判定する。
:::proposition[余弦を含む一次分母]{#prop-cpx-cosine-integral}
実数 $a>|b|$ に対し $\int_0^{2\pi}\dd t/(a+b\cos t)=2\pi/\sqrt{a^2-b^2}$ である。
:::

:::proof
$b=0$ なら直接積分できる。$b\ne0$ では積分は $(2/i)\int_{|z|=1}\dd z/(bz^2+2az+b)$ となる。根は $r=(-a+\sqrt{a^2-b^2})/b=-b/(a+\sqrt{a^2-b^2})$ と $1/r$。従って $|r|<1$ であり、その根での分母の微分は $2\sqrt{a^2-b^2}$。:ref[thm-cpx-residue-theorem] を用いて結論を得る。
:::

:::example[分母の二乗]
$a=5,b=3$ なら上式は $\pi/2$ である。$a$ の小近傍で分母に一様な正の下限があるので $a$ で微分できる。負号を取り除けば $\int_0^{2\pi}(a+b\cos t)^{-2}\dd t=2\pi a/(a^2-b^2)^{3/2}$、同じ数値で $5\pi/32$ である。
:::

:::example[Fourier 係数]
$b\ne0$、整数 $n\ge0$ で $e^{int}$ を分子に加えると変換後に $z^n$ が掛かる。内点の留数が $r^n$ 倍になるから $\int_0^{2\pi}\cos(nt)/(a+b\cos t)\dd t=2\pi r^n/\sqrt{a^2-b^2}$。虚部は $t\mapsto2\pi-t$ による対称性で消える。
:::

### 正弦を含む分母と原点の極

$a+b\cos t+c\sin t$ を扱うとき、まず $d=\sqrt{b^2+c^2}$ と置く。$d>0$ なら $\cos\phi=b/d,\sin\phi=c/d$ を満たす実数 $\phi$ を選べる。すると分母は $a+d\cos(t-\phi)$ となる。一周積分では出発点を変えても値は同じだが、分子が $t$ を含めばそちらにも角の移動が要る。特に $a>d$ なら逆数の積分は $2\pi/\sqrt{a^2-d^2}$ である。$a<-d$ の場合は分母全体から負号をくくる。結果は同じ平方根の前に負号を付けたものになる。

:::example[混合項を含む具体例]
$\int_0^{2\pi}\dd t/(13+3\cos t+4\sin t)$ では $d=5$ で、分母の最小値は $8$。従って積分は $2\pi/\sqrt{169-25}=\pi/6$ である。次に同じ分母の二乗に対して分子に $\sin t$ を置く。$c$ に関する微分を使えば、$-\int\sin t\,(a+b\cos t+c\sin t)^{-2}\dd t=2\pi c(a^2-b^2-c^2)^{-3/2}$。指定値で積分は $-8\pi/12^3=-\pi/216$ である。分子が正の半周と負の半周で分母も変わるため、正弦が入っても零にならない。
:::

:::example[$z=0$ の見落としを防ぐ]
$a>1$ で $\int_0^{2\pi}\cos t/(a-\cos t)\dd t$ を計算する。実変数で先に割れば $a\int\dd t/(a-\cos t)-2\pi=2\pi(a/\sqrt{a^2-1}-1)$ となる。一方 $z=e^{it}$ を直接使うと、有理型の微分は
$$
-\frac{z^2+1}{iz(z^2-2az+1)}\,\dd z
$$
であり、二次式の内根 $a-\sqrt{a^2-1}$ のほかに原点にも単純極がある。原点の留数は $i$ で、積分に $-2\pi$ を加える。三角式の分母の根だけを数えた結果と実変数の式との差は、この測度と分子に由来する極である。
:::

:::example[境界に触れる極]
$a>1$ から $a\downarrow1$ とすると $\int_0^{2\pi}(a-\cos t)^{-1}\dd t=2\pi/\sqrt{(a-1)(a+1)}$ は発散する。円内の極と円外の極が $z=1$ に合流するからである。実変数でも $1-\cos t\sim t^2/2$ により、$a=1$ での被積分関数は $2/t^2$ に等価となる。一方分子が $\sin^2t$ なら $a=1$ で $(1-\cos^2t)/(1-\cos t)=1+\cos t$ まで約せるので、積分は $2\pi$ である。境界の点が分母を消すときは、約分と局所展開を先に行う。
:::

:::example[分子の一次の余弦がもたらす補正]
分母は $[1,3]$ に収まり、原点の極の寄与も残る。

$$
I=\int_0^{2\pi}\frac{\cos{\left(t \right)}}{2 - \cos{\left(t \right)}}\,\dd t.
$$

単位円への置換で次式を得る。

$$
I=\int_{|z|=1} F(z)\,\dd z,\qquad F(z)=\frac{i \left(z - i\right) \left(z + i\right)}{z \left(z^{2} - 4 z + 1\right)}.
$$

各極 $a$ で $(z-a)^mF(z)$ の Taylor 展開を作る。以下では左に円内の極と位数、右に第 $m-1$ 次係数である留数を示す。

$$
\begin{aligned}
z=2 - \sqrt{3}\quad(m=1) &\quad - \frac{2 \sqrt{3} i}{3}\\z=0\quad(m=1) &\quad i
\end{aligned}
$$

従って $I=2\pi i(\frac{i \left(3 - 2 \sqrt{3}\right)}{3})=\frac{2 \pi \left(-3 + 2 \sqrt{3}\right)}{3}$ である。分子を $2-(2-\cos t)$ に分けると、同じ答えを逆数の平均と定数の差として得る。
:::

:::example[奇数の調波を伴う積分]
分子は $\cos(3t)$ である。角度を三倍するため、単位円置換の分子にも負の次数が現れる。

$$
I=\int_0^{2\pi}\frac{4 \cos^{3}{\left(t \right)} - 3 \cos{\left(t \right)}}{5 - 4 \cos{\left(t \right)}}\,\dd t.
$$

単位円への置換で次式を得る。

$$
I=\int_{|z|=1} F(z)\,\dd z,\qquad F(z)=\frac{i \left(z - i\right) \left(z + i\right) \left(z^{2} - i z - 1\right) \left(z^{2} + i z - 1\right)}{4 z^{3} \left(z - 2\right) \left(z - \frac{1}{2}\right)}.
$$

単位円内での極と留数は次の組になる。

$$
\begin{aligned}
z=\frac{1}{2}\quad(m=1) &\quad - \frac{65 i}{48}\\z=0\quad(m=3) &\quad \frac{21 i}{16}
\end{aligned}
$$

従って $I=2\pi i(- \frac{i}{24})=\frac{\pi}{12}$ である。分母の逆数の第三 Fourier 係数を測ったことになる。ここでの正の値は分母が小さくなる $t=0$ 付近の寄与と整合する。
:::

:::example[正弦の二乗で定まる重み]
$\sin^2t$ は非負なので答えの符号を積分前に予測できる。

$$
I=\int_0^{2\pi}\frac{\sin^{2}{\left(t \right)}}{\left(\cos{\left(t \right)} + 2\right)^{2}}\,\dd t.
$$

単位円への置換で次式を得る。

$$
I=\int_{|z|=1} F(z)\,\dd z,\qquad F(z)=\frac{i \left(z - 1\right)^{2} \left(z + 1\right)^{2}}{z \left(z^{2} + 4 z + 1\right)^{2}}.
$$

単位円内での極と留数は次の組になる。

$$
\begin{aligned}
z=-2 + \sqrt{3}\quad(m=2) &\quad - \frac{2 \sqrt{3} i}{3}\\z=0\quad(m=1) &\quad i
\end{aligned}
$$

従って $I=2\pi i(\frac{i \left(3 - 2 \sqrt{3}\right)}{3})=\frac{2 \pi \left(-3 + 2 \sqrt{3}\right)}{3}$ である。実変数では $\sin^2t=1-\cos^2t$ として分母の逆べきに分解できる。
:::

:::example[三角式の約分で消える境界零点]
分母と分子は $t=0$ で同時に零になる。そこは連続延長した値で積分する。

$$
I=\int_0^{2\pi}\frac{1 - \cos^{2}{\left(t \right)}}{1 - \cos{\left(t \right)}}\,\dd t.
$$

単位円への置換で次式を得る。

$$
I=\int_{|z|=1} F(z)\,\dd z,\qquad F(z)=- \frac{i \left(z + 1\right)^{2}}{2 z^{2}}.
$$

単位円内での極と留数は次の組になる。

$$
\begin{aligned}
z=0\quad(m=2) &\quad - i
\end{aligned}
$$

従って $I=2\pi i(- i)=2 \pi$ である。約分後は $1+\cos t$ である。一周の余弦の平均が零であることからも値を確認できる。
:::

:::example[余弦の二乗を含む正の分母]
$1+3\cos^2t$ は一周期の中で最小値を二回取る。二倍角で処理する方法と極の数を比較する。

$$
I=\int_0^{2\pi}\frac{1}{3 \cos^{2}{\left(t \right)} + 1}\,\dd t.
$$

単位円への置換で次式を得る。

$$
I=\int_{|z|=1} F(z)\,\dd z,\qquad F(z)=- \frac{4 i z}{3 \left(z^{2} + \frac{1}{3}\right) \left(z^{2} + 3\right)}.
$$

単位円内での極と留数は次の組になる。

$$
\begin{aligned}
z=- \frac{\sqrt{3} i}{3}\quad(m=1) &\quad - \frac{i}{4}\\z=\frac{\sqrt{3} i}{3}\quad(m=1) &\quad - \frac{i}{4}
\end{aligned}
$$

従って $I=2\pi i(- \frac{i}{2})=\pi$ である。実変数では $1+3\cos^2t=(5+3\cos 2t)/2$。$u=2t$ と置くと二周になるため、測度の $1/2$ と周期の個数が相殺する。
:::

:::example[二つの異なる因子をもつ円周平均]
二つの一次分母を持ち、単位円内には異なる二つの正の実極が現れる。

$$
I=\int_0^{2\pi}\frac{1}{\left(2 - \cos{\left(t \right)}\right) \left(3 - \cos{\left(t \right)}\right)}\,\dd t.
$$

単位円への置換で次式を得る。

$$
I=\int_{|z|=1} F(z)\,\dd z,\qquad F(z)=- \frac{4 i z}{\left(z^{2} - 6 z + 1\right) \left(z^{2} - 4 z + 1\right)}.
$$

単位円内での極と留数は次の組になる。

$$
\begin{aligned}
z=3 - 2 \sqrt{2}\quad(m=1) &\quad \frac{\sqrt{2} i}{4}\\z=2 - \sqrt{3}\quad(m=1) &\quad - \frac{\sqrt{3} i}{3}
\end{aligned}
$$

従って $I=2\pi i(i \left(- \frac{\sqrt{3}}{3} + \frac{\sqrt{2}}{4}\right))=\frac{\pi \left(- 3 \sqrt{2} + 4 \sqrt{3}\right)}{6}$ である。実変数でも $(2-C)^{-1}-(3-C)^{-1}$ に分解できる。二因子の積は正なので差の符号を確認できる。
:::

:::example[平行移動した角の効果]
正弦の一次式が分母にあるときは有理式の係数に虚数が現れる。

$$
I=\int_0^{2\pi}\frac{1}{3 \sin{\left(t \right)} + 5}\,\dd t.
$$

単位円への置換で次式を得る。

$$
I=\int_{|z|=1} F(z)\,\dd z,\qquad F(z)=\frac{2}{3 \left(z + \frac{i}{3}\right) \left(z + 3 i\right)}.
$$

単位円内での極と留数は次の組になる。

$$
\begin{aligned}
z=- \frac{i}{3}\quad(m=1) &\quad - \frac{i}{4}
\end{aligned}
$$

従って $I=2\pi i(- \frac{i}{4})=\frac{\pi}{2}$ である。角を $\pi/2$ 平行移動すると余弦の場合に一致する。積分値が実数になったことは、留数計算の向きを確かめる手掛かりである。
:::

## 有理関数の広義積分

:::proposition[有理関数を半円で閉じる]{#prop-cpx-rational-integral}
実軸上に極をもたない有理関数 $R=P/Q$ が $\deg Q\ge\deg P+2$ を満たせば、実軸全体で絶対可積分であり、その積分は上半平面にある留数の和の $2\pi i$ 倍である。
:::

:::proof
十分大きい $|z|$ で $|R(z)|\le C|z|^{-2}$。従って上半円の積分は絶対値で $\pi R_0 C/R_0^2\to0$ となる。ここで $R_0$ は半径である。実軸の両端でも $x^{-2}$ と比較して絶対収束する。半円と直径からなる正向き境界に留数定理を使い、半径を無限大へ送ればよい。
:::

:::example[四次式]
$\int_0^\infty\dd x/(1+x^4)$ を求める。上半平面の根を $u=e^{\pi i/4},v=e^{3\pi i/4}$ とすると各留数は $1/(4z^3)=-z/4$。和は $-i\sqrt2/4$、実軸積分は $\pi/\sqrt2$ である。偶関数なので半直線では $\pi/(2\sqrt2)$。
:::

:::example[重極を含む積分]
$a>0$ なら $(z^2+a^2)^{-2}$ の上半平面の極は $ia$ の二位の極である。留数は $[(z+ia)^{-2}]'_{z=ia}=-2/(2ia)^3=1/(4ia^3)$。従って $\int_{\R}(x^2+a^2)^{-2}\dd x=\pi/(2a^3)$。$\pi/a$ をパラメータで微分しても同じ結果となる。
:::

### 四次分母の実パラメータを残す

:::example[共役な極が重なっても値が続く積分]
$A>-1$ について
$$
I(A)=\int_0^\infty \frac{\dd x}{x^4+2Ax^2+1}
$$
を一つの式で求める。$b=\sqrt{2(A+1)}$ とすれば分母は $(z^2+ibz-1)(z^2-ibz-1)$。後者の根は $(ib\pm\sqrt{4-b^2})/2$ で、$b<2$ なら実部が逆符号の二つの上半平面の点である。$b>2$ なら二つとも正の虚軸上にある。分母の正負ではなく虚部の符号により閉路内の因子を判定できる。$b=2$、すなわち $A=1$ のときは二根が $i$ に重なる。

$Q_\pm=z^2\pm ibz-1$ と置くと
$$
\frac{1+z^2}{Q_+Q_-}=\frac1{2ib}
\left(\frac{Q_+'}{Q_+}-\frac{Q_-'}{Q_-}\right)(-1).
$$
実際、左の差の分子は $(2z+ib)(z^2-ibz-1)-(2z-ib)(z^2+ibz-1)=-2ib(z^2+1)$。従って上半平面における $(1+z^2)/(Q_+Q_-)$ の留数和は $2/(2ib)=1/(ib)$ で、二重根の場合も対数微分の留数を重複度で数えれば同じである。実軸全体での積分は $2\pi/b$。半直線で $x=1/t$ と変えると分子 $1$ と分子 $x^2$ の二つの積分は等しい。偶性も合わせれば $4I(A)=2\pi/b$、すなわち $I(A)=\pi/(2\sqrt{2(A+1)})$ である。$A=1$ では $\pi/4$、$A=0$ では $\pi/(2\sqrt2)$ が再現される。
:::

:::example[減衰次数が一つ不足するとき]
$z/(z^2+1)$ の上側の留数は $1/2$ である。この関数を直径と上半円上で積分しても、実軸の主値をそのまま $\pi i$ としてはならない。十分遠くで $z/(1+z^2)=1/z+O(|z|^{-3})$ だから半円の寄与は $i\pi+O(R^{-2})$ である。これを閉積分 $\pi i$ から引けば対称主値は零となる。通常の広義積分では正負の端が別々に対数発散する。
:::

### 局所展開を使った三位の極の計算

:::example[$\int_{\mathbb R}(1+x^2)^{-3}\dd x$]
上側の極は $i$ である。$h=z-i$ として
$$
(z+i)^{-3}=(2i)^{-3}
\left(1-\frac{3h}{2i}+\frac{6h^2}{(2i)^2}+O(h^3)\right)
$$
と書けば $h^{-3}$ を掛けた後の留数は $6/(2i)^5=3/(16i)$。積分は $3\pi/8$ である。別法では $J_n(a)=\int_{\R}(x^2+a^2)^{-n}\dd x$ と置き、$J_{n+1}(a)=-J_n'(a)/(2na)$ を用いる。$J_1(a)=\pi/a$ から $J_2(a)=\pi/(2a^3)$、$J_3(a)=3\pi/(8a^5)$ を得る。$a$ を正のコンパクト区間に限れば微分後も $(1+x^2)^{-n-1}$ の定数倍で抑えられるので、パラメータ微分は正当である。
:::

:::example[左右対称でない実極なしの分母]
頂点が $x=1$ にある Lorentz 型の関数である。正区間と負区間には異なる面積があるため、全体から半直線の値を半分にして求めることはできない。

$$
I=\int_{-\infty}^\infty \frac{1}{\left(x - 1\right)^{2} + 4}\,\dd x.
$$

分母と分子の次数差は $2$ である。半径 $R$ の半円上では、長さまで含めた積分の評価が $O(R^{-1})$ となって消える。中心 $a=1 + 2 i$ で $v=z-a$ とすれば、負次数までの局所式は

$$
F(a+v)=- \frac{i}{4 v}+O(1),\qquad \operatorname{Res}_{a}F=- \frac{i}{4}.
$$

正向きに閉じて $I=2\pi i\sum_{\Im a>0}\operatorname{Res}_{a}F=\frac{\pi}{2}$ を得る。実変数 $u=(x-1)/2$ で積分は $(1/2)\int_{\R}(1+u^2)^{-1}\dd u$ となる。
:::

:::example[正の山と負の山が打ち消す実積分]
分子は実軸上で符号を変える。正の関数を積分したときの直観で答えの符号を決めない。

$$
I=\int_{-\infty}^\infty \frac{x^{2} - 1}{\left(x^{2} + 1\right)^{2}}\,\dd x.
$$

分母と分子の次数差は $2$ である。半径 $R$ の半円上では、長さまで含めた積分の評価が $O(R^{-1})$ となって消える。中心 $a=i$ で $v=z-a$ とすれば、負次数までの局所式は

$$
F(a+v)=\frac{1}{2 v^{2}}+O(1),\qquad \operatorname{Res}_{a}F=0.
$$

正向きに閉じて $I=2\pi i\sum_{\Im a>0}\operatorname{Res}_{a}F=0$ を得る。実変数の原始関数 $-x/(1+x^2)$ は両端で零に戻る。上半平面の二位の極でも留数が零となりうる。
:::

:::example[平行移動の異なる二つの山]
二つの分母因子の中心は $0,2$。それぞれの極は同じ高さにあるが、個別の留数は純虚数に限らない。

$$
I=\int_{-\infty}^\infty \frac{1}{\left(x^{2} + 1\right) \left(\left(x - 2\right)^{2} + 1\right)}\,\dd x.
$$

分母と分子の次数差は $4$ である。半径 $R$ の半円上では、長さまで含めた積分の評価が $O(R^{-3})$ となって消える。中心 $a=2 + i$ で $v=z-a$ とすれば、負次数までの局所式は

$$
F(a+v)=\frac{- \frac{1}{16} - \frac{i}{16}}{v}+O(1),\qquad \operatorname{Res}_{a}F=- \frac{1}{16} - \frac{i}{16}.
$$

中心 $a=i$ で $v=z-a$ とすれば、負次数までの局所式は

$$
F(a+v)=\frac{\frac{1}{16} - \frac{i}{16}}{v}+O(1),\qquad \operatorname{Res}_{a}F=\frac{1}{16} - \frac{i}{16}.
$$

正向きに閉じて $I=2\pi i\sum_{\Im a>0}\operatorname{Res}_{a}F=\frac{\pi}{4}$ を得る。変換 $x\mapsto2-x$ で実関数は不変であり、二つの極の寄与の実部が相殺する。
:::

:::example[二つの異なる減衰長]
分子の $x^2$ が原点付近の寄与を弱める。両端では $x^{-2}$ に等価なので絶対収束する。

$$
I=\int_{-\infty}^\infty \frac{x^{2}}{\left(x^{2} + 1\right) \left(x^{2} + 9\right)}\,\dd x.
$$

分母と分子の次数差は $2$ である。半径 $R$ の半円上では、長さまで含めた積分の評価が $O(R^{-1})$ となって消える。中心 $a=i$ で $v=z-a$ とすれば、負次数までの局所式は

$$
F(a+v)=\frac{i}{16 v}+O(1),\qquad \operatorname{Res}_{a}F=\frac{i}{16}.
$$

中心 $a=3 i$ で $v=z-a$ とすれば、負次数までの局所式は

$$
F(a+v)=- \frac{3 i}{16 v}+O(1),\qquad \operatorname{Res}_{a}F=- \frac{3 i}{16}.
$$

正向きに閉じて $I=2\pi i\sum_{\Im a>0}\operatorname{Res}_{a}F=\frac{\pi}{4}$ を得る。分子を比較すると部分分数の係数は $-1/8$ と $9/8$。二つの Lorentz 積分にもどせる。
:::

:::example[円周上に並ぶ四つの単純極]
四次方程式の根は偏角が $\pi/4$ ずつではなく、最初が $\pi/4$ でそこから $\pi/2$ おきに現れる。

$$
I=\int_{-\infty}^\infty \frac{1}{x^{4} + 1}\,\dd x.
$$

分母と分子の次数差は $4$ である。半径 $R$ の半円上では、長さまで含めた積分の評価が $O(R^{-3})$ となって消える。中心 $a=- \frac{\sqrt{2}}{2} + \frac{\sqrt{2} i}{2}$ で $v=z-a$ とすれば、負次数までの局所式は

$$
F(a+v)=\frac{\sqrt{2} \left(1 - i\right)}{8 v}+O(1),\qquad \operatorname{Res}_{a}F=\frac{\sqrt{2} \left(1 - i\right)}{8}.
$$

中心 $a=\frac{\sqrt{2}}{2} + \frac{\sqrt{2} i}{2}$ で $v=z-a$ とすれば、負次数までの局所式は

$$
F(a+v)=\frac{\sqrt{2}}{2 v \left(1 + i\right)^{3}}+O(1),\qquad \operatorname{Res}_{a}F=\frac{\sqrt{2} \left(-1 - i\right)}{8}.
$$

正向きに閉じて $I=2\pi i\sum_{\Im a>0}\operatorname{Res}_{a}F=\frac{\sqrt{2} \pi}{2}$ を得る。同じ積分は四分円を用いて実軸と虚軸を比較しても得られる。
:::

:::example[六次分母と四次分子]
正の半直線で $x\mapsto1/x$ とすれば、分子 $x^4$ が分子 $1$ に移る。上半円は三つの極を含む。

$$
I=\int_{-\infty}^\infty \frac{x^{4}}{x^{6} + 1}\,\dd x.
$$

分母と分子の次数差は $2$ である。半径 $R$ の半円上では、長さまで含めた積分の評価が $O(R^{-1})$ となって消える。中心 $a=- \frac{\sqrt{3}}{2} + \frac{i}{2}$ で $v=z-a$ とすれば、負次数までの局所式は

$$
F(a+v)=\frac{\left(\sqrt{3} - i\right)^{4}}{48 v \left(\sqrt{3} + i\right)}+O(1),\qquad \operatorname{Res}_{a}F=- \frac{\sqrt{3}}{12} - \frac{i}{12}.
$$

中心 $a=i$ で $v=z-a$ とすれば、負次数までの局所式は

$$
F(a+v)=- \frac{i}{6 v}+O(1),\qquad \operatorname{Res}_{a}F=- \frac{i}{6}.
$$

中心 $a=\frac{\sqrt{3}}{2} + \frac{i}{2}$ で $v=z-a$ とすれば、負次数までの局所式は

$$
F(a+v)=\frac{1 - \sqrt{3} i}{6 v \left(\sqrt{3} - i\right)}+O(1),\qquad \operatorname{Res}_{a}F=\frac{\sqrt{3}}{12} - \frac{i}{12}.
$$

正向きに閉じて $I=2\pi i\sum_{\Im a>0}\operatorname{Res}_{a}F=\frac{2 \pi}{3}$ を得る。単位円上の三つの上側の根のうち、虚軸上の根の寄与を忘れると正の実積分の値が変わる。
:::

:::example[三位の極を伴う偶モーメント]
同じ三位の極でも $x^2$ による重みが留数を変える。局所座標に移したときの分子も展開する必要がある。

$$
I=\int_{-\infty}^\infty \frac{x^{2}}{\left(x^{2} + 1\right)^{3}}\,\dd x.
$$

分母と分子の次数差は $4$ である。半径 $R$ の半円上では、長さまで含めた積分の評価が $O(R^{-3})$ となって消える。中心 $a=i$ で $v=z-a$ とすれば、負次数までの局所式は

$$
F(a+v)=- \frac{i}{16 v} - \frac{1}{16 v^{2}} - \frac{i}{8 v^{3}}+O(1),\qquad \operatorname{Res}_{a}F=- \frac{i}{16}.
$$

正向きに閉じて $I=2\pi i\sum_{\Im a>0}\operatorname{Res}_{a}F=\frac{\pi}{8}$ を得る。実軸上で $x^2=(1+x^2)-1$ と分ければ、二位と三位の分母の既知の積分の差でも検算できる。
:::

:::example[偶性だけでは消せない一次の分子]
実軸の原点について奇関数ではないが、中心をずらすと一方の山がもう一方を打ち消す。

$$
I=\int_{-\infty}^\infty \frac{x + 1}{\left(\left(x + 1\right)^{2} + 1\right)^{2}}\,\dd x.
$$

分母と分子の次数差は $3$ である。半径 $R$ の半円上では、長さまで含めた積分の評価が $O(R^{-2})$ となって消える。中心 $a=-1 + i$ で $v=z-a$ とすれば、負次数までの局所式は

$$
F(a+v)=- \frac{i}{4 v^{2}}+O(1),\qquad \operatorname{Res}_{a}F=0.
$$

正向きに閉じて $I=2\pi i\sum_{\Im a>0}\operatorname{Res}_{a}F=0$ を得る。実際、元の関数は $-1/[2((x+1)^2+1)]$ の導関数で、全積分は両端値の差である。
:::

### 無限遠留数による符号の独立検算

:::example[非対称な実分母でも留数の実部が消える理由]
実係数有理関数 $R$ が実軸上で正則、無限遠で $O(|z|^{-2})$ とする。非実極 $a$ の Laurent 展開を共役すれば、$\bar a$ は同じ位数の極で $\operatorname{Res}_{\bar a}R=\overline{\operatorname{Res}_{a}R}$ である。上側の留数和を $S$ と置けば下側の和は $\bar S$。無限遠の局所座標では $-R(1/w)/w^2=O(1)$ であり、その留数は零である。球面上の総留数の消失により $S+\bar S=0$、従って $S$ は純虚数である。実積分 $2\pi iS$ が必ず実数となる理由が分かる。

例えば $R(z)=1/[(z^2+1)((z-2)^2+1)]$ の上側の二点では、留数の実部は個別には消えない。$i$ での留数は $1/[2i(4-4i)]=(1-i)/16$、$2+i$ では $1/[2i(4+4i)]=(-1-i)/16$ である。従って $S=-i/8$、積分は $\pi/4$。もし二番目の分母の符号を間違えると $S$ に実部が残り、この独立な条件を破る。この検算は「答えが実数」の条件までであり、純虚数部分の係数まで正しいことは各極での計算によって確かめる必要がある。
:::

## Jordan の補題と Fourier 型積分

:::lemma[Jordan の補題]{#lem-cpx-jordan}
$b>0$、$g$ が十分大きい半径の上半円上で連続で、$M_R=\sup_{0\le\theta\le\pi}|g(Re^{i\theta})|\to0$ なら、上半円上で $\int e^{ibz}g(z)\dd z\to0$ である。
:::

:::proof
指数の絶対値は $e^{-bR\sin\theta}$ である。$0\le\theta\le\pi/2$ で $\sin\theta\ge2\theta/\pi$、後半は前半と対称なので、積分の絶対値は $2RM_R\int_0^{\pi/2}e^{-2bR\theta/\pi}\dd\theta\le\pi M_R/b$。$b<0$ では下半円を使い、閉曲線の向きが負になる。
:::

:::example[Lorentz 型の Fourier 積分]
$a>0$、$b>0$ として $e^{ibz}/(z^2+a^2)$ の上側の留数は $e^{-ab}/(2ia)$。従って $\int_{\R}e^{ibx}/(x^2+a^2)\dd x=(\pi/a)e^{-a|b|}$ となる。負の $b$ では下側の $-ia$ と負向きの符号が相殺し、$b=0$ では有理関数の積分である。虚部は奇関数で零である。
:::

:::example[絶対収束しない Fourier 積分]
$a,b>0$ で $x\sin(bx)/(x^2+a^2)$ の広義積分は両端で存在する。有理因子が各端で単調に零に近づき、振動因子 $e^{ibx}$ の原始関数が有界なので、部分積分で複素積分の片側ごとの収束も従う。一方、正の実軸で $|\sin(bx)|\ge1/2$ となる一定長の区間を周期ごとに取ると、絶対値の積分は調和級数との比較で発散する。$ze^{ibz}/(z^2+a^2)$ の留数は $e^{-ab}/2$。虚部を取ると全実軸の積分は $\pi e^{-ab}$ である。
:::

### 因果的な Fourier 積分と周波数零

:::example[極が片側だけにある場合]
$a>0$ として $F(t)=\int_{\R} e^{itx}/(a+ix)^2\,\dd x$ を考える。分母の零点は $ia$ であり、符号は $a+ix=i(x-ia)$ から確認できる。$t>0$ では上側で閉じ、留数は $[-e^{itz}]'_{z=ia}=-it e^{-at}$ であるから $F(t)=2\pi t e^{-at}$。$t<0$ では下側に極はなく $F(t)=0$。$t=0$ では二位の極の留数は零で、やはり積分は零である。元の関数は絶対可積分なので $F$ の連続性は優収束からも従う。結果が原点で連続に接続することは分母の次数と一致する。
:::

:::example[二位の極と周波数の多項式]
$a>0,t>0$ で $e^{itz}/(z^2+a^2)^2$ の留数は
$$
\left.\frac{\dd}{\dd z}\left(\frac{e^{itz}}{(z+ia)^2}\right)\right|_{z=ia}
=e^{-at}\left(\frac{it}{(2ia)^2}-\frac{2}{(2ia)^3}\right)
=\frac{e^{-at}(1+at)}{4ia^3}.
$$
負の $t$ は共役対称性から求まり、実軸全体での余弦積分は $\pi(1+a|t|)e^{-a|t|}/(2a^3)$ である。$t=0$ の値と前節を照合できる。パラメータ $a$ について $\pi e^{-a|t|}/a$ を微分して $-1/(2a)$ 倍しても同じ結果になる。
:::

:::example[直線を横にずらす]
実数 $c$ を用い $[(x-c)^2+a^2]^{-1}$ の Fourier 積分を求めると、実変数の平行移動から前節の値に $e^{itc}$ が掛かる。留数計算でも上側の極は $c+ia$、そこでの指数は $e^{itc}e^{-at}$ となる。例えば $c=\pi,\ t=1,\ a=2$ では余弦積分は $-\pi e^{-2}/2$ である。留数から読み取れる二つの要素のうち、極の高さは減衰率を、実部は振動の位相を決めている。
:::

### 長方形の上下の辺を関係付ける

半円の上で評価しにくい関数でも、虚方向の周期、または反周期があれば積分の上辺を下辺と比較できる。例えば $\cosh(z+i\pi)=-\cosh z$ である。この変換の下では、分母の符号と Fourier 因子の倍率 $e^{-\pi t}$ が同時に変わる。長方形を使えば上辺を無限遠に追いやる必要はない。

:::example[双曲線正割の Fourier 変換]
実数 $t$ に対して
$$
H(t)=\int_{\R}\frac{e^{itx}}{\cosh x}\,\dd x
$$
を考える。$\cosh x\ge e^{|x|}/2$ により絶対収束し、$t$ は正負どちらでもよい。底辺 $[-R,R]$、高さ $\pi$ の正向き長方形上で $F(z)=e^{itz}/\cosh z$ を積分する。極は内部の $z=i\pi/2$ 一つである。$\sinh(i\pi/2)=i$ なので留数は $e^{-\pi t/2}/i$。上辺では $F(x+i\pi)=-e^{-\pi t}F(x)$ で、積分方向も右から左へ変わるから、下辺との和は $(1+e^{-\pi t})\int_{-R}^R F(x)\dd x$ となる。

垂直辺を落とせることを確かめる。$|\cosh(R+iy)|^2=\sinh^2R+\cos^2y\ge\sinh^2R$、分子の絶対値は $e^{-ty}\le e^{\pi|t|}$ である。二本の辺の長さを掛けた上界は $2\pi e^{\pi|t|}/\sinh R\to0$。従って
$$
(1+e^{-\pi t})H(t)=2\pi e^{-\pi t/2},\qquad
H(t)=\frac{\pi}{\cosh(\pi t/2)} .
$$
答えが偶関数かつ実数となるのは、元の実数値の重みが偶関数で正弦成分の積分が零であるためである。$t=0$ なら $u=e^x$ と置いて $H(0)=2\int_0^\infty(1+u^2)^{-1}\dd u=\pi$ とも確かめられる。
:::

:::example[二乗した双曲線正割の変換]
$K(t)=\int_{\R} e^{itx}\cosh^{-2}x\,\dd x$ では高さ $i\pi$ の移動で分母の符号変化が消える。そのため上下辺が与える係数は今度は $1-e^{-\pi t}$ である。$z=i\pi/2+h$ では $\cosh z=i\sinh h=ih(1+h^2/6+\cdots)$ なので
$$
\frac{e^{itz}}{\cosh^2z}=-e^{-\pi t/2}h^{-2}(1+ith+O(h^2)).
$$
留数は $-it e^{-\pi t/2}$。垂直辺は前例より速く減衰するから、$t\ne0$ で
$$
K(t)=\frac{2\pi t e^{-\pi t/2}}{1-e^{-\pi t}}
=\frac{\pi t}{\sinh(\pi t/2)} .
$$
$t=0$ では上下辺が打ち消し合い、閉積分は $0=0$ を与えるだけである。ここでは $K$ の連続性と $\sinh u/u\to1$ を用い $K(0)=2$ となる。直接には $\tanh x$ の両端値の差でもある。変換式で割る係数が零になる点を、無条件に元の等式へ代入してはならない。
:::

## 多価関数を含む積分と鍵穴型積分路

:::theorem[鍵穴型積分]{#thm-cpx-keyhole}
$0<p<1$ のとき $\int_0^\infty x^{p-1}/(1+x)\dd x=\pi/\sin(\pi p)$ である。
:::

:::proof
正の実軸を切り $0<\arg z<2\pi$ を選ぶ。$z^{p-1}/(1+z)$ の積分路は切れ目の上側を $\eps$ から $R$、外円を正向き、下側を $R$ から $\eps$、内円を負向きとする。厳密には切れ目と小角をなす二本の半直線で積分し、固定した両半径で角を零へ送る。上下の差は $(1-e^{2\pi ip})$ 倍の実積分になる。外円は $O(R^{p-1})$、内円は $O(\eps^p)$ なので消える。$-1$ の留数は $e^{i\pi(p-1)}$。従って求める値は $2\pi i e^{i\pi(p-1)}/(1-e^{2\pi ip})=\pi/\sin(\pi p)$ である。
:::

:::example[平方根と対数]
$p=1/2$ では $\int_0^\infty \dd x/(\sqrt{x}(1+x))=\pi$。$p$ が $0,1$ から正距離離れていれば $\log x$ を一つ、二つ掛けても両端で支配可能なので微分できる。$p=1/2$ で一次の対数モーメントは零、二次の対数モーメントは $\pi^3$ である。一次の零は $x\mapsto1/x$ でも確かめられる。
:::

:::example[尺度とべき]
$0<s<q,\ a>0$ とすると $u=x^q/a$ により $\int_0^\infty x^{s-1}/(a+x^q)\dd x=\pi a^{s/q-1}/[q\sin(\pi s/q)]$。例えば $s=2,q=6,a=1$ で $\pi/(3\sqrt3)$。端点 $s=q$ では尾部が $1/x$ に等価なので発散する。
:::

### 切れ目の跳びと対数の二乗

:::example[符号の変わる対数積分]
$\int_0^\infty \log x/(1+x)^2\,\dd x$ は絶対収束する。$x=1/u$ で実積分がちょうどその負になるから答えは零である。鍵穴積分では $L(z)=\log|z|+i\arg z$、$0<\arg z<2\pi$ を取る。この対数の上側と下側の差は一定の $2\pi i$ にとどまり、欲しい $\log x$ が残らない。そこで $L(z)^2/(1+z)^2$ を選ぶ。上下の和の分子は $(\log x)^2-(\log x+2\pi i)^2=4\pi^2-4\pi i\log x$。$-1$ での留数は $2L(-1)/(-1)=-2\pi i$。よって閉積分は $4\pi^2$ であり、$\int_0^\infty\dd x/(1+x)^2=1$ を使えば対数積分が零と分かる。外円の上界は $O((\log R)^2/R)$、内円の上界は $O(\eps|\log\eps|^2)$ なのでいずれも消える。
:::

:::example[一見消えない端点の差を取る]
$a,b>0$ について
$$
\int_0^\infty\frac{\log(x+a)-\log(x+b)}{(x+a)(x+b)}\dd x
$$
を求める。$a\ne b$ なら $y=\log[(x+a)/(x+b)]$ の微分は $(b-a)\dd x/[(x+a)(x+b)]$。$x=0$ で $y=\log(a/b)$、無限遠で $y=0$ と境界値を決めれば、結果は $(\log(a/b))^2/[2(a-b)]$ となる。複素平面で二つの対数の切れ目を考える前に、差の微分が選んだ有理因子に一致するかを調べると積分路を省ける。$a=b$ のときは元の積分が恒等的に零であり、上式の極限とも一致する。
:::

### Mellin 積分の帯と対数モーメント

$q>0$ とし、複素パラメータの冪を $x^{s-1}=\exp((s-1)\log x)$ とする。分母が $1+x^q$ の場合、$x<1$ での上界と $x>1$ での上界を別々に取ると絶対収束する帯は $0<\Re s<q$ である。帯の閉部分では両端で指数に余裕が残るため、任意の固定回数だけ $s$ で微分しても支配関数を置ける。$|\log x|^j$ がどの小さい正のべき $x^{-\delta}$ よりも遅く原点で増大すること、無限遠でも $x^\delta$ より遅いことを使えばよい。従って以下の計算では、一度鍵穴型の積分を求めれば、その同じ導出から複数の実積分を得られる。

$a>0$ を入れて $M(s)=\int_0^\infty x^{s-1}/(a+x^q)\dd x$ とする。$u=x^q/a$ と置換すれば、まず実数 $0<s<q$ で次の $M(s)$ の式を得る。両辺は帯 $0<\Re s<q$ で正則なので、一致の定理により帯全体で成り立つ。後の二式は対数微分から従う。

$$
M(s)=\frac{\pi a^{s/q-1}}{q\sin(\pi s/q)},\quad
\frac{M'}{M}=\frac{\log a-\pi\cot(\pi s/q)}q,\quad
\left(\frac{M'}{M}\right)'=\frac{\pi^2}{q^2\sin^2(\pi s/q)}.
$$
実積分で一次の対数が正負の値を取り、二次の対数が常に非負になることを、これらの式から点検する。$\log a$ は積分変数を拡大する際に増えた定数であり、$a=1$ の公式に単に倍率を掛けるだけでは対数モーメントは得られない。

:::example[三乗分母の一次対数]
原点側の重みが比較的大きい。一次対数の符号は負、二次対数の符号は正になる。 ここでは $q=3,\ a=1,\ s=1$ を取る。対数を付ける前の関数は

$$
M(s)=\int_0^\infty\frac{x^{s-1}}{x^{3} + 1}\,\dd x.
$$

積分前に両端を調べると、原点では $x^{0}$ の定数倍、無限遠では $x^{-3}$ に等価である。従って対数の二乗まで付けても収束する。正弦と余接を指定点で評価すれば

$$
\begin{aligned}
\sin(\pi s/q)&=\frac{\sqrt{3}}{2}, \\
M^\prime/M&=- \frac{\sqrt{3} \pi}{9}.
\end{aligned}
$$

三つのモーメントを別々に書くと

$$
\begin{aligned}
\int_0^\infty \frac{1}{x^{3} + 1}\,\dd x&=\frac{2 \sqrt{3} \pi}{9}, \\
\int_0^\infty \frac{1}{x^{3} + 1}\log x\,\dd x&=- \frac{2 \pi^{2}}{27}, \\
\int_0^\infty \frac{1}{x^{3} + 1}(\log x)^2\,\dd x&=\frac{10 \sqrt{3} \pi^{3}}{243}.
\end{aligned}
$$
:::

:::example[逆数変換の固定指数]
中心指数 $s=q/2$ では $x$ と $1/x$ の寄与が釣り合い、奇数次の対数モーメントが消える。 ここでは $q=4,\ a=1,\ s=2$ を取る。対数を付ける前の関数は

$$
M(s)=\int_0^\infty\frac{x^{s-1}}{x^{4} + 1}\,\dd x.
$$

積分前に両端を調べると、原点では $x^{1}$ の定数倍、無限遠では $x^{-3}$ に等価である。従って対数の二乗まで付けても収束する。正弦と余接を指定点で評価すれば

$$
\begin{aligned}
\sin(\pi s/q)&=1, \\
M^\prime/M&=0.
\end{aligned}
$$

三つのモーメントを別々に書くと

$$
\begin{aligned}
\int_0^\infty \frac{x}{x^{4} + 1}\,\dd x&=\frac{\pi}{4}, \\
\int_0^\infty \frac{x}{x^{4} + 1}\log x\,\dd x&=0, \\
\int_0^\infty \frac{x}{x^{4} + 1}(\log x)^2\,\dd x&=\frac{\pi^{3}}{64}.
\end{aligned}
$$
:::

:::example[尺度による対称点の移動]
無次元の積分変数を使えば中央で消える一次対数が、もとの $x$ では正の値になる。 ここでは $q=2,\ a=4,\ s=1$ を取る。対数を付ける前の関数は

$$
M(s)=\int_0^\infty\frac{x^{s-1}}{x^{2} + 4}\,\dd x.
$$

積分前に両端を調べると、原点では $x^{0}$ の定数倍、無限遠では $x^{-2}$ に等価である。従って対数の二乗まで付けても収束する。正弦と余接を指定点で評価すれば

$$
\begin{aligned}
\sin(\pi s/q)&=1, \\
M^\prime/M&=\log{\left(2 \right)}.
\end{aligned}
$$

三つのモーメントを別々に書くと

$$
\begin{aligned}
\int_0^\infty \frac{1}{x^{2} + 4}\,\dd x&=\frac{\pi}{4}, \\
\int_0^\infty \frac{1}{x^{2} + 4}\log x\,\dd x&=\frac{\pi \log{\left(2 \right)}}{4}, \\
\int_0^\infty \frac{1}{x^{2} + 4}(\log x)^2\,\dd x&=\frac{\pi \left(4 \log{\left(2 \right)}^{2} + \pi^{2}\right)}{16}.
\end{aligned}
$$
:::

:::example[四分の三の指数が現れる場合]
$s$ を $q-s$ に取り替えると尺度一の基本積分は同じ値をもつが、その一次微分の符号は反転する。 ここでは $q=4,\ a=1,\ s=3$ を取る。対数を付ける前の関数は

$$
M(s)=\int_0^\infty\frac{x^{s-1}}{x^{4} + 1}\,\dd x.
$$

積分前に両端を調べると、原点では $x^{2}$ の定数倍、無限遠では $x^{-2}$ に等価である。従って対数の二乗まで付けても収束する。正弦と余接を指定点で評価すれば

$$
\begin{aligned}
\sin(\pi s/q)&=\frac{\sqrt{2}}{2}, \\
M^\prime/M&=\frac{\pi}{4}.
\end{aligned}
$$

三つのモーメントを別々に書くと

$$
\begin{aligned}
\int_0^\infty \frac{x^{2}}{x^{4} + 1}\,\dd x&=\frac{\sqrt{2} \pi}{4}, \\
\int_0^\infty \frac{x^{2}}{x^{4} + 1}\log x\,\dd x&=\frac{\sqrt{2} \pi^{2}}{16}, \\
\int_0^\infty \frac{x^{2}}{x^{4} + 1}(\log x)^2\,\dd x&=\frac{3 \sqrt{2} \pi^{3}}{64}.
\end{aligned}
$$
:::

:::example[二つの側で異なる可積分性]
原点では冪が零となり、無限遠では速く減衰する。$s$ を零まで外挿する前に原点での対数発散を調べる必要がある。 ここでは $q=6,\ a=1,\ s=1$ を取る。対数を付ける前の関数は

$$
M(s)=\int_0^\infty\frac{x^{s-1}}{x^{6} + 1}\,\dd x.
$$

積分前に両端を調べると、原点では $x^{0}$ の定数倍、無限遠では $x^{-6}$ に等価である。従って対数の二乗まで付けても収束する。正弦と余接を指定点で評価すれば

$$
\begin{aligned}
\sin(\pi s/q)&=\frac{1}{2}, \\
M^\prime/M&=- \frac{\sqrt{3} \pi}{6}.
\end{aligned}
$$

三つのモーメントを別々に書くと

$$
\begin{aligned}
\int_0^\infty \frac{1}{x^{6} + 1}\,\dd x&=\frac{\pi}{3}, \\
\int_0^\infty \frac{1}{x^{6} + 1}\log x\,\dd x&=- \frac{\sqrt{3} \pi^{2}}{18}, \\
\int_0^\infty \frac{1}{x^{6} + 1}(\log x)^2\,\dd x&=\frac{7 \pi^{3}}{108}.
\end{aligned}
$$
:::

:::example[二次分母と平方根の分子]
$x$ の冪は整数でなく、実軸上で実対数を使うという約束が鍵穴の両側の因子を一意にする。 ここでは $q=2,\ a=1,\ s=\frac{3}{2}$ を取る。対数を付ける前の関数は

$$
M(s)=\int_0^\infty\frac{x^{s-1}}{x^{2} + 1}\,\dd x.
$$

積分前に両端を調べると、原点では $x^{\frac{1}{2}}$ の定数倍、無限遠では $x^{- \frac{3}{2}}$ に等価である。従って対数の二乗まで付けても収束する。正弦と余接を指定点で評価すれば

$$
\begin{aligned}
\sin(\pi s/q)&=\frac{\sqrt{2}}{2}, \\
M^\prime/M&=\frac{\pi}{2}.
\end{aligned}
$$

三つのモーメントを別々に書くと

$$
\begin{aligned}
\int_0^\infty \frac{\sqrt{x}}{x^{2} + 1}\,\dd x&=\frac{\sqrt{2} \pi}{2}, \\
\int_0^\infty \frac{\sqrt{x}}{x^{2} + 1}\log x\,\dd x&=\frac{\sqrt{2} \pi^{2}}{4}, \\
\int_0^\infty \frac{\sqrt{x}}{x^{2} + 1}(\log x)^2\,\dd x&=\frac{3 \sqrt{2} \pi^{3}}{8}.
\end{aligned}
$$
:::

### 非整数べきがある分母を二方向から評価する

:::example[$1+2x\cos\theta+x^2$ の Mellin 積分]
$0<\theta<\pi$、$0<p<2$ とする。$x>0$ では二次分母は $(x+\cos\theta)^2+\sin^2\theta>0$ である。原点と無限遠における冪指数は $p-1$ と $p-3$ だから
$$
L_\theta(p)=\int_0^\infty\frac{x^{p-1}}{1+2x\cos\theta+x^2}\dd x
$$
はこの範囲で絶対収束する。前と同じ正実軸の切断を取り、二極 $-e^{i\theta}$、$-e^{-i\theta}$ の偏角をそれぞれ $\pi+\theta,\pi-\theta$ とする。$0<\theta<\pi$ のため両方とも指定範囲 $(0,2\pi)$ に入る。この二極での留数和は
$$
\frac{e^{i(p-1)(\pi-\theta)}-e^{i(p-1)(\pi+\theta)}}{2i\sin\theta}
=-e^{i\pi(p-1)}\frac{\sin((p-1)\theta)}{\sin\theta}.
$$
外弧の評価は $O(R^{p-2})$、内弧は $O(\eps^p)$ である。$p\ne1$ なら上下辺の係数 $1-e^{2\pi ip}$ は非零なので
$$
L_\theta(p)=\frac{\pi\sin((1-p)\theta)}{\sin(\pi p)\sin\theta}.
$$
$x=1/u$ の変換により $L_\theta(p)=L_\theta(2-p)$ であり、右辺も同じ対称性を満たす。$p=1$ は相殺による $0/0$ であって特異な実積分ではない。$p\to1$ の極限を取ると $L_\theta(1)=\theta/\sin\theta$。平方完成して実積分を逆正接で求めても同じである。$\theta=\pi/2$ では $\int_0^\infty x^{p-1}/(1+x^2)\dd x=\pi/[2\sin(\pi p/2)]$ が戻る。
:::

## Cauchy の主値積分

:::definition[Cauchy の主値]{#def-cpx-principal-value}
実軸上 $a$ に特異点があるとき、両側から同じ幅 $\eps$ を除き、左右の積分の和を取ってから $\eps\downarrow0$ とした極限を :term[Cauchy の主値]{en="Cauchy principal value" yomi="こーしーのしゅち"} という。無限区間の両端を相殺させる場合は $[-R,R]$ という切り方も指定する。
:::

:::proposition[実軸上の単純極]{#prop-cpx-indentation}
上半平面で閉じ、実軸を左から右へ進み、実軸上の単純極を小上半円で避ける。大円弧の寄与が消えるとき、実軸の主値は上半平面内部の留数和の $2\pi i$ 倍に実軸上の留数和の $\pi i$ 倍を足した値である。ただし関数は該当する閉半平面の近傍で有理型で、極は有限個とする。
:::

:::proof
局所表示は $c/(z-a)+h(z)$、$h$ は有界である。小円弧 $z=a+\eps e^{it}$ では $t$ は $\pi$ から $0$ に減少するので、その積分は $-\pi i c+O(\eps)$。極は閉路の外にあるため内部の留数には加えない。小円弧の積分を右辺に移すと符号が正になる。まず有限半径で各小弧の半径を零に送ると、直線部分は局所主値に収束する。次に外半径を無限大へ送り、大円弧の項を消す。無限遠で得るのは対称な切り方による極限である。
:::

:::example[除き方が異なる場合]
$\int_{-1}^1\dd x/x$ の主値は零である。負側を $-\eps$ まで、正側を $2\eps$ からとすると和は $-\log2$ に変わる。従って通常の広義積分では、独立な片側極限の存在を要求する。
:::

### 広義積分、主値、有限部分の違い

:::example[二位の極の両側は相殺しない]
$\int_{-1}^{-\eps}x^{-2}\dd x+\int_\eps^1 x^{-2}\dd x=2/\eps-2$ である。正負の側を同じ幅で除いても極限は存在しない。発散項 $2/\eps$ を引いた後の定数 $-2$ を Hadamard の有限部分と呼ぶことがあるが、これは引く項を指定した別の操作である。留数が零の二位の極にもこの発散は生じる。
:::

:::example[半円を下側に回したときの帳尻]
$e^{ix}/x$ の上半平面の閉路で、原点だけを実軸の下にふくらませて通る。この小弧は $z=\eps e^{it}$、$-\pi\le t\le0$ なので積分は $+\pi i+o(1)$。今回は原点が閉路の内側にあり、留数定理の右辺も $2\pi i$ 増える。従って直線部分の主値は $2\pi i-\pi i=\pi i$ と、上にくぼませた場合と等しい。小弧の符号を変えるときは内部の極の集合も同時に変更する必要がある。
:::

:::example[Sokhotski–Plemelj の直線上の形]
実直線上 $C^1$ 級でコンパクト台をもつ $\phi$ と実数 $a$ に対し、$\eps\downarrow0$ のとき
$$
\int_\R\frac{\phi(x)}{x-a\mp i\eps}\dd x\longrightarrow
\operatorname{PV}\int_\R\frac{\phi(x)}{x-a}\dd x\ \pm i\pi\phi(a).
$$
$u=x-a$ とおくと分母の逆数は $u/(u^2+\eps^2)\pm i\eps/(u^2+\eps^2)$。前半は $|u|\le1$ で $\phi(a)$ を引いても奇性で変わらず、残りは平均値定理で一様に有界となり優収束を使える。後半は $u=\eps v$ と置いて $\phi(a+\eps v)/(1+v^2)$ の積分であり、上界 $\|\phi\|_\infty/(1+v^2)$ と $\int_\R\dd v/(1+v^2)=\pi$ を使う。極をどちらから実軸に近づけるかで、半留数の符号が決まることを実積分だけでも証明できた。
:::

### 高位の極で対称主値が存在する条件

孤立した極 $a\in\R$ の局所変数を $u=x-a$ とし、主部を $\sum_{k=1}^m c_{-k}u^{-k}$ とする。固定した小区間 $[-\delta,\delta]$ から $(-\eps,\eps)$ を除く。奇数 $k$ の項は奇関数なので左右の積分が相殺する。偶数 $k\ge2$ の項からは
$$
\frac{2c_{-k}}{k-1}\left(\eps^{1-k}-\delta^{1-k}\right)
$$
が残る。最も大きい偶数 $k$ で $c_{-k}\ne0$ なら、その発散はそれより低い次数の項で打ち消せない。正則部分は普通に積分できる。従って局所対称主値が存在するための必要十分条件は、負の偶数次数の係数がすべて零になることである。積分区間の左右の端点が異なる場合も、極を中心とする共通の小区間を切り出せば同じ判定になる。

例えば $\int_{-1}^2 x^{-3}\dd x$ の通常の広義積分では両側が別々に発散する。対称主値では近くの三次の項が相殺し、残る $\int_1^2 x^{-3}\dd x=3/8$ が答えとなる。高位の極があるときは、小円弧を計算する前にこのように主値の存在を確認する。

### 実軸上に複数の極があるときの主値

実軸上の極を一つずつ除くとき、各極の左右の幅は同じにし、異なる極の間では同じ幅である必要はない。有限個なら各点の Laurent 展開から単純極を引いた連続部分を別々に評価できるためである。上側で閉じると内部の寄与は $2\pi i$ 倍、実軸上の寄与は $\pi i$ 倍となる。実係数有理関数を用いた以下の例では最終値は実数であり、虚数部分が消えるかどうかを符号の検算に使える。

:::example[原点の奇特異性と共役極]
経路を上に閉じたときの複素留数の寄与は個々には虚数である。しかし元の実関数は奇関数で、左右対称な主値は零になるはずである。

$$
I=\operatorname{PV}\int_\R \frac{1}{x \left(x^{2} + 1\right)}\,\dd x.
$$

分母を因数分解し、実軸上と上半平面の極に分ける。

$$
\begin{aligned}
a=0:\quad \operatorname{Res}_a F&=1\quad(a\in\R), \\
a=i:\quad \operatorname{Res}_a F&=- \frac{1}{2}\quad(\Im a>0),
\end{aligned}
$$

大円弧は分母と分子の次数差が二以上なので消える。実極を上に避ける場合、直線部分以外に残る小円弧の和は $- i \pi$。従って

$$
I=2\pi i\left(- \frac{1}{2}\right)+\pi i\left(1\right)=0.
$$

実変数からも、特異な部分を分離した $F(z)-(1)/(z-(0))=- \frac{z}{z^{2} + 1}$ が $z=0$ で有限であることを確かめられる。
:::

:::example[中心をずらした重み]
実極の位置を右へ移したため全体の奇性が失われる。原点の回りだけで機械的に左右を切る方法では特異点を除けない。

$$
I=\operatorname{PV}\int_\R \frac{1}{\left(x - 1\right) \left(x^{2} + 1\right)}\,\dd x.
$$

分母を因数分解し、実軸上と上半平面の極に分ける。

$$
\begin{aligned}
a=1:\quad \operatorname{Res}_a F&=\frac{1}{2}\quad(a\in\R), \\
a=i:\quad \operatorname{Res}_a F&=- \frac{1}{4} + \frac{i}{4}\quad(\Im a>0),
\end{aligned}
$$

大円弧は分母と分子の次数差が二以上なので消える。実極を上に避ける場合、直線部分以外に残る小円弧の和は $- \frac{i \pi}{2}$。従って

$$
I=2\pi i\left(- \frac{1}{4} + \frac{i}{4}\right)+\pi i\left(\frac{1}{2}\right)=- \frac{\pi}{2}.
$$

実変数からも、特異な部分を分離した $F(z)-(\frac{1}{2})/(z-(1))=- \frac{z + 1}{2 \left(z^{2} + 1\right)}$ が $z=1$ で有限であることを確かめられる。
:::

:::example[二つの実極と一対の虚極]
実極の左右にはそれぞれ対数発散がある。各極の局所留数を計算し、同じ向きの小弧を一つずつ入れる。

$$
I=\operatorname{PV}\int_\R \frac{1}{\left(x - 1\right) \left(x + 2\right) \left(x^{2} + 1\right)}\,\dd x.
$$

分母を因数分解し、実軸上と上半平面の極に分ける。

$$
\begin{aligned}
a=1:\quad \operatorname{Res}_a F&=\frac{1}{6}\quad(a\in\R), \\
a=-2:\quad \operatorname{Res}_a F&=- \frac{1}{15}\quad(a\in\R), \\
a=i:\quad \operatorname{Res}_a F&=- \frac{1}{20} + \frac{3 i}{20}\quad(\Im a>0),
\end{aligned}
$$

大円弧は分母と分子の次数差が二以上なので消える。実極を上に避ける場合、直線部分以外に残る小円弧の和は $- \frac{i \pi}{10}$。従って

$$
I=2\pi i\left(- \frac{1}{20} + \frac{3 i}{20}\right)+\pi i\left(\frac{1}{10}\right)=- \frac{3 \pi}{10}.
$$
:::

:::example[実極で正則成分が消える例]
分子は実極で非零である。実極を閉路から除いた場合と内部に入れた場合の双方で同じ主値が得られる。

$$
I=\operatorname{PV}\int_\R \frac{x + 2}{\left(x - 2\right) \left(x^{2} + 4\right)}\,\dd x.
$$

分母を因数分解し、実軸上と上半平面の極に分ける。

$$
\begin{aligned}
a=2:\quad \operatorname{Res}_a F&=\frac{1}{2}\quad(a\in\R), \\
a=2 i:\quad \operatorname{Res}_a F&=- \frac{1}{4}\quad(\Im a>0),
\end{aligned}
$$

大円弧は分母と分子の次数差が二以上なので消える。実極を上に避ける場合、直線部分以外に残る小円弧の和は $- \frac{i \pi}{2}$。従って

$$
I=2\pi i\left(- \frac{1}{4}\right)+\pi i\left(\frac{1}{2}\right)=0.
$$

実変数からも、特異な部分を分離した $F(z)-(\frac{1}{2})/(z-(2))=- \frac{z}{2 \left(z^{2} + 4\right)}$ が $z=2$ で有限であることを確かめられる。
:::

:::example[複素極が二位の場合]
実軸上で半分になるのは単純極の寄与だけで、上半平面の二位の極には通常の留数定理を用いる。

$$
I=\operatorname{PV}\int_\R \frac{1}{\left(x - 1\right) \left(x^{2} + 1\right)^{2}}\,\dd x.
$$

分母を因数分解し、実軸上と上半平面の極に分ける。

$$
\begin{aligned}
a=1:\quad \operatorname{Res}_a F&=\frac{1}{4}\quad(a\in\R), \\
a=i:\quad \operatorname{Res}_a F&=- \frac{1}{8} + \frac{i}{4}\quad(\Im a>0),
\end{aligned}
$$

大円弧は分母と分子の次数差が二以上なので消える。実極を上に避ける場合、直線部分以外に残る小円弧の和は $- \frac{i \pi}{4}$。従って

$$
I=2\pi i\left(- \frac{1}{8} + \frac{i}{4}\right)+\pi i\left(\frac{1}{4}\right)=- \frac{\pi}{2}.
$$

実変数からも、特異な部分を分離した $F(z)-(\frac{1}{4})/(z-(1))=- \frac{\left(z + 1\right) \left(z^{2} + 3\right)}{4 \left(z^{2} + 1\right)^{2}}$ が $z=1$ で有限であることを確かめられる。
:::

## 留数による級数の和の計算

:::proposition[整数に極をもつ核]{#prop-cpx-cot-sum}
整数上に極をもたない有理関数 $g=O(|z|^{-2})$ に対し、
$$
\sum_{n\in\Z}g(n)=-\sum_{a:\,g\text{ の極}}\operatorname{Res}_{a}\bigl(\pi\cot(\pi z)g(z)\bigr).
$$
:::

:::proof
$\pi\cot\pi z$ は各整数で留数 $1$ をもつ。辺が $\Re z,\Im z=\pm(N+1/2)$ の正方形で積分する。垂直辺では $\cot\pi z$ の絶対値は $|\tanh(\pi\Im z)|\le1$、水平辺では指数表示から $(1+e^{-2\pi|\Im z|})/(1-e^{-2\pi|\Im z|})$ 以下である。従って境界全体の積分は $O(N^{-1})$。留数定理で境界内の全留数和を零に近づければ結論となる。
:::

:::example[双曲線余接の和]
$g(z)=1/(z^2+a^2)$、$a>0$ では $\cot(iu)=-i\coth u$。$ia$ での積の留数は $-\pi\coth(\pi a)/(2a)$、$-ia$ でも同じ値である。従って $\sum_{n\in\Z}1/(n^2+a^2)=\pi\coth(\pi a)/a$。$n=0$ の項を除き $a\to0$ とすれば $2\sum_{n\ge1}1/n^2=\pi^2/3$ を得る。極限交換は $n^{-2}$ による支配で正当化される。
:::

### 交代符号をもつ和

:::example[$\sum_{n\in\Z}(-1)^n/(n^2+a^2)$]
整数 $n$ で留数が $(-1)^n$ になる核は $\pi/\sin\pi z$ である。前と同じ正方形の垂直辺では $|\sin\pi z|=\cosh(\pi\Im z)\ge1$、水平辺では $|\sin\pi z|\ge\sinh(\pi(N+1/2))$。従って核は境界全体で一様有界であり、$1/(z^2+a^2)$ を掛ければ境界積分は消える。$ia$ での積核の留数は $\pi/[(2ia)i\sinh(\pi a)]=-\pi/[2a\sinh(\pi a)]$、$-ia$ でも同じである。従って求める全整数の和は $\pi/[a\sinh(\pi a)]$ となる。交代しない場合の $\coth(\pi a)$ が $1/\sinh(\pi a)$ に替わる。$a\downarrow0$ では $n=0$ の項 $a^{-2}$ を先に除けば $\sum_{n\ge1}(-1)^{n-1}/n^2=\pi^2/12$ である。
:::

:::example[半整数を取る核]
半整数 $n+1/2$ で留数が $1$ になる核は $-\pi\tan\pi z$ である。$g(z)=1/(z^2+a^2)$ とし、今度は半整数に当たらないよう、辺を $\Re z,\Im z=\pm N$ とする。十分大きい $N>a$ では $g$ の極も境界を避ける。垂直辺では $|\tan(\pi z)|=|\tanh(\pi\Im z)|\le1$、水平辺では指数表示により $|\tan(\pi z)|\le(1+e^{-2\pi N})/(1-e^{-2\pi N})$ である。従って境界積分は $O(N^{-1})$ で消える。$ia$ での留数は $-\pi\tanh(\pi a)/(2a)$、$-ia$ でも同じである。従って $\sum_{n\in\Z}[(n+1/2)^2+a^2]^{-1}=\pi\tanh(\pi a)/a$。$a\to0$ で左辺は $\sum(n+1/2)^{-2}$、右辺は $\pi^2$。正の奇数だけに書き直すと $\sum_{m\ge0}(2m+1)^{-2}=\pi^2/8$ となる。格子を半分ずらすと、正負の添字の対応は $n\leftrightarrow -n-1$ になる点に注意する。
:::

### 整数格子の和と重極の相互作用

$\sum_{n\in\mathbb Z}(n^2+a^2)^{-1}$ を出発点にして、指数の大きい分母の和を作る。$a$ が正実軸のコンパクト区間にあれば、任意の固定回数の $a$ 微分後にも $|n|\ge1$ の各項は定数倍の $n^{-2}$ で抑えられる。そこで
$$
S_j(a)=\sum_{n\in\mathbb Z}(n^2+a^2)^{-j},\qquad
S_{j+1}(a)=-\frac{1}{2ja}\frac{\dd}{\dd a}S_j(a)
$$
が使える。同じ操作は符号 $(-1)^n$ を付けた場合や半整数格子にも使えるが、最初に選ぶ積分核が変わる。

:::example[整数上の二乗分母]
正の項だけの級数であるため、和が正になるかを最終式の符号からも確かめられる。 一乗分母の和と微分から

$$
\begin{aligned}
S_1(a)&=\frac{\pi \coth{\left(\pi a \right)}}{a}, \\
S_2(a)&=-\frac{1}{2a} S_1^\prime(a).
\end{aligned}
$$

右辺を積の微分で評価すれば

$$
\sum_{n\in\mathbb Z}\frac{1}{(n^2+a^2)^2}=\frac{\pi^{2}}{2 a^{2} \sinh^{2}{\left(\pi a \right)}} + \frac{\pi \coth{\left(\pi a \right)}}{2 a^{3}}.
$$

特に $a=1$ の全格子和は $\frac{\pi \left(\frac{\pi}{\sinh^{2}{\left(\pi \right)}} + \frac{1}{\tanh{\left(\pi \right)}}\right)}{2}$ である。正整数部分を取り出す際は、原点の項 $1$ を差し引いてから二で割る。
:::

:::example[交代する二乗分母]
交代符号により大きい $a$ での和が指数的に小さくなる。 一乗分母の和と微分から

$$
\begin{aligned}
S_1(a)&=\frac{\pi}{a \sinh{\left(\pi a \right)}}, \\
S_2(a)&=-\frac{1}{2a} S_1^\prime(a).
\end{aligned}
$$

右辺を積の微分で評価すれば

$$
\sum_{n\in\mathbb Z}\frac{(-1)^n}{(n^2+a^2)^2}=\frac{\pi^{2} \cosh{\left(\pi a \right)}}{2 a^{2} \sinh^{2}{\left(\pi a \right)}} + \frac{\pi}{2 a^{3} \sinh{\left(\pi a \right)}}.
$$

特に $a=1$ の全格子和は $\frac{\pi \left(1 + \frac{\pi}{\tanh{\left(\pi \right)}}\right)}{2 \sinh{\left(\pi \right)}}$ である。正整数部分を取り出す際は、原点の項 $1$ を差し引いてから二で割る。
:::

:::example[半整数上の二乗分母]
原点に格子点がないので $a$ が零に近づいても原点の項を差し引く必要がない。 一乗分母の和と微分から

$$
\begin{aligned}
S_1(a)&=\frac{\pi \tanh{\left(\pi a \right)}}{a}, \\
S_2(a)&=-\frac{1}{2a} S_1^\prime(a).
\end{aligned}
$$

右辺を積の微分で評価すれば

$$
\sum_{n\in\mathbb Z}\frac{1}{((n+1/2)^2+a^2)^2}=\frac{\pi^{2} \tanh^{2}{\left(\pi a \right)}}{2 a^{2}} - \frac{\pi^{2}}{2 a^{2}} + \frac{\pi \tanh{\left(\pi a \right)}}{2 a^{3}}.
$$

特に $a=1$ の全格子和は $- \frac{\pi \left(- \tanh{\left(\pi \right)} + \frac{\pi}{\cosh^{2}{\left(\pi \right)}}\right)}{2}$ である。非負整数部分の和はこの半分である。負の添字側の $-n-1$ が正側の同じ分母に対応するためである。
:::

## 例：Fresnel 積分と Dirichlet 積分

:::example[Dirichlet 積分]
$b>0$ で $e^{ibz}/z$ を上半円で閉じ、原点を上側に避ける。Jordan の補題と半留数の公式より $\operatorname{PV}\int_\R e^{ibx}\dd x/x=\pi i$。虚部 $\sin(bx)/x$ は原点で連続に延長でき、無限遠では部分積分により片側の尾部が $O(1/R)$ となる。偶性から $\int_0^\infty\sin(bx)\dd x/x=\pi/2$。$b=0$ では零、$b<0$ なら $-\pi/2$ である。
:::

:::example[Fresnel 積分]
$e^{-z^2}$ を角 $0\le\arg z\le\pi/4$ の扇形で積分する。円弧上で絶対値は $e^{-R^2\cos 2\theta}$。$u=\pi/4-\theta$ と置けば $\cos2\theta=\sin2u\ge4u/\pi$、従って円弧積分の絶対値は $R\int_0^{\pi/4}e^{-4R^2u/\pi}\dd u\le\pi/(4R)$ である。二本の半径の積分を比較し
$$
\frac{\sqrt\pi}{2}=\int_0^\infty e^{-x^2}\dd x=e^{i\pi/4}\int_0^\infty e^{-it^2}\dd t
$$
を得る。左の実 Gauss 積分は二乗して第一象限の二重積分に移し、極座標で $\int_0^{\pi/2}\dd\theta\int_0^\infty re^{-r^2}\dd r=\pi/4$ と求める。従って余弦と正弦の各積分は $\sqrt{\pi/8}$。円弧の評価は振動積分の極限の存在も同時に与える。
:::

### 減衰を付けてから外す別証

:::example[Dirichlet 積分の Abel 正則化]
$\eta>0$、$b\ge0$ とし $D_\eta(b)=\int_0^\infty e^{-\eta x}\sin(bx)\dd x/x$ とする。$b$ で微分して $D_\eta'(b)=\Re[(\eta-ib)^{-1}]=\eta/(\eta^2+b^2)$、初期値 $D_\eta(0)=0$ より $D_\eta(b)=\arctan(b/\eta)$。固定 $b>0$ で $\eta\downarrow0$ とすれば $\pi/2$ である。これが元の広義積分になるためには、広がる無限区間で極限を勝手に交換せず尾部を制御する。$e^{-\eta x}/x$ は単調に減少し、$\sin(bx)$ の任意の区間での積分の絶対値は $2/b$ 以下。部分積分により $R$ より先の尾部は $2/(bR)$ 以下で、$\eta$ に依存しない。有限区間では優収束が使えるから主張が従う。
:::

:::example[振動の次数を変えた積分]
$p>1$、$0<s<p$ とし、$\int_0^\infty x^{s-1}e^{ix^p}\dd x$ を求める。$u=x^p$ と置くと $(1/p)\int_0^\infty u^{s/p-1}e^{iu}\dd u$。原点では $s/p>0$、無限遠では $s/p-1<0$ なので Dirichlet 判定により通常の広義積分として存在する。$\eta>0$ を入れると Laplace 積分は $\Gamma(s/p)(\eta-i)^{-s/p}/p$ である。ここでは正の実数 $\lambda$ に対して $\Gamma(\lambda)=\int_0^\infty v^{\lambda-1}e^{-v}\dd v$ と定義する。Laplace 変数が正実数のときは実変数の置換で式が従う。右半平面では小コンパクト集合ごとに指数減衰による可積分上界があるので積分は正則で、右辺も主対数で正則である。よって一致の定理で右半平面へ延びる。前例と同様に単調に減少する振幅の尾部を $\eta$ に無関係に抑えられるので、極限は $\Gamma(s/p)e^{i\pi s/(2p)}/p$。特に $p=2,s=1$ で Fresnel 積分が再現される。[ガンマ関数](/math/complex-analysis/gamma-zeta/)の Euler 積分との接点である。
:::

### 扇形積分の端の評価を一般化する

:::example[回転する角を小さくして極限を取る]
$0\le\alpha<\pi/4$ の直線に沿って $\int_0^\infty e^{-z^2}\dd z$ を計算する。$z=te^{i\alpha}$ なら $|e^{-z^2}|=\exp(-t^2\cos2\alpha)$ であり、$\cos2\alpha>0$ だから直線上の積分は絶対収束する。角 $\alpha$ の扇形では円弧上の積分の絶対値が $\alpha R e^{-R^2\cos2\alpha}$ 以下であり、従って
$$
\int_0^\infty e^{-t^2e^{2i\alpha}}\dd t
=e^{-i\alpha}\frac{\sqrt\pi}{2}
$$
である。$\alpha\uparrow\pi/4$ ではこの指数減衰の定数が零になり、無限区間の優収束はこのままでは使えない。そこで積分尾部を部分積分する。$c=e^{2i\alpha}$ として $(e^{-ct^2})'=-2ct e^{-ct^2}$ だから、$T>R>0$ では
$$
\int_R^T e^{-ct^2}\dd t
=\left[-\frac{e^{-ct^2}}{2ct}\right]_R^T-\frac1{2c}\int_R^T\frac{e^{-ct^2}}{t^2}\dd t .
$$
$|c|=1,\Re c\ge0$ であるから境界項と積分を合わせて $1/R+1/(2T)$ 以下に抑えられる。これなら $\alpha$ に依存しない。有限区間で極限を取り、次いで $R\to\infty$ とすることで、境界の振動積分の値を正当化できる。直接の円弧評価と異なり、この方法は角度を動かすときにも同じ尾部評価を使える。
:::

:::example[半直線で片側だけを取り出すとき]
Fresnel の実部と虚部はともに正で等しいが、$e^{it^2}$ の実部そのものは正とは限らない。積分の符号を示したいだけなら、$u=t^2$ として $\int_0^\infty\sin(t^2)\dd t=\frac12\int_0^\infty u^{-1/2}\sin u\,\dd u$ と置く方法もある。長さ $\pi$ の各区間で符号は交代し、絶対値の積分は正の単調減少列となる。実際、$u$ を $\pi$ 平行移動して二つの区間を同じ区間に載せると、$|\sin u|$ は変わらず重みだけが減る。従って正の最初の項が次の負の項を上回り、以後も二項ずつまとめて正となる。一方、絶対値の積分は発散する。各区間の中央で $|\sin u|\ge1/2$ となる一定長の部分を取れば、下界が $\sum n^{-1/2}$ の定数倍となるからである。この条件収束を確認してから複素積分で値を求めると、積分順序や絶対値を交換できる場所が明確になる。
:::

## 演習問題

::::exercise{level=basic}
$\int_0^{2\pi}\dd t/(3+2\cos t)$ を求めよ。

:::solution
$3>2$ だから単位円上に極はない。余弦の公式で $2\pi/\sqrt{9-4}=2\pi/\sqrt5$。
:::
::::

::::exercise{level=basic}
$\int_{\R}\dd x/[(x^2+1)(x^2+4)]$ を求めよ。

:::solution
被積分関数は $[(x^2+1)^{-1}-(x^2+4)^{-1}]/3$ に等しい。各項の上側の単純極を取ると $(\pi-\pi/2)/3=\pi/6$。
:::
::::

::::exercise{level=basic}
$\int_0^\infty x^2/(1+x^4)\dd x$ を求めよ。

:::solution
$x=1/t$ と置けば微分の $t^{-2}$ と分子分母のべきが消え、$\int_0^\infty\dd t/(1+t^4)=\pi/(2\sqrt2)$ である。
:::
::::

::::exercise{level=basic}
$\int_{\R}\cos(3x)/(x^2+4)\dd x$ を求めよ。

:::solution
上半円で $e^{3iz}$ を使い、極 $2i$ での留数は $e^{-6}/(4i)$。実部を取れば $\pi e^{-6}/2$。
:::
::::

::::exercise{level=advanced}
$\int_0^\infty x^{-2/3}/(1+x)\dd x$ を求めよ。

:::solution
$p-1=-2/3$ より $p=1/3$ で収束域内にある。鍵穴型公式から $\pi/\sin(\pi/3)=2\pi/\sqrt3$。
:::
::::

::::exercise{level=advanced}
$\operatorname{PV}\int_{-2}^3 \dd x/x$ を求めよ。

:::solution
負側が $\log\eps-\log2$、正側が $\log3-\log\eps$。和の極限は $\log(3/2)$ である。
:::
::::

::::exercise{level=advanced}
$\sum_{n\ge1}(n^2+1)^{-1}$ を求めよ。

:::solution
全整数の和 $\pi\coth\pi$ から原点の $1$ を引き、正負の対称性で割ると $(\pi\coth\pi-1)/2$。
:::
::::

::::exercise{level=advanced}
$\int_0^\infty(1-\cos bx)/x^2\,\dd x$ を実数 $b$ に対して求めよ。

:::solution
原点では $1-\cos bx=O(x^2)$、遠方では有界なので絶対収束する。部分積分の境界項 $-(1-\cos bx)/x$ は両端で消え、残りは $b\int_0^\infty\sin(bx)/x\,\dd x=\pi|b|/2$ である。$b=0$ でも式は成立する。
:::
::::

## 確認テスト

::::quiz{answer="0.5" tol="1e-8"}
$\int_{\R}\dd x/(x^2+4)$ を $\pi$ で割った値はいくらか。

:::explanation
上側の極 $2i$ の留数は $1/(4i)$ なので積分は $\pi/2$。これを $\pi$ で割ると $1/2$ である。
:::
::::

::::quiz
$b<0$ で $e^{ibz}$ を含む Fourier 積分を閉じる側を選べ。

- [ ] 上半平面
- [x] 下半平面
- [ ] 常に両側

:::explanation
絶対値 $e^{-b\Im z}$ が減衰するには $b\Im z>0$ が必要である。
:::
::::

::::quiz{answer="-1" tol="1e-8"}
$1/z$ を原点の上で左から右へ避ける円弧積分を $\pi i$ で割るといくらか。

:::explanation
偏角が $\pi$ から $0$ に減るので、積分は $i(0-\pi)=-\pi i$ である。
:::
::::

::::quiz
級数を足すために $\pi\cot\pi z$ を使う理由を選べ。

- [x] 各整数で留数が $1$ になる。
- [ ] 整数で零になる。
- [ ] 全平面で有界である。

:::explanation
整数 $n$ で $\sin\pi z$ の導関数は $\pi\cos\pi n$ であり、分子の値と等しいので留数は $1$。
:::
::::

---
title: "複素積分の応用"
description: "留数・切断・境界値を使い、実積分、線形応答、二次元の境界値問題を計算する。"
notation:
  - symbol: '\operatorname{Res}_{z=a} f'
    meaning: "点 a における f の留数"
  - symbol: '\operatorname{PV}'
    meaning: "Cauchy の主値"
  - symbol: '\chi(\omega)'
    meaning: "正の指数を用いた因果応答の周波数表示"
---

## この章の目標

- 被積分関数に合わせて閉曲線と分岐を選び、積分を評価する。
- 実軸上の極を主値と境界値に分ける。
- 因果律から応答の実部と虚部の関係を導く。
- 等角写像で調和関数の境界値を移す。

## 留数による定積分（三角関数・有理関数・Fourier 型積分）

### 単位円への変数変換

角度が一周すると $z=e^{i\theta}$ が単位円を一周する。$\cos\theta=(z+z^{-1})/2$、$\sin\theta=(z-z^{-1})/(2i)$、$\dd\theta=\dd z/(iz)$ を代入すれば、有理式の周期積分を :ref[thm-cpx-residue-theorem] に帰着できる。元の被積分関数が有限でも、変数変換後の $z=0$ の極は調べる必要がある。

:::example[余弦を含む分母]
$a>|b|>0$ とする。$\int_0^{2\pi}(a+b\cos\theta)^{-1}\dd\theta$ では $\dd z$ の係数が $2/[i(bz^2+2az+b)]$ となる。二根は $z_\pm=(-a\pm\sqrt{a^2-b^2})/b$ で、積は $1$、$|z_+|=|b|/(a+\sqrt{a^2-b^2})<1$ である。分母の微分から $z_+$ での留数は $1/(i\sqrt{a^2-b^2})$ となる。従って値は $2\pi/\sqrt{a^2-b^2}$ であり、$a=5,b=3$ なら $\pi/2$ である。
:::

### 円弧の寄与

実軸に極を持たない有理関数 $R$ が無限遠で $O(|z|^{-2})$ なら、半径 $L$ の半円の寄与は長さ $\pi L$ と最大値の積で $O(L^{-1})$ となる。そこで実軸を左から右へ、上半円を反時計回りに閉じて上半平面の留数和を取る。

Fourier 型の積分では指数因子の符号を先に調べる。$e^{ikz}$ は $k>0$ のとき上半平面で減衰する。$R(z)=O(1/z)$ でも $\sin\theta\ge2\theta/\pi$（$0\le\theta\le\pi/2$）より
$$
\int_0^\pi e^{-kL\sin\theta}\dd\theta\le
2\int_0^{\pi/2}e^{-2kL\theta/\pi}\dd\theta\le\frac{\pi}{kL}.
$$
$\dd z$ の大きさ $L\dd\theta$ を掛けても寄与は $O(L^{-1})$ となる。これがこの場合の Jordan の補題である。$k<0$ のときは下で閉じ、向きは時計回りとなる。$k=0$ では別に円弧を評価する。

:::example[指数減衰するスペクトル]
$a>0$ とする。$k>0$ に対する $\int_\R e^{ikx}/(x^2+a^2)\dd x$ は極 $ia$ の留数 $e^{-ka}/(2ia)$ から $\pi e^{-ka}/a$ となる。$k<0$ では $-ia$ での留数 $-e^{ka}/(2ia)$ を $-2\pi i$ 倍するので $\pi e^{-|k|a}/a$ にまとまる。$a=2,k=1$ では $\pi e^{-2}/2$ である。
:::

:::example[四次の分母]
$\int_\R(1+x^4)^{-1}\dd x$ を上半円で閉じる。極は $e^{i\pi/4},e^{3i\pi/4}$ であり、各留数は $1/(4z^3)=-z/4$ となる。和は $-i\sqrt2/4$ なので実積分は $\pi/\sqrt2$ である。もし半径を $a$ で尺度調整して $\int(a^4+x^4)^{-1}$ とすれば、$x=au$ より全体は $a^{-3}$ 倍になる。
:::

### 極の位置と円弧の評価を別に行う

有理関数で分子次数が分母次数より二つ以上小さければ、半円の積分は一様な評価で消える。差が一つしかない場合、振動因子を含むかどうかが重要である。例えば $e^{ikz}$ を付けて $k>0$ とすれば、上半円でその絶対値は $e^{-kR\sin\theta}$ となる。$\sin\theta\ge2\theta/\pi$ を $0\le\theta\le\pi/2$ で使い、反対側は $\theta\mapsto\pi-\theta$ で移すと円弧積分を $O(R^{-1})$ に抑えられる。この評価では $k=0$ に一様な極限は取れない。実部、虚部を最後に取ると余弦、正弦の積分を分離できる。

:::example[極の前の多項式（1）]
実軸積分を複素化し、分母を因数分解して上半平面の極の寄与を選ぶ。 一次の模型から始める。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

実軸上の分母を定める。

$$
\begin{aligned}
f(x)\big|_{m=1} &= \frac{1}{x^{2} + 1}
\end{aligned}
$$

上半平面の極の留数は二つの根の差の逆数になり、$2\pi i$ を掛けると実数になる。

$$
\begin{aligned}
2\pi i\,\mathrm{Res}\,f\big|_{m=1} &= \pi
\end{aligned}
$$

指数因子は極の位置で評価する。この例では変換側の変数も $x$ と書く。全変換は $|x|$ を含むので、原点での微分には右極限を使う。

$$
\begin{aligned}
\widehat f(x)\quad (x>0)\big|_{m=1} &= \pi e^{- x}
\end{aligned}
$$

観測点 $x$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
f(x)\big|_{x=0} &= 1 \\
\frac{\partial^1 (f(x))}{\partial x^1} &= - \frac{2 x}{\left(x^{2} + 1\right)^{2}} \\
\left.\frac{\partial^1 (f(x))}{\partial x^1}\right|_{x=0} &= 0
\end{aligned}
$$

観測点 $x$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\widehat f(x)\big|_{x=0^+} &= \pi \\
\frac{\partial^1 (\widehat f(x)\quad (x>0))}{\partial x^1} &= - \pi e^{- x} \\
\left.\frac{\partial^1 (\widehat f(x)\quad (x>0))}{\partial x^1}\right|_{x=0^+} &= - \pi
\end{aligned}
$$
:::

:::example[極の前の多項式（2）]
実軸積分を複素化し、分母を因数分解して上半平面の極の寄与を選ぶ。 二次の模型を考える。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

実軸上の分母を定める。

$$
\begin{aligned}
f(x)\big|_{m=2} &= \frac{1}{x^{2} + 4}
\end{aligned}
$$

上半平面の極の留数は二つの根の差の逆数になり、$2\pi i$ を掛けると実数になる。

$$
\begin{aligned}
2\pi i\,\mathrm{Res}\,f\big|_{m=2} &= \frac{\pi}{2}
\end{aligned}
$$

指数因子は極の位置で評価する。この例では変換側の変数も $x$ と書く。全変換は $|x|$ を含むので、原点での微分には右極限を使う。

$$
\begin{aligned}
\widehat f(x)\quad (x>0)\big|_{m=2} &= \frac{\pi e^{- 2 x}}{2}
\end{aligned}
$$

観測点 $x$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
f(x)\big|_{x=0} &= \frac{1}{4} \\
\frac{\partial^1 (f(x))}{\partial x^1} &= - \frac{2 x}{\left(x^{2} + 4\right)^{2}} \\
\left.\frac{\partial^1 (f(x))}{\partial x^1}\right|_{x=0} &= 0
\end{aligned}
$$

観測点 $x$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\widehat f(x)\big|_{x=0^+} &= \frac{\pi}{2} \\
\frac{\partial^1 (\widehat f(x)\quad (x>0))}{\partial x^1} &= - \pi e^{- 2 x} \\
\left.\frac{\partial^1 (\widehat f(x)\quad (x>0))}{\partial x^1}\right|_{x=0^+} &= - \pi
\end{aligned}
$$
:::

:::example[極の前の多項式（3）]
実軸積分を複素化し、分母を因数分解して上半平面の極の寄与を選ぶ。 三次の模型を考える。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

実軸上の分母を定める。

$$
\begin{aligned}
f(x)\big|_{m=3} &= \frac{1}{x^{2} + 9}
\end{aligned}
$$

上半平面の極の留数は二つの根の差の逆数になり、$2\pi i$ を掛けると実数になる。

$$
\begin{aligned}
2\pi i\,\mathrm{Res}\,f\big|_{m=3} &= \frac{\pi}{3}
\end{aligned}
$$

指数因子は極の位置で評価する。この例では変換側の変数も $x$ と書く。全変換は $|x|$ を含むので、原点での微分には右極限を使う。

$$
\begin{aligned}
\widehat f(x)\quad (x>0)\big|_{m=3} &= \frac{\pi e^{- 3 x}}{3}
\end{aligned}
$$

観測点 $x$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
f(x)\big|_{x=0} &= \frac{1}{9} \\
\frac{\partial^1 (f(x))}{\partial x^1} &= - \frac{2 x}{\left(x^{2} + 9\right)^{2}} \\
\left.\frac{\partial^1 (f(x))}{\partial x^1}\right|_{x=0} &= 0
\end{aligned}
$$

観測点 $x$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\widehat f(x)\big|_{x=0^+} &= \frac{\pi}{3} \\
\frac{\partial^1 (\widehat f(x)\quad (x>0))}{\partial x^1} &= - \pi e^{- 3 x} \\
\left.\frac{\partial^1 (\widehat f(x)\quad (x>0))}{\partial x^1}\right|_{x=0^+} &= - \pi
\end{aligned}
$$
:::

:::example[異なる二つの尺度 1 と 2]
実軸上の四次分母を積分する。半径 $R$ を二つの尺度より大きく取り、上半平面で閉じる。

実軸上で分母は正なので、経路をよける小円弧は要らない。複素平面では分母を四つの一次因子に分ける。

$$
\begin{aligned}
f(z) &= \frac{1}{\left(z^{2} + 1\right) \left(z^{2} + 4\right)} \\
z_a &= i \\
z_b &= 2 i
\end{aligned}
$$

上半平面で拾う二つの極は単純である。分母から対応する因子を除き、それぞれの点で評価する。

$$
\begin{aligned}
\operatorname{Res}_{z=z_a}f &= - \frac{i}{6} \\
\operatorname{Res}_{z=z_b}f &= \frac{i}{12}
\end{aligned}
$$

大円弧上では分母が $R^4$ の大きさになり、円弧の長さは $\pi R$ なので、円弧積分は $O(R^{-3})$ になる。

$$
\begin{aligned}
2\pi i(\operatorname{Res}_{z_a}f+\operatorname{Res}_{z_b}f) &= \frac{\pi}{6} \\
\frac{\pi}{b^2-a^2}\left(\frac{1}{a}-\frac{1}{b}\right) &= \frac{\pi}{6}
\end{aligned}
$$

部分分数に分けて実軸上で積分しても同じ値になる。この例では二つの留数に異なる符号が付き、和を先に整理できる。
:::
:::example[異なる二つの尺度 1 と 3]
実軸上の四次分母を積分する。半径 $R$ を二つの尺度より大きく取り、上半平面で閉じる。

実軸上で分母は正なので、経路をよける小円弧は要らない。複素平面では分母を四つの一次因子に分ける。

$$
\begin{aligned}
f(z) &= \frac{1}{\left(z^{2} + 1\right) \left(z^{2} + 9\right)} \\
z_a &= i \\
z_b &= 3 i
\end{aligned}
$$

上半平面で拾う二つの極は単純である。分母から対応する因子を除き、それぞれの点で評価する。

$$
\begin{aligned}
\operatorname{Res}_{z=z_a}f &= - \frac{i}{16} \\
\operatorname{Res}_{z=z_b}f &= \frac{i}{48}
\end{aligned}
$$

大円弧上では分母が $R^4$ の大きさになり、円弧の長さは $\pi R$ なので、円弧積分は $O(R^{-3})$ になる。

$$
\begin{aligned}
2\pi i(\operatorname{Res}_{z_a}f+\operatorname{Res}_{z_b}f) &= \frac{\pi}{12} \\
\frac{\pi}{b^2-a^2}\left(\frac{1}{a}-\frac{1}{b}\right) &= \frac{\pi}{12}
\end{aligned}
$$

部分分数に分けて実軸上で積分しても同じ値になる。この例では二つの留数に異なる符号が付き、和を先に整理できる。
:::
:::example[異なる二つの尺度 2 と 3]
実軸上の四次分母を積分する。半径 $R$ を二つの尺度より大きく取り、上半平面で閉じる。

実軸上で分母は正なので、経路をよける小円弧は要らない。複素平面では分母を四つの一次因子に分ける。

$$
\begin{aligned}
f(z) &= \frac{1}{\left(z^{2} + 4\right) \left(z^{2} + 9\right)} \\
z_a &= 2 i \\
z_b &= 3 i
\end{aligned}
$$

上半平面で拾う二つの極は単純である。分母から対応する因子を除き、それぞれの点で評価する。

$$
\begin{aligned}
\operatorname{Res}_{z=z_a}f &= - \frac{i}{20} \\
\operatorname{Res}_{z=z_b}f &= \frac{i}{30}
\end{aligned}
$$

大円弧上では分母が $R^4$ の大きさになり、円弧の長さは $\pi R$ なので、円弧積分は $O(R^{-3})$ になる。

$$
\begin{aligned}
2\pi i(\operatorname{Res}_{z_a}f+\operatorname{Res}_{z_b}f) &= \frac{\pi}{30} \\
\frac{\pi}{b^2-a^2}\left(\frac{1}{a}-\frac{1}{b}\right) &= \frac{\pi}{30}
\end{aligned}
$$

部分分数に分けて実軸上で積分しても同じ値になる。この例では二つの留数に異なる符号が付き、和を先に整理できる。
:::
:::example[異なる二つの尺度 2 と 5]
実軸上の四次分母を積分する。半径 $R$ を二つの尺度より大きく取り、上半平面で閉じる。

実軸上で分母は正なので、経路をよける小円弧は要らない。複素平面では分母を四つの一次因子に分ける。

$$
\begin{aligned}
f(z) &= \frac{1}{\left(z^{2} + 4\right) \left(z^{2} + 25\right)} \\
z_a &= 2 i \\
z_b &= 5 i
\end{aligned}
$$

上半平面で拾う二つの極は単純である。分母から対応する因子を除き、それぞれの点で評価する。

$$
\begin{aligned}
\operatorname{Res}_{z=z_a}f &= - \frac{i}{84} \\
\operatorname{Res}_{z=z_b}f &= \frac{i}{210}
\end{aligned}
$$

大円弧上では分母が $R^4$ の大きさになり、円弧の長さは $\pi R$ なので、円弧積分は $O(R^{-3})$ になる。

$$
\begin{aligned}
2\pi i(\operatorname{Res}_{z_a}f+\operatorname{Res}_{z_b}f) &= \frac{\pi}{70} \\
\frac{\pi}{b^2-a^2}\left(\frac{1}{a}-\frac{1}{b}\right) &= \frac{\pi}{70}
\end{aligned}
$$

部分分数に分けて実軸上で積分しても同じ値になる。この例では二つの留数に異なる符号が付き、和を先に整理できる。
:::

## 主値積分と Sokhotski–Plemelj の公式

### 極をどちらから避けるか

:::definition[主値]{#def-mp-complex-pv}
有限の点 $a$ に単純極を持つ実積分で、極の両側から同じ幅 $\epsilon$ の区間を除き、その幅を零にした極限を :term[コーシーの主値]{en="Cauchy principal value" yomi="こーしーのしゅち"} と呼ぶ。無限遠でも必要なら対称に切る。切り方は極限とともに指定する。
:::

:::theorem[Sokhotski–Plemelj の公式]{#thm-mp-plemelj}
滑らかで台がコンパクトな $\varphi$ に対して
$$
\lim_{\epsilon\downarrow0}\int_\R\frac{\varphi(x)}{x-a\pm i\epsilon}\dd x
=\operatorname{PV}\int_\R\frac{\varphi(x)}{x-a}\dd x\mp i\pi\varphi(a).
\label{eq-mp-plemelj}
$$
:::

:::proof
$u=x-a$ とし、分母を実部 $u/(u^2+\epsilon^2)$ と虚部 $\mp\epsilon/(u^2+\epsilon^2)$ に分ける。原点近くで $\varphi(a)$ に対称なカットオフを掛けた分を引くと、実部に残る分は $O(u)$ だから優収束で主値に収束する。引いた分は奇関数の積分で零である。虚部で $u=\epsilon t$ と置くと、$\varphi(a+\epsilon t)/(1+t^2)$ を積分する。優収束と $\int_\R (1+t^2)^{-1}\dd t=\pi$ より結論を得る。
:::

従って Cauchy 型積分 $C(z)=(2\pi i)^{-1}\int \varphi(x)/(x-z)\dd x$ では上側境界値から下側境界値を引くと $\varphi$ になる。輪郭を使う場合、上に小円で迂回するとその小円は時計回りで、単純極の留数の $-i\pi$ 倍に収束する。実軸部分だけの値を求めるときには、その小円の寄与を閉曲線積分から差し引く。

ここでは二種類の上下を区別する。$C(a\pm i0)$ で動くのは $z$ なので、被積分関数の分母は $x-a\mp i0$ になる。従って
$$
C(a\pm i0)=\frac1{2\pi i}\operatorname{PV}\int_\R\frac{\varphi(x)}{x-a}\dd x
\left(\pm\frac12\right)\varphi(a).
$$
元の公式は $\varphi$ が $a$ の近くで $C^1$ で、$\int_{|x-a|\ge1}|\varphi(x)|\,|x-a|^{-1}\dd x<\infty$ を満たす局所可積分な場合にも成り立つ。$\varphi$ を $a$ の近くで $1$ となる滑らかなカットオフで分け、近傍外では $1/|x-a|$ を使って優収束を適用すればよい。

:::example[Dirichlet 積分]{#exm-mp-dirichlet-integral}
$k>0$ で $e^{ikz}/z$ を上半平面で積分し、原点を上に迂回する。輪郭内には極がなく、大円は上の評価で消える。小円の値は $-i\pi$ なので主値は $i\pi$ である。実部の積分は奇関数で零、虚部は偶関数なので $\int_0^\infty\sin(kx)/x\,\dd x=\pi/2$ を得る。無限遠の収束は部分積分でも確認できる。
:::

:::example[輪郭内の極と輪郭上の極が共存する積分]
$k>0$ とし、$e^{ikz}/[z(1+z^2)]$ を原点の上で迂回する。上半平面内の極 $i$ の留数は $-e^{-k}/2$、迂回する原点での留数は $1$ である。大円弧は $|e^{ikz}|\le1$ と三次の減衰により消える。残る輪郭の積分は主値から $i\pi$ を引いた値だから、留数定理により主値は $i\pi(1-e^{-k})$ となる。虚部を取れば $\int_\R \sin(kx)/[x(1+x^2)]\,\dd x=\pi(1-e^{-k})$ である。実積分の被積分関数は原点で連続延長でき、無限遠では絶対可積分である。
:::

### 主値と境界値の実部・虚部

主値は同じ距離だけ除く極限である。例えば連続微分可能な $f$ に対し、原点近傍で $f(x)=f(0)+x g(x)$ と書くと、$f(0)/x$ は対称区間で消える。$g$ には通常の積分を使える。一方 $1/(x+i\epsilon)=(x-i\epsilon)/(x^2+\epsilon^2)$ の虚部には、幅が $\epsilon$ 程度の負のピークがある。ピークの面積は $-\pi$ であり、縮小しても消えない。したがって極を迂回する積分と主値を等置するには、半円の向きと半留数の項を指定する必要がある。

:::example[対称除去と半留数（1）]
原点に極を持つ例を取る。原点付近で定数成分と残りを分ける。 一次の模型から始める。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

原点に掛かる滑らかな因子を選ぶ。

$$
\begin{aligned}
f(x)\big|_{m=1} &= \frac{x + 1}{x^{2} + 1}
\end{aligned}
$$

定数分の寄与を対称区間で消すと、残った偶関数が有限の積分を与える。

$$
\begin{aligned}
\operatorname{PV}\int_\R f(x)/x\,\dd x\big|_{m=1} &= \pi
\end{aligned}
$$

上からの境界値は、$f(0)=1$ なので虚部を $-\pi$ だけ加える。

$$
\begin{aligned}
\int_\R f(x)/(x+i0)\,\dd x\big|_{m=1} &= \pi \left(1 - i\right)
\end{aligned}
$$

観測点 $x$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
f(x)\big|_{x=0} &= 1 \\
\frac{\partial^1 (f(x))}{\partial x^1} &= \frac{x^{2} - 2 x \left(x + 1\right) + 1}{\left(x^{2} + 1\right)^{2}} \\
\left.\frac{\partial^1 (f(x))}{\partial x^1}\right|_{x=0} &= 1
\end{aligned}
$$
:::

:::example[対称除去と半留数（2）]
原点に極を持つ例を取る。原点付近で定数成分と残りを分ける。 二次の模型を考える。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

原点に掛かる滑らかな因子を選ぶ。

$$
\begin{aligned}
f(x)\big|_{m=2} &= \frac{2 x + 1}{x^{2} + 1}
\end{aligned}
$$

定数分の寄与を対称区間で消すと、残った偶関数が有限の積分を与える。

$$
\begin{aligned}
\operatorname{PV}\int_\R f(x)/x\,\dd x\big|_{m=2} &= 2 \pi
\end{aligned}
$$

上からの境界値は、$f(0)=1$ なので虚部を $-\pi$ だけ加える。

$$
\begin{aligned}
\int_\R f(x)/(x+i0)\,\dd x\big|_{m=2} &= \pi \left(2 - i\right)
\end{aligned}
$$

観測点 $x$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
f(x)\big|_{x=0} &= 1 \\
\frac{\partial^1 (f(x))}{\partial x^1} &= \frac{2 \left(- x^{2} - x + 1\right)}{x^{4} + 2 x^{2} + 1} \\
\left.\frac{\partial^1 (f(x))}{\partial x^1}\right|_{x=0} &= 2
\end{aligned}
$$
:::

:::example[対称除去と半留数（3）]
原点に極を持つ例を取る。原点付近で定数成分と残りを分ける。 三次の模型を考える。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

原点に掛かる滑らかな因子を選ぶ。

$$
\begin{aligned}
f(x)\big|_{m=3} &= \frac{3 x + 1}{x^{2} + 1}
\end{aligned}
$$

定数分の寄与を対称区間で消すと、残った偶関数が有限の積分を与える。

$$
\begin{aligned}
\operatorname{PV}\int_\R f(x)/x\,\dd x\big|_{m=3} &= 3 \pi
\end{aligned}
$$

上からの境界値は、$f(0)=1$ なので虚部を $-\pi$ だけ加える。

$$
\begin{aligned}
\int_\R f(x)/(x+i0)\,\dd x\big|_{m=3} &= \pi \left(3 - i\right)
\end{aligned}
$$

観測点 $x$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
f(x)\big|_{x=0} &= 1 \\
\frac{\partial^1 (f(x))}{\partial x^1} &= \frac{- 3 x^{2} - 2 x + 3}{x^{4} + 2 x^{2} + 1} \\
\left.\frac{\partial^1 (f(x))}{\partial x^1}\right|_{x=0} &= 3
\end{aligned}
$$
:::

## 分岐切断と多価関数の積分

### 切断の両岸で位相を読む

非整数冪 $z^\alpha$ は $\log z=\log|z|+i\arg z$ の枝を固定して定義する。正の実軸に切断を置くとき、$0<\arg z<2\pi$ を選べば上の岸では位相が $1$、下の岸では $e^{2\pi i\alpha}$ となる。上下の向きも逆である。

:::example[鍵穴型輪郭]
$0<a<1$ とし $z^{a-1}/(1+z)$ を正の実軸で切る。原点の小円は $O(\epsilon^a)$、外円は $O(R^{a-1})$ で消える。実軸からの和は $(1-e^{2\pi ia})I$、$I=\int_0^\infty x^{a-1}/(1+x)\dd x$ である。極 $-1$ の留数は $e^{i\pi(a-1)}$ なので
$$
I=\frac{2\pi i e^{i\pi(a-1)}}{1-e^{2\pi ia}}=\frac{\pi}{\sin\pi a}.
\label{eq-mp-keyhole-beta}
$$
$a=1/2$ では $\pi$ となる。$x=t^2$ なら $I=2\int_0^\infty(1+t^2)^{-1}\dd t$ になる。
:::

対数因子を得るにはパラメータで微分してもよい。$a$ が $(0,1)$ のコンパクトな部分区間を動く限り、両端で $|\log x|$ を掛けた関数の積分が収束して優収束が使える。従って $\int_0^\infty x^{a-1}\log(x)/(1+x)\dd x=-\pi^2\cos(\pi a)/\sin^2(\pi a)$ である。$a=1/2$ で零となることは $x\mapsto1/x$ で符号が入れ替わることからも分かる。

:::example[負の実軸を切った平方根]
$\arg z\in(-\pi,\pi)$ を選ぶと $(-x+i0)^{1/2}=i\sqrt x$、$(-x-i0)^{1/2}=-i\sqrt x$（$x>0$）である。従って $\sqrt z$ の境界値の差は $2i\sqrt x$ となる。例えば $x=4$ なら両境界値は $\pm2i$ で差は $4i$ である。これは正実軸の切断のときの上岸・下岸の比と組になり方が異なる。どちらの枝でも一周した後の倍率は $-1$ である。
:::

### 分岐点のまわりでどこまで切断するか

$z^\alpha$ を使う前に $\log z$ の偏角の区間を決める。正の実軸を切り、$0<\arg z<2\pi$ を選ぶと、上岸の $x^\alpha$ は下岸で $e^{2\pi i\alpha}x^\alpha$ になる。輪郭が下岸を反対向きに進むので、両者の差に因子 $1-e^{2\pi i\alpha}$ が現れる。$z^{\alpha-1}/(1+z)^b$ で $b$ が正整数なら $-1$ は $b$ 位の極であり、小円弧が $O(\epsilon^{\Re\alpha})$、大円弧が $O(R^{\Re\alpha-b})$ となる。従って両円弧が消える範囲は $0<\Re\alpha<b$ である。この範囲でも $\alpha$ が整数なら係数 $1-e^{2\pi i\alpha}$ が零になるため、非整数の $\alpha$ で計算した式を連続延長する。$b$ を非整数にする場合は $-1$ の枝も指定し、その切断に沿う寄与を含める。

:::example[同じ切断に対する分母の次数（1）]
平方根には正実軸で切断を入れ、$-1$ の極の次数を変える。 一次の模型から始める。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

正の実軸の上岸に現れる被積分関数を $u$ で書く。

$$
\begin{aligned}
q(u)\big|_{m=1} &= \frac{1}{\sqrt{u} \left(u + 1\right)}
\end{aligned}
$$

切断の下岸では平方根の符号が変わるため両岸の差は上岸の二倍になる。

$$
\begin{aligned}
q_+(u)-q_-(u)\big|_{m=1} &= \frac{2}{\sqrt{u} \left(u + 1\right)}
\end{aligned}
$$

実軸上で $u=s/(1-s)$ と置けば Beta 積分になり、正の値を与える。

$$
\begin{aligned}
\int_0^\infty q(u)\,\dd u\big|_{m=1} &= \pi
\end{aligned}
$$

観測点 $u$=$1$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
q(u)\big|_{u=1} &= \frac{1}{2} \\
\frac{\partial^1 (q(u))}{\partial u^1} &= \frac{- 3 u - 1}{2 u^{\frac{3}{2}} \left(u^{2} + 2 u + 1\right)} \\
\left.\frac{\partial^1 (q(u))}{\partial u^1}\right|_{u=1} &= - \frac{1}{2}
\end{aligned}
$$

観測点 $u$=$1$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
q_+(u)-q_-(u)\big|_{u=1} &= 1 \\
\frac{\partial^1 (q_+(u)-q_-(u))}{\partial u^1} &= \frac{- 3 u - 1}{u^{\frac{3}{2}} \left(u^{2} + 2 u + 1\right)} \\
\left.\frac{\partial^1 (q_+(u)-q_-(u))}{\partial u^1}\right|_{u=1} &= -1
\end{aligned}
$$
:::

:::example[同じ切断に対する分母の次数（2）]
平方根には正実軸で切断を入れ、$-1$ の極の次数を変える。 二次の模型を考える。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

正の実軸の上岸に現れる被積分関数を $u$ で書く。

$$
\begin{aligned}
q(u)\big|_{m=2} &= \frac{1}{\sqrt{u} \left(u + 1\right)^{2}}
\end{aligned}
$$

切断の下岸では平方根の符号が変わるため両岸の差は上岸の二倍になる。

$$
\begin{aligned}
q_+(u)-q_-(u)\big|_{m=2} &= \frac{2}{\sqrt{u} \left(u + 1\right)^{2}}
\end{aligned}
$$

実軸上で $u=s/(1-s)$ と置けば Beta 積分になり、正の値を与える。

$$
\begin{aligned}
\int_0^\infty q(u)\,\dd u\big|_{m=2} &= \frac{\pi}{2}
\end{aligned}
$$

観測点 $u$=$1$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
q(u)\big|_{u=1} &= \frac{1}{4} \\
\frac{\partial^1 (q(u))}{\partial u^1} &= \frac{- 5 u - 1}{2 u^{\frac{3}{2}} \left(u + 1\right)^{3}} \\
\left.\frac{\partial^1 (q(u))}{\partial u^1}\right|_{u=1} &= - \frac{3}{8}
\end{aligned}
$$

観測点 $u$=$1$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
q_+(u)-q_-(u)\big|_{u=1} &= \frac{1}{2} \\
\frac{\partial^1 (q_+(u)-q_-(u))}{\partial u^1} &= \frac{- 5 u - 1}{u^{\frac{3}{2}} \left(u + 1\right)^{3}} \\
\left.\frac{\partial^1 (q_+(u)-q_-(u))}{\partial u^1}\right|_{u=1} &= - \frac{3}{4}
\end{aligned}
$$
:::

:::example[同じ切断に対する分母の次数（3）]
平方根には正実軸で切断を入れ、$-1$ の極の次数を変える。 三次の模型を考える。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

正の実軸の上岸に現れる被積分関数を $u$ で書く。

$$
\begin{aligned}
q(u)\big|_{m=3} &= \frac{1}{\sqrt{u} \left(u + 1\right)^{3}}
\end{aligned}
$$

切断の下岸では平方根の符号が変わるため両岸の差は上岸の二倍になる。

$$
\begin{aligned}
q_+(u)-q_-(u)\big|_{m=3} &= \frac{2}{\sqrt{u} \left(u + 1\right)^{3}}
\end{aligned}
$$

実軸上で $u=s/(1-s)$ と置けば Beta 積分になり、正の値を与える。

$$
\begin{aligned}
\int_0^\infty q(u)\,\dd u\big|_{m=3} &= \frac{3 \pi}{8}
\end{aligned}
$$

観測点 $u$=$1$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
q(u)\big|_{u=1} &= \frac{1}{8} \\
\frac{\partial^1 (q(u))}{\partial u^1} &= \frac{- 7 u - 1}{2 u^{\frac{3}{2}} \left(u + 1\right)^{4}} \\
\left.\frac{\partial^1 (q(u))}{\partial u^1}\right|_{u=1} &= - \frac{1}{4}
\end{aligned}
$$

観測点 $u$=$1$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
q_+(u)-q_-(u)\big|_{u=1} &= \frac{1}{4} \\
\frac{\partial^1 (q_+(u)-q_-(u))}{\partial u^1} &= \frac{- 7 u - 1}{u^{\frac{3}{2}} \left(u + 1\right)^{4}} \\
\left.\frac{\partial^1 (q_+(u)-q_-(u))}{\partial u^1}\right|_{u=1} &= - \frac{1}{2}
\end{aligned}
$$
:::

## 解析接続と Mittag-Leffler 展開

### 局所の表示をつなぐ

:term[解析接続]{en="analytic continuation" yomi="かいせきせつぞく"} は、重なりで正則関数が一致するように定義域を延ばす操作である。連結な重なりで開集合上に一致があれば一致の定理で一意である。ただし穴の周りを通る経路の違いで枝が変わり得る。たとえば $\log z$ を原点の周りに一周接続すると、導関数の積分 $\int \dd z/z=2\pi i$ だけ変わる。

:::proposition[余接の Mittag–Leffler 展開]{#prop-mp-cotangent-expansion}
$z\notin\Z$ に対し
$$
\pi\cot\pi z=\frac1z+\sum_{n=1}^\infty\frac{2z}{z^2-n^2}.
\label{eq-mp-cotangent-expansion}
$$
:::

:::proof
$\pi\cot(\pi\zeta)/(\zeta^2-z^2)$ を辺の座標が $\pm(N+1/2)$ の正方形で積分する。縦辺の余接は絶対値 $|\tanh(\pi\Im\zeta)|\le1$、横辺でも指数表示から一様有界である。分母は境界上で $O(N^{-2})$ の逆数を与えるので積分は零に近づく。整数 $n$ での留数は $1/(n^2-z^2)$、$\pm z$ の留数の和は $\pi\cot(\pi z)/z$ である。和を整理すれば結論を得る。各項がコンパクト集合上一様に $O(n^{-2})$ なので局所一様収束する。
:::

一般の Mittag–Leffler 構成でも主部をただ足すと発散し得る。有限平面内で離散な極 $a_n$ の主部を $P_n$ とし、$|a_n|\to\infty$ の順に並べる。原点を含み半径が無限に増える円板を各 $n$ に対して $|a_n|$ より小さく選び、その円板で $P_n$ を Taylor 多項式 $Q_n$ により誤差 $2^{-n}$ 以下で近似する。有限個の初項はそのまま置く。$\sum(P_n-Q_n)$ は極を避けて局所一様収束し、指定の主部を持つ。別の解との差はすべての極で可除だから整関数である。

:::example[逆二乗和への応用]
$\pi\cot\pi z=1/z-(\pi^2/3)z+O(z^3)$ は $\sin$ と $\cos$ の Taylor 展開の商から得る。式 :eqref[eq-mp-cotangent-expansion] の $z$ の係数は $-2\sum_{n\ge1}1/n^2$ である。項別微分は $|z|<1$ の部分円板で一様収束するため許され、和は $\pi^2/6$ となる。
:::

:::example[実軸上に極のない格子和]
式 :eqref[eq-mp-cotangent-expansion] に $z=ia$（$a>0$）を入れる。指数表示から $\cot(i\pi a)=-i\coth(\pi a)$ である。両辺から $-ia$ を因子として外すと
$$
\sum_{n\in\Z}\frac1{n^2+a^2}=\frac\pi a\coth(\pi a)
$$
となる。左の $n=0$ は $a^{-2}$ で、正負の対の寄与は $2/(n^2+a^2)$ である。例えば $a=1$ では和が $\pi\coth\pi$ になる。$a$ が小さいと $n=0$ の項が支配し、$\coth(\pi a)\sim1/(\pi a)$ と整合する。$a>0$ のコンパクトな区間ではこの和を項別微分でき、二乗した分母を持つ格子和も計算できる。
:::

### 主部を引いて正規収束させる

無限個の極を足す場合、各極で正しい主部を持つだけでは和の収束が保証されない。離れた点 $n$ での単純極について $1/(z-n)+1/n=z/[n(z-n)]$ とすれば、任意の固定したコンパクト集合上で大きい $n$ の項は $O(n^{-2})$ になる。これは極近傍の特異項を変えずに遠方の寄与を抑えるための差し引きである。極を除いたコンパクト集合で一様収束するので、和は正則であり微分も局所で項別にできる。与えられた関数と和の差は整関数になる。その整関数が零と示すには、周期性と有界性などの追加条件が必要である。

:::example[一周期に並ぶ極の和（1）]
実数パラメータをずらして格子和を求める。差し引き済みの展開を微分して二乗分母を作る。 単位の周期を取る。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

周期を縮めた余接を取る。$x$ の極は間隔 $1/m$ で並ぶ。

$$
\begin{aligned}
C(x)\big|_{m=1} &= \pi \cot{\left(\pi x \right)}
\end{aligned}
$$

各単純極の主部を微分し、全体の符号を反転すると絶対収束する二乗分母の和になる。

$$
\begin{aligned}
-C^\prime(x)\big|_{m=1} &= \frac{\pi^{2}}{\sin^{2}{\left(\pi x \right)}}
\end{aligned}
$$

観測点 $x$=$\frac{1}{4}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
C(x)\big|_{x=\frac{1}{4}} &= \pi \\
\frac{\partial^1 (C(x))}{\partial x^1} &= - \frac{\pi^{2}}{\sin^{2}{\left(\pi x \right)}} \\
\left.\frac{\partial^1 (C(x))}{\partial x^1}\right|_{x=\frac{1}{4}} &= - 2 \pi^{2}
\end{aligned}
$$

観測点 $x$=$\frac{1}{4}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
-C^\prime(x)\big|_{x=\frac{1}{4}} &= 2 \pi^{2} \\
\frac{\partial^1 (-C^\prime(x))}{\partial x^1} &= - \frac{2 \pi^{3} \cos{\left(\pi x \right)}}{\sin^{3}{\left(\pi x \right)}} \\
\left.\frac{\partial^1 (-C^\prime(x))}{\partial x^1}\right|_{x=\frac{1}{4}} &= - 4 \pi^{3}
\end{aligned}
$$
:::

:::example[一周期に並ぶ極の和（2）]
実数パラメータをずらして格子和を求める。差し引き済みの展開を微分して二乗分母を作る。 周期を半分にする。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

周期を縮めた余接を取る。$x$ の極は間隔 $1/m$ で並ぶ。

$$
\begin{aligned}
C(x)\big|_{m=2} &= 2 \pi \cot{\left(2 \pi x \right)}
\end{aligned}
$$

各単純極の主部を微分し、全体の符号を反転すると絶対収束する二乗分母の和になる。

$$
\begin{aligned}
-C^\prime(x)\big|_{m=2} &= \frac{4 \pi^{2}}{\sin^{2}{\left(2 \pi x \right)}}
\end{aligned}
$$

観測点 $x$=$\frac{1}{4}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
C(x)\big|_{x=\frac{1}{4}} &= 0 \\
\frac{\partial^1 (C(x))}{\partial x^1} &= - \frac{4 \pi^{2}}{\sin^{2}{\left(2 \pi x \right)}} \\
\left.\frac{\partial^1 (C(x))}{\partial x^1}\right|_{x=\frac{1}{4}} &= - 4 \pi^{2}
\end{aligned}
$$

観測点 $x$=$\frac{1}{4}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
-C^\prime(x)\big|_{x=\frac{1}{4}} &= 4 \pi^{2} \\
\frac{\partial^1 (-C^\prime(x))}{\partial x^1} &= - \frac{16 \pi^{3} \cos{\left(2 \pi x \right)}}{\sin^{3}{\left(2 \pi x \right)}} \\
\left.\frac{\partial^1 (-C^\prime(x))}{\partial x^1}\right|_{x=\frac{1}{4}} &= 0
\end{aligned}
$$
:::

:::example[一周期に並ぶ極の和（3）]
実数パラメータをずらして格子和を求める。差し引き済みの展開を微分して二乗分母を作る。 格子の間隔を三分の一にする。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

周期を縮めた余接を取る。$x$ の極は間隔 $1/m$ で並ぶ。

$$
\begin{aligned}
C(x)\big|_{m=3} &= 3 \pi \cot{\left(3 \pi x \right)}
\end{aligned}
$$

各単純極の主部を微分し、全体の符号を反転すると絶対収束する二乗分母の和になる。

$$
\begin{aligned}
-C^\prime(x)\big|_{m=3} &= \frac{9 \pi^{2}}{\sin^{2}{\left(3 \pi x \right)}}
\end{aligned}
$$

観測点 $x$=$\frac{1}{4}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
C(x)\big|_{x=\frac{1}{4}} &= - 3 \pi \\
\frac{\partial^1 (C(x))}{\partial x^1} &= - \frac{9 \pi^{2}}{\sin^{2}{\left(3 \pi x \right)}} \\
\left.\frac{\partial^1 (C(x))}{\partial x^1}\right|_{x=\frac{1}{4}} &= - 18 \pi^{2}
\end{aligned}
$$

観測点 $x$=$\frac{1}{4}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
-C^\prime(x)\big|_{x=\frac{1}{4}} &= 18 \pi^{2} \\
\frac{\partial^1 (-C^\prime(x))}{\partial x^1} &= - \frac{54 \pi^{3} \cos{\left(3 \pi x \right)}}{\sin^{3}{\left(3 \pi x \right)}} \\
\left.\frac{\partial^1 (-C^\prime(x))}{\partial x^1}\right|_{x=\frac{1}{4}} &= 108 \pi^{3}
\end{aligned}
$$
:::

### 与えられた主部を持つ関数を構成する

:::theorem[定めた主部を持つ有理型関数]{#thm-mp-mittag-leffler}
複素平面内に集積点のない相異なる点 $a_j$ と、それぞれ有限個の負の冪からなる主部 $P_j(z)=\sum_{\ell=1}^{d_j}c_{j\ell}(z-a_j)^{-\ell}$ を与える。この主部をすべての $a_j$ に持ち、それ以外に極を持たない有理型関数が存在する。同じ主部を持つ二関数の差は整関数である。
:::

:::proof
有限個の場合は主部を足せばよい。無限個の場合、有限個の点を先に取り除き、残りを $|a_j|\to\infty$ となる順に並べる。$r_j=|a_j|/2$ とすれば、$P_j$ の原点まわりの Taylor 級数は閉円板 $|z|\le r_j$ で一様に収束する。従って Taylor 多項式 $Q_j$ を十分高い次数まで取り、$\sup_{|z|\le r_j}|P_j(z)-Q_j(z)|\le2^{-j}$ とできる。任意の閉円板は十分大きい $j$ のすべてについて $|z|\le r_j$ に入るため、和 $\sum_j(P_j-Q_j)$ の尾はその閉円板上で一様収束する。そこで、極を除いたコンパクト集合では各項の和が正則である。各 $a_j$ ではその点に対応する項以外は近傍で正則に足せるので、負の冪は $P_j$ と一致する。最後に取り除いた有限個の主部を足せばよい。二関数の差はすべての極で負の冪が消えるため、除去可能特異点を埋めると整関数になる。
:::

例えば $a_j=2^j$、主部 $P_j(z)=2^j/(z-2^j)$ とする。各項は固定した $z$ で $-1$ に近づくので、主部をそのまま足すと級数の必要条件すら満たさない。定数項を引くと $P_j(z)+1=z/(z-2^j)$ になる。$2^j\ge2R$ なら $|z|\le R$ 上で絶対値は $2R/2^j$ 以下であり、幾何級数で抑えられる。したがってこの例では各極で Taylor の定数項を引くだけで定理の構成ができる。$z=0$ ではすべての補正後の項が零となるので、選んだ整関数部分の規格化も確認できる。

## 分散関係と Kramers–Kronig の関係

### 因果応答の境界値

入力 $E(t)$ に対し、線形で時間並進不変な応答を $P(t)=\int_0^\infty h(s)E(t-s)\dd s$ とする。$E(t)=e^{-i\omega t}$ を用いると、周波数の乗数は $\chi(z)=\int_0^\infty h(t)e^{izt}\dd t$ となる。これは負指数 Fourier 変換の周波数 $k=-z$ に当たる。$h\in L^1(0,\infty)$ なら上半平面で正則で、閉上半平面上で連続である。正則性は任意の内部コンパクト集合上で $t e^{-(\Im z)t}$ が有界であることから従う。

:::proposition[Kramers–Kronig 関係]{#prop-mp-dispersion}
$\chi$ が上半平面で正則、閉半平面上で連続で、実軸上で局所 $C^1$、無限遠で一様に $O(1/|z|)$ とする。$\chi=\chi_R+i\chi_I$ と書くと
$$
\begin{aligned}
\chi_R(\omega)&=\frac1\pi\operatorname{PV}\int_\R\frac{\chi_I(v)}{v-\omega}\dd v,\\
\chi_I(\omega)&=-\frac1\pi\operatorname{PV}\int_\R\frac{\chi_R(v)}{v-\omega}\dd v.
\end{aligned}
$$
:::

:::proof
上半円を閉じて Cauchy の公式を使うと $\chi(z)=(2\pi i)^{-1}\int \chi(v)/(v-z)\dd v$ となる。円弧の寄与は $O(1/R)$ である。$z=\omega+i\epsilon$ を実軸に近づけ、式 :eqref[eq-mp-plemelj] を $\chi$ に適用する。遠方の尾は $O(v^{-2})$ なのでコンパクト台の条件は切り落とせる。右辺の局所項は $\chi(\omega)/2$ である。移項して実部と虚部を比較すればよい。
:::

:::example[指数緩和]
$h(t)=e^{-\gamma t}$（$t\ge0,\ \gamma>0$）では $\chi(\omega)=1/(\gamma-i\omega)$ である。極は $-i\gamma$ で上半平面を避ける。$\omega=0$ の分散積分は $\pi^{-1}\int_\R (\gamma^2+v^2)^{-1}\dd v=1/\gamma$ で、$\chi_R(0)$ と一致する。実時間応答が実数なら $\chi(-\omega)=\overline{\chi(\omega)}$ なので偶奇性も確認できる。
:::

定数 $\chi_\infty$ を差し引いた応答が閉上半平面で一様に $O(1/|z|)$ となる場合は、$\chi-\chi_\infty$ に公式を適用する。例えば入力に比例する瞬時応答は周波数に依らないため、この定数に含めて扱う。

### 減衰条件を引き算で改善する

$\chi$ は閉上半平面で連続、上半平面で正則かつ有界とする。実軸上で局所 $C^1$ とし、基準周波数の近くでは上半平面側の一階微分も境界まで連続とする。基準周波数 $\omega_0$ を選び、
$F(z)=[\chi(z)-\chi(\omega_0)]/(z-\omega_0)$ と置く。境界の値は $\omega_0$ で導関数により連続延長できる。$F$ は無限遠で $O(1/|z|)$ である。$\omega\ne\omega_0$ で上半平面の Cauchy の式を使って境界値を取ると
$$
\chi(\omega)-\chi(\omega_0)
=\frac{\omega-\omega_0}{\pi i}\operatorname{PV}\int_\R
\frac{\chi(v)-\chi(\omega_0)}{(v-\omega_0)(v-\omega)}\dd v
$$
を得る。基準周波数の近くでは、$\chi(v)-\chi(\omega_0)$ が一階の零を持つため、そこでの分母は可除である。例えば $\chi(z)=2+(1-iz)^{-1}$ のとき $\omega_0=0$ とすれば $F(z)=i/(1-iz)$ なので、定数の瞬時応答 $2$ は差から自動的に消える。反対に、差だけの関係から絶対的な静的応答は決まらない。

### 実応答では正の周波数だけを積分できる

$h(t)$ が実なら $\chi_R$ は偶、$\chi_I$ は奇である。式の負の部分を $v\mapsto-v$ で持ち上げると、$\omega>0$ で
$$
\chi_R(\omega)=\frac2\pi\operatorname{PV}\int_0^\infty
\frac{v\chi_I(v)}{v^2-\omega^2}\dd v,\qquad
\chi_I(\omega)=-\frac{2\omega}{\pi}\operatorname{PV}\int_0^\infty
\frac{\chi_R(v)}{v^2-\omega^2}\dd v .
$$
例えば静的応答は、$\chi_I(v)=O(v)$ が原点で成り立てば $\chi_R(0)=(2/\pi)\int_0^\infty\chi_I(v)\dd v/v$ である。減衰係数 $\gamma=2\,\mathrm{s^{-1}}$ の単純な緩和核なら、積分の中身は $(4+v^2)^{-1}$ で、結果は $1/2\,\mathrm{s}$ になる。実際の物理量の単位は $h$ への比例係数も含めて決まる。

### 低周波数の引き算と測定量

実時間の応答が可積分で、因果的であるとき、正の虚部への延長が得られる。Hilbert 変換で実部を復元するにはさらに大円弧を評価する。実軸の高周波で定数に近いと分かっているなら、その定数を引いた量で先に考えられる。引く定数が未知の場合は $\chi(z)-\chi(z_0)$ の商を考えれば、$z-z_0$ が分母に加わる。この差を取ると、高周波の影響が一段小さい重みになる。その代わり $\chi(z_0)$ は独立の入力として必要になる。静的応答では $z_0=0$ が便利であるが、原点での連続性を確認する。

:::example[時間の一次因子を含む緩和（1）]
核 $h(t)=te^{-mt}$（$t\ge0$）を周波数へ移す。以下では実周波数を $x$ とし、$\Im\chi(x)/x$ は原点に連続に延長する。 一次の模型から始める。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

積分で $t$ が掛かるため、複素周波数の極は二位になる。

$$
\begin{aligned}
\chi(x)\big|_{m=1} &= \frac{1}{\left(i x - 1\right)^{2}}
\end{aligned}
$$

実部を取り出す。零周波数では分母の尺度だけが残る。

$$
\begin{aligned}
\Re\chi(x)\big|_{m=1} &= \frac{1 - x^{2}}{\left(x^{2} + 1\right)^{2}}
\end{aligned}
$$

虚部は奇関数であり、$x$ で割れば偶関数になる。

$$
\begin{aligned}
\Im\chi(x)/x\big|_{m=1} &= \frac{2}{\left(x^{2} + 1\right)^{2}}
\end{aligned}
$$

観測点 $x$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\chi(x)\big|_{x=0} &= 1 \\
\frac{\partial^1 (\chi(x))}{\partial x^1} &= - \frac{2 i}{\left(i x - 1\right)^{3}} \\
\left.\frac{\partial^1 (\chi(x))}{\partial x^1}\right|_{x=0} &= 2 i \\
\frac{\partial^2 (\chi(x))}{\partial x^2} &= - \frac{6}{\left(i x - 1\right)^{4}} \\
\left.\frac{\partial^2 (\chi(x))}{\partial x^2}\right|_{x=0} &= -6
\end{aligned}
$$

観測点 $x$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\Re\chi(x)\big|_{x=0} &= 1 \\
\frac{\partial^1 (\Re\chi(x))}{\partial x^1} &= \frac{2 x \left(x^{2} - 3\right)}{\left(x^{2} + 1\right)^{3}} \\
\left.\frac{\partial^1 (\Re\chi(x))}{\partial x^1}\right|_{x=0} &= 0 \\
\frac{\partial^2 (\Re\chi(x))}{\partial x^2} &= \frac{6 \left(- x^{4} + 6 x^{2} - 1\right)}{x^{8} + 4 x^{6} + 6 x^{4} + 4 x^{2} + 1} \\
\left.\frac{\partial^2 (\Re\chi(x))}{\partial x^2}\right|_{x=0} &= -6
\end{aligned}
$$

観測点 $x$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\Im\chi(x)/x\big|_{x=0} &= 2 \\
\frac{\partial^1 (\Im\chi(x)/x)}{\partial x^1} &= - \frac{8 x}{\left(x^{2} + 1\right)^{3}} \\
\left.\frac{\partial^1 (\Im\chi(x)/x)}{\partial x^1}\right|_{x=0} &= 0 \\
\frac{\partial^2 (\Im\chi(x)/x)}{\partial x^2} &= \frac{8 \left(5 x^{2} - 1\right)}{\left(x^{2} + 1\right)^{4}} \\
\left.\frac{\partial^2 (\Im\chi(x)/x)}{\partial x^2}\right|_{x=0} &= -8
\end{aligned}
$$
:::

:::example[時間の一次因子を含む緩和（2）]
核 $h(t)=te^{-mt}$（$t\ge0$）を周波数へ移す。以下では実周波数を $x$ とし、$\Im\chi(x)/x$ は原点に連続に延長する。 二次の模型を考える。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

積分で $t$ が掛かるため、複素周波数の極は二位になる。

$$
\begin{aligned}
\chi(x)\big|_{m=2} &= \frac{1}{\left(i x - 2\right)^{2}}
\end{aligned}
$$

実部を取り出す。零周波数では分母の尺度だけが残る。

$$
\begin{aligned}
\Re\chi(x)\big|_{m=2} &= \frac{4 - x^{2}}{\left(x^{2} + 4\right)^{2}}
\end{aligned}
$$

虚部は奇関数であり、$x$ で割れば偶関数になる。

$$
\begin{aligned}
\Im\chi(x)/x\big|_{m=2} &= \frac{4}{\left(x^{2} + 4\right)^{2}}
\end{aligned}
$$

観測点 $x$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\chi(x)\big|_{x=0} &= \frac{1}{4} \\
\frac{\partial^1 (\chi(x))}{\partial x^1} &= - \frac{2 i}{\left(i x - 2\right)^{3}} \\
\left.\frac{\partial^1 (\chi(x))}{\partial x^1}\right|_{x=0} &= \frac{i}{4} \\
\frac{\partial^2 (\chi(x))}{\partial x^2} &= - \frac{6}{\left(i x - 2\right)^{4}} \\
\left.\frac{\partial^2 (\chi(x))}{\partial x^2}\right|_{x=0} &= - \frac{3}{8}
\end{aligned}
$$

観測点 $x$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\Re\chi(x)\big|_{x=0} &= \frac{1}{4} \\
\frac{\partial^1 (\Re\chi(x))}{\partial x^1} &= \frac{2 x \left(x^{2} - 12\right)}{\left(x^{2} + 4\right)^{3}} \\
\left.\frac{\partial^1 (\Re\chi(x))}{\partial x^1}\right|_{x=0} &= 0 \\
\frac{\partial^2 (\Re\chi(x))}{\partial x^2} &= \frac{6 \left(- x^{4} + 24 x^{2} - 16\right)}{x^{8} + 16 x^{6} + 96 x^{4} + 256 x^{2} + 256} \\
\left.\frac{\partial^2 (\Re\chi(x))}{\partial x^2}\right|_{x=0} &= - \frac{3}{8}
\end{aligned}
$$

観測点 $x$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\Im\chi(x)/x\big|_{x=0} &= \frac{1}{4} \\
\frac{\partial^1 (\Im\chi(x)/x)}{\partial x^1} &= - \frac{16 x}{\left(x^{2} + 4\right)^{3}} \\
\left.\frac{\partial^1 (\Im\chi(x)/x)}{\partial x^1}\right|_{x=0} &= 0 \\
\frac{\partial^2 (\Im\chi(x)/x)}{\partial x^2} &= \frac{16 \left(5 x^{2} - 4\right)}{\left(x^{2} + 4\right)^{4}} \\
\left.\frac{\partial^2 (\Im\chi(x)/x)}{\partial x^2}\right|_{x=0} &= - \frac{1}{4}
\end{aligned}
$$
:::

:::example[時間の一次因子を含む緩和（3）]
核 $h(t)=te^{-mt}$（$t\ge0$）を周波数へ移す。以下では実周波数を $x$ とし、$\Im\chi(x)/x$ は原点に連続に延長する。 三次の模型を考える。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

積分で $t$ が掛かるため、複素周波数の極は二位になる。

$$
\begin{aligned}
\chi(x)\big|_{m=3} &= \frac{1}{\left(i x - 3\right)^{2}}
\end{aligned}
$$

実部を取り出す。零周波数では分母の尺度だけが残る。

$$
\begin{aligned}
\Re\chi(x)\big|_{m=3} &= \frac{9 - x^{2}}{\left(x^{2} + 9\right)^{2}}
\end{aligned}
$$

虚部は奇関数であり、$x$ で割れば偶関数になる。

$$
\begin{aligned}
\Im\chi(x)/x\big|_{m=3} &= \frac{6}{\left(x^{2} + 9\right)^{2}}
\end{aligned}
$$

観測点 $x$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\chi(x)\big|_{x=0} &= \frac{1}{9} \\
\frac{\partial^1 (\chi(x))}{\partial x^1} &= - \frac{2 i}{\left(i x - 3\right)^{3}} \\
\left.\frac{\partial^1 (\chi(x))}{\partial x^1}\right|_{x=0} &= \frac{2 i}{27} \\
\frac{\partial^2 (\chi(x))}{\partial x^2} &= - \frac{6}{\left(i x - 3\right)^{4}} \\
\left.\frac{\partial^2 (\chi(x))}{\partial x^2}\right|_{x=0} &= - \frac{2}{27}
\end{aligned}
$$

観測点 $x$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\Re\chi(x)\big|_{x=0} &= \frac{1}{9} \\
\frac{\partial^1 (\Re\chi(x))}{\partial x^1} &= \frac{2 x \left(x^{2} - 27\right)}{\left(x^{2} + 9\right)^{3}} \\
\left.\frac{\partial^1 (\Re\chi(x))}{\partial x^1}\right|_{x=0} &= 0 \\
\frac{\partial^2 (\Re\chi(x))}{\partial x^2} &= \frac{6 \left(- x^{4} + 54 x^{2} - 81\right)}{x^{8} + 36 x^{6} + 486 x^{4} + 2916 x^{2} + 6561} \\
\left.\frac{\partial^2 (\Re\chi(x))}{\partial x^2}\right|_{x=0} &= - \frac{2}{27}
\end{aligned}
$$

観測点 $x$=$0$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\Im\chi(x)/x\big|_{x=0} &= \frac{2}{27} \\
\frac{\partial^1 (\Im\chi(x)/x)}{\partial x^1} &= - \frac{24 x}{\left(x^{2} + 9\right)^{3}} \\
\left.\frac{\partial^1 (\Im\chi(x)/x)}{\partial x^1}\right|_{x=0} &= 0 \\
\frac{\partial^2 (\Im\chi(x)/x)}{\partial x^2} &= \frac{24 \left(5 x^{2} - 9\right)}{\left(x^{2} + 9\right)^{4}} \\
\left.\frac{\partial^2 (\Im\chi(x)/x)}{\partial x^2}\right|_{x=0} &= - \frac{8}{243}
\end{aligned}
$$
:::

:::example[緩和率が 1 のときの吸収と静的応答]
因果応答に対し、通常の Fourier 変換と逆符号の $e^{i\omega t}$ を用いて、正の虚部に解析的な周波数応答を作る。

時間の積分を直接計算する。収束を決める実部は複素周波数に正の虚部を加えるほど大きくなる。

$$
\begin{aligned}
h(t) &= e^{- t} \\
\chi(\omega) &= \frac{1}{- i \omega + 1}
\end{aligned}
$$

分母に複素共役を掛けて、分散を表す実部と吸収を表す虚部を分ける。以下の $x$ は実周波数である。

$$
\begin{aligned}
\Re\chi(x) &= \frac{1}{x^{2} + 1} \\
\Im\chi(x) &= \frac{x}{x^{2} + 1} \\
\Im\chi(x)/x &= \frac{1}{x^{2} + 1}
\end{aligned}
$$

最後の比は $x=0$ に連続に延びる。両側で積分し、$\pi$ で割ると静的な実部に一致する。

$$
\begin{aligned}
\frac1\pi\int_\R\frac{\Im\chi(x)}x\,\dd x &= 1 \\
\chi(0) &= 1
\end{aligned}
$$
:::
:::example[緩和率が 2 のときの吸収と静的応答]
因果応答に対し、通常の Fourier 変換と逆符号の $e^{i\omega t}$ を用いて、正の虚部に解析的な周波数応答を作る。

時間の積分を直接計算する。収束を決める実部は複素周波数に正の虚部を加えるほど大きくなる。

$$
\begin{aligned}
h(t) &= e^{- 2 t} \\
\chi(\omega) &= \frac{1}{- i \omega + 2}
\end{aligned}
$$

分母に複素共役を掛けて、分散を表す実部と吸収を表す虚部を分ける。以下の $x$ は実周波数である。

$$
\begin{aligned}
\Re\chi(x) &= \frac{2}{x^{2} + 4} \\
\Im\chi(x) &= \frac{x}{x^{2} + 4} \\
\Im\chi(x)/x &= \frac{1}{x^{2} + 4}
\end{aligned}
$$

最後の比は $x=0$ に連続に延びる。両側で積分し、$\pi$ で割ると静的な実部に一致する。

$$
\begin{aligned}
\frac1\pi\int_\R\frac{\Im\chi(x)}x\,\dd x &= \frac{1}{2} \\
\chi(0) &= \frac{1}{2}
\end{aligned}
$$
:::
:::example[緩和率が 4 のときの吸収と静的応答]
因果応答に対し、通常の Fourier 変換と逆符号の $e^{i\omega t}$ を用いて、正の虚部に解析的な周波数応答を作る。

時間の積分を直接計算する。収束を決める実部は複素周波数に正の虚部を加えるほど大きくなる。

$$
\begin{aligned}
h(t) &= e^{- 4 t} \\
\chi(\omega) &= \frac{1}{- i \omega + 4}
\end{aligned}
$$

分母に複素共役を掛けて、分散を表す実部と吸収を表す虚部を分ける。以下の $x$ は実周波数である。

$$
\begin{aligned}
\Re\chi(x) &= \frac{4}{x^{2} + 16} \\
\Im\chi(x) &= \frac{x}{x^{2} + 16} \\
\Im\chi(x)/x &= \frac{1}{x^{2} + 16}
\end{aligned}
$$

最後の比は $x=0$ に連続に延びる。両側で積分し、$\pi$ で割ると静的な実部に一致する。

$$
\begin{aligned}
\frac1\pi\int_\R\frac{\Im\chi(x)}x\,\dd x &= \frac{1}{4} \\
\chi(0) &= \frac{1}{4}
\end{aligned}
$$
:::
:::example[緩和率が 6 のときの吸収と静的応答]
因果応答に対し、通常の Fourier 変換と逆符号の $e^{i\omega t}$ を用いて、正の虚部に解析的な周波数応答を作る。

時間の積分を直接計算する。収束を決める実部は複素周波数に正の虚部を加えるほど大きくなる。

$$
\begin{aligned}
h(t) &= e^{- 6 t} \\
\chi(\omega) &= \frac{1}{- i \omega + 6}
\end{aligned}
$$

分母に複素共役を掛けて、分散を表す実部と吸収を表す虚部を分ける。以下の $x$ は実周波数である。

$$
\begin{aligned}
\Re\chi(x) &= \frac{6}{x^{2} + 36} \\
\Im\chi(x) &= \frac{x}{x^{2} + 36} \\
\Im\chi(x)/x &= \frac{1}{x^{2} + 36}
\end{aligned}
$$

最後の比は $x=0$ に連続に延びる。両側で積分し、$\pi$ で割ると静的な実部に一致する。

$$
\begin{aligned}
\frac1\pi\int_\R\frac{\Im\chi(x)}x\,\dd x &= \frac{1}{6} \\
\chi(0) &= \frac{1}{6}
\end{aligned}
$$
:::

## 等角写像と 2 次元ポテンシャル問題

### 調和関数と境界を同時に移す

正則写像の導関数がある点で零でなければ、その点で向きと角度を保つ。これを :term[等角写像]{en="conformal map" yomi="とうかくしゃぞう"} と呼ぶ。Cauchy–Riemann 方程式によって微分行列は回転と一様伸縮の積になる。領域全体で境界値問題を移すときは、像への全単射と境界の対応も確かめる。写像 $w=F(z)$ に対し実関数 $u(w)$ の合成は $\Delta_z(u\circ F)=|F'|^2(\Delta_w u)\circ F$ を満たす。連鎖律の二階の項では $F$ の実部・虚部が調和なので一階微分に掛かる項が消え、Hessian を含む項では二方向の伸縮が同じで交差項が消えるためである。

:::example[くさび内の電位]
$0<\arg z<\beta$ の両辺の電位を $0,V_0$ とする。$w=z^{\pi/\beta}$ は上半平面へ写し、正実軸を $0$、負実軸を $V_0$ と対応させる。$V=V_0\arg w/\pi=V_0\arg z/\beta$ は調和かつ境界値を満たす。電場は $-(V_0/(\beta r))\hat{\vb e}_\theta$ である。例えば $\beta=\pi/2$、$V_0=10\,\mathrm V$、$r=0.1\,\mathrm m$ では大きさは $200/\pi\,\mathrm{V\,m^{-1}}$ である。頂点で異なる境界値が接するため、頂点近傍では理想幾何の電場が発散する。
:::

:::example[円板から上半平面へ]
$w=i(1+z)/(1-z)$ では $\Im w=(1-|z|^2)/|1-z|^2>0$ である。逆写像は $z=(w-i)/(w+i)$ なので上半平面への全単射である。円板の調和関数 $\Re z$ を写すと $\Re[(w-i)/(w+i)]$ となり、$w=2i$ では $1/3$ となる。円周の点 $z=1$ は無限遠に移る。
:::

### 上半平面の Dirichlet 核

十分滑らかでコンパクト台を持つ実境界値 $f$ に対して
$$
u(x,y)=\frac1\pi\int_\R \frac{y f(t)}{(x-t)^2+y^2}\dd t,\qquad y>0
$$
とする。積分核は $(t-z)^{-1}$ の虚部である。従って積分の下で微分すると調和である。$y$ を固定して $x-t=ys$ と置けば核の全質量は $1$ であり、$y\downarrow0$ で任意の原点近傍の外の質量が零になる。$f$ の連続性と有界性で積分を近傍と外に分ければ $u(x,y)\to f(x)$ を得る。この式は有界な区分連続の境界値でも連続点で同じ極限を持つ。

例えば $f=\mathbf1_{[-1,1]}$ なら $u(0,y)=2\arctan(1/y)/\pi$ であり、$y=1$ で $1/2$ になる。円板への変換を組み合わせれば、対応する弧にだけ値 $1$ を与えた円板の解も得られる。

有界な調和関数の範囲では、有限個の点を除いて境界値を指定すれば解は一意である。二つの解の差を円板へ移すと、円周の有限個の点を除き境界極限が零となる。上半平面での無限遠も円周上の一点として除いてよい。円板の任意の内点が中心へ移るよう、円板の正則自己同型をさらに合成する。こうして得た調和関数の中心での値は、半径 $r<1$ の円周での平均に等しい。$r\uparrow1$ で有界性により優収束を使うと平均は零になる。任意の内点で差が零なので、一意性が従う。例えばくさびの両辺で零となる $r^{\pi/\beta}\sin(\pi\theta/\beta)$ は無限遠で増大するため、有界性を課さない場合にはこれを足した解も存在する。

### Dirichlet 問題で境界データを移す

$w=f(z)$ が正則で微分が零でなければ、局所的に長さは $|f'(z)|$ 倍、角度は $\arg f'(z)$ だけ回る。実部 $U$ と虚部 $V$ の勾配は等長直交する。従って合成関数の Laplacian を連鎖律で展開すると、混合項が消え、二階項は $|f'|^2$ 倍の Laplacian になる。さらに一階項は $U,V$ が調和であるため消える。したがって調和な解を合成して移せるが、法線微分を与える問題では長さの倍率を含める必要がある。くさびでは $z^{\pi/\alpha}$ で開き角 $\alpha$ を $\pi$ に合わせる。

:::example[くさびの多項式解（1）]
$w=z^m$ の虚部を取り、境界の角度が $\pi/m$ 離れた二辺で零になる解を作る。 一次の模型から始める。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

複素数の冪を直交座標で展開して虚部を取る。

$$
\begin{aligned}
U(x,y)\big|_{m=1} &= y
\end{aligned}
$$

連鎖律では二つの二階項が互いに反対符号となる。

$$
\begin{aligned}
U_{xx}\big|_{m=1} &= 0
\end{aligned}
$$

二階の $y$ 微分は、複素数 $i$ が二回掛かることに対応する。

$$
\begin{aligned}
U_{yy}\big|_{m=1} &= 0
\end{aligned}
$$
:::

:::example[くさびの多項式解（2）]
$w=z^m$ の虚部を取り、境界の角度が $\pi/m$ 離れた二辺で零になる解を作る。 二次の模型を考える。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

複素数の冪を直交座標で展開して虚部を取る。

$$
\begin{aligned}
U(x,y)\big|_{m=2} &= 2 x y
\end{aligned}
$$

連鎖律では二つの二階項が互いに反対符号となる。

$$
\begin{aligned}
U_{xx}\big|_{m=2} &= 0
\end{aligned}
$$

二階の $y$ 微分は、複素数 $i$ が二回掛かることに対応する。

$$
\begin{aligned}
U_{yy}\big|_{m=2} &= 0
\end{aligned}
$$

観測点 $x$=$1$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
U(x,y)\big|_{x=1} &= 2 y
\end{aligned}
$$
:::

:::example[くさびの多項式解（3）]
$w=z^m$ の虚部を取り、境界の角度が $\pi/m$ 離れた二辺で零になる解を作る。 三次の模型を考える。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

複素数の冪を直交座標で展開して虚部を取る。

$$
\begin{aligned}
U(x,y)\big|_{m=3} &= y \left(3 x^{2} - y^{2}\right)
\end{aligned}
$$

連鎖律では二つの二階項が互いに反対符号となる。

$$
\begin{aligned}
U_{xx}\big|_{m=3} &= 6 y
\end{aligned}
$$

二階の $y$ 微分は、複素数 $i$ が二回掛かることに対応する。

$$
\begin{aligned}
U_{yy}\big|_{m=3} &= - 6 y
\end{aligned}
$$

観測点 $x$=$1$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
U(x,y)\big|_{x=1} &= y \left(3 - y^{2}\right)
\end{aligned}
$$
:::

## 演習問題

::::exercise{level=basic}
$\int_{-\infty}^{\infty}(x^2+1)^{-2}\dd x$ を留数で求めよ。
:::solution
上半平面の極は $i$ の二位極である。留数は $(z+i)^{-2}$ の $i$ での微分 $-2/(2i)^3=1/(4i)$ である。円弧は $O(R^{-3})$ で消えるので、値は $2\pi i/(4i)=\pi/2$ である。
:::
::::

::::exercise{level=basic}
$\int_\R x\sin(kx)/(x^2+a^2)\dd x$ を $k>0,a>0$ で求めよ。
:::solution
$ze^{ikz}/(z^2+a^2)$ に約 $1/z$ の減衰があり、Jordan の評価で上の円弧が消える。$ia$ の留数は $e^{-ka}/2$ なので、複素積分は $i\pi e^{-ka}$ である。虚部が求める値 $\pi e^{-ka}$ である。実積分の両側の尾は部分積分で収束する。
:::
::::

::::exercise{level=basic}
$\int_0^{2\pi}(5+3\cos\theta)^{-1}\dd\theta$ を単位円に移し、内部の極を求めよ。
:::solution
分母は $i(3z^2+10z+3)$ で分子は $2$ である。根は $-1/3,-3$ なので内部は $-1/3$ だけである。留数は $2/(8i)$ で、積分値は $\pi/2$ となる。
:::
::::

::::exercise{level=basic}
$w=z^2$ が第 1 象限の角度をどう移すかを述べ、$z=1+i$ の像を求めよ。
:::solution
極形式の偏角は二倍され、$0<\arg z<\pi/2$ は $0<\arg w<\pi$ になる。半径の二乗は狭義単調なので一対一である。$(1+i)^2=2i$ である。
:::
::::

::::exercise{level=advanced}
$\int_0^\infty x^{-1/2}(1+x)^{-2}\dd x$ を式 :eqref[eq-mp-keyhole-beta] から求めよ。
:::solution
$J=\int x^{a-1}/(1+x)^2$、$K=\int x^a/(1+x)^2$ とし、両方とも積分は $0$ から $\infty$ までとする。$J+K=I$ である。端点で消える $x^a/(1+x)$ の導関数を積分すると $aI-K=0$ なので $J=(1-a)I$ を得る。$a=1/2$ を入れると $\pi/2$ である。
:::
::::

::::exercise{level=advanced}
$\operatorname{PV}\int_\R (x^2+1)^{-1}x^{-1}\dd x$ と $1/(x+i0)$ に置き換えた積分を比較せよ。
:::solution
最初の被積分関数は奇関数で、対称な切断ごとに積分は零である。境界値では式 :eqref[eq-mp-plemelj] によりさらに $-i\pi$ が加わる。重みの原点での値が $1$ であるためである。
:::
::::

::::exercise{level=advanced}
式 :eqref[eq-mp-cotangent-expansion] から $\sum_{n=1}^\infty(n^2-1/4)^{-1}$ を求めよ。
:::solution
$z=1/2$ で余接は零である。$1/z=2$、$2z=1$ なので式は $0=2-\sum_{n\ge1}(n^2-1/4)^{-1}$ となり、和は $2$ である。部分分数 $1/(n-1/2)-1/(n+1/2)$ を用いても同じ値になる。
:::
::::

::::exercise{level=advanced}
積分核 $h(t)=t e^{-2t}$（$t\ge0$）の $\chi$ を求め、零周波数で分散関係を確認せよ。
:::solution
部分積分で $\int_0^\infty t e^{-(2-i\omega)t}\dd t=(2-i\omega)^{-2}$ となる。虚部は $4\omega/(4+\omega^2)^2$ である。$\omega=0$ の右辺は $\pi^{-1}\int_\R4/(4+v^2)^2\dd v$ で、$v=2u$ により $\pi^{-1}\cdot \tfrac12\int_\R(1+u^2)^{-2}\dd u=1/4$ となる。左辺の $\Re\chi(0)$ も $1/4$ である。
:::
::::

## 確認テスト

::::quiz{answer="2" tol="1e-9"}
$R(z)=1/(z-2i)$ の極の虚部を答えよ。

:::explanation
分母は $z=2i$ でのみ零になり、分子は非零なので単純極である。座標 $0+2i$ の虚部は $2$ である。
:::
::::

::::quiz
実軸上で $e^{-ikx}$、$k>0$ を含む輪郭積分で減衰する向きを選べ。

- [ ] 上半平面を反時計回りに閉じる
- [x] 下半平面を時計回りに閉じる
- [ ] 上下どちらでも指数は同じ絶対値になる

:::explanation
$z=x+iy$ で $|e^{-ikz}|=e^{ky}$ なので $y<0$ が減衰する。実軸を左から右に進むため、下側の閉曲線の向きは時計回りである。
:::
::::

::::quiz{answer="0.25" tol="1e-9"}
$h(t)=e^{-4t}$、$t\ge0$ の静的応答 $\chi(0)$ を求めよ。

:::explanation
$\omega=0$ では指数変調が消えるため、$\chi(0)=\int_0^\infty e^{-4t}\dd t=1/4=0.25$ である。
:::
::::

::::quiz
正の実軸を切り $0<\arg z<2\pi$ としたとき、$z^{1/2}$ の下岸の値はどれか。

- [ ] $+\sqrt x$
- [x] $-\sqrt x$
- [ ] $i\sqrt x$

:::explanation
下岸の偏角は $2\pi$ に近づき、平方根の位相因子は $e^{i\pi}=-1$ になる。絶対値は $\sqrt x$ である。
:::
::::

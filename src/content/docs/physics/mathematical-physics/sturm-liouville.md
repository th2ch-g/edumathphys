---
title: "Sturm–Liouville 理論と固有関数展開"
description: "境界条件を含む自己共役問題として Sturm–Liouville 方程式を扱い、振動・完全性・変分原理を導く。"
notation:
  - symbol: 'L'
    meaning: "重み付き Sturm–Liouville 作用素"
  - symbol: 'p,q,w'
    meaning: "微分係数、ポテンシャル、重み"
  - symbol: '\lambda_n'
    meaning: "昇順に並べた固有値"
  - symbol: '\mathcal R[u]'
    meaning: "Rayleigh 商"
---

## この章の目標

- 境界形式から許される自己共役境界条件を判定する。
- 固有関数の零点と固有値の順序を結び付ける。
- 逆作用素のコンパクト性から完全性を示す。
- 特異な端点での境界条件と連続スペクトルを見分ける。

[境界値問題](/math/ode/boundary-value/) では端点の条件が解を制約する。本章では同じ条件を課した作用素全体を考え、固有関数が源や初期データの展開に使えることを示す。

## 自己共役微分作用素と境界条件

### 部分積分で残る端点の情報

有限区間 $[a,b]$ で $p\in C^1$、$q,w\in C$ は実数値、$p,w>0$ とする。$L u=w^{-1}[-(pu')'+qu]$ とし、内積 $\inner{u}{v}_w=\int_a^b u\overline v\,w\dd x$ を第 1 引数で線形とする。部分積分により
$$
\inner{Lu}{v}_w-\inner{u}{Lv}_w
=[u\overline{pv'}-pu'\overline v]_a^b
\label{eq-mp-sl-boundary}
$$
となる。$u,pu'$ が絶対連続で $Lu\in L^2(w)$ となる定義域をまず取り、端点で $\alpha_c u(c)+\beta_c p(c)u'(c)=0$ を課す（$c=a,b$、実 $(\alpha_c,\beta_c)\ne(0,0)$）。$\beta=0$ は Dirichlet、$\alpha=0$ は Neumann、それ以外は Robin 条件である。

この領域で $L$ は :term[自己共役作用素]{en="self-adjoint operator" yomi="じこきょうやくさようそ"} である。まず境界条件により二成分 $(u,pu')$ は各端点で同じ一次元空間にあるため境界形式が零となる。逆に随伴の候補 $v$ に対し台が内部にある試験関数を入れると、超関数として $-(pv')'+qv\in L^2$ である。これを一回積分すると $pv'$ は絶対連続な代表を持ち、$p^{-1}$ を掛けて再度積分すると $v$ も絶対連続である。残る式 :eqref[eq-mp-sl-boundary] は端点ごとに同じ一次元空間の全ベクトルと零にならねばならない。許された端データは一点ずつ Hermite 補間で実現できるので、$v$ にも同じ境界条件が必要である。例えば周期条件 $u(a)=u(b),p(a)u'(a)=p(b)u'(b)$ でも二端点を相殺でき、同様の任意性から自己共役である。

### 二つの端を結ぶ自己共役条件

端の値と流束を $\Gamma_0u=(u(a),u(b))^{\mathsf T}$、$\Gamma_1u=(p(a)u'(a),-p(b)u'(b))^{\mathsf T}$ と並べる。このとき境界形式は $\inner{\Gamma_1u}{\Gamma_0v}-\inner{\Gamma_0u}{\Gamma_1v}$ である。行列 $A,B\in\C^{2\times2}$ による条件 $A\Gamma_0u+B\Gamma_1u=0$ が自己共役条件となるための必要十分条件は
$$
 \rank[\,A\ \ B\,]=2,\qquad AB^*=BA^*
$$
である。実際、$M=[\,A\ \ B\,]$、$J=\begin{pmatrix}0&I\\-I&0\end{pmatrix}$ とすると境界形式は $\inner{J(\Gamma_0u,\Gamma_1u)^{\mathsf T}}{(\Gamma_0v,\Gamma_1v)^{\mathsf T}}$ である。$J$ は可逆なので境界形式の直交補の次元は $4-\dim\Ker M$ である。従って $\Ker M$ が自分の直交補と一致するには $\rank M=2$ が必要である。このときその直交補は $\Img(JM^*)$ に等しいので、$\Ker M$ と一致することと $MJM^*=AB^*-BA^*=0$ は同値である。端の値は独立に指定できるから、境界形式での直交補が随伴の端条件を与える。

例えば両端の値と流束を同じ位相 $e^{i\vartheta}$ で結ぶ $u(b)=e^{i\vartheta}u(a)$、$p(b)u'(b)=e^{i\vartheta}p(a)u'(a)$ は、この条件を満たす。実際、両端での境界形式は $|e^{i\vartheta}|^2=1$ により相殺し、許容端データの次元は二である。周期条件は $\vartheta=0$ に当たる。分離条件のように各端の値を一つの直線に制限する場合と違い、ここでは一方の端の二成分が他方の端を決める。

### 混合条件と境界項

正則な有限区間で $wLu=-(p u^\prime)^\prime+q u$ とする。Green の境界形式は端での二成分 $(u,p u^\prime)$ の交代形式である。各端で $(\cos\alpha)u+(\sin\alpha)p u^\prime=0$ を課すと、許される二成分は一次元の直線へ入る。このため任意の二つの許容関数の境界形式は零になる。逆に端の値と導関数の組は三次の補間多項式で独立に指定できる。従って随伴領域でも同じ直線の条件が必要になる。

:::example[Dirichlet と Neumann の組（1）]
$[0,1]$ で左の値を零、右の導関数を零とした $-u^{\prime\prime}$ のモードを作る。

$$
m=1
$$

左端の条件から正弦を選ぶ。

$$
\begin{aligned}
u_m(x)\big|_{m=1} &= \sin{\left(\frac{\pi x}{2} \right)}
\end{aligned}
$$

右端で余弦が零になる半整数を波数にする。

$$
\begin{aligned}
\lambda_m\big|_{m=1} &= \frac{\pi^{2}}{4}
\end{aligned}
$$

この境界の正弦の二乗積分を計算する。

$$
\begin{aligned}
\|u_m\|^2 &= \frac{1}{2}
\end{aligned}
$$

$x=\frac{1}{3}$ で値と導関数を評価する。

$$
\begin{aligned}
u_m(x)\big|_{x=\frac{1}{3}} &= \frac{1}{2} \\
\frac{\partial^1 (u_m(x))}{\partial x^1} &= \frac{\pi \cos{\left(\frac{\pi x}{2} \right)}}{2} \\
\left.\frac{\partial^1 (u_m(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{\sqrt{3} \pi}{4}
\end{aligned}
$$
:::

:::example[Dirichlet と Neumann の組（2）]
$[0,1]$ で左の値を零、右の導関数を零とした $-u^{\prime\prime}$ のモードを作る。値を変えて依存性を調べる。

$$
m=2
$$

左端の条件から正弦を選ぶ。

$$
\begin{aligned}
u_m(x)\big|_{m=2} &= \sin{\left(\frac{3 \pi x}{2} \right)}
\end{aligned}
$$

右端で余弦が零になる半整数を波数にする。

$$
\begin{aligned}
\lambda_m\big|_{m=2} &= \frac{9 \pi^{2}}{4}
\end{aligned}
$$

この境界の正弦の二乗積分を計算する。

$$
\begin{aligned}
\|u_m\|^2 &= \frac{1}{2}
\end{aligned}
$$

$x=\frac{1}{3}$ で値と導関数を評価する。

$$
\begin{aligned}
u_m(x)\big|_{x=\frac{1}{3}} &= 1 \\
\frac{\partial^1 (u_m(x))}{\partial x^1} &= \frac{3 \pi \cos{\left(\frac{3 \pi x}{2} \right)}}{2} \\
\left.\frac{\partial^1 (u_m(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= 0
\end{aligned}
$$
:::

:::example[Dirichlet と Neumann の組（3）]
$[0,1]$ で左の値を零、右の導関数を零とした $-u^{\prime\prime}$ のモードを作る。さらに値を変えて比較する。

$$
m=3
$$

左端の条件から正弦を選ぶ。

$$
\begin{aligned}
u_m(x)\big|_{m=3} &= \sin{\left(\frac{5 \pi x}{2} \right)}
\end{aligned}
$$

右端で余弦が零になる半整数を波数にする。

$$
\begin{aligned}
\lambda_m\big|_{m=3} &= \frac{25 \pi^{2}}{4}
\end{aligned}
$$

この境界の正弦の二乗積分を計算する。

$$
\begin{aligned}
\|u_m\|^2 &= \frac{1}{2}
\end{aligned}
$$

$x=\frac{1}{3}$ で値と導関数を評価する。

$$
\begin{aligned}
u_m(x)\big|_{x=\frac{1}{3}} &= \frac{1}{2} \\
\frac{\partial^1 (u_m(x))}{\partial x^1} &= \frac{5 \pi \cos{\left(\frac{5 \pi x}{2} \right)}}{2} \\
\left.\frac{\partial^1 (u_m(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - \frac{5 \sqrt{3} \pi}{4}
\end{aligned}
$$
:::

:::example[Dirichlet と Neumann の組（4）]
$[0,1]$ で左の値を零、右の導関数を零とした $-u^{\prime\prime}$ のモードを作る。最後にこの値でも計算する。

$$
m=4
$$

左端の条件から正弦を選ぶ。

$$
\begin{aligned}
u_m(x)\big|_{m=4} &= \sin{\left(\frac{7 \pi x}{2} \right)}
\end{aligned}
$$

右端で余弦が零になる半整数を波数にする。

$$
\begin{aligned}
\lambda_m\big|_{m=4} &= \frac{49 \pi^{2}}{4}
\end{aligned}
$$

この境界の正弦の二乗積分を計算する。

$$
\begin{aligned}
\|u_m\|^2 &= \frac{1}{2}
\end{aligned}
$$

$x=\frac{1}{3}$ で値と導関数を評価する。

$$
\begin{aligned}
u_m(x)\big|_{x=\frac{1}{3}} &= - \frac{1}{2} \\
\frac{\partial^1 (u_m(x))}{\partial x^1} &= \frac{7 \pi \cos{\left(\frac{7 \pi x}{2} \right)}}{2} \\
\left.\frac{\partial^1 (u_m(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - \frac{7 \sqrt{3} \pi}{4}
\end{aligned}
$$
:::

### Robin 条件での負の固有値

外向きの導関数に $h u$ を足す境界条件を考える。単位区間で $-u$ を取り、左 Dirichlet、右 $u'(1)+h u(1)=0$ とする。負の固有値 $-\mu^2$ の解は $u(x)=C\sinh(\mu x)$ なので、右条件は $-h=\mu\coth\mu$ になる。$g(\mu)=\mu\coth\mu$ は $\mu>0$ で $1$ から無限大まで狭義単調に増える。実際、導関数の分子は $\sinh\mu\cosh\mu-\mu$ であり、この微分は $2\sinh^2\mu>0$ で、原点での値は零である。従って $h<-1$ なら負の固有値が一つ、$h=-1$ なら線形の解による零固有値、$h>-1$ なら負と零のどちらもない。符号の違いは微分作用素の形式が同じでも、二次形式の端点項を通じて固有値に影響することを示す。
:::example[負の固有値を逆に設計する]
左 Dirichlet、右 Robin で双曲正弦の指数を二に固定する。境界係数を $h=-2\coth 2$ とすれば右条件を満たす。係数を読み、試行関数を微分して方程式へ戻す。

この例の係数を左から二階、一階、零階の順に書く。右辺がある場合は三つの項の和と比較する。

$$
\left( 1, \  0, \  -4\right)
$$

関数とその二つの導関数を並べる。特に一つ下の階数の導関数には、一階項の係数を掛ける。

$$
\begin{aligned}u&=\sinh{\left(2 x \right)}\\u^{\prime}&=2 \cosh{\left(2 x \right)}\\u^{\prime\prime}&=4 \sinh{\left(2 x \right)}\end{aligned}
$$

もとの演算子を作用させると、次のようになる。

$$
\begin{aligned}Lu&=4 \sinh{\left(2 x \right)}+\left(0\right)+\left(- 4 \sinh{\left(2 x \right)}\right)\\&=0\end{aligned}
$$

$x$ が $1$ のところで値を評価する。変換した方程式を使う場合も、初期値や境界値は元の変数で照合する。

$$
\begin{aligned}u&=\sinh{\left(2 \right)}\\ u^\prime&=2 \cosh{\left(2 \right)}\\ Lu&=0\end{aligned}
$$
:::

:::example[指数率 1 の負固有値]
左 Dirichlet と右 Robin の条件を順に照合し、固有関数を規格化する。

零を左端に置くため双曲正弦を選ぶ。$-u^{\prime\prime}$ と $u$ の比が固有値である。

$$
\begin{aligned}
u &= \sinh{\left(x \right)} \\
u(0) &= 0 \\
-u^{\prime\prime}/u &= -1
\end{aligned}
$$

右端に代入する。$h$ は無次元化した Robin 係数である。

$$
\begin{aligned}
u(1) &= \sinh{\left(1 \right)} \\
u^\prime(1) &= \cosh{\left(1 \right)} \\
h &= - \coth{\left(1 \right)} \\
u^\prime(1)+h u(1) &= 0
\end{aligned}
$$

二乗積分の正の平方根で割る。$u$ の符号はこの例では正側に固定する。

$$
\begin{aligned}
\|u\|^2 &= \frac{- 4 e^{2} - 1 + e^{4}}{8 e^{2}} \\
\phi &= \frac{2 \sqrt{2} e \sinh{\left(x \right)}}{\sqrt{- 4 e^{2} - 1 + e^{4}}}
\end{aligned}
$$
:::


:::example[指数率 3 の負固有値]
左 Dirichlet と右 Robin の条件を順に照合し、固有関数を規格化する。

零を左端に置くため双曲正弦を選ぶ。$-u^{\prime\prime}$ と $u$ の比が固有値である。

$$
\begin{aligned}
u &= \sinh{\left(3 x \right)} \\
u(0) &= 0 \\
-u^{\prime\prime}/u &= -9
\end{aligned}
$$

右端に代入する。$h$ は無次元化した Robin 係数である。

$$
\begin{aligned}
u(1) &= \sinh{\left(3 \right)} \\
u^\prime(1) &= 3 \cosh{\left(3 \right)} \\
h &= - 3 \coth{\left(3 \right)} \\
u^\prime(1)+h u(1) &= 0
\end{aligned}
$$

二乗積分の正の平方根で割る。$u$ の符号はこの例では正側に固定する。

$$
\begin{aligned}
\|u\|^2 &= \frac{- 12 e^{6} - 1 + e^{12}}{24 e^{6}} \\
\phi &= \frac{2 \sqrt{6} e^{3} \sinh{\left(3 x \right)}}{\sqrt{- 12 e^{6} - 1 + e^{12}}}
\end{aligned}
$$
:::


:::example[指数率 4 の負固有値]
左 Dirichlet と右 Robin の条件を順に照合し、固有関数を規格化する。

零を左端に置くため双曲正弦を選ぶ。$-u^{\prime\prime}$ と $u$ の比が固有値である。

$$
\begin{aligned}
u &= \sinh{\left(4 x \right)} \\
u(0) &= 0 \\
-u^{\prime\prime}/u &= -16
\end{aligned}
$$

右端に代入する。$h$ は無次元化した Robin 係数である。

$$
\begin{aligned}
u(1) &= \sinh{\left(4 \right)} \\
u^\prime(1) &= 4 \cosh{\left(4 \right)} \\
h &= - 4 \coth{\left(4 \right)} \\
u^\prime(1)+h u(1) &= 0
\end{aligned}
$$

二乗積分の正の平方根で割る。$u$ の符号はこの例では正側に固定する。

$$
\begin{aligned}
\|u\|^2 &= \frac{- 16 e^{8} - 1 + e^{16}}{32 e^{8}} \\
\phi &= \frac{4 \sqrt{2} e^{4} \sinh{\left(4 x \right)}}{\sqrt{- 16 e^{8} - 1 + e^{16}}}
\end{aligned}
$$
:::

## 固有値の実数性と固有関数の直交性

### 内積から固有値の性質を読む

$Lu=\lambda u$、$u\ne0$ なら $\lambda\|u\|_w^2=\inner{Lu}{u}_w=\inner{u}{Lu}_w=\overline\lambda\|u\|_w^2$ なので $\lambda$ は実数である。$Lv=\mu v$ なら $(\lambda-\mu)\inner{u}{v}_w=0$ であり、異なる固有値では直交する。分離境界条件で固有値が単純なのは、左端の許された初期値が一次元で、二階方程式の初期値解が一意だからである。周期条件ではこの議論は使えない。例えば $[0,2\pi]$ の $-u''$ では $\sin x$ と $\cos x$ がともに周期固有値 $1$ を持つ。

### 境界条件による縮退の違い

分離条件では、一方の端で許される初期データが一次元になり、同じ固有値で独立な二解は取れない。連結条件ではこの議論が適用できない。例えば周期境界では $u(1)=u(0),u^\prime(1)=u^\prime(0)$ により両端の寄与が相殺する。正弦と余弦が同じ周期で戻るため、定数以外の自由モードには重複がある。

:::example[周期問題の二方向（1）]
区間 $[0,1]$ の正の固有値について、固有空間の二つの基底を取る。

$$
m=1
$$

偶の成分を作る。

$$
\begin{aligned}
c_m(x)\big|_{m=1} &= \cos{\left(2 \pi x \right)}
\end{aligned}
$$

奇の成分を作る。

$$
\begin{aligned}
s_m(x)\big|_{m=1} &= \sin{\left(2 \pi x \right)}
\end{aligned}
$$

共通の固有値と内積を求める。

$$
\begin{aligned}
\lambda_m\big|_{m=1} &= 4 \pi^{2}
\end{aligned}
$$

直交性は正弦と余弦の積を一周積分して確認できる。

$$
\begin{aligned}
\langle c_m,s_m\rangle &= 0
\end{aligned}
$$

$x=\frac{1}{3}$ で値と導関数を評価する。

$$
\begin{aligned}
c_m(x)\big|_{x=\frac{1}{3}} &= - \frac{1}{2} \\
\frac{\partial^1 (c_m(x))}{\partial x^1} &= - 2 \pi \sin{\left(2 \pi x \right)} \\
\left.\frac{\partial^1 (c_m(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - \sqrt{3} \pi
\end{aligned}
$$

$x=\frac{1}{3}$ で値と導関数を評価する。

$$
\begin{aligned}
s_m(x)\big|_{x=\frac{1}{3}} &= \frac{\sqrt{3}}{2} \\
\frac{\partial^1 (s_m(x))}{\partial x^1} &= 2 \pi \cos{\left(2 \pi x \right)} \\
\left.\frac{\partial^1 (s_m(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - \pi
\end{aligned}
$$
:::

:::example[周期問題の二方向（2）]
区間 $[0,1]$ の正の固有値について、固有空間の二つの基底を取る。値を変えて依存性を調べる。

$$
m=2
$$

偶の成分を作る。

$$
\begin{aligned}
c_m(x)\big|_{m=2} &= \cos{\left(4 \pi x \right)}
\end{aligned}
$$

奇の成分を作る。

$$
\begin{aligned}
s_m(x)\big|_{m=2} &= \sin{\left(4 \pi x \right)}
\end{aligned}
$$

共通の固有値と内積を求める。

$$
\begin{aligned}
\lambda_m\big|_{m=2} &= 16 \pi^{2}
\end{aligned}
$$

直交性は正弦と余弦の積を一周積分して確認できる。

$$
\begin{aligned}
\langle c_m,s_m\rangle &= 0
\end{aligned}
$$

$x=\frac{1}{3}$ で値と導関数を評価する。

$$
\begin{aligned}
c_m(x)\big|_{x=\frac{1}{3}} &= - \frac{1}{2} \\
\frac{\partial^1 (c_m(x))}{\partial x^1} &= - 4 \pi \sin{\left(4 \pi x \right)} \\
\left.\frac{\partial^1 (c_m(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= 2 \sqrt{3} \pi
\end{aligned}
$$

$x=\frac{1}{3}$ で値と導関数を評価する。

$$
\begin{aligned}
s_m(x)\big|_{x=\frac{1}{3}} &= - \frac{\sqrt{3}}{2} \\
\frac{\partial^1 (s_m(x))}{\partial x^1} &= 4 \pi \cos{\left(4 \pi x \right)} \\
\left.\frac{\partial^1 (s_m(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - 2 \pi
\end{aligned}
$$
:::

:::example[周期問題の二方向（3）]
区間 $[0,1]$ の正の固有値について、固有空間の二つの基底を取る。さらに値を変えて比較する。

$$
m=3
$$

偶の成分を作る。

$$
\begin{aligned}
c_m(x)\big|_{m=3} &= \cos{\left(6 \pi x \right)}
\end{aligned}
$$

奇の成分を作る。

$$
\begin{aligned}
s_m(x)\big|_{m=3} &= \sin{\left(6 \pi x \right)}
\end{aligned}
$$

共通の固有値と内積を求める。

$$
\begin{aligned}
\lambda_m\big|_{m=3} &= 36 \pi^{2}
\end{aligned}
$$

直交性は正弦と余弦の積を一周積分して確認できる。

$$
\begin{aligned}
\langle c_m,s_m\rangle &= 0
\end{aligned}
$$

$x=\frac{1}{3}$ で値と導関数を評価する。

$$
\begin{aligned}
c_m(x)\big|_{x=\frac{1}{3}} &= 1 \\
\frac{\partial^1 (c_m(x))}{\partial x^1} &= - 6 \pi \sin{\left(6 \pi x \right)} \\
\left.\frac{\partial^1 (c_m(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= 0
\end{aligned}
$$

$x=\frac{1}{3}$ で値と導関数を評価する。

$$
\begin{aligned}
s_m(x)\big|_{x=\frac{1}{3}} &= 0 \\
\frac{\partial^1 (s_m(x))}{\partial x^1} &= 6 \pi \cos{\left(6 \pi x \right)} \\
\left.\frac{\partial^1 (s_m(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= 6 \pi
\end{aligned}
$$
:::

:::example[周期問題の二方向（4）]
区間 $[0,1]$ の正の固有値について、固有空間の二つの基底を取る。最後にこの値でも計算する。

$$
m=4
$$

偶の成分を作る。

$$
\begin{aligned}
c_m(x)\big|_{m=4} &= \cos{\left(8 \pi x \right)}
\end{aligned}
$$

奇の成分を作る。

$$
\begin{aligned}
s_m(x)\big|_{m=4} &= \sin{\left(8 \pi x \right)}
\end{aligned}
$$

共通の固有値と内積を求める。

$$
\begin{aligned}
\lambda_m\big|_{m=4} &= 64 \pi^{2}
\end{aligned}
$$

直交性は正弦と余弦の積を一周積分して確認できる。

$$
\begin{aligned}
\langle c_m,s_m\rangle &= 0
\end{aligned}
$$

$x=\frac{1}{3}$ で値と導関数を評価する。

$$
\begin{aligned}
c_m(x)\big|_{x=\frac{1}{3}} &= - \frac{1}{2} \\
\frac{\partial^1 (c_m(x))}{\partial x^1} &= - 8 \pi \sin{\left(8 \pi x \right)} \\
\left.\frac{\partial^1 (c_m(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - 4 \sqrt{3} \pi
\end{aligned}
$$

$x=\frac{1}{3}$ で値と導関数を評価する。

$$
\begin{aligned}
s_m(x)\big|_{x=\frac{1}{3}} &= \frac{\sqrt{3}}{2} \\
\frac{\partial^1 (s_m(x))}{\partial x^1} &= 8 \pi \cos{\left(8 \pi x \right)} \\
\left.\frac{\partial^1 (s_m(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - 4 \pi
\end{aligned}
$$
:::

## 振動定理と比較定理

### 位相で節を数える

Dirichlet 条件で左から $u(a)=0,pu'(a)=1$ の解を出す。$(u,pu')=\rho(\sin\theta,\cos\theta)$、$\theta(a)=0$ と置くと $\rho$ は零にならない。この連続な偏角を :term[プリューファー角]{en="Prüfer angle" yomi="ぷりゅーふぁーかく"} という。二成分の一次系から
$$
\theta'=\frac{\cos^2\theta}{p}+(\lambda w-q)\sin^2\theta.
\label{eq-mp-pruefer}
$$
従って $\theta$ が整数倍の $\pi$ を通る方向は常に増加側である。$\partial_\lambda\theta$ は右辺の $\theta$ 微分を係数とする一次方程式を満たし、強制項は $w\sin^2\theta$ である。積分因子で解くと $x>a$ で正になる。$\lambda w-q\le0$ なら $(pu')'=(q-\lambda w)u$ により初め正の $u$ が増加をやめないため、右端の角は $\pi$ 未満である。大きい正 $\lambda$ では $(u,pu'/\sqrt\lambda)$ の角を同様に計算する。角の速度は $\sqrt\lambda[\cos^2\theta/p+(w-q/\lambda)\sin^2\theta]$ で $c\sqrt\lambda$ 以上になる。尺度を戻しても整数倍の $\pi$ の通過回数は同じである。従って各 $n=1,2,\ldots$ に対し $\theta(b)=n\pi$ を満たす固有値が一つあり、内部の零点は $n-1$ 個である。これが Dirichlet 問題の振動定理である。

同じ $p$ で係数 $\lambda w-q$ を点ごとに増やすと、同じ左角からの右角は減少しない。実際、二つの係数を線形補間して補間変数で角方程式を微分すると、上の強制項が係数の差に $\sin^2\theta$ を掛けた形になる。差が正な開区間を通った後では、比較した二つの角の大小は厳密になる。従って同じ初期零点から比べると、下側の係数を持つ解に次の零点があるなら、上側の係数を持つ解は同じ点までに次の零点を持つ。比較の途中で係数差が正の区間を通過していれば、後者の零点は厳密に手前にある。例えば $u''+\lambda u=0$ の零点間隔は $\pi/\sqrt\lambda$ なので、有界な正の係数を上下の定数で挟めば零点間隔も挟める。

### 分離条件での端点位相

一般の分離条件では $(u(a),pu'(a))$ を零でない許容ベクトルに固定し、その角 $\theta_a$ を $[0,\pi)$ に選ぶ。右端条件 $\alpha_bu+\beta_bpu'=0$ を満たす角も $\pi$ の整数倍の間隔で並ぶ。最初の正の許容角を $\beta\in(0,\pi]$ とする。右が Dirichlet の場合は $\beta=\pi$ である。

$\partial_\lambda\theta(b;\lambda)>0$ は左の初期角が $\lambda$ によらないため Dirichlet と同じ積分因子で示せる。さらに $\lambda\to-\infty$ では $\theta(b;\lambda)\to0$ である。実際、任意の小さい $\epsilon>0$ に対し、式 :eqref[eq-mp-pruefer] の右辺は $\theta=\epsilon$ で負となる。初期角が $\epsilon$ 以下なら上からこの障壁を越えられない。初期角が $\epsilon$ より大きければ、閉区間 $[\epsilon,\theta_a]\subset(0,\pi)$ では $\sin^2\theta$ に正の下限があるため、十分負の $\lambda$ で角は $-c_\epsilon|\lambda|$ 以下の速度で減少し、$b$ より前に $\epsilon$ に達する。いずれの場合も $\theta=0$ を下向きには通れない。従って $0<\theta(b;\lambda)\le\epsilon$ である。$\lambda\to+\infty$ では先ほどの尺度付きの角の速度の下限が左端の角に依存せず、通過回数が発散するので $\theta(b;\lambda)\to\infty$ である。

従って $\theta(b;\lambda)=\beta+(n-1)\pi$ を満たす各固有値は一つで、第 $n$ 固有関数の内部零点は $n-1$ 個である。左端の角は $[0,\pi)$ にあり、整数倍の $\pi$ は増加方向にしか横切れない。右端の角がちょうど $n\pi$ の場合は、その右端自体を内部零点から除く。例えば $[0,1]$ で $u(0)=0,u'(1)=0$ なら固有波数は $(n-1/2)\pi$ で、内部零点は $x=\ell/(n-1/2)$（$\ell=1,\ldots,n-1$）である。

### 零点での位相の増加

$u=R\sin\theta,\ pu^\prime=R\cos\theta$ とする。$R^2=u^2+(pu^\prime)^2$ は非零解では零にならない。二成分の微分を代入すると $\theta^\prime=p^{-1}\cos^2\theta+(\lambda w-q)\sin^2\theta$ になる。零点では $\theta$ は $\pi$ の整数倍を正の速度 $1/p$ で横切る。従って振動を調べるには、零点付近だけの位相の向きと端点までの全回転を追えばよい。

:::example[区間を縮めたときの零点（1）]
長さ $2$ の両端 Dirichlet 問題で第 $m$ モードの最初の零点を調べる。

$$
m=1
$$

波数を決める。

$$
\begin{aligned}
\omega_m\big|_{m=1} &= \frac{\pi}{2}
\end{aligned}
$$

固有関数の次の零点は位相が $\pi$ 進んだ点である。$m=1$ ではここは右端に当たる。

$$
\begin{aligned}
x_1\big|_{m=1} &= 2
\end{aligned}
$$

そこでの導関数は負であり、振幅が非零なら単純零点である。

$$
\begin{aligned}
u_m^\prime(x_1)\big|_{m=1} &= - \frac{\pi}{2}
\end{aligned}
$$

内部で同じ手順を繰り返す回数を得る。

$$
\begin{aligned}
N_{\rm interior}\big|_{m=1} &= 0
\end{aligned}
$$
:::

:::example[区間を縮めたときの零点（2）]
長さ $2$ の両端 Dirichlet 問題で第 $m$ モードの最初の零点を調べる。値を変えて依存性を調べる。

$$
m=2
$$

波数を決める。

$$
\begin{aligned}
\omega_m\big|_{m=2} &= \pi
\end{aligned}
$$

固有関数の次の零点は位相が $\pi$ 進んだ点である。$m=1$ ではここは右端に当たる。

$$
\begin{aligned}
x_1\big|_{m=2} &= 1
\end{aligned}
$$

そこでの導関数は負であり、振幅が非零なら単純零点である。

$$
\begin{aligned}
u_m^\prime(x_1)\big|_{m=2} &= - \pi
\end{aligned}
$$

内部で同じ手順を繰り返す回数を得る。

$$
\begin{aligned}
N_{\rm interior}\big|_{m=2} &= 1
\end{aligned}
$$
:::

:::example[区間を縮めたときの零点（3）]
長さ $2$ の両端 Dirichlet 問題で第 $m$ モードの最初の零点を調べる。さらに値を変えて比較する。

$$
m=3
$$

波数を決める。

$$
\begin{aligned}
\omega_m\big|_{m=3} &= \frac{3 \pi}{2}
\end{aligned}
$$

固有関数の次の零点は位相が $\pi$ 進んだ点である。$m=1$ ではここは右端に当たる。

$$
\begin{aligned}
x_1\big|_{m=3} &= \frac{2}{3}
\end{aligned}
$$

そこでの導関数は負であり、振幅が非零なら単純零点である。

$$
\begin{aligned}
u_m^\prime(x_1)\big|_{m=3} &= - \frac{3 \pi}{2}
\end{aligned}
$$

内部で同じ手順を繰り返す回数を得る。

$$
\begin{aligned}
N_{\rm interior}\big|_{m=3} &= 2
\end{aligned}
$$
:::

:::example[区間を縮めたときの零点（4）]
長さ $2$ の両端 Dirichlet 問題で第 $m$ モードの最初の零点を調べる。最後にこの値でも計算する。

$$
m=4
$$

波数を決める。

$$
\begin{aligned}
\omega_m\big|_{m=4} &= 2 \pi
\end{aligned}
$$

固有関数の次の零点は位相が $\pi$ 進んだ点である。$m=1$ ではここは右端に当たる。

$$
\begin{aligned}
x_1\big|_{m=4} &= \frac{1}{2}
\end{aligned}
$$

そこでの導関数は負であり、振幅が非零なら単純零点である。

$$
\begin{aligned}
u_m^\prime(x_1)\big|_{m=4} &= - 2 \pi
\end{aligned}
$$

内部で同じ手順を繰り返す回数を得る。

$$
\begin{aligned}
N_{\rm interior}\big|_{m=4} &= 3
\end{aligned}
$$
:::

### 比較の不等式を積分から示す

:::proposition[節の間での Sturm 比較]{#prop-mp-sl-zero-comparison}
$u,v$ はそれぞれ $-(pu')'=q_1u$、$-(pv')'=q_2v$ の実非零解とする。$p\in C^1$ は正、$q_i$ は連続、$q_2>q_1$ が区間全体で成り立つとする。$u$ の隣接する二零点の間には、$v$ の零点が少なくとも一つある。
:::

:::proof
$u$ の隣接零点を $a,b$ とする。$v$ に内部の零点がなければ、符号を選び直して内部で $u,v>0$ とできる。端では $v(a),v(b)\ge0$ であり、$u'(a)>0$ と $u'(b)<0$ である。後の二つは端まで正の解が近づく向きと初期値の一意性から出る。$W=p(u'v-uv')$ とすると $W'=(q_2-q_1)uv>0$ である。ところが $W(a)=p(a)u'(a)v(a)\ge0$、$W(b)=p(b)u'(b)v(b)\le0$ なので、積分して得る $W(b)-W(a)>0$ に反する。
:::

例えば $u=\sin(2\pi x)$ と $v=\cos(3\pi x)$ では、$u$ の最初の節の区間 $(0,1/2)$ の中で $v$ が $1/6$ に零点を持つ。右端 $1/2$ にも $v$ の零点があるが、比較定理の主張は内部の零点についてであり、端だけを数えて証明を置き換えることはできない。
:::example[Dirichlet の節での微分]
右端も零になる三つ目のモードを取り、固有関数を直接確認する。

固有関数の形は正弦であり、区間の両端での値を調べる。

$$
\begin{aligned}
f(x) &= \sin{\left(\frac{3 \pi x}{2} \right)} \\
f(0) &= 0 \\
f(L) &= 0
\end{aligned}
$$

一階導関数を使い、自然境界条件または内部の節での傾きを照合する。

$$
\begin{aligned}
f^\prime(x) &= \frac{3 \pi \cos{\left(\frac{3 \pi x}{2} \right)}}{2} \\
f^\prime(L) &= - \frac{3 \pi}{2}
\end{aligned}
$$

二階微分の比から固有値を求める。併せて二乗を区間全体で積分して規格化に使う。

$$
\begin{aligned}
\lambda &= \frac{9 \pi^{2}}{4} \\
\|f\|^2 &= 1
\end{aligned}
$$
:::

## 完全性と固有関数展開

### 逆作用素から展開を作る

Dirichlet の場合、$C$ を十分大きくして $q+Cw>0$ とする。$(L+C)u=0$ の斉次境界解は、$\inner{(L+C)u}{u}_w=\int[p|u'|^2+(q+Cw)|u|^2]$ により零しかない。左と右を満たす二つの解 $u_a,u_b$ の加重 Wronskian $W=p(u_a u_b'-u_a'u_b)$ は一定で非零である。$G(x,t)=-u_a(\min(x,t))u_b(\max(x,t))/W$ と置くと $G$ は連続、$pG_x$ の跳びが $-1$ で、$(L+C)^{-1}f=\int G(x,t)f(t)w(t)\dd t$ となる。区間上で連続な核は一様に長方形ごとの定数に近似できる。対応する有限階数作用素との差のノルムは Cauchy–Schwarz によって核の重み付き $L^2$ ノルム以下なので、逆作用素はコンパクトである。

:ref[thm-fa-compact-spectral]{full} を逆作用素に用いる。その逆作用素の零空間は自明なので、全空間の正規直交基底が固有ベクトルで得られる。逆の固有値を $\mu_n$ とすると元の固有値は $\mu_n^{-1}-C$ である。従って $f=\sum_n\inner{f}{u_n}_w u_n$ が $L^2(w)$ で収束し、$(L-\lambda)^{-1}$ の各係数は $\inner{f}{u_n}_w/(\lambda_n-\lambda)$ となる。例えば区間 $[0,1]$ で $p=w=1,q=0$ なら $u_n=\sqrt2\sin(n\pi x)$、$\lambda_n=n^2\pi^2$ である。

### Robin 条件での逆も有限区間ではコンパクトになる

分離 Robin 条件では部分積分に $-[pu'\bar u]_a^b$ が加わる。条件で $pu'$ が端点値の定数倍で表されるので、この項は $|u(a)|^2,|u(b)|^2$ の実線形結合である。端点値は任意の $\epsilon>0$ で
$$
|u(a)|^2+|u(b)|^2\le
\epsilon\int_a^b|u'|^2\dd x+C_\epsilon\int_a^b|u|^2\dd x
$$
と抑えられる。例えば $a$ で $|u(a)|^2=|u(x)|^2-\int_a^x (|u|^2)'$ として $x$ で平均し、$2|u'||u|\le\eta|u'|^2+\eta^{-1}|u|^2$ を使う。$b$ についても同様である。$p$ の正の下限より $\epsilon$ を小さく取れば負になり得る境界項を微分項で吸収できる。残りを $Cw$ の増加で吸収すれば $(L+C)$ の斉次境界解は零だけになる。後は左解・右解で作った連続核とその対称性を使う前と同じ議論で、正規直交基底を得る。

例えば $-u''$ の区間 $[0,1]$ で $u'(0)=hu(0)$、$u'(1)=-hu(1)$、$h\ge0$ なら境界項は $h(|u(0)|^2+|u(1)|^2)\ge0$ である。$h>0$ では $u'=0$ でも非零定数の分子は正になる。$h=0$ では定数が零モードになる。このように Robin の係数の符号は外向き法線と対応させて確認できる。

### スペクトル展開での残差の見積もり

直交射影 $\Pi_N f$ の誤差の二乗は残りの係数の二乗和である。$f$ が微分作用素の領域に入り、$Lf$ も二乗可積分なら、Green の式により $(Lf)_j=\lambda_j f_j$ が成り立つ。従って正の固有値が単調に増える部分では、$\sum_{j>N}|f_j|^2\le\lambda_{N+1}^{-2}\|Lf\|^2$ である。この評価に境界条件が要るのは、移した微分の境界項を消すためである。

:::example[定数の正弦展開（1）]
正規直交化した Dirichlet 基底へ定数を射影する。

$$
m=1
$$

定数の係数は正弦を積分して求まる。

$$
\begin{aligned}
f_m\big|_{m=1} &= \frac{2 \sqrt{2}}{\pi}
\end{aligned}
$$

その係数の二乗が射影によって回収されたノルムに寄与する。

$$
\begin{aligned}
|f_m|^2\big|_{m=1} &= \frac{8}{\pi^{2}}
\end{aligned}
$$

係数を固有値で割ると $-u^{\prime\prime}=1$ の同じ正規基底での解の係数になる。

$$
\begin{aligned}
u_m\big|_{m=1} &= \frac{2 \sqrt{2}}{\pi^{3}}
\end{aligned}
$$
:::

:::example[定数の正弦展開（2）]
正規直交化した Dirichlet 基底へ定数を射影する。値を変えて依存性を調べる。

$$
m=2
$$

定数の係数は正弦を積分して求まる。

$$
\begin{aligned}
f_m\big|_{m=2} &= 0
\end{aligned}
$$

その係数の二乗が射影によって回収されたノルムに寄与する。

$$
\begin{aligned}
|f_m|^2\big|_{m=2} &= 0
\end{aligned}
$$

係数を固有値で割ると $-u^{\prime\prime}=1$ の同じ正規基底での解の係数になる。

$$
\begin{aligned}
u_m\big|_{m=2} &= 0
\end{aligned}
$$
:::

:::example[定数の正弦展開（3）]
正規直交化した Dirichlet 基底へ定数を射影する。さらに値を変えて比較する。

$$
m=3
$$

定数の係数は正弦を積分して求まる。

$$
\begin{aligned}
f_m\big|_{m=3} &= \frac{2 \sqrt{2}}{3 \pi}
\end{aligned}
$$

その係数の二乗が射影によって回収されたノルムに寄与する。

$$
\begin{aligned}
|f_m|^2\big|_{m=3} &= \frac{8}{9 \pi^{2}}
\end{aligned}
$$

係数を固有値で割ると $-u^{\prime\prime}=1$ の同じ正規基底での解の係数になる。

$$
\begin{aligned}
u_m\big|_{m=3} &= \frac{2 \sqrt{2}}{27 \pi^{3}}
\end{aligned}
$$
:::

:::example[定数の正弦展開（4）]
正規直交化した Dirichlet 基底へ定数を射影する。最後にこの値でも計算する。

$$
m=4
$$

定数の係数は正弦を積分して求まる。

$$
\begin{aligned}
f_m\big|_{m=4} &= 0
\end{aligned}
$$

その係数の二乗が射影によって回収されたノルムに寄与する。

$$
\begin{aligned}
|f_m|^2\big|_{m=4} &= 0
\end{aligned}
$$

係数を固有値で割ると $-u^{\prime\prime}=1$ の同じ正規基底での解の係数になる。

$$
\begin{aligned}
u_m\big|_{m=4} &= 0
\end{aligned}
$$
:::

### 周期条件での逆作用素

周期条件では $-[pu'\bar u]_a^b=0$ なので、$q+Cw>0$ なら $(L+C)u=0$ の周期解は零である。任意の $f\in L^2(w)$ に対して Dirichlet の逆による解 $u_D$ をまず作る。次に $(L+C)v=0$ の二次元解空間から $v$ を取り、$u_D+v$ が周期条件を満たすよう二つの係数を選ぶ。周期条件への残差を取る写像は、同次問題に非零の周期解が存在しないため単射であり、両側の次元が二なので全射である。よって補正は一意で、Dirichlet の逆との差は階数二以下である。この差は有界である。実際、Dirichlet の核の端での $p\partial_xG$ は $t$ の有界関数であり、補正の二係数はその端での積分値の線形結合になる。従って周期条件でも逆はコンパクトで、重複固有値の固有空間ごとに正規直交基底を選べば同じ $L^2$ 展開が得られる。

### 逆作用素の固有関数と核の対称性

有限区間での実対称核 $G(x,y)$ が $L^2$ に入るとき、固有値零の空間の基底も含めた正規直交固有基底を $\phi_j$ とする。固有対 $(\mu_j,\phi_j)$ に対して核も $\sum_j\mu_j\phi_j(x)\overline{\phi_j(y)}$ と $L^2$ で一致する。実際、直積の正規直交基底 $\phi_i(x)\overline{\phi_j(y)}$ への係数は $\inner{K\phi_j}{\phi_i}_w$（$K$ は核 $G$ の積分作用素） であり、固有値問題により $\mu_j\delta_{ij}$ になる。この基底が完全であることは、直積の関数に一方ずつ直交展開を適用し、Fubini と二乗ノルムの等式を使えば示せる。したがって核の $L^2$ ノルムの二乗は $\sum|\mu_j|^2$ に一致する。

例えば Dirichlet の $-\partial_x^2$ の逆で $G(x,y)=x(1-y)$ が $0<x<y<1$ に成り立つ。対称側も積分して $\iint G^2=2\int_0^1\int_0^y x^2(1-y)^2\dd x\dd y=1/90$ になる。一方固有値からは $\sum_{j=1}^\infty 1/(\pi j)^4$ である。比較により $\sum j^{-4}=\pi^4/90$ が得られる。この計算は級数の点ごとの収束を仮定せず、逆作用素の核の二乗可積分性で済む。

## Rayleigh 商と固有値の変分特性付け

### 試行関数で固有値を挟む

Dirichlet の :term[レイリー商]{en="Rayleigh quotient" yomi="れいりーしょう"} を $\mathcal R[u]=\int(p|u'|^2+q|u|^2)/\int w|u|^2$ とする。作用素の定義域では分子は $\inner{Lu}{u}_w$ である。展開により $\mathcal R=\sum\lambda_n|c_n|^2/\sum|c_n|^2$ となるため、最小値は $\lambda_1$ で、等号は第一固有関数の倍数である。さらに $j$ 次元部分空間 $S$ の中には $u_1,\ldots,u_{j-1}$ と直交する非零元があるので $\max_S\mathcal R\ge\lambda_j$ である。$S=\spn(u_1,\ldots,u_j)$ では最大が $\lambda_j$ なので
$$
\lambda_j=\min_{\dim S=j}\max_{u\in S\setminus\{0\}}\mathcal R[u].
\label{eq-mp-sl-minimax}
$$
例えば $-u''$ の単位区間で $x(1-x)$ を入れると分子 $1/3$、分母 $1/30$ なので上界 $10$ を得る。一般に $q$ を $q+cw$ に置き換えると全商が定数 $c$ だけ上がり、全固有値が同じ量だけ移る。

### 形式領域での最小最大原理

分離条件のうち Dirichlet 端だけで値を零とした $H^1(a,b)$ の閉部分空間を $V$ とする。Robin または Neumann の端では $pu'(c)=\gamma_cu(c)$ と書ける。このとき分子を
$$
\mathfrak a[u,v]=\int_a^b (p u'\overline{v'}+q u\overline v)\dd x
 +\gamma_a u(a)\overline{v(a)}-\gamma_b u(b)\overline{v(b)}
$$
とする。Dirichlet 端の項は省く。$u\in D(L)$ なら $\mathfrak a[u,v]=\inner{Lu}{v}_w$ は任意の $v\in V$ で成立する。一方 $V$ の関数は $H^1$ でさえあればよく、端の導関数を指定する必要はない。これは試行関数を選ぶ自由度を増やす。

先の端点評価により、十分大きい $C$ に対して $\mathfrak a_C[u,v]=\mathfrak a[u,v]+C\inner{u}{v}_w$ が定めるノルムは $H^1$ のノルムと同値である。固有関数 $u_n$ はこの内積でも直交し、$\mathfrak a_C[u_n,u_n]=\lambda_n+C$ である。もし $v\in V$ が $\mathfrak a_C$ で全ての $u_n$ と直交すれば、$\mathfrak a_C[v,u_n]=(\lambda_n+C)\inner{v}{u_n}_w=0$ だから $L^2$ の完全性から $v=0$ である。従って $u_n/\sqrt{\lambda_n+C}$ は $V$ の正規直交基底であり、
$$
 \mathfrak a_C[u,u]=\sum_{n\ge1}(\lambda_n+C)|\inner{u}{u_n}_w|^2 .
$$
重み付きノルムを引けば Rayleigh 商の展開が形式領域にも延びる。このため式 :eqref[eq-mp-sl-minimax] では部分空間 $S$ を $V$ の任意の $j$ 次元部分空間に取れる。例えば Neumann 問題では、端で導関数が零にならない試行関数でも上界の計算に使える。

停留条件を直接使うなら $u\in V$、$\inner{u}{u}_w=1$ として $\mathfrak a[u,v]=\lambda\inner{u}{v}_w$ を全ての $v\in V$ に課す。内部に台を持つ $v$ を選ぶと分布として $-(pu')'+qu=\lambda wu$ となり、右辺の可積分性から $pu'$ は絶対連続である。次に部分積分して端の値を独立に動かせば、Dirichlet でない端に元の Robin または Neumann 条件が戻る。有限試行空間ではこの等式のテスト関数もその空間に制限されるので、基底 $\psi_i$ に対して $Ac=\lambda Mc$（$A_{ij}=\mathfrak a[\psi_j,\psi_i]$、$M_{ij}=\inner{\psi_j}{\psi_i}_w$）を解くことになる。

### 試行関数の微分と商

Rayleigh 商の分母を $1$ に規格化すると、二次形式の停留条件が固有値式になる。試行関数を有限の部分空間に制限すると、その空間の中での停留条件しか得られない。第一固有値の上界は得られるが、各点で微分方程式が満たされるとは限らない。例えば非規格化の試行関数に定数を掛けても、分子と分母の二乗が相殺して商は変わらない。

:::example[同じ境界に従う多項式の商（1）]
試行関数 $u=x^m(1-x)$ を取り、$[0,1]$ で Dirichlet の商を計算する。

$$
m=1
$$

まず端点で零になる関数を選ぶ。

$$
\begin{aligned}
u\big|_{m=1} &= x \left(1 - x\right)
\end{aligned}
$$

分母では二乗を三つの積分に展開できる。

$$
\begin{aligned}
\int_0^1u^2\dd x\big|_{m=1} &= \frac{1}{30}
\end{aligned}
$$

分子では微分した式を平方してから積分する。

$$
\begin{aligned}
\int_0^1(u^\prime)^2\dd x\big|_{m=1} &= \frac{1}{3}
\end{aligned}
$$

比を取ると高次で右端付近に局在する試行関数ほど商が大きくなる。

$$
\begin{aligned}
R[u]\big|_{m=1} &= 10
\end{aligned}
$$

$x=\frac{1}{3}$ で値と導関数を評価する。

$$
\begin{aligned}
u\big|_{x=\frac{1}{3}} &= \frac{2}{9} \\
\frac{\partial^1 (u)}{\partial x^1} &= 1 - 2 x \\
\left.\frac{\partial^1 (u)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{1}{3}
\end{aligned}
$$
:::

:::example[同じ境界に従う多項式の商（2）]
試行関数 $u=x^m(1-x)$ を取り、$[0,1]$ で Dirichlet の商を計算する。値を変えて依存性を調べる。

$$
m=2
$$

まず端点で零になる関数を選ぶ。

$$
\begin{aligned}
u\big|_{m=2} &= x^{2} \left(1 - x\right)
\end{aligned}
$$

分母では二乗を三つの積分に展開できる。

$$
\begin{aligned}
\int_0^1u^2\dd x\big|_{m=2} &= \frac{1}{105}
\end{aligned}
$$

分子では微分した式を平方してから積分する。

$$
\begin{aligned}
\int_0^1(u^\prime)^2\dd x\big|_{m=2} &= \frac{2}{15}
\end{aligned}
$$

比を取ると高次で右端付近に局在する試行関数ほど商が大きくなる。

$$
\begin{aligned}
R[u]\big|_{m=2} &= 14
\end{aligned}
$$

$x=\frac{1}{3}$ で値と導関数を評価する。

$$
\begin{aligned}
u\big|_{x=\frac{1}{3}} &= \frac{2}{27} \\
\frac{\partial^1 (u)}{\partial x^1} &= x \left(2 - 3 x\right) \\
\left.\frac{\partial^1 (u)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{1}{3}
\end{aligned}
$$
:::

:::example[同じ境界に従う多項式の商（3）]
試行関数 $u=x^m(1-x)$ を取り、$[0,1]$ で Dirichlet の商を計算する。さらに値を変えて比較する。

$$
m=3
$$

まず端点で零になる関数を選ぶ。

$$
\begin{aligned}
u\big|_{m=3} &= x^{3} \left(1 - x\right)
\end{aligned}
$$

分母では二乗を三つの積分に展開できる。

$$
\begin{aligned}
\int_0^1u^2\dd x\big|_{m=3} &= \frac{1}{252}
\end{aligned}
$$

分子では微分した式を平方してから積分する。

$$
\begin{aligned}
\int_0^1(u^\prime)^2\dd x\big|_{m=3} &= \frac{3}{35}
\end{aligned}
$$

比を取ると高次で右端付近に局在する試行関数ほど商が大きくなる。

$$
\begin{aligned}
R[u]\big|_{m=3} &= \frac{108}{5}
\end{aligned}
$$

$x=\frac{1}{3}$ で値と導関数を評価する。

$$
\begin{aligned}
u\big|_{x=\frac{1}{3}} &= \frac{2}{81} \\
\frac{\partial^1 (u)}{\partial x^1} &= x^{2} \left(3 - 4 x\right) \\
\left.\frac{\partial^1 (u)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{5}{27}
\end{aligned}
$$
:::

:::example[同じ境界に従う多項式の商（4）]
試行関数 $u=x^m(1-x)$ を取り、$[0,1]$ で Dirichlet の商を計算する。最後にこの値でも計算する。

$$
m=4
$$

まず端点で零になる関数を選ぶ。

$$
\begin{aligned}
u\big|_{m=4} &= x^{4} \left(1 - x\right)
\end{aligned}
$$

分母では二乗を三つの積分に展開できる。

$$
\begin{aligned}
\int_0^1u^2\dd x\big|_{m=4} &= \frac{1}{495}
\end{aligned}
$$

分子では微分した式を平方してから積分する。

$$
\begin{aligned}
\int_0^1(u^\prime)^2\dd x\big|_{m=4} &= \frac{4}{63}
\end{aligned}
$$

比を取ると高次で右端付近に局在する試行関数ほど商が大きくなる。

$$
\begin{aligned}
R[u]\big|_{m=4} &= \frac{220}{7}
\end{aligned}
$$

$x=\frac{1}{3}$ で値と導関数を評価する。

$$
\begin{aligned}
u\big|_{x=\frac{1}{3}} &= \frac{2}{243} \\
\frac{\partial^1 (u)}{\partial x^1} &= x^{3} \left(4 - 5 x\right) \\
\left.\frac{\partial^1 (u)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{7}{81}
\end{aligned}
$$
:::

:::example[二次元の試行空間の固有値上界]
この例では各試行関数が作用素の領域に入り、第一と第二の Ritz 値をそのまま変分原理で比較できる。

各関数が端点で零であることを確かめ、二つを列の順に固定する。

$$
\begin{aligned}
\phi_{1}(x) &= x \left(1 - x\right) \\
\phi_{2}(x) &= x^{2} \left(1 - x\right)
\end{aligned}
$$

微分して、積を積分すると剛性行列を得る。

$$
\begin{aligned}
\phi_{1}^\prime &= 1 - 2 x \\
\phi_{2}^\prime &= - x^{2} + 2 x \left(1 - x\right) \\
A &= \begin{pmatrix}\frac{1}{3} & \frac{1}{6}\\\frac{1}{6} & \frac{2}{15}\end{pmatrix}
\end{aligned}
$$

関数自体の積を積分して、質量行列を作る。ここで零でない対角外成分も残して計算する。

$$
\begin{aligned}
B &= \begin{pmatrix}\frac{1}{30} & \frac{1}{60}\\\frac{1}{60} & \frac{1}{105}\end{pmatrix} \\
\det B &= \frac{1}{25200}
\end{aligned}
$$

一般化固有方程式を行列式にして根を取る。第一の根はこの試行空間での最小の Rayleigh 商である。

$$
\begin{aligned}
\det(tB-A) &= \frac{\left(t - 42\right) \left(t - 10\right)}{25200} \\
(\lambda_j)_j &= \left( 10, \  42\right)
\end{aligned}
$$
:::


:::example[同じ対称性に限定した比較]
第一固有値の比較を空間の包含から読む。一方、偶対称の部分空間だけを使うと奇対称の固有モードは近似されない。

各関数が端点で零であることを確かめ、二つを列の順に固定する。

$$
\begin{aligned}
\phi_{1}(x) &= x \left(1 - x\right) \\
\phi_{2}(x) &= x^{2} \left(1 - x\right)^{2}
\end{aligned}
$$

微分して、積を積分すると剛性行列を得る。

$$
\begin{aligned}
\phi_{1}^\prime &= 1 - 2 x \\
\phi_{2}^\prime &= x^{2} \left(2 x - 2\right) + 2 x \left(1 - x\right)^{2} \\
A &= \begin{pmatrix}\frac{1}{3} & \frac{1}{15}\\\frac{1}{15} & \frac{2}{105}\end{pmatrix}
\end{aligned}
$$

関数自体の積を積分して、質量行列を作る。ここで零でない対角外成分も残して計算する。

$$
\begin{aligned}
B &= \begin{pmatrix}\frac{1}{30} & \frac{1}{140}\\\frac{1}{140} & \frac{1}{630}\end{pmatrix} \\
\det B &= \frac{1}{529200}
\end{aligned}
$$

一般化固有方程式を行列式にして根を取る。第一の根はこの試行空間での最小の Rayleigh 商である。

$$
\begin{aligned}
\det(tB-A) &= \frac{t^{2} - 112 t + 1008}{529200} \\
(\lambda_j)_j &= \left( 56 - 4 \sqrt{133}, \  4 \sqrt{133} + 56\right)
\end{aligned}
$$
:::

## 特異 Sturm–Liouville 問題と連続スペクトル

### 端点の可積分性と無限区間

開区間の内部で $p,w>0$、$1/p,q,w$ が局所可積分であるとする。端点が有限で $1/p,q,w$ がその近くで可積分なら、その端点を正則と呼ぶ。無限端点、またはこれらの可積分性が破れる端点は特異端点である。本章の最初の連続係数の仮定は、両端が正則となる十分条件である。標準の Weyl の分類では、非実のスペクトル変数で端点近くの全解が $L^2(w)$ に入るなら極限円、一つの独立解だけが入るなら極限点である。分離条件で自己共役化する場合は、極限円端点ごとに一つの境界条件を課し、極限点端点には課さない。極限円が両端にあるときは両端を結ぶ自己共役条件も選べる。これを非実パラメータに依存しない分類として示すには特異端点の随伴領域と境界形式の極限を扱う必要があり、証明は :cite[coddington-levinson] を参照する。

Bessel の動径問題では $p=w=r,q=\nu^2/r$ であり、零のスペクトル変数での原点の二解は $r^\nu,r^{-\nu}$（$\nu>0$）、$\nu=0$ では $1,\log r$ である。重み $r\dd r$ では $\nu\ge1$ だと減少冪が可積分でなく、$0\le\nu<1$ では両方とも可積分である。分類は、極限円なら実のスペクトル変数でも全解が可積分であるという Weyl の結果を含む。ここでは非実の変数でも原点の指数が同じで、階数低下により特異解の先頭は $r^{-\nu}$ または $\log r$ になるため、直接同じ結論が出る。従って $0\le\nu<1$ では円板の有界な枝を選ぶ境界条件が別に必要である。この正則側と外端 Dirichlet を選んだ場合、零固有値を避けて逆核を直接作れる。$\nu>0$ の核は $r_<^\nu(r_>^{-\nu}-a^{-2\nu}r_>^\nu)/(2\nu)$、$\nu=0$ では $\log(a/r_>)$ である。二解の接続と $pG_x$ の跳びを調べれば得られ、前者は有界、後者の二乗も $r\dd r\,t\dd t$ で可積分である。核を $L^2$ 近似して同じコンパクト性を得るため、円板の Fourier–Bessel 基底は完全である。

半直線 $x>0$ の Dirichlet $-u''$ では、奇拡張して Fourier 変換すると正弦変換 $\sqrt{2/\pi}\int_0^\infty u(x)\sin(kx)\dd x$ になる。Plancherel からこれは $L^2(0,\infty)$ のユニタリ変換で、$-u''$ は $k^2$ の乗算に移る。$k^2=\lambda$ の集合は測度零なので非零の $L^2$ 固有ベクトルはない。一方、任意の $\lambda\ge0$ の近くの有限幅の帯域に台を持つ単位ベクトル $\widehat u$ を取れば、幅を縮めて $(k^2-\lambda)\widehat u$ のノルムを零に近づけられる。従って $[0,\infty)$ はスペクトルであり、負の $\lambda$ は有界な逆乗算を持つ。この例では連続スペクトルが現れる。

### Liouville の正規形

正の $p,w\in C^2$ について、変数 $s$ と振幅 $\rho$ を
$$
 \frac{\dd s}{\dd x}=\sqrt{w/p},\qquad
 \rho=(pw)^{1/4},\qquad v(s)=\rho(x)u(x)
$$
とする。すると $pu'=\rho v_s-\rho_s v$ であり、もう一度微分して $(pu')'=(\dd s/\dd x)(\rho v_{ss}-\rho_{ss}v)$ となる。従って $Lu=\lambda u$ は
$$
 -v_{ss}+\left(\frac q w+\frac{\rho_{ss}}{\rho}\right)v=\lambda v
$$
に移る。測度は $\rho^2\dd s=w\dd x$ なので、同時に $\int|v|^2\dd s=\int |u|^2w\dd x$ が成り立つ。正則端の分離条件は、$u=v/\rho$ と $pu'=\rho v_s-\rho_s v$ を代入して移せる。特異端では端点値が存在するという仮定を置かず、先に二乗可積分な枝を調べる。

例えば Bessel のモデルでは $p=w=r$、$q=\nu^2/r$ なので $s=r$、$\rho=\sqrt r$ である。変換後のポテンシャルは $(\nu^2-1/4)/r^2$ になる。原点での二つの枝は $r^{1/2+\nu}$ と $r^{1/2-\nu}$ で、$\nu=0$ では二つ目は $\sqrt r\log r$ になる。こうして $L^2(\dd r)$ で判定しても、元の重み $r\dd r$ で判定したときと同じ境界条件の数が得られる。

### 無限区間の近似固有関数

正の $k$ を取り、$\eta\in C_c^\infty(1,2)$ を選ぶ。$u_R(x)=R^{-1/2}\eta(x/R)e^{ikx}$ とおけば半直線の端点で零になり、ノルムは $R$ に依存しない。二回微分して $k^2u_R$ を引くと、$\eta^\prime$ の項は $O(R^{-1})$、$\eta^{\prime\prime}$ の項は $O(R^{-2})$ である。もし $-D^2-k^2$ に有界な逆があればノルム一定の列を零へ送ることになり矛盾する。このため $k^2$ はスペクトルに属する。

:::example[局在した波束の残差（1）]
正弦変換で構成する連続スペクトルの模型として、無限遠で切り取った指数の残差を調べる。

$$
m=1
$$

正の軸上での波数を指定する。

$$
\begin{aligned}
k\big|_{m=1} &= 1
\end{aligned}
$$

カットオフの一階微分の項の係数は $2k/R$ である。

$$
\begin{aligned}
|-2ik/R|\big|_{m=1} &= \frac{2}{u}
\end{aligned}
$$

二階微分の項はさらに $1/R$ 小さい。

$$
\begin{aligned}
|-1/R^2| &= \frac{1}{u^{2}}
\end{aligned}
$$

非零のカットオフを定数倍し、$\|\eta'\|_2,\|\eta''\|_2\le1$ となるように取る。このときの上界を求める。

$$
\begin{aligned}
\|(L-k^2)u_R\|\ \text{の上界}\big|_{m=1} &= \frac{2 u + 1}{u^{2}}
\end{aligned}
$$

$u=10$ で値と導関数を評価する。

$$
\begin{aligned}
|-2ik/R|\big|_{u=10} &= \frac{1}{5} \\
\frac{\partial^1 (|-2ik/R|)}{\partial u^1} &= - \frac{2}{u^{2}} \\
\left.\frac{\partial^1 (|-2ik/R|)}{\partial u^1}\right|_{u=10} &= - \frac{1}{50} \\
\frac{\partial^2 (|-2ik/R|)}{\partial u^2} &= \frac{4}{u^{3}} \\
\left.\frac{\partial^2 (|-2ik/R|)}{\partial u^2}\right|_{u=10} &= \frac{1}{250}
\end{aligned}
$$

$u=10$ で値と導関数を評価する。

$$
\begin{aligned}
|-1/R^2|\big|_{u=10} &= \frac{1}{100} \\
\frac{\partial^1 (|-1/R^2|)}{\partial u^1} &= - \frac{2}{u^{3}} \\
\left.\frac{\partial^1 (|-1/R^2|)}{\partial u^1}\right|_{u=10} &= - \frac{1}{500} \\
\frac{\partial^2 (|-1/R^2|)}{\partial u^2} &= \frac{6}{u^{4}} \\
\left.\frac{\partial^2 (|-1/R^2|)}{\partial u^2}\right|_{u=10} &= \frac{3}{5000}
\end{aligned}
$$

$u=10$ で値と導関数を評価する。

$$
\begin{aligned}
\|(L-k^2)u_R\|\ \text{の上界}\big|_{u=10} &= \frac{21}{100} \\
\frac{\partial^1 (\|(L-k^2)u_R\|\ \text{の上界})}{\partial u^1} &= \frac{2 \left(- u - 1\right)}{u^{3}} \\
\left.\frac{\partial^1 (\|(L-k^2)u_R\|\ \text{の上界})}{\partial u^1}\right|_{u=10} &= - \frac{11}{500} \\
\frac{\partial^2 (\|(L-k^2)u_R\|\ \text{の上界})}{\partial u^2} &= \frac{2 \left(2 u + 3\right)}{u^{4}} \\
\left.\frac{\partial^2 (\|(L-k^2)u_R\|\ \text{の上界})}{\partial u^2}\right|_{u=10} &= \frac{23}{5000}
\end{aligned}
$$
:::

:::example[局在した波束の残差（2）]
正弦変換で構成する連続スペクトルの模型として、無限遠で切り取った指数の残差を調べる。値を変えて依存性を調べる。

$$
m=2
$$

正の軸上での波数を指定する。

$$
\begin{aligned}
k\big|_{m=2} &= 2
\end{aligned}
$$

カットオフの一階微分の項の係数は $2k/R$ である。

$$
\begin{aligned}
|-2ik/R|\big|_{m=2} &= \frac{4}{u}
\end{aligned}
$$

二階微分の項はさらに $1/R$ 小さい。

$$
\begin{aligned}
|-1/R^2| &= \frac{1}{u^{2}}
\end{aligned}
$$

非零のカットオフを定数倍し、$\|\eta'\|_2,\|\eta''\|_2\le1$ となるように取る。このときの上界を求める。

$$
\begin{aligned}
\|(L-k^2)u_R\|\ \text{の上界}\big|_{m=2} &= \frac{4 u + 1}{u^{2}}
\end{aligned}
$$

$u=10$ で値と導関数を評価する。

$$
\begin{aligned}
|-2ik/R|\big|_{u=10} &= \frac{2}{5} \\
\frac{\partial^1 (|-2ik/R|)}{\partial u^1} &= - \frac{4}{u^{2}} \\
\left.\frac{\partial^1 (|-2ik/R|)}{\partial u^1}\right|_{u=10} &= - \frac{1}{25} \\
\frac{\partial^2 (|-2ik/R|)}{\partial u^2} &= \frac{8}{u^{3}} \\
\left.\frac{\partial^2 (|-2ik/R|)}{\partial u^2}\right|_{u=10} &= \frac{1}{125}
\end{aligned}
$$

$u=10$ で値と導関数を評価する。

$$
\begin{aligned}
|-1/R^2|\big|_{u=10} &= \frac{1}{100} \\
\frac{\partial^1 (|-1/R^2|)}{\partial u^1} &= - \frac{2}{u^{3}} \\
\left.\frac{\partial^1 (|-1/R^2|)}{\partial u^1}\right|_{u=10} &= - \frac{1}{500} \\
\frac{\partial^2 (|-1/R^2|)}{\partial u^2} &= \frac{6}{u^{4}} \\
\left.\frac{\partial^2 (|-1/R^2|)}{\partial u^2}\right|_{u=10} &= \frac{3}{5000}
\end{aligned}
$$

$u=10$ で値と導関数を評価する。

$$
\begin{aligned}
\|(L-k^2)u_R\|\ \text{の上界}\big|_{u=10} &= \frac{41}{100} \\
\frac{\partial^1 (\|(L-k^2)u_R\|\ \text{の上界})}{\partial u^1} &= \frac{2 \left(- 2 u - 1\right)}{u^{3}} \\
\left.\frac{\partial^1 (\|(L-k^2)u_R\|\ \text{の上界})}{\partial u^1}\right|_{u=10} &= - \frac{21}{500} \\
\frac{\partial^2 (\|(L-k^2)u_R\|\ \text{の上界})}{\partial u^2} &= \frac{2 \left(4 u + 3\right)}{u^{4}} \\
\left.\frac{\partial^2 (\|(L-k^2)u_R\|\ \text{の上界})}{\partial u^2}\right|_{u=10} &= \frac{43}{5000}
\end{aligned}
$$
:::

:::example[局在した波束の残差（3）]
正弦変換で構成する連続スペクトルの模型として、無限遠で切り取った指数の残差を調べる。さらに値を変えて比較する。

$$
m=3
$$

正の軸上での波数を指定する。

$$
\begin{aligned}
k\big|_{m=3} &= 3
\end{aligned}
$$

カットオフの一階微分の項の係数は $2k/R$ である。

$$
\begin{aligned}
|-2ik/R|\big|_{m=3} &= \frac{6}{u}
\end{aligned}
$$

二階微分の項はさらに $1/R$ 小さい。

$$
\begin{aligned}
|-1/R^2| &= \frac{1}{u^{2}}
\end{aligned}
$$

非零のカットオフを定数倍し、$\|\eta'\|_2,\|\eta''\|_2\le1$ となるように取る。このときの上界を求める。

$$
\begin{aligned}
\|(L-k^2)u_R\|\ \text{の上界}\big|_{m=3} &= \frac{6 u + 1}{u^{2}}
\end{aligned}
$$

$u=10$ で値と導関数を評価する。

$$
\begin{aligned}
|-2ik/R|\big|_{u=10} &= \frac{3}{5} \\
\frac{\partial^1 (|-2ik/R|)}{\partial u^1} &= - \frac{6}{u^{2}} \\
\left.\frac{\partial^1 (|-2ik/R|)}{\partial u^1}\right|_{u=10} &= - \frac{3}{50} \\
\frac{\partial^2 (|-2ik/R|)}{\partial u^2} &= \frac{12}{u^{3}} \\
\left.\frac{\partial^2 (|-2ik/R|)}{\partial u^2}\right|_{u=10} &= \frac{3}{250}
\end{aligned}
$$

$u=10$ で値と導関数を評価する。

$$
\begin{aligned}
|-1/R^2|\big|_{u=10} &= \frac{1}{100} \\
\frac{\partial^1 (|-1/R^2|)}{\partial u^1} &= - \frac{2}{u^{3}} \\
\left.\frac{\partial^1 (|-1/R^2|)}{\partial u^1}\right|_{u=10} &= - \frac{1}{500} \\
\frac{\partial^2 (|-1/R^2|)}{\partial u^2} &= \frac{6}{u^{4}} \\
\left.\frac{\partial^2 (|-1/R^2|)}{\partial u^2}\right|_{u=10} &= \frac{3}{5000}
\end{aligned}
$$

$u=10$ で値と導関数を評価する。

$$
\begin{aligned}
\|(L-k^2)u_R\|\ \text{の上界}\big|_{u=10} &= \frac{61}{100} \\
\frac{\partial^1 (\|(L-k^2)u_R\|\ \text{の上界})}{\partial u^1} &= \frac{2 \left(- 3 u - 1\right)}{u^{3}} \\
\left.\frac{\partial^1 (\|(L-k^2)u_R\|\ \text{の上界})}{\partial u^1}\right|_{u=10} &= - \frac{31}{500} \\
\frac{\partial^2 (\|(L-k^2)u_R\|\ \text{の上界})}{\partial u^2} &= \frac{6 \left(2 u + 1\right)}{u^{4}} \\
\left.\frac{\partial^2 (\|(L-k^2)u_R\|\ \text{の上界})}{\partial u^2}\right|_{u=10} &= \frac{63}{5000}
\end{aligned}
$$
:::

:::example[局在した波束の残差（4）]
正弦変換で構成する連続スペクトルの模型として、無限遠で切り取った指数の残差を調べる。最後にこの値でも計算する。

$$
m=4
$$

正の軸上での波数を指定する。

$$
\begin{aligned}
k\big|_{m=4} &= 4
\end{aligned}
$$

カットオフの一階微分の項の係数は $2k/R$ である。

$$
\begin{aligned}
|-2ik/R|\big|_{m=4} &= \frac{8}{u}
\end{aligned}
$$

二階微分の項はさらに $1/R$ 小さい。

$$
\begin{aligned}
|-1/R^2| &= \frac{1}{u^{2}}
\end{aligned}
$$

非零のカットオフを定数倍し、$\|\eta'\|_2,\|\eta''\|_2\le1$ となるように取る。このときの上界を求める。

$$
\begin{aligned}
\|(L-k^2)u_R\|\ \text{の上界}\big|_{m=4} &= \frac{8 u + 1}{u^{2}}
\end{aligned}
$$

$u=10$ で値と導関数を評価する。

$$
\begin{aligned}
|-2ik/R|\big|_{u=10} &= \frac{4}{5} \\
\frac{\partial^1 (|-2ik/R|)}{\partial u^1} &= - \frac{8}{u^{2}} \\
\left.\frac{\partial^1 (|-2ik/R|)}{\partial u^1}\right|_{u=10} &= - \frac{2}{25} \\
\frac{\partial^2 (|-2ik/R|)}{\partial u^2} &= \frac{16}{u^{3}} \\
\left.\frac{\partial^2 (|-2ik/R|)}{\partial u^2}\right|_{u=10} &= \frac{2}{125}
\end{aligned}
$$

$u=10$ で値と導関数を評価する。

$$
\begin{aligned}
|-1/R^2|\big|_{u=10} &= \frac{1}{100} \\
\frac{\partial^1 (|-1/R^2|)}{\partial u^1} &= - \frac{2}{u^{3}} \\
\left.\frac{\partial^1 (|-1/R^2|)}{\partial u^1}\right|_{u=10} &= - \frac{1}{500} \\
\frac{\partial^2 (|-1/R^2|)}{\partial u^2} &= \frac{6}{u^{4}} \\
\left.\frac{\partial^2 (|-1/R^2|)}{\partial u^2}\right|_{u=10} &= \frac{3}{5000}
\end{aligned}
$$

$u=10$ で値と導関数を評価する。

$$
\begin{aligned}
\|(L-k^2)u_R\|\ \text{の上界}\big|_{u=10} &= \frac{81}{100} \\
\frac{\partial^1 (\|(L-k^2)u_R\|\ \text{の上界})}{\partial u^1} &= \frac{2 \left(- 4 u - 1\right)}{u^{3}} \\
\left.\frac{\partial^1 (\|(L-k^2)u_R\|\ \text{の上界})}{\partial u^1}\right|_{u=10} &= - \frac{41}{500} \\
\frac{\partial^2 (\|(L-k^2)u_R\|\ \text{の上界})}{\partial u^2} &= \frac{2 \left(8 u + 3\right)}{u^{4}} \\
\left.\frac{\partial^2 (\|(L-k^2)u_R\|\ \text{の上界})}{\partial u^2}\right|_{u=10} &= \frac{83}{5000}
\end{aligned}
$$
:::

### 半直線の Robin 条件にできる束縛状態

自由作用素 $-u''$ を $x>0$ に取る。端の条件を $u'(0)=h u(0)$ とすると、負の固有値 $-\kappa^2$ に対する二つの枝のうち $e^{+\kappa x}$ は二乗可積分でない。許される枝 $e^{-\kappa x}$ を境界条件へ入れれば $-\kappa=h$ となる。従って $h<0$ のときだけ負固有値 $-h^2$ があり、規格化した関数は $\sqrt{-2h}e^{hx}$ である。$h=0$ で得られる定数は半直線上で二乗可積分でない。このため有限区間の Neumann 零モードを半直線へそのまま延長できない。

例えば $h=-2$ では固有値は $-4$、規格化した解は $2e^{-2x}$ となる。境界で値が $2$、導関数が $-4$ なので右辺の $h u(0)$ と一致する。この離散固有値が存在する場合でも、遠方で切り取った正の波数の近似固有関数は同じ境界条件を満たす。従って境界が束縛状態を加えても、自由作用素の正の連続帯は消えない。

### 半直線の連続帯を Robin 条件に合わせる

$h$ は実数とする。$A=\partial_x-h$ を $D(A)=H^1(0,\infty)$ 上で定める。部分積分により $A^*=-\partial_x-h$ で、その定義域は $H_0^1(0,\infty)$ である。内部の試験関数から随伴の導関数が $L^2$ であることが分かり、残る原点の境界項が $v(0)=0$ を要求する。従って $A^*A=-\partial_x^2+h^2$ は $u'(0)=hu(0)$ を満たす $H^2$ 上に、$AA^*=-\partial_x^2+h^2$ は Dirichlet 条件を満たす $H^2$ 上に作用する。

両作用素の自己共役性は、閉作用素のグラフ内積からも確かめられる。$A$ のグラフノルムは $H^1$ のノルムと同値なので $A$ は閉じており、$A^*$ も閉じている。$B=A$ または $A^*$ とすると、任意の $f\in L^2$ に対し Riesz の表現定理から $u\in D(B)$ が存在し、$\inner{u}{v}+\inner{Bu}{Bv}=\inner{f}{v}$ が全ての $v\in D(B)$ で成り立つ。従って $Bu\in D(B^*)$、$(I+B^*B)u=f$ である。$T=B^*B$ は対称で、$I+T$ が全射だから自己共役である。実際、$v\in D(T^*)$ に対し $(I+T)u=(I+T^*)v$ を解くと $v-u$ は $\Ker(I+T^*)=\Img(I+T)^\perp=\{0\}$ に属する。

正弦変換で $AA^*$ を $k^2+h^2$ の乗算に直す。正規化した連続モード $\sqrt{2/\pi}\sin(kx)$ に $-A^*/\sqrt{k^2+h^2}$ を作用させると
$$
\varphi_k(x)=\sqrt{\frac2\pi}
  \frac{k\cos(kx)+h\sin(kx)}{\sqrt{k^2+h^2}},\qquad k>0
$$
を得る。分子の原点での値は $k$、導関数は $hk$ なので、このモードは確かに Robin 条件に合う。等長性は次のように示せる。正弦変換で $k$ の台が $0$ と $\infty$ から離れたデータ $g$ に対し、$v=(AA^*)^{-1/2}g$ と置くと $\|A^*v\|^2=\inner{AA^*v}{v}=\|g\|^2$ である。この集合が稠密なので $-A^*(AA^*)^{-1/2}$ は等長作用素に延びる。その値域の直交補は $Av=0$ の $L^2$ 解である。実際、使用した $v$ は正弦変換の帯域を拡張すれば $H_0^1$ のノルムで稠密であり、全ての $A^*v$ と直交すれば随伴の定義からこの結論を得る。

よって $h<0$ のときだけ束縛状態を直交補に加えればよく、残りは $\varphi_k$ による連続展開になる。微分作用素の積 $A^*(AA^*)=(A^*A)A^*$ は帯域を絞った稠密集合で成り立つため、この等長変換は Robin 作用素の連続部分を $k^2$ の乗算へ移す。従って連続帯は $[0,\infty)$ であり、それ以外のスペクトルは $h<0$ の単純固有値 $-h^2$ だけである。$h=0$ のときは $\varphi_k=\sqrt{2/\pi}\cos(kx)$ となる。

## 演習問題

::::exercise{level=basic}
$[0,2]$ の Dirichlet $-u''$ の第一固有値を求めよ。
:::solution
左を満たす解は $\sin(kx)$、右で $2k=\pi$ が最初なので $\lambda_1=\pi^2/4$ である。零以下は部分積分の正性と境界条件から固有値にならない。
:::
::::

::::exercise{level=basic}
$[0,1]$ の Neumann $-u''$ の零固有空間を求めよ。
:::solution
$u''=0$ より $u=A+Bx$、端点の導関数が零なので $B=0$ となる。固有空間は定数全体である。
:::
::::

::::exercise{level=basic}
Dirichlet で $\int_0^1\sin(\pi x)\sin(2\pi x)\dd x$ を計算せよ。
:::solution
積を $\tfrac12[\cos(\pi x)-\cos(3\pi x)]$ と書くと両項の積分が零になる。固有値が異なることからも同じ結論である。
:::
::::

::::exercise{level=basic}
$p=w=1,q=3$ の $[0,1]$ Dirichlet の第二固有値を求めよ。
:::solution
定数ポテンシャルは全固有値に $3$ を足すので $4\pi^2+3$ である。
:::
::::

::::exercise{level=advanced}
$p=w=1,q=0$、$[0,1]$ で $u=x(1-x)$ の商を計算せよ。
:::solution
分子は $\int_0^1(1-2x)^2=1/3$、分母は $\int_0^1(x^2-2x^3+x^4)=1/30$ である。従って $10$ で、第一固有値 $\pi^2$ の上界である。
:::
::::

::::exercise{level=advanced}
$r^{-3/4}$ と $r^{-5/4}$ のどちらが原点近くで $L^2(r\dd r)$ に入るか。
:::solution
二乗し重みを掛けるとそれぞれ $r^{-1/2}$ と $r^{-3/2}$ である。原点での冪積分は指数 $>-1$ の場合だけ有限なので前者のみ入る。
:::
::::

::::exercise{level=advanced}
$q_2\ge q_1$ で $p,w$ が同じ二つの Dirichlet 問題の固有値を比較せよ。
:::solution
任意の同じ試行関数に対して分子の差は $\int(q_2-q_1)|u|^2\ge0$ である。各次元の部分空間で最大を取ってから最小を取っても順序が保たれ、$\lambda_j(q_2)\ge\lambda_j(q_1)$ となる。
:::
::::

::::exercise{level=advanced}
Dirichlet の第 4 固有関数の内部零点数を求め、その零点が単純であることを示せ。
:::solution
右端の位相が $4\pi$、内部で $\pi,2\pi,3\pi$ を一度ずつ通るので零点は三つである。もし $u$ と $u'$ がともに零なら初期値一意性で自明解になるため、各零点は単純である。
:::
::::

## 確認テスト

::::quiz{answer="3" tol="1e-9"}
第 4 Dirichlet 固有関数の内部零点数を答えよ。

:::explanation
左端の角 $0$ と右端の $4\pi$ の間に整数倍の $\pi$ は三つある。各点で増加方向に通るため、内部零点は三つである。
:::
::::

::::quiz
一般の直交性に使う内積の重みはどれか。

- [ ] $p$
- [ ] $q$
- [x] $w$

:::explanation
方程式は $-(pu')'+qu=\lambda wu$ であり、二つの解を交差して引くと $(\lambda-\mu)\int wu\overline v$ が残る。
:::
::::

::::quiz{answer="0" tol="1e-9"}
Neumann の $-u''$ で定数関数の固有値を答えよ。

:::explanation
導関数と二階導関数がすべて零である。非零定数なので固有値は $0$ である。
:::
::::

::::quiz
半直線の自由 Dirichlet 作用素のスペクトルについて正しいものはどれか。

- [ ] 正弦は正規化できる $L^2$ 固有関数である
- [x] 正の波数は連続スペクトルを与える
- [ ] 正の波数は必ず離散化する

:::explanation
正弦変換は連続変数 $k$ を与え、作用素は $k^2$ 倍となる。個々の波数の正弦は半直線全体では二乗可積分でない。
:::
::::

---
title: "Green 関数"
description: "点源の応答から、境界、放射条件、時間方向を指定した Green 関数を構成する。"
notation:
  - symbol: 'G(x,y)'
    meaning: "点 y の源に対する点 x の応答"
  - symbol: 'G_{\mathrm{ret}}'
    meaning: "遅延 Green 関数"
  - symbol: 'G_{\mathrm{adv}}'
    meaning: "先進 Green 関数"
  - symbol: 'R(z)'
    meaning: "(L-z) の逆作用素"
---

## この章の目標

- 二つの斉次解を連結して境界値問題の核を作る。
- 固有関数展開と鏡像法で同じ逆作用素を表す。
- 放射する波、因果的な波、拡散の基本解を区別する。
- 作用素の摂動を Dyson 方程式で扱う。

[ポテンシャル論](/math/pde/potential-theory/) で扱う調和関数の境界値を、点源の核で表す。時間依存の問題では、同じ空間作用素でも時刻の条件によって核が変わる。

## 常微分方程式の Green 関数の構成

### 点源で導関数が跳ぶ

有限区間 $[a,b]$ で $Lu=-(pu')'+qu$ とし、$p\in C^1([a,b])$ は正、$q$ は実連続とする。両端でそれぞれ $u$ と $pu'$ の実係数の同次一次結合を指定する。各端の係数は同時には零としない。:term[グリーン関数]{en="Green function" yomi="ぐりーんかんすう"} $G(x,y)$ はこの条件で点源 $L_xG=\delta_y$ に応答する核である。零が固有値でないとし、左の条件を満たす非零の斉次解 $u_a$ と右を満たす非零の斉次解 $u_b$ を選ぶ。$W=p(u_au_b'-u_a'u_b)$ とすれば、$W'=u_a(pu_b')'-(pu_a')'u_b=0$ である。$W=0$ なら、ある点で両方の初期データが比例するため二解が従属し、両端を満たす非零解になる。従って $W\ne0$ である。そこで
$$
G(x,y)=-\frac{u_a(\min(x,y))u_b(\max(x,y))}{W}
\label{eq-mp-one-dimensional-green}
$$
は $x=y$ で連続で、$p\partial_xG$ の跳びは $-1$ である。よって式 :eqref[eq-mp-distribution-jump] により $L_xG=\delta_y$ となる。解は $\int_a^bG(x,y)f(y)\dd y$ で、両端の条件も満たす。例えば $L=-\partial_x^2$、$[0,1]$ の Dirichlet では $u_a=x,u_b=1-x$、$W=-1$ なので $G=x_<(1-x_>)$ である。$f=1$ を積分して $u=x(1-x)/2$ を得る。

### 二端のデータを満たす枝を選ぶ

二階作用素の逆を作るとき、左枝は左端の条件を満たす同次解、右枝は右端の条件を満たす同次解から作る。点源で二枝を等しくしなければ、二階微分でデルタの微分が出てしまう。等しくした後、一階微分の跳びを先頭係数に合わせる。二枝が比例してしまう場合は同次問題に非零解があるため、通常の逆は作れない。

:::example[混合端条件と核の積分]
$Lu=-u''$ を $[0,1]$ で考え、$u(0)=u'(1)=0$ とする。$u_a=x$、$u_b=1$、$W=-1$ なので $G(x,y)=\min(x,y)$ である。$f(y)=y$ に作用させると
$$
\int_0^x y^2\,\dd y+x\int_x^1 y\,\dd y=\frac x2-\frac{x^3}{6}
$$
となる。二階微分に負号を付ければ $x$ に戻り、右端では一次微分が消える。
:::

同次でない分離境界条件は、与えられた二つの端条件を満たす滑らかな関数 $v$ を一つ選んで $u=v+G(f-Lv)$ とすれば扱える。例えば $u(a),u(b)$ を指定する場合は線形補間を使える。Robin 条件を含む場合は、各端で必要な値と微分を決め、両端の値と微分を補間する三次多項式を使えばよい。

積分測度に関する対称性も確認できる。実係数では式 :eqref[eq-mp-one-dimensional-green] が $x,y$ に関して対称なので、$G(x,y)=\overline{G(y,x)}$ である。一般の $a_2(x)u''+a_1(x)u'+a_0(x)u$（$a_2$ は連続で零を持たない）では、微分の跳びは $1/a_2(y)$ である。二端を満たす斉次解が従属しないなら核は
$$
\frac{u_a(\min(x,y))u_b(\max(x,y))}
{a_2(y)[u_a(y)u_b'(y)-u_a'(y)u_b(y)]}
$$
となる。分母の $y$ 依存性から、元の測度 $\dd y$ に対する核は一般には対称でない。

:::example[正のポテンシャルを入れた区間の核（1）]
単位区間に $-\partial_x^2+m^2$ と両端 Dirichlet 条件をとる。源は $y=1/2$ に置く。

$$
m=1
$$

左側は原点で零になる双曲正弦を使う。

$$
\begin{aligned}
G_<(x,1/2)\big|_{m=1} &= \frac{\sinh{\left(x \right)}}{2 \cosh{\left(\frac{1}{2} \right)}}
\end{aligned}
$$

右側は右端で零になる双曲正弦に切り替える。

$$
\begin{aligned}
G_>(x,1/2)\big|_{m=1} &= - \frac{\sinh{\left(\frac{1}{2} \right)} \sinh{\left(x - 1 \right)}}{\sinh{\left(1 \right)}}
\end{aligned}
$$

接続点では二枝の値が一致する。

$$
\begin{aligned}
G(1/2,1/2)\big|_{m=1} &= \frac{\sinh^{2}{\left(\frac{1}{2} \right)}}{\sinh{\left(1 \right)}}
\end{aligned}
$$

観測点 $x=\frac{1}{2}$ での値と導関数を計算する。

$$
\begin{aligned}
G_<(x,1/2)\big|_{x=\frac{1}{2}} &= \frac{\tanh{\left(\frac{1}{2} \right)}}{2} \\
\frac{\partial^1 (G_<(x,1/2))}{\partial x^1} &= \frac{\cosh{\left(x \right)}}{2 \cosh{\left(\frac{1}{2} \right)}} \\
\left.\frac{\partial^1 (G_<(x,1/2))}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{1}{2}
\end{aligned}
$$

観測点 $x=\frac{1}{2}$ での値と導関数を計算する。

$$
\begin{aligned}
G_>(x,1/2)\big|_{x=\frac{1}{2}} &= \frac{\sinh^{2}{\left(\frac{1}{2} \right)}}{\sinh{\left(1 \right)}} \\
\frac{\partial^1 (G_>(x,1/2))}{\partial x^1} &= - \frac{\sinh{\left(\frac{1}{2} \right)} \cosh{\left(x - 1 \right)}}{\sinh{\left(1 \right)}} \\
\left.\frac{\partial^1 (G_>(x,1/2))}{\partial x^1}\right|_{x=\frac{1}{2}} &= - \frac{1}{2}
\end{aligned}
$$
:::

:::example[正のポテンシャルを入れた区間の核（2）]
単位区間に $-\partial_x^2+m^2$ と両端 Dirichlet 条件をとる。源は $y=1/2$ に置く。

$$
m=2
$$

左側は原点で零になる双曲正弦を使う。

$$
\begin{aligned}
G_<(x,1/2)\big|_{m=2} &= \frac{\sinh{\left(2 x \right)}}{4 \cosh{\left(1 \right)}}
\end{aligned}
$$

右側は右端で零になる双曲正弦に切り替える。

$$
\begin{aligned}
G_>(x,1/2)\big|_{m=2} &= - \frac{\sinh{\left(1 \right)} \sinh{\left(2 x - 2 \right)}}{2 \sinh{\left(2 \right)}}
\end{aligned}
$$

接続点では二枝の値が一致する。

$$
\begin{aligned}
G(1/2,1/2)\big|_{m=2} &= \frac{\tanh{\left(1 \right)}}{4}
\end{aligned}
$$

観測点 $x=\frac{1}{2}$ での値と導関数を計算する。

$$
\begin{aligned}
G_<(x,1/2)\big|_{x=\frac{1}{2}} &= \frac{\tanh{\left(1 \right)}}{4} \\
\frac{\partial^1 (G_<(x,1/2))}{\partial x^1} &= \frac{\cosh{\left(2 x \right)}}{2 \cosh{\left(1 \right)}} \\
\left.\frac{\partial^1 (G_<(x,1/2))}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{1}{2}
\end{aligned}
$$

観測点 $x=\frac{1}{2}$ での値と導関数を計算する。

$$
\begin{aligned}
G_>(x,1/2)\big|_{x=\frac{1}{2}} &= \frac{\tanh{\left(1 \right)}}{4} \\
\frac{\partial^1 (G_>(x,1/2))}{\partial x^1} &= - \frac{\cosh{\left(2 x - 2 \right)}}{2 \cosh{\left(1 \right)}} \\
\left.\frac{\partial^1 (G_>(x,1/2))}{\partial x^1}\right|_{x=\frac{1}{2}} &= - \frac{1}{2}
\end{aligned}
$$
:::

:::example[正のポテンシャルを入れた区間の核（3）]
単位区間に $-\partial_x^2+m^2$ と両端 Dirichlet 条件をとる。源は $y=1/2$ に置く。

$$
m=3
$$

左側は原点で零になる双曲正弦を使う。

$$
\begin{aligned}
G_<(x,1/2)\big|_{m=3} &= \frac{\sinh{\left(3 x \right)}}{6 \cosh{\left(\frac{3}{2} \right)}}
\end{aligned}
$$

右側は右端で零になる双曲正弦に切り替える。

$$
\begin{aligned}
G_>(x,1/2)\big|_{m=3} &= - \frac{\sinh{\left(\frac{3}{2} \right)} \sinh{\left(3 x - 3 \right)}}{3 \sinh{\left(3 \right)}}
\end{aligned}
$$

接続点では二枝の値が一致する。

$$
\begin{aligned}
G(1/2,1/2)\big|_{m=3} &= \frac{\sinh^{2}{\left(\frac{3}{2} \right)}}{3 \sinh{\left(3 \right)}}
\end{aligned}
$$

観測点 $x=\frac{1}{2}$ での値と導関数を計算する。

$$
\begin{aligned}
G_<(x,1/2)\big|_{x=\frac{1}{2}} &= \frac{\tanh{\left(\frac{3}{2} \right)}}{6} \\
\frac{\partial^1 (G_<(x,1/2))}{\partial x^1} &= \frac{\cosh{\left(3 x \right)}}{2 \cosh{\left(\frac{3}{2} \right)}} \\
\left.\frac{\partial^1 (G_<(x,1/2))}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{1}{2}
\end{aligned}
$$

観測点 $x=\frac{1}{2}$ での値と導関数を計算する。

$$
\begin{aligned}
G_>(x,1/2)\big|_{x=\frac{1}{2}} &= \frac{\sinh^{2}{\left(\frac{3}{2} \right)}}{3 \sinh{\left(3 \right)}} \\
\frac{\partial^1 (G_>(x,1/2))}{\partial x^1} &= - \frac{\sinh{\left(\frac{3}{2} \right)} \cosh{\left(3 x - 3 \right)}}{\sinh{\left(3 \right)}} \\
\left.\frac{\partial^1 (G_>(x,1/2))}{\partial x^1}\right|_{x=\frac{1}{2}} &= - \frac{1}{2}
\end{aligned}
$$
:::

:::example[正のポテンシャルを入れた区間の核（4）]
単位区間に $-\partial_x^2+m^2$ と両端 Dirichlet 条件をとる。源は $y=1/2$ に置く。

$$
m=4
$$

左側は原点で零になる双曲正弦を使う。

$$
\begin{aligned}
G_<(x,1/2)\big|_{m=4} &= \frac{\sinh{\left(4 x \right)}}{8 \cosh{\left(2 \right)}}
\end{aligned}
$$

右側は右端で零になる双曲正弦に切り替える。

$$
\begin{aligned}
G_>(x,1/2)\big|_{m=4} &= - \frac{\sinh{\left(2 \right)} \sinh{\left(4 x - 4 \right)}}{4 \sinh{\left(4 \right)}}
\end{aligned}
$$

接続点では二枝の値が一致する。

$$
\begin{aligned}
G(1/2,1/2)\big|_{m=4} &= \frac{\tanh{\left(2 \right)}}{8}
\end{aligned}
$$

観測点 $x=\frac{1}{2}$ での値と導関数を計算する。

$$
\begin{aligned}
G_<(x,1/2)\big|_{x=\frac{1}{2}} &= \frac{\tanh{\left(2 \right)}}{8} \\
\frac{\partial^1 (G_<(x,1/2))}{\partial x^1} &= \frac{\cosh{\left(4 x \right)}}{2 \cosh{\left(2 \right)}} \\
\left.\frac{\partial^1 (G_<(x,1/2))}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{1}{2}
\end{aligned}
$$

観測点 $x=\frac{1}{2}$ での値と導関数を計算する。

$$
\begin{aligned}
G_>(x,1/2)\big|_{x=\frac{1}{2}} &= \frac{\tanh{\left(2 \right)}}{8} \\
\frac{\partial^1 (G_>(x,1/2))}{\partial x^1} &= - \frac{\sinh{\left(2 \right)} \cosh{\left(4 x - 4 \right)}}{\sinh{\left(4 \right)}} \\
\left.\frac{\partial^1 (G_>(x,1/2))}{\partial x^1}\right|_{x=\frac{1}{2}} &= - \frac{1}{2}
\end{aligned}
$$
:::

## 固有関数展開による Green 関数

### 逆作用素の係数を読む

重み $w>0$ が連続な正則有限区間の自己共役問題 $Lu=w^{-1}[-(pu')'+qu]$ を考える。固有値と正規直交基底を $(\lambda_n,u_n)$ とし、$\inner{f}{g}_w=\int f\bar g\,w\,\dd x$ とする。$z$ がスペクトル外なら各成分を $\lambda_n-z$ で割るので、測度 $w(y)\dd y$ に対する逆核は
$$
G_z(x,y)=\sum_n\frac{u_n(x)\overline{u_n(y)}}{\lambda_n-z}.
\label{eq-mp-spectral-green}
$$
作用素としての収束は $\sum|\inner{f}{u_n}_w|^2<\infty$ と分母の下限による。核を得るには、前節で $q$ を $q-zw$ に替えて得る連続核を、積の正規直交基底 $u_m(x)\overline{u_n(y)}$ へ射影する。$y$ の積分で核を $u_n$ に作用させると $(\lambda_n-z)^{-1}u_n(x)$ になるため、係数は $\delta_{mn}/(\lambda_n-z)$ である。従って上の級数は積測度で $L^2$ 収束し、連続核とほとんど至る所で一致する。分布の式は $w(x)(L_x-z)G_z=\delta_y$ である。$z$ が固有値 $\lambda$ に近づくと、レゾルベントの極の留数は固有空間への直交射影 $-P_\lambda$ であり、単純固有値ではその核は $-u_n(x)\overline{u_n(y)}$ である。

### 小固有値と射影を扱う

核に固有値零があると $L^{-1}$ を全空間で定められない。固有関数へ射影すると、右辺の零モード成分が零でなければ解がないことが分かる。零モードを除く逆を作り、解にも零モードと直交する規格化を課すと一意になる。このとき点源をそのまま入れる代わりに、点源から零モードへの射影を引いたものを右辺に置く。

:::example[単位区間の Neumann 擬逆]
$L=-\partial_x^2$、$u'(0)=u'(1)=0$ とする。定数 $1$ が規格化された零固有関数であり、:term[縮約グリーン関数]{en="reduced Green function" yomi="しゅくやくぐりーんかんすう"} は $-\partial_x^2G_\perp=\delta_y-1$ を満たす。両枝で二階微分を $1$ に取り、Neumann 条件から左枝を $x^2/2+A(y)$、右枝を $x^2/2-x+B(y)$ とする。連続性から $B=A+y$、平均零から $A=1/3+y^2/2-y$ を得る。従って
$$
G_\perp(x,y)=\frac13+\frac{x^2+y^2}{2}-\max(x,y)
$$
である。$\int_0^1 f=0$ のとき、この核を $f$ に作用させれば平均零の解を得る。例えば $f(x)=x-1/2$ では、核を積分すると $u(x)=-x^3/6+x^2/4-1/24$ となり、微分と積分で両端条件と平均零を確認できる。
:::

:::example[Neumann の非零モードの逆（1）]
単位区間の $-\partial_x^2$ を平均零の空間で反転する。

$$
m=1
$$

非零モードは正規化された余弦である。

$$
\begin{aligned}
\phi_m(x)\big|_{m=1} &= \sqrt{2} \cos{\left(\pi x \right)}
\end{aligned}
$$

固有値で割った核の一項を求める。

$$
\begin{aligned}
G_m(x,y)\big|_{m=1} &= \frac{2 \cos{\left(\pi x \right)} \cos{\left(\pi y \right)}}{\pi^{2}}
\end{aligned}
$$

点源の代わりに、$y=0$ に近づけたときの成分を境界からの極限として評価する。

$$
\begin{aligned}
G_m(x,0)\big|_{m=1} &= \frac{2 \cos{\left(\pi x \right)}}{\pi^{2}}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
\phi_m(x)\big|_{x=\frac{1}{3}} &= \frac{\sqrt{2}}{2} \\
\frac{\partial^1 (\phi_m(x))}{\partial x^1} &= - \sqrt{2} \pi \sin{\left(\pi x \right)} \\
\left.\frac{\partial^1 (\phi_m(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - \frac{\sqrt{6} \pi}{2}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
G_m(x,y)\big|_{x=\frac{1}{3}} &= \frac{\cos{\left(\pi y \right)}}{\pi^{2}} \\
\frac{\partial^1 (G_m(x,y))}{\partial x^1} &= - \frac{2 \sin{\left(\pi x \right)} \cos{\left(\pi y \right)}}{\pi} \\
\left.\frac{\partial^1 (G_m(x,y))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - \frac{\sqrt{3} \cos{\left(\pi y \right)}}{\pi}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
G_m(x,0)\big|_{x=\frac{1}{3}} &= \frac{1}{\pi^{2}} \\
\frac{\partial^1 (G_m(x,0))}{\partial x^1} &= - \frac{2 \sin{\left(\pi x \right)}}{\pi} \\
\left.\frac{\partial^1 (G_m(x,0))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - \frac{\sqrt{3}}{\pi}
\end{aligned}
$$
:::

:::example[Neumann の非零モードの逆（2）]
単位区間の $-\partial_x^2$ を平均零の空間で反転する。

$$
m=2
$$

非零モードは正規化された余弦である。

$$
\begin{aligned}
\phi_m(x)\big|_{m=2} &= \sqrt{2} \cos{\left(2 \pi x \right)}
\end{aligned}
$$

固有値で割った核の一項を求める。

$$
\begin{aligned}
G_m(x,y)\big|_{m=2} &= \frac{\cos{\left(2 \pi x \right)} \cos{\left(2 \pi y \right)}}{2 \pi^{2}}
\end{aligned}
$$

点源の代わりに、$y=0$ に近づけたときの成分を境界からの極限として評価する。

$$
\begin{aligned}
G_m(x,0)\big|_{m=2} &= \frac{\cos{\left(2 \pi x \right)}}{2 \pi^{2}}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
\phi_m(x)\big|_{x=\frac{1}{3}} &= - \frac{\sqrt{2}}{2} \\
\frac{\partial^1 (\phi_m(x))}{\partial x^1} &= - 2 \sqrt{2} \pi \sin{\left(2 \pi x \right)} \\
\left.\frac{\partial^1 (\phi_m(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - \sqrt{6} \pi
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
G_m(x,y)\big|_{x=\frac{1}{3}} &= - \frac{\cos{\left(2 \pi y \right)}}{4 \pi^{2}} \\
\frac{\partial^1 (G_m(x,y))}{\partial x^1} &= - \frac{\sin{\left(2 \pi x \right)} \cos{\left(2 \pi y \right)}}{\pi} \\
\left.\frac{\partial^1 (G_m(x,y))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - \frac{\sqrt{3} \cos{\left(2 \pi y \right)}}{2 \pi}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
G_m(x,0)\big|_{x=\frac{1}{3}} &= - \frac{1}{4 \pi^{2}} \\
\frac{\partial^1 (G_m(x,0))}{\partial x^1} &= - \frac{\sin{\left(2 \pi x \right)}}{\pi} \\
\left.\frac{\partial^1 (G_m(x,0))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - \frac{\sqrt{3}}{2 \pi}
\end{aligned}
$$
:::

:::example[Neumann の非零モードの逆（3）]
単位区間の $-\partial_x^2$ を平均零の空間で反転する。

$$
m=3
$$

非零モードは正規化された余弦である。

$$
\begin{aligned}
\phi_m(x)\big|_{m=3} &= \sqrt{2} \cos{\left(3 \pi x \right)}
\end{aligned}
$$

固有値で割った核の一項を求める。

$$
\begin{aligned}
G_m(x,y)\big|_{m=3} &= \frac{2 \cos{\left(3 \pi x \right)} \cos{\left(3 \pi y \right)}}{9 \pi^{2}}
\end{aligned}
$$

点源の代わりに、$y=0$ に近づけたときの成分を境界からの極限として評価する。

$$
\begin{aligned}
G_m(x,0)\big|_{m=3} &= \frac{2 \cos{\left(3 \pi x \right)}}{9 \pi^{2}}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
\phi_m(x)\big|_{x=\frac{1}{3}} &= - \sqrt{2} \\
\frac{\partial^1 (\phi_m(x))}{\partial x^1} &= - 3 \sqrt{2} \pi \sin{\left(3 \pi x \right)} \\
\left.\frac{\partial^1 (\phi_m(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= 0
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
G_m(x,y)\big|_{x=\frac{1}{3}} &= - \frac{2 \cos{\left(3 \pi y \right)}}{9 \pi^{2}} \\
\frac{\partial^1 (G_m(x,y))}{\partial x^1} &= - \frac{2 \sin{\left(3 \pi x \right)} \cos{\left(3 \pi y \right)}}{3 \pi} \\
\left.\frac{\partial^1 (G_m(x,y))}{\partial x^1}\right|_{x=\frac{1}{3}} &= 0
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
G_m(x,0)\big|_{x=\frac{1}{3}} &= - \frac{2}{9 \pi^{2}} \\
\frac{\partial^1 (G_m(x,0))}{\partial x^1} &= - \frac{2 \sin{\left(3 \pi x \right)}}{3 \pi} \\
\left.\frac{\partial^1 (G_m(x,0))}{\partial x^1}\right|_{x=\frac{1}{3}} &= 0
\end{aligned}
$$
:::

:::example[Neumann の非零モードの逆（4）]
単位区間の $-\partial_x^2$ を平均零の空間で反転する。

$$
m=4
$$

非零モードは正規化された余弦である。

$$
\begin{aligned}
\phi_m(x)\big|_{m=4} &= \sqrt{2} \cos{\left(4 \pi x \right)}
\end{aligned}
$$

固有値で割った核の一項を求める。

$$
\begin{aligned}
G_m(x,y)\big|_{m=4} &= \frac{\cos{\left(4 \pi x \right)} \cos{\left(4 \pi y \right)}}{8 \pi^{2}}
\end{aligned}
$$

点源の代わりに、$y=0$ に近づけたときの成分を境界からの極限として評価する。

$$
\begin{aligned}
G_m(x,0)\big|_{m=4} &= \frac{\cos{\left(4 \pi x \right)}}{8 \pi^{2}}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
\phi_m(x)\big|_{x=\frac{1}{3}} &= - \frac{\sqrt{2}}{2} \\
\frac{\partial^1 (\phi_m(x))}{\partial x^1} &= - 4 \sqrt{2} \pi \sin{\left(4 \pi x \right)} \\
\left.\frac{\partial^1 (\phi_m(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= 2 \sqrt{6} \pi
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
G_m(x,y)\big|_{x=\frac{1}{3}} &= - \frac{\cos{\left(4 \pi y \right)}}{16 \pi^{2}} \\
\frac{\partial^1 (G_m(x,y))}{\partial x^1} &= - \frac{\sin{\left(4 \pi x \right)} \cos{\left(4 \pi y \right)}}{2 \pi} \\
\left.\frac{\partial^1 (G_m(x,y))}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{\sqrt{3} \cos{\left(4 \pi y \right)}}{4 \pi}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
G_m(x,0)\big|_{x=\frac{1}{3}} &= - \frac{1}{16 \pi^{2}} \\
\frac{\partial^1 (G_m(x,0))}{\partial x^1} &= - \frac{\sin{\left(4 \pi x \right)}}{2 \pi} \\
\left.\frac{\partial^1 (G_m(x,0))}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{\sqrt{3}}{4 \pi}
\end{aligned}
$$
:::

## Poisson 方程式の Green 関数と鏡像法

### 境界を打ち消す像

$-\Delta$ の全空間基本解は :ref[exm-mp-fundamental-laplace] の $g(x-y)=1/(4\pi|x-y|)$ である。上半空間の Dirichlet で $y^*$ を $y$ の鏡映像とすると
$$
G_D(x,y)=\frac1{4\pi|x-y|}-\frac1{4\pi|x-y^*|}
$$
となる。像は領域の外なので内部の源を変えず、境界面では距離が等しいため零である。零 Neumann では符号を正に替える。実際、法線方向で互いの微分符号が逆だからである。球 $|x|<a$ の Dirichlet では $0<|y|<a$ に対して $y^*=a^2 y/|y|^2$、像の係数を $-a/|y|$ とする。境界では $|x-y^*|=(a/|y|)|x-y|$ が、両辺を二乗して確かめられる。これにより源、境界の両方を満たし、差の調和関数に最大値原理を使うと一意である。中心の源の極限では $G_D=1/(4\pi r)-1/(4\pi a)$ となる。

境界値も非零のとき、$x$ を積分変数として Green の第 2 恒等式を使う。$-\Delta u=f$ とし、滑らかな有界領域から $y$ のまわりの小球を除く。小球の外側で $G_D$ は調和なので、$\int G_Df=\int_{\partial}(u\partial_nG_D-G_D\partial_nu)$ となる。小球面で領域の外向き法線が $-\partial_r$ を向くため、$u\partial_nG_D$ の積分は $u(y)$ に収束する。$G_D\partial_nu$ の積分は $O(r)$ で零に収束する。外側の境界では $G_D=0$ なので、極限で
$$
u(y)=\int_\Omega G_D(x,y)f(x)\,\dd x
-\int_{\partial\Omega}u(x)\partial_{n_x}G_D(x,y)\,\dd S_x
$$
を得る。第二項の核 $-\partial_{n_x}G_D$ を :term[ポアソン核]{en="Poisson kernel" yomi="ぽあそんかく"} と呼ぶ。

:::example[球の Poisson 核と中心値]
上の Kelvin 像を外向き $\partial_{n_x}$ で微分する。$x=a\widehat x$ では $|x-y^*|=(a/|y|)|x-y|$ なので、二項を同じ分母にそろえると
$$
P(y,x)=\frac{a^2-|y|^2}{4\pi a|x-y|^3},\qquad |x|=a
$$
となる。実際、第一項の分子は $a-\widehat x\cdot y$、像から差し引く分子は $|y|^2/a-\widehat x\cdot y$ である。$y=0$ では $P=1/(4\pi a^2)$ となり、調和関数の中心値は球面平均に等しい。全ての $y$ で $\int_{\partial\Omega}P(y,x)\dd S_x=1$ であることは、境界値 $1$ の調和解が定数 $1$ であることからも分かる。
:::

### 平面の点源は対数になる

二次元では $-\Delta$ の基本解を $g(r)=-(2\pi)^{-1}\log(r/r_0)$ とできる。$r>0$ で $r g'$ は定数なので調和である。原点での質量は小円の法線微分から $-\int_{r=\epsilon}\partial_r g\,\dd S=1$ となる。残りの境界項は $\epsilon|\log\epsilon|$ に比例して消えるので、試験関数での作用も点源である。$r_0$ は対数の引数を無次元にする任意の基準長であり、変えると定数だけ変わる。例えば半径 $a$ の円板で中心に点源を置き外周を零とする解は $(2\pi)^{-1}\log(a/r)$ である。$r=a/e$ で $1/(2\pi)$ となる。三次元の $1/r$ をそのまま平面に使うと、原点以外でも Laplacian が零にならない。

### 平面境界に反対符号を置く

三次元の半空間 $z>0$ で Dirichlet 条件を課す。高さ $a$ の正点源に対し、高さ $-a$ に反対符号の像を置けば平面上で逆距離が一致して相殺する。像の位置は領域外なので実際の方程式の点源は変わらない。面上の垂直微分は、二点源の距離の変化が逆向きであるため逆に加算される。

:::example[半空間の軸上の像（1）]
点源を高さ $m$ に置き、軸上の点を高さ $x$ で表す。

$$
m=1
$$

源より下の領域 $0<x<m$ では距離が $m-x$ と $m+x$ になる。

$$
\begin{aligned}
G(x)\big|_{m=1} &= - \frac{x}{2 \pi \left(x^{2} - 1\right)}
\end{aligned}
$$

境界での内向き法線微分を求める。

$$
\begin{aligned}
\partial_xG(0)\big|_{m=1} &= \frac{1}{2 \pi}
\end{aligned}
$$

軸から水平方向に距離 $\rho$ 離れた境界での同じ微分を求める。

$$
\begin{aligned}
\partial_zG(\rho,0)\big|_{m=1} &= \frac{1}{2 \pi \left(\rho^{2} + 1\right)^{\frac{3}{2}}}
\end{aligned}
$$

観測点 $x=\frac{1}{2}$ での値と導関数を計算する。

$$
\begin{aligned}
G(x)\big|_{x=\frac{1}{2}} &= \frac{1}{3 \pi} \\
\frac{\partial^1 (G(x))}{\partial x^1} &= \frac{x^{2} + 1}{2 \pi \left(x^{2} - 1\right)^{2}} \\
\left.\frac{\partial^1 (G(x))}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{10}{9 \pi}
\end{aligned}
$$
:::

:::example[半空間の軸上の像（2）]
点源を高さ $m$ に置き、軸上の点を高さ $x$ で表す。

$$
m=2
$$

源より下の領域 $0<x<m$ では距離が $m-x$ と $m+x$ になる。

$$
\begin{aligned}
G(x)\big|_{m=2} &= - \frac{x}{2 \pi \left(x^{2} - 4\right)}
\end{aligned}
$$

境界での内向き法線微分を求める。

$$
\begin{aligned}
\partial_xG(0)\big|_{m=2} &= \frac{1}{8 \pi}
\end{aligned}
$$

軸から水平方向に距離 $\rho$ 離れた境界での同じ微分を求める。

$$
\begin{aligned}
\partial_zG(\rho,0)\big|_{m=2} &= \frac{1}{\pi \left(\rho^{2} + 4\right)^{\frac{3}{2}}}
\end{aligned}
$$

観測点 $x=\frac{1}{2}$ での値と導関数を計算する。

$$
\begin{aligned}
G(x)\big|_{x=\frac{1}{2}} &= \frac{1}{15 \pi} \\
\frac{\partial^1 (G(x))}{\partial x^1} &= \frac{x^{2} + 4}{2 \pi \left(x^{2} - 4\right)^{2}} \\
\left.\frac{\partial^1 (G(x))}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{34}{225 \pi}
\end{aligned}
$$
:::

:::example[半空間の軸上の像（3）]
点源を高さ $m$ に置き、軸上の点を高さ $x$ で表す。

$$
m=3
$$

源より下の領域 $0<x<m$ では距離が $m-x$ と $m+x$ になる。

$$
\begin{aligned}
G(x)\big|_{m=3} &= - \frac{x}{2 \pi \left(x^{2} - 9\right)}
\end{aligned}
$$

境界での内向き法線微分を求める。

$$
\begin{aligned}
\partial_xG(0)\big|_{m=3} &= \frac{1}{18 \pi}
\end{aligned}
$$

軸から水平方向に距離 $\rho$ 離れた境界での同じ微分を求める。

$$
\begin{aligned}
\partial_zG(\rho,0)\big|_{m=3} &= \frac{3}{2 \pi \left(\rho^{2} + 9\right)^{\frac{3}{2}}}
\end{aligned}
$$

観測点 $x=\frac{1}{2}$ での値と導関数を計算する。

$$
\begin{aligned}
G(x)\big|_{x=\frac{1}{2}} &= \frac{1}{35 \pi} \\
\frac{\partial^1 (G(x))}{\partial x^1} &= \frac{x^{2} + 9}{2 \pi \left(x^{2} - 9\right)^{2}} \\
\left.\frac{\partial^1 (G(x))}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{74}{1225 \pi}
\end{aligned}
$$
:::

:::example[半空間の軸上の像（4）]
点源を高さ $m$ に置き、軸上の点を高さ $x$ で表す。

$$
m=4
$$

源より下の領域 $0<x<m$ では距離が $m-x$ と $m+x$ になる。

$$
\begin{aligned}
G(x)\big|_{m=4} &= - \frac{x}{2 \pi \left(x^{2} - 16\right)}
\end{aligned}
$$

境界での内向き法線微分を求める。

$$
\begin{aligned}
\partial_xG(0)\big|_{m=4} &= \frac{1}{32 \pi}
\end{aligned}
$$

軸から水平方向に距離 $\rho$ 離れた境界での同じ微分を求める。

$$
\begin{aligned}
\partial_zG(\rho,0)\big|_{m=4} &= \frac{2}{\pi \left(\rho^{2} + 16\right)^{\frac{3}{2}}}
\end{aligned}
$$

観測点 $x=\frac{1}{2}$ での値と導関数を計算する。

$$
\begin{aligned}
G(x)\big|_{x=\frac{1}{2}} &= \frac{1}{63 \pi} \\
\frac{\partial^1 (G(x))}{\partial x^1} &= \frac{x^{2} + 16}{2 \pi \left(x^{2} - 16\right)^{2}} \\
\left.\frac{\partial^1 (G(x))}{\partial x^1}\right|_{x=\frac{1}{2}} &= \frac{130}{3969 \pi}
\end{aligned}
$$
:::

## Helmholtz 方程式の Green 関数と放射条件

### 波を外へ出す逆作用素

三次元で $k>0$ とし、時間因子に $e^{-i\omega t}$（$\omega>0$）を取る。$-\Delta-k^2$ の外向き基本解は $g_k(r)=e^{ikr}/(4\pi r)$ である。$r>0$ では $\Delta(v/r)=v''/r$ であり、$v''=-k^2v$ を満たす。原点で $1/(4\pi r)$ と同じ $1/r$ 特異性を持ち、差の微分は有界なので、小球の境界項の極限は同じ $\delta$ である。遠方では $(\partial_r-ik)g_k=-e^{ikr}/(4\pi r^2)$ なので Sommerfeld 条件を満たす。$k+i\epsilon$（$\epsilon>0$）を使うと空間無限遠で指数減衰する。この核を Fourier 変換すれば $(|\xi|^2-(k+i\epsilon)^2)^{-1}$ であり、$\epsilon\downarrow0$ が外向きの選択に対応する。内向き核 $e^{-ikr}/(4\pi r)$ は、減衰する $e^{-i(k-i\epsilon)r}/(4\pi r)$ の極限である。

### 時間の核から周波数の核へ

時間にも空間と同じ Fourier の符号規約を用いる。時間因子 $e^{-i\omega t}$ に対応する係数は、時間 Fourier 変換を周波数 $-\omega$ で評価した $\int e^{i\omega t}G_{\rm ret}(t,r)\dd t$ である。三次元の核は次節で導くように時刻 $r/c$ に台を持つから、結果は $e^{i\omega r/c}/(4\pi c^2r)$ である。作用素 $\partial_t^2-c^2\Delta$ も変換で $c^2(-\Delta-\omega^2/c^2)$ となるため、$\omega>0$ では $1/c^2$ を除けば外向き核と一致する。負周波数では正周波数の核の複素共役が現れ、時間因子も共役になる。例えば $r/c=1\,\mathrm{s}$ で $\omega=2\pi\,\mathrm{s^{-1}}$ なら位相は $1$、$\omega=\pi\,\mathrm{s^{-1}}$ なら $-1$ である。

### 放射条件による一意性

三次元で滑らかなコンパクト障害物の外部の連結成分を考える。全ての方向で一様に $u=O(r^{-1})$、$(\partial_r-ik)u=o(r^{-1})$ とする。同じ源と Dirichlet データを持つ二解の差 $v$ は領域で $(\Delta+k^2)v=0$、障害物表面で零である。障害物を包む半径 $R$ の球で領域を切り、式に $\bar v$ を掛けて部分積分する。体積積分が実数なので
$$
\Im\int_{|x|=R}\bar v\,\partial_rv\,\dd S=0
$$
である。これに放射条件を代入すると $k\int_{|x|=R}|v|^2\,\dd S=o(1)$ となる。一方、障害物を囲む球の外側で各球面調和関数に射影すれば、動径係数は球 Hankel 関数の線形結合になる。射影した放射条件のもとで内向き係数は零である。各外向き係数 $a_{\ell m}$ は、$R^2|a_{\ell m}h_\ell^{(1)}(kR)|^2\to |a_{\ell m}|^2/k^2$ を満たす。角方向の直交性からこの量は $\int_{|x|=R}|v|^2\,\dd S$ 以下なので全係数が零である。従って $v$ は十分遠方で零で、連結な外部領域では Helmholtz 方程式の解の実解析性により全域で零である。ここで実解析性は定数係数の楕円方程式に対する結果である。$(\Delta+k^2)v=0$ なら $v(x)e^{ks}$ が一変数加えた空間の調和関数になるため、この場合は調和関数の実解析性だけでも従う。

### 二次元の放射核

二次元で中心対称な放射解は $H_0^{(1)}(kr)$ の定数倍である。$J_0(r)=1+O(r^2)$ と $W(J_0,Y_0)=2/(\pi r)$ を使うと $(Y_0/J_0)'=2/(\pi r J_0^2)$ で、$Y_0(r)=(2/\pi)\log r+O(1)$ が分かる。従って $(i/4)H_0^{(1)}(kr)$ の対数特異性は $-(2\pi)^{-1}\log r$ である。小円の計算により $(-\Delta-k^2)$ を作用させると $\delta^{(2)}$ を得る。二次元の遠方振幅は $r^{-1/2}$、三次元では $r^{-1}$ である。これらの二乗を円周または球面の大きさと掛けると遠方で次数が零となり、定常放射で面全体の流束を保てる。

### 放射条件の差を直接調べる

時間因子を $e^{-i\omega t}$ に取る。球面波 $g=e^{ikr}/(4\pi r)$ に $\partial_r-ik$ を作用させると、指数の微分は消えて $-e^{ikr}/(4\pi r^2)$ が残る。このため $r(\partial_r-ik)g\to0$ となる。内向きの $e^{-ikr}$ を入れると $-2ik$ の振動成分が残る。規約を $e^{+i\omega t}$ に変える場合は放射条件の符号も反対になる。

:::example[球面波の二つの選択（1）]
波数を $m$ として放射作用素の結果を比較する。

$$
m=1
$$

外向き波に残る項を計算する。

$$
\begin{aligned}
(\partial_r-im)g_{\rm out}\big|_{m=1} &= - \frac{e^{i r}}{4 \pi r^{2}}
\end{aligned}
$$

内向き波では $r^{-1}$ と $r^{-2}$ の両方が残る。

$$
\begin{aligned}
(\partial_r-im)g_{\rm in}\big|_{m=1} &= \frac{\left(- 2 i r - 1\right) e^{- i r}}{4 \pi r^{2}}
\end{aligned}
$$

外向き波の二乗を半径 $r$ の球面全体で積分する。距離とともに場が減る一方、面積は増える。

$$
\begin{aligned}
4\pi r^2|g_{\rm out}|^2 &= \frac{1}{4 \pi}
\end{aligned}
$$

観測点 $r=2$ での値と導関数を計算する。

$$
\begin{aligned}
(\partial_r-im)g_{\rm out}\big|_{r=2} &= - \frac{e^{2 i}}{16 \pi} \\
\frac{\partial^1 ((\partial_r-im)g_{\rm out})}{\partial r^1} &= \frac{\left(- i r + 2\right) e^{i r}}{4 \pi r^{3}} \\
\left.\frac{\partial^1 ((\partial_r-im)g_{\rm out})}{\partial r^1}\right|_{r=2} &= \frac{\left(1 - i\right) e^{2 i}}{16 \pi}
\end{aligned}
$$

観測点 $r=2$ での値と導関数を計算する。

$$
\begin{aligned}
(\partial_r-im)g_{\rm in}\big|_{r=2} &= \frac{\left(-1 - 4 i\right) e^{- 2 i}}{16 \pi} \\
\frac{\partial^1 ((\partial_r-im)g_{\rm in})}{\partial r^1} &= \frac{\left(- 2 r^{2} + 3 i r + 2\right) e^{- i r}}{4 \pi r^{3}} \\
\left.\frac{\partial^1 ((\partial_r-im)g_{\rm in})}{\partial r^1}\right|_{r=2} &= \frac{3 \left(-1 + i\right) e^{- 2 i}}{16 \pi}
\end{aligned}
$$
:::

:::example[球面波の二つの選択（2）]
波数を $m$ として放射作用素の結果を比較する。

$$
m=2
$$

外向き波に残る項を計算する。

$$
\begin{aligned}
(\partial_r-im)g_{\rm out}\big|_{m=2} &= - \frac{e^{2 i r}}{4 \pi r^{2}}
\end{aligned}
$$

内向き波では $r^{-1}$ と $r^{-2}$ の両方が残る。

$$
\begin{aligned}
(\partial_r-im)g_{\rm in}\big|_{m=2} &= \frac{\left(- i r - \frac{1}{4}\right) e^{- 2 i r}}{\pi r^{2}}
\end{aligned}
$$

外向き波の二乗を半径 $r$ の球面全体で積分する。距離とともに場が減る一方、面積は増える。

$$
\begin{aligned}
4\pi r^2|g_{\rm out}|^2 &= \frac{1}{4 \pi}
\end{aligned}
$$

観測点 $r=2$ での値と導関数を計算する。

$$
\begin{aligned}
(\partial_r-im)g_{\rm out}\big|_{r=2} &= - \frac{e^{4 i}}{16 \pi} \\
\frac{\partial^1 ((\partial_r-im)g_{\rm out})}{\partial r^1} &= \frac{\left(- i r + 1\right) e^{2 i r}}{2 \pi r^{3}} \\
\left.\frac{\partial^1 ((\partial_r-im)g_{\rm out})}{\partial r^1}\right|_{r=2} &= \frac{\left(1 - 2 i\right) e^{4 i}}{16 \pi}
\end{aligned}
$$

観測点 $r=2$ での値と導関数を計算する。

$$
\begin{aligned}
(\partial_r-im)g_{\rm in}\big|_{r=2} &= \frac{\left(-1 - 8 i\right) e^{- 4 i}}{16 \pi} \\
\frac{\partial^1 ((\partial_r-im)g_{\rm in})}{\partial r^1} &= \frac{\left(- 4 r^{2} + 3 i r + 1\right) e^{- 2 i r}}{2 \pi r^{3}} \\
\left.\frac{\partial^1 ((\partial_r-im)g_{\rm in})}{\partial r^1}\right|_{r=2} &= \frac{3 \left(-5 + 2 i\right) e^{- 4 i}}{16 \pi}
\end{aligned}
$$
:::

:::example[球面波の二つの選択（3）]
波数を $m$ として放射作用素の結果を比較する。

$$
m=3
$$

外向き波に残る項を計算する。

$$
\begin{aligned}
(\partial_r-im)g_{\rm out}\big|_{m=3} &= - \frac{e^{3 i r}}{4 \pi r^{2}}
\end{aligned}
$$

内向き波では $r^{-1}$ と $r^{-2}$ の両方が残る。

$$
\begin{aligned}
(\partial_r-im)g_{\rm in}\big|_{m=3} &= \frac{\left(- 6 i r - 1\right) e^{- 3 i r}}{4 \pi r^{2}}
\end{aligned}
$$

外向き波の二乗を半径 $r$ の球面全体で積分する。距離とともに場が減る一方、面積は増える。

$$
\begin{aligned}
4\pi r^2|g_{\rm out}|^2 &= \frac{1}{4 \pi}
\end{aligned}
$$

観測点 $r=2$ での値と導関数を計算する。

$$
\begin{aligned}
(\partial_r-im)g_{\rm out}\big|_{r=2} &= - \frac{e^{6 i}}{16 \pi} \\
\frac{\partial^1 ((\partial_r-im)g_{\rm out})}{\partial r^1} &= \frac{\left(- 3 i r + 2\right) e^{3 i r}}{4 \pi r^{3}} \\
\left.\frac{\partial^1 ((\partial_r-im)g_{\rm out})}{\partial r^1}\right|_{r=2} &= \frac{\left(1 - 3 i\right) e^{6 i}}{16 \pi}
\end{aligned}
$$

観測点 $r=2$ での値と導関数を計算する。

$$
\begin{aligned}
(\partial_r-im)g_{\rm in}\big|_{r=2} &= \frac{\left(-1 - 12 i\right) e^{- 6 i}}{16 \pi} \\
\frac{\partial^1 ((\partial_r-im)g_{\rm in})}{\partial r^1} &= \frac{\left(- 18 r^{2} + 9 i r + 2\right) e^{- 3 i r}}{4 \pi r^{3}} \\
\left.\frac{\partial^1 ((\partial_r-im)g_{\rm in})}{\partial r^1}\right|_{r=2} &= \frac{\left(-35 + 9 i\right) e^{- 6 i}}{16 \pi}
\end{aligned}
$$
:::

:::example[球面波の二つの選択（4）]
波数を $m$ として放射作用素の結果を比較する。

$$
m=4
$$

外向き波に残る項を計算する。

$$
\begin{aligned}
(\partial_r-im)g_{\rm out}\big|_{m=4} &= - \frac{e^{4 i r}}{4 \pi r^{2}}
\end{aligned}
$$

内向き波では $r^{-1}$ と $r^{-2}$ の両方が残る。

$$
\begin{aligned}
(\partial_r-im)g_{\rm in}\big|_{m=4} &= \frac{\left(- 8 i r - 1\right) e^{- 4 i r}}{4 \pi r^{2}}
\end{aligned}
$$

外向き波の二乗を半径 $r$ の球面全体で積分する。距離とともに場が減る一方、面積は増える。

$$
\begin{aligned}
4\pi r^2|g_{\rm out}|^2 &= \frac{1}{4 \pi}
\end{aligned}
$$

観測点 $r=2$ での値と導関数を計算する。

$$
\begin{aligned}
(\partial_r-im)g_{\rm out}\big|_{r=2} &= - \frac{e^{8 i}}{16 \pi} \\
\frac{\partial^1 ((\partial_r-im)g_{\rm out})}{\partial r^1} &= \frac{\left(- i r + \frac{1}{2}\right) e^{4 i r}}{\pi r^{3}} \\
\left.\frac{\partial^1 ((\partial_r-im)g_{\rm out})}{\partial r^1}\right|_{r=2} &= \frac{\left(1 - 4 i\right) e^{8 i}}{16 \pi}
\end{aligned}
$$

観測点 $r=2$ での値と導関数を計算する。

$$
\begin{aligned}
(\partial_r-im)g_{\rm in}\big|_{r=2} &= \frac{\left(- \frac{1}{16} - i\right) e^{- 8 i}}{\pi} \\
\frac{\partial^1 ((\partial_r-im)g_{\rm in})}{\partial r^1} &= \frac{\left(- 16 r^{2} + 6 i r + 1\right) e^{- 4 i r}}{2 \pi r^{3}} \\
\left.\frac{\partial^1 ((\partial_r-im)g_{\rm in})}{\partial r^1}\right|_{r=2} &= \frac{3 \left(-21 + 4 i\right) e^{- 8 i}}{16 \pi}
\end{aligned}
$$
:::

## 波動方程式の遅延・先進 Green 関数

### 波が届く時刻に台を持たせる

$L=\partial_t^2-c^2\Delta$ とする。空間 Fourier 変換で $Lg=\delta(t)\delta(x)$ を解くと、遅延条件 $t<0$ で零の :term[遅延グリーン関数]{en="retarded Green function" yomi="ちえんぐりーんかんすう"} は $\widehat G_{\rm ret}(t,k)=H(t)\sin(c|k|t)/(c|k|)$ である。$t=0$ の値は零、時間微分の跳びは $1$ だからである。三次元で固定 $t>0$ の球面 $r=ct$ 上に
$$
G_{\rm ret}(t,r)=\frac{\delta(t-r/c)}{4\pi c^2r},\qquad
G_{\rm adv}(t,r)=\frac{\delta(t+r/c)}{4\pi c^2r}
\label{eq-mp-wave-green}
$$
と置くと、角度の Fourier 平均は $\sin(|k|r)/(|k|r)$、動径積分で $r=ct$ に取り出すから、遅延核の変換は上式と一致する。先進核は時間反転である。$t=0$ の点も含めた超関数は、$4\pi c^2r$ を分母にした球面積分の式で定まり、原点での積分は局所有限である。一次元では遅延核は $H(t-|x|/c)/(2c)$ で、区間の Fourier 変換で同じ乗数を与える。三次元では光円錐の表面、一次元では内部にも応答がある。

### 波面上のデルタを時間で積分する

三次元の遅延核では正の距離 $r$ の点へ源の変化が時刻 $r/c$ だけ遅れて届く。時刻 $t$ の観測では、時間畳み込みのデルタが $s=t-r/c$ を選ぶ。外力が $s=0$ より前に零なら、到着時刻まで応答は零である。源に時間微分が入る場合も、先に因果的に延長してから微分する。初期の跳びがあるとその微分にデルタが現れる。

### 次元降下と初期値の寄与

三次元の遅延核を一つの空間座標 $z$ で積分すると、二次元の遅延核になる。任意のコンパクトな時刻区間では三次元核の空間の台が有界なので、試験関数を $z$ 方向へ台の外まで定数で延長して三次元の方程式を適用できる。$r$ を二次元の動径距離とすると、$ct>r$ では $\sqrt{r^2+z^2}=ct$ の根 $z=\pm\sqrt{c^2t^2-r^2}$ が寄与する。変数変換の係数はそれぞれ $c^2t/|z|$ だから
$$
G_{\mathrm{ret}}^{(2)}(t,r)=
\frac{H(t-r/c)}{2\pi c\sqrt{c^2t^2-r^2}}
$$
となる。右辺は $ct>r$ でこの値、$ct<r$ で零とする局所可積分関数である。二次元では源が止まった後も円錐の内部に応答が残る。

外力 $F$ のほかに初期変位 $a(x)$、速度 $b(x)$ がある場合、正時刻の解を $t<0$ で零に延長すると
$$
(\partial_t^2-c^2\Delta)(H u)=H F+\delta(t)b+\delta'(t)a
$$
となる。実際、一回目の微分で $\delta(t)a$ が現れ、二回目で $\delta'(t)a$ と $\delta(t)b$ が現れる。この式に遅延核を畳み込めば、$t>0$ の解は $G_{\mathrm{ret}}*(HF)+G_{\mathrm{ret}}(t,\cdot)*_x b+\partial_tG_{\mathrm{ret}}(t,\cdot)*_x a$ となる。ここで $*_x$ は空間だけの畳み込みである。

:::example[調和源を時刻零で入れる（1）]
波速を $1$、原点の点源を $H(t)\cos(mt)\delta^{(3)}(\vb x)$ とする。本文の波動演算子で $t>r$ の場を求める。

$$
m=1
$$

遅延時刻を計算する。

$$
\begin{aligned}
t_{\rm emit} &= - r + t
\end{aligned}
$$

デルタで時間源を評価する。

$$
\begin{aligned}
u(t,r)\big|_{m=1} &= \frac{\cos{\left(r - t \right)}}{4 \pi r}
\end{aligned}
$$

到達直後の初期値には球面の広がりの因子が残る。

$$
\begin{aligned}
u(r+0,r) &= \frac{1}{4 \pi r}
\end{aligned}
$$

半径 $0<r<3$ に対して $t=3$ の値と、到着後の時間微分を求める。

$$
\begin{aligned}
t_{\rm emit}\big|_{t=3} &= 3 - r \\
\frac{\partial^1 (t_{\rm emit})}{\partial t^1} &= 1 \\
\left.\frac{\partial^1 (t_{\rm emit})}{\partial t^1}\right|_{t=3} &= 1
\end{aligned}
$$

半径 $0<r<3$ に対して $t=3$ の値と、到着後の時間微分を求める。

$$
\begin{aligned}
u(t,r)\big|_{t=3} &= \frac{\cos{\left(r - 3 \right)}}{4 \pi r} \\
\frac{\partial^1 (u(t,r))}{\partial t^1} &= \frac{\sin{\left(r - t \right)}}{4 \pi r} \\
\left.\frac{\partial^1 (u(t,r))}{\partial t^1}\right|_{t=3} &= \frac{\sin{\left(r - 3 \right)}}{4 \pi r}
\end{aligned}
$$
:::

:::example[調和源を時刻零で入れる（2）]
波速を $1$、原点の点源を $H(t)\cos(mt)\delta^{(3)}(\vb x)$ とする。本文の波動演算子で $t>r$ の場を求める。

$$
m=2
$$

遅延時刻を計算する。

$$
\begin{aligned}
t_{\rm emit} &= - r + t
\end{aligned}
$$

デルタで時間源を評価する。

$$
\begin{aligned}
u(t,r)\big|_{m=2} &= \frac{\cos{\left(2 r - 2 t \right)}}{4 \pi r}
\end{aligned}
$$

到達直後の初期値には球面の広がりの因子が残る。

$$
\begin{aligned}
u(r+0,r) &= \frac{1}{4 \pi r}
\end{aligned}
$$

半径 $0<r<3$ に対して $t=3$ の値と、到着後の時間微分を求める。

$$
\begin{aligned}
t_{\rm emit}\big|_{t=3} &= 3 - r \\
\frac{\partial^1 (t_{\rm emit})}{\partial t^1} &= 1 \\
\left.\frac{\partial^1 (t_{\rm emit})}{\partial t^1}\right|_{t=3} &= 1
\end{aligned}
$$

半径 $0<r<3$ に対して $t=3$ の値と、到着後の時間微分を求める。

$$
\begin{aligned}
u(t,r)\big|_{t=3} &= \frac{\cos{\left(2 r - 6 \right)}}{4 \pi r} \\
\frac{\partial^1 (u(t,r))}{\partial t^1} &= \frac{\sin{\left(2 r - 2 t \right)}}{2 \pi r} \\
\left.\frac{\partial^1 (u(t,r))}{\partial t^1}\right|_{t=3} &= \frac{\sin{\left(2 r - 6 \right)}}{2 \pi r}
\end{aligned}
$$
:::

:::example[調和源を時刻零で入れる（3）]
波速を $1$、原点の点源を $H(t)\cos(mt)\delta^{(3)}(\vb x)$ とする。本文の波動演算子で $t>r$ の場を求める。

$$
m=3
$$

遅延時刻を計算する。

$$
\begin{aligned}
t_{\rm emit} &= - r + t
\end{aligned}
$$

デルタで時間源を評価する。

$$
\begin{aligned}
u(t,r)\big|_{m=3} &= \frac{\cos{\left(3 r - 3 t \right)}}{4 \pi r}
\end{aligned}
$$

到達直後の初期値には球面の広がりの因子が残る。

$$
\begin{aligned}
u(r+0,r) &= \frac{1}{4 \pi r}
\end{aligned}
$$

半径 $0<r<3$ に対して $t=3$ の値と、到着後の時間微分を求める。

$$
\begin{aligned}
t_{\rm emit}\big|_{t=3} &= 3 - r \\
\frac{\partial^1 (t_{\rm emit})}{\partial t^1} &= 1 \\
\left.\frac{\partial^1 (t_{\rm emit})}{\partial t^1}\right|_{t=3} &= 1
\end{aligned}
$$

半径 $0<r<3$ に対して $t=3$ の値と、到着後の時間微分を求める。

$$
\begin{aligned}
u(t,r)\big|_{t=3} &= \frac{\cos{\left(3 r - 9 \right)}}{4 \pi r} \\
\frac{\partial^1 (u(t,r))}{\partial t^1} &= \frac{3 \sin{\left(3 r - 3 t \right)}}{4 \pi r} \\
\left.\frac{\partial^1 (u(t,r))}{\partial t^1}\right|_{t=3} &= \frac{3 \sin{\left(3 r - 9 \right)}}{4 \pi r}
\end{aligned}
$$
:::

:::example[調和源を時刻零で入れる（4）]
波速を $1$、原点の点源を $H(t)\cos(mt)\delta^{(3)}(\vb x)$ とする。本文の波動演算子で $t>r$ の場を求める。

$$
m=4
$$

遅延時刻を計算する。

$$
\begin{aligned}
t_{\rm emit} &= - r + t
\end{aligned}
$$

デルタで時間源を評価する。

$$
\begin{aligned}
u(t,r)\big|_{m=4} &= \frac{\cos{\left(4 r - 4 t \right)}}{4 \pi r}
\end{aligned}
$$

到達直後の初期値には球面の広がりの因子が残る。

$$
\begin{aligned}
u(r+0,r) &= \frac{1}{4 \pi r}
\end{aligned}
$$

半径 $0<r<3$ に対して $t=3$ の値と、到着後の時間微分を求める。

$$
\begin{aligned}
t_{\rm emit}\big|_{t=3} &= 3 - r \\
\frac{\partial^1 (t_{\rm emit})}{\partial t^1} &= 1 \\
\left.\frac{\partial^1 (t_{\rm emit})}{\partial t^1}\right|_{t=3} &= 1
\end{aligned}
$$

半径 $0<r<3$ に対して $t=3$ の値と、到着後の時間微分を求める。

$$
\begin{aligned}
u(t,r)\big|_{t=3} &= \frac{\cos{\left(4 r - 12 \right)}}{4 \pi r} \\
\frac{\partial^1 (u(t,r))}{\partial t^1} &= \frac{\sin{\left(4 r - 4 t \right)}}{\pi r} \\
\left.\frac{\partial^1 (u(t,r))}{\partial t^1}\right|_{t=3} &= \frac{\sin{\left(4 r - 12 \right)}}{\pi r}
\end{aligned}
$$
:::

## 拡散方程式の Green 関数

### 拡散核の質量を保つ

$\partial_t-\kappa\Delta$ の遅延基本解は
$$
G(t,x)=H(t)(4\pi\kappa t)^{-d/2} e^{-|x|^2/(4\kappa t)}
\label{eq-mp-diffusion-green}
$$
である。空間変換後は $H(t)e^{-\kappa t|k|^2}$ となり、時間微分で $\delta(t)$ が一つ現れる。逆変換の Gaussian 積分は各方向に分けられ、時刻 $t>0$ で全質量は $1$ である。空間次元 $d$ で核の単位は $\mathrm{m^{-d}}$、一方向の分散は $2\kappa t$ である。半直線の零 Dirichlet では $G(t,x-y)-G(t,x+y)$ を使う。境界で二項が打ち消し、$t\downarrow0$ には $y>0$ の点源だけが内部に残る。零 Neumann では和になり半直線内の質量が保存される。

### 反射した熱核と質量

半直線で零 Dirichlet 条件なら全直線核の奇延長を取り、零 Neumann 条件なら偶延長を取る。半直線だけで積分すると、正の像では移した二つの半区間が合わさって全直線積分になる。負の像では境界により失われた分が引かれる。これは熱方程式を積分したときに境界から熱流が出ることと整合する。

:::example[有限区間の熱核を像と固有関数で表す]
長さ $a$ の区間の Dirichlet 熱核を考える。全直線の核を $K(t,z)$ とすると、奇延長を周期 $2a$ で繰り返すことで
$$
K_D(t,x,y)=\sum_{m\in\Z}
\big[K(t,x-y+2ma)-K(t,x+y+2ma)\big]
$$
を得る。$x=0,a$ では添字の付け替えと偶性で二つの和が相殺する。$t>0$ で Gaussian 減衰により項別微分でき、区間の内部で $t\downarrow0$ に残る点源は $y$ だけである。周期和 $\sum_m K(t,z+2ma)$ の第 $n$ Fourier 係数は、一周期の積分を全実線へ展開することで $(2a)^{-1}e^{-\kappa(n\pi/a)^2t}$ となる。$z=x-y$ と $z=x+y$ の差を取ると
$$
K_D(t,x,y)=\frac2a\sum_{n\ge1}e^{-\kappa(n\pi/a)^2t}
\sin\frac{n\pi x}{a}\sin\frac{n\pi y}{a}
$$
を得る。短時間では像の和、長時間では最低モードが残る後者の表示から挙動を読みやすい。Neumann 条件では二つの周期和を加え、余弦の積と定数項 $1/a$ が残る。
:::

:::example[吸収境界と熱量の減少（1）]
拡散係数を $1$、初期点源を $y=m$ とする。

$$
m=1
$$

全直線の核で正の位置にある源を観測する。

$$
\begin{aligned}
K(t,x-m)\big|_{m=1} &= \frac{e^{- \frac{\left(x - 1\right)^{2}}{4 u}}}{2 \sqrt{\pi} \sqrt{u}}
\end{aligned}
$$

反対符号の像を足す。

$$
\begin{aligned}
K_D(t,x,m)\big|_{m=1} &= \frac{\left(e^{\frac{x}{u}} - 1\right) e^{- \frac{x^{2} + 2 x + 1}{4 u}}}{2 \sqrt{\pi} \sqrt{u}}
\end{aligned}
$$

半直線で積分すると誤差関数が出る。

$$
\begin{aligned}
\int_0^\infty K_D(t,x,m)\dd x\big|_{m=1} &= \operatorname{erf}{\left(\frac{1}{2 \sqrt{u}} \right)}
\end{aligned}
$$

観測点 $u=1$ での値と導関数を計算する。

$$
\begin{aligned}
K(t,x-m)\big|_{u=1} &= \frac{e^{- \frac{\left(x - 1\right)^{2}}{4}}}{2 \sqrt{\pi}} \\
\frac{\partial^1 (K(t,x-m))}{\partial u^1} &= \frac{\left(- 2 u + \left(x - 1\right)^{2}\right) e^{- \frac{\left(x - 1\right)^{2}}{4 u}}}{8 \sqrt{\pi} u^{\frac{5}{2}}} \\
\left.\frac{\partial^1 (K(t,x-m))}{\partial u^1}\right|_{u=1} &= \frac{\left(\left(x - 1\right)^{2} - 2\right) e^{- \frac{\left(x - 1\right)^{2}}{4}}}{8 \sqrt{\pi}}
\end{aligned}
$$

観測点 $u=1$ での値と導関数を計算する。

$$
\begin{aligned}
K_D(t,x,m)\big|_{u=1} &= \frac{\left(e^{x} - 1\right) e^{- \frac{x^{2}}{4} - \frac{x}{2} - \frac{1}{4}}}{2 \sqrt{\pi}} \\
\frac{\partial^1 (K_D(t,x,m))}{\partial u^1} &= \frac{\left(- 2 u e^{\frac{x}{u}} + 2 u + x^{2} e^{\frac{x}{u}} - x^{2} - 2 x e^{\frac{x}{u}} - 2 x + e^{\frac{x}{u}} - 1\right) e^{- \frac{x^{2} + 2 x + 1}{4 u}}}{8 \sqrt{\pi} u^{\frac{5}{2}}} \\
\left.\frac{\partial^1 (K_D(t,x,m))}{\partial u^1}\right|_{u=1} &= \frac{\left(x^{2} e^{x} - x^{2} - 2 x e^{x} - 2 x - e^{x} + 1\right) e^{- \frac{x^{2}}{4} - \frac{x}{2} - \frac{1}{4}}}{8 \sqrt{\pi}}
\end{aligned}
$$

観測点 $u=1$ での値と導関数を計算する。

$$
\begin{aligned}
\int_0^\infty K_D(t,x,m)\dd x\big|_{u=1} &= \operatorname{erf}{\left(\frac{1}{2} \right)} \\
\frac{\partial^1 (\int_0^\infty K_D(t,x,m)\dd x)}{\partial u^1} &= - \frac{e^{- \frac{1}{4 u}}}{2 \sqrt{\pi} u^{\frac{3}{2}}} \\
\left.\frac{\partial^1 (\int_0^\infty K_D(t,x,m)\dd x)}{\partial u^1}\right|_{u=1} &= - \frac{1}{2 \sqrt{\pi} e^{\frac{1}{4}}}
\end{aligned}
$$
:::

:::example[吸収境界と熱量の減少（2）]
拡散係数を $1$、初期点源を $y=m$ とする。

$$
m=2
$$

全直線の核で正の位置にある源を観測する。

$$
\begin{aligned}
K(t,x-m)\big|_{m=2} &= \frac{e^{- \frac{\left(x - 2\right)^{2}}{4 u}}}{2 \sqrt{\pi} \sqrt{u}}
\end{aligned}
$$

反対符号の像を足す。

$$
\begin{aligned}
K_D(t,x,m)\big|_{m=2} &= \frac{\left(e^{\frac{2 x}{u}} - 1\right) e^{- \frac{\frac{x^{2}}{4} + x + 1}{u}}}{2 \sqrt{\pi} \sqrt{u}}
\end{aligned}
$$

半直線で積分すると誤差関数が出る。

$$
\begin{aligned}
\int_0^\infty K_D(t,x,m)\dd x\big|_{m=2} &= \operatorname{erf}{\left(\frac{1}{\sqrt{u}} \right)}
\end{aligned}
$$

観測点 $u=1$ での値と導関数を計算する。

$$
\begin{aligned}
K(t,x-m)\big|_{u=1} &= \frac{e^{- \frac{\left(x - 2\right)^{2}}{4}}}{2 \sqrt{\pi}} \\
\frac{\partial^1 (K(t,x-m))}{\partial u^1} &= \frac{\left(- 2 u + \left(x - 2\right)^{2}\right) e^{- \frac{\left(x - 2\right)^{2}}{4 u}}}{8 \sqrt{\pi} u^{\frac{5}{2}}} \\
\left.\frac{\partial^1 (K(t,x-m))}{\partial u^1}\right|_{u=1} &= \frac{\left(\left(x - 2\right)^{2} - 2\right) e^{- \frac{\left(x - 2\right)^{2}}{4}}}{8 \sqrt{\pi}}
\end{aligned}
$$

観測点 $u=1$ での値と導関数を計算する。

$$
\begin{aligned}
K_D(t,x,m)\big|_{u=1} &= \frac{\left(e^{2 x} - 1\right) e^{- \frac{x^{2}}{4} - x - 1}}{2 \sqrt{\pi}} \\
\frac{\partial^1 (K_D(t,x,m))}{\partial u^1} &= \frac{\left(- 2 u e^{\frac{2 x}{u}} + 2 u + x^{2} e^{\frac{2 x}{u}} - x^{2} - 4 x e^{\frac{2 x}{u}} - 4 x + 4 e^{\frac{2 x}{u}} - 4\right) e^{- \frac{\frac{x^{2}}{4} + x + 1}{u}}}{8 \sqrt{\pi} u^{\frac{5}{2}}} \\
\left.\frac{\partial^1 (K_D(t,x,m))}{\partial u^1}\right|_{u=1} &= \frac{\left(x^{2} e^{2 x} - x^{2} - 4 x e^{2 x} - 4 x + 2 e^{2 x} - 2\right) e^{- \frac{x^{2}}{4} - x - 1}}{8 \sqrt{\pi}}
\end{aligned}
$$

観測点 $u=1$ での値と導関数を計算する。

$$
\begin{aligned}
\int_0^\infty K_D(t,x,m)\dd x\big|_{u=1} &= \operatorname{erf}{\left(1 \right)} \\
\frac{\partial^1 (\int_0^\infty K_D(t,x,m)\dd x)}{\partial u^1} &= - \frac{e^{- \frac{1}{u}}}{\sqrt{\pi} u^{\frac{3}{2}}} \\
\left.\frac{\partial^1 (\int_0^\infty K_D(t,x,m)\dd x)}{\partial u^1}\right|_{u=1} &= - \frac{1}{e \sqrt{\pi}}
\end{aligned}
$$
:::

:::example[吸収境界と熱量の減少（3）]
拡散係数を $1$、初期点源を $y=m$ とする。

$$
m=3
$$

全直線の核で正の位置にある源を観測する。

$$
\begin{aligned}
K(t,x-m)\big|_{m=3} &= \frac{e^{- \frac{\left(x - 3\right)^{2}}{4 u}}}{2 \sqrt{\pi} \sqrt{u}}
\end{aligned}
$$

反対符号の像を足す。

$$
\begin{aligned}
K_D(t,x,m)\big|_{m=3} &= \frac{\left(e^{\frac{3 x}{u}} - 1\right) e^{- \frac{x^{2} + 6 x + 9}{4 u}}}{2 \sqrt{\pi} \sqrt{u}}
\end{aligned}
$$

半直線で積分すると誤差関数が出る。

$$
\begin{aligned}
\int_0^\infty K_D(t,x,m)\dd x\big|_{m=3} &= \operatorname{erf}{\left(\frac{3}{2 \sqrt{u}} \right)}
\end{aligned}
$$

観測点 $u=1$ での値と導関数を計算する。

$$
\begin{aligned}
K(t,x-m)\big|_{u=1} &= \frac{e^{- \frac{\left(x - 3\right)^{2}}{4}}}{2 \sqrt{\pi}} \\
\frac{\partial^1 (K(t,x-m))}{\partial u^1} &= \frac{\left(- 2 u + \left(x - 3\right)^{2}\right) e^{- \frac{\left(x - 3\right)^{2}}{4 u}}}{8 \sqrt{\pi} u^{\frac{5}{2}}} \\
\left.\frac{\partial^1 (K(t,x-m))}{\partial u^1}\right|_{u=1} &= \frac{\left(\left(x - 3\right)^{2} - 2\right) e^{- \frac{\left(x - 3\right)^{2}}{4}}}{8 \sqrt{\pi}}
\end{aligned}
$$

観測点 $u=1$ での値と導関数を計算する。

$$
\begin{aligned}
K_D(t,x,m)\big|_{u=1} &= \frac{\left(e^{3 x} - 1\right) e^{- \frac{x^{2}}{4} - \frac{3 x}{2} - \frac{9}{4}}}{2 \sqrt{\pi}} \\
\frac{\partial^1 (K_D(t,x,m))}{\partial u^1} &= \frac{\left(- 2 u e^{\frac{3 x}{u}} + 2 u + x^{2} e^{\frac{3 x}{u}} - x^{2} - 6 x e^{\frac{3 x}{u}} - 6 x + 9 e^{\frac{3 x}{u}} - 9\right) e^{- \frac{x^{2} + 6 x + 9}{4 u}}}{8 \sqrt{\pi} u^{\frac{5}{2}}} \\
\left.\frac{\partial^1 (K_D(t,x,m))}{\partial u^1}\right|_{u=1} &= \frac{\left(x^{2} e^{3 x} - x^{2} - 6 x e^{3 x} - 6 x + 7 e^{3 x} - 7\right) e^{- \frac{x^{2}}{4} - \frac{3 x}{2} - \frac{9}{4}}}{8 \sqrt{\pi}}
\end{aligned}
$$

観測点 $u=1$ での値と導関数を計算する。

$$
\begin{aligned}
\int_0^\infty K_D(t,x,m)\dd x\big|_{u=1} &= \operatorname{erf}{\left(\frac{3}{2} \right)} \\
\frac{\partial^1 (\int_0^\infty K_D(t,x,m)\dd x)}{\partial u^1} &= - \frac{3 e^{- \frac{9}{4 u}}}{2 \sqrt{\pi} u^{\frac{3}{2}}} \\
\left.\frac{\partial^1 (\int_0^\infty K_D(t,x,m)\dd x)}{\partial u^1}\right|_{u=1} &= - \frac{3}{2 \sqrt{\pi} e^{\frac{9}{4}}}
\end{aligned}
$$
:::

:::example[吸収境界と熱量の減少（4）]
拡散係数を $1$、初期点源を $y=m$ とする。

$$
m=4
$$

全直線の核で正の位置にある源を観測する。

$$
\begin{aligned}
K(t,x-m)\big|_{m=4} &= \frac{e^{- \frac{\left(x - 4\right)^{2}}{4 u}}}{2 \sqrt{\pi} \sqrt{u}}
\end{aligned}
$$

反対符号の像を足す。

$$
\begin{aligned}
K_D(t,x,m)\big|_{m=4} &= \frac{\left(e^{\frac{4 x}{u}} - 1\right) e^{- \frac{\frac{x^{2}}{4} + 2 x + 4}{u}}}{2 \sqrt{\pi} \sqrt{u}}
\end{aligned}
$$

半直線で積分すると誤差関数が出る。

$$
\begin{aligned}
\int_0^\infty K_D(t,x,m)\dd x\big|_{m=4} &= \operatorname{erf}{\left(\frac{2}{\sqrt{u}} \right)}
\end{aligned}
$$

観測点 $u=1$ での値と導関数を計算する。

$$
\begin{aligned}
K(t,x-m)\big|_{u=1} &= \frac{e^{- \frac{\left(x - 4\right)^{2}}{4}}}{2 \sqrt{\pi}} \\
\frac{\partial^1 (K(t,x-m))}{\partial u^1} &= \frac{\left(- 2 u + \left(x - 4\right)^{2}\right) e^{- \frac{\left(x - 4\right)^{2}}{4 u}}}{8 \sqrt{\pi} u^{\frac{5}{2}}} \\
\left.\frac{\partial^1 (K(t,x-m))}{\partial u^1}\right|_{u=1} &= \frac{\left(\left(x - 4\right)^{2} - 2\right) e^{- \frac{\left(x - 4\right)^{2}}{4}}}{8 \sqrt{\pi}}
\end{aligned}
$$

観測点 $u=1$ での値と導関数を計算する。

$$
\begin{aligned}
K_D(t,x,m)\big|_{u=1} &= \frac{\left(e^{4 x} - 1\right) e^{- \frac{x^{2}}{4} - 2 x - 4}}{2 \sqrt{\pi}} \\
\frac{\partial^1 (K_D(t,x,m))}{\partial u^1} &= \frac{\left(- \frac{u e^{\frac{4 x}{u}}}{4} + \frac{u}{4} + \frac{x^{2} e^{\frac{4 x}{u}}}{8} - \frac{x^{2}}{8} - x e^{\frac{4 x}{u}} - x + 2 e^{\frac{4 x}{u}} - 2\right) e^{- \frac{\frac{x^{2}}{4} + 2 x + 4}{u}}}{\sqrt{\pi} u^{\frac{5}{2}}} \\
\left.\frac{\partial^1 (K_D(t,x,m))}{\partial u^1}\right|_{u=1} &= \frac{\left(x^{2} e^{4 x} - x^{2} - 8 x e^{4 x} - 8 x + 14 e^{4 x} - 14\right) e^{- \frac{x^{2}}{4} - 2 x - 4}}{8 \sqrt{\pi}}
\end{aligned}
$$

観測点 $u=1$ での値と導関数を計算する。

$$
\begin{aligned}
\int_0^\infty K_D(t,x,m)\dd x\big|_{u=1} &= \operatorname{erf}{\left(2 \right)} \\
\frac{\partial^1 (\int_0^\infty K_D(t,x,m)\dd x)}{\partial u^1} &= - \frac{2 e^{- \frac{4}{u}}}{\sqrt{\pi} u^{\frac{3}{2}}} \\
\left.\frac{\partial^1 (\int_0^\infty K_D(t,x,m)\dd x)}{\partial u^1}\right|_{u=1} &= - \frac{2}{\sqrt{\pi} e^{4}}
\end{aligned}
$$
:::

## Green 関数による摂動展開と Dyson 方程式

### 摂動のたびに粒子または場が伝播する

閉作用素 $L_0$ と有界作用素 $V$ に対し、$L=L_0+V$ を $L_0$ と同じ定義域で考える。$R_0=(L_0-z)^{-1}$、$R=(L-z)^{-1}$ が存在するとする。$(L_0-z)R=I-VR$ の両辺に左から $R_0$ を掛けると
$$
R=R_0-R_0VR .
\label{eq-mp-dyson}
$$
これがこの分母の規約での Dyson 方程式である。$\|R_0V\|<1$ なら $(I+R_0V)^{-1}=\sum_{n\ge0}(-R_0V)^n$ であり、有限和を掛けた残りがノルムで零になるため
$R=R_0-R_0VR_0+R_0VR_0VR_0-\cdots$ が成立する。例えば $L_0=2I,V=\epsilon I,z=0$ では $R=(2+\epsilon)^{-1}I=(1/2-\epsilon/4+\epsilon^2/8+\cdots)I$、収束は $|\epsilon|<2$ である。外向き核で実スペクトル上へ極限を取る場合、ノルムの条件は用いた関数空間で確認する。周波数の符号や逆作用素を $(z-L)^{-1}$ で定義する流儀と混用すると、摂動項の符号が変わる。

### 階数一の摂動を正確に足す

$L=L_0+\alpha\ket{v}\bra{v}$ とする。ここでは Dirac 記法を用い、$\braket{v}{u}$ は $u$ に線形である。$G_0=L_0^{-1}$ がある場合、$LG=I$ を $G=G_0-\alpha G_0\ket{v}\bra{v}G$ と書く。左から $\bra{v}$ を掛けて $(1+\alpha\mel{v}{G_0}{v})\bra{v}G=\bra{v}G_0$ を得る。従って分母が零でないときは
$$
G=G_0-\frac{\alpha G_0\ket{v}\bra{v}G_0}{1+\alpha\mel{v}{G_0}{v}}.
$$
両側からの掛け算でこれが逆になることを確認できる。分母が零なら $L G_0 v=(1+\alpha\mel{v}{G_0}{v})v=0$ であり、$v\ne0$ では $G_0v$ が零固有ベクトルになる。この式は Neumann 級数が収束しない範囲でも、分母が零でなければ成り立つ。

### 外向き核を用いた散乱と Born 近似

$(-\Delta+V-k^2)u=0$ で、有界な $V$ の台が半径 $a$ の球の内部にあるとする。入射波 $u_{\rm in}(x)=e^{ik\widehat d\cdot x}$ を指定し、散乱部分には外向き条件を課す。すると $(-\Delta-k^2)(u-u_{\rm in})=-Vu$ なので、外向き核の畳み込みで
$$
u(x)=u_{\rm in}(x)-\int g_k(|x-y|)V(y)u(y)\,\dd y
$$
を得る。この表示を :term[リップマン・シュウィンガー方程式]{en="Lippmann–Schwinger equation" yomi="りっぷまんしゅうぃんがーほうていしき"} という。逆向きには $-\Delta-k^2$ を作用させれば微分方程式に戻り、積分項が外向き条件を満たす。$x$ を球の内部に限ると、積分作用素の一様ノルムは $\|V\|_\infty\int_{|z|\le2a}(4\pi|z|)^{-1}\dd z=2a^2\|V\|_\infty$ 以下である。この上界が $1$ 未満なら球内で反復は収束し、得られた値から球外での場も積分で決まる。

$|x|=r\to\infty$ で $\widehat x=x/r$ を固定すると、台の上で $|x-y|=r-\widehat x\cdot y+O(a^2/r)$ である。よって $u-u_{\rm in}\sim e^{ikr}f(\widehat x)/r$ の係数は
$$
f(\widehat x)=-\frac1{4\pi}\int e^{-ik\widehat x\cdot y}V(y)u(y)\,\dd y
$$
となる。ここで $u$ を $u_{\rm in}$ に置いたものが第一 Born 近似である。$V$ の次元は長さの逆二乗、$g_k$ の次元は長さの逆数であり、この積分で $f$ は長さの次元を持つ。$2a^2\|V\|_\infty<1$ のとき、核の反復で $u-u_{\rm in}$ は $V$ の一階以上だから、振幅で省いた項は $V$ の二階以上である。

:::example[球内で一定の弱い散乱ポテンシャル]
$|y|<a$ で $V=V_0$、外で零とする。$q=k|\widehat x-\widehat d|$ と置けば、角積分の平均は $\sin(qr)/(qr)$ である。従って第一 Born 近似の振幅は
$$
f_{\rm B}=-V_0\int_0^a r^2\frac{\sin(qr)}{qr}\,\dd r
=-V_0\frac{\sin(qa)-qa\cos(qa)}{q^3}.
$$
前方の $q=0$ では $\sin z-z\cos z=z^3/3+O(z^5)$ なので $f_{\rm B}=-V_0 a^3/3$ となる。これは $-1/(4\pi)$ 倍した $V$ の全空間積分と一致する。
:::

:::example[二成分の応答と一方向の摂動（1）]
非摂動の作用素を単位行列とし、方向を $(1,1)^{\mathsf T}$ とする。

$$
m=1
$$

摂動の階数一行列を求める。

$$
\begin{aligned}
L\big|_{m=1} &= \begin{pmatrix}2 & 1\\1 & 2\end{pmatrix}
\end{aligned}
$$

分母は $1+2m$ になる。

$$
\begin{aligned}
G\big|_{m=1} &= \begin{pmatrix}\frac{2}{3} & - \frac{1}{3}\\- \frac{1}{3} & \frac{2}{3}\end{pmatrix}
\end{aligned}
$$

摂動に平行な方向と直交する方向で固有値を比べる。

$$
\begin{aligned}
(\lambda_\parallel(G),\lambda_\perp(G))\big|_{m=1} &= \left( \frac{1}{3}, \  1\right)
\end{aligned}
$$
:::

:::example[二成分の応答と一方向の摂動（2）]
非摂動の作用素を単位行列とし、方向を $(1,1)^{\mathsf T}$ とする。

$$
m=2
$$

摂動の階数一行列を求める。

$$
\begin{aligned}
L\big|_{m=2} &= \begin{pmatrix}3 & 2\\2 & 3\end{pmatrix}
\end{aligned}
$$

分母は $1+2m$ になる。

$$
\begin{aligned}
G\big|_{m=2} &= \begin{pmatrix}\frac{3}{5} & - \frac{2}{5}\\- \frac{2}{5} & \frac{3}{5}\end{pmatrix}
\end{aligned}
$$

摂動に平行な方向と直交する方向で固有値を比べる。

$$
\begin{aligned}
(\lambda_\parallel(G),\lambda_\perp(G))\big|_{m=2} &= \left( \frac{1}{5}, \  1\right)
\end{aligned}
$$
:::

:::example[二成分の応答と一方向の摂動（3）]
非摂動の作用素を単位行列とし、方向を $(1,1)^{\mathsf T}$ とする。

$$
m=3
$$

摂動の階数一行列を求める。

$$
\begin{aligned}
L\big|_{m=3} &= \begin{pmatrix}4 & 3\\3 & 4\end{pmatrix}
\end{aligned}
$$

分母は $1+2m$ になる。

$$
\begin{aligned}
G\big|_{m=3} &= \begin{pmatrix}\frac{4}{7} & - \frac{3}{7}\\- \frac{3}{7} & \frac{4}{7}\end{pmatrix}
\end{aligned}
$$

摂動に平行な方向と直交する方向で固有値を比べる。

$$
\begin{aligned}
(\lambda_\parallel(G),\lambda_\perp(G))\big|_{m=3} &= \left( \frac{1}{7}, \  1\right)
\end{aligned}
$$
:::

:::example[二成分の応答と一方向の摂動（4）]
非摂動の作用素を単位行列とし、方向を $(1,1)^{\mathsf T}$ とする。

$$
m=4
$$

摂動の階数一行列を求める。

$$
\begin{aligned}
L\big|_{m=4} &= \begin{pmatrix}5 & 4\\4 & 5\end{pmatrix}
\end{aligned}
$$

分母は $1+2m$ になる。

$$
\begin{aligned}
G\big|_{m=4} &= \begin{pmatrix}\frac{5}{9} & - \frac{4}{9}\\- \frac{4}{9} & \frac{5}{9}\end{pmatrix}
\end{aligned}
$$

摂動に平行な方向と直交する方向で固有値を比べる。

$$
\begin{aligned}
(\lambda_\parallel(G),\lambda_\perp(G))\big|_{m=4} &= \left( \frac{1}{9}, \  1\right)
\end{aligned}
$$
:::

## 演習問題

::::exercise{level=basic}
単位区間の Dirichlet $-\partial_x^2$ の $G(1/4,3/4)$ を求めよ。
:::solution
核は $x_<(1-x_>)$ なので $(1/4)(1/4)=1/16$ である。
:::
::::

::::exercise{level=basic}
同じ区間で $f=2$ の解を求めよ。
:::solution
$-u''=2$ から $u=-x^2+Ax+B$、端点条件から $B=0,A=1$ で $x(1-x)$ となる。核を定数 $2$ に作用させても同じ値になる。
:::
::::

::::exercise{level=basic}
$g_k$ の $(\partial_r-ik)g_k$ を求めよ。
:::solution
微分して $ik e^{ikr}/(4\pi r)-e^{ikr}/(4\pi r^2)$ を得る。第一項が打ち消され、残りは $-e^{ikr}/(4\pi r^2)$ である。
:::
::::

::::exercise{level=basic}
一次元遅延波動核を $x$ で積分せよ。
:::solution
$t>0$ では長さ $2ct$ の区間に高さ $1/(2c)$ があるので積分は $t$ である。$t<0$ では零であり、全体で $tH(t)$ となる。
:::
::::

::::exercise{level=advanced}
三次元遅延波動核の空間積分を一次元の値と比較せよ。
:::solution
球面積を掛けると $\int_0^\infty r\delta(t-r/c)\dd r/c^2$ になる。$t>0$ では $r=ct$ と Jacobian $c$ を得るので $t$、$t<0$ では零である。従って $tH(t)$ で一致する。
:::
::::

::::exercise{level=advanced}
球の中心に点源があるときの法線微分を $r=a$ で求めよ。
:::solution
$G_D=1/(4\pi r)-1/(4\pi a)$ で外向きは $+\partial_r$ だから $\partial_nG_D=-1/(4\pi a^2)$ である。その球面積分は $-1$ となり、$-\Delta G_D$ の質量 $1$ と整合する。
:::
::::

::::exercise{level=advanced}
半直線の Neumann 熱核の質量が $1$ であることを示せ。
:::solution
二項の $x$ 積分を変数変換すると $\int_{-y}^\infty G(t,u)\dd u+\int_y^\infty G(t,u)\dd u$ となる。偶性により第 2 項は $\int_{-\infty}^{-y}$ と等しいので全実線の質量 $1$ である。
:::
::::

::::exercise{level=advanced}
$R_0=\diag(1,1/2)$、$V=\epsilon\begin{pmatrix}0&1\\1&0\end{pmatrix}$ の逆の一次の非対角成分を求めよ。
:::solution
一次項 $-R_0VR_0$ の非対角は $-\epsilon/2$ である。直接逆行列を作ると非対角は $-\epsilon/(2-\epsilon^2)$ となり、一次項が一致する。
:::
::::

## 確認テスト

::::quiz{answer="1" tol="1e-9"}
$t>0$ の全空間熱核の質量を答えよ。

:::explanation
各座標で Gaussian 積分が分母の $\sqrt{4\pi\kappa t}$ を打ち消すので、全積は $1$ である。
:::
::::

::::quiz
$L=-\partial_x^2$ の Green 核の $x=y$ での導関数の跳びはどれか。

- [ ] $+1$
- [x] $-1$
- [ ] $0$

:::explanation
連続な核の二階微分には一階微分の跳びがデルタ係数として出る。作用素の負符号を掛けて $+1$ になるには跳びが $-1$ である。
:::
::::

::::quiz{answer="0.0625" tol="1e-9"}
単位区間の Dirichlet 核 $G(1/4,3/4)$ を答えよ。

:::explanation
$x_<=1/4$、$1-x_>=1/4$ なので積は $1/16=0.0625$ である。
:::
::::

::::quiz
遅延 Green 関数の負時刻での値はどれか。

- [x] 零
- [ ] 常に先進核と等しい
- [ ] 一様な非零定数

:::explanation
原点の時刻に源を置いたとき、その源より前の応答を零とするのが遅延条件である。先進核では正時刻側が零になる。
:::
::::

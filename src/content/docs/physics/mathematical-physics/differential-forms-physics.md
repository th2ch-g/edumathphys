---
title: "微分形式と物理"
description: "微分形式の演算と積分を導き、電磁場および熱力学の局所構造を記述する。"
notation:
  - symbol: '\Omega^p(M)'
    meaning: "滑らかな p 形式の空間"
  - symbol: '\wedge'
    meaning: "外積"
  - symbol: '\mathrm d'
    meaning: "外微分"
  - symbol: '\star'
    meaning: "Hodge 双対"
  - symbol: '\delta'
    meaning: "余微分"
---

## この章の目標

- 外積、引き戻し、外微分を座標で計算する。
- Hodge 双対と Stokes の定理を使い、ベクトル解析の公式を導く。
- Maxwell 方程式を電磁場と電流の微分形式で表す。
- 熱の Pfaff 形式の可積分条件を調べる。

:ref[def-mfd-differential-form]{text="微分形式"} を流束や状態空間の積分に適用する。外積代数を使うと、線に沿う循環と面を通る流束をそれぞれ積分で表せる。

## 外積代数と微分形式

### 向き付き面に値を与える

流速は一本の方向への成分で表せるが、流束を知るには面の向きと広さが必要である。二本の接ベクトルを入れ替えると符号が反転する量を使えば、同じ面を別の座標で記述しても値を保てる。

:::definition[微分形式]{#def-mp-differential-form}
滑らかな多様体 $M$ 上で、各点の接空間の $p$ 個のベクトルに対して交代的な多重線形写像を滑らかに割り当てたものを :term[微分形式]{en="differential form" yomi="びぶんけいしき"} という。$p=0$ の形式は関数である。$p$ 形式と $q$ 形式の :term[外積]{en="wedge product" yomi="がいせき"} は、$\omega\otimes\eta$ を全置換で交代化し、$p!q!$ で割った $(p+q)$ 形式である。ここでは交代化の和に $1/(p+q)!$ を含めない。
:::

座標 $x^i$ では $\dd x^{i_1}\wedge\cdots\wedge\dd x^{i_p}$ を $i_1<\cdots<i_p$ に限ったものが基底となる。評価すると、各 $\dd x^{i_a}$ の入力値を並べた行列の行列式になる。従って基底の交換則は一対ごとの符号反転であり、$\omega\wedge\eta=(-1)^{pq}\eta\wedge\omega$ となる。三つの単項式の外積は括弧の位置にかかわらず同じ順の列を連結した形なので結合的である。一般の形式では双線形性でこの性質が延びる。外積代数の基底に関しては :ref[thm-linalg-exterior-basis] も参照せよ。

写像 $\Phi:N\to M$ の :term[引き戻し]{en="pullback" yomi="ひきもどし"} $\Phi^*\omega$ は、各入力に微分 $\dd\Phi$ を作用させてから $\omega$ を評価する。外積の定義に代入すると $\Phi^*(\omega\wedge\eta)=\Phi^*\omega\wedge\Phi^*\eta$ となる。

:::example[斜めの面積]
$\Phi(u,v)=(u,v,2u+3v)$ とすると $\Phi^*(\dd y\wedge\dd z)=\dd v\wedge(2\dd u+3\dd v)=-2\dd u\wedge\dd v$ である。単位正方形で積分した値は $-2$ で、$yz$ 平面への符号付き射影面積を表す。
:::

:::example[円板の面積形式]
$x=r\cos\theta$、$y=r\sin\theta$ の微分を外積に代入すると $\dd x\wedge\dd y=r\dd r\wedge\dd\theta$ である。従って半径 $2$ の円板の面積は $\int_0^{2\pi}\int_0^2 r\,\dd r\dd\theta=4\pi$ となる。写像が座標にならない $r=0$ は積分上の測度零の集合である。
:::

### 差し込みによる面の向き

ベクトル $X$ の $p$ 形式への差し込み $\iota_X\omega$ は $(\iota_X\omega)(v_1,\ldots,v_{p-1})=\omega(X,v_1,\ldots,v_{p-1})$ とする。零形式への差し込みは零である。一形式 $\alpha_1,\ldots,\alpha_p$ の外積に差し込むと、行列式の第 1 列を展開することで
$$
\iota_X(\alpha_1\wedge\cdots\wedge\alpha_p)
=\sum_{j=1}^p(-1)^{j-1}\alpha_j(X)\alpha_1\wedge\cdots\widehat{\alpha_j}\cdots\wedge\alpha_p
$$
を得る。外す位置を $j$ とすると、それより前を通過する回数が $j-1$ だからである。従って $\iota_X(\alpha\wedge\beta)=\iota_X\alpha\wedge\beta+(-1)^{\deg\alpha}\alpha\wedge\iota_X\beta$ である。例えば $X=2\partial_x-\partial_y$ なら $\iota_X(\dd x\wedge\dd y\wedge\dd z)=2\dd y\wedge\dd z+\dd x\wedge\dd z$ である。後者の正号は $\partial_y$ の係数と列の符号の両方が負であるためである。

### 行列式から引き戻しを計算する

一形式を座標写像 $F$ で引き戻すと、各 $\dd x^i$ を $\sum_j\partial_jF^i\,\dd u^j$ で置き換える。二形式では二つの置換を外積に入れ、同じ $\dd u^j$ が並ぶ項を消す。交換した項に負号を付けて集めると、Jacobian の小行列式になる。二次元の最高次形式ではこの小行列式が全行列式である。同じ変数変換を向きを無視する面積積分に使うときは絶対値が要る。

:::example[平面の曲線座標での面積形式（1）]
写像 $(x,y)=(u,v+u^m)$ で一形式と二形式を引き戻す。基準の値で計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

二番目の座標の微分に、$u$ 方向の寄与が加わる。

$$
\begin{aligned}
\partial_u y\big|_{m=1} &= 1
\end{aligned}
$$

元の一形式を $y\,\dd x+x\,\dd y$ とする。$\dd u$ の係数をまとめる。

$$
\begin{aligned}
[F^*(y\,\dd x+x\,\dd y)]_u\big|_{m=1} &= 2 u + v
\end{aligned}
$$

二形式の $\dd u\wedge\dd v$ の係数には小行列式が入る。ここでは座標のせん断が行列式を変えない。

$$
\begin{aligned}
[F^*(\dd x\wedge\dd y)]_{uv} &= 1
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
[F^*(y\,\dd x+x\,\dd y)]_u\big|_{u=2} &= v + 4 \\
\frac{\partial^1 ([F^*(y\,\dd x+x\,\dd y)]_u)}{\partial u^1} &= 2 \\
\left.\frac{\partial^1 ([F^*(y\,\dd x+x\,\dd y)]_u)}{\partial u^1}\right|_{u=2} &= 2
\end{aligned}
$$
:::

:::example[平面の曲線座標での面積形式（2）]
写像 $(x,y)=(u,v+u^m)$ で一形式と二形式を引き戻す。値を変えて依存性を調べる。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

二番目の座標の微分に、$u$ 方向の寄与が加わる。

$$
\begin{aligned}
\partial_u y\big|_{m=2} &= 2 u
\end{aligned}
$$

元の一形式を $y\,\dd x+x\,\dd y$ とする。$\dd u$ の係数をまとめる。

$$
\begin{aligned}
[F^*(y\,\dd x+x\,\dd y)]_u\big|_{m=2} &= 3 u^{2} + v
\end{aligned}
$$

二形式の $\dd u\wedge\dd v$ の係数には小行列式が入る。ここでは座標のせん断が行列式を変えない。

$$
\begin{aligned}
[F^*(\dd x\wedge\dd y)]_{uv} &= 1
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\partial_u y\big|_{u=2} &= 4 \\
\frac{\partial^1 (\partial_u y)}{\partial u^1} &= 2 \\
\left.\frac{\partial^1 (\partial_u y)}{\partial u^1}\right|_{u=2} &= 2
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
[F^*(y\,\dd x+x\,\dd y)]_u\big|_{u=2} &= v + 12 \\
\frac{\partial^1 ([F^*(y\,\dd x+x\,\dd y)]_u)}{\partial u^1} &= 6 u \\
\left.\frac{\partial^1 ([F^*(y\,\dd x+x\,\dd y)]_u)}{\partial u^1}\right|_{u=2} &= 12
\end{aligned}
$$
:::

:::example[平面の曲線座標での面積形式（3）]
写像 $(x,y)=(u,v+u^m)$ で一形式と二形式を引き戻す。さらに値を変えて比較する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

二番目の座標の微分に、$u$ 方向の寄与が加わる。

$$
\begin{aligned}
\partial_u y\big|_{m=3} &= 3 u^{2}
\end{aligned}
$$

元の一形式を $y\,\dd x+x\,\dd y$ とする。$\dd u$ の係数をまとめる。

$$
\begin{aligned}
[F^*(y\,\dd x+x\,\dd y)]_u\big|_{m=3} &= 4 u^{3} + v
\end{aligned}
$$

二形式の $\dd u\wedge\dd v$ の係数には小行列式が入る。ここでは座標のせん断が行列式を変えない。

$$
\begin{aligned}
[F^*(\dd x\wedge\dd y)]_{uv} &= 1
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\partial_u y\big|_{u=2} &= 12 \\
\frac{\partial^1 (\partial_u y)}{\partial u^1} &= 6 u \\
\left.\frac{\partial^1 (\partial_u y)}{\partial u^1}\right|_{u=2} &= 12
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
[F^*(y\,\dd x+x\,\dd y)]_u\big|_{u=2} &= v + 32 \\
\frac{\partial^1 ([F^*(y\,\dd x+x\,\dd y)]_u)}{\partial u^1} &= 12 u^{2} \\
\left.\frac{\partial^1 ([F^*(y\,\dd x+x\,\dd y)]_u)}{\partial u^1}\right|_{u=2} &= 48
\end{aligned}
$$
:::

:::example[平面の曲線座標での面積形式（4）]
写像 $(x,y)=(u,v+u^m)$ で一形式と二形式を引き戻す。最後にこの値でも計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=4
$$

二番目の座標の微分に、$u$ 方向の寄与が加わる。

$$
\begin{aligned}
\partial_u y\big|_{m=4} &= 4 u^{3}
\end{aligned}
$$

元の一形式を $y\,\dd x+x\,\dd y$ とする。$\dd u$ の係数をまとめる。

$$
\begin{aligned}
[F^*(y\,\dd x+x\,\dd y)]_u\big|_{m=4} &= 5 u^{4} + v
\end{aligned}
$$

二形式の $\dd u\wedge\dd v$ の係数には小行列式が入る。ここでは座標のせん断が行列式を変えない。

$$
\begin{aligned}
[F^*(\dd x\wedge\dd y)]_{uv} &= 1
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\partial_u y\big|_{u=2} &= 32 \\
\frac{\partial^1 (\partial_u y)}{\partial u^1} &= 12 u^{2} \\
\left.\frac{\partial^1 (\partial_u y)}{\partial u^1}\right|_{u=2} &= 48
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
[F^*(y\,\dd x+x\,\dd y)]_u\big|_{u=2} &= v + 80 \\
\frac{\partial^1 ([F^*(y\,\dd x+x\,\dd y)]_u)}{\partial u^1} &= 20 u^{3} \\
\left.\frac{\partial^1 ([F^*(y\,\dd x+x\,\dd y)]_u)}{\partial u^1}\right|_{u=2} &= 160
\end{aligned}
$$
:::

## 外微分と Poincaré の補題

### 外微分と局所的なポテンシャル

閉曲線での循環はその内側の微小な循環を合計したものである。この局所的な循環を作る演算は、一変数の微分を高次の向き付き面へ拡張する。

:::definition[外微分]{#def-mp-exterior-derivative}
$\omega=\sum_I\omega_I\dd x^I$ と書き、:term[外微分]{en="exterior derivative" yomi="がいびぶん"} を $\dd\omega=\sum_I \dd\omega_I\wedge\dd x^I$ とする。ここで $\dd f=\sum_i(\partial_i f)\dd x^i$ である。$\dd\omega=0$ の形式を閉形式、$\omega=\dd\beta$ と書ける形式を完全形式という。
:::

座標式に積の微分を使うと $\dd(\alpha\wedge\beta)=\dd\alpha\wedge\beta+(-1)^{\deg\alpha}\alpha\wedge\dd\beta$ である。また $\dd^2=0$ は二階微分の対称性と外積の反対称性から従う。引き戻しは関数の微分と連鎖律で交換し、$\dd^2\Phi^i=0$ だから座標の一形式でも交換する。積の微分により任意の形式で $\Phi^*\dd=\dd\Phi^*$ が成り立ち、外微分の定義は座標に依存しない。

:::theorem[星形領域での Poincaré の補題]{#thm-mp-poincare}
原点に関して星形な開集合 $U\subset\R^n$ で $p>0$ とする。閉 $p$ 形式 $\omega$ は完全である。
:::

:::proof
$H:[0,1]\times U\to U$ を $H(t,x)=tx$ とする。$H^*\omega=\alpha(t)+\dd t\wedge\beta(t)$ と分ける。$x$ 方向だけの外微分を $\dd_x$ とすれば、$\dd H^*\omega$ の $\dd t$ 成分は $\partial_t\alpha-\dd_x\beta$ となる。閉性によりこれは零である。積分すると $\alpha(1)-\alpha(0)=\dd_x\int_0^1\beta(t)\dd t$ であり、$\alpha(1)=\omega$、$\alpha(0)=0$ だから結論を得る。具体的に $\beta(t)_x(v_1,\ldots,v_{p-1})=t^{p-1}\omega_{tx}(x,v_1,\ldots,v_{p-1})$ である。積分は $t=0$ まで滑らかである。
:::

:::example[二形式の原始形式]
$\omega=3x\,\dd x\wedge\dd y$ では $\beta=3t^2x(x\dd y-y\dd x)$ なので原始形式は $x(x\dd y-y\dd x)$ になる。微分すると $2x\dd x\wedge\dd y-x\dd y\wedge\dd x=3x\dd x\wedge\dd y$ と確認できる。
:::

:::example[穴の周りの閉形式]
$\R^2\setminus\{0\}$ で $\eta=(-y\dd x+x\dd y)/(x^2+y^2)$ とする。係数の微分から $\dd\eta=0$ であるが、単位円に引き戻すと $\dd\theta$ となるので積分は $2\pi$ である。完全な一形式の閉曲線積分は端点での関数値の差で零になる。従って $\eta$ は全領域では完全でない。
:::

### 流れの向きに沿った微分

滑らかな $X$ が作る流れ $\Phi_t$ に対し $\mathcal L_X\omega=(\dd/\dd t)_{t=0}\Phi_t^*\omega$ と定める。引き戻しは外積に分配するので $\mathcal L_X$ は次数零の積の微分則を満たす。さらに $\dd\Phi_t^*=\Phi_t^*\dd$ なので $\mathcal L_X\dd f=\dd\mathcal L_X f$ である。$\dd\iota_X+\iota_X\dd$ も二つの積の微分則を組み合わせると次数零の積の微分則を満たす。関数では両方とも $Xf$、$\dd x^i$ では $\dd X^i$ だから、座標の外積基底に適用して
$$
\mathcal L_X\omega=\dd\iota_X\omega+\iota_X\dd\omega
$$
を得る。例えば $X=x\partial_x+2y\partial_y$ の流れは $(e^tx,e^{2t}y)$ で、$\dd x\wedge\dd y$ は引き戻しで $e^{3t}$ 倍になる。差し込みで $x\dd y-2y\dd x$ を得て外微分すると $3\dd x\wedge\dd y$ となり、倍率の微分と一致する。

### 外微分と相殺を成分で確かめる

二次元の一形式 $\omega=A\dd x+B\dd y$ の外微分は $(B_x-A_y)\dd x\wedge\dd y$ となる。ここで $A,B$ の交差微分が等しければ局所ポテンシャルの候補を作れる。実際、原点を含む星形領域なら $f(x,y)=\int_0^1[x A(tx,ty)+y B(tx,ty)]\dd t$ に微分を入れる。閉じている条件を使うと各偏微分の被積分関数は $tA(tx,ty)$ または $tB(tx,ty)$ の $t$ 微分になる。従って $f_x=A,f_y=B$ である。

:::example[斉次な閉形式の原始関数（1）]
一形式 $\omega=mx^{m-1}y\,\dd x+x^m\dd y$ を星形領域で積分する。基準の値で計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

二つの交差微分は一致する。

$$
\begin{aligned}
\partial_y A=\partial_x B\big|_{m=1} &= 1
\end{aligned}
$$

放射方向に一形式を適用すると、$t$ に関するモーメントになる。

$$
\begin{aligned}
xA(tx,ty)+yB(tx,ty)\big|_{m=1} &= 2 t x y
\end{aligned}
$$

積分した関数をもう一度微分すると原始形式を確認できる。

$$
\begin{aligned}
f(x,y)\big|_{m=1} &= x y
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
xA(tx,ty)+yB(tx,ty)\big|_{x=2} &= 4 t y \\
\frac{\partial^1 (xA(tx,ty)+yB(tx,ty))}{\partial x^1} &= 2 t y \\
\left.\frac{\partial^1 (xA(tx,ty)+yB(tx,ty))}{\partial x^1}\right|_{x=2} &= 2 t y \\
\frac{\partial^2 (xA(tx,ty)+yB(tx,ty))}{\partial x^2} &= 0 \\
\left.\frac{\partial^2 (xA(tx,ty)+yB(tx,ty))}{\partial x^2}\right|_{x=2} &= 0
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
f(x,y)\big|_{x=2} &= 2 y \\
\frac{\partial^1 (f(x,y))}{\partial x^1} &= y \\
\left.\frac{\partial^1 (f(x,y))}{\partial x^1}\right|_{x=2} &= y \\
\frac{\partial^2 (f(x,y))}{\partial x^2} &= 0 \\
\left.\frac{\partial^2 (f(x,y))}{\partial x^2}\right|_{x=2} &= 0
\end{aligned}
$$
:::

:::example[斉次な閉形式の原始関数（2）]
一形式 $\omega=mx^{m-1}y\,\dd x+x^m\dd y$ を星形領域で積分する。値を変えて依存性を調べる。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

二つの交差微分は一致する。

$$
\begin{aligned}
\partial_y A=\partial_x B\big|_{m=2} &= 2 x
\end{aligned}
$$

放射方向に一形式を適用すると、$t$ に関するモーメントになる。

$$
\begin{aligned}
xA(tx,ty)+yB(tx,ty)\big|_{m=2} &= 3 t^{2} x^{2} y
\end{aligned}
$$

積分した関数をもう一度微分すると原始形式を確認できる。

$$
\begin{aligned}
f(x,y)\big|_{m=2} &= x^{2} y
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\partial_y A=\partial_x B\big|_{x=2} &= 4 \\
\frac{\partial^1 (\partial_y A=\partial_x B)}{\partial x^1} &= 2 \\
\left.\frac{\partial^1 (\partial_y A=\partial_x B)}{\partial x^1}\right|_{x=2} &= 2 \\
\frac{\partial^2 (\partial_y A=\partial_x B)}{\partial x^2} &= 0 \\
\left.\frac{\partial^2 (\partial_y A=\partial_x B)}{\partial x^2}\right|_{x=2} &= 0
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
xA(tx,ty)+yB(tx,ty)\big|_{x=2} &= 12 t^{2} y \\
\frac{\partial^1 (xA(tx,ty)+yB(tx,ty))}{\partial x^1} &= 6 t^{2} x y \\
\left.\frac{\partial^1 (xA(tx,ty)+yB(tx,ty))}{\partial x^1}\right|_{x=2} &= 12 t^{2} y \\
\frac{\partial^2 (xA(tx,ty)+yB(tx,ty))}{\partial x^2} &= 6 t^{2} y \\
\left.\frac{\partial^2 (xA(tx,ty)+yB(tx,ty))}{\partial x^2}\right|_{x=2} &= 6 t^{2} y
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
f(x,y)\big|_{x=2} &= 4 y \\
\frac{\partial^1 (f(x,y))}{\partial x^1} &= 2 x y \\
\left.\frac{\partial^1 (f(x,y))}{\partial x^1}\right|_{x=2} &= 4 y \\
\frac{\partial^2 (f(x,y))}{\partial x^2} &= 2 y \\
\left.\frac{\partial^2 (f(x,y))}{\partial x^2}\right|_{x=2} &= 2 y
\end{aligned}
$$
:::

:::example[斉次な閉形式の原始関数（3）]
一形式 $\omega=mx^{m-1}y\,\dd x+x^m\dd y$ を星形領域で積分する。さらに値を変えて比較する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

二つの交差微分は一致する。

$$
\begin{aligned}
\partial_y A=\partial_x B\big|_{m=3} &= 3 x^{2}
\end{aligned}
$$

放射方向に一形式を適用すると、$t$ に関するモーメントになる。

$$
\begin{aligned}
xA(tx,ty)+yB(tx,ty)\big|_{m=3} &= 4 t^{3} x^{3} y
\end{aligned}
$$

積分した関数をもう一度微分すると原始形式を確認できる。

$$
\begin{aligned}
f(x,y)\big|_{m=3} &= x^{3} y
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\partial_y A=\partial_x B\big|_{x=2} &= 12 \\
\frac{\partial^1 (\partial_y A=\partial_x B)}{\partial x^1} &= 6 x \\
\left.\frac{\partial^1 (\partial_y A=\partial_x B)}{\partial x^1}\right|_{x=2} &= 12 \\
\frac{\partial^2 (\partial_y A=\partial_x B)}{\partial x^2} &= 6 \\
\left.\frac{\partial^2 (\partial_y A=\partial_x B)}{\partial x^2}\right|_{x=2} &= 6
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
xA(tx,ty)+yB(tx,ty)\big|_{x=2} &= 32 t^{3} y \\
\frac{\partial^1 (xA(tx,ty)+yB(tx,ty))}{\partial x^1} &= 12 t^{3} x^{2} y \\
\left.\frac{\partial^1 (xA(tx,ty)+yB(tx,ty))}{\partial x^1}\right|_{x=2} &= 48 t^{3} y \\
\frac{\partial^2 (xA(tx,ty)+yB(tx,ty))}{\partial x^2} &= 24 t^{3} x y \\
\left.\frac{\partial^2 (xA(tx,ty)+yB(tx,ty))}{\partial x^2}\right|_{x=2} &= 48 t^{3} y
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
f(x,y)\big|_{x=2} &= 8 y \\
\frac{\partial^1 (f(x,y))}{\partial x^1} &= 3 x^{2} y \\
\left.\frac{\partial^1 (f(x,y))}{\partial x^1}\right|_{x=2} &= 12 y \\
\frac{\partial^2 (f(x,y))}{\partial x^2} &= 6 x y \\
\left.\frac{\partial^2 (f(x,y))}{\partial x^2}\right|_{x=2} &= 12 y
\end{aligned}
$$
:::

:::example[斉次な閉形式の原始関数（4）]
一形式 $\omega=mx^{m-1}y\,\dd x+x^m\dd y$ を星形領域で積分する。最後にこの値でも計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=4
$$

二つの交差微分は一致する。

$$
\begin{aligned}
\partial_y A=\partial_x B\big|_{m=4} &= 4 x^{3}
\end{aligned}
$$

放射方向に一形式を適用すると、$t$ に関するモーメントになる。

$$
\begin{aligned}
xA(tx,ty)+yB(tx,ty)\big|_{m=4} &= 5 t^{4} x^{4} y
\end{aligned}
$$

積分した関数をもう一度微分すると原始形式を確認できる。

$$
\begin{aligned}
f(x,y)\big|_{m=4} &= x^{4} y
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\partial_y A=\partial_x B\big|_{x=2} &= 32 \\
\frac{\partial^1 (\partial_y A=\partial_x B)}{\partial x^1} &= 12 x^{2} \\
\left.\frac{\partial^1 (\partial_y A=\partial_x B)}{\partial x^1}\right|_{x=2} &= 48 \\
\frac{\partial^2 (\partial_y A=\partial_x B)}{\partial x^2} &= 24 x \\
\left.\frac{\partial^2 (\partial_y A=\partial_x B)}{\partial x^2}\right|_{x=2} &= 48
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
xA(tx,ty)+yB(tx,ty)\big|_{x=2} &= 80 t^{4} y \\
\frac{\partial^1 (xA(tx,ty)+yB(tx,ty))}{\partial x^1} &= 20 t^{4} x^{3} y \\
\left.\frac{\partial^1 (xA(tx,ty)+yB(tx,ty))}{\partial x^1}\right|_{x=2} &= 160 t^{4} y \\
\frac{\partial^2 (xA(tx,ty)+yB(tx,ty))}{\partial x^2} &= 60 t^{4} x^{2} y \\
\left.\frac{\partial^2 (xA(tx,ty)+yB(tx,ty))}{\partial x^2}\right|_{x=2} &= 240 t^{4} y
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
f(x,y)\big|_{x=2} &= 16 y \\
\frac{\partial^1 (f(x,y))}{\partial x^1} &= 4 x^{3} y \\
\left.\frac{\partial^1 (f(x,y))}{\partial x^1}\right|_{x=2} &= 32 y \\
\frac{\partial^2 (f(x,y))}{\partial x^2} &= 12 x^{2} y \\
\left.\frac{\partial^2 (f(x,y))}{\partial x^2}\right|_{x=2} &= 48 y
\end{aligned}
$$
:::

## Hodge 双対と余微分・Laplacian

### 計量が与える双対

外微分は滑らかな構造から定まる。一方、向き付き線分とそれに垂直な面を対応させるには長さと角度が要る。そこで計量と向きを使って $p$ 形式と $(n-p)$ 形式を対応させる。

:::definition[Hodge 双対と余微分]{#def-mp-hodge}
向き付き Riemann 多様体の体積形式を $\mathrm{vol}$ とする。正向き正規直交余標構の増加添字の外積を正規直交とする内積を $\inner{\alpha}{\beta}$ と書く。:term[ホッジ双対]{en="Hodge star" yomi="ほっじそうつい"} $\star$ を $\alpha\wedge\star\beta=\inner{\alpha}{\beta}\mathrm{vol}$ で定める。$p$ 形式上の :term[余微分]{en="codifferential" yomi="よびぶん"} を $\delta=(-1)^{n(p+1)+1}\star\dd\star$ とする。$p=0$ では $\delta=0$ とする。
:::

基底の形式と重ならない添字を埋めれば $\star$ を構成できる。たとえば $\star(\dd x\wedge\dd z)=-\dd y$ であり、二回作用すると $p$ 個と $n-p$ 個の移動の符号から $\star^2=(-1)^{p(n-p)}$ となる。負符号方向が $s$ 本ある擬 Riemann 計量では、直交基底の各形式の内積にその方向の符号を掛けるため $\star^2=(-1)^{p(n-p)+s}$ となる。

:::proposition[形式的随伴]{#prop-mp-codifferential-adjoint}
滑らかな実形式 $\alpha\in\Omega^{p-1}(M)$ と $\beta\in\Omega^p(M)$ のうち少なくとも一方が内部に含まれるコンパクトな台を持つなら、$\int\inner{\dd\alpha}{\beta}\mathrm{vol}=\int\inner{\alpha}{\delta\beta}\mathrm{vol}$ である。
:::

:::proof
積の外微分を $\dd(\alpha\wedge\star\beta)$ に適用し、次節の積分公式を使うと $\int\dd\alpha\wedge\star\beta=(-1)^p\int\alpha\wedge\dd\star\beta$ を得る。$\star^2$ の符号から $\star\delta\beta=(-1)^p\dd\star\beta$ である。これを戻すと示す式になる。
:::

:term[ホッジ Laplacian]{en="Hodge Laplacian" yomi="ほっじらぷらしあん"} を $\Delta_H=\dd\delta+\delta\dd$ とする。上の関係から 内部に含まれるコンパクトな台を持つ $\omega$ に対して $\int\inner{\omega}{\Delta_H\omega}\mathrm{vol}=\|\dd\omega\|^2+\|\delta\omega\|^2$ である。Euclid 空間で $\delta=-\sum_i\iota_{\partial_i}\partial_i$ となる。基底の外積を一つずつ調べると $\iota_{\partial_i}(\dd x^j\wedge\omega)+\dd x^j\wedge\iota_{\partial_i}\omega=\delta_i^j\omega$ となるため、$\Delta_H=-\sum_i\partial_i^2$ が各成分に作用する。スカラーに関する符号は前章の $\Delta$ と逆である。

:::example[二次関数と一形式]
$\R^3$ で $f=x^2+y^2+z^2$ なら $\Delta_H f=-6$ である。また $\alpha=x^2\dd y$ に対して $\delta\alpha=0$、$\Delta_H\alpha=-2\dd y$ となる。従って $\delta\dd\alpha=-2\dd y$ である。
:::

:::example[極座標の Hodge 双対]
平面を $\dd r\wedge r\dd\theta$ で向き付けると $\star\dd r=r\dd\theta$、$\star\dd\theta=-\dd r/r$ となる。例えば $f=r^3$ は $\star\dd f=3r^3\dd\theta$ を与える。微分は $9r^2\dd r\wedge\dd\theta$、さらに双対を取ると $9r$ である。$\Delta_Hf=-9r$ となり、$\Delta=r^{-1}\partial_r(r\partial_r)$ と符号まで照合できる。
:::

### 座標成分と体積因子

向き付き Riemann 多様体の正向き座標で、計量の行列式を $g=\det(g_{ij})$ とすると $\mathrm{vol}=\sqrt g\,\dd x^1\wedge\cdots\wedge\dd x^n$ である。これは計量で測る単位座標箱の体積が Gram 行列式の平方根になるためである。座標一形式 $\alpha=\alpha_i\dd x^i$ に対し $A^i=g^{ij}\alpha_j$ と置けば $\star\alpha=\iota_A\mathrm{vol}$ となる。これは任意の一形式 $\beta$ を左から外積すると、両辺が $\beta(A)\mathrm{vol}$ となることから分かる。

差し込みの展開を微分すると $\dd\star\alpha=\partial_i(\sqrt g\,g^{ij}\alpha_j)\,\dd x^1\wedge\cdots\wedge\dd x^n$ である。$\star\mathrm{vol}=1$ だから
$$
\delta\alpha=-\frac1{\sqrt g}\partial_i(\sqrt g\,g^{ij}\alpha_j),\qquad
\Delta_H f=-\frac1{\sqrt g}\partial_i(\sqrt g\,g^{ij}\partial_j f)
$$
となる。例えば単位球面の局所座標 $(\theta,\phi)$ では計量が $\dd\theta^2+\sin^2\theta\,\dd\phi^2$ であり、$\cos\theta$ は $\Delta_H(\cos\theta)=2\cos\theta$ を満たす。実際、$\partial_\theta(\sin\theta\,\partial_\theta\cos\theta)=-2\sin\theta\cos\theta$ である。極ではこの座標の体積係数が零になるが、$\cos\theta$ は球面上の滑らかな高さの関数なので、等式は連続性により極にも延びる。

### 余微分の符号を試す

向き付き Euclid 平面で $*\dd x=\dd y$、$*\dd y=-\dd x$ である。従って $**=-1$ が一形式に作用する。さらに $\delta=-*\dd*$ とすると $\delta(A\dd x+B\dd y)=-(A_x+B_y)$ であり、$d\delta+\delta d$ は直交成分ごとの負の Laplacian になる。最後の式では混合微分が相殺する。部分積分による随伴と比較すると符号を確かめやすい。

:::example[勾配一形式の余微分（1）]
関数 $f=x^m y^2$ の微分と Laplacian を比較する。基準の値で計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

一形式の第一成分を求める。

$$
\begin{aligned}
(\dd f)_x\big|_{m=1} &= y^{2}
\end{aligned}
$$

第二成分を求める。

$$
\begin{aligned}
(\dd f)_y\big|_{m=1} &= 2 x y
\end{aligned}
$$

二成分をそれぞれの座標で再度微分して負号を付ける。

$$
\begin{aligned}
\delta\dd f\big|_{m=1} &= - 2 x
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
(\dd f)_y\big|_{x=2} &= 4 y \\
\frac{\partial^1 ((\dd f)_y)}{\partial x^1} &= 2 y \\
\left.\frac{\partial^1 ((\dd f)_y)}{\partial x^1}\right|_{x=2} &= 2 y
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\delta\dd f\big|_{x=2} &= -4 \\
\frac{\partial^1 (\delta\dd f)}{\partial x^1} &= -2 \\
\left.\frac{\partial^1 (\delta\dd f)}{\partial x^1}\right|_{x=2} &= -2
\end{aligned}
$$
:::

:::example[勾配一形式の余微分（2）]
関数 $f=x^m y^2$ の微分と Laplacian を比較する。値を変えて依存性を調べる。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

一形式の第一成分を求める。

$$
\begin{aligned}
(\dd f)_x\big|_{m=2} &= 2 x y^{2}
\end{aligned}
$$

第二成分を求める。

$$
\begin{aligned}
(\dd f)_y\big|_{m=2} &= 2 x^{2} y
\end{aligned}
$$

二成分をそれぞれの座標で再度微分して負号を付ける。

$$
\begin{aligned}
\delta\dd f\big|_{m=2} &= - 2 x^{2} - 2 y^{2}
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
(\dd f)_x\big|_{x=2} &= 4 y^{2} \\
\frac{\partial^1 ((\dd f)_x)}{\partial x^1} &= 2 y^{2} \\
\left.\frac{\partial^1 ((\dd f)_x)}{\partial x^1}\right|_{x=2} &= 2 y^{2}
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
(\dd f)_y\big|_{x=2} &= 8 y \\
\frac{\partial^1 ((\dd f)_y)}{\partial x^1} &= 4 x y \\
\left.\frac{\partial^1 ((\dd f)_y)}{\partial x^1}\right|_{x=2} &= 8 y
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\delta\dd f\big|_{x=2} &= - 2 y^{2} - 8 \\
\frac{\partial^1 (\delta\dd f)}{\partial x^1} &= - 4 x \\
\left.\frac{\partial^1 (\delta\dd f)}{\partial x^1}\right|_{x=2} &= -8
\end{aligned}
$$
:::

:::example[勾配一形式の余微分（3）]
関数 $f=x^m y^2$ の微分と Laplacian を比較する。さらに値を変えて比較する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

一形式の第一成分を求める。

$$
\begin{aligned}
(\dd f)_x\big|_{m=3} &= 3 x^{2} y^{2}
\end{aligned}
$$

第二成分を求める。

$$
\begin{aligned}
(\dd f)_y\big|_{m=3} &= 2 x^{3} y
\end{aligned}
$$

二成分をそれぞれの座標で再度微分して負号を付ける。

$$
\begin{aligned}
\delta\dd f\big|_{m=3} &= 2 x \left(- x^{2} - 3 y^{2}\right)
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
(\dd f)_x\big|_{x=2} &= 12 y^{2} \\
\frac{\partial^1 ((\dd f)_x)}{\partial x^1} &= 6 x y^{2} \\
\left.\frac{\partial^1 ((\dd f)_x)}{\partial x^1}\right|_{x=2} &= 12 y^{2}
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
(\dd f)_y\big|_{x=2} &= 16 y \\
\frac{\partial^1 ((\dd f)_y)}{\partial x^1} &= 6 x^{2} y \\
\left.\frac{\partial^1 ((\dd f)_y)}{\partial x^1}\right|_{x=2} &= 24 y
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\delta\dd f\big|_{x=2} &= - 12 y^{2} - 16 \\
\frac{\partial^1 (\delta\dd f)}{\partial x^1} &= - 6 x^{2} - 6 y^{2} \\
\left.\frac{\partial^1 (\delta\dd f)}{\partial x^1}\right|_{x=2} &= - 6 y^{2} - 24
\end{aligned}
$$
:::

:::example[勾配一形式の余微分（4）]
関数 $f=x^m y^2$ の微分と Laplacian を比較する。最後にこの値でも計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=4
$$

一形式の第一成分を求める。

$$
\begin{aligned}
(\dd f)_x\big|_{m=4} &= 4 x^{3} y^{2}
\end{aligned}
$$

第二成分を求める。

$$
\begin{aligned}
(\dd f)_y\big|_{m=4} &= 2 x^{4} y
\end{aligned}
$$

二成分をそれぞれの座標で再度微分して負号を付ける。

$$
\begin{aligned}
\delta\dd f\big|_{m=4} &= 2 x^{2} \left(- x^{2} - 6 y^{2}\right)
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
(\dd f)_x\big|_{x=2} &= 32 y^{2} \\
\frac{\partial^1 ((\dd f)_x)}{\partial x^1} &= 12 x^{2} y^{2} \\
\left.\frac{\partial^1 ((\dd f)_x)}{\partial x^1}\right|_{x=2} &= 48 y^{2}
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
(\dd f)_y\big|_{x=2} &= 32 y \\
\frac{\partial^1 ((\dd f)_y)}{\partial x^1} &= 8 x^{3} y \\
\left.\frac{\partial^1 ((\dd f)_y)}{\partial x^1}\right|_{x=2} &= 64 y
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\delta\dd f\big|_{x=2} &= - 48 y^{2} - 32 \\
\frac{\partial^1 (\delta\dd f)}{\partial x^1} &= 8 x \left(- x^{2} - 3 y^{2}\right) \\
\left.\frac{\partial^1 (\delta\dd f)}{\partial x^1}\right|_{x=2} &= - 48 y^{2} - 64
\end{aligned}
$$
:::

## 一般化された Stokes の定理

### 境界を向き付ける

区間の端点の符号を $b$ で $+$、$a$ で $-$ とする約束は、外向きの方向を先に置く規則の一次元版である。高次元では外向き法線に続けた境界の基底が内部と同じ向きになるようにする。

:::theorem[一般化された Stokes の定理]{#thm-mp-general-stokes}
向き付き滑らかな $n$ 次元多様体 $M$ とその境界への包含 $\jmath:\partial M\hookrightarrow M$ に対し、台がコンパクトな滑らかな $(n-1)$ 形式 $\omega$ は $\int_M\dd\omega=\int_{\partial M}\jmath^*\omega$ を満たす。境界は外向き先頭の規約で向き付ける。
:::

:::proof
単位の分割によって $\omega$ を各座標の中に台が収まる有限個の形式の和に分ける。分割の和は台の近くで $1$ なので、和を外微分しても元と一致する。内点の正向き座標で $\omega=\sum_i(-1)^{i-1}a_i\dd x^1\wedge\cdots\widehat{\dd x^i}\cdots\wedge\dd x^n$ と書くと、$\dd\omega=(\sum_i\partial_i a_i)\dd x^1\wedge\cdots\wedge\dd x^n$ である。実軸全体の積分は微積分の基本定理で零となる。

境界点では $x^n\ge0$ の半空間をとる。同じ計算で $i<n$ の寄与は零で、$i=n$ は $-\int a_n(x',0)\dd x'$ となる。外向き方向は $-\partial_n$ であり、座標順の境界基底には $(-1)^n$ の向き補正が必要になる。引き戻しに残る係数 $(-1)^{n-1}a_n$ と掛け合わせると $-a_n$ を得る。分割した各項を足せばよい。
:::

:::example[円周と円板]
$\omega=(-y\dd x+x\dd y)/2$ では $\dd\omega=\dd x\wedge\dd y$ である。半径 $3$ の円周上で $\omega=\tfrac92\dd\theta$ なので積分は $9\pi$ である。これは円板の面積に等しい。
:::

:::example[境界の相対的な向き]
$1\le r\le2$ の環で外微分の節の $\eta=(-y\dd x+x\dd y)/(x^2+y^2)$ を使うと、外側の円では $2\pi$、内側の円は時計回りに向き付けるので積分は $-2\pi$ である。従って境界全体では零となり、$\dd\eta=0$ と一致する。
:::

### 境界の向きと積分の相殺

平面で正の面積形式を $\dd x\wedge\dd y$ としたとき、長方形の境界は反時計回りになる。隣接する面素をつないだ場合、内部の境界は二つの向きで一度ずつ出るため相殺する。この相殺が局所的な微分式を領域全体に足せる理由である。曲線を逆向きに走ると一形式の積分だけで符号を追える。面積の絶対値を使う計算では、境界の向きは別に指定する。

:::example[長方形の周りの循環（1）]
長方形 $0\le x\le a,0\le y\le2$ に $\omega=x^m y\,\dd y$ を取る。基準の値で計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

水平辺では $\dd y=0$ であり、左辺では係数が零である。右辺の積分だけが残る。

$$
\begin{aligned}
\int_{\partial R}\omega\big|_{m=1} &= 2 a
\end{aligned}
$$

外微分では $x$ の微分だけが残る。

$$
\begin{aligned}
[\dd\omega]_{xy}\big|_{m=1} &= y
\end{aligned}
$$

面積積分を二つの独立な積分に分ける。

$$
\begin{aligned}
\int_R\dd\omega\big|_{m=1} &= 2 a
\end{aligned}
$$

観測点 $a$=$3$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\int_{\partial R}\omega\big|_{a=3} &= 6 \\
\frac{\partial^1 (\int_{\partial R}\omega)}{\partial a^1} &= 2 \\
\left.\frac{\partial^1 (\int_{\partial R}\omega)}{\partial a^1}\right|_{a=3} &= 2
\end{aligned}
$$

観測点 $a$=$3$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\int_R\dd\omega\big|_{a=3} &= 6 \\
\frac{\partial^1 (\int_R\dd\omega)}{\partial a^1} &= 2 \\
\left.\frac{\partial^1 (\int_R\dd\omega)}{\partial a^1}\right|_{a=3} &= 2
\end{aligned}
$$
:::

:::example[長方形の周りの循環（2）]
長方形 $0\le x\le a,0\le y\le2$ に $\omega=x^m y\,\dd y$ を取る。値を変えて依存性を調べる。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

水平辺では $\dd y=0$ であり、左辺では係数が零である。右辺の積分だけが残る。

$$
\begin{aligned}
\int_{\partial R}\omega\big|_{m=2} &= 2 a^{2}
\end{aligned}
$$

外微分では $x$ の微分だけが残る。

$$
\begin{aligned}
[\dd\omega]_{xy}\big|_{m=2} &= 2 x y
\end{aligned}
$$

面積積分を二つの独立な積分に分ける。

$$
\begin{aligned}
\int_R\dd\omega\big|_{m=2} &= 2 a^{2}
\end{aligned}
$$

観測点 $a$=$3$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\int_{\partial R}\omega\big|_{a=3} &= 18 \\
\frac{\partial^1 (\int_{\partial R}\omega)}{\partial a^1} &= 4 a \\
\left.\frac{\partial^1 (\int_{\partial R}\omega)}{\partial a^1}\right|_{a=3} &= 12
\end{aligned}
$$

観測点 $a$=$3$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\int_R\dd\omega\big|_{a=3} &= 18 \\
\frac{\partial^1 (\int_R\dd\omega)}{\partial a^1} &= 4 a \\
\left.\frac{\partial^1 (\int_R\dd\omega)}{\partial a^1}\right|_{a=3} &= 12
\end{aligned}
$$
:::

:::example[長方形の周りの循環（3）]
長方形 $0\le x\le a,0\le y\le2$ に $\omega=x^m y\,\dd y$ を取る。さらに値を変えて比較する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

水平辺では $\dd y=0$ であり、左辺では係数が零である。右辺の積分だけが残る。

$$
\begin{aligned}
\int_{\partial R}\omega\big|_{m=3} &= 2 a^{3}
\end{aligned}
$$

外微分では $x$ の微分だけが残る。

$$
\begin{aligned}
[\dd\omega]_{xy}\big|_{m=3} &= 3 x^{2} y
\end{aligned}
$$

面積積分を二つの独立な積分に分ける。

$$
\begin{aligned}
\int_R\dd\omega\big|_{m=3} &= 2 a^{3}
\end{aligned}
$$

観測点 $a$=$3$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\int_{\partial R}\omega\big|_{a=3} &= 54 \\
\frac{\partial^1 (\int_{\partial R}\omega)}{\partial a^1} &= 6 a^{2} \\
\left.\frac{\partial^1 (\int_{\partial R}\omega)}{\partial a^1}\right|_{a=3} &= 54
\end{aligned}
$$

観測点 $a$=$3$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\int_R\dd\omega\big|_{a=3} &= 54 \\
\frac{\partial^1 (\int_R\dd\omega)}{\partial a^1} &= 6 a^{2} \\
\left.\frac{\partial^1 (\int_R\dd\omega)}{\partial a^1}\right|_{a=3} &= 54
\end{aligned}
$$
:::

:::example[長方形の周りの循環（4）]
長方形 $0\le x\le a,0\le y\le2$ に $\omega=x^m y\,\dd y$ を取る。最後にこの値でも計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=4
$$

水平辺では $\dd y=0$ であり、左辺では係数が零である。右辺の積分だけが残る。

$$
\begin{aligned}
\int_{\partial R}\omega\big|_{m=4} &= 2 a^{4}
\end{aligned}
$$

外微分では $x$ の微分だけが残る。

$$
\begin{aligned}
[\dd\omega]_{xy}\big|_{m=4} &= 4 x^{3} y
\end{aligned}
$$

面積積分を二つの独立な積分に分ける。

$$
\begin{aligned}
\int_R\dd\omega\big|_{m=4} &= 2 a^{4}
\end{aligned}
$$

観測点 $a$=$3$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\int_{\partial R}\omega\big|_{a=3} &= 162 \\
\frac{\partial^1 (\int_{\partial R}\omega)}{\partial a^1} &= 8 a^{3} \\
\left.\frac{\partial^1 (\int_{\partial R}\omega)}{\partial a^1}\right|_{a=3} &= 216
\end{aligned}
$$

観測点 $a$=$3$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\int_R\dd\omega\big|_{a=3} &= 162 \\
\frac{\partial^1 (\int_R\dd\omega)}{\partial a^1} &= 8 a^{3} \\
\left.\frac{\partial^1 (\int_R\dd\omega)}{\partial a^1}\right|_{a=3} &= 216
\end{aligned}
$$
:::

## ベクトル解析の公式の微分形式表示

### 循環形式と流束形式

$\vb A^\flat$ をベクトル $\vb A$ と計量で対応する一形式 $\vb A^\flat(\vb v)=g(\vb A,\vb v)$とする。$\iota_{\vb A}\mathrm{vol}=\star\vb A^\flat$ は流束の二形式である。等式は正規直交基底の各ベクトルに差し込みを行えば分かる。

:::proposition[ベクトル解析との対応]{#prop-mp-vector-forms}
右手系の Euclid 3 次元空間で $(\grad f)^\flat=\dd f$、$(\rot\vb A)^\flat=\star\dd\vb A^\flat$、$\divg\vb A=\star\dd\star\vb A^\flat$ である。
:::

:::proof
第一式は $\dd f(\vb v)=\grad f\cdot\vb v$ による。第二式の $\dd x^1$ 成分は $\dd x^2\wedge\dd x^3$ の係数 $\partial_2 A_3-\partial_3 A_2$ であり、巡回順の残りも同様である。第三式は $\star\vb A^\flat=A_1\dd x^2\wedge\dd x^3+A_2\dd x^3\wedge\dd x^1+A_3\dd x^1\wedge\dd x^2$ を微分すれば従う。
:::

従って $\dd^2=0$ の二つの現れが $\rot\grad f=0$ と $\divg\rot\vb A=0$ である。また :ref[thm-mp-general-stokes] を一形式に用いると循環の Stokes の定理、二形式に用いると発散定理が出る。

:::example[角度の成分と接線成分]
円筒座標で $\vb A=r^2\hat{\vb e}_\phi$ とする。$\vb A^\flat=r^3\dd\phi$ なので $\dd\vb A^\flat=3r^2\dd r\wedge\dd\phi$ である。正規直交余標構は $(\dd r,r\dd\phi,\dd z)$ だから $\star\dd\vb A^\flat=3r\dd z$ となり、式 :eqref[eq-mp-orthogonal-operators] と一致する。
:::

:::example[面を流れる膨張場]
$X=x\partial_x+y\partial_y+z\partial_z$ では $\iota_X\mathrm{vol}=x\dd y\wedge\dd z+y\dd z\wedge\dd x+z\dd x\wedge\dd y$ で、その微分は $3\mathrm{vol}$ である。辺の長さが $2\,\mathrm m$ の立方体で、各成分の比例定数を $1\,\mathrm{s^{-1}}$ とすれば正味流量は $24\,\mathrm{m^3\,s^{-1}}$ になる。立方体 $[0,2]^3$ では座標が零の三面の寄与が零、反対の三面が各 $8\,\mathrm{m^3\,s^{-1}}$ であるため、面ごとの計算でも確認できる。
:::

### 流束形式と回転を区別する

平坦な三次元空間で、ベクトル $V$ の循環形式 $\alpha=V^\flat$ は曲線に沿って積分する。流束形式 $*\alpha$ は面上で積分する。前者に $\dd$ を作用させると回転の流束になり、後者に作用させると発散の体積形式になる。$\dd^2=0$ を二回使うと、勾配の回転と回転の発散が零になる。これは座標ごとに混合微分を比べる方法と同じ内容である。

:::example[軸の周囲の流れ（1）]
速度 $V=(-y x^m,x^{m+1},0)$ の局所回転と発散を求める。基準の値で計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

循環形式の $y$ 成分は $x^{m+1}$ である。$x$ 成分の $y$ 微分を引くと回転の軸成分になる。

$$
\begin{aligned}
(\operatorname{curl}V)_z\big|_{m=1} &= 3 x
\end{aligned}
$$

発散に入るのは各成分を同じ変数で微分した項である。

$$
\begin{aligned}
\operatorname{div}V\big|_{m=1} &= - y
\end{aligned}
$$

回転は $x$ に依存するが、唯一の成分は $z$ 向きなのでその発散は消える。

$$
\begin{aligned}
\operatorname{div}\operatorname{curl}V &= 0
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
(\operatorname{curl}V)_z\big|_{x=2} &= 6 \\
\frac{\partial^1 ((\operatorname{curl}V)_z)}{\partial x^1} &= 3 \\
\left.\frac{\partial^1 ((\operatorname{curl}V)_z)}{\partial x^1}\right|_{x=2} &= 3
\end{aligned}
$$
:::

:::example[軸の周囲の流れ（2）]
速度 $V=(-y x^m,x^{m+1},0)$ の局所回転と発散を求める。値を変えて依存性を調べる。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

循環形式の $y$ 成分は $x^{m+1}$ である。$x$ 成分の $y$ 微分を引くと回転の軸成分になる。

$$
\begin{aligned}
(\operatorname{curl}V)_z\big|_{m=2} &= 4 x^{2}
\end{aligned}
$$

発散に入るのは各成分を同じ変数で微分した項である。

$$
\begin{aligned}
\operatorname{div}V\big|_{m=2} &= - 2 x y
\end{aligned}
$$

回転は $x$ に依存するが、唯一の成分は $z$ 向きなのでその発散は消える。

$$
\begin{aligned}
\operatorname{div}\operatorname{curl}V &= 0
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
(\operatorname{curl}V)_z\big|_{x=2} &= 16 \\
\frac{\partial^1 ((\operatorname{curl}V)_z)}{\partial x^1} &= 8 x \\
\left.\frac{\partial^1 ((\operatorname{curl}V)_z)}{\partial x^1}\right|_{x=2} &= 16
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\operatorname{div}V\big|_{x=2} &= - 4 y \\
\frac{\partial^1 (\operatorname{div}V)}{\partial x^1} &= - 2 y \\
\left.\frac{\partial^1 (\operatorname{div}V)}{\partial x^1}\right|_{x=2} &= - 2 y
\end{aligned}
$$
:::

:::example[軸の周囲の流れ（3）]
速度 $V=(-y x^m,x^{m+1},0)$ の局所回転と発散を求める。さらに値を変えて比較する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

循環形式の $y$ 成分は $x^{m+1}$ である。$x$ 成分の $y$ 微分を引くと回転の軸成分になる。

$$
\begin{aligned}
(\operatorname{curl}V)_z\big|_{m=3} &= 5 x^{3}
\end{aligned}
$$

発散に入るのは各成分を同じ変数で微分した項である。

$$
\begin{aligned}
\operatorname{div}V\big|_{m=3} &= - 3 x^{2} y
\end{aligned}
$$

回転は $x$ に依存するが、唯一の成分は $z$ 向きなのでその発散は消える。

$$
\begin{aligned}
\operatorname{div}\operatorname{curl}V &= 0
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
(\operatorname{curl}V)_z\big|_{x=2} &= 40 \\
\frac{\partial^1 ((\operatorname{curl}V)_z)}{\partial x^1} &= 15 x^{2} \\
\left.\frac{\partial^1 ((\operatorname{curl}V)_z)}{\partial x^1}\right|_{x=2} &= 60
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\operatorname{div}V\big|_{x=2} &= - 12 y \\
\frac{\partial^1 (\operatorname{div}V)}{\partial x^1} &= - 6 x y \\
\left.\frac{\partial^1 (\operatorname{div}V)}{\partial x^1}\right|_{x=2} &= - 12 y
\end{aligned}
$$
:::

:::example[軸の周囲の流れ（4）]
速度 $V=(-y x^m,x^{m+1},0)$ の局所回転と発散を求める。最後にこの値でも計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=4
$$

循環形式の $y$ 成分は $x^{m+1}$ である。$x$ 成分の $y$ 微分を引くと回転の軸成分になる。

$$
\begin{aligned}
(\operatorname{curl}V)_z\big|_{m=4} &= 6 x^{4}
\end{aligned}
$$

発散に入るのは各成分を同じ変数で微分した項である。

$$
\begin{aligned}
\operatorname{div}V\big|_{m=4} &= - 4 x^{3} y
\end{aligned}
$$

回転は $x$ に依存するが、唯一の成分は $z$ 向きなのでその発散は消える。

$$
\begin{aligned}
\operatorname{div}\operatorname{curl}V &= 0
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
(\operatorname{curl}V)_z\big|_{x=2} &= 96 \\
\frac{\partial^1 ((\operatorname{curl}V)_z)}{\partial x^1} &= 24 x^{3} \\
\left.\frac{\partial^1 ((\operatorname{curl}V)_z)}{\partial x^1}\right|_{x=2} &= 192
\end{aligned}
$$

観測点 $x$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\operatorname{div}V\big|_{x=2} &= - 32 y \\
\frac{\partial^1 (\operatorname{div}V)}{\partial x^1} &= - 12 x^{2} y \\
\left.\frac{\partial^1 (\operatorname{div}V)}{\partial x^1}\right|_{x=2} &= - 48 y
\end{aligned}
$$
:::

## Maxwell 方程式の微分形式表示

### 電磁場と電荷保存

SI 単位系で真空の構成関係を用いる。光速を $c$、誘電率を $\epsilon_0$、透磁率を $\mu_0$ とし、空間の Hodge 双対を $\star_3$、空間の外微分を $\dd_3$ とする。電場と磁場を一形式 $E=\vb E^\flat$ と二形式 $B=\star_3\vb B^\flat$ で表す。$D=\epsilon_0\star_3 E$、$H=\mu_0^{-1}\star_3 B$、$j=\star_3\vb J^\flat$ と置く。$E$ の電場成分の単位は $\mathrm{V\,m^{-1}}$、$B$ の磁束密度成分の単位は $\mathrm T$ である。$E$ の線積分は $\mathrm V$ の単位を持ち、静電場では始点と終点の電位の差を与える。一方、$B$ の面積分は磁束を与える。

:::proposition[時空の電磁方程式]{#prop-mp-maxwell-forms}
$F=B+E\wedge\dd t$、$G=D-H\wedge\dd t$、$\mathcal J=\rho\,\mathrm{vol}_3-\dd t\wedge j$ と定めると、Maxwell 方程式は次の二式にまとまる。

$$
\dd F=0,\qquad\dd G=\mathcal J .
\label{eq-mp-maxwell-forms}
$$
:::

:::proof
$\dd F=\dd_3 B+\dd t\wedge(\partial_t B+\dd_3 E)$ である。空間成分は磁場の発散が零という式、時間成分は Faraday の法則である。また $\dd G=\dd_3D+\dd t\wedge(\partial_tD-\dd_3H)$ である。これを $\mathcal J$ と比較すると、電気の Gauss の法則と変位電流を含む Ampère の法則を得る。各導出は逆向きにも成立する。
:::

時空の向きを $c\dd t\wedge\mathrm{vol}_3$、計量を $-c^2\dd t^2+\dd x^2+\dd y^2+\dd z^2$ とする。このとき $\star B=c\dd t\wedge\star_3B$、$\star(c\dd t\wedge E)=-\star_3E$ である。後者の負符号は時間方向の内積による。$\epsilon_0\mu_0c^2=1$ より $G=\epsilon_0c\star F$ となる。この構成関係が真空の計量と単位を式に入れている。

$\dd^2G=0$ により $\dd\mathcal J=\dd t\wedge(\partial_t\rho\,\mathrm{vol}_3+\dd_3j)=0$ であり、電荷保存 $\partial_t\rho+\divg\vb J=0$ を得る。星形領域では $F=\dd A$ と書ける。$A=\vb a^\flat-\varphi\dd t$ とすれば $B=\dd_3\vb a^\flat$、$E=-\dd_3\varphi-\partial_t\vb a^\flat$ である。$A$ に $\dd\chi$ を足しても $F$ は変わらない。

:::example[一様な軸方向磁場]
$B=b(t)\dd x\wedge\dd y$、$\vb a^\flat=b(t)(x\dd y-y\dd x)/2$、$\varphi=0$ とする。$E=-b'(t)(x\dd y-y\dd x)/2$ なので半径 $R$ の円周で $\int E=-\pi R^2 b'(t)$ となる。$R=0.1\,\mathrm m$、$b'=2\,\mathrm{T\,s^{-1}}$ なら起電力は $-0.02\pi\,\mathrm V$ である。この例は $\dd F=0$ を満たす場を構成しており、必要な電流はもう一方の方程式で決まる。
:::

真空で $\mathcal J=0$ のとき、$\dd F=0$ と $\dd\star F=0$ の両方が成り立つ。Lorentz 計量の二形式では $\star^2=-1$ なので、$F$ を $F_\vartheta=\cos\vartheta\, F+\sin\vartheta\,\star F$ に置き換えても両式を満たす。実際、$\dd F_\vartheta$ と $\dd\star F_\vartheta$ はそれぞれ二つの零の線形結合である。例えば $\vartheta=\pi/2$ のとき、空間二形式は $B_\vartheta=\star_3 E/c$、電場一形式は $E_\vartheta=-c\star_3 B$ となる。ここで $c$ は電場と磁束密度の単位を換算する。電荷や電流があると $\dd\star F$ が零でなくなり、この変換は一般に $\dd F_\vartheta=0$ を保たない。

### 電磁場の横波と保存則

真空で電荷と電流を零にし、波を $z$ 方向へ進める。$\xi=z-ct$ の関数なら $E_x=f(\xi)$ と $B_y=f(\xi)/c$ の組が Faraday の式と Ampère–Maxwell の式を満たす。例えば $\partial_z E_x=f^\prime$、$-\partial_tB_y=f^\prime$ となる。もう一方でも $-\partial_zB_y=-f^\prime/c$ と $\partial_t E_x/c^2=-f^\prime/c$ が一致する。源のある式を使う場合、$\dd^2=0$ から電荷の連続の式も従う。

:::example[真空の単色平面波（1）]
$\xi=z-ct$ とし、位置と時間を数値で無次元化して $E_x=\cos(m\xi)$ を調べる。基準の値で計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

電場と磁場の比は $c$ である。

$$
\begin{aligned}
E_x(\xi)\big|_{m=1} &= \cos{\left(\xi \right)}
\end{aligned}
$$

$\partial_z\xi=1$ なので、$z$ 微分は $\xi$ 微分に等しい。

$$
\begin{aligned}
\partial_zE_x\big|_{m=1} &= - \sin{\left(\xi \right)}
\end{aligned}
$$

電場の二乗と $c^2$ を掛けた磁場の二乗は等しい。エネルギー密度を $\varepsilon_0$ で割るとこの二乗が残る。

$$
\begin{aligned}
\mathcal E/\varepsilon_0\big|_{m=1} &= \cos^{2}{\left(\xi \right)}
\end{aligned}
$$

観測点 $\xi$=$\frac{\pi}{4}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
E_x(\xi)\big|_{\xi=\frac{\pi}{4}} &= \frac{\sqrt{2}}{2} \\
\frac{\partial^1 (E_x(\xi))}{\partial \xi^1} &= - \sin{\left(\xi \right)} \\
\left.\frac{\partial^1 (E_x(\xi))}{\partial \xi^1}\right|_{\xi=\frac{\pi}{4}} &= - \frac{\sqrt{2}}{2}
\end{aligned}
$$

観測点 $\xi$=$\frac{\pi}{4}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\partial_zE_x\big|_{\xi=\frac{\pi}{4}} &= - \frac{\sqrt{2}}{2} \\
\frac{\partial^1 (\partial_zE_x)}{\partial \xi^1} &= - \cos{\left(\xi \right)} \\
\left.\frac{\partial^1 (\partial_zE_x)}{\partial \xi^1}\right|_{\xi=\frac{\pi}{4}} &= - \frac{\sqrt{2}}{2}
\end{aligned}
$$

観測点 $\xi$=$\frac{\pi}{4}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\mathcal E/\varepsilon_0\big|_{\xi=\frac{\pi}{4}} &= \frac{1}{2} \\
\frac{\partial^1 (\mathcal E/\varepsilon_0)}{\partial \xi^1} &= - \sin{\left(2 \xi \right)} \\
\left.\frac{\partial^1 (\mathcal E/\varepsilon_0)}{\partial \xi^1}\right|_{\xi=\frac{\pi}{4}} &= -1
\end{aligned}
$$
:::

:::example[真空の単色平面波（2）]
$\xi=z-ct$ とし、位置と時間を数値で無次元化して $E_x=\cos(m\xi)$ を調べる。値を変えて依存性を調べる。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

電場と磁場の比は $c$ である。

$$
\begin{aligned}
E_x(\xi)\big|_{m=2} &= \cos{\left(2 \xi \right)}
\end{aligned}
$$

$\partial_z\xi=1$ なので、$z$ 微分は $\xi$ 微分に等しい。

$$
\begin{aligned}
\partial_zE_x\big|_{m=2} &= - 2 \sin{\left(2 \xi \right)}
\end{aligned}
$$

電場の二乗と $c^2$ を掛けた磁場の二乗は等しい。エネルギー密度を $\varepsilon_0$ で割るとこの二乗が残る。

$$
\begin{aligned}
\mathcal E/\varepsilon_0\big|_{m=2} &= \cos^{2}{\left(2 \xi \right)}
\end{aligned}
$$

観測点 $\xi$=$\frac{\pi}{4}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
E_x(\xi)\big|_{\xi=\frac{\pi}{4}} &= 0 \\
\frac{\partial^1 (E_x(\xi))}{\partial \xi^1} &= - 2 \sin{\left(2 \xi \right)} \\
\left.\frac{\partial^1 (E_x(\xi))}{\partial \xi^1}\right|_{\xi=\frac{\pi}{4}} &= -2
\end{aligned}
$$

観測点 $\xi$=$\frac{\pi}{4}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\partial_zE_x\big|_{\xi=\frac{\pi}{4}} &= -2 \\
\frac{\partial^1 (\partial_zE_x)}{\partial \xi^1} &= - 4 \cos{\left(2 \xi \right)} \\
\left.\frac{\partial^1 (\partial_zE_x)}{\partial \xi^1}\right|_{\xi=\frac{\pi}{4}} &= 0
\end{aligned}
$$

観測点 $\xi$=$\frac{\pi}{4}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\mathcal E/\varepsilon_0\big|_{\xi=\frac{\pi}{4}} &= 0 \\
\frac{\partial^1 (\mathcal E/\varepsilon_0)}{\partial \xi^1} &= - 2 \sin{\left(4 \xi \right)} \\
\left.\frac{\partial^1 (\mathcal E/\varepsilon_0)}{\partial \xi^1}\right|_{\xi=\frac{\pi}{4}} &= 0
\end{aligned}
$$
:::

:::example[真空の単色平面波（3）]
$\xi=z-ct$ とし、位置と時間を数値で無次元化して $E_x=\cos(m\xi)$ を調べる。さらに値を変えて比較する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

電場と磁場の比は $c$ である。

$$
\begin{aligned}
E_x(\xi)\big|_{m=3} &= \cos{\left(3 \xi \right)}
\end{aligned}
$$

$\partial_z\xi=1$ なので、$z$ 微分は $\xi$ 微分に等しい。

$$
\begin{aligned}
\partial_zE_x\big|_{m=3} &= - 3 \sin{\left(3 \xi \right)}
\end{aligned}
$$

電場の二乗と $c^2$ を掛けた磁場の二乗は等しい。エネルギー密度を $\varepsilon_0$ で割るとこの二乗が残る。

$$
\begin{aligned}
\mathcal E/\varepsilon_0\big|_{m=3} &= \cos^{2}{\left(3 \xi \right)}
\end{aligned}
$$

観測点 $\xi$=$\frac{\pi}{4}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
E_x(\xi)\big|_{\xi=\frac{\pi}{4}} &= - \frac{\sqrt{2}}{2} \\
\frac{\partial^1 (E_x(\xi))}{\partial \xi^1} &= - 3 \sin{\left(3 \xi \right)} \\
\left.\frac{\partial^1 (E_x(\xi))}{\partial \xi^1}\right|_{\xi=\frac{\pi}{4}} &= - \frac{3 \sqrt{2}}{2}
\end{aligned}
$$

観測点 $\xi$=$\frac{\pi}{4}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\partial_zE_x\big|_{\xi=\frac{\pi}{4}} &= - \frac{3 \sqrt{2}}{2} \\
\frac{\partial^1 (\partial_zE_x)}{\partial \xi^1} &= - 9 \cos{\left(3 \xi \right)} \\
\left.\frac{\partial^1 (\partial_zE_x)}{\partial \xi^1}\right|_{\xi=\frac{\pi}{4}} &= \frac{9 \sqrt{2}}{2}
\end{aligned}
$$

観測点 $\xi$=$\frac{\pi}{4}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\mathcal E/\varepsilon_0\big|_{\xi=\frac{\pi}{4}} &= \frac{1}{2} \\
\frac{\partial^1 (\mathcal E/\varepsilon_0)}{\partial \xi^1} &= - 3 \sin{\left(6 \xi \right)} \\
\left.\frac{\partial^1 (\mathcal E/\varepsilon_0)}{\partial \xi^1}\right|_{\xi=\frac{\pi}{4}} &= 3
\end{aligned}
$$
:::

:::example[真空の単色平面波（4）]
$\xi=z-ct$ とし、位置と時間を数値で無次元化して $E_x=\cos(m\xi)$ を調べる。最後にこの値でも計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=4
$$

電場と磁場の比は $c$ である。

$$
\begin{aligned}
E_x(\xi)\big|_{m=4} &= \cos{\left(4 \xi \right)}
\end{aligned}
$$

$\partial_z\xi=1$ なので、$z$ 微分は $\xi$ 微分に等しい。

$$
\begin{aligned}
\partial_zE_x\big|_{m=4} &= - 4 \sin{\left(4 \xi \right)}
\end{aligned}
$$

電場の二乗と $c^2$ を掛けた磁場の二乗は等しい。エネルギー密度を $\varepsilon_0$ で割るとこの二乗が残る。

$$
\begin{aligned}
\mathcal E/\varepsilon_0\big|_{m=4} &= \cos^{2}{\left(4 \xi \right)}
\end{aligned}
$$

観測点 $\xi$=$\frac{\pi}{4}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
E_x(\xi)\big|_{\xi=\frac{\pi}{4}} &= -1 \\
\frac{\partial^1 (E_x(\xi))}{\partial \xi^1} &= - 4 \sin{\left(4 \xi \right)} \\
\left.\frac{\partial^1 (E_x(\xi))}{\partial \xi^1}\right|_{\xi=\frac{\pi}{4}} &= 0
\end{aligned}
$$

観測点 $\xi$=$\frac{\pi}{4}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\partial_zE_x\big|_{\xi=\frac{\pi}{4}} &= 0 \\
\frac{\partial^1 (\partial_zE_x)}{\partial \xi^1} &= - 16 \cos{\left(4 \xi \right)} \\
\left.\frac{\partial^1 (\partial_zE_x)}{\partial \xi^1}\right|_{\xi=\frac{\pi}{4}} &= 16
\end{aligned}
$$

観測点 $\xi$=$\frac{\pi}{4}$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
\mathcal E/\varepsilon_0\big|_{\xi=\frac{\pi}{4}} &= 1 \\
\frac{\partial^1 (\mathcal E/\varepsilon_0)}{\partial \xi^1} &= - 8 \sin{\left(4 \xi \right)} \cos{\left(4 \xi \right)} \\
\left.\frac{\partial^1 (\mathcal E/\varepsilon_0)}{\partial \xi^1}\right|_{\xi=\frac{\pi}{4}} &= 0
\end{aligned}
$$
:::

## 熱力学と微分形式（Carathéodory の定理）

### 断熱面と積分因子

平衡状態を座標で表し、可逆過程の微小な熱を一形式 $\omega$ と書く。準静的な体積仕事だけがある単純圧縮系では $\omega=\dd U+p\dd V$ である。可逆断熱過程に沿う接ベクトル $X$ は $\omega(X)=0$ を満たす。各点でこの条件を満たす方向が、近くの状態を通る面に接しているかを調べる。

:::theorem[Carathéodory の局所条件]{#thm-mp-caratheodory}
非零で滑らかな一形式 $\omega$ が定める $\Ker\omega$ に沿う区分滑らかな曲線だけを許す。各点 $p$ について、十分小さい任意の近傍 $U$ の内部に留まる許容曲線で $p$ から到達できる点の集合が、$p$ の近傍を含まないという条件を考える。この条件は $\omega\wedge\dd\omega=0$ と同値であり、そのとき各点の近くで $\omega=T\dd S$ と書ける。$T$ は零でなく、$\dd S\ne0$ である。
:::

:::proof
座標の最後を $z$ とし、係数を割って $\omega/T_0=\dd z-\sum_{i=1}^{n-1}a_i\dd x^i$ と書く。$\Ker\omega$ は $X_i=\partial_i+a_i\partial_z$ で張られる。$[X_i,X_j]$ は縦方向であり、外微分の定義から $\omega([X_i,X_j])=-\dd\omega(X_i,X_j)$ である。$\omega\wedge\dd\omega=0$ とこれらの括弧がすべて零であることは同値である。実際、$\omega(Y)=1$ となる $Y$ を選ぶと $(\omega\wedge\dd\omega)(Y,X_i,X_j)=\dd\omega(X_i,X_j)$ であり、すべて横方向に入れると零である。

ある点で $[X_i,X_j]\ne0$ とする。$X_i$ の流れを $\Phi^i_s$ と書き、$\Phi^i_{-s}\circ\Phi^j_t\circ\Phi^i_s$ を $t$ で零の周りに動かす。この移動は許された曲線の連結であり、$t=0$ で元の点を返す。接ベクトルは $X_j+s[X_i,X_j]+O(s^2)$ である。これは引き戻したベクトル場の $s=0$ での微分を連鎖律で計算して得る。固定した小さい $s\ne0$ に対して縦方向の成分が横方向の張る空間から外れる。各 $X_k$ による短い流れを続けると終点写像の微分が全階数になる。逆写像定理により近傍がすべて到達可能になり、仮定に反する。

逆にすべての括弧が零なら、$X_i$ の流れは互いに相手のベクトル場を保つ。実際、流れで引き戻した場の時間微分は括弧の引き戻しであり零である。$z$ 軸上の初期点をパラメータ $S$ で指定し、各 $X_i$ の流れの時間を残りの座標に取る。初期点では接ベクトルが独立なので逆写像定理が使える。流れが交換するため $X_iS=0$ となり $\omega=T\dd S$ を得る。この座標近傍の内部に留まる許容曲線では $S$ が一定である。従って、同じ近傍内で $S$ の値が異なる任意に近い点には到達できない。
:::

:::example[理想気体の積分因子]
物質量一定の理想気体で $pV=nRT$、$U=C_VT$、$C_V$ 一定とする。$T>0$、$V>0$ で $\omega/T=C_V\dd T/T+nR\dd V/V$ なので $S-S_0=C_V\log(T/T_0)+nR\log(V/V_0)$ と取れる。体積を一定温度で二倍にすると $S$ は $nR\log2$（$\mathrm{J\,K^{-1}}$）増加する。局所的な因子の存在だけから温度の絶対尺度は定まらない。
:::

:::example[可積分性を失う三変数の形式]
$\omega=\dd z-x\dd y$ とすると $\omega\wedge\dd\omega=-\dd x\wedge\dd y\wedge\dd z$ である。$z$ の変化は $\dd z=x\dd y$ に従う。$(x,y)$ 平面で一辺 $a,b$ の長方形を反時計回りに一周すると、$z$ が $ab$ だけ変わる。平面への射影が同じ点に戻っても元の状態とは異なる。
:::

### 積分因子の自由度

$\omega=T\dd S$ が局所的に得られたとし、$S$ は近傍の一つの座標にできるとする。別の表示 $\omega=\widetilde T\dd\widetilde S$ があれば、$\widetilde S$ の微分は $\dd S$ の倍数なので、他の座標と $S$ の直積に縮めた近傍で $\widetilde S=f(S)$ と書ける。従って $\widetilde T=T/f'(S)$ である。逆に $f'\ne0$ ならこの置き換えで同じ $\omega$ が得られる。積分因子の存在から直ちに単位まで定まらない理由である。

例えば $T=1$、$S=x^2+y$ に対して $\widetilde S=e^{2S}$ を取れば $\widetilde T=\tfrac12e^{-2S}$ になる。原点の近くではどちらの係数も有限で、断熱面 $x^2+y=\mathrm{const.}$ は変わらない。変数に物理単位を持たせた場合は、指数の引数を無次元にする基準値も必要となる。

### 積分因子の非一意性と温度

$\omega=T\,\dd S$ が成り立つとき、$S$ を単調な関数 $F(S)$ に取り直せる。$F^\prime\ne0$ なら $\omega=[T/F^\prime(S)]\,\dd F(S)$ である。従って局所可積分性だけからは積分因子の数値的な規格化は決まらない。温度の尺度とエントロピーの加法性を合わせるには、熱接触などの追加の物理条件が要る。非完全微分を二変数で扱う場合も、原始関数にする前に積分因子を選ぶ。

:::example[可逆な理想気体の断熱線（1）]
物質量を $1\,\mathrm{mol}$ とし、モル定積熱容量と気体定数の比を $m=C_{V,\mathrm{mol}}/R$ とする。温度と体積はそれぞれ初期値で割って無次元化する。$T V^{1/m}$ が一定の断熱線を調べる。計算式では $u=V$ とする。基準の値で計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=1
$$

初期点 $(T,V)=(1,1)$ を通る断熱線を求める。

$$
\begin{aligned}
T(V)\big|_{m=1} &= \frac{1}{u}
\end{aligned}
$$

定積熱容量の項と体積仕事の項を $R$ で割り、温度の変化を代入する。

$$
\begin{aligned}
m\,\dd T/\dd V\big|_{m=1} &= - \frac{1}{u^{2}}
\end{aligned}
$$

残った体積項は上の値と反対である。

$$
\begin{aligned}
T/V\big|_{m=1} &= \frac{1}{u^{2}}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
T(V)\big|_{u=2} &= \frac{1}{2} \\
\frac{\partial^1 (T(V))}{\partial u^1} &= - \frac{1}{u^{2}} \\
\left.\frac{\partial^1 (T(V))}{\partial u^1}\right|_{u=2} &= - \frac{1}{4} \\
\frac{\partial^2 (T(V))}{\partial u^2} &= \frac{2}{u^{3}} \\
\left.\frac{\partial^2 (T(V))}{\partial u^2}\right|_{u=2} &= \frac{1}{4}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
m\,\dd T/\dd V\big|_{u=2} &= - \frac{1}{4} \\
\frac{\partial^1 (m\,\dd T/\dd V)}{\partial u^1} &= \frac{2}{u^{3}} \\
\left.\frac{\partial^1 (m\,\dd T/\dd V)}{\partial u^1}\right|_{u=2} &= \frac{1}{4} \\
\frac{\partial^2 (m\,\dd T/\dd V)}{\partial u^2} &= - \frac{6}{u^{4}} \\
\left.\frac{\partial^2 (m\,\dd T/\dd V)}{\partial u^2}\right|_{u=2} &= - \frac{3}{8}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
T/V\big|_{u=2} &= \frac{1}{4} \\
\frac{\partial^1 (T/V)}{\partial u^1} &= - \frac{2}{u^{3}} \\
\left.\frac{\partial^1 (T/V)}{\partial u^1}\right|_{u=2} &= - \frac{1}{4} \\
\frac{\partial^2 (T/V)}{\partial u^2} &= \frac{6}{u^{4}} \\
\left.\frac{\partial^2 (T/V)}{\partial u^2}\right|_{u=2} &= \frac{3}{8}
\end{aligned}
$$
:::

:::example[可逆な理想気体の断熱線（2）]
物質量を $1\,\mathrm{mol}$ とし、モル定積熱容量と気体定数の比を $m=C_{V,\mathrm{mol}}/R$ とする。温度と体積はそれぞれ初期値で割って無次元化する。$T V^{1/m}$ が一定の断熱線を調べる。計算式では $u=V$ とする。値を変えて依存性を調べる。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=2
$$

初期点 $(T,V)=(1,1)$ を通る断熱線を求める。

$$
\begin{aligned}
T(V)\big|_{m=2} &= \frac{1}{\sqrt{u}}
\end{aligned}
$$

定積熱容量の項と体積仕事の項を $R$ で割り、温度の変化を代入する。

$$
\begin{aligned}
m\,\dd T/\dd V\big|_{m=2} &= - \frac{1}{u^{\frac{3}{2}}}
\end{aligned}
$$

残った体積項は上の値と反対である。

$$
\begin{aligned}
T/V\big|_{m=2} &= \frac{1}{u^{\frac{3}{2}}}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
T(V)\big|_{u=2} &= \frac{\sqrt{2}}{2} \\
\frac{\partial^1 (T(V))}{\partial u^1} &= - \frac{1}{2 u^{\frac{3}{2}}} \\
\left.\frac{\partial^1 (T(V))}{\partial u^1}\right|_{u=2} &= - \frac{\sqrt{2}}{8} \\
\frac{\partial^2 (T(V))}{\partial u^2} &= \frac{3}{4 u^{\frac{5}{2}}} \\
\left.\frac{\partial^2 (T(V))}{\partial u^2}\right|_{u=2} &= \frac{3 \sqrt{2}}{32}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
m\,\dd T/\dd V\big|_{u=2} &= - \frac{\sqrt{2}}{4} \\
\frac{\partial^1 (m\,\dd T/\dd V)}{\partial u^1} &= \frac{3}{2 u^{\frac{5}{2}}} \\
\left.\frac{\partial^1 (m\,\dd T/\dd V)}{\partial u^1}\right|_{u=2} &= \frac{3 \sqrt{2}}{16} \\
\frac{\partial^2 (m\,\dd T/\dd V)}{\partial u^2} &= - \frac{15}{4 u^{\frac{7}{2}}} \\
\left.\frac{\partial^2 (m\,\dd T/\dd V)}{\partial u^2}\right|_{u=2} &= - \frac{15 \sqrt{2}}{64}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
T/V\big|_{u=2} &= \frac{\sqrt{2}}{4} \\
\frac{\partial^1 (T/V)}{\partial u^1} &= - \frac{3}{2 u^{\frac{5}{2}}} \\
\left.\frac{\partial^1 (T/V)}{\partial u^1}\right|_{u=2} &= - \frac{3 \sqrt{2}}{16} \\
\frac{\partial^2 (T/V)}{\partial u^2} &= \frac{15}{4 u^{\frac{7}{2}}} \\
\left.\frac{\partial^2 (T/V)}{\partial u^2}\right|_{u=2} &= \frac{15 \sqrt{2}}{64}
\end{aligned}
$$
:::

:::example[可逆な理想気体の断熱線（3）]
物質量を $1\,\mathrm{mol}$ とし、モル定積熱容量と気体定数の比を $m=C_{V,\mathrm{mol}}/R$ とする。温度と体積はそれぞれ初期値で割って無次元化する。$T V^{1/m}$ が一定の断熱線を調べる。計算式では $u=V$ とする。さらに値を変えて比較する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=3
$$

初期点 $(T,V)=(1,1)$ を通る断熱線を求める。

$$
\begin{aligned}
T(V)\big|_{m=3} &= \frac{1}{\sqrt[3]{u}}
\end{aligned}
$$

定積熱容量の項と体積仕事の項を $R$ で割り、温度の変化を代入する。

$$
\begin{aligned}
m\,\dd T/\dd V\big|_{m=3} &= - \frac{1}{u^{\frac{4}{3}}}
\end{aligned}
$$

残った体積項は上の値と反対である。

$$
\begin{aligned}
T/V\big|_{m=3} &= \frac{1}{u^{\frac{4}{3}}}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
T(V)\big|_{u=2} &= \frac{2^{\frac{2}{3}}}{2} \\
\frac{\partial^1 (T(V))}{\partial u^1} &= - \frac{1}{3 u^{\frac{4}{3}}} \\
\left.\frac{\partial^1 (T(V))}{\partial u^1}\right|_{u=2} &= - \frac{2^{\frac{2}{3}}}{12} \\
\frac{\partial^2 (T(V))}{\partial u^2} &= \frac{4}{9 u^{\frac{7}{3}}} \\
\left.\frac{\partial^2 (T(V))}{\partial u^2}\right|_{u=2} &= \frac{2^{\frac{2}{3}}}{18}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
m\,\dd T/\dd V\big|_{u=2} &= - \frac{2^{\frac{2}{3}}}{4} \\
\frac{\partial^1 (m\,\dd T/\dd V)}{\partial u^1} &= \frac{4}{3 u^{\frac{7}{3}}} \\
\left.\frac{\partial^1 (m\,\dd T/\dd V)}{\partial u^1}\right|_{u=2} &= \frac{2^{\frac{2}{3}}}{6} \\
\frac{\partial^2 (m\,\dd T/\dd V)}{\partial u^2} &= - \frac{28}{9 u^{\frac{10}{3}}} \\
\left.\frac{\partial^2 (m\,\dd T/\dd V)}{\partial u^2}\right|_{u=2} &= - \frac{7 \cdot 2^{\frac{2}{3}}}{36}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
T/V\big|_{u=2} &= \frac{2^{\frac{2}{3}}}{4} \\
\frac{\partial^1 (T/V)}{\partial u^1} &= - \frac{4}{3 u^{\frac{7}{3}}} \\
\left.\frac{\partial^1 (T/V)}{\partial u^1}\right|_{u=2} &= - \frac{2^{\frac{2}{3}}}{6} \\
\frac{\partial^2 (T/V)}{\partial u^2} &= \frac{28}{9 u^{\frac{10}{3}}} \\
\left.\frac{\partial^2 (T/V)}{\partial u^2}\right|_{u=2} &= \frac{7 \cdot 2^{\frac{2}{3}}}{36}
\end{aligned}
$$
:::

:::example[可逆な理想気体の断熱線（4）]
物質量を $1\,\mathrm{mol}$ とし、モル定積熱容量と気体定数の比を $m=C_{V,\mathrm{mol}}/R$ とする。温度と体積はそれぞれ初期値で割って無次元化する。$T V^{1/m}$ が一定の断熱線を調べる。計算式では $u=V$ とする。最後にこの値でも計算する。

まずパラメータを固定する。ここでは次の値で計算する。

$$
m=4
$$

初期点 $(T,V)=(1,1)$ を通る断熱線を求める。

$$
\begin{aligned}
T(V)\big|_{m=4} &= \frac{1}{\sqrt[4]{u}}
\end{aligned}
$$

定積熱容量の項と体積仕事の項を $R$ で割り、温度の変化を代入する。

$$
\begin{aligned}
m\,\dd T/\dd V\big|_{m=4} &= - \frac{1}{u^{\frac{5}{4}}}
\end{aligned}
$$

残った体積項は上の値と反対である。

$$
\begin{aligned}
T/V\big|_{m=4} &= \frac{1}{u^{\frac{5}{4}}}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
T(V)\big|_{u=2} &= \frac{2^{\frac{3}{4}}}{2} \\
\frac{\partial^1 (T(V))}{\partial u^1} &= - \frac{1}{4 u^{\frac{5}{4}}} \\
\left.\frac{\partial^1 (T(V))}{\partial u^1}\right|_{u=2} &= - \frac{2^{\frac{3}{4}}}{16} \\
\frac{\partial^2 (T(V))}{\partial u^2} &= \frac{5}{16 u^{\frac{9}{4}}} \\
\left.\frac{\partial^2 (T(V))}{\partial u^2}\right|_{u=2} &= \frac{5 \cdot 2^{\frac{3}{4}}}{128}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
m\,\dd T/\dd V\big|_{u=2} &= - \frac{2^{\frac{3}{4}}}{4} \\
\frac{\partial^1 (m\,\dd T/\dd V)}{\partial u^1} &= \frac{5}{4 u^{\frac{9}{4}}} \\
\left.\frac{\partial^1 (m\,\dd T/\dd V)}{\partial u^1}\right|_{u=2} &= \frac{5 \cdot 2^{\frac{3}{4}}}{32} \\
\frac{\partial^2 (m\,\dd T/\dd V)}{\partial u^2} &= - \frac{45}{16 u^{\frac{13}{4}}} \\
\left.\frac{\partial^2 (m\,\dd T/\dd V)}{\partial u^2}\right|_{u=2} &= - \frac{45 \cdot 2^{\frac{3}{4}}}{256}
\end{aligned}
$$

観測点 $u$=$2$ で値を照合する。微分も必要なら、パラメータを固定してから求める。

$$
\begin{aligned}
T/V\big|_{u=2} &= \frac{2^{\frac{3}{4}}}{4} \\
\frac{\partial^1 (T/V)}{\partial u^1} &= - \frac{5}{4 u^{\frac{9}{4}}} \\
\left.\frac{\partial^1 (T/V)}{\partial u^1}\right|_{u=2} &= - \frac{5 \cdot 2^{\frac{3}{4}}}{32} \\
\frac{\partial^2 (T/V)}{\partial u^2} &= \frac{45}{16 u^{\frac{13}{4}}} \\
\left.\frac{\partial^2 (T/V)}{\partial u^2}\right|_{u=2} &= \frac{45 \cdot 2^{\frac{3}{4}}}{256}
\end{aligned}
$$
:::

## 演習問題

::::exercise{level=basic}
$\alpha=x\dd y$ と $\beta=y\dd z$ の外積とその外微分を求めよ。
:::solution
$\alpha\wedge\beta=xy\dd y\wedge\dd z$ である。係数の微分は $y\dd x+x\dd y$ であり、二つ目の項は $\dd y$ が重複して消える。従って $\dd(\alpha\wedge\beta)=y\dd x\wedge\dd y\wedge\dd z$ となる。
:::
::::

::::exercise{level=basic}
$\omega=2xy\dd x+x^2\dd y$ の原始関数を求めよ。
:::solution
$x$ 方向の係数を積分すると $f=x^2y+C(y)$ になる。$y$ 方向に微分すると $x^2+C'(y)$ だから $C'=0$ であり、$f=x^2y+C$ が解である。
:::
::::

::::exercise{level=basic}
$\R^3$ で $\star(2\dd x+\dd z)$ と $\star^2(2\dd x+\dd z)$ を求めよ。
:::solution
$\dd x$ の補集合は $(y,z)$、$\dd z$ の補集合は $(x,y)$ であり、どちらも巡回順なので第一の結果は $2\dd y\wedge\dd z+\dd x\wedge\dd y$ である。二回の符号は $(-1)^2=1$ で元に戻る。
:::
::::

::::exercise{level=basic}
$E=x^2\dd y$、$B=0$ で $\dd F$ を求めよ。
:::solution
$F=x^2\dd y\wedge\dd t$ なので、$\dd F=2x\dd x\wedge\dd y\wedge\dd t$ である。$x\ne0$ では零でなく、静的な電場と零磁場の組は Faraday の法則を満たさない。
:::
::::

::::exercise{level=advanced}
$\R^3\setminus\{0\}$ で $\alpha=\iota_{\vb r/r^3}\mathrm{vol}$ が閉で完全でないことを示せ。
:::solution
直交座標で発散は $3r^{-3}-3r^2r^{-5}=0$ であり、:ref[prop-mp-vector-forms] より閉である。半径 $R$ の球面では $\vb r/r^3$ の法線成分は $1/R^2$ なので積分は $4\pi$ である。もし $\alpha=\dd\beta$ なら球面へ引き戻した形式も完全であり、境界のない球面上の積分は :ref[thm-mp-general-stokes] により零である。この矛盾から結論を得る。
:::
::::

::::exercise{level=advanced}
$\omega=\dd z+2x\dd x+\dd y$ に積分因子があることを示し、断熱面を求めよ。
:::solution
$\dd(z+x^2+y)=\omega$ であるから、$T=1$、$S=z+x^2+y$ と取れる。従って断熱面は $z+x^2+y$ が一定の面であり、たとえば原点を通るものは $z=-x^2-y$ である。
:::
::::

::::exercise{level=advanced}
$\omega=\dd z+y\dd x$ の可積分性を判定し、平面の長方形 $(0,0),(2,0),(2,3),(0,3)$ に沿う状態の変化を求めよ。
:::solution
$\dd\omega=\dd y\wedge\dd x$ だから $\omega\wedge\dd\omega=-\dd x\wedge\dd y\wedge\dd z\ne0$ である。順に歩くと $\dd z=-y\dd x$ は第三辺だけで寄与し、$-3(0-2)=6$ となる。従って積分因子は局所的にも存在しない。
:::
::::

::::exercise{level=advanced}
滑らかな $\omega$ に対し $\Delta_H\dd\omega=\dd\Delta_H\omega$ を示せ。さらに直交座標で $\omega=x^3\dd y$ を代入して確かめよ。
:::solution
両辺を展開すると $\dd^2=0$ からいずれも $\dd\delta\dd\omega$ になる。左辺は $\Delta_H(3x^2\dd x\wedge\dd y)=-6\dd x\wedge\dd y$、右辺は $\dd(-6x\dd y)=-6\dd x\wedge\dd y$ である。
:::
::::

## 確認テスト

::::quiz{answer="-1" tol="1e-9"}
$\R^3$ で $\star(\dd x\wedge\dd z)=a\dd y$ の係数 $a$ を求めよ。

:::explanation
$\dd x\wedge\dd z\wedge\dd y$ は $\dd x\wedge\dd y\wedge\dd z$ と反対向きである。内積のノルムは $1$ なので $a=-1$ が必要である。
:::
::::

::::quiz
閉形式が完全であると保証される条件はどれか。

- [x] 星形開領域上の正の次数の閉形式
- [ ] 穴のある任意の領域上の閉一形式
- [ ] 閉零形式

:::explanation
正次数と星形性から :ref[thm-mp-poincare] が使える。穴があれば円周積分の非零な閉一形式が存在し得る。零形式ではそもそもこの複体に $-1$ 形式はなく、零でない定数は完全でない。
:::
::::

::::quiz{answer="6" tol="1e-9"}
$f=-x^2-y^2-z^2$ について $\R^3$ で $\Delta_H f$ を求めよ。

:::explanation
三方向の二階微分がそれぞれ $-2$ なので、$\Delta_H=-\Delta$ は和の符号を変え、$6$ となる。
:::
::::

::::quiz
$\dd G=\mathcal J$ にもう一度外微分を作用させるとどの物理法則を得るか。

- [ ] 磁束の単位の定義
- [ ] 真空中の光速の値
- [x] 電荷の局所保存

:::explanation
$\dd^2G=0$ なので $\dd\mathcal J=0$ である。時間を含む四形式の係数を読むと $\partial_t\rho+\divg\vb J=0$ となる。
:::
::::

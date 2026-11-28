---
title: 双曲型方程式
description: 二階波動方程式のエネルギー弱解と、対称双曲型系の有限伝播および Fourier 表示を扱う。
notation:
  - symbol: 'E(t)'
    meaning: 運動項と変形項を合わせたエネルギー
  - symbol: 'A^0'
    meaning: 対称双曲型系の時間微分の係数
  - symbol: 'B(\xi)'
    meaning: 一階定数係数系の空間記号
---

## この章の目標

- 二階双曲型方程式の弱形式とエネルギー空間を定める。
- 弱解の存在、一意性、連続依存性を証明する。
- 移動球上のエネルギーで有限伝播速度を求める。
- 対称双曲型系と一般の定数係数系で Fourier 乗数を比較する。

## 2 階双曲型方程式の弱解

[波動方程式](/math/pde/wave/) では基本解を用いた。係数が場所で変わる場合でも、変位の勾配と速度に関するエネルギーは残る。そこで明示的な解の表示を使わずに方程式を定義する。

:::definition[二階のエネルギー弱解]{#def-pde-hyperbolic-weak}
有界領域 $\Omega$ で $V=H_0^1(\Omega)$、$H=L^2(\Omega)$ とする。$a(u,v)=\int_\Omega A(x)Du\cdot Dv$ は有界かつ強圧的で、$A$ は実対称とする。$g\in V$、$h\in H$、$f\in L^1(0,T;H)$ に対し、$u\in C([0,T];V)\cap C^1([0,T];H)$ が

$$
-\int_0^T(u',v')_H\,\dd t+\int_0^Ta(u,v)\,\dd t
=\int_0^T(f,v)_H\,\dd t+(h,v(0))_H
\label{eq-pde-hyperbolic-weak}
$$

を $v\in C_c^1([0,T);V)$ に対して満たし、$u(0)=g$ であるとき $u''-\divg(A Du)=f$ の :term[エネルギー弱解]{en="energy weak solution" yomi="えねるぎーじゃっかい"} と呼ぶ。
:::

空間部分積分で境界値 $0$ が試験空間に移り、時間部分積分では初速度だけが端点項に残る。個々の時刻について二階時間微分を $H$ で取る必要はない。方程式から $u''\in L^1(0,T;V^*)$ である。

:::example[二階微分を持たない初期変位]
$(0,\pi)$ の Dirichlet スペクトルで $g_j=1/j^2$、$h=0$ とする。$\sum_j j^2|g_j|^2=\sum_jj^{-2}<\infty$ なので初期変位は $V$ に属する。一方、$\sum_jj^4|g_j|^2=\infty$ なので $-\Delta g$ は $L^2$ に入らない。このデータを扱うには $u''$ を双対空間で読む必要がある。
:::

:::example[弱形式の初速度の端点項]
$(0,\pi)$ で $u=\sin t\sin x$ とする。$v(t,x)=\eta(t)\sin x$ なら :eqref[eq-pde-hyperbolic-weak] の左辺は $(\pi/2)\int(-\cos t\,\eta'+\sin t\,\eta)$ である。時間部分積分で内部項が消え、値は $(\pi/2)\eta(0)$ である。これは初速度 $h=\sin x$ と $v(0)$ の内積である。
:::

### 時間に関する弱式の二つの初期値

$u''+Au=f$ とし、$v\in V$ と時間試験関数 $\eta$ を取る。$\eta(T)=\eta'(T)=0$ とする。二回時間積分すれば
$$
\int_0^T[(u,v)_H\eta''+a(u,v)\eta]\dd t
=\int_0^T(f,v)_H\eta\dd t+(h,v)_H\eta(0)-(g,v)_H\eta'(0).
$$
端点の符号は $[(u',v)_H\eta-(u,v)_H\eta']_0^T$ を左辺から右辺へ移したものである。$\eta(0)$ と $\eta'(0)$ は独立に選べるので、変位の跡と速度の跡を同じ初期項にまとめることはできない。逆に $u\in C(V)\cap C^1(H)$ がこの式を満たせば、まず時刻内部に台を持つ試験で分布方程式を得る。境界項を比較して両方の初期値を回復する。

:::example[弱式の二つの端点項に入る多項式]
一次元で $P=4\partial_{xx}$、$g=x^3+1$、$h=2x^4$ とする。局所弱式で時刻零のどの係数を使うかを計算する。

時間二階の方程式なので偶数階と奇数階を分ける。$g$ からは $g,Pg,P^2g,\ldots$、$h$ からは $h,Ph,P^2h,\ldots$ を得る。空間の多項式に定数係数の二階作用素を掛けると次数が少なくとも二つ下がる。従ってこの場合には両方の列が有限回で零になり、級数の収束を別に仮定する必要はない。

$$
\begin{aligned}
P^{0}g&=x^{3} + 1\\
P^{1}g&=24 x\\
P^{0}h&=2 x^{4}\\
P^{1}h&=96 x^{2}\\
P^{2}h&=768\\
P^{3}g&=0
\end{aligned}
$$

$P^j g$ を $t^{2j}/(2j)!$、$P^j h$ を $t^{2j+1}/(2j+1)!$ と掛けて足す。時間で二回微分すると各項の添字が一つ下がり、空間作用素を一回掛けた列と一致する。時刻零の変位と速度は別々に確認する。

$$
\begin{aligned}
u(t,x)&=\frac{32 t^{5}}{5} + 16 t^{3} x^{2} + 12 t^{2} x + 2 t x^{4} + x^{3} + 1\\
u_{tt}&=\frac{\partial^{2}}{\partial t^{2}} \left(\frac{32 t^{5}}{5} + 16 t^{3} x^{2} + 12 t^{2} x + 2 t x^{4} + x^{3} + 1\right)=8 \left(16 t^{3} + 12 t x^{2} + 3 x\right)\\
Pu&=8 \left(16 t^{3} + 12 t x^{2} + 3 x\right)\\
u(0,x)&=x^{3} + 1,\qquad t=0\\
u_t(0,x)&=2 x^{4},\qquad t=0
\end{aligned}
$$

空間の導関数を測ると、時間の多項式の各係数が局所でどう効くかが分かる。初期関数との差は単に平行移動で置き換えず、以下では導関数と時刻を代入して計算する。

$$
\begin{aligned}
\partial_{x} u&=\frac{\partial}{\partial x} \left(\frac{32 t^{5}}{5} + 16 t^{3} x^{2} + 12 t^{2} x + 2 t x^{4} + x^{3} + 1\right)=\left(4 t^{2} + x^{2}\right) \left(8 t x + 3\right)\\
g&=2,\qquad x=1\\
u&=\frac{41}{5},\qquad t=\frac{1}{2},x=1\\
u_t&=28,\qquad t=\frac{1}{2},x=1\\
h&=2,\qquad x=1
\end{aligned}
$$

この解を全空間の弱式で使うとき、空間試験関数はコンパクトな台を持つものに取る。$g$ の項は $-\eta^{\prime}(0)$、$h$ の項は $\eta(0)$ と組み合わされる。従って例えば $\eta(0)=0$、$\eta^{\prime}(0)=1$ なら、初期速度の項は消え、変位の項が負号で残る。
:::

:::example[有限次元の試験空間で弱形式を読む]
$(0,1)$ で $A=2$、$c=0$ とし、時刻一定の解 $u=x^2(1-x)$ を支える外力を置く。波動方程式でも定常解は同じ空間方程式を満たす。

剛性行列を作る前に、方程式と試験空間を決める。ここでは全ての基底が両端で消える。厳密解を含まない試験空間に制限すると、近似は割り当てた方向の弱形式を満たすように決まる。厳密解は誤差を測るときに使う。

$$
\begin{aligned}
A&=2\\
u&=x^{2} \left(1 - x\right)\\
f&=4 \left(3 x - 1\right)\\
\psi_1&=x \left(1 - x\right)\\
\psi_2&=x \left(1 - x\right) \left(2 x - 1\right)
\end{aligned}
$$

剛性の成分は $K_{ij}=\int A\psi^{\prime}_j\psi^{\prime}_i+\int c\psi_j\psi_i$、右辺は $F_i=\int f\psi_i$ である。基底を正規直交には取っていないので、内積の質量行列 $M_{ij}=\int\psi_j\psi_i$ も別に求める。

$$
\begin{aligned}
K&=\left[\begin{matrix}\frac{2}{3} & 0\\0 & \frac{2}{5}\end{matrix}\right]\\
M&=\left[\begin{matrix}\frac{1}{30} & 0\\0 & \frac{1}{210}\end{matrix}\right]\\
F&=\left[\begin{matrix}\frac{1}{3}\\\frac{1}{5}\end{matrix}\right]\\
\det(K-\lambda M)&=\frac{\lambda^{2}}{6300} - \frac{26 \lambda}{1575} + \frac{4}{15}=\frac{\left(\lambda - 84\right) \left(\lambda - 20\right)}{6300}
\end{aligned}
$$

連立方程式 $Kc=F$ を解き、元の関数に戻す。残差 $r=f-Lu_N$ を計算すると一般に零ではないが、選んだ試験関数に対する積分は全て零になる。PDE を点で満たすことと有限次元の弱形式を満たすことを区別する。

$$
\begin{aligned}
c&=\left[\begin{matrix}\frac{1}{2}\\\frac{1}{2}\end{matrix}\right]\\
u_N&=- x^{2} \left(x - 1\right)\\
e&=u-u_N=0\\
r&=0\\
(\int r\psi_i)_i&=\left[\begin{matrix}0\\0\end{matrix}\right]
\end{aligned}
$$

最後に誤差の二乗ノルムを二通りに計算する。近似空間へのエネルギー直交性により、誤差の勾配エネルギーは全体のエネルギーから近似の分を引いた値になる。ここでは二次形式そのものを用いており、エネルギー汎関数の二分の一は付けていない。

$$
\begin{aligned}
a(u,u)&=\frac{4}{15}\\
a(u_N,u_N)&=\frac{4}{15}\\
a(e,e)&=0\\
\norm e_2^2&=\int\limits_{0}^{1} 0\, dx=0
\end{aligned}
$$

時間を戻すと式は $M c^{\prime\prime}+Kc=F$ である。ここで求めた $K^{-1}F$ を初期変位に、零を初期速度に選ぶと係数は一定になる。初期変位または速度を変えると、ずれは同じ $M,K$ で定まる振動子になる。弱式の二つの初期項が必要なのは、このずれを決めるためでもある。
:::

:::example[二つのモードを持つ弱式の初期端点]
波動方程式の独立な二固有方向で $M=I$、剛性固有値を $1,9$ とし、初期変位 $(2,0)^\top$、初期速度 $(0,6)^\top$ を与える。

各行の質量係数と剛性係数を区別して、候補を連立方程式に代入する。$B$ の非対角成分は別の成分を結び付ける。従って個々の微分の値を求めるだけではなく、全成分を掛け合わせた結果を確認する。

$$
\begin{aligned}
M&=\left[\begin{matrix}1 & 0\\0 & 1\end{matrix}\right]\\
B&=\left[\begin{matrix}1 & 0\\0 & 9\end{matrix}\right]\\
U&=\left[\begin{matrix}2 \cos{\left(t \right)}\\2 \sin{\left(3 t \right)}\end{matrix}\right]\\
U_{tt}&=\left[\begin{matrix}- 2 \cos{\left(t \right)}\\- 18 \sin{\left(3 t \right)}\end{matrix}\right]\\
F&=\left[\begin{matrix}0\\0\end{matrix}\right]
\end{aligned}
$$

初期時刻へ戻して、先に指定したデータと比較する。次に正の時刻における成分の値を読む。成分が混ざっていても、初期データは方程式を積分する前の元の座標で指定される。

$$
\begin{aligned}
U(0)&=\left[\begin{matrix}2\\0\end{matrix}\right],\qquad t=0\\
U_t(0)&=\left[\begin{matrix}0\\6\end{matrix}\right],\qquad t=0\\
U&=\left[\begin{matrix}\sqrt{3}\\2\end{matrix}\right],\qquad t=\frac{\pi}{6}\\
U_t&=\left[\begin{matrix}-1\\0\end{matrix}\right],\qquad t=\frac{\pi}{6}
\end{aligned}
$$

二次量を具体的に求め、微分した値と流束または外力の寄与を照合する。質量行列が単位行列でない場合は係数の普通の二乗和がエネルギーになるとは限らない。

$$
\begin{aligned}
e&=20\\
e_t&=\frac{d}{d t} 20=0\\
e_t&=0
\end{aligned}
$$

試験方向を第一成分だけに取ると、時刻零の項は変位の $-2\eta^{\prime}(0)$ のみである。第二成分だけを取ると速度の $6\eta(0)$ のみである。合成した試験関数では両方が残り、それぞれを独立に指定できる。
:::

## エネルギー評価と一意性

:::theorem[一意なエネルギー解と評価]{#thm-pde-hyperbolic-energy}
:ref[def-pde-hyperbolic-weak] のデータに対して弱解は一つ存在する。$E(t)=\frac12(\norm{u'}_H^2+a(u,u))$ と置くと

$$
E(t)=E(0)+\int_0^t(f,u')_H\,\dd s,\qquad
\sqrt{2E(t)}\le\sqrt{2E(0)}+\int_0^t\norm f_H\,\dd s
\label{eq-pde-hyperbolic-energy}
$$

である。
:::

:::proof
:ref[thm-pde-elliptic-spectrum] の正規直交固有系を使い、係数 $u_j$ を $u_j''+\lambda_j u_j=f_j$ と初期値 $(g_j,h_j)$ で解く。有限和 $u_N$ の式に $u_N'$ を掛けて足せば $E_N'=(P_N f,u_N')_H$ である。$F_\delta=\sqrt{2E_N+\delta}$ とすると $F_\delta'\le\norm f_H$ であり、積分して $\delta\downarrow0$ で評価を得る。二つの部分和の差にも同じ評価を適用する。$g$ の $V$ 展開と $h$ の $H$ 展開は強収束し、$\int\norm{(P_M-P_N)f}_H\to0$ は優収束で従う。よって $u_N$ は $C(V)$、$u_N'$ は $C(H)$ で Cauchy である。各近似の弱式を得るには $P_Nv$ を用いればよく、$P_N$ は $V$ でも強収束するので極限から弱式と初期条件を得る。エネルギー恒等式は一様収束と $f\in L^1(H)$ で極限に移せる。一意性は任意の弱解の係数 $(u,e_j)_H$ が同じ ODE と初期条件を分布の意味で満たすことによる。各係数は一意で基底が完全である。
:::

時間依存の対称行列 $A(t,x)$ について、滑らかな解には $E'=(f,u')+\frac12\int (\partial_t A)Du\cdot Du$ が成立する。$A\ge\lambda I$、$\partial_tA$ が有界なら追加項は $(\norm{\partial_t A}_\infty/\lambda)E$ 以下である。従って積分因子で指数定数を許す評価が得られる。時間に依存する固有系は変化するため、弱解の構成には固定した試験基底を使い、一意性には時間で積分した試験関数を使う。

:::proposition[時間依存の対称形式]{#prop-pde-time-dependent-wave}
$V=H^1_0(\Omega)$、$H=L^2(\Omega)$ とし、$a_t$ は一様有界、対称、強圧的な形式とする。さらに $t\mapsto a_t$ は形式ノルムについて Lipschitz とし、各 $v,w$ に対する a.e. の微分を $\dot a_t(v,w)$ とする。$g\in V$、$h\in H$、$f\in L^1(0,T;H)$ に対し、$u\in L^\infty(0,T;V)$、$u'\in L^\infty(0,T;H)$ で :eqref[eq-pde-hyperbolic-weak] の $a$ を $a_t$ に替えた弱式と初期値を満たす解が存在し、このクラスで一意である。また $uin C([0,T];V)$、$u'in C([0,T];H)$ であり、二つの初期値はそれぞれの強い跡として成り立つ。
:::

:::proof
固定した Dirichlet Laplace 基底の有限部分に射影し、連立 ODE を解く。有限次元の解では $E_N'=(f,u_N')+\dot a_t(u_N,u_N)/2$ である。強圧性から $|\dot a_t(u_N,u_N)|\le C E_N$ なので、$\sqrt{2E_N+\delta}$ の微分を積分因子で抑えると $u_N$ は $L^\infty(V)$、$u_N'$ は $L^\infty(H)$ で一様有界である。射影の $V$ 有界性により $u_N''$ は $L^1(V^*)$ で一様有界である。前二者から弱星収束を取り、有限個の基底と滑らかな時間関数で試験すれば弱式を得る。時間微分の定義と初期項を含めた部分積分を取れば初期変位 $g$ と初速度 $h$ が極限で保たれる。極限の $u$ は $u'\in L^\infty(H)$ により $H$ で強連続であり、$u''=f-A(t)u\in L^1(V^*)$ により $u'$ は $V^*$ で連続である。後者と $L^\infty(H)$ 上界を合わせると $u'$ は弱 $H$ で連続であるため、初期項を比較できる。

一意性には差 $w$ と固定した時刻 $s$ に対し $v(t)=\int_t^s w(r)\,\dd r$（$0\le t\le s$）、$v=0$（$t>s$）を使う。$v\in H^1(0,T;V)$ で終端側で消えるため、平滑化と台の切断で弱式に許される試験関数へ近似できる。$v'=-w$ が成り立つ区間 $s<T$ で計算し、必要なら $s\uparrow T$ とする。初期値はともに零であり、$\partial_t a_t(v,v)=\dot a_t(v,v)+2a_t(v',v)$ なので

$$
\norm{w(s)}_H^2+a_0(v(0),v(0))=-\int_0^s\dot a_t(v(t),v(t))\,\dd t .
$$

形式の積の法則は、有限基底への射影で成り立つ法則を $H^1_tV$ の収束で移せる。$z(t)=\int_0^t w(r)\,\dd r$ と置くと $v(t)=z(s)-z(t)$ である。強圧定数を $\alpha$、$\norm{\dot a_t}\le L$ とすれば

$$
\alpha\norm{z(s)}_V^2\le 2Ls\norm{z(s)}_V^2+2L\int_0^s\norm{z(t)}_V^2\,\dd t .
$$

$2Ls\le\alpha/2$ の区間では第一項を左へ移せる。Grönwall の不等式により $z=0$ なので $w=0$ である。この区間の終点では連続な跡が両方とも零なので、同じ長さの区間を順に適用すれば全区間で一意である。

強い跡を確認する。まず $u$ は $H$ で強連続かつ $V$ で一様有界なので、各時刻に $V$ 値を持つ弱連続代表が得られる。実際、時刻を近づけたとき $V$ の任意の弱極限は、$H$ での極限に一致する。Galerkin 列も各時刻で $u$ に $V$ 弱収束し、$u'$ に $H$ 弱収束する部分列を取れる。$u_N$ は $H$ で等連続であり、$u_N'$ の $V^*$ での等連続性は方程式から
$$
\norm{u_N'(t)-u_N'(s)}_{V^*}
\le C|t-s|+C\int_s^t\norm f_H
$$
で得られる。稠密な可算個の試験ベクトルとの組を Ascoli の定理で一様収束させ、時刻積分で得た弱極限と照合すればよい。

任意の時刻 $s$ から $u(s),u'(s)$ をデータとして同じ近似を作る。有限次元のエネルギー評価を各時刻の弱下半連続性で移し、一意性で元の解に戻すと
$$
\sqrt{2E(t)}\le e^{C(t-s)}\left(\sqrt{2E(s)}+\int_s^t\norm f_H\right)
\qquad(t\ge s)
$$
となる。初期時刻の項は射影が $V$ と $H$ で強収束するので極限で等しい。時刻を反転した方程式も同じ仮定を満たし、速度の符号を反転して同じ議論を適用できる。よって $t\le s$ でも、指数の差を $s-t$、積分を $\int_t^s$ とした不等式が成り立つ。これにより $\limsup_{t\to s}E(t)\le E(s)$ を得る。一方、弱連続性と $a_t$ の形式ノルム連続性から $E(s)\le\liminf_{t\to s}E(t)$ である。従ってノルムが収束する。$a_s$ を $V$ の同値な内積とみなすと、弱収束とノルムの収束により $(u(t),u'(t))\to(u(s),u'(s))$ が $V\times H$ で成り立つ。
:::

:::example[仕事率と符号]
$(0,\pi)$ で $u=\sin t\sin x$、$f=0$ のとき $E=\pi/4$ である。$u=t\sin x$ と保つには $f=t\sin x$ が必要で、$E=\pi(1+t^2)/4$、$E'=\pi t/2$ は $(f,u')$ と一致する。
:::

:::example[初期速度だけを近似する場合]
初期変位を固定し、初速度の誤差が $H$ ノルムで $1/10$ とする。外力が同じなら差の全エネルギーは $1/200$ で保たれる。従って任意の時刻で速度の誤差は $1/10$ 以下、変位の勾配の誤差は $\lambda^{-1/2}/10$ 以下である。時間とともに運動項から変形項へ移っても、この二つの上界は変わらない。
:::

### エネルギーが強圧でない場合のずらし

$a(u,u)+\gamma\norm u_H^2\ge\alpha\norm u_V^2$ とし、$a$ は対称とする。$E_\gamma=(\norm{u'}_H^2+a(u,u)+\gamma\norm u_H^2)/2$ は非負で、滑らかな解について
$$
E_\gamma'=(f,u')_H+\gamma(u,u')_H\le\sqrt{2E_\gamma}\norm f_H+C E_\gamma .
$$
最後の項は $V\hookrightarrow H$ と強圧性による。$\sqrt{E_\gamma+\eps}$ を使って割り、積分因子で有限時間の評価を得る。従って負の固有値が有限個あって波が成長しても、有限時間の連続依存性は失われない。例えば $A=-\partial_{xx}-2$ の第一正弦モードは $u_1''-u_1=0$ であり、初期値 $(1,0)$ に対し $\cosh t$ になる。元の符号付きエネルギー $((u_1')^2-u_1^2)/2=-1/2$ は保存されるが、その値だけでは大きさを上から押さえられない。

:::example[正と負の固有値を同じ初期データに含める]
$(0,\pi)$ で $u_{tt}-u_{xx}-2u=0$、初期変位 $\sin x+\sin2x$、初期速度零とする。符号付きエネルギーと成長モードを比べる。

まず初期値を固有関数ごとに分け、対応する時間方程式を解く。外力がある場合も、同じ固有関数に属する外力だけを時間方程式に入れる。以下では得られた表示の微分を取って、式の各項の関係を確かめる。

$$
\begin{aligned}
u&=\sin{\left(x \right)} \cosh{\left(t \right)} + \sin{\left(2 x \right)} \cos{\left(\sqrt{2} t \right)}\\
u_{tt}&=\frac{\partial^{2}}{\partial t^{2}} \left(\sin{\left(x \right)} \cosh{\left(t \right)} + \sin{\left(2 x \right)} \cos{\left(\sqrt{2} t \right)}\right)=- \left(4 \cos{\left(x \right)} \cos{\left(\sqrt{2} t \right)} - \cosh{\left(t \right)}\right) \sin{\left(x \right)}\\
u_{xx}&=\frac{\partial^{2}}{\partial x^{2}} \left(\sin{\left(x \right)} \cosh{\left(t \right)} + \sin{\left(2 x \right)} \cos{\left(\sqrt{2} t \right)}\right)=- \left(8 \cos{\left(x \right)} \cos{\left(\sqrt{2} t \right)} + \cosh{\left(t \right)}\right) \sin{\left(x \right)}\\
f&=0
\end{aligned}
$$

境界条件を各モードで確かめてから、時刻零のデータに戻る。端点値がすべて零でも時間の係数は任意に選べない。時刻を零に戻して初めて各係数の初期条件が確認できる。

$$
\begin{aligned}
u(0,x)&=\sin{\left(x \right)} + \sin{\left(2 x \right)},\qquad t=0\\
u_t(0,x)&=0,\qquad t=0\\
u(t,0)&=0,\qquad x=0\\
u(t,L)&=0,\qquad x=\pi
\end{aligned}
$$

エネルギーを空間積分で測る。異なるモードの積の積分は零であり、一つずつ積分して足し合わせられる。その時間微分を取り、解の式を使って積分する方法と比較する。

$$
\begin{aligned}
e(t,x)&=\frac{\left(\sin{\left(x \right)} \sinh{\left(t \right)} - \sqrt{2} \sin{\left(2 x \right)} \sin{\left(\sqrt{2} t \right)}\right)^{2}}{2} - \left(\sin{\left(x \right)} \cosh{\left(t \right)} + \sin{\left(2 x \right)} \cos{\left(\sqrt{2} t \right)}\right)^{2} + \frac{\left(\cos{\left(x \right)} \cosh{\left(t \right)} + 2 \cos{\left(2 x \right)} \cos{\left(\sqrt{2} t \right)}\right)^{2}}{2}\\
\mathcal E(t)&=\frac{\pi}{4}\\
\mathcal E^{\prime}&=\frac{d}{d t} \frac{\pi}{4}=0\\
\text{入力と散逸の和}&=\int\limits_{0}^{\pi} 0\, dx=0
\end{aligned}
$$

最後に一つの観測点で変位または温度と勾配を計算する。点で値が相殺されることと、区間全体のエネルギーが零になることは区別する。

$$
\begin{aligned}
u&=\cosh{\left(\frac{\sqrt{2} \pi}{2} \right)},\qquad t=\frac{\sqrt{2} \pi}{2},x=\frac{\pi}{2}\\
u_x&=2,\qquad t=\frac{\sqrt{2} \pi}{2},x=\frac{\pi}{2}\\
\mathcal E&=\frac{\pi}{4},\qquad t=\frac{\sqrt{2} \pi}{2}
\end{aligned}
$$

第1モードは負固有値 $-1$ に対応し、$\cosh t$ で成長する。第2モードは固有値 $2$ により振動する。符号付きエネルギーは保存されても第1モードの負のポテンシャル項と正の速度項が打ち消し合い、振幅の上界にはならない。有限時間の評価では正のずらしを入れる。
:::

### 強圧性を用いた Galerkin 極限での初期値の確認

時間で Lipschitz な対称形式 $a_t$ に対し、有限次元近似が $L^\infty(0,T;V)$ および一階微分について $L^\infty(0,T;H)$ で一様有界とする。近似空間への射影は $V$ で一様有界に取り、近似式から $u_N''$ は $L^1(0,T;V^*)$ で一様有界である。さらに $a_t(u_N,\cdot)$ のノルムは $L^\infty$ で一様有界なので、時間微分を積分した差は $C|t-s|+\int_s^t\|f\|_{V^*}$ 以下となり等連続である。ここでは外力が固定で $L^1$ に入ることを使う。

稠密な可算個の試験ベクトルに対して初期値を含む時間積分を取り、固定した試験ベクトルについて極限へ移す。$a_t$ は有界なので、時間試験関数を掛けた $a_t(\cdot,v)$ は $L^1(0,T;V^*)$ に属し、$L^\infty(0,T;V)$ での弱星極限を取れる。時刻零の境界項は $g_N,h_N$ の強収束で元のデータの項になる。ここで得られる弱い跡から強い跡へ進むには、:ref[prop-pde-time-dependent-wave] のように各時刻のエネルギーを両側から比較する。時間一定の場合は固有係数の差の評価により同じ結論がより短く従う。

:::example[速度比例の減衰が正の仕事を取り除く]
一つの波動モードで $a^{\prime\prime}+4a^{\prime}+4a=0$、$a(0)=1$、$a^{\prime}(0)=-1$ とする。重根の時間解と散逸を計算する。

各行の質量係数と剛性係数を区別して、候補を連立方程式に代入する。$B$ の非対角成分は別の成分を結び付ける。従って個々の微分の値を求めるだけではなく、全成分を掛け合わせた結果を確認する。

$$
\begin{aligned}
M&=\left[\begin{matrix}1\end{matrix}\right]\\
B&=\left[\begin{matrix}4\end{matrix}\right]\\
C&=\left[\begin{matrix}4\end{matrix}\right]\\
U&=\left[\begin{matrix}\left(t + 1\right) e^{- 2 t}\end{matrix}\right]\\
U_{tt}&=\left[\begin{matrix}4 t e^{- 2 t}\end{matrix}\right]\\
F&=\left[\begin{matrix}0\end{matrix}\right]
\end{aligned}
$$

初期時刻へ戻して、先に指定したデータと比較する。次に正の時刻における成分の値を読む。成分が混ざっていても、初期データは方程式を積分する前の元の座標で指定される。

$$
\begin{aligned}
U(0)&=\left[\begin{matrix}1\end{matrix}\right],\qquad t=0\\
U_t(0)&=\left[\begin{matrix}-1\end{matrix}\right],\qquad t=0\\
U&=\left[\begin{matrix}\frac{2}{e^{2}}\end{matrix}\right],\qquad t=1\\
U_t&=\left[\begin{matrix}- \frac{3}{e^{2}}\end{matrix}\right],\qquad t=1
\end{aligned}
$$

二次量を具体的に求め、微分した値と流束または外力の寄与を照合する。質量行列が単位行列でない場合は係数の普通の二乗和がエネルギーになるとは限らない。

$$
\begin{aligned}
e&=\frac{\left(8 t^{2} + 12 t + 5\right) e^{- 4 t}}{2}\\
e_t&=\frac{d}{d t} \frac{\left(8 t^{2} + 12 t + 5\right) e^{- 4 t}}{2}=- 4 \left(2 t + 1\right)^{2} e^{- 4 t}\\
e_t&=\left(- 2 \left(t + 1\right) e^{- 2 t} + e^{- 2 t}\right) \left(8 \left(t + 1\right) e^{- 2 t} - 4 e^{- 2 t}\right)=- 4 \left(2 t + 1\right)^{2} e^{- 4 t}
\end{aligned}
$$

特性方程式は $(\lambda+2)^2=0$ であり、時間の独立解は $e^{-2t},te^{-2t}$ である。初期条件に合わせると両係数は $1$。二次エネルギーの微分は $-4(a^{\prime})^2\le0$ で、正のばね定数を保ったまま速度比例項だけが仕事を失わせる。
:::

:::example[一定の外力を受ける減衰モード]
速度比例の減衰を持つモード $a^{\prime\prime}+2a^{\prime}+2a=2$、初期値は変位も速度も零とする。

時間変換の変数を $s$ とし、$\mathcal L[f](s)=\int_0^\infty e^{-st}f(t)\dd t$ とする。大きい $s$ の範囲で積分が収束する。時間微分を部分積分すると、$t=0$ の値が分子の項になる。まず外力の変換、微分から来る分母、初期データから来る分子を分ける。

$$
\begin{aligned}
f(t)&=2\\
\widetilde f(s)&=\frac{2}{s}\\
\text{微分の多項式}&=s^{2} + 2 s + 2\\
\text{初期項}&=0\\
\widetilde a(s)&=\frac{2}{s^{3} + 2 s^{2} + 2 s}
\end{aligned}
$$

次に分母を因数分解する。一つの極が重なるときは単純な指数の係数に加え、時間のべきが掛かる項を求める。例えば $(s+\lambda)^{-2}$ は $te^{-\lambda t}$、分母の二次式 $(s+\lambda)^2+\omega^2$ は減衰する正弦または余弦に対応する。

$$
\begin{aligned}
\widetilde a(s)&=- \frac{s + 2}{s^{2} + 2 s + 2} + \frac{1}{s}\\
a(t)&=\left(e^{t} - \sqrt{2} \sin{\left(t + \frac{\pi}{4} \right)}\right) e^{- t}\\
a^{\prime}&=\frac{d}{d t} \left(e^{t} - \sqrt{2} \sin{\left(t + \frac{\pi}{4} \right)}\right) e^{- t}=2 e^{- t} \sin{\left(t \right)}\\
a^{\prime\prime}&=\frac{d^{2}}{d t^{2}} \left(e^{t} - \sqrt{2} \sin{\left(t + \frac{\pi}{4} \right)}\right) e^{- t}=2 \sqrt{2} e^{- t} \cos{\left(t + \frac{\pi}{4} \right)}
\end{aligned}
$$

最後に時刻を戻して初期条件を調べ、次に得た時間関数を元の左辺に代入する。変換領域の式を満たしても、逆変換で定数や初期項を落とすとこの二つの確認のどちらかに失敗する。

$$
\begin{aligned}
a(0)&=0,\qquad t=0\\
a^{\prime}(0)&=0,\qquad t=0\\
\text{時間方程式の左辺}&=2 \left(e^{t} - \sqrt{2} \cos{\left(t + \frac{\pi}{4} \right)}\right) e^{- t} + 2 \sqrt{2} e^{- t} \cos{\left(t + \frac{\pi}{4} \right)}=2\\
a(t_*)&=\frac{-1 + e^{\frac{\pi}{2}}}{e^{\frac{\pi}{2}}},\qquad t=\frac{\pi}{2}\\
a^{\prime}(t_*)&=- \frac{-1 + e^{\frac{\pi}{2}}}{e^{\frac{\pi}{2}}} + \frac{1 + e^{\frac{\pi}{2}}}{e^{\frac{\pi}{2}}}=\frac{2}{e^{\frac{\pi}{2}}},\qquad t=\frac{\pi}{2}
\end{aligned}
$$

一定の外力は静止平衡 $a=1$ を作る。差 $z=a-1$ のエネルギー $((z^{\prime})^2+2z^2)/2$ を微分すれば $-2(z^{\prime})^2$ となる。従って入力がずっとあることと、平衡から測ったエネルギーが減ることは同時に成り立つ。
:::

## 擾乱の有限伝播速度

:::proposition[変数係数の円錐評価]{#prop-pde-variable-wave-cone}
全空間で $u''-\divg(A Du)=0$ とし、$A=A(x)$ は滑らかな実対称行列で $0<\lambda I\le A\le c_*^2 I$ とする。滑らかな解について、$B_{c_*(T-t)}(x_0)$ 内のエネルギーは非増加である。
:::

:::proof
$e=(u_t^2+Du\cdot A Du)/2$ について $\partial_t e=\divg(u_t A Du)$ が成り立つ。移動球の境界で

$$
|u_t\,\nu^\top A Du|
\le |u_t|(\nu^\top A\nu)^{1/2}(Du\cdot A Du)^{1/2}
\le c_* e .
\label{eq-pde-variable-wave-flux}
$$

最初に $A$ 内積の Cauchy–Schwarz、次に $2ab\le a^2+b^2$ を使った。球の境界は速度 $-c_*$ で内向きに移動するため、エネルギー微分の被積分量 $u_t\nu^\top A Du-c_* e$ は非正になる。
:::

初期の球上で二解が一致し円錐内で外力も一致すれば、差のエネルギーは円錐内で零であり、差の変位と速度はそこで消える。このように領域境界周辺を使わずに結論が出る。係数が $A(t,x)$ で $\partial_tA$ 有界でも、移動球上で右辺に $C E(t)$ を加えるだけなので、初期エネルギーが $0$ なら円錐内は引き続き $0$ である。

:::example[異方性での方向別上界]
$A=\diag(9,1)$ では $c_*=3$ を選べる。法線 $e_2$ の平面を動かすときは $\nu^\top A\nu=1$ なので流束は $e$ 以下であり、同じ評価でその方向の速さは $1$ で十分である。定数行列なら $x=A^{1/2}y$ で速度 $1$ の等方的波となるので、依存領域の断面は半軸 $3t,t$ の楕円である。
:::

:::example[外力が到着できる時刻]
$A=4I$ とし、初期値はゼロ、外力は時刻 $s\ge1$ かつ空間 $B_1(0)$ に台を持つとする。点 $(4,0)$ までの最短距離は $3$ なので、$t<1+3/2$ ではそこに後ろ向き円錐が外力の台と交わらない。従って、その時刻の値はゼロである。
:::

### 異方的な円錐を座標変換で得る

定数対称正定値 $A$ に対し $u_{tt}-\divg(A Du)=0$ とする。$x=A^{1/2}y$、$v(t,y)=u(t,A^{1/2}y)$ と置けば $v_{tt}-\Delta_y v=0$ である。従って依存領域の初期断面は
$$
(x-x_0)^\top A^{-1}(x-x_0)\le t^2
$$
で与えられる。この楕円体の法線方向 $\nu$ への支持関数は $t\sqrt{\nu^\top A\nu}$ である。実際、$x-x_0=A^{1/2}y$、$|y|\le t$ として $\nu\cdot A^{1/2}y$ の Cauchy–Schwarz の等号を取ればよい。例えば $A=\diag(9,1)$ なら $t=2$ で半軸は $6,2$ である。点 $(0,4)$ は半径 $6$ の円に入るが楕円体の外なので、一般の最大速度 $3$ の球評価より強い結論が得られる。

:::example[移動境界に必要な速度]
対称一階系で $A^0=I$ とし、ある法線方向の速度行列を調べる。

法線記号の行列を $B$ とする。対角成分は座標方向での値を与える。固有方程式から直交する方向を取り出すため、特性多項式 $\det(\lambda I-B)$ を計算する。

$$
\begin{aligned}
B&=\left[\begin{matrix}\frac{5053}{4225} & - \frac{13296}{4225} & \frac{144}{169}\\- \frac{13296}{4225} & - \frac{2703}{4225} & - \frac{108}{169}\\\frac{144}{169} & - \frac{108}{169} & \frac{244}{169}\end{matrix}\right]\\
\det(\lambda I-B)&=\lambda^{3} - 2 \lambda^{2} - 11 \lambda + 12=\left(\lambda - 4\right) \left(\lambda - 1\right) \left(\lambda + 3\right)\\
d_1&=-3,\quad v_1=\left[\begin{matrix}\frac{3}{5}\\\frac{4}{5}\\0\end{matrix}\right]\\
d_2&=1,\quad v_2=\left[\begin{matrix}- \frac{4}{13}\\\frac{3}{13}\\\frac{12}{13}\end{matrix}\right]\\
d_3&=4,\quad v_3=\left[\begin{matrix}\frac{48}{65}\\- \frac{36}{65}\\\frac{5}{13}\end{matrix}\right]
\end{aligned}
$$

$V$ の列に固有ベクトルを並べる。$V^\top V=I$ と $V^\top BV=D$ を使えば、元の座標 $X$ と固有座標 $Y=V^\top X$ の双方で状態ベクトルを記述できる。二次形式は固有座標にすると交差項を持たない。

$$
\begin{aligned}
V&=\left[\begin{matrix}\frac{3}{5} & - \frac{4}{13} & \frac{48}{65}\\\frac{4}{5} & \frac{3}{13} & - \frac{36}{65}\\0 & \frac{12}{13} & \frac{5}{13}\end{matrix}\right]\\
V^\top BV&=\left[\begin{matrix}-3 & 0 & 0\\0 & 1 & 0\\0 & 0 & 4\end{matrix}\right]\\
X^\top BX&=\frac{5053 x_{1}^{2} - 26592 x_{1} x_{2} + 7200 x_{1} x_{3} - 2703 x_{2}^{2} - 5400 x_{2} x_{3} + 6100 x_{3}^{2}}{4225}\\
\sum_j d_jY_j^2&=- 3 Y_{1}^{2} + Y_{2}^{2} + 4 Y_{3}^{2}
\end{aligned}
$$

一つの状態ベクトルで固有係数を計算し、流束の二次形式を求める。直交変換なので二乗ノルムは変わらない。二次形式を二乗ノルムで割ると、固有値を非負の重みで平均した値になる。

$$
\begin{aligned}
X&=\left[\begin{matrix}1\\2\\1\end{matrix}\right]\\
Y&=\left[\begin{matrix}\frac{11}{5}\\\frac{14}{13}\\\frac{1}{65}\end{matrix}\right]\\
|X|^2&=6\\
X^\top BX&=- \frac{56443}{4225}\\
\frac{X^\top BX}{|X|^2}&=- \frac{56443}{25350}\\
\min_j d_j&=-3\\
\max_j d_j&=4
\end{aligned}
$$

$q\cdot\nu-ce$ はエネルギーを固定すると、速度 $4$ の固有ベクトルで最大になる。$c=4$ なら非正、$c=3$ ならその固有ベクトル上で正になる。反対方向では行列が $-B$ となって最大速度は $3$ になる。
:::

:::example[異方的な波の同じ次数に異なる係数が出る]
二次元で $P=9\partial_{xx}+\partial_{yy}$、$g=x^4+y^4$、$h=xy$ とする。基準となる楕円体の主軸は座標軸に平行である。

時間二階の方程式なので偶数階と奇数階を分ける。$g$ からは $g,Pg,P^2g,\ldots$、$h$ からは $h,Ph,P^2h,\ldots$ を得る。空間の多項式に定数係数の二階作用素を掛けると次数が少なくとも二つ下がる。従ってこの場合には両方の列が有限回で零になり、級数の収束を別に仮定する必要はない。

$$
\begin{aligned}
P^{0}g&=x^{4} + y^{4}\\
P^{1}g&=12 \left(9 x^{2} + y^{2}\right)\\
P^{2}g&=1968\\
P^{0}h&=x y\\
P^{3}g&=0
\end{aligned}
$$

$P^j g$ を $t^{2j}/(2j)!$、$P^j h$ を $t^{2j+1}/(2j+1)!$ と掛けて足す。時間で二回微分すると各項の添字が一つ下がり、空間作用素を一回掛けた列と一致する。時刻零の変位と速度は別々に確認する。

$$
\begin{aligned}
u(t,x)&=82 t^{4} + 6 t^{2} \left(9 x^{2} + y^{2}\right) + t x y + x^{4} + y^{4}\\
u_{tt}&=\frac{\partial^{2}}{\partial t^{2}} \left(82 t^{4} + 6 t^{2} \left(9 x^{2} + y^{2}\right) + t x y + x^{4} + y^{4}\right)=12 \left(82 t^{2} + 9 x^{2} + y^{2}\right)\\
Pu&=12 \left(82 t^{2} + 9 x^{2} + y^{2}\right)\\
u(0,x)&=x^{4} + y^{4},\qquad t=0\\
u_t(0,x)&=x y,\qquad t=0
\end{aligned}
$$

空間の導関数を測ると、時間の多項式の各係数が局所でどう効くかが分かる。初期関数との差は単に平行移動で置き換えず、以下では導関数と時刻を代入して計算する。

$$
\begin{aligned}
\partial_{x} u&=\frac{\partial}{\partial x} \left(82 t^{4} + 6 t^{2} \left(9 x^{2} + y^{2}\right) + t x y + x^{4} + y^{4}\right)=108 t^{2} x + t y + 4 x^{3}\\
\partial_{y} u&=\frac{\partial}{\partial y} \left(82 t^{4} + 6 t^{2} \left(9 x^{2} + y^{2}\right) + t x y + x^{4} + y^{4}\right)=12 t^{2} y + t x + 4 y^{3}\\
g&=2,\qquad x=1,y=1\\
u&=145,\qquad t=1,x=1,y=1\\
u_t&=449,\qquad t=1,x=1,y=1\\
h&=1,\qquad x=1,y=1
\end{aligned}
$$

変数を $X=x/3$、$Y=y$ と置くと標準波動方程式である。しかし初期多項式は $81X^4+Y^4$ に変わり、観測点は $(1/3,1)$ へ移る。この変換で半径 $t$ の依存領域を元に戻すと、$x$ 方向では半径 $3t$、$y$ 方向では $t$ である。
:::

### 時間に依存する一様な波速の到達範囲

$c=c(t)>0$ を $C^1$ とし、$u_{tt}-c(t)^2\Delta u=f$ を考える。観測点 $(T,x_0)$ に対し半径 $R(t)=\int_t^T c(s)\dd s$ の球を用いる。$e=(u_t^2+c(t)^2|Du|^2)/2$ とすれば
$$
\partial_t e-\divg(c(t)^2 u_tDu)=fu_t+c(t)c'(t)|Du|^2 .
$$
動く球での境界項は $c(t)^2u_t\partial_\nu u+R'(t)e=c(t)^2u_t\partial_\nu u-c(t)e\le0$ である。$c$ が正の下限を持つ時間区間では係数微分の項は $2|c'|e/c$ 以下なので、初期エネルギーと外力が円錐にない場合は積分因子で解の差が零になる。上限の一定速度を使った球より、時間積分で半径を決める方が狭い。

:::example[拡大する速度の履歴]
$c(t)=1+t$、$T=2$ なら依存領域の半径は $\int_0^2(1+s)\dd s=4$ である。最大速度 $3$ を全区間に使う半径 $6$ より小さい。時刻 $1$ からの入力が到達し得る初期断面は $R(1)=\int_1^2(1+s)\dd s=5/2$ である。源の時刻が異なると使う半径も変わるため、時刻 $0$ での最大の距離だけでは入力の到達時刻を決められない。
:::

## 対称双曲型系

:::definition[対称双曲型系]{#def-pde-symmetric-hyperbolic}
実ベクトル場 $U$ に対する一階系 $A^0 U_t+\sum_{j=1}^n A^j\partial_jU+BU=F$ において、$A^0$ が対称で一様に正定値、$A^j$ が対称であるとき :term[対称双曲型]{en="symmetric hyperbolic" yomi="たいしょうそうきょくがた"} という。
:::

:::proposition[一階系のエネルギー恒等式]{#prop-pde-symmetric-flux}
滑らかな解について $e=U^\top A^0U/2$、$q_j=U^\top A^jU/2$ とおくと

$$
\partial_t e+\sum_j\partial_j q_j
=U^\top F+U^\top\left(\tfrac12\partial_tA^0+\tfrac12\sum_j\partial_jA^j-B\right)U .
\label{eq-pde-symmetric-flux}
$$

:::

:::proof
方程式に左から $U^\top$ を掛ける。$A^0$ の対称性で $U^\top A^0 U_t$ が $\partial_te-U^\top(\partial_t A^0)U/2$ に等しい。同様に各空間項を流束の微分と係数微分に分ければ式を得る。
:::

$A^0$ の最小固有値が一様に正であり、係数微分と $B$ が有界なら、右辺の下位項は $Ce$ で抑えられる。解が $L^\infty(0,T;L^2)$、外力が $L^1(0,T;L^2)$ に属し、主係数も有界とする。球 $B_R$ で $1$、$B_{2R}$ の外で $0$ の滑らかな切断を掛けて積分すると、流束の誤差は $CR^{-1}\norm U_2^2$ 以下である。$R\to\infty$ でこの誤差を消し、全空間で $\frac{d}{dt}\sqrt{\int e}\le C\sqrt{\int e}+C\norm F_2$ を得る。法線 $\nu$ に沿う流束の絶対値は $c_* e$ 以下である。ここで

$$
c_*=\sup_{t,x,|\nu|=1}\left\|(A^0)^{-1/2}
\left(\sum_j\nu_j A^j\right)(A^0)^{-1/2}\right\|.
$$

従って円錐上でも同じエネルギー評価が成立する。$A^j$ の対称性により主部の積分で空間導関数の損失が起きない。

### 全空間での近似の組み立て

$A^0=I$ とし、$A^j$ が有界で空間 Lipschitz、$B$ が有界な場合を考える。係数は時刻について可測で、これらの上界は一様とする。空間の偶な平滑化作用素を $J_\eps$ として、平滑化方程式

$$
U_\eps'=-\sum_jJ_\eps(A^j\partial_j J_\eps U_\eps)-BU_\eps+F
$$

を解く。これは固定 $\eps$ では $L^2$ 上の有界 ODE であり、短区間での積分方程式の縮小を繰り返せば可解である。$J_\eps$ が自己共役なのでエネルギー積分において主部の両方に $J_\eps$ が作用し、部分積分で残るのは $(\partial_j A^j)$ と $J_\eps U_\eps$ の二次式である。従って得られる $L^\infty_t L^2_x$ 評価は $\eps$ に依存せず、$U_\eps'$ も $L^1_tH^{-1}_x$ で有界である。弱収束を取り、滑らかな試験関数に平滑化と微分を移せば極限が分布解で初期跡を $H^{-1}$ に持つことが従う。

一意性には分布解そのものを $J_\eps$ で平滑化する。このとき主部との交換子 $[J_\eps,A^j]\partial_j U$ が残る。そのノルムは $C\operatorname{Lip}(A^j)\norm U_2$ 以下である。実際、核に $A^j(y)-A^j(x)$ が掛かり、$y$ で部分積分すると核の微分には $|x-y|$ の因子が残る。従って $\int |z||D J_\eps(z)|\,\dd z$ と $\int J_\eps$ が一様に有界なら Young の畳み込み評価が使える。滑らかな $U$ では交換子が強く $0$ に収束し、一般の $U$ でも $L^2$ の密度と一様評価から強収束する。$B$ の交換子も強収束し、$F\in L^1_tL^2_x$ なら強収束を時間積分に移せる。従って平滑化解のエネルギー恒等式を極限に移せる。弱連続な $L^2$ 代表は、$H^{-1}$ 連続性と一様 $L^2$ 上界から得られる。各時刻の弱連続代表に $J_\eps$ を作用させれば $L^2$ で強収束し、極限の恒等式は $L^2$ ノルムの連続性を与える。弱連続性と合わせ強連続性が従う。一意性と局所エネルギー評価も同様に従う。局所評価では球の特性関数の代わりに滑らかな移動切断を用い、交換子の $L^1_tL^2_x$ 収束でその誤差が消える。

:::example[波を速度と勾配の系にする]
$c$ が一定なら $U=(u_t,cDu)$ と置き、$v=u_t$、$w=cDu$ と書くと $v_t-c\divg w=f$、$w_t-cDv=0$ である。$A^0=I$、$A^j$ の非零ブロックは $-ce_j^\top$ と $-ce_j$ なので対称である。法線記号の固有値は $\pm c$ と、接方向があれば $0$ であり、円錐上界は $c$ となる。
:::

### 時間成分の行列を単位行列にする

$A^0$ を含め主係数が時空間で有界 Lipschitz、$A^0\ge\lambda I$ とし、$B$ は有界とする。$C=(A^0)^{1/2}$、$W=CU$ とおいて元の式に左から $C^{-1}$ を掛ける。すると

$$
W_t+\sum_j (C^{-1}A^j C^{-1})\partial_jW+\widetilde B W=C^{-1}F,
\qquad
\widetilde B=C^{-1}BC^{-1}+C\partial_t C^{-1}
+\sum_j C^{-1}A^j\partial_j C^{-1}
$$

を得る。新しい空間係数は対称であり、$\widetilde B$ は有界である。平方根の Lipschitz 性は、二つの正行列の平方根 $C_1,C_2$ の差 $D$ が $C_1D+DC_2=A^0_1-A^0_2$ を満たすことから分かる。実際 $\int_0^\infty e^{-sC_1}(A^0_1-A^0_2)e^{-sC_2}\,\dd s=D$ であり、$\norm D\le (2\sqrt\lambda)^{-1}\norm{A^0_1-A^0_2}$ である。逆行列の差にも同様に有界性が従う。従って前の平滑化の構成が適用でき、初期値 $U(0)\in L^2$ と $F\in L^1_tL^2_x$ に対して $C([0,T];L^2)$ の分布解が一意に得られる。元の密度 $U^\top A^0U/2$ は $|W|^2/2$ なので、局所エネルギー評価もこの変換を介して適用できる。

:::example[場所による速度変化と保存型の係数]
一次元で $U_t+b(x)U_x=0$、$b(x)=\sin x$ とする。$A^0=1$、$A^1=b$ なので $e=U^2/2$、$q=bU^2/2$ に対して $e_t+q_x=(\cos x)e$ である。従って全空間で $\int e(t)\le e^t\int e(0)$ である。速度の絶対値は $1$ 以下なので局所円錐の上界には $c_*=1$ を使える。
:::

### エネルギー弱解でも成立する円錐評価

全空間で $u_{tt}-\divg(A Du)=f$ を考える。$A=A(t,x)$ は実対称で一様楕円、時空間について有界 Lipschitz とする。$R=A^{1/2}$、$v=u_t$、$w=RDu$ とおけば

$$
v_t-\divg(Rw)=f,\qquad
w_t-RDv=(\partial_t R)R^{-1}w
$$

である。発散を展開すると、第 $j$ 空間微分の係数の非対角ブロックは $-e_j^\top R$ と $-Re_j$ である。いずれも互いに転置なので対称であり、係数の微分は有界な下位項に入る。従って全空間の平滑化の構成で、初期データ $(h,R(0)Dg)$ に対する $C_tL^2_x$ の解が得られる。ここで $g\in H^1(\R^n)$、$h\in L^2(\R^n)$、$f\in L^1_tL^2_x$ とする。第二式から分布として $\partial_t(R^{-1}w)=Dv$ なので、$u=g+\int_0^t v(s)\,\dd s$ が $Du=R^{-1}w$ を満たす。こうして $u\in C_tH^1_x$、$u_t\in C_tL^2_x$ の弱解になる。

密度は $(v^2+|w|^2)/2$、外向き流束は $-v\,\nu^\top Rw$ である。その絶対値は $\norm{R\nu}$ に密度を掛けたもの以下なので、$A\le c_*^2 I$ なら速度 $c_*$ の円錐で局所エネルギー評価が成り立つ。二解のデータが初期の断面内で一致し、$f$ の差が円錐内で零なら、差のエネルギーは初め零である。有界な下位項を積分因子で処理すると円錐の各断面で零になり、$v$ と $Du$ の差が消える。各内点から初期断面へ垂直にたどって $u=g+\int v$ を使えば、変位の差も零である。これにより有限伝播速度は滑らかなデータから弱解のデータへ移る。

### 音響方程式での対称化

$p_t+K\divg v=0$、$\rho v_t+Dp=0$ とし、$K,\rho>0$ は時間に依存しない空間関数である。第一式を $K$ で割ると $A^0=\diag(1/K,\rho I)$ になり、各 $A^j$ の非対角ブロックは $e_j^\top,e_j$ で対称である。従ってエネルギー $e=p^2/(2K)+\rho|v|^2/2$ と流束 $q=pv$ で $e_t+\divg q=0$ となる。位置依存の $K,\rho$ は $A^0$ に集めており、空間微分の係数行列は一定なので追加の生成項はない。法線方向の正規化行列には非対角の係数 $\sqrt{K/\rho}$ が出るため、局所波速はこの値である。例えば $K=8,\rho=2$ なら $c=2$ で、$p=4,v=1$ の一次元状態では $e=2$、$q=4=ce$ となる。

:::example[入流と流出に対応する成分]
半直線の端における法線記号が異なる符号を持つ場合を調べる。

法線記号の行列を $B$ とする。対角成分は座標方向での値を与える。固有方程式から直交する方向を取り出すため、特性多項式 $\det(\lambda I-B)$ を計算する。

$$
\begin{aligned}
B&=\left[\begin{matrix}\frac{6}{5} & - \frac{12}{5}\\- \frac{12}{5} & - \frac{1}{5}\end{matrix}\right]\\
\det(\lambda I-B)&=\lambda^{2} - \lambda - 6=\left(\lambda - 3\right) \left(\lambda + 2\right)\\
d_1&=-2,\quad v_1=\left[\begin{matrix}\frac{3}{5}\\\frac{4}{5}\end{matrix}\right]\\
d_2&=3,\quad v_2=\left[\begin{matrix}- \frac{4}{5}\\\frac{3}{5}\end{matrix}\right]
\end{aligned}
$$

$V$ の列に固有ベクトルを並べる。$V^\top V=I$ と $V^\top BV=D$ を使えば、元の座標 $X$ と固有座標 $Y=V^\top X$ の双方で状態ベクトルを記述できる。二次形式は固有座標にすると交差項を持たない。

$$
\begin{aligned}
V&=\left[\begin{matrix}\frac{3}{5} & - \frac{4}{5}\\\frac{4}{5} & \frac{3}{5}\end{matrix}\right]\\
V^\top BV&=\left[\begin{matrix}-2 & 0\\0 & 3\end{matrix}\right]\\
X^\top BX&=\frac{6 x_{1}^{2} - 24 x_{1} x_{2} - x_{2}^{2}}{5}\\
\sum_j d_jY_j^2&=- 2 Y_{1}^{2} + 3 Y_{2}^{2}
\end{aligned}
$$

一つの状態ベクトルで固有係数を計算し、流束の二次形式を求める。直交変換なので二乗ノルムは変わらない。二次形式を二乗ノルムで割ると、固有値を非負の重みで平均した値になる。

$$
\begin{aligned}
X&=\left[\begin{matrix}2\\1\end{matrix}\right]\\
Y&=\left[\begin{matrix}2\\-1\end{matrix}\right]\\
|X|^2&=5\\
X^\top BX&=-5\\
\frac{X^\top BX}{|X|^2}&=-1\\
\min_j d_j&=-2\\
\max_j d_j&=3
\end{aligned}
$$

外向きの流束を正とすれば、$d=-2$ のモードは入流、$d=3$ のモードは流出である。零入流条件 $v_1^\top U=0$ のもとでは流束が非負になる。すべての成分を零に指定すると流出モードまで拘束する。
:::

### 音響の界面における透過条件

一次元で $p_t/K+v_x=0$、$\rho v_t+p_x=0$ とし、正の係数 $K,\rho$ が $x=0$ の左右で定数を取る。$p,v$ が各半直線で滑らかで有界であり、界面に片側跡を持つとする。弱形式に界面近くの試験関数を使うと、$v_x,p_x$ のデルタが消えるため $[p]=[v]=0$ が必要である。逆にこの条件があれば、両側の古典式を合わせた弱形式からデルタは出ない。エネルギー流束 $pv$ も連続なので、界面自体で散逸しない。

各側の波速を $c_i=\sqrt{K_i/\rho_i}$、インピーダンスを $Z_i=\rho_i c_i$ とする。一定の背景を引き、右向きの固有モード $p=P(x-c_it)$ では $v=P/Z_i$、左向きの固有モードでは $v=-P/Z_i$ と取ると双方の式を満たす。左からの入射圧力 $A$、反射圧力 $B$、右への透過圧力 $C$ なら
$$
A+B=C,\qquad (A-B)/Z_1=C/Z_2
$$
を得る。従って $B/A=(Z_2-Z_1)/(Z_1+Z_2)$、$C/A=2Z_2/(Z_1+Z_2)$ である。エネルギー流束では速度の係数も掛かるので、圧力の透過係数が $1$ より大きくても保存則には反しない。実際 $(A^2-B^2)/Z_1=C^2/Z_2$ は二つの界面式を掛け合わせると従う。

:::example[圧力の反射率とエネルギー率の違い]
$Z_1=1,Z_2=3,A=2$ なら $B=1$、$C=3$ である。入射の流束は $4$、反射は左向きに $1$、透過は右向きに $3$。従って正味の流束は界面の両側で $3$ である。反射・透過のエネルギーの比は $1/4,3/4$ となり、圧力比 $1/2,3/2$ をそのまま流束比として読むことはできない。
:::

:::example[二成分が逆方向へ移す圧力と速度]
一次元音響で $K=8,\rho=2$ とする。右向きの圧力波を $4\cos(x-2t)$、左向きを $2\cos(x+2t)$ として重ね合わせる。

時間微分に質量行列 $M$ を掛け、空間微分に $B$ を掛ける。以下の表示で二つの微分を別々に計算する。和が所定の外力と一致することにより、成分ごとの方程式を検証できる。

$$
\begin{aligned}
M&=\left[\begin{matrix}\frac{1}{8} & 0\\0 & 2\end{matrix}\right]\\
B&=\left[\begin{matrix}0 & 1\\1 & 0\end{matrix}\right]\\
U&=\left[\begin{matrix}4 \cos{\left(2 t - x \right)} + 2 \cos{\left(2 t + x \right)}\\\cos{\left(2 t - x \right)} - \frac{\cos{\left(2 t + x \right)}}{2}\end{matrix}\right]\\
U_t&=\left[\begin{matrix}- 8 \sin{\left(2 t - x \right)} - 4 \sin{\left(2 t + x \right)}\\- 2 \sin{\left(2 t - x \right)} + \sin{\left(2 t + x \right)}\end{matrix}\right]\\
U_x&=\left[\begin{matrix}4 \sin{\left(2 t - x \right)} - 2 \sin{\left(2 t + x \right)}\\\sin{\left(2 t - x \right)} + \frac{\sin{\left(2 t + x \right)}}{2}\end{matrix}\right]\\
F&=\left[\begin{matrix}0\\0\end{matrix}\right]
\end{aligned}
$$

初期時刻へ戻して、先に指定したデータと比較する。次に正の時刻における成分の値を読む。成分が混ざっていても、初期データは方程式を積分する前の元の座標で指定される。

$$
\begin{aligned}
U(0)&=\left[\begin{matrix}6 \cos{\left(x \right)}\\\frac{\cos{\left(x \right)}}{2}\end{matrix}\right],\qquad t=0\\
U&=\left[\begin{matrix}2\\\frac{3}{2}\end{matrix}\right],\qquad t=\frac{\pi}{4},x=\frac{\pi}{2}\\
U_t&=\left[\begin{matrix}0\\0\end{matrix}\right],\qquad t=\frac{\pi}{4},x=\frac{\pi}{2}
\end{aligned}
$$

二次量を具体的に求め、微分した値と流束または外力の寄与を照合する。質量行列が単位行列でない場合は係数の普通の二乗和がエネルギーになるとは限らない。

$$
\begin{aligned}
e&=\frac{4 \cos^{2}{\left(2 t - x \right)} + \cos^{2}{\left(2 t + x \right)}}{2}\\
e_t&=\frac{\partial}{\partial t} \frac{4 \cos^{2}{\left(2 t - x \right)} + \cos^{2}{\left(2 t + x \right)}}{2}=- 4 \sin{\left(4 t - 2 x \right)} - \sin{\left(4 t + 2 x \right)}\\
q&=- \left(- 2 \cos{\left(2 t - x \right)} + \cos{\left(2 t + x \right)}\right) \left(2 \cos{\left(2 t - x \right)} + \cos{\left(2 t + x \right)}\right)\\
q_x&=\frac{\partial}{\partial x} \left(- \left(- 2 \cos{\left(2 t - x \right)} + \cos{\left(2 t + x \right)}\right) \left(2 \cos{\left(2 t - x \right)} + \cos{\left(2 t + x \right)}\right)\right)=4 \sin{\left(4 t - 2 x \right)} + \sin{\left(4 t + 2 x \right)}\\
e_t+q_x&=0\\
E(t)&=\int\limits_{0}^{2 \pi} \frac{4 \cos^{2}{\left(2 t - x \right)} + \cos^{2}{\left(2 t + x \right)}}{2}\, dx=\frac{5 \pi \sin^{2}{\left(2 t \right)}}{2} + \frac{5 \pi \cos^{2}{\left(2 t \right)}}{2}
\end{aligned}
$$

圧力の和に対し速度では左向きの寄与だけ符号が反対である。エネルギーを展開すると二方向の交差項は相殺し、流束は右向き圧力の二乗を $4$ で割った量から左向き圧力の二乗を $4$ で割った量を引く。この場合の周期積分は時刻に依らない。
:::

### 境界条件に必要な成分の数

滑らかな解が有界領域にあり、エネルギー恒等式の境界積分が存在するとする。外法線に対し $B_\nu=(A^0)^{-1/2}(\sum_j\nu_j A^j)(A^0)^{-1/2}$ とし、$W=(A^0)^{1/2}U$ を固有ベクトルに分解する。固有値が負の成分は外向き流束 $q\cdot\nu=\tfrac12 W^\top B_\nu W$ を負にするので、領域へ入る向きである。そこで負固有空間への射影を零にする境界条件なら $q\cdot\nu\ge0$、エネルギーの時間微分への寄与は $-q\cdot\nu\le0$ となる。二解の差に同じ条件を適用すれば、均一有界な係数の下で初期値零からの一意性を得る。この議論は一意性のための境界条件を示すものであり、曲がった境界や零固有値が交わる場合の存在まで自動的に与えない。

一次元定係数で $B_\nu$ が非特異なら、対角化した各成分は輸送方程式であり、入流の値と初期値を指定して実際に解ける。例えば $x>0$ の外法線は $-1$ なので、正の空間速度を持つ成分が入流になる。右向きの速度 $2$、左向きの速度 $-3$ の二成分では、左端で右向き成分だけに境界値を指定する。左向き成分の端点値は時刻 $t$ に初期位置 $3t$ から来る。従ってその値まで独立に決めると、もとの輸送問題の適合性を失う。

## 定数係数双曲型系と Fourier 解析

[Fourier 解析](/math/fourier-analysis/fourier-transform/) の規約 $\widehat{\partial_j f}=2\pi i\xi_j\widehat f$ に従う。$U_t+\sum_jB^j\partial_jU=0$ なら $B(\xi)=\sum_j\xi_jB^j$ に対して

$$
\widehat U(t,\xi)=e^{-2\pi it B(\xi)}\widehat U(0,\xi)
\label{eq-pde-hyperbolic-fourier}
$$

である。$B(\xi)$ が対称なら指数はユニタリで、Plancherel により $L^2$ ノルムを保つ。$A^0$ が一定正定値で対称化する場合は $(A^0)^{1/2}$ で共役すると同じ結論が $A^0$ 内積で成り立つ。一般の $B(\xi)$ でも実固有値を持ち、対角化行列 $R(\xi)$ の条件数 $\norm R\norm{R^{-1}}$ が方向について一様に有界なら、指数の作用素ノルムが一様有界である。よって Sobolev 空間でも同じ階数の連続依存性を得る。逆に実固有値だけでは Jordan ブロックを排除できない。

:::example[二つの実波]
一次元の $B=\begin{pmatrix}0&2\\2&0\end{pmatrix}$ は固有ベクトル $(1,1),(1,-1)$、固有値 $2,-2$ を持つ。左へ進む組合せと右へ進む組合せが分離し、$p+q=(p_0+q_0)(x-2t)$、$p-q=(p_0-q_0)(x+2t)$ となる。
:::

:::example[Jordan ブロックによる微分の損失]
$B=\begin{pmatrix}0&1\\0&0\end{pmatrix}$ に対して $B^2=0$ なので指数は $I-2\pi it\xi B$ である。従って $(p,q)=(p_0-tq_0',q_0)$ となる。$q_0$ の Fourier 台を $|\xi|\ge N$ に限定して $\norm{q_0}_2=1$ とすれば $\norm{q_0'}_2\ge2\pi N$ である。$p_0=0$ とし、$N\le|\xi|\le 2N$ に台を限定すれば滑らかな $L^2$ の初期値で試せる。初期ノルムが一定でも $t>0$ での第一成分のノルムは無界になり、同じ $L^2$ 階数で適切性が失われる。
:::

### 一様な Fourier 有界性から対称化を取り出す

定数係数行列の記号 $B(\theta)=\sum_j\theta_j B^j$ とする。ある固定の正時刻で解作用素が $L^2$ 有界と仮定する。Plancherel と周波数の小球に台を持つデータで試せば、乗数の行列ノルムが一様に有界であることが分かる。具体的には、ある周波数と単位ベクトルで上界を超えれば、連続性により同じベクトルを使って小球全体で上界を超えるため、$L^2$ 有界性に矛盾する。周波数の大きさと符号を動かすと $U_\theta(s)=e^{-isB(\theta)}$ が $s\in\R$ と $\theta$ に一様に有界、例えば $\norm U\le M$ である。

各 $\theta$ について $P_T=(2T)^{-1}\int_{-T}^T U_\theta(s)^*U_\theta(s)\dd s$ を取る。$U_\theta(s)^{-1}=U_\theta(-s)$ なので $M^{-2}I\le P_T\le M^2I$。有限次元で部分列極限 $P$ を取る。$U_\theta(r)^*P_TU_\theta(r)$ の積分区間は $r$ だけずれ、差はノルムで $M^2|r|/T$ 以下である。従って $U_\theta(r)^* P U_\theta(r)=P$。$r=0$ で微分して $PB=B^*P$ を得る。従って $P^{1/2}BP^{-1/2}$ は Hermite で、$B$ は実固有値を持ち、条件数が $M^2$ 以下の基底で対角化できる。こうして一つの正時刻での $L^2$ 有界性から、方向ごとの一様な対角化の必要性が得られる。

方向ごとに実固有値と一様な条件数の対角化を持つ性質を :term[強双曲性]{en="strong hyperbolicity" yomi="きょうそうきょくせい"} という。この定数係数の場合は上と合わせて、同じ Sobolev 階数で Cauchy 問題が適切であることと同値である。実際、行列指数そのものは周波数について連続であり、固有ベクトルの連続な選択は不要である。$H^s$ の Fourier ノルムに同じ乗数上界を使うと $S(t)g$ が定義でき、$t=0$ での強連続性は優収束から従う。外力が $F\in L^1(0,T;H^s)$ なら、各周波数で定数変化法を使い

$$
U(t)=S(t)g+\int_0^t S(t-r)F(r)\,\dd r,
\qquad
\norm{U(t)}_{H^s}\le M\left(\norm g_{H^s}+\int_0^t\norm F_{H^s}\right)
$$

を得る。滑らかなデータで式を確認してから、右の評価で一般のデータへ移せる。一方、対称双曲型では単一の正行列が全ての空間係数を対称化する。この二つの条件は区別して使う。

:::example[固有方向での Fourier 乗数]
一次元の $U_t+BU_x=0$ を一つの Fourier モードで考える。

法線記号の行列を $B$ とする。対角成分は座標方向での値を与える。固有方程式から直交する方向を取り出すため、特性多項式 $\det(\lambda I-B)$ を計算する。

$$
\begin{aligned}
B&=\left[\begin{matrix}\frac{62}{25} & - \frac{84}{25}\\- \frac{84}{25} & \frac{13}{25}\end{matrix}\right]\\
\det(\lambda I-B)&=\lambda^{2} - 3 \lambda - 10=\left(\lambda - 5\right) \left(\lambda + 2\right)\\
d_1&=-2,\quad v_1=\left[\begin{matrix}\frac{3}{5}\\\frac{4}{5}\end{matrix}\right]\\
d_2&=5,\quad v_2=\left[\begin{matrix}- \frac{4}{5}\\\frac{3}{5}\end{matrix}\right]
\end{aligned}
$$

$V$ の列に固有ベクトルを並べる。$V^\top V=I$ と $V^\top BV=D$ を使えば、元の座標 $X$ と固有座標 $Y=V^\top X$ の双方で状態ベクトルを記述できる。二次形式は固有座標にすると交差項を持たない。

$$
\begin{aligned}
V&=\left[\begin{matrix}\frac{3}{5} & - \frac{4}{5}\\\frac{4}{5} & \frac{3}{5}\end{matrix}\right]\\
V^\top BV&=\left[\begin{matrix}-2 & 0\\0 & 5\end{matrix}\right]\\
X^\top BX&=\frac{62 x_{1}^{2} - 168 x_{1} x_{2} + 13 x_{2}^{2}}{25}\\
\sum_j d_jY_j^2&=- 2 Y_{1}^{2} + 5 Y_{2}^{2}
\end{aligned}
$$

一つの状態ベクトルで固有係数を計算し、流束の二次形式を求める。直交変換なので二乗ノルムは変わらない。二次形式を二乗ノルムで割ると、固有値を非負の重みで平均した値になる。

$$
\begin{aligned}
X&=\left[\begin{matrix}1\\3\end{matrix}\right]\\
Y&=\left[\begin{matrix}3\\1\end{matrix}\right]\\
|X|^2&=10\\
X^\top BX&=-13\\
\frac{X^\top BX}{|X|^2}&=- \frac{13}{10}\\
\min_j d_j&=-2\\
\max_j d_j&=5
\end{aligned}
$$

周波数 $\xi$ で初期係数 $Y_j$ は $e^{-2\pi i t\xi d_j}Y_j$ に写る。この因子の絶対値は $1$ なのでノルムを保つ。実空間では二成分の固有射影がそれぞれ速度 $-2,5$ で異なる向きへ移動する。
:::

### Jordan ブロックで失う微分の階数

$B=\lambda I+N$、$N^{k+1}=0$、$N^k\ne0$ とする。Fourier 変換の乗数は
$$
e^{-2\pi i t\xi B}=e^{-2\pi i t\xi\lambda}\sum_{j=0}^k\frac{(-2\pi i t\xi)^j}{j!}N^j .
$$
$\lambda$ が実数なら指数因子の絶対値は $1$ である。固定した正の時刻では多項式の絶対値が $C_t(1+|\xi|^k)$ 以下なので $H^{s+k}\to H^s$ の有界性を得る。逆に $N^k e\ne0$ となる $e$ を取り、周波数 $R$ の近傍だけに Fourier 台を持つデータを $e$ 方向に選ぶ。一番高い項の大きさは $c_t R^k$ であり、低い項の総和は $O(R^{k-1})$ なので打ち消せない。従って $r<k$ なら $H^{s+r}\to H^s$ の一様有界性はない。固有値がすべて実という条件だけでは微分の損失が残り、その階数を Jordan の長さで読める。

:::example[三成分で二階の損失が起きる]
$u_t+v_x=0$、$v_t+w_x=0$、$w_t=0$ とする。初期値を $(0,0,g)$ とすれば順に $w=g$、$v=-tg'$、$u=t^2g''/2$ である。周期区間で $g(x)=\sin(nx)$ を入れると、第一成分の振幅は $t^2n^2/2$ になる。例えば $t=1$、$n=4$ なら $8$ であり、第三成分の初期振幅は $1$ である。三成分の全固有値は零だが、行列指数には二次項まで残る。
:::

:::example[冪零な結合の解は有限多項式になる]
$U_t+BU_x=0$、$B$ は上隣接成分が $1$ の三次 Jordan 行列とする。初期に第3成分だけ $\cos(2x)$ を置く。

時間微分に質量行列 $M$ を掛け、空間微分に $B$ を掛ける。以下の表示で二つの微分を別々に計算する。和が所定の外力と一致することにより、成分ごとの方程式を検証できる。

$$
\begin{aligned}
M&=\left[\begin{matrix}1 & 0 & 0\\0 & 1 & 0\\0 & 0 & 1\end{matrix}\right]\\
B&=\left[\begin{matrix}0 & 1 & 0\\0 & 0 & 1\\0 & 0 & 0\end{matrix}\right]\\
U&=\left[\begin{matrix}- 2 t^{2} \cos{\left(2 x \right)}\\2 t \sin{\left(2 x \right)}\\\cos{\left(2 x \right)}\end{matrix}\right]\\
U_t&=\left[\begin{matrix}- 4 t \cos{\left(2 x \right)}\\2 \sin{\left(2 x \right)}\\0\end{matrix}\right]\\
U_x&=\left[\begin{matrix}4 t^{2} \sin{\left(2 x \right)}\\4 t \cos{\left(2 x \right)}\\- 2 \sin{\left(2 x \right)}\end{matrix}\right]\\
F&=\left[\begin{matrix}0\\0\\0\end{matrix}\right]
\end{aligned}
$$

初期時刻へ戻して、先に指定したデータと比較する。次に正の時刻における成分の値を読む。成分が混ざっていても、初期データは方程式を積分する前の元の座標で指定される。

$$
\begin{aligned}
U(0)&=\left[\begin{matrix}0\\0\\\cos{\left(2 x \right)}\end{matrix}\right],\qquad t=0\\
U&=\left[\begin{matrix}-2\\0\\1\end{matrix}\right],\qquad t=1,x=0\\
U_t&=\left[\begin{matrix}-4\\0\\0\end{matrix}\right],\qquad t=1,x=0
\end{aligned}
$$

二次量を具体的に求め、微分した値と流束または外力の寄与を照合する。質量行列が単位行列でない場合は係数の普通の二乗和がエネルギーになるとは限らない。

$$
\begin{aligned}
e&=\frac{4 t^{4} \cos^{2}{\left(2 x \right)} + 4 t^{2} \sin^{2}{\left(2 x \right)} + \cos^{2}{\left(2 x \right)}}{2}\\
e_t&=\frac{\partial}{\partial t} \frac{4 t^{4} \cos^{2}{\left(2 x \right)} + 4 t^{2} \sin^{2}{\left(2 x \right)} + \cos^{2}{\left(2 x \right)}}{2}=4 t \left(2 t^{2} \cos^{2}{\left(2 x \right)} + \sin^{2}{\left(2 x \right)}\right)\\
E(t)&=\int\limits_{0}^{2 \pi} \frac{4 t^{4} \cos^{2}{\left(2 x \right)} + 4 t^{2} \sin^{2}{\left(2 x \right)} + \cos^{2}{\left(2 x \right)}}{2}\, dx=2 \pi t^{4} + 2 \pi t^{2} + \frac{\pi}{2}
\end{aligned}
$$

$B$ は対称でないので、ここで示した $|U|^2/2$ に対する対称系の流束等式は使えない。第3成分は変わらないが、第2成分は時間に一次、第1成分は二次で増える。空間周波数を増やすと第1成分の係数が周波数の二乗で増大する。
:::

## 演習問題

::::exercise{level=basic}
$(0,\pi)$ の $u_{tt}=u_{xx}$ で $g=0$、$h=3\sin2x$ の解を求めよ。

:::solution
係数は $a''+4a=0$、$a(0)=0$、$a'(0)=3$ なので $a=3\sin2t/2$ である。従って $u=3\sin2t\sin2x/2$ である。
:::
::::

::::exercise{level=basic}
$A=\diag(9,4)$ の一様な円錐速度を一つ求めよ。

:::solution
最大固有値は $9$ であり、$\nu^\top A\nu\le9$ である。:eqref[eq-pde-variable-wave-flux] で $c_*=3$ を選べる。
:::
::::

::::exercise{level=basic}
斉次波で $g=0$、$\norm h_H=2$ のとき $E(t)$ を求めよ。

:::solution
初期エネルギーは $\norm h_H^2/2=2$ である。外力の仕事がないため $E(t)=2$ で保存される。
:::
::::

::::exercise{level=basic}
$B=\diag(2,-1)$ の一次元系を初期値 $(p_0,q_0)$ から解け。

:::solution
成分ごとに $p_t+2p_x=0$、$q_t-q_x=0$ である。特性を逆にたどると $p=p_0(x-2t)$、$q=q_0(x+t)$ となる。
:::
::::

::::exercise{level=advanced}
$A$ の定数下界が $\lambda I$、$f\in L^1_tH$ のとき $\norm{Du(t)}_2$ を初期エネルギーで抑えよ。

:::solution
$\lambda\norm{Du(t)}_2^2\le a(u,u)\le2E(t)$ である。:eqref[eq-pde-hyperbolic-energy] により上界は $\lambda^{-1/2}(\sqrt{2E(0)}+\int_0^t\norm f_H)$ である。
:::
::::

::::exercise{level=advanced}
$A^0=\diag(4,1)$、$A^1=\begin{pmatrix}0&2\\2&0\end{pmatrix}$ の一階系で伝播速度を決めよ。

:::solution
正規化行列 $(A^0)^{-1/2}A^1(A^0)^{-1/2}$ は $\begin{pmatrix}0&1\\1&0\end{pmatrix}$ となる。固有値は $\pm1$ であり、最大絶対値 $1$ が円錐速度である。
:::
::::

::::exercise{level=advanced}
$2\times2$ Jordan ブロックの例で $q_0=\sin nx$、$p_0=0$ を周期領域で使う。$p(t)$ の $L^2$ ノルムを求めよ。

:::solution
$q=q_0$ なので $p_t=-n\cos nx$ から $p=-tn\cos nx$ である。$[0,2\pi]$ で $\int\cos^2nx=\pi$ だからノルムは $tn\sqrt\pi$（$t\ge0$、$n\ge1$）である。
:::
::::

::::exercise{level=advanced}
$A(t)= (1+t)I$ の滑らかな波で $f=0$ とする。$E(t)\le (1+t)E(0)$ を証明せよ。

:::solution
$E'=\frac12\norm{Du}_2^2\le E/(1+t)$ である。従って $(E/(1+t))'\le0$ であり、時間 $0$ から積分して結論を得る。
:::
::::

## 確認テスト

::::quiz
斉次固定端の波で自然な初期値空間はどれか。

- [x] 変位 $H_0^1$、速度 $L^2$
- [ ] 変位 $L^2$、速度 $H_0^1$ のみ
- [ ] 変位と速度が常に $C^\infty$

:::explanation
エネルギーは変位の一階弱微分と速度そのものの二乗積分からなる。従って変位に $H^1_0$、速度に $L^2$ を課す。
:::
::::

::::quiz{answer="3" tol="1e-6"}
$A=\diag(9,4)$ に対する円錐速度の一様上界 $\sqrt{\lambda_{\max}(A)}$ を求めよ。

:::explanation
最大固有値は $9$ であるため平方根は $3$ である。これを境界の移動速度にすると、流束から境界の移動による項を引いた寄与は非正になる。
:::
::::

::::quiz
実固有値を持つ一階定数係数系でも $L^2$ の適切性を失う理由はどれか。

- [x] Jordan ブロックで Fourier 指数に周波数のべきが現れる
- [ ] 全固有値が実なら常に指数減衰が起こる
- [ ] Fourier 変換が線形でない

:::explanation
非自明な冪零成分 $N$ について $e^{-2\pi it\xi N}$ は多項式となる。例えば $N^2=0$ では一次の項が現れ、周波数とともに乗数が無界になる。
:::
::::

::::quiz{answer="1" tol="1e-6"}
$A^0=\diag(4,1)$、$A^1_{12}=A^1_{21}=2$、他の成分は $0$ の系の最大伝播速度を求めよ。

:::explanation
正規化した行列は対角が $0$、非対角が $1$ なので固有値は $\pm1$ である。最大絶対値は $1$ である。
:::
::::

::::quiz{answer="4" tol="1e-6"}
波速 $c(t)=1+t$ の空間一様な波で、時刻 $2$ の観測点に対する時刻零の依存球の半径を求めよ。

:::explanation
境界を内側へ縮める速度を各時刻で $c(t)$ と取れば、半径は前向きに $R^{\prime}(t)=-c(t)$ を満たす。$R(2)=0$ から $R(0)=\int_0^2(1+t)\dd t=2+2=4$ である。
:::
::::


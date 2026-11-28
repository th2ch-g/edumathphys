---
title: 放物型方程式
description: 時間依存の Sobolev 空間と Galerkin 近似を用いて弱解を構成し、正則性と半群表示を導く。
notation:
  - symbol: 'V\subset H\subset V^*'
    meaning: 試験空間、基準 Hilbert 空間、双対空間
  - symbol: 'W(0,T)'
    meaning: 空間エネルギーと双対値の時間微分を持つ関数空間
  - symbol: 'e^{-tA}'
    meaning: 作用素 A に対応する熱半群
---

## この章の目標

- 双対空間に値を取る時間微分と初期値の意味を説明する。
- Galerkin 近似とエネルギー評価で弱解を作る。
- 自己共役な作用素の平滑化と時間微分を評価する。
- 弱最大値原理を証明し、半群表示を弱解と結び付ける。

## Bochner 空間と時間依存の Sobolev 空間

この章では $\Omega$ を有界開集合とし、関数は実数値とする。$H_0^1(\Omega)$ のノルムには $\norm{Du}_2$ を用いる。Poincaré の不等式 :ref[thm-fa-poincare] により通常の $H^1$ ノルムと同値である。空間微分の正則性と時間微分の正則性は一般に異なる。例えば $u\in H_0^1$ なら $-\Delta u$ は分布として $H^{-1}$ に入るが、必ずしも $L^2$ には入らない。従って時間微分を空間と同じノルムで測ると、自然な初期値が排除される。

:::definition[時間依存のエネルギー空間]{#def-pde-bochner-energy}
$V=H_0^1(\Omega)$、$H=L^2(\Omega)$ とし、$H$ を $(h,v)_H$ によって $V^*$ の部分空間とみなす。$1\le p<\infty$ とする。Banach 空間 $X$ に値を取る強可測関数で $\int_0^T\norm{u(t)}_X^p\,\dd t<\infty$ のものを、ほとんど至る所での一致で同一視して :term[Bochner 空間]{en="Bochner space" yomi="ぼほなーくうかん"} $L^p(0,T;X)$ とする。$p=\infty$ では $\norm{u(t)}_X$ の本質的上限が有限であるものを取る。時間微分 $u'$ は $V^*$ 値の分布微分として定義し、

$$
W(0,T)=\{u\in L^2(0,T;V):u'\in L^2(0,T;V^*)\}
$$

と置く。
:::

:::theorem[時間の跡とエネルギー恒等式]{#thm-pde-evolution-trace}
$u\in W(0,T)$ は $C([0,T];H)$ の代表を持ち、任意の $0\le s\le t\le T$ に対し

$$
\tfrac12\norm{u(t)}_H^2-\tfrac12\norm{u(s)}_H^2
=\int_s^t\langle u',u\rangle\,\dd r
\label{eq-pde-evolution-trace}
$$

が成り立つ。$u,v$ の内積に対する積の法則も成立する。
:::

:::proof
$u'\in L^2(V^*)$ なので $u$ は $V^*$ 値の絶対連続な代表を持つ。実際、$u'$ の Bochner 積分を引けば分布微分が零となり、双対との組ごとに定数である。従って端点で反射延長しても時間微分に跳躍項は出ず、$u\in L^2(V)$、$u'\in L^2(V^*)$ の性質を少し広い区間に延ばせる。時間方向の平滑化 $u_\eps$ は $L^2(V)$ で収束し、微分も $L^2(V^*)$ で収束する。滑らかな二つの近似の差 $w$ に対し、$t,s$ の順序によらず

$$
\norm{w(t)}_H^2\le\norm{w(s)}_H^2+2\norm{w'}_{L^2(V^*)}\norm w_{L^2(V)}
$$

である。$s$ で平均すると右辺の第一項は $T^{-1}\norm w_{L^2(H)}^2$ になる。よって近似は $C([0,T];H)$ で Cauchy となり、その極限は $u$ の代表である。積の恒等式は微分と関数の双対積の $L^1$ 収束を使えば極限に移せる。二つの関数の式は $u+v,u-v$ の二乗を引いて得られる。
:::

:::example[初期値がエネルギー空間の外にある場合]
$(0,\pi)$ の熱方程式で初期値 $g=1$ とする。定数は $L^2$ に属するが零トレースではない。解は正の時刻で $H_0^1$ に入り、初期値へは $L^2$ で近づく。実際、固有展開により $\int_0^T\sum_j\lambda_j e^{-2\lambda_jt}|g_j|^2\dd t\le\frac12\sum_j|g_j|^2$ であり、$L^2_tV$ は時刻 $0$ の $V$ 値を要求しない。
:::

:::example[時間微分が双対ノルムで収束する例]
$u_j(t)=e^{-j^2t}\sin jx$ を区間 $(0,\pi)$ で考える。$\sin jx$ の $H$ ノルムは $\sqrt{\pi/2}$ で、$V^*$ ノルムはこれを $j$ で割ったものである。従って $\int_0^T\norm{u_j'}_{V^*}^2\,\dd t=\pi(1-e^{-2j^2T})/4$ は一様に有界である。同じ計算で $H$ の時間微分の二乗積分にはさらに $j^2$ が掛かり、高周波で発散する。
:::

### 時間方向に失うノルムを具体化する

$H=L^2(0,\pi)$ で正規直交正弦基底を $e_j$ とする。$V$ での二乗ノルムは $\sum_j j^2|v_j|^2$、$V^*$ での二乗ノルムは $\sum_j j^{-2}|f_j|^2$ である。後者は Cauchy–Schwarz で上から評価し、有限個の係数について $v_j=j^{-2}f_j$ を選ぶと逆の評価が得られる。$u(t)=\sum_j e^{-j^2t}g_j e_j$ なら
$$
\int_0^T\norm u_V^2\dd t=\frac12\sum_j(1-e^{-2j^2T})|g_j|^2
=\int_0^T\norm{u'}_{V^*}^2\dd t.
$$
非負の和と積分を交換して各指数を積分した。例えば $g_j=1/j$ は $H$ に入るので両積分は有限だが、$u'$ の $L^2_t H$ ノルムを計算すると $j^2|g_j|^2$ に比例する和が発散する。正時刻ごとに $u'$ が滑らかであることと、零時刻まで時間に積分したノルムが有限であることは異なる。

### 時間端点の強い跡を近似で測る

$u,v$ が時間のエネルギー空間に属し、滑らかな時間近似で得た双線形恒等式を持つとする。任意の長さ $T$ の区間で、$w=u-v$ の $H$ ノルムの二乗の平均値を取る点 $s$ を選ぶ代わりに、$s$ での恒等式を $s$ について平均すれば
$$
\sup_t\norm{w(t)}_H^2
\le T^{-1}\int_0^T\norm w_H^2
+2\norm{w'}_{L^2(V^*)}\norm w_{L^2(V)}
$$
となる。実際、任意の $t,s$ について $\norm{w(t)}_H^2\le\norm{w(s)}_H^2+2\int_0^T|\langle w',w\rangle|$ とし、Cauchy–Schwarz を使う。この評価は滑らかな二つの近似の差にそのまま適用できるので、近似が $L^2(V)$ で収束し、その時間微分が $L^2(V^*)$ で収束すれば $C(H)$ で Cauchy となる。従って跡は弱収束の値としてだけでなく、強い $H$ 連続代表として定まる。

:::example[一点で跳ぶ係数は Bochner の時間微分に入らない]
$e\in V$ を非零とし、$u(t)=\mathbf1_{t>s}e$（$0<s<T$）を考える。$u$ は $L^2(V)$ に属するが、分布時間微分は $e\delta_s$ である。$v\in V$ を $(e,v)_H\ne0$ と選び、時刻 $s$ で $1$、台の幅 $\eps$ の滑らかな関数 $\eta_\eps$ を試験すると、$\langle u',v\eta_\eps\rangle=(e,v)_H$ は幅に依らない。もし $u'\in L^2(V^*)$ なら右辺は $\norm{u'}_{L^2(V^*)}\norm v_V\norm{\eta_\eps}_{L^2}$ 以下で、$\eps\downarrow0$ で零へ行くので矛盾する。従って空間的に滑らかな成分を一つ取っても、時刻の跳躍があれば跡の定理の仮定を満たさない。
:::

## Galerkin 近似

$a_t$ を $V$ 上の可測な実双線形形式とし、ほとんど全ての $t$ で

$$
|a_t(u,v)|\le M\norm u_V\norm v_V,\qquad
a_t(u,u)\ge \alpha\norm u_V^2-\beta\norm u_H^2,\quad \alpha>0,\ \beta\ge0
\label{eq-pde-parabolic-coercivity}
$$

とする。これは時間依存の発散形の楕円型作用素の変分形式を含む。$\{e_j\}$ を Dirichlet Laplace 作用素の $H$ 正規直交固有系とし、$V_N=\mathrm{span}(e_1,\dots,e_N)$ とする。射影 $P_N$ は $V$ でもノルムを増やさない。実際、:ref[thm-pde-elliptic-spectrum] のエネルギー展開で係数を切り詰めるだけである。

:::definition[Galerkin 近似]{#def-pde-galerkin}
$g\in H$、$f\in L^2(0,T;V^*)$ に対し $u_N(t)\in V_N$ を

$$
(u_N',v)_H+a_t(u_N,v)=\langle f,v\rangle\quad(v\in V_N),\qquad u_N(0)=P_Ng
\label{eq-pde-parabolic-galerkin}
$$

で定める。この有限次元化を :term[Galerkin 近似]{en="Galerkin approximation" yomi="がれーるきんきんじ"} と呼ぶ。
:::

係数については可測有界な行列の線形 ODE である。短区間の積分方程式を連続関数に作用させ、区間長と行列ノルムの積を $1$ より小さく取れば縮小写像になる。右辺は可積分なので解は絶対連続であり、有限個の区間をつなげば $[0,T]$ 上で得られる。$v=u_N$ と置けば

$$
\tfrac12(\norm{u_N}_H^2)'+\tfrac\alpha2\norm {u_N}_V^2
\le \beta\norm{u_N}_H^2+\tfrac{1}{2\alpha}\norm f_{V^*}^2
$$

を得る。積分因子 $e^{-2\beta t}$ を掛けると $\norm{u_N}_{L^\infty(H)}+\norm{u_N}_{L^2(V)}$ が $N$ によらず抑えられる。さらに $\langle u_N',v\rangle=(u_N',P_Nv)_H$ から $\norm{u_N'}_{V^*}\le \norm f_{V^*}+M\norm{u_N}_V$ である。これは時間微分の評価で射影の $V$ 有界性が必要な理由である。

:::example[時間変化する拡散率]
$(0,\pi)$ で $A(t)= (1+t)(-\partial_{xx})$、$g=2\sin x$、$f=0$ とする。どの $N\ge1$ でも第 1 モードだけがあり、$a'+(1+t)a=0$、$a(0)=2$ である。従って $u_N=2e^{-t-t^2/2}\sin x$ になり、この例では有限次元近似が厳密である。
:::

:::example[定数係数でないモード間の結合]
$(0,\pi)$ で $a(u,v)=\int(u'v'+(1+\cos x)uv)$ とする。正規化した $e_j=(2/\pi)^{1/2}\sin jx$ を使うと、$V_2$ 上で行列は $\begin{pmatrix}2&1/2\\1/2&5\end{pmatrix}$ である。実際、$\cos x\sin x=\sin2x/2$ なので非対角項は $1/2$ で、$\cos x$ を含む対角積分は $0$ である。初期値 $e_1$ でも $a_2'(0)=-1/2$ になり、第二モードが直ちに生じる。
:::

### 有限次元の質量行列

Galerkin 基底を一般の独立な $\psi_1,\ldots,\psi_N$ とすると、質量行列 $M_{ij}=(\psi_j,\psi_i)_H$ と剛性行列 $K_{ij}=a_t(\psi_j,\psi_i)$ により $M\dot c+Kc=F$ となる。$c^\top M c=\norm{\sum c_j\psi_j}_H^2>0$ なので $M$ は可逆である。$c^\top$ を掛けると
$$
\tfrac12\frac{\dd}{\dd t}(c^\top M c)+c^\top Kc=c^\top F .
$$
基底が正規直交でなくてもエネルギーは $c$ の Euclid ノルムではなく質量行列で測れば同じになる。$V^*$ 上で時間微分の一様上界まで得るには、近似空間への $H$ 射影が $V$ 上でも一様有界でなければならない。そのため存在証明では Laplace 基底が扱いやすい。

:::example[二つの基底が直交していない場合]
$\psi_1=e_1,\psi_2=e_1+e_2$ とすると $M=\left(\begin{smallmatrix}1&1\\1&2\end{smallmatrix}\right)$、$a(u,v)=\int Du\cdot Dv$ なら $K=\left(\begin{smallmatrix}1&1\\1&5\end{smallmatrix}\right)$ である。$\det(K-\lambda M)=(1-\lambda)(4-\lambda)$ なので時間の減衰率は $1,4$ のままである。$K$ 自体の固有値を減衰率として読むと、基底の取り方によって結果が変わってしまう。
:::

:::example[非対角の剛性を持つ二モード]
時間一定な Galerkin 方程式 $c_t+Kc=0$ の初期係数を固有方向に分ける。

Galerkin の剛性行列を $K$ とする。対角成分は座標方向での値を与える。固有方程式から直交する方向を取り出すため、特性多項式 $\det(\lambda I-K)$ を計算する。

$$
\begin{aligned}
K&=\left[\begin{matrix}\frac{73}{25} & - \frac{36}{25}\\- \frac{36}{25} & \frac{52}{25}\end{matrix}\right]\\
\det(\lambda I-K)&=\lambda^{2} - 5 \lambda + 4=\left(\lambda - 4\right) \left(\lambda - 1\right)\\
d_1&=1,\quad v_1=\left[\begin{matrix}\frac{3}{5}\\\frac{4}{5}\end{matrix}\right]\\
d_2&=4,\quad v_2=\left[\begin{matrix}- \frac{4}{5}\\\frac{3}{5}\end{matrix}\right]
\end{aligned}
$$

$V$ の列に固有ベクトルを並べる。$V^\top V=I$ と $V^\top KV=D$ を使えば、元の座標 $X$ と固有座標 $Y=V^\top X$ の双方で初期係数を記述できる。二次形式は固有座標にすると交差項を持たない。

$$
\begin{aligned}
V&=\left[\begin{matrix}\frac{3}{5} & - \frac{4}{5}\\\frac{4}{5} & \frac{3}{5}\end{matrix}\right]\\
V^\top KV&=\left[\begin{matrix}1 & 0\\0 & 4\end{matrix}\right]\\
X^\top KX&=\frac{73 x_{1}^{2} - 72 x_{1} x_{2} + 52 x_{2}^{2}}{25}\\
\sum_j d_jY_j^2&=Y_{1}^{2} + 4 Y_{2}^{2}
\end{aligned}
$$

一つの初期係数で固有係数を計算し、モードごとの散逸を求める。直交変換なので二乗ノルムは変わらない。二次形式を二乗ノルムで割ると、固有値を非負の重みで平均した値になる。

$$
\begin{aligned}
X&=\left[\begin{matrix}1\\3\end{matrix}\right]\\
Y&=\left[\begin{matrix}3\\1\end{matrix}\right]\\
|X|^2&=10\\
X^\top KX&=13\\
\frac{X^\top KX}{|X|^2}&=\frac{13}{10}\\
\min_j d_j&=1\\
\max_j d_j&=4
\end{aligned}
$$

固有座標での解は $Y_1e^{-t},Y_2e^{-4t}$ であり、元の座標には $V$ を掛けて戻す。$K$ の非対角項は元の基底での結合であるが、固定した直交変換でこの例の結合を解ける。
:::

:::example[空間近似を作ってから時刻依存の方程式に移す]
$(0,1)$ で $A=1$、$c=0$ とし、定常解が $u=x^2(1-x)$ となる源を入れる。有限次元の定常値から時間依存の方程式を組み立てる。

剛性行列を作る前に、方程式と試験空間を決める。ここでは全ての基底が両端で消える。厳密解を含まない試験空間に制限すると、近似は割り当てた方向の弱形式を満たすように決まる。厳密解は誤差を測るときに使う。

$$
\begin{aligned}
A&=1\\
u&=x^{2} \left(1 - x\right)\\
f&=2 \left(3 x - 1\right)\\
\psi_1&=x \left(1 - x\right)\\
\psi_2&=x^{2} \left(1 - x\right)^{2}
\end{aligned}
$$

剛性の成分は $K_{ij}=\int A\psi^{\prime}_j\psi^{\prime}_i+\int c\psi_j\psi_i$、右辺は $F_i=\int f\psi_i$ である。基底を正規直交には取っていないので、内積の質量行列 $M_{ij}=\int\psi_j\psi_i$ も別に求める。

$$
\begin{aligned}
K&=\left[\begin{matrix}\frac{1}{3} & \frac{1}{15}\\\frac{1}{15} & \frac{2}{105}\end{matrix}\right]\\
M&=\left[\begin{matrix}\frac{1}{30} & \frac{1}{140}\\\frac{1}{140} & \frac{1}{630}\end{matrix}\right]\\
F&=\left[\begin{matrix}\frac{1}{6}\\\frac{1}{30}\end{matrix}\right]\\
\det(K-\lambda M)&=\frac{\lambda^{2}}{529200} - \frac{\lambda}{4725} + \frac{1}{525}=\frac{\lambda^{2} - 112 \lambda + 1008}{529200}
\end{aligned}
$$

連立方程式 $Kc=F$ を解き、元の関数に戻す。残差 $r=f-Lu_N$ を計算すると一般に零ではないが、選んだ試験関数に対する積分は全て零になる。PDE を点で満たすことと有限次元の弱形式を満たすことを区別する。

$$
\begin{aligned}
c&=\left[\begin{matrix}\frac{1}{2}\\0\end{matrix}\right]\\
u_N&=- \frac{x \left(x - 1\right)}{2}\\
e&=u-u_N=- \frac{x \left(x - 1\right) \left(2 x - 1\right)}{2}\\
r&=3 \left(2 x - 1\right)\\
(\int r\psi_i)_i&=\left[\begin{matrix}0\\0\end{matrix}\right]
\end{aligned}
$$

最後に誤差の二乗ノルムを二通りに計算する。近似空間へのエネルギー直交性により、誤差の勾配エネルギーは全体のエネルギーから近似の分を引いた値になる。ここでは二次形式そのものを用いており、エネルギー汎関数の二分の一は付けていない。

$$
\begin{aligned}
a(u,u)&=\frac{2}{15}\\
a(u_N,u_N)&=\frac{1}{12}\\
a(e,e)&=\frac{1}{20}\\
\norm e_2^2&=\int\limits_{0}^{1} \frac{x^{2} \left(x - 1\right)^{2} \left(2 x - 1\right)^{2}}{4}\, dx=\frac{1}{840}
\end{aligned}
$$

熱方程式の二係数は $M c^{\prime}+Kc=F$ である。ここで得た $K^{-1}F$ は時間依存系の定常値である。これを引いた係数のエネルギー $c^\top M c/2$ は $-c^\top Kc$ の割合で減衰する。従って定常解へ近づく速さを求めるには $M$ と $K$ の両方が必要である。
:::

:::example[非直交な基底で同じ二モードを読む]
基底 $\psi_1=e_1,\psi_2=e_1+e_2$ を使い、初期係数 $c=(0,1)^\top$ とする。これは物理的には $e_1+e_2$ であり、固有値 $1,4$ がそれぞれに作用する。

各行の質量係数と剛性係数を区別して、候補を連立方程式に代入する。$B$ の非対角成分は別の成分を結び付ける。従って個々の微分の値を求めるだけではなく、全成分を掛け合わせた結果を確認する。

$$
\begin{aligned}
M&=\left[\begin{matrix}1 & 1\\1 & 2\end{matrix}\right]\\
K&=\left[\begin{matrix}1 & 1\\1 & 5\end{matrix}\right]\\
U&=\left[\begin{matrix}e^{- t} - e^{- 4 t}\\e^{- 4 t}\end{matrix}\right]\\
U_t&=\left[\begin{matrix}- e^{- t} + 4 e^{- 4 t}\\- 4 e^{- 4 t}\end{matrix}\right]\\
F&=\left[\begin{matrix}0\\0\end{matrix}\right]
\end{aligned}
$$

初期時刻へ戻して、先に指定したデータと比較する。次に正の時刻における成分の値を読む。成分が混ざっていても、初期データは方程式を積分する前の元の座標で指定される。

$$
\begin{aligned}
U(0)&=\left[\begin{matrix}0\\1\end{matrix}\right],\qquad t=0\\
U&=\left[\begin{matrix}\frac{7}{16}\\\frac{1}{16}\end{matrix}\right],\qquad t=\log{\left(2 \right)}\\
U_t&=\left[\begin{matrix}- \frac{1}{4}\\- \frac{1}{4}\end{matrix}\right],\qquad t=\log{\left(2 \right)}
\end{aligned}
$$

二次量を具体的に求め、微分した値と流束または外力の寄与を照合する。質量行列が単位行列でない場合は係数の普通の二乗和がエネルギーになるとは限らない。

$$
\begin{aligned}
e&=\frac{\left(e^{2 t} + 1\right) \left(e^{4 t} - e^{2 t} + 1\right) e^{- 8 t}}{2}\\
e_t&=\frac{d}{d t} \frac{\left(e^{2 t} + 1\right) \left(e^{4 t} - e^{2 t} + 1\right) e^{- 8 t}}{2}=- \left(e^{6 t} + 4\right) e^{- 8 t}\\
e_t&=\left(- e^{- t} - 4 e^{- 4 t}\right) e^{- 4 t} - \left(e^{- t} - e^{- 4 t}\right) e^{- t}=- \left(e^{6 t} + 4\right) e^{- 8 t}
\end{aligned}
$$

元の基底で第一係数は初期に零なのに正時刻で増える。これは物理量のエネルギー増加を意味しない。エネルギーには質量行列が入っており、固有空間の係数に戻せば二つとも減少する。
:::

## 弱解の存在と一意性

:::theorem[可測な時間依存形式の弱解]{#thm-pde-parabolic-existence}
:eqref[eq-pde-parabolic-coercivity] の条件のもと、$g\in H$、$f\in L^2(V^*)$ に対して、$u\in W(0,T)$ で $u(0)=g$ を満たし、

$$
\langle u',v\rangle+a_t(u,v)=\langle f,v\rangle\quad(v\in V)
\label{eq-pde-parabolic-weak}
$$

をほとんど全ての時刻で満たす解が一つ存在する。
:::

:::proof
前節の評価により部分列で $u_N\rightharpoonup u$ が $L^2(V)$、$u_N'\rightharpoonup w$ が $L^2(V^*)$ で成り立つ。時間の試験関数で積分すれば $w=u'$ である。固定した $e_j$ に時間試験関数を掛けて :eqref[eq-pde-parabolic-galerkin] を積分する。$N\ge j$ で式が成り立ち、$\int a_t(u_N,e_j)\eta(t)\dd t$ は $L^2(V)$ 上の有界線形汎関数を通じて収束する。稠密な可算個の試験ベクトルから全ての $v\in V$ に延ばせる。時間試験関数を $t=0$ で非零にすれば、部分積分で初期項 $(P_Ng,e_j)_H$ が残る。極限の端点値は :ref[thm-pde-evolution-trace] で定義でき、$(u(0),e_j)_H=(g,e_j)_H$ なので $u(0)=g$ である。一意性では差 $w$ に対して試験関数を $w(t)$ とし、$(\norm w_H^2)'\le 2\beta\norm w_H^2$ を得る。初期値が $0$ なので積分因子で $w=0$ が従う。
:::

同じ議論は連続依存性を与える。解の差 $w$ の初期値を $h$、右辺の差を $F$ とすれば $\norm w_{C(H)}^2+\norm w_{L^2(V)}^2+\norm{w'}_{L^2(V^*)}^2\le C(\norm h_H^2+\norm F_{L^2(V^*)}^2)$ である。定数は $T,\alpha,\beta,M$ に依存する。例えば $f=0$、$\beta=0$ なら $\norm{u(t)}_H$ は非増加である。

:::example[高周波初期値の差]
$\kappa=1$ の区間 $(0,\pi)$ で初期値の差を $\sin(3x)$ とする。時刻 $t$ の差は $e^{-9t}\sin(3x)$ であり、$H$ ノルムの比は $e^{-9t}$ になる。初期値の差が微分を持たない場合にも、定理は $H$ ノルムでの安定性を与える。
:::

:::example[時間に可測な外力と初期微分]
初期係数 $0$、固有値 $1$ のモードに $f(t)=\mathbf1_{t>1}$ を加える。$t\le1$ で $a=0$、$t>1$ で $a=1-e^{-(t-1)}$ である。$a$ は時刻 $1$ で連続で時間弱微分は有界だが、古典的時間微分はそこで跳ぶ。弱解の存在が外力の不連続な時刻での古典微分まで保証するわけではない。
:::

### 外力と初期値を同時に変えた場合

二解の差 $w$ に対し $\tfrac12(\norm w_H^2)'+\alpha\norm w_V^2\le\beta\norm w_H^2+\norm F_{V^*}\norm w_V$ である。右の積を $\alpha\norm w_V^2/2+\norm F_{V^*}^2/(2\alpha)$ と抑え、$e^{-2\beta t}$ を掛けて積分すると
$$
e^{-2\beta t}\norm{w(t)}_H^2+\alpha\int_0^t e^{-2\beta s}\norm w_V^2\dd s
\le\norm{w(0)}_H^2+\alpha^{-1}\int_0^t e^{-2\beta s}\norm F_{V^*}^2\dd s .
$$
また形式自体が異なる場合、$F=f_1-f_2+(A_2-A_1)u_2$ とすれば同じ評価が使える。作用素ノルム差を $\eta$ とすれば最後の項の $L^2(V^*)$ ノルムは $\eta\norm{u_2}_{L^2(V)}$ 以下である。従って存在証明で扱った弱解で、係数の一様近似に対する強い連続依存性が得られる。

:::example[切り替えの時刻に解は連続である]
$u'+a(t)u=0$、$a=1$（$t<1$）、$a=3$（$t>1$）、$u(0)=2$ なら、$u=2e^{-t}$ と $2e^{-1}e^{-3(t-1)}$ を $t=1$ でつなぐ。弱式に時刻 $1$ のデルタがないため値の跳躍はできないが、時間微分は $-2/e$ から $-6/e$ へ変わる。可測な時間係数の存在定理はこの差を許す。
:::

:::example[定常の源項を受ける二つのモード]
$(0,\pi)$ で $u_t-u_{xx}=\sin x+4\sin2x$、初期値零とする。弱解を二つのモードで構成し、エネルギーの収支を比較する。

まず初期値を固有関数ごとに分け、対応する時間方程式を解く。外力がある場合も、同じ固有関数に属する外力だけを時間方程式に入れる。以下では得られた表示の微分を取って、式の各項の関係を確かめる。

$$
\begin{aligned}
u&=\left(1 - e^{- 4 t}\right) \sin{\left(2 x \right)} + \left(1 - e^{- t}\right) \sin{\left(x \right)}\\
u_t&=\frac{\partial}{\partial t} \left(\left(1 - e^{- 4 t}\right) \sin{\left(2 x \right)} + \left(1 - e^{- t}\right) \sin{\left(x \right)}\right)=\left(e^{3 t} + 8 \cos{\left(x \right)}\right) e^{- 4 t} \sin{\left(x \right)}\\
u_{xx}&=\frac{\partial^{2}}{\partial x^{2}} \left(\left(1 - e^{- 4 t}\right) \sin{\left(2 x \right)} + \left(1 - e^{- t}\right) \sin{\left(x \right)}\right)=- \left(e^{t} - 1\right) \left(e^{3 t} \sin{\left(x \right)} + 4 e^{3 t} \sin{\left(2 x \right)} + 4 e^{2 t} \sin{\left(2 x \right)} + 4 e^{t} \sin{\left(2 x \right)} + 4 \sin{\left(2 x \right)}\right) e^{- 4 t}\\
f&=4 \left(1 - e^{- 4 t}\right) \sin{\left(2 x \right)} + \left(1 - e^{- t}\right) \sin{\left(x \right)} + e^{- t} \sin{\left(x \right)} + 4 e^{- 4 t} \sin{\left(2 x \right)}=\sin{\left(x \right)} + 4 \sin{\left(2 x \right)}
\end{aligned}
$$

境界条件を各モードで確かめてから、時刻零のデータに戻る。端点値がすべて零でも時間の係数は任意に選べない。時刻を零に戻して初めて各係数の初期条件が確認できる。

$$
\begin{aligned}
u(0,x)&=0,\qquad t=0\\
u(t,0)&=0,\qquad x=0\\
u(t,L)&=0,\qquad x=\pi
\end{aligned}
$$

エネルギーを空間積分で測る。異なるモードの積の積分は零であり、一つずつ積分して足し合わせられる。その時間微分を取り、解の式を使って積分する方法と比較する。

$$
\begin{aligned}
e(t,x)&=\frac{\left(\left(1 - e^{- 4 t}\right) \sin{\left(2 x \right)} + \left(1 - e^{- t}\right) \sin{\left(x \right)}\right)^{2}}{2}\\
\mathcal E(t)&=\frac{\pi \left(e^{t} - 1\right)^{2} \left(2 e^{6 t} + 2 e^{5 t} + 3 e^{4 t} + 4 e^{3 t} + 3 e^{2 t} + 2 e^{t} + 1\right) e^{- 8 t}}{4}\\
\mathcal E^{\prime}&=\frac{d}{d t} \frac{\pi \left(e^{t} - 1\right)^{2} \left(2 e^{6 t} + 2 e^{5 t} + 3 e^{4 t} + 4 e^{3 t} + 3 e^{2 t} + 2 e^{t} + 1\right) e^{- 8 t}}{4}=\frac{\pi \left(e^{t} - 1\right) \left(e^{6 t} + 4 e^{3 t} + 4 e^{2 t} + 4 e^{t} + 4\right) e^{- 8 t}}{2}\\
\text{入力と散逸の和}&=\int\limits_{0}^{\pi} \left(\left(\left(1 - e^{- 4 t}\right) \sin{\left(2 x \right)} + \left(1 - e^{- t}\right) \sin{\left(x \right)}\right) \left(\sin{\left(x \right)} + 4 \sin{\left(2 x \right)}\right) - \left(2 \left(1 - e^{- 4 t}\right) \cos{\left(2 x \right)} + \left(1 - e^{- t}\right) \cos{\left(x \right)}\right)^{2}\right)\, dx=\frac{\pi \left(e^{t} - 1\right) \left(e^{6 t} + 4 e^{3 t} + 4 e^{2 t} + 4 e^{t} + 4\right) e^{- 8 t}}{2}
\end{aligned}
$$

最後に一つの観測点で変位または温度と勾配を計算する。点で値が相殺されることと、区間全体のエネルギーが零になることは区別する。

$$
\begin{aligned}
u&=\frac{1}{2},\qquad t=\log{\left(2 \right)},x=\frac{\pi}{2}\\
u_x&=- \frac{15}{8},\qquad t=\log{\left(2 \right)},x=\frac{\pi}{2}\\
\mathcal E&=\frac{289 \pi}{1024},\qquad t=\log{\left(2 \right)}
\end{aligned}
$$

定常解は $\sin x+\sin2x$ である。各係数は定常値から指数的にずれ、$t=\log2$ では第1モードのずれが $1/2$、第2モードでは $1/16$ である。外力が可測なだけの場合にも積分因子の式は使え、有限次元近似では初期値を連続に保つ。
:::

:::example[可測係数を時間積分して初期値を保つ]
$(0,\pi)$ で時間だけに依存する拡散率 $\alpha(t)=1+2t$、第1モードの外力 $f_1=2(1+2t)$、初期係数 $0$ を用いる。

まず $u(t,x)=a(t)\sin(mx)$ を代入する。空間二階の微分は $-m^2$ 倍であり、時間微分は $a^{\prime}$ だけに作用する。従って時間依存係数を含む場合でも、他の空間モードを結合させなければ一次の方程式に帰着する。ここで $\mu(t)$ を空間作用素のこのモードでの係数とする。

$$
\begin{aligned}
m&=1\\
\mu(t)&=2 t + 1\\
f_m&=4 t + 2\\
M(t) &= \int_0^t \mu(s)\,\mathrm{d}s=\int\limits_{0}^{t} \left(2 z + 1\right)\, dz=t^{2} + t\\
I(t)&=e^{M(t)}=e^{t^{2} + t}
\end{aligned}
$$

積分因子を掛けると $(Ia)^{\prime}=If_m$ である。零時刻から積分し、初期係数を加えてから因子で割る。外力を単に時間積分するのと違い、積分の中で指数の重みが掛かる。

$$
\begin{aligned}
If_m&=2 \left(2 t + 1\right) e^{t} e^{t^{2}}\\
\int_0^t I(s)f_m(s)\,\mathrm{d}s&=\int\limits_{0}^{t} 2 \left(2 z + 1\right) e^{z} e^{z^{2}}\, dz=2 e^{t} e^{t^{2}} - 2\\
a(t)&=- 2 \left(-1 + e^{- t} e^{- t^{2}}\right)\\
a^{\prime}&=\frac{d}{d t} \left(- 2 \left(-1 + e^{- t} e^{- t^{2}}\right)\right)=2 \left(2 t + 1\right) e^{- t} e^{- t^{2}}\\
a^{\prime}+\mu a&=4 t e^{- t} e^{- t^{2}} - 2 \left(-1 + e^{- t} e^{- t^{2}}\right) \left(2 t + 1\right) + 2 e^{- t} e^{- t^{2}}=4 t + 2\\
a(0)&=0,\qquad t=0
\end{aligned}
$$

正弦の二乗積分は $\pi/2$ である。解の二乗ノルムの半分を時間微分し、散逸および入力から直接求めた変化率と比べる。係数が時間に依存しても、ここで微分する二乗ノルム自体には時間依存の作用素が入らない。

$$
\begin{aligned}
u(t,x)&=- 2 \left(-1 + e^{- t} e^{- t^{2}}\right) \sin{\left(x \right)}\\
E(t)&=\pi \left(e^{t} e^{t^{2}} - 1\right)^{2} e^{- 2 t} e^{- 2 t^{2}}\\
E^{\prime}&=\frac{d}{d t} \pi \left(e^{t} e^{t^{2}} - 1\right)^{2} e^{- 2 t} e^{- 2 t^{2}}=2 \pi \left(2 t + 1\right) \left(e^{t} e^{t^{2}} - 1\right) e^{- 2 t} e^{- 2 t^{2}}\\
\frac{\pi}2 a(f_m-\mu a)&=2 \pi \left(2 t + 1\right) \left(e^{t} e^{t^{2}} - 1\right) e^{- 2 t} e^{- 2 t^{2}}\\
u(t_*,x_*)&=2 - \frac{2}{e^{2}}=\frac{2 \left(-1 + e\right) \left(1 + e\right)}{e^{2}},\qquad t=1,x=\frac{\pi}{2}\\
E(t_*)&=\frac{\pi \left(-1 + e^{2}\right)^{2}}{e^{4}}=\frac{\pi \left(-1 + e\right)^{2} \left(1 + e\right)^{2}}{e^{4}},\qquad t=1
\end{aligned}
$$

積分因子の指数は $t+t^2$ なので、定常目標の値 $2$ に対し差は $-2e^{-(t+t^2)}$ である。一般の可測な正の係数でも原始関数は絶対連続で、同じ乗数の変形はほとんど至る所で成立する。ここで明示解が滑らかになるのは係数も滑らかだからである。
:::

## 正則性

:::proposition[対称かつ時間一定な形式での正則性]{#prop-pde-parabolic-regularity}
$a$ は対称で強圧的とし、その $H$ 上の自己共役作用素を $A$ とする。$g\in V$、$f\in L^2(0,T;H)$ なら $u'\in L^2(H)$、$Au\in L^2(H)$、$u\in C([0,T];V)$ である。
:::

:::proof
:ref[thm-pde-elliptic-spectrum] の $A$ の基底で Galerkin 近似を作る。このとき各係数は $u_j'+\lambda_j u_j=f_j$ に従う。$u_N'$ で試験すると

$$
\norm{u_N'}_H^2+\frac12(a(u_N,u_N))'=(P_N f,u_N')_H
\le\tfrac12\norm f_H^2+\tfrac12\norm{u_N'}_H^2 .
$$

積分し $N$ によらない $L^2(H)$ の時間微分と $C(V)$ の上界を得る。$Au_N=P_N f-u_N'$ も $L^2(H)$ で有界である。さらに二つの部分和の差に同じ評価を適用すると、初期値の $V$ 展開と右辺の $L^2(H)$ 展開の尾が $0$ に収束するため、$u_N$ は $C(V)$ で Cauchy であり、$u_N'$ と $Au_N$ は $L^2(H)$ で Cauchy である。各係数の極限を取ると最後の極限は $Au$ である。
:::

### 時間に Lipschitz な対称形式

$g\in V$、$f\in L^2(H)$ とし、$a_t$ は対称で一様強圧、形式ノルムについて時間に Lipschitz とする。この場合にも $u'\in L^2(H)$、$A(t)u\in L^2(H)$ を証明できる。ただし前節の固有基底は $t$ ごとに変わり得るため、近似には固定した Dirichlet Laplace 基底を用いる。近似の式に $u_N'$ を入れると a.e. に

$$
\norm{u_N'}_H^2+\tfrac12\partial_t a_t(u_N,u_N)
=(f,u_N')_H+\tfrac12\dot a_t(u_N,u_N).
$$

$L$ を形式の Lipschitz 定数とする。稠密な可算部分空間でスカラー関数 $a_t(v,w)$ を微分し、共通の充満測度集合上で $|\dot a_t(v,w)|\le L\norm v_V\norm w_V$ を使って $V$ 上へ延長できる。この微分は形式ノルムでの微分可能性を仮定しない。有限次元の積の法則と Young の不等式を使えば

$$
\int_0^t\norm{u_N'}_H^2+a_t(u_N(t),u_N(t))
\le a_0(P_Ng,P_Ng)+\int_0^t\norm f_H^2+L\int_0^t\norm{u_N}_V^2 .
$$

最後の項は前節のエネルギー評価で一様有界である。従って $u_N'$ は $L^2(H)$ で弱収束部分列を持ち、その極限は $u'$ である。:eqref[eq-pde-parabolic-weak] から $A(t)u=f-u'$ は $L^2(H)$ に入る。この証明では固定基底への射影と $A(t)$ は可換と仮定していない。強圧性と近似の一様 $L^\infty(V)$ 上界から $u$ も $L^\infty(V)$ に入る。

例えば $a_t(u,v)=(1+t)\int Du\cdot Dv$ なら $\dot a_t(u,u)=\norm{Du}_2^2$ であり、時間依存項は主形式以下である。一方、単に時刻について可測な形式ではこの時間微分の上界を使えない。その場合に前節で得た $u'\in L^2(V^*)$ と区別する。

$C^{1,1}$ 領域で、空間について一様な Lipschitz ノルムを持つ発散形の係数なら :ref[thm-pde-elliptic-boundary-h2] を a.e. の時刻に適用して $L^2_tH^2_x$ が従う。逆に時間依存の一般の形式では $a_t(u,u)$ の微分に $\partial_ta_t$ の項が現れるので、上の証明から同じ評価は出ない。自律な対称作用素の斉次問題では任意の $g\in H$ について $\norm{A^m u(t)}_H\le (m/(et))^m\norm g_H$（$m\ge1$）である。これは $\lambda^m e^{-t\lambda}$ の最大が $\lambda=m/t$ で達成されるためである。

:::example[モード別の平滑化]
$(0,\pi)$ の Dirichlet 熱作用素で初期値 $g_j=1/j$ を選ぶと $g\in H$ だが $\sum_j j^2|g_j|^2=\infty$ なので $g\notin V$ である。任意の $t>0$ では $\sum_j j^{4m-2}e^{-2j^2t}$ が有限で、$A^m u(t)\in H$ となる。
:::

:::example[初期値と生成項の正則性の差]
固有値 $\lambda$ の一つのモードに $f_j=1$ を加え零初期値を課すと、$u_j=(1-e^{-\lambda t})/\lambda$ である。$Au$ のその係数は $1-e^{-\lambda t}$ で、外力を打ち消さない時刻には $\lambda$ に関して指数減衰しない。斉次問題で初期値に適用した任意の高階の平滑化を、粗い外力のある問題にそのまま適用できない。
:::

### 自律解の初期時刻での正則性の逆

正の自己共役作用素 $A$ の固有値 $\lambda_j>0$、$g\in H$ に対し $u(t)=\sum e^{-t\lambda_j}g_j e_j$ とする。単調収束で
$$
\int_0^T\norm{Au(t)}_H^2\dd t
=\frac12\sum_j\lambda_j(1-e^{-2T\lambda_j})|g_j|^2 .
$$
最低固有値が正なので $1-e^{-2T\lambda_j}$ は零から離れる。従って積分が有限であることは $g\in D(A^{1/2})$ と同値である。$u'=-Au$ にも同じ判定が使える。拡散の平滑化は $t>0$ で任意階に及ぶが、時刻零まで含める最大正則性は初期値に半階分を要求する。

:::example[初期値の指数で境界が分かれる]
Dirichlet Laplace 作用素で $\lambda_j=j^2$、$g_j=j^{-s}$ なら $g\in H$ の条件は $s>1/2$、$Au\in L^2_tH$ の条件は $2-2s<-1$、すなわち $s>3/2$ である。$s=3/2$ のとき右辺は調和級数に比較できるので発散する。$s=2$ なら和は $\tfrac12\sum j^{-2}(1-e^{-2Tj^2})$ で有限である。
:::

### 初期時刻から離れた部分での評価

対称な形式が時間について Lipschitz で、$g\in H$、$f\in L^2(0,T;H)$ とする。初期値の $V$ 正則性を仮定しなくても、任意の $\tau>0$ に対して $u',A(t)u\in L^2(\tau,T;H)$ である。実際、$t=0$ の近くで零、$t\ge\tau$ で $1$ となる滑らかな $\eta$ を取り、$v=\eta u$ と置く。$v(0)=0\in V$ であり、弱式の積の法則から $v'+A(t)v=\eta f+\eta' u$ である。$u\in C([0,T];H)$ なので右辺は $L^2(H)$ に入る。従って上で証明した正則性を $v$ に使える。

さらに $\Omega$ が $C^{1,1}$、主係数が空間で一様に Lipschitz なら、:ref[thm-pde-elliptic-boundary-h2] により $u\in L^2(\tau,T;H^2(\Omega))$ となる。例えば区間の熱方程式では $g_j=1/j$ のような初期値にも使える。実際、この場合の $u'$ の二乗ノルムを $(\tau,T)$ で積分すると $\frac12\sum_j(e^{-2\tau j^2}-e^{-2Tj^2})$ となる。$\tau>0$ なら有限だが、$\tau\downarrow0$ では発散する。

:::example[作用素を一回掛けた後の平滑化]
三基底の自律問題で、固有値が大きく異なる模型を考える。

Galerkin の剛性行列を $K$ とする。対角成分は座標方向での値を与える。固有方程式から直交する方向を取り出すため、特性多項式 $\det(\lambda I-K)$ を計算する。

$$
\begin{aligned}
K&=\left[\begin{matrix}\frac{49601}{4225} & - \frac{34032}{4225} & \frac{720}{169}\\- \frac{34032}{4225} & \frac{29749}{4225} & - \frac{540}{169}\\\frac{720}{169} & - \frac{540}{169} & \frac{1220}{169}\end{matrix}\right]\\
\det(\lambda I-K)&=\lambda^{3} - 26 \lambda^{2} + 125 \lambda - 100=\left(\lambda - 20\right) \left(\lambda - 5\right) \left(\lambda - 1\right)\\
d_1&=1,\quad v_1=\left[\begin{matrix}\frac{3}{5}\\\frac{4}{5}\\0\end{matrix}\right]\\
d_2&=5,\quad v_2=\left[\begin{matrix}- \frac{4}{13}\\\frac{3}{13}\\\frac{12}{13}\end{matrix}\right]\\
d_3&=20,\quad v_3=\left[\begin{matrix}\frac{48}{65}\\- \frac{36}{65}\\\frac{5}{13}\end{matrix}\right]
\end{aligned}
$$

$V$ の列に固有ベクトルを並べる。$V^\top V=I$ と $V^\top KV=D$ を使えば、元の座標 $X$ と固有座標 $Y=V^\top X$ の双方で初期係数を記述できる。二次形式は固有座標にすると交差項を持たない。

$$
\begin{aligned}
V&=\left[\begin{matrix}\frac{3}{5} & - \frac{4}{13} & \frac{48}{65}\\\frac{4}{5} & \frac{3}{13} & - \frac{36}{65}\\0 & \frac{12}{13} & \frac{5}{13}\end{matrix}\right]\\
V^\top KV&=\left[\begin{matrix}1 & 0 & 0\\0 & 5 & 0\\0 & 0 & 20\end{matrix}\right]\\
X^\top KX&=\frac{49601 x_{1}^{2} - 68064 x_{1} x_{2} + 36000 x_{1} x_{3} + 29749 x_{2}^{2} - 27000 x_{2} x_{3} + 30500 x_{3}^{2}}{4225}\\
\sum_j d_jY_j^2&=Y_{1}^{2} + 5 Y_{2}^{2} + 20 Y_{3}^{2}
\end{aligned}
$$

一つの初期係数で固有係数を計算し、モードごとの散逸を求める。直交変換なので二乗ノルムは変わらない。二次形式を二乗ノルムで割ると、固有値を非負の重みで平均した値になる。

$$
\begin{aligned}
X&=\left[\begin{matrix}1\\1\\2\end{matrix}\right]\\
Y&=\left[\begin{matrix}\frac{7}{5}\\\frac{23}{13}\\\frac{62}{65}\end{matrix}\right]\\
|X|^2&=6\\
X^\top KX&=\frac{151286}{4225}\\
\frac{X^\top KX}{|X|^2}&=\frac{75643}{12675}\\
\min_j d_j&=1\\
\max_j d_j&=20
\end{aligned}
$$

時刻 $t>0$ で $Kc(t)$ の第 $j$ 係数は $d_je^{-td_j}Y_j$ であり、$d_je^{-td_j}\le1/(et)$ である。$z e^{-tz}$ の最大は $z=1/t$ にあるためである。同じ上界を全ての近似空間で使えるので、高周波数と無関係な作用素ノルムの評価になる。
:::

:::example[移り変わる係数と時間微分の強い評価]
$(0,\pi)$ で $\alpha(t)=1/(1+t)$、第2モードに $f_2(t)=4/(1+t)$ を加え、初期係数 $3$ とする。

まず $u(t,x)=a(t)\sin(mx)$ を代入する。空間二階の微分は $-m^2$ 倍であり、時間微分は $a^{\prime}$ だけに作用する。従って時間依存係数を含む場合でも、他の空間モードを結合させなければ一次の方程式に帰着する。ここで $\mu(t)$ を空間作用素のこのモードでの係数とする。

$$
\begin{aligned}
m&=2\\
\mu(t)&=\frac{4}{t + 1}\\
f_m&=\frac{4}{t + 1}\\
M(t) &= \int_0^t \mu(s)\,\mathrm{d}s=\int\limits_{0}^{t} \frac{4}{z + 1}\, dz=4 \log{\left(t + 1 \right)}\\
I(t)&=e^{M(t)}=\left(t + 1\right)^{4}
\end{aligned}
$$

積分因子を掛けると $(Ia)^{\prime}=If_m$ である。零時刻から積分し、初期係数を加えてから因子で割る。外力を単に時間積分するのと違い、積分の中で指数の重みが掛かる。

$$
\begin{aligned}
If_m&=4 \left(t + 1\right)^{3}\\
\int_0^t I(s)f_m(s)\,\mathrm{d}s&=\int\limits_{0}^{t} 4 \left(z + 1\right)^{3}\, dz=t^{4} + 4 t^{3} + 6 t^{2} + 4 t\\
a(t)&=\frac{t^{4} + 4 t^{3} + 6 t^{2} + 4 t + 3}{\left(t + 1\right)^{4}}\\
a^{\prime}&=\frac{d}{d t} \frac{t^{4} + 4 t^{3} + 6 t^{2} + 4 t + 3}{\left(t + 1\right)^{4}}=- \frac{8}{\left(t + 1\right)^{5}}\\
a^{\prime}+\mu a&=\frac{4 t^{3} + 12 t^{2} + 12 t + 4}{\left(t + 1\right)^{4}}=\frac{4}{t + 1}\\
a(0)&=3,\qquad t=0
\end{aligned}
$$

正弦の二乗積分は $\pi/2$ である。解の二乗ノルムの半分を時間微分し、散逸および入力から直接求めた変化率と比べる。係数が時間に依存しても、ここで微分する二乗ノルム自体には時間依存の作用素が入らない。

$$
\begin{aligned}
u(t,x)&=\frac{\left(t^{4} + 4 t^{3} + 6 t^{2} + 4 t + 3\right) \sin{\left(2 x \right)}}{\left(t + 1\right)^{4}}\\
E(t)&=\frac{\pi \left(t^{4} + 4 t^{3} + 6 t^{2} + 4 t + 3\right)^{2}}{4 \left(t + 1\right)^{8}}\\
E^{\prime}&=\frac{d}{d t} \frac{\pi \left(t^{4} + 4 t^{3} + 6 t^{2} + 4 t + 3\right)^{2}}{4 \left(t + 1\right)^{8}}=- \frac{4 \pi \left(t^{4} + 4 t^{3} + 6 t^{2} + 4 t + 3\right)}{\left(t + 1\right)^{9}}\\
\frac{\pi}2 a(f_m-\mu a)&=- \frac{4 \pi \left(t^{4} + 4 t^{3} + 6 t^{2} + 4 t + 3\right)}{\left(t + 1\right)^{9}}\\
u(t_*,x_*)&=\frac{9}{8},\qquad t=1,x=\frac{\pi}{4}\\
E(t_*)&=\frac{81 \pi}{256},\qquad t=1
\end{aligned}
$$

有限時間では $\alpha$ が零から離れ、Lipschitz である。実際 $a^{\prime}=-8/(1+t)^5$、$A(t)u$ の係数は $4a/(1+t)$ で、ともに有限時間の $L^2$ に入る。$a_t(u,u)$ を微分する場合は、この計算の $E(t)$ と異なり拡散率の微分 $\alpha^{\prime}$ が追加される。
:::

## 放物型方程式の最大値原理

:::proposition[弱解の正部分による比較]{#prop-pde-parabolic-positive-part}
$A(t,x)$ は有界可測かつ一様楕円、$c$ は有界可測で非負とし、$u\in W(0,T)$ が $u'-\divg(A Du)+cu=f$ を満たすとする。$u(0)\le0$ で $f(t)$ が全ての非負試験関数に非正なら、$u(t)\le0$ が全時刻に従う。
:::

:::proof
時刻で平滑化した関数では $\frac12(\norm{u_+}_H^2)'=(u',u_+)_H$ である。平滑化を $W$ で戻すと連続埋め込みにより端点の正部分は $H$ で収束する。$u\mapsto u_+$ は $L^2(V)$ でも強収束を保つ。実際 $Du_+=\mathbf1_{\{u>0\}}Du$ であり、$u_j\to u$ の a.e. 部分列に対して指示関数は $u\ne0$ の範囲で収束し、$\{u=0\}$ では $Du=0$ a.e. である。勾配の差を $\mathbf1_{\{u_j>0\}}(Du_j-Du)+(\mathbf1_{\{u_j>0\}}-\mathbf1_{\{u>0\}})Du$ と分ければ結論を得る。このため双対積を極限に移せる。弱式を $u_+(t)$ で試験すると $\langle u',u_+\rangle\le-\lambda\norm{Du_+}_2^2$ である。積分恒等式で $\norm{u_+(t)}_H^2\le\norm{u_+(0)}_H^2=0$ となる。
:::

二解の差に適用すれば、初期値と外力の順序を解が保つ。例えば $f\ge0$、$g\ge0$ の問題では $-u$ に適用して $u\ge0$ である。古典解が強最大値原理を満たす場合とは異なり、この議論では特定の内点での厳密な正値性は主張しない。全空間の定数係数の解では熱核の正値性により、非負かつ零でない初期値から任意の $t>0$ で各点の正値性が従う。

:::example[源と減衰がある一つのモード]
$(0,\pi)$ で $u_t-u_{xx}+2u=\sin x$、初期値 $0$ とする。係数は $a'+3a=1$ であり $a=(1-e^{-3t})/3$ である。従って $u$ は正の時刻の領域内部で正となり、弱比較が与える非負性とも一致する。
:::

:::example[正部分と負の外力]
$(0,\pi)$ で零初期値、零境界値、$f=-2\sin x$ とすると係数は $a=-2(1-e^{-t})$ である。従って $u_+=0$ だが $\norm{u(t)}_H$ は増える。上からの比較とノルムの単調性は別の結論であり、後者には外力の仕事項がないことを使う。
:::

### 反応項があるときの上下解

$u_t-\divg(A Du)+b\cdot Du+c u=f$ で $b,c$ は有界とする。劣解と優解の差を $w$ とし、差の右辺と初期値を非正とする。零境界の差なら $w\in W(0,T)$ とする。非零の境界値を比べる場合は有界 Lipschitz 領域で $w\in L^2(0,T;H^1(\Omega))\cap C([0,T];H)$、$w'\in L^2(0,T;V^*)$、境界のトレースが非正とする。この場合も $w_+\in L^2(0,T;V)$ である。前の証明で時間平滑化に非負な核と端点での反射延長を使えば、平滑化した差のトレースも非正であり、その正部分を試験に使える。$H^1$ での正部分の強収束によって同じ連鎖律が従う。

$w_+$ を試験すると、ドリフトは絶対値で $\norm b_\infty\norm{Dw_+}_2\norm{w_+}_2$ 以下である。Young の不等式で勾配の強圧項の半分を吸収し、負の $c$ も $\norm{c_-}_\infty\norm{w_+}_2^2$ で抑える。従って $(\norm{w_+}_2^2)'\le C\norm{w_+}_2^2$ であり、積分因子と初期値の順序から比較が従う。

:::example[減衰率より増幅率の方が強い場合]
$(0,\pi)$ で $u_t-u_{xx}-3u=0$、$g=\sin x$ とすると $u=e^{2t}\sin x$ である。初期値が非負なので順序は保つが、最大値は増える。時間因子を付けた $w=e^{-3t}u=e^{-t}\sin x$ は同次熱方程式になり、一様ノルムが非増加である。比較と一様ノルムの縮小が別の性質であることが分かる。
:::

:::example[反応が正値性と時間増大へ与える違い]
$(0,\pi)$ で $u_t-u_{xx}-6u=0$、初期値 $2\sin x$ とする。零との比較と振幅の増大をあわせて調べる。

まず初期値を固有関数ごとに分け、対応する時間方程式を解く。外力がある場合も、同じ固有関数に属する外力だけを時間方程式に入れる。以下では得られた表示の微分を取って、式の各項の関係を確かめる。

$$
\begin{aligned}
u&=2 e^{5 t} \sin{\left(x \right)}\\
u_t&=\frac{\partial}{\partial t} 2 e^{5 t} \sin{\left(x \right)}=10 e^{5 t} \sin{\left(x \right)}\\
u_{xx}&=\frac{\partial^{2}}{\partial x^{2}} 2 e^{5 t} \sin{\left(x \right)}=- 2 e^{5 t} \sin{\left(x \right)}\\
f&=0
\end{aligned}
$$

境界条件を各モードで確かめてから、時刻零のデータに戻る。端点値がすべて零でも時間の係数は任意に選べない。時刻を零に戻して初めて各係数の初期条件が確認できる。

$$
\begin{aligned}
u(0,x)&=2 \sin{\left(x \right)},\qquad t=0\\
u(t,0)&=0,\qquad x=0\\
u(t,L)&=0,\qquad x=\pi
\end{aligned}
$$

エネルギーを空間積分で測る。異なるモードの積の積分は零であり、一つずつ積分して足し合わせられる。その時間微分を取り、解の式を使って積分する方法と比較する。

$$
\begin{aligned}
e(t,x)&=2 e^{10 t} \sin^{2}{\left(x \right)}\\
\mathcal E(t)&=\pi e^{10 t}\\
\mathcal E^{\prime}&=\frac{d}{d t} \pi e^{10 t}=10 \pi e^{10 t}\\
\text{入力と散逸の和}&=\int\limits_{0}^{\pi} \left(24 e^{10 t} \sin^{2}{\left(x \right)} - 4 e^{10 t} \cos^{2}{\left(x \right)}\right)\, dx=10 \pi e^{10 t}
\end{aligned}
$$

最後に一つの観測点で変位または温度と勾配を計算する。点で値が相殺されることと、区間全体のエネルギーが零になることは区別する。

$$
\begin{aligned}
u&=2,\qquad t=\frac{\log{\left(2 \right)}}{5},x=\frac{\pi}{6}\\
u_x&=2 \sqrt{3},\qquad t=\frac{\log{\left(2 \right)}}{5},x=\frac{\pi}{6}\\
\mathcal E&=4 \pi,\qquad t=\frac{\log{\left(2 \right)}}{5}
\end{aligned}
$$

反応係数を引いた $w=e^{-6t}u$ は通常の熱方程式に従う。従って元の解は正値性を保つが、定数の上界が時間によらず保たれるという結論は使えない。エネルギーの増加は $6\norm u_2^2$ の反応が $\norm{u_x}_2^2$ の散逸を上回ることに対応する。
:::

:::example[系では成分ごとの最大値原理が自動ではない]
$M=I$、$K=\left(\begin{smallmatrix}2&1\\1&2\end{smallmatrix}\right)$ とし、初期係数 $(1,0)^\top$ を考える。正の初期成分から負の成分が生じるか調べる。

各行の質量係数と剛性係数を区別して、候補を連立方程式に代入する。$B$ の非対角成分は別の成分を結び付ける。従って個々の微分の値を求めるだけではなく、全成分を掛け合わせた結果を確認する。

$$
\begin{aligned}
M&=\left[\begin{matrix}1 & 0\\0 & 1\end{matrix}\right]\\
K&=\left[\begin{matrix}2 & 1\\1 & 2\end{matrix}\right]\\
U&=\left[\begin{matrix}\frac{e^{- t}}{2} + \frac{e^{- 3 t}}{2}\\- \frac{e^{- t}}{2} + \frac{e^{- 3 t}}{2}\end{matrix}\right]\\
U_t&=\left[\begin{matrix}- \frac{e^{- t}}{2} - \frac{3 e^{- 3 t}}{2}\\\frac{e^{- t}}{2} - \frac{3 e^{- 3 t}}{2}\end{matrix}\right]\\
F&=\left[\begin{matrix}0\\0\end{matrix}\right]
\end{aligned}
$$

初期時刻へ戻して、先に指定したデータと比較する。次に正の時刻における成分の値を読む。成分が混ざっていても、初期データは方程式を積分する前の元の座標で指定される。

$$
\begin{aligned}
U(0)&=\left[\begin{matrix}1\\0\end{matrix}\right],\qquad t=0\\
U&=\left[\begin{matrix}\frac{5}{16}\\- \frac{3}{16}\end{matrix}\right],\qquad t=\log{\left(2 \right)}\\
U_t&=\left[\begin{matrix}- \frac{7}{16}\\\frac{1}{16}\end{matrix}\right],\qquad t=\log{\left(2 \right)}
\end{aligned}
$$

二次量を具体的に求め、微分した値と流束または外力の寄与を照合する。質量行列が単位行列でない場合は係数の普通の二乗和がエネルギーになるとは限らない。

$$
\begin{aligned}
e&=\frac{\left(e^{4 t} + 1\right) e^{- 6 t}}{4}\\
e_t&=\frac{d}{d t} \frac{\left(e^{4 t} + 1\right) e^{- 6 t}}{4}=- \frac{\left(e^{4 t} + 3\right) e^{- 6 t}}{2}\\
e_t&=\left(- \frac{e^{- t}}{2} - \frac{3 e^{- 3 t}}{2}\right) \left(\frac{e^{- t}}{2} + \frac{e^{- 3 t}}{2}\right) + \left(- \frac{e^{- t}}{2} + \frac{e^{- 3 t}}{2}\right) \left(\frac{e^{- t}}{2} - \frac{3 e^{- 3 t}}{2}\right)=- \frac{\left(e^{4 t} + 3\right) e^{- 6 t}}{2}
\end{aligned}
$$

正時刻では第2成分が負であるが、二乗エネルギーは減少する。単一のスカラー PDE の正部分の証明を各成分に適用すると、他成分を含む項の符号が決まらない。このような連立反応系で順序保存を示すには、成分の結合に別の符号条件が要る。
:::

### 古典解における正値性の伝播

連結な開集合 $\Omega$ と時間区間を考える。$P=\partial_t-A:D^2-b\cdot D+c$ の係数は有界、$A$ は一様楕円とする。非負の古典的優解 $u$ が $Pu\ge0$ と $u(s,y)>0$ を満たせば、任意の $t_1>s$ と $x\in\Omega$ で $u(t_1,x)>0$ となる。$u$ と必要な微分は時空間の内部で連続とする。

実際、$y$ から $x$ へ至る滑らかな経路 $\gamma(t)$（$s\le t\le t_1$）を取り、全経路のまわりに半径 $R$ の閉球を置いても $\Omega$ に収まるようにする。さらに $R$ を小さくして $u(s,\cdot)$ が初めの球で正の下限を持つようにする。$\rho(t,z)=R^2-|z-\gamma(t)|^2$ とし、

$$
w(t,z)=\delta e^{-K(t-s)}(\rho(t,z)_+)^3
$$

を比較関数にする。この関数は時間一階、空間二階まで連続である。球内で $Pw/(\delta e^{-K(t-s)})$ は

$$
(c-K)\rho^3
+3\rho^2(\rho_t-A:D^2\rho-b\cdot D\rho)
-6\rho\,D\rho^\top A D\rho
$$

となる。球面に近い所では $|D\rho|$ が零から離れるため、最後の負の項が中項を抑える。球の残りでは $\rho$ が正の定数以上なので、$K$ を大きくすれば全体が非正になる。球の外では $w=0$ なので、全体で $Pw\le0$ である。$\delta>0$ を小さく選んで初期値を $w(s,\cdot)\le u(s,\cdot)$ とする。

経路と球の台を含む有界な部分領域で $w-u$ を比較する。側面では $w=0$、$u\ge0$ である。$c$ の下界が負でも、$w-u$ に減衰する時間因子を掛けて零次係数を正にすればよい。正の最大値が側面や初期面以外にあれば、そこでは時間微分が非負、空間の勾配は零、Hessian は非正であり、作用素の値は正になる。これは劣解の不等式と矛盾する。従って $u(t_1,x)\ge w(t_1,x)=\delta e^{-K(t_1-s)}R^6>0$ である。この結果を $M-u$ に適用すると、$c\ge0$ のとき非負な上限 $M$ に内部で達した劣解は、それより前の時刻の空間断面で $M$ に等しい。これが古典的な放物型の強最大値原理である。

## 半群理論による定式化

:::definition[強連続半群]{#def-pde-heat-semigroup}
$H$ の有界作用素族 $S(t)$ が $S(0)=I$、$S(t+s)=S(t)S(s)$、$\lim_{t\downarrow0}S(t)g=g$ を $H$ で満たすとき :term[強連続半群]{en="strongly continuous semigroup" yomi="きょうれんぞくはんぐん"} という。
:::

対称強圧的な $A$ について $S(t)g=\sum_j e^{-t\lambda_j}g_j e_j$ と置く。半群性は係数の積から、強連続性は $\ell^2$ の優収束から従う。$\norm{S(t)}\le1$ である。生成作用素は $\lim_{t\downarrow0}(S(t)g-g)/t$ が $H$ に存在するときの極限である。成分ごとの極限は $-\lambda_j g_j$ なので、その存在には $\sum_j\lambda_j^2|g_j|^2<\infty$ が必要である。これは十分でもある。$|(e^{-t\lambda_j}-1)/t|\le\lambda_j$ により $\ell^2$ の優収束を適用できるからである。この定義域上の生成作用素が $-A$ である。

$f\in L^2(0,T;H)$ では

$$
u(t)=S(t)g+\int_0^t S(t-s)f(s)\,\dd s
\label{eq-pde-parabolic-semigroup}
$$

を :term[緩和解]{en="mild solution" yomi="かんわかい"} の表示と呼ぶ。積分は $H$ で収束する。各固有係数は Galerkin の ODE を解く。従って $H$ の基底で比較すると先の弱解と一致する。時間依存の作用素では $S(t-s)$ だけで過程を表すことは通常できず、始点と終点を別々に記す二時刻の発展作用素が必要である。一般の生成定理は [作用素半群](/math/functional-analysis/semigroups/) を参照する。

:::example[二時刻の依存性]
拡散率 $1+t$ で固有値 $4$ のモードを始時刻 $s$ から終時刻 $t$ に運ぶと、倍率は $\exp[-4(t-s)-2(t^2-s^2)]$ である。長さ $1$ でも $s=0$ なら $e^{-6}$、$s=1$ なら $e^{-10}$ となる。
:::

:::example[半群の Laplace 変換]
$s>0$ で $\int_0^\infty e^{-st}S(t)g\,\dd t$ を取ると、第 $j$ 係数は $g_j\int_0^\infty e^{-(s+\lambda_j)t}\,\dd t=g_j/(s+\lambda_j)$ である。半群の縮小性で積分は $H$ で収束し、この係数は $(sI+A)^{-1}g$ のものである。例えば $\lambda_j=4$、$s=2$ なら、そのモードは $1/6$ 倍になる。
:::

### 非対称な自律形式が定める半群

$a_t=a$ が時間に依らず、:eqref[eq-pde-parabolic-coercivity] を満たすとする。$A: D(A)\subset H\to H$ を、$u\in V$ について $a(u,v)=(Au,v)_H$ が全 $v\in V$ で成り立つときに定義する。$A+\beta I$ の形式を $\widetilde a$ とすると強圧性があり、$G=-(A+\beta I)$ は散逸的である。任意の $\lambda>0$、$h\in H$ に対し

$$
\widetilde a(u,v)+\lambda(u,v)_H=(h,v)_H\qquad(v\in V)
$$

は Lax–Milgram で可解である。これにより $u\in D(A)$ で $(\lambda I-G)u=h$ が成り立つ。定義域の稠密性も同じ式から従う。実際、$D(A)$ と直交する $h$ にこの式を解き、$v=u$ を入れると右辺は零なので $u=0$ となる。そのとき右辺は全 $v\in V$ で零であり、$h=0$ である。従って :ref[thm-fa-lumer-phillips] により $G$ は縮小半群を生成する。元の $-A$ の半群 $S(t)$ はこれに $e^{\beta t}$ を掛けたものである。

まず $g\in D(A)$ の斉次解を考える。半群性から $S(t)g$ は $D(A)$ に入り、時間微分は $-AS(t)g$ である。これは定義域を定める差商に $S(t)$ を作用させれば分かる。グラフノルムから $V$ ノルムは強圧性で抑えられるので、この解は弱解の空間にも属する。一意性により Galerkin 解と一致する。次に $g\in H$ は $D(A)$ の元で近似する。半群の作用素ノルムの上界と弱解のエネルギー評価によって、両方を $C([0,T];H)$ で極限に移せる。従って任意の $H$ 初期値について $S(t)$ が弱解を与える。

$f\in L^2(0,T;H)$ の場合も、時間に滑らかで $D(A)$ に値を取る外力から始めると、半群の畳み込みは微分できて $u'+Au=f$ を満たす。そのような外力は $L^2(H)$ で稠密である。畳み込みは $\norm{S(t)}\le e^{\beta t}$ により $C([0,T];H)$ で収束し、弱解もエネルギー評価により同じ空間で収束する。従って :eqref[eq-pde-parabolic-semigroup] は非対称な自律形式にも使える。

### 生成作用素の定義域を係数から読む

$S(t)e_j=e^{-t\lambda_j}e_j$ とする。$g$ が $\sum \lambda_j^2|g_j|^2<\infty$ を満たせば、$\left|(e^{-t\lambda_j}-1)/t\right|\le\lambda_j$ なので優収束で $(S(t)g-g)/t\to-\sum\lambda_jg_je_j$ が $H$ で成立する。逆に $H$ で差商が収束すれば各係数の極限は $-\lambda_jg_j$ なので、その二乗和は有限である。従って半群の生成作用素は定義域 $D(A)$ を持つ $-A$ である。

$h>0$ について $R_h=(I+hA)^{-1}$ とすると、$R_{t/N}^N g$ の係数は $(1+t\lambda_j/N)^{-N}g_j$ である。各係数は $e^{-t\lambda_j}g_j$ へ収束し、その絶対値は $|g_j|$ 以下である。従って $H$ で $R_{t/N}^N g\to S(t)g$ である。時間刻みを解の存在証明に用いる際、この計算は自律対称の場合に陰的離散化が同じ半群を選ぶことを示す。

:::example[陰的刻みのレゾルベント]
二基底の系で陰的 Euler 法を半群表示と比べる。刻みは $1/2$ とする。

Galerkin の剛性行列を $K$ とする。対角成分は座標方向での値を与える。固有方程式から直交する方向を取り出すため、特性多項式 $\det(\lambda I-K)$ を計算する。

$$
\begin{aligned}
K&=\left[\begin{matrix}\frac{114}{25} & - \frac{48}{25}\\- \frac{48}{25} & \frac{86}{25}\end{matrix}\right]\\
\det(\lambda I-K)&=\lambda^{2} - 8 \lambda + 12=\left(\lambda - 6\right) \left(\lambda - 2\right)\\
d_1&=2,\quad v_1=\left[\begin{matrix}\frac{3}{5}\\\frac{4}{5}\end{matrix}\right]\\
d_2&=6,\quad v_2=\left[\begin{matrix}- \frac{4}{5}\\\frac{3}{5}\end{matrix}\right]
\end{aligned}
$$

$V$ の列に固有ベクトルを並べる。$V^\top V=I$ と $V^\top KV=D$ を使えば、元の座標 $X$ と固有座標 $Y=V^\top X$ の双方で初期係数を記述できる。二次形式は固有座標にすると交差項を持たない。

$$
\begin{aligned}
V&=\left[\begin{matrix}\frac{3}{5} & - \frac{4}{5}\\\frac{4}{5} & \frac{3}{5}\end{matrix}\right]\\
V^\top KV&=\left[\begin{matrix}2 & 0\\0 & 6\end{matrix}\right]\\
X^\top KX&=\frac{2 \left(57 x_{1}^{2} - 48 x_{1} x_{2} + 43 x_{2}^{2}\right)}{25}\\
\sum_j d_jY_j^2&=2 Y_{1}^{2} + 6 Y_{2}^{2}
\end{aligned}
$$

一つの初期係数で固有係数を計算し、モードごとの散逸を求める。直交変換なので二乗ノルムは変わらない。二次形式を二乗ノルムで割ると、固有値を非負の重みで平均した値になる。

$$
\begin{aligned}
X&=\left[\begin{matrix}3\\0\end{matrix}\right]\\
Y&=\left[\begin{matrix}\frac{9}{5}\\- \frac{12}{5}\end{matrix}\right]\\
|X|^2&=9\\
X^\top KX&=\frac{1026}{25}\\
\frac{X^\top KX}{|X|^2}&=\frac{114}{25}\\
\min_j d_j&=2\\
\max_j d_j&=6
\end{aligned}
$$

一歩で固有係数はそれぞれ $1/(1+d_1/2)=1/2$ と $1/(1+d_2/2)=1/4$ 倍になる。半時刻の厳密な乗数は $e^{-1}$ と $e^{-3}$。両方とも縮小するが一歩の値は違う。刻みを小さくし総時刻を固定すると、陰的刻みの積が各指数へ収束する。
:::

### 陰的 Euler 法での一様な時間評価

等間隔の刻み $h$、時刻 $t_j=jh$ に対し、$V^*$ で
$$
(u^j-u^{j-1})/h+A_j u^j=f_j
$$
を解く。$A_j$ の形式が一様強圧、$f_j$ は区間平均とし、$u^0\in H$ を許す。$A_j+h^{-1}I:V\to V^*$ を Lax–Milgram で解ける。$u^j$ で試験し、
$$
2(u^j-u^{j-1},u^j)_H=\norm{u^j}_H^2-\norm{u^{j-1}}_H^2+\norm{u^j-u^{j-1}}_H^2
$$
を使う。Young の不等式で外力を吸収し、$j$ について和を取れば、$\max_j\norm{u^j}_H^2$、$\sum_j\norm{u^j-u^{j-1}}_H^2$、$h\sum_j\norm{u^j}_V^2$ は $C(\norm{u^0}_H^2+\int\norm f_{V^*}^2)$ 以下である。時間の差商 $D_hu^j=(u^j-u^{j-1})/h$ は離散方程式から $h\sum\norm{D_hu^j}_{V^*}^2$ で一様有界にできる。

階段状の補間を $\bar u_h(t)=u^j$、区分線形の補間を $\widehat u_h$ とする。両者の $L^2(0,T;H)$ での差の二乗は $h\sum\norm{u^j-u^{j-1}}_H^2/3$ で、零へ行く。これは各小区間で差が $((t_j-t)/h)(u^{j-1}-u^j)$ であることから直接分かる。階段補間の $L^2(V)$ 弱極限と、線形補間の時間微分の $L^2(V^*)$ 弱極限は、時間の部分積分で同じ極限関数とその微分になる。$A_j=A$ が時間に依らない場合は、$A\bar u_h$ の極限を弱式でそのまま取れる。時刻零で試験関数を消さなければ、$\widehat u_h(0)=u^0$ が初期項に入る。従って弱解が得られ、一意性により Galerkin 法と同じ解を選ぶ。時間依存の場合は $A_j$ の取り方も必要であり、例えば作用素ノルムで時間に連続な形式なら最大幅 $h$ での一様近似から同じ極限操作ができる。

:::example[一歩の離散エネルギーの不足分]
スカラー $D_hu^j+2u^j=0$、$h=1/2$、$u^0=4$ とする。最初の値は $u^1=2$ である。二乗ノルムの変化は $4-16=-12$、差の二乗は $(2-4)^2=4$ であるから、積の式の左辺 $2(-2)2=-8$ と一致する。これに拡散項 $2h\cdot2\cdot(u^1)^2=8$ を足すと零になる。差の二乗を落とすと、連続のエネルギー式と違って等号を維持できない。
:::

:::example[自律な半群の外力畳み込み]
$(0,\pi)$ で $\alpha=1$、初期値は $2\sin3x$、外力は $3e^{-9t}\sin3x$ とする。減衰率と外力の指数が一致する場合を調べる。

まず $u(t,x)=a(t)\sin(mx)$ を代入する。空間二階の微分は $-m^2$ 倍であり、時間微分は $a^{\prime}$ だけに作用する。従って時間依存係数を含む場合でも、他の空間モードを結合させなければ一次の方程式に帰着する。ここで $\mu(t)$ を空間作用素のこのモードでの係数とする。

$$
\begin{aligned}
m&=3\\
\mu(t)&=9\\
f_m&=3 e^{- 9 t}\\
M(t) &= \int_0^t \mu(s)\,\mathrm{d}s=\int\limits_{0}^{t} 9\, dz=9 t\\
I(t)&=e^{M(t)}=e^{9 t}
\end{aligned}
$$

積分因子を掛けると $(Ia)^{\prime}=If_m$ である。零時刻から積分し、初期係数を加えてから因子で割る。外力を単に時間積分するのと違い、積分の中で指数の重みが掛かる。

$$
\begin{aligned}
If_m&=3\\
\int_0^t I(s)f_m(s)\,\mathrm{d}s&=\int\limits_{0}^{t} 3\, dz=3 t\\
a(t)&=\left(3 t + 2\right) e^{- 9 t}\\
a^{\prime}&=\frac{d}{d t} \left(3 t + 2\right) e^{- 9 t}=- 3 \left(9 t + 5\right) e^{- 9 t}\\
a^{\prime}+\mu a&=3 e^{- 9 t}\\
a(0)&=2,\qquad t=0
\end{aligned}
$$

正弦の二乗積分は $\pi/2$ である。解の二乗ノルムの半分を時間微分し、散逸および入力から直接求めた変化率と比べる。係数が時間に依存しても、ここで微分する二乗ノルム自体には時間依存の作用素が入らない。

$$
\begin{aligned}
u(t,x)&=\left(3 t + 2\right) e^{- 9 t} \sin{\left(3 x \right)}\\
E(t)&=\frac{\pi \left(3 t + 2\right)^{2} e^{- 18 t}}{4}\\
E^{\prime}&=\frac{d}{d t} \frac{\pi \left(3 t + 2\right)^{2} e^{- 18 t}}{4}=- \frac{3 \pi \left(3 t + 2\right) \left(9 t + 5\right) e^{- 18 t}}{2}\\
\frac{\pi}2 a(f_m-\mu a)&=- \frac{3 \pi \left(3 t + 2\right) \left(9 t + 5\right) e^{- 18 t}}{2}\\
u(t_*,x_*)&=\frac{3}{e^{3}},\qquad t=\frac{1}{3},x=\frac{\pi}{6}\\
E(t_*)&=\frac{9 \pi}{4 e^{6}},\qquad t=\frac{1}{3}
\end{aligned}
$$

半群を経由した外力の寄与は $\int_0^t e^{-9(t-s)}3e^{-9s}\dd s$ で、指数の積が $e^{-9t}$ にそろうため $3t e^{-9t}$ になる。これが積分因子で得た項と同じである。指数が一致しても全体の解は有界であり、波動の共鳴による長時間成長とは異なる。
:::

:::example[生成作用素のレゾルベントが初期値と外力を運ぶ]
一つの放物型モードで $a^{\prime}+3a=e^{-3t}$、$a(0)=1$ とする。変換領域で $(s+3)^{-1}$ の役割を見る。

時間変換の変数を $s$ とし、$\mathcal L[f](s)=\int_0^\infty e^{-st}f(t)\dd t$ とする。大きい $s$ の範囲で積分が収束する。時間微分を部分積分すると、$t=0$ の値が分子の項になる。まず外力の変換、微分から来る分母、初期データから来る分子を分ける。

$$
\begin{aligned}
f(t)&=e^{- 3 t}\\
\widetilde f(s)&=\frac{1}{s + 3}\\
\text{微分の多項式}&=s + 3\\
\text{初期項}&=1\\
\widetilde a(s)&=\frac{s + 4}{s^{2} + 6 s + 9}
\end{aligned}
$$

次に分母を因数分解する。一つの極が重なるときは単純な指数の係数に加え、時間のべきが掛かる項を求める。例えば $(s+\lambda)^{-2}$ は $te^{-\lambda t}$、分母の二次式 $(s+\lambda)^2+\omega^2$ は減衰する正弦または余弦に対応する。

$$
\begin{aligned}
\widetilde a(s)&=\frac{1}{s + 3} + \frac{1}{\left(s + 3\right)^{2}}\\
a(t)&=\left(t + 1\right) e^{- 3 t}\\
a^{\prime}&=\frac{d}{d t} \left(t + 1\right) e^{- 3 t}=- \left(3 t + 2\right) e^{- 3 t}
\end{aligned}
$$

最後に時刻を戻して初期条件を調べ、次に得た時間関数を元の左辺に代入する。変換領域の式を満たしても、逆変換で定数や初期項を落とすとこの二つの確認のどちらかに失敗する。

$$
\begin{aligned}
a(0)&=1,\qquad t=0\\
\text{時間方程式の左辺}&=e^{- 3 t}\\
a(t_*)&=\frac{4}{3 e},\qquad t=\frac{1}{3}\\
a^{\prime}(t_*)&=- \frac{3}{e},\qquad t=\frac{1}{3}
\end{aligned}
$$

初期項は $(s+3)^{-1}$、外力の項は $(s+3)^{-2}$ である。後者はレゾルベントが二度掛かった形で、逆変換の $te^{-3t}$ に対応する。半群経由の外力積分でも同じ時間因子を得る。
:::

## 演習問題

::::exercise{level=basic}
$u_t=u_{xx}$、$(0,\pi)$ の固定境界で初期値 $2\sin2x$ の解を求めよ。

:::solution
固有値は $4$ であり、Galerkin 係数は $a'+4a=0$、$a(0)=2$ を満たす。よって $u=2e^{-4t}\sin2x$ である。
:::
::::

::::exercise{level=basic}
$(0,\pi)$ の零境界で $u_t-2u_{xx}=2\sin x$ とする。初期値 $u(0,x)=\sin x$ に対して $\sin x$ の係数を求めよ。

:::solution
係数方程式は $a'+2a=2$ である。初期係数 $1$ はこの ODE の定常値なので $a=1$ が一意な解である。
:::
::::

::::exercise{level=basic}
$u_t=(1+t)u_{xx}$ の $\sin x$ モードを時刻 $0$ から $2$ に運ぶ倍率を求めよ。

:::solution
指数は $-\int_0^2(1+s)\dd s=-(2+2)=-4$ である。したがって倍率は $e^{-4}$ である。
:::
::::

::::exercise{level=basic}
$c\ge0$ で外力 $f=0$ の弱解の初期値が非負なら、解の下界を求めよ。

:::solution
$-u$ は同じ斉次方程式を満たし初期値が非正である。:ref[prop-pde-parabolic-positive-part] により $-u\le0$ だから $u\ge0$ である。
:::
::::

::::exercise{level=advanced}
対称作用素の半群で $\norm{AS(t)}$ を $t$ だけで上から抑えよ。

:::solution
固有展開で各倍率は $\lambda_j e^{-t\lambda_j}$ である。連続変数 $\lambda>0$ に広げて微分すると $e^{-t\lambda}(1-t\lambda)$ で、最大値は $\lambda=1/t$ の $1/(et)$ である。Parseval から同じものが作用素ノルムの上界である。
:::
::::

::::exercise{level=advanced}
弱解のエネルギー評価に必要な $v=u_N$ と $v=u_N'$ の違いを説明せよ。

:::solution
前者では $(u_N',u_N)_H=\frac12(\norm{u_N}_H^2)'$ を常に使える。後者で $a(u_N,u_N')=\frac12(a(u_N,u_N))'$ とするには $a$ の対称性と時間独立性を使う。時間依存なら右辺の微分に $\frac12\partial_t a_t(u_N,u_N)$ が加わり別の評価が要る。
:::
::::

::::exercise{level=advanced}
$(0,\pi)$ の Dirichlet 熱作用素で $g_j=j^{-2}$ とする。$g\in V$ と $g\in D(A)$ を判定せよ。

:::solution
$V$ の和は $\sum_j j^2j^{-4}=\sum_j j^{-2}<\infty$ である。$D(A)$ の和は $\sum_j j^4j^{-4}=\sum_j1=\infty$ である。よって $g\in V$ だが $g\notin D(A)$ である。
:::
::::

::::exercise{level=advanced}
$\lambda=1$ のモードに $f(t)=e^{-t}$ を加え、初期係数 $0$ とする。Duhamel 積分を計算せよ。

:::solution
$a(t)=\int_0^t e^{-(t-s)}e^{-s}\dd s=t e^{-t}$ である。微分すると $a'+a=e^{-t}$、$a(0)=0$ が確かめられる。
:::
::::

## 確認テスト

::::quiz
$W(0,T)$ の弱解の時間微分は初めどこに属すると仮定するか。

- [x] $L^2(0,T;V^*)$
- [ ] 常に $L^2(0,T;V)$
- [ ] 常に $C([0,T];H)$

:::explanation
二階の空間作用素は一般に $V$ から $V^*$ へ作用するため、時間微分も $V^*$ で測る。解自体はエネルギー恒等式により $H$ 連続な代表を持つ。
:::
::::

::::quiz{answer="4" tol="1e-6"}
$u_t=u_{xx}$ の $\sin2x$ モードの減衰率を求めよ。

:::explanation
負の二階微分が $4$ 倍なので係数は $e^{-4t}$ であり、減衰率は $4$ である。
:::
::::

::::quiz
弱比較で $u$ の上からの符号を抑える試験関数はどれか。

- [x] $u_+$
- [ ] $u_-$ だけ
- [ ] 常に空間定数 $1$

:::explanation
正の部分で値と勾配を取り出し、$\frac12(\norm{u_+}_H^2)'$ と楕円なエネルギーに結び付ける。境界値 $0$ の $u$ では $u_+$ も同じ試験空間に入る。
:::
::::

::::quiz{answer="3" tol="1e-6"}
$a'+3a=1$、$a(0)=0$ の解の定常値の逆数を求めよ。

:::explanation
定常値は $3a=1$ より $a=1/3$ であり、時間解も $(1-e^{-3t})/3$ に収束する。その逆数は $3$ である。
:::
::::

::::quiz{answer="0.25" tol="1e-6"}
$a^{\prime}+6a=0$、$a(0)=1$ に陰的 Euler 法を刻み $1/2$ で一回適用する。$a^1$ を求めよ。

:::explanation
離散式は $(a^1-1)/(1/2)+6a^1=0$ である。整理して $4a^1=1$ を得るので $a^1=1/4$ である。陽的法なら $1-(1/2)6$ となるが、陰的法では右辺の減衰項を新しい値で評価する。
:::
::::


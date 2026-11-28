---
title: 作用素半群
description: 強連続半群の生成作用素をレゾルベントで特徴付け、解析性と発展方程式に応用する。
notation:
  - symbol: 'T(t)'
    meaning: 時刻 t の半群作用素
  - symbol: 'A_\lambda'
    meaning: 生成作用素の Yosida 近似
  - symbol: 'R(\lambda,A)'
    meaning: λI−A の有界逆
---

## この章の目標

- 強連続半群から閉じた生成作用素を求める。
- Hille–Yosida の条件を用いて半群を構成する。
- 散逸性から縮小半群の生成条件を判定する。
- 解析半群の平滑化を調べ、発展方程式の解を表す。

## 強連続半群と生成作用素

時間発展を繰り返せる場合、解を求める写像には時刻の加法に対応する積の法則がある。初期値を動かしたときの連続性に加え、固定した初期値で時刻を動かしたときの連続性も扱う。

:::definition[強連続半群]{#def-fa-semigroup}
Banach 空間 $X$ 上の有界作用素族 $T(t)$、$t\ge0$ が $T(0)=I$、$T(s+t)=T(s)T(t)$ を満たし、各 $x\in X$ で $T(t)x\to x$（$t\downarrow0$）となるとき :term[強連続半群]{en="strongly continuous semigroup" yomi="きょうれんぞくはんぐん"} と呼ぶ。:term[生成作用素]{en="generator" yomi="せいせいさようそ"} $A$ の定義域は $\lim_{h\downarrow0}(T(h)x-x)/h$ が存在する $x$ 全体であり、$Ax$ をその極限とする。
:::

:::lemma[指数上界と領域上の微分]{#lem-fa-generator}
強連続半群には $M\ge1,\omega\in\R$ で $\norm{T(t)}\le Me^{\omega t}$ となるものがある。全ての軌道は連続であり、$A$ は稠密定義で閉じている。$x\in D(A)$ なら $T(t)x\in D(A)$、$\frac{d}{dt}T(t)x=T(t)Ax=AT(t)x$ である。
:::

:::proof
$0<h\le1/n$ の $T(h)$ の族のうち、$\norm{T(h)x}\le n$ が全 $h$ で成立する $x$ の集合を考える。各集合は閉じ、原点での強連続性からその可算和は $X$ を覆う。Baire の定理と球内の二点の差から、ある $n$ について $\sup_{0\le h\le1/n}\norm{T(h)}<\infty$ である。半群則で任意の $t$ を長さ $1/n$ 以下の区間に分ければ指数上界となる。この上界と積の法則から右連続性が従い、左からは $T(t)x-T(t-h)x=T(t-h)(T(h)x-x)$ を使う。

平均 $x_h=h^{-1}\int_0^hT(s)x\,ds\to x$ は差商の変数変換により $Ax_h=(T(h)x-x)/h$ を満たし、定義域の稠密性を示す。$x\in D(A)$ では $T(t)$ と差商が可換なので右微分は $T(t)Ax$ となる。左微分は $T(t-h)(T(h)x-x)/h$ の極限で同じになる。従って $T(t)x-x=\int_0^tT(s)Ax\,ds$ である。$x_n\to x,Ax_n\to y$ のとき、局所有界性でこの積分の極限を取り、零で右微分すれば $Ax=y$ が得られ、閉性が従う。
:::

:::example[列の指数減衰]
$\ell^2$ で $T(t)e_n=e^{-nt}e_n$ は半群でノルムが $e^{-t}$ である。尾部の平方和で支配収束を使うと強連続性が従う。差商は各座標で $-nx_n$ に近づく。$\sum n^2\abs{x_n}^2<\infty$ なら $\abs{(e^{-nh}-1)/h}\le n$ で平方和の極限を取れ、逆に差商が収束すれば座標極限が $\ell^2$ に属さなければならない。
:::

:::example[上三角行列の強連続半群を微分から作る]
$A=\lambda I+N$、$\lambda=-2$、$N=\begin{pmatrix}0 & 4 & 0\\0 & 0 & 1\\0 & 0 & 0\end{pmatrix}$ を考える。$N$ は座標を前へ一つずつ送り、$N^3=0$ である。従って $z\ne\lambda$ で $(zI-A)^{-1}=\sum_{j=0}^{2}(z-\lambda)^{-j-1}N^j$ となる。実際、$zI-A=(z-\lambda)I-N$ を掛けると隣り合う項が消え、末尾の $N^3$ が零なので $I$ が残る。$z=\lambda$ では $e_1$ が核に入るため、スペクトルは一点 $\{\lambda\}$ である。

| 次数 $j$ | $N^j$ | 逆の係数 | $f^{(j)}(\lambda)/j!$ |
| --- | --- | --- | --- |
| 0 | $\begin{pmatrix}1 & 0 & 0\\0 & 1 & 0\\0 & 0 & 1\end{pmatrix}$ | $\frac{1}{3}$ | $-27$ |
| 1 | $\begin{pmatrix}0 & 4 & 0\\0 & 0 & 1\\0 & 0 & 0\end{pmatrix}$ | $\frac{1}{9}$ | $27$ |
| 2 | $\begin{pmatrix}0 & 0 & 4\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix}$ | $\frac{1}{27}$ | $-9$ |

表の第 3 列で足すと $R_A(1)=\begin{pmatrix}\frac{1}{3} & \frac{4}{9} & \frac{4}{27}\\0 & \frac{1}{3} & \frac{1}{9}\\0 & 0 & \frac{1}{3}\end{pmatrix}$ になる。第 4 列は $f(z)=(z-1)^3$ の Taylor 係数であり、同じ $N^j$ と掛けて足すと $f(A)=\begin{pmatrix}-27 & 108 & -36\\0 & -27 & 27\\0 & 0 & -27\end{pmatrix}$ である。多項式の場合は二項定理で、一般の正則関数の場合は逆の各項を Cauchy の導関数公式に代入することで同じ展開を得る。$f(A)$ を作るには $\lambda$ での $2$ 階までの情報があればよい。

$U(\tau)=e^{\tau\lambda}\sum_{j=0}^{2}\tau^jN^j/j!$ とする。$U(0)=I$ で微分すると $U'(\tau)=AU(\tau)$ になる。さらに二つの有限和の積で $N^j$ の係数は $\sum_{k=0}^j\tau^k s^{j-k}/(k!(j-k)!)=(\tau+s)^j/j!$ なので半群則を満たす。$e_3$ での第 3 座標は $e^{\tau\lambda}$ だが、先行する座標にも多項式を含む成分を生じる。例えば時刻 $1$ の作用は $U(1)x=e^{-2}\begin{pmatrix}1 & 4 & 2\\0 & 1 & 1\\0 & 0 & 1\end{pmatrix}x$ であり、$\lambda=-2$ でも有限時間における増幅は上三角成分に依存する。
:::

:::example[モードごとの生成作用素を差商で読む]
正規直交方向 $e_1,\ldots,e_3$ を取り、$Ae_j=-d_je_j$ とする。初期値 $\sum a_j e_j$ と時間に依存しない右辺 $\sum h_j e_j$ に対し、$u'(t)=Au(t)+\sum_j h_j e_j$ を解く。時刻は $t=\frac{1}{2}$ で比較する。$u_j'+d_ju_j=h_j$ に $e^{d_jt}$ を掛けて積分すると $u_j(t)=e^{-d_jt}a_j+(1-e^{-d_jt})h_j/d_j$ になる。

| モード | 減衰率 | 初期値 | 外力 | 解の係数 | 指定時刻の値 |
| --- | --- | --- | --- | --- | --- |
| 1 | $1$ | $1$ | $0$ | $e^{- t}$ | $e^{- \frac{1}{2}}$ |
| 2 | $3$ | $2$ | $0$ | $2 e^{- 3 t}$ | $\frac{2}{e^{\frac{3}{2}}}$ |
| 3 | $5$ | $-1$ | $0$ | $- e^{- 5 t}$ | $- \frac{1}{e^{\frac{5}{2}}}$ |

外力が一定なら定常解は $u_\infty=\sum_j(h_j/d_j)e_j$ である。差 $u_j-u_{\infty,j}=e^{-d_jt}(a_j-h_j/d_j)$ を平方和に入れると
$$
\norm{u(t)-u_\infty}\le e^{-t}\norm{u(0)-u_\infty}
$$
となる。入力によって最初のモードの差が消えていれば、表に出てくる次の減衰率で改善できる。これは全入力に共通のノルム評価と、特定の入力に対して得られる評価の違いである。

斉次項の差商はモードごとに $a_j(e^{-d_jt}-1)/t\to-d_j a_j$ となる。例えば表の初期値で生成作用素を作用させると係数は $(-1,-6,5)$ である。有限モードならどの入力も生成作用素の定義域にある。これを無限個へ広げると差商のノルム収束が成立するために $\sum_j d_j^2\abs{a_j}^2<\infty$ が必要となる。係数ごとの収束だけではこの和の収束は保証されない。
:::

:::proposition[平均と積分が定義域に入る]{#prop-fa-semigroup-integrated-domain}
強連続半群 $T$ の生成作用素を $A$ とする。任意の $x\in X$ と $t>0$ で $v(t)=\int_0^t T(s)x\,ds$ は $D(A)$ に属し、$Av(t)=T(t)x-x$ である。
:::

:::proof
$h>0$ に対し $T(h)v(t)=\int_h^{t+h}T(s)x\,ds$ なので、差は $\int_t^{t+h}T(s)x\,ds-\int_0^hT(s)x\,ds$ である。これを $h$ で割り、各端の連続性を使うと差商は $T(t)x-x$ へ収束する。$v(t)/t\to x$ から定義域の稠密性も再び得る。例えば $T(t)e_j=e^{-jt}e_j$ では $v(t)_j=(1-e^{-jt})x_j/j$ で、元の $x$ が $D(A)$ になくても重み $j$ が積分の分母で相殺される。
:::

:::example[減衰半群が強連続でノルム連続でない例]
$\ell^2$ で $T(t)e_j=e^{-j^2t}e_j$ とする。$t>0$ で $0\le1-e^{-j^2t}<1$、$j\to\infty$ で差が $1$ に近づくので、$\norm{T(t)-I}=1$ である。一方、任意の固定した $x$ では差の平方は $\sum_j\abs{e^{-j^2t}-1}^2\abs{x_j}^2$ であり、可算和の支配収束で零に近づく。例えば $x_j=j^{-2}$ では強連続だが生成作用素の定義域には入らない。実際、差商の各成分は $-j^2x_j=-1$ に近づくので、全体のノルム収束は不可能である。

しかし時刻 $t>0$ を先に取った後では、$\sum_jj^4e^{-2j^2t}\abs{x_j}^2$ が有限であり、$T(t)x\in D(A)$ となる。指数の減衰によって $r e^{-tr}$ が $r\ge0$ 上で有界であるためである。時刻零で差商が収束しない入力でも、正時刻で微分できる場合がある。
:::

:::example[固有方向が直交しない場合の指数上界]
$P=\begin{pmatrix}1 & 4 & 1\\0 & 1 & 2\\0 & 0 & 1\end{pmatrix}$ は可逆であり、$\Lambda=\diag(-1,-2,-4)$ を取って $A=P\Lambda P^{-1}=\begin{pmatrix}-1 & -4 & 5\\0 & -2 & -4\\0 & 0 & -4\end{pmatrix}$ とする。行列 $P$ の列はそれぞれの固有値の固有ベクトルである。入力 $x=\begin{pmatrix}2\\-1\\3\end{pmatrix}$ をこれらに分解するには $Pc=x$ を解いて $c=P^{-1}x$ を得る。標準内積で各列が正規直交しなければ、各列との内積だけで $c$ は求まらない。

| $j$ | 固有値 | 固有ベクトル | 展開係数 | 入力のその成分 |
| --- | --- | --- | --- | --- |
| 1 | $-1$ | $\begin{pmatrix}1\\0\\0\end{pmatrix}$ | $27$ | $\begin{pmatrix}27\\0\\0\end{pmatrix}$ |
| 2 | $-2$ | $\begin{pmatrix}4\\1\\0\end{pmatrix}$ | $-7$ | $\begin{pmatrix}-28\\-7\\0\end{pmatrix}$ |
| 3 | $-4$ | $\begin{pmatrix}1\\2\\1\end{pmatrix}$ | $3$ | $\begin{pmatrix}3\\6\\3\end{pmatrix}$ |

この表で入力の各成分を足せば $x$ に戻る。各成分に対応する固有値を掛けて足すと $Ax=\begin{pmatrix}17\\-10\\-12\end{pmatrix}$ である。従って固有座標の操作と標準座標での行列積が一致する。$\Lambda$ が対角であることは作用の仕組みを簡単にするが、$P$ によるノルムの変化は残る。

斉次方程式 $u'=Au$、$u(0)=x$ は新しい座標で $c_j'=\lambda_jc_j$ となる。従って $u(t)=P\diag(e^{t\lambda_j})c$ である。時刻 $1$ の解は $u(1)=\begin{pmatrix}- \frac{28}{e^{2}} + \frac{3}{e^{4}} + \frac{27}{e}\\- \frac{7}{e^{2}} + \frac{6}{e^{4}}\\\frac{3}{e^{4}}\end{pmatrix}$ になる。差商を各指数について微分すると $u'(0)=Ax=\begin{pmatrix}17\\-10\\-12\end{pmatrix}$ が得られる。作用素全体については $T(t)=P\diag(e^{t\lambda_j})P^{-1}$ であり、同じ $P^{-1}P=I$ から半群則を検算できる。

$\ell^1$ の作用素ノルムを用いれば $\norm{T(t)}\le \norm P_1 \norm{P^{-1}}_1 e^{t\max_j\lambda_j}=50e^{-t}$ である。逆に最大の固有値の固有ベクトルでは長さの比がちょうど $e^{t\max_j\lambda_j}$ なので、指数の値はこれ以上下げられない。対角表示の出口で再び $P$ を掛けることから、指数前の係数も評価に現れる。レゾルベント $R_A(2)=\begin{pmatrix}\frac{1}{3} & - \frac{1}{3} & \frac{1}{2}\\0 & \frac{1}{4} & - \frac{1}{6}\\0 & 0 & \frac{1}{6}\end{pmatrix}$ は同じ問題に Laplace 変換を適用した結果である。
:::

:::proposition[生成作用素が半群を一意に定める]{#prop-fa-semigroup-uniqueness}
二つの強連続半群 $T,S$ が同じ生成作用素 $A$ を持てば、各 $t\ge0$ で $T(t)=S(t)$ である。
:::

:::proof
$x\in D(A)$ とする。$s\mapsto T(t-s)S(s)x$ は $[0,t]$ 上で連続である。$0<s<t$ では積の差商を二項に分け、$S(s)x\in D(A)$ と $AS(s)x=S(s)Ax$ を使うと導関数は $-T(t-s)AS(s)x+T(t-s)S(s)Ax=0$ である。従って両端の値 $T(t)x$ と $S(t)x$ は一致する。$D(A)$ の稠密性と両作用素の有界性により全ての $x$ で一致する。
:::

## Hille–Yosida の定理

:::lemma[半群の Laplace 変換]{#lem-fa-semigroup-resolvent}
$\norm{T(t)}\le Me^{\omega t}$ なら $\lambda>\omega$ で
$$
R(\lambda,A)x=\int_0^\infty e^{-\lambda t}T(t)x\,dt,\qquad
\norm{R(\lambda,A)^n}\le\frac{M}{(\lambda-\omega)^n}
$$
である。
:::

:::proof
積分 $Ry$ は指数上界によりノルムで収束する。$T(h)Ry=e^{\lambda h}\int_h^\infty e^{-\lambda s}T(s)y\,ds$ を差商にすると $ARy=\lambda Ry-y$ である。$x\in D(A)$ では $e^{-\lambda t}T(t)x$ の導関数を積分し、無限遠で零へ向かうことから $R(\lambda I-A)x=x$ を得る。$n$ 回積分の変数の和を $t$ とすると、半群則と Fubini の定理により
$$
R(\lambda,A)^nx=\frac1{(n-1)!}\int_0^\infty t^{n-1}e^{-\lambda t}T(t)x\,dt.
$$
重みの式は $n=1$ から始め、畳み込みで $\int_0^t s^{n-1}/(n-1)!\,ds=t^n/n!$ として帰納的に得られる。指数上界を積分して主張のノルムを得る。
:::

:::theorem[Hille–Yosida の定理]{#thm-fa-hille-yosida}
稠密定義閉作用素 $A$ が $\norm{T(t)}\le Me^{\omega t}$ の強連続半群を生成するための必要十分条件は、全ての実 $\lambda>\omega$ がレゾルベント集合に属し、全ての $n\ge1$ で $\norm{R(\lambda,A)^n}\le M/(\lambda-\omega)^n$ となることである。
:::

:::proof
必要性は直前の補題である。十分性では $A-\omega I$ を扱って $\omega=0$ としてよい。$R_\lambda=R(\lambda,A)$、$A_\lambda=\lambda^2R_\lambda-\lambda I$ とおく。$x\in D(A)$ では $\lambda R_\lambda x-x=R_\lambda Ax$ だから $\lambda R_\lambda x\to x$ である。作用素が一様有界で $D(A)$ が稠密なので全 $x$ でも収束する。$A_\lambda x=\lambda R_\lambda Ax\to Ax$ が $D(A)$ で成立する。指数級数と冪の上界で
$$
\norm{e^{tA_\lambda}}
\le e^{-\lambda t}\left(1+\sum_{n\ge1}\frac{(\lambda t)^nM}{n!}\right)\le M.
$$

レゾルベント恒等式で全ての $A_\lambda$ は可換である。有界作用素の指数級数を微分し、$e^{(t-s)A_\lambda}e^{sA_\mu}x$ を $s$ で積分すると
$$
\norm{e^{tA_\lambda}x-e^{tA_\mu}x}\le M^2t\norm{(A_\lambda-A_\mu)x}.
$$
従って $x\in D(A)$ でコンパクトな時刻区間上一様 Cauchy であり、一様上界 $M$ と稠密性で全 $x$ でも同じとなる。極限 $T(t)x$ は連続で、積を固定ベクトルの強極限と一様上界で移すと半群則を満たす。$x\in D(A)$ で
$$
e^{tA_\lambda}x-x=\int_0^t e^{sA_\lambda}A_\lambda x\,ds
$$
を極限に送れば $T(t)x-x=\int_0^tT(s)Ax\,ds$ となる。よって生成作用素 $G$ は $A$ を拡大する。任意の $\lambda>0$ で $\lambda I-A$ が全射、直前の補題で $\lambda I-G$ が単射なので、$y\in D(G)$ の $(\lambda I-G)y$ と同じ像を持つ $x\in D(A)$ は $y=x$ を満たす。よって $G=A$ である。$\omega$ を戻すには半群に $e^{\omega t}$ を掛ける。
:::

:::example[中央の近似の係数]
$Ae_n=-ne_n$ なら $R(\lambda,A)e_n=(\lambda+n)^{-1}e_n$ であり、その $k$ 乗のノルムは $(\lambda+1)^{-k}$ である。$A_\lambda e_n=-\lambda n/(\lambda+n)e_n$ なので、例えば $\lambda=6,n=3$ での近似固有値は $-2$ で、$\lambda\to\infty$ では $-3$ に近づく。
:::

:::example[Laplace 変換と Yosida 近似を同じモードで計算する]
正規直交方向 $e_1,\ldots,e_3$ を取り、$Ae_j=-d_je_j$ とする。初期値 $\sum a_j e_j$ と時間に依存しない右辺 $\sum h_j e_j$ に対し、$u'(t)=Au(t)+\sum_j h_j e_j$ を解く。時刻は $t=\frac{1}{3}$ で比較する。$u_j'+d_ju_j=h_j$ に $e^{d_jt}$ を掛けて積分すると $u_j(t)=e^{-d_jt}a_j+(1-e^{-d_jt})h_j/d_j$ になる。

| モード | 減衰率 | 初期値 | 外力 | 解の係数 | 指定時刻の値 |
| --- | --- | --- | --- | --- | --- |
| 1 | $2$ | $1$ | $0$ | $e^{- 2 t}$ | $e^{- \frac{2}{3}}$ |
| 2 | $5$ | $0$ | $0$ | $0$ | $0$ |
| 3 | $7$ | $2$ | $0$ | $2 e^{- 7 t}$ | $\frac{2}{e^{\frac{7}{3}}}$ |

外力が一定なら定常解は $u_\infty=\sum_j(h_j/d_j)e_j$ である。差 $u_j-u_{\infty,j}=e^{-d_jt}(a_j-h_j/d_j)$ を平方和に入れると
$$
\norm{u(t)-u_\infty}\le e^{-2t}\norm{u(0)-u_\infty}
$$
となる。入力によって最初のモードの差が消えていれば、表に出てくる次の減衰率で改善できる。これは全入力に共通のノルム評価と、特定の入力に対して得られる評価の違いである。

$\lambda>0$ で斉次半群の Laplace 変換を取ると各係数は $\int_0^\infty e^{-(\lambda+d_j)t}\,dt=1/(\lambda+d_j)$ となる。作用素ノルムは $1/(\lambda+2)$ で、最小の減衰率の方向で達成する。同じ方向が冪の最大も与えるので、第 $k$ 冪は $(\lambda+2)^{-k}$ になる。$\lambda=2$ ではレゾルベントの対角係数が $(\frac{1}{4},\frac{1}{7},\frac{1}{9})$ である。これを $2^2$ 倍して $2I$ を引くと Yosida 近似の係数 $(-1,- \frac{10}{7},- \frac{14}{9})$ が出る。
:::

:::proposition[Yosida 近似と生成作用素の誤差]{#prop-fa-yosida-error}
$A$ が縮小半群を生成し、$A_\lambda=\lambda^2R(\lambda,A)-\lambda I$ とする。$x\in D(A^2)$ なら
$$
\norm{(A_\lambda-A)x}\le\lambda^{-1}\norm{A^2x},\qquad
\norm{e^{tA_\lambda}x-T(t)x}\le(t/\lambda)\norm{A^2x}.
$$
:::

:::proof
$A_\lambda x=\lambda R(\lambda,A)Ax$ と $\lambda R(\lambda,A)Ax-Ax=R(\lambda,A)A^2x$ を使えば第一式である。$s\mapsto e^{(t-s)A_\lambda}T(s)x$ を微分すると $e^{(t-s)A_\lambda}T(s)(A-A_\lambda)x$ になる。ここで $R(\lambda,A)$ と $T(s)$ は半群の Laplace 表示から可換である。両側の作用のノルムは $1$ 以下なので、積分して $t\norm{(A_\lambda-A)x}$ を得る。例えば $A e_j=-j e_j$、$x=e_3$ なら近似生成値の差は $9/(\lambda+3)$ であり、上界は $9/\lambda$ である。
:::

:::example[冪の一様な評価が一次の評価を改良する]
$T(t)=e^{-t}(I+tN)$、$N^2=0$、$N=\left(\begin{smallmatrix}0&2\\0&0\end{smallmatrix}\right)$ とする。最大行和で $\norm{T(t)}=e^{-t}(1+2t)$ であり、その最大値は $t=1/2$ で $2e^{-1/2}$ である。例えば $M=2,\omega=0$ の上界が使える。生成作用素は $A=-I+N$ なので $R(\lambda,A)=(\lambda+1)^{-1}I+(\lambda+1)^{-2}N$ である。冪零性から
$$
R(\lambda,A)^n=(\lambda+1)^{-n}I+n(\lambda+1)^{-n-1}N
$$
であり、ノルムは $(\lambda+1)^{-n}(1+2n/(\lambda+1))$ である。単純に一回分のノルム上界を $n$ 回掛けると $2^n/\lambda^n$ になるが、半群の Laplace 表示で得る定数は $2/\lambda^n$ である。

例えば $\lambda=1$ では冪のノルムは $2^{-n}(1+n)$ になる。$n=2$ で $3/4$、$n=3$ で $1/2$ である。共通の積の上界に必要なのは、次数 $n$ とともに増えない先頭の定数である。
:::

:::proposition[半群の Euler 近似]{#prop-fa-semigroup-euler}
強連続半群 $T$ の生成作用素を $A$ とする。$t>0$ で $(I-tA/n)^{-n}x\to T(t)x$ が全ての $x$ に対して成り立つ。逆作用素の存在は十分大きい $n$ に対して保証される。縮小半群なら $x\in D(A)$ で誤差は $t\norm{Ax}/\sqrt n$ 以下である。
:::

:::proof
まず $\norm{T(t)}\le M$ がすべての $t\ge0$ で成り立つ場合を考える。Laplace の冪の公式に $\lambda=n/t$ を入れると、Euler 近似は $T(s)x$ を密度
$$
p_n(s)=\frac{(n/t)^n}{(n-1)!}s^{n-1}e^{-ns/t}\qquad(s\ge0)
$$
で平均したものになる。部分積分で $\int_0^\infty s^k e^{-\lambda s}ds=k!/\lambda^{k+1}$ を得る。従ってこの密度の積分は $1$、平均は $t$、分散は $t^2/n$ である。$x\in D(A)$ なら軌道の微分式から $\norm{T(s)x-T(t)x}\le M\abs{s-t}\norm{Ax}$ であり、Cauchy–Schwarz で誤差を $Mt\norm{Ax}/\sqrt n$ 以下に抑える。平均作用素のノルムは $M$ 以下なので稠密性により全入力で収束する。

一般に $\norm{T(t)}\le Me^{\omega t}$ の場合、$B=A-\omega I$ の半群 $S(s)=e^{-\omega s}T(s)$ に前の結果を使う。$n>\omega t$ に対し $I-tA/n=(1-\omega t/n)(I-t_nB/n)$、$t_n=t/(1-\omega t/n)$ と書ける。$t_n\to t$ であり、有界半群に対する上の誤差は $t_n$ を入れても零へ向かう。同じ稠密性と一様上界でこの結論は全入力に移る。外の係数は $(1-\omega t/n)^{-n}\to e^{\omega t}$ であり、$S(t_n)x\to S(t)x$ と合わせて元の公式を得る。
:::

:::example[Euler の係数を指数関数と比較する]
$Ae_j=-je_j$、$t=1$、$x=e_2+e_3$ とする。Euler の第 $n$ 近似は $(1+2/n)^{-n}e_2+(1+3/n)^{-n}e_3$ になる。$n=2$ なら係数は $1/4,4/25$、$n=4$ なら $16/81,256/2401$ である。真の値は $e^{-2},e^{-3}$ であり、$\log(1+x)<x$ から有限 $n$ では各近似係数が真値より大きい。$n\log(1+j/n)\to j$ なので各座標で真値へ近づく。

$x\in D(A)$ に対する Euler 近似の評価では、時間幅に基づく係数 $Mt/\sqrt n$ が全てのモードに共通である。固定した $j$ だけを見ると $n\log(1+j/n)=j-j^2/(2n)+O(n^{-2})$ なので係数の誤差は $e^{-j}j^2/(2n)+O(n^{-2})$ となる。固定モードで得られる改善は、$j$ に依存する Taylor 評価による。
:::

:::proposition[レゾルベントの強収束から半群の強収束へ]{#prop-fa-semigroup-resolvent-convergence}
$A_n,A$ がいずれも縮小半群 $T_n,T$ を生成するとする。ある $\lambda>0$ について $R_n=(\lambda I-A_n)^{-1}$ が $R=(\lambda I-A)^{-1}$ に強収束するなら、各 $x$ で $T_n(t)x\to T(t)x$ は有界な時刻区間上一様である。
:::

:::proof
まず $x\in D(A^2)$ を取る。$y(s)=(\lambda I-A)T(s)x$ は連続微分可能であり、$y'(s)=A T(s)(\lambda I-A)x$ である。$[0,b]$ 上の $y,y'$ の像はコンパクトなので、$R_n-R$ の強収束はその像上で一様である。これには、共通のノルム上界 $2/\lambda$ と像の有限網を用いればよい。$w_n(s)=R_ny(s)\in D(A_n)$ は $T(s)x$ に一様に近づく。$A_nR_n=\lambda R_n-I$ から
$$
w_n'(s)-A_n w_n(s)=(R_n-R)y'(s)-\lambda(R_n-R)y(s)
$$
であり、右辺は $[0,b]$ 上で一様に零へ近づく。ここで $Ry'(s)=A T(s)x$ を使った。任意の $t\le b$ について $s\mapsto T_n(t-s)w_n(s)$ を微分して積分すると、$w_n(t)-T_n(t)w_n(0)$ のノルムは $b\sup_{s\le b}\norm{w_n'(s)-A_n w_n(s)}$ 以下である。従って $T_n(t)x\to T(t)x$ が一様に成り立つ。

$D(A^2)$ は稠密である。実際、$\lambda R(\lambda,A)x\to x$ は $D(A)$ 上の公式 $\lambda R(\lambda,A)x-x=R(\lambda,A)Ax$ と一様なノルム上界から従い、同じ上界により $(\lambda R(\lambda,A))^2x\to x$ である。その像は $D(A^2)$ に入る。従ってこの稠密性と $T_n,T$ の共通縮小性を使い、任意の $x$ に結論を延長できる。
:::

:::example[高いモードだけを変更した近似半群の判定]
$\ell^2$ で $Ae_j=-j e_j$ とする。$A_n$ を最初の $n$ 個で同じ係数、それ以外では係数 $0$ とした有界作用素にする。半群は最初の $n$ 成分だけを指数で減衰させ、尾部は固定する。$\lambda>0$ でのレゾルベント差は尾部で $1/\lambda-1/(\lambda+j)$ なので、任意の固定した $x$ の像の平方ノルムは $\lambda^{-2}\sum_{j>n}\abs{x_j}^2$ 以下で零に近づく。従って縮小半群の収束をレゾルベントから得られる。

ここでは直接の半群計算でも、時刻区間に依存せず $\norm{(T_n(t)-T(t))x}^2\le\sum_{j>n}\abs{x_j}^2$ が従う。例えば $x_j=2^{-j}$ で誤差は $2^{-n}/\sqrt3$ 以下である。一方 $t>0$ では尾部の $1-e^{-jt}$ が $1$ に近づくため、$\norm{T_n(t)-T(t)}=1$ であり、作用素ノルムでの収束は起きない。近似半群の生成作用素が全域で有界でも、その上界は $n$ で増大して極限は非有界になりうる。
:::

## Lumer–Phillips の定理と縮小半群

:::definition[散逸性]{#def-fa-dissipative}
Banach 空間の作用素 $A$ が全ての $\lambda>0,x\in D(A)$ について $\norm{(\lambda I-A)x}\ge\lambda\norm x$ を満たすとき、:term[散逸的]{en="dissipative" yomi="さんいつてき"} であるという。Hilbert 空間ではこれは $\Re\inner{Ax}{x}\le0$ と同値である。実際、平方を展開するとこの条件は十分であり、逆に展開した不等式を $2\lambda$ で割って $\lambda\to\infty$ とすれば必要性も従う。
:::

:::theorem[Lumer–Phillips の定理]{#thm-fa-lumer-phillips}
稠密に定義された作用素 $A$ が縮小半群を生成することと、$A$ が散逸的であること及びある $\lambda_0>0$ で $\Img(\lambda_0I-A)=X$ であることは同値である。
:::

:::proof
必要性は Laplace 変換の評価である。十分性では、全射性と散逸性から $\lambda_0I-A$ は全単射で、その逆は有界である。その交換グラフは閉じているので $A$ は閉じている。全射となる正の $\lambda$ の集合は、Neumann 級数で開いている。$\lambda_n\to\lambda>0$ で $\lambda_n$ が全射点なら、$\abs{\lambda-\lambda_n}\norm{R(\lambda_n,A)}\le\abs{\lambda-\lambda_n}/\lambda_n<1$ となる十分大きな $n$ を使い、$\lambda$ も全射点だと分かる。従って全ての $\lambda>0$ で全射である。逆のノルムは $1/\lambda$ 以下で、その冪は $1/\lambda^n$ 以下である。:ref[thm-fa-hille-yosida] に $M=1,\omega=0$ を適用する。
:::

:::example[片側の平行移動]
$L^2(0,\infty)$ の $T(t)f(s)=f(s+t)$ はノルムを増やさない。試験関数では一様連続性と有界な台で強連続性を示せ、稠密性で一般の場合に延長できる。生成作用素は弱微分 $Af=f'$、定義域は $H^1(0,\infty)$ である。実際、$H^1$ では差を導関数の積分で書き、平行移動の強連続性から差商の収束を得る。逆に差商が収束すれば試験関数で変数をずらし、極限が弱微分であることが分かる。左向きに観測する半群なので、原点の値を零にする条件は付かない。
:::

:::example[散逸性はあるが全射性が失われる例]
$L^2(0,1)$ で $A u=u'$ を両端零の $H^1$ に定める。$2\Re\inner{Au}{u}=\abs{u(1)}^2-\abs{u(0)}^2=0$ なので散逸的である。しかし $(\lambda I-A)u=f$ を左端零で解くと $u(t)=-e^{\lambda t}\int_0^t e^{-\lambda s}f(s)\,ds$ であり、右端も零にするには $\int_0^1e^{-\lambda s}f(s)\,ds=0$ が必要になる。$f=1$ ではこの積分が正であるため、全射ではない。

右端だけを零とすると $u(t)=\int_t^1 e^{-\lambda(s-t)}f(s)\,ds$ が全ての $f\in L^2$ に対する解になる。この式は絶対連続で $u'=\lambda u-f\in L^2$、端点は零である。散逸性も $2\Re\inner{Au}{u}=-\abs{u(0)}^2\le0$ となって残る。従って右端だけを零とした作用素は縮小半群を生成する。定義域を拡大すると可解条件が消える過程を、この二つの積分式で確認できる。
:::

:::example[右シフト半群で境界条件を確認する]
$L^2(0,\infty)$ で $T(t)f(s)=f(s-t)$ を $s>t$ に定め、$s\le t$ では零とする。積分変数を移すと $\norm{T(t)f}_2=\norm f_2$ であり、半群則も台を二度移すことで従う。零拡張と平行移動の強連続性から強連続である。$f\in H^1(0,\infty)$ のうち $f(0)=0$ なら零拡張は弱導関数が $L^2$ のままなので、差商は $-f'$ に収束する。逆に差商が $L^2$ で収束すれば、内部の試験関数から弱導関数が $L^2$ と分かる。さらに $0<s<t$ の部分での差商の平方ノルムは $t^{-2}\int_0^t\abs{f(s)}^2ds$ である。境界値が零でなければ、この値は $t^{-1}$ の次数で発散する。従って定義域は零境界値の $H^1$ である。

例えば $f(s)=e^{-s}$ は $H^1$ に属するが右シフトの生成作用素の定義域には入らない。前に扱った左シフトでは同じ関数に $f'=-e^{-s}$ を作用させられる。移動の方向によって、零を補う端点が生成作用素の条件に入るかどうかが変わる。
:::

:::proposition[縮小半群のエネルギー評価]{#prop-fa-semigroup-energy}
Hilbert 空間上で $\Re\inner{Ax}{x}\le-\kappa\norm x^2$、$x\in D(A)$ とし、$A$ が強連続半群 $T$ を生成するとする。このとき $\norm{T(t)}\le e^{-\kappa t}$ である。
:::

:::proof
$x\in D(A)$ では $(d/dt)\norm{T(t)x}^2=2\Re\inner{AT(t)x}{T(t)x}\le-2\kappa\norm{T(t)x}^2$ である。$e^{2\kappa t}$ を掛けて微分すると非正になるので、値は初期値の平方ノルム以下となる。平方根を取って評価を得る。$D(A)$ は稠密で $T(t)$ は有界なので全入力へ移る。例えば $A=-2I+K$、$K^*=-K$ の場合は二次形式で等号が成り立ち、$\norm{T(t)x}=e^{-2t}\norm x$ である。
:::

:::example[零境界値の熱作用素で散逸性を積分から計算する]
区間 $(0,1)$ で $-u''=f$、$f(t)=t^{2} - t + 1$、$u(0)=u(1)=0$ とする。左端で零の原始関数 $F_1(t)=\int_0^t f$、$F_2(t)=\int_0^t F_1$ を置くと候補は $u(t)=\beta t-F_2(t)$ である。右端条件だけで $\beta=F_2(1)/1$ が決まり、$u(t)=- \frac{t^{4}}{12} + \frac{t^{3}}{6} - \frac{t^{2}}{2} + \frac{5 t}{12}$ となる。二解の差は導関数の差を二回積分して一次関数となり、両端で零なので消える。

| 評価点 | $u$ | $u'$ | $-u''$ |
| --- | --- | --- | --- |
| $0$ | $0$ | $\frac{5}{12}$ | $1$ |
| $\frac{1}{4}$ | $\frac{77}{1024}$ | $\frac{37}{192}$ | $\frac{13}{16}$ |
| $\frac{1}{2}$ | $\frac{19}{192}$ | $0$ | $\frac{3}{4}$ |
| $\frac{3}{4}$ | $\frac{77}{1024}$ | $- \frac{37}{192}$ | $\frac{13}{16}$ |
| $1$ | $0$ | $- \frac{5}{12}$ | $1$ |

積分による検算として、$\int_0^1\abs{u'}^2=\int_0^1f\overline u=\frac{269}{5040}$ である。これは $-u''=f$ と端点での値零を使った等式であり、解の導関数の絶対値の平方を積分しても同じ値を得る。

形式 $q(v,w)=\int v'\conj{w'}$ を考える。零境界値の $w$ で部分積分すると $q(u,w)=\inner{f}{w}$ であるから、この形式に対応する非負の Dirichlet 作用素 $L$ は $Lu=f$ を満たす。熱半群の生成作用素は $A=-L$ であり、$Au=-f$ なので $\Re\inner{Au}{u}=-q(u,u)$ である。平方ノルムを順に計算すると $\norm{u}_2^2=\frac{481}{90720}$、$q(u,u)=\frac{269}{5040}$、$\norm{Au}_2^2=\frac{7}{10}$ である。

従って形式ノルム $\sqrt{\norm u_2^2+q(u,u)}$ は $\frac{\sqrt{372610}}{2520}$、作用素のグラフノルムは $\frac{\sqrt{179158}}{504}$ となる。前者は一次導関数で、後者は二次導関数で測っている。
:::

:::corollary[正自己共役作用素の負号は縮小半群を生成する]{#cor-fa-positive-selfadjoint-semigroup}
$L\ge0$ が自己共役なら $A=-L$ は縮小半群を生成する。下界が $L\ge\kappa I$、$\kappa>0$ であればノルムは $e^{-\kappa t}$ 以下である。
:::

:::proof
$\Re\inner{-Lx}{x}\le0$ が散逸性を与える。$(\lambda I+L)^{-1}$ はスペクトル上の $1/(\lambda+s)$ を掛ける作用素であり、$s/(\lambda+s)$ も有界なのでその像は $D(L)$ に属する。従って $\lambda I+L$ は全射である。Lumer–Phillips から縮小半群が存在する。同じ生成作用素を持つ半群の一意性により、これがスペクトル積分 $e^{-tL}$ に等しいと分かる。実際 $e^{-tL}$ は各モードの積と支配収束で強連続半群であり、差商は $D(L)$ 上で $-L$ に収束する。逆に差商がノルム収束すれば、Fatou の補題から入力のスペクトル測度の二次モーメントが有限であり、入力は $D(L)$ に属する。指数上界は $\sup_{s\ge\kappa}e^{-ts}=e^{-\kappa t}$ から得られる。
:::

:::example[強い散逸性とレゾルベントの上界]
$\R^2$ に Euclid ノルムを入れ、$A=\left(\begin{smallmatrix}-2&-3\\3&-2\end{smallmatrix}\right)$ とする。交差項は打ち消し合うため $\inner{Ax}{x}=-2\norm x^2$ である。$\lambda>0$ で $(\lambda I-A)^*(\lambda I-A)=((\lambda+2)^2+9)I$ となり、逆のノルムは $((\lambda+2)^2+9)^{-1/2}$ である。散逸性から出る $1/\lambda$ の上界より小さく、例えば $\lambda=2$ のとき真の値は $1/5$ である。指数写像は $e^{-2t}$ と角度 $3t$ の回転の積で、ノルムは $e^{-2t}$ である。散逸性が減衰を与え、反自己共役部分が長さを変えずに方向を回すことを行列から直接読める。
:::

:::corollary[共役作用素の散逸性による生成判定]{#cor-fa-dual-dissipativity}
Banach 空間上の稠密定義閉作用素 $A$ と、双対空間上の共役 $A^*$ がともに散逸的なら、$A$ は縮小半群を生成する。ここで $A^*$ の定義域のノルム稠密性は仮定しない。
:::

:::proof
$\lambda>0$ に対し、散逸性から $\lambda I-A$ は単射で、その像上の逆は $\lambda^{-1}$ 以下である。値域の Cauchy 列に逆を適用すれば前像も Cauchy 列となり、$A$ の閉性によって値域は閉じる。真の閉部分空間なら Hahn–Banach の定理により値域を消す非零の $\ell\in X^*$ がある。すると $\ell(Ax)=\lambda\ell(x)$ が全 $x\in D(A)$ で成立し、$A^*\ell=\lambda\ell$ となる。これは $A^*$ の散逸性に反する。従って値域は全空間で、Lumer–Phillips の定理が適用できる。
:::

## 解析半群

:::theorem[解析延長と微分評価]{#thm-fa-analytic-semigroup}
複素 Banach 空間の一様有界な強連続半群 $T$ を考える。以下の正則延長を持つ半群を :term[解析半群]{en="analytic semigroup" yomi="かいせきはんぐん"} という。$T$ が開扇形 $\{z\ne0:\abs{\arg z}<\theta\}$（$0<\theta<\pi/2$）へ作用素ノルムで正則に延長され、より狭い各扇形で有界なら、ある $C$ があって $T(t)X\subset D(A)$、$\norm{AT(t)}\le C/t$ が全ての $t>0$ で成立する。逆にこの評価が成り立てば、ある $\theta>0$ でそのような延長がある。
:::

:::proof
正の $t$ を中心とする半径 $ct$ の円が狭い扇形に入るよう $c>0$ を選ぶ。Cauchy 評価により $\norm{T'(t)}\le C/t$ である。$T(t+h)-T(t)=(T(h)-I)T(t)$ から $T'(t)=AT(t)$ である。

逆向きでは、はじめに $t>0$ で各ベクトルに対する強微分が $AT(t)$ であることを示す。$T(t)x=T(t/2)T(t/2)x$ と書けば内部の $T(t/2)x$ が $D(A)$ に属し、生成作用素の微分公式が使える。同じ分割で $n$ 階の微分が $A^n T(t)$ であることを帰納的に得る。各 $T(t/n)$ が定義域に写し、そこ上で $A$ と可換なので $\norm{A^nT(t)}\le(Cn/t)^n$ である。$n!\ge(n/e)^n$ より Taylor 級数は $\abs{z-t}<t/(eC)$ でノルム収束する。$C=0$ なら $T(t)=I$ なので直ちに延長できる。半径を例えば $t/(2+4eC)$ に縮めると、実軸上での Taylor の積分剰余のノルムは、区間上の $s\ge t/2$ を用いて幾何級数型に零へ向かう。従って級数は元の半群に一致する。

各円板で作った級数は、円板が交わるとき実軸上にも共通区間を持つため、恒等定理で一致する。これらの円板は原点を頂点とするある扇形を覆い、さらに狭い扇形では $t=\Re z$ における級数のノルムを幾何級数で一様に抑えられる。半群則も恒等定理で延長される。実際、まず $w>0$ を固定して $z$ を動かし、次に $z$ を固定して $w$ を動かす。角度を $\pi/2$ 未満に取れば和は同じ扇形に入る。
:::

:::example[指数の平滑化]
$Ae_n=-n e_n$ なら $T(z)e_n=e^{-nz}e_n$ は $\Re z>0$ でノルム正則である。実際、コンパクト部分集合で $\Re z$ の正の下界を取ると、有限階打切りが一様ノルムで収束する。$t>0$ で $\sup_{n\ge1}ne^{-tn}\le1/(et)$ は実変数 $s\mapsto se^{-ts}$ の最大値から従う。これは $T(t)$ が生成作用素の領域へ像を送る定量的評価である。
:::

### 扇形のレゾルベント評価

:::proposition[解析半群のレゾルベント]{#prop-fa-analytic-resolvent}
角度 $\theta$ の正則延長がより狭い扇形で有界とする。$0<\theta'<\theta$ に対し $\abs{\arg\lambda}\le\pi/2+\theta'$、$\lambda\ne0$ で $\lambda$ はレゾルベント集合に属し、$\norm{R(\lambda,A)}\le C_{\theta'}/\abs\lambda$ となる。
:::

:::proof
$\theta'<\theta_0<\theta$ を取り、$\abs{\varphi}\le\theta_0$ で $\Re(\lambda e^{i\varphi})\ge c\abs\lambda$ を満たす角 $\varphi$ を選ぶ。正の実軸に近ければ $\varphi=0$ とし、レゾルベント変数の偏角が大きければ反対向きに $\theta_0$ まで回せるため、右辺の定数を全ての $\lambda$ に共通にできる。半群を射線上で $S_\varphi(t)=T(e^{i\varphi}t)$ と書く。これは有界で強連続であり、生成作用素は $e^{i\varphi}A$ であることを以下に示す。

$x\in D(A)$ では $T'(z)x=T(z)Ax$ が恒等定理から成立する。$t$ と $te^{i\varphi}$ を円弧で結んで積分すると、$x$ に作用させた差のノルムが $Ct\abs\varphi\norm{Ax}$ 以下になる。従って $T(e^{i\varphi}t)x\to x$ であり、扇形での有界性と稠密性から全 $x$ でも強連続である。射線上で微分を積分して $D(A)$ での生成作用素の値が $e^{i\varphi}Ax$ だと分かる。逆の定義域包含は、$x$ が射線半群の生成作用素の定義域に入る場合、同様の議論を射線を基準として実軸へ適用すれば従う。具体的にその射線での微分は生成作用素 $G_\varphi$ と可換であり、恒等定理から $T'(z)x=e^{-i\varphi}T(z)G_\varphi x$ となるので、実軸の差商を取れる。

Laplace 変換の証明は実指数を複素に替えても実部が正なら成立する。$\lambda e^{i\varphi}$ で適用すると
$$
R(\lambda,A)=e^{i\varphi}\int_0^\infty e^{-\lambda e^{i\varphi}t}T(e^{i\varphi}t)\,dt
$$
を得る。ノルムは射線上の上界を $c\abs\lambda$ で割ったもの以下である。
:::

### 扇形のレゾルベント評価から解析性を得る

:::proposition[レゾルベントからの微分評価]{#prop-fa-sectorial-converse}
$T$ を有界な強連続半群とし、その生成作用素を $A$ とする。ある $0<\theta<\pi/2$ について $0\ne\lambda$、$\abs{\arg\lambda}<\pi/2+\theta$ でレゾルベントが存在し、$\norm{R(\lambda,A)}\le M/\abs\lambda$ とする。このとき $T(t)X\subset D(A)$、$\norm{AT(t)}\le C/t$ であり、:ref[thm-fa-analytic-semigroup] により $T$ は解析延長を持つ。
:::

:::proof
角度 $\pi/2<\varphi<\pi/2+\theta$ の二射線を、下側の無限遠から原点を経て上側の無限遠へ向けた経路 $\Gamma$ とする。$\norm{AR(\lambda,A)}=\norm{\lambda R(\lambda,A)-I}\le M+1$ である。よって
$$
B(t)=\frac{1}{2\pi i}\int_\Gamma e^{t\lambda}AR(\lambda,A)\,d\lambda
$$
はノルムで収束し、各射線で $\abs{e^{t\lambda}}=e^{-t\abs\lambda\abs{\cos\varphi}}$ を積分すると $\norm{B(t)}\le(M+1)/(\pi t\abs{\cos\varphi})$ となる。

$y\in D(A^2)$ を固定する。$\gamma,\nu>0$ に対し、$t\ge0$ で $g(t)=e^{-\gamma t}(T(t)Ay-e^{-\nu t}Ay)$、$t<0$ で零とする。これは連続で可積分である。角周波数 $\tau$ による Fourier 変換 $\int e^{-i\tau t}g(t)\,dt$ は、Laplace 変換の公式で、$\lambda=\gamma+i\tau$ として
$$
H_y(\lambda)=R(\lambda,A)Ay-\frac{Ay}{\lambda+\nu}
=\frac{R(\lambda,A)A^2y}{\lambda}
+\frac{\nu Ay}{\lambda(\lambda+\nu)}
$$
である。右端により無限遠で $O(\abs\lambda^{-2})$ なので変換も可積分である。スカラー汎関数を適用して :ref[thm-four-euclid-inversion] を使う。角周波数 $\tau$ では反転積分に $1/(2\pi)$ が付く。汎関数による点の分離に戻すと
$$
T(t)Ay-e^{-\nu t}Ay
=\frac{1}{2\pi i}\int_{\gamma-i\infty}^{\gamma+i\infty}e^{t\lambda}H_y(\lambda)\,d\lambda.
$$
積分路は $\Gamma$ へ変形できる。大円弧の区間では $\Re\lambda\le\gamma$ であり、長さが $O(\abs\lambda)$、被積分関数のノルムが $O(\abs\lambda^{-2})$ なので寄与は消える。原点を避ける小円弧では $R(\lambda,A)Ay=AR(\lambda,A)y$ が有界であるため寄与は消える。間の領域で $H_y$ は正則である。さらに $\Gamma$ を左側の大円弧で閉じれば、$e^{t\lambda}/(\lambda+\nu)$ の積分は $2\pi i\,e^{-\nu t}$ である。極 $-\nu$ が一つ入り、円弧で $\Re\lambda\le R\cos\varphi<0$ なのでその寄与は消える。従って移項すると $AT(t)y=T(t)Ay=B(t)y$ である。

$D(A^2)$ は稠密である。実際、$R(\lambda,A)$ は $D(A)$ の元を $D(A^2)$ に送り、$\lambda R(\lambda,A)\to I$ が強収束するため $\lambda^2 R(\lambda,A)^2x\to x$ となる。任意の $x$ を $D(A^2)$ の元で近似すれば、$T(t)$ と $B(t)$ の有界性及び $A$ の閉性から $T(t)x\in D(A)$、$AT(t)x=B(t)x$ が従う。
:::

:::example[複素時間で減衰率の実部を用いる]
正規直交方向 $e_1,\ldots,e_3$ を取り、$Ae_j=-d_je_j$ とする。初期値 $\sum a_j e_j$ と時間に依存しない右辺 $\sum h_j e_j$ に対し、$u'(t)=Au(t)+\sum_j h_j e_j$ を解く。時刻は $t=\frac{1}{2}$ で比較する。$u_j'+d_ju_j=h_j$ に $e^{d_jt}$ を掛けて積分すると $u_j(t)=e^{-d_jt}a_j+(1-e^{-d_jt})h_j/d_j$ になる。

| モード | 減衰率 | 初期値 | 外力 | 解の係数 | 指定時刻の値 |
| --- | --- | --- | --- | --- | --- |
| 1 | $1$ | $0$ | $0$ | $0$ | $0$ |
| 2 | $4$ | $2$ | $0$ | $2 e^{- 4 t}$ | $\frac{2}{e^{2}}$ |
| 3 | $16$ | $1$ | $0$ | $e^{- 16 t}$ | $e^{-8}$ |

外力が一定なら定常解は $u_\infty=\sum_j(h_j/d_j)e_j$ である。差 $u_j-u_{\infty,j}=e^{-d_jt}(a_j-h_j/d_j)$ を平方和に入れると
$$
\norm{u(t)-u_\infty}\le e^{-t}\norm{u(0)-u_\infty}
$$
となる。入力によって最初のモードの差が消えていれば、表に出てくる次の減衰率で改善できる。これは全入力に共通のノルム評価と、特定の入力に対して得られる評価の違いである。

斉次部分は複素時刻 $z=s+it$ に拡張され、$\abs{e^{-d_jz}}=e^{-d_js}$ である。従って $s>0$ では各モードに減衰があり、$s<0$ で時間を戻すと逆に高いモードほど増える。各有限モードでは複素指数関数の微分を取れる。無限個なら $s\ge\delta>0$ に限って $\sup_{r\ge0}r e^{-r\delta}=1/(e\delta)$ のような上界を使い、差商の一様制御を行う。この最大値は微分 $e^{-r\delta}(1-r\delta)$ の符号が $r=1/\delta$ で変わることから得る。
:::

:::proposition[解析半群の高階平滑化]{#prop-fa-analytic-higher-smoothing}
有界解析半群について、扇形の内部で作用素ノルムが有界とする。任意の整数 $k\ge1$ と $t>0$ で $T(t)X\subset D(A^k)$ かつ $\norm{A^kT(t)}\le C_k t^{-k}$ となる。
:::

:::proof
扇形内の中心 $t$、半径 $ct$ の円を取る。小さな $c>0$ を角度だけで決めれば全ての $t>0$ で取れる。作用素値 Cauchy の公式は各汎関数にスカラーの公式を適用すれば従い、$\norm{T^{(k)}(t)}\le k! M/(ct)^k$ である。実時刻で $T(h)T(t)=T(t+h)$ なので差商から $AT(t)=T'(t)$ を得る。$T(h)T^{(k-1)}(t)=T^{(k-1)}(t+h)$ は半群則を $t$ で微分した式である。

従って帰納的に $AT^{(k-1)}(t)=T^{(k)}(t)$ となる。二つを組み合わせて定義域と評価を得る。例えば $Ae_j=-j^2e_j$ の場合には $\norm{A^kT(t)}\le (k/(et))^k$ が成り立つ。$r^ke^{-tr}$ を $r\ge0$ 上で微分すると $r=k/t$ が最大点である。
:::

:::example[ユニタリ群が有界な解析拡張を持たない場合]
$T(t)e_j=e^{ijt}e_j$、$j\ge1$ とする。強連続性は係数の支配収束による。正の実軸の両側を含む扇形に有界作用素値の正則拡張があると仮定する。成分 $\inner{T(z)e_j}{e_k}$ はスカラー正則で、正の実軸上で $\delta_{jk}e^{ijt}$ と等しい。一致の定理から扇形内で $\delta_{jk}e^{ijz}$ と等しくなる。任意の $\Im z<0$ の点では $\abs{e^{ijz}}=e^{-j\Im z}$ が無限大に向かうので、$T(z)$ は有界作用素になれない。

この例では実時刻で全てのモードの大きさが $1$ のままである。複素時間の半平面のうち $\Im z>0$ だけなら指数は減衰するが、実軸を内部とした扇形では逆側の増大を避けられない。正の時刻での軌道の連続性と、扇形への解析拡張を別々に判定する必要がある。
:::

## 発展方程式への応用

:::theorem[定数変化公式]{#thm-fa-variation-constants}
$A$ が強連続半群 $T$ の生成作用素で、$x\in D(A)$、$f:[0,b]\to X$ が連続微分可能なら
$$
u(t)=T(t)x+\int_0^t T(t-s)f(s)\,ds
$$
は $u'=Au+f$、$u(0)=x$ の一意な古典解である。$x\in X$、$f$ が連続の場合もこの式は連続な関数を与え、それを :term[緩和解]{en="mild solution" yomi="かんわかい"} と呼ぶ。
:::

:::proof
積分項を $v(t)=\int_0^tT(r)f(t-r)\,dr$ と書くと、連続微分可能性と局所有界性により $v'(t)=T(t)f(0)+\int_0^tT(r)f'(t-r)\,dr$ である。一方、積分を時刻 $t+h$ の前後に分けると
$$
v(t+h)=T(h)v(t)+\int_t^{t+h}T(t+h-s)f(s)\,ds.
$$
最後の積分を $h$ で割った極限は $f(t)$ である。従って $v(t)\in D(A)$、$Av(t)=v'(t)-f(t)$ が従う。$t=b$ の主張は $f$ の連続微分可能な延長を用いるか、左から $Av(t)$ の極限に閉性を使えばよい。斉次項については生成作用素の微分公式を使う。一意性について、古典解 $w$ に対し $s\mapsto T(t-s)w(s)$ を $0<s<t$ で微分すると $T(t-s)(w'(s)-Aw(s))$ である。積の差商を分けて局所有界性を使えばこの微分が正当化できる。積分すると表示を得る。
:::

:::example[熱方程式の正弦モード]
$L^2[0,1]$ の Dirichlet Laplace 作用素 $L$ は正自己共役であり、固有関数 $e_n=\sqrt2\sin(n\pi t)$ の固有値は $n^2\pi^2$ である。熱方程式に対する $A=-L$ の半群は係数に $e^{-n^2\pi^2 s}$ を掛ける。例えば初期値 $e_1+2e_2$ なら時刻 $s$ の解は $e^{-\pi^2s}e_1+2e^{-4\pi^2s}e_2$ である。$s>0$ では任意の整数 $k$ について $\sup_{\lambda\ge0}\lambda^ke^{-s\lambda}<\infty$ なので、初期値が一般の $L^2$ でも解は $D(L^k)$ に入る。
:::

:::example[強制熱方程式の三つの正弦成分]
正規直交方向 $e_1,\ldots,e_3$ を取り、$Ae_j=-d_je_j$ とする。初期値 $\sum a_j e_j$ と時間に依存しない右辺 $\sum h_j e_j$ に対し、$u'(t)=Au(t)+\sum_j h_j e_j$ を解く。時刻は $t=\frac{1}{4}$ で比較する。$u_j'+d_ju_j=h_j$ に $e^{d_jt}$ を掛けて積分すると $u_j(t)=e^{-d_jt}a_j+(1-e^{-d_jt})h_j/d_j$ になる。

| モード | 減衰率 | 初期値 | 外力 | 解の係数 | 指定時刻の値 |
| --- | --- | --- | --- | --- | --- |
| 1 | $\pi^{2}$ | $0$ | $2$ | $\frac{2}{\pi^{2}} - \frac{2 e^{- \pi^{2} t}}{\pi^{2}}$ | $- \frac{2}{\pi^{2} e^{\frac{\pi^{2}}{4}}} + \frac{2}{\pi^{2}}$ |
| 2 | $4 \pi^{2}$ | $1$ | $-1$ | $- \frac{1}{4 \pi^{2}} + \frac{e^{- 4 \pi^{2} t}}{4 \pi^{2}} + e^{- 4 \pi^{2} t}$ | $\frac{- e^{\pi^{2}} + 1 + 4 \pi^{2}}{4 \pi^{2} e^{\pi^{2}}}$ |
| 3 | $9 \pi^{2}$ | $2$ | $0$ | $2 e^{- 9 \pi^{2} t}$ | $\frac{2}{e^{\frac{9 \pi^{2}}{4}}}$ |

外力が一定なら定常解は $u_\infty=\sum_j(h_j/d_j)e_j$ である。差 $u_j-u_{\infty,j}=e^{-d_jt}(a_j-h_j/d_j)$ を平方和に入れると
$$
\norm{u(t)-u_\infty}\le e^{-\pi^{2}t}\norm{u(0)-u_\infty}
$$
となる。入力によって最初のモードの差が消えていれば、表に出てくる次の減衰率で改善できる。これは全入力に共通のノルム評価と、特定の入力に対して得られる評価の違いである。

同じ係数を定数変化公式で作ると、斉次項は初期値 $a_j$ に $e^{-d_jt}$ を掛けたもの、積分項は $\int_0^t e^{-d_j(t-s)}h_j\,ds=h_j(1-e^{-d_jt})/d_j$ である。両項は表の微分方程式の解と一致する。初期値と外力にゼロの座標を入れたモードは時間発展でも零を保つ。一般の入力を有限切断した場合には、切り落とした初期値のノルムと外力の積分を別々に評価する。
:::

:::example[各モードの定常値と減衰を分ける]
正規直交方向 $e_1,\ldots,e_4$ を取り、$Ae_j=-d_je_j$ とする。初期値 $\sum a_j e_j$ と時間に依存しない右辺 $\sum h_j e_j$ に対し、$u'(t)=Au(t)+\sum_j h_j e_j$ を解く。時刻は $t=\frac{1}{3}$ で比較する。$u_j'+d_ju_j=h_j$ に $e^{d_jt}$ を掛けて積分すると $u_j(t)=e^{-d_jt}a_j+(1-e^{-d_jt})h_j/d_j$ になる。

| モード | 減衰率 | 初期値 | 外力 | 解の係数 | 指定時刻の値 |
| --- | --- | --- | --- | --- | --- |
| 1 | $2$ | $1$ | $1$ | $\frac{1}{2} + \frac{e^{- 2 t}}{2}$ | $\frac{1}{2 e^{\frac{2}{3}}} + \frac{1}{2}$ |
| 2 | $8$ | $2$ | $0$ | $2 e^{- 8 t}$ | $\frac{2}{e^{\frac{8}{3}}}$ |
| 3 | $18$ | $1$ | $1$ | $\frac{1}{18} + \frac{17 e^{- 18 t}}{18}$ | $\frac{17 + e^{6}}{18 e^{6}}$ |
| 4 | $32$ | $0$ | $2$ | $\frac{1}{16} - \frac{e^{- 32 t}}{16}$ | $\frac{1}{16} - \frac{1}{16 e^{\frac{32}{3}}}$ |

外力が一定なら定常解は $u_\infty=\sum_j(h_j/d_j)e_j$ である。差 $u_j-u_{\infty,j}=e^{-d_jt}(a_j-h_j/d_j)$ を平方和に入れると
$$
\norm{u(t)-u_\infty}\le e^{-2t}\norm{u(0)-u_\infty}
$$
となる。入力によって最初のモードの差が消えていれば、表に出てくる次の減衰率で改善できる。これは全入力に共通のノルム評価と、特定の入力に対して得られる評価の違いである。

同じ係数を定数変化公式で作ると、斉次項は初期値 $a_j$ に $e^{-d_jt}$ を掛けたもの、積分項は $\int_0^t e^{-d_j(t-s)}h_j\,ds=h_j(1-e^{-d_jt})/d_j$ である。両項は表の微分方程式の解と一致する。初期値と外力にゼロの座標を入れたモードは時間発展でも零を保つ。一般の入力を有限切断した場合には、切り落とした初期値のノルムと外力の積分を別々に評価する。
:::

:::proposition[緩和解の初期値と右辺への連続依存性]{#prop-fa-mild-stability}
$\norm{T(t)}\le Me^{\omega t}$ のとき、初期値と連続な右辺 $(x,f)$、$(y,g)$ に対する緩和解の差は
$$
\norm{u(t)-v(t)}
\le Me^{\omega t}\norm{x-y}
+M\int_0^t e^{\omega(t-s)}\norm{f(s)-g(s)}\,ds
$$
である。
:::

:::proof
二つの定数変化公式を引くと、初期差に $T(t)$ を掛けた項と右辺の差の畳み込みに分かれる。それぞれに半群の上界を適用し、積分のノルムはノルムの積分以下であることを使えばよい。例えば $M=1,\omega=-2$ のとき、初期差が $1$、右辺の差が常に $3$ 以下なら、全誤差は $e^{-2t}+(3/2)(1-e^{-2t})$ 以下である。時刻が大きくなるにつれて初期値の影響は消えるが、右辺が異なる限り定常状態の差は残りうる。
:::

:::proposition[解析半群と Hölder 連続な外力]{#prop-fa-analytic-semigroup-forcing}
$T$ を有界解析半群とし、$f:[0,b]\to X$ が指数 $\alpha>0$ の Hölder 連続関数とする。初期値 $x\in X$ に対する定数変化公式の緩和解 $u$ は各 $0<t<b$ で $D(A)$ に属し、$u'(t)=Au(t)+f(t)$ を満たす。
:::

:::proof
積分項 $v(t)$ を
$$
v(t)=\int_0^tT(r)(f(t-r)-f(t))\,dr+\int_0^tT(r)f(t)\,dr
$$
と分ける。第二項は :ref[prop-fa-semigroup-integrated-domain] により定義域に入り、その作用像は $(T(t)-I)f(t)$ である。第一項の積分下端を $\eps>0$ にすると、被積分関数は $D(A)$ に入り、作用像はノルムで連続である。従ってグラフノルムで積分できる。$\norm{AT(r)}\le C/r$ と Hölder の評価により、その作用像のノルムは $C' r^{\alpha-1}$ 以下だから $\eps\downarrow0$ で積分が収束する。$A$ の閉性で元の第一項に $A$ を作用させられる。以上により $Au(t)$ は $t>0$ で連続である。実際、積分区間を $[0,b]$ にそろえて指示関数で切れば、被積分関数の点ごとの連続性と共通の可積分上界 $C' r^{\alpha-1}$ が使える。

定数変化公式を時刻 $t+h$ で分けると $u(t+h)=T(h)u(t)+\int_t^{t+h}T(t+h-s)f(s)\,ds$ となる。$h\downarrow0$ の差商は $Au(t)+f(t)$ へ収束する。左微分では $u(t)-u(t-h)$ に同じ分解を使う。$u(t-h)\in D(A)$ なので $(T(h)-I)u(t-h)=\int_0^h T(r)Au(t-h)\,dr$ と書ける。$Au$ の連続性から左の差商も同じ極限になる。
:::

:::proposition[有界摂動の積分表示と生成作用素]{#prop-fa-semigroup-bounded-perturbation}
$A$ が強連続半群 $T$ を生成し、$\norm{T(t)}\le Me^{\omega t}$ とする。$B\in\mathcal B(X)$ なら $A+B$ は $D(A)$ 上で強連続半群 $S$ を生成し、
$$
S(t)x=T(t)x+\int_0^t T(t-s)BS(s)x\,ds,\qquad
\norm{S(t)}\le Me^{(\omega+M\norm B)t}.
$$
:::

:::proof
$S_0=T$、$S_{j+1}(t)x=\int_0^t T(t-s)BS_j(s)x\,ds$ と繰り返す。帰納的に $\norm{S_j(t)}\le Me^{\omega t}(M\norm B\,t)^j/j!$ となる。実際、前段の $s^j$ を積分すると $t^{j+1}/(j+1)$ が付く。各軌道は積分と局所上界から連続であり、和 $S(t)=\sum_j S_j(t)$ は有界な時刻区間で作用素ノルムについて一様収束する。従って積分式と指数上界を得る。

同じ積分式の連続解は一意である。差が斉次の積分式を満たすなら、反復して時刻区間での差の上界に $C^jt^j/j!$ を掛けた評価を得るためである。$h\ge0$ を固定し、$S(t+h)x$ の積分を区間 $[0,h]$ と $[h,t+h]$ に分けると、前者と初期項は $T(t)S(h)x$、後者は $\int_0^t T(t-r)BS(r+h)x\,dr$ となる。

従って $S(t+h)x$ と $S(t)S(h)x$ は同じ初期値を持つ解であり、一意性から半群則を得る。$S(0)=I$ である。時刻零で積分項を $t$ で割ると $Bx$ に収束する。実際 $BS(s)x\to Bx$ と $T(r)Bx\to Bx$ は $s,r\in[0,t]$ で一様に成立し、$T$ の局所有界性を使える。従って $S$ の差商の極限が存在することと $T$ の差商の極限が存在することは同値で、値の差は $Bx$ である。これが $D(A)$ 上の $A+B$ を生成作用素として同定する。
:::

:::example[有界摂動が異なるモードを混ぜる]
$A=\diag(-1,-3)$、$B=\left(\begin{smallmatrix}0&2\\0&0\end{smallmatrix}\right)$ とする。第 2 成分は $u_2(t)=e^{-3t}x_2$ である。これを第 1 成分の右辺に入れ、積分因子で
$$
u_1(t)=e^{-t}x_1+2\int_0^t e^{-(t-s)}e^{-3s}x_2\,ds
=e^{-t}x_1+(e^{-t}-e^{-3t})x_2
$$
を得る。従って摂動後の半群は $\left(\begin{smallmatrix}e^{-t}&e^{-t}-e^{-3t}\\0&e^{-3t}\end{smallmatrix}\right)$ である。例えば初期値が $e_2$ でも、正時刻では第 1 成分が現れる。積分式の循環は $B$ が第 2 成分から第 1 成分への一方向にしか作用しないため、$S_j=0$ が $j\ge2$ で成立する。指数型の上界より、この例では行和を直接読む方が鋭い。
:::

:::example[非直交モードの時間発展とレゾルベントの照合]
$P=\begin{pmatrix}1 & 2 & 0\\0 & 1 & 3\\0 & 0 & 1\end{pmatrix}$ は可逆であり、$\Lambda=\diag(-2,-4,-6)$ を取って $A=P\Lambda P^{-1}=\begin{pmatrix}-2 & -4 & 12\\0 & -4 & -6\\0 & 0 & -6\end{pmatrix}$ とする。行列 $P$ の列はそれぞれの固有値の固有ベクトルである。入力 $x=\begin{pmatrix}1\\2\\1\end{pmatrix}$ をこれらに分解するには $Pc=x$ を解いて $c=P^{-1}x$ を得る。標準内積で各列が正規直交しなければ、各列との内積だけで $c$ は求まらない。

| $j$ | 固有値 | 固有ベクトル | 展開係数 | 入力のその成分 |
| --- | --- | --- | --- | --- |
| 1 | $-2$ | $\begin{pmatrix}1\\0\\0\end{pmatrix}$ | $3$ | $\begin{pmatrix}3\\0\\0\end{pmatrix}$ |
| 2 | $-4$ | $\begin{pmatrix}2\\1\\0\end{pmatrix}$ | $-1$ | $\begin{pmatrix}-2\\-1\\0\end{pmatrix}$ |
| 3 | $-6$ | $\begin{pmatrix}0\\3\\1\end{pmatrix}$ | $1$ | $\begin{pmatrix}0\\3\\1\end{pmatrix}$ |

この表で入力の各成分を足せば $x$ に戻る。各成分に対応する固有値を掛けて足すと $Ax=\begin{pmatrix}2\\-14\\-6\end{pmatrix}$ である。従って固有座標の操作と標準座標での行列積が一致する。$\Lambda$ が対角であることは作用の仕組みを簡単にするが、$P$ によるノルムの変化は残る。

斉次方程式 $u'=Au$、$u(0)=x$ は新しい座標で $c_j'=\lambda_jc_j$ となる。従って $u(t)=P\diag(e^{t\lambda_j})c$ である。時刻 $1$ の解は $u(1)=\begin{pmatrix}- \frac{2}{e^{4}} + \frac{3}{e^{2}}\\- \frac{1}{e^{4}} + \frac{3}{e^{6}}\\e^{-6}\end{pmatrix}$ になる。差商を各指数について微分すると $u'(0)=Ax=\begin{pmatrix}2\\-14\\-6\end{pmatrix}$ が得られる。作用素全体については $T(t)=P\diag(e^{t\lambda_j})P^{-1}$ であり、同じ $P^{-1}P=I$ から半群則を検算できる。

$\ell^1$ の作用素ノルムを用いれば $\norm{T(t)}\le \norm P_1 \norm{P^{-1}}_1 e^{t\max_j\lambda_j}=40e^{-2t}$ である。逆に最大の固有値の固有ベクトルでは長さの比がちょうど $e^{t\max_j\lambda_j}$ なので、指数の値はこれ以上下げられない。対角表示の出口で再び $P$ を掛けることから、指数前の係数も評価に現れる。レゾルベント $R_A(2)=\begin{pmatrix}\frac{1}{4} & - \frac{1}{6} & \frac{1}{2}\\0 & \frac{1}{6} & - \frac{1}{8}\\0 & 0 & \frac{1}{8}\end{pmatrix}$ は同じ問題に Laplace 変換を適用した結果である。
:::

:::example[積分で時刻依存の外力に共鳴する場合]
実数値の方程式 $u'=-2u+3e^{-2t}$、$u(0)=1$ を考える。定数変化公式は
$$
u(t)=e^{-2t}+3\int_0^t e^{-2(t-s)}e^{-2s}\,ds=e^{-2t}(1+3t)
$$
となる。積分の因子が $s$ に依らないため、時間因子 $t$ が出る。微分すると $u'=e^{-2t}(1-6t)$ であり、右辺 $-2e^{-2t}(1+3t)+3e^{-2t}$ と一致する。最大値は $t=1/6$ で $\frac32e^{-1/3}$ である。初期値からいったん増加しても、後には指数減衰する。外力付きの軌道を調べるとき、斉次半群の縮小性だけでは軌道の単調性を言えない。

外力を $3e^{-\alpha t}$、$\alpha\ne2$ に変えると積分項は $3(e^{-\alpha t}-e^{-2t})/(2-\alpha)$ になる。$\alpha\to2$ ではこの商は $3te^{-2t}$ へ収束し、最初の計算に戻る。異なる指数の公式の分母が零になる場合は、積分を先に実行するか、この極限で取り直せる。
:::

### スペクトル写像と有限時刻で消える半群

:::proposition[半群へのスペクトルの移り方]{#prop-fa-semigroup-spectral-inclusion}
複素 Banach 空間上の強連続半群 $T$ と生成作用素 $A$ に対し、$t>0$ で $e^{t\sigma(A)}\subset\sigma(T(t))$ である。
:::

:::proof
$Kx=\int_0^t e^{\lambda(t-s)}T(s)x\,ds$ とおく。$T(h)Kx$ を変数変換し、積分区間を $[h,t+h]$ に移すと、$Kx\in D(A)$、$(\lambda I-A)Kx=e^{\lambda t}x-T(t)x$ を得る。$x\in D(A)$ では $e^{\lambda(t-s)}T(s)x$ を微分して積分すると $K(\lambda I-A)x=e^{\lambda t}x-T(t)x$ でもある。$e^{\lambda t}I-T(t)$ が可逆なら、その逆と $K$ の積が $\lambda I-A$ の右逆になる。後者は二番目の式により単射なので、この右逆は両側逆であり、有界でもある。従って $\lambda\in\sigma(A)$ なら $e^{\lambda t}$ は $\sigma(T(t))$ に属する。
:::

:::example[生成作用素のスペクトルが空となる平行移動]
複素 $L^2(0,1)$ で、$T(t)f(s)=f(s+t)$ とし、$s+t\ge1$ では零とする。零拡張の平行移動から強連続な縮小半群になる。$A f=f'$ の定義域は $f(1)=0$ を満たす $H^1(0,1)$ である。実際、端点値が零なら区間の右側への零拡張に跳びがなく、導関数の平行移動から差商が収束する。逆に差商が収束すれば、内側の試験関数から弱微分が $L^2$ に属し、$1-t<s<1$ での差商平方ノルム $t^{-2}\int_{1-t}^1|f(s)|^2ds$ の有界性から $f(1)=0$ が必要である。

任意の $\lambda\in\C$ に対し $(\lambda I-A)u=f$ の解は $u(s)=\int_s^1 e^{-\lambda(r-s)}f(r)\,dr$ である。この式は絶対連続で $u'=\lambda u-f$、$u(1)=0$ である。積分区間の長さは高々 $1$ なので $\norm u_2\le e^{|\Re\lambda|}\norm f_2$ を得る。斉次式の解 $ce^{\lambda s}$ は右端条件で零になり、従って $\rho(A)=\C$ である。一方 $T(t)^n=T(nt)=0$ は $nt\ge1$ で成り立つから、$t>0$ で $\sigma(T(t))=\{0\}$ である。この例では生成作用素のスペクトルからは、半群のスペクトルに現れる零を読めない。
:::

## 演習問題

::::exercise{level=basic}
実数上の $T(t)x=e^{-3t}x$ の生成作用素と $R(2,A)$ を求めよ。
:::solution
零での差商は $-3x$ なので $A=-3I$ であり全域定義である。$2I-A=5I$ なので逆は $I/5$ である。
:::
::::

::::exercise{level=basic}
$A=-4I$ の $\lambda=4$ での Yosida 近似を求めよ。
:::solution
$R(4,A)=I/8$ である。従って $A_4=16I/8-4I=-2I$ である。
:::
::::

::::exercise{level=basic}
$\lambda>0$ に対し $\int_0^\infty t e^{-\lambda t}\,dt$ を求めよ。
:::solution
部分積分で境界項は零で、残る積分は $\lambda^{-1}\int_0^\infty e^{-\lambda t}dt=\lambda^{-2}$ である。これはレゾルベントの二乗の Laplace 積分の重みである。
:::
::::

::::exercise{level=basic}
$u'=-2u+1,\ u(0)=0$ を定数変化公式で解け。
:::solution
$T(t)=e^{-2t}$ なので $u(t)=\int_0^t e^{-2(t-s)}\,ds=(1-e^{-2t})/2$ である。微分して代入すると $u'=e^{-2t}=-2u+1$ であり、初期値も満たす。
:::
::::

::::exercise{level=advanced}
同じ生成作用素を持つ二つの強連続半群は一致することを示せ。
:::solution
$x\in D(A)$ の $s\mapsto T(t-s)S(s)x$ を微分すると $-T(t-s)AS(s)x+T(t-s)S(s)Ax=0$ となる。従って両端で $T(t)x=S(t)x$ であり、定義域の稠密性と有界性により全 $x$ で一致する。
:::
::::

::::exercise{level=advanced}
非負自己共役 $L$ で $A=-L$ が縮小半群を生成することを Lumer–Phillips で示せ。
:::solution
$\Re\inner{Ax}{x}=-\inner{Lx}{x}\le0$ なので散逸的である。$\lambda>0$ で $(\lambda+L)^{-1}$ は有界なスペクトル関数 $(\lambda+s)^{-1}$ であり、その像が $D(L)$ に入ることは $s/(\lambda+s)\le1$ による。従って全射条件が満たされる。
:::
::::

::::exercise{level=advanced}
複素 Banach 空間上の有界線形作用素 $A$ が生成する半群を示し、任意の複素時刻で正則であることを証明せよ。
:::solution
$e^{zA}=\sum_{n\ge0}z^nA^n/n!$ はコンパクト部分集合上でノルム一様収束する。Cauchy 積で加法法則が従い、零での微分は $A$ である。従って実非負時刻では強連続半群である。正則性はべき級数から従うが、全扇形で有界とは限らない。
:::
::::

::::exercise{level=advanced}
$T(t)e_n=e^{-n t}e_n$ のとき $t>0$ で $T(t)$ がコンパクトであることを示せ。
:::solution
先頭 $N$ 座標で打ち切った作用素との差のノルムは $e^{-(N+1)t}$ で零へ近づく。有限階作用素のノルム極限だからコンパクトである。$t=0$ では恒等作用素なので無限次元性によりコンパクトでない。
:::
::::

::::exercise{level=advanced}
$A=\diag(-2,-5,-8)$、$x=(1,2,0)$、定数外力 $f=(2,1,-1)$ とする。定数変化公式を計算し、定常値、時刻 $1/2$ の解、定常値への収束の評価を求めよ。
:::solution
**係数ごとの初期値と外力**

正規直交方向 $e_1,\ldots,e_3$ を取り、$Ae_j=-d_je_j$ とする。初期値 $\sum a_j e_j$ と時間に依存しない右辺 $\sum h_j e_j$ に対し、$u'(t)=Au(t)+\sum_j h_j e_j$ を解く。時刻は $t=\frac{1}{2}$ で比較する。$u_j'+d_ju_j=h_j$ に $e^{d_jt}$ を掛けて積分すると $u_j(t)=e^{-d_jt}a_j+(1-e^{-d_jt})h_j/d_j$ になる。

| モード | 減衰率 | 初期値 | 外力 | 解の係数 | 指定時刻の値 |
| --- | --- | --- | --- | --- | --- |
| 1 | $2$ | $1$ | $2$ | $1$ | $1$ |
| 2 | $5$ | $2$ | $1$ | $\frac{1}{5} + \frac{9 e^{- 5 t}}{5}$ | $\frac{9}{5 e^{\frac{5}{2}}} + \frac{1}{5}$ |
| 3 | $8$ | $0$ | $-1$ | $- \frac{1}{8} + \frac{e^{- 8 t}}{8}$ | $\frac{1 - e^{4}}{8 e^{4}}$ |

外力が一定なら定常解は $u_\infty=\sum_j(h_j/d_j)e_j$ である。差 $u_j-u_{\infty,j}=e^{-d_jt}(a_j-h_j/d_j)$ を平方和に入れると
$$
\norm{u(t)-u_\infty}\le e^{-2t}\norm{u(0)-u_\infty}
$$
となる。入力によって最初のモードの差が消えていれば、表に出てくる次の減衰率で改善できる。これは全入力に共通のノルム評価と、特定の入力に対して得られる評価の違いである。

同じ係数を定数変化公式で作ると、斉次項は初期値 $a_j$ に $e^{-d_jt}$ を掛けたもの、積分項は $\int_0^t e^{-d_j(t-s)}h_j\,ds=h_j(1-e^{-d_jt})/d_j$ である。両項は表の微分方程式の解と一致する。初期値と外力にゼロの座標を入れたモードは時間発展でも零を保つ。一般の入力を有限切断した場合には、切り落とした初期値のノルムと外力の積分を別々に評価する。
:::
::::

## 確認テスト

::::quiz
強連続半群の生成作用素について正しいものを選べ。

- [x] 稠密に定義され閉じている。
- [ ] 常に有界である。
- [ ] 定義域は常に全空間である。

:::explanation
時間平均で定義域の稠密性が得られ、軌道積分の恒等式の極限から閉性が従う。対角半群 $T(t)e_n=e^{-nt}e_n$ では生成作用素は有界でなく、定義域は重み付き平方和の条件を満たす真部分空間である。
:::
::::

::::quiz{answer="0.2" tol="1e-8"}
$A=-2I$ に対する $R(3,A)$ のノルムを求めよ。
:::explanation
$3I-A=5I$ であるから逆は $I/5$ であり、零でない空間でノルムは $1/5$ である。
:::
::::

::::quiz
縮小半群の Hille–Yosida 条件に使うノルム評価を選べ。

- [x] 全 $\lambda>0$ について $\norm{R(\lambda,A)}\le1/\lambda$
- [ ] 全 $\lambda>0$ について $\norm{R(\lambda,A)}\le\lambda$
- [ ] ある $\lambda>0$ について可逆であることだけ

:::explanation
縮小半群の Laplace 変換は $\int_0^\infty e^{-\lambda t}\,dt=1/\lambda$ によって評価される。逆方向でも各 $\lambda$ でこの評価があれば劣乗法性により冪の評価が従う。
:::
::::

::::quiz{answer="1" tol="1e-8"}
$u'=0u+1,\ u(0)=0$ の $u(1)$ を求めよ。
:::explanation
斉次半群は恒等写像なので定数変化公式は $u(t)=\int_0^t1\,ds=t$ を与える。従って $u(1)=1$ である。
:::
::::

::::quiz{answer="-2" tol="1e-8"}
$A e_j=-j^2e_j$ とし、$A_4=16R(4,A)-4I$ とする。$A_4e_2$ の係数を求めよ。
:::explanation
$R(4,A)e_2=(4+4)^{-1}e_2$ である。従って係数は $16/8-4=-2$ となる。
:::
::::

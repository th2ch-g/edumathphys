---
title: "ℝⁿ 上の Fourier 変換"
description: "Euclid 空間上の Fourier 変換を定義し、反転、Plancherel、補間、畳み込みの公式を導く。"
notation:
  - symbol: '\mathcal Ff'
    meaning: "指数の符号を負に取る Fourier 変換"
  - symbol: '\mathcal F^{-1}g'
    meaning: "指数の符号を正に取る逆変換"
  - symbol: '\mathcal S'
    meaning: "Schwartz 空間"
  - symbol: '\partial^\alpha'
    meaning: "多重指数による偏微分"
  - symbol: 'p'''
    meaning: "p の共役指数"
---

## この章の目標

- $L^1$ と Schwartz 空間上の Fourier 変換を計算する。
- Gauss 核で正則化して反転公式を証明する。
- Plancherel の定理と Hausdorff–Young の不等式を導く。
- 畳み込み、球対称性、周波数の尺度を計算に使う。

## L¹ 上の Fourier 変換

円周の周波数は整数だった。実空間では周期境界を課さないため、周波数を任意の実ベクトルに取る。$\xi$ は単位長さ当たりの周期数とし、指数に $2\pi$ を含める。

:::definition[Euclid 空間上の Fourier 変換]{#def-four-euclid-transform}
$f\in L^1(\R^n)$ の :term[Fourier 変換]{en="Fourier transform" yomi="ふーりえへんかん"} を

$$
\widehat f(\xi)=\mathcal Ff(\xi)=\int_{\R^n} f(x)e^{-2\pi i x\cdot\xi}\,dx
\label{eq-four-euclid-normalization}
$$

とする。$g\in L^1$ に対して $\check g(x)=\int g(\xi)e^{2\pi i x\cdot\xi}\,d\xi$ と書く。
:::

:::theorem[Riemann–Lebesgue と連続性]{#thm-four-euclid-rl}
$\widehat f\in C_0(\R^n)$ であり、$\|\widehat f\|_\infty\le\|f\|_1$ である。
:::

:::proof
絶対値を取ればノルム評価が得られる。連続性には優収束を使える。さらに $\widehat f(\xi+h)-\widehat f(\xi)$ の絶対値は $\int|f(x)||e^{-2\pi ix\cdot h}-1|\,dx$ 以下なので一様連続である。無限遠での消滅には $h=\xi/(2|\xi|^2)$ を取る。$\tau_h f(x)=f(x-h)$ の変換は $e^{-2\pi i h\cdot\xi}\widehat f(\xi)=-\widehat f(\xi)$ だから $2|\widehat f(\xi)|\le\|\tau_h f-f\|_1$ となる。$|\xi|\to\infty$ で $|h|\to0$ となり、:ref[prop-meas-translation] が使える。
:::

:::proposition[平行移動・変調・線形変換]{#prop-four-euclid-symmetries}
$a,b\in\R^n$、$A$ を可逆実行列とすると

$$
\begin{aligned}
\mathcal F(f(\cdot-a))(\xi)&=e^{-2\pi ia\cdot\xi}\widehat f(\xi),\\
\mathcal F(e^{2\pi ib\cdot x} f(x))(\xi)&=\widehat f(\xi-b),\\
\mathcal F(f(Ax))(\xi)&=|\det A|^{-1}\widehat f((A^{-1})^{\mathsf T}\xi),\\
\mathcal F(\bar f)(\xi)&=\overline{\widehat f(-\xi)}.
\end{aligned}
$$
:::

:::proof
順に $y=x-a$ の置換、指数の結合、$y=Ax$ の置換、積分の複素共役を用いる。線形変換では $A^{-1}y\cdot\xi=y\cdot(A^{-1})^{\mathsf T}\xi$ であり、Jacobian は $|\det A|^{-1}$ である。
:::

:::example[区間と指数減衰]
$1_{[-a,a]}$ の変換は $\sin(2\pi a\xi)/(\pi\xi)$ で、原点の値は $2a$ である。また $e^{-a|x|}$（$a>0$）の変換は、正負の半軸の積分を足して

$$
\frac1{a+2\pi i\xi}+\frac1{a-2\pi i\xi}
=\frac{2a}{a^2+4\pi^2\xi^2}
$$

となる。前者は点ごとの符号振動を持ち、後者は正で可積分である。
:::

:::example[箱の各軸の尺度]
$1_{[-1,1]\times[-2,2]}$ の変換は $\sin(2\pi\xi_1)\sin(4\pi\xi_2)/(\pi^2\xi_1\xi_2)$ である。積の形に Fubini を適用した。第 2 軸の最初の正の零点は $1/4$ で、第 1 軸の $1/2$ の半分である。
:::

### 変換が持つモーメントの情報

絶対積分の定義では、$x$ に大きな重みを掛けても積分できるほど、周波数側を多く微分できる。逆に、周波数側の減衰には入力の導関数が必要である。二つの条件は異なる方向の情報を与える。

:::proposition[モーメントと導関数]{#prop-four-moments-transform}
整数 $m\ge0$ とし $(1+|x|)^mf\in L^1$ なら、$\widehat f$ は $C^m$、各 $|\alpha|\le m$ で $\partial^\alpha\widehat f=\mathcal F((-2\pi ix)^\alpha f)$ である。これらの導関数は無限遠でゼロに近づく。さらに $f$ が原点を中心とする半径 $R$ の球に台を持つなら、一次元の場合の $\widehat f$ は整関数に延長され $|\widehat f(\xi+i\eta)|\le e^{2\pi R|\eta|}\|f\|_1$ を満たす。
:::

:::proof
一回の差商を取ると、指数の差商は $2\pi|x_j|$ 以下である。これを $|\alpha|$ 回繰り返すと $(2\pi |x|)^{|\alpha|}|f(x)|$ で支配できるため、微分と積分が交換できる。無限遠での消滅は各入力 $x^\alpha f$ に :ref[thm-four-euclid-rl] を使う。コンパクト台なら $e^{-2\pi ix(\xi+i\eta)}$ の絶対値は $e^{2\pi R|\eta|}$ 以下であり、局所的な複素差商も同じ型で支配されるので整関数になる。
:::

:::example[平均と一次モーメント]
$f=1_{[0,2]}$ では $\widehat f(0)=2$、$\widehat f'(0)=-2\pi i\int_0^2x\,dx=-4\pi i$ である。同じ質量を持つ $1_{[-1,1]}$ は奇な一次モーメントがゼロなので、周波数原点で傾きがゼロである。二つの入力は移動の位相で結び付く。
:::

### 変換の大きさに等号が成り立つ場合

積分で絶対値を外へ出す評価が鋭くなるには、各点から寄与する複素数の偏角がそろわなければならない。$f\in L^1$ が非零ならば :ref[thm-four-euclid-rl] により $|\widehat f|$ は無限遠で消える。したがって $\|\widehat f\|_\infty=\|f\|_1$ なら、その最大を実現する $\xi_0$ が存在する。$e^{-i\theta}\widehat f(\xi_0)=\|f\|_1$ となる $\theta$ に対して

$$
\int_{\R^n}\bigl(|f(x)|-\operatorname{Re}(e^{-i\theta}e^{-2\pi ix\cdot\xi_0}f(x))\bigr)\,dx=0
$$

である。被積分関数は非負だから a.e. ゼロになり、$f(x)=e^{i\theta}e^{2\pi ix\cdot\xi_0}|f(x)|$ となる。逆にこの形の $f$ を $\xi_0$ で積分すれば等号になる。非負入力は $\xi_0=0,\theta=0$ の例である。さらに非零の非負入力で $\xi\ne0$ なら $|\widehat f(\xi)|<\widehat f(0)$ である。等号なら $f>0$ の集合が $x\cdot\xi$ を法 $1$ で固定した可算個の超平面に含まれることになり、正の測度を持てないからである。

:::example[複素指数で絶対値の最大を移動する]
$f(x)=e^{-2|x|}e^{6\pi ix}$ とすると
$\widehat f(\xi)=4/(4+4\pi^2(\xi-3)^2)$ である。空間側で $|f|=e^{-2|x|}$、$\|f\|_1=1$ であり、周波数 $\xi=3$ だけで変換の絶対値は $1$ になる。$f$ にさらに空間移動を施すと周波数に依存する単位絶対値の因子が掛かる。このため、変換の絶対値だけでは入力の位置を復元できない。
:::

### 偶性と因子分解を積分の前に使う

一次元で $f$ の偶部分と奇部分をそれぞれ $f_e,f_o$ とすれば、

$$
\widehat f(\xi)=2\int_0^\infty f_e(x)\cos(2\pi x\xi)\,dx
-2i\int_0^\infty f_o(x)\sin(2\pi x\xi)\,dx .
$$

$\R$ の負の半軸を $x\mapsto-x$ で移せば従う。実数値の偶関数の変換は実偶関数、実数値の奇関数の変換は純虚奇関数である。一方、片側だけに台を持つ実数値関数はこの二つの成分をともに持つことがある。例えば $a>0$ で $f=1_{[0,\infty)}e^{-ax}$ なら $\widehat f=(a+2\pi i\xi)^{-1}$ であり、実部が $a/(a^2+4\pi^2\xi^2)$、虚部が $-2\pi\xi/(a^2+4\pi^2\xi^2)$ である。特に $a=2\pi$ として $\xi=1$ に代入すると $(1-i)/(4\pi)$ を得る。

:::example[座標のせん断と周波数のせん断]
$f(x_1,x_2)=g(x_1+x_2)h(x_2)$ とし $g,h\in L^1(\R)$ とする。$A=\begin{pmatrix}1&1\\0&1\end{pmatrix}$ に対して $(A^{-1})^{\mathsf T}\xi=(\xi_1,\xi_2-\xi_1)$ であり行列式は $1$ である。よって $\widehat f(\xi)=\widehat g(\xi_1)\widehat h(\xi_2-\xi_1)$ となる。例えば $g=1_{[-1/2,1/2]}$、$h(x)=e^{-2\pi|x|}$ なら $\xi=(1/2,1/2)$ では $\widehat g(1/2)=2/\pi$、$\widehat h(0)=1/\pi$ なので変換の値は $2/\pi^2$ である。周波数の第 2 座標に現れる差は、空間変換をそのまま周波数へ代入しても得られない。
:::

### 非整数のモーメントと一様な Hölder 評価

$0<\alpha\le1$、$\int|x|^\alpha |f(x)|\,dx<\infty$ なら

$$
|\widehat f(\xi+h)-\widehat f(\xi)|
\le 2^{1-\alpha}(2\pi|h|)^\alpha\int |x|^\alpha|f(x)|\,dx .
$$

実際、$|e^{iu}-1|\le\min(2,|u|)\le2^{1-\alpha}|u|^\alpha$ を積分に適用する。$x^\beta f$ に同じ評価を使えば $\partial^\beta\widehat f$ にも Hölder 評価を与えられる。例えば $f(x)=(1+|x|)^{-7/4}$ なら $\int|x|^\alpha f<\infty$ が $\alpha<3/4$ で成り立つので、各 $\alpha<3/4$ の一様 Hölder 連続性を得る。一方、$|x| f$ は可積分でなく、このモーメントによる方法では一般の一階微分を正当化できない。必要なモーメントを調べずに積分の中で微分すると、原点のように位相の相殺がない点で特に問題となる。

## Schwartz 空間と Fourier 変換

微分と多項式倍を繰り返しても無限遠で十分速く減る空間を用意する。この空間では積分と微分が交換でき、反転の証明後に他の空間への拡張に使える。

:::definition[Schwartz 空間]{#def-four-schwartz}
多重指数 $\alpha=(\alpha_1,\ldots,\alpha_n)\in\Z_{\ge0}^n$ について $x^\alpha=\prod x_j^{\alpha_j}$、$\partial^\alpha=\prod\partial_{x_j}^{\alpha_j}$、$|\alpha|=\sum\alpha_j$ とする。:term[Schwartz 空間]{en="Schwartz space" yomi="しゅわるつくうかん"} $\mathcal S(\R^n)$ は、すべての $\alpha,\beta$ で $p_{\alpha,\beta}(f)=\sup_x|x^\alpha\partial^\beta f(x)|<\infty$ となる $C^\infty$ 関数の空間である。これらのセミノルムで位相を定める。
:::

:::proposition[微分と多項式倍]{#prop-four-schwartz-operations}
$f\in\mathcal S$ に対し

$$
\mathcal F(\partial^\alpha f)=(2\pi i\xi)^\alpha\widehat f,\qquad
\partial_\xi^\beta\widehat f=\mathcal F((-2\pi ix)^\beta f).
$$

$\mathcal F$ は $\mathcal S\to\mathcal S$ で連続である。
:::

:::proof
座標ごとの部分積分で第 1 式を得る。各境界項は多項式より速い減衰で消え、残る積分は絶対可積分である。第 2 式は各階の微分を $C|x|^{|\beta|}|f(x)|$ で支配して積分内で行う。両式から $|\xi^\alpha\partial^\beta\widehat f(\xi)|\le C_{\alpha,\beta}\|\partial^\alpha(x^\beta f)\|_1$ となる。積の微分は有限個の $x^\gamma\partial^\delta f$ の和である。それぞれに $(1+|x|)^{-n-1}$ を掛け、残りを有限個のセミノルムで抑えられるので連続性を得る。
:::

:::example[微分が作る周波数の零点]
$g(x)=e^{-\pi x^2}$ は各導関数が多項式と $g$ の積なので $\mathcal S$ に属する。$g''$ の変換は $-4\pi^2\xi^2\widehat g(\xi)$ で、原点で二次の零点を持つ。特に $\int g''=0$ である。
:::

:::example[積分可能性と Schwartz 条件]
$f(x)=(1+x^2)^{-2}$ は $L^1(\R)$ に属するが、$x^5 f(x)$ は有界でない。一方コンパクト台を持つ $C^\infty$ 関数は台の外で全導関数が消えるので $\mathcal S$ に属する。
:::

### Schwartz 位相での収束

各導関数について急減少を見るには、座標ごとの単項式より $1+|x|$ を使う方が便利な場合がある。$q_{m,k}(f)=\max_{|\beta|\le k}\sup_x(1+|x|)^m|\partial^\beta f(x)|$（$m,k\in\Z_{\ge0}$）とする。$|x^\alpha|\le (1+|x|)^{|\alpha|}$ なので $p_{\alpha,\beta}\le q_{|\alpha|,|\beta|}$ である。逆に $(1+|x|)^m\le (1+\sum_j|x_j|)^m$ を展開すると、$q_{m,k}$ は $|\alpha|\le m,|\beta|\le k$ の $p_{\alpha,\beta}$ の有限和の定数倍で抑えられる。従って両者は同じ位相を定める。$\|f\|_p\le C q_{m,0}(f)$ は $mp>n$（$p<\infty$）なら $(1+|x|)^{-mp}$ の積分から従い、$p=\infty$ なら $m=0$ で成り立つ。

:::proposition[Schwartz 空間の完備性と切断近似]{#prop-four-schwartz-topology}
Schwartz 空間は可算個のセミノルムを持つ完備な距離化可能局所凸空間であり、$C_c^\infty(\R^n)$ がその位相で稠密である。
:::

:::proof
距離は例えば $\sum_{j=0}^\infty 2^{-j-1}\min(1,q_{j,j}(f-g))$ で与えられ、この距離での Cauchy 性は全ての $q_{m,k}$ についての Cauchy 性と同値である。$f_\ell$ が Cauchy 列なら、各 $\partial^\beta f_\ell$ は一様に連続関数 $g_\beta$ へ近づく。線分上で

$$
\partial^\beta f_\ell(x+he_j)-\partial^\beta f_\ell(x)
=\int_0^h\partial^{\beta+e_j}f_\ell(x+se_j)\,ds
$$

と書いて極限を取れば、$\partial_jg_\beta=g_{\beta+e_j}$ である。従って $g_0$ は滑らかで、全ての $g_\beta$ がその導関数である。各重み付き一様ノルムでも Cauchy だから、点ごとの極限 $(1+|x|)^m g_\beta(x)$ は有界で、一様収束もする。このため $g_0\in\mathcal S$ であり $f_\ell\to g_0$ が Schwartz 位相で成り立つ。

稠密性には、$|x|\le1$ で $1$、$|x|\ge2$ で $0$ である滑らかな関数 $\chi$ を取り $f_R(x)=\chi(x/R)f(x)$ と置く。$\partial^\beta(f-f_R)$ の積の微分の各項は $|x|\ge R$ に台を持つ。$\chi(x/R)$ の階数 $j$ の導関数は $C R^{-j}$ 以下であり、$f$ の導関数には $(1+|x|)^{-m-1}q_{m+1,k}(f)$ の評価を使える。$R\ge1$ で $q_{m,k}(f-f_R)\le C_{m,k,\chi}(1+R)^{-1}q_{m+1,k}(f)$ となり、結論が従う。
:::

:::example[指数で消えるが滑らかでない入力]
$e^{-|x|}$ は無限遠で任意の多項式より速く減る。しかし原点の左右で導関数は $+1,-1$ となるので $\mathcal S$ には属さない。一方、$e^{-x^2}$ は $\partial^j e^{-x^2}=P_j(x)e^{-x^2}$ と書け、$P_0=1,P_{j+1}=P_j'-2xP_j$ である。例えば $P_2=4x^2-2$ である。すべての $P_j$ は多項式なので $q_{m,k}$ は有限になる。無限遠での入力の減衰だけでなく、各導関数について同じ性質が必要である。
:::

### 変換と微分作用素の交換順序

位置の乗算 $X_jf=x_jf$ と微分 $\partial_j$ は交換しない。実際、積の微分により $(\partial_j X_k-X_k\partial_j)f=\delta_{jk}f$ である。Fourier 変換後は $\partial_j$ が $2\pi i\xi_j$ の乗算、$X_k$ が $-(2\pi i)^{-1}\partial_{\xi_k}$ になる。右辺が失われないのは、移った先で乗算と微分の順序が逆になるためである。例えば一次元で

$$
\begin{aligned}
\widehat{(xf)'}(\xi)&=-\xi(\widehat f)'(\xi),\\
\widehat{xf'}(\xi)&=-\widehat f(\xi)-\xi(\widehat f)'(\xi)
\end{aligned}
$$

である。最初の式では先に $xf$ を変換してから乗数を掛け、第二式では先に $f'$ を変換してから全体を微分している。$f=e^{-\pi x^2}$ なら右辺はそれぞれ $2\pi\xi^2e^{-\pi\xi^2}$ と $(2\pi\xi^2-1)e^{-\pi\xi^2}$ になる。この順序は後で微分作用素の固有関数を求める際にも使う。

### 変換の減衰を微分の積分から読む

$f\in C_c^\infty(\R^n)$ なら
$(1+4\pi^2|\xi|^2)^k\widehat f=\widehat{(1-\Delta)^k f}$
なので、絶対値は $\|(1-\Delta)^k f\|_1$ 以下である。$k>n/2$ を取れば変換は可積分になる。実際、$|\xi|\ge1$ で逆の重みは $C|\xi|^{-2k}$ 以下であり、極座標で積分できる。有限回の微分しか持たない場合も、境界項を消せる回数まで同じ計算ができる。例えば一次元で $f$ がコンパクト台を持つ $C^2$ 関数なら $|\widehat f(\xi)|\le(1+4\pi^2\xi^2)^{-1}\|f-f''\|_1$ となり、$C^\infty$ を仮定せず反転の絶対積分が得られる。一方、区間の指示関数では最初の部分積分で端の項が残るため、減衰は $O(|\xi|^{-1})$ に止まる。

## 反転公式

$\widehat f$ が可積分なら指数の符号を反転した積分が定義できる。その積分が元へ戻ることを示すため、両方の積分を交換できるように周波数側に Gauss 因子を一度挿入する。

:::lemma[正規化した Gauss 関数]{#lem-four-gaussian}
$g(x)=e^{-\pi |x|^2}$ なら $\widehat g=g$ である。$g_\eps(x)=\eps^{-n}g(x/\eps)$ なら $\widehat{g_\eps}(\xi)=e^{-\pi\eps^2|\xi|^2}$ である。
:::

:::proof
一次元で $G=\widehat g$ と置く。$g'=-2\pi x g$ と :ref[prop-four-schwartz-operations] より $G'=-2\pi\xi G$ である。$G(0)=1$ は、積分の二乗を極座標にして $\int_{\R^2}e^{-\pi(x^2+y^2)}dxdy=2\pi\int_0^\infty re^{-\pi r^2}dr=1$ と計算して得る。この初期値の微分方程式の解は $G=e^{-\pi\xi^2}$ である。各座標への分解で $n$ 次元に移り、拡大縮小で最後の式を得る。
:::

:::theorem[$L^1$ 反転公式]{#thm-four-euclid-inversion}
$f\in L^1$ かつ $\widehat f\in L^1$ なら $f=\check{\widehat f}$ a.e. である。右辺は連続であり、$f$ の連続点の値と一致する。さらに $\mathcal F^2 f(x)=f(-x)$ が $\mathcal S$ 上で成り立ち、$\mathcal F$ は $\mathcal S$ の位相同型である。
:::

:::proof
Fubini と :ref[lem-four-gaussian] により

$$
\begin{aligned}
\int \widehat f(\xi)e^{-\pi\eps^2|\xi|^2}e^{2\pi ix\cdot\xi}\,d\xi
&=\int f(y)g_\eps(x-y)\,dy.
\end{aligned}
$$

右辺は近似単位により $L^1$ で $f$ に収束する（:ref[prop-meas-approx-error]）。これは $g_\eps$ の質量が $1$ で固定半径の外側の質量が消えることによる。左辺は優収束を $\int|\widehat f||e^{-\pi\eps^2|\xi|^2}-1|$ に適用して $\check{\widehat f}$ に一様収束する。従ってコンパクト集合上で両極限が $L^1$ で一致する。連続点で異なる値を取るなら近傍で差が正の下限を持ち、a.e. の一致に矛盾する。Schwartz 関数とその変換は可積分かつ連続であり、反転を全点で使える。$\mathcal F^{-1}$ は $\mathcal F$ と反射の合成なので連続である。
:::

$f\in L^1$ の変換がゼロなら :ref[thm-four-euclid-inversion] により $f=0$ a.e. である。さらに $\widehat f$ の可積分性がなくても、上の Gauss 因子付き積分は $L^1$ で $f$ に戻る。

:::example[Poisson 核と両側指数]
$e^{-2\pi y|\xi|}$（$y>0$）の逆変換は $y/(\pi(y^2+x^2))$ である。先の両側指数の式に $a=2\pi y$ を入れ、偶性で正負の符号の違いを消した。両方が可積分なので反転からこの有理関数の変換は $e^{-2\pi y|\xi|}$ となる。
:::

:::example[Gauss 正則化による復元の幅]
一次元で $f=e^{-\pi x^2}$ を Gauss 因子付きで反転すると、積の変換は $e^{-\pi(1+\eps^2)\xi^2}$ である。従って正則化された元の関数は $(1+\eps^2)^{-1/2}e^{-\pi x^2/(1+\eps^2)}$ になる。$\eps=1$ では原点の値が $1/\sqrt2$ に下がるが、尺度変換により積分値は $1$ のままである。$\eps\downarrow0$ で点の値と幅がともに戻る。
:::

### Lebesgue 点と跳躍での Gauss 正則化

Gauss 因子を掛けた反転は、$\widehat f\in L^1$ を仮定せず各点で考えられる。$f\in L^1$ なら $\widehat f$ が有界であり、Gauss 因子は可積分だからである。$x$ が $f$ の Lebesgue 点なら、その積分は $f(x)$ に近づく。

:::proof[Gauss 正則化の各点収束]
$h(t)=|f(x-t)-f(x)|$ とし、$\eta>0$ に対して $0<a\le2\delta$ で $\int_{|t|\le a}h(t)\,dt\le\eta a^n$ となる $\delta>0$ を選ぶ。$|t|\le\eps$ で $g_\eps\le\eps^{-n}$ なので積分は $\eta$ 以下である。$2^j\eps<|t|\le2^{j+1}\eps$ かつ $|t|<\delta$ で $g_\eps\le \eps^{-n}e^{-\pi 4^j}$ であり、$h$ の積分は外半径を使って $\eta(2^{j+1}\eps)^n$ 以下になる。そのため全内側寄与は $\eta$ の定数倍以下である。$|t|\ge\delta$ では $f(x-t)$ の項を $\|f\|_1\sup_{|t|\ge\delta} g_\eps(t)$ で、$f(x)$ の項を $|f(x)|\int_{|t|\ge\delta}g_\eps(t)dt$ で抑える。それぞれ $\eps\downarrow0$ でゼロになる。最後に $\eta\downarrow0$ とする。
:::

:::example[矩形の端点と Gauss 誤差関数]
$f=1_{[0,1]}$ とし、$x=\eps y$ で正則化した値を取る。$u=(x-t)/\eps$ と置けば

$$
(f*g_\eps)(\eps y)=\int_{y-1/\eps}^y e^{-\pi u^2}\,du
\longrightarrow \int_{-\infty}^y e^{-\pi u^2}\,du .
$$

$y=0$ で極限は $1/2$、拡大した波形は導関数 $e^{-\pi y^2}>0$ を持つ。Gauss 正則化では平均が $[0,1]$ に保たれ、跳躍の両側へはみ出さない。$x\in(0,1)$ を固定すると積分区間は全実軸を覆うように広がり、復元値は $1$ へ近づく。
:::

### 平滑化の二次の主項

$f$ が $C^2$、有界で、二階導関数が有界かつ一様連続なら、Gauss で平滑化したときの最初の誤差を一様に計算できる。$g=e^{-\pi|x|^2}$ に対し

$$
f*g_\eps=f+\frac{\eps^2}{4\pi}\Delta f+o(\eps^2)
$$

が上限ノルムで成り立つ。実際、$t=\eps u$ として $f(x-\eps u)$ を $x$ の周りで二次まで展開する。一次の項は偶性により消える。$\int u_j u_k g(u)\,du=0$（$j\ne k$）も座標ごとの奇性から従う。$\int u_j^2g(u)\,du=1/(2\pi)$ は一変数で $(e^{-\pi u_j^2})'=-2\pi u_j e^{-\pi u_j^2}$ を使う部分積分による。積分型剰余を $\eps^2|u|^2\omega(\eps|u|)/2$ で抑えると、$\omega$ は二階導関数の一様連続度なので、$u$ について優収束を適用できる。例えば $f=e^{-\pi x^2}$ の原点では $f''(0)=-2\pi$ だから平滑化は $1-\eps^2/2+o(\eps^2)$ であり、$(1+\eps^2)^{-1/2}$ の展開と一致する。

### 有限測度に対する一意性

絶対積分を有限複素 Borel 測度 $\mu$ に対して $\widehat\mu(\xi)=\int e^{-2\pi ix\cdot\xi}\,d\mu(x)$ と定められる。優収束で連続になり、絶対値は全変動 $\|\mu\|$ 以下である。ただし無限遠で消えるとは限らない。例えば $\mu=\delta_a$ では指数そのものが変換になり、絶対値は全周波数で $1$ である。

:::proposition[有限測度を変換から区別する]{#prop-four-measure-transform-injectivity}
有限複素 Borel 測度 $\mu,\nu$ の変換が一致するなら $\mu=\nu$ である。
:::

:::proof
差を取って $\widehat\mu=0$ とする。Gauss 因子を挿入した反転と同じ Fubini の計算により、$\mu*g_\eps(x)=\int g_\eps(x-y)\,d\mu(y)$ は全ての $x$ でゼロになる。$g_\eps\in L^1$ と $\|\mu\|<\infty$ が交換を保証する。任意の $\varphi\in C_c(\R^n)$ に対し、この関数を $\varphi(x)$ と掛け積分する。Fubini と Gauss の偶性により $0=\int(\varphi*g_\eps)(y)\,d\mu(y)$ となる。$\varphi*g_\eps\to\varphi$ が一様に成り立つのは、$\varphi$ の一様連続性と $g_\eps$ の質量集中による。従って $\int\varphi\,d\mu=0$ である。コンパクト台の連続関数は有限 Borel 測度を区別する。例えばコンパクト集合 $K$ の指示関数を、$K$ で $1$、距離 $1/m$ 以上で $0$ となる $[0,1]$ 値の連続関数で点ごとに近似することで $\mu(K)=0$ が従い、全変動の正則性で全 Borel 集合に移る。
:::

:::example[有限個の点質量を区別する]
$\mu=(\delta_{-1}+\delta_1)/2$ では $\widehat\mu(\xi)=\cos(2\pi\xi)$ である。正則化した密度は $(g_\eps(x+1)+g_\eps(x-1))/2$ になる。それぞれの峰の幅が $\eps$ に比例して小さくなり、コンパクト台の連続関数で積分すれば、その関数の $-1,1$ での値の平均に近づく。有限測度の復元ではこの積分での収束を使う。各点の密度の有限値への収束は、もとの測度が密度を持たない場合には要求できない。
:::

## Plancherel の定理

$L^2$ 関数の絶対積分が存在するとは限らない。そこで、先に Schwartz 関数で変換が二乗ノルムを保存することを示し、ノルム極限によって定義を拡張する。

:::theorem[Plancherel]{#thm-four-plancherel}
$\mathcal F$ は $L^2(\R^n)$ 上のユニタリ作用素に一意に拡張される。$f\in L^1\cap L^2$ では積分で定めた変換と a.e. 一致する。
:::

:::proof
$f,g\in\mathcal S$ に対し、$g$ の反転式を共役して $\int f\bar g$ に入れる。$\iint|f(x)||\widehat g(\xi)|\,dx\,d\xi$ は有限なので

$$
\int f(x)\overline{g(x)}\,dx=\int \widehat f(\xi)\overline{\widehat g(\xi)}\,d\xi
\label{eq-four-plancherel}
$$

が従う。:ref[thm-meas-smooth-density] で $\mathcal S$ は $L^2$ に稠密である。等長性から変換の極限は近似列によらず、像は完備なので閉じている。:ref[thm-four-euclid-inversion] により像は $\mathcal S$ を含むため全射である。

$f\in L^1\cap L^2$ では、先に空間と振幅を切断し、次に滑らかなコンパクト台の近似単位で畳み込めば両ノルムで収束する $\mathcal S$ の近似列を作れる。変換は $L^1$ に関する評価で一様収束し、等長性で $L^2$ 収束する。コンパクト集合上では両方が $L^2$ で同じ極限になる。
:::

:::example[矩形の二乗エネルギー]
$f=1_{[-1/2,1/2]}$ の二乗ノルムは $1$ である。従って $\int_{\R}(\sin(\pi\xi)/(\pi\xi))^2\,d\xi=1$ である。被積分関数は非負で、振動する積分を直接評価せずに値を得た。
:::

:::example[絶対積分を持たない入力]
$f(x)=(1+|x|)^{-3/4}$ は $L^2$ だが $L^1$ に属さない。切断した $f\,1_{[-R,R]}$ の積分変換の $L^2$ 極限で $\widehat f$ を定義できる。その二乗ノルムは $2\int_0^\infty(1+x)^{-3/2}dx=4$ である。
:::

### 二つの定義域の和に変換を延長する

$f=f_1+f_2$、$f_1\in L^1$、$f_2\in L^2$ と書けるとき、$\mathcal Ff=\widehat f_1+\widehat f_2$ と置ける。この和は分解によらない。実際、二つの分解の第一項の差は、第二項の差でもあるため $L^1\cap L^2$ に属し、:ref[thm-four-plancherel] によりどちらの変換を使っても同じである。同様に、$L^1+L^2$ への逆変換も分解によらず定義できる。例えば $f(x)=|x|^{-3/4}1_{\{|x|\le1\}}+(1+|x|)^{-3/4}1_{\{|x|>1\}}$ は $L^1$ にも $L^2$ にも入らないが、第一項は $L^1$、第二項は $L^2$ に入るので、この変換を持つ。

### 空間側と周波数側の切断

$f\in L^2$ では $\widehat f$ の値はほとんど至る所で定まる。積分式から $\widehat f$ を計算したいときは、まず $f_R=f\,1_{\{|x|\le R\}}$ を取る。Cauchy–Schwarz により $f_R\in L^1\cap L^2$ であり、$\|f_R-f\|_2\to0$ である。従って $\widehat f_R$ は積分で求まり、$\widehat f$ に $L^2$ で近づく。逆に $\widehat f\,1_{\{|\xi|\le R\}}$ は $L^1\cap L^2$ であるため

$$
P_R f(x)=\int_{|\xi|\le R}\widehat f(\xi)e^{2\pi ix\cdot\xi}\,d\xi
$$

は連続関数を定める。Plancherel により $\|P_R f-f\|_2^2=\int_{|\xi|>R}|\widehat f(\xi)|^2d\xi$ である。$P_R$ は帯域が $|\xi|\le R$ の閉部分空間への直交射影である。実際、周波数側の指示関数の乗算が直交射影であり、変換でその内積が保存される。同じ帯域の $g$ を用いる任意の近似に対し $\|f-g\|_2^2=\|f-P_Rf\|_2^2+\|P_Rf-g\|_2^2$ である。切断で得た誤差が最小であることと、各点で復元できることは別の主張である。

:::example[帯域の端と空間の裾]
一次元で $\widehat f=1_{[-B,B]}$ とすると、$f(x)=\sin(2\pi Bx)/(\pi x)$（原点で $2B$）は $L^2$ の逆変換である。$0<R<B$ なら $P_R f$ は $B$ を $R$ に替えた同じ形で、$\|f-P_R f\|_2^2=2(B-R)$ である。例えば $B=2,R=1/2$ では二乗誤差は $3$ だが、原点での誤差も $4-1=3$ となる。この数値の一致は選んだ $B,R$ によるもので、前者は周波数領域の測度、後者は積分の振幅の差から出る。
:::

### 平行移動と自己相関

$f,g\in L^2$ に対し $C_{f,g}(y)=\int f(x+y)\bar g(x)\,dx$ と置く。$L^1\cap L^2$ の稠密性と等長性により $\widehat{f(\cdot+y)}(\xi)=e^{2\pi i y\cdot\xi}\widehat f(\xi)$ が $L^2$ 上でも成り立つ。従って Parseval により

$$
C_{f,g}(y)=\int \widehat f(\xi)\overline{\widehat g(\xi)}e^{2\pi i y\cdot\xi}\,d\xi .
$$

積は Cauchy–Schwarz で可積分である。従って $C_{f,g}$ は連続で無限遠でゼロに近づき、絶対値は $\|f\|_2\|g\|_2$ 以下になる。$g=f$ では $\|f\|_2^2=C_{f,f}(0)$ で、自己相関の逆変換側の密度は $|\widehat f|^2$ である。この関係には $f\in L^1$ を仮定する必要がない。

:::example[帯域制限関数の自己相関]
$\widehat f=1_{[-1,1]}$ とすると $|\widehat f|^2=\widehat f$ なので $C_{f,f}(y)=\sin(2\pi y)/(\pi y)=f(y)$ となる。原点では $C_{f,f}(0)=2$、$y=1/2$ では $0$ である。後者は $f(\cdot+1/2)$ と $f$ が直交することを意味する。どちらも空間では全実軸に裾を持つが、位相を含めた積分がゼロになる。
:::

### 二乗ノルムでの微分の判定

通常の導関数が存在しない場合にも、変換側の $\xi$ の重みで微分の有無を測れる。$f,g\in L^2$ に対し、すべての $\varphi\in C_c^\infty$ について $\int g\bar\varphi=-\int f\overline{\partial_j\varphi}$ が成り立つとき、$g$ を $f$ の第 $j$ 方向の :term[弱微分]{en="weak derivative" yomi="じゃくびぶん"} という。

:::proposition[弱微分と周波数の重み]{#prop-four-l2-weak-derivative}
$f\in L^2$ が $L^2$ に属する第 $j$ 弱微分を持つことと、$\xi_j\widehat f\in L^2$ であることは同値である。そのとき
$\widehat{\partial_j f}=2\pi i\xi_j\widehat f$
である。
:::

:::proof
弱微分の等式は、Schwartz 関数を滑らかな切断で近似すると $\varphi\in\mathcal S$ に対しても成り立つ。実際、近似は $\varphi$ と $\partial_j\varphi$ の両方について $L^2$ で成り立ち、積分の極限は Cauchy–Schwarz で取れる。Plancherel と Schwartz の微分則で等式を変換すると

$$
\int\widehat g\,\overline{\widehat\varphi}\,d\xi
=\int(2\pi i\xi_j\widehat f)\overline{\widehat\varphi}\,d\xi
$$

となる。右辺は $\widehat f\in L^2$ と $\xi_j\widehat\varphi\in L^2$ により絶対収束する。任意の $\psi\in C_c^\infty$ に対し $\widehat\varphi=\psi$ と選べる。従って局所可積分関数 $\widehat g-2\pi i\xi_j\widehat f$ の積分が全試験関数でゼロになる。これは、その関数を滑らかな近似単位で正則化するとゼロになり、局所 $L^1$ で元へ戻ることから a.e. の一致を意味する。ゆえに $\xi_j\widehat f\in L^2$ である。

逆に $\xi_j\widehat f\in L^2$ なら、$2\pi i\xi_j\widehat f$ の $L^2$ 逆変換 $g$ を取る。上の等式を逆に辿れば、各 $\varphi\in C_c^\infty$ に対し弱微分の等式が得られる。
:::

:::example[角と段差が弱微分に与える差]
$f(x)=(1-|x|)_+$ の弱微分は $(-1,0)$ 上で $1$、$(0,1)$ 上で $-1$、外側で $0$ である。区間別の部分積分で端の値が相殺するので、この関数が弱微分の等式を満たし、ノルムの二乗は $2$ である。従って $4\pi^2\int \xi^2(\sin\pi\xi/(\pi\xi))^4d\xi=2$ である。一方 $g=1_{[-1/2,1/2]}$ では $\xi\widehat g(\xi)=\sin(\pi\xi)/\pi$ は二乗可積分でない。角では値がつながり、各区間の傾きが弱微分を与えるが、値そのものが跳ぶと $L^2$ の弱微分は得られない。
:::

### 周波数に有限の台を持つ正の正則化

Gauss 因子は全周波数で正である。厳密に有限の周波数だけを使いながら、空間で非負の重みを得るには、一次元の矩形を周波数側で畳み込む。$R>0$ とし、$m_R(\xi)=(1-|\xi|/R)_+$ と置く。$m_R=R^{-1}1_{[-R/2,R/2]}*1_{[-R/2,R/2]}$ であることは二つの区間の重なりの長さから従う。従って逆変換は

$$
K_R(x)=\frac1R\left(\frac{\sin(\pi R x)}{\pi x}\right)^2
=R\left(\frac{\sin(\pi R x)}{\pi R x}\right)^2 .
$$

積への変換は、二つの積分を Fubini で交換して得られる。非負な $K_R$ の質量は、Plancherel で二乗積分を矩形の長さ $R$ に帰着すれば $1$ である。空間側の評価は $K_R(x)\le\min(R,1/(\pi^2Rx^2))$ であるため、固定された $\delta>0$ の外側質量は $2/(\pi^2R\delta)$ 以下になる。よって $K_R$ は $R\to\infty$ の近似単位である。$m_R$ と $K_R$ がともに可積分なので、反転公式により $\widehat K_R=m_R$ でもある。

$f\in L^1$ に対して、$m_R\widehat f$ の逆変換と $K_R*f$ は Fubini により全点で一致する。$f\in L^p$（$1\le p<\infty$）に対して $K_R*f\to f$ が $L^p$ で成り立つ。これは $L^p$ での移動連続性を使い、畳み込みの誤差を $|x|<\delta$ と外側に分けることで従う。$L^1$ 関数の Lebesgue 点での各点収束には、まず $|x|\le R^{-1}$ の誤差を $K_R\le R$ で、次に $2^j/R<|x|\le2^{j+1}/R$ の誤差を $K_R\le R/(\pi^2 4^j)$ で抑える。Lebesgue 点の差の積分が半径に対し $o(r)$ であるため、十分小さな $\delta$ までの環の和は任意に小さくできる。外側では $f$ の項は $O_\delta(R^{-1})\|f\|_1$、選んだ点の値の項は外側質量で消える。このように、帯域が有限という条件を保った正則化でも一つの点を復元できる。

:::example[鋭い切断と線形な切断の違い]
$\widehat f=1_{[-B,B]}$、$R\ge B$ とする。鋭い切断では $P_R f=f$ である。線形な切断では $1-m_R(\xi)=|\xi|/R$ が帯域内でも正なので、

$$
\|K_R*f-f\|_2^2
=\frac1{R^2}\int_{-B}^B\xi^2\,d\xi=\frac{2B^3}{3R^2}
$$

である。$B=R=1$ では二乗誤差は $2/3$ になる。一方、空間での $K_R$ は非負で質量 $1$ だから、例えば $0\le g\le1$ の $L^1$ 入力には $0\le K_R*g\le1$ が成り立つ。帯域内の周波数を正確に保存する性質と、空間で入力の値域を保存する性質には違いがある。
:::

## Hausdorff–Young の不等式

二つの端点で、変換先はそれぞれ $L^\infty$ と $L^2$ であった。入力の指数を $1$ から $2$ に動かすと、出力は共役指数へ動く。

:::theorem[Hausdorff–Young]{#thm-four-hausdorff-young}
$1\le p\le2$、$p'=p/(p-1)$（$p=1$ で $p'=\infty$）とする。Fourier 変換は $L^p\to L^{p'}$ へ拡張され、

$$
\|\widehat f\|_{p'}\le\|f\|_p
$$

が成り立つ。
:::

:::proof
:ref[thm-four-euclid-rl] と :ref[thm-four-plancherel] の作用素ノルムは高々 $1$ である。共通定義域で一致するので :ref[thm-meas-riesz-thorin] を用い、$1/p=(1-\theta)+\theta/2$、$1/p'=\theta/2$ と置く。コンパクト台の有界関数は $L^p$ に稠密なので拡張は一意である。$L^1\cap L^p$ との一致は切断による両ノルム近似から、$L^1$ と $L^2$ の場合と同様に従う。
:::

:::example[指数の確認]
$p=4/3$ なら $p'=4$ であり、$f=1_{[-1/2,1/2]}$ では $\|\widehat f\|_4\le1$ である。実際 $f*f(x)=(1-|x|)_+$ なので次節の積公式と Plancherel で $\|\widehat f\|_4^4=2/3$ となる。上界 $1$ はこの例で等号を与えない。
:::

:::example[Gauss 関数のノルム比]
一次元の $g=e^{-\pi x^2}$ では $\|g\|_p=p^{-1/(2p)}$ である。$p=4/3$ で $\|\widehat g\|_{4}/\|g\|_{4/3}=4^{-1/8}(4/3)^{3/8}<1$ になる。補間の定数は端点の情報から得られ、内部の最良定数の決定は別の問題である。
:::

### 入力指数を二より大きくしたとき

$p=2$ を越えて同じ共役指数への評価が成り立つと仮定すると、離れた場所の入力を並べた例で矛盾が生じる。空間でのノルムは各項を独立に数えるが、変換では位相の多項式になり、その大きさを下から評価できる。

:::proposition[Hausdorff–Young の指数範囲]{#prop-four-hausdorff-young-range}
$2<p\le\infty$ では、$f\in C_c^\infty(\R^n)$ に対する評価 $\|\widehat f\|_{p'}\le C\|f\|_p$ を満たす共通の定数 $C$ は存在しない。
:::

:::proof
一変数で二つの多項式 $P_0=Q_0=1$ を取り
$$
P_{k+1}(z)=P_k(z)+z^{2^k}Q_k(z),\qquad
Q_{k+1}(z)=P_k(z)-z^{2^k}Q_k(z)
$$
と置く。次数 $2^k-1$ までの係数はいずれも $\pm1$ である。前半と後半が重ならないためである。また $|z|=1$ では平行四辺形の法則で $|P_k|^2+|Q_k|^2=2^{k+1}$ である。$N=2^k$ と書くと $\|P_k\|_\infty\le\sqrt{2N}$、$\int_0^1|P_k(e^{2\pi it})|^2dt=N$ となる。$q=p'<2$ では $|P_k|^2\le(2N)^{(2-q)/2}|P_k|^q$ なので $\|P_k(e^{2\pi it})\|_{L^q(0,1)}\ge c_q\sqrt N$ である。

台が $(-1/4,1/4)$ に入る非負非零の滑らかな関数 $g$ を選ぶ。$\widehat g(0)>0$ なので $[-b,b]$ で $|\widehat g|\ge c>0$ となる $b>0$ がある。$A>\max(1,1/b)$ を固定し、$f_k(x)=\sum_{j=0}^{N-1}a_jg(x-jA)$ と置く。ただし $a_j$ は $P_k$ の係数である。台が互いに素なので $\|f_k\|_p=N^{1/p}\|g\|_p$ である。変換は $\widehat g(\xi)P_k(e^{-2\pi iA\xi})$ であり、$[0,1/A]\subset[-b,b]$ 上で $L^q$ ノルムは $cA^{-1/q}c_q\sqrt N$ 以上である。比は $N^{1/2-1/p}$ の定数倍以上なので発散する。高次元では他の座標に固定された非零の滑らかなコンパクト台関数を掛ければ、両ノルムに固定因子が掛かるだけである。
:::

:::example[四項の符号の選び方]
$P_1=1+z$、$Q_1=1-z$ なので $P_2=1+z+z^2-z^3$ である。係数の二乗和は $4$、単位円上の絶対値は $\sqrt8$ 以下となる。この上限は全ての係数を $1$ にして $z=1$ で値 $4$ を取る多項式より小さい。符号の選択により変換のエネルギーを狭い周波数だけに集めずに済む。
:::

### 尺度が共役指数を強制する

$1\le p,q\le\infty$ で未知の評価 $\|\widehat f\|_q\le C\|f\|_p$ が全ての $f\in\mathcal S$ に成り立つとする。非零の $f$ を固定し、$f_a(x)=f(ax)$（$a>0$）を入れる。変数変換により $\|f_a\|_p=a^{-n/p}\|f\|_p$、$\|\widehat{f_a}\|_q=a^{-n+n/q}\|\widehat f\|_q$ である。従って $a^{n(1/p+1/q-1)}$ が全ての $a>0$ で同じ定数により抑えられる必要がある。$a\downarrow0$ と $a\uparrow\infty$ の両側を考えると $1/p+1/q=1$ が必要である。:ref[prop-four-hausdorff-young-range] は、その上で $p\le2$ も必要なことを示している。片方の指数を多少変えても同じ形の評価が成り立つとは限らないのは、$\R^n$ では拡大と縮小の両方が可能だからである。

:::example[有限測度の空間にある包含を使えない]
一次元で $L^1\to L^2$ の評価を試みる。$f_a(x)=a e^{-\pi (ax)^2}$ の質量は $1$ だが、$\widehat f_a=e^{-\pi \xi^2/a^2}$ の二乗ノルムは $\|\widehat f_a\|_2^2=a/\sqrt2$ である。$a\to\infty$ にすると $L^1$ ノルムは変わらず $L^2$ ノルムだけが増える。この例の空間では測度が無限であるため、$L^\infty$ の有界性から $L^2$ の有界性を導けない。
:::

### 周期の場合と最良定数の違い

円周では対応する不等式が $\|\widehat f\|_{\ell^{p'}}\le\|f\|_{L^p(dm)}$ となる。$L^1\to\ell^\infty$ は係数の積分の三角不等式、$L^2\to\ell^2$ は円周の Parseval によるから、同じ Riesz–Thorin 補間が使える。円周では $f(t)=ce^{ikt}$ を入れると $\|\widehat f\|_{\ell^{p'}}=|c|=\|f\|_p$ なので、定数 $1$ は最良である。

Euclid 空間で $1<p<2$ の最良定数は Beckner の定理によって
$(p^{1/p}/(p')^{1/p'})^{n/2}$
となる。鋭い上界には Gauss 測度に関する複素超縮小性の議論が必要であり、本章の二端点補間の範囲を超える。証明は :cite[grafakos-classical] の第 2 章の歴史的注記に挙げられた Beckner の原論文を参照。Gauss 入力でこの比が達成されることは本章の計算で確かめられる。実際、$f=e^{-\pi a|x|^2}$ なら $\|f\|_p=(pa)^{-n/(2p)}$、$\|\widehat f\|_{p'}=a^{-n/2}(p'/a)^{-n/(2p')}$ である。共役指数の関係により商から $a$ が消え、表示した定数に一致する。例えば $n=2,p=4/3$ で商は $(4/3)^{3/4}/4^{1/4}$ で、四乗すると $16/27$ である。この計算で最良定数の下界を得る。すべての入力で同じ値を上界にできることは別途必要である。

### 補間した変換と空間での切断

$1<p<2$ では $f$ が可積分でなくても、$f_R=f\,1_{\{|x|\le R,\ |f(x)|\le R\}}$ は $L^1\cap L^p$ に属し、$f_R\to f$ が $L^p$ で成り立つ。従って $\widehat f_R\to\widehat f$ が $L^{p'}$ である。これで定義されるのは $L^{p'}$ の同値類であり、各 $\xi$ での値を取り出すには別の条件が必要になる。例えば $p=3/2$ で $f(x)=(1+|x|)^{-5/6}$ とすれば $\|f\|_{3/2}^{3/2}=2\int_0^\infty(1+x)^{-5/4}dx=8$ であるが $\|f\|_1=\infty$ である。補間変換の $L^3$ ノルムは $8^{2/3}=4$ 以下であり、原点で絶対積分による値は定まらない。

逆変換にも同じ評価が成り立つ。$\check g(x)=\widehat g(-x)$ であり、反射はすべての $L^q$ ノルムを保つためである。従って $f\in L^2$ かつ $\widehat f\in L^p$（$1\le p\le2$）なら $f\in L^{p'}$ で $\|f\|_{p'}\le\|\widehat f\|_p$ である。切断で近似した $\widehat f$ の逆変換は $L^2$ で元の $f$ に近づき、同時に $L^{p'}$ で逆変換の拡張値に近づく。コンパクト集合上では両方が測度収束を含意するので、同じ関数を与える。この形は周波数情報から空間のノルムを制御する際に使う。

## 畳み込みと Fourier 変換

:::theorem[畳み込みと積]{#thm-four-euclid-convolution}
$f,g\in L^1$ なら $\widehat{f*g}=\widehat f\,\widehat g$ である。$f,g\in\mathcal S$ なら $\widehat{fg}=\widehat f*\widehat g$ も成り立つ。さらに $f\in L^1$、$g\in L^2$ なら第 1 式は $L^2$ の変換で成り立つ。
:::

:::proof
畳み込みの定義は :ref[def-meas-convolution] に従う。$u=x-y$ と変数変換し、二重積分の絶対値が $\|f\|_1\|g\|_1$ であることから Fubini を使う。Schwartz の積では $g$ の反転式を代入し、$\int|f|\int|\widehat g|<\infty$ で順序を交換すると $\int\widehat f(\xi-\eta)\widehat g(\eta)\,d\eta$ となる。$g\in L^2$ のときは $L^1\cap L^2$ の関数で近似する。:ref[thm-meas-young] が畳み込み側を、$\|\widehat f\|_\infty\le\|f\|_1$ が積の側を制御する。
:::

:::example[区間を二回平均する]
$f=1_{[-1/2,1/2]}$ とすると重なりの長さから $f*f(x)=(1-|x|)_+$ である。よって三角形の変換は $(\sin(\pi\xi)/(\pi\xi))^2$ となる。変換が可積分になったので、反転式を原点で使ってもその積分が $1$ となる。
:::

:::example[平均の分散が足し合わさる]
$G_s(x)=s^{-n/2}e^{-\pi|x|^2/s}$ とすれば $\widehat G_s=e^{-\pi s|\xi|^2}$ である。$s,t>0$ で変換を掛け、:ref[thm-four-euclid-inversion] の一意性を使うと $G_s*G_t=G_{s+t}$ となる。例えば一次元で $s=1,t=3$ なら原点の値は $1/2$ である。
:::

### $L^2$ 関数の積も変換できる

$f,g\in L^2$ の積は Cauchy–Schwarz により $L^1$ に入る。周波数側では $\widehat f*\widehat g$ も各点で絶対積分により定まり、
$|(\widehat f*\widehat g)(\xi)|\le\|\widehat f\|_2\|\widehat g\|_2$
である。よって Schwartz の積公式をこの場合に延長できる。実際 $f_m,g_m\in\mathcal S$ が $f,g$ に $L^2$ で近づくように取ると

$$
\|f_mg_m-fg\|_1\le\|f_m-f\|_2\|g_m\|_2+\|f\|_2\|g_m-g\|_2\longrightarrow0
$$

である。その変換は一様に $\widehat{fg}$ に近づく。周波数側でも畳み込み差の上限を同様の Cauchy–Schwarz で評価でき、Plancherel によって $\widehat f_m*\widehat g_m$ は $\widehat f*\widehat g$ に一様収束する。よって $\widehat{fg}=\widehat f*\widehat g$ が全点で成り立つ。例えば $\widehat f=1_{[-1,1]}$、$\widehat g=1_{[-1/2,1/2]}$ なら変換側の畳み込みは、$|\xi|\le1/2$ で $1$、$1/2<|\xi|<3/2$ で $3/2-|\xi|$、残りで $0$ である。積によって二つの帯域の半径が加わるが、その中での高さは積分区間の重なりに依存する。

### 畳み込みの台と平滑性

$f,g\in L^1$ がそれぞれコンパクト集合 $A,B$ に台を持つならば、$f*g$ は $A+B=\{a+b:a\in A,b\in B\}$ に台を持つ。実際 $x\notin A+B$ なら、$y\in B$ と $x-y\in A$ は同時に成り立たない。台の外でゼロとなる代表元を選んで畳み込みの積分を見れば結論を得る。さらに $g\in C_c^\infty$ のとき、任意の階の差商は $|f(y)|$ と $g$ の適切な導関数の上限で支配される。従って $\partial^\alpha(f*g)=f*\partial^\alpha g$ であり、畳み込みは滑らかになる。このとき周波数側で $\widehat f$ は単に有界でも、$\widehat g$ が急減少するから、その積は急減少する。実空間でも周波数空間でも、同じ正則化を異なる式で説明できる。

:::example[三回の区間平均]
$r=1_{[-1/2,1/2]}$ とし $r*r=(1-|x|)_+$ をもう一度 $r$ と畳み込む。$0\le x\le1/2$ の積分区間は $[x-1/2,x+1/2]$ であり、$0$ で二つに分けて $\int(1-|u|)du$ を計算すると $3/4-x^2$ となる。$1/2\le x\le3/2$ では $[x-1/2,1]$ 上の $1-u$ の積分で $(3/2-x)^2/2$ となる。外側はゼロであり、負の $x$ へは偶に延長する。したがって最終的な関数は $C^1$ であり、変換は $(\sin(\pi\xi)/(\pi\xi))^3$ である。三回の畳み込みで台の長さは $3$ になり、周波数の原点の値は質量 $1$ のままである。
:::

### 畳み込みの逆問題

$k\in L^1$ とし、測定値 $h=k*f$ から $f\in L^2$ を求める。$L^2$ の積公式により、変換を取ると $\widehat h=\widehat k\,\widehat f$ である。$|\widehat k|>0$ が a.e. で成り立てば割り算で $\widehat f$ が a.e. 一意に定まる。逆に $\{\widehat k=0\}$ の測度が正なら、その集合の有界な部分で正測度のものを取り、指示関数を変換とする $f\ne0$ を Plancherel で作れる。すると $k*f=0$ なので復元は一意でない。

核が一意性を保っても、復元の $L^2$ 安定性は別に調べる。例えば $k(x)=e^{-\pi x^2}$ とし、$\widehat f_j=1_{[j,j+1]}$ とする。$\|f_j\|_2=1$ だが $\|k*f_j\|_2\le e^{-\pi j^2}$ であり、測定値はゼロへ近づく。復元の作用素を変換側で考えると、$1/\widehat k$ の乗算が高周波を増幅することが分かる。一方、$f$ の周波数が $[-B,B]$ に限られると既知なら $|\widehat k|\ge e^{-\pi B^2}$ であるから、$\|f\|_2\le e^{\pi B^2}\|k*f\|_2$ が成り立つ。例えば $B=1$ の安定性定数は $e^\pi$ であり、入力の帯域が広がると定数も増大する。

## 例：Gauss 関数と球対称関数の Fourier 変換

球対称な関数の変換は、周波数ベクトルの向きに依存しない。この性質により角度の積分を先に行い、半径に関する一次元積分へ変えられる。

:::proposition[球対称変換の積分式]{#prop-four-radial-transform}
$n\ge2$、$f(x)=F(|x|)\in L^1$ とする。$s=|\xi|$ について

$$
\widehat f(\xi)=|S^{n-2}|\int_0^\infty F(r) r^{n-1}
\int_{-1}^1 e^{-2\pi irsu}(1-u^2)^{(n-3)/2}\,du\,dr.
$$

$n=3$ では $\widehat f(\xi)=4\pi\int_0^\infty F(r)r^2\,\sin(2\pi rs)/(2\pi rs)\,dr$ となる。
:::

:::proof
直交行列について :ref[prop-four-euclid-symmetries] を使えば $\xi=se_1$ に帰着する。$x=r\omega$ と置き、$\omega=(u,\sqrt{1-u^2}\eta)$、$\eta\in S^{n-2}$ と書く。$u$ 方向の微小長さに $(1-u^2)^{-1/2}$、残りの $n-2$ 方向に $\sqrt{1-u^2}$ が掛かるので、球面測度は $(1-u^2)^{(n-3)/2}\,du\,d\sigma(\eta)$ である。$n=2$ では $S^0$ の二点を数える。$n=3$ では $|S^1|=2\pi$ で、内側は $2\sin(2\pi rs)/(2\pi rs)$ である。$rs=0$ では極限値を取る。
:::

角度積分を展開すると Bessel 関数の表示が得られる。$\nu=n/2-1$ とし、

$$
J_\nu(z)=\sum_{k\ge0}\frac{(-1)^k(z/2)^{2k+\nu}}{k!\,\Gamma(k+\nu+1)}
\qquad(z>0)
$$

と定める。ここで $\Gamma(a)=\int_0^\infty e^{-t}t^{a-1}dt$ である。角度積分で奇数冪は消え、偶数冪の積分は $v=u^2$ による置換で

$$
\int_{-1}^1 u^{2k}(1-u^2)^{\nu-1/2}du
=\frac{\Gamma(k+1/2)\Gamma(\nu+1/2)}{\Gamma(k+\nu+1)}
$$

となる。右辺は二つの Gamma 積分を和と比の変数に移すことで得られる。$\Gamma(k+1/2)/(2k)!=\sqrt{\pi}/(4^k k!)$ と球面積 $|S^{n-2}|=2\pi^{(n-1)/2}/\Gamma((n-1)/2)$ を使うと

$$
\widehat f(\xi)=2\pi\,s^{-\nu}\int_0^\infty F(r)J_\nu(2\pi rs)r^{\nu+1}\,dr\quad(s>0).
$$

球面積は $\int_{\R^{n-1}}e^{-|x|^2}dx=\pi^{(n-1)/2}$ を極座標表示すれば従う。級数の交換は角度積分の有限区間で各 $rs$ を固定して行った。

:::example[三次元の単位球]
$F=1_{[0,1]}$ とする。$b=2\pi s$ と置き、$\int_0^1r\sin(br)\,dr=(\sin b-b\cos b)/b^2$ を部分積分で求めると

$$
\widehat f(\xi)=4\pi\,\frac{\sin b-b\cos b}{b^3}.
$$

原点では分子が $b^3/3+O(b^5)$ で、極限は球の体積 $4\pi/3$ である。
:::

:::example[異方的な Gauss 関数]
$f(x_1,x_2)=e^{-\pi(4x_1^2+x_2^2)}$ なら $\widehat f=(1/2)e^{-\pi(\xi_1^2/4+\xi_2^2)}$ である。原点の値は $\int f=1/2$ である。第 1 軸の空間幅が半分になったため、変換側の幅は 2 倍になる。
:::

### 反転を使う場合の積分順序

$f*g$ の積分と $f$ の逆変換は、絶対収束がある場合に同じように見える。ところが矩形波では変換が $L^1$ に属さず、反転の積分を通常の絶対積分とは読めない。そこで周波数で対称に切断し、点を固定した極限を調べる。

:::proposition[対称切断に対する Dini の判定法]{#prop-four-euclid-dini-inversion}
$f\in L^1(\R)$、$x\in\R$、$s\in\C$ とする。ある $\delta>0$ に対して $\int_0^\delta |f(x+t)+f(x-t)-2s|\,dt/t<\infty$ なら、$\int_{-R}^R\widehat f(\xi)e^{2\pi ix\xi}d\xi\to s$ である。
:::

:::proof
有限区間で Fubini を使うと、逆積分は $\int_\R f(x-t)\sin(2\pi Rt)/(\pi t)\,dt$ となる。正負の $t$ を合わせ、$F_x(t)=f(x+t)+f(x-t)$ と置く。$0<t<\delta$ の振幅から $2s$ を引けば、Dini 条件により $(F_x(t)-2s)/t$ は可積分である。$t\ge\delta$ では $|F_x(t)|/t\le |F_x(t)|/\delta$ なので可積分である。従って両区間で振幅をゼロ延長して Riemann–Lebesgue の補題を適用できる。残る項は $(2s/\pi)\int_0^\delta\sin(2\pi Rt)\,dt/t$ であり、$v=2\pi Rt$ と置くと :ref[prop-four-sine-integral-limit] により $s$ に近づく。
:::

:::proposition[区間指示関数の条件収束による反転]{#prop-four-box-symmetric-inversion}
$f=1_{[-a,a]}$（$a>0$）なら
$$
\lim_{R\to\infty}\int_{-R}^R\widehat f(\xi)e^{2\pi ix\xi}d\xi
=\begin{cases}
1&(|x|<a),\\
1/2&(|x|=a),\\
0&(|x|>a).
\end{cases}
$$
:::

:::proof
$\widehat f$ は実偶関数なので指数の実部だけでよい。積和公式で有限積分は
$$
\frac1\pi\int_0^R\frac{\sin(2\pi(a+x)\xi)+\sin(2\pi(a-x)\xi)}{\xi}\,d\xi
$$
である。各周波数因子が正なら :ref[prop-four-sine-integral-limit] で極限は $\pi/2$、負なら $-\pi/2$、ゼロなら積分はゼロである。$a+x,a-x$ の符号を三つの場合で数えれば結果が従う。
:::

:::example[二軸で端点が重なる]
$f(x_1,x_2)=1_{[-1,1]}(x_1)1_{[-1,1]}(x_2)$ の逆積分を $[-R,R]^2$ で切ると、積は一次元の二つの有限積分に分かれる。角 $(1,1)$ での極限は $(1/2)^2=1/4$、辺 $(1,0)$ では $1/2$ である。対称な正方形の切断を使っており、ここから球による切断をただちに同一視することはできない。
:::

### 球による周波数切断

三次元で $R>0$ とし、$K_R(x)=\int_{|\xi|\le R}e^{2\pi ix\cdot\xi}\,d\xi$ と置く。球の変換式を半径 $R$ に拡大すると、$r=|x|>0$ に対して

$$
K_R(x)=\frac{\sin(2\pi Rr)-2\pi Rr\cos(2\pi Rr)}{2\pi^2r^3}
$$

である。$f\in L^1(\R^3)$ なら、有限領域での Fubini により $\int_{|\xi|\le R}\widehat f(\xi)e^{2\pi ix\cdot\xi}\,d\xi=(f*K_R)(x)$ である。球面平均 $M_x(r)=(4\pi)^{-1}\int_{S^2}f(x-r\omega)\,d\sigma(\omega)$ を使うと、この値は

$$
\frac2\pi\int_0^\infty M_x(r)
\left(\frac{\sin(2\pi Rr)}r-2\pi R\cos(2\pi Rr)\right)\,dr
$$

となる。例えば $f\in C_c^1$ なら $M_x$ はコンパクト台を持つ $C^1$ 関数であり、$M_x(0)=f(x)$ である。余弦の項を部分積分すると、境界項は消え、$M_x'(r)\sin(2\pi Rr)$ の積分は Riemann–Lebesgue の補題によりゼロへ向かう。正弦の項では $M_x(r)-f(x)1_{[0,1]}(r)$ を $r$ で割ったものが可積分である。Riemann–Lebesgue の補題と :ref[prop-four-sine-integral-limit] を使えば極限は $f(x)$ となる。

:::example[原点の周囲で一定でも球切断が振動する例]
$f=1_{\{|x|\le1\}}$ の原点での球面平均は $M_0(r)=1_{[0,1]}(r)$ である。従って球切断した逆積分は
$$
\frac2\pi\left(\int_0^{2\pi R}\frac{\sin v}{v}\,dv-\sin(2\pi R)\right)
$$
となる。第一項は $1$ に近づくが、第二項は振動する。例えば $R$ を整数で取ると極限は $1$、整数に $1/4$ を足した列では $1-2/\pi$ である。原点から離れた半径 $1$ の跳躍が振動を残す点で、一次元の Dini 判定とは条件が異なる。
:::

### 正定値二次形式の Gauss 積分

各座標の幅が異なる場合も、直交座標へ移れば同じ計算になる。$A$ を実対称正定値行列とし $g_A(x)=e^{-\pi x^{\mathsf T}Ax}$ と置く。直交行列で $A$ を対角化し、正の固有値を $\lambda_1,\ldots,\lambda_n$ とする。その座標で積分を分け、一次元 Gauss 関数の拡大縮小を使うと

$$
\widehat g_A(\xi)=(\det A)^{-1/2}e^{-\pi\xi^{\mathsf T}A^{-1}\xi}
$$

である。行列式の平方根は正のものを取る。$x$ の変数変換で Jacobian が現れ、周波数の指数には逆行列が現れる。例えば $A=\begin{pmatrix}2&1\\1&2\end{pmatrix}$ なら $\det A=3$、$A^{-1}=\frac13\begin{pmatrix}2&-1\\-1&2\end{pmatrix}$ で、$\xi=(1,1)$ での値は $3^{-1/2}e^{-2\pi/3}$ である。第 $(1,1)$ 成分の逆数を個別に取ると、非対角項から来る相関を失う。

### 二次の位相を含む Gauss 関数

幅を複素数へ移した計算は、実部の正値性を保てば積分で行える。$c=a+ib$、$a>0$、$b\in\R$ とし、$c^{1/2}$ は実部が正の平方根を取る。すると一次元で

$$
\int_{\R}e^{-\pi c x^2}e^{-2\pi ix\xi}\,dx
=c^{-1/2}e^{-\pi \xi^2/c}
$$

である。実際、$G(\xi)$ を左辺として $x$ で部分積分すると、$\partial_x e^{-\pi c x^2}=-2\pi c x e^{-\pi c x^2}$ から $G'=-(2\pi/c)\xi G$ を得る。境界項は実部 $a>0$ により消える。初期値 $I(b)=G(0)$ は $b$ で微分でき、その積分は $\pi x^2 e^{-\pi a x^2}$ により支配される。$\int(xe^{-\pi c x^2})'\,dx=0$ を使うと

$$
I'(b)=-\pi i\int x^2e^{-\pi c x^2}\,dx=-\frac{i}{2c}I(b).
$$

$c^{-1/2}$ は同じ微分方程式を満たし、$b=0$ でどちらも $a^{-1/2}$ なので一致する。$G$ の微分方程式を解けば求める式を得る。

:::example[絶対値と位相を別々に追う]
$c=1+i$ なら $c^{-1/2}=2^{-1/4}e^{-i\pi/8}$、$1/c=(1-i)/2$ である。従って変換は
$2^{-1/4} e^{-\pi \xi^2/2} e^{i(\pi \xi^2/2-\pi/8)}$
となる。入力の絶対値は $e^{-\pi x^2}$ のままだが、変換の絶対値は $e^{-\pi\xi^2/2}$ の幅を持つ。二乗積分は空間側で $1/\sqrt2$、周波数側でも $2^{-1/2}\int e^{-\pi \xi^2}\,d\xi=1/\sqrt2$ であり、Plancherel に合う。
:::

### 球と球面を変換する

半径 $a$ の球の指示関数に球対称変換を適用する。$s=|\xi|>0$、$\nu=n/2-1$ とする。Bessel の級数から、$\kappa>0$ について
$\frac{d}{dr}(r^{\nu+1}J_{\nu+1}(\kappa r))=\kappa r^{\nu+1}J_\nu(\kappa r)$
が従う。実際、左辺を項別微分すると、各項のべきに由来する $2k+2\nu+2$ が $\Gamma(k+\nu+2)$ の因子 $k+\nu+1$ を打ち消し、右の冪級数になる。級数と導関数の収束は $r$ の閉有界区間で項比から従う。従って $n\ge2$ で

$$
\widehat{1_{\{|x|\le a\}}}(\xi)=(a/s)^{n/2}J_{n/2}(2\pi as).
$$

$r=0$ で原始関数 $r^{\nu+1}J_{\nu+1}(2\pi rs)$ はゼロになるため、この式では $r=a$ の項だけが残る。$J_{n/2}(z)\sim(z/2)^{n/2}/\Gamma(n/2+1)$ より、$s\downarrow0$ で球の体積 $\pi^{n/2}a^n/\Gamma(n/2+1)$ に近づく。

同じ半径 $a$ の球面測度 $\sigma_a$ の変換では $r$ の積分を取らない。球面の変数を $x=a\omega$ と置くと $d\sigma_a=a^{n-1}d\sigma_1$ であり、角度積分から

$$
\widehat{\sigma_a}(\xi)=2\pi a^{n/2}s^{1-n/2}J_{n/2-1}(2\pi as)
$$

を得る。球の変換を $a$ で微分しても同じ式になる。実際、球の積分を半径と角度に分けた後、上端の微分が球面の積分を与える。球面の質量は $n\pi^{n/2}a^{n-1}/\Gamma(n/2+1)$ で、球の体積とは $a$ の指数が一つ違う。

:::example[円板の中心と球面変換の零点]
$n=2,a=1$ の円板では $\widehat{1_{\{|x|\le1\}}}(\xi)=J_1(2\pi s)/s$ である。$J_1(z)=z/2-z^3/16+O(z^5)$ より原点近くで $\pi-(\pi^3/2)s^2+O(s^4)$ となる。$n=3,a=1$ の球面では角度積分を直接行うと $4\pi\sin(2\pi s)/(2\pi s)$ であり、$s=1/2$ が最初の正の零点である。同じ半径の三次元の球では $s=1/2$ における分子が $\sin\pi-\pi\cos\pi=\pi$ なので、そこで変換はゼロにならない。
:::

### Poisson 核を Gauss 積分から求める

Gauss 核を幅の異なるものの重み付き和として積分すると、周波数で $e^{-2\pi y|\xi|}$ となる核が得られる。この計算では一次元と高次元の正規化を同じ式から求められる。

:::lemma[平方根の指数の積分表示]{#lem-four-gaussian-subordination}
$a>0,b\ge0$ について

$$
e^{-a\sqrt b}=\frac{a}{2\sqrt\pi}\int_0^\infty t^{-3/2}e^{-a^2/(4t)}e^{-bt}\,dt .
$$
:::

:::proof
$J(c)=\int_0^\infty e^{-u^2-c^2/u^2}\,du$ とする。$c>0$ では微分の被積分関数を原点と無限遠で積分可能な関数で支配できる。$J'(c)=-2c\int_0^\infty u^{-2}e^{-u^2-c^2/u^2}du$ であり、$v=c/u$ と置くと右辺は $-2J(c)$ になる。一方 $c\downarrow0$ で優収束より $J(c)\to\sqrt\pi/2$ である。従って $J(c)=(\sqrt\pi/2)e^{-2c}$ である。求める積分で $u=a/(2\sqrt t)$ と置けば $(2/\sqrt\pi)J(a\sqrt b/2)$ となる。
:::

:::proposition[高次元の Poisson 核の変換]{#prop-four-euclid-poisson-transform}
$y>0$ とし

$$
P_y(x)=\frac{\Gamma((n+1)/2)}{\pi^{(n+1)/2}}\frac{y}{(y^2+|x|^2)^{(n+1)/2}}
$$

と置く。$P_y$ は質量 $1$ で $\widehat P_y(\xi)=e^{-2\pi y|\xi|}$ である。
:::

:::proof
:ref[lem-four-gaussian-subordination] に $a=y,b=4\pi^2|\xi|^2$ を使う。$e^{-4\pi^2t|\xi|^2}$ の逆変換は Gauss の尺度変換で $G_t(x)=(4\pi t)^{-n/2}e^{-|x|^2/(4t)}$ である。よって $G_t$ を $\frac{y}{2\sqrt\pi}t^{-3/2}e^{-y^2/(4t)}dt$ で平均する。この重みの積分は補題の $b=0$ から $1$ で、各 $G_t$ の質量も $1$ なので、Tonelli により核の平均は質量 $1$ の $L^1$ 関数になる。変換と $t$ 積分を交換して $e^{-2\pi y|\xi|}$ を得る。直接 $x$ 側の式を求めるには

$$
\frac{y}{2\sqrt\pi(4\pi)^{n/2}}
\int_0^\infty t^{-(n+3)/2}e^{-(y^2+|x|^2)/(4t)}\,dt
$$

に $u=(y^2+|x|^2)/(4t)$ を入れる。残る Gamma 積分により表示式を得る。
:::

:::example[変換と質量を二次元で確認する]
$n=2$ では $\Gamma(3/2)=\sqrt\pi/2$ なので $P_y(x)=y/(2\pi(y^2+|x|^2)^{3/2})$ である。極座標で質量は $y\int_0^\infty r(y^2+r^2)^{-3/2}dr=1$ となる。例えば $y=2$ で原点の値は $1/(8\pi)$ である。$P_y*P_z$ の変換は $e^{-2\pi(y+z)|\xi|}$ だから、一意性により $P_y*P_z=P_{y+z}$ である。この半群では $y$ がそのまま加わる。先の Gauss による正則化で $\eps$ を用いた場合は、畳み込みで $\eps^2$ が加わった。
:::

## 演習問題

::::exercise{level=basic}
$1_{[0,2]}$ の Fourier 変換を求めよ。
:::solution
$\xi\ne0$ では指数を積分すると $(1-e^{-4\pi i\xi})/(2\pi i\xi)$ になる。原点では長さ $2$ である。$[-1,1]$ を右に $1$ 移した結果として $e^{-2\pi i\xi}\sin(2\pi\xi)/(\pi\xi)$ とも書ける。
:::
::::

::::exercise{level=basic}
$e^{-\pi(x-2)^2}e^{4\pi ix}$ の変換を求めよ。
:::solution
先に Gauss 関数を $2$ だけ移すと $e^{-4\pi i\xi}e^{-\pi\xi^2}$ となる。さらに周波数を $2$ だけ移すので $e^{-4\pi i(\xi-2)}e^{-\pi(\xi-2)^2}$ である。
:::
::::

::::exercise{level=basic}
$g=e^{-\pi x^2}$ として $\widehat{xg}$ を求めよ。
:::solution
:ref[prop-four-schwartz-operations] から $\widehat{xg}=-(2\pi i)^{-1}(\widehat g)'$ である。$\widehat g=g$ を微分して $-i\xi e^{-\pi\xi^2}$ を得る。
:::
::::

::::exercise{level=basic}
$f=1_{[-2,2]}$ の変換の二乗積分を求めよ。
:::solution
$f^2=f$ で積分は $4$ である。:ref[thm-four-plancherel] により $\int|\widehat f|^2=4$ である。
:::
::::

::::exercise{level=advanced}
$f\in L^1(\R)$ とし、$xf(x)\in L^1$ を仮定する。$\widehat f$ が $C^1$ で $(\widehat f)'=-2\pi i\widehat{xf}$ となることを証明せよ。
:::solution
差商の被積分関数は $f(x)e^{-2\pi ix\xi}(e^{-2\pi ixh}-1)/h$ である。$|e^{iu}-1|\le|u|$ より $2\pi|xf(x)|$ で支配される。優収束で極限を積分に入れれば微分式になり、:ref[thm-four-euclid-rl] で導関数の連続性が従う。
:::
::::

::::exercise{level=advanced}
$f\in L^1$ について $\widehat f(\xi)=0$ が原点の近傍で成り立つとする。$x^\alpha f\in L^1$ が全ての多重指数で成立するなら全てのモーメントが消えることを示せ。
:::solution
前問の支配を各方向に繰り返すと $\partial_\xi^\alpha\widehat f(0)=(-2\pi i)^{|\alpha|}\int x^\alpha f(x)dx$ になる。原点の近傍で $\widehat f$ が恒等的にゼロなので左辺はゼロであり、非零定数で割れる。
:::
::::

::::exercise{level=advanced}
$\int_{\R}(\sin(\pi\xi)/(\pi\xi))^4\,d\xi$ を求めよ。
:::solution
$f=1_{[-1/2,1/2]}$ の畳み込みの変換は $(\sin(\pi\xi)/(\pi\xi))^2$ である。畳み込みは $(1-|x|)_+$ なので、その二乗積分は $2\int_0^1(1-x)^2dx=2/3$ である。Plancherel で求める積分と一致する。
:::
::::

::::exercise{level=advanced}
$g\in\mathcal S(\R)$ とする。$\widehat{\widehat g\,*\,\widehat g}(x)=g(-x)^2$ を証明せよ。
:::solution
:ref[thm-four-euclid-convolution] で左辺は $\mathcal F^2g(x)\mathcal F^2g(x)$ となる。:ref[thm-four-euclid-inversion] によって各因子は $g(-x)$ に等しい。
:::
::::

::::exercise{level=basic}
$f(x_1,x_2)=e^{-\pi(3x_1^2+2x_1x_2+3x_2^2)}$ とする。$\widehat f(1,-1)$ を求めよ。
:::solution
二次形式の行列は $A=\begin{pmatrix}3&1\\1&3\end{pmatrix}$ である。固有値は $4,2$ で正定値、行列式は $8$ である。$(1,-1)$ は固有値 $2$ の固有ベクトルなので、逆行列の二次形式は $(1^2+(-1)^2)/2=1$ である。よって $\widehat f(1,-1)=e^{-\pi}/\sqrt8$ である。
:::
::::

::::exercise{level=advanced}
一次元で $\widehat f=1_{[-B,B]}$、$B>0$ とする。$f\in L^2$ の弱微分が $L^2$ に入ることを示し、そのノルムの二乗を求めよ。
:::solution
$|\xi|\le B$ に台があるので $\xi\widehat f\in L^2$ である。:ref[prop-four-l2-weak-derivative] から弱微分が存在し、その変換は $2\pi i\xi\,1_{[-B,B]}$ である。Plancherel により $\|f'\|_2^2=4\pi^2\int_{-B}^B\xi^2d\xi=8\pi^2B^3/3$ となる。一方 $f$ は有限帯域の逆積分で滑らかであり、その式を一回微分しても同じ変換の逆積分になる。
:::
::::

## 確認テスト

::::quiz{answer="4" tol="1e-6"}
$\widehat{1_{[-2,2]}}(0)$ を求めよ。
:::explanation
指数因子は $1$ なので区間の長さ $2-(-2)=4$ である。
:::
::::

::::quiz
$L^1(\R^n)$ 関数の変換について必ず成り立つものはどれか。

- [x] 連続で無限遠でゼロに近づく。
- [ ] 可積分である。
- [ ] 全点で微分可能である。

:::explanation
連続性と無限遠での消滅は :ref[thm-four-euclid-rl] に従う。区間指示関数の変換は正弦を $\xi$ で割った形で、各半周期の絶対積分を下から調和級数で抑えられるため可積分でない。可積分な Poisson 核の変換 $e^{-2\pi y|\xi|}$ は原点で微分できない。
:::
::::

::::quiz{answer="4" tol="1e-6"}
Hausdorff–Young の不等式で $p=4/3$ のとき $p'$ を求めよ。
:::explanation
$1/p+1/p'=1$ であり、$1/p=3/4$ だから $1/p'=1/4$、従って $p'=4$ である。
:::
::::

::::quiz
$f,g\in\mathcal S$ のとき $\mathcal F(fg)$ はどれか。

- [ ] $\widehat f\,\widehat g$
- [x] $\widehat f*\widehat g$
- [ ] $2\pi\,\widehat f*\widehat g$

:::explanation
$g(x)=\int\widehat g(\eta)e^{2\pi ix\cdot\eta}d\eta$ を積の変換に代入する。絶対可積分性で積分を交換すると、内側の積分が $\widehat f(\xi-\eta)$ になり、結果は $\widehat f*\widehat g$ である。
:::
::::

::::quiz{answer="0.6666666666666667" tol="1e-6"}
一次元で $\widehat f=1_{[-1,1]}$ とする。線形な切断 $m_1(\xi)=(1-|\xi|)_+$ を使った復元の二乗誤差を求めよ。
:::explanation
$[-1,1]$ 上で $(1-m_1)\widehat f=|\xi|$ であり、台の外ではゼロである。Plancherel により誤差の二乗は $\int_{-1}^1|\xi|^2d\xi=2/3$ である。
:::
::::

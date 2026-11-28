---
title: Laplace 方程式と調和関数
description: Laplace 方程式の基本解、平均値、境界積分による表示、Dirichlet エネルギーを導く。
notation:
  - symbol: '\Delta'
    meaning: 各 Euclid 座標についての二階微分の和
  - symbol: '\sigma_{n-1}'
    meaning: 単位球面の面積
  - symbol: '\Gamma'
    meaning: 負の Laplace 作用素の基本解
---

## この章の目標

- 基本解の定数を球面上の流束から決める。
- 球面と球の平均値公式から最大値原理を導く。
- 導関数評価から Liouville の定理と解析性を説明する。
- 球と半空間の Green 関数を構成して Poisson 積分を求める。
- 境界値問題を Dirichlet エネルギーの最小化に結び付ける。

## 基本解と Newton ポテンシャル

:::definition[調和関数]{#def-pde-harmonic}
開集合 $\Omega\subset\R^n$ 上の $C^2$ 関数 $u$ が $\Delta u=\sum_{j=1}^n u_{x_jx_j}=0$ を満たすとき :term[調和関数]{en="harmonic function" yomi="ちょうわかんすう"} という。$-\Delta u=f$ を :term[Poisson 方程式]{en="Poisson equation" yomi="ぽあそんほうていしき"} という。
:::

点に集中した入力を表すため、試験関数に $\delta_0(\phi)=\phi(0)$ を作用させる。[超関数](/math/fourier-analysis/distributions/) の意味で $-\Delta\Gamma=\delta_0$ を満たす $\Gamma$ を :term[基本解]{en="fundamental solution" yomi="きほんかい"} と呼ぶ。符号は Poisson 方程式の $-\Delta$ に対応させる。

:::proposition[基本解の規格化]{#prop-pde-laplace-fundamental}
$n\ge2$ で $\sigma_{n-1}=|\mathbb S^{n-1}|$ とする。基本解は

$$
\Gamma(x)=\begin{cases}
\displaystyle\frac{|x|^{2-n}}{(n-2)\sigma_{n-1}} & n\ge3,\\
-\dfrac{1}{2\pi}\log |x| & n=2 .
\end{cases}
\label{eq-pde-laplace-fundamental}
$$

:::

:::proof
$r=|x|>0$ で $u(x)=h(r)$ とすると、$\partial_i u=h'(r)x_i/r$ をもう一度微分して足せば $\Delta u=h''+(n-1)h'/r=r^{1-n}(r^{n-1}h')'$ である。よって原点以外で調和であり、$-\int_{\partial B_r}\partial_r\Gamma\,\dd S=1$ とするには $\Gamma'=-1/(\sigma_{n-1}r^{n-1})$ と選べばよい。これを積分すると表示を得る。原点で局所可積分であり、半径 $\eps$ の穴の外で Green の公式を用いると

$$
\int_{|x|>\eps}\Gamma(-\Delta\phi)\,\dd x
=\int_{\partial B_\eps}\big(\Gamma\partial_r\phi-\phi\partial_r\Gamma\big)\,\dd S .
$$

最初の項は $O(\eps)$、二次元では $O(\eps|\log\eps|)$ で消え、次の項は $\phi$ の球面平均であるから $\phi(0)$ に収束する。
:::

$f\in C_c^2(\R^n)$ なら $Nf=\Gamma*f$ を :term[Newton ポテンシャル]{en="Newtonian potential" yomi="にゅーとんぽてんしゃる"} という。$Nf(x)=\int\Gamma(z)f(x-z)\,\dd z$ と書くと、$x$ のコンパクト集合に対し $z$ の積分領域は有界に取れるので二階微分を $f$ に作用させられる。$-\Delta_xNf(x)=f(x)$ は :ref[prop-pde-laplace-fundamental] を $z\mapsto f(x-z)$ に適用すると従う。

:::example[三次元の一様球の外側]
原点を中心とする半径 $1$ の球内に一定密度 $1$ を置く。球対称のポテンシャルは外側で $A/r$ と書ける。半径 $r>1$ を通る負の流束が全量 $4\pi/3$ であるため $4\pi A=4\pi/3$、従って $Nf=1/(3r)$ である。内側で微分方程式 $-(r^2u')'=r^2$ と原点での有限性から $u'= -r/3$ となり、$r=1$ で接続すると $u=1/2-r^2/6$ になる。値と法線微分がつながるので境界に面密度はない。
:::

:::example[一次元の基本解]
$n=1$ では $\Gamma=-|x|/2$ である。傾きは左の $1/2$ から右の $-1/2$ に変わるので、$-\Gamma''=\delta_0$ となる。直接には二つの半直線で二度部分積分すると境界項の和が $\phi(0)$ である。
:::

### 有界な体積密度と遠方展開

有界でコンパクトな台を持つ $f$ にも $Nf=\Gamma*f$ を定義できる。$\Gamma$ と $D\Gamma$ は原点で局所可積分なので、固定した有界領域で積分の特異性を半径 $\eps$ の球で切ると、捨てた部分は $\norm f_\infty\int_{B_\eps}(|\Gamma|+|D\Gamma|)$ 以下である。これは一様に零へ収束する。切った核では連続性が得られるため、$Nf$ と $D\Gamma*f$ は連続である。任意の試験関数に対して Fubini と基本解の等式を使うと、$D(Nf)=D\Gamma*f$、$-\Delta Nf=f$ が分布の意味で成立する。従って $Nf$ は $C^1$ である。密度が段差を持つ球でも一階微分は連続するが、$-\Delta Nf$ には密度そのものの段差が残る。

$n\ge3$、$\operatorname{supp} f\subset B_a$ の場合、$|x|>2a$ なら $|x-sy|\ge|x|/2$（$0\le s\le1$）である。Taylor 展開と $|D^2\Gamma(z)|\le C_n|z|^{-n}$ から
$$
Nf(x)=M\Gamma(x)-D\Gamma(x)\cdot p+O\left(|x|^{-n}\int|y|^2|f(y)|\dd y\right),
\qquad M=\int f,\quad p=\int yf.
$$
従って総量が零のときは一次モーメントが最初の候補になる。密度に符号がある場合、単に $\int|f|$ で押さえた上界ではこの相殺は見えない。

:::example[距離が大きい正負二つの点源]
三次元で $e=(1,0,0)$ とし、点源 $\delta_e-\delta_{-e}$ を考える。基本解の平行移動から解は $(|x-e|^{-1}-|x+e|^{-1})/(4\pi)$ である。$x=Re$、$R>1$ では $1/(2\pi(R^2-1))$ となり、$R^2 u\to1/(2\pi)$ である。$M=0$、$p=2e$ として遠方展開の $-D\Gamma(Re)\cdot p$ を計算しても $1/(2\pi R^2)$ となる。垂直二等分面 $x_1=0$ では距離が等しく、すべての点で相殺する。
:::

:::example[一定でない密度を球殻の中と外でつなぐ]
三次元の単位球に密度 $f(r)=6+20r^2$ を入れ、無限遠で零になる Newton ポテンシャルを求める。まず球内の式を作る。

放射方向だけに依存する解を考える。半径 $r$ の球面積が $\sigma_{n-1}r^{n-1}$ であることから、発散を計算すると $-(r^{n-1}u_r)^{\prime}=r^{n-1}f$ となる。これを一度積分し、中心を含む場合には中心に余分な点源がないよう流束の極限を零にする。環状領域では内側の境界条件で積分定数を決める。

$$
\begin{aligned}
u&=- r^{4} - r^{2} + 8\\
u_r&=\frac{d}{d r} \left(- r^{4} - r^{2} + 8\right)=- 2 r \left(2 r^{2} + 1\right)\\
r^{n-1}u_r&=- 2 r^{3} \left(2 r^{2} + 1\right)\\
-r^{n-1}f&=\frac{d}{d r} \left(- 2 r^{3} \left(2 r^{2} + 1\right)\right)=- 2 r^{2} \left(10 r^{2} + 3\right)\\
-\Delta u&=2 \left(10 r^{2} + 3\right)=20 r^{2} + 6
\end{aligned}
$$

次に値と流束を境界で読む。外側では法線は半径が増える向きである。内側にも境界がある場合は法線が反対向きになる。源項の体積積分を角度方向に積分すると面積係数が現れ、流束との差を確かめられる。

$$
\begin{aligned}
u&=8,\qquad r=0\\
u&=6,\qquad r=1\\
u_r&=-6,\qquad r=1\\
\sigma_{n-1}&=4 \pi\\
\text{全源量}&=\int\limits_{0}^{1} 4 \pi r^{2} \left(20 r^{2} + 6\right)\, dr=24 \pi\\
\text{負の全流束}&=24 \pi
\end{aligned}
$$

観測半径で値と半径微分を計算する。勾配の二乗を用いてエネルギーも計算できる。観測点で勾配が小さい場合でも、球面積の因子があるため、そこを通る全流束が小さいとは限らない。

$$
\begin{aligned}
u&=\frac{123}{16},\qquad r=\frac{1}{2}\\
u_r&=- \frac{3}{2},\qquad r=\frac{1}{2}\\
\sigma_{n-1}r^{n-1}u_r&=- \frac{3 \pi}{2},\qquad r=\frac{1}{2}\\
\frac12\int |Du|^2&=\int\limits_{0}^{1} 2 \pi r^{2} \left(- 4 r^{3} - 2 r\right)^{2}\, dr=\frac{3064 \pi}{315}
\end{aligned}
$$

外部では源がないので $u=B/r$ である。球内の積分定数は、外部の $6/r$ と境界値を一致させるように $8$ に選んだ。値は $6$、半径微分は両側とも $-6$ であり、継ぎ目の面測度は生じない。
:::

## 平均値の性質

球面平均を $M(r)=\sigma_{n-1}^{-1}\int_{\mathbb S^{n-1}}u(x+r\theta)\,\dd S_\theta$ とおく。球面の半径を増やして平均がどう変わるかを見ると、変化率が球内の源の総量に比例する。

:::theorem[平均値公式と逆]{#thm-pde-harmonic-mean}
$\overline B_r(x)\subset\Omega$ を満たす球について、調和関数は

$$
u(x)=\frac{1}{|\partial B_r|}\int_{\partial B_r(x)}u\,\dd S
=\frac{1}{|B_r|}\int_{B_r(x)}u\,\dd y
\label{eq-pde-harmonic-mean}
$$

を満たす。逆に連続関数が全ての内球の球体平均に等しければ、調和である。
:::

:::proof
球面平均を微分し :ref[thm-anl-divergence] を使うと $M'(r)=|\partial B_r|^{-1}\int_{B_r(x)}\Delta u\,\dd y=0$ である。$r\to0$ で $M(r)\to u(x)$ となり、半径で積分すれば球体平均も得られる。逆を示す。球体積で両辺を掛けて半径について微分すると、連続性により球面平均も同じになる。$\rho$ を単位球に台を持つ $C_c^\infty$ の非負な放射対称関数で $\int\rho=1$ とし、$\rho_r(y)=r^{-n}\rho(y/r)$ とおく。球面平均を半径で重ねると $u=u*\rho_r$ が境界から距離 $r$ より離れた場所で成り立つ。よって $u$ は滑らかである。球面上で Taylor 展開すると一次項と非対角の二次項は対称性で消える。また $\theta_i^2$ の球面平均は全て等しく、和は $1$ なので $1/n$ である。従って $M(r)=u(x)+r^2\Delta u(x)/(2n)+o(r^2)$ であり、平均値の等式から $\Delta u=0$ を得る。
:::

:::example[二次式の平均]
$u=x_1^2-x_2^2$ について中心 $(2,1)$ の円周の平均は $3$ である。周上で二乗を展開すると一次項は平均 $0$ で、二方向の二次の寄与 $r^2/2$ が消し合う。逆に $v=x_1^2+x_2^2$ の同じ円周平均は $5+r^2$ であり、$\Delta v=4$ を反映する。
:::

:::example[Hölder の不等式による局所値評価]
$u$ が $B_r(x)$ 上で調和なら、$1\le p<\infty$ について平均値公式と Hölder より $|u(x)|\le |B_r|^{-1/p}\norm u_{L^p(B_r(x))}$ である。例えば二次元の半径 $2$ の円盤で $\norm u_2=4$ なら中心値の絶対値は $2/\sqrt\pi$ 以下である。
:::

### 球面平均から得る Poisson 型の誤差

$\Delta u=f$、$f$ は連続とする。半径 $r$ の球面平均を $M(r)$ とし、同じ球の体積平均を $A_r f$ とすると、発散定理から $M'(r)=rA_r f/n$ となる。従って
$$
M(r)-u(x)=\frac1n\int_0^r s A_s f\dd s .
$$
特に $m\le f\le M_0$ なら $mr^2/(2n)\le M(r)-u(x)\le M_0r^2/(2n)$ である。$f$ の符号が一定なら平均の増減も決まる。$\Delta u$ の符号が反対になる $-\Delta u=f$ では不等号をそのまま写せない。

:::example[二次式と四次式の平均の差]
二次元で $u=x_1^2+2x_2^2$ なら $\Delta u=6$ である。原点での円平均は $3r^2/2$、円板平均は $3r^2/4$ で、後者は前者に面素 $2s\dd s/r^2$ を掛けて積分した値である。次に $u=|x|^4$ なら $\Delta u=16|x|^2$ であり、$A_s f=8s^2$ である。従って平均差は $\int_0^r4s^3\dd s=r^4$。点での Laplacian が零でも、すべての半径で平均が値に一致するとは限らない。
:::

:::example[円周の平均で消える異なる振動数]
半径 $2$ の円板で $g=4+3\cos\theta-2\sin(3\theta)$ を与える。内側で有界な調和拡張の平均と、半径 $1$ での値を比べる。

円周データを定数と三角関数に分ける。角振動数 $j$ について半径の方程式は $R^{\prime\prime}+r^{-1}R^{\prime}-j^2r^{-2}R=0$ である。元の条件に合う根を選び、境界半径で係数を規格化する。この例での関数は次のとおりである。

$$
\begin{aligned}
u(r,\theta)&=- \frac{r^{3} \sin{\left(3 \theta \right)}}{4} + \frac{3 r \cos{\left(\theta \right)}}{2} + 4\\
g&=- 2 \sin{\left(3 \theta \right)} + 3 \cos{\left(\theta \right)} + 4,\qquad r=2\\
u_r&=\frac{\partial}{\partial r} \left(- \frac{r^{3} \sin{\left(3 \theta \right)}}{4} + \frac{3 r \cos{\left(\theta \right)}}{2} + 4\right)=- \frac{3 \left(r^{2} \sin{\left(3 \theta \right)} - 2 \cos{\left(\theta \right)}\right)}{4}\\
u_\theta&=\frac{\partial}{\partial \theta} \left(- \frac{r^{3} \sin{\left(3 \theta \right)}}{4} + \frac{3 r \cos{\left(\theta \right)}}{2} + 4\right)=- \frac{3 r \left(r^{2} \cos{\left(3 \theta \right)} + 2 \sin{\left(\theta \right)}\right)}{4}
\end{aligned}
$$

二階微分を別々に計算する。半径方向の二階微分だけを角度方向の二階微分と足しても相殺しない。極座標の Laplacian の一階項も入れると各振動数で相殺する。

$$
\begin{aligned}
u_{rr}&=\frac{\partial^{2}}{\partial r^{2}} \left(- \frac{r^{3} \sin{\left(3 \theta \right)}}{4} + \frac{3 r \cos{\left(\theta \right)}}{2} + 4\right)=- \frac{3 r \sin{\left(3 \theta \right)}}{2}\\
u_{\theta\theta}&=\frac{\partial^{2}}{\partial \theta^{2}} \left(- \frac{r^{3} \sin{\left(3 \theta \right)}}{4} + \frac{3 r \cos{\left(\theta \right)}}{2} + 4\right)=\frac{3 r \left(3 r^{2} \sin{\left(3 \theta \right)} - 2 \cos{\left(\theta \right)}\right)}{4}\\
\Delta u&=- \frac{3 r \sin{\left(3 \theta \right)}}{2} + \frac{- \frac{3 r^{2} \sin{\left(3 \theta \right)}}{4} + \frac{3 \cos{\left(\theta \right)}}{2}}{r} + \frac{3 \left(3 r^{2} \sin{\left(3 \theta \right)} - 2 \cos{\left(\theta \right)}\right)}{4 r}=0
\end{aligned}
$$

境界における値の平均と法線微分を読む。$j\ne0$ の角方向の平均は零であり、全体の平均は定数項である。またエネルギーは体積積分の方法と Green の境界積分の方法の両方で計算できる。

$$
\begin{aligned}
\text{円周の平均}&=\int\limits_{- \pi}^{\pi} \frac{- 2 \sin{\left(3 \theta \right)} + 3 \cos{\left(\theta \right)} + 4}{2 \pi}\, d\theta=4\\
\partial_\nu u&=- 3 \sin{\left(3 \theta \right)} + \frac{3 \cos{\left(\theta \right)}}{2}=\frac{3 \left(- 2 \sin{\left(3 \theta \right)} + \cos{\left(\theta \right)}\right)}{2}\\
\int_{-\pi}^{\pi}|Du|^2\,\mathrm{d}\theta&=\frac{9 \pi r^{4}}{8} + \frac{9 \pi}{2}=\frac{9 \pi \left(r^{2} - 2 r + 2\right) \left(r^{2} + 2 r + 2\right)}{8}\\
\int |Du|^2&=\int\limits_{0}^{2} r \left(\frac{9 \pi r^{4}}{8} + \frac{9 \pi}{2}\right)\, dr=21 \pi\\
\int_{\partial\Omega} u\partial_\nu u&=\int\limits_{- \pi}^{\pi} \left(- 3 \sin{\left(3 \theta \right)} + \frac{3 \cos{\left(\theta \right)}}{2}\right) \left(- 4 \sin{\left(3 \theta \right)} + 6 \cos{\left(\theta \right)} + 8\right)\, d\theta=21 \pi
\end{aligned}
$$

次に同じ振動数の半径因子を比較し、角度を固定した一点で値と半径方向の傾きを調べる。角度の取り方によって、振幅のあるモードでも点での値が零になることがある。

$$
\begin{aligned}
u&=\frac{3 \sqrt{3}}{4} + \frac{15}{4}=\frac{3 \left(\sqrt{3} + 5\right)}{4},\qquad r=1,\theta=\frac{\pi}{6}\\
u_r&=- \frac{3}{4} + \frac{3 \sqrt{3}}{4}=\frac{3 \left(-1 + \sqrt{3}\right)}{4},\qquad r=1,\theta=\frac{\pi}{6}\\
u_\theta&=- \frac{3}{4},\qquad r=1,\theta=\frac{\pi}{6}
\end{aligned}
$$

どの同心円上でも平均は $4$ である。$r=1$ の境界に比べた振幅比は、一次モードが $1/2$、三次モードが $1/8$ である。平均値の等式は振動数零だけを残すが、中心以外の点では各振動数が異なる比で残る。
:::

## 最大値原理と一意性

:::theorem[調和関数の強最大値原理]{#thm-pde-harmonic-maximum}
連結開集合 $\Omega$ 上の調和関数が内点で $\sup_\Omega u$ を達成すれば定数である。$\Omega$ が有界で $u\in C(\overline\Omega)$ なら、最大値と最小値は境界で達成される。
:::

:::proof
最大値 $M$ をとる点を中心とした小球上で $M-u$ は連続かつ非負で、平均が $0$ である。球内の一点で正なら小近傍でも正となり平均が正になるので、球全体で $u=M$ である。最大点集合はこうして開であり、連続性により相対閉でもある。連結性により全体に等しい。有界の場合は閉包上で最大点があり、内点なら定数である。最小値には $-u$ を用いる。
:::

:::corollary[Dirichlet 問題の安定性]{#cor-pde-laplace-dirichlet}
有界領域で $-\Delta u_i=f$、$u_i|_{\partial\Omega}=g_i$ を満たす古典解があれば、$\norm{u_1-u_2}_{C(\overline\Omega)}\le\norm{g_1-g_2}_{C(\partial\Omega)}$ である。
:::

:::proof
差は調和であるので、:ref[thm-pde-harmonic-maximum] を差とその負に適用する。
:::

:::example[境界値から全域の符号を決める]
単位円盤で境界値が $2+\cos 3\theta$ の調和関数は $1\le u\le3$ である。境界値が一定でないので内点では狭義の不等号になる。実際 $u=2+r^3\cos3\theta$ は候補であり、差の最大値原理で一意である。
:::

:::example[右辺の差もあるとき]
$\Omega\subset B_R(0)$、$-\Delta w=f$、$w|_{\partial\Omega}=0$、$|f|\le F$ とする。$b(x)=F(R^2-|x|^2)/(2n)$ は $-\Delta b=F$ を満たす。$\Delta(w-b)\ge0$ の場合、$w-b+\eps|x|^2$ は内部最大点を持てない。最大点では Laplace 作用素が非正のはずだが、追加項により正になるからである。境界評価から $\eps\downarrow0$ として $w\le b$ であり、$-w$ に同じ議論を適用すれば $|w|\le b$ を得る。
:::

### Poisson 方程式に対する障壁評価

$\Omega\subset B_R(a)$ を有界領域とし、$u\in C^2(\Omega)\cap C(\overline\Omega)$ が $-\Delta u=f$、$u|_{\partial\Omega}=g$ を満たすとする。$w(x)=\norm f_\infty(R^2-|x-a|^2)/(2n)$ は非負で $-\Delta w=\norm f_\infty$ を満たす。$u-\norm g_\infty-w$ の Laplacian は非負で境界値は非正である。球面平均の単調性を用いた弱最大値原理からこれは非正である。$-u$ にも同じ議論を適用すると
$$
|u(x)|\le\norm g_\infty+\frac{R^2-|x-a|^2}{2n}\norm f_\infty .
$$
球に定数外力を置き $g=0$ とすれば $u=w$ が成立するので、球全体についての係数 $R^2/(2n)$ は改善できない。

:::example[源項が異なる二つの解]
二次元の半径 $2$ の円板で境界値が共通の二つの解を考える。源項の差が定数 $3$ なら解の差は $3(4-r^2)/4$ である。最大差は中心の $3$ で、境界から距離が $1$ の点では $9/4$ である。源項の差が $3$ 以下としか分からない場合にも同じ式が点ごとの上界になる。境界値も異なるなら境界差の一様ノルムを足して評価する。
:::

:::example[球上の障壁と点ごとの誤差]
三次元の半径 $2$ の球で $-\Delta u=6$、零境界とする。最大値原理から得た上界がどこで達せられるかを求める。

放射方向だけに依存する解を考える。半径 $r$ の球面積が $\sigma_{n-1}r^{n-1}$ であることから、発散を計算すると $-(r^{n-1}u_r)^{\prime}=r^{n-1}f$ となる。これを一度積分し、中心を含む場合には中心に余分な点源がないよう流束の極限を零にする。環状領域では内側の境界条件で積分定数を決める。

$$
\begin{aligned}
u&=4 - r^{2}\\
u_r&=\frac{d}{d r} \left(4 - r^{2}\right)=- 2 r\\
r^{n-1}u_r&=- 2 r^{3}\\
-r^{n-1}f&=\frac{d}{d r} \left(- 2 r^{3}\right)=- 6 r^{2}\\
-\Delta u&=6
\end{aligned}
$$

次に値と流束を境界で読む。外側では法線は半径が増える向きである。内側にも境界がある場合は法線が反対向きになる。源項の体積積分を角度方向に積分すると面積係数が現れ、流束との差を確かめられる。

$$
\begin{aligned}
u&=4,\qquad r=0\\
u&=0,\qquad r=2\\
u_r&=-4,\qquad r=2\\
\sigma_{n-1}&=4 \pi\\
\text{全源量}&=\int\limits_{0}^{2} 24 \pi r^{2}\, dr=64 \pi\\
\text{負の全流束}&=64 \pi
\end{aligned}
$$

観測半径で値と半径微分を計算する。勾配の二乗を用いてエネルギーも計算できる。観測点で勾配が小さい場合でも、球面積の因子があるため、そこを通る全流束が小さいとは限らない。

$$
\begin{aligned}
u&=3,\qquad r=1\\
u_r&=-2,\qquad r=1\\
\sigma_{n-1}r^{n-1}u_r&=- 8 \pi,\qquad r=1\\
\frac12\int |Du|^2&=\int\limits_{0}^{2} 8 \pi r^{4}\, dr=\frac{256 \pi}{5}
\end{aligned}
$$

境界で零、中心で $4$ となる。比較に使う $w=(R^2-r^2)\norm f_\infty/(2n)$ はこの解そのものである。時刻を含まない Poisson 方程式では、定数の源が球全体にあると中心まで値を増やす。
:::

## 調和関数の正則性と導関数の評価

平均値は平滑化しても値が変わらないことを意味する。そこで微分は未知関数ではなく既知の平均核に作用させる。

:::proposition[内部導関数評価]{#prop-pde-harmonic-derivative}
$u$ が $\overline B_R(x_0)$ の近傍で調和であれば任意の多重指数 $\alpha$ について

$$
|D^\alpha u(x_0)|\le C_{n,\alpha}R^{-n-|\alpha|}\int_{B_R(x_0)}|u(y)|\,\dd y .
\label{eq-pde-harmonic-derivative}
$$

また $m=|\alpha|\ge1$ に対し $|D^\alpha u(x_0)|\le (nm/R)^m\sup_{B_R(x_0)}|u|$ である。
:::

:::proof
:ref[thm-pde-harmonic-mean] の放射状平滑化で $u(x)=\int\rho_R(x-y)u(y)\,\dd y$ とする。$D^\alpha\rho_R$ の上界は $R^{-n-|\alpha|}\norm{D^\alpha\rho}_\infty$ なので最初の評価が従う。さらに平滑性により $\partial_j u$ も調和である。球体平均を用いて $\partial_j u(x)=|B_r|^{-1}\int_{B_r(x)}\partial_j u$ とし、発散定理で右辺を $|B_r|^{-1}\int_{\partial B_r(x)}u\nu_j\,\dd S$ に直すと $|\partial_j u(x)|\le(n/r)\sup_{B_r(x)}|u|$ である。半径 $R$ の球を幅 $R/m$ ずつ縮めながら、各段階で一階微分の評価を適用すると第 2 の式を得る。
:::

この評価から解析性も従う。$x\in B_{R/2}(x_0)$ なら $x$ を中心とする半径 $R/2$ の球で評価できる。$m^m\le e^m m!$ は $\sum_{j=1}^m\log j\ge\int_1^m\log s\,\dd s$ から従う。$h$ 方向の $m$ 階微分の絶対値は $\norm h_1^m$ と各偏微分の最大値の積以下である。線分上の Taylor 余項は $C(2en\norm h_1/R)^m$ 以下になるため、$h$ が十分小さければ Taylor 級数は $u(x_0+h)$ に収束する。

:::example[調和多項式の勾配]
$u=x_1^3-3x_1x_2^2$ は $\partial_{11}u=6x_1$、$\partial_{22}u=-6x_1$ なので調和である。$(1,0)$ では $\partial_1u=3$ である。単位円盤の中心では全ての一階、二階導関数が $0$ でも、三階の $\partial_{111}u=6$ が残る。
:::

:::example[境界への距離の役割]
単位円盤の $u_m(r,\theta)=r^m\cos(m\theta)$ は $\norm{u_m}_\infty\le1$ だが境界 $(1,0)$ の法線微分は $m$ である。したがって内部評価の定数を、境界からの距離に依存させず閉円盤に適用することはできない。
:::

### 分布解から調和代表を取り出す

:::proposition[Weyl の補題の局所形]{#prop-pde-weyl-harmonic}
$u\in L^1_{\mathrm{loc}}(\Omega)$ が分布の意味で $\Delta u=0$ を満たすなら、$u$ と a.e. に等しい滑らかな調和関数が存在する。
:::

:::proof
$u_\eps=\rho_\eps*u$ とする。切断の影響がない内側では $\Delta u_\eps=0$ である。$K\Subset\Omega$ に対し、ある $r>0$ で $K$ の $2r$ 近傍が $\Omega$ に含まれるとする。十分小さい $\eps,\delta$ なら $u_\eps-u_\delta$ に球体平均を使えて
$$
\sup_{x\in K}|u_\eps(x)-u_\delta(x)|
\le |B_r|^{-1}\norm{u_\eps-u_\delta}_{L^1(K+B_r)} .
$$
右辺は局所 $L^1$ での平滑化の収束により零へ行く。従って局所一様な極限 $v$ は連続で、平均値の等式を引き継ぐ。:ref[thm-pde-harmonic-mean] の逆から $v$ は調和であり、$L^1$ 極限の一意性から $v=u$ a.e. である。
:::

同じ議論で $\Delta u$ が滑らかな関数 $f$ の場合、局所で $f$ を切断して $-Nf$ を引くと残りが分布調和である。従って $u$ は古典的な $C^\infty$ 代表を持つ。既知の源項をポテンシャルに分けることで、調和関数の評価を一般の Poisson 解へ持ち込める。

:::example[高い角振動数で導関数がどこまで残るか]
単位円板の境界に $g=\cos(4\theta)+2\sin(5\theta)$ を課す。中心から離れて導関数を求める。

円周データを定数と三角関数に分ける。角振動数 $j$ について半径の方程式は $R^{\prime\prime}+r^{-1}R^{\prime}-j^2r^{-2}R=0$ である。元の条件に合う根を選び、境界半径で係数を規格化する。この例での関数は次のとおりである。

$$
\begin{aligned}
u(r,\theta)&=2 r^{5} \sin{\left(5 \theta \right)} + r^{4} \cos{\left(4 \theta \right)}\\
g&=2 \sin{\left(5 \theta \right)} + \cos{\left(4 \theta \right)},\qquad r=1\\
u_r&=\frac{\partial}{\partial r} \left(2 r^{5} \sin{\left(5 \theta \right)} + r^{4} \cos{\left(4 \theta \right)}\right)=2 r^{3} \left(5 r \sin{\left(5 \theta \right)} + 2 \cos{\left(4 \theta \right)}\right)\\
u_\theta&=\frac{\partial}{\partial \theta} \left(2 r^{5} \sin{\left(5 \theta \right)} + r^{4} \cos{\left(4 \theta \right)}\right)=2 r^{4} \left(5 r \cos{\left(5 \theta \right)} - 2 \sin{\left(4 \theta \right)}\right)
\end{aligned}
$$

二階微分を別々に計算する。半径方向の二階微分だけを角度方向の二階微分と足しても相殺しない。極座標の Laplacian の一階項も入れると各振動数で相殺する。

$$
\begin{aligned}
u_{rr}&=\frac{\partial^{2}}{\partial r^{2}} \left(2 r^{5} \sin{\left(5 \theta \right)} + r^{4} \cos{\left(4 \theta \right)}\right)=4 r^{2} \left(10 r \sin{\left(5 \theta \right)} + 3 \cos{\left(4 \theta \right)}\right)\\
u_{\theta\theta}&=\frac{\partial^{2}}{\partial \theta^{2}} \left(2 r^{5} \sin{\left(5 \theta \right)} + r^{4} \cos{\left(4 \theta \right)}\right)=- 2 r^{4} \left(25 r \sin{\left(5 \theta \right)} + 8 \cos{\left(4 \theta \right)}\right)\\
\Delta u&=4 r^{2} \left(10 r \sin{\left(5 \theta \right)} + 3 \cos{\left(4 \theta \right)}\right) - 2 r^{2} \left(25 r \sin{\left(5 \theta \right)} + 8 \cos{\left(4 \theta \right)}\right) + \frac{10 r^{4} \sin{\left(5 \theta \right)} + 4 r^{3} \cos{\left(4 \theta \right)}}{r}=0
\end{aligned}
$$

境界における値の平均と法線微分を読む。$j\ne0$ の角方向の平均は零であり、全体の平均は定数項である。またエネルギーは体積積分の方法と Green の境界積分の方法の両方で計算できる。

$$
\begin{aligned}
\text{円周の平均}&=\int\limits_{- \pi}^{\pi} \frac{2 \sin{\left(5 \theta \right)} + \cos{\left(4 \theta \right)}}{2 \pi}\, d\theta=0\\
\partial_\nu u&=10 \sin{\left(5 \theta \right)} + 4 \cos{\left(4 \theta \right)}=2 \left(5 \sin{\left(5 \theta \right)} + 2 \cos{\left(4 \theta \right)}\right)\\
\int_{-\pi}^{\pi}|Du|^2\,\mathrm{d}\theta&=200 \pi r^{8} + 32 \pi r^{6}=8 \pi r^{6} \left(25 r^{2} + 4\right)\\
\int |Du|^2&=\int\limits_{0}^{1} r \left(200 \pi r^{8} + 32 \pi r^{6}\right)\, dr=24 \pi\\
\int_{\partial\Omega} u\partial_\nu u&=\int\limits_{- \pi}^{\pi} \left(2 \sin{\left(5 \theta \right)} + \cos{\left(4 \theta \right)}\right) \left(10 \sin{\left(5 \theta \right)} + 4 \cos{\left(4 \theta \right)}\right)\, d\theta=24 \pi
\end{aligned}
$$

次に同じ振動数の半径因子を比較し、角度を固定した一点で値と半径方向の傾きを調べる。角度の取り方によって、振幅のあるモードでも点での値が零になることがある。

$$
\begin{aligned}
u&=\frac{1}{16},\qquad r=\frac{1}{2},\theta=0\\
u_r&=\frac{1}{2},\qquad r=\frac{1}{2},\theta=0\\
u_\theta&=\frac{5}{16},\qquad r=\frac{1}{2},\theta=0
\end{aligned}
$$

この調和多項式の中心での三階までの微分は零になる。例えば放射軸 $\theta=0$ に沿えば最初の寄与は $r^4$ であり、四階微分で $24$ になる。一方、中心から半分の円における値の各振幅は $1/16$ と $2/32$ である。内部導関数の評価で境界との距離が関係することが、この違いに現れる。
:::

### 調和成分と源の成分を分ける局所計算

平面で $-\Delta u=x^2-y^2$ とし、内部の多項式解を探す。$P=x^2-y^2$ は二次の同次調和多項式なので、積の微分から $\Delta(|x|^2 P)=4P+4(x\partial_x+y\partial_y)P=12P$ である。従って $u_0=-(x^2+y^2)(x^2-y^2)/12$ は特解であり、任意の $u\in L^1_{\mathrm{loc}}$ の分布解との差は Weyl の補題により調和な代表を持つ。例えば $u_0=-x^4/12+y^4/12$ であり、$u_{0,xx}=-x^2$、$u_{0,yy}=y^2$ と直接確かめられる。境界値をまだ指定していなければ、$1$、$x$、$xy$ などを加えたものは同じ源を持つ。正則性を求める際は特解の滑らかさと調和部分の評価を別々に扱える。

## Liouville の定理と Harnack の不等式

:::theorem[Liouville の定理]{#thm-pde-harmonic-liouville}
全空間で調和かつ一様有界な関数は定数である。より一般に $|u(x)|\le C(1+|x|^k)$、$k$ が非負整数なら $u$ は次数 $k$ 以下の調和多項式である。
:::

:::proof
任意の $x$ に対し次数 $m=k+1$ の導関数の評価を半径 $R$ の球で行うと、右辺は固定した $x$ に対し $O(R^{k-m})$ である。$R\to\infty$ でその導関数が消える。任意の線分上で Taylor の公式を使えば $u$ は次数 $k$ 以下の多項式である。$k=0$ が最初の主張である。
:::

:::theorem[球上の Harnack 評価]{#thm-pde-harmonic-harnack}
$u\ge0$ が $B_{4r}(x_0)$ 上で調和とする。$x,y\in B_r(x_0)$ に対し $u(x)\le3^n u(y)$ である。
:::

:::proof
$B_r(x)\subset B_{2r}(x_0)\subset B_{3r}(y)\subset B_{4r}(x_0)$ である。非負性により球体平均の分子は各包含で増える。$u(x)=|B_r|^{-1}\int_{B_r(x)}u\le |B_r|^{-1}\int_{B_{3r}(y)}u=3^n u(y)$ である。
:::

正値の調和関数は全空間で定数である。実際、非負な $u$ に対し :eqref[eq-pde-harmonic-derivative] の一階評価を使えば $|Du(x)|\le C R^{-1}u(x)$ である。$R\to\infty$ で勾配が消える。したがって片側だけ有界でも、定数を足し符号を変えることで Liouville の結論が従う。

:::example[小球を鎖でつなぐ]
二つの点を結ぶ折れ線が境界から距離 $4r$ より離れ、各線分の長さが $r/2$ 以下なら、各頂点で半径 $r$ の Harnack 評価を適用して両端の比を $3^{nN}$ で抑えられる。ここで $N$ は線分数である。二次元で $N=2$ なら上界は $81$ である。領域の連結性は点をつなげるが、境界に近づくと必要な球の半径は小さくなる。
:::

:::example[符号条件の意味]
$u(x)=x_1$ は全空間で調和であるが原点で $0$ となり、任意の小球内で正負の値をとる。Harnack 評価を符号の条件なしに当てはめると正の値が $0$ 以下になるため、非負性は不可欠である。
:::

### 増大度と調和多項式

$u$ が全空間で調和であり $|u(x)|\le C(1+|x|^m)$、$m$ は非負整数とする。任意の $a$ と $|\alpha|=m+1$ について、半径 $R$ の内部導関数評価は
$$
|D^\alpha u(a)|\le C_{\alpha,n}R^{-(m+1)}C(1+(|a|+R)^m)
$$
を与える。$R\to\infty$ で右辺は零である。従って原点から点 $x$ までの線分で Taylor 公式を使うと余項が零になり、$u$ は次数 $m$ 以下の調和多項式である。逆に各調和多項式はその次数の多項式増大度を満たす。例えば $m=1$ では全てのアフィン関数が出る。$m=2$ で二次項を $x^\top A x$ と置けば $\Delta(x^\top A x)=2\operatorname{tr} A$ なので、二次の制約は跡が零であることだけである。

:::example[一点で極小を取れる一次調和多項式]
$u$ が全平面で $|u(x)|\le 1+2|x|$ を満たし、原点で大域的最小値を取るとする。上の議論から $u=a+b\cdot x$ であり、$b\ne0$ なら $x=-tb$ で $t\to\infty$ により $u\to-\infty$ である。従って $b=0$ である。さらに原点で $u(0)=2$ を要求すると増大度の原点での上界 $1$ と矛盾し、解はない。
:::

:::example[正の調和関数が境界近くで大きく変化する]
単位円板で $g=4+3\cos\theta$、すなわち拡張が正の線形関数である例を調べる。

円周データを定数と三角関数に分ける。角振動数 $j$ について半径の方程式は $R^{\prime\prime}+r^{-1}R^{\prime}-j^2r^{-2}R=0$ である。元の条件に合う根を選び、境界半径で係数を規格化する。この例での関数は次のとおりである。

$$
\begin{aligned}
u(r,\theta)&=3 r \cos{\left(\theta \right)} + 4\\
g&=3 \cos{\left(\theta \right)} + 4,\qquad r=1\\
u_r&=\frac{\partial}{\partial r} \left(3 r \cos{\left(\theta \right)} + 4\right)=3 \cos{\left(\theta \right)}\\
u_\theta&=\frac{\partial}{\partial \theta} \left(3 r \cos{\left(\theta \right)} + 4\right)=- 3 r \sin{\left(\theta \right)}
\end{aligned}
$$

二階微分を別々に計算する。半径方向の二階微分だけを角度方向の二階微分と足しても相殺しない。極座標の Laplacian の一階項も入れると各振動数で相殺する。

$$
\begin{aligned}
u_{rr}&=\frac{\partial^{2}}{\partial r^{2}} \left(3 r \cos{\left(\theta \right)} + 4\right)=0\\
u_{\theta\theta}&=\frac{\partial^{2}}{\partial \theta^{2}} \left(3 r \cos{\left(\theta \right)} + 4\right)=- 3 r \cos{\left(\theta \right)}\\
\Delta u&=0
\end{aligned}
$$

境界における値の平均と法線微分を読む。$j\ne0$ の角方向の平均は零であり、全体の平均は定数項である。またエネルギーは体積積分の方法と Green の境界積分の方法の両方で計算できる。

$$
\begin{aligned}
\text{円周の平均}&=\int\limits_{- \pi}^{\pi} \frac{3 \cos{\left(\theta \right)} + 4}{2 \pi}\, d\theta=4\\
\partial_\nu u&=3 \cos{\left(\theta \right)}\\
\int_{-\pi}^{\pi}|Du|^2\,\mathrm{d}\theta&=18 \pi\\
\int |Du|^2&=\int\limits_{0}^{1} 18 \pi r\, dr=9 \pi\\
\int_{\partial\Omega} u\partial_\nu u&=\int\limits_{- \pi}^{\pi} 3 \left(3 \cos{\left(\theta \right)} + 4\right) \cos{\left(\theta \right)}\, d\theta=9 \pi
\end{aligned}
$$

次に同じ振動数の半径因子を比較し、角度を固定した一点で値と半径方向の傾きを調べる。角度の取り方によって、振幅のあるモードでも点での値が零になることがある。

$$
\begin{aligned}
u&=\frac{5}{2},\qquad r=\frac{1}{2},\theta=\pi\\
u_r&=-3,\qquad r=\frac{1}{2},\theta=\pi\\
u_\theta&=0,\qquad r=\frac{1}{2},\theta=\pi
\end{aligned}
$$

半径 $s$ の円板における最大・最小の比は $(4+3s)/(4-3s)$ であり、$s=1/2$ では $11/5$ である。この比は中心値 $4$ だけからは出ず、半径を固定して初めて内側の範囲が決まる。境界に零を許した $4+4r\cos\theta$ なら $s\uparrow1$ で比は発散する。
:::

### 正値だけから得る全空間での定数性

全空間で非負な調和関数 $u$ を考える。中心 $x$ の球体平均で一階微分を表すと、$r>0$ に対して $|\partial_j u(x)|\le |B_r|^{-1}\int_{\partial B_r(x)}u\,\dd S$ である。ここでは $u\ge0$ なので球面上の絶対値を外せる。球面平均の等式をあわせると
$$
|\partial_j u(x)|\le \frac{n}{r}u(x).
$$
$u(x)$ は固定された有限値なので、$r\to\infty$ で $\partial_j u(x)=0$ となる。従って $u$ は定数である。$u$ が下に有界なら、下界を引いて同じ議論を使える。上下の両方を押さえた一様ノルムは必要ない。例えば $\R^2$ 上の $e^x\cos y$ は調和で非定数だが、$y=0,\pi$ の両直線を比較すると上にも下にも有界でないことが分かる。

調和多項式を同次部分に分けると、Laplacian が各部分の次数をちょうど二つ下げるため、異なる次数からの寄与は相殺できない。従って各同次部分も調和である。三変数で次数二の項は対称行列の六成分で表せ、跡が零という一つの独立な条件が付く。従って空間の次元は五であり、例えば $xy,yz,zx,x^2-y^2,y^2-z^2$ が基底である。$x^2+y^2+z^2$ を加えると Laplacian が $6$ だけ変化し、調和性は失われる。

## Green 関数と Poisson 積分

### Green の表示

:::definition[Dirichlet Green 関数]{#def-pde-dirichlet-green}
固定した $x\in\Omega$ に対して $G(x,y)=\Gamma(y-x)-h_x(y)$ とし、$h_x$ が $\Omega$ で調和で、境界値 $\Gamma(y-x)$ を持つとする。このとき $G$ を :term[Green 関数]{en="Green function" yomi="ぐりーんかんすう"} と呼ぶ。$G(x,\cdot)$ は境界で $0$、内部で $-\Delta_yG=\delta_x$ を満たす。
:::

有界な連結 $C^2$ 領域で $-\Delta u=f$、$u|_{\partial\Omega}=g$ の古典解を考え、$u$ と $h_x$ は境界まで $C^1$ とする。点 $x$ の周囲を除いて Green の第二公式を用いる。小球側の極限は基本解の規格化で $u(x)$ を与え、外部境界では $G=0$ である。従って

$$
u(x)=\int_\Omega G(x,y)f(y)\,\dd y
-\int_{\partial\Omega}\partial_{\nu_y}G(x,y)g(y)\,\dd S_y
\label{eq-pde-green-representation}
$$

となる。$G$ の対称性は $G(x,\cdot)$ と $G(z,\cdot)$ に同じ公式を二つの穴の外で適用すると得られる。外部境界項はどちらも $0$ で、二つの小球から $G(x,z)$ と $-G(z,x)$ が残る。$G$ は内部で正である。小球境界で $\Gamma$ が正に発散し $h_x$ は有界なので、小球を除いた領域の最小値原理を使って $G\ge0$ を得る。強最大値原理で内点では正になる。

### 反転と孤立特異点

:::proposition[Kelvin 変換]{#prop-pde-kelvin}
$n\ge2$ とし、原点を含まない領域で定義した $C^2$ 関数 $u$ に対し $v(x)=|x|^{2-n}u(x/|x|^2)$ とする。反転した領域上では $\Delta v(x)=|x|^{-n-2}(\Delta u)(x/|x|^2)$ である。
:::

:::proof
極座標の Laplace 作用素を $\partial_{rr}+(n-1)\partial_r/r+r^{-2}\Delta_{\mathbb S^{n-1}}$ と書く。これは勾配を動径成分と接成分に分け、動径の発散を $r^{1-n}\partial_r(r^{n-1}\partial_r)$、接成分の発散を $r^{-2}\Delta_{\mathbb S^{n-1}}$ と書いたものである。$s=1/r$ とし $v=r^{2-n}u(s,\theta)$ を動径で二度微分すると

$$
v_{rr}+\tfrac{n-1}{r}v_r=r^{-n-2}u_{ss}+(n-1)r^{-n-1}u_s .
$$

ゼロ階の項は $(2-n)(1-n)+(n-1)(2-n)=0$ で消える。接方向の項は $r^{-n}\Delta_{\mathbb S^{n-1}}u$ である。ここで $s^{-1}=r$、$s^{-2}=r^2$ と見れば主張が従う。
:::

例えば三次元の $u(x)=x_1$ は $x_1/|x|^3$ に移り、これは原点以外で調和である。反転前に定数だった関数は $r^{2-n}$ に移るため、この変換は基本解と境界の像を結び付ける。

:::proposition[有界な孤立特異点の除去]{#prop-pde-removable-harmonic}
$n\ge2$ とする。$B_R\setminus\{0\}$ 上で調和で原点近くで有界な関数は $B_R$ に調和に延長できる。
:::

:::proof
半径 $r$ の球面平均を $m(r)$ とする。二つの球面の間で発散定理を使うと、流束 $\int_{\partial B_r}\partial_r u$ は $r$ によらない。従って $m'(r)$ は $r^{1-n}$ の定数倍である。$n\ge2$ ではこの関数の原始関数が原点で発散するが、$m(r)$ は有界である。従って定数は零であり、流束は零になる。

半径 $r$ の球面上では :ref[prop-pde-harmonic-derivative] を半径 $r/2$ の内球で使うと $Du=O(r^{-1})$ である。小球を除いて試験関数 $\phi$ と Green の第二公式を使うと、境界項は $\int_{\partial B_r}(u\partial_r\phi-\phi\partial_r u)$ である。流束が零なので、第二項で $\phi$ を $\phi-\phi(0)$ に替えられる。両項ともに $O(r^{n-1})$ で消え、$u$ は球全体で分布調和となる。:ref[prop-pde-weyl-harmonic] が調和な代表を与える。この代表は原点以外で元の $u$ と連続性により一致する。
:::

:::example[境界に一致しても特異性が残る場合]
二次元で $u=\log r$ は単位円周で $0$ であり、原点外では調和だが有界でない。$-\Delta u=-2\pi\delta_0$ なので、値を一点で選び直しても球上の調和関数にはならない。境界値の一意性は領域全体で調和性を要求する。
:::

### 球と半空間での像

$B_R(0)$、$n\ge3$、$x\ne0$ とし、$x^*=R^2x/|x|^2$ と置く。$|y|=R$ では $|y-x^*|=(R/|x|)|y-x|$ なので、

$$
G(x,y)=\Gamma(y-x)-(R/|x|)^{n-2}\Gamma(y-x^*)
$$

は境界で $0$ になる。$x^*$ は外部にあるので引いた項は球内で調和である。$n=2$ では $\Gamma(y-x)-\Gamma(y-x^*)-(2\pi)^{-1}\log(R/|x|)$ とする。$x=0$ では $\Gamma(y)-\Gamma(R)$ である。法線方向に微分すると

$$
P_R(x,y)=-\partial_{\nu_y}G(x,y)=\frac{R^2-|x|^2}{\sigma_{n-1}R|x-y|^n},\qquad |y|=R .
\label{eq-pde-ball-poisson}
$$

実際、最初の基本解の負の法線微分は $(R^2-y\cdot x)/(\sigma_{n-1}R|x-y|^n)$ である。像の項では距離の関係と $y\cdot x^*=R^2y\cdot x/|x|^2$ により、差の分子が $R^2-y\cdot x-(|x|^2-y\cdot x)$ になる。

:::proposition[連続境界値の Poisson 積分]{#prop-pde-poisson-ball}
$g\in C(\partial B_R)$ に対し $u(x)=\int_{\partial B_R}P_R(x,y)g(y)\,\dd S_y$ は球内で調和、閉球に連続に延長でき、境界値 $g$ を持つ。
:::

:::proof
核は $x$ の関数として調和である。直接には $|x-y|^{-n+2}$ が調和であることと $R^2-|x|^2=-|x-y|^2-2y\cdot(x-y)$ を使うと確認できる。$n=2$ では $y\cdot(x-y)|x-y|^{-2}$ は $\log|x-y|$ の方向微分である。よって内球上で微分を積分に入れられる。定数解 $1$ に :eqref[eq-pde-green-representation] を適用すると $\int P_R\,\dd S=1$ である。$P_R\ge0$ であり、$x\to y_0$ のとき $|y-y_0|\ge\eta$ の範囲で分子は $0$ に収束し分母は下から有界である。この範囲の積分は $0$ に収束する。残りの範囲では $g(y)-g(y_0)$ が一様に小さいため、$u(x)\to g(y_0)$ である。
:::

:::example[円盤の Fourier モード]
半径 $2$ の円で境界値 $g=3\sin(2\theta)$ を与える。$3(r/2)^2\sin(2\theta)=3x_1x_2/2$ が同じ境界値を持つ調和関数なので、Poisson 積分と一致する。点 $(1,1)$ では $3/2$ である。
:::

:::example[半空間の反射]
$x=(x',h)$、$h>0$ を境界面で反射した点を $x^*$ とする。$G=\Gamma(y-x)-\Gamma(y-x^*)$ は平面で $0$ であり、外法線 $-e_n$ で微分した Poisson 核は $2h/(\sigma_{n-1}(|x'-y'|^2+h^2)^{n/2})$ である。二次元では $h/(\pi(s^2+h^2))$ となるので、境界の $[-1,1]$ に値 $1$、外側に値 $0$ を置くと点 $(0,1)$ の値は $2\arctan(1)/\pi=1/2$ である。
:::

半空間の Poisson 核 $P_h(s)=2h/(\sigma_{n-1}(|s|^2+h^2)^{n/2})$ を $s\in\R^{n-1}$ 上で積分すると $1$ になる。実際、点 $(s,h)$ の方向を単位上半球へ写すと、球面面積要素は $h(|s|^2+h^2)^{-n/2}\dd s$ である。この係数は平面上の面積を半径方向に垂直な面へ投影し、距離の $n-1$ 乗で割ることで得られる。上半球の面積は $\sigma_{n-1}/2$ である。また $P_h(s)=h^{1-n}P_1(s/h)$ なので、任意の $\delta>0$ に対して $\int_{|s|>\delta}P_h(s)\dd s\to0$ である。したがって有界連続な境界値 $g$ に対する $P_h*g$ は、境界に近づくと $g$ に収束する。核は $h>0$ で境界点を極とする基本解の法線微分の定数倍であり、調和である。各内点の近傍で、二階までの微分は遠方で可積分な関数に支配されるので、積分で得た関数も調和である。

この半空間の解は有界な解の中で一意である。差 $w$ は境界で零であり、境界面を横切って奇関数として反射する。反射した関数が調和であることは、境界面を中心とする任意の小球で確かめられる。その球面上の奇対称な値を球の Poisson 積分で延長すると、一意性により球内でも奇対称である。上半球では球面部分と平面部分の双方で $w$ と値が一致するので、最大値原理により上半球全体で一致する。下半球でも奇対称性により一致する。したがって反射した関数は全空間で有界かつ調和であり、Liouville の定理と奇対称性により零である。例えば $w(x)=x_n$ は同じ零境界値を持つ無界な調和関数である。

球の Poisson 核から、中心と内点を比較する Harnack 評価も得られる。$u\ge0$ が閉球の近傍で調和とし、$r=|x|<R$ とする。球面上では $R-r\le|x-y|\le R+r$ なので
$$
\frac{(R-r)R^{n-2}}{(R+r)^{n-1}}u(0)
\le u(x)\le
\frac{(R+r)R^{n-2}}{(R-r)^{n-1}}u(0).
$$
ここでは $P_R(x,y)$ を $P_R(0,y)=1/(\sigma_{n-1}R^{n-1})$ で割り、非負な境界値に対して積分した。両係数は $x$ から最遠または最近の球面点に核が取る比である。その点の近くに台を持つ非負な連続境界値を取り、台を縮めるとこの比に近づく。

### Green 関数の対称性と領域の包含

滑らかな有界領域で $G(x,y)=\Gamma(x-y)-h_x(y)$ とし、$h_x$ は境界で $\Gamma(x-y)$ と一致する調和関数とする。異なる $a,b$ の周りの半径 $\eps$ の球を除いた領域に Green の第二恒等式を適用する。$G(a,\cdot),G(b,\cdot)$ は外側の境界で零なので外側の積分は零である。小球の法線は領域から見ると中心向きであり、その法線に沿う基本解の微分の積分は $1$ である。そこで $a$ 側の積分では滑らかな因子を $G(a,b)$ に置き換えられ、$b$ 側では $G(b,a)$ になる。他方、基本解そのものと滑らかな因子の微分の積は面積を掛けると $O(\eps)$（二次元で $O(\eps|\log\eps|)$）になる。従って差の極限が零になり $G(a,b)=G(b,a)$ である。

$G(y,\cdot)$ は小球の境界で正、外境界で零なので比較により非負である。$\Omega_1\subset\Omega_2$ に対して同じ極 $y\in\Omega_1$ を取れば、$G_2-G_1$ では特異性が打ち消されて $\Omega_1$ で調和となる。境界で $G_1=0$、$G_2\ge0$ なので $G_1\le G_2$。従って非負の源に対する零 Dirichlet 解は、領域を広げると小さくならない。例えば三次元の同心球で中心の極に対し $G_R(x,0)=(|x|^{-1}-R^{-1})/(4\pi)$ であり、$R=1$ から $R=2$ に広げると共通部分の差は一定 $1/(8\pi)$ である。

:::example[同心環での調和測度の重み]
三次元の環 $1<r<4$ において内面の値を $1$、外面の値を $0$ とする。

放射方向だけに依存する解を考える。半径 $r$ の球面積が $\sigma_{n-1}r^{n-1}$ であることから、発散を計算すると $-(r^{n-1}u_r)^{\prime}=r^{n-1}f$ となる。これを一度積分し、中心を含む場合には中心に余分な点源がないよう流束の極限を零にする。環状領域では内側の境界条件で積分定数を決める。

$$
\begin{aligned}
u&=- \frac{1}{3} + \frac{4}{3 r}\\
u_r&=\frac{d}{d r} \left(- \frac{1}{3} + \frac{4}{3 r}\right)=- \frac{4}{3 r^{2}}\\
r^{n-1}u_r&=- \frac{4}{3}\\
-r^{n-1}f&=\frac{d}{d r} \left(- \frac{4}{3}\right)=0\\
-\Delta u&=0
\end{aligned}
$$

次に値と流束を境界で読む。外側では法線は半径が増える向きである。内側にも境界がある場合は法線が反対向きになる。源項の体積積分を角度方向に積分すると面積係数が現れ、流束との差を確かめられる。

$$
\begin{aligned}
u&=1,\qquad r=1\\
u&=0,\qquad r=4\\
u_r&=- \frac{1}{12},\qquad r=4\\
\sigma_{n-1}&=4 \pi\\
\text{全源量}&=\int\limits_{1}^{4} 0\, dr=0\\
\text{負の全流束}&=0
\end{aligned}
$$

観測半径で値と半径微分を計算する。勾配の二乗を用いてエネルギーも計算できる。観測点で勾配が小さい場合でも、球面積の因子があるため、そこを通る全流束が小さいとは限らない。

$$
\begin{aligned}
u&=\frac{1}{3},\qquad r=2\\
u_r&=- \frac{1}{3},\qquad r=2\\
\sigma_{n-1}r^{n-1}u_r&=- \frac{16 \pi}{3},\qquad r=2\\
\frac12\int |Du|^2&=\int\limits_{1}^{4} \frac{32 \pi}{9 r^{2}}\, dr=\frac{8 \pi}{3}
\end{aligned}
$$

半径 $2$ では内側の境界データの重みが $1/3$ である。外側の値を $b$、内側の値を $a$ に変えれば解は $b+(a-b)u$。全源量は零なので外面と内面の流束は相殺し、各面の半径微分の大きさが違っても全流束は等しい。
:::

:::example[Poisson 積分をモードで評価する]
半径 $3$ の円板、$g=1+4\cos2\theta+2\sin\theta$ とする。Poisson 核の円周上での Fourier 係数 $(r/3)^j$ を使う。

円周データを定数と三角関数に分ける。角振動数 $j$ について半径の方程式は $R^{\prime\prime}+r^{-1}R^{\prime}-j^2r^{-2}R=0$ である。元の条件に合う根を選び、境界半径で係数を規格化する。この例での関数は次のとおりである。

$$
\begin{aligned}
u(r,\theta)&=\frac{4 r^{2} \cos{\left(2 \theta \right)}}{9} + \frac{2 r \sin{\left(\theta \right)}}{3} + 1\\
g&=2 \sin{\left(\theta \right)} + 4 \cos{\left(2 \theta \right)} + 1,\qquad r=3\\
u_r&=\frac{\partial}{\partial r} \left(\frac{4 r^{2} \cos{\left(2 \theta \right)}}{9} + \frac{2 r \sin{\left(\theta \right)}}{3} + 1\right)=\frac{2 \left(4 r \cos{\left(2 \theta \right)} + 3 \sin{\left(\theta \right)}\right)}{9}\\
u_\theta&=\frac{\partial}{\partial \theta} \left(\frac{4 r^{2} \cos{\left(2 \theta \right)}}{9} + \frac{2 r \sin{\left(\theta \right)}}{3} + 1\right)=- \frac{2 r \left(8 r \sin{\left(\theta \right)} - 3\right) \cos{\left(\theta \right)}}{9}
\end{aligned}
$$

二階微分を別々に計算する。半径方向の二階微分だけを角度方向の二階微分と足しても相殺しない。極座標の Laplacian の一階項も入れると各振動数で相殺する。

$$
\begin{aligned}
u_{rr}&=\frac{\partial^{2}}{\partial r^{2}} \left(\frac{4 r^{2} \cos{\left(2 \theta \right)}}{9} + \frac{2 r \sin{\left(\theta \right)}}{3} + 1\right)=\frac{8 \cos{\left(2 \theta \right)}}{9}\\
u_{\theta\theta}&=\frac{\partial^{2}}{\partial \theta^{2}} \left(\frac{4 r^{2} \cos{\left(2 \theta \right)}}{9} + \frac{2 r \sin{\left(\theta \right)}}{3} + 1\right)=- \frac{2 r \left(8 r \cos{\left(2 \theta \right)} + 3 \sin{\left(\theta \right)}\right)}{9}\\
\Delta u&=\frac{8 \cos{\left(2 \theta \right)}}{9} + \frac{\frac{8 r \cos{\left(2 \theta \right)}}{9} + \frac{2 \sin{\left(\theta \right)}}{3}}{r} - \frac{2 \left(8 r \cos{\left(2 \theta \right)} + 3 \sin{\left(\theta \right)}\right)}{9 r}=0
\end{aligned}
$$

境界における値の平均と法線微分を読む。$j\ne0$ の角方向の平均は零であり、全体の平均は定数項である。またエネルギーは体積積分の方法と Green の境界積分の方法の両方で計算できる。

$$
\begin{aligned}
\text{円周の平均}&=\int\limits_{- \pi}^{\pi} \frac{2 \sin{\left(\theta \right)} + 4 \cos{\left(2 \theta \right)} + 1}{2 \pi}\, d\theta=1\\
\partial_\nu u&=\frac{2 \sin{\left(\theta \right)}}{3} + \frac{8 \cos{\left(2 \theta \right)}}{3}=\frac{2 \left(\sin{\left(\theta \right)} + 4 \cos{\left(2 \theta \right)}\right)}{3}\\
\int_{-\pi}^{\pi}|Du|^2\,\mathrm{d}\theta&=\frac{128 \pi r^{2}}{81} + \frac{8 \pi}{9}=\frac{8 \pi \left(16 r^{2} + 9\right)}{81}\\
\int |Du|^2&=\int\limits_{0}^{3} r \left(\frac{128 \pi r^{2}}{81} + \frac{8 \pi}{9}\right)\, dr=36 \pi\\
\int_{\partial\Omega} u\partial_\nu u&=\int\limits_{- \pi}^{\pi} \left(\frac{2 \sin{\left(\theta \right)}}{3} + \frac{8 \cos{\left(2 \theta \right)}}{3}\right) \left(6 \sin{\left(\theta \right)} + 12 \cos{\left(2 \theta \right)} + 3\right)\, d\theta=36 \pi
\end{aligned}
$$

次に同じ振動数の半径因子を比較し、角度を固定した一点で値と半径方向の傾きを調べる。角度の取り方によって、振幅のあるモードでも点での値が零になることがある。

$$
\begin{aligned}
u&=1,\qquad r=\frac{3}{2},\theta=\frac{\pi}{2}\\
u_r&=- \frac{2}{3},\qquad r=\frac{3}{2},\theta=\frac{\pi}{2}\\
u_\theta&=0,\qquad r=\frac{3}{2},\theta=\frac{\pi}{2}
\end{aligned}
$$

$q=r/3$ とすれば核の角部分は $1+2\sum_{j\ge1}q^j\cos(j(\theta-\phi))$ である。これは $\sum q^je^{ij(\theta-\phi)}$ の幾何級数と共役を足して得られる。これを $g$ と積分し、三角関数の直交性を使うと表の半径因子になる。定数の係数が $1$ であることから境界核の全量も $1$ である。
:::

### 非負の孤立特異点の強さ

$n\ge3$、$u\ge0$ が $B_1\setminus\{0\}$ で調和であるとする。球面平均 $m(r)$ は放射方向の調和式 $(r^{n-1}m')'=0$ を満たす。これは球面平均の微分の式と、二つの球の間で発散定理を使って $\int_{\partial B_r}\partial_r u$ が一定になることから分かる。従って $m(r)=a+b r^{2-n}$、非負性から $b\ge0$ である。半径 $r$ の球面は環 $r/2<|x|<2r$ 内で有限個の球によりつなげられる。尺度を $1/r$ 倍した固定環で Harnack 不等式を有限回使えば、球面上の最大値は最小値の定数倍以下、従って $C m(r)$ 以下である。従って $u=O(r^{2-n})$、導関数評価から $Du=O(r^{1-n})$ である。

この評価で $u$ は原点まで局所可積分である。試験関数 $\phi\in C_c^\infty(B_1)$ に対し、小球を除いて Green の公式を使う。境界法線が中心向きなので
$$
\int_{|x|>r}u(-\Delta\phi)
=\int_{|x|=r}(u\partial_r\phi-\phi\partial_r u).
$$
第一項は $O(r)$、第二項で $\phi(x)$ を $\phi(0)$ に替えた誤差も $O(r)$ である。平均の式から $\int_{|x|=r}\partial_r u=(2-n)\sigma_{n-1}b$。従って $-\Delta u=(n-2)\sigma_{n-1}b\delta_0$ となる。基本解をこの係数で引けば差は分布調和であり、:ref[prop-pde-weyl-harmonic] により原点まで滑らかに調和延長できる。すなわち $u=b|x|^{2-n}+v$ で、$v$ は球全体で調和である。

二次元では $m(r)=a+b\log(1/r)$ となり、非負性から $b\ge0$ である。同じ環上の Harnack 評価で $u=O(1+|\log r|)$、内部導関数評価で $Du=O((1+|\log r|)/r)$ になる。小円上の二つの誤差はともに $O(r(1+|\log r|))$ で消え、全流束は $\int_{|x|=r}\partial_r u=-2\pi b$ である。したがって $u=b\log(1/|x|)+v$ と分解できる。これらの分解を :term[Bôcher の定理]{en="Bôcher's theorem" yomi="ぼしぇーのていり"} という。例えば三次元で $u=2/r+3+x_1$ なら特異成分の係数は $2$、分布の源は $8\pi\delta_0$ である。非負性を外すと $\partial_{x_1}\Gamma$ のように強い特異点を持つ例があり、球面平均だけでは分類できない。

## エネルギー法と Dirichlet の原理

:::proposition[Dirichlet の原理]{#prop-pde-dirichlet-energy}
$\Omega$ を有界連結な滑らかな領域とし、$u\in C^2(\overline\Omega)$ が $-\Delta u=f$、$u|_{\partial\Omega}=g$ を満たすとする。$E(v)=\int_\Omega(\frac12|Dv|^2-fv)\,\dd x$ とおけば、同じ境界値を持つ $C^1$ 関数について

$$
E(v)-E(u)=\frac12\int_\Omega|D(v-u)|^2\,\dd x .
\label{eq-pde-dirichlet-energy}
$$

従って $u$ は一意な最小化元である。逆に $C^2$ の最小化元は Poisson 方程式を満たす。
:::

:::proof
$w=v-u$ とおき展開すると、一次項は $\int(Du\cdot Dw-fw)$ である。$w=0$ が境界で成り立つため部分積分で一次項は $0$ になる。二次項が消える場合は $Dw=0$ なので連結性と境界値から $w=0$ である。逆に $u+s\phi$ を $\phi\in C_c^\infty(\Omega)$ とともに代入すると $s=0$ で一次微分が消える。部分積分で $\int(-\Delta u-f)\phi=0$ になり、連続な右辺なら点ごとの方程式が従う。
:::

この議論は滑らかな最小化元が存在するときの特徴付けである。最小化列の極限に必要な関数空間と、極限から PDE を得る方法は [楕円型方程式の弱解](/math/pde/weak-solutions/) で扱う。

:::example[線分上の最小化]
$[0,1]$ 上で $f=0$、端点値を $0,2$ とすると $u=2x$ である。$E(u)=2$ である。同じ端点値の $v=2x^2$ では $E(v)=8/3$ であり、差 $2/3$ は $\frac12\int_0^1(4x-2)^2\,\dd x$ に等しい。
:::

:::example[正の源と負の最小値]
半径 $1$ の球で $f=2n$、境界値 $0$ とすると $u=1-|x|^2$ である。$\int|Du|^2=2n\int u$ なので $E(u)=-n\int u<0$ である。境界と源を固定した最小化では、勾配の二乗だけでなく源との積も費用に含まれる。
:::

### 境界データが固定された場合の構成

$g$ がある $G\in H^1(\Omega)$ の跡であるとする。$f\in H^{-1}$ のときのエネルギーは $\norm{Du}_2^2/2-\langle f,u-G\rangle$ と書く。$H_0^1(\Omega)$ 上で $v=u-G$ を求める。すべての $\phi\in H_0^1$ について
$$
\int Dv\cdot D\phi=\langle f,\phi\rangle-\int DG\cdot D\phi
$$
を課す。右辺は勾配ノルムに関する有界線形汎関数であり、左辺は Poincaré の不等式により $H_0^1$ の完全な内積になる。Riesz の表現で $v$ が一意に得られ、エネルギー差は $\norm{D\phi}_2^2/2$ である。従って得た解はアフィン集合 $G+H_0^1$ の唯一の最小化元である。延長を $G+\psi$（$\psi\in H_0^1$）に取り替えても同じ集合と弱方程式になるため、$u$ は延長の選択によらない。

:::example[区間の端点と一定の源項]
$(0,2)$ で $f=3$、$u(0)=1$、$u(2)=5$ とする。アフィン延長 $G=1+2x$ に対し $\int G'\phi'=2[\phi]_0^2=0$ なので、補正は $-v''=3$ の零境界解 $v=3x(2-x)/2$ である。従って $u(1)=1+2+3/2=9/2$ である。任意の $\phi\in H_0^1(0,2)$ で $u+\phi$ のエネルギーとの差は $\int|\phi'|^2/2$ であり、非零の補正では正になる。
:::

:::example[定数外力のエネルギー最小化]
長さ $2$ の区間で $-u_{xx}=4$、$u(0)=0$、$u(2)=2$ とする。流束と値を求め、エネルギーの一次項が消えることを確認する。

まず流束 $q=Au_x$ を未知量にする。$q_x=-f$ を積分すると、任意の定数 $q(0)$ について $q(x)=q(0)-\int_0^x f$ である。さらに $q/A$ を積分し、指定された両端の値または平均に合わせて積分定数を決める。以下では計算後の $u$ について関係を逆に確かめる。

$$
\begin{aligned}
u(x)&=- 2 x^{2} + 5 x\\
u_x&=\frac{d}{d x} \left(- 2 x^{2} + 5 x\right)=5 - 4 x\\
q&=5 - 4 x\\
-f&=\frac{d}{d x} \left(5 - 4 x\right)=-4
\end{aligned}
$$

次に値の条件と流出量を両端で区別して評価する。区間の外法線は左端で $-1$、右端で $1$ である。このため源項の積分は左側の流束から右側の流束を引いたものになる。流束を係数で割った積分は両端の値の差に等しい。

$$
\begin{aligned}
u(0)&=0,\qquad x=0\\
u(L)&=2,\qquad x=2\\
q(0)&=5,\qquad x=0\\
q(L)&=-3,\qquad x=2\\
\int_0^L f\,\mathrm{d}x&=\int\limits_{0}^{2} 4\, dx=8\\
\int_0^L (q/A)\,\mathrm{d}x&=\int\limits_{0}^{2} \left(5 - 4 x\right)\, dx=2
\end{aligned}
$$

弱形式も具体的に比べる。試験関数 $\phi=x(L-x)$ は両端で零であり、積分の間に境界項は出ない。以下の値が一致するのは、任意の零トレース試験関数に対しても部分積分できることの一例である。

$$
\begin{aligned}
\phi&=x \left(2 - x\right)\\
\phi^{\prime}&=\frac{d}{d x} x \left(2 - x\right)=- 2 \left(x - 1\right)\\
\int_0^L q\phi^{\prime}\,\mathrm{d}x&=\int\limits_{0}^{2} \left(2 - 2 x\right) \left(5 - 4 x\right)\, dx=\frac{16}{3}\\
\int_0^L f\phi\,\mathrm{d}x&=\int\limits_{0}^{2} 4 x \left(2 - x\right)\, dx=\frac{16}{3}
\end{aligned}
$$

最後に内部での値、勾配、流束を同じ点で測る。区間全体での勾配の寄与には $A$ の重みが付くので、点での値だけでは全体のエネルギーは分からない。

$$
\begin{aligned}
u&=2,\qquad x=\frac{1}{2}\\
u_x&=3,\qquad x=\frac{1}{2}\\
q&=3,\qquad x=\frac{1}{2}\\
\frac12\int_0^L A(u_x)^2\,\mathrm{d}x&=\int\limits_{0}^{2} \frac{\left(5 - 4 x\right)^{2}}{2}\, dx=\frac{19}{3}
\end{aligned}
$$

この例では $J(u)=\int_0^2((u_x)^2/2-4u)=-37/3$ である。$u$ に $\eps\phi$ を加えると差は $\eps^2\int_0^2|\phi^{\prime}|^2/2=4\eps^2/3$ となる。従って、この方向では $\eps=0$ においてだけ最小となる。任意の零境界の方向でも同じ平方の形が成立する。
:::

## 演習問題

::::exercise{level=basic}
$u=x_1^2+x_2^2-2x_3^2$ が調和であるか判定し、中心 $(1,1,1)$ の球での平均を求めよ。

:::solution
$\Delta u=2+2-4=0$ であり、中心値は $1+1-2=0$ である。球体平均も球面平均も $0$ になる。
:::
::::

::::exercise{level=basic}
二次元の基本解について半径 $2$ の円を通る負の流束を求めよ。

:::solution
$\Gamma_r=-1/(2\pi r)$、円周長は $2\pi r$ なので $-\int\Gamma_r\,\dd S=1$ である。$r=2$ でも同じである。
:::
::::

::::exercise{level=basic}
単位円盤上の境界値 $x_1^2$ の調和拡張を求めよ。

:::solution
境界で $x_1^2=(1+x_1^2-x_2^2)/2$ である。右辺は内部で調和なのでこれが解であり、差の最大値原理で一意である。
:::
::::

::::exercise{level=basic}
$-\Delta u=6$、$u|_{\partial B_2}=0$ を三次元で解け。

:::solution
$\Delta |x|^2=6$ であるから $u=4-|x|^2$ である。正しい右辺と境界値を持ち、差が調和なので他の古典解と一致する。
:::
::::

::::exercise{level=advanced}
$\Omega=B_1(0)\subset\R^n$ で $u\in C^2(\overline\Omega)$、$\Delta u=0$、$\partial_\nu u=0$ なら定数であることを示せ。

:::solution
発散定理で $\int|Du|^2=\int_{\partial\Omega}u\partial_\nu u-\int u\Delta u=0$ である。$Du$ は連続なので各点で消え、球の連結性により定数である。従って Neumann データでは定数分の不定性が残る。
:::
::::

::::exercise{level=advanced}
全空間で調和で、$|u(x)|\le C(1+|x|^{3/2})$ とする。$u$ はアフィンであることを示せ。

:::solution
次数 $2$ の内部評価を任意の中心で半径 $R$ に対して使うと $|D^2u(x)|\le C_x R^{-1/2}$（$R\ge1$）となる。$R\to\infty$ で Hessian が消えるので各一次微分は定数であり、主張が従う。
:::
::::

::::exercise{level=advanced}
単位円盤で $y=(1,0)$ と $x=(1/2,0)$ の Poisson 核を求め、中心値での核との比を求めよ。

:::solution
:eqref[eq-pde-ball-poisson] で分子 $3/4$、分母 $2\pi\cdot1/4$ なので核は $3/(2\pi)$ である。中心の核は $1/(2\pi)$ なので比は $3$ である。
:::
::::

::::exercise{level=advanced}
区間 $[0,1]$ で $-u''=2$、$u(0)=u(1)=0$ を解き、最小エネルギーを求めよ。

:::solution
二度積分して $u=x(1-x)$ を得る。$\int_0^1(u')^2\,\dd x=1/3$、$\int_0^1 2u\,\dd x=1/3$ なので $E(u)=1/6-1/3=-1/6$ である。:ref[prop-pde-dirichlet-energy] により最小値である。
:::
::::

## 確認テスト

::::quiz
本章の基本解の符号規約はどれか。

- [ ] $\Delta\Gamma=\delta_0$
- [x] $-\Delta\Gamma=\delta_0$
- [ ] $\Delta\Gamma=1$

:::explanation
試験関数と組み合わせると、原点を囲む球面で $-\partial_r\Gamma$ の積分が $1$ になる。これが $-\Delta\Gamma=\delta_0$ の規格化である。
:::
::::

::::quiz{answer="3" tol="1e-6"}
調和関数 $u=x_1^2-x_2^2$ を中心 $(2,1)$ の円周で平均した値を求めよ。

:::explanation
二階微分の和は $2-2=0$ なので平均値公式を適用できる。中心値 $2^2-1^2=3$ が答えである。
:::
::::

::::quiz
有界な連結領域で、境界まで連続な調和関数が内点に最大値を持つ。この関数について何が従うか。

- [x] 領域全体で定数である
- [ ] 境界で必ず $0$ である
- [ ] その一点以外は任意に選べる

:::explanation
最大点を中心とした球の平均が最大値に等しいため、球内の全ての値が最大値になる。最大点集合は開かつ閉なので連結領域全体に一致する。
:::
::::

::::quiz{answer="0.5" tol="1e-6"}
上半平面で境界の $[-1,1]$ に値 $1$、外側に値 $0$ を置いた Poisson 積分を $(0,1)$ で評価せよ。

:::explanation
値は $\pi^{-1}\int_{-1}^1(1+s^2)^{-1}\,\dd s=2\arctan 1/\pi=1/2$ である。
:::
::::

::::quiz{answer="5" tol="1e-6"}
単位円板における調和関数の境界値を $g(\theta)=2+3\cos(2\theta)$ とする。点 $(1,0)$ における境界値を $A$、中心の値を $B$ とし、$A$ を求めよ。

:::explanation
$r=1,\theta=0$ で $\cos(2\theta)=1$ なので $A=5$ である。なお内部調和拡張は $2+3r^2\cos(2\theta)$ で中心の値は $2$ になる。平均値の等式は中心の値を決めるが、それを境界の一点の値として使うことはできない。
:::
::::


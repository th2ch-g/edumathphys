---
title: "共役関数と Hilbert 変換"
description: "円周上の共役関数を核と乗数で記述し、強型・弱型評価から Fourier 部分和の収束を導く。"
notation:
  - symbol: 'Q_r'
    meaning: "共役 Poisson 核"
  - symbol: 'H'
    meaning: "円周の Hilbert 変換"
  - symbol: 'P_{\ge0}'
    meaning: "非負の周波数への射影"
  - symbol: 'S_N'
    meaning: "対称 Fourier 部分和"
---

## この章の目標

- 正則関数の共役調和成分から Hilbert 変換を導く。
- Fourier 乗数と主値核による二つの表示を証明する。
- $L^p$ 有界性と $L^1$ の弱型評価を区別する。
- 周波数への射影で $L^p$ ノルム収束を証明し、上半平面の Hardy 空間と結び付ける。

## 共役関数と共役 Poisson 核

$dm=d\theta/(2\pi)$ とする。実数値の境界関数 $f$ の Poisson 積分を正則関数の実部にするには、負の周波数を正の周波数の共役としてまとめればよい。虚部の定数は円板の中心でゼロに選ぶ。

:::proposition[正則関数と共役 Poisson 核]{#prop-four-conjugate-poisson}
$f\in L^1(\mathbb T)$ が実数値なら

$$
F(z)=\int_{\mathbb T}\frac{e^{it}+z}{e^{it}-z}f(t)\,dm(t)
=\widehat f(0)+2\sum_{k\ge1}\widehat f(k)z^k
$$

は円板内で正則で、$z=re^{ix}$ における実部は $P_r*f(x)$、虚部は $Q_r*f(x)$ である。ここで :term[共役 Poisson 核]{en="conjugate Poisson kernel" yomi="きょうやくぽあそんかく"} は

$$
Q_r(t)=\frac{2r\sin t}{1-2r\cos t+r^2},
\qquad \widehat Q_r(k)=-i\sgn(k)r^{|k|}.
$$
:::

:::proof
商を等比級数で展開し、小円板上で一様収束することから項別積分する。$z/e^{it}=re^{i(x-t)}$ を代入し、分母の共役を掛けると実部が $P_r(x-t)$、虚部が $Q_r(x-t)$ になる。虚部は $2\sum_{k\ge1}r^k\sin kt$ で、正負の指数に分ければ係数の式を得る。
:::

:::example[一つの共役モード]
$f=5+2\cos(3x)$ では $F(z)=5+2z^3$ であり、虚部は $2r^3\sin(3x)$ である。定数成分は虚部に寄与しない。$f=\sin(3x)$ に対しては $F(z)=-iz^3$ なので虚部は $-r^3\cos(3x)$ となる。
:::

:::example[原点で交差する核]
$Q_r$ は奇関数で積分はゼロである。$x=\pi/2$ で $Q_{1/2}(x)=4/5$ だが、$Q_{1/2}(-x)=-4/5$ となる。正の核による平均の考え方はそのまま適用できない。
:::

### 境界の実部と虚部を相互に決める条件

正則関数 $F=u+iv$ が与えられると、$u,v$ は Cauchy–Riemann 方程式で結ばれる。従って実部を固定すると虚部は加法定数を除いて決まる。円板で境界から $u$ を作るときには、中心での $v$ を指定すれば境界の共役を一意に選べる。

:::proposition[共役調和関数の一意性と極座標微分]{#prop-four-circle-conjugate-cr}
単位円板内で正則な $F=u+iv$ について、$r>0$ なら
$$
r\partial_r u=\partial_x v,\qquad r\partial_r v=-\partial_x u
$$
である。二つの正則関数の実部が一致すれば虚部の差は実定数である。従って $v(0)=0$ は共役を一意にする。
:::

:::proof
$z=re^{ix}$ とすると $\partial_r F=e^{ix}F'(z)$、$\partial_xF=ir e^{ix}F'(z)$ である。$\partial_xF=ir\partial_rF$ の実部と虚部を比較すれば二式を得る。実部の差がゼロなら複素微分可能性に伴う直交座標での Cauchy–Riemann 方程式により虚部の二つの一次微分はゼロである。円板は連結なので虚部の差は定数である。
:::

:::example[円板の中心での定数を修正する]
$F(z)=3z+2iz^2+4i$ の実部は $3r\cos x-2r^2\sin(2x)$ である。虚部は $3r\sin x+2r^2\cos(2x)+4$ だが、共役 Poisson 核で得る虚部には最後の $4$ がない。実際、$F(0)=4i$ であり、:ref[prop-four-conjugate-poisson] の積分では $F(0)$ は実数である。境界では $f=3\cos x-2\sin2x$ から $3\sin x+2\cos2x$ が得られる。原点での正規化を変えても、円周での非零周波数には影響しない。
:::

### 正の実部を持つ正則関数の測度表示

積分核 $\frac{e^{it}+z}{e^{it}-z}$ の実部は正である。したがって正測度を与えると正の実部を持つ正則関数を得る。逆に、実部が正であれば各半径におけるその積分の質量が中心値に等しいため、境界で弱い極限を取れる。

:::theorem[Herglotz の表示]{#thm-four-herglotz}
単位円板上で正則な $F$ が $\Re F\ge0$ を満たすとする。このとき一意な有限正測度 $\mu$ と実数 $b$ があって
$$
F(z)=ib+\int_{\mathbb T}\frac{e^{it}+z}{e^{it}-z}\,d\mu(t).
$$
このとき $\mu(\mathbb T)=\Re F(0)$、$b=\Im F(0)$ である。
:::

:::proof
$r<1$ に対し $d\mu_r(t)=\Re F(re^{it})dm(t)$ と置く。円板における平均値公式または Taylor 級数の定数項から $\mu_r(\mathbb T)=\Re F(0)$ である。円周上の連続関数空間は可分なので、$r_j\uparrow1$ から可算稠密な関数に対して積分が収束する部分列を取る。一様な質量上界で全ての連続関数に収束を拡張でき、極限は正の有界線形汎関数である。Riesz 表現により正測度 $\mu$ を得る。

固定 $z$ で $r>|z|$ なら、$w\mapsto F(rw)$ の実部の境界値は $\Re F(re^{it})$ である。この正則関数の Taylor 級数とその共役を展開すると、実部の境界値から :ref[prop-four-conjugate-poisson] の積分で元の関数を虚定数を除いて復元できる。従って
$$
F(z)-i\Im F(0)=\int_{\mathbb T}
   \frac{e^{it}+z/r}{e^{it}-z/r}\,d\mu_r(t).
$$
$r\to1$ で被積分関数は $t$ に一様に収束する。弱い極限と一様な質量上界を使えば表示が得られる。

一意性については、積分を Taylor 展開すると定数項は $\mu(\mathbb T)$、次数 $k\ge1$ の係数は $2\int e^{-ikt}d\mu(t)$ である。正測度なので負の係数はこれらの共役で決まり、全三角多項式上の積分が定まる。三角多項式は :ref[thm-four-fejer] により連続関数に一様稠密であるから、有限測度の一意性が従う。
:::

:::example[離れた二点の質量から有理関数を作る]
$\mu=\delta_0+2\delta_\pi$、$b=0$ とすると
$$
F(z)=\frac{1+z}{1-z}+2\frac{1-z}{1+z}.
$$
$F(0)=3$ である。$z=1/2$ では最初の項が $3$、次が $2/3$ なので $F(1/2)=11/3$ となる。半径 $r$ の実部は $\theta=0$ と $\theta=\pi$ を中心とする Poisson 核の和である。ほとんど全ての境界点では各核の実部がゼロへ近づく一方、各半径での積分は常に $3$ である。点ごとの境界値だけでは、表示に現れる測度の質量は回収できない。
:::

### 共役核の積分ノルム

$Q_r$ の積分はゼロである。だが奇関数なので、絶対値の積分では正負を足すことになる。半径を $1$ に近づけたときの符号付きの打ち消しと、絶対値の大きさを区別する。

:::proposition[共役 Poisson 核の $L^1$ ノルム]{#prop-four-conjugate-poisson-l1}
$0<r<1$ では
$$
\|Q_r\|_1=\frac2\pi\log\frac{1+r}{1-r}.
$$
:::

:::proof
$Q_r(t)$ は $0<t<\pi$ で正で、その分子 $2r\sin t$ は分母 $1-2r\cos t+r^2$ の微分である。負の半周では奇性で絶対積分が同じだから
$$
\|Q_r\|_1=\frac1\pi[\log(1-2r\cos t+r^2)]_0^\pi
=\frac1\pi\log\frac{(1+r)^2}{(1-r)^2}.
$$
:::

:::example[半径を一段上げるときの絶対積分]
$r=1/2$ で $\|Q_r\|_1=2\log3/\pi$、$r=3/4$ で $2\log7/\pi$ である。$Q_r$ が最高値に達する正の角は、分子と分母を微分して $\cos t=2r/(1+r^2)$ を解けば求まる。例えば $r=1/2$ の角では $\cos t=4/5$、$\sin t=3/5$ で最大値は $4/3$ である。$r$ を $1$ に近づけると角が原点に集まり、核の両側の絶対積分が増える。
:::

## Hilbert 変換の L² 理論

:::definition[円周の Hilbert 変換]{#def-four-hilbert-circle}
$f\in L^2(\mathbb T)$ に対し、$\widehat{Hf}(k)=-i\sgn(k)\widehat f(k)$ で定まる $L^2$ 関数を :term[Hilbert 変換]{en="Hilbert transform" yomi="ひるべるとへんかん"} という。$\sgn(0)=0$ である。実数値の $f$ について $Hf$ を共役関数とも呼ぶ。
:::

:::proposition[二乗ノルムと主値]{#prop-four-hilbert-l2}
$H^*=-H$、$H^2f=-f+\widehat f(0)$、$\|Hf\|_2^2=\|f\|_2^2-|\widehat f(0)|^2$ であり、$P_r*Hf=Q_r*f\to Hf$ が $L^2$ で成り立つ。また $f\in C^1$ なら

$$
Hf(x)=\lim_{\eps\downarrow0}\int_{\eps<|t|\le\pi} f(x-t)\cot(t/2)\,dm(t).
\label{eq-four-circle-hilbert-pv}
$$

この切断は任意の $f\in L^2$ についても $L^2$ で $Hf$ に収束する。
:::

:::proof
:ref[thm-four-circle-riesz-fischer] で定義が可能であり、Parseval により随伴は乗数の共役、二乗は乗数の二乗で求まる。$P_r$ と $Q_r$ の係数および :ref[thm-four-good-kernel] から径方向の収束が従う。

$t\ne0$ では $Q_r(t)\to\cot(t/2)$ であり、$r\ge1/2$ で $|Q_r(t)|\le C/|t|$ である。奇性で $f(x)$ を積分から引けるため、$C^1$ では差を $O(|t|)$ で抑えて優収束を使える。主値切断の乗数は奇性により $(-i/\pi)\int_\eps^\pi\sin(kt)\cot(t/2)dt$ である。$\cot(t/2)-2/t$ は可積分なのでその寄与は一様有界であり、$2/t$ の寄与は $u=|k|t$ の置換と :ref[lem-four-series-sine-primitive] で一様に抑えられる。各 $k$ の乗数の極限は $f=e_k$ の場合に先の議論で求めた。係数平方和へ優収束を適用すれば結論を得る。
:::

:::example[定数を残す二重変換]
$f=2+\cos(2x)-3\sin x$ なら $Hf=\sin(2x)+3\cos x$ であり、$H^2 f=-\cos(2x)+3\sin x$ である。差 $f+H^2 f$ は平均 $2$ だけを持つ。
:::

:::example[二つのモードの主値を原点で評価する]
$f(t)=\cos t+\sin(2t)$ とする。$f(-t)=\cos t-\sin(2t)$ だから、奇核との積分で $\cos t$ の項は消える。残る積分は $-(1/\pi)\int_0^\pi\sin(2t)\cot(t/2)dt$ である。$\sin(2t)\cot(t/2)=2\cos t(1+\cos t)$ なので積分値は $-1$ になる。乗数で計算した $Hf=\sin t-\cos(2t)$ の原点の値も $-1$ である。
:::

### 実軸での $L^2$ 主値

円周では定数モードが残った。実軸では周波数 $0$ は Lebesgue 零集合であり、二乗可積分な非零の定数関数も存在しない。このため実軸の Hilbert 変換は全 $L^2$ 上で等長となる。

:::theorem[実軸の Hilbert 変換の主値表示]{#thm-four-real-hilbert-l2}
$H_{\R}$ を乗数 $-i\sgn(\xi)$ の $L^2(\R)$ 作用素とする。$0<\eps<R$ の切断
$$
H_{\eps,R}^{\R}f(x)=\frac1\pi\int_{\eps<|t|<R}\frac{f(x-t)}t dt
$$
は $\eps\downarrow0, R\uparrow\infty$ で $H_\R f$ に $L^2$ 収束する。各 $f\in L^2$ では $|t|>\eps$ の積分を $R=\infty$ として定義でき、その $L^2$ 極限も同じである。また $H_\R^*=-H_\R$、$H_\R^2=-I$、$\|H_\R f\|_2=\|f\|_2$ である。
:::

:::proof
有限切断核は $L^1$ であり、奇性から乗数は
$$
m_{\eps,R}(\xi)=-\frac{2i}\pi\int_\eps^R\frac{\sin(2\pi t\xi)}t dt
$$
である。$\xi=0$ ではゼロ、$\xi\ne0$ では $u=2\pi|\xi|t$ と置く。$u<1$ では $|\sin u|/u\le1$、$u>1$ では部分積分で $\int_a^b\sin u/u$ の絶対値が $3/a$ 以下となる。従って全切断で乗数は一様有界である。極限は :ref[prop-four-sine-integral-limit] により $-i\sgn(\xi)$ であり、Plancherel と優収束から $L^2$ 強収束を得る。

固定 $\eps$ に対し $|t|>\eps$ では $t^{-1}\in L^2$ なので、$f$ との積の絶対積分は Cauchy–Schwarz で有限である。また切断の $R\to\infty$ の乗数極限は存在して一様有界だから $L^2$ 極限があり、a.e. 収束部分列で点ごとの積分と同定できる。この乗数で $\eps\downarrow0$ とすれば先と同じ極限となる。三つの作用素の関係は、乗数の共役、二乗、絶対値をそれぞれ計算して得る。
:::

:::example[実軸の区間指示関数]
$f=1_{[-1,1]}$ とすると、$x\ne\pm1$ で
$$
H_\R f(x)=\frac1\pi\log\left|\frac{x+1}{x-1}\right|.
$$
実際、積分を $t=x-y$ へ移すと範囲は $x-1<t<x+1$ である。$|x|<1$ の場合は原点の両側から $\eps$ までを除き、二つの $\log\eps$ が消える。$x=0$ はゼロ、$x=2$ では $\log3/\pi$ である。点ごとの切断の極限が存在し、:ref[thm-four-real-hilbert-l2] の $L^2$ 極限の部分列から乗数による変換とも一致する。
:::

### 共役で作る周波数射影

:::proposition[正周波数への直交射影]{#prop-four-hilbert-orthogonal-projections}
円周上の $L^2$ では $P_+=(I+iH-P_0)/2$、$P_-=(I-iH-P_0)/2$ がそれぞれ正、負の周波数への直交射影である。$P_0f=\widehat f(0)$ とする。実軸上では $P_\pm^\R=(I\pm iH_\R)/2$ が正負の半軸への直交射影である。
:::

:::proof
円周では $k>0$ のとき $iH$ の乗数は $1$、$k<0$ で $-1$、$k=0$ で $0$ である。$P_0$ の乗数は原点だけで $1$ だから、$P_\pm$ の乗数はそれぞれの半軸の指示関数になる。指示関数の乗数は値が実で二乗して自分自身であるため、Parseval により自己随伴かつ冪等である。実軸では Plancherel と $0$ が零集合であることを使って同様に従う。
:::

:::example[実信号から正周波数を作る]
$f(x)=3+\cos(2x)+\sin x$ とする。$Hf=\sin(2x)-\cos x$ なので $f+iHf=3+e^{2ix}-ie^{ix}$ である。$(f+iHf-P_0 f)/2=(e^{2ix}-ie^{ix})/2$ は $f$ の正の Fourier 係数だけを持つ。$\|P_+f\|_2^2=1/2$、$\|P_-f\|_2^2=1/2$、$\|f\|_2^2=10$ であり、残りの $9$ は定数成分の二乗である。共役で作った複素信号の正周波数の係数は、元の正周波数の係数の二倍となる。
:::

### 一次微分と符号乗数の合成

:::proposition[絶対周波数の取り出し]{#prop-four-hilbert-absolute-derivative}
三角多項式 $f$ について $H(\partial_x f)$ の係数は $|k|\widehat f(k)$ である。実軸で $f\in\mathcal S$ なら $\widehat{H_\R(\partial_x f)}(\xi)=2\pi|\xi|\widehat f(\xi)$ であり、:ref[prop-four-half-plane-dirichlet-neumann] の法線微分と同じである。
:::

:::proof
円周では微分の乗数 $ik$ に $-i\sgn(k)$ を掛け、積は $|k|$ になる。実軸では微分の乗数が $2\pi i\xi$ に変わるため積は $2\pi|\xi|$ になる。いずれも二つの乗数の積で順序を交換できる。実軸の法線微分は同じ乗数を持ち、Plancherel により二つの $L^2$ 関数は一致する。
:::

:::example[一方向の位相微分と絶対周波数]
$f=e^{-3ix}+2e^{ix}$ では $\partial_x f=-3ie^{-3ix}+2ie^{ix}$ である。負の成分に $i$、正の成分に $-i$ を掛けると $Hf'=3e^{-3ix}+2e^{ix}$ となる。普通の微分では符号が変わる二つのモードが、合成後には両方とも正の重みを受ける。
:::

## M. Riesz の定理

二乗ノルムの評価だけでは他の $p$ は分からない。共役調和関数を正則関数として組み合わせると、正則関数の二乗に対しても実部と虚部が再び共役になる。この関係で指数を倍にできる。

:::theorem[M. Riesz]{#thm-four-riesz-hilbert}
$1<p<\infty$ に対し $H$ は $L^p(\mathbb T)$ に有界に拡張される。各 $f\in L^p$ で $Q_r*f\to Hf$ が $L^p$ で成り立つ。
:::

:::proof
実三角多項式 $f$ について $f+iHf$ は正の周波数と実定数係数のみを持つ。その二乗の実部は $f^2-(Hf)^2$、虚部は $2fHf$ で、定数係数は $\widehat f(0)^2$ である。よって $H(2fHf)=-(f^2-(Hf)^2)+\widehat f(0)^2$ である。$H$ が $L^p$ で定数 $C_p$ を持つと仮定し、この式を $(Hf)^2$ について解くと

$$
\|Hf\|_{2p}^2\le 2C_p\|f\|_{2p}\|Hf\|_{2p}+2\|f\|_{2p}^2.
$$

積は Hölder、定数は $|\widehat f(0)|\le\|f\|_{2p}$ で評価した。二次不等式により $\|Hf\|_{2p}\le(C_p+\sqrt{C_p^2+2})\|f\|_{2p}$ を得る。複素多項式は実部と虚部を分けて定数を高々二倍にすればよい。

$p=2$ から反復し、:ref[thm-four-fejer] の稠密性で各 $p=2^j$ に拡張する。拡張は共通定義域で係数が同じなので :ref[thm-four-series-uniqueness] により一致する。:ref[thm-meas-riesz-thorin] で $2\le p<\infty$ の間を埋める。$1<p<2$ では、多項式上で $\int Hf\,\bar g=-\int f\,\overline{Hg}$ を使い、$p'>2$ の評価を適用する。多項式 $g$ は $L^{p'}$ で稠密なので双対ノルムの上限を取ることができる。これで $L^p$ 評価が従い、拡張できる。最後に $P_r*Hf$ と $Q_r*f$ の係数は一致し、:ref[thm-four-good-kernel] により $L^p$ 極限を得る。
:::

:::example[$p=4$ の計算]
$p=2$ の定数 $1$ を上の二次不等式に入れると実多項式について $C_4\le1+\sqrt3$ を得る。例えば $f=\cos t$ なら $Hf=\sin t$ で両者の四乗平均は $3/8$ である。単一モードの比 $1$ と一般の上界は異なる。
:::

### 主値と共役 Poisson 平均の比較

$L^p$ の極限だけでは、切断半径を任意に小さくしたときに各点で極限があるか分からない。主値切断と同じ尺度の共役 Poisson 平均との差を、局所の極大関数で抑える。

:::proposition[円周の主値の概収束]{#prop-four-circle-hilbert-ae}
$1<p<\infty$ なら :eqref[eq-four-circle-hilbert-pv] の右辺はほとんど至る所で収束し、その値は $Hf$ である。
:::

:::proof
円周上の中心付き極大関数 $M_{\mathbb T}f$ を弧の絶対平均の上限とする。これは周期関数を少なくとも三周期だけ残して実軸上の $M$ に掛ければ、中心を一周期内に制限したとき上から抑えられる。従って $L^p$ で有界である。実軸での有界性は :ref[thm-meas-maximal-weak] と $L^\infty$ 評価の補間に従う。

$0<\eps<1/2$ とし $r=1-\eps$ と置く。$|t|\le\eps$ で $|Q_r(t)|\le C|t|/\eps^2\le C/\eps$ である。一方、$|t|>\eps$ では
$$
\cot(t/2)-Q_r(t)=\cot(t/2)\frac{(1-r)^2}{(1-r)^2+2r(1-\cos t)}
$$
で、絶対値は $C\eps^2|t|^{-3}$ 以下である。内側と外側の二進環を積分すると、主値切断と $Q_r*f$ の差の絶対値は $C M_{\mathbb T}f$ 以下である。さらに $Q_r*f=P_r*Hf$ であり、$P_r(t)\le C\eps/(\eps^2+t^2)$ を環で積分すると $\sup_{r\ge1/2}|P_r*Hf|\le C M_{\mathbb T}(Hf)$ である。このため切断の極大関数は $L^p$ で有界となる。三角多項式では主値の極限は :ref[prop-four-hilbert-l2] で一致する。多項式による $L^p$ 近似を使い、:ref[thm-four-maximal-convergence] の誤差の証明で極限を $Hf$ へ拡張する。
:::

:::example[共役核の尺度を読み取る]
$t=\eps$ に置くと $Q_{1-\eps}(\eps)\sim1/\eps$ である。実際、分子は $2\eps+O(\eps^2)$、分母は $2\eps^2+O(\eps^3)$ である。一方で $\cot(\eps/2)\sim2/\eps$ であり、同じ尺度でも値は一致しない。その差はこの証明では、原点からの距離を含む積分評価で扱った。
:::

### 円周を拡大して実軸へ移す

実軸の $L^2$ 乗数を得たが、他の指数での評価を改めて複素関数から証明する必要はない。大きい円周の共役核が、固定した有界範囲で実軸の核に近づくことを利用する。

:::theorem[実軸の M. Riesz 評価]{#thm-four-real-hilbert-lp}
$1<p<\infty$ で $H_\R$ は $L^p(\R)$ に有界に拡張する。円周の評価 $\|Hf\|_p\le C_p\|f\|_p$ と同じ定数を上界として使える。
:::

:::proof
$f\in C_c^\infty(\R)$ を取り、台を $[-A,A]$ に含める。$L>2A$ に対し $2L$ 周期へ延長し、$x\mapsto \pi x/L$ の角度変換で円周上の作用素を適用する。$x$ を固定し $L>|x|+A$ とすれば、積分で被積分関数が非零となるのは同じ一周期の $y\in[-A,A]$ だけである。このとき変換は
$$
H_L f(x)=\frac1{2L}\operatorname{pv}\int_{-A}^A f(y)\cot\left(\frac{\pi(x-y)}{2L}\right)dy.
$$
$|u|\le\pi/4$ で $\cot u-1/u=O(u)$ である。固定 $x$ で $L$ が十分大きければ、核の差 $(2L)^{-1}\cot(\pi(x-y)/(2L))-1/(\pi(x-y))$ は $C|x-y|/L^2$ 以下であり、原点で連続に延長できる。従って $H_L f(x)\to H_\R f(x)$ となる。

円周の $L^p$ 評価では、入力と出力の積分の両方に付く $1/(2L)$ が相殺する。従って $\int_{-L}^L|H_Lf|^pdx\le C_p^p\int_\R|f|^pdx$ である。$L\to\infty$ の列で $H_L f\,1_{[-L,L]}$ に Fatou を適用すれば実軸の評価になる。$C_c^\infty$ の稠密性から全 $L^p$ に拡張する。同時に $L^2$ に属する入力では $L^2$ と $L^p$ の双方で近似できるため、測度での極限の一意性により拡張は一致する。
:::

:::example[周期を広げたときの一点の値]
実軸の $1_{[-1,1]}$ を長さ $2L$ で周期化し、$x=2$（$L>3$）で変換する。そこでは主値を取る必要がなく、
$$
H_L f(2)=\frac1\pi\log\frac{\sin(3\pi/(2L))}{\sin(\pi/(2L))}.
$$
$L=6$ では正弦の比は $\sin(\pi/4)/\sin(\pi/12)=1+\sqrt3$ である。$L\to\infty$ では $\sin v/v\to1$ により比は $3$ となり、実軸の $\log3/\pi$ に戻る。台の近くを固定したまま周期を大きくすると、核の違いが小さくなる。
:::

### 四乗の平均から得られる追加評価

:::proposition[実多項式の四乗平均]{#prop-four-hilbert-fourth-moment}
実三角多項式 $f$ の平均を $m$、$g=Hf$ とする。このとき
$$
\|f\|_4^4+\|g\|_4^4=6\int_{\mathbb T}f^2g^2dm+m^4.
$$
従って $\|g\|_4\le\sqrt6\,\|f\|_4$ であり、$m=0$ の場合は $\|g\|_4\le(1+\sqrt2)\|f\|_4$ である。
:::

:::proof
$f+ig$ は負の周波数を持たず、定数係数が $m$ である。これを四乗しても定数係数は $m^4$ であり、他の項は正の周波数を持つ。平均を取り、実部を二項展開すると等式を得る。$|m|\le\|f\|_4$ と Hölder より
$$
\|g\|_4^4\le 6\|f\|_4^2\|g\|_4^2
$$
だから前半の評価となる。$m=0$、$\|f\|_4>0$ なら $r=\|g\|_4^2/\|f\|_4^2$ が $1+r^2\le6r$ を満たす。上の根は $3+2\sqrt2=(1+\sqrt2)^2$ であり、平方根を取って結論を得る。$f=0$ は自明である。
:::

:::example[正弦と余弦を同じ実部に入れる]
$f(x)=\cos x+\sin2x$、$g(x)=\sin x-\cos2x$ とする。$f^4$ を展開すると $\cos^4x$ と $\sin^4(2x)$ の平均は各 $3/8$、奇数個を含む混合項の平均はゼロである。残る $\cos^2x\sin^2(2x)$ の平均は $1/4$ なので $\|f\|_4^4=9/4$ となる。同じように $\|g\|_4^4=9/4$ であり、恒等式から $\int f^2g^2 dm=3/4$ である。例えば $\cos^3x\sin2x$ は奇関数であり、$\sin^3x\cos2x$ も奇関数なので消える。もう一方の奇数項は周波数の直交性でも消える。平均の反復積には、各周波数を単独で評価するだけでは分からない混合項が残る。
:::

## Kolmogorov の弱型評価

$L^1$ では入力が一箇所に集中すると出力に $1/x$ 型の尾ができる。尾の絶対積分を制御する代わりに、出力が閾値を超える集合の測度を調べる。

:::theorem[Kolmogorov の弱型評価]{#thm-four-kolmogorov}
ある絶対定数 $C$ により

$$
m\{|Hf|>\lambda\}\le C\|f\|_1/\lambda\qquad(\lambda>0)
$$

が $f\in L^1\cap L^2$ で成り立つ。これにより $H$ は $L^1$ 全体に測度収束の意味で拡張できる。
:::

:::proof
$\lambda<\|f\|_1$ なら測度が高々 $1$ なのでよい。そうでなければ :ref[thm-four-cz-decomposition] と同じ手順で、円周の弧を二分し続ける dyadic 分解を取る。全周の平均は閾値以下なので各選択弧 $I$ の親の平均が $\lambda$ 以下である。従って $f=g+\sum b_I$、$\|g\|_2^2\le2\lambda\|f\|_1$、$\int b_I=0$、$\sum\|b_I\|_1\le2\|f\|_1$ を得る。$Hg$ は二乗評価と Chebyshev で処理する。三倍に拡大した選択弧の和集合 $E$ の測度は $3\|f\|_1/\lambda$ 以下である。

$I$ の中心を $c_I$ とし $K(t)=\cot(t/2)$ とする。$x\notin 3I$ では、主値切断の $L^2$ 極限から

$$
Hb_I(x)=\int_I [K(x-y)-K(x-c_I)]b_I(y)\,dm(y)
$$

が a.e. 成り立つ。円周距離を $d$ とすれば $|K'(t)|\le C/d(t,0)^2$ である。従って差は $C|I|/d(x,c_I)^2$ 以下で、$x\notin3I$ の絶対積分は $C'\|b_I\|_1$ 以下である。弧を足す前に有限部分和を考える。$f\in L^2$ なら $b$ の和は互いに素な台上で $L^2$ 収束し、$H$ を掛けた和も収束する。その a.e. 収束部分列に Fatou を使えば $\int_{E^c}|Hb|\le C'\sum\|b_I\|_1$ である。これに $L^1$ の Chebyshev を適用し、$Hg$ と $E$ を合わせる。

$L^1$ 近似列に評価を差に適用すると出力は測度で Cauchy になる。閾値 $2^{-j}$ に対し隣接差がそれを超える集合の測度が $2^{-j}$ 以下の部分列を選べば、a.e. で出力の差が絶対可和になり極限が存在する。全列は同じ極限に測度収束し、近似列の選択にも依存しない。閾値を少し下げた評価に a.e. 収束部分列を使い、その後閾値を戻せば拡張にも評価が成り立つ。
:::

:::example[$L^1$ の強型が失われる例]
正の Fejér 核 $F_N$ は $\|F_N\|_1=1$ である。ところが $HF_N=N^{-1}\sum_{j=0}^{N-1}2\sum_{k=1}^j\sin kt$ であり、$t\ne0$ では等比級数から

$$
HF_N(t)=\cot(t/2)-\frac{\sin(Nt)}{2N\sin^2(t/2)}
$$

となる。実際 $\sum_{k=1}^j\sin kt=(\cos(t/2)-\cos((j+1/2)t))/(2\sin(t/2))$ と $\sum_{j=0}^{N-1}\cos((j+1/2)t)=\sin(Nt)/(2\sin(t/2))$ を使った。固定 $t\ne0$ で極限は $\cot(t/2)$ であり、その絶対積分は無限である。Fatou により $\|HF_N\|_1$ は無限大に向かう。
:::

:::example[円周の小弧で起こる対数]
$f=1_{[0,\pi/2]}$ とする。弧の外側では核に特異点がないため
$$
Hf(x)=\frac1\pi\log\left|\frac{\sin(x/2)}{\sin((x-\pi/2)/2)}\right|.
$$
これは $\cot(t/2)$ の原始関数 $2\log|\sin(t/2)|$ で積分した結果である。$x=\pi$ で値は $(\log 2)/(2\pi)$、$x\to(\pi/2)+$ では対数的に発散する。端点近くで対数は可積分なので、この例そのものでは強型 $L^1$ の失敗は現れない。
:::

### 切断半径に一様な弱型評価

変換そのものを弱型評価しても、任意の半径での切断を同時には制御できない。そこで最初に、切断を共役 Poisson 核に替える。この核の導関数は半径に一様な評価を持つため、悪い部分の平均ゼロを半径全体に一度に使える。

:::theorem[共役 Poisson 極大関数と主値の $L^1$ 端点]{#thm-four-hilbert-maximal-weak}
$f\in L^1(\mathbb T)$、$\lambda>0$ とする。
$$
m\left\{\sup_{0<r<1}|Q_r*f|>\lambda\right\}
+m\left\{\sup_{0<\eps<\pi}|H_\eps f|>\lambda\right\}
\le C\|f\|_1/\lambda ,
$$
ここで $H_\eps f$ は :eqref[eq-four-circle-hilbert-pv] の切断積分である。$H_\eps f$ は $\eps\downarrow0$ で a.e. に存在し、:ref[thm-four-kolmogorov] で得た $Hf$ に一致する。$Q_r*f$ も $r\uparrow1$ で a.e. に $Hf$ に収束する。
:::

:::proof
各固定 $r<1$ の畳み込みは有界で滑らかな核で定義される。$r$ について連続なので上限は可算添字で取れる。$r\ge1/2$ では
$$
|\partial_t Q_r(t)|\le
C\left(\frac1{(1-r)^2+t^2}
+\frac{t^2}{((1-r)^2+t^2)^2}\right)
\le C/|t|^2
\quad(0<|t|\le\pi).
$$
分母を $(1-r)^2+2r(1-\cos t)$ と書き、$1-\cos t\ge c t^2$ と $|\sin t|\le|t|$ を用いた。$r<1/2$ は分母に下界があるので同じ右辺で抑えられる。円周上で $|t|$ は原点までの最短距離ともみなせる。

$\lambda\ge\|f\|_1$ として :ref[thm-four-kolmogorov] の証明と同じ二進分解をとる。$\lambda<\|f\|_1$ は全体の測度 $1$ で評価できる。良い部分 $g$ は $\|g\|_2^2\le2\lambda\|f\|_1$ である。$Q_r*g=P_r*(Hg)$ が成り立ち、Poisson 核の二進環評価により $\sup_r|P_r*(Hg)|\le C M_{\mathbb T}(Hg)$ である。$r<1/2$ でも核の有界性と全周での平均を使えば同様に成り立つ。従って極大関数の $L^2$ 評価と $H$ の $L^2$ 評価から、良い部分が $\lambda/2$ を超える集合の測度は $C\|f\|_1/\lambda$ 以下である。

各悪い部分 $b_I$ では $x\notin3I$ で平均ゼロを用いる。円周距離を $d$ とし、中心を $c_I$ とすると
$$
\sup_r |Q_r*b_I(x)|\le C\frac{|I|}{d(x,c_I)^2}\|b_I\|_1.
$$
特に拡大弧の外側で積分すると $C\|b_I\|_1$ 以下である。固定 $r$ では $\sum b_I$ の $L^1$ 収束と $\|Q_r\|_\infty<\infty$ で畳み込みの和が収束する。従って全ての $r$ に対する上限も各項の上限の和で抑えられ、Tonelli で遠方での $L^1$ 評価を得る。拡大弧の和集合の測度と Chebyshev を合わせれば $Q_r$ の弱型評価が従う。

次に $r=1-\eps$ とした比較の証明 :ref[prop-four-circle-hilbert-ae] は、差の核を絶対値で積分しており $f\in L^1$ でも使える。従って $0<\eps<1/2$ で $|H_\eps f|\le|Q_{1-\eps}*f|+CM_{\mathbb T}f$ である。$\eps\ge1/2$ は核が有界である。このため切断極大関数の弱型評価も得られる。

最後に三角多項式では両方の極限が $Hf$ に一致する。一般の $f$ を多項式 $g$ で $L^1$ 近似し、各近似について尾の極限上限を
$$
\limsup |H_\eps f-Hf|
\le\sup_\eps |H_\eps(f-g)|+|H(f-g)|
$$
と抑える。$Q_r$ も同じ形で評価する。二つの弱型評価により、右辺が任意の固定閾値を超える集合の測度は $C\|f-g\|_1$ をその閾値で割った値以下である。近似誤差をゼロへ送れば両方の a.e. 収束を得る。切断積分は $\eps>0$ で連続だから可算な稠密半径を使える。
:::

:::example[一点への集中と最大絶対積分]
正の核を $F_N$、$H F_N$ を既出の式で表す。$0<t\le1$ では $\cos(t/2)\ge1/2$ と $\sin(t/2)\le t/2$ により $\cot(t/2)\ge1/t$ であり、第 2 項の絶対値は $\pi^2/(2Nt^2)$ 以下である。最後は $\sin(t/2)\ge t/\pi$ から従う。例えば $10/N\le t\le1$ に制限すると $\pi^2/(2Nt^2)\le 1/(2t)$ なので $|HF_N(t)|\ge 1/(2t)$ である。従って $\|HF_N\|_1\ge (4\pi)^{-1}\log(N/10)$（$N\ge10$）である。弱型評価は全ての $N$ に一様だが、その分布関数を全閾値について積分すると一様な $L^1$ 評価は得られない。
:::

### 実軸の弱型評価と区間の分布関数

:::proposition[実軸の Hilbert 変換の弱型 $(1,1)$]{#prop-four-real-hilbert-weak}
$f\in L^1\cap L^2(\R)$ に対し $|\{|H_\R f|>\lambda\}|\le C\lambda^{-1}\|f\|_1$ である。これにより $H_\R$ は全 $L^1$ に測度で拡張できる。
:::

:::proof
コンパクト台を持つ滑らかな $f$ について、:ref[thm-four-real-hilbert-lp] と同じ周期化 $H_L f$ を使う。円周の弱型評価を空間変数 $x$ に戻すと、両辺に付く $1/(2L)$ が消え
$$
|\{x\in[-L,L]: |H_L f(x)|>\lambda\}|\le C\lambda^{-1}\|f\|_1
$$
となる。固定 $x$ で $H_L f(x)\to H_\R f(x)$ である。$|H_\R f(x)|>\lambda$ なら十分大きい $L$ において $x$ は左辺の集合に入るので、Fatou によって実軸の評価を得る。一般の $f\in L^1\cap L^2$ は $L^1,L^2$ の両方で滑らかなコンパクト台関数に近似できる。$H_\R$ の $L^2$ 有界性から a.e. 収束部分列を取り、同じく閾値の厳密な不等号で Fatou を使えばよい。全 $L^1$ では、差への弱型評価が出力を測度で Cauchy にし、:ref[thm-four-kolmogorov] と同じ可算部分列の議論で拡張を得る。この議論は無限測度でも、差の上位集合の測度の和が有限になる部分列を使うので成立する。
:::

:::example[区間が作る分布関数の厳密値]
$f=1_{[-1,1]}$ に対し $x>0$ では $H_\R f=\pi^{-1}\log|(x+1)/(x-1)|$ は非負である。$\lambda>0$、$E=e^{\pi\lambda}$ と置く。$0<x<1$ で閾値を超える条件は $x>(E-1)/(E+1)$、$x>1$ では $x<(E+1)/(E-1)$ である。負の側は奇性で同じ長さであり、
$$
|\{|H_\R f|>\lambda\}|
=2\left(\frac{E+1}{E-1}-\frac{E-1}{E+1}\right)
=\frac4{\sinh(\pi\lambda)}
$$
を得る。$\lambda\downarrow0$ で $\lambda|\{|H_\R f|>\lambda\}|\to4/\pi$ である。これに対して大きな $\lambda$ では測度は指数的に減る。遠方の $1/x$ の尾は小さな閾値に、端点の対数的な発散は大きな閾値に対応する。
:::

:::example[区間の長さを変えた弱型係数]
$f=A1_{[b-\ell,b+\ell]}$、$A,\ell>0$ とする。核の変数を $t=\ell v$ と置くと $\ell$ は分子と分母で消えるため、$H_\R f(b+\ell x)=A H_\R1_{[-1,1]}(x)$ である。従って分布関数は $4\ell/\sinh(\pi\lambda/A)$ となる。例えば $A=2,\ell=3,\lambda=2\log2/\pi$ なら $\sinh(\pi\lambda/A)=(2-1/2)/2=3/4$ なので、閾値を超える集合の長さは $16$ である。$\sinh t\ge t$ であり、$t\downarrow0$ で $\sinh t/t\to1$ なので、この入力での $\sup_{\lambda>0}\lambda|\{|H_\R f|>\lambda\}|$ は $4\ell A/\pi$ である。元の $L^1$ ノルム $2\ell A$ との比は $2/\pi$ となる。
:::

:::proposition[実軸の最大切断と主値]{#prop-four-real-hilbert-maximal}
$1\le p<\infty$、$f\in L^p(\R)$ とし、
$$
H_\eps^\R f(x)=\frac1\pi\int_{|t|>\eps}\frac{f(x-t)}{t}\,dt,\qquad
H_*^\R f(x)=\sup_{\eps>0}|H_\eps^\R f(x)|
$$
と置く。$H_*^\R$ は弱型 $(1,1)$ および $1<p<\infty$ の強型 $(p,p)$ を持つ。$H_\eps^\R f\to H_\R f$ が a.e. に成り立つ。
:::

:::proof
各切断は Hölder で絶対収束し、半径に関して連続である。従って上限は正の有理数から取れる。まず $f$ が滑らかで $[-A,A]$ に台を持つとする。:ref[thm-four-real-hilbert-lp] と同じように $2L$ 周期へ延長する。$x$ と $\eps$ を固定して $L>|x|+A+\eps$ とすれば、円周の半径 $\pi\eps/L$ の切断は $[-A,A]$ 上で $|x-y|>\eps$ の部分を積分する。正規化した余接核はそこで $1/(\pi(x-y))$ へ収束する。従って $H_*^\R f(x)$ は周期切断の最大値の $\liminf_{L\to\infty}$ 以下である。周期の最大切断の弱型評価 :ref[thm-four-hilbert-maximal-weak] で分母の $2L$ を相殺し、$[-L,L]$ 上の上位集合に Fatou を適用すると実軸の評価になる。同じ手順で、:ref[prop-four-circle-hilbert-ae] の証明で得た周期の最大切断の強型評価も移る。そこでは $\eps\ge1/2$ の角度の切断は、核の有界性と全周平均で抑えられる。

一般の $f\in L^p$ を滑らかなコンパクト台関数 $f_j$ で近似する。固定した切断の核は $L^{p'}$ に属するので $H_\eps^\R f_j\to H_\eps^\R f$ が一様に成り立つ。半径の上限に $\liminf$ を使い、弱型と強型の評価を移す。滑らかな入力では対称切断で $f(x)$ を差し引けるため主値が存在し、:ref[thm-four-real-hilbert-l2] から $H_\R$ と一致する。一般の入力では近似 $g=f_j$ に対し $\limsup_{\eps\downarrow0}|H_\eps^\R f-H_\R f|\le H_*^\R(f-g)+|H_\R(f-g)|$ である。両項の上位集合を同じ指数で評価し、$\|f-g\|_p\to0$ とすれば a.e. 収束を得る。
:::

## Fourier 級数の Lᵖ ノルム収束

符号で分ける乗数は、変調で境界の周波数をずらせる。二つの半直線への射影の差を取れば、有限区間への射影になる。

:::theorem[$L^p$ ノルム収束]{#thm-four-series-lp}
$1<p<\infty$ の $f\in L^p(\mathbb T)$ なら $\|S_N f-f\|_p\to0$ である。
:::

:::proof
$P_{\ge0}=(I+iH+P_0)/2$ とし、$P_0f=\widehat f(0)$ とする。係数で見ると $k\ge0$ への射影である。変調 $M_j f(x)=e^{ijx}f(x)$ は等長であり、$M_jP_{\ge0}M_{-j}$ は $k\ge j$ への射影である。従って

$$
S_N=M_{-N}P_{\ge0}M_N-M_{N+1}P_{\ge0}M_{-N-1}.
$$

等式は係数の一意性で得られる。:ref[thm-four-riesz-hilbert] により右辺の作用素ノルムは $N$ によらず有界である。三角多項式に対しては十分大きな $N$ で一致する。任意の $f$ を多項式で $L^p$ 近似し、一様なノルム評価で誤差を抑える。
:::

:::example[周波数の切断をずらす]
$f=1+e^{-3ix}+2e^{4ix}$ に $M_{-2}P_{\ge0}M_2$ を作用させると、$k\ge-2$ を残して $1+2e^{4ix}$ を得る。$S_2$ を作るには、ここから $k\ge3$ への射影 $2e^{4ix}$ を引く。
:::

:::example[部分和に対する確実な誤差評価]
$f(t)=1+\sum_{k\ge1}2^{-k}\cos(kt)$ とする。各項が可和なので $S_N$ はその尾を除いた有限和になる。尾の上限は $\sum_{k>N}2^{-k}=2^{-N}$ であり、$t=0$ で等号になる。例えば $N=3$ なら一様誤差は $1/8$、各 $p$ の平均ノルムの誤差も $1/8$ 以下である。一般の $L^p$ データではこのような係数の和の評価を仮定せず、射影の有界性を使う。
:::

### 区間への射影と最良近似

同じ証明によって、対称でないどの整数区間にも一様有界な射影を作れる。また、部分和が与える近似誤差を、任意の三角多項式で達成できる誤差と比較できる。

:::proposition[任意の周波数区間への射影]{#prop-four-series-interval-projection}
$1<p<\infty$ とし、$a,b$ を整数で $a\le b$ とする。周波数 $a,\dots,b$ への射影 $S_{[a,b]}$ は
$$
S_{[a,b]}=M_aP_{\ge0}M_{-a}-M_{b+1}P_{\ge0}M_{-b-1},
\qquad \|S_{[a,b]}\|_{p\to p}\le2\|P_{\ge0}\|_{p\to p}
$$
を満たす。特に $E_{N,p}(f)=\inf_{\deg P\le N}\|f-P\|_p$ に対して
$$
E_{N,p}(f)\le\|S_Nf-f\|_p\le(1+2\|P_{\ge0}\|_{p\to p})E_{N,p}(f).
$$
ここで $\deg P\le N$ は $P$ が絶対値 $N$ 以下の周波数だけを持つことを表す。
:::

:::proof
最初の項の乗数は $1_{\{k\ge a\}}$、第 2 項は $1_{\{k\ge b+1\}}$ であり、差が求める指示関数である。各変調は等長なので上界は三角不等式から出る。$S_Nf$ 自身が次数 $N$ 以下だから左の近似不等式を得る。任意の $P$ がその次数以下なら $S_NP=P$ なので
$$
\|S_N f-f\|_p=\|(S_N-I)(f-P)\|_p\le(1+\|S_N\|_{p\to p})\|f-P\|_p.
$$
最後に右辺の下限を取る。
:::

:::example[残す次数が小さいときの最良近似]
$f(x)=e^{5ix}$ とし、$N=4$ とする。$S_N f=0$ で誤差は $1$ である。どの次数 $4$ 以下の $P$ でも $\widehat{f-P}(5)=1$ であり、Hölder と $m(\mathbb T)=1$ から $1\le\|f-P\|_1\le\|f-P\|_p$ となる。したがって $E_{4,p}(f)=1$ である。$N=5$ では $f$ 自身を取れるので誤差はゼロである。この例では最良近似の境界が周波数の切断位置に一致する。
:::

:::proposition[一次導関数からの部分和誤差]{#prop-four-series-lp-derivative-rate}
$1<p<\infty$、$f$ は絶対連続な周期関数で $f'\in L^p$ とする。$N\ge2$ に対し
$$
\|S_Nf-f\|_p\le C_p N^{-1}\|f'\|_p
$$
である。
:::

:::proof
周期を使うと $f(x-t)-f(x)=-\int_0^t f'(x-s)ds$ であり、Minkowski と平行移動の等長性で $\|f(\cdot-t)-f\|_p\le |t|\|f'\|_p$ が得られる。:ref[prop-four-jackson] で使う Jackson 核を $J_M$ とする。そこにある構成では $\int J_Mdm=1$、次数 $\le2(M-1)$、$\int|t|J_M(t)dm(t)\le C/M$ である。$M=\lfloor N/2\rfloor+1$ と選ぶと、$P=J_M*f$ は次数 $N$ 以下である。また
$$
\|f-P\|_p
\le\int J_M(t)\|f-f(\cdot-t)\|_pdm(t)
\le (C/M)\|f'\|_p.
$$
これを :ref[prop-four-series-interval-projection] の近似誤差に入れる。
:::

:::example[三角波を滑らかさで評価する]
$f(x)=|x|$（$-\pi\le x\le\pi$）を周期延長する。$f'=\sgn x$ a.e. なので $\|f'\|_p=1$ であり、上の評価は $C_p/N$ となる。さらに $n\ne0$ で部分積分を二回（$0$ で区間を分ける）行うと
$$
\widehat f(n)=\frac{(-1)^n-1}{\pi n^2}
$$
である。従って $p=2$ では Parseval の尾より $\|S_Nf-f\|_2^2=(8/\pi^2)\sum_{n>N,\ n\ {\rm odd}}n^{-4}\le8/(3\pi^2 N^3)$ となる。導関数のノルムだけを使う一般の評価より、この具体例では速い。特定の関数の係数が計算できるときは、その級数を用いると精度が上がる。
:::

### 端点でのノルム収束

:::proposition[円周の $L^1$ 端点での失敗]{#prop-four-series-l1-divergence}
$L^1(\mathbb T)$ には $\sup_N\|S_Nf\|_1=\infty$ となる $f$ が存在する。また $\|S_N\|_{L^1\to L^1}=\|D_N\|_1$ であり、$D_N$ は対称部分和の Dirichlet 核である。
:::

:::proof
Young の不等式から $\|S_N\|_{1\to1}\le\|D_N\|_1$ である。正の Fejér 核 $F_M$ の $L^1$ ノルムは $1$ であり、固定 $N$ について $D_N*F_M\to D_N$ が $L^1$ で成り立つ。実際、有限個の Fourier 係数の重み $1-|k|/M$ が全て $1$ に収束するためである。よって $\|S_NF_M\|_1\to\|D_N\|_1$ であり、逆向きの不等式も得られる。:ref[prop-four-lebesgue-constant] によればこの核ノルムは $N$ とともに無限へ行く。$L^1$ が Banach 空間なので一様有界性原理を線形作用素列 $S_N$ に適用すると、少なくとも一つの $f$ で出力ノルムの上限が無限である。
:::

:::example[有限の切断に端点ノルムの増加を見せる]
$D_N=\sum_{|k|\le N}e^{ikx}$ に $F_M$ を畳み込み、$M>N$ とする。$D_N*F_M-D_N$ の $L^1$ ノルムは各係数の絶対値の和により $2\sum_{k=1}^N k/M=N(N+1)/M$ 以下である。例えば $N=9,M=900$ では誤差は $1/10$ 以下であり、入力の $L^1$ ノルム $1$ のまま Dirichlet 核のノルムをこの精度で実現できる。ここでは入力が $N$ に依存するが、一様有界性原理は単一の入力について発散を与える。
:::

:::proposition[$L^\infty$ 端点での失敗]{#prop-four-series-uniform-divergence}
$C(\mathbb T)$ 上でも $L^\infty(\mathbb T)$ 上でも、$\|S_N\|_{\infty\to\infty}=\|D_N\|_1$ である。従って $\sup_N\|S_N f\|_\infty=\infty$ となる有界関数があり、その関数は連続に選べる。
:::

:::proof
畳み込みの絶対値を取ると $L^\infty$ 上での上界が従う。逆向きには :ref[prop-four-lebesgue-constant] の $f\mapsto S_Nf(0)$ のノルムを使う。$S_Nf$ は三角多項式なので本質上限は原点での絶対値以上であり、連続な入力に制限しても $\|D_N\|_1$ が下界になる。:ref[thm-four-continuous-divergence] の連続関数は原点の値が非有界だから、この本質上限も非有界になる。
:::

### 実軸での帯域切断

:::proposition[実軸の区間乗数の一様評価]{#prop-four-real-interval-projection}
$1<p<\infty$、$a<b$ とする。$L^2$ で乗数 $1_{[a,b]}$ を持つ作用素 $T_{a,b}$ は $L^p$ に有界に拡張し、上界は $a,b$ によらない。特に $T_{-R,R}f\to f$ が $L^p$ で成り立つ。
:::

:::proof
$P_+^\R=(I+iH_\R)/2$ は :ref[thm-four-real-hilbert-lp] により $L^p$ 有界である。変調 $M_a^\R f=e^{2\pi iax}f$ で共役すると半軸 $\xi>a$ への射影になる。$a,b$ の半軸射影の差を取れば $T_{a,b}$ であり、端点は Lebesgue 零集合なので選び方は問わない。よって一様な上界を得る。

収束にはまず $f\in\mathcal S$ を取る。Schwartz 関数 $\phi$ を $\widehat\phi=1$ が $[-1,1]$ で成り立ち、$\widehat\phi$ が $[-2,2]$ に台を持つように選べる。$f*\phi_\eps$ は周波数のコンパクト台を持ち、$\int\phi=1$ の近似単位により $L^p$ で $f$ に収束する。従って帯域制限された $L^p\cap L^2$ の関数が $L^p$ で稠密である。これらの関数は $R$ が十分大きければ変化しないので、作用素ノルムの一様性で結論が従う。
:::

:::example[有限帯域からはみ出した部分を読む]
$\widehat f=1_{[-2,2]}$ とし、$T_{-1,1}$ を考える。反転により $f(x)=\sin(4\pi x)/(\pi x)$、切断後は $\sin(2\pi x)/(\pi x)$ である。いずれも $x=0$ では連続極限を取る。$L^2$ 誤差の二乗は Plancherel で落とした周波数区間の長さ $2$ である。$p>1$ ならこれらの空間関数は遠方で $O(1/|x|)$ なので $L^p$ であり、帯域切断の一様評価を適用できる。
:::

## 上半平面の Hardy 空間（概観）

円板の正則関数を作ると負の周波数が消えた。上半平面では周波数が連続となり、正の周波数だけを持つデータが高さとともに $e^{-2\pi y\xi}$ で減衰する。

:::definition[上半平面の Hardy 空間]{#def-four-half-plane-hardy}
$1\le p<\infty$ に対し、上半平面で正則で $\sup_{y>0}\|F(\cdot+iy)\|_{L^p(\R)}<\infty$ の関数の空間を :term[Hardy 空間]{en="Hardy space" yomi="はーでぃくうかん"} $H^p(\C_+)$ とする。
:::

:::proposition[$H^2$ の境界値]{#prop-four-half-plane-hardy}
$H^2(\C_+)$ の関数は、$(-\infty,0)$ 上で $\widehat f=0$ がほとんど至る所成り立つ $f\in L^2(\R)$ の Poisson 積分に一意に表される。$\sup_y\|F(\cdot+iy)\|_2=\|f\|_2$ で、境界値へは $L^2$ と非接 a.e. の双方で近づく。
:::

:::proof
正の周波数だけを持つ $f$ に対して

$$
F(z)=\int_0^\infty \widehat f(\xi)e^{2\pi iz\xi}\,d\xi
$$

は Cauchy–Schwarz により $y>0$ で絶対収束し、同じ支配で正則である。各高さの変換は $e^{-2\pi y|\xi|}\widehat f$ であり、:ref[prop-four-upper-poisson-n] により Poisson 積分である。Plancherel と優収束でノルムおよび $L^2$ 境界値を得る。非接極限は :ref[thm-four-fatou-poisson] に従う。

逆に $F\in H^2$ とする。Cauchy–Riemann 方程式 $\partial_y F=i\partial_x F$ により、水平変換 $V_y$ は $\partial_y V_y=-2\pi\xi V_y$ を満たす。これを正当化するには、水平 Schwartz 試験関数とコンパクト台の $y$ 試験関数を取り、$x$ で切断した後に部分積分する。$L^2$ の一様上界と試験関数の急減少で切断の誤差は消える。また各水平試験関数に対する作用は $y$ に連続である。有限の $x$ 区間で正則関数の連続性を使い、外側は $L^2$ の一様上界で任意に小さくできる。従ってコンパクトな周波数区間上で $e^{2\pi y\xi}V_y$ が時刻によらない。これは $\partial_y(e^{2\pi y\xi}V_y)=0$ を試験関数に作用させ、:ref[def-four-distribution-derivative] を使って得る。負の周波数の閉区間上では $y\to\infty$ で増大するため、一様な $L^2$ 上界によりそこでは $V_y=0$ となる。正の周波数で $V_y=e^{-2\pi y\xi}v(\xi)$ と書くと、$y\downarrow0$ の単調収束から $\int_0^\infty|v|^2\le\sup_y\|V_y\|_2^2$ である。$v$ の逆変換を境界値とすれば表示を得る。
:::

:::example[極が下半平面にある正則関数]
$\widehat f(\xi)=1_{(0,\infty)}e^{-2\pi\xi}$ とすると $F(z)=1/(2\pi(1-iz))$ である。極は $z=-i$ で、$H^2$ ノルムの二乗は $\int_0^\infty e^{-4\pi\xi}d\xi=1/(4\pi)$ である。
:::

### 実軸と円周での共役核

円周の共役核は原点付近で $\cot(t/2)\sim2/t$ となる。正規化測度の $1/(2\pi)$ を含めれば、実軸の共役核に現れる定数と一致する。そこで上半平面の高さを残したまま対応を計算する。

:::proposition[上半平面の共役 Poisson 核]{#prop-four-half-plane-conjugate-poisson}
$y>0$ に対し $Q_y(x)=x/(\pi(x^2+y^2))$ と置く。その $L^2$ 変換は $-i\sgn(\xi)e^{-2\pi y|\xi|}$ である。$f\in L^2(\R)$ が実数値なら
$$
(P_y*f)(x)+i(Q_y*f)(x)=2\int_0^\infty \widehat f(\xi)e^{2\pi i(x+iy)\xi}d\xi
$$
であり、これは $H^2(\C_+)$ に属する。$Q_y*f$ の $L^2$ 境界値は乗数 $-i\sgn(\xi)$ による実軸の Hilbert 変換である。
:::

:::proof
:ref[prop-four-tempered-pv-heaviside] の正則化に現れた $x/(x^2+y^2)$ の変換を $\pi$ で割る。$Q_y\in L^2$ だが $L^1$ には属さないので、$L^2$ データとの畳み込みは Cauchy–Schwarz で点ごとに定義する。$L^2$ で近似し、変換の積は $L^1$ で収束するから逆積分と同じになる。Poisson 核の変換に $i$ 倍した共役核の変換を加えると正の半軸では $2e^{-2\pi y\xi}$、負ではゼロである。これで積分式と :ref[prop-four-half-plane-hardy] の仮定を得る。境界値の収束は Plancherel と乗数の優収束に従う。
:::

:::example[Poisson 核の共役を別の高さで読む]
$f=P_a$（$a>0$）なら $\widehat f=e^{-2\pi a|\xi|}$ である。$P_y*f=P_{y+a}$ となり、$Q_y*f=Q_{y+a}$ である。例えば $a=y=1$、$x=2$ では実部と虚部はともに $1/(4\pi)$ となる。
:::

### $1<p<\infty$ の境界値と解析的射影

$H^2$ では変換した境界値が二乗可積分な関数であった。それ以外の $p$ では Fourier 変換を通常の関数として扱えるとは限らない。そこで半軸の指示関数を掛ける操作を、Hilbert 変換から作った有界射影によって表す。

:::theorem[上半平面の $H^p$ 境界値]{#thm-four-half-plane-hardy-lp}
$1<p<\infty$ とする。$F\in H^p(\C_+)$ と同値な条件は、$f\in L^p(\R)$ が存在して $P_+^\R f=f$ かつ $F(x+iy)=P_y*f(x)$ となることである。$f$ は一意で、$F(\cdot+iy)\to f$ は $L^p$ と非接 a.e. の両方で成り立つ。
:::

:::proof
まず $F\in H^p$ ならその実部、虚部、または複素値調和関数 $F$ 自身に :ref[prop-four-harmonic-lp-boundary] を適用し、一意な $f$ と Poisson 表示を得る。:ref[thm-four-fatou-poisson] は非接収束を与える。固定高さ $t>0$ に対し $F_t=F(\cdot+it)$ と置く。

$1<p\le2$ なら $F_t\in L^2$ である。実際 $F_t$ は調和平均と Hölder により $\|F_t\|_\infty\le C t^{-1/p}\sup_y\|F_y\|_p$ を満たし、$\int |F_t|^2\le \|F_t\|_\infty^{2-p}\int|F_t|^p$ である。任意の $t>0$ について $F(z+it)$ は $H^2$ に属するので、その境界 $F_t$ は正周波数を持つ。よって $P_+^\R F_t=F_t$ が $L^p$ で成り立ち、$t\downarrow0$ の $L^p$ 極限で $P_+^\R f=f$ となる。

全ての $p$ に共通の証明として、$N>0$ に対する正則関数 $b_N(z)=(1-iz/N)^{-1}$ を使う。$|b_N|\le1$ であり、$F(\cdot+i(t+y))$ は $y\ge0$ について一様に有界だから
$$
\int |b_N(x+iy) F(x+i(t+y))|^2dx
\le \sup_{s\ge t}\|F_s\|_\infty^2
  \int \frac{dx}{(1+y/N)^2+(x/N)^2}
\le \pi N\sup_{s\ge t}\|F_s\|_\infty^2 .
$$
したがって $b_N(z)F(z+it)$ は $H^2$ に属する。その $H^2$ 境界値は連続性により $b_N(x)F_t(x)$ である。$b_N(x)\to1$ と $|b_N|\le1$ から $L^p$ で $F_t$ に戻るため、$P_+^\R$ の有界性によって $P_+^\R F_t=F_t$ が得られる。$p>2$ でもこれで上と同様に境界値が射影の像に入る。

逆に $P_+^\R f=f$ とする。$f_j\in L^p\cap L^2$ が $L^p$ で $f$ に近づくように選び、$g_j=P_+^\R f_j$ と置く。$g_j$ は $L^2$ で正周波数を持つので $P_y*g_j$ は正則である。$g_j\to f$ が $L^p$ で成り立ち、$y$ を正のコンパクト範囲に制限すれば $\|P_y\|_{p'}$ は一様有界である。従って Hölder により $P_y*g_j\to P_y*f$ が局所一様であり、極限は正則である。断面の $L^p$ ノルムは Young の不等式で $\|f\|_p$ 以下なので $H^p$ に入る。
:::

:::example[正則関数の境界がどの指数に入るか]
$F(z)=(a-iz)^{-1}$、$a>0$ とすると
$$
\|F(\cdot+iy)\|_p^p=(a+y)^{1-p}\int_\R(1+x^2)^{-p/2}dx
$$
である。$p>1$ では積分は有限であり、断面ノルムは $y\downarrow0$ で最大に近づく。$p=2$ のノルムの二乗は $\pi/a$ である。$p=1$ では積分の尾が $1/|x|$ のため発散する。そこで $F_2(z)=(a-iz)^{-2}$ に替えると $p=1$ でもノルムの上限は $\pi/a$ となる。極が上半平面の外にあるだけでは $H^1$ に属さず、水平の無限遠での減衰も必要になる。
:::

### $H^1$ の境界値と非接極大関数

$F\in H^1$ の実部と虚部に調和関数としての測度表示を個別に使うと、特異測度が残る可能性がある。正則性をさらに利用するには、$|F|$ の小さい冪を劣調和関数として扱う。そこで $L^1$ 端点の問題を、冪を取った後の $L^2$ 極大評価へ移す。

:::lemma[劣調和な冪の Poisson による上界]{#lem-four-hardy-subharmonic-majorant}
$F$ を上半平面で正則とし、各 $t>0$ で $F$ が高さ $s\ge t$ 上有界であるとする。$q>0$ なら
$$
|F(x+i(y+t))|^q\le (P_y*|F_t|^q)(x),\qquad F_t=F(\cdot+it).
$$
:::

:::proof
$\eps>0$ とし $h_\eps=(|F|^2+\eps)^{q/2}$ と置く。$s=|F|^2$ と書けば、正則性から $\Delta s=4|F'|^2$、$|\nabla s|^2=4s|F'|^2$ である。従って
$$
\Delta h_\eps= 2q(s+\eps)^{q/2-2}(\eps+(q/2)s)|F'|^2\ge0.
$$
高さを $t$ だけずらし、$h_\eps(x,y+t)$ と境界 $h_\eps(x,t)$ の Poisson 積分との差を取る。差 $w$ は有界で調和関数を引いても劣調和であり、連続なゼロ境界値を持つ。$M=\|w\|_\infty$ とする。帯 $0<y<H$ で $w$ の上界として
$$
v_\delta(x,y)=My/H+\delta(x^2+H^2-y^2)
$$
を置く。これは調和的で上下で $w$ 以上、十分遠い垂直辺でも $w$ 以上である。差 $w-v_\delta$ に劣調和関数の最大値原理を適用できる。実際、長方形で $\Delta w\ge0$ の滑らかな関数の最大値は境界上にあり、この主張は $\eta(x^2+y^2)$ を加えて厳密に $\Delta>0$ とし、内点最大での Hessian の符号から示せる。固定点で $\delta\to0$、続いて $H\to\infty$ とすれば $w\le0$ である。最後に $\eps\downarrow0$ として、有界な境界データと核の質量 $1$ による優収束を用いる。
:::

:::theorem[解析的 $H^1$ の強い境界値]{#thm-four-half-plane-hardy-l1}
$F\in H^1(\C_+)$ とし $A=\sup_{t>0}\|F_t\|_1$ とする。一意な $f\in L^1$ が存在し、$F=P_y*f$、$F_t\to f$ が $L^1$ と非接 a.e. で成り立つ。任意の開き $a>0$ について $\|N_a F\|_1\le C_a A$ である。$\widehat f$ は負の半軸上でゼロであり、逆にこの条件を満たす任意の $f\in L^1$ の Poisson 積分は $H^1(\C_+)$ に属する。この対応で $A=\|f\|_1$ である。
:::

:::proof
調和平均を使えば、固定 $t>0$ 以上の高さでは $|F|\le C t^{-1}A$ である。従って :ref[lem-four-hardy-subharmonic-majorant] を $q=1/2$ で適用できる。高さを $t$ だけずらした関数を $F^{(t)}(z)=F(z+it)$ とする。非接 Poisson 評価 :ref[thm-four-fatou-poisson] と合わせて
$$
N_a F^{(t)}(x)\le C_a \big(M(|F_t|^{1/2})(x)\big)^2.
$$
$|F_t|^{1/2}\in L^2$ なので、:ref[prop-four-maximal-strong] により $\|N_a F^{(t)}\|_1\le C_a'\|F_t\|_1\le C_a'A$ である。$t_j\downarrow0$ とすると、固定した $(z,y)\in\Gamma_a(x)$ の値 $|F(z+i(y+t_j))|$ は $|F(z+iy)|$ に収束する。このため上限を取っても $N_a F(x)\le\liminf_j N_a F^{(t_j)}(x)$ である。Fatou により $N_a F$ の $L^1$ 評価を得る。

:ref[thm-four-harmonic-l1-boundary] により $F=P_y*\mu$ の測度表示を持つ。:ref[prop-four-poisson-measure-boundary] から、非接 a.e. 極限は $\mu$ の絶対連続部分の密度 $f$ である。鉛直線で $|F_t|\le N_a F\in L^1$ なので優収束により $F_t\to f$ は $L^1$ でも成り立つ。この強収束は、$F_tdx$ の弱 * 極限が $f\,dx$ であることを意味する。測度表示の一意性から $\mu=f\,dx$ である。

最後に、固定 $t>0$ で $F(z+it)$ は一様に有界かつ断面の $L^1$ ノルムが有界だから $H^2$ に属する。従って $F_t$ の $L^2$ 変換は負の半軸でゼロである。$F_t\in L^1\cap L^2$ のため連続な $L^1$ 変換でも同じことが成り立つ。$L^1$ 収束は変換の一様収束を与えるので、極限の $\widehat f$ も負の半軸でゼロである。

逆に $f\in L^1$ の $\widehat f$ が負の半軸で消えるとする。$e^{-2\pi y|\xi|}\widehat f(\xi)$ は $y>0$ で可積分なので、$P_y*f$ を反転すれば $\int_0^\infty e^{2\pi i(x+iy)\xi}\widehat f(\xi)d\xi$ になる。この積分は正の高さのコンパクト範囲で、任意の $\xi$ の冪を掛けても可積分である。従って微分でき、Cauchy–Riemann 方程式を満たす。各断面の $L^1$ ノルムは Young の不等式から $\|f\|_1$ 以下であり、近似単位による $L^1$ 収束を合わせれば上限は $\|f\|_1$ である。
:::

:::example[二乗した Cauchy 核の境界]
$F_2(z)=(1-iz)^{-2}$ の $H^1$ ノルムは $\pi$ である。$F_2(x)=(1-x^2+2ix)/(1+x^2)^2$ なので、その実部の積分は $[x/(1+x^2)]_{-\infty}^{\infty}=0$、虚部は奇性で積分がゼロである。境界の Fourier 変換を直接得るには
$$
F_2(x)=4\pi^2\int_0^\infty \xi e^{-2\pi\xi}e^{2\pi ix\xi}d\xi
$$
と書く。積分は $\int_0^\infty \xi e^{-z\xi}d\xi=z^{-2}$（$\Re z>0$）で計算した。よって負の半軸だけでなく連続性により周波数 $0$ もゼロである。一般の $H^1$ でも $\widehat f$ の連続性と負の半軸での消失から $\int f=0$ が従う。
:::

### $H^2$ の再生核

正の周波数だけを積分する表示からは、正の高さの値を境界に対する有界線形汎関数として計算できる。その作用を内積で表すベクトルは、Cauchy 核と同じ有理式になる。

:::proposition[上半平面 $H^2$ の再生核と点評価]{#prop-four-half-plane-reproducing}
$w=x_0+iy_0\in\C_+$ とし、$H^2$ の内積は境界関数の $L^2$ 内積（第 1 変数に線形）とする。
$$
K_w(z)=\frac{i}{2\pi(z-\bar w)},\qquad
F(w)=\inner{F}{K_w}_{H^2},\qquad
\|K_w\|_{H^2}^2=\frac1{4\pi y_0}.
$$
従って $|F(w)|\le (4\pi y_0)^{-1/2}\|F\|_{H^2}$ であり、定数は最良である。
:::

:::proof
:ref[prop-four-half-plane-hardy] で境界の正周波数の変換を使う。$\xi>0$ で $k_w(\xi)=e^{-2\pi i\bar w\xi}$ と置くと
$$
\int_0^\infty \widehat f(\xi)\overline{k_w(\xi)}d\xi
=\int_0^\infty\widehat f(\xi)e^{2\pi i w\xi}d\xi=F(w).
$$
$k_w$ の二乗ノルムは $\int_0^\infty e^{-4\pi y_0\xi}d\xi=1/(4\pi y_0)$ である。対応する正則関数は $\int_0^\infty e^{2\pi i(z-\bar w)\xi}d\xi=i/(2\pi(z-\bar w))$ なので、Plancherel で内積表示を得る。Cauchy–Schwarz で点評価が従い、$F=K_w$ で等号となる。
:::

:::example[異なる高さの再生核の内積]
$w=i$、$z=3i$ では $K_i(3i)=1/(8\pi)$ である。再生公式でこれは $\inner{K_i}{K_{3i}}_{H^2}$ であり、各核のノルムの積は $1/(4\sqrt3\pi)$ である。したがって正規化した二つのベクトルの内積は $\sqrt3/2$ となる。二点が同じ高さで水平にずれている場合には、分母にその水平距離の項が現れるため、距離が増すとこの内積は小さくなる。
:::

## 演習問題

::::exercise{level=basic}
$H(3+\sin(2x)+4\cos x)$ を求めよ。
:::solution
定数は消え、$\sin(2x)$ は $-\cos(2x)$、$4\cos x$ は $4\sin x$ に移る。従って $-\cos(2x)+4\sin x$ である。
:::
::::

::::exercise{level=basic}
$f=2+\cos x$ の $\|Hf\|_2^2$ を求めよ。
:::solution
$Hf=\sin x$ である。その二乗平均は $1/2$ である。Parseval でも二つの係数の絶対値が $1/2$ なので同じ値になる。
:::
::::

::::exercise{level=basic}
$Q_r*\cos(4x)$ を求めよ。
:::solution
$Q_r=P_r H$ の乗数を持つ。まず $\cos(4x)$ が $\sin(4x)$ に移り、周波数 $\pm4$ の Poisson 係数は $r^4$ なので $r^4\sin(4x)$ である。
:::
::::

::::exercise{level=basic}
$f=e^{-2ix}+3+e^{ix}$ の $P_{\ge0}f$ を求めよ。
:::solution
周波数 $-2$ は消え、定数と正の周波数 $1$ は残る。よって $P_{\ge0}f=3+e^{ix}$ である。
:::
::::

::::exercise{level=advanced}
$H$ は三角多項式の微分と交換することを示せ。
:::solution
微分の乗数は $ik$、$H$ の乗数は $-i\sgn(k)$ である。各 $k$ で複素数の積の順序を交換できるので、二つの順序の係数が一致する。有限反転で多項式も一致する。
:::
::::

::::exercise{level=advanced}
$f$ を実三角多項式とする。$f$ の平均がゼロなら $\int f^2 Hf\,dm=\frac13\int(Hf)^3\,dm$ を示せ。
:::solution
$f+iHf$ は正の周波数だけを持つので、その三乗の定数係数もゼロである。虚部を展開すると $3 f^2Hf-(Hf)^3$ となり、その平均がゼロである。
:::
::::

::::exercise{level=advanced}
Kolmogorov の評価から、$0<q<1$ で $\int_{\mathbb T}|Hf|^qdm\le C_q\|f\|_1^q$ を示せ。
:::solution
$\int|Hf|^qdm=q\int_0^\infty \lambda^{q-1}m\{|Hf|>\lambda\}d\lambda$ と書く。$\|f\|_1=A>0$ で分け、内側では測度 $1$、外側では $CA/\lambda$ を使う。二つの積分は $A^q$ と $qCA^q/(1-q)$ である。$A=0$ なら $Hf=0$ である。
:::
::::

::::exercise{level=advanced}
$F(z)=(2-iz)^{-1}$ が $H^2(\C_+)$ に属することを示し、水平境界値の変換を求めよ。
:::solution
$\int_0^\infty e^{-4\pi\xi}e^{2\pi iz\xi}d\xi=1/(2\pi(2-iz))$ である。従って境界の変換は $2\pi e^{-4\pi\xi}1_{(0,\infty)}$ である。この関数は $L^2$ で正の半軸に台を持つので :ref[prop-four-half-plane-hardy] により属する。
:::
::::

::::exercise{level=basic}
$r=1/3$ で共役 Poisson 核 $Q_r$ の $L^1$ ノルムと、$t=\pi/2$ での値を求めよ。
:::solution
:ref[prop-four-conjugate-poisson-l1] で $(1+r)/(1-r)=2$ なのでノルムは $2\log2/\pi$ である。$t=\pi/2$ では $\sin t=1,\cos t=0$ なので値は $(2/3)/(1+1/9)=3/5$ である。
:::
::::

::::exercise{level=advanced}
$F\in H^2(\C_+)$ が $F(i)=1$ となるとき、取り得る $H^2$ ノルムの最小値と達成する関数を求めよ。
:::solution
:ref[prop-four-half-plane-reproducing] により $1\le(4\pi)^{-1/2}\|F\|_{H^2}$ なので下界は $2\sqrt\pi$ である。等号を得るには $K_i$ の定数倍を使う。$K_i(i)=1/(4\pi)$ であるから $F(z)=4\pi K_i(z)=2i/(z+i)$ とすれば条件を満たす。この関数のノルムは $4\pi(4\pi)^{-1/2}=2\sqrt\pi$ である。
:::
::::

## 確認テスト

::::quiz{answer="0" tol="1e-6"}
$H$ を定数関数 $7$ に適用した値を求めよ。
:::explanation
定数は周波数ゼロだけを持ち、そこへの乗数は $-i\sgn(0)=0$ である。従って変換はゼロである。
:::
::::

::::quiz
円周上で $\cos x$ の共役関数はどれか。

- [x] $\sin x$
- [ ] $-\sin x$
- [ ] $\cos x$

:::explanation
$\cos x=(e^{ix}+e^{-ix})/2$ の両項にそれぞれ $-i$ と $i$ を掛けると $(e^{ix}-e^{-ix})/(2i)=\sin x$ となる。
:::
::::

::::quiz{answer="0.5" tol="1e-6"}
$\|H\sin x\|_2^2$ を求めよ。
:::explanation
$H\sin x=-\cos x$ であり、正規化した二乗平均は $(2\pi)^{-1}\int_{-\pi}^{\pi}\cos^2x\,dx=1/2$ である。
:::
::::

::::quiz
上半平面の $H^2$ 境界値の変換はどこでゼロになるか。

- [x] 負の周波数
- [ ] 正の周波数
- [ ] 全周波数

:::explanation
高さの微分は変換後に $-2\pi\xi$ を掛ける。$\xi<0$ では高さとともに指数的に増えるので、一様 $L^2$ 上界によりその周波数成分はゼロになる。
:::
::::

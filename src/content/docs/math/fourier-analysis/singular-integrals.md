---
title: "Calderón–Zygmund の特異積分作用素"
description: "核の打ち消しと二乗評価から特異積分の有界性を証明し、Riesz 変換と楕円型方程式に応用する。"
notation:
  - symbol: 'K(x,y)'
    meaning: "対角線を除いた積分核"
  - symbol: 'T_\eps'
    meaning: "原点から距離 ε の内側を除いた積分"
  - symbol: 'R_j'
    meaning: "第 j Riesz 変換"
  - symbol: 'A'
    meaning: "定数係数の正定値対称行列"
---

## この章の目標

- 特異核の大きさと積分した変動を評価する。
- Calderón–Zygmund 分解から弱型評価を証明する。
- Riesz 変換を斉次核と Fourier 乗数の双方で表す。
- 二階楕円型方程式の導関数を $L^p$ で評価する。

## 特異積分核と Hörmander 条件

原点の近くで $|x|^{-n}$ 程度の核は絶対可積分でない。一方、入力のある部分の平均がゼロなら、遠い位置では核の値そのものを差し引き、核の変化だけを評価できる。

:::definition[核の Hörmander 条件]{#def-four-hormander-kernel}
対角線 $x=y$ の外で定義された可測関数 $K(x,y)$ について

$$
\sup_{y\ne z}\int_{|x-z|>2|y-z|}|K(x,y)-K(x,z)|\,dx<\infty
$$

を第 2 変数に関する :term[Hörmander 条件]{en="Hörmander condition" yomi="へるまんだーじょうけん"} と呼ぶ。随伴核は $K^*(x,y)=\overline{K(y,x)}$ である。
:::

:::proposition[微分評価から Hörmander 条件]{#prop-four-kernel-smoothness}
$K(x,y)$ が第 2 変数で微分可能で $|\nabla_y K(x,y)|\le A|x-y|^{-n-1}$ なら Hörmander 条件を満たす。畳み込み核 $K(x-y)$ に対しては $|\nabla K(w)|\le A|w|^{-n-1}$ で十分であり、その随伴核も条件を満たす。
:::

:::proof
$|x-z|>2|y-z|$ では線分 $[y,z]$ の全ての点が $x$ から距離 $|x-z|/2$ 以上離れる。従って差は $C_n A|y-z||x-z|^{-n-1}$ 以下である。極座標で積分すると半径積分は $|y-z|\int_{2|y-z|}^\infty r^{-2}dr=1/2$ になる。畳み込み核の随伴は $\overline{K(-w)}$ なので微分の大きさは同じである。
:::

:::example[一次元の原型]
$K(w)=1/(\pi w)$ に対して $|K'|=1/(\pi |w|^2)$ である。例えば $y=1,z=0$ で $|x|>2$ なら $|K(x-1)-K(x)|\le2/(\pi x^2)$ で、その積分は $2/\pi$ 以下である。原点での絶対積分の発散は、この遠方の評価に影響しない。
:::

:::example[球対称な大きさと第一成分の符号]
二次元で $K(w)=w_1/|w|^3$ とする。$\partial_j K=\delta_{1j}|w|^{-3}-3w_1w_j|w|^{-5}$ なので各成分の絶対値は $4|w|^{-3}$ 以下である。従って $|x-z|>2|y-z|$ では差の絶対値は $C|y-z||x-z|^{-3}$ 以下となる。距離 $d=|y-z|=2$ なら半径積分は $2\int_4^\infty r^{-2}dr=1/2$ となり、近い二点でも遠い二点でも同じ上界の桁になる。
:::

### 核から作用素を読むための条件

被積分関数の特異性が現れるのは対角線上である。二つの互いに離れた台を使えば、核積分には特異性が現れず、作用素がどう働くかを試験できる。ただしその試験だけでは、入力と出力の位置が同じところで作用する項を区別できない。

:::definition[標準核]{#def-four-cz-standard-kernel}
$K(x,y)$ を $x\ne y$ で定義された可測関数とする。ある $A>0$ と $0<\delta\le1$ について
$$
\begin{aligned}
|K(x,y)|&\le A|x-y|^{-n},\\
|K(x,y)-K(x,z)|&\le A|y-z|^\delta |x-z|^{-n-\delta}
\quad (|x-z|>2|y-z|)
\end{aligned}
$$
が成り立ち、第 1 変数についても同じ滑らかさの評価が成り立つとき、これを :term[標準核]{en="standard kernel" yomi="ひょうじゅんかく"} と呼ぶ。核に対応する $L^2$ 作用素とは、互いに離れたコンパクト台を持つ $f,g$ について
$$
\inner{Tf}{g}=\int\!\!\int K(x,y)f(y)\overline{g(x)}\,dy\,dx
$$
を満たすものをいう。ここでは内積の第 1 変数を線形とする。
:::

固定の入力 $f$ に対する台の外の核積分は、台から離れた各コンパクト集合上で有界である。従って上の内積の公式を、その集合内部に台を持つ全ての滑らかな $g$ に使うと、$Tf$ は台の外で核積分に a.e. に一致する。逆向きは台の外での式に $\bar g$ を掛けて積分するだけである。これは次節の核の仮定をペアリングとして確認する方法である。

:::example[同じ核を持つ二つの作用素]
$T$ が核 $K$ を持ち、$a\in L^\infty$ とする。$Tf$ と $Tf+af$ は作用素として異なり得るが、入力の台を出た点では $af$ がゼロなので、どちらも同じ核 $K$ を持つ。例えば $K=0$ なら $T=0$ と $T=3I$ の両方が対応する。従って核表示を台の外で示しただけでは、「切断極限そのものに一致する」という主張は得られない。主値を扱う場合には、滑らかな入力上の極限も独立に確認する。
:::

### Dini 型の連続度

核が微分できなくても、相対距離に対する連続度が十分速くゼロに近づけば Hörmander 条件を得られる。半径が大きくなるにつれて、二つの入力点の間隔が相対的に小さくなる点を利用する。

:::proposition[Dini 連続度からの Hörmander 条件]{#prop-four-cz-dini-kernel}
$\omega:[0,1]\to[0,\infty)$ は単調非減少で $\omega(0)=0$ とする。$|x-z|>2|y-z|$ で
$$
|K(x,y)-K(x,z)|\le A|x-z|^{-n}\omega(|y-z|/|x-z|)
$$
を満たし、$\int_0^1\omega(t)dt/t<\infty$ とする。このとき第 2 変数に関する Hörmander 定数は $A|S^{n-1}|\int_0^{1/2}\omega(t)dt/t$ 以下である。第 1 変数にも同じ条件があれば随伴核にも同じ結論が成り立つ。
:::

:::proof
$d=|y-z|>0$ とする。定義域 $|x-z|>2d$ で右辺を極座標で積分すると
$$
A|S^{n-1}|\int_{2d}^\infty \omega(d/r)\frac{dr}{r}
=A|S^{n-1}|\int_0^{1/2}\omega(t)\frac{dt}{t}
$$
である。極座標の体積 $r^{n-1}dr$ と大きさの上界 $r^{-n}$ の積が $dr/r$ を与える。第 1 変数に同じ手順を施すと、複素共役は絶対値を変えないので随伴核を扱える。
:::

:::example[対数の連続度でちょうど生じる境界]
$\omega(t)=(\log(e/t))^{-\beta}$（$0<t\le1$、$\omega(0)=0$）を考える。$u=\log(e/t)$ と置くと
$$
\int_0^{1/2}\omega(t)\frac{dt}{t}=\int_{1+\log2}^\infty u^{-\beta}du.
$$
$\beta>1$ では $(1+\log2)^{1-\beta}/(\beta-1)$、$\beta\le1$ では無限である。例えば $\beta=2$ で $1/(1+\log2)$ を得る。単に $\omega(t)\to0$ だけでは遠方の環を足し合わせられず、その収束の速さが必要である。
:::

:::example[導関数を持たない角度を含む核]
$n=2$、$0<\delta<1$ とし、
$$
K(w)=|w|^{-2}\bigl(|w_1/|w||^\delta-a_\delta\bigr),\qquad
a_\delta=\frac1{2\pi}\int_{S^1}|\omega_1|^\delta d\sigma(\omega)
$$
とする。$||s|^\delta-|t|^\delta|\le|s-t|^\delta$ である。実際、$s,t$ の絶対値を $u\ge v\ge0$ として $(v+(u-v))^\delta\le v^\delta+(u-v)^\delta$ を用いる。これは $s\mapsto s^\delta$ の増加性と導関数の減少性から従う。有向単位ベクトル $w/|w|$ の差は $|h|\le|w|/2$ で $C|h|/|w|$ 以下であり、半径の因子 $|w|^{-2}$ の差は $C|h||w|^{-3}$ 以下である。従って核の差は $C|w|^{-2}(|h|/|w|)^\delta$ 以下となる。$\delta=1/2$ の場合の Dini 積分は $\int_0^{1/2}t^{-1/2}dt=\sqrt2$ である。
:::

### 切断領域で変わる有限な項

:::proposition[球と星形領域による切断の差]{#prop-four-cz-truncation-shape}
$\Omega\in L^\infty(S^{n-1})$、$\int_{S^{n-1}}\Omega\,d\sigma=0$ とし、$K(w)=|w|^{-n}\Omega(w/|w|)$ と置く。$a:S^{n-1}\to(0,\infty)$ は可測で $0<c\le a\le C<\infty$ とする。$f\in C_c^\infty$ について、切断 $r>\eps a(\omega)$ と球切断 $r>\eps$ の差は $\eps\downarrow0$ で
$$
-f(x)\int_{S^{n-1}}\Omega(\omega)\log a(\omega)\,d\sigma(\omega)
$$
へ収束する。
:::

:::proof
差の積分は極座標で
$$
-\int_{S^{n-1}}\Omega(\omega)\int_\eps^{\eps a(\omega)}
 f(x-r\omega)\,\frac{dr}{r}\,d\sigma(\omega)
$$
である。内側は向き付き積分である。$f(x-r\omega)=f(x)+O(r)$ を入れ、$O(r)$ の定数を $\|\nabla f\|_\infty$ とする。誤差は $C'\eps\|\Omega\|_{L^1(S^{n-1})}\|\nabla f\|_\infty$ 以下でゼロに近づく。定数項の動径積分は $\log a(\omega)$ である。球での主値は、原点付近で平均ゼロを使って $f(x)$ を差し引けば絶対収束するので、星形切断にも極限が存在する。
:::

:::example[左右で除く長さを変える]
$n=1$ では球面 $S^0=\{-1,1\}$ の面積測度を計数測度とする。$\Omega(1)=1/\pi$、$\Omega(-1)=-1/\pi$ は $K(w)=1/(\pi w)$ を与える。正側で $a(1)=2$、負側で $a(-1)=1$ とすると、切断差は $-(\log2)f(x)/\pi$ となる。正負の核の係数は逆符号だが、切断で除いた長さが違うため $\log2$ が残る。この項は $f(x)$ を掛ける作用素であり、台の外での核を変えない。
:::

## Calderón–Zygmund の定理

線形作用素に対して最初に $L^2$ 有界性を仮定する。離れた点では積分核で作用が分かるとき、平均ゼロの部分に対して上の差の評価を適用できる。

:::theorem[Calderón–Zygmund の有界性]{#thm-four-cz}
$T:L^2\to L^2$ は有界な線形作用素とする。コンパクト台の $L^2$ 関数 $f$ について、台の外では $Tf(x)=\int K(x,y)f(y)dy$ が a.e. 成り立つとする。$K$ は互いに離れた二つのコンパクト集合の積上で有界とする。$K$ が第 2 変数の Hörmander 条件を満たせば $T$ は弱型 $(1,1)$ と $1<p\le2$ の強型 $(p,p)$ を持つ。$K^*$ も条件を満たせば $2<p<\infty$ でも強型である。
:::

:::proof
$f\in L^1\cap L^2$ に :ref[thm-four-cz-decomposition] を高さ $\lambda$ で適用する。$\|g\|_2^2\le2^n\lambda\|f\|_1$ だから

$$
|\{|Tg|>\lambda/2\}|\le C\lambda^{-1}\|f\|_1.
$$

選択立方体を中心を保って $4\sqrt n$ 倍に拡大した和集合を $E$ とすると、$|E|\le C_n\|f\|_1/\lambda$ である。$x\notin E$ では各中心 $c_j$ について $|x-c_j|>2|y-c_j|$ が $y\in Q_j$ に対して成り立つ。$\int b_j=0$ より

$$
Tb_j(x)=\int_{Q_j}(K(x,y)-K(x,c_j))b_j(y)dy .
$$

Hörmander 条件と Tonelli により $\int_{E^c}|Tb_j|\le C\|b_j\|_1$ である。$\sum_jb_j$ は、互いに素な台と $f,g\in L^2$ により $L^2$ で収束する。従って $T$ 後の部分和の a.e. 収束部分列に Fatou を使え、$\int_{E^c}|T\sum_j b_j|\le C\sum_j\|b_j\|_1$ である。これを Chebyshev で評価すると弱型評価が得られる。稠密な部分空間から測度での拡張は :ref[thm-four-kolmogorov] と同じ差の議論による。:ref[thm-meas-marcinkiewicz] により弱型 $(1,1)$ と強型 $(2,2)$ を補間すれば $1<p<2$ を得る。

$T^*$ は、台が互いに離れた二つのコンパクト台関数に対して積分の順序を交換すれば随伴核を持つ。実際片方を台の外の任意の試験関数に選び、局所可積分関数としての一意性を使う。従って同じ証明で $T^*$ は $1<p'<2$ に有界である。双対ノルムをコンパクト台の有界関数で取れば、$|\langle Tf,g\rangle|\le\|f\|_p\|T^*g\|_{p'}$ から $p>2$ の評価が従う。
:::

:::proposition[滑らかな畳み込み核の主値収束]{#prop-four-cz-principal-value}
さらに $T$ が畳み込み核 $K(w)$ を持ち、$|K(w)|\le A|w|^{-n}$、$|\nabla K(w)|\le A|w|^{-n-1}$、$\int_{a<|w|<b}K(w)dw=0$（$0<a<b$）を満たすとする。$T$ が滑らかなコンパクト台関数上で主値積分に一致するとき、全ての $f\in L^p$（$1<p<\infty$）で主値が a.e. に存在し $Tf$ に一致する。
:::

:::proof
$1<q<p$ を選び、$M_qf=(M(|f|^q))^{1/q}$ とする。$B=B(x,\eps/4)$、$h=f1_{B(x,2\eps)}$ と置く。$y\in B$ では $T(f-h)(y)$ は離れた台の核積分である。非コンパクト台の場合も、切断して $L^p$ と Hölder による遠方の核積分の極限を取れば表式は保たれる。$T_\eps f(x)=\int_{|x-z|>\eps}K(x-z)f(z)dz$ との差には、$\eps<|x-z|\le2\eps$ の部分と、外側の核の差がある。前者は大きさの評価で $CMf(x)$ 以下、後者は半径 $2^j\eps$ の環に分けると $C2^{-j}Mf(x)$ 以下の和になる。よって

$$
|T_\eps f(x)|\le |Tf(y)|+|Th(y)|+CMf(x)
$$

が a.e. の $y\in B$ で成立する。$y$ について平均し、$Th$ の項を $\|Th\|_q|B|^{-1/q}$ と :ref[thm-four-cz] によって評価すると $\sup_\eps|T_\eps f|\le M(Tf)+CM_qf$ になる。$M$ は $L^r$（$r>1$）に有界である。これは :ref[thm-meas-maximal-weak] の弱型 $(1,1)$ と明らかな $L^\infty$ 評価を :ref[thm-meas-marcinkiewicz] で補間すれば従う。$r=p/q>1$ に適用して $M_q$ の $L^p$ 有界性を得る。

滑らかな関数では原点近くで入力から定数を引けば、環の平均ゼロと $O(|w|)$ の差で主値が存在する。$\mathcal D$ の $L^p$ 稠密性と :ref[thm-four-maximal-convergence] の証明を、目標作用素を恒等作用素から $T$ に替えて適用する。$T$ の $L^p$ 有界性が差の第 2 項を制御する。上限は $\eps$ について連続なので有理数だけで取れる。
:::

:::example[大きな値を持つ領域の見積もり]
$\|f\|_1=6$、閾値 $\lambda=3$ なら、分解で選ぶ立方体の体積和は $2$ 以下である。$f$ の高さが非常に大きくなる場合でも、平均を閾値で制御することで良い部分の二乗ノルムは $2^n\cdot18$ 以下に保たれる。
:::

:::example[平均ゼロの対で遠方の次数が一つ下がる]
一次元で $b=1_{[0,1]}-1_{[1,2]}$ とし、$x>4$ で $K(w)=1/(\pi w)$ の積分を考える。平均ゼロなので $1/(\pi x)$ を引ける。$0\le y\le2$ では $|1/(x-y)-1/x|\le 2y/x^2$ である。このため $|Tb(x)|\le (2/\pi x^2)\int_0^2 y\,dy=4/(\pi x^2)$ であり、$x>4$ の絶対積分は $1/\pi$ 以下である。
:::

### 最大切断を変換の極大関数で抑える

元の主値収束の証明では、$L^q$ の強型評価を使うため $q>1$ を選んだ。弱型端点を扱うには小さい冪 $0<\delta<1$ を平均する。その冪については弱型 $(1,1)$ の上位集合の評価を局所積分に変換できる。

:::lemma[弱型評価から小さい冪の局所平均へ]{#lem-four-cz-kolmogorov-local}
$0<\delta<1$、$E$ は $0<|E|<\infty$ を満たし、$|\{|v|>\lambda\}|\le A/\lambda$ とする。すると
$$
\left(\frac1{|E|}\int_E |v|^\delta\right)^{1/\delta}
\le (1-\delta)^{-1/\delta}A/|E|.
$$
:::

:::proof
$A=0$ なら $v=0$ a.e. であり自明とする。分布関数を用いて $\int_E|v|^\delta=\delta\int_0^\infty\lambda^{\delta-1}|E\cap\{|v|>\lambda\}|d\lambda$ と書く。$A/|E|$ 以下では $|E|$ で、上では $A/\lambda$ で抑える。二つを積分すると合計は $A^\delta |E|^{1-\delta}/(1-\delta)$ となる。
:::

:::theorem[最大切断の Cotlar 型評価]{#thm-four-cz-cotlar}
:ref[prop-four-cz-principal-value] の核と作用素を仮定し、$T_*f=\sup_{\eps>0}|T_\eps f|$ とする。任意の $0<\delta<1$ と $f\in L^1\cap L^2$ に対し
$$
T_*f(x)\le C_\delta\big((M(|Tf|^\delta)(x))^{1/\delta}+Mf(x)\big)
$$
が a.e. で成り立つ。特に $T_*$ は弱型 $(1,1)$ を持ち、全ての $f\in L^1$ に対して $\lim_{\eps\downarrow0}T_\eps f$ は a.e. に存在して $T$ の弱型拡張と一致する。
:::

:::proof
$B=B(x,\eps/4)$、$h=f1_{B(x,2\eps)}$ とする。:ref[prop-four-cz-principal-value] と同じ核差の評価で、a.e. の $y\in B$ に対し
$$
|T_\eps f(x)|\le |Tf(y)|+|Th(y)|+CMf(x).
$$
ここで $\delta<1$ なので非負数の和を $\delta$ 乗すると各 $\delta$ 乗の和以下である。この上界を $B$ 上で平均する。$Th$ は弱型 $(1,1)$ であるため :ref[lem-four-cz-kolmogorov-local] により、その局所平均を $1/\delta$ 乗したものは $C_\delta\|h\|_1/|B|$ 以下である。この比は体積比 $8^n$ と球平均から $C_\delta Mf(x)$ 以下となる。$Tf$ の項は出発点 $x$ の球極大関数で評価できる。以上は各 $\eps$ について成り立ち、まず有理な $\eps$ で共通にし、切断積分の連続性から全ての $\eps$ に広げる。

$|Tf|^\delta=V$ は分布評価 $|\{V>s\}|\le A s^{-r}$ を満たす。ただし $r=1/\delta>1$、$A=C\|f\|_1$ である。この種の評価は $M$ を掛けても保たれる。実際 $V=V1_{\{V>s/2\}}+V1_{\{V\le s/2\}}$ と分け、後者の極大関数は $s/2$ 以下である。前者は
$$
\int_{\{V>s/2\}}V
=(s/2)|\{V>s/2\}|+\int_{s/2}^\infty|\{V>\lambda\}|d\lambda
\le A(s/2)^{1-r}\left(1+\frac1{r-1}\right)
$$
なので弱型 $(1,1)$ を使って $|\{MV>s\}|\le C_r A s^{-r}$ を得る。$s=\lambda^\delta$ と置けば $(MV)^{1/\delta}$ の弱型評価となる。$Mf$ の評価も合わせると $T_*$ の弱型評価を得る。

$L^1$ の一般の入力では、固定した $\eps$ の核は $A\eps^{-n}$ 以下なので、$L^1\cap L^2$ で近似すれば切断積分は各点で収束する。各有理 $\eps$ の上限について
$$
T_*f(x)\le\liminf_{j\to\infty}T_*f_j(x)
$$
が成り立つから、Fatou で弱型評価が移る。滑らかなコンパクト台関数では主値が存在し $Tf$ に等しいので、:ref[prop-four-maximal-oscillation] の証明を $L^1$ とこの上限に用いると一般の $f$ でも存在する。さらに $L^1$ 近似で切断誤差を $T_*(f-g)$、極限作用素の誤差を $|T(f-g)|$ で抑え、両者の弱型評価により a.e. 極限は弱型拡張と一致する。
:::

:::example[弱型の分布を局所平均へ移す]
一次元で $v(x)=1/x$、$E=(0,1)$ とする。$|\{|v|>\lambda\}|=2/\lambda$ である。$\delta=1/2$ で $\int_0^1|v(x)|^{1/2}dx=2$ なので、平均の二乗は $4$ になる。一方 $v$ の $E$ 上の絶対積分は無限である。変換の遠方での $L^1$ 端点評価と、その一部の小さい冪の局所平均が両立する理由をこの計算で確認できる。
:::

### 一つの核に複数の入力を与える

:::proposition[有限列に対するベクトル値評価]{#prop-four-cz-vector}
:ref[thm-four-cz] で $K,K^*$ が双方の変数について該当する条件を満たすとする。$1<p<\infty$ と有限列 $(f_1,\dots,f_N)$ に対し
$$
\left\|\left(\sum_{j=1}^N|Tf_j|^2\right)^{1/2}\right\|_p
\le C_p\left\|\left(\sum_{j=1}^N|f_j|^2\right)^{1/2}\right\|_p
$$
が成り立ち、定数は $N$ によらない。左辺の関数の弱型 $(1,1)$ も同じ形で成り立つ。
:::

:::proof
$\C^N$ の Euclid ノルムを $|\cdot|_2$ とし、$F=(f_j)_{j=1}^N$ と置く。$TF=(Tf_j)_{j=1}^N$ は各成分の Plancherel 型の $L^2$ 有界性を足すことで、元の $T$ と同じ $L^2(\C^N)$ 上界を持つ。二進分解では $\int_Q|F|_2$ を使って立方体を選び、$F_Q$ を成分別平均とする。積分の三角不等式 $|F_Q|_2\le |Q|^{-1}\int_Q|F|_2$ から良い部分と悪い部分のノルムはスカラー時と同じ評価になる。悪い部分は $\int_Q (F-F_Q)=0$ をベクトルとして満たす。台の外で積分に核の差を入れた後、Euclid ノルムの三角不等式と Tonelli によって
$$
\int_{E^c}|T((F-F_Q)1_Q)|_2dx
\le C\int_Q|F-F_Q|_2dy
$$
となる。良い部分に二乗評価、悪い部分にこの遠方評価を使えば、$|TF|_2$ の分布関数に弱型評価を得る。$F\mapsto|TF|_2$ は劣線形なので、Marcinkiewicz 補間で $1<p<2$ を扱える。$p>2$ では $T^*$ に同じ結果を適用し、Euclid 内積と $L^p$ 双対ノルムを組み合わせる。どの段階でも成分の個数を定数に入れていない。
:::

:::example[局所的にキャンセルする二成分]
$n=1$ とし、$F(x)=(3,4)1_{[0,1/2)}(x)$、閾値 $\lambda=2$ とする。$|F|_2=5\,1_{[0,1/2)}$ だから選択区間は $[0,1)$ であり、$F_Q=(3/2,2)$ である。悪い部分は左右で逆符号の二成分を持ち、両方の積分がゼロである。各点でその Euclid ノルムは $5/2$ である。もし各成分を別々の閾値で選ぶと、共通の区間が得られるとは限らない。ベクトルの大きさで選ぶことで、同じ拡大区間の外側を二成分に対して使える。
:::

## Riesz 変換

一次元の Hilbert 変換の乗数 $-i\sgn(\xi)$ は、$-i\xi/|\xi|$ と書ける。高次元ではその方向ベクトルの成分を使う。

:::definition[Riesz 変換]{#def-four-riesz-transform}
$L^2(\R^n)$ 上で $\widehat{R_j f}(\xi)=-i\xi_j|\xi|^{-1}\widehat f(\xi)$ と定める。原点での値は $L^2$ では任意である。この作用素を :term[Riesz 変換]{en="Riesz transform" yomi="りーすへんかん"} という。
:::

:::theorem[Riesz 変換の核と評価]{#thm-four-riesz-transform}
$R_j$ は主値核

$$
K_j(x)=c_n\frac{x_j}{|x|^{n+1}},\qquad c_n=\frac{\Gamma((n+1)/2)}{\pi^{(n+1)/2}}
$$

を持ち、弱型 $(1,1)$ および $1<p<\infty$ の強型 $(p,p)$ である。また $R_j^*=-R_j$、$\sum_{j=1}^n R_j^2=-I$ が $L^2$ 上で成立する。
:::

:::proof
$n=1$ では :ref[prop-four-tempered-pv-heaviside] に $\pi^{-1}$ を掛ければ核を得る。$n>1$ では $\alpha=n-1$ として :ref[prop-four-riesz-kernel] を使うと、$\mathcal F(|x|^{1-n})=d_n|\xi|^{-1}$、$d_n=\pi^{(n-1)/2}/\Gamma((n-1)/2)$ となる。従って $-(2\pi d_n)^{-1}\partial_j |x|^{1-n}$ の変換が $-i\xi_j/|\xi|$ である。

その微分を球を除いて部分積分すると、通常微分 $(1-n)x_j|x|^{-n-1}$ に加え、球面上で $x_j/|x|$ に比例する境界積分が現れる。球面測度と $|x|^{1-n}$ の積は $\eps$ に依存しないが、定数 $\phi(0)$ の部分は奇性でゼロである。残りは $O(\eps)$ なので微分は主値である。$(n-1)/(2\pi d_n)=c_n$ は Gamma の漸化式による。Schwartz 関数との畳み込みの変換は双対と :ref[thm-four-euclid-convolution] で確かめられる。具体的には、試験関数 $\phi$ に対し $\widehat{T*f}$ の作用は $\langle T,\widehat{\widehat f\,\phi}\rangle$ であり、$x$ の積分を先に行うとこの式になる。急減少性から試験関数のセミノルムでも交換できる。従って主値表示は乗数と一致する。

二乗有界性、随伴、平方和は各乗数を見れば従う。核の大きさと微分の評価は冪関数を微分すれば得られ、:ref[prop-four-kernel-smoothness] と :ref[thm-four-cz] を適用できる。さらに核は奇関数なので各環で積分がゼロである。Schwartz 関数との分布畳み込みは、この打ち消しで得られる各点での球主値に一致する。従って :ref[prop-four-cz-principal-value] と :ref[thm-four-cz-cotlar] により、$1\le p<\infty$ の $L^p$ 入力でも主値が a.e. に存在し、対応する拡張に一致する。
:::

:::example[二成分の方向]
二次元では $K_1(x)=x_1/(2\pi|x|^3)$、$K_2(x)=x_2/(2\pi|x|^3)$ である。周波数 $(3,4)$ における乗数はそれぞれ $-3i/5,-4i/5$ で、二乗の和は $-1$ になる。軸上では一方の乗数だけが大きさ $1$ になる。
:::

:::example[方向を入れ替えたときの乗数]
$f$ の変換が小球 $B((0,2),r)$（$r<1$）に台を持つなら、その球で $\xi_2>1$、$|\xi|<3$ なので $\|R_2 f\|_2\ge\|f\|_2/3$ である。一方 $|\xi_1|<r$ と $|\xi|>1$ から $\|R_1 f\|_2\le r\|f\|_2$ である。周波数の方向が第 2 軸に揃うほど、一次元の符号乗数に近い成分とゼロに近い成分に分かれる。
:::

### Poisson 積分に付随する共役系

上半空間の Poisson 積分に対し、各水平座標の向きに共役を作る。二次元の Cauchy–Riemann 方程式は、一般の次元では発散と回転に関する方程式になる。

:::proposition[上半空間での共役 Riesz 核]{#prop-four-conjugate-riesz-poisson}
$f\in L^2(\R^n)$、$u=P_y*f$、$v_j=P_y*(R_jf)$ とする。$v_j=Q_{j,y}*f$ であり、
$$
Q_{j,y}(x)=c_n\frac{x_j}{(|x|^2+y^2)^{(n+1)/2}}.
$$
これらの関数は調和的で、
$$
\partial_y v_j=\partial_j u,\qquad
\partial_j v_k=\partial_k v_j,\qquad
\partial_y u+\sum_j\partial_j v_j=0
$$
を満たす。$y\downarrow0$ で $v_j$ は $R_jf$ に $L^2$ と非接 a.e. で収束する。
:::

:::proof
:ref[prop-four-upper-poisson-n] の定数を $c_n$ とする。$P_s$ を $x_j$ で微分すると
$$
\partial_jP_s(x)=-(n+1)c_n sx_j(|x|^2+s^2)^{-(n+3)/2}.
$$
従って $-\int_y^\infty\partial_jP_s(x)ds=Q_{j,y}(x)$ である。有限上端 $S>y$ で積分すれば $Q_{j,y}-Q_{j,S}$ となる。$Q_{j,S}(x)=S^{-n}Q_{j,1}(x/S)$ であり、$Q_{j,1}$ は無限遠で $O(|x|^{-n})$ なので $L^2$ である。よって $\|Q_{j,S}\|_2=S^{-n/2}\|Q_{j,1}\|_2\to0$ である。

一方、有限上端では $\partial_jP_s$ の変換を積分できる。乗数は
$$
-\int_y^S 2\pi i\xi_j e^{-2\pi s|\xi|}ds
= -i\frac{\xi_j}{|\xi|}(e^{-2\pi y|\xi|}-e^{-2\pi S|\xi|})
$$
であり、$L^2$ の極限を取ると $Q_{j,y}$ の変換になる。$Q_{j,y}$ と $f$ の畳み込みは Cauchy–Schwarz で定義され、この乗数と $\widehat f$ の積は $L^1$ にある。$L^1\cap L^2$ の滑らかな関数で近似すれば、逆積分と畳み込みが一致する。従って $v_j=P_y*R_jf$ である。

水平微分に $2\pi i\xi_j$、高さの微分に $-2\pi|\xi|$ を掛けると、表示した三式がそれぞれ成り立つ。第 3 式では $\sum \xi_j^2=|\xi|^2$ を用いる。任意の正の高さで積分を微分できるのは指数の減衰による。調和性と境界の二種類の収束は Poisson 積分の性質 :ref[prop-four-upper-poisson-n] と :ref[thm-four-fatou-poisson] から従う。
:::

:::example[点状に集中した境界値の共役]
$n=2$、$f=P_1$ とする。半群則から $u=P_{y+1}$、変換の積から $v_j=Q_{j,y+1}$ である。$y=1$、$x=(3,4)$ では分母の平方根の二乗が $29$ であり、
$$
u=\frac{2}{2\pi\,29^{3/2}},\quad
v_1=\frac{3}{2\pi\,29^{3/2}},\quad
v_2=\frac{4}{2\pi\,29^{3/2}}.
$$
共役のベクトルは水平方向 $(3,4)$ に平行で、その長さと $u$ の比は $5/2$ である。原点では全ての $v_j$ がゼロになる。
:::

### 一方向の Hilbert 変換を回転平均する

$\omega\in S^{n-1}$ を固定し、各直線 $x+\R\omega$ 上の一次元 Hilbert 変換を $H_\omega$ と書く。すなわち $\omega^\perp$ で横の座標を固定し、直線の変数に $H_{\R}$ を作用させる。Fubini と :ref[thm-four-real-hilbert-lp] により $L^p$ ノルムは一次元と同じ定数で抑えられる。直交座標による Fourier 変換を取れば乗数は $-i\sgn(\xi\cdot\omega)$ である。

:::theorem[回転法と次元に依らない成分評価]{#thm-four-riesz-method-rotations}
$n\ge2$ とし、$a_n=\int_{S^{n-1}}|\omega_1|d\sigma(\omega)>0$ と置く。$L^2$ 上で
$$
R_jf=a_n^{-1}\int_{S^{n-1}}\omega_j H_\omega f\,d\sigma(\omega)
$$
である。$1<p<\infty$ で実軸の Hilbert 変換の上界 $C_p$ を用いれば $\|R_jf\|_p\le C_p\|f\|_p$ となり、定数は次元にも成分 $j$ にもよらない。
:::

:::proof
$\theta\in S^{n-1}$ に対し $\int \omega\sgn(\theta\cdot\omega)d\sigma(\omega)$ を考える。$\theta$ を第 1 軸とする座標では、第 1 成分の被積分関数は $|\omega_1|$、それ以外は座標反射で奇なので積分はゼロである。元の座標へ戻すとベクトルは $a_n\theta$ である。従って $\theta=\xi/|\xi|$ を入れれば、積分の乗数が $-i\xi_j/|\xi|$ になる。球面上の作用素積分は、まず $L^2$ の乗数を使えば可測であり、各方向の二乗上界と球面の有限測度で積分できる。Fubini と Plancherel が許されるので等式を得る。

$f$ を Schwartz とすれば $\omega\mapsto H_\omega f$ は $L^p$ 値で強連続である。実際、近い方向の変換は恒等に近い回転で共役になる。回転作用の $L^p$ 強連続性はコンパクト台の連続関数での優収束と密度から従い、$H_\omega$ の回転前の作用素は $L^p$ 有界である。従って積分に Minkowski の不等式を適用でき、$\int |\omega_j|=a_n$ と一次元の上界を使うと $C_p\|f\|_p$ を得る。最後に Schwartz 関数の稠密性で拡張する。
:::

:::example[円周の回転平均を直接計算する]
$n=2$ では $\omega=(\cos\phi,\sin\phi)$、$a_2=\int_0^{2\pi}|\cos\phi|d\phi=4$ である。$\xi=(0,2)$ とすると、$\int_0^{2\pi}\cos\phi\sgn(\sin\phi)d\phi=0$、$\int_0^{2\pi}\sin\phi\sgn(\sin\phi)d\phi=4$ となる。従って $R_1$ の乗数はゼロ、$R_2$ の乗数は $-i$ である。直接の式 $-i\xi_j/|\xi|$ からも同じ値を得る。
:::

### 勾配成分と発散ゼロ成分を分ける

:::proposition[勾配への Fourier 射影]{#prop-four-riesz-helmholtz}
ベクトル場 $F=(F_1,\dots,F_n)\in L^2(\R^n;\C^n)$ に対し $P_{\rm grad}F$ の第 $j$ 成分を $-\sum_kR_jR_k F_k$ と定める。これは各非零周波数で $\widehat F(\xi)$ を $\xi$ が張る直線に直交射影し、$F-P_{\rm grad}F$ は発散ゼロである。また $1<p<\infty$ で $L^p(\R^n;\C^n)$ に有界に拡張する。
:::

:::proof
乗数行列は $(\xi_j\xi_k/|\xi|^2)_{j,k}$ である。この行列は $\xi/|\xi|$ の直交射影で、実対称かつ冪等である。$\sum_j\xi_j((I-P_{\rm grad})\widehat F)_j=0$ なので、分布の意味でも補成分の発散はゼロになる。逆に発散ゼロなら、この内積が周波数の局所 $L^2$ 関数としてゼロであり、$P_{\rm grad}F=0$ となる。各 $R_jR_k$ の $L^p$ 有界性と有限和により $L^p$ 拡張が得られる。$F=\nabla \phi$、$\phi\in\mathcal S$ では変換は $2\pi i\xi\widehat\phi$ だから射影で保たれる。
:::

:::example[斜めの周波数で場を二成分に分ける]
$\xi=(3,4)$ で振幅ベクトル $(1,0)$ を考える。勾配成分は $(3,4)(3/25)=(9/25,12/25)$、発散ゼロ成分は $(16/25,-12/25)$ である。後者と $(3,4)$ の内積は $(48-48)/25=0$ であり、二つの成分同士の内積もゼロである。この計算は単一の周波数における行列の作用であり、$L^2$ 全体では各周波数の二乗を積分する。
:::

## 斉次核をもつ特異積分

角度ごとの密度を指定すると、半径方向の積分は $dr/r$ である。この積分は無限遠と原点に問題を持つので、球面上の平均ゼロを用いて半径の基準を差し引く。

:::theorem[平均ゼロの滑らかな斉次核]{#thm-four-homogeneous-cz}
$n\ge2$、$\Omega\in C^1(S^{n-1})$、$\int_{S^{n-1}}\Omega\,d\sigma=0$ とし、$K(x)=|x|^{-n}\Omega(x/|x|)$ と置く。環 $\eps<|x|<R$ で切断した畳み込みは、$\eps\downarrow0,R\uparrow\infty$ で $L^2$ の強極限 $T$ を持つ。$T$ は弱型 $(1,1)$ および $1<p<\infty$ の強型 $(p,p)$ である。
:::

:::proof
$\xi\ne0$ について $\theta=\xi/|\xi|$、$s=2\pi|\xi|r$ に置き換える。乗数の角度積分内の動径積分は、$\Omega$ の平均ゼロにより

$$
\int_{2\pi|\xi|\eps}^{2\pi|\xi|R}
\frac{e^{-is\omega\cdot\theta}-1_{\{s\le1\}}}{s}\,ds
$$

に置き換えられる。$u=\omega\cdot\theta\ne0$ とする。$s\le1$ は $|e^{-isu}-1|\le s|u|$ で一様に有界で、$1<s\le1/|u|$ は $\log(1/|u|)$ 以下である。$s>1/|u|$ の振動積分は部分積分で一様に有界であり、上限を無限大に送った極限も存在する。このため全体は $C(1+|\log|u||)$ 以下である。球面で $u$ の密度は $(1-u^2)^{(n-3)/2}$ に比例し、$u=0$ 付近では対数を積分できる。$\Omega$ は有界なので優収束で乗数が $\xi\ne0$ で収束し、一様有界となる。Plancherel により $L^2$ 強極限を得る。

滑らかなコンパクト台関数で原点近くの定数を引けば、$\int_{\eps<|w|<b}K(w)dw=0$ によって原点の発散が消える。台から離れた点では通常積分である。このため極限の核表示が従う。大きさと微分の評価は斉次性と球面での $C^1$ 条件から従い、:ref[thm-four-cz] を用いる。
:::

:::example[角度の二次成分]
二次元で $\Omega(\omega)=\omega_1^2-\omega_2^2$ は座標交換で符号が反転するので平均ゼロである。対応する核は $(x_1^2-x_2^2)/|x|^4$ で、偶関数である。Riesz 変換の核は奇関数だが、斉次核の評価に必要なのは球面平均がゼロになることである。
:::

:::example[三次の角度成分]
$n=2$ で $\Omega(\cos\theta,\sin\theta)=\cos(3\theta)$ を取る。$0\le\theta\le2\pi$ の積分は $\sin(3\theta)/3$ の端点差でゼロになる。このため半径 $1$ から $3$ の環でも核の積分はゼロである。一方、下半円と上半円での値を比べると $\Omega(\theta+\pi)=-\Omega(\theta)$ であり、この例は奇核の条件も満たす。
:::

### 動径積分が与える乗数の角度依存性

環で切断した乗数の収束では、振動積分を対数で支配した。極限の動径積分を実際に評価すると、偶核と奇核で異なる成分が残る。

:::proposition[球面核の対数と符号による乗数]{#prop-four-homogeneous-symbol}
:ref[thm-four-homogeneous-cz] の仮定の下、$T$ の乗数は $\xi\ne0$ で
$$
m(\xi)=\int_{S^{n-1}}\Omega(\omega)
\left(-\log|\omega\cdot\theta|-\frac{\pi i}{2}\sgn(\omega\cdot\theta)\right)d\sigma(\omega),
\qquad \theta=\xi/|\xi|.
$$
偶関数の $\Omega$ では符号項、奇関数の $\Omega$ では対数項が消える。
:::

:::proof
動径変数を $s=2\pi|\xi|r$ に変え、$u=\omega\cdot\theta$ とする。平均ゼロを使って $e^{-s}$ を差し引き、$\eta>0$ で $e^{-\eta s}$ を掛けて
$$
\int_0^\infty \frac{e^{-(\eta+iu)s}-e^{-(\eta+1)s}}s ds
=\log(\eta+1)-\log(\eta+iu)
$$
と求める。等式の証明には、右半平面にある二つの複素数を線分で結び、そのパラメータで被積分関数を微分する。$s$ が約分されて $\int_0^\infty e^{-zs}ds=1/z$ となるから、線分積分が対数の差を与える。右半平面では主値の偏角を選ぶ。

$u\ne0$ で $\log(\eta+iu)\to\log|u|+i(\pi/2)\sgn(u)$ である。$0<\eta\le1$ では対数の絶対値が $C+|\log|u||$ 以下であり、球面上で可積分なので優収束が適用できる。動径積分による元の乗数と同じであることは、各 $u\ne0$ について $e^{-\eta s}$ を外す手順から分かる。$0<s<1$ では分子が $O(s)$、有限範囲 $1<s<A$ では連続性で極限を取れる。$s>A$ の振動項は部分積分して $2/(A|u|)$ 以下で、これは $\eta$ によらない。減衰項の尾も一様に消える。このため減衰を外した積分に収束し、:ref[thm-four-homogeneous-cz] の切断極限と同じである。最後に球面で $\omega\mapsto-\omega$ を使うと、偶と奇の組合せの積分がゼロになる。
:::

:::corollary[可積分な奇角度に対する強型評価]{#cor-four-odd-rough-rotation}
$\Omega\in L^1(S^{n-1})$ が奇関数なら、$K=|x|^{-n}\Omega(x/|x|)$ の球切断は滑らかなコンパクト台関数で主値を持ち、その作用素は $1<p<\infty$ で
$$
\|Tf\|_p\le \tfrac\pi2 C_p\|\Omega\|_{L^1(S^{n-1})}\|f\|_p
$$
と評価できる。
:::

:::proof
有限環での積分を極座標に直す。$\Omega(-\omega)=-\Omega(\omega)$ によって
$$
\begin{aligned}
\int_{\eps<|w|<R}K(w)f(x-w)dw
&=\frac12\int_{S^{n-1}}\Omega(\omega)\int_\eps^R
   \frac{f(x-r\omega)-f(x+r\omega)}r\,dr\,d\sigma(\omega)\\
&=\frac\pi2\int_{S^{n-1}}\Omega(\omega)H_{\omega,\eps,R}f(x)d\sigma(\omega),
\end{aligned}
$$
ここで $H_{\omega,\eps,R}$ は直線に沿う一次元の有限切断である。滑らかなコンパクト台関数と固定 $x$ について、$r<1$ の差は $2r\|\nabla f\|_\infty$ 以下である。また $r>|x|+A$（台を $B(0,A)$ に含めたとき）は両項がゼロである。従って $\Omega$ に関する優収束を使って主値極限を取れ、$\frac{\pi}{2}\int\Omega(\omega)H_\omega f(x)d\sigma$ が得られる。Minkowski と直線に沿う $L^p$ 有界性により結論を得る。拡張は滑らかなコンパクト台関数の稠密性で決まる。この証明は角度の連続性を用いず、奇性により一次元の変換だけに帰着した。
:::

:::example[半円だけで符号を切り替える核]
$n=2$ で $\Omega(\omega)=\sgn(\omega_1)$ とすると、不連続な奇核であり $\|\Omega\|_1=2\pi$ である。上の評価の定数は $\pi^2C_p$ となる。$\xi=(1,0)$ では $\Omega(\omega)\sgn(\omega\cdot\xi)=1$ a.e. なので乗数は $-i\pi^2$ である。$\xi=(0,1)$ では積が第 1 座標の反射で奇になりゼロである。回転平均の乗数を使えば、角度に微分がない場合にも二つの向きでの差が計算できる。
:::

### 球面調和成分の Fourier 変換

球面上で多項式の成分を一つ選ぶと、角度を表す関数は変換後にも同じ形で現れる。次数だけで変わる定数と位相を求めるには、まず Gauss 関数に多項式を掛けてから半径方向を積分する。

:::proposition[斉次調和多項式の核と位相]{#prop-four-harmonic-polynomial-kernel}
$P$ は次数 $k\ge1$ の斉次調和多項式とし、$K_P=\operatorname{pv}(P(x)/|x|^{n+k})$ を球切断で定める。その Fourier 変換は非零周波数で
$$
\widehat{K_P}(\xi)=
(-i)^k\pi^{n/2}\frac{\Gamma(k/2)}{\Gamma((n+k)/2)}
     \frac{P(\xi)}{|\xi|^k}
$$
であり、右辺の有界関数が表す緩増加超関数に一致する。
:::

:::proof
$P$ の球面平均は調和関数の平均値公式により $P(0)=0$ である。従って原点の近くで試験関数からその原点値を引けば、主値が定義できる。$0<\alpha<n$ を固定して $K_{P,\alpha}(x)=P(x)|x|^{-n-k+\alpha}$ とし、$a=(n+k-\alpha)/2$ と置く。Gamma 積分で
$$
K_{P,\alpha}(x)=\frac{\pi^a}{\Gamma(a)}
           \int_0^\infty t^{a-1}P(x)e^{-\pi t|x|^2}dt
$$
と書ける。この積分を Schwartz 試験関数に作用させた絶対積分は、$\int |P(x)||x|^{-2a}|\phi(x)|dx$ の定数倍以下で有限である。

$P(\partial_\xi)e^{-\pi|\xi|^2/t}=(-2\pi/t)^k P(\xi)e^{-\pi|\xi|^2/t}$ である。これには訳された Gauss 関数を $h$ について微分する。すなわち
$$
e^{-\pi|\xi+h|^2/t}
=e^{-\pi|\xi|^2/t}e^{-2\pi\xi\cdot h/t}e^{-\pi|h|^2/t}.
$$
次数 $k$ の項のうち $|h|^{2j}$ を含むものに $P(\partial_h)|_{h=0}$ を作用させると、$j\ge1$ では $\Delta^j P$ が現れてゼロになる。実際、多項式の係数ペアリング $P(\partial_h)Q(h)|_{h=0}$ では $h_j$ を $Q$ に掛けることは $P$ の $\partial_j$ に対応し、これは単項式で階乗を照合すれば得られる。残る指数の項から所望の公式を得る。Fourier 変換での $x_j$ の作用を $(i/(2\pi))\partial_{\xi_j}$ とすれば
$$
\mathcal F(P(x)e^{-\pi t|x|^2})(\xi)
=(-i)^k t^{-n/2-k}P(\xi)e^{-\pi|\xi|^2/t}.
$$

$t$ 積分を代入して $b=(k+\alpha)/2$ と置くと、$t$ の冪は $t^{-b-1}$ である。$s=\pi|\xi|^2/t$ に置き換えれば
$$
\widehat{K_{P,\alpha}}(\xi)=
(-i)^k\pi^{n/2-\alpha}
  \frac{\Gamma((k+\alpha)/2)}{\Gamma((n+k-\alpha)/2)}
  \frac{P(\xi)}{|\xi|^{k+\alpha}}.
$$
交換の絶対収束は、右辺の絶対値の計算に現れる $|\xi|^{-\alpha}|\phi(\xi)|$ が可積分であることから従う。最後に $\alpha\downarrow0$ とする。空間側では球面平均ゼロにより $|x|<1$ で $\phi(x)-\phi(0)$ を使い、その絶対値を $C|x|$ で抑えれば優収束で主値に近づく。遠方は $\alpha$ を固定の $\alpha_0<n$ 以下にして Schwartz 減衰を使う。周波数側でも原点近くを $|\xi|^{-\alpha_0}$、遠方を有界関数で支配できる。こうして緩増加超関数として式が得られる。
:::

:::example[二次の核と Beurling–Ahlfors 変換]
$n=2$、$P(x)=(x_1-ix_2)^2$ とする。$\partial_1^2P+\partial_2^2P=0$ なので調和であり、上の定数は $-\pi$ である。従って複素座標 $z=x_1+ix_2$ で核 $-1/(\pi z^2)=-(x_1-ix_2)^2/(\pi |x|^4)$ は、$\xi=(\rho\cos\phi,\rho\sin\phi)$ における乗数 $e^{-2i\phi}$ を持つ。この作用素を :term[Beurling–Ahlfors 変換]{en="Beurling–Ahlfors transform" yomi="ばーりんぐあーるふぉるすへんかん"} という。例えば $\xi=(1,1)$ では $\phi=\pi/4$ で乗数は $-i$ である。乗数の絶対値は $1$ なので Plancherel により $L^2$ で等長となる。$\Omega$ は滑らかで平均ゼロなので :ref[thm-four-homogeneous-cz] の $L^p$ 評価も適用できる。
:::

### 角度に滑らかさを仮定しない二乗評価

乗数の式で絶対値を取りにくいのは、$\omega\cdot\theta=0$ 付近の対数項である。角度密度の大きな値と、赤道に近い点の寄与を一緒に評価すると、微分を使わずに二乗評価を得られる。

:::proposition[角度密度の $L\log L$ 条件]{#prop-four-rough-angular-l2}
$n\ge2$ とし、球面上の $\Omega$ は平均ゼロで $\int |\Omega|\log(e+|\Omega|)\,d\sigma<\infty$ とする。核 $|x|^{-n}\Omega(x/|x|)$ の有限環での切断は、原点側と無限遠側の切断を外すとき $L^2$ 強極限を持つ。極限の乗数は :ref[prop-four-homogeneous-symbol] の角度積分で与えられる。
:::

:::proof
$u=\omega\cdot\theta$ と置く。:ref[thm-four-homogeneous-cz] の動径積分の評価は、切断の端点に一様に $C(1+\log(1/|u|))$ である。$A\ge0$、$B\ge0$ に対し $AB\le A\log(e+A)+e^B$ である。実際、$A>0$ のとき $AB-A\log A+A\le e^B$ は、左辺を $A$ の関数として最大化して得られる。この不等式を $A=|\Omega|$、$B=\frac12\log(1/|u|)$ に使うと、対数を掛けた積は
$$
|\Omega|\log(1/|u|)\le 2|\Omega|\log(e+|\Omega|)+2|u|^{-1/2}
$$
で抑えられる。後者の積分は回転で $\theta$ に依存せず、赤道の近くでは $\int_{-1/2}^{1/2}|u|^{-1/2}du<\infty$、赤道から離れれば球面の有限測度で有限である。従って角度積分の優収束が使え、その上限は $\theta$ に一様である。動径積分の極限は :ref[prop-four-homogeneous-symbol] の計算の通りであり、差し引いた基準項の違いは $\Omega$ の平均ゼロで消える。最後に Plancherel と優収束で $L^2$ 強極限を得る。
:::

:::example[赤道に向かって発散する角度密度]
$n=2$ とし、$0<\beta<1$ に対して $\Omega(\omega)=|\omega_1|^{-\beta}-c_\beta$ とする。$c_\beta$ は球面平均である。$\omega_1=u$ がゼロに近い場所の積分は $\int_0^{1/2}u^{-\beta}(1+|\log u|)du$ で制御できる。$u=e^{-v}$ と置くと $v\ge\log2$ での $(1+v)e^{-(1-\beta)v}$ の積分になり、有限である。従ってこの非有界な偶核でも二乗極限を得る。角度に微分の仮定を置いた弱型評価を適用するには、別途その仮定を確認する必要がある。
:::

:::theorem[角度密度の $L\log L$ 条件による強型評価]{#thm-four-rough-angular-lp}
:ref[prop-four-rough-angular-l2] の仮定の下で、切断の極限 $T$ は $1<p<\infty$ の $L^p$ 作用素に延長される。また各有限環での切断は、この拡張に $L^p$ で強収束する。
:::

:::proof
まず $\Omega\in L^2(S^{n-1})$ として $A=\|\Omega\|_1$、$B=\|\Omega\|_2$ と置く。$K_j=K1_{\{2^j<|x|\le2^{j+1}\}}$ とする。$d=2^j|\xi|$ が $1$ 以下なら球面平均ゼロで指数から $1$ を引き、$|\widehat K_j(\xi)|\le CAd$ を得る。$d>1$ では $\theta=\xi/|\xi|$ とし、$|\omega\cdot\theta|\le d^{-1/2}$ の帯を分ける。その測度は $C_n d^{-1/2}$ 以下なので、Cauchy–Schwarz で帯の寄与は $CBd^{-1/4}$ 以下である。帯の外では動径積分を部分積分すると $C/(d|\omega\cdot\theta|)\le C d^{-1/2}$ であり、この寄与も $CBd^{-1/4}$ 以下になる。従って
$$
|\widehat K_j(\xi)|\le C B\min(d,d^{-1/4}).
$$

次に $\Delta_k$ を、乗数 $\varphi(2^{-k}\xi)$ を掛ける滑らかな切断とする。$\varphi=\eta-\eta(2\cdot)$、$\eta$ は半径 $1$ の球で $1$、半径 $2$ の球外で $0$ となる滑らかな球対称関数とする。すると $\sum_k\varphi(2^{-k}\xi)=1$（$\xi\ne0$）で、各周波数で重なりは有界である。$\widetilde\Delta_k$ は $\varphi$ の台の近傍で $1$ となる実数値の滑らかな環の切断を同じように拡大したものとする。いずれの切断族 $D_k$ にも
$$
\left\|(\textstyle\sum_k |D_k f|^2)^{1/2}\right\|_r\le C_r\|f\|_r
\qquad (1<r<\infty)
$$
が成り立つ。実際、有限列を取って $f\mapsto(D_kf)_k$ を考える。$L^2$ 上では台の重なりと Plancherel で二乗和を抑えられる。核の各成分は $2^{kn}\check\varphi(2^k x)$ の形である。$2^k|x|\le1$ と $2^k|x|>1$ に分けて Schwartz 減衰の等比級数を足すと、成分の絶対和は $C|x|^{-n}$、勾配の絶対和は $C|x|^{-n-1}$ 以下となる。従って Euclid ノルムでも同じ評価であり、随伴の核にも同じ上界がある。:ref[prop-four-cz-vector] の証明では核の絶対値を作用素ノルムに替えれば、スカラー入力から有限ベクトルへの写像とその随伴にも、同じ二進分解が使える。ここでは個数に依らない核の差の積分上界と $L^2$ 評価があるので、その弱型評価と補間・双対化から上の不等式が従う。

有限な添字集合 $F$ と整数 $l$ について $S_l f=\sum_{j\in F}K_j*\Delta_{l-j}f$ とする。この各項で $d$ は $2^l$ と同程度なので、Plancherel によって $\|S_l\|_{2\to2}\le C B2^{-|l|/4}$ となる。これとは別に、全ての $1<r<\infty$ で $\|S_l\|_{r\to r}\le C_r A$ を示す。$U_jh=|K_j|*h$ とすると、非負の $h$ に対し
$$
\sup_j U_j h(x)\le C\int_{S^{n-1}}|\Omega(\omega)| M_\omega h(x)d\sigma(\omega).
$$
ここで $M_\omega$ は方向 $\omega$ の直線上の中心付き一次元極大関数である。動径積分で $dr/r\le2^{-j}dr$ を使い、積分を $[-2^{j+1},2^{j+1}]$ に広げれば従う。従って一次元極大関数の強型評価と Minkowski で、左辺の $L^s$ ノルムは $C_s A\|h\|_s$ 以下である。反射核に替えてもこの評価は保たれる。

$r>2$ のとき、Cauchy–Schwarz により $|K_j*g_j|^2\le \|K_j\|_1U_j(|g_j|^2)$ である。両辺を足して非負の $h\in L^{(r/2)'}$ と積分すると、反射核の上限の評価により
$$
\int h\sum_j|K_j*g_j|^2
\le C_r A^2\|h\|_{(r/2)'}\left\|\sum_j |g_j|^2\right\|_{r/2}
$$
となる。$r=2$ は Young の不等式を各成分で用いればよい。$r<2$ は対角作用素 $(g_j)\mapsto(K_j*g_j)$ の随伴に $r'>2$ の評価を使う。このベクトル評価に $g_j=\Delta_{l-j}f$ を入れる。双対試験では $\langle K_j*\Delta_{l-j} f,h\rangle=\langle K_j*\Delta_{l-j}f,\widetilde\Delta_{l-j}h\rangle$ なので、二つの二乗和に Cauchy–Schwarz と Hölder を使えば $\|S_l\|_{r\to r}\le C_r A$ を得る。

$1<p<2$ なら $1<r<p$、$p>2$ なら $r>p$ を一つ選ぶ。$2$ と $r$ の間で補間すると $\|S_l\|_{p\to p}\le C_p A^{1-\theta}B^\theta 2^{-\theta|l|/4}$ となる。ただし $0<\theta<1$ は二乗評価の重みである。$p=2$ では $\theta=1$ でよい。一方、$p$ での直接の上界は $C_p A$ である。この二つの小さい方を $l$ について足し、$|l|\le4\log_2^+(B/A)+1$ とその外側に分ければ、全体は $C_p A(1+\log^+(B/A))$ 以下になる。ここで $\log^+ v=\max(0,\log v)$ とする。$A=0$ なら作用素はゼロである。$L^2$ 上で $\sum_l S_l f=\sum_{j\in F}K_j*f$ は切断の復元で従うため、有限環の和はこの $L^p$ 上界を持つ。

一般の $\Omega$ を $|\Omega|\le2$ と $2^m<|\Omega|\le2^{m+1}$（$m\ge1$）に分け、各部分からその球面平均を引く。引いた平均を足した値は元の平均ゼロである。元の各部分の $L^1$ ノルムを $a_m$ とし、$|\Omega|\le2$ の部分は $m=0$ とする。平均を引いた部分の $L^1$ ノルムは $2a_m$ 以下、$L^2$ ノルムは $C2^{m/2}\sqrt{a_m}$ 以下である。前段の和の評価にこれらの上界を直接入れると、各部分の上界は $C_p a_m(1+m+\log^+(1/a_m))$ 以下である。$a_m=0$ の項はゼロと読む。$\sum_m(1+m)a_m<\infty$ は仮定による。さらに $a_m\ge e^{-(m+1)}$ なら $a_m\log^+(1/a_m)\le(m+1)a_m$ である。$0<a_m<e^{-(m+1)}$ なら $a\log(1/a)$ はこの範囲で増加するので $a_m\log(1/a_m)\le(m+1)e^{-(m+1)}$ である。従って対数の項も可和である。有限環では角度密度の和を $L^1$ で取って核積分に戻せる。これで一般の場合も二進環の有限和について一様な $L^p$ 評価が得られる。

任意の有限環は二進環の有限和と両端の高々二つの部分環に分かれる。部分環の核の $L^1$ ノルムは $A\log2$ 以下なので、Young の不等式で同じ一様評価が成り立つ。$L^2$ 強収束は :ref[prop-four-rough-angular-l2] による。$p\ne2$ に対しては、Schwartz 入力上で $L^2$ の Cauchy 評価と、$p$ より同じ側にある $L^r$ の一様評価を Hölder で補間すると $L^p$ でも Cauchy となる。稠密性と切断の一様評価で $L^p$ 全体に収束が移り、共通部分では $L^2$ の極限に一致する。
:::

## 楕円型方程式の Lᵖ 評価への応用

$-\Delta u=f$ の変換は $4\pi^2|\xi|^2\widehat u=\widehat f$ である。$u$ 自体を復元する際は原点の特異性があるが、二階微分は同じ次数の分子を持つため、比が有界な乗数になる。

:::theorem[定数係数楕円作用素の二階評価]{#thm-four-elliptic-lp}
実対称正定値行列 $A$ と $L_A=-\sum_{j,k}a_{jk}\partial_j\partial_k$ に対し、$1<p<\infty$ なら

$$
\sum_{j,k}\|\partial_j\partial_k u\|_p\le C_{n,p,A}\|L_Au\|_p\qquad(u\in C_c^\infty(\R^n)).
$$
:::

:::proof
$A=I$ の場合、$f=-\Delta u$ と置くと $\widehat{\partial_j\partial_k u}=-\xi_j\xi_k|\xi|^{-2}\widehat f$、従って $\partial_j\partial_k u=R_jR_k f$ である。:ref[thm-four-riesz-transform] の有界性を二度適用する。一般の $A$ では対称正定値平方根 $B=A^{1/2}$ を取り、$v(y)=u(By)$ と置く。連鎖律で $D^2v(y)=B^{\mathsf T}D^2u(By)B$、従って $-\Delta v=(L_Au)(By)$ である。$D^2u$ は $D^2v$ に $B^{-1}$ を両側から掛けて戻り、両ノルムに共通な Jacobian $|\det B|^{-1/p}$ は相殺される。定数に $B^{-1}$ の成分を含めれば結論となる。
:::

:::example[異なる二軸の尺度]
$L_A=-4\partial_1^2-\partial_2^2$ では $v(y)=u(2y_1,y_2)$ とする。$\partial_{y_1}^2v=4\partial_{x_1}^2u(By)$、$\partial_{y_1}\partial_{y_2}v=2\partial_{x_1}\partial_{x_2}u(By)$ である。等方的な評価から各成分を戻すと、二つの微分の前置係数が異なる理由が分かる。
:::

### 局所化した楕円評価

全空間での評価を局所で使うには、調べたい領域の近くでだけ $u$ を残す。切断関数で微分をすると余分な一次項と零次項が出るため、それらを分けて評価する。

:::corollary[二つの球の間での二階評価]{#cor-four-elliptic-local}
$r>0$、$1<p<\infty$ とし、$u$ が $B(x_0,2r)$ の近傍で滑らかで $L_Au=f$ を満たすとする。すると
$$
\|D^2u\|_{L^p(B(x_0,r))}
\le C\left(\|f\|_{L^p(B(x_0,2r))}
+r^{-1}\|\nabla u\|_{L^p(B(x_0,2r))}
+r^{-2}\|u\|_{L^p(B(x_0,2r))}\right)
$$
である。$D^2$ のノルムは成分ノルムの和としてよい。
:::

:::proof
$\chi$ を $B(x_0,r)$ の近くで $1$、$B(x_0,2r)$ に台を持つ滑らかな切断とし、$\|\partial^\alpha\chi\|_\infty\le C_\alpha r^{-|\alpha|}$ となるよう固定切断を拡大する。積の微分により
$$
L_A(\chi u)=\chi f-2\sum_{j,k}a_{jk}(\partial_j\chi)(\partial_k u)+(L_A\chi)u .
$$
右辺を各ノルムで抑え、:ref[thm-four-elliptic-lp] を $\chi u$ に適用する。左辺の二階微分は内側の球で $u$ の二階微分と一致する。
:::

:::example[半径で生じる係数]
$A=I$、$r=1/2$ なら一次項の前には $2$、零次項の前には $4$ が掛かる。例えば $u(x)=x_1^2-x_2^2$ は調和関数なので $f=0$ だが、二階微分はゼロでない。局所化の評価で右辺に $u$ または一次導関数の情報を含める理由がこの例で確認できる。
:::

### Newton ポテンシャルの二階微分と局所項

$-\Delta$ の基本解を二度微分すると、原点以外では $|x|^{-n}$ の核になる。ただし二度目の分布微分では球面上の境界項が残る。この項を含めると、対角成分の和が元の Poisson 方程式に戻る。

:::proposition[二階 Riesz 変換の主値核]{#prop-four-second-riesz-kernel}
$n\ge2$、$s_{n-1}=|S^{n-1}|$ とし、
$$
K_{jk}(x)=\frac{n x_jx_k-\delta_{jk}|x|^2}{s_{n-1}|x|^{n+2}}.
$$
このとき $R_jR_k f=\operatorname{pv}(K_{jk})*f-(\delta_{jk}/n)f$ が Schwartz 関数上で成り立つ。従って $R_jR_k$ は弱型 $(1,1)$ を持つ。$L^1$ 入力で核の主値は a.e. に収束し、この局所項を含めた式が $R_jR_k$ の弱型拡張を与える。
:::

:::proof
$n>2$ で $\Phi(x)=|x|^{2-n}/((n-2)s_{n-1})$、$n=2$ で $\Phi(x)=-\log|x|/(2\pi)$ と置く。両方の $\partial_j\Phi=-x_j/(s_{n-1}|x|^n)$ は通常関数として局所可積分である。実際、分布微分を求めるために $|x|>\eps$ で部分積分すると、最初の境界項は $n>2$ で $O(\eps)$、$n=2$ で $O(\eps|\log\eps|)$ なので消える。これを分布としてさらに微分し、$|x|>\eps$ で部分積分すると、原点の外側にある微分 $K_{jk}$ と境界項
$$
-\frac1{s_{n-1}}\int_{S^{n-1}}\omega_j\omega_k\phi(\eps\omega)d\sigma(\omega)
$$
が出る。ここで球を除いた領域の外向き法線は $-\omega$ である。$\int_{S^{n-1}}\omega_j\omega_k\,d\sigma/s_{n-1}=\delta_{jk}/n$ なので、$\eps\downarrow0$ で境界項は $-(\delta_{jk}/n)\phi(0)$ となる。球面平均の式は $j\ne k$ では反射、$j=k$ では対称性と $\sum_j\omega_j^2=1$ から従う。対角成分を足すと $\sum_jK_{jj}=0$ であり、局所項の和は $-\delta_0$ である。これが $-\Delta\Phi=\delta_0$ を与える。基本解の境界値問題への応用は [Laplace 方程式](/math/pde/laplace/) で扱う。

Fourier 乗数を直接確認するには $P(x)=n x_jx_k-\delta_{jk}|x|^2$ を使う。$P$ は斉次調和二次式であり、:ref[prop-four-harmonic-polynomial-kernel] の定数を $s_{n-1}=2\pi^{n/2}/\Gamma(n/2)$ で割ると $-1/n$ である。このため $\operatorname{pv}(K_{jk})$ の乗数は $-\xi_j\xi_k/|\xi|^2+\delta_{jk}/n$ であり、局所項を引くと $R_jR_k$ の乗数になる。$K_{jk}$ は滑らかな平均ゼロの斉次核なので :ref[thm-four-homogeneous-cz] の弱型評価および :ref[thm-four-cz-cotlar] の主値の結論を持つ。恒等作用素との有限和に Chebyshev を使えば弱型評価は保たれる。
:::

:::example[二次元で二つの対角成分を足す]
$n=2$ では $K_{11}=(x_1^2-x_2^2)/(2\pi|x|^4)$、$K_{22}=-K_{11}$ である。主値を足すと消えるが、局所項は各成分 $-f/2$ なので和は $-f$ である。$K_{12}=x_1x_2/(\pi|x|^4)$ には局所項がなく、例えば $(x_1,x_2)=(1,1)$ で $K_{12}=1/(4\pi)$ である。原点外の計算だけから $R_1^2+R_2^2=0$ と考えると、この局所項を落とすことになる。
:::

:::example[二階微分の強型 $(1,1)$ が破れる]
$n=2$ とし、滑らかな非負の $f$ を $B(0,1)$ に台を持つ質量 $1$ の関数に取る。$|x|>2$ では
$$
R_1^2f(x)=K_{11}(x)+O(|x|^{-3}),
$$
なぜなら $\sup_{0\le t\le1}|\nabla K_{11}(x-ty)|\le C|x|^{-3}$ を使い、$\int |y|f(y)dy\le1$ とできるからである。角度 $\theta$ が $|\theta|\le\pi/6$ なら $K_{11}(x)=\cos(2\theta)/(2\pi|x|^2)\ge1/(4\pi|x|^2)$ である。十分遠いところでは誤差はその半分より小さい。この扇形で絶対値を積分すると半径積分は $\int^\infty dr/r$ となる。従って滑らかな $L^1$ 入力でも出力が $L^1$ に入らないことがある。
:::

### 内部評価から一次導関数を除く

切断を一回使った評価には、右辺に $\nabla u$ が現れた。一次差分でこの項を $u$ と $D^2u$ に分け、後者に小さい係数を掛ける。球を少しずつ広げながら上限を取ると、その係数を左辺に吸収できる。

:::lemma[局所の一次差分評価]{#lem-four-elliptic-difference-interpolation}
$u$ を $B(x_0,r+h)$ の近傍で滑らかとし、$1\le p\le\infty$ とする。$h>0$ に対し
$$
\|\partial_j u\|_{L^p(B(x_0,r))}
\le 2h^{-1}\|u\|_{L^p(B(x_0,r+h))}
+\frac h2\|\partial_j^2 u\|_{L^p(B(x_0,r+h))}.
$$
:::

:::proof
線分 $x+te_j$（$0\le t\le h$）に線形の Taylor 公式を積分余項付きで使い、
$$
\partial_j u(x)=\frac{u(x+he_j)-u(x)}h-\frac1h\int_0^h(h-t)\partial_j^2u(x+te_j)dt
$$
とする。$x\in B(x_0,r)$ の移動は全て大きい球の中に入る。移動後の関数のノルムをその球で評価し、時間積分に Minkowski を使えば、最後の係数は $h^{-1}\int_0^h(h-t)dt=h/2$ である。$p=\infty$ では上限を積分の中へ移せる。
:::

:::theorem[内部 $W^{2,p}$ 評価]{#thm-four-elliptic-interior-w2p}
$A$ を正定値実対称定数行列、$1<p<\infty$ とする。$r>0$、$u,f\in L^p(B(x_0,2r))$ が分布の意味で $L_Au=f$ を満たすなら $u\in W^{2,p}(B(x_0,r))$ であり、
$$
\|D^2u\|_{L^p(B(x_0,r))}
\le C_{n,p,A}\bigl(\|f\|_{L^p(B(x_0,2r))}+r^{-2}\|u\|_{L^p(B(x_0,2r))}\bigr)
$$
である。$W^{2,p}$ は次数 $0,1,2$ の弱微分がすべて $L^p$ にある空間である。
:::

:::proof
まず $u$ は閉球の近傍で滑らかとし、中心を $x_0$ に固定して半径を添字にする。$0<\rho<R$、$t=(R-\rho)/3$ とする。$B_\rho$ で $1$、$B_{\rho+t}$ に台を持つ切断が作れる。例えば球の指示関数 $1_{B_{\rho+t/2}}$ に半径 $t/4$ の滑らかな質量 $1$ の核を畳み込めばよい。階数 $j$ の導関数の上限は $C_jt^{-j}$ 以下である。積の公式と :ref[thm-four-elliptic-lp] から
$$
\|D^2u\|_{p,B_\rho}\le C\big(\|f\|_{p,B_R}
+t^{-1}\|\nabla u\|_{p,B_{\rho+t}}
+t^{-2}\|u\|_{p,B_R}\big)
$$
を得る。:ref[lem-four-elliptic-difference-interpolation] を $h=\eta t$（$0<\eta<1$）で使う。一次項に掛かる $t^{-1}$ を含めると、導関数の寄与は $C\eta\|D^2u\|_{p,B_{\rho+2t}}$、$u$ の寄与は $C\eta^{-1}t^{-2}\|u\|_{p,B_R}$ 以下である。$M=\sup_{0<\rho<R}(R-\rho)^2\|D^2u\|_{p,B_\rho}$ と置く。$\rho+2t$ と $R$ の差は $(R-\rho)/3$ なので、前者を $(R-\rho)^2$ 倍すると $9C\eta M$ 以下になる。よって
$$
M\le C_\eta(R^2\|f\|_{p,B_R}+\|u\|_{p,B_R})+9C\eta M.
$$
滑らかな場合は $M<\infty$ である。$\eta$ を小さく固定して最後の項を左に移し、$\rho=R/2$ を選ぶ。さらに同じ上限より $\rho/R$ を任意の固定比にしても定数が変わるだけである。

一般の場合は滑らかな質量 $1$ の非負核 $\phi_\eps$ を使い、$u_\eps=\phi_\eps*u$、$f_\eps=\phi_\eps*f$ を内側の球で定義する。$u,f$ は大きい球の外でゼロとして積分できるが、方程式は境界から $\eps$ 以上離れた球だけで使う。定数係数なのでそこでは $L_Au_\eps=f_\eps$ である。$R=3r/2$ を選び、滑らかな場合の評価と Young の不等式から $D^2u_\eps$ は $B_r$ で一様に有界である。例えば中間半径 $5r/4$ でも一様に有界で、一次差分評価を $h=r/4$ で使うと一次導関数も $B_r$ で有界となる。$L^p$ の反射性で各成分から弱収束部分列を選び、試験関数で部分積分して $u_\eps\to u$ の $L^p$ 収束を使うと、それらが $u$ の弱微分であることが分かる。ノルムの弱下半連続性から二階の評価も保たれる。
:::

:::example[区間で低次の項が必要な理由]
一次元、$L_A=-\partial_x^2$、$u(x)=x$ とすると $f=0$、$D^2u=0$ だが一次導関数は $1$ である。よって二階の評価で $\nabla u$ を右辺から除いても、一次導関数を $f$ のみで抑えられるわけではない。また二次元では $u(x)=x_1^2-x_2^2$、$f=0$ でも $D^2u\ne0$ である。$B(0,r)$ で $\|\partial_{11}u\|_2=2\sqrt\pi\,r$ であり、$r^{-2}\|u\|_{L^2(B(0,2r))}$ は同じ $r$ の冪に比例する。実際 $u(\rho,\theta)=\rho^2\cos2\theta$ なのでノルムの二乗は $\pi(2r)^6/6$ である。この尺度で低次の項が必要になる。
:::

### 連続係数を小球ごとに固定する

:::proposition[係数の固定による局所評価]{#prop-four-elliptic-freeze}
$A(x)$ は $B(x_0,2r)$ 上で実対称、$a_0 I\le A(x)\le a_1 I$（$0<a_0\le a_1$）とする。$L u=-\sum_{j,k}a_{jk}(x)\partial_{jk}u$ とし、$\max_{j,k}|a_{jk}(x)-a_{jk}(x_0)|\le\delta$ が球内で成り立つとする。$1<p<\infty$ を固定する。$\delta$ が $n,p,a_0,a_1$ のみに依存する値以下なら、滑らかな $u$ に対し
$$
\|D^2u\|_{p,B_r}\le C\bigl(\|Lu\|_{p,B_{2r}}+r^{-2}\|u\|_{p,B_{2r}}\bigr)
$$
である。連続係数なら各点の十分小さい球でこの振動の仮定を満たす。
:::

:::proof
$A_0=A(x_0)$ を固定する。:ref[thm-four-elliptic-lp] の定数は $A_0$ の最小固有値を下から抑えれば一様となる。実際、証明で使う座標変換の後の両ノルムに共通の Jacobian が付き、元へ戻す際の $\|A_0^{-1/2}\|^2$ は $1/a_0$ 以下である。$v$ の台を球内に取ると
$$
\|D^2v\|_p\le C_0\|L_{A_0}v\|_p
\le C_0\|Lv\|_p+C_0\delta\sum_{j,k}\|\partial_{jk}v\|_p .
$$
二階ノルムを成分ノルムの和で取っているので $\delta\le1/(2C_0)$ なら最後の項を左辺に移せる。従って球内に台を持つ $v$ に対し $\|D^2v\|_p\le2C_0\|Lv\|_p$ である。

ここで $v=\chi u$ とする。作用素は非発散形なので、積の微分には係数の導関数が現れない。$Lu$ の項、$-2\sum a_{jk}(\partial_j\chi)(\partial_k u)$ の項、$(L\chi)u$ の項が出る。$|a_{jk}|\le a_1$ と切断の微分評価を使えば、:ref[thm-four-elliptic-interior-w2p] の最初の滑らかな場合の評価と同じ形になる。そこにある一次差分と $M$ の吸収を、$f=Lu$、外側の半径 $R=3r/2$ として使えば結論となる。この半径では $u$ の滑らかさを閉球の近傍で使え、固定した比 $r/R=2/3$ に対して評価が得られる。連続性は、所定の $\delta$ に対して球の半径を十分小さく選べることを保証する。
:::

:::example[係数の変化の大きさと球の半径]
$L=-(1+x_1^2)\partial_1^2-\partial_2^2$ を原点の周りで考える。半径 $2r$ の球上で $A(0)=I$ との差は高々 $4r^2$ である。従って $\delta$ を上の吸収定数とすると、$r\le \sqrt{\delta}/2$ を選べばよい。この例で係数が微分できることは評価に使っておらず、$x_1^2$ を $|x_1|$ に替えても振動の上界 $2r$ を使って小球で同じ手順を取れる。
:::

### 零次項を加えて全空間で逆を作る

:ref[thm-four-elliptic-lp] は二階導関数の大きさを抑えるが、$-\Delta$ の逆の乗数 $1/(4\pi^2|\xi|^2)$ は原点で発散する。正の $\lambda$ を足すとこの特異性が消える。熱半群の時刻積分を用いると、零次項と微分項を同時に制御できる。

:::proposition[$\lambda-\Delta$ の $L^p$ 逆作用素]{#prop-four-elliptic-resolvent}
$\lambda>0$、$1<p<\infty$、$f\in L^p(\R^n)$ とする。
$$
u=\int_0^\infty e^{-\lambda t}G_t*f\,dt
$$
は $L^p$ で収束し、$u\in W^{2,p}$ で $(\lambda-\Delta)u=f$ を満たす。一意性は $L^p$ の分布解全体で成立し、
$$
\lambda\|u\|_p+\lambda^{1/2}\|\nabla u\|_p+\|D^2u\|_p\le C_{n,p}\|f\|_p
$$
である。
:::

:::proof
熱核 $G_t$ は質量 $1$ なので $\lambda\|u\|_p\le\|f\|_p$ である。:ref[prop-four-heat-smoothing] の核評価から
$$
\|\partial_j u\|_p\le C_n\int_0^\infty e^{-\lambda t}t^{-1/2}dt\,\|f\|_p
=C_n\sqrt{\pi/\lambda}\|f\|_p.
$$
従って弱微分としても積分と微分を交換できる。Schwartz 関数の場合は、変換で $u$ の乗数が $(\lambda+4\pi^2|\xi|^2)^{-1}$ である。よって $(\lambda-\Delta)u=f$、さらに $\partial_{jk}u=R_jR_k(f-\lambda u)$ となる。右辺の $L^p$ 有界性から二階評価を得る。一般の $f$ は Schwartz 関数で $L^p$ 近似し、上の各評価を近似の差に使うと、$u$ と一次、二次の導関数は全て $L^p$ で Cauchy である。部分積分で極限が弱微分であることを確かめれば存在が従う。

一意性のため $(\lambda-\Delta)v=0$ の $v\in L^p$ を Fourier 変換する。$v$ は Hölder により緩増加超関数とみなせる。変換後に $(\lambda+4\pi^2|\xi|^2)\widehat v=0$ である。任意の Schwartz 試験関数 $\phi$ をこの正の多項式で割ったものも Schwartz に属する。従って $\widehat v$ は全ての $\phi$ に対してゼロであり、逆変換から $v=0$ となる。
:::

:::example[一次元の逆作用素を空間で計算する]
$n=1$、$\lambda=4$ では核は $E(x)=e^{-2|x|}/4$ である。$x\ne0$ で $4E-E''=0$、$E'(0+)=-1/2$、$E'(0-)=1/2$ なので分布の二階微分に $-\delta_0$ が現れる。従って $(4-\partial_x^2)E=\delta_0$ である。質量は $2(1/4)\int_0^\infty e^{-2x}dx=1/4$ となり、$f\in L^p$ に対する畳み込み $E*f$ は上の一意性により熱半群の積分と同じになる。
:::

## 演習問題

::::exercise{level=basic}
$\xi=(3,4)$ で $R_1R_2$ の乗数を求めよ。
:::solution
二つの乗数は $-3i/5,-4i/5$ なので積は $-12/25$ である。
:::
::::

::::exercise{level=basic}
一次元で $|x|>2|y|$ のとき $|1/(x-y)-1/x|$ を評価せよ。
:::solution
差は $y/(x(x-y))$ である。$|x-y|>|x|/2$ なので絶対値は $2|y|/|x|^2$ 以下である。
:::
::::

::::exercise{level=basic}
$\int_{1<|x|<3}x_1/|x|^{n+1}\,dx$ を求めよ。
:::solution
領域は $x_1\mapsto-x_1$ で不変で、被積分関数は符号を変える。積分は絶対収束する環上で行うので答えはゼロである。
:::
::::

::::exercise{level=basic}
$L_A=-9\partial_1^2-4\partial_2^2$ を Laplace 作用素へ移す座標変換を求めよ。
:::solution
$A^{1/2}=\diag(3,2)$ なので $v(y)=u(3y_1,2y_2)$ と置く。$\Delta v=9u_{11}+4u_{22}$ を $(3y_1,2y_2)$ で評価した値となる。
:::
::::

::::exercise{level=advanced}
$\sum_j\|R_j f\|_2^2=\|f\|_2^2$ を証明せよ。
:::solution
Plancherel から左辺は $\int(\sum_j|\xi_j|^2/|\xi|^2)|\widehat f|^2$ である。原点以外で和は $1$、原点は零測度である。
:::
::::

::::exercise{level=advanced}
$T$ の核が Hörmander 条件を満たす。台が $B(c,r)$ に含まれ、積分がゼロの $b\in L^2$ について $\int_{|x-c|>2r}|Tb|\le C\|b\|_1$ を示せ。
:::solution
台の外では核表示を使える。$\int b=0$ なので $K(x,c)$ を引く。各 $y$ に対して積分範囲は $|x-c|>2|y-c|$ に含まれる。Tonelli と Hörmander 条件で結論を得る。
:::
::::

::::exercise{level=advanced}
$n=2$ の斉次核で $\Omega(\cos\theta,\sin\theta)=1+\cos(2\theta)$ とする。定数成分を除く必要を、原点で一定の試験関数を用いて示せ。
:::solution
$\phi=1$ が $|x|\le1$ で成り立つよう選ぶ。$\eps<|x|<1$ 上の積分は $\int_\eps^1dr/r\int_0^{2\pi}(1+\cos(2\theta))d\theta=2\pi\log(1/\eps)$ となる。外側は有限であるため主値は発散する。平均 $1$ を引けばこの発散項は消える。
:::
::::

::::exercise{level=advanced}
$u\in\mathcal S$ に対し $\sum_{j,k}\|\partial_j\partial_k u\|_2^2=\|\Delta u\|_2^2$ を証明せよ。
:::solution
左辺を変換すると重みは $(2\pi)^4\sum_{j,k}\xi_j^2\xi_k^2=(2\pi)^4|\xi|^4$ である。これは $\Delta$ の乗数の絶対値二乗に等しいので Plancherel で結論を得る。
:::
::::

::::exercise{level=advanced}
二次元で核 $K_{12}(x)=x_1x_2/(\pi|x|^4)$ の球面平均を求め、切断領域を $a(\omega)=e^{4\omega_1\omega_2}$ で変更したときに加わる局所項を求めよ。
:::solution
$K_{12}$ の角度の係数は $\Omega(\omega)=\omega_1\omega_2/\pi$ である。第 1 座標の反射により平均はゼロである。:ref[prop-four-cz-truncation-shape] で $\Omega\log a=(4/\pi)\omega_1^2\omega_2^2$ を積分する。$\omega=(\cos\theta,\sin\theta)$ でその二乗積は $(1-\cos4\theta)/8$ であり、積分値は $\pi/4$ だから係数は $1$ である。式の符号により新しい切断極限は元の主値から $f(x)$ を引いたものになる。
:::
::::

::::exercise{level=advanced}
$f\in L^2(\R^2)$ に対して、Beurling–Ahlfors 変換を $R_1,R_2$ の積の和で表せ。
:::solution
:ref[prop-four-harmonic-polynomial-kernel] の例で乗数は $(\xi_1-i\xi_2)^2/|\xi|^2$ である。$R_1-iR_2$ の乗数は $-i(\xi_1-i\xi_2)/|\xi|$ だから、変換は $-(R_1-iR_2)^2=-R_1^2+2iR_1R_2+R_2^2$ となる。乗数の積は可換なので中間項の係数は $2i$ である。
:::
::::

## 確認テスト

::::quiz{answer="-1" tol="1e-6"}
一次元で $R_1^2$ の非零周波数における乗数を求めよ。
:::explanation
$R_1$ の乗数は $-i\sgn(\xi)$ であり、二乗して $-1$ になる。
:::
::::

::::quiz
Calderón–Zygmund 分解で悪い部分の遠方を評価する際に使う性質はどれか。

- [x] 各立方体上の積分がゼロである。
- [ ] 各立方体上で常に非負である。
- [ ] 全点で元の関数と同じである。

:::explanation
積分がゼロなので、核の値から立方体中心での核の値を引ける。Hörmander 条件は、この差の遠方の絶対積分を一様に抑える。
:::
::::

::::quiz{answer="1" tol="1e-6"}
$\xi=(0,5)$ で $R_2$ の乗数の絶対値を求めよ。
:::explanation
$\xi_2/|\xi|=5/5=1$ なので乗数は $-i$ であり、絶対値は $1$ である。
:::
::::

::::quiz
滑らかな斉次核 $K(x)=|x|^{-n}\Omega(x/|x|)$ で打ち消しに使う条件はどれか。

- [x] $\Omega$ の球面平均がゼロである。
- [ ] $\Omega$ は必ず奇関数である。
- [ ] $\Omega$ は必ず非負である。

:::explanation
球面平均がゼロなら、環の積分は $\int dr/r$ にそのゼロを掛けた値になる。偶関数 $\Omega(\omega)=\omega_1^2-\omega_2^2$ でも平均ゼロである。
:::
::::

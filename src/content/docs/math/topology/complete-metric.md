---
title: 完備距離空間
description: Cauchy 列がすべて収束する距離空間を調べ、完備化の構成、Baire のカテゴリー定理とその応用、縮小写像の原理と常微分方程式の解の存在、Arzelà–Ascoli の定理を証明する。
notation:
  - symbol: 'D(x,r)'
    meaning: '中心 x、半径 r の閉球 {y : d(x,y) ≤ r}'
  - symbol: '\widehat X'
    meaning: 距離空間 X の完備化
  - symbol: '\iota'
    meaning: 完備化への標準的な等長埋め込み
  - symbol: 'B(S,Y)'
    meaning: 集合 S から距離空間 Y への有界写像全体
  - symbol: 'd_\infty'
    meaning: 写像の一様距離（各点の距離の上限）
  - symbol: 'C(K,Y)'
    meaning: 位相空間 K から距離空間 Y への連続写像全体
  - symbol: '\ell^\infty'
    meaning: 有界な実数列全体
  - symbol: 'c_0'
    meaning: 0 に収束する実数列全体
---

## この章の目標

- Cauchy 列と入れ子の閉集合の両面から完備性を判定し、積空間・有界写像の空間・連続写像の空間の完備性を証明する。
- 任意の距離空間の完備化を Cauchy 列の同値類で構成し、一様連続写像の延長定理から一意性を導く。
- Baire のカテゴリー定理を証明し、連続点の集合、完全距離化可能な部分空間、至る所微分不可能な連続関数の存在に応用する。
- 縮小写像の原理を証明し、常微分方程式の初期値問題の解の存在と一意性（Picard–Lindelöf の定理）を導く。
- Arzelà–Ascoli の定理によって連続関数族の相対コンパクト性を特徴付け、Peano の存在定理に応用する。

[距離空間](/math/topology/metric-spaces/) の章で Cauchy 列を定義し、[コンパクト性](/math/topology/compactness/) の章では距離空間のコンパクト性が「完備かつ全有界」と同値であることを証明した。本章では完備性そのものを主題とする。完備性は位相ではなく距離に依存する性質であるが、それが保証する「極限の存在」は解析学のほとんどすべての存在定理の土台になっている。本章の後半で扱う三つの定理（Baire のカテゴリー定理、縮小写像の原理、Arzelà–Ascoli の定理）は、いずれも完備性を使って「求める対象を近似列の極限として作る」ための道具である。

本章を通じて、距離空間 $(X,d)$ の開球を $B(x,r)=\{y:d(x,y)<r\}$、閉球を $D(x,r)=\{y:d(x,y)\le r\}$ と書く（:ref[def-top-metric-open]）。

## 完備性

実数の連続性は「Cauchy 列は収束する」という形でも述べられた（[実数の連続性](/math/analysis/real-numbers/)、[数列と級数](/math/analysis/sequences-series/)）。$\sqrt2$ を小数で近似する有理数列 $1,\ 1.4,\ 1.41,\ 1.414,\ \dots$ は、先へ進むほど項どうしの差がいくらでも小さくなる。しかし $\Q$ の中に極限はなく、極限を置くには $\Q$ の外に出なければならない。Cauchy 列という概念の利点は、極限の値をあらかじめ知らなくても、項どうしの距離だけで「収束すべき列」を判定できることにある。極限をもたない Cauchy 列が存在するかどうかが、空間に「穴」があるかどうかを表す。

完備とは、任意の Cauchy 列が空間内に極限をもつことである（:ref[def-top-complete]）。例えば実数は完備であり、その閉部分集合も完備である。一方、有理数の中では上の列は収束しない。閉部分集合と完備性の対応は :ref[prop-top-complete-closed] で示した。

### 完備性は距離に依存する

一様同値な距離は完備性を共有する（:ref[prop-top-uniform-equivalence-complete]）。一方、:ref[exm-top-log-metric] で見たように、一つの開区間にも完備な距離と不完備な距離が同じ位相を与える。

:::definition[完全距離化可能性]{#def-top-completely-metrizable}
位相空間 $X$ の位相を定める完備な距離が存在するとき、$X$ は :term[完全距離化可能]{en="completely metrizable" yomi="かんぜんきょりかかのう"} であるという。
:::

上の例により、通常の位相の開区間 $(0,1)$ は（通常の距離では完備でないが）完全距離化可能である。完全距離化可能性は位相的性質であり、同相写像で保たれる。どのような部分空間が完全距離化可能であるかは、Baire の定理の節で完全に決定する（:ref[prop-top-gdelta-complete]、:ref[prop-top-complete-subspace-gdelta]）。

### 入れ子の閉集合による判定

実数の連続性は区間縮小法としても述べられた。距離空間では、区間の代わりに直径が $0$ に近づく閉集合の列を使う。

:::theorem[Cantor の共通部分定理]{#thm-top-cantor-intersection}
距離空間 $X$ について次は同値である。

1. $X$ は完備である。
2. 空でない閉集合の減少列 $F_1\supset F_2\supset\cdots$ で $\operatorname{diam}F_n\to0$ を満たすものは、必ずただ一点を共有する。
:::

:::proof
**(1 ⇒ 2)** 各 $n$ で $x_n\in F_n$ を選ぶ。$m,n\ge k$ なら $x_m\in F_m\subset F_k$、$x_n\in F_n\subset F_k$ なので

$$
d(x_m,x_n)\le\operatorname{diam}F_k
$$

である。$\eps>0$ に対し $\operatorname{diam}F_k<\eps$ となる $k$ を選べば、$m,n\ge k$ で $d(x_m,x_n)<\eps$ となるから、$(x_n)$ は Cauchy 列であり、完備性から極限 $x$ をもつ。$k$ を固定すると、尾部 $(x_n)_{n\ge k}$ は閉集合 $F_k$ の点列で $x$ に収束するから、:ref[prop-top-metric-closed-sequence] により $x\in F_k$ である。これがすべての $k$ で成り立つので $x\in\bigcap_kF_k$ である。共通部分の任意の点 $y$ について $d(x,y)\le\operatorname{diam}F_k\to0$ なので $y=x$ であり、共通点はただ一つである。

**(2 ⇒ 1)** $(x_n)$ を Cauchy 列とし、尾部 $T_n=\{x_m:m\ge n\}$ の閉包を $F_n=\overline{T_n}$ と置く。$F_n$ は空でない閉集合で、$T_n\supset T_{n+1}$ より減少列である。閉包を取っても直径は変わらない（:ref[prop-top-closure-diameter]）から $\operatorname{diam}F_n=\operatorname{diam}T_n$ である。$\eps>0$ に対し、$m,l\ge N$ で $d(x_m,x_l)<\eps$ となる $N$ を取れば、$n\ge N$ で $\operatorname{diam}T_n\le\eps$ となる。従って $\operatorname{diam}F_n\to0$ であり、仮定から共通点 $x$ がある。$x$ と $x_n$ はともに $F_n$ に属するので $d(x_n,x)\le\operatorname{diam}F_n\to0$、すなわち $x_n\to x$ である。
:::

:::example[定理の仮定を外すと]{#exm-top-cantor-hypotheses}
1. 完備な $\R$ で $F_n=[n,\infty)$ は空でない閉集合の減少列だが、共通部分は空である。直径は $\infty$ で、条件 $\operatorname{diam}F_n\to0$ が欠けている。
2. 完備でない $X=(0,1]$ で $F_n=(0,1/n]$ は $X$ の閉集合（$X\cap[0,1/n]$ と書ける）で、直径 $1/n\to0$ だが共通部分は空である。完備性の仮定は外せない。
3. 有界な閉集合でも直径が $0$ に近づかなければ共通点は保証されない。具体例は下の :ref[exm-top-linfty] で与える。
:::

### 積空間と写像空間の完備性

座標ごとの極限から積全体の極限を作る方法は :ref[prop-top-complete-product] で証明した。写像空間では同じ方法に一様な誤差評価を加える。

関数空間の完備性は、一様収束の議論に帰着する。次の補題は Arzelà–Ascoli の定理でも用いる。

:::lemma[一様 Cauchy 列と各点での極限]{#lem-top-uniform-cauchy}
$S$ を集合、$(Y,\rho)$ を距離空間とし、写像の列 $f_n:S\to Y$ が**一様 Cauchy 条件**、すなわち

> 各 $\eps>0$ に対しある $N$ があり、$m,n\ge N$ とすべての $s\in S$ で $\rho(f_m(s),f_n(s))<\eps$

を満たすとする。各点 $s$ で $f_n(s)$ が $Y$ に極限 $f(s)$ をもつなら（$Y$ が完備なら、$(f_n(s))_n$ が Cauchy 列なので常にそうである）、$f_n$ は $f$ に一様収束する。さらに $S$ が位相空間で各 $f_n$ が連続なら、$f$ も連続である。
:::

:::proof
$\eps>0$ に対し、$m,n\ge N$ とすべての $s$ で $\rho(f_m(s),f_n(s))<\eps/2$ となる $N$ を取る。$s$ と $n\ge N$ を固定して $m\to\infty$ とすると、距離の連続性（:ref[lem-top-reverse-triangle]）から

$$
\rho(f(s),f_n(s))=\lim_{m\to\infty}\rho(f_m(s),f_n(s))\le\frac\eps2<\eps
$$

となる。$N$ は $s$ によらないので、$f_n$ は $f$ に一様収束する。連続性は一様極限の連続性（:ref[prop-top-uniform-limit]）による。
:::

:::proposition[有界写像の空間と連続写像の空間]{#prop-top-bounded-maps-complete}
$S$ を集合、$(Y,\rho)$ を距離空間とし、像が有界な写像 $f:S\to Y$ の全体を $B(S,Y)$ と書く。

$$
d_\infty(f,g)=\sup_{s\in S}\rho\bigl(f(s),g(s)\bigr)
$$

（$S=\varnothing$ のときは $0$）は $B(S,Y)$ 上の距離であり、:term[一様距離]{en="uniform metric" yomi="いちようきょり"} と呼ばれる。

1. $Y$ が完備なら $(B(S,Y),d_\infty)$ は完備である。
2. $K$ が位相空間なら、有界な連続写像全体 $BC(K,Y)$ は $B(K,Y)$ の閉集合である。従って $Y$ が完備なら $BC(K,Y)$ も完備である。
3. $K$ がコンパクトなら、すべての連続写像 $K\to Y$ は有界である。従って連続写像全体 $C(K,Y)$ に一様距離 $d_\infty$ が入り、$Y$ が完備なら $(C(K,Y),d_\infty)$ は完備である。
:::

:::proof
**距離であること。** $f,g\in B(S,Y)$ とし、$S\ne\varnothing$ として $s_0\in S$ を固定する。像が有界なので、ある $r,r'$ で $f(S)\subset B(f(s_0),r)$、$g(S)\subset B(g(s_0),r')$ となる。すると任意の $s$ で

$$
\rho(f(s),g(s))\le\rho(f(s),f(s_0))+\rho(f(s_0),g(s_0))+\rho(g(s_0),g(s))<r+\rho(f(s_0),g(s_0))+r'
$$

なので上限は有限である。非負性と対称性は明らかで、$d_\infty(f,g)=0$ ならすべての $s$ で $f(s)=g(s)$ である。三角不等式は、各 $s$ で $\rho(f(s),h(s))\le\rho(f(s),g(s))+\rho(g(s),h(s))\le d_\infty(f,g)+d_\infty(g,h)$ が成り立つので、左辺の上限を取ればよい。

**1.** $(f_n)$ を $d_\infty$ に関する Cauchy 列とすると、$\rho(f_m(s),f_n(s))\le d_\infty(f_m,f_n)$ より一様 Cauchy 条件が成り立つ。$Y$ は完備なので各点の極限 $f(s)$ があり、:ref[lem-top-uniform-cauchy] により一様収束する。特にある $N$ ですべての $s$ について $\rho(f(s),f_N(s))\le1$ となり、$f_N$ の像が有界だから $f$ の像も有界、すなわち $f\in B(S,Y)$ である。さらに一様収束から、任意の $\eps>0$ にある $N'$ があって $n\ge N'$ で $d_\infty(f_n,f)\le\eps$ となる。よって $f_n\to f$ である。

**2.** $f_n\in BC(K,Y)$ が $d_\infty$ で $f\in B(K,Y)$ に収束すれば一様収束であり、:ref[prop-top-uniform-limit] により $f$ は連続である。従って $BC(K,Y)$ は閉集合で、$Y$ が完備なら 1 と :ref[prop-top-complete-closed] により完備である。

**3.** $f:K\to Y$ が連続なら $f(K)$ はコンパクトである（:ref[thm-top-compact-image]）。像が空なら有界であり、空でなければ $y_0\in f(K)$ を固定すると開球 $B(y_0,n)$（$n\ge1$）が $f(K)$ を覆うので、有限個で覆え、$f(K)$ は有界である。従って $C(K,Y)=BC(K,Y)$ であり、2 から結論を得る。
:::

:::example[有界数列の空間 ℓ^∞ と c_0]{#exm-top-linfty}
有界な実数列 $a=(a_i)_{i\ge1}$ の全体 $\ell^\infty=B(\N_{\ge1},\R)$ に $d(a,b)=\sup_i\abs{a_i-b_i}$ を入れると、前命題により完備である。

この空間の閉単位球 $\{a:\sup_i\abs{a_i}\le1\}$ は有界な閉集合だがコンパクトではない。第 $n$ 項だけ $1$ で他が $0$ の列 $e_n$ はすべて閉単位球に属し、$m\ne n$ なら $d(e_m,e_n)=1$ なので、$(e_n)$ は Cauchy 部分列をもたず、収束部分列をもたない。

$F_n=\{e_k:k\ge n\}$ と置くと、$F_n$ の点列が収束すればその列は Cauchy 列で、異なる二点の距離が $1$ だからやがて定数になる。従って極限は $F_n$ に属し、$F_n$ は閉である。$F_n$ は空でない有界閉集合の減少列で直径はすべて $1$ だが、共通部分は空である（:ref[exm-top-cantor-hypotheses] の 3）。

$0$ に収束する数列の全体 $c_0$ は $\ell^\infty$ の閉集合である。実際、$a^{(k)}\in c_0$ が $a\in\ell^\infty$ に収束するとし、$\eps>0$ を取る。$d(a,a^{(k)})<\eps/2$ となる $k$ を固定し、$i\ge N$ で $\abs{a^{(k)}_i}<\eps/2$ となる $N$ を取れば、$i\ge N$ で $\abs{a_i}\le\abs{a_i-a^{(k)}_i}+\abs{a^{(k)}_i}<\eps$ となる。よって $a\in c_0$ であり、:ref[prop-top-complete-closed] から $c_0$ も完備である。
:::

積分距離での対照例は :ref[exm-top-l1-incomplete] を参照せよ。

## 完備化

$\Q$ を $\R$ に埋め込めば、$\Q$ の Cauchy 列はすべて $\R$ で極限をもつ。同じことを任意の距離空間で行いたい。ただし一般の距離空間には $\R$ のような「外側の空間」が用意されていない。そこで、収束すべき Cauchy 列そのものを新しい点とみなし、同じ点に収束すべき列を同一視する。

:::definition[完備化]{#def-top-completion}
距離空間 $(X,d)$ の :term[完備化]{en="completion" yomi="かんびか"} とは、完備距離空間 $(\widehat X,\widehat d)$ と等長写像 $\iota:X\to\widehat X$ の組であって、$\iota(X)$ が $\widehat X$ で稠密なものをいう。
:::

元の点は等長に保たれ、加わるのは $\iota(X)$ の点の極限として必要な点だけである。次の補題は、新しい空間の完備性を稠密部分集合の Cauchy 列だけで確かめてよいことを述べる。

:::lemma[稠密部分集合の Cauchy 列による判定]{#lem-top-dense-cauchy}
$E$ を距離空間 $Z$ の稠密部分集合とする。$E$ の点からなる Cauchy 列がすべて $Z$ で収束するなら、$Z$ は完備である。
:::

:::proof
$(z_k)$ を $Z$ の Cauchy 列とする。稠密性から $e_k\in E$ を $d(z_k,e_k)<1/k$ となるように選ぶ。

$$
d(e_k,e_l)\le d(e_k,z_k)+d(z_k,z_l)+d(z_l,e_l)<\frac1k+d(z_k,z_l)+\frac1l
$$

なので $(e_k)$ は Cauchy 列であり、仮定からある $z\in Z$ に収束する。$d(z_k,z)\le1/k+d(e_k,z)\to0$ より $z_k\to z$ である。
:::

:::theorem[距離空間の完備化の存在]{#thm-top-completion-exists}
任意の距離空間 $(X,d)$ は完備化をもつ。$X$ の Cauchy 列の同値類を新しい点とすることで、標準的な完備化を構成できる。
:::

:::proof
$X=\varnothing$ なら空空間自身が完備化である。以下 $X\ne\varnothing$ とし、$X$ の Cauchy 列全体を $\mathcal C$ と書く。

**第 1 段（列の間の距離）。** $a=(a_n),\ b=(b_n)\in\mathcal C$ に対し $s_n=d(a_n,b_n)$ と置く。逆三角不等式（:ref[lem-top-reverse-triangle]）から

$$
\abs{s_m-s_n}\le d(a_m,a_n)+d(b_m,b_n)
$$

であり、右辺は $m,n\to\infty$ で $0$ に近づくので、$(s_n)$ は $\R$ の Cauchy 列である。$\R$ は完備なので極限

$$
\rho(a,b)=\lim_{n\to\infty}d(a_n,b_n)
$$

が存在する。$\rho\ge0$、$\rho(a,a)=0$、$\rho(a,b)=\rho(b,a)$ は明らかであり、各 $n$ での三角不等式 $d(a_n,c_n)\le d(a_n,b_n)+d(b_n,c_n)$ の極限を取れば $\rho(a,c)\le\rho(a,b)+\rho(b,c)$ を得る。

**第 2 段（同一視と距離）。** $\rho(a,b)=0$ のとき $a\sim b$ と書く。反射律と対称律は明らかで、推移律は $\rho$ の三角不等式から従うので、$\sim$ は $\mathcal C$ 上の同値関係である。商集合を $\widehat X=\mathcal C/{\sim}$、$a$ の類を $[a]$ と書く。$a\sim a'$、$b\sim b'$ なら、三角不等式から $\rho(a,b)\le\rho(a,a')+\rho(a',b')+\rho(b',b)=\rho(a',b')$ であり、$a,b$ と $a',b'$ の役割を入れ替えれば逆向きの不等式も得る。従って

$$
\widehat d([a],[b])=\rho(a,b)
$$

は代表元の取り方によらずに定まる。$\widehat d$ は第 1 段の性質を受け継ぎ、さらに $\widehat d([a],[b])=0$ は $a\sim b$、すなわち $[a]=[b]$ を意味するから、$\widehat X$ 上の距離である。

**第 3 段（元の点の埋め込みと稠密性）。** $x\in X$ に定数列 $(x,x,\dots)$ の類を対応させる写像を $\iota$ とする。$\widehat d(\iota(x),\iota(y))=\lim_nd(x,y)=d(x,y)$ なので $\iota$ は等長写像である。次に $p=[a]\in\widehat X$ と $\eps>0$ を取る。$m,n\ge N$ で $d(a_m,a_n)<\eps/2$ となる $N$ を選ぶと、$k\ge N$ のとき

$$
\widehat d(\iota(a_k),p)=\lim_{n\to\infty}d(a_k,a_n)\le\frac\eps2<\eps
$$

である。従って $\iota(X)$ は $\widehat X$ で稠密であり、同時に

$$
\iota(a_k)\longrightarrow[a]\qquad(k\to\infty) \label{eq-top-completion-density}
$$

が示された。

**第 4 段（完備性）。** :ref[lem-top-dense-cauchy] により、$\iota(X)$ の Cauchy 列が収束することを示せばよい。そのような列は $\iota(a_k)$ と書け、$\iota$ が等長なので $(a_k)$ は $X$ の Cauchy 列である。:eqref[eq-top-completion-density] により $\iota(a_k)\to[(a_k)]$ であるから、列は $\widehat X$ で収束する。よって $\widehat X$ は完備であり、$(\widehat X,\iota)$ は $X$ の完備化である。
:::

:::remark[実数の構成との関係]
この構成では距離の値が実数であり、第 1 段で $\R$ の完備性を用いた。従ってこの定理を $X=\Q$ に適用して $\R$ を構成することはできない（循環論法になる）。$\Q$ から $\R$ を構成するには、距離の値も有理数の Cauchy 列で扱う必要がある。その構成は [実数と複素数の構成](/math/sets-logic/number-systems/) で行う。
:::

完備化の一意性と、完備化への写像の延長は、次の延長定理から従う。

:::theorem[稠密部分からの一様連続延長]{#thm-top-completion-extension}
$E$ を距離空間 $(X,d)$ の稠密部分集合、$(Y,\rho)$ を完備距離空間とする。一様連続写像 $f:E\to Y$ は連続写像 $F:X\to Y$ にただ一通りに延長でき、$F$ も一様連続である。$f$ が等長なら $F$ も等長である。
:::

:::proof
**値の定義。** $x\in X$ を取る。稠密性から $a_n\in E\cap B(x,1/n)$ を選べるので、$x$ に収束する $E$ の点列が存在する。そのような列 $(a_n)$ は $X$ の Cauchy 列であり、:ref[prop-top-uniform-equivalence-complete] から $(f(a_n))$ は $Y$ の Cauchy 列である。$Y$ は完備なので $\lim f(a_n)$ が存在する。

この極限は列の選び方によらない。実際、$(b_n)$ も $x$ に収束する $E$ の点列なら、交互に並べた列 $(a_1,b_1,a_2,b_2,\dots)$ も $x$ に収束する $E$ の点列なので、その像は収束する。収束列の部分列は同じ極限をもつ（:ref[prop-top-metric-limit-laws]）ので、$\lim f(a_n)=\lim f(b_n)$ である。そこで $F(x)=\lim f(a_n)$ と定める。$x\in E$ なら定数列 $a_n=x$ を使えるので $F(x)=f(x)$ である。

**一様連続性。** $\eps>0$ に対し、$a,a'\in E$、$d(a,a')<\delta$ なら $\rho(f(a),f(a'))<\eps/2$ となる $\delta>0$ を取る。$x,x'\in X$、$d(x,x')<\delta/3$ とし、$x,x'$ に収束する $E$ の点列 $(a_n),(a'_n)$ を取る。十分大きい $n$ で $d(a_n,x)<\delta/3$、$d(a'_n,x')<\delta/3$ となるので、$d(a_n,a'_n)<\delta$、従って $\rho(f(a_n),f(a'_n))<\eps/2$ である。$n\to\infty$ とすれば距離の連続性から $\rho(F(x),F(x'))\le\eps/2<\eps$ となる。よって $F$ は一様連続である。

**一意性。** $G:X\to Y$ が連続で $G|_E=f$ なら、$x$ に収束する $E$ の点列 $(a_n)$ について $G(x)=\lim G(a_n)=\lim f(a_n)=F(x)$ である。

**等長性。** $f$ が等長なら、$x,x'$ に収束する $E$ の点列について $\rho(F(x),F(x'))=\lim\rho(f(a_n),f(a'_n))=\lim d(a_n,a'_n)=d(x,x')$ である。最後の等号は距離の連続性による。
:::

:::corollary[完備化の一意性]{#cor-top-completion-unique}
$(Z,i)$ と $(W,j)$ が $X$ の二つの完備化なら、$\Phi\circ i=j$ を満たす全射等長写像 $\Phi:Z\to W$ がただ一つ存在する。
:::

:::proof
$i(X)$ 上の写像 $i(x)\mapsto j(x)$ は等長であり（$i$ は単射なので well-defined である）、$i(X)$ は $Z$ で稠密、$W$ は完備である。前定理により等長写像 $\Phi:Z\to W$ に延長され、$\Phi\circ i=j$ を満たす。同様に $j(x)\mapsto i(x)$ を延長した等長写像 $\Psi:W\to Z$ がある。$\Psi\circ\Phi:Z\to Z$ は連続で、稠密部分集合 $i(X)$ 上で恒等写像に一致するから、延長の一意性により $\Psi\circ\Phi=\id_Z$ である。同様に $\Phi\circ\Psi=\id_W$ なので、$\Phi$ は全射である。$\Phi'$ も同じ条件を満たせば、$\Phi'$ と $\Phi$ は $i(X)$ 上で一致する連続写像なので等しい。
:::

この一意性があるので、完備化を具体的に求めるには、条件を満たす組を一つ見つければよい。

:::proposition[完備な空間の中での完備化]{#prop-top-completion-closure}
$X$ が完備距離空間 $Z$ の部分空間なら、$X$ の $Z$ での閉包 $\overline X$ と包含写像の組は $X$ の完備化である。
:::

:::proof
$\overline X$ は完備な $Z$ の閉集合なので完備である（:ref[prop-top-complete-closed]）。包含写像は等長で、$X$ は $\overline X$ で稠密である。
:::

:::example[完備化の具体例]{#exm-top-completion-examples}
1. $\Q$ の通常の距離による完備化は $\R$ である。$\Q\subset\R$ で $\overline\Q=\R$ だからである。同様に $\Q\cap[0,1]$ の完備化は $[0,1]$ である。
2. 平面の開円板 $\{z:\norm z<1\}$ の完備化は閉円板 $\{z:\norm z\le1\}$ である。境界の点 $p$ は $(1-1/n)p$ の極限として閉包に入る。
3. 完備な空間の完備化はそれ自身である。例えば無限離散距離空間には一点も追加されない。$(0,1)$ に :ref[exm-top-log-metric] の完備距離 $\rho$ を入れた空間も、それ自身が完備化である。一方、同じ $(0,1)$ に通常の距離を入れると完備化は $[0,1]$ である。完備化は位相ではなく距離によって決まる。
4. $\R$ に $d(x,y)=\abs{\arctan x-\arctan y}$ を入れると、$\arctan$ は $(\R,d)$ から $(-\pi/2,\pi/2)$（通常の距離）への全射等長写像である。従って完備化は $[-\pi/2,\pi/2]$ と等長であり、二つの点が加わる。列 $n$ と $-n$ は $d$ に関する Cauchy 列で、それぞれ異なる追加点に収束する。$\R$ の一点コンパクト化（:ref[thm-top-one-point]）では両方向の無限遠が一点になるので、二つの構成は異なる。
5. 有限個を除いて $0$ である実数列の全体 $c_{00}$ に上限距離を入れる。$c_{00}\subset c_0$ で、$c_0$ は完備である（:ref[exm-top-linfty]）。$a\in c_0$ の第 $N$ 項より後を $0$ に置き換えた列 $a^{[N]}\in c_{00}$ について $d(a,a^{[N]})=\sup_{i>N}\abs{a_i}\to0$ なので、$c_{00}$ は $c_0$ で稠密である。前命題により $c_{00}$ の完備化は $c_0$ である。定数列 $(1,1,\dots)$ は $\ell^\infty$ に属するが、どの $c_{00}$ の元とも距離が $1$ 以上なので完備化に入らない。
:::

:::example[連続性だけでは延長できない]
$f:(0,1]\to\R$、$f(x)=1/x$ は連続であるが、完備化 $[0,1]$ への連続延長をもたない。$1/n\to0$ に対し像 $n$ は Cauchy 列ですらないからである。これは $f$ が一様連続でないことと対応している。一方 $g(x)=x^2$ は $\abs{x^2-y^2}=\abs{x+y}\abs{x-y}\le2\abs{x-y}$ を満たすので一様連続であり、延長は $G(0)=0$ で与えられる。
:::

:::remark[ノルム空間の完備化]
ノルム空間 $V$ の完備化では $[u_n]+[v_n]=[u_n+v_n]$、$a[u_n]=[au_n]$、$\norm{[u_n]}=\lim_n\norm{u_n}$ と定める。和は三角不等式、固定した $a$ による積は $\norm{au-av}=\abs a\norm{u-v}$、ノルムは逆三角不等式により代表元によらない。各演算の公理は項ごとに成り立ち、$\norm{[u_n]}=0$ は $[u_n]=[0]$ と同値である。従って完備化は Banach 空間になる。なお二変数の積 $(a,u)\mapsto au$ は一般に一様連続でない。$\R^2$ の点 $(n,0)$ と $(n,1/n)$ の距離は $0$ に近づくが、積の差は $1$ である。
:::

## Baire のカテゴリー定理

集合の「小ささ」にはいくつかの測り方がある。濃度による測り方（可算集合は小さい）や、長さによる測り方（総長がいくらでも小さい区間の族で覆える集合は小さい）のほかに、位相による測り方がある。それは「どの開集合を取っても、その中に集合を避ける開集合がある」という意味での小ささである。Baire のカテゴリー定理は、完備距離空間がこの意味で小さい集合の可算個では覆えないことを主張する。

疎集合、第一類集合、Baire 空間の定義には :ref[def-top-category] を用いる。疎集合を避けるには、任意の空でない開集合の中でさらに空でない開集合を選び直せばよい（:ref[lem-top-nowhere-dense]）。

:::proposition[Baire 空間の言い換え]{#prop-top-baire-equivalent}
位相空間 $X$ について次は同値である。

1. $X$ は Baire 空間である。
2. 内部が空である閉集合の可算和は、内部が空である。
3. 空でない開集合は第一類でない。
:::

:::proof
**(1 ⇔ 2)** 閉集合 $F_n$ の補集合 $G_n=X\setminus F_n$ を考える。$\operatorname{int}F_n=\varnothing$ は $G_n$ が稠密であることと同値であり、$X\setminus\bigcup F_n=\bigcap G_n$ である。$\bigcup F_n$ の内部が空であることは、$\bigcap G_n$ が稠密であることと同値である。この対応で 1 と 2 は互いに言い換えになる。

**(2 ⇒ 3)** 空でない開集合 $O$ が疎集合 $A_n$ の可算和であるとする。$O\subset\bigcup\overline{A_n}$ であり、各 $\overline{A_n}$ は内部が空な閉集合だから、2 により $\bigcup\overline{A_n}$ の内部は空である。しかし空でない開集合 $O$ を含むので矛盾である。

**(3 ⇒ 2)** 内部が空な閉集合 $F_n$ の和が空でない開集合 $O$ を含むとすると、$O=\bigcup(O\cap F_n)$ であり、各 $O\cap F_n$ は疎集合 $F_n$ の部分集合なので疎である（:ref[lem-top-nowhere-dense]）。従って $O$ は第一類となり、3 に反する。
:::

:::theorem[Baire のカテゴリー定理]{#thm-top-baire}
次の空間は Baire 空間である。

1. 完備距離空間。
2. 局所コンパクト Hausdorff 空間。

従ってこれらの空間の空でない開集合（特に空でない空間全体）は第一類でない。
:::

:::proof
$(G_n)_{n\ge1}$ を稠密な開集合の列、$O$ を空でない開集合とし、$O\cap\bigcap_nG_n\ne\varnothing$ を示せばよい。

**完備距離空間の場合。** $G_1$ は稠密なので $O\cap G_1$ は空でない開集合である。その点 $x_1$ を取り、$B(x_1,2r_1)\subset O\cap G_1$ となる $0<r_1<1/2$ を選ぶ。$D(x_1,r_1)\subset B(x_1,2r_1)$ なので $D(x_1,r_1)\subset O\cap G_1$ である。

$x_n,r_n$ まで選んだとする。$G_{n+1}$ は稠密なので $B(x_n,r_n)\cap G_{n+1}$ は空でない開集合であり、同様にその点 $x_{n+1}$ と $0<r_{n+1}<2^{-(n+1)}$ を

$$
D(x_{n+1},r_{n+1})\subset B(x_n,r_n)\cap G_{n+1}
$$

となるように選べる。閉球 $D_n=D(x_n,r_n)$ は閉集合で（:ref[prop-top-ball-open]）、$D_{n+1}\subset D_n$、$\operatorname{diam}D_n\le2r_n<2^{1-n}\to0$ である。Cantor の共通部分定理（:ref[thm-top-cantor-intersection]）により共通点 $z$ がある。$z\in D_1\subset O\cap G_1$ であり、$n\ge1$ で $z\in D_{n+1}\subset G_{n+1}$ なので、$z\in O\cap\bigcap_nG_n$ である。

**局所コンパクト Hausdorff 空間の場合。** $O\cap G_1$ は空でない開集合なので、その点 $x_1$ を取ると、:ref[lem-top-compact-neighborhood] により $x_1\in V_1\subset\overline{V_1}\subset O\cap G_1$ を満たし $\overline{V_1}$ がコンパクトな開集合 $V_1$ がある。$V_n$ まで選んだとき、$V_n\cap G_{n+1}$ は空でない開集合なので、その点 $x_{n+1}$ に同じ補題を使い、

$$
x_{n+1}\in V_{n+1}\subset\overline{V_{n+1}}\subset V_n\cap G_{n+1}
$$

で $\overline{V_{n+1}}$ がコンパクトな開集合 $V_{n+1}$ を選ぶ。$\overline{V_1}\supset\overline{V_2}\supset\cdots$ はコンパクト空間 $\overline{V_1}$ の空でない閉集合の減少列であり、有限交叉性をもつから、:ref[thm-top-compact-fip] により共通点 $z$ がある。$z\in\overline{V_1}\subset O\cap G_1$、$z\in\overline{V_{n+1}}\subset G_{n+1}$ なので $z\in O\cap\bigcap_nG_n$ である。

最後の主張は :ref[prop-top-baire-equivalent] による。
:::

:::remark[証明で使われたもの]
完備距離空間の場合の証明では、各段で点と半径を一つずつ選んだ。これは前の選択に依存する可算回の選択であり、従属選択公理を用いている。完備性は最後に Cantor の共通部分定理で一回だけ使われる。$\Q$ では同じ手順で入れ子の閉球を作れても共通点が無理数になりうるため、結論が成り立たない（次の例）。
:::

:::corollary[閉集合による可算被覆]{#cor-top-baire-closed-cover}
完備距離空間 $X$ が閉集合の可算和 $X=\bigcup_nF_n$ で表されるとき、$\bigcup_n\operatorname{int}F_n$ は $X$ で稠密な開集合である。特に $X\ne\varnothing$ なら、ある $F_n$ は内部をもつ。
:::

:::proof
空でない開集合 $O$ が $\bigcup_n\operatorname{int}F_n$ と交わらないとする。$O=\bigcup_n(O\cap F_n)$ であり、各 $O\cap F_n$ が疎であることを示す。$W$ を $\overline{O\cap F_n}$ に含まれる開集合とする。$W\subset F_n$ なので $W\subset\operatorname{int}F_n$ であり、仮定から $W\cap O=\varnothing$ である。一方 $W\subset\overline O$ であるから、$W\ne\varnothing$ なら $W$ は $O$ と交わるはずである。よって $W=\varnothing$ であり、$O\cap F_n$ は疎である。すると空でない開集合 $O$ が第一類となり、:ref[thm-top-baire] に反する。
:::

:::example[Baire の定理の最初の帰結]{#exm-top-baire-first}
1. **$\Q$ は Baire 空間でない**：:ref[exm-top-nowhere-dense] の 5 により、空でない空間 $\Q$ は自分自身の中で第一類である。従って $\Q$ の通常の位相を定める完備距離は存在せず、$\Q$ は完全距離化可能でない。
2. **$\R$ は非可算である**：$\R$ が可算なら一点集合（疎集合）の可算和となり、完備距離空間 $\R$ が第一類になってしまう。同じ議論で、孤立点をもたない空でない完備距離空間は非可算である（演習問題）。
3. **無理数全体は第一類でない**：$\R\setminus\Q$ が第一類なら、$\R=\Q\cup(\R\setminus\Q)$ も第一類となって矛盾する。
:::

:::example[各点で有界な連続関数列]{#exm-top-pointwise-bounded}
空でない完備距離空間 $X$ 上の連続関数列 $f_n:X\to\R$ が、各点 $x$ で $\sup_n\abs{f_n(x)}<\infty$ を満たすとする。

$$
F_m=\bigcap_{n\ge1}\{x:\abs{f_n(x)}\le m\}
$$

は閉集合の共通部分なので閉であり、各点での有界性から $\bigcup_{m\ge1}F_m=X$ である。:ref[cor-top-baire-closed-cover] により、ある $F_m$ は空でない開集合 $W$ を含む。すなわち $W$ 上ではすべての $n$ で $\abs{f_n}\le m$ である。各点で異なってよかった上界から、ある開集合上で一斉に成り立つ上界を取り出したことになる。この議論は関数解析の一様有界性原理の核心である（[関数解析](/math/functional-analysis/)）。
:::

### 残留集合と G_δ 集合

Baire 空間では、第一類集合の補集合が「大きい」集合として振る舞う。

:::definition[残留集合]
可算個の稠密な開集合の共通部分を含む集合を :term[残留集合]{en="residual set" yomi="ざんりゅうしゅうごう"} という。可算個の開集合の共通部分を $G_\delta$ 集合と呼んだ（:ref[prop-top-zero-set]）ので、稠密な $G_\delta$ 集合を含む集合と言っても同じである。
:::

:::proposition[残留集合の性質]{#prop-top-residual}
位相空間 $X$ の部分集合 $A$ について次が成り立つ。

1. $A$ が残留集合であることと、$X\setminus A$ が第一類であることは同値である。
2. 可算個の残留集合の共通部分は残留集合である。
3. $X$ が空でない Baire 空間なら、残留集合は稠密であり、第一類でない。
:::

:::proof
1. $A\supset\bigcap_nG_n$（$G_n$ は稠密な開集合）とする。$X\setminus A\subset\bigcup_n(X\setminus G_n)$ であり、各 $X\setminus G_n$ は内部が空な閉集合なので疎である（:ref[lem-top-nowhere-dense]）。疎集合の部分集合は疎であるから、$X\setminus A=\bigcup_n\bigl((X\setminus A)\cap(X\setminus G_n)\bigr)$ は第一類である。逆に $X\setminus A\subset\bigcup_nA_n$（$A_n$ は疎）なら、$G_n=X\setminus\overline{A_n}$ は稠密な開集合で（:ref[lem-top-nowhere-dense]）、

$$
\bigcap_nG_n=X\setminus\bigcup_n\overline{A_n}\subset X\setminus\bigcup_nA_n\subset A
$$

である。
2. 可算個の稠密な開集合の族を可算個集めても可算個なので、それらの共通部分は可算個の稠密開集合の共通部分を含む。
3. $A\supset\bigcap G_n$ で、Baire 性から右辺は稠密なので $A$ も稠密である。もし $A$ が第一類なら、1 から $X\setminus A$ も第一類なので、$X=A\cup(X\setminus A)$ が第一類となり、:ref[prop-top-baire-equivalent] の 3 に反する。
:::

:::remark[長さによる小ささとは別物である]
有理数を $q_1,q_2,\dots$ と並べ、$n\ge1$ に対し

$$
U_n=\bigcup_{k\ge1}\bigl(q_k-2^{-k-n},\ q_k+2^{-k-n}\bigr)
$$

と置く。$U_n$ は $\Q$ を含む開集合なので稠密であり、$G=\bigcap_nU_n$ は稠密な $G_\delta$ 集合、従って残留集合である。一方 $G\subset U_n$ は総長 $\sum_k2\cdot2^{-k-n}=2^{1-n}$ の区間の族で覆われ、$n$ は任意である。つまり $\R$ は、位相の意味で小さい集合（$\R\setminus G$、第一類）と長さの意味で小さい集合（$G$）の和に分けられる。二つの小ささは互いに独立な概念である。
:::

完全距離化可能な部分空間は、$G_\delta$ 集合として特徴付けられる。

:::proposition[G_δ 部分集合の完全距離化]{#prop-top-gdelta-complete}
完備距離空間 $(X,d)$ の $G_\delta$ 部分集合 $A=\bigcap_{n\ge1}U_n$（$U_n$ は開）は完全距離化可能である。すなわち、$A$ の相対位相を定める完備な距離が存在する。
:::

:::proof
$F_n=X\setminus U_n$ と置く。$F_n\ne\varnothing$ のとき $A$ 上の関数 $h_n(x)=1/d(x,F_n)$ を考え、$F_n=\varnothing$ のときは $h_n=0$ とする。$d(\cdot,F_n)$ は $1$-Lipschitz で連続であり（:ref[prop-top-distance-set]）、$F_n$ が閉なので $d(x,F_n)=0$ と $x\in F_n$ は同値である。$A\cap F_n=\varnothing$ だから分母は $A$ 上で正であり、$h_n$ は $A$ 上の連続関数である。

**埋め込み。** 完備空間の可算積 $P=X\times\R\times\R\times\cdots$ に :ref[prop-top-complete-product] の距離 $D$ を入れると、$P$ は完備で、$D$ は積位相を定める。$\varphi:A\to P$ を

$$
\varphi(x)=\bigl(x,h_1(x),h_2(x),\dots\bigr)
$$

で定める。各座標が連続なので $\varphi$ は連続であり（:ref[thm-top-product-property]）、第一座標を見れば単射である。像の上で逆写像は第一射影の制限なので連続であり、$\varphi$ は $A$ から像 $\varphi(A)$ への同相写像である。

**像が閉であること。** $\varphi(x_k)\to(x,(t_n)_n)$ とする。積位相での収束は座標ごとの収束なので、$x_k\to x$ かつ各 $n$ で $h_n(x_k)\to t_n$ である。もしある $n$ で $x\in F_n$ なら、$0<d(x_k,F_n)\le d(x_k,x)\to0$ となり、$h_n(x_k)\to\infty$ であって実数 $t_n$ に収束できない。従ってすべての $n$ で $x\in U_n$、すなわち $x\in A$ である。$h_n$ の連続性から $t_n=\lim h_n(x_k)=h_n(x)$ となり、極限は $\varphi(x)\in\varphi(A)$ である。距離空間で点列について閉じているので、$\varphi(A)$ は閉集合である（:ref[prop-top-metric-closed-sequence]）。

**完備な距離。** $\varphi(A)$ は完備な $P$ の閉集合なので完備である。$\rho(x,y)=D(\varphi(x),\varphi(y))$ と置くと、$\varphi$ は $(A,\rho)$ から $\varphi(A)$ への全射等長写像なので $\rho$ は完備であり（:ref[prop-top-uniform-equivalence-complete]）、$\varphi$ が同相写像であることから $\rho$ は $A$ の相対位相を定める。
:::

:::proposition[完全距離化可能な部分空間は G_δ 集合である]{#prop-top-complete-subspace-gdelta}
距離空間 $(X,d)$ の部分空間 $A$ が完全距離化可能なら、$A$ は $X$ の $G_\delta$ 集合である。従って完備距離空間の部分空間が完全距離化可能であることと、$G_\delta$ 集合であることは同値である。
:::

:::proof
$\rho$ を $A$ の相対位相を定める完備な距離とし、$\overline A$ を $X$ での閉包とする。$n\ge1$ に対し

$$
V_n=\bigl\{x\in\overline A:\ x\text{ のある開近傍 }W\text{（}X\text{ で開）について }\operatorname{diam}_\rho(W\cap A)<1/n\bigr\}
$$

と置く。$x\in V_n$ に対する $W$ を一つ取れば、$\overline A\cap W$ の各点も同じ $W$ によって $V_n$ に属する。従って $V_n$ は $\overline A$ の相対開集合であり、$X$ の開集合 $O_n$ を用いて $V_n=\overline A\cap O_n$ と書ける。

**$A\subset V_n$。** $a\in A$ とする。$\rho$ の開球 $B_\rho(a,1/(3n))$ は $A$ の相対開集合なので、$X$ の開集合 $W\ni a$ で $W\cap A=B_\rho(a,1/(3n))$ となるものがある。その $\rho$ 直径は $2/(3n)$ 以下なので $a\in V_n$ である。

**$\bigcap_nV_n\subset A$。** $x\in\bigcap_nV_n$ とし、各 $n$ で $\operatorname{diam}_\rho(W'_n\cap A)<1/n$ となる $x$ の開近傍 $W'_n$ を取る。$W_n=W'_1\cap\cdots\cap W'_n\cap B(x,1/n)$ と置けば、$W_n$ は減少する $x$ の開近傍で、$\operatorname{diam}_\rho(W_n\cap A)<1/n$ を満たす。$x\in\overline A$ なので $W_n\cap A\ne\varnothing$ であり、$a_n\in W_n\cap A$ を選べる。$d(a_n,x)<1/n$ なので $a_n\to x$（距離 $d$ で）である。また $m\ge n$ なら $a_m,a_n\in W_n\cap A$ なので $\rho(a_m,a_n)<1/n$ であり、$(a_n)$ は $\rho$ の Cauchy 列である。$\rho$ は完備なので、ある $a\in A$ に $\rho$ で収束する。$\rho$ は $A$ の相対位相を定めるから、$a_n\to a$ は相対位相での収束、すなわち $d$ での収束でもある。極限の一意性から $x=a\in A$ である。

以上より $A=\overline A\cap\bigcap_nO_n$ である。最後に閉集合 $\overline A$ は $\overline A=\bigcap_k\{x:d(x,\overline A)<1/k\}$ と開集合の可算共通部分で書ける（$A=\varnothing$ なら $\overline A=\varnothing$ は開集合である）。従って $A$ は $G_\delta$ 集合である。後半は前命題と合わせればよい。
:::

:::example[無理数の空間と有理数の空間]{#exm-top-irrationals}
有理数を $q_1,q_2,\dots$ と並べると、無理数全体は $\bigcap_n(\R\setminus\{q_n\})$ と書けるので $G_\delta$ 集合であり、:ref[prop-top-gdelta-complete] により完全距離化可能、従って Baire 空間である。通常の距離では有理数に近づく無理数列があるので完備でないが、位相を変えずに完備な距離に取り替えられる。

一方 $\Q$ は完全距離化可能でない（:ref[exm-top-baire-first]）ので、:ref[prop-top-complete-subspace-gdelta] により $\Q$ は $\R$ の $G_\delta$ 集合でない。可算稠密集合を取り除いた残りと、取り除かれた $\Q$ 自身とでは、Baire 性がまったく異なる。
:::

### 各点極限の連続点

連続関数列の各点極限は連続とは限らない。しかし完備距離空間上では、不連続点が「小さい」ことが Baire の定理から従う。

:::corollary[各点収束する連続関数列の連続点]{#cor-top-baire-one}
完備距離空間 $X$ 上の連続関数列 $f_n:X\to\R$ がすべての点で $f$ に収束するとする。このとき $f$ が連続である点全体は稠密な $G_\delta$ 集合である。
:::

:::proof
$\eps>0$ に対し、「任意の二点 $z,w\in U$ で $\abs{f(z)-f(w)}<\eps$」を満たす開集合 $U$ すべての和を $G_\eps$ とする。$G_\eps$ は開集合である。

**$G_\eps$ は稠密である。** 空でない開集合 $O$ を取る。$B(x,2r)\subset O$ となる開球を取り、$V=B(x,r)$、$K=\overline V$ と置くと $K\subset D(x,r)\subset O$ である。$K$ は完備な $X$ の閉集合なので完備であり、空でない。$K$ の閉集合

$$
E_N=\bigcap_{m,n\ge N}\{x\in K:\abs{f_m(x)-f_n(x)}\le\eps/4\}
$$

を考える（$f_m-f_n$ の連続性から各集合は閉である）。各点で $f_n(x)$ は収束するので Cauchy 列であり、どの $x\in K$ もある $E_N$ に属する。$K=\bigcup_NE_N$ に :ref[cor-top-baire-closed-cover] を適用すると、ある $E_N$ は $K$ での相対内部 $W'$ が空でない。$W'=K\cap O'$（$O'$ は $X$ で開）と書ける。$V$ は $K$ で稠密であり $W'$ は $K$ の空でない相対開集合なので、$W=V\cap O'=V\cap W'$ は空でない。$W$ は $X$ の開集合で、$W\subset E_N$ である。

$W$ 上では、$m\ge N$ で $\abs{f_m-f_N}\le\eps/4$ なので、$m\to\infty$ として $\abs{f-f_N}\le\eps/4$ である。$w\in W$ を一つ取り、$f_N$ の $w$ での連続性から、$w$ の開近傍 $W_0\subset W$ ですべての $z\in W_0$ について $\abs{f_N(z)-f_N(w)}<\eps/8$ となるものを取る。$z,z'\in W_0$ なら $\abs{f_N(z)-f_N(z')}<\eps/4$ であり、

$$
\abs{f(z)-f(z')}\le\abs{f(z)-f_N(z)}+\abs{f_N(z)-f_N(z')}+\abs{f_N(z')-f(z')}<\frac{3\eps}4<\eps
$$

となる。よって $W_0\subset G_\eps\cap O$ であり、$G_\eps$ は $O$ と交わる。

**連続点の集合。** $C=\bigcap_{k\ge1}G_{1/k}$ は稠密な開集合の可算共通部分なので、Baire の定理により稠密な $G_\delta$ 集合である。$x\in C$ なら、各 $k$ で $x$ を含む開集合 $U$ 上の値の差が $1/k$ 未満なので、$f$ は $x$ で連続である。逆に $f$ が $x$ で連続なら、$\abs{f(z)-f(x)}<1/(2k)$ を満たす $x$ の開近傍 $U$ 上では二点の差が $1/k$ 未満なので $x\in G_{1/k}$ である。従って $C$ はちょうど $f$ の連続点の集合である。
:::

:::example[各点極限になれる関数となれない関数]
1. 有理数の特性関数 $\chi_\Q$ はどの開区間でも $0$ と $1$ をともに取るので、連続点をもたない。前系により、$\chi_\Q$ は $\R$ 上の連続関数列の各点極限として表せない。
2. $g_n(x)=\max(0,1-n\abs x)$ は連続で、原点では常に $1$、$x\ne0$ では $n>1/\abs x$ 以後 $0$ なので、各点極限は $\chi_{\{0\}}$ である。その連続点の集合 $\R\setminus\{0\}$ は確かに稠密な開集合である。
3. $f:\R\to\R$ が至る所微分可能なら、$f'(x)=\lim_{n\to\infty}n\bigl(f(x+1/n)-f(x)\bigr)$ は連続関数列の各点極限である。従って導関数 $f'$ の連続点は稠密な $G_\delta$ 集合をなす。微分可能な関数の導関数は不連続でありうるが、不連続点は第一類の集合に限られる。
:::

### 至る所微分不可能な連続関数

Baire の定理の最も有名な応用の一つは、「典型的な」連続関数がどこでも微分不可能であることの証明である。具体的な関数を構成する代わりに、微分可能な点をもつ関数の集合が第一類であることを示す。

以下 $C([0,1])$ に一様距離 $d_\infty$ を入れる（:ref[prop-top-bounded-maps-complete] により完備）。$f$ が $x\in[0,1]$ で微分可能であるとは、端点では片側微分係数が有限値として存在することとする。

:::theorem[至る所微分不可能な関数の存在]{#thm-top-nowhere-differentiable}
$C([0,1])$ において、$[0,1]$ のどの点でも微分可能でない関数の全体は残留集合である。すなわち、それは稠密な $G_\delta$ 集合を含む。特に、至る所微分不可能な連続関数が存在する。
:::

:::proof
$n\ge1$ に対し

$$
E_n=\Bigl\{f\in C([0,1]):\ \text{ある }x\in[0,1-\tfrac1n]\text{ があり、すべての }h\in(0,\tfrac1n]\text{ で }\abs{f(x+h)-f(x)}\le nh\Bigr\}
$$

と置く。

**第 1 段（右微分係数をもつ関数は $\bigcup E_n$ に属する）。** $f$ が $x\in[0,1)$ で有限な右微分係数 $c$ をもつとする。ある $\delta\in(0,1-x]$ があって、$0<h<\delta$ で $\abs{f(x+h)-f(x)}\le(\abs c+1)h$ が成り立つ。$\delta\le h\le1-x$ では $\abs{f(x+h)-f(x)}\le2d_\infty(f,0)\le(2d_\infty(f,0)/\delta)\,h$ である。そこで $n\ge\max\{\abs c+1,\ 2d_\infty(f,0)/\delta,\ 1/(1-x)\}$ を満たす整数 $n$ を取れば、$x\le1-1/n$ であり、$h\in(0,1/n]\subset(0,1-x]$ に対しどちらの場合も $\abs{f(x+h)-f(x)}\le nh$ となる。よって $f\in E_n$ である。

**第 2 段（$E_n$ は閉集合である）。** $f_k\in E_n$ が $f$ に一様収束するとし、$f_k$ に対する点 $x_k\in[0,1-1/n]$ を取る。$[0,1-1/n]$ はコンパクトな距離空間なので（:ref[thm-top-heine-borel]、:ref[thm-top-metric-compactness]）、部分列を取り直して $x_k\to x\in[0,1-1/n]$ としてよい。$h\in(0,1/n]$ を固定すると $x_k+h,\ x+h\in[0,1]$ であり、

$$
\begin{aligned}
\abs{f(x+h)-f(x)}\le{}&\abs{f(x+h)-f(x_k+h)}+\abs{f(x_k+h)-f_k(x_k+h)}\\
&+\abs{f_k(x_k+h)-f_k(x_k)}+\abs{f_k(x_k)-f(x_k)}+\abs{f(x_k)-f(x)}\\
\le{}&\abs{f(x+h)-f(x_k+h)}+2d_\infty(f_k,f)+nh+\abs{f(x_k)-f(x)}
\end{aligned}
$$

である。$k\to\infty$ とすると $f$ の連続性と一様収束から右辺は $nh$ に収束するので、$\abs{f(x+h)-f(x)}\le nh$ となる。従って $f\in E_n$ であり、$E_n$ は閉集合である。

**第 3 段（$E_n$ は内部をもたない）。** $f\in C([0,1])$ と $\eps>0$ に対し、$d_\infty(f,g)<\eps$ で $g\notin E_n$ となる $g$ を作る。

まず $f$ を折れ線で近似する。$f$ は一様連続なので（:ref[cor-top-heine-cantor]）、$\abs{s-t}\le1/N$ なら $\abs{f(s)-f(t)}<\eps/4$ となる $N$ がある。$p$ を、点 $k/N$（$0\le k\le N$）で $f$ と一致し各区間 $[k/N,(k+1)/N]$ で一次関数である折れ線とする。$t\in[k/N,(k+1)/N]$ で $p(t)$ は $f(k/N)$ と $f((k+1)/N)$ の凸結合なので、$\abs{p(t)-f(t)}\le\max\{\abs{f(k/N)-f(t)},\abs{f((k+1)/N)-f(t)}\}<\eps/4$ である。各区間での $p$ の傾きの絶対値の最大値を $L$ とすると、$p$ は $L$-Lipschitz である。

次に鋸歯状関数 $\phi(t)=\min_{k\in\Z}\abs{t-k}$ を用いる。$\phi$ は周期 $1$、値は $[0,1/2]$ にあり、各区間 $[j/2,(j+1)/2]$（$j\in\Z$）上で傾き $\pm1$ の一次関数である。整数 $M$ を

$$
M\ge n,\qquad \frac{\eps M}{6}>n+L
$$

となるように選び、$\sigma(t)=\tfrac\eps2\phi(Mt)$、$g=p+\sigma$ と置く。$0\le\sigma\le\eps/4$ なので $d_\infty(f,g)\le d_\infty(f,p)+d_\infty(p,\sigma+p)<\eps/4+\eps/4<\eps$ である。$\ell=1/(2M)$ と置くと、$\sigma$ は各区間 $I_j=[j\ell,(j+1)\ell]$ 上で傾き $\pm\eps M/2$ の一次関数であり、$\sigma((j+2)\ell)=\sigma(j\ell)$ である。

$x\in[0,1-1/n]$ を任意に取り、$x\in[j\ell,(j+1)\ell)$ となる $j$ を取る。

- $(j+1)\ell-x\ge\ell/2$ のとき、$h=\ell/2$ と置く。$x$ と $x+h$ は同じ区間 $I_j$ に属するので $\abs{\sigma(x+h)-\sigma(x)}=\tfrac{\eps M}2h$ である。
- $(j+1)\ell-x<\ell/2$、すなわち $x-j\ell>\ell/2$ のとき、$h=(j+2)\ell-x$ と置く。$\ell<h<3\ell/2$ であり、$\sigma(x+h)=\sigma(j\ell)$ なので、$I_j$ 上の傾きから $\abs{\sigma(x+h)-\sigma(x)}=\tfrac{\eps M}2(x-j\ell)>\tfrac{\eps M}2\cdot\tfrac\ell2>\tfrac{\eps M}2\cdot\tfrac h3$ である。

いずれの場合も $0<h<3\ell/2=3/(4M)\le1/n$ なので $x+h\le1$ であり、

$$
\abs{\sigma(x+h)-\sigma(x)}\ge\frac{\eps M}{6}\,h
$$

が成り立つ。従って

$$
\abs{g(x+h)-g(x)}\ge\abs{\sigma(x+h)-\sigma(x)}-\abs{p(x+h)-p(x)}\ge\Bigl(\frac{\eps M}6-L\Bigr)h>nh
$$

となる。$x$ は任意なので $g\notin E_n$ である。よって $E_n$ のどの点の $\eps$ 近傍も $E_n$ の外の点を含み、閉集合 $E_n$ の内部は空、すなわち $E_n$ は疎である。

**第 4 段（左側の微分係数と結論）。** $R(f)(t)=f(1-t)$ は $C([0,1])$ の全射等長写像である。$E'_n=R(E_n)$ は閉で内部が空な集合であり、

$$
E'_n=\Bigl\{f:\ \text{ある }x\in[\tfrac1n,1]\text{ があり、すべての }h\in(0,\tfrac1n]\text{ で }\abs{f(x-h)-f(x)}\le nh\Bigr\}
$$

である。第 1 段と同じ議論（$R$ で移したもの）により、$x\in(0,1]$ で有限な左微分係数をもつ $f$ はある $E'_n$ に属する。

$G=C([0,1])\setminus\bigcup_n(E_n\cup E'_n)=\bigcap_n\bigl(C([0,1])\setminus(E_n\cup E'_n)\bigr)$ は稠密な開集合の可算共通部分であり（:ref[lem-top-nowhere-dense]）、Baire の定理により稠密な $G_\delta$ 集合である。$f\in G$ が点 $x\in[0,1)$ で微分可能なら右微分係数をもつので $f\in\bigcup E_n$、$x=1$ で微分可能なら左微分係数をもつので $f\in\bigcup E'_n$ となり、いずれも矛盾する。従って $G$ の元はどの点でも微分可能でない。$C([0,1])\ne\varnothing$ なので $G\ne\varnothing$ である。
:::

:::remark[構成的な例との比較]
この構成により、任意の連続関数を至る所微分不可能な連続関数で一様近似できる。なお定理は「微分不可能な関数の全体が稠密な $G_\delta$ 集合を含む」と述べており、その全体自体が $G_\delta$ 集合であるとは主張していない。
:::

## 縮小写像の原理

方程式 $x=Tx$ を解く最も素朴な方法は、初期値 $x_0$ から $x_{n+1}=Tx_n$ と反復することである。反復が収束するためには、$T$ が二点の距離を一様な割合で縮めればよい。そうすれば反復列は Cauchy 列になり、完備性が極限の存在を保証する。この単純な原理が、常微分方程式の解の存在や逆関数定理（[逆関数定理と陰関数定理](/math/analysis/implicit-function/)）の証明の基礎になる。

:::definition[縮小写像]
距離空間 $(X,d)$ の自己写像 $T:X\to X$ に対し、定数 $0\le q<1$ が存在してすべての $x,y$ で $d(Tx,Ty)\le q\,d(x,y)$ が成り立つとき、$T$ を :term[縮小写像]{en="contraction mapping" yomi="しゅくしょうしゃぞう"}、$q$ をその縮小率という。$Tx=x$ を満たす点 $x$ を $T$ の :term[不動点]{en="fixed point" yomi="ふどうてん"} という。
:::

縮小写像は縮小率を Lipschitz 定数とする Lipschitz 写像（:ref[def-top-metric-continuity] の後の定義）であり、特に一様連続である。

:::theorem[Banach の不動点定理]{#thm-top-contraction}
$(X,d)$ を空でない完備距離空間、$T:X\to X$ を縮小率 $q$ の縮小写像とする。$T$ はただ一つの不動点 $p$ をもつ。任意の初期値 $x_0\in X$ から $x_{n+1}=Tx_n$ で定めた列は $p$ に収束し、次の誤差評価が成り立つ。

$$
d(x_n,p)\le\frac{q^n}{1-q}\,d(x_0,x_1),\qquad d(x_n,p)\le\frac{d(x_n,x_{n+1})}{1-q} \label{eq-top-contraction-estimate}
$$
:::

:::proof
**反復列は Cauchy 列である。** 帰納法で $d(x_{k+1},x_k)\le q^kd(x_1,x_0)$ を示す。$k=0$ では等号であり、$k$ で成り立てば $d(x_{k+2},x_{k+1})=d(Tx_{k+1},Tx_k)\le q\,d(x_{k+1},x_k)\le q^{k+1}d(x_1,x_0)$ である。$m>n$ に対し三角不等式を繰り返して

$$
d(x_n,x_m)\le\sum_{k=n}^{m-1}d(x_k,x_{k+1})\le\sum_{k=n}^{m-1}q^kd(x_0,x_1)\le\frac{q^n}{1-q}d(x_0,x_1)
$$

を得る。$q^n\to0$ なので $(x_n)$ は Cauchy 列であり、完備性から極限 $p$ をもつ。

**極限は不動点である。** $T$ は Lipschitz 連続なので、$Tp=\lim_nTx_n=\lim_nx_{n+1}=p$ である。

**一意性。** $p'$ も不動点なら $d(p,p')=d(Tp,Tp')\le q\,d(p,p')$ で、$1-q>0$ より $d(p,p')=0$ である。

**誤差評価。** 上の不等式で $m\to\infty$ とすれば第一式を得る。また $d(x_n,p)\le d(x_n,x_{n+1})+d(x_{n+1},p)=d(x_n,x_{n+1})+d(Tx_n,Tp)\le d(x_n,x_{n+1})+q\,d(x_n,p)$ を整理すれば第二式を得る。$q=0$ のときは $T$ が定数写像で、$x_1=p$ である（$q^0=1$ と約束する）。
:::

第一式は反復を始める前に必要な回数を見積もる事前評価、第二式は計算した二項の差から誤差を見積もる事後評価である。

:::example[反復に必要な回数]
$T(x)=1+x/3$ は $\R$ 上の縮小率 $q=1/3$ の縮小写像で、不動点は $p=3/2$ である。$x_0=0$ なら $d(x_0,x_1)=1$ なので、事前評価は $d(x_n,p)\le\tfrac32\cdot3^{-n}$ である。誤差を $10^{-3}$ 未満にするには、$n=6$ で $1/486>10^{-3}$、$n=7$ で $1/1458<10^{-3}$ なので $n=7$ で足りる。実際、帰納法で $x_n=\tfrac32(1-3^{-n})$ が示されるので $p-x_n=\tfrac32\cdot3^{-n}$ であり、この例では事前評価が等号で成り立つ。
:::

:::example[方程式 x = cos x]
$T(x)=\cos x$ を $X=[0,1]$ 上で考える。$\cos$ は $[0,1]$ で減少するので $T([0,1])=[\cos1,1]\subset[0,1]$ であり（$\cos1\approx0.5403$）、$T$ は $[0,1]$ の自己写像である。平均値の定理（[1 変数関数の微分](/math/analysis/differentiation/)）から、ある $\xi\in[0,1]$ で $\abs{\cos x-\cos y}=\abs{\sin\xi}\abs{x-y}\le\sin1\cdot\abs{x-y}$ なので、$T$ は縮小率 $q=\sin1\approx0.8415$ の縮小写像である。$[0,1]$ は完備なので、$x=\cos x$ は $[0,1]$ にただ一つの解 $p\approx0.739085$ をもつ。

$x_0=1$ なら $d(x_0,x_1)=1-\cos1\approx0.4597$ であり、事前評価から誤差を $10^{-6}$ 未満にするには $n\ge87$ で足りる。実際の反復では $n=32$ で誤差が $10^{-6}$ 未満になる。収束の速さは $p$ の近くでの $\abs{T'(p)}=\sin p\approx0.674$ で決まり、全区間での縮小率 $\sin1$ より小さいからである。事前評価は安全側の見積もりであり、事後評価を併用すると実際の誤差に近い判定ができる。
:::

:::example[仮定を一つでも外すと]{#exm-top-contraction-hypotheses}
1. **完備性**：$(0,1)$ 上の $T(x)=x/2$ は縮小写像だが、不動点 $0$ は空間の外にある。
2. **一様な縮小率**：完備な $[1,\infty)$ 上の $T(x)=x+1/x$ は自己写像で（$x+1/x\ge2$）、$x\ne y$ なら

$$
\abs{T(x)-T(y)}=\abs{x-y}\Bigl\lvert1-\frac1{xy}\Bigr\rvert<\abs{x-y}
$$

と距離を真に縮めるが、不動点をもたない（$x=x+1/x$ は解をもたない）。$x,y$ がともに大きいと比が $1$ に近づき、一つの $q<1$ を選べない。
3. **自己写像であること**：$[0,1]$ 上の $T(x)=(x+3)/2$ は縮小率 $1/2$ で距離を縮めるが、$T(1)=2\notin[0,1]$ であり、不動点 $3$ は $[0,1]$ の外にある。
:::

:::corollary[反復が縮小写像になる場合]{#cor-top-contraction-power}
$X$ を空でない完備距離空間、$T:X\to X$ を写像とし、ある正整数 $m$ について $T^m$ が縮小写像であるとする。このとき $T$ はただ一つの不動点をもつ。
:::

:::proof
$T^m$ の一意な不動点を $p$ とする。$T^m(Tp)=T(T^mp)=Tp$ なので $Tp$ も $T^m$ の不動点であり、一意性から $Tp=p$ である。また $T$ の不動点は $T^m$ の不動点なので $p$ 以外にない。$T$ 自身の連続性は使っていない。
:::

:::proposition[不動点のパラメータへの連続依存性]{#prop-top-contraction-parameter}
$\Lambda$ を位相空間、$X$ を空でない完備距離空間とし、$T:\Lambda\times X\to X$ が次を満たすとする。

1. ある $0\le q<1$ について、すべての $\lambda$ と $x,y$ で $d(T(\lambda,x),T(\lambda,y))\le q\,d(x,y)$。
2. 各 $x$ について $\lambda\mapsto T(\lambda,x)$ は連続。

このとき $T(\lambda,\cdot)$ の不動点 $p(\lambda)$ は $\lambda$ について連続である。
:::

:::proof
$\lambda,\mu\in\Lambda$ に対し

$$
\begin{aligned}
d(p(\lambda),p(\mu))&=d\bigl(T(\lambda,p(\lambda)),T(\mu,p(\mu))\bigr)\\
&\le d\bigl(T(\lambda,p(\lambda)),T(\lambda,p(\mu))\bigr)+d\bigl(T(\lambda,p(\mu)),T(\mu,p(\mu))\bigr)\\
&\le q\,d(p(\lambda),p(\mu))+d\bigl(T(\lambda,p(\mu)),T(\mu,p(\mu))\bigr)
\end{aligned}
$$

なので、$d(p(\lambda),p(\mu))\le d\bigl(T(\lambda,p(\mu)),T(\mu,p(\mu))\bigr)/(1-q)$ である。$\mu$ を固定して $\lambda\to\mu$ とすると、仮定 2 を点 $x=p(\mu)$ に使って右辺は $0$ に近づく。
:::

:::example[距離を真に縮めるコンパクト空間の自己写像]
空でないコンパクト距離空間 $X$ の写像 $T$ が、$x\ne y$ なら $d(Tx,Ty)<d(x,y)$ を満たすとする。$T$ は $1$-Lipschitz なので連続であり、$g(x)=d(x,Tx)$ も連続である。コンパクト性から $g$ は最小点 $p$ をもつ（:ref[thm-top-compact-image]）。$p\ne Tp$ なら $g(Tp)=d(Tp,T^2p)<d(p,Tp)=g(p)$ となって最小性に反するので、$p$ は不動点である。不動点が二つあれば、その二点で厳密な不等式に反するので一意である。:ref[exm-top-contraction-hypotheses] の 2 の $[1,\infty)$ では $g(x)=1/x$ が下限 $0$ を達成しないので、この議論は使えない。
:::

### 常微分方程式の初期値問題

縮小写像の原理の代表的な応用として、常微分方程式の初期値問題

$$
y'(t)=f(t,y(t)),\qquad y(t_0)=y_0 \label{eq-top-ivp}
$$

の解の存在と一意性を証明する。本節では $\R^m$ に最大値ノルム $\norm y=\max_{1\le i\le m}\abs{y_i}$ を入れる。$\R^m$ 値関数の微分と積分は成分ごとに行い、連続関数 $g:[\alpha,\beta]\to\R^m$ について各成分で $\abs{\int g_i}\le\int\abs{g_i}\le\int\norm g$ なので

$$
\Bigl\lVert\int_\alpha^\beta g(s)\,\dd s\Bigr\rVert\le\int_\alpha^\beta\norm{g(s)}\,\dd s
$$

が成り立つ。

:::lemma[積分方程式への書き換え]{#lem-top-ode-integral}
$I$ を $t_0$ を含む区間、$f$ を $\R\times\R^m$ の部分集合 $\Omega$ 上の連続写像とし、$u:I\to\R^m$ を連続で、すべての $t\in I$ で $(t,u(t))\in\Omega$ とする。$u$ が $C^1$ 級で :eqref[eq-top-ivp] を満たすことと、

$$
u(t)=y_0+\int_{t_0}^tf(s,u(s))\,\dd s\qquad(t\in I) \label{eq-top-picard-integral}
$$

が成り立つことは同値である。
:::

:::proof
$s\mapsto f(s,u(s))$ は連続写像の合成なので連続である。$u$ が $C^1$ 級の解なら、微積分学の基本定理（[Riemann 積分](/math/analysis/riemann-integral/)）により $u(t)-u(t_0)=\int_{t_0}^tu'(s)\,\dd s=\int_{t_0}^tf(s,u(s))\,\dd s$ である。逆に :eqref[eq-top-picard-integral] が成り立てば、右辺は連続関数の不定積分なので $C^1$ 級で導関数は $f(t,u(t))$ であり、$u(t_0)=y_0$ である。
:::

:::theorem[Picard–Lindelöf の定理]{#thm-top-picard-lindelof}
$a,b>0$ とし、$\Omega=\{(t,y):\abs{t-t_0}\le a,\ \norm{y-y_0}\le b\}$ 上の連続写像 $f:\Omega\to\R^m$ が $y$ について Lipschitz 条件

$$
\norm{f(t,y)-f(t,z)}\le L\norm{y-z}\qquad((t,y),(t,z)\in\Omega)
$$

を満たすとする。$M=\max_\Omega\norm f$ と置き、$h=\min\{a,b/M\}$（$M=0$ なら $h=a$）、$I=[t_0-h,t_0+h]$ とする。

1. $I$ 上で定義され、$\norm{y(t)-y_0}\le b$ を満たす :eqref[eq-top-ivp] の $C^1$ 級の解 $y$ がただ一つ存在する。
2. $y_0(t)\equiv y_0$、$y_{k+1}(t)=y_0+\int_{t_0}^tf(s,y_k(s))\,\dd s$ で定める Picard 反復 $(y_k)$ は $I$ 上で $y$ に一様収束する。
3. $f$ が $[t_0-a,t_0+a]\times\R^m$ 上で連続かつ $y$ について大域的に Lipschitz なら、解は区間 $[t_0-a,t_0+a]$ 全体でただ一つ存在する。
:::

:::proof
$\Omega$ はコンパクトなので（:ref[thm-top-heine-borel]）連続関数 $\norm f$ は最大値 $M$ をもつ。$L=0$ の場合も以下の議論はそのまま通る（そのとき $T$ は値の差を $0$ に写す）ので、$L>0$ としてよい。

**空間。** $X=\{u\in C(I,\R^m):\ \norm{u(t)-y_0}\le b\ (t\in I)\}$ と置く。:ref[prop-top-bounded-maps-complete] により $C(I,\R^m)$ は一様距離 $d_\infty$ で完備であり、$X$ はその閉集合（各点での閉条件の共通部分）なので完備である。$X$ 上に重み付きの距離

$$
d_L(u,v)=\max_{t\in I}e^{-2L\abs{t-t_0}}\norm{u(t)-v(t)}
$$

を入れる。$e^{-2Lh}d_\infty\le d_L\le d_\infty$ なので $d_L$ は $d_\infty$ と双 Lipschitz 同値であり、:ref[prop-top-uniform-equivalence-complete] により $(X,d_L)$ も完備である。

**写像。** $u\in X$ に対し $Tu(t)=y_0+\int_{t_0}^tf(s,u(s))\,\dd s$ と置く。被積分関数は連続なので $Tu$ は連続であり、

$$
\norm{Tu(t)-y_0}\le\Bigl\lvert\int_{t_0}^t\norm{f(s,u(s))}\,\dd s\Bigr\rvert\le M\abs{t-t_0}\le Mh\le b
$$

なので $Tu\in X$ である。

**縮小性。** $u,v\in X$、$t\ge t_0$ とすると、$\norm{u(s)-v(s)}\le e^{2L(s-t_0)}d_L(u,v)$ より

$$
\norm{Tu(t)-Tv(t)}\le\int_{t_0}^tL\norm{u(s)-v(s)}\,\dd s\le Ld_L(u,v)\int_{t_0}^te^{2L(s-t_0)}\,\dd s=\frac{e^{2L(t-t_0)}-1}2\,d_L(u,v)
$$

である。両辺に $e^{-2L(t-t_0)}$ を掛けると $e^{-2L\abs{t-t_0}}\norm{Tu(t)-Tv(t)}\le\tfrac12d_L(u,v)$ を得る。$t<t_0$ の場合も積分の向きを変えて同じ評価を得る。従って $d_L(Tu,Tv)\le\tfrac12d_L(u,v)$ であり、$T$ は縮小率 $1/2$ の縮小写像である。

**結論。** Banach の不動点定理により、$T$ は $X$ にただ一つの不動点 $y$ をもつ。:ref[lem-top-ode-integral] により、$X$ の元が $T$ の不動点であることと、$\norm{y(t)-y_0}\le b$ を満たす $I$ 上の $C^1$ 級の解であることは同値なので 1 が従う。Picard 反復は $y_0\in X$ からの反復列 $T^ky_0$ であり、$d_L$ で $y$ に収束する。$d_\infty\le e^{2Lh}d_L$ なので一様収束でもあり、2 が従う。

3 では $X=C([t_0-a,t_0+a],\R^m)$ 全体を取り、同じ距離 $d_L$（$h$ を $a$ に置き換える）を用いる。値の範囲の制限がないので $T$ は $X$ の自己写像であり、上と同じ計算で縮小率 $1/2$ の縮小写像となる。
:::

:::remark[区間の取り方]
重みのない一様距離 $d_\infty$ を使うと、評価は $d_\infty(Tu,Tv)\le hL\,d_\infty(u,v)$ となり、縮小性のために $hL<1$ という余分な条件が必要になる。重み $e^{-2L\abs{t-t_0}}$ は、$t_0$ から離れるほど誤差の積み重なりを許容する距離を与え、この制限を取り除く。また解の一意性は、$t_0$ を含む部分区間 $J\subset I$ 上で定義された解どうしでも成り立つ（任意の $t\in J$ に対し、$t_0$ と $t$ を両端とする閉区間に二解を制限して同じ縮小評価を適用する）。
:::

:::example[y′ = y の Picard 反復]
$0\le t\le1/2$ で $y'=y$、$y(0)=1$ を考える。$f(t,y)=y$ は大域的に Lipschitz（$L=1$）なので、定理の 3 により解はただ一つである。この例では重みなしの一様距離でも、

$$
\sup_{0\le t\le1/2}\abs{Tu(t)-Tv(t)}\le\sup_{0\le t\le1/2}\int_0^t\abs{u(s)-v(s)}\,\dd s\le\frac12d_\infty(u,v)
$$

より $T$ は縮小率 $1/2$ の縮小写像である。$y_0=1$ から始めると $y_1=1+t$、$y_2=1+t+t^2/2$ であり、帰納的に

$$
y_k(t)=\sum_{j=0}^k\frac{t^j}{j!}
$$

となる。実際 $1+\int_0^t\sum_{j\le k}s^j/j!\,\dd s=1+\sum_{j\le k}t^{j+1}/(j+1)!$ である。$d_\infty(y_0,y_1)=\sup t=1/2$ なので事前評価は $d_\infty(y_k,y)\le(1/2)^k\cdot2\cdot\tfrac12=2^{-k}$ となる。極限は指数関数 $e^t$ であり、実際の誤差 $e^{1/2}-\sum_{j\le k}2^{-j}/j!$ は例えば $k=4$ で約 $2.8\times10^{-4}$ と、評価 $2^{-4}=0.0625$ より十分小さい。この議論は解の式を先に知らなくても、存在と一意性と近似の方法を同時に与える。
:::

:::example[解は局所的にしか存在しないことがある]
$y'=y^2$、$y(0)=1$ の解は $y(t)=1/(1-t)$ であり、$t\to1$ で発散するので $t<1$ でしか存在しない。定理を $a=1$、$b>0$ で適用すると、$\Omega=[-1,1]\times[1-b,1+b]$ 上で $M=(1+b)^2$、$L=2(1+b)$ であり、$h=\min\{1,\ b/(1+b)^2\}$ である。$b/(1+b)^2$ は $b=1$ で最大値 $1/4$ を取るので、定理が保証する区間は高々 $[-1/4,1/4]$ である。定理は解の存在を局所的にしか主張せず、その区間も実際の存在区間よりずっと小さいことがある。解を延長する方法は [解の存在と一意性](/math/ode/existence-uniqueness/) で扱う。
:::

:::example[Lipschitz 条件がないと一意性が崩れる]
$y'=3y^{2/3}$、$y(0)=0$（$y^{2/3}=(\sqrt[3]y)^2$）には、$y\equiv0$ と $y=t^3$ という二つの解がある。実際 $(t^3)'=3t^2=3(t^3)^{2/3}$ である。$f(y)=3y^{2/3}$ は原点の近くで Lipschitz 条件を満たさない。$\abs{f(y)-f(0)}/\abs{y-0}=3\abs y^{-1/3}\to\infty$（$y\to0$）だからである。この方程式でも解の存在は保証される（:ref[thm-top-peano]）が、一意性には Lipschitz 条件のような追加の仮定が必要である。
:::

## Arzelà–Ascoli の定理

関数列の各点での値を個別に抑えても、一様収束する部分列が取れるとは限らない。例えば細い山の位置が動いていく関数列は、各点で $0$ に収束しても一様には収束しない。Arzelà–Ascoli の定理は、近い二点での値の差を族のすべての関数について同時に制御すれば（同程度連続性）、各点での値の集合の閉包がコンパクトなら、一様収束部分列が取れることを主張する。これは有限次元の Bolzano–Weierstrass の定理を関数空間に拡張したものである。

:::definition[同程度連続性と相対コンパクト性]{#def-top-equicontinuity}
$K$ を位相空間、$(Y,\rho)$ を距離空間とする。$\mathcal H\subset C(K,Y)$ が点 $x\in K$ で :term[同程度連続]{en="equicontinuous" yomi="どうていどれんぞく"} であるとは、各 $\eps>0$ に対し $x$ の開近傍 $U$ が存在し、

$$
\rho\bigl(h(x),h(z)\bigr)<\eps\qquad(h\in\mathcal H,\ z\in U)
$$

が成り立つことをいう。すべての点で同程度連続であるとき、単に同程度連続という。位相空間の部分集合で、その閉包がコンパクトであるものを :term[相対コンパクト]{en="relatively compact" yomi="そうたいこんぱくと"} という。
:::

同程度連続性では、近傍 $U$ が点 $x$ と $\eps$ だけで決まり、族の関数 $h$ によらない。$K$ がコンパクトなら $C(K,Y)$ には一様距離 $d_\infty$ が入る（:ref[prop-top-bounded-maps-complete]）。以下、$C(K,Y)$ の位相はこの一様距離によるものとする。

:::example[同程度連続な族とそうでない族]{#exm-top-equicontinuous}
1. 有限個の連続写像の族は同程度連続である。各関数について近傍を取り、その有限個の共通部分を取ればよい。
2. 距離空間 $K$ 上の $L$-Lipschitz 写像の族は同程度連続である。$\delta=\eps/(L+1)$ とすれば、$d(x,z)<\delta$ で $\rho(h(x),h(z))\le Ld(x,z)<\eps$ がすべての $h$ について成り立つ。
3. $[0,1]$ 上の $f_n(t)=t^n$ の族は $t=1$ で同程度連続でない。$1$ のどの近傍もある $s<1$ を含み、$s^n\to0$ なので十分大きい $n$ で $\abs{f_n(1)-f_n(s)}=1-s^n>1/2$ となるからである。
:::

:::theorem[Arzelà–Ascoli の定理]{#thm-top-arzela-ascoli}
$K$ をコンパクト空間、$(Y,\rho)$ を距離空間とする。$\mathcal H\subset C(K,Y)$ が $d_\infty$ について相対コンパクトであるための必要十分条件は、次の二条件が成り立つことである。

1. $\mathcal H$ は同程度連続である。
2. 各 $x\in K$ で、$\mathcal H(x)=\{h(x):h\in\mathcal H\}$ は $Y$ で相対コンパクトである。

$Y=\R^m$ のとき、条件 2 は各点での有界性 $\sup_{h\in\mathcal H}\norm{h(x)}<\infty$ と同値である。一般の距離空間 $Y$ でも、条件 1、2 がそろえば $\mathcal H$ の任意の列は一様収束する部分列をもつ。
:::

:::proof
各 $x\in K$ について評価写像 $e_x:C(K,Y)\to Y$、$e_x(f)=f(x)$ を考える。$\rho(f(x),g(x))\le d_\infty(f,g)$ なので $e_x$ は $1$-Lipschitz で連続である。$\mathcal H$ の $C(K,Y)$ での閉包を $\mathcal C$ と書く。

**必要性：各点での相対コンパクト性。** $\mathcal C$ がコンパクトなら、$e_x(\mathcal C)$ はコンパクトで（:ref[thm-top-compact-image]）、距離空間 $Y$ の中で閉である（:ref[prop-top-compact-closed]）。$\mathcal H(x)\subset e_x(\mathcal C)$ なので $\overline{\mathcal H(x)}\subset e_x(\mathcal C)$ であり、$\overline{\mathcal H(x)}$ はコンパクト集合の閉部分集合としてコンパクトである。

**必要性：同程度連続性。** $x\in K$ と $\eps>0$ を取る。$\mathcal H=\varnothing$ なら主張は自明なので空でないとする。コンパクト距離空間 $\mathcal C$ は全有界なので（:ref[thm-top-metric-compactness]）、$g_1,\dots,g_r\in\mathcal C$ で $\mathcal C\subset\bigcup_jB(g_j,\eps/3)$ となるものがある。各 $g_j$ は $x$ で連続なので、$z\in U_j$ なら $\rho(g_j(z),g_j(x))<\eps/3$ となる開近傍 $U_j$ がある。$U=U_1\cap\cdots\cap U_r$ と置く。$h\in\mathcal H$ に対し $d_\infty(h,g_j)<\eps/3$ となる $j$ を選べば、$z\in U$ のとき

$$
\rho(h(z),h(x))\le\rho(h(z),g_j(z))+\rho(g_j(z),g_j(x))+\rho(g_j(x),h(x))<\eps
$$

である。

**十分性：全有界性。** $\mathcal H\ne\varnothing$ とし、$\eps>0$ を固定する。同程度連続性により、各 $x\in K$ に開近傍 $U_x$ を、$z\in U_x$ とすべての $h\in\mathcal H$ で $\rho(h(z),h(x))<\eps/5$ となるように取る。$K$ はコンパクトなので有限個 $U_{x_1},\dots,U_{x_s}$ が $K$ を覆う。各 $i$ について $\overline{\mathcal H(x_i)}$ はコンパクトなので全有界であり、半径 $\eps/5$ の有限個の開球 $B_{i,1},\dots,B_{i,r_i}$ で覆える。各 $h\in\mathcal H$ に、$h(x_i)\in B_{i,\nu_i}$ となる番号の組 $\nu(h)=(\nu_1,\dots,\nu_s)$ を一つ対応させる。組の種類は有限個である。

$\nu(h)=\nu(g)$ なら、$h(x_i),g(x_i)$ は同じ半径 $\eps/5$ の球に属するので $\rho(h(x_i),g(x_i))<2\eps/5$ である。任意の $z\in K$ は、ある $U_{x_i}$ に属するので

$$
\rho(h(z),g(z))\le\rho(h(z),h(x_i))+\rho(h(x_i),g(x_i))+\rho(g(x_i),g(z))<\frac\eps5+\frac{2\eps}5+\frac\eps5=\frac{4\eps}5
$$

となり、$d_\infty(h,g)\le4\eps/5<\eps$ である。実際に現れる各組から代表 $h$ を一つずつ選べば、有限個の $h$ を中心とする半径 $\eps$ の球で $\mathcal H$ が覆われる。$\eps$ は任意なので $\mathcal H$ は全有界である。

さらに閉包 $\mathcal C$ も全有界である。半径 $\eps/2$ の有限個の球 $B(h_j,\eps/2)$（$h_j\in\mathcal H$）が $\mathcal H$ を覆えば、閉球の和 $\bigcup_jD(h_j,\eps/2)$ は $\mathcal H$ を含む閉集合なので $\mathcal C$ を含み、$\mathcal C\subset\bigcup_jB(h_j,\eps)$ となるからである。

**十分性：閉包の完備性。** $(g_n)$ を $\mathcal C$ の Cauchy 列とする。各 $x$ で $\rho(g_n(x),g_m(x))\le d_\infty(g_n,g_m)$ なので $(g_n(x))$ は $Y$ の Cauchy 列である。$e_x$ の連続性から $e_x(\mathcal C)\subset\overline{e_x(\mathcal H)}=\overline{\mathcal H(x)}$ であり、$\overline{\mathcal H(x)}$ はコンパクト距離空間なので完備である（:ref[thm-top-metric-compactness]）。従って $g_n(x)$ は $\overline{\mathcal H(x)}$ の点 $g(x)$ に収束する。$(g_n)$ は一様 Cauchy 条件を満たすので、:ref[lem-top-uniform-cauchy] により $g_n\to g$ は一様収束で、$g$ は連続である。$K$ はコンパクトなので $g\in C(K,Y)$ で $d_\infty(g_n,g)\to0$ であり、$\mathcal C$ は閉なので $g\in\mathcal C$ である。

**結論。** $\mathcal C$ は完備かつ全有界なのでコンパクトである（:ref[thm-top-metric-compactness]）。コンパクト距離空間は点列コンパクトなので、$\mathcal H$ の任意の列は $\mathcal C$ の中で収束する部分列、すなわち一様収束する部分列をもつ。$Y=\R^m$ のとき、部分集合が相対コンパクトであることは有界であることと同値である（:ref[thm-top-heine-borel]）。
:::

定理の二条件はそれぞれ独立に必要である。

:::example[二条件の役割]
1. $[0,1]$ 上で $f(0)=0$ かつ $\abs{f(s)-f(t)}\le2\abs{s-t}$ を満たす関数全体 $\mathcal L$ を考える。$\abs{f(t)}=\abs{f(t)-f(0)}\le2$ なので各点で有界であり、:ref[exm-top-equicontinuous] の 2 により同程度連続である。従って相対コンパクトである。さらに一様極限では $f(0)=0$ と Lipschitz 不等式が保たれるので $\mathcal L$ は閉であり、$\mathcal L$ 自身がコンパクトである。
2. $f_n(t)=t^n$（$0\le t\le1$）は各点で有界だが、$t=1$ で同程度連続でない（:ref[exm-top-equicontinuous] の 3）。どの部分列も各点では不連続な関数（$t<1$ で $0$、$t=1$ で $1$）に収束するので、一様収束すれば極限は連続でなければならないことに反し、一様収束する部分列はない。
3. 定数関数の列 $g_n(t)=n$ は同程度連続だが、一点での値すら有界でない。$d_\infty(g_m,g_n)=\abs{m-n}\ge1$ なので一様収束する部分列はない。
:::

コンパクト距離空間上では、同程度連続性の近傍を共通の半径にそろえられる。

:::proposition[コンパクト距離空間での一様同程度連続性]{#prop-top-uniform-equicontinuity}
$K$ がコンパクト距離空間なら、$\mathcal H\subset C(K,Y)$ の同程度連続性は次の条件と同値である。各 $\eps>0$ に対し一つの $\delta>0$ があり、

$$
d(x,y)<\delta\quad\Longrightarrow\quad\rho(h(x),h(y))<\eps
$$

がすべての $h\in\mathcal H$ とすべての $x,y\in K$ で成り立つ。
:::

:::proof
同程度連続であるとし、$\eps>0$ を取る。各 $a\in K$ に対し、すべての $h$ と $z\in U_a$ で $\rho(h(a),h(z))<\eps/2$ となる開近傍 $U_a$ を取る。開被覆 $(U_a)_{a\in K}$ の Lebesgue 数（:ref[lem-top-lebesgue-number]）を $\delta>0$ とする。$d(x,y)<\delta$ なら二点集合 $\{x,y\}$ の直径は $\delta$ 未満なので、ある $U_a$ に含まれ、

$$
\rho(h(x),h(y))\le\rho(h(x),h(a))+\rho(h(a),h(y))<\eps
$$

となる。逆に共通の $\delta$ があれば、各点 $x$ の近傍として $B(x,\delta)$ を選べば同程度連続性が従う。
:::

:::example[定義域がコンパクトでないとき]
$\R$ 上で $h_n(x)=\max(0,1-\abs{x-n})$（$n\ge1$）とする。各関数は Lipschitz 定数 $1$ をもち、値は $[0,1]$ に入るので、族は同程度連続で各点で有界である。しかし $m\ne n$ なら $h_n(n)=1$、$h_m(n)=0$ なので $\sup_x\abs{h_n(x)-h_m(x)}=1$ であり、$\R$ 全体で一様収束する部分列はない。一方、任意の有界閉区間 $[-R,R]$ 上では、$n>R+1$ で $h_n$ は恒等的に $0$ なので、一様に $0$ に収束する。移動していく山は、コンパクトな定義域の上でしか捕まえられない。Arzelà–Ascoli の証明で定義域を有限個の近傍で覆ったことが、この違いに対応している。
:::

距離化可能なコンパクト空間上では、Arzelà–Ascoli の結論を対角線論法で直接示すこともできる。この方法は収束部分列の構成を具体的に与える。

:::proposition[対角線論法による部分列の選択]{#prop-top-ascoli-diagonal}
$K$ をコンパクト距離空間、$(Y,\rho)$ を距離空間とし、関数列 $(f_n)$ が同程度連続で、各点 $x$ で $\{f_n(x)\}$ が相対コンパクトであるとする。可算稠密集合上で順に収束部分列を選ぶ方法で、一様収束する部分列を構成できる。
:::

:::proof
$K=\varnothing$ なら明らかなので $K\ne\varnothing$ とする。$K$ は全有界なので、各 $k$ について有限 $1/k$ ネットを取り、それらを合わせれば可算稠密集合 $\{a_1,a_2,\dots\}$ が得られる（有限なら同じ点を繰り返して並べる）。

**対角線論法。** $\overline{\{f_n(a_1)\}}$ はコンパクト距離空間なので点列コンパクトであり、$(f_n(a_1))$ が収束する部分列 $(f_{1,n})_n$ を選べる。その部分列から $a_2$ での値も収束する部分列 $(f_{2,n})_n$ を選ぶ。これを続けて、第 $i$ 段の列 $(f_{i,n})_n$ が第 $i-1$ 段の列の部分列で、$a_i$ での値が収束するようにする。$g_k=f_{k,k}$ と置くと、$(g_k)$ は元の列の部分列であり（元の添字は狭義増加である）、各 $i$ について $(g_k)_{k\ge i}$ は第 $i$ 段の列の部分列なので、$g_k(a_i)$ は収束する。

**一様 Cauchy 性。** $\eps>0$ を取る。同程度連続性から各 $x$ に開近傍 $U_x$ を、$z\in U_x$ とすべての $n$ で $\rho(f_n(z),f_n(x))<\eps/6$ となるように取る。有限個 $U_{x_1},\dots,U_{x_s}$ で $K$ を覆い、各 $U_{x_j}$ から稠密集合の点 $b_j$ を一つ選ぶ。$z,b_j\in U_{x_j}$ なら $\rho(f_n(z),f_n(b_j))<\eps/3$ である。有限個の $b_j$ で $(g_k(b_j))_k$ は収束するので Cauchy 列であり、$k,l\ge N$ ならすべての $j$ で $\rho(g_k(b_j),g_l(b_j))<\eps/6$ となる $N$ がある。任意の $z$ について、$z\in U_{x_j}$ となる $j$ を取れば、$k,l\ge N$ で

$$
\rho(g_k(z),g_l(z))\le\rho(g_k(z),g_k(b_j))+\rho(g_k(b_j),g_l(b_j))+\rho(g_l(b_j),g_l(z))<\frac\eps3+\frac\eps6+\frac\eps3=\frac{5\eps}6
$$

である。従って $(g_k)$ は一様 Cauchy 条件を満たす。

**収束。** 各点 $x$ で $(g_k(x))$ はコンパクト集合 $\overline{\{f_n(x)\}}$ の Cauchy 列なので収束する。:ref[lem-top-uniform-cauchy] により $(g_k)$ は連続関数に一様収束する。
:::

### 応用：導関数の評価と Peano の存在定理

:::corollary[導関数が一様に有界な関数列]{#cor-top-ascoli-derivative}
$[a,b]$ 上の $C^1$ 級関数列 $f_n$ が、ある定数 $C$ について $\abs{f_n(a)}\le C$ と $\abs{f'_n(t)}\le C$（すべての $n$ と $t$）を満たすとする。このとき $(f_n)$ は一様収束する部分列をもち、極限は $C$-Lipschitz である。
:::

:::proof
平均値の定理（[1 変数関数の微分](/math/analysis/differentiation/)）により $\abs{f_n(s)-f_n(t)}\le C\abs{s-t}$ なので、族は $C$-Lipschitz で同程度連続である。また $\abs{f_n(t)}\le\abs{f_n(a)}+C(t-a)\le C(1+b-a)$ なので各点で有界である。Arzelà–Ascoli の定理により一様収束する部分列 $f_{n_k}\to f$ がある。$\abs{f_{n_k}(s)-f_{n_k}(t)}\le C\abs{s-t}$ で $k\to\infty$ とすれば、$f$ も $C$-Lipschitz である。
:::

Picard–Lindelöf の定理では Lipschitz 条件を用いた。連続性だけを仮定すると一意性は失われるが、解の存在は Arzelà–Ascoli の定理から従う。

:::theorem[Peano の存在定理]{#thm-top-peano}
$a,b>0$ とする。$\Omega=\{(t,y):t_0\le t\le t_0+a,\ \norm{y-y_0}\le b\}$ 上の連続写像 $f:\Omega\to\R^m$ を考え、$M=\max_\Omega\norm f$、$h=\min\{a,b/M\}$（$M=0$ なら $h=a$）と置く。このとき $[t_0,t_0+h]$ 上で定義され、$\norm{y(t)-y_0}\le b$ を満たす :eqref[eq-top-ivp] の $C^1$ 級の解が少なくとも一つ存在する。
:::

:::proof
**Euler の折れ線。** $n\ge1$ に対し $t_k=t_0+kh/n$（$0\le k\le n$）と置き、$\varphi_n(t_0)=y_0$、

$$
\varphi_n(t)=\varphi_n(t_k)+(t-t_k)\,f\bigl(t_k,\varphi_n(t_k)\bigr)\qquad(t_k\le t\le t_{k+1})
$$

で $\varphi_n$ を順に定める。これが定義できること、すなわち $(t_k,\varphi_n(t_k))\in\Omega$ であることを帰納法で示す。$[t_0,t_k]$ 上で $\norm{\varphi_n(t)-y_0}\le M(t-t_0)$ が成り立つとすると、$\norm{\varphi_n(t_k)-y_0}\le Mh\le b$ なので $f(t_k,\varphi_n(t_k))$ が定義され、その最大値ノルムは $M$ 以下である。従って $[t_k,t_{k+1}]$ 上で $\norm{\varphi_n(t)-y_0}\le M(t_k-t_0)+M(t-t_k)=M(t-t_0)$ となる。こうして $\varphi_n$ は $[t_0,t_0+h]$ 上の連続な折れ線で、各区間での傾きのノルムが $M$ 以下なので $M$-Lipschitz である。

**部分列の収束。** 族 $\{\varphi_n\}$ は $M$-Lipschitz なので同程度連続であり、値は有界集合 $\{\norm{y-y_0}\le b\}$ に入る。Arzelà–Ascoli の定理により、一様収束する部分列 $\varphi_{n_j}\to\varphi$ がある。$\varphi$ は連続で、$\norm{\varphi(t)-y_0}\le b$ を満たす。

**積分方程式。** 階段関数 $g_n(s)=f(t_k,\varphi_n(t_k))$（$t_k\le s<t_{k+1}$）、$g_n(t_0+h)=f(t_{n-1},\varphi_n(t_{n-1}))$ を考えると、各区間で積分すれば

$$
\varphi_n(t)=y_0+\int_{t_0}^tg_n(s)\,\dd s
$$

である。$\Omega$ はコンパクトなので $f$ は一様連続であり（:ref[cor-top-heine-cantor]）、$\eps>0$ に対し $\abs{t-t'}<\delta$、$\norm{y-y'}<\delta$ なら $\norm{f(t,y)-f(t',y')}<\eps$ となる $\delta>0$ がある。$s\in[t_k,t_{k+1})$ なら $\abs{s-t_k}\le h/n$ であり、

$$
\norm{\varphi_n(t_k)-\varphi(s)}\le\norm{\varphi_n(t_k)-\varphi_n(s)}+\norm{\varphi_n(s)-\varphi(s)}\le\frac{Mh}n+d_\infty(\varphi_n,\varphi)
$$

である。$n=n_j$ が十分大きければ $h/n<\delta$ かつ右辺が $\delta$ 未満になるので、すべての $s$ で $\norm{g_n(s)-f(s,\varphi(s))}<\eps$ となる（$s=t_0+h$ でも同様である）。従って

$$
\Bigl\lVert\int_{t_0}^tg_{n_j}(s)\,\dd s-\int_{t_0}^tf(s,\varphi(s))\,\dd s\Bigr\rVert\le\eps h
$$

であり、$j\to\infty$ として $\varphi(t)=y_0+\int_{t_0}^tf(s,\varphi(s))\,\dd s$ を得る。:ref[lem-top-ode-integral] により $\varphi$ は $C^1$ 級の解である。
:::

:::remark[Peano の定理について]
$y'=3y^{2/3}$、$y(0)=0$ には零解と $y(t)=t^3$ がある。Euler の折れ線は各節点で傾きが零なので、どの刻み幅でも零解そのものになる。従ってこの構成で全解が得られるわけではない。また、$f$ が $|t-t_0|\le a$ の両側で同じ条件を満たす場合には、$u'(s)=-f(2t_0-s,u(s))$ を右向きに解けば左側の解も得られる。両解の $t_0$ での導関数はともに $f(t_0,y_0)$ なので、貼り合わせた解は $C^1$ 級である。
:::

:::example[正則性は極限にどこまで残るか]
$[-1,1]$ 上の $f_n(t)=\sqrt{t^2+1/n}$ はどれも微分可能で、$f_n'(t)=t/\sqrt{t^2+1/n}$ より $\abs{f'_n}\le1$ である。$f_n(0)=1/\sqrt n\le1$ なので :ref[cor-top-ascoli-derivative] の仮定を満たす。実際

$$
0\le f_n(t)-\abs t=\frac{1/n}{\sqrt{t^2+1/n}+\abs t}\le\frac1{\sqrt n}
$$

だから $\abs t$ に一様収束する。極限は $1$-Lipschitz だが原点で微分可能でない。Arzelà–Ascoli の定理が与える収束は関数値の一様収束であり、導関数の収束まで必要なら導関数の族にも同程度連続性などの評価を課す必要がある。
:::

## 演習問題

::::exercise{level=basic}
$\Z$ に通常の距離を入れた空間が完備であることを示せ。
:::solution
$(n_k)$ を $\Z$ の Cauchy 列とする。$\eps=1$ に対する $N$ を取ると、$k,l\ge N$ で $\abs{n_k-n_l}<1$ である。異なる二整数の差の絶対値は $1$ 以上なので $n_k=n_l$ となり、列は $N$ 番目以降定数 $n_N$ である。従って $n_N$ に収束する。（$\Z$ が完備な $\R$ の閉集合であることから :ref[prop-top-complete-closed] を使ってもよい。）
:::
::::

::::exercise{level=basic}
$\Q\cap[0,1]$ に通常の距離を入れた空間の完備化を求めよ。
:::solution
$\Q\cap[0,1]$ は完備な $\R$ の部分空間であり、その閉包は $[0,1]$ である。実際 $[0,1]$ は閉で $\Q\cap[0,1]$ を含み、$[0,1]$ の各点 $x$ の近くには有理数の稠密性から $[0,1]$ 内の有理数がある（$x$ が端点でも片側から近づける）。:ref[prop-top-completion-closure] により、$[0,1]$ と包含写像の組が完備化であり、一意性（:ref[cor-top-completion-unique]）から完備化は $[0,1]$ と等長である。
:::
::::

::::exercise{level=basic}
$T(x)=(x+4)/2$ を $\R$ 上で考え、$x_0=0$ から反復する。不動点、$x_1,x_2,x_3$、$x_3$ の実際の誤差、および事前評価による $x_3$ の誤差上界を求めよ。
:::solution
$\abs{T(x)-T(y)}=\abs{x-y}/2$ なので縮小率 $q=1/2$ の縮小写像である。$p=(p+4)/2$ から $p=4$ である。$x_1=2$、$x_2=3$、$x_3=7/2$ なので実際の誤差は $4-7/2=1/2$ である。$d(x_0,x_1)=2$ だから事前評価は $q^3/(1-q)\cdot2=(1/8)\cdot2\cdot2=1/2$ であり、この例では評価が等号で成り立つ。一般に $x_n=4(1-2^{-n})$ で $p-x_n=4\cdot2^{-n}=q^n/(1-q)\cdot d(x_0,x_1)$ である。
:::
::::

::::exercise{level=basic}
有界数列 $a=(1/i)_{i\ge1}$ の第 $N$ 項より後を $0$ にした数列 $a^{[N]}$ は、上限距離で $a$ からどれだけ離れているか。また定数列 $(1,1,\dots)$ に同じ操作をするとどうなるか。
:::solution
$d(a,a^{[N]})=\sup_{i>N}1/i=1/(N+1)$ であり、$N\to\infty$ で $0$ に近づく。定数列 $b=(1,1,\dots)$ では $d(b,b^{[N]})=\sup_{i>N}1=1$ がすべての $N$ で成り立ち、打ち切りでは近似できない。これは $a\in c_0$、$b\notin c_0$ に対応し、$c_{00}$ の完備化が $c_0$ であって $\ell^\infty$ でないことを表している（:ref[exm-top-completion-examples] の 5）。
:::
::::

::::exercise{level=basic}
初期値問題 $y'=2ty$、$y(0)=1$ の Picard 反復 $y_k$ を求め、$\abs t\le1/2$ での縮小率と、$y_2$ の誤差の事前評価および実際の誤差を比較せよ。
:::hint
$y_k(t)=\sum_{j=0}^kt^{2j}/j!$ を帰納法で示す。
:::
:::solution
$y_0=1$ とすると $y_1=1+\int_0^t2s\,\dd s=1+t^2$ である。$y_k=\sum_{j\le k}t^{2j}/j!$ と仮定すると

$$
y_{k+1}(t)=1+\int_0^t2s\sum_{j=0}^k\frac{s^{2j}}{j!}\,\dd s=1+\sum_{j=0}^k\frac{2t^{2j+2}}{(2j+2)\,j!}=1+\sum_{j=0}^k\frac{t^{2(j+1)}}{(j+1)!}
$$

となり、帰納法が完成する。極限は $e^{t^2}$ である。$I=[-1/2,1/2]$ 上で $Tu(t)=1+\int_0^t2su(s)\,\dd s$ とすると、$\abs{Tu(t)-Tv(t)}\le\abs{\int_0^t2\abs s\,\dd s}\,d_\infty(u,v)=t^2d_\infty(u,v)\le\tfrac14d_\infty(u,v)$ なので縮小率は $q=1/4$ である。$d_\infty(y_0,y_1)=\sup t^2=1/4$ なので、事前評価は $d_\infty(y_2,y)\le q^2/(1-q)\cdot\tfrac14=\tfrac1{48}\approx0.0208$ である。実際の誤差は $t=\pm1/2$ で最大となり、$e^{1/4}-(1+\tfrac14+\tfrac1{32})\approx0.00278$ である。
:::
::::


::::exercise{level=advanced}
距離空間 $X$ が全有界であることと、その完備化 $\widehat X$ が全有界であることが同値であることを示せ。これを用いて、$\widehat X$ がコンパクトであるための条件を与えよ。
:::solution
$X$ を $\iota(X)$ と同一視する。$X$ が全有界なら、$\eps>0$ に対し半径 $\eps/2$ の有限個の球 $B(x_j,\eps/2)$ で $X$ が覆われる。$\widehat X$ の任意の点 $p$ には $d(p,x)<\eps/2$ となる $x\in X$ があり、$x$ はある $B(x_j,\eps/2)$ に属するので $\widehat d(p,x_j)<\eps$ である。従って $\widehat X$ は全有界である。

逆に $\widehat X$ が全有界なら、半径 $\eps/2$ の有限個の球で $\widehat X$ を覆う。$X$ と交わる各球から $X$ の点 $x_j$ を一つ選ぶ（交わらない球は捨てる）。$X$ の任意の点 $x$ はある球に属し、その球から選んだ $x_j$ との距離は球の中心を経由する三角不等式で $\eps$ 未満なので、$X\subset\bigcup_jB(x_j,\eps)$ となり、$X$ は全有界である。

$\widehat X$ は完備なので、:ref[thm-top-metric-compactness] により $\widehat X$ がコンパクトであることは $\widehat X$ が全有界であること、すなわち $X$ が全有界であることと同値である。
:::
::::

::::exercise{level=advanced}
空でない完備距離空間 $X$ が可算集合なら、$X$ は孤立点をもつことを示せ。これを用いて、孤立点をもたない空でない完備距離空間は非可算であることを示し、$[0,1]$ に適用せよ。
:::solution
$X=\{x_1,x_2,\dots\}$（有限でもよい）と書くと、$X=\bigcup_n\{x_n\}$ は閉集合の可算和である。:ref[cor-top-baire-closed-cover] により、ある $\{x_n\}$ は内部をもつ。一点集合の内部が空でなければ $\{x_n\}$ 自身が開集合であり、$x_n$ は孤立点である。対偶を取れば、孤立点をもたない空でない完備距離空間は可算でない。$[0,1]$ は完備（$\R$ の閉集合）で、任意の点 $x$ の任意の近傍が $x$ 以外の点を含むので孤立点をもたない。従って $[0,1]$ は非可算である。
:::
::::

::::exercise{level=advanced}
$T(x)=\cos x$ を $\R$ 上で考える。$T$ は縮小写像でないが、$T^2$ は縮小写像であることを示し、方程式 $x=\cos x$ が $\R$ でただ一つの解をもつことを結論せよ。
:::solution
$\abs{\cos x-\cos y}\le\abs{x-y}$ であるが、$h\to0$ のとき $\abs{\cos(\pi/2+h)-\cos(\pi/2)}/\abs h=\abs{\sin h}/\abs h\to1$ なので、$1$ 未満の縮小率は取れない。一方 $T^2(x)=\cos(\cos x)$ の導関数は $\sin(\cos x)\sin x$ である。$\abs{\cos x}\le1$ と $\sin$ が $[-1,1]$ で増加することから $\abs{\sin(\cos x)}\le\sin1$ であり、$\abs{\sin x}\le1$ なので、平均値の定理により $\abs{T^2x-T^2y}\le\sin1\cdot\abs{x-y}$ である。$\sin1<1$ なので $T^2$ は完備な $\R$ 上の縮小写像であり、:ref[cor-top-contraction-power] により $T$ はただ一つの不動点をもつ。本文の例により、その値は約 $0.739085$ である。
:::
::::

::::exercise{level=advanced}
$[0,1]$ 上の関数族 $h_a(x)=\sin(ax)$（$0\le a\le3$）が $C([0,1])$ のコンパクト集合であることを二通りに示せ。
:::solution
**第一の方法。** $\abs{\sin u-\sin v}\le\abs{u-v}$ より $\abs{h_a(x)-h_b(x)}\le\abs{a-b}\,x\le\abs{a-b}$ なので $d_\infty(h_a,h_b)\le\abs{a-b}$ である。写像 $[0,3]\to C([0,1])$、$a\mapsto h_a$ は連続であり、コンパクトな $[0,3]$ の連続像としてコンパクトである。

**第二の方法。** 各点で $\abs{h_a(x)}\le1$ であり、$\abs{h_a(x)-h_a(y)}\le a\abs{x-y}\le3\abs{x-y}$ なので族は同程度連続である。Arzelà–Ascoli の定理により相対コンパクトである。族が閉であることを示す。$h_{a_n}\to g$ と一様収束するとき、$[0,3]$ のコンパクト性から部分列 $a_{n_k}\to a\in[0,3]$ を取れ、$d_\infty(h_{a_{n_k}},h_a)\le\abs{a_{n_k}-a}\to0$ なので $g=h_a$ である。従って族は閉で相対コンパクト、すなわちコンパクトである。
:::
::::


::::exercise{level=advanced}
$[0,\pi]$ 上の $f_n(x)=\sin(nx)$ は各点で有界であるが、同程度連続でないことを示せ。さらに、どの部分列も一様収束しないことを示せ。
:::solution
$x=0$ で同程度連続でないことを示す。$0$ の任意の近傍は、十分大きい $n$ について点 $\pi/(2n)$ を含み、$\abs{f_n(\pi/(2n))-f_n(0)}=1$ である。従って $\eps=1$ に対して条件を満たす近傍はない。

部分列 $(f_{n_k})$ が一様収束すれば、$\{f_{n_k}\}$ とその極限からなる集合はコンパクトであり（収束列とその極限の集合はコンパクトである）、$\{f_{n_k}\}$ は相対コンパクトである。Arzelà–Ascoli の定理の必要性により $\{f_{n_k}\}$ は同程度連続でなければならない。しかし $n_k\to\infty$ なので、上と同じ議論で $\{f_{n_k}\}$ も $x=0$ で同程度連続でない。これは矛盾である。
:::
::::

::::exercise{level=advanced}
$\Q$ が $\R$ の $G_\delta$ 集合でないことを、Baire の定理から直接示せ。
:::solution
$\Q=\bigcap_nU_n$（$U_n$ は開）と仮定する。各 $U_n$ は $\Q$ を含むので稠密である。有理数を $q_1,q_2,\dots$ と並べると、$V_k=\R\setminus\{q_k\}$ も稠密な開集合である。可算個の稠密開集合 $U_n,V_k$ の共通部分は

$$
\bigcap_nU_n\cap\bigcap_kV_k=\Q\cap(\R\setminus\Q)=\varnothing
$$

となるが、Baire の定理（:ref[thm-top-baire]）により完備距離空間 $\R$ ではこれは稠密、特に空でないはずである。これは矛盾である。
:::
::::

::::exercise{level=advanced}
$K$ をコンパクト距離空間とし、空でない閉集合 $A\subset K$ に対する関数 $d_A(x)=d(x,A)$ 全体の族 $\mathcal D$ を考える。$\mathcal D$ が $C(K)$ で相対コンパクトであることを示し、閉集合の列 $(A_n)$ に対して $(d_{A_n})$ が一様収束する部分列をもつことを結論せよ。
:::solution
:ref[prop-top-distance-set] により各 $d_A$ は $1$-Lipschitz なので、$\mathcal D$ は同程度連続である。また $0\le d_A(x)\le\operatorname{diam}K<\infty$（コンパクト距離空間は有界である）なので各点で有界である。Arzelà–Ascoli の定理により $\mathcal D$ は相対コンパクトであり、特に任意の列 $(d_{A_n})$ は一様収束する部分列をもつ。
:::
::::

## 確認テスト

::::quiz
有限個の成分を除いて $0$ である実数列全体 $c_{00}$ に上限距離を入れた空間の完備化はどれか。

- [ ] 有界な実数列全体 $\ell^\infty$
- [x] $0$ に収束する実数列全体 $c_0$
- [ ] 最終的に定数になる実数列全体

:::explanation
$c_0$ は完備な $\ell^\infty$ の閉集合なので完備であり、$a\in c_0$ は打ち切り $a^{[N]}\in c_{00}$ で $d(a,a^{[N]})=\sup_{i>N}\abs{a_i}\to0$ と近似できるので $c_{00}$ は $c_0$ で稠密である。:ref[prop-top-completion-closure] により完備化は $c_0$ である。定数列 $(1,1,\dots)$ はどの $c_{00}$ の元とも距離 $1$ 以上離れているので、$\ell^\infty$ 全体や最終的に定数になる列の全体は完備化にならない。
:::
::::

::::quiz
$\R$（通常の距離）の部分集合のうち、第一類であるものをすべて選べ。

- [x] 有理数全体 $\Q$
- [x] Cantor 集合
- [ ] 無理数全体 $\R\setminus\Q$
- [ ] 閉区間 $[0,1]$

:::explanation
$\Q$ は一点集合（疎集合）の可算和なので第一類である。Cantor 集合は内部が空な閉集合なので疎集合であり、第一類である。$\R\setminus\Q$ が第一類なら $\R=\Q\cup(\R\setminus\Q)$ も第一類となり、Baire の定理に反する。$[0,1]$ は空でない開区間 $(0,1)$ を含むので、Baire の定理により第一類でない（第一類集合の部分集合は第一類だからである）。
:::
::::

::::quiz
完備距離空間で、空でない閉集合の減少列であって共通部分が空になりうるものをすべて選べ。

- [x] $\R$ で $F_n=[n,\infty)$
- [ ] $\R$ で $F_n=[-1/n,1/n]$
- [x] $\ell^\infty$ で $F_n=\{e_k:k\ge n\}$（$e_k$ は第 $k$ 成分だけ $1$ の数列）
- [ ] $\R$ で $F_n=[0,1+1/n]$

:::explanation
$[n,\infty)$ の共通部分は空であり、直径は $\infty$ である。$\{e_k:k\ge n\}$ は閉集合で（距離 $1$ の点からなるので収束列はやがて定数）、共通部分は空であり、直径は $1$ のまま $0$ に近づかない。一方 $[-1/n,1/n]$ は直径が $0$ に近づくので Cantor の共通部分定理により一点 $\{0\}$ を共有し、$[0,1+1/n]$ は $[0,1]$ を共通部分にもつ。直径が $0$ に近づかない場合や有界でない場合には、完備でも共通点は保証されない。
:::
::::

::::quiz{answer="0.25" tol="1e-6"}
縮小率 $q=1/2$ の縮小写像の反復で、連続する二項 $x_n,x_{n+1}$ の距離が $1/8$ であった。$x_n$ から不動点までの距離の、事後評価による上界を求めよ。

:::explanation
事後評価 $d(x_n,p)\le d(x_n,x_{n+1})/(1-q)$ に代入すると $(1/8)/(1-1/2)=1/4=0.25$ である。この評価は $d(x_n,p)\le d(x_n,x_{n+1})+d(Tx_n,Tp)\le d(x_n,x_{n+1})+q\,d(x_n,p)$ を整理して得られる。
:::
::::

::::quiz{answer="0.2222222222" tol="1e-6"}
初期値問題 $y'=y^2$、$y(0)=1$ に Picard–Lindelöf の定理を $a=1$、$b=1/2$ で適用する。定理が解の存在を保証する区間 $[-h,h]$ の $h=\min\{a,b/M\}$ を小数で求めよ。

:::explanation
$\Omega=[-1,1]\times[1/2,3/2]$ 上で $\abs{f(t,y)}=y^2\le(3/2)^2=9/4$ なので $M=9/4$ である。$h=\min\{1,\ (1/2)/(9/4)\}=2/9\approx0.2222$ である。真の解 $1/(1-t)$ は $t<1$ で存在するので、保証される区間は実際の存在区間よりかなり小さい。
:::
::::

::::quiz
コンパクト区間上の実数値連続関数の族が一様距離について相対コンパクトであるための条件として正しいものを選べ。

- [x] 同程度連続であり、各点で有界である。
- [ ] 各関数が連続であり、各点で有界である。
- [ ] 同程度連続である。

:::explanation
Arzelà–Ascoli の定理により、条件は同程度連続性と各点での相対コンパクト性であり、$\R$ では後者は各点での有界性と同値である（Heine–Borel の定理）。二番目の選択肢には $f_n(t)=t^n$ が反例である（各点で有界だが $t=1$ で同程度連続でなく、一様収束する部分列をもたない）。三番目には定数関数の列 $g_n=n$ が反例である（同程度連続だが各点で有界でない）。
:::
::::

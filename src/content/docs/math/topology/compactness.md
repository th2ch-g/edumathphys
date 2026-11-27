---
title: コンパクト性
description: 開被覆の有限部分被覆からコンパクト性を定義し、Heine–Borel の定理、距離空間でのコンパクト性の同値条件、超フィルターによる Tychonoff の定理、局所コンパクト空間の一点コンパクト化を証明する。
notation:
  - symbol: '\mathcal U'
    meaning: 開被覆
  - symbol: '\operatorname{diam} A'
    meaning: 部分集合 A の直径
  - symbol: 'B(x,r)'
    meaning: 中心 x、半径 r の開球
  - symbol: '\mathscr F,\ \mathscr U'
    meaning: 集合上のフィルター・超フィルター
  - symbol: '\pi_i'
    meaning: 積空間の第 i 射影
  - symbol: '\alpha X'
    meaning: 局所コンパクト Hausdorff 空間 X の一点コンパクト化
  - symbol: '\infty'
    meaning: 一点コンパクト化で加える点
  - symbol: 'Q(x)'
    meaning: 点 x の準成分
---

## この章の目標

- 開被覆の有限部分被覆と閉集合族の有限交叉性という二つの言葉でコンパクト性を扱い、閉集合・Hausdorff 空間・連続写像との関係を証明する。
- 閉区間のコンパクト性から Heine–Borel の定理を導き、コンパクト空間上の連続関数の最大値・最小値の存在と一様連続性を示す。
- 距離空間ではコンパクト・極限点コンパクト・点列コンパクト・完備かつ全有界の四条件が同値であることを、Lebesgue 数の補題とともに証明する。
- 超フィルターを用いて、任意個のコンパクト空間の積がコンパクトであること（Tychonoff の定理）を証明する。
- 局所コンパクト Hausdorff 空間に一点を加えてコンパクト化し、その一意性と固有写像との関係を調べる。

[解析学](/math/analysis/) では、有界閉区間 $[a,b]$ 上の連続関数が最大値・最小値をとること、一様連続であることを学んだ。
その証明を振り返ると、区間の各点の近くで連続性から得た局所的な情報を、有限個に絞って区間全体に広げる操作が本質的であった。
有限個ならば誤差の許容値の最小値や有界性の定数の最大値を取れるが、無限個では取れないからである。
この「局所的に選んだ開集合を有限個に絞れる」という性質を取り出したものがコンパクト性である。
本章ではコンパクト性を一般の位相空間で定義し、距離空間での具体的な判定法、任意個の積での保存、局所的な版とそこからのコンパクト化を順に扱う。

本章では Hausdorff 空間を頻繁に用いる。
詳しい扱いは [分離公理](/math/topology/separation/) で行うが、ここで定義と最小限の性質を述べておく。

:::definition[Hausdorff 空間]
位相空間 $X$ の異なる任意の二点 $x\ne y$ に対し、$x\in U$、$y\in V$、$U\cap V=\varnothing$ を満たす開集合 $U,V$ が存在するとき、$X$ を Hausdorff 空間という。
:::

距離空間は Hausdorff である。
実際 $x\ne y$ なら $r=d(x,y)/2>0$ と置くと $B(x,r)\cap B(y,r)=\varnothing$ となる。
両方に属する $z$ があれば $d(x,y)\le d(x,z)+d(z,y)<2r=d(x,y)$ となって矛盾するからである。
Hausdorff 空間の部分空間 $A$ も Hausdorff である。
$A$ の二点を $X$ で分離する開集合 $U,V$ を $A$ と交わせば、$A\cap U$ と $A\cap V$ が $A$ でそれらを分離する。
また Hausdorff 空間では一点集合は閉である。
$y\ne x$ ごとに $x$ を含まない $y$ の開近傍があり、それらの和が $X\setminus\{x\}$ になるからである。

## コンパクト空間

閉区間 $[0,1]$ 上の連続関数 $f$ が有界であることを、局所的な情報から示してみる。
各点 $x$ で連続性を $\eps=1$ に使うと、$x$ の開近傍 $U_x$ 上で $|f|\le|f(x)|+1$ となる。
近傍は無数にあるが、有限個 $U_{x_1},\dots,U_{x_n}$ で区間を覆えるなら、$|f|\le\max_j|f(x_j)|+1$ という一つの上界が得られる。
これに対し $(0,1]$ 上の $f(x)=1/x$ では、各点 $x$ の近傍 $(x/2,2x)\cap(0,1]$ 上で $f$ は有界であるのに、どの有限個を選んでも $0$ の近くが覆われず、実際 $f$ は有界でない。
二つの空間の違いは、開集合で覆ったときに有限個を選び出せるかどうかにある。

:::definition[開被覆とコンパクト性]{#def-top-compact}
$X$ を位相空間とし、$K\subset X$ とする。
$X$ の開集合の族 $\mathcal U=(U_\lambda)_{\lambda\in\Lambda}$ が $K\subset\bigcup_{\lambda\in\Lambda}U_\lambda$ を満たすとき、$\mathcal U$ を $K$ の:term[開被覆]{en="open cover" yomi="かいひふく"}という。
$\Lambda'\subset\Lambda$ について $(U_\lambda)_{\lambda\in\Lambda'}$ もまた $K$ を覆うとき、これを部分被覆といい、$\Lambda'$ が有限のとき有限部分被覆という。

$X$ のすべての開被覆が有限部分被覆をもつとき、$X$ は:term[コンパクト]{en="compact" yomi="こんぱくと"}であるという。
部分集合 $K\subset X$ は、相対位相を入れた部分空間としてコンパクトであるとき、コンパクト集合であるという。
:::

部分集合のコンパクト性は、周囲の空間の開集合だけで判定できる。

:::lemma[部分集合のコンパクト性の判定]{#lem-top-compact-subspace}
$K\subset X$ がコンパクトであるための必要十分条件は、$X$ の開集合からなる $K$ の任意の開被覆が有限部分被覆をもつことである。
:::

:::proof
$K$ がコンパクトとし、$X$ の開集合の族 $(U_\lambda)$ が $K$ を覆うとする。
各 $K\cap U_\lambda$ は $K$ の相対開集合であり（:ref[def-top-subspace]）、その和は $K$ に等しい。
コンパクト性から有限個 $\lambda_1,\dots,\lambda_n$ で $K=\bigcup_j(K\cap U_{\lambda_j})$ となり、従って $K\subset\bigcup_jU_{\lambda_j}$ である。

逆に条件を仮定し、$(V_\lambda)$ を $K$ の相対開集合による被覆とする。
相対位相の定義により各 $\lambda$ に $X$ の開集合 $U_\lambda$ を選んで $V_\lambda=K\cap U_\lambda$ と書ける。
$(U_\lambda)$ は $K$ を覆うので、仮定から有限個 $U_{\lambda_1},\dots,U_{\lambda_n}$ が $K$ を覆う。
このとき $K=\bigcup_j(K\cap U_{\lambda_j})=\bigcup_jV_{\lambda_j}$ である。
:::

この補題により、以後はコンパクト集合の被覆を周囲の空間の開集合で考えてよい。
また $K\subset A\subset X$ のとき、$K$ に $A$ から入る相対位相と $X$ から入る相対位相は一致するので、$K$ のコンパクト性は $A$ と $X$ のどちらで考えても同じである。

:::remark[Hausdorff 性を定義に含める流儀]
Bourbaki の流儀では上の性質を準コンパクトと呼び、Hausdorff 性を合わせたものをコンパクトと呼ぶ。
このサイトでは Hausdorff 性を定義に含めず、必要なときは仮定として明記する。
空集合は空な部分族で覆えるのでコンパクトである。
:::

:::example[コンパクトな空間とそうでない空間]{#exm-top-compact-basic}
1. 有限集合はどの位相でもコンパクトである。開被覆が与えられたら、各点についてその点を含む元を一つずつ選べばよい。
2. 無限集合に離散位相を入れた空間はコンパクトでない。一点集合全体 $\{\{x\}:x\in X\}$ は開被覆であり、その有限部分族は有限個の点しか覆わないからである。
3. $\R$ はコンパクトでない。開被覆 $\{(-n,n):n\ge1\}$ の有限部分族は最大の添字 $N$ の $(-N,N)$ に含まれ、$N$ を覆わない。
4. $(0,1]$ はコンパクトでない。開被覆 $\{(1/n,2):n\ge2\}$ の有限部分族の最大の添字を $N$ とすると、その和は $(1/N,2)$ であり $1/N$ が覆われない。
5. 密着空間はコンパクトである。開集合は $\varnothing$ と $X$ だけなので、被覆には $X$ 自身が含まれ（$X\ne\varnothing$ のとき）、それ一つで覆える。
:::

:::example[補有限位相ではすべての部分集合がコンパクトである]{#exm-top-cofinite-compact}
無限集合 $X$ に補有限位相を入れる。
空でない $K\subset X$ と、$X$ の開集合による $K$ の被覆 $\mathcal U$ を取る。
$K$ の一点を含む $U_0\in\mathcal U$ は空でない開集合なので、$X\setminus U_0$ は有限集合である。
従って $K\setminus U_0$ は有限個の点 $y_1,\dots,y_m$ からなり、各 $y_j$ を含む $U_j\in\mathcal U$ を選べば $U_0,U_1,\dots,U_m$ が $K$ を覆う。
よってすべての部分集合がコンパクトである。
一方、閉集合は有限集合と $X$ 自身だけなので、$X$ の無限の真部分集合はコンパクトであるが閉でない。
後で示すように Hausdorff 空間ではコンパクト集合は閉であり、この例はその仮定を外せないことを示している。
実際この空間は Hausdorff でない。空でない二つの開集合は補集合が有限なので必ず交わる。
:::

### 有限交叉性による言い換え

開被覆の条件を補集合に移すと、閉集合族の言葉になる。

:::definition[有限交叉性]
集合族 $\mathcal F$ の任意の有限個の元 $F_1,\dots,F_n$（$n\ge1$）について $F_1\cap\dots\cap F_n\ne\varnothing$ が成り立つとき、$\mathcal F$ は:term[有限交叉性]{en="finite intersection property" yomi="ゆうげんこうさせい"}をもつという。
:::

:::theorem[閉集合の有限交叉性によるコンパクト性]{#thm-top-compact-fip}
位相空間 $X$ がコンパクトであるための必要十分条件は、有限交叉性をもつ閉集合の空でない任意の族 $\mathcal F$ について $\bigcap_{F\in\mathcal F}F\ne\varnothing$ となることである。
:::

:::proof
De Morgan の法則により、閉集合族 $\mathcal F$ について

$$
\bigcup_{F\in\mathcal F}(X\setminus F)=X\setminus\bigcap_{F\in\mathcal F}F
$$

が成り立つ。

$X$ がコンパクトであるとし、有限交叉性をもつ閉集合の空でない族 $\mathcal F$ の共通部分が空であると仮定する。
上の式から $\{X\setminus F:F\in\mathcal F\}$ は $X$ の開被覆である。
コンパクト性によって有限個 $F_1,\dots,F_n$ で $X=\bigcup_j(X\setminus F_j)=X\setminus\bigcap_jF_j$ となり、$F_1\cap\dots\cap F_n=\varnothing$ が従う。
$\mathcal F$ が空でないので $n\ge1$ と取れ、これは有限交叉性に反する。

逆に条件を仮定し、$\mathcal U$ を $X$ の開被覆とする。
$X=\varnothing$ なら空な部分族で覆えるので、$X\ne\varnothing$ とする。このとき $\mathcal U\ne\varnothing$ である。
$\mathcal U$ が有限部分被覆をもたないと仮定し、$\mathcal F=\{X\setminus U:U\in\mathcal U\}$ と置く。
任意の $U_1,\dots,U_n\in\mathcal U$ について $U_1\cup\dots\cup U_n\ne X$ だから $\bigcap_j(X\setminus U_j)\ne\varnothing$ であり、$\mathcal F$ は有限交叉性をもつ閉集合の空でない族である。
仮定から $\bigcap\mathcal F\ne\varnothing$ であるが、上の式からこれは $\bigcup\mathcal U\ne X$ を意味し、$\mathcal U$ が被覆であることに反する。
:::

:::corollary[減少閉集合列の共通部分]{#cor-top-cantor-compact}
1. コンパクト空間 $X$ の空でない閉集合の減少列 $F_1\supset F_2\supset\cdots$ について $\bigcap_nF_n\ne\varnothing$ である。
2. Hausdorff 空間の空でないコンパクト集合の減少列 $K_1\supset K_2\supset\cdots$ について $\bigcap_nK_n\ne\varnothing$ である。
:::

:::proof
1. 減少列では $F_{n_1}\cap\dots\cap F_{n_k}=F_{\max_jn_j}\ne\varnothing$ なので有限交叉性があり、:ref[thm-top-compact-fip] を適用すればよい。
2. 後で示す :ref[prop-top-compact-closed] により各 $K_n$ は閉である。従って $K_n=K_n\cap K_1$ はコンパクト空間 $K_1$ の空でない閉部分集合であり、1 を $K_1$ に適用すればよい。
:::

コンパクト性の仮定は外せない。
$\R$ の閉集合 $[n,\infty)$ は空でなく減少するが共通点をもたない。
$(0,1]$ の相対閉集合 $(0,1/n]$ も同様である。

### Hausdorff 空間のコンパクト集合

Hausdorff 空間では、コンパクト集合は点と同じように開集合で分離できる。
有限個の点ならば分離に使う開集合を有限回交わせばよいが、コンパクト集合でも有限個に帰着できるからである。

:::lemma[コンパクト集合の開集合による分離]{#lem-top-compact-separation}
$X$ を Hausdorff 空間とする。

1. コンパクト集合 $K\subset X$ と点 $x\notin K$ に対し、$K\subset U$、$x\in V$、$U\cap V=\varnothing$ を満たす開集合 $U,V$ が存在する。
2. 交わらないコンパクト集合 $K,L\subset X$ に対し、$K\subset U$、$L\subset V$、$U\cap V=\varnothing$ を満たす開集合 $U,V$ が存在する。
:::

:::proof
1. $K=\varnothing$ なら $U=\varnothing$、$V=X$ とすればよいので、$K\ne\varnothing$ とする。
各 $y\in K$ は $x$ と異なるので、Hausdorff 性から $y\in U_y$、$x\in V_y$、$U_y\cap V_y=\varnothing$ となる開集合を選べる。
$(U_y)_{y\in K}$ は $K$ の開被覆なので、:ref[lem-top-compact-subspace] により有限個 $y_1,\dots,y_n$ で $K\subset U_{y_1}\cup\dots\cup U_{y_n}$ となる。
$U=\bigcup_jU_{y_j}$、$V=\bigcap_jV_{y_j}$ と置くと、両者は開であり（$V$ は有限個の開集合の共通部分）、$K\subset U$、$x\in V$ である。
$z\in U\cap V$ があれば、ある $j$ で $z\in U_{y_j}$ かつ $z\in V\subset V_{y_j}$ となり $U_{y_j}\cap V_{y_j}=\varnothing$ に反する。

2. $L=\varnothing$ なら $U=X$、$V=\varnothing$ とすればよい。
$L\ne\varnothing$ とし、各 $x\in L$ に 1 を使って $K\subset U_x$、$x\in V_x$、$U_x\cap V_x=\varnothing$ となる開集合を選ぶ。
$(V_x)_{x\in L}$ は $L$ の開被覆なので有限個 $x_1,\dots,x_m$ で $L\subset\bigcup_kV_{x_k}$ となる。
$U=\bigcap_kU_{x_k}$、$V=\bigcup_kV_{x_k}$ と置けば、1 と同じ議論で $K\subset U$、$L\subset V$、$U\cap V=\varnothing$ である。
:::

:::proposition[閉部分集合・有限和・Hausdorff 空間のコンパクト集合]{#prop-top-compact-closed}
1. コンパクト空間の閉部分集合はコンパクトである。
2. 有限個のコンパクト部分集合の和はコンパクトである。
3. Hausdorff 空間のコンパクト部分集合は閉である。
4. Hausdorff 空間では、空でない任意個のコンパクト集合の共通部分はコンパクトである。
:::

:::proof
1. $X$ をコンパクト、$F\subset X$ を閉とし、$X$ の開集合による $F$ の被覆 $\mathcal U$ を取る。
$X\setminus F$ は開なので、$\mathcal U\cup\{X\setminus F\}$ は $X$ の開被覆である。
その有限部分被覆を取り、$X\setminus F$ が含まれていれば取り除く。
残った $\mathcal U$ の有限個の元は、$F$ の点を覆うのに $X\setminus F$ を必要としないので $F$ を覆う。
:ref[lem-top-compact-subspace] により $F$ はコンパクトである。

2. $K_1,\dots,K_m$ をコンパクトとし、$\mathcal U$ を $\bigcup_jK_j$ の $X$ の開集合による被覆とする。
各 $K_j$ を覆う有限部分族 $\mathcal U_j$ を選べば、$\mathcal U_1\cup\dots\cup\mathcal U_m$ は和集合を覆う有限部分族である。

3. $K$ をコンパクトとし、$x\notin K$ とする。
:ref[lem-top-compact-separation] の 1 により、$x$ の開近傍 $V$ で $K$ と交わらないものがある（$V\cap U=\varnothing$、$K\subset U$）。
従って $X\setminus K$ は各点の開近傍を含み、開である。

4. $(K_\lambda)_{\lambda\in\Lambda}$（$\Lambda\ne\varnothing$）をコンパクト集合の族とし、$\lambda_0\in\Lambda$ を一つ固定する。
3 により各 $K_\lambda$ は閉なので $C=\bigcap_\lambda K_\lambda$ は閉であり、$C=C\cap K_{\lambda_0}$ は $K_{\lambda_0}$ の閉部分集合である。
$K_{\lambda_0}$ はそれ自身コンパクト空間なので、1 により $C$ はコンパクトである。
:::

:::remark[よくある誤り]
「コンパクト集合は閉」「コンパクト集合の共通部分はコンパクト」はどちらも Hausdorff 性を使っている。
:ref[exm-top-cofinite-compact] では無限の真部分集合がコンパクトで閉でなかった。
共通部分の主張が非 Hausdorff 空間で失敗する例は演習問題で扱う。
また逆向きの「閉集合はコンパクト」はコンパクト空間の中でしか成り立たない。$\R$ 自身は閉だがコンパクトでない。
:::

### 閉区間と Heine–Borel の定理

$\R$ の部分集合のコンパクト性は実数の連続性に支えられている。
まず閉区間を直接扱う。

:::theorem[閉区間のコンパクト性]{#thm-top-interval-compact}
$a\le b$ を実数とする。閉区間 $[a,b]$ はコンパクトである。
:::

:::proof
$\mathcal U$ を $\R$ の開集合による $[a,b]$ の被覆とし、:ref[lem-top-compact-subspace] を用いる。

$$
S=\{t\in[a,b]:[a,t]\text{ は }\mathcal U\text{ の有限個の元で覆える}\}
$$

と置く。
$a$ を含む $U\in\mathcal U$ が一つで $[a,a]=\{a\}$ を覆うので $a\in S$ である。
$S$ は空でなく $b$ で上に有界なので、[実数の連続性](/math/analysis/real-numbers/) により上限 $c=\sup S$ が存在し、$a\le c\le b$ である。

$c\in[a,b]$ を含む $U_0\in\mathcal U$ を取る。$U_0$ は開なので、ある $\eps>0$ で $(c-\eps,c+\eps)\subset U_0$ となる。
$c$ は $S$ の上限なので、$c-\eps<t\le c$ を満たす $t\in S$ がある。
$[a,t]$ を覆う有限部分族を $\mathcal U'$ とし、$t'=\min(c+\eps/2,\,b)$ と置く。
$t\le c<c+\eps/2$ かつ $t\le b$ なので $t\le t'$ であり、$[t,t']\subset(c-\eps,c+\eps)\subset U_0$ である。
従って $\mathcal U'\cup\{U_0\}$ は $[a,t']$ を覆い、$t'\in S$ である。
上限の性質から $t'\le c$ であるが、$c+\eps/2>c$ なので $t'=b$ でなければならない。
よって $b\in S$、すなわち $[a,b]$ は有限個の元で覆える。
:::

実数の連続性が本質的であることは、有理数の中で同じことを考えると分かる。
$\Q\cap[0,1]$ はコンパクトでない。
無理数 $s=1/\sqrt2\in[0,1]$ について、$\Q\cap[0,1]$ の相対開集合

$$
U_n=\{q\in\Q\cap[0,1]:|q-s|>1/n\}\qquad(n\ge1)
$$

は $s$ が有理数でないので $\Q\cap[0,1]$ を覆うが、有限個の和は最大の添字 $N$ の $U_N$ であり、$|q-s|\le1/N$ を満たす有理数 $q\in[0,1]$（有理数の稠密性により存在する）を覆わない。

:::theorem[Heine–Borel の定理]{#thm-top-heine-borel}
$\R^n$ の部分集合 $K$ がコンパクトであるための必要十分条件は、$K$ が閉かつ有界であることである。
:::

:::proof
空集合はコンパクトで閉かつ有界なので、以下 $K\ne\varnothing$ とする。

**必要性。** $\R^n$ は距離空間なので Hausdorff であり、:ref[prop-top-compact-closed] の 3 により $K$ は閉である。
また開球 $B(0,m)$（$m=1,2,\dots$）の全体は $\R^n$ を覆うので、$K$ を覆う有限個 $B(0,m_1),\dots,B(0,m_k)$ がある。
これらは最大の $m=\max_jm_j$ の球 $B(0,m)$ に含まれるので、$K\subset B(0,m)$ であり $K$ は有界である。

**十分性。** まず $M\ge0$ に対し立方体 $Q=[-M,M]^n$ がコンパクトであることを示す。
$M=0$ なら $Q$ は一点なのでよい。$M>0$ とし、$Q$ の開被覆 $\mathcal U$ で有限部分被覆をもたないものがあると仮定する。
$Q_0=Q$ と置く。
$Q_k=\prod_{i=1}^n[a_i^{(k)},b_i^{(k)}]$ が一辺 $2M/2^k$ の閉立方体で、$\mathcal U$ の有限個では覆えないとする。
各辺を中点 $m_i=(a_i^{(k)}+b_i^{(k)})/2$ で $[a_i^{(k)},m_i]$ と $[m_i,b_i^{(k)}]$ に分け、各座標でどちらかを選ぶと、一辺 $2M/2^{k+1}$ の閉立方体が $2^n$ 個でき、その和は $Q_k$ である。
すべてが有限個で覆えるなら、それらの有限部分族を合わせて $Q_k$ も有限個で覆えてしまう。
従って有限個で覆えないものが少なくとも一つあり、それを $Q_{k+1}$ とする。

こうして得た $Q_0\supset Q_1\supset\cdots$ について、各座標 $i$ の区間 $[a_i^{(k)},b_i^{(k)}]$ は入れ子で長さ $2M/2^k\to0$ である。
[区間縮小法](/math/analysis/real-numbers/) により、すべての区間に属する実数 $x_i$ がある。
$x=(x_1,\dots,x_n)$ はすべての $Q_k$ に属し、特に $x\in Q$ である。
$x$ を含む $U\in\mathcal U$ を取ると、$U$ は $Q$ の相対開集合なので $Q\cap B(x,r)\subset U$ となる $r>0$ がある。
$Q_k$ の任意の二点の距離は対角線の長さ $\sqrt n\,2M/2^k$ 以下である。
$\sqrt n\,2M/2^k<r$ となる $k$ を取ると、$x\in Q_k$ なので $Q_k\subset Q\cap B(x,r)\subset U$ となり、$Q_k$ は $\mathcal U$ の一つの元で覆える。
これは $Q_k$ の選び方に反する。

最後に $K$ が閉かつ有界なら、ある $M$ で $K\subset[-M,M]^n=Q$ である。
$K$ は $\R^n$ で閉なので $K=K\cap Q$ は $Q$ で閉であり、:ref[prop-top-compact-closed] の 1 によりコンパクトである。
:::

:::example[Heine–Borel の定理の適用]{#exm-top-heine-borel}
1. 単位球面 $S^n=\{x\in\R^{n+1}:\norm x=1\}$ は、連続関数 $x\mapsto\norm x$ による閉集合 $\{1\}$ の逆像なので閉であり、有界だからコンパクトである。閉球 $\{\norm x\le r\}$ も同様である。
2. 実 $n$ 次正方行列全体を $\R^{n^2}$ と同一視する。直交群 $O(n)=\{A:\transp AA=I\}$ は、成分が多項式で与えられる連続写像 $A\mapsto\transp AA$ による一点 $\{I\}$ の逆像なので閉である。$A\in O(n)$ の各列は単位ベクトルなので $\sum_{i,j}a_{ij}^2=n$ であり、$O(n)$ は半径 $\sqrt n$ の球面に含まれて有界である。よって $O(n)$ はコンパクトであり、その閉部分集合 $SO(n)=O(n)\cap\{\det A=1\}$ もコンパクトである。
3. 一般線形群 $GL_n(\R)=\{\det A\ne0\}$ はコンパクトでない。$tI$（$t>0$）は非有界であり、また $t\to0$ のとき $tI\to0\notin GL_n(\R)$ なので閉でもない。
4. 無限集合 $X$ に離散距離（異なる二点の距離を $1$ とする）を入れると、$X$ 自身は閉で直径 $1$ の有界集合であるが、:ref[exm-top-compact-basic] の 2 によりコンパクトでない。従って「閉かつ有界」からコンパクト性を導くには、距離の公理に加えて仮定が要る。距離空間の正しい判定条件は次々節で与える。
:::

開被覆の有限部分被覆があることは、すべてのネットに収束部分ネットがあることとも同値である。開被覆からの反例ネットの作り方は :ref[prop-top-compact-net] による。

## コンパクト性と連続写像

コンパクト性は連続写像で保たれる。
これが最大値・最小値の定理の位相的な理由である。

:::theorem[コンパクト集合の連続像]{#thm-top-compact-image}
$f:X\to Y$ を連続写像とする。

1. $K\subset X$ がコンパクトなら $f(K)$ はコンパクトである。
2. （最大値・最小値の定理）$X$ が空でないコンパクト空間で $Y=\R$ なら、$f$ は有界であり、最大値と最小値をとる。
3. $X$ がコンパクトで $Y$ が Hausdorff なら $f$ は閉写像である。特に連続全単射 $f$ は同相写像であり、連続単射 $f$ は像の上への同相写像（埋め込み）である。
:::

:::proof
1. $(V_\lambda)$ を $Y$ の開集合による $f(K)$ の被覆とする。
$f$ は連続なので $f^{-1}(V_\lambda)$ は $X$ で開であり、$x\in K$ なら $f(x)\in V_\lambda$ となる $\lambda$ があるので $(f^{-1}(V_\lambda))$ は $K$ を覆う。
:ref[lem-top-compact-subspace] により有限個 $\lambda_1,\dots,\lambda_n$ で $K\subset\bigcup_jf^{-1}(V_{\lambda_j})$ となる。
任意の $y=f(x)\in f(K)$ について $x\in f^{-1}(V_{\lambda_j})$ となる $j$ があり $y\in V_{\lambda_j}$ である。
よって $V_{\lambda_1},\dots,V_{\lambda_n}$ が $f(K)$ を覆い、$f(K)$ はコンパクトである。

2. 1 により $f(X)$ は $\R$ の空でないコンパクト集合であり、:ref[thm-top-heine-borel] により閉かつ有界である。
有界なので上限 $s=\sup f(X)$ が実数として存在する。
任意の $\eps>0$ に対し上限の性質から $(s-\eps,s]$ に $f(X)$ の点があるので、$s\in\overline{f(X)}=f(X)$ である。
従って $s=f(x_0)$ となる $x_0$ があり、これが最大値である。最小値も同様である。

3. $F\subset X$ を閉とする。:ref[prop-top-compact-closed] の 1 により $F$ はコンパクト、1 により $f(F)$ はコンパクト、:ref[prop-top-compact-closed] の 3 により $f(F)$ は Hausdorff 空間 $Y$ で閉である。
よって $f$ は閉写像である。
連続全単射が閉写像なら同相写像である（:ref[prop-top-homeomorphism-test]）。
連続単射 $f$ については、$f(X)$ は Hausdorff 空間の部分空間なので Hausdorff であり、$X\to f(X)$ は連続全単射だから同相写像である。
:::

この閉写像判定を連続全射に適用すると :ref[cor-top-compact-quotient] を得る。特に円周への端点同一視が :ref[exm-top-circle-quotient] の同相を与える。

:::example[最大値・最小値の計算]
1. 円周 $S^1$ 上で $f(x,y)=x+y$ は最大値と最小値をとる。$x^2+y^2=1$ のとき $(x+y)^2+(x-y)^2=2(x^2+y^2)=2$ なので $(x+y)^2\le2$ であり、最大値 $\sqrt2$ は $(1/\sqrt2,1/\sqrt2)$、最小値 $-\sqrt2$ は $(-1/\sqrt2,-1/\sqrt2)$ でとる。
2. 距離空間 $X$ の空でないコンパクト集合 $K$ と点 $p$ に対し、$y\mapsto d(p,y)$ は連続なので $K$ 上で最小値をとる。すなわち $K$ には $p$ に最も近い点がある。$\R^n$ の空でない閉集合 $F$ でも同じことが成り立つ。$y_0\in F$ を一つ取り $r=d(p,y_0)$ とすると、$F\cap\{y:d(p,y)\le r\}$ は閉かつ有界なのでコンパクトであり、そこでの最小値が $F$ 全体での最小値になる。
3. 仮定を外すと結論は失敗する。$(0,1)$ 上の $f(x)=x$ は有界であるが最大値も最小値もとらず、$(0,1]$ 上の $1/x$ は有界でない。
:::

:::remark[閉写像になるための条件]
3 の Hausdorff 性も外せない。
$[0,1]$ の通常位相から密着位相への恒等写像は連続全単射であるが、同相写像でない。
また連続写像は一般に閉集合を閉集合に写さない。例えば $x\mapsto e^x$ による閉集合 $\R$ の像は $(0,\infty)$ である。
:::

### 管状近傍の補題と有限積

各点の上の切片を含む開集合が管を含むことは、積位相の開長方形から :ref[lem-top-tube] で証明した。ここではその結論を有限積に用いる。

:::theorem[有限積のコンパクト性]{#thm-top-finite-product-compact}
コンパクト空間 $X_1,\dots,X_n$ の積 $X_1\times\dots\times X_n$ は積位相でコンパクトである。
:::

:::proof
まず二つのコンパクト空間 $X,Y$ について示す。
$\mathcal W$ を $X\times Y$ の開被覆とする。
$x\in X$ を固定する。写像 $Y\to X\times Y$、$y\mapsto(x,y)$ は座標がともに連続なので連続であり（:ref[thm-top-product-property]）、その像 $\{x\}\times Y$ はコンパクトである（:ref[thm-top-compact-image]）。
従って $\mathcal W$ の有限個 $W^x_1,\dots,W^x_{m_x}$ が $\{x\}\times Y$ を覆う。
その和 $W^x$ は $\{x\}\times Y$ を含む開集合なので、:ref[lem-top-tube] により $U_x\times Y\subset W^x$ となる $x$ の開近傍 $U_x$ がある。
$(U_x)_{x\in X}$ は $X$ の開被覆なので、有限個 $x_1,\dots,x_k$ で $X=\bigcup_lU_{x_l}$ となる。
このとき

$$
X\times Y=\bigcup_{l=1}^kU_{x_l}\times Y\subset\bigcup_{l=1}^k\bigcup_{j=1}^{m_{x_l}}W^{x_l}_j
$$

であり、右辺は $\mathcal W$ の有限個の元の和である。

一般の $n$ については帰納法による。
$n=0$ の積は一点、$n=1$ は仮定そのものである。
$n\ge2$ のとき、写像 $(x_1,\dots,x_n)\mapsto((x_1,\dots,x_{n-1}),x_n)$ は積の普遍性（:ref[thm-top-product-property]）により両方向に連続な全単射であり、$X_1\times\dots\times X_n$ は $(X_1\times\dots\times X_{n-1})\times X_n$ と同相である。
帰納法の仮定と二因子の場合からこれはコンパクトである。
:::

:::remark[Heine–Borel の定理の別証明]
$\R^n$ の Euclid 距離と最大値ノルム $\norm x_\infty=\max_i|x_i|$ の距離は同値であり（:ref[thm-top-norm-equivalence]）、後者の開球は開区間の積である。
従って $\R^n$ の位相は $\R$ の $n$ 個の積の積位相と一致する。
:ref[thm-top-interval-compact] と :ref[thm-top-finite-product-compact] により $[-M,M]^n$ はコンパクトであり、Heine–Borel の定理の十分性の証明の前半を二分法なしに得る。
:::

## 点列コンパクト性

[解析学](/math/analysis/) の Bolzano–Weierstrass の定理は、$\R$ の有界な点列が収束部分列をもつことを主張する。
これは「コンパクト集合の中の点列は収束部分列をもつ」という形に一般化できそうである。
実際、距離空間ではコンパクト性はこの点列の条件と同値になる。
一方、一般の位相空間では両者は独立な条件であり、その違いを本節と次節で調べる。

:::definition[点列コンパクト性・極限点コンパクト性・全有界性]{#def-top-total-bounded}
1. 位相空間 $X$ の任意の点列が、$X$ の点に収束する部分列をもつとき、$X$ は:term[点列コンパクト]{en="sequentially compact" yomi="てんれつこんぱくと"}であるという。
2. 位相空間 $X$ の任意の無限部分集合が $X$ の中に集積点をもつとき、$X$ は:term[極限点コンパクト]{en="limit point compact" yomi="きょくげんてんこんぱくと"}であるという（Bolzano–Weierstrass の性質をもつともいう）。
3. 距離空間 $X$ が:term[全有界]{en="totally bounded" yomi="ぜんゆうかい"}であるとは、任意の $r>0$ に対し有限個の点 $a_1,\dots,a_m\in X$ があって $X=\bigcup_{j}B(a_j,r)$ となることである。このような $\{a_1,\dots,a_m\}$ を有限 $r$ ネットという。
4. 距離空間 $X$ のすべての Cauchy 列が $X$ の点に収束するとき、$X$ は完備であるという。完備性は [完備距離空間](/math/topology/complete-metric/) で詳しく扱う。
:::

全有界性と通常の有界性の関係は :ref[lem-top-total-bounded-sequence] の前で調べた。

:::proposition[コンパクト空間の極限点コンパクト性]{#prop-top-compact-limit-point}
コンパクト空間は極限点コンパクトである。
:::

:::proof
$X$ をコンパクトとし、集積点をもたない無限部分集合 $A\subset X$ があると仮定する。
各 $x\in X$ は $A$ の集積点でないので、$x$ の開近傍 $U_x$ で $U_x\cap(A\setminus\{x\})=\varnothing$、すなわち $U_x\cap A\subset\{x\}$ となるものがある。
$(U_x)_{x\in X}$ は開被覆なので、有限個 $x_1,\dots,x_n$ で $X=\bigcup_jU_{x_j}$ となる。
すると $A=\bigcup_j(U_{x_j}\cap A)\subset\{x_1,\dots,x_n\}$ となり、$A$ が無限であることに反する。
:::

### 全有界性と Lebesgue 数

全有界性はすべての点列に Cauchy 部分列があることと同値である（:ref[lem-top-total-bounded-sequence]）。

次の補題の数 $\lambda$ は、被覆のどの元に入るかを「十分小さい集合」に対して一様に保証する。

:::lemma[Lebesgue 数の補題]{#lem-top-lebesgue-number}
点列コンパクトな距離空間 $X$ の任意の開被覆 $\mathcal U$ に対し、ある $\lambda>0$ が存在して、各 $x\in X$ について球 $B(x,\lambda)$ が $\mathcal U$ のある一つの元に含まれる。
特に直径が $\lambda$ 未満の空でない部分集合は $\mathcal U$ のある一つの元に含まれる。
このような $\lambda$ を被覆 $\mathcal U$ の :term[Lebesgue 数]{en="Lebesgue number" yomi="るべーぐすう"} という。
:::

:::proof
$X=\varnothing$ なら任意の $\lambda>0$ でよいので、$X\ne\varnothing$ とする。
結論が成り立たないと仮定すると、各 $n\ge1$ に対し $\lambda=1/n$ が条件を満たさず、ある点 $x_n$ の球 $B(x_n,1/n)$ は $\mathcal U$ のどの元にも含まれない。
点列コンパクト性から部分列 $x_{n_k}$ がある点 $x$ に収束する。
$x$ を含む $U\in\mathcal U$ を取り、$B(x,r)\subset U$ となる $r>0$ を選ぶ。
$k$ を十分大きく取れば $d(x_{n_k},x)<r/2$ かつ $1/n_k<r/2$ となる。
このとき $y\in B(x_{n_k},1/n_k)$ なら $d(y,x)\le d(y,x_{n_k})+d(x_{n_k},x)<r$ なので $B(x_{n_k},1/n_k)\subset B(x,r)\subset U$ となり、$x_{n_k}$ の選び方に反する。

後半について、$A$ を直径 $\operatorname{diam}A<\lambda$ の空でない集合とし、$a\in A$ を取る。
任意の $y\in A$ について $d(a,y)\le\operatorname{diam}A<\lambda$ なので $A\subset B(a,\lambda)$ であり、前半から $A$ は $\mathcal U$ の一つの元に含まれる。
:::

### 距離空間のコンパクト性

:::theorem[距離空間のコンパクト性の同値条件]{#thm-top-metric-compactness}
距離空間 $X$ について、次の四条件は同値である。

1. $X$ はコンパクトである。
2. $X$ は極限点コンパクトである。
3. $X$ は点列コンパクトである。
4. $X$ は完備かつ全有界である。
:::

:::proof
**1 ⇒ 2。** :ref[prop-top-compact-limit-point] による。

**2 ⇒ 3。** $(x_n)$ を点列とし、値の集合 $S=\{x_n:n\ge1\}$ を考える。
$S$ が有限なら、ある値 $a$ について $x_n=a$ となる $n$ が無限個あり、それらの添字を並べた部分列は $a$ に収束する定数列である。
$S$ が無限なら、仮定から $S$ の集積点 $x$ がある。
各 $\eps>0$ について $B(x,\eps)$ は $S$ の点を無限個含む。
実際、$B(x,\eps)\cap(S\setminus\{x\})$ が有限集合 $\{s_1,\dots,s_m\}$ なら、$\eps'=\min(\eps,d(x,s_1),\dots,d(x,s_m))>0$ について $B(x,\eps')\cap(S\setminus\{x\})=\varnothing$ となり、$x$ が集積点であることに反する（$m=0$ なら $\eps'=\eps$ とする）。
従って $x_n\in B(x,\eps)$ となる添字 $n$ も無限個ある。
そこで $n_1<n_2<\cdots$ を $x_{n_k}\in B(x,1/k)$ となるように順に選べば、$x_{n_k}\to x$ である。

**3 ⇒ 4。** 収束列は Cauchy 列なので、任意の点列は Cauchy 部分列をもち、:ref[lem-top-total-bounded-sequence] により $X$ は全有界である。
Cauchy 列の収束部分列が得られれば、:ref[lem-top-cauchy-subsequence] により全列も同じ点に収束する。従って $X$ は完備である。

**4 ⇒ 3。** :ref[lem-top-total-bounded-sequence] により任意の点列は Cauchy 部分列をもち、完備性によりそれは収束する。

**3 ⇒ 1。** $\mathcal U$ を開被覆とし、:ref[lem-top-lebesgue-number] によりその Lebesgue 数 $\lambda>0$ を取る。
3 ⇒ 4 で示したように $X$ は全有界なので、有限 $\lambda$ ネット $a_1,\dots,a_m$ がある。
各 $B(a_j,\lambda)$ は $\mathcal U$ のある元 $U_j$ に含まれ、$X=\bigcup_jB(a_j,\lambda)\subset\bigcup_jU_j$ である。
よって $U_1,\dots,U_m$ が有限部分被覆である。
:::

:::corollary[コンパクト距離空間の可算性]{#cor-top-compact-metric-separable}
コンパクト距離空間は有界であり、可分かつ第二可算である。
:::

:::proof
:ref[thm-top-metric-compactness] により $X$ は全有界なので、上で述べたように有界である。
各 $k\ge1$ について有限 $1/k$ ネット $A_k$ を取り、$D=\bigcup_kA_k$ と置くと $D$ は可算である。
任意の $x\in X$ と $r>0$ について、$1/k<r$ となる $k$ を取れば $x\in B(a,1/k)$ となる $a\in A_k$ があり、$a\in B(x,r)$ である。
よって $D$ は稠密であり、$X$ は可分である。
距離空間では可分性と第二可算性は同値である（:ref[prop-top-separable-metric]）。
:::

:::corollary[Heine–Cantor の定理]{#cor-top-heine-cantor}
コンパクト距離空間 $(X,d)$ から距離空間 $(Y,\rho)$ への連続写像 $f$ は一様連続である。
:::

:::proof
$\eps>0$ を固定する。
各 $y\in Y$ について開球 $B(y,\eps/2)$ の逆像 $f^{-1}(B(y,\eps/2))$ は開であり、それらの全体は $X$ の開被覆である。
:ref[thm-top-metric-compactness] により $X$ は点列コンパクトなので、この被覆の Lebesgue 数 $\lambda>0$ がある。
$d(x,x')<\lambda$ なら二点集合 $\{x,x'\}$ の直径は $\lambda$ 未満なので、ある $y$ について $x,x'\in f^{-1}(B(y,\eps/2))$ となる。
従って $\rho(f(x),f(x'))\le\rho(f(x),y)+\rho(y,f(x'))<\eps$ である。
$\lambda$ は $x,x'$ によらないので $f$ は一様連続である。
:::

:::example[全有界性と完備性はそれぞれ必要である]
1. $(0,1)$ は全有界である。$r>0$ に対し $1/m<r$ となる $m$ を取れば、点 $j/m$（$1\le j\le m-1$）が有限 $r$ ネットになる。しかし Cauchy 列 $1/(n+1)$ は $(0,1)$ に極限をもたないので完備でなく、コンパクトでない。
2. 無限離散距離空間は完備である。Cauchy 条件を $\eps=1$ で使うと、ある番号以後の項はすべて等しいからである。しかし上で見たように全有界でない。
3. $\R$ も完備だが全有界でない。$\R^n$ の部分集合では、全有界性は有界性と同値なので、定理の条件 4 は「閉かつ有界」と同値になる（完備な $\R^n$ の部分集合が完備であることと閉であることは同値である。[完備距離空間](/math/topology/complete-metric/) を参照）。こうして Heine–Borel の定理の三つ目の証明が得られる。
:::

:::example[Heine–Cantor の定理の仮定]
$f(x)=x^2$ は $\R$ 上で一様連続でない。$x_n=n$、$y_n=n+1/n$ は $|x_n-y_n|\to0$ を満たすが $|f(y_n)-f(x_n)|=2+1/n^2\ge2$ である（:ref[prop-top-uniform-pairs]）。
$(0,1]$ 上の $\sin(1/x)$ も一様連続でない。$x_n=1/(2\pi n+\pi/2)$、$y_n=1/(2\pi n)$ は距離が $0$ に近づくが、値は $1$ と $0$ である。
いずれも定義域がコンパクトでない。
:::

:::example[Lebesgue 数を具体的に求める]
$[0,1]$ の相対開集合 $U=[0,2/3)$、$V=(1/3,1]$ からなる被覆を考える。
直径が $1/3$ 未満の空でない集合 $A$ は $U$ か $V$ に含まれる。
$A\not\subset U$ なら $A$ は $2/3$ 以上の点 $a$ を含み、$A$ の任意の点 $y$ は $y>a-1/3\ge1/3$ を満たすので $A\subset V$ だからである。
直径の上界 $1/3$ はこれ以上大きくできない。
二点集合 $\{1/3,2/3\}$ の直径はちょうど $1/3$ で、$2/3\notin U$、$1/3\notin V$ なのでどちらにも含まれない。

一方、球の形の条件「すべての $x$ で $B(x,r)$ が一方に含まれる」では最良の値は $r=1/6$ である。
中心 $x\le1/2$ なら $B(x,1/6)\subset[0,x+1/6)\subset[0,2/3)$、$x\ge1/2$ なら $B(x,1/6)\subset(x-1/6,1]\subset(1/3,1]$ である。
$r>1/6$ なら $B(1/2,r)$ は $1/3$ と $2/3$ をともに含み、どちらにも含まれない。
このように直径による定式化と球による定式化では最良の数が異なるが、どちらの Lebesgue 数も存在することが重要である。
:::

```tikz caption="区間を覆う二つの相対開集合。白丸は除いた端点で、両端点の距離が直径版の最良の Lebesgue 数を与える。"
\begin{tikzpicture}[x=6cm]
\draw[gray] (0,-.4) -- (1,-.4);
\foreach \x/\t in {0/0,.3333/{\frac{1}{3}},.6667/{\frac{2}{3}},1/1}{
  \draw (\x,-.36) -- (\x,-.44) node[below] {$\t$};
}
\draw[very thick] (0,.65) -- (.6667,.65);
\fill (0,.65) circle (1.6pt);
\draw[fill=white] (.6667,.65) circle (1.9pt);
\draw[very thick] (.3333,.12) -- (1,.12);
\draw[fill=white] (.3333,.12) circle (1.9pt);
\fill (1,.12) circle (1.6pt);
\node[left] at (0,.65) {$U$};
\node[right] at (1,.12) {$V$};
\end{tikzpicture}
```

### 一般の位相空間での比較

距離空間の外では、定理の四条件の関係は崩れる。
まず極限点コンパクト性からコンパクト性が従わない例を見る。

:::example[極限点コンパクトだがコンパクトでない空間]
$Y=\{a,b\}$ に密着位相を入れ、$X=\N\times Y$ に積位相を入れる。
$X$ の開集合は $A\times Y$（$A\subset\N$）の形である。
空でない任意の部分集合 $S\subset X$ は集積点をもつ。
$(n,a)\in S$ なら $(n,b)$ の任意の開近傍は $\{n\}\times Y$ を含むので、$(n,b)$ と異なる $S$ の点 $(n,a)$ を含む。従って $(n,b)$ は $S$ の集積点である（$(n,b)\in S$ の場合は $(n,a)$ が集積点になる）。
よって $X$ は極限点コンパクトである。
一方、開被覆 $\{\{n\}\times Y:n\in\N\}$ は有限部分被覆をもたないので $X$ はコンパクトでない。
点列 $(n,a)$ はどの点にも収束する部分列をもたないので、点列コンパクトでもない。
この空間は $T_1$ でない（一点集合 $\{(n,a)\}$ は閉でない）ことに注意する。
:::

次に点列の条件と被覆の条件の中間にある条件を導入する。

:::definition[可算コンパクト性]{#def-top-countably-compact}
位相空間 $X$ の任意の可算開被覆が有限部分被覆をもつとき、$X$ は:term[可算コンパクト]{en="countably compact" yomi="かさんこんぱくと"}であるという。
:::

:::proposition[可算コンパクト性の同値条件]{#prop-top-countably-compact}
位相空間 $X$ について次は同値である。

1. $X$ は可算コンパクトである。
2. 空でない閉集合の任意の減少列 $F_1\supset F_2\supset\cdots$ について $\bigcap_nF_n\ne\varnothing$ である。
3. 任意の点列 $(x_n)$ が集積点をもつ。ここで $x$ が点列の集積点であるとは、$x$ の任意の近傍 $U$ について $x_n\in U$ となる $n$ が無限個あることをいう。

さらに $X$ が第一可算なら、これらは点列コンパクト性とも同値である。
またコンパクト空間と点列コンパクト空間は常に可算コンパクトである。
:::

:::proof
**1 ⇔ 2。** 減少列 $F_n$ の共通部分が空なら、補集合は増加する可算開被覆になる。有限個で覆えれば最大の添字 $N$ で $F_N=\varnothing$ となる。逆に、有限個では覆えない可算開被覆 $(U_n)$ からは、$F_n=X\setminus\bigcup_{j\le n}U_j$ として共通点のない空でない閉集合の減少列を得る。

**2 ⇔ 3。** 点列の集積点は、すべての末尾の閉包 $F_n=\overline{\{x_m:m\ge n\}}$ に属する点にほかならない。実際、どの近傍にも任意の $n$ 以後の項があることが、まさにこの閉包の条件である。従って 2 から 3 が従う。逆に、空でない閉集合の減少列 $F_n$ から $x_n\in F_n$ を選ぶ。固定した $k$ に対して閉集合 $F_k$ は第 $k$ 項以後を含むので、その点列の集積点は $F_k$ に入る。これはすべての $k$ に成り立つ。

**第一可算な場合。** 集積点 $x$ の可算基本近傍系 $V_j$ を取り、$W_j=\bigcap_{i\le j}V_i$ とする。各 $W_j$ に無限個の項が入るので、$x_{n_j}\in W_j$、$n_j>n_{j-1}$ と順に選べる。任意の近傍がある $W_j$ を含むため $x_{n_j}\to x$ である。逆に収束部分列の極限は常に集積点である。コンパクト空間では可算開被覆も有限部分被覆をもつので、最後の主張も従う。
:::

:::remark[極限点コンパクト性との関係]
可算コンパクト空間は極限点コンパクトである。
無限集合 $A$ から相異なる点の列 $a_1,a_2,\dots$ を取り、その集積点 $x$ を取ると、$x$ の任意の近傍は無限個の $a_n$ を含むので $x$ と異なる $A$ の点を含むからである。
逆は $T_1$ 空間では成り立つ。
点列 $(x_n)$ の値の集合 $S$ が有限なら無限回現れる値が集積点である。
$S$ が無限ならその集積点 $x$ を取る。$x$ のある開近傍 $U$ に入る $x_n$ が有限個しかなければ、それらの値のうち $x$ と異なるものの有限集合 $E$ は $T_1$ 性から閉であり、$U\setminus E$ は $S\setminus\{x\}$ と交わらない $x$ の開近傍となって矛盾する。
上の $\N\times\{a,b\}$ の例は、$T_1$ 性を外すとこの逆が成り立たないことを示している。
:::

:::example[点列コンパクトだがコンパクトでない空間]{#exm-top-omega-one}
最小の非可算順序数を $\omega_1$ とし、$\omega_1$ より小さい順序数全体 $[0,\omega_1)$ に順序位相を入れる。
順序数と超限帰納法については [整列集合と順序数](/math/sets-logic/ordinals/) を参照する。
各 $\beta<\omega_1$ は可算集合であり、可算個の可算順序数の上限は可算順序数である。

**各 $[0,\alpha]$ はコンパクトである。** $\alpha<\omega_1$ と開被覆を固定する。$[0,\beta]$ が有限個で覆えない最小の $\beta\le\alpha$ があると仮定する。$0$ は一つの被覆要素で覆える。後続順序数なら直前までの有限被覆に $\beta$ を含む元を加えればよい。極限順序数なら $\beta$ を含む被覆要素はある $\gamma<\beta$ について $(\gamma,\beta]$ を含むので、$[0,\gamma]$ の有限被覆にその元を加えればよい。どの場合も矛盾する。

**点列コンパクトである。** 任意の点列 $(\beta_n)$ は、可算順序数 $\alpha=\sup_n\beta_n$ までのコンパクト区間に入るので、:ref[prop-top-countably-compact] によりその中に集積点をもつ。各 $\beta>0$ の基本近傍 $(\gamma,\beta]$（$\gamma<\beta$）は可算個で、$0$ では $\{0\}$ が開である。従って第一可算であり、同じ命題により集積点に収束する部分列を選べる。

**コンパクトでない。** $\{[0,\beta):\beta<\omega_1\}$ は開被覆である（$\alpha\in[0,\alpha+1)$）。
有限個 $[0,\beta_1),\dots,[0,\beta_k)$ の和は $[0,\max_j\beta_j)$ であり、点 $\max_j\beta_j$ を覆わない。
:::

本節の関係をまとめておく。
コンパクト ⇒ 可算コンパクト ⇐ 点列コンパクトであり、可算コンパクト ⇒ 極限点コンパクトである。
距離空間ではこれらはすべて同値である。
$[0,\omega_1)$ は点列コンパクトだがコンパクトでなく、次節の $\{0,1\}^{\mathcal P(\N)}$ はコンパクトだが点列コンパクトでない。

## Tychonoff の定理

有限個のコンパクト空間の積がコンパクトであることは :ref[thm-top-finite-product-compact] で示した。
可算個の積では、点列の対角線論法が使える。

:::proposition[可算積の点列コンパクト性]{#prop-top-countable-product-sequential}
点列コンパクト空間の列 $X_1,X_2,\dots$ の積 $\prod_{m\ge1}X_m$ は点列コンパクトである。
特に可算個のコンパクト距離空間の積はコンパクト距離化可能であり、Hilbert 立方体 $[0,1]^{\N}$ はコンパクトである。
:::

:::proof
$(x^{(k)})_{k\ge1}$ を積の点列とする。
$N_0=\{1,2,\dots\}$ と置き、$m\ge1$ について無限集合 $N_m\subset N_{m-1}$ と点 $p_m\in X_m$ を、$k\in N_m$ を増加順に並べた列 $x^{(k)}_m$ が $p_m$ に収束するように選ぶ。
これは $N_{m-1}$ に沿った第 $m$ 座標の点列に $X_m$ の点列コンパクト性を使えば可能である。
次に $k_1\in N_1$ を選び、$k_m\in N_m$ を $k_m>k_{m-1}$ となるように選ぶ。
各 $m$ について $j\ge m$ なら $k_j\in N_j\subset N_m$ なので、$(x^{(k_j)}_m)_{j\ge m}$ は $N_m$ に沿った収束列の部分列であり $p_m$ に収束する。
従って $x^{(k_j)}$ は各座標で $p=(p_m)$ に収束し、:ref[thm-top-product-property] により積位相で $p$ に収束する。

各 $X_m$ がコンパクト距離空間なら、:ref[thm-top-metric-compactness] により点列コンパクトであり、積は :ref[prop-top-product-metric] の距離で距離化される。
積は点列コンパクトな距離空間なので、再び :ref[thm-top-metric-compactness] によりコンパクトである。
$[0,1]$ はコンパクト距離空間なので $[0,1]^{\N}$ はコンパクトである。
:::

非可算個の積ではこの方法は使えない。
実際、後の例で見るようにコンパクトな積が点列コンパクトでないことがある。
そこで点列の代わりに、有限個の条件を同時に扱う「フィルター」を用いる。

### 超フィルターの極限

フィルターの収束には :ref[def-top-ultrafilter] を用いる。有限交叉性をもつ族から超フィルターを作る手順と補集合に関する二択性は :ref[lem-top-ultrafilter-extension] で証明した。

:::lemma[コンパクト性と超フィルターの収束]{#lem-top-ultrafilter-limit}
空でない位相空間 $X$ がコンパクトであるための必要十分条件は、$X$ 上のすべての超フィルターが少なくとも一点に収束することである。
:::

:::proof
$X$ をコンパクトとし、$\mathscr U$ を超フィルターとする。
閉集合族 $\{\overline F:F\in\mathscr U\}$ は有限交叉性をもつ。
$F_1,\dots,F_n\in\mathscr U$ なら $F_1\cap\dots\cap F_n\in\mathscr U$ は空でなく、$\overline{F_1}\cap\dots\cap\overline{F_n}$ に含まれるからである。
:ref[thm-top-compact-fip] により共通点 $x$ がある。
$x$ の開近傍 $V$ が $\mathscr U$ に属さないと仮定すると、:ref[lem-top-ultrafilter-extension] の 2 により $X\setminus V\in\mathscr U$ である。
$X\setminus V$ は閉なので $x\in\overline{X\setminus V}=X\setminus V$ となり、$x\in V$ に反する。
従って $x$ のすべての開近傍が $\mathscr U$ に属し、上方への包含からすべての近傍が属する。
すなわち $\mathscr U$ は $x$ に収束する。

逆に、すべての超フィルターが収束すると仮定し、$\mathcal F$ を有限交叉性をもつ閉集合の空でない族とする。
:ref[lem-top-ultrafilter-extension] の 1 により $\mathcal F$ を含む超フィルター $\mathscr U$ があり、その極限の一つを $x$ とする。
ある $F\in\mathcal F$ が $x$ を含まなければ、$X\setminus F$ は $x$ の開近傍なので $\mathscr U$ に属し、$F\in\mathscr U$ と合わせて $\varnothing\in\mathscr U$ となり矛盾する。
よって $x\in\bigcap\mathcal F$ であり、:ref[thm-top-compact-fip] により $X$ はコンパクトである。
:::

Hausdorff 空間ではフィルターの極限は高々一つである。
$\mathscr F$ が異なる二点 $x,y$ に収束すれば、それらを分離する開集合 $U\ni x$、$V\ni y$ がともに $\mathscr F$ に属し、$\varnothing=U\cap V\in\mathscr F$ となるからである。
非 Hausdorff 空間では極限が複数ありうる。たとえば密着空間では任意のフィルターがすべての点に収束する。

### Tychonoff の定理の証明

:::theorem[Tychonoff の定理]{#thm-top-tychonoff}
任意の集合 $I$ で添字付けられたコンパクト空間の族 $(X_i)_{i\in I}$ について、積位相を入れた $P=\prod_{i\in I}X_i$ はコンパクトである。
因子の Hausdorff 性や添字集合の可算性は仮定しない。
:::

:::proof
ある $X_i$ が空なら $P$ は空でありコンパクトである。
$I=\varnothing$ なら $P$ は一点でありコンパクトである。
$P$ が空の場合は既に結論が成り立つので、以下 $P\ne\varnothing$、$I\ne\varnothing$ とする。
:ref[lem-top-ultrafilter-limit] により、$P$ 上の任意の超フィルター $\mathscr U$ が収束することを示せばよい。

**射影した超フィルター。** 各 $i\in I$ について

$$
\mathscr U_i=\{A\subset X_i:\pi_i^{-1}(A)\in\mathscr U\}
$$

と置く。逆像は全体集合、空集合、有限交叉および包含を保つので、$\mathscr U_i$ は真のフィルターである。また $\pi_i^{-1}(X_i\setminus A)=P\setminus\pi_i^{-1}(A)$ だから、$\mathscr U_i$ は各 $A$ とその補集合のちょうど一方を含む。:ref[lem-top-ultrafilter-extension] によりこれは超フィルターである。

**極限の選択。** $X_i$ はコンパクトなので、:ref[lem-top-ultrafilter-limit] により $\mathscr U_i$ の極限全体 $L_i\subset X_i$ は空でない。
選択公理により各 $i$ について $x_i\in L_i$ を選び、$x=(x_i)_{i\in I}\in P$ と置く。

**$\mathscr U$ は $x$ に収束する。** $W$ を $x$ の近傍とする。
積位相の定義（:ref[def-top-product]）から、有限集合 $J\subset I$ と開集合 $V_i\ni x_i$（$i\in J$）で

$$
x\in\bigcap_{i\in J}\pi_i^{-1}(V_i)\subset W
$$

となるものがある。
各 $i\in J$ について $V_i$ は $x_i$ の近傍であり、$\mathscr U_i$ は $x_i$ に収束するので $V_i\in\mathscr U_i$、すなわち $\pi_i^{-1}(V_i)\in\mathscr U$ である。
フィルターは有限個の共通部分で閉じているので $\bigcap_{i\in J}\pi_i^{-1}(V_i)\in\mathscr U$ であり（$J=\varnothing$ なら共通部分は $P\in\mathscr U$）、上方への包含により $W\in\mathscr U$ である。
従って $\mathscr U$ は $x$ に収束し、$P$ はコンパクトである。
:::

:::remark[証明の要点と選択公理]
証明が有効に働いたのは、積位相の基本近傍が有限個の座標しか制限しないこと、そしてフィルターが有限個の共通部分で閉じていることが噛み合うからである。
選択公理は二か所で用いた。超フィルターの存在（Zorn の補題）と、各座標の極限 $x_i$ の選択である。
因子が Hausdorff なら極限は一意なので、後者の選択は不要になる。
:::

:::example[Tychonoff の定理から得られるコンパクト空間]
1. 二点離散空間 $\{0,1\}$ の任意個の積 $\{0,1\}^I$（Cantor 立方体）はコンパクトである。
2. $[0,1]^{[0,1]}$、すなわち $[0,1]$ から $[0,1]$ への写像全体に各点収束の位相（積位相）を入れた空間はコンパクトである。
3. 実数列 $(a_n)$ で $|a_n|\le c_n$（$c_n\ge0$ は固定）を満たすもの全体は、$\prod_n[-c_n,c_n]$ として各点収束の位相でコンパクトである。
:::

:::example[コンパクトだが点列コンパクトでない積]
$P=\{0,1\}^{\mathcal P(\N)}$、すなわち $\N$ の部分集合全体を添字集合とする二点離散空間の積を考える。
Tychonoff の定理により $P$ はコンパクトである。
点列 $x_n\in P$ を、座標 $A\subset\N$ で $x_n(A)=1$（$n\in A$ のとき）、$x_n(A)=0$（$n\notin A$ のとき）と定める。
任意の部分列 $x_{n_1},x_{n_2},\dots$（$n_1<n_2<\cdots$）について、座標 $A=\{n_2,n_4,n_6,\dots\}$ を見ると、$x_{n_k}(A)$ は $k$ が奇数なら $0$、偶数なら $1$ である。
この座標で収束しないので、:ref[thm-top-product-property] により部分列は積位相で収束しない。
よって $P$ は点列コンパクトでない。
添字が非可算なので、可算積で使った対角線論法は使えない。
:::

:::example[箱位相ではコンパクト性が失われる]
$\{0,1\}^{\N}$ に箱位相（各座標の開集合の任意個の積を開基とする位相）を入れると、各一点 $\{x\}=\prod_n\{x_n\}$ が開になり、離散空間になる。
この空間は無限集合（実際は非可算集合）なので、:ref[exm-top-compact-basic] の 2 によりコンパクトでない。
Tychonoff の定理の証明では、基本近傍が有限個の座標しか制限しないことを使ったので、箱位相では成り立たない。
:::

### Alexander の準開基定理

:::theorem[Alexander の準開基定理]{#thm-top-alexander}
$\mathcal S$ を位相空間 $X$ の位相の準開基とする。
$\mathcal S$ の元だけからなる $X$ の開被覆がすべて有限部分被覆をもつなら、$X$ はコンパクトである。
:::

:::proof
$X=\varnothing$ ならコンパクトなので、$X\ne\varnothing$ とする。
$X$ がコンパクトでないと仮定すると、:ref[lem-top-ultrafilter-limit] により、どの点にも収束しない超フィルター $\mathscr U$ がある。

各点 $x\in X$ について、$x\in S$ かつ $S\notin\mathscr U$ となる $S\in\mathcal S$ が存在することを示す。
そうでなければ、$x$ を含む $\mathcal S$ の元はすべて $\mathscr U$ に属する。
準開基の定義から、$x$ の任意の近傍はある有限個の $S_1,\dots,S_n\in\mathcal S$（$x\in S_j$）の共通部分を含み、その共通部分は $\mathscr U$ に属するので、近傍も $\mathscr U$ に属する。
すると $\mathscr U$ は $x$ に収束してしまい、仮定に反する。

従って $\mathcal C=\{S\in\mathcal S:S\notin\mathscr U\}$ は $X$ を覆い、仮定から有限個 $S_1,\dots,S_n\in\mathcal C$ で $X=S_1\cup\dots\cup S_n$ となる（$X\ne\varnothing$ なので $n\ge1$）。
各 $S_j\notin\mathscr U$ なので、:ref[lem-top-ultrafilter-extension] の 2 により $X\setminus S_j\in\mathscr U$ である。
その共通部分 $\bigcap_j(X\setminus S_j)=X\setminus\bigcup_jS_j=\varnothing$ も $\mathscr U$ に属することになり、矛盾する。
:::

Alexander の定理を使うと、Tychonoff の定理の別証明が得られる。
積位相の準開基 $\mathcal S=\{\pi_i^{-1}(U):i\in I,\ U\text{ は }X_i\text{ で開}\}$ の元からなる被覆 $\mathcal C$ で、有限部分被覆をもたないものがあるとする（$P\ne\varnothing$ の場合を考えればよい）。
各 $i$ について $\mathcal C_i=\{U:\pi_i^{-1}(U)\in\mathcal C\}$ と置く。
もし $\mathcal C_i$ が $X_i$ を覆えば、$X_i$ のコンパクト性から有限個 $U_1,\dots,U_n\in\mathcal C_i$ が $X_i$ を覆い、$\pi_i^{-1}(U_1),\dots,\pi_i^{-1}(U_n)$ が $P$ を覆ってしまう。
従って各 $i$ について $X_i\setminus\bigcup\mathcal C_i$ は空でなく、選択公理によりその点 $x_i$ を選べる。
点 $x=(x_i)$ はある $\pi_i^{-1}(U)\in\mathcal C$ に属するが、それは $x_i\in U\in\mathcal C_i$ を意味し、$x_i$ の選び方に反する。

:::example[Tychonoff の定理は選択公理を含む]
選択公理を仮定せずに Tychonoff の定理を仮定すると、選択公理が導かれる。
非空集合の族 $(A_i)_{i\in I}$ が与えられたとする。
新しい点 $*_i$ を加えた $Y_i=A_i\amalg\{*_i\}$ に、開集合が $\varnothing$、$\{*_i\}$、$Y_i$ の三つだけの位相を入れる。
各 $Y_i$ は開集合が有限個なのでコンパクトであり、積 $\prod_iY_i$ には点 $(*_i)_i$ があるので空でない（ここでは選択を使わない）。
$\{*_i\}$ は開なので $A_i$ は $Y_i$ で閉であり、$C_i=\pi_i^{-1}(A_i)$ は積の閉集合である。
族 $\{C_i\}$ は有限交叉性をもつ。
有限個の $i_1,\dots,i_n$ について $a_{i_j}\in A_{i_j}$ を（有限個なので選択公理なしに）選び、残りの座標に $*_i$ を入れた点が $C_{i_1}\cap\dots\cap C_{i_n}$ に属するからである。
Tychonoff の定理から積はコンパクトなので、:ref[thm-top-compact-fip] により $\bigcap_iC_i$ に点 $a$ がある。
$a_i\in A_i$ がすべての $i$ で成り立つので、$a$ は選択関数である。
選択公理と Zorn の補題の同値性は [選択公理と Zorn の補題](/math/sets-logic/choice/) で扱う。
この事実を含む Tychonoff の定理の位置付けについては :cite[uchida-topology] も参照するとよい。
:::

### 無限積でのコンパクト性の応用

:::example[有限個ずつ解ける二値制約を同時に解く]
集合 $I$ の各元 $i$ に $0$ または $1$ を割り当てる問題を考える。
一つの制約は有限集合 $J\subset I$ の座標だけを参照し、許される値の組の集合 $E\subset\{0,1\}^J$ を指定するものとする。
制約の族が与えられ、そのうちのどの有限個も同時に満たす割り当てがあるならば、すべての制約を同時に満たす割り当てが存在する。

実際、割り当て全体 $P=\{0,1\}^I$ は Tychonoff の定理によりコンパクトである。
制約 $(J,E)$ を満たす割り当ての集合は、射影 $p_J:P\to\{0,1\}^J$、$x\mapsto(x_i)_{i\in J}$ による $E$ の逆像である。
$\{0,1\}^J$ は有限個の離散空間の積なので離散であり、$E$ は閉、$p_J$ は連続なので、この逆像は閉である。
仮定によりこれらの閉集合の族は有限交叉性をもち、:ref[thm-top-compact-fip] により共通点がある。
それが求める割り当てである。

例えばグラフの頂点全体を $I$ とし、各辺 $\{u,v\}$ について「$u$ と $v$ の値が異なる」という制約を置く。
有限部分グラフがすべて二色で塗り分けられるなら、有限個の辺に関する制約は、それらの辺の端点を塗り分け、残りの頂点に $0$ を割り当てれば同時に満たされる。
従ってグラフ全体も二色で塗り分けられる。
ここでは各制約が有限個の座標しか見ないことが閉性を保証した。
無限個の座標について「どこか一つは $1$」と要求する集合は開であるが、$I$ が無限なら閉でない。
すべての座標が $0$ の点のどの基本近傍も、制限されていない座標の一つだけを $1$ にした点を含むからである。
:::

## 局所コンパクト空間と一点コンパクト化

$\R^n$ はコンパクトでないが、各点の周りに閉球というコンパクトな近傍をもつ。
このように局所的にはコンパクト性が使える空間は、解析学・幾何学で最も頻繁に現れる。
本節ではそのような空間の基本的性質を調べ、一点を付け加えてコンパクト空間にする方法を述べる。

:::definition[局所コンパクト性]{#def-top-locally-compact}
位相空間 $X$ の各点がコンパクトな近傍をもつとき、すなわち各 $x\in X$ に対しコンパクト集合 $K$ と開集合 $U$ で $x\in U\subset K$ となるものがあるとき、$X$ は:term[局所コンパクト]{en="locally compact" yomi="きょくしょこんぱくと"}であるという。
局所コンパクトな Hausdorff 空間を局所コンパクト Hausdorff 空間という。
:::

コンパクトな近傍そのものが開である必要はない。

:::example[局所コンパクトな空間とそうでない空間]{#exm-top-locally-compact}
1. コンパクト空間は局所コンパクトである。$X$ 自身が各点のコンパクトな近傍である。
2. $\R^n$ は局所コンパクト Hausdorff である。各点 $x$ の閉球 $\{y:\norm{y-x}\le1\}$ は Heine–Borel の定理によりコンパクトで、開球 $B(x,1)$ を含む。$\R^n$ の開集合 $U$ も、$x\in U$ に対し $U$ に含まれる閉球を取れば局所コンパクトである。
3. 離散空間は局所コンパクトである。$\{x\}$ がコンパクトな開近傍である。
4. $\Q$ は局所コンパクトでない。ある $q\in\Q$ がコンパクトな近傍 $K$ をもつとすると、$K$ はある $\Q\cap(q-\eps,q+\eps)$ を含む。$K$ は $\R$ でもコンパクト（包含写像による像）なので $\R$ で閉であり、$\Q\cap(q-\eps,q+\eps)$ の $\R$ での閉包 $[q-\eps,q+\eps]$ を含む。これは無理数を含み $K\subset\Q$ に反する。
5. 可算積 $\R^{\N}$（積位相）は局所コンパクトでない。点 $0$ がコンパクトな近傍 $K$ をもつとすると、$K$ は基本開集合 $B=\bigcap_{i\in J}\pi_i^{-1}(V_i)\ni0$（$J$ は有限）を含む。$i\notin J$ を一つ取ると、$B$ は第 $i$ 座標を制限しないので $\pi_i(B)=\R$ であり、$\pi_i(K)\supset\R$ となる。ところが $\pi_i(K)$ はコンパクト集合の連続像なのでコンパクトであり、$\pi_i(K)=\R$ がコンパクトでないことに反する。
:::

:::lemma[コンパクトな閉包をもつ開近傍]{#lem-top-compact-neighborhood}
$X$ を局所コンパクト Hausdorff 空間とし、$x\in U$、$U$ を開集合とする。
このとき $x\in V\subset\overline V\subset U$ を満たし、$\overline V$ がコンパクトである開集合 $V$ が存在する。
:::

:::proof
$x$ のコンパクトな近傍 $K$ と、$x\in O\subset K$ となる開集合 $O$ を取る。
Hausdorff 空間なので $K$ は閉であり（:ref[prop-top-compact-closed]）、$F=K\setminus U=K\cap(X\setminus U)$ はコンパクト空間 $K$ の閉部分集合なのでコンパクトである。
$x\notin F$ なので、:ref[lem-top-compact-separation] により開集合 $H\ni x$、$G\supset F$ で $H\cap G=\varnothing$ となるものがある（$F=\varnothing$ なら $H=X$、$G=\varnothing$ とする）。
$V=H\cap O\cap U$ と置くと、$V$ は $x$ を含む開集合である。
$V\subset O\subset K$ で $K$ は閉なので $\overline V\subset K$ である。
また $V\subset H\subset X\setminus G$ で $X\setminus G$ は閉なので $\overline V\subset X\setminus G$ である。
従って $\overline V\subset K\setminus G$ であり、$K\setminus U=F\subset G$ より $K\setminus G\subset U$ なので $\overline V\subset U$ である。
$\overline V$ はコンパクトな $K$ の閉部分集合なのでコンパクトである。
:::

:::corollary[コンパクト集合の周りのコンパクトな閉包をもつ開集合]{#cor-top-lch-compact-shrink}
局所コンパクト Hausdorff 空間でコンパクト集合 $K$ が開集合 $U$ に含まれるなら、$K\subset V\subset\overline V\subset U$ を満たし、$\overline V$ がコンパクトである開集合 $V$ が存在する。
:::

:::proof
各 $x\in K$ について :ref[lem-top-compact-neighborhood] で開集合 $V_x$（$x\in V_x\subset\overline{V_x}\subset U$、$\overline{V_x}$ はコンパクト）を取る。
有限個 $V_{x_1},\dots,V_{x_n}$ が $K$ を覆い、$V=\bigcup_jV_{x_j}$ と置けば $\overline V=\bigcup_j\overline{V_{x_j}}\subset U$ であり、:ref[prop-top-compact-closed] の 2 により $\overline V$ はコンパクトである（$K=\varnothing$ なら $V=\varnothing$）。
:::

### 一点コンパクト化

$\R$ に一点 $\infty$ を加えて、$|x|\to\infty$ を「$x\to\infty$」と読めるようにすると、円周が得られる。
同じことを一般の局所コンパクト Hausdorff 空間で行う。
$\infty$ の近傍は「コンパクト集合の外側」であるべきである。

:::theorem[一点コンパクト化]{#thm-top-one-point}
$X$ をコンパクトでない局所コンパクト Hausdorff 空間とし、$X$ に属さない一点 $\infty$ を加えた集合 $\alpha X=X\cup\{\infty\}$ を考える。
$\alpha X$ の部分集合で、次のいずれかの形のもの全体を $\mathcal T^*$ とする。

1. $X$ の開集合 $U$。
2. $\{\infty\}\cup(X\setminus K)$、ただし $K\subset X$ はコンパクト。

このとき次が成り立つ。

- $\mathcal T^*$ は $\alpha X$ の位相であり、$X$ は $\alpha X$ の開部分空間で、その相対位相は $X$ のもとの位相に一致する。
- $\alpha X$ はコンパクト Hausdorff 空間であり、$X$ は $\alpha X$ で稠密である。
- （一意性）$Y$ がコンパクト Hausdorff 空間で、$X\subset Y$ が部分空間（相対位相がもとの位相に一致する）であり、$Y\setminus X$ がただ一点 $p$ からなるならば、$X$ 上で恒等写像、$\infty\mapsto p$ で定まる全単射 $h:\alpha X\to Y$ は同相写像である。

$\alpha X$ を $X$ の:term[一点コンパクト化]{en="one-point compactification" yomi="いってんこんぱくとか"}（Alexandroff コンパクト化）という。
:::

:::proof
**位相と相対位相。** 第二の形の集合は補集合が $X$ のコンパクト集合であるような $\infty$ の近傍候補である。Hausdorff 性から各コンパクト集合 $K$ は $X$ で閉であり、どちらの形の集合も $X$ と交わすと開になる。空集合は第一の形、全体は $K=\varnothing$ とした第二の形である。二つの候補の交わりが $\infty$ を含まなければ $X$ の開集合であり、含めばその補集合は二つのコンパクト集合の和でコンパクトである。任意和 $W$ についても $W\cap X$ は開である。$\infty\in W$ のときは第二の形の一つを選べば $X\setminus W$ はその補集合 $K$ に含まれる。これは $X$ で閉なので $K$ の閉部分集合としてコンパクトである。従って $W$ は第二の形であり、位相の公理がそろう。第一の形には $X$ の全開集合が含まれるので、$X$ はもとの位相をもつ開部分空間である。

**コンパクト性と分離。** 任意の開被覆から $\infty$ を含む元を一つ選ぶ。その補集合 $K$ はコンパクトなので残りは有限個の被覆要素で覆える。$X$ 内の二点はもとの開集合で分離できる。$x\in X$ と $\infty$ を分離するには局所コンパクト性によりコンパクト近傍 $K$ を取り、$x\in U\subset K$ となる開集合 $U$ と $\alpha X\setminus K$ を用いればよい。さらに $X$ が非コンパクトならコンパクト集合 $K$ は $X$ を尽くさず、$\infty$ の各開近傍は $X$ に交わる。従って $X$ は稠密である。

**一意性。** $Y$ は Hausdorff なので $\{p\}$ は閉、従って $X$ は $Y$ で開である。$p$ を含まない開集合はちょうど $X$ の開集合である。$p$ を含む開集合の補集合は $Y$ の閉部分集合としてコンパクトであり、$X$ に含まれる。逆に $X$ のコンパクト集合は包含写像により $Y$ でもコンパクトであり、Hausdorff 性から $Y$ で閉である。こうして $p$ の開近傍も第二の形と正確に対応する。指定の全単射は両方向に連続である。
:::

:::remark[仮定の役割]
$X$ がすでにコンパクトでも、同じ $\mathcal T^*$ は位相になり、$\alpha X$ はコンパクト Hausdorff である。
しかし $K=X$ と取れるので $\{\infty\}$ が開になり、$\infty$ は孤立点で $X$ は稠密でない。
非コンパクト性は稠密性のために使った。
一方、局所コンパクト性は Hausdorff 性の証明で使った。
後の :ref[prop-top-lch-subspace] で示すように、Hausdorff 空間 $X$ について $\alpha X$ が Hausdorff になることと $X$ が局所コンパクトであることは同値である。
:::

:::example[一点コンパクト化の具体例]{#exm-top-one-point-examples}
1. $n\ge1$ のとき $\alpha\R^n\cong S^n$ である。北極 $N=(0,\dots,0,1)$ からの立体射影は $S^n\setminus\{N\}$ と $\R^n$ の同相写像を与える（[連続写像と同相](/math/topology/continuous-maps/) の立体射影の例を参照）。逆写像は $v\mapsto\bigl(2v/(1+\norm v^2),\ (\norm v^2-1)/(1+\norm v^2)\bigr)$ である。$S^n$ は Heine–Borel の定理によりコンパクトな距離空間であり、$S^n\setminus\{N\}$ はその部分空間で補集合は一点である。一意性により $\alpha\R^n\cong S^n$ となる。特に $\R$ の「両側の無限遠」は円周の同じ一点 $N$ に対応する。
2. 正整数全体の離散空間 $\N$ では、コンパクト集合はちょうど有限集合である（無限部分集合は一点集合による被覆が有限化できない）。従って $\infty$ の開近傍は有限個を除くすべての整数と $\infty$ からなる集合を含む。写像 $n\mapsto1/n$、$\infty\mapsto0$ は $\alpha\N$ から $\{0\}\cup\{1/n:n\ge1\}$ への全単射であり、後者はコンパクト（Heine–Borel）で $\{1/n\}$ は孤立点からなる開部分空間なので、一意性により同相である。
3. $\alpha(0,1)\cong S^1$、$\alpha[0,1)\cong[0,1]$ である。前者は $(0,1)\cong\R$ と 1 から、後者は $[0,1]$ がコンパクト Hausdorff で $[0,1)$ がその部分空間、補集合が一点 $\{1\}$ であることから、一意性による。
4. 列 $x_n\in\R$ が $\alpha\R$ で $\infty$ に収束することは、任意のコンパクト集合 $K$ について十分後の項が $K$ の外にあることであり、$K$ として $[-M,M]$ を取れば $|x_n|\to\infty$ と同値である。
:::

:::proposition[局所コンパクト Hausdorff 空間の部分空間]{#prop-top-lch-subspace}
$X$ を局所コンパクト Hausdorff 空間とする。

1. $X$ の開部分空間は局所コンパクト Hausdorff である。
2. $X$ の閉部分空間は局所コンパクト Hausdorff である。
3. 位相空間 $Z$ が局所コンパクト Hausdorff であるための必要十分条件は、$Z$ があるコンパクト Hausdorff 空間の開部分空間と同相であることである。
:::

:::proof
Hausdorff 空間の部分空間は Hausdorff なので、局所コンパクト性を示せばよい。

1. $U\subset X$ を開とし、$x\in U$ とする。
:ref[lem-top-compact-neighborhood] により $x\in V\subset\overline V\subset U$、$\overline V$ コンパクトとなる開集合 $V$ がある。
$V$ は $U$ でも開であり、$\overline V\subset U$ は $U$ の部分集合としてコンパクトなので、$\overline V$ は $x$ の $U$ におけるコンパクトな近傍である。

2. $F\subset X$ を閉とし、$x\in F$ とする。
$x$ のコンパクトな近傍 $K$ と開集合 $O$（$x\in O\subset K$）を取る。
$F\cap K$ はコンパクト空間 $K$ の閉部分集合なのでコンパクトであり、$F\cap O$ は $F$ の相対開集合で $x\in F\cap O\subset F\cap K$ である。
よって $F\cap K$ は $x$ の $F$ におけるコンパクトな近傍である。

3. コンパクト Hausdorff 空間 $Y$ は局所コンパクト Hausdorff であり、1 によりその開部分空間も局所コンパクト Hausdorff である。この性質は同相写像で保たれる。
逆に $Z$ が局所コンパクト Hausdorff なら、$Z$ がコンパクトのときは $Z$ 自身を、そうでないときは :ref[thm-top-one-point] の $\alpha Z$ を取れば、$Z$ はコンパクト Hausdorff 空間の開部分空間である。
:::

1 と 2 から、開集合と閉集合の共通部分 $U\cap F$ も局所コンパクト Hausdorff である（$U\cap F$ は局所コンパクト Hausdorff 空間 $F$ の開部分空間）。
一方、一般の部分空間では成り立たない。局所コンパクトな $\R$ の部分空間 $\Q$ は局所コンパクトでない。
また 3 の必要性から、:ref[thm-top-one-point] の注意で述べた主張が従う。
Hausdorff 空間 $X$ の $\alpha X$ が Hausdorff なら、$X$ はコンパクト Hausdorff 空間 $\alpha X$ の開部分空間なので局所コンパクトである。

### 無限遠での連続性と固有写像

:::proposition[一点コンパクト化への延長と固有写像]{#prop-top-proper-extension}
$X,Y$ をコンパクトでない局所コンパクト Hausdorff 空間とし、$f:X\to Y$ を連続写像とする。
$f$ を $\infty_X\mapsto\infty_Y$ で延長した写像 $\bar f:\alpha X\to\alpha Y$ が連続であるための必要十分条件は、$Y$ の任意のコンパクト集合 $K$ について $f^{-1}(K)$ が $X$ でコンパクトであることである。
この条件を満たす連続写像を:term[固有写像]{en="proper map" yomi="こゆうしゃぞう"}という。
:::

:::proof
$\alpha Y$ の開集合のうち 1 の形のもの、すなわち $Y$ の開集合 $V$ については、$\bar f^{-1}(V)=f^{-1}(V)$ は $X$ で開なので $\alpha X$ でも開である。
2 の形 $\{\infty_Y\}\cup(Y\setminus K)$ については

$$
\bar f^{-1}\bigl(\{\infty_Y\}\cup(Y\setminus K)\bigr)=\{\infty_X\}\cup\bigl(X\setminus f^{-1}(K)\bigr)
$$

である。
$f^{-1}(K)$ がコンパクトなら右辺は 2 の形で開である。
逆に右辺が $\alpha X$ で開なら、$\infty_X$ を含むので 2 の形 $\{\infty_X\}\cup(X\setminus C)$（$C$ はコンパクト）であり、$X$ の部分を比べて $f^{-1}(K)=C$ はコンパクトである。
従って $\bar f$ が連続であることと $f$ が固有であることは同値である。
:::

:::example[固有写像と固有でない写像]
1. $f(x)=x^2$ は $\R$ から $\R$ への固有写像である。コンパクト集合 $K$ はある $[-M,M]$（$M\ge0$）に含まれ、$f^{-1}(K)$ は閉で $f^{-1}([-M,M])=[-\sqrt M,\sqrt M]$ に含まれるので有界閉、すなわちコンパクトである。従って $f$ は $\alpha\R\cong S^1$ の自己写像に延長できる。
2. $\sin:\R\to\R$ は固有でない。コンパクト集合 $\{0\}$ の逆像 $\pi\Z$ は非有界だからである。
3. 定数でない複素係数多項式 $p(z)=a_nz^n+\dots+a_0$（$n\ge1$、$a_n\ne0$）は $\C$ から $\C$ への固有写像である。$R=\max\bigl(1,\ 2\sum_{k<n}|a_k|/|a_n|\bigr)$ と置くと、$|z|\ge R$ のとき

$$
\Bigl|\sum_{k<n}a_kz^k\Bigr|\le\Bigl(\sum_{k<n}|a_k|\Bigr)|z|^{n-1}\le\frac{|a_n|}2|z|^n
$$

なので $|p(z)|\ge|a_n||z|^n/2$ である。
コンパクト集合 $K\subset\{|w|\le M\}$ の逆像は閉であり、$|z|\le\max\bigl(R,(2M/|a_n|)^{1/n}\bigr)$ に含まれるので有界閉、すなわちコンパクトである。
従って $p$ は $\alpha\C\cong S^2$（Riemann 球面）の連続な自己写像に延長される。
:::

### 局所コンパクトな因子と商写像の積

一般に商写像と恒等写像の積は商写像とは限らないが、局所コンパクト Hausdorff 空間との積なら商写像になる。
これは商空間の積を扱うときに有用である。
局所コンパクト性を外したときの反例は :ref[exm-top-quotient-product-failure] で構成した。

:::proposition[局所コンパクトな因子との商写像の積]{#prop-top-quotient-product}
$q:X\to Y$ を商写像、$L$ を局所コンパクト Hausdorff 空間とする。
このとき $q\times\id_L:X\times L\to Y\times L$、$(x,t)\mapsto(q(x),t)$ は商写像である。
:::

:::proof
$p=q\times\id_L$ は座標がともに連続なので連続であり、$q$ が全射なので全射である。
$W\subset Y\times L$ について $G=p^{-1}(W)$ が開であると仮定し、$W$ が開であることを示す。
$(y,t)\in W$ を固定し、$q(x_0)=y$ となる $x_0\in X$ を一つ選ぶ。

**切り口が開であること。** $W_y=\{s\in L:(y,s)\in W\}$ は、連続写像 $s\mapsto(x_0,s)$ による $G$ の逆像に等しい（$(x_0,s)\in G\iff(q(x_0),s)=(y,s)\in W$）ので、$L$ で開であり $t$ を含む。
:ref[lem-top-compact-neighborhood] により、$t\in O\subset\overline O\subset W_y$、$\overline O$ コンパクトとなる開集合 $O$ がある。$K=\overline O$ と置く。

**$K$ に沿った管。** $V=\{z\in Y:\{z\}\times K\subset W\}$ と置くと、$K\subset W_y$ なので $y\in V$ である。
$x\in X$ について

$$
x\notin q^{-1}(V)\iff\exists s\in K,\ (q(x),s)\notin W\iff\exists s\in K,\ (x,s)\notin G
$$

なので、$X\setminus q^{-1}(V)$ は閉集合 $(X\times K)\setminus G$ の射影 $X\times K\to X$ による像に等しい。
$K$ はコンパクトなので、:ref[lem-top-tube] によりこの射影は閉写像であり、$q^{-1}(V)$ は $X$ で開である。
$q$ は商写像なので $V$ は $Y$ で開である。
$V\times O$ は $(y,t)$ の開近傍で、$z\in V$、$s\in O\subset K$ なら $(z,s)\in W$ なので $V\times O\subset W$ である。
$(y,t)$ は任意だったので $W$ は開である。
:::

例えば $[0,1]$ の両端を同一視して円周を作り（:ref[exm-top-circle-quotient]）、それに $[0,1]$ を掛けた円筒 $S^1\times[0,1]$ は、先に正方形 $[0,1]^2$ を作ってから左右の辺を同一視した商と同じ位相をもつ。
命題により $q\times\id_{[0,1]}$ は商写像だからである。

### コンパクト Hausdorff 空間の成分と準成分

連結成分は「内側から」連結集合の和として定義した（:ref[def-top-component]）。
「外側から」開かつ閉な集合で削っていくと別の集合が得られる。

:::definition[準成分]
位相空間 $X$ の点 $x$ について、$x$ を含む開かつ閉な集合すべての共通部分 $Q(x)$ を $x$ の:term[準成分]{en="quasicomponent" yomi="じゅんせいぶん"}という。
:::

連結成分は常に準成分に含まれる。
$x$ を含む開かつ閉な集合 $H$ について、連結集合 $C(x)$ が $H$ と $X\setminus H$ の両方と交われば、$C(x)\cap H$ と $C(x)\setminus H$ が $C(x)$ の分離になるからである。
一般には等号は成り立たないが、コンパクト Hausdorff 空間では成り立つ。

:::proposition[コンパクト Hausdorff 空間の成分と準成分]{#prop-top-quasicomponent}
コンパクト Hausdorff 空間 $K$ の各点 $x$ について $Q(x)=C(x)$ である。
:::

:::proof
$C(x)\subset Q(x)$ は上で示したので、$Q(x)$ が連結であることを示せばよい。
そうすれば $Q(x)$ は $x$ を含む連結集合であり、連結成分の最大性から $Q(x)\subset C(x)$ となる。

$Q(x)$ が連結でないと仮定し、$Q(x)=A\cup B$ を分離とする（$A,B$ は空でなく交わらない $Q(x)$ の相対開集合）。$x\in A$ としてよい。
$Q(x)$ は閉集合の共通部分なので $K$ で閉であり、$A,B$ は $Q(x)$ で相対閉なので $K$ でも閉である。
従って $A,B$ はコンパクトであり、:ref[lem-top-compact-separation] により交わらない開集合 $U\supset A$、$V\supset B$ がある。

$L=K\setminus(U\cup V)$ は閉なのでコンパクトであり、$Q(x)\subset U\cup V$ なので $L\cap Q(x)=\varnothing$ である。
各 $z\in L$ は $Q(x)$ に属さないので、$x$ を含み $z$ を含まない開かつ閉な集合 $H_z$ がある。
$\{K\setminus H_z:z\in L\}$ は $L$ の開被覆なので有限個 $z_1,\dots,z_m$ で覆え、$H=H_{z_1}\cap\dots\cap H_{z_m}$ と置く（$L=\varnothing$ なら $H=K$）。
$H$ は開かつ閉で、$x\in H$、$H\cap L=\varnothing$、すなわち $H\subset U\cup V$ である。

$H\cap U$ を考える。これは開である。
また $H\subset U\cup V$、$U\cap V=\varnothing$ より $H\setminus(H\cap U)=H\cap V$ は開なので、$K\setminus(H\cap U)=(K\setminus H)\cup(H\cap V)$ は開であり、$H\cap U$ は閉でもある。
$x\in A\subset U$ かつ $x\in H$ なので $x\in H\cap U$ であり、準成分の定義から $Q(x)\subset H\cap U$ である。
ところが $B\subset Q(x)$ は空でなく $B\subset V$ なので $B\cap U=\varnothing$ となり、矛盾する。
:::

コンパクト性がないと準成分が成分より大きくなる例は :ref[exm-top-quasicomponent-teeth] で計算した。

## 演習問題

::::exercise{level=basic}
$K=\{0\}\cup\{1/n:n\ge1\}\subset\R$ がコンパクトであることを、開被覆の定義から直接示せ。
また $\{1/n:n\ge1\}$ がコンパクトでないことを示せ。

:::solution
$\mathcal U$ を $\R$ の開集合による $K$ の被覆とする。
$0$ を含む $U_0\in\mathcal U$ を取ると、$U_0$ は開なので $(-\eps,\eps)\subset U_0$ となる $\eps>0$ がある。
$1/N<\eps$ となる $N$ を取れば、$n\ge N$ の $1/n$ はすべて $U_0$ に属する。
残りの有限個の点 $1,1/2,\dots,1/(N-1)$ について、それぞれを含む $\mathcal U$ の元を一つずつ選べば、$U_0$ と合わせて有限部分被覆になる。

一方 $A=\{1/n:n\ge1\}$ では、各 $1/n$ を含む開区間 $U_n=(1/n-\delta_n,1/n+\delta_n)$、$\delta_n=1/(2n(n+1))$ を取ると、$U_n$ は $A$ の点を $1/n$ 一つしか含まない（隣の点 $1/(n+1)$ との距離は $1/(n(n+1))=2\delta_n$ で、$1/(n-1)$ との距離はそれより大きい）。
$(U_n)$ は $A$ を覆うが、有限部分族は有限個の点しか覆わない。
:::
::::

::::exercise{level=basic}
$[0,1]$ の相対開被覆 $U=[0,2/3)$、$V=(1/3,1]$ について、「直径が $\lambda$ 未満の空でない集合は $U$ か $V$ に含まれる」を満たす $\lambda$ の最大値を求めよ。

:::solution
最大値は $1/3$ である。
$A$ が $U$ にも $V$ にも含まれなければ、$A\not\subset U$ から $a\ge2/3$ となる $a\in A$ が、$A\not\subset V$ から $b\le1/3$ となる $b\in A$ がある。
従って $\operatorname{diam}A\ge a-b\ge1/3$ である。
対偶を取れば、直径が $1/3$ 未満の空でない集合は $U$ か $V$ に含まれ、$\lambda=1/3$ は条件を満たす。
$\lambda>1/3$ では、$A=\{1/3,2/3\}$ の直径 $1/3$ は $\lambda$ 未満であるが、$2/3\notin U$、$1/3\notin V$ なので条件を満たさない。
:::
::::

::::exercise{level=basic}
次の部分集合がコンパクトかどうか判定せよ。

1. $\{(x,y)\in\R^2:xy=1,\ 1\le x\le2\}$
2. $\{(x,y)\in\R^2:xy=1\}$
3. $\{(x,y)\in\R^2:x^2+4y^2\le4\}$
4. $\Q\cap[0,1]$（$\R$ の部分集合として）

:::solution
1. 連続関数 $xy$ と $x$ で定まる閉条件の共通部分なので閉であり、$1\le x\le2$ のとき $y=1/x\in[1/2,1]$ なので有界である。Heine–Borel の定理によりコンパクトである。
2. 閉であるが、$(n,1/n)$ を含み有界でないのでコンパクトでない。
3. 連続関数 $x^2+4y^2$ による閉半直線 $(-\infty,4]$ の逆像なので閉であり、$|x|\le2$、$|y|\le1$ なので有界である。よってコンパクトである。
4. $q_n=10^{-n}\lfloor10^n\sqrt2/2\rfloor$ は $\Q\cap[0,1]$ の点列で、その極限 $\sqrt2/2$ は集合に属さない。従って閉でない（:ref[prop-top-metric-closed-sequence]）。従ってコンパクトでない。
:::
::::

::::exercise{level=basic}
$X$ をコンパクト空間、$f:X\to\R$ を連続関数とし、すべての $x$ で $f(x)>0$ とする。
ある $c>0$ があってすべての $x$ で $f(x)\ge c$ となることを示せ。
$X=(0,1]$、$f(x)=x$ ではどうなるか。

:::solution
$X=\varnothing$ なら $c=1$ でよい。
$X\ne\varnothing$ なら :ref[thm-top-compact-image] の 2 により $f$ は最小値 $c=f(x_0)$ をとり、仮定から $c>0$ である。
別証明として、開集合 $U_n=\{x:f(x)>1/n\}$ は増大列で $X$ を覆うので、有限部分被覆を取ると一つの $U_N=X$ となり、$c=1/N$ が取れる。
$X=(0,1]$、$f(x)=x$ では $f>0$ だが $\inf f=0$ であり、そのような $c$ は存在しない。$(0,1]$ がコンパクトでないからである。
:::
::::


::::exercise{level=advanced}
[部分空間・積空間・商空間](/math/topology/constructions/) の原点を二つもつ直線 $X$（:ref[exm-top-double-origin]）において、二つのコピーの閉区間 $[-1,1]$ の像 $K_a,K_b$ はそれぞれコンパクトであるが、$K_a\cap K_b$ はコンパクトでないことを示せ。
これは :ref[prop-top-compact-closed] の 4 で Hausdorff 性が必要であることを示す。

:::solution
商写像を $q$ とし、第一のコピーの包含を $i_a:\R\to X$ とする。$i_a$ は連続なので、$K_a=i_a([-1,1])$ はコンパクト集合の連続像としてコンパクトである。$K_b$ も同様である。
$K_a\cap K_b$ は二つの原点を除いた $[-1,1]\setminus\{0\}$ の像 $E$ である。
$E$ の点と $[-1,1]\setminus\{0\}$ の点は $i_a$ で一対一に対応し、$X$ の開集合 $q(\text{二つのコピーの }(c,d)\setminus\{0\})$ を考えれば、$E$ の相対位相は $[-1,1]\setminus\{0\}$ の通常の位相と一致する（$i_a$ は $\R\setminus\{0\}$ から $X\setminus\{0_a,0_b\}$ の上への同相写像である）。
$[-1,1]\setminus\{0\}$ は $\R$ で閉でないのでコンパクトでない。例えば $\{(1/n,2)\cup(-2,-1/n):n\ge2\}$ は有限化できない被覆である。
従って $K_a\cap K_b$ はコンパクトでない。
:::
::::


::::exercise{level=advanced}
コンパクト距離空間 $(X,d)$ の自己写像 $f$ が等長（$d(f(x),f(y))=d(x,y)$）なら、$f$ は全射であることを示せ。
$X=[0,\infty)$、$f(x)=x+1$ では成り立たないことも確かめよ。

:::solution
$f$ は連続なので $f(X)$ はコンパクトであり、閉である。
$f(X)\ne X$ と仮定し、$x\in X\setminus f(X)$ を取る。
$X\setminus f(X)$ は開なので $B(x,\delta)\cap f(X)=\varnothing$ となる $\delta>0$ がある。
$x_n=f^n(x)$（$f^0=\id$）と置く。$m>n\ge0$ なら、等長性を $n$ 回使って

$$
d(x_m,x_n)=d(f^n(f^{m-n}(x)),f^n(x))=d(f^{m-n}(x),x)\ge\delta
$$

である。$f^{m-n}(x)\in f(X)$ だからである。
従って $(x_n)$ のどの部分列も Cauchy 列でなく収束しないが、これは :ref[thm-top-metric-compactness] の点列コンパクト性に反する。

$X=[0,\infty)$ では $f(x)=x+1$ は等長であるが、像は $[1,\infty)$ で全射でない。$X$ はコンパクトでない。
:::
::::


::::exercise{level=advanced}
距離空間 $X$ のコンパクト集合 $A$ と閉集合 $B$ が交わらないとき、$d(A,B)=\inf\{d(a,b):a\in A,\ b\in B\}>0$ であることを示せ（$A,B\ne\varnothing$）。
両方が閉集合なだけでは成り立たない例を $\R^2$ で挙げよ。

:::solution
関数 $g(a)=d(a,B)=\inf_{b\in B}d(a,b)$ は $|g(a)-g(a')|\le d(a,a')$ を満たすので連続である。
$B$ は閉で $a\notin B$ なので $g(a)>0$ である（$g(a)=0$ なら $a\in\overline B=B$）。
$A$ はコンパクトなので $g$ は $A$ 上で最小値 $g(a_0)>0$ をとり、$d(A,B)=\min_Ag=g(a_0)>0$ である。

反例として $A=\{(x,0):x\in\R\}$、$B=\{(x,y):xy=1,\ x>0\}$ を取る。どちらも閉で交わらないが、$d((n,0),(n,1/n))=1/n\to0$ なので $d(A,B)=0$ である。
:::
::::

::::exercise{level=advanced}
コンパクト空間 $X$ 上の連続関数列 $f_n:X\to\R$ が各点で単調減少し、連続関数 $f$ に各点収束するとき、一様収束することを示せ（Dini の定理）。
$X=[0,1)$、$f_n(x)=x^n$ ではどうなるか。

:::solution
$g_n=f_n-f$ は連続で、各点で単調減少して $0$ に収束するので $g_n\ge0$ である。
$\eps>0$ を固定し、$U_n=\{x:g_n(x)<\eps\}$ と置くと、$U_n$ は開であり、単調性から $U_1\subset U_2\subset\cdots$ である。
各点 $x$ で $g_n(x)\to0$ なので、ある $n$ で $x\in U_n$ となり、$(U_n)$ は $X$ を覆う。
有限部分被覆を取ると、増大列なので最大の添字 $N$ について $U_N=X$ である。
$n\ge N$ なら $U_n\supset U_N=X$ なので、すべての $x$ で $0\le g_n(x)<\eps$ である。
よって $\sup_x|f_n(x)-f(x)|\le\eps$（$n\ge N$）であり、一様収束する。

$X=[0,1)$ では $x^n$ は各点で単調減少して $0$ に収束するが、$\sup_{[0,1)}x^n=1$ なので一様収束しない。$U_n=\{x<\eps^{1/n}\}$ は有限個で $[0,1)$ を覆えない。
:::
::::

::::exercise{level=advanced}
$X$ を Hausdorff 空間、$K_1\supset K_2\supset\cdots$ をコンパクト集合の減少列とし、開集合 $U$ が $\bigcap_nK_n\subset U$ を満たすとする。
ある $N$ で $K_N\subset U$ となることを示せ。
これを用いて、$\bigcap_nK_n$ が一点 $\{x\}$ なら、$x$ の任意の近傍が十分大きい $n$ の $K_n$ を含むことを導け。

:::solution
$F_n=K_n\setminus U=K_n\cap(X\setminus U)$ と置く。
$K_n$ は Hausdorff 空間のコンパクト集合なので閉であり、$F_n$ はコンパクト空間 $K_1$ の閉集合の減少列である。
$\bigcap_nF_n=\bigl(\bigcap_nK_n\bigr)\setminus U=\varnothing$ である。
すべての $F_n$ が空でなければ、:ref[cor-top-cantor-compact] の 1（$K_1$ に適用）により共通部分は空でないので矛盾する。
従ってある $N$ で $F_N=\varnothing$、すなわち $K_N\subset U$ である。
後半は、近傍 $W$ に含まれる開近傍 $U\ni x$ に前半を適用すればよい。減少列なので $n\ge N$ で $K_n\subset K_N\subset U\subset W$ である。
:::
::::

::::exercise{level=advanced}
任意の $a\in\R$ について $\{x:f(x)<a\}$ が開であるような関数 $f:X\to\R$ を上半連続という。
空でないコンパクト空間上の上半連続関数は最大値をとることを示せ。
最小値は必ずしもとらないことを例で示せ。

:::solution
$U_n=\{x:f(x)<n\}$（$n\ge1$）は開で $X$ を覆うので、有限個で覆え、増大列なので一つの $U_N=X$ となる。従って $f<N$ で、$s=\sup f$ は有限である。
$s$ が値として取られないと仮定すると、$V_m=\{x:f(x)<s-1/m\}$（$m\ge1$）は開で、各 $x$ で $f(x)<s$ なので $X$ を覆う。
有限部分被覆を取ると一つの $V_M=X$ となり、$f<s-1/M$ がすべての点で成り立って $s$ が上限であることに反する。
従って $f$ は最大値 $s$ をとる。

最小値については、$X=[0,1]$ 上で $f(x)=x$（$x>0$）、$f(0)=1$ と置くと、$\{f<a\}$ は $a\le0$ で空、$0<a\le1$ で $(0,a)$、$a>1$ で $[0,1]$ となりすべて開なので上半連続である。$\inf f=0$ であるが値 $0$ はとらない。
:::
::::

## 確認テスト

::::quiz
可算コンパクトであるがコンパクトでない空間の例として、正しい根拠をもつものを選べ。

- [x] $[0,\omega_1)$：各点列は一つのコンパクトな初切片 $[0,\alpha]$ に入るので集積点をもつが、初切片 $[0,\beta)$ による開被覆は有限化できない。
- [ ] $\R$：任意の点列が有界だから集積点をもつ。
- [ ] 有限集合に離散位相を入れた空間：一点集合による被覆が有限化できない。

:::explanation
正解は最初の選択肢である。
$[0,\omega_1)$ の点列 $(\beta_n)$ の上限 $\alpha$ は可算順序数であり、全項はコンパクトな $[0,\alpha]$ に入るので、:ref[prop-top-countably-compact] により集積点をもち、空間は可算コンパクトである。
一方 $\{[0,\beta):\beta<\omega_1\}$ は開被覆で、有限個の和は $[0,\max\beta_j)$ となり点 $\max\beta_j$ を覆わない。
$\R$ の点列 $x_n=n$ は集積点をもたないので $\R$ は可算コンパクトでない。
有限空間は常にコンパクトであり、一点集合による被覆は最初から有限である。
:::
::::

::::quiz{answer="3" tol="1e-9"}
直交群 $O(3)$ を $\R^9$ の部分集合とみたとき、その元 $A=(a_{ij})$ について $\sum_{i,j}a_{ij}^2$ の値はいくらか（この値が一定であることが $O(3)$ の有界性を与える）。

:::explanation
$\transp AA=I$ の対角成分を見ると、各列 $a_{\cdot j}$ は $\sum_ia_{ij}^2=1$ を満たす単位ベクトルである。
列は三本なので $\sum_{i,j}a_{ij}^2=3$ である。
従って $O(3)$ は半径 $\sqrt3$ の球面に含まれて有界であり、閉でもあるので Heine–Borel の定理によりコンパクトである。
:::
::::

::::quiz
距離空間 $X$ について、コンパクト性と同値な条件をすべて選べ。

- [x] 点列コンパクトである。
- [x] 完備かつ全有界である。
- [ ] 完備かつ有界である。
- [x] 任意の無限部分集合が集積点をもつ。
- [ ] 全有界である。

:::explanation
:ref[thm-top-metric-compactness] により、コンパクト・極限点コンパクト（任意の無限部分集合が集積点をもつ）・点列コンパクト・完備かつ全有界は同値である。
完備かつ有界だけでは不十分で、無限離散距離空間が反例である（完備、直径 $1$ だが、半径 $1/2$ の球は一点なので全有界でない）。
全有界だけでも不十分で、$(0,1)$ が反例である（全有界だが完備でない）。
:::
::::

::::quiz{answer="1" tol="1e-9"}
離散空間 $\N$ の一点コンパクト化 $\alpha\N$ において、孤立点でない点はいくつあるか。

:::explanation
各自然数 $n$ について $\{n\}$ は $\N$ の開集合なので $\alpha\N$ でも開であり、$n$ は孤立点である。
$\infty$ の開近傍は $\{\infty\}\cup(\N\setminus K)$（$K$ は有限集合）の形なので、無限個の自然数を含み、$\{\infty\}$ は開でない。
従って孤立点でない点は $\infty$ の一つだけである。
実際 $\alpha\N$ は $\{0\}\cup\{1/n:n\ge1\}$ と同相であり、$0$ だけが孤立点でない。
:::
::::

::::quiz
Tychonoff の定理（任意個のコンパクト空間の積はコンパクト）について正しい記述を選べ。

- [ ] 添字集合が可算であることが必要である。
- [ ] 各因子が Hausdorff であることが必要である。
- [x] 積には積位相を入れる必要があり、箱位相では一般に成り立たない。

:::explanation
正解は三つ目である。
超フィルターによる証明では、各座標で射影した超フィルターの極限を一つずつ選び、積位相の基本近傍が有限個の座標しか制限しないことと、フィルターが有限個の共通部分で閉じていることを組み合わせた。
添字の可算性も極限の一意性（Hausdorff 性）も使っていない。
箱位相では $\{0,1\}^{\N}$ が無限離散空間になり、コンパクトでない。
:::
::::

::::quiz
局所コンパクトな空間をすべて選べ。

- [x] $\R^n$
- [ ] $\Q$（通常の位相）
- [ ] $\R^{\N}$（積位相）
- [x] 無限集合に離散位相を入れた空間
- [x] $\R^2$ の開円板

:::explanation
$\R^n$ とその開集合（開円板など）は、各点の周りに閉球というコンパクトな近傍をもつ。
離散空間では一点集合がコンパクトな開近傍である。
$\Q$ のコンパクトな近傍は $\R$ で閉なので有理数の開区間の閉包、すなわち無理数を含む閉区間を含むことになり矛盾する。
$\R^{\N}$ では、基本開集合が制限しない座標への射影の像が $\R$ 全体になり、コンパクトな近傍の連続像がコンパクトでない $\R$ を含むことになって矛盾する。
:::
::::

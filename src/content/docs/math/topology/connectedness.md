---
title: 連結性
description: 開集合による分離と道による接続を比較し、区間の連結性、成分、積、Cantor 集合、位相幾何学者の正弦曲線、局所連結性を調べる。
notation:
  - symbol: 'C(x)'
    meaning: 点 x の連結成分
  - symbol: 'P(x)'
    meaning: 点 x の弧状連結成分
  - symbol: '\gamma\ast\eta'
    meaning: 終点と始点が一致する二つの道の連結
  - symbol: '\gamma^{-}'
    meaning: 道 γ を逆向きにたどる道
  - symbol: 'K'
    meaning: Cantor 集合
  - symbol: '\GL_n^{+}(\R)'
    meaning: 行列式が正の n 次正則行列全体
---

## この章の目標

- 分離、開かつ閉な集合、二点離散空間への連続写像の三つの方法で連結性を判定し、実数の連結部分集合がちょうど区間であることを証明する。
- 中間値の定理を連結性の帰結として捉え、連結性・切断点・成分の個数を同相の判別に使う。
- 連結成分を構成し、連続写像・閉包・任意積に対する振る舞いを証明し、Cantor 集合が $\{0,1\}^{\N}$ と同相な全不連結空間であることを示す。
- 道による接続（弧状連結性）を導入し、位相幾何学者の正弦曲線と辞書式順序の正方形で連結性との違いを確かめる。
- 局所連結性・局所弧状連結性を開集合の成分で特徴付け、商空間と積空間での振る舞いを決定する。

前章までに、位相空間・連続写像・部分空間・積・商を準備した（[部分空間・積空間・商空間](/math/topology/constructions/)）。
本章で扱う連結性は「空間が一続きであるか」を開集合だけで表す性質であり、
解析学の中間値の定理の背後にある位相的な事実である。
連続写像で保たれるので、同相でない空間を見分ける最初の道具にもなる。

## 連結空間

### 分離と連結性の定義

直線から原点を除くと、残りは負の側 $(-\infty,0)$ と正の側 $(0,\infty)$ に分かれる。
どちらの側も開集合であり、一方の点から少し動いただけでは他方に移れない。
逆に区間 $[0,1]$ をこのように二つの開集合に分けることはできない、というのが直観である。
この「全体が交わらない二つの空でない開集合の和に分かれる」という現象を、
距離や線分を使わずに開集合だけで表したものが分離である。

:::definition[分離と連結性]{#def-top-connected}
位相空間 $X$ の部分集合 $U,V$ が

$$
X=U\cup V,\qquad U\cap V=\varnothing,\qquad U\ne\varnothing,\qquad V\ne\varnothing
$$

を満たし、$U,V$ がともに開集合であるとき、組 $(U,V)$ を $X$ の分離という。
分離をもたない空間を:term[連結空間]{en="connected space" yomi="れんけつくうかん"}という。
部分集合 $A\subset X$ が連結であるとは、$A$ に相対位相を入れた部分空間が連結であることをいう。
:::

定義から、空集合と一点空間は連結である（空でない二つの集合に分けることができない）。
文献によっては空集合を連結と呼ばない流儀もあるが、本章では上の定義に従う。
部分集合の連結性は相対位相で判定することに注意する。
例えば $A=(-1,0)\cup(0,1)\subset\R$ では、$(-1,0)$ と $(0,1)$ は $A$ の相対開集合であり、
$A$ の分離を与えるので $A$ は連結でない。

:::definition[開かつ閉な集合]
$X$ で開集合でもあり閉集合でもある部分集合を:term[開かつ閉な集合]{en="clopen set" yomi="かいかつへいなしゅうごう"}という。
$\varnothing$ と $X$ は常に開かつ閉であり、これらを自明な開かつ閉集合という。
:::

:::proposition[連結性の言い換えと連続像]{#prop-top-connected-image}
位相空間 $X$ について、次の四条件は同値である。

1. $X$ は連結である。
2. $X$ の開かつ閉な集合は $\varnothing$ と $X$ だけである。
3. 任意の離散空間 $D$ への任意の連続写像 $f:X\to D$ の像は高々一点である。
4. 二点離散空間 $\{0,1\}$ への任意の連続写像の像は高々一点である。

さらに、連結空間 $X$ から位相空間 $Y$ への連続写像 $f$ について、像 $f(X)$ は連結である。
:::

:::proof
**(1)⇒(2)。** $A$ が開かつ閉で $A\ne\varnothing,\ A\ne X$ とする。
補集合 $X\setminus A$ も開かつ閉で空でない。
すると $(A,X\setminus A)$ は分離であり、$X$ が連結であることに反する。

**(2)⇒(1)。** 分離 $(U,V)$ があったとする。
$U=X\setminus V$ であり $V$ は開だから $U$ は閉である。
$U$ はもともと開なので開かつ閉であり、$U\ne\varnothing$、$V\ne\varnothing$ より $U\ne X$ である。
これは (2) に反する。

**(2)⇒(3)。** $f:X\to D$ が連続で、相異なる二値 $y,z$ を取るとする。
離散空間では $\{y\}$ は開かつ閉であるから、$A=f^{-1}(\{y\})$ は開であり、
$X\setminus A=f^{-1}(D\setminus\{y\})$ も開なので $A$ は閉でもある。
$A$ は $f$ が $y$ を取る点を含むので空でなく、$f$ が $z$ を取る点を含まないので $X$ でもない。
これは (2) に反する。

**(3)⇒(4)。** $\{0,1\}$ は離散空間なので特別な場合である。

**(4)⇒(2)。** $A$ を開かつ閉な集合とし、特性関数 $\chi_A:X\to\{0,1\}$（$A$ 上で $1$、外で $0$）を考える。
$\{0,1\}$ の部分集合 $\varnothing,\{1\},\{0\},\{0,1\}$ の逆像はそれぞれ $\varnothing,A,X\setminus A,X$ で、いずれも開である。
従って $\chi_A$ は連続であり、(4) より値は高々一種類である。
値が $1$ だけなら $A=X$、$0$ だけなら $A=\varnothing$ である（$X=\varnothing$ なら両方）。

**連続像。** $f:X\to Y$ を連続とし、$f$ を値域を $f(X)$ に制限した写像 $f_1:X\to f(X)$ とみなす。
:ref[prop-top-subspace-laws] により $f_1$ は連続である。
連続写像 $g:f(X)\to\{0,1\}$ を任意に取ると、$g\circ f_1:X\to\{0,1\}$ は連続だから、
$X$ の連結性と (4) より $g\circ f_1$ の値は高々一種類である。
$f_1$ は全射なので $g$ の値も高々一種類であり、(4) から $f(X)$ は連結である。
:::

特に連結性は同相写像で保たれる位相的性質である。
条件 (4) は以後の証明で最もよく使う形であり、
「連結空間上の連続な $\{0,1\}$ 値関数は定数である」と覚えておくとよい。

:::lemma[連結集合は分離の片側に入る]{#lem-top-connected-in-separation}
$Y\subset X$ が連結で、$U,V$ が $X$ の交わらない開集合で $Y\subset U\cup V$ を満たすならば、
$Y\subset U$ または $Y\subset V$ である。
:::

:::proof
$Y\cap U$ と $Y\cap V$ は $Y$ の相対開集合で、交わらず、和は $Y$ である。
両方が空でなければ $Y$ の分離となり矛盾するので、一方は空である。
$Y\cap V=\varnothing$ なら $Y\subset U$、$Y\cap U=\varnothing$ なら $Y\subset V$ である。
:::

### 周囲の空間で分離を表す

部分集合の分離は相対開集合で定義したが、周囲の空間 $X$ の閉包を使っても表せる。

:::lemma[分離された集合]{#lem-top-separated-sets}
$Y\subset X$ とする。$Y$ が連結でないための必要十分条件は、
空でない $A,B\subset X$ で

$$
Y=A\cup B,\qquad \overline A\cap B=\varnothing,\qquad A\cap\overline B=\varnothing
$$

を満たすものが存在することである（閉包は $X$ で取る）。
この条件を満たす $A,B$ を互いに分離された集合という。
:::

:::proof
:ref[prop-top-subspace-laws] より、$C\subset Y$ の $Y$ での閉包は $Y\cap\overline C$ である。

$(A,B)$ が $Y$ の分離であるとする。
$A=Y\setminus B$ は $Y$ で閉だから $A=Y\cap\overline A$ である。
従って $\overline A\cap B=\overline A\cap Y\cap B=A\cap B=\varnothing$ であり、
同様に $A\cap\overline B=\varnothing$ である。

逆に条件を満たす $A,B$ があるとする。
$A\cap B\subset\overline A\cap B=\varnothing$ なので $A,B$ は交わらない。
$Y\cap\overline A=(A\cup B)\cap\overline A=A\cup(B\cap\overline A)=A$ だから $A$ は $Y$ で閉であり、
同様に $B$ も $Y$ で閉である。
互いに $Y$ での補集合なので両方とも $Y$ で開でもあり、$(A,B)$ は $Y$ の分離である。
:::

:::example[閉包どうしが接しても分離できる]
$A=(-1,0)$、$B=(0,1)$ とすると、閉包 $[-1,0]$ と $[0,1]$ は点 $0$ を共有する。
それでも $\overline A\cap B=[-1,0]\cap(0,1)=\varnothing$、$A\cap\overline B=\varnothing$ であるから、
$A\cup B$ は連結でない。
一方 $B=[0,1)$ とすると $\overline A\cap B=\{0\}\ne\varnothing$ であり、
実際 $A\cup B=(-1,1)$ は次の定理により連結である。
どちらの集合に接点そのものが属するかで結論が変わる点に注意する。
「分離された」は周囲の空間で両方が開であることを要求しない。
:::

### 実数の連結部分集合

連結性の最も基本的な例は実数の区間である。
まず「区間」を正確に定義し、それが二点の間を埋める集合と一致することを確かめる。

:::definition[区間と凸集合]
$I\subset\R$ が凸であるとは、$a,b\in I$、$a<c<b$ ならば $c\in I$ となることをいう。
$a\le b$ を実数として、

$$
\varnothing,\ [a,b],\ [a,b),\ (a,b],\ (a,b),\ [a,\infty),\ (a,\infty),\ (-\infty,b],\ (-\infty,b),\ \R
$$

の形の集合を :term[区間]{en="interval" yomi="くかん"} という（$[a,a]=\{a\}$ も区間である）。
:::

:::proposition[区間の特徴付け]{#prop-top-real-convex}
$\R$ の部分集合が凸であるための必要十分条件は、それが区間であることである。
:::

:::proof
各区間が凸であることは不等式から直ちに従う。
例えば $I=[a,b)$ で $x,y\in I$、$x<c<y$ なら $a\le x<c<y<b$ である。

逆に $I$ を空でない凸集合とする。
$I$ が下に有界なら $\alpha=\inf I\in\R$、そうでなければ $\alpha=-\infty$ とし、
上に有界なら $\beta=\sup I\in\R$、そうでなければ $\beta=+\infty$ とする。
上限・下限の存在は [実数の連続性](/math/analysis/real-numbers/) による。
定義から $I$ の点 $x$ は $\alpha\le x\le\beta$ を満たす。

$\alpha<c<\beta$ を満たす実数 $c$ を取る。
$c$ は $I$ の下界でないので $a\in I$ で $a<c$ となるものがあり、
$c$ は $I$ の上界でないので $b\in I$ で $c<b$ となるものがある。
凸性から $c\in I$ である。
よって

$$
(\alpha,\beta)\subset I\subset[\alpha,\beta]\cap\R
$$

である。
$I$ と $(\alpha,\beta)$ の差は、端点 $\alpha,\beta$ のうち実数で $I$ に属するものだけである。
$\alpha=\beta$ なら $I=\{\alpha\}$ であり、$\alpha<\beta$ なら、
各端点が有限か、有限なら $I$ に属するかで場合を分けると、$I$ は上に列挙した区間のどれかに一致する。
:::

:::theorem[実数の連結部分集合]{#thm-top-interval-connected}
$\R$ の部分集合が連結であるための必要十分条件は、それが凸であることである。
従って $\R$ の連結部分集合はちょうど区間（空集合と一点を含む）である。
:::

:::proof
**連結なら凸。** $A$ が凸でないとし、$a,b\in A$、$a<c<b$、$c\notin A$ とする。

$$
U=A\cap(-\infty,c),\qquad V=A\cap(c,\infty)
$$

は $A$ の相対開集合であり、$a\in U$、$b\in V$ より空でない。
$c\notin A$ だから $A=U\cup V$ であり、$U\cap V=\varnothing$ である。
従って $(U,V)$ は $A$ の分離で、$A$ は連結でない。

**凸なら連結。** $I$ を凸集合とし、分離 $(U,V)$ があると仮定する。
$a\in U$、$b\in V$ を取る。$a\ne b$ であり、$U$ と $V$ の役割は対称なので、必要なら名前を入れ替えて $a<b$ としてよい。
凸性から $[a,b]\subset I$ である。

$$
U'=U\cap[a,b],\qquad V'=V\cap[a,b]
$$

と置くと、$U',V'$ は $[a,b]$ を交わらない二つの集合に分け、$a\in U'$、$b\in V'$ である。
また $U=I\cap O$（$O$ は $\R$ の開集合）と書けるので $U'=[a,b]\cap O$ は $[a,b]$ の相対開集合であり、$V'$ も同様である。
$U'$ は空でなく $b$ で上から抑えられるので、上限 $s=\sup U'$ が存在し、$a\le s\le b$ である。
従って $s\in U'$ または $s\in V'$ である。

$s\in U'$ の場合、$b\in V'$ だから $s\ne b$、すなわち $s<b$ である。
$U'$ は $[a,b]$ で相対開なので、ある $\delta>0$ で $(s-\delta,s+\delta)\cap[a,b]\subset U'$ となる。
$t=\min(s+\delta/2,\,b)$ と置くと $s<t\le b$、$t<s+\delta$ だから $t\in U'$ である。
これは $t>s=\sup U'$ に反する。

$s\in V'$ の場合、$a\in U'$ だから $s\ne a$、すなわち $s>a$ である。
$V'$ は相対開なので、ある $\delta>0$ で $(s-\delta,s+\delta)\cap[a,b]\subset V'$ となる。
$s$ は $U'$ の上限なので、$s-\delta$ は上界でなく、$u\in U'$ で $s-\delta<u\le s$ となるものがある。
$u\in[a,b]$ だから $u\in V'$ となり、$U'\cap V'=\varnothing$ に反する。

いずれの場合も矛盾するので $I$ は連結である。
後半は :ref[prop-top-real-convex] による。
:::

この証明で本質的に使ったのは実数の上限の存在（連続性）である。
同じ議論を有理数で行うと $\sup$ が存在しないことがあり、実際 $\Q$ は連結でない。
$\sqrt2$ が無理数なので $\Q\cap(-\infty,\sqrt2)$ と $\Q\cap(\sqrt2,\infty)$ が $\Q$ の分離になるからである。

:::corollary[中間値の定理]{#cor-top-intermediate-value}
$X$ を連結空間、$f:X\to\R$ を連続写像とする。
$a,b\in X$ と実数 $c$ が $f(a)<c<f(b)$ を満たすならば、$f(x)=c$ となる $x\in X$ が存在する。
特に区間 $I$ 上の連続関数 $f:I\to\R$ の像 $f(I)$ は区間である。
:::

:::proof
:ref[prop-top-connected-image] より $f(X)$ は $\R$ の連結部分集合であり、
:ref[thm-top-interval-connected] より凸である。
$f(a),f(b)\in f(X)$ と $f(a)<c<f(b)$ から $c\in f(X)$ である。
区間は連結なので後半も同じ理由による。
:::

:::example[方程式の解の存在と二分法]
$f(t)=t^3-t-1$ は $\R$ 上連続で、$f(1)=-1<0<5=f(2)$ である。
$[1,2]$ は連結だから、:ref[cor-top-intermediate-value] により $[1,2]$ に零点がある。
中点で符号を調べると $f(3/2)=7/8>0$、$f(5/4)=-19/64<0$、$f(11/8)=115/512>0$ なので、
零点は $[5/4,11/8]$ にある（実際の零点は $1.3247\ldots$ である）。
このように区間を半分にしながら符号の変わる側を残す二分法は、
上の証明の「上限 $s$」を数値的に近似する手続きである。
:::

:::example[中間値の定理が使えない場合]
定義域が連結でないと結論は成り立たない。
$X=[0,1]\cup[2,3]$ 上で $f(x)=x$ は連続で、値 $1$ と $2$ を取るが $3/2$ を取らない。
また値域の側を $\Q$ に変えると、$g:\Q\to\Q$、$g(q)=q^2-2$ は連続で $g(0)<0<g(2)$ だが、
$g(q)=0$ となる有理数は存在しない。
これは $\Q$ が連結でないことの言い換えである。
:::

### 連結集合の和と閉包

:::proposition[和と閉包]{#prop-top-connected-union}
$X$ を位相空間とする。

1. 連結部分集合の族 $(A_\lambda)_{\lambda\in\Lambda}$ がある共通点 $p$ をもつ（$p\in\bigcap_\lambda A_\lambda$）なら、$\bigcup_\lambda A_\lambda$ は連結である。
2. 連結部分集合の列 $A_1,A_2,\ldots$ が各 $n$ で $A_n\cap A_{n+1}\ne\varnothing$ を満たすなら、$\bigcup_n A_n$ は連結である。
3. $C$ が連結で $C\subset D\subset\overline C$ なら $D$ も連結である。特に連結集合の閉包は連結である。
:::

:::proof
いずれも :ref[prop-top-connected-image] の条件 (4) で示す。

(1) $Y=\bigcup_\lambda A_\lambda$ とし、$g:Y\to\{0,1\}$ を連続とする。
各 $A_\lambda$ への制限は連続であり、$A_\lambda$ は連結なので $g|_{A_\lambda}$ は定数である。
その値は $p\in A_\lambda$ での値 $g(p)$ に等しい。
各点 $y\in Y$ はある $A_\lambda$ に属するので $g(y)=g(p)$ であり、$g$ は定数である。

(2) $B_n=A_1\cup\cdots\cup A_n$ が連結であることを $n$ についての帰納法で示す。
$B_1=A_1$ は連結である。
$B_n$ が連結なら、$A_n\cap A_{n+1}\subset B_n\cap A_{n+1}$ は空でないので、
(1) を二つの集合 $B_n,A_{n+1}$ に使えば $B_{n+1}$ は連結である。
すべての $B_n$ は $A_1$ の一点（$A_1\cap A_2\ne\varnothing$ なので $A_1$ は空でない）を共有するから、
(1) により $\bigcup_nB_n=\bigcup_nA_n$ は連結である。

(3) $C=\varnothing$ なら $D\subset\overline\varnothing=\varnothing$ なので明らかである。
$C\ne\varnothing$ とし、$g:D\to\{0,1\}$ を連続とする。
$g|_C$ は連続で $C$ は連結なので、ある値 $c$ を取る定数である。
ある $d\in D$ で $g(d)\ne c$ だったとする。
$W=g^{-1}(\{g(d)\})$ は $D$ の相対開集合で $d$ を含むので、$X$ の開集合 $O$ を用いて $W=D\cap O$ と書ける。
$d\in O$ かつ $d\in\overline C$ だから $O\cap C\ne\varnothing$ である。
$O\cap C\subset D\cap O=W$ の点で $g$ の値は $g(d)\ne c$ となり、$g|_C\equiv c$ に反する。
:::

例えば二つの開区間 $(-1,0)$、$(0,1)$ は個々に連結だが、共通点がなく和は連結でない。
これに点 $0$ を加えた $(-1,0)\cup\{0\}$ は $(-1,0)$ とその閉包の間にあるので (3) により連結であり、
これと $(0,1)\cup\{0\}$ は $0$ を共有するので和 $(-1,1)$ は (1) により連結である。

:::proposition[連結集合は境界を横切る]{#prop-top-connected-boundary}
$C\subset X$ が連結で、$A\subset X$ について $C\cap A\ne\varnothing$ かつ $C\setminus A\ne\varnothing$ ならば、
$C\cap\partial A\ne\varnothing$ である。
:::

:::proof
$X$ は $A^\circ$、$\partial A=\overline A\setminus A^\circ$、$X\setminus\overline A$ の三つに交わらずに分かれる。
$C\cap\partial A=\varnothing$ と仮定すると $C\subset A^\circ\cup(X\setminus\overline A)$ である。
右辺の二つは交わらない開集合なので、:ref[lem-top-connected-in-separation] により $C$ は一方に含まれる。
$C\cap A$ の点は $\partial A$ にないので $A^\circ$ に属し（$A\subset\overline A$ だから $X\setminus\overline A$ には属さない）、
従って $C\cap A^\circ\ne\varnothing$ である。
$C\setminus A$ の点は $A^\circ\subset A$ に属さず、$\partial A$ にもないので $X\setminus\overline A$ に属する。
従って $C$ は両方と交わり、矛盾する。
:::

:::example[円板の内から外への曲線]
$\R^2$ で $A=\{x:\norm x<1\}$ とすると $\partial A=\{x:\norm x=1\}$ である。
実際、$\norm x<1$ の点は $A$ の内点、$\norm x>1$ の点は閉集合 $\{\norm x\le1\}$ の外にあるので $\overline A$ に属さない。
$\norm x=1$ の点 $x$ は、$(1-1/n)x\in A$ と $(1+1/n)x\notin A$ がともに $x$ に収束するので境界点である。
従って原点と点 $(2,0)$ を含む連結集合は必ず単位円周と交わる。
連続関数 $x\mapsto\norm x$ と :ref[cor-top-intermediate-value] を使っても同じ結論を得る。
:::

### 順序位相と線形連続体

:ref[thm-top-interval-connected] の証明は、実数の順序の二つの性質、
すなわち「二点の間に点があること」と「上に有界な集合が上限をもつこと」だけを使った。
これを抽象化すると、実数以外の連結な順序空間を作れる。
全順序集合の順序位相は [位相空間](/math/topology/topological-spaces/) で定義した（開半直線を準開基とする位相）。

:::definition[線形連続体]
全順序集合 $L$ が二点以上をもち、次の二条件を満たすとき :term[線形連続体]{en="linear continuum" yomi="せんけいれんぞくたい"} という。

1. $x<y$ ならば $x<z<y$ となる $z\in L$ が存在する（稠密性）。
2. 空でなく上に有界な部分集合は上限をもつ（上限性）。

$C\subset L$ が順序凸であるとは、$a,b\in C$、$a<c<b$ ならば $c\in C$ となることをいう。
:::

順序位相の端点と片側近傍には :ref[lem-top-order-neighborhood] を使う。

:::proposition[順序連続体の連結性]{#prop-top-linear-continuum}
線形連続体 $L$ に順序位相を入れると、$L$ の順序凸部分集合はすべて連結である。
特に $L$ 自身と、$L$ の開区間・閉区間・半直線は連結である。
:::

:::proof
$C$ を順序凸集合とし、相対位相での分離 $(U,V)$ があると仮定する。
$a\in U$、$b\in V$ を取り、必要なら $U,V$ を入れ替えて $a<b$ とする。
$[a,b]=\{x:a\le x\le b\}$ は凸性により $C$ に含まれる。
$U'=U\cap[a,b]$ は空でなく $b$ で上から抑えられるので、上限性から $s=\sup U'$ が存在し、
$a\le s\le b$ だから $s\in C$ である。

$s\in U$ の場合、$b\in V$ より $s<b$ であり、$s$ は $L$ の最大元でない。
$U=C\cap O$（$O$ は $L$ の開集合）とすると $s\in O$ であり、
:ref[lem-top-order-neighborhood] からある $u>s$ で $[s,u)\subset O$ となる。
$u'=\min(u,b)$ とすると $s<u'$ であり、稠密性から $s<z<u'$ となる $z$ がある。
$z\in[s,u)\subset O$、$z\in[a,b]\subset C$ だから $z\in U'$ であり、$z>s$ は上限性に反する。

$s\in V$ の場合、$a\in U$ より $a<s$ であり、$s$ は最小元でない。
$V=C\cap O'$ とすると、ある $t<s$ で $(t,s]\subset O'$ となる。
$t$ は $U'$ の上界でないので $x\in U'$ で $x>t$ となるものがあり、$x\le s$ である。
すると $x\in(t,s]\cap C\subset V$ となり、$U\cap V=\varnothing$ に反する。

どちらの場合も矛盾するので $C$ は連結である。
:::

$\R$ は線形連続体なので、この命題は :ref[thm-top-interval-connected] の一般化である。
二つの条件はどちらも外せない。
$\Z$ の順序位相は離散位相であり、隣り合う点の対 $\{0,1\}$ は凸だが連結でない（稠密性の欠如）。
$\Q$ は稠密だが $\{q\in\Q:q>0,\ q^2<2\}$ が上限をもたず、実際 $\Q$ は連結でない（上限性の欠如）。

:::example[長い半直線]{#exm-top-long-ray}
$\omega_1$ を最小の非可算順序数とし、$[0,\omega_1)$ を可算順序数全体とする。
順序数の性質（整列性、上に有界な順序数の集合が上限をもつこと、可算個の可算順序数の上限が可算であること）は
[整列集合と順序数](/math/sets-logic/ordinals/) による。
$L=[0,\omega_1)\times[0,1)$ に辞書式順序

$$
(\alpha,s)<(\beta,t)\iff \alpha<\beta\ \text{または}\ (\alpha=\beta\ \text{かつ}\ s<t)
$$

と順序位相を入れたものを長い半直線という。
これは $[0,1)$ を非可算個（$\omega_1$ 個）縦につないだ空間である。

**稠密性。** $(\alpha,s)<(\beta,t)$ とする。
$\alpha=\beta$ なら $(\alpha,(s+t)/2)$ が間にある。
$\alpha<\beta$ なら $s<(s+1)/2<1$ なので $(\alpha,(s+1)/2)$ が間にある。

**上限性。** $A\subset L$ を空でなく、$(\gamma,u)$ で上から抑えられているとする。
$A$ の第一座標全体の上限を $\beta$ とすると $\beta\le\gamma<\omega_1$ である。
(i) $\beta$ が $A$ の第一座標として現れないとき、$(\beta,0)$ が上限である。
実際、$A$ の元の第一座標は $\beta$ 未満なので $(\beta,0)$ は上界であり、
$(\delta,r)<(\beta,0)$ なら $\delta<\beta$ なので、第一座標が $\delta$ より大きい $A$ の元があり、それは $(\delta,r)$ より大きい。
(ii) $\beta$ が現れるとき、$t^*=\sup\{t:(\beta,t)\in A\}\in(0,1]$ または $0$ と置く。
$t^*<1$ なら $(\beta,t^*)$ が上限である。
$t^*=1$ なら（各 $t<1$ なので上限は達成されず）$(\beta+1,0)$ が上限である。
いずれも、より小さい元 $(\delta,r)$ を取ると、$\delta<\beta$ なら第一座標 $\beta$ の元が、
$\delta=\beta$ なら $r<t^*$ により $t>r$ となる元 $(\beta,t)\in A$ がそれを超えることから確かめられる。

従って $L$ は線形連続体であり、:ref[prop-top-linear-continuum] により連結である。

**可分でないこと。** 可算集合 $D\subset L$ の第一座標は可算個の可算順序数なので、その上限 $\gamma$ も可算順序数である。
空でない開集合 $\{x:x>(\gamma+1,0)\}$（$(\gamma+1,1/2)$ を含む）は $D$ と交わらない。
従って $L$ は可分でなく、:ref[prop-top-countable-basis] により第二可算でもない。
一方、各点の基本近傍は順序凸である。各基本近傍も順序稠密であり、その中に上界をもつ空でない部分集合の上限は順序凸性からその近傍に属する。従って同じ線形連続体の判定で基本近傍は連結である。
:::

### 細かい鎖による必要条件

距離空間では、連結性から「任意の精度の小刻みな歩みで二点を結べる」ことが従う。

:::definition[ε 鎖]
距離空間 $X$ の点列 $a=x_0,x_1,\ldots,x_m=b$ が、すべての $j$ で $d(x_{j-1},x_j)<\eps$ を満たすとき、
これを $a$ から $b$ への $\eps$ 鎖という。
途中の点はすべて $X$ に属するものとする。
:::

:::proposition[連結距離空間の鎖]{#prop-top-connected-chain}
連結な距離空間では、任意の $\eps>0$ と任意の二点 $a,b$ に対し、$a$ から $b$ への有限な $\eps$ 鎖が存在する。
:::

:::proof
$a$ を固定し、$a$ から $\eps$ 鎖で到達できる点全体を $C$ とする。
$m=0$ の鎖を許すと $a\in C$ である。

$C$ は開である。
$x\in C$、$y\in B(x,\eps)$ なら、$x$ に至る鎖の末尾に $y$ を付け加えると $y$ への $\eps$ 鎖になるので $B(x,\eps)\subset C$ である。

$X\setminus C$ も開である。
$x\notin C$ とし、$B(x,\eps)$ が $C$ の点 $y$ を含むとすると、$y$ に至る鎖に $x$ を付け加えれば
$d(y,x)<\eps$ より $x$ への $\eps$ 鎖となり、$x\notin C$ に反する。
従って $B(x,\eps)\subset X\setminus C$ である。

$C$ は空でない開かつ閉な集合なので、:ref[prop-top-connected-image] により $C=X$ である。
:::

:::example[鎖があっても連結とは限らない]
逆は成り立たない。
$\Q$ で有理数 $a<b$ と $\eps>0$ を与え、自然数 $m>(b-a)/\eps$ を選ぶと、

$$
x_j=a+\frac jm(b-a)\qquad(0\le j\le m)
$$

はすべて有理数で、隣接距離 $(b-a)/m$ は $\eps$ 未満である。
例えば $a=0$、$b=1$、$\eps=1/4$ なら $m=5$ で歩幅は $1/5$ である。
それでも $\Q$ は連結でない。
有限個の点を小刻みに並べられることと、その間を埋める連結な集合があることとは別である。

一方 $X=[-1,0]\cup[1,2]$ では、$-1$ から $2$ への $\eps$ 鎖は $\eps>1$ でなければ存在しない。
初めて $[1,2]$ に入る一歩は $0$ 以下の点から $1$ 以上の点への移動であり、歩幅が $1$ 以上だからである。
:::

## 連結成分

空間全体が連結でないときは、各点の周りで「つながっている最大の部分」に分けて調べる。

:::definition[連結成分]{#def-top-component}
位相空間 $X$ の点 $x$ を含む連結部分集合すべての和を $C(x)$ と書き、
$x$ の:term[連結成分]{en="connected component" yomi="れんけつせいぶん"}（または単に成分）という。
すべての連結成分が一点である空間を :term[全不連結]{en="totally disconnected" yomi="ぜんふれんけつ"} という。
:::

:::proposition[成分による分割]{#prop-top-component-closed}
1. $C(x)$ は $x$ を含む連結部分集合であり、$x$ を含む連結部分集合はすべて $C(x)$ に含まれる（最大性）。
2. $C(x)$ は閉集合である。
3. 二つの成分 $C(x),C(y)$ は一致するか交わらない。従って成分全体は $X$ の分割をなし、$y\in C(x)$ と $C(x)=C(y)$ は同値である。
4. $X$ の空でない連結部分集合は、ちょうど一つの成分に含まれる。
:::

:::proof
(1) $x$ を含む連結集合の族は $\{x\}$ を含むので空でなく、すべて $x$ を共有する。
:ref[prop-top-connected-union] (1) により和 $C(x)$ は連結である。
最大性は定義から明らかである。

(2) $\overline{C(x)}$ は :ref[prop-top-connected-union] (3) により連結で $x$ を含むので、
(1) の最大性から $\overline{C(x)}\subset C(x)$ となり、$C(x)$ は閉である。

(3) $z\in C(x)\cap C(y)$ とすると、$C(x)\cup C(y)$ は $z$ を共有する二つの連結集合の和なので連結であり、
$x$ を含むので $C(x)$ に含まれ、$y$ を含むので $C(y)$ に含まれる。
よって $C(x)=C(x)\cup C(y)=C(y)$ である。
各点 $x$ は $C(x)$ に属するので、成分全体は $X$ を覆う。
$y\in C(x)$ なら $y\in C(x)\cap C(y)$ より $C(x)=C(y)$ であり、逆は $y\in C(y)$ から従う。

(4) 空でない連結集合 $A$ の一点 $a$ を取ると、(1) から $A\subset C(a)$ である。
$A$ が成分 $C(y)$ にも含まれれば $a\in C(a)\cap C(y)$ より $C(y)=C(a)$ である。
:::

$X$ が空でないとき、$X$ が連結であることと成分がただ一つであることは同値である。
成分は常に閉だが、開とは限らない。
ただし成分が有限個なら、各成分の補集合は残りの有限個の閉集合の和なので閉であり、各成分は開でもある。

:::example[成分が開にならない例]
$\Q$ の相異なる二点 $p<q$ の間には無理数 $r$ がある（[実数の連続性](/math/analysis/real-numbers/)）。
$\Q\cap(-\infty,r)$ と $\Q\cap(r,\infty)$ は $\Q$ の分離であり、$p,q$ を別々の側に含む。
:ref[lem-top-connected-in-separation] より $p,q$ を含む連結集合はないので、$\Q$ の成分はすべて一点である。
しかし一点は $\Q$ の開集合でない（どの開区間も無限個の有理数を含む）。
従って $\Q$ は全不連結だが離散ではない。

また $A=\{0\}\cup\{1/n:n\ge1\}\subset\R$ でも成分はすべて一点である。
$1/n$ の成分 $\{1/n\}$ は $A$ で開である（$1/n$ の十分小さい近傍は $A$ の他の点を含まない）が、
$0$ の成分 $\{0\}$ は開でない。
:::

:::example[補有限位相の空間]
無限集合 $X$ に補有限位相を入れると、空でない二つの開集合 $U,V$ はそれぞれ有限個の点しか除かないので、
$U\cap V$ は無限集合であり空でない。
従って分離は存在せず $X$ は連結である。
無限部分集合 $Y\subset X$ に誘導される位相も補有限位相なので $Y$ も連結である。
一方、二点以上の有限部分集合 $F$ には離散位相が誘導される（各点の補集合が有限で開だから）ので連結でない。
:::

:::example[整数点を除いた直線]
$X=\R\setminus\Z$ の成分は開区間 $(n,n+1)$（$n\in\Z$）である。
各区間は連結であり、$X\cap(n,n+1)$ は $X$ で開、その補集合 $X\cap\bigl((-\infty,n)\cup(n+1,\infty)\bigr)$ も開なので開かつ閉である。
$(n,n+1)$ を含む連結集合は :ref[lem-top-connected-in-separation] によりこの開かつ閉な集合からはみ出せないので、最大性から成分になる。
:::

### 連続写像と成分

:::proposition[連続写像による成分の移り方]{#prop-top-component-map}
$f:X\to Y$ が連続なら、すべての $x\in X$ で $f(C(x))\subset C(f(x))$ である。
$f$ が同相写像なら等号が成り立ち、$f$ は $X$ の成分全体と $Y$ の成分全体の間の全単射を引き起こす。
:::

:::proof
$f(C(x))$ は連結集合の連続像として連結であり、$f(x)$ を含むので、最大性から $C(f(x))$ に含まれる。
$f$ が同相なら逆写像 $f^{-1}$ に同じ議論を使うと $f^{-1}(C(f(x)))\subset C(x)$、
すなわち $C(f(x))\subset f(C(x))$ となる。
従って $f$ は成分を成分に写し、異なる成分は（交わらないので）異なる成分に写る。
全射性は $C(y)=f(C(f^{-1}(y)))$ から従う。
:::

連続全射でも、成分の像が値域の成分より真に小さい場合がある。
実数に離散位相を入れた空間から通常の $\R$ への恒等写像は連続全射だが、
定義域の成分は一点、値域の成分は $\R$ 全体である。

:::example[成分と準成分が異なる平面の部分集合]{#exm-top-quasicomponent-teeth}
点を含む開かつ閉な集合すべての共通部分を、その点の準成分 $Q(x)$ と呼ぶ。成分 $C(x)$ はそこに含まれる。各開閉集合とその補集合にまたがる連結集合は存在しないからである。包含が真となる例を調べよう。

平面の部分空間

$$
X=\bigcup_{n\ge1}\{1/n\}\times[0,1]\ \cup\ \{(0,0),(0,1)\}
$$

を考える。各線分 $I_n=\{1/n\}\times[0,1]$ は $X$ で開かつ閉である（$I_n$ は閉集合 $\{1/n\}\times[0,1]$ との交わりで、各 $1/n$ を他の第一座標から隔てる小開区間と $\R$ の積との交わりでもある）。
$(0,0)$ を含む連結集合 $D$ がある $I_n$ と交われば、$I_n$ が開かつ閉なので $D\subset I_n$ となり矛盾するから、$D\subset\{(0,0),(0,1)\}$ である。この二点集合は離散なので $C((0,0))=\{(0,0)\}$ である。
一方、$(0,0)$ を含む開かつ閉な集合 $H$ は開なので、ある $N$ 以上の $n$ について $(1/n,0)\in H$ となる。$H\cap I_n$ は連結な $I_n$ の空でない開かつ閉な部分集合なので $I_n$ 全体であり、$(1/n,1)\in H$（$n\ge N$）となる。$H$ は閉なので極限 $(0,1)$ も $H$ に属する。
従って $Q((0,0))=\{(0,0),(0,1)\}\ne C((0,0))$ である（$I_n$ の補集合は開かつ閉なので、$Q((0,0))$ に線分の点は入らない）。例えば $(1/n,1/2)$ の平面内の極限 $(0,1/2)$ は $X$ に含まれない。
:::

### 切断点と同相の判別

一点を除いたときに空間がいくつに分かれるかは、同相で保たれる量である。
これを使うと、連結性だけでは区別できない空間を見分けられる。

:::definition[切断点]{#def-top-cut-point}
連結空間 $X$ の点 $p$ について、$X\setminus\{p\}$ が連結でないとき $p$ を :term[切断点]{en="cut point" yomi="せつだんてん"} という。
:::

:::proposition[切断点の位相不変性]{#prop-top-cut-point-invariant}
$h:X\to Y$ を同相写像とし、$p\in X$ とする。
$h$ の制限は同相 $X\setminus\{p\}\to Y\setminus\{h(p)\}$ を与える。
従って $X\setminus\{p\}$ と $Y\setminus\{h(p)\}$ の成分の個数は等しい。
特に $h$ は切断点を切断点に、切断点でない点を切断点でない点に写し、
両者の個数（濃度）は $X$ と $Y$ で等しい。
:::

:::proof
$h$ は全単射なので $h(X\setminus\{p\})=Y\setminus\{h(p)\}$ である。
制限 $X\setminus\{p\}\to Y$ は連続であり、値が部分空間 $Y\setminus\{h(p)\}$ に入るので、
:ref[prop-top-subspace-laws] によりこの部分空間への写像として連続である。
逆写像は $h^{-1}$ の $Y\setminus\{h(p)\}$ への制限であり、同じ理由で連続である。
成分の個数については :ref[prop-top-component-map] を適用すればよい。
:::

:::example[三種類の区間]{#exm-top-three-intervals}
区間 $I$ の点 $p$ が切断点であるのは、$p$ が $I$ の最小元でも最大元でもないときに限る。
実際、そのような $p$ では $I\setminus\{p\}$ が $p$ の両側の点を含み凸でないので、:ref[thm-top-interval-connected] から連結でない。
逆に $p$ が最小元（または最大元）なら $I\setminus\{p\}$ は凸であり連結である。
従って切断点でない点の個数は

$$
[0,1]\colon 2,\qquad [0,1)\colon 1,\qquad (0,1)\colon 0
$$

であり、:ref[prop-top-cut-point-invariant] により三つの区間は互いに同相でない。
どれも連結なので、連結性だけではこの区別はできない。
:::

:::example[直線と平面]{#exm-top-line-plane}
$n\ge2$ のとき $\R^n$ から一点 $p$ を除いた空間は連結である。
証明のため、二点 $x,y$ を結ぶ線分 $[x,y]=\{(1-t)x+ty:0\le t\le1\}$ が $[0,1]$ の連続像として連結であることに注意する。

点 $x_0\ne p$ を固定し、各 $y\ne p$ に対して $x_0$ と $y$ を含む連結集合 $A_y\subset\R^n\setminus\{p\}$ を作る。
$p\notin[x_0,y]$ なら $A_y=[x_0,y]$ とする。
$p\in[x_0,y]$ なら $x_0\ne y$ であり、$x_0,y$ を通る直線 $\ell$ を考える。
$n\ge2$ なので $\ell$ 上にない点 $z$ がある。
線分 $[x_0,z]$ の点 $(1-t)x_0+tz$ が $t>0$ で $\ell$ に属すれば、$z=\bigl(((1-t)x_0+tz)-(1-t)x_0\bigr)/t$ は $\ell$ の二点のアフィン結合として $\ell$ に属してしまう。
従って $[x_0,z]$ と $\ell$ の共通点は $x_0$ だけであり、同様に $[z,y]$ と $\ell$ の共通点は $y$ だけである。
$p\in\ell$、$p\ne x_0,y$ なので $A_y=[x_0,z]\cup[z,y]$ は $p$ を含まず、$z$ を共有するので連結である。

すべての $A_y$ は $x_0$ を含み、和は $\R^n\setminus\{p\}$ に一致するので、
:ref[prop-top-connected-union] (1) により $\R^n\setminus\{p\}$ は連結である。
一方 $\R$ のすべての点は切断点である。
従って $n\ge2$ のとき $\R$ と $\R^n$ は同相でない。
$\R^2$ と $\R^3$ のように、どちらにも切断点がない場合の区別にはホモロジー論などのより精密な道具が必要になる。
:::

:::example[円周と閉区間]
円周 $S^1=\{(\cos\theta,\sin\theta)\}$ から点 $p=(\cos\theta_0,\sin\theta_0)$ を除いた集合は、
開区間 $(\theta_0,\theta_0+2\pi)$ の連続写像 $\theta\mapsto(\cos\theta,\sin\theta)$ による像なので連結である。
従って $S^1$ に切断点はないが、$[0,1]$ の内点はすべて切断点である。
よって $S^1$ と $[0,1]$ は同相でない。
:::

### 積の連結性

積空間の連結性を調べるため、まず積の中の「一部の座標だけを動かした切片」を確認する。

:::lemma[積の切片]{#lem-top-product-slice}
$P=\prod_{i\in I}X_i$ の点 $p$ と部分集合 $J\subset I$ を取り、

$$
E_J=\{x\in P:\ i\notin J\ \text{ならば}\ x_i=p_i\}
$$

と置く。
写像 $\iota_J:\prod_{i\in J}X_i\to P$ を、$i\in J$ では $y_i$、$i\notin J$ では $p_i$ を座標とする点に $(y_i)_{i\in J}$ を写すものとすると、
$\iota_J$ は $\prod_{i\in J}X_i$ から部分空間 $E_J$ への同相写像である。
:::

:::proof
$\iota_J$ の各座標は、$i\in J$ では射影 $(y_j)_{j\in J}\mapsto y_i$、$i\notin J$ では定数写像なので連続であり、
:ref[thm-top-product-property] により $\iota_J$ は連続である。
$\iota_J$ は単射で像は $E_J$ である。
逆写像は射影 $P\to\prod_{i\in J}X_i$、$x\mapsto(x_i)_{i\in J}$ の $E_J$ への制限であり、
この射影は各座標が $\pi_i$ なので再び :ref[thm-top-product-property] により連続である。
:::

例えば $X\times Y$ で $J=\{1\}$ とすれば、横の切片 $X\times\{q\}$ が $X$ と同相であることを表す。

:::theorem[積の連結性]{#thm-top-connected-product}
連結空間の族 $(X_i)_{i\in I}$ の積 $P=\prod_iX_i$ は連結である。
逆に空でない積 $P$ が連結なら、各因子 $X_i$ は連結である。
:::

:::proof
ある因子が空なら $P$ は空で連結であり、$I=\varnothing$ なら $P$ は一点で連結である。
以下すべての因子が空でないとする。

**二因子の場合。** $X,Y$ を空でない連結空間とし、$a\in X$ を固定する。
各 $y\in Y$ に対し $T_y=(\{a\}\times Y)\cup(X\times\{y\})$ と置く。
:ref[lem-top-product-slice] により $\{a\}\times Y\cong Y$、$X\times\{y\}\cong X$ は連結で、
点 $(a,y)$ を共有するので $T_y$ は連結である。
すべての $T_y$ は $\{a\}\times Y$ を含み、特に一点 $(a,y_0)$ を共有する。
任意の $(x,y)$ は $X\times\{y\}\subset T_y$ に属するので $\bigcup_yT_y=X\times Y$ であり、
:ref[prop-top-connected-union] (1) により $X\times Y$ は連結である。

**有限積。** 有限集合 $J$ の大きさに関する帰納法で $\prod_{i\in J}X_i$ が連結であることを示す。
$k\in J$ を一つ選び $J'=J\setminus\{k\}$ とする。
写像 $(x_i)_{i\in J}\mapsto\bigl((x_i)_{i\in J'},x_k\bigr)$ は全単射で、両方向とも各座標が射影なので
:ref[thm-top-product-property] により連続であり、同相 $\prod_{i\in J}X_i\cong\bigl(\prod_{i\in J'}X_i\bigr)\times X_k$ を与える。
帰納法の仮定と二因子の場合から右辺は連結である。

**一般の積。** 基点 $p\in P$ を固定し、有限個を除く座標が $p$ と一致する点全体を

$$
E=\bigcup_{J\subset I\ \text{有限}}E_J
$$

とする。
有限な $J$ ごとに $E_J\cong\prod_{i\in J}X_i$（:ref[lem-top-product-slice]）は有限積として連結で、
すべての $E_J$ は $p$ を含むので、$E$ は連結である。

$E$ は $P$ で稠密である。
空でない基本開集合 $B=\bigcap_{i\in J}\pi_i^{-1}(U_i)$（$J$ は有限、$U_i\ne\varnothing$）を取り、
$i\in J$ では $y_i\in U_i$、$i\notin J$ では $y_i=p_i$ とした点 $y$ を作ると、$y\in B\cap E_J\subset B\cap E$ である。
従って $\overline E=P$ であり、:ref[prop-top-connected-union] (3) により $P$ は連結である。

**逆。** $P$ が空でないとき、各射影 $\pi_i$ は全射である（基点の第 $i$ 座標だけを任意の値に取り替えればよい）。
連結空間の連続像として $X_i=\pi_i(P)$ は連結である。
:::

この証明で、積位相の基本開集合が有限個の座標しか制限しないことが本質的に使われた。
箱位相では結論が成り立たない（演習問題）。

:::proposition[積空間の成分]{#prop-top-product-components}
積空間 $P=\prod_iX_i$ で点 $x=(x_i)$ を含む成分は、各因子の成分の積 $\prod_iC_i(x_i)$ に等しい。
ここで $C_i(x_i)$ は $X_i$ における $x_i$ の成分である。
:::

:::proof
各射影 $\pi_i$ は連続なので、:ref[prop-top-component-map] により $\pi_i(C(x))\subset C_i(x_i)$ であり、
$C(x)\subset\prod_iC_i(x_i)$ を得る。

逆向きのため、$A=\prod_iC_i(x_i)$ の相対位相が各 $C_i(x_i)$ の積位相と一致することを確かめる。
相対位相の準開基は $A\cap\pi_i^{-1}(U_i)$（$U_i$ は $X_i$ で開）であり、
これは $\prod_jA_j$ のうち第 $i$ 座標が $C_i(x_i)\cap U_i$ に入る点全体、
すなわち積位相の準開基の元そのものである。
従って :ref[thm-top-connected-product] により $A$ は連結で、$x$ を含むので $A\subset C(x)$ である。
:::

:::example[縦線が成分になる積]
$\Q\times[0,1]$ では、$\Q$ の成分は一点、$[0,1]$ は連結なので、
点 $(q,t)$ の成分は縦線分 $\{q\}\times[0,1]$ である。
各縦線分は（成分なので）閉だが開でない。
積の任意の基本近傍 $(U\cap\Q)\times V$ は、$q$ 以外の有理数を第一座標にもつ点を含むからである。
また $\{0,1\}^{\N}$ は各因子が二点離散で全不連結なので、成分はすべて一点であり全不連結である。
:::

### Cantor 集合

全不連結な空間の重要な例が Cantor 集合である。
閉区間 $[0,1]$ を三等分して中央の開区間を除き、残った各閉区間で同じ操作を繰り返す。
これを二進の「左右の選択」で記述し、$\{0,1\}^{\N}$ との同相を証明する。

:::definition[Cantor 集合]
$0,1$ からなる長さ $n$ の語 $s=(s_1,\ldots,s_n)$ に対し

$$
l_s=\sum_{k=1}^n\frac{2s_k}{3^k},\qquad I_s=[\,l_s,\ l_s+3^{-n}\,]
$$

とする（空の語では $I_\varnothing=[0,1]$）。
$K_n=\bigcup_{|s|=n}I_s$、$K=\bigcap_{n\ge0}K_n$ と置き、$K$ を :term[Cantor 集合]{en="Cantor set" yomi="かんとーるしゅうごう"} という。
:::

長さ $n$ の語 $s$ について $l_{s0}=l_s$、$l_{s1}=l_s+2\cdot3^{-n-1}$ なので

$$
I_s=I_{s0}\cup\bigl(l_s+3^{-n-1},\ l_s+2\cdot3^{-n-1}\bigr)\cup I_{s1}
$$

であり、$K_{n+1}$ はちょうど $K_n$ の各閉区間から中央の三分の一の開区間を除いたものである。
特に $I_{t}\subset I_s$（$t$ が $s$ の延長のとき）であり、$K_n$ は $2^n$ 個の長さ $3^{-n}$ の閉区間からなる。

:::lemma[段階区間の間隔]{#lem-top-cantor-gap}
長さ $n$ の相異なる語 $s,t$ について、$I_s$ と $I_t$ の距離は $3^{-n}$ 以上である。
特に $I_s$ と $I_t$ は交わらない。
:::

:::proof
$s,t$ が初めて異なる位置を $k\le n$ とし、共通の頭部を $r=(s_1,\ldots,s_{k-1})$ とする。
必要なら入れ替えて $s_k=0$、$t_k=1$ とする。
$I_s\subset I_{r0}=[l_r,\ l_r+3^{-k}]$、$I_t\subset I_{r1}=[l_r+2\cdot3^{-k},\ l_r+3^{-k+1}]$ であり、
二つの閉区間の間には長さ $3^{-k}$ の開区間がある。
従って距離は $3^{-k}\ge3^{-n}$ 以上である。
:::

:::proposition[Cantor 集合と二進列空間]{#prop-top-cantor-set}
$\{0,1\}$ に離散位相を入れ、$\{0,1\}^{\N}$（$\N=\{1,2,\ldots\}$）に積位相を入れる。

$$
\varphi(a)=\sum_{k=1}^\infty\frac{2a_k}{3^k}
$$

は $\{0,1\}^{\N}$ から Cantor 集合 $K$ への同相写像である。
従って $K$ は全不連結であり、孤立点をもたない。
:::

:::proof
$a\in\{0,1\}^{\N}$ の初めの $n$ 項からなる語を $a|_n$ と書く。

**像が $K$ に入る。** $\varphi(a)=l_{a|_n}+\sum_{k>n}2a_k3^{-k}$ であり、尾部は

$$
0\le\sum_{k>n}\frac{2a_k}{3^k}\le\sum_{k>n}\frac{2}{3^k}=3^{-n}
$$

を満たすので $\varphi(a)\in I_{a|_n}\subset K_n$ である。$n$ は任意なので $\varphi(a)\in K$ である。

**単射。** $a\ne b$ が初めて異なる位置を $k$ とすると、$\varphi(a)\in I_{a|_k}$、$\varphi(b)\in I_{b|_k}$ で、
:ref[lem-top-cantor-gap] よりこの二区間は交わらない。

**全射。** $x\in K$ とする。
各 $n$ で $x\in K_n$ なので、$x\in I_s$ となる長さ $n$ の語 $s$ がある。
:ref[lem-top-cantor-gap] よりそれはただ一つであり、$s^{(n)}$ と書く。
$x\in I_{s^{(n+1)}}\subset I_{s^{(n+1)}|_n}$ と一意性から $s^{(n+1)}|_n=s^{(n)}$ である。
従って $a_k=s^{(n)}_k$（$n\ge k$）は $n$ によらずに定まり、$a|_n=s^{(n)}$ である。
$x$ と $\varphi(a)$ はともに長さ $3^{-n}$ の区間 $I_{a|_n}$ に属するので $|x-\varphi(a)|\le3^{-n}$ であり、
$n$ は任意だから $x=\varphi(a)$ である。

**連続性。** $a,b$ の初めの $N$ 項が一致すれば

$$
|\varphi(a)-\varphi(b)|\le\sum_{k>N}\frac{2|a_k-b_k|}{3^k}\le3^{-N}
$$

である。$\eps>0$ に対し $3^{-N}<\eps$ となる $N$ を選ぶと、
$a$ の基本近傍 $\{b:b_k=a_k\ (k\le N)\}$ は $\varphi$ によって $(\varphi(a)-\eps,\varphi(a)+\eps)$ に写る。

**開写像であること。** 長さ $n$ の語 $s$ に対し柱状集合 $Z_s=\{a:a|_n=s\}$ を考える。
$\{0,1\}^{\N}$ の基本開集合は有限個の座標の値を指定する集合であり、
指定された座標の最大番号を $n$ とすれば、条件に合う長さ $n$ の語 $s$ にわたる $Z_s$ の和である。
従って柱状集合全体は開基であり、各 $\varphi(Z_s)$ が $K$ の相対開集合であることを示せばよい。

まず $\varphi(Z_s)=K\cap I_s$ である。
$\subset$ は最初の段落による。
$x\in K\cap I_s$ なら、全射性の証明で作った $a$ は $a|_n=s^{(n)}=s$ を満たすので $x=\varphi(a)\in\varphi(Z_s)$ である。

次に $O=(l_s-3^{-n}/2,\ l_s+3^{-n}+3^{-n}/2)$ と置く。
:ref[lem-top-cantor-gap] より、長さ $n$ の他の語 $t$ の区間 $I_t$ の点は $I_s$ から $3^{-n}$ 以上離れているので $O$ に属さない。
$K\subset K_n$ だから $K\cap O=K\cap I_s$ であり、これは $K$ の相対開集合である。

以上より $\varphi$ は連続な開全単射であり、:ref[prop-top-homeomorphism-test] により同相写像である。

**帰結。** :ref[prop-top-product-components] により $\{0,1\}^{\N}$ の成分は一点なので、$K$ は全不連結である。
また各柱状集合 $Z_s$ は無限個の点を含むので $\{0,1\}^{\N}$ に孤立点はなく、$K$ にもない。
:::

$K_n$ の長さの総和は $(2/3)^n\to0$ だが、$K$ は $\{0,1\}^{\N}$ と全単射で結ばれるので非可算である。
「長さがない」「孤立点がない」「全不連結」という性質が両立する例として、Cantor 集合は位相・解析のいたるところに現れる。

## 弧状連結性

連結性は「分けられない」という否定的な定義であった。
より直観に近い「任意の二点を連続な曲線で結べる」という性質を定義し、両者を比較する。

:::definition[道と弧状連結性]{#def-top-path}
位相空間 $X$ の :term[道]{en="path" yomi="みち"} とは連続写像 $\gamma:[0,1]\to X$ であり、
$\gamma(0)$ を始点、$\gamma(1)$ を終点という。
任意の二点 $x,y$ に対し $x$ を始点、$y$ を終点とする道が存在するとき、
$X$ は:term[弧状連結]{en="path-connected" yomi="こじょうれんけつ"}であるという。
点 $x$ と道で結べる点全体を $x$ の :term[弧状連結成分]{en="path component" yomi="こじょうれんけつせいぶん"} といい $P(x)$ と書く。
本章では道の単射性は要求しない。
:::

道 $\gamma$ に対し、逆向きの道 $\gamma^{-}(t)=\gamma(1-t)$ を定める。
$\gamma(1)=\eta(0)$ を満たす道 $\gamma,\eta$ に対し、連結 $\gamma\ast\eta$ を

$$
(\gamma\ast\eta)(t)=
\begin{cases}
\gamma(2t) & (0\le t\le 1/2),\\
\eta(2t-1) & (1/2\le t\le 1)
\end{cases}
$$

で定める。
$t=1/2$ では両式の値は $\gamma(1)=\eta(0)$ で一致し、二つの閉集合 $[0,1/2]$、$[1/2,1]$ 上でそれぞれ連続なので、
:ref[lem-top-pasting] により $\gamma\ast\eta$ は連続である。

:::proposition[道の合成と成分]{#prop-top-path-component}
1. 「$x$ と $y$ を道で結べる」は $X$ 上の同値関係であり、その同値類 $P(x)$ は弧状連結である。
2. 弧状連結空間は連結である。さらに $P(x)\subset C(x)$ である。
3. 弧状連結空間の連続像は弧状連結である。
:::

:::proof
(1) 定数道 $\gamma(t)=x$ により反射律が、$\gamma^{-}$ により対称律が、$\gamma\ast\eta$ により推移律が成り立つ。
$y,z\in P(x)$ なら、推移律から $y$ と $z$ は道 $\gamma$ で $X$ 内で結べる。
その道の途中の点 $\gamma(t)$ も道 $s\mapsto\gamma(st)$ によって $y$ と結べるので $P(x)$ に属し、
$\gamma$ は $P(x)$ 内の道である。
従って $P(x)$ は弧状連結である。

(2) $X$ を弧状連結とし、$X\ne\varnothing$ として点 $x$ を固定する。
各 $y$ に対し $x$ から $y$ への道 $\gamma_y$ を取ると、像 $\gamma_y([0,1])$ は連結な $[0,1]$ の連続像なので連結であり、$x$ を含む。
これらの和は $X$ なので、:ref[prop-top-connected-union] (1) により $X$ は連結である。
(1) と合わせると、$P(x)$ は $x$ を含む連結集合なので $P(x)\subset C(x)$ である。

(3) $f:X\to Y$ が連続で $X$ が弧状連結なら、$f(X)$ の二点 $f(x),f(y)$ に対し、
$x$ から $y$ への道 $\gamma$ を取ると $f\circ\gamma$ が $f(x)$ から $f(y)$ への道である。
:::

:::example[凸集合・星形集合・積]
$\R^n$ の凸集合 $A$（$x,y\in A$ なら線分 $[x,y]\subset A$）では $\gamma(t)=(1-t)x+ty$ が道である。
ある点 $c$ からすべての点へ線分を引ける集合（$c$ に関する星形集合）では、$c$ を経由して $\gamma\ast\eta$ で二点を結べる。
例えば二本の座標軸の和 $\{xy=0\}$ は凸でないが原点に関する星形であり、弧状連結である。
弧状連結な空間の族 $(X_i)$ の積では、各座標で道 $\gamma_i$ を選ぶと
:ref[thm-top-product-property] により $t\mapsto(\gamma_i(t))_i$ が連続なので、積も弧状連結である。
:::

:::example[球面と穴あき空間]
$n\ge2$ のとき $\R^n\setminus\{p\}$ は弧状連結である。
:ref[exm-top-line-plane] で作った集合 $A_y$ は線分一本または二本の和であり、
線分をたどる道を $\ast$ でつなげば $x_0$ から $y$ への道が得られる。
$n\ge1$ のとき球面 $S^n=\{x\in\R^{n+1}:\norm x=1\}$ は弧状連結である。
実際、$x\mapsto x/\norm x$ は $\R^{n+1}\setminus\{0\}$ から $S^n$ への連続全射であり、
$n+1\ge2$ だから定義域は弧状連結なので、(3) により像も弧状連結である。
同様に円環 $\{x\in\R^2:1\le\norm x\le2\}$ は $S^1\times[1,2]$ の連続像 $(u,r)\mapsto ru$ として弧状連結である。
:::

:::proposition[一般線形群の成分]{#prop-top-gl-components}
$n\ge1$ とし、$\GL_n(\R)$ に $M_n(\R)\cong\R^{n^2}$ からの相対位相を入れる。

$$
\GL_n^{+}(\R)=\{A:\det A>0\},\qquad \GL_n^{-}(\R)=\{A:\det A<0\}
$$

はそれぞれ弧状連結であり、$\GL_n(\R)$ の連結成分（かつ弧状連結成分）はちょうどこの二つである。
:::

:::proof
行列式は成分の多項式なので連続である。
$\GL_n^{\pm}(\R)$ は開集合 $(0,\infty)$、$(-\infty,0)$ の逆像なので $\GL_n(\R)$ の開集合であり、
互いに交わらず、和が $\GL_n(\R)$ で、$I$ と $Q_n(-1)=\diag(1,\ldots,1,-1)$ をそれぞれ含むので空でない。
従って $\GL_n(\R)$ は連結でなく、:ref[lem-top-connected-in-separation] によりどの連結集合も一方に含まれる。
あとは各々が弧状連結であることを示せばよい。

基本行列 $R_{ij}(c)=I+cE_{ij}$（$i\ne j$、:ref[def-linalg-elementary-matrix]）を左から掛けると、第 $i$ 行に第 $j$ 行の $c$ 倍が加わる。
**主張：** 任意の $A\in\GL_n(\R)$ に対し、このような行列を有限個掛けて $R_mR_{m-1}\cdots R_1A=Q_n(d)=\diag(1,\ldots,1,d)$ とでき、$d=\det A$ である。

$n$ に関する帰納法で示す。$n=1$ なら $A=(d)$ でありそのままでよい。
$n\ge2$ とする。$A$ は正則なので第 $1$ 列は零でない。
$a_{21}=0$ なら、$a_{k1}\ne0$ となる $k\ne2$ を取り第 $k$ 行を第 $2$ 行に加えて $a_{21}\ne0$ とする。
次に第 $2$ 行の $(1-a_{11})/a_{21}$ 倍を第 $1$ 行に加えると $a_{11}=1$ となる。
さらに各 $i\ge2$ について第 $1$ 行の $-a_{i1}$ 倍を第 $i$ 行に加えると、第 $1$ 列は $e_1$ になる。
こうして $A$ は $\begin{pmatrix}1&*\\0&A'\end{pmatrix}$ の形になり、:ref[prop-linalg-det-operations] より行列式は変わらないので $\det A'=\det A\ne0$ である。
第 $2$ 行以降どうしの操作は第 $1$ 列の零を保つので、帰納法の仮定を $A'$ に適用して右下を $\diag(1,\ldots,1,d)$ にできる。
最後に第 $1$ 行の第 $j$ 列（$j\ge2$）の成分を、第 $j$ 行（$j<n$ では $\transp{e_j}$、$j=n$ では $d\,\transp{e_n}$）の適当な倍数を引いて消すと $Q_n(d)$ を得る。
各操作で行列式は変わらないので $d=\det A$ である。

$R_{ij}(c)^{-1}=R_{ij}(-c)$ だから、主張は $A=R_{i_1j_1}(c_1)\cdots R_{i_mj_m}(c_m)\,Q_n(d)$ と書けることを意味する。
$\sigma=\sgn d\in\{1,-1\}$ として

$$
\Gamma(t)=R_{i_1j_1}\bigl((1-t)c_1\bigr)\cdots R_{i_mj_m}\bigl((1-t)c_m\bigr)\,Q_n\bigl((1-t)d+t\sigma\bigr)\qquad(0\le t\le1)
$$

と置く。
各成分は $t$ の多項式なので $\Gamma$ は $M_n(\R)$ への連続写像であり、
:ref[thm-linalg-det-product] より $\det\Gamma(t)=(1-t)d+t\sigma$ は $d$ と同符号で零にならない。
従って $\Gamma$ は $\GL_n^{\pm}(\R)$ 内の道で、$A=\Gamma(0)$ と $\Gamma(1)=Q_n(\sigma)$ を結ぶ。
$\GL_n^{+}(\R)$ の任意の元は $I$ と、$\GL_n^{-}(\R)$ の任意の元は $Q_n(-1)$ と結べるので、両者は弧状連結である。
:::

:::example[二次の行列での計算]
$A=\begin{pmatrix}0&2\\1&3\end{pmatrix}$ は $\det A=-2$ である。
$a_{21}=1\ne0$ なので、第 $2$ 行の $1$ 倍を第 $1$ 行に加えて $\begin{pmatrix}1&5\\1&3\end{pmatrix}$、
第 $1$ 行を第 $2$ 行から引いて $\begin{pmatrix}1&5\\0&-2\end{pmatrix}$、
第 $2$ 行の $5/2$ 倍を第 $1$ 行に加えて $\diag(1,-2)$ を得る。
すなわち $R_{12}(5/2)R_{21}(-1)R_{12}(1)A=\diag(1,-2)$ であり、

$$
A=R_{12}(-1)\,R_{21}(1)\,R_{12}(-5/2)\,\diag(1,-2)
$$

である。
対応する道 $\Gamma(t)=R_{12}(-(1-t))R_{21}(1-t)R_{12}(-\tfrac52(1-t))\diag(1,t-2)$ の行列式は $t-2<0$ であり、
$\Gamma(0)=A$ と $\Gamma(1)=\diag(1,-1)$ を $\GL_2^{-}(\R)$ 内で結ぶ。
:::

### 位相幾何学者の正弦曲線

連結だが弧状連結でない空間の標準的な例を調べる。
$x\to0$ で $\sin(1/x)$ は $-1$ と $1$ の間を無限回振動するので、グラフは縦線分にいくらでも近づくが、
そこへ連続にたどり着くことはできない。

:::theorem[位相幾何学者の正弦曲線]{#thm-top-sine-curve}
平面で

$$
G=\{(x,\sin(1/x)):0<x\le1\},\qquad
S=\{0\}\times[-1,1],\qquad
T=G\cup S
$$

と置く。
$T=\overline G$ であり、$T$ は連結だが弧状連結でない。
$T$ の弧状連結成分はちょうど $G$ と $S$ の二つである。
:::

:::proof
**$G$ と $S$ は弧状連結。** $G$ は区間 $(0,1]$ の連続写像 $x\mapsto(x,\sin(1/x))$ による像であり、
区間は凸なので弧状連結だから、:ref[prop-top-path-component] (3) により $G$ も弧状連結である。
$S$ は線分なので弧状連結である。

**$\overline G=T$。**
$y\in[-1,1]$ に対し $\sin\theta=y$ となる $\theta\in[-\pi/2,\pi/2]$ を取り、$x_n=1/(\theta+2\pi n)$（$n\ge1$）と置く。
$\theta+2\pi n\ge2\pi-\pi/2>1$ なので $x_n\in(0,1)$ であり、$\sin(1/x_n)=\sin\theta=y$ である。
$(x_n,y)\in G$ は $(0,y)$ に収束するので $S\subset\overline G$ である。
逆に $(x,y)\in\overline G$ とする。
$G\subset[0,1]\times[-1,1]$ で右辺は閉なので $0\le x\le1$ である。
$x>0$ なら、$G$ の点 $(x_k,\sin(1/x_k))$ で $(x,y)$ に収束するものを取ると（:ref[prop-top-metric-closed-sequence]）、
$x_k\to x>0$ と $\sin(1/\cdot)$ の $x$ での連続性から $y=\sin(1/x)$ となり $(x,y)\in G$ である。
$x=0$ なら $|y|\le1$ より $(x,y)\in S$ である。

**連結性。** $T=\overline G$ は連結集合 $G$ の閉包なので、:ref[prop-top-connected-union] (3) により連結である。

**$S$ と $G$ を結ぶ道はない。** 道 $\gamma=(u,v):[0,1]\to T$ で $\gamma(0)\in S$、$\gamma(1)\in G$ となるものがあると仮定する。
$Z=u^{-1}(0)$ は $[0,1]$ の閉集合で $0$ を含むので、最大元 $s=\max Z$ をもつ。
$u(1)>0$ より $s<1$ であり、$t\in(s,1]$ では $u(t)>0$、従って $\gamma(t)\in G$、すなわち $v(t)=\sin(1/u(t))$ である。

各 $k\ge1$ に対し $t_k=s+(1-s)/k$、$c_k=u(t_k)>0$ と置く。
十分大きい自然数 $m$ を選べば

$$
a_k=\frac1{\pi/2+2\pi m},\qquad b_k=\frac1{3\pi/2+2\pi m}
$$

はともに $(0,c_k)$ に入り、$\sin(1/a_k)=1$、$\sin(1/b_k)=-1$ である。
$u$ は $[s,t_k]$ で連続で $u(s)=0<a_k,b_k<c_k=u(t_k)$ だから、
:ref[cor-top-intermediate-value] により $u(\tau_k)=a_k$、$u(\sigma_k)=b_k$ となる $\tau_k,\sigma_k\in(s,t_k)$ がある。
このとき $v(\tau_k)=1$、$v(\sigma_k)=-1$ である。
$\tau_k,\sigma_k\to s$ なので $v$ の連続性から $v(s)=\lim_kv(\tau_k)=1$ かつ $v(s)=\lim_kv(\sigma_k)=-1$ となり、矛盾する。
$G$ から $S$ への道があれば逆向きの道が $S$ から $G$ への道になるので、これも存在しない。

従って $G$ と $S$ はそれぞれ一つの弧状連結成分であり、$T$ は弧状連結でない。
:::

```tikz caption="位相幾何学者の正弦曲線。横軸を拡大し、グラフの一部と閉包に加わる縦線分を示す。"
\begin{tikzpicture}[x=5.5cm,y=.9cm]
\draw[->,gray] (-.08,0) -- (1.1,0) node[right] {$x$};
\draw[->,gray] (0,-1.35) -- (0,1.48) node[above] {$y$};
\draw[thick] (0,-1) -- (0,1);
\draw[thick] plot[domain=1:32,samples=200,variable=\t] ({1/\t},{sin(deg(\t))});
\node[left] at (0,1) {$1$};
\node[left] at (0,-1) {$-1$};
\node[below] at (1,0) {$1$};
\end{tikzpicture}
```

:::remark[閉包は弧状連結性を保たない]
:ref[prop-top-connected-union] (3) により連結集合の閉包は連結であったが、
弧状連結集合 $G$ の閉包 $T$ は弧状連結でない。
また $G\cup\{(0,0)\}$ のように $G$ と閉包の中間の集合も、連結だが弧状連結でない（演習問題）。
弧状連結性は「点列の極限を取る」操作と相性が悪いのである。
:::

:::example[区間の連続像にはなれない]
$T$ は区間からの連続全射の像になれない。
区間は弧状連結であり、:ref[prop-top-path-component] (3) から像も弧状連結でなければならないからである。
一方、$(0,1]$ と $[-1,1]$ の直和からの写像 $x\mapsto(x,\sin(1/x))$、$y\mapsto(0,y)$ は連続な全単射 $(0,1]\amalg[-1,1]\to T$ を与える。
しかし逆写像は連続でない。
直和では各成分が開であるのに、$S$ は $T$ で開でない（$S$ の点の近くに $G$ の点がある）からである。
「一つのパラメータで連続にたどれること」と「各部分を別々にパラメータ付けでき、和が連結であること」は異なる。
:::

### Euclid 空間の連結開集合

正弦曲線のような病理は、開集合では起こらない。

:::proposition[連結開集合と折れ線]{#prop-top-open-connected-polygonal}
$\Omega\subset\R^n$ を連結な開集合とすると、$\Omega$ の任意の二点は $\Omega$ 内の折れ線（有限個の線分の和）で結べる。
特に $\Omega$ は弧状連結である。
:::

:::proof
$\Omega=\varnothing$ なら明らかなので、$a\in\Omega$ を固定する。
$a$ から $\Omega$ 内の有限個の線分を順につないで到達できる点の集合を $P$ とする。
$a\in P$ である。

$P$ は開である。
$x\in P$ なら、$\Omega$ が開なので $B(x,r)\subset\Omega$ となる $r>0$ がある。
球は凸なので、$y\in B(x,r)$ に対して線分 $[x,y]\subset\Omega$ を最後に付け加えれば $y\in P$ となる。

$\Omega\setminus P$ も開である。
$x\in\Omega\setminus P$ で $B(x,r)\subset\Omega$ とする。
$B(x,r)$ が $P$ の点 $y$ を含めば、$y$ に至る折れ線に線分 $[y,x]\subset B(x,r)$ を付け加えて $x\in P$ となり矛盾する。
従って $B(x,r)\subset\Omega\setminus P$ である。

$P$ は $\Omega$ の空でない開かつ閉な集合なので、:ref[prop-top-connected-image] により $P=\Omega$ である。
折れ線は線分をたどる道の連結なので、$\Omega$ は弧状連結である。
:::

例えば $\R^2$ の開集合 $\{(x,y):1<x^2+y^2<4\}$ は連結な円環の内部であり、
任意の二点を有限個の線分（例えば中心を避けて回り込む三、四本の線分）で結べる。
一方、開集合でない円周 $S^1$ は弧状連結だが、二点を結ぶ折れ線は $S^1$ 内に一本もない。

### 位相群の単位元の成分

:::example[位相群の単位元の成分]
群 $G$ に位相が入り、積 $G\times G\to G$、$(x,y)\mapsto xy$ と逆元 $x\mapsto x^{-1}$ がともに連続であるとき、
$G$ を:term[位相群]{en="topological group" yomi="いそうぐん"}という。
単位元 $e$ の成分 $H=C(e)$ は閉正規部分群であり、各成分はその剰余類 $gH$ である。

実際、$H\times H$ は :ref[thm-top-connected-product] により連結で、
連続写像 $(x,y)\mapsto xy^{-1}$ による像は連結で $e=ee^{-1}$ を含むので $H$ に含まれる。
従って $x,y\in H$ なら $xy^{-1}\in H$ であり、$e\in H$ と合わせて $H$ は部分群である。
各 $g\in G$ に対し共役 $x\mapsto gxg^{-1}$ は連続で、逆写像 $x\mapsto g^{-1}xg$ も連続なので同相である。
これは $e$ を固定するので :ref[prop-top-component-map] により $gHg^{-1}=C(e)=H$ となり、$H$ は正規である。
左移動 $x\mapsto gx$ も同相で $e$ を $g$ に写すので $gH=C(g)$ である。
成分は閉なので $H$ は閉である。

乗法群 $\R\setminus\{0\}$ では $H=(0,\infty)$ であり、もう一つの成分は $(-1)H=(-\infty,0)$ である。
$\GL_n(\R)$ では :ref[prop-top-gl-components] により $H=\GL_n^{+}(\R)$ であり、$\GL_n^{-}(\R)=Q_n(-1)\,\GL_n^{+}(\R)$ である。
:::

### 中間値性による連結性の特徴付け

:::example[実数値関数で連結性を検査する]
空間 $X$ のすべての連続関数 $f:X\to\R$ が中間値性（像が $u<v$ を含めば $[u,v]$ を含む）をもつと仮定すると、$X$ は連結である。
非自明な開かつ閉な集合 $A$ があれば、$A$ 上で $0$、補集合で $1$ を取る関数は、
任意の開集合の逆像が $\varnothing,A,X\setminus A,X$ のいずれかなので連続である。
その像は $\{0,1\}$ で中間の値 $1/2$ を欠き、仮定に反する。
逆は :ref[cor-top-intermediate-value] である。

連続性が欠けると中間値性は失われる。
$\gamma:[0,1]\to\R$ を $t<1/2$ で $2$、$t\ge1/2$ で $-1$ とすると、値 $0$ を取らない。
この写像では開集合 $(-2,0)$ の逆像 $[1/2,1]$ が $[0,1]$ の相対開集合でなく、確かに連続でない。
:::

## 局所連結性

連結性は大域的な性質であり、各点の近くでの振る舞いについては何も言わない。
例えば二つの離れた開区間の和は連結でないが、各点の周りは区間のように見える。
逆に正弦曲線 $T$ は連結だが、縦線分の点の近くでは無数の断片に分かれている。
この局所的な振る舞いを捉えるのが局所連結性である。

:::definition[局所連結性と局所弧状連結性]{#def-top-locally-connected}
位相空間 $X$ が点 $x$ で局所連結であるとは、$x$ の任意の開近傍 $U$ に対し、
$x\in V\subset U$ を満たす連結な開集合 $V$ が存在することをいう。
すべての点で局所連結なとき、$X$ は:term[局所連結]{en="locally connected" yomi="きょくしょれんけつ"}であるという。
「連結」を「弧状連結」に置き換えた条件を :term[局所弧状連結]{en="locally path-connected" yomi="きょくしょこじょうれんけつ"} という。
:::

弧状連結な集合は連結なので、局所弧状連結な空間は局所連結である。
連結な開集合からなる開基をもつことは局所連結性と同値であり、
例えば $\R^n$ とその開集合は、開球が凸なので局所弧状連結である。

:::theorem[開集合の成分による判定]{#thm-top-local-components}
1. $X$ が局所連結であるための必要十分条件は、$X$ のすべての開集合 $U$ について、$U$（の相対位相）の連結成分がすべて $X$ で開であることである。
2. $X$ が局所弧状連結であるための必要十分条件は、すべての開集合の弧状連結成分が $X$ で開であることである。
3. $X$ が局所弧状連結なら、$X$ の連結成分と弧状連結成分は一致する。
:::

:::proof
(1) $X$ を局所連結とし、開集合 $U$ の成分 $C$ と $x\in C$ を取る。
局所連結性から $x\in V\subset U$ となる連結開集合 $V$ がある。
$V$ は $U$ の連結部分集合で $x$ を含むので、:ref[prop-top-component-closed] の最大性から $V\subset C$ である。
従って $C$ は各点の開近傍を含み、$X$ で開である。
逆に条件が成り立つとし、$x$ とその開近傍 $U$ を取る。
$U$ における $x$ の成分 $V$ は条件から $X$ で開であり、連結で $x\in V\subset U$ を満たす。

(2) (1) の証明で「連結」を「弧状連結」、成分を弧状連結成分に置き換えればよい。
$V\subset U$ が $x$ を含む弧状連結開集合なら、$V$ の各点は $V\subset U$ 内の道で $x$ と結べるので、
$U$ における $x$ の弧状連結成分に含まれる。

(3) $X$ を局所弧状連結とする。
(2) を $U=X$ に使うと、各弧状連結成分は開である。
各弧状連結成分の補集合は他の弧状連結成分の和なので開であり、各弧状連結成分は開かつ閉である。
連結成分 $C(x)$ は連結なので、:ref[lem-top-connected-in-separation] により開かつ閉な $P(x)$ と補集合のうち $x$ を含む側、すなわち $P(x)$ に含まれる。
:ref[prop-top-path-component] の $P(x)\subset C(x)$ と合わせて $P(x)=C(x)$ である。
:::

:::corollary[連結な局所弧状連結空間]{#cor-top-lpc-path-connected}
連結かつ局所弧状連結な空間は弧状連結である。
特に $\R^n$ の連結開集合、連結な位相多様体は弧状連結である。
:::

:::proof
:ref[thm-top-local-components] (3) より弧状連結成分は連結成分と一致し、連結空間の連結成分は全体ただ一つである。
:::

:::proposition[開部分空間の局所連結性]{#prop-top-local-open-subspace}
局所連結（局所弧状連結）な空間の開部分空間は局所連結（局所弧状連結）である。
:::

:::proof
$U\subset X$ を開とし、$x\in U$ と $U$ での $x$ の開近傍 $W$ を取る。
$U$ は $X$ で開なので $W$ も $X$ で開である。
局所連結性から $x\in V\subset W$ となる連結な $X$ の開集合 $V$ があり、$V\subset U$ なので $V$ は $U$ でも開である。
弧状連結の場合も同様である。
:::

開でない部分空間には局所連結性は受け継がれない。
$\R$ は局所連結だが、部分空間 $\Q$ や $\{0\}\cup\{1/n:n\ge1\}$ は局所連結でない。
後者では $0$ を含む連結集合は $\{0\}$ だけであり、これは開でないからである。

:::example[局所と大域の違い]
二つの離れた開円板の和は、各点に凸な小円板の近傍があるので局所弧状連結であるが、連結でない。
離散空間も一点集合が連結開集合なので局所連結であり、二点以上あれば連結でない。
一方、:ref[exm-top-long-ray] の長い半直線 $L$ では、開区間と開半直線 $\{x<b\}$ は順序凸であり、
:ref[prop-top-linear-continuum] によって連結である。
これらは順序位相の開基をなす（有限個の開半直線の共通部分は開区間か開半直線である）ので、$L$ は局所連結である。
$L$ は連結・局所連結で各点の近くは実数の区間と同じ順序構造をもつが、可分でない。
:::

:::example[正弦曲線の局所での障害]
:ref[thm-top-sine-curve] の $T$ は $(0,0)$ で局所連結でない。
$(0,0)$ の開近傍として帯との共通部分 $W=T\cap\bigl(\R\times(-1/2,1/2)\bigr)$ を取り、
$(0,0)\in V\subset W$ となる $T$ の連結な開集合 $V$ があると仮定する。
$V$ は $T$ で開なので、ある $r>0$ で $T\cap B((0,0),r)\subset V$ となる。
$m$ を大きく取れば $x_0=1/(\pi m)<r/2$ であり、$(x_0,0)\in G$ は $V$ に属する。
第一座標への射影 $\pi_1$ は連続なので $\pi_1(V)$ は $0$ と $x_0$ を含む $\R$ の連結集合、すなわち区間であり、$[0,x_0]$ を含む。
$c=1/(\pi/2+2\pi m')$ を $c<x_0$ となるように選ぶと、$V$ の点で第一座標が $c$ のものがある。
しかし $T$ で第一座標が $c>0$ の点は $(c,\sin(1/c))=(c,1)$ だけであり、これは $W$ に属さない。
矛盾するので、$W$ の中に $(0,0)$ を含む連結開集合は取れない。

$G$ の点では $T$ は局所弧状連結である（$G$ は $T$ で開であり、$(0,1]$ と同相だから）。
局所連結性が破れるのは縦線分 $S$ の上だけである。
:::

### 弧状連結でも局所連結でない櫛形空間

:::example[櫛形空間]
$D=\{0\}\cup\{1/n:n\ge1\}$ とし、

$$
H=([0,1]\times\{0\})\cup(D\times[0,1])
$$

を考える。
各点から自分の縦線（歯）を底辺まで降り、底辺に沿って原点へ進めるので、$H$ は弧状連結である。

一方 $H$ は $(0,1/2)$ で局所連結でない。
開近傍 $W=H\cap\bigl(\R\times(1/4,3/4)\bigr)$ を取り、$(0,1/2)\in V\subset W$ となる連結な開集合 $V$ があるとする。
$W$ の点の第一座標は $D$ に属するので、$\pi_1(V)$ は $D$ の連結部分集合である。
$D$ の相異なる二点の間には $D$ に属さない実数があるので、$D$ の連結部分集合は高々一点であり、$\pi_1(V)=\{0\}$ である。
ところが $V$ は $H$ で開なので $(0,1/2)$ のある近傍の $H$ の点をすべて含み、
十分大きい $n$ では $(1/n,1/2)\in V$ となって $\pi_1(V)\ni1/n$ となる。
これは矛盾である。

遠く離れた歯の点どうしを結ぶ道は必ず底辺を通るので、道は存在しても「近くの点を近くで結ぶ」ことができない。
弧状連結性は局所弧状連結性を、さらには局所連結性すら含意しない。
:::

### 局所連結でも道を引けない順序空間

:::example[辞書式順序を入れた正方形]
$L=[0,1]\times[0,1]$ を辞書式順序（先に第一座標、等しければ第二座標で比べる）で並べ、順序位相を入れる。
この位相は平面からの相対位相と異なる。
例えば縦線の内部 $\{s\}\times(0,1)$ は順序の開区間 $((s,0),(s,1))$ であり $L$ で開だが、平面では開でない。

**線形連続体であること。** $(s,t)<(s',t')$ のとき、$s<s'$ なら $((s+s')/2,t)$、$s=s'$ なら $(s,(t+t')/2)$ が間にある。
空でない $A\subset L$ の第一座標の上限を $s_0$ とする。
$s_0$ が $A$ の第一座標として現れるなら、$t_0=\sup\{t:(s_0,t)\in A\}$ として $(s_0,t_0)$ が上限であり、
現れないなら $(s_0,0)$ が上限である。
実際、どちらも上界であり、それより小さい点 $(p,q)$ を取ると、$p<s_0$ なら第一座標が $p$ より大きい $A$ の点が、
$p=s_0$ なら $q<t_0$ より第二座標が $q$ より大きい点 $(s_0,t)\in A$ がそれを超える。
従って :ref[prop-top-linear-continuum] により $L$ は連結であり、
開区間と端を含む半直線 $[(0,0),b)$、$(a,(1,1)]$ も順序凸で連結である。
これらは開基をなすので $L$ は局所連結である。

**第一座標の異なる二点は道で結べない。** 道の像 $P$ が第一座標 $a<b$ の点 $z_1<z_2$ を含むとする。
$P$ は連結である。
$z_1<z<z_2$ を満たす $z$ が $P$ に属さなければ、$P\cap\{x<z\}$ と $P\cap\{x>z\}$ が $P$ の分離になるので、$P$ は $z_1$ と $z_2$ の間の点をすべて含む。
特に $a<s<b$ を満たす各 $s$ について、空でない開区間 $J_s=\{s\}\times(0,1)$ は $P$ に含まれ、$P$ の相対開集合である。
$s\ne s'$ なら $J_s\cap J_{s'}=\varnothing$ である。

一方 $P$ は $[0,1]$ の連続像なので可分である。
実際、$\gamma:[0,1]\to P$ を全射とし、$\Q\cap[0,1]$ の像を $E$ とすると、
$P$ の空でない開集合 $W$ の逆像 $\gamma^{-1}(W)$ は $[0,1]$ の空でない開集合で有理数を含むので、$W\cap E\ne\varnothing$ である。
各 $J_s$ は $E$ の点を含み、$J_s$ たちは互いに交わらないので、$s\mapsto$（$J_s$ に含まれる $E$ の点の一つ）は非可算集合 $(a,b)$ から可算集合 $E$ への単射になる。
これは不可能である。

**弧状連結成分は縦線である。** 写像 $t\mapsto(s,t)$ は $[0,1]$ から縦線 $\{s\}\times[0,1]$ への同相である。
実際、開半直線 $\{x<(a,b)\}$ の逆像は $a>s$ なら $[0,1]$、$a<s$ なら $\varnothing$、$a=s$ なら $[0,b)$ であり、
$\{x>(a,b)\}$ についても同様なので連続である。
逆向きも、$[0,1]$ の相対開区間の像が $L$ の開区間または半直線と縦線との共通部分として表せるので連続である。
従って縦線は弧状連結であり、弧状連結成分はちょうど縦線 $\{s\}\times[0,1]$ である。

この空間は連結・局所連結だが弧状連結でない。
局所連結性と連結性を合わせても、局所「弧状」連結性がなければ道の存在は保証されない。
:::

### 局所連結性を商に移す

:::proposition[局所連結空間の商]{#prop-top-local-quotient}
$q:X\to Y$ を商写像とする。
$X$ が局所連結なら $Y$ も局所連結であり、$X$ が局所弧状連結なら $Y$ も局所弧状連結である。
:::

:::proof
:ref[thm-top-local-components] により、$Y$ の開集合 $V$ の成分 $C$ が $Y$ で開であることを示せばよい。
商位相の定義から、そのためには $q^{-1}(C)$ が $X$ で開であることを示せば十分である。

$q^{-1}(V)$ は $X$ の開集合なので、:ref[thm-top-local-components] によりその成分はすべて $X$ で開である。
$x\in q^{-1}(C)$ を取り、$q^{-1}(V)$ における $x$ の成分を $D$ とする。
$q(D)$ は連結集合の連続像として連結で、$V$ に含まれ、$q(x)\in C$ を含むので、
成分の最大性から $q(D)\subset C$、すなわち $D\subset q^{-1}(C)$ である。
従って $q^{-1}(C)$ は開集合 $D$ たちの和であり開である。

局所弧状連結の場合は、成分を弧状連結成分に置き換え、弧状連結集合の連続像が弧状連結であること
（:ref[prop-top-path-component] (3)）を使えば同じ議論が成り立つ。
:::

:::example[円周とトーラス]
$[0,1]$ は局所弧状連結であり、両端を同一視する商写像 $[0,1]\to S^1$ を通して $S^1$ も局所弧状連結である。
同様に、正方形の向かい合う辺の同一視による商として得られるトーラス $S^1\times S^1$（:ref[prop-top-torus]）も局所弧状連結である。
商写像を経由すれば、局所座標を具体的に書かなくても局所連結性が分かる。
:::

### 無限積が局所連結になる条件

連結性は任意積で保たれたが、局所連結性は無限積で失われうる。

:::example[二点離散空間の無限積]
$\{0,1\}^{\N}$ では各因子が離散なので局所弧状連結であるのに、積は局所連結でない。
:ref[thm-top-connected-product] の後の例で見たように成分はすべて一点であり、
積の基本開集合は有限個の座標しか制限しないので一点集合は開でない。
:ref[thm-top-local-components] (1) から局所連結でない。
:::

:::theorem[積の局所連結性]{#thm-top-locally-connected-product}
空でない積 $P=\prod_{i\in I}X_i$ が局所連結であるための必要十分条件は、
すべての因子 $X_i$ が局所連結で、しかも有限個を除くすべての因子が連結であることである。
:::

:::proof
**十分性。** 非連結な因子の添字は有限集合 $F$ をなすとする。
点 $p\in P$ とその開近傍 $W$ を取り、$W$ に含まれる基本開集合 $\bigcap_{i\in J_0}\pi_i^{-1}(U_i)\ni p$（$J_0$ は有限）を選ぶ。
$J=J_0\cup F$ とし、$i\in J\setminus J_0$ では $U_i=X_i$ と置く。
各 $i\in J$ で、$X_i$ の局所連結性から $p_i\in V_i\subset U_i$ となる連結開集合 $V_i$ を選ぶ。
$i\notin J$ では $V_i=X_i$ と置くと、これは連結である（$i\notin F$）。
すると

$$
V=\prod_{i\in I}V_i=\bigcap_{i\in J}\pi_i^{-1}(V_i)
$$

は $P$ の開集合で $p\in V\subset W$ を満たし、:ref[thm-top-connected-product] により連結空間の積として連結である
（部分空間としての位相が積位相に一致することは :ref[prop-top-product-components] の証明と同じである）。
従って $P$ は局所連結である。

**必要性。** $P$ は空でないので各射影 $\pi_i$ は全射であり、連続な開写像である（基本開集合の像は開集合なので）。
:ref[prop-top-surjective-test] により $\pi_i$ は商写像であり、:ref[prop-top-local-quotient] から各 $X_i$ は局所連結である。
次に一点 $p\in P$ を固定し、局所連結性を開近傍 $P$ に適用して $p$ を含む連結開集合 $W$ を取る。
$W$ は $p$ を含む基本開集合 $B=\bigcap_{i\in J}\pi_i^{-1}(U_i)$（$J$ は有限）を含む。
$i\notin J$ なら、$p$ の第 $i$ 座標だけを任意の $y\in X_i$ に取り替えた点が $B$ に属するので $\pi_i(B)=X_i$、従って $\pi_i(W)=X_i$ である。
$W$ は連結なので $X_i=\pi_i(W)$ は連結である。
よって非連結な因子の添字は $J$ に含まれ、有限個である。
:::

:::example[一つだけ非連結な因子を許す]
$\{0,1\}\times\R^{\N}$ は局所連結である。
例えば点 $p$ の基本近傍に対し、第一座標を $p$ の値に固定し、指定された有限個の実座標に開区間を選べば連結開近傍ができる。
一方 $\{0,1\}^{\N}\times\R$ は局所連結でない。
連結開集合 $W$ が含む基本開集合は制限していない離散座標 $i$ を残すので $\pi_i(W)=\{0,1\}$ となるが、
連結集合の像は連結なので一点でなければならないからである。
:::

## 演習問題

::::exercise{level=basic}
$X=(-2,-1]\cup[0,1]\cup\{3\}\subset\R$ の連結成分を求めよ。

:::solution
$A=(-2,-1]$、$B=[0,1]$、$C=\{3\}$ はそれぞれ区間なので連結である。
$X\cap(-\infty,-1/2)=A$、$X\cap(-1/2,2)=B$、$X\cap(2,\infty)=C$ はいずれも $X$ の相対開集合で、
互いに交わらず和が $X$ なので、各々は開かつ閉である（補集合が他の二つの和で開）。
:ref[lem-top-connected-in-separation] により、連結集合はこのうちの一つに含まれる。
$A,B,C$ はそれぞれ連結で、それを含む連結集合は同じ片側からはみ出せないので、最大性から成分はちょうど $A,B,C$ の三つである。
:::
::::

::::exercise{level=basic}
二点以上をもつ密着空間 $X$ は弧状連結であることを示せ。

:::solution
$x,y\in X$ に対し、$\gamma(t)=x$（$0\le t<1$）、$\gamma(1)=y$ と置く。
$X$ の開集合は $\varnothing$ と $X$ だけで、その逆像は $\varnothing$ と $[0,1]$ なので開であり、$\gamma$ は連続である。
従って $\gamma$ は $x$ から $y$ への道である。
密着空間ではどんな写像も連続なので、道の存在はほとんど自明になる。
:::
::::

::::exercise{level=basic}
$\R\setminus\{1,2,3\}$ の弧状連結成分を求めよ。

:::solution
$\R\setminus\{1,2,3\}=(-\infty,1)\cup(1,2)\cup(2,3)\cup(3,\infty)$ であり、各区間は凸なので弧状連結である。
異なる区間の点 $x<y$ を結ぶ道 $\gamma$ があれば、$\gamma([0,1])$ は連結な $\R$ の部分集合なので区間であり、$x,y$ の間にある除かれた点を含んでしまう。
従って弧状連結成分はこの四つの開区間である。
:::
::::

::::exercise{level=basic}
$[-1,1]$ と円周 $S^1$ が同相でないことを、一点を除く議論で示せ。

:::solution
$[-1,1]\setminus\{0\}=[-1,0)\cup(0,1]$ は凸でないので連結でない。
一方 $S^1$ からどの一点 $(\cos\theta_0,\sin\theta_0)$ を除いても、残りは開区間 $(\theta_0,\theta_0+2\pi)$ の連続像なので連結である。
同相 $h:[-1,1]\to S^1$ があれば、:ref[prop-top-cut-point-invariant] により $[-1,1]\setminus\{0\}$ と $S^1\setminus\{h(0)\}$ が同相になり、矛盾する。
:::
::::

::::exercise{level=basic}
連続写像 $f:[0,1]\to[0,1]$ は不動点（$f(c)=c$ となる $c$）をもつことを示せ。
また $(0,1)$ では同じ主張が成り立たないことを例で示せ。

:::solution
$g(x)=f(x)-x$ は連続で、$g(0)=f(0)\ge0$、$g(1)=f(1)-1\le0$ である。
$g(0)=0$ または $g(1)=0$ なら端点が不動点である。
そうでなければ $g(0)>0>g(1)$ であり、$[0,1]$ は連結なので :ref[cor-top-intermediate-value] により $g(c)=0$ となる $c$ がある。

$(0,1)$ 上の $f(x)=x/2$ は $(0,1)$ への連続写像で、$f(x)=x$ は $x=0$ だけを解にもつが $0\notin(0,1)$ なので不動点はない。
:::
::::

::::exercise{level=basic}
$X$ を連結空間、$f:X\to\R$ を連続とし、$f$ の値がすべて整数であるとする。
$f$ は定数であることを示せ。

:::solution
$f(X)$ は連結なので区間である。
相異なる二つの整数 $m<n$ を含む区間は $m+1/2$ を含むが、これは整数でない。
従って $f(X)$ は高々一点であり、$f$ は定数である（$X=\varnothing$ なら自明）。
:::
::::

::::exercise{level=basic}
平面の「T 字形」$Y=([-1,1]\times\{0\})\cup(\{0\}\times[0,1])$ は $[0,1]$ と同相でないことを示せ。

:::solution
$p=(0,0)$ を除くと、$Y\setminus\{p\}$ は $A=[-1,0)\times\{0\}$、$B=(0,1]\times\{0\}$、$E=\{0\}\times(0,1]$ の和である。
各々は区間と同相で連結である。
$\overline A=[-1,0]\times\{0\}$ は $Y\setminus\{p\}$ と $A$ でしか交わらないので $A$ は $Y\setminus\{p\}$ で閉であり、$B,E$ も同様である。
三つとも閉で互いの補集合が他の二つの和なので開でもあり、:ref[lem-top-connected-in-separation] により成分はちょうど $A,B,E$ の三つである。
一方 $[0,1]$ から一点を除いた空間の成分は高々二つである（:ref[exm-top-three-intervals]）。
:ref[prop-top-cut-point-invariant] より $Y$ と $[0,1]$ は同相でない。
:::
::::

::::exercise{level=advanced}
連結集合 $A,B\subset X$ が $A\cap\overline B\ne\varnothing$ を満たせば、$A\cup B$ は連結であることを示せ（$A\cap B=\varnothing$ でもよい）。

:::solution
$p\in A\cap\overline B$ を取る。
$B\subset B\cup\{p\}\subset\overline B$ なので、:ref[prop-top-connected-union] (3) により $B\cup\{p\}$ は連結である。
$A$ と $B\cup\{p\}$ は $p$ を共有するので、同じ命題の (1) により和 $A\cup B$ は連結である。
:::
::::

::::exercise{level=advanced}
連結空間 $X$ から $\Q$ への連続写像の像を決定せよ。

:::solution
像 $f(X)$ は連結である。
部分空間 $f(X)\subset\Q$ の位相は $\R$ の部分空間としての位相と一致するので、$f(X)$ は $\R$ の連結部分集合、すなわち区間で、しかも $\Q$ に含まれる。
二点以上を含む区間は無理数を含むので、$f(X)$ は高々一点である。
従って $X\ne\varnothing$ なら $f$ は定数、$X=\varnothing$ なら像は空である。
:::
::::

::::exercise{level=advanced}
Sorgenfrey 直線（半開区間 $[a,b)$ を開基とする $\R$）が全不連結であることを示し、離散空間でないことを確かめよ。

:::solution
$x<y$ とする。
$(-\infty,y)=\bigcup_n[y-n,y)$ と $[y,\infty)=\bigcup_n[y,y+n)$ は半開区間の和なので開であり、
交わらず和が全体なので Sorgenfrey 直線の分離を与える。
$x$ は前者、$y$ は後者に属するので、:ref[lem-top-connected-in-separation] により $x,y$ をともに含む連結集合はない。
従って成分はすべて一点である。
しかし $x$ の任意の基本近傍 $[x,x+\eps)$ は $x$ 以外の点を含むので、一点集合 $\{x\}$ は開でなく、離散空間ではない。
:::
::::

::::exercise{level=advanced}
$E=\{(x,\sin(1/x)):0<x\le1\}\cup\{(0,0)\}$ の連結成分と弧状連結成分を求めよ。

:::solution
$G=\{(x,\sin(1/x)):0<x\le1\}$ について $G\subset E\subset\overline G$ なので、:ref[prop-top-connected-union] (3) により $E$ は連結であり、成分は $E$ 自身である。
$E$ 内で原点から $G$ の点への道があれば、それは $T\supset E$ 内の道でもあり、:ref[thm-top-sine-curve] に反する。
$G$ は弧状連結なので、弧状連結成分は $G$ と $\{(0,0)\}$ の二つである。
:::
::::

::::exercise{level=advanced}
開集合 $U=\{(x,y)\in\R^2:xy\ne0\}$ の連結成分と弧状連結成分を求めよ。

:::solution
$U$ は四つの開象限 $Q_{\pm\pm}=\{\pm x>0,\ \pm y>0\}$ の和である。
各象限は凸なので弧状連結であり、$\R^2$ で開なので $U$ でも開である。
各象限の $U$ での補集合は他の三つの和なので開であり、各象限は $U$ で開かつ閉である。
:ref[lem-top-connected-in-separation] により連結集合は一つの象限に含まれるので、連結成分も弧状連結成分も四つの象限である。
:::
::::

::::exercise{level=advanced}
$X$ の点 $x$ について、$x$ を含む開かつ閉な集合すべての共通部分を $Q(x)$ とする（準成分）。
$C(x)\subset Q(x)$ を示し、$X=\Q$ では $Q(x)=C(x)=\{x\}$ であることを確かめよ。

:::solution
$A$ を $x$ を含む開かつ閉な集合とする。
$A$ と $X\setminus A$ は交わらない開集合で、$C(x)$ を覆うので、:ref[lem-top-connected-in-separation] により $C(x)$ は一方に含まれ、$x\in A$ より $C(x)\subset A$ である。
$A$ は任意なので $C(x)\subset Q(x)$ である。

$X=\Q$ で $y\ne x$ なら、$x,y$ の間の無理数 $r$ について $\Q\cap(-\infty,r)$ と $\Q\cap(r,\infty)$ は開かつ閉で、$x$ を含む方は $y$ を含まない。
従って $y\notin Q(x)$ であり、$Q(x)=\{x\}$ である。
$\{x\}\subset C(x)\subset Q(x)$ より三つは一致する。
コンパクト Hausdorff 空間では常に $C(x)=Q(x)$ となることを [コンパクト性](/math/topology/compactness/) で示す。
:::
::::

::::exercise{level=advanced}
連続関数 $f:S^1\to\R$ に対し、$f(x)=f(-x)$ を満たす $x\in S^1$ が存在することを示せ。

:::solution
$g(x)=f(x)-f(-x)$ は連続で $g(-x)=-g(x)$ を満たす。
$S^1$ は弧状連結なので連結である。
ある $x_0$ で $g(x_0)\ne0$ なら、$g(x_0)$ と $g(-x_0)=-g(x_0)$ は異符号なので、
:ref[cor-top-intermediate-value] により $g(x)=0$ となる $x$ がある。
$g(x_0)=0$ ならそのまま $x=x_0$ でよい。
例えば赤道上の気温が連続なら、ちょうど反対側と同じ気温の地点が存在する。
:::
::::

::::exercise{level=advanced}
連続な単射 $f:\R^2\to\R$ は存在しないことを示せ。

:::solution
$f$ を連続な単射と仮定する。
$I=f(\R^2)$ は連結なので区間であり、$f$ が単射なので二点以上を含む。
任意の点 $p\in\R^2$ について、$\R^2\setminus\{p\}$ は :ref[exm-top-line-plane] により連結なので、
$f(\R^2\setminus\{p\})=I\setminus\{f(p)\}$（単射性による）は連結、すなわち凸である。
$f(p)$ が $I$ の最小元でも最大元でもなければ $I\setminus\{f(p)\}$ は $f(p)$ の両側の点を含んで凸でないので、
$f(p)$ は $I$ の最小元か最大元である。
これはすべての $p$ について成り立つが、$I$ の最小元・最大元は合わせて高々二つであり、
相異なる三点 $p_1,p_2,p_3$ の像が相異なること（単射性）に反する。
:::
::::

::::exercise{level=advanced}
二点以上をもつ連結な距離空間は非可算であることを示せ。

:::solution
$a\ne b$ を取り、$f(x)=d(x,a)$ とする。
三角不等式から $|f(x)-f(y)|\le d(x,y)$ なので $f$ は連続である。
$f(a)=0$、$f(b)=d(a,b)>0$ であり、$f(X)$ は連結なので区間で $[0,d(a,b)]$ を含む。
$f(X)$ は非可算集合 $[0,d(a,b)]$ を含むので、定義域 $X$ も非可算である。
二点以上をもつ可算な距離空間（例えば $\Q$）が連結でないことはこの特別な場合である。
:::
::::

::::exercise{level=advanced}
$\R^{\N}$ に箱位相（すべての $\prod_nU_n$、$U_n$ は開、を開基とする位相）を入れると連結でないことを示せ。

:::solution
有界な数列全体を $B$ とする。
$x\in B$ なら $\prod_n(x_n-1,x_n+1)$ は箱位相の開集合で、その元はすべて有界なので $B$ は開である。
$x\notin B$ なら $\prod_n(x_n-1,x_n+1)$ の元はすべて非有界なので、補集合も開である。
$0\in B$、$(n)_n\notin B$ なので $(B,\R^{\N}\setminus B)$ は分離である。
積位相では :ref[thm-top-connected-product] により $\R^{\N}$ は連結であり、
その証明が基本開集合の座標制限の有限性に依存していたことが分かる。
:::
::::

::::exercise{level=advanced}
空間 $X$ の各点 $x$ と $x$ の各開近傍 $U$ に対し、$x$ の連結な近傍 $N\subset U$（$N$ は開でなくてよい）が存在するとする。
このとき $X$ は局所連結であることを示せ。

:::solution
:ref[thm-top-local-components] (1) により、開集合 $U$ の成分 $C$ が開であることを示せばよい。
$x\in C$ を取ると、仮定から $x$ の連結な近傍 $N\subset U$ がある。
$N$ は $U$ の連結部分集合で $x$ を含むので、成分の最大性から $N\subset C$ である。
$N$ は $x$ の近傍なので $x\in N^\circ\subset C$ となり、$x$ は $C$ の内点である。
$x$ は任意なので $C$ は開である。
:::
::::

## 確認テスト

::::quiz
平面の櫛形空間 $H=([0,1]\times\{0\})\cup(D\times[0,1])$（$D=\{0\}\cup\{1/n:n\ge1\}$）について正しいものを選べ。

- [x] 弧状連結である。
- [ ] 各点のどの開近傍の中にも、その点を含む連結開近傍がある。
- [ ] 弧状連結成分が無限個ある。

:::explanation
正解は (A) である。
どの点からも自分の歯を底辺まで降り、底辺に沿って原点に行けるので、$H$ は弧状連結であり、弧状連結成分は一つである（(C) は誤り）。
(B) は誤りである。
点 $(0,1/2)$ の近傍 $H\cap(\R\times(1/4,3/4))$ に含まれる連結集合の第一座標は $D$ の連結部分集合なので一点であり、
一方 $(0,1/2)$ のどの開近傍も $(1/n,1/2)$ を含むので、連結な開近傍は取れない。
:::
::::

::::quiz{answer="3" tol="1e-6"}
実数直線から相異なる二点を除いた空間の連結成分はいくつか。

:::explanation
二点を $a<b$ とすると、残りは $(-\infty,a)$、$(a,b)$、$(b,\infty)$ の三つの開区間である。
各区間は連結で、互いに相手の補集合の一部として開かつ閉なので、除いた点をまたぐ連結集合はない。
答えは $3$ である。
:::
::::

::::quiz
連結集合 $A\subset X$ について、常に連結であるものをすべて選べ。

- [x] その閉包 $\overline A$
- [x] 連続写像 $f:X\to Y$ による像 $f(A)$
- [ ] 任意の別の連結集合 $B$ との共通部分 $A\cap B$

:::explanation
正解は (A) と (B) である。
閉包は :ref[prop-top-connected-union] (3)、連続像は :ref[prop-top-connected-image] による。
(C) は誤りで、円周の上半分 $A=\{(x,y)\in S^1:y\ge0\}$ と下半分 $B=\{(x,y)\in S^1:y\le0\}$ はそれぞれ連結だが、
$A\cap B=\{(1,0),(-1,0)\}$ は二点の離散空間で連結でない。
:::
::::

::::quiz{answer="2" tol="1e-6"}
位相幾何学者の正弦曲線 $T=\{(x,\sin(1/x)):0<x\le1\}\cup(\{0\}\times[-1,1])$ の弧状連結成分はいくつか。

:::explanation
グラフ部分 $G$ は $(0,1]$ の連続像、縦線分は線分なので、それぞれ弧状連結である。
縦線分から $G$ への道 $(u,v)$ があると、最後に $u=0$ となる時刻の直後に $u$ が小さな正の値を取り続け、
$v=\sin(1/u)$ が $1$ と $-1$ を無限回取るので、$v$ がその時刻で連続になれない。
従って弧状連結成分は二つである（:ref[thm-top-sine-curve]）。
:::
::::

::::quiz
局所弧状連結な空間について常に成り立つものを選べ。

- [ ] 空間は連結である。
- [x] 連結成分と弧状連結成分が一致する。
- [ ] 連結成分はすべて一点である。

:::explanation
正解は (B) である。
:ref[thm-top-local-components] により、局所弧状連結空間の弧状連結成分は開であり、補集合が他の弧状連結成分の和として開なので閉でもある。
従って連結成分は弧状連結成分からはみ出せず、両者は一致する。
二つの離れた開区間の和は局所弧状連結だが連結でなく、成分は一点でもないので、(A) と (C) は誤りである。
:::
::::

::::quiz{answer="2" tol="1e-6"}
$\GL_3(\R)$（$M_3(\R)\cong\R^9$ の部分空間）の連結成分の個数を答えよ。

:::explanation
行列式は連続なので、$\GL_3^{+}(\R)=\{\det>0\}$ と $\GL_3^{-}(\R)=\{\det<0\}$ は交わらない空でない開集合であり、$\GL_3(\R)$ は連結でない。
:ref[prop-top-gl-components] により、行基本変形 $R_{ij}(c)$ の係数を $0$ へ縮める道で、
各元は $I$ または $\diag(1,1,-1)$ と行列式の符号を保ったまま結べるので、二つの部分はそれぞれ弧状連結である。
従って連結成分はちょうど $2$ 個である。
:::
::::

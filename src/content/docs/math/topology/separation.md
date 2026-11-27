---
title: 分離公理
description: 点や閉集合を開集合と連続関数で分離する条件を比較し、Urysohn の補題、Tietze の拡張定理、Tychonoff の埋め込み定理、Urysohn の距離化定理を証明する。
notation:
  - symbol: 'T_0,\ T_1,\ T_2,\ T_3,\ T_{3\frac12},\ T_4'
    meaning: 分離公理（T_3, T_{3½}, T_4 は正則・完全正則・正規に T_1 を加えたもの）
  - symbol: '\Delta_X'
    meaning: 積空間 X×X の対角集合 {(x,x)}
  - symbol: 'D'
    meaning: 区間 [0,1] 内の二進有理数 k/2^n 全体
  - symbol: 'Z(g)'
    meaning: 実数値関数 g の零点集合 g^{-1}(0)
  - symbol: '\operatorname{supp}\phi'
    meaning: 関数 φ の台（非零集合の閉包）
  - symbol: 'C(X,[0,1])'
    meaning: X から [0,1] への連続関数全体
  - symbol: '[0,1]^J'
    meaning: 単位区間の J 個の積（立方体）に積位相を入れた空間
---

## この章の目標

- 分離公理 $T_0$〜$T_4$ と完全正則性の定義を固定し、含意関係と、逆が成り立たないことを示す反例を整理する。
- Hausdorff 性を対角集合とネットの極限の一意性で特徴付け、部分空間・積・商での振る舞いを調べる。
- 正規空間で Urysohn の補題を証明し、閉集合上の連続関数を全体へ延長する Tietze の拡張定理を導く。
- 完全正則空間を立方体へ埋め込み（Tychonoff の埋め込み定理）、第二可算な正則空間が距離化可能であること（Urysohn の距離化定理）を証明する。

[距離空間](/math/topology/metric-spaces/) では、異なる二点 $x,y$ は半径 $d(x,y)/2$ の開球で囲めば互いに交わらない近傍に入る。
さらに閉集合 $A$ までの距離 $d(x,A)$ を使えば、点や閉集合を連続関数の値で区別できる。
一般の位相空間ではこのような「分ける」操作が保証されない。
[位相空間](/math/topology/topological-spaces/) の章で見た密着位相や補有限位相では、開集合が少なすぎて点を区別できないことがある。
本章では「何を何から分離できるか」を段階的な公理として整理する。
後半では、閉集合どうしを開集合で分離できるという条件（正規性）から、
距離関数の代わりとなる連続関数を構成する。
その結果、どのような位相空間が距離から来るかという問題（距離化問題）に一つの答えが得られる。

## T₀〜T₄ 空間

### 分離公理の定義

分離公理は「区別する」「囲む」「関数で分ける」という三種類の操作を、
二点、点と閉集合、閉集合どうしという三種類の対象に適用したものである。
対象が大きくなるほど条件は強くなり、得られる道具も強力になる。

:::definition[分離公理]{#def-top-separation}
$X$ を位相空間とする。
二つの部分集合 $A,B\subset X$ を :term[分離する]{en="separate" yomi="ぶんりする"} とは、
$A\subset U$、$B\subset V$、$U\cap V=\varnothing$ を満たす開集合 $U,V$ が存在することをいう。
点 $x$ については $A=\{x\}$ として同じ言葉を使う。

- $T_0$（Kolmogorov）：異なる二点 $x,y$ に対し、一方を含み他方を含まない開集合が存在する。
- $T_1$：異なる二点 $x,y$ に対し、$x$ を含み $y$ を含まない開集合と、$y$ を含み $x$ を含まない開集合がともに存在する。
- $T_2$：異なる二点を分離できる。$T_2$ 空間を :term[Hausdorff 空間]{en="Hausdorff space" yomi="はうすどるふくうかん"} という。
- $X$ が :term[正則]{en="regular" yomi="せいそく"} であるとは、閉集合 $F$ と $F$ に属さない点 $x$ を常に分離できることである。正則な $T_1$ 空間を $T_3$ 空間という。
- $X$ が :term[正規]{en="normal" yomi="せいき"} であるとは、交わらない二つの閉集合を常に分離できることである。正規な $T_1$ 空間を $T_4$ 空間という。
:::

:::remark[流儀の違いについて]
正則・正規に $T_1$ を含めるかどうか、また $T_3$・$T_4$ という記号がどちらを指すかは文献によって異なる。
例えば :cite[munkres-topology] では正則・正規の定義に「一点集合が閉である」ことを含めている。
本章では上の定義に固定し、$T_1$ を含む条件には必ず $T_3$・$T_4$ の記号を用いる。
文献を読むときは、まず定義がどちらの流儀かを確認する必要がある。
:::

$T_1$ の条件は、一点集合の閉性という扱いやすい形に言い換えられる。
これにより、点を「閉集合の特別な場合」として扱え、公理の間の含意が見通しよくなる。

:::proposition[公理の含意と一点集合]{#prop-top-separation-chain}
$X$ が $T_1$ であることは、すべての一点集合が閉であることと同値である。
またこの約束の下で

$$
T_4\Longrightarrow T_3\Longrightarrow T_2\Longrightarrow T_1\Longrightarrow T_0
$$

が成り立つ。
:::

:::proof
$X$ を $T_1$ とし、$x\in X$ を固定する。
各 $y\ne x$ に対し、$y$ を含み $x$ を含まない開集合 $V_y$ がある。
すると $X\setminus\{x\}=\bigcup_{y\ne x}V_y$ である。
実際、右辺の各 $V_y$ は $x$ を含まないので右辺は左辺に含まれ、左辺の各点 $y$ は $V_y$ に属する。
開集合の和は開なので $\{x\}$ は閉である。
逆にすべての一点集合が閉なら、$x\ne y$ に対し $X\setminus\{y\}$ は $x$ を含み $y$ を含まない開集合であり、
$X\setminus\{x\}$ はその逆の役割をもつ開集合である。
よって $X$ は $T_1$ である。

$T_4\Rightarrow T_3$：$X$ を $T_4$ とし、閉集合 $F$ と $x\notin F$ を取る。
$T_1$ 性から $\{x\}$ は閉であり、$\{x\}$ と $F$ は交わらない二つの閉集合である。
正規性によりこれらを分離できるので $X$ は正則であり、$T_1$ と合わせて $T_3$ である。

$T_3\Rightarrow T_2$：$x\ne y$ とすると $\{y\}$ は閉で $x\notin\{y\}$ である。
正則性により $x$ と $\{y\}$、すなわち $x$ と $y$ を分離できる。

$T_2\Rightarrow T_1$：$x\in U$、$y\in V$、$U\cap V=\varnothing$ となる開集合を取ると、
$U$ は $x$ を含み $y$ を含まず、$V$ は $y$ を含み $x$ を含まない。

$T_1\Rightarrow T_0$ は定義から明らかである。
:::

### 公理が成り立たない典型例

各公理の強さを実感するため、それが成り立たない最小限の例を調べる。
次の三つの空間は後で反例として何度も使う。

:::example[分離の失敗の典型例]{#exm-top-separation-failures}
**二点密着空間。**
$X=\{a,b\}$ に開集合 $\varnothing,X$ だけを入れる。
二点を含む開集合は $X$ しかないので、どちらか一方だけを含む開集合はなく、$T_0$ でない。
一方、閉集合も $\varnothing,X$ だけである。
点 $x$ を含まない閉集合は $\varnothing$ だけなので、$x$ と $\varnothing$ は開集合 $X$ と $\varnothing$ で分離できる。
交わらない二閉集合の組も一方が空集合であり、同様に分離できる。
従ってこの空間は正則かつ正規であるが $T_0$ ですらない。

**Sierpiński 空間。**
$X=\{0,1\}$ に開集合 $\varnothing,\{1\},X$ を入れる（[位相空間](/math/topology/topological-spaces/) の例）。
開集合 $\{1\}$ が $1$ を含み $0$ を含まないので $T_0$ である。
しかし $0$ を含み $1$ を含まない開集合はないから $T_1$ でない。
閉集合は $\varnothing,\{0\},X$ であり、交わらない二閉集合の組は一方が空なので正規である。
ところが点 $1$ と閉集合 $\{0\}$ を分離しようとすると、$\{0\}$ を含む開集合は $X$ だけで、必ず $1$ を含む。
従って正規であるが正則でない。

**無限集合の補有限位相。**
$X$ を無限集合とし、空集合と補集合が有限な集合を開とする。
一点集合の補集合は開なので一点集合は閉であり、$T_1$ である。
しかし空でない開集合 $U,V$ の補集合はともに有限だから、
$X\setminus(U\cap V)=(X\setminus U)\cup(X\setminus V)$ も有限であり、$X$ が無限なので $U\cap V\ne\varnothing$ である。
従って異なる二点を分離できず $T_2$ でない。
同じ理由で、点 $x$ と $x$ を含まない空でない有限集合（閉集合）も分離できず、正則でも正規でもない。
:::

最初の二例が示すように、$T_1$ を仮定しない正則性・正規性は非常に弱い条件であり、
Hausdorff 性すら導かない。
「正規 $\Rightarrow$ 正則 $\Rightarrow$ Hausdorff」という連鎖は $T_1$ を付け加えて初めて成立する。
これが $T_3,T_4$ の定義に $T_1$ を含める理由である。

### 一点の閉包による判定

$T_0$ と $T_1$ は、各点の閉包を調べるだけで判定できる。
$y\in\overline{\{x\}}$ は「$y$ のどの開近傍も $x$ を含む」ことと同値であり、
これは $y$ の近くに常に $x$ がいる、つまり開集合では $y$ を $x$ から引き離せないことを意味する。

:::proposition[一点の閉包と識別可能性]{#prop-top-tzero-closure}
$X$ が $T_0$ であることは、異なる点 $x,y$ に対して常に $\overline{\{x\}}\ne\overline{\{y\}}$ となることと同値である。
特に有限空間が $T_1$ であるための必要十分条件は、離散空間であることである。
:::

:::proof
まず一般に、$y\in\overline{\{x\}}$ と「$y$ を含むすべての開集合が $x$ を含む」ことは同値である。
実際、閉包の点とは、そのどの開近傍も集合と交わる点のことだからである。

$x\ne y$ で $\overline{\{x\}}=\overline{\{y\}}$ とすると、$y\in\overline{\{y\}}=\overline{\{x\}}$ より $y$ を含む開集合はすべて $x$ を含み、
対称に $x$ を含む開集合はすべて $y$ を含む。
よって一方だけを含む開集合はなく、$X$ は $T_0$ でない。

逆に $X$ が $T_0$ でないとし、どちらか一方だけを含む開集合をもたない $x\ne y$ を取る。
上の言い換えにより $y\in\overline{\{x\}}$ かつ $x\in\overline{\{y\}}$ である。
$\overline{\{x\}}$ は $y$ を含む閉集合なので $\{y\}$ の閉包を含み、$\overline{\{y\}}\subset\overline{\{x\}}$ である。
対称に逆の包含も成り立つので $\overline{\{x\}}=\overline{\{y\}}$ となる。

最後に $X$ を有限な $T_1$ 空間とする。
:ref[prop-top-separation-chain] により一点集合は閉であり、任意の部分集合 $A$ は有限個の一点集合の和として閉である。
よって $X\setminus A$ も閉、すなわち $A$ は開であり、$X$ は離散空間である。
逆に離散空間ではすべての部分集合が閉なので $T_1$ である。
:::

例えば三点の鎖で上方閉集合を開とすると、各点の閉包はそれ以下の点の集合である。閉包は三つとも異なるが $\{2\}$ は閉でないので、この空間は $T_0$ であって $T_1$ ではない。

開近傍で区別できない点をまとめると :ref[prop-top-kolmogorov-quotient] の商を得る。

### 近傍を縮める

正則性と正規性を実際に使うときは、「二つの開集合で囲む」形よりも
「開集合の中に、閉包ごと収まる小さな開集合を取る」形が便利である。
距離空間では、点 $x$ の開近傍 $G$ が $B(x,r)$ を含むとき、
$U=B(x,r/2)$ の閉包は閉球 $\{y:d(x,y)\le r/2\}$ に含まれ、従って $G$ に含まれる。
この「一回り小さい近傍」を一般の空間で取れることが正則性にほかならない。

:::lemma[近傍を縮める条件]{#lem-top-shrinking}
$X$ が正則であることは、任意の開集合 $G$ と $x\in G$ に対して
$x\in U\subset\overline U\subset G$ となる開集合 $U$ が存在することと同値である。
$X$ が正規であることは、任意の開集合 $G$ と閉集合 $A\subset G$ に対して
$A\subset U\subset\overline U\subset G$ となる開集合 $U$ が存在することと同値である。
:::

:::proof
正則の場合を示す。
$X$ を正則とし、$x\in G$（$G$ は開）とする。
$F=X\setminus G$ は $x$ を含まない閉集合だから、$x\in U$、$F\subset V$、$U\cap V=\varnothing$ となる開集合 $U,V$ がある。
$U\subset X\setminus V$ で右辺は閉なので $\overline U\subset X\setminus V$ である。
$F\subset V$ より $X\setminus V\subset X\setminus F=G$ だから $\overline U\subset G$ を得る。

逆に縮小条件が成り立つとし、閉集合 $F$ と $x\notin F$ を取る。
$G=X\setminus F$ は $x$ を含む開集合なので、$x\in U\subset\overline U\subset G$ となる開集合 $U$ がある。
$V=X\setminus\overline U$ は開で、$\overline U\subset X\setminus F$ より $F\subset V$ である。
$U\cap V\subset\overline U\cap(X\setminus\overline U)=\varnothing$ だから、$U,V$ が $x$ と $F$ を分離する。

正規の場合は、上の議論で点 $x$ を閉集合 $A$ に置き換えればよい。
:::

### 部分空間への遺伝

:::proposition[部分空間への遺伝]{#prop-top-separation-hereditary}
$T_0$、$T_1$、$T_2$、正則性、$T_3$ は任意の部分空間に受け継がれる。
正規性は閉部分空間に受け継がれる。
:::

:::proof
$A\subset X$ を部分空間とする。
$x,y\in A$ を区別または分離する $X$ の開集合を $A$ と交わすと、:ref[def-top-subspace] により $A$ の開集合となる。
交わしても $x,y$ の所属は変わらず、交わらない二集合は交わらないままなので、$T_0,T_1,T_2$ は $A$ に受け継がれる。

正則性を示す。
$F$ を $A$ の閉集合、$a\in A\setminus F$ とする。
:ref[prop-top-subspace-laws] により $F=A\cap F'$（$F'$ は $X$ の閉集合）と書ける。
$a\notin F'$ である。実際 $a\in F'$ なら $a\in A\cap F'=F$ となる。
$X$ の正則性から $a\in U$、$F'\subset V$、$U\cap V=\varnothing$ となる開集合 $U,V$ があり、
$A\cap U$ と $A\cap V$ が $A$ で $a$ と $F$ を分離する。
$T_1$ も受け継がれるので $T_3$ も受け継がれる。

最後に $A$ を $X$ の閉部分空間、$X$ を正規とする。
$A$ の交わらない閉集合 $E_1,E_2$ は、$A$ の閉集合が $A\cap F$（$F$ は $X$ で閉）の形であることと $A$ が閉であることから、$X$ でも閉である。
$X$ の正規性で $E_1\subset U_1$、$E_2\subset U_2$、$U_1\cap U_2=\varnothing$ とし、$A\cap U_1,A\cap U_2$ を取ればよい。
:::

正規性は一般の部分空間には受け継がれない。
これは本章の最後の節で、完全正則空間の立方体への埋め込みを使って示す（:ref[exm-top-normal-not-hereditary]）。
閉部分空間に限られるのは、上の証明で $A$ の閉集合が $X$ でも閉であることを使ったからである。

### Hausdorff だが正則でない直線

$T_2$ と $T_3$ の間には本質的な差がある。
次の例では、通常の直線の位相に一つの閉集合を「人工的に」付け加えることで、
点と閉集合を分離できなくする。

:::example[K 位相の直線]{#exm-top-k-topology}
$K=\{1/n:n\ge1\}$ とし、$\R$ 上で

$$
\mathcal B=\{(a,b):a<b\}\cup\{(a,b)\setminus K:a<b\}
$$

を開基とする位相を考える。

**開基の条件。**
$\mathcal B$ の和は $\R$ である。
二つの元の交わりは、$(a,b)\cap(c,d)$ が開区間（または空）であり、
$(a,b)\cap\bigl((c,d)\setminus K\bigr)=\bigl((a,b)\cap(c,d)\bigr)\setminus K$、
$\bigl((a,b)\setminus K\bigr)\cap\bigl((c,d)\setminus K\bigr)=\bigl((a,b)\cap(c,d)\bigr)\setminus K$ なので、
再び $\mathcal B$ の元または空集合である。
よって :ref[prop-top-basis-test] により $\mathcal B$ はある位相の開基である。

**Hausdorff 性と $K$ の閉性。**
この位相は開区間をすべて開集合として含むので、通常の位相より細かい。
通常の位相で異なる二点を分離する開区間はこの位相でも開なので、$\R$ は Hausdorff である。
また $\R\setminus K=\bigcup_{n\ge1}\bigl((-n,n)\setminus K\bigr)$ は開だから $K$ は閉である。

**正則でないこと。**
$0\notin K$ である。
$0$ と $K$ を分離する開集合 $U\ni0$、$V\supset K$ があったと仮定する。
$U$ は $0$ を含む開基の元を含むが、それは $(a,b)$ か $(a,b)\setminus K$（$a<0<b$）の形である。
いずれの場合も $(a,b)\setminus K\subset U$ である。
$1/n<b$ となる $n$ を取ると、$1/n\in K\subset V$ なので $V$ は $1/n$ を含む開基の元を含む。
$(c,d)\setminus K$ の形の元は $1/n$ を含まないから、それは開区間 $(c,d)$（$c<1/n<d$）である。
ここで $\max\{c,1/(n+1)\}<z<1/n$ となる $z$ を取ると、$z$ は二つの連続する点 $1/(n+1),1/n$ の間にあるので $z\notin K$ であり、
$0<z<1/n<b$ だから $z\in(a,b)\setminus K\subset U$、また $z\in(c,d)\subset V$ である。
これは $U\cap V=\varnothing$ に反する。
従ってこの空間は $T_2$ であるが正則でない（$T_3$ でない）。
:::

この例で $0$ の近傍は $K$ を避けられるが、$K$ の近傍は各 $1/n$ の両側の小区間を含まざるをえない。
$1/n$ が $0$ に集積するため、その小区間たちが $0$ の近傍に必ず食い込む。
正則性の失敗は、このように閉集合の近傍が点の近傍から逃げられない状況として現れる。

## Hausdorff 空間

### 極限の一意性と対角集合

距離空間では、列の極限はただ一つである。
解析学で極限を「$\lim$」という記号で書けるのはこの一意性のおかげである。
一般の位相空間では、例えば密着空間で定数列がすべての点に収束するように、極限は一意とは限らない。
Hausdorff 性はちょうど極限の一意性を保証する条件である。
ただし点列では足りず、[位相空間](/math/topology/topological-spaces/) の章で導入したネットを使う必要がある。

:::theorem[対角集合と極限]{#thm-top-hausdorff-diagonal}
位相空間 $X$ について、次の三条件は同値である。

1. $X$ は Hausdorff である。
2. 対角集合 $\Delta_X=\{(x,x):x\in X\}$ は $X\times X$ で閉である。
3. $X$ の任意のネットは高々一点にしか収束しない。
:::

:::proof
**1 と 2 の同値。**
$\Delta_X$ が閉であることは、補集合 $X\times X\setminus\Delta_X$ が開であること、
すなわち各点 $(x,y)$（$x\ne y$）が $\Delta_X$ と交わらない開近傍をもつことである。
積位相では各開近傍が基本開集合 $U\times V$（$U\ni x$、$V\ni y$ は開）を含むので、
これは「$U\times V$ が $\Delta_X$ と交わらない $U,V$ が存在する」ことと同値である。
$(U\times V)\cap\Delta_X=\{(z,z):z\in U\cap V\}$ なので、この条件は $U\cap V=\varnothing$ と同値であり、
まさに $x$ と $y$ の分離である。

**1 から 3。**
ネット $(x_\alpha)_{\alpha\in A}$ が異なる二点 $x,y$ に収束したとする。
$x\in U$、$y\in V$、$U\cap V=\varnothing$ となる開集合を取る。
収束の定義から、ある $\alpha_1$ 以後のすべての項が $U$ に、ある $\alpha_2$ 以後のすべての項が $V$ に入る。
有向集合なので $\alpha_1,\alpha_2$ の共通の上界 $\alpha_3$ があり、$x_{\alpha_3}\in U\cap V=\varnothing$ となって矛盾する。

**3 から 1。**
対偶を示す。
$x\ne y$ が分離できないとする。
$x$ の開近傍 $U$ と $y$ の開近傍 $V$ の組 $(U,V)$ 全体を $\Lambda$ とし、
$(U,V)\le(U',V')$ を $U'\subset U$ かつ $V'\subset V$ で定める。
これは反射的・推移的であり、$(U\cap U',V\cap V')$ が $(U,V),(U',V')$ の共通の上界なので、$\Lambda$ は有向集合である。
仮定より各組で $U\cap V\ne\varnothing$ なので、点 $z_{(U,V)}\in U\cap V$ を選ぶ（選択公理を使う）。
$x$ の任意の開近傍 $H$ に対し、$(U,V)\ge(H,X)$ なら $z_{(U,V)}\in U\subset H$ である。
よってこのネットは $x$ に収束し、同様に $y$ にも収束する。
極限が一意でないので 3 は成り立たない。
:::

点列で同じ判定ができるかどうかは空間による。
第一可算空間では点列で十分である。
実際、$x\ne y$ が分離できず、$x,y$ が減少する可算基本近傍系 $U_1\supset U_2\supset\cdots$、$V_1\supset V_2\supset\cdots$ をもつなら、
$z_n\in U_n\cap V_n$ を選んだ点列は $x$ にも $y$ にも収束する。
しかし一般には点列の極限の一意性から Hausdorff 性は従わない。

:::example[補可算位相：点列の極限は一意だが Hausdorff でない]{#exm-top-cocountable}
$X$ を非可算集合とし、空集合と補集合が高々可算な集合を開とする補可算位相を入れる。
空でない二つの開集合 $U,V$ について $X\setminus(U\cap V)$ は高々可算なので、$X$ が非可算であることから $U\cap V\ne\varnothing$ である。
よって $X$ は Hausdorff でない（一点集合の補集合は開なので $T_1$ ではある）。

一方、点列 $x_n\to x$ を取る。
$C=\{x_n:x_n\ne x\}$ は高々可算なので $X\setminus C$ は $x$ を含む開集合であり、
ある番号以後のすべての項がこれに入る。
その項は $C$ に属さないので $x$ に等しい。
従って収束する点列は最終的に定数であり、極限はその定数値ただ一つである。
:::

### 稠密集合で決まる連続写像

対角集合の閉性は、「二つの写像が一致する点の集合は閉である」という形でよく使われる。

:::corollary[稠密集合で決まる連続写像]{#cor-top-dense-uniqueness}
$Y$ を Hausdorff 空間とする。
二つの連続写像 $f,g:X\to Y$ が $X$ の稠密部分集合 $D$ 上で一致すれば、$X$ 全体で $f=g$ である。
:::

:::proof
$h:X\to Y\times Y$、$h(x)=(f(x),g(x))$ は :ref[thm-top-product-property] により連続である。
一致点の集合 $E=\{x:f(x)=g(x)\}$ は $h^{-1}(\Delta_Y)$ に等しく、
:ref[thm-top-hausdorff-diagonal] により $\Delta_Y$ が閉なので $E$ は閉である。
$D\subset E$ より $X=\overline D\subset E$ となる。
:::

例えば $\R$ 上の連続関数 $f$ が有理数上で $f(q)=q^2$ を満たすなら、$\Q$ の稠密性から $f(x)=x^2$ がすべての実数で成り立つ。
値域の Hausdorff 性は省けない。

:::example[二原点の直線での一意性の破れ]
:ref[exm-top-double-origin] の二原点をもつ直線 $X$ を考える。
二本の直線 $\R_a,\R_b$ から $X$ への写像として、それぞれの包含と商写像の合成 $\iota_a,\iota_b:\R\to X$ を取る。
両者は連続であり、稠密な $\R\setminus\{0\}$ 上で一致するが、$\iota_a(0)=0_a\ne0_b=\iota_b(0)$ である。
:ref[cor-top-dense-uniqueness] が成り立たないので、:ref[thm-top-hausdorff-diagonal] により $X$ の対角集合は閉でなく、
点列 $1/n$ の像は $0_a$ にも $0_b$ にも収束する。
:::

### 部分空間・積・商

:::proposition[遺伝性と積]{#prop-top-hausdorff-product}
Hausdorff 空間の部分空間と、Hausdorff 空間の任意個の積は Hausdorff である。
コンパクト Hausdorff 空間は正規であり、従って $T_4$ 空間である。
:::

:::proof
部分空間については :ref[prop-top-separation-hereditary] で示した。

$P=\prod_{i\in I}X_i$ で $x\ne y$ とすると、ある座標 $i$ で $x_i\ne y_i$ である。
$X_i$ の Hausdorff 性で $x_i\in U_i$、$y_i\in V_i$、$U_i\cap V_i=\varnothing$ となる開集合を取る。
$\pi_i^{-1}(U_i)$ と $\pi_i^{-1}(V_i)$ は積位相で開であり、それぞれ $x,y$ を含む。
$z$ が両方に属すれば $z_i\in U_i\cap V_i$ となるので、両者は交わらない。

最後に $X$ をコンパクト Hausdorff 空間とし、$A,B$ を交わらない閉集合とする。
:ref[prop-top-compact-closed] により $A,B$ はコンパクトである。
:ref[lem-top-compact-separation] により、Hausdorff 空間の交わらない二つのコンパクト集合は交わらない開集合で分離できる。
よって $X$ は正規であり、$T_2\Rightarrow T_1$ と合わせて $T_4$ である。
:::

Hausdorff 性は商空間には受け継がれない。
二原点の直線は Hausdorff な $\R\amalg\R$ の商であるが Hausdorff でない。
また Hausdorff 空間では有限集合は閉であり、:ref[prop-top-compact-closed] によりコンパクト集合も閉である。
これを使うと、:ref[thm-top-compact-image] のとおり、コンパクト空間から Hausdorff 空間への連続写像は閉写像になる。

:::example[連続写像のグラフ]
$Y$ が Hausdorff なら、連続写像 $f:X\to Y$ のグラフ $\Gamma_f=\{(x,f(x))\}$ は $X\times Y$ で閉である。
実際、連続写像 $X\times Y\to Y\times Y$、$(x,y)\mapsto(f(x),y)$ による $\Delta_Y$ の逆像が $\Gamma_f$ である。
逆は一般に成り立たない。
$f(x)=1/x$（$x\ne0$）、$f(0)=0$ で定まる $f:\R\to\R$ のグラフは、閉集合である双曲線 $\{xy=1\}$ と一点 $(0,0)$ の和で閉であるが、$f$ は原点で連続でない。
値域がコンパクト Hausdorff の場合には逆が成り立つ（演習問題）。
:::

### 同値関係が閉になる商

商空間が Hausdorff になるかどうかは、構成した空間を幾何学的に使えるかを左右する。
コンパクト Hausdorff 空間の商では、同値関係のグラフの閉性だけで判定できる。

:::proposition[コンパクト Hausdorff 空間の Hausdorff 商]{#prop-top-closed-relation}
$K$ をコンパクト Hausdorff 空間、$q:K\to Q$ を同値関係 $\sim$ による商写像とする。
$Q$ が Hausdorff であるための必要十分条件は

$$
R=\{(s,t)\in K\times K:q(s)=q(t)\}
$$

が $K\times K$ で閉であることである。
このとき $q$ は閉写像である。
:::

:::proof
$Q$ が Hausdorff なら、$R$ は連続写像 $q\times q:K\times K\to Q\times Q$ による閉集合 $\Delta_Q$ の逆像なので閉である。

逆に $R$ が閉であるとする。

**$q$ は閉写像である。**
$C\subset K$ を閉集合とする。
:ref[thm-top-finite-product-compact] により $K\times K$ はコンパクトであり、
$R\cap(K\times C)$ はその閉部分集合なので :ref[prop-top-compact-closed] によりコンパクトである。
第一射影 $p_1:K\times K\to K$ は連続なので、:ref[thm-top-compact-image] により $p_1\bigl(R\cap(K\times C)\bigr)$ はコンパクトである。
この集合は $\{s:\text{ある }t\in C\text{ で }s\sim t\}=q^{-1}(q(C))$ に等しい。
$K$ は Hausdorff なのでコンパクト集合 $q^{-1}(q(C))$ は閉であり、商位相の定義により $q(C)$ は閉である。

**$Q$ は Hausdorff である。**
$Q$ の異なる二点 $a,b$ を取り、$A=q^{-1}(a)$、$B=q^{-1}(b)$ とする。
$K$ は $T_1$ なので一点 $\{s\}$ は閉であり、上で示したことから $A=q^{-1}(q(\{s\}))$（$s\in A$）は閉である。
同様に $B$ も閉で、$A\cap B=\varnothing$ である。
:ref[prop-top-hausdorff-product] により $K$ は正規なので、$A\subset U$、$B\subset V$、$U\cap V=\varnothing$ となる開集合がある。
$W_a=Q\setminus q(K\setminus U)$ と置くと、$q$ が閉写像なので $W_a$ は開である。
$a\in W_a$ である。
実際 $a=q(k)$（$k\notin U$）とすると $k\in q^{-1}(a)=A\subset U$ となって矛盾する。
また $q^{-1}(W_a)=K\setminus q^{-1}(q(K\setminus U))\subset K\setminus(K\setminus U)=U$ である。
同様に $W_b=Q\setminus q(K\setminus V)$ は $b$ を含む開集合で $q^{-1}(W_b)\subset V$ を満たす。
$q^{-1}(W_a\cap W_b)\subset U\cap V=\varnothing$ で $q$ は全射なので $W_a\cap W_b=\varnothing$ である。
:::

:::example[閉集合を一点に潰す商と射影空間]
コンパクト Hausdorff 空間 $K$ の空でない閉集合 $C$ を一点に潰す同値関係のグラフは
$R=\Delta_K\cup(C\times C)$ であり、二つの閉集合の和として閉である。
従って商 $K/C$ は Hausdorff である。
例えば :ref[prop-top-sphere-quotient] の $D^n/S^{n-1}$ がそうである。

同様に、単位球面 $S^n$ で $u\sim-u$ と同一視した射影空間 $\mathbb{RP}^n$ の関係は

$$
R=\{(u,u):u\in S^n\}\cup\{(u,-u):u\in S^n\}
$$

であり、恒等写像と $u\mapsto-u$ という二つの連続写像のグラフの和として閉である。
$S^n$ はコンパクト Hausdorff なので $\mathbb{RP}^n$ は Hausdorff である。

反対に $[0,1]$ の部分集合 $(0,1)$ を一点 $p$ に潰し、端点は別々に残す商を考える。
$q^{-1}(p)=(0,1)$ は閉でないので $\{p\}$ は閉でなく、この商は $T_1$ ですらない。
対応する $R=\Delta\cup(0,1)^2$ の閉包には $(1/n,1/2)\to(0,1/2)$（$n\ge2$）の極限が入る。この極限は対角線にも $(0,1)^2$ にも属さないので、$R$ は閉でない。
:::

### 正則性と積

Hausdorff 性と同様に、正則性も任意個の積に受け継がれる。
証明の要点は、積の基本近傍が有限個の座標しか制限しないことである。

:::proposition[正則性と積位相]{#prop-top-regular-product}
正則空間の任意個の積は正則である。
因子がすべて $T_3$ なら積も $T_3$ である。
:::

:::proof
$P=\prod_{i\in I}X_i$ とし、点 $p\in P$ とその開近傍 $G$ を取る。
:ref[lem-top-shrinking] により、$p\in V\subset\overline V\subset G$ となる開集合 $V$ を見つければよい。
積位相の定義により、有限集合 $J\subset I$ と各 $i\in J$ の開集合 $G_i\ni p_i$ を選んで

$$
p\in\bigcap_{i\in J}\pi_i^{-1}(G_i)\subset G
$$

とできる（$J=\varnothing$ のときは左辺の共通部分を $P$ と読む）。
$i\in J$ ごとに $X_i$ の正則性と :ref[lem-top-shrinking] により

$$
p_i\in V_i\subset\overline{V_i}\subset G_i
$$

となる開集合 $V_i$ を選ぶ。
$V=\bigcap_{i\in J}\pi_i^{-1}(V_i)$ は $p$ の開近傍である。
一方 $\bigcap_{i\in J}\pi_i^{-1}(\overline{V_i})$ は閉集合の逆像の有限個の共通部分なので閉であり、$V$ を含む。
従って

$$
\overline V\subset\bigcap_{i\in J}\pi_i^{-1}(\overline{V_i})\subset\bigcap_{i\in J}\pi_i^{-1}(G_i)\subset G
$$

となり、$P$ は正則である。
全因子が $T_1$ なら、各点について $\{p\}=\bigcap_{i\in I}\pi_i^{-1}(\{p_i\})$ は閉集合の共通部分として閉であり、
:ref[prop-top-separation-chain] により $P$ も $T_1$ である。
:::

この証明では、点の近傍を有限個の座標に関して縮めれば十分だった。
二つの閉集合全体を同時に分けようとすると、その集合がどのように座標にまたがるかを有限個の座標で制御できない。
実際、正規空間の積は正規とは限らない。
次節の Sorgenfrey 平面がその例である。

## 正規空間と Urysohn の補題

距離空間の交わらない閉集合 $A,B$ は、関数

$$
h(x)=\frac{d(x,A)}{d(x,A)+d(x,B)}
$$

で分けられる。
$h$ は $A$ 上で $0$、$B$ 上で $1$ を取る連続関数であり、
$\{h<1/2\}$ と $\{h>1/2\}$ が $A,B$ を分離する開集合を与える。
一般の位相空間には距離がないので、このような関数を書き下せない。
本節の目標は、正規性という開集合についての条件だけから、同じ性質の関数を作り出すことである。

### 正規空間の例

:::proposition[距離空間の正規性]{#prop-top-metric-normal}
距離空間は $T_4$ 空間である。
さらに任意の閉集合 $A$ はある連続関数 $g:X\to[0,1]$ の零点集合 $Z(g)=g^{-1}(0)$ である。
:::

:::proof
一点集合 $\{x\}$ の補集合の各点 $y$ は開球 $B(y,d(x,y))$ をもち、これは $x$ を含まないので、$\{x\}$ は閉であり $X$ は $T_1$ である。

$A,B$ を交わらない閉集合とする。
一方が空なら、空集合と $X$ が分離を与える。
両方が空でないとき、:ref[prop-top-distance-set] により $x\mapsto d(x,A)$ は $1$-Lipschitz で連続であり、
$d(x,A)=0$ は $x\in\overline A=A$ と同値である。
$x\in A$ なら $x\notin B$ なので $d(x,B)>0$ であり、分母 $d(x,A)+d(x,B)$ はどこでも正である。
よって上の $h$ は連続で、$A$ 上で $0$、$B$ 上で $1$ を取り、値は $[0,1]$ に入る。
$U=\{h<1/2\}$、$V=\{h>1/2\}$ は開区間の逆像として開で、交わらず、それぞれ $A,B$ を含む。

最後の主張は、$A\ne\varnothing$ なら $g(x)=\min\{1,d(x,A)\}$、$A=\varnothing$ なら $g\equiv1$ とすればよい。
:::

例えば $\R$ で $A=(-\infty,0]$、$B=[2,\infty)$ なら、$0\le x\le2$ で $d(x,A)=x$、$d(x,B)=2-x$ なので $h(x)=x/2$ である。
$x\le0$ では $h=0$、$x\ge2$ では $h=1$ となり、$h$ は折れ線である。

コンパクト Hausdorff 空間が $T_4$ であることは :ref[prop-top-hausdorff-product] で示した。
もう一つの重要な十分条件は可算性から来る。

:::lemma[正則 Lindelöf 空間の正規性]{#lem-top-regular-lindelof}
正則で Lindelöf な空間は正規である。
:::

:::proof
**閉部分集合も Lindelöf である。**
$A$ を閉集合とし、$A$ を覆う $X$ の開集合族 $(O_\lambda)$ を取る。
これに開集合 $X\setminus A$ を加えると $X$ の開被覆となり、可算部分被覆が取れる。
そこから $X\setminus A$ を除いたものは $A$ を覆う可算部分族である。

**可算個の縮小近傍を取る。**
$A,B$ を交わらない閉集合とする。
一方が空なら空集合と $X$ で分離できるので、両方が空でないとする。
各 $a\in A$ に対し、$a\in X\setminus B$ と :ref[lem-top-shrinking] から、$a\in U_a\subset\overline{U_a}\subset X\setminus B$ となる開集合 $U_a$ がある。
$(U_a)_{a\in A}$ は $A$ を覆うので、可算個 $U_1,U_2,\ldots$ で $A$ を覆える（有限個で済むときは最後の集合を繰り返して列にする）。
各 $n$ で $\overline{U_n}\cap B=\varnothing$ である。
同様に $B$ を覆う開集合列 $V_1,V_2,\ldots$ で $\overline{V_n}\cap A=\varnothing$ となるものを取る。

**重なりを取り除く。**
次の二つの集合を考える。

$$
U=\bigcup_{n\ge1}\Bigl(U_n\setminus\bigcup_{i\le n}\overline{V_i}\Bigr),
\qquad
V=\bigcup_{n\ge1}\Bigl(V_n\setminus\bigcup_{i\le n}\overline{U_i}\Bigr)
$$

各項は開集合から有限個の閉集合の和を除いたものなので開であり、$U,V$ は開である。
$a\in A$ はある $U_n$ に属し、どの $\overline{V_i}$ にも属さないので $a\in U$ である。
よって $A\subset U$、同様に $B\subset V$ である。
$x\in U\cap V$ とし、$x$ が $U$ の第 $n$ 項と $V$ の第 $m$ 項に属するとする。
$m\le n$ なら、$x\in U_n\setminus\bigcup_{i\le n}\overline{V_i}$ より $x\notin\overline{V_m}\supset V_m$ であり、$x\in V_m$ に反する。
$n\le m$ なら、$x\in V_m\setminus\bigcup_{i\le m}\overline{U_i}$ より $x\notin U_n$ であり、やはり矛盾する。
従って $U\cap V=\varnothing$ である。
:::

:::corollary[第二可算な正則空間の正規性]{#cor-top-second-countable-normal}
第二可算な正則空間は正規である。
特に第二可算な $T_3$ 空間は $T_4$ である。
:::

:::proof
:ref[prop-top-countable-basis] により第二可算空間は Lindelöf であるから、:ref[lem-top-regular-lindelof] を適用すればよい。
:::

:::example[Sorgenfrey 直線の正規性]{#exm-top-sorgenfrey-normal}
Sorgenfrey 直線 $S$（半開区間 $[a,b)$ を開基とする $\R$）を考える。
$[a,b)$ は閉でもある。
実際、補集合 $(-\infty,a)\cup[b,\infty)=\bigcup_{n\ge1}[a-n,a)\cup\bigcup_{n\ge1}[b,b+n)$ は開だからである。
$x\in G$（$G$ は開）なら $x\in[x,b)\subset G$ となる $b$ があり、$[x,b)$ は閉なので $\overline{[x,b)}=[x,b)\subset G$ である。
よって :ref[lem-top-shrinking] により $S$ は正則である。
$S$ は通常の位相より細かいので $T_1$ である。
:ref[prop-top-sorgenfrey-lindelof] により $S$ は Lindelöf なので、:ref[lem-top-regular-lindelof] から $S$ は $T_4$ である。
ただし $S$ は第二可算でないので、:ref[cor-top-second-countable-normal] は使えない。
:::

### 正規でない正則空間

正規性を否定するには、交わらない閉集合の対で分離できないものを見つける必要がある。
直接見つけるのは難しいが、次の数え上げの議論が有効である。

:::lemma[Jones の補題]{#lem-top-jones}
$X$ を可分な正規空間とし、$\Lambda\subset X$ を閉集合で、部分空間として離散であるとする。
このとき $\Lambda$ の部分集合全体の濃度は、可算集合の部分集合全体の濃度 $2^{\aleph_0}$ 以下である。
:::

:::proof
$Q\subset X$ を可算稠密集合とする。
$\Lambda$ は離散なので任意の $E\subset\Lambda$ は $\Lambda$ で閉であり、$\Lambda$ が閉なので $X$ でも閉である。
$\Lambda\setminus E$ も同様に閉なので、正規性から $E\subset U_E$、$\Lambda\setminus E\subset V_E$、$U_E\cap V_E=\varnothing$ となる開集合を選べる。
写像 $E\mapsto U_E\cap Q$ が単射であることを示す。
$E\ne E'$ とし、例えば $x\in E\setminus E'$ とする（逆の場合は役割を入れ替える）。
$x\in U_E$ かつ $x\in\Lambda\setminus E'\subset V_{E'}$ なので、$U_E\cap V_{E'}$ は空でない開集合であり、稠密集合の点 $q\in Q$ を含む。
$q\in U_E$ であり、$q\in V_{E'}$ より $q\notin U_{E'}$ である。
従って $U_E\cap Q\ne U_{E'}\cap Q$ となる。
よって $\Lambda$ の部分集合全体から $Q$ の部分集合全体への単射が得られた。
:::

:::example[Sorgenfrey 平面は正規でない]{#exm-top-sorgenfrey-plane}
Sorgenfrey 直線 $S$ の積 $S^2$ を考える。
開基は長方形 $[a,b)\times[c,d)$ である。
:ref[exm-top-sorgenfrey-normal] と :ref[prop-top-regular-product] により $S^2$ は $T_3$ である。
[位相空間](/math/topology/topological-spaces/) の章で見たとおり、$\Q^2$ は $S^2$ で稠密であり、$S^2$ は可分である。

反対角線 $\Lambda=\{(x,-x):x\in\R\}$ は閉である。
実際、$(x,y)\notin\Lambda$ とすると、$x+y>0$ なら長方形 $[x,x+1)\times[y,y+1)$ の点の座標和は $x+y$ 以上で正のままであり、
$x+y<0$ なら $\eps=-(x+y)/2$ として $[x,x+\eps)\times[y,y+\eps)$ の点の座標和は $x+y+2\eps=0$ 未満である。
いずれも $\Lambda$ と交わらない開近傍である。
また $(x,-x)$ の近傍 $[x,x+1)\times[-x,-x+1)$ と $\Lambda$ の交わりは、$t\ge x$ かつ $-t\ge-x$ より $t=x$ の一点だけであり、
$\Lambda$ は部分空間として離散である。

$S^2$ が正規なら、:ref[lem-top-jones] により $\Lambda$ の部分集合全体の濃度は $2^{\aleph_0}$ 以下である。
しかし $\Lambda$ は $\R$ と対等で、$\R$ の濃度は $2^{\aleph_0}$ なので、$\Lambda$ の部分集合全体の濃度は $2^{2^{\aleph_0}}$ であり、
Cantor の定理（[濃度](/math/sets-logic/cardinality/)）により $2^{2^{\aleph_0}}>2^{\aleph_0}$ である。
これは矛盾であり、$S^2$ は正規でない。
:::

この例から次の事実が分かる。

- 正規空間（$S$）どうしの積は正規とは限らない。
- Lindelöf 空間どうしの積は Lindelöf とは限らない。もし $S^2$ が Lindelöf なら、正則なので :ref[lem-top-regular-lindelof] により正規になってしまう。
- $T_3$ 空間は $T_4$ とは限らない。

### Urysohn の補題

正規空間で $A$ 上 $0$、$B$ 上 $1$ となる関数を作るには、関数の等高線 $\{u<r\}$ を先に作るのが自然である。
連続関数 $u$ があれば、$U_r=\{u<r\}$ は $r$ とともに増大する開集合で、
$r<s$ なら $\overline{U_r}\subset\{u\le r\}\subset U_s$ を満たす。
逆に、このような開集合の族を二進有理数 $r$ について作れれば、$u(x)$ を「$x$ が初めて入る $U_r$ の添字」として定義できる。
正規性は、閉集合の近傍をその閉包ごと収まるように一段ずつ挿入する操作（:ref[lem-top-shrinking]）を保証する。

$D=\bigcup_{n\ge0}\{k/2^n:0\le k\le2^n\}$ を $[0,1]$ 内の二進有理数全体とする。
$D$ は $[0,1]$ で稠密である。

:::lemma[二進有理数による開集合の挿入]{#lem-top-dyadic-shrinking}
$X$ を正規空間、$A,B$ を交わらない閉集合とする。
このとき開集合の族 $(U_r)_{r\in D}$ で

$$
A\subset U_0,\qquad U_1=X\setminus B,\qquad
r<s\ \Longrightarrow\ \overline{U_r}\subset U_s
$$

を満たすものが存在する。
:::

:::proof
$D_n=\{k/2^n:0\le k\le2^n\}$ とおくと $D_0\subset D_1\subset\cdots$、$D=\bigcup_nD_n$ である。
$n$ に関する帰納法で、$D_n$ を添字とする開集合 $U_r$ を、$r<s$（$r,s\in D_n$）なら $\overline{U_r}\subset U_s$ となるように定め、
$n$ を増やすときには既に定めた集合を変えないようにする。

$n=0$：$U_1=X\setminus B$ とおく。
$A\subset X\setminus B$ で $A$ は閉なので、:ref[lem-top-shrinking] により $A\subset U_0\subset\overline{U_0}\subset U_1$ となる開集合 $U_0$ がある。

$n\to n+1$：$D_{n+1}\setminus D_n$ の元は $t=(2k+1)/2^{n+1}$ の形であり、
$D_n$ の隣り合う二元 $r=k/2^n<s=(k+1)/2^n$ の中点である。
帰納法の仮定より $\overline{U_r}\subset U_s$ なので、閉集合 $\overline{U_r}$ と開集合 $U_s$ に :ref[lem-top-shrinking] を適用して

$$
\overline{U_r}\subset U_t\subset\overline{U_t}\subset U_s
$$

となる開集合 $U_t$ を選ぶ。
これで $D_{n+1}$ 上の族が定まる。
$r'<s'$（$r',s'\in D_{n+1}$）の条件を確かめる。
$D_{n+1}$ の隣り合う二元 $r'<s'$ については、一方が新しい中点で他方が $D_n$ の元なので、上の選び方から $\overline{U_{r'}}\subset U_{s'}$ である。
一般の $r'<s'$ では、その間の $D_{n+1}$ の元を $r'=t_0<t_1<\cdots<t_m=s'$ と並べると、

$$
\overline{U_{t_0}}\subset U_{t_1}\subset\overline{U_{t_1}}\subset U_{t_2}\subset\cdots\subset U_{t_m}
$$

となって $\overline{U_{r'}}\subset U_{s'}$ を得る。

任意の $r<s$（$r,s\in D$）はある $D_n$ に同時に属するので、全体で条件が成り立つ。
:::

:::theorem[Urysohn の補題]{#thm-top-urysohn}
$X$ を正規空間、$A,B$ を交わらない閉集合とする。
このとき連続関数 $u:X\to[0,1]$ で、$A$ 上で $u=0$、$B$ 上で $u=1$ となるものが存在する。
逆に、交わらない任意の閉集合の対に対してこのような関数が存在すれば、$X$ は正規である。
:::

:::proof
:ref[lem-top-dyadic-shrinking] の族 $(U_r)_{r\in D}$ を取り、

$$
u(x)=\inf\bigl(\{r\in D:x\in U_r\}\cup\{1\}\bigr)
$$

と定める。
$u(x)\in[0,1]$ である。

**境界値。**
$x\in A$ なら $x\in U_0$ なので $u(x)=0$ である。
$x\in B$ なら $x\notin U_1=X\setminus B$ であり、$r\le1$ なら $U_r\subset\overline{U_r}\subset U_1$（$r<1$ のとき）または $U_r=U_1$ なので、$x$ はどの $U_r$ にも属さない。
よって $u(x)=\inf\{1\}=1$ である。

**二つの基本的な観察。**
(i) $x\in\overline{U_r}$ なら $u(x)\le r$ である。
実際 $r<s$ となる任意の $s\in D$ について $x\in\overline{U_r}\subset U_s$ なので $u(x)\le s$ であり、
$D$ の稠密性から $s$ を $r$ にいくらでも近く取れる（$r=1$ のときは $u\le1$ が自明である）。
(ii) $x\notin U_r$ なら $u(x)\ge r$ である。
実際 $s<r$ なら $U_s\subset\overline{U_s}\subset U_r$ なので $x\notin U_s$ であり、
$x$ が属する $U_s$ の添字はすべて $r$ 以上だからである。

**$\{u<t\}$ が開であること。**
$0<t\le1$ に対して

$$
\{x:u(x)<t\}=\bigcup_{r\in D,\ r<t}U_r
$$

を示す。
$u(x)<t$ なら、下限の定義から $\{r\in D:x\in U_r\}\cup\{1\}$ に $t$ 未満の元がある。
$t\le1$ なのでそれは $1$ 以外、すなわち $x\in U_r$ となる $r\in D$、$r<t$ である。
逆に $x\in U_r$（$r<t$）なら $u(x)\le r<t$ である。
右辺は開集合の和なので開である。
$t\le0$ なら $\{u<t\}=\varnothing$、$t>1$ なら $\{u<t\}=X$ であり、これらも開である。

**$\{u>t\}$ が開であること。**
$0\le t<1$ に対して

$$
\{x:u(x)>t\}=\bigcup_{r\in D,\ t<r<1}\bigl(X\setminus\overline{U_r}\bigr)
$$

を示す。
右辺の点 $x$ は、ある $r>t$ について $x\notin\overline{U_r}\supset U_r$ を満たすので、(ii) から $u(x)\ge r>t$ である。
逆に $u(x)>t$ なら、$D$ の稠密性から $t<r<s<u(x)$ となる $r,s\in D$ を選べる（$u(x)\le1$ なので $s<1$）。
(i) の対偶から $x\notin\overline{U_s}$、従って $x\notin U_s\supset\overline{U_r}$ であり、$x$ は右辺に属する。
右辺は閉集合の補集合の和なので開である。
$t<0$ なら $\{u>t\}=X$、$t\ge1$ なら $\{u>t\}=\varnothing$ である。

**連続性。**
開半直線 $(-\infty,t)$、$(t,\infty)$ は $\R$ の準開基をなし、その逆像がすべて開であることを示したので、
:ref[prop-top-subbase-continuity] により $u$ は連続である。

**逆。**
そのような $u$ があれば、$\{u<1/2\}$ と $\{u>1/2\}$ が $A,B$ を分離する開集合である。
:::

この関数を :term[Urysohn 関数]{en="Urysohn function" yomi="うりそーんかんすう"} と呼ぶことがある。
構成の仕組みを具体例で確かめておく。

:::example[区間での Urysohn 関数の計算]
$X=[0,1]$、$A=\{0\}$、$B=\{1\}$ とし、二進有理数 $r<1$ に対して $U_r=[0,(1+r)/2)$、$U_1=[0,1)=X\setminus B$ と置く。
$U_0=[0,1/2)\ni0$ であり、$r<s<1$ なら $\overline{U_r}=[0,(1+r)/2]\subset[0,(1+s)/2)=U_s$、
$r<1$ なら $(1+r)/2<1$ より $\overline{U_r}\subset U_1$ なので、補題の条件を満たす。
$x\in[0,1/2)$ なら $x\in U_0$ なので $u(x)=0$ である。
$1/2\le x<1$ なら $x\in U_r$ は $r>2x-1$ と同値であり、
$2x-1$ より大きい二進有理数の下限は $D$ の稠密性から $2x-1$ に等しいので $u(x)=2x-1$ である。
$x=1$ はどの $U_r$ にも属さないので $u(1)=1$ である。
まとめると $u(x)=\max\{0,2x-1\}$ であり、確かに連続である。
集合 $U_r$ の取り方を変えれば別の Urysohn 関数が得られ、Urysohn 関数は一意ではない。
:::

:::example[値が飛ぶように見える Urysohn 関数]
Sorgenfrey 直線 $S$ で $A=(-\infty,0)$、$B=[0,\infty)$ とする。
$A$ の補集合 $B=\bigcup_n[0,n)$ は開なので $A$ は閉であり、$B$ の補集合 $A=\bigcup_n[-n,0)$ も開なので $B$ は閉である。
$B$ の特性関数 $\chi_B$ は、$\R$ の任意の開集合の逆像が $\varnothing,A,B,S$ のいずれかなので $S$ 上で連続であり、Urysohn 関数になる。
通常の直線では閉集合 $(-\infty,0]$ と $[1,\infty)$ を分ける関数は途中のすべての値を取る（[連結性](/math/topology/connectedness/) の中間値の定理）が、
$S$ は連結でないので値の飛びが許される。
:::

:::example[距離が零に近づく閉集合]
平面の $A=\{(x,0):x\ge1\}$ と $B=\{(x,1/x):x\ge1\}$ は交わらない閉集合である。
$x\to\infty$ で両者の点の距離 $1/x$ は $0$ に近づくので、二集合の間の距離の下限は $0$ である。
それでも :ref[prop-top-metric-normal] の関数 $h$ は使える。
各点で $d(\cdot,A)+d(\cdot,B)>0$ であることだけが必要だからである。
正規性は、固定した幅の帯で二つの集合を隔てることを要求していない。
:::

Urysohn の補題の正規性の仮定は外せない。
:ref[exm-top-sorgenfrey-plane] の $S^2$ では、反対角線の部分集合 $E$ と $\Lambda\setminus E$ を分ける連続関数が存在しない $E$ がある
（存在すれば $\{u<1/2\}$、$\{u>1/2\}$ が分離を与えてしまう）。

### 零点集合と完全正規性

Urysohn の補題は $A$ 上で $u=0$ を保証するが、$A$ の外でも $0$ を取ることを許す。
零点集合をちょうど $A$ にできるかどうかは、補集合の各点を可算個の条件で検出できるかによる。

:::proposition[閉集合を正確に零点集合にする条件]{#prop-top-zero-set}
$X$ を正規空間、$A$ を閉集合とする。
ある連続関数 $g:X\to[0,1]$ で $A=Z(g)=g^{-1}(0)$ となるものが存在するための必要十分条件は、
$A$ が可算個の開集合の共通部分（$G_\delta$ 集合）であることである。
:::

:::proof
$A=g^{-1}(0)$ なら $A=\bigcap_{n\ge1}\{x:g(x)<1/n\}$ であり、右辺の各集合は開である。

逆に $A=\bigcap_{n\ge1}U_n$（$U_n$ は開）とする。
各 $n$ で $A$ と $X\setminus U_n$ は交わらない閉集合なので、Urysohn の補題により
$A$ 上 $0$、$X\setminus U_n$ 上 $1$ となる連続関数 $g_n:X\to[0,1]$ がある。
$g=\sum_{n\ge1}2^{-n}g_n$ と置く。
部分和 $s_N=\sum_{n\le N}2^{-n}g_n$ は連続であり、すべての $x$ で

$$
0\le g(x)-s_N(x)\le\sum_{n>N}2^{-n}=2^{-N}
$$

だから $s_N$ は $g$ に一様収束する。
:ref[prop-top-uniform-limit] により $g$ は連続で、値は $[0,\sum_n2^{-n}]=[0,1]$ に入る。
$x\in A$ ならすべての $g_n(x)=0$ なので $g(x)=0$ である。
$x\notin A$ ならある $n$ で $x\notin U_n$ であり、$g_n(x)=1$ から $g(x)\ge2^{-n}>0$ である。
よって $Z(g)=A$ となる。
:::

:::definition[完全正規]{#def-top-perfectly-normal}
正規空間で、すべての閉集合が $G_\delta$ 集合（同値な言い方では、ある連続関数の零点集合）であるものを
:term[完全正規]{en="perfectly normal" yomi="かんぜんせいき"} であるという。
:::

:ref[prop-top-metric-normal] により距離空間は完全正規である。

:::example[完全正規でない $T_4$ 空間]
非可算な離散空間 $X$ の一点コンパクト化 $\alpha X=X\cup\{\infty\}$（:ref[thm-top-one-point]）を考える。
$\alpha X$ はコンパクト Hausdorff なので :ref[prop-top-hausdorff-product] により $T_4$ である。
$X$ のコンパクト集合は有限集合なので、$\infty$ の開近傍は有限集合の補集合を含む。
可算個の開近傍 $W_1,W_2,\ldots$ の共通部分の補集合は、各 $W_n$ の補集合（有限集合）の可算和なので高々可算であり、
$\bigcap_nW_n$ は $X$ の非可算個の点を含む。
よって閉集合 $\{\infty\}$ は $G_\delta$ でなく、:ref[prop-top-zero-set] により $\{\infty\}$ をちょうど零点集合とする連続関数は存在しない。
:::

### 有限開被覆に従う一の分割

局所的に定義された量を、重みを付けて平均することで大域的な量に貼り合わせる操作は、
解析学や多様体論の基本的な技法である。
重みとして使う関数を正規性から作る。

:::proposition[正規空間での有限な一の分割]{#prop-top-finite-partition}
正規空間 $X$ を有限個の開集合 $U_1,\ldots,U_m$ で覆うとき、連続関数 $\phi_i:X\to[0,1]$ で

$$
\sum_{i=1}^m\phi_i=1,\qquad
\operatorname{supp}\phi_i\subset U_i\quad(1\le i\le m)
$$

を満たすものが存在する。
ここで $\operatorname{supp}\phi=\overline{\{x:\phi(x)\ne0\}}$ は $\phi$ の :term[台]{en="support" yomi="だい"} である。
このような族を被覆 $(U_i)$ に従う :term[一の分割]{en="partition of unity" yomi="いちのぶんかつ"} という。
:::

:::proof
**被覆の縮小。**
$i=1,\ldots,m$ の順に、開集合 $V_i$ で $\overline{V_i}\subset U_i$ かつ
$V_1,\ldots,V_i,U_{i+1},\ldots,U_m$ が $X$ を覆うものを選ぶ。
$V_1,\ldots,V_{i-1},U_i,\ldots,U_m$ が $X$ を覆っているとき、

$$
F_i=X\setminus\bigl(V_1\cup\cdots\cup V_{i-1}\cup U_{i+1}\cup\cdots\cup U_m\bigr)
$$

は閉集合であり、被覆の仮定から $F_i\subset U_i$ である。
:ref[lem-top-shrinking] により $F_i\subset V_i\subset\overline{V_i}\subset U_i$ となる開集合 $V_i$ を取る。
$F_i$ の点は $V_i$ に、それ以外の点は $V_1,\ldots,V_{i-1},U_{i+1},\ldots,U_m$ のどれかに属するので、
$U_i$ を $V_i$ に置き換えても $X$ を覆う。
$m$ 回で $V_1,\ldots,V_m$ が $X$ を覆い、$\overline{V_i}\subset U_i$ となる。

**関数の構成。**
各 $i$ でもう一度 :ref[lem-top-shrinking] を使い、$\overline{V_i}\subset W_i\subset\overline{W_i}\subset U_i$ となる開集合 $W_i$ を取る。
Urysohn の補題により、$\overline{V_i}$ 上で $1$、$X\setminus W_i$ 上で $0$ となる連続関数 $g_i:X\to[0,1]$ がある。
$\{g_i\ne0\}\subset W_i$ なので $\operatorname{supp}g_i\subset\overline{W_i}\subset U_i$ である。
各点はある $V_j$ に属し、そこで $g_j=1$ なので $\sum_jg_j\ge1$ である。
$\phi_i=g_i/\sum_jg_j$ は分母が正の連続関数の商として連続で、値は $[0,1]$ に入り、$\sum_i\phi_i=1$ を満たす。
$\{\phi_i\ne0\}=\{g_i\ne0\}$ なので $\operatorname{supp}\phi_i=\operatorname{supp}g_i\subset U_i$ である。
:::

:::example[重なる二区間で重みを渡す]
$[0,1]$ の開被覆 $U=[0,3/4)$、$V=(1/4,1]$ を取る。
$A=[0,1/3]$、$B=[2/3,1]$ に対する :ref[prop-top-metric-normal] の関数を $\phi_V$ とすると、

$$
\phi_V(x)=\begin{cases}0&(0\le x\le1/3)\\3x-1&(1/3\le x\le2/3)\\1&(2/3\le x\le1)\end{cases}
$$

である。
$\phi_U=1-\phi_V$ と置くと、$\operatorname{supp}\phi_U=[0,2/3]\subset U$、$\operatorname{supp}\phi_V=[1/3,1]\subset V$ である。
例えば $x=1/2$ では $\phi_U=\phi_V=1/2$、$x=0.6$ では $\phi_V=0.8$、$\phi_U=0.2$ となる。
$U$ 上の量 $a$ と $V$ 上の量 $b$ を $\phi_Ua+\phi_Vb$ で貼り合わせれば、重なりの領域で連続的に切り替わる。折れ点では一般に微分できない。
:::

## Tietze の拡張定理

Urysohn の補題は「二つの閉集合の上で値 $0,1$ を指定した関数」を作った。
本節ではこれを、閉集合 $A$ 上で与えられた任意の連続関数 $f$ を全空間の連続関数に延長する問題へ一般化する。
一度に一致させるのは難しいので、Urysohn 関数で $f$ を粗く近似し、
残った誤差を毎回一定の割合で小さくしていく。
誤差の幾何級数の和として延長を得るのである。

:::lemma[誤差を三分の二にする近似]{#lem-top-tietze-step}
$X$ を正規空間、$A$ を閉集合、$M\ge0$ とし、$g:A\to[-M,M]$ を連続とする。
このとき連続関数 $h:X\to[-M/3,M/3]$ で、$A$ 上で $\abs{g-h}\le2M/3$ を満たすものが存在する。
:::

:::proof
$M=0$ なら $h=0$ でよい。
$M>0$ とし、

$$
L=\{a\in A:g(a)\le-M/3\},\qquad H=\{a\in A:g(a)\ge M/3\}
$$

と置く。
$L,H$ は連続関数 $g$ による閉集合の逆像なので $A$ で閉であり、$A$ が閉なので $X$ でも閉である。
また明らかに交わらない（空でもよい）。
Urysohn の補題で $L$ 上 $0$、$H$ 上 $1$ となる連続関数 $u:X\to[0,1]$ を取り、

$$
h=\frac{M}{3}(2u-1)
$$

と置く。
$h$ は連続で $[-M/3,M/3]$ に値を取る。
$a\in L$ なら $h(a)=-M/3$、$g(a)\in[-M,-M/3]$ なので $\abs{g(a)-h(a)}\le2M/3$ である。
$a\in H$ なら同様である。
$a\in A\setminus(L\cup H)$ なら $\abs{g(a)}<M/3$ かつ $\abs{h(a)}\le M/3$ なので $\abs{g(a)-h(a)}<2M/3$ である。
:::

:::theorem[Tietze の拡張定理]{#thm-top-tietze}
$X$ を正規空間、$A$ を閉集合とする。

1. 連続関数 $f:A\to[a,b]$（$a\le b$）は、連続関数 $F:X\to[a,b]$ に延長できる。
2. 連続関数 $f:A\to\R$ は、連続関数 $F:X\to\R$ に延長できる。
3. 連続関数 $f:A\to(-1,1)$ は、連続関数 $F:X\to(-1,1)$ に延長できる。
:::

:::proof
**有界な場合。**
まず $f:A\to[-M,M]$ とする。
$g_0=f$、$M_k=M(2/3)^k$ と置き、$k=0,1,2,\ldots$ について次を帰納的に行う。
$g_k:A\to[-M_k,M_k]$ が連続のとき、:ref[lem-top-tietze-step] により連続関数 $h_k:X\to[-M_k/3,M_k/3]$ で
$A$ 上 $\abs{g_k-h_k}\le2M_k/3=M_{k+1}$ となるものを取り、$g_{k+1}=g_k-h_k|_A$ と置く。
$g_{k+1}$ は連続で $[-M_{k+1},M_{k+1}]$ に値を取るので、帰納法が続く。

$F_n=\sum_{k=0}^{n-1}h_k$ と置く。
$n<m$ ならすべての $x\in X$ で

$$
\abs{F_m(x)-F_n(x)}\le\sum_{k=n}^{m-1}\frac{M_k}{3}\le\sum_{k=n}^{\infty}\frac{M}{3}\Bigl(\frac23\Bigr)^k=M\Bigl(\frac23\Bigr)^n
$$

なので、各 $x$ で $(F_n(x))$ は Cauchy 列であり極限 $F(x)$ をもつ。
上の不等式で $m\to\infty$ とすると $\abs{F(x)-F_n(x)}\le M(2/3)^n$ がすべての $x$ で成り立ち、$F_n$ は $F$ に一様収束する。
:ref[prop-top-uniform-limit] により $F$ は連続である。
$A$ 上では $g_n=f-F_n|_A$ が帰納的に成り立つので $\abs{f-F_n}\le M_n\to0$ となり、$F|_A=f$ である。
また $\abs{F}\le\sum_{k\ge0}M_k/3=M$ である。

**主張 1。**
$a=b$ なら定数関数でよい。
$a<b$ なら $c=(a+b)/2$、$M=(b-a)/2$ として $f-c:A\to[-M,M]$ に上の結果を適用し、
得た延長に $c$ を加えれば $[a,b]$ に値を取る延長が得られる。

**主張 3。**
$f:A\to(-1,1)$ とする。
主張 1 により $f$ を $H:X\to[-1,1]$ に延長する。
$E=H^{-1}(\{-1,1\})$ は閉集合であり、$A$ 上では $H=f\in(-1,1)$ なので $E\cap A=\varnothing$ である。
Urysohn の補題で $A$ 上 $1$、$E$ 上 $0$ となる連続関数 $v:X\to[0,1]$ を取り、$F=vH$ と置く。
$F$ は連続で $A$ 上 $F=f$ である。
$x\in E$ なら $F(x)=0$、$x\notin E$ なら $\abs{F(x)}\le\abs{H(x)}<1$ なので、$F$ は $(-1,1)$ に値を取る。

**主張 2。**
$\psi(t)=t/(1+\abs t)$ は $\R$ から $(-1,1)$ への同相写像であり、逆写像は $s\mapsto s/(1-\abs s)$ である。
主張 3 により $\psi\circ f:A\to(-1,1)$ を $G:X\to(-1,1)$ に延長し、$F=\psi^{-1}\circ G$ と置けばよい。
:::

:::corollary[立方体や Euclid 空間への写像の延長]{#cor-top-tietze-product}
$X$ を正規空間、$A$ を閉集合とし、各 $j\in J$ について $I_j$ を $\R$ または有界閉区間とする。
連続写像 $f:A\to\prod_{j\in J}I_j$ は連続写像 $F:X\to\prod_{j\in J}I_j$ に延長できる。
特に $A$ 上の $\R^n$ 値または $[0,1]^n$ 値の連続写像は、同じ空間に値をもつ連続写像に延長できる。
:::

:::proof
:ref[thm-top-product-property] により各座標 $f_j=\pi_j\circ f:A\to I_j$ は連続である。
:ref[thm-top-tietze] により各 $j$ で延長 $F_j:X\to I_j$ を選ぶ（$J$ が無限なら選択公理を使う）。
$F=(F_j)_{j\in J}$ は再び :ref[thm-top-product-property] により連続で、$A$ 上で $f$ に一致する。
:::

:::example[具体的な延長]
$X=\R$、$A=(-\infty,-1]\cup[1,\infty)$ とし、$f$ を $(-\infty,-1]$ 上で $-2$、$[1,\infty)$ 上で $3$ とする。
$A$ は閉で、$f$ は二つの閉集合の上で連続なので :ref[lem-top-pasting] により $A$ 上連続である。
延長として

$$
F(x)=\begin{cases}-2&(x\le-1)\\\dfrac{5x+1}{2}&(-1\le x\le1)\\3&(x\ge1)\end{cases}
$$

が取れる。
$x=\pm1$ で両側の式が一致する（$(5\cdot(-1)+1)/2=-2$、$(5+1)/2=3$）ので、$F$ は貼り合わせの補題により連続である。
$F$ は $[-2,3]$ に値を取り、定理の主張 1 と整合する。
:::

延長定理の仮定はどれも省けない。

:::example[延長できない状況]
**定義域が閉でない場合。**
$X=[0,1]$、$A=(0,1]$、$f(x)=\sin(1/x)$ とする。
$f$ は有界連続であるが、$x_n=1/(\pi/2+2\pi n)$、$y_n=1/(3\pi/2+2\pi n)$ はともに $0$ に収束し、$f(x_n)=1$、$f(y_n)=-1$ である。
連続な延長 $F$ があれば $F(0)=\lim F(x_n)=1$ かつ $F(0)=\lim F(y_n)=-1$ となって矛盾する。
定理の証明では、$\{g\le-M/3\}$ と $\{g\ge M/3\}$ が全空間で閉であることに $A$ の閉性を使っていた。

**値域が区間でない場合。**
$X=[0,1]$、$A=\{0,1\}$ とし、$f:A\to\{0,1\}$ を恒等写像とする（$\{0,1\}$ は $\R$ の部分空間として離散）。
連続な延長 $F:[0,1]\to\{0,1\}$ があれば、:ref[prop-top-connected-image] により連結な $[0,1]$ の像 $\{0,1\}$ は連結でなければならないが、これは誤りである。
Tietze の定理は値域が $\R$ や区間であることに本質的に依存している。
:::

:::example[延長性から正規性を取り戻す]
逆に、$X$ の任意の閉集合上の任意の連続関数 $A\to[-1,1]$ が $X$ 上に延長できると仮定する。
交わらない閉集合 $A,B$ を取り、$A\cup B$ 上で $A$ では $-1$、$B$ では $1$ を取る関数 $f$ を考える。
$A,B$ は $A\cup B$ の閉集合で $A\cup B$ を覆うので、:ref[lem-top-pasting] により $f$ は連続である。
延長 $F:X\to[-1,1]$ を取れば $\{F<0\}$、$\{F>0\}$ が $A,B$ を分離する。
従って Tietze の拡張定理の結論は正規性と同値である。
:::

## 距離化可能性

:::definition[距離化可能性]{#def-top-metrizable}
位相空間 $(X,\mathcal T)$ に対し、ある距離 $d$ の定める位相が $\mathcal T$ と一致するとき、
$X$ は :term[距離化可能]{en="metrizable" yomi="きょりかかのう"} であるという。
:::

距離化可能性は位相的性質であり、同相写像で保たれる。
どの空間が距離化可能かを判定するには、まず距離空間が必ずもつ性質を調べ、
次にそれらから距離を構成する方法を探すのが自然である。

:::proposition[距離化可能であるための必要条件]{#prop-top-metrizable-necessary}
距離化可能な空間 $X$ は次の性質をもつ。

1. 完全正規な $T_4$ 空間である。
2. 第一可算である。
3. 任意の部分空間も距離化可能である。
4. 可分であること、第二可算であること、Lindelöf であることは互いに同値である。
:::

:::proof
1 は :ref[prop-top-metric-normal] である。
2 と、4 のうち可分性と第二可算性の同値は :ref[prop-top-separable-metric] による。
3 は距離を部分集合に制限すればよい（相対位相は制限した距離の位相に一致する）。
4 の残りを示す。
第二可算なら :ref[prop-top-countable-basis] により Lindelöf である。
逆に Lindelöf なら、各 $n$ について開被覆 $\{B(x,1/n):x\in X\}$ の可算部分被覆の中心の集合 $C_n$ を取る。
$C=\bigcup_nC_n$ は可算であり、任意の開球 $B(y,r)$ に対し $1/n<r$ とすれば $y\in B(c,1/n)$ となる $c\in C_n$ があり、$c\in B(y,r)$ である。
よって $C$ は稠密であり、$X$ は可分である。
:::

:::example[Sorgenfrey 直線は距離化可能でない]
Sorgenfrey 直線 $S$ は :ref[exm-top-sorgenfrey-normal] により $T_4$ であり、
$[x,x+1/n)$ を基本近傍系にもつので第一可算、$\Q$ が稠密なので可分、さらに Lindelöf である。
しかし :ref[prop-top-sorgenfrey-lindelof] のとおり第二可算でない。
距離化可能なら :ref[prop-top-metrizable-necessary] の 4 により第二可算になるはずなので、$S$ は距離化可能でない。
上の必要条件のうち 1 と 2 をすべて満たしても距離化可能とは限らないことが分かる。
:::

### 完全正則空間

距離化の基本的なアイデアは、点を区別する連続関数を座標として並べ、
空間を距離空間（例えば立方体 $[0,1]^{\N}$）の部分空間として実現することである。
そのためには「点と閉集合を連続関数で分ける」ことができなければならない。

:::definition[完全正則空間]{#def-top-completely-regular}
$X$ が :term[完全正則]{en="completely regular" yomi="かんぜんせいそく"} であるとは、
任意の閉集合 $F$ と $x\notin F$ に対し、連続関数 $f:X\to[0,1]$ で $f(x)=1$、$F$ 上で $f=0$ となるものが存在することである。
完全正則な $T_1$ 空間を $T_{3\frac12}$ 空間または :term[Tychonoff 空間]{en="Tychonoff space" yomi="ちほのふくうかん"} という。
:::

:::proposition[完全正則性の基本性質]{#prop-top-completely-regular-basic}
1. $T_4\Rightarrow T_{3\frac12}\Rightarrow T_3$ が成り立つ。特に距離化可能な空間とコンパクト Hausdorff 空間は Tychonoff 空間である。
2. 完全正則空間の部分空間は完全正則であり、Tychonoff 空間の部分空間は Tychonoff 空間である。
3. 完全正則空間の任意個の積は完全正則であり、Tychonoff 空間の任意個の積は Tychonoff 空間である。
:::

:::proof
**1。**
$X$ を $T_4$ とし、閉集合 $F$ と $x\notin F$ を取る。
$\{x\}$ は閉なので、Urysohn の補題を閉集合 $F$ と $\{x\}$ に使えば $F$ 上 $0$、$x$ で $1$ の連続関数が得られる。
次に $X$ を完全正則とし、$f(x)=1$、$F$ 上 $f=0$ とすると、$\{f>1/2\}$ と $\{f<1/2\}$ が $x$ と $F$ を分離する。
よって完全正則なら正則であり、$T_1$ を付ければ $T_{3\frac12}\Rightarrow T_3$ である。
距離化可能な空間とコンパクト Hausdorff 空間は :ref[prop-top-metric-normal] と :ref[prop-top-hausdorff-product] により $T_4$ である。

**2。**
$A\subset X$、$F$ を $A$ の閉集合、$a\in A\setminus F$ とする。
$F=A\cap F'$（$F'$ は $X$ で閉）と書くと $a\notin F'$ である。
$X$ の連続関数 $f$ で $f(a)=1$、$F'$ 上 $0$ となるものを取れば、制限 $f|_A$ が求める関数である。
$T_1$ も部分空間に受け継がれる。

**3。**
$P=\prod_{i\in I}X_i$、$p\in P$、$F$ を $p$ を含まない閉集合とする。
$G=P\setminus F$ は $p$ の開近傍なので、有限集合 $J\subset I$ と開集合 $G_i\ni p_i$（$i\in J$）で
$p\in\bigcap_{i\in J}\pi_i^{-1}(G_i)\subset G$ となるものがある。
各 $i\in J$ で $X_i\setminus G_i$ は $p_i$ を含まない閉集合なので、
連続関数 $f_i:X_i\to[0,1]$ で $f_i(p_i)=1$、$X_i\setminus G_i$ 上 $0$ となるものがある。

$$
f(y)=\min_{i\in J}f_i(y_i)\qquad(J=\varnothing\text{ なら }f\equiv1)
$$

と置く。
各 $y\mapsto f_i(y_i)=f_i(\pi_i(y))$ は連続で、有限個の連続関数の最小値は $\min\{s,t\}=(s+t-\abs{s-t})/2$ から連続である。
$f(p)=1$ である。
$y\in F$ なら $y\notin G$ なので、ある $i\in J$ で $y_i\notin G_i$ となり $f_i(y_i)=0$、従って $f(y)=0$ である。
$J=\varnothing$ の場合は $G=P$ なので $F=\varnothing$ であり、$f\equiv1$ で条件を満たす。
$T_1$ 性の積への遺伝は :ref[prop-top-regular-product] の証明で示した。
:::

局所コンパクト Hausdorff 空間は一般に正規とは限らないが、常に Tychonoff 空間である。

:::proposition[局所コンパクト Hausdorff 空間の完全正則性]{#prop-top-lch-completely-regular}
局所コンパクト Hausdorff 空間は Tychonoff 空間である。
:::

:::proof
$X$ を局所コンパクト Hausdorff 空間とする。
$X$ がコンパクトなら :ref[prop-top-completely-regular-basic] の 1 による。
$X$ がコンパクトでないとき、:ref[thm-top-one-point] により一点コンパクト化 $\alpha X$ はコンパクト Hausdorff 空間であり、$X$ はその部分空間である。
$\alpha X$ は Tychonoff 空間なので、:ref[prop-top-completely-regular-basic] の 2 により $X$ も Tychonoff 空間である。
（:ref[prop-top-lch-subspace] のとおり、局所コンパクト Hausdorff 空間はちょうどコンパクト Hausdorff 空間の開部分空間であり、この証明はその一方向を使っている。）
:::

### 立方体への埋め込み

座標として使う関数の族が十分豊富なら、評価写像は埋め込みになる。
この補題は Tychonoff の埋め込み定理と Urysohn の距離化定理の共通の核である。

:::lemma[関数族による埋め込み]{#lem-top-embedding-functions}
$X$ を $T_1$ 空間とし、連続関数の族 $f_j:X\to[0,1]$（$j\in J$）が次の条件を満たすとする。

> 任意の開集合 $G$ と $x\in G$ に対し、ある $j\in J$ で $f_j(x)>0$ かつ $X\setminus G$ 上 $f_j=0$ となる。

このとき $e:X\to[0,1]^J$、$e(x)=(f_j(x))_{j\in J}$ は埋め込みである。
:::

:::proof
**連続性。**
各座標 $\pi_j\circ e=f_j$ が連続なので、:ref[thm-top-product-property] により $e$ は連続である。

**単射性。**
$x\ne y$ とすると、$T_1$ 性から $G=X\setminus\{y\}$ は $x$ を含む開集合である。
条件から $f_j(x)>0$、$f_j(y)=0$ となる $j$ があり、$e(x)\ne e(y)$ である。

**像への開写像であること。**
$G\subset X$ を開集合とし、$e(G)$ が $e(X)$ で開であることを示す。
$x\in G$ に対し条件の $j=j(x)$ を取り、$O_x=\pi_j^{-1}\bigl((0,1]\bigr)$ と置く。
$(0,1]$ は $[0,1]$ の開集合なので $O_x$ は $[0,1]^J$ の開集合であり、$e(x)\in O_x$ である。
$y\in X$ で $e(y)\in O_x$ なら $f_j(y)>0$ なので $y\in G$ である。
従って $O_x\cap e(X)\subset e(G)$ であり、

$$
e(G)=\bigcup_{x\in G}\bigl(O_x\cap e(X)\bigr)
$$

は $e(X)$ の開集合である。
$e:X\to e(X)$ は連続な全単射で開写像だから、:ref[prop-top-homeomorphism-test] により同相写像である。
:::

:::theorem[Tychonoff の埋め込み定理]{#thm-top-tychonoff-embedding}
位相空間 $X$ について、次は同値である。

1. $X$ は Tychonoff 空間である。
2. $X$ はある集合 $J$ に対する立方体 $[0,1]^J$ の部分空間と同相である。
3. $X$ はあるコンパクト Hausdorff 空間の部分空間と同相である。
:::

:::proof
**1 から 2。**
$J=C(X,[0,1])$（$X$ から $[0,1]$ への連続関数全体）とし、各 $f\in J$ を座標関数とする。
$x\in G$（$G$ は開）なら、$X\setminus G$ は $x$ を含まない閉集合なので、完全正則性から $f(x)=1$、$X\setminus G$ 上 $f=0$ となる $f\in J$ がある。
:ref[lem-top-embedding-functions] により $e(x)=(f(x))_{f\in J}$ は埋め込みである。

**2 から 3。**
$[0,1]$ は :ref[thm-top-interval-compact] によりコンパクトで、距離空間なので Hausdorff である。
:ref[thm-top-tychonoff] と :ref[prop-top-hausdorff-product] により $[0,1]^J$ はコンパクト Hausdorff 空間である。

**3 から 1。**
コンパクト Hausdorff 空間は :ref[prop-top-completely-regular-basic] の 1 により Tychonoff 空間であり、
その部分空間は同じ命題の 2 により Tychonoff 空間である。
Tychonoff 空間であることは同相で保たれる。
:::

この定理により、Tychonoff 空間は「コンパクト Hausdorff 空間の中に住む空間」として特徴付けられる。
一方で、コンパクト Hausdorff 空間は正規であるのに、その部分空間は正規とは限らない。

:::example[正規性は部分空間に遺伝しない]{#exm-top-normal-not-hereditary}
Sorgenfrey 平面 $S^2$ は Tychonoff 空間である。
実際、Sorgenfrey 直線 $S$ は :ref[exm-top-sorgenfrey-normal] により $T_4$ なので :ref[prop-top-completely-regular-basic] の 1 により Tychonoff 空間であり、
同じ命題の 3 により $S^2$ も Tychonoff 空間である。
:ref[thm-top-tychonoff-embedding] により $S^2$ は立方体 $[0,1]^J$ の部分空間 $e(S^2)$ と同相である。
$[0,1]^J$ はコンパクト Hausdorff なので :ref[prop-top-hausdorff-product] により正規であるが、
$e(S^2)\cong S^2$ は :ref[exm-top-sorgenfrey-plane] により正規でない。
従って正規空間の部分空間は正規とは限らない。
:::

### Urysohn の距離化定理

立方体 $[0,1]^J$ は $J$ が可算なら距離化可能である（:ref[prop-top-product-metric]）。
従って可算個の関数で :ref[lem-top-embedding-functions] の条件を満たせれば、距離が得られる。
可算個で済ませるために、可算開基の元の対を使う。

:::theorem[Urysohn の距離化定理]{#thm-top-urysohn-metrization}
第二可算な $T_3$ 空間は距離化可能であり、可算積 $[0,1]^{\N}$（Hilbert 立方体）の部分空間と同相である。
:::

:::proof
$X$ を第二可算な $T_3$ 空間とする。
$X=\varnothing$ なら空の距離でよいので、$X\ne\varnothing$ とする。

**正規性。**
:ref[cor-top-second-countable-normal] により $X$ は $T_4$ である。

**可算個の関数。**
可算開基 $\mathcal B=\{B_1,B_2,\ldots\}$ を取る。
$\overline{B_i}\subset B_j$ を満たす添字の対 $(i,j)$ 全体を $\Pi$ とする。
$(i,j)\in\Pi$ ごとに、交わらない閉集合 $\overline{B_i}$ と $X\setminus B_j$ に Urysohn の補題を使い、
$\overline{B_i}$ 上 $1$、$X\setminus B_j$ 上 $0$ となる連続関数 $f_{ij}:X\to[0,1]$ を選ぶ。
$\Pi\subset\N\times\N$ は高々可算である。

**埋め込みの条件。**
開集合 $G$ と $x\in G$ を取る。
開基の性質から $x\in B_j\subset G$ となる $j$ がある。
:ref[lem-top-shrinking] により $x\in W\subset\overline W\subset B_j$ となる開集合 $W$ を取り、さらに $x\in B_i\subset W$ となる $i$ を取る。
$\overline{B_i}\subset\overline W\subset B_j$ なので $(i,j)\in\Pi$ であり、
$f_{ij}(x)=1>0$、$X\setminus G\subset X\setminus B_j$ 上で $f_{ij}=0$ である。
特に $\Pi$ は空でない。

**埋め込みと距離。**
$\Pi$ の元を一列に並べ（有限なら最後の元を繰り返して）、対応する関数を $f_1,f_2,\ldots$ とする。
上で示した条件は :ref[lem-top-embedding-functions] の仮定そのものなので、
$e:X\to[0,1]^{\N}$、$e(x)=(f_n(x))_{n\ge1}$ は埋め込みである。
:ref[prop-top-product-metric] により $[0,1]^{\N}$ の積位相は距離 $D(s,t)=\sum_n2^{-n}\min\{1,\abs{s_n-t_n}\}$ で与えられ、
各座標の差は $1$ 以下なので $\min$ は不要である。
$e$ は単射なので

$$
d(x,y)=D(e(x),e(y))=\sum_{n=1}^{\infty}2^{-n}\abs{f_n(x)-f_n(y)}
$$

は $X$ 上の距離である。
$e$ が $X$ から $e(X)$ への同相写像で、$e(X)$ の相対位相は $D$ を制限した距離の位相に一致するので、$d$ は $X$ の位相を与える。
:::

この定理の仮定のうち、$T_3$ を正規性に強めることはできても弱めることはできない（距離空間は正則である）。
第二可算性は十分条件であるが必要条件ではない。
例えば非可算な離散空間は離散距離で距離化できるが、第二可算でない。

:::corollary[コンパクト Hausdorff 空間の距離化]{#cor-top-compact-metrizable}
コンパクト Hausdorff 空間 $X$ が距離化可能であるための必要十分条件は、$X$ が第二可算であることである。
:::

:::proof
$X$ が第二可算なら、:ref[prop-top-hausdorff-product] により $T_4$、従って $T_3$ なので、Urysohn の距離化定理により距離化可能である。
逆に $X$ がコンパクトで距離 $d$ をもつとする。
:ref[thm-top-metric-compactness] により $X$ は全有界なので、各 $n$ に有限集合 $C_n$ で $X=\bigcup_{c\in C_n}B(c,1/n)$ となるものがある。
$C=\bigcup_nC_n$ は可算であり、:ref[prop-top-metrizable-necessary] の 4 の証明と同様に稠密である。
よって $X$ は可分であり、:ref[prop-top-separable-metric] により第二可算である。
:::

:::example[距離化できないコンパクト Hausdorff 空間]
**非可算個の積。**
$J$ を非可算集合とすると、$[0,1]^J$ はコンパクト Hausdorff であるが、どの点でも第一可算でない。
実際、点 $x$ の可算個の開近傍 $W_1,W_2,\ldots$ を任意に取る。
各 $W_n$ は基本開集合 $B_n\ni x$ を含み、$B_n$ が制限する座標は有限集合 $J_n$ である。
$\bigcup_nJ_n$ は可算なので、それに属さない $j\in J$ がある。
開近傍 $O=\pi_j^{-1}\bigl((x_j-1/2,x_j+1/2)\cap[0,1]\bigr)$ を考える。
$B_n$ は座標 $j$ を制限しないので、$x$ の座標 $j$ だけを $\abs{y_j-x_j}\ge1/2$ となる値（$x_j\le1/2$ なら $1$、そうでなければ $0$）に替えた点 $y$ は $B_n\subset W_n$ に属し、$O$ には属さない。
従ってどの $W_n$ も $O$ に含まれず、$x$ は可算基本近傍系をもたない。
:ref[prop-top-metrizable-necessary] の 2 により $[0,1]^J$ は距離化可能でない。

**非可算離散空間の一点コンパクト化。**
非可算離散空間 $X$ の $\alpha X$ はコンパクト Hausdorff であるが、各点 $x\in X$ について $\{x\}$ が開であるため、
任意の開基は非可算個の一点集合をすべて含み、第二可算でない。
:ref[cor-top-compact-metrizable] により $\alpha X$ は距離化可能でない（完全正規でないことからも分かる）。
:::

:::theorem[Nagata–Smirnov の距離化定理]{#thm-top-nagata-smirnov}
$T_3$ 空間 $X$ が距離化可能であることと、局所有限な開集合族 $\mathcal B_n$（$n\ge1$）の和を開基にもつことは同値である。
:::

:::proof
距離がある場合は :ref[prop-top-metric-discrete-base] により、各段を離散な族にさえ取れる。逆に仮定の開基を取り、任意の開集合 $G$ に対し $P_n(G)$ を、$\mathcal B_n$ の元 $B$ で $\overline B\subset G$ を満たすものの和とする。正則性により $x\in G$ を含む開近傍を閉包ごと $G$ の中に取り、その中に開基の元を取れるので $G=\bigcup_n P_n(G)$ である。局所有限性と :ref[prop-top-locally-finite-closure] により $\overline{P_n(G)}\subset G$ である。従ってすべての開集合が可算個の閉集合の和になる。

次に、互いに素な閉集合 $A,F$ に対し $P_n=P_n(X\setminus F)$、$Q_n=P_n(X\setminus A)$ と置く。開集合
$$
U=\bigcup_n\left(P_n\setminus\bigcup_{j\le n}\overline{Q_j}\right),\qquad
V=\bigcup_n\left(Q_n\setminus\bigcup_{j\le n}\overline{P_j}\right)
$$
はそれぞれ $A,F$ を含む。$U$ の第 $n$ 項と $V$ の第 $m$ 項は、$n\ge m$ なら前者が $Q_m$ を避け、$m\ge n$ なら後者が $P_n$ を避けるので交わらない。これで $X$ は正規である。

各開基の元 $B$ の補集合は上の結果から $G_\delta$ である。:ref[prop-top-zero-set] により連続関数 $f_B:X\to[0,1]$ を $\{f_B>0\}=B$ となるように選べる。空の上限をここでは零と約束し、
$$
p_n(x,y)=\sup_{B\in\mathcal B_n}|f_B(x)-f_B(y)|,\qquad
d(x,y)=\sum_{n\ge1}2^{-n}p_n(x,y)
$$
と置く。各 $p_n$ は三角不等式と対称性を満たし、$0\le p_n\le1$ である。$x\ne y$ なら $T_1$ 性から $x\in B\subset X\setminus\{y\}$ なる開基の元が取れるので $d(x,y)>0$、従って $d$ は距離である。

固定した $x,n$ に対し、$x$ の近傍で値が零以外になりうる $f_B$ は局所有限性により有限個である。そのため $p_n(x,y)$ はその近傍で有限個の連続関数と零の最大値になり、$y=x$ で連続である。級数の末尾は $2^{-N}$ 以下なので、$x$ のどの $d$-球も元の位相での $x$ の近傍を含む。逆に $x\in B\subset G$、$B\in\mathcal B_n$ と取れば、$y\notin B$ では $d(x,y)\ge2^{-n}f_B(x)>0$ である。その半径の開球は $G$ 内にあり、二つの位相は一致する。
:::

### T₃ だが完全正則でない空間

$T_{3\frac12}$ は $T_3$ と $T_4$ の間の条件である。
$T_4$ でない $T_{3\frac12}$ 空間は Sorgenfrey 平面であった。
$T_3$ でありながら $T_{3\frac12}$ でない空間も存在する。
ここでは Mysior による比較的簡単な例を、すべての主張を確かめながら構成する。
この例は「点と閉集合を開集合で分離できても、連続関数では分離できない」という現象を示す。

:::example[Mysior の平面]{#exm-top-mysior}
上半平面 $H=\{(x,y)\in\R^2:y\ge0\}$ に新しい一点 $a$ を加えた $X=H\cup\{a\}$ を考える。
$x$ 軸上の点 $z=(s,0)$ に対し、縦線分と斜め線分

$$
I_s=\{(s,y):0\le y<2\},\qquad J_s=\{(s+y,y):0\le y<2\}
$$

を定める。
次の集合を各点の基本近傍とする。

- $y>0$ の点 $p$：$\{p\}$（孤立点とする）。
- $z=(s,0)$：$(I_s\cup J_s)\setminus E$（$E$ は $z$ を含まない有限集合）。
- $a$：$N_n=\{a\}\cup\{(x,y)\in H:x>n\}$（$n\in\N$）。

$U\subset X$ が開であるとは、$U$ の各点について、その基本近傍の一つが $U$ に含まれることと定める。
これは位相である。
各基本近傍は開集合である。
実際、$(I_s\cup J_s)\setminus E$ の $z$ 以外の点は $y>0$ の孤立点であり、
$N_n$ の点 $(t,0)$（$t>n$）については $I_t\cup J_t\subset\{x\ge t\}\subset N_n$ である。
また二つの基本近傍の共通部分は同じ点の基本近傍を含む。

**$T_1$ であること。**
$p\ne q$ とする。
$p$ が孤立点なら $\{p\}$ が $q$ を含まない開集合である。
$p=(s,0)$ なら、$E=\{q\}\cap(I_s\cup J_s)$ とした $(I_s\cup J_s)\setminus E$ が $q$ を含まない。
$p=a$ なら、$q$ の第一座標より大きい $n$ を取れば $N_n$ は $q$ を含まない。

**$x$ 軸上の点の基本近傍は閉であること。**
$N=(I_s\cup J_s)\setminus E$ の補集合が開であることを示す。
孤立点については明らかである。
$x$ 軸上の点 $w=(t,0)\ne z$ について、$I_t\cap I_s=\varnothing$（異なる縦線）、$J_t\cap J_s=\varnothing$（平行な異なる斜線）であり、
$I_t\cap J_s$ と $J_t\cap I_s$ はそれぞれ高々一点で、その $y$ 座標 $\abs{t-s}$ は正である。
従って $I_t\cup J_t$ から高々二点を除いた $w$ の基本近傍は $N$ と交わらない。
$a$ については、$n>s+2$ なら $N_n$ は $\{x<s+2\}$ に含まれる $I_s\cup J_s$ と交わらない。

**正則であること。**
:ref[lem-top-shrinking] の縮小条件を確かめる。
孤立点 $p$ では $\{p\}$ が開かつ閉である（$X\setminus\{p\}$ の各点の基本近傍から $p$ を除けばよい）。
$x$ 軸上の点では上で示したとおり基本近傍が開かつ閉である。
点 $a$ では、$\overline{N_n}\subset\{a\}\cup\{(x,y):x\ge n-2\}$ を示す。
右辺に属さない点は、$x<n-2$ の孤立点（$N_n$ に属さない開集合 $\{p\}$ をもつ）か、
$x$ 軸上の点 $(t,0)$（$t<n-2$）であり、後者の近傍 $I_t\cup J_t$ は $\{x<t+2<n\}$ に含まれて $N_n$ と交わらない。
よって $\overline{N_{n+3}}\subset\{a\}\cup\{x\ge n+1\}\subset N_n$ である。
従って $X$ は $T_3$ 空間である。

**完全正則でないこと。**
$F=[0,1]\times\{0\}$ と置く。
$F$ は閉である。
実際、孤立点は開集合 $\{p\}$ をもち、$x$ 軸上の点 $(t,0)\notin F$ の近傍 $I_t\cup J_t$ は $x$ 軸と $(t,0)$ だけで交わり、$a$ の近傍 $N_2$ は $F$ と交わらない。
$a\notin F$ である。
連続関数 $f:X\to[0,1]$ が $F$ 上で $0$ ならば $f(a)=0$ となることを示す。
これにより $a$ と $F$ は連続関数で分離できない。

(a) $f(s,0)=0$ なら、$B_s=\{w\in I_s\cup J_s:f(w)\ne0\}$ は高々可算である。
実際、各 $k$ について連続性から基本近傍 $(I_s\cup J_s)\setminus E_k$ で $\abs f<1/k$ となるものがあり、
$\abs f\ge1/k$ となる $I_s\cup J_s$ の点は有限集合 $E_k$ に含まれる。
$B_s=\bigcup_k\{w:\abs{f(w)}\ge1/k\}$ は可算個の有限集合の和である。

(b) $I_t$ の無限個の点で $f=0$ なら $f(t,0)=0$ である。
$f(t,0)=c\ne0$ とすると、連続性から基本近傍 $(I_t\cup J_t)\setminus E$ 上で $\abs{f-c}<\abs c/2$、従って $f\ne0$ となるが、
$f=0$ となる $I_t$ の無限個の点のうち有限集合 $E$ の外にあるものが存在して矛盾する。

(c) $Z_n=\{t\in[n-1,n]:f(t,0)=0\}$ と置く。
$Z_n$ が無限なら $[n,n+1]\setminus Z_{n+1}$ は高々可算であることを示す。
可算無限集合 $A\subset Z_n$ を取り、$B=\bigcup_{s\in A}B_s$ と置くと、(a) により $B$ は高々可算である。
縦線分 $I_t$ は $t$ ごとに互いに交わらないので、$T=\{t:I_t\cap B\ne\varnothing\}$ は高々可算である。
$t\in[n,n+1]\setminus T$ を取る。
$s\in A$ に対して $0\le t-s\le2$ であり、$t-s\in\{0,2\}$ となる $s$ は高々二つ（$s=t$ または $s=n-1,t=n+1$）である。
残りの $s\in A$ について点 $(t,t-s)$ は $J_s\cap I_t$ に属し、$s$ が異なれば $y$ 座標が異なる。
$t\notin T$ よりこれらの点は $B_s$ に属さず、$f(t,t-s)=0$ である。
従って $I_t$ の無限個の点で $f=0$ であり、(b) により $f(t,0)=0$、すなわち $t\in Z_{n+1}$ である。

(d) 仮定より $Z_1=[0,1]$ は無限である。
(c) を繰り返すと、すべての $n$ で $Z_n$ は $[n-1,n]$ から高々可算集合を除いたもの（$n\ge2$）であり、特に空でない。
$f(a)=c>0$ と仮定すると、連続性からある $N_m$ 上で $f>c/2$ である。
しかし $n>m+1$ と $t\in Z_n$ を取れば $(t,0)\in N_m$ かつ $f(t,0)=0$ であり、矛盾する。
よって $f(a)=0$ である。
:::

### 分離公理の全体像

本章で扱った条件の含意と、逆が成り立たないことを示す例をまとめる。
各行の左の性質は右隣の行の性質を導くが、その逆は表の右端の例で否定される。

| 性質 | 導く性質 | 導かれる性質をもつが元の性質をもたない例 |
|---|---|---|
| 距離化可能 | $T_4$ | Sorgenfrey 直線（:ref[exm-top-sorgenfrey-normal]） |
| $T_4$ | $T_{3\frac12}$ | Sorgenfrey 平面（:ref[exm-top-sorgenfrey-plane]） |
| $T_{3\frac12}$ | $T_3$ | Mysior の平面（:ref[exm-top-mysior]） |
| $T_3$ | $T_2$ | K 位相の直線（:ref[exm-top-k-topology]） |
| $T_2$ | $T_1$ | 無限集合の補有限位相（:ref[exm-top-separation-failures]） |
| $T_1$ | $T_0$ | Sierpiński 空間（:ref[exm-top-separation-failures]） |

$T_1$ を仮定しない場合、正規 $\Rightarrow$ 正則も成り立たない（Sierpiński 空間）。
また正則でも $T_0$ とは限らない（二点密着空間）。
ただし正則な $T_0$ 空間は自動的に Hausdorff になる（演習問題）。

遺伝性と積についても次のようにまとめられる。

| 性質 | 任意の部分空間 | 閉部分空間 | 任意個の積 |
|---|---|---|---|
| $T_0,T_1,T_2$ | 遺伝する | 遺伝する | 保たれる |
| $T_3$, $T_{3\frac12}$ | 遺伝する | 遺伝する | 保たれる |
| $T_4$ | 遺伝しない | 遺伝する | 保たれない |
| 距離化可能 | 遺伝する | 遺伝する | 可算積は保たれる |

最後の行の可算積については :ref[prop-top-product-metric] を、非可算積が距離化可能でないことは上の例を参照せよ。
正規性の振る舞いが他と異なることは、正規性が「二つの閉集合」という大域的な対象を扱う条件であることを反映している。

## 演習問題

::::exercise{level=basic}
二点密着空間は正規か。また $T_4$ か。

:::solution
閉集合は $\varnothing$ と全体だけであり、交わらない閉集合の組は少なくとも一方が空集合である。
空集合と全体を開集合として取れば分離できるので正規である。
しかし一点集合は閉でないので $T_1$ でなく、本章の約束では $T_4$ ではない。
:::
::::

::::exercise{level=basic}
Sierpiński 空間 $X=\{0,1\}$（開集合は $\varnothing,\{1\},X$）が $T_0$、$T_1$、正則、正規のそれぞれを満たすか判定せよ。

:::solution
$\{1\}$ が $1$ を含み $0$ を含まないので $T_0$ である。
$0$ を含み $1$ を含まない開集合はないので $T_1$ でない。
閉集合は $\varnothing,\{0\},X$ であり、交わらない閉集合の組は一方が空なので正規である。
点 $1$ と閉集合 $\{0\}$ について、$\{0\}$ を含む開集合は $X$ だけで $1$ を含むので分離できず、正則でない。
:::
::::

::::exercise{level=basic}
$\R$ で $A=\{0\}$、$B=\{3\}$ とし、$h(x)=d(x,A)/(d(x,A)+d(x,B))$ の $x=1$ での値と、$0\le x\le3$ での一般形を求めよ。

:::solution
$d(1,A)=1$、$d(1,B)=2$ なので $h(1)=1/(1+2)=1/3$ である。
$0\le x\le3$ では $d(x,A)=x$、$d(x,B)=3-x$ なので $h(x)=x/(x+3-x)=x/3$ である。
:::
::::

::::exercise{level=basic}
$A=\{-1,0,1\}\subset\R$ 上で値が順に $1,0,1$ の関数を、値域 $[0,1]$ を保って $\R$ に連続に延長せよ。

:::solution
$F(x)=\min\{1,\abs x\}$ と置く。
絶対値は連続で、$\min\{s,t\}=(s+t-\abs{s-t})/2$ も連続なので $F$ は連続である。
$F(-1)=F(1)=1$、$F(0)=0$ であり、値は $[0,1]$ に入る。
:::
::::

::::exercise{level=basic}
有界な場合の Tietze の拡張定理の証明で $M=27$ とする。
近似項 $h_0,h_1,h_2$ の大きさの上界と、三項を足した直後の誤差 $\abs{f-F_3}$ の上界を求めよ。

:::solution
$M_k=27(2/3)^k$ なので $M_0=27$、$M_1=18$、$M_2=12$、$M_3=8$ である。
$\abs{h_k}\le M_k/3$ より、上界は順に $9,6,4$ である。
$A$ 上で $\abs{f-F_3}=\abs{g_3}\le M_3=8$ である。
:::
::::

::::exercise{level=advanced}
正規空間の閉部分空間は正規であることを、定義から直接示せ。

:::solution
$A$ を正規空間 $X$ の閉部分空間とし、$E_1,E_2$ を $A$ の交わらない閉集合とする。
$E_k=A\cap F_k$（$F_k$ は $X$ で閉）と書けるので、$E_k$ は二つの閉集合の共通部分として $X$ でも閉である。
正規性から $E_1\subset U_1$、$E_2\subset U_2$、$U_1\cap U_2=\varnothing$ となる $X$ の開集合を取ると、
$A\cap U_1$、$A\cap U_2$ は $A$ の開集合で、$E_1,E_2$ を分離する。
:::
::::

::::exercise{level=advanced}
正則な $T_0$ 空間は Hausdorff であることを示せ。

:::solution
$x\ne y$ とする。
$T_0$ 性から、一方だけを含む開集合 $U$ がある。
必要なら $x,y$ を入れ替えて $x\in U$、$y\notin U$ とする。
$F=X\setminus U$ は $y$ を含み $x$ を含まない閉集合である。
正則性から $x\in V$、$F\subset W$、$V\cap W=\varnothing$ となる開集合が取れ、$y\in F\subset W$ なので $V,W$ が $x,y$ を分離する。
:::
::::

::::exercise{level=advanced}
$X,Y$ を空でない位相空間とする。
$X\times Y$ が Hausdorff であることと、$X$ と $Y$ がともに Hausdorff であることが同値であることを示せ。

:::solution
$X,Y$ が Hausdorff なら :ref[prop-top-hausdorff-product] により積も Hausdorff である。
逆に $X\times Y$ が Hausdorff とする。
$Y$ は空でないので $y_0\in Y$ を取ると、$X\to X\times\{y_0\}$、$x\mapsto(x,y_0)$ は連続な全単射で、逆写像は第一射影の制限なので連続である。
よって $X$ は Hausdorff 空間 $X\times Y$ の部分空間 $X\times\{y_0\}$ と同相であり、:ref[prop-top-separation-hereditary] により Hausdorff である。
$Y$ についても同様である。
:::
::::

::::exercise{level=advanced}
どの位相空間 $X$ とどの連続写像 $f:X\to Y$ に対しても、グラフが $X\times Y$ で閉となるための、空間 $Y$ に関する必要十分条件を求めよ。さらに $Y=\R$ で、グラフが閉なのに不連続な実関数を一つ構成せよ。

:::solution
恒等写像 $Y\to Y$ のグラフは対角集合である。対角集合が閉ならば、異なる点 $a,b$ に対し $(a,b)$ を含み対角集合を避ける開長方形 $U\times V$ が取れ、$U\cap V=\varnothing$ なので $Y$ は Hausdorff である。逆に $Y$ が Hausdorff なら対角集合は閉である。連続な $f$ のグラフは連続写像 $(x,y)\mapsto(f(x),y)$ による対角集合の逆像なので閉である。

最後に $f(0)=0$、$f(x)=1/x$（$x\ne0$）とする。グラフは閉集合 $\{(x,y):xy=1\}$ と一点 $\{(0,0)\}$ の和として閉であるが、$1/n\to0$ の像 $n$ は $0$ に収束しない。値域までコンパクトな場合には :ref[lem-top-tube] により、閉グラフと $X\times F$（$F\subset Y$ は閉）の共通部分を $X$ に射影して得る $f^{-1}(F)$ が閉となり、この反例は成立しない。
:::
::::

::::exercise{level=advanced}
可分な空間 $X$ から $\R$ への連続関数全体の濃度は、$\R$ の濃度以下であることを示せ。

:::solution
$Q$ を可算稠密集合とする。
二つの連続関数 $f,g$ が $Q$ 上で一致すれば、$\R$ は Hausdorff なので :ref[cor-top-dense-uniqueness] により $f=g$ である。
よって制限写像 $f\mapsto f|_Q$ は単射であり、連続関数全体の濃度は $\R^{Q}$ の濃度、すなわち $(2^{\aleph_0})^{\aleph_0}$ 以下である。
$(2^{\aleph_0})^{\aleph_0}=2^{\aleph_0\cdot\aleph_0}=2^{\aleph_0}$（$\N\times\N$ が可算であることによる、[濃度の演算](/math/sets-logic/cardinal-arithmetic/) を参照）であり、これは $\R$ の濃度に等しい。
:::
::::

::::exercise{level=advanced}
可算個の連続関数が点を区別するだけでは、それらから作った距離が元の位相を与えるとは限らない。
Sorgenfrey 直線上の関数 $f(x)=\tfrac12+\arctan(x)/\pi$ を使ってこのことを示せ。

:::solution
$f$ は通常の位相で連続な狭義単調増加関数であり、Sorgenfrey 位相は通常の位相より細かいので（:ref[prop-top-identity-topology]）、Sorgenfrey 直線上でも連続である。
単射なので点を区別する。
この一つの関数から作った距離 $\rho(x,y)=\abs{f(x)-f(y)}$ は、$f$ が $\R$ から $(0,1)$ への通常位相での同相写像なので、通常の位相を与える。
しかし Sorgenfrey 位相で開な $[0,1)$ は通常の位相で開でないので、$\rho$ は Sorgenfrey 位相を与えない。
距離化定理の証明では、各点 $x$ と開近傍 $G$ に対し、$x$ で正で $G$ の外で消える関数まで用意した（:ref[lem-top-embedding-functions] の条件）。
$f$ はこの条件を満たさない。
:::
::::

::::exercise{level=advanced}
局所コンパクト Hausdorff 空間が正則であることを、:ref[lem-top-compact-neighborhood] から直接示せ。

:::solution
$X$ を局所コンパクト Hausdorff 空間、$G$ を開集合、$x\in G$ とする。
:ref[lem-top-compact-neighborhood] により、$x\in V\subset\overline V\subset G$ で $\overline V$ がコンパクトな開集合 $V$ がある。
これは :ref[lem-top-shrinking] の縮小条件そのものなので $X$ は正則である。
Hausdorff なので $T_1$ でもあり、$X$ は $T_3$ である（:ref[prop-top-lch-completely-regular] によればさらに $T_{3\frac12}$ である）。
:::
::::

## 確認テスト

::::quiz
本章の約束で、常に成り立つ含意をすべて選べ。

- [x] $T_4\Rightarrow T_{3\frac12}$
- [x] $T_{3\frac12}\Rightarrow T_3$
- [ ] 正規 $\Rightarrow$ 正則（$T_1$ を仮定しない）
- [ ] $T_3\Rightarrow T_4$

:::explanation
(A) は Urysohn の補題を一点集合（$T_1$ 性から閉）と閉集合に適用すれば得られる。
(B) は、$f(x)=1$、$F$ 上 $f=0$ の関数に対して $\{f>1/2\}$、$\{f<1/2\}$ が $x$ と $F$ を分離することから従う。
(C) は Sierpiński 空間が反例である。正規であるが点 $1$ と閉集合 $\{0\}$ を分離できない。
(D) は Sorgenfrey 平面が反例である。$T_3$ であるが、反対角線と Jones の補題による数え上げから正規でない。
:::
::::

::::quiz
非可算離散空間の一点コンパクト化で、追加点 $\infty$ だけをちょうど零点とする連続実数値関数は存在するか。

- [ ] 正規空間なので必ず存在する。
- [x] $\{\infty\}$ が可算個の開集合の共通部分にならないので存在しない。
- [ ] $\{\infty\}$ が閉でないので存在しない。

:::explanation
一点コンパクト化はコンパクト Hausdorff なので正規であり、$\{\infty\}$ は閉である。
$\infty$ の開近傍の補集合は有限集合なので、可算個の開近傍の共通部分の補集合は高々可算であり、非可算個の元の点を除けない。
よって $\{\infty\}$ は $G_\delta$ 集合でない。
連続関数 $g$ の零点集合は $\bigcap_n\{\abs g<1/n\}$ と $G_\delta$ 集合なので、(B) が正しい。
:::
::::

::::quiz
正規性だけから保証されるものを選べ。

- [x] 交わらない二つの閉集合を、連続な $[0,1]$ 値関数で分けられる。
- [ ] 可算な開基が存在する。
- [ ] 距離化可能である。

:::explanation
(A) は Urysohn の補題そのものである。
非可算な離散空間は正規であるが可算開基をもたないので (B) は誤りである。
二点密着空間は正規であるが Hausdorff でないので距離化できず、(C) も誤りである。
:::
::::

::::quiz{answer="6" tol="1e-6"}
Tietze の拡張定理の一段階の補題で、$A$ 上の残差 $g$ の上界が $M=9$ のとき、修正後の残差 $g-h$ に保証される上界はいくらか。

:::explanation
補題により $\abs h\le M/3=3$ の連続関数 $h$ で $\abs{g-h}\le2M/3$ とできる。
$M=9$ なら $2\cdot9/3=6$ である。
:::
::::

::::quiz{answer="0.6666666667" tol="1e-6"}
$\{0,1\}^{\N}$ に距離 $d(x,y)=\sum_{n\ge1}2^{-n}\abs{x_n-y_n}$ を入れる。
奇数番目の桁だけが $1$ である点 $x=(1,0,1,0,\ldots)$ と零点 $0=(0,0,\ldots)$ の距離を小数で答えよ。

:::explanation
異なる桁は $n=1,3,5,\ldots$ なので
$d(x,0)=\sum_{k\ge1}2^{-(2k-1)}=2\sum_{k\ge1}4^{-k}=2\cdot\frac{1/4}{1-1/4}=\frac23\approx0.6666666667$ である。
:::
::::

::::quiz
Urysohn の距離化定理の証明で、可算個の座標関数で足りる理由を選べ。

- [ ] 空間の点が可算個だから。
- [x] 可算開基の元の順序対 $(i,j)$ が高々可算個だから。
- [ ] 開集合がすべて可算集合だから。

:::explanation
証明では、$\overline{B_i}\subset B_j$ を満たす開基の元の対 $(i,j)$ ごとに Urysohn 関数 $f_{ij}$ を一つ選んだ。
開基が可算なら、その対の全体は $\N\times\N$ の部分集合として高々可算である。
空間自体や各開集合の点の個数には制限がない（例えば $\R$ は非可算だが第二可算である）。
:::
::::

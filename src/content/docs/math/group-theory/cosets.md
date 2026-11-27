---
title: 剰余類と Lagrange の定理
description: 部分群による群の分割（剰余類）から Lagrange の定理と指数の公式を導き、正規部分群の判定法と剰余群の構成を与える。
notation:
  - symbol: 'gH,\ Hg'
    meaning: 部分群 H の、g を代表元とする左剰余類と右剰余類
  - symbol: 'G/H'
    meaning: H の左剰余類全体の集合（H が正規部分群なら剰余群）
  - symbol: 'H\backslash G'
    meaning: H の右剰余類全体の集合
  - symbol: '[G:H]'
    meaning: 部分群 H の G における指数
  - symbol: 'HK'
    meaning: 積集合 {hk | h∈H, k∈K}
  - symbol: 'HgK,\ H\backslash G/K'
    meaning: 両側剰余類とその全体
  - symbol: 'gHg^{-1}'
    meaning: H の g による共役部分群
  - symbol: 'N\trianglelefteq G'
    meaning: N が G の正規部分群であること
  - symbol: 'N_G(H)'
    meaning: 部分群 H の正規化群
  - symbol: 'G/N'
    meaning: 正規部分群 N による剰余群
  - symbol: '\pi_N'
    meaning: 自然な射影 G → G/N, g ↦ gN
---

## この章の目標

- 部分群 $H$ による左剰余類・右剰余類への分割を理解し、代表元を取り替える条件 $g^{-1}k\in H$ を使いこなす。
- Lagrange の定理と指数の乗法性を証明し、元の位数、Euler の定理、位数の小さい群の分類に応用する。Lagrange の定理の逆が成り立たないことを交代群 $A_4$ で確かめる。
- 正規部分群の同値な特徴付けを証明し、具体的な部分群の正規性を判定する。部分群の積、正規化群、正規核を扱う。
- 正規部分群による剰余群を構成し、剰余類の積が代表元によらずに定まることと正規性が同値であることを示す。剰余群の元の位数と構造を計算する。

前章 [群の定義と例](/math/group-theory/groups/) では、群とその部分群、巡回群と元の位数を調べた。本章では部分群 $H\le G$ を一つ固定し、$G$ を「$H$ を平行移動した集合」に分割する。整数を $3$ で割った余りで分類すると、各類は部分群 $3\Z$ を平行移動した $a+3\Z$ である。これを一般の群で行ったものが剰余類であり、分割の各部分がすべて $H$ と同じ大きさをもつことから、有限群の部分群の位数は群の位数を割る（Lagrange の定理）。

後半では、剰余類の全体そのものを群にすることを考える。整数の場合は余りどうしを足せたが、一般の群では「代表元を取り替えても積の類が変わらない」とは限らない。これが成り立つ条件がちょうど正規性であり、そのとき剰余類の全体は剰余群をなす。剰余群は次章 [準同型と同型定理](/math/group-theory/homomorphisms/) で準同型の像を記述する道具となる。

本章を通して、群の演算は乗法的に書き、単位元を $e$ とする（加法群では $0$、行列群では $I$、$Q_8$ では $1$）。置換の積 $\sigma\tau$ は合成 $\sigma\circ\tau$ であり、右の $\tau$ を先に施す（:ref[def-grp-symmetric]）。二面体群 $D_n$ の元は標準形 $r^is^j$ で書く（:ref[prop-grp-dihedral-form]）。

## 左剰余類と右剰余類

### 部分群が定める同値関係

整数 $a,b$ が「$3$ で割った余りが等しい」ことは $b-a\in3\Z$ と書ける。差が部分群に入るという条件は、乗法的に書けば $a^{-1}b\in H$ となる。ところが非可換な群では $a^{-1}b$ と $ba^{-1}$ は一般に異なるので、左右二通りの関係が生じる。

:::proposition[部分群が定める同値関係]{#prop-grp-coset-relation}
$H\le G$ とする。$G$ 上の関係を
$$
g\sim_H k\ :\Longleftrightarrow\ g^{-1}k\in H,\qquad
g\ {}_H{\sim}\ k\ :\Longleftrightarrow\ kg^{-1}\in H
$$
で定めると、いずれも同値関係である。$g$ を含む同値類は、前者では $\{gh\mid h\in H\}$、後者では $\{hg\mid h\in H\}$ である。
:::

:::proof
前者について、$g^{-1}g=e\in H$ より反射律が成り立つ。$g^{-1}k\in H$ なら、その逆元 $(g^{-1}k)^{-1}=k^{-1}g$ も $H$ に属するので対称律が成り立つ。$g^{-1}k,\ k^{-1}l\in H$ なら、積 $(g^{-1}k)(k^{-1}l)=g^{-1}l$ も $H$ に属するので推移律が成り立つ。$g$ の同値類は $\{k\mid g^{-1}k\in H\}$ であり、$g^{-1}k=h$ と $k=gh$ は同値だから、これは $\{gh\mid h\in H\}$ に等しい。

後者も同様である。$gg^{-1}=e$、$(kg^{-1})^{-1}=gk^{-1}$、$(lk^{-1})(kg^{-1})=lg^{-1}$ から三つの律が従い、$kg^{-1}=h$ と $k=hg$ は同値なので $g$ の類は $\{hg\mid h\in H\}$ である。
:::

:::definition[剰余類と指数]{#def-grp-coset}
$H\le G$、$g\in G$ とする。
$$
gH=\{gh\mid h\in H\},\qquad Hg=\{hg\mid h\in H\}
$$
をそれぞれ $H$ の :term[左剰余類]{en="left coset" yomi="ひだりじょうよるい"}、:term[右剰余類]{en="right coset" yomi="みぎじょうよるい"} といい、両者を総称して :term[剰余類]{en="coset" yomi="じょうよるい"} という。$g$ をその剰余類の**代表元**という。左剰余類全体の集合を $G/H$、右剰余類全体の集合を $H\backslash G$ と書く。$G/H$ の元の個数（濃度）を $H$ の $G$ における :term[指数]{en="index" yomi="しすう"} といい、$[G:H]$ と書く。$G/H$ が無限集合のときは $[G:H]=\infty$ とも書く。
:::

加法的に書かれた群では、左剰余類を $g+H=\{g+h\mid h\in H\}$ と書く。可換群では $gH=Hg$ なので左右の区別はない。また、文献によっては $gH$ を右剰余類と呼ぶものもあるので、他書を読むときは定義を確かめる必要がある。この段階で $G/H$ は単なる集合であり、群の構造は入っていない。群の構造が入る条件は本章の後半で調べる。

:::proposition[剰余類の基本性質]{#prop-grp-coset-test}
$H\le G$、$g,k\in G$ とする。

1. $gH=kH$ であることと $g^{-1}k\in H$ は同値であり、これはさらに $k\in gH$ と同値である。特に $gH=H$ と $g\in H$ は同値である。
2. 二つの左剰余類は一致するか交わらないかのいずれかであり、$G$ は左剰余類の互いに交わらない和集合である。
3. $Hg=Hk$ であることと $kg^{-1}\in H$ は同値であり、これはさらに $k\in Hg$ と同値である。右剰余類についても 2 と同じことが成り立つ。
4. $h\mapsto gh$ は $H$ から $gH$ への全単射であり、$h\mapsto hg$ は $H$ から $Hg$ への全単射である。特に、すべての左剰余類と右剰余類は $H$ と同じ濃度をもつ。
:::

:::proof
1. :ref[prop-grp-coset-relation] により $gH$ と $kH$ は同値関係 $\sim_H$ の同値類である。同値類が一致することは代表元が同値であること、すなわち $g^{-1}k\in H$ と同値である。また $k\in gH$ は $k=gh$（$h\in H$）と書けることであり、これも $g^{-1}k\in H$ と同値である。$k=e$ とおけば、$gH=H$ と $g^{-1}\in H$、すなわち $g\in H$ の同値性を得る。

2. 同値関係の同値類は $G$ の分割をなす（[同値関係と商集合](/math/sets-logic/relations/)）。直接示すなら、$x\in gH\cap kH$ のとき 1 により $xH=gH$ かつ $xH=kH$ だから $gH=kH$ である。また各 $g$ は $g=ge\in gH$ を満たすので、左剰余類の和集合は $G$ である。

3. 右側の同値関係 ${}_H{\sim}$ について 1、2 と同じ議論をすればよい。$k\in Hg$ は $k=hg$、すなわち $kg^{-1}\in H$ と同値である。

4. $x\mapsto g^{-1}x$ は $gH$ を $H$ に写し、$h\mapsto gh$ の逆写像である。実際 $g^{-1}(gh)=h$、$g(g^{-1}x)=x$ である。右剰余類についても $x\mapsto xg^{-1}$ が逆写像になる。
:::

:::remark
1 の条件 $g^{-1}k\in H$ を $gk^{-1}\in H$ と取り違えやすい。$gk^{-1}\in H$ は 3 により $Hk=Hg$ を意味し、右剰余類についての条件である。可換な加法群では二つの差 $k-g$ と $g-k$ が互いに逆元であり、一方が $H$ に入れば他方も入るため、整数の合同式の感覚で計算すると誤る。また $gH=kH$ でも $g=k$ とは限らない。剰余類のどの元も代表元になりうる。
:::

左剰余類と右剰余類は一般に異なる分割を与えるが、その個数は常に等しい。

:::proposition[左右の剰余類の対応]{#prop-grp-left-right-cosets}
$H\le G$ とする。部分集合 $A\subseteq G$ に対して $A^{-1}=\{a^{-1}\mid a\in A\}$ とおくと、$(gH)^{-1}=Hg^{-1}$ が成り立ち、$gH\mapsto Hg^{-1}$ は $G/H$ から $H\backslash G$ への全単射である。特に、左剰余類の個数と右剰余類の個数は等しい。
:::

:::proof
$(gh)^{-1}=h^{-1}g^{-1}$ であり、$h$ が $H$ を動くと $h^{-1}$ も $H$ 全体を動くから $(gH)^{-1}=Hg^{-1}$ である。写像 $A\mapsto A^{-1}$ は部分集合 $A$ だけで定まるので、$gH\mapsto Hg^{-1}$ は代表元の選び方によらない。同様に $Hg\mapsto(Hg)^{-1}=g^{-1}H$ は $H\backslash G$ から $G/H$ への写像であり、$(A^{-1})^{-1}=A$ だから二つの写像は互いに逆である。
:::

ここで $gH\mapsto Hg$ とすると、代表元によって行き先が変わることがある。次の例で確かめる。

### 剰余類の例

:::example[整数と実数]{#exm-grp-cosets-integers}
$n\ge1$ とし、$G=\Z$、$H=n\Z$ とする。:ref[prop-grp-coset-test] により $a+n\Z=b+n\Z$ は $b-a\in n\Z$、すなわち $a\equiv b\pmod n$ と同値である。各整数は $a=qn+r$（$0\le r<n$）と割れるので $a+n\Z=r+n\Z$ である。また $0,1,\ldots,n-1$ のうち相異なる二数の差は絶対値が $n$ 未満の非零整数で $n$ の倍数でないから、これらは互いに異なる類に属する。よって
$$
\Z/n\Z=\{\,0+n\Z,\ 1+n\Z,\ \ldots,\ (n-1)+n\Z\,\},\qquad [\Z:n\Z]=n
$$
である。この集合は前章で $\bar a$ と書いた類の集合（:ref[exm-grp-integers-mod]）にほかならない。なお $n=0$ なら $H=\{0\}$ であり、各剰余類は一点からなるので $[\Z:\{0\}]=\infty$ である。

実数の加法群 $\R$ と部分群 $\Z$ では、$x+\Z=y+\Z$ は $x-y\in\Z$ と同値である。床関数 $\lfloor x\rfloor$ を用いると $x+\Z=(x-\lfloor x\rfloor)+\Z$ であり、区間 $[0,1)$ の相異なる二数の差は整数でないから、各剰余類は $[0,1)$ にちょうど一つの代表元をもつ。従って $\R/\Z$ は集合として $[0,1)$ と一対一に対応し、指数 $[\R:\Z]$ は非可算である。
:::

:::example[三文字の置換群]{#exm-grp-cosets-s3}
$G=S_3$、$H=\langle(1\,2)\rangle=\{e,(1\,2)\}$ とする。積は右から施すので、$(1\,2\,3)(1\,2)$ は $1\mapsto2\mapsto3$、$2\mapsto1\mapsto2$、$3\mapsto3\mapsto1$ と動かし $(1\,3)$ に等しい。一方 $(1\,2)(1\,2\,3)$ は $1\mapsto2\mapsto1$、$2\mapsto3\mapsto3$、$3\mapsto1\mapsto2$ と動かし $(2\,3)$ に等しい。同様に計算すると、左剰余類と右剰余類はそれぞれ
$$
\begin{aligned}
G/H&:\ \{e,(1\,2)\},\quad \{(1\,2\,3),(1\,3)\},\quad \{(1\,3\,2),(2\,3)\},\\
H\backslash G&:\ \{e,(1\,2)\},\quad \{(1\,2\,3),(2\,3)\},\quad \{(1\,3\,2),(1\,3)\}
\end{aligned}
$$
である。$(1\,2\,3)H$ と $H(1\,2\,3)$ はともに $(1\,2\,3)$ を含む二元集合だが一致しない。「交われば一致する」のは同じ側の剰余類どうしに限る。

$gH\mapsto Hg$ が定義できないことも分かる。$(1\,2\,3)H=(1\,3)H$ であるが、$H(1\,2\,3)=\{(1\,2\,3),(2\,3)\}$ と $H(1\,3)=\{(1\,3),(1\,3\,2)\}$ は異なる。これに対し :ref[prop-grp-left-right-cosets] の対応では、$(1\,2\,3)H$ と $(1\,3)H$ はともに $H(1\,3\,2)=H(1\,3)=\{(1\,3\,2),(1\,3)\}$ に写る。

部分群を $A_3=\langle(1\,2\,3)\rangle=\{e,(1\,2\,3),(1\,3\,2)\}$ に替えると、左剰余類は $A_3$ と $(1\,2)A_3=\{(1\,2),(2\,3),(1\,3)\}$ であり、右剰余類も $A_3$ と $A_3(1\,2)=\{(1\,2),(1\,3),(2\,3)\}$ である。この場合は左右の分割が一致する。
:::

:::example[一文字を固定する部分群]{#exm-grp-cosets-symmetric}
$n\ge2$ とし、$H=\{\sigma\in S_n\mid\sigma(n)=n\}$ とおく。$\sigma,\tau\in H$ なら $\tau^{-1}(n)=n$ より $\sigma\tau^{-1}(n)=n$ なので、$H$ は部分群である（:ref[prop-grp-subgroup-test]）。$\sigma\in H$ を $\{1,\ldots,n-1\}$ に制限する写像は $H$ から $S_{n-1}$ への全単射なので $|H|=(n-1)!$ である。

$\sigma H=\tau H$ は $\sigma^{-1}\tau\in H$、すなわち $\sigma^{-1}\tau(n)=n$ と同値であり、両辺に $\sigma$ を施せば $\tau(n)=\sigma(n)$ と同値である。従って左剰余類は「$n$ の行き先」で分類され、$\sigma H$ は $n$ を $\sigma(n)$ に送る置換全体である。$n$ の行き先は $1,\ldots,n$ のどれにもなりうる（$i\ne n$ なら互換 $(i\,n)$ が $n$ を $i$ に送り、$i=n$ なら恒等置換が $n$ を固定する）ので $[S_n:H]=n$ である。

右剰余類では $H\sigma=H\tau$ が $\tau\sigma^{-1}(n)=n$、すなわち $\sigma^{-1}(n)=\tau^{-1}(n)$ と同値である。従って右剰余類は「$n$ に送られる文字」で分類される。左右で分類の基準が異なることに注意する。

各剰余類が $(n-1)!$ 個の元をもち、それらが $n$ 個あるので $|S_n|=n\cdot(n-1)!$ となる。これは $|S_n|=n!$（:ref[prop-linalg-symmetric-group]）の帰納的な証明にもなっている。
:::

:::example[二面体群の剰余類]{#exm-grp-cosets-dihedral}
$D_n$（$n\ge3$）の頂点を $\Z/n\Z$ で番号づけ、$r(i)=i+1$、$s(i)=-i$ とする（:ref[def-grp-dihedral]）。$r^is^j$ は頂点 $0$ を $i+(-1)^j\cdot0=i$ に送る。

回転全体 $R=\langle r\rangle$ について、$sR=\{sr^i\}=\{r^{-i}s\}$ と $Rs=\{r^is\}$ はどちらも鏡映全体である。従って左剰余類も右剰余類も「回転全体」と「鏡映全体」の二つであり、$[D_n:R]=2$ である。

$H=\langle s\rangle=\{e,s\}$ については事情が異なる。$r^is^j(0)=i$ だから、頂点 $0$ を固定する元は $e$ と $s$ だけであり $H$ に等しい。従って $gH=kH$ は $g^{-1}k(0)=0$、すなわち $g(0)=k(0)$ と同値である。左剰余類 $r^iH=\{r^i,\ r^is\}$ は頂点 $0$ を $i$ に送る元全体であり、左剰余類と頂点が一対一に対応する。特に $[D_n:H]=n$ である。右剰余類は $Hr^i=\{r^i,\ sr^i\}=\{r^i,\ r^{-i}s\}$ であり、$n\ge3$ なら $r^2\ne e$ より $rs\ne r^{-1}s$ だから $rH\ne Hr$ である。

この例と :ref[exm-grp-cosets-symmetric] は、「ある点を固定する部分群の左剰余類は、その点の行き先と対応する」という同じ構造をもつ。この対応は [群の作用](/math/group-theory/group-actions/) で軌道・固定群定理として一般化される。
:::

:::example[行列式による分類]{#exm-grp-cosets-linear}
体 $K$ 上の $\SL_n(K)=\{A\in\GL_n(K)\mid\det A=1\}$（:ref[exm-grp-matrix-subgroups]）を考える。行列式の積公式（:ref[thm-linalg-det-product]）より
$$
A\,\SL_n(K)=B\,\SL_n(K)\iff\det(A^{-1}B)=1\iff\det A=\det B
$$
である。右剰余類でも $\SL_n(K)A=\SL_n(K)B\iff\det(BA^{-1})=1\iff\det A=\det B$ となる。従って左剰余類も右剰余類も「行列式の値が $c$ である行列全体」（$c\in K^\times$）であり、両者は一致する。$c\in K^\times$ に対し $\diag(c,1,\ldots,1)$ の行列式は $c$ なので、どの $c$ も現れ、$[\GL_n(K):\SL_n(K)]=|K^\times|$ である。

$K=\F_p$ なら剰余類は $p-1$ 個であり、各剰余類は $\SL_n(\F_p)$ と同じ個数の元をもつ（:ref[prop-grp-coset-test] の 4）。よって :ref[prop-grp-gl-order] から
$$
|\SL_2(\F_3)|=\frac{(3^2-1)(3^2-3)}{2}=24,\qquad |\SL_2(\F_5)|=\frac{(5^2-1)(5^2-5)}{4}=120
$$
を得る。
:::

:::example[平面の中の剰余類]{#exm-grp-cosets-geometric}
乗法群 $\C^\times$ と単位円 $S^1=\{z\mid |z|=1\}$ を考える。$zS^1=wS^1$ は $|z^{-1}w|=1$、すなわち $|z|=|w|$ と同値である。従って左剰余類は原点を中心とする円であり、半径 $\rho>0$ と一対一に対応する。

加法群 $\R^2$ と、零でないベクトル $v$ の張る直線 $L=\R v$ を考える。$x+L$ は $x$ を通り $L$ に平行な直線であり、$x+L=y+L$ は $y-x\in L$ と同値である。$v$ と一次独立な $w$ を一つとると、各平行線は直線 $\R w$ とちょうど一点で交わるので、$\R w$ の点が剰余類の代表元の系をなす。部分空間による剰余類はアフィン部分空間であり、その全体は商空間（[双対空間と商空間](/math/linear-algebra/dual-and-quotient/)）として再びベクトル空間になる。後で述べる剰余群はその群版である。
:::

## Lagrange の定理

有限集合を同じ大きさの部分に分割すれば、全体の個数は「部分の個数 × 部分の大きさ」である。剰余類はまさにこの状況を与えるので、部分群の位数についての強い制約が得られる。

### Lagrange の定理と元の位数

:::theorem[Lagrange の定理]{#thm-grp-lagrange}
有限群 $G$ とその部分群 $H$ に対して
$$
|G|=[G:H]\,|H| \label{eq-grp-lagrange}
$$
が成り立つ。特に $|H|$ と $[G:H]$ はともに $|G|$ を割り、$[G:H]=|G|/|H|$ である。
:::

:::proof
:ref[prop-grp-coset-test] の 2 により、$G$ は $[G:H]$ 個の左剰余類の互いに交わらない和集合である。同じ命題の 4 により各左剰余類はちょうど $|H|$ 個の元をもつ。元の個数を剰余類ごとに数えて加えれば $|G|=[G:H]\,|H|$ を得る。
:::

:::remark{#rem-grp-lagrange-infinite}
$G$ が無限群の場合も、各左剰余類から代表元を一つずつ選んだ集合 $R$ をとると（一般には選択公理を用いる。[選択公理と Zorn の補題](/math/sets-logic/choice/) 参照）、$R\times H\to G$、$(g,h)\mapsto gh$ は全単射である。実際、$x\in G$ の属する剰余類 $xH$ の代表元を $g$ とすると $x\in gH$ なので全射であり、$gh=g'h'$ なら $gH\cap g'H\ne\emptyset$ より $gH=g'H$、代表元の一意性から $g=g'$、消去して $h=h'$ となるので単射である。従って濃度の等式 $|G|=[G:H]\cdot|H|$ が成り立つ。特に $[G:H]$ と $|H|$ がともに有限なら $G$ も有限である。ただし $n\ge1$ のとき、$\Z$ と $n\Z$ はどちらも可算無限であり、この濃度の等式だけから個々の指数 $n$ を決めることはできない。
:::

:::corollary[元の位数と群の位数]{#cor-grp-order-divides}
$G$ を有限群とする。

1. 任意の $a\in G$ の位数は有限で $|G|$ を割り、$a^{|G|}=e$ である。
2. $|G|=p$ が素数なら $G$ は巡回群であり、単位元以外のどの元も $G$ を生成する。$G$ の部分群は $\{e\}$ と $G$ だけである。
:::

:::proof
1. $a,a^2,a^3,\ldots$ は有限集合 $G$ の元なので、ある $i<j$ で $a^i=a^j$ となり、$a^{j-i}=e$ である。よって $a$ の位数 $m=\operatorname{ord}(a)$ は有限であり、:ref[thm-grp-cyclic-order] より $|\langle a\rangle|=m$ である。:ref[thm-grp-lagrange] を部分群 $\langle a\rangle$ に適用すれば $m\mid|G|$ を得る。$|G|=mq$ と書けば $a^{|G|}=(a^m)^q=e$ である（:ref[prop-grp-power-laws]）。

2. 部分群の位数は $p$ の約数 $1$ か $p$ であり、それぞれ $\{e\}$ と $G$ に対応する。$a\ne e$ なら $\langle a\rangle\ne\{e\}$ だから $\langle a\rangle=G$ である。
:::

Lagrange の定理を乗法群 $(\Z/n\Z)^\times$ に適用すると、初等整数論の基本定理が得られる。$n\ge1$ に対し、$1$ 以上 $n$ 以下で $n$ と互いに素な整数の個数を $\varphi(n)$ と書く（Euler の関数）。

:::corollary[Euler の定理と Fermat の小定理]{#cor-grp-euler-fermat}
1. （Euler）$n\ge1$ と $\gcd(a,n)=1$ を満たす整数 $a$ に対し $a^{\varphi(n)}\equiv1\pmod n$ である。
2. （Fermat）$p$ が素数で $p\nmid a$ なら $a^{p-1}\equiv1\pmod p$ である。また任意の整数 $a$ に対し $a^p\equiv a\pmod p$ である。
:::

:::proof
1. :ref[exm-grp-integers-mod] により、$n$ と互いに素な整数の類全体 $(\Z/n\Z)^\times$ は乗法について位数 $\varphi(n)$ の群をなし、$\bar a$ はその元である。:ref[cor-grp-order-divides] より $\bar a^{\varphi(n)}=\bar1$、すなわち $a^{\varphi(n)}\equiv1\pmod n$ である。

2. $p$ が素数なら $1,\ldots,p-1$ はすべて $p$ と互いに素なので $\varphi(p)=p-1$ であり、前半は 1 の特別な場合である。$p\nmid a$ なら前半の両辺に $a$ を掛けて $a^p\equiv a$ を得る。$p\mid a$ なら両辺とも $p$ の倍数である。
:::

:::example[冪の余りの計算]{#exm-grp-euler-computation}
法 $15$ では $\varphi(15)=8$ なので、$15$ と互いに素な $a$ について $a^8\equiv1$ である。実際には $2^4=16\equiv1$ であり、$2$ の位数は $4$ で $8$ を割る。Lagrange の定理から分かるのは位数が $\varphi(n)$ を割ることであり、位数そのものは個別に計算する必要がある。

$3^{100}$ を $7$ で割った余りを求める。Fermat の小定理より $3^6\equiv1\pmod7$ であり、$100=6\cdot16+4$ だから $3^{100}\equiv3^4=81\equiv4\pmod 7$ である。同様に $\varphi(10)=4$ より $7^{4}\equiv1\pmod{10}$ であり、$7^{222}=(7^4)^{55}\cdot7^2\equiv49\equiv9\pmod{10}$ を得る。
:::

:::remark
Fermat の小定理の逆は成り立たない。$341=11\cdot31$ は合成数であるが $2^{340}\equiv1\pmod{341}$ が成り立つ。実際、$2^{10}=1024=3\cdot341+1$ より $2^{10}\equiv1\pmod{341}$ であり、$2^{340}=(2^{10})^{34}\equiv1$ である。「$a^{n-1}\equiv1\pmod n$ なら $n$ は素数」という判定は誤りを含みうる。
:::

位数が偶数の群には位数 $2$ の元がある。これは Lagrange の定理の逆の最も簡単な場合であり、逆元をとる操作による組分けから示される。

:::proposition[位数 2 の元の存在]{#prop-grp-order-two-element}
$|G|$ が偶数の有限群 $G$ には位数 $2$ の元が存在する。さらに、位数 $2$ の元の個数は奇数である。
:::

:::proof
$X=\{x\in G\mid x\ne x^{-1}\}$ とおく。$x\in X$ なら $x^{-1}\in X$ であり、$x\mapsto x^{-1}$ は $X$ 上の写像で、二回施すと恒等写像になり、固定点をもたない。従って $X$ は二元集合 $\{x,x^{-1}\}$ たちに分割され、$|X|$ は偶数である。補集合 $G\setminus X=\{x\mid x^2=e\}$ は単位元と位数 $2$ の元からなる。$|G\setminus X|=|G|-|X|$ は偶数なので、位数 $2$ の元の個数 $|G\setminus X|-1$ は奇数であり、特に $1$ 以上である。
:::

例えば位数 $2$ の元は、$S_3$ に三つ（互換）、$D_4$ に五つ（四つの鏡映と $r^2$）、$Q_8$ に一つ（$-1$）、$C_6$ に一つ（$\bar3$）あり、いずれも奇数である。

### 位数の小さい群

Lagrange の定理は、位数の小さい群の構造を決めるのに役立つ。素数位数の群は :ref[cor-grp-order-divides] により巡回群である。位数 $4$ の群は次の二種類に限る。

:::proposition[位数 4 の群]{#prop-grp-order-four}
位数 $4$ の群は巡回群 $C_4$ または Klein の四元群 $V_4$ のいずれかと同型であり、$C_4$ と $V_4$ は同型でない。特に位数 $4$ の群は可換である。
:::

:::proof
$|G|=4$ とする。Lagrange の定理より各元の位数は $1,2,4$ のいずれかである。

位数 $4$ の元 $a$ があれば $\langle a\rangle$ は $4$ 元からなり $G$ に等しい。よって $G$ は位数 $4$ の巡回群であり、:ref[thm-grp-cyclic-classification] により $C_4$ と同型である。

位数 $4$ の元がなければ、単位元以外の三元 $a,b,c$ はすべて位数 $2$ である。積 $ab$ を考えると、$ab=e$ なら $b=a^{-1}=a$、$ab=a$ なら $b=e$、$ab=b$ なら $a=e$ となりいずれも矛盾するので $ab=c$ である。$a,b,c$ の役割は対称なので、相異なる二つの非単位元の積は（順序によらず）残りの一つである。さらに $x^2=e$、$ex=xe=x$ なので、積表は完全に決まる：
$$
x^2=e,\qquad ab=ba=c,\quad bc=cb=a,\quad ca=ac=b .
$$
$V_4=\{e,(1\,2)(3\,4),(1\,3)(2\,4),(1\,4)(2\,3)\}$（:ref[exm-grp-klein-four]）は位数 $4$ の群で、非単位元はいずれも互いに素な二つの互換の積なので位数 $2$ である（:ref[prop-grp-permutation-order]）。従って $V_4$ にも同じ議論が適用でき、その積表も上と同じ形になる。よって
$$
e\mapsto e,\quad a\mapsto(1\,2)(3\,4),\quad b\mapsto(1\,3)(2\,4),\quad c\mapsto(1\,4)(2\,3)
$$
は積を保つ全単射、すなわち同型（:ref[def-grp-isomorphism]）である。

最後に、同型写像は元の位数を保つ（:ref[prop-grp-isomorphism-basic]）。$C_4$ は位数 $4$ の元をもち $V_4$ はもたないので、両者は同型でない。どちらも可換だから、位数 $4$ の群は可換である。
:::

:::example[位数 4 の群の判別]{#exm-grp-order-four-examples}
- $(\Z/8\Z)^\times=\{\bar1,\bar3,\bar5,\bar7\}$ では $3^2=9$、$5^2=25$、$7^2=49$ がいずれも $1$ と合同なので、非単位元はすべて位数 $2$ であり $V_4$ と同型である。$(\Z/12\Z)^\times=\{\bar1,\bar5,\bar7,\bar{11}\}$ も同様である。
- $(\Z/5\Z)^\times$ では $2,2^2,2^3,2^4$ の類が $2,4,3,1$ となるので $\bar2$ の位数は $4$ であり、$C_4$ と同型である。
- $Q_8$ の部分群 $\langle i\rangle=\{\pm1,\pm i\}$ は $C_4$ と、$D_4$ の部分群 $\{e,r^2,s,r^2s\}$ は $V_4$ と同型である。後者では $r^2$ も $s$ も $r^2s$ も位数 $2$ だからである。
- 成分ごとの演算による $C_2\times C_2$（:ref[exm-grp-direct-product-basic]）は $V_4$ と同型である。
:::

素数位数の群と位数 $4$ の群が可換であるから、**位数 $5$ 以下の群はすべて可換である**。位数 $6$ の群は $C_6$ または $S_3$ に限ることが演習問題で示され、$S_3$ が最小の非可換群であることが分かる。

### Lagrange の定理の逆は成り立たない

有限群 $G$ と $|G|$ の約数 $d$ が与えられたとき、位数 $d$ の部分群は存在するだろうか。巡回群では各約数に対してちょうど一つ存在した（:ref[thm-grp-cyclic-subgroups]）。しかし一般には成り立たない。その証明には次の簡単な補題を用いる。

:::lemma[指数 2 の部分群は平方を含む]{#lem-grp-index-two-squares}
$[G:H]=2$ なら、すべての $g\in G$ に対して $g^2\in H$ である。
:::

:::proof
$g\in H$ なら $g^2\in H$ である。$g\notin H$ とすると、左剰余類は $H$ と $gH$ の二つであり、$gH=G\setminus H$ である。もし $g^2\notin H$ なら $g^2\in gH$ であり、$g^2=gh$（$h\in H$）と書ける。左から $g^{-1}$ を掛けると $g=h\in H$ となり矛盾する。よって $g^2\in H$ である。
:::

:::example[位数 6 の部分群をもたない A₄]{#exm-grp-a4-no-six}
$4$ 文字の偶置換全体 $A_4$（:ref[exm-grp-alternating-subgroup]）は位数 $4!/2=12$ の群である。その元を数えると、単位元が $1$ 個、互いに素な二つの互換の積 $(1\,2)(3\,4),(1\,3)(2\,4),(1\,4)(2\,3)$ が $3$ 個、$3$ 巡回置換が $8$ 個（固定する文字の選び方 $4$ 通りに、残り三文字の巡回の向き $2$ 通り）である。これらはいずれも偶置換であり（:ref[cor-linalg-cycle-sign]）、合計 $1+3+8=12$ 個なので $A_4$ の元はこれで尽きる。

$|A_4|=12$ の約数 $6$ に対し、位数 $6$ の部分群 $H\le A_4$ は存在しない。実際、存在したとすると $[A_4:H]=12/6=2$ である。$3$ 巡回置換 $c$ は位数 $3$ なので $c=c^4=(c^2)^2$ と書け、:ref[lem-grp-index-two-squares] により $c\in H$ である。すると $H$ は $8$ 個の $3$ 巡回置換と単位元を含み、$|H|\ge9>6$ となって矛盾する。

一方、他の約数に対しては部分群が存在する。位数 $1,12$ は自明であり、位数 $2$ は $\langle(1\,2)(3\,4)\rangle$ など、位数 $3$ は $\langle(1\,2\,3)\rangle$ など、位数 $4$ は $V_4$ である。欠けている約数は $6$ だけである。
:::

:::remark
Lagrange の定理の部分的な逆として、素数 $p$ が $|G|$ を割れば位数 $p$ の元が存在すること（Cauchy の定理）、さらに $|G|$ を割る最大の $p$ 冪を位数とする部分群が存在すること（Sylow の定理）が [Sylow の定理](/math/group-theory/sylow/) で示される。上の例は、約数が素数冪でないときにはこのような存在定理が期待できないことを示している。
:::

### 指数の乗法性

$K\le H\le G$ のとき、$G$ を $H$ の剰余類に分け、さらに各剰余類を $K$ の剰余類に分けることができる。これを数えると指数が掛け算になる。$G$ が無限群でも指数が有限になることがあるので、Lagrange の定理を経由しない証明を与える。

:::proposition[指数の乗法性]{#prop-grp-index-tower}
$K\le H\le G$ とする。$G/H$ の代表元の系を $\{g_i\}_{i\in I}$、$H/K$ の代表元の系を $\{h_j\}_{j\in J}$ とすると、$(i,j)\mapsto g_ih_jK$ は $I\times J$ から $G/K$ への全単射である。従って
$$
[G:K]=[G:H]\,[H:K]
$$
が濃度として成り立つ。特に $[G:K]$ が有限であることと、$[G:H]$ と $[H:K]$ がともに有限であることは同値であり、そのとき上の等式は整数の等式である。
:::

:::proof
全射性：$g\in G$ とする。$g$ はある $g_iH$ に属するので $g=g_ih$（$h\in H$）と書ける。$h$ はある $h_jK$ に属するので $h=h_jk$（$k\in K$）と書ける。よって $g=g_ih_jk\in g_ih_jK$、すなわち $gK=g_ih_jK$ である。

単射性：$g_ih_jK=g_{i'}h_{j'}K$ とする。$h_j,h_{j'}\in H$、$K\subseteq H$ なので、左辺は $g_iH$ に、右辺は $g_{i'}H$ に含まれる。二つの左 $H$ 剰余類が交わるので $g_iH=g_{i'}H$ であり、代表元の系の性質から $i=i'$ である。等式の両辺に左から $g_i^{-1}$ を掛けると $h_jK=h_{j'}K$ となり、同様に $j=j'$ である。

$I,J$ が空でないことに注意すると、$I\times J$ が有限であることと $I,J$ がともに有限であることは同値である。
:::

:::example[指数の乗法性の例]{#exm-grp-index-tower}
$12\Z\le3\Z\le\Z$ では $[\Z:3\Z]=3$、$[3\Z:12\Z]=4$ である。代表元 $g_i=i$（$0\le i<3$）と $h_j=3j$（$0\le j<4$）から得る $12$ 個の和 $i+3j$ は $0,1,\ldots,11$ を一度ずつ尽くし、$[\Z:12\Z]=12$ の代表元の系になる。$\Z$ は無限群であり、Lagrange の定理を直接使うことはできない。

$V_4\le A_4\le S_4$ では $[S_4:A_4]=2$、$[A_4:V_4]=3$ である。代表元として $S_4/A_4$ に $e,(1\,2)$、$A_4/V_4$ に $e,(1\,2\,3),(1\,3\,2)$ をとる（$(1\,2\,3)$、$(1\,3\,2)$ と $(1\,2\,3)^{-1}(1\,3\,2)=(1\,3\,2)^2=(1\,2\,3)$ はいずれも $V_4$ に属さないので、三つは相異なる左 $V_4$ 剰余類を代表する。同様に $(1\,2)$ は奇置換なので $e$ と $(1\,2)$ は相異なる左 $A_4$ 剰余類を代表する）。積をとると
$$
e,\ (1\,2\,3),\ (1\,3\,2),\ (1\,2),\ (1\,2)(1\,2\,3)=(2\,3),\ (1\,2)(1\,3\,2)=(1\,3)
$$
の $6$ 個が $S_4/V_4$ の代表元の系になる。これはちょうど文字 $4$ を固定する置換全体（$S_3$ と同一視できる部分群）である。
:::

### 二つの部分群の積集合と共通部分

二つの部分群 $H,K\le G$ に対し、積集合を
$$
HK=\{hk\mid h\in H,\ k\in K\}
$$
と定める。$HK$ は一般には部分群ではない（部分群となる条件は次節で扱う）が、その大きさは常に計算できる。

:::proposition[積集合の大きさ]{#prop-grp-product-size}
$H,K$ を群 $G$ の有限部分群とすると
$$
|HK|=\frac{|H|\,|K|}{|H\cap K|}
$$
である。
:::

:::proof
写像 $\mu:H\times K\to HK$、$(h,k)\mapsto hk$ は全射である。$hk\in HK$ の逆像を求める。$(h',k')\in H\times K$ が $h'k'=hk$ を満たすとき、$t=h^{-1}h'=kk'^{-1}$ は $H$ と $K$ の両方に属し、$h'=ht$、$k'=t^{-1}k$ である。逆に任意の $t\in H\cap K$ に対して $(ht,t^{-1}k)\in H\times K$ であり、$(ht)(t^{-1}k)=hk$ である。$t$ が異なれば $ht$ も異なるので、逆像 $\mu^{-1}(hk)$ はちょうど $|H\cap K|$ 個の元からなる。$H\times K$ を逆像ごとに数えれば $|H||K|=|HK|\,|H\cap K|$ を得る。
:::

別の見方をすると、$HK=\bigcup_{h\in H}hK$ は左 $K$ 剰余類の和集合であり、$hK=h'K$ は $h^{-1}h'\in H\cap K$ と同値だから、現れる左 $K$ 剰余類は $[H:H\cap K]$ 個である。各剰余類が $|K|$ 個の元をもつので、Lagrange の定理を $H\cap K\le H$ に適用しても同じ式を得る。

:::example[積集合の計算]{#exm-grp-product-size}
- $S_3$ で $H=\langle(1\,2)\rangle$、$K=\langle(2\,3)\rangle$ とすると $H\cap K=\{e\}$ なので $|HK|=4$ である。実際 $HK=\{e,(1\,2),(2\,3),(1\,2)(2\,3)\}$ であり、$(1\,2)(2\,3)=(1\,2\,3)$ である。$4$ は $6$ を割らないので、$HK$ は部分群ではない。
- $D_n$ で $\langle s\rangle\cap\langle r\rangle=\{e\}$ だから $|\langle s\rangle\langle r\rangle|=2\cdot n=|D_n|$ であり、$D_n=\{s^jr^i\}$ となる。標準形 $r^is^j$ と順序を逆にした表示も可能である。
- $S_4$ で文字 $4$ を固定する部分群 $S_3$ と $V_4$ を考える。$V_4$ の非単位元はどれも $4$ を動かすので $S_3\cap V_4=\{e\}$ であり、$|S_3V_4|=6\cdot4=24$ より $S_4=S_3V_4$ である。
:::

剰余類の個数についても、共通部分をとったときの評価が得られる。

:::proposition[共通部分の指数（Poincaré）]{#prop-grp-index-intersection}
$H,K\le G$ とする。

1. $h(H\cap K)\mapsto hK$ は $H/(H\cap K)$ から $G/K$ への単射であり、その像は $HK$ に含まれる左 $K$ 剰余類全体である。従って $[H:H\cap K]\le[G:K]$ であり、$[G:K]$ が有限なら等号は $HK=G$ と同値である。
2. $[G:H\cap K]\le[G:H]\,[G:K]$ である。特に $[G:H]$、$[G:K]$ が有限なら $[G:H\cap K]$ も有限であり、等号は $HK=G$ と同値である。
3. $[G:H]$ と $[G:K]$ が有限で互いに素なら、$[G:H\cap K]=[G:H][G:K]$ かつ $G=HK$ である。
:::

:::proof
1. $h,h'\in H$ について、$hK=h'K$ は $h^{-1}h'\in K$ と同値であり、$h^{-1}h'$ は最初から $H$ に属するので、これは $h^{-1}h'\in H\cap K$、すなわち $h(H\cap K)=h'(H\cap K)$ と同値である。定義域で二つの類が等しければ像も等しいので、写像は代表元によらずに定まる。逆に像が等しければ定義域の類も等しいので単射である。像は $\{hK\mid h\in H\}$ であり、その和集合が $HK$ である。$[G:K]$ が有限なら、等号は「すべての左 $K$ 剰余類が $hK$（$h\in H$）の形である」こと、すなわち $G$ の各元がある $hK$ に属すること、つまり $G=HK$ と同値である。

2. $H\cap K\le H\le G$ に :ref[prop-grp-index-tower] を適用し、1 を用いると
$$
[G:H\cap K]=[G:H]\,[H:H\cap K]\le[G:H]\,[G:K]
$$
である。両辺が有限のとき、$[G:H]\ge1$ で割れば等号は $[H:H\cap K]=[G:K]$ と同値であり、1 により $HK=G$ と同値である。

3. $H\cap K\le H$ と $H\cap K\le K$ に指数の乗法性を用いると、$[G:H]$ と $[G:K]$ はともに $[G:H\cap K]$ を割る。互いに素なので積 $[G:H][G:K]$ も $[G:H\cap K]$ を割り、2 の不等式と合わせて等号が成り立つ。2 により $G=HK$ である。
:::

この命題は、写像 $g(H\cap K)\mapsto(gH,gK)$ で見ると分かりやすい。$gH=g'H$ かつ $gK=g'K$ は $g^{-1}g'\in H\cap K$ と同値なので、この写像は単射である。左剰余類の組 $(aH,bK)$ が像に入る、すなわち $aH\cap bK\ne\emptyset$ となる条件は $a^{-1}b\in HK$ である。実際、$ah=bk$ なら $a^{-1}b=hk^{-1}$ であり、逆に $a^{-1}b=hk'$ と書ければ $g=ah=bk'^{-1}$ が両方に属する。従ってすべての組が実現されるのは $HK=G$ のときに限る。

:::example[合同式の連立]{#exm-grp-index-intersection}
$\Z$ で $H=4\Z$、$K=6\Z$ とすると $H\cap K=12\Z$ であり、$[\Z:12\Z]=12<4\cdot6$ である。実際 $H+K=2\Z\ne\Z$ であり、$a+4\Z$ と $b+6\Z$ が交わるのは $a-b\in2\Z$ のときに限る。例えば法 $4$ で $1$、法 $6$ で $3$ となる整数は $9+12\Z$ であるが、法 $4$ で $1$、法 $6$ で $2$ となる整数は存在しない。

$H=4\Z$、$K=9\Z$ なら指数 $4$ と $9$ は互いに素なので、3 により $[\Z:36\Z]=36$ かつ $4\Z+9\Z=\Z$ であり、任意の余りの組が実現される。これは中国剰余定理（[直積・半直積と有限 Abel 群](/math/group-theory/products-abelian/)）の群論的な形である。

$S_3$ で $H=\langle(1\,2)\rangle$（指数 $3$）と $K=A_3$（指数 $2$）をとると、3 から $H\cap K=\{e\}$ と $S_3=HA_3$ が従う。
:::

### 両側剰余類の大きさ

左右の両側から代表元を取り替えることを許すと、剰余類をさらに大きな組にまとめられる。二つの部分群は異なってもよい。

:::definition[両側剰余類]{#def-grp-double-coset}
$H,K\le G$、$g\in G$ に対し $HgK=\{hgk\mid h\in H,\ k\in K\}$ を :term[両側剰余類]{en="double coset" yomi="りょうがわじょうよるい"}（二重剰余類）という。両側剰余類全体の集合を $H\backslash G/K$ と書く。
:::

$g'=hgk$（$h\in H$、$k\in K$）と書けるとき $g\approx g'$ と定めると、これは同値関係である。実際、$g=ege$ より反射律、$g=h^{-1}g'k^{-1}$ より対称律、$g''=h'g'k'=(h'h)g(kk')$ より推移律が成り立つ。$g$ の同値類が $HgK$ なので、両側剰余類も $G$ の分割を与える。各両側剰余類 $HgK=\bigcup_{h\in H}hgK$ は左 $K$ 剰余類の和集合であり、同時に右 $H$ 剰余類の和集合でもある。片側の剰余類と異なり、両側剰余類の大きさは一定とは限らない。

:::proposition[両側剰余類の分解]{#prop-grp-double-size}
$H,K\le G$、$g\in G$ とする。$gKg^{-1}=\{gkg^{-1}\mid k\in K\}$ は部分群であり、$HgK$ に含まれる左 $K$ 剰余類の個数は $[H:H\cap gKg^{-1}]$ である。特に $H,K$ が有限なら
$$
|HgK|=\frac{|H|\,|K|}{|H\cap gKg^{-1}|} \label{eq-grp-double-size}
$$
である。
:::

:::proof
$gkg^{-1}(gk'g^{-1})^{-1}=g(kk'^{-1})g^{-1}$ なので :ref[prop-grp-subgroup-test] により $gKg^{-1}$ は部分群であり、$H\cap gKg^{-1}$ も部分群である。写像 $h\mapsto hgK$ は $H$ から「$HgK$ に含まれる左 $K$ 剰余類全体」への全射である。$h_1,h_2\in H$ について
$$
h_1gK=h_2gK\iff g^{-1}h_1^{-1}h_2g\in K\iff h_1^{-1}h_2\in gKg^{-1}\iff h_1^{-1}h_2\in H\cap gKg^{-1}
$$
であり（最後の同値は $h_1^{-1}h_2\in H$ による）、これは $h_1(H\cap gKg^{-1})=h_2(H\cap gKg^{-1})$ と同値である。従って $h(H\cap gKg^{-1})\mapsto hgK$ は $H/(H\cap gKg^{-1})$ から $HgK$ に含まれる左 $K$ 剰余類全体への全単射である。有限の場合、各左 $K$ 剰余類は $|K|$ 個の元をもち互いに交わらないので、Lagrange の定理と合わせて $|HgK|=[H:H\cap gKg^{-1}]\,|K|=|H||K|/|H\cap gKg^{-1}|$ を得る。
:::

$g=e$ とすれば :ref[prop-grp-product-size] の式 $|HK|=|H||K|/|H\cap K|$ が再び得られる。

:::example[三文字の両側剰余類]{#exm-grp-double-coset-s3}
$S_3$ で $H=K=\langle(1\,2)\rangle$ とする。$HeH=H$ は $2$ 元である。$g=(1\,3)$ とすると $gHg^{-1}=\{e,(1\,3)(1\,2)(1\,3)\}=\langle(2\,3)\rangle$（:ref[prop-grp-cycle-conjugation]）であり、$H$ と自明に交わる。従って $|HgH|=2\cdot2/1=4$ であり、$HgH$ は $H$ 以外の四元 $(1\,3),(2\,3),(1\,2\,3),(1\,3\,2)$ をすべて含む。$S_3$ は大きさ $2$ と $4$ の二つの両側剰余類に分かれる。大きさが一定でないので、Lagrange の定理のような「個数 × 大きさ」の公式は両側剰余類には成り立たない（$2\cdot2\ne6$）。
:::

:::example[GL₂ の上三角群による両側剰余類]{#exm-grp-bruhat-gl2}
体 $K$ 上で、上三角正則行列全体 $B=\left\{\begin{pmatrix}a&b\\0&d\end{pmatrix}\middle|\ ad\ne0\right\}$ は $\GL_2(K)$ の部分群であり（:ref[exm-grp-matrix-subgroups]）、$w=\begin{pmatrix}0&1\\1&0\end{pmatrix}$ とおく。両側剰余類は $B$ と $BwB$ の二つだけであることを示す。

$A=\begin{pmatrix}a&b\\c&d\end{pmatrix}\in\GL_2(K)$ について、$c=0$ なら $A\in B$ である。$c\ne0$ なら、$\Delta=ad-bc\ne0$ として
$$
\begin{pmatrix}1&a/c\\0&1\end{pmatrix}
\begin{pmatrix}0&1\\1&0\end{pmatrix}
\begin{pmatrix}c&d\\0&-\Delta/c\end{pmatrix}
=\begin{pmatrix}1&a/c\\0&1\end{pmatrix}
\begin{pmatrix}0&-\Delta/c\\c&d\end{pmatrix}
=\begin{pmatrix}a&b\\c&d\end{pmatrix}
$$
である（右上成分は $-\Delta/c+ad/c=b$）。左右の行列は $B$ に属するので $A\in BwB$ である。$w\notin B$ なので $BwB\ne B$ であり、両側剰余類は交わらないから、$\GL_2(K)=B\sqcup BwB$ である。この分解は Bruhat 分解の最も簡単な場合である。

$K=\F_p$ のとき大きさを確かめる。$|B|=(p-1)^2p$ である。$wBw^{-1}$ は $w\begin{pmatrix}a&b\\0&d\end{pmatrix}w=\begin{pmatrix}d&0\\b&a\end{pmatrix}$ より下三角正則行列全体なので、$B\cap wBw^{-1}$ は正則対角行列全体で位数 $(p-1)^2$ である。:eqref[eq-grp-double-size] より $|BwB|=|B|^2/(p-1)^2=p^2(p-1)^2$ であり、
$$
|B|+|BwB|=p(p-1)^2(1+p)=(p^2-1)(p^2-p)=|\GL_2(\F_p)|
$$
となって :ref[prop-grp-gl-order] と一致する。$p=3$ では $12+36=48$ である。
:::

## 正規部分群

前節までの例では、$A_3\le S_3$、$\SL_n(K)\le\GL_n(K)$、回転全体 $R\le D_n$ のように左剰余類と右剰余類の分割が一致する場合と、$\langle(1\,2)\rangle\le S_3$ や $\langle s\rangle\le D_n$ のように一致しない場合があった。$gH=Hg$ は $gHg^{-1}=H$ と書き直せる。すなわち「$H$ の元を $g$ で挟んでも $H$ から出ない」という条件である。次節で見るように、この条件はちょうど剰余類どうしの積が定義できる条件になる。

### 正規部分群の定義と判定法

:::definition[共役部分群と正規部分群]{#def-grp-normal}
$H\le G$、$g\in G$ に対して $gHg^{-1}=\{ghg^{-1}\mid h\in H\}$ とおく。これは $G$ の部分群であり、$H$ の $g$ による :term[共役部分群]{en="conjugate subgroup" yomi="きょうやくぶぶんぐん"} という。元 $ghg^{-1}$ を $h$ の $g$ による**共役**という。すべての $g\in G$ に対して $gHg^{-1}=H$ が成り立つとき、$H$ を $G$ の :term[正規部分群]{en="normal subgroup" yomi="せいきぶぶんぐん"} といい、$H\trianglelefteq G$ と書く。
:::

$gHg^{-1}$ が部分群であることは :ref[prop-grp-double-size] の証明で確かめた。写像 $h\mapsto ghg^{-1}$ は $x\mapsto g^{-1}xg$ を逆写像にもつ全単射なので、$|gHg^{-1}|=|H|$ である。また共役は積を保つ：
$$
(gxg^{-1})(gyg^{-1})=g(xy)g^{-1},\qquad (gxg^{-1})^{-1}=gx^{-1}g^{-1}. \label{eq-grp-conjugation-product}
$$
$\{e\}$ と $G$ 自身は常に正規部分群である。

:::proposition[正規性の判定法]{#prop-grp-normal-criteria}
$H\le G$ について、次は同値である。

1. $H\trianglelefteq G$。
2. すべての $g\in G$、$h\in H$ について $ghg^{-1}\in H$。すなわち、すべての $g$ について $gHg^{-1}\subseteq H$。
3. すべての $g\in G$ について $gH=Hg$。
4. $H$ の左剰余類はどれも $H$ の右剰余類である。
5. $G=\langle S\rangle$、$H=\langle T\rangle$ とするとき、すべての $s\in S$、$t\in T$ について $sts^{-1}\in H$ かつ $s^{-1}ts\in H$。

さらに $G$ が有限なら、5 の条件のうち $sts^{-1}\in H$ だけを確かめれば十分である。
:::

:::proof
$1\Rightarrow2$ は明らかである。

$2\Rightarrow1$：2 を $g^{-1}$ に適用すると $g^{-1}Hg\subseteq H$ である。この包含の両辺を $g$ で共役すると $H\subseteq gHg^{-1}$ を得るので、2 と合わせて $gHg^{-1}=H$ である。

$1\Rightarrow3$：$gH=(gHg^{-1})g=Hg$ である。$3\Rightarrow1$：$gHg^{-1}=(Hg)g^{-1}=H$ である。

$3\Rightarrow4$ は明らかである。$4\Rightarrow3$：$g\in G$ に対し、仮定により $gH=Hk$ となる $k$ がある。$g=ge\in gH=Hk$ なので、:ref[prop-grp-coset-test] の 3 により $Hk=Hg$ であり、$gH=Hg$ を得る。

$1\Rightarrow5$ は明らかである（$S,T$ の元は $G,H$ の元である）。$5\Rightarrow1$：
$$
M=\{g\in G\mid gHg^{-1}=H\}
$$
とおく。$e\in M$ であり、$g,g'\in M$ なら $(gg')H(gg')^{-1}=g(g'Hg'^{-1})g^{-1}=gHg^{-1}=H$、また $g^{-1}Hg=g^{-1}(gHg^{-1})g=H$ なので、$M$ は部分群である。次に $s\in S$ を固定する。:eqref[eq-grp-conjugation-product] より、$H$ の元を $T$ の元とその逆元の積として書けば（:ref[prop-grp-generated-words]）、その $s$ による共役は $sts^{-1}$ とその逆元の積になる。仮定によりこれは $H$ に属するので $sHs^{-1}\subseteq H$ であり、同様に $s^{-1}Hs\subseteq H$ である。後者を $s$ で共役すると $H\subseteq sHs^{-1}$ なので $sHs^{-1}=H$、すなわち $s\in M$ である。$M$ は $S$ を含む部分群なので $M\supseteq\langle S\rangle=G$ であり、$H\trianglelefteq G$ である。

$G$ が有限の場合、$sts^{-1}\in H$（$t\in T$）から上と同様に $sHs^{-1}\subseteq H$ が従い、$|sHs^{-1}|=|H|$ なので $sHs^{-1}=H$ となる。従って $s^{-1}$ についての条件は不要である。
:::

:::remark{#rem-grp-normal-pitfalls}
3 の $gH=Hg$ は集合としての等式であり、各 $h\in H$ について $gh=hg$ を意味するのではない。例えば $(1\,2)A_3=A_3(1\,2)$ であるが、$(1\,2)(1\,2\,3)=(2\,3)$ と $(1\,2\,3)(1\,2)=(1\,3)$ は異なる。

また、無限群では一つの $g$ について $gHg^{-1}\subseteq H$ が成り立っても $gHg^{-1}=H$ とは限らない。$G=\GL_2(\Q)$、$H=\left\{\begin{pmatrix}1&m\\0&1\end{pmatrix}\middle|\ m\in\Z\right\}$、$g=\diag(2,1)$ とすると
$$
g\begin{pmatrix}1&m\\0&1\end{pmatrix}g^{-1}=\begin{pmatrix}1&2m\\0&1\end{pmatrix}
$$
なので $gHg^{-1}$ は $m$ が偶数の行列全体であり、$H$ の真部分集合である。実際、$u=\begin{pmatrix}1&1\\0&1\end{pmatrix}$ とし、$G$ を $\langle g,u\rangle$ に制限すると $H=\langle u\rangle$ である。生成元 $g,u$ はそれぞれ $u$ を $u^2,u$ に共役するが、$g^{-1}ug$ の右上成分は $1/2$ なので $H$ に入らない。これで判定法 5 の逆向きの共役も必要であることが分かる。
:::

### 正規部分群の例

:::example[正規部分群の例]{#exm-grp-normal-examples}
1. **可換群**：$G$ が可換なら $ghg^{-1}=h$ なので、すべての部分群が正規である。
2. **中心**：中心 $Z(G)=\{z\in G\mid zg=gz\ (\forall g\in G)\}$（:ref[def-grp-center]）について $gzg^{-1}=z$ なので $Z(G)\trianglelefteq G$ である。同じ理由で $Z(G)$ に含まれる任意の部分群は $G$ で正規である。
3. **特殊線形群**：$A\in\GL_n(K)$、$B\in\SL_n(K)$ なら、行列式の積公式（:ref[thm-linalg-det-product]）より $\det(ABA^{-1})=\det A\det B(\det A)^{-1}=1$ なので $\SL_n(K)\trianglelefteq\GL_n(K)$ である。
4. **交代群**：偶置換全体 $A_n$ について、符号の積法則（:ref[thm-linalg-sign-product]）より $\sgn(\sigma\tau\sigma^{-1})=\sgn\sigma\,\sgn\tau\,\sgn\sigma^{-1}=\sgn\tau$ なので $A_n\trianglelefteq S_n$ である。
5. **Klein の四元群**：$\sigma\in S_4$ と $V_4$ の非単位元 $(a\,b)(c\,d)$ について、:ref[prop-grp-cycle-conjugation] と :eqref[eq-grp-conjugation-product] より
   $$
   \sigma(a\,b)(c\,d)\sigma^{-1}=\bigl(\sigma(a)\ \sigma(b)\bigr)\bigl(\sigma(c)\ \sigma(d)\bigr)
   $$
   である。$\{\sigma(a),\sigma(b),\sigma(c),\sigma(d)\}=\{1,2,3,4\}$ なので、右辺は互いに素な二つの互換の積であり、$V_4$ の非単位元である（$S_4$ の中でこの形の置換はちょうど三つあり、すべて $V_4$ に属する）。従って $V_4\trianglelefteq S_4$ である。
6. **二面体群の回転部分群**：$D_n=\langle r,s\rangle$ の部分群 $\langle r^k\rangle$ について、$rr^kr^{-1}=r^k$、$sr^ks^{-1}=r^{-k}\in\langle r^k\rangle$ である（:ref[prop-grp-dihedral-form] の $srs=r^{-1}$ より $sr^ks=r^{-k}$）。$D_n$ は有限なので判定法 5 により $\langle r^k\rangle\trianglelefteq D_n$ である。特に回転全体 $R=\langle r\rangle$ は正規である。
7. **四元数群**：$Q_8$ の元の位数は、$1$ が $1$、$-1$ が $2$、$\pm i,\pm j,\pm k$ が $4$ である。部分群 $L\le Q_8$ が位数 $4$ の元、例えば $i$ を含めば $\langle i\rangle=\{\pm1,\pm i\}\subseteq L$ であり、$L$ が $\langle i\rangle$ の外の元も含めば $|L|>4$ かつ $|L|\mid8$ より $L=Q_8$ である。位数 $4$ の元を含まない部分群は $\{\pm1\}$ に含まれる。従って $Q_8$ の部分群は
   $$
   \{1\},\quad\{\pm1\},\quad\langle i\rangle,\quad\langle j\rangle,\quad\langle k\rangle,\quad Q_8
   $$
   で尽きる。$\{\pm1\}=Z(Q_8)$ は 2 により正規であり、$\langle i\rangle,\langle j\rangle,\langle k\rangle$ は指数 $2$ なので次の :ref[prop-grp-index-two] により正規である。こうして $Q_8$ のすべての部分群は正規であるが、$ij=k\ne-k=ji$ なので $Q_8$ は可換でない。すべての部分群が正規であることは可換性よりも真に弱い条件である。
:::

:::example[正規でない部分群]{#exm-grp-non-normal}
1. $S_3$ の $H=\langle(1\,2)\rangle$ について、$(1\,3)(1\,2)(1\,3)^{-1}=(3\,2)\notin H$ なので $H$ は正規でない。:ref[exm-grp-cosets-s3] で左右の剰余類が異なったのはこのためである。
2. $D_n$（$n\ge3$）の $\langle s\rangle$ について、$sr^{-1}=rs$ より $rsr^{-1}=r^2s$ であり、$r^2\ne e$ なので $r^2s\ne s$ である。従って $\langle s\rangle$ は正規でない。
3. $\GL_2(K)$ の上三角群 $B$ について、$w=\begin{pmatrix}0&1\\1&0\end{pmatrix}$ とすると $w\begin{pmatrix}1&1\\0&1\end{pmatrix}w^{-1}=\begin{pmatrix}1&0\\1&1\end{pmatrix}\notin B$ なので $B$ は正規でない。
4. $S_4$ の $\langle(1\,2)(3\,4)\rangle$ について、$(2\,3)$ で共役すると $(1\,3)(2\,4)$ になるので正規でない。
:::

指数 $2$ の部分群は常に正規である。これは正規部分群を見つける最も簡単な方法の一つである。

:::proposition[指数 2 の部分群]{#prop-grp-index-two}
$[G:H]=2$ なら $H\trianglelefteq G$ である。
:::

:::proof
$g\in H$ なら $gH=H=Hg$ である。$g\notin H$ なら、左剰余類は $H$ と $gH$ の二つで $G$ を分割するので $gH=G\setminus H$ である。:ref[prop-grp-left-right-cosets] により右剰余類も二つであり、$g\notin H$ より $Hg\ne H$ なので $Hg=G\setminus H$ である。いずれの場合も $gH=Hg$ であり、判定法 3 により $H\trianglelefteq G$ である。
:::

例えば $A_n\trianglelefteq S_n$（$n\ge2$）、$R\trianglelefteq D_n$、$\langle i\rangle\trianglelefteq Q_8$ はすべてこの命題からも従う。指数 $3$ では同様の主張は成り立たない。$\langle(1\,2)\rangle\le S_3$ は指数 $3$ で正規でない。また :ref[exm-grp-a4-no-six] は「$A_4$ には指数 $2$ の部分群がない」と言い換えられる。

### 正規性と部分群の演算

:::remark{#rem-grp-normal-not-transitive}
正規性はどの群の中で考えるかに依存し、推移的でない。$D_4$ で $K=\{e,r^2,s,r^2s\}$ とすると、$K$ は指数 $2$ なので $K\trianglelefteq D_4$ であり、$\langle s\rangle$ は $K$ の中で指数 $2$ なので $\langle s\rangle\trianglelefteq K$ である。しかし :ref[exm-grp-non-normal] の 2 により $\langle s\rangle$ は $D_4$ で正規でない。同様に $\langle(1\,2)(3\,4)\rangle\trianglelefteq V_4\trianglelefteq S_4$ であるが、$\langle(1\,2)(3\,4)\rangle$ は $S_4$ で正規でない。「正規部分群の正規部分群」を扱うときにはこの点に注意する必要があり、[可解群と冪零群](/math/group-theory/solvable-nilpotent/) の組成列はこの状況を扱う。
:::

二つの部分群の積集合 $HK$ が部分群になる条件を調べる。

:::proposition[部分群の積]{#prop-grp-product-subgroup}
$H,K\le G$ とする。

1. $HK$ が部分群であることと $HK=KH$ は同値である。このとき $HK=\langle H\cup K\rangle$ である。
2. すべての $h\in H$ について $hKh^{-1}=K$ ならば $HK=KH$ である。特に $H$ または $K$ が正規なら $HK$ は部分群である。
:::

:::proof
1. $HK$ が部分群なら逆元をとる操作で閉じているので $HK=(HK)^{-1}=K^{-1}H^{-1}=KH$ である。逆に $HK=KH$ とする。$HK$ は $e=ee$ を含む。$h,h'\in H$、$k,k'\in K$ に対し、$k(k')^{-1}(h')^{-1}\in KH=HK$ なので、これを $h''k''$ と書くと
$$
(hk)(h'k')^{-1}=h\,k(k')^{-1}(h')^{-1}=hh''k''\in HK
$$
となり、:ref[prop-grp-subgroup-test] により $HK$ は部分群である。$HK$ は $H$ と $K$ を含み、$H\cup K$ を含む部分群は積 $hk$ をすべて含むので、$HK=\langle H\cup K\rangle$ である。

2. $hk=(hkh^{-1})h\in KH$ かつ $kh=h(h^{-1}kh)\in HK$ なので $HK=KH$ である。$K$ が正規ならこの仮定が成り立つ。$H$ が正規なら $H$ と $K$ の役割を入れ替えて $KH=HK$ を得る。
:::

$HK=KH$ は集合としての等式であり、$H$ の元と $K$ の元が可換であることを意味しない。例えば $S_3=A_3\langle(1\,2)\rangle$ では $A_3$ が正規なので積は部分群（全体）だが、$(1\,2\,3)$ と $(1\,2)$ は可換でない。一方 :ref[exm-grp-product-size] の $\langle(1\,2)\rangle\langle(2\,3)\rangle$ は $(1\,2\,3)$ を含むが $(1\,3\,2)=(2\,3)(1\,2)$ を含まないので $HK\ne KH$ である。

:::proposition[正規部分群と部分群の演算]{#prop-grp-normal-operations}
1. 正規部分群の族 $\{N_\lambda\}_{\lambda\in\Lambda}$ の共通部分 $\bigcap_\lambda N_\lambda$ は正規部分群である。
2. $N\trianglelefteq G$、$H\le G$ なら $H\cap N\trianglelefteq H$ である。
3. $N\trianglelefteq G$、$N\le H\le G$ なら $N\trianglelefteq H$ である。
4. $N,M\trianglelefteq G$ なら $NM=MN$ は $G$ の正規部分群である。
5. $N,M\trianglelefteq G$ かつ $N\cap M=\{e\}$ なら、すべての $n\in N$、$m\in M$ について $nm=mn$ である。
:::

:::proof
1. 部分群の共通部分は部分群である。$x\in\bigcap_\lambda N_\lambda$、$g\in G$ なら、各 $\lambda$ で $gxg^{-1}\in N_\lambda$ なので $gxg^{-1}\in\bigcap_\lambda N_\lambda$ である。

2. $H\cap N$ は部分群であり、$h\in H$、$x\in H\cap N$ なら $hxh^{-1}$ は $H$ の元の積なので $H$ に属し、$N$ の正規性から $N$ にも属する。

3. 判定法 2 の条件を $g\in H$ に制限すればよい。

4. :ref[prop-grp-product-subgroup] の 2 により $NM=MN$ は部分群である。$g\in G$ について $g(nm)g^{-1}=(gng^{-1})(gmg^{-1})\in NM$ なので正規である。

5. $x=nmn^{-1}m^{-1}$ とおく。$x=(nmn^{-1})m^{-1}$ と見れば $M$ の正規性から $x\in M$、$x=n(mn^{-1}m^{-1})$ と見れば $N$ の正規性から $x\in N$ である。よって $x\in N\cap M=\{e\}$ であり、$nm=mn$ を得る。
:::

5 は、二つの正規部分群が自明に交わるとき群がそれらの「直積」に近い形になることを示唆する。これは [直積・半直積と有限 Abel 群](/math/group-theory/products-abelian/) で内部直積として扱う。

### 正規化群

正規でない部分群 $H$ についても、$H$ を正規部分群として含む最大の部分群を考えることができる。

:::proposition[正規化群]{#prop-grp-normalizer}
$H\le G$ に対し
$$
N_G(H)=\{g\in G\mid gHg^{-1}=H\}
$$
とおき、$H$ の**正規化群**という（[群の作用](/math/group-theory/group-actions/) で共役作用の固定部分群として再び現れる）。

1. $N_G(H)$ は $H$ を含む $G$ の部分群であり、$H\trianglelefteq N_G(H)$ である。
2. $H\le K\le G$ かつ $H\trianglelefteq K$ なら $K\le N_G(H)$ である。すなわち $N_G(H)$ は $H$ を正規部分群として含む最大の部分群である。
3. $H\trianglelefteq G$ であることと $N_G(H)=G$ は同値である。
4. $gN_G(H)\mapsto gHg^{-1}$ は、$G/N_G(H)$ から $H$ の共役部分群全体の集合への全単射である。特に $H$ の相異なる共役部分群の個数は $[G:N_G(H)]$ に等しく、$G$ が有限ならそれは $[G:H]$ を割る。
5. $K\le N_G(H)$ なら $HK=KH$ は部分群である。
:::

:::proof
1. $N_G(H)$ は判定法の証明に現れた $M$ と同じものであり、部分群であることはそこで示した。$h\in H$ なら $hHh^{-1}=H$ なので $H\le N_G(H)$ であり、$N_G(H)$ の元はすべて $H$ を共役で保つから $H\trianglelefteq N_G(H)$ である。

2. $H\trianglelefteq K$ なら、各 $k\in K$ について $kHk^{-1}=H$ なので $k\in N_G(H)$ である。

3. 定義から明らかである。

4. $g,g'\in G$ について
$$
gHg^{-1}=g'Hg'^{-1}\iff (g^{-1}g')H(g^{-1}g')^{-1}=H\iff g^{-1}g'\in N_G(H)\iff gN_G(H)=g'N_G(H)
$$
であるから、写像は代表元によらずに定まり、単射である。全射は明らかである。$G$ が有限なら $H\le N_G(H)\le G$ に :ref[prop-grp-index-tower] を適用して $[G:H]=[G:N_G(H)]\,[N_G(H):H]$ を得る。

5. 各 $k\in K$ が $kHk^{-1}=H$ を満たすので、:ref[prop-grp-product-subgroup] の 2（$H$ と $K$ の役割を入れ替えたもの）による。
:::

:::example[正規化群の計算]{#exm-grp-normalizer}
1. $S_3$ で $H=\langle(1\,2)\rangle$ とする。:ref[prop-grp-cycle-conjugation] より $\sigma(1\,2)\sigma^{-1}=(\sigma(1)\,\sigma(2))$ はすべての互換を実現するので、$H$ の共役部分群は $\langle(1\,2)\rangle,\langle(1\,3)\rangle,\langle(2\,3)\rangle$ の三つである。よって $[S_3:N_{S_3}(H)]=3$ であり、$N_{S_3}(H)=H$ である。
2. $D_4$ で $H=\langle s\rangle$ とする。$r^ksr^{-k}=r^{2k}s$、$(r^ks)s(r^ks)^{-1}=r^ksr^{-k}=r^{2k}s$ なので、$s$ の共役は $s,r^2s$ の二つであり、$H$ の共役部分群は $\langle s\rangle,\langle r^2s\rangle$ の二つである。従って $N_{D_4}(H)$ は指数 $2$、位数 $4$ であり、$r^2sr^{-2}=s$ より $r^2\in N_{D_4}(H)$ なので $N_{D_4}(H)=\{e,r^2,s,r^2s\}$ である。
3. $S_4$ で $H=\langle(1\,2\,3\,4)\rangle$ とする。$4$ 巡回置換は $3!=6$ 個あり（$(1\,a\,b\,c)$ と $1$ から書き始めると、$a,b,c$ は残り三文字の並べ方で決まる）、各位数 $4$ の巡回部分群はそのうち二つ（$c$ と $c^{-1}=c^3$）を含むので、$4$ 巡回置換が生成する部分群は三つある。:ref[prop-grp-cycle-conjugation] より $4$ 巡回置換の共役は $4$ 巡回置換であり、$(a\,b\,c\,d)$ は $1\mapsto a,2\mapsto b,3\mapsto c,4\mapsto d$ となる $\sigma$ による $(1\,2\,3\,4)$ の共役なので、三つの部分群はすべて $H$ と共役である。従って $[S_4:N_{S_4}(H)]=3$、$|N_{S_4}(H)|=8$ である。実際、$(1\,3)(1\,2\,3\,4)(1\,3)^{-1}=(3\,2\,1\,4)=(1\,4\,3\,2)\in H$ なので $(1\,3)\in N_{S_4}(H)$ である。$\langle H,(1\,3)\rangle$ は $N_{S_4}(H)$ の部分群で、$H$ と $H$ の外の元 $(1\,3)$ を含むから位数は $4$ より大きい $8$ の約数、すなわち $8$ であり、$N_{S_4}(H)=\langle(1\,2\,3\,4),(1\,3)\rangle$ である。これは頂点に $1,2,3,4$ と番号をつけた正方形の対称群にあたる。
:::

### 正規核と正規閉包

一般の部分群 $H\le G$ から正規部分群を作るには、$H$ を縮める方向と拡げる方向がある。

:::proposition[正規核と正規閉包]{#prop-grp-normal-core-closure}
$H\le G$ とする。

1. $\operatorname{Core}_G(H)=\bigcap_{g\in G}gHg^{-1}$ は $H$ に含まれる正規部分群であり、$H$ に含まれる $G$ の正規部分群はすべてこれに含まれる。
2. $H^G=\langle ghg^{-1}\mid g\in G,\ h\in H\rangle$ は $H$ を含む正規部分群であり、$H$ を含む $G$ の正規部分群はすべてこれを含む。$H^G$ は $H$ を含む正規部分群すべての共通部分に等しい。
:::

:::proof
1. $g=e$ の項があるので $\operatorname{Core}_G(H)\subseteq H$ であり、部分群の共通部分なので部分群である。$x\in G$ による共役 $y\mapsto xyx^{-1}$ は全単射なので共通部分と交換し、
$$
x\Bigl(\bigcap_{g}gHg^{-1}\Bigr)x^{-1}=\bigcap_{g}(xg)H(xg)^{-1}=\bigcap_{g'}g'Hg'^{-1}
$$
である（$g'=xg$ は $g$ とともに $G$ 全体を動く）。よって正規である。$L\trianglelefteq G$、$L\le H$ なら、すべての $g$ について $L=gLg^{-1}\le gHg^{-1}$ なので $L\le\operatorname{Core}_G(H)$ である。

2. 生成系 $X=\{ghg^{-1}\}$ について $x(ghg^{-1})x^{-1}=(xg)h(xg)^{-1}\in X$ である。$H^G$ の元は $X$ の元とその逆元の積なので、:eqref[eq-grp-conjugation-product] によりその共役も $X$ の元とその逆元の積であり、$H^G$ に属する。よって $H^G\trianglelefteq G$ であり、$g=e$ の項から $H\subseteq H^G$ である。$H$ を含む正規部分群は $X$ のすべての元を含むので $H^G$ を含む。$H^G$ 自身が $H$ を含む正規部分群なので、それらすべての共通部分は $H^G$ に等しい。
:::

$\operatorname{Core}_G(H)$ を $H$ の**正規核**、$H^G$ を $H$ の**正規閉包**という。正規閉包は [自由群と群の表示](/math/group-theory/free-groups/) で関係式から群を定める際に中心的な役割を果たす。

:::example[正規核と正規閉包の計算]{#exm-grp-normal-core-closure}
$S_3$ で $H=\langle(1\,2)\rangle$ とすると、その共役は三つの互換部分群である（:ref[exm-grp-normalizer]）。共通部分は $\{e\}$、三つの互換が生成する部分群は $S_3$ なので、$\operatorname{Core}_{S_3}(H)=\{e\}$、$H^{S_3}=S_3$ である。

$D_4$ で $H=\langle s\rangle$ とすると、$s$ の共役は $s,r^2s$ なので $H^{D_4}=\langle s,r^2s\rangle=\{e,r^2,s,r^2s\}$ であり、$\operatorname{Core}_{D_4}(H)=\langle s\rangle\cap\langle r^2s\rangle=\{e\}$ である。正規閉包は全群になるとは限らない。
:::

## 剰余群

$\Z/n\Z$ では、類 $a+n\Z$ と $b+n\Z$ の和を代表元の和の類 $(a+b)+n\Z$ と定めた。代表元を $a'=a+kn$、$b'=b+ln$ に取り替えても $a'+b'=(a+b)+(k+l)n$ なので結果は変わらない。同じことを一般の群の左剰余類で行うには、$(gN)(hN)=(gh)N$ が代表元によらないことが必要である。

部分集合 $A,B\subseteq G$ に対して $AB=\{ab\mid a\in A,\ b\in B\}$ と定める。$N\trianglelefteq G$ なら、判定法 3 と $NN=N$ を用いて
$$
(gN)(hN)=g(Nh)N=g(hN)N=gh(NN)=ghN
$$
となる。つまり正規部分群の二つの左剰余類の積集合は再び左剰余類であり、代表元を用いずに積が定まる。正規でない場合は、$(gH)(hH)$ が一つの左剰余類になるとは限らない。

### 剰余群の構成

:::theorem[剰余群の構成]{#thm-grp-quotient}
$N\le G$ とする。

1. $G/N$ 上の規則 $(gN)(hN)=(gh)N$ が代表元の選び方によらずに定まるための必要十分条件は $N\trianglelefteq G$ である。
2. $N\trianglelefteq G$ のとき、$G/N$ はこの積について群をなす。単位元は $N=eN$、$gN$ の逆元は $g^{-1}N$ であり、$|G/N|=[G:N]$ である。
3. 写像 $\pi_N:G\to G/N$、$\pi_N(g)=gN$ は全射であり、$\pi_N(gh)=\pi_N(g)\pi_N(h)$ を満たす。$\pi_N(g)$ が $G/N$ の単位元であることと $g\in N$ は同値である。

2 の群 $G/N$ を $N$ による :term[剰余群]{en="quotient group" yomi="じょうよぐん"}（商群）といい、$\pi_N$ を**自然な射影**という。
:::

:::proof
1. $N\trianglelefteq G$ とし、$g'N=gN$、$h'N=hN$ とする。$g'=gn$、$h'=hm$（$n,m\in N$）と書けるので
$$
g'h'=gnhm=gh\,(h^{-1}nh)\,m
$$
である。正規性より $h^{-1}nh\in N$ なので $g'h'\in ghN$、すなわち $g'h'N=ghN$ である。

逆に積が代表元によらないとする。$n\in N$、$h\in G$ をとる。$nN=eN$ なので、$(nN)(hN)$ と $(eN)(hN)$ は等しく、$nhN=hN$ である。:ref[prop-grp-coset-test] の 1 より $h^{-1}nh\in N$ である。これがすべての $h,n$ で成り立つので、$h$ を $h^{-1}$ に替えて判定法 2 の条件を得る。よって $N\trianglelefteq G$ である。

2. 結合法則は $\bigl((gN)(hN)\bigr)(kN)=(gh)kN=g(hk)N=(gN)\bigl((hN)(kN)\bigr)$ による。$(eN)(gN)=gN=(gN)(eN)$、$(gN)(g^{-1}N)=eN=(g^{-1}N)(gN)$ である。$G/N$ の元の個数は定義により $[G:N]$ である。

3. 全射であることと積を保つことは定義から明らかである。$gN=N$ と $g\in N$ の同値は :ref[prop-grp-coset-test] の 1 による。
:::

:::remark
剰余群の元は $G$ の部分集合（剰余類）であり、$G$ の元ではない。単位元は $N$ という集合であり、「$N$ の元をすべて単位元とみなした」群と考えることができる。$\pi_N$ の性質は次章の言葉では「全射準同型であり核は $N$」と表される（[準同型と同型定理](/math/group-theory/homomorphisms/)）。
:::

:::proposition[剰余群の基本性質]{#prop-grp-quotient-basic}
$N\trianglelefteq G$ とする。

1. すべての整数 $k$ について $(gN)^k=g^kN$ である。
2. $gN$ の位数は $g^k\in N$ を満たす最小の正整数 $k$ である（存在しなければ無限位数）。$g$ の位数が有限なら $\operatorname{ord}(gN)$ は $\operatorname{ord}(g)$ を割り、$G/N$ が有限なら $\operatorname{ord}(gN)$ は $[G:N]$ を割る。
3. $G$ が可換なら $G/N$ も可換である。$G=\langle S\rangle$ なら $G/N=\langle sN\mid s\in S\rangle$ であり、特に $G$ が巡回群なら $G/N$ も巡回群である。
4. $G/N$ が可換であることと、すべての $g,h\in G$ について $ghg^{-1}h^{-1}\in N$ が成り立つことは同値である。
:::

:::proof
1. $k=0$ では両辺とも $N$ である。$k\ge0$ で成り立てば $(gN)^{k+1}=(g^kN)(gN)=g^{k+1}N$ なので、帰納法により $k\ge0$ で成り立つ。$k<0$ なら $(gN)^k=\bigl((gN)^{-k}\bigr)^{-1}=(g^{-k}N)^{-1}=g^kN$ である。

2. 1 より $(gN)^k=N$ は $g^k\in N$ と同値なので前半を得る。$g^m=e$ なら $(gN)^m=N$ であり、:ref[thm-grp-cyclic-order] より $\operatorname{ord}(gN)\mid m$ である。最後の主張は :ref[cor-grp-order-divides] を $G/N$ に適用したものである。

3. $(gN)(hN)=ghN=hgN=(hN)(gN)$ である。$G$ の元は $S$ の元とその逆元の有限積 $s_1^{\epsilon_1}\cdots s_m^{\epsilon_m}$（:ref[prop-grp-generated-words]）なので、その類は $(s_1N)^{\epsilon_1}\cdots(s_mN)^{\epsilon_m}$ である。

4. $ghN=hgN$ は $(hg)^{-1}gh=g^{-1}h^{-1}gh\in N$ と同値である。これがすべての $g,h$ で成り立つことは、$g,h$ をそれぞれ $g^{-1},h^{-1}$ に置き換えれば、すべての $g,h$ について $ghg^{-1}h^{-1}\in N$ が成り立つことと同値である。
:::

4 の元 $ghg^{-1}h^{-1}$ を $g,h$ の**交換子**という。交換子全体が生成する部分群は「割ると可換になる最小の正規部分群」であり、[可解群と冪零群](/math/group-theory/solvable-nilpotent/) で詳しく扱う。

### 剰余群の例

:::example[剰余群の例]{#exm-grp-quotient-examples}
1. **整数の剰余類**：$\Z/n\Z$ の剰余群としての和 $(a+n\Z)+(b+n\Z)=(a+b)+n\Z$ は、前章で定めた $\bar a+\bar b=\overline{a+b}$（:ref[exm-grp-integers-mod]）にほかならない。巡回群 $C_n$ は $\Z$ の剰余群として得られる。
2. **巡回群の剰余群**：$G=C_{12}$、$N=\langle\bar4\rangle=\{\bar0,\bar4,\bar8\}$ とすると $|G/N|=12/3=4$ であり、$G/N$ は $\bar1+N$ で生成される巡回群である（:ref[prop-grp-quotient-basic] の 3）。$\bar2+N$ の位数は、$2\cdot\bar2=\bar4\in N$ より $2$ である。$\bar1$ は $G$ で位数 $12$ だが、$\bar1+N$ の位数は $4$ に下がる。
3. **回転と鏡映**：$D_n/R$（$R=\langle r\rangle$）は $R$ と $sR$ の二元からなる位数 $2$ の群である。積は「回転と回転、鏡映と鏡映の積は回転、回転と鏡映の積は鏡映」という規則を表す。
4. **行列式**：$\GL_n(K)/\SL_n(K)$ の元は「行列式が $c$ の行列全体」であり（:ref[exm-grp-cosets-linear]）、積は $\det(AB)=\det A\det B$ に従う。剰余群の積は行列式の値の積と対応する。この対応は次章で第一同型定理として定式化される。
5. **円周**：$\R/\Z$ の和 $(x+\Z)+(y+\Z)=(x+y)+\Z$ は、代表元を $[0,1)$ にとれば「和をとって整数部分を捨てる」演算である。$n(x+\Z)=\Z$ は $nx\in\Z$ と同値である。従って $x=a/b$ を既約分数（$b>0$）で表すと、この条件を満たす最小の正整数は $b$ であり、$x+\Z$ の位数は $b$ である。逆に有限位数ならある正整数 $n$ に対して $nx\in\Z$ となるので $x\in\Q$ である。従って有限位数の元全体は部分群 $\Q/\Z$ をなす。$\frac1n+\Z$ の位数は $n$ なので、$\R/\Z$ はあらゆる位数の元をもつ。無理数 $x$ の類は無限位数である。
6. **四元数群を符号で割る**：$Z(Q_8)=\{\pm1\}$ は正規であり、$Q_8/\{\pm1\}$ は $\{\pm1\},\{\pm i\},\{\pm j\},\{\pm k\}$ の四元からなる。$i^2=-1\in\{\pm1\}$ などより非単位元の類はすべて位数 $2$ なので、:ref[prop-grp-order-four] により $Q_8/\{\pm1\}\cong V_4$ である。同型は $\{\pm i\}\mapsto(1\,2)(3\,4)$、$\{\pm j\}\mapsto(1\,3)(2\,4)$、$\{\pm k\}\mapsto(1\,4)(2\,3)$ で与えられる。$Q_8$ 自身は非可換だが、その剰余群は可換になる。
7. **$D_4$ を中心で割る**：$Z(D_4)=\{e,r^2\}$ である。実際、$r^2$ は $r$ と可換で、$sr^2=r^{-2}s=r^2s$ である。一方 $r^{\pm1}$ は $s$ と可換でなく（$sr=r^{-1}s\ne rs$）、鏡映 $r^is$ は $r$ と可換でない（$r\cdot r^is=r^{i+1}s\ne r^{i-1}s=r^is\cdot r$）。$D_4/Z(D_4)$ は位数 $4$ であり、回転の平方 $r^{2i}$ は $e$ か $r^2$ で $Z(D_4)$ に属し、鏡映の平方は $(r^is)^2=e$ なので、すべての元の平方が $Z(D_4)$ に属する。従って非単位元の類はすべて位数 $2$ であり、$D_4/Z(D_4)\cong V_4$ である。
8. **$S_4$ を $V_4$ で割る**：$S_4/V_4$ は位数 $6$ の群である。$(1\,2)(1\,3)=(1\,3\,2)$、$(1\,3)(1\,2)=(1\,2\,3)$ であり、$(1\,3\,2)^{-1}(1\,2\,3)=(1\,2\,3)^2=(1\,3\,2)\notin V_4$ なので $(1\,2)V_4$ と $(1\,3)V_4$ は可換でない。従って $S_4/V_4$ は非可換である。:ref[exm-grp-index-tower] より文字 $4$ を固定する六つの置換が代表元の系をなし、演習問題の結果（位数 $6$ の群の分類）から $S_4/V_4\cong S_3$ となる。
:::

:::example[格子を周期で割る]{#exm-grp-lattice-quotient}
平面格子 $\Z^2$（加法群）の部分群 $N=\langle(2,0),(1,2)\rangle=\{(2a+b,\,2b)\mid a,b\in\Z\}$ を考える。$(u,v)\in\Z^2$ から $(1,2)$ の整数倍を引いて第二成分を $0$ または $1$ にし、さらに $(2,0)$ の整数倍を引いて第一成分も $0$ または $1$ にできる。従って $(0,0),(1,0),(0,1),(1,1)$ が代表元の系を含む。二つの代表元の差 $(x,y)$ が $N$ に入るなら、$y=2b$ は偶数かつ $|y|\le1$ なので $y=0$、すると $b=0$ で $x=2a$ も偶数かつ $|x|\le1$ なので $x=0$ である。よって四元は相異なる類を代表し、$[\Z^2:N]=4$ である。

$x=(0,1)+N$ とおくと、$2x=(0,2)+N=(-1,0)+N$ は零でなく（$(-1,0)$ の第二成分 $0$ から $b=0$、第一成分 $-1$ は偶数でない）、$4x=(0,4)+N=(-2,0)+N=N$ である。$3x=(0,3)+N=(-1,1)+N$ も零でないので $x$ の位数は $4$ であり、$\Z^2/N$ は $x$ で生成される位数 $4$ の巡回群である。二つの座標をそれぞれ二通りに制限できたことだけからは、$C_2\times C_2$ か $C_4$ かは判断できない。
:::

### 中心による剰余群

剰余群 $G/Z(G)$ は、$G$ が「どれだけ可換から遠いか」を測る。次の命題は、この剰余群が非自明な巡回群にはならないことを示す。

:::proposition[中心による剰余群が巡回なら可換]{#prop-grp-center-quotient-cyclic}
$G/Z(G)$ が巡回群なら $G$ は可換群であり、従って $G=Z(G)$、$G/Z(G)$ は自明な群である。
:::

:::proof
$Z=Z(G)$ とおき、$G/Z=\langle aZ\rangle$ とする。任意の $x,y\in G$ について、$xZ=(aZ)^i=a^iZ$、$yZ=a^jZ$ となる整数 $i,j$ があるので、$x=a^iz$、$y=a^jw$（$z,w\in Z$）と書ける。$z,w$ はすべての元と可換なので
$$
xy=a^iza^jw=a^{i+j}zw=a^ja^iwz=a^jwa^iz=yx
$$
である。よって $G$ は可換であり、$Z(G)=G$ である。
:::

この命題から、$|G/Z(G)|$ は素数になりえない（素数位数の群は巡回群だから）。例えば位数 $8$ の非可換群 $D_4$、$Q_8$ では $Z(G)\ne G$ であり、$|Z(G)|$ は $1,2,4$ のいずれかである。$|Z(G)|=4$ なら $G/Z(G)$ は位数 $2$ の巡回群となって矛盾するので、$|Z(G)|\in\{1,2\}$ である。実際どちらも $|Z(G)|=2$ であり、$G/Z(G)\cong V_4$ であった（:ref[exm-grp-quotient-examples] の 6、7）。

### 積と両立する同値関係

剰余群は「群の演算と両立する同値関係による商」として特徴付けられる。

:::proposition[積と両立する同値関係]{#prop-grp-compatible-relation}
$G$ 上の同値関係 $\sim$ が、$a\sim a'$ かつ $b\sim b'$ ならば $ab\sim a'b'$ を満たすとき、$\sim$ は**積と両立する**という。

1. $N\trianglelefteq G$ に対し、$a\equiv b\pmod N:\Leftrightarrow a^{-1}b\in N$ は積と両立する同値関係であり、その同値類は左剰余類 $aN$ である。
2. 積と両立する同値関係 $\sim$ に対し、$N=\{x\in G\mid x\sim e\}$ は正規部分群であり、$a\sim b$ と $a^{-1}b\in N$ は同値である。
3. 1 と 2 の対応は、$G$ の正規部分群全体と、$G$ 上の積と両立する同値関係全体との間の全単射を与える。このとき商集合 $G/{\sim}$ 上に $[a][b]=[ab]$ で定まる演算は、剰余群 $G/N$ の演算と一致する。
:::

:::proof
1. 同値関係であることと同値類が $aN$ であることは :ref[prop-grp-coset-relation] による。両立性は :ref[thm-grp-quotient] の 1 の前半、すなわち $aN=a'N$、$bN=b'N$ なら $abN=a'b'N$ であることにほかならない。

2. 反射律より $e\in N$ である。$x,y\in N$ なら、両立性より $xy\sim ee=e$ なので $xy\in N$ である。$x\in N$ なら、$x^{-1}\sim x^{-1}$ と $x\sim e$ から $x^{-1}x\sim x^{-1}e$、すなわち $e\sim x^{-1}$ であり、対称律から $x^{-1}\in N$ である。従って $N$ は部分群である。$g\in G$、$x\in N$ なら、$g\sim g$、$x\sim e$、$g^{-1}\sim g^{-1}$ から $gxg^{-1}\sim geg^{-1}=e$ となるので、$N$ は正規である。

$a\sim b$ なら、$a^{-1}\sim a^{-1}$ と合わせて $e=a^{-1}a\sim a^{-1}b$ となり、$a^{-1}b\in N$ である。逆に $a^{-1}b\sim e$ なら、$a\sim a$ と合わせて $b=aa^{-1}b\sim ae=a$ である。

3. $N$ から作った同値関係 $\equiv\pmod N$ に 2 を適用すると、$\{x\mid x\equiv e\}=\{x\mid e^{-1}x\in N\}=N$ に戻る。逆に $\sim$ から 2 で $N$ を作ると、2 の後半により $\equiv\pmod N$ は $\sim$ と一致する。従って二つの対応は互いに逆である。$[a]=aN$ なので、$[a][b]=[ab]$ は剰余群の積 $(aN)(bN)=abN$ そのものである。
:::

:::remark
左からの積と両立する同値関係（$a\sim b$ ならすべての $g$ について $ga\sim gb$）は、正規とは限らない部分群 $H$ の左剰余類への分割 $\sim_H$ と一対一に対応する。実際、$H$ を $e$ の類とすると、$x,y\in H$ なら $x\sim y$ に左から $x^{-1}$ を掛けて $e\sim x^{-1}y$ を得るので $H$ は部分群である。また $a\sim b$ は $e\sim a^{-1}b$ と同値であり、元の同値類は $aH$ と一致する。逆に左剰余類が定める関係は $(ga)^{-1}(gb)=a^{-1}b$ により左からの積と両立する。$H$ が正規でない場合にも、$G$ の元を左から掛ける操作 $g\cdot(kH)=(gk)H$ は定義できる。これが [群の作用](/math/group-theory/group-actions/) の典型例である。
:::

### 正規性が破れると何が起こるか

正規でない部分群で剰余類の積を定義しようとすると、同じ類どうしを掛けたはずなのに異なる答えが出る。$S_3$ で $H=\langle(1\,2)\rangle$ とし、類 $H$ と類 $(1\,3)H$ を掛ける。$H$ の代表元として $e$ を使えば $e(1\,3)H=(1\,3)H$ であり、$(1\,2)$ を使えば $(1\,2)(1\,3)H$ である。両者が等しいためには $(1\,3)^{-1}(1\,2)(1\,3)=(2\,3)$ が $H$ に入る必要があるが、実際には入らない。失敗したのは $G$ の分割ではなく、その上に新しく定めようとした積の規則である。:ref[thm-grp-quotient] の証明の後半は、この計算を一般化したものである。

## 演習問題

::::exercise{level=basic}
$C_{12}=\Z/12\Z$ の部分群 $H=\langle\bar3\rangle$ の剰余類をすべて列挙し、$[C_{12}:H]$ を求めよ。

:::solution
$H=\{\bar0,\bar3,\bar6,\bar9\}$ である。$\bar1+H=\{\bar1,\bar4,\bar7,\bar{10}\}$、$\bar2+H=\{\bar2,\bar5,\bar8,\bar{11}\}$ であり、$H$ と合わせて $C_{12}$ の $12$ 元をちょうど一度ずつ含むので、剰余類はこの三つである。各元はその代表元を $3$ で割った余りによってただ一つの類に入る。$[C_{12}:H]=3=12/4$ であり、Lagrange の定理と一致する。
:::
::::

::::exercise{level=basic}
$S_3$ に位数 $4$ の部分群は存在するか。位数 $3$ の部分群はいくつあるか。

:::solution
$4$ は $|S_3|=6$ を割らないので、Lagrange の定理により位数 $4$ の部分群は存在しない。位数 $3$ の部分群は素数位数なので巡回群であり（:ref[cor-grp-order-divides]）、位数 $3$ の元で生成される。$S_3$ の位数 $3$ の元は $(1\,2\,3)$ と $(1\,3\,2)=(1\,2\,3)^2$ だけであり、両者は同じ部分群 $A_3=\{e,(1\,2\,3),(1\,3\,2)\}$ を生成する。従って位数 $3$ の部分群はちょうど一つである。
:::
::::

::::exercise{level=basic}
$D_6$ において $N=\langle r^3\rangle$ が正規部分群であることを示し、$D_6/N$ の位数と、$rN$、$sN$ の位数を求めよ。また $D_6/N$ が可換でないことを示せ。

:::solution
$sr^3s^{-1}=r^{-3}=r^3$ であり、$r$ は $r^3$ と可換なので、判定法（:ref[prop-grp-normal-criteria] の 5）により $N=\{e,r^3\}\trianglelefteq D_6$ である（実際 $r^3$ は中心に属する）。$|D_6/N|=12/2=6$ である。$r^k\in N$ となる最小の正整数は $k=3$ なので $rN$ の位数は $3$、$s\notin N$ かつ $s^2=e$ なので $sN$ の位数は $2$ である。$(sN)(rN)=srN=r^{-1}sN$、$(rN)(sN)=rsN$ であり、両者が等しいことは $(rs)^{-1}r^{-1}s=sr^{-2}s=r^2\in N$ と同値であるが、$r^2\notin N$ なので等しくない。従って $D_6/N$ は可換でない。
:::
::::

::::exercise{level=basic}
指数の乗法性を用いて $[\Z:18\Z]$ を $6\Z$ を経由して求めよ。

:::solution
$18\Z\le6\Z\le\Z$ である。$[\Z:6\Z]=6$ である（:ref[exm-grp-cosets-integers]）。$6\Z$ の元 $6k$ は $6k-6k'\in18\Z$、すなわち $k\equiv k'\pmod3$ のときに同じ類に属するので、$6\Z/18\Z$ は $0,6,12$ の類からなり $[6\Z:18\Z]=3$ である。:ref[prop-grp-index-tower] より $[\Z:18\Z]=6\cdot3=18$ である。
:::
::::

::::exercise{level=basic}
Euler の定理または Fermat の小定理を用いて、$5^{123}$ を $11$ で割った余りと、$2^{1000}$ を $21$ で割った余りを求めよ。

:::solution
$11$ は素数で $11\nmid5$ なので $5^{10}\equiv1\pmod{11}$ である。$123=10\cdot12+3$ より $5^{123}\equiv5^3=125=11\cdot11+4\equiv4\pmod{11}$ である。

$21=3\cdot7$ と互いに素な $1$ 以上 $21$ 以下の整数は、$3$ の倍数 $7$ 個と $7$ の倍数 $3$ 個（$21$ は両方に数えられる）を除いた $21-7-3+1=12$ 個なので $\varphi(21)=12$ であり、$2^{12}\equiv1\pmod{21}$ である。$1000=12\cdot83+4$ より $2^{1000}\equiv2^4=16\pmod{21}$ である。なお実際には $2^6=64=3\cdot21+1$ より $2$ の位数は $6$ であり、$1000\equiv4\pmod6$ からも同じ答えを得る。
:::
::::

::::exercise{level=basic}
$H$ が有限群 $G$ の位数 $m$ の部分群のうちただ一つのものであれば、$H\trianglelefteq G$ であることを示せ。これを用いて、$C_n$ の任意の部分群と、$Q_8$ の部分群 $\{\pm1\}$ が正規であることを再確認せよ。

:::solution
$g\in G$ に対し $gHg^{-1}$ は部分群であり、$|gHg^{-1}|=|H|=m$ である（:ref[def-grp-normal] の直後の注意）。位数 $m$ の部分群は $H$ だけなので $gHg^{-1}=H$ であり、$H\trianglelefteq G$ である。

$C_n$ では各約数 $d$ に対して位数 $d$ の部分群がただ一つなので（:ref[thm-grp-cyclic-subgroups]）、すべての部分群が正規である（可換性からも分かる）。$Q_8$ の位数 $2$ の元は $-1$ だけなので、位数 $2$ の部分群は $\{\pm1\}$ だけであり、正規である。
:::
::::

::::exercise{level=basic}
有限群 $G$ の部分群 $H,K$ の位数が互いに素なら $H\cap K=\{e\}$ であることを示せ。

:::solution
$H\cap K$ は $H$ と $K$ の両方の部分群なので、Lagrange の定理によりその位数は $|H|$ と $|K|$ の公約数である。仮定から公約数は $1$ だけなので $|H\cap K|=1$、すなわち $H\cap K=\{e\}$ である。
:::
::::

::::exercise{level=basic}
$H\le G$、$x\in G$ とする。左剰余類 $xH$ が $G$ の部分群であるための必要十分条件は $x\in H$ であることを示せ。

:::solution
$xH$ が部分群なら単位元を含むので $e=xh$（$h\in H$）と書け、$x=h^{-1}\in H$ である。逆に $x\in H$ なら :ref[prop-grp-coset-test] の 1 により $xH=H$ であり、部分群である。
:::
::::

::::exercise{level=advanced}
$n\ge2$ とし、$S_n$ の部分群 $H=\{\sigma\mid\sigma(n)=n\}$ を考える。両側剰余類 $H\backslash S_n/H$ をすべて求め、それぞれの元の個数を $n$ で表せ。

:::solution
$g(n)=n$ なら $g\in H$ であり $HgH=H$ である。$g(n)\ne n$ とする。右から $H$ の元を掛けても $n$ の行き先は変わらない（$(gk)(n)=g(n)$）。左から $h\in H$ を掛けると $n$ の行き先は $h(g(n))$ になり、$h$ は $n$ を固定して $\{1,\ldots,n-1\}$ を任意に並べ替えるので、行き先は $n$ 以外の任意の文字にできる。従って $g(n)\ne n$、$k(n)\ne n$ なら $h(g(n))=k(n)$ となる $h\in H$ があり、$(k^{-1}hg)(n)=n$ より $k^{-1}hg\in H$、すなわち $g\in HkH$ である。

よって両側剰余類は $H$ と $S_n\setminus H$ の二つであり、元の個数は $(n-1)!$ と $n!-(n-1)!=(n-1)\cdot(n-1)!$ である。後者は :ref[prop-grp-double-size] からも確かめられる。$g=(n-1\ \ n)$ とすると $gHg^{-1}$ は $n-1$ を固定する置換全体であり、$H\cap gHg^{-1}$ は $n-1,n$ をともに固定する置換全体で位数 $(n-2)!$ である。従って $|HgH|=(n-1)!^2/(n-2)!=(n-1)\cdot(n-1)!$ である。$n=3$ では $2$ 元と $4$ 元になり、:ref[exm-grp-double-coset-s3] と一致する。
:::
::::

::::exercise{level=advanced}
$N\trianglelefteq G$ の位数を $m$（有限）とし、$g\in G$ について $gN$ の位数が $k$（有限）であるとする。$g$ の位数 $d$ は有限で、$k$ の倍数かつ $km$ の約数であることを示せ。さらに両端 $d=k$ と $d=km$ がともに起こる例を与えよ。

:::solution
$(gN)^k=N$ より $g^k\in N$ であり、:ref[cor-grp-order-divides] を $N$ に適用すると $(g^k)^m=e$ である。従って $g$ の位数 $d$ は有限で $km$ を割る（:ref[thm-grp-cyclic-order]）。また $g^d=e\in N$ より $(gN)^d=N$ なので $k\mid d$ である。

例えば $G=C_6$、$N=\langle\bar3\rangle=\{\bar0,\bar3\}$ とすると $m=2$ である。$\bar1+N$ と $\bar2+N$ はどちらも位数 $3$（$\bar1,\ 2\cdot\bar1=\bar2\notin N$、$3\cdot\bar1=\bar3\in N$、および $\bar2,\ 2\cdot\bar2=\bar4\notin N$、$3\cdot\bar2=\bar0\in N$）なので $k=3$ であるが、$\bar1$ の位数は $6=km$、$\bar2$ の位数は $3=k$ である。
:::
::::

::::exercise{level=advanced}
位数 $6$ の群 $G$ は $C_6$ または $S_3$ と同型であることを示せ。

:::hint
位数 $2$ の元は :ref[prop-grp-order-two-element] から得られる。位数 $3$ の元がないと仮定して矛盾を導き、位数 $3$ の元 $a$ が生成する部分群が正規であることを使う。
:::

:::solution
$G$ に位数 $6$ の元があれば $G$ は巡回群であり、:ref[thm-grp-cyclic-classification] により $C_6$ と同型である。以下、位数 $6$ の元はないとする。Lagrange の定理より非単位元の位数は $2$ または $3$ である。

**位数 $3$ の元の存在**：すべての非単位元の位数が $2$ だとすると、すべての $x$ について $x=x^{-1}$ なので、$ab=(ab)^{-1}=b^{-1}a^{-1}=ba$ となり $G$ は可換である。相異なる非単位元 $a,b$ をとると、$ab\ne e,a,b$ であり（:ref[prop-grp-order-four] の証明と同じ議論）、$\{e,a,b,ab\}$ は積と逆元で閉じた位数 $4$ の部分群になる。$4\nmid6$ なのでこれは Lagrange の定理に反する。よって位数 $3$ の元 $a$ が存在する。

**位数 $2$ の元**：:ref[prop-grp-order-two-element] により位数 $2$ の元 $b$ が存在する。

**構造の決定**：$N=\langle a\rangle=\{e,a,a^2\}$ は指数 $2$ なので正規である（:ref[prop-grp-index-two]）。$b$ の位数は $2$ で $N$ の非単位元の位数は $3$ なので $b\notin N$ であり、$N$ の正規性から $bN=Nb$ なので $G=N\sqcup Nb=\{e,a,a^2\}\sqcup\{b,ab,a^2b\}$ である。従って $G$ の元は $a^ib^j$（$0\le i<3$、$0\le j<2$）の六つで、これらは相異なる。$N$ の正規性より $bab^{-1}\in N$ であり、共役は位数を保つので（:eqref[eq-grp-conjugation-product] より $(bab^{-1})^k=ba^kb^{-1}$）、$bab^{-1}\in\{a,a^2\}$ である。$bab^{-1}=a$ なら $a,b$ は可換であり、位数 $3$ と $2$ は互いに素なので :ref[prop-grp-coprime-orders] により $ab$ の位数は $6$ となって仮定に反する。従って $bab^{-1}=a^{-1}$、すなわち $ba=a^{-1}b$ であり、帰納的に $ba^k=a^{-k}b$ を得る。これから
$$
(a^ib^j)(a^kb^l)=a^{i+(-1)^jk}b^{j+l}
$$
（指数はそれぞれ法 $3$、法 $2$）が成り立ち、$G$ の積表は完全に決まる。

$S_3$ で $\alpha=(1\,2\,3)$、$\beta=(1\,2)$ とすると、$\alpha$ は位数 $3$、$\beta$ は位数 $2$ で、:ref[prop-grp-cycle-conjugation] より $\beta\alpha\beta^{-1}=(2\,1\,3)=(1\,3\,2)=\alpha^{-1}$ である。$S_3$ も上の仮定（位数 $6$ の元をもたない位数 $6$ の群）を満たすので、同じ議論により元は $\alpha^i\beta^j$ の六つで、同じ積の規則に従う。従って $a^ib^j\mapsto\alpha^i\beta^j$ は積を保つ全単射であり、$G\cong S_3$ である。最後に、$C_6$ は可換で $S_3$ は可換でないので、両者は同型でない。
:::
::::

::::exercise{level=advanced}
素数 $p$ に対して $(p-1)!\equiv-1\pmod p$ が成り立つこと（Wilson の定理）を、群 $(\Z/p\Z)^\times$ を用いて示せ。

:::solution
$p=2$ なら $1!=1\equiv-1\pmod2$ である。以下 $p$ は奇素数とする。$(\Z/p\Z)^\times=\{\bar1,\ldots,\overline{p-1}\}$ は位数 $p-1$ の可換群である。

$x=x^{-1}$、すなわち $x^2=\bar1$ となる元を求める。$x=\bar a$ とすると条件は $p\mid(a-1)(a+1)$ である。素数 $p$ が積 $uv$ を割るなら $p\mid u$ または $p\mid v$ である（$\Z$ の部分群 $u\Z+p\Z$ は :ref[thm-grp-cyclic-subgroups] により $d\Z$（$d\ge1$）と書け、$p\in d\Z$ より $d$ は $1$ か $p$ である。$p\nmid u$ なら $u\in d\Z$ より $d\ne p$、従って $d=1$ であり、$uu'+pp'=1$ となる整数 $u',p'$ がある。すると $v=uvu'+pvp'$ は $p$ で割り切れる）。従って $a\equiv1$ または $a\equiv-1\pmod p$ であり、$p$ が奇数なので $\bar1\ne\overline{-1}$ である。

残りの $p-3$ 個の元は $x\ne x^{-1}$ を満たし、組 $\{x,x^{-1}\}$ に分かれる。各組の積は $\bar1$ である。可換群なので全元の積は順序によらず、
$$
\overline{(p-1)!}=\prod_{x\in(\Z/p\Z)^\times}x=\bar1\cdot\overline{-1}\cdot\bar1^{(p-3)/2}=\overline{-1}
$$
である。例えば $p=7$ では $6!=720=7\cdot103-1$ である。
:::
::::

::::exercise{level=advanced}
$H\le G$ の指数 $n=[G:H]$ が有限なら、$H$ に含まれ指数が $n^n$ 以下の正規部分群が存在することを示せ。

:::solution
**共役は高々 $n$ 個**：$g'=gh$（$h\in H$）なら $g'Hg'^{-1}=g(hHh^{-1})g^{-1}=gHg^{-1}$ なので、共役部分群 $gHg^{-1}$ は左剰余類 $gH$ だけで決まる。従って相異なる共役部分群は高々 $n$ 個であり、それらを $H_1,\ldots,H_m$（$m\le n$）とする。

**共役の指数は $n$**：$g\in G$ を固定すると、$xH\mapsto(gxg^{-1})(gHg^{-1})=gxHg^{-1}$ は $G/H$ から $G/(gHg^{-1})$ への全単射である。実際、$xH=yH\iff x^{-1}y\in H\iff (gxg^{-1})^{-1}(gyg^{-1})=gx^{-1}yg^{-1}\in gHg^{-1}$ なので代表元によらず単射であり、$G/(gHg^{-1})$ の元 $z(gHg^{-1})$ は $x=g^{-1}zg$ の像なので全射である。よって $[G:H_i]=n$ である。

**共通部分**：:ref[prop-grp-index-intersection] の 2 を繰り返し用いると $[G:H_1\cap\cdots\cap H_m]\le n^m\le n^n$ である。この共通部分はすべての共役 $gHg^{-1}$ の共通部分、すなわち正規核 $\operatorname{Core}_G(H)$ であり、:ref[prop-grp-normal-core-closure] により $H$ に含まれる正規部分群である。

なお [群の作用](/math/group-theory/group-actions/) を用いると、正規核の指数は $n!$ を割ることが示される。
:::
::::

::::exercise{level=advanced}
$G$ を有限 Abel 群、$p$ を $|G|$ を割る素数とする。$G$ に位数 $p$ の元が存在することを、$|G|$ に関する帰納法と剰余群を用いて示せ。

:::solution
$|G|$ に関する帰納法で示す。$p\mid|G|$ より $G\ne\{e\}$ なので、$a\ne e$ をとり $m=\operatorname{ord}(a)\ge2$ とする。

$p\mid m$ の場合、:ref[thm-grp-cyclic-order] より $a^{m/p}$ の位数は $m/\gcd(m,m/p)=p$ である。

$p\nmid m$ の場合、$G$ は可換なので $N=\langle a\rangle$ は正規であり、$G/N$ は位数 $|G|/m$ の Abel 群である（:ref[prop-grp-quotient-basic] の 3）。$p\mid|G|=m\cdot|G/N|$ かつ $p\nmid m$ なので $p\mid|G/N|$ である（素数が積を割れば一方の因子を割る。Wilson の定理の演習の解答を参照）。また $m\ge2$ より $|G/N|<|G|$ なので、帰納法の仮定により位数 $p$ の元 $bN\in G/N$ が存在する。$t=\operatorname{ord}(b)$ とおくと $(bN)^t=N$ なので、:ref[prop-grp-quotient-basic] の 2 より $p\mid t$ である。従って $b^{t/p}$ の位数は $p$ である。

（非可換群でも同じ結論が成り立つ。これが Cauchy の定理であり、[Sylow の定理](/math/group-theory/sylow/) で証明する。）
:::
::::

::::exercise{level=advanced}
$n\ge3$ とする。二面体群 $D_n$ の中心 $Z(D_n)$ を求めよ。さらに $n$ が偶数で $n\ge6$ のとき、$D_n/Z(D_n)\cong D_{n/2}$ であることを示せ。

:::solution
**中心**：回転 $r^k$ はすべての回転と可換であり、$sr^k=r^{-k}s$ なので、$r^k$ が $s$ と可換であることは $r^{-k}=r^k$、すなわち $n\mid2k$ と同値である。鏡映 $r^is$ は $r$ と可換でない。実際 $r\cdot r^is=r^{i+1}s$、$r^is\cdot r=r^{i-1}s$ であり、$r^2\ne e$ より両者は異なる。$D_n=\langle r,s\rangle$ なので、中心の元は $r$ と $s$ の両方と可換な元である。従って
$$
Z(D_n)=\begin{cases}\{e\}&(n\text{ が奇数}),\\ \{e,\ r^{n/2}\}&(n\text{ が偶数})\end{cases}
$$
である。

**剰余群**：$n=2m$、$m\ge3$ とし、$Z=Z(D_n)=\{e,r^m\}$、$\bar r=rZ$、$\bar s=sZ$ とおく。$r^k\in Z$ は $m\mid k$ と同値なので $\bar r$ の位数は $m$ であり、$s\notin Z$、$s^2=e$ より $\bar s$ の位数は $2$ である。$D_n$ の元 $r^is^j$ の類は $\bar r^i\bar s^j$ であり、$\bar r^m=\bar e$ なので $0\le i<m$、$0\le j<2$ の $2m$ 個の記号 $\bar r^i\bar s^j$ が $D_n/Z$ のすべての元を表す。$|D_n/Z|=2n/2=2m$ なので、これらは相異なる。さらに $D_n$ の積の規則 :ref[prop-grp-dihedral-form] を $\pi_Z$ で写すと
$$
(\bar r^i\bar s^j)(\bar r^k\bar s^l)=\bar r^{\,i+(-1)^jk}\,\bar s^{\,j+l}
$$
であり、これは $D_m$ の積の規則と同じ形である（$\bar r$ の指数は法 $m$ で考えればよい）。従って $\bar r^i\bar s^j\mapsto r^is^j\in D_m$ は積を保つ全単射であり、$D_n/Z(D_n)\cong D_m$ である。$n=4$ の場合は :ref[exm-grp-quotient-examples] の 7 で見たように $V_4$ になる。
:::
::::

## 確認テスト

::::quiz
$H\le G$、$a,b\in G$ とする。左剰余類について $aH=bH$ が成り立つための必要十分条件はどれか。

- [x] $a^{-1}b\in H$
- [ ] $ab\in H$
- [ ] $ab^{-1}\in H$
- [ ] $a\in H$ かつ $b\in H$

:::explanation
$aH=bH$ は $b\in aH$、すなわち $b=ah$（$h\in H$）と書けることと同値であり、これは $a^{-1}b\in H$ と同値である（:ref[prop-grp-coset-test]）。$ab^{-1}\in H$ は右剰余類についての条件 $Ha=Hb$ であり、非正規な部分群では一般に異なる。$a,b$ は $H$ の外の元でもよいので、最後の選択肢は必要条件ではない。$ab\in H$ は例えば $a=b\notin H$、$a^2\notin H$ のとき成り立たないのに $aH=bH$ となる。
:::
::::

::::quiz{answer="120" tol="0"}
$\SL_2(\F_5)$ の位数を求めよ。

:::explanation
$|\GL_2(\F_5)|=(5^2-1)(5^2-5)=24\cdot20=480$ である（:ref[prop-grp-gl-order]）。$\SL_2(\F_5)$ の左剰余類は行列式の値 $c\in\F_5^\times$ ごとの行列全体であり（:ref[exm-grp-cosets-linear]）、$[\GL_2(\F_5):\SL_2(\F_5)]=|\F_5^\times|=4$ である。Lagrange の定理より $|\SL_2(\F_5)|=480/4=120$ である。
:::
::::

::::quiz{answer="2" tol="0"}
剰余群 $Q_8/\{\pm1\}$ における $i$ の類 $\{\pm i\}$ の位数を求めよ。

:::explanation
$i\notin\{\pm1\}$ なので類は単位元でなく、$i^2=-1\in\{\pm1\}$ なので $(i\{\pm1\})^2=i^2\{\pm1\}$ は単位元である。:ref[prop-grp-quotient-basic] の 2 により位数は $2$ である。$i$ 自身の位数 $4$ の約数になっていることにも注意する。
:::
::::

::::quiz
群 $G$ の部分群 $H$ について、常に $H\trianglelefteq G$ となる条件はどれか。

- [ ] $|H|=2$
- [x] $[G:H]=2$
- [ ] $H$ が巡回群である
- [ ] $[G:H]=3$

:::explanation
指数 $2$ なら、$H$ 以外の左剰余類と右剰余類はどちらも補集合 $G\setminus H$ なので $gH=Hg$ が成り立つ（:ref[prop-grp-index-two]）。$S_3$ の $\langle(1\,2)\rangle$ は位数 $2$ の巡回群で指数 $3$ だが、$(1\,3)$ による共役で $\langle(2\,3)\rangle$ に移るので正規でない。従ってほかの三つの条件はいずれも十分でない。
:::
::::

::::quiz
任意の群 $G$ とその部分群 $H$ について常に成り立つものをすべて選べ。

- [x] $gH\mapsto Hg^{-1}$ は左剰余類全体から右剰余類全体への全単射である
- [ ] $gH$ と $Hg$ が交われば $gH=Hg$ である
- [x] $gH$ が単位元を含むなら $gH=H$ である
- [x] $H$ の共役部分群 $gHg^{-1}$ は $H$ と同じ位数をもつ

:::explanation
$(gH)^{-1}=Hg^{-1}$ なので一つ目の写像は代表元によらず、逆元をもう一度とる写像が逆写像である（:ref[prop-grp-left-right-cosets]）。$e\in gH$ なら $g^{-1}\in H$、従って $g\in H$ であり $gH=H$ である。$h\mapsto ghg^{-1}$ は全単射なので $|gHg^{-1}|=|H|$ である。一方 $gH$ と $Hg$ は常に $g$ を共有するが、$S_3$ の $H=\langle(1\,2)\rangle$、$g=(1\,2\,3)$ では $gH=\{(1\,2\,3),(1\,3)\}\ne Hg=\{(1\,2\,3),(2\,3)\}$ である。「交われば一致する」のは同じ側の剰余類どうしに限る。
:::
::::

::::quiz
交代群 $A_4$（位数 $12$）の部分群の位数として現れるものをすべて選べ。

- [x] $2$
- [x] $3$
- [x] $4$
- [ ] $6$

:::explanation
$\langle(1\,2)(3\,4)\rangle$ は位数 $2$、$\langle(1\,2\,3)\rangle$ は位数 $3$、$V_4$ は位数 $4$ の部分群である。位数 $6$ の部分群は指数 $2$ になるので、:ref[lem-grp-index-two-squares] によりすべての平方を含む。$3$ 巡回置換 $c$ は $c=(c^2)^2$ と平方なので、$8$ 個の $3$ 巡回置換と単位元がすべて含まれ、位数が $6$ を超えて矛盾する（:ref[exm-grp-a4-no-six]）。これは Lagrange の定理の逆が成り立たない例である。
:::
::::

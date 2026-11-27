---
title: 可解群と冪零群
description: 組成列と Jordan–Hölder の定理・Schreier の細分定理、交換子と可換化、導来列による可解群、中心列による冪零群と有限冪零群の構造を扱う。
notation:
  - symbol: 'G=G_0\trianglerighteq G_1\trianglerighteq\cdots\trianglerighteq G_r=1'
    meaning: 劣正規列（各項が直前の項の正規部分群である部分群の列）
  - symbol: '\ell(G)'
    meaning: 群 G の組成列の長さ（組成長）
  - symbol: '[x,y]'
    meaning: 交換子 xyx^{-1}y^{-1}
  - symbol: '[H,K]'
    meaning: H の元と K の元の交換子全体が生成する部分群
  - symbol: "G'=[G,G]"
    meaning: G の交換子群（導来群）
  - symbol: 'G_{\mathrm{ab}}'
    meaning: G の可換化 G/G'
  - symbol: 'G^{(i)}'
    meaning: 第 i 導来群
  - symbol: '\gamma_i(G)'
    meaning: 下中心列の第 i 項
  - symbol: 'Z_i(G)'
    meaning: 上中心列の第 i 項
  - symbol: '\mathrm{UT}_n(K)'
    meaning: 体 K 上の n 次単位上三角行列群
---

## この章の目標

- 組成列を構成し、Jordan–Hölder の定理と Schreier の細分定理を証明して、組成因子が群の不変量であることを理解する。
- 交換子の計算規則を使って交換子群と可換化を求め、Hall–Witt の恒等式から三部分群補題を導く。
- 導来列によって可解群を定義し、部分群・剰余群・拡大で閉じること、有限可解群の組成因子が素数位数の巡回群であることを示す。
- 下中心列・上中心列によって冪零群を定義し、有限 $p$ 群が冪零であること、有限冪零群が Sylow 部分群の直積であることを証明する。

[準同型と同型定理](/math/group-theory/homomorphisms/) で見たように、正規部分群 $N\trianglelefteq G$ があると、群 $G$ は「部分 $N$」と「商 $G/N$」の二つの小さい群に分かれる。$N$ や $G/N$ にさらに正規部分群があれば、同じ操作を繰り返せる。有限群ならこの分解はいつか止まり、それ以上分解できない群、すなわち単純群（:ref[def-grp-simple]）が残る。整数を素因数に分解する操作と似ており、本章の最初の主題は「分解の仕方によらず、最後に残る単純群の組は一定である」という Jordan–Hölder の定理である。

ただし整数と違い、因子の組から元の群は復元できない。例えば $C_6$ と $S_3$ はどちらも $C_2$ と $C_3$ に分解される。そこで、因子の種類によって群の大まかな性質を分類する。すべての因子を Abel 群にできる群が**可解群**であり、方程式が根号で解けるかどうかという Galois 理論の問いから生まれた概念である。さらに強く、各段が「中心に入る」ように分解できる群が**冪零群**であり、有限 $p$ 群がその典型例である。

本章でも置換の積 $\sigma\tau$ は右の $\tau$ を先に施す。群の単位元を $e$、自明群を $1$ と書く。

## 組成列と Jordan–Hölder の定理

### 劣正規列と組成列

:::definition[劣正規列と組成列]{#def-grp-composition}
群 $G$ の部分群の有限列

$$
G=G_0\trianglerighteq G_1\trianglerighteq\cdots\trianglerighteq G_r=1
$$

で、各 $i$ について $G_{i+1}$ が $G_i$ の正規部分群であるものを :term[劣正規列]{en="subnormal series" yomi="れつせいきれつ"} という。$G_{i+1}$ が $G$ 自身の正規部分群である必要はない。すべての $G_i$ が $G$ の正規部分群であるとき :term[正規列]{en="normal series" yomi="せいきれつ"} という。剰余群 $G_i/G_{i+1}$（$0\le i<r$）を列の**因子**、$r$ を列の**長さ**という。

すべての因子が単純群である劣正規列を :term[組成列]{en="composition series" yomi="そせいれつ"} といい、その因子を :term[組成因子]{en="composition factor" yomi="そせいいんし"} という。
:::

単純群は非自明と約束したので、組成列では $G_{i+1}\ne G_i$ であり、包含はすべて真である。文献によっては劣正規列を「正規列」と呼ぶので、読むときは各項が直前の項で正規なのか全体で正規なのかを確認する必要がある。

:::remark[正規性は推移的でない]
劣正規列と正規列を区別するのは、正規部分群の正規部分群が正規とは限らないからである。$A_4$ の中で $V_4=\{e,(12)(34),(13)(24),(14)(23)\}$ は正規であり（共役は巡回型を保つ）、$V_4$ は Abel 群なので $K=\langle(12)(34)\rangle$ は $V_4$ で正規である。しかし

$$
(123)\,(12)(34)\,(123)^{-1}=(23)(14)\notin K
$$

であり、$K$ は $A_4$ で正規でない。したがって $A_4\trianglerighteq V_4\trianglerighteq K\trianglerighteq1$ は劣正規列だが正規列ではない。ここで共役の計算には、$\tau(a_1\,\ldots\,a_k)\tau^{-1}=(\tau(a_1)\,\ldots\,\tau(a_k))$ を用いた。
:::

列を比べるために、項の挿入と因子の一致を定義する。

:::definition[細分と同値]{#def-grp-refinement}
劣正規列 $(G_i)_{0\le i\le r}$ に部分群を挿入して得られる劣正規列を、元の列の :term[細分]{en="refinement" yomi="さいぶん"} という。挿入した項が隣の項と等しくてもよい（すなわち項の重複を許す）。

二つの劣正規列 $(G_i)_{0\le i\le r}$ と $(H_j)_{0\le j\le s}$ が**同値**であるとは、$r=s$ であり、$\{0,\ldots,r-1\}$ の全単射 $\pi$ で、すべての $i$ について

$$
G_i/G_{i+1}\cong H_{\pi(i)}/H_{\pi(i)+1}
$$

となるものが存在することをいう。すなわち、因子を並べ替えれば同型なものが重複度込みで一対一に対応する。
:::

ある段 $G_{i+1}\trianglelefteq G_i$ の間に真に挿入できる項は、商 $G_i/G_{i+1}$ の正規部分群に対応する。

:::lemma[単純な商と極大正規部分群]{#lem-grp-simple-quotient}
$N\trianglelefteq G$ とする。次は同値である。

1. $G/N$ は単純群である。
2. $N\ne G$ であり、$N\subsetneq M\subsetneq G$ を満たす正規部分群 $M\trianglelefteq G$ は存在しない。

2 を満たす $N$ を $G$ の**極大正規部分群**という。
:::

:::proof
対応定理 :ref[thm-grp-correspondence] により、$G/N$ の正規部分群は $N$ を含む $G$ の正規部分群と一対一に対応する。両端はそれぞれ $1\leftrightarrow N$、$G/N\leftrightarrow G$ である。端が異なり、中間に項がない条件がちょうど 1 と 2 である。
:::

この補題から、劣正規列が組成列であることは「各段が真の包含で、どの段の間にも真に新しい項を挿入できない」ことと同値である。実際、$G_{i+1}\trianglelefteq K\trianglelefteq G_i$ となる $K$ を挿入できることは、$K/G_{i+1}$ が $G_i/G_{i+1}$ の正規部分群であることと同値だからである。特に、組成列の細分は項を重複させることしかできない。

:::proposition[有限群の組成列の存在]{#prop-grp-composition-exists}
有限群は組成列をもつ。
:::

:::proof
$|G|=1$ なら長さ $0$ の列をとる。それ以外では真の正規部分群のうち位数最大の $N$ を選ぶ。候補に $1$ があり有限個なので選択できる。最大性と :ref[lem-grp-simple-quotient] より $G/N$ は単純である。$|N|<|G|$ だから位数についての帰納法で $N$ の組成列を得る。その先頭に $G$ を置けばよい。
:::

:::example[組成列をもたない群]{#exm-grp-z-no-composition}
$\Z$ は組成列をもたない。仮に $\Z=G_0\trianglerighteq\cdots\trianglerighteq G_r=0$ が組成列なら、$\Z\ne0$ より $r\ge1$ であり、最後の非零項 $G_{r-1}\cong G_{r-1}/G_r$ は単純である。一方 :ref[thm-grp-cyclic-subgroups] より $G_{r-1}=m\Z$（$m\ge1$）であり、$2m\Z$ は $m\Z$ の非自明な真部分群で、可換性から正規である。これは単純性に反する。

より一般に、組成列をもつ Abel 群は有限である。実際、Abel 群の部分群はすべて正規なので、単純な Abel 群は非自明な真部分群をもたない。非単位元 $a$ をとると $\langle a\rangle$ が全体であり、$\Z$ と同型なら $\langle a^2\rangle$ が非自明な真部分群になるので、有限巡回群である。位数 $n$ が合成数 $n=dm$（$1<d<n$）なら $\langle a^d\rangle$ は位数 $m$ の非自明な真部分群である。よって単純 Abel 群は素数位数の巡回群であり、組成列をもつ Abel 群の位数は組成因子の位数の積として有限である。
:::

:::example[巡回群の組成列]{#exm-grp-cyclic-composition}
$C_{12}=\Z/12\Z$ において、位数 $12,6,3,1$ の部分群からなる列

$$
\Z/12\Z\trianglerighteq\langle\bar2\rangle\trianglerighteq\langle\bar4\rangle\trianglerighteq0
$$

の因子の位数は $2,2,3$ であり、素数位数の群は巡回群（:ref[cor-grp-order-divides]）なので因子は $C_2,C_2,C_3$ である。一方、位数 $12,4,2,1$ の部分群からなる列 $\Z/12\Z\trianglerighteq\langle\bar3\rangle\trianglerighteq\langle\bar6\rangle\trianglerighteq0$ の因子は $C_3,C_2,C_2$ である。両者は中間の部分群も因子の順序も違うが、因子は重複度込みで一致する。

一般に $n=p_1p_2\cdots p_r$（$p_i$ は素数、重複を許す）と書けば、:ref[thm-grp-cyclic-subgroups] により $C_n$ には位数 $p_{i+1}\cdots p_r$ の部分群 $G_i$ がただ一つあり、$G_{i+1}$ は $G_i$ の（唯一の）位数 $p_{i+2}\cdots p_r$ の部分群であるから $G_{i+1}\le G_i$ である。因子 $G_i/G_{i+1}$ は位数 $p_{i+1}$ の巡回群で単純である。したがって $C_n$ の組成因子は $C_{p_1},\ldots,C_{p_r}$ である。
:::

:::example[対称群の組成列]{#exm-grp-symmetric-composition}
1. $S_3\trianglerighteq A_3\trianglerighteq1$ の因子は $C_2,C_3$ である。$C_6\trianglerighteq C_3\trianglerighteq1$ の因子も $C_2,C_3$ であり、組成因子が同じでも群は同型とは限らない。同様に $C_4$ と $V_4=C_2\times C_2$ の組成因子はどちらも $C_2,C_2$ である。
2. $S_4\trianglerighteq A_4\trianglerighteq V_4\trianglerighteq\langle(12)(34)\rangle\trianglerighteq1$ の因子は位数から順に $C_2,C_3,C_2,C_2$ である。$A_4\trianglelefteq S_4$ は指数 $2$ であることから（:ref[prop-grp-index-two]）、$V_4\trianglelefteq A_4$ は共役が巡回型を保つことから、$\langle(12)(34)\rangle\trianglelefteq V_4$ は可換性から従う。
3. $n\ge5$ なら $S_n\trianglerighteq A_n\trianglerighteq1$ は組成列である。$S_n/A_n\cong C_2$ であり、$A_n$ は :ref[thm-grp-alternating-simple] により単純だからである。組成因子は $C_2$ と $A_n$ であり、非可換な単純群が現れる。
:::

### Jordan–Hölder の定理

組成列の一意性を示すために、まず正規部分群と剰余群に組成列を「制限」「射影」する方法を準備する。

:::lemma[正規部分群と剰余群の組成列]{#lem-grp-composition-normal}
$G=G_0\trianglerighteq G_1\trianglerighteq\cdots\trianglerighteq G_r=1$ を組成列、$N\trianglelefteq G$ とする。$N_i=N\cap G_i$、$Q_i=G_iN/N\le G/N$ とおく。

1. $N_{i+1}\trianglelefteq N_i$ であり、$N_i/N_{i+1}$ は自明群であるか $G_i/G_{i+1}$ と同型である。
2. $Q_{i+1}\trianglelefteq Q_i$ であり、$Q_i/Q_{i+1}$ は自明群であるか $G_i/G_{i+1}$ と同型である。
3. 各 $i$ について、$N_i/N_{i+1}$ と $Q_i/Q_{i+1}$ のちょうど一方だけが非自明である。

したがって $(N_i)$ と $(Q_i)$ からそれぞれ重複する項を除けば $N$ と $G/N$ の組成列が得られ、両者の長さの和は $r$ であり、両者の組成因子を合わせたものは $(G_i)$ の組成因子と重複度込みで一致する。
:::

:::proof
$S_i=G_i/G_{i+1}$ とする。$G_i$ の元を $G_iN/G_{i+1}N$ に送る自然な射影を考える。後者の商ができるのは、自然射影 $\pi:G\to G/N$ の下で $\pi(G_{i+1})\trianglelefteq\pi(G_i)$ で、その逆像でも正規性が成り立つからである。$\pi(G_i)=Q_i$ なので、合成 $G_i\to Q_i/Q_{i+1}$ は全射であり、$G_{i+1}$ を核に含む。従って全射 $S_i\to Q_i/Q_{i+1}$ を誘導する。

この全射の核は $(G_i\cap G_{i+1}N)/G_{i+1}$ である。$x=hn\in G_i$（$h\in G_{i+1},\ n\in N$）なら $n=h^{-1}x\in N_i$ であり、逆向きの包含も成り立つから $G_i\cap G_{i+1}N=G_{i+1}N_i$ となる。一方 $N_i\to S_i$ の核は $N_i\cap G_{i+1}=N_{i+1}$ で、像はちょうどこの核 $G_{i+1}N_i/G_{i+1}$ である。これで $N_{i+1}\trianglelefteq N_i$ と、その因子が $S_i$ の正規部分群に同型であることが分かる。

$S_i$ は非自明な単純群なので、この正規部分群は $1$ と $S_i$ のいずれかである。前者なら $Q_i/Q_{i+1}\cong S_i$、後者なら $Q_i=Q_{i+1}$ である。従って各 $i$ の因子 $S_i$ は制限した列と射影した列のちょうど一方に現れる。両列で重複する項を除けば主張を得る。
:::

:::theorem[Jordan–Hölder の定理]{#thm-grp-jordan-holder}
群 $G$ が組成列をもつならば、$G$ の任意の二つの組成列は同値である。すなわち長さが等しく、組成因子は順序と同型を除いて重複度込みで一致する。
:::

:::proof
長さ $r$ の組成列を一つもつ任意の群について、他のすべての組成列がそれと同値になることを $r$ に関する強い帰納法で示す。$r=0$ なら自明群であり、真の包含だけからなる組成列は長さ $0$ に限られる。

$r>0$ とし、二列の最初の真の項を $M,N$ とする。どちらも :ref[lem-grp-simple-quotient] により $G$ の極大正規部分群である。$M=N$ なら残りの二列を $M$ で比較する。一方は長さ $r-1$ なので帰納法を使え、共通の因子 $G/M$ を戻せばよい。

$M\ne N$ なら $MN$ は $M$ を真に含む正規部分群だから $MN=G$ である。実際 $N\le M$ であれば $N$ の極大性から $N=M$ となる。$D=M\cap N$ とおくと第二同型定理により $M/D\cong G/N$、$N/D\cong G/M$ はともに単純である。ここで $D$ が組成列をもつことも確かめる。$D\trianglelefteq M$ なので $M$ にある長さ $r-1$ の組成列に :ref[lem-grp-composition-normal] を適用すれば、$D$ の組成列が得られる。その長さを $t$ とする。

この同じ列の先頭に $M$ または $N$ を付けると、各商が単純なので、いずれも長さ $t+1$ の組成列となる。まず $M$ で帰納法を適用すれば $t+1=r-1$ と分かる。従って $N$ にも長さ $r-1$ の組成列があり、ここでも帰納法で他の列と比較できる。元の二列の因子は、$D$ の共通の列の因子に、それぞれ $G/M,M/D$ と $G/N,N/D$ を加えたものである。先に得た二つの同型でこれらは順を交換して対応し、重複度も長さも一致する。
:::

この定理により、組成列をもつ群 $G$ に対して、組成列の長さ $\ell(G)$（**組成長**）と組成因子の多重集合は $G$ だけで定まる。:ref[lem-grp-composition-normal] と合わせると次を得る。

:::corollary[組成長の加法性]{#cor-grp-length-additive}
$G$ が組成列をもち $N\trianglelefteq G$ ならば、$N$ と $G/N$ も組成列をもち

$$
\ell(G)=\ell(N)+\ell(G/N)
$$

である。さらに $G$ の組成因子は、$N$ の組成因子と $G/N$ の組成因子を合わせたものに重複度込みで一致する。
:::

:::proof
$G$ の組成列を一つとり :ref[lem-grp-composition-normal] を適用すれば、$N$ と $G/N$ の組成列で、長さの和が $\ell(G)$、因子を合わせると $G$ の組成因子になるものが得られる。Jordan–Hölder の定理により、$\ell(N)$, $\ell(G/N)$ と組成因子は組成列の選び方によらないから、主張が従う。
:::

整数の素因数分解の一意性も、Jordan–Hölder の定理の特別な場合として得られる。

:::corollary[素因数分解の一意性]{#cor-grp-unique-factorization}
$2$ 以上の整数 $n$ は素数の積に表され、その表し方は因数の順序を除いて一意である。
:::

:::proof
存在は $n$ に関する帰納法で示す。$n$ の $1$ より大きい約数のうち最小のもの $p$ は素数である（$p$ の $1<d<p$ なる約数は $n$ の約数になり最小性に反する）。$n=p$ なら終わり、そうでなければ $n/p<n$ に帰納法の仮定を用いる。

一意性を示す。$n=p_1\cdots p_r=q_1\cdots q_s$ を二つの素因数分解とする。:ref[exm-grp-cyclic-composition] により、$C_n$ には組成因子が $C_{p_1},\ldots,C_{p_r}$ である組成列と $C_{q_1},\ldots,C_{q_s}$ である組成列がある。Jordan–Hölder の定理より $r=s$ であり、並べ替えれば $C_{p_i}\cong C_{q_{\pi(i)}}$ である。同型な有限群は位数が等しいので $p_i=q_{\pi(i)}$ である。
:::

この応用では巡回群の部分群を指数の整除性から求めた。:ref[thm-grp-cyclic-subgroups] の証明には素因数分解の一意性を必要としない。

### 細分と Schreier の定理

Jordan–Hölder の定理では因子が単純群であった。因子がまだ分解できる二つの劣正規列を比較するには、両方に項を挿入して同じ細かさにそろえればよい。そのための鍵が、二段ずつの列を交差させる次の補題である。

:::lemma[Zassenhaus の補題]{#lem-grp-butterfly}
群 $G$ の部分群 $A'\trianglelefteq A$ と $B'\trianglelefteq B$ をとる。$D=A\cap B$, $E=(A'\cap B)(A\cap B')$ とおく。このとき $A'(A\cap B')\trianglelefteq A'(A\cap B)$, $B'(A'\cap B)\trianglelefteq B'(A\cap B)$, $E\trianglelefteq D$ であり、

$$
\frac{A'(A\cap B)}{A'(A\cap B')}\ \cong\ \frac{D}{E}\ \cong\ \frac{B'(A\cap B)}{B'(A'\cap B)} . \label{eq-grp-zassenhaus}
$$

この補題は、四つの群の包含関係を図に描くと蝶の形になることから**蝶の補題**とも呼ばれる。
:::

:::proof
左側の同型を示す。右側は $A$ と $B$ の役割を入れ替えた同じ議論で得られる。

1. **$E\trianglelefteq D$。** $A'\cap B=A'\cap D$ である（$A'\subseteq A$ より $A'\cap B=A'\cap A\cap B$）。$d\in D\subseteq A$ に対して $dA'd^{-1}=A'$, $dDd^{-1}=D$ だから $d(A'\cap D)d^{-1}=A'\cap D$ である。よって $A'\cap B\trianglelefteq D$ であり、同様に $A\cap B'=B'\cap D\trianglelefteq D$ である。正規部分群の積は正規部分群なので $E\trianglelefteq D$ である。
2. **$A'D$ は部分群。** $D\subseteq A$ は $A'$ を正規化するので $A'D=DA'$ であり、$A'D$ は部分群である。$A'D=A'(A\cap B)$ に注意する。
3. **写像の定義。** $\varphi:A'D\to D/E$ を $\varphi(ad)=dE$（$a\in A'$, $d\in D$）で定める。表し方によらないことを確かめる。$a_1d_1=a_2d_2$ なら $a_2^{-1}a_1=d_2d_1^{-1}\in A'\cap D\subseteq E$ であり、$E\trianglelefteq D$ なので $d_2E=d_1E$ である。
4. **準同型。** $a_1,a_2\in A'$, $d_1,d_2\in D$ に対して

$$
(a_1d_1)(a_2d_2)=a_1(d_1a_2d_1^{-1})\,d_1d_2,\qquad d_1a_2d_1^{-1}\in A'
$$

だから $\varphi((a_1d_1)(a_2d_2))=d_1d_2E=\varphi(a_1d_1)\varphi(a_2d_2)$ である。$\varphi$ は明らかに全射である。

5. **核。** $\varphi(ad)=E$ と $d\in E$ は同値なので $\Ker\varphi=A'E$ である。$A'\cap B\subseteq A'$ だから $A'E=A'(A'\cap B)(A\cap B')=A'(A\cap B')$ である。

:ref[thm-grp-first-iso] により $A'(A\cap B')=\Ker\varphi$ は $A'D$ の正規部分群で、$A'D/A'(A\cap B')\cong D/E$ である。
:::

:::theorem[Schreier の細分定理]{#thm-grp-schreier-refinement}
群 $G$ の任意の二つの劣正規列は、互いに同値な細分をもつ。
:::

:::proof
片方の列の長さが $0$ なら $G=1$ であり、他方と同じ長さまで $1$ を重ねて細分すればよい。以下、長さがどちらも正である二つの劣正規列を

$$
G=A_0\trianglerighteq A_1\trianglerighteq\cdots\trianglerighteq A_r=1,\qquad
G=B_0\trianglerighteq B_1\trianglerighteq\cdots\trianglerighteq B_s=1
$$

とする。$0\le i<r$, $0\le j\le s$ に対して

$$
A_{i,j}=A_{i+1}(A_i\cap B_j),\qquad
B_{j,i}=B_{j+1}(B_j\cap A_i)\quad(0\le j<s,\ 0\le i\le r)
$$

とおく。

1. **端点。** $B_0=G$ より $A_{i,0}=A_{i+1}A_i=A_i$、$B_s=1$ より $A_{i,s}=A_{i+1}$ である。同様に $B_{j,0}=B_j$, $B_{j,r}=B_{j+1}$ である。
2. **正規性。** :ref[lem-grp-butterfly] を $A=A_i$, $A'=A_{i+1}$, $B=B_j$, $B'=B_{j+1}$ に適用すると、$A_{i,j+1}=A'(A\cap B')\trianglelefteq A'(A\cap B)=A_{i,j}$ である。したがって

$$
A_i=A_{i,0}\trianglerighteq A_{i,1}\trianglerighteq\cdots\trianglerighteq A_{i,s}=A_{i+1}
$$

であり、これを各 $i$ について並べた列は $(A_i)$ の長さ $rs$ の細分である。同様に $B_{j,0}\trianglerighteq\cdots\trianglerighteq B_{j,r}$ を並べた列は $(B_j)$ の長さ $rs$ の細分である。

3. **因子の対応。** 同じ補題により

$$
A_{i,j}/A_{i,j+1}\cong (A\cap B)/E\cong B_{j,i}/B_{j,i+1}
$$

である。第一の細分の因子は添字の対 $(i,j)$ で、第二の細分の因子は対 $(j,i)$ で番号付けられており、$(i,j)\mapsto(j,i)$ は添字集合の全単射である。よって二つの細分は同値である。
:::

細分では項の重複を許したが、重複した項（自明な因子）を除いても同値性は保たれる。上の対応で自明な因子は自明な因子に対応するからである。

:::proof[Jordan–Hölder の定理の別証明]
組成列 $(G_i)$ と $(H_j)$ に Schreier の細分定理を適用し、同値な細分を得る。組成列の細分は項を重複させることしかできない（:ref[lem-grp-simple-quotient] の直後の注意）。したがって細分の非自明な因子は元の組成列の因子そのものである。同値な細分から自明な因子を除くと、$(G_i)$ の因子と $(H_j)$ の因子が重複度込みで一対一に対応する。
:::

:::example[無限群での細分]{#exm-grp-schreier-z}
組成列をもたない $\Z$ でも Schreier の定理は使える。$A$: $\Z\trianglerighteq2\Z\trianglerighteq0$ と $B$: $\Z\trianglerighteq3\Z\trianglerighteq0$ に上の構成を適用する（加法群なので積は和である）。

$$
\begin{aligned}
A_{0,j}&=2\Z+(\Z\cap B_j):\quad \Z,\ \Z,\ 2\Z,\\
A_{1,j}&=0+(2\Z\cap B_j):\quad 2\Z,\ 6\Z,\ 0
\end{aligned}
$$

であり、$A$ の細分は $\Z\trianglerighteq\Z\trianglerighteq2\Z\trianglerighteq6\Z\trianglerighteq0$、因子は順に $1,\ C_2,\ C_3,\ \Z$ である。同様に $B$ の細分は $\Z\trianglerighteq\Z\trianglerighteq3\Z\trianglerighteq6\Z\trianglerighteq0$、因子は $1,\ C_3,\ C_2,\ \Z$ である。対応 $(i,j)\leftrightarrow(j,i)$ は $(0,1)\leftrightarrow(1,0)$ で $C_2$ どうし、$(1,0)\leftrightarrow(0,1)$ で $C_3$ どうし、$(1,1)\leftrightarrow(1,1)$ で $6\Z\cong\Z$ どうしを結ぶ。
:::

### 主組成列と単純群の分類

正規列（全項が $G$ で正規）だけを考えても同様の理論が成り立つ。$G$ の正規列で、$G$ の正規部分群をそれ以上真に挿入できないものを :term[主組成列]{en="chief series" yomi="しゅそせいれつ"} といい、その因子を**主因子**という。$A_i$, $B_j$ がすべて $G$ で正規なら、Schreier の証明に現れる $A_{i,j}=A_{i+1}(A_i\cap B_j)$ も $G$ で正規である（正規部分群の共通部分と積は正規）。主組成列の正規な細分は項を重複させることしかできないので、Jordan–Hölder の定理の別証明と同じ議論により、二つの主組成列は同値である。有限群 $G$ が主組成列をもつことは、:ref[prop-grp-composition-exists] と同様に示せる。$G_0=G$ とし、$G_i\ne1$ が定まったら、$G_i$ に真に含まれる $G$ の正規部分群（$1$ がその一つである）のうち位数最大のものを $G_{i+1}$ とする。位数は真に減るので有限回で $1$ に達し、最大性から各段の間に $G$ の正規部分群を真に挿入することはできない。

主因子は単純とは限らない。例えば $S_4\trianglerighteq A_4\trianglerighteq V_4\trianglerighteq1$ は主組成列である。実際 $V_4$ の位数 $2$ の部分群は三つとも $S_4$ で共役なので $S_4$ の正規部分群ではなく、$V_4$ と $1$ の間には $S_4$ の正規部分群を挿入できない。主因子は $C_2,C_3,V_4$ である。

:::remark[Hölder のプログラム]
Jordan–Hölder の定理により、有限群の研究は (1) 有限単純群をすべて決定すること、(2) 与えられた単純群を組成因子にもつ群をすべて決定すること（拡大問題）の二つに分けられる。この方針を Hölder のプログラムという。(1) は多数の数学者による膨大な研究の結果、有限単純群の分類定理として完成した：有限単純群は、素数位数の巡回群、$5$ 次以上の交代群、Lie 型の群の系列、および $26$ 個の散在型単純群のいずれかに同型である :cite[dummit-foote]。その証明は本サイトの範囲を大きく超える。(2) は一般には非常に難しく、次節以降の可解群・冪零群は、因子に条件を課して扱いやすくした群の族である。
:::

## 交換子群

群を Abel 群へ写す準同型を考えると、$xy$ と $yx$ の像は必ず一致する。その差 $xy(yx)^{-1}=xyx^{-1}y^{-1}$ を集めたものが、Abel 群へ移るときに必ず失われる情報である。逆にこの差だけを消せば、情報をできるだけ残したまま Abel 群が得られる。

### 交換子と計算規則

:::definition[交換子と交換子群]{#def-grp-commutator}
群 $G$ の元 $x,y$ に対して

$$
[x,y]=xyx^{-1}y^{-1}
$$

を $x$ と $y$ の :term[交換子]{en="commutator" yomi="こうかんし"} という。部分群 $H,K\le G$ に対して、$h\in H$, $k\in K$ の交換子全体が生成する部分群を

$$
[H,K]=\langle\,[h,k]\mid h\in H,\ k\in K\,\rangle
$$

と書く。特に $G'=[G,G]$ を $G$ の :term[交換子群]{en="commutator subgroup" yomi="こうかんしぐん"}（導来群）という。
:::

交換子の集合そのものが積で閉じることは定義からは保証されないので、生成する部分群をとる。$[H,K]$ の元は、交換子とその逆元の有限積である。

:::remark[交換子の記法の違い]
$[x,y]=x^{-1}y^{-1}xy$ と定義する文献も多い（例えば :cite[dummit-foote]）。この元は本章の記法では $[x^{-1},y^{-1}]$ に等しいので、$x,y$ が群全体（または部分群全体）を動くとき、生成される部分群 $[H,K]$ や $G'$ はどちらの定義でも同じである。ただし下に述べる個々の恒等式は記法によって形が変わるので、他書の公式を引用するときは注意を要する。
:::

:::lemma[交換子の計算規則]{#lem-grp-commutator-identities}
群 $G$ の任意の元 $x,y,z,g$ について次が成り立つ。

1. $[x,y]=e$ と $xy=yx$ は同値である。また $[x,y]^{-1}=[y,x]$ であり、$xy=[x,y]\,yx$ である。
2. 任意の準同型 $f:G\to H$ について $f([x,y])=[f(x),f(y)]$ である。特に $g[x,y]g^{-1}=[gxg^{-1},gyg^{-1}]$ である。
3. $[xy,z]=x[y,z]x^{-1}\,[x,z]$ である。
4. $[x,yz]=[x,y]\,y[x,z]y^{-1}$ である。
5. （**Hall–Witt の恒等式**）

$$
y\bigl[[y^{-1},x],z\bigr]y^{-1}\cdot z\bigl[[z^{-1},y],x\bigr]z^{-1}\cdot x\bigl[[x^{-1},z],y\bigr]x^{-1}=e . \label{eq-grp-hall-witt}
$$
:::

:::proof
1. $[x,y]=e$ の両辺に右から $yx$ を掛ければ $xy=yx$ を得る。逆も同様である。$[x,y]^{-1}=(xyx^{-1}y^{-1})^{-1}=yxy^{-1}x^{-1}=[y,x]$ であり、$[x,y]\,yx=xyx^{-1}y^{-1}yx=xy$ である。

2. $f$ は積と逆元を保つので $f(xyx^{-1}y^{-1})=f(x)f(y)f(x)^{-1}f(y)^{-1}$ である。共役 $c_g(x)=gxg^{-1}$ は準同型なので後半を得る。

3. 右辺を展開すると

$$
x(yzy^{-1}z^{-1})x^{-1}\,(xzx^{-1}z^{-1})=xyzy^{-1}(z^{-1}z)x^{-1}z^{-1}=(xy)z(xy)^{-1}z^{-1}
$$

となり、左辺に等しい。

4. 同様に

$$
(xyx^{-1}y^{-1})\,y(xzx^{-1}z^{-1})y^{-1}=xy(x^{-1}x)zx^{-1}z^{-1}y^{-1}=x(yz)x^{-1}(yz)^{-1}
$$

である。

5. $u(x,y,z)=xyx^{-1}zx$ とおく。左辺の第一因子を展開すると、$[y^{-1},x]=y^{-1}xyx^{-1}$, $[y^{-1},x]^{-1}=xy^{-1}x^{-1}y$ より

$$
y\bigl[[y^{-1},x],z\bigr]y^{-1}
=y\,(y^{-1}xyx^{-1})\,z\,(xy^{-1}x^{-1}y)\,z^{-1}y^{-1}
=(xyx^{-1}zx)\,(y^{-1}x^{-1}yz^{-1}y^{-1})
$$

である。右側の括弧は $u(y,z,x)^{-1}=(yzy^{-1}xy)^{-1}$ に等しいので、第一因子は $u(x,y,z)\,u(y,z,x)^{-1}$ である。第二・第三因子はこの式で $(x,y,z)$ を $(y,z,x)$, $(z,x,y)$ に置き換えたものなので、それぞれ $u(y,z,x)u(z,x,y)^{-1}$, $u(z,x,y)u(x,y,z)^{-1}$ である。三つを掛けると中央で打ち消し合い、

$$
u(x,y,z)u(y,z,x)^{-1}\cdot u(y,z,x)u(z,x,y)^{-1}\cdot u(z,x,y)u(x,y,z)^{-1}=e
$$

を得る。
:::

3 と 4 は「交換子の双線形性」が共役の分だけずれることを表している。通常の積のように $[xy,z]=[x,z][y,z]$ と展開するのは誤りであり、一般には成り立たない。ただし交換子がすべて中心に入る群（後で見る冪零類 $2$ 以下の群）では共役が外れて、双線形性が成り立つ（演習問題）。

:::proposition[交換子部分群の基本性質]{#prop-grp-commutator-subgroups}
$H,K,H_1,K_1\le G$ とする。

1. $[H,K]=[K,H]$ である。また $H\le H_1$, $K\le K_1$ なら $[H,K]\le[H_1,K_1]$ である。
2. 準同型 $f:G\to G_1$ に対して $f([H,K])=[f(H),f(K)]$ である。
3. $[H,K]=1$ と、$H$ のすべての元が $K$ のすべての元と可換であることは同値である。
4. $K$ が $H$ を正規化する（すべての $k\in K$ で $kHk^{-1}=H$）なら $[H,K]\le H$ である。
5. $H,K\trianglelefteq G$ なら $[H,K]\trianglelefteq G$ かつ $[H,K]\le H\cap K$ である。$H,K$ がともに $G$ の特性部分群なら $[H,K]$ も特性部分群である。
6. $N\trianglelefteq G$ とし、$\pi:G\to G/N$ を射影とする。$[H,K]\le N$ と、$\pi(H)$ の元と $\pi(K)$ の元がすべて可換であることは同値である。
:::

:::proof
1. :ref[lem-grp-commutator-identities] の 1 により、$[H,K]$ の生成元 $[h,k]$ の逆元は $[K,H]$ の生成元 $[k,h]$ である。部分群は生成元の逆元で生成されても同じなので $[H,K]=[K,H]$ である。後半は生成系の包含から従う。

2. 準同型の像について $f(\langle S\rangle)=\langle f(S)\rangle$ である（$\langle S\rangle$ の元は $S$ の元とその逆元の積であり、$f$ はそれを $f(S)$ の元とその逆元の積に写す）。$S=\{[h,k]\}$ に適用し、$f([h,k])=[f(h),f(k)]$ を使えばよい。

3. :ref[lem-grp-commutator-identities] の 1 による。

4. $[h,k]=h\,(kh^{-1}k^{-1})$ であり、$kh^{-1}k^{-1}\in H$ だから $[h,k]\in H$ である。

5. $g\in G$ による共役 $c_g$ は自己同型なので、2 より $c_g([H,K])=[c_g(H),c_g(K)]=[H,K]$ である。特性部分群の場合も、任意の $\phi\in\Aut(G)$ について同様に $\phi([H,K])=[\phi(H),\phi(K)]=[H,K]$ である。$K$ は $H$ を正規化するので 4 より $[H,K]\le H$、同様に 1 と 4 から $[H,K]=[K,H]\le K$ である。

6. $\pi([h,k])=[\pi(h),\pi(k)]$ なので、すべての生成元が $N=\Ker\pi$ に入ることと、$\pi(h)$ と $\pi(k)$ がすべて可換であることは同値である。$N$ は部分群なので、生成元がすべて $N$ に入ることと $[H,K]\le N$ は同値である。
:::

### 可換化

:::proposition[交換子群と可換化]{#prop-grp-abelianization}
群 $G$ について次が成り立つ。

1. $G'$ は $G$ の特性部分群であり、特に正規部分群である。
2. $N\trianglelefteq G$ に対して、$G/N$ が Abel 群であることと $G'\le N$ は同値である。したがって $G'$ は、剰余群が Abel 群になる正規部分群のうち最小のものである。
3. $G_{\mathrm{ab}}=G/G'$ は Abel 群である。これを $G$ の :term[可換化]{en="abelianization" yomi="かかんか"} という。射影を $\pi:G\to G_{\mathrm{ab}}$ とすると、任意の Abel 群 $A$ と準同型 $f:G\to A$ に対し、$f=\bar f\circ\pi$ を満たす準同型 $\bar f:G_{\mathrm{ab}}\to A$ がただ一つ存在する。
:::

:::proof
1. :ref[prop-grp-commutator-subgroups] の 5 を $H=K=G$ に適用すればよい。

2. :ref[prop-grp-commutator-subgroups] の 6 を $H=K=G$ に適用すると、$G'\le N$ と「$G/N$ の任意の二元が可換」が同値である。

3. 2 で $N=G'$ とすれば $G/G'$ は Abel 群である。$f:G\to A$ が Abel 群への準同型なら $f([x,y])=[f(x),f(y)]=e$ なので、すべての交換子が $\Ker f$ に含まれ、$G'\le\Ker f$ である。:ref[prop-grp-quotient-universal] により $f=\bar f\circ\pi$ となる準同型 $\bar f$ がただ一つ存在する。
:::

3 は、$G$ から Abel 群 $A$ への準同型と $G_{\mathrm{ab}}$ から $A$ への準同型が $f\leftrightarrow\bar f$ で一対一に対応することを意味する。例えば $G$ の一次元表現（$\C^\times$ への準同型）はすべて $G_{\mathrm{ab}}$ を経由する。

交換子群を具体的に求めるときは、生成元どうしの交換子から出発すればよい。部分集合 $S\subseteq G$ を含む最小の正規部分群を $S$ の**正規閉包**という。これは $S$ の元の共役 $gsg^{-1}$（$g\in G$, $s\in S$）全体が生成する部分群に等しい。実際、その部分群は共役で閉じた集合から生成されるので、:ref[lem-grp-commutator-identities] の 2 と同様に共役で保たれて正規であり、$S$ を含む任意の正規部分群は各 $gsg^{-1}$ を含むからである。

:::proposition[生成元による交換子群の計算]{#prop-grp-commutator-generators}
$G=\langle S\rangle$ ならば、$G'$ は $\{[s,t]\mid s,t\in S\}$ の正規閉包に等しい。
:::

:::proof
この正規閉包を $N$ とする。$G'$ は各 $[s,t]$ を含む正規部分群なので $N\le G'$ である。

逆を示す。$\bar G=G/N$ とし、$s\in S$ の像を $\bar s$ と書く。$[s,t]\in N$ なので、任意の $s,t\in S$ について $\bar s\bar t=\bar t\bar s$ である。$s\in S$ を固定すると、$\bar s$ の中心化群 $C_{\bar G}(\bar s)$ はすべての $\bar t$ を含む部分群であり、$\bar G$ は $\{\bar t\}$ で生成されるから $C_{\bar G}(\bar s)=\bar G$ である。すなわちすべての $\bar s$ は中心 $Z(\bar G)$ に属する。中心は部分群で、生成系を含むから $Z(\bar G)=\bar G$ であり、$\bar G$ は Abel 群である。:ref[prop-grp-abelianization] の 2 より $G'\le N$ である。
:::

:::example[二面体群の交換子群]{#exm-grp-dihedral-commutator}
$D_n=\langle r,s\rangle$（$n\ge3$）では :ref[prop-grp-dihedral-form] より $srs^{-1}=r^{-1}$ だから

$$
[r,s]=r\,(sr^{-1}s^{-1})=r\cdot r=r^2,\qquad [s,r]=r^{-2},\qquad [r,r]=[s,s]=e
$$

である。$\langle r^2\rangle$ は $r$ と可換で、$s r^2s^{-1}=r^{-2}$ だから正規である。よって :ref[prop-grp-commutator-generators] より $D_n'=\langle r^2\rangle$ である。

- $n$ が奇数なら $\gcd(2,n)=1$ なので $\langle r^2\rangle=\langle r\rangle$（:ref[thm-grp-cyclic-order]）であり、$D_n'$ の位数は $n$、可換化は位数 $2$ の群 $C_2$ である。
- $n$ が偶数なら $\langle r^2\rangle$ の位数は $n/2$ で、可換化は位数 $4$ である。その元 $\bar r,\bar s$ は $\bar r^2=\bar s^2=e$ を満たし可換なので、すべての元の平方が単位元である。位数 $4$ の元がないので巡回群でなく、$(D_n)_{\mathrm{ab}}\cong C_2\times C_2$ である。
:::

:::example[四元数群の交換子群]{#exm-grp-quaternion-commutator}
$Q_8=\{\pm1,\pm i,\pm j,\pm k\}$（:ref[exm-grp-quaternion]）では $i^{-1}=-i$, $j^{-1}=-j$ より

$$
[i,j]=ij(-i)(-j)=(ij)^2=k^2=-1
$$

である。$\{\pm1\}$ は中心なので正規であり、$Q_8=\langle i,j\rangle$ に :ref[prop-grp-commutator-generators] を適用すると $Q_8'=\{\pm1\}$ である。可換化 $Q_8/\{\pm1\}$ は位数 $4$ で、$\bar i^2=\bar j^2=\overline{-1}=e$ だから $C_2\times C_2$ である。
:::

:::example[対称群と交代群の交換子群]{#exm-grp-symmetric-commutator}
1. **$S_n'=A_n$。** $n\le2$ なら両辺は自明群である。$n\ge2$ では符号 $\sgn:S_n\to\{\pm1\}$ が Abel 群への準同型で核が $A_n$ だから（:ref[def-grp-alternating]）、:ref[prop-grp-abelianization] より $S_n'\le A_n$ である。$n\ge3$ なら

$$
[(12),(23)]=\bigl((12)(23)\bigr)^2=(123)^2=(132)
$$

である（互換は自身の逆元である）。:ref[thm-grp-cycle-classes] よりすべての $3$ 巡回は $S_n$ で $(132)$ と共役であり、$S_n'$ は正規だからすべての $3$ 巡回を含む。:ref[prop-grp-three-generate] により $A_n\le S_n'$ である。特に $(S_n)_{\mathrm{ab}}\cong C_2$（$n\ge2$）である。

2. **$A_4'=V_4$。** $V_4\trianglelefteq A_4$ で $A_4/V_4$ は位数 $3$ の巡回群なので $A_4'\le V_4$ である。$A_4$ は可換でない（例えば $(123)(12)(34)=(134)$, $(12)(34)(123)=(243)$）から $A_4'\ne1$ であり、$A_4'$ は $V_4$ の非単位元、すなわち二重互換を一つ含む。$(123)$ による共役は三つの二重互換を巡回的に入れ替える（$(12)(34)\mapsto(23)(14)\mapsto(31)(24)$）ので、正規部分群 $A_4'$ は三つすべてを含み、$A_4'=V_4$ である。
3. **$n\ge5$ なら $A_n'=A_n$。** $A_n'$ は $A_n$ の正規部分群で、$A_n$ は非可換だから非自明である。:ref[thm-grp-alternating-simple] により $A_n'=A_n$ である。具体的にも

$$
[(124),(135)]=(123)
$$

が成り立つ（右から順に $1\mapsto5\mapsto5\mapsto1\mapsto2$ などと追えば確かめられる）。
:::

$G'=G$ を満たす群を :term[完全群]{en="perfect group" yomi="かんぜんぐん"} という。非可換単純群は完全群である（$G'$ は非自明な正規部分群だから）。$A_5$ はその最小の例である。

:::example[Heisenberg 群]{#exm-grp-heisenberg}
体 $K$ に対し

$$
u(a,b,c)=\begin{pmatrix}1&a&c\\0&1&b\\0&0&1\end{pmatrix}\qquad(a,b,c\in K)
$$

の全体を $U$ とする。行列の積を計算すると

$$
\begin{aligned}
u(a,b,c)\,u(x,y,z)&=u(a+x,\ b+y,\ c+z+ay),\\
u(a,b,c)^{-1}&=u(-a,\ -b,\ ab-c)
\end{aligned} \label{eq-grp-heisenberg-product}
$$

である（第二式は第一式に代入して $u(0,0,c+(ab-c)+a(-b))=u(0,0,0)$ となることで確かめられる）。よって $U$ は $\GL_3(K)$ の部分群であり、:term[Heisenberg 群]{en="Heisenberg group" yomi="はいぜんべるくぐん"} と呼ばれる。

1. **中心。** $u(a,b,c)$ と $u(x,y,z)$ が可換であることは、積の第三成分を比べて $ay=xb$ と同値である。これがすべての $x,y$ で成り立つのは $a=b=0$ のときに限るので、$Z(U)=\{u(0,0,c)\mid c\in K\}\cong(K,+)$ である。
2. **交換子。** $[u(a,b,c),u(x,y,z)]=u(0,0,ay-bx)$ である。実際、$g=u(a,b,c)$, $h=u(x,y,z)$ とすると $gh$ と $hg$ は第一・第二成分が等しく第三成分の差が $ay-xb$ なので、$gh=u(0,0,ay-xb)\,hg$ である（$u(0,0,\cdot)$ は中心にあり、積の公式で第三成分を加えるだけである）。$[g,h]=gh(hg)^{-1}$ より主張を得る。
3. **交換子群。** 2 より $U'\le Z(U)$ である。$a=t$, $y=1$, $b=x=0$ とすれば $u(0,0,t)$ が交換子として現れるので $U'=Z(U)$ である。可換化は $u(a,b,c)\mapsto(a,b)$ により $(K^2,+)$ と同型である。この写像は積の公式から準同型で、核が $Z(U)=U'$ だからである。

$K=\F_p$ なら $|U|=p^3$, $|Z(U)|=p$ である。
:::

### 三部分群補題

Hall–Witt の恒等式は、三つの部分群の交換子についての次の補題を導く。後で冪零群の下中心列を調べるときに用いる。

:::lemma[三部分群補題]{#lem-grp-three-subgroups}
$H,K,L\le G$, $N\trianglelefteq G$ とする。$[[H,K],L]\le N$ かつ $[[K,L],H]\le N$ ならば $[[L,H],K]\le N$ である。
:::

:::proof
まず $N=1$ の場合を示す。仮定は、$[H,K]$ の元が $L$ の元とすべて可換であり、$[K,L]$ の元が $H$ の元とすべて可換であることを意味する（:ref[prop-grp-commutator-subgroups] の 3）。

$l\in L$, $h\in H$, $k\in K$ を任意にとり、Hall–Witt の恒等式 :eqref[eq-grp-hall-witt] で $x=h$, $y=l^{-1}$, $z=k$ とおく。第二因子の内側 $[[k^{-1},l^{-1}],h]$ は $[k^{-1},l^{-1}]\in[K,L]$ と $h\in H$ の交換子なので $e$ であり、第三因子の内側 $[[h^{-1},k],l^{-1}]$ は $[h^{-1},k]\in[H,K]$ と $l^{-1}\in L$ の交換子なので $e$ である。よって第一因子 $l^{-1}[[l,h],k]\,l$ も $e$ であり、$[[l,h],k]=e$ を得る。すなわち各 $[l,h]$ は $K$ のすべての元と可換である。$K$ の中心化群 $C_G(K)=\{g\mid gk=kg\ (k\in K)\}$ は部分群で、生成元 $[l,h]$ をすべて含むので $[L,H]\le C_G(K)$ であり、$[[L,H],K]=1$ である。

一般の場合は射影 $\pi:G\to G/N$ を用いる。:ref[prop-grp-commutator-subgroups] の 2 により $[[\pi(H),\pi(K)],\pi(L)]=\pi([[H,K],L])=1$ であり、同様に $[[\pi(K),\pi(L)],\pi(H)]=1$ である。$N=1$ の場合を $G/N$ に適用すると $\pi([[L,H],K])=[[\pi(L),\pi(H)],\pi(K)]=1$、すなわち $[[L,H],K]\le N$ である。
:::

## 可解群

可換化は群から非可換な部分 $G'$ を一度に取り除く。取り除いた $G'$ にも同じ操作を施し、これを繰り返したときに有限回で自明群に達する群が可解群である。

### 導来列

:::definition[導来列と可解群]{#def-grp-solvable}
群 $G$ に対し $G^{(0)}=G$, $G^{(i+1)}=\bigl(G^{(i)}\bigr)'=[G^{(i)},G^{(i)}]$ と定める。列

$$
G=G^{(0)}\trianglerighteq G^{(1)}\trianglerighteq G^{(2)}\trianglerighteq\cdots
$$

を $G$ の :term[導来列]{en="derived series" yomi="どうらいれつ"} という。ある $d\ge0$ について $G^{(d)}=1$ となるとき、$G$ を :term[可解群]{en="solvable group" yomi="かかいぐん"} といい、そのような最小の $d$ を $G$ の :term[導来長]{en="derived length" yomi="どうらいちょう"} という。
:::

導来長 $0$ の群は自明群、導来長 $1$ 以下の群は Abel 群（$G'=1$）である。:ref[prop-grp-abelianization] により $G^{(i+1)}$ は $G^{(i)}$ の特性部分群で、$G^{(i)}/G^{(i+1)}$ は Abel 群である。さらに各 $G^{(i)}$ は $G$ の特性部分群である。実際 $i$ に関する帰納法で、$G^{(i)}$ が $G$ の特性部分群なら $\phi\in\Aut(G)$ の制限は $G^{(i)}$ の自己同型であり、:ref[prop-grp-commutator-subgroups] の 2 より $\phi(G^{(i+1)})=[\phi(G^{(i)}),\phi(G^{(i)})]=G^{(i+1)}$ となる。特に導来列は正規列である。

:::lemma[導来列と準同型]{#lem-grp-derived-hom}
1. $H\le G$ なら、すべての $i$ について $H^{(i)}\le G^{(i)}$ である。
2. 全射準同型 $f:G\to Q$ について、すべての $i$ で $f(G^{(i)})=Q^{(i)}$ である。
:::

:::proof
いずれも $i$ に関する帰納法である。$i=0$ では明らかである。

1. $H^{(i)}\le G^{(i)}$ なら、:ref[prop-grp-commutator-subgroups] の 1 より $H^{(i+1)}=[H^{(i)},H^{(i)}]\le[G^{(i)},G^{(i)}]=G^{(i+1)}$ である。
2. $f(G^{(i)})=Q^{(i)}$ なら、:ref[prop-grp-commutator-subgroups] の 2 より $f(G^{(i+1)})=[f(G^{(i)}),f(G^{(i)})]=[Q^{(i)},Q^{(i)}]=Q^{(i+1)}$ である。
:::

:::theorem[可解群の特徴付けと閉性]{#thm-grp-solvable-extension}
1. $G$ が可解であることと、因子がすべて Abel 群である劣正規列

$$
G=H_0\trianglerighteq H_1\trianglerighteq\cdots\trianglerighteq H_d=1
$$

をもつことは同値である。このとき $G^{(i)}\le H_i$ であり、導来長は $d$ 以下である。

2. 可解群の部分群は可解である。
3. 可解群の剰余群（準同型による像）は可解である。
4. $N\trianglelefteq G$ について、$N$ と $G/N$ がともに可解なら $G$ も可解であり、導来長について $\mathrm{dl}(G)\le\mathrm{dl}(N)+\mathrm{dl}(G/N)$ が成り立つ。
5. 有限個の可解群の直積は可解である。
:::

:::proof
1. 可解なら導来列 $G^{(0)}\trianglerighteq\cdots\trianglerighteq G^{(d)}=1$ がそのような列である。逆にそのような列があるとする。$H_i/H_{i+1}$ が Abel 群なので :ref[prop-grp-abelianization] より $H_i'\le H_{i+1}$ である。$G^{(i)}\le H_i$ を帰納法で示す。$i=0$ は $G^{(0)}=G=H_0$ である。$G^{(i)}\le H_i$ なら $G^{(i+1)}=[G^{(i)},G^{(i)}]\le[H_i,H_i]=H_i'\le H_{i+1}$ である。よって $G^{(d)}\le H_d=1$ である。

2. $G^{(d)}=1$, $H\le G$ なら :ref[lem-grp-derived-hom] の 1 より $H^{(d)}\le G^{(d)}=1$ である。

3. 全射 $f:G\to Q$ について、:ref[lem-grp-derived-hom] の 2 より $Q^{(d)}=f(G^{(d)})=1$ である。

4. $\mathrm{dl}(N)=a$, $\mathrm{dl}(G/N)=b$ とし、射影を $\pi:G\to G/N$ とする。:ref[lem-grp-derived-hom] の 2 より $\pi(G^{(b)})=(G/N)^{(b)}=1$、すなわち $G^{(b)}\le N$ である。導来列の定義から $G^{(b+i)}=(G^{(b)})^{(i)}$ なので、:ref[lem-grp-derived-hom] の 1 より $G^{(a+b)}=(G^{(b)})^{(a)}\le N^{(a)}=1$ である。

5. $H,K$ が可解なら $H\times K$ は正規部分群 $H\times1\cong H$ をもち、剰余群は $K$ と同型である。4 より可解である。有限個の場合は帰納法による。
:::

1 の条件は「Abel 群を有限回の拡大で積み重ねて得られる」ことを意味し、可解群の定義としてこちらを採用する文献も多い。導来列は、このような列のうち各項が最も小さいものである。

:::corollary[組成因子による可解性の判定]{#cor-grp-solvable-factors}
有限群 $G$ が可解であるための必要十分条件は、$G$ の組成因子がすべて素数位数の巡回群であることである。特に有限 $p$ 群は可解である。
:::

:::proof
$G$ が可解とし、組成列 $(G_i)$ をとる（:ref[prop-grp-composition-exists]）。組成因子 $G_i/G_{i+1}$ は可解群の部分群 $G_i$ の剰余群なので、定理の 2, 3 により可解である。非自明な可解単純群 $S$ では $S'$ は正規部分群なので $1$ か $S$ であり、$S'=S$ なら導来列が $S$ のまま止まらず可解性に反するから $S'=1$、すなわち $S$ は Abel 群である。:ref[exm-grp-z-no-composition] で見たように、単純 Abel 群は素数位数の巡回群である。

逆に組成因子がすべて素数位数の巡回群なら、組成列は因子が Abel 群の劣正規列であり、定理の 1 より $G$ は可解である。

位数 $p^a$ の群の組成因子の位数は $p^a$ を割る（:ref[thm-grp-lagrange] を各段に用いる）ので、上の議論の前半を使わずとも次のように示せる。$a$ に関する帰納法で、$a=0$ は自明である。$a\ge1$ なら :ref[cor-grp-p-center] により中心に位数 $p$ の部分群 $C$ があり、$C$ は正規で Abel 群、$G/C$ は位数 $p^{a-1}$ の群なので帰納法の仮定より可解である。定理の 4 より $G$ は可解である。
:::

### 可解群の例

:::example[小さい群の導来列]{#exm-grp-derived-examples}
1. **二面体群。** :ref[exm-grp-dihedral-commutator] より $D_n'=\langle r^2\rangle$ は巡回群なので $D_n''=1$ であり、$D_n$（$n\ge3$）は導来長 $2$ の可解群である（$D_n$ は非可換なので導来長は $1$ ではない）。
2. **四元数群。** $Q_8'=\{\pm1\}$, $Q_8''=1$ なので導来長は $2$ である。
3. **$S_4$。** :ref[exm-grp-symmetric-commutator] より $S_4'=A_4$, $S_4''=A_4'=V_4$, $S_4'''=V_4'=1$ である。導来列 $S_4\trianglerighteq A_4\trianglerighteq V_4\trianglerighteq1$ の因子は $C_2,\ C_3,\ C_2\times C_2$ であり、導来長は $3$ である。導来列の因子は Abel 群だが単純とは限らないことに注意する。
4. **Heisenberg 群。** :ref[exm-grp-heisenberg] より $U'=Z(U)$ は Abel 群なので $U''=1$ であり、導来長は $2$ である。
:::

:::proposition[対称群の可解性]{#prop-grp-sn-solvable}
$S_n$ が可解であるための必要十分条件は $n\le4$ である。
:::

:::proof
$n\le2$ なら $S_n$ は Abel 群である。$S_3\trianglerighteq A_3\trianglerighteq1$ は因子 $C_2,C_3$ の劣正規列であり、$S_4$ は上の例で見たように導来長 $3$ の可解群である。

$n\ge5$ とする。:ref[exm-grp-symmetric-commutator] より $S_n'=A_n$, $A_n'=A_n$ なので、すべての $i\ge1$ で $S_n^{(i)}=A_n\ne1$ であり、$S_n$ は可解でない。同じことは、可解群の部分群 $A_n$ が非可換単純群になることはない（:ref[cor-grp-solvable-factors] の証明の前半）ことからも分かる。
:::

有限群の可解性を位数だけから判定できる場合がある。

:::proposition[位数 pq と p²q の群]{#prop-grp-pq-solvable}
$p,q$ を素数とする。位数が $pq$ または $p^2q$ の群は可解である。
:::

:::proof
$p=q$ なら有限 $p$ 群の可解性を使う。$p\ne q$ のときは、一つの Sylow 部分群が正規なら、それ自身も商も素数冪位数の可解群になるため、拡大の可解性から結論が出る。位数 $pq$ で $p<q$ とすれば $n_q\mid p,\ n_q\equiv1\bmod q$ より $n_q=1$ である。

位数 $p^2q$ の場合に両 Sylow 部分群が唯一でないと仮定する。$n_p\mid q$ より $n_p=q\equiv1\bmod p$、従って $q>p$ である。$n_q$ は $p,p^2$ のいずれかだが $p<q$ なので $p$ は除かれ、$q\mid p^2-1$ となる。$q$ は $p-1$ を割れず、$p+1$ を割るから $q=p+1$ であり、連続する素数 $p=2,\ q=3$ だけが残る。このとき四個の位数 $3$ の部分群が単位元を除いて互いに交わらず、八個の非単位元を使う。位数 $4$ の Sylow 部分群は残り四元をちょうど占めるので唯一となり矛盾する。
:::

素因数が二種類の場合の一般論は、表現論を用いて初めて証明された深い定理である。

:::theorem[Burnside の pᵃqᵇ 定理]{#thm-grp-burnside-paqb}
$p,q$ を素数、$a,b\ge0$ を整数とする。位数 $p^aq^b$ の群は可解である。
:::

標準的な証明は複素指標の直交関係と代数的整数の性質を必要とするため、ここでは :cite[dummit-foote] に譲る。これらの道具は[表現論の群論への応用](/math/representation-theory/applications-groups/)で扱う。奇数位数の有限群がすべて可解であるという Feit–Thompson の定理も知られている。その証明は有限群の詳細な局所解析を必要とする。:cite[dummit-foote]

:::example[位数 60 未満の群はすべて可解]{#exm-grp-order-less-60}
三種類以上の素数で割れる $60$ 未満の正整数は $30,42$ だけである。実際、$2,3$ のどちらかを欠けば最小でも $2\cdot5\cdot7=70$、両方を含めば第三の素数は $5$ または $7$ で、もう一因子を掛ければ $60$ 以上になる。他の位数なら :ref[thm-grp-burnside-paqb] を適用できる。位数 $30$ の場合は :ref[prop-grp-order-30] により位数 $15$ の巡回正規部分群があり、商は位数 $2$ なので可解である。位数 $42$ の場合は $n_7\mid6,\ n_7\equiv1\bmod7$ より正規な $C_7$ があり、商の位数は $6$ で可解である。いずれも可解群の拡大だから可解であり、位数 $60$ の非可換単純群 $A_5$ が位数最小の非可解群となる。
:::

### 極小正規部分群と主因子

有限可解群の「最小の構成要素」は特別な形をしている。$G$ の正規部分群 $N\ne1$ で、$1$ と $N$ の間に $G$ の正規部分群を真に挟めないものを :term[極小正規部分群]{en="minimal normal subgroup" yomi="きょくしょうせいきぶぶんぐん"} という。素数 $p$ について、Abel 群 $A$ ですべての元が $x^p=e$ を満たす非自明なものを :term[基本可換 p 群]{en="elementary abelian p-group" yomi="きほんかかんぴーぐん"} という。有限な基本可換 $p$ 群は、整数 $k$ の作用 $x\mapsto x^k$ が $k\bmod p$ だけで決まることから $\F_p$ 上の有限次元ベクトル空間とみなせ、基底をとれば $C_p\times\cdots\times C_p$ と同型である（[ベクトル空間・基底・次元](/math/linear-algebra/vector-spaces/)）。

:::proposition[可解群の極小正規部分群]{#prop-grp-minimal-normal}
有限可解群 $G$ の極小正規部分群は、ある素数 $p$ に対する基本可換 $p$ 群である。
:::

:::proof
$N$ を極小正規部分群とする。

1. **$N$ は Abel 群。** $N'$ は $N$ の特性部分群で $N\trianglelefteq G$ だから、$G$ の元による共役（$N$ の自己同型を引き起こす）で保たれ、$N'\trianglelefteq G$ である（[準同型と同型定理](/math/group-theory/homomorphisms/) の特性部分群の項）。$N$ は可解群の部分群なので可解であり、$N\ne1$ だから $N'\ne N$ である（そうでなければ導来列が止まらない）。極小性から $N'=1$ である。
2. **$N$ は基本可換 $p$ 群。** $|N|$ を割る素数 $p$ を一つとり、$N[p]=\{x\in N\mid x^p=e\}$ とおく。$N$ は Abel 群なので $N[p]$ は部分群であり（$(xy)^p=x^py^p$）、自己同型は位数を保つので $N$ の特性部分群である。よって 1 と同様に $N[p]\trianglelefteq G$ である。:ref[thm-grp-cauchy] より $N$ は位数 $p$ の元をもつので $N[p]\ne1$ であり、極小性から $N[p]=N$ である。
:::

:::corollary[可解群の主因子]{#cor-grp-chief-solvable}
有限可解群の主組成列の因子はすべて基本可換 $p$ 群（$p$ は因子ごとに異なりうる）である。
:::

:::proof
主組成列 $(G_i)$ の因子 $G_i/G_{i+1}$ は、:ref[thm-grp-correspondence] により $G/G_{i+1}$ の極小正規部分群である。$G/G_{i+1}$ は可解群の剰余群なので可解であり、命題を適用すればよい。
:::

例えば $S_4$ の主因子 $C_2,C_3,V_4$ はいずれも基本可換である。

### 方程式の可解性との関係

標数 $0$ の体上では、多項式の根が根号を用いて求められる条件は、分解体の Galois 群が可解であることである。独立な不定元を係数とする一般の $n$ 次方程式の Galois 群は $S_n$ である。従って :ref[prop-grp-sn-solvable] より $n\ge5$ では、係数の四則演算と有限回の冪根で全根を与える一般公式は存在しない。個別の方程式では Galois 群がより小さく、根号で解ける場合がある。体の拡大との対応は[方程式の可解性](/math/fields-galois/solvability/)で扱う。

## 冪零群

導来列は残った部分群どうしの交換子をとる。下中心列では一方を毎回 $G$ 全体に戻すので、交換子を作る相手が増え、列が自明群に達しないこともある。例えば $S_3$ の導来列は $S_3,A_3,1$ と進むが、下中心列は $S_3,A_3,A_3,\ldots$ と続く。このより強い停止条件によって冪零性を定義する。

### 中心列

:::definition[下中心列・上中心列と冪零群]{#def-grp-nilpotent}
群 $G$ に対し、$\gamma_1(G)=G$, $\gamma_{i+1}(G)=[G,\gamma_i(G)]$ と定めた列

$$
G=\gamma_1(G)\ge\gamma_2(G)\ge\gamma_3(G)\ge\cdots
$$

を :term[下中心列]{en="lower central series" yomi="かちゅうしんれつ"} という。また $Z_0(G)=1$ とし、$Z_{i+1}(G)$ を $G/Z_i(G)$ の中心の、射影 $G\to G/Z_i(G)$ による逆像と定めた列

$$
1=Z_0(G)\le Z_1(G)\le Z_2(G)\le\cdots
$$

を :term[上中心列]{en="upper central series" yomi="じょうちゅうしんれつ"} という。$Z_1(G)=Z(G)$ である。ある $c\ge0$ について $\gamma_{c+1}(G)=1$ となるとき、$G$ を :term[冪零群]{en="nilpotent group" yomi="べきれいぐん"} といい、そのような最小の $c$ を $G$ の :term[冪零類]{en="nilpotency class" yomi="べきれいるい"} という。
:::

$\gamma_2(G)=G'$ である。冪零類 $0$ の群は自明群、冪零類 $1$ 以下の群は Abel 群である。定義が意味をもつこと（$\gamma_i$ が減少し、$Z_i$ が正規部分群で増加すること）を確かめておく。

:::lemma[中心列の基本性質]{#lem-grp-central-basic}
1. 各 $\gamma_i(G)$ は $G$ の特性部分群であり、$\gamma_{i+1}(G)\le\gamma_i(G)$ である。
2. 各 $Z_i(G)$ は $G$ の特性部分群であり、$Z_i(G)\le Z_{i+1}(G)$ である。
3. $x\in Z_{i+1}(G)$ であるための必要十分条件は、すべての $g\in G$ について $[g,x]\in Z_i(G)$ となることである。すなわち $[G,Z_{i+1}(G)]\le Z_i(G)$ であり、$Z_{i+1}(G)$ はこの性質をもつ部分群のうち最大のものである。
:::

:::proof
1. $i$ に関する帰納法で示す。$\gamma_1=G$ は特性部分群である。$\gamma_i$ が特性部分群なら、:ref[prop-grp-commutator-subgroups] の 5 により $\gamma_{i+1}=[G,\gamma_i]$ も特性部分群であり、同じく 5 により $\gamma_{i+1}=[G,\gamma_i]\le\gamma_i$ である。

2. $Z_0=1$ は特性部分群である。$Z_i$ が特性部分群なら、$\phi\in\Aut(G)$ は $G/Z_i$ の自己同型 $\bar\phi(gZ_i)=\phi(g)Z_i$ を引き起こす（$\phi(Z_i)=Z_i$ なので矛盾なく定まり、逆写像も同様に作れる）。自己同型は中心を中心に写すので $\bar\phi(Z(G/Z_i))=Z(G/Z_i)$ であり、逆像をとって $\phi(Z_{i+1})=Z_{i+1}$ を得る。$Z_{i+1}$ は単位元 $Z_i$ の逆像 $Z_i$ を含む。

3. $xZ_i\in Z(G/Z_i)$ は、すべての $g$ について $(gZ_i)(xZ_i)=(xZ_i)(gZ_i)$、すなわち $[g,x]Z_i=Z_i$ であることと同値である。後半は、$[G,X]\le Z_i$ を満たす部分群 $X$ の元はすべてこの条件を満たすことから従う。
:::

:::proposition[中心列による特徴付け]{#prop-grp-central-series}
部分群の列 $1=B_0\le B_1\le\cdots\le B_c=G$ で、すべての $1\le i\le c$ について $[G,B_i]\le B_{i-1}$ を満たすものを長さ $c$ の :term[中心列]{en="central series" yomi="ちゅうしんれつ"} という。$c\ge0$ について次は同値である。

1. $\gamma_{c+1}(G)=1$。
2. $Z_c(G)=G$。
3. $G$ は長さ $c$ の中心列をもつ。

さらに、任意の長さ $c$ の中心列 $(B_i)$ について

$$
\gamma_{c+1-i}(G)\le B_i\le Z_i(G)\qquad(0\le i\le c) \label{eq-grp-central-sandwich}
$$

が成り立つ。したがって冪零群の冪零類は、下中心列が $1$ に達するまでの段数とも、上中心列が $G$ に達するまでの段数とも等しい。また中心列の各項は $G$ の正規部分群で、隣り合う商は Abel 群であり、冪零群は可解である。
:::

:::proof
まず 3 を仮定して :eqref[eq-grp-central-sandwich] を示す。

- **$B_i\le Z_i$。** $i$ に関する帰納法で示す。$B_0=1=Z_0$ である。$B_{i-1}\le Z_{i-1}$ なら $[G,B_i]\le B_{i-1}\le Z_{i-1}$ なので、:ref[lem-grp-central-basic] の 3 より $B_i\le Z_i$ である。
- **$\gamma_{c+1-i}\le B_i$。** $i=c,c-1,\ldots,0$ の順に示す。$i=c$ では $\gamma_1=G=B_c$ である。$\gamma_{c+1-i}\le B_i$（$i\ge1$）なら、:ref[prop-grp-commutator-subgroups] の 1 より $\gamma_{c+2-i}=[G,\gamma_{c+1-i}]\le[G,B_i]\le B_{i-1}$ である。

$i=0$ とすれば $\gamma_{c+1}\le B_0=1$、$i=c$ とすれば $G=B_c\le Z_c$ なので、3 から 1 と 2 が従う。

逆に 1 を仮定すると $B_i=\gamma_{c+1-i}$ とおいた列は、$B_0=\gamma_{c+1}=1$, $B_c=G$, $[G,B_i]=\gamma_{c+2-i}=B_{i-1}$ を満たす中心列である（増加列であることは :ref[lem-grp-central-basic] の 1 による）。2 を仮定すると $B_i=Z_i$ とおいた列が、:ref[lem-grp-central-basic] の 3 により中心列である。

最後に、中心列の項 $B_i$ は $g\in G$, $b\in B_i$ について $gbg^{-1}=[g,b]\,b\in B_{i-1}B_i=B_i$ を満たすので正規である。また $[B_i,B_i]\le[G,B_i]\le B_{i-1}$ なので :ref[prop-grp-commutator-subgroups] の 6 により $B_i/B_{i-1}$ は Abel 群である。冪零群は中心列 $G=B_c\trianglerighteq\cdots\trianglerighteq B_0=1$ をもつので、:ref[thm-grp-solvable-extension] の 1 により可解である。
:::

中心列の条件 $[G,B_i]\le B_{i-1}$ は、$B_i/B_{i-1}\le Z(G/B_{i-1})$ と同値である（:ref[prop-grp-commutator-subgroups] の 6）。すなわち冪零群とは「中心による拡大を有限回繰り返して得られる群」である。可解群の定義（Abel 群による拡大の繰り返し）より強い条件であり、各段で $B_i/B_{i-1}$ がその段の商 $G/B_{i-1}$ の中心に入ることを要求する。

### 冪零群の例

:::example[二面体群の冪零性]{#exm-grp-dihedral-nilpotent}
$D_n$（$n\ge3$）が冪零であるための必要十分条件は $n$ が $2$ の冪であることであり、$n=2^k$（$k\ge2$）のとき冪零類は $k$ である。

$i\ge1$ について $\gamma_{i+1}(D_n)=\langle r^{2^i}\rangle$ を帰納法で示す。$i=1$ は :ref[exm-grp-dihedral-commutator] である。$\gamma_{i+1}=\langle r^m\rangle$（$m=2^i$）とすると、$\gamma_{i+2}=[D_n,\langle r^m\rangle]$ は交換子 $[g,r^{mk}]$ で生成される。$g=r^j$ なら交換子は $e$ であり、$g=r^js$ なら $s r^{mk}s^{-1}=r^{-mk}$ より

$$
[r^js,r^{mk}]=r^js\,r^{mk}\,s^{-1}r^{-j}\,r^{-mk}=r^j r^{-mk}r^{-j}r^{-mk}=r^{-2mk}
$$

である。よって $\gamma_{i+2}=\langle r^{2m}\rangle=\langle r^{2^{i+1}}\rangle$ である。

$\langle r^{2^i}\rangle=1$ は $n\mid2^i$ と同値である。$n$ が奇素数で割り切れればこれは決して成り立たず、$D_n$ は冪零でない。$n=2^k$ なら $\gamma_{i+1}=1$ と $i\ge k$ は同値なので冪零類は $k$ である。例えば $D_4$（位数 $8$）の冪零類は $2$、$D_8$（位数 $16$）は $3$ である。$D_3\cong S_3$ では $\gamma_2=\gamma_3=\cdots=A_3$ であり、上中心列も $Z(S_3)=1$ のため $Z_i(S_3)=1$ のまま止まる。
:::

:::example[四元数群と Heisenberg 群]{#exm-grp-quaternion-heisenberg-nilpotent}
$Q_8'=\{\pm1\}$ は中心にあるので $\gamma_3(Q_8)=[Q_8,\{\pm1\}]=1$ であり、$Q_8$ の冪零類は $2$ である。同様に Heisenberg 群 $U$ では $\gamma_2(U)=U'=Z(U)$（:ref[exm-grp-heisenberg]）なので $\gamma_3(U)=1$ であり、冪零類は $2$ である。$K$ が無限体でも同じなので、冪零性は有限性とは無関係な概念である。

$K=\F_p$ のときの $U$ を調べる。:eqref[eq-grp-heisenberg-product] から $m$ に関する帰納法で

$$
u(a,b,c)^m=u\Bigl(ma,\ mb,\ mc+\binom m2ab\Bigr)\qquad(m\ge1)
$$

が成り立つ。実際、右辺に $u(a,b,c)$ を右から掛けると第三成分は $mc+\binom m2ab+c+ma\cdot b$ となり、$\binom m2+m=\binom{m+1}2$ である。$p$ が奇素数なら $\binom p2=p\cdot\frac{p-1}2$ は $p$ の倍数なので、$U$ のすべての元は $u^p=e$ を満たす。一方 $p=2$, $a=b=1$, $c=0$ では $u(1,1,0)^2=u(0,0,1)\ne e$ である。$r=u(1,1,0)$, $s=u(1,0,0)$ とおくと、$r$ の位数は $4$、$s^2=e$ であり、積の公式から

$$
srs=u(0,1,1)\,s=u(1,1,1)=u(-1,-1,1)=r^{-1}
$$

である（$\F_2$ では $-1=1$）。$s\notin\langle r\rangle=\{e,\,r,\,u(0,0,1),\,u(1,1,1)\}$ なので、$\langle r,s\rangle$ は位数 $8$ の群全体であり、:ref[prop-grp-dihedral-form] と同じ関係式を満たす。よって $\mathrm{UT}_3(\F_2)\cong D_4$ である。
:::

Heisenberg 群は次の行列群の特別な場合である。体 $K$ 上の $n$ 次上三角行列で対角成分がすべて $1$ のもの全体を $\mathrm{UT}_n(K)$ と書き、:term[単位上三角行列群]{en="unitriangular group" yomi="たんいうえさんかくぎょうれつぐん"} という。$E_{ij}$ を $(i,j)$ 成分だけが $1$ の行列単位とし、$x_{ij}(t)=I+tE_{ij}$（$i<j$, $t\in K$）とおく。

:::proposition[単位上三角行列群の冪零類]{#prop-grp-unitriangular}
$n\ge2$ のとき $\mathrm{UT}_n(K)$ は群であり、冪零類 $n-1$ の冪零群である。
:::

:::proof
$1\le m\le n$ に対し、$(i,j)$ 成分が $j-i<m$ のとき $0$ である行列全体を $\mathfrak n_m$ とする。$\mathfrak n_1$ は狭義上三角行列全体、$\mathfrak n_n=0$ であり、$U_m=I+\mathfrak n_m$ とおくと $U_1=\mathrm{UT}_n(K)$, $U_n=\{I\}$ である。

1. **$\mathfrak n_a\mathfrak n_b\subseteq\mathfrak n_{a+b}$。** $A\in\mathfrak n_a$, $B\in\mathfrak n_b$ なら $(AB)_{ij}=\sum_kA_{ik}B_{kj}$ の各項は $k-i\ge a$ かつ $j-k\ge b$ のときしか非零でないので、$j-i<a+b$ なら $0$ である（$a+b>n-1$ のときは $AB=0$）。
2. **$U_m$ は部分群。** $A,B\in\mathfrak n_m$ なら $(I+A)(I+B)=I+A+B+AB\in U_m$ である。また 1 より $A^n\in\mathfrak n_{mn}=0$ なので

$$
(I+A)\sum_{k=0}^{n-1}(-A)^k=I-(-A)^n=I
$$

であり、$(I+A)^{-1}=\sum_{k=0}^{n-1}(-A)^k\in U_m$ である。特に $\mathrm{UT}_n(K)=U_1$ は群である。

3. **$[U_1,U_m]\le U_{m+1}$。** $g=I+A\in U_1$, $h=I+B\in U_m$ とする。$gh-hg=AB-BA\in\mathfrak n_{m+1}$ である。一方

$$
[g,h]-I=ghg^{-1}h^{-1}-hgg^{-1}h^{-1}=(gh-hg)\,g^{-1}h^{-1}
$$

であり、$g^{-1}h^{-1}\in I+\mathfrak n_1$ と 1 から $(gh-hg)g^{-1}h^{-1}\in\mathfrak n_{m+1}$ である。よって $[g,h]\in U_{m+1}$ であり、$U_{m+1}$ は部分群なので生成される部分群も $U_{m+1}$ に含まれる。

4. **$\gamma_m\le U_m$。** $\gamma_1=U_1$ であり、$\gamma_m\le U_m$ なら 3 より $\gamma_{m+1}=[U_1,\gamma_m]\le[U_1,U_m]\le U_{m+1}$ である。よって $\gamma_n\le U_n=\{I\}$ であり、冪零類は $n-1$ 以下である。
5. **冪零類はちょうど $n-1$。** $i<j<k$ なら $E_{ij}E_{jk}=E_{ik}$, $E_{jk}E_{ij}=0$, $E_{ik}E_{ij}=E_{ik}E_{jk}=0$ である。$g=x_{ij}(a)$, $h=x_{jk}(b)$ に 3 の式を使うと $gh-hg=abE_{ik}$、$g^{-1}h^{-1}=(I-aE_{ij})(I-bE_{jk})$ であり、

$$
[x_{ij}(a),x_{jk}(b)]=I+abE_{ik}(I-aE_{ij})(I-bE_{jk})=x_{ik}(ab) . \label{eq-grp-unitriangular-commutator}
$$

これを用いて、すべての $i$ と $t\in K$ について $x_{i,i+m}(t)\in\gamma_m$ であることを $m$ に関する帰納法で示す。$m=1$ は明らかである。$m$ で成り立てば $x_{i,i+1+m}(t)=[x_{i,i+1}(1),x_{i+1,i+1+m}(t)]\in[G,\gamma_m]=\gamma_{m+1}$ である。特に $x_{1n}(1)\in\gamma_{n-1}$ は単位行列でないので、$\gamma_{n-1}\ne1$ である。
:::

$n=3$ のとき $\mathrm{UT}_3(K)$ は Heisenberg 群であり、冪零類 $2$ である。$\mathrm{UT}_n(\F_p)$ は対角より上の $n(n-1)/2$ 個の成分を自由に選べるので、位数 $p^{n(n-1)/2}$ の有限 $p$ 群である。一方、正則行列の第 $k$ 列は最初の $k-1$ 列が張る $p^{k-1}$ 元の部分空間の外から選ぶので $|\GL_n(\F_p)|=(p^n-1)(p^n-p)\cdots(p^n-p^{n-1})=p^{n(n-1)/2}\prod_{k=1}^{n}(p^k-1)$ であり、その $p$ 部分はちょうど $p^{n(n-1)/2}$ である。よって $\mathrm{UT}_n(\F_p)$ は $\GL_n(\F_p)$ の Sylow $p$ 部分群である。

### 冪零群の閉性

:::proposition[冪零群の閉性]{#prop-grp-nilpotent-closure}
1. 冪零類 $c$ の群の部分群は、冪零類 $c$ 以下の冪零群である。
2. 冪零類 $c$ の群の剰余群は、冪零類 $c$ 以下の冪零群である。
3. 冪零群 $H,K$ の直積 $H\times K$ は冪零であり、冪零類は両者の冪零類の大きい方である。
4. $N\le Z(G)$ で $G/N$ が冪零類 $d$ の冪零群なら、$G$ は冪零類 $d+1$ 以下の冪零群である。特に $G/Z(G)$ が冪零なら $G$ も冪零である。
5. $G$ が冪零類 $c\ge1$ の冪零群なら $G/Z(G)$ の冪零類は $c-1$ である。
:::

:::proof
1. $H\le G$ について $\gamma_i(H)\le\gamma_i(G)$ を帰納法で示す：$\gamma_{i+1}(H)=[H,\gamma_i(H)]\le[G,\gamma_i(G)]=\gamma_{i+1}(G)$（:ref[prop-grp-commutator-subgroups] の 1）。

2. 全射 $f:G\to Q$ について $f(\gamma_i(G))=\gamma_i(Q)$ を帰納法で示す：$f(\gamma_{i+1}(G))=[f(G),f(\gamma_i(G))]=[Q,\gamma_i(Q)]$（同 2）。

3. 直積では $[(h,k),(h',k')]=([h,h'],[k,k'])$ である。$\gamma_i(H\times K)=\gamma_i(H)\times\gamma_i(K)$ を帰納法で示す。$\subseteq$ は、生成元 $[(h,k),(x,y)]$（$(x,y)\in\gamma_i(H)\times\gamma_i(K)$）が $\gamma_{i+1}(H)\times\gamma_{i+1}(K)$ に入ることから従う。$\supseteq$ は、$h\in H$, $x\in\gamma_i(H)$ について $([h,x],e)=[(h,e),(x,e)]$ であり $(x,e)\in\gamma_i(H\times K)$ だから $\gamma_{i+1}(H)\times1\le\gamma_{i+1}(H\times K)$ となること、および $K$ についての同様の包含から従う。

4. 2 の議論を射影 $\pi:G\to G/N$ に用いると $\pi(\gamma_{d+1}(G))=\gamma_{d+1}(G/N)=1$、すなわち $\gamma_{d+1}(G)\le N\le Z(G)$ である。よって $\gamma_{d+2}(G)=[G,\gamma_{d+1}(G)]\le[G,Z(G)]=1$ である。

5. $\bar G=G/Z(G)$ とし、$Z_i(\bar G)=Z_{i+1}(G)/Z(G)$ を $i$ に関する帰納法で示す。$i=0$ では両辺が自明群である。$i$ で成り立てば、:ref[lem-grp-central-basic] の 3 により

$$
\bar x\in Z_{i+1}(\bar G)\iff[\bar g,\bar x]\in Z_i(\bar G)\ (\forall g)\iff[g,x]\in Z_{i+1}(G)\ (\forall g)\iff x\in Z_{i+2}(G)
$$

である（中央の同値には $Z(G)\le Z_{i+1}(G)$ を用いた）。よって $Z_{c-1}(\bar G)=\bar G$ と $Z_c(G)=G$ は同値であり、:ref[prop-grp-central-series] により主張を得る。
:::

冪零群の拡大は冪零とは限らない。$S_3$ の正規部分群 $A_3\cong C_3$ も剰余群 $S_3/A_3\cong C_2$ も Abel 群（冪零類 $1$）だが、$S_3$ は冪零でない（:ref[exm-grp-dihedral-nilpotent]）。4 のように、核が中心に入れば冪零性を保証できる。これは十分条件である。例えば $D_4$ は冪零群であり、回転群を核とする商は $C_2$ だが、回転群は中心に含まれない。この点で冪零群は可解群（:ref[thm-grp-solvable-extension] の 4）と異なる。

### 有限 p 群

:::lemma[中心による剰余群が巡回なら可換]{#lem-grp-center-cyclic-quotient}
$G/Z(G)$ が巡回群なら $G$ は Abel 群である。
:::

:::proof
$G/Z(G)=\langle aZ(G)\rangle$ とすると、$G$ の元はすべて $a^iz$（$i\in\Z$, $z\in Z(G)$）の形である。$(a^iz)(a^jw)=a^{i+j}zw=(a^jw)(a^iz)$ なので $G$ は Abel 群である。
:::

:::theorem[有限 p 群の冪零性]{#thm-grp-p-nilpotent}
位数 $p^a$（$p$ は素数）の群 $P$ は冪零であり、上中心列は高々 $a$ 段で $P$ に達する。さらに $a\ge2$ なら $P$ の冪零類は $a-1$ 以下である。
:::

:::proof
$Z_i=Z_i(P)\ne P$ ならば、$P/Z_i$ は非自明な $p$ 群なので :ref[cor-grp-p-center] によりその中心の位数は $p$ 以上であり、$|Z_{i+1}|\ge p\,|Z_i|$ である。したがって上中心列は $P$ に達するまで位数が段ごとに $p$ 倍以上になり、高々 $a$ 段で $P$ に達する。:ref[prop-grp-central-series] により $P$ は冪零で、冪零類は $a$ 以下である。

$a\ge2$ とし、冪零類を $c$ とする。$c\le1$ なら $c\le a-1$ である。$c\ge2$ とする。冪零類が $c$ なので、$i<c$ では $Z_i\ne P$ であり、上の議論から $|Z_{c-1}|\ge p^{c-1}$ である。次に $P/Z_{c-1}$ は巡回群でない。実際 $\bar P=P/Z_{c-2}$ とおくと、上中心列の定義から $Z(\bar P)=Z_{c-1}/Z_{c-2}$ であり、:ref[thm-grp-third-iso] より $\bar P/Z(\bar P)\cong P/Z_{c-1}$ である。これが巡回群なら :ref[lem-grp-center-cyclic-quotient] より $\bar P$ は Abel 群で $Z(\bar P)=\bar P$、すなわち $Z_{c-1}=P$ となり、冪零類が $c$ であることに反する。位数 $p$ の群は巡回群なので $|P/Z_{c-1}|\ge p^2$ であり、

$$
p^a=|P|=|Z_{c-1}|\,|P/Z_{c-1}|\ge p^{c-1}\cdot p^2=p^{c+1}
$$

を得る。よって $c\le a-1$ である。
:::

例えば位数 $8$ の非可換群 $D_4$ と $Q_8$ の冪零類はちょうど $2=3-1$ である。位数 $2^{k+1}$ の $D_{2^k}$ の冪零類は $k$（:ref[exm-grp-dihedral-nilpotent]）であり、この上限を達成している。

### 冪零群の性質

冪零群では、下中心列の項どうしの交換子が「次数の和」の項に入る。この性質から、冪零群の導来長が冪零類の対数程度で抑えられる。

:::proposition[下中心列の交換子]{#prop-grp-gamma-commutator}
任意の群 $G$ と $i,j\ge1$ について $[\gamma_i(G),\gamma_j(G)]\le\gamma_{i+j}(G)$ である。特に $G^{(k)}\le\gamma_{2^k}(G)$ であり、冪零類 $c\ge1$ の群の導来長は $\lceil\log_2(c+1)\rceil$ 以下である。
:::

:::proof
$\gamma_i=\gamma_i(G)$ と書き、「すべての $i\ge1$ について $[\gamma_i,\gamma_j]\le\gamma_{i+j}$」を $j$ に関する帰納法で示す。$j=1$ では :ref[prop-grp-commutator-subgroups] の 1 より $[\gamma_i,G]=[G,\gamma_i]=\gamma_{i+1}$ である。

$j$ で成り立つとし、$[\gamma_i,\gamma_{j+1}]=[\gamma_i,[G,\gamma_j]]=[[\gamma_j,G],\gamma_i]$ を調べる（:ref[prop-grp-commutator-subgroups] の 1 を二回使った）。:ref[lem-grp-three-subgroups] を $L=\gamma_j$, $H=G$, $K=\gamma_i$, $N=\gamma_{i+j+1}$（正規部分群）に適用する。仮定の二つを確かめる。

- $[[H,K],L]=[[G,\gamma_i],\gamma_j]=[\gamma_{i+1},\gamma_j]\le\gamma_{i+j+1}$（帰納法の仮定を $i+1$ と $j$ に適用）。
- $[[K,L],H]=[[\gamma_i,\gamma_j],G]\le[\gamma_{i+j},G]=\gamma_{i+j+1}$（帰納法の仮定と単調性）。

よって $[[L,H],K]=[[\gamma_j,G],\gamma_i]\le\gamma_{i+j+1}$ であり、$j+1$ でも成り立つ。

後半は $k$ に関する帰納法である。$G^{(0)}=G=\gamma_1$ であり、$G^{(k)}\le\gamma_{2^k}$ なら $G^{(k+1)}=[G^{(k)},G^{(k)}]\le[\gamma_{2^k},\gamma_{2^k}]\le\gamma_{2^{k+1}}$ である。冪零類が $c$ なら $2^k\ge c+1$ で $G^{(k)}\le\gamma_{2^k}\le\gamma_{c+1}=1$ となる。
:::

例えば冪零類 $2$ または $3$ の群の導来長は $2$ 以下である。$D_8$（冪零類 $3$）の導来列は $D_8\trianglerighteq\langle r^2\rangle\trianglerighteq1$ であり、確かに導来長 $2$ である。

:::proposition[正規化群の増大]{#prop-grp-normalizer-condition}
$G$ を冪零群、$H$ を $G$ の真部分群とすると、$H\subsetneq N_G(H)$ である。
:::

:::proof
$Z_0(G)=1\le H$ であり、$Z_c(G)=G\not\le H$ である（$c$ は冪零類）。$Z_i(G)\not\le H$ となる最小の $i$ をとると $i\ge1$ で $Z_{i-1}(G)\le H$ である。$z\in Z_i(G)\setminus H$ を一つとる。$h\in H$ について、:ref[lem-grp-central-basic] の 3 より $[h,z]\in Z_{i-1}(G)\le H$ なので $[z,h]=[h,z]^{-1}\in H$ であり、

$$
zhz^{-1}=[z,h]\,h\in H
$$

である。よって $zHz^{-1}\subseteq H$ である。$z^{-1}\in Z_i(G)\setminus H$ にも同じ議論を適用すれば $z^{-1}Hz\subseteq H$、すなわち $H\subseteq zHz^{-1}$ であり、$zHz^{-1}=H$ である。したがって $z\in N_G(H)\setminus H$ である。
:::

:::proposition[正規部分群と中心]{#prop-grp-normal-meets-center}
$G$ を冪零群、$N\ne1$ を $G$ の正規部分群とすると、$N\cap Z(G)\ne1$ である。特に非自明な冪零群の中心は非自明である。
:::

:::proof
$N\cap Z_c(G)=N\ne1$ なので、$N\cap Z_i(G)\ne1$ となる最小の $i$ が存在し、$i\ge1$ である。$i\ge2$ と仮定し、$x\in N\cap Z_i(G)$, $x\ne e$ をとる。任意の $g\in G$ について $[g,x]=(gxg^{-1})x^{-1}\in N$（$N$ は正規）であり、:ref[lem-grp-central-basic] の 3 より $[g,x]\in Z_{i-1}(G)$ である。$i$ の最小性から $N\cap Z_{i-1}(G)=1$ なので $[g,x]=e$ であり、$x\in Z(G)=Z_1(G)$ となる。これは $N\cap Z_1(G)=1$（$i\ge2$ の最小性）に反する。よって $i=1$ である。
:::

### 有限冪零群の構造

有限冪零群は、Sylow 部分群に分解して調べることができる。

:::theorem[有限冪零群の特徴付け]{#thm-grp-finite-nilpotent}
有限群 $G$ について次は同値である。

1. $G$ は冪零である。
2. $G$ のすべての真部分群 $H$ について $H\subsetneq N_G(H)$ である。
3. $G$ の極大部分群（真部分群で、それを真に含む真部分群がないもの）はすべて正規である。
4. $G$ のすべての Sylow 部分群は正規である。
5. $G$ は、$|G|$ を割る各素数についての Sylow 部分群の内部直積である。
6. 位数が互いに素な $G$ の二元は可換である。
:::

:::proof
**1 ⇒ 2。** :ref[prop-grp-normalizer-condition] である。

**2 ⇒ 3。** 極大部分群 $M$ について $M\subsetneq N_G(M)\le G$ であり、極大性から $N_G(M)=G$、すなわち $M\trianglelefteq G$ である。

**3 ⇒ 4。** Sylow $p$ 部分群 $P$ が正規でないと仮定すると $N_G(P)\ne G$ である。$G$ は有限なので、$N_G(P)$ を含む真部分群のうち位数最大のもの $M$ は極大部分群であり、3 より $M\trianglelefteq G$ である。$P\le M$ であり、$|M|$ の $p$ 部分は $|G|$ の $p$ 部分 $|P|$ 以下だから $P\in\operatorname{Syl}_p(M)$ である。:ref[prop-grp-frattini-argument] より $G=MN_G(P)=M$ となり、矛盾する。

**4 ⇒ 5。** $|G|=p_1^{a_1}\cdots p_k^{a_k}$（$p_i$ は相異なる素数）とし、$P_i$ を Sylow $p_i$ 部分群とする。正規な Sylow 部分群は（すべての Sylow $p_i$ 部分群が共役なので :ref[thm-grp-sylow-conjugacy]）ただ一つである。$Q_m=P_1\cdots P_m$ が位数 $p_1^{a_1}\cdots p_m^{a_m}$ の正規部分群であることを $m$ に関する帰納法で示す。$m=1$ は仮定である。$Q_m$ と $P_{m+1}$ はともに正規で、位数が互いに素なので :ref[thm-grp-lagrange] より $Q_m\cap P_{m+1}=1$ である。正規部分群の積 $Q_{m+1}=Q_mP_{m+1}$ は正規部分群で、:ref[prop-grp-product-size] より位数は $|Q_m||P_{m+1}|$ である。各段で $Q_m$ と $P_{m+1}$ は正規で $Q_m\cap P_{m+1}=1$ だから、:ref[prop-grp-internal-direct] により $Q_{m+1}\cong Q_m\times P_{m+1}$ であり、帰納的に $Q_m\cong P_1\times\cdots\times P_m$ である。$m=k$ のとき $|Q_k|=|G|$ なので $G=Q_k\cong P_1\times\cdots\times P_k$ である。

**5 ⇒ 1。** 各 $P_i$ は :ref[thm-grp-p-nilpotent] により冪零であり、:ref[prop-grp-nilpotent-closure] の 3 により直積も冪零である。

**5 ⇒ 6。** $G=P_1\times\cdots\times P_k$ の元を $x=(x_1,\ldots,x_k)$ と書くと、$x_i$ の位数は $p_i$ の冪であり、$x$ の位数はそれらの最小公倍数である。$x,y$ の位数が互いに素なら、各 $i$ について $x_i=e$ または $y_i=e$ である（両方が非単位元なら $p_i$ が両方の位数を割る）。よって各成分で $x_iy_i=y_ix_i$ であり、$xy=yx$ である。

**6 ⇒ 4。** $P\in\operatorname{Syl}_p(G)$ とする。$p$ と異なる素数 $q$ の Sylow $q$ 部分群 $Q$ の元の位数は $q$ の冪、$P$ の元の位数は $p$ の冪なので、6 より $Q\le C_G(P)\le N_G(P)$ である。また $P\le N_G(P)$ である。Sylow 部分群は各素数について存在する（:ref[thm-grp-sylow-existence]）から、$|N_G(P)|$ は $|G|$ を割る各素数の最高冪で割り切れ、$|N_G(P)|=|G|$、すなわち $P\trianglelefteq G$ である。
:::

:::corollary[p 群の極大部分群]{#cor-grp-p-maximal}
有限 $p$ 群 $P$ の極大部分群 $M$ は正規で、指数は $p$ である。
:::

:::proof
$P$ は冪零なので、定理の 3 により $M\trianglelefteq P$ である。:ref[thm-grp-correspondence] により $P/M$ の部分群は $M$ と $P$ の間の部分群に対応するので、$P/M$ は非自明で、自明群と自分自身以外の部分群をもたない。$P/M$ は $p$ 群なので :ref[thm-grp-cauchy] により位数 $p$ の部分群をもち、それが $P/M$ 全体である。よって $[P:M]=p$ である。
:::

:::remark[有限性の役割]
1 ⇒ 2 の証明は有限性を使っておらず、任意の冪零群で成り立つ。一方 3 ⇒ 4 では極大部分群の存在と Sylow の定理に有限性を用いた。また 5 から、有限冪零群 $G$ は $|G|$ の任意の約数 $d$ に対して位数 $d$ の部分群をもつことが分かる（演習問題）。$A_4$ は位数 $6$ の部分群をもたないので冪零でなく、実際 $A_4$ の Sylow $3$ 部分群は $4$ 個あって正規でない。
:::

## 演習問題

::::exercise{level=basic}
巡回群 $C_{18}$ の組成列をすべて求め、それぞれの組成因子を答えよ。
:::hint
:ref[thm-grp-cyclic-subgroups] により、$C_{18}$ の部分群は位数で決まる。
:::
:::solution
$C_{18}$ の部分群は各約数 $d\mid18$ に対して位数 $d$ のものがただ一つあり、位数 $d'\mid d$ の部分群は位数 $d$ の部分群に含まれる（:ref[exm-grp-cyclic-composition] の議論）。したがって組成列は、$18$ から $1$ まで各段で素数で割っていく位数の列と一対一に対応する。$18=2\cdot3\cdot3$ なので、割る素数の順序は $(2,3,3)$, $(3,2,3)$, $(3,3,2)$ の三通りであり、組成列は

$$
\begin{aligned}
&C_{18}\trianglerighteq\langle\bar2\rangle\trianglerighteq\langle\bar6\rangle\trianglerighteq0 &&(\text{位数 }18,9,3,1;\ \text{因子 }C_2,C_3,C_3),\\
&C_{18}\trianglerighteq\langle\bar3\rangle\trianglerighteq\langle\bar6\rangle\trianglerighteq0 &&(\text{位数 }18,6,3,1;\ \text{因子 }C_3,C_2,C_3),\\
&C_{18}\trianglerighteq\langle\bar3\rangle\trianglerighteq\langle\bar9\rangle\trianglerighteq0 &&(\text{位数 }18,6,2,1;\ \text{因子 }C_3,C_3,C_2)
\end{aligned}
$$

の三つである。可換なので各項は正規であり、因子は素数位数で単純である。Jordan–Hölder の定理のとおり、どの列でも組成因子は $C_2,C_3,C_3$ である。
:::
::::

::::exercise{level=basic}
$D_6$（位数 $12$）の導来列と、組成列を一つ求めよ。
:::solution
:ref[exm-grp-dihedral-commutator] より $D_6'=\langle r^2\rangle$ であり、これは位数 $3$ の巡回群なので $D_6''=1$ である。導来列は $D_6\trianglerighteq\langle r^2\rangle\trianglerighteq1$、因子は $(D_6)_{\mathrm{ab}}\cong C_2\times C_2$ と $C_3$ である。

組成列としては、例えば

$$
D_6\trianglerighteq\langle r\rangle\trianglerighteq\langle r^2\rangle\trianglerighteq1
$$

がとれる。$\langle r\rangle$ は指数 $2$ なので正規（:ref[prop-grp-index-two]）、$\langle r^2\rangle$ は巡回群 $\langle r\rangle$ の部分群なので正規である。因子の位数は $2,2,3$ で、組成因子は $C_2,C_2,C_3$ である。導来列の第一因子 $C_2\times C_2$ は単純でないので、導来列そのものは組成列ではない。
:::
::::

::::exercise{level=basic}
$Q_8$ の下中心列と上中心列を求めよ。
:::solution
:ref[exm-grp-quaternion-commutator] より $\gamma_2(Q_8)=Q_8'=\{\pm1\}$ であり、$-1$ は中心にあるので $\gamma_3(Q_8)=[Q_8,\{\pm1\}]=1$ である。

中心は $Z(Q_8)=\{\pm1\}$ である（$\pm i,\pm j,\pm k$ はそれぞれ別の単位と反可換で中心に入らない）。$Q_8/\{\pm1\}$ は位数 $4$ の Abel 群なので、その中心は全体であり $Z_2(Q_8)=Q_8$ である。よって下中心列は $Q_8\ge\{\pm1\}\ge1$、上中心列は $1\le\{\pm1\}\le Q_8$ であり、冪零類は $2$ である。
:::
::::

::::exercise{level=basic}
$C_6\times Q_8$ の冪零類を求めよ。
:::solution
:ref[prop-grp-nilpotent-closure] の 3 の証明より $\gamma_i(C_6\times Q_8)=\gamma_i(C_6)\times\gamma_i(Q_8)$ である。$C_6$ は可換なので $\gamma_2(C_6)=1$、前問より $\gamma_2(Q_8)=\{\pm1\}$, $\gamma_3(Q_8)=1$ である。よって $\gamma_2=1\times\{\pm1\}\ne1$, $\gamma_3=1$ であり、冪零類は $2$ である。
:::
::::

::::exercise{level=basic}
$S_3$ の下中心列と上中心列を求め、$S_3$ が可解だが冪零でないことを確かめよ。
:::solution
$S_3'=A_3$（:ref[exm-grp-symmetric-commutator]）なので $\gamma_2=A_3=\langle(123)\rangle$ である。$\gamma_3=[S_3,A_3]$ は $[(12),(123)]=(12)(123)(12)(132)=(132)(132)=(123)$ を含む（$(12)(123)(12)=(213)=(132)$）。よって $\gamma_3=A_3$ であり、以後 $\gamma_i=A_3$ のまま $1$ に達しない。

中心は $Z(S_3)=1$ である。実際、$(12)$ と可換な元は $\{e,(12)\}$ のみで（$(13)(12)=(123)\ne(132)=(12)(13)$ など）、同様に $(23)$ と可換な元は $\{e,(23)\}$ のみだから、両方と可換なのは $e$ だけである。よって $Z_1=1$ で、以後 $Z_i=1$ のまま $S_3$ に達しない。一方、導来列 $S_3\trianglerighteq A_3\trianglerighteq1$ は $1$ に達するので $S_3$ は可解である。
:::
::::

::::exercise{level=basic}
素数 $p$ について、位数 $p^2$ の群は Abel 群であることを本章の結果から示せ。
:::solution
$|G|=p^2$ とする。:ref[cor-grp-p-center] より $Z(G)\ne1$ なので、$|Z(G)|$ は $p$ または $p^2$ である。$|Z(G)|=p$ なら $G/Z(G)$ は位数 $p$ で巡回群なので、:ref[lem-grp-center-cyclic-quotient] より $G$ は Abel 群となり $Z(G)=G$ で $|Z(G)|=p^2$ となって矛盾する。よって $Z(G)=G$ であり、$G$ は Abel 群である。
:::
::::

::::exercise{level=advanced}
$G\ne1$ を可解群とする。可換化 $G_{\mathrm{ab}}$ は非自明であることを示せ。さらに $G$ が有限なら、$G$ は素数指数の正規部分群をもつことを示せ。
:::solution
$G_{\mathrm{ab}}=1$ なら $G'=G$ であり、帰納的に $G^{(i)}=G\ne1$ がすべての $i$ で成り立つので、可解性に反する。よって $G'\ne G$ である。

$G$ が有限とする。$A=G/G'$ は非自明な有限 Abel 群なので組成列をもち（:ref[prop-grp-composition-exists]）、その最初の段 $A\trianglerighteq A_1$ の因子 $A/A_1$ は単純 Abel 群、すなわち素数位数 $p$ の巡回群である（:ref[exm-grp-z-no-composition]）。:ref[thm-grp-correspondence] により $A_1=N/G'$ となる $N\trianglelefteq G$ があり、:ref[thm-grp-third-iso] より $G/N\cong A/A_1$ の位数は $p$ である。よって $N$ は指数 $p$ の正規部分群である。
:::
::::

::::exercise{level=advanced}
$D_8$（位数 $16$）の下中心列と上中心列を求め、両者が逆順に一致することを確かめよ。
:::solution
:ref[exm-grp-dihedral-nilpotent] より $\gamma_2=\langle r^2\rangle$, $\gamma_3=\langle r^4\rangle$, $\gamma_4=1$ である。

上中心列を求める。$r^k$ が中心にあるのは $sr^ks^{-1}=r^{-k}=r^k$、すなわち $8\mid2k$ のときであり、$r^ks$ は $r(r^ks)=r^{k+1}s\ne r^{k-1}s=(r^ks)r$ なので中心にない。よって $Z_1=\langle r^4\rangle=\{e,r^4\}$ である。:ref[lem-grp-central-basic] の 3 により $x\in Z_2$ と「すべての $g$ で $[g,x]\in\{e,r^4\}$」は同値である。$x=r^k$ なら $[r^j,r^k]=e$, $[r^js,r^k]=r^{-2k}$（:ref[exm-grp-dihedral-nilpotent] の計算）なので、条件は $r^{-2k}\in\{e,r^4\}$、すなわち $k$ が偶数であることである。$x=r^ks$ なら $[r,r^ks]=r\,r^ks\,r^{-1}\,s^{-1}r^{-k}=r^{k+1}\cdot r\cdot r^{-k}=r^2\notin\{e,r^4\}$ である（$sr^{-1}s^{-1}=r$ を用いた）。よって $Z_2=\langle r^2\rangle$ である。最後に $[g,x]\in D_8'=\langle r^2\rangle=Z_2$ はすべての $x$ で成り立つので $Z_3=D_8$ である。

まとめると、下中心列は $D_8\ge\langle r^2\rangle\ge\langle r^4\rangle\ge1$、上中心列は $1\le\langle r^4\rangle\le\langle r^2\rangle\le D_8$ であり、逆順に一致する。冪零類は $3$ である。一般には二つの列は一致するとは限らず、:eqref[eq-grp-central-sandwich] の不等式が成り立つだけである（例えば $C_2\times D_4$ では $\gamma_2=1\times\langle r^2\rangle$ だが $Z_1=C_2\times\langle r^2\rangle$ である）。
:::
::::

::::exercise{level=advanced}
$p$ を素数とし、$G$ を位数 $p^3$ の非可換群とする。$Z(G)=G'$ であり、これは位数 $p$ であること、および $G/Z(G)\cong C_p\times C_p$ であることを示せ。
:::solution
$Z(G)\ne1$（:ref[cor-grp-p-center]）で、$G$ は非可換なので $Z(G)\ne G$ である。$|Z(G)|=p^2$ なら $G/Z(G)$ は位数 $p$ の巡回群となり、:ref[lem-grp-center-cyclic-quotient] より $G$ が可換になって矛盾する。よって $|Z(G)|=p$, $|G/Z(G)|=p^2$ である。

位数 $p^2$ の群は Abel 群（前の演習問題）なので、:ref[prop-grp-abelianization] の 2 より $G'\le Z(G)$ である。$G$ は非可換なので $G'\ne1$ であり、$|Z(G)|=p$ だから $G'=Z(G)$ である。

$G/Z(G)$ が巡回群なら再び :ref[lem-grp-center-cyclic-quotient] に反するので、位数 $p^2$ の Abel 群 $G/Z(G)$ は位数 $p^2$ の元をもたず、非単位元の位数はすべて $p$ である。非単位元 $a$ と $b\notin\langle a\rangle$ をとれば $\langle a\rangle\cap\langle b\rangle=1$ で、:ref[prop-grp-internal-direct] より $G/Z(G)=\langle a\rangle\langle b\rangle\cong C_p\times C_p$ である。例えば $D_4$, $Q_8$ では $Z=G'=\{e,r^2\}$, $\{\pm1\}$ であり、Heisenberg 群 $\mathrm{UT}_3(\F_p)$ では $Z=U'=\{u(0,0,c)\}$ である。
:::
::::

::::exercise{level=advanced}
$A_5$ と $C_{60}$ はどちらも位数 $60$ である。両者の組成因子を比較し、組成因子が群の位数だけでは決まらないことを説明せよ。
:::solution
$A_5$ は単純なので（:ref[lem-grp-a-five]）、組成列は $A_5\trianglerighteq1$ だけであり、組成因子は $A_5$ 一つである。$C_{60}$ では位数 $60,30,15,5,1$ の部分群の列が組成列で、組成因子は $C_2,C_2,C_3,C_5$ である。組成長は $1$ と $4$ で異なる。前者は非可換単純群を因子にもつので可解でなく、後者は可換なので可解である。:ref[cor-grp-solvable-factors] のとおり、可解性は組成因子が素数位数かどうかで判定され、位数 $60=2^2\cdot3\cdot5$ だけでは判定できない。
:::
::::

::::exercise{level=advanced}
位数 $105=3\cdot5\cdot7$ の群は可解であることを示せ。
:::solution
$n_7\mid15$, $n_7\equiv1\pmod7$ より $n_7\in\{1,15\}$、$n_5\mid21$, $n_5\equiv1\pmod5$ より $n_5\in\{1,21\}$ である（:ref[thm-grp-sylow-count]）。$n_7=15$ かつ $n_5=21$ と仮定する。位数 $7$ の相異なる部分群は自明に交わるので位数 $7$ の元は $15\cdot6=90$ 個、同様に位数 $5$ の元は $21\cdot4=84$ 個あり、合計 $174>105$ となって矛盾する。よって $n_7=1$ または $n_5=1$ であり、正規な Sylow 部分群 $P$（位数 $7$ または $5$）がある。$G/P$ の位数は $15=3\cdot5$ または $21=3\cdot7$ なので :ref[prop-grp-pq-solvable] により可解であり、$P$ は巡回群で可解である。:ref[thm-grp-solvable-extension] の 4 より $G$ は可解である。
:::
::::

::::exercise{level=advanced}
群 $G$ で $G'\le Z(G)$（冪零類 $2$ 以下）とする。任意の $x,y,z\in G$ と整数 $n\ge1$ について次を示せ。

1. $[xy,z]=[x,z][y,z]$ および $[x^n,y]=[x,y]^n$。
2. $(xy)^n=x^ny^n[y,x]^{n(n-1)/2}$。
:::hint
交換子は中心にあるので、:ref[lem-grp-commutator-identities] の 3 で共役が消える。2 は $ab=[a,b]ba$ を使って帰納法で示す。
:::
:::solution
1. :ref[lem-grp-commutator-identities] の 3 より $[xy,z]=x[y,z]x^{-1}[x,z]$ であり、$[y,z]$ は中心にあるので $x[y,z]x^{-1}=[y,z]$ である。中心の元どうしは可換なので $[xy,z]=[x,z][y,z]$ を得る。$x^{n+1}=x\cdot x^n$ に適用すると $[x^{n+1},y]=[x,y][x^n,y]$ となり、帰納法により $[x^n,y]=[x,y]^n$ である。

2. $n=1$ では明らかである。$n$ で成り立つとすると

$$
(xy)^{n+1}=(xy)^n\,xy=x^ny^n[y,x]^{n(n-1)/2}\,xy .
$$

$[y,x]^{n(n-1)/2}$ は中心にあるので右端に移せる。また :ref[lem-grp-commutator-identities] の 1 と 1 の結果より $y^nx=[y^n,x]\,xy^n=[y,x]^n\,xy^n$ である。よって

$$
(xy)^{n+1}=x^n[y,x]^nxy^n\,y\,[y,x]^{n(n-1)/2}=x^{n+1}y^{n+1}[y,x]^{n+n(n-1)/2}
$$

であり、$n+\frac{n(n-1)}2=\frac{(n+1)n}2$ だから $n+1$ でも成り立つ。

例えば $p$ が奇素数のとき Heisenberg 群 $\mathrm{UT}_3(\F_p)$ では $x^p=y^p=e$ かつ $p\mid p(p-1)/2$ なので $(xy)^p=e$ となり、:ref[exm-grp-quaternion-heisenberg-nilpotent] の計算と一致する。
:::
::::

::::exercise{level=advanced}
有限冪零群 $G$ は、$|G|$ の任意の正の約数 $d$ に対して位数 $d$ の部分群をもつことを示せ。また $A_4$ が位数 $6$ の部分群をもたないことを示し、$A_4$ が冪零でないことを確かめよ。
:::solution
:ref[thm-grp-finite-nilpotent] より $G=P_1\times\cdots\times P_k$（$P_i$ は Sylow $p_i$ 部分群、$|P_i|=p_i^{a_i}$）である。$d=p_1^{b_1}\cdots p_k^{b_k}$（$0\le b_i\le a_i$）とする。:ref[thm-grp-sylow-existence] を $p$ 群 $P_i$ に適用すると位数 $p_i^{b_i}$ の部分群 $H_i\le P_i$ がある。$H_1\times\cdots\times H_k$ は $G$ の位数 $d$ の部分群である。

$A_4$ に位数 $6$ の部分群 $H$ があるとする。指数 $2$ なので $H\trianglelefteq A_4$（:ref[prop-grp-index-two]）であり、$A_4/H$ は位数 $2$ だから任意の $g$ について $(gH)^2=H$、すなわち $g^2\in H$ である。$3$ 巡回 $\sigma$ については $\sigma=\sigma^4=(\sigma^2)^2\in H$ なので、$A_4$ の $8$ 個の $3$ 巡回がすべて $H$ に入り、$|H|\ge9$ となって矛盾する。よって $A_4$ は位数 $6$ の部分群をもたず、前半より冪零でない。実際、$A_4$ の Sylow $3$ 部分群は $\langle(123)\rangle,\langle(124)\rangle,\langle(134)\rangle,\langle(234)\rangle$ の $4$ 個あり、正規でない。
:::
::::

::::exercise{level=advanced}
体 $K$ 上の可逆な上三角行列全体のなす群 $B_n(K)$ は可解であることを示せ。
:::solution
上三角行列の積の対角成分は、対角成分どうしの積である（$(AB)_{ii}=\sum_kA_{ik}B_{ki}$ で、$A_{ik}\ne0$ には $k\ge i$、$B_{ki}\ne0$ には $k\le i$ が必要）。したがって

$$
\delta:B_n(K)\to(K^\times)^n,\qquad A\mapsto(A_{11},\ldots,A_{nn})
$$

は準同型である（可逆な上三角行列の対角成分は行列式 $\prod A_{ii}\ne0$ より非零）。その核は対角成分がすべて $1$ の上三角行列全体 $\mathrm{UT}_n(K)$ であり、:ref[prop-grp-unitriangular] により冪零、したがって可解である（:ref[prop-grp-central-series]）。像は Abel 群 $(K^\times)^n$ の部分群で可解であり、:ref[thm-grp-first-iso] により $B_n(K)/\mathrm{UT}_n(K)$ はこれと同型である。:ref[thm-grp-solvable-extension] の 4 より $B_n(K)$ は可解である。$B_n(K)$ は一般に冪零ではない。例えば $K=\F_3$, $n=2$ のとき、対角行列 $\diag(-1,1)$ と $x_{12}(1)$ の交換子は $x_{12}(-2)=x_{12}(1)$ であり、$\gamma_i(B_2(\F_3))$ は $x_{12}(1)$ を含み続ける。
:::
::::

## 確認テスト

::::quiz
組成列をもつ群について、Jordan–Hölder の定理により組成列の選び方によらず定まるものはどれか。

- [ ] 組成列の各項の部分群
- [x] 組成因子の同型類の多重集合
- [ ] 組成因子が現れる順序
- [ ] 群の同型類

:::explanation
Jordan–Hölder の定理は、二つの組成列の長さが等しく、因子が順序と同型を除いて重複度込みで一致することを主張する。$C_{12}$ の二つの組成列（:ref[exm-grp-cyclic-composition]）では中間の部分群も因子の順序も異なる。また $C_6$ と $S_3$ は組成因子 $C_2,C_3$ を共有するが同型でないので、組成因子から群は定まらない。
:::
::::

::::quiz{answer="3" tol="0"}
$S_4$ の導来長はいくらか。

:::explanation
$S_4'=A_4$, $A_4'=V_4$, $V_4'=1$ なので導来列は $S_4\trianglerighteq A_4\trianglerighteq V_4\trianglerighteq1$ であり、初めて自明群になるのは $S_4^{(3)}$ である。よって導来長は $3$ である（:ref[exm-grp-derived-examples]）。
:::
::::

::::quiz{answer="4" tol="0"}
二面体群 $D_{16}$（位数 $32$）の冪零類はいくらか。

:::explanation
:ref[exm-grp-dihedral-nilpotent] より $\gamma_{i+1}(D_n)=\langle r^{2^i}\rangle$ である。$n=16=2^4$ では $\gamma_2,\gamma_3,\gamma_4,\gamma_5$ の位数が $8,4,2,1$ となり、$\gamma_5=1\ne\gamma_4$ なので冪零類は $4$ である。位数 $2^5$ の群の冪零類の上限 $5-1=4$（:ref[thm-grp-p-nilpotent]）を達成している。
:::
::::

::::quiz
有限群 $G$ について、「$G$ が冪零である」ことと同値な条件をすべて選べ。

- [ ] $G$ は可解である
- [x] $G$ のすべての Sylow 部分群は正規である
- [x] $G$ のすべての極大部分群は正規である
- [ ] $G$ の中心は非自明である

:::explanation
第二・第三は :ref[thm-grp-finite-nilpotent] の条件 4, 3 である。可解性は冪零性より弱く、$S_3$ は可解だが冪零でない。中心が非自明であることも十分でなく、$S_3\times C_2$ は中心 $1\times C_2$ をもつが、部分群 $S_3\times1$ が冪零でないので冪零でない（:ref[prop-grp-nilpotent-closure] の 1）。
:::
::::

::::quiz
正しいものをすべて選べ。

- [x] 有限 $p$ 群の上中心列は有限回で群全体に達する
- [x] 可解群の剰余群は可解である
- [ ] 正規部分群と剰余群がともに冪零な群は冪零である
- [x] 冪零類 $3$ の群の導来長は $2$ 以下である

:::explanation
第一は非自明な $p$ 群の中心が非自明であることを各剰余群に適用して得られる（:ref[thm-grp-p-nilpotent]）。第二は全射が導来列を導来列に写すことによる（:ref[thm-grp-solvable-extension]）。第三は誤りで、$A_3\trianglelefteq S_3$ では $A_3$ も $S_3/A_3$ も Abel 群だが $S_3$ は冪零でない。第四は :ref[prop-grp-gamma-commutator] より $G''\le\gamma_4(G)=1$ となることから正しい。
:::
::::

::::quiz{answer="6" tol="0"}
巡回群 $C_{360}$ の組成長（組成列の長さ）はいくらか。

:::explanation
$360=2^3\cdot3^2\cdot5$ であり、:ref[exm-grp-cyclic-composition] より $C_{360}$ の組成因子は素因数分解の各素数に対応する $C_2,C_2,C_2,C_3,C_3,C_5$ である。よって組成長は $3+2+1=6$ である。
:::
::::

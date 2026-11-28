---
title: テンソル積と核型 C*-環
description: 空間的テンソル積と普遍的テンソル積を構成し、核型性と単射性の関係を調べる。
notation:
  - symbol: 'A\odot B'
    meaning: 代数的テンソル積
  - symbol: 'A\otimes_{\min}B'
    meaning: 最小 C* テンソル積
  - symbol: 'A\otimes_{\max}B'
    meaning: 最大 C* テンソル積
  - symbol: 'M\overline\otimes N'
    meaning: 空間的 von Neumann テンソル積
---

## この章の目標

- Hilbert 空間のテンソル積とその内積を構成する。
- von Neumann 環のテンソル積と可換子の対応を説明する。
- 最小・最大 C* ノルムの定義と違いを理解する。
- 有限次元近似から核型性と単射性を調べる。

## Hilbert 空間のテンソル積

:::definition[Hilbert テンソル積]{#def-oa-hilbert-tensor}
代数的テンソル積 $H\odot K$ は、形式的な $\xi\otimes\eta$ を双線形関係で割った空間である。単純テンソルの間に $\inner{\xi\otimes\eta}{\xi'\otimes\eta'}=\inner{\xi}{\xi'}\inner{\eta}{\eta'}$ と定め、線形性と共役線形性で延長する。この内積で完備化した $H\otimes K$ を :term[Hilbert 空間のテンソル積]{en="Hilbert space tensor product" yomi="ひるべるとくうかんのてんそるせき"} という。
:::

内積の正定値性は、有限和の第一因子の張る空間に正規直交基底 $e_j$ を選び、$u=\sum_je_j\otimes\eta_j$ と直すと $\norm u^2=\sum_j\norm{\eta_j}^2$ となることから従う。$\eta_j$ は第一因子の座標汎関数で取り出せるため、表示の一意性も分かる。第二因子の基底で同じ分解を使うと $\norm{S\otimes1}\le\norm S$ である。従って $(S\otimes1)(1\otimes T)$ は有界で $\norm{S\otimes T}\le\norm S\norm T$ である。単位単純テンソルで両因子のノルム上限に近づけば等号である。随伴と積は単純テンソル上で計算して全空間に延長する。

:::example[積のノルム]
$S=\diag(2,-1)$、$T=\diag(3,4)$ とするとテンソルの対角値は $6,8,-3,-4$ で、ノルムは $8=2\cdot4$ である。一方 $\xi=e_1\otimes e_1+e_2\otimes e_2$ は階数二の係数行列を持つため、単純テンソル一つでは表せない。ノルムは $\sqrt2$ である。
:::

### Schmidt 分解と最良の階数一近似

有限表示 $u=\sum_{ij}a_{ij}e_i\otimes f_j$ に対して係数行列 $(a_{ij})$ を特異値分解すれば $u=\sum_{j=1}^rs_j\xi_j\otimes\eta_j$ と書ける。右特異ベクトルの成分を共役して第 2 因子に置く点に注意する。$\xi_j,\eta_j$ はそれぞれ正規直交系、$s_j>0$ は特異値である。任意の単位単純テンソル $\xi\otimes\eta$ に対して $|\inner{u}{\xi\otimes\eta}|\le s_1$ である。特異値の最大を $s_1$ とし、Cauchy–Schwarz を二つの座標列に適用すればよい。等号は $\xi_1,\eta_1$ で出る。そこで各単位方向に $u$ を直交射影してから誤差を最小化すると、最良の単純テンソルによる近似の誤差は $(\sum_{j>1}s_j^2)^{1/2}$ になる。

:::example[各列で長さを測る積ベクトル]
有限次元の単純テンソル $u=a\otimes b$ と、第二因子だけに作用する $1\otimes D$ を計算する。基底の順序は第一因子の添字を先に固定する。

座標 $u_{ij}=a_i b_j$ を添字の順に並べる。第 1 因子の一つの成分が、第 2 因子の全成分をスケールする。

$$
\begin{aligned}
a&=\begin{pmatrix}1\\2\end{pmatrix}\\
b&=\begin{pmatrix}3\\1\end{pmatrix}\\
u&=\begin{pmatrix}3\\1\\6\\2\end{pmatrix}
\end{aligned}
$$

ノルムの二乗を取ると $|a_i|^2|b_j|^2$ を全て足す形になる。有限和を二つに分ければノルムの積の二乗になる。

$$
\begin{aligned}
\norm a^2&=5\\
\norm b^2&=10\\
\norm u^2&=50
\end{aligned}
$$

$1\otimes D$ は第一因子の各成分を保ったまま、第二因子だけに対角行列 $D$ を掛ける。二つの方法で得るベクトルは同じである。

$$
\begin{aligned}
D&=\begin{pmatrix}1 & 0\\0 & 2\end{pmatrix}\\
(1\otimes D)u&=\begin{pmatrix}3\\2\\6\\4\end{pmatrix}\\
\norm{(1\otimes D)u}^2&=65
\end{aligned}
$$

この入力ではノルムが各因子のノルムから読める。一方、一般の和のノルムには別の単純テンソルとの内積も寄与する。
:::

:::example[第二因子の三座標]
有限次元の単純テンソル $u=a\otimes b$ と、第二因子だけに作用する $1\otimes D$ を計算する。基底の順序は第一因子の添字を先に固定する。

座標 $u_{ij}=a_i b_j$ を添字の順に並べる。第 1 因子の一つの成分が、第 2 因子の全成分をスケールする。

$$
\begin{aligned}
a&=\begin{pmatrix}1\\0\end{pmatrix}\\
b&=\begin{pmatrix}2\\1\\3\end{pmatrix}\\
u&=\begin{pmatrix}2\\1\\3\\0\\0\\0\end{pmatrix}
\end{aligned}
$$

ノルムの二乗を取ると $|a_i|^2|b_j|^2$ を全て足す形になる。有限和を二つに分ければノルムの積の二乗になる。

$$
\begin{aligned}
\norm a^2&=1\\
\norm b^2&=14\\
\norm u^2&=14
\end{aligned}
$$

$1\otimes D$ は第一因子の各成分を保ったまま、第二因子だけに対角行列 $D$ を掛ける。二つの方法で得るベクトルは同じである。

$$
\begin{aligned}
D&=\begin{pmatrix}1 & 0 & 0\\0 & -1 & 0\\0 & 0 & 2\end{pmatrix}\\
(1\otimes D)u&=\begin{pmatrix}2\\-1\\6\\0\\0\\0\end{pmatrix}\\
\norm{(1\otimes D)u}^2&=41
\end{aligned}
$$

この入力ではノルムが各因子のノルムから読める。一方、一般の和のノルムには別の単純テンソルとの内積も寄与する。
:::

:::example[符号を含むテンソル]
有限次元の単純テンソル $u=a\otimes b$ と、第二因子だけに作用する $1\otimes D$ を計算する。基底の順序は第一因子の添字を先に固定する。

座標 $u_{ij}=a_i b_j$ を添字の順に並べる。第 1 因子の一つの成分が、第 2 因子の全成分をスケールする。

$$
\begin{aligned}
a&=\begin{pmatrix}2\\1\end{pmatrix}\\
b&=\begin{pmatrix}-1\\1\end{pmatrix}\\
u&=\begin{pmatrix}-2\\2\\-1\\1\end{pmatrix}
\end{aligned}
$$

ノルムの二乗を取ると $|a_i|^2|b_j|^2$ を全て足す形になる。有限和を二つに分ければノルムの積の二乗になる。

$$
\begin{aligned}
\norm a^2&=5\\
\norm b^2&=2\\
\norm u^2&=10
\end{aligned}
$$

$1\otimes D$ は第一因子の各成分を保ったまま、第二因子だけに対角行列 $D$ を掛ける。二つの方法で得るベクトルは同じである。

$$
\begin{aligned}
D&=\begin{pmatrix}2 & 0\\0 & 3\end{pmatrix}\\
(1\otimes D)u&=\begin{pmatrix}-4\\6\\-2\\3\end{pmatrix}\\
\norm{(1\otimes D)u}^2&=65
\end{aligned}
$$

この入力ではノルムが各因子のノルムから読める。一方、一般の和のノルムには別の単純テンソルとの内積も寄与する。
:::

:::example[一つの座標に集中した入力]
有限次元の単純テンソル $u=a\otimes b$ と、第二因子だけに作用する $1\otimes D$ を計算する。基底の順序は第一因子の添字を先に固定する。

座標 $u_{ij}=a_i b_j$ を添字の順に並べる。第 1 因子の一つの成分が、第 2 因子の全成分をスケールする。

$$
\begin{aligned}
a&=\begin{pmatrix}0\\2\end{pmatrix}\\
b&=\begin{pmatrix}1\\2\end{pmatrix}\\
u&=\begin{pmatrix}0\\0\\2\\4\end{pmatrix}
\end{aligned}
$$

ノルムの二乗を取ると $|a_i|^2|b_j|^2$ を全て足す形になる。有限和を二つに分ければノルムの積の二乗になる。

$$
\begin{aligned}
\norm a^2&=4\\
\norm b^2&=5\\
\norm u^2&=20
\end{aligned}
$$

$1\otimes D$ は第一因子の各成分を保ったまま、第二因子だけに対角行列 $D$ を掛ける。二つの方法で得るベクトルは同じである。

$$
\begin{aligned}
D&=\begin{pmatrix}3 & 0\\0 & 2\end{pmatrix}\\
(1\otimes D)u&=\begin{pmatrix}0\\0\\6\\8\end{pmatrix}\\
\norm{(1\otimes D)u}^2&=100
\end{aligned}
$$

この入力ではノルムが各因子のノルムから読める。一方、一般の和のノルムには別の単純テンソルとの内積も寄与する。
:::

## von Neumann 環のテンソル積と可換子定理

:::definition[空間的テンソル積]{#def-oa-vn-tensor}
$M\subset B(H)$、$N\subset B(K)$ を単位的 von Neumann 環とする。作用素 $x\otimes y$（$x\in M,y\in N$）の線形包の $H\otimes K$ 上での二重可換子を :term[空間的テンソル積]{en="spatial tensor product" yomi="くうかんてきてんそるせき"} $M\overline\otimes N$ と書く。
:::

:::lemma[巡回ベクトル上で可換子を判定する]{#lem-oa-tensor-real-commutant}
$A,B\subset B(H)$ を互いに可換な単位的 von Neumann 環とし、$\xi$ を $A$ の巡回ベクトルとする。実内積を $\Re\inner{\cdot}{\cdot}$ とする。このとき $B=A'$ であることと $A_{\rm sa}\xi+iB_{\rm sa}\xi$ が実 Hilbert 空間として稠密であることは同値である。
:::

:::proof
$a\in A_{\rm sa}$、$b\in A'_{\rm sa}$ なら $ab$ が自己共役なので $\inner{a\xi}{b\xi}$ は実数である。従って $A_{\rm sa}\xi$ と $iA'_{\rm sa}\xi$ は実直交する。稠密性を仮定すると $b\in A'_{\rm sa}$ に対し $b_j\xi\to b\xi$ となる $b_j\in B_{\rm sa}$ を取れる。$c\in B'$ と $a,d\in A$ に対し
$$
\inner{cb a\xi}{d\xi}
=\lim_j\inner{cb_j a\xi}{d\xi}
=\lim_j\inner{c a\xi}{d b_j\xi}
=\inner{bc a\xi}{d\xi}.
$$
最初と最後の極限では $b_j$ と $A$ の可換性を使う。巡回性から $b\in B''=B$ となり、自己共役成分に分けて $A'=B$ を得る。

逆に $B=A'$ とする。$\zeta$ が $A_{\rm sa}\xi+iB_{\rm sa}\xi$ に実直交すると仮定する。$H\oplus H$ で $(\xi,\zeta)$ の対角な $A$ による巡回空間への直交射影 $P$ を取る。成分で書くと $P=\left(\begin{smallmatrix}p&r\\ r^*&q\end{smallmatrix}\right)$ であり、$p,q,r\in A'$ で $0\le p\le1$ である。まず $P(\xi,\zeta)=(\xi,\zeta)$ より $p\xi+r\zeta=\xi$ である。さらに $\zeta$ と $A_{\rm sa}\xi$ の実直交性から、$\inner{\zeta}{a\xi}+\inner{\xi}{a\zeta}=0$ が全 $a\in A$ で成り立つ。従って $P(\zeta,\xi)=0$ であり、$p\zeta+r\xi=0$ である。一方、$iB_{\rm sa}\xi$ との直交性から $\inner{b\xi}{\zeta}=\inner{b\zeta}{\xi}$ が自己共役 $b\in B$ で成立し、線形性により全 $b\in B$ に広がる。特に $r\in B$ を代入すると
$$
0\le\inner{p\zeta}{\zeta}
=-\inner{r\xi}{\zeta}
=-\inner{r\zeta}{\xi}
=-\inner{(1-p)\xi}{\xi}\le0 .
$$
従って $(1-p)\xi=0$ であり、$1-p\in A'$ と巡回性から $p=1$ である。すると $\inner{p\zeta}{\zeta}=0$ により $\zeta=0$ となる。
:::

:::lemma[実部分空間のテンソル積の稠密性]{#lem-oa-tensor-real-density}
$U\subset H,V\subset K$ を閉実部分空間とし、$U+iU$ と $V+iV$ がそれぞれ稠密とする。$U^\perp,V^\perp$ を実内積についての直交補とする。$u\otimes v$（$u\in U,v\in V$）と $i u'\otimes v'$（$u'\in U^\perp,v'\in V^\perp$）の実線形包は $H\otimes K$ に稠密である。
:::

:::proof
この線形包に実直交する $\zeta$ を取る。$F:H\to K$ を $\inner{Fh}{k}=\inner{\zeta}{h\otimes k}$ で定めると、有界な反線形作用素となる。二種類の直交条件は $F(U)\subset V^\perp$ と $F(U^\perp)\subset iV$ を与える。実 Hilbert 空間の随伴 $F^*$ も反線形であり、$F^*(V)\subset U^\perp$ が成り立つ。よって $D=F^*F$ は複素線形な正作用素であって、$D(U^\perp)\subset iU^\perp$ となる。従って $D^2(U^\perp)\subset U^\perp$ である。実係数多項式で平方根を近似すると $D=(D^2)^{1/2}$ 自体も $U^\perp$ を保つ。$U+iU$ の稠密性から $U^\perp\cap iU^\perp=\{0\}$ なので $D(U^\perp)=0$、従って $F(U^\perp)=0$ である。すると $F^*(V)\subset U$ でもあり、既に $F^*(V)\subset U^\perp$ なので $F^*(V)=0$ となる。反線形性と $V+iV$ の稠密性から $F=0$ であり、単純テンソルの全行列係数が消えるので $\zeta=0$ である。
:::

:::theorem[von Neumann 環のテンソル積の可換子定理]{#thm-oa-tensor-commutation}
$(M\overline\otimes N)'=M'\overline\otimes N'$ である。
:::
:::proof
右辺は各 $x\otimes y$ と可換であり、WOT 閉性から左辺に含まれる。$M=B(H)$ の場合、左辺の作用素は $B(H)\otimes1$ の行列単位と可換なので $1\otimes T$ と書ける。さらに $1\otimes N$ と可換となる条件は $T\in N'$ であり、逆包含が従う。

$M,N$ がそれぞれ巡回ベクトル $\xi,\eta$ を持つとする。実閉部分空間 $U=\overline{M_{\rm sa}\xi}$ と $V=\overline{N_{\rm sa}\eta}$ を取ると、:ref[lem-oa-tensor-real-commutant] と実直交性により $U^\perp=i\overline{M'_{\rm sa}\xi}$、$V^\perp=i\overline{N'_{\rm sa}\eta}$ である。$U,V$ の単純テンソルは $(M\overline\otimes N)_{\rm sa}(\xi\otimes\eta)$ の閉包に入り、$U^\perp,V^\perp$ の単純テンソルに $i$ を掛けたものは $i(M'\overline\otimes N')_{\rm sa}(\xi\otimes\eta)$ の閉包に入る。後者では二つの因子からそれぞれ $i$ が出て、積は実数 $-1$ になることを使った。また $U+iU$ と $V+iV$ はそれぞれ巡回性で稠密である。従って :ref[lem-oa-tensor-real-density] により二種類の実部分空間の和が稠密になる。$\xi\otimes\eta$ は積の環に巡回的なので、:ref[lem-oa-tensor-real-commutant] を積の環と $M'\overline\otimes N'$ に適用して等号を得る。

まず忠実正規状態を持つ環は、その GNS 表現では巡回分離表現になる。 その表現は正規である。実際、有界増大ネットの各行列係数を稠密な巡回ベクトルの軌道上で取ると、元の正規状態が上限を保存するため、表現も上限を保存する。これを別の任意の正規表現に移すため、次の増幅の事実を使う。$M\subset B(H)$ の忠実正規表現と正規正値汎関数 $\varphi$ に対し、前双対がトレース級の商であることから $\varphi=\Tr(K\,\cdot)|_M$ の持ち上げを取れる。実部を取れば $K=K^*$ とできる。$K$ の正部分が与える $\psi$ は $\varphi$ を支配し、$\psi$ は $H$ の可算直和上のベクトル汎関数である。実際 $K_+^{1/2}$ の列を並べたベクトルを用いればよい。稠密な $M\xi_\psi$ 上の形式 $\inner{T a\xi_\psi}{b\xi_\psi}=\varphi(b^*a)$ は $0\le T\le1$ の作用素を定め、左作用と可換である。従って $\varphi$ は同じ直和上の $T^{1/2}\xi_\psi$ によるベクトル汎関数になる。任意の正規表現を巡回部分空間の直交和に分けると、この結果から、もとの忠実正規表現の適当な増幅の還元部分空間に埋め込める。

二つの因子をそれぞれ増幅しても交換子定理の等式は保たれる。増幅の可換子は、増幅添字の各行列係数が元の可換子に属することから、元の可換子と増幅側の全作用素環の空間的テンソル積になる。有限個の添字への切断が有界強収束するので、この行列係数による判定は無限個の添字でも成立する。さらに $e$ が還元射影なら可換子は $eM'e$ に圧縮される。角の可換子の作用素を補空間で零に延長すれば元の可換子に入ることからである。この性質を両因子に適用すれば、忠実正規状態を持つ二つの環の任意の忠実正規表現で等式を得る。

一般の場合は $H$ を $M$ の巡回部分空間に直交分解する。各空間は自己共役な作用で還元されるため、極大な直交族を取ったときの補空間が非零なら、そこから別の巡回空間を選べて矛盾する。各空間への射影を $p_i\in M'$ とする。各角 $p_iM'p_i$ では、その巡回部分空間を生成した単位ベクトルが分離的なので忠実正規状態を与える。$N'$ でも同様に $q_j$ を取る。有限個の $p_i$ の和 $p_F$ による角も、各角の忠実正規状態を対角圧縮と合成した平均を取れば忠実正規状態を持つ。正元の対角圧縮が全て零なら、その平方根が全ての $p_iH$ を消すことから忠実性が従う。

そこで $p_F M'p_F$ と $q_G N'q_G$ について、既に示した交換子定理を適用して両辺の交換子を取る。圧縮の可換子の記述により、これは $p_F M$ と $q_G N$ の交換子定理である。元の $(M\overline\otimes N)'$ の任意の作用素を $p_F\otimes q_G$ で両側から圧縮すると、この有限部分の等式により $M'\overline\otimes N'$ に属する。これらの圧縮は元の作用素へ強収束するので、逆包含が一般に従う。
:::

$T\in B(H\otimes K)$、$\eta,\eta'\in K$ に対して $\inner{T_{\eta,\eta'}\xi}{\xi'}=\inner{T(\xi\otimes\eta)}{\xi'\otimes\eta'}$ と定めると $T_{\eta,\eta'}\in B(H)$ である。この作用素を $K$ の行列係数によるスライスと呼ぶ。$T$ が $1\otimes N$ と可換なら、$H$ の行列係数による $K$ 側のスライスは $N'$ に属する。同様に $M\otimes1$ との可換性から $H$ 側のスライスは $M'$ に属する。片側が全作用素環であれば上の行列単位による証明のようにスライスを組み立て直せる。

### スライスによる係数回収

:::example[行列側の各成分を抜き出す]
$M\overline\otimes M_n$ の元 $T$ に $K=\C^n$ の標準基底を使う。$T_{ij}$ を $T(\xi\otimes e_j)$ の第 $i$ 成分とすると、$T=\sum_{i,j}T_{ij}\otimes E_{ij}$ である。各スライスは WOT 連続なので、$T$ をテンソルの有限和の WOT 極限として取れば $T_{ij}\in M$ である。逆に $T_{ij}\in M$ ならこの和は代数的テンソル積に属する。従って片側が有限行列の場合、積のすべての元を有限個の係数で表せる。例えば $M=L^\infty(X)$ では $T_{ij}$ は可測関数であり、$T$ は $X$ 上の本質的有界な行列関数になる。
:::

:::example[二つの行列単位に分解する作用素]
$M_2\overline\otimes M_2$ を $\C^2\otimes\C^2$ 上に置く。第二因子の行列係数によるスライスをとって、作用素 $x$ から第一因子の行列を回収する。

$x=A\otimes E_{12}+B\otimes E_{21}$ とする。テンソルの順序は $(1,1),(1,2),(2,1),(2,2)$ である。従って片方の添字を固定すると対応する行と列は飛び飛びになる。

$$
\begin{aligned}
A&=\begin{pmatrix}2 & 0\\0 & 1\end{pmatrix}\\
B&=\begin{pmatrix}0 & 1\\1 & 0\end{pmatrix}\\
x&=\begin{pmatrix}0 & 2 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\\1 & 0 & 0 & 0\end{pmatrix}
\end{aligned}
$$

第二因子の第 $j$ ベクトルを入力、第 $i$ を出力に取ると $x_{ij}$ を得る。例えば $x_{12}$ は全行列行の第 $1,3$ 行と第 $2,4$ 列から読む。

$$
\begin{aligned}
x_{12}&=\begin{pmatrix}2 & 0\\0 & 1\end{pmatrix}\\
x_{21}&=\begin{pmatrix}0 & 1\\1 & 0\end{pmatrix}\\
x_{11}&=\begin{pmatrix}0 & 0\\0 & 0\end{pmatrix}\\
x_{22}&=\begin{pmatrix}0 & 0\\0 & 0\end{pmatrix}
\end{aligned}
$$

この二つの係数を使って作用を再構成できる。さらに第一因子の行列 $C$ と $x$ が可換である条件は、$[C,A]=[C,B]=0$ である。積の交換子のそれぞれのスライスを取れば必要性、係数が全て零なら十分性が分かる。

片方が有限行列環なら、スライスは有限表で作用を完全に決める。可換子定理の一般の場合には、二つの因子の強閉包について別の議論も必要になる。
:::

:::example[二つの向きの異なるスライス]
$M_2\overline\otimes M_2$ を $\C^2\otimes\C^2$ 上に置く。第二因子の行列係数によるスライスをとって、作用素 $x$ から第一因子の行列を回収する。

$x=A\otimes E_{12}+B\otimes E_{21}$ とする。テンソルの順序は $(1,1),(1,2),(2,1),(2,2)$ である。従って片方の添字を固定すると対応する行と列は飛び飛びになる。

$$
\begin{aligned}
A&=\begin{pmatrix}1 & 0\\0 & 2\end{pmatrix}\\
B&=\begin{pmatrix}0 & 1\\1 & 0\end{pmatrix}\\
x&=\begin{pmatrix}0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 2\\1 & 0 & 0 & 0\end{pmatrix}
\end{aligned}
$$

第二因子の第 $j$ ベクトルを入力、第 $i$ を出力に取ると $x_{ij}$ を得る。例えば $x_{12}$ は全行列行の第 $1,3$ 行と第 $2,4$ 列から読む。

$$
\begin{aligned}
x_{12}&=\begin{pmatrix}1 & 0\\0 & 2\end{pmatrix}\\
x_{21}&=\begin{pmatrix}0 & 1\\1 & 0\end{pmatrix}\\
x_{11}&=\begin{pmatrix}0 & 0\\0 & 0\end{pmatrix}\\
x_{22}&=\begin{pmatrix}0 & 0\\0 & 0\end{pmatrix}
\end{aligned}
$$

この二つの係数を使って作用を再構成できる。さらに第一因子の行列 $C$ と $x$ が可換である条件は、$[C,A]=[C,B]=0$ である。積の交換子のそれぞれのスライスを取れば必要性、係数が全て零なら十分性が分かる。

片方が有限行列環なら、スライスは有限表で作用を完全に決める。可換子定理の一般の場合には、二つの因子の強閉包について別の議論も必要になる。
:::

:::example[負の係数を持つスライス]
$M_2\overline\otimes M_2$ を $\C^2\otimes\C^2$ 上に置く。第二因子の行列係数によるスライスをとって、作用素 $x$ から第一因子の行列を回収する。

$x=A\otimes E_{12}+B\otimes E_{21}$ とする。テンソルの順序は $(1,1),(1,2),(2,1),(2,2)$ である。従って片方の添字を固定すると対応する行と列は飛び飛びになる。

$$
\begin{aligned}
A&=\begin{pmatrix}3 & 0\\0 & -1\end{pmatrix}\\
B&=\begin{pmatrix}0 & 1\\1 & 0\end{pmatrix}\\
x&=\begin{pmatrix}0 & 3 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & -1\\1 & 0 & 0 & 0\end{pmatrix}
\end{aligned}
$$

第二因子の第 $j$ ベクトルを入力、第 $i$ を出力に取ると $x_{ij}$ を得る。例えば $x_{12}$ は全行列行の第 $1,3$ 行と第 $2,4$ 列から読む。

$$
\begin{aligned}
x_{12}&=\begin{pmatrix}3 & 0\\0 & -1\end{pmatrix}\\
x_{21}&=\begin{pmatrix}0 & 1\\1 & 0\end{pmatrix}\\
x_{11}&=\begin{pmatrix}0 & 0\\0 & 0\end{pmatrix}\\
x_{22}&=\begin{pmatrix}0 & 0\\0 & 0\end{pmatrix}
\end{aligned}
$$

この二つの係数を使って作用を再構成できる。さらに第一因子の行列 $C$ と $x$ が可換である条件は、$[C,A]=[C,B]=0$ である。積の交換子のそれぞれのスライスを取れば必要性、係数が全て零なら十分性が分かる。

片方が有限行列環なら、スライスは有限表で作用を完全に決める。可換子定理の一般の場合には、二つの因子の強閉包について別の議論も必要になる。
:::

:::example[一つのスライスに核がある]
$M_2\overline\otimes M_2$ を $\C^2\otimes\C^2$ 上に置く。第二因子の行列係数によるスライスをとって、作用素 $x$ から第一因子の行列を回収する。

$x=A\otimes E_{12}+B\otimes E_{21}$ とする。テンソルの順序は $(1,1),(1,2),(2,1),(2,2)$ である。従って片方の添字を固定すると対応する行と列は飛び飛びになる。

$$
\begin{aligned}
A&=\begin{pmatrix}0 & 0\\0 & 2\end{pmatrix}\\
B&=\begin{pmatrix}0 & 1\\1 & 0\end{pmatrix}\\
x&=\begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 2\\1 & 0 & 0 & 0\end{pmatrix}
\end{aligned}
$$

第二因子の第 $j$ ベクトルを入力、第 $i$ を出力に取ると $x_{ij}$ を得る。例えば $x_{12}$ は全行列行の第 $1,3$ 行と第 $2,4$ 列から読む。

$$
\begin{aligned}
x_{12}&=\begin{pmatrix}0 & 0\\0 & 2\end{pmatrix}\\
x_{21}&=\begin{pmatrix}0 & 1\\1 & 0\end{pmatrix}\\
x_{11}&=\begin{pmatrix}0 & 0\\0 & 0\end{pmatrix}\\
x_{22}&=\begin{pmatrix}0 & 0\\0 & 0\end{pmatrix}
\end{aligned}
$$

この二つの係数を使って作用を再構成できる。さらに第一因子の行列 $C$ と $x$ が可換である条件は、$[C,A]=[C,B]=0$ である。積の交換子のそれぞれのスライスを取れば必要性、係数が全て零なら十分性が分かる。

片方が有限行列環なら、スライスは有限表で作用を完全に決める。可換子定理の一般の場合には、二つの因子の強閉包について別の議論も必要になる。
:::

## C*-環の代数的テンソル積と C*-ノルム

$A\odot B$ には $(a\otimes b)(c\otimes d)=ac\otimes bd$、$(a\otimes b)^*=a^*\otimes b^*$ を入れる。以下ではまず両環が単位的な場合を扱う。忠実表現 $\pi,\rho$ により $\sum_i a_i\otimes b_i$ を $\sum_i\pi(a_i)\otimes\rho(b_i)$ で表すと単射である。実際、$a_i$ を一次独立に取り直し、第一因子の行列係数を使う。それらの制限の線形包は有限次元空間 $\spn\{\pi(a_i)\}$ の双対全体になる。そうでなければ非零な作用素が全行列係数で消える。各 $a_i$ を取り出す係数で $b_i=0$ が従う。この作用素ノルムは従って C* ノルムである。

:::example[行列環による増幅]
$M_n\odot B$ の元は $\sum_{i,j}E_{ij}\otimes b_{ij}$ と一意に書け、その積と随伴は $B$ の元を成分に持つ行列の積と随伴である。空間的ノルムは $K^n$ 上のブロック作用素のノルムになる。$M_2\odot M_3$ では行列単位のテンソル $36$ 個が $M_6$ の行列単位となる。
:::

### 代数的テンソル積の表示の一意性

有限和 $x=\sum_{j=1}^m a_j\otimes b_j$ で第一因子を一次独立に取り直す。各 $j$ に対する座標汎関数 $\ell_j$ を $\spn\{a_i\}$ 上で $\ell_j(a_i)=\delta_{ij}$ とすると、$(\ell_j\otimes\id)x=b_j$ である。従って $x=0$ なら全係数 $b_j$ が零であり、逆も成り立つ。第一因子を忠実に表現したとき、その有限次元部分空間をベクトル行列係数は分離する。実際、全行列係数が零なら作用素も零である。有限次元なので、制限した行列係数の線形包が双対全体になる。従って空間的表示で $x$ が零なら、行列係数を使って各 $b_j$ の忠実像が零と分かり、$x=0$ となる。

:::example[対角に二成分が残る係数行列]
$\C^2\odot\C^2$ を座標ごとの積を持つ可換代数の代数的テンソル積とする。係数行列 $X$ で $x=\sum_{i,j}X_{ij}e_i\otimes e_j$ を表し、最小の単純テンソルの個数を求める。

$e_1,e_2$ が基底なので、第一因子の二つの座標汎関数を作用させると係数の行を取り出せる。そのうち非零の行を使った表示が直ちに得られる。

$$
\begin{aligned}
X&=\begin{pmatrix}1 & 0\\0 & 2\end{pmatrix}\\
(e_1^*\otimes\id)x&=\begin{pmatrix}1\\0\end{pmatrix}\\
(e_2^*\otimes\id)x&=\begin{pmatrix}0\\2\end{pmatrix}
\end{aligned}
$$

単純テンソル $u\otimes v$ の係数は $uv^{\mathsf T}$ で階数が高々一である。従って $r$ 個の単純テンソルで表示できるなら係数行列の階数は高々 $r$ である。逆に独立な列から列空間の基底を取り、それぞれの座標を第二因子に置けば階数個の表示が得られる。

$$
\begin{aligned}
\det X&=2\\
\rank X&=2
\end{aligned}
$$

積と対合は係数ごとの積と共役である。従って $x^*x$ の係数は各成分の絶対値の二乗になり、テンソルの個数と正性は別の計算になる。

$$
\begin{aligned}
\operatorname{coeff}(x^*x)&=\begin{pmatrix}1 & 0\\0 & 4\end{pmatrix}
\end{aligned}
$$

ここでいう階数は係数行列の行列階数である。可換な表現で作用素としての階数を調べるなら、各点での非零値を数えることになる。
:::

:::example[一次独立性が消える係数行列]
$\C^2\odot\C^2$ を座標ごとの積を持つ可換代数の代数的テンソル積とする。係数行列 $X$ で $x=\sum_{i,j}X_{ij}e_i\otimes e_j$ を表し、最小の単純テンソルの個数を求める。

$e_1,e_2$ が基底なので、第一因子の二つの座標汎関数を作用させると係数の行を取り出せる。そのうち非零の行を使った表示が直ちに得られる。

$$
\begin{aligned}
X&=\begin{pmatrix}1 & 2\\2 & 4\end{pmatrix}\\
(e_1^*\otimes\id)x&=\begin{pmatrix}1\\2\end{pmatrix}\\
(e_2^*\otimes\id)x&=\begin{pmatrix}2\\4\end{pmatrix}
\end{aligned}
$$

単純テンソル $u\otimes v$ の係数は $uv^{\mathsf T}$ で階数が高々一である。従って $r$ 個の単純テンソルで表示できるなら係数行列の階数は高々 $r$ である。逆に独立な列から列空間の基底を取り、それぞれの座標を第二因子に置けば階数個の表示が得られる。

$$
\begin{aligned}
\det X&=0\\
\rank X&=1
\end{aligned}
$$

積と対合は係数ごとの積と共役である。従って $x^*x$ の係数は各成分の絶対値の二乗になり、テンソルの個数と正性は別の計算になる。

$$
\begin{aligned}
\operatorname{coeff}(x^*x)&=\begin{pmatrix}1 & 4\\4 & 16\end{pmatrix}
\end{aligned}
$$

ここでいう階数は係数行列の行列階数である。可換な表現で作用素としての階数を調べるなら、各点での非零値を数えることになる。
:::

:::example[同じ第二因子をまとめる]
$\C^2\odot\C^2$ を座標ごとの積を持つ可換代数の代数的テンソル積とする。係数行列 $X$ で $x=\sum_{i,j}X_{ij}e_i\otimes e_j$ を表し、最小の単純テンソルの個数を求める。

$e_1,e_2$ が基底なので、第一因子の二つの座標汎関数を作用させると係数の行を取り出せる。そのうち非零の行を使った表示が直ちに得られる。

$$
\begin{aligned}
X&=\begin{pmatrix}1 & 3\\1 & 3\end{pmatrix}\\
(e_1^*\otimes\id)x&=\begin{pmatrix}1\\3\end{pmatrix}\\
(e_2^*\otimes\id)x&=\begin{pmatrix}1\\3\end{pmatrix}
\end{aligned}
$$

単純テンソル $u\otimes v$ の係数は $uv^{\mathsf T}$ で階数が高々一である。従って $r$ 個の単純テンソルで表示できるなら係数行列の階数は高々 $r$ である。逆に独立な列から列空間の基底を取り、それぞれの座標を第二因子に置けば階数個の表示が得られる。

$$
\begin{aligned}
\det X&=0\\
\rank X&=1
\end{aligned}
$$

積と対合は係数ごとの積と共役である。従って $x^*x$ の係数は各成分の絶対値の二乗になり、テンソルの個数と正性は別の計算になる。

$$
\begin{aligned}
\operatorname{coeff}(x^*x)&=\begin{pmatrix}1 & 9\\1 & 9\end{pmatrix}
\end{aligned}
$$

ここでいう階数は係数行列の行列階数である。可換な表現で作用素としての階数を調べるなら、各点での非零値を数えることになる。
:::

:::example[二つの独立な係数を持つ場合]
$\C^2\odot\C^2$ を座標ごとの積を持つ可換代数の代数的テンソル積とする。係数行列 $X$ で $x=\sum_{i,j}X_{ij}e_i\otimes e_j$ を表し、最小の単純テンソルの個数を求める。

$e_1,e_2$ が基底なので、第一因子の二つの座標汎関数を作用させると係数の行を取り出せる。そのうち非零の行を使った表示が直ちに得られる。

$$
\begin{aligned}
X&=\begin{pmatrix}1 & 2\\3 & 4\end{pmatrix}\\
(e_1^*\otimes\id)x&=\begin{pmatrix}1\\2\end{pmatrix}\\
(e_2^*\otimes\id)x&=\begin{pmatrix}3\\4\end{pmatrix}
\end{aligned}
$$

単純テンソル $u\otimes v$ の係数は $uv^{\mathsf T}$ で階数が高々一である。従って $r$ 個の単純テンソルで表示できるなら係数行列の階数は高々 $r$ である。逆に独立な列から列空間の基底を取り、それぞれの座標を第二因子に置けば階数個の表示が得られる。

$$
\begin{aligned}
\det X&=-2\\
\rank X&=2
\end{aligned}
$$

積と対合は係数ごとの積と共役である。従って $x^*x$ の係数は各成分の絶対値の二乗になり、テンソルの個数と正性は別の計算になる。

$$
\begin{aligned}
\operatorname{coeff}(x^*x)&=\begin{pmatrix}1 & 4\\9 & 16\end{pmatrix}
\end{aligned}
$$

ここでいう階数は係数行列の行列階数である。可換な表現で作用素としての階数を調べるなら、各点での非零値を数えることになる。
:::

:::example[複素成分で階数一になる場合]
$\C^2\odot\C^2$ を座標ごとの積を持つ可換代数の代数的テンソル積とする。係数行列 $X$ で $x=\sum_{i,j}X_{ij}e_i\otimes e_j$ を表し、最小の単純テンソルの個数を求める。

$e_1,e_2$ が基底なので、第一因子の二つの座標汎関数を作用させると係数の行を取り出せる。そのうち非零の行を使った表示が直ちに得られる。

$$
\begin{aligned}
X&=\begin{pmatrix}1 & i\\i & -1\end{pmatrix}\\
(e_1^*\otimes\id)x&=\begin{pmatrix}1\\i\end{pmatrix}\\
(e_2^*\otimes\id)x&=\begin{pmatrix}i\\-1\end{pmatrix}
\end{aligned}
$$

単純テンソル $u\otimes v$ の係数は $uv^{\mathsf T}$ で階数が高々一である。従って $r$ 個の単純テンソルで表示できるなら係数行列の階数は高々 $r$ である。逆に独立な列から列空間の基底を取り、それぞれの座標を第二因子に置けば階数個の表示が得られる。

$$
\begin{aligned}
\det X&=0\\
\rank X&=1
\end{aligned}
$$

積と対合は係数ごとの積と共役である。従って $x^*x$ の係数は各成分の絶対値の二乗になり、テンソルの個数と正性は別の計算になる。

$$
\begin{aligned}
\operatorname{coeff}(x^*x)&=\begin{pmatrix}1 & 1\\1 & 1\end{pmatrix}
\end{aligned}
$$

ここでいう階数は係数行列の行列階数である。可換な表現で作用素としての階数を調べるなら、各点での非零値を数えることになる。
:::

## 最小テンソル積と最大テンソル積

:::definition[空間的ノルムと普遍ノルム]{#def-oa-min-max-tensor}
忠実表現のテンソル作用素のノルムを $\norm{\cdot}_{\min}$ とする。共通の Hilbert 空間上で、像が互いに可換な全ての表現 $\alpha,\beta$ に対する $\norm{\sum_i\alpha(a_i)\beta(b_i)}$ の上限を $\norm{\cdot}_{\max}$ とする。対応する完備化を :term[最小テンソル積]{en="minimal tensor product" yomi="さいしょうてんそるせき"}、:term[最大テンソル積]{en="maximal tensor product" yomi="さいだいてんそるせき"} という。
:::

:::proposition[表現の選択と比較]{#prop-oa-tensor-norms}
空間的ノルムは忠実表現の選択に依存せず、任意のテンソル表現のノルムを抑える。また $\norm x_{\min}\le\norm x_{\max}\le\sum_i\norm{a_i}\norm{b_i}$ である。両者は C* ノルムで、恒等写像は $A\otimes_{\max}B$ から $A\otimes_{\min}B$ への全射になる。
:::
:::proof
忠実な単位的表現 $\pi$ について、有限個のベクトル汎関数の凸結合は状態空間で弱 * 稠密である。もし状態をその弱 * 閉凸包から分離できれば、ある自己共役 $a$ で状態の値が $\sup_{\norm\xi=1}\inner{\pi(a)\xi}{\xi}$ より大きくなる。右辺はスペクトルの上端であり、全状態がその以下の値を取るため矛盾である。同じ結論は行列増幅にも成り立つ。

別の表現 $\alpha$ と有限ベクトル族 $(\xi_j)$ に対し、$[a_{ij}]\mapsto\sum_{i,j}\inner{\alpha(a_{ij})\xi_j}{\xi_i}$ は $M_n(A)$ の正値汎関数である。それを上の凸結合で近似すると、$\alpha$ の有限個の行列係数を $\pi$ の直和のベクトル族で任意に近似できる。$x^*x$ と $\sum_j\xi_j\otimes\eta_j$ に適用すると、$x$ のテンソル表現のノルムが忠実表現を第一因子に使ったノルム以下と分かる。右因子も同様である。非退化でない表現の零作用部分は無視できる。非単位的な環には単位化の中で同じ証明を適用する。

最大ノルムの上界は三角不等式であり、空間表現は $\pi(a)\otimes1$ と $1\otimes\rho(b)$ を用いてその上限に含まれる。各表現で C* 条件が成り立つため二乗と上限を交換できる。よって完備化は C*-環である。恒等写像の延長は *-準同型で、像は代数的テンソル積を含み、閉性から全射となる。
:::

:::example[行列因子のノルムは一つ]
$M_n$ の非零単位的表現は $a\mapsto a\otimes1_L$ の形である。その可換子は $1\otimes B(L)$ なので、可換なもう一方の表現は $1\otimes\rho(b)$ の形である。従って最大ノルムは最小ノルムに等しい。有限直和についても中心射影で各成分に分け、ノルムの最大を取れば同じことが分かる。
:::

:::proposition[任意の C* ノルムとの比較]{#prop-oa-tensor-extremal}
$A,B$ が単位的とする。$A\odot B$ 上の任意の C* ノルム $\norm{\cdot}_\gamma$ に対し $\norm x_{\min}\le\norm x_\gamma\le\norm x_{\max}$ である。
:::

:::proof
完備化 $C_\gamma$ の忠実表現を取り、二つの因子 $a\otimes1$ と $1\otimes b$ に制限した表現を使えば上界が出る。下界には、二つの状態 $\varphi,\psi$ が与えられると $\varphi(a)\psi(b)$ が $C_\gamma$ 上の状態に延びることを示せばよい。

自己共役 $h\in A,k\in B$ とし、$h\otimes1$ と $1\otimes k$ が生成する可換部分環を取る。その指標空間は各指標を二つの座標で識別するため $\sigma(h)\times\sigma(k)$ の閉部分集合と見なせる。この集合は全ての点を含む。実際、外の点があればその点の積近傍も集合から離れて選べる。各座標のスペクトル上でその近傍に台を持ちその点で非零となる関数 $f,g$ を取ると、$f(h)\otimes g(k)$ は全指標で消えるため零となる。しかし両因子は非零で、単純テンソルは代数的テンソル積で非零である。$\norm{\cdot}_\gamma$ はノルムなので矛盾する。従って $h\otimes1+1\otimes k$ のスペクトルは $\sigma(h)+\sigma(k)$ である。

次に自己共役な和の部分空間に $\omega(h\otimes1+1\otimes k)=\varphi(h)+\psi(k)$ と置く。和の表示のずれは、二因子の交わりがスカラーであるため一方から他方へ移した定数だけである。従って値は一意に定まる。スペクトルの端点の足し算から $\omega$ の値の絶対値は和のノルム以下である。実 Hahn–Banach でノルムを保って延長し、複素化すれば単位元の値 $1$ を持つ正値汎関数を得る。正値性は正縮小元 $x$ に $\norm{1-x}\le1$ を使って分かる。

もし $\varphi$ が純粋なら、$b\ge0$ に対し $a\mapsto\omega(a\otimes b)$ は $0$ と $\norm b\,\varphi$ の間の正値汎関数である。これは積の可換性から $0\le a\otimes b\le\norm b(a\otimes1)$ が $a\ge0$ で成り立つためである。:ref[thm-oa-pure-gns] の証明で示した純粋性の判定により、この汎関数は $\varphi$ の倍数である。$a=1$ の値を取ると係数は $\psi(b)$ となる。正元で線形展開すれば全 $b$ で積状態を得る。一般の $\varphi$ は純粋状態の凸結合の弱 * 極限として近似できる。対応する積状態の弱 * クラスター点は積状態を延長する。

積状態の GNS 内積は、両因子の巡回ベクトルのテンソルから得る内積と等しい。従ってその表現は $\pi_\varphi\otimes\pi_\psi$ であり、どの組についても $\gamma$ ノルムで縮小的である。全状態の GNS 直和を両因子の忠実表現として取った空間的ノルムは、これらのテンソル表現のノルムの上限である。従って下界も得られる。
:::

### 可換因子をもつテンソルノルム

コンパクト Hausdorff 空間 $X$ と C*-環 $B$ に対し、$x=\sum_i f_i\otimes b_i$ は $x(t)=\sum_i f_i(t)b_i$ という連続 $B$ 値関数を与える。最小ノルムは全点評価の忠実直和と $B$ の忠実表現を取り、$\sup_t\norm{x(t)}$ になる。最大ノルムにも同じ上界がある。各 $f_i$ の変動を $\eps$ 以下にする開被覆を取り、その従属する分割 $h_j$ と各開集合の点 $t_j$ で $x$ を $\sum_jh_j\otimes x(t_j)$ に近似する。誤差の最大ノルムは $\eps\sum_i\norm{b_i}$ 以下である。任意の可換表現の組 $(\pi,\rho)$ について、後者の像は
$$
V^*\diag(\rho(x(t_j)))V,\qquad V\xi=(\pi(h_j)^{1/2}\xi)_j
$$
となる。$V$ は縮小的で、単位的表現なら等長である。従って最大ノルムも $\sup_t\norm{x(t)}$ 以下となる。ここで使った分割の存在は :ref[thm-top-urysohn]{text="Urysohn の補題"}から、有限被覆を縮めて得た関数の和で正規化すればよい。

:::example[最大値が一つの点で出る]
有限可換代数 $\C^2$ の二つの因子を、二点上の乗算として表現する。$x=f\otimes g+h\otimes k$ の値を並べ、空間的ノルムを求める。

点の対 $(i,j)$ での値は $f_i g_j+h_i k_j$ である。積の二項を外積の行列で表示すると、座標ごとの寄与を分けて確認できる。

$$
\begin{aligned}
(f_i g_j)_{ij}&=\begin{pmatrix}2 & 1\\-4 & -2\end{pmatrix}\\
(h_i k_j)_{ij}&=\begin{pmatrix}-1 & 2\\-3 & 6\end{pmatrix}\\
(x_{ij})&=\begin{pmatrix}1 & 3\\-7 & 4\end{pmatrix}
\end{aligned}
$$

表現先では $(i,j)$ に対応する基底ベクトルは固有ベクトルである。従って全ての値の絶対値の最大が作用素ノルムになる。

$$
\begin{aligned}
\sigma(x)&=\left\{-7, 1, 3, 4\right\}\\
\norm{x}_{\min}&=7
\end{aligned}
$$

可換表現の各因子は二つの直交射影の和である。両因子の射影が可換なのでその四つの積で切ると、各部分空間でも同じ四つのスカラーのどれかが作用する。従って最大ノルムにも同じ上界がある。

$$
\begin{aligned}
\norm{x}_{\max}&=7
\end{aligned}
$$

二項それぞれのノルムを足した値は三角不等式による上界である。実際のノルムでは、同じ点の対で二項が加わった後の値を調べる必要がある。
:::

:::example[符号の打ち消しがある場合]
有限可換代数 $\C^2$ の二つの因子を、二点上の乗算として表現する。$x=f\otimes g+h\otimes k$ の値を並べ、空間的ノルムを求める。

点の対 $(i,j)$ での値は $f_i g_j+h_i k_j$ である。積の二項を外積の行列で表示すると、座標ごとの寄与を分けて確認できる。

$$
\begin{aligned}
(f_i g_j)_{ij}&=\begin{pmatrix}2 & 2\\-2 & -2\end{pmatrix}\\
(h_i k_j)_{ij}&=\begin{pmatrix}2 & -2\\2 & -2\end{pmatrix}\\
(x_{ij})&=\begin{pmatrix}4 & 0\\0 & -4\end{pmatrix}
\end{aligned}
$$

表現先では $(i,j)$ に対応する基底ベクトルは固有ベクトルである。従って全ての値の絶対値の最大が作用素ノルムになる。

$$
\begin{aligned}
\sigma(x)&=\left\{-4, 0, 4\right\}\\
\norm{x}_{\min}&=4
\end{aligned}
$$

可換表現の各因子は二つの直交射影の和である。両因子の射影が可換なのでその四つの積で切ると、各部分空間でも同じ四つのスカラーのどれかが作用する。従って最大ノルムにも同じ上界がある。

$$
\begin{aligned}
\norm{x}_{\max}&=4
\end{aligned}
$$

二項それぞれのノルムを足した値は三角不等式による上界である。実際のノルムでは、同じ点の対で二項が加わった後の値を調べる必要がある。
:::

:::example[同じ第一因子をまとめる]
有限可換代数 $\C^2$ の二つの因子を、二点上の乗算として表現する。$x=f\otimes g+h\otimes k$ の値を並べ、空間的ノルムを求める。

点の対 $(i,j)$ での値は $f_i g_j+h_i k_j$ である。積の二項を外積の行列で表示すると、座標ごとの寄与を分けて確認できる。

$$
\begin{aligned}
(f_i g_j)_{ij}&=\begin{pmatrix}1 & 0\\3 & 0\end{pmatrix}\\
(h_i k_j)_{ij}&=\begin{pmatrix}0 & 1\\0 & 3\end{pmatrix}\\
(x_{ij})&=\begin{pmatrix}1 & 1\\3 & 3\end{pmatrix}
\end{aligned}
$$

表現先では $(i,j)$ に対応する基底ベクトルは固有ベクトルである。従って全ての値の絶対値の最大が作用素ノルムになる。

$$
\begin{aligned}
\sigma(x)&=\left\{1, 3\right\}\\
\norm{x}_{\min}&=3
\end{aligned}
$$

可換表現の各因子は二つの直交射影の和である。両因子の射影が可換なのでその四つの積で切ると、各部分空間でも同じ四つのスカラーのどれかが作用する。従って最大ノルムにも同じ上界がある。

$$
\begin{aligned}
\norm{x}_{\max}&=3
\end{aligned}
$$

二項それぞれのノルムを足した値は三角不等式による上界である。実際のノルムでは、同じ点の対で二項が加わった後の値を調べる必要がある。
:::

:::example[最大絶対値の符号が変わる]
有限可換代数 $\C^2$ の二つの因子を、二点上の乗算として表現する。$x=f\otimes g+h\otimes k$ の値を並べ、空間的ノルムを求める。

点の対 $(i,j)$ での値は $f_i g_j+h_i k_j$ である。積の二項を外積の行列で表示すると、座標ごとの寄与を分けて確認できる。

$$
\begin{aligned}
(f_i g_j)_{ij}&=\begin{pmatrix}-2 & 3\\-4 & 6\end{pmatrix}\\
(h_i k_j)_{ij}&=\begin{pmatrix}4 & 0\\1 & 0\end{pmatrix}\\
(x_{ij})&=\begin{pmatrix}2 & 3\\-3 & 6\end{pmatrix}
\end{aligned}
$$

表現先では $(i,j)$ に対応する基底ベクトルは固有ベクトルである。従って全ての値の絶対値の最大が作用素ノルムになる。

$$
\begin{aligned}
\sigma(x)&=\left\{-3, 2, 3, 6\right\}\\
\norm{x}_{\min}&=6
\end{aligned}
$$

可換表現の各因子は二つの直交射影の和である。両因子の射影が可換なのでその四つの積で切ると、各部分空間でも同じ四つのスカラーのどれかが作用する。従って最大ノルムにも同じ上界がある。

$$
\begin{aligned}
\norm{x}_{\max}&=6
\end{aligned}
$$

二項それぞれのノルムを足した値は三角不等式による上界である。実際のノルムでは、同じ点の対で二項が加わった後の値を調べる必要がある。
:::

### 二つのノルムが異なる自由群の例

:::example[左正則表現と右正則表現の積]
$G=\mathbb F_2=\langle a,b\rangle$ とし、$S=\{a,a^{-1},b,b^{-1}\}$ とする。$A$ は左正則作用素 $\lambda(g)$ が生成する単位的 C*-環である。$x=\sum_{g\in S}\lambda(g)\otimes\lambda(g)\in A\odot A$ について
$$
\norm x_{\min}=2\sqrt3,\qquad \norm x_{\max}=4
$$
である。

右正則表現 $\rho(g)\delta_h=\delta_{hg^{-1}}$ は左正則表現とユニタリ同値であり、$A$ の元の左作用と可換なもう一つの $A$ の表現を与える。実際、基底の反転 $\delta_h\mapsto\delta_{h^{-1}}$ は $\lambda(g)$ を $\rho(g)$ へ移す。$\lambda(g)\rho(g)$ は $\delta_e$ を固定するため、最大ノルムは $4$ 以上である。各項がユニタリなので上界も $4$ である。

空間的表現では $(s,t)\mapsto(s,s^{-1}t)$ の基底置換により、$\lambda(g)\otimes\lambda(g)$ は $\lambda(g)\otimes1$ へ移る。従って最小ノルムは Cayley グラフの隣接作用素 $T=\sum_{g\in S}\lambda(g)$ のノルムである。縮約語の長さを $|h|$ とし、正の重み $w(h)=3^{-|h|/2}$ を置く。非単位元には語長を一つ減らす隣点が一つ、一つ増やす隣点が三つある。従って隣点の重みの和は $2\sqrt3\,w(h)$ であり、単位元では $4/\sqrt3\le2\sqrt3$ である。各辺に対して $2|f(h)f(k)|\le|f(h)|^2w(k)/w(h)+|f(k)|^2w(h)/w(k)$ を使い、有限台の $f$ 上で和を取ると $|\inner{Tf}{f}|\le2\sqrt3\norm f^2$ を得る。$T$ は自己共役なのでこれがノルムの上界となる。

逆に $f_R=w\,1_{\{|h|\le R\}}$ を取る。語長 $j\ge1$ の点は $4\cdot3^{j-1}$ 個なので、$\norm{f_R}^2=1+4R/3$ である。深さ $j-1$ と $j$ の間の辺も $4\cdot3^{j-1}$ 本で、二次形式には各辺が二回現れる。従って $\inner{Tf_R}{f_R}=8R/\sqrt3$ である。その比が $2\sqrt3$ に収束するので、下界も一致する。
:::

### 最大テンソル積と商

:::proposition[最大テンソル積は短完全列を保つ]{#prop-oa-max-tensor-exact}
$I$ を C*-環 $A$ の閉両側イデアルとし、$B$ を C*-環とする。このとき自然な写像による列
$$
0\longrightarrow I\otimes_{\max}B
\longrightarrow A\otimes_{\max}B
\longrightarrow (A/I)\otimes_{\max}B
\longrightarrow 0
$$
は完全である。
:::

:::proof
$I,B$ の可換な表現の対 $\pi,\rho$ を考える。$\pi$ が零となる直交補を除けば $\pi$ は非退化である。$I$ の縮小近似単位元 $(e_\lambda)$ に対して $\pi(ae_\lambda)$ は $\norm a$ 以下であり、$\pi(i)\xi$ 上で $\pi(ai)\xi$ に収束する。従って強極限は $A$ の有界表現 $\widetilde\pi$ を定める。積と随伴は稠密な $\pi(I)H$ 上で確認できる。また各 $\pi(ae_\lambda)$ は $\rho(B)$ と可換なので、極限も可換である。こうして $I\odot B$ の最大ノルムと $A\odot B$ から来るノルムが一致し、第一の写像は単射となる。その像 $J$ は代数的な積と近似から閉イデアルである。

次に $(A\otimes_{\max}B)/J$ の忠実非退化表現を取り、商写像と合成する。最大テンソル積の表現は二つの可換な非退化表現の積である。例えば $B$ が単位的でなければ、$a\otimes f_\lambda$ を $B$ の近似単位元 $(f_\lambda)$ で強極限に送って、$A$ の表現を取り出す。極限は積表現の稠密な像に作用させれば確認できる。$I\odot B$ が消えるため、$B$ の表現の非退化性から $I$ の像は零である。従って $A$ の表現が $A/I$ の表現に降り、この商ノルムは $(A/I)\otimes_{\max}B$ のノルム以下となる。逆の評価は $A/I$ と $B$ の可換表現を $A$ に引き戻せば得られる。従って商環は $(A/I)\otimes_{\max}B$ と等長に同型であり、列は完全である。
:::

## 核型 C*-環と単射的 von Neumann 環

:::definition[核型性と単射性]{#def-oa-nuclearity}
任意の C*-環 $B$ に対して $\norm{\cdot}_{\min}=\norm{\cdot}_{\max}$ となる $A$ を :term[核型 C* 環]{en="nuclear C-star algebra" yomi="かくがたしーすたーかん"} という。単位を共有する任意の作用素系の包含 $S\subset T$ について、$S$ から von Neumann 環 $M$ への任意の UCP 写像が $T$ まで $M$ 値の UCP 写像として延長できるとき、$M$ は :term[単射的]{en="injective" yomi="たんしゃてき"} であるという。
:::

:::proposition[単射性と射影]{#prop-oa-injective-expectation}
$M\subset B(H)$ が単射的であることと、$B(H)$ から $M$ への条件付き期待値が存在することは同値である。
:::
:::proof
単射性を $\id_M$ に適用すると UCP 射影を得る。:ref[thm-oa-tomiyama] で期待値になる。逆に拡張すべき写像を $B(H)$ 値とみなし :ref[thm-oa-arveson-extension] を適用した後、期待値と合成すれば $M$ 値となる。
:::

単位的 C*-環 $A$ の完全正近似性とは、恒等写像を $A\xrightarrow{\phi}M_n\xrightarrow{\psi}A$ の UCP 合成で点ごとのノルム近似できることをいう。非単位的な環では完全正縮小写像を用いる。

:::proposition[完全正近似から核型性へ]{#prop-oa-cpap-nuclear}
完全正近似性を持つ C*-環 $A$ は核型である。
:::

:::proof
CP 縮小写像 $\phi:A\to D$ を一方の因子に作用させると最小ノルムで縮小的になる。忠実表現に合成した $\phi$ を Stinespring の形にし、もう一方の表現とテンソルを取れば、$\phi\odot\id$ のテンソル表現は縮小作用素 $V\otimes1$ による圧縮である。膨張した表現のノルムは忠実テンソル表現のノルム以下である。

最大ノルムについては、値域の $D$ と別の環 $B$ の可換な表現 $\alpha,\beta$ を取る。$\alpha\phi$ を $V^*\pi(\cdot)V$ に膨張させると、$\beta$ の作用を膨張空間に移せる。Stinespring の構成で $[a\otimes\xi]$ に $b$ を作用させる式を $[a\otimes\beta(b)\xi]$ とする。有限和の半内積を測る正行列 $[\alpha\phi(a_i^*a_j)]$ は $\beta(b)$ の対角増幅と可換である。この正行列の平方根を通して $\beta(b)$ を移動すれば、定めた作用のノルムは $\norm b$ 以下となり、零空間にも抵触しない。積と随伴は式で検査でき、得られた表現は $\pi$ と可換である。$V$ は $\beta$ とこの表現を絡める。単位元がない場合も :ref[cor-oa-nonunital-stinespring] の $V\xi$ の極限の式で同じことが分かる。従って $\sum_i\alpha(\phi(a_i))\beta(b_i)$ は可換な膨張表現の圧縮であり、ノルムは $\norm{\sum_i a_i\otimes b_i}_{\max}$ 以下である。最大ノルムで上限を取れば縮小性が出る。

$x=\sum_i a_i\otimes b_i$ に対し $x_n=(\psi\phi\odot\id)(x)$ は $x$ に最大ノルムで近づく。実際、差のノルムは $\sum_i\norm{\psi\phi(a_i)-a_i}\norm{b_i}$ 以下である。行列環の二つのノルムが一致することと両写像の縮小性から $\norm{x_n}_{\max}\le\norm{(\phi\odot\id)(x)}_{\min}\le\norm x_{\min}$ である。極限で $\norm x_{\max}\le\norm x_{\min}$ となる。
:::

逆に核型なら完全正近似性を持ち、核型 $A$ の各非退化表現の二重可換子は単射的である。前者は Choi–Effros の特徴付けであり、逆方向の有限行列を通る近似の抽出には二重双対上の半離散性が必要である。ここではその抽出の証明を扱わない。テンソルノルムによる核型性の議論は :cite[murphy-cstar]{loc="Chapter 6"}、von Neumann 環の単射性との関係は :cite[takesaki-structure] を参照する。

具体的に $C(X)$ がコンパクト Hausdorff $X$ でこの近似性を持つことは、有限個の関数の振幅が小さい開被覆と従属分割 $h_j$ を用いて分かる。$x_j$ を各開集合から取り、$f\mapsto\diag(f(x_j))$ と $b\mapsto\sum_j b_{jj}h_j$ を合成する。後者は各座標圧縮と正係数 $h_j$ の和なので CP であり、単位を保つ。各点で $\abs{\sum_j h_j(x)(f(x_j)-f(x))}$ が振幅の上界以下になる。

### 期待値の正規性は別に調べる

:::example[拡散 MASA の期待値は正規でない]
$N=L^\infty([0,1])\subset B(L^2([0,1]))$ とする。期待値 $E$ があるとして、階数一正元 $T=\theta_{\xi,\xi}$ を取る。積分 $\int_E|\xi|^2$ は絶対連続なので、$[0,1]$ を十分短い有限個の区間に分け、各部分で積分を $\eps$ 以下にできる。各区間の射影を $p_j$ とすれば $\norm{p_jTp_j}=\norm{p_j\xi}^2\le\eps$ である。両加群性より $p_jE(T)p_j=E(p_jTp_j)$ であり、$E(T)$ は $N$ にあるのでそのノルムは各ブロックのノルムの最大である。したがって $\norm{E(T)}\le\eps$、$\eps$ が任意なので $E(T)=0$ になる。有限階数射影もすべて消え、$1$ へ強収束する増大列でも上限が保たれない。従って期待値は正規でない。期待値自体の存在は可換環の単射性、または有限区間分割での状態平均の弱極限から得られる。
:::

:::example[コヒーレンスを完全に消す]
$d=\diag(1,-1)$、$p=(1+k)/2$、$q=1-p$ とし、$\Phi(a)=pa+qdad$ を考える。有限次元の表示で完全正性と非対角成分の変化を同時に調べる。

代数の表現を $\pi(a)=a\oplus a$ とし、$V\xi=(\sqrt p\,\xi,\sqrt q\,d\xi)$ と取る。$V^*V=(p+q)I=I$ なので、これは等長による圧縮になる。

$$
\begin{aligned}
(k,p,q)&=\left( 0, \  \frac{1}{2}, \  \frac{1}{2}\right)\\
V&=\begin{pmatrix}\frac{\sqrt{2}}{2} & 0\\0 & \frac{\sqrt{2}}{2}\\\frac{\sqrt{2}}{2} & 0\\0 & - \frac{\sqrt{2}}{2}\end{pmatrix}\\
a&=\begin{pmatrix}2 & 1 + i\\1 - i & 3\end{pmatrix}
\end{aligned}
$$

$V^*\pi(a)V$ の二項は、同じブロックからの寄与 $pa$ と $qdad$ である。$dad$ では対角は同じで、非対角だけ符号が変わる。従って非対角には $p-q=k$ が掛かる。

$$
\begin{aligned}
dad&=\begin{pmatrix}2 & -1 - i\\-1 + i & 3\end{pmatrix}\\
\Phi(a)&=\begin{pmatrix}2 & 0\\0 & 3\end{pmatrix}
\end{aligned}
$$

この表示を行列増幅しても $I_m\otimes V$ を使った圧縮になるため正性を保つ。端点 $p=0$ または $q=0$ では対応するブロックが不要になる。それ以外では $d$ と $I$ は一次独立なので、行列単位を適用して生成される空間は両ブロック全体である。

$$
\begin{aligned}
\dim H_{\min}&=4
\end{aligned}
$$

この場合、対角の確率を全て保っても非対角の内積は変化し得る。写像が状態を変えるかどうかを対角だけから判定できない。
:::

:::example[コヒーレンスを半分にする]
$d=\diag(1,-1)$、$p=(1+k)/2$、$q=1-p$ とし、$\Phi(a)=pa+qdad$ を考える。有限次元の表示で完全正性と非対角成分の変化を同時に調べる。

代数の表現を $\pi(a)=a\oplus a$ とし、$V\xi=(\sqrt p\,\xi,\sqrt q\,d\xi)$ と取る。$V^*V=(p+q)I=I$ なので、これは等長による圧縮になる。

$$
\begin{aligned}
(k,p,q)&=\left( \frac{1}{2}, \  \frac{3}{4}, \  \frac{1}{4}\right)\\
V&=\begin{pmatrix}\frac{\sqrt{3}}{2} & 0\\0 & \frac{\sqrt{3}}{2}\\\frac{1}{2} & 0\\0 & - \frac{1}{2}\end{pmatrix}\\
a&=\begin{pmatrix}2 & 1 + i\\1 - i & 3\end{pmatrix}
\end{aligned}
$$

$V^*\pi(a)V$ の二項は、同じブロックからの寄与 $pa$ と $qdad$ である。$dad$ では対角は同じで、非対角だけ符号が変わる。従って非対角には $p-q=k$ が掛かる。

$$
\begin{aligned}
dad&=\begin{pmatrix}2 & -1 - i\\-1 + i & 3\end{pmatrix}\\
\Phi(a)&=\begin{pmatrix}2 & \frac{1}{2} + \frac{i}{2}\\\frac{1}{2} - \frac{i}{2} & 3\end{pmatrix}
\end{aligned}
$$

この表示を行列増幅しても $I_m\otimes V$ を使った圧縮になるため正性を保つ。端点 $p=0$ または $q=0$ では対応するブロックが不要になる。それ以外では $d$ と $I$ は一次独立なので、行列単位を適用して生成される空間は両ブロック全体である。

$$
\begin{aligned}
\dim H_{\min}&=4
\end{aligned}
$$

この場合、対角の確率を全て保っても非対角の内積は変化し得る。写像が状態を変えるかどうかを対角だけから判定できない。
:::

:::example[コヒーレンスの符号を反転する]
$d=\diag(1,-1)$、$p=(1+k)/2$、$q=1-p$ とし、$\Phi(a)=pa+qdad$ を考える。有限次元の表示で完全正性と非対角成分の変化を同時に調べる。

代数の表現を $\pi(a)=a\oplus a$ とし、$V\xi=(\sqrt p\,\xi,\sqrt q\,d\xi)$ と取る。$V^*V=(p+q)I=I$ なので、これは等長による圧縮になる。

$$
\begin{aligned}
(k,p,q)&=\left( - \frac{1}{2}, \  \frac{1}{4}, \  \frac{3}{4}\right)\\
V&=\begin{pmatrix}\frac{1}{2} & 0\\0 & \frac{1}{2}\\\frac{\sqrt{3}}{2} & 0\\0 & - \frac{\sqrt{3}}{2}\end{pmatrix}\\
a&=\begin{pmatrix}2 & 1 + i\\1 - i & 3\end{pmatrix}
\end{aligned}
$$

$V^*\pi(a)V$ の二項は、同じブロックからの寄与 $pa$ と $qdad$ である。$dad$ では対角は同じで、非対角だけ符号が変わる。従って非対角には $p-q=k$ が掛かる。

$$
\begin{aligned}
dad&=\begin{pmatrix}2 & -1 - i\\-1 + i & 3\end{pmatrix}\\
\Phi(a)&=\begin{pmatrix}2 & - \frac{1}{2} - \frac{i}{2}\\- \frac{1}{2} + \frac{i}{2} & 3\end{pmatrix}
\end{aligned}
$$

この表示を行列増幅しても $I_m\otimes V$ を使った圧縮になるため正性を保つ。端点 $p=0$ または $q=0$ では対応するブロックが不要になる。それ以外では $d$ と $I$ は一次独立なので、行列単位を適用して生成される空間は両ブロック全体である。

$$
\begin{aligned}
\dim H_{\min}&=4
\end{aligned}
$$

この場合、対角の確率を全て保っても非対角の内積は変化し得る。写像が状態を変えるかどうかを対角だけから判定できない。
:::

:::example[恒等写像に一致する]
$d=\diag(1,-1)$、$p=(1+k)/2$、$q=1-p$ とし、$\Phi(a)=pa+qdad$ を考える。有限次元の表示で完全正性と非対角成分の変化を同時に調べる。

代数の表現を $\pi(a)=a\oplus a$ とし、$V\xi=(\sqrt p\,\xi,\sqrt q\,d\xi)$ と取る。$V^*V=(p+q)I=I$ なので、これは等長による圧縮になる。

$$
\begin{aligned}
(k,p,q)&=\left( 1, \  1, \  0\right)\\
V&=\begin{pmatrix}1 & 0\\0 & 1\\0 & 0\\0 & 0\end{pmatrix}\\
a&=\begin{pmatrix}2 & 1 + i\\1 - i & 3\end{pmatrix}
\end{aligned}
$$

$V^*\pi(a)V$ の二項は、同じブロックからの寄与 $pa$ と $qdad$ である。$dad$ では対角は同じで、非対角だけ符号が変わる。従って非対角には $p-q=k$ が掛かる。

$$
\begin{aligned}
dad&=\begin{pmatrix}2 & -1 - i\\-1 + i & 3\end{pmatrix}\\
\Phi(a)&=\begin{pmatrix}2 & 1 + i\\1 - i & 3\end{pmatrix}
\end{aligned}
$$

この表示を行列増幅しても $I_m\otimes V$ を使った圧縮になるため正性を保つ。端点 $p=0$ または $q=0$ では対応するブロックが不要になる。それ以外では $d$ と $I$ は一次独立なので、行列単位を適用して生成される空間は両ブロック全体である。

$$
\begin{aligned}
\dim H_{\min}&=2
\end{aligned}
$$

この場合、対角の確率を全て保っても非対角の内積は変化し得る。写像が状態を変えるかどうかを対角だけから判定できない。
:::

:::example[対角ユニタリ共役に一致する]
$d=\diag(1,-1)$、$p=(1+k)/2$、$q=1-p$ とし、$\Phi(a)=pa+qdad$ を考える。有限次元の表示で完全正性と非対角成分の変化を同時に調べる。

代数の表現を $\pi(a)=a\oplus a$ とし、$V\xi=(\sqrt p\,\xi,\sqrt q\,d\xi)$ と取る。$V^*V=(p+q)I=I$ なので、これは等長による圧縮になる。

$$
\begin{aligned}
(k,p,q)&=\left( -1, \  0, \  1\right)\\
V&=\begin{pmatrix}0 & 0\\0 & 0\\1 & 0\\0 & -1\end{pmatrix}\\
a&=\begin{pmatrix}2 & 1 + i\\1 - i & 3\end{pmatrix}
\end{aligned}
$$

$V^*\pi(a)V$ の二項は、同じブロックからの寄与 $pa$ と $qdad$ である。$dad$ では対角は同じで、非対角だけ符号が変わる。従って非対角には $p-q=k$ が掛かる。

$$
\begin{aligned}
dad&=\begin{pmatrix}2 & -1 - i\\-1 + i & 3\end{pmatrix}\\
\Phi(a)&=\begin{pmatrix}2 & -1 - i\\-1 + i & 3\end{pmatrix}
\end{aligned}
$$

この表示を行列増幅しても $I_m\otimes V$ を使った圧縮になるため正性を保つ。端点 $p=0$ または $q=0$ では対応するブロックが不要になる。それ以外では $d$ と $I$ は一次独立なので、行列単位を適用して生成される空間は両ブロック全体である。

$$
\begin{aligned}
\dim H_{\min}&=2
\end{aligned}
$$

この場合、対角の確率を全て保っても非対角の内積は変化し得る。写像が状態を変えるかどうかを対角だけから判定できない。
:::

## Connes の分類定理（概観）

:::theorem[単射因子環と超有限性]{#thm-oa-connes-hyperfinite}
可分な前双対を持つ因子環で、$\mathrm{III}_1$ 型を除くと、単射性と超有限性が同値である。単射的な $\mathrm{II}_1$ 因子環は :ref[prop-oa-hyperfinite-two] の $R$ と同型であり、単射的な $\mathrm{II}_\infty$ 因子環は $R\overline\otimes B(\ell^2(\N))$ と同型である。また $0<\lambda<1$ ごとに単射的な $\mathrm{III}_\lambda$ 因子環は同型を除いて一つであり、$\mathrm{III}_0$ 型の単射因子環は重みの流れの共役類で区別される。
:::
:::proof
超有限環は、有限次元部分環の増大合併のノルム閉包 $A$ を使えば単射性を示せる。実際、有限集合を一つの有限次元部分環 $F$ に近似し、単射的な $F$ への UCP 射影と包含が与える合成写像で恒等写像を近似できる。$F$ は全行列環の有限直和なので、ブロック対角に埋め込み、対角ブロックだけを取る UCP 射影を使えば全行列環を通る近似になる。よって $A$ は核型であり、その二重可換子は単射的である。逆方向は有限次元作用素系の CP 近似を有限次元部分環での近似に替える段階で、半離散性と因子環の構造論を要する。III 型の同型を決めるには、さらにモジュラー群から重みの流れを構成する必要がある。これらは本章では概観にとどめる。関連する超有限模型については :cite[takesaki-structure] を参照する。
:::

後に Haagerup が可分な前双対を持つ単射的な $\mathrm{III}_1$ 因子環も同型を除いて一つであることを示した。従ってこの型でも単射性と超有限性は同値である。III 型の区別に使うスペクトルと重みの流れについては [Tomita–Takesaki 理論](/math/operator-algebras/tomita-takesaki/) を参照する。

有限次元部分環への近似から単射性が現れる仕組みは、二進行列模型 $R$ で直接見える。$R$ 上の規格化トレースを $\tau$、有限段を $A_n=M_d$ とし、その行列単位を $e_{ij}$ とする。$E_n(x)=\sum_{i,j}d\,\tau(e_{ji}x)e_{ij}$ は $L^2$ で $A_n$ への直交射影である。$E_n$ は UCP である。実際、$\sum_{i,j}\bar c_i c_jd\,\tau(e_{ji}x)$ は $d\,\tau(v^*xv)$（$v=\sum_j c_je_{j1}$）であり、任意の増幅でもベクトルごとに $v$ を選ぶと正性を確認できる。右作用との可換性と一様有界性から $E_n(x)\to x$ は強収束する。各 $E_n$ を Arveson の定理で $B(H)$ に拡張し、点ごとの WOT クラスター点を取ると $R$ に値を取る UCP 射影となる。

### 有限段近似で確認できる向き

$M$ が有限次元部分環の増大合併の σ 弱閉包であり、忠実トレースを持つ場合を考える。各有限次元部分環 $F$ にはトレースを保つ期待値を式で作れる。行列単位 $e_{ij}^{(k)}$ のある各成分について、$x$ の対応する係数を $\tau(e_{ji}^{(k)}x)/\tau(e_{jj}^{(k)})$ と置く。トレースの循環性により各対角射影のトレースは成分内で等しい。この写像は $F$ を固定し、両加群的である。正性は各ブロックの二次形式を $\tau(v^*xv)$ として書けば従う。さらに行列増幅にも同じ方法を使え、完全正になる。$L^2$ 内積で直交射影となるため、有限段の合併がノルムで稠密な $L^2$ 部分空間ではこれらの射影が恒等に収束する。この構成は、有限次元近似が既にある場合にその近似を完全正写像として表す。

:::example[有限段の外をトレースで平均する]
超有限因子環の二つの行列段 $M_2\subset M_2\otimes M_2$ を取り、$E=\id\otimes(\tr/2)$ で小さい段へ戻す。二乗平均での近似を具体的に計算する。

$x=E_{11}\otimes D+(E_{12}+E_{21})\otimes F$ とし、$D,F$ は第二因子の対角行列である。$E$ は最初の因子を残し、第二因子をその規格化トレースで置き換える。

$$
\begin{aligned}
D&=\begin{pmatrix}2 & 0\\0 & 1\end{pmatrix}\\
F&=\begin{pmatrix}1 & 0\\0 & 2\end{pmatrix}\\
x&=\begin{pmatrix}2 & 0 & 1 & 0\\0 & 1 & 0 & 2\\1 & 0 & 0 & 0\\0 & 2 & 0 & 0\end{pmatrix}\\
E(x)&=\begin{pmatrix}\frac{3}{2} & \frac{3}{2}\\\frac{3}{2} & 0\end{pmatrix}
\end{aligned}
$$

差を比較するときには $E(x)$ をもう一度大きい段に埋め込む。$x-E(x)\otimes I$ の各第二因子は規格化トレースが零になる。したがって小さい段の任意の $a\otimes I$ とトレース内積で直交する。

$$
\begin{aligned}
x-E(x)\otimes I&=\begin{pmatrix}\frac{1}{2} & 0 & - \frac{1}{2} & 0\\0 & - \frac{1}{2} & 0 & \frac{1}{2}\\- \frac{1}{2} & 0 & 0 & 0\\0 & \frac{1}{2} & 0 & 0\end{pmatrix}
\end{aligned}
$$

規格化トレースでの二乗誤差は次の通りである。第二因子がスカラーに近いほど誤差は小さくなり、両方がスカラーなら零になる。

$$
\begin{aligned}
\norm{x-E(x)\otimes I}_2^2&=\frac{3}{8}
\end{aligned}
$$

有限段を増やして同じ期待値を作ると、各固定の有限段に入る元は十分後から正確に復元される。増大合併の閉包から出発する場合は、これが超有限性による近似の具体的な向きである。
:::

:::example[一方のテンソル項の平均が縮む]
超有限因子環の二つの行列段 $M_2\subset M_2\otimes M_2$ を取り、$E=\id\otimes(\tr/2)$ で小さい段へ戻す。二乗平均での近似を具体的に計算する。

$x=E_{11}\otimes D+(E_{12}+E_{21})\otimes F$ とし、$D,F$ は第二因子の対角行列である。$E$ は最初の因子を残し、第二因子をその規格化トレースで置き換える。

$$
\begin{aligned}
D&=\begin{pmatrix}3 & 0\\0 & 2\end{pmatrix}\\
F&=\begin{pmatrix}0 & 0\\0 & 1\end{pmatrix}\\
x&=\begin{pmatrix}3 & 0 & 0 & 0\\0 & 2 & 0 & 1\\0 & 0 & 0 & 0\\0 & 1 & 0 & 0\end{pmatrix}\\
E(x)&=\begin{pmatrix}\frac{5}{2} & \frac{1}{2}\\\frac{1}{2} & 0\end{pmatrix}
\end{aligned}
$$

差を比較するときには $E(x)$ をもう一度大きい段に埋め込む。$x-E(x)\otimes I$ の各第二因子は規格化トレースが零になる。したがって小さい段の任意の $a\otimes I$ とトレース内積で直交する。

$$
\begin{aligned}
x-E(x)\otimes I&=\begin{pmatrix}\frac{1}{2} & 0 & - \frac{1}{2} & 0\\0 & - \frac{1}{2} & 0 & \frac{1}{2}\\- \frac{1}{2} & 0 & 0 & 0\\0 & \frac{1}{2} & 0 & 0\end{pmatrix}
\end{aligned}
$$

規格化トレースでの二乗誤差は次の通りである。第二因子がスカラーに近いほど誤差は小さくなり、両方がスカラーなら零になる。

$$
\begin{aligned}
\norm{x-E(x)\otimes I}_2^2&=\frac{3}{8}
\end{aligned}
$$

有限段を増やして同じ期待値を作ると、各固定の有限段に入る元は十分後から正確に復元される。増大合併の閉包から出発する場合は、これが超有限性による近似の具体的な向きである。
:::

:::example[両方の項がそのまま残る]
超有限因子環の二つの行列段 $M_2\subset M_2\otimes M_2$ を取り、$E=\id\otimes(\tr/2)$ で小さい段へ戻す。二乗平均での近似を具体的に計算する。

$x=E_{11}\otimes D+(E_{12}+E_{21})\otimes F$ とし、$D,F$ は第二因子の対角行列である。$E$ は最初の因子を残し、第二因子をその規格化トレースで置き換える。

$$
\begin{aligned}
D&=\begin{pmatrix}2 & 0\\0 & 2\end{pmatrix}\\
F&=\begin{pmatrix}2 & 0\\0 & 2\end{pmatrix}\\
x&=\begin{pmatrix}2 & 0 & 2 & 0\\0 & 2 & 0 & 2\\2 & 0 & 0 & 0\\0 & 2 & 0 & 0\end{pmatrix}\\
E(x)&=\begin{pmatrix}2 & 2\\2 & 0\end{pmatrix}
\end{aligned}
$$

差を比較するときには $E(x)$ をもう一度大きい段に埋め込む。$x-E(x)\otimes I$ の各第二因子は規格化トレースが零になる。したがって小さい段の任意の $a\otimes I$ とトレース内積で直交する。

$$
\begin{aligned}
x-E(x)\otimes I&=\begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{pmatrix}
\end{aligned}
$$

規格化トレースでの二乗誤差は次の通りである。第二因子がスカラーに近いほど誤差は小さくなり、両方がスカラーなら零になる。

$$
\begin{aligned}
\norm{x-E(x)\otimes I}_2^2&=0
\end{aligned}
$$

有限段を増やして同じ期待値を作ると、各固定の有限段に入る元は十分後から正確に復元される。増大合併の閉包から出発する場合は、これが超有限性による近似の具体的な向きである。
:::

:::example[符号の異なる二つのトレース]
超有限因子環の二つの行列段 $M_2\subset M_2\otimes M_2$ を取り、$E=\id\otimes(\tr/2)$ で小さい段へ戻す。二乗平均での近似を具体的に計算する。

$x=E_{11}\otimes D+(E_{12}+E_{21})\otimes F$ とし、$D,F$ は第二因子の対角行列である。$E$ は最初の因子を残し、第二因子をその規格化トレースで置き換える。

$$
\begin{aligned}
D&=\begin{pmatrix}1 & 0\\0 & -2\end{pmatrix}\\
F&=\begin{pmatrix}1 & 0\\0 & -1\end{pmatrix}\\
x&=\begin{pmatrix}1 & 0 & 1 & 0\\0 & -2 & 0 & -1\\1 & 0 & 0 & 0\\0 & -1 & 0 & 0\end{pmatrix}\\
E(x)&=\begin{pmatrix}- \frac{1}{2} & 0\\0 & 0\end{pmatrix}
\end{aligned}
$$

差を比較するときには $E(x)$ をもう一度大きい段に埋め込む。$x-E(x)\otimes I$ の各第二因子は規格化トレースが零になる。したがって小さい段の任意の $a\otimes I$ とトレース内積で直交する。

$$
\begin{aligned}
x-E(x)\otimes I&=\begin{pmatrix}\frac{3}{2} & 0 & 1 & 0\\0 & - \frac{3}{2} & 0 & -1\\1 & 0 & 0 & 0\\0 & -1 & 0 & 0\end{pmatrix}
\end{aligned}
$$

規格化トレースでの二乗誤差は次の通りである。第二因子がスカラーに近いほど誤差は小さくなり、両方がスカラーなら零になる。

$$
\begin{aligned}
\norm{x-E(x)\otimes I}_2^2&=\frac{17}{8}
\end{aligned}
$$

有限段を増やして同じ期待値を作ると、各固定の有限段に入る元は十分後から正確に復元される。増大合併の閉包から出発する場合は、これが超有限性による近似の具体的な向きである。
:::

## 演習問題

::::exercise{level=basic}
$\norm{\xi}=2$、$\norm{\eta}=3$ のとき $\norm{\xi\otimes\eta}$ を求めよ。
:::solution
単純テンソルの内積は内積の積なので、ノルムの二乗は $4\cdot9=36$ であり、ノルムは $6$ となる。
:::
::::
::::exercise{level=basic}
$\diag(1,-2)\otimes\diag(3,0)$ の固有値とノルムを求めよ。
:::solution
対角基底で値は $3,0,-6,0$ である。絶対値の最大は $6$ である。
:::
::::
::::exercise{level=basic}
$M_2\otimes_{\min}M_4$ の複素次元を求めよ。
:::solution
行列単位のテンソルが $M_8$ の行列単位に対応する。したがって次元は $8^2=64$ である。
:::
::::
::::exercise{level=basic}
$M_2\otimes1$ が $\C^2\otimes\C^3$ に作用するとき可換子を求めよ。
:::solution
行列単位と可換なら非対角ブロックは零で、二つの対角ブロックは等しい。よって $1\otimes M_3$ である。
:::
::::
::::exercise{level=advanced}
$M_n\odot A$ について空間ノルムと普遍ノルムが一致する理由を述べよ。
:::solution
$M_n$ の非退化表現は多重度付きの通常表現である。可換な $A$ の表現は可換子の形から多重度空間だけに作用する。最大ノルムを作る各表現がテンソル表現となり、:ref[prop-oa-tensor-norms] によりそのノルムは最小ノルム以下になる。
:::
::::
::::exercise{level=advanced}
二つの von Neumann 環が因子環なら空間的テンソル積も因子環であることを示せ。
:::solution
中心元 $T$ は $M\overline\otimes N$ とその可換子の双方に属する。$K$ 側でスライスすると $H$ 側は $M$（$M'\otimes1$ と可換）かつ $M'$（$M\otimes1$ と可換）なのでスカラーである。行列係数から $T=1\otimes S$ となる。同様に $S$ は $Z(N)=\C1$ に属するため $T$ もスカラーである。
:::
::::
::::exercise{level=advanced}
有限個の連続関数の完全正近似について、近似誤差が共通に $\eps$ 以下となる被覆を作る理由を説明せよ。
:::solution
各点 $x$ に対し全関数について $\abs{f_i(y)-f_i(x)}<\eps/2$ となる共通近傍を取る。有限部分被覆を選ぶと、同じ近傍の任意の二点の値の差は $\eps$ 未満である。従属分割の非負係数の和が $1$ なので、この誤差が合成写像にも保たれる。
:::
::::
::::exercise{level=advanced}
単射的 $M$ に UCP 射影を持つ von Neumann 部分環 $N$ は単射的であることを示せ。
:::solution
$N$ 値の UCP 写像をまず $M$ 値とみなして任意の大きな作用素系へ延長する。これと射影を合成すれば元の写像を保つ $N$ 値の UCP 拡張になる。
:::
::::

## 確認テスト

::::quiz{answer="8" tol="1e-6"}
$\diag(2,-1)\otimes \diag(3,4)$ のノルムを答えよ。
:::explanation
単純テンソルの作用素ノルムはノルムの積であり、$2\cdot4=8$ となる。
:::
::::
::::quiz
最大テンソルノルムを計算するときの表現の条件はどれか。

- [x] 同じ空間上で像が可換である。
- [ ] 各像が必ず可換環である。
- [ ] 両因子が必ず同じ環である。

:::explanation
異なる因子から来る作用素同士が可換なため、$a\otimes b$ の像を両者の積にできる。各因子の内部まで可換である必要はない。例えば全行列環の二つのテンソル因子は各々非可換である。
:::
::::
::::quiz{answer="6" tol="1e-6"}
$M_2\otimes M_3$ の全行列環としての行列サイズを答えよ。
:::explanation
$\C^2\otimes\C^3$ の次元は $2\cdot3=6$ なので $M_6$ である。
:::
::::
::::quiz
どの写像が単射性の具体的な判定になるか。

- [x] $B(H)$ から $M$ への条件付き期待値
- [ ] $M$ からスカラーへの任意の状態
- [ ] $M$ を $B(H)$ に入れる任意の忠実表現

:::explanation
Arveson 拡張で得る $B(H)$ 値の写像を $M$ 値に戻せる条件付き期待値が、単射性を与える。忠実表現やスカラー値状態の存在は一般の非零 C*-環にも成り立つため、$M$ 値の拡張を与えない。
:::
::::

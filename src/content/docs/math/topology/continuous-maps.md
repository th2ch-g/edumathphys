---
title: 連続写像と同相
description: 逆像・閉包・内部・近傍・ネットによる連続性の同値な特徴付け、開写像と閉写像、同相写像と位相的性質、位相の強弱と始位相・終位相を扱う。
notation:
  - symbol: 'f^{-1}(V)'
    meaning: 集合 V の写像 f による逆像
  - symbol: 'f(A)'
    meaning: 集合 A の写像 f による像
  - symbol: 'f|_A'
    meaning: 写像 f の部分集合 A への制限
  - symbol: 'X\cong Y'
    meaning: 空間 X と Y が同相であること
  - symbol: '\chi_A'
    meaning: 集合 A の特性関数（A 上で 1、外で 0）
  - symbol: 'S^n'
    meaning: ℝⁿ⁺¹ の単位球面
  - symbol: 'B^n,\ D^n'
    meaning: ℝⁿ の開単位球と閉単位球（Euclid ノルム）
  - symbol: '\mathcal T_1\subset\mathcal T_2'
    meaning: 位相 𝒯₂ が 𝒯₁ より細かいこと
  - symbol: 'X_s'
    meaning: 空間 X の点列開集合全体を位相とする空間
---

## この章の目標

- 連続写像を開集合の逆像で定義し、閉集合・閉包・内部・近傍・開基・ネット・点列による同値な言い換えを証明する。
- 実数値連続関数の和・積・商・最大・最小が連続であることを位相空間上で示し、貼り合わせの補題と一様極限の定理を使いこなす。
- 開写像・閉写像を連続性と区別し、像についての判定法と典型的な例・反例を理解する。
- 同相写像を明示的な式で構成し（区間、球と Euclid 空間、立体射影、球と立方体）、位相的性質とそうでない性質を区別する。
- 同じ集合上の位相を比較し、写像族から始位相・終位相を構成してその普遍性を証明する。

[距離空間](/math/topology/metric-spaces/) では、連続性を $\eps$-$\delta$ 論法で定義し、それが「開集合の逆像が開であること」と同値であることを示した（:ref[thm-top-metric-continuity-tests]）。後者は距離を用いずに述べられるので、[位相空間](/math/topology/topological-spaces/) の上でもそのまま定義として採用できる。本章ではこの定義から出発し、連続性をさまざまな形に言い換える。次に、像の振る舞いに関する開写像・閉写像を調べ、二つの空間を「位相的に同じ」とみなす同相の概念を導入する。最後に、写像を連続にするという要請から逆に位相を定める方法（始位相・終位相）を述べる。これは次章 [部分空間・積空間・商空間](/math/topology/constructions/) の構成の共通の土台となる。

本章を通じて $X,Y,Z$ は位相空間を表し、部分集合 $A\subset X$ に対する閉包を $\overline A$、内部を $A^\circ$ と書く（:ref[def-top-closure]）。距離空間には常に距離から定まる位相を入れる。

## 連続写像の同値な特徴付け

### 連続性の定義

距離空間における連続性は「入力を十分近づければ出力がいくらでも近くなる」ことであった。「出力が $f(a)$ の $\eps$ 近く」という条件は $f(a)$ の開球で、「入力が $a$ の $\delta$ 近く」という条件は $a$ の開球で表される。開球を一般の近傍に置き換えれば、距離を持たない空間でも同じ考え方ができる。

:::definition[位相空間の連続写像]{#def-top-continuous}
$f:X\to Y$ を位相空間の間の写像とする。

1. $x\in X$ とする。$f(x)$ の任意の近傍 $N$ に対し逆像 $f^{-1}(N)$ が $x$ の近傍であるとき、$f$ は $x$ で:term[連続]{en="continuous" yomi="れんぞく"}であるという。
2. 任意の開集合 $V\subset Y$ に対して $f^{-1}(V)=\{x\in X:f(x)\in V\}$ が $X$ の開集合であるとき、$f$ を:term[連続写像]{en="continuous map" yomi="れんぞくしゃぞう"}という。

$X$ から $Y$ への連続写像全体を $C(X,Y)$ と書き、$C(X)=C(X,\R)$ とおく。
:::

点 $x$ での連続性は、近傍の定義（:ref[def-top-neighborhood]）を展開すると次のように言い換えられる：$f(x)$ を含む任意の開集合 $V$ に対し、$x$ を含む開集合 $U$ で $f(U)\subset V$ となるものが存在する。実際、$f^{-1}(V)$ が $x$ の近傍であることは、$x\in U\subset f^{-1}(V)$ となる開集合 $U$ があること、すなわち $f(U)\subset V$ となる開集合 $U\ni x$ があることにほかならない。また $f(x)$ の任意の近傍は $f(x)$ を含む開集合を含み、逆像は包含を保つので、開近傍についてだけ確かめれば十分である。

逆像は写像の単射性や全射性を必要としない操作であり、和・共通部分・補集合と交換する：

$$
f^{-1}\Bigl(\bigcup_i V_i\Bigr)=\bigcup_i f^{-1}(V_i),\qquad
f^{-1}\Bigl(\bigcap_i V_i\Bigr)=\bigcap_i f^{-1}(V_i),\qquad
f^{-1}(Y\setminus V)=X\setminus f^{-1}(V).
$$

連続性の議論はこの三つの等式を繰り返し使う。像 $f(A)$ は和とは交換するが、共通部分や補集合とは一般に交換しない。これが連続性（逆像の性質）と後で扱う開写像（像の性質）の扱いやすさの違いになる。

:::proposition[距離空間の連続性との一致]{#prop-top-continuity-metric-agree}
$(X,d)$, $(Y,\rho)$ を距離空間とし、それぞれに距離から定まる位相を入れる。写像 $f:X\to Y$ と $a\in X$ について、$f$ が :ref[def-top-continuous] の意味で $a$ で連続であることと、:ref[def-top-metric-continuity] の $\eps$-$\delta$ の意味で $a$ で連続であることは同値である。したがって二つの意味の「連続写像」も一致する。
:::

:::proof
距離空間では、$N$ が $a$ の近傍であることと、ある $r>0$ で $B_d(a,r)\subset N$ となることは同値である。実際、$N$ が開集合 $U\ni a$ を含めば、開集合の定義（:ref[def-top-metric-open]）により $B_d(a,r)\subset U$ となる $r>0$ がある。逆に $B_d(a,r)\subset N$ なら、開球は開集合である（:ref[prop-top-ball-open]）から $N$ は近傍である。

位相的な意味で連続とする。$\eps>0$ に対し $B_\rho(f(a),\eps)$ は $f(a)$ の近傍なので、その逆像は $a$ の近傍であり、ある $\delta>0$ で $B_d(a,\delta)\subset f^{-1}(B_\rho(f(a),\eps))$ となる。これは $d(x,a)<\delta\Rightarrow\rho(f(x),f(a))<\eps$ ということである。

逆に $\eps$-$\delta$ の意味で連続とし、$N$ を $f(a)$ の近傍とする。$B_\rho(f(a),\eps)\subset N$ となる $\eps>0$ をとり、対応する $\delta$ をとれば $B_d(a,\delta)\subset f^{-1}(B_\rho(f(a),\eps))\subset f^{-1}(N)$ なので、$f^{-1}(N)$ は $a$ の近傍である。

写像全体の連続性については、次の :ref[thm-top-continuity-tests] で「各点で連続」と「開集合の逆像が開」が同値であることを示すので、上の点ごとの一致から従う（:ref[thm-top-metric-continuity-tests] とも整合する）。
:::

:::example[位相の両極端と連続性]{#exm-top-continuity-extremes}
1. **定数写像・恒等写像**：定数写像 $c:X\to Y$, $c(x)=y_0$ による開集合 $V$ の逆像は、$y_0\in V$ なら $X$、そうでなければ $\varnothing$ であり、いずれも開である。恒等写像 $\id_X$ による逆像は集合そのものである。したがって両者は連続である。
2. **離散空間からの写像**：$X$ が離散位相をもてば、$X$ のすべての部分集合が開なので、任意の写像 $f:X\to Y$ は連続である。
3. **密着空間への写像**：$Y$ が密着位相をもてば、$Y$ の開集合は $\varnothing,Y$ だけであり、その逆像は $\varnothing,X$ なので、任意の写像 $f:X\to Y$ は連続である。
4. **離散空間への写像**：$Y$ が離散位相をもつとき、$f:X\to Y$ が連続であるための必要十分条件は、すべての $y\in Y$ について $f^{-1}(\{y\})$ が $X$ で開であることである。必要性は一点集合が開であることから、十分性は $f^{-1}(V)=\bigcup_{y\in V}f^{-1}(\{y\})$ から従う。例えば $X=\R$（通常の位相）から二点離散空間 $\{0,1\}$ への連続写像は定数写像に限ることが [連結性](/math/topology/connectedness/) の章でわかる。
:::

:::example[特性関数で開集合を読み取る]{#exm-top-sierpinski-characteristic}
$S=\{0,1\}$ に開集合 $\varnothing,\{1\},S$ を与えた Sierpiński 空間を考える。$A\subset X$ の特性関数 $\chi_A:X\to S$（$A$ 上で $1$、外で $0$）が連続であるための条件を調べる。開集合の逆像は

$$
\chi_A^{-1}(\varnothing)=\varnothing,\qquad \chi_A^{-1}(\{1\})=A,\qquad \chi_A^{-1}(S)=X
$$

なので、$\chi_A$ が連続であることは $A$ が開であることと同値である。こうして $X$ の開集合全体は $C(X,S)$ と一対一に対応する。一方、値域 $\{0,1\}$ に離散位相を入れると、$\{0\}$ の逆像 $X\setminus A$ も開でなければならず、$\chi_A$ が連続であることは $A$ が開かつ閉であることと同値になる。同じ写像でも、値域の位相によって連続性の意味が変わる。
:::

:::example[有限空間の間の写像をすべて調べる]
Sierpiński 空間 $S$ から自身への写像は 4 つある。二つの定数写像と恒等写像は連続である。残る交換写像 $\sigma(0)=1,\ \sigma(1)=0$ では $\sigma^{-1}(\{1\})=\{0\}$ が開でないので連続でない。したがって $C(S,S)$ はちょうど 3 個の元をもつ。定義域を離散な二点空間に取り替えれば、例 :ref[exm-top-continuity-extremes] の 2 により 4 個すべてが連続になる。
:::

### 連続性の同値な特徴付け

開集合ではなく閉集合・閉包・内部・近傍で連続性を表すと、状況に応じて検査しやすい形を選べる。

:::theorem[連続性の特徴付け]{#thm-top-continuity-tests}
写像 $f:X\to Y$ について次の条件は同値である。

1. $f$ は連続である（$Y$ の開集合の逆像は $X$ で開）。
2. $f$ は $X$ の各点で連続である。
3. $Y$ の閉集合の逆像は $X$ で閉である。
4. 任意の $A\subset X$ について $f(\overline A)\subset\overline{f(A)}$ である。
5. 任意の $B\subset Y$ について $\overline{f^{-1}(B)}\subset f^{-1}(\overline B)$ である。
6. 任意の $B\subset Y$ について $f^{-1}(B^\circ)\subset\bigl(f^{-1}(B)\bigr)^\circ$ である。
:::

:::proof
**1 ⇒ 2.** $x\in X$ とし、$N$ を $f(x)$ の近傍とする。$f(x)\in V\subset N$ となる開集合 $V$ をとると、$f^{-1}(V)$ は $x$ を含む開集合で $f^{-1}(N)$ に含まれる。よって $f^{-1}(N)$ は $x$ の近傍である。

**2 ⇒ 1.** $V\subset Y$ を開集合とし、$x\in f^{-1}(V)$ とする。$V$ は $f(x)$ の近傍なので $f^{-1}(V)$ は $x$ の近傍であり、$x\in U_x\subset f^{-1}(V)$ となる開集合 $U_x$ がある。すると $f^{-1}(V)=\bigcup_{x\in f^{-1}(V)}U_x$ は開集合の和であり、開である。

**1 ⇔ 3.** $F\subset Y$ が閉であることは $Y\setminus F$ が開であることであり、$f^{-1}(Y\setminus F)=X\setminus f^{-1}(F)$ である。したがって「開集合の逆像が開」と「閉集合の逆像が閉」は補集合を取ることで互いに移り合う。

**3 ⇒ 4.** $A\subset X$ とする。$\overline{f(A)}$ は閉なので、3 により $f^{-1}(\overline{f(A)})$ は閉であり、$A\subset f^{-1}(f(A))\subset f^{-1}(\overline{f(A)})$ である。閉包は $A$ を含む最小の閉集合である（:ref[prop-top-closure-laws]）から $\overline A\subset f^{-1}(\overline{f(A)})$、すなわち $f(\overline A)\subset\overline{f(A)}$ を得る。

**4 ⇒ 5.** $B\subset Y$ とし、4 を $A=f^{-1}(B)$ に適用する。$f(f^{-1}(B))\subset B$ と閉包の単調性から

$$
f\bigl(\overline{f^{-1}(B)}\bigr)\subset\overline{f(f^{-1}(B))}\subset\overline B
$$

となり、両辺の逆像をとれば $\overline{f^{-1}(B)}\subset f^{-1}(\overline B)$ である。

**5 ⇒ 3.** $F\subset Y$ を閉集合とすると $\overline F=F$ なので、5 から $\overline{f^{-1}(F)}\subset f^{-1}(F)$ である。逆の包含は常に成り立つから $f^{-1}(F)$ は閉である。

**1 ⇒ 6.** $B^\circ$ は開なので $f^{-1}(B^\circ)$ は開であり、$f^{-1}(B)$ に含まれる。内部は含まれる最大の開集合であるから $f^{-1}(B^\circ)\subset(f^{-1}(B))^\circ$ である。

**6 ⇒ 1.** $V\subset Y$ を開とすると $V^\circ=V$ なので $f^{-1}(V)\subset(f^{-1}(V))^\circ\subset f^{-1}(V)$ となり、$f^{-1}(V)$ は自身の内部に等しく、開である。
:::

条件 4 の包含は一般に等号にならない。例えば $f=\exp:\R\to\R$ と $A=\R$ では $f(\overline A)=(0,\infty)$ だが $\overline{f(A)}=[0,\infty)$ である。連続写像は閉包の点を閉包の点に移すが、像の閉包に新しく現れる点（ここでは $0$）を捉えるとは限らない。条件 5 も同様で、$f$ を定数写像 $f\equiv0:\R\to\R$、$B=(0,1)$ とすると $\overline{f^{-1}(B)}=\varnothing$ だが $f^{-1}(\overline B)=\R$ である。

点ごとの連続性も閉包で表せる。

:::proposition[点での連続性の閉包による判定]{#prop-top-continuity-point-closure}
$f:X\to Y$ と $x\in X$ について、$f$ が $x$ で連続であることと、「$x\in\overline A$ となる任意の $A\subset X$ について $f(x)\in\overline{f(A)}$」が成り立つことは同値である。
:::

:::proof
$f$ が $x$ で連続で $x\in\overline A$ とする。$V$ を $f(x)$ を含む開集合とすると、$f(U)\subset V$ となる開集合 $U\ni x$ がある。閉包の近傍判定（:ref[prop-top-closure-test]）により $U\cap A$ は空でなく、その点 $a$ について $f(a)\in V\cap f(A)$ である。$V$ は任意なので、再び :ref[prop-top-closure-test] により $f(x)\in\overline{f(A)}$ である。

逆に $f$ が $x$ で連続でないとする。すると $f(x)$ を含む開集合 $V$ で、$x$ を含むどの開集合 $U$ も $f(U)\not\subset V$、すなわち $U\cap(X\setminus f^{-1}(V))\ne\varnothing$ となるものがある。$A=X\setminus f^{-1}(V)$ とおくと、:ref[prop-top-closure-test] により $x\in\overline A$ である。一方 $f(A)\subset Y\setminus V$ であり、$V$ は $f(x)$ を含み $f(A)$ と交わらない開集合なので $f(x)\notin\overline{f(A)}$ である。これは条件に反する。
:::

### 開基・準開基による検査と合成

すべての開集合の逆像を調べるのは大変であるが、開集合を生成する小さな族で調べれば足りる。

:::proposition[準開基による検査と合成]{#prop-top-subbase-continuity}
1. $\mathcal S$ を $Y$ の位相の準開基（特に開基）とする（:ref[def-top-basis]）。写像 $f:X\to Y$ が連続であるための必要十分条件は、すべての $S\in\mathcal S$ について $f^{-1}(S)$ が開であることである。
2. 恒等写像と定数写像は連続である。$f:X\to Y$ と $g:Y\to Z$ が連続なら、合成 $g\circ f:X\to Z$ も連続である。
:::

:::proof
1. 必要性は明らかである。十分性を示す。:ref[prop-top-basis-test] により、$\mathcal S$ の元の有限個の共通部分（空の共通部分は $Y$ とする）の全体 $\mathcal B$ は $Y$ の位相の開基であり、$Y$ の任意の開集合は $\mathcal B$ の元の和で表される。$V=\bigcup_\lambda(S_{\lambda,1}\cap\dots\cap S_{\lambda,k_\lambda})$ と書けば

$$
f^{-1}(V)=\bigcup_\lambda\bigl(f^{-1}(S_{\lambda,1})\cap\dots\cap f^{-1}(S_{\lambda,k_\lambda})\bigr)
$$

であり、右辺は開集合の有限個の共通部分の和なので開である（空の共通部分については $f^{-1}(Y)=X$ が開である）。

2. 恒等写像と定数写像は例 :ref[exm-top-continuity-extremes] で示した。合成については、$W\subset Z$ が開なら $(g\circ f)^{-1}(W)=f^{-1}(g^{-1}(W))$ であり、$g^{-1}(W)$ が開、さらにその $f$ による逆像も開である。
:::

:::proposition[合成写像の点ごとの連続性]{#prop-top-continuity-composition}
$f:X\to Y$ が $x$ で連続、$g:Y\to Z$ が $f(x)$ で連続なら、$g\circ f$ は $x$ で連続である。
:::

:::proof
$N$ を $g(f(x))$ の近傍とすると、$g^{-1}(N)$ は $f(x)$ の近傍であり、さらに $f^{-1}(g^{-1}(N))=(g\circ f)^{-1}(N)$ は $x$ の近傍である。
:::

:::example[実数値関数の連続性を開半直線で調べる]
$\R$ の通常の位相では、開区間 $(a,b)=(-\infty,b)\cap(a,\infty)$ が開基をなすので、開半直線全体 $\{(-\infty,b)\}\cup\{(a,\infty)\}$ は準開基である。したがって $f:X\to\R$ が連続であるための必要十分条件は、すべての $a\in\R$ について

$$
\{x:f(x)>a\}=f^{-1}((a,\infty)),\qquad \{x:f(x)<a\}=f^{-1}((-\infty,a))
$$

がともに開集合であることである。この二つの条件を別々に要求したものが、後で述べる下半連続性と上半連続性である。
:::

### ネットと点列による判定

距離空間では連続性を点列の収束で判定できた。一般の位相空間では点列だけでは足りず、[位相空間](/math/topology/topological-spaces/) で導入したネットが必要になる。ネット $(x_\alpha)_{\alpha\in D}$ が $x$ に収束するとは、$x$ の任意の近傍 $N$ に対してある $\alpha_0\in D$ があり、$\alpha\ge\alpha_0$ ならば $x_\alpha\in N$ となることであった。

:::proposition[収束を保つ写像]{#prop-top-net-continuity}
1. $f:X\to Y$ が $x$ で連続であるための必要十分条件は、$x$ に収束する任意のネット $(x_\alpha)$ に対して $(f(x_\alpha))$ が $f(x)$ に収束することである。
2. $X$ が第一可算（:ref[def-top-countability]）ならば、$f$ が $x$ で連続であるための必要十分条件は、$x$ に収束する任意の点列 $(x_n)$ に対して $f(x_n)\to f(x)$ となることである。
:::

:::proof
1. $f$ が $x$ で連続で $x_\alpha\to x$ とする。$f(x)$ の近傍 $N$ に対して $f^{-1}(N)$ は $x$ の近傍なので、ある $\alpha_0$ 以降 $x_\alpha\in f^{-1}(N)$、すなわち $f(x_\alpha)\in N$ となる。よって $f(x_\alpha)\to f(x)$ である。

逆に $f$ が $x$ で連続でないとする。$f(x)$ を含む開集合 $V$ で、$x$ を含むどの開集合 $U$ についても $f(U)\not\subset V$ となるものがある。$x$ を含む開集合全体 $D$ に、$U\le U'\iff U\supset U'$ で順序を入れると、$U\cap U'$ が $U,U'$ の共通の上界となるので $D$ は有向集合である。各 $U\in D$ に対し $x_U\in U$ で $f(x_U)\notin V$ となるものを選ぶ（選択公理を用いる）。$x$ の任意の近傍 $N$ に対し、$N$ に含まれる開集合 $U_0\ni x$ をとれば、$U\ge U_0$ のとき $x_U\in U\subset U_0\subset N$ なので $x_U\to x$ である。しかし $f(x_U)$ は一度も $V$ に入らないので $f(x)$ に収束しない。

2. 必要性は 1 の特別な場合（$D=\N$）である。十分性を対偶で示す。$f$ が $x$ で連続でないとし、上の $V$ をとる。$x$ の可算な基本近傍系 $N_1,N_2,\dots$ をとり、各 $N_k$ に含まれる開集合 $W_k\ni x$ を選んで $U_n=W_1\cap\dots\cap W_n$ とおく。$U_n$ は $x$ を含む開集合の減少列で、$x$ の任意の近傍はある $N_k\supset W_k\supset U_k$ を含む。各 $n$ に $x_n\in U_n$ で $f(x_n)\notin V$ となるものを選ぶと、$x$ の任意の近傍 $N\supset U_k$ に対して $n\ge k$ なら $x_n\in U_n\subset U_k\subset N$ なので $x_n\to x$ であるが、$f(x_n)\not\to f(x)$ である。
:::

距離空間は第一可算なので、2 は :ref[thm-top-metric-continuity-tests] の点列条件を一般化している。第一可算性を外すと 2 は成り立たない。

:::example[点列の極限を保つが連続でない写像]{#exm-top-cocountable-sequential-continuity}
$X$ を非可算集合とし、補可算位相（空集合と、補集合が高々可算な集合を開とする位相）を入れる。$Y$ は同じ集合に離散位相を入れた空間とし、恒等写像 $\id:X\to Y$ を考える。

まず $X$ で $x_n\to x$ ならば、ある番号以降 $x_n=x$ であることを示す。$C=\{x_n:n\ge1\}\setminus\{x\}$ は高々可算なので、$X\setminus C$ は $x$ を含む開集合である。したがってある $N$ 以降 $x_n\in X\setminus C$ であるが、$x_n\in C\cup\{x\}$ なので $x_n=x$ となる。このような点列は離散空間 $Y$ でも $x$ に収束する。よって $\id$ は点列の極限を保つ。

一方、$\{x\}$ は $Y$ で開であるが、その逆像 $\{x\}$ は補集合 $X\setminus\{x\}$ が非可算なので $X$ で開でない。よって $\id$ は連続でない。:ref[prop-top-net-continuity] の 2 により、$X$ は第一可算でないこともわかる。
:::

### 制限と値域の取り替え

部分集合 $A\subset X$ には、$X$ の開集合と $A$ との共通部分 $A\cap U$ 全体を開集合とする位相が入る。$A\cap\varnothing=\varnothing$, $A\cap X=A$ であり、$\bigcup_i(A\cap U_i)=A\cap\bigcup_iU_i$, $(A\cap U)\cap(A\cap U')=A\cap(U\cap U')$ なので、これは位相の公理を満たす。これを**相対位相**という。詳しくは [部分空間・積空間・商空間](/math/topology/constructions/) で扱うが、本章では次の性質だけを使う：$A$ の閉集合は $X$ の閉集合 $F$ によって $A\cap F$ と書ける集合である（$A\setminus(A\cap U)=A\cap(X\setminus U)$ による）。特に $A$ 自身が $X$ で開（閉）ならば、$A$ で開（閉）な集合は $X$ でも開（閉）である。

:::proposition[制限と値域の取り替え]{#prop-top-restriction}
$f:X\to Y$ を写像とする。部分集合にはすべて相対位相を入れる。

1. $f$ が連続ならば、任意の $A\subset X$ に対し制限 $f|_A:A\to Y$ は連続である。
2. $f(X)\subset B\subset Y$ とし、$f$ を $B$ への写像 $f_B:X\to B$ とみなす。$f$ が連続であることと $f_B$ が連続であることは同値である。
3. $Y\subset Z$ を部分空間とし、$\iota:Y\to Z$ を包含写像とする。$\iota$ は連続であり、$f:X\to Y$ が連続であることと $\iota\circ f:X\to Z$ が連続であることは同値である。
:::

:::proof
1. $V\subset Y$ が開なら $(f|_A)^{-1}(V)=A\cap f^{-1}(V)$ であり、$f^{-1}(V)$ は $X$ で開なので、これは $A$ で開である。

2. $B$ の開集合は $Y$ の開集合 $V$ によって $B\cap V$ と書け、$f(X)\subset B$ より $f_B^{-1}(B\cap V)=f^{-1}(V)$ である。したがって「$Y$ の開集合の逆像がすべて開」と「$B$ の開集合の逆像がすべて開」は同じ条件である。

3. $Z$ の開集合 $W$ について $\iota^{-1}(W)=Y\cap W$ は $Y$ で開なので $\iota$ は連続である。$\iota\circ f$ は $f$ を値域 $Z$ の写像とみなしたものであるから、2 を $Y\subset Z$ に適用すればよい。
:::

値域を部分集合に制限するときには必ず相対位相を入れる。値域に別の位相を入れると連続性は変わりうる。例えば恒等写像 $\R\to\R$ で値域だけに離散位相を入れると、$\{0\}$ の逆像 $\{0\}$ が開でないので連続でなくなる。

### 実数値連続関数の演算

実数値関数 $f:X\to\R$ が $x$ で連続であることは、開区間 $(f(x)-\eps,f(x)+\eps)$ $(\eps>0)$ が $f(x)$ の基本近傍系をなすことから、次のように言い換えられる：任意の $\eps>0$ に対し $x$ の近傍 $U$ があって、$y\in U$ ならば $\abs{f(y)-f(x)}<\eps$ となる。この形を使えば、微分積分学で学んだ連続関数の演算の議論を、定義域の位相だけを用いてそのまま再現できる。

:::proposition[実数値連続関数の演算]{#prop-top-real-function-algebra}
$f,g:X\to\R$ が $x\in X$ で連続であり、$c\in\R$ とする。

1. $f+g$, $f-g$, $cf$, $fg$ は $x$ で連続である。
2. $\abs f$, $\max(f,g)$, $\min(f,g)$ は $x$ で連続である。
3. $g(x)\ne0$ ならば、$G=\{y:g(y)\ne0\}$ 上で定義された $f/g$ は（$G$ に相対位相を入れて）$x$ で連続である。

特に $f,g\in C(X)$ ならば上の関数はすべて連続であり、$C(X)$ は関数の和と積について可換環をなす。
:::

:::proof
$\eps>0$ を固定する。$\delta>0$ を後で定め、$x$ の近傍 $U_f,U_g$ を $y\in U_f\Rightarrow\abs{f(y)-f(x)}<\delta$、$y\in U_g\Rightarrow\abs{g(y)-g(x)}<\delta$ となるようにとる。$U=U_f\cap U_g$ は $x$ の近傍であり、以下 $y\in U$ とする。

1. 和と差：$\abs{(f\pm g)(y)-(f\pm g)(x)}\le\abs{f(y)-f(x)}+\abs{g(y)-g(x)}<2\delta$ なので $\delta=\eps/2$ とすればよい。スカラー倍は $g$ を定数関数 $c$ とした積の特別な場合である。積については、$\delta\le1$ なら $\abs{g(y)}\le\abs{g(x)}+1$ であり、

$$
\begin{aligned}
\abs{f(y)g(y)-f(x)g(x)}
&\le\abs{f(y)-f(x)}\,\abs{g(y)}+\abs{f(x)}\,\abs{g(y)-g(x)}\\
&\le\delta\bigl(\abs{g(x)}+1\bigr)+\abs{f(x)}\,\delta
=\delta\bigl(\abs{f(x)}+\abs{g(x)}+1\bigr)
\end{aligned}
$$

となる。$\delta=\min\bigl(1,\ \eps/(2(\abs{f(x)}+\abs{g(x)}+1))\bigr)$ とすれば右辺は $\eps/2<\eps$ 以下である。

2. 三角不等式から $\bigl|\abs{f(y)}-\abs{f(x)}\bigr|\le\abs{f(y)-f(x)}<\delta$ なので、$\delta=\eps$ とすれば $\abs f$ は $x$ で連続である。さらに実数 $s,t$ について

$$
\max(s,t)=\frac{s+t+\abs{s-t}}2,\qquad \min(s,t)=\frac{s+t-\abs{s-t}}2
$$

が成り立つ（$s\ge t$ と $s<t$ に分けて確かめられる）ので、$\max(f,g)$, $\min(f,g)$ は 1 と $\abs{\cdot}$ の連続性を組み合わせて得られる。ここで $\abs{f-g}$ が $x$ で連続であることは、$f-g$ が $x$ で連続であることと、$\abs{\cdot}$ についての上の議論を $f-g$ に適用することから従う。

3. まず $x$ が $G$ の内点であることを示す。$g$ は $x$ で連続なので、$\delta\le\abs{g(x)}/2$ ととれば $y\in U_g$ について $\abs{g(y)}\ge\abs{g(x)}-\abs{g(y)-g(x)}>\abs{g(x)}/2>0$ となり、$U_g\subset G$ である。したがって $x$ は $G$ の内点であり、$U_g$ の中で $1/g$ は定義される。$y\in U_g$ について

$$
\Bigl|\frac1{g(y)}-\frac1{g(x)}\Bigr|=\frac{\abs{g(x)-g(y)}}{\abs{g(x)}\,\abs{g(y)}}
<\frac{2\delta}{\abs{g(x)}^2}
$$

であるから、$\delta=\min\bigl(\abs{g(x)}/2,\ \eps\abs{g(x)}^2/2\bigr)$ とすれば $\abs{1/g(y)-1/g(x)}<\eps$ となる。$G\cap U_g$ は $G$ における $x$ の近傍なので、$1/g:G\to\R$ は $x$ で連続である。$f/g=f\cdot(1/g)$ は 1 の積の主張（$X$ を $G$ に置き換え、$f|_G$ は :ref[prop-top-restriction] の議論と同様に $x$ で連続）から $x$ で連続である。

$f,g$ がすべての点で連続なら、上の関数はすべての点で連続であり、:ref[thm-top-continuity-tests] により連続写像である。$g$ が連続なら $G=g^{-1}(\R\setminus\{0\})$ は開集合である。環の公理は各点の値の等式から従う。
:::

:::corollary[連続関数で定まる集合]{#cor-top-real-function-sets}
$f,g\in C(X)$ ならば、$\{x:f(x)<g(x)\}$ は開集合であり、$\{x:f(x)\le g(x)\}$ と $\{x:f(x)=g(x)\}$ は閉集合である。特に零点集合 $f^{-1}(0)$ は閉集合である。
:::

:::proof
$h=g-f$ は連続である。$\{f<g\}=h^{-1}((0,\infty))$ は開集合の逆像なので開であり、$\{f\le g\}=h^{-1}([0,\infty))$, $\{f=g\}=h^{-1}(\{0\})$ は閉集合の逆像なので閉である（:ref[thm-top-continuity-tests]）。
:::

Euclid 空間 $\R^m$ には Euclid 距離 $d_2(x,y)=\norm{x-y}_2$ を入れる。Euclid 空間に値をとる写像の連続性は成分ごとに判定できる。

:::proposition[Euclid 空間への写像の成分による判定]{#prop-top-euclidean-components}
$f=(f_1,\dots,f_m):X\to\R^m$ が連続であるための必要十分条件は、各成分 $f_i:X\to\R$ が連続であることである。特に座標関数 $\pi_i:\R^m\to\R$, $\pi_i(y)=y_i$ と Euclid ノルム $y\mapsto\norm y_2$ は連続である。
:::

:::proof
$\abs{\pi_i(y)-\pi_i(z)}=\abs{y_i-z_i}\le\norm{y-z}_2$ であり、逆三角不等式（:ref[lem-top-reverse-triangle]）から $\bigl|\norm y_2-\norm z_2\bigr|\le\norm{y-z}_2$ である。したがって $\pi_i$ とノルムは Lipschitz 連続であり、連続である。よって $f$ が連続なら $f_i=\pi_i\circ f$ は連続である（:ref[prop-top-subbase-continuity]）。

逆に各 $f_i$ が連続とし、$x\in X$, $\eps>0$ をとる。各 $i$ に対し $x$ の近傍 $U_i$ で $y\in U_i\Rightarrow\abs{f_i(y)-f_i(x)}<\eps/\sqrt m$ となるものがある。$U=U_1\cap\dots\cap U_m$ は $x$ の近傍であり、$y\in U$ ならば

$$
\norm{f(y)-f(x)}_2=\Bigl(\sum_{i=1}^m\abs{f_i(y)-f_i(x)}^2\Bigr)^{1/2}<\Bigl(m\cdot\frac{\eps^2}m\Bigr)^{1/2}=\eps
$$

である。よって $f^{-1}(B_{d_2}(f(x),\eps))$ は $x$ の近傍であり、開球が $f(x)$ の基本近傍系をなすことから $f$ は $x$ で連続である。
:::

:::example[多項式で定まる集合]{#exm-top-polynomial-sets}
:ref[prop-top-euclidean-components] と :ref[prop-top-real-function-algebra] により、$\R^m$ 上の実係数多項式関数はすべて連続である（座標関数から和・積・スカラー倍を有限回繰り返して得られる）。したがって :ref[cor-top-real-function-sets] から

- 単位球面 $S^{m-1}=\{y:y_1^2+\dots+y_m^2=1\}$ と閉単位球 $D^m=\{y:\norm y_2\le1\}$ は閉集合、開単位球 $B^m=\{y:\norm y_2<1\}$ は開集合である。
- 双曲線 $\{(x,y)\in\R^2:xy=1\}$ は閉集合、開上半平面 $\{(x,y):y>0\}$ は開集合である。
- 実係数多項式 $p$ について $\{p\ne0\}$ は開集合である。例えば $\{(x,y):x^2+y^2\ne1\}$ は開である。

これらの事実は本章の同相写像の構成や次章の例で繰り返し用いる。
:::

:::remark[連続関数の演算と積位相]
$f+g$ は、写像 $x\mapsto(f(x),g(x))\in\R^2$ と加法 $\R^2\to\R$ の合成である。前者の連続性は :ref[prop-top-euclidean-components]、後者の連続性は $\abs{(s+t)-(s'+t')}\le\sqrt2\,\norm{(s,t)-(s',t')}_2$ からわかるので、この合成として示すこともできる。積や商も同様である。この見方は、次章で積位相を導入すると一般の位相群などへ拡張される。
:::

### 貼り合わせの補題

関数を領域ごとに異なる式で定義することはよくある。その関数が連続であるかどうかは、各部分での連続性と、部分の取り方に依存する。

:::proposition[開被覆と有限閉被覆による貼り合わせ]{#prop-top-continuity-pasting}
$f:X\to Y$ を写像とし、各部分集合には相対位相を入れる。

1. $X=\bigcup_{i\in I}U_i$ が開集合による被覆であり、各制限 $f|_{U_i}$ が連続ならば $f$ は連続である。特に、各点が開近傍 $U$ で $f|_U$ が連続となるものをもてば $f$ は連続である（連続性は局所的な性質である）。
2. $X=F_1\cup\dots\cup F_k$ が有限個の閉集合による被覆であり、各制限 $f|_{F_j}$ が連続ならば $f$ は連続である。

また、被覆 $(A_i)$ の各部分に連続写像 $f_i:A_i\to Y$ が与えられ、$A_i\cap A_j$ 上で $f_i=f_j$ が成り立つならば、$f|_{A_i}=f_i$ となる写像 $f:X\to Y$ がただ一つ存在する。この $f$ は、$(A_i)$ が開被覆または有限閉被覆なら連続である。
:::

:::proof
1. $V\subset Y$ を開集合とする。$(f|_{U_i})^{-1}(V)=U_i\cap f^{-1}(V)$ は $U_i$ で開なので、$X$ の開集合 $W_i$ によって $U_i\cap W_i$ と書ける。これは開集合の共通部分なので $X$ で開である。したがって

$$
f^{-1}(V)=\bigcup_i\bigl(U_i\cap f^{-1}(V)\bigr)
$$

は開集合の和であり開である。後半は、各点に選んだ開近傍の全体が開被覆をなすことから従う。

2. $H\subset Y$ を閉集合とする。$(f|_{F_j})^{-1}(H)=F_j\cap f^{-1}(H)$ は $F_j$ で閉なので、$X$ の閉集合 $C_j$ によって $F_j\cap C_j$ と書け、$X$ で閉である。$f^{-1}(H)=\bigcup_{j=1}^k(F_j\cap f^{-1}(H))$ は有限個の閉集合の和なので閉である。:ref[thm-top-continuity-tests] により $f$ は連続である。

最後の主張について、$x\in X$ に対し $x\in A_i$ となる $i$ を選んで $f(x)=f_i(x)$ と定める。別の $j$ で $x\in A_j$ なら $x\in A_i\cap A_j$ なので $f_i(x)=f_j(x)$ であり、値は $i$ の選び方によらない。一意性は $f|_{A_i}=f_i$ の要請から明らかである。連続性は 1, 2 による。
:::

:::example[絶対値と区分的な定義]
$\R=(-\infty,0]\cup[0,\infty)$ は二つの閉集合による被覆であり、$f_1(x)=-x$, $f_2(x)=x$ は重なり $\{0\}$ で一致する。:ref[prop-top-continuity-pasting] の 2 により、貼り合わせた関数 $\abs x$ は連続である。同様に、閉区間 $[0,1]$ 上で $t\le1/2$ では $2t$、$t\ge1/2$ では $2-2t$ と定めたテント関数も連続である。
:::

:::example[貼り合わせの補題が使えない場合]
1. **閉でない部分**：$A=(-\infty,0)$, $B=[0,\infty)$ とし、$f$ を $A$ 上で $0$、$B$ 上で $1$ とする。$f|_A$, $f|_B$ は定数なので連続であるが、$f$ は $0$ で連続でない（$0$ の任意の近傍に値 $0$ の点がある）。$A\cup B=\R$ は被覆であり重なりもないが、開集合 $A$ は閉でなく、閉集合 $B$ は開でないため、1 も 2 も使えない。
2. **無限個の閉集合**：$\R$ を一点集合 $\{x\}$ $(x\in\R)$ で覆うと、各部分は閉で、任意の関数の各制限は連続である。それでも不連続な関数はいくらでもある。有限性の仮定は外せない。
:::

局所有限な閉被覆への拡張も同じ逆像の等式による。閉集合の逆像は各部分で閉な集合の局所有限和だから閉である。

### 一様極限

連続関数の各点収束極限は連続とは限らないが、一様収束極限は連続である。この事実は定義域に距離を必要としない。

:::proposition[一様極限の連続性]{#prop-top-uniform-limit}
$X$ を位相空間、$(Y,\rho)$ を距離空間とする。連続写像 $f_n:X\to Y$ が $f:X\to Y$ に一様収束する、すなわち任意の $\eps>0$ に対しある $N$ があって、$n\ge N$ と全 $x\in X$ について $\rho(f_n(x),f(x))<\eps$ となるならば、$f$ は連続である。
:::

:::proof
$x\in X$ と $\eps>0$ を固定する。一様収束より、ある $n$ で全 $y\in X$ について $\rho(f_n(y),f(y))<\eps/3$ となる。$f_n$ は $x$ で連続なので、$x$ の近傍 $U$ で $y\in U\Rightarrow\rho(f_n(y),f_n(x))<\eps/3$ となるものがある。$y\in U$ ならば

$$
\rho(f(y),f(x))\le\rho(f(y),f_n(y))+\rho(f_n(y),f_n(x))+\rho(f_n(x),f(x))<\eps
$$

である。よって $f^{-1}(B_\rho(f(x),\eps))\supset U$ は $x$ の近傍であり、$f$ は $x$ で連続である（:ref[prop-top-continuity-metric-agree] の証明と同様に、開球が基本近傍系をなすことを用いた）。証明では一つの $n$ を固定してから近傍を選んでいるので、$X$ に距離は不要である。
:::

:::example[各点収束では連続性が失われる]
$[0,1]$ 上で $f_n(x)=x^n$ は各点で $f(x)=0$ $(0\le x<1)$, $f(1)=1$ に収束する。$f^{-1}((1/2,2))=\{1\}$ は $[0,1]$ で開でないので $f$ は連続でない。実際この収束は一様でない：$x_n=2^{-1/n}\in[0,1)$ とおくと $f(x_n)=0$ だが $f_n(x_n)=1/2$ であり、どの $n$ でも $\sup_x\abs{f_n(x)-f(x)}\ge1/2$ である。
:::

### 半連続性

連続性の準開基による判定で、二種類の開半直線の片方だけを要求すると、最大値・最小値の問題や関数の上限・下限の操作に適した弱い連続性が得られる。

:::definition[半連続関数]{#def-top-semicontinuous}
$f:X\to\R$ とする。すべての $a\in\R$ について $\{x:f(x)>a\}$ が開集合であるとき $f$ は :term[下半連続]{en="lower semicontinuous" yomi="かはんれんぞく"}であるといい、すべての $a$ について $\{x:f(x)<a\}$ が開集合であるとき :term[上半連続]{en="upper semicontinuous" yomi="じょうはんれんぞく"}であるという。
:::

$f$ が上半連続であることと $-f$ が下半連続であることは同値である（$\{-f>a\}=\{f<-a\}$）。開半直線が $\R$ の準開基をなすので、:ref[prop-top-subbase-continuity] により $f$ が連続であることは、下半連続かつ上半連続であることと同値である。

点ごとには次のように言い換えられる：$f$ が下半連続であることは、各 $x$ と各 $a<f(x)$ に対して $x$ の近傍 $U$ で $U$ 上 $f>a$ となるものが存在することと同値である。実際、前者から後者は $U=\{f>a\}$ をとればよく、後者を仮定すると $\{f>a\}$ は自身の各点の近傍を含むので開である。直観的には、下半連続関数は近くの点で値が「下に跳ぶ」ことがない。

:::example[特性関数と床関数]
1. $U\subset X$ の特性関数 $\chi_U:X\to\R$ について、$\{\chi_U>a\}$ は $a<0$ なら $X$、$0\le a<1$ なら $U$、$a\ge1$ なら $\varnothing$ である。したがって $\chi_U$ が下半連続であることは $U$ が開であることと同値である。同様に $\{\chi_U<a\}$ は $a\le0$ なら $\varnothing$、$0<a\le1$ なら $X\setminus U$、$a>1$ なら $X$ なので、上半連続であることは $U$ が閉であることと同値である。例えば $\chi_{(0,\infty)}$ は下半連続だが上半連続でない。
2. 床関数 $\lfloor x\rfloor$（$x$ 以下の最大の整数）について、$\lfloor x\rfloor<a\iff\lfloor x\rfloor\le\lceil a\rceil-1\iff x<\lceil a\rceil$ なので $\{\lfloor x\rfloor<a\}=(-\infty,\lceil a\rceil)$ は開であり、床関数は上半連続である。一方 $\lfloor x\rfloor>a\iff\lfloor x\rfloor\ge\lfloor a\rfloor+1\iff x\ge\lfloor a\rfloor+1$ なので $\{\lfloor x\rfloor>a\}=[\lfloor a\rfloor+1,\infty)$ は開でなく、下半連続でない。
:::

:::proposition[点ごとの上限と半連続性]{#prop-top-semicontinuity}
$(f_i)_{i\in I}$ を下半連続関数 $X\to\R$ の空でない族とし、各点で $f(x)=\sup_if_i(x)$ が有限であるとする。このとき $f$ は下半連続である。同様に、上半連続関数の空でない族の点ごとの下限が有限値ならば、それは上半連続である。
:::

:::proof
$f(x)>a$ であることは、上限の定義により、ある $i$ で $f_i(x)>a$ となることと同値である。したがって

$$
\{x:f(x)>a\}=\bigcup_{i\in I}\{x:f_i(x)>a\}
$$

であり、右辺は開集合の和なので開である。下限については $\{\inf_if_i<a\}=\bigcup_i\{f_i<a\}$ を用いる。
:::

連続関数は下半連続なので、連続関数族の上限は下半連続である。しかし連続とは限らない。例えば $f_n(x)=\min(1,\max(0,nx))$ は $\R$ 上の連続関数であり、その上限は $\chi_{(0,\infty)}$ である。$\{\chi_{(0,\infty)}<1/2\}=(-\infty,0]$ は開でないので、上限は上半連続でなく、連続でもない。

## 開写像と閉写像

### 定義と基本性質

連続性は逆像についての条件であった。逆に、像についての条件を課すと別の種類の写像が得られる。これらは連続性とは独立な概念であるが、連続写像の逆写像の連続性を調べるとき（同相写像）や、商空間の位相を調べるとき（次章）に本質的な役割を果たす。

:::definition[開写像と閉写像]{#def-top-open-closed-map}
写像 $f:X\to Y$ が、$X$ の任意の開集合 $U$ に対して $f(U)$ を $Y$ の開集合に写すとき :term[開写像]{en="open map" yomi="かいしゃぞう"}、$X$ の任意の閉集合 $F$ に対して $f(F)$ を $Y$ の閉集合に写すとき :term[閉写像]{en="closed map" yomi="へいしゃぞう"}という。連続性は仮定しない。
:::

像は和と交換する（$f(\bigcup_iA_i)=\bigcup_if(A_i)$）が、共通部分とは交換しない。例えば $f(t)=t^2$, $A=\{-1\}$, $B=\{1\}$ とすると $f(A\cap B)=\varnothing$ だが $f(A)\cap f(B)=\{1\}$ である。補集合とも交換しない：$f(\R\setminus\{-1\})=[0,\infty)$ は $\R\setminus f(\{-1\})=\R\setminus\{1\}$ と異なる。このため、開写像と閉写像は（逆像の場合と異なり）補集合を取って移り合う関係にはなく、互いに独立な概念である。ただし全単射なら像と補集合は交換するので、次が成り立つ。

:::proposition[開写像の判定と合成]{#prop-top-open-map-base-test}
1. $\mathcal B$ を $X$ の開基とする。$f:X\to Y$ が開写像であるための必要十分条件は、すべての $B\in\mathcal B$ について $f(B)$ が開であることである。
2. $f:X\to Y$ が開写像であるための必要十分条件は、各 $x\in X$ と $x$ の各近傍 $N$ について $f(N)$ が $f(x)$ の近傍であることである。
3. 全単射 $f$ について、$f$ が開写像であることと閉写像であることは同値である。
:::

:::proof
1. 必要性は開基の元が開であることから明らかである。逆に $U$ が開なら $U=\bigcup_\lambda B_\lambda$ $(B_\lambda\in\mathcal B)$ と書け、$f(U)=\bigcup_\lambda f(B_\lambda)$ は開集合の和なので開である。

2. $f$ が開写像で $N$ が $x$ の近傍なら、$x\in U\subset N$ となる開集合 $U$ について $f(U)$ は $f(x)$ を含む開集合で $f(N)$ に含まれる。よって $f(N)$ は $f(x)$ の近傍である。逆に条件を仮定し、$U$ を開集合とする。各 $y\in f(U)$ について $y=f(x)$, $x\in U$ とすると、$U$ は $x$ の近傍なので $f(U)$ は $y$ の近傍である。自身の各点の近傍であるような集合は、各点を含む開集合の和として開である。

3. 全単射では $f(X\setminus A)=Y\setminus f(A)$ が成り立つ。したがって開集合 $U$ の像が開であることと、閉集合 $X\setminus U$ の像 $Y\setminus f(U)$ が閉であることは同値である。
:::

:::proposition[開写像・閉写像の合成と制限]{#prop-top-open-closed-map-composition}
1. $f:X\to Y$, $g:Y\to Z$ がともに開写像（閉写像）ならば、$g\circ f$ も開写像（閉写像）である。
2. $f:X\to Y$ が開写像で $U\subset X$ が開集合ならば、$f|_U:U\to Y$ は開写像である。$f$ が閉写像で $F\subset X$ が閉集合ならば、$f|_F:F\to Y$ は閉写像である。
3. $f:X\to Y$ が開写像（閉写像）で $f(X)\subset B\subset Y$ ならば、$f$ を $B$ への写像とみなしても開写像（閉写像）である。
4. $f:X\to Y$ が連続な全射で $g\circ f$ が開写像（閉写像）ならば、$g:Y\to Z$ は開写像（閉写像）である。
:::

:::proof
1. $(g\circ f)(A)=g(f(A))$ による。

2. $U$ の開集合は $U\cap W$（$W$ は $X$ の開集合）の形であり、$U$ が開なので $X$ の開集合である。よってその像は開である。閉の場合も、$F$ の閉集合は $X$ の閉集合なので同様である。

3. $V=f(U)$ が $Y$ で開なら、$V=B\cap V$ は $B$ で開である。閉の場合も同様である。

4. $V\subset Y$ を開集合とすると、全射性から $V=f(f^{-1}(V))$ であり、$f^{-1}(V)$ は開なので、$g(V)=(g\circ f)(f^{-1}(V))$ は開である。閉の場合は閉集合 $H\subset Y$ に $H=f(f^{-1}(H))$ を用いる。
:::

### 例と反例

連続性・開写像性・閉写像性の 3 つの性質は互いに独立であり、どの組み合わせも起こりうる。以下で典型例をまとめる。

:::example[射影は開だが閉でない]{#exm-top-projection-open}
$p:\R^2\to\R$, $p(x,y)=x$ とする。$\abs{p(z)-p(w)}\le\norm{z-w}_2$ なので $p$ は連続である。

$p$ が開写像であることを示す。$U\subset\R^2$ を開集合、$(x,y)\in U$ とし、$B((x,y),r)\subset U$ となる $r>0$ をとる。$\abs{t-x}<r$ なら $(t,y)\in B((x,y),r)\subset U$ なので $(x-r,x+r)\subset p(U)$ である。よって $p(U)$ は開である。

一方、双曲線 $H=\{(x,y):xy=1\}$ は閉集合である（例 :ref[exm-top-polynomial-sets]）が、$p(H)=\R\setminus\{0\}$ は閉でない。よって $p$ は閉写像でない。
:::

:::example[平方写像は閉だが開でない]{#exm-top-square-closed}
$f:\R\to\R$, $f(t)=t^2$ は連続であるが、$f((-1,1))=[0,1)$ は開でないので開写像でない。

$f$ は閉写像である。$F\subset\R$ を閉集合とする。:ref[prop-top-metric-closed-sequence] により、$f(F)$ の点列 $f(x_n)=x_n^2$ $(x_n\in F)$ が $y\in\R$ に収束するとき $y\in f(F)$ となることを示せばよい。$x_n^2\ge0$ なので $y\ge0$ であり、平方根の連続性により $\abs{x_n}=\sqrt{x_n^2}\to\sqrt y$ である。$x_n\ge0$ となる $n$ が無限個あればその部分列は $\sqrt y$ に、そうでなければ $x_n<0$ となる部分列が $-\sqrt y$ に収束する。$F$ は閉なので、その極限 $\pm\sqrt y$ は $F$ に属し、$y=f(\pm\sqrt y)\in f(F)$ である。

なお値域を $[0,\infty)$ に制限した $f:\R\to[0,\infty)$ は開写像になる。開区間 $(a,b)$ の像は、$0\in(a,b)$ なら $[0,\max(a^2,b^2))$、そうでなければ $(\min(a^2,b^2),\max(a^2,b^2))$ であり（中間値の定理による。[関数の極限と連続性](/math/analysis/continuity/) を参照）、いずれも $[0,\infty)$ で開だからである。:ref[prop-top-open-map-base-test] の 1 を用いた。値域の取り方で開写像かどうかが変わる。
:::

:::example[指数関数・定数写像・包含写像]
1. $\exp:\R\to\R$ は狭義単調増加で連続であり、中間値の定理により $\exp((a,b))=(e^a,e^b)$ なので、開区間の開基に :ref[prop-top-open-map-base-test] を適用して開写像である。しかし閉集合 $\R$ の像 $(0,\infty)$ は閉でないので閉写像でない。
2. 定数写像 $c:\R\to\R$, $c(x)=3$ は連続で、閉集合の像は $\varnothing$ または $\{3\}$ でいずれも閉なので閉写像であるが、開集合 $\R$ の像 $\{3\}$ が開でないので開写像でない。
3. 包含写像 $[0,1]\to\R$ は連続な閉写像（$[0,1]$ の閉集合は $\R$ の閉集合 $F$ によって $[0,1]\cap F$ と書け、閉集合の共通部分として閉）であるが開写像でない。包含写像 $(0,1)\to\R$ は連続な開写像であるが閉写像でない（$(0,1)$ 自身の像が閉でない）。
:::

:::example[連続でない開かつ閉な全単射と、その逆]
$X=\{0,1\}$ とし、$\mathcal T_i$ を密着位相、$\mathcal T_d$ を離散位相とする。恒等写像 $\id:(X,\mathcal T_i)\to(X,\mathcal T_d)$ は、定義域の開集合・閉集合が $\varnothing,X$ だけなので開かつ閉であるが、値域の開集合 $\{0\}$ の逆像が開でないので連続でない。逆向きの恒等写像 $(X,\mathcal T_d)\to(X,\mathcal T_i)$ は連続な全単射であるが、$\{0\}$ の像が開でも閉でもないので開写像でも閉写像でもない。
:::

### 像による判定と全射への応用

:::proposition[開写像と閉写像の集合による判定]{#prop-top-image-criteria}
任意の写像 $f:X\to Y$ について次が成り立つ。

$$
\begin{aligned}
f\text{ が開写像}&\iff f(A^\circ)\subset\bigl(f(A)\bigr)^\circ\quad(\forall A\subset X),\\
f\text{ が閉写像}&\iff \overline{f(A)}\subset f(\overline A)\quad(\forall A\subset X).
\end{aligned}
$$
:::

:::proof
$f$ が開写像なら、$f(A^\circ)$ は $f(A)$ に含まれる開集合なので $f(A)$ の内部に含まれる。逆に条件を仮定し $A$ を開集合とすると $f(A)=f(A^\circ)\subset(f(A))^\circ\subset f(A)$ なので $f(A)$ は開である。

$f$ が閉写像なら、$f(\overline A)$ は $f(A)$ を含む閉集合なので $\overline{f(A)}$ を含む。逆に条件を仮定し $A$ を閉集合とすると $\overline{f(A)}\subset f(\overline A)=f(A)$ なので $f(A)$ は閉である。
:::

閉写像の包含は、連続性の判定 :ref[thm-top-continuity-tests] の 4 と逆向きである。したがって連続な閉写像では $f(\overline A)=\overline{f(A)}$ がすべての $A$ で成り立ち、逆にこれが成り立てば $f$ は連続な閉写像である。例 :ref[exm-top-projection-open] の射影と双曲線 $H$ では $p(\overline H)=\R\setminus\{0\}\subsetneq\R=\overline{p(H)}$ であり、閉写像の条件が原点一つで破れている。

:::proposition[開全射・閉全射による判定]{#prop-top-surjective-test}
$q:X\to Y$ を連続な全射とする。$q$ が開写像または閉写像ならば、$V\subset Y$ について

$$
V\text{ が }Y\text{ で開}\iff q^{-1}(V)\text{ が }X\text{ で開}
$$

が成り立つ。同様に $H\subset Y$ が閉であることと $q^{-1}(H)$ が閉であることは同値である。
:::

:::proof
$\Rightarrow$ は連続性である。逆を示す。全射性より任意の $V\subset Y$ について $q(q^{-1}(V))=V$ である。$q$ が開写像で $q^{-1}(V)$ が開なら、$V=q(q^{-1}(V))$ は開である。$q$ が閉写像で $q^{-1}(V)$ が開なら、$X\setminus q^{-1}(V)=q^{-1}(Y\setminus V)$ は閉であり、その像 $q(q^{-1}(Y\setminus V))=Y\setminus V$ は閉なので $V$ は開である。閉集合についての主張は補集合をとれば得られる。
:::

この性質をもつ全射を **商写像** といい、[部分空間・積空間・商空間](/math/topology/constructions/) で詳しく調べる。開写像でも閉写像でもない連続全射では一般に成り立たない。例えば恒等写像 $(\R,\text{離散位相})\to(\R,\text{通常の位相})$ は連続全射であり、任意の $V$ について $q^{-1}(V)=V$ は開であるが、$V=\{0\}$ は値域で開でない。

### 閉写像の近傍による判定

閉写像であることは、ファイバー（一点の逆像）を囲む開集合の振る舞いで特徴付けられる。この形は、次章以降で商空間や積空間の性質を調べるときに使われる。

:::proposition[閉写像の近傍による判定]{#prop-top-closed-fiber-test}
写像 $q:X\to Y$ について次は同値である。

1. $q$ は閉写像である。
2. 任意の $B\subset Y$ と $q^{-1}(B)$ を含む任意の開集合 $U\subset X$ に対し、$B$ を含む開集合 $V\subset Y$ で $q^{-1}(V)\subset U$ となるものが存在する。
3. 任意の $y\in Y$ と $q^{-1}(y)$ を含む任意の開集合 $U\subset X$ に対し、$y$ を含む開集合 $V$ で $q^{-1}(V)\subset U$ となるものが存在する。
:::

:::proof
**1 ⇒ 2.** $V=Y\setminus q(X\setminus U)$ とおく。$X\setminus U$ は閉なのでその像は閉であり、$V$ は開である。$y\in B$ なら $q^{-1}(y)\subset U$ なので $y\notin q(X\setminus U)$、すなわち $y\in V$ である。また $q(x)\in V$ なら $x\notin X\setminus U$ なので $q^{-1}(V)\subset U$ である。

**2 ⇒ 3.** $B=\{y\}$ とすればよい。

**3 ⇒ 1.** $F\subset X$ を閉集合とし、$y\notin q(F)$ とする。$q^{-1}(y)\cap F=\varnothing$ なので $q^{-1}(y)$ は開集合 $U=X\setminus F$ に含まれる。3 による $V\ni y$ は $q^{-1}(V)\cap F=\varnothing$ を満たすので $V\cap q(F)=\varnothing$ である。よって $Y\setminus q(F)$ は自身の各点の開近傍を含み、開である。
:::

$q^{-1}(y)=\varnothing$（$y$ が像に入らない）の場合も 3 は意味をもち、$U=\varnothing$ として「$y$ の開近傍で像と交わらないものがある」すなわち $Y\setminus q(X)$ が開であることを要求する。

:::example[射影で縦線の幅をそろえられない]
例 :ref[exm-top-projection-open] の射影 $p(x,y)=x$ が閉写像でないことを、:ref[prop-top-closed-fiber-test] の 3 で見直す。開集合 $U=\{(x,y):\abs{xy}<1\}$（例 :ref[exm-top-polynomial-sets] により開）はファイバー $p^{-1}(0)=\{0\}\times\R$ を含む。しかし $0$ のどの開近傍 $V$ にも $t\ne0$ が含まれ、点 $(t,1/t)\in p^{-1}(V)$ は $U$ に属さない。したがって縦線全体を囲む縦の帯 $p^{-1}(V)$ を $U$ の中にとることはできない。ファイバー上の各点 $(0,y)$ には $U$ に含まれる開正方形がとれるが、その横幅は $\abs y$ が大きくなるにつれて小さくしなければならず、一様にとれないのである。ファイバーが有界閉集合のときにはこの現象が起きないことを、[コンパクト性](/math/topology/compactness/) の章で示す。
:::

## 同相写像

### 定義と同値条件

二つの位相空間の開集合の構造が、点の対応によって完全に移り合うならば、位相に関する限り両者を区別する理由はない。このような対応が同相写像である。

:::definition[同相写像と埋め込み]{#def-top-homeomorphism}
全単射 $h:X\to Y$ で、$h$ も逆写像 $h^{-1}$ も連続であるものを :term[同相写像]{en="homeomorphism" yomi="どうそうしゃぞう"}という。同相写像が存在するとき $X$ と $Y$ は :term[同相]{en="homeomorphic" yomi="どうそう"}であるといい、$X\cong Y$ と書く。

単射 $e:X\to Y$ が、$X$ から像 $e(X)$（相対位相を入れる）への同相写像を与えるとき、$e$ を :term[埋め込み]{en="embedding" yomi="うめこみ"}という。
:::

$h$ が同相写像であることは、$h$ が全単射で、$U\subset X$ が開であることと $h(U)\subset Y$ が開であることが同値であること、すなわち $h$ が開集合全体の間の全単射 $\mathcal T_X\to\mathcal T_Y$, $U\mapsto h(U)$ を引き起こすことと言い換えられる。実際、$h^{-1}$ が連続であることは $(h^{-1})^{-1}(U)=h(U)$ が開であることにほかならない。

:::proposition[同相は同値関係である]{#prop-top-homeomorphism-equivalence-relation}
恒等写像は同相写像であり、同相写像の逆写像と合成は同相写像である。したがって同相 $\cong$ は位相空間の間の同値関係である。
:::

:::proof
$\id_X$ は自身が逆写像であり連続である。$h$ が同相写像なら $h^{-1}$ は連続な全単射でその逆 $h$ も連続である。$h:X\to Y$, $k:Y\to Z$ が同相写像なら $k\circ h$ は全単射で、その逆 $h^{-1}\circ k^{-1}$ とともに連続写像の合成として連続である（:ref[prop-top-subbase-continuity]）。
:::

連続な全単射の逆写像は連続とは限らない（後の例 :ref[exm-top-angle-parametrization]）。逆写像の連続性は次のように開写像・閉写像の言葉で判定できる。

:::proposition[連続全単射の逆写像]{#prop-top-homeomorphism-test}
連続な全単射 $h:X\to Y$ について、次の 3 条件は同値である。

1. $h$ は同相写像である。
2. $h$ は開写像である。
3. $h$ は閉写像である。
:::

:::proof
$g=h^{-1}$ とすると、$U\subset X$ について $g^{-1}(U)=h(U)$ である。したがって $g$ が連続であること（開集合の逆像が開）は $h$ が開写像であることと同値であり、:ref[thm-top-continuity-tests] の 3 により $g$ が連続であること（閉集合の逆像が閉）は $h$ が閉写像であることとも同値である。
:::

:::proposition[同相写像の制限]{#prop-top-homeomorphism-restriction}
$h:X\to Y$ が同相写像で $A\subset X$ ならば、$h|_A:A\to h(A)$ は（相対位相に関して）同相写像である。
:::

:::proof
$h|_A:A\to Y$ は連続（:ref[prop-top-restriction] の 1）であり、像が $h(A)$ に含まれるので $h|_A:A\to h(A)$ も連続である（同 2）。これは全単射で、逆写像は $h^{-1}|_{h(A)}:h(A)\to A$ であり、同じ理由で連続である。
:::

### 同相写像の例

同相であることを示すには、多くの場合、逆写像を明示的な式で書き、両方の連続性を確かめるのが確実である。以下では Euclid 空間とその部分集合に Euclid 距離の位相（およびその相対位相）を入れる。連続性の確認には :ref[prop-top-real-function-algebra], :ref[prop-top-euclidean-components], :ref[prop-top-restriction] を用いる。

:::example[区間の同相]{#exm-top-interval-homeomorphisms}
1. $a<b$, $c<d$ とする。$\ell(t)=c+\frac{d-c}{b-a}(t-a)$ は一次関数で、逆も一次関数 $\ell^{-1}(s)=a+\frac{b-a}{d-c}(s-c)$ である。$\ell$ は狭義増加なので $(a,b)$ を $(c,d)$ に、$[a,b]$ を $[c,d]$ に、$[a,b)$ を $[c,d)$ に全単射に写す。:ref[prop-top-homeomorphism-restriction] により、端点の種類が同じ有界区間はすべて同相である。
2. $h:(-1,1)\to\R$, $h(t)=t/(1-\abs t)$ は同相写像である。$g(s)=s/(1+\abs s)$ とおくと $\abs{g(s)}<1$ であり、$1-\abs{g(s)}=1/(1+\abs s)$ から $h(g(s))=s$、同様に $1+\abs{h(t)}=1/(1-\abs t)$ から $g(h(t))=t$ である。$h,g$ は分母が零にならない連続関数の商なので連続である。1 と合わせて、すべての開区間 $(a,b)$ は $\R$ と同相である。
3. $\log:(0,\infty)\to\R$ は $\exp$ を逆にもつ連続全単射なので、$(0,\infty)\cong\R$ である。平行移動と $t\mapsto-t$ を合成すれば $(a,\infty)\cong(-\infty,a)\cong\R$ を得る。
4. $k:[0,1)\to[0,\infty)$, $k(t)=t/(1-t)$ は逆 $s\mapsto s/(1+s)$ をもつ同相写像である。したがって半開区間 $[a,b)$, 半直線 $[a,\infty)$ はすべて同相である。

こうして $\R$ の区間は、空集合と一点を除くと、「開区間型」「閉区間型」「半開区間型」の 3 種類に同相で分類される。この 3 種類が互いに同相でないことは、点を一つ除いたときの連結性の違いから [連結性](/math/topology/connectedness/) の章で示される。
:::

:::example[開球と Euclid 空間]{#exm-top-open-ball-euclidean}
$B^n=\{x\in\R^n:\norm x_2<1\}$ とする。$h:B^n\to\R^n$, $g:\R^n\to B^n$ を

$$
h(x)=\frac{x}{1-\norm x_2},\qquad g(y)=\frac{y}{1+\norm y_2}
$$

で定める。$\norm{g(y)}_2=\norm y_2/(1+\norm y_2)<1$ なので $g$ は $B^n$ に値をとる。$1-\norm{g(y)}_2=1/(1+\norm y_2)$ より $h(g(y))=y$、$1+\norm{h(x)}_2=1/(1-\norm x_2)$ より $g(h(x))=x$ である。ノルムは連続（:ref[prop-top-euclidean-components]）で分母は正なので、各成分は連続関数の商として連続であり、$h,g$ は互いに逆な連続写像である。よって $B^n\cong\R^n$ である。平行移動と正のスカラー倍 $x\mapsto a+rx$ は同相写像なので、$\R^n$ の任意の開球は $\R^n$ と同相である。$n=1$ の場合が例 :ref[exm-top-interval-homeomorphisms] の 2 である。
:::

:::example[立体射影]{#exm-top-stereographic-projection}
$S^n=\{(u,s)\in\R^n\times\R:\norm u_2^2+s^2=1\}$ とし、北極 $N=(0,1)$ をとる。$N$ と $(u,s)\in S^n\setminus\{N\}$ を結ぶ直線と超平面 $s=0$ の交点を対応させる写像

$$
\sigma(u,s)=\frac{u}{1-s}\in\R^n
$$

を **立体射影** という。逆写像の候補として

$$
\tau(v)=\Bigl(\frac{2v}{1+\norm v_2^2},\ \frac{\norm v_2^2-1}{\norm v_2^2+1}\Bigr)
$$

をとる。$\norm{2v}_2^2+(\norm v_2^2-1)^2=(\norm v_2^2+1)^2$ なので $\tau(v)\in S^n$ であり、第 2 成分は $1$ より小さいので $\tau(v)\ne N$ である。$1-s=2/(1+\norm v_2^2)$ から $\sigma(\tau(v))=v$ が従う。逆に $(u,s)\in S^n\setminus\{N\}$ と $v=u/(1-s)$ について、$\norm u_2^2=1-s^2$ を用いると

$$
\norm v_2^2=\frac{1-s^2}{(1-s)^2}=\frac{1+s}{1-s},\qquad 1+\norm v_2^2=\frac{2}{1-s}
$$

であるから、$\frac{2v}{1+\norm v_2^2}=u$, $\frac{\norm v_2^2-1}{\norm v_2^2+1}=s$ となり $\tau(\sigma(u,s))=(u,s)$ である。$\sigma$ は座標関数の商で分母 $1-s$ が $S^n\setminus\{N\}$ 上で正、$\tau$ は分母が正の有理式なので、ともに連続である。よって $S^n\setminus\{N\}\cong\R^n$ である。例えば $n=2$ では点 $(2/3,2/3,1/3)\in S^2$ は $\sigma$ により $(1,1)$ に写る。
:::

単位球と閉立方体を対応させる原点からの放射写像は :ref[exm-top-ball-cube] で構成した。

:::example[角度による表示は同相でない]{#exm-top-angle-parametrization}
$e:[0,2\pi)\to S^1$, $e(t)=(\cos t,\sin t)$ は連続な全単射である（三角関数の周期性による。[解析学](/math/analysis/) を参照）。しかし同相写像ではない。$z_n=e(2\pi-1/n)$ とおくと、$\cos,\sin$ の連続性により $z_n\to(1,0)=e(0)$ であるが、$e^{-1}(z_n)=2\pi-1/n$ は $[0,2\pi)$ のどの点にも収束しない（もし $t\in[0,2\pi)$ に収束すれば、$\R$ で $2\pi$ と $t$ の両方に収束することになり極限の一意性に反する）。:ref[thm-top-metric-continuity-tests] により $e^{-1}$ は $(1,0)$ で連続でない。

:ref[prop-top-homeomorphism-test] と照合すると、$e$ は開写像でも閉写像でもないはずである。実際、$[0,1)$ は $[0,2\pi)$ で開であるが、その像は $(1,0)$ を含み、$(1,0)$ に収束する点列 $z_n$ $(n\ge1)$ は像に入らないので、像は開でない。また $[\pi,2\pi)$ は $[0,2\pi)$ で閉であるが、その像は $z_n$ を含み極限 $(1,0)$ を含まないので閉でない。
:::

### 位相的性質

:::definition[位相的性質]{#def-top-topological-property}
位相空間についての性質 $P$ が、$X$ が $P$ をもち $X\cong Y$ ならば $Y$ も $P$ をもつ、という条件を満たすとき、$P$ を :term[位相的性質]{en="topological property" yomi="いそうてきせいしつ"}（位相不変な性質）という。同相な空間に同じ値を与える量を :term[位相不変量]{en="topological invariant" yomi="いそうふへんりょう"}という。
:::

開集合族だけを用いて定義される性質は、同相写像が開集合族の間の全単射を引き起こすので、すべて位相的性質である。具体的に確かめておく。

:::proposition[同相が保つ集合演算と性質]{#prop-top-homeomorphic-invariants}
$h:X\to Y$ を同相写像、$A\subset X$ とする。

1. $h(A^\circ)=h(A)^\circ$, $h(\overline A)=\overline{h(A)}$, $h(\partial A)=\partial h(A)$ である。特に $A$ が開（閉、稠密）であることと $h(A)$ が開（閉、稠密）であることは同値である。
2. $x$ の基本近傍系は $h$ によって $h(x)$ の基本近傍系に、$X$ の開基は $Y$ の開基に写る。
3. 次は位相的性質である：離散であること、密着であること、第一可算性、第二可算性、可分性、Lindelöf 性。また濃度と孤立点の個数は位相不変量である。
:::

:::proof
1. $h$ は開写像なので :ref[prop-top-image-criteria] により $h(A^\circ)\subset h(A)^\circ$ である。同じことを同相写像 $h^{-1}$ と $h(A)$ に適用すると $h^{-1}(h(A)^\circ)\subset A^\circ$、すなわち $h(A)^\circ\subset h(A^\circ)$ である。閉包については、$h$ の連続性から $h(\overline A)\subset\overline{h(A)}$（:ref[thm-top-continuity-tests]）、$h$ が閉写像であることから $\overline{h(A)}\subset h(\overline A)$（:ref[prop-top-image-criteria]）である。全単射は差集合を保つので、$\partial A=\overline A\setminus A^\circ$ から $h(\partial A)=\partial h(A)$ を得る。最後の主張は $A^\circ=A$, $\overline A=A$, $\overline A=X$ という特徴付けと $h(X)=Y$ から従う。

2. $\mathcal B_x$ を $x$ の基本近傍系とする。$h$ は開写像なので $h(N)$ $(N\in\mathcal B_x)$ は $h(x)$ の近傍である（:ref[prop-top-open-map-base-test] の 2）。$h(x)$ の任意の近傍 $M$ について、$h^{-1}(M)$ は $x$ の近傍なので、ある $N\in\mathcal B_x$ が $N\subset h^{-1}(M)$ を満たし、$h(N)\subset M$ となる。開基 $\mathcal B$ については、$Y$ の開集合 $V$ について $h^{-1}(V)=\bigcup_\lambda B_\lambda$ と書けば $V=\bigcup_\lambda h(B_\lambda)$ であり、$h(B_\lambda)$ は開である。

3. 離散・密着は開集合族がそれぞれ冪集合全体・$\{\varnothing,\text{全体}\}$ であることで、全単射 $U\mapsto h(U)$ によって保たれる。第一可算性と第二可算性は 2 と、全単射が集合族の濃度を保つことから従う。可分性は 1 により可算稠密集合が可算稠密集合に写ることから従う。Lindelöf 性については、$Y$ の開被覆 $(V_i)$ の逆像 $(h^{-1}(V_i))$ が $X$ の開被覆なので可算部分被覆 $(h^{-1}(V_{i_k}))_k$ があり、$(V_{i_k})_k$ が $Y$ を覆う。濃度は全単射で保たれ、$x$ が孤立点（$\{x\}$ が開）であることと $h(x)$ が孤立点であることは同値である。
:::

一方、距離を用いて定義される性質の多くは位相的性質ではない。

:::example[位相的でない性質]{#exm-top-non-invariants}
1. **有界性・直径**：$(-1,1)\cong\R$（例 :ref[exm-top-interval-homeomorphisms]）であり、前者は有界、後者は有界でない。
2. **Cauchy 列であること**：$h:(0,1]\to[1,\infty)$, $h(x)=1/x$ は自身と同じ式の逆をもつ同相写像である。$x_n=1/n$ は $(0,1]$ の Cauchy 列であるが、像 $h(x_n)=n$ は Cauchy 列でない。
3. **完備性**（すべての Cauchy 列が収束すること。[完備距離空間](/math/topology/complete-metric/) で扱う）：$[1,\infty)$ の Cauchy 列は実数の完備性（[数列と級数](/math/analysis/sequences-series/)）により $\R$ で収束し、極限は閉集合 $[1,\infty)$ に属する（:ref[prop-top-metric-closed-sequence]）ので完備である。一方 $(0,1]$ では $1/n$ が $(0,1]$ に極限をもたないので完備でない。両者は 2 の $h$ で同相である。
4. **角の有無**：例 :ref[exm-top-ball-cube] により正方形の周と円周は同相である。

これらは距離の構造（一様構造）には依存するが、開集合族だけでは決まらない性質である。
:::

二つの空間が同相でないことを示すには、一方だけがもつ位相的性質を見つければよい。例えば、:ref[prop-top-finite-preorder] で扱った三点集合上の二つの位相、すなわち $a\preceq b$, $a\preceq c$ から定まる位相（開集合は $\varnothing,\{b\},\{c\},\{b,c\},X$）とその順序を逆にした位相（開集合は $\varnothing,\{a\},\{a,b\},\{a,c\},X$）は、ともに開集合を 5 個もつが、孤立点の個数がそれぞれ 2 と 1 なので同相でない。

:::remark[より深い位相不変量]
$\R$ と $\R^2$、あるいは $[0,1)$ と $(0,1)$ が同相でないことは直観的には明らかであるが、ここまでの不変量（濃度・可算公理・可分性など）では区別できない。これらは「一点を除いた空間が連結かどうか」によって [連結性](/math/topology/connectedness/) の章で区別される。さらに一般に、Brouwer の **領域不変性定理**（$U\subset\R^n$ が開集合で $f:U\to\R^n$ が連続な単射ならば $f(U)$ は開集合であり、$f$ は $U$ から $f(U)$ への同相写像である）から、$m\ne n$ なら $\R^m\not\cong\R^n$ であることが従う。この定理の証明にはホモロジー論などの代数的トポロジーの道具が必要であり、本講座の一般位相の範囲を超えるので、ここでは証明しない。
:::

### 埋め込み

部分空間と同相な「コピー」を作る写像が埋め込みである。

:::proposition[埋め込みの判定]{#prop-top-embedding-criteria}
連続な単射 $e:X\to Y$ について次は同値である。

1. $e$ は埋め込みである。
2. $X$ の任意の開集合 $U$ に対し、$Y$ の開集合 $V$ で $e(U)=e(X)\cap V$ となるものが存在する。
3. $e$ を像 $e(X)$ への写像とみなしたものは開写像である。

特に、連続な単射で開写像（閉写像）であるものは、像が開集合（閉集合）である埋め込みである。これを :term[開埋め込み]{en="open embedding" yomi="かいうめこみ"}（:term[閉埋め込み]{en="closed embedding" yomi="へいうめこみ"}）という。
:::

:::proof
$e':X\to e(X)$ を値域を像に制限した写像とすると、:ref[prop-top-restriction] の 2 により $e'$ は連続な全単射である。:ref[prop-top-homeomorphism-test] により、$e'$ が同相写像であることと開写像であることは同値であるから 1 ⇔ 3 である。3 ⇔ 2 は相対位相の定義（$e(X)$ の開集合は $e(X)\cap V$ の形）そのものである。

$e$ が $Y$ への開写像なら、$e(U)$ は $Y$ で開なので $e(X)\cap e(U)=e(U)$ によって 2 が成り立ち、$e(X)$ 自身も開である。$e$ が閉写像なら、:ref[prop-top-open-closed-map-composition] の 3 により $e'$ は閉写像であり、:ref[prop-top-homeomorphism-test] により同相写像である。$e(X)$ は閉である。
:::

:::example[埋め込みの例と反例]
1. 部分空間 $A\subset X$ の包含写像は埋め込みである（$A$ の開集合はちょうど $A\cap V$ の形）。
2. $\R\to\R^2$, $x\mapsto(x,0)$ は閉埋め込みである。逆写像は射影 $(x,y)\mapsto x$ の制限で連続であり、像 $\{y=0\}$ は連続関数 $(x,y)\mapsto y$ の零点集合なので閉である（:ref[cor-top-real-function-sets]）。
3. $\exp:\R\to\R$ は連続な単射で開写像なので、像 $(0,\infty)$ への開埋め込みである。
4. 例 :ref[exm-top-angle-parametrization] の $e:[0,2\pi)\to\R^2$, $t\mapsto(\cos t,\sin t)$ は連続な単射であるが、像 $S^1$ への同相でないので埋め込みでない。単射で連続でも埋め込みとは限らない。
:::

部分空間への引き戻しを使う埋め込みの例は [部分空間・積空間・商空間](/math/topology/constructions/) で扱う。

### 局所同相写像

:::definition[局所同相写像]{#def-top-local-homeomorphism}
写像 $f:X\to Y$ が :term[局所同相写像]{en="local homeomorphism" yomi="きょくしょどうそうしゃぞう"}であるとは、各点 $x\in X$ に開近傍 $U$ が存在して、$f(U)$ が $Y$ の開集合であり、制限 $f|_U:U\to f(U)$ が同相写像となることをいう。
:::

:::proposition[局所同相は連続かつ開である]{#prop-top-local-homeomorphism}
局所同相写像 $f:X\to Y$ は連続な開写像である。さらに $f$ が単射ならば、$f$ は開埋め込みである。
:::

:::proof
各 $x$ に定義の開近傍 $U_x$ を選ぶ。$V\subset Y$ を開集合とすると、$f(U_x)\cap V$ は $f(U_x)$ で開なので、$(f|_{U_x})^{-1}(V)=U_x\cap f^{-1}(V)$ は $U_x$ で開、$U_x$ が開なので $X$ でも開である。$f^{-1}(V)=\bigcup_x(U_x\cap f^{-1}(V))$ は開である（:ref[prop-top-continuity-pasting] の 1 と同じ議論）。

$W\subset X$ を開集合とすると $W\cap U_x$ は $U_x$ で開であり、$f|_{U_x}$ は同相写像なので $f(W\cap U_x)$ は $f(U_x)$ で開、$f(U_x)$ が $Y$ で開なので $Y$ でも開である。$f(W)=\bigcup_xf(W\cap U_x)$ は開である。最後の主張は :ref[prop-top-embedding-criteria] による。
:::

:::example[円周への巻き付け]{#exm-top-circle-covering}
$p:\R\to S^1$, $p(t)=(\cos t,\sin t)$ が局所同相写像であることを示す。$\theta_0\in\R$ を固定し、$J=(\theta_0-\pi,\theta_0+\pi)$, $W=S^1\setminus\{p(\theta_0+\pi)\}$ とおく。$W$ は一点の補集合なので $S^1$ で開である。$(x,y)\in W$ に対し、$\theta_0$ だけ逆回転した点

$$
x'=x\cos\theta_0+y\sin\theta_0,\qquad y'=-x\sin\theta_0+y\cos\theta_0
$$

は $S^1$ 上にあり $(x',y')\ne(-1,0)$ である。そこで

$$
\varphi(x,y)=\theta_0+2\arctan\frac{y'}{1+x'}
$$

とおく。$\varphi$ は分母が正の連続関数で、値は $J$ に入る。$\tan(\theta/2)=\sin\theta/(1+\cos\theta)$ と加法定理から、$t\in J$ について $\varphi(p(t))=t$ が成り立ち、逆に $u=y'/(1+x')$ とおくと $x'^2+y'^2=1$ より $u^2=(1-x')/(1+x')$, $\cos(2\arctan u)=(1-u^2)/(1+u^2)=x'$, $\sin(2\arctan u)=2u/(1+u^2)=y'$ となるので $p(\varphi(x,y))=(x,y)$ である。したがって $p|_J:J\to W$ は連続な逆 $\varphi$ をもつ同相写像であり、$p$ は局所同相写像、特に開写像である。しかし $p$ は単射でないので埋め込みではない。

例 :ref[exm-top-angle-parametrization] の $e=p|_{[0,2\pi)}$ が同相でないこととの違いは、$[0,2\pi)$ が $\R$ の開集合でないことにある。局所同相写像の制限は開集合上でしか局所同相にならない。
:::

## 位相の強弱

### 細かい位相と粗い位相

同じ集合に複数の位相を入れて比べることがある。開集合が多いほど、点を細かく区別でき、収束は起こりにくく、そこから出る写像は連続になりやすい。

:::definition[細かい位相と粗い位相]{#def-top-finer}
集合 $X$ 上の位相 $\mathcal T_1,\mathcal T_2$ が $\mathcal T_1\subset\mathcal T_2$ を満たすとき、$\mathcal T_2$ は $\mathcal T_1$ より :term[細かい]{en="finer" yomi="こまかい"}（強い）、$\mathcal T_1$ は $\mathcal T_2$ より :term[粗い]{en="coarser" yomi="あらい"}（弱い）という。等しい場合も含める。どちらの包含も成り立たないとき、二つの位相は比較できないという。
:::

離散位相は最も細かく、密着位相は最も粗い。「強い・弱い」の用語は文献によって逆の意味に使われることがあるので、包含関係で明記するのが安全である。

:::proposition[恒等写像による判定と連続性の単調性]{#prop-top-identity-topology}
1. 恒等写像 $\id:(X,\mathcal T_2)\to(X,\mathcal T_1)$ が連続であるための必要十分条件は $\mathcal T_1\subset\mathcal T_2$ である。
2. $f:(X,\mathcal T_X)\to(Y,\mathcal T_Y)$ が連続ならば、定義域の位相をより細かい $\mathcal T_X'\supset\mathcal T_X$ に、値域の位相をより粗い $\mathcal T_Y'\subset\mathcal T_Y$ に取り替えても連続である。
3. $\mathcal T_1\subset\mathcal T_2$ のとき、$\mathcal T_2$ で収束するネット（点列）は $\mathcal T_1$ でも同じ点に収束する。
:::

:::proof
1. 恒等写像では $V$ の逆像は $V$ 自身なので、連続性は「$\mathcal T_1$ の元はすべて $\mathcal T_2$ の元である」ことと同じである。
2. $V\in\mathcal T_Y'$ なら $V\in\mathcal T_Y$ なので $f^{-1}(V)\in\mathcal T_X\subset\mathcal T_X'$ である。
3. 1 の恒等写像は連続なので、:ref[prop-top-net-continuity] により収束を保つ。
:::

開基による比較は :ref[lem-top-basis-comparison-finer] による。左右の半開区間が与える位相は比較できない（:ref[exm-top-real-line-topologies]）。

:::example[関数空間上の二つの距離]
$[0,1]$ 上の実連続関数全体 $C[0,1]$ に、上限距離 $d_\infty(f,g)=\max_t\abs{f(t)-g(t)}$ と積分距離 $d_1(f,g)=\int_0^1\abs{f(t)-g(t)}\,\dd t$ を入れる。$d_1\le d_\infty$ なので恒等写像 $(C[0,1],d_\infty)\to(C[0,1],d_1)$ は Lipschitz 連続であり、:ref[prop-top-identity-topology] により $d_\infty$ の位相は $d_1$ の位相より細かい。逆向きの恒等写像は連続でない：$n\ge2$ で $g_n(t)=\max(0,1-n\abs{t-1/2})$ とおくと $d_1(g_n,0)=1/n\to0$ だが $d_\infty(g_n,0)=1$ である。したがって $d_\infty$ の位相は真に細かい（:ref[prop-top-metric-comparison] の点列による判定とも一致する）。
:::

位相の共通部分と、和集合が生成する位相はそれぞれ包含順序の下限と上限である（:ref[prop-top-topology-lattice]）。

### 始位相と終位相

写像の族を与えて、それらが連続になるように位相を定めることができる。写像の定義域に位相を入れる場合は、連続にするために開集合が必要なので「最も粗い」位相を、値域に入れる場合は「最も細かい」位相を選ぶのが自然である。

:::proposition[始位相と終位相]{#prop-top-initial-final}
1. 集合 $X$ と位相空間への写像族 $f_i:X\to Y_i$ $(i\in I)$ を考える。$\{f_i^{-1}(V):i\in I,\ V\text{ は }Y_i\text{ の開集合}\}$ を準開基とする位相 $\mathcal T_{\mathrm{in}}$ は、すべての $f_i$ を連続にする最も粗い位相である。これを :term[始位相]{en="initial topology" yomi="しいそう"}という。写像 $g:Z\to(X,\mathcal T_{\mathrm{in}})$ が連続であるための必要十分条件は、すべての $f_i\circ g$ が連続であることである。
2. 位相空間から集合 $X$ への写像族 $h_i:Y_i\to X$ を考える。$\mathcal T_{\mathrm{fin}}=\{U\subset X:\text{すべての }i\text{ で }h_i^{-1}(U)\text{ が開}\}$ は位相であり、すべての $h_i$ を連続にする最も細かい位相である。これを :term[終位相]{en="final topology" yomi="しゅういそう"}という。写像 $k:(X,\mathcal T_{\mathrm{fin}})\to Z$ が連続であるための必要十分条件は、すべての $k\circ h_i$ が連続であることである。
3. 1（2）の普遍性は始位相（終位相）を特徴付ける。すなわち $X$ 上の位相 $\mathcal T$ が、任意の $g:Z\to X$ について「$g$ が連続 $\iff$ すべての $f_i\circ g$ が連続」を満たすならば $\mathcal T=\mathcal T_{\mathrm{in}}$ である（終位相についても同様）。
:::

:::proof
1. 各 $f_i$ が連続であることは、準開基の元 $f_i^{-1}(V)$ が開であることにほかならない。$\mathcal T$ がすべての $f_i$ を連続にするなら、$\mathcal T$ は準開基のすべての元を含むので、それが生成する $\mathcal T_{\mathrm{in}}$ を含む。普遍性について、$g$ が連続なら $f_i\circ g$ は合成として連続である。逆にすべての $f_i\circ g$ が連続なら、準開基の元の逆像 $g^{-1}(f_i^{-1}(V))=(f_i\circ g)^{-1}(V)$ は開であり、:ref[prop-top-subbase-continuity] により $g$ は連続である。

2. 逆像は和と共通部分を保ち、$h_i^{-1}(\varnothing)=\varnothing$, $h_i^{-1}(X)=Y_i$ であるから、$\mathcal T_{\mathrm{fin}}$ は位相の公理を満たす。定義から各 $h_i$ は連続であり、すべての $h_i$ を連続にする位相の元 $U$ は $h_i^{-1}(U)$ がすべて開なので $\mathcal T_{\mathrm{fin}}$ に属する。普遍性について、$k$ が連続なら $k\circ h_i$ も連続である。逆にすべての $k\circ h_i$ が連続なら、$Z$ の開集合 $W$ について $h_i^{-1}(k^{-1}(W))=(k\circ h_i)^{-1}(W)$ がすべて開なので、$k^{-1}(W)\in\mathcal T_{\mathrm{fin}}$ である。

3. $\mathcal T$ が普遍性をもつとする。$g=\id:(X,\mathcal T)\to(X,\mathcal T)$ は連続なので、性質からすべての $f_i$ は $\mathcal T$ に関して連続であり、1 より $\mathcal T_{\mathrm{in}}\subset\mathcal T$ である。また $g=\id:(X,\mathcal T_{\mathrm{in}})\to(X,\mathcal T)$ について $f_i\circ g=f_i$ は $\mathcal T_{\mathrm{in}}$ に関して連続なので、性質から $g$ は連続であり、:ref[prop-top-identity-topology] により $\mathcal T\subset\mathcal T_{\mathrm{in}}$ である。終位相についても矢印の向きを逆にして同様に示せる。
:::

始位相の典型例は座標射影で定まる積位相、終位相の典型例は全射で定まる商位相である。これらの例は [部分空間・積空間・商空間](/math/topology/constructions/) で計算する。

:::example[実数値関数による始位相]
$X$ 上の実数値関数の族 $\mathcal F$ による始位相を $\mathcal F$ の定める **弱位相** という。

1. $\R$ に、通常の位相で連続な実数値関数すべてによる始位相を入れる。各準開基元 $f^{-1}(V)$ は通常の位相で開なので始位相は通常の位相より粗く、恒等関数 $\id_\R$ が族に含まれるので通常の開集合 $V=\id^{-1}(V)$ はすべて始位相で開である。よって始位相は通常の位相に一致する。
2. Sierpiński 空間 $S$ から $\R$ への連続関数は定数に限る。実際 $f(0)\ne f(1)$ とすると、$f(0)$ を含み $f(1)$ を含まない開区間 $I$ について $f^{-1}(I)=\{0\}$ が開になるが、$S$ で $\{0\}$ は開でない。したがって $C(S)$ による始位相は密着位相であり、$S$ の位相を復元しない。実数値連続関数で位相を復元できる空間（完全正則空間）は [分離公理](/math/topology/separation/) の章で扱う。
:::

### 点列空間

距離空間や第一可算空間では、位相を点列の収束だけで記述できた。この性質を取り出したものが点列空間である。

:::definition[点列開集合と点列空間]{#def-top-sequential-space}
$X$ を位相空間とする。$A\subset X$ が :term[点列開集合]{en="sequentially open set" yomi="てんれつかいしゅうごう"}であるとは、$A$ の点に収束する任意の点列が、ある番号以降すべて $A$ に属することをいう。点列開集合がすべて開集合であるような空間を :term[点列空間]{en="sequential space" yomi="てんれつくうかん"}という。
:::

開集合 $U$ は点列開である（$x\in U$ に収束する点列は、近傍 $U$ にある番号以降入る）。

:::proposition[点列開集合の位相]{#prop-top-sequentialization}
$(X,\mathcal T)$ の点列開集合全体 $\mathcal T_s$ は $\mathcal T$ より細かい位相である。$X_s=(X,\mathcal T_s)$ とおくと、次が成り立つ。

1. $X$ と $X_s$ で、収束する点列とその極限は一致する。
2. 写像 $f:X\to Y$ が点列の極限を保つ（$x_n\to x$ ならば $f(x_n)\to f(x)$）ことと、$f:X_s\to Y$ が連続であることは同値である。
3. $X$ が点列空間であることと $X_s=X$ は同値であり、このとき点列の極限を保つ写像 $f:X\to Y$ は連続である。
:::

:::proof
$\varnothing$ と $X$ は点列開である。$A_i$ が点列開で $x_n\to x\in\bigcup_iA_i$ なら、$x\in A_j$ となる $j$ について $x_n$ はある番号以降 $A_j\subset\bigcup_iA_i$ に属する。$A,B$ が点列開で $x_n\to x\in A\cap B$ なら、$n\ge N_1$ で $x_n\in A$、$n\ge N_2$ で $x_n\in B$ となるので、$n\ge\max(N_1,N_2)$ で $x_n\in A\cap B$ である。よって $\mathcal T_s$ は位相であり、上の注意から $\mathcal T\subset\mathcal T_s$ である。

1. $X_s$ で収束すれば、$\mathcal T\subset\mathcal T_s$ と :ref[prop-top-identity-topology] の 3 により $X$ でも収束する。逆に $X$ で $x_n\to x$ とし、$A\in\mathcal T_s$ が $x$ を含むとすると、点列開の定義により $x_n$ はある番号以降 $A$ に属する。$x$ の $X_s$ での任意の近傍はそのような $A$ を含むので、$X_s$ でも $x_n\to x$ である。

2. $f$ が点列の極限を保つとし、$V\subset Y$ を開集合とする。$x_n\to x\in f^{-1}(V)$ なら $f(x_n)\to f(x)\in V$ なので、ある番号以降 $f(x_n)\in V$、すなわち $x_n\in f^{-1}(V)$ である。よって $f^{-1}(V)$ は点列開であり、$f:X_s\to Y$ は連続である。逆に $f:X_s\to Y$ が連続なら、$X$ で $x_n\to x$ のとき 1 により $X_s$ でも $x_n\to x$ であり、:ref[prop-top-net-continuity] により $f(x_n)\to f(x)$ である。

3. 前半は定義そのものである。後半は 2 による。
:::

:::proposition[第一可算空間は点列空間である]{#prop-top-first-countable-sequential-space}
第一可算空間は点列空間である。特に距離空間は点列空間である。
:::

:::proof
$X$ を第一可算とし、$A\subset X$ が点列開であるが開でないとする。すると $A$ の点 $x$ で、$x$ のどの近傍も $A$ に含まれないもの、すなわち $x$ のどの近傍も $X\setminus A$ と交わるものがある。よって $x\in\overline{X\setminus A}$ であり、:ref[prop-top-countable-basis] により $X\setminus A$ の点列 $(x_n)$ で $x$ に収束するものがある。これは一つの項も $A$ に属さないので、$A$ が点列開であることに反する。距離空間は各点で半径 $1/n$ の開球が可算な基本近傍系をなすので第一可算である。
:::

:::example[補可算位相の点列化]
例 :ref[exm-top-cocountable-sequential-continuity] の補可算位相をもつ非可算集合 $X$ では、収束する点列はある番号以降極限に等しい。したがって任意の部分集合が点列開であり、$X_s$ は離散空間である。$X$ 自身は離散でないので点列空間ではなく、:ref[prop-top-first-countable-sequential-space] により第一可算でもない。:ref[prop-top-sequentialization] の 2 に従えば、$X$ から任意の空間への任意の写像が点列の極限を保つ。これが例 :ref[exm-top-cocountable-sequential-continuity] の現象の説明である。
:::

## 演習問題

::::exercise{level=basic}
$f:\R\to\R$, $f(x)=x^2$ による $(1,4)$ の逆像と $[1,4]$ の逆像を求め、それぞれ開・閉であることを確かめよ。

:::solution
$1<x^2<4$ は $1<\abs x<2$ と同値なので、$f^{-1}((1,4))=(-2,-1)\cup(1,2)$ であり、開区間の和として開である。同様に $f^{-1}([1,4])=[-2,-1]\cup[1,2]$ であり、閉区間二つの和として閉である。$f$ は連続（:ref[prop-top-real-function-algebra]）なので :ref[thm-top-continuity-tests] からも従う。
:::
::::

::::exercise{level=basic}
三点離散空間 $D_3$ から二点密着空間 $I_2$ への連続写像の個数、および $I_2$ から $D_3$ への連続写像の個数を求めよ。

:::solution
$I_2$ の開集合は $\varnothing,I_2$ だけであり、その逆像は $\varnothing,D_3$ なので、$D_3\to I_2$ の写像はすべて連続である。各点の像を 2 通りから選べるので $2^3=8$ 個である。

逆向きでは、$f:I_2\to D_3$ が連続なら各一点集合 $\{y\}$ の逆像が $\varnothing$ または $I_2$ でなければならない。$f$ が二点を異なる点 $y\ne y'$ に写すと $f^{-1}(\{y\})$ が一点集合になって開でないので、$f$ は定数写像である。定数写像は連続なので、答えは $3$ 個である。
:::
::::

::::exercise{level=basic}
$f:\R\to\{0,1\}$（値域は離散）を $x\ge0$ で $1$、$x<0$ で $0$ と定める。$f$ が連続となる点と不連続となる点をすべて求めよ。

:::solution
$x>0$ なら開区間 $(0,\infty)$ 上で $f$ は定数なので、$f(x)$ の任意の近傍の逆像は $(0,\infty)$ を含み、$x$ の近傍である。$x<0$ も同様に $(-\infty,0)$ を用いて連続である。$x=0$ では $f(0)=1$ の近傍 $\{1\}$ の逆像 $[0,\infty)$ が $0$ の近傍でない（$0$ を含むどの開区間も負の数を含む）ので不連続である。不連続点は $0$ だけである。
:::
::::

::::exercise{level=basic}
定数写像 $f:\R\to\R$, $f(x)=3$ は連続か、開写像か、閉写像か。

:::solution
開集合の逆像は $\varnothing$ または $\R$ なので連続である。閉集合の像は $\varnothing$ または $\{3\}$ であり、いずれも閉なので閉写像である。開集合 $\R$ の像 $\{3\}$ は開でないので開写像でない。
:::
::::

::::exercise{level=basic}
$(0,1)$ から $(2,5)$ への同相写像と、$[0,1)$ から $[0,\infty)$ への同相写像を、逆写像とともに具体的に与えよ。

:::solution
$\ell(t)=2+3t$ は $(0,1)$ を $(2,5)$ に全単射に写し、逆は $\ell^{-1}(s)=(s-2)/3$ である。ともに一次関数で連続である。

$k(t)=t/(1-t)$ は $[0,1)$ 上で分母が正の連続関数で値は $[0,\infty)$ に入り、$m(s)=s/(1+s)$ は $[0,\infty)$ 上の連続関数で値は $[0,1)$ に入る。$k(m(s))=\dfrac{s/(1+s)}{1/(1+s)}=s$, $m(k(t))=\dfrac{t/(1-t)}{1/(1-t)}=t$ であるから、$k$ は同相写像である。
:::
::::

::::exercise{level=basic}
$f:X\to Y$ を連続写像、$D\subset X$ を稠密な部分集合とする。$f(D)$ は $f(X)$（相対位相）で稠密であることを示せ。

:::solution
:ref[thm-top-continuity-tests] の 4 により $f(X)=f(\overline D)\subset\overline{f(D)}$ である。$f(X)$ における $f(D)$ の閉包は $f(X)\cap\overline{f(D)}$ である（相対位相の閉集合は $f(X)\cap F$ の形であり、$f(D)$ を含む最小のものは $F=\overline{f(D)}$ のときに得られる）から、これは $f(X)$ に等しい。よって $f(D)$ は $f(X)$ で稠密である。例えば $\Q$ は $\R$ で稠密なので、連続関数 $f:\R\to\R$ の値は有理数での値で決まる。
:::
::::

::::exercise{level=advanced}
距離空間の間の写像 $f:(X,d)\to(Y,\rho)$ が開写像であるための必要十分条件は、任意の $x\in X$ と $r>0$ に対して $s>0$ が存在し、$B_\rho(f(x),s)\subset f(B_d(x,r))$ となることであることを示せ。

:::solution
:ref[prop-top-open-map-base-test] の 2 により、$f$ が開写像であることは、各 $x$ の各近傍の像が $f(x)$ の近傍であることと同値である。距離空間では開球 $B_d(x,r)$ $(r>0)$ が $x$ の基本近傍系をなし、近傍の像は含まれる開球の像を含むので、条件は各 $x,r$ について $f(B_d(x,r))$ が $f(x)$ の近傍であることと同値である。$f(B_d(x,r))$ が $f(x)$ の近傍であることは、ある $s>0$ で $B_\rho(f(x),s)\subset f(B_d(x,r))$ となることと同値である（:ref[prop-top-continuity-metric-agree] の証明の第 1 段落）。
:::
::::

::::exercise{level=advanced}
Sorgenfrey 直線 $\R_S$（$[a,b)$ を開基とする位相）から通常の $\R$ への写像 $f$ が連続であるための必要十分条件は、$f$ が各点で右連続（$\lim_{y\downarrow x}f(y)=f(x)$）であることを示せ。さらに床関数 $\lfloor x\rfloor$ は $\R_S\to\R$ として連続であるが、天井関数 $\lceil x\rceil$ は連続でないことを示せ。

:::solution
$[x,x+\delta)$ $(\delta>0)$ は $\R_S$ における $x$ の基本近傍系である（$x$ を含む基本開集合 $[a,b)$ は $[x,b)$ を含む）。したがって $f$ が $x$ で連続であることは、任意の $\eps>0$ に対して $\delta>0$ があり、$x\le y<x+\delta$ ならば $\abs{f(y)-f(x)}<\eps$ となることと同値であり、これは右連続性の定義そのものである。:ref[thm-top-continuity-tests] により、各点で右連続であることと連続性は同値である。

床関数は $k\le x<k+1$ で値 $k$ をとるので、$x\in[k,k+1)$ に対し $[x,k+1)$ 上で定数であり、右連続である。天井関数は $x=0$ で $\lceil0\rceil=0$ だが $0<y<1$ で $\lceil y\rceil=1$ なので右連続でなく、$\R_S\to\R$ として連続でない。
:::
::::

::::exercise{level=advanced}
$f:X\to Y$ を連続な開写像とする。任意の $B\subset Y$ について $f^{-1}(\overline B)=\overline{f^{-1}(B)}$ が成り立つことを示せ。また、開写像の仮定を外すと成り立たない例を挙げよ。

:::solution
$\supset$ は :ref[thm-top-continuity-tests] の 5 である。$\subset$ を示す。$x\in f^{-1}(\overline B)$ とし、$U$ を $x$ を含む開集合とする。$f(U)$ は $f(x)$ を含む開集合であり、$f(x)\in\overline B$ なので :ref[prop-top-closure-test] により $f(U)\cap B\ne\varnothing$ である。$f(u)\in B$ となる $u\in U$ をとれば $u\in U\cap f^{-1}(B)$ である。$U$ は任意なので $x\in\overline{f^{-1}(B)}$ である。

反例：定数写像 $f\equiv0:\R\to\R$ と $B=(0,1)$ について $f^{-1}(\overline B)=f^{-1}([0,1])=\R$ だが $\overline{f^{-1}(B)}=\overline\varnothing=\varnothing$ である。
:::
::::

::::exercise{level=advanced}
$f:X\to Y$ を写像とし、$X\times Y$ には $U\times V$（$U,V$ は開集合）の形の集合全体を開基とする位相を入れる（積位相。[部分空間・積空間・商空間](/math/topology/constructions/) を参照）。$f$ が連続ならば $j:X\to X\times Y$, $j(x)=(x,f(x))$ は埋め込みであることを示せ。

:::solution
$U\times V$ の形の集合の族は $X\times Y$ を覆い、$(U\times V)\cap(U'\times V')=(U\cap U')\times(V\cap V')$ なので開基の条件（:ref[prop-top-basis-test]）を満たす。$j$ は第 1 成分から $x$ が読めるので単射である。開基の元の逆像は $j^{-1}(U\times V)=U\cap f^{-1}(V)$ であり、開なので $j$ は連続である（:ref[prop-top-subbase-continuity]）。$p(x,y)=x$ とすると $p^{-1}(U)=U\times Y$ は開なので $p$ は連続であり、$p\circ j=\id_X$ である。像 $j(X)$ 上の逆写像は連続写像 $p$ の制限なので、$j$ は埋め込みである。
:::
::::

::::exercise{level=advanced}
$f,g:X\to\R$ が下半連続ならば $f+g$ も下半連続であることを示せ。また $f$ が下半連続で $c>0$ なら $cf$ は下半連続、$c<0$ なら $cf$ は上半連続であることを示せ。

:::solution
$f(x)+g(x)>a$ であることは、$a-g(x)<q<f(x)$ を満たす有理数 $q$ が存在することと同値である（有理数の稠密性）。したがって

$$
\{f+g>a\}=\bigcup_{q\in\Q}\bigl(\{f>q\}\cap\{g>a-q\}\bigr)
$$

であり、右辺の各項は二つの開集合の共通部分なので開、その和も開である。$c>0$ なら $\{cf>a\}=\{f>a/c\}$ は開である。$c<0$ なら $\{cf<a\}=\{f>a/c\}$ は開なので $cf$ は上半連続である。
:::
::::

::::exercise{level=advanced}
$X$ を無限集合とし、補有限位相を入れる。写像 $f:X\to X$ が連続であるための必要十分条件は、$f$ が定数写像であるか、すべての $y\in X$ について $f^{-1}(y)$ が有限集合であることを示せ。また $X$ の任意の全単射は同相写像であることを示せ。

:::solution
閉集合は有限集合と $X$ なので、$f$ が連続であることは、任意の有限集合 $F$ について $f^{-1}(F)$ が有限または $X$ であることと同値である（:ref[thm-top-continuity-tests] の 3）。

$f$ が定数 $c$ なら $f^{-1}(F)$ は $c\in F$ のとき $X$、そうでなければ $\varnothing$ で、連続である。すべてのファイバーが有限なら、$f^{-1}(F)=\bigcup_{y\in F}f^{-1}(y)$ は有限個の有限集合の和として有限であり、連続である。

逆に $f$ が連続で定数でないとする。$y\in X$ について $f^{-1}(y)$ は有限または $X$ であるが、$X$ なら $f$ は定数 $y$ となって仮定に反するので、有限である。

全単射 $h$ はファイバーが一点なので連続であり、逆写像 $h^{-1}$ も全単射なので連続である。よって同相写像である。
:::
::::

## 確認テスト

::::quiz{answer="3" tol="1e-6"}
Sierpiński 空間 $S=\{0,1\}$（開集合は $\varnothing,\{1\},S$）から自身への連続写像はいくつあるか。

:::explanation
$S$ から $S$ への写像は 4 つある。二つの定数写像と恒等写像は連続である。交換写像 $\sigma$ は $\sigma^{-1}(\{1\})=\{0\}$ が開でないので連続でない。したがって答えは $3$ である。
:::
::::

::::quiz
任意の連続写像 $f:X\to Y$ について常に成り立つものをすべて選べ。

- [ ] 閉集合の像は閉である。
- [x] 閉集合の逆像は閉である。
- [ ] 開集合の像は開である。
- [x] 任意の $A\subset X$ について $f(\overline A)\subset\overline{f(A)}$ である。

:::explanation
逆像は補集合と交換するので、閉集合の逆像は閉である（:ref[thm-top-continuity-tests] の 3）。閉包の包含は同定理の 4 である。像についての主張は誤りであり、平面の第 1 射影と閉集合 $\{xy=1\}$（像 $\R\setminus\{0\}$ は閉でない）、平方写像と開区間 $(-1,1)$（像 $[0,1)$ は開でない）が反例である。正解は (B) と (D) である。
:::
::::

::::quiz{answer="2" tol="1e-6"}
$S^2\subset\R^3$ の北極 $N=(0,0,1)$ からの立体射影 $\sigma(u,s)=u/(1-s)$ による点 $(2/3,2/3,1/3)$ の像を $v\in\R^2$ とする。$\norm v_2^2$ の値を求めよ。

:::explanation
$(2/3)^2+(2/3)^2+(1/3)^2=1$ なので点は $S^2$ 上にある。$1-s=2/3$ なので $v=(2/3,2/3)/(2/3)=(1,1)$ であり、$\norm v_2^2=2$ である。逆写像 $\tau(v)=\bigl(2v/(1+\norm v_2^2),(\norm v_2^2-1)/(\norm v_2^2+1)\bigr)=\bigl((2/3,2/3),1/3\bigr)$ で元の点に戻ることも確かめられる。
:::
::::

::::quiz
集合 $X$ 上の位相が $\mathcal T_1\subsetneq\mathcal T_2$ を満たすとき、連続な恒等写像の向きとして正しいものを選べ。

- [x] $(X,\mathcal T_2)\to(X,\mathcal T_1)$ だけが連続である。
- [ ] $(X,\mathcal T_1)\to(X,\mathcal T_2)$ だけが連続である。
- [ ] 両方向とも連続である。

:::explanation
恒等写像 $(X,\mathcal T)\to(X,\mathcal T')$ が連続であることは $\mathcal T'\subset\mathcal T$ と同値である（:ref[prop-top-identity-topology]）。$\mathcal T_1\subset\mathcal T_2$ なので細かい方から粗い方への恒等写像は連続であり、包含が真なので逆向きは連続でない。正解は (A) である。
:::
::::

::::quiz
位相的性質（同相写像で保たれる性質）であるものをすべて選べ。

- [x] 可分であること
- [x] 第一可算であること
- [ ] 距離空間として有界であること
- [ ] 距離空間として完備であること
- [x] 孤立点をちょうど 2 個もつこと

:::explanation
可分性・第一可算性・孤立点の個数は開集合族だけで定義されるので、:ref[prop-top-homeomorphic-invariants] により保たれる。有界性は $(-1,1)\cong\R$ によって、完備性は $(0,1]\cong[1,\infty)$（$x\mapsto1/x$）によって保たれないことがわかる（例 :ref[exm-top-non-invariants]）。正解は (A), (B), (E) である。
:::
::::

::::quiz{answer="1" tol="1e-6"}
$X=\{a,b,c\}$ に開集合 $\varnothing,\{a\},\{a,b\},X$ を与える。$X$ から自身への同相写像はいくつあるか。

:::explanation
同相写像 $h$ は開集合の全体を自身に写す全単射を引き起こす。一点の開集合は $\{a\}$ だけなので $h(\{a\})=\{a\}$、すなわち $h(a)=a$ である。二点の開集合は $\{a,b\}$ だけなので $h(\{a,b\})=\{a,b\}$ となり $h(b)=b$、したがって $h(c)=c$ である。よって同相写像は恒等写像だけであり、答えは $1$ である。
:::
::::

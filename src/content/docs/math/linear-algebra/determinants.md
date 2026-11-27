---
title: 行列式
description: 置換の符号、多重線形交代形式としての行列式の特徴付けと Leibniz の公式、積公式・転置・基本変形、小行列式と階数、余因子展開と逆行列、Cramer の公式
notation:
  - symbol: 'S_n'
    meaning: n 文字の置換全体（n 次対称群）
  - symbol: '\sgn(\sigma)'
    meaning: 置換 σ の符号
  - symbol: '(a_1\ a_2\ \cdots\ a_k)'
    meaning: 長さ k の巡回置換
  - symbol: '\det A'
    meaning: 正方行列 A の行列式
  - symbol: '\det T'
    meaning: 線形変換 T の行列式
  - symbol: 'A_{\widehat i,\widehat j}'
    meaning: A から第 i 行と第 j 列を除いた行列
  - symbol: 'C_{ij}'
    meaning: A の (i, j) 余因子
  - symbol: '\operatorname{adj}(A)'
    meaning: A の余因子行列（余因子を並べた行列の転置）
---

## この章の目標

- 置換の符号を定義し、互換・巡回置換への分解と符号の関係を理解する。
- 列についての多重線形性・交代性・正規化条件によって行列式を特徴付け、Leibniz の公式を導く。
- 行列式の積公式・転置不変性・基本変形による変化を証明し、行列式が $0$ でないことと正則性が同値であることを示す。
- 余因子展開によって行列式を計算し、余因子行列による逆行列の公式と Cramer の公式を導く。
- 小行列式による階数の特徴付け、Vandermonde の行列式、線形変換の行列式などの応用を身につける。

2 元連立一次方程式 $ax+by=e$, $cx+dy=f$ から $y$ を消去すると $(ad-bc)x=ed-bf$ を得る。したがって $ad-bc\ne0$ ならば解はただ 1 つに定まる。この量 $ad-bc$ を $n$ 元の場合に一般化したものが行列式である。

$ad-bc$ には次の性質がある。列ベクトル $a_1=\transp{(a,c)}$, $a_2=\transp{(b,d)}$ の関数 $D(a_1,a_2)=ad-bc$ とみると、(i) 一方を固定すれば他方について線形であり、(ii) $a_1=a_2$ ならば $0$ であり、(iii) $D(e_1,e_2)=1$ である。$K=\R$ のとき、$\abs{D(a_1,a_2)}$ は $a_1,a_2$ の張る平行四辺形の面積であり、符号は $a_1$ から $a_2$ への回転の向きを表す。本章では、性質 (i)–(iii) をもつ関数が $n$ 次元でもただ 1 つ存在することを示し、それを行列式と定義する。その準備として、まず列の並べ替えに伴う符号を調べる。

本章の結果はすべて任意の体 $K$ 上で成り立つ。さらに、多くの結果は $\Z$ や多項式環 $K[t]$ を成分とする行列に対しても成り立つ（:ref[prop-linalg-det-ring]）。これは次章で特性多項式を扱う際に用いる。

## 置換と符号

### 置換

:::definition[置換と互換]{#def-linalg-permutation}
$\{1,\ldots,n\}$ からそれ自身への全単射を $n$ 文字の :term[置換]{en="permutation" yomi="ちかん"} といい、その全体を $S_n$ と書く。置換の積 $\sigma\tau$ は合成 $\sigma\circ\tau$ とし、右の $\tau$ を先に作用させる：$(\sigma\tau)(i)=\sigma(\tau(i))$。恒等写像を $\id$ と書く。相異なる $p,q$ を入れ替え、他の文字を固定する置換を :term[互換]{en="transposition" yomi="ごかん"} といい $(p\ q)$ と書く。
:::

置換 $\sigma$ は、値を並べた列 $[\sigma(1),\ldots,\sigma(n)]$（**一行表示**）や、上下 2 段に並べた表示

$$
\sigma=\begin{pmatrix}1&2&\cdots&n\\ \sigma(1)&\sigma(2)&\cdots&\sigma(n)\end{pmatrix}
$$

で表される。右から互換 $(p\ q)$ を掛けると、$(\sigma(p\ q))(p)=\sigma(q)$, $(\sigma(p\ q))(q)=\sigma(p)$ なので、一行表示の第 $p$ 項と第 $q$ 項が入れ替わる。左から掛けると、値 $p$ と値 $q$ が入れ替わる。

:::proposition[置換全体の性質]{#prop-linalg-symmetric-group}
1. $\sigma,\tau\in S_n$ ならば $\sigma\tau\in S_n$, $\sigma^{-1}\in S_n$ であり、積は結合法則を満たし、$\id\sigma=\sigma\id=\sigma$, $\sigma\sigma^{-1}=\sigma^{-1}\sigma=\id$ である。
2. $S_n$ の元の個数は $n!$ である。
3. $\tau\in S_n$ を固定すると、$\sigma\mapsto\sigma\tau$, $\sigma\mapsto\tau\sigma$, $\sigma\mapsto\sigma^{-1}$ はいずれも $S_n$ からそれ自身への全単射である。
:::

:::proof
1. 全単射の合成と逆写像は全単射であり、写像の合成は結合的である（[集合と論理](/math/sets-logic/)）。

2. $\{1,\dots,k\}$ から $\{1,\dots,n\}$ への単射の個数が $n(n-1)\cdots(n-k+1)$ であることを $k$ に関する帰納法で示す。$k=0$ のときは空写像 1 個である。$\{1,\dots,k\}$ からの単射は、$\{1,\dots,k-1\}$ からの単射 $f$ と、$f$ の値以外の $n-(k-1)$ 個の値から選んだ $k$ の行き先との組に 1 対 1 に対応するので、個数は $n(n-1)\cdots(n-k+2)\cdot(n-k+1)$ となる。有限集合 $\{1,\dots,n\}$ からそれ自身への単射は全単射なので、$k=n$ として $\abs{S_n}=n!$ を得る。

3. それぞれ $\sigma\mapsto\sigma\tau^{-1}$, $\sigma\mapsto\tau^{-1}\sigma$, $\sigma\mapsto\sigma^{-1}$ が逆写像である。
:::

1 は $S_n$ が合成について**群**をなすことを述べている。$S_n$ を $n$ 次**対称群**という（群としての詳しい性質は [群論](/math/group-theory/) で扱う）。

### 巡回置換への分解

:::definition[巡回置換]{#def-linalg-cycle}
相異なる $a_1,\dots,a_k\in\{1,\dots,n\}$（$k\ge2$）に対し、$a_1\mapsto a_2\mapsto\cdots\mapsto a_k\mapsto a_1$ と順に送り、他の文字を固定する置換を長さ $k$ の :term[巡回置換]{en="cycle" yomi="じゅんかいちかん"} といい、$(a_1\ a_2\ \cdots\ a_k)$ と書く。互換は長さ $2$ の巡回置換である。2 つの巡回置換が共通の文字を含まないとき、**互いに素**であるという。
:::

互いに素な巡回置換 $\gamma,\delta$ は可換である。実際、各文字は高々一方によってしか動かされず、$\gamma$ が動かす文字の集合は $\gamma$ で自分自身に移るので、$\gamma\delta$ と $\delta\gamma$ はどの文字でも同じ値をとる。

:::theorem[巡回置換分解]{#thm-linalg-cycle-decomposition}
任意の $\sigma\in S_n$ は、互いに素な巡回置換の積として表される（$\sigma=\id$ のときは空積とする）。この表し方は、因子の順序を除いて一意である。
:::

:::proof
文字 $a$ に対し、列 $a,\sigma(a),\sigma^2(a),\dots$ は有限集合の元なので、どこかで同じ値が現れる。$\sigma^i(a)=\sigma^j(a)$（$0\le i<j$）となる組のうち $j$ が最小のものをとると $i=0$ である（$i\ge1$ ならば $\sigma$ の単射性から $\sigma^{i-1}(a)=\sigma^{j-1}(a)$ となり、$j$ の最小性に反する）。この $j$ を $\ell_a$ と書くと、$a,\sigma(a),\dots,\sigma^{\ell_a-1}(a)$ は相異なり、$\sigma^{\ell_a}(a)=a$ である。この集合 $O_a\defeq\{\sigma^i(a)\mid i\ge0\}$ を $a$ の**軌道**と呼ぶ。

$b\in O_a$ ならば $O_b=O_a$ である（$b=\sigma^i(a)$ ならば $O_b\subset O_a$ は明らかで、$a=\sigma^{\ell_a-i}(b)$ より $O_a\subset O_b$）。よって 2 つの軌道は一致するか交わらないかのどちらかであり、$\{1,\dots,n\}$ は互いに交わらない軌道に分割される。長さ $2$ 以上の各軌道 $O_a$ について巡回置換 $\gamma_a\defeq(a\ \sigma(a)\ \cdots\ \sigma^{\ell_a-1}(a))$ を作ると、これらは互いに素であり、その積はどの文字 $b$ でも $\sigma(b)$ に一致する（$b$ を含む軌道の巡回置換だけが $b$ を動かす）。

一意性：$\sigma=\gamma_1\cdots\gamma_m$ を互いに素な巡回置換への任意の分解とする。$\gamma_i$ に現れる文字 $b$ について、$\sigma$ は $b$ に $\gamma_i$ だけが作用するように働くので $\sigma(b)=\gamma_i(b)$ であり、$\gamma_i$ の文字は $\sigma$ の 1 つの軌道をなす。よって各 $\gamma_i$ は上で作った $\gamma_a$ のいずれかに一致する。
:::

例えば一行表示 $[3,5,4,1,2]$ の置換は $1\mapsto3\mapsto4\mapsto1$, $2\mapsto5\mapsto2$ なので $(1\ 3\ 4)(2\ 5)$ である。

:::proposition[巡回置換の互換分解]{#prop-linalg-cycle-transpositions}
$(a_1\ a_2\ \cdots\ a_k)=(a_1\ a_k)(a_1\ a_{k-1})\cdots(a_1\ a_2)$ が成り立つ。したがって、任意の置換は互換の積で表される。
:::

:::proof
$k$ に関する帰納法で示す。$k=2$ は自明である。$k\ge3$ とし、$\gamma=(a_1\ \cdots\ a_{k-1})$ について成り立つとする。$(a_1\ a_k)\gamma$ を計算すると、$a_i$（$i\le k-2$）は $\gamma$ で $a_{i+1}\ne a_1,a_k$ に移り、その後固定されるので $a_{i+1}$ へ、$a_{k-1}$ は $a_1$ に移った後 $a_k$ へ、$a_k$ は $\gamma$ で固定された後 $a_1$ へ移る。他の文字は固定される。よって $(a_1\ a_k)\gamma=(a_1\ \cdots\ a_k)$ である。後半は :ref[thm-linalg-cycle-decomposition] による。
:::

互換への分解は一意ではない。例えば $(1\ 3)=(1\ 2)(2\ 3)(1\ 2)$ であり、右辺は $1\mapsto2\mapsto3$, $3\mapsto2\mapsto1$, $2\mapsto1\mapsto2$ と確かめられる。しかし、次に示すように、用いる互換の個数の**偶奇**は一定である。

### 転倒数と符号

:::definition[転倒数と符号]{#def-linalg-permutation-sign}
$\sigma\in S_n$ について、$i<j$ かつ $\sigma(i)>\sigma(j)$ となる組 $(i,j)$ を $\sigma$ の :term[転倒]{en="inversion" yomi="てんとう"} と呼び、その個数を $N(\sigma)$ とする。置換の :term[符号]{en="sign of a permutation" yomi="ふごう"} を

$$
\sgn(\sigma)\defeq(-1)^{N(\sigma)}
$$

と定める。$\sgn(\sigma)=1$ のとき $\sigma$ を :term[偶置換]{en="even permutation" yomi="ぐうちかん"}、$-1$ のとき :term[奇置換]{en="odd permutation" yomi="きちかん"} という。
:::

例えば $[3,1,4,2]$ の転倒は、位置の組で書くと $(1,2),(1,4),(3,4)$ の 3 個なので、符号は $-1$ である。転倒数による定義は、置換がどのように作られたかによらず値の列だけから計算できる点で便利である。

:::lemma[隣り合う位置の交換]{#lem-linalg-adjacent-inversions}
$s_i=(i\ \ i+1)$（$1\le i<n$）とおく。$N(\sigma s_i)$ は、$\sigma(i)<\sigma(i+1)$ のとき $N(\sigma)+1$、$\sigma(i)>\sigma(i+1)$ のとき $N(\sigma)-1$ である。また任意の置換は、このような隣接互換 $s_i$ の積で表される。
:::

:::proof
$\sigma s_i$ の一行表示は、$\sigma$ の一行表示の第 $i$ 項と第 $i+1$ 項を入れ替えたものである。$k<i$ に対し、組 $(k,i),(k,i+1)$ のうち転倒であるものの個数は、位置 $i,i+1$ にある 2 つの値の集合だけで決まるので、入れ替えで変わらない。$k>i+1$ に対する組 $(i,k),(i+1,k)$ も同様である。位置 $i,i+1$ の両方を含まない組は影響を受けない。変わるのは組 $(i,i+1)$ だけであり、入れ替えによってこの組の大小が逆転する。これで前半を得る。

恒等置換以外の一行表示には、隣り合う 2 項が減少している箇所がある。実際、そのような箇所がなければ一行表示は狭義単調増加であり、$1,\dots,n$ をこの順に並べたもの、すなわち恒等置換しかない。減少箇所 $i$ を選んで右から $s_i$ を掛けると、前半より転倒数は $1$ 減る。転倒数は $0$ 以上なので、有限回で $\sigma s_{i_1}\cdots s_{i_h}=\id$ となる。$s_i^{-1}=s_i$ なので $\sigma=s_{i_h}\cdots s_{i_1}$ である。
:::

:::theorem[符号の積法則と互換分解の偶奇]{#thm-linalg-sign-product}
任意の $\sigma,\tau\in S_n$ に対し

$$
\sgn(\sigma\tau)=\sgn(\sigma)\sgn(\tau),\qquad
\sgn(\sigma^{-1})=\sgn(\sigma)
$$

が成り立つ。また互換の符号は $-1$ であり、$\sigma$ を $k$ 個の互換の積で表したとき $\sgn(\sigma)=(-1)^k$ である。特に、互換の個数の偶奇は分解の選び方によらない。
:::

:::proof
:ref[lem-linalg-adjacent-inversions] により $\tau=s_{i_1}\cdots s_{i_h}$ と書ける。同じ補題により、右から隣接互換を 1 つ掛けるたびに転倒数は $\pm1$ 変化し、符号は反転する。よって $\sigma$ に右から $s_{i_1},\dots,s_{i_h}$ を順に掛けると

$$
\sgn(\sigma\tau)=(-1)^h\sgn(\sigma)
$$

である。$\sigma=\id$（転倒数 $0$）とすれば $\sgn(\tau)=(-1)^h$ なので、積法則が従う。$1=\sgn(\id)=\sgn(\sigma)\sgn(\sigma^{-1})$ と $\sgn(\sigma)=\pm1$ から $\sgn(\sigma^{-1})=\sgn(\sigma)$ である。

次に互換 $(p\ q)$（$p<q$）の転倒を数える。一行表示は $1,\dots,n$ の第 $p$ 項と第 $q$ 項を入れ替えたもの $[\dots,q,\dots,j,\dots,p,\dots]$ である。転倒となる組は、$(p,q)$ と、$p<j<q$ を満たす各 $j$ についての $(p,j)$ と $(j,q)$ だけである（他の組では値の大小が位置の大小と一致する）。よって転倒数は $1+2(q-p-1)$ で奇数であり、互換の符号は $-1$ である。積法則を $k$ 回用いれば最後の主張を得る。
:::

:::corollary[巡回置換の符号・偶置換の個数]{#cor-linalg-cycle-sign}
1. 長さ $k$ の巡回置換の符号は $(-1)^{k-1}$ である。
2. $\sigma$ を互いに素な巡回置換に分解したとき、固定点も長さ $1$ の巡回とみなして巡回の個数を $c(\sigma)$ とすると $\sgn(\sigma)=(-1)^{n-c(\sigma)}$ である。
3. $n\ge2$ ならば、$S_n$ の偶置換と奇置換はそれぞれ $n!/2$ 個ずつある。
:::

:::proof
1. :ref[prop-linalg-cycle-transpositions] より長さ $k$ の巡回置換は $k-1$ 個の互換の積である。
2. 長さ $\ell_1,\dots,\ell_c$ の巡回（$\sum\ell_i=n$）に分解されるならば、1 と積法則より $\sgn(\sigma)=\prod_i(-1)^{\ell_i-1}=(-1)^{n-c}$。
3. :ref[prop-linalg-symmetric-group] の 3 より $\sigma\mapsto\sigma(1\ 2)$ は $S_n$ の全単射であり、積法則により偶置換と奇置換を入れ替える。よって両者の個数は等しい。
:::

偶置換全体を $A_n$ と書く（[群論](/math/group-theory/) の交代群）。例えば $[3,5,4,1,2]=(1\ 3\ 4)(2\ 5)$ については、転倒が $(1,4),(1,5),(2,3),(2,4),(2,5),(3,4),(3,5)$ の 7 個で符号は $-1$ であり、巡回置換分解からも $(+1)\cdot(-1)=-1$、あるいは $(-1)^{5-2}=-1$ と求まる。

## 多重線形交代形式としての行列式

### 多重線形交代形式

:::definition[多重線形形式と交代性]{#def-linalg-alternating-form}
$V$ を $K$ 上のベクトル空間とし、$k\ge1$ とする。写像 $F:V^k=V\times\cdots\times V\to K$ が、他の $k-1$ 個の引数を固定するごとに残る 1 つの引数について線形であるとき、$F$ を（$k$ 重の）:term[多重線形形式]{en="multilinear form" yomi="たじゅうせんけいけいしき"} という。多重線形形式 $F$ が、2 つの引数が等しいとき（$x_p=x_q$, $p\ne q$）常に $F(x_1,\dots,x_k)=0$ となる性質を :term[交代性]{en="alternating property" yomi="こうたいせい"} といい、このとき $F$ を**交代形式**（交代的な多重線形形式）という。
:::

$k$ 個のベクトルの組 $x=(x_1,\dots,x_k)$ と $\sigma\in S_k$ に対し、$x^\sigma\defeq(x_{\sigma(1)},\dots,x_{\sigma(k)})$ と書く（引数の並べ替え）。$x^\sigma$ の第 $i$ 成分は $x_{\sigma(i)}$ なので、$(x^\sigma)^\tau$ の第 $i$ 成分は $x^\sigma$ の第 $\tau(i)$ 成分 $x_{\sigma(\tau(i))}$ であり、

$$
(x^\sigma)^\tau=x^{\sigma\tau} \label{eq-linalg-permute-arguments}
$$

が成り立つ。特に互換 $\tau=(p\ q)$ に対し $y^\tau$ は組 $y$ の第 $p$ 成分と第 $q$ 成分を入れ替えたものである。

:::lemma[交代形式と引数の並べ替え]{#lem-linalg-alternating-permute}
$F:V^k\to K$ を交代形式とする。

1. 2 つの引数を入れ替えると、$F$ の値は $-1$ 倍になる。
2. 任意の $\sigma\in S_k$ と $x\in V^k$ に対し $F(x^\sigma)=\sgn(\sigma)F(x)$。
3. $x_1,\dots,x_k$ が一次従属ならば $F(x_1,\dots,x_k)=0$。
:::

:::proof
1. 第 $p$ 引数と第 $q$ 引数（$p<q$）以外を固定し、$G(u,v)\defeq F(\dots,u,\dots,v,\dots)$（$u$ が第 $p$ 引数、$v$ が第 $q$ 引数）とおく。$G$ は双線形で $G(w,w)=0$ を満たすので

$$
0=G(u+v,u+v)=G(u,u)+G(u,v)+G(v,u)+G(v,v)=G(u,v)+G(v,u)
$$

となり、$G(v,u)=-G(u,v)$ である。

2. :ref[prop-linalg-cycle-transpositions] により $\sigma=\tau_1\tau_2\cdots\tau_h$（$\tau_i$ は互換）と書ける。:eqref[eq-linalg-permute-arguments] を繰り返し用いると $x^\sigma=(\cdots((x^{\tau_1})^{\tau_2})\cdots)^{\tau_h}$ であり、各段階は 2 つの成分の入れ替えなので、1 により $F$ の値は $-1$ 倍される。よって $F(x^\sigma)=(-1)^hF(x)=\sgn(\sigma)F(x)$（:ref[thm-linalg-sign-product]）。

3. ある $x_i$ が他の一次結合 $x_i=\sum_{j\ne i}c_jx_j$ で書ける。第 $i$ 引数についての線形性により $F(x_1,\dots,x_k)=\sum_{j\ne i}c_jF(\dots,x_j,\dots)$ であり、右辺の各項は第 $i$ 引数と第 $j$ 引数がともに $x_j$ なので $0$ である。
:::

:::remark
1 の性質（**歪対称性**）を交代性の定義とする本もある。交代性から歪対称性が従うことは上で示した。逆に $1+1\ne0$ ならば、歪対称性から $F(\dots,w,\dots,w,\dots)=-F(\dots,w,\dots,w,\dots)$ より $2F(\dots)=0$、したがって交代性が従う。しかし $K=\F_2$ では $-1=1$ なので、例えば $\F_2^2$ 上の $F(x,y)=x_1y_1$ は歪対称（すなわち対称）だが $F(e_1,e_1)=1\ne0$ で交代的でない。任意の体で議論するため、本サイトでは交代性を「等しい引数で $0$」として定義する。
:::

### 行列式の存在と一意性

:::theorem[行列式の存在と一意性]{#thm-linalg-determinant-unique}
交代形式 $F:(K^n)^n\to K$ で $F(e_1,\ldots,e_n)=1$ を満たすものがただ 1 つ存在する。$A=(a_1\ \cdots\ a_n)=(a_{ij})\in M_n(K)$ の列を代入した値 $F(a_1,\dots,a_n)$ を $A$ の :term[行列式]{en="determinant" yomi="ぎょうれつしき"} といい、$\det A$ または $\abs{A}$ と書く。その値は次の **Leibniz の公式**で与えられる。

$$
\det A=\sum_{\sigma\in S_n}\sgn(\sigma)\prod_{j=1}^n a_{\sigma(j),j}=\sum_{\sigma\in S_n}\sgn(\sigma)\,a_{\sigma(1),1}a_{\sigma(2),2}\cdots a_{\sigma(n),n}.
\label{eq-linalg-leibniz}
$$
:::

:::proof
**一意性**：$F$ を条件を満たす交代形式とする。第 $j$ 列を $a_j=\sum_{i}a_{ij}e_i$ と展開し、各引数についての線形性を順に用いると

$$
F(a_1,\ldots,a_n)=
\sum_{i_1=1}^n\cdots\sum_{i_n=1}^n
a_{i_1,1}a_{i_2,2}\cdots a_{i_n,n}\,F(e_{i_1},\ldots,e_{i_n})
$$

となる。添字 $i_1,\dots,i_n$ に重複があれば、交代性により $F(e_{i_1},\dots,e_{i_n})=0$ である。重複がなければ $j\mapsto i_j$ は置換 $\sigma\in S_n$ であり、$(e_{\sigma(1)},\dots,e_{\sigma(n)})=(e_1,\dots,e_n)^\sigma$ なので、:ref[lem-linalg-alternating-permute] の 2 より

$$
F(e_{\sigma(1)},\ldots,e_{\sigma(n)})=\sgn(\sigma)F(e_1,\ldots,e_n)=\sgn(\sigma)
$$

である。よって $F(a_1,\dots,a_n)$ は :eqref[eq-linalg-leibniz] の右辺に等しく、$F$ は一意である。

**存在**：:eqref[eq-linalg-leibniz] の右辺を $D(a_1,\dots,a_n)$ とおき、条件を確かめる。

- 多重線形性：各項 $\sgn(\sigma)\prod_ja_{\sigma(j),j}$ は第 $k$ 列の成分をちょうど 1 つ（$a_{\sigma(k),k}$）だけ含むので、他の列を固定すると第 $k$ 列の線形関数である。その和も線形である。
- 交代性：第 $p$ 列と第 $q$ 列（$p<q$）が等しいとし、$\tau=(p\ q)$ とおく。:ref[cor-linalg-cycle-sign] の証明のとおり $\sigma\mapsto\sigma\tau$ は偶置換全体 $A_n$ から奇置換全体への全単射なので、$D=\sum_{\sigma\in A_n}\bigl(\prod_ja_{\sigma(j),j}-\prod_ja_{\sigma\tau(j),j}\bigr)$ と書ける。$j\ne p,q$ では $\sigma\tau(j)=\sigma(j)$、また $a_{\sigma\tau(p),p}=a_{\sigma(q),p}=a_{\sigma(q),q}$, $a_{\sigma\tau(q),q}=a_{\sigma(p),q}=a_{\sigma(p),p}$（第 $p$ 列と第 $q$ 列が等しい）なので、2 つの積は一致し、$D=0$ である。
- 正規化：$A=I_n$ のとき $a_{\sigma(j),j}=\delta_{\sigma(j),j}$ なので、$\sigma=\id$ の項だけが $1$ で、他の項は $0$ である。

$n=1$ のときは $S_1=\{\id\}$ であり、$\det(a)=a$ である。
:::

:::corollary[交代形式は行列式の定数倍]{#cor-linalg-alternating-forms}
任意の交代形式 $G:(K^n)^n\to K$ について、すべての $A=(a_1\ \cdots\ a_n)$ に対し

$$
G(a_1,\ldots,a_n)=G(e_1,\ldots,e_n)\det A \label{eq-linalg-alternating-multiple}
$$

が成り立つ。すなわち $(K^n)^n$ 上の交代形式全体は 1 次元のベクトル空間であり、$\det$ がその基底である。
:::

:::proof
一意性の証明の展開をそのまま $G$ に適用すると、$G(e_{\sigma(1)},\dots,e_{\sigma(n)})=\sgn(\sigma)G(e_1,\dots,e_n)$ より $G(a_1,\dots,a_n)=G(e_1,\dots,e_n)\sum_\sigma\sgn(\sigma)\prod_ja_{\sigma(j),j}$ を得る。$G(e_1,\dots,e_n)=0$ の場合も同じ計算が成り立つ。
:::

:::example[2 次と 3 次の行列式]{#exm-linalg-det-small}
$n=2$ の置換は恒等置換（偶）と $(1\ 2)$（奇）だけなので

$$
\det\begin{pmatrix}a&b\\c&d\end{pmatrix}=ad-bc.
$$

$n=3$ では偶置換 3 個、奇置換 3 個の計 6 項があり、各積を行番号の順に並べ直すと

$$
\begin{aligned}
\det A={}&a_{11}a_{22}a_{33}+a_{12}a_{23}a_{31}+a_{13}a_{21}a_{32}\\
         &-a_{13}a_{22}a_{31}-a_{12}a_{21}a_{33}-a_{11}a_{23}a_{32}
\end{aligned}
$$

となる（**Sarrus の方法**）。例えば $A=\begin{pmatrix}1&2&0\\0&1&1\\1&0&1\end{pmatrix}$ では、零でない項は第 1 項の $1$ と第 2 項の $2$ だけなので $\det A=3$ である。Sarrus の方法のような「斜めの積の和」は $n=3$ までしか正しくない。$n=4$ では $24$ 項が必要である。
:::

:::remark
$K=\R$ のとき、$\abs{\det A}$ は列ベクトル $a_1,\dots,a_n$ の張る平行体 $\{\sum_jt_ja_j\mid 0\le t_j\le1\}$ の $n$ 次元体積に等しく、$\det A$ の符号は $(a_1,\dots,a_n)$ が標準基底と同じ「向き」をもつかどうかを表す。符号付き体積は、1 つの辺を $c$ 倍すれば $c$ 倍になり、1 つの辺に他の辺の定数倍を加える（平行移動的にずらす）操作で変わらず、2 つの辺が一致すれば潰れて $0$ になり、単位立方体で $1$ である。これらは行列式を特徴付ける性質（多重線形性・交代性・正規化）そのものである。体積の厳密な扱いは積分論（[解析学](/math/analysis/) の重積分と変数変換公式）で行う。
:::

## 行列式の基本性質

### 積と転置

:::theorem[積の行列式]{#thm-linalg-det-product}
$A,B\in M_n(K)$ に対し、

$$
\det(AB)=\det A\det B. \label{eq-linalg-det-product}
$$
:::

:::proof
$A$ を固定し、$B$ の列 $b_1,\ldots,b_n$ の関数として

$$
G(b_1,\ldots,b_n)\defeq\det(Ab_1\ \cdots\ Ab_n)
$$

を考える。$AB$ の第 $j$ 列は $Ab_j$ なので $G(b_1,\dots,b_n)=\det(AB)$ である。$b_j\mapsto Ab_j$ が線形であり、$\det$ が各列について線形であることから $G$ は多重線形である。2 つの $b_j$ が等しければ対応する $Ab_j$ も等しいので、$G$ は交代的である。さらに $Ae_j$ は $A$ の第 $j$ 列だから $G(e_1,\ldots,e_n)=\det A$ である。:eqref[eq-linalg-alternating-multiple] を $G$ に適用すれば $\det(AB)=\det A\det B$ を得る。
:::

:::proposition[転置の行列式]{#prop-linalg-det-transpose}
任意の $A\in M_n(K)$ に対して $\det(\transp{A})=\det A$ である。したがって

$$
\det A=\sum_{\sigma\in S_n}\sgn(\sigma)\prod_{i=1}^na_{i,\sigma(i)} \label{eq-linalg-leibniz-rows}
$$

であり、行列式は**行**についても多重線形かつ交代的である。
:::

:::proof
$\transp{A}$ の $(i,j)$ 成分は $a_{ji}$ なので、Leibniz の公式より $\det\transp{A}=\sum_\sigma\sgn(\sigma)\prod_ja_{j,\sigma(j)}$ である。積の因子を $k=\sigma(j)$ で並べ直すと $\prod_ja_{j,\sigma(j)}=\prod_ka_{\sigma^{-1}(k),k}$ である。$\sgn(\sigma^{-1})=\sgn(\sigma)$（:ref[thm-linalg-sign-product]）であり、$\sigma\mapsto\sigma^{-1}$ は $S_n$ の全単射（:ref[prop-linalg-symmetric-group]）なので、和を $\rho=\sigma^{-1}$ についてとり直せば

$$
\det\transp{A}=\sum_{\rho\in S_n}\sgn(\rho)\prod_ka_{\rho(k),k}=\det A
$$

である。:eqref[eq-linalg-leibniz-rows] はこの計算の途中の式である。$A$ の行は $\transp{A}$ の列なので、行についての性質は列についての性質から従う。
:::

### 基本変形と正則性

:::proposition[基本変形と行列式]{#prop-linalg-det-operations}
$A\in M_n(K)$, $c\in K$ とする。

1. 2 つの行（または列）を入れ替えると、行列式は $-1$ 倍になる。
2. 1 つの行（または列）を $c$ 倍すると、行列式は $c$ 倍になる。特に $\det(cA)=c^n\det A$。
3. 1 つの行（または列）に、他の行（または列）の定数倍を加えても、行列式は変わらない。
4. 基本行列について $\det P_{ij}=-1$, $\det Q_i(c)=c$, $\det R_{ij}(c)=1$ である。
:::

:::proof
列について示せば、行については :ref[prop-linalg-det-transpose] から従う。1 は :ref[lem-linalg-alternating-permute] の 1、2 はその列についての線形性による（$cA$ はすべての列を $c$ 倍したもの）。3：第 $j$ 列を $a_j+ca_k$（$k\ne j$）に替えると、第 $j$ 列についての線形性により、元の行列式と「第 $j$ 列を $a_k$ に替えた行列の行列式の $c$ 倍」との和になる。後者は等しい列を 2 つもつので $0$ である。4 は $I_n$ に 1–3 を適用すればよい。
:::

$E$ を基本行列とすると、:ref[lem-linalg-elementary-left-mult] と 4 から $\det(EA)=\det E\det A$ である（これは積公式の特別な場合でもある）。

:::corollary[行列式と正則性]{#cor-linalg-det-invertible}
$A\in M_n(K)$ について、$A$ が正則であることと $\det A\ne0$ は同値である。$A$ が正則ならば $\det(A^{-1})=(\det A)^{-1}$ である。したがって、斉次方程式 $Ax=0$ が自明でない解をもつための必要十分条件は $\det A=0$ である。
:::

:::proof
$A$ が正則ならば、積公式から $1=\det I_n=\det(AA^{-1})=\det A\det(A^{-1})$ であり、$\det A\ne0$ と逆数の公式を得る。$A$ が正則でなければ、:ref[thm-linalg-invertible-equiv] より $\rank A<n$ であり、:ref[prop-linalg-independence-rank] より $A$ の列は一次従属である。:ref[lem-linalg-alternating-permute] の 3 より $\det A=0$ である。最後の主張は :ref[thm-linalg-invertible-equiv] の (1)$\Leftrightarrow$(4) による。
:::

:::remark
行列式は**線形ではない**。一般に $\det(A+B)\ne\det A+\det B$ であり（例：$A=B=I_2$ で $4\ne2$）、$\det(cA)=c^n\det A\ne c\det A$ である。行列式が線形なのは、1 つの列（または行）以外を固定したときに限る。
:::

:::proposition[三角行列とブロック三角行列]{#prop-linalg-det-block-triangular}
上三角行列（$i>j$ ならば $a_{ij}=0$）および下三角行列の行列式は、対角成分の積である。また $B\in M_r(K)$, $C\in M_{n-r}(K)$ と $r\times(n-r)$ 行列 $H$ に対して

$$
\det\begin{pmatrix}B&H\\O&C\end{pmatrix}=\det B\det C
$$

が成り立ち、下ブロック三角 $\begin{pmatrix}B&O\\H'&C\end{pmatrix}$ についても同じ公式が成り立つ。
:::

:::proof
上三角行列 $A$ の Leibniz の公式で $0$ でない項 $\prod_ja_{\sigma(j),j}$ は、すべての $j$ について $\sigma(j)\le j$ を満たす。両辺の総和はどちらも $1+\cdots+n$ なので、すべて等号でなければならず、$\sigma=\id$ である。よって $\det A=\prod_ia_{ii}$。

ブロックの場合、$M=\begin{pmatrix}B&H\\O&C\end{pmatrix}$ の成分 $m_{ij}$ は $i>r\ge j$ のとき $0$ である。$0$ でない項では、$j\le r$ について $\sigma(j)\le r$ でなければならない。$\sigma$ は単射なので $\sigma$ は $\{1,\dots,r\}$ をそれ自身の上に全単射に移し、したがって $\{r+1,\dots,n\}$ もそれ自身の上に移す。このような $\sigma$ は $\sigma_1\in S_r$ と $\{r+1,\dots,n\}$ の置換 $\sigma_2$ の組に 1 対 1 に対応する。$i\le r<j$ の組 $(i,j)$ は $\sigma(i)\le r<\sigma(j)$ なので転倒でなく、$N(\sigma)=N(\sigma_1)+N(\sigma_2)$、すなわち $\sgn(\sigma)=\sgn(\sigma_1)\sgn(\sigma_2)$ である。成分の積も $B$ の成分の積と $C$ の成分の積に分かれるので、和をとれば $\det B\det C$ を得る。下ブロック三角の場合は転置して :ref[prop-linalg-det-transpose] を用いる。
:::

帰納法により、対角ブロックが $B_1,\dots,B_k$ のブロック上三角行列の行列式は $\det B_1\cdots\det B_k$ である。

:::example[消去法による計算]{#exm-linalg-det-elimination}
$A=\begin{pmatrix}0&2&1\\1&1&0\\2&0&3\end{pmatrix}$ の行列式を求める。第 1 行と第 2 行を入れ替え、新しい第 3 行から第 1 行の $2$ 倍を引き、次に第 3 行に第 2 行を加えると

$$
\begin{pmatrix}1&1&0\\0&2&1\\2&0&3\end{pmatrix}
\longrightarrow
\begin{pmatrix}1&1&0\\0&2&1\\0&-2&3\end{pmatrix}
\longrightarrow
\begin{pmatrix}1&1&0\\0&2&1\\0&0&4\end{pmatrix}
$$

となる。最後の行列式は $1\cdot2\cdot4=8$ であり、行の入れ替えを 1 回行ったので $\det A=-8$ である。

もう 1 つの例として、Pascal の三角形を並べた行列を考える。各行から前の行を引く操作を繰り返すと

$$
\det\begin{pmatrix}1&1&1&1\\1&2&3&4\\1&3&6&10\\1&4&10&20\end{pmatrix}
=\det\begin{pmatrix}1&1&1&1\\0&1&2&3\\0&2&5&9\\0&3&9&19\end{pmatrix}
=\det\begin{pmatrix}1&1&1&1\\0&1&2&3\\0&0&1&3\\0&0&3&10\end{pmatrix}
=\det\begin{pmatrix}1&1&1&1\\0&1&2&3\\0&0&1&3\\0&0&0&1\end{pmatrix}=1
$$

となる（第 1 段は第 2–4 行から第 1 行を、第 2 段は第 3, 4 行から第 2 行の $2$ 倍、$3$ 倍を、第 3 段は第 4 行から第 3 行の $3$ 倍を引いた）。
:::

:::remark
Leibniz の公式は $n!$ 個の項からなり、各項に $n-1$ 回の掛け算が必要である。$n=20$ で $20!\approx2.4\times10^{18}$ 項となり、実用的な計算はできない。一方、行基本変形で上三角行列にする方法（Gauss の消去法）は約 $n^3/3$ 回の乗除算で済む。Leibniz の公式は理論的な道具（成分の多項式であること、性質の証明）として使い、数値の計算には消去法を用いるのが原則である。
:::

### 線形変換の行列式

:::definition[線形変換の行列式]{#def-linalg-operator-determinant}
$V$ を $n$ 次元ベクトル空間（$n\ge1$）、$T\in\End(V)$ とする。$V$ の基底 $\mathcal B$ を 1 つ選び、

$$
\det T\defeq\det[T]_{\mathcal B}
$$

と定める。これは基底の選び方によらない。実際、別の基底 $\mathcal B'$ について $[T]_{\mathcal B'}=P^{-1}[T]_{\mathcal B}P$（:ref[prop-linalg-similarity-equivalence]）なので、積公式より $\det[T]_{\mathcal B'}=\det(P^{-1})\det[T]_{\mathcal B}\det P=\det[T]_{\mathcal B}$ である。
:::

:::proposition[線形変換の行列式の性質]{#prop-linalg-operator-det}
$S,T\in\End(V)$ について、$\det(S\circ T)=\det S\det T$, $\det\id_V=1$ であり、$T$ が同型写像であることと $\det T\ne0$ は同値である。
:::

:::proof
同じ基底 $\mathcal B$ を用いれば $[S\circ T]_{\mathcal B}=[S]_{\mathcal B}[T]_{\mathcal B}$, $[\id_V]_{\mathcal B}=I_n$ である（:ref[cor-linalg-invertible-map-matrix]）。積公式と、$T$ が同型写像 $\iff$ $[T]_{\mathcal B}$ が正則（:ref[cor-linalg-invertible-map-matrix]）$\iff$ $\det[T]_{\mathcal B}\ne0$（:ref[cor-linalg-det-invertible]）から従う。
:::

例えば $M_2(K)$ 上の転置 $X\mapsto\transp X$ の行列式は、:ref[exm-linalg-transpose-matrix] の表現行列（$I_4$ の第 2 列と第 3 列を入れ替えたもの）の行列式 $-1$ である。相似な行列は同じ行列式をもつので、行列式は（階数・トレースと同様に）相似の不変量である。

### 小行列式と階数

:::definition[小行列式]{#def-linalg-minor}
$A\in M_{m,n}(K)$ と $1\le k\le\min\{m,n\}$ に対し、行番号 $i_1<\cdots<i_k$ と列番号 $j_1<\cdots<j_k$ を選び、それらの行と列の交差する成分からなる $k$ 次正方行列の行列式を、$A$ の $k$ 次の :term[小行列式]{en="minor" yomi="しょうぎょうれつしき"} という。
:::

:::theorem[小行列式による階数の特徴付け]{#thm-linalg-rank-minor}
$A\in M_{m,n}(K)$, $A\ne O$ とする。$\rank A$ は、$0$ でない $k$ 次小行列式が存在するような $k$ の最大値に等しい。
:::

:::proof
$r=\rank A$ とおく。

$0$ でない $k$ 次小行列式があるとし、それが行 $I=\{i_1,\dots,i_k\}$ と列 $J=\{j_1,\dots,j_k\}$ によるものとする。この $k$ 次正方行列 $A_{I,J}$ は正則（:ref[cor-linalg-det-invertible]）なので、その列は一次独立である。$A$ の第 $j_1,\dots,j_k$ 列の間に一次関係 $\sum_lc_la_{j_l}=0$ があれば、その第 $i_1,\dots,i_k$ 成分だけを見ると $A_{I,J}$ の列の一次関係になるので $c_l=0$ である。よって $A$ の $k$ 個の列が一次独立であり、:ref[thm-linalg-row-column-rank] より $k\le r$ である。

逆に、:ref[thm-linalg-row-column-rank] より $A$ の列空間は $r$ 次元なので、:ref[thm-linalg-basis-exists] により $A$ の列から一次独立な $r$ 個を選べる。それらを並べた $m\times r$ 行列 $B$ は $\rank B=r$ をもつので、:ref[thm-linalg-row-column-rank] より $B$ の行空間も $r$ 次元であり、$B$ の行から一次独立な $r$ 個を選べる。それらを並べた $r$ 次正方行列 $C$ は階数 $r$ なので正則であり、$\det C\ne0$ である。$\det C$ は $A$ の $r$ 次小行列式である。
:::

例えば $A=\begin{pmatrix}1&2&3\\2&4&6\\1&0&1\end{pmatrix}$ では、第 2 行が第 1 行の 2 倍なので $\det A=0$ であるが、第 1, 3 行と第 1, 2 列からなる小行列式は $\det\begin{pmatrix}1&2\\1&0\end{pmatrix}=-2\ne0$ なので、$\rank A=2$ である。小行列式は成分の多項式なので、この定理は「$K=\R$ のとき、$\rank A\ge k$ という条件は成分を少し動かしても保たれる」ことなどを示すのに使われる。

## 余因子展開と逆行列

### 余因子展開

:::definition[余因子]{#def-linalg-cofactor}
$A\in M_n(K)$（$n\ge2$）から第 $i$ 行と第 $j$ 列を除いた $(n-1)$ 次正方行列を $A_{\widehat i,\widehat j}$ と書く。

$$
C_{ij}\defeq (-1)^{i+j}\det A_{\widehat i,\widehat j}
$$

を $A$ の $(i,j)$ :term[余因子]{en="cofactor" yomi="よいんし"} という。$n=1$ のときは $C_{11}\defeq1$ とする（$0\times0$ 行列の行列式を $1$ と約束する）。
:::

符号 $(-1)^{i+j}$ は、$(1,1)$ 成分を $+$ として市松模様に $+,-$ が並ぶ。

:::theorem[余因子展開]{#thm-linalg-cofactor-expansion}
$A\in M_n(K)$ の任意の第 $i$ 行、および任意の第 $j$ 列に沿って、次の :term[余因子展開]{en="cofactor expansion" yomi="よいんしてんかい"}（Laplace 展開）が成り立つ。

$$
\det A=\sum_{k=1}^n a_{ik}C_{ik},\qquad
\det A=\sum_{k=1}^n a_{kj}C_{kj}. \label{eq-linalg-cofactor-expansion}
$$
:::

:::proof
$n=1$ では両式とも $a_{11}=a_{11}\cdot1$ なので、$n\ge2$ とする。まず列展開を示す。第 $j$ 列を $a_j=\sum_{i}a_{ij}e_i$ と展開すると、第 $j$ 列についての線形性より

$$
\det A=\sum_{i=1}^na_{ij}\det A^{(i,j)}
$$

である。ここで $A^{(i,j)}$ は $A$ の第 $j$ 列を $e_i$ に置き換えた行列である。$\det A^{(i,j)}=C_{ij}$ を示せばよい。

$A^{(i,j)}$ において、第 $j$ 列を隣り合う列との入れ替えを $j-1$ 回行って第 1 列に移し（他の列の順序は保たれる）、次に第 $i$ 行を隣り合う行との入れ替えを $i-1$ 回行って第 1 行に移す。得られる行列 $M$ は、第 1 列が $e_1$（元の第 $j$ 列 $e_i$ の $1$ が第 1 行に移る）であり、第 1 行・第 1 列を除いた部分は、$A$ から第 $i$ 行と第 $j$ 列を除いた行列 $A_{\widehat i,\widehat j}$ に等しい（行と列の順序は元のまま）。すなわち

$$
M=\begin{pmatrix}1&*\\ O&A_{\widehat i,\widehat j}\end{pmatrix}
$$

であり、:ref[prop-linalg-det-block-triangular] より $\det M=\det A_{\widehat i,\widehat j}$ である。入れ替えを計 $(i-1)+(j-1)$ 回行ったので、:ref[prop-linalg-det-operations] の 1 より $\det A^{(i,j)}=(-1)^{i+j-2}\det M=C_{ij}$ である。

行展開は、$\transp{A}$ の第 $i$ 列に沿った列展開から得られる。$\transp{A}$ の $(k,i)$ 成分は $a_{ik}$ であり、$(\transp{A})_{\widehat k,\widehat i}=\transp{(A_{\widehat i,\widehat k})}$ なので、$\transp{A}$ の $(k,i)$ 余因子は $(-1)^{k+i}\det A_{\widehat i,\widehat k}=C_{ik}$ である。よって $\det A=\det\transp{A}=\sum_ka_{ik}C_{ik}$ である。
:::

$0$ の多い行や列に沿って展開すれば、計算すべき小行列式の個数が減る。

:::example[零の多い列に沿う展開]{#exm-linalg-cofactor-expansion}
次の行列では第 2 列に $0$ でない成分が 1 つしかないので、第 2 列に沿って展開すると

$$
\det\begin{pmatrix}1&0&4\\2&0&1\\3&2&0\end{pmatrix}
=2\cdot(-1)^{3+2}\det\begin{pmatrix}1&4\\2&1\end{pmatrix}
=-2\cdot(1-8)=14
$$

である。
:::

:::example[三重対角行列]{#exm-linalg-tridiagonal}
対角成分が $2$、そのすぐ上と下の成分が $-1$、他の成分が $0$ の $n$ 次正方行列 $A$ の行列式を $D_n$ とする（この行列は 2 階微分 $-\dv{^2}{x^2}$ を差分で近似するときに現れる）。$n\ge3$ とし、第 1 行に沿って展開する。第 1 行で $0$ でない成分は $a_{11}=2$ と $a_{12}=-1$ だけであり、$(1,1)$ 余因子は $C_{11}=D_{n-1}$ である。$A_{\widehat1,\widehat2}$ の第 1 列は $\transp{(-1,0,\dots,0)}$ であり、その第 1 行・第 1 列を除いた部分は $n-2$ 次の同じ形の行列なので、第 1 列に沿って展開すると $\det A_{\widehat1,\widehat2}=-D_{n-2}$、すなわち $C_{12}=(-1)^{1+2}(-D_{n-2})=D_{n-2}$ である。よって

$$
D_n=2C_{11}+(-1)C_{12}=2D_{n-1}-D_{n-2}\qquad(n\ge3)
$$

である。$D_1=2$, $D_2=3$ から、帰納法により $D_n=n+1$ を得る（$2n-(n-1)=n+1$）。
:::

### Vandermonde の行列式

:::theorem[Vandermonde の行列式]{#thm-linalg-vandermonde}
$x_1,\dots,x_n\in K$ に対し

$$
\det\begin{pmatrix}1&1&\cdots&1\\ x_1&x_2&\cdots&x_n\\ x_1^2&x_2^2&\cdots&x_n^2\\ \vdots&\vdots&&\vdots\\ x_1^{n-1}&x_2^{n-1}&\cdots&x_n^{n-1}\end{pmatrix}=\prod_{1\le i<j\le n}(x_j-x_i). \label{eq-linalg-vandermonde}
$$

特にこの行列が正則であるための必要十分条件は、$x_1,\dots,x_n$ が相異なることである。
:::

:::proof
左辺を $V(x_1,\dots,x_n)$ とおき、$n$ に関する帰納法で示す。$n=1$ では両辺とも $1$（空積）である。$n\ge2$ とする。第 $n$ 行から第 $n-1$ 行の $x_1$ 倍を引き、次に第 $n-1$ 行から第 $n-2$ 行の $x_1$ 倍を引き、…、最後に第 2 行から第 1 行の $x_1$ 倍を引く（下の行から順に行うので、各段階で引く行はまだ変更されていない）。行列式は変わらず（:ref[prop-linalg-det-operations] の 3）、第 $i$ 行（$i\ge2$）第 $j$ 列の成分は $x_j^{i-1}-x_1x_j^{i-2}=x_j^{i-2}(x_j-x_1)$ となる。第 1 列は $e_1$ になるので第 1 列に沿って展開し、第 $j$ 列から共通因子 $x_j-x_1$ をくくり出すと（列についての線形性）

$$
V(x_1,\dots,x_n)=\prod_{j=2}^n(x_j-x_1)\cdot V(x_2,\dots,x_n)
$$

となる。帰納法の仮定より主張が従う。最後の主張は :ref[cor-linalg-det-invertible] と、体では積が $0$ になるのはいずれかの因子が $0$ のときに限ることによる。
:::

例えば $V(a,b,c)=(b-a)(c-a)(c-b)$ である。Vandermonde 行列の転置は、:ref[exm-linalg-lagrange-interpolation] の写像 $p\mapsto(p(c_0),\dots,p(c_n))$ の基底 $(1,t,\dots,t^n)$ と標準基底に関する表現行列である。この定理は、相異なる点での値から多項式がただ 1 つ決まることの、行列式による別証明を与える。

### 余因子行列と逆行列

:::definition[余因子行列]{#def-linalg-adjugate}
$A\in M_n(K)$ の余因子を並べた行列 $(C_{ij})$ の転置

$$
\operatorname{adj}(A)\defeq\transp{(C_{ij})}
$$

を $A$ の :term[余因子行列]{en="adjugate matrix" yomi="よいんしぎょうれつ"} という。$\operatorname{adj}(A)$ の $(i,j)$ 成分は $C_{ji}$ である（添字の順序に注意）。
:::

:::theorem[余因子行列と逆行列]{#thm-linalg-adjugate-identity}
すべての $A\in M_n(K)$ に対して

$$
A\operatorname{adj}(A)=\operatorname{adj}(A)A=(\det A)I_n. \label{eq-linalg-adjugate}
$$

したがって $\det A\ne0$ ならば $A^{-1}=\dfrac{1}{\det A}\operatorname{adj}(A)$ である。
:::

:::proof
$A\operatorname{adj}(A)$ の $(i,j)$ 成分は $\sum_ka_{ik}C_{jk}$ である。$i=j$ ならば、これは第 $j$ 行に沿った余因子展開なので $\det A$ に等しい。$i\ne j$ のとき、$A$ の第 $j$ 行を第 $i$ 行で置き換えた行列 $A'$ を考える。$A'$ の第 $j$ 行に関する余因子 $C'_{jk}$ は第 $j$ 行を除いて作るので $C_{jk}$ と等しく、$A'$ の第 $j$ 行に沿った展開は $\sum_ka_{ik}C_{jk}$ である。$A'$ は等しい 2 行をもつので $\det A'=0$ であり、この和は $0$ である。

$\operatorname{adj}(A)A$ の $(i,j)$ 成分は $\sum_kC_{ki}a_{kj}$ である。$i=j$ ならば第 $i$ 列に沿った展開で $\det A$ に等しい。$i\ne j$ ならば、第 $i$ 列を第 $j$ 列で置き換えた行列の第 $i$ 列に沿った展開であり、等しい 2 列をもつので $0$ である。

$\det A\ne0$ ならば、両辺を $\det A$ で割れば $A\cdot\frac{1}{\det A}\operatorname{adj}(A)=\frac{1}{\det A}\operatorname{adj}(A)\cdot A=I_n$ となる。
:::

例えば $A=\begin{pmatrix}a&b\\c&d\end{pmatrix}$ では $C_{11}=d$, $C_{12}=-c$, $C_{21}=-b$, $C_{22}=a$ なので

$$
\operatorname{adj}(A)=\begin{pmatrix}d&-b\\-c&a\end{pmatrix},\qquad
A^{-1}=\frac{1}{ad-bc}\begin{pmatrix}d&-b\\-c&a\end{pmatrix}\quad(ad-bc\ne0).
$$

:eqref[eq-linalg-adjugate] 自体は $ad-bc=0$ でも成り立つ。逆行列の公式は、逆行列の各成分が $A$ の成分の有理式（分母は $\det A$）であることを示しており、理論的に重要である。一方、数値計算には :ref[exm-linalg-inverse] の掃き出し法のほうが効率的である。

### 多項式や整数を成分とする行列式

次章では、成分に変数 $t$ を含む行列 $tI_n-A$ の行列式を扱う。また整数を成分とする行列の逆行列が整数行列になる条件も興味深い。そこで、行列式の理論がどこまで一般の「数の体系」で成り立つかを確認しておく。

加法と乗法をもち、体の公理（:ref[def-linalg-field]）のうち「$0$ でない元の乗法逆元の存在」を除いたすべてを満たす集合 $R$ を**可換環**という（[環と加群](/math/rings-modules/)）。整数全体 $\Z$ や、多項式全体 $K[t]$ は可換環であるが体ではない。可換環 $R$ の元を成分とする $n$ 次正方行列全体を $M_n(R)$ と書き、$A\in M_n(R)$ の行列式を Leibniz の公式 :eqref[eq-linalg-leibniz] で定義する。

:::proposition[可換環上の行列式]{#prop-linalg-det-ring}
$R$ を可換環とする。$M_n(R)$ の行列式について、次が成り立つ。

1. 行列式は列についても行についても多重線形かつ交代的であり、$\det I_n=1$ である。
2. $\det(AB)=\det A\det B$, $\det\transp{A}=\det A$ であり、:ref[prop-linalg-det-operations] の 1–3 と :ref[prop-linalg-det-block-triangular] が成り立つ。
3. 余因子展開 :eqref[eq-linalg-cofactor-expansion] と :eqref[eq-linalg-adjugate] が成り立つ。
4. $A\in M_n(R)$ に対し、$AB=BA=I_n$ を満たす $B\in M_n(R)$ が存在するための必要十分条件は、$\det A$ が $R$ で乗法逆元をもつことである。
5. $\varphi:R\to R'$ が和・積・$1$ を保つ写像（環準同型）ならば、成分ごとに $\varphi$ を施した行列 $\varphi(A)$ について $\det\varphi(A)=\varphi(\det A)$ である。
:::

:::proof
1–3：:ref[thm-linalg-determinant-unique] の存在部分の証明、:ref[lem-linalg-alternating-permute] の 1, 2、:ref[cor-linalg-alternating-forms]、:ref[thm-linalg-det-product]、:ref[prop-linalg-det-transpose]、:ref[prop-linalg-det-operations] の 1–3、:ref[prop-linalg-det-block-triangular]、:ref[thm-linalg-cofactor-expansion]、:ref[thm-linalg-adjugate-identity] の証明を見直すと、用いたのは和と積の結合法則・交換法則・分配法則、$0$ と $1$ の性質、加法の逆元だけであり、$0$ でない元で割る操作も、一次独立性や階数も用いていない（ここで $R^n$ の元の「線形性」は $R$ の元を係数とする一次結合について考える）。よってこれらはそのまま $R$ 上で成り立つ。用いることができないのは、:ref[lem-linalg-alternating-permute] の 3、:ref[cor-linalg-det-invertible]、および以下の Cramer の公式のように、割り算や一次独立性に基づく結果である。

4. $AB=I_n$ ならば $\det A\det B=1$ なので $\det A$ は逆元 $\det B$ をもつ。逆に $\det A$ が逆元 $d$ をもてば、:eqref[eq-linalg-adjugate] より $B=d\operatorname{adj}(A)\in M_n(R)$ が $AB=BA=I_n$ を満たす。

5. Leibniz の公式の右辺は成分の和と積だけからなる（$\sgn(\sigma)=\pm1$ も含め）ので、$\varphi$ と可換である。
:::

例えば $R=K[t]$ で $\varphi$ を $t=c$ の代入 $p\mapsto p(c)$ とすれば、5 は「行列式を計算してから $t=c$ を代入しても、先に代入してから行列式を計算しても同じ」ことを述べている。また 4 により、整数行列 $A$ の逆行列が整数行列であるための必要十分条件は $\det A=\pm1$ である（$\Z$ で逆元をもつのは $\pm1$ だけ）。

## Cramer の公式

### 公式とその証明

:::theorem[Cramer の公式]{#thm-linalg-cramer}
$A=(a_1\ \cdots\ a_n)\in M_n(K)$ を正則行列とし、$b\in K^n$ とする。$A$ の第 $j$ 列を $b$ に置き換えた行列を $A_j(b)$ と書くと、$Ax=b$ のただ 1 つの解 $x=\transp{(x_1,\dots,x_n)}$ は

$$
x_j=\frac{\det A_j(b)}{\det A}\qquad(1\le j\le n) \label{eq-linalg-cramer}
$$

で与えられる。
:::

:::proof
$A$ は正則なので、解 $x=A^{-1}b$ はただ 1 つ存在し、$b=Ax=\sum_{k=1}^n x_ka_k$ である。$A_j(b)$ の第 $j$ 列にこの和を代入し、第 $j$ 列についての線形性で展開すると

$$
\det A_j(b)=\sum_{k=1}^n x_k\det(a_1\ \cdots\ a_{j-1}\ a_k\ a_{j+1}\ \cdots\ a_n)
$$

となる。$k\ne j$ の項では列 $a_k$ が 2 回現れるので $0$ であり、$k=j$ の項は $x_j\det A$ である。よって $\det A_j(b)=x_j\det A$ であり、:ref[cor-linalg-det-invertible] より $\det A\ne0$ なので両辺を割ればよい。
:::

:::example[Cramer の公式による解法]{#exm-linalg-cramer}
$2x+y=5$, $x-y=1$ の係数行列を $A$、右辺を $b$ とすると

$$
\begin{aligned}
\det A&=\det\begin{pmatrix}2&1\\1&-1\end{pmatrix}=-3,\\
\det A_1(b)&=\det\begin{pmatrix}5&1\\1&-1\end{pmatrix}=-6,\\
\det A_2(b)&=\det\begin{pmatrix}2&5\\1&1\end{pmatrix}=-3
\end{aligned}
$$

なので $x=2$, $y=1$ である。どの分子も同じ右辺 $b$ を用いるが、置き換える列は求めたい未知数に応じて変わる。
:::

### Cramer の公式の意味と注意

$b=e_i$ とおくと、$A^{-1}$ の第 $i$ 列 $x=A^{-1}e_i$ の第 $j$ 成分は $\det A_j(e_i)/\det A$ である。$A_j(e_i)$ を第 $j$ 列に沿って展開すると $\det A_j(e_i)=C_{ij}$ なので、$(A^{-1})_{ji}=C_{ij}/\det A$ となる。これは :ref[thm-linalg-adjugate-identity] の逆行列の公式と一致する。

$K=\R$ のとき、Cramer の公式には幾何学的な意味がある。$\det A_j(b)$ は、平行体の第 $j$ 辺 $a_j$ を $b=\sum_kx_ka_k$ に置き換えたものの符号付き体積であり、$b$ の $a_j$ 以外の成分は他の辺に平行なずらしなので体積を変えない。したがって体積は $x_j$ 倍になる。

Cramer の公式は解を係数の有理式として明示するので、解が係数に連続的に依存することや、係数が整数で $\det A=\pm1$ ならば整数の右辺に対して整数解をもつことなどを示すのに有用である。ただし計算量の点では、$n+1$ 個の $n$ 次行列式が必要なので、掃き出し法に比べて効率が悪い。

:::remark
$\det A=0$ の場合にも、解 $x$ が存在すれば、証明中の等式 $\det A_j(b)=x_j\det A=0$ が成り立つ。したがって、ある $j$ について $\det A_j(b)\ne0$ ならば解は存在しない。しかし逆は成り立たない。すなわち、すべての $j$ について $\det A_j(b)=0$ であっても、解が存在するとは限らない。例えば $A=\begin{pmatrix}1&1\\1&1\end{pmatrix}$, $b=\transp{(1,2)}$ では $\det A=0$ で、$\det A_1(b)=1-2=-1\ne0$ なので解はない。また $A=O_{2,2}$, $b\ne0$ では置き換えた後の行列にも零列が残るので分子はすべて $0$ だが、$Ax=b$ は解をもたない。$\det A=0$ の場合の解の存在の判定には :ref[thm-linalg-rank-criterion] を用いる。
:::

## 演習問題

::::exercise{level=basic}
$\sigma=[2,4,1,3]\in S_4$ の転倒をすべて求め、$\sigma$ と $\sigma^{-1}$ の符号を求めよ。さらに $\sigma$ を互いに素な巡回置換の積と、隣接互換の積で表せ。
:::solution
転倒は $(1,3)$（$2>1$）、$(2,3)$（$4>1$）、$(2,4)$（$4>3$）の 3 個なので $\sgn(\sigma)=-1$ であり、:ref[thm-linalg-sign-product] より $\sgn(\sigma^{-1})=-1$ である。$1\mapsto2\mapsto4\mapsto3\mapsto1$ なので $\sigma=(1\ 2\ 4\ 3)$ であり、長さ 4 の巡回置換の符号 $(-1)^3=-1$ とも一致する。$s_i=(i\ \ i+1)$ とし、:ref[lem-linalg-adjacent-inversions] の証明のとおり減少箇所を順に入れ替えると

$$
[2,4,1,3]\ \xrightarrow{s_2}\ [2,1,4,3]\
\xrightarrow{s_1}\ [1,2,4,3]\ \xrightarrow{s_3}\ [1,2,3,4]
$$

となる。各矢印は右から掛ける操作なので $\sigma s_2s_1s_3=\id$ であり、$\sigma=s_3s_1s_2$ を得る。実際 $s_3s_1s_2$ は $1\mapsto1\mapsto2\mapsto2$, $2\mapsto3\mapsto3\mapsto4$, $3\mapsto2\mapsto1\mapsto1$, $4\mapsto4\mapsto4\mapsto3$ と確かめられる。
:::
::::

::::exercise{level=basic}
実数 $t$ に対して次の行列の行列式を因数分解し、正則でない $t$ をすべて求めよ。

$$
A(t)=\begin{pmatrix}t&1&0\\1&t&0\\2&3&t-2\end{pmatrix}
$$
:::solution
第 3 列に沿って展開する。第 3 列で $0$ でない成分は $(3,3)$ 成分だけで、その余因子の符号は $(-1)^{3+3}=+1$ である。したがって

$$
\det A(t)=(t-2)\det\begin{pmatrix}t&1\\1&t\end{pmatrix}=(t-2)(t^2-1)=(t-2)(t-1)(t+1).
$$

:ref[cor-linalg-det-invertible] より、正則でないのは $t=-1,1,2$ のときに限る。
:::
::::

::::exercise{level=basic}
$A=\begin{pmatrix}1&2&0\\0&1&1\\1&0&1\end{pmatrix}$ の余因子行列と逆行列を求めよ。
:::solution
各余因子を計算すると、例えば $C_{12}=-\det\begin{pmatrix}0&1\\1&1\end{pmatrix}=1$, $C_{21}=-\det\begin{pmatrix}2&0\\0&1\end{pmatrix}=-2$, $C_{23}=-\det\begin{pmatrix}1&2\\1&0\end{pmatrix}=2$ であり、全体では

$$
(C_{ij})=\begin{pmatrix}1&1&-1\\-2&1&2\\2&-1&1\end{pmatrix},
\qquad
\operatorname{adj}(A)=\begin{pmatrix}1&-2&2\\1&1&-1\\-1&2&1\end{pmatrix}
$$

となる。第 1 行に沿って展開すると $\det A=1\cdot1+2\cdot1+0\cdot(-1)=3$ である。よって $A^{-1}=\tfrac13\operatorname{adj}(A)$ であり、:ref[exm-linalg-inverse] の掃き出し法による結果と一致する。
:::
::::

::::exercise{level=basic}
$2x+y-z=1$, $x+y+z=4$, $x-y+2z=2$ を Cramer の公式で解け。
:::solution
係数行列を $A$、右辺を $b=\transp{(1,4,2)}$ とする。それぞれ第 1 行に沿って展開すると

$$
\begin{aligned}
\det A &=2(2+1)-(2-1)+(-1)(-1-1)=7,\\
\det A_1(b)&=1(2+1)-(8-2)+(-1)(-4-2)=3,\\
\det A_2(b)&=2(8-2)-(2-1)+(-1)(2-4)=13,\\
\det A_3(b)&=2(2+4)-(2-4)+1(-1-1)=12
\end{aligned}
$$

である。したがって解は $(x,y,z)=\bigl(\tfrac37,\tfrac{13}7,\tfrac{12}7\bigr)$ であり、元の方程式に代入して確かめられる（例えば $2\cdot\tfrac37+\tfrac{13}7-\tfrac{12}7=1$）。
:::
::::

::::exercise{level=basic}
次の行列式を、行基本変形で上三角行列に帰着させて計算せよ。

$$
\det\begin{pmatrix}1&0&2&1\\2&1&0&3\\0&3&1&2\\1&1&1&1\end{pmatrix}
$$
:::solution
第 2 行から第 1 行の 2 倍を、第 4 行から第 1 行を引くと（行列式は変わらない）

$$
\begin{pmatrix}1&0&2&1\\0&1&-4&1\\0&3&1&2\\0&1&-1&0\end{pmatrix}.
$$

さらに第 3 行から第 2 行の 3 倍を、第 4 行から第 2 行を引くと

$$
\begin{pmatrix}1&0&2&1\\0&1&-4&1\\0&0&13&-1\\0&0&3&-1\end{pmatrix}.
$$

:ref[prop-linalg-det-block-triangular] より、求める行列式は $1\cdot1\cdot\det\begin{pmatrix}13&-1\\3&-1\end{pmatrix}=-13+3=-10$ である。
:::
::::

::::exercise{level=basic}
$1+1\ne0$ を満たす体 $K$ 上の $n$ 次交代行列 $A$（$\transp{A}=-A$）について、$n$ が奇数ならば $\det A=0$ であることを示せ。$n$ が偶数のときはどうか。
:::solution
:ref[prop-linalg-det-transpose] と :ref[prop-linalg-det-operations] の 2 より $\det A=\det\transp{A}=\det(-A)=(-1)^n\det A$ である。$n$ が奇数ならば $\det A=-\det A$、すなわち $2\det A=0$ であり、$2\ne0$ より $\det A=0$ である。$n$ が偶数のときはこの議論は何も与えず、例えば $\begin{pmatrix}0&1\\-1&0\end{pmatrix}$ の行列式は $1\ne0$ である。
:::
::::

::::exercise{level=advanced}
$\sigma\in S_n$ について、有理数の範囲で

$$
\sgn(\sigma)=\prod_{1\le i<j\le n}\frac{\sigma(j)-\sigma(i)}{j-i}
$$

が成り立つことを示せ。さらにこの式から、:ref[lem-linalg-adjacent-inversions] を用いずに積法則 $\sgn(\sigma\tau)=\sgn(\sigma)\sgn(\tau)$ を導け。
:::solution
分子 $\prod_{i<j}(\sigma(j)-\sigma(i))$ を考える。$\sigma$ は全単射なので、$i<j$ の組 $(i,j)$ と、値の 2 元集合 $\{k,l\}$（$k<l$）は $\{\sigma(i),\sigma(j)\}=\{k,l\}$ によって 1 対 1 に対応する。対応する因子 $\sigma(j)-\sigma(i)$ は、$(i,j)$ が転倒でなければ $l-k$、転倒ならば $-(l-k)$ である。よって分子は $(-1)^{N(\sigma)}\prod_{k<l}(l-k)$ であり、分母 $\prod_{i<j}(j-i)$ と比べて主張を得る。

積法則：$f(a,b)\defeq\frac{\sigma(b)-\sigma(a)}{b-a}$ は $a,b$ について対称（$f(a,b)=f(b,a)$）である。

$$
\prod_{i<j}\frac{\sigma\tau(j)-\sigma\tau(i)}{j-i}
=\prod_{i<j}f(\tau(i),\tau(j))\cdot\prod_{i<j}\frac{\tau(j)-\tau(i)}{j-i}.
$$

第 1 の積では、$i<j$ のとき $\{\tau(i),\tau(j)\}$ はすべての 2 元集合をちょうど 1 回ずつ動き、$f$ は対称なので、積は $\prod_{a<b}f(a,b)=\sgn(\sigma)$ に等しい。第 2 の積は $\sgn(\tau)$ である。
:::
::::

::::exercise{level=advanced}
$u,v\in K^n$ を列ベクトルとする。

1. $\det(I_n+u\transp{v})=1+\transp{v}u$ を示せ。
2. $A\in M_n(K)$ が正則ならば $\det(A+u\transp{v})=\det A\,(1+\transp{v}A^{-1}u)$ を示せ。
:::hint
1 では $I_n+u\transp{v}$ の第 $j$ 列を $e_j+v_ju$ と書いて、列について多重線形に展開する。
:::
:::solution
1. 第 $j$ 列は $e_j+v_ju$ である。各列について線形に展開すると、各列から $e_j$ か $v_ju$ のどちらかを選んだ $2^n$ 項の和になる。すべての列から $e_j$ を選んだ項は $\det I_n=1$ である。2 つ以上の列から $v_ju$ を選んだ項は、スカラーをくくり出すと等しい列 $u$ を 2 つもつので $0$ である。第 $j$ 列だけで $v_ju$ を選んだ項は $v_j\det(e_1\ \cdots\ u\ \cdots\ e_n)$ であり、$u=\sum_iu_ie_i$ と展開すると $i\ne j$ の項は $e_i$ が重複して $0$ なので、$v_ju_j$ である。以上を加えて $\det(I_n+u\transp{v})=1+\sum_jv_ju_j=1+\transp{v}u$ を得る（$1\times1$ 行列 $\transp{v}u$ をその成分と同一視した）。

2. $A+u\transp{v}=A(I_n+(A^{-1}u)\transp{v})$ なので、積公式と 1 より $\det(A+u\transp{v})=\det A\,(1+\transp{v}A^{-1}u)$ である。
:::
::::

::::exercise{level=advanced}
$A\in M_k(K)$ を正則行列、$B\in M_{k,l}(K)$, $C\in M_{l,k}(K)$, $D\in M_l(K)$ とする。

1. $\det\begin{pmatrix}A&B\\C&D\end{pmatrix}=\det A\,\det(D-CA^{-1}B)$ を示せ（$D-CA^{-1}B$ を **Schur 補行列**という）。
2. $k=l$ で $AC=CA$ ならば、$\det\begin{pmatrix}A&B\\C&D\end{pmatrix}=\det(AD-CB)$ を示せ。
:::solution
1. ブロックの積を計算すると

$$
\begin{pmatrix}I_k&O\\-CA^{-1}&I_l\end{pmatrix}\begin{pmatrix}A&B\\C&D\end{pmatrix}=\begin{pmatrix}A&B\\O&D-CA^{-1}B\end{pmatrix}
$$

である。左の因子は対角成分がすべて $1$ の下三角行列なので行列式は $1$ である（:ref[prop-linalg-det-block-triangular]）。積公式と :ref[prop-linalg-det-block-triangular] より、求める行列式は $\det A\det(D-CA^{-1}B)$ である。

2. 1 と積公式より、行列式は $\det(A(D-CA^{-1}B))=\det(AD-ACA^{-1}B)$ である。$AC=CA$ より $ACA^{-1}=CAA^{-1}=C$ なので、これは $\det(AD-CB)$ に等しい。
:::
::::

::::exercise{level=advanced}
$n\ge2$, $A\in M_n(K)$ とする。

1. $\det\operatorname{adj}(A)=(\det A)^{n-1}$ を示せ。
2. $\operatorname{adj}(A)$ の階数は、$\rank A=n$ ならば $n$、$\rank A=n-1$ ならば $1$、$\rank A\le n-2$ ならば $0$ であることを示せ。
:::solution
1. $A$ が正則ならば、:eqref[eq-linalg-adjugate] の行列式をとると $\det A\det\operatorname{adj}(A)=\det((\det A)I_n)=(\det A)^n$ なので、$\det A$ で割って主張を得る。$A$ が正則でないとき、$\det A=0$ なので $\det\operatorname{adj}(A)=0$ を示せばよい。$\operatorname{adj}(A)$ が正則であると仮定すると、$A\operatorname{adj}(A)=O$ に右から $\operatorname{adj}(A)^{-1}$ を掛けて $A=O$ となるが、$n\ge2$ では $O$ の余因子はすべて $0$ なので $\operatorname{adj}(O)=O$ となり、正則性に反する。

2. $\rank A=n$ ならば、1 より $\det\operatorname{adj}(A)\ne0$ なので階数は $n$ である。$\rank A\le n-2$ ならば、:ref[thm-linalg-rank-minor] より $(n-1)$ 次小行列式はすべて $0$ なので $\operatorname{adj}(A)=O$ である。$\rank A=n-1$ ならば、:ref[thm-linalg-rank-minor] より $0$ でない $(n-1)$ 次小行列式があるので $\operatorname{adj}(A)\ne O$、すなわち階数は $1$ 以上である。一方、$A\operatorname{adj}(A)=(\det A)I_n=O$ より $\operatorname{adj}(A)$ の各列は $\Ker L_A$ に属し、:ref[cor-linalg-solution-space-dim] より $\dim\Ker L_A=1$ なので、$\operatorname{adj}(A)$ の列空間は 1 次元以下である。よって階数は $1$ である。
:::
::::

::::exercise{level=advanced}
（**Cauchy–Binet の公式**）$m\le n$ とし、$A\in M_{m,n}(K)$, $B\in M_{n,m}(K)$ とする。$\{1,\dots,n\}$ の $m$ 元部分集合 $S=\{s_1<\cdots<s_m\}$ に対し、$A$ の第 $s_1,\dots,s_m$ 列からなる $m$ 次正方行列を $A_{:,S}$、$B$ の第 $s_1,\dots,s_m$ 行からなる $m$ 次正方行列を $B_{S,:}$ と書く。

1. $\det(AB)=\sum_S\det A_{:,S}\det B_{S,:}$（和は $m$ 元部分集合 $S$ 全体にわたる）を示せ。
2. $K=\R$ のとき $\det(A\transp{A})\ge0$ を示し、$A=\begin{pmatrix}1&2&0\\0&1&1\end{pmatrix}$ について両辺を計算して公式を確かめよ。
3. $m>n$ ならば $\det(AB)=0$ であることを示せ。
:::solution
1. $A$ の列を $a_1,\dots,a_n$ とすると、$AB$ の第 $k$ 列は $\sum_{j=1}^nb_{jk}a_j$ である。各列について多重線形に展開すると

$$
\det(AB)=\sum_{j_1,\dots,j_m=1}^nb_{j_1,1}\cdots b_{j_m,m}\det(a_{j_1}\ \cdots\ a_{j_m})
$$

となる。$j_1,\dots,j_m$ に重複があれば、行列式は等しい列をもつので $0$ である。重複がなければ、$S=\{j_1,\dots,j_m\}$ を小さい順に $s_1<\cdots<s_m$ と並べ、$j_k=s_{\pi(k)}$ となる $\pi\in S_m$ がただ 1 つ定まる。:ref[lem-linalg-alternating-permute] の 2 より $\det(a_{s_{\pi(1)}}\ \cdots\ a_{s_{\pi(m)}})=\sgn(\pi)\det A_{:,S}$ である。逆に $S$ と $\pi$ の組から $(j_1,\dots,j_m)$ が決まるので

$$
\det(AB)=\sum_S\det A_{:,S}\sum_{\pi\in S_m}\sgn(\pi)\prod_{k=1}^mb_{s_{\pi(k)},k}
$$

である。$B_{S,:}$ の $(i,k)$ 成分は $b_{s_i,k}$ なので、内側の和は Leibniz の公式による $\det B_{S,:}$ である。

2. $B=\transp{A}$ とすると $B_{S,:}=\transp{(A_{:,S})}$ なので、1 と :ref[prop-linalg-det-transpose] より $\det(A\transp{A})=\sum_S(\det A_{:,S})^2\ge0$ である。例の $A$ では、$S=\{1,2\},\{1,3\},\{2,3\}$ に対する小行列式は $1,1,2$ なので右辺は $1+1+4=6$。一方 $A\transp{A}=\begin{pmatrix}5&2\\2&2\end{pmatrix}$ で、行列式は $10-4=6$ である。

3. 1 の証明の展開は $m>n$ でも成り立つが、$m$ 個の添字 $j_1,\dots,j_m\in\{1,\dots,n\}$ には必ず重複があるので、すべての項が $0$ である（$\rank(AB)\le n<m$ と :ref[thm-linalg-rank-minor] からも従う）。
:::
::::

::::exercise{level=advanced}
$A\in M_n(\Z)$（整数を成分とする行列）とする。

1. $A^{-1}$ が存在して整数行列となるための必要十分条件は $\det A=\pm1$ であることを示せ。
2. $\det A=\pm1$ ならば、任意の $b\in\Z^n$ に対して $Ax=b$ は整数解をもつことを示せ。$A=\begin{pmatrix}2&1\\1&1\end{pmatrix}$, $b=\transp{(3,5)}$ の整数解を求めよ。
:::solution
1. :ref[prop-linalg-det-ring] の 4 を $R=\Z$ に適用する。$\Z$ で乗法逆元をもつ元は $\pm1$ だけである（$ab=1$ ならば $\abs{a}\abs{b}=1$ で $\abs{a}=1$）。

2. 1 より $A^{-1}=\pm\operatorname{adj}(A)$ は整数行列なので、$x=A^{-1}b\in\Z^n$ である。例では $\det A=2-1=1$, $A^{-1}=\operatorname{adj}(A)=\begin{pmatrix}1&-1\\-1&2\end{pmatrix}$ なので $x=A^{-1}b=\transp{(-2,7)}$ であり、$2\cdot(-2)+7=3$, $-2+7=5$ と確かめられる。
:::
::::

## 確認テスト

::::quiz{answer="-1"}
転倒数が $5$ の置換の符号はいくつか。

:::explanation
符号は $(-1)^{N(\sigma)}$ で定義されるので、転倒数 $5$（奇数）ならば $(-1)^5=-1$ である。
:::
::::

::::quiz
次の操作のうち、すべての正方行列について行列式を変えないものをすべて選べ。ただし行や列の番号は異なるものとする。

- [x] 1 つの列に他の列の $3$ 倍を加える
- [ ] 2 つの行を入れ替える
- [x] 行列を転置する
- [ ] 行列全体を $2$ 倍する

:::explanation
列の加算は :ref[prop-linalg-det-operations] の 3、転置は :ref[prop-linalg-det-transpose] により行列式を変えない。行の入れ替えでは符号が反転し（例えば $I_2$ に施すと $-1$）、行列全体を $2$ 倍すると $n$ 次正方行列では $2^n$ 倍になる。
:::
::::

::::quiz{answer="-6"}
$\det A=2$, $\det B=-3$ のとき、$\det(AB)$ を求めよ。$A,B$ は同じサイズの正方行列とする。

:::explanation
積公式 :eqref[eq-linalg-det-product] から $\det(AB)=2\cdot(-3)=-6$ である。
:::
::::

::::quiz{answer="40"}
$A$ を $3$ 次正方行列とし、$\det A=5$ とする。$\det(2A)$ を求めよ。

:::explanation
$2A$ はすべての列（3 つ）を $2$ 倍したものなので、:ref[prop-linalg-det-operations] の 2 より $\det(2A)=2^3\det A=8\cdot5=40$ である。$2\det A=10$ とする誤りに注意する。
:::
::::

::::quiz{answer="12"}
$S_4$ に含まれる偶置換の個数はいくつか。

:::explanation
$\abs{S_4}=4!=24$ であり、:ref[cor-linalg-cycle-sign] の 3 より偶置換と奇置換は同数なので、偶置換は $24/2=12$ 個である。
:::
::::

::::quiz
$A\in M_n(K)$ が $\det A=0$ を満たすとき、必ず成り立つものをすべて選べ。

- [ ] $A^{-1}=\operatorname{adj}(A)/\det A$
- [x] $A\operatorname{adj}(A)=O$
- [ ] すべての $b$ に対して $Ax=b$ は解をもつ
- [x] $Ax=0$ は自明でない解をもつ

:::explanation
:eqref[eq-linalg-adjugate] は $\det A$ で割らずに成り立つので、$\det A=0$ ならば $A\operatorname{adj}(A)=O$ である。$\det A=0$ のとき $A$ は正則でない（:ref[cor-linalg-det-invertible]）ので、逆行列は存在せず、$\Ker L_A\ne\{0\}$ であり、次元定理より $L_A$ は全射でないので、像の外の $b$ に対して $Ax=b$ は解をもたない。
:::
::::

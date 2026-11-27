---
title: 準同型と同型定理
description: 群準同型の基本性質と構成法、核と像、自己同型群と特性部分群を調べ、第一・第二・第三同型定理と対応定理を証明して剰余群を既知の群と同定する。
notation:
  - symbol: '\Hom(G,H)'
    meaning: G から H への群準同型全体
  - symbol: '\Ker f,\ \Img f'
    meaning: 準同型 f の核と像
  - symbol: '\Aut(G)'
    meaning: G の自己同型群
  - symbol: 'c_g'
    meaning: g による内部自己同型 x ↦ gxg^{-1}
  - symbol: '\operatorname{Inn}(G)'
    meaning: G の内部自己同型群
  - symbol: '\operatorname{Out}(G)'
    meaning: G の外部自己同型群 Aut(G)/Inn(G)
  - symbol: '\pi_N'
    meaning: 正規部分群 N による自然射影 G → G/N
  - symbol: 'S^1'
    meaning: 絶対値 1 の複素数全体のなす乗法群
---

## この章の目標

- 準同型の基本性質を証明し、生成元の行き先と関係式から準同型を構成・数え上げる。
- 核と像を計算し、核が正規部分群であること、ファイバーが核の剰余類であることを用いて単射性と像の大きさを判定する。
- 自己同型群・内部自己同型群・特性部分群を調べ、$\Aut(C_n)\cong(\Z/n\Z)^\times$ などを導く。
- 第一・第二・第三同型定理を証明し、具体的な剰余群を既知の群と同定する。
- 対応定理によって、剰余群の部分群・正規部分群・指数を元の群で計算する。

前章 [剰余類と Lagrange の定理](/math/group-theory/cosets/) では、正規部分群 $N$ による剰余群 $G/N$ を構成した。しかし $G/N$ は剰余類の集合であり、それがどのような群なのか（例えば $\GL_n(K)/\SL_n(K)$ が何であるか）は、定義を眺めるだけでは分からない。本章では、群の間の「積を保つ写像」である準同型を導入する。線形代数で線形写像の核と像を調べると次元定理が得られたように（:ref[thm-linalg-rank-nullity]）、準同型の核と像を調べると「剰余群 $G/\Ker f$ は像と同型である」という第一同型定理が得られる。これによって、剰余群を既知の群として具体的に記述できるようになる。

本章を通じて、群の単位元を $e$（必要なら $e_G$）と書き、置換の積 $\sigma\tau$ は $\tau$ を先に施す合成 $\sigma\circ\tau$ とする。二面体群 $D_n$（$n\ge3$）の元は標準形 $r^is^j$（$0\le i<n$, $0\le j<2$）で表し、$r^n=s^2=e$, $srs=r^{-1}$ である（:ref[prop-grp-dihedral-form]）。

## 準同型写像

整数を $n$ で割った余りだけを記録しても、和の余りは正しく計算できる。この記録法では情報の一部が失われる（$1$ と $n+1$ は区別されない）。一方、正三角形の対称操作を頂点の置換で記録すると、情報を失わずに合成を計算できる。二つの記録法には「記録してから計算しても、計算してから記録しても結果が同じ」という共通点がある。この性質を取り出したものが準同型である。

### 定義と基本性質

:::definition[準同型と同型]{#def-grp-hom}
$G,H$ を群とする。写像 $f:G\to H$ がすべての $a,b\in G$ について

$$
f(ab)=f(a)f(b)
$$

を満たすとき、$f$ を :term[準同型]{en="group homomorphism" yomi="じゅんどうけい"}（群準同型）という。左辺の積は $G$ の積、右辺の積は $H$ の積である。$G$ から $H$ への準同型全体を $\Hom(G,H)$ と書く。

1. 単射な準同型を**単射準同型**（埋め込み）、全射な準同型を**全射準同型**という。
2. 全単射な準同型は、:ref[def-grp-isomorphism] の意味の同型写像にほかならない（条件がまったく同じである）。
3. $G=H$ のとき、準同型 $G\to G$ を**自己準同型**、同型 $G\to G$ を $G$ の :term[自己同型]{en="automorphism" yomi="じこどうけい"} という。
:::

加法で書かれた群（加法群）の間では条件は $f(a+b)=f(a)+f(b)$ となり、乗法群から加法群への写像なら $f(ab)=f(a)+f(b)$ となる。どの演算で書くかは表記の問題にすぎない。

:::proposition[準同型の計算規則]{#prop-grp-hom-rules}
$f:G\to H$, $k:H\to L$ を準同型とする。

1. $f(e_G)=e_H$ であり、任意の $a\in G$ と $m\in\Z$ について $f(a^{-1})=f(a)^{-1}$, $f(a^m)=f(a)^m$ が成り立つ。
2. 合成 $k\circ f:G\to L$ は準同型である。$f$ が全単射ならば、逆写像 $f^{-1}:H\to G$ も準同型である。
3. $S\subseteq G$ が $G$ を生成するとき、$S$ の各元で値が一致する二つの準同型 $f,f':G\to H$ は等しい。
4. $a\in G$ が有限位数ならば $f(a)$ も有限位数で、$\operatorname{ord}f(a)$ は $\operatorname{ord}a$ を割り切る。$f$ が単射ならば $\operatorname{ord}f(a)=\operatorname{ord}a$ である。
:::

:::proof
1. $f(e)=f(ee)=f(e)f(e)$ の両辺に左から $f(e)^{-1}$ を掛けると $e_H=f(e)$ を得る（:ref[prop-grp-cancellation]）。次に $f(a)f(a^{-1})=f(aa^{-1})=f(e)=e_H$, $f(a^{-1})f(a)=f(a^{-1}a)=e_H$ なので、逆元の一意性から $f(a^{-1})=f(a)^{-1}$ である。$m\ge0$ のときは $m$ に関する帰納法で示す。$m=0$ は $f(e)=e_H$ そのものであり、$f(a^m)=f(a)^m$ ならば $f(a^{m+1})=f(a^ma)=f(a^m)f(a)=f(a)^mf(a)=f(a)^{m+1}$ である。$m<0$ のときは $a^m=(a^{-m})^{-1}$ であるから、すでに示したことより $f(a^m)=f(a^{-m})^{-1}=(f(a)^{-m})^{-1}=f(a)^m$ となる。

2. $(k\circ f)(ab)=k(f(a)f(b))=k(f(a))\,k(f(b))$ である。$f$ が全単射のとき、$u,v\in H$ に対し $a=f^{-1}(u)$, $b=f^{-1}(v)$ とおくと $f(ab)=f(a)f(b)=uv$ なので、$f^{-1}(uv)=ab=f^{-1}(u)f^{-1}(v)$ である。

3. :ref[prop-grp-generated-words] により、$G$ の任意の元は $g=s_1^{\eps_1}\cdots s_m^{\eps_m}$（$m\ge0$, $s_i\in S$, $\eps_i=\pm1$）と書ける。$m$ に関する帰納法と 1 により $f(g)=f(s_1)^{\eps_1}\cdots f(s_m)^{\eps_m}$ であり、$f'(g)$ も同じ式で与えられる。仮定より $f(s_i)=f'(s_i)$ なので $f(g)=f'(g)$ である。

4. $n=\operatorname{ord}a$ とすると $f(a)^n=f(a^n)=f(e)=e_H$ である。:ref[thm-grp-cyclic-order] により、ある元の $m$ 乗が単位元になることと、その元の位数が $m$ を割り切ることは同値である。よって $\operatorname{ord}f(a)\mid n$ である。$f$ が単射なら、$f(a)^m=e_H$ は $f(a^m)=f(e)$、すなわち $a^m=e$ と同値である。冪が単位元になる指数の集合が $a$ と $f(a)$ で一致するので、位数も一致する。
:::

:::remark
$f(e)=e$ は準同型であるための必要条件だが十分条件ではない。例えば $S_3$ 上の写像 $x\mapsto x^2$ は $e\mapsto e$ を満たすが、後の例で見るように準同型ではない。また、1 により準同型が逆元を保つことは自動的に従うので、定義に「逆元を保つ」を課す必要はない。
:::

:::example[準同型の例]{#exm-grp-hom-examples}
1. **自明な準同型・恒等写像・包含**：すべての元を $e_H$ に送る写像 $G\to H$ を**自明な準同型**という。恒等写像 $\id_G$ は自己同型である。$H\le G$ のとき包含写像 $\iota:H\to G$ は単射準同型である。
2. **法 $n$ への還元**：$n\ge1$ に対し、$\Z\to\Z/n\Z$, $a\mapsto\bar a$ は準同型である。$\overline{a+b}=\bar a+\bar b$ は $\Z/n\Z$ の加法の定義そのものである（:ref[exm-grp-integers-mod]）。この写像は全射である。
3. **冪写像**：$A$ を可換群とし $m\in\Z$ とする。$x\mapsto x^m$ は $A$ の自己準同型である。実際、可換な二元 $x,y$ について $(xy)^m=x^my^m$ が成り立つ。$m\ge0$ では帰納法により、$y^mx=xy^m$（$y$ と $x$ が可換なので $y$ の冪とも可換）を用いて $(xy)^{m+1}=(xy)^mxy=x^my^mxy=x^{m+1}y^{m+1}$ となり、$m<0$ では $(xy)^m=((xy)^{-m})^{-1}=(x^{-m}y^{-m})^{-1}=y^mx^m=x^my^m$ となる。加法で書けば $x\mapsto mx$ である。非可換群ではこれは一般に成り立たない。$S_3$ で $x\mapsto x^2$ とすると、$(1\,2)(2\,3)=(1\,2\,3)$ の像は $(1\,2\,3)^2=(1\,3\,2)$ だが、$(1\,2)^2(2\,3)^2=e$ である。
4. **行列式**：体 $K$ について $\det:\GL_n(K)\to K^\times$ は準同型である（:ref[thm-linalg-det-product]）。$\det\diag(t,1,\ldots,1)=t$ なので全射である。
5. **符号**：$\sgn:S_n\to\{\pm1\}$ は準同型である（:ref[thm-linalg-sign-product]）。ここで $\{\pm1\}$ は乗法について位数 $2$ の群である。$n\ge2$ なら $\sgn(1\,2)=-1$ なので全射である。
6. **指数関数**：$\exp:(\R,+)\to(\R_{>0},\times)$ は $e^{x+y}=e^xe^y$ により準同型である。$\exp$ は狭義単調増加で値域が $\R_{>0}$ 全体であり（[解析学](/math/analysis/)）、逆写像 $\log$ をもつので同型である。足し算の群と掛け算の群が同型であることは、対数表による掛け算の計算の原理にほかならない。複素数でも $\exp:(\C,+)\to\C^\times$, $z\mapsto e^z$ は $e^{z+w}=e^ze^w$ により準同型である。$w\in\C^\times$ を極形式で $w=\rho e^{i\theta}$（$\rho>0$）と書けば $w=e^{\log\rho+i\theta}$ なので全射であるが、$e^{2\pi i}=1=e^0$ なので単射ではない。
7. **円周への巻き付け**：$S^1=\{z\in\C:\abs z=1\}$ は $\abs{zw}=\abs z\abs w$, $\abs{z^{-1}}=\abs z^{-1}$ により $\C^\times$ の部分群である。$\R\to S^1$, $t\mapsto e^{2\pi it}$ は準同型である。また絶対値 $\C^\times\to\R_{>0}$, $z\mapsto\abs z$ も準同型である。
8. **二面体群の頂点置換**：$D_n$ の各元は正 $n$ 角形の頂点を置換する。頂点を $\Z/n\Z$ で番号づけ、$g\in D_n$ に頂点上の置換 $\rho(g)$ を対応させる。写像の合成を頂点に制限したものは制限の合成なので、$\rho:D_n\to S_{\Z/n\Z}$ は準同型である。$r(x)=x+1$, $s(x)=-x$ より $\rho(r^is^j)(x)=i+(-1)^jx$ である。$\rho(r^is^j)=\rho(r^ks^l)$ とすると、$x=0$ を代入して $i\equiv k$、$x=1$ を代入して $(-1)^j\equiv(-1)^l\pmod n$ を得る。$n\ge3$ では $1\not\equiv-1\pmod n$ なので $j=l$ である。よって $\rho$ は単射準同型である。頂点を $1,\ldots,n$ と番号づけ直せば、$D_n$ は $S_n$ に埋め込まれる。
9. **行列による実現**：$D_n$ の元を平面の $2$ 次直交行列で表す準同型 $D_n\to\GL_2(\R)$ もある。これは後の :ref[prop-grp-hom-dihedral] を用いて構成する。
10. **共役**：$g\in G$ に対し $c_g(x)=gxg^{-1}$ は $G$ の自己同型である。これは後の :ref[def-grp-inner-automorphism] で詳しく扱う。
:::

準同型でない写像の例も挙げておく。$\Z\to\Z$, $x\mapsto x+1$ は $0$ を $1$ に送るので準同型でない。行列式 $\det:M_n(K)\to K$ は $\det(AB)=\det A\det B$ を満たすが、$M_n(K)$ は乗法について群でない（$0$ に逆元がない）ので、群準同型という概念の外にある。

### 同型で保たれる性質

:ref[prop-grp-isomorphism-basic] で見たように、同型は位数、可換性、各元の位数を保つ。従って、ある群論的な性質を一方の群だけがもてば、二つの群は同型でない。逆に同型を示すには、具体的に全単射な準同型を作る必要がある。

- $(\R,+)\cong(\R_{>0},\times)$ は $\exp$ が与える（:ref[exm-grp-hom-examples] の 6）。
- $(\R,+)$ と $(\R^\times,\times)$ は同型でない。$\R^\times$ には位数 $2$ の元 $-1$ があるが、$(\R,+)$ で $2x=0$ となるのは $x=0$ だけである。
- $(\Q,+)$ と $(\Q_{>0},\times)$ は同型でない。同型 $\phi:\Q\to\Q_{>0}$ があったとすると、$\phi(x)=2$ となる $x$ をとれば $\phi(x/2)^2=\phi(x)=2$ となり、$\sqrt2$ が有理数となって矛盾する。$\Q$ では任意の元が「$2$ で割れる」が、$\Q_{>0}$ では $2$ の平方根がとれない。
- $C_4$ と $C_2\times C_2$（:ref[exm-grp-direct-product-basic]）は位数 $4$ の可換群だが、前者には位数 $4$ の元があり後者にはない。
- $Q_8$ と $D_4$ は位数 $8$ の非可換群だが、位数 $2$ の元は $Q_8$ では $-1$ だけ、$D_4$ では $r^2$ と四つの鏡映の五つである。

### 巡回群からの準同型

巡回群からの準同型は、生成元の行き先だけで決まる。行き先として許される元を正確に特定する。

:::proposition[巡回群からの準同型]{#prop-grp-hom-from-cyclic}
$H$ を群とする。

1. 任意の $h\in H$ に対し、$f(1)=h$ を満たす準同型 $f:\Z\to H$ がただ一つ存在し、$f(m)=h^m$ で与えられる。従って $f\mapsto f(1)$ は全単射 $\Hom(\Z,H)\to H$ である。
2. $n\ge1$ とする。$h\in H$ に対し、$f(\bar1)=h$ を満たす準同型 $f:C_n=\Z/n\Z\to H$ が存在するための必要十分条件は $h^n=e$ である。このとき $f$ は一意で $f(\bar m)=h^m$ であり、像は $\langle h\rangle$、$f$ が単射である条件は $\operatorname{ord}h=n$ である。従って $f\mapsto f(\bar1)$ は $\Hom(C_n,H)$ から $\{h\in H:h^n=e\}$ への全単射である。
3. $m,n\ge1$ に対し、$\Hom(C_m,C_n)$ の元の個数は $\gcd(m,n)$ である。
:::

:::proof
1. 指数法則 $h^{m+k}=h^mh^k$（:ref[prop-grp-power-laws]）により $m\mapsto h^m$ は $(\Z,+)$ からの準同型で、$1\mapsto h$ である。一意性は :ref[prop-grp-hom-rules] の 3（$\Z$ は $1$ で生成される）から従う。

2. 必要性：$f$ が準同型なら $h^n=f(\bar1)^n=f(n\cdot\bar1)=f(\bar0)=e$ である（加法群 $C_n$ での $n\cdot\bar1$ は $\bar1$ の $n$ 個の和）。十分性：$h^n=e$ とし、$f(\bar m)=h^m$ と定めたい。そのためには右辺が代表元 $m$ のとり方によらないことを示す必要がある。$\bar m=\bar k$ なら $m-k=qn$ と書け、$h^m=h^kh^{qn}=h^k(h^n)^q=h^k$ である。こうして定義された $f$ は、$f(\bar m+\bar k)=f(\overline{m+k})=h^{m+k}=h^mh^k$ により準同型である。一意性は 1 と同様である。像は $\{h^m:m\in\Z\}=\langle h\rangle$ である。$f(\bar m)=f(\bar k)$ は $h^{m-k}=e$、すなわち $\operatorname{ord}h\mid m-k$ と同値である（:ref[thm-grp-cyclic-order]）。$f$ が単射であることは「$\operatorname{ord}h\mid m-k$ ならば $n\mid m-k$」と同値であり、$m-k=\operatorname{ord}h$ とおけばこれは $n\mid\operatorname{ord}h$ を意味する。一方 $h^n=e$ より $\operatorname{ord}h\mid n$ なので、単射性は $\operatorname{ord}h=n$ と同値である。

3. 2 により、数えるべきは $mh\equiv0\pmod n$ を満たす $h\in\{0,1,\ldots,n-1\}$ の個数である。$d=\gcd(m,n)$, $m=dm'$, $n=dn'$ とおくと $\gcd(m',n')=1$ である。$n\mid mh$ は $n'\mid m'h$ と同値である。$um'+vn'=1$ となる整数 $u,v$ をとると $h=um'h+vn'h$ なので、$n'\mid m'h$ は $n'\mid h$ と同値である。$0\le h<n$ の範囲にある $n'$ の倍数は $0,n',\ldots,(d-1)n'$ の $d$ 個である。
:::

:::example[準同型を数える]
1. $C_6\to C_4$ の準同型は $\gcd(6,4)=2$ 個で、$6h\equiv0\pmod 4$ を満たす $h=\bar0,\bar2$ に対応する。非自明な方は $\bar m\mapsto\overline{2m}$ である。
2. $C_{12}\to C_{18}$ の準同型は $\gcd(12,18)=6$ 個で、$\bar1$ の行き先は $\bar0,\bar3,\bar6,\bar9,\overline{12},\overline{15}$ である。
3. $\Hom(C_2,S_3)$ の元は $h^2=e$ となる $h\in S_3$、すなわち $e,(1\,2),(1\,3),(2\,3)$ の四つに対応する。$\Hom(C_3,S_3)$ の元は $e,(1\,2\,3),(1\,3\,2)$ の三つに対応する。
4. 「$C_4\to C_6$, $\bar m\mapsto\bar m$」は写像として定義されない。$C_4$ では $\bar0=\bar4$ だが、行き先の $\bar0$ と $\bar4$ は $C_6$ で異なる。:ref[prop-grp-hom-from-cyclic] の第 2 項にある条件 $4h=\bar0$ は $h=\bar1$ で破れている。
:::

### 生成元と関係式からの準同型

巡回群でない群でも、生成元の行き先が「生成元の間の関係式」を満たせば準同型が得られることがある。二面体群で正確に述べる。一般の群の表示による扱いは [自由群と群の表示](/math/group-theory/free-groups/) で行う。

:::proposition[二面体群からの準同型]{#prop-grp-hom-dihedral}
$n\ge3$ とし、$H$ を群、$A,B\in H$ とする。$f(r)=A$, $f(s)=B$ を満たす準同型 $f:D_n\to H$ が存在するための必要十分条件は

$$
A^n=e,\qquad B^2=e,\qquad BAB^{-1}=A^{-1} \label{eq-grp-dihedral-relations}
$$

が成り立つことである。このとき $f$ は一意で、$f(r^is^j)=A^iB^j$ である。従って $f\mapsto(f(r),f(s))$ は、$\Hom(D_n,H)$ から :eqref[eq-grp-dihedral-relations] を満たす組 $(A,B)$ 全体への全単射である。
:::

:::proof
必要性：$D_n$ で $r^n=e$, $s^2=e$, $srs^{-1}=srs=r^{-1}$ が成り立つので、両辺に $f$ を施して :ref[prop-grp-hom-rules] の 1 を用いれば :eqref[eq-grp-dihedral-relations] を得る。

十分性：:eqref[eq-grp-dihedral-relations] を仮定する。まず任意の整数 $k$ について

$$
BA^k=A^{-k}B
$$

が成り立つ。実際、写像 $x\mapsto BxB^{-1}$ は $BxyB^{-1}=(BxB^{-1})(ByB^{-1})$ により $H$ の自己準同型なので、:ref[prop-grp-hom-rules] の 1 から $BA^kB^{-1}=(BAB^{-1})^k=(A^{-1})^k=A^{-k}$ であり、右から $B$ を掛ければよい。これを $j\in\{0,1\}$ について繰り返せば $B^jA^k=A^{(-1)^jk}B^j$ である。

:ref[prop-grp-dihedral-form] により $D_n$ の各元は $r^is^j$（$0\le i<n$, $0\le j<2$）と一意に書けるので、$f(r^is^j)=A^iB^j$ によって写像 $f:D_n\to H$ が定まる。さらに任意の整数 $i,j$ に対しても $f(r^is^j)=A^iB^j$ が成り立つ。実際、$i=qn+i'$, $j=2t+j'$（$0\le i'<n$, $0\le j'<2$）と書くと、$r^n=s^2=e$ より $r^is^j=r^{i'}s^{j'}$ であり、$A^n=B^2=e$ より $A^iB^j=A^{i'}B^{j'}$ である。

準同型であることを確かめる。:ref[prop-grp-dihedral-form] の積の公式 $(r^is^j)(r^ks^l)=r^{i+(-1)^jk}s^{j+l}$ と、上で示した二つの事実から

$$
f\bigl((r^is^j)(r^ks^l)\bigr)=A^{i+(-1)^jk}B^{j+l}=A^i\bigl(B^jA^k\bigr)B^l=A^iB^jA^kB^l=f(r^is^j)\,f(r^ks^l)
$$

である。一意性は、$D_n$ が $r,s$ で生成されることと :ref[prop-grp-hom-rules] の 3 から従う。
:::

この命題の要点は、関係式 :eqref[eq-grp-dihedral-relations] を確かめるだけで、$2n\times2n$ 通りの積をすべて調べる必要がなくなることである。関係式は必要条件として自然に現れるが、それが十分でもあることは、標準形の積の公式が関係式だけから導かれていることによる。

:::example[二面体群の行列表示]
$\theta=2\pi/n$ とし、

$$
R=\begin{pmatrix}\cos\theta&-\sin\theta\\ \sin\theta&\cos\theta\end{pmatrix},\qquad S=\begin{pmatrix}1&0\\0&-1\end{pmatrix}
$$

とおく。$R^k$ は角 $k\theta$ の回転行列なので $R^n=I$ であり、$S^2=I$ である。また

$$
SRS^{-1}=SRS=\begin{pmatrix}\cos\theta&\sin\theta\\ -\sin\theta&\cos\theta\end{pmatrix}=R^{-1}
$$

である。よって :ref[prop-grp-hom-dihedral] により、$r\mapsto R$, $s\mapsto S$ を満たす準同型 $f:D_n\to\GL_2(\R)$ がただ一つ存在し、$f(r^is^j)=R^iS^j$ である。$f(r^is^j)=I$ とすると、行列式をとって $(-1)^j=1$ より $j=0$、続いて $R^i=I$ より $n\mid i$ を得る。従って $f$ は単射であり、$D_n$ は $\GL_2(\R)$ の位数 $2n$ の部分群 $\{R^iS^j\}$ と同型である。
:::

:::example[二面体群から位数 2 の群へ]
$C_2=\{\bar0,\bar1\}$ を加法で書く。$D_n\to C_2$ の準同型は、:eqref[eq-grp-dihedral-relations] を加法で書いた $nA=\bar0$, $2B=\bar0$, $B+A-B=-A$ を満たす $(A,B)$ に対応する。第三式は $2A=\bar0$ となり、$C_2$ では常に成り立つ。第二式も常に成り立つ。$n$ が奇数なら第一式から $A=\bar0$ であり、準同型は $B$ の選び方の二つである。$n$ が偶数なら $A,B$ を独立に選べるので四つある。

$D_4$ の非自明な三つの準同型を $(A,B)=(\bar0,\bar1),(\bar1,\bar0),(\bar1,\bar1)$ と書くと、$r^is^j$ の行き先はそれぞれ $\bar j$, $\bar i$, $\overline{i+j}$ である。つまり「鏡映か」「回転指数が奇数か」「その二条件の一方だけが成り立つか」を記録している。
:::

:::example[関係式が破れる場合]
$D_n\to C_n$ で $r\mapsto\bar1$, $s\mapsto\bar0$ とする準同型を作ろうとすると、第三の関係式は $\bar1=-\bar1$、すなわち $2\equiv0\pmod n$ を要求する。$n\ge3$ ではこれは成り立たないので、そのような準同型は存在しない。可換群への準同型では $f(r)=f(srs^{-1})=f(r^{-1})=f(r)^{-1}$ となるため、回転 $r$ の像の位数は常に $2$ 以下である。
:::

:::example[二面体群の間の準同型]{#exm-grp-dihedral-quotient-map}
$m\ge3$ が $n$ を割り切るとし、$D_m$ の生成元を $r',s'$ と書く。$A=r'$, $B=s'$ とおくと $r'^n=(r'^m)^{n/m}=e$, $s'^2=e$, $s'r's'^{-1}=r'^{-1}$ なので、$r\mapsto r'$, $s\mapsto s'$ を満たす準同型 $D_n\to D_m$ が存在する。$r^is^j\mapsto r'^is'^j$ が単位元になるのは $j=0$ かつ $m\mid i$ のときなので、単位元に写る元は $\langle r^m\rangle$ の $n/m$ 個である。例えば $D_6\to D_3$ では $e$ と $r^3$ が単位元に写る。正六角形の向かい合う頂点の組（三組）を正三角形の頂点とみなすと、この準同型は組の置換を記録するものである。
:::

### 準同型全体の演算

値域 $B$ が可換群ならば、二つの準同型 $f,g:G\to B$ の値ごとの積 $(fg)(x)=f(x)g(x)$ はまた準同型である。実際

$$
(fg)(xy)=f(x)f(y)g(x)g(y)=f(x)g(x)f(y)g(y)=(fg)(x)\,(fg)(y)
$$

であり、二番目の等号で $B$ の可換性を用いた。自明な準同型が単位元、$x\mapsto f(x)^{-1}$ が $f$ の逆元となり（$f(x)^{-1}$ は可換群の冪写像 $y\mapsto y^{-1}$ と $f$ の合成）、$\Hom(G,B)$ は Abel 群をなす。値域が非可換ならこの入れ替えはできない。例えば非可換群 $G$ の恒等写像を自分自身と値ごとに掛けると平方写像 $x\mapsto x^2$ になるが、これが準同型なら $xyxy=x^2y^2$ の両端を消去して $yx=xy$ が全元で成り立ってしまう。

:ref[prop-grp-hom-from-cyclic] の全単射 $\Hom(C_m,C_n)\to\{h:mh=\bar0\}$ はこの演算を保つ。従って $\Hom(C_m,C_n)$ は $C_n$ の部分群 $\{h:mh=\bar0\}$ と同型である。後者は巡回群 $C_n$ の部分群なので巡回群であり（:ref[thm-grp-cyclic-subgroups]）、位数は $\gcd(m,n)$ である。よって $\Hom(C_m,C_n)\cong C_{\gcd(m,n)}$ である。

### 自己同型群と内部自己同型

群 $G$ の自己同型は「元の名前の付け替えで積の表を保つもの」であり、群の対称性を表す。自己同型全体もまた群をなす。

:::definition[自己同型群と内部自己同型]{#def-grp-inner-automorphism}
群 $G$ の自己同型全体を $\Aut(G)$ と書き、合成を積として $G$ の :term[自己同型群]{en="automorphism group" yomi="じこどうけいぐん"} という。$g\in G$ に対し

$$
c_g:G\to G,\qquad c_g(x)=gxg^{-1}
$$

を $g$ による :term[内部自己同型]{en="inner automorphism" yomi="ないぶじこどうけい"}（共役）といい、内部自己同型全体を $\operatorname{Inn}(G)$ と書く。
:::

:::proposition[自己同型群の基本性質]{#prop-grp-aut-basic}
1. $\Aut(G)$ は写像の合成について群をなす。
2. 各 $c_g$ は $G$ の自己同型であり、$c_gc_h=c_{gh}$, $c_e=\id_G$ が成り立つ。$\operatorname{Inn}(G)$ は $\Aut(G)$ の部分群である。また $c_g=c_h$ と $h^{-1}g\in Z(G)$ は同値である。
3. 任意の $\phi\in\Aut(G)$ と $g\in G$ について $\phi\circ c_g\circ\phi^{-1}=c_{\phi(g)}$ である。特に $\operatorname{Inn}(G)\trianglelefteq\Aut(G)$ である。
:::

:::proof
1. :ref[prop-grp-hom-rules] の 2 により、自己同型の合成と逆写像は自己同型である。恒等写像は自己同型であり、写像の合成は結合的である。

2. $c_g(xy)=gxyg^{-1}=(gxg^{-1})(gyg^{-1})=c_g(x)c_g(y)$ なので $c_g$ は準同型である。$c_g(c_h(x))=ghxh^{-1}g^{-1}=(gh)x(gh)^{-1}=c_{gh}(x)$ であり、$c_e=\id_G$ は明らかである。従って $c_gc_{g^{-1}}=c_e=c_{g^{-1}}c_g$ となり、$c_g$ は全単射で逆写像 $c_{g^{-1}}$ をもつ。$\operatorname{Inn}(G)$ は $\id_G$ を含み、$c_gc_h^{-1}=c_gc_{h^{-1}}=c_{gh^{-1}}\in\operatorname{Inn}(G)$ なので部分群である（:ref[prop-grp-subgroup-test]）。最後に、$c_g=c_h$ は両辺に左から $c_{h^{-1}}$ を合成して $c_{h^{-1}g}=\id_G$ と同値であり、これは全 $x$ について $h^{-1}gx(h^{-1}g)^{-1}=x$、すなわち $h^{-1}g$ が全元と可換であることと同値である。これが $h^{-1}g\in Z(G)$ の定義である（:ref[def-grp-center]）。

3. $x\in G$ について $\phi(c_g(\phi^{-1}(x)))=\phi(g)\,\phi(\phi^{-1}(x))\,\phi(g^{-1})=\phi(g)\,x\,\phi(g)^{-1}=c_{\phi(g)}(x)$ である。従って $\phi\operatorname{Inn}(G)\phi^{-1}\subseteq\operatorname{Inn}(G)$ がすべての $\phi$ について成り立ち、:ref[prop-grp-normal-criteria] により $\operatorname{Inn}(G)$ は正規部分群である。
:::

剰余群 $\operatorname{Out}(G)=\Aut(G)/\operatorname{Inn}(G)$ を $G$ の :term[外部自己同型群]{en="outer automorphism group" yomi="がいぶじこどうけいぐん"} という。内部自己同型でない自己同型を**外部自己同型**ということもあるが、$\operatorname{Out}(G)$ の元は自己同型そのものではなく内部自己同型を法とした類である。$G$ が可換なら $c_g=\id_G$ なので $\operatorname{Inn}(G)=\{\id_G\}$, $\operatorname{Out}(G)\cong\Aut(G)$ である。

:::example[内部自己同型と外部自己同型]
1. :ref[prop-grp-cycle-conjugation] により $c_{(1\,2\,3)}((1\,2))=(2\,3)$ である。内部自己同型は置換の「文字の名前の付け替え」を表す。
2. $\GL_n(K)$ の写像 $\tau(A)=(\transp A)^{-1}$ は自己同型である。実際 $\tau(AB)=(\transp B\,\transp A)^{-1}=(\transp A)^{-1}(\transp B)^{-1}=\tau(A)\tau(B)$ であり、$\tau\circ\tau=\id$ なので全単射である。$K=\R$ のとき $\tau$ は内部自己同型でない。内部自己同型は行列式を保つ（$\det(gAg^{-1})=\det A$）が、$A=\diag(2,1,\ldots,1)$ について $\det\tau(A)=1/2\ne2=\det A$ だからである。
3. $(\Z,+)$ の自己同型 $x\mapsto-x$ は外部自己同型である。可換群では恒等写像以外の自己同型はすべて内部自己同型でない。
:::

巡回群の自己同型群は完全に決定できる。

:::theorem[巡回群の自己同型群]{#thm-grp-aut-cyclic}
1. $\Aut(\Z)=\{\id,-\id\}$ であり、$\Aut(\Z)\cong\{\pm1\}$ である。
2. $n\ge1$ とする。$\bar u\in(\Z/n\Z)^\times$ に $\mu_u(x)=ux$ を対応させる写像 $(\Z/n\Z)^\times\to\Aut(C_n)$ は同型である。特に $\abs{\Aut(C_n)}=\varphi(n)$ である。
:::

:::proof
1. :ref[prop-grp-hom-from-cyclic] の 1 により、$\Z$ の自己準同型は $h=f(1)$ を用いて $f(m)=mh$ と書ける。その像は $h\Z$ であり、全射となるのは $h\Z=\Z$、すなわち $h=\pm1$ のときに限る。$h=\pm1$ なら $f=\pm\id$ は全単射である。$\pm\id$ の合成は符号の積に対応するので、$\Aut(\Z)\cong\{\pm1\}$ である。

2. :ref[prop-grp-hom-from-cyclic] の 2 で $H=C_n$ とすると、条件 $nh=\bar0$ は常に満たされるので、$C_n$ の自己準同型は $h=\bar u=f(\bar1)$ を任意に選んで $f(\bar m)=\overline{mu}$、すなわち $f=\mu_u$ と書ける。$\mu_u$ の像は $\langle\bar u\rangle$ であり、$C_n$ は有限集合なので、$\mu_u$ が全単射であることは全射であること、すなわち $\langle\bar u\rangle=C_n$ と同値である。:ref[thm-grp-cyclic-order] により $\bar u$ の位数は $n/\gcd(n,u)$ なので、これは $\gcd(u,n)=1$、すなわち $\bar u\in(\Z/n\Z)^\times$ と同値である（:ref[exm-grp-integers-mod]）。

従って $\bar u\mapsto\mu_u$ は $(\Z/n\Z)^\times$ から $\Aut(C_n)$ への全射である。$\mu_u(\bar1)=\bar u$ なので単射でもある。最後に $\mu_u\circ\mu_v(x)=uvx=\mu_{uv}(x)$ なので準同型である。
:::

:::example[巡回群の自己同型群の具体例]
1. $\Aut(C_8)\cong(\Z/8\Z)^\times=\{\bar1,\bar3,\bar5,\bar7\}$ である。$3^2=9$, $5^2=25$, $7^2=49$ はすべて法 $8$ で $1$ なので、非単位元の位数はすべて $2$ である。従って $\Aut(C_8)$ は巡回群でなく、:ref[prop-grp-order-four] により $C_2\times C_2$ と同型である。
2. $\Aut(C_9)\cong(\Z/9\Z)^\times$ は位数 $6$ である。$2$ の冪は法 $9$ で $2,4,8,7,5,1$ と巡るので、$\bar2$ が生成元であり $\Aut(C_9)\cong C_6$ である。対応する自己同型 $\mu_2$ は $x\mapsto2x$ である。
3. $p$ が素数なら $(\Z/p\Z)^\times=\F_p^\times$ は巡回群である（:ref[thm-grp-field-units-cyclic]）。従って $\Aut(C_p)\cong C_{p-1}$ である。例えば $\Aut(C_7)$ は $x\mapsto3x$ で生成される（$3$ の冪は法 $7$ で $3,2,6,4,5,1$）。
:::

:::example[小さい群の自己同型群]{#exm-grp-aut-small}
1. **$\Aut(V_4)\cong S_3$**：Klein の四元群 $V_4=\{e,a,b,c\}$（:ref[exm-grp-klein-four]）では、非単位元はすべて位数 $2$ で、異なる二つの非単位元の積は残りの一つである（$ab=c$ など）。自己同型は $e$ を固定し位数を保つので、$\{a,b,c\}$ を置換する。この制限によって単射準同型 $\Aut(V_4)\to S_{\{a,b,c\}}$ が得られる（:ref[prop-grp-hom-rules] の 3 より、$\{a,b,c\}$ 上の値で自己同型は決まる）。逆に $\{a,b,c\}$ の任意の置換 $\pi$ を $e\mapsto e$ で延長した写像 $\phi$ は準同型である。実際、$x$ または $y$ が $e$ なら明らかであり、$x=y\ne e$ なら $\phi(x^2)=\phi(e)=e=\phi(x)^2$、$x\ne y$ がともに非単位元なら $xy$ は第三の元 $z$ で、$\phi(x)\phi(y)$ は $\phi(x),\phi(y)$ 以外の非単位元、すなわち $\phi(z)$ である（$\pi$ は全単射）。従って $\Aut(V_4)\cong S_3$ であり、$V_4$ は可換なので $\operatorname{Out}(V_4)\cong S_3$ である。
2. **$\Aut(S_3)=\operatorname{Inn}(S_3)\cong S_3$**：自己同型は位数を保つので、$S_3$ の位数 $2$ の元全体、すなわち三つの互換の集合 $T=\{(1\,2),(1\,3),(2\,3)\}$ を置換する。$T$ は $S_3$ を生成するので、制限によって単射準同型 $\Aut(S_3)\to S_T$ が得られ、$\abs{\Aut(S_3)}\le6$ である。一方 $Z(S_3)=\{e\}$ である。実際、$z\in Z(S_3)$ なら :ref[prop-grp-cycle-conjugation] より $(z(1)\,z(2))=z(1\,2)z^{-1}=(1\,2)$ なので $z$ は $\{1,2\}$ を保ち、同様に $\{1,3\}$ も保つので $z(1)=1$、従って $z=e$ である。:ref[prop-grp-aut-basic] の 2 により $g\mapsto c_g$ は単射なので $\abs{\operatorname{Inn}(S_3)}=6$ である。よって $\Aut(S_3)=\operatorname{Inn}(S_3)$ であり、$g\mapsto c_g$ は同型 $S_3\to\Aut(S_3)$ を与える。
:::

## 核と像

二つの元が同じ値に写る条件を、単位元に写る元だけで判定する。これは線形写像の核で単射性を判定した議論（:ref[prop-linalg-kernel-injective]）の乗法版である。

:::definition[核と像]{#def-grp-kernel}
準同型 $f:G\to H$ に対して

$$
\Ker f=\{g\in G:f(g)=e_H\},\qquad \Img f=f(G)=\{f(g):g\in G\}
$$

をそれぞれ $f$ の :term[核]{en="kernel" yomi="かく"}、:term[像]{en="image" yomi="ぞう"} という。$y\in H$ の逆像 $f^{-1}(\{y\})$ を $y$ 上の**ファイバー**ともいう。
:::

:::proposition[部分群の像と逆像]{#prop-grp-image-preimage}
$f:G\to H$ を準同型とし、$A\le G$, $B\le H$ とする。

1. $f(A)\le H$, $f^{-1}(B)\le G$ である。特に $\Img f\le H$, $\Ker f\le G$ である。
2. $B\trianglelefteq H$ ならば $f^{-1}(B)\trianglelefteq G$ である。特に $\Ker f\trianglelefteq G$ である。
3. $A\trianglelefteq G$ ならば $f(A)\trianglelefteq f(G)$ である。特に $f$ が全射なら $f(A)\trianglelefteq H$ である。
:::

:::proof
1. $f(e_G)=e_H$ なので $f(A)$ は空でない。$f(a),f(b)\in f(A)$ に対して $f(a)f(b)^{-1}=f(ab^{-1})\in f(A)$ である。$f^{-1}(B)$ は $e_G$ を含み、$x,y\in f^{-1}(B)$ なら $f(xy^{-1})=f(x)f(y)^{-1}\in B$ である。いずれも :ref[prop-grp-subgroup-test] により部分群である。$A=G$, $B=\{e_H\}$ とすれば後半を得る。

2. $g\in G$, $x\in f^{-1}(B)$ とすると $f(gxg^{-1})=f(g)f(x)f(g)^{-1}\in B$ である（$B$ の正規性）。従って $g\,f^{-1}(B)\,g^{-1}\subseteq f^{-1}(B)$ が全 $g$ で成り立ち、:ref[prop-grp-normal-criteria] により $f^{-1}(B)$ は正規である。$B=\{e_H\}$ は $H$ の正規部分群なので $\Ker f$ も正規である。

3. $f(G)$ の元は $f(g)$ の形である。$a\in A$ について $f(g)f(a)f(g)^{-1}=f(gag^{-1})\in f(A)$ である（$A$ の正規性）。従って $f(A)$ は $f(G)$ の中で共役で閉じ、正規である。
:::

:::proposition[ファイバー]{#prop-grp-fiber}
$f:G\to H$ を準同型とし、$K=\Ker f$ とおく。$a,b\in G$ について $f(a)=f(b)$ と $a^{-1}b\in K$ は同値である。特に $f(a)$ 上のファイバーは左剰余類 $aK$ であり（右剰余類 $Ka$ でもある）、$f$ が単射であるための必要十分条件は $K=\{e\}$ である。
:::

:::proof
$f(a)^{-1}f(b)=f(a^{-1}b)$ なので、$f(a)=f(b)$ は $f(a^{-1}b)=e_H$、すなわち $a^{-1}b\in K$ と同値である。$a^{-1}b\in K$ と $b\in aK$ は同値なので、$f(a)$ 上のファイバーは $aK$ である。同様に $f(b)f(a)^{-1}=f(ba^{-1})$ からファイバーは $Ka$ でもある（$K$ が正規なので $aK=Ka$ と整合する）。核が自明なら各ファイバー $aK=\{a\}$ は一点なので $f$ は単射である。逆に単射なら、$e$ と同じ値 $e_H$ をとる元は $e$ のみなので $K=\{e\}$ である。
:::

従って準同型の単射性は、すべての二元の組を比べなくても、「$f(g)=e$ なら $g=e$」を確かめるだけで判定できる。また、準同型 $f$ は $G$ を核の剰余類に分割し、各剰余類を $H$ の一点に潰す写像である。この描像が次節の第一同型定理の内容である。

:::corollary[核と像の位数]{#cor-grp-kernel-image-order}
$G$ を有限群、$f:G\to H$ を準同型とする。このとき

$$
\abs G=\abs{\Ker f}\cdot\abs{\Img f} \label{eq-grp-kernel-image-order}
$$

である。特に $\abs{\Img f}$ は $\abs G$ を割り切り、$H$ も有限なら $\abs{\Img f}$ は $\gcd(\abs G,\abs H)$ を割り切る。
:::

:::proof
$G$ は $\Img f$ の各点上のファイバーの非交和である。:ref[prop-grp-fiber] により各ファイバーは $aK$（$K=\Ker f$）の形であり、$x\mapsto ax$ は $K$ から $aK$ への全単射（逆写像は $y\mapsto a^{-1}y$）なので、各ファイバーはちょうど $\abs K$ 個の元をもつ。ファイバーは $\abs{\Img f}$ 個あるので :eqref[eq-grp-kernel-image-order] を得る。$\Img f\le H$ なので、$H$ が有限なら Lagrange の定理（:ref[thm-grp-lagrange]）により $\abs{\Img f}$ は $\abs H$ も割り切る。
:::

:eqref[eq-grp-kernel-image-order] は線形代数の次元定理「$\dim V=\dim\Ker f+\dim\Img f$」の群版である。次元の和が位数の積に置き換わっている。

:::example[核と像の計算]{#exm-grp-kernel-examples}
1. **行列式**：$\det:\GL_2(\F_3)\to\F_3^\times$ は全射準同型（:ref[exm-grp-hom-examples] の 4）で、核は $\SL_2(\F_3)$ である。:ref[prop-grp-gl-order] より $\abs{\GL_2(\F_3)}=(9-1)(9-3)=48$ なので、:eqref[eq-grp-kernel-image-order] から $\abs{\SL_2(\F_3)}=48/2=24$ である。一般に $q$ 元からなる有限体 $K$ 上で $\abs{\SL_n(K)}=\abs{\GL_n(K)}/(q-1)$ であり、例えば $\abs{\SL_2(\F_5)}=(25-1)(25-5)/4=120$ である。
2. **符号**：$\sgn:S_n\to\{\pm1\}$（$n\ge2$）は全射で、核は偶置換全体 $A_n$（:ref[exm-grp-alternating-subgroup]）である。従って $\abs{A_n}=n!/2$ であり、$A_n\trianglelefteq S_n$ である。これは :ref[cor-linalg-cycle-sign] の数え上げと一致する。
3. **複素指数関数**：$\exp:\C\to\C^\times$ の核は $2\pi i\Z$ である。実際 $z=x+iy$ について $e^z=e^x(\cos y+i\sin y)=1$ は、絶対値を比べて $e^x=1$、すなわち $x=0$ であり、さらに $\cos y=1$, $\sin y=0$、すなわち $y\in2\pi\Z$ と同値である。同様に $\R\to S^1$, $t\mapsto e^{2\pi it}$ の核は $\Z$ である。
4. **巡回群の倍写像**：$C_n$ の自己準同型 $x\mapsto kx$ の像は $\langle\bar k\rangle$ で位数は $n/\gcd(n,k)$ である（:ref[thm-grp-cyclic-order]）。従って :eqref[eq-grp-kernel-image-order] より核の位数は $\gcd(n,k)$ である。例えば $C_{12}$ で $x\mapsto8x$ の像は $\{\bar0,\bar4,\bar8\}$、核は $\{\bar0,\bar3,\bar6,\bar9\}$ である。
5. **平方写像**：$\C^\times$ の自己準同型 $z\mapsto z^2$（:ref[exm-grp-hom-examples] の 3）の核は $\{\pm1\}$ であり、$\rho e^{i\theta}=(\sqrt\rho\,e^{i\theta/2})^2$ より全射である。$\R^\times$ 上では核は同じく $\{\pm1\}$ だが、像は $\R_{>0}$ である。同じ式で定義された準同型でも、像は定義域によって大きく変わる。
6. **二面体群**：$D_n\to C_2$, $r^is^j\mapsto\bar j$ の核は回転部分群 $\langle r\rangle$ である。行列表示 $D_n\to\GL_2(\R)$ の核は自明である。
:::

:::example[正規な像の範囲]
$C_2\to S_3$ で生成元を $(1\,2)$ に送る単射準同型を考える。$C_2$ は自分自身の中で正規だが、その像 $\langle(1\,2)\rangle$ は $S_3$ の中で正規でない（$(1\,3)(1\,2)(1\,3)^{-1}=(2\,3)$）。:ref[prop-grp-image-preimage] の 3 で像の正規性が保証されるのは $f(G)$ の中であり、終域 $H$ 全体の中ではない。一方、$H$ の正規部分群を逆像に取ると、$f$ の全射性によらず $G$ で正規である。
:::

:::example[核から分かる準同型の制約]
1. 素数位数 $p$ の群 $G$ からの準同型 $f$ は自明か単射である。$\Ker f$ は $G$ の部分群なので、Lagrange の定理により位数は $1$ か $p$ だからである。
2. $G,H$ が有限群で $\gcd(\abs G,\abs H)=1$ ならば $G\to H$ の準同型は自明なものに限る。:ref[cor-grp-kernel-image-order] より像の位数が $1$ だからである。例えば $C_5\to S_3$ の準同型は自明なものしかない。
3. $S_3\to C_3$ の準同型は自明なものに限る。像の位数は $\gcd(6,3)=3$ の約数だが、像が $C_3$ 全体なら核は位数 $2$ の正規部分群である。$S_3$ の位数 $2$ の部分群は互換で生成される三つで、どれも正規でないので、像は自明である。
:::

### 特性部分群

正規部分群は内部自己同型で不変な部分群である。すべての自己同型で不変という、より強い条件を考える。

:::definition[特性部分群]{#def-grp-characteristic}
部分群 $L\le G$ が、すべての $\phi\in\Aut(G)$ について $\phi(L)=L$ を満たすとき、$L$ を $G$ の :term[特性部分群]{en="characteristic subgroup" yomi="とくせいぶぶんぐん"} という。
:::

$\phi(L)\subseteq L$ がすべての $\phi$ で成り立てば、$\phi^{-1}$ にも適用して $\phi^{-1}(L)\subseteq L$、すなわち $L\subseteq\phi(L)$ となるので等号が従う。

:::proposition[特性部分群の性質]{#prop-grp-characteristic-basic}
1. 特性部分群は正規部分群である。
2. $L$ が $K$ の特性部分群で $K\trianglelefteq G$ ならば、$L\trianglelefteq G$ である。
3. $L$ が $K$ の特性部分群で、$K$ が $G$ の特性部分群ならば、$L$ は $G$ の特性部分群である。
4. 中心 $Z(G)$ は $G$ の特性部分群である。また、有限群 $G$ が位数 $m$ の部分群をただ一つしかもたないとき、その部分群は特性部分群である。
:::

:::proof
1. 内部自己同型 $c_g$ も自己同型なので $gLg^{-1}=c_g(L)=L$ である。

2. $g\in G$ とする。$K$ は正規なので $c_g(K)=K$ であり、$c_g$ の $K$ への制限 $c_g|_K:K\to K$ は準同型で、$c_{g^{-1}}|_K$ を逆写像にもつ。従って $c_g|_K\in\Aut(K)$ であり、$L$ が $K$ の特性部分群なので $gLg^{-1}=c_g|_K(L)=L$ である。

3. $\phi\in\Aut(G)$ とすると $\phi(K)=K$ なので、2 と同様に $\phi|_K\in\Aut(K)$ であり、$\phi(L)=\phi|_K(L)=L$ である。

4. $z\in Z(G)$, $\phi\in\Aut(G)$ とする。任意の $y\in G$ は $\phi$ の全射性により $y=\phi(x)$ と書け、$\phi(z)y=\phi(zx)=\phi(xz)=y\phi(z)$ である。よって $\phi(z)\in Z(G)$ である。後半は、$\phi(L)$ が $L$ と同じ位数の部分群になる（$\phi$ は単射準同型）ことから従う。
:::

:::example[特性部分群と正規部分群の違い]
1. **巡回群**：有限巡回群の部分群は位数ごとにただ一つなので（:ref[thm-grp-cyclic-subgroups]）、すべて特性部分群である。例えば $C_{12}$ の位数 $3$ の部分群 $\langle\bar4\rangle$ は特性部分群である。$\Z$ の部分群 $d\Z$ も $\Aut(\Z)=\{\pm\id\}$ で保たれる。
2. **$V_4$**：$C_2\times C_2$ の部分群 $C_2\times\{\bar0\}$ は可換性により正規だが、二つの座標を入れ替える自己同型で $\{\bar0\}\times C_2$ に移るので、特性部分群でない。
3. **$Q_8$ の外部自己同型**：:ref[exm-grp-quaternion] の行列表示 $u=\begin{pmatrix}\mathrm i&0\\0&-\mathrm i\end{pmatrix}$, $v=\begin{pmatrix}0&1\\-1&0\end{pmatrix}$, $w=uv$ で $Q_8=\{\pm I,\pm u,\pm v,\pm w\}$ を表す（$\mathrm i=\sqrt{-1}$ であり、$u,v,w$ が四元数の $i,j,k$ に当たる）。$u^2=v^2=-I$, $vu=-uv$ を用いる。$P=u+v=\begin{pmatrix}\mathrm i&1\\-1&-\mathrm i\end{pmatrix}$ とおくと $\det P=2$ で $P$ は正則であり、

$$
Pu=u^2+vu=-I+vu=v(u+v)=vP,\qquad Pv=uv+v^2=uv-I=u(u+v)=uP
$$

なので $PuP^{-1}=v$, $PvP^{-1}=u$、従って $PwP^{-1}=(PuP^{-1})(PvP^{-1})=vu=-w$ である。写像 $\kappa(x)=PxP^{-1}$ は $\GL_2(\C)$ の上で積を保ち、$\pm I\mapsto\pm I$, $\pm u\mapsto\pm v$, $\pm v\mapsto\pm u$, $\pm w\mapsto\mp w$ と $Q_8$ を自分自身の上に全単射に写す。よって $\kappa$ の $Q_8$ への制限は $Q_8$ の自己同型であり、四元数の記号では $i\mapsto j$, $j\mapsto i$, $k\mapsto-k$ である。$\langle i\rangle$ は指数 $2$ なので正規だが（:ref[prop-grp-index-two]）、$\kappa(\langle i\rangle)=\langle j\rangle\ne\langle i\rangle$ なので特性部分群でない。

さらに $\kappa$ は内部自己同型でない。$Q_8$ の元 $g$ による共役は $i$ を $\pm i$ に送るからである。実際 $g=\pm1,\pm i$ なら $gig^{-1}=i$ であり、$g=\pm j$ なら $jij^{-1}=ji(-j)=-(ji)j=kj=-i$、$g=\pm k$ なら $kik^{-1}=ki(-k)=-(ki)k=-jk=-i$ である。従って $\operatorname{Out}(Q_8)$ は自明でない。

4. **推移性の破れ**：正規部分群の正規部分群は正規とは限らない。$D_4$ で $\langle s\rangle\trianglelefteq\langle r^2,s\rangle\trianglelefteq D_4$ だが $\langle s\rangle$ は $D_4$ で正規でない（[剰余類と Lagrange の定理](/math/group-theory/cosets/) の例）。:ref[prop-grp-characteristic-basic] の 2 は、中間の段で「正規」を「特性」に強めれば推移性が回復することを述べている。$\langle s\rangle$ は $\langle r^2,s\rangle\cong V_4$ の特性部分群でない。
:::

特性部分群の概念は、可解群の研究で、導来群や中心などの「群の積だけから定まる部分群」が次々に正規部分群になることを示すのに用いる（[可解群と冪零群](/math/group-theory/solvable-nilpotent/)）。

## 第一同型定理

核の一つの剰余類を一点にまとめれば、準同型が失う情報をちょうど取り除ける。これを定式化するため、まず剰余群への標準的な準同型を用意し、次に別の群への準同型を剰余群経由で表す条件を調べる。

### 自然射影

:::definition[自然射影]{#def-grp-natural-projection}
$N\trianglelefteq G$ とする。$g\in G$ に剰余類 $gN$ を対応させる写像

$$
\pi_N:G\to G/N,\qquad \pi_N(g)=gN
$$

を :term[自然射影]{en="canonical projection" yomi="しぜんしゃえい"}（標準射影）という。
:::

:::proposition[自然射影の性質]{#prop-grp-natural-projection}
$\pi_N$ は全射準同型であり、$\Ker\pi_N=N$ である。従って、部分群 $N\le G$ が正規であることと、$N$ がある準同型の核であることは同値である。
:::

:::proof
剰余群の積の定義（:ref[thm-grp-quotient]）により $\pi_N(gh)=ghN=(gN)(hN)=\pi_N(g)\pi_N(h)$ である。$G/N$ の元はすべて $gN$ の形なので全射である。剰余群の単位元は $N=eN$ であり、$gN=eN$ は $g\in N$ と同値である（:ref[prop-grp-coset-test]）から $\Ker\pi_N=N$ である。後半について、正規部分群は $\pi_N$ の核であり、逆に核は正規である（:ref[prop-grp-image-preimage]）。
:::

正規部分群は「ある準同型で単位元とみなされる元の集まり」として特徴付けられた。これが正規性の本当の意味である。

### 剰余群の普遍性

:::proposition[剰余群の普遍性]{#prop-grp-quotient-universal}
$N\trianglelefteq G$ とし、$f:G\to H$ を準同型とする。準同型 $\bar f:G/N\to H$ で $f=\bar f\circ\pi_N$ を満たすものが存在するための必要十分条件は $N\subseteq\Ker f$ である。このとき $\bar f$ は一意で $\bar f(gN)=f(g)$ により与えられ、

$$
\Img\bar f=\Img f,\qquad \Ker\bar f=\Ker f/N=\{gN:g\in\Ker f\}
$$

が成り立つ。特に $\bar f$ が単射であることと $N=\Ker f$ は同値である。
:::

:::proof
必要性：$f=\bar f\circ\pi_N$ なら、$n\in N$ について $f(n)=\bar f(nN)=\bar f(N)=e_H$ である（$\bar f$ は準同型で $N$ は $G/N$ の単位元）。

十分性：$N\subseteq\Ker f$ と仮定する。$gN=g'N$ なら $g^{-1}g'\in N\subseteq\Ker f$ なので、:ref[prop-grp-fiber] により $f(g)=f(g')$ である。よって $\bar f(gN)=f(g)$ は代表元のとり方によらず定まる。$\bar f((gN)(hN))=\bar f(ghN)=f(gh)=f(g)f(h)=\bar f(gN)\bar f(hN)$ なので準同型であり、定義から $\bar f\circ\pi_N=f$ である。

一意性：$f=\bar f\circ\pi_N$ を満たす準同型は、$\bar f(gN)=\bar f(\pi_N(g))=f(g)$ を満たさなければならず、$G/N$ の元はすべて $gN$ の形なので値が決まる。

像と核：$\Img\bar f=\{f(g):g\in G\}=\Img f$ である。$\bar f(gN)=e_H$ は $f(g)=e_H$ と同値なので、$\Ker\bar f=\{gN:g\in\Ker f\}=\pi_N(\Ker f)$ であり、$N\le\Ker f$ なのでこれは $\Ker f/N$ と書ける。$\bar f$ が単射であることは $\Ker\bar f$ が単位元 $N$ だけからなること、すなわち $\Ker f\subseteq N$ と同値であり、仮定 $N\subseteq\Ker f$ と合わせて $N=\Ker f$ と同値である。
:::

この命題は「$N$ を潰す準同型は $G/N$ を経由する」と言い換えられる。$G/N$ は、$N$ の元を単位元とみなす群のうちで「最も情報を残したもの」である。

### 第一同型定理

:::theorem[第一同型定理]{#thm-grp-first-iso}
任意の群準同型 $f:G\to H$ について

$$
G/\Ker f\cong\Img f
$$

であり、同型は $g\Ker f\mapsto f(g)$ で与えられる。特に $f$ が全射なら $G/\Ker f\cong H$ である。
:::

:::proof
$N=\Ker f$ は正規部分群である（:ref[prop-grp-image-preimage]）。:ref[prop-grp-quotient-universal] を $N=\Ker f$ に適用すると、$\bar f(gN)=f(g)$ で定まる準同型 $\bar f:G/N\to H$ は単射で、像は $\Img f$ である。値域を $\Img f$ に制限すれば全単射な準同型、すなわち同型 $G/\Ker f\to\Img f$ を得る。
:::

まとめると、任意の準同型は

$$
G\xrightarrow{\ \pi\ }G/\Ker f\xrightarrow{\ \bar f\ }\Img f\xrightarrow{\ \iota\ }H
$$

と、全射（自然射影）、同型、単射（包含）の合成に分解される。例えば $G/N$ を群 $H$ と同定するには、核が $N$ である全射準同型 $G\to H$ を見つければよい。これが第一同型定理の典型的な使い方である。

:::example[第一同型定理による剰余群の同定]{#exm-grp-first-iso-examples}
1. **$\Z/n\Z$**：還元写像 $\Z\to C_n$, $a\mapsto\bar a$ は全射で核は $n\Z$ である。従って剰余群 $\Z/n\Z$（:ref[thm-grp-quotient] の意味）は $C_n$ と同型であり、実際に集合としても同じものである。
2. **$\GL_n(K)/\SL_n(K)\cong K^\times$**：$\det$ は全射で核は $\SL_n(K)$ である。$\GL_n(K)$ の二つの行列が同じ剰余類に入ることは、行列式が等しいことと同値である。
3. **$S_n/A_n\cong\{\pm1\}$**（$n\ge2$）：$\sgn$ は全射で核は $A_n$ である。$S_n$ は偶置換と奇置換の二つの剰余類に分かれる。
4. **$\R/\Z\cong S^1$**：$t\mapsto e^{2\pi it}$ は全射（$S^1$ の元は $e^{i\theta}$ と書ける）で、核は $\Z$ である（:ref[exm-grp-kernel-examples] の 3）。実数を整数部分を忘れて記録すると円周上の点になる。同様に $\C/2\pi i\Z\cong\C^\times$ である。
5. **$\C^\times/S^1\cong\R_{>0}$**：$z\mapsto\abs z$ は全射（$\rho>0$ は $\rho$ 自身の像）で核は $S^1$ である。$S^1$ の剰余類は原点を中心とする円である。
6. **$D_n/\langle r\rangle\cong C_2$**：$r^is^j\mapsto\bar j$ は全射で核は $\langle r\rangle$ である。
7. **$D_n/\langle r^m\rangle\cong D_m$**（$m\ge3$, $m\mid n$）：:ref[exm-grp-dihedral-quotient-map] の準同型 $D_n\to D_m$ は生成元 $r',s'$ を像に含むので全射であり、核は $\langle r^m\rangle$ である。
:::

:::example[倍写像の行き先]
$f:C_{12}\to C_{18}$, $\bar x\mapsto\overline{3x}$ は、:ref[prop-grp-hom-from-cyclic] の条件 $12\cdot\bar3=\overline{36}=\bar0$（法 $18$）を満たすので準同型である。$3x\equiv0\pmod{18}$ は $x\equiv0\pmod6$ と同値なので核は $\{\bar0,\bar6\}$、像は法 $18$ の $3$ の倍数全体 $\langle\bar3\rangle$（位数 $6$）である。従って $C_{12}/\langle\bar6\rangle\cong\langle\bar3\rangle\cong C_6$ であり、:eqref[eq-grp-kernel-image-order] は $12=2\cdot6$ となる。像は終域 $C_{18}$ 全体とは位数が異なる。第一同型定理が同定するのは像であり、終域ではない。
:::

### 中心と内部自己同型

:::proposition[内部自己同型群と中心]{#prop-grp-inn-center}
任意の群 $G$ について $G/Z(G)\cong\operatorname{Inn}(G)$ であり、同型は $gZ(G)\mapsto c_g$ で与えられる。
:::

:::proof
:ref[prop-grp-aut-basic] の 2 により、$\Phi:G\to\Aut(G)$, $g\mapsto c_g$ は準同型（$c_{gh}=c_gc_h$）で、像は定義により $\operatorname{Inn}(G)$ である。$c_g=\id_G=c_e$ は $g=e^{-1}g\in Z(G)$ と同値なので、$\Ker\Phi=Z(G)$ である。第一同型定理を適用すればよい。
:::

:::example[内部自己同型の個数]
1. $Q_8$ の中心は $\{\pm1\}$ である（$i$ と可換な元は $\pm1,\pm i$、$j$ と可換な元は $\pm1,\pm j$ で、両方と可換なのは $\pm1$）。従って $\abs{\operatorname{Inn}(Q_8)}=8/2=4$ である。$Q_8/\{\pm1\}$ の非単位元 $\{\pm i\},\{\pm j\},\{\pm k\}$ はいずれも平方が $\{\pm1\}$ なので位数 $2$ であり、$\operatorname{Inn}(Q_8)\cong C_2\times C_2$ である（:ref[prop-grp-order-four]）。
2. $D_4$ の中心は $\{e,r^2\}$ である。$r^2$ は $r$ と可換で、$sr^2s^{-1}=r^{-2}=r^2$ より $s$ とも可換である。$r,r^3$ は $srs^{-1}=r^{-1}\ne r$ より $s$ と可換でなく、鏡映 $r^js$ は $r(r^js)r^{-1}=r^{j+2}s\ne r^js$ より $r$ と可換でない。従って $\abs{\operatorname{Inn}(D_4)}=4$ である。
3. $Z(S_3)=\{e\}$ なので $\operatorname{Inn}(S_3)\cong S_3$ である（:ref[exm-grp-aut-small]）。
4. $G$ が可換なら $Z(G)=G$ で $\operatorname{Inn}(G)$ は自明である。
:::

:ref[prop-grp-inn-center] と :ref[prop-grp-center-quotient-cyclic] を合わせると、$\operatorname{Inn}(G)$ が巡回群なら $G$ は可換、従って $\operatorname{Inn}(G)$ は自明である。特に $\operatorname{Inn}(G)$ は位数 $2$ や $3$ にはなりえない。

### Cayley の定理

第一同型定理の単射版として、どの群も置換群として実現できることを示す。

:::example[Cayley の定理]{#exm-grp-cayley-embedding}
群 $G$ の各元 $g$ に、左からの掛け算 $\lambda_g:G\to G$, $\lambda_g(x)=gx$ を対応させる。$\lambda_g$ は $\lambda_{g^{-1}}$ を逆写像にもつ全単射であり、$\lambda_g\lambda_h(x)=ghx=\lambda_{gh}(x)$ なので、$\lambda:G\to S_G$ は準同型である。$\lambda_g=\id$ なら $g=\lambda_g(e)=e$ なので核は自明であり、第一同型定理により $G\cong\lambda(G)\le S_G$ である。$\abs G=n$ なら、$G$ の元に $1,\ldots,n$ と番号をつけることで $G$ は $S_n$ の部分群と同型になる。例えば $C_4$ の $\bar1$ は $\bar0\mapsto\bar1\mapsto\bar2\mapsto\bar3\mapsto\bar0$ と動かすので、四文字の巡回置換に対応する。群作用の言葉による一般化は [群の作用](/math/group-theory/group-actions/) で扱う。
:::

### $S_4$ から $S_3$ への全射

第一同型定理の効果がよく分かる例として、$S_4$ の剰余群 $S_4/V_4$ を同定する。

:::example[四次対称群を Klein の四元群で割る]{#exm-grp-s4-klein-quotient}
集合 $\{1,2,3,4\}$ を二元ずつの二組に分ける方法は三通りある。

$$
\pi_1=\{\{1,2\},\{3,4\}\},\qquad \pi_2=\{\{1,3\},\{2,4\}\},\qquad \pi_3=\{\{1,4\},\{2,3\}\}.
$$

$\sigma\in S_4$ は組分け $\{\{a,b\},\{c,d\}\}$ を組分け $\{\{\sigma(a),\sigma(b)\},\{\sigma(c),\sigma(d)\}\}$ に移す。これは $P=\{\pi_1,\pi_2,\pi_3\}$ 上の置換 $\Phi(\sigma)$ を定め（$\sigma^{-1}$ が逆の置換を与える）、$\tau$ で移してから $\sigma$ で移すことは $\sigma\tau$ で移すことと同じなので、$\Phi:S_4\to S_P\cong S_3$ は準同型である。

**全射性**：$(1\,2)$ は $\pi_1$ を固定し、$\pi_2=\{\{1,3\},\{2,4\}\}$ を $\{\{2,3\},\{1,4\}\}=\pi_3$ に移すので、$\Phi((1\,2))$ は $\pi_2,\pi_3$ の互換である。$(2\,3)$ は $\pi_1$ を $\{\{1,3\},\{2,4\}\}=\pi_2$ に移し、$\pi_3$ を固定するので、$\Phi((2\,3))$ は $\pi_1,\pi_2$ の互換である。$S_P$ の二つの異なる互換の積は長さ $3$ の巡回置換なので、像は位数 $2$ と $3$ の元を含み、Lagrange の定理により位数は $6$ の倍数、すなわち $S_P$ 全体である。

**核**：:ref[cor-grp-kernel-image-order] より $\abs{\Ker\Phi}=24/6=4$ である。$(1\,2)(3\,4)$ は $\{1,2\}$ と $\{3,4\}$ をそれぞれ保つので $\pi_1$ を固定し、$\pi_2$ を $\{\{2,4\},\{1,3\}\}=\pi_2$ に、$\pi_3$ を $\{\{2,3\},\{1,4\}\}=\pi_3$ に移す。同様に $(1\,3)(2\,4)$ と $(1\,4)(2\,3)$ も三つの組分けをすべて固定する（例えば $(1\,3)(2\,4)$ は $\pi_1$ を $\{\{3,4\},\{1,2\}\}$ に移す）。従って $V_4=\{e,(1\,2)(3\,4),(1\,3)(2\,4),(1\,4)(2\,3)\}\subseteq\Ker\Phi$ であり（:ref[exm-grp-klein-four]）、位数が等しいので $\Ker\Phi=V_4$ である。

第一同型定理により

$$
S_4/V_4\cong S_3 \label{eq-grp-s4-klein}
$$

である。副産物として $V_4\trianglelefteq S_4$ も分かる。また $\Phi((1\,2\,3\,4))$ は $\pi_1\leftrightarrow\pi_3$ の互換（$\pi_2$ を固定）であり、$\Phi((1\,2\,3))$ は $\pi_1\mapsto\pi_3\mapsto\pi_2\mapsto\pi_1$ と巡る。$\Phi(A_4)$ は $\abs{A_4}/\abs{V_4}=3$ 個の元からなるので（$V_4\subseteq A_4$ に注意）、$S_P$ の位数 $3$ の部分群、すなわち偶置換全体である。
:::

この準同型は、$S_4$ が「三つの組分けの対称性」を通じて $S_3$ に写ることを示している。$n\ge5$ では $S_n$ から $S_{n-1}$ への全射準同型は存在しない。$S_n$ の正規部分群は $\{e\},A_n,S_n$ に限り（[対称群](/math/group-theory/symmetric-groups/) で示す $A_n$ の単純性から従う）、剰余群の位数は $n!,2,1$ のいずれかで $(n-1)!$ にならないからである。$S_4$ のこの特殊性は、四次方程式の三次分解式の背後にある構造でもある（[体と Galois 理論](/math/fields-galois/)）。

## 第二・第三同型定理

第一同型定理は一つの準同型についての定理であった。ここでは、部分群に制限してから割る操作と、段階的に割る操作を比較する。どちらも、適切な準同型を作って第一同型定理を適用するだけで証明できる。

### 第二同型定理

大きな群で割った後に部分群だけを見る操作と、部分群の中で直接割る操作を比べる。

:::theorem[第二同型定理]{#thm-grp-second-iso}
$H\le G$, $N\trianglelefteq G$ とする。このとき $HN\le G$, $N\trianglelefteq HN$, $H\cap N\trianglelefteq H$ であり、

$$
H/(H\cap N)\ \cong\ HN/N,\qquad h(H\cap N)\longmapsto hN
$$

が成り立つ。特に $[HN:N]=[H:H\cap N]$ である。
:::

:::proof
$N$ が正規なので $HN$ は部分群である（:ref[prop-grp-product-subgroup]）。$N=eN\subseteq HN$ であり、$HN$ の元による共役は $G$ の元による共役なので $N$ を保つ。従って $N\trianglelefteq HN$ であり、剰余群 $HN/N$ が定まる。

写像 $\phi:H\to HN/N$, $\phi(h)=hN$ は、自然射影 $\pi_N:HN\to HN/N$ を $H$ に制限したものなので準同型である。$HN/N$ の元は $(hn)N=hN=\phi(h)$（$h\in H$, $n\in N$）と書けるので $\phi$ は全射である。$\phi(h)=N$ は $h\in N$ と同値なので、$\Ker\phi=H\cap N$ である。従って $H\cap N$ は $H$ の正規部分群であり（:ref[prop-grp-image-preimage]）、第一同型定理により $H/(H\cap N)\cong HN/N$ である。同型な群の位数は等しいので、指数についての等式も従う。
:::

有限群では $\abs{HN}/\abs N=\abs H/\abs{H\cap N}$ となり、:ref[prop-grp-product-size] の公式 $\abs{HN}=\abs H\abs N/\abs{H\cap N}$ が再び得られる。第二同型定理の内容は、$H$ を $G/N$ に写したときに失われるのは $H\cap N$ だけであり、像は $HN/N$ である、ということである。$H\cap N$ は $H$ では正規だが、$G$ で正規とは限らない点に注意する。

:::example[二つの倍数の群]
$\Z$ で $H=4\Z$, $N=6\Z$ とする。$H+N$ は $4$ と $6$ を含むので $6-4=2$ を含み、$2\Z$ に含まれるので $H+N=2\Z$ である。$H\cap N$ は $4$ と $6$ の公倍数全体で $12\Z$ である。従って

$$
4\Z/12\Z\cong2\Z/6\Z
$$

を得る。左辺の代表 $0,4,8$ は右辺の類 $0,4,2$ にそれぞれ移り、対応は加法を保つ。どちらも位数 $3$ の巡回群である。
:::

:::example[四次対称群の一点固定部分群と Klein の四元群]
$S_4$ で $4$ を固定する置換全体を $H$ とすると、$H$ は $\{1,2,3\}$ の置換全体と同一視でき $H\cong S_3$ である。$N=V_4$（:ref[exm-grp-s4-klein-quotient] により正規）の非単位元はすべて $4$ を動かすので $H\cap N=\{e\}$ である。:ref[prop-grp-product-size] より $\abs{HN}=6\cdot4/1=24$ なので $HN=S_4$ であり、第二同型定理から

$$
S_3\cong H/\{e\}\cong S_4/V_4
$$

となる。これは :eqref[eq-grp-s4-klein] の別証明であるとともに、$\Phi$ を $H$ に制限したものが同型であることを示している。実際、組分け $\pi_1,\pi_2,\pi_3$ は $4$ の相手 $3,2,1$ でそれぞれ決まり、$4$ を固定する $\sigma$ は $4$ の相手 $a$ を $\sigma(a)$ に移すので、$\Phi|_H$ は $\{1,2,3\}$ 上の置換そのものと対応する。
:::

:::example[射影特殊線形群]
$G=\GL_2(\F_5)$ とし、$H=\SL_2(\F_5)$、$N=\{aI:a\in\F_5^\times\}$ をスカラー行列全体とする。スカラー行列はすべての行列と可換なので $N\trianglelefteq G$ である。$H\cap N=\{aI:a^2=1\}=\{\pm I\}$ なので、第二同型定理により

$$
\SL_2(\F_5)/\{\pm I\}\cong\SL_2(\F_5)N/N
$$

である。:ref[exm-grp-kernel-examples] の 1 より $\abs{\SL_2(\F_5)}=120$ なので、左辺の位数は $60$ である。$\abs{\SL_2(\F_5)N}=120\cdot4/2=240$ であり、これは $\abs{\GL_2(\F_5)}=480$ の半分である。左辺の群を $\mathrm{PSL}_2(\F_5)$ と書き、**射影特殊線形群**という。第二同型定理は、それが「行列のスカラー倍を同一視した群」$\GL_2(\F_5)/N$ の部分群として実現されることを示している。
:::

### 第三同型定理

$N$ で割った後にさらに大きな $M$ で割ることは、はじめから $M$ で割ることと同じである。

:::theorem[第三同型定理]{#thm-grp-third-iso}
$N\subseteq M$ をともに $G$ の正規部分群とする。このとき $M/N\trianglelefteq G/N$ であり、

$$
(G/N)/(M/N)\cong G/M,\qquad (gN)(M/N)\longmapsto gM
$$

が成り立つ。特に $G$ が有限なら $[G:M]=[G:N]/[M:N]$ である。
:::

:::proof
自然射影 $\pi_M:G\to G/M$ の核 $M$ は $N$ を含むので、:ref[prop-grp-quotient-universal] により $\psi(gN)=gM$ で定まる準同型 $\psi:G/N\to G/M$ が存在する。$\psi$ は明らかに全射である。$\psi(gN)=M$ は $g\in M$ と同値であり、$g\in M$ と $gN\in M/N=\{mN:m\in M\}$ は同値である（$gN=mN$ なら $g\in mN\subseteq M$）。従って $\Ker\psi=M/N$ であり、$M/N$ は $G/N$ の正規部分群である。第一同型定理を $\psi$ に適用すれば同型を得る。位数の等式は $\abs{G/M}=\abs{G/N}/\abs{M/N}$ から従う。
:::

:::example[段階的な合同算術]
$12\Z\subseteq4\Z\subseteq\Z$ について、第三同型定理は

$$
(\Z/12\Z)/(4\Z/12\Z)\cong\Z/4\Z
$$

を与える。$C_{12}$ から $C_4$ への自然な写像 $\bar x\mapsto\bar x$ の核は $\{\bar0,\bar4,\bar8\}=4\Z/12\Z$ である。最初から法 $4$ で見ることと、法 $12$ で見た後に $\bar4$ を零とみなすことは同じである。
:::

:::example[四次対称群の正規部分群の列]
$S_4\supseteq A_4\supseteq V_4\supseteq\{e\}$ を考える。$A_4=\Ker\sgn$ と $V_4=\Ker\Phi$（:ref[exm-grp-s4-klein-quotient]）はともに $S_4$ の正規部分群であり、$V_4$ の非単位元は二つの互換の積なので偶置換、すなわち $V_4\subseteq A_4$ である。第三同型定理により

$$
(S_4/V_4)/(A_4/V_4)\cong S_4/A_4\cong\{\pm1\}
$$

である。$A_4/V_4$ は位数 $12/4=3$ なので巡回群であり（:ref[cor-grp-order-divides]）、同型 $S_4/V_4\cong S_3$ のもとで $S_3$ の位数 $3$ の部分群 $\langle(1\,2\,3)\rangle$ に対応する。左辺は $S_3/\langle(1\,2\,3)\rangle\cong C_2$ と読める。こうして $S_4$ は、順に割ると位数 $2,3,4$ の「層」に分解される。この見方は可解群の理論で組成列として発展する（[可解群と冪零群](/math/group-theory/solvable-nilpotent/)）。
:::

### 商の間の写像

自然射影に限らず、一般の準同型 $f:G\to H$ を両側の剰余群の間に下ろしたい場合がある。

:::proposition[剰余群の間に誘導される準同型]{#prop-grp-induced-map}
$f:G\to H$ を準同型とし、$N\trianglelefteq G$, $M\trianglelefteq H$ とする。

1. $\bar f(gN)=f(g)M$ によって写像 $\bar f:G/N\to H/M$ が（代表元によらずに）定まるための必要十分条件は $f(N)\subseteq M$ である。このとき $\bar f$ は準同型である。
2. このとき $\Ker\bar f=f^{-1}(M)/N$, $\Img\bar f=f(G)M/M$ である。従って $\bar f$ が単射であることは $f^{-1}(M)=N$ と、全射であることは $f(G)M=H$ と同値である。
:::

:::proof
$F=\pi_M\circ f:G\to H/M$ とおくと、$F(g)=M$ は $f(g)\in M$ と同値なので $\Ker F=f^{-1}(M)$ である。式が代表元によらず定まるなら、$n\in N$ に対して $nN=eN$ なので $f(n)M=f(e)M=M$ となり、$f(N)\subseteq M$ が必要である。逆に誘導したい準同型は $\bar f\circ\pi_N=F$ を満たすものなので、:ref[prop-grp-quotient-universal] によりそれが存在する条件は $N\subseteq f^{-1}(M)$、すなわち $f(N)\subseteq M$ であり、そのとき $\bar f$ は準同型である。同じ命題により $\Ker\bar f=\Ker F/N=f^{-1}(M)/N$, $\Img\bar f=\Img F=\{f(g)M:g\in G\}$ である。$M$ が正規なので $f(G)M$ は $H$ の部分群であり（:ref[prop-grp-product-subgroup]）、$\{f(g)M\}$ は $f(G)M$ に含まれる $M$ の剰余類全体、すなわち $f(G)M/M$ に等しい。単射性は $f^{-1}(M)/N$ が自明であること、すなわち $f^{-1}(M)\subseteq N$ と同値であり、$N\subseteq f^{-1}(M)$ と合わせて等号と同値である。全射性は、任意の $hM$ が $f(g)M$ と書けること、すなわち $H=f(G)M$ と同値である。
:::

:::example[整数の倍写像を剰余群に下ろす]
$f:\Z\to\Z$, $x\mapsto6x$ を $C_{10}\to C_{15}$ に下ろす。$f(10\Z)=60\Z\subseteq15\Z$ なので $\bar f(\bar x)=\overline{6x}$ が定まる。核は $15\mid6x$、すなわち $5\mid2x$、すなわち $5\mid x$ となる類で、法 $10$ の $\bar0,\bar5$ である。像は法 $15$ の $\bar0,\bar3,\bar6,\bar9,\overline{12}$ であり、:eqref[eq-grp-kernel-image-order] の通り $10=2\cdot5$ となる。元の整数上の写像は単射だが、$f^{-1}(15\Z)=5\Z\ne10\Z$ なので誘導写像は単射でない。
:::

### 二つの剰余群への同時の写像

二つの正規部分群 $N,M\trianglelefteq G$ で同じ元 $g$ を観測し、記録の組 $(gN,gM)$ を作る。どのような組が実際に現れるかを決定する。直積 $G/N\times G/M$ は成分ごとの積で群をなす（:ref[exm-grp-direct-product-basic]）。

:::proposition[二つの剰余群への埋め込み]{#prop-grp-two-quotients}
$N,M\trianglelefteq G$ とする。$\psi:G\to G/N\times G/M$, $\psi(g)=(gN,gM)$ は準同型で $\Ker\psi=N\cap M$ である。従って

$$
G/(N\cap M)\cong\Img\psi=\{(aN,bM):a^{-1}b\in NM\}
$$

であり、条件 $a^{-1}b\in NM$ は $aNM=bNM$ と同値である。特に $NM=G$ ならば $G/(N\cap M)\cong G/N\times G/M$ である。
:::

:::proof
$\psi(gh)=(ghN,ghM)=(gN\cdot hN,\ gM\cdot hM)=\psi(g)\psi(h)$ である。$\psi(g)$ が単位元 $(N,M)$ になるのは $g\in N$ かつ $g\in M$ のときなので、$\Ker\psi=N\cap M$ である。第一同型定理より $G/(N\cap M)\cong\Img\psi$ である。

像を決定する。$(aN,bM)\in\Img\psi$ は、$gN=aN$ かつ $gM=bM$ となる $g$、すなわち $g\in aN\cap bM$ が存在することと同値である。そのような $g$ があれば $g=an=bm$（$n\in N$, $m\in M$）と書け、$a^{-1}b=nm^{-1}\in NM$ である。逆に $a^{-1}b=nm'$（$n\in N$, $m'\in M$）なら $g=an=bm'^{-1}$ は $aN\cap bM$ に属する。$NM$ は正規部分群なので（:ref[prop-grp-normal-operations]）、:ref[prop-grp-coset-test] により $a^{-1}b\in NM$ は $aNM=bNM$ と同値である。$NM=G$ なら条件は常に満たされ、$\psi$ は全射である。
:::

有限群では、像の位数を比べて $[G:N\cap M]=[G:N][G:M]/[G:NM]$ が得られる。像は直積 $G/N\times G/M$ のうち「共通の剰余群 $G/NM$ で一致する組」全体であり、$NM=G$ ならこの制約が消える。

:::corollary[中国剰余定理]{#cor-grp-chinese-remainder}
$m,n\ge1$ を互いに素な整数とする。このとき $a\bmod mn\mapsto(a\bmod m,\ a\bmod n)$ は同型

$$
\Z/mn\Z\cong\Z/m\Z\times\Z/n\Z
$$

を与える。
:::

:::proof
:ref[prop-grp-two-quotients] を $G=\Z$, $N=m\Z$, $M=n\Z$ に適用する。$um+vn=1$ となる整数 $u,v$ をとると、任意の $x\in\Z$ は $x=(xu)m+(xv)n\in m\Z+n\Z$ なので $N+M=\Z$ である。次に $x\in m\Z\cap n\Z$ とすると $x=umx+vnx$ であり、$n\mid x$ より $mn\mid umx$、$m\mid x$ より $mn\mid vnx$ なので $mn\mid x$ である。逆に $mn$ の倍数は両方に属するので $N\cap M=mn\Z$ である。従って $\Z/mn\Z\cong\Z/m\Z\times\Z/n\Z$ であり、同型は $\psi$ から誘導される写像 $a+mn\Z\mapsto(a+m\Z,\,a+n\Z)$ である。
:::

例えば $\Z/12\Z\cong\Z/3\Z\times\Z/4\Z$ で、法 $3$ で $2$、法 $4$ で $1$ となる組 $(\bar2,\bar1)$ は法 $12$ の $\bar5$ に対応する。有限 Abel 群の分解における役割は [直積・半直積と有限 Abel 群](/math/group-theory/products-abelian/) で扱う。

:::example[互いに素でない法の場合]
$G=\Z$, $N=8\Z$, $M=12\Z$ とする。$N\cap M=24\Z$, $N+M=4\Z$ なので、:ref[prop-grp-two-quotients] により $C_{24}$ は $C_8\times C_{12}$ の中の「法 $4$ の余りが等しい組」全体と同型である。位数は $8\cdot12/4=24$ で整合する。例えば $(\bar5,\bar9)$ は $5\equiv9\equiv1\pmod4$ なので像に入る。実際 $5+8t\equiv9\pmod{12}$ を解くと $8t\equiv4$、すなわち $2t\equiv1\pmod3$ より $t\equiv2\pmod3$ であり、$5+16=21$ が法 $8$ で $5$、法 $12$ で $9$ となる。この組を与える元は法 $24$ の $\overline{21}$ である。一方 $(\bar5,\bar6)$ は法 $4$ の余りが $1$ と $2$ で異なるので、逆像は空である。
:::

## 対応定理

剰余群 $G/N$ の部分群を調べるとき、剰余類の集合として扱う代わりに、元の群 $G$ の部分群として計算できる。ただし、一つの剰余類を含めるならその類に属する元をすべて含めることになるので、対応するのは $N$ を含む部分群だけである。

:::theorem[部分群の対応定理]{#thm-grp-correspondence}
$N\trianglelefteq G$ とし、$\pi=\pi_N$ を自然射影とする。$N$ を含む $G$ の部分群全体を $\mathcal S(G,N)$、$G/N$ の部分群全体を $\mathcal S(G/N)$ と書く。

1. $\alpha(H)=H/N=\pi(H)$ と $\beta(L)=\pi^{-1}(L)$ は互いに逆な全単射 $\alpha:\mathcal S(G,N)\to\mathcal S(G/N)$, $\beta:\mathcal S(G/N)\to\mathcal S(G,N)$ を与える。
2. $H_1,H_2\in\mathcal S(G,N)$ について、$H_1\subseteq H_2$ と $H_1/N\subseteq H_2/N$ は同値である。また

$$
(H_1\cap H_2)/N=H_1/N\cap H_2/N,\qquad \langle H_1\cup H_2\rangle/N=\langle H_1/N\cup H_2/N\rangle
$$

が成り立つ。

3. $N\le H\le K\le G$ ならば $[K:H]=[K/N:H/N]$ である。
4. $N\le H\le K\le G$ ならば、$H\trianglelefteq K$ と $H/N\trianglelefteq K/N$ は同値である。このとき $(K/N)/(H/N)\cong K/H$ である。
:::

:::proof
1. $H\in\mathcal S(G,N)$ なら $\pi(H)=\{hN:h\in H\}=H/N$ は $G/N$ の部分群である（:ref[prop-grp-image-preimage]）。$L\in\mathcal S(G/N)$ なら $\pi^{-1}(L)$ は $G$ の部分群であり、$L$ が単位元 $N$ を含むので $\pi^{-1}(L)\supseteq\pi^{-1}(\{N\})=\Ker\pi=N$ である。よって $\alpha,\beta$ は写像として定まる。

$\pi$ は全射なので、$L$ の各元は $\pi(g)$ と書け、その $g$ は $\pi^{-1}(L)$ に属する。従って $\alpha(\beta(L))=\pi(\pi^{-1}(L))=L$ である。逆向きについて、$\beta(\alpha(H))=\pi^{-1}(\pi(H))\supseteq H$ は明らかである。$x\in\pi^{-1}(\pi(H))$ なら $xN=hN$ となる $h\in H$ があり、$x\in hN\subseteq H$（$N\subseteq H$ を用いた）である。よって $\beta(\alpha(H))=H$ である。

2. $H_1\subseteq H_2$ なら $\pi(H_1)\subseteq\pi(H_2)$ である。逆に $\pi(H_1)\subseteq\pi(H_2)$ なら、1 より $H_1=\beta(\pi(H_1))\subseteq\beta(\pi(H_2))=H_2$ である。

共通部分について：逆像は共通部分と交換するので $\beta(L_1\cap L_2)=\beta(L_1)\cap\beta(L_2)$ である。$L_i=\alpha(H_i)$ とおくと $\beta(\alpha(H_1)\cap\alpha(H_2))=H_1\cap H_2$ であり、$H_1\cap H_2\in\mathcal S(G,N)$ なので両辺に $\alpha$ を施して $\alpha(H_1)\cap\alpha(H_2)=\alpha(H_1\cap H_2)$ を得る。

生成される部分群について：$J=\langle H_1\cup H_2\rangle$ は $N$ を含むので $J\in\mathcal S(G,N)$ である。$\alpha(J)\supseteq\alpha(H_1),\alpha(H_2)$ なので $\alpha(J)\supseteq\langle\alpha(H_1)\cup\alpha(H_2)\rangle$ である。逆に $L=\langle\alpha(H_1)\cup\alpha(H_2)\rangle$ とすると、前半より $\beta(L)\supseteq\beta(\alpha(H_i))=H_i$ なので $\beta(L)\supseteq J$ であり、$L=\alpha(\beta(L))\supseteq\alpha(J)$ である。よって $\alpha(J)=L$ である。

3. $K/N$ は $G/N$ の部分群で、2 より $H/N\le K/N$ である。写像 $\kappa:K/H\to(K/N)/(H/N)$（左剰余類の集合の間の写像）を $\kappa(kH)=(kN)(H/N)$ で定めたい。$k,k'\in K$ について

$$
(kN)(H/N)=(k'N)(H/N)\iff(kN)^{-1}(k'N)=k^{-1}k'N\in H/N\iff k^{-1}k'\in\beta(H/N)=H\iff kH=k'H
$$

である（:ref[prop-grp-coset-test]）。$kH=k'H$ なら像が等しくなるので $\kappa$ は代表元によらずに定まる。逆に像が等しければ $kH=k'H$ となるので単射である。$K/N$ の元はすべて $kN$（$k\in K$）の形なので $\kappa$ は全射である。従って二つの剰余類の集合は同じ濃度をもち、$[K:H]=[K/N:H/N]$ である。

4. $\pi$ を $K$ に制限した $\pi|_K:K\to K/N$ は全射準同型である。$H\trianglelefteq K$ なら、:ref[prop-grp-image-preimage] の 3 により $\pi|_K(H)=H/N$ は $\pi|_K(K)=K/N$ の正規部分群である。逆に $H/N\trianglelefteq K/N$ なら、同じ命題の 2 により $(\pi|_K)^{-1}(H/N)=K\cap\pi^{-1}(H/N)=K\cap H=H$ は $K$ の正規部分群である。このとき $N\subseteq H$ はともに $K$ の正規部分群なので（$N\trianglelefteq G$ より $N\trianglelefteq K$）、群 $K$ に第三同型定理を適用して $(K/N)/(H/N)\cong K/H$ を得る。
:::

対応定理は、$G/N$ の部分群の「格子」（包含関係による図）が、$G$ の部分群の格子のうち $N$ より上の部分と同じ形をしていることを述べている。$N$ を含まない部分群 $H$ については、$\pi(H)=HN/N$ に対応するのは $H$ ではなく $HN$ である（第二同型定理）。

:::corollary[正規部分群の対応と単純な剰余群]{#cor-grp-correspondence-normal}
$N\trianglelefteq G$ とする。

1. $H\mapsto H/N$ は、$N$ を含む $G$ の正規部分群全体から $G/N$ の正規部分群全体への全単射であり、$G/H\cong(G/N)/(H/N)$ が成り立つ。
2. $G/N$ が自明でなく、正規部分群として $\{e\}$ と自分自身しかもたないこと（このような群を**単純群**という。[対称群](/math/group-theory/symmetric-groups/) を参照）は、$N\ne G$ かつ $N\subsetneq M\subsetneq G$ を満たす正規部分群 $M$ が存在しないことと同値である。
:::

:::proof
1. :ref[thm-grp-correspondence] の 1 と、4 を $K=G$ に適用したものから従う。

2. 1 の全単射は包含を保つので、$N$ と $G$ の間に真に挟まる正規部分群は、$G/N$ の自明でも全体でもない正規部分群に対応する。$G/N$ が自明であることは $N=G$ と同値である。
:::

:::example[巡回群の部分群]
$G=\Z$, $N=12\Z$ とする。$\Z$ の部分群は $d\Z$（$d\ge0$）の形であり（:ref[thm-grp-cyclic-subgroups]）、$12\Z\subseteq d\Z$ は $d\mid12$ と同値である。従って $C_{12}=\Z/12\Z$ の部分群は $d\Z/12\Z$（$d=1,2,3,4,6,12$）の六つであり、:ref[thm-grp-correspondence] の 3 よりその指数は $[\Z:d\Z]=d$、位数は $12/d$ である。包含は整数の大小ではなく整除で決まり、$d\Z\subseteq d'\Z$ は $d'\mid d$ と同値である。例えば $4\Z/12\Z$ と $6\Z/12\Z$ の共通部分は $\operatorname{lcm}(4,6)\Z/12\Z=12\Z/12\Z$、合わせて生成する部分群は $\gcd(4,6)\Z/12\Z=2\Z/12\Z$ である。
:::

:::example[対応から外れる部分群]
$G=C_2\times C_2$ で第二成分を忘れる射影 $G\to C_2$, $(x,y)\mapsto x$ を考える。核は $N=\{\bar0\}\times C_2$ である。部分群 $C_2\times\{\bar0\}$ と対角部分群 $\{(\bar0,\bar0),(\bar1,\bar1)\}$ はどちらも $N$ を含まず、どちらも像は $C_2$ 全体である。対応定理の条件 $N\le H$ を外すと、像から元の部分群を復元できない。どちらの $H$ についても $HN=G$ であり、像 $C_2$ に対応定理で対応する $N$ を含む部分群は $G$ 自身である。
:::

:::example[四次対称群で Klein の四元群を含む部分群]{#exm-grp-correspondence-s4}
:ref[exm-grp-s4-klein-quotient] の全射 $\Phi:S_4\to S_3$（核 $V_4$）に対応定理を適用する。$S_3$ の部分群は、自明群、互換が生成する位数 $2$ の三つ、$\langle(1\,2\,3)\rangle$、$S_3$ の六つである。従って $V_4$ を含む $S_4$ の部分群もちょうど六つあり、位数はそれぞれ $4,8,8,8,12,24$ である。

- 自明群には $V_4$ が、全体には $S_4$ が対応する。
- 位数 $3$ の部分群には位数 $12$ の部分群が対応する。$A_4\supseteq V_4$ で $\Phi(A_4)$ は位数 $3$（:ref[exm-grp-s4-klein-quotient]）なので、これは $A_4$ である。
- 互換 $\Phi((1\,2))$（$\pi_1$ を固定）が生成する部分群には、$\pi_1$ を固定する置換全体が対応する。それは $V_4\cup(1\,2)V_4$ であり、

$$
\{e,\ (1\,2)(3\,4),\ (1\,3)(2\,4),\ (1\,4)(2\,3),\ (1\,2),\ (3\,4),\ (1\,3\,2\,4),\ (1\,4\,2\,3)\}
$$

である。例えば $(1\,2)\cdot(1\,3)(2\,4)$ は $1\mapsto3$, $3\mapsto2$, $2\mapsto4$, $4\mapsto1$ と送るので $(1\,3\,2\,4)$ である。これは頂点を $1,3,2,4$ の順に並べた正方形の対称性の群で、$D_4$ と同型である。$\pi_2,\pi_3$ を固定する置換全体も同様に位数 $8$ の部分群をなす。

$S_3$ の正規部分群は自明群、$\langle(1\,2\,3)\rangle$、$S_3$ の三つなので、:ref[cor-grp-correspondence-normal] により $V_4$ を含む $S_4$ の正規部分群は $V_4,A_4,S_4$ の三つであり、位数 $8$ の三つは正規でない。
:::

### どの定理で部分群を比較するか

$G$ からの全射が一つ与えられたとき、像そのものを同定したければ第一同型定理を使う。同じ写像の定義域を部分群 $H$ に絞ると、第二同型定理の $H\cap N$ が現れる。$N$ より多くの元を単位元とみなす $M$ への変更は、第三同型定理になる。これらの違いを一つの群の中で追ってみる。

$G=D_6$ とし、:ref[exm-grp-dihedral-quotient-map] の全射 $f:D_6\to D_3$, $r^is^j\mapsto r'^is'^j$ を考える。核は $N=\langle r^3\rangle=\{e,r^3\}$ なので、第一同型定理により $D_6/\langle r^3\rangle\cong D_3$ である。幾何学的には、正六角形の向かい合う頂点を一組にすると三つの組が正三角形の頂点のように置換される。番号を法 $3$ に落とせば $r$ は一つ進み、$s$ は符号を反転する。

次に $H=\langle r^2,s\rangle=\{r^{2k},r^{2k}s:0\le k<3\}$ に制限する。$H$ の回転指数は偶数なので $H\cap N=\{e\}$ であり、第二同型定理により $H\cong HN/N$ である。$\abs{HN}=6\cdot2=12$ なので $HN=G$ であり、右辺は $G/N\cong D_3$ そのものになる。つまり $H$ は正六角形の頂点 $0,2,4$ がなす正三角形の対称性の群で、$f|_H$ は同型である。

最後に $M=\langle r\rangle$ をとると $N\le M$ であり、$M/N$ は $D_3$ の回転部分群 $\langle r'\rangle$ に対応する。第三同型定理により $(G/N)/(M/N)\cong G/M\cong C_2$ であり、回転か鏡映かだけを記録する群が残る。これははじめから $G/M$ をとった結果と一致する。こうして、同じ像への制限（第二同型定理）と、像をさらに割る操作（第三同型定理）を区別できる。対応定理により、$D_3$ の六つの部分群は $D_6$ の $N=\langle r^3\rangle$ を含む六つの部分群に一対一に対応する。

## 演習問題

::::exercise{level=basic}
$C_{10}$ から $C_{15}$ への準同型をすべて求めよ。各準同型の核と像も答えよ。

:::solution
:ref[prop-grp-hom-from-cyclic] により、準同型は $\bar1$ の行き先 $y\in C_{15}$ で $10y=\bar0$ を満たすものと一対一に対応し、$f(\bar x)=xy$ で与えられる。$15\mid10y$ は $3\mid2y$、すなわち $3\mid y$ と同値なので、$y=\bar0,\bar3,\bar6,\bar9,\overline{12}$ の五通りであり、$\gcd(10,15)=5$ と一致する。

$y=\bar0$ なら自明な準同型で、核は $C_{10}$、像は $\{\bar0\}$ である。$y\ne\bar0$ のとき $y=\overline{3c}$（$c=1,2,3,4$）と書ける。$15\mid3cx$ は $5\mid cx$、すなわち $5\mid x$ と同値なので（$\gcd(c,5)=1$）、核は $\{\bar0,\bar5\}$ である。像は位数 $10/2=5$ の部分群で（:ref[cor-grp-kernel-image-order]）、$C_{15}$ の位数 $5$ の部分群はただ一つなので（:ref[thm-grp-cyclic-subgroups]）$\langle\bar3\rangle=\{\bar0,\bar3,\bar6,\bar9,\overline{12}\}$ である。
:::
::::

::::exercise{level=basic}
$f:\Z^2\to\Z$, $(a,b)\mapsto2a+3b$ の核と像を求め、$\Z^2/\Ker f$ を同定せよ。ただし $\Z^2$ は成分ごとの和で群とみなす。

:::solution
$f((a,b)+(c,d))=2(a+c)+3(b+d)=f(a,b)+f(c,d)$ なので準同型である。$f(-1,1)=1$ より像は $1$ を含み、従って $\Z$ 全体である。$2a+3b=0$ なら $3\mid2a$ より $3\mid a$ なので $a=3k$ とおけて $b=-2k$ である。逆に $(3k,-2k)$ は核に入る。よって $\Ker f=\langle(3,-2)\rangle$ であり、第一同型定理により $\Z^2/\langle(3,-2)\rangle\cong\Z$ である。
:::
::::

::::exercise{level=basic}
$D_5$ の元 $r^is^j$ を $\bar j\in C_2$ に送る写像が準同型であることを示し、その核を求めよ。

:::solution
:ref[prop-grp-hom-dihedral] を $A=\bar0$, $B=\bar1$ に適用する（加法で書く）。$5A=\bar0$, $2B=\bar0$, $B+A-B=\bar0=-A$ なので条件を満たし、$r^is^j\mapsto iA+jB=\bar j$ は準同型である。核は $j=0$ の元全体、すなわち回転部分群 $\langle r\rangle$ であり、位数 $5$ である。像は $C_2$ 全体なので、:eqref[eq-grp-kernel-image-order] の通り $10=5\cdot2$ である。
:::
::::

::::exercise{level=basic}
自然射影 $\Z\to\Z/4\Z=C_4$ について、$C_4$ の位数 $2$ の部分群の逆像を求めよ。

:::solution
$C_4$ の位数 $2$ の部分群は $\{\bar0,\bar2\}$ である。この二つの類に入る整数は $4$ で割った余りが $0$ または $2$ の整数、すなわち偶数全体であり、逆像は $2\Z$ である。対応定理（:ref[thm-grp-correspondence]）の通り $2\Z\supseteq4\Z$ であり、$[\Z:2\Z]=2=[C_4:\{\bar0,\bar2\}]$ となる。
:::
::::

::::exercise{level=basic}
$\R^\times/\{\pm1\}\cong\R_{>0}$ を示せ。ここで $\R^\times$, $\R_{>0}$ は乗法群である。

:::solution
$f:\R^\times\to\R_{>0}$, $x\mapsto\abs x$ は $\abs{xy}=\abs x\abs y$ により準同型である。$y>0$ は $f(y)=y$ なので全射であり、$\abs x=1$ は $x=\pm1$ と同値なので $\Ker f=\{\pm1\}$ である。第一同型定理により $\R^\times/\{\pm1\}\cong\R_{>0}$ である。各剰余類 $x\{\pm1\}=\{x,-x\}$ は二つの元をもち、正の方は $\abs x$ である。従って各 $y>0$ に剰余類 $\{y,-y\}$ が一つずつ対応する。
:::
::::

::::exercise{level=basic}
1. 位数が素数 $p$ の群 $G$ からの準同型は、自明であるか単射であることを示せ。
2. $C_3$ から $S_3$ への準同型をすべて求め、それぞれ自明か単射かを判定せよ。
3. $C_7$ から $S_4$ への準同型をすべて求めよ。

:::solution
1. $\Ker f$ は $G$ の部分群なので、Lagrange の定理（:ref[thm-grp-lagrange]）によりその位数は $1$ か $p$ である。$1$ なら単射（:ref[prop-grp-fiber]）、$p$ なら $\Ker f=G$ で自明である。

2. :ref[prop-grp-hom-from-cyclic] により、準同型は $h^3=e$ を満たす $h\in S_3$ と対応する。$S_3$ の元の位数は $e$ が $1$、互換が $2$、$(1\,2\,3)$, $(1\,3\,2)$ が $3$ なので、$h=e,(1\,2\,3),(1\,3\,2)$ の三つである。$h=e$ は自明であり、残る二つは $\operatorname{ord}h=3$ なので単射である。

3. :ref[cor-grp-kernel-image-order] より像の位数は $\gcd(7,24)=1$ を割り切るので、自明な準同型だけである。
:::
::::

::::exercise{level=basic}
準同型 $f:G\to H$ について、$G$ が巡回群なら $f(G)$ も巡回群であり、$G$ が可換なら $f(G)$ も可換であることを示せ。これを用いて、巡回群の剰余群は巡回群であることを示せ。

:::solution
$G=\langle a\rangle$ なら、:ref[prop-grp-hom-rules] の 1 より $f(G)=\{f(a^m):m\in\Z\}=\{f(a)^m\}=\langle f(a)\rangle$ である。$G$ が可換なら、$f(a)f(b)=f(ab)=f(ba)=f(b)f(a)$ である。$N\le G$ による剰余群 $G/N$ は自然射影 $\pi_N$ の像なので（:ref[prop-grp-natural-projection]）、$G$ が巡回群なら $G/N=\langle aN\rangle$ も巡回群である。
:::
::::

::::exercise{level=advanced}
$\Aut(C_8)$, $\Aut(C_9)$, $\Aut(C_{10})$ のそれぞれについて、巡回群であるかどうかを判定せよ。

:::solution
:ref[thm-grp-aut-cyclic] により $\Aut(C_n)\cong(\Z/n\Z)^\times$ である。

- $(\Z/8\Z)^\times=\{\bar1,\bar3,\bar5,\bar7\}$ では $3^2,5^2,7^2$ がすべて法 $8$ で $1$ なので、位数 $4$ の元がない。従って巡回群でなく、:ref[prop-grp-order-four] により $C_2\times C_2$ と同型である。
- $(\Z/9\Z)^\times=\{\bar1,\bar2,\bar4,\bar5,\bar7,\bar8\}$ では、$2$ の冪が法 $9$ で $2,4,8,7,5,1$ となり、$\bar2$ の位数は $6$ である。従って $C_6$ と同型な巡回群である。
- $(\Z/10\Z)^\times=\{\bar1,\bar3,\bar7,\bar9\}$ では、$3$ の冪が法 $10$ で $3,9,7,1$ となり、$\bar3$ の位数は $4$ である。従って $C_4$ と同型な巡回群である。

位数が同じ $4$ でも、$\Aut(C_8)$ と $\Aut(C_{10})$ は同型でない。
:::
::::

::::exercise{level=advanced}
準同型 $f:G\to H$ の核を $N$ とし、$K\le G$ とする。$f^{-1}(f(K))=KN$ を示せ。

:::solution
$x\in f^{-1}(f(K))$ なら $f(x)=f(k)$ となる $k\in K$ がある。:ref[prop-grp-fiber] より $k^{-1}x\in N$ なので $x=k(k^{-1}x)\in KN$ である。逆に $x=kn$（$k\in K$, $n\in N$）なら $f(x)=f(k)f(n)=f(k)\in f(K)$ である。両方の包含から等号を得る。特に $KN$ は部分群の逆像なので部分群である（$N$ の正規性からも従う）。
:::
::::

::::exercise{level=advanced}
内部自己同型群 $\operatorname{Inn}(G)$ が巡回群ならば、$\operatorname{Inn}(G)$ は自明であることを示せ。

:::solution
:ref[prop-grp-inn-center] により $\operatorname{Inn}(G)\cong G/Z(G)$ なので、$G/Z(G)$ が巡回群である。:ref[prop-grp-center-quotient-cyclic] により $G$ は可換であり、$Z(G)=G$ となる。従って $G/Z(G)$、すなわち $\operatorname{Inn}(G)$ は自明である。例えば $\operatorname{Inn}(G)$ の位数は $2,3,5$ などの素数にはなりえない（素数位数の群は巡回群であるから）。
:::
::::

::::exercise{level=advanced}
二つの準同型 $f:G\to H$, $k:G\to K$ に対し、$G/(\Ker f\cap\Ker k)$ が $H\times K$ のある部分群と同型であることを示せ。

:::solution
$\psi:G\to H\times K$, $g\mapsto(f(g),k(g))$ は $\psi(gg')=(f(g)f(g'),k(g)k(g'))=\psi(g)\psi(g')$ により準同型である。$\psi(g)=(e_H,e_K)$ は $g\in\Ker f$ かつ $g\in\Ker k$ と同値なので、$\Ker\psi=\Ker f\cap\Ker k$ である。第一同型定理により $G/(\Ker f\cap\Ker k)\cong\Img\psi\le H\times K$ である。$f,k$ を自然射影 $G\to G/N$, $G\to G/M$ にとれば :ref[prop-grp-two-quotients] の前半になる。
:::
::::

::::exercise{level=advanced}
$S_3$ から $C_6$ への準同型をすべて求めよ。また、それぞれの準同型について第一同型定理が与える同型を書け。

:::solution
$r=(1\,2\,3)$, $s=(1\,2)$ とおくと、$S_3$ は $D_3$ と同一視でき（$srs^{-1}=(2\,1\,3)=r^{-1}$）、:ref[prop-grp-hom-dihedral] により準同型は $3A=\bar0$, $2B=\bar0$, $B+A-B=-A$ を満たす $(A,B)\in C_6^2$ に対応する。第三式は $2A=\bar0$ であり、$3A=\bar0$ と合わせて $A=3A-2A=\bar0$ を得る。$2B=\bar0$ より $B=\bar0$ または $\bar3$ である。

一つ目は自明な準同型で、核は $S_3$、第一同型定理は $S_3/S_3\cong\{\bar0\}$ を与える。二つ目は $r^is^j\mapsto\overline{3j}$ で、核は $\langle r\rangle=A_3$、像は $\{\bar0,\bar3\}$ であり、$S_3/\langle(1\,2\,3)\rangle\cong\langle\bar3\rangle\cong C_2$ を与える。これは符号 $\sgn$ を $\{\pm1\}\cong\{\bar0,\bar3\}$ と読み替えたものである。
:::
::::

::::exercise{level=advanced}
$f:G\to H$ を全射準同型、$K\trianglelefteq G$ とする。$f(K)\trianglelefteq H$ であり、$H/f(K)\cong G/(K\Ker f)$ であることを示せ。さらに $G=C_{24}$、$f$ を自然射影 $C_{24}\to C_8$、$K=\langle\bar6\rangle$ として両辺を具体的に求めよ。

:::solution
全射性と :ref[prop-grp-image-preimage] の 3 から $f(K)\trianglelefteq H$ である。合成 $\pi_{f(K)}\circ f:G\to H/f(K)$ は全射の合成なので全射であり、その核は $f^{-1}(f(K))$ である。前の演習により $f^{-1}(f(K))=K\Ker f$ なので、第一同型定理から $G/(K\Ker f)\cong H/f(K)$ である。

例では $\Ker f=\langle\bar8\rangle=\{\bar0,\bar8,\overline{16}\}$、$K=\{\bar0,\bar6,\overline{12},\overline{18}\}$ である。$K+\Ker f$ は $\bar8-\bar6=\bar2$ を含み、すべての元は偶数の類なので $K+\Ker f=\langle\bar2\rangle$ である。一方 $C_8$ で $f(K)=\langle\bar6\rangle=\{\bar0,\bar6,\bar4,\bar2\}=\langle\bar2\rangle$ である。両辺とも偶奇だけを記録する位数 $2$ の群であり、同型は $\bar1$ の類を $\bar1$ の類へ送る。
:::
::::

::::exercise{level=advanced}
$\GL_2(\F_2)\cong S_3$ を示せ。

:::hint
$\F_2^2$ の非零ベクトルは三つある。
:::

:::solution
$X=\F_2^2\setminus\{0\}=\{e_1,e_2,e_1+e_2\}$ とおく。$A\in\GL_2(\F_2)$ は単射な線形写像なので非零ベクトルを非零ベクトルに移し、$X$ 上の置換 $\rho(A)$ を定める。$(AB)x=A(Bx)$ より $\rho(AB)=\rho(A)\rho(B)$ なので、$\rho:\GL_2(\F_2)\to S_X$ は準同型である。$\rho(A)=\id$ なら $Ae_1=e_1$, $Ae_2=e_2$ より $A=I$ なので、核は自明で $\rho$ は単射である。:ref[prop-grp-gl-order] より $\abs{\GL_2(\F_2)}=(4-1)(4-2)=6=\abs{S_X}$ なので $\rho$ は全単射であり、$\GL_2(\F_2)\cong S_X\cong S_3$ である。
:::
::::

::::exercise{level=advanced}
1. $(\Q,+)$ から $(\Z,+)$ への準同型は零写像に限ることを示せ。
2. $\mu_\infty=\{z\in\C:\text{ある }n\ge1\text{ について }z^n=1\}$ は $\C^\times$ の部分群であり、$\Q/\Z\cong\mu_\infty$ であることを示せ。

:::solution
1. $f:\Q\to\Z$ を準同型とし、$x\in\Q$ とする。任意の $n\ge1$ について $f(x)=f(n\cdot\tfrac xn)=nf(\tfrac xn)$ なので、整数 $f(x)$ はすべての正整数で割り切れる。従って $f(x)=0$ である。

2. $z^n=1$, $w^m=1$ なら $(zw^{-1})^{nm}=1$ なので（$\C^\times$ は可換）、$\mu_\infty$ は部分群である。準同型 $\R\to S^1$, $t\mapsto e^{2\pi it}$ を $\Q$ に制限した $g:\Q\to\C^\times$ を考える。$t=a/b$（$b\ge1$）なら $g(t)^b=e^{2\pi ia}=1$ なので $g(\Q)\subseteq\mu_\infty$ である。逆に $z\in\mu_\infty$, $z^n=1$ とすると $\abs z^n=1$ より $\abs z=1$ なので $z=e^{2\pi it}$（$0\le t<1$）と書け、$e^{2\pi int}=1$ より $nt\in\Z$、従って $t\in\Q$ である。よって $g$ の像は $\mu_\infty$ であり、核は $\Q\cap\Z=\Z$ である（:ref[exm-grp-kernel-examples] の 3）。第一同型定理により $\Q/\Z\cong\mu_\infty$ である。特に $\Q/\Z$ のすべての元は有限位数である。
:::
::::

::::exercise{level=advanced}
$D_4$ の中心 $Z=\{e,r^2\}$ について、$D_4/Z\cong C_2\times C_2$ を示せ。さらに対応定理を用いて、$Z$ を含む $D_4$ の正規部分群をすべて求めよ。

:::solution
$\abs{D_4/Z}=8/2=4$ である。$D_4$ の任意の元 $g$ について $g^2\in Z$ である。実際、回転 $r^i$ なら $r^{2i}\in\{e,r^2\}$、鏡映なら $g^2=e$ である。従って $D_4/Z$ の非単位元はすべて位数 $2$ であり、:ref[prop-grp-order-four] により $D_4/Z\cong C_2\times C_2$ である。

$C_2\times C_2$ は可換なので部分群はすべて正規であり、自明群、位数 $2$ の三つ、全体の五つである。$D_4/Z$ の位数 $2$ の部分群は $\{Z,rZ\}$, $\{Z,sZ\}$, $\{Z,rsZ\}$ である。:ref[cor-grp-correspondence-normal] により、$Z$ を含む $D_4$ の正規部分群はその逆像

$$
Z,\quad \langle r\rangle=\{e,r,r^2,r^3\},\quad \{e,r^2,s,r^2s\},\quad \{e,r^2,rs,r^3s\},\quad D_4
$$

の五つである。位数 $4$ の三つはいずれも指数 $2$ であり、正規であることは :ref[prop-grp-index-two] とも整合する。
:::
::::

## 確認テスト

::::quiz
準同型 $f:G\to H$ が単射であるための必要十分条件はどれか。

- [ ] 像が可換である
- [x] 核が自明である
- [ ] 核が正規部分群である
- [ ] $f(e_G)=e_H$ である

:::explanation
正解は「核が自明である」。:ref[prop-grp-fiber] により $f(a)=f(b)$ は $a^{-1}b\in\Ker f$ と同値であり、各ファイバーは核の剰余類 $a\Ker f$ である。核が自明なら各ファイバーは一点なので単射であり、逆に単射なら $e_H$ の逆像は $\{e_G\}$ である。核はどんな準同型でも正規部分群であり、$f(e_G)=e_H$ もすべての準同型で成り立つので、これらは単射性の判定に使えない。像の可換性も単射性とは無関係である（例えば $C_2$ からの自明な準同型の像は可換だが単射でない）。
:::
::::

::::quiz{answer="12" tol="0"}
位数 $60$ の群 $G$ からの準同型 $f$ の核の位数が $5$ であるとき、像 $\Img f$ の位数はいくらか。

:::explanation
:ref[cor-grp-kernel-image-order] より $\abs G=\abs{\Ker f}\cdot\abs{\Img f}$ なので、$\abs{\Img f}=60/5=12$ である。第一同型定理の言葉では、$\Img f\cong G/\Ker f$ であり、その位数は指数 $[G:\Ker f]=12$ である。
:::
::::

::::quiz{answer="8" tol="0"}
巡回群 $C_{15}$ の自己同型群 $\Aut(C_{15})$ の位数はいくらか。

:::explanation
:ref[thm-grp-aut-cyclic] により $\Aut(C_{15})\cong(\Z/15\Z)^\times$ であり、位数は $15$ と互いに素な $1\le u\le15$ の個数 $\varphi(15)=\varphi(3)\varphi(5)=2\cdot4=8$ である。具体的には $u=1,2,4,7,8,11,13,14$ に対応する $x\mapsto ux$ が自己同型である。
:::
::::

::::quiz
$N\trianglelefteq G$ のとき、$G/N$ の部分群全体と一対一に対応するものはどれか。

- [ ] $G$ の正規部分群全体
- [x] $N$ を含む $G$ の部分群全体
- [ ] $N$ に含まれる $G$ の部分群全体
- [ ] $G$ の部分群全体

:::explanation
正解は「$N$ を含む $G$ の部分群全体」（:ref[thm-grp-correspondence]）。自然射影の逆像 $\pi_N^{-1}(L)$ は常に核 $N$ を含み、逆に $N\subseteq H$ なら $\pi_N^{-1}(H/N)=H$ となって元に戻る。$N$ を含まない部分群 $H$ では $\pi_N^{-1}(\pi_N(H))=HN\ne H$ となり、異なる部分群が同じ像をもちうる（$C_2\times C_2$ の例）。正規部分群に限ると、対応するのは $N$ を含む正規部分群と $G/N$ の正規部分群である。
:::
::::

::::quiz
$f:G\to H$ を任意の群準同型とする。常に成り立つものをすべて選べ。

- [x] $H$ の正規部分群の逆像は $G$ の正規部分群である
- [ ] $G$ の正規部分群の像は $H$ の正規部分群である
- [x] 像に属する任意の点のファイバーは核の左剰余類である
- [x] $a\in G$ が有限位数なら $\operatorname{ord}f(a)$ は $\operatorname{ord}a$ を割り切る

:::explanation
正解は一つ目、三つ目、四つ目である。逆像の正規性は $f(gxg^{-1})=f(g)f(x)f(g)^{-1}$ から従う（:ref[prop-grp-image-preimage]）。ファイバーについては、$f(x)=f(a)$ と $a^{-1}x\in\Ker f$ が同値なので $f(a)$ 上のファイバーは $a\Ker f$ である（:ref[prop-grp-fiber]）。位数については $f(a)^{\operatorname{ord}a}=f(e)=e$ による（:ref[prop-grp-hom-rules]）。二つ目は誤りで、$C_2$ を $S_3$ の部分群 $\langle(1\,2)\rangle$ に埋め込むと、定義域全体（正規）の像が $S_3$ で正規にならない。正規性が保証されるのは像 $f(G)$ の中だけである。
:::
::::

::::quiz{answer="6" tol="0"}
$S_4$ の部分群のうち、Klein の四元群 $V_4=\{e,(1\,2)(3\,4),(1\,3)(2\,4),(1\,4)(2\,3)\}$ を含むものはいくつあるか。

:::explanation
:ref[exm-grp-s4-klein-quotient] により全射準同型 $\Phi:S_4\to S_3$ で核が $V_4$ のものがあり、$S_4/V_4\cong S_3$ である。対応定理（:ref[thm-grp-correspondence]）により、$V_4$ を含む $S_4$ の部分群は $S_3$ の部分群と一対一に対応する。$S_3$ の部分群は自明群、互換が生成する三つ、$\langle(1\,2\,3)\rangle$、$S_3$ の六つなので、答えは $6$ である。具体的には $V_4$、位数 $8$ の三つ、$A_4$、$S_4$ である（:ref[exm-grp-correspondence-s4]）。
:::
::::

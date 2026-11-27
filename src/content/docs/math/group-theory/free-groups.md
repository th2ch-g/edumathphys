---
title: 自由群と群の表示
description: 簡約語から自由群を構成し、Schreier グラフでその部分群を調べる。さらに群の表示の普遍性と Tietze の定理を証明し、有限群と無限群の具体的な表示を求める。
notation:
  - symbol: 'X^{\pm}'
    meaning: 集合 X の文字とその形式的な逆の文字全体
  - symbol: 'F(X)'
    meaning: 集合 X を自由基底とする自由群
  - symbol: 'F_r'
    meaning: 階数 r の自由群
  - symbol: '\ell(w)'
    meaning: 語 w の文字数（長さ）
  - symbol: 'R(w)'
    meaning: 語 w を簡約して得られる簡約語
  - symbol: '\langle\!\langle R\rangle\!\rangle'
    meaning: 関係語集合 R の正規閉包
  - symbol: '\langle X\mid R\rangle'
    meaning: 生成元 X と関係語 R による群の表示
  - symbol: '\Z^{(X)}'
    meaning: X を基底とする自由 Abel 群（有限個を除き 0 の整数列全体）
---

## この章の目標

- 簡約語の一意性を証明し、簡約語全体が自由群をなすこと、その普遍性、階数の不変性を示す。
- 自由群の元の長さと Cayley グラフを用いて、自由群が捩れをもたないこと、可換な二元が同じ巡回群に入ることなどを調べる。
- 剰余類のグラフを用いて Nielsen–Schreier の定理（自由群の部分群は自由群）と Schreier の指数公式を証明する。
- 正規閉包で割ることにより群の表示を定義し、von Dyck の定理と「標準形による上からの評価」と「モデルへの全射による下からの評価」で表示群を同定する方法を身につける。
- Tietze 変換が群を変えないこと、および同型な群の有限表示が Tietze 変換で結ばれること（Tietze の定理）を証明し、具体的な表示の変形に使う。

二面体群 $D_n$ の回転と鏡映は $r^n=s^2=e,\ srs=r^{-1}$ を満たす（:ref[prop-grp-dihedral-form]）。この三式だけから積を計算できるのはなぜか。まず文字の間に関係式を課さずに自由群を作り、次に関係語を単位元にする最小の正規部分群で割る。こうして得た群と $D_n$ を比較すれば、この問いに答えられる。

本章では群の積 $\sigma\tau$ は右の $\tau$ を先に施す合成とし（:ref[def-grp-symmetric]）、交換子は $[x,y]=xyx^{-1}y^{-1}$ とする（:ref[def-grp-commutator]）。

## 自由群

### 語と簡約

語 $abb^{-1}$ と $a$ はどの群に代入しても同じ値になるが、$ab$ と $ba$ は $S_3$ に $a=(12),b=(23)$ を代入すれば $(123)\ne(132)$ となる。この違いを、隣接する逆の文字の対を消す操作で記述する。

:::definition[語と簡約語]{#def-grp-reduced-word}
集合 $X$ に対し、$X$ と交わらない集合 $X^{-1}=\{x^{-1}:x\in X\}$ を、$x\mapsto x^{-1}$ が全単射となるように用意する。$X^{\pm}=X\cup X^{-1}$ の元を **文字** という。文字 $s$ に対して $s^{-1}$ を、$s=x\in X$ なら $x^{-1}$、$s=x^{-1}$ なら $x$ と定める。したがって $(s^{-1})^{-1}=s$ であり、$s^{-1}\ne s$ である。

1. 文字の有限列 $w=s_1s_2\cdots s_k$ を $X$ 上の **語** といい、$k=\ell(w)$ をその長さという。長さ $0$ の語を空語といい $1$ と書く。二つの語の積は連結 $uv$ とし、語 $w=s_1\cdots s_k$ の逆を $w^{-1}=s_k^{-1}\cdots s_1^{-1}$ と定める。
2. 語の中で隣り合う二文字 $ss^{-1}$ を取り除く操作を **基本簡約**、その逆の操作（任意の位置に $ss^{-1}$ を挿入すること）を **基本挿入** という。基本簡約と基本挿入を有限回繰り返して $w$ から $w'$ に移れるとき $w\sim w'$ と書く。$\sim$ は同値関係である。
3. 基本簡約を一度も施せない語、すなわち $s_{i+1}=s_i^{-1}$ となる $i$ がない語を :term[簡約語]{en="reduced word" yomi="かんやくご"} という。
:::

例えば $aa^{-1}a$ は前の二文字を消しても後の二文字を消しても $a$ になる。長い語で消す順序を自由に選んでも結果が変わらないことを、左から一文字ずつ読む手続きで示す。

:::lemma[簡約の一意性]{#lem-grp-reduction}
1. 任意の語 $w$ から基本簡約を可能な限り繰り返すと有限回で簡約語に達し、その簡約語は基本簡約を施す順序によらない。これを $R(w)$ と書く。
2. 語 $w,w'$ について、$w\sim w'$ であることと $R(w)=R(w')$ は同値である。特に各同値類はちょうど一つの簡約語を含む。
:::

:::proof
**手続きの定義。** 簡約語 $r$ と文字 $s$ に対して
$$
r\ast s=\begin{cases}r' & (r=r's^{-1}\ \text{と書けるとき}),\\ rs & (\text{それ以外}),\end{cases}
$$
と定める。前者では簡約語の末尾を削り、後者では逆でない文字を追加するので、どちらも簡約語になる。語 $w=s_1\cdots s_k$ に対し
$$
\rho(w)=(\cdots((1\ast s_1)\ast s_2)\cdots)\ast s_k
$$
とおく。また簡約語 $r$ と語 $v=t_1\cdots t_m$ に対して $r\ast v=(\cdots(r\ast t_1)\cdots)\ast t_m$ と書く。定義から $\rho(w)$ は簡約語であり、任意の語 $u,v$ について
$$
\rho(uv)=\rho(u)\ast v \label{eq-grp-rho-concat}
$$
が成り立つ。

**ステップ 1。** $w$ は基本簡約だけで $\rho(w)$ になる。空語では自明である。$w=w's$ では、帰納法で $w'$ を $\rho(w')$ に簡約する。末尾が $s^{-1}$ ならもう一対を消し、そうでなければそのままで $\rho(w)$ を得る。

**ステップ 2。** $w$ が簡約語なら、その先頭から何文字読んでも消去は起きない。従って $\rho(w)=w$ である。

**ステップ 3。** $(r\ast s)\ast s^{-1}=r$ が成り立つ。$r=r's^{-1}$ なら $r\ast s=r'$ であり、$r'$ の末尾は $s$ でないので $s^{-1}$ を再び付ける。それ以外の場合は付けた $s$ を直ちに削る。

**ステップ 4：$\rho$ は基本簡約・基本挿入で変わらない。** 語 $u,v$ と文字 $s$ について、:eqref[eq-grp-rho-concat] とステップ 3 から
$$
\rho(uss^{-1}v)=\bigl((\rho(u)\ast s)\ast s^{-1}\bigr)\ast v=\rho(u)\ast v=\rho(uv)
$$
となる。従って $w\sim w'$ なら $\rho(w)=\rho(w')$ である。

**結論。** 基本簡約は長さを $2$ 減らすので必ず停止する。得た簡約語を $r$ とすると、ステップ 4 と 2 から $\rho(w)=\rho(r)=r$ となり、消去順序によらない。$R=\rho$ とすれば、ステップ 4 より同値な語の $R$ は等しく、逆に $R(w)=R(w')$ ならステップ 1 より $w\sim R(w)=R(w')\sim w'$ である。
:::

証明中の $\rho$ は、計算機で語を簡約するときの標準的な手順（スタックによる簡約）そのものである。例えば $w=ab^{-1}ba^{-1}b$ に対して
$$
1\xrightarrow{\ast a}a\xrightarrow{\ast b^{-1}}ab^{-1}\xrightarrow{\ast b}a\xrightarrow{\ast a^{-1}}1\xrightarrow{\ast b}b
$$
となり、$R(w)=b$ である。

:::remark
簡約語 $ab$ と $b^{-1}a^{-1}a^{-1}$ を連結すると、境目の $bb^{-1}$、次に $aa^{-1}$ が消えて $a^{-1}$ が残る。このように境目での消去は一対で終わるとは限らない。
:::

### 自由群の構成と普遍性

:::theorem[自由群の存在と普遍性]{#thm-grp-free-universal}
$X$ 上の簡約語全体を $F(X)$ とし、積を $u\cdot v=R(uv)$ で定める。

1. $F(X)$ はこの積について群をなす。単位元は空語、$u$ の逆元は逆語 $u^{-1}$ である。
2. 一文字の語を対応させる写像 $\iota:X\to F(X)$ は単射である。
3. 任意の群 $G$ と任意の写像 $f:X\to G$ に対し、$\hat f\circ\iota=f$ を満たす準同型 $\hat f:F(X)\to G$ がただ一つ存在する。

この群 $F(X)$ を $X$ 上の :term[自由群]{en="free group" yomi="じゆうぐん"} という。
:::

:::proof
1. **結合法則。** 語 $u,v,w$ について、ステップ 1（:ref[lem-grp-reduction] の証明）より $uv\sim R(uv)$ なので $uvw\sim R(uv)w$、同様に $uvw\sim uR(vw)$ である。:ref[lem-grp-reduction] の 2 より
$$
R\bigl(R(uv)w\bigr)=R(uvw)=R\bigl(uR(vw)\bigr)
$$
となり、$(u\cdot v)\cdot w=u\cdot(v\cdot w)$ を得る。
**単位元。** 簡約語 $u$ に対し $R(1u)=R(u1)=R(u)=u$ である。
**逆元。** $u=s_1\cdots s_k$ が簡約語なら $u^{-1}=s_k^{-1}\cdots s_1^{-1}$ も簡約語である。実際、$u^{-1}$ で隣り合う $s_{i+1}^{-1}s_i^{-1}$ が消去できるのは $s_i^{-1}=s_{i+1}$ のときで、それは $u$ の中の $s_is_{i+1}$ が消去できることを意味する。また $uu^{-1}=s_1\cdots s_ks_k^{-1}\cdots s_1^{-1}$ は中央から順に $s_ks_k^{-1}$, $s_{k-1}s_{k-1}^{-1}$, … と消えて空語になるので $u\cdot u^{-1}=1$、同様に $u^{-1}\cdot u=1$ である。

2. 異なる $x,y\in X$ の一文字の語 $x,y$ は異なる簡約語である。

3. **存在。** 写像 $\tilde f:X^{\pm}\to G$ を $\tilde f(x)=f(x)$, $\tilde f(x^{-1})=f(x)^{-1}$ で定め、語 $w=s_1\cdots s_k$ の値を $\mathrm{ev}(w)=\tilde f(s_1)\cdots\tilde f(s_k)$（空語では $e$）とする。任意の文字 $s$ について $\tilde f(s^{-1})=\tilde f(s)^{-1}$ だから、基本簡約・基本挿入で $\mathrm{ev}$ の値は変わらない。また $\mathrm{ev}(uv)=\mathrm{ev}(u)\mathrm{ev}(v)$ は定義から明らかである。そこで $\hat f(u)=\mathrm{ev}(u)$（$u\in F(X)$）とおくと、
$$
\hat f(u\cdot v)=\mathrm{ev}\bigl(R(uv)\bigr)=\mathrm{ev}(uv)=\mathrm{ev}(u)\,\mathrm{ev}(v)=\hat f(u)\hat f(v)
$$
となり、$\hat f$ は準同型で $\hat f(\iota(x))=f(x)$ を満たす。
**一意性。** 簡約語 $u=s_1\cdots s_k$ は $F(X)$ の中で一文字の語の積 $s_1\cdot s_2\cdots s_k$ に等しい。実際、簡約語の先頭部分 $s_1\cdots s_{i}$ に $s_{i+1}$ を連結しても消去は起こらない。さらに $F(X)$ の中で $x^{-1}$ は $\iota(x)$ の逆元である。従って $\iota(X)$ は $F(X)$ を生成し、:ref[prop-grp-hom-rules] により $\iota(X)$ 上で $f$ に一致する準同型は一つしかない。
:::

以後、$x\in X$ と $\iota(x)$ を同一視する。文字に $f$ の値を代入して $G$ で掛ける $\hat f$ を **代入準同型** と呼ぶ。

:::corollary[自由群の標準形]{#cor-grp-free-normal-form}
$F(X)$ の各元は
$$
x_1^{n_1}x_2^{n_2}\cdots x_m^{n_m}\qquad(m\ge0,\ x_i\in X,\ n_i\in\Z\setminus\{0\},\ x_i\ne x_{i+1})
$$
の形にただ一通りに書ける。
:::

:::proof
簡約語を同じ基底文字に由来する連続区間に切る。各区間では符号が一定なので一つの非零冪となり、区間を極大に取れば隣の基底文字は異なる。逆にこの形の式を一文字ずつ展開すると元の簡約語に戻る。従って :ref[lem-grp-reduction] より表し方は一意である。
:::

普遍性に現れる性質を、抽象的な群の言葉で定義しておく。

:::definition[自由基底]{#def-grp-free-basis}
群 $F$ とその部分集合 $S$ について、任意の群 $G$ と任意の写像 $f:S\to G$ に対し、$f$ を延長する準同型 $F\to G$ がただ一つ存在するとき、$S$ を $F$ の :term[自由基底]{en="free basis" yomi="じゆうきてい"} といい、$F$ を $S$ 上の自由群という。
:::

:ref[thm-grp-free-universal] は、$X$ が $F(X)$ の自由基底であることを述べている。

:::proposition[自由群の一意性]{#prop-grp-free-unique}
1. $S$ が群 $F$ の自由基底なら、$S$ は $F$ を生成する。
2. $S$ が $F$ の自由基底なら、包含写像 $S\hookrightarrow F$ を延長する準同型 $\varphi:F(S)\to F$ は同型である。特に自由群は自由基底の集合だけで同型を除いて決まる。
3. $X,Y$ の間に全単射 $\beta:X\to Y$ があれば、$\beta$ を延長する準同型 $F(X)\to F(Y)$ は同型である。
:::

:::proof
1. $H=\langle S\rangle$ とする。$S$ から $H$ への包含写像を延長する準同型 $\pi:F\to H$ がある。$H$ から $F$ への包含 $i$ との合成 $i\circ\pi:F\to F$ は、恒等写像 $\id_F$ と同じく $S$ 上で包含写像に一致する準同型である。自由基底の定義における一意性から $i\circ\pi=\id_F$ であり、$F=i(\pi(F))\subseteq H$ を得る。

2. $S$ が $F$ の自由基底であることから、$S\to F(S)$（文字を一文字の語に送る写像）を延長する準同型 $\psi:F\to F(S)$ がある。$\psi\circ\varphi:F(S)\to F(S)$ は $S$ 上で恒等的であり、:ref[thm-grp-free-universal] の一意性から $\id_{F(S)}$ に等しい。$\varphi\circ\psi:F\to F$ も $S$ 上で恒等的なので、自由基底の定義の一意性から $\id_F$ に等しい。よって $\varphi$ は同型である。

3. $\beta^{-1}$ を延長する準同型 $F(Y)\to F(X)$ をとれば、二つの合成は生成元上で恒等写像なので、:ref[prop-grp-hom-rules] より恒等写像である。
:::

:::corollary[任意の群は自由群の剰余群]{#cor-grp-quotient-of-free}
群 $G$ の生成系 $S$ をとると、包含写像 $S\hookrightarrow G$ を延長する準同型 $F(S)\to G$ は全射であり、$G\cong F(S)/N$（$N$ はその核）となる。特に、$S=G$ とすれば任意の群は自由群の剰余群である。
:::

:::proof
準同型の像は部分群であり、$S$ を含むので $\langle S\rangle=G$ を含む。後半は :ref[thm-grp-first-iso] による。
:::

:::example[小さい自由群]{#exm-grp-small-free}
1. $F(\varnothing)$ は空語だけからなる自明群である。
2. 一元集合 $\{a\}$ 上の簡約語は $a^n$（$n\in\Z$）であり、負の $n$ は $a^{-1}$ を $-n$ 個並べる。$a^ma^n=a^{m+n}$ だから $n\mapsto a^n$ は同型 $\Z\cong F(\{a\})$ を与える。
3. $F(\{a,b\})$ では $ab$ と $ba$ は異なる簡約語なので $ab\ne ba$ であり、$F(\{a,b\})$ は非可換である。$|X|\ge2$ なら同様に $F(X)$ は非可換である。
4. $f:\{a,b\}\to S_3$ を $f(a)=(12),\ f(b)=(23)$ とすると、代入準同型により $\hat f(aba)=(12)(23)(12)=(13)$、$\hat f(abab)=((12)(23))^2=(123)^2=(132)$ である。$\hat f$ は全射であるが、$\hat f(a^2)=e$ なので単射ではない。
:::

以下、$X=\{x_1,\ldots,x_r\}$ のとき $F(X)$ を $F(x_1,\ldots,x_r)$ とも書く。

### 自由群の階数

自由 Abel 群 $\Z^r$ の階数 $r$ が群だけから決まったように（:ref[thm-grp-fg-abelian]）、自由群の自由基底の大きさも群だけから決まる。

:::proposition[階数の不変性]{#prop-grp-free-rank}
$F(X)\cong F(Y)$ であるための必要十分条件は、$X$ と $Y$ が対等（同じ濃度）であることである。
:::

:::proof
十分性は :ref[prop-grp-free-unique] の 3 である。必要性を示す。$C_2=\{0,1\}$ を位数 $2$ の群とする。:ref[thm-grp-free-universal] の 3 により、準同型 $F(X)\to C_2$ は $X$ 上の値で決まり、どのような値の組も実現される。従って $\Hom(F(X),C_2)$ は写像全体 $\{0,1\}^X$ と全単射に対応する。

**$X$ が有限の場合。** $|\Hom(F(X),C_2)|=2^{|X|}$ は有限である。同型 $F(X)\cong F(Y)$ は全単射 $\Hom(F(Y),C_2)\to\Hom(F(X),C_2)$ を引き起こす（同型との合成）から、$\{0,1\}^Y$ も有限で $2^{|Y|}=2^{|X|}$ となり、$|Y|=|X|$ を得る。

**$X$ が無限の場合。** 上の議論より $Y$ も無限である。$X\subseteq F(X)$ であり、逆に $F(X)$ は $\bigcup_{k\ge0}(X^\pm)^k$ に含まれる。選択公理のもとでは無限集合とその有限直積（正の個数の因子）は対等であり、可算個のそのような集合の和も元と同じ濃度をもつ（[濃度の演算](/math/sets-logic/cardinal-arithmetic/)）。よって $|X|=|F(X)|=|F(Y)|=|Y|$ である。
:::

:::definition[自由群の階数]{#def-grp-free-rank}
自由群 $F$ の自由基底の濃度を $F$ の :term[階数]{en="rank of a free group" yomi="かいすう"} といい、$\rank F$ と書く。:ref[prop-grp-free-unique] と :ref[prop-grp-free-rank] により、これは自由基底の選び方によらない。階数 $r$ の自由群を $F_r$ と書く。
:::

自由基底は一つの群に対して一通りではない。

:::example[自由基底の取り替え]{#exm-grp-free-basis-change}
$F=F(a,b)$ で、代入 $a\mapsto ab,\ b\mapsto b$ が定める準同型を $\alpha$、$a\mapsto ab^{-1},\ b\mapsto b$ が定める準同型を $\beta$ とする。$\beta(\alpha(a))=\beta(a)\beta(b)=ab^{-1}b=a$、$\alpha(\beta(a))=abb^{-1}=a$ であり、両者は $b$ を固定する。生成元上で恒等的な準同型は恒等写像だから、$\alpha$ は自己同型で $\alpha^{-1}=\beta$ である。

一般に $\alpha$ が $F$ の自己同型で $S$ が自由基底なら、$\alpha(S)$ も自由基底である。実際、写像 $f:\alpha(S)\to G$ が与えられたとき、$f\circ\alpha|_S$ を延長する準同型 $g:F\to G$ をとれば、$g\circ\alpha^{-1}$ は $f$ を延長する。延長が二つあれば、それらに $\alpha$ を合成したものは $S$ 上で一致するので等しい。従って $\{ab,b\}$ も $F(a,b)$ の自由基底である。このように、自由基底の元を別の元との積で置き換える操作を **Nielsen 変換** という。
:::

自由群を可換化すると自由 Abel 群になる。$X$ を基底とする自由 Abel 群を $\Z^{(X)}=\bigoplus_{x\in X}\Z$ と書き、$x$ に対応する基底ベクトルを $e_x$ とする。

:::example[自由群の可換化と指数和]{#exm-grp-free-abelianization}
:ref[thm-grp-free-universal] により $x\mapsto e_x$ を延長する準同型 $\varepsilon:F(X)\to\Z^{(X)}$ がある。$\varepsilon(w)$ の $x$ 成分は、語 $w$ に現れる $x$ の個数から $x^{-1}$ の個数を引いた **指数和** である。$\Z^{(X)}$ は可換なので、:ref[prop-grp-abelianization] により $\varepsilon$ は可換化 $F(X)_{\mathrm{ab}}$ を経由する準同型 $\bar\varepsilon$ を引き起こす。逆向きに $\kappa:\Z^{(X)}\to F(X)_{\mathrm{ab}}$ を $\kappa\bigl(\sum_xn_xe_x\bigr)=\prod_x\bar x^{\,n_x}$（有限積、$\bar x$ は $x$ の像）で定めると、可換性によりこれは順序によらず準同型である。$\bar\varepsilon\circ\kappa$ と $\kappa\circ\bar\varepsilon$ は生成元 $e_x$, $\bar x$ を固定するので恒等写像であり、
$$
F(X)_{\mathrm{ab}}\cong\Z^{(X)}
$$
を得る。有限の $X$ では :ref[thm-grp-fg-abelian] からも階数の不変性が従う。例えば $aba^{-1}b^{-1}$ と $a^2ba^{-2}b^{-1}$ は両座標の指数和が $0$ なので交換子群に属する。
:::

### 巡回簡約語と捩れ

$aba^{-1}$ は簡約語だが、冪を取ると $(aba^{-1})^k=ab^ka^{-1}$ となる。一般の語でも冪の長さを調べるには、両端で互いに逆になる文字を先に取り分ける。

:::definition[巡回簡約語]{#def-grp-cyclically-reduced}
簡約語 $w=s_1\cdots s_k$ が、$k\le1$ であるか $s_k\ne s_1^{-1}$ を満たすとき、**巡回簡約語** という。すなわち、$w$ を円周上に並べても消去箇所がない語である。
:::

:::lemma[巡回簡約語への分解]{#lem-grp-cyclic-decomposition}
空でない任意の簡約語 $w$ は、簡約語 $u$ と空でない巡回簡約語 $v$ を用いて、文字列として $w=uvu^{-1}$（連結しても消去が起こらない形）と書ける。このとき $\ell(w)=2\ell(u)+\ell(v)$ である。
:::

:::proof
$\ell(w)$ に関する帰納法で示す。$w$ が巡回簡約語なら $u=1,\ v=w$ とすればよい。そうでなければ $\ell(w)\ge2$ で、$w=sw's^{-1}$（$s$ は文字）と書ける。$w'$ は簡約語 $w$ の中間部分なので簡約語であり、空ではない。実際、$w'$ が空なら $w=ss^{-1}$ となり簡約性に反する。帰納法の仮定から $w'=u'vu'^{-1}$ と書けるので、$u=su'$ とおけば文字列として $w=uvu^{-1}$ である。
:::

:::proposition[自由群は捩れをもたない]{#prop-grp-free-torsion-free}
$w\ne1$ を $F(X)$ の元とし、:ref[lem-grp-cyclic-decomposition] のように $w=uvu^{-1}$ と書く。任意の $k\ge1$ について $w^k$ の簡約語は文字列 $uv^ku^{-1}$ であり、
$$
\ell(w^k)=2\ell(u)+k\,\ell(v)
$$
である。特に $w^k\ne1$ であり、$F(X)$ の単位元以外の元はすべて無限位数である。
:::

:::proof
$w^k$ を連結すると因子の間の $u^{-1}u$ が消え、$uv^ku^{-1}$ が残る。両端の $u,v$ および $v,u^{-1}$ の境目は元の簡約語と同じである。さらに $v$ は巡回簡約されているので $v$ と $v$ の境目でも消えない。従って表示した文字数が簡約後の長さであり、$\ell(v)\ge1$ より $w^k\ne1$ である。負の冪についても逆元を取れば同じ結論を得る。
:::

$w=ab^2a^{-1}$ なら $w^3=ab^6a^{-1}$ の長さは $8$、巡回簡約語 $w=aba^{-1}b^{-1}$ でも $w^2$ の長さは $8$ となる。これに対して $\Z/n\Z$ や $S_n$ では、付加した関係式によって有限位数の非単位元が生じうる。

### Cayley グラフ

群 $G$ と生成系 $S$ に対し、$G$ の各元を頂点とし、各 $g\in G$ と $s\in S$ について $g$ から $gs$ へラベル $s$ の有向辺を引いたものを $(G,S)$ の **Cayley グラフ** という。辺を逆向きにたどるときは $s^{-1}$ を読むと約束すると、頂点 $g$ から出発して語 $w=s_1\cdots s_k$（$s_i\in S^{\pm}$）を読む道がただ一つ定まり、その道は $g,\ gs_1,\ gs_1s_2,\ \ldots$ を通って $g\bar w$ に到着する。ここで $\bar w$ は $w$ を $G$ で評価した元である。従って、頂点 $g$ から出る閉じた道で読む語は、ちょうど $G$ で単位元に評価される語である。

道が **後戻りしない** とは、ある辺を通った直後に同じ辺を逆向きに通ることがないことをいう。後戻りしない閉じた道（長さ $\ge1$）を **閉路** と呼び、連結で閉路のないグラフを **木** という。

:::proposition[自由群の Cayley グラフは木]{#prop-grp-cayley-tree}
$G=F(X)$、$S=X$ の Cayley グラフ $\Gamma$ について次が成り立つ。

1. $\Gamma$ の二頂点の間には高々一本の辺しかなく、辺の端点は相異なる。
2. 頂点 $g$ から出る道が後戻りしないことと、その道で読む語が簡約語であることは同値である。
3. $\Gamma$ は木である。さらに任意の二頂点 $g,h$ を結ぶ後戻りしない道はただ一つであり、それは $g^{-1}h$ の簡約語を読む道である。
4. グラフの距離（最短の道の辺の数）は $d(g,h)=\ell(g^{-1}h)$ である。
:::

:::proof
1. 辺の端点は $g$ と $gx$（$x\in X$）で、$x\ne1$ だから相異なる。頂点 $g,h$ を結ぶ辺が二本あるとすると、$h=gx$ と $h=gy$、または $h=gx$ と $g=hy$（$x,y\in X$）のいずれかの形に書ける。前者なら $x=y$ で同じ辺である。後者なら $xy=1$ となるが、$xy$ は長さ $2$ の簡約語であり（$y\ne x^{-1}$ は文字として明らか）、空語と異なるので矛盾である。

2. 1 より、隣り合う頂点の間の辺は一本しかない。従ってある辺を通った直後に同じ辺を逆向きに通ることは、読む語に隣り合う $ss^{-1}$ が現れることと同値である。

3. $g$ から $g'$ へは $g^{-1}g'$ の簡約語を読めば到達するので $\Gamma$ は連結である。$g$ から出る閉路があれば、2 よりそれは長さ $\ge1$ の簡約語 $w$ を読み、$g\bar w=g$、すなわち $F(X)$ で $w=1$ となる。これは異なる簡約語が異なる元を表すこと（:ref[thm-grp-free-universal]）に反する。$g$ から $h$ への後戻りしない道が読む語 $w$ は簡約語で $\bar w=g^{-1}h$ を満たすから、$w$ は $g^{-1}h$ の簡約語に等しく、道はこの語で一意に決まる。

4. $g$ から $h$ への任意の道が読む語 $w$ は $g^{-1}h$ を表す。:ref[lem-grp-reduction] より $w$ は基本簡約で $g^{-1}h$ の簡約語 $R(w)$ になるので $\ell(w)\ge\ell(R(w))=\ell(g^{-1}h)$ である。一方 $R(w)$ を読む道の長さはちょうど $\ell(g^{-1}h)$ である。
:::

例えば $F(a,b)$ で $g=ab,\ h=abb$ なら $g^{-1}h=b$ なので距離は $1$ である。$g=ab,\ h=ba$ なら $g^{-1}h=b^{-1}a^{-1}ba$ は既に簡約されており、距離は $4$ である。これに対し $\Z^2=\langle a,b\rangle$ の Cayley グラフは平面の格子であり、$a,b,a^{-1},b^{-1}$ を順に読む一辺 $1$ の正方形が閉路をなす。この閉路が関係式 $aba^{-1}b^{-1}=1$ に対応している。一般に $(G,S)$ の Cayley グラフの閉じた道は、自然な全射 $F(S)\to G$ の核に属する語と対応する。本章の後半で扱う関係式は、この核を記述する手段である。

### 自由群の部分群

自由 Abel 群の部分群の階数は元の階数以下である（:ref[lem-grp-lattice]）。自由群では $F(a,b)$ の部分群が階数 $3$ にもなる（:ref[exm-grp-index-two-subgroup]）。その自由基底を求めるため、右剰余類を頂点として、部分群の元を基点に戻る道で表す。

:::definition[Schreier グラフ]{#def-grp-schreier-graph}
$F=F(X)$ とし、$H\le F$ とする。右剰余類全体 $H\backslash F=\{Hg:g\in F\}$ を頂点の集合とし、各頂点 $v=Hg$ と各 $x\in X$ について、$v$ から $vx=Hgx$ へラベル $x$ の有向辺 $(v,x)$ を引く。$vx$ は代表元 $g$ の選び方によらない（$Hg=Hg'$ なら $Hgx=Hg'x$）。このグラフを $H$ の **Schreier グラフ** といい、頂点 $H$ を **基点** と呼ぶ。辺 $(v,x)$ を **正の辺** と呼び、Cayley グラフと同じく逆向きにたどるときは $x^{-1}$ を読む。
:::

頂点 $v$ から語 $w=s_1\cdots s_k$ を読む道は一通りに定まる。実際、$s=x$ なら現在の頂点 $u$ から正の辺 $(u,x)$ を進み、$s=x^{-1}$ なら正の辺 $(ux^{-1},x)$（これは $ux^{-1}$ から $u$ への辺である）を逆向きに進む。どちらの場合も到着点は $us$ である。従って $v$ から $w$ を読む道は $v\bar w$ に到着する。特に次が成り立つ。

:::lemma[閉じた道と部分群の元]{#lem-grp-schreier-paths}
基点から語 $w$ を読む道は頂点 $H\bar w$ に到着する。従ってその道が閉じていることと $\bar w\in H$ は同値である。また Schreier グラフは連結である。
:::

:::proof
前半は直前の説明による。$H\bar w=H$ と $\bar w\in H$ は同値である。任意の頂点 $Hg$ には、基点から $g$ の簡約語を読めば到達する。
:::

次に、各頂点から基点へ戻る「決まった道」を選ぶ。基点からの距離（最短の道の辺の数）を $\delta(v)$ と書く。

:::lemma[親の辺の選択]{#lem-grp-parent-edges}
基点以外の各頂点 $v$ に対し、$v$ と、$\delta(p(v))=\delta(v)-1$ を満たす頂点 $p(v)$ とを結ぶ正の辺 $e_v$ を一本ずつ選ぶ（これを $v$ の **親の辺** と呼ぶ）。$T=\{e_v\}$ とおくと次が成り立つ。

1. 各頂点 $v$ から $p$ を $\delta(v)$ 回繰り返すと基点に達する。
2. $v\mapsto e_v$ は単射である。
3. 語 $t_v$ を、$t_H=1$ とし、$e_v=(p(v),x)$ なら $t_v=t_{p(v)}x$、$e_v=(v,x)$ なら $t_v=t_{p(v)}x^{-1}$ として $\delta(v)$ に関して帰納的に定める。このとき基点から $t_v$ を読む道は $T$ の辺だけを通って $v$ に到着する。特に $Ht_v=v$ である。
:::

:::proof
親の辺は存在する。実際、$\delta(v)=n\ge1$ なら長さ $n$ の道で基点から $v$ に到達でき、その最後の辺の手前の頂点 $u$ は $\delta(u)\le n-1$ を満たし、$v$ と隣り合うので $\delta(u)\ge n-1$ でもある。各頂点について一本ずつ選ぶには一般に選択公理を用いる（$X$ が可算なら $F$ も可算で、語を一列に並べて最初のものを選べばよい）。

1. $\delta$ が一回ごとに $1$ ずつ減るので $\delta(v)$ 回で $0$、すなわち基点に達する。

2. $e_v=e_w$（$v\ne w$）とすると、この辺は $v,p(v)$ を結び、かつ $w,p(w)$ を結ぶので、$w=p(v)$ かつ $v=p(w)$ である。すると $\delta(w)=\delta(v)-1$ かつ $\delta(v)=\delta(w)-1$ となり矛盾する。

3. $\delta(v)$ に関する帰納法による。基点から $t_{p(v)}$ を読む道は $T$ の辺だけを通って $p(v)$ に着く。そこから $e_v=(p(v),x)$ なら $x$ を読んで正の辺 $e_v$ を進み、$e_v=(v,x)$ なら $x^{-1}$ を読んで $e_v$ を逆向きに進む。いずれも $T$ の辺で、到着点は $v$ である。
:::

親をたどると必ず基点に着くので $T$ は連結である。閉路があれば最短のものを取る。その始点でも最初と最後の辺は異なる。同じ辺なら両端の往復を取り除いて、より短い閉路が残るからである。そこでこの閉路上で $\delta$ が最大の頂点を見る。そこから両側へ出る二辺はともに距離を減らすが、そのような親の辺は一本しかない。従って $T$ は全頂点を結ぶ木（全域木）である。木の基点から $v$ までの語 $t_v$ は、先頭部分も代表系に属する。このような右剰余類の代表系を **Schreier 代表系** という。

:::theorem[Nielsen–Schreier の定理]{#thm-grp-free-subgroup}
$H\le F(X)$ とし、Schreier グラフと親の辺の集合 $T$ を上のようにとる。$T$ に属さない正の辺全体を $E'$ とし、$e=(v,x)\in E'$（$v$ から $w=vx$ への辺）に対して
$$
s_e=t_v\,x\,t_w^{-1}\in F(X) \label{eq-grp-schreier-generator}
$$
とおく。このとき $H$ は自由群であり、$\{s_e:e\in E'\}$ は $H$ の自由基底である（$e\mapsto s_e$ は単射である）。特に、自由群の部分群は自由群である。さらに $|X|=r<\infty$、$[F(X):H]=d<\infty$ なら
$$
\rank H=1+d(r-1) \label{eq-grp-schreier-index}
$$
である（**Schreier の指数公式**）。
:::

:::proof
:eqref[eq-grp-schreier-generator] の式で、$T$ に属する正の辺 $e$ に対しても $s_e$ を定めておく。

**ステップ 1：$s_e\in H$、かつ $e\in T$ なら $s_e=1$。** 基点から $s_e$ を読む道は、$t_v$ で $v$ へ行き、辺 $e$ を進んで $w$ に着き、$t_w^{-1}$ で $t_w$ の道を逆にたどって基点に戻る。閉じた道なので :ref[lem-grp-schreier-paths] より $s_e\in H$ である。$e\in T$ なら、$e$ は $w$ の親の辺（$p(w)=v$）か $v$ の親の辺（$p(v)=w$）である。前者なら $t_w=t_vx$ で $s_e=t_vxx^{-1}t_v^{-1}=1$、後者なら $t_v=t_wx^{-1}$ で $s_e=t_wx^{-1}xt_w^{-1}=1$ である。

**ステップ 2：$\{s_e:e\in E'\}$ は $H$ を生成する。** $h\in H$ を表す語 $w=s_1\cdots s_k$ をとり、基点からこの語を読む道の頂点を $v_0=H,v_1,\ldots,v_k=H$ とする（$v_i=Hs_1\cdots s_i$）。$t_{v_0}=t_{v_k}=1$ だから、$F(X)$ の中で
$$
h=\prod_{i=1}^{k}\bigl(t_{v_{i-1}}\,s_i\,t_{v_i}^{-1}\bigr)
$$
が成り立つ（隣り合う因子の間で $t_{v_i}^{-1}t_{v_i}$ が消える）。第 $i$ 因子は、$s_i=x$ で第 $i$ 歩が正の辺 $e=(v_{i-1},x)$ を進むなら $s_e$ に等しい。$s_i=x^{-1}$ で第 $i$ 歩が正の辺 $e=(v_i,x)$ を逆向きに進むなら $t_{v_{i-1}}x^{-1}t_{v_i}^{-1}=(t_{v_i}xt_{v_{i-1}}^{-1})^{-1}=s_e^{-1}$ である。ステップ 1 より $T$ の辺に対応する因子は $1$ なので、$h$ は $s_e^{\pm1}$（$e\in E'$）の積である。

**ステップ 3：記録写像。** $E'$ の元に対応する記号 $z_e$ を用意し、自由群 $F(E')$（基底 $\{z_e\}$）を考える。道 $\pi$ が辺 $e_1,\ldots,e_k$ を向き $\epsilon_1,\ldots,\epsilon_k\in\{\pm1\}$（正の向きなら $+1$）で通るとき、その **記録** を
$$
\operatorname{rec}(\pi)=\prod_{i:\ e_i\in E'}z_{e_i}^{\epsilon_i}\in F(E')
$$
（$i$ の小さい順の積、$T$ の辺は読み飛ばす）と定める。$h\in H$ に対し、$h$ を表す語 $w$ を基点から読む道の記録を $\psi(h)$ とおく。

- **代表の語によらないこと。** $h$ を表す二つの語は :ref[lem-grp-reduction] より基本挿入・基本簡約で移り合う。語の途中に $ss^{-1}$ を挿入すると、道はその位置である辺を通ってすぐ同じ辺を戻り、同じ頂点から元の道を続ける。記録には $z_e^{\epsilon}z_e^{-\epsilon}$ が挿入されるか（$e\in E'$ のとき）、何も変わらない（$e\in T$ のとき）。いずれにせよ $F(E')$ の元としての記録は変わらない。従って $\psi:H\to F(E')$ は well-defined である。
- **準同型であること。** $h,h'\in H$ を表す語 $w,w'$ について、$ww'$ は $hh'$ を表し、それを読む道は $w$ の閉じた道に $w'$ の閉じた道を続けたものである。従って $\psi(hh')=\psi(h)\psi(h')$ である。
- **$\psi(s_e)=z_e$（$e\in E'$）。** $s_e=t_vxt_w^{-1}$ を読む道は、$T$ の辺だけを通って $v$ へ行き、$e$ を正の向きに通り、$T$ の辺だけを通って基点に戻る（:ref[lem-grp-parent-edges] の 3）。記録は $z_e$ である。

**ステップ 4：同型。** :ref[thm-grp-free-universal] により $z_e\mapsto s_e$ を延長する準同型 $\varphi:F(E')\to H$ があり、ステップ 2 より全射である。$\psi\circ\varphi:F(E')\to F(E')$ は各 $z_e$ を固定するので恒等写像である（:ref[prop-grp-hom-rules]）。よって $\varphi$ は単射でもあり同型である。特に $e\ne e'$ なら $s_e=\varphi(z_e)\ne\varphi(z_{e'})=s_{e'}$ であり、$\{s_e\}$ は $H$ の自由基底である（:ref[prop-grp-free-unique] の 2 の議論を $\varphi$ に適用すればよい）。

**ステップ 5：指数公式。** 頂点は $d$ 個、各頂点から各 $x\in X$ について正の辺がちょうど一本出るので正の辺は $dr$ 本ある。:ref[lem-grp-parent-edges] の 2 により $T$ の辺は基点以外の頂点と一対一に対応し、$d-1$ 本である。従って $|E'|=dr-(d-1)=1+d(r-1)$ である。
:::

最短路から選んだ木以外でも同じ計算ができる。任意の全域木には基点から各頂点への単純な道がただ一つあり、それを読む語を $t_v$ とすれば木の辺で $s_e=1$、それ以外の辺で記録は $z_e$ となる。従って上の証明はそのまま通る。位相的には全域木を一点に縮めたグラフの基本群を計算している。被覆空間を用いる証明は :cite[hatcher-at] を参照。

:::example[指数 2 の部分群]{#exm-grp-index-two-subgroup}
$F=F(a,b)$ から $C_2=\{\bar0,\bar1\}$ への代入準同型 $a\mapsto\bar1,\ b\mapsto\bar0$ の核を $H$ とする。$H$ は「$a$ の指数和が偶数の元」の全体で、指数は $2$ である。剰余類 $H,Ha$ を $0,1$ と書くと、Schreier グラフの正の辺は

- $(0,a):0\to1$、$(1,a):1\to0$、
- $(0,b):0\to0$、$(1,b):1\to1$（ループ）

の四本である。頂点 $1$ の親の辺として $(0,a)$ を選ぶと $t_0=1,\ t_1=a$ である。$E'=\{(1,a),(0,b),(1,b)\}$ に対応する自由基底は
$$
s_{(1,a)}=t_1\,a\,t_0^{-1}=a^2,\qquad s_{(0,b)}=b,\qquad s_{(1,b)}=t_1\,b\,t_1^{-1}=aba^{-1}
$$
であり、$\rank H=3=1+2(2-1)$ は指数公式と一致する。

$h=abab^{-1}\in H$ をこの基底で書くには、ステップ 2 の手順に従う。基点から $a,b,a,b^{-1}$ を読むと、辺 $(0,a)$（木）、$(1,b)$、$(1,a)$、$(0,b)$ の逆向きを順に通る。従って
$$
abab^{-1}=s_{(1,b)}\,s_{(1,a)}\,s_{(0,b)}^{-1}=(aba^{-1})(a^2)(b^{-1})
$$
である。右辺を簡約すると確かに $abab^{-1}$ に戻る。
:::

:::example[交換子群は無限階数の自由群]{#exm-grp-commutator-subgroup-free}
$F=F(a,b)$ の交換子群 $F'$ を考える。:ref[exm-grp-free-abelianization] により $F/F'\cong\Z^2$ であり、剰余類 $F'g$ は $g$ の指数和 $(m,n)\in\Z^2$ で区別される。従って Schreier グラフの頂点は格子点 $(m,n)$ であり、正の辺は $a$ による $(m,n)\to(m+1,n)$ と $b$ による $(m,n)\to(m,n+1)$ である。基点からの距離は $\delta(m,n)=|m|+|n|$ である。

親の辺を、$n\ne0$ なら縦の辺（$(m,n)$ と $(m,n\mp1)$ を結ぶ $b$ の辺）、$n=0,\ m\ne0$ なら横の辺（$(m,0)$ と $(m\mp1,0)$ を結ぶ $a$ の辺）と選ぶ。すると $T$ は「横軸全体とすべての縦の辺」からなり、$t_{(m,n)}=a^mb^n$ である。$T$ に属さない正の辺は $n\ne0$ の行にある横の辺 $(m,n)\to(m+1,n)$ であり、
$$
s_{(m,n)}=a^mb^n\,a\,b^{-n}a^{-(m+1)}\qquad(m\in\Z,\ n\in\Z\setminus\{0\})
$$
が $F'$ の自由基底をなす。特に $F'$ は可算無限階数の自由群である。例えば交換子 $[a,b]=aba^{-1}b^{-1}$ を読む道は $(0,0)\to(1,0)\to(1,1)\to(0,1)\to(0,0)$ で、第三歩だけが木にない辺 $(0,1)\to(1,1)$ を逆向きに通る。よって $[a,b]=s_{(0,1)}^{-1}=(bab^{-1}a^{-1})^{-1}$ である。

従って、有限生成群の部分群は無限個の生成元を必要とすることもある。
:::

:::corollary[自由群の中の自由部分群]{#cor-grp-free-subgroups-ranks}
$B$ が自由群 $F$ の自由基底で $B'\subseteq B$ なら、$\langle B'\rangle$ は $B'$ を自由基底とする自由群である。特に $F_2$ は、任意の有限階数および可算無限階数の自由群を部分群として含む。
:::

:::proof
:ref[prop-grp-free-unique] により $F=F(B)$ としてよい。$B'$ 上の簡約語は $B$ 上でも簡約語であり、異なる簡約語は異なる元を表す。従って包含 $B'\hookrightarrow F(B)$ を延長する準同型 $F(B')\to F(B)$ は単射で、像は $\langle B'\rangle$ である。後半は :ref[exm-grp-commutator-subgroup-free] の可算無限個の自由基底から、有限個または可算無限個を選べばよい。
:::

Nielsen–Schreier の定理を使うと、自由群の可換性について強い結論が得られる。

:::corollary[自由群の可換な元]{#cor-grp-free-commuting}
$F=F(X)$ とする。

1. $g,h\in F$ が $gh=hg$ を満たせば、ある $c\in F$ があって $g,h\in\langle c\rangle$ となる。
2. $F$ の元 $x\in X$ は真の冪ではない。すなわち $x=c^k$ なら $k=\pm1$ である。
3. $|X|\ge2$ なら $F$ の中心は自明である。
4. $g\ne1$ の中心化群 $C_F(g)$ は無限巡回群である。
:::

:::proof
1. $K=\langle g,h\rangle$ は可換であり、:ref[thm-grp-free-subgroup] により自由基底 $B$ をもつ。$B$ が異なる二元 $p,q$ を含めば、$pq$ と $qp$ は $B$ 上の異なる簡約語なので $pq\ne qp$ となり可換性に反する。従って $|B|\le1$ で $K$ は巡回群である。

2. $x=c^k$ とする。$k=0$ はありえない。$k<0$ なら $c$ を $c^{-1}$ に替えて $k>0$ としてよい。$c=uvu^{-1}$ と :ref[lem-grp-cyclic-decomposition] のように書けば、:ref[prop-grp-free-torsion-free] より $1=\ell(x)=2\ell(u)+k\ell(v)\ge k$ なので $k=1$ である。

3. $z$ を中心の元とし、$X$ の異なる二元 $x,y$ をとる。1 より $\langle z,x\rangle=\langle c\rangle$ と書け、$x=c^k$ だから 2 より $c=x^{\pm1}$、従って $z\in\langle x\rangle$ である。同様に $z\in\langle y\rangle$ であり、$z=x^m=y^n$ と書ける。標準形の一意性（:ref[cor-grp-free-normal-form]）から $m=n=0$、すなわち $z=1$ である。

4. $C=C_F(g)$ は :ref[thm-grp-free-subgroup] により自由群で、その自由基底を $B$ とする。$g$ は $C$ の中心に属し $g\ne1$ だから、3 を自由群 $C\cong F(B)$ に適用すると $|B|\le1$ である。$g\in C$ より $C\ne1$ なので $|B|=1$、すなわち $C$ は無限巡回群である（:ref[exm-grp-small-free]）。
:::

例えば $C_{F(a,b)}(a)=\langle a\rangle$ であり、共役を取ると $C_{F(a,b)}(aba^{-1})=\langle aba^{-1}\rangle$ である。

## 生成元と関係式

正方形の四分の一回転 $r$ について $r^4=e$ を課したいとする。$F(r,s)$ の中で $r^4$ を単位元にする剰余群を作るには、$r^4$ を含む **正規** 部分群で割らなければならない（:ref[thm-grp-quotient]）。$r^4$ が単位元なら、その共役 $ur^4u^{-1}$ も、それらの積や逆元も単位元になる。従って割るべきものは、関係語のすべての共役から生成される部分群である。

:::definition[正規閉包と群の表示]{#def-grp-presentation}
群 $G$ の部分集合 $R$ を含む正規部分群すべての共通部分を $R$ の :term[正規閉包]{en="normal closure" yomi="せいきへいほう"} といい、$\langle\!\langle R\rangle\!\rangle$（群を明示するときは $\langle\!\langle R\rangle\!\rangle_G$）と書く。これは $R$ を含む最小の正規部分群である。

集合 $X$ と $R\subseteq F(X)$ に対し、剰余群
$$
\langle X\mid R\rangle=F(X)/\langle\!\langle R\rangle\!\rangle
$$
を、:term[生成元]{en="generator" yomi="せいせいげん"} $X$ と :term[関係語]{en="relator" yomi="かんけいご"} $R$ による :term[群の表示]{en="group presentation" yomi="ぐんのひょうじ"} という。群 $G$ と同型 $G\cong\langle X\mid R\rangle$ の組を $G$ の表示ともいう。語 $u,v$ について **関係式** $u=v$ を課すとは、$uv^{-1}$ を関係語に加えることを意味する。$X,R$ がともに有限集合である表示を **有限表示** という。
:::

表示を書くときは $\langle x_1,\ldots,x_n\mid r_1,\ldots,r_m\rangle$ のように集合の括弧を省き、関係語の代わりに関係式を書くことも多い。例えば $\langle r,s\mid r^n=s^2=1,\ srs=r^{-1}\rangle$ は $\langle r,s\mid r^n,s^2,srsr\rangle$ と同じ群である（$srs(r^{-1})^{-1}=srsr$）。自然な射影を $\pi:F(X)\to\langle X\mid R\rangle$ とし、$\pi(x)$ も単に $x$ と書く。ただし異なる $x,y\in X$ の像が一致することもある（例えば $\langle a,b\mid a=b\rangle$）ので注意する。

:::proposition[正規閉包の元]{#prop-grp-normal-closure}
群 $G$ と $R\subseteq G$ について、$\langle\!\langle R\rangle\!\rangle$ は
$$
\prod_{i=1}^{k}u_i\,r_i^{\epsilon_i}\,u_i^{-1}\qquad(k\ge0,\ u_i\in G,\ r_i\in R,\ \epsilon_i=\pm1) \label{eq-grp-normal-closure}
$$
の形の元全体に等しい。
:::

:::proof
:eqref[eq-grp-normal-closure] の形の元全体を $N$ とする。

1. $N$ は部分群である。空積として $e\in N$ であり、二つの積を並べたものはまた同じ形である。逆元は $\prod_{i=k}^{1}u_ir_i^{-\epsilon_i}u_i^{-1}$（順序を逆にした積）で、同じ形である。
2. $N$ は正規である。$g\in G$ について $g\Bigl(\prod_iu_ir_i^{\epsilon_i}u_i^{-1}\Bigr)g^{-1}=\prod_i(gu_i)r_i^{\epsilon_i}(gu_i)^{-1}$ であり、これは $N$ に属する（隣り合う因子の間に $g^{-1}g$ を挿入した）。
3. $R\subseteq N$ は明らかである。
4. $R$ を含む任意の正規部分群 $M$ は各因子 $u_ir_i^{\epsilon_i}u_i^{-1}$ を含むので $N\subseteq M$ である。

2, 3 より $N$ は $R$ を含む正規部分群なので $\langle\!\langle R\rangle\!\rangle\subseteq N$ であり、4 より逆の包含も成り立つ。
:::

従って表示群で $w=1$ となる条件は、自由群で $w$ が上の有限積に書けることである。このとき $w=1$ を関係語の **帰結** という。例えば $a^6=a^{10}=1$ から $a^4=a^{10}(a^6)^{-1}=1$、さらに $a^2=a^6(a^4)^{-1}=1$ を得る。

### von Dyck の定理

表示群から他の群への準同型は、生成元の行き先を指定し、関係語が単位元に写ることを確かめれば得られる。これが表示を扱うときの基本的な道具である。

:::theorem[von Dyck の定理（表示群の普遍性）]{#thm-grp-presented-universal}
$P=\langle X\mid R\rangle$ とし、$\pi:F(X)\to P$ を射影とする。群 $G$ と写像 $f:X\to G$ に対し、$\varphi\circ\pi|_X=f$ を満たす準同型 $\varphi:P\to G$ が存在するための必要十分条件は、すべての $r\in R$ について代入準同型の値 $\hat f(r)$ が $e_G$ に等しいことである。このとき $\varphi$ は一意であり、その像は $f(X)$ が生成する部分群である。特に $G$ が $f(X)$ で生成され、$G$ で関係式がすべて成り立つならば、$\varphi$ は全射であり $G$ は $P$ の剰余群である。
:::

:::proof
**必要性。** $\varphi\circ\pi$ と $\hat f$ はどちらも $X$ 上で $f$ に一致する $F(X)$ からの準同型なので、:ref[thm-grp-free-universal] の一意性から等しい。$r\in R$ なら $\pi(r)=1$ だから $\hat f(r)=\varphi(\pi(r))=e_G$ である。

**十分性。** $\Ker\hat f$ は $F(X)$ の正規部分群であり、仮定から $R$ を含む。正規閉包の最小性から $\langle\!\langle R\rangle\!\rangle\subseteq\Ker\hat f$ である。:ref[prop-grp-quotient-universal] により、$\hat f=\varphi\circ\pi$ を満たす準同型 $\varphi:P\to G$ が存在する。

**一意性と像。** $P$ は $\pi(X)$ で生成されるので、:ref[prop-grp-hom-rules] より $\varphi$ は $\pi(X)$ 上の値で決まる。像は $\varphi(P)=\hat f(F(X))$ であり、$F(X)$ が $X$ で生成されることから、これは $f(X)$ で生成される部分群である。
:::

:::remark
von Dyck の定理が与えるのは準同型であって同型ではない。関係式を満たす具体的な群（**モデル**）を一つ見つけても、表示群がそのモデルと同型であるとは限らない。例えば $\langle a,b\mid aba^{-1}b^{-1}\rangle$ の関係式は $C_2\times C_2$ の二つの生成元でも成り立つが、この表示群は $\Z^2$ である（:ref[exm-grp-free-abelian-presentation]）。同型を示すには、表示群の大きさを「上から」抑える議論が別に必要である。
:::

### 表示の存在と同定の方法

:::proposition[表示の存在]{#prop-grp-presentation-exists}
1. 任意の群 $G$ と生成系 $S$ について、$F(S)\to G$ の核を $K$ とすると $G\cong\langle S\mid K\rangle$ である。
2. 有限群 $G$ は有限表示をもつ。実際、記号 $x_g$（$g\in G$）を生成元とし、$x_gx_hx_{gh}^{-1}$（$g,h\in G$）を関係語とする表示（**乗積表による表示**）は $G$ を与える。
:::

:::proof
1. :ref[cor-grp-quotient-of-free] の全射の核 $K$ は正規なので $\langle\!\langle K\rangle\!\rangle=K$ であり、:ref[thm-grp-first-iso] から $F(S)/K\cong G$ である。

2. $x_g\mapsto g$ は乗積表を保つので表示群 $P$ から $G$ への全射を与える。逆に $g\mapsto x_g$ は $x_gx_h=x_{gh}$ により積を保ち、$x_e^2=x_e$ から $x_e=1$ である。この準同型との二つの合成は各 $g$、各生成元 $x_g$ を固定するので恒等写像である。
:::

乗積表を使うと関係式は $|G|^2$ 本に及ぶ。より少数の関係で足りることは、次の方法で確かめられる。

:::proposition[標準形による表示の同定]{#prop-grp-identify-presentation}
$P=\langle X\mid R\rangle$、$\pi:F(X)\to P$ を射影とする。群 $H$ と写像 $f:X\to H$ が次を満たすとする。

1. （モデル）すべての $r\in R$ について $\hat f(r)=e$ であり、$f(X)$ は $H$ を生成する。
2. （上からの評価）語の集合 $W\subseteq F(X)$ があって、$P$ の任意の元は $\pi(w)$（$w\in W$）の形に書ける。
3. （区別）$\hat f$ は $W$ 上で単射である。

このとき von Dyck の定理の準同型 $\varphi:P\to H$ は同型であり、$\pi|_W:W\to P$ は全単射である。すなわち $W$ は $P$ の元の **標準形** の集合である。
:::

:::proof
1 と von Dyck の定理から全射準同型 $\varphi$ があり、$\varphi\circ\pi=\hat f$ である。$p,p'\in P$ が $\varphi(p)=\varphi(p')$ を満たすとし、2 により $p=\pi(w),\ p'=\pi(w')$（$w,w'\in W$）と書く。すると $\hat f(w)=\varphi(p)=\varphi(p')=\hat f(w')$ であり、3 より $w=w'$、従って $p=p'$ である。よって $\varphi$ は単射である。$\pi|_W$ は 2 より全射で、$\varphi\circ\pi|_W=\hat f|_W$ が単射だから単射でもある。
:::

$P$ が有限なら、2 は $|P|\le|W|$ を、1 と 3 は $|P|\ge|\varphi(P)|=|H|\ge|W|$ を意味する。「上からと下から挟む」議論である。

:::example[巡回群の表示]{#exm-grp-cyclic-presentation}
$n\ge1$ とし、$P=\langle a\mid a^n\rangle$ とする。$H=C_n=\Z/n\Z$、$f(a)=\bar1$ とすると、$\hat f(a^n)=\bar n=\bar0$ なのでモデルである。$P$ の元は $a^k$（$k\in\Z$）であり、$a^n=1$ より $a^k=a^{k\bmod n}$ だから $W=\{1,a,\ldots,a^{n-1}\}$ が条件 2 を満たす。$\hat f(a^k)=\bar k$ は $0\le k<n$ で相異なる。従って $\langle a\mid a^n\rangle\cong C_n$ である。関係語のない $\langle a\mid\ \rangle=F(a)$ は $\Z$ である。

上で求めた $a^2=1$ から逆に $a^6=a^{10}=1$ も従うので、$\langle a\mid a^6,a^{10}\rangle=\langle a\mid a^2\rangle\cong C_2$ である。一般に正整数 $m,n$ では $\langle a\mid a^m,a^n\rangle\cong C_{\gcd(m,n)}$ となる（演習問題）。
:::

:::example[自由 Abel 群の表示]{#exm-grp-free-abelian-presentation}
$P=\langle a,b\mid aba^{-1}b^{-1}\rangle$ とする。$P$ の中で $ab=ba$ であり、生成元の逆元どうしも可換になる（$ab=ba$ の両辺に左右から $a^{-1}$ を掛けると $ba^{-1}=a^{-1}b$、など）。従って $P$ は可換群で、任意の元は $a^ib^j$（$i,j\in\Z$）と書ける。モデルを $\Z^2$、$f(a)=(1,0),\ f(b)=(0,1)$ とすると関係語は $(0,0)$ に写り、$\hat f(a^ib^j)=(i,j)$ は相異なる。従って $P\cong\Z^2$ であり、$a^ib^j$ が標準形である。$n$ 個の生成元についても、すべての対の交換子を関係語とする表示は $\Z^n$ を与える（:ref[prop-grp-product-presentation] の後の注意）。

関係式を一つ変えるだけで群は大きく変わる。$\langle a,b\mid ab\rangle$ では $b=a^{-1}$ が強制され、$a\mapsto1$ で $\Z$ と同型になる（$W=\{a^k\}$、モデル $\Z$）。また $\langle a,b\mid a=b\rangle$ も $\Z$ である。
:::

:::example[自明群になる表示]{#exm-grp-trivial-presentation}
表示から群の大きさを見抜くことは一般には難しい。$P=\langle a,b\mid aba^{-1}=b^2,\ bab^{-1}=a^2\rangle$ を考える。第一式から $ab=b^2a$、第二式から $ba=a^2b$ である。従って
$$
ab=b^2a=b(ba)=b(a^2b)=ba^2b
$$
となり、右から $b^{-1}$ を掛けて $a=ba^2$、さらに右から $a^{-1}$ を掛けて $1=ba$、すなわち $b=a^{-1}$ を得る。これを第二式に代入すると $a^{-1}aa=a^2$、すなわち $a=a^2$ で $a=1$、従って $b=1$ である。よって $P$ は自明群である。二つの関係式は一見無害だが、組み合わせると生成元を消してしまう。
:::

### 関係式の追加と直積

:::proposition[関係式の追加]{#prop-grp-add-relations}
$P=\langle X\mid R\rangle$、$\pi:F(X)\to P$ を射影とし、$S\subseteq F(X)$ とする。このとき
$$
\langle X\mid R\cup S\rangle\cong P\big/\langle\!\langle\pi(S)\rangle\!\rangle_P
$$
である。すなわち、関係式を追加することは、その像の正規閉包で割ることに等しい。
:::

:::proof
$N=\langle\!\langle R\rangle\!\rangle$、$M=\langle\!\langle R\cup S\rangle\!\rangle$（$F(X)$ の中の正規閉包）とすると $N\subseteq M$ であり、:ref[thm-grp-third-iso] より $F(X)/M\cong(F(X)/N)/(M/N)=P/\pi(M)$ である。$\pi(M)=\langle\!\langle\pi(S)\rangle\!\rangle_P$ を示せばよい。$\pi$ は全射なので $\pi(M)$ は $P$ の正規部分群であり、$\pi(S)$ を含む。従って $\pi(M)\supseteq\langle\!\langle\pi(S)\rangle\!\rangle_P$ である。逆に $\pi^{-1}\bigl(\langle\!\langle\pi(S)\rangle\!\rangle_P\bigr)$ は $F(X)$ の正規部分群で $N=\Ker\pi$ と $S$ を含むので $M$ を含む。両辺を $\pi$ で写して $\pi(M)\subseteq\langle\!\langle\pi(S)\rangle\!\rangle_P$ を得る。
:::

例えば位数 $12$ の $D_6$ に $r^3=1$ を加えると、中心にある $\langle r^3\rangle$ で割るため位数は $6$ になる。残る関係は $r^3=s^2=srsr=1$ であり、得られる群は $D_3$ である（:ref[exm-grp-dihedral-presentation]）。これに対し $s=1$ を加えると $srsr=1$ から $r^2=1$ が従い、$\langle r\mid r^6,r^2\rangle\cong C_2$ だけが残る。

:::proposition[直積の表示]{#prop-grp-product-presentation}
$P=\langle X\mid R\rangle$、$Q=\langle Y\mid S\rangle$（$X\cap Y=\varnothing$）とすると
$$
P\times Q\cong\bigl\langle X\cup Y\ \big|\ R,\ S,\ [x,y]\ (x\in X,\ y\in Y)\bigr\rangle
$$
である。
:::

:::proof
右辺を $G$ とする。$x\mapsto(x,e),\ y\mapsto(e,y)$ は $P\times Q$ で関係式をすべて満たす（$R$ は第一成分、$S$ は第二成分で成り立ち、$(x,e)$ と $(e,y)$ は可換）ので、von Dyck の定理により準同型 $\Phi:G\to P\times Q$ がある。

逆向きを作る。$x\mapsto x$ は $G$ の中で $R$ を満たすので準同型 $\alpha:P\to G$ があり、同様に $\beta:Q\to G$ がある。$\alpha(P)$ の元と $\beta(Q)$ の元は可換である。実際、各 $y$ の像の中心化群は部分群であり $\alpha(X)$ を含むので $\alpha(P)$ を含む。従って $\alpha(P)$ の各元の中心化群は $\beta(Y)$ を含み、$\beta(Q)$ を含む。よって $\Psi(p,q)=\alpha(p)\beta(q)$ は
$$
\Psi\bigl((p,q)(p',q')\bigr)=\alpha(p)\alpha(p')\beta(q)\beta(q')=\alpha(p)\beta(q)\alpha(p')\beta(q')
$$
を満たし、準同型である。$\Psi\circ\Phi$ と $\Phi\circ\Psi$ は生成元を固定するので恒等写像であり、$\Phi$ は同型である。
:::

帰納法により、$\langle x_1,\ldots,x_n\mid[x_i,x_j]\ (i<j)\rangle\cong\Z^n$、$\langle a,b\mid a^m,b^n,[a,b]\rangle\cong C_m\times C_n$ が従う。

### 表示から可換化を求める

表示群そのものの構造は分かりにくくても、その可換化は整数行列の計算で求まる。語 $w$ の指数和ベクトルを $\varepsilon(w)\in\Z^{(X)}$ とする（:ref[exm-grp-free-abelianization]）。

:::proposition[表示群の可換化]{#prop-grp-presentation-abelianization}
$P=\langle X\mid R\rangle$ とし、$L\subseteq\Z^{(X)}$ を $\{\varepsilon(r):r\in R\}$ が生成する部分群とする。このとき
$$
P_{\mathrm{ab}}\cong\Z^{(X)}/L
$$
であり、同型は $x$ の像を $e_x+L$ に送る。
:::

:::proof
$A=\Z^{(X)}/L$ とおく。写像 $x\mapsto e_x+L$ の代入準同型は $w\mapsto\varepsilon(w)+L$ であり、関係語 $r$ を $\varepsilon(r)+L=0$ に送る。von Dyck の定理により準同型 $\theta:P\to A$ があり、$A$ は可換なので :ref[prop-grp-abelianization] により $\bar\theta:P_{\mathrm{ab}}\to A$ を引き起こす。

逆向きに、$\kappa:\Z^{(X)}\to P_{\mathrm{ab}}$ を $\kappa\bigl(\sum n_xe_x\bigr)=\prod\bar x^{\,n_x}$ で定める（$\bar x$ は $x$ の $P_{\mathrm{ab}}$ での像、可換なので順序によらない）。$\kappa\circ\varepsilon$ と、$F(X)\to P\to P_{\mathrm{ab}}$ の合成はともに $x\mapsto\bar x$ を満たす準同型なので一致する。従って $\kappa(\varepsilon(r))$ は $r$ の $P_{\mathrm{ab}}$ での像、すなわち単位元であり、$\kappa(L)=0$ である。よって $\kappa$ は $\bar\kappa:A\to P_{\mathrm{ab}}$ を引き起こす。$\bar\theta\bar\kappa$ と $\bar\kappa\bar\theta$ は生成元を固定するので恒等写像である。
:::

$X=\{x_1,\ldots,x_n\}$ と $R=\{r_1,\ldots,r_m\}$ が有限なら、$\varepsilon(r_i)$ を行とする $m\times n$ 整数行列（**関係行列**）を作り、行と列の整数基本変形で Smith 標準形（:ref[thm-grp-smith]）に対角化すれば、$\Z^n/L$ の不変因子が読み取れる（[直積・半直積と有限 Abel 群](/math/group-theory/products-abelian/) の「生成元と関係式から群を決定する」、:ref[thm-grp-fg-abelian]）。

:::example[可換化の計算]{#exm-grp-abelianization-computation}
1. $P=\langle a,b\mid a^2b^{-3}\rangle$（関係式 $a^2=b^3$）の関係行列は $(2\ \ {-3})$ である。$(i,j)\mapsto3i+2j$ は $\Z^2\to\Z$ の全射で、その核は $(2,-3)$ で生成される（$3i+2j=0$ なら $i=2k,\ j=-3k$）。従って $P_{\mathrm{ab}}\cong\Z$ であり、$P$ は無限群である。一方 $a\mapsto(12),\ b\mapsto(123)$ は $(12)^2=e=(123)^3$ を満たすので $P$ は $S_3$ に全射し、非可換である。
2. $P=\langle a,b\mid a^2,b^3,(ab)^3\rangle$ の関係行列の行は $(2,0),(0,3),(3,3)$ である。$(3,3)-(2,0)=(1,3)$、$(1,3)-(0,3)=(1,0)$ より $L$ は $(1,0),(0,3)$ で生成され、$P_{\mathrm{ab}}\cong C_3$ である。後で見るように $P\cong A_4$ であり（:ref[exm-grp-a4-presentation]）、$A_4/V_4\cong C_3$ と整合する。
3. $\langle r,s\mid r^n,s^2,srsr\rangle$ の行は $(n,0),(0,2),(2,2)$ で、$L$ は $(n,0),(0,2),(2,0)$ で生成される。従って可換化は $n$ が奇数なら $C_2$、偶数なら $C_2\times C_2$ である。
:::

## Tietze 変換

同じ群に多くの表示がある。計算しやすい生成元を付け加えたり、既に従う関係式を取り除いたりして表示を変形したい。どのような変形が群を変えないかを定め、それで十分であることを示す。

:::definition[Tietze 変換]{#def-grp-tietze}
表示 $\langle X\mid R\rangle$ に対する次の四種類の操作を :term[Tietze 変換]{en="Tietze transformation" yomi="てぃーつぇへんかん"} という。

- **(T1) 関係語の追加：** $r\in\langle\!\langle R\rangle\!\rangle$ を満たす語 $r\in F(X)$ を $R$ に加え、$\langle X\mid R\cup\{r\}\rangle$ とする。
- **(T2) 関係語の削除：** $r\in R$ が $r\in\langle\!\langle R\setminus\{r\}\rangle\!\rangle$ を満たすとき、$\langle X\mid R\setminus\{r\}\rangle$ とする。
- **(T3) 生成元の追加：** $X$ に属さない新しい記号 $y$ と語 $w\in F(X)$ をとり、$\langle X\cup\{y\}\mid R\cup\{yw^{-1}\}\rangle$ とする（$y$ を $w$ と **定義** する）。
- **(T4) 生成元の削除：** 表示が $\langle X\cup\{y\}\mid R\cup\{yw^{-1}\}\rangle$ の形で、$w\in F(X)$ かつ $R\subseteq F(X)$（すなわち $y$ は関係語 $yw^{-1}$ 以外に現れない）のとき、$\langle X\mid R\rangle$ とする。

(T2), (T4) はそれぞれ (T1), (T3) の逆操作である。
:::

:::proposition[Tietze 変換は群を変えない]{#prop-grp-tietze-sound}
各 Tietze 変換の前後の表示群は同型である。より詳しく、(T1), (T2) の前後では正規閉包が等しく、(T3) の前後では $x\mapsto x$（$x\in X$）が同型 $\langle X\mid R\rangle\to\langle X\cup\{y\}\mid R\cup\{yw^{-1}\}\rangle$ を与え、その逆は $x\mapsto x,\ y\mapsto w$ である。
:::

:::proof
**(T1), (T2)。** $R$ の正規閉包に既に入っている $r$ を加えても、それを含む正規部分群の族は変わらない。従って共通部分である正規閉包は同じである。(T2) はこの操作を逆にたどる。

**(T3), (T4)。** 旧表示群を $P$、新しいものを $P'$ とする。von Dyck の定理により $x\mapsto x$ は $\Phi:P\to P'$ を与え、$x\mapsto x,\ y\mapsto w$ は $\Psi:P'\to P$ を与える。後者は旧関係を保ち、新関係 $yw^{-1}$ を $1$ に写すからである。$\Psi\Phi$ は旧生成元を固定する。$\Phi\Psi$ はさらに $y$ を $w$ に送るが、$P'$ では $y=w$ である。従って両方の合成が恒等写像であり、(T4) には逆の同型を使える。
:::

実際の計算では、次の操作を一つの手順として扱うと便利である。

:::lemma[派生する変形]{#lem-grp-tietze-derived}
次の操作は (T1), (T2) の組合せで実現でき、従って群を変えない。

1. 関係語 $r$ を、その逆 $r^{-1}$、共役 $uru^{-1}$（$u\in F(X)$）、または別の関係語 $r'$ との積 $rr'$ に置き換える。
2. 関係語 $yw^{-1}$（$y$ は生成元、$w$ は $y$ を含まない語）があるとき、ほかの関係語 $r$ の中の $y$ をすべて $w$ に置き換えた語 $r'$ で $r$ を置き換える。
:::

:::proof
1. 新しい語は $\langle\!\langle R\rangle\!\rangle$ に属するので (T1) で加えられる。その後、元の $r$ は新しい語と残りの関係語から従う（$r=(r^{-1})^{-1}$、$r=u^{-1}(uru^{-1})u$、$r=(rr')r'^{-1}$）ので (T2) で削除できる。

2. $N=\langle\!\langle R\rangle\!\rangle$ とすると、$F(X)/N$ の中では $y=w$ だから $r$ と $r'$ の像は等しく、$r'r^{-1}\in N$、従って $r'\in N$ である。(T1) で $r'$ を加える。逆に、関係語 $r'$ と $yw^{-1}$ から $r$ が従う（同じ理由で $r$ と $r'$ は $yw^{-1}$ を含む正規閉包を法として等しい）ので、(T2) で $r$ を削除できる。
:::

:::theorem[Tietze の定理]{#thm-grp-tietze}
二つの有限表示 $\langle X\mid R\rangle$ と $\langle Y\mid S\rangle$ が同型な群を与えるための必要十分条件は、一方から他方へ有限回の Tietze 変換で移れることである（生成元の記号の名前の付け替えは許す）。
:::

:::proof
十分性は :ref[prop-grp-tietze-sound] による。必要性を示す。記号の名前を付け替えて $X\cap Y=\varnothing$ としてよい。$P=\langle X\mid R\rangle$、$Q=\langle Y\mid S\rangle$ とし、同型 $\theta:P\to Q$ を固定する。各 $y\in Y$ について $\theta(\bar u_y)=\bar y$ を満たす語 $u_y\in F(X)$ を、各 $x\in X$ について $\theta(\bar x)=\bar v_x$ を満たす語 $v_x\in F(Y)$ を選ぶ（バーは表示群での像）。

**ステップ 1。** $P$ の表示に (T3) を $|Y|$ 回施し、
$$
P_1=\bigl\langle X\cup Y\ \big|\ R,\ yu_y^{-1}\ (y\in Y)\bigr\rangle
$$
を得る。:ref[prop-grp-tietze-sound] の同型 $\Psi:P_1\to P$（$x\mapsto\bar x,\ y\mapsto\bar u_y$）と $\theta$ の合成 $\theta\circ\Psi:P_1\to Q$ は同型で、$y\mapsto\theta(\bar u_y)=\bar y$、$x\mapsto\bar v_x$ を満たす。

**ステップ 2。** $Q$ で各 $s\in S$ と $xv_x^{-1}$ は単位元に写る。$\theta\Psi$ が単射だから、これらは $P_1$ でも単位元であり、関係語の正規閉包に属する。従って (T1) で不足する語を加えて
$$
P_2=\bigl\langle X\cup Y\ \big|\ R,\ yu_y^{-1}\ (y\in Y),\ S,\ xv_x^{-1}\ (x\in X)\bigr\rangle
$$
を得る。

**ステップ 3。** 逆方向から $Q$ に (T3) で $x=v_x$ を加えた表示を $Q_1$ とする。$\theta^{-1}$ に同じ議論を適用すると、$R$ の元と $yu_y^{-1}$ は $Q_1$ の関係語の正規閉包に属する。そこで $P_2$ から、$Q_1$ に残す語をすべて保ちながら余分な関係語だけを (T2) で取り除けば $Q_1$ になる。

**ステップ 4。** $Q_1$ では、各 $x$ は関係語 $xv_x^{-1}$ にしか現れず、$v_x$ は $Y$ 上の語である。従って (T4) を $|X|$ 回施して $x$ をすべて削除でき、$\langle Y\mid S\rangle$ に達する。

$P$ から $Q$ へ進む経路では (T3) は $|Y|$ 回、(T4) は $|X|$ 回である。(T1) は高々 $|S|+|X|$ 回、(T2) は高々 $|R|+|Y|$ 回であり、重複する関係語は一度だけ扱えばよい。有限表示だからすべて有限回で済む。
:::

無限表示でも関係語の族や新生成元とその定義の族を一度に扱えば、同じ証明が通る。ただし証明は同型 $\theta$ を用いて語を選んでおり、二つの表示の同型性を判定する手続きを与えるわけではない。

:::example[定義した生成元の消去]{#exm-grp-tietze-eliminate}
$P=\langle a,b\mid b=a^2,\ b^3=1\rangle$ とする。:ref[lem-grp-tietze-derived] の 2 で関係語 $b^3$ の $b$ を $a^2$ に置き換えると $\langle a,b\mid ba^{-2},\ a^6\rangle$ となる。$b$ は関係語 $ba^{-2}$ にしか現れないので (T4) で削除でき、$P\cong\langle a\mid a^6\rangle\cong C_6$ を得る。

この例で $b^3$ への代入を省くと (T4) の条件を満たさない。$b$ を含む式をすべて捨てれば $\Z$ が残り、位数さえ変わる。
:::

:::example[生成元を積で取り替える]{#exm-grp-tietze-product}
$P=\langle a,b\mid a^2,\ b^3,\ aba^{-1}b^{-1}\rangle\cong C_2\times C_3$（:ref[prop-grp-product-presentation]）に、(T3) で $t=ab$ を加える。$P$ では $a,b$ が可換なので $t^3=a^3b^3=a$、$t^4=a^4b^4=b$、$t^6=1$ が成り立ち、(T1) で関係語 $at^{-3},\ bt^{-4},\ t^6$ を加えられる。逆に $t^6=1,\ a=t^3,\ b=t^4$ からは $a^2=t^6=1$、$b^3=t^{12}=1$、$ab=t^7=t=ba$ が従うので、元の三つの関係語と $t(ab)^{-1}$ を (T2) で削除できる。最後に (T4) で $a,b$ を削除すると $\langle t\mid t^6\rangle\cong C_6$ を得る。これは $C_2\times C_3\cong C_6$（:ref[exm-grp-crt]）の表示による証明である。
:::

:::example[三葉結び目の群]{#exm-grp-trefoil}
$P=\langle x,y\mid xyx=yxy\rangle$ と $Q=\langle a,b\mid a^2=b^3\rangle$ が同型であることを、Tietze 変換を明示して示す。関係語を $r=xyx(yxy)^{-1}$ と書く。

1. (T3) を二回施し、$a=xyx,\ b=xy$ を定義する：$\langle x,y,a,b\mid r,\ a(xyx)^{-1},\ b(xy)^{-1}\rangle$。
2. この群で $a=xyx=bx$ だから $x=b^{-1}a$、また $y=x^{-1}b=a^{-1}b\cdot b=a^{-1}b^2$ である。さらに $b^3=(xy)(xy)(xy)=(xyx)(yxy)=a\cdot a$ である（$yxy=xyx=a$ を用いた）。従って (T1) で関係語 $x(b^{-1}a)^{-1},\ y(a^{-1}b^2)^{-1},\ a^2b^{-3}$ を加えられる。
3. 逆に $x=b^{-1}a,\ y=a^{-1}b^2,\ a^2=b^3$ を仮定すれば $xy=b$、$xyx=a$、$yxy=y(xy)=a^{-1}b^3=a$ である。従って (T2) で $r,\ a(xyx)^{-1},\ b(xy)^{-1}$ を順に削除できる。
4. $x,y$ はそれぞれの定義関係にしか現れず、定義する語は $a,b$ だけでできているので、(T4) で削除して $\langle a,b\mid a^2b^{-3}\rangle=Q$ を得る。

:ref[exm-grp-abelianization-computation] の 1 より、この群は無限で非可換である。さらに $z=a^2=b^3$ は $a$ とも $b$ とも可換なので中心に属する。この群は位相幾何学では三葉結び目の補空間の基本群として現れ :cite[hatcher-at]、また $x,y$ の表示は 3 本の紐の組紐群の標準的な表示でもある。
:::

## 例

本節では、これまでの方法（:ref[prop-grp-identify-presentation] による上下からの評価と、Tietze 変換）を使って、よく現れる群の表示を同定する。

### 二面体群

:::example[二面体群の表示]{#exm-grp-dihedral-presentation}
$n\ge3$ とし、$P=\langle r,s\mid r^n,\ s^2,\ srsr\rangle$（関係式 $r^n=s^2=1,\ srs=r^{-1}$）とする。このとき $P\cong D_n$ であり、$r^is^j$（$0\le i<n,\ j=0,1$）が標準形である。

**上からの評価。** $P$ の中では $s^{-1}=s$、$r^{-1}=r^{n-1}$ であり、$srs=r^{-1}$ から $sr=r^{-1}s$ である。$W=\{r^is^j:0\le i<n,\ j=0,1\}$（指数 $i$ は法 $n$ で読む）の像の集合 $U$ が、右から $r$ と $s$ を掛ける操作で閉じることを確かめる。$r^is^j\cdot s=r^is^{j+1}$ で、$s^2=1$ により $U$ に戻る。$r^i\cdot r=r^{i+1}$、$r^is\cdot r=r^ir^{-1}s=r^{i-1}s$ である。$U$ は単位元を含み、生成元 $r,s$ の逆元も $r,s$ の冪なので、生成元とその逆元の任意の積、すなわち $P$ 全体が $U$ に含まれる。

**モデル。** :ref[def-grp-dihedral] の回転 $r$ と鏡映 $s$ は :ref[prop-grp-dihedral-form] により三つの関係式を満たし、$D_n$ を生成する。同じ命題により $r^is^j$（$0\le i<n,\ j=0,1$）は相異なる。

:ref[prop-grp-identify-presentation] により $P\cong D_n$ である。
:::

同じ群を、二つの鏡映を生成元とする表示に書き換える。Tietze 変換を順に施す。

1. (T3) で $t=sr$ を加える：$\langle r,s,t\mid r^n,s^2,srsr,\ t(sr)^{-1}\rangle$。
2. この群では $t^2=srsr=1$、$st=s\cdot sr=r$、$(st)^n=r^n=1$ である。(T1) で $t^2,\ r(st)^{-1},\ (st)^n$ を加える。
3. 逆に $s^2=t^2=(st)^n=1$ と $r=st$ からは、$r^n=(st)^n=1$、$srsr=s\,st\,s\,st=t\cdot t=1$、$sr=s\,st=t$ が従う。(T2) で $r^n,\ srsr,\ t(sr)^{-1}$ を削除する。
4. $r$ は定義関係 $r=st$ にしか現れないので (T4) で削除する。

結果は
$$
D_n\cong\langle s,t\mid s^2,\ t^2,\ (st)^n\rangle \label{eq-grp-dihedral-coxeter}
$$
である。幾何学的には、$s$ と $t=sr$ は軸のなす角が $\pi/n$ の二つの鏡映であり（$t(i)=s(i+1)=-i-1$ は頂点の番号を反転する）、二つの鏡映の積は角 $2\pi/n$ の回転である。$n=3$ のとき、$S_3$ の互換 $s=(12),\ t=(23)$ は $s^2=t^2=e$、$(st)^3=(123)^3=e$ を満たして $S_3$ を生成するので、位数を比べて $S_3\cong D_3$ も分かる。

### 無限二面体群

:::example[無限二面体群]{#exm-grp-infinite-dihedral}
:eqref[eq-grp-dihedral-coxeter] から $(st)^n$ を取り除いた $P=\langle s,t\mid s^2,t^2\rangle$ を考える。

**モデル。** 整数全体 $\Z$ の全単射 $\sigma(x)=-x,\ \tau(x)=1-x$ は $\sigma^2=\tau^2=\id$ を満たす。$\rho=\sigma\tau$ は $\rho(x)=\sigma(1-x)=x-1$ で、$1$ だけずらす平行移動である。$\sigma,\tau$ が生成する群を $D_\infty$ と書く。

**上からの評価。** $P$ で $r=st$ とおくと、$s^{-1}=s,\ t^{-1}=t$ より $srs=ts=(st)^{-1}=r^{-1}$、従って $sr=r^{-1}s$ であり、また $t=sr$ である。$U=\{r^is^j:i\in\Z,\ j=0,1\}$ は右から $s$ を掛けても $t$ を掛けても閉じる：$r^i\cdot t=r^isr=r^{i-1}s$、$r^is\cdot t=r^is\,sr=r^{i+1}$。従って $P=U$ である。

**区別。** $r^i$ の像は $x\mapsto x-i$、$r^is$ の像は $x\mapsto-x-i$ である。前者は向きを保つ平行移動、後者は向きを逆にする写像であり、$0$ の像から $i$ が決まるので、$(i,j)$ が異なれば像も異なる。

:ref[prop-grp-identify-presentation] により $P\cong D_\infty$ であり、$st$ は無限位数である。$D_\infty$ は $\Z$ の距離を保つ全単射（$x\mapsto\pm x+c$）全体に一致し、平行移動の部分群 $\langle\rho\rangle\cong\Z$ と $\langle\sigma\rangle\cong C_2$ の半直積 $\Z\rtimes C_2$（$C_2$ は $\Z$ に $-1$ 倍で作用する、:ref[def-grp-semidirect]）である。また :ref[prop-grp-add-relations] により、$D_n$ は $D_\infty$ を $(st)^n$ の正規閉包で割った群である。
:::

二つの群 $G=\langle X\mid R\rangle$, $H=\langle Y\mid S\rangle$（$X\cap Y=\varnothing$）に対して、交換子 $[x,y]$ を加えない $\langle X\cup Y\mid R\cup S\rangle$ を $G$ と $H$ の **自由積** $G*H$ という。直積（:ref[prop-grp-product-presentation]）と違い、二つの因子の間に関係を何も課さない。von Dyck の定理により、$G*H$ から群 $K$ への準同型は、$G\to K$ と $H\to K$ の準同型の組と一対一に対応する。この性質は表示によらない言葉で述べられているので、:ref[prop-grp-free-unique] と同じ議論により、$G*H$ は $G,H$ の表示の選び方によらず同型を除いて定まる。$D_\infty=C_2*C_2$ はその最も簡単な例である。

### 四元数群

:::example[四元数群の表示]{#exm-grp-quaternion-presentation}
$P=\langle a,b\mid a^4=1,\ b^2=a^2,\ bab^{-1}=a^{-1}\rangle$ とする。

**上からの評価。** $P$ では $ba=a^{-1}b$、$b^2=a^2$、$a^{-1}=a^3$、$b^{-1}=b^3=a^2b$ である（$b^4=a^4=1$）。$U=\{a^ib^j:0\le i<4,\ j=0,1\}$ について、$a^ib\cdot a=a^ia^{-1}b=a^{i-1}b$、$a^ib\cdot b=a^ia^2=a^{i+2}$ であり、$a^i\cdot a,\ a^i\cdot b$ は明らかに $U$ に入る。生成元の逆元は生成元の冪なので $P=U$、$|P|\le8$ である。

**モデル。** :ref[exm-grp-quaternion] の行列 $u,v$ は $u^4=I$、$v^2=-I=u^2$、$vuv^{-1}=u^{-1}$ を満たす（$vu=-uv$ と $v^{-1}=-v$ から $vuv^{-1}=-vuv=uv^2=-u=u^{-1}$）。$a\mapsto u,\ b\mapsto v$ の像 $u^iv^j$ は $I,u,-I,-u,v,uv,-v,-uv$ で相異なり、$Q_8$ 全体である。

従って $P\cong Q_8$ である。関係式 $b^2=a^2$ を $b^2=1$ に替えると、:ref[exm-grp-dihedral-presentation] により $D_4$ になる。二つの位数 $8$ の非可換群は、$a^4=1$ と $bab^{-1}=a^{-1}$ を共有し、$b^2$ の値だけで区別される。
:::

### 交代群 A₄

:::example[四次交代群の表示]{#exm-grp-a4-presentation}
$P=\langle a,b\mid a^2,\ b^3,\ (ab)^3\rangle$ は四次交代群 $A_4$（:ref[def-grp-alternating]）と同型である。

**ステップ 1：二つの共役の積。** $c=bab^{-1},\ d=b^{-1}ab$ とおく。$(ab)^3=1$ の左から $a^{-1}=a$ を掛けると $babab=a$、従って
$$
bab=a(ab)^{-1}=ab^{-1}a \label{eq-grp-a4-key}
$$
である。$b^{-1}=b^2$ を用いると $c=bab\cdot b=ab^{-1}ab$ なので $ac=b^{-1}ab=d$ である。また :eqref[eq-grp-a4-key] を使うと $ca=bab^{-1}a=b\,(ab^{-1}a)=b\,(bab)=b^{-1}ab=d$ である。従って $ac=ca=d$ である。

**ステップ 2：正規な Klein 部分群。** $a^2=1$、$c^2=ba^2b^{-1}=1$ で $a,c$ は可換だから、$V=\{1,a,c,d\}$（$d=ac$）は積で閉じ（例えば $ad=c,\ cd=a,\ d^2=1$）、各元は自分自身が逆元なので部分群である。さらに
$$
bab^{-1}=c,\qquad bcb^{-1}=b^2ab^{-2}=b^{-1}ab=d,\qquad bdb^{-1}=a
$$
（$b^3=1$ より $b^2=b^{-1}$、$b^{-2}=b$）であり、$V$ は可換で $a\in V$ だから $aVa^{-1}=V$ である。正規化群 $N_P(V)$ は $a,b$ を含む部分群なので $P$ 全体であり、$V\trianglelefteq P$ である。

**ステップ 3：上からの評価。** $P/V$ は $a,b$ の像で生成され、$a$ の像は単位元、$b$ の像は位数 $3$ 以下なので $|P/V|\le3$ である。$|V|\le4$ と合わせて $|P|\le12$ である（:ref[thm-grp-lagrange]）。

**ステップ 4：モデル。** $A_4$ で $a\mapsto(12)(34)$、$b\mapsto(123)$ とすると、$((12)(34))^2=e$、$(123)^3=e$ であり、$(12)(34)(123)=(243)$（右の置換を先に施す）は位数 $3$ なので関係式が成り立つ。像の部分群 $K$ は $(12)(34)$ とその $(123)$ による共役 $(14)(23)$, $(13)(24)$ を含むので $V_4$ を含み、$(123)$ も含む。従って $|K|$ は $4$ と $3$ で割り切れ、$|K|=12$、$K=A_4$ である。

von Dyck の定理により全射 $P\to A_4$ があるので $|P|\ge12$、従って $P\cong A_4$ であり、$V$ は $A_4$ の $V_4$ に対応する。
:::

### 対称群

隣接互換 $\sigma_i=(i\ \ i+1)$（$1\le i<n$）は $S_n$ を生成し（:ref[lem-linalg-adjacent-inversions]）、次の関係を満たす。

- $\sigma_i^2=e$。
- $|i-j|\ge2$ なら $\sigma_i$ と $\sigma_j$ は互いに素な互換なので可換である。
- $\sigma_i\sigma_{i+1}\sigma_i$ は互換 $\sigma_{i+1}$ を $\sigma_i$ で共役したものなので $(i\ \ i+2)$ に等しく、$\sigma_{i+1}\sigma_i\sigma_{i+1}$ も $\sigma_i$ を $\sigma_{i+1}$ で共役したものなので $(i\ \ i+2)$ に等しい。

これらの関係だけで $S_n$ の積がすべて決まることを示す。

:::theorem[対称群の Coxeter 表示]{#thm-grp-symmetric-presentation}
$n\ge2$ とし、生成元 $s_1,\ldots,s_{n-1}$ と関係式
$$
\begin{aligned}
s_i^2 &=1 && (1\le i\le n-1),\\
s_is_j &=s_js_i && (|i-j|\ge2),\\
s_is_{i+1}s_i &=s_{i+1}s_is_{i+1} && (1\le i\le n-2)
\end{aligned}
$$
による表示群を $W_n$ とする。$s_i\mapsto(i\ \ i+1)$ は同型 $W_n\cong S_n$ を与える。
:::

:::proof
**ステップ 0：全射。** 上で確かめたように隣接互換は関係式を満たし $S_n$ を生成するので、von Dyck の定理により全射 $\varphi:W_n\to S_n$ がある。あとは $|W_n|\le n!$ を $n$ に関する帰納法で示せばよい。$n=2$ では $W_2=\langle s_1\mid s_1^2\rangle$ なので $|W_2|\le2$ である。

**ステップ 1：部分群 $H$。** $n\ge3$ とし、$H=\langle s_1,\ldots,s_{n-2}\rangle\le W_n$ とおく。$s_i\mapsto s_i$（$i\le n-2$）は $W_{n-1}$ の関係式を $W_n$ の中で満たすので、von Dyck の定理により $W_{n-1}$ から $H$ への全射がある。帰納法の仮定から $|H|\le|W_{n-1}|\le(n-1)!$ である。

**ステップ 2：剰余類の代表。** $u_n=1$、$1\le j<n$ では $u_j=s_{n-1}\cdots s_j$ とおく。$j=n$ なら $u_ns_i=s_i$ で、$i\le n-2$ のとき $H=Hu_n$ に入り、$i=n-1$ のとき $u_{n-1}$ に等しい。残る $1\le j<n$ について次の四つに分ける。

- $i\le j-2$ のとき：$s_i$ は $s_{n-1},\ldots,s_j$ のすべてと添字が $2$ 以上離れているので可換であり、$u_js_i=s_iu_j\in Hu_j$ である（$i\le n-3$ より $s_i\in H$）。
- $i=j-1$ のとき：$u_js_{j-1}=u_{j-1}$ である。
- $i=j$ のとき：$s_j^2=1$ より $u_js_j=s_{n-1}\cdots s_{j+1}=u_{j+1}$ である（$j=n-1$ なら $u_n=1$）。
- $i\ge j+1$ のとき：$u_j=s_{n-1}\cdots s_{i+1}\,s_is_{i-1}\,s_{i-2}\cdots s_j$ と書ける（$i=j+1$ なら $s_{i-2}\cdots s_j$ は空）。末尾に掛けた $s_i$ は $s_{i-2},\ldots,s_j$ と可換なので左へ移動でき、$s_is_{i-1}s_i=s_{i-1}s_is_{i-1}$ を使うと
$$
u_js_i=s_{n-1}\cdots s_{i+1}\,s_is_{i-1}s_i\,s_{i-2}\cdots s_j=s_{n-1}\cdots s_{i+1}\,s_{i-1}s_is_{i-1}\,s_{i-2}\cdots s_j
$$
となる。左に現れた $s_{i-1}$ は $s_{n-1},\ldots,s_{i+1}$ と可換なので先頭まで移動でき、$u_js_i=s_{i-1}u_j\in Hu_j$ である（$i-1\le n-2$ より $s_{i-1}\in H$）。


**ステップ 3：個数の評価。** $U=\bigcup_{j=1}^{n}Hu_j$ は $1=u_n$ を含み、ステップ 2 により各 $s_i$ を右から掛ける操作で閉じる（$hu_js_i=h(u_js_i)\in Hu_{j'}$）。$s_i^{-1}=s_i$ なので、$U$ は生成元とその逆元の任意の積を含み、$U=W_n$ である。従って $|W_n|\le n|H|\le n\cdot(n-1)!=n!$ である。

$\varphi$ は全射で $|W_n|\le|S_n|$ だから、$\varphi$ は全単射である。
:::

例えば $n=4$ で $u_1=s_3s_2s_1$ に右から $s_3$ を掛けると、$s_3s_2s_1s_3=s_3s_2s_3s_1=s_2s_3s_2s_1=s_2u_1$ となる。証明の第四の場合の操作である。$s_i^2=1$ のもとで第三の関係式は $(s_is_{i+1})^3=1$ と同値であり、第二の関係式は $(s_is_j)^2=1$ と同値である。従って $W_n=\langle s_1,\ldots,s_{n-1}\mid(s_is_j)^{m_{ij}}\rangle$（$m_{ii}=1$、$|i-j|=1$ なら $m_{ij}=3$、$|i-j|\ge2$ なら $m_{ij}=2$）と書ける。この形の表示で与えられる群を一般に **Coxeter 群** という。$n=3$ の場合は :eqref[eq-grp-dihedral-coxeter] の $n=3$ の場合と一致し、ここでも $S_3\cong D_3$ が確かめられる。

### Baumslag–Solitar 群

有限表示をもつ群でも、その部分群や元の振る舞いは自由群や有限群とはかなり異なりうる。

:::example[Baumslag–Solitar 群 BS(1,2)]{#exm-grp-baumslag-solitar}
$P=\langle a,b\mid bab^{-1}=a^2\rangle$ とする。$\Z[1/2]=\{k/2^m:k\in\Z,\ m\ge0\}$ とし、$\Z$ が $\Z[1/2]$ に $j\cdot q=2^jq$ で作用する半直積を $\Z[1/2]\rtimes\Z$ とすると、$P\cong\Z[1/2]\rtimes\Z$ である。

**上からの評価。** $P$ の中では $ba=a^2b$、従って帰納的に $ba^k=a^{2k}b$、$b^na^kb^{-n}=a^{2^nk}$（$n\ge0,\ k\in\Z$）である。また $bab^{-1}=a^2$ の左から $b^{-1}$ を掛けると $ab^{-1}=b^{-1}a^2$、従って $a^kb^{-1}=b^{-1}a^{2k}$ である。集合 $U=\{b^{-m}a^kb^n:m,n\ge0,\ k\in\Z\}$ が右からの $a^{\pm1},b^{\pm1}$ の積で閉じることを確かめる。

- $b^{-m}a^kb^n\cdot b=b^{-m}a^kb^{n+1}$。
- $n\ge1$ なら $b^{-m}a^kb^n\cdot b^{-1}=b^{-m}a^kb^{n-1}$、$n=0$ なら $b^{-m}a^k\cdot b^{-1}=b^{-(m+1)}a^{2k}$。
- $b^{-m}a^kb^n\cdot a^{\pm1}=b^{-m}a^k\,(b^na^{\pm1}b^{-n})\,b^n=b^{-m}a^{k\pm2^n}b^n$。

$U$ は単位元を含むので $P=U$ である。

**モデルと単射性。** 実数直線のアフィン変換 $\alpha(x)=x+1$、$\beta(x)=2x$ は $\beta\alpha\beta^{-1}(x)=2(x/2+1)=x+2=\alpha^2(x)$ を満たすので、準同型 $\varphi:P\to\mathrm{Aff}$、$a\mapsto\alpha,\ b\mapsto\beta$ がある。
$$
\varphi(b^{-m}a^kb^n)(x)=\frac{2^nx+k}{2^m}=2^{n-m}x+\frac{k}{2^m}
$$
である。$g=b^{-m}a^kb^n\in\Ker\varphi$ なら $n=m$ かつ $k=0$ で、$g=b^{-m}b^m=1$ である。従って $\varphi$ は単射である。

**像の記述。** 像は $x\mapsto2^jx+q$（$j\in\Z,\ q\in\Z[1/2]$）の形の変換全体である。実際、上の式から像はこの形であり、逆に $q=k/2^m$ と書いて $m\ge\max(0,-j)$ を選び $n=j+m\ge0$ とおけば実現される。この変換を組 $(q,j)$ で表すと、合成は
$$
(q,j)(q',j')=(q+2^jq',\ j+j')
$$
となり、これは :ref[def-grp-semidirect] の $\Z[1/2]\rtimes\Z$ の積である。

$bab^{-1}=a^2$ より $a$ と $a^2$ は共役である。平行移動の全体は傾きの指数 $(q,j)\mapsto j$ の核だから正規である。逆に平行移動 $x\mapsto x+k/2^m$ は $b^{-m}a^kb^m$ の像なので、平行移動全体が $a$ の正規閉包に一致する。その加法群 $\Z[1/2]$ では、有限個の元が生成する群はある $2^{-N}\Z$ に含まれ、$2^{-N-1}$ を含まない。従って二元で生成される $P$ に、有限生成でない正規部分群がある。
:::

### 関係式から読める情報と語の問題

与えられた語が表示群の単位元を表すかを判定する問題を :term[語の問題]{en="word problem" yomi="ごのもんだい"} という。

- 自由群では、語を簡約して空語になるかを見ればよい（:ref[lem-grp-reduction]）。
- 本節の例のように標準形とモデルが得られた群では、語を標準形に書き直せば判定できる。例えば $D_n$ では $sr^k=r^{-k}s$ を使って $r^is^j$ に直す。
- 一般の有限表示 $\langle X\mid R\rangle$ でも、:ref[prop-grp-normal-closure] により、$w=1$ ならば $w$ は関係語の共役とその逆の有限積に自由群の中で等しい。そこで「因子の個数と共役に用いる語の長さの合計」が $0,1,2,\ldots$ の積を順にすべて作り、簡約して $w$ の簡約語と比べる。各段階の候補は有限個なので、$w=1$ なら必ず有限時間で一致が見つかる。しかし $w\ne1$ のときこの探索は終わらない。
- 逆に $w\ne1$ を示すには、関係式を満たす群（例えば置換群）への準同型で $w$ の像が単位元でないものを見つければよい。:ref[exm-grp-abelianization-computation] で $\langle a,b\mid a^2=b^3\rangle$ が非可換であることを $S_3$ への写像で示したのはこの方法である。

実際、語の問題を解くアルゴリズムが存在しない有限表示群がある（Novikov–Boone の定理）。証明には計算の停止問題を群の関係式で表す構成が必要なので、ここでは扱わない。:cite[dummit-foote]

## 演習問題

::::exercise{level=basic}
語 $a^{-1}bb^{-1}aac^{-1}c$ に :ref[lem-grp-reduction] の証明の手続き $\rho$ を適用し、簡約語を求めよ。
:::solution
左から一文字ずつ読むと
$$
1\xrightarrow{\ast a^{-1}}a^{-1}\xrightarrow{\ast b}a^{-1}b\xrightarrow{\ast b^{-1}}a^{-1}\xrightarrow{\ast a}1\xrightarrow{\ast a}a\xrightarrow{\ast c^{-1}}ac^{-1}\xrightarrow{\ast c}a
$$
となり、簡約語は $a$ である。:ref[lem-grp-reduction] により、$bb^{-1}$ と $c^{-1}c$ を先に消してから $a^{-1}a$ を消すなど、別の順序で簡約しても同じ $a$ が得られる。
:::
::::

::::exercise{level=basic}
$F(a,b)$ からの代入準同型 $f:F(a,b)\to C_6$（$a\mapsto\bar2,\ b\mapsto\bar3$）と $g:F(a,b)\to S_3$（$a\mapsto(12),\ b\mapsto(123)$）について、$w=a^2b^{-1}a$ の像を求めよ。
:::solution
$C_6$ は加法で書くので $f(w)=2\cdot\bar2-\bar3+\bar2=\bar3$ である。$S_3$ では $g(w)=(12)^2(123)^{-1}(12)=(132)(12)$ である。右の $(12)$ を先に施すと $1\mapsto2\mapsto1$、$2\mapsto1\mapsto3$、$3\mapsto3\mapsto2$ なので $g(w)=(23)$ である。いずれも :ref[thm-grp-free-universal] により一意に定まる準同型の値である。
:::
::::

::::exercise{level=basic}
有限な基底 $X$ の大きさが $r\ge1$ のとき、長さ $k\ge1$ の簡約語の個数が $2r(2r-1)^{k-1}$ であることを示せ。また $F(a,b)$ で $\ell(w)\le2$ の元はいくつあるか。
:::solution
最初の文字は $X^{\pm}$ の $2r$ 個から自由に選べる。二文字目以降は、直前の文字の逆だけが禁止されるので $2r-1$ 通りずつある。簡約語であることは隣接する二文字の条件だけで決まるから、個数は $2r(2r-1)^{k-1}$ である。

$F(a,b)$（$r=2$）では長さ $0,1,2$ の簡約語がそれぞれ $1,4,12$ 個あり、異なる簡約語は異なる元を表すので、$1+4+12=17$ 個である。一般に長さ $k$ 以下の元は $1+\sum_{j=1}^k4\cdot3^{j-1}=2\cdot3^k-1$ 個であり、これは Cayley グラフ（:ref[prop-grp-cayley-tree]）で単位元から距離 $k$ 以内の頂点の数である。
:::
::::

::::exercise{level=basic}
$m,n\ge1$ について $\langle a\mid a^m,a^n\rangle\cong C_{\gcd(m,n)}$ を示せ。
:::solution
$d=\gcd(m,n)$ とし、Bézout の等式 $d=um+vn$（$u,v\in\Z$）をとる。$F(a)$ の中で $a^d=(a^m)^u(a^n)^v$ なので $a^d\in\langle\!\langle a^m,a^n\rangle\!\rangle$ である。逆に $a^m=(a^d)^{m/d}$、$a^n=(a^d)^{n/d}$ は $\langle\!\langle a^d\rangle\!\rangle$ に属する。従って二つの正規閉包は一致し、$\langle a\mid a^m,a^n\rangle=\langle a\mid a^d\rangle\cong C_d$ である（:ref[exm-grp-cyclic-presentation]）。例えば $\langle a\mid a^6,a^{10}\rangle\cong C_2$ である。
:::
::::

::::exercise{level=basic}
$\langle x,y\mid x^2,\ y^2,\ (xy)^2\rangle\cong C_2\times C_2$ を示せ。
:::solution
$F(x,y)$ の中で次の等式が成り立つ（右辺を展開して確かめられる）。
$$
(xy)^2=[x,y]\cdot\bigl(yx^2y^{-1}\bigr)\cdot y^2,\qquad
[x,y]=(xy)^2\cdot y^{-2}\cdot\bigl(yx^{-2}y^{-1}\bigr).
$$
例えば第一式の右辺は $xyx^{-1}y^{-1}\cdot yx^2y^{-1}\cdot y^2=xyx^{-1}x^2y=xyxy$ である。従って $\{x^2,y^2,(xy)^2\}$ と $\{x^2,y^2,[x,y]\}$ の正規閉包は互いに含み合い、一致する（:ref[prop-grp-normal-closure]）。後者の表示は :ref[prop-grp-product-presentation] により $\langle x\mid x^2\rangle\times\langle y\mid y^2\rangle\cong C_2\times C_2$ を与える。
:::
::::

::::exercise{level=basic}
$P=\langle a,b\mid a^3b^2,\ a^2b^{-2}\rangle$ の可換化を求めよ。
:::solution
:ref[prop-grp-presentation-abelianization] により $P_{\mathrm{ab}}\cong\Z^2/L$、$L$ は $(3,2)$ と $(2,-2)$ で生成される。$(3,2)-(2,-2)=(1,4)$、$(2,-2)-2(1,4)=(0,-10)$ なので $L$ は $(1,4),(0,10)$ で生成される。準同型 $\Z^2\to\Z/10\Z$、$(i,j)\mapsto j-4i$ は全射であり、核は $j\equiv4i\pmod{10}$ を満たす $(i,j)=i(1,4)+(0,j-4i)$ 全体、すなわち $L$ である。従って $P_{\mathrm{ab}}\cong C_{10}$ である。関係行列 $\begin{pmatrix}3&2\\2&-2\end{pmatrix}$ の成分の最大公約数が $1$、行列式が $-10$ であることからも不変因子 $1,10$ が読み取れる。
:::
::::

::::exercise{level=advanced}
$F(a,b)$ の自己準同型 $\varphi$ を $a\mapsto a^2,\ b\mapsto b$ で定める。$\varphi$ は単射であるが全射でないことを示せ。
:::solution
**単射性。** $w\ne1$ を標準形 $x_1^{n_1}\cdots x_m^{n_m}$（:ref[cor-grp-free-normal-form]）で書く。$\varphi(w)$ は各 $a^{n_i}$ を $a^{2n_i}$ に替えたものであり、$x_i\ne x_{i+1}$ と $2n_i\ne0$ が保たれるので、再び $m\ge1$ 個の因子をもつ標準形である。従って $\varphi(w)\ne1$ であり、$\varphi$ は単射である。

**全射でないこと。** 指数和の準同型 $\varepsilon:F(a,b)\to\Z^2$（:ref[exm-grp-free-abelianization]）を合成すると、$\varepsilon(\varphi(w))$ の第一成分は常に偶数である。$\varepsilon(a)=(1,0)$ なので $a$ は $\varphi$ の像に入らない。

従って $F_2$ は真部分群 $\langle a^2,b\rangle$ と同型である。単射な自己準同型の全射性が保証される有限群や有限次元ベクトル空間とは事情が異なる。
:::
::::

::::exercise{level=advanced}
$F=F(a,b)$ における $a$ の正規閉包 $N=\langle\!\langle a\rangle\!\rangle$ について、(1) $N\ne\langle a\rangle$ を示し、(2) $F/N\cong\Z$ を示し、(3) $N$ の自由基底を一つ求めよ。
:::solution
(1) $bab^{-1}\in N$ であるが、簡約語 $bab^{-1}$ は $b$ を含むので $a$ の冪 $a^k$（簡約語 $a^k$）とは異なる。

(2) $F/N=\langle a,b\mid a\rangle$ である。関係語 $a$ は「$a$ を空語で定義する」関係なので (T4) で $a$ を削除でき、$\langle b\mid\ \rangle\cong\Z$ を得る。同型は $a\mapsto0,\ b\mapsto1$ で与えられる。

(3) $N$ は全射 $F\to\Z$（$a\mapsto0,\ b\mapsto1$）の核である。剰余類 $Nb^n$ を $n\in\Z$ と書くと、Schreier グラフの正の辺は $b$ による $n\to n+1$ と、$a$ による各頂点のループ $n\to n$ である。親の辺を $b$ の辺（$n>0$ なら $n-1\to n$、$n<0$ なら $n\to n+1$）に選ぶと $t_n=b^n$ である。木にない正の辺は $a$ のループだけであり、:ref[thm-grp-free-subgroup] より
$$
\{\,b^nab^{-n}:n\in\Z\,\}
$$
が $N$ の自由基底である。特に $N$ は可算無限階数の自由群である。
:::
::::

::::exercise{level=advanced}
$P=\langle a,b\mid a^3,\ b^4,\ bab^{-1}=a^{-1}\rangle$ の位数が $12$ であり、$P$ が非可換であることを示せ。
:::solution
**上から。** $P$ では $ba=a^{-1}b$、従って $ba^i=a^{-i}b$ である。$U=\{a^ib^j:0\le i<3,\ 0\le j<4\}$ は右から $a$ を掛けると $a^ib^ja=a^{i+(-1)^j}b^j$、$b$ を掛けると $a^ib^{j+1}$ となって閉じ、生成元の逆元はその冪なので $P=U$、$|P|\le12$ である。

**下から。** $C_4=\langle\bar1\rangle$ が $C_3$ に $\bar1\cdot x=-x$ で作用する半直積 $C_3\rtimes C_4$（:ref[def-grp-semidirect]）を考える。$a\mapsto(\bar1,\bar0),\ b\mapsto(\bar0,\bar1)$ とすると $a^3,b^4$ は単位元に写り、$bab^{-1}\mapsto(-\bar1,\bar0)$ は $a^{-1}$ の像に等しい。$a^ib^j$ の像は $(\bar i,\bar j)$ で $12$ 個すべて異なる。:ref[prop-grp-identify-presentation] により $P\cong C_3\rtimes C_4$、$|P|=12$ である。

$ba=a^{-1}b$ で $a\ne a^{-1}$（$a$ の位数は $3$）なので $ab\ne ba$ であり、$P$ は非可換である。
:::
::::

::::exercise{level=advanced}
$P=\langle s,t\mid s^2,\ t^2,\ sts=tst\rangle$ を同定せよ。
:::solution
関係語を $r=sts(tst)^{-1}$ とし、$N=\langle\!\langle s^2,t^2\rangle\!\rangle\le F(s,t)$、$\bar G=F(s,t)/N=\langle s,t\mid s^2,t^2\rangle$ とおく。$\bar G$ の中では $s^{-1}=s,\ t^{-1}=t$ なので
$$
r=sts\,t^{-1}s^{-1}t^{-1}=sts\,tst=(st)^3\qquad(\bar G\ \text{の中で})
$$
である。すなわち $r$ と $(st)^3$ の $\bar G$ での像は等しく、$r\bigl((st)^3\bigr)^{-1}\in N$ である。従って $(st)^3\in\langle\!\langle s^2,t^2,r\rangle\!\rangle$ かつ $r\in\langle\!\langle s^2,t^2,(st)^3\rangle\!\rangle$ であり、二つの正規閉包は一致する。(T1) で $(st)^3$ を加え、(T2) で $r$ を削除すれば
$$
P=\langle s,t\mid s^2,\ t^2,\ (st)^3\rangle\cong D_3
$$
を得る（:eqref[eq-grp-dihedral-coxeter]）。$D_3\cong S_3$ であり、具体的には $s\mapsto(12),\ t\mapsto(23)$ が同型を与える。これは :ref[thm-grp-symmetric-presentation] の $n=3$ の場合でもある。
:::
::::

::::exercise{level=advanced}
$H$ を代入準同型 $F(a,b)\to C_3$（$a\mapsto\bar1,\ b\mapsto\bar0$）の核とする。Schreier グラフを使って $H$ の自由基底を一つ求め、指数公式と比較せよ。
:::solution
剰余類 $H,Ha,Ha^2$ を $0,1,2$ と書く。正の辺は $a$ による三角形 $i\to i+1\pmod3$ と、$b$ による各頂点のループである。全域木として二辺 $(0,a),(1,a)$ を選ぶ。この木で基点から読む語は $t_0=1,t_1=a,t_2=a^2$ である。木に入らないのは $(2,a):2\to0$ と三本のループであり、対応する元は
$$
s_{(2,a)}=a^2\cdot a\cdot1=a^3,\quad s_{(0,b)}=b,\quad s_{(1,b)}=aba^{-1},\quad s_{(2,b)}=a^2ba^{-2}
$$
である。:ref[thm-grp-free-subgroup] によりこれらが $H$ の自由基底であり、$\rank H=4=1+3(2-1)$ は指数公式に一致する。（親の辺として $(2,a)$ を選び $t_2=a^{-1}$ とすれば、別の自由基底 $a^3,\ b,\ aba^{-1},\ a^{-1}ba$ が得られる。）
:::
::::

::::exercise{level=advanced}
自由群 $F$ の元 $g,h$ と整数 $k\ne0$ について、$g^k=h^k$ なら $g=h$ であることを示せ。$S_3$ ではこれが成り立たない例を挙げよ。
:::solution
$g=1$ なら $h^k=1$ で、:ref[prop-grp-free-torsion-free] より $h=1=g$ である。$g\ne1$ とし、$c=g^k$ とおくと同じ命題より $c\ne1$ である。$g$ は $g^k=c$ と可換で、$h$ は $h^k=c$ と可換なので、$g,h\in C_F(c)$ である。:ref[cor-grp-free-commuting] の 4 より $C_F(c)=\langle z\rangle$ は無限巡回群であり、$g=z^m,\ h=z^n$ と書ける。$z^{mk}=z^{nk}$ と $z$ が無限位数であることから $mk=nk$、$k\ne0$ より $m=n$、すなわち $g=h$ である。

$S_3$ では $(12)^2=e=(13)^2$ であるが $(12)\ne(13)$ である。
:::
::::

::::exercise{level=advanced}
$P=\langle a,b\mid a^2,\ b^3,\ (ab)^4\rangle$ から $S_4$ への全射が存在することを示し、$P$ の可換化を求めよ。
:::solution
$a\mapsto(12),\ b\mapsto(234)$ とすれば $(12)(234)=(1234)$ なので関係式を満たす。像は $(12)$ と $(1234)$ を含むので、$(12)$ をこの四巡回の一乗、二乗で共役した $(23),(34)$ も含む。三つの隣接互換は $S_4$ を生成するから、von Dyck の定理が与える準同型は全射である。

可換化は :ref[prop-grp-presentation-abelianization] により $\Z^2/L$、$L$ は $(2,0),(0,3),(4,4)$ で生成される。$(4,4)-2(2,0)=(0,4)$、$(0,4)-(0,3)=(0,1)$ より $L$ は $(2,0),(0,1)$ で生成され、$P_{\mathrm{ab}}\cong C_2$ である。これは $S_4$ の可換化 $S_4/A_4\cong C_2$ と整合する。
:::
::::

## 確認テスト

::::quiz
$F(X)$ の関係語の集合 $R$ をすべて単位元にする剰余群を作るとき、割るべき最小の正規部分群はどれか。

- [ ] $R$ が生成する部分群 $\langle R\rangle$
- [x] $R$ の正規閉包 $\langle\!\langle R\rangle\!\rangle$
- [ ] $F(X)$ の交換子群
- [ ] $F(X)$ の中心

:::explanation
剰余群 $F(X)/N$ で $R$ の元が単位元になる条件は $R\subseteq N$ であり、$N$ は正規でなければならない（:ref[thm-grp-quotient]）。$R$ を含む最小の正規部分群が正規閉包であり、:ref[prop-grp-normal-closure] によりその元は関係語の共役とその逆の有限積である。$\langle R\rangle$ は一般に正規でない（例えば $F(a,b)$ の $\langle a\rangle$ は $bab^{-1}$ を含まない）。交換子群は可換化のための部分群で $R$ と無関係であり、$|X|\ge2$ なら中心は自明である（:ref[cor-grp-free-commuting]）。
:::
::::

::::quiz{answer="6" tol="0"}
$\langle a,b\mid b=a^2,\ b^3=1\rangle$ の位数を答えよ。

:::explanation
:ref[lem-grp-tietze-derived] の 2 で $b^3$ の $b$ を $a^2$ に置き換えると $\langle a,b\mid ba^{-2},\ a^6\rangle$ となり、(T4) で $b$ を削除すると $\langle a\mid a^6\rangle\cong C_6$ を得る（:ref[exm-grp-tietze-eliminate]、:ref[exm-grp-cyclic-presentation]）。従って位数は $6$ である。$b^3=1$ を書き換えずに $b$ を消すと $\Z$ になってしまう点に注意する。
:::
::::

::::quiz{answer="36" tol="0"}
階数 $2$ の自由群 $F(a,b)$ の元で、簡約語の長さがちょうど $3$ であるものはいくつあるか。

:::explanation
最初の文字は $a,a^{-1},b,b^{-1}$ の $4$ 通り、二文字目と三文字目はそれぞれ直前の文字の逆以外の $3$ 通りなので、$4\cdot3\cdot3=36$ 個の簡約語がある。異なる簡約語は異なる元を表す（:ref[lem-grp-reduction]）ので、元の個数も $36$ である。
:::
::::

::::quiz{answer="11" tol="0"}
階数 $3$ の自由群の指数 $5$ の部分群の階数を答えよ。

:::explanation
Schreier の指数公式 :eqref[eq-grp-schreier-index] により $1+d(r-1)=1+5\cdot(3-1)=11$ である。Schreier グラフには頂点が $5$ 個、正の辺が $5\cdot3=15$ 本あり、親の辺 $4$ 本を除いた $11$ 本が自由基底に対応する。
:::
::::

::::quiz
Tietze 変換（またはその組合せ）として、表示群の同型類を常に保つ操作をすべて選べ。

- [x] 既存の関係語の正規閉包に属する語を関係語に追加する
- [ ] 既存の生成元だけを用いた任意の語を関係語に追加する
- [x] 新しい生成元 $y$ を、旧生成元の語 $w$ による定義関係 $y=w$ とともに追加する
- [x] 関係語 $r$ を共役 $uru^{-1}$ に置き換える
- [ ] 生成元 $y$ とその定義関係 $y=w$ を、他の関係語に $y$ が残ったまま削除する

:::explanation
第一は (T1)、第三は (T3)、第四は :ref[lem-grp-tietze-derived] の 1 で、いずれも群を変えない（:ref[prop-grp-tietze-sound]）。第二は一般に剰余群への移行であり、例えば $\langle a\mid\ \rangle\cong\Z$ に $a^2$ を加えると $C_2$ になる。第五は (T4) の条件に反し、:ref[exm-grp-tietze-eliminate] のように $C_6$ が $\Z$ に変わることがある。
:::
::::

::::quiz
自由群について正しいものをすべて選べ。

- [x] 自由群の任意の部分群は自由群である
- [ ] $F_2$ の部分群の階数は $2$ 以下である
- [ ] $F_2$ には位数 $2$ の元がある
- [x] $F_3$ の可換化は $\Z^3$ と同型である

:::explanation
第一は Nielsen–Schreier の定理（:ref[thm-grp-free-subgroup]）である。第二は誤りで、指数 $2$ の部分群は階数 $3$ であり（:ref[exm-grp-index-two-subgroup]）、交換子群は無限階数である（:ref[exm-grp-commutator-subgroup-free]）。第三は誤りで、自由群は捩れをもたない（:ref[prop-grp-free-torsion-free]）。第四は :ref[exm-grp-free-abelianization] による。
:::
::::

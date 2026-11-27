---
title: 実数の連続性
description: 体と順序の公理および上限の公理から出発し、自然数と帰納法、上限・下限、Dedekind の切断と区間縮小法による連続性の諸形式、Archimedes の性質、有理数の稠密性、n 乗根と 10 進展開の存在を証明する。
notation:
  - symbol: '\N,\ \Z,\ \Q'
    meaning: 自然数（1 から始まる）、整数、有理数の全体
  - symbol: '\sup A'
    meaning: 集合 A の最小の上界（上限）
  - symbol: '\inf A'
    meaning: 集合 A の最大の下界（下限）
  - symbol: '\overline{\R}'
    meaning: 拡張実数 ℝ ∪ {−∞, +∞}
  - symbol: '\abs{x}'
    meaning: 実数 x の絶対値
  - symbol: '\lfloor x\rfloor'
    meaning: x 以下の最大の整数
  - symbol: '\sqrt[n]{t}'
    meaning: t ≥ 0 の非負の n 乗根
---

## この章の目標

- 体と順序の公理から四則演算と不等式の規則を導き、自然数・整数・有理数と数学的帰納法を公理の上に組み立てる。
- 上限と下限を定義し、誤差 $\eps$ を用いた特徴付けによってその値を証明する。拡張実数での約束を理解する。
- 上限の公理、Dedekind の切断の性質、区間縮小法と Archimedes の性質の組が、順序体において互いに同値であることを証明する。
- Archimedes の性質から整数部分と有理数の稠密性を導き、$n$ 乗根の存在、$\sqrt2$ の無理性、10 進展開の存在、実数の非可算性を証明する。

長さが数で表せるとすれば、一辺の長さが $1$ の正方形の対角線にも、平方が $2$ になる正の数が対応するはずである。
しかし、有理数の四則演算を何度繰り返しても、その数には届かない。
そこで本章では、数を計算する規則に加えて、集合の境界に数が存在することを公理に採る。
以後の解析学のすべての定理、すなわち極限の存在、中間値の定理、最大値の存在、積分の存在は、最終的にこの公理に帰着する。

## 実数の公理

### 体と順序の公理

:::definition[体と順序体]{#def-anl-ordered-field}
集合 $K$ に加法 $+$ と乗法 $\cdot$ という二つの演算が与えられ、すべての $a,b,c\in K$ について次の (F1)〜(F9) が成り立つとき、$K$ を :term[体]{en="field" yomi="たい"} という。

- (F1) $(a+b)+c=a+(b+c)$。
- (F2) $a+b=b+a$。
- (F3) すべての $a$ について $a+0=a$ となる元 $0\in K$ が存在する。
- (F4) 各 $a$ に対して $a+a'=0$ となる元 $a'\in K$ が存在する。
- (F5) $(ab)c=a(bc)$。
- (F6) $ab=ba$。
- (F7) すべての $a$ について $a\cdot1=a$ となる元 $1\in K$ で、$1\ne0$ を満たすものが存在する。
- (F8) 各 $a\ne0$ に対して $aa^*=1$ となる元 $a^*\in K$ が存在する。
- (F9) $a(b+c)=ab+ac$。

さらに $K$ 上の二項関係 $\le$ が、すべての $a,b,c\in K$ について次の (O1)〜(O6) を満たすとき、$(K,\le)$ を :term[順序体]{en="ordered field" yomi="じゅんじょたい"} という。

- (O1) $a\le a$。
- (O2) $a\le b$ かつ $b\le a$ ならば $a=b$。
- (O3) $a\le b$ かつ $b\le c$ ならば $a\le c$。
- (O4) $a\le b$ または $b\le a$ が成り立つ。
- (O5) $a\le b$ ならば $a+c\le b+c$。
- (O6) $0\le a$ かつ $0\le b$ ならば $0\le ab$。

$a\le b$ かつ $a\ne b$ のとき $a<b$ と書き、$b\ge a$、$b>a$ も同じ意味に用いる。$a>0$ となる元を正、$a<0$ となる元を負という。
:::

(O1)〜(O4) は $\le$ が全順序であることをいう。(O5) と (O6) は順序と演算の両立条件である。
例えば有理数全体 $\Q$ は通常の演算と大小関係について順序体である。
複素数全体 $\C$ は体であるが、どのように順序を入れても順序体にはならない（下の命題の後で示す）。

:::proposition[体の基本性質]{#prop-anl-field-basics}
体 $K$ において次が成り立つ。

1. (F3) の $0$ と (F7) の $1$ はただ一つである。各 $a$ に対する (F4) の $a'$ と、$a\ne0$ に対する (F8) の $a^*$ もただ一つである。これらをそれぞれ $-a$、$a^{-1}$ と書き、$b-a=b+(-a)$、$b/a=ba^{-1}$ と定める。
2. $a+c=b+c$ ならば $a=b$ である。$ac=bc$ かつ $c\ne0$ ならば $a=b$ である。
3. $a\cdot0=0$ である。
4. $-(-a)=a$、$(-1)a=-a$、$(-a)b=-(ab)$、$(-a)(-b)=ab$ である。
5. $ab=0$ ならば $a=0$ または $b=0$ である。
:::

:::proof
1. $0'$ も (F3) を満たせば、(F2) により $0'=0'+0=0+0'=0$ である。$1$ についても同様である。
$a'$ と $a''$ がともに (F4) を満たせば
$a''=a''+0=a''+(a+a')=(a''+a)+a'=0+a'=a'$ となる。$a^*$ についても乗法で同じ計算をすればよい。
2. 両辺に $-c$ を加えると、(F1) により $a=a+(c+(-c))=(a+c)+(-c)=(b+c)+(-c)=b$ である。乗法の場合は $c^{-1}$ を掛ける。
3. $a\cdot0=a(0+0)=a\cdot0+a\cdot0$ であるから、$0+a\cdot0=a\cdot0+a\cdot0$ に 2 を用いて $0=a\cdot0$ を得る。
4. $(-a)+a=0$ は、$a$ が $-a$ の加法の逆元であることを意味するから、1 の一意性により $-(-a)=a$ である。
次に $a+(-1)a=(1+(-1))a=0\cdot a=0$ より $(-1)a=-a$ である。
これを用いると $(-a)b=((-1)a)b=(-1)(ab)=-(ab)$、さらに $(-a)(-b)=-(a(-b))=-(-(ab))=ab$ となる。
5. $a\ne0$ ならば $b=a^{-1}(ab)=a^{-1}\cdot0=0$ である。
:::

:::proposition[順序体の基本性質]{#prop-anl-order-basics}
順序体 $K$ の元 $a,b,c,d$ について次が成り立つ。

1. $a<b$、$a=b$、$b<a$ のうちちょうど一つが成り立つ。
2. $a\le b\iff0\le b-a$ である。特に $a\le0\iff0\le-a$ である。
3. $a<b$ ならば $a+c<b+c$ である。$a\le b$ かつ $c\le d$ ならば $a+c\le b+d$ であり、どちらかが狭義なら和も狭義である。
4. $a<b$ かつ $c>0$ ならば $ac<bc$ である。$a<b$ かつ $c<0$ ならば $bc<ac$ である。
5. $a\ne0$ ならば $a^2>0$ である。特に $1>0$ である。
6. $a>0$ ならば $a^{-1}>0$ である。$0<a<b$ ならば $0<b^{-1}<a^{-1}$ である。
7. $2=1+1$ とおくと $2>0$ であり、$a<b$ ならば $a<(a+b)/2<b$ である。
:::

:::proof
1. (O4) により $a\le b$ か $b\le a$ の少なくとも一方が成り立ち、両方なら (O2) により $a=b$ である。
従って $a\ne b$ ならば $a<b$ と $b<a$ のちょうど一方が成り立ち、$a=b$ ならどちらも成り立たない。
2. $a\le b$ の両辺に (O5) で $-a$ を加えれば $0\le b-a$、逆に $0\le b-a$ に $a$ を加えれば $a\le b$ である。後半は $b=0$ の場合である。
3. (O5) で $a+c\le b+c$ であり、等号なら :ref[prop-anl-field-basics] の 2 により $a=b$ となるので狭義である。
後半は (O5) を二度用いて $a+c\le b+c\le b+d$ とし、(O3) を使う。狭義の場合も前半から従う。
4. 2 により $0<b-a$ である。(O6) から $0\le(b-a)c$ であり、:ref[prop-anl-field-basics] の 5 により $(b-a)c\ne0$ である。
従って $0<bc-ac$、すなわち $ac<bc$ である。$c<0$ なら 2 により $-c>0$ なので $a(-c)<b(-c)$ を得て、両辺に $ac+bc$ を加えれば $bc<ac$ となる。
5. $a>0$ なら 4 を $0<a$ と $c=a$ に用いて $0=0\cdot a<a\cdot a$ である。
$a<0$ なら $-a>0$ であり、:ref[prop-anl-field-basics] の 4 より $a^2=(-a)^2>0$ である。$1=1^2$ かつ $1\ne0$ なので $1>0$ である。
6. $a>0$ とする。$a^{-1}=0$ なら $1=aa^{-1}=0$ となるので $a^{-1}\ne0$ である。
もし $a^{-1}<0$ なら、4 を $c=a>0$ に用いて $a^{-1}a<0\cdot a$、すなわち $1<0$ となり 5 に反する。
後半では $a^{-1},b^{-1}>0$ より $a^{-1}b^{-1}>0$ であり、$a<b$ にこれを掛けると $b^{-1}<a^{-1}$ である。
7. 3 により $2=1+1>0+1=1>0$ である。$a<b$ なら 3 により $2a=a+a<a+b<b+b=2b$ であり、$2^{-1}>0$ を掛ければよい。
:::

$\C$ が順序体にならないことは 5 から従う。$i\ne0$ なら $i^2>0$ でなければならないが、$i^2=-1$ と $1>0$ から $-1<0$ である。

### 絶対値と区間

:::definition[絶対値と区間]
順序体の元 $x$ の :term[絶対値]{en="absolute value" yomi="ぜったいち"} を、$x\ge0$ なら $\abs{x}=x$、$x<0$ なら $\abs{x}=-x$ と定める。これは $\max\{x,-x\}$ に等しい。
$a\le b$ に対して

$$
[a,b]=\{x\mid a\le x\le b\},\quad (a,b)=\{x\mid a<x<b\},\quad [a,b)=\{x\mid a\le x<b\},\quad (a,b]=\{x\mid a<x\le b\}
$$

をそれぞれ :term[閉区間]{en="closed interval" yomi="へいくかん"}、:term[開区間]{en="open interval" yomi="かいくかん"}、半開区間という。
また $(a,\infty)=\{x\mid x>a\}$、$[a,\infty)=\{x\mid x\ge a\}$、$(-\infty,b)$、$(-\infty,b]$、$(-\infty,\infty)=K$ を無限区間という。
記号 $\infty$ はここでは単に「片側に制限がない」ことを表し、$K$ の元ではない。
:::

:::proposition[絶対値の性質と三角不等式]{#prop-anl-triangle}
順序体の任意の元 $x,y$ と $c\ge0$ について次が成り立つ。

1. $\abs{x}\ge0$ であり、$\abs{x}=0\iff x=0$ である。また $\abs{-x}=\abs{x}$、$-\abs{x}\le x\le\abs{x}$ である。
2. $\abs{xy}=\abs{x}\abs{y}$ である。
3. $\abs{x}\le c\iff -c\le x\le c$、$\abs{x}<c\iff -c<x<c$ である。
4. $\abs{x+y}\le\abs{x}+\abs{y}$ および $\bigl|\abs{x}-\abs{y}\bigr|\le\abs{x-y}$ が成り立つ。
:::

:::proof
1. $x\ge0$ なら $\abs{x}=x\ge0$、$x<0$ なら :ref[prop-anl-order-basics] の 2 により $\abs{x}=-x>0$ である。従って $\abs{x}=0$ となるのは $x=0$ のときに限る。
$\abs{-x}=\max\{-x,x\}=\abs{x}$ であり、$x$ と $-x$ はともに最大値 $\abs{x}$ 以下なので $-\abs{x}\le x\le\abs{x}$ である。
2. $\abs{x}=\sigma x$、$\abs{y}=\tau y$ となる $\sigma,\tau\in\{1,-1\}$ を選ぶ。:ref[prop-anl-field-basics] の 4 より $\abs{x}\abs{y}=(\sigma\tau)xy$ であり、$\sigma\tau\in\{1,-1\}$ である。
左辺は (O6) により非負であるから、これは $xy$ と $-xy$ のうち非負の方、すなわち $\abs{xy}$ に等しい。
3. $\abs{x}\le c$ は $x\le c$ かつ $-x\le c$ と同値であり、後者は $-c\le x$ と同値である。狭義の場合も同じである。
4. 1 より $-\abs{x}\le x\le\abs{x}$、$-\abs{y}\le y\le\abs{y}$ であり、辺々加えると $-(\abs{x}+\abs{y})\le x+y\le\abs{x}+\abs{y}$ となる。3 より最初の不等式を得る。
次に $x=(x-y)+y$ に適用すれば $\abs{x}-\abs{y}\le\abs{x-y}$ であり、$x,y$ を交換すれば $\abs{y}-\abs{x}\le\abs{x-y}$ も成り立つ。3 により後半を得る。
:::

3 により、$\eps>0$ に対して $\abs{x-c}<\eps$ は $c-\eps<x<c+\eps$ と同じ条件である。
以後、近似の誤差を述べるときはこの二つの書き方を使い分ける。

### 自然数・整数・有理数

順序体 $K$ の中で $1$ を順に加えると $0<1<1+1<1+1+1<\cdots$ となる（:ref[prop-anl-order-basics] の 3, 5）。
「このようにして得られる元全体」を厳密に定義するために、$1$ を含み $1$ を加える操作で閉じた集合のうち最小のものを考える。

:::definition[帰納的集合と自然数]{#def-anl-inductive}
順序体 $K$ の部分集合 $S$ が $1\in S$ を満たし、$x\in S$ ならば $x+1\in S$ となるとき、$S$ を :term[帰納的集合]{en="inductive set" yomi="きのうてきしゅうごう"} という。
$K$ の帰納的部分集合すべての共通部分を $\N$ と書き、その元を自然数という。さらに

$$
\Z=\N\cup\{0\}\cup\{-n\mid n\in\N\},\qquad \Q=\{m/n\mid m\in\Z,\ n\in\N\}
$$

とおき、その元をそれぞれ整数、有理数という。$\Z_{\ge0}=\N\cup\{0\}$ とも書く。
:::

$K$ 自身が帰納的なので、共通部分をとる族は空でない。$1$ はすべての帰納的集合に属し、$x$ がすべての帰納的集合に属せば $x+1$ もそうであるから、$\N$ 自身も帰納的であり、どの帰納的集合にも含まれる。
自然数を集合論から構成する方法は [自然数・整数・有理数](/math/sets-logic/natural-numbers/) で扱う。本章の定義は、与えられた順序体の中に自然数の複製を見つけるものである。

:::theorem[数学的帰納法]{#thm-anl-induction}
$S\subset\N$ が $1\in S$ を満たし、$n\in S$ ならば $n+1\in S$ となるとき、$S=\N$ である。
:::

:::proof
仮定により $S$ は帰納的であるから、$\N$ の定義により $\N\subset S$ である。$S\subset\N$ と合わせて $S=\N$ を得る。
:::

自然数に関する主張 $P(n)$ を示すには、$S=\{n\in\N\mid P(n)\text{ が成り立つ}\}$ にこの定理を適用する。これが :term[数学的帰納法]{en="mathematical induction" yomi="すうがくてききのうほう"} である。

:::proposition[自然数の基本性質]{#prop-anl-natural-basics}
1. すべての $n\in\N$ について $n\ge1$ である。
2. $n\in\N$、$n\ne1$ ならば $n-1\in\N$ である。
3. $m,n\in\N$ ならば $m+n\in\N$、$mn\in\N$ である。
4. $m,n\in\N$、$m<n$ ならば $n-m\in\N$ である。
5. $m,n\in\N$、$m>n$ ならば $m\ge n+1$ である。すなわち $n<m<n+1$ となる自然数 $m$ は存在しない。
:::

:::proof
1. $\{x\in K\mid x\ge1\}$ は $1$ を含み、$x\ge1$ なら $x+1>x\ge1$ なので帰納的である。従ってこの集合は $\N$ を含む。
2. $S=\{1\}\cup\{n\in\N\mid n-1\in\N\}$ とおく。$1\in S$ であり、$n\in S$ なら $(n+1)-1=n\in\N$ なので $n+1\in S$ である。:ref[thm-anl-induction] により $S=\N$ である。
3. $m$ を固定し、$S=\{n\in\N\mid m+n\in\N\}$ とおく。$m+1\in\N$ より $1\in S$ であり、$n\in S$ なら $m+(n+1)=(m+n)+1\in\N$ である。従って $S=\N$ である。
積については $T=\{n\in\N\mid mn\in\N\}$ とおくと、$m\cdot1=m\in\N$ であり、$n\in T$ なら和の結果から $m(n+1)=mn+m\in\N$ である。
4. $S=\{m\in\N\mid n>m\text{ を満たすすべての }n\in\N\text{ について }n-m\in\N\}$ とおく。
$n>1$ なら $n\ne1$ なので 2 から $n-1\in\N$ であり、$1\in S$ である。
$m\in S$ とし、$n>m+1$ とする。$n>m$ なので $k=n-m\in\N$ であり、$k>1$ だから 2 により $k-1=n-(m+1)\in\N$ である。従って $m+1\in S$ である。
5. 4 により $m-n\in\N$ であり、1 により $m-n\ge1$ である。
:::

:::theorem[自然数の整列性]{#thm-anl-well-ordering}
$\N$ の空でない部分集合は最小元をもつ。
:::

:::proof
$S\subset\N$ が最小元をもたないと仮定し、$S=\varnothing$ を示す。
$T=\{n\in\N\mid k\le n\text{ を満たすすべての }k\in\N\text{ について }k\notin S\}$ とおく。

$1\in S$ ならば、:ref[prop-anl-natural-basics] の 1 により $1$ は $S$ の最小元となって仮定に反する。$k\le1$ となる自然数は $k=1$ だけなので $1\in T$ である。

$n\in T$ とする。$S$ の元 $k$ はすべて $k>n$ を満たすので、:ref[prop-anl-natural-basics] の 5 により $k\ge n+1$ である。
従って $n+1\in S$ ならば $n+1$ は $S$ の最小元になり、仮定に反する。よって $n+1\notin S$ である。
$k\le n+1$ となる自然数は、再び同じ性質により $k\le n$ または $k=n+1$ であるから、$n+1\in T$ を得る。

:ref[thm-anl-induction] により $T=\N$ である。各 $n\in S$ は $n\in T$ により $n\notin S$ となるので、$S=\varnothing$ である。
:::

整列性からは「:term[完全帰納法]{en="strong induction" yomi="かんぜんきのうほう"}」も従う。すなわち $S\subset\N$ が「$k<n$ となるすべての自然数 $k$ が $S$ に属するならば $n\in S$」を満たせば $S=\N$ である。
実際 $\N\setminus S$ が空でなければ、その最小元 $n$ より小さい自然数はすべて $S$ に属するので $n\in S$ となり、矛盾する。

:::proposition[整数と有理数]{#prop-anl-integers}
1. $\Z$ は加法、減法、乗法について閉じている。$\Q$ は $K$ の部分体であり、$K$ の演算と順序の制限によって順序体となる。
2. 正の整数全体は $\N$ に等しい。$k\in\Z$ ならば $k<j<k+1$ を満たす $j\in\Z$ は存在しない。
3. $\Z$ の空でない部分集合 $S$ が、ある整数 $l$ によって下から抑えられる（すべての $s\in S$ で $l\le s$）ならば、$S$ は最小元をもつ。ある整数によって上から抑えられるならば最大元をもつ。
:::

:::proof
1. まず $m,n\in\N$ に対して $m-n\in\Z$ である。実際、:ref[prop-anl-natural-basics] の 4 により、$m>n$ なら $m-n\in\N$、$m<n$ なら $m-n=-(n-m)$ で $n-m\in\N$、$m=n$ なら $0$ である。
逆にすべての整数は自然数の差で書ける。$n=(n+1)-1$、$0=1-1$、$-n=1-(n+1)$ だからである。
そこで $x=m-n$、$y=m'-n'$ と書くと

$$
x+y=(m+m')-(n+n'),\quad x-y=(m+n')-(n+m'),\quad xy=(mm'+nn')-(mn'+nm')
$$

であり、:ref[prop-anl-natural-basics] の 3 により各括弧は自然数であるから、いずれも整数である。
次に $r=m/n$、$r'=m'/n'$ ならば

$$
r\pm r'=\frac{mn'\pm m'n}{nn'},\qquad rr'=\frac{mm'}{nn'}
$$

であり $nn'\in\N$ なので、$\Q$ は加法・減法・乗法で閉じる。$r\ne0$ なら $m\ne0$ であり、$m\in\N$ のとき $r^{-1}=n/m$、$-m\in\N$ のとき $r^{-1}=(-n)/(-m)$ が $\Q$ に属する。
$0,1\in\Q$ であり、体と順序体の公理は $K$ で成り立つので $\Q$ でも成り立つ。
2. 整数は自然数、$0$、または負の数 $-n$（$n\in\N$）のいずれかであるから、正の整数は自然数である。
$k<j<k+1$ となる整数 $j$ があれば、$j-k$ は正の整数なので自然数であり、しかも $j-k<1$ となって :ref[prop-anl-natural-basics] の 1 に反する。
3. $S'=\{s-l+1\mid s\in S\}$ は 1 より $\Z$ に含まれ、その元は $1$ 以上であるから、2 により $S'\subset\N$ である。
:ref[thm-anl-well-ordering] により $S'$ は最小元 $s_0-l+1$ をもち、$s_0$ が $S$ の最小元である。上から抑えられる場合は $\{-s\mid s\in S\}$ に適用すればよい。
:::

自然数 $n$ に対するべき $a^n$ は $a^1=a$、$a^{n+1}=a^na$ で帰納的に定め、$a\ne0$ のとき $a^0=1$、$a^{-n}=(a^{-1})^n$ とする（便宜上 $0^0=1$ とも約束する）。
指数法則 $a^{m+n}=a^ma^n$、$(ab)^n=a^nb^n$、$(a^m)^n=a^{mn}$ は $n$ についての帰納法で示される。例えば $a^{m+(n+1)}=a^{m+n}a=a^ma^na=a^ma^{n+1}$、$(ab)^{n+1}=a^nb^nab=a^{n+1}b^{n+1}$、$(a^m)^{n+1}=a^{mn}a^m=a^{mn+m}$ である。
有限和 $\sum_{k=1}^n a_k$ と階乗 $n!$（$0!=1$、$(n+1)!=n!\,(n+1)$）も同様に帰納的に定める。

:::proposition[Bernoulli の不等式]{#prop-anl-bernoulli}
順序体において、$h\ge-1$ と $n\in\N$ に対し $(1+h)^n\ge1+nh$ が成り立つ。
:::

:::proof
$n=1$ では等号である。$(1+h)^n\ge1+nh$ を仮定すると、$1+h\ge0$ を両辺に掛けて

$$
(1+h)^{n+1}\ge(1+nh)(1+h)=1+(n+1)h+nh^2\ge1+(n+1)h
$$

を得る。:ref[thm-anl-induction] により結論が従う。
:::

特に $h=1$ とすると $2^n\ge1+n$ であり、$n\ge2$ で $2^{n-1}\ge1+(n-1)=n$ となる。$n=1$ でも $2^0=1$ なので、すべての $n$ で $2^{n-1}\ge n$ である。
これは二分法で区間の長さを評価するときに繰り返し使う。

:::proposition[二項定理]{#prop-anl-binomial}
$0\le k\le n$ に対して $\binom nk=\dfrac{n!}{k!\,(n-k)!}$ とおく。体の元 $a,b$ と $n\in\N$ について

$$
(a+b)^n=\sum_{k=0}^n\binom nk a^{n-k}b^k .
$$

また $\binom nk$ はすべて自然数である。
:::

:::proof
$1\le k\le n$ について、通分により

$$
\binom n{k-1}+\binom nk=\frac{n!\,\{k+(n-k+1)\}}{k!\,(n-k+1)!}=\binom{n+1}k
\label{eq-anl-pascal}
$$

が成り立つ（Pascal の関係）。$\binom n0=\binom nn=1$ と :eqref[eq-anl-pascal] から、$n$ についての帰納法で $\binom nk\in\N$ がわかる。
展開式は $n=1$ で明らかである。$n$ で成り立つとすると

$$
\begin{aligned}
(a+b)^{n+1}
&=\sum_{k=0}^n\binom nk a^{n+1-k}b^k+\sum_{k=0}^n\binom nk a^{n-k}b^{k+1}\\
&=a^{n+1}+\sum_{k=1}^n\left\{\binom nk+\binom n{k-1}\right\}a^{n+1-k}b^k+b^{n+1}
\end{aligned}
$$

であり、:eqref[eq-anl-pascal] により $n+1$ の場合の式になる。
:::

### 上限の公理

$A\subset K$ とする。すべての $a\in A$ に対して $a\le u$ となる $u\in K$ を $A$ の
:term[上界]{en="upper bound" yomi="じょうかい"} といい、上界が存在するとき $A$ は上に有界であるという。
すべての $a\in A$ に対して $l\le a$ となる $l\in K$ を
:term[下界]{en="lower bound" yomi="かかい"} といい、下界が存在するとき下に有界であるという。上にも下にも有界な集合を有界という。

上界の中に最小の元があれば、それを $A$ の :term[上限]{en="supremum" yomi="じょうげん"} といい、$\sup A$ と書く。
実数に特有の条件を、次の形で公理として採用する。

:::axiom[上限の公理]{#ax-anl-supremum}
実数体 $\R$ は順序体であり、その空でない上に有界な部分集合は、必ず $\R$ の中に上限をもつ。
この性質を実数の :term[完備性]{en="completeness" yomi="かんびせい"}、あるいは実数の連続性という。
:::

:::remark[完備順序体の存在と一意性]
上限の公理を満たす順序体は存在し、そのような順序体が二つあれば、その間に演算と順序を保つ全単射がただ一つ存在する。
存在は Dedekind の切断や Cauchy 列による構成で、一意性は有理数の稠密性を用いて証明される。いずれも集合論の主題であり、[実数と複素数の構成](/math/sets-logic/number-systems/) で証明する。
本書の解析学はこの公理系だけを出発点とし、実数の具体的な構成には依存しない。以下 $\N,\Z,\Q$ は $K=\R$ として上の定義で定めたものとする。
:::

順序と連続性から解析を組み立てる扱いは :cite[takagi-kaiseki] にもある。
$\Q$ は順序体であるが上限の公理を満たさない。このことは $\sqrt2$ の無理性を示した後、本章の最後の節で確かめる。

## 上限と下限

### 最大値との違い

:::definition[上限と下限の順序による特徴付け]{#def-anl-bounds}
空でない $A\subset\R$ に対し、$\alpha=\sup A$ は次の二条件を意味する。

1. 任意の $a\in A$ に対して $a\le\alpha$ である。
2. $u$ が $A$ の上界ならば $\alpha\le u$ である。

下界の中の最大の元を $A$ の :term[下限]{en="infimum" yomi="かげん"} といい、$\inf A$ と書く。
従って $\beta=\inf A$ とは、各 $a\in A$ について $\beta\le a$ であり、各下界 $l$ について $l\le\beta$ となることである。
:::

$\alpha,\alpha'$ が同じ集合の上限の条件を満たすとする。
どちらも上界なので、第 2 条件を両方に適用すると $\alpha\le\alpha'$ かつ $\alpha'\le\alpha$ となる。
従って上限は存在すれば一意であり、下限も同様である。

上界であると同時に $A$ の元でもある数を $A$ の**最大値**といい、$\max A$ と書く。
最大値があれば、どの上界もその元以上なので最大値は上限に一致する。
最小値 $\min A$ とは $A$ に属する下界のことであり、存在すれば下限に一致する。
空でない有限集合は常に最大値と最小値をもつ。元の個数についての帰納法により、一点集合 $\{a\}$ の最大値は $a$ であり、$A$ の最大値が $m$ なら $A\cup\{b\}$ の最大値は $\max\{m,b\}$ だからである（最小値も同様）。

:::example[端点が集合に属するか]
$A=(0,1]$ の最大値は $1$ である。一方、$0$ は下界であるが $A$ に属さない。
どの $x\in A$ についても $x/2\in A$ かつ $x/2<x$ だから、最小値は存在しない。
それでも下限は $0$ である。実際、$l>0$ を任意にとると $x=\min\{l/2,1\}$ は $A$ の元で $x<l$ を満たすので、正の数は下界でない。
:::

:::proposition[下限の存在]{#prop-anl-infimum}
空でない下に有界な $A\subset\R$ は下限をもち、$-A=\{-a\mid a\in A\}$ とおくと

$$
\inf A=-\sup(-A).
$$
:::

:::proof
$l$ が $A$ の下界であるための必要十分条件は、$-l$ が $-A$ の上界であることである。
従って $-A$ は空でなく上に有界であり、:ref[ax-anl-supremum] より $s=\sup(-A)$ が存在する。
$s$ が上界であることから $-s$ は $A$ の下界である。
さらに任意の下界 $l$ について $s\le-l$、従って $l\le-s$ であるから、$-s$ が最大の下界である。
:::

### 誤差を指定して上限に近づく

次の定理は上限の存在を主張するものではなく、与えられた元が上限であるかどうかを判定する。
証明は順序体の公理しか使わないので、任意の順序体 $K$ で成り立つ。この一般性は次節で上限の公理と他の形式の同値性を示すときに用いる。

:::theorem[上限・下限の ε による特徴付け]{#thm-anl-bound-epsilon}
順序体 $K$ の空でない部分集合 $A$ と $s\in K$ に対して

$$
s=\sup A
\quad\Longleftrightarrow\quad
\begin{cases}
a\le s & (\text{すべての }a\in A),\\
\text{ある }a\in A\text{ が }s-\eps<a & (\text{各 }\eps>0)
\end{cases}
\label{eq-anl-sup-epsilon}
$$

である。同様に $t=\inf A$ であるための必要十分条件は、$t$ が下界であり、各 $\eps>0$ に対して $a<t+\eps$ となる $a\in A$ があることである。
:::

:::proof
$s=\sup A$ とする。$\eps>0$ なら $s-\eps<s$ は上界になれないので、ある $a\in A$ が $s-\eps<a$ を満たす。
逆に右辺の二条件を仮定する。もし $u<s$ が上界なら、$\eps=s-u>0$ を用いて $u=s-\eps<a$ となる $a\in A$ が得られ、$u$ が上界であることに反する。
よって $s$ は最小の上界である。下限については $-A$ と $-t$ に同じ議論を適用し、符号を戻せばよい。
:::

式 :eqref[eq-anl-sup-epsilon] で $a$ は $\eps$ を決めた後で選ぶ。
すべての $\eps$ に同じ $a$ を使えるのは、$s$ 自身が $A$ に属するときに限る。

:::example[開区間の上限]
$a<b$ として $A=(a,b)$ を考える。$b$ は上界であり、任意の $\eps>0$ に対し

$$
d=\tfrac12\min\{\eps,b-a\},\qquad x=b-d
$$

と選ぶと $a<x<b$ かつ $b-\eps<x$ である。従って $\sup A=b$ である。
同じ $d$ に対して $x=a+d$ を選ぶと $\inf A=a$ も示せる。
端点はいずれも $A$ に属さないので、どちらの値も $A$ の元では達しない。
:::

### 上限の計算

:::proposition[集合の和の上限]{#prop-anl-sup-sum}
$A,B\subset\R$ を空でない上に有界な集合とする。
$A+B=\{a+b\mid a\in A,\ b\in B\}$ とおくと $\sup(A+B)=\sup A+\sup B$ である。
:::

:::proof
$s=\sup A,\ t=\sup B$ とおく。各 $a+b$ に対して $a+b\le s+t$ であり、$s+t$ は上界である。
これが最小であることを示すため、$\eps>0$ をとる。
:ref[thm-anl-bound-epsilon] を誤差 $\eps/2$ で二度用い、

$$
s-\eps/2<a\le s,\qquad t-\eps/2<b\le t
$$

を満たす $a\in A,\ b\in B$ を選ぶ。すると $a+b\in A+B$ かつ $s+t-\eps<a+b$ となるから、同じ定理で結論を得る。
:::

:::proposition[上限の基本的な性質]{#prop-anl-sup-properties}
$A,B\subset\R$ は空でないとする。

1. $A\subset B$ かつ $B$ が上に有界ならば $\sup A\le\sup B$ である。$B$ が下に有界ならば $\inf B\le\inf A$ である。
2. $A$ が上に有界で $c>0$ ならば $\sup(cA)=c\sup A$ である。ここで $cA=\{ca\mid a\in A\}$ とする。
3. $A,B$ が上に有界ならば $\sup(A\cup B)=\max\{\sup A,\sup B\}$ である。
4. すべての $a\in A$、$b\in B$ について $a\le b$ ならば、$A$ は上に、$B$ は下に有界であり、$\sup A\le\inf B$ である。
5. $A$ が有界ならば $\inf A\le\sup A$ であり、等号は $A$ が一点集合のときに限る。
:::

:::proof
1. $\sup B$ は $B$ の上界であるから $A$ の上界でもあり、上限の最小性から $\sup A\le\sup B$ である。下限も同様である。
2. $a\le\sup A$ に $c>0$ を掛ければ $c\sup A$ は $cA$ の上界である。逆に $u$ が $cA$ の上界なら、$ca\le u$ に $c^{-1}>0$ を掛けて $u/c$ が $A$ の上界になるので、$\sup A\le u/c$、すなわち $c\sup A\le u$ である。
3. $m=\max\{\sup A,\sup B\}$ は $A\cup B$ の上界である。$A\cup B$ の上界 $u$ は $A$ と $B$ の上界でもあるから $u\ge\sup A$ かつ $u\ge\sup B$、従って $u\ge m$ である。
4. 各 $b\in B$ は $A$ の上界なので $\sup A\le b$ である。従って $\sup A$ は $B$ の下界であり、$\sup A\le\inf B$ となる。
5. $a\in A$ を一つとれば $\inf A\le a\le\sup A$ である。等号が成り立てば、すべての $a\in A$ が $\inf A\le a\le\sup A=\inf A$ を満たすので $A$ は一点である。
:::

関数 $f:D\to\R$ に対しては $\sup_{x\in D}f(x)=\sup f(D)$、$\inf_{x\in D}f(x)=\inf f(D)$ と書く。$f(D)$ が上に有界なとき $f$ は上に有界であるという。
$f,g$ が上に有界なら、各 $x$ で $f(x)+g(x)\le\sup f+\sup g$ なので

$$
\sup_{x\in D}\bigl(f(x)+g(x)\bigr)\le\sup_{x\in D}f(x)+\sup_{x\in D}g(x)
$$

である。:ref[prop-anl-sup-sum] と異なり、ここで等号は一般に成り立たない。
$D=[0,1]$ で $f(x)=x$、$g(x)=1-x$ とすると右辺は $2$ であるが、$f+g$ は常に $1$ である。
集合の和 $A+B$ では $a$ と $b$ を独立に選べたのに対し、関数の和では同じ $x$ で二つの値を同時に選ばなければならないからである。

### 拡張実数

上に有界でない集合には上限がない。しかし極限を扱う際には「上限が無限大である」と言えると便利である。

:::definition[拡張実数]{#def-anl-extended-reals}
$\R$ に二つの記号 $-\infty,+\infty$ を付け加えた集合を $\overline{\R}=\R\cup\{-\infty,+\infty\}$ と書き、:term[拡張実数]{en="extended real numbers" yomi="かくちょうじっすう"} という。
順序は $\R$ の順序に、すべての $x\in\R$ について $-\infty<x<+\infty$ を加えて定める。
$A\subset\overline{\R}$ の上限・下限は、この全順序における最小の上界・最大の下界と定める。
:::

:::proposition[拡張実数における上限]{#prop-anl-extended-sup}
$\overline{\R}$ の任意の部分集合 $E$ は $\overline{\R}$ の中に上限と下限をもつ。特に $A\subset\R$ について

$$
\sup A=
\begin{cases}
\text{$\R$ における上限} & (A\ne\varnothing\text{ が上に有界}),\\
+\infty & (A\text{ が上に有界でない}),\\
-\infty & (A=\varnothing)
\end{cases}
$$

であり、下限についても符号を逆にした等式が成り立つ。
:::

:::proof
$+\infty\in E$ なら $+\infty$ が最大元であり、上限である。以下 $+\infty\notin E$ とし、$A=E\cap\R$ とおく。
$A=\varnothing$ なら $E\subset\{-\infty\}$ であり、すべての元が上界になるので最小の上界は $-\infty$ である。
$A$ が空でなく $\R$ で上に有界なら、$s=\sup A$ は $E$ の上界であり（$-\infty\le s$）、$E$ の上界で $s$ より小さいものは $A$ の上界になるので存在しない。
$A$ が $\R$ で上に有界でなければ、実数は $E$ の上界にならず、上界は $+\infty$ だけである。下限は $-E$ に同じ議論を適用すればよい。
:::

この約束の下で、$\sup A=+\infty$ は「任意の $M\in\R$ に対し $a>M$ となる $a\in A$ がある」ことと同値である。
$\sup\varnothing=-\infty$、$\inf\varnothing=+\infty$ となり、:ref[prop-anl-sup-properties] の 4 の不等式は空集合を含めても成り立つ。
演算については $x\in\R$ に対して $x+(+\infty)=+\infty$、$x+(-\infty)=-\infty$、$x>0$ のとき $x\cdot(\pm\infty)=\pm\infty$ と約束するが、$(+\infty)+(-\infty)$ や $0\cdot(+\infty)$ は定義しない。
$\overline{\R}$ は体ではないので、等式変形には注意が要る。上極限・下極限を扱う [数列と級数](/math/analysis/sequences-series/) でこの約束を用いる。

## 実数の連続性の諸形式

上限の公理は「上に有界な集合の上界全体には最小元がある」と述べる。
同じ内容を、数直線の切断、区間の列、数列の言葉で言い換えることができる。どの形式を用いるかは問題によって使い分ける。
本節と次節では、上限の公理を仮定しない一般の順序体 $K$ で議論する場合を明示する。

### Dedekind の切断

数直線を一点で左右に切り分けたとき、その切り口には必ず数がある。これを公理の形にしたのが Dedekind の切断である。

:::definition[切断]{#def-anl-dedekind-cut}
順序体 $K$ の部分集合の組 $(L,U)$ が次の三条件を満たすとき、$K$ の :term[切断]{en="Dedekind cut" yomi="せつだん"} という。

1. $L\ne\varnothing$、$U\ne\varnothing$ である。
2. $L\cup U=K$、$L\cap U=\varnothing$ である。
3. $x\in L$、$y\in U$ ならば $x<y$ である。

どの切断についても「$L$ が最大元をもつか、$U$ が最小元をもつ」とき、$K$ は **Dedekind の性質** をもつという。
:::

切断で $L$ の最大元 $a$ と $U$ の最小元 $b$ が同時に存在することはない。実際 $a<b$ であり、:ref[prop-anl-order-basics] の 7 により $c=(a+b)/2$ は $a<c<b$ を満たすので、$c$ は $L$ にも $U$ にも属さず $L\cup U=K$ に反する。
従って Dedekind の性質の下では、各切断に対し、ちょうど一方の端に境目の数がある。

:::theorem[上限の公理と Dedekind の性質]{#thm-anl-dedekind}
順序体 $K$ について、次の二条件は同値である。

1. $K$ の空でない上に有界な部分集合は、$K$ の中に上限をもつ。
2. $K$ は Dedekind の性質をもつ。
:::

:::proof
**1 ⇒ 2。** 切断 $(L,U)$ をとる。$L$ は空でなく、$U$ の任意の元を上界にもつので、上限 $s=\sup L$ が存在する。
$s\in L$ ならば $s$ は $L$ の最大元である。
$s\in U$ ならば、任意の $y\in U$ は $L$ の上界なので $s\le y$ であり、$s$ は $U$ の最小元である。

**2 ⇒ 1。** 空でない上に有界な $A\subset K$ をとり、$U$ を $A$ の上界全体、$L=K\setminus U$ とする。
$U$ は仮定により空でない。$a\in A$ をとれば $a-1<a$ なので $a-1$ は上界でなく、$L$ も空でない。
$x\in L$、$y\in U$ ならば、$x$ は上界でないので $x<a$ となる $a\in A$ があり、$x<a\le y$ である。よって $(L,U)$ は切断である。

$L$ が最大元 $m$ をもつと仮定する。$m$ は上界でないので $m<a$ となる $a\in A$ がある。
すると $m<(m+a)/2<a$ であり、$(m+a)/2$ も上界でないから $L$ に属する。これは $m$ の最大性に反する。
従って Dedekind の性質により $U$ が最小元をもち、それは最小の上界、すなわち $\sup A$ である。
:::

$\Q$ では、$L=\{q\in\Q\mid q\le0\text{ または }q^2<2\}$ と $U=\Q\setminus L$ が切断になり、$L$ に最大元も $U$ に最小元もない（演習問題）。
切り口にあるべき数 $\sqrt2$ が $\Q$ に欠けているのである。

### 入れ子になった閉区間

$[a_{n+1},b_{n+1}]\subset[a_n,b_n]$ がすべての $n\in\N$ で成り立つとき、閉区間の列は**入れ子**になっているという。
どの段階でも区間を空にしなければ、無限に続けた後にも点が残るだろうか。
共通点を探すには、すべての左端以上で、すべての右端以下の数を求めればよい。

:::theorem[閉区間の共通点]{#thm-anl-nested-intervals}
空でない有界閉区間 $I_n=[a_n,b_n]$ が入れ子ならば、すべての $I_n$ に属する実数が存在する。
さらに各 $\eps>0$ に対してある $n$ が $b_n-a_n<\eps$ を満たすなら、共通点はただ一つである。
この後半を含む主張を :term[区間縮小法]{en="nested interval principle" yomi="くかんしゅくしょうほう"} という。
:::

:::proof
任意の $j,k$ について $a_j\le b_k$ である。実際、$j\le k$ なら $a_j\le a_k\le b_k$ であり、$j>k$ なら $a_j\le b_j\le b_k$ である。
従って各 $b_k$ は $A=\{a_j\mid j\in\N\}$ の上界である。
:ref[ax-anl-supremum] により上限 $s=\sup A$ をとると、$s$ は上界なので $a_k\le s$、最小の上界なので $s\le b_k$ であり、すべての $k$ について $s\in I_k$ を得る。

もし二つの共通点 $x<y$ があれば、すべての $n$ について $0<y-x\le b_n-a_n$ である。
ところが追加条件を $\eps=y-x$ に適用すると、これに反する $n$ が得られる。
:::

この証明は上限の公理だけを用いるので、上限の公理を満たす任意の順序体で成り立つ。

:::example[端点を除いた場合]
開区間 $J_n=(0,1/n)$ も入れ子であるが、共通点はない。
実際、共通点の候補 $x$ は正でなければならず、次節で示す Archimedes の性質により $n>1/x$ と選ぶと $x\notin J_n$ となる。
一方、閉区間 $[0,1/n]$ なら、すべてに属する点は $0$ である。
有界性も外せない。$[n,\infty)$ は入れ子であるが、共通点 $x$ があれば $n>x$ となる $n$ で矛盾する。
:::

区間縮小法の典型的な応用として、実数は自然数で番号付けし尽くせないことを示す。

:::theorem[実数の非可算性]{#thm-anl-uncountable}
任意の実数列 $(x_n)_{n\in\N}$ に対し、すべての $n$ について $x\ne x_n$ となる $x\in[0,1]$ が存在する。
従って $\N$ から $[0,1]$ への全射も、$\N$ から $\R$ への全射も存在しない。
:::

:::proof
$I_0=[0,1]$ とし、長さ $3^{-(n-1)}$ の閉区間 $I_{n-1}=[c,c+3^{-(n-1)}]$ が定まったとする。
その左の三分の一 $[c,c+3^{-n}]$ と右の三分の一 $[c+2\cdot3^{-n},c+3^{-(n-1)}]$ は共通部分をもたないので、少なくとも一方は $x_n$ を含まない。
左が $x_n$ を含まなければ左を、そうでなければ右を $I_n$ とする。

こうして得た $I_n$ は入れ子の有界閉区間であるから、:ref[thm-anl-nested-intervals] によりすべての $I_n$ に属する $x$ がある。
$x\in I_n$ かつ $x_n\notin I_n$ なので $x\ne x_n$ である。
全射 $f:\N\to\R$ があれば $x_n=f(n)$ はすべての実数を尽くすが、これは今示したことに反する。$[0,1]$ についても同じである。
:::

一方、$\Q$ は自然数で番号付けできる（[可算集合](/math/sets-logic/cardinality/) を参照）。可算集合二つの和集合は可算であるから、無理数全体も可算でない。無理数は存在するだけでなく、有理数よりはるかに多い。

### 二分法で上限を取り戻す

逆に、順序体で区間縮小法と「どの元よりも大きい自然数がある」という条件を仮定してみよう。
後者を :term[Archimedes の性質]{en="Archimedean property" yomi="あるきめですのせいしつ"} という。この二条件から上限の公理を導くには、区間を半分に分け、上界となる位置を絞ればよい。

:::proposition[区間縮小法から上限の公理へ]{#prop-anl-bisection-sup}
順序体 $K$ が Archimedes の性質をもち、長さが任意に小さくなる入れ子の有界閉区間には必ず共通点があるとする。
このとき、$K$ の空でない上に有界な部分集合は上限をもつ。
:::

:::proof
集合を $A$ とし、$a\in A$ と上界 $u$ を選ぶ。$l_1=a-1$ は上界でなく、$u_1=u$ は上界である。
左端が上界でなく右端が上界である区間 $[l_n,u_n]$ が得られたとき、中央の点 $m=(l_n+u_n)/2$ を調べる。

- $m$ が上界なら $[l_{n+1},u_{n+1}]=[l_n,m]$ とする。
- $m$ が上界でなければ $[l_{n+1},u_{n+1}]=[m,u_n]$ とする。

これで左端は常に上界でなく、右端は常に上界である。$l_n$ は上界でないので $l_n<a'$ となる $a'\in A$ があり、$u_n$ は上界なので $a'\le u_n$、従って $l_n<u_n$ である。
中点は $l_n<m<u_n$ を満たすので（:ref[prop-anl-order-basics] の 7）、区間は入れ子である。しかも

$$
u_n-l_n=\frac{u_1-l_1}{2^{n-1}}\le\frac{u_1-l_1}{n}.
$$

最後の不等式には :ref[prop-anl-bernoulli] の後に示した $2^{n-1}\ge n$ を用いた。
$\eps>0$ に対し、Archimedes の性質により $n>(u_1-l_1)/\eps$ と選べば長さは $\eps$ 未満になる。仮定により共通点 $s$ が存在する。

もし $x\in A$ が $s<x$ を満たすなら、長さが $x-s$ 未満の段階では $u_n\le s+(u_n-l_n)<x$ となり、$u_n$ が上界であることに反する。
従って $s$ は上界である。さらに任意の $\eps>0$ に対し、長さが $\eps$ 未満の段階では $l_n\ge s-(u_n-l_n)>s-\eps$ となる。
$l_n$ は上界でないので $l_n<x$ となる $x\in A$ があり、$s-\eps<x$ である。
:ref[thm-anl-bound-epsilon] は任意の順序体で成り立つので、$s=\sup A$ である。
:::

例えば $A=(0,2/3)$ に対し、上界でない左端 $0$、上界である右端 $1$ から始めると、中点 $1/2$ は上界でなく、次の中点 $3/4$ は上界である。
続いて $5/8$ は上界でないので、三度の二分で上限は $[5/8,3/4]$ に絞られる。

```tikz caption="上界でない左端と上界である右端を保ちながら、上限を挟む区間を半分にする。"
\begin{tikzpicture}[x=0.85cm,y=0.64cm]
\draw[dashed] ({16/3},0.75) -- ({16/3},4.35) node[above] {$s=2/3$};
\foreach \row/\lo/\hi/\idx in {4/0/8/1,3/4/8/2,2/4/6/3,1/5/6/4}{
  \node[left] at (-0.15,\row) {$I_{\idx}$};
  \draw[thick,|-|] (\lo,\row) -- (\hi,\row);
}
\draw[->] (-0.15,0) -- (8.45,0);
\foreach \x/\v in {0/0,4/{1/2},5/{5/8},6/{3/4},8/1}
  \draw (\x,0.09) -- (\x,-0.09) node[below] {$\v$};
\end{tikzpicture}
```

### 数列を用いる形式

順序体 $K$ の中の数列とは、自然数 $n$ に元 $x_n\in K$ を一つずつ対応させたものである。
各 $\eps>0$（$\eps\in K$）に対してある $N\in\N$ があり、$n\ge N$ なら $\abs{x_n-x}<\eps$ となるとき、数列は $x$ に収束するという。
これを用いると、上限の公理に相当する条件を次のようにも書ける。

| 形式 | 条件 |
| --- | --- |
| 単調収束 | 上に有界で $x_n\le x_{n+1}$ を満たす数列は、$K$ の元に収束する。 |
| Bolzano–Weierstrass | 有界な数列から、添字を狭義増加に選んで収束する部分列を取り出せる。 |
| Cauchy 完備性と Archimedes の性質 | 各 $\eps>0$ に対して十分後の任意の二項の距離が $\eps$ 未満になる数列は収束し、かつ $\N$ は上に有界でない。 |

例えば単調収束は上限から直接わかる。$s=\sup\{x_n\mid n\in\N\}$ とし、$\eps>0$ に対して $s-\eps<x_N$ となる項を選ぶ。
増加性から $n\ge N$ のすべてで $s-\eps<x_N\le x_n\le s$ となるので、$x_n$ は $s$ に収束する。
上限近くで一項を見つければ、それ以降のすべての項の位置が決まる。

これらの条件から上限を導く逆向きの議論と、部分列および Cauchy 列に関する証明は [数列と級数](/math/analysis/sequences-series/) で行う。
そこでは、数列の極限の一意性や順序の保存のように、順序体の公理だけから従う性質を一般の順序体で用いる。
区間縮小法や Cauchy 完備性を使う形式では、上の二分法に必要な Archimedes の性質を明示的に仮定しておく。次節の例が示すように、順序体では Archimedes の性質は自動的には成り立たないからである。

## Archimedes の性質

### 自然数がいくらでも大きくなること

:::theorem[Archimedes の性質]{#thm-anl-archimedes}
任意の $x\in\R$ に対して $n>x$ となる $n\in\N$ が存在する。
従って $a,b>0$ なら $na>b$ となる自然数 $n$ がある。
:::

:::proof
もし $\N$ に上界があれば、:ref[ax-anl-supremum] により $\alpha=\sup\N$ が存在する。
$\alpha-1$ は上界になれないから $\alpha-1<m$ を満たす自然数 $m$ がある。
しかし $m+1$ も自然数であり、$\alpha<m+1$ となって $\alpha$ が上界であることに反する。
従って $\N$ は上に有界でない。$x=b/a$ と置いて $n>x$ の両辺に $a>0$ を掛ければ後半を得る。
:::

この証明は上限の公理だけを用いるので、上限の公理を満たす任意の順序体は Archimedes の性質をもつ。

:::corollary[任意に小さい単位分数]{#cor-anl-unit-fractions}
各 $\eps>0$ に対して $1/N<\eps$ を満たす $N\in\N$ がある。
また、すべての自然数 $n$ について $0\le d\le1/n$ ならば $d=0$ である。
:::

:::proof
:ref[thm-anl-archimedes] を $x=1/\eps$ に用いれば $N>1/\eps$ となり、:ref[prop-anl-order-basics] の 6 により $1/N<\eps$ である。
もし $d>0$ ならば $\eps=d$ に対して $1/N<d$ を得るので、仮定に反する。
:::

以上で、本章の前半の結果を次のようにまとめることができる。

:::corollary[実数の連続性の同値な形式]{#cor-anl-completeness-forms}
順序体 $K$ について、次の三条件は同値である。

1. 上限の公理：空でない上に有界な部分集合は上限をもつ。
2. Dedekind の性質：どの切断についても、下組に最大元があるか上組に最小元がある。
3. $K$ は Archimedes の性質をもち、長さが任意に小さくなる入れ子の有界閉区間は共通点をもつ。
:::

:::proof
1 と 2 の同値は :ref[thm-anl-dedekind] である。
1 を仮定すると、:ref[thm-anl-archimedes] と :ref[thm-anl-nested-intervals] の証明は上限の公理しか使わないので、$K$ でも成り立ち 3 が従う。
3 から 1 は :ref[prop-anl-bisection-sup] である。
:::

数列による形式（単調収束、Bolzano–Weierstrass、Cauchy 完備性と Archimedes の性質）もこれらと同値であることは、[数列と級数](/math/analysis/sequences-series/) で証明する。

### Archimedes 的でない順序体

Archimedes の性質は順序体の公理からは従わない。それを示す例を構成する。

:::example[有理関数体]{#exm-anl-non-archimedean}
実係数の多項式 $p,q$（$q\ne0$）の商 $p/q$ 全体を $K=\R(t)$ とする。通常の分数の計算によってこれが体になることは代数学で扱う（[体と Galois 理論](/math/fields-galois/)）。
零でない多項式 $p$ の最高次の係数を $\operatorname{lc}(p)$ と書き、$f=p/q\ne0$ が **正** であることを $\operatorname{lc}(p)\operatorname{lc}(q)>0$ で定める。

*定義が表示によらないこと。* $p/q=p'/q'$ なら $pq'=p'q$ であり、最高次係数は積に関して乗法的なので $\operatorname{lc}(p)\operatorname{lc}(q')=\operatorname{lc}(p')\operatorname{lc}(q)$ である。
両辺に $\operatorname{lc}(q)\operatorname{lc}(q')$ を掛けると $\operatorname{lc}(p)\operatorname{lc}(q)\cdot\operatorname{lc}(q')^2=\operatorname{lc}(p')\operatorname{lc}(q')\cdot\operatorname{lc}(q)^2$ となり、二つの積の符号は一致する。

*正の元全体 $P$ の性質。* $f\ne0$ なら $f$ と $-f$ のちょうど一方が $P$ に属する。$f,g\in P$ なら $fg\in P$ である（最高次係数は乗法的）。
和については、分母と分子に $-1$ を掛けて $f=p/q$、$g=r/w$ の分母の最高次係数を正にしておくと、$p,r$ の最高次係数も正である。
$f+g=(pw+rq)/(qw)$ の分子は最高次係数が正の多項式二つの和であり、次数が異なれば高い方の、等しければ二つの正の数の和が最高次係数になるので、いずれも正である。

*順序。* $f\le g$ を「$g-f\in P\cup\{0\}$」で定める。(O1) は明らか、(O2) は $g-f$ と $f-g$ が同時に正にならないこと、(O3) は $P$ が和で閉じること、(O4) は $g-f$ か $f-g$ が $P\cup\{0\}$ に属すること、(O5) は $(g+c)-(f+c)=g-f$、(O6) は $P$ が積で閉じることから従う。こうして $K$ は順序体になる。

*Archimedes の性質が成り立たないこと。* 実数の定数 $1,2,3,\dots$ 全体は $K$ の帰納的集合なので、$K$ の $\N$ はこれに含まれる。
各自然数 $n$ について $t-n$ の最高次係数は $1>0$ であるから $t>n$ である。従って $\N$ は $t$ を上界にもち、:ref[thm-anl-archimedes] により $K$ は上限の公理を満たさない。
同様に $1/n-1/t=(t-n)/(nt)$ は正なので $0<1/t<1/n$ がすべての $n$ で成り立ち、:ref[cor-anl-unit-fractions] の後半も $K$ では成り立たない。
:::

### 整数で実数を挟む

実数 $x$ を超える最初の整数がわかれば、間隔 $1$ の目盛りのどこに $x$ があるかが決まる。
$x$ が負の場合も含めて、この整数が存在することを確かめよう。

:::proposition[整数部分]{#prop-anl-floor}
各 $x\in\R$ に対し、$k\le x<k+1$ となる整数 $k$ がただ一つ存在する。
この $k$ を :term[床関数]{en="floor function" yomi="ゆかかんすう"} の $x$ での値といい、$\lfloor x\rfloor$ と書く。$\lfloor x\rfloor$ は $x$ 以下の整数のうち最大のものである。
:::

:::proof
$S=\{j\in\Z\mid j>x\}$ とおく。:ref[thm-anl-archimedes] により $n>x$ となる $n\in\N$ があるので $S\ne\varnothing$ である。
また $M>-x$ となる $M\in\N$ をとれば、$S$ の元 $j$ は $j>x>-M$ を満たすので、$S$ は整数 $-M$ で下から抑えられる。
:ref[prop-anl-integers] の 3 により最小元 $j$ が存在する。$j-1\notin S$ なので $j-1\le x<j$ であり、$k=j-1$ が求める整数である。

別の整数 $k'$ も同じ条件を満たし $k<k'$ とすると、:ref[prop-anl-integers] の 2 により $k'\ge k+1>x$ となり矛盾する。$k'<k$ の場合も同様である。
最後に、整数 $j\le x$ が $j>k$ を満たせば $j\ge k+1>x$ となるので、$k$ は $x$ 以下の最大の整数である。
:::

例えば $\lfloor-7/3\rfloor=-3$ である。負の数を扱うときも、目盛りの左側の整数を選ぶ。
$\eps>0$ に対して $N=\lfloor1/\eps\rfloor+1$ とすれば、$\eps$ の逆数が整数である場合を含め、必ず $N>1/\eps$、すなわち $1/N<\eps$ となる。

## 有理数の稠密性

### 開区間の中の有理数

:::theorem[有理数の稠密性]{#thm-anl-rational-density}
任意の実数 $a<b$ に対して $a<r<b$ を満たす $r\in\Q$ が存在する。
この性質を、$\Q$ は $\R$ において :term[稠密]{en="dense" yomi="ちゅうみつ"} であるという。
:::

:::proof
:ref[cor-anl-unit-fractions] により $1/n<b-a$ となる自然数 $n$ を選ぶ。間隔 $1/n$ の目盛りで $a$ のすぐ右にある点をとるために、$m=\lfloor na\rfloor+1$ とする。
$na<m\le na+1$ であるから

$$
a<\frac mn\le a+\frac1n<b.
$$

よって $r=m/n$ が求める有理数である。
:::

区間 $(a,b)$ に有理数が一つあれば、$(a,r)$ と $(r,b)$ にも有理数があるので、どの開区間にも有理数は無限に多く含まれる。

:::example[端点が負である区間]
$-7/5<r<-4/3$ を満たす有理数を分母 $16$ の形で探す。
区間の幅は $1/15>1/16$ なので、上の証明で $n=16$ とできる。
$16(-7/5)=-112/5$ の床は $-23$ であり、一つ右の整数は $-22$ となる。
従って $r=-22/16=-11/8$ である。
:::

:::proposition[実数は有理数の上限である]{#prop-anl-sup-rationals}
各 $x\in\R$ について $x=\sup\{r\in\Q\mid r<x\}$ である。
:::

:::proof
集合 $Q_x=\{r\in\Q\mid r<x\}$ は :ref[thm-anl-rational-density] により $(x-1,x)$ の有理数を含むので空でなく、$x$ を上界にもつ。
任意の $\eps>0$ に対し、同じ定理で $x-\eps<r<x$ となる有理数 $r\in Q_x$ を選べる。:ref[thm-anl-bound-epsilon] により $x=\sup Q_x$ である。
:::

このように各実数は有理数から上限の操作で復元される。実数の構成では、逆にこの事実を出発点として有理数の集合から実数を定義する。

### n 乗根の存在と無理数

有理数はどの開区間にも入るが、実数の値を厳密に指定したとき、その値が有理数になるとは限らない。
有理数でない実数を :term[無理数]{en="irrational number" yomi="むりすう"} という。
まず、正の数の $n$ 乗根が実数の中に存在することを上限の公理から示す。

:::theorem[n 乗根の存在]{#thm-anl-nth-root}
$n\in\N$ と $t\ge0$ に対し、$s^n=t$ を満たす $s\ge0$ がただ一つ存在する。
この $s$ を $t$ の :term[n 乗根]{en="n-th root" yomi="えぬじょうこん"} といい、$\sqrt[n]{t}$ または $t^{1/n}$ と書く。$n=2$ のときは $\sqrt t$ と書く。
:::

:::proof
**段階 1（単調性と一意性）。** $0\le x<y$ ならば $x^n<y^n$ である。$n=1$ では明らかであり、$0\le x^n<y^n$ を仮定すると

$$
x^{n+1}=x^nx\le x^ny<y^ny=y^{n+1}
$$

となる（最初の不等号は $x\le y$ に $x^n\ge0$ を掛け、次の不等号は $x^n<y^n$ に $y>0$ を掛けた）。
従って非負の解が二つあれば、それらは等しい。

**段階 2（差の評価）。** $0\le x\le y\le c$ ならば、展開して隣り合う項を打ち消すと

$$
y^n-x^n=(y-x)\sum_{k=0}^{n-1}y^{n-1-k}x^k\le nc^{n-1}(y-x)
\label{eq-anl-power-difference}
$$

である。和の各項は $c^{n-1}$ 以下だからである。

**段階 3（存在）。** $t=0$ なら $s=0$ である。$t>0$ とし、$A=\{x\ge0\mid x^n<t\}$ とおく。
$0\in A$ である。$x\ge1+t$ なら $x\ge1$ より $x^{n-1}\ge1$（段階 1）なので $x^n\ge x\ge1+t>t$ であり、$1+t$ は $A$ の上界である。
:ref[ax-anl-supremum] により $s=\sup A\ge0$ が存在する。$c=s+1$ とおく。

$s^n<t$ と仮定し、$h=\min\{1,(t-s^n)/(2nc^{n-1})\}>0$ とおく。$0\le s<s+h\le c$ なので :eqref[eq-anl-power-difference] より

$$
(s+h)^n-s^n\le nc^{n-1}h\le\frac{t-s^n}2
$$

となり、$(s+h)^n<t$、すなわち $s+h\in A$ である。これは $s$ が上界であることに反する。

$s^n>t$ と仮定すると $s>0$ である。$h=\min\{s,(s^n-t)/(2nc^{n-1})\}>0$ とおくと $0\le s-h<s\le c$ であり

$$
s^n-(s-h)^n\le nc^{n-1}h\le\frac{s^n-t}2
$$

から $(s-h)^n\ge(s^n+t)/2>t$ を得る。$x\in A$ が $x\ge s-h$ を満たせば段階 1 により $x^n\ge(s-h)^n>t$ となるので、$A$ の元はすべて $s-h$ より小さい。
すると $s-h$ は $s$ より小さい上界になり、矛盾する。:ref[prop-anl-order-basics] の 1 により $s^n=t$ である。
:::

一意性から、$a,b\ge0$ に対し $\sqrt[n]{ab}=\sqrt[n]a\,\sqrt[n]b$ が成り立つ。右辺も非負で、$n$ 乗すると $ab$ になるからである。

:::theorem[√2 の存在と無理性]{#thm-anl-sqrt-two}
$s^2=2$ となる正の実数 $s$ がただ一つ存在する。この数 $\sqrt{2}$ は無理数である。
:::

:::proof
存在と一意性は :ref[thm-anl-nth-root] の $n=2$、$t=2$ の場合である。$1^2<2<2^2$ と段階 1 の単調性から $1<\sqrt2<2$ である。

$\sqrt2\in\Q$ と仮定する。$\sqrt2=m/n$ なら $n\sqrt2=m\in\Z$ なので、$S=\{q\in\N\mid q\sqrt2\in\Z\}$ は空でない。
:ref[thm-anl-well-ordering] により最小元 $q$ があり、$p=q\sqrt2\in\Z$ とおく。
$q'=p-q=q(\sqrt2-1)$ は整数であり、$0<\sqrt2-1<1$ から $0<q'<q$ を満たすので $q'\in\N$ である。
ところが $q'\sqrt2=q\cdot2-q\sqrt2=2q-p\in\Z$ であるから $q'\in S$ となり、$q$ の最小性に反する。
:::

この証明は約分や素因数分解を使わず、整列性だけを用いている。同じ方法で、平方数でない自然数の平方根が無理数であることも示せる（演習問題）。

:::corollary[無理数の稠密性]{#cor-anl-irrational-density}
どの開区間 $(a,b)$、$a<b$ にも無理数がある。
:::

:::proof
:ref[thm-anl-rational-density] を区間 $(a-\sqrt{2},b-\sqrt{2})$ に用い、有理数 $r$ を一つとる。
すると $a<r+\sqrt{2}<b$ である。もし $r+\sqrt{2}$ が有理数なら、その数から $r$ を引いた $\sqrt{2}$ も有理数となり、:ref[thm-anl-sqrt-two] に反する。
:::

:::example[有理数体で上限が存在しない集合]
$B=\{q\in\Q\mid 0\le q,\ q^2<2\}$ を考える。$\R$ では $\sup B=\sqrt{2}$ である。
実際、$q\in B$ なら段階 1 の単調性から $q<\sqrt2$ であり、$\eps>0$ に対し $\max\{0,\sqrt{2}-\eps\}<q<\sqrt{2}$ を満たす有理数は $B$ に属する。

一方、$\Q$ の中で $B$ の上界となる数 $u$ は $\sqrt{2}$ より大きい。等号は無理性により生じず、それより小さい数は今示した近似により上界になれないからである。
そこで $\sqrt{2}<v<u$ となる有理数 $v$ を選ぶと、$v$ も $B$ の上界である。
どの有理数の上界も最小ではないので、$B$ は $\Q$ の中に上限をもたない。従って $\Q$ は上限の公理を満たさない。
:::

### 10 進展開

日常的に用いる小数表示も、床関数と上限から厳密に基礎付けられる。

:::proposition[10 進展開]{#prop-anl-decimal}
$x\in\R$ と $n\in\Z_{\ge0}$ に対し $k_n=\lfloor10^nx\rfloor$、$d_n=k_n/10^n$ とおく。

1. $d_n\le x<d_n+10^{-n}$、$d_0\le d_1\le d_2\le\cdots$ であり、$x=\sup_n d_n$ である。
2. $n\ge1$ について $a_n=k_n-10k_{n-1}$ は $\{0,1,\dots,9\}$ に属し、$d_n=k_0+\sum_{j=1}^n a_j10^{-j}$ である。
3. ある $m$ 以降のすべての $j>m$ で $a_j=9$ となることはない。
:::

:::proof
1. 床関数の定義から $k_n\le10^nx<k_n+1$ であり、$10^n$ で割れば最初の不等式を得る。
$10k_{n-1}\le10^nx$ で $10k_{n-1}$ は整数だから、:ref[prop-anl-floor] により $10k_{n-1}\le k_n$、すなわち $d_{n-1}\le d_n$ である。
$x$ は $\{d_n\}$ の上界である。$\eps>0$ に対し、:ref[prop-anl-bernoulli] より $10^n\ge1+9n>n$ なので、$n>1/\eps$ と選べば $10^{-n}<\eps$ となり、$x-\eps<x-10^{-n}<d_n$ である。:ref[thm-anl-bound-epsilon] により $x=\sup_nd_n$ である。
2. $k_{n-1}\le10^{n-1}x<k_{n-1}+1$ に $10$ を掛けると $10^nx<10k_{n-1}+10$ であるから、整数 $k_n$ は $k_n<10k_{n-1}+10$、:ref[prop-anl-integers] の 2 により $k_n\le10k_{n-1}+9$ を満たす。1 と合わせて $0\le a_n\le9$ である。
また $d_j-d_{j-1}=(k_j-10k_{j-1})/10^j=a_j10^{-j}$ を $j=1,\dots,n$ について加えればよい。
3. $j>m$ で $a_j=9$ と仮定する。$n>m$ について

$$
d_n=d_m+9\sum_{j=m+1}^n10^{-j}=d_m+10^{-m}-10^{-n}
$$

である（等比数列の有限和）。1 により $\delta=d_m+10^{-m}-x>0$ である。1 の証明と同様に $10^{-n}<\delta$ となる $n>m$ をとると
$d_n=d_m+10^{-m}-10^{-n}>x$ となり、1 の $d_n\le x$ に反する。
:::

$k_0$ と数字の列 $a_1a_2a_3\cdots$ を並べて $x=k_0.a_1a_2a_3\cdots$ と書く。例えば $1/3=0.333\cdots$ である。
級数の言葉を用いれば $x=k_0+\sum_{j=1}^\infty a_j10^{-j}$ であり、逆に $9$ が続かない任意の数字列はこの形で一つの実数を定める（[数列と級数](/math/analysis/sequences-series/) の単調収束定理による）。
$0.999\cdots$ のように $9$ が続く列の和は $1$ であるが、上の手続きは $1$ に対して $1.000\cdots$ を与え、$9$ が続く表示を生まない。

## 演習問題

::::exercise{level=basic}
$A=\{2-3/n\mid n\in\N\}$ の上限・下限を求め、それぞれが $A$ に属するかを判定せよ。
:::solution
$n\ge1$ より $3/n\le3$ なので $-1\le2-3/n<2$ である。左の等号は $n=1$ で成り立つから、下限は最小値 $-1$ であり $A$ に属する。

上限の候補 $2$ は上界である。$\eps>0$ が与えられたとき $n=\lfloor3/\eps\rfloor+1$ ととれば $n>3/\eps$、すなわち $3/n<\eps$ である。
従って $2-\eps<2-3/n$ であり、:ref[thm-anl-bound-epsilon] より上限は $2$ である。
すべての $n$ で $3/n>0$ なので、$2$ は $A$ に属さない。
:::
::::

::::exercise{level=basic}
空でない有界集合 $A\subset\R$ と $c<0$ に対して $cA=\{ca\mid a\in A\}$ とおく。$\sup(cA)$ を $\inf A$ で表せ。
:::solution
$t=\inf A$ とおく。$t\le a$ に負の数 $c$ を掛けると $ca\le ct$ になるので、$ct$ は $cA$ の上界である。
任意の $\eps>0$ に対し、$t\le a<t+\eps/\abs{c}$ となる $a\in A$ を選ぶ。$c<0$ を掛けると不等号が逆転し

$$
ct-\eps<ca\le ct
$$

となる。:ref[thm-anl-bound-epsilon] により $\sup(cA)=c\inf A$ である。負の係数を掛けると上側と下側の役割が入れ替わる。
:::
::::

::::exercise{level=basic}
$1.41<q<1.42$ を満たす有理数を、分母を $128$ とする分数の形で一つ求めよ。
:::solution
幅は $1/100>1/128$ であるから、:ref[thm-anl-rational-density] の証明で $n=128$ とできる。左端に $128$ を掛けると $180.48$、右端に掛けると $181.76$ となる。
$m=\lfloor180.48\rfloor+1=181$ は $180.48<181<181.76$ を満たすので、$q=181/128=1.4140625$ が条件を満たす。
:::
::::

::::exercise{level=basic}
$a,b\ge0$ に対して $\sqrt{ab}\le(a+b)/2$ を示し、等号が成り立つのは $a=b$ のときに限ることを示せ。
:::solution
:ref[thm-anl-nth-root] により $x=\sqrt a$、$y=\sqrt b$ が存在し、その後の注意により $\sqrt{ab}=xy$ である。

$$
\frac{a+b}2-\sqrt{ab}=\frac{x^2-2xy+y^2}2=\frac{(x-y)^2}2\ge0
$$

であるから不等式が成り立つ。等号は $(x-y)^2=0$、すなわち :ref[prop-anl-order-basics] の 5 により $x=y$ のときに限る。
$x,y\ge0$ と平方根の一意性から、$x=y$ は $a=x^2=y^2=b$ と同値である。
:::
::::

::::exercise{level=basic}
$L=\{q\in\Q\mid q\le0\text{ または }q^2<2\}$、$U=\Q\setminus L$ とおく。$(L,U)$ が $\Q$ の切断であり、$L$ に最大元がなく $U$ に最小元がないことを示せ。
:::solution
$U=\{q\in\Q\mid q>0,\ q^2\ge2\}$ であり、:ref[thm-anl-sqrt-two] により $q^2\ne2$ なので $U=\{q\in\Q\mid q>0,\ q^2>2\}$ である。
$0\in L$、$2\in U$ なので両方とも空でなく、定義から $L\cup U=\Q$、$L\cap U=\varnothing$ である。
$x\in L$、$y\in U$ とする。$x\le0$ なら $x<y$ である。$x>0$ なら $x^2<2<y^2$ であり、:ref[thm-anl-nth-root] の段階 1 から $x\ge y$ はありえないので $x<y$ である。従って切断である。

$x\in L$ が最大元でないことを示す。$x\le0$ なら $1\in L$ が $x$ より大きい。$x>0$ なら $x^2<2$ より $x<\sqrt2$ であり、:ref[thm-anl-rational-density] で $x<r<\sqrt2$ となる有理数 $r$ をとると $r>0$、$r^2<2$ なので $r\in L$ である。
同様に $y\in U$ なら $y>\sqrt2$ であり、$\sqrt2<r<y$ となる有理数 $r$ は $U$ に属するので、$y$ は最小元でない。
よって $\Q$ は Dedekind の性質をもたず、:ref[thm-anl-dedekind] により上限の公理も満たさない。
:::
::::

::::exercise{level=advanced}
空でない有界集合 $A$ の直径を $D=\sup\{\abs{x-y}\mid x,y\in A\}$ とする。$D=\sup A-\inf A$ を示せ。
:::hint
上限に近い元と下限に近い元を別々に選ぶ。
:::
:::solution
$s=\sup A,\ t=\inf A$ とおく。$t\le x,y\le s$ より $-(s-t)\le x-y\le s-t$、すなわち $\abs{x-y}\le s-t$ であるから、直径は存在し $D\le s-t$ である。
他方、$\eps>0$ に対し $x>s-\eps/2$ と $y<t+\eps/2$ を満たす元をそれぞれ選べば

$$
s-t-\eps<x-y\le\abs{x-y}\le D
$$

となる。もし $D<s-t$ ならば $\eps=s-t-D$ として $D<D$ となり矛盾する。従って $D=s-t$ である。
:::
::::

::::exercise{level=advanced}
任意の $x\in\R$ に対して $q_n=\lfloor 2^n x\rfloor/2^n$ とおく。
区間 $[q_n,q_n+2^{-n}]$ は入れ子であり、共通点は $x$ だけであることを示せ。
:::solution
$k_n=\lfloor2^n x\rfloor$ とおけば $k_n\le2^n x<k_n+1$ であるから、$q_n\le x<q_n+2^{-n}$ であり、$x$ はすべての区間に属する。
$k_n\le2^nx<k_n+1$ に $2$ を掛けると $2k_n\le2^{n+1}x<2k_n+2$ を得る。
:ref[prop-anl-floor] により $k_{n+1}$ は $2^{n+1}x$ 以下の最大の整数なので、$k_{n+1}$ は $2k_n$ または $2k_n+1$ である。
従って次の区間 $[k_{n+1}/2^{n+1},(k_{n+1}+1)/2^{n+1}]$ は前の区間 $[2k_n/2^{n+1},(2k_n+2)/2^{n+1}]$ の左半分か右半分であり、入れ子である。

共通点 $y$ があれば、$x,y$ はともに長さ $2^{-n}$ の区間に属するので、各 $n$ で $\abs{x-y}\le2^{-n}\le1/n$ となる（$2^n\ge n$）。
よって :ref[cor-anl-unit-fractions] により $\abs{x-y}=0$、すなわち $y=x$ である。
:::
::::

::::exercise{level=advanced}
$A,B\subset(0,\infty)$ を空でない上に有界な集合とし、$AB=\{ab\mid a\in A,\ b\in B\}$ とおく。$\sup(AB)=\sup A\cdot\sup B$ を示せ。また正の数という仮定を外すと成り立たない例を挙げよ。
:::solution
$s=\sup A$、$t=\sup B$ とおくと $s,t>0$ である。$0<a\le s$、$0<b\le t$ から $ab\le sb\le st$ なので、$st$ は上界である。
$\eps>0$ とし、$\eta=\min\{s,t,\eps/(s+t)\}>0$ とおく。$a>s-\eta\ge0$、$b>t-\eta\ge0$ となる $a\in A$、$b\in B$ を選ぶと

$$
ab>(s-\eta)b\ge(s-\eta)(t-\eta)=st-\eta(s+t)+\eta^2>st-\eta(s+t)\ge st-\eps
$$

である。最初の不等号は $b>0$、次は $s-\eta\ge0$ による。:ref[thm-anl-bound-epsilon] により $\sup(AB)=st$ である。

正という仮定がないと、$A=B=\{-1,0\}$ では $\sup A=\sup B=0$ であるが、$AB=\{0,1\}$ の上限は $1$ である。
:::
::::

::::exercise{level=advanced}
$m\in\N$ が平方数でない（$m=k^2$ となる $k\in\N$ がない）ならば $\sqrt m$ は無理数であることを示せ。
:::hint
$k=\lfloor\sqrt m\rfloor$ を用い、:ref[thm-anl-sqrt-two] の証明を真似る。
:::
:::solution
$k=\lfloor\sqrt m\rfloor$ とおくと $k\le\sqrt m<k+1$ である。$\sqrt m=k$ なら $m=k^2$ となるので $k<\sqrt m<k+1$、すなわち $0<\sqrt m-k<1$ である。
$\sqrt m\in\Q$ と仮定すると、$S=\{q\in\N\mid q\sqrt m\in\Z\}$ は空でなく、:ref[thm-anl-well-ordering] により最小元 $q$ をもつ。$p=q\sqrt m\in\Z$ とおく。
$q'=q(\sqrt m-k)=p-kq$ は整数で $0<q'<q$ なので $q'\in\N$ である。さらに

$$
q'\sqrt m=q m-k\,q\sqrt m=qm-kp\in\Z
$$

であるから $q'\in S$ となり、$q$ の最小性に反する。
:::
::::

::::exercise{level=advanced}
:ref[exm-anl-non-archimedean] の順序体 $K=\R(t)$ について次を示せ。(a) $\N$ は $K$ の中に上限をもたない。(b) 区間 $[0,1/n]$（$n\in\N$）の共通部分は $0$ と $1/t$ を含む。(c) (b) は :ref[thm-anl-nested-intervals] の一意性の部分と矛盾しない。
:::solution
(a) $s\in K$ が $\N$ の上界なら、すべての $n$ で $n+1\le s$、すなわち $n\le s-1$ なので $s-1$ も上界である。$s-1<s$ だから、どの上界も最小ではない。

(b) $0\in[0,1/n]$ は明らかである。例で示したように $1/n-1/t=(t-n)/(nt)$ の分子と分母の最高次係数はそれぞれ $1$ と $n$ で正なので $1/t<1/n$、また $1/t>0$ である。

(c) 一意性の部分は「各 $\eps>0$ について長さが $\eps$ 未満になる区間がある」ことを仮定する。$K$ で $\eps=1/t$ をとると、すべての $n$ で $1/n>1/t$ なのでこの仮定は満たされない。
実際、二つの共通点の差 $1/t$ がちょうどこの反例を与えている。
:::
::::

## 確認テスト

::::quiz
空でない上に有界な $A\subset\R$ について、必ず成り立つものをすべて選べ。

- [x] すべての $a\in A$ は $\sup A$ 以下である。
- [ ] $\sup A$ は $A$ の元である。
- [x] $u<\sup A$ ならば $u<a$ となる $a\in A$ がある。

:::explanation
(A) は上限が上界であることそのものである。(C) は、上限が最小の上界なので、それより小さい $u$ は上界になれないことによる（:ref[thm-anl-bound-epsilon] で $\eps=\sup A-u$ とした場合）。
(B) は誤りで、例えば $(0,1)$ の上限は $1$ だが、$1$ はこの集合に属さない。
:::
::::

::::quiz{answer="-4"}
$\lfloor-17/5\rfloor$ を求めよ。

:::explanation
$-17/5=-3.4$ であり、$-4\le-3.4<-3$ なので、:ref[prop-anl-floor] の条件 $k\le x<k+1$ を満たす整数は $k=-4$ である。負の数では小数部分を切り捨てるのではなく、左側の整数を選ぶことに注意する。
:::
::::

::::quiz{answer="101"}
$1/N<1/100$ を満たす最小の自然数 $N$ を求めよ。

:::explanation
正の数の逆数をとると不等号が逆になる（:ref[prop-anl-order-basics] の 6）。従って条件は $N>100$ と同値であり、最小値は $101$ である。$N=100$ では等号になり条件を満たさない。
:::
::::

::::quiz
空でない有界閉区間が入れ子になっているとき、共通点の一意性まで保証する追加条件はどれか。

- [ ] すべての区間が同じ長さをもつ。
- [x] 任意の $\eps>0$ について長さが $\eps$ 未満になる段階がある。
- [ ] すべての端点が有理数である。

:::explanation
(B) が :ref[thm-anl-nested-intervals] の後半の条件である。二つの異なる共通点の距離は正であり、その距離より短い区間には両方の点は入らない。
同じ区間 $[0,1]$ を繰り返す列は、端点が有理数で長さも一定であるが、$[0,1]$ のすべての点が共通点となる。従って (A) と (C) では一意性は従わない。
:::
::::

::::quiz
次の順序体のうち、Archimedes の性質をもつものをすべて選べ。

- [x] 有理数体 $\Q$
- [x] 実数体 $\R$
- [ ] 有理関数体 $\R(t)$（:ref[exm-anl-non-archimedean] の順序）

:::explanation
$\R$ は :ref[thm-anl-archimedes] による。$\Q$ では、有理数 $x$ を実数とみて同じ定理を適用すれば $n>x$ となる自然数が得られる（$\Q$ の自然数と $\R$ の自然数は同じ集合である）。
$\R(t)$ では、すべての自然数 $n$ について $t>n$ となるので Archimedes の性質は成り立たない。
なお $\Q$ は Archimedes の性質をもつが上限の公理は満たさない。従って :ref[cor-anl-completeness-forms] の条件 3 で欠けているのは区間縮小法の方である。
:::
::::

---
title: 最小多項式と Jordan 標準形
description: 多項式の整除と最小多項式による空間の直和分解、冪零作用素の鎖の構成から Jordan 標準形の存在と一意性を示し、相似の判定・行列の冪・線形漸化式・実 Jordan 標準形に応用する。
notation:
  - symbol: 'm_T(t)'
    meaning: 線形変換 T の最小多項式
  - symbol: 'm_{T,v}(t)'
    meaning: ベクトル v を消す T の最小多項式
  - symbol: '\gcd(p,q)'
    meaning: 多項式 p, q の最大公約元
  - symbol: 'C(f)'
    meaning: モニック多項式 f の同伴行列
  - symbol: '\widetilde E_\lambda(T)'
    meaning: 固有値 λ に対する T の広義固有空間
  - symbol: 'K_j'
    meaning: 冪零作用素の j 乗の核
  - symbol: 'J_s(\lambda)'
    meaning: 大きさ s、固有値 λ の Jordan ブロック
---

## この章の目標

- 最小多項式を定義し、特性多項式との関係（同じ既約因子をもつこと）と対角化可能性の判定を証明する。
- 多項式の最大公約元と Bézout の等式を用いて、空間を不変部分空間の直和に分解する（一次分解定理）。
- 冪零作用素の鎖による基底を構成し、核の次元列から Jordan ブロックの大きさを復元する。
- Jordan 標準形の存在と一意性を証明し、相似の判定、実 Jordan 標準形、Jordan–Chevalley 分解を導く。
- 標準形を用いて行列の冪、線形漸化式、可換な行列などを計算する。

## 最小多項式

本章では、特に断らない限り $K$ は任意の体、$V$ は $K$ 上の $n$ 次元ベクトル空間（$n\ge1$）、$T\in\End(V)$ とする。
多項式 $p\in K[t]$ に $T$ を代入した $p(T)$ は :ref[lem-linalg-polynomial-substitution] の意味であり、$(pq)(T)=p(T)q(T)$ などの規則を断りなく用いる。
$q\mid p$ は、ある $h\in K[t]$ により $p=hq$ と書けること（$q$ が $p$ を割り切ること）を表す。

[固有値と対角化](/math/linear-algebra/eigenvalues/) の章で見たように、固有ベクトルを集めても基底を作れない作用素がある。
それでも例えば $(T-I)^2=0$ が分かっていれば、$T$ のすべての冪は $I,T$ の一次結合に帰着する。
任意の体で成り立つ Cayley–Hamilton の定理（:ref[thm-linalg-cayley-hamilton]）は $\Phi_T(T)=0$ を保証するが、もっと次数の低い多項式で $T$ が消えることも多い。
そこで $T$ が満たす多項式関係のうち、次数が最も低いものを調べる。
数値を指定した例・演習・確認テストでは、別に体を指定しない限り係数を実数とする。

:::definition[消去多項式と最小多項式]{#def-linalg-minimal-polynomial}
$p\in K[t]$ が $p(T)=0$ を満たすとき、$p$ を $T$ の :term[消去多項式]{en="annihilating polynomial" yomi="しょうきょたこうしき"} という。
最高次係数が $1$ の多項式を :term[モニック]{en="monic" yomi="もにっく"} な多項式という。
零でない消去多項式のうち次数が最小のモニックなものを $T$ の :term[最小多項式]{en="minimal polynomial" yomi="さいしょうたこうしき"} といい、$m_T$ と書く。
正方行列 $A\in M_n(K)$ についても同様に $m_A$ を定める。
:::

:::proposition[最小多項式の基本性質]{#prop-linalg-minimal-polynomial-divides}
1. $m_T$ は存在して一意である。
2. $p\in K[t]$ について、$p(T)=0$ であることと $m_T\mid p$ であることは同値である。
3. $m_T$ は特性多項式 $\Phi_T$ を割り切る。特に $1\le\deg m_T\le n$ である。
4. 任意の基底 $\mathcal B$ について $m_T=m_{[T]_{\mathcal B}}$ である。特に相似な行列は同じ最小多項式をもつ。
:::

:::proof
:ref[thm-linalg-cayley-hamilton] により $\Phi_T$ は零でない消去多項式である。零でない消去多項式のうち次数最小のものをモニックにして $m$ とする。:ref[lem-linalg-polynomial-division] により $p=qm+r$ と割ると $p(T)=r(T)$ である。消去多項式 $p$ に対して $r\ne0$ なら、余りをモニックにしたものが最小性に反する。よって消去多項式はすべて $m$ の倍数である。逆に倍数は :ref[lem-linalg-polynomial-substitution] により $T$ を消す。等次数のモニックな候補は互いに割り切るから一致する。

特に $m\mid\Phi_T$ である。非零定数は $V\ne0$ 上の非零のスカラー倍を表すから $1\le\deg m\le n$ となる。表現行列を取る操作は和と積と恒等写像を保ち、零写像だけを零行列に移す。したがって任意の基底で消去多項式の集合がそのまま保たれ、相似でも $m$ は変わらない。
:::

2 により、$m_A$ を求めるには $I,A,A^2,\ldots$ を順に並べ、初めて一次従属になる $A^k$ を見つければよい。
$A^k=-\sum_{j<k}c_jA^j$ ならば $m_A=t^k+\sum_{j<k}c_jt^j$ である。
実際、次数 $k$ 未満の消去多項式は $I,\ldots,A^{k-1}$ の非自明な一次関係を与えてしまう。

:::example[重根の有無]
二つの行列を比べる。

$$
A=\diag(2,2),\qquad B=\begin{pmatrix}2&1\\0&2\end{pmatrix}
$$

特性多項式はともに $(t-2)^2$ である。
$A-2I=O$ だが $B-2I\ne O,\ (B-2I)^2=O$ なので、最小多項式は順に $t-2,\ (t-2)^2$ となる。
:::

:::example[回転の最小多項式]
実平面の $90^\circ$ 回転を考える。

$$
T=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\qquad T^2=-I.
$$

$T$ はスカラー行列ではないので $I,T$ は一次独立であり、一次式では消えない。
よって最小多項式は $t^2+1$ である。
これは実数上では既約だが、複素数上では $(t-i)(t+i)$ に分かれる。
:::

### 多項式の最大公約元と既約分解

最小多項式を使って空間を分解するには、整数の場合と同様の多項式の整除の理論が必要になる。
零でない定数は $K[t]$ で乗法逆元をもつので、多項式の約元は定数倍を除いてしか決まらない。
そこでモニックな多項式に正規化して扱う。

まず、モニックな $d,d'$ が $d\mid d'$ かつ $d'\mid d$ を満たせば $d=d'$ である。
実際 $d'=hd,\ d=h'd'$ と書くと $d=h'hd$ であり、:ref[lem-linalg-polynomial-division] の 1 で $d$ を約分して $h'h=1$ を得る。
次数を比べると $h,h'$ は零でない定数であり、最高次係数を比べて $h=1$ である。

:::definition[最大公約元と互いに素]{#def-linalg-polynomial-gcd}
$p,q\in K[t]$ の少なくとも一方は零でないとする。
モニック多項式 $d$ が $p,q$ の公約元（$d\mid p$ かつ $d\mid q$）であり、$p,q$ の任意の公約元が $d$ を割り切るとき、$d$ を $p,q$ の :term[最大公約元]{en="greatest common divisor" yomi="さいだいこうやくげん"} といい、$\gcd(p,q)$ と書く。
$\gcd(p,q)=1$ のとき $p,q$ は :term[互いに素]{en="coprime" yomi="たがいにそ"} であるという。
:::

:::lemma[最大公約元の存在と Bézout の等式]{#lem-linalg-bezout}
$p,q\in K[t]$ の少なくとも一方は零でないとする。

1. $\gcd(p,q)$ はただ一つ存在し、ある $a,b\in K[t]$ により $\gcd(p,q)=ap+bq$ と書ける（:term[Bézout の等式]{en="Bézout's identity" yomi="べずーのとうしき"}）。
2. $p,q$ が互いに素であるための必要十分条件は、$ap+bq=1$ を満たす $a,b\in K[t]$ が存在することである。
3. $\gcd(p,q)=\gcd(p,r)=1$ ならば $\gcd(p,qr)=1$ である。
:::

:::proof
1. $I=\{ap+bq: a,b\in K[t]\}$ とおく。
$I$ は和について閉じ、任意の多項式を掛ける操作についても閉じている。
また $p,q\in I$ なので $I$ は零でない元を含む。
$I$ の零でない元のうち次数最小のものを一つ選び、最高次係数で割ってモニックな $d=ap+bq\in I$ を得る。
任意の $f\in I$ を :ref[lem-linalg-polynomial-division] の 2 で $f=sd+r$ と割ると、$r=f-sd\in I$ であり、$r\ne0$ ならば $\deg r<\deg d$ となって $d$ の選び方に反する。
よって $d$ は $I$ のすべての元を割り切り、特に $p$ と $q$ を割り切る。
一方、$c$ が $p,q$ の公約元で $p=cp',\ q=cq'$ ならば $d=c(ap'+bq')$ なので $c\mid d$ である。
したがって $d$ は最大公約元である。
二つの最大公約元は互いに割り切り合うモニック多項式なので、上の注意により一致する。

2. 互いに素ならば 1 から $1=ap+bq$ と書ける。
逆に $ap+bq=1$ ならば、公約元 $\gcd(p,q)$ は $1$ を割り切るので定数であり、モニックだから $1$ である。

3. $1=ap+bq$, $1=a'p+b'r$ と書いて掛け合わせると
$$
1=(aa'p+ab'r+ba'q)\,p+(bb')\,qr
$$
となるので、2 により $p$ と $qr$ は互いに素である。
:::

最大公約元は割り算の繰り返し（:term[Euclid の互除法]{en="Euclidean algorithm" yomi="ゆーくりっどのごじょほう"}）で計算できる。
$p=sq+r$ ならば、$p,q$ の公約元と $q,r$ の公約元は一致する（$r=p-sq$ と $p=sq+r$ による）ので、$\gcd(p,q)=\gcd(q,r)$ である。
余りの次数は真に減少するから、有限回で余りが零になり、最後の零でない余りをモニックにしたものが最大公約元である。
割り算の式を逆にたどれば Bézout の等式の係数 $a,b$ も求まる。

:::example[互除法の計算]
$p=t^3-1,\ q=t^2-1$ とする。
$t^3-1=t\,(t^2-1)+(t-1)$、$t^2-1=(t+1)(t-1)$ なので、余りの列は $t-1,0$ で終わり、$\gcd(p,q)=t-1$ である。
第一式から
$$
t-1=1\cdot(t^3-1)+(-t)\cdot(t^2-1)
$$
であり、これが Bézout の等式である。
:::

:::definition[既約多項式]{#def-linalg-irreducible}
次数 $1$ 以上の $f\in K[t]$ が、次数 $1$ 以上の二つの多項式の積に書けないとき、$f$ は $K$ 上 :term[既約]{en="irreducible" yomi="きやく"} であるという。
:::

既約多項式 $f$ の約元は、零でない定数と $f$ の定数倍だけである（$f=gh$ ならば $g,h$ の一方は定数）。
一次式はつねに既約である。
次数 $2$ または $3$ の多項式が可約ならば一次の因子をもつので、:ref[lem-linalg-factor-theorem] の 1 により $K$ に根をもつ。
したがって次数 $2,3$ の多項式が既約であることと $K$ に根をもたないことは同値である。
例えば $t^2+1$ は $\R$ 上既約だが、$\C$ 上では $(t-i)(t+i)$ と分解する。
既約性は係数体に依存する概念である。

:::theorem[既約分解の存在と一意性]{#thm-linalg-unique-factorization}
1. （Euclid の補題）$f$ が既約で $f\mid gh$ ならば、$f\mid g$ または $f\mid h$ である。より一般に $f\mid g_1\cdots g_k$ ならば、ある $i$ について $f\mid g_i$ である。
2. 次数 $1$ 以上のモニック多項式はモニック既約多項式の積に書け、その表示は因子の順序を除いて一意である。
3. $f_1,\ldots,f_r$ を相異なるモニック既約多項式とし、$p=\prod_{i=1}^rf_i^{s_i}$（$s_i\ge0$）とする。$p$ のモニックな約元は $\prod_if_i^{k_i}$（$0\le k_i\le s_i$）の形のものに限る。特に $f_i^{s_i}$ と $\prod_{j\ne i}f_j^{s_j}$ は互いに素である。
4. 零でない $p_1,\ldots,p_k$ が $p_j=c_j\prod_if_i^{s_{ij}}$（$c_j\in K$）と書けているとき、$l=\prod_if_i^{\max_js_{ij}}$ はすべての $p_j$ の倍元であり、$p_1,\ldots,p_k$ の任意の公倍元を割り切る。このようなモニックな $l$ は一意であり、:term[最小公倍元]{en="least common multiple" yomi="さいしょうこうばいげん"} と呼んで $\operatorname{lcm}(p_1,\ldots,p_k)$ と書く。
:::

:::proof
1. $f\nmid g$ とする。既約性から $\gcd(f,g)=1$ なので :ref[lem-linalg-bezout] により $af+bg=1$ と書ける。$h=afh+bgh$ の右辺は $f$ の倍数だから $f\mid h$ である。因子が三つ以上あるときはこの議論を繰り返す。

2. 既約ならそれ自身を取る。可約なら二つの正次数の因子に分け、両方をモニックにする。もとの最高次係数は $1$ なので積は変わらない。各因子の次数は下がるから、帰納的に既約分解が得られる。二つの分解があれば、一方の最初の因子は 1 によって他方のどれかを割る。どちらもモニックかつ既約だから一致する。:ref[lem-linalg-polynomial-division] でその因子を約分し、次数が零になるまで繰り返せば一意性が従う。

3. $p=de$ をモニックな因子に分ける。$d,e$ を既約分解したときの指数の和は 2 により $s_i$ である。逆に指定の範囲の $k_i$ には、残りの指数 $s_i-k_i$ で商を作れる。互いに素という結論は、共通する既約因子がないことから従う。

4. 各 $f_i$ の最大指数を取れば各 $p_j$ で割り切れる。どの非零公倍元にも少なくともこの指数が必要なので $l$ はそれを割る。零も $l$ の倍数である。同じ条件を満たす二つのモニック多項式は互いに割り切るから一致する。
:::

:::proposition[複素数体と実数体上の既約多項式]{#prop-linalg-irreducible-real-complex}
1. $\C$ 上のモニック既約多項式は一次式 $t-c$ に限る。したがって次数 $1$ 以上の複素係数多項式は一次式の積に分解する。
2. $\R$ 上のモニック既約多項式は、一次式 $t-c$ と、判別式が負の二次式 $t^2+\beta t+\gamma$（$\beta^2-4\gamma<0$）に限る。
:::

:::proof
代数学の基本定理により、正次数の複素多項式には複素根がある。その解析的な証明は [複素解析](/math/complex-analysis/) の Liouville の定理から得られる（:cite[ahlfors-complex]）。:ref[lem-linalg-factor-theorem] より根 $c$ があれば $t-c$ を取り出せるので、複素既約多項式の次数は $1$ である。

実係数なら非実根 $c$ と $\bar c$ は対で現れる。その積因子 $g=(t-c)(t-\bar c)$ は実係数で、判別式は $-4(\Im c)^2<0$ である。実係数で $f$ を $g$ で割った余りは次数が $2$ 未満で $c,\bar c$ の両方を根にもつので零である。したがって実既約多項式は実根をもてば一次、もたなければこの二次因子に等しい。逆に負の判別式をもつ実二次式は実根をもたず、一次式の積にはならない。
:::

### ベクトルを消す多項式と巡回部分空間

$v\in V$ を一つ固定し、$p(T)v=0$ を満たす $p$ を考える。
:ref[prop-linalg-minimal-polynomial-divides] の 1, 2 の証明は、$p(T)=0$ を $p(T)v=0$ に置き換えてもそのまま通用する（$\Phi_T(T)v=0$ なので零でないそのような $p$ が存在する）。
したがって、$p(T)v=0$ と $m_{T,v}\mid p$ が同値になるモニック多項式 $m_{T,v}$ がただ一つ存在する。
これを $v$ を消す $T$ の最小多項式と呼ぶ。
$m_T(T)v=0$ だから $m_{T,v}\mid m_T$ であり、$v=0$ のときに限り $m_{T,v}=1$ である。

モニック多項式 $f=t^d+a_{d-1}t^{d-1}+\cdots+a_1t+a_0$（$d\ge1$）に対し、$d$ 次正方行列

$$
C(f)=\begin{pmatrix}
0&0&\cdots&0&-a_0\\
1&0&\cdots&0&-a_1\\
0&1&\ddots&\vdots&-a_2\\
\vdots&\ddots&\ddots&0&\vdots\\
0&\cdots&0&1&-a_{d-1}
\end{pmatrix}
$$

を $f$ の :term[同伴行列]{en="companion matrix" yomi="どうはんぎょうれつ"} という。
対角のすぐ下の成分がすべて $1$、最後の列が $-a_0,\ldots,-a_{d-1}$ で、他は零である。

:::lemma[同伴行列の特性多項式と最小多項式]{#lem-linalg-companion}
$\Phi_{C(f)}=m_{C(f)}=f$ である。
:::

:::proof
$C=C(f)$ とおく。$C^je_1=e_{j+1}$（$0\le j<d$）なので次数 $d$ 未満の非零多項式は $e_1$ を消せない。一方、最後の列より $C^de_1=-\sum_{j<d}a_jC^je_1$、すなわち $f(C)e_1=0$ である。$f(C)$ は $C$ と可換だから全基底 $C^je_1$ を消す。したがって $m_C=f$ である。:ref[thm-linalg-cayley-hamilton] より $\Phi_C$ は $m_C$ の倍数であり、ともに次数 $d$ のモニック多項式だから一致する。
:::

:::proposition[巡回部分空間]{#prop-linalg-cyclic-subspace}
$v\ne0$ とし、$m_{T,v}=t^d+a_{d-1}t^{d-1}+\cdots+a_0$ とする。

1. $v,Tv,\ldots,T^{d-1}v$ は一次独立であり、$Z(v)\defeq\spn\{T^kv: k\ge0\}$ の基底である。
2. $Z(v)$ は $v$ を含む最小の $T$ 不変部分空間である。
3. 基底 $(v,Tv,\ldots,T^{d-1}v)$ に関する $T|_{Z(v)}$ の表現行列は $C(m_{T,v})$ であり、$T|_{Z(v)}$ の特性多項式と最小多項式はともに $m_{T,v}$ に等しい。
:::

:::proof
1. 非自明な関係 $\sum_{j<d}c_jT^jv=0$ があれば、$p=\sum_{j<d}c_jt^j\ne0$ は $p(T)v=0$ を満たすから $m_{T,v}\mid p$ となるが、$\deg p<d$ なので不可能である。
また $T^dv=-\sum_{j<d}a_jT^jv$ は $U=\spn\{v,\ldots,T^{d-1}v\}$ に属する。
$T(U)=\spn\{Tv,\ldots,T^dv\}\subset U$ なので、帰納的にすべての $T^kv$ が $U$ に属し、$Z(v)=U$ である。

2. 1 の証明から $Z(v)$ は $T$ 不変である。
$v$ を含む $T$ 不変部分空間はすべての $T^kv$ を含むので $Z(v)$ を含む。

3. $T(T^jv)=T^{j+1}v$（$j<d-1$）と $T(T^{d-1}v)=T^dv=-\sum_ja_jT^jv$ から、表現行列は $C(m_{T,v})$ である。
後半は :ref[lem-linalg-companion] による。
:::

$Z(v)$ を $v$ の生成する :term[巡回部分空間]{en="cyclic subspace" yomi="じゅんかいぶぶんくうかん"} という。
$Z(v)=V$ となる $v$ を $T$ の :term[巡回ベクトル]{en="cyclic vector" yomi="じゅんかいべくとる"} という。

:::proposition[最小多項式と最小公倍元]{#prop-linalg-minpoly-lcm}
1. $(v_1,\ldots,v_n)$ を $V$ の基底とすると $m_T=\operatorname{lcm}(m_{T,v_1},\ldots,m_{T,v_n})$ である。
2. $V=W_1\oplus\cdots\oplus W_k$ で各 $W_i$ が $T$ 不変ならば $m_T=\operatorname{lcm}(m_{T|_{W_1}},\ldots,m_{T|_{W_k}})$ である。
:::

:::proof
1. $p(T)=0$ はすべての $i$ について $p(T)v_i=0$ であることと同値であり、それは $m_{T,v_i}\mid p$（すべての $i$）、すなわち :ref[thm-linalg-unique-factorization] の 4 により $\operatorname{lcm}\mid p$ と同値である。
$p=m_T$ と $p=\operatorname{lcm}$ に適用すれば両者は互いに割り切り合うので一致する。

2. $p(T)$ は各 $W_i$ を保ち、その制限は $p(T|_{W_i})$ である。
任意のベクトルは各 $W_i$ の成分の和なので、$p(T)=0$ とすべての $i$ について $p(T|_{W_i})=0$ は同値である。
あとは 1 と同じ議論である。
:::

:::example[一つのベクトルだけでは足りない場合]
$T=\diag(1,2,2)$ に対し、$v=e_1$ なら $m_{T,v}=t-1$、$v=e_2$ なら $m_{T,v}=t-2$ である。
$v=e_1+e_2$ なら両成分を消す必要があるので $m_{T,v}=(t-1)(t-2)$ であり、これは $m_T=\operatorname{lcm}(t-1,t-2,t-2)$ に一致する。
一方 $\deg m_T=2<3$ なので、$T$ に巡回ベクトルは存在しない（:ref[prop-linalg-cyclic-subspace] の 1 から $\dim Z(v)\le2$）。
:::

:::example[一列の計算から行列全体の多項式を得る]
$T(x,y,z)=(z,x+z,y)$、$v=e_1$ とする。
$Tv=e_2,\ T^2v=e_3,\ T^3v=Te_3=e_1+e_2=v+Tv$ である。
よって $m_{T,v}=t^3-t-1$ である。
$m_{T,v}\mid m_T\mid\Phi_T$ で $\deg\Phi_T=3$ なので、三つはすべて等しく、$e_1$ は巡回ベクトルである。
特に $T^3=T+I$ が全空間で成り立ち、例えば $T^5=T^2\cdot T^3=T^3+T^2=T^2+T+I$ となる。
多項式の複素根を求める必要はない。
:::

### 特性多項式との関係

最小多項式は特性多項式を割り切るが、逆向きにも強い関係がある。

:::proposition[特性多項式と最小多項式の既約因子]{#prop-linalg-charpoly-minpoly-factors}
$\Phi_T\mid m_T^{\,n}$ である。
したがって $\Phi_T$ と $m_T$ は同じモニック既約因子をもつ。
すなわち $\Phi_T=\prod_{i=1}^rf_i^{e_i}$（$f_i$ は相異なるモニック既約多項式、$e_i\ge1$）ならば、$m_T=\prod_{i=1}^rf_i^{s_i}$、$1\le s_i\le e_i$ である。
特に $\Phi_T$ と $m_T$ は $K$ において同じ根をもつ。
:::

:::proof
基底を選び $A=[T]_{\mathcal B}$ とする（:ref[prop-linalg-minimal-polynomial-divides] の 4 より $m_T=m_A$）。
$m=m_A=\sum_{k=0}^dc_kt^k$ と書く。
各 $k\ge1$ について、$M_n(K[t])$ において $tI$ と $A$ は可換なので

$$
t^kI-A^k=(tI-A)\bigl(t^{k-1}I+t^{k-2}A+\cdots+A^{k-1}\bigr)
$$

が成り立つ（右辺を展開すると隣り合う項が打ち消し合う）。
$m(A)=O$ を引けば

$$
m(t)I=m(t)I-m(A)=\sum_{k=1}^dc_k\bigl(t^kI-A^k\bigr)=(tI-A)B(t)
$$

となる $B(t)\in M_n(K[t])$ がある。
両辺の行列式を取ると、:ref[prop-linalg-det-ring] の 2（$K[t]$ 上の乗法性）により

$$
m(t)^n=\det(tI-A)\det B(t)=\Phi_A(t)\det B(t)
$$

であり、$\Phi_A\mid m^n$ を得る。
後半：$f$ が $\Phi_T$ の既約因子ならば $f\mid m_T^n$ であり、:ref[thm-linalg-unique-factorization] の 1 から $f\mid m_T$ である。
逆に $m_T\mid\Phi_T$ なので、$m_T$ の既約因子は $\Phi_T$ の既約因子である。
指数の不等式は :ref[thm-linalg-unique-factorization] の 3 による。
根については :ref[lem-linalg-factor-theorem] の 1 により「$t-c$ が既約因子であること」と言い換えればよい。
:::

:::example[最小多項式の候補を絞る]
$\Phi_A=(t-1)^2(t-2)$ ならば、$m_A$ は $(t-1)(t-2)$ か $(t-1)^2(t-2)$ のいずれかである。
どちらであるかは $(A-I)(A-2I)=O$ かどうかを一度計算すれば判定できる。
$\Phi_A=(t^2+1)^2$ の実行列ならば、$m_A$ は $t^2+1$ か $(t^2+1)^2$ である。
:::

### 多項式の逆と固有値

$m_T$ による割り算は $p(T)$ の逆写像も与える。

:::proposition[多項式に代入した作用素の可逆性]{#prop-linalg-polynomial-invertible}
$p\in K[t]$ について、$p(T)$ が可逆であることと $\gcd(p,m_T)=1$ は同値である。
このとき $p(T)^{-1}$ は $T$ の多項式として書ける。
特に $\lambda\in K$ が $T$ の固有値であることと $m_T(\lambda)=0$ は同値である。
:::

:::proof
$\gcd(p,m_T)=1$ ならば :ref[lem-linalg-bezout] により $ap+bm_T=1$ と書け、$T$ を代入すると $a(T)p(T)=p(T)a(T)=I$ である。
逆に $d=\gcd(p,m_T)$ が次数 $1$ 以上とし、$m_T=dh$ と書く。
$\deg h<\deg m_T$ なので $h(T)\ne0$ であり、$h(T)x\ne0$ となる $x$ がある。
$p=dp'$ と書けば $p(T)h(T)=p'(T)m_T(T)=0$ なので $p(T)\bigl(h(T)x\bigr)=0$ であり、$p(T)$ は単射でない。
最後の主張：$t-\lambda$ のモニックな約元は $1$ と $t-\lambda$ だけなので、$\gcd(t-\lambda,m_T)=1$ は $(t-\lambda)\nmid m_T$、すなわち :ref[lem-linalg-factor-theorem] の 1 により $m_T(\lambda)\ne0$ と同値である。
一方、$T-\lambda I$ が可逆でないことは、:ref[cor-linalg-equal-dimension-map] により単射でないこと、すなわち $\lambda$ が固有値であることと同値である。
:::

:::example[割り算で逆を作る]
$T^2=-I$ のとき $(I+T)(I-T)=I-T^2=2I$ であり、$K=\R$ 上では $(I+T)^{-1}=(I-T)/2$ である。
多項式で書けば $(1+t)\cdot\frac{1-t}{2}+\frac12(t^2+1)=1$ であり、これは $p=1+t$ と $m_T=t^2+1$ に対する Bézout の等式である。
:::

## 不変部分空間と広義固有空間

$T$ の作用を小さい空間ごとに調べるには、:ref[def-linalg-invariant-subspace] の $T$ 不変部分空間を用いる。
$V$ が $T$ 不変部分空間の直和 $W_1\oplus\cdots\oplus W_k$ に分かれれば、:ref[prop-linalg-invariant-block] により、各 $W_i$ の基底を並べた基底で $T$ の表現行列はブロック対角になり、$T$ の研究は各制限 $T|_{W_i}$ の研究に帰着する。
任意の $p\in K[t]$ について $p(T)$ は $T$ と可換なので、$\Ker p(T)$ と $\Img p(T)$ は $T$ 不変である（:ref[def-linalg-invariant-subspace] の直後の説明）。
多項式の核は、このような分解を組織的に与える。

以下、零空間上の作用素の最小多項式は $1$ と約束する。

:::lemma[制限の最小多項式]{#lem-linalg-restriction-minpoly}
$W$ が $T$ 不変部分空間ならば $m_{T|_W}\mid m_T$ である。
:::

:::proof
$(T|_W)^k=T^k|_W$ なので、任意の $p$ について $p(T|_W)=p(T)|_W$ である。
特に $m_T(T|_W)=m_T(T)|_W=0$ であり、:ref[prop-linalg-minimal-polynomial-divides] の 2 から主張が従う（$W=\{0\}$ なら $m_{T|_W}=1$ で明らか）。
:::

特性多項式についても :ref[cor-linalg-invariant-charpoly] により $\Phi_{T|_W}\mid\Phi_T$ である。

:::definition[広義固有空間]{#def-linalg-generalized-eigenspace}
$\lambda\in K$ に対し

$$
\widetilde E_\lambda(T)=\{v\in V: \text{ある } j\ge1 \text{ について } (T-\lambda I)^jv=0\}=\bigcup_{j\ge1}\Ker(T-\lambda I)^j
$$

を $T$ の $\lambda$ に対する :term[広義固有空間]{en="generalized eigenspace" yomi="こうぎこゆうくうかん"} という。
その零でない元を :term[広義固有ベクトル]{en="generalized eigenvector" yomi="こうぎこゆうべくとる"} という。
:::

和集合で定義したが、次の命題により、これは一つの冪の核に等しく、特に部分空間である。

:::proposition[核の列と像の列の安定]{#prop-linalg-kernel-stabilization}
$S\in\End(V)$ とする。

1. $\Ker S^j\subset\Ker S^{j+1}$、$\Img S^j\supset\Img S^{j+1}$ である。
2. ある $k$ で $\Ker S^k=\Ker S^{k+1}$ ならば、すべての $j\ge0$ で $\Ker S^{k+j}=\Ker S^k$ かつ $\Img S^{k+j}=\Img S^k$ である。
3. $\Ker S^n=\Ker S^{n+j}$、$\Img S^n=\Img S^{n+j}$（$j\ge0$）である。特に $\widetilde E_\lambda(T)=\Ker(T-\lambda I)^n$ であり、これは固有空間 $E_\lambda(T)$ を含む $T$ 不変部分空間である。
:::

:::proof
1. $S^jx=0$ ならば $S^{j+1}x=0$、また $S^{j+1}x=S^j(Sx)$ である。

2. $x\in\Ker S^{k+2}$ ならば $Sx\in\Ker S^{k+1}=\Ker S^k$ なので $x\in\Ker S^{k+1}$ である。
よって $\Ker S^{k+2}=\Ker S^{k+1}$ であり、帰納法で核についての主張を得る。
:ref[thm-linalg-rank-nullity] により $\dim\Img S^j=n-\dim\Ker S^j$ なので、像の次元も $k$ 以降一定であり、1 の包含から像も一定である。

3. 核の次元列 $0=\dim\Ker S^0\le\dim\Ker S^1\le\cdots$ は $n$ 以下であり、2 によりいったん等号が成り立てば以後一定である。
したがって最初に等号が成り立つまでは真に増加し、その段階は $n$ 以下である。
最後の主張は $S=T-\lambda I$ に適用し、$T$ と $(T-\lambda I)^n$ が可換であることによる。
:::

:::proposition[Fitting 分解]{#prop-linalg-fitting}
$S\in\End(V)$ について $V=\Ker S^n\oplus\Img S^n$ である。
両成分は $S$ 不変であり、$S$ の制限は第一成分上で冪零、第二成分上で可逆である。
:::

:::proof
$x\in\Ker S^n\cap\Img S^n$ とし、$x=S^ny$ と書く。
$S^{2n}y=S^nx=0$ なので、:ref[prop-linalg-kernel-stabilization] の 3 より $y\in\Ker S^{2n}=\Ker S^n$、したがって $x=S^ny=0$ である。
共通部分が零なので :ref[prop-linalg-direct-sum-criterion] の 1 により和は直和である。
:ref[thm-linalg-rank-nullity] から二つの次元の和は $n$ なので、:ref[cor-linalg-subspace-dimension] によりこの直和は $V$ 全体である。
両空間は $S^n$ の核と像なので $S$ 不変であり、核への制限の $n$ 乗は零である。
像の側では $S(\Img S^n)=\Img S^{n+1}=\Img S^n$ なので、制限は全射であり、:ref[cor-linalg-equal-dimension-map] により可逆である。
:::

:::example[スペクトルを求めずに一成分だけ分ける]
次の行列から、固有値 $2$ に対応する成分を取り出す。

$$
A=\begin{pmatrix}2&1&0\\0&2&0\\0&0&-1\end{pmatrix},\qquad S=A-2I.
$$

$S^2=\diag(0,0,9)$ であり、$S^3=\diag(0,0,-27)$ と比べて核も像もすでに安定している。
核は最初の二座標の平面、像は第三座標の直線である。
$S$ は平面上で二乗すると零になり、直線上では $-3$ 倍である。
広義固有空間を一つだけ取り出すには、他の固有値に対する計算をしなくてよい。
:::

### 互いに素な因子による分解

:::lemma[互いに素な因子による分解]{#lem-linalg-coprime-splitting}
$p,q\in K[t]$ が互いに素ならば

$$
\Ker\,(pq)(T)=\Ker p(T)\oplus\Ker q(T)
$$

であり、この直和分解の各成分への射影は $T$ の多項式の制限で与えられる。
特に $p(T)q(T)=0$ ならば $V=\Ker p(T)\oplus\Ker q(T)$ である。
この補題の証明は $V$ の次元を用いないので、無限次元の空間でも成り立つ。
:::

:::proof
:ref[lem-linalg-bezout] により $ap+bq=1$ となる $a,b$ がある。
$T$ を代入すると、すべての $x$ について

$$
x=b(T)q(T)x+a(T)p(T)x
$$

である。
$x\in\Ker(pq)(T)$ ならば、$p(T)\bigl(b(T)q(T)x\bigr)=b(T)(pq)(T)x=0$ なので第一項は $\Ker p(T)$ に、同様に第二項は $\Ker q(T)$ に属する。
逆に $\Ker p(T)$ と $\Ker q(T)$ はともに $\Ker(pq)(T)$ に含まれる。
$x$ が両方の核に属すれば、上の式の両項が零になるので $x=0$ である。
よって直和であり、$\Ker p(T)$ への射影は $b(T)q(T)$ の制限、$\Ker q(T)$ への射影は $a(T)p(T)$ の制限である。
:::

:::theorem[一次分解定理]{#thm-linalg-primary-decomposition}
$m_T=\prod_{i=1}^rf_i^{s_i}$ を、相異なるモニック既約多項式 $f_i$ の冪（$s_i\ge1$）への分解とし、$V_i=\Ker f_i(T)^{s_i}$ とおく。
また $\Phi_T=\prod_{i=1}^rf_i^{e_i}$ とする（:ref[prop-linalg-charpoly-minpoly-factors]）。

1. $V=V_1\oplus\cdots\oplus V_r$ であり、各 $V_i$ は $T$ 不変である。
2. 1 の直和分解に関する $V_i$ への射影 $P_i$（$V_i$ 成分を取り出し、他の成分を消す写像）は $T$ の多項式である。
3. $m_{T|_{V_i}}=f_i^{s_i}$ であり、特に $V_i\ne\{0\}$ である。
4. $\Phi_{T|_{V_i}}=f_i^{e_i}$ であり、$\dim V_i=e_i\deg f_i$ である。
:::

:::proof
$p_i=f_i^{s_i}$、$h_i=\prod_{j\ne i}p_j=m_T/p_i$ とおく。
:ref[thm-linalg-unique-factorization] の 3 により $p_i$ と $h_i$ は互いに素なので、:ref[lem-linalg-bezout] により $a_ih_i+b_ip_i=1$ となる $a_i,b_i$ がある。
$\eta_i=a_ih_i$ とおく。

1, 2. $j\ne i$ ならば $p_j\mid h_i$ なので $p_j\mid \eta_i$ であり、また $\eta_i-1=-b_ip_i$ は $p_i$ で割り切れる。
したがって $\sum_i\eta_i-1$ はすべての $p_j$ で割り切れ、:ref[thm-linalg-unique-factorization] の 4 によりそれらの最小公倍元 $\prod_jp_j=m_T$ で割り切れる。
よって $\sum_i\eta_i(T)=I$ である。
また $p_i(T)\eta_i(T)=a_i(T)m_T(T)=0$ より $\Img \eta_i(T)\subset V_i$ である。
したがって任意の $x$ は $x=\sum_i\eta_i(T)x$（$\eta_i(T)x\in V_i$）と書け、$V=\sum_iV_i$ である。
表示の一意性を示す。
$x_i\in V_i$ が $\sum_ix_i=0$ を満たすとする。
$\eta_j(T)$ を作用させると、$i\ne j$ では $p_i\mid \eta_j$ より $\eta_j(T)x_i=0$、$i=j$ では $\eta_j(T)x_j=x_j-b_j(T)p_j(T)x_j=x_j$ なので、$x_j=0$ を得る。
よって直和であり、$V_i$ への射影は $P_i=\eta_i(T)$ である。
$V_i$ は多項式の核なので $T$ 不変である。

3. $\mu_i=m_{T|_{V_i}}$ とおく（$V_i=\{0\}$ なら $\mu_i=1$）。
$p_i(T|_{V_i})=0$ なので $\mu_i\mid p_i$ である。
一方、$x=\sum_ix_i$（$x_i\in V_i$）に対し $\bigl(\prod_j\mu_j\bigr)(T)x_i=\bigl(\prod_{j\ne i}\mu_j\bigr)(T)\,\mu_i(T)x_i=0$ なので、$\prod_j\mu_j$ は $T$ を消し、$m_T\mid\prod_j\mu_j$ である。
次数を比べると
$$
\deg m_T\le\sum_j\deg\mu_j\le\sum_j\deg p_j=\deg m_T
$$
であり、すべての $j$ で $\deg\mu_j=\deg p_j$ となる。
$\mu_j$ は $p_j$ のモニックな約元で次数が等しいので $\mu_j=p_j$ である。
$\deg p_j\ge1$ なので $V_j\ne\{0\}$ である。

4. 1 と :ref[prop-linalg-invariant-block] により $T$ の表現行列はブロック対角であり、:ref[cor-linalg-invariant-charpoly] を繰り返し用いて $\Phi_T=\prod_i\Phi_{T|_{V_i}}$ を得る。
:ref[prop-linalg-charpoly-minpoly-factors] を $T|_{V_i}$ に適用すると、$\Phi_{T|_{V_i}}$ の既約因子は $\mu_i=f_i^{s_i}$ の既約因子 $f_i$ だけなので、$\Phi_{T|_{V_i}}=f_i^{c_i}$ と書ける。
$\prod_if_i^{c_i}=\prod_if_i^{e_i}$ と既約分解の一意性から $c_i=e_i$ であり、次数を比べて $\dim V_i=e_i\deg f_i$ である。
:::

この分解を :term[一次分解]{en="primary decomposition" yomi="いちじぶんかい"}（準素分解）という。
実数上で既約因子が二次式になる場合も含めて、任意の体で成り立つことに注意する。

:::corollary[広義固有空間分解]{#cor-linalg-generalized-eigenspace-decomposition}
$\Phi_T=\prod_{i=1}^r(t-\lambda_i)^{m_i}$（$\lambda_i$ は相異なる固有値、$m_i$ は代数的重複度）と分解するとき、$m_T=\prod_i(t-\lambda_i)^{s_i}$、$1\le s_i\le m_i$ であり、次が成り立つ。

1. $\widetilde E_{\lambda_i}(T)=\Ker(T-\lambda_iI)^{s_i}=\Ker(T-\lambda_iI)^{m_i}$ であり、$\dim\widetilde E_{\lambda_i}(T)=m_i$ である。
2. $V=\widetilde E_{\lambda_1}(T)\oplus\cdots\oplus\widetilde E_{\lambda_r}(T)$ である。
3. $\widetilde E_{\lambda_i}(T)$ 上で $T-\lambda_iI$ は冪零で、$(T-\lambda_iI)^{s_i}=0$、$(T-\lambda_iI)^{s_i-1}\ne0$ である。$\mu\ne\lambda_i$ ならば $T-\mu I$ は $\widetilde E_{\lambda_i}(T)$ 上で可逆である。
4. $\mu$ が固有値でなければ $\widetilde E_\mu(T)=\{0\}$ である。
:::

:::proof
$f_i=t-\lambda_i$ として :ref[thm-linalg-primary-decomposition] を用いる。
$V_i=\Ker(T-\lambda_iI)^{s_i}\subset\widetilde E_{\lambda_i}(T)$ は明らかである。
3 の後半を先に示す。
$T|_{V_k}$ の最小多項式は $(t-\lambda_k)^{s_k}$ であり、$\mu\ne\lambda_k$ ならば $t-\mu$ と互いに素なので、:ref[prop-linalg-polynomial-invertible] により $T-\mu I$ は $V_k$ 上で可逆である。
さて $x\in\widetilde E_{\lambda_i}(T)$、$(T-\lambda_iI)^jx=0$ とし、$x=\sum_kx_k$（$x_k\in V_k$）と分ける。
$(T-\lambda_iI)^j$ は各 $V_k$ を保つので $(T-\lambda_iI)^jx_k=0$ であり、$k\ne i$ では可逆性から $x_k=0$ である。
よって $x=x_i\in V_i$ であり、$\widetilde E_{\lambda_i}(T)=V_i$ を得る。
$s_i\le m_i$ なので $V_i\subset\Ker(T-\lambda_iI)^{m_i}\subset\widetilde E_{\lambda_i}(T)$ であり、1 の等式が従う。
次元は定理の 4 による。
2 は定理の 1 である。
3 の前半は $m_{T|_{V_i}}=(t-\lambda_i)^{s_i}$ の言い換えである。
4：$\mu$ が固有値でなければ $T-\mu I$ は単射なので、その冪も単射である。
:::

:ref[prop-linalg-eigenspaces-direct-sum] では固有空間の和が直和になることを示したが、それが全空間になるのは対角化可能な場合に限られた。
広義固有空間に広げれば、特性多項式が分解する限り全空間が埋め尽くされる。

:::example[射影多項式を作る]
$m_T=(t-1)^2(t-3)$ とする。
$(t-1)^2$ は $t=3$ で $4$ だから、$P_3=(T-I)^2/4$ は固有値 $3$ の成分で恒等、固有値 $1$ の成分で零になる（定理の証明の $\eta_i$ を直接見つけたことになる）。
よって $P_1=I-P_3$ であり、各ベクトルの広義固有成分を、基底を求めずに取り出せる。

例えば次の行列に対し $e_3$ を分解する。

$$
T=\begin{pmatrix}1&1&1\\0&1&0\\0&0&3\end{pmatrix}
$$

$(T-I)(T-3I)$ の $(1,2)$ 成分は $-2$ なので $m_T\ne(t-1)(t-3)$ であり、:ref[prop-linalg-charpoly-minpoly-factors] から $m_T=\Phi_T=(t-1)^2(t-3)$ である。
$(T-I)e_3=(1,0,2)$ にもう一度 $T-I$ を掛けると $(2,0,4)$ になる。
よって $P_3e_3=(1/2,0,1)$、$P_1e_3=(-1/2,0,0)$ となる。
実際 $T(1/2,0,1)=(3/2,0,3)$ であり、前者は固有値 $3$ の固有ベクトルである。
:::

### 最小多項式による対角化可能性の判定

:::theorem[最小多項式による対角化可能性の判定]{#thm-linalg-minpoly-diagonalizable}
1. 次の三条件は同値である。
   (a) $T$ は $K$ 上対角化可能である。
   (b) $m_T$ は $K$ 上の相異なる一次式の積である。
   (c) 相異なる $c_1,\ldots,c_k\in K$ について $\prod_{i=1}^k(T-c_iI)=0$ が成り立つ。
2. $T$ が $K$ 上三角化可能であることと、$m_T$ が $K$ 上一次式の積に分解することは同値である。
:::

:::proof
固有基底があれば、異なる固有値 $\lambda_i$ について $p=\prod_i(t-\lambda_i)$ を取る。:ref[lem-linalg-polynomial-substitution] により各基底ベクトル $v$ で $p(T)v=p(\lambda_i)v=0$ だから (c) が従う。(c) の消去多項式を $m_T$ が割るので、:ref[lem-linalg-factor-theorem] の約元の記述から (b) が従う。(b) のとき :ref[thm-linalg-primary-decomposition] の各成分は $\Ker(T-\lambda_iI)$ そのものであり、基底を合わせれば (a) を得る。

三角化可能性は :ref[thm-linalg-triangularization] により $\Phi_T$ が一次式の積になることと同値である。:ref[prop-linalg-charpoly-minpoly-factors] より $m_T$ と $\Phi_T$ は同じ既約因子をもつので、一方が一次式の積なら他方もそうである。
:::

この判定は、固有ベクトルを一本も求めずに使えるのが利点である。
例えば $T^2=T$ ならば $T$ は $t(t-1)$ で消えるので、任意の体で対角化可能である。
$K$ の標数が $2$ でなければ（$1\ne-1$）、$T^2=I$ を満たす $T$ も $(t-1)(t+1)$ で消えるので対角化可能である。
標数 $2$ の体 $\F_2$ では $1=-1$ であり、$\begin{pmatrix}1&1\\0&1\end{pmatrix}^2=I$ だがこの行列は対角化できない。
また :ref[lem-linalg-restriction-minpoly] と :ref[lem-linalg-factor-theorem] の 5 により、対角化可能な作用素の不変部分空間への制限は対角化可能である。

### 可換な対角化可能作用素の同時対角化

[正規作用素とスペクトル定理](/math/linear-algebra/spectral-theorem/) の :ref[thm-linalg-simultaneous] では、正規直交基底で対角化するために正規性を仮定した。
一般の基底でよければ、可換性だけが本質的な条件である。

:::theorem[可換な対角化可能作用素の同時対角化]{#thm-linalg-simultaneous-diagonalization}
$T_1,\ldots,T_k\in\End(V)$ がそれぞれ $K$ 上対角化可能であるとき、すべての $T_j$ を同時に対角化する基底が存在するための必要十分条件は、どの二つも可換であることである。
:::

:::proof
同じ基底で対角なら可換である。逆に互いに可換とし、作用素の個数 $k$ で帰納する。一本なら仮定そのものである。$T_1$ の固有空間 $E$ は $T_1T_jx=T_jT_1x$ により全 $T_j$ で不変である。各制限の最小多項式はもとの最小多項式を割るので、直前の判定により制限も対角化可能である。$T_2,\ldots,T_k$ の制限は互いに可換だから、帰納仮定で $E$ の共通固有基底を選べる。これらをすべての $E$ で合わせれば $V$ の基底となり、$T_1$ も各 $E$ ではスカラー倍である。
:::

:::example[共通固有基底が直交しない場合]
次の行列とその二乗は可換である。

$$
A=\begin{pmatrix}1&1\\0&2\end{pmatrix},\qquad B=A^2=\begin{pmatrix}1&3\\0&4\end{pmatrix}
$$

$e_1$ と $(1,1)$ は $A$ では順に固有値 $1,2$、$B$ では $1,4$ に属する固有ベクトルで、両者を同時対角化する。
しかし $A$ の固有空間はどちらも一次元で、この二方向の内積は零でない。
したがって標準内積に関する共通の正規直交固有基底は存在しない。
$A$ が正規でないことと整合している。
:::

## 冪零作用素の構造

:ref[cor-linalg-generalized-eigenspace-decomposition] により、特性多項式が分解するとき、各広義固有空間の上で $N=T-\lambda I$ は冪零になる。
したがって $T$ の構造を決めるには、冪零作用素の構造を決めればよい。
$N$ を繰り返し作用させたとき、どの段階でどれだけのベクトルが消えるかを記録すれば、基底を選ぶ手順が見えてくる。

:::definition[冪零作用素と鎖]{#def-linalg-nilpotent-chain}
$N\in\End(V)$ について、$N^s=0$ となる正整数 $s$ が存在するとき $N$ を :term[冪零作用素]{en="nilpotent operator" yomi="べきれいさようそ"} といい、そのような最小の $s$ を :term[冪零指数]{en="index of nilpotency" yomi="べきれいしすう"} という。
$v\in V$ が $N^kv=0,\ N^{k-1}v\ne0$ を満たすとき、列 $(N^{k-1}v,\ldots,Nv,v)$ を $v$ を頂点とする長さ $k$ の :term[鎖]{en="chain" yomi="くさり"} と呼び、$N^{k-1}v$ をその末端という。
:::

鎖では、末端は $N$ で消え、それ以外の各ベクトルはすぐ左のベクトルに移る。
末端は $\Ker N$ の零でない元である。

:::proposition[冪零作用素の特徴づけ]{#prop-linalg-nilpotent-basic}
$N\in\End(V)$ について次は同値である。

1. $N$ は冪零である。
2. ある $s\ge1$ について $m_N=t^s$ である。このとき $s$ は冪零指数に等しい。
3. $\Phi_N=t^n$ である。
4. $N^n=0$ である。

特に冪零指数は $n$ 以下であり、冪零作用素の固有値は $0$ だけである。
:::

:::proof
1⇒2：$N^k=0$ ならば $m_N\mid t^k$ であり、:ref[thm-linalg-unique-factorization] の 3 により $m_N=t^s$（$s\ge1$）である。
$t^j$ が消去多項式であることは $m_N=t^s\mid t^j$、すなわち $s\le j$ と同値なので、$s$ は冪零指数である。
2⇒3：:ref[prop-linalg-charpoly-minpoly-factors] により $\Phi_N$ の既約因子は $t$ だけであり、次数 $n$ のモニック多項式なので $\Phi_N=t^n$ である。
3⇒4 は :ref[thm-linalg-cayley-hamilton]、4⇒1 は定義による。
固有値は $\Phi_N=t^n$ の根である。
:::

一本の鎖は一次独立である。実際、非自明な関係 $\sum_j a_jN^jv=0$ で $a_j\ne0$ の最小の添字を $h$ とする。長さを $k$ とすると $N^{k-1-h}$ を作用させて $a_hN^{k-1}v=0$ となり、末端が非零であることに反する。複数の鎖を合わせるには末端どうしの独立性も確保する必要がある。

:::theorem[冪零鎖による基底]{#thm-linalg-nilpotent-chains}
$V\ne\{0\}$ 上の冪零作用素 $N$ は、いくつかの鎖を並べた基底をもつ。
$K_j=\Ker N^j$、$d_j=\dim K_j$（$d_0=0$）とおくと、鎖を並べたどの基底においても、長さが $j$ 以上の鎖の本数は $d_j-d_{j-1}$、長さがちょうど $j$ の鎖の本数は $2d_j-d_{j-1}-d_{j+1}$ である。
特に鎖の本数は $\dim\Ker N$ に等しい。
:::

:::proof
$s$ を冪零指数とし、$1\le j\le s$ に対して $L_j=N^{j-1}K_j$、また $L_{s+1}=\{0\}$ とおく。
$N\cdot N^{j-1}K_j=N^jK_j=\{0\}$ なので $L_j\subset K_1$ であり、$L_j$ は「長さ $j$ 以上の鎖の末端になれるベクトル」の空間である。
$NK_{j+1}\subset K_j$ より $L_{j+1}=N^{j-1}(NK_{j+1})\subset N^{j-1}K_j=L_j$ である。
$N^{j-1}$ の $K_j$ への制限は $K_j$ から $L_j$ への全射で、その核は $K_j\cap\Ker N^{j-1}=K_{j-1}$ なので、:ref[thm-linalg-rank-nullity] から

$$
\dim L_j=d_j-d_{j-1}
$$

である。

$L_s$ の基底から始め、:ref[thm-linalg-basis-extension] により $L_{s-1},\ldots,L_1$ の基底へ順に延長する。
$L_j$ の段階で新たに加えた各末端 $u$ に対し、$v\in K_j$ を $N^{j-1}v=u$ となるよう選び、長さ $j$ の鎖 $(N^{j-1}v,\ldots,v)$ を作る（$u\ne0$ なので $N^{j-1}v\ne0$、$v\in K_j$ なので $N^jv=0$）。
こうして得たすべての末端を合わせたものは $L_1$ の基底なので一次独立である。

すべての鎖のベクトルを合わせた列が一次独立であることを示す。
各鎖において、末端から数えて $i$ 番目のベクトル（長さ $k$ の鎖なら $N^{k-i}v$）の高さを $i$ と定める。
非自明な一次関係があったとし、零でない係数が付いたベクトルの高さの最大値を $h$ とする。
関係に $N^{h-1}$ を作用させると、高さ $i<h$ のベクトル $N^{k-i}v$ は $N^{k-i+h-1}v=0$ となって消え、高さ $h$ のベクトルは $N^{k-1}v$、すなわち自分の鎖の末端に移る。
異なる鎖の末端は異なるので、末端たちの非自明な一次関係が得られ、末端の一次独立性に反する。

最後に本数を数える。
長さが $j$ 以上の鎖の末端は、ちょうど $L_s,\ldots,L_j$ の段階で選んだもの、すなわち $L_j$ の基底であるから、その本数は $\dim L_j=d_j-d_{j-1}$ である。
各鎖は長さの数だけの段階 $j=1,2,\ldots$ で一度ずつ数えられるので、ベクトルの総数は

$$
\sum_{j=1}^s(d_j-d_{j-1})=d_s=\dim V
$$

となる。
よって一次独立な列は基底である。
隣り合う二つの本数の差を取れば、長さがちょうど $j$ の本数の公式を得る（$j\ge s$ では $d_j=d_s$ と読む）。
任意の鎖基底でも、長さ $k$ の鎖で $N^j$ の核に入る基底ベクトルは末端から $\min(j,k)$ 本である。残りの像は互いに異なる基底ベクトルなので、それらの一次結合が新たに核に入ることもない。したがって $d_j=\sum_{\text{鎖}}\min(j,k)$ であり、隣り合う項の差を取ると同じ本数の公式を得る。
:::

長さが $j$ 以上の鎖の本数 $p_j=d_j-d_{j-1}$ は $L_{j+1}\subset L_j$ から $p_1\ge p_2\ge\cdots$ を満たす。
鎖の長さを大きい順に $k_1\ge k_2\ge\cdots$ と並べると、これは $n$ の分割（和が $n$ の正整数の非増加列）であり、$p_j$ は「長さが $j$ 以上の項の個数」である。
各鎖を一列の箱として並べた図（Young 図形）で言えば、$p_j$ は第 $j$ 段にある箱の個数である。

:::example[核の増分から読む]
次元 $7$ で $d_0,d_1,d_2,d_3,d_4=0,3,5,6,7$ とする。
増分 $3,2,1,1$ は、長さが $1,2,3,4$ 以上の鎖の本数である。
したがって鎖は三本で、長さは $4,2,1$ となる。
和が元の次元 $7$ に戻ることも確かめられる。
:::

:::example[鎖を手で選ぶ]
$Ne_3=e_1+e_2,\ Ne_2=e_1,\ Ne_1=0$ とする。
$N^2e_3=N(e_1+e_2)=e_1$ なので $(e_1,e_1+e_2,e_3)$ が長さ $3$ の鎖である。
$e_1,e_2,e_3$ 自身の順では $Ne_3$ が二項を含むが、この鎖を基底にすると $N$ は各ベクトルを左隣だけに送る。
:::

### 複数の鎖を同時に選ぶ計算

:::example[長い鎖の末端を二度使わない]
$Ne_1=0,\ Ne_2=e_1,\ Ne_3=e_2,\ Ne_4=e_1$ とする。$Nx=(x_2+x_4)e_1+x_3e_2,\ N^2x=x_3e_1,\ N^3x=0$ なので核の次元は $2,3,4$ である。増分 $2,1,1$ を読むと鎖の長さは $3,1$ である。

長さ $3$ の鎖に $(e_1,e_2,e_3)$ を取る。$e_4$ から始めると $Ne_4=e_1$ がすでに使った末端と重なる。そこで $N(e_4-e_2)=0$ に着目する。$\Ker N=\spn\{e_1,e_4-e_2\}$ なので $e_4-e_2$ を長さ $1$ の鎖に加えればよい。得られた四本は $e_4=(e_4-e_2)+e_2$ によって元の基底を復元でき、確かに基底である。
:::

## Jordan 標準形の存在と一意性

広義固有空間への分解と冪零作用素の鎖を組み合わせれば、特性多項式が分解する任意の作用素の標準形が得られる。

:::definition[Jordan ブロックと Jordan 標準形]{#def-linalg-jordan-block}
$\lambda\in K$ と $k\ge1$ に対し、対角成分がすべて $\lambda$、そのすぐ上の成分（$(i,i+1)$ 成分）がすべて $1$、他の成分が零の $k$ 次正方行列を $J_k(\lambda)$ と書き、:term[Jordan ブロック]{en="Jordan block" yomi="じょるだんぶろっく"} という。

$$
J_k(\lambda)=\begin{pmatrix}
\lambda&1&&\\
&\lambda&\ddots&\\
&&\ddots&1\\
&&&\lambda
\end{pmatrix}
$$

Jordan ブロックを対角に並べたブロック対角行列 $J_{k_1}(\lambda_1)\oplus\cdots\oplus J_{k_m}(\lambda_m)$ を :term[Jordan 標準形]{en="Jordan normal form" yomi="じょるだんひょうじゅんけい"}（Jordan 行列）という。
:::

$J_k(\lambda)=\lambda I+J_k(0)$ であり、$N=J_k(0)$ は標準基底について $Ne_1=0,\ Ne_i=e_{i-1}$（$i\ge2$）を満たす。
すなわち $(e_1,\ldots,e_k)$ は $e_k$ を頂点とする $N$ の長さ $k$ の鎖である。
$N^j$ は $(i,i+j)$ 成分だけが $1$ の行列であり、$\dim\Ker N^j=\min(j,k)$、$N^{k-1}\ne0,\ N^k=0$ である。
超対角に $1$ が来るのは、鎖の末端を最初の基底ベクトルに置いたためである（逆順に並べると下の副対角に $1$ が来る）。

:::theorem[Jordan 標準形]{#thm-linalg-jordan}
$\Phi_T$ が $K$ 上で一次式の積に分解するならば、ある基底に関する $T$ の表現行列は Jordan 標準形になる。
この Jordan 標準形はブロックの並べ方を除いて一意であり、固有値 $\lambda$ について $d_{\lambda,j}=\dim\Ker(T-\lambda I)^j$（$d_{\lambda,0}=0$）とおけば、大きさ $j$ のブロック $J_j(\lambda)$ の個数は

$$
2d_{\lambda,j}-d_{\lambda,j-1}-d_{\lambda,j+1}
$$

である。
逆に、ある基底で表現行列が Jordan 標準形になるならば、$\Phi_T$ は一次式の積に分解する。
行列で言えば、$\Phi_A$ が $K$ 上で分解する $A\in M_n(K)$ に対し、$P^{-1}AP$ が Jordan 標準形となる正則行列 $P\in M_n(K)$ が存在する。
特に複素正方行列はすべて $\C$ 上で Jordan 標準形に相似である。
:::

:::proof
:ref[cor-linalg-generalized-eigenspace-decomposition] により $V$ は広義固有空間の直和であり、各成分で $N=T-\lambda I$ は冪零である。:ref[thm-linalg-nilpotent-chains] の鎖基底を各成分に取る。長さ $k$ の鎖を末端から $w_1,\ldots,w_k$ と並べると $Tw_h=\lambda w_h+w_{h-1}$（$w_0=0$）であり、鎖の上で行列は $J_k(\lambda)$ である。:ref[prop-linalg-invariant-block] により全体の行列はこれらの直和となる。

逆にブロックの直和が与えられたとする。$\mu\ne\lambda$ のブロックでは $J_k(\mu)-\lambda I$ は対角が非零の三角行列なので核に寄与しない。$\mu=\lambda$ のブロックでは $j$ 乗の核は先頭の $\min(j,k)$ 本で張られる。よって $d_{\lambda,j}$ はこれらの $\min(j,k)$ の和である。差 $d_{\lambda,j}-d_{\lambda,j-1}$ は大きさ $j$ 以上の個数であり、もう一度差を取れば主張の式になる。核の次元は基底によらないためブロックの組も一意である。

また Jordan 標準形は三角行列なので :ref[prop-linalg-triangular-eigenvalues] より特性多項式は分解する。行列については鎖基底を列に並べれば $P$ を得る。複素数上では :ref[prop-linalg-irreducible-real-complex] より分解の仮定は自動的に成立する。
:::

:::corollary[Jordan 標準形から読める不変量]{#cor-linalg-jordan-invariants}
$\Phi_T$ が分解するとき、固有値 $\lambda$ について次が成り立つ。

1. $\lambda$ の Jordan ブロックの個数は幾何学的重複度 $g_\lambda=\dim E_\lambda(T)$ に等しい。
2. $\lambda$ の Jordan ブロックの大きさの総和は代数的重複度 $m_\lambda$ に等しい。
3. $\lambda$ の Jordan ブロックの大きさの最大値は、$m_T$ における $t-\lambda$ の指数 $s_\lambda$ に等しい。

特に $T$ が対角化可能であることと、すべての Jordan ブロックの大きさが $1$ であることは同値である。
:::

:::proof
1. 大きさ $1$ 以上のブロックの個数は $d_{\lambda,1}-d_{\lambda,0}=\dim\Ker(T-\lambda I)=g_\lambda$ である。
2. $\Phi_T=\prod_b(t-\mu_b)^{k_b}$ と :ref[lem-linalg-factor-theorem] の 4 による。
3. $J_k(\lambda)$ の最小多項式は $(t-\lambda)^k$ である（$J_k(0)^{k-1}\ne0=J_k(0)^k$ と :ref[thm-linalg-unique-factorization] の 3）。
異なる固有値のブロックの最小多項式は互いに素であり、:ref[prop-linalg-minpoly-lcm] の 2 により $m_T$ はそれらの最小公倍元 $\prod_\lambda(t-\lambda)^{\max k}$ である。
最後の主張は 3 と :ref[thm-linalg-minpoly-diagonalizable] による（あるいは 1, 2 と :ref[thm-linalg-diagonalization-criterion] による）。
:::

一つの固有値についても、ブロックの個数、大きさの総和、最大の大きさはそれぞれ別の量であり、混同しないことが大切である。
これら三つだけでは中間の大きさのブロックは決まらない。

:::example[同じ特性多項式と最小多項式]
$J_3(0)\oplus J_3(0)$ と $J_3(0)\oplus J_2(0)\oplus J_1(0)$ はいずれも特性多項式 $t^6$、最小多項式 $t^3$ をもつ。
しかし核の次元はそれぞれ $2,3$ で異なるから、相似ではない。
特性多項式と最小多項式だけでは、行列の相似類は一般に決まらない。
ただし $n\le3$ の場合は、各固有値のブロックの大きさの組が代数的重複度（$\le3$）と最大の大きさから決まるので、二つで足りる。
:::

:::example[基底変換を確かめる]
次の行列の固有値 $2$ に対して鎖を求める。

$$
T=\begin{pmatrix}2&1&1\\0&2&1\\0&0&2\end{pmatrix}
$$

$N=T-2I$ とおくと $Ne_3=e_1+e_2,\ N^2e_3=Ne_1+Ne_2=e_1$ である。
$S=(e_1,\ e_1+e_2,\ e_3)$（列として並べた行列）は正則で、$TS=SJ_3(2)$ を満たす。
行列を掛ける際は、基底を並べた $S$ が新しい座標を元の座標に戻す側であり、$S^{-1}TS=J_3(2)$ である。
:::

### 相似の判定

Jordan 標準形の一意性から、相似の完全な判定法が得られる。

:::corollary[階数による相似の判定]{#cor-linalg-jordan-similarity}
$A,B\in M_n(K)$ とし、$\Phi_A$ が $K$ 上分解するとする。
$A$ と $B$ が相似であるための必要十分条件は、$\Phi_A=\Phi_B$ であり、かつ $A$ の各固有値 $\lambda$ と $1\le j\le n$ について

$$
\rank(A-\lambda I)^j=\rank(B-\lambda I)^j
$$

が成り立つことである。
:::

:::proof
相似ならば $\Phi_A=\Phi_B$ であり（:ref[prop-linalg-charpoly-similarity]）、$P^{-1}(A-\lambda I)^jP=(B-\lambda I)^j$ なので階数も等しい。
逆を示す。
$\Phi_B=\Phi_A$ も分解するので、:ref[thm-linalg-jordan] により $A,B$ はそれぞれ Jordan 標準形 $J,J'$ に相似である。
両者の固有値は同じであり、:ref[thm-linalg-rank-nullity] により $d_{\lambda,j}=n-\rank(A-\lambda I)^j$ は $A$ と $B$ で一致する。
$j>n$ では :ref[prop-linalg-kernel-stabilization] の 3 により $d_{\lambda,j}=d_{\lambda,n}$ である。
ブロック数の公式から $J$ と $J'$ はブロックの並べ方を除いて一致し、ブロックの並べ替えは基底の並べ替え（置換行列による相似）で実現できる。
よって $A\sim J\sim J'\sim B$ である。
:::

:::corollary[転置との相似]{#cor-linalg-transpose-similar}
$\Phi_A$ が $K$ 上分解するならば、$A$ と $\transp{A}$ は相似である。
:::

:::proof
:ref[prop-linalg-det-ring] の 2（$K[t]$ 上で $\det\transp{M}=\det M$）より $\Phi_{\transp{A}}=\det(tI-\transp{A})=\det\transp{(tI-A)}=\Phi_A$ である。
また $(\transp{A}-\lambda I)^j=\transp{\bigl((A-\lambda I)^j\bigr)}$ であり、転置は階数を変えない（:ref[thm-linalg-row-column-rank]）。
:ref[cor-linalg-jordan-similarity] から主張が従う。
:::

具体的な相似行列は反転置換行列を用いて明示できる（演習問題）。

### Jordan 標準形を求める手順

定理の証明は、そのまま計算手順になる。
$A\in M_n(K)$ の $\Phi_A$ が分解するとき、次のように進める。

1. $\Phi_A$ を求めて因数分解し、相異なる固有値 $\lambda$ と代数的重複度 $m_\lambda$ を得る。
2. 各 $\lambda$ について $r_j=\rank(A-\lambda I)^j$ を $r_j=n-m_\lambda$ となるまで計算する（:ref[cor-linalg-generalized-eigenspace-decomposition] によりここで安定する）。$d_j=n-r_j$ からブロック数の公式で Jordan ブロックの大きさを決める。
3. 大きい鎖から順に、頂点 $v\in\Ker(A-\lambda I)^j$ を、その末端 $(A-\lambda I)^{j-1}v$ がそれまでに選んだ末端と一次独立になるように選ぶ。これは :ref[thm-linalg-nilpotent-chains] の証明で $L_j$ の基底を延長する操作である。
4. 各鎖を末端から頂点の順に並べ、鎖ごとに列として並べた行列を $P$ とすると $P^{-1}AP$ が Jordan 標準形になる。検算には逆行列を求めず $AP=PJ$ を確かめればよい。

:::example[二つの固有値をもつ四次行列]{#exm-linalg-jordan-4x4}
次の行列の Jordan 標準形と変換行列を求める。

$$
A=\begin{pmatrix}2&0&1&-1\\0&2&1&-1\\3&-3&5&-6\\3&-3&3&-4\end{pmatrix}
$$

**特性多項式**　:ref[prop-linalg-det-ring] を使う。$tI-A$ の第 2 行から第 1 行を、第 4 行から第 3 行を引くと、それぞれ $(t-2)(-1,1,0,0)$、$(t-2)(0,0,-1,1)$ になる。
この二つの因子 $t-2$ をくくり出し、さらに第 1 列に第 2 列を、第 3 列に第 4 列を加えると

$$
\Phi_A=(t-2)^2\det\begin{pmatrix}t-2&0&0&1\\0&1&0&0\\0&3&t+1&6\\0&0&0&1\end{pmatrix}=(t-2)^3(t+1)
$$

を得る（最後は第 1 列で展開する）。
固有値は $2$（代数的重複度 $3$）と $-1$（重複度 $1$）である。

**ブロックの大きさ**：$N=A-2I$ とおくと

$$
N=\begin{pmatrix}0&0&1&-1\\0&0&1&-1\\3&-3&3&-6\\3&-3&3&-6\end{pmatrix},\qquad
N^2=\begin{pmatrix}0&0&0&0\\0&0&0&0\\-9&9&-9&18\\-9&9&-9&18\end{pmatrix}
$$

であり、$\rank N=2$、$\rank N^2=1=4-3$ である。
よって $d_1=2,\ d_2=3$ で、ここで安定する。
大きさ $1$ 以上のブロックは $2$ 個、$2$ 以上は $1$ 個、$3$ 以上は $0$ 個なので、固有値 $2$ のブロックは $J_2(2)$ と $J_1(2)$ である。
固有値 $-1$ のブロックは $J_1(-1)$ である。
したがって

$$
J=J_2(2)\oplus J_1(2)\oplus J_1(-1),\qquad m_A=(t-2)^2(t+1)
$$

となる。

**鎖の選択**：$\Ker N$ は $x_3=x_4,\ x_1=x_2+x_4$ で定まり、$(1,1,0,0),(1,0,1,1)$ が基底である。
$\Ker N^2$ は $x_1-x_2+x_3-2x_4=0$ で定まる。
長さ $2$ の鎖の頂点として $\Ker N^2$ に属し $\Ker N$ に属さない $v_2=(-1,0,1,0)$ を選ぶと、末端は $v_1=Nv_2=(1,1,0,0)$ である。
次に長さ $1$ の鎖として、$\Ker N$ の元で $v_1$ と一次独立な $v_3=(1,0,1,1)$ を選ぶ。
ここで $v_3$ の代わりに $v_1$ の定数倍を選ぶと基底にならない。
最後に $A+I$ の核から $v_4=(0,0,1,1)$ を得る。

**検算**：$P=(v_1,v_2,v_3,v_4)$ とすると $\det P=1$ であり、

$$
P=\begin{pmatrix}1&-1&1&0\\1&0&0&0\\0&1&1&1\\0&0&1&1\end{pmatrix},\qquad
AP=PJ .
$$

例えば第 2 列は $Av_2=(-1,1,2,0)=2v_2+v_1$ である。
:::

### 実 Jordan 標準形

実行列の特性多項式は実数上で分解するとは限らない。
複素数上の Jordan 標準形を、共役な固有値を組にして実数の基底で書き直す。

:::theorem[実 Jordan 標準形]{#thm-linalg-real-jordan}
$A\in M_n(\R)$ とする。
$a,b\in\R$、$b>0$、$C=\begin{pmatrix}a&-b\\b&a\end{pmatrix}$ に対し、$2k$ 次の実行列

$$
R_k(a,b)=\begin{pmatrix}
C&I_2&&\\
&C&\ddots&\\
&&\ddots&I_2\\
&&&C
\end{pmatrix}
$$

を考える。
ある $P\in\GL_n(\R)$ について、$P^{-1}AP$ は、$A$ の実固有値 $\lambda$ に対する $J_k(\lambda)$ と、非実固有値の対 $a\pm ib$（$b>0$）に対する $R_k(a,b)$ のブロック対角行列になる。
$R_k(a,b)$ の個数は、$A$ の複素 Jordan 標準形における $J_k(a+ib)$ の個数に等しく、この実標準形はブロックの並べ方を除いて一意である。
:::

:::proof
実固有値 $\lambda$ については実空間 $\Ker_{\R}(A-\lambda I)^n$ 上で :ref[thm-linalg-nilpotent-chains] を用いる。実行列の核の複素解は実解の実部・虚部を組み合わせたものだから、ここで選んだ実鎖は複素広義固有空間の基底にもなる。

非実固有値 $\lambda=a+ib$（$b>0$）に対しては複素広義固有空間に鎖基底を取る。$A$ は実行列なので共役は $\lambda$ の鎖基底を $\bar\lambda$ の鎖基底に移す。一本の鎖を $(z_1,\ldots,z_k)$、$z_0=0$、$Az_h=\lambda z_h+z_{h-1}$ とする。$z_h=u_h-iw_h$ と書けば、$u_h=(z_h+\bar z_h)/2,\ w_h=(\bar z_h-z_h)/(2i)$ は実ベクトルである。この基底変換は各対で可逆であり、異なる広義固有空間は直和をなす。したがって実固有値の鎖と各対の $u_h,w_h$ を合わせた $n$ 本は複素基底をなし、実ベクトルなので実基底でもある。

鎖の式で実部と虚部を比較すると
$$
Au_h=au_h+bw_h+u_{h-1},\qquad Aw_h=-bu_h+aw_h+w_{h-1}.
$$
各組を $(u_h,w_h)$ の順に並べれば、対角の作用が $C$、一つ前の組への作用が $I_2$ であり、$R_k(a,b)$ を得る。逆にこのブロックの各座標対を $u_h\mp iw_h$ に戻せば $J_k(\lambda)$ と $J_k(\bar\lambda)$ が一つずつ現れる。したがって各実ブロックの個数は :ref[thm-linalg-jordan] の複素ブロック数から一意に復元できる。
:::

$C$ の形 $\begin{pmatrix}a&-b\\b&a\end{pmatrix}$ は、[正規作用素とスペクトル定理](/math/linear-algebra/spectral-theorem/) で実正規行列の標準形に用いた二次ブロックと同じ規約である。
$C$ は平面の $\sqrt{a^2+b^2}$ 倍の伸縮と回転の合成であり、固有値は $a\pm ib$ である。
基底を $(u_h,-w_h)$ の順に取れば、対角ブロックは $\begin{pmatrix}a&b\\-b&a\end{pmatrix}$ となる。
例えば $a=0,b=1$ のブロック $C$ は $90^\circ$ 回転であり、実数の範囲では二本の実直線に分解できない。

### 拡大体での最小多項式

最小多項式の次数 $d$ は、$I,A,\ldots,A^d$ の成分を縦に並べて初めて列が従属になる段である。$K$ の拡大体に移っても消去で使った非零のピボットはそのまま非零であり、階数は変わらない。したがって次数は同じである。もとのモニックな消去多項式は拡大体でも作用素を消すので、一意性から多項式そのものが一致する。特性多項式も :ref[prop-linalg-det-ring] の Leibniz の公式から変わらない。例えば実行列で変わるのは、既約な二次因子 $(t-a)^2+b^2$ が $\C$ 上で $(t-a-ib)(t-a+ib)$ に分かれる点である。

### 任意の体上の標準形について

$K$ 上で $\Phi_T$ が分解しない場合にも、相似類を完全に決める標準形がある。
:ref[lem-linalg-companion] の同伴行列を用いると、次が成り立つ。
任意の $T\in\End(V)$ に対し、モニック多項式 $q_1\mid q_2\mid\cdots\mid q_k$（次数 $1$ 以上）がただ一組存在して、ある基底に関する $T$ の表現行列が $C(q_1)\oplus\cdots\oplus C(q_k)$ となる。
このとき $q_k=m_T$、$q_1\cdots q_k=\Phi_T$ である。
これを :term[有理標準形]{en="rational canonical form" yomi="ゆうりひょうじゅんけい"}（Frobenius 標準形）といい、$q_i$ を :term[不変因子]{en="invariant factor" yomi="ふへんいんし"} という。
各 $q_i$ を既約分解して現れる冪 $f^e$ を :term[単因子]{en="elementary divisor" yomi="たんいんし"}（単因子は文献により不変因子を指すこともある）といい、二つの正方行列が相似であることと、不変因子（あるいは単因子）が一致することは同値である。
同値な条件は $tI-A$ と $tI-B$ が多項式の基本変形で移り合うことである。
$\Phi_T$ が分解する場合、単因子 $(t-\lambda)^k$ はちょうど Jordan ブロック $J_k(\lambda)$ に対応する。

$R=K[t]$ とおき、ベクトル多項式 $\sum_jt^jx_j$ を $\sum_jA^jx_j$ に送る。この写像は定数列だけですでに全射であり、$t$ 倍を $A$ の作用に移す。$(tI-A)R^n$ の元はすべて零に移る。逆に $t^jx-A^jx=(tI-A)\sum_{\ell=0}^{j-1}t^{j-1-\ell}A^\ell x$ なので、各ベクトル多項式はこの像を法としてその評価値に等しい。こうして $V$ は $R^n/(tI-A)R^n$ と同一視できる。

零でない多項式行列に可逆な行・列基本変形を施して、左上に非零な成分 $d$ を置く。同じ列に $d$ で割り切れない成分があれば、割り算で余りに変えて行を交換する。同じ行の成分には列操作を使う。どちらも新しい左上成分の次数が下がる。左上の成分がその行と列を割り切れば、他の成分をその行と列から消せる。残った右下ブロックに $d$ で割れない成分があれば、その成分を含む行を第 1 行に加え、再び列操作で次数を下げる。次数は負にならないので、有限回で残りの全成分を割る $d$ に達する。単元倍して $d$ をモニックにし、残りに同じ操作を施す。その後の操作も $d$ の倍数どうしの演算なので、順に各対角成分が次を割る。ここでは $\det(tI-A)\ne0$ なので零の対角成分は残らない。

行操作は分子の座標を、列操作は分母の生成元を替えるので、商は $\bigoplus_i R/(q_i)$ となる。対角成分 $1$ は零の商を与えるので省く。各 $R/(q_i)$ で $1,t,\ldots,t^{\deg q_i-1}$ は割り算の一意性から基底であり、$t$ 倍は $C(q_i)$ で表される。同伴行列の特性多項式は $q_i$ なので、全体ではその積になる。直和を消す多項式はすべての $q_i$ の公倍数であり、整除の順序から最小多項式は $q_k$ である。単因子から戻すには、既約因子ごとに指数を昇順に並べ、最長の列に合わせて左を零で埋め、位置ごとに冪を掛ける。これは $q_i\mid q_{i+1}$ が各指数を非減少にするからである。分解体では、各同伴ブロックの最小多項式と特性多項式が等しい。従ってその一次分解で $(t-\lambda)^e$ に属する空間は次元 $e$、最大鎖長も $e$ となり、$J_e(\lambda)$ 一個を与える。

一意性には、$j$ 次小行列式全部のモニック最大公約元を用いる。行の定数倍（非零の体の元）、行交換、他の行の多項式倍の加算によって、変形後の各小行列式は変形前のものの $R$ 線形結合になる。これは :ref[prop-linalg-det-ring] の多重線形性と交代性による。列操作も同じで、逆操作も可能なので最大公約元は変わらない。対角形では先頭 $j$ 個の積がすべての $j$ 次小行列式を割り、自身も小行列式である。したがって隣り合う次数の最大公約元の商から、$1$ も含めた対角成分がただ一通りに戻る。逆に各 $tI-C(q_i)$ には下の副対角から作る次数 $\deg q_i-1$ の単元小行列式があるので、最後以外の対角成分は $1$ となる。最後の成分は行列式 $q_i$ である。したがって候補の同伴ブロックを並べても同じ手順でその $q_i$ を回収できる。また、相似変換 $A\mapsto P^{-1}AP$ は $tI-A$ の可逆な定数行・列操作である。逆に同じ対角形をもてば同じ同伴ブロックの直和に相似である。体を広げてもこの対角形は有効で、転置でも小行列式は変わらない。

## Jordan 標準形の応用

対角化できない行列でも、各ブロックに分ければ、スカラー部分と冪零部分の二項展開で計算できる。
行列の指数関数など冪級数で定まる関数への応用は [行列解析](/math/linear-algebra/matrix-analysis/) で扱う。

:::proposition[ブロックへの多項式の代入]{#prop-linalg-jordan-polynomial}
$p\in K[t]$ について $p(\lambda+u)=\sum_jc_ju^j$ と展開し、$N=J_s(0)$ とおくと

$$
p(J_s(\lambda))=\sum_{j=0}^{s-1}c_jN^j
$$

となる。
特に非負整数 $m$ について

$$
J_s(\lambda)^m=\sum_{j=0}^{\min(m,s-1)}\binom mj\lambda^{m-j}N^j
$$

である（$\lambda^0=1$ と読む）。
$K$ の標数が零ならば $c_j=p^{(j)}(\lambda)/j!$ である。
:::

:::proof
$p(\lambda+u)=\sum_jc_ju^j$ は $K[u]$ における多項式の等式であり、:ref[lem-linalg-polynomial-substitution] の 1 により両辺に可換な行列を代入してよい。
$\lambda I$ と $N$ は可換なので $u=N$ を代入すると $p(\lambda I+N)=\sum_jc_jN^j$ を得て、$N^s=0$ により $j\ge s$ の項は消える。
$p(t)=t^m$ なら二項定理が係数 $c_j=\binom mj\lambda^{m-j}$ を与える。
標数零では、展開の両辺を $u$ で $j$ 回微分して $u=0$ とおけば $c_j=p^{(j)}(\lambda)/j!$ である。
:::

$J_s(\lambda)^m$ では同じ超対角上の成分は等しい。第 $j$ 超対角の成分は $j\le m$ なら $\binom mj\lambda^{m-j}$、$j>m$ なら零である。

:::example[反復による増加]
$J_3(1)^m$ の対角は $1$、第一超対角は $m$、第二超対角は $\binom m2$ である。
したがって初期値 $e_3$ は三回の反復で $(3,3,1)$ に移る。
固有値の絶対値がすべて $1$ でも、反復が有界とは限らない。
:::

### Jordan–Chevalley 分解

:::theorem[Jordan–Chevalley 分解]{#thm-linalg-jordan-chevalley}
$\Phi_T$ が $K$ 上分解するとする。
対角化可能な $D$ と冪零な $N$ で

$$
T=D+N,\qquad DN=ND
$$

を満たすものがただ一組存在する。
さらに $D,N$ はともに $T$ の多項式である。
:::

:::proof
:ref[thm-linalg-primary-decomposition] の射影 $P_i$ を用いて $D=\sum_i\lambda_iP_i$、$N=T-D$ とおく。両者は $T$ の多項式なので可換である。各広義固有空間上では $D=\lambda_iI$、$N=T-\lambda_iI$ なので、前者は全体で対角化可能、後者は $N^n=0$ を満たす。

別の分解 $T=D'+N'$ でも $D',N'$ は和 $T$ と可換であるから、その多項式 $D,N$ とも可換である。そこで $D-D'=N'-N$ の両辺を比べる。左辺は :ref[thm-linalg-simultaneous-diagonalization] により対角化可能である。右辺の $2n$ 乗を二項展開すると、各項に $N^n$ か $(N')^n$ が含まれて零になる。冪零かつ対角化可能な作用素は固有基底のすべてを零に送るから零である。よって $D=D',N=N'$ となる。
:::

$D$ を $T$ の :term[半単純部分]{en="semisimple part" yomi="はんたんじゅんぶぶん"}、$N$ を :term[冪零部分]{en="nilpotent part" yomi="べきれいぶぶん"} という。
Jordan 標準形の基底では、$D$ は対角部分、$N$ は超対角の $1$ たちである。
しかし任意の基底で「対角部分 + 狭義上三角部分」と分けても、一般にはこの分解にならない（演習問題）。
可換性の条件が本質的である。

この分解から冪の公式が得られる。
$D$ と $N$ は可換なので二項定理が使え、$N^s=0$（$s$ は冪零指数）として

$$
T^m=\sum_{j=0}^{\min(m,s-1)}\binom mjD^{m-j}N^j,\qquad D^{m-j}=\sum_i\lambda_i^{m-j}P_i
$$

である。

:::example[四次行列の冪の公式]
:ref[exm-linalg-jordan-4x4] の $A$ では $m_A=(t-2)^2(t+1)$ である。
$S=A-2I$ は固有値 $2$ の広義固有空間上で $S^2=0$、固有値 $-1$ の固有空間上で $-3$ 倍なので、$P_{-1}=S^2/9$、$P_2=I-P_{-1}$ である。
したがって

$$
P_{-1}=\begin{pmatrix}0&0&0&0\\0&0&0&0\\-1&1&-1&2\\-1&1&-1&2\end{pmatrix},\qquad
N=A-2P_2+P_{-1}=\begin{pmatrix}0&0&1&-1\\0&0&1&-1\\0&0&0&0\\0&0&0&0\end{pmatrix}
$$

であり、$N^2=O$、$NP_{-1}=O$ である。
よって $m\ge1$ について

$$
A^m=2^mP_2+(-1)^mP_{-1}+m\,2^{m-1}N
$$

を得る。
例えば $m=1$ では右辺は $D+N=A$ に戻る。
:::

### 線形漸化式

Jordan ブロックが現れる典型例は、特性方程式が重根をもつ線形漸化式である。
以下 $K=\R$ または $\C$ とする。

:::theorem[定数係数線形漸化式の解]{#thm-linalg-linear-recurrence}
$c(t)=t^d+a_{d-1}t^{d-1}+\cdots+a_1t+a_0$（$a_0\ne0$）が $K$ 上で $c(t)=\prod_{i=1}^r(t-\lambda_i)^{e_i}$（$\lambda_i$ は相異なる）と分解するとする。
数列 $(x_m)_{m\ge0}$ で

$$
x_{m+d}+a_{d-1}x_{m+d-1}+\cdots+a_0x_m=0\qquad(m\ge0)
$$

を満たすものの全体 $X$ は $d$ 次元の部分空間であり、$d$ 個の数列 $\bigl(m^k\lambda_i^m\bigr)_{m\ge0}$（$1\le i\le r$、$0\le k<e_i$）がその基底をなす。
:::

:::proof
数列の全体 $K^{\N}$ はベクトル空間であり、その上の :term[シフト作用素]{en="shift operator" yomi="しふとさようそ"} $S$ を $(Sx)_m=x_{m+1}$ で定める。
漸化式は $c(S)x=0$ と書けるので $X=\Ker c(S)$ は部分空間である。
初期値 $(x_0,\ldots,x_{d-1})$ を任意に与えると漸化式により数列がただ一つ定まるので、$x\mapsto(x_0,\ldots,x_{d-1})$ は $X$ から $K^d$ への同型写像であり、$\dim X=d$ である。

:ref[thm-linalg-unique-factorization] の 3 により $(t-\lambda_i)^{e_i}$ たちは対ごとに互いに素であり、さらに一つと残りの積も互いに素である。
:ref[lem-linalg-coprime-splitting]（無限次元でも成り立つ）を繰り返し用いて

$$
X=\Ker c(S)=\bigoplus_{i=1}^r\Ker(S-\lambda_i)^{e_i}
$$

を得る。

次に $\lambda\ne0$ と多項式 $q$ に対し $x_m=q(m)\lambda^m$ とおくと

$$
\bigl((S-\lambda)x\bigr)_m=q(m+1)\lambda^{m+1}-\lambda q(m)\lambda^m=\lambda^{m+1}(\Delta q)(m),\qquad(\Delta q)(u)=q(u+1)-q(u)
$$

である。
$\Delta q$ は $q$ より次数が一つ低い（$q$ が定数なら零）ので、$\deg q<e$ ならば $(S-\lambda)^ex=0$ である。
よって $U_i=\spn\{(m^k\lambda_i^m)_m: 0\le k<e_i\}$ は $\Ker(S-\lambda_i)^{e_i}$ に含まれる。
$\lambda_i\ne0$（$a_0\ne0$ による）なので、$\sum_kc_km^k\lambda_i^m=0$ がすべての $m$ で成り立てば、多項式 $\sum_kc_ku^k$ は無限個の点 $u=0,1,2,\ldots$ で零になり、:ref[lem-linalg-factor-theorem] の 2 により零多項式である。
したがって $\dim U_i=e_i$ である。
$U_i$ は直和をなす部分空間 $\Ker(S-\lambda_i)^{e_i}$ に含まれるので $\sum_iU_i$ も直和であり、その次元は $\sum_ie_i=d=\dim X$ である。
よって $X=\bigoplus_iU_i$ であり、主張の $d$ 個の数列は $X$ の基底である。
:::

$X$ は $S$ 不変であり、:ref[prop-linalg-charpoly-minpoly-factors] の観点から言えば、$S|_X$ の最小多項式と特性多項式はともに $c(t)$ である（初期値の座標で $S|_X$ を表すと $C(c)$ の転置になる）。
各固有値に Jordan ブロックが一つずつあり、大きさ $e_i$ のブロックが数列 $m^k\lambda_i^m$（$k<e_i$）を生む。

:::example[重根をもつ漸化式]
次の初期条件をもつ二階の漸化式を解く。

$$
a_{m+2}=4a_{m+1}-4a_m,\qquad a_0=1,\ a_1=4
$$

$c(t)=t^2-4t+4=(t-2)^2$ なので、定理により $a_m=(\alpha+\beta m)2^m$ と書ける。
$a_0=\alpha=1$、$a_1=2(1+\beta)=4$ から $\beta=1$ であり、$a_m=(m+1)2^m$ である。
直接にも、$b_m=a_{m+1}-2a_m$ とおけば $b_{m+1}=2b_m$、$b_0=2$ なので $a_{m+1}/2^{m+1}-a_m/2^m=1$ となり、同じ答えを得る。
$(a_m,b_m)$ を $(a_{m+1},b_{m+1})$ に移す行列は $\begin{pmatrix}2&1\\0&2\end{pmatrix}=J_2(2)$ であり、因子 $m+1$ は Jordan ブロックの超対角から生じている。
:::

### 巡回ベクトルと可換な作用素

:::proposition[最小多項式を実現するベクトルの存在]{#prop-linalg-cyclic-vector}
$m_{T,v}=m_T$ を満たす $v\in V$ が存在する。
:::

:::proof
$m_T=\prod_if_i^{s_i}$ と既約分解し、:ref[thm-linalg-primary-decomposition] の成分 $V_i$ を考える。
$V_i$ の基底 $(w_1,\ldots,w_l)$ について、各 $m_{T,w_j}$ は $f_i^{s_i}$ を割り切るので $f_i^{k_j}$ の形であり（:ref[thm-linalg-unique-factorization] の 3）、:ref[prop-linalg-minpoly-lcm] の 1 を $T|_{V_i}$ に適用すると $f_i^{s_i}=m_{T|_{V_i}}=f_i^{\max_jk_j}$ である。
よって $m_{T,v_i}=f_i^{s_i}$ となる $v_i\in V_i$ がある。
$v=\sum_iv_i$ とおく。
$p(T)v=0$ ならば、$p(T)v_i\in V_i$ と直和性からすべての $i$ で $p(T)v_i=0$、すなわち $f_i^{s_i}\mid p$ である。
:ref[thm-linalg-unique-factorization] の 4 により $m_T\mid p$ であり、$m_{T,v}=m_T$ を得る。
:::

:::corollary[巡回ベクトルをもつ条件]{#cor-linalg-cyclic-criterion}
次は同値である。

1. $T$ は巡回ベクトルをもつ。
2. $m_T=\Phi_T$ である。
3. ある基底に関する $T$ の表現行列は $C(\Phi_T)$ である。

$\Phi_T$ が分解するときは、さらに次とも同値である。

4. 各固有値の Jordan ブロックはただ一つである（すなわち $g_\lambda=1$）。
:::

:::proof
1⇒3：巡回ベクトル $v$ について $Z(v)=V$ なので、:ref[prop-linalg-cyclic-subspace] の 3 により表現行列は $C(m_{T,v})$ であり、その特性多項式 $m_{T,v}$ は $\Phi_T$ に等しい。
3⇒2 は :ref[lem-linalg-companion] による。
2⇒1：:ref[prop-linalg-cyclic-vector] の $v$ について $\deg m_{T,v}=\deg\Phi_T=n$ なので $\dim Z(v)=n$、$Z(v)=V$ である。
2⇔4：:ref[cor-linalg-jordan-invariants] により $\Phi_T=\prod_\lambda(t-\lambda)^{(\text{大きさの和})}$、$m_T=\prod_\lambda(t-\lambda)^{(\text{大きさの最大値})}$ なので、両者が等しいことはブロックが固有値ごとに一つであることと同値である。
:::

:::proposition[巡回ベクトルをもつ作用素と可換な作用素]{#prop-linalg-cyclic-commutant}
$T$ が巡回ベクトルをもつならば、$S\in\End(V)$ が $T$ と可換であるための必要十分条件は、$S=p(T)$ となる $p\in K[t]$ が存在することである。
$\deg p<n$ とすれば $p$ は一意である。
:::

:::proof
多項式 $p(T)$ が $T$ と可換なのは明らかである。
逆に $ST=TS$ とし、$v$ を巡回ベクトルとする。
$(v,Tv,\ldots,T^{n-1}v)$ は基底なので $Sv=\sum_{k<n}c_kT^kv=p(T)v$ と書ける。
すると $S(T^kv)=T^kSv=T^kp(T)v=p(T)(T^kv)$ であり、$S$ と $p(T)$ は基底上で一致するので等しい。
$I,T,\ldots,T^{n-1}$ は $\deg m_T=n$ より一次独立なので、$p$ は一意である。
:::

巡回ベクトルがない場合には、可換な作用素はもっと多い。
冪零な場合を Jordan 基底で調べる。

:::proposition[二つの鎖の間で作用を保つ写像]{#prop-linalg-chain-intertwiner}
$N_s=J_s(0)$、$N_r=J_r(0)$ とする。
$X\in M_{r\times s}(K)$ のうち $XN_s=N_rX$ を満たすもの全体は $\min(r,s)$ 次元の部分空間である。
:::

:::proof
定義域の鎖の頂点を $v=e_s$ とすると、$N_s^jv=e_{s-j}$ である。
$XN_s=N_rX$ ならば、$y=Xv$ を指定すれば残りの基底の像は $X(N_s^jv)=N_r^jy$（$j=0,\ldots,s-1$）で決まり、$N_r^sy=XN_s^sv=0$ が必要である。
逆に $N_r^sy=0$ を満たす $y$ を任意に取り、$X(e_{s-j})=N_r^jy$ で $X$ を定めると、$j\le s-2$ では $XN_se_{s-j}=Xe_{s-j-1}=N_r^{j+1}y=N_rXe_{s-j}$、$j=s-1$ では $XN_se_1=0=N_r^sy=N_rXe_1$ なので $XN_s=N_rX$ である。
したがって求める空間は $\Ker N_r^s$ と同型であり、その次元は $\min(r,s)$ である。
:::

:::example[可換性を保つ長方形行列]
$s=3,\ r=2$ とする。
頂点の像を $y=(a,b)$ とすれば $N_2y=(b,0),\ N_2^2y=0$ である。
定義域の標準基底は $(N_3^2v,N_3v,v)=(e_1,e_2,e_3)$ なので、求める行列は

$$
X=\begin{pmatrix}0&b&a\\0&0&b\end{pmatrix}
$$

である。
実際、$XN_3$ と $N_2X$ はともに $(1,3)$ 成分だけが $b$ の行列である。
:::

$N$ が大きさ $s_1,\ldots,s_m$ の冪零 Jordan ブロックの直和のとき、$X$ を対応する長方形ブロック $X_{ij}$ に分けると、$XN=NX$ は各 $X_{ij}$ に $X_{ij}N_{s_j}=N_{s_i}X_{ij}$ という独立な条件を課す。
したがって $N$ と可換な作用素全体の次元は

$$
\sum_{i=1}^m\sum_{j=1}^m\min(s_i,s_j)
$$

である。
例えば $J_3(0)\oplus J_1(0)$ では $3+1+1+1=6$ である。
特性多項式が $K$ 上で一次式の積に分解する $T$ では、$T$ と可換な $X$ は $T$ の多項式である射影 $P_i$ と可換なので各広義固有空間を保ち、その上で $T$ との可換性は $T-\lambda_iI$ との可換性に等しい。
したがって次元は固有値ごとの上の和の総和である。
単一のブロック $J_n(\lambda)$ では次元は $n$ であり、:ref[prop-linalg-cyclic-commutant] と整合する。

### 多項式を作用させたあとのブロック

$J_s(\lambda)$ に多項式 $p$ を代入すると、固有値は $p(\lambda)$ に移る（:ref[cor-linalg-spectral-mapping]）が、ブロックの大きさも変わりうる。

:::proposition[多項式によるブロックの分裂]{#prop-linalg-polynomial-block-splitting}
$J=J_s(\lambda)$、$p\in K[t]$ を定数でない多項式とし、:ref[lem-linalg-factor-theorem] の 3 により $p(t)-p(\lambda)=(t-\lambda)^hg(t)$、$g(\lambda)\ne0$ と書く。
$s=qh+r$（$0\le r<h$）とすると、$p(J)$ の Jordan 標準形は、大きさ $q+1$ のブロック $J_{q+1}(p(\lambda))$ が $r$ 個と、大きさ $q$ のブロック $J_q(p(\lambda))$ が $h-r$ 個（$q=0$ なら無し）の直和である。
特に $h=1$（任意の標数で $p'(\lambda)\ne0$ と同値）ならば $p(J)\sim J_s(p(\lambda))$ である。
:::

:::proof
$N=J-\lambda I$ とおく。
$g(J)$ は対角成分が $g(\lambda)\ne0$ の上三角行列なので正則で、$N$ と可換である。
よって $(p(J)-p(\lambda)I)^j=N^{hj}g(J)^j$ であり、$\dim\Ker(p(J)-p(\lambda)I)^j=\dim\Ker N^{hj}=\min(hj,s)$ である。
:ref[thm-linalg-jordan] のブロック数の公式によれば、大きさ $j$ 以上のブロックの個数は $\min(hj,s)-\min(h(j-1),s)$ であり、これは $j\le q$ で $h$、$j=q+1$ で $r$、$j\ge q+2$ で $0$ である。
したがって大きさ $q+1$ のものが $r$ 個、ちょうど $q$ のものが $h-r$ 個である。
形式微分による言い換えは、$p'(t)=h(t-\lambda)^{h-1}g(t)+(t-\lambda)^hg'(t)$ から $p'(\lambda)\ne0\iff h=1$ による。
:::

:::example[鎖の途中を飛ぶとどうなるか]
$J=J_5(0)$、$p(t)=t^2$ とすると $h=2$、$5=2\cdot2+1$ なので、$J^2\sim J_3(0)\oplus J_2(0)$ である。
実際、標準基底で $J^2e_j=e_{j-2}$（$j\le2$ なら零）なので、二つの鎖は $(e_1,e_3,e_5)$ と $(e_2,e_4)$ に分かれる。
:::

### 実 Jordan ブロックに対応する振動

:::example[共役固有値と二段の鎖]
$C=\begin{pmatrix}0&-1\\1&0\end{pmatrix}$ とし、四次実行列を

$$
A=\begin{pmatrix}C&I_2\\O&C\end{pmatrix}=R_2(0,1)
$$

とする。
$C^2=-I_2$ なので $A^2+I_4$ の対角ブロックは零、右上は $2C$ である。
したがって $A^2+I_4\ne O$ だが $(A^2+I_4)^2=O$ となる。
:ref[thm-linalg-unique-factorization] の 3 により $(t^2+1)^2$ のモニックな真の約元は $1,\ t^2+1$ だけであり、どちらでも $A$ は消えないので、実最小多項式は $(t^2+1)^2$ である。

複素数上では $C(1,-i)=i(1,-i)$ である。
$w=(1,-i)$ として $(w,0),(0,w)$ を並べると、$A(w,0)=i(w,0)$、$A(0,w)=(w,0)+i(0,w)$ なので、$A$ はこの二本の上で $J_2(i)$ として作用する。
共役をとった二本の上では $J_2(-i)$ となる。
:ref[thm-linalg-real-jordan] のとおり、一つの実ブロック $R_2(0,1)$ に長さ $2$ の複素鎖が二本含まれている。
:::

### 数値計算での注意

正確な係数の上では、核の次元列が相似類を区別する。
一方、計算した固有値が近いという理由だけで同じ Jordan ブロックにまとめることはできない。
例えば

$$
A_h=\begin{pmatrix}1&1\\0&1+h\end{pmatrix}
$$

は $h\ne0$ なら固有値が異なるので対角化できるが、$h=0$ では $J_2(1)$ になる。
$h\ne0$ のときの固有ベクトルを $(1,0),(1,h)$ と選ぶと、$h\to0$ で両者はほぼ平行になり、二本を並べた行列の行列式は $h$ である。
このように Jordan 標準形は行列の成分に対して連続ではなく、固有基底の条件の悪さが摂動に対する固有値の感度を左右する。
その定量的な評価は [行列解析](/math/linear-algebra/matrix-analysis/) で扱う。

## 演習問題

::::exercise{level=basic}
$A=\diag(1,1,-1)$ の特性多項式と最小多項式を求めよ。

:::solution
特性多項式は対角成分に対応する因子の積 $(t-1)^2(t+1)$ である。
$(A-I)(A+I)$ は対角成分が $0\cdot2,\ 0\cdot2,\ (-2)\cdot0$ の対角行列なので零であり、$A$ はスカラー行列でないので一次式では消えない。
:ref[prop-linalg-charpoly-minpoly-factors] から最小多項式の候補は $(t-1)(t+1)$ と $(t-1)^2(t+1)$ に限られ、前者が消去多項式なので $m_A=(t-1)(t+1)$ である。
:::
::::

::::exercise{level=basic}
$p=t^4-1,\ q=t^3+t^2-t-1$ の最大公約元を互除法で求め、$\gcd(p,q)=ap+bq$ となる $a,b$ を一組求めよ。

:::solution
$p$ を $q$ で割ると $p=(t-1)q+(2t^2-2)$ である（$(t-1)q=t^4-2t^2+1$）。
次に $q=\bigl(\tfrac12t+\tfrac12\bigr)(2t^2-2)$ と余りなく割り切れる。
よって最後の零でない余り $2t^2-2$ をモニックにした $t^2-1$ が最大公約元である。
第一式から

$$
t^2-1=\tfrac12(2t^2-2)=\tfrac12\,p-\tfrac12(t-1)\,q
$$

なので、$a=\tfrac12,\ b=-\tfrac12(t-1)$ が一組の解である。
実際 $p=(t-1)(t+1)(t^2+1)$、$q=(t-1)(t+1)^2$ と分解すれば、共通因子は $(t-1)(t+1)$ である。
:::
::::

::::exercise{level=basic}
$J_3(2)^4$ を求めよ。

:::solution
$N=J_3(0)$ とすると、:ref[prop-linalg-jordan-polynomial] により $(2I+N)^4=16I+\binom41 8N+\binom42 4N^2=16I+32N+24N^2$ である（$N^3=O$）。
したがって対角、第一超対角、第二超対角が順に $16,32,24$、それ以外が零の上三角行列である。
:::
::::

::::exercise{level=basic}
六次元空間の冪零作用素 $N$ について、$\Ker N,\Ker N^2,\Ker N^3$ の次元が $3,5,6$ であるとする。
$N$ の Jordan 標準形を求めよ。

:::solution
増分 $3,2,1$ は、長さが $1,2,3$ 以上の鎖の本数である（:ref[thm-linalg-nilpotent-chains]）。
したがって鎖は三本で、そのうち二本が長さ $2$ 以上、一本が長さ $3$ である。
よって長さは $3,2,1$ であり、標準形は $J_3(0)\oplus J_2(0)\oplus J_1(0)$ である。
:::
::::

::::exercise{level=basic}
$A=J_2(1)\oplus J_1(3)$ に対し、固有値 $3$ の広義固有空間への一次分解の射影を $A$ の多項式として求めよ。

:::solution
$m_A=(t-1)^2(t-3)$ である。
$e(t)=(t-1)^2/4$ は $(t-1)^2$ で割り切れ、$e(t)-1=(t-3)(t+1)/4$ は $t-3$ で割り切れるので、:ref[thm-linalg-primary-decomposition] の証明の $\eta_i$ の条件を満たす。
よって射影は $(A-I)^2/4$ であり、$J_2(0)^2=O$ と $(3-1)^2/4=1$ から $\diag(0,0,1)$ となる。
これは任意の列の第三成分だけを残す。
:::
::::

::::exercise{level=basic}
次の行列の Jordan 標準形 $J$ と、$P^{-1}BP=J$ となる正則行列 $P$ を求めよ。

$$
B=\begin{pmatrix}1&1&0\\-1&3&0\\1&1&2\end{pmatrix}
$$

:::solution
第 3 列で展開すると $\Phi_B=(t-2)\bigl((t-1)(t-3)+1\bigr)=(t-2)^3$ である。
$N=B-2I$ とおくと

$$
N=\begin{pmatrix}-1&1&0\\-1&1&0\\1&1&0\end{pmatrix},\qquad
N^2=\begin{pmatrix}0&0&0\\0&0&0\\-2&2&0\end{pmatrix},\qquad N^3=O
$$

なので、$\rank N=2,\ \rank N^2=1$ であり、$d_1,d_2,d_3=1,2,3$ である。
鎖は長さ $3$ の一本で、$J=J_3(2)$、$m_B=(t-2)^3$ である。
頂点として $N^2e_1\ne0$ を満たす $e_1$ を選ぶと、$Ne_1=(-1,-1,1)$、$N^2e_1=(0,0,-2)$ である。
よって

$$
P=\begin{pmatrix}0&-1&1\\0&-1&0\\-2&1&0\end{pmatrix}
$$

とすれば $\det P=-2\ne0$ であり、$BP=PJ_3(2)$ が成り立つ。
例えば第 2 列は $B(-1,-1,1)=(-2,-2,0)=2(-1,-1,1)+(0,0,-2)$ である。
:::
::::

::::exercise{level=basic}
$a_{m+3}=3a_{m+2}-3a_{m+1}+a_m$、$a_0=0,\ a_1=1,\ a_2=4$ を満たす数列の一般項を求めよ。

:::solution
特性多項式は $t^3-3t^2+3t-1=(t-1)^3$ なので、:ref[thm-linalg-linear-recurrence] により $a_m=\alpha+\beta m+\gamma m^2$ と書ける。
$a_0=\alpha=0$、$a_1=\beta+\gamma=1$、$a_2=2\beta+4\gamma=4$ から $\gamma=1,\ \beta=0$ である。
よって $a_m=m^2$ であり、実際 $9=3\cdot4-3\cdot1+0$ などが成り立つ。
:::
::::

::::exercise{level=basic}
$W$ を $k$ 次の反転置換行列（$(i,k+1-i)$ 成分が $1$、他は零）とする。
$W^2=I$ かつ $WJ_k(\lambda)W=\transp{J_k(\lambda)}$ を示し、これを用いて Jordan 標準形 $J$ と $P^{-1}AP=J$ から $\transp{A}$ と $A$ の相似を与える行列を明示せよ。

:::solution
$We_j=e_{k+1-j}$ なので $W^2e_j=e_j$ である。
$WMW$ の $(i,j)$ 成分は $M$ の $(k+1-i,k+1-j)$ 成分である。
$M=J_k(\lambda)$ のとき、これは $i=j$ で $\lambda$、$k+1-j=(k+1-i)+1$ すなわち $j=i-1$ で $1$、他は零なので、$WJ_k(\lambda)W=\transp{J_k(\lambda)}$ である。
$J=\bigoplus_bJ_{k_b}(\lambda_b)$ に対し $\widetilde W=\bigoplus_bW_{k_b}$ とすれば $\widetilde W^2=I$、$\widetilde WJ\widetilde W=\transp{J}$ である。
$A=PJP^{-1}$ から $\transp{A}=\transp{(P^{-1})}\transp{J}\transp{P}=\transp{(P^{-1})}\widetilde WJ\widetilde W\transp{P}$ であり、$Q=P\widetilde W\transp{P}$ とおくと $Q^{-1}=\transp{(P^{-1})}\widetilde WP^{-1}$ なので

$$
Q^{-1}AQ=\transp{(P^{-1})}\widetilde W\,J\,\widetilde W\transp{P}=\transp{A}
$$

である。
:::
::::

::::exercise{level=advanced}
$\Phi_T$ が分解するとき、固有値 $\lambda$ の大きさ $j$ の Jordan ブロックの個数を $r_j=\rank(T-\lambda I)^j$（$r_0=n$）だけで表せ。
これを用いて、$T$ が対角化可能であることと、すべての固有値 $\lambda$ について $\rank(T-\lambda I)=\rank(T-\lambda I)^2$ であることが同値であることを示せ。

:::solution
:ref[thm-linalg-rank-nullity] から $d_j=n-r_j$ なので、:ref[thm-linalg-jordan] の公式は

$$
2(n-r_j)-(n-r_{j-1})-(n-r_{j+1})=r_{j-1}-2r_j+r_{j+1}
$$

となる。
後半：$r_1=r_2$ ならば $d_1=d_2$ であり、:ref[prop-linalg-kernel-stabilization] の 2 から $d_j=d_1$（$j\ge1$）となる。
大きさ $2$ 以上のブロックの個数は $d_2-d_1=0$ である。
これがすべての固有値で成り立てば、すべてのブロックの大きさは $1$ であり、:ref[cor-linalg-jordan-invariants] により対角化可能である。
逆に対角化可能ならブロックの大きさはすべて $1$ なので $d_2-d_1=0$、すなわち $r_1=r_2$ である。
:::
::::

::::exercise{level=advanced}
$A\in M_n(\C)$ が $A^k=I$（$k\ge1$）を満たすならば、$A$ は対角化可能であることを示せ。
また $\R$ 上では対角化可能とは限らない例と、標数 $2$ の体上で対角化できない例を挙げよ。

:::solution
$\omega=e^{2\pi i/k}$ とすると $1,\omega,\ldots,\omega^{k-1}$ は相異なる $k$ 個の $t^k-1$ の根である。
:ref[lem-linalg-factor-theorem] の 1 を繰り返し用いると、次数とモニック性から $t^k-1=\prod_{j=0}^{k-1}(t-\omega^j)$ である。
$A$ はこの相異なる一次式の積で消えるので、:ref[thm-linalg-minpoly-diagonalizable] の (c) により $\C$ 上対角化可能である。
$\R$ 上では、$90^\circ$ 回転 $\begin{pmatrix}0&-1\\1&0\end{pmatrix}$ は $A^4=I$ を満たすが、実固有値をもたないので実数上では対角化できない。
$\F_2$ 上では $\begin{pmatrix}1&1\\0&1\end{pmatrix}^2=\begin{pmatrix}1&2\\0&1\end{pmatrix}=I$ だが、これは $J_2(1)$ であり対角化できない。
標数 $2$ では $t^2-1=(t-1)^2$ が重根をもつことが原因である。
:::
::::

::::exercise{level=advanced}
複素行列 $T$ について、$m\to\infty$ で $T^m\to O$（各成分の収束）となることと、すべての固有値が $\abs\lambda<1$ を満たすことが同値であることを示せ。

:::solution
必要性：固有ベクトル $v$ について $T^mv=\lambda^mv$ であり、$T^m\to O$ ならば $\lambda^mv\to0$、したがって $\abs\lambda<1$ である。
十分性：$T=PJP^{-1}$ とすると $T^m=PJ^mP^{-1}$ であり、$P,P^{-1}$ は固定されているので $J^m\to O$ を示せばよい。
:ref[prop-linalg-jordan-polynomial] により、大きさ $s$ のブロック $J_s(\lambda)^m$ の成分は $\binom mj\lambda^{m-j}$（$0\le j<s$）である。
$\lambda=0$ ならば $m\ge s$ で零である。
$0<\abs\lambda<1$ ならば $\binom mj\le m^j$ より成分の絶対値は $m^j\abs\lambda^{m-j}$ 以下である。
$c_m=m^j\abs\lambda^m$ とおくと $c_{m+1}/c_m=(1+1/m)^j\abs\lambda\to\abs\lambda<1$ なので、ある $r<1$ と $M$ について $m\ge M$ で $c_{m+1}\le rc_m$ となり、$c_m\le c_Mr^{m-M}\to0$ である。
よって $J^m$ の全成分は零に収束する。
:::
::::

::::exercise{level=advanced}
六次の複素行列 $A$ の特性多項式が $(t-2)^4(t+1)^2$ であり、$\rank(A-2I)=4,\ \rank(A-2I)^2=3,\ \rank(A+I)=5$ とする。
Jordan 標準形と最小多項式を求めよ。
また $\rank(A-2I)$ だけでは Jordan 標準形が決まらない理由を述べよ。

:::solution
固有値 $2$ について $d_1=6-4=2,\ d_2=6-3=3$ である。
ブロックは $d_1=2$ 個で、大きさ $2$ 以上のものは $d_2-d_1=1$ 個である。
大きさの総和は代数的重複度 $4$ なので、一方が $1$、他方が $3$ である（$d_3=4$ と整合する）。
固有値 $-1$ については $d_1=6-5=1$ なので、ブロックは一つで大きさは $2$ である。
よって

$$
A\sim J_3(2)\oplus J_1(2)\oplus J_2(-1),\qquad m_A=(t-2)^3(t+1)^2
$$

である（:ref[cor-linalg-jordan-invariants] の 3）。
固有値 $2$ のブロックが二つであることだけなら、大きさは $3+1$ と $2+2$ の両方がありうる。
後者なら $d_2=4$、すなわち $\rank(A-2I)^2=2$ となるので、二乗の階数を調べて初めて区別できる。
:::
::::

::::exercise{level=advanced}
$J_7(0)^3$ の Jordan 標準形を、:ref[prop-linalg-polynomial-block-splitting] と核の次元の直接計算の両方で求めよ。

:::solution
$p(t)=t^3$、$\lambda=0$ では $h=3$ であり、$7=2\cdot3+1$ なので、大きさ $3$ のブロックが $1$ 個、大きさ $2$ のブロックが $3-1=2$ 個である。
よって $J_7(0)^3\sim J_3(0)\oplus J_2(0)\oplus J_2(0)$ である。
直接には、$M=J_7(0)^3$ について $\dim\Ker M^j=\dim\Ker J_7(0)^{3j}=\min(3j,7)$ なので、$d_1,d_2,d_3=3,6,7$ である。
増分 $3,3,1$ から、長さ $1,2,3$ 以上の鎖がそれぞれ $3,3,1$ 本であり、長さは $3,2,2$ となる。
鎖は $(e_1,e_4,e_7)$、$(e_2,e_5)$、$(e_3,e_6)$ である。
:::
::::

::::exercise{level=advanced}
$T=\begin{pmatrix}1&1&1\\0&1&0\\0&0&3\end{pmatrix}$ の Jordan–Chevalley 分解 $T=D+N$ を求め、$D,N$ を $T$ の多項式として表せ。
また $T$ を対角部分 $\diag(1,1,3)$ と狭義上三角部分に分けたものはこの分解でないことを確かめよ。

:::solution
$m_T=(t-1)^2(t-3)$ であり、本文の例で見たように固有値 $3$ への射影は $P_3=(T-I)^2/4$、$P_1=I-P_3$ である。
よって $D=P_1+3P_3=I+2P_3=I+\tfrac12(T-I)^2$ である。
$(T-I)^2=\begin{pmatrix}0&0&2\\0&0&0\\0&0&4\end{pmatrix}$ なので

$$
D=\begin{pmatrix}1&0&1\\0&1&0\\0&0&3\end{pmatrix},\qquad
N=T-D=\begin{pmatrix}0&1&0\\0&0&0\\0&0&0\end{pmatrix}
$$

である。
$N=T-I-\tfrac12(T-I)^2$ も $T$ の多項式であり、$N^2=O$、$DN=ND$ が成り立つ。
$D$ は最小多項式 $(t-1)(t-3)$ で消えるので対角化可能である。
一方、$D_0=\diag(1,1,3)$、$N_0=T-D_0=\begin{pmatrix}0&1&1\\0&0&0\\0&0&0\end{pmatrix}$ とすると、$D_0N_0-N_0D_0$ の $(1,3)$ 成分は $1-3=-2$ で零でない。
可換でないので、一意性からこれは Jordan–Chevalley 分解ではない。
:::
::::

## 確認テスト

::::quiz
$m_T$ が $K$ 上で相異なる一次因子の積であるとき、常に成り立つものを選べ。

- [x] $T$ は $K$ 上で対角化可能である
- [ ] すべての固有空間が一次元である
- [ ] $m_T=\Phi_T$ である

:::explanation
:ref[thm-linalg-minpoly-diagonalizable] により、一次分解の各成分が固有空間になるので固有基底が存在する。
例えば $I_2$ の最小多項式は $t-1$ だが、固有空間は二次元で、$\Phi=(t-1)^2\ne m$ である。
したがって残りの二つは一般には成り立たない。
:::
::::

::::quiz{answer="2" tol="1e-6"}
$J_4(0)\oplus J_2(0)\oplus J_1(0)$ で長さ $2$ 以上の鎖は何本あるか。

:::explanation
長さ $4$ と $2$ の二本が該当する。
核の次元で計算しても、$d_1=3$、$d_2=2+2+1=5$ より $d_2-d_1=2$ である。
:::
::::

::::quiz
Jordan 標準形から一意に決まるものをすべて選べ。

- [x] 最小多項式
- [ ] 元の座標での変換行列
- [x] 各固有値の幾何学的重複度

:::explanation
:ref[cor-linalg-jordan-invariants] により、各固有値の最大のブロックの大きさが最小多項式の指数を決め、ブロックの個数が固有空間の次元を決める。
一方、鎖を零でない定数倍しても標準形は変わらないので、変換行列は一意ではない。
:::
::::

::::quiz{answer="6" tol="1e-6"}
$J_3(1)^4$ の $(1,3)$ 成分を求めよ。

:::explanation
$J_3(1)^4=\sum_j\binom4jN^j$（$N=J_3(0)$）であり、$(1,3)$ 成分に寄与するのは $N^2$ の項だけなので、値は $\binom42=6$ である。
:::
::::

::::quiz{answer="5" tol="1e-6"}
$4$ 次の冪零複素行列の相似類はいくつあるか。

:::explanation
冪零行列の Jordan 標準形は固有値 $0$ のブロックだけからなり、:ref[thm-linalg-jordan] によりブロックの大きさの組（和が $4$ の分割）で相似類が決まる。
分割は $4,\ 3+1,\ 2+2,\ 2+1+1,\ 1+1+1+1$ の $5$ 通りである。
:::
::::

::::quiz{answer="5" tol="1e-6"}
$J_2(0)\oplus J_1(0)$ と可換な $3$ 次正方行列全体の空間の次元を求めよ。

:::explanation
ブロックの大きさ $2,1$ について $\sum_{i,j}\min(s_i,s_j)=\min(2,2)+\min(2,1)+\min(1,2)+\min(1,1)=2+1+1+1=5$ である。
直接にも、$XN=NX$ を成分で書くと独立な条件が $4$ 個あり、$9-4=5$ となる。
:::
::::

---
title: 直積・半直積と有限 Abel 群
description: 直積と半直積による群の組み立てと分解を調べて小さい位数の群の分類に用い、Smith 標準形を通して有限生成 Abel 群の基本定理を不変量の一意性まで含めて証明する。
notation:
  - symbol: 'H\times K,\ A\oplus B'
    meaning: 群の直積と Abel 群の直和
  - symbol: 'N\rtimes_\alpha H'
    meaning: 準同型 α による半直積
  - symbol: '\Aut(N)'
    meaning: 群 N の自己同型群
  - symbol: '(\Z/n\Z)^\times'
    meaning: 法 n の単数群（乗法群）
  - symbol: 'T(A)'
    meaning: Abel 群 A の捩れ部分群
  - symbol: 'A_p'
    meaning: 有限 Abel 群 A の p 準素成分
  - symbol: 'A[m]'
    meaning: m 倍すると零になる元全体
  - symbol: '\GL_n(\Z)'
    meaning: 整数成分で逆行列も整数成分の n 次正方行列全体
---

## この章の目標

- 外部直積と内部直積を対応させ、群が直積に分解する条件を判定する。
- 半直積を構成し、分裂する完全列・補群・自己同型群の計算と結びつけて、位数 $18,20,30$ の群を同型を除いて分類する。
- 整数行列の Smith 標準形を証明し、有限生成 Abel 群の基本定理（存在と一意性）を導く。
- 不変因子と単因子を相互に変換し、有限 Abel 群の元・部分群・生成系を数え、体の乗法群の有限部分群が巡回群であることを示す。

これまでの章では、与えられた群の部分群や剰余群を調べてきた。本章では逆に、小さな群を組み合わせて大きな群を作る二つの方法、直積と半直積を扱う。直積は二つの群を「互いに干渉せずに」並べる構成であり、半直積は一方の群がもう一方に自己同型として作用しながら並ぶ構成である。どちらも、群の中に適切な部分群の組を見つけることで内部的に認識できる。この認識法と前章の Sylow の定理を組み合わせると、位数の小さい群を系統的に分類できる。

後半では Abel 群に限定し、有限個の元で生成される Abel 群がすべて巡回群の直和に分解され、その分解の型が群だけで決まることを示す。証明の中心は、整数行列を整数係数の基本変形で対角化する Smith 標準形の理論である。これは線形代数の基本変形を、割り算のできない整数環 $\Z$ の上で行うものであり、生成元と関係式で与えられた Abel 群を具体的に計算する手段にもなる。同じ定理を単項イデアル整域上の加群に一般化したものは[単項イデアル整域上の加群](/math/rings-modules/pid-modules/)で扱う。

## 直積

### 外部直積

二つの独立な時計の針を一緒に動かすと、両方の針が同時に出発点へ戻るのは、それぞれの周期の最小公倍数の時刻である。一方の針の動きはもう一方に影響しない。このような「独立に並べた操作」を記述するのが直積である。

:::definition[直積]{#def-grp-direct-product}
群 $H_1,\ldots,H_n$ に対し、集合としての直積 $H_1\times\cdots\times H_n$ に成分ごとの積

$$
(h_1,\ldots,h_n)(h_1',\ldots,h_n')=(h_1h_1',\ldots,h_nh_n')
$$

を入れたものを $H_1,\ldots,H_n$ の :term[直積]{en="direct product" yomi="ちょくせき"} といい、同じ記号 $H_1\times\cdots\times H_n$ で表す。各 $H_i$ が加法で書かれた Abel 群のときは、同じ構成を :term[直和]{en="direct sum" yomi="ちょくわ"} とも呼び $H_1\oplus\cdots\oplus H_n$ と書く。$H$ の $n$ 個の直積を $H^n$ と書く。
:::

結合法則は各成分で $H_i$ の結合法則から従う。単位元は $(e,\ldots,e)$、$(h_1,\ldots,h_n)$ の逆元は $(h_1^{-1},\ldots,h_n^{-1})$ である。従って直積は群をなし、その位数は $|H_1|\cdots|H_n|$ である。$H_1\times H_2\times H_3$ と $(H_1\times H_2)\times H_3$ は $(h_1,h_2,h_3)\leftrightarrow((h_1,h_2),h_3)$ で同型であり、因子の順序の入れ替えも同型を与える。

任意の添字集合 $I$ についても成分ごとの演算で直積 $\prod_{i\in I}H_i$ を作れる。Abel 群の場合、非零成分が有限個だけの元からなる部分群が直和 $\bigoplus_{i\in I}H_i$ である。例えば $\prod_{i\ge1}C_2$ の全成分が $\bar1$ の元は直和に入らない。非自明な因子が有限個だけなら、添字集合が無限でも両者は一致する。

:::proposition[直積の普遍性]{#prop-grp-product-universal}
$P=H_1\times\cdots\times H_n$ とし、$\pi_i:P\to H_i$ を第 $i$ 成分への射影、$\iota_i:H_i\to P$ を第 $i$ 成分に $h$ を置き他の成分を $e$ とする写像とする。

1. $\pi_i$ は全射準同型、$\iota_i$ は単射準同型であり、$\pi_i\iota_i=\id_{H_i}$ である。
2. 任意の群 $B$ と準同型 $f_i:B\to H_i$（$1\le i\le n$）に対し、$\pi_i\circ f=f_i$ をすべての $i$ で満たす準同型 $f:B\to P$ がただ一つ存在する。それは $f(b)=(f_1(b),\ldots,f_n(b))$ である。
3. $\iota_i(H_i)$ は $P$ の正規部分群であり、$i\ne j$ なら $\iota_i(H_i)$ の元と $\iota_j(H_j)$ の元は可換である。
:::

:::proof
射影と埋め込みは成分ごとの積を保つ。$\pi_i\iota_i=\id$ より前者は全射、後者は単射である。$f(b)$ の各成分を指定すると組 $f(b)=(f_i(b))_i$ が一意に決まり、各 $f_i$ が積を保つのでこの $f$ も準同型である。$\iota_i(H_i)$ は残る全成分への射影の核だから正規である。異なる座標に埋め込んだ二元は、各成分で少なくとも一方が単位元なので可換である。
:::

2 は「$B$ から直積への準同型を与えることは、$B$ から各因子への準同型を一つずつ与えることと同じである」と言い換えられる。例えば $\Z\to H\times K$ という準同型は $1$ の像 $(h,k)$ で決まり、これは $\Z\to H$ と $\Z\to K$ の組 $m\mapsto h^m$、$m\mapsto k^m$ に対応する。

:::proposition[直積の元の位数]{#prop-grp-product-order}
$(h_1,\ldots,h_n)\in H_1\times\cdots\times H_n$ が有限位数であるための必要十分条件は、各 $h_i$ が有限位数であることであり、そのとき

$$
\operatorname{ord}(h_1,\ldots,h_n)=\operatorname{lcm}\bigl(\operatorname{ord}(h_1),\ldots,\operatorname{ord}(h_n)\bigr)
$$

である。
:::

:::proof
$(h_1,\ldots,h_n)^k=(h_1^k,\ldots,h_n^k)$ であるから、この元が単位元になるのは、すべての $i$ で $h_i^k=e$ となるときである。どれかの $h_i$ が無限位数なら、正の $k$ で $h_i^k=e$ とはならないので組も無限位数である。すべて有限位数なら、:ref[thm-grp-cyclic-order] により $h_i^k=e$ は $\operatorname{ord}(h_i)\mid k$ と同値である。従って $(h_1,\ldots,h_n)^k=e$ は $k$ がすべての $\operatorname{ord}(h_i)$ の公倍数であることと同値であり、そのような最小の正整数は最小公倍数である。
:::

例えば $C_4\times C_6$ の元 $(\bar1,\bar1)$ の位数は $\operatorname{lcm}(4,6)=12$ である。どの元の位数も $4$ の約数と $6$ の約数の最小公倍数なので $12$ を超えない。位数 $24$ の元がないから、$C_4\times C_6$ は巡回群ではない。同様に $C_2\times C_2$ の非単位元はすべて位数 $2$ なので、$C_2\times C_2\not\cong C_4$ である。

### 巡回群の直積と中国剰余定理

:::proposition[巡回群の直積]{#prop-grp-cyclic-product}
$m,n\ge1$ とする。$C_m\times C_n$ が巡回群であるための必要十分条件は $\gcd(m,n)=1$ である。このとき

$$
\rho:\Z/mn\Z\to\Z/m\Z\times\Z/n\Z,\qquad x\bmod mn\mapsto(x\bmod m,\ x\bmod n)
$$

は群の同型であり、$um+vn=1$ を満たす整数 $u,v$ を用いると、逆写像は $(a\bmod m,\ b\bmod n)\mapsto(avn+bum)\bmod mn$ で与えられる。
:::

:::proof
直積の全元の位数は $\operatorname{lcm}(m,n)$ を割るので、$\gcd(m,n)>1$ なら位数 $mn$ の元はない。両者が互いに素なら還元 $\rho$ は和を保ち、その核は $m,n$ の双方で割り切れる類、すなわち $\bar0$ だけである。両辺の位数は $mn$ だから $\rho$ は同型である。$um+vn=1$ より $avn+bum$ は法 $m$ で $a$、法 $n$ で $b$ に等しく、指定の逆写像を得る。
:::

:::example[中国剰余定理の計算]{#exm-grp-crt}
$x\equiv2\pmod3$、$x\equiv1\pmod4$ を満たす $x$ を法 $12$ で求める。$(-1)\cdot3+1\cdot4=1$ なので $u=-1,\ v=1$ とおき、$x=2\cdot1\cdot4+1\cdot(-1)\cdot3=5$ を得る。実際 $5\equiv2\pmod3$、$5\equiv1\pmod4$ である。これは「互いに素な法についての連立合同式は、積の法でただ一つの解をもつ」という中国剰余定理の、群の同型 $C_{12}\cong C_3\times C_4$ による言い換えである。

:ref[prop-grp-cyclic-product] を繰り返し用いると、$n=p_1^{a_1}\cdots p_k^{a_k}$（$p_i$ は相異なる素数）に対して

$$
C_n\cong C_{p_1^{a_1}}\times\cdots\times C_{p_k^{a_k}}
$$

となる。実際、$p_1^{a_1}$ と残りの積は互いに素なので一因子ずつ切り離せる。例えば $C_{60}\cong C_4\times C_3\times C_5$ である。一方 $C_{12}\not\cong C_2\times C_6$ である。後者の元の位数は $6$ 以下だからである。巡回群の直積がどう分解するかは、素因数ごとに見なければ判断できない。
:::

還元 $\rho$ は成分ごとに乗法も保つので、単数群にも適用できる。

:::proposition[単数群の分解]{#prop-grp-units-crt}
$\gcd(m,n)=1$ なら、$\rho$ は群の同型

$$
(\Z/mn\Z)^\times\cong(\Z/m\Z)^\times\times(\Z/n\Z)^\times
$$

を誘導する。従って Euler の関数について $\varphi(mn)=\varphi(m)\varphi(n)$ であり、$n=\prod_ip_i^{a_i}$ と素因数分解すると

$$
\varphi(n)=\prod_ip_i^{a_i-1}(p_i-1). \label{eq-grp-euler-phi}
$$

:::

:::proof
$\rho$ は積と $1$ を保つ全単射なので、$xy=1$ は $\rho(x)\rho(y)=(1,1)$ と同値である。従って単数だけに制限しても全単射である。位数を比べれば $\varphi$ の乗法性が従う。法 $p^a$ の非単元は $p$ の倍数の $p^{a-1}$ 類なので、$\varphi(p^a)=p^{a-1}(p-1)$ を掛け合わせればよい。
:::

例えば $(\Z/15\Z)^\times\cong(\Z/3\Z)^\times\times(\Z/5\Z)^\times\cong C_2\times C_4$ である。実際、法 $15$ で $2$ の冪は $2,4,8,1$ となり位数 $4$ をもつが、位数 $8$ の元はない。また $(\Z/12\Z)^\times\cong(\Z/4\Z)^\times\times(\Z/3\Z)^\times\cong C_2\times C_2$ であり、全元の平方が $1$ になることと一致する。$\varphi(60)=\varphi(4)\varphi(3)\varphi(5)=2\cdot2\cdot4=16$ である。

### 内部直積

群が与えられたとき、それが直積に分解するかどうかを部分群の言葉で判定したい。:ref[prop-grp-product-universal] の 3 により、直積 $H\times K$ の中では $\iota_1(H)$ と $\iota_2(K)$ が正規で、互いに可換で、自明に交わり、全体を生成する。この性質が直積を特徴づける。

:::proposition[内部直積の判定]{#prop-grp-internal-direct}
$H,K\le G$ とする。積写像

$$
\mu:H\times K\to G,\qquad (h,k)\mapsto hk
$$

が群の同型であるための必要十分条件は、$H\trianglelefteq G$、$K\trianglelefteq G$、$H\cap K=\{e\}$、$HK=G$ が成り立つことである。このとき $G$ は $H$ と $K$ の :term[内部直積]{en="internal direct product" yomi="ないぶちょくせき"} であるという。
:::

:::proof
$\mu$ が同型であるとする。$\mu$ は $\iota_1(H)=H\times\{e\}$ を $H$ に、$\iota_2(K)$ を $K$ に移すので、:ref[prop-grp-product-universal] の 3 から $H,K$ は正規である。$\mu$ の全射性は $HK=G$ を意味する。$x\in H\cap K$ なら $\mu(x,e)=x=\mu(e,x)$ で、$\mu$ の単射性から $(x,e)=(e,x)$、すなわち $x=e$ である。

逆に四条件を仮定する。まず $h\in H,\ k\in K$ が可換であることを示す。交換子 $c=hkh^{-1}k^{-1}$ は、$(hkh^{-1})k^{-1}$ と見れば $K$ の正規性から $K$ に属し、$h(kh^{-1}k^{-1})$ と見れば $H$ の正規性から $H$ に属する。よって $c\in H\cap K=\{e\}$ であり $hk=kh$ である。すると

$$
\mu\bigl((h,k)(h',k')\bigr)=hh'kk'=hkh'k'=\mu(h,k)\,\mu(h',k')
$$

となり、$\mu$ は準同型である。$HK=G$ から全射である。$\mu(h,k)=e$ なら $h=k^{-1}\in H\cap K$ なので $h=k=e$ となり、核が自明だから単射である。
:::

:::remark
二つの部分群の正規性はどちらも外せない。$S_3$ では $N=\langle(1\,2\,3)\rangle$ と $H=\langle(1\,2)\rangle$ が $N\cap H=\{e\}$、$|NH|=6$ を満たし、$N$ は正規である。しかし $S_3$ は非可換なので可換群 $C_3\times C_2$ と同型でない。このように一方だけが正規な場合が、次節の半直積である。
:::

:::proposition[有限個の因子の内部直積]{#prop-grp-internal-direct-n}
$N_1,\ldots,N_n\trianglelefteq G$ が $G=N_1N_2\cdots N_n$ を満たし、各 $i$ について

$$
N_i\cap\bigl(N_1\cdots N_{i-1}N_{i+1}\cdots N_n\bigr)=\{e\}
$$

であるとする。このとき $(x_1,\ldots,x_n)\mapsto x_1\cdots x_n$ は $N_1\times\cdots\times N_n$ から $G$ への同型である。
:::

:::proof
$i\ne j$ なら $N_i\cap N_j=1$ である。$x\in N_i,\ y\in N_j$ の交換子は両群の正規性からこの共通部分に入るので $xy=yx$ となる。従って積写像は準同型であり、仮定より全射である。$x_1\cdots x_n=e$ なら各 $x_i$ は他の因子の元の積の逆元だから、仮定の共通部分に入って $e$ となる。よって核は自明である。
:::

条件を「どの二つも自明に交わる」に弱めることはできない。$G=C_2\times C_2$ の三つの位数 $2$ の部分群はどの二つも自明に交わり、三つを合わせると全体を生成する。しかし三つの直積は位数 $8$ であり、位数 $4$ の $G$ と同型ではない。実際、第三の部分群は最初の二つの積 $G$ に含まれ、上の条件が破れている。

:::example[内部直積の例]{#exm-grp-internal-direct}
1. $C_{30}$ で $\bar{15},\bar{10},\bar6$ の位数は $2,3,5$ で、和は $\bar1$ である。どの一成分の位数も残る二成分の位数の積と互いに素なので、共通部分の元は単位元に限られる。よって $C_{30}\cong C_2\times C_3\times C_5$。
2. $x=\sgn(x)|x|$ の一意性と積公式より $\R^\times=\{\pm1\}\times\R_{>0}$ と内部直積に分かれる。同じ議論を $z=(z/|z|)|z|$ に適用すれば $\C^\times$ の単位円と正の実数への分解を得る。
3. 奇数 $m\ge3$ について $D_{2m}\cong D_m\times C_2$ である。$H=\langle r^2,s\rangle$ は標準形より位数 $2m$ で $D_m$ と同じ積をもち、指数 $2$ なので正規である。半回転 $r^m$ は中心に入る。$m$ が奇数なので $r^m\notin H$ であり、二群の積は $4m$ 個の元をもつ。
4. $D_4$ の非自明な分解があれば因子の位数は $2,4$ で、どちらも可換群のため矛盾する。$Q_8$ の非自明な部分群はすべて $-1$ を含む（位数 $4$ の元の平方は $-1$）。従ってこの二群は直積に分解できない。
:::

### 直積の部分群と中心

:::proposition[直積の中心と剰余群]{#prop-grp-product-center}
1. $Z(H\times K)=Z(H)\times Z(K)$ である。特に $H\times K$ が可換であることと、$H,K$ がともに可換であることは同値である。
2. $N\trianglelefteq H,\ M\trianglelefteq K$ なら $N\times M\trianglelefteq H\times K$ であり、$(H\times K)/(N\times M)\cong H/N\times K/M$ である。
:::

:::proof
1. $(z,w)$ がすべての $(h,k)$ と可換であることは、$zh=hz$ と $wk=kw$ がすべての $h,k$ で成り立つことと同値である（一方の成分を $e$ にとればよい）。

2. $(h,k)\mapsto(hN,kM)$ は :ref[prop-grp-product-universal] の 2 により準同型 $H\times K\to H/N\times K/M$ であり、明らかに全射である。核は $h\in N$ かつ $k\in M$ となる組の全体 $N\times M$ なので正規であり、:ref[thm-grp-first-iso] から同型を得る。
:::

直積の部分群がすべて「部分群の直積」の形をしているわけではない。$H$ の **対角部分群** $\Delta_H=\{(h,h):h\in H\}\le H\times H$ はその典型である。$\Delta_H$ が正規であるのは $H$ が可換のときに限る。実際、$(a,e)(h,h)(a,e)^{-1}=(aha^{-1},h)$ が $\Delta_H$ に入る条件は $aha^{-1}=h$ である。

$H=K=C_p$（$p$ は素数）の場合を数えると、$C_p\times C_p$ の位数 $p$ の部分群は、$p^2-1$ 個の非零元を、各部分群が $p-1$ 個ずつ含むので $(p^2-1)/(p-1)=p+1$ 個ある。自明部分群と全体を合わせて部分群は $p+3$ 個であるが、そのうち部分群の直積の形のものは $\{0\}\times\{0\}$、$C_p\times\{0\}$、$\{0\}\times C_p$、全体の $4$ 個だけである。

## 半直積

$S_3$ や二面体群 $D_n$ は、回転の部分群 $N=\langle r\rangle$ と一つの鏡映の生成する部分群 $H=\langle s\rangle$ によって $G=NH$、$N\cap H=\{e\}$ と分解する。しかし $H$ は正規でなく、$rs\ne sr$ である。その代わり、鏡映を回転の右から左へ移すと、$sr^k=r^{-k}s$ のように回転の向きが反転する。つまり $H$ の元は共役によって $N$ に自己同型として作用し、積はこの作用で決まる。この構造を抽象化したものが半直積である。

### 外部半直積

:::definition[半直積]{#def-grp-semidirect}
群 $N,H$ と準同型 $\alpha:H\to\Aut(N)$、$h\mapsto\alpha_h$ が与えられたとする。集合 $N\times H$ に

$$
(n,h)(m,k)=\bigl(n\,\alpha_h(m),\ hk\bigr) \label{eq-grp-semidirect-product}
$$

で積を入れたものを、$\alpha$ による $N$ と $H$ の :term[半直積]{en="semidirect product" yomi="はんちょくせき"} といい、$N\rtimes_\alpha H$ と書く。作用 $\alpha$ が文脈から明らかなときは $N\rtimes H$ と略記する。
:::

記号 $\rtimes$ の三角形の開いた側の $N$ が、次の命題で見るように正規部分群になる因子である。

:::proposition[半直積の基本性質]{#prop-grp-semidirect-basic}
$G=N\rtimes_\alpha H$ とする。

1. $G$ は群であり、単位元は $(e,e)$、$(n,h)$ の逆元は $(\alpha_{h^{-1}}(n^{-1}),\ h^{-1})$ である。
2. $\hat N=\{(n,e)\}$、$\hat H=\{(e,h)\}$ はそれぞれ $N,H$ と同型な部分群であり、$\hat N$ は射影 $(n,h)\mapsto h$ の核として正規である。さらに $\hat N\cap\hat H=\{e\}$、$G=\hat N\hat H$ で、$(n,h)=(n,e)(e,h)$ と書ける。
3. $(e,h)(n,e)(e,h)^{-1}=(\alpha_h(n),e)$ である。
4. $\hat H$ が正規であることと $\alpha$ が自明（すべての $h$ で $\alpha_h=\id_N$）であることは同値であり、そのとき $N\rtimes_\alpha H=N\times H$ である。
:::

:::proof
1. 三元 $(n,h),(m,k),(l,j)$ をどちらに括っても第一成分は $n\alpha_h(m)\alpha_{hk}(l)$、第二成分は $hkj$ になる。ここで $\alpha_h$ が積を保つことと $\alpha_h\alpha_k=\alpha_{hk}$ を使った。$\alpha_e=\id$ より $(e,e)$ は単位元である。$(n,h)$ と $(\alpha_{h^{-1}}(n^{-1}),h^{-1})$ の積はどちらの順でも第一成分が $e$、第二成分も $e$ となる。
2. 各座標への埋め込みは積の式から単射準同型である。第二成分への射影の核が $\hat N$ なのでこれは正規であり、共通部分は $(e,e)$ のみである。$(n,e)(e,h)=(n,h)$ より全体が得られる。
3. $(e,h)(n,e)(e,h^{-1})=(\alpha_h(n),e)$ である。
4. $(n,e)(e,h)(n^{-1},e)=(n\alpha_h(n)^{-1},h)$ が $\hat H$ に入る条件は $\alpha_h(n)=n$ である。従って $\hat H$ が正規なら作用は自明で、逆に自明なら積は直積の式となり両因子が正規である。
:::

3 は、半直積の中では $H$ の元による共役が $N$ 上で指定された自己同型 $\alpha_h$ として実現されることを示している。逆に、群の中で一方が正規な分解が見つかれば、群は半直積として再構成できる。

:::theorem[内部半直積]{#thm-grp-semidirect-internal}
$N\trianglelefteq G,\ H\le G$ が $N\cap H=\{e\}$、$NH=G$ を満たすとする。$\alpha_h(n)=hnh^{-1}$ とおくと $\alpha:H\to\Aut(N)$ は準同型であり、

$$
\Phi:N\rtimes_\alpha H\to G,\qquad (n,h)\mapsto nh
$$

は同型である。このとき $G$ は $N$ と $H$ の **内部半直積** であるといい、$H$ を $N$ の :term[補群]{en="complement" yomi="ほぐん"} という。逆に任意の $N\rtimes_\alpha H$ は、:ref[prop-grp-semidirect-basic] の $\hat N$ と $\hat H$ の内部半直積である。
:::

:::proof
$N$ が正規なので $\alpha_h(N)\subseteq N$ であり、$\alpha_h$ は積を保ち、$\alpha_{h^{-1}}$ を逆写像にもつから $\alpha_h\in\Aut(N)$ である。$\alpha_h\alpha_k(n)=hknk^{-1}h^{-1}=\alpha_{hk}(n)$ なので $\alpha$ は準同型である。

$\Phi$ は $G=NH$ より全射である。$nh=n'h'$ なら $n'^{-1}n=h'h^{-1}\in N\cap H=\{e\}$ なので $n=n'$、$h=h'$ となり単射である（:ref[prop-grp-internal-semidirect]）。さらに

$$
\Phi\bigl((n,h)(m,k)\bigr)=n\,(hmh^{-1})\,hk=nh\,mk=\Phi(n,h)\,\Phi(m,k)
$$

だから準同型である。後半は :ref[prop-grp-semidirect-basic] の 2, 3 と、前半を $\hat N,\hat H$ に適用した結果から従う。
:::

:::example[二面体群とアフィン群]{#exm-grp-semidirect-basic}
1. $D_n$ で $N=\langle r\rangle$、$H=\langle s\rangle$ とすると、$N$ は指数 $2$ で正規、$N\cap H=\{e\}$、$|NH|=2n$ である。$srs^{-1}=r^{-1}$ なので、$s$ は $N\cong C_n$ に反転 $x\mapsto-x$ として作用し、$D_n\cong C_n\rtimes_{-1}C_2$ である。ここで $C_n\rtimes_{-1}C_2$ は $C_2$ の生成元が $-1$ 倍で作用する半直積を表す。
2. 素数 $p$ について、$\F_p$ 上のアフィン変換 $x\mapsto ax+b$（$a\in\F_p^\times,\ b\in\F_p$）全体は合成について群 $\mathrm{Aff}(\F_p)$ をなす。$(x\mapsto ax+b)\circ(x\mapsto cx+d)$ は $x\mapsto acx+(ad+b)$ なので、変換を $(b,a)$ で表すと積は $(b,a)(d,c)=(b+ad,\ ac)$ となる。これは :eqref[eq-grp-semidirect-product] で $N=\F_p$（加法群）、$H=\F_p^\times$、$\alpha_a(d)=ad$ としたものであり、$\mathrm{Aff}(\F_p)\cong\F_p\rtimes\F_p^\times$、位数 $p(p-1)$ である。平行移動全体が正規部分群 $\hat N$ に、原点を固定する変換全体が補群 $\hat H$ に対応する。倍率 $a$ の変換で平行移動 $x\mapsto x+d$ を共役すると $x\mapsto x+ad$ になる。$p=3$ では位数 $6$ の非可換群なので、$\mathrm{Aff}(\F_3)\cong S_3$ である（$\F_3$ の三点の置換として忠実に作用し、位数が等しい）。
3. 任意の群 $N$ に対し、恒等写像 $\Aut(N)\to\Aut(N)$ による半直積 $N\rtimes\Aut(N)$ を $N$ の **holomorph** といい $\operatorname{Hol}(N)$ と書く。2 の群は $\operatorname{Hol}(C_p)$ にほかならない（次の :ref[prop-grp-aut-cyclic] を参照）。
:::

### 分裂する完全列

半直積を、全射準同型の側から見直す。

:::definition[短完全列と分裂]{#def-grp-split-sequence}
群と準同型の列 $1\to N\xrightarrow{\iota}G\xrightarrow{\pi}H\to1$ で、$\iota$ が単射、$\pi$ が全射、$\Img\iota=\Ker\pi$ を満たすものを :term[短完全列]{en="short exact sequence" yomi="たんかんぜんれつ"} という（$1$ は自明群）。このとき $G$ を $N$ の $H$ による**拡大**という。準同型 $s:H\to G$ で $\pi\circ s=\id_H$ を満たすものを $\pi$ の :term[切断]{en="section" yomi="せつだん"} といい、切断が存在するとき短完全列は :term[分裂する]{en="split" yomi="ぶんれつする"} という。
:::

全射 $\pi:G\to H$ が与えられれば、$N=\Ker\pi$ として包含写像により短完全列が得られる。各 $h\in H$ に対して $\pi(g)=h$ となる $g$ を一つずつ選ぶ写像は常に存在するが、それが準同型になるとは限らない。分裂するかどうかは、この選び方を積と両立させられるかという問題である。

:::theorem[分裂と補群]{#thm-grp-split}
全射準同型 $\pi:G\to H$ の核を $N$ とする。次は同値である。

1. $\pi$ は切断をもつ。
2. $N$ は補群をもつ。すなわち $NK=G$、$N\cap K=\{e\}$ を満たす部分群 $K$ が存在する。

このとき $s\mapsto s(H)$ は切断全体から補群全体への全単射であり、切断 $s$ について $\alpha_h(n)=s(h)\,n\,s(h)^{-1}$ とおけば $G\cong N\rtimes_\alpha H$ である。
:::

:::proof
切断 $s$ があれば $K=s(H)$ は部分群である。$\pi s=\id$ から $K\cap N=1$、また $g\,s(\pi(g))^{-1}\in N$ なので $G=NK$ である。逆に補群 $K$ があれば $\pi|_K$ は核が自明で、$\pi(NK)=\pi(K)=H$ より全射でもある。従って $s=(\pi|_K)^{-1}$ が切断である。この二つの構成は $s(H)=K$ と $\pi s=\id$ から互いに逆である。:ref[thm-grp-semidirect-internal] の共役作用を $K=s(H)$ で書き直せば指定の半直積になる。
:::

商が巡回群のときは、分裂の判定がとくに簡単になる。

:::proposition[巡回群への全射の分裂]{#prop-grp-cyclic-split}
$H=\langle h_0\rangle$ を位数 $m$ の巡回群、$\pi:G\to H$ を全射準同型とする。$\pi$ が切断をもつための必要十分条件は、$\pi(g)=h_0$ かつ $g^m=e$ を満たす $g\in G$ が存在することである。
:::

:::proof
切断 $s$ があれば $g=s(h_0)$ について $\pi(g)=h_0$、$g^m=s(h_0^m)=e$ である。逆にそのような $g$ があれば、$g^m=e$ なので $s(h_0^k)=g^k$ は指数の選び方によらず定まり、準同型 $H\to G$ を与える（[準同型と同型定理](/math/group-theory/homomorphisms/)の巡回群からの準同型の例）。$\pi(s(h_0^k))=\pi(g)^k=h_0^k$ だから $s$ は切断である。
:::

:::example[分裂する例としない例]{#exm-grp-split-examples}
1. $D_n$ の回転群による商 $C_2$ の生成元は鏡映 $s$ に持ち上がり、$s^2=e$ なので分裂する。
2. $S_3/A_3$ の生成元の持ち上げは三つの互換である。いずれも位数 $2$ なので切断は三つある。直積 $C_2^2$ でも第一成分の補群に第二成分と対角部分群の二つを選べる。
3. $C_4\to C_2$ で $\bar1$ の逆像は $\bar1,\bar3$ とも位数 $4$ であり、切断はない。
4. $Q_8/\langle i\rangle$ の非自明な類の全元も平方が $-1$ なので切断はない。どの非自明な部分群も $-1$ を含むため、別の二群を選んでも半直積にならない。
5. $\Z\to C_2$ は $\Z$ に位数 $2$ の元がないため分裂しない。
:::

### 自己同型群の計算

半直積 $N\rtimes_\alpha H$ を作るには準同型 $H\to\Aut(N)$ を列挙する必要がある。そのため基本的な群の自己同型群を求めておく。

:::proposition[巡回群の自己同型群]{#prop-grp-aut-cyclic}
$n\ge1$ とし、$C_n=\Z/n\Z$ を加法で書く。$u\in(\Z/n\Z)^\times$ に $\mu_u(x)=ux$ を対応させる写像は同型 $(\Z/n\Z)^\times\cong\Aut(C_n)$ を与える。特に $|\Aut(C_n)|=\varphi(n)$ であり、$\Aut(C_n)$ は可換である。
:::

:::proof
$f(\bar1)=u$ とすれば $f(\bar k)=ku$ であり、逆にこの式は任意の $u$ で準同型を与える。像は $\langle u\rangle$ だから、有限性と :ref[thm-grp-cyclic-order] により全単射となる条件は $\gcd(u,n)=1$ である。$\mu_u\mu_v=\mu_{uv}$ と $\mu_u(\bar1)=u$ より対応は合成を保つ全単射である。
:::

例えば法 $5,7,9$ ではそれぞれ $2,3,2$ の位数が $4,6,6$ なので、対応する自己同型群は巡回群である。法 $8$ の四つの単元はすべて平方が $1$ であり、$\Aut(C_8)\cong C_2^2$ となる。

:::proposition[基本 Abel 群の自己同型群]{#prop-grp-aut-elementary}
$p$ を素数とし、$C_p^n$ を $\F_p^n$ の加法群と同一視する。$C_p^n$ の自己準同型はちょうど $\F_p$ 線形写像であり、$\Aut(C_p^n)\cong\GL_n(\F_p)$ である。特に

$$
|\Aut(C_p^n)|=\prod_{i=0}^{n-1}(p^n-p^i).
$$

:::

:::proof
線形写像は和を保つので自己準同型である。逆に $f$ が和を保てば、$c\in\{0,1,\ldots,p-1\}$ について $f(cx)=f(x+\cdots+x)=cf(x)$ であり、$\F_p$ の元はこのような $c$ の類で尽くされるから $f$ は $\F_p$ 線形である。従って自己同型は可逆な線形写像であり、標準基底に関する行列を対応させれば $\GL_n(\F_p)$ と同型である。位数の式は[群の定義と例](/math/group-theory/groups/)で数えたとおりである（:eqref[eq-grp-gl-order] とその一般化）。
:::

$n=2,\ p=2$ では $|\Aut(C_2\times C_2)|=6$ である。自己同型は三つの非零元を置換し、二つの非零元の行き先を決めれば決まるので、$\Aut(C_2\times C_2)$ は三つの非零元の対称群 $S_3$ と同型である。$n=2,\ p=3$ では $|\Aut(C_3\times C_3)|=48$ である。

分類では $\Aut(N)$ の位数 $2$ の元（対合）がよく現れる。

:::lemma[対合の対角化]{#lem-grp-involution-diag}
体 $K$ で $2\ne0$ とし、$A\in\GL_n(K)$ が $A^2=I$ を満たすとする。$V_\pm=\{v\in K^n:Av=\pm v\}$ とおくと $K^n=V_+\oplus V_-$ であり、$A$ は $\diag(1,\ldots,1,-1,\ldots,-1)$（$-1$ が $\dim V_-$ 個）と $\GL_n(K)$ で共役である。二つの対合が共役であるための必要十分条件は $\dim V_-$ が等しいことである。
:::

:::proof
任意の $v$ について

$$
v=\tfrac12(v+Av)+\tfrac12(v-Av),\qquad A(v+Av)=Av+v,\qquad A(v-Av)=-(v-Av)
$$

だから $K^n=V_++V_-$ である。$v\in V_+\cap V_-$ なら $v=Av=-v$ で $2v=0$、従って $v=0$ である。$V_+$ と $V_-$ の基底を並べた正則行列 $P$ について $P^{-1}AP$ は上の対角行列になる。共役 $B=QAQ^{-1}$ は $V_-$ を $QV_-$ に移すので $\dim V_-$ は共役で変わらず、逆に $\dim V_-$ が等しい二つの対合は同じ対角行列と共役である。
:::

$K=\F_3,\ n=2$ では、対合は $\dim V_-=2$ の $-I$ と、$\dim V_-=1$ の $\diag(1,-1)$ の共役類（:ref[exm-grp-gl2f3-classes] によれば $12$ 個）であり、合計 $13$ 個ある。

### 半直積が同型になる条件

記号 $N\rtimes H$ だけでは群は決まらない。例えば $C_4\rtimes C_2$ で自明な作用を選べば可換群 $C_4\times C_2$、生成元を反転させる作用を選べば $D_4$ になる。一方、作用の式が異なって見えても、因子の生成元や座標を取り替えただけなら同型な群になる。この現象を定式化する。

:::proposition[作用の取り替えによる同型]{#prop-grp-semidirect-iso}
$\alpha:H\to\Aut(N)$、$\beta:H'\to\Aut(N')$ を準同型とし、同型 $u:N\to N'$、$v:H\to H'$ が

$$
u\circ\alpha_h\circ u^{-1}=\beta_{v(h)}\qquad(h\in H) \label{eq-grp-semidirect-compatible}
$$

を満たすとする。このとき $\Psi(n,h)=(u(n),v(h))$ は同型 $N\rtimes_\alpha H\to N'\rtimes_\beta H'$ である。特に次が成り立つ。

1. $w\in\Aut(H)$ について $N\rtimes_{\alpha}H\cong N\rtimes_{\alpha\circ w}H$ である。
2. $\sigma\in\Aut(N)$ について、$\beta_h=\sigma\alpha_h\sigma^{-1}$ で定まる $\beta$ に対し $N\rtimes_\alpha H\cong N\rtimes_\beta H$ である。
:::

:::proof
$\Psi$ は成分ごとに全単射なので全単射である。積を保つことは、:eqref[eq-grp-semidirect-compatible] から $u(\alpha_h(m))=\beta_{v(h)}(u(m))$ であることを用いて

$$
\Psi\bigl((n,h)(m,k)\bigr)=\bigl(u(n)\,u(\alpha_h(m)),\ v(h)v(k)\bigr)=\bigl(u(n)\,\beta_{v(h)}(u(m)),\ v(h)v(k)\bigr)=\Psi(n,h)\,\Psi(m,k)
$$

と確かめられる。1 は $u=\id_N$、$v=w^{-1}$、$\beta=\alpha\circ w$ とおけば条件 $\alpha_h=(\alpha\circ w)_{w^{-1}(h)}$ が成り立つ場合であり、2 は $u=\sigma$、$v=\id_H$ の場合である。
:::

$H$ が有限巡回群なら、作用の像が共役であることから半直積の同型が従う。そのために次の初等的な補題を用意する。

:::lemma[単数の持ち上げ]{#lem-grp-unit-lift}
正整数 $d,m$ は $d\mid m$ を満たすとし、整数 $k$ が $\gcd(k,d)=1$ を満たすとする。このとき $k'\equiv k\pmod d$ かつ $\gcd(k',m)=1$ を満たす整数 $k'$ が存在する。すなわち還元写像 $(\Z/m\Z)^\times\to(\Z/d\Z)^\times$ は全射である。
:::

:::proof
$m$ を割る素数のうち $k$ を割らないものすべての積を $t$ とし（そのような素数がなければ $t=1$）、$k'=k+dt$ とおく。$m$ を割る素数 $p$ をとる。$p\mid k$ の場合、$\gcd(k,d)=1$ より $p\nmid d$ であり、$t$ の定義から $p\nmid t$ なので $p\nmid dt$、従って $p\nmid k+dt$ である。$p\nmid k$ の場合、$p\mid t$ なので $p\mid dt$ であり、やはり $p\nmid k+dt$ である。よって $\gcd(k',m)=1$ である。
:::

:::proposition[巡回群による半直積の同型判定]{#prop-grp-semidirect-cyclic}
$H$ を有限巡回群とし、$\alpha,\beta:H\to\Aut(N)$ を準同型とする。像 $\alpha(H)$ と $\beta(H)$ が $\Aut(N)$ の共役な部分群であれば、$N\rtimes_\alpha H\cong N\rtimes_\beta H$ である。
:::

:::proof
$h$ を $H$ の生成元とする。像の共役を与える自己同型で座標を替えれば、:ref[prop-grp-semidirect-iso] により $\alpha(H)=\beta(H)$ の場合に帰着する。この共通の像は位数 $d\mid |H|$ の巡回群で、$\beta(h)=\alpha(h)^k$、$\gcd(k,d)=1$ と書ける。:ref[lem-grp-unit-lift] より $k$ と法 $d$ で等しく $|H|$ と互いに素な $k'$ をとれる。従って $h\mapsto h^{k'}$ は $H$ の自己同型 $w$ を定め、生成元での一致から $\beta=\alpha\circ w$ である。再び :ref[prop-grp-semidirect-iso] を適用する。
:::

$\Aut(N)$ が可換なら部分群の共役は自分自身だけなので、二つの作用の像が等しければ半直積は同型である。

:::example[位数 21 の非可換群]{#exm-grp-order-21}
$\Aut(C_7)\cong(\Z/7\Z)^\times$ は位数 $6$ の巡回群なので、位数 $3$ の部分群は $\{1,2,4\}$ だけである。$C_3$ の生成元を $2$ 倍で作用させる準同型と $4$ 倍で作用させる準同型は、像がともに $\{1,2,4\}$ なので、:ref[prop-grp-semidirect-cyclic] により同型な半直積 $C_7\rtimes C_3$ を与える。具体的には $C_3$ の生成元 $h$ を $h^2=h^{-1}$ に取り替える自己同型 $w$ によって、$2$ 倍の作用は $2^2=4$ 倍の作用に移る。従って非自明な作用による $C_7\rtimes C_3$ は同型を除いて一つである。同様に $C_5\rtimes C_4$ で生成元が $2$ 倍で作用するものと $3$ 倍で作用するものは、像がともに $\Aut(C_5)$ 全体なので同型である。
:::

:::remark
$H$ が巡回群でないと、像が等しくても同型とは限らない。$H=C_4\times C_2$、$N=C_3$ とし、$\alpha$ は $(\bar1,\bar0)$ を $-1$ 倍、$(\bar0,\bar1)$ を恒等写像に、$\beta$ は $(\bar1,\bar0)$ を恒等写像、$(\bar0,\bar1)$ を $-1$ 倍に送るとする。像はどちらも $\Aut(C_3)=\{\pm1\}$ である。$N\rtimes_\alpha H$ では $(0,(\bar0,\bar1))$ が中心に属し、:ref[prop-grp-internal-direct] により $N\rtimes_\alpha H\cong(C_3\rtimes_{-1}C_4)\times C_2$ となる。$C_3\rtimes_{-1}C_4$ の元 $(n,k)$ の平方は $(n+(-1)^kn,\ 2k)$ なので、位数 $2$ の元は $(0,\bar2)$ だけであり、$(C_3\rtimes_{-1}C_4)\times C_2$ の位数 $2$ の元は $3$ 個である。一方 $N\rtimes_\beta H\cong S_3\times C_4$ の位数 $2$ の元は、$x^2=e$ を満たす $S_3$ の $4$ 元と $C_4$ の $2$ 元の組から単位元を除いた $7$ 個である。従って二つの群は同型でない。
:::

### 半直積による分類

位数 $18,20,30$ について、正規部分群、補群、共役作用の順に求める。位数 $p^2$ の群の候補は $C_{p^2}$ と $C_p^2$ である（[Sylow の定理](/math/group-theory/sylow/)）。

:::proposition[位数 18 の群]{#prop-grp-order-18}
位数 $18$ の群は、次の五つのいずれか一つと同型である。

$$
C_{18},\qquad C_6\times C_3,\qquad D_9,\qquad S_3\times C_3,\qquad (C_3\times C_3)\rtimes_{-1}C_2.
$$

ここで最後の群は $C_2$ の生成元が $C_3\times C_3$ に $-1$ 倍で作用する半直積である。
:::

:::proof
Sylow $3$ 部分群 $P$ の個数は $2$ を割り $1\bmod3$ なので一つである。位数 $2$ の元 $b$ をとると $P\cap\langle b\rangle=1$、$|P\langle b\rangle|=18$ より $G=P\rtimes\langle b\rangle$ である。
$P=C_9$ なら $b$ の作用は $u^2\equiv1\bmod9$ を満たす $u$ 倍である。$u-1,u+1$ は同時に $3$ で割れないので、一方が $9$ の倍数である。従って作用は恒等または反転であり、直積 $C_{18}$ と二面体群 $D_9$ を得る。

$P=C_3^2$ なら作用の行列は平方が $I$ である。:ref[lem-grp-involution-diag] により固有値 $-1$ の空間の次元は $0,1,2$ のいずれかで、この次元で共役類が決まる。$0$ の場合は直積 $C_6\times C_3$、$2$ の場合は表題の反転作用である。$1$ なら $+1$ 固有直線は中心に入り、もう一つの直線と $b$ は $D_3\cong S_3$ を作る。固有直線の直和分解よりこれら二群の積は全体で交わりは自明だから、$C_3\times S_3$ を得る。

五群を区別するには、まず可換な二群のうち位数 $9$ の元があるのは $C_{18}$ だけと分かる。非可換な三群のうちでも位数 $9$ の元をもつのは $D_9$ だけである。実際、残る二群では任意の元の平方が $P=C_3^2$ に入り、その三乗が単位元になる。さらに $C_3\times S_3$ には位数 $2$ の元が三つあるのに対し、反転作用では $(v,b)^2=(v-v,e)$ となるため $P$ の外の九元全部が位数 $2$ である。
:::

:::proposition[位数 20 の群]{#prop-grp-order-20}
位数 $20$ の群は、次の五つのいずれか一つと同型である。

$$
C_{20},\qquad C_{10}\times C_2,\qquad D_{10},\qquad C_5\rtimes_{-1}C_4,\qquad \mathrm{Aff}(\F_5)\cong C_5\rtimes_2C_4.
$$

ここで $C_5\rtimes_uC_4$ は $C_4$ の生成元が $C_5$ に $u$ 倍で作用する半直積である。
:::

:::proof
Sylow $5$ 部分群は個数が $4$ を割り $1\bmod5$ なので唯一である。これを $N=C_5$ とし、位数 $4$ の Sylow 部分群を $H$ とすれば $N\cap H=1$ と $|NH|=20$ から $G=N\rtimes H$ である。$\Aut(N)$ は法 $5$ の単元 $2$ が生成する巡回群 $C_4$ である。

$H=C_4$ なら生成元の作用は $1,-1,2,3$ 倍のいずれかである。$2$ と $3$ は逆なので生成元を逆に取り替えれば同じ型となる。この三通りから順に $C_{20}$、$C_5\rtimes_{-1}C_4$、全アフィン群 $x\mapsto ax+b$（$a\in\F_5^\times,\ b\in\F_5$）を得る。
$H=C_2^2$ なら像は $\{1,-1\}$ に含まれる。自明な作用は $C_{10}\times C_2$ を与える。非自明なら核の生成元と、その外の一点を $H$ の基底に選べる。前者は中心に入り、後者は $N$ と $D_5$ を作る。従って群は $C_2\times D_5\cong D_{10}$ である（:ref[exm-grp-internal-direct] の 3）。

可換な二群では位数 $20$ の元の有無が異なる。非可換な三群は位数 $2$ の元を数えて区別できる。$D_{10}$ では半回転と十個の鏡映、$C_5\rtimes_{-1}C_4$ では $(v,j)^2=(v+(-1)^jv,2j)$ より $(0,2)$ のみである。アフィン群では二乗が $x\mapsto a^2x+(1+a)b$ だから、非恒等で平方が恒等となるのは $a=-1$ の五個である。
:::

:::proposition[位数 30 の群]{#prop-grp-order-30}
位数 $30$ の群は位数 $15$ の巡回正規部分群をもち、次の四つのいずれか一つと同型である。

$$
C_{30},\qquad C_3\times D_5,\qquad C_5\times S_3,\qquad D_{15}.
$$

:::

:::proof
Sylow の個数条件から $n_3=1,10$、$n_5=1,6$ が候補となる。両方とも $1$ より大きければ、相異なる素数位数の部分群は自明に交わるため、位数 $3,5$ の元が $10\cdot2+6\cdot4=44$ 個必要になり矛盾する。一方が正規なので Sylow 部分群の積 $M=P_3P_5$ は位数 $15$ の部分群であり、指数 $2$ だから正規である。$M$ の中では $n_3\mid5,\ n_5\mid3$ と合同条件から両方とも $1$ となる。二つの部分群は正規で自明に交わるので $M=C_3\times C_5\cong C_{15}$ である。

位数 $2$ の元が $M$ の補群を生成する。その共役作用は、位数 $3,5$ の特性部分群をそれぞれ保つ。各成分上で平方が恒等となる単元は $+1,-1$ のみだから、二つの符号を選ぶ四通りである。両方 $+$ なら $C_{30}$、両方 $-$ なら $D_{15}$ となる。一方だけ反転する場合、固定する因子は中心に入り、反転する因子は補群と二面体群を作るので、$C_3\times D_5$ または $C_5\times D_3$ となる。奇数 $m\ge3$ に対し $Z(D_m)=1$（:ref[def-grp-center] に続く二面体群の計算）なので、四群の中心の位数は順に $30,1,3,5$ と異なる。
:::

:::remark
同じ方法で、素数 $p<q$ に対する位数 $pq$ の群も分類できる（[Sylow の定理](/math/group-theory/sylow/)の章で扱った）。Sylow $q$ 部分群が正規で $G\cong C_q\rtimes C_p$ となり、$\Aut(C_q)\cong(\Z/q\Z)^\times$ の位数 $p$ の部分群は高々一つなので、:ref[prop-grp-semidirect-cyclic] により非自明な作用はすべて同型な群を与える。位数 $21$ はその例である（:ref[exm-grp-order-21]）。
:::

## 有限生成 Abel 群の基本定理

本節と次節では Abel 群を加法で書き、単位元を $0$、$x$ の $k$ 倍を $kx$ と書く。部分集合 $S$ の生成する部分群 $\langle S\rangle$ は $S$ の元の整数係数の有限和全体である。有限集合で生成される Abel 群を :term[有限生成]{en="finitely generated" yomi="ゆうげんせいせい"} という。Abel 群なので部分群はすべて正規であり、剰余群が常に作れる。

本節の目標は、有限生成 Abel 群が巡回群の直和に分解し、その分解の型が群だけで決まることを示すことである。証明は三段階からなる。まず生成元を選んで群を $\Z^n$ の剰余群として表す。次に、割る部分群の生成系を並べた整数行列を基本変形で対角化する（Smith 標準形）。最後に、得られた分解の型が群の不変量で読み取れることを示す。

### 捩れ部分群と自由 Abel 群

:::definition[捩れ部分群]{#def-grp-torsion}
Abel 群 $A$ において、ある正整数 $n$ で $nx=0$ となる元（有限位数の元）全体を $A$ の :term[捩れ部分群]{en="torsion subgroup" yomi="ねじれぶぶんぐん"} といい、$T(A)$ と書く。$T(A)=0$ の群を**捩れのない群**、$T(A)=A$ の群を**捩れ群**という。
:::

$T(A)$ は実際に部分群である。$0\in T(A)$ であり、$mx=0,\ ny=0$ なら $mn(x-y)=n(mx)-m(ny)=0$ だからである。ここで可換性を用いた。非可換群では有限位数の元の全体は部分群とは限らない（演習問題を参照）。例えば $T(\Z\oplus C_6)=\{0\}\oplus C_6$、$T(\R/\Z)=\Q/\Z$ であり、$\C^\times$ の捩れ部分群は 1 の冪根全体である。$\Z^n$ や $\Q$ は捩れのない群である。

:::definition[自由 Abel 群と基底]{#def-grp-free-abelian}
Abel 群 $F$ の部分集合 $B$ が**基底**であるとは、各元が $B$ の元の整数係数の和として一意に表されることをいう。係数は有限個を除いて零とする。基底をもつ群を :term[自由 Abel 群]{en="free abelian group" yomi="じゆうあーべるぐん"} という。以下の構造定理で扱う自由部分には有限基底がある。
:::

$\Z^r$ の標準基底 $e_i=(0,\ldots,1,\ldots,0)$ は基底であり、基底 $e_1,\ldots,e_r$ をもつ $F$ は座標写像 $\sum k_ie_i\mapsto(k_1,\ldots,k_r)$ によって $\Z^r$ と同型である。以下、$\Z^n$ の元は**行ベクトル**として書く。

:::proposition[自由 Abel 群の普遍性]{#prop-grp-free-universal}
$F$ を基底 $e_1,\ldots,e_r$ をもつ自由 Abel 群とする。任意の Abel 群 $B$ と $b_1,\ldots,b_r\in B$ に対し、$f(e_i)=b_i$ を満たす準同型 $f:F\to B$ がただ一つ存在する。特に、$n$ 個の元 $a_1,\ldots,a_n$ で生成される Abel 群 $A$ は、$K=\{(k_i)\in\Z^n:\sum k_ia_i=0\}$ とおくと $A\cong\Z^n/K$ である。
:::

:::proof
$f(\sum k_ie_i)=\sum k_ib_i$ とおく。表示が一意なのでこれは写像として定まり、$B$ の可換性から和を保つ。準同型は整数倍と和を保つので、$e_i$ での値が決まれば全体で決まり、一意性が従う。後半は $e_i\mapsto a_i$ で定まる $\Z^n\to A$ が全射で核が $K$ であることに :ref[thm-grp-first-iso] を適用すればよい。
:::

$\Z^n$ の部分群 $K$ は、生成元 $a_i$ の間に成り立つ**関係式** $\sum k_ia_i=0$ の全体を表している。

:::proposition[自由 Abel 群の階数]{#prop-grp-rank-invariance}
$\Z^m\cong\Z^n$ なら $m=n$ である。従って有限階数の自由 Abel 群の基底の元の個数は基底の選び方によらない。これを自由 Abel 群の :term[階数]{en="rank" yomi="かいすう"} という。
:::

:::proof
同型 $f:\Z^m\to\Z^n$ は $f(2x)=2f(x)$ と全射性により $2\Z^m$ を $2\Z^n$ の上に移すので、同型 $\Z^m/2\Z^m\cong\Z^n/2\Z^n$ を誘導する。$\Z^m/2\Z^m$ は各座標を $2$ を法として見ることで $(\Z/2\Z)^m$ と同型であり、$2^m$ 個の元をもつ。従って $2^m=2^n$ である。
:::

基底の取り替えは整数行列で記述される。$n$ 次の整数正方行列 $P$ で、逆行列も整数行列であるもの全体を $\GL_n(\Z)$ と書く。:ref[prop-linalg-det-ring] の 4 により、整数正方行列 $P$ が $\GL_n(\Z)$ に属するための必要十分条件は $\det P=\pm1$ である。

:::lemma[Z^n の基底と GL_n(Z)]{#lem-grp-integer-basis}
$P\in M_n(\Z)$ の行ベクトル $f_1,\ldots,f_n$ が $\Z^n$ の基底であるための必要十分条件は、$P\in\GL_n(\Z)$ であることである。
:::

:::proof
行ベクトルで書くと、整数係数の一次結合 $\sum y_if_i$ は $yP$（$y\in\Z^n$）である。$P\in\GL_n(\Z)$ で逆行列を $Q$ とすれば、任意の $x$ は $x=(xQ)P$ と表され、$yP=y'P$ なら右から $Q$ を掛けて $y=y'$ である。逆に $f_i$ が基底なら、各 $e_i$ を $e_i=q_iP$ と表して $q_i$ を行とする整数行列 $Q$ を作ると $QP=I$ である。行列式をとって $\det Q\det P=1$、従って $\det P=\pm1$ であり、$P\in\GL_n(\Z)$ となる。
:::

:::lemma[整数格子の部分群]{#lem-grp-lattice}
$\Z^r$ の任意の部分群 $B$ は、高々 $r$ 個の元からなる基底をもつ自由 Abel 群である。
:::

:::proof
$r=0$ なら空基底をとる。$r\ge1$ とし、第一座標への射影の像を $d\Z$（$d\ge0$）と書く。核 $B_0$ は $\Z^{r-1}$ の部分群だから、帰納法により高々 $r-1$ 個の元からなる基底をもつ。$d=0$ なら $B=B_0$ である。$d>0$ なら第一座標が $d$ の $b\in B$ を選ぶ。第一座標が $kd$ の $x\in B$ は $x=kb+(x-kb)$ と書け、括弧内は $B_0$ に入る。第一座標から $k$ が一意に決まるので、$b$ を $B_0$ の基底に加えれば $B$ の基底となる。
:::

:::corollary[有限生成 Abel 群の部分群]{#cor-grp-fg-subgroup}
$n$ 個の元で生成される Abel 群の部分群は、$n$ 個以下の元で生成される。
:::

:::proof
:ref[prop-grp-free-universal] の全射 $\varphi:\Z^n\to A$ をとる。部分群 $B\le A$ の逆像 $\varphi^{-1}(B)$ は $\Z^n$ の部分群なので、:ref[lem-grp-lattice] により $n$ 個以下の元で生成される。その像が $B$ を生成する。
:::

### Smith 標準形

整数行列 $M\in M_{m,n}(\Z)$ の行ベクトルの整数係数の一次結合全体を $R(M)=\{yM:y\in\Z^m\}\le\Z^n$ と書く。生成元 $x_1,\ldots,x_n$ の間の関係式 $\sum_ja_{ij}x_j=0$（$1\le i\le m$）を課した Abel 群とは $\Z^n/R(M)$、$M=(a_{ij})$ のことであり、$M$ をその**関係行列**という。

整数行列に対する次の操作を**整数基本変形**という。

1. 二つの行（列）を入れ替える。
2. 一つの行（列）を $-1$ 倍する。
3. 一つの行（列）に、別の行（列）の整数倍を加える。

体上の基本変形との違いは、2 で $\pm1$ 以外の数を掛けないことである。$\Z$ では $\pm1$ 以外の数は逆数をもたないからである。行の変形は左から、列の変形は右から、それぞれ対応する基本行列を掛けることで実現され、基本行列の行列式は $\pm1$、逆行列も同じ型の基本行列である。従って整数基本変形を繰り返すことは、$M$ を $PMQ$（$P\in\GL_m(\Z),\ Q\in\GL_n(\Z)$）に変えることである。

:::lemma[関係行列の変形]{#lem-grp-relation-change}
$P\in\GL_m(\Z)$、$Q\in\GL_n(\Z)$ なら $\Z^n/R(M)\cong\Z^n/R(PMQ)$ である。
:::

:::proof
$yM=(yP^{-1})(PM)$ と $y(PM)=(yP)M$ から $R(PM)=R(M)$ である。次に $\phi(v)=vQ$ は $\Z^n$ の自己同型（逆写像は $v\mapsto vQ^{-1}$）であり、$\phi(R(M))=\{yMQ\}=R(MQ)$ である。合成 $\Z^n\xrightarrow{\phi}\Z^n\to\Z^n/R(MQ)$ は全射で、核は $\phi^{-1}(R(MQ))=R(M)$ である。:ref[thm-grp-first-iso] より $\Z^n/R(M)\cong\Z^n/R(MQ)=\Z^n/R(PMQ)$ を得る。
:::

行の変形は関係式の組を取り替えること（関係の全体は変わらない）、列の変形は生成元を取り替えることに当たる。例えば第 $j$ 列の $c$ 倍を第 $l$ 列に加える変形は、新しい生成元 $x_j'=x_j-cx_l$ を用いて関係式を書き直すことに対応する。

:::theorem[Smith 標準形]{#thm-grp-smith}
任意の $M\in M_{m,n}(\Z)$ に対し、整数基本変形の積として得られる $P\in\GL_m(\Z)$、$Q\in\GL_n(\Z)$ があって、$PMQ$ は $(i,i)$ 成分（$1\le i\le k$）が $d_i$、その他の成分が $0$ の行列になる。ここで $k\le\min(m,n)$ であり、$d_i$ は正整数で

$$
d_1\mid d_2\mid\cdots\mid d_k
$$

を満たす。この行列を $M$ の :term[Smith 標準形]{en="Smith normal form" yomi="すみすひょうじゅんけい"} という。
:::

:::proof
$m+n$ に関する帰納法で示す。行または列の数が $0$ の行列、および零行列は、$k=0$ としてすでに標準形である。$M\ne0$ とし、次の手順を考える。

(a) 零でない成分のうち絶対値が最小のものを、行と列の入れ替えで $(1,1)$ 成分に移し、必要なら第 $1$ 行を $-1$ 倍して正にする。これを $a$ とする。

(b) 第 $1$ 行のある成分 $a_{1j}$ が $a$ で割り切れなければ、$a_{1j}=qa+\rho$（$0<\rho<a$）と割り、第 $j$ 列から第 $1$ 列の $q$ 倍を引く。すると $(1,j)$ 成分は $\rho$ になり、零でない成分の絶対値の最小値は $a$ より真に小さくなる。このときは (a) に戻る。第 $1$ 列の成分についても行の変形で同様にする。

(c) 第 $1$ 行と第 $1$ 列の成分がすべて $a$ で割り切れるなら、第 $1$ 列の整数倍を各列から、第 $1$ 行の整数倍を各行から引いて、$(1,1)$ 以外の第 $1$ 行と第 $1$ 列の成分をすべて $0$ にする。

(d) さらに、ある $(i,j)$ 成分（$i,j\ge2$）が $a$ で割り切れなければ、第 $i$ 行を第 $1$ 行に加える。第 $i$ 行の第 $1$ 成分は $0$ なので、新しい第 $1$ 行は $(a,a_{i2},\ldots,a_{in})$ となり、$(1,j)$ 成分が $a$ で割り切れない。そこで (b) を行い、(a) に戻る。

(a) に戻るたびに、零でない成分の絶対値の最小値は真に減少する正整数である。従って (a) に戻るのは有限回であり、手順はいずれ (c) を終えた時点で (d) の条件が満たされない状態に達する。そのとき行列は

$$
\begin{pmatrix}a&0\\0&M'\end{pmatrix},\qquad M'\in M_{m-1,n-1}(\Z)
$$

の形で、$M'$ の成分はすべて $a$ の倍数である。帰納法の仮定により $P'M'Q'$ が標準形 $\diag(d_2,\ldots,d_k)$ となる $P',Q'$ がある。$M'$ に施す整数基本変形は、成分を元の成分の整数係数の和に変えるだけなので、「すべての成分が $a$ の倍数である」という性質を保つ。従って $a\mid d_2$ である（$k\ge2$ の場合）。$\diag(1,P')$ と $\diag(1,Q')$ は整数基本変形の積であり（$M'$ に対する変形を第 $2$ 行・列以降に施せばよい）、それらを掛けて $d_1=a$ とすれば求める形を得る。
:::

Smith 標準形の対角成分 $d_i$ は $M$ だけで決まる。このことは次の基本定理の一意性から従う（:ref[thm-grp-fg-abelian] の証明の後の注意を参照）。小さい行列では、$d_1$ が全成分の最大公約数、$d_1d_2\cdots d_j$ が $j$ 次小行列式全体の最大公約数に等しいことを使って検算できる（演習問題）。

:::corollary[部分群に適合する基底]{#cor-grp-stacked-basis}
$K\le\Z^n$ を部分群とする。$\Z^n$ の基底 $f_1,\ldots,f_n$ と、$d_1\mid d_2\mid\cdots\mid d_k$ を満たす正整数 $d_1,\ldots,d_k$（$k\le n$）があって、$d_1f_1,\ldots,d_kf_k$ は $K$ の基底である。
:::

:::proof
:ref[lem-grp-lattice] により $K$ は有限個の元で生成される。それらを行とする行列 $M$ をとると $K=R(M)$ である（$K=0$ なら $M$ を零行列にとる）。:ref[thm-grp-smith] の $PMQ=D$ について、$f_i$ を $Q^{-1}$ の第 $i$ 行とすると、:ref[lem-grp-integer-basis] により $f_1,\ldots,f_n$ は基底である。$R(M)=R(PM)=R(DQ^{-1})$ であり、$DQ^{-1}$ の第 $i$ 行は $i\le k$ なら $d_if_i$、それ以外は $0$ である。よって $K$ は $d_1f_1,\ldots,d_kf_k$ で生成される。これらの一次結合 $\sum c_id_if_i$ が $0$ なら、$f_i$ の基底性から $c_id_i=0$、従って $c_i=0$ なので、表示は一意である。
:::

### 基本定理

:::theorem[有限生成 Abel 群の基本定理]{#thm-grp-fg-abelian}
$A$ を有限生成 Abel 群とする。

1. 整数 $r\ge0$ と、$1<d_1\mid d_2\mid\cdots\mid d_s$ を満たす整数 $d_1,\ldots,d_s$（$s\ge0$）があって

$$
A\cong\Z^r\oplus\Z/d_1\Z\oplus\cdots\oplus\Z/d_s\Z \label{eq-grp-invariant-factor-form}
$$

となる。$r$ と $d_1,\ldots,d_s$ は $A$ の同型類だけで決まる。

2. 同値な形として、素数冪 $q_1,\ldots,q_t>1$ があって $A\cong\Z^r\oplus\Z/q_1\Z\oplus\cdots\oplus\Z/q_t\Z$ と表され、$r$ と多重集合 $\{q_1,\ldots,q_t\}$ は $A$ の同型類だけで決まる。

$r$ を $A$ の**階数**（自由階数）、$d_i$ を :term[不変因子]{en="invariant factor" yomi="ふへんいんし"}、$q_j$ を :term[単因子]{en="elementary divisor" yomi="たんいんし"} という。$A$ が有限であることと $r=0$ は同値である。
:::

:::proof[存在の証明]
$A$ が $n$ 個の元で生成されるとする。:ref[prop-grp-free-universal] により $A\cong\Z^n/K$ である。:ref[cor-grp-stacked-basis] の基底 $f_1,\ldots,f_n$ と $d_1\mid\cdots\mid d_k$ をとる。座標写像

$$
\Z^n\to\Z/d_1\Z\oplus\cdots\oplus\Z/d_k\Z\oplus\Z^{n-k},\qquad \sum_ic_if_i\mapsto(c_1\bmod d_1,\ \ldots,\ c_k\bmod d_k,\ c_{k+1},\ldots,c_n)
$$

は全射準同型であり、核は $c_i\in d_i\Z$（$i\le k$）かつ $c_i=0$（$i>k$）となる元、すなわち $K$ である。従って $A\cong\bigoplus_{i\le k}\Z/d_i\Z\oplus\Z^{n-k}$ である。$d_i=1$ の因子は零群なので除いてよく、残りは整除の列をなす。これで 1 の分解を得る（$r=n-k$）。各 $\Z/d_i\Z$ を :ref[exm-grp-crt] のように素数冪位数の巡回群の直和に分解すれば 2 の分解を得る。一意性は以下の補題で示す。
:::

:::lemma[自由部分の分離]{#lem-grp-free-splitting}
有限生成 Abel 群 $A$ について、$T(A)$ は有限群であり、$A/T(A)$ は自由 Abel 群で、$A\cong T(A)\oplus A/T(A)$ である。:eqref[eq-grp-invariant-factor-form] のような任意の分解 $A\cong\Z^r\oplus F$（$F$ 有限）について、$F\cong T(A)$ かつ $r$ は $A/T(A)$ の階数に等しい。特に $r$ は $A$ だけで決まる。
:::

:::proof
存在の証明で得た $\Z^r\oplus F$ において、有限位数の元はちょうど $(0,f)$ である。実際、$nv=0$（$n>0$）を満たす整数ベクトルは $v=0$ に限られ、$F$ の全元は有限位数である。同型は元の位数を保つので $T(A)\cong F$、射影の核をとれば $A/T(A)\cong\Z^r$ となる。これで直和の主張が従い、:ref[prop-grp-rank-invariance] より $r$ は群だけで決まる。
:::

$T(A)$ の補群（$A/T(A)$ と同型な部分群で $T(A)$ と自明に交わり、合わせて全体になるもの）は一意ではない。$A=\Z\oplus C_2$ で $A/T(A)\cong\Z$ の生成元を $(1,\bar0)$ に持ち上げても $(1,\bar1)$ に持ち上げてもよく、それぞれが生成する $\Z\oplus\{\bar0\}$ と $\{(k,\bar k)\}$ はともに補群であるが、異なる部分集合である。基本定理の一意性は、部分群としての分解の一意性ではなく、因子の**型**の一意性である。

:::lemma[準素分解]{#lem-grp-primary}
有限 Abel 群 $T$ と素数 $p$ に対し、$T_p=\{x\in T:\text{ある }k\ge0\text{ で }p^kx=0\}$ を $T$ の :term[p 準素成分]{en="p-primary component" yomi="ぴーじゅんそせいぶん"} という。$T_p$ は位数が $p$ の冪の部分群であり、$|T|$ を割る素数の全体を $p_1,\ldots,p_l$ とすると

$$
T=T_{p_1}\oplus\cdots\oplus T_{p_l}
$$

（内部直和）となる。同型 $\theta:T\to T'$ は $T_p$ を $T'_p$ の上に移す。また $T=\bigoplus_j\Z/q_j\Z$（$q_j$ は素数冪）なら、$T_p$ は $q_j$ が $p$ の冪である因子の直和である。
:::

:::proof
$p^kx=0,\ p^ly=0$ なら $p^{\max(k,l)}(x-y)=0$ なので $T_p$ は部分群である。$|T_p|$ を割る素数 $q\ne p$ があれば、:ref[thm-grp-cauchy] により $T_p$ は位数 $q$ の元 $x$ をもつが、$p^kx=0$ と $qx=0$ から（$up^k+vq=1$ をとって）$x=0$ となり矛盾する。従って $|T_p|$ は $p$ の冪である。

$|T|=1$ なら空直和である。それ以外では $M=|T|=\prod_i p_i^{a_i}$、$M_i=M/p_i^{a_i}$ とおく。$\gcd(M_1,\ldots,M_l)=1$ より $\sum_i c_iM_i=1$ と書ける。$Mx=0$ だから $x=\sum_i c_iM_ix$ の各項は $p_i^{a_i}$ 倍で消える。他の成分の和は $M_i$ 倍で消えるため、$T_{p_i}$ との共通部分の元は互いに素な $M_i,p_i^{a_i}$ の双方で消え、零である。従ってこの和は直和である。

同型は $p^kx=0$ という条件を保つので、$T_p$ を $T'_p$ に移す。最後に、$\Z/\ell^e\Z$（$\ell$ は素数）の元 $x$ が $p^kx=0$ を満たすとき、$\ell\ne p$ なら $\ell^ex=0$ と合わせて $x=0$ であり、$\ell=p$ なら全元が条件を満たす。直和の元が $T_p$ に属することは各成分が属することと同値なので、最後の主張を得る。
:::

有限 Abel $p$ 群については、$p^j$ 倍で消える元の個数が分解の型を決める。Abel 群 $A$ と正整数 $m$ に対し $A[m]=\{x\in A:mx=0\}$ と書く。これは $m$ 倍写像の核なので部分群であり、同型 $\theta:A\to A'$ は $A[m]$ を $A'[m]$ の上に移す。

:::proposition[p 群の型の決定]{#prop-grp-p-type}
$P\cong\Z/p^{e_1}\Z\oplus\cdots\oplus\Z/p^{e_t}\Z$（$e_i\ge1$）とすると、$|P[p^j]|=p^{\sum_i\min(j,e_i)}$ であり、

$$
\log_p\frac{|P[p^j]|}{|P[p^{j-1}]|}=\#\{i:e_i\ge j\}\qquad(j\ge1). \label{eq-grp-primary-layers}
$$

従って多重集合 $\{e_1,\ldots,e_t\}$ は $P$ の同型類だけで決まる。
:::

:::proof
$\Z/p^e\Z$ の元 $x$（$0\le x<p^e$）が $p^jx\equiv0\pmod{p^e}$ を満たす条件は、$j\ge e$ なら常に成り立ち、$j<e$ なら $p^{e-j}\mid x$ である。いずれの場合も解は $p^{\min(j,e)}$ 個である。直和では $(\bigoplus C_i)[m]=\bigoplus C_i[m]$ なので個数の式を得る。$\min(j,e_i)-\min(j-1,e_i)$ は $e_i\ge j$ なら $1$、そうでなければ $0$ だから :eqref[eq-grp-primary-layers] が従う。左辺は $P$ の同型類で決まるので、$\#\{i:e_i\ge j\}$ がすべての $j$ について決まり、その差 $\#\{i:e_i=j\}=\#\{i:e_i\ge j\}-\#\{i:e_i\ge j+1\}$ も決まる。
:::

:::proof[基本定理の一意性の証明]
$r$ の一意性は :ref[lem-grp-free-splitting] で示した。同じ補題により、分解の有限部分 $\bigoplus_j\Z/q_j\Z$ はいつも $T(A)$ と同型である。:ref[lem-grp-primary] により、$q_j$ のうち $p$ の冪であるものの直和は $T(A)_p$ と同型であり、:ref[prop-grp-p-type] によってその指数の多重集合は $T(A)_p$ で決まる。従って単因子の多重集合は一意である。

捩れ部分が自明なら不変因子の列は空である。そうでなければ $d_1>1$ の素因数はすべての $d_i$ を割る。従って不変因子の個数 $s$ は、一つの素数に属する単因子の個数の最大値である。各素数 $p$ の指数を昇順に並べ、左端を $0$ で補って長さ $s$ にする。$d_i\mid d_{i+1}$ よりこの配列が各 $d_i$ の $p$ 進付値を与える唯一の順序であり、各列で素数冪を掛ければ $d_i$ が復元される。
:::

:::remark
一意性から Smith 標準形の対角成分も $M$ で決まる。実際、$\Z^n/R(M)\cong\bigoplus_{i\le k}\Z/d_i\Z\oplus\Z^{n-k}$ で、$1$ でない $d_i$ は不変因子、$n-k$ は階数として一意であり、$1$ に等しい $d_i$ の個数は $k$ から $1$ でないものの個数を引いた数である。また「単因子」という語は、行列の理論では Smith 標準形の $d_i$（本章の不変因子）を指すことが多い。文献を読む際には定義を確かめる必要がある。
:::

### 別証明：最大位数の巡回因子

有限 Abel 群に限れば、行列を使わずに分解の存在を示すこともできる。:ref[lem-grp-primary] により有限 Abel $p$ 群を扱えば十分である。

:::lemma[最大位数の巡回因子を取り出す]{#lem-grp-cyclic-summand}
$P$ を非自明な有限 Abel $p$ 群、$a\in P$ を位数最大の元とし、その位数を $p^e$ とする。このとき $P=\langle a\rangle\oplus K$ を満たす部分群 $K$ が存在する。
:::

:::proof
$C=\Z/p^e\Z$ とし、$F(a)=\bar1$ を満たす準同型 $F:P\to C$ を構成する。$B_0=\langle a\rangle$ 上では $f_0(ca)=\bar c$ とおく。$a$ の位数が $p^e$ なので、これは定義でき準同型である。

部分群 $B\supseteq\langle a\rangle$ 上の準同型 $f:B\to C$ で $f(a)=\bar1$ を満たすものが与えられ、$x\in P\setminus B$ とする。$x+B$ の $P/B$ における位数を $p^k$ とすると $k\ge1$ であり、$x$ の位数は $p^e$ 以下なので $k\le e$ である。$p^kx\in B$ であり、$p^{e-k}f(p^kx)=f(p^ex)=0$ である。$C$ の元 $z$ で $p^{e-k}z=0$ を満たすものはちょうど $p^k$ の倍数なので、$f(p^kx)=p^ky$ となる $y\in C$ がある。

$B'=B+\langle x\rangle$ 上で $f'(b+tx)=f(b)+ty$ と定める。$b+tx=b'+t'x$ なら $(t-t')x=b'-b\in B$ なので $p^k\mid t-t'$ であり、$t-t'=lp^k$ と書けば $b'-b=lp^kx$ である。すると

$$
\bigl(f(b)+ty\bigr)-\bigl(f(b')+t'y\bigr)=-f(lp^kx)+lp^ky=-l\bigl(f(p^kx)-p^ky\bigr)=0
$$

なので $f'$ は定義でき、明らかに準同型で $f$ を延長する。$B'$ は $B$ より真に大きく、$P$ は有限なので、この操作を有限回繰り返して $F:P\to C$ を得る。

$K=\Ker F$ とおく。任意の $v\in P$ について $F(v)=\bar c$ とすれば $F(v-ca)=0$ なので $v\in\langle a\rangle+K$ である。$ca\in K$ なら $\bar c=F(ca)=0$、すなわち $p^e\mid c$ で $ca=0$ だから $\langle a\rangle\cap K=0$ である。:ref[prop-grp-internal-direct] より $P=\langle a\rangle\oplus K$ である。
:::

$|K|<|P|$ なので、$K$ に同じ補題を繰り返し適用すると、有限 Abel $p$ 群が巡回群の直和になることが位数に関する帰納法で示される。

位数最大という仮定は外せない。$P=C_8\oplus C_4$ で $a=(\bar1,\bar1)$（位数 $8$）をとると、第一成分への射影 $F$ が $F(a)=\bar1$ を満たし、各 $(x,y)$ は $(x,y)=xa+(0,y-x)$ と分解される（第二成分で $x$ を法 $4$ に落とせるのは $4\mid8$ による）。一方 $a'=(\bar2,\bar0)$ は位数 $4$ だが、$\langle a'\rangle$ は直和因子でない。直和因子なら射影 $P\to\langle a'\rangle\cong C_4$ で $a'\mapsto\bar1$ となるものがあるが、$a'=2(\bar1,\bar0)$ の像は $2F'((\bar1,\bar0))$ で偶数の類であり、$\bar1$ にはならないからである。同様に $C_4$ の位数 $2$ の部分群は直和因子でない。

## 有限 Abel 群の分類

基本定理により、有限 Abel 群は単因子の多重集合、あるいは不変因子の列で完全に分類される。本節では、この分類を具体的な数え上げと計算に用いる。

### 位数を与えたときの分類

正整数 $k$ を正整数の和として表す方法（順序を区別しない）の個数を $k$ の**分割数** $\mathrm{p}(k)$ という。例えば $3=3=2+1=1+1+1$ なので $\mathrm{p}(3)=3$ である。小さい値は

| $k$ | $1$ | $2$ | $3$ | $4$ | $5$ | $6$ |
| --- | --- | --- | --- | --- | --- | --- |
| $\mathrm{p}(k)$ | $1$ | $2$ | $3$ | $5$ | $7$ | $11$ |

である。

:::proposition[有限 Abel 群の個数]{#prop-grp-abelian-count}
$n=p_1^{a_1}\cdots p_l^{a_l}$（$p_i$ は相異なる素数）とする。位数 $n$ の Abel 群の同型類は、各 $i$ について $a_i$ の分割を一つずつ選ぶことと一対一に対応する。特にその個数は $\mathrm{p}(a_1)\cdots\mathrm{p}(a_l)$ である。
:::

:::proof
位数 $n$ の Abel 群 $A$ の単因子のうち $p_i$ の冪であるものを $p_i^{e_1},\ldots,p_i^{e_t}$ とすると、:ref[lem-grp-primary] により $|A_{p_i}|=p_i^{e_1+\cdots+e_t}$ であり、$|A|=\prod_i|A_{p_i}|$ だから $e_1+\cdots+e_t=a_i$ である。つまり単因子は各 $a_i$ の分割を定める。逆に分割を任意に選べば、対応する巡回群の直和は位数 $n$ の Abel 群である。:ref[thm-grp-fg-abelian] の一意性から、異なる選び方は同型でない群を与える。
:::

例えば位数 $p^4$ の Abel 群は $4=4=3+1=2+2=2+1+1=1+1+1+1$ に対応して

$$
C_{p^4},\quad C_{p^3}\oplus C_p,\quad C_{p^2}\oplus C_{p^2},\quad C_{p^2}\oplus C_p\oplus C_p,\quad C_p^4
$$

の $\mathrm{p}(4)=5$ 個である。

:::example[位数 72 の Abel 群]{#exm-grp-order-72}
$72=2^3\cdot3^2$ なので $\mathrm{p}(3)\mathrm{p}(2)=6$ 個ある。$2$ 部分は $C_8,\ C_2\oplus C_4,\ C_2\oplus C_2\oplus C_2$、$3$ 部分は $C_9,\ C_3\oplus C_3$ である。不変因子に直すには、:ref[thm-grp-fg-abelian] の一意性の証明のように、各素数の指数を大きい順に並べて同じ順位どうしを掛ける。

| $2$ 部分 | $C_9$ と合わせる | $C_3\oplus C_3$ と合わせる |
| --- | --- | --- |
| $C_8$ | $C_{72}$ | $C_3\oplus C_{24}$ |
| $C_2\oplus C_4$ | $C_2\oplus C_{36}$ | $C_6\oplus C_{12}$ |
| $C_2\oplus C_2\oplus C_2$ | $C_2\oplus C_2\oplus C_{18}$ | $C_2\oplus C_6\oplus C_6$ |

例えば右下は、$2$ 部分の指数 $(1,1,1)$ と $3$ 部分の指数 $(1,1)$ を大きい順に揃えて $(2\cdot3,\ 2\cdot3,\ 2\cdot1)$ とし、整除の順に並べ直したものである。
:::

逆方向の変換は素因数分解で行う。例えば $C_6\oplus C_{10}$ と $C_2\oplus C_{30}$ はどちらも単因子 $2,2,3,5$ をもつので同型であり、不変因子は $(2,30)$ である。$C_{12}\oplus C_{18}$ の単因子は $4,3,2,9$ なので、不変因子は $(2\cdot3,\ 4\cdot9)=(6,36)$ である。

### 元と部分群を数える

:::proposition[倍写像の核と方程式の解]{#prop-grp-abelian-kernel-count}
$A=\Z/d_1\Z\oplus\cdots\oplus\Z/d_s\Z$ と $m\ge1$ について、

$$
|A[m]|=\prod_{i=1}^s\gcd(m,d_i),\qquad |mA|=\frac{|A|}{|A[m]|}
$$

である。また $b=(b_1,\ldots,b_s)\in A$ について、方程式 $mx=b$ が解をもつ条件はすべての $i$ で $\gcd(m,d_i)\mid b_i$ が成り立つことであり、そのとき解はちょうど $|A[m]|$ 個ある。
:::

:::proof
成分ごとに考えればよい。$\Z/d\Z$ で $g=\gcd(m,d)$ とする。$mx\equiv b\pmod d$ は整数 $y$ を用いて $mx+dy=b$ と書け、Bézout の等式により解があるのは $g\mid b$ のときに限る。$b=0$ のとき、$mx\equiv0\pmod d$ は $(d/g)\mid x$ と同値なので解は $g$ 個である。解が一つ $x_0$ あれば、解全体は $x_0+(\Z/d\Z)[m]$ である。直和では $A[m]=\bigoplus(\Z/d_i\Z)[m]$ であり、$m$ 倍写像 $A\to mA$ に :ref[thm-grp-first-iso] を適用すれば $|mA|=|A|/|A[m]|$ を得る。
:::

例えば $A=C_6\oplus C_{12}$ で $4x=(\bar2,\bar8)$ を解く。第一成分は $4u\equiv2\pmod6$、すなわち $u\equiv2\pmod3$ で $u=2,5$、第二成分は $4v\equiv8\pmod{12}$、すなわち $v\equiv2\pmod3$ で $v=2,5,8,11$ である。解は $2\cdot4=8=\gcd(4,6)\gcd(4,12)$ 個である。右辺を $(\bar2,\bar6)$ に替えると、第二成分で $\gcd(4,12)=4\nmid6$ なので解はない。

特に $m=2$ とすると、有限 Abel 群の位数 $2$ の元の個数は $2^t-1$ である。ここで $t$ は偶数の不変因子の個数である。

有限 Abel $p$ 群 $P$ では、位数がちょうど $p^j$ の元は $|P[p^j]|-|P[p^{j-1}]|$ 個ある。位数 $p^j$ の巡回部分群は $\varphi(p^j)=p^{j-1}(p-1)$ 個の生成元をもち、各元はただ一つの巡回部分群を生成するので、位数 $p^j$ の巡回部分群の個数はこれを $\varphi(p^j)$ で割ったものである。例えば $C_4\oplus C_4$ では $|P[2]|=4$、$|P[4]|=16$ なので、位数 $4$ の元は $12$ 個、位数 $4$ の巡回部分群は $12/2=6$ 個、位数 $2$ の部分群は $3$ 個ある。

:::proposition[有限 Abel 群の部分群の位数]{#prop-grp-abelian-subgroup-orders}
有限 Abel 群 $A$ と $|A|$ の任意の正の約数 $m$ に対し、位数 $m$ の部分群が存在する。
:::

:::proof
:ref[lem-grp-primary] により $A=\bigoplus_pA_p$ であり、$m=\prod_pp^{b_p}$ とすると $p^{b_p}\le|A_p|$ である。各 $A_p$ に位数 $p^{b_p}$ の部分群 $B_p$ があれば、$\bigoplus_pB_p$ は位数 $m$ である。そこで $P=A_p\cong\bigoplus_{i=1}^t\Z/p^{e_i}\Z$、$b\le\sum_ie_i$ とする。$0\le b_i\le e_i$、$\sum_ib_i=b$ を満たす $b_i$ を選び（先頭から順に可能な限り大きくとればよい）、$\Z/p^{e_i}\Z$ の部分群 $\langle\overline{p^{e_i-b_i}}\rangle$（位数 $p^{b_i}$、:ref[thm-grp-cyclic-subgroups]）の直和をとれば、位数 $p^b$ の部分群を得る。
:::

非可換群ではこの主張は成り立たない。位数 $12$ の交代群 $A_4$ は位数 $6$ の部分群をもたない（[剰余類と Lagrange の定理](/math/group-theory/cosets/)、[対称群](/math/group-theory/symmetric-groups/)）。

### 生成系と極大部分群

有限 Abel $p$ 群 $P$ では $P/pP$ の元は $p$ 倍で零になるので、整数 $c$ の類 $\bar c\in\F_p$ によるスカラー倍 $\bar c(x+pP)=cx+pP$ が代表のとり方によらず定まり、$P/pP$ は $\F_p$ 上のベクトル空間になる。$P\cong\bigoplus_{i=1}^t\Z/p^{e_i}\Z$（$e_i\ge1$）なら、各因子で $(\Z/p^e\Z)/p(\Z/p^e\Z)\cong\Z/p\Z$ なので（:ref[prop-grp-product-center] の 2 を繰り返し用いて）$P/pP\cong\F_p^t$ である。従って $\dim_{\F_p}P/pP$ は巡回因子の個数 $t$ に等しい。

:::proposition[法 p での生成判定]{#prop-grp-abelian-generators}
$P$ を有限 Abel $p$ 群、$S\subseteq P$ とする。$S$ が $P$ を生成するための必要十分条件は、$S$ の像が $P/pP$ を $\F_p$ ベクトル空間として生成することである。特に $P$ を生成する元の最小個数は $t=\dim_{\F_p}P/pP$ である。
:::

:::proof
$S$ が $P$ を生成すれば、像は $P/pP$ を整数係数の和で生成し、従ってベクトル空間としても生成する。逆に $H=\langle S\rangle$ とすると、像が $P/pP$ を生成することは $P=H+pP$ を意味する。$P$ の全元を消す最小の $p$ 冪を $p^e$ とする。任意の $x$ を $x=h_0+px_1$（$h_0\in H$）と書き、$x_1$ についても $x_1=h_1+px_2$ と書くことを繰り返すと

$$
x=h_0+ph_1+\cdots+p^{e-1}h_{e-1}+p^ex_e
$$

となり、末項は $0$ だから $x\in H$ である。よって $H=P$ である。最後の主張は、$t$ 次元空間の生成系は $t$ 個以上の元を含み、標準的な分解の生成元 $t$ 個が実際に $P$ を生成することから従う。
:::

:::corollary[有限生成 Abel 群の最小生成元数]{#cor-grp-min-generators}
$A\cong\Z^r\oplus\Z/d_1\Z\oplus\cdots\oplus\Z/d_s\Z$（$1<d_1\mid\cdots\mid d_s$）を生成する元の最小個数は $r+s$ である。
:::

:::proof
標準的な生成元 $r+s$ 個で生成される。逆向きの不等式を示す。$s\ge1$ なら $d_1$ を割る素数 $p$ を、$s=0$ なら $p=2$ をとる。$p$ はすべての $d_i$ を割るので $(\Z/d_i\Z)/p(\Z/d_i\Z)\cong\Z/p\Z$、自由部分の各座標も法 $p$ で一次元になるから、$A/pA\cong\F_p^{r+s}$ である。$A$ の生成系の像は $A/pA$ を生成するので、$r+s$ 個以上の元を含む。
:::

例えば $\Z\oplus C_6\oplus C_{12}$ を生成する元の最小個数は $3$ である。$C_6\oplus C_{10}\cong C_2\oplus C_{30}$ は二元で生成されるが一元では生成されない。

部分群 $M\subsetneq P$ で、$M\subsetneq L\subsetneq P$ となる部分群 $L$ が存在しないものを**極大部分群**という。

:::proposition[極大部分群の個数]{#prop-grp-abelian-maximal}
有限 Abel $p$ 群 $P$ の極大部分群は $pP$ を含む指数 $p$ の部分群であり、$P/pP$ の超平面（余次元 $1$ の部分空間）と一対一に対応する。その個数は $t=\dim_{\F_p}P/pP$ として $(p^t-1)/(p-1)$ である。
:::

:::proof
$M$ を極大部分群とし、$x\notin M$ をとる。$M+\langle x\rangle$ は $M$ を真に含むので $P$ に等しく、$P/M$ は $x+M$ で生成される巡回群である。その位数 $n$ が合成数 $n=ab$（$1<a<n$）なら、$\langle a(x+M)\rangle$ は位数 $b$ の非自明な真部分群であり、その逆像は $M$ と $P$ の真に間にある部分群となって矛盾する。$n$ は $|P|$ を割るので $p$ の冪であり、従って $n=p$ である。すると $p(P/M)=0$、すなわち $pP\subseteq M$ である。

$pP$ を含む部分群 $M$ と $P/pP$ の部分群 $M/pP$ は、逆像をとる操作によって一対一に対応し、包含関係も保たれる（[準同型と同型定理](/math/group-theory/homomorphisms/)の対応定理）。$P/pP$ の部分群は整数倍で閉じているので部分空間であり、極大部分群は極大な真部分空間、すなわち超平面に対応する。超平面は零でない線形形式 $\lambda:\F_p^t\to\F_p$ の核であり、二つの形式が同じ核をもつのは一方が他方の零でないスカラー倍であるときに限る（核で割った一次元空間上で値が比例する）。零でない形式は $p^t-1$ 個あり、同じ核をもつものが $p-1$ 個ずつあるので、超平面は $(p^t-1)/(p-1)$ 個である。
:::

例えば $C_4\oplus C_2$ では $t=2$、$p=2$ なので極大部分群は $3$ 個である。二つの座標の偶奇を $u,v\in\F_2$ と書けば、それぞれ $u=0$、$v=0$、$u+v=0$ で定まる部分群（位数 $4$）である。

### 体の乗法群の有限部分群

基本定理の応用として、体の乗法群の有限部分群が巡回群であることを示す。[群の定義と例](/math/group-theory/groups/)では群の指数を使って同じ結果を示したが、ここでは不変因子から直接読み取る。

:::lemma[多項式の根の個数]{#lem-grp-root-count}
$K$ を体とし、$f\in K[t]$ を次数 $n\ge0$ の零でない多項式とする。$f$ の $K$ における根は高々 $n$ 個である。
:::

:::proof
$n$ に関する帰納法で示す。$n=0$ なら $f$ は零でない定数で根をもたない。$n\ge1$ とし、$f=\sum_{k=0}^na_kt^k$ が根 $c$ をもつとする。$t^k-c^k=(t-c)(t^{k-1}+t^{k-2}c+\cdots+c^{k-1})$ を用いると

$$
f(t)=f(t)-f(c)=\sum_{k=1}^na_k(t^k-c^k)=(t-c)\,g(t)
$$

と書け、$g$ は次数 $n-1$（最高次係数 $a_n\ne0$）の多項式である。$c'\ne c$ も根なら $0=(c'-c)g(c')$ であり、体では $c'-c\ne0$ で割れるので $g(c')=0$ である。よって $f$ の根は $c$ と $g$ の根に限られ、帰納法の仮定により高々 $1+(n-1)=n$ 個である。
:::

:::theorem[体の乗法群の有限部分群]{#thm-grp-finite-mult-cyclic}
体 $K$ の乗法群 $K^\times$ の有限部分群は巡回群である。
:::

:::proof
$U=\{1\}$ は巡回群である。$U\ne\{1\}$ のとき不変因子を $1<d_1\mid\cdots\mid d_s$ とすると、すべての元は $d_s$ 乗して $1$ になる。:ref[lem-grp-root-count] を $t^{d_s}-1$ に用いれば $d_1\cdots d_s=|U|\le d_s$ であり、$s=1$ が従う。
:::

:::corollary[原始根の存在]{#cor-grp-primitive-root}
素数 $p$ に対し $(\Z/p\Z)^\times=\F_p^\times$ は位数 $p-1$ の巡回群である。その生成元を法 $p$ の :term[原始根]{en="primitive root" yomi="げんしこん"} といい、原始根はちょうど $\varphi(p-1)$ 個ある。
:::

:::proof
$\F_p=\Z/p\Z$ は体なので :ref[thm-grp-finite-mult-cyclic] を適用すればよい。位数 $p-1$ の巡回群の生成元は $\varphi(p-1)$ 個である（[群の定義と例](/math/group-theory/groups/)の法 $18$ の例）。
:::

例えば法 $7$ の原始根は $3,5$ の $\varphi(6)=2$ 個、法 $13$ の原始根は $2,6,7,11$ の $\varphi(12)=4$ 個である。一方、合成数を法とする単数群は巡回群とは限らない。$(\Z/8\Z)^\times\cong C_2\times C_2$、$(\Z/15\Z)^\times\cong C_2\times C_4$ はいずれも巡回群でない。ただし $(\Z/9\Z)^\times$ は $2$ で生成される巡回群である。

### 生成元と関係式から群を決定する

生成元と関係式で与えられた Abel 群は、次の手順で基本定理の形に書き直せる。

1. 関係式の係数を行とする関係行列 $M$ を作る。群は $\Z^n/R(M)$ である。
2. 整数基本変形で Smith 標準形 $\diag(d_1,\ldots,d_k)$ に変形する（:ref[thm-grp-smith]）。
3. :ref[lem-grp-relation-change] により、群は $\Z/d_1\Z\oplus\cdots\oplus\Z/d_k\Z\oplus\Z^{n-k}$ と同型である。$d_i=1$ の因子は除く。

:::example[二つの生成元と二つの関係式]{#exm-grp-relations-2x2}
生成元 $x,y$ に関係式 $6x+4y=0$、$4x+6y=0$ だけを課した Abel 群 $A$ を求める。「だけ」とは、$\Z^2$ を関係ベクトル $(6,4),(4,6)$ の生成する部分群で割るという意味である。関係行列に次の変形を施す。

$$
\begin{pmatrix}6&4\\4&6\end{pmatrix}
\xrightarrow{R_1-R_2}
\begin{pmatrix}2&-2\\4&6\end{pmatrix}
\xrightarrow{R_2-2R_1}
\begin{pmatrix}2&-2\\0&10\end{pmatrix}
\xrightarrow{C_2+C_1}
\begin{pmatrix}2&0\\0&10\end{pmatrix}.
$$

従って $A\cong C_2\oplus C_{10}$ であり、位数は $|\det M|=20$ に一致する。列の変形 $C_2+C_1$ は生成元の取り替えに当たる。$Q=\begin{pmatrix}1&1\\0&1\end{pmatrix}$ について $Q^{-1}$ の行は $(1,-1),(0,1)$ なので、新しい生成元は $u=x-y$、$v=y$ である。実際、$x=u+v$ を代入すると関係式は $6u+10v=0$、$4u+10v=0$ となり、その差から $2u=0$、さらに $10v=0$ を得る。逆に $2u=10v=0$ なら元の二式が成り立つ。
:::

途中で「$2u=0$ だから $u$ の位数は $2$」と断定してはいけない。残りの関係式が $u=0$ を強制する可能性があるからである。この例で位数がちょうど $2$ と分かるのは、関係の全体が新しい座標で $2\Z\oplus10\Z$ と完全に記述できた段階である。また、実数の連立方程式のように関係式を $2$ で割ると、有限位数の因子を失ってしまう。

:::example[三つの生成元の例]{#exm-grp-relations-3x3}
生成元 $x,y,z$ と関係式

$$
2x+4y+4z=0,\qquad -6x+6y+12z=0,\qquad 10x-4y-16z=0
$$

で定まる群を求める。$(1,1)$ 成分 $2$ はすべての成分を割るので、第 $2$、$3$ 列から第 $1$ 列の $2$ 倍を引き、第 $2$ 行に第 $1$ 行の $3$ 倍を足し、第 $3$ 行から第 $1$ 行の $5$ 倍を引くと

$$
\begin{pmatrix}2&4&4\\-6&6&12\\10&-4&-16\end{pmatrix}
\longrightarrow
\begin{pmatrix}2&0&0\\0&18&24\\0&-24&-36\end{pmatrix}
$$

となる。右下の $2\times2$ 部分の成分はすべて $2$ の倍数である。これをさらに

$$
\begin{pmatrix}18&24\\-24&-36\end{pmatrix}
\xrightarrow{C_2-C_1}
\begin{pmatrix}18&6\\-24&-12\end{pmatrix}
\xrightarrow{C_1\leftrightarrow C_2}
\begin{pmatrix}6&18\\-12&-24\end{pmatrix}
\xrightarrow{C_2-3C_1}
\begin{pmatrix}6&0\\-12&12\end{pmatrix}
\xrightarrow{R_2+2R_1}
\begin{pmatrix}6&0\\0&12\end{pmatrix}
$$

と変形する。Smith 標準形は $\diag(2,6,12)$ であり、群は $C_2\oplus C_6\oplus C_{12}$、位数 $144=|\det M|$ である。単因子は $2,2,4,3,3$ である。検算として、全成分の最大公約数 $2=d_1$、$2$ 次小行列式全体の最大公約数 $12=d_1d_2$ も確かめられる。
:::

:::example[自由部分が残る例]{#exm-grp-relations-free}
生成元 $x,y,z$ に関係式 $x+2y+3z=0$、$4x+5y+6z=0$ を課す。

$$
\begin{pmatrix}1&2&3\\4&5&6\end{pmatrix}
\xrightarrow{C_2-2C_1,\ C_3-3C_1}
\begin{pmatrix}1&0&0\\4&-3&-6\end{pmatrix}
\xrightarrow{R_2-4R_1}
\begin{pmatrix}1&0&0\\0&-3&-6\end{pmatrix}
\xrightarrow{C_3-2C_2,\ -R_2}
\begin{pmatrix}1&0&0\\0&3&0\end{pmatrix}
$$

より、群は $\Z/1\Z\oplus\Z/3\Z\oplus\Z\cong C_3\oplus\Z$ であり、階数は $1$ である。関係式が一つ $2x+4y-6z=0$ だけなら、$u=x+2y-3z,\ v=y,\ w=z$ が $\Z^3$ の新しい基底であり（$x=u-2v+3w$ と戻せる）、関係式は $2u=0$ となるので、群は $C_2\oplus\Z^2$ である。
:::

## 演習問題

::::exercise{level=basic}
$C_{12}\times C_{18}$ の元 $(\bar4,\bar3)$ の位数を求めよ。

:::solution
:ref[thm-grp-cyclic-order] により、$C_{12}$ での $\bar4$ の位数は $12/\gcd(12,4)=3$、$C_{18}$ での $\bar3$ の位数は $18/\gcd(18,3)=6$ である。:ref[prop-grp-product-order] より、組の位数は $\operatorname{lcm}(3,6)=6$ である。
:::
::::

::::exercise{level=basic}
$C_{12}\oplus C_{18}\oplus C_{10}$ の単因子と不変因子を求めよ。

:::solution
$12=4\cdot3$、$18=2\cdot9$、$10=2\cdot5$ と分解すると（:ref[exm-grp-crt]）、単因子は $2$ 部分が $2,2,4$、$3$ 部分が $3,9$、$5$ 部分が $5$ である。各素数の指数を大きい順に並べて同順位を掛けると、最大の不変因子は $4\cdot9\cdot5=180$、次は $2\cdot3=6$、最小は $2$ である。従って不変因子は $(2,6,180)$ であり、$C_{12}\oplus C_{18}\oplus C_{10}\cong C_2\oplus C_6\oplus C_{180}$ となる。位数 $2\cdot6\cdot180=2160=12\cdot18\cdot10$ も一致する。
:::
::::

::::exercise{level=basic}
位数 $360$ の Abel 群を同型を除いてすべて求め、不変因子で表せ。

:::solution
$360=2^3\cdot3^2\cdot5$ なので、:ref[prop-grp-abelian-count] により $\mathrm{p}(3)\mathrm{p}(2)\mathrm{p}(1)=3\cdot2\cdot1=6$ 個ある。$2$ 部分の指数は $(3),(2,1),(1,1,1)$、$3$ 部分は $(2),(1,1)$、$5$ 部分は $(1)$ である。同順位を掛けると

| $2$ 部分 | $3$ 部分 $(2)$ | $3$ 部分 $(1,1)$ |
| --- | --- | --- |
| $(3)$ | $C_{360}$ | $C_3\oplus C_{120}$ |
| $(2,1)$ | $C_2\oplus C_{180}$ | $C_6\oplus C_{60}$ |
| $(1,1,1)$ | $C_2\oplus C_2\oplus C_{90}$ | $C_2\oplus C_6\oplus C_{30}$ |

を得る。:ref[thm-grp-fg-abelian] の一意性から、これらは互いに同型でない。
:::
::::

::::exercise{level=basic}
$\Aut(C_{12})$ と $\Aut(C_{16})$ の構造を決定せよ。

:::solution
:ref[prop-grp-aut-cyclic] により $\Aut(C_{12})\cong(\Z/12\Z)^\times=\{1,5,7,11\}$ であり、:ref[prop-grp-units-crt] から $(\Z/4\Z)^\times\times(\Z/3\Z)^\times\cong C_2\times C_2$ である。

$\Aut(C_{16})\cong(\Z/16\Z)^\times=\{1,3,5,7,9,11,13,15\}$ は位数 $8$ の Abel 群である。$3$ の冪は $3,9,11,1$ なので位数 $4$ の元があり、$7^2=49\equiv1$、$9^2=81\equiv1$、$15^2\equiv1$ なので位数 $2$ の元は $7,9,15$ の $3$ 個である。位数 $8$ の Abel 群は $C_8$、$C_2\oplus C_4$、$C_2^3$ のいずれかで、位数 $2$ の元の個数はそれぞれ $1,3,7$ である（:ref[prop-grp-abelian-kernel-count] のあとの注意）。従って $\Aut(C_{16})\cong C_2\times C_4$ である。
:::
::::

::::exercise{level=basic}
半直積 $C_4\rtimes_\alpha C_2$ を同型を除いてすべて求めよ。

:::solution
$\Aut(C_4)\cong(\Z/4\Z)^\times=\{1,3\}$ は位数 $2$ なので、準同型 $\alpha:C_2\to\Aut(C_4)$ は自明なものと、生成元を $3=-1$ 倍に送るものの二つである。前者は $C_4\times C_2$、後者は :ref[exm-grp-semidirect-basic] と同様に $D_4$ を与える。前者は可換、後者は非可換なので同型でない。従って二つである。
:::
::::

::::exercise{level=basic}
$A$ を Abel 群とし、$C_2$ の生成元が $-1$ 倍で作用する半直積 $G=A\rtimes_{-1}C_2$ を考える（**一般二面体群**）。$A$ の外の元 $(a,\bar1)$ はすべて位数 $2$ であることを示し、$G$ が可換であるための条件と、$2A\ne0$ のときの中心を求めよ。

:::solution
$A$ は可換なので $a\mapsto-a$ は自己同型である。:eqref[eq-grp-semidirect-product] より $(a,\bar1)(a,\bar1)=(a+(-a),\bar0)=(0,\bar0)$ であり、$(a,\bar1)\ne e$ なので位数 $2$ である。:ref[prop-grp-semidirect-basic] の 3 から $(0,\bar1)$ による共役は $A$ 上で $-1$ 倍として働くので、$G$ が可換なら $-a=a$、すなわち $2A=0$ である。逆に $2A=0$ なら作用は自明で $G=A\times C_2$ は可換である。

$2A\ne0$ とする。$(a,\bar0)$ が $(0,\bar1)$ と可換である条件は $-a=a$ であり、$(a,\bar0)$ は $(b,\bar0)$ とは常に可換である。$(a,\bar1)(b,\bar0)=(a-b,\bar1)$、$(b,\bar0)(a,\bar1)=(a+b,\bar1)$ なので、$(a,\bar1)$ が中心に属するには全 $b$ で $2b=0$ が必要であり、これは仮定に反する。従って $Z(G)=A[2]\times\{\bar0\}$ である。$A=C_n$、$n\ge3$ なら $G\cong D_n$ であり、中心は $n$ が奇数なら自明、偶数なら $\{e,r^{n/2}\}$ である。
:::
::::

::::exercise{level=advanced}
有限 Abel $2$ 群 $P$ の位数が $32$ で、$|P[2]|=8$、$|P[4]|=16$ である。$P$ の型を決定せよ。

:::solution
$P\cong\bigoplus_i\Z/2^{e_i}\Z$ とする。:eqref[eq-grp-primary-layers] により $\#\{i:e_i\ge1\}=\log_28=3$、$\#\{i:e_i\ge2\}=\log_2(16/8)=1$ である。従って因子は $3$ 個で、そのうち $2$ 個は $e_i=1$ である。$\sum e_i=5$ なので残る因子の指数は $5-1-1=3$ である。よって $P\cong C_2\oplus C_2\oplus C_8$ である。実際この群では $|P[2]|=2\cdot2\cdot2=8$、$|P[4]|=2\cdot2\cdot4=16$ となり条件を満たす。
:::
::::

::::exercise{level=advanced}
有理数の加法群 $\Q$ は捩れのない群であるが、有限生成でなく、自由 Abel 群でもないことを示せ。

:::solution
$nq=0$（$n\ge1$）なら $q=0$ なので捩れはない。有限個の有理数 $q_1,\ldots,q_k$ の分母の公倍数を $D$ とすると、それらの整数係数の和はすべて $\frac1D\Z$ に属する。$\frac1{2D}\notin\frac1D\Z$ なので、有限個の元では $\Q$ を生成できない。

自由 Abel 群でないことを示す。零でない二つの有理数 $a/b,\ c/d$ について $(cb)\cdot\frac ab-(ad)\cdot\frac cd=0$ であり、係数 $cb,ad$ は零でない。従って基底は高々一つの元しか含めず、一元からなる基底があれば $\Q\cong\Z$ は有限生成となって矛盾する。有限生成という仮定を外すと、「捩れのない群は自由」という主張は成り立たない。
:::
::::

::::exercise{level=advanced}
正整数 $m,n$ について、準同型全体 $\Hom(C_m,C_n)$ は値ごとの和について位数 $\gcd(m,n)$ の巡回群と同型であることを示せ。

:::solution
準同型 $f:C_m\to C_n$ は $f(\bar1)$ で決まり、$h\in C_n$ が $f(\bar1)$ の値になれる条件は $mh=0$ である（[準同型と同型定理](/math/group-theory/homomorphisms/)の巡回群からの準同型の例）。$(f+g)(\bar1)=f(\bar1)+g(\bar1)$ なので、$f\mapsto f(\bar1)$ は $\Hom(C_m,C_n)$ から $C_n[m]$ への群の同型である。:ref[prop-grp-abelian-kernel-count] より $|C_n[m]|=\gcd(m,n)$ であり、$C_n[m]$ は巡回群 $C_n$ の部分群なので巡回群である（:ref[thm-grp-cyclic-subgroups]）。具体的には $\overline{n/\gcd(m,n)}$ が生成する。
:::
::::

::::exercise{level=advanced}
$A=C_2\oplus C_4$ を生成する順序対 $(x,y)$ の個数と、$|\Aut(A)|$ を求めよ。

:::solution
$2A=\{(\bar0,\bar0),(\bar0,\bar2)\}$ で $A/2A\cong\F_2^2$ である。:ref[prop-grp-abelian-generators] により、$(x,y)$ が $A$ を生成することと、像 $(\bar x,\bar y)$ が $\F_2^2$ の基底であることは同値である。順序つき基底は $3\cdot2=6$ 個あり、$A/2A$ の各元の逆像は $|2A|=2$ 個なので、生成する順序対は $6\cdot2\cdot2=24$ 個である。

自己同型 $\phi$ は $e=(\bar1,\bar0)$、$f=(\bar0,\bar1)$ の像で決まる。逆に生成対 $(x,y)$ が $2x=0$ を満たせば（$4y=0$ は常に成り立つ）、$(a,b)\mapsto ax+by$ は定義できる準同型で、全射なので有限性から全単射である。$2x=0$ を満たす $x$ は $A[2]=\{(\bar0,\bar0),(\bar1,\bar0),(\bar0,\bar2),(\bar1,\bar2)\}$ のうち $2A$ に属さない $(\bar1,\bar0),(\bar1,\bar2)$ の $2$ 個である。$\bar y$ は $\bar0,\bar x$ 以外の $2$ 類のいずれかで、各類に $2$ 元ずつあるので $y$ は $4$ 通りである。よって $|\Aut(A)|=2\cdot4=8$ である。
:::
::::

::::exercise{level=advanced}
有限生成 Abel 群 $A,B$ が $A\cong A\oplus B$ を満たすなら $B=0$ であることを示せ。また、有限生成の仮定を外すと成り立たない例を与えよ。

:::solution
$B$ は有限生成群 $A\oplus B$（$A$ と同型）の射影像なので有限生成である。$T(A\oplus B)=T(A)\oplus T(B)$ であり、$(A\oplus B)/T(A\oplus B)\cong A/T(A)\oplus B/T(B)$ なので、:ref[lem-grp-free-splitting] の階数は和になる。$A\cong A\oplus B$ の両辺の階数を比べると $B$ の階数は $0$ で、$B=T(B)$ は有限である。さらに両辺の捩れ部分群の位数を比べると $|T(A)|=|T(A)|\,|B|$ なので $|B|=1$ である。

可算無限個の直和 $A=\bigoplus_{i\ge1}C_2$ をとると、$(x_1,x_2,\ldots)\mapsto(x_2,x_3,\ldots;\ x_1)$ は $A\cong A\oplus C_2$ を与える同型である。
:::
::::

::::exercise{level=advanced}
位数 $28$ の群を同型を除いて分類せよ。

:::solution
:ref[thm-grp-sylow-count] により Sylow $7$ 部分群の個数は $4$ を割り $1\bmod7$ なので $1$ であり、$N\cong C_7$ は正規である。Sylow $2$ 部分群 $H$ は位数 $4$ で $C_4$ か $C_2\times C_2$ であり、$G\cong C_7\rtimes_\alpha H$ である（:ref[thm-grp-semidirect-internal]）。$\Aut(C_7)\cong(\Z/7\Z)^\times$ は位数 $6$ の巡回群なので、位数 $2$ の元は $-1$ だけである。

$H=C_4=\langle b\rangle$ なら $\alpha(b)$ の位数は $\gcd(4,6)=2$ を割るので $\alpha(b)=\pm1$ であり、$C_{28}$ と $C_7\rtimes_{-1}C_4$ を得る。$H=C_2\times C_2$ なら像は $\{\pm1\}$ に含まれ、自明なら $C_{14}\times C_2$、自明でなければ :ref[prop-grp-order-20] の証明と同様に $D_7\times C_2\cong D_{14}$（:ref[exm-grp-internal-direct] の 3）である。

四つは互いに同型でない。$C_{28}$ と $C_{14}\times C_2$ は可換で、前者のみ位数 $28$ の元をもつ。非可換な二つのうち $D_{14}$ の位数 $2$ の元は $15$ 個、$C_7\rtimes_{-1}C_4$ の位数 $2$ の元は $(0,\bar2)$ の $1$ 個である（$(n,k)^2=(n+(-1)^kn,2k)$）。従って位数 $28$ の群は $C_{28},\ C_{14}\times C_2,\ D_{14},\ C_7\rtimes_{-1}C_4$ の四つである。
:::
::::

::::exercise{level=advanced}
零行列でない整数行列 $M$ の全成分の最大公約数が Smith 標準形の $d_1$ に等しいことを示せ。さらに$\det M\ne0$ の $2$ 次整数行列では $d_2=|\det M|/d_1$ となることを示し、関係行列 $\begin{pmatrix}12&18\\8&30\end{pmatrix}$ で定まる Abel 群を求めよ。

:::solution
行の入れ替えと $-1$ 倍は成分の集合を符号を除いて変えない。第 $i$ 行に第 $j$ 行の $c$ 倍を加えると、新しい成分は元の成分の整数係数の和なので、元の最大公約数 $g$ の倍数であり、新しい最大公約数 $g'$ について $g\mid g'$ である。逆の変形（$-c$ 倍を加える）も整数基本変形なので $g'\mid g$、従って $g=g'$ である。列についても同様である。Smith 標準形の成分の最大公約数は、$d_1$ がすべての $d_i$ を割ることから $d_1$ である。

$2$ 次正方行列では、$P,Q$ の行列式は $\pm1$ なので $|\det(PMQ)|=|\det M|\ne0$ である。従って標準形は $\diag(d_1,d_2)$ で $d_1d_2=|\det M|$ となる。

与えられた行列では成分 $12,18,8,30$ の最大公約数は $2$、$\det=12\cdot30-18\cdot8=216$ なので $d_1=2$、$d_2=108$ である。群は $C_2\oplus C_{108}$ である。
:::
::::

::::exercise{level=advanced}
$\GL_2(\Z)$ の元 $A=\begin{pmatrix}0&-1\\1&0\end{pmatrix}$、$B=\begin{pmatrix}0&1\\-1&-1\end{pmatrix}$ はともに有限位数であるが、$AB$ は無限位数であることを示せ。従って非可換群では有限位数の元の全体は部分群をなすとは限らない。

:::solution
直接計算で $A^2=-I$ なので $A^4=I$ である。$B^2=\begin{pmatrix}-1&-1\\1&0\end{pmatrix}$、$B^3=B^2B=I$ である。一方

$$
AB=\begin{pmatrix}1&1\\0&1\end{pmatrix},\qquad (AB)^k=\begin{pmatrix}1&k\\0&1\end{pmatrix}
$$

であり（$k$ に関する帰納法）、$k\ne0$ なら $(AB)^k\ne I$ である。従って有限位数の元 $A,B$ の積が無限位数となり、有限位数の元全体は積で閉じていない。:ref[def-grp-torsion] で可換性を用いたのはこのためである。
:::
::::

::::exercise{level=advanced}
$k\ge3$ のとき $(\Z/2^k\Z)^\times$ は巡回群でないことを示せ。

:::solution
$x=2^{k-1}\pm1$ とすると $x^2=2^{2k-2}\pm2^k+1$ で、$2k-2\ge k$ だから $x^2\equiv1\pmod{2^k}$ である。$-1$ も平方が $1$ である。$k\ge3$ なら $1,\ -1,\ 2^{k-1}-1,\ 2^{k-1}+1$ は法 $2^k$ で相異なる（$3\le2^{k-1}-1<2^{k-1}+1<2^k-1$）。従って位数 $2$ の元が $3$ 個以上ある。一方、巡回群には位数 $2$ の部分群が高々一つしかなく（:ref[thm-grp-cyclic-subgroups]）、位数 $2$ の元は高々一つである。よって $(\Z/2^k\Z)^\times$ は巡回群でない。
:::
::::

## 確認テスト

::::quiz
半直積 $N\rtimes_\alpha H$ について常に正しいものはどれか。

- [x] $N\times\{e\}$ は正規部分群である
- [ ] $\{e\}\times H$ は正規部分群である
- [ ] $N\times\{e\}$ と $\{e\}\times H$ はともに中心に含まれる

:::explanation
$N\times\{e\}$ は第二成分への射影 $(n,h)\mapsto h$ の核なので常に正規である（:ref[prop-grp-semidirect-basic]）。$\{e\}\times H$ が正規であるのは作用 $\alpha$ が自明なときに限る。例えば $D_n\cong C_n\rtimes_{-1}C_2$ で鏡映の生成する部分群は正規でない。中心に含まれるという主張も、$D_3\cong S_3$ の中心が自明であることから誤りである。
:::
::::

::::quiz{answer="6" tol="0"}
位数 $72$ の Abel 群の同型類の個数はいくらか。

:::explanation
$72=2^3\cdot3^2$ である。:ref[prop-grp-abelian-count] により、同型類は $3$ の分割（$3,\ 2+1,\ 1+1+1$ の三通り）と $2$ の分割（$2,\ 1+1$ の二通り）の組に対応するので、$3\cdot2=6$ 個である。具体的な一覧は :ref[exm-grp-order-72] にある。
:::
::::

::::quiz
$C_4\oplus C_6$ の不変因子はどれか。

- [ ] $(4,6)$
- [x] $(2,12)$
- [ ] $(24)$

:::explanation
$C_6\cong C_2\oplus C_3$ なので単因子は $2,4,3$ である。各素数の指数を大きい順に並べて同順位を掛けると $4\cdot3=12$ と $2$ になり、不変因子は $(2,12)$ である。$(4,6)$ は $4\nmid6$ なので整除の列でなく、$C_{24}$ は位数 $24$ の元をもつが $C_4\oplus C_6$ の元の位数は $\operatorname{lcm}(4,6)=12$ 以下である。
:::
::::

::::quiz{answer="5" tol="0"}
位数 $20$ の群は同型を除いていくつあるか。

:::explanation
:ref[prop-grp-order-20] により、$C_{20}$、$C_{10}\times C_2$、$D_{10}$、$C_5\rtimes_{-1}C_4$、$\mathrm{Aff}(\F_5)$ の $5$ 個である。Sylow $5$ 部分群が正規なので群は $C_5\rtimes H$（$|H|=4$）の形であり、$H=C_4$ からは作用の像 $\{1\},\{\pm1\},\Aut(C_5)$ に応じて三つ、$H=C_2\times C_2$ からは作用が自明か否かに応じて二つの群が得られ、それらは可換性・位数 $2$ の元の個数で区別される。
:::
::::

::::quiz
有限生成 Abel 群について、同型で保たれる量をすべて選べ。

- [x] 自由階数 $r$
- [x] 各素数 $p$ と正整数 $j$ に対する方程式 $p^jx=0$ の解の個数
- [ ] 選んだ生成系に含まれる元の個数
- [x] 群を生成する元の最小個数

:::explanation
自由階数は $A/T(A)$ の階数として群から定まる（:ref[lem-grp-free-splitting]）。$p^jx=0$ の解全体 $A[p^j]$ は同型によって全単射に移る。生成元の最小個数も同型で保たれ、:ref[cor-grp-min-generators] により $r+s$ に等しい。一方、生成系には既存の生成元の和などをいくらでも追加できるので、任意に選んだ生成系の大きさは不変量ではない。
:::
::::

::::quiz{answer="10" tol="0"}
関係行列 $\begin{pmatrix}6&4\\4&6\end{pmatrix}$ の Smith 標準形を $\diag(d_1,d_2)$（$d_1\mid d_2$）とするとき、$d_2$ はいくらか。

:::explanation
:ref[exm-grp-relations-2x2] の変形により標準形は $\diag(2,10)$ である。検算として、$d_1$ は全成分の最大公約数 $2$、$d_1d_2$ は $|\det M|=36-16=20$ なので $d_2=10$ である。対応する群は $C_2\oplus C_{10}$ である。
:::
::::

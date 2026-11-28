---
title: "BMO と Hardy 空間"
description: "平均振動と原子を使って、特異積分の端点評価と H¹–BMO 双対性を証明する。"
notation:
  - symbol: 'b_Q'
    meaning: "立方体 Q 上の平均"
  - symbol: '\|b\|_{\mathrm{BMO}}'
    meaning: "平均振動の上限"
  - symbol: '\mathcal M f'
    meaning: "試験関数全体で測る極大関数"
  - symbol: 'H^1(\R^n)'
    meaning: "実変数 Hardy 空間"
  - symbol: 'a_j'
    meaning: "平均ゼロの原子"
---

## この章の目標

- BMO の平均振動と John–Nirenberg の指数型評価を証明する。
- 特異積分を $L^\infty$ から BMO へ延長する。
- Hardy 空間を極大関数で定め、原子分解を構成する。
- 局所的な $L^2$ 双対から $H^1$ と BMO の双対性を導く。

## BMO 空間

点の値を一様に抑えられない関数でも、各立方体で平均値からの揺れだけなら制御できる。立方体ごとに平均値を選び直せるため、積分で生じる加法定数を扱いやすくなる。

:::definition[有界平均振動]{#def-four-bmo}
$b\in L^1_{\mathrm{loc}}(\R^n)$、$b_Q=|Q|^{-1}\int_Q b$ とする。座標軸に平行な全立方体について

$$
\|b\|_{\mathrm{BMO}}=\sup_Q |Q|^{-1}\int_Q |b-b_Q|
$$

が有限なら $b$ は :term[有界平均振動]{en="bounded mean oscillation" yomi="ゆうかいへいきんしんどう"} を持つという。この半ノルムで定数を同一視した空間を $\mathrm{BMO}$ と書く。
:::

:::proposition[別の中心で測る振動]{#prop-four-bmo-centers}
任意の定数 $c$ に対し $\int_Q|b-b_Q|\le2\int_Q|b-c|$ である。半ノルムがゼロの関数は a.e. 定数である。
:::

:::proof
$|b_Q-c|\le |Q|^{-1}\int_Q|b-c|$ と三角不等式で最初の主張を得る。半ノルムがゼロなら各 $Q$ 上で a.e. 定数であり、空間を覆う増加立方体列の重なりでその定数は全て同じである。
:::

:::example[跳躍を持つ有界関数]
$n=1$ で $b=1_{(0,\infty)}$ とする。区間の正側の長さの割合を $\theta$ とすると、振動の平均は $2\theta(1-\theta)$ である。従って BMO 半ノルムは $1/2$ で、原点を中心とする区間で達成される。
:::

:::example[有界でない BMO 関数]
$b(x)=\log|x|$ とする。立方体 $Q$ の中心から原点までの距離が $2\sqrt n\,\ell(Q)$ より大きければ、$\nabla b=x/|x|^2$ と線分での積分から $Q$ 内の振動は次元だけの定数で抑えられる。中心がそれより近ければ $c=\log\ell(Q)$ を取り、$x=\ell(Q)z$ と置くと $Q$ の振動積分は、尺度を除いて固定球での $|\log|z||$ の積分以下である。これは原点近くでも球座標で有限である。:ref[prop-four-bmo-centers] により BMO である。
:::

### 加法定数をそろえた完備性

BMO のノルムは代表元の定数部分を測らない。収束を調べるには、一つの立方体の平均をゼロに選び、そこを基準にして他の立方体上の平均を比べる。

:::proposition[BMO の完備性]{#prop-four-bmo-completeness}
定数を同一視した $\mathrm{BMO}(\R^n)$ はノルム $\|\cdot\|_{\mathrm{BMO}}$ に関して完備である。
:::

:::proof
$Q_0$ を固定し、各代表元の平均を $Q_0$ 上でゼロにする。$Q_0\subset R$ なら $|b_R|=|b_R-b_{Q_0}|\le |Q_0|^{-1}\int_{Q_0}|b-b_R|\le(|R|/|Q_0|)\|b\|_{\mathrm{BMO}}$ である。従って $R$ 上の $L^1$ ノルムは $C_R\|b\|_{\mathrm{BMO}}$ 以下である。BMO の Cauchy 列は、$Q_0$ を含む増加立方体列のそれぞれで $L^1$ の極限を持つ。極限は包含上で一致するので $b\in L^1_{\rm loc}$ を定める。任意の $Q$ 上で、この局所 $L^1$ 収束は平均値と平均振動の収束を与える。Cauchy 条件 $\|(b_j-b_k)\|_{\mathrm{BMO}}\le\eps$ に対し $k\to\infty$ として各 $Q$ 上の振動を $\eps$ 以下とし、次に $Q$ の上限を取れば $\|b_j-b\|_{\mathrm{BMO}}\le\eps$ を得る。
:::

### 立方体と球で測った振動

平均振動の定義には座標軸に平行な立方体を使った。核を球殻ごとに積分する場合は球を使う方が簡単なので、ここで二種類の測り方を比べる。包含する集合の体積比が一定なら、その二つの平均の差も制御できる。

:::lemma[入れ子の集合における平均の比較]{#lem-four-bmo-nested-means}
$E\subset Q$ が正の測度を持つ可測集合なら
$$
|b_E-b_Q|\le\frac{|Q|}{|E|}\|b\|_{\mathrm{BMO}},
\qquad
\frac1{|E|}\int_E|b-b_E|\le2\frac{|Q|}{|E|}\|b\|_{\mathrm{BMO}}.
$$
立方体の代わりに球全体を用いて定義した半ノルムは $\|\cdot\|_{\mathrm{BMO}}$ と同値である。
:::

:::proof
第 1 式は $b_E-b_Q$ を $E$ 上の $b-b_Q$ の平均として評価すれば得られる。第 2 式は $|\int_E(b-b_Q)|\le\int_E|b-b_Q|$ と三角不等式によって得られる。同じ証明は $Q$ を球としても成り立つ。各球を外接する座標軸に平行な立方体に含めると体積比は次元だけに依存する。逆に各立方体を外接球に含めても体積比は次元だけに依存する。両方に第 2 式を使えば半ノルムが比較できる。
:::

:::example[同じ有界関数でも立方体の位置で平均が変わる]
$b=1_{[0,1]}$ とする。$Q=[-1,1]$ では $b_Q=1/2$、平均振動は $1/2$ である。$R=[-2,2]$ では $b_R=1/4$ なので、平均振動は $(1/4)(3/4)+(3/4)(1/4)=3/8$ になる。一方、任意の区間で $b=1$ の割合 $\theta$ に対し平均振動は $2\theta(1-\theta)\le1/2$ である。従って $Q$ の値が半ノルムを達成する。大きい区間に移しても振動が必ず増えるわけではない。
:::

:::example[半ノルムが原点の移動や拡大に依存しない理由]
$g(x)=b(a+rx)$、$r>0$ と置くと $g_Q=b_{a+rQ}$ であり、振動積分でも Jacobian の $r^n$ が体積で相殺される。$Q$ を動かすと $a+rQ$ が全立方体を動くので、$\|g\|_{\mathrm{BMO}}=\|b\|_{\mathrm{BMO}}$ である。例えば $b(x)=\log|x|$ に対し $g(x)=\log|2+4x|$ は $\log4+\log|x+1/2|$ と分けられる。共通の半ノルムは対数の原点の位置にも、加わった定数にも依存しない。
:::

### 一つの二進格子と複数の格子

John–Nirenberg の証明では、一つの立方体の内部を二進分割する。一方、空間全体で一つの二進格子だけに立方体を制限すると、その格子の境界をまたぐ振動を捉えられないことがある。

:::definition[二進 BMO]{#def-four-dyadic-bmo}
各尺度 $2^k$ の立方体が分割をなし、異なる尺度で交わる立方体が入れ子となる格子 $\mathcal D$ を取る。$\sup_{Q\in\mathcal D}\operatorname{avg}_Q |b-b_Q|$ を $\|b\|_{\mathrm{BMO}(\mathcal D)}$ と書く。ここで $\operatorname{avg}_Q$ は $|Q|^{-1}\int_Q$ である。
:::

:::example[格子の境界に隠れる振動]
$\R$ の通常の格子 $[2^k j,2^k(j+1))$ では、どの区間も正と負の両側を同時には含まない。$b(x)=1_{(0,\infty)}(x)\log|x|$ とすると、正側の格子区間では対数の BMO 評価を使え、負側では振動がゼロである。従って二進 BMO 半ノルムは有限である。ところが $Q=(-r,r)$ では負側の平均は $0$、正側の平均は $\log r-1$ である。全区間の平均は $(\log r-1)/2$ なので、負側だけの寄与から $\operatorname{avg}_Q|b-b_Q|\ge|\log r-1|/4$ を得る。$r\downarrow0$ とすれば全立方体での半ノルムは無限大である。
:::

:::proposition[隣接する二進格子による判定]{#prop-four-bmo-adjacent-grids}
$\R^n$ には $3^n$ 個の二進格子 $\mathcal D_\alpha$ があり
$$
\max_\alpha\|b\|_{\mathrm{BMO}(\mathcal D_\alpha)}
\le \|b\|_{\mathrm{BMO}}
\le 2\cdot6^n\max_\alpha\|b\|_{\mathrm{BMO}(\mathcal D_\alpha)}
$$
が成り立つ。
:::

:::proof
格子には :ref[prop-four-adjacent-dyadic-grids] の $\mathcal D^\alpha$ を用いる。辺長 $\ell$ の立方体 $Q$ に対し $3\ell<2^k\le6\ell$ と取る。同命題の証明を各座標への長さ $\ell$ の射影に適用すると、$Q$ を含む立方体 $R\in\mathcal D^\alpha$ が得られ、$|R|/|Q|\le6^n$ である。$b_R$ を中心に振動を測り直せば
$$
\operatorname{avg}_Q|b-b_Q|\le2\operatorname{avg}_Q|b-b_R|
\le2\frac{|R|}{|Q|}\operatorname{avg}_R|b-b_R|.
$$
これが上側を与える。下側は各格子を全立方体の一部として考えれば従う。
:::

:::example[原点をまたぐ区間を別の格子に入れる]
$Q=[-1,1]$ では $\ell=2$ なので $2^k=8$ を選べる。$\alpha=1/3$ とすると $k=3$ で格子は $8(\Z-1/3)$ から始まる。$Q$ は $[-8/3,16/3)$ に含まれ、長さの比は $4$ である。通常の格子では原点が全ての尺度の境界なので、この区間を含む格子区間は存在しない。
:::

### 平均の増加と緩増加性

BMO で点の大きさを直接抑えることはできない。ただし原点を中心とする立方体を順に広げると、平均値の増加は尺度の対数までに抑えられる。この評価は緩増加超関数として Fourier 変換するときにも使える。

:::proposition[BMO 代表元の重み付き可積分性]{#prop-four-bmo-weight-integrable}
原点を中心とする辺長 $2$ の立方体を $Q_0$ とし、$b_{Q_0}=0$ とする。任意の $\delta>0$ に対して
$$
\int_{\R^n}|b(x)|(1+|x|)^{-n-\delta}\,dx
\le C_{n,\delta}\|b\|_{\mathrm{BMO}}
$$
である。従って全ての BMO 代表元は緩増加超関数を定める。
:::

:::proof
$Q_k=2^kQ_0$ とすると :ref[lem-four-bmo-nested-means] により $|b_{Q_{k+1}}-b_{Q_k}|\le2^n\|b\|_{\mathrm{BMO}}$ である。差を足し、$\operatorname{avg}_{Q_k}|b-b_{Q_k}|\le\|b\|_{\mathrm{BMO}}$ を使うと $\int_{Q_k}|b|\le |Q_k|(1+2^nk)\|b\|_{\mathrm{BMO}}$ となる。$Q_k\setminus Q_{k-1}$ では $1+|x|\ge 2^{k-1}$ である。この殻の重み付き積分は $C_n(1+k)2^{-k\delta}\|b\|_{\mathrm{BMO}}$ 以下であり、$k\ge1$ で可和である。内部 $Q_0$ は振動の定義で評価できる。一般の定数部分もこの重みで可積分である。積分汎関数を Schwartz 関数に作用させると、$\sup_x (1+|x|)^{n+\delta}|\phi(x)|$ の定数倍で抑えられ、連続性が従う。
:::

:::example[対数に対する尺度差の直接計算]
$\R$ で $b(x)=\log|x|$、$Q_k=(-2^k,2^k)$ とする。$\operatorname{avg}_{Q_k}\log|x|=k\log2-1$ なので隣接する平均値の差は $\log2$ である。$b_{Q_0}$ を引くと $Q_k$ 上の平均は $k\log2$ になり、上の評価に現れた $k$ の増加が実際に起きる。これを $|x|$ で測ると対数の増加に相当する。
:::

## John–Nirenberg の不等式

平均振動の制約は各立方体に課される。平均から遠い場所を小立方体に集め、そこで局所平均を引き直すと、各段で悪い集合の体積は一定比で縮む。距離が大きいほど多くの段を通るため指数減衰が得られる。

:::theorem[John–Nirenberg]{#thm-four-john-nirenberg}
次元だけに依存する正定数 $C,c$ があり、$B=\|b\|_{\mathrm{BMO}}>0$ のとき

$$
|\{x\in Q:|b(x)-b_Q|>\lambda\}|
\le C e^{-c\lambda/B}|Q|\quad(\lambda>0)
$$

である。特に $1\le q<\infty$ について $(|Q|^{-1}\int_Q|b-b_Q|^q)^{1/q}\le C_{n,q}B$ である。
:::

:::proof
$Q$ を二進分割し、$|b-b_Q|$ の平均が $2B$ を超える極大部分立方体 $Q_i$ を取る。根での平均は $B$ 以下である。:ref[thm-four-cz-decomposition] と同じ親の評価で $\sum_i|Q_i|\le |Q|/2$、$|b_{Q_i}-b_Q|\le2^{n+1}B$ を得る。補集合では Lebesgue 点において $|b-b_Q|\le2B$ である。各 $Q_i$ でも引く中心を $b_{Q_i}$ に替えて繰り返せる。第 $k$ 世代の和集合は $2^{-k}|Q|$ 以下であり、これを出た点では平均の差を世代ごとに足して $|b-b_Q|\le k\,2^{n+1}B$ である。世代を $\lfloor\lambda/(2^{n+1}B)\rfloor$ に取ると指数型評価を得る。$q$ 乗積分は $\int_0^\infty q\lambda^{q-1}$ に分布関数を掛けて積分し、$t=c\lambda/B$ に置き換えると求める評価になる。
:::

:::example[対数の尾と二乗振動]
$Q=(0,1)$ で $b(x)=-\log x$ は平均 $1$、二乗平均 $2$ を持つ。$x=e^{-t}$ と置けばそれぞれ $\int_0^\infty te^{-t}dt=1$、$\int_0^\infty t^2e^{-t}dt=2$ である。$\lambda>1$ に対し $|\{|b-b_Q|>\lambda\}|=e^{-1-\lambda}$ であり、指数型の尾が実際に現れる。
:::

:::example[指数型評価をモーメントへ換算する]
もしある $Q$ で $|Q|^{-1}|\{|b-b_Q|>\lambda\}|\le2e^{-\lambda/3}$ が得られたなら、二乗平均は $4\int_0^\infty\lambda e^{-\lambda/3}d\lambda=36$ 以下である。従って二乗平均の平方根は $6$ 以下である。これは $q=2$ の層別積分を数値で追った例である。
:::

### 指数可積分性とモーメントの定数

分布関数の指数型上界から、平均振動の高い冪を積分できる。高い冪のノルムが指数 $q$ に対してどの程度増加するかを追っておくと、指数関数の積分との対応が明らかになる。

:::corollary[指数平均と線形のモーメント増加]{#cor-four-bmo-exponential-moments}
$B=\|b\|_{\mathrm{BMO}}>0$ とする。ある $a_n>0$ に対して全ての立方体で
$$
\operatorname{avg}_Q e^{a_n|b-b_Q|/B}\le2
$$
が成り立つ。さらに全ての $q\ge1$ に対して $(\operatorname{avg}_Q|b-b_Q|^q)^{1/q}\le C_n q B$ となる。逆に、ある $D>0$ で $(\operatorname{avg}_Q|b-b_Q|^q)^{1/q}\le Dq$ が全ての正整数 $q$ と全立方体に対して成り立てば、$\operatorname{avg}_Q\exp(|b-b_Q|/(2eD))\le2$ である。
:::

:::proof
$Y=|b-b_Q|/B$ と置き、$Q$ 上の一様確率で考える。各点で $e^{aY}=1+\int_0^Y ae^{at}dt$ なので Tonelli により
$$
\operatorname{avg}_Q e^{aY}\le1+Ca\int_0^\infty e^{-(c-a)t}dt=1+\frac{Ca}{c-a}.
$$
$a_n\le c/(1+C)$ かつ $a_n<c$ に取れば $2$ 以下になる。任意の $u\ge0$ に対し $u^q e^{-a_nu}$ は $u=q/a_n$ で最大となり、最大値は $(q/(e a_n))^q$ である。従って $\operatorname{avg}_Q Y^q\le 2(q/(ea_n))^q$ から所望の上界を得る。

逆向きでは指数を冪級数に展開し、非負項を Tonelli で積分する。$k!\ge(k/e)^k$ は $\log k!\ge\int_1^k\log t\,dt=k\log k-k+1$ から従う。従って $k\ge1$ の項は $(Dk)^k/((2eD)^k k!)\le 2^{-k}$ である。定数項を含めて足すと $2$ 以下になる。
:::

:::example[対数の指数可積分性には閾値がある]
$Q=(0,1)$、$b(x)=-\log x$ では $b_Q=1$ である。$t=-\log x$ により
$$
\operatorname{avg}_Q e^{a|b-b_Q|}
=e^a\int_0^1 e^{-(a+1)t}dt+e^{-a}\int_1^\infty e^{(a-1)t}dt.
$$
第 2 項が有限なのは $a<1$ のときに限る。$a=1/2$ では第 1 項が $\frac23(e^{1/2}-e^{-1})$、第 2 項が $2e^{-1}$ なので和は $\frac23 e^{1/2}+\frac43e^{-1}$ である。また $q$ が正整数なら $\operatorname{avg}_Q|b|^q=q!$ であり、$(q!)^{1/q}\ge q/e$ が成り立つ。中心化した $b-b_Q$ も三角不等式で少なくとも $q/e-1$ なので、線形増加は一般には改良できない。
:::

### 指数を取った重みの平均

小さい実数 $t$ では $e^{tb}$ とその逆冪を同時に積分できる。平均 $b_Q$ の大きさは立方体によって変わるが、二つの平均の積を取るとその定数部分は消える。

:::definition[$A_p$ 重み]{#def-four-bmo-ap-weight}
$1<p<\infty$ とする。正の a.e. 有限な関数 $w$ について
$$
[w]_{A_p}=\sup_Q\left(\operatorname{avg}_Q w\right)
\left(\operatorname{avg}_Q w^{-1/(p-1)}\right)^{p-1}<\infty
$$
なら $w$ は :term[$A_p$ 重み]{en="Muckenhoupt weight" yomi="えーぴーおもみ"} という。各平均は Lebesgue 測度で取る。
:::

:::proposition[BMO の小さい指数は $A_p$ 重みを与える]{#prop-four-bmo-exponential-weight}
実数値の $b\in\mathrm{BMO}$ に対し $|t|\|b\|_{\mathrm{BMO}}\le a_n\min(1,p-1)$ なら、$w=e^{tb}$ は $[w]_{A_p}\le2^p$ を満たす。逆に $w\in A_p$ なら $\log w\in\mathrm{BMO}$ であり、$\|\log w\|_{\mathrm{BMO}}\le 2[w]_{A_p}$ である。
:::

:::proof
$B=0$ なら $w$ は定数で $[w]_{A_p}=1$ である。$B>0$ の場合、$e^{t b_Q}$ を第 1 平均から、$e^{-t b_Q/(p-1)}$ を第 2 平均から外に出すと両者が打ち消される。残る二つの平均は :ref[cor-four-bmo-exponential-moments] によりそれぞれ $2$ 以下である。

逆向きでは、$\log^+ w\le w$ と $\log^- w\le(p-1)w^{-1/(p-1)}$ により $b=\log w$ は局所可積分である。Jensen により $\operatorname{avg}_Q e^{b-b_Q}\ge1$、$\operatorname{avg}_Q e^{-(b-b_Q)/(p-1)}\ge1$ である。$A_p$ 条件と第 2 平均の下界から $\operatorname{avg}_Q e^{b-b_Q}\le[w]_{A_p}$ を得る。$u_+\le e^u$ と $\operatorname{avg}_Q(b-b_Q)=0$ を用いると $\operatorname{avg}_Q|b-b_Q|=2\operatorname{avg}_Q(b-b_Q)_+\le2[w]_{A_p}$ となる。
:::

:::example[べき重みと積分できる範囲]
$\R$ で $w=|x|^\alpha$ とすると、原点近くで $w$ と $w^{-1/(p-1)}$ が局所可積分となる条件は $-1<\alpha<p-1$ である。$Q=(0,r)$ の $A_p$ の積は
$$
\frac{1}{1+\alpha}\left(\frac{1}{1-\alpha/(p-1)}\right)^{p-1}
$$
で、$r$ に依存しない。この範囲では他の区間でも一様である。長さ $\ell$ に対して原点までの距離が $\ell$ 以上なら、区間内の $|x|$ は最小値の $2$ 倍以下なので二つの平均の積は $2^{|\alpha|}$ 以下に抑えられる。距離が $\ell$ 未満なら区間は $[-2\ell,2\ell]$ に含まれ、正の被積分関数の二つの平均はそれぞれこの区間の平均の $4$ 倍以下である。従って全区間で条件が成立する。これは $b=\log|x|$ による $e^{\alpha b}$ の例であり、閾値に達すれば一方の平均が無限大になる。
:::

## 特異積分の L^∞ から BMO への有界性

有界関数は無限遠でも減衰するとは限らない。そこで遠くの入力による出力の定数部分を、各立方体の中心の値で取り除く。近い入力は有限体積に台を持つので二乗評価が使える。

:::theorem[特異積分の BMO 端点]{#thm-four-cz-bmo}
$T$ が :ref[thm-four-cz] の二乗有界性と核表示を持ち、第 1 変数で $|\nabla_xK(x,y)|\le A|x-y|^{-n-1}$ とする。$T$ は $L^\infty$ から $\mathrm{BMO}$ へ、$L^\infty\cap L^2$ 上の作用と定数を除いて一致する有界な写像を定める。
:::

:::proof
$Q$ の中心を $c_Q$ とし、$Q^\dagger$ を同じ中心で辺長を $4\sqrt n$ 倍にした立方体とする。$Q$ 上の局所出力を

$$
u_Q(x)=T(f1_{Q^\dagger})(x)+\int_{(Q^\dagger)^c}[K(x,y)-K(c_Q,y)]f(y)dy
$$

と定める。近い部分の絶対平均は $|Q|^{-1/2}\|T(f1_{Q^\dagger})\|_2\le C\|f\|_\infty$ である。外側の差は $C|x-c_Q||y-c_Q|^{-n-1}$ で、積分して $C'\|f\|_\infty$ 以下である。よって :ref[prop-four-bmo-centers] で平均振動が抑えられる。

$R$ を大きくして $Q^\dagger\subset B(0,R)$ とすると、$T(f1_{B(0,R)})$ の $Q$ 上の値から $\int_{B(0,R)\setminus Q^\dagger}K(c_Q,y)f(y)dy$ を引いたものは $u_Q$ に $L^1(Q)$ で収束する。遠方の差の積分に優収束を使った。従って重なる二立方体で局所出力の差は重なり上の定数である。増加立方体列でこれらの定数を調整して全空間の代表元を作る。任意の $Q$ はその列の一つに含まれるため、平均振動は $u_Q$ で決まり、上の評価が全 $Q$ に一様に成り立つ。$f\in L^2$ なら切断した出力は元の $Tf$ に $L^2$ 収束し、局所表示との差は定数である。
:::

:::example[半直線の Hilbert 変換]
$n=1$、$f=H=1_{(0,\infty)}$ を $1_{(0,R)}$ で切断する。主値核 $1/(\pi x)$ で出力は $\pi^{-1}(\log|x|-\log|x-R|)$ である。固定区間上で $\log|x-R|-\log R\to0$ なので、定数を調整した極限は $\pi^{-1}\log|x|$ である。これは BMO に属し、原点近くでは有界でない。
:::

:::example[有界で短い入力に対する局所出力]
一次元で $f=1_{[0,1]}$ とすれば $R_1 f(x)=\pi^{-1}(\log|x|-\log|x-1|)$ である。$Q=(0,1)$ の平均は $x\mapsto1-x$ の反対称性でゼロになる。平均振動は
$$
\frac{2}{\pi}\int_0^{1/2}\log\frac{1-x}{x}dx=\frac{2\log2}{\pi}
$$
である。点 $0,1$ の近くでは有界性を失っても、立方体で平均した振動はこの有限値である。
:::

### 鋭い極大関数と局所振動からのノルム評価

平均振動の上限を取ると BMO 半ノルムになる。一方、立方体を一点を含むものに限ってから上限を取ると、場所ごとの振動が得られる。平均の大きさと振動の大きさを同時に調べれば、$L^p$ ノルムも復元できる。

:::definition[鋭い極大関数]{#def-four-bmo-sharp-maximal}
局所可積分な $u$ に対し、$M^\#u(x)=\sup_{Q\ni x}\operatorname{avg}_Q|u-u_Q|$ と置く。これを :term[鋭い極大関数]{en="sharp maximal function" yomi="するどいきょくだいかんすう"} という。二進格子に限定したものは :ref[def-four-dyadic-sharp-maximal] の $M^\#_{\mathcal D}$ である。同じ格子で $|u|$ の平均を測る $M_{\mathcal D}$ と比較する。
:::

:::lemma[振動の小さい部分における分布関数の比較]{#lem-four-sharp-good-lambda}
$u\in L^r(\R^n)$、$1\le r<\infty$ とする。$A>2^n$、$\lambda,\gamma>0$ に対し
$$
|\{M_{\mathcal D}u>A\lambda,\ M^\#_{\mathcal D}u\le\gamma\lambda\}|
\le\frac{\gamma}{A-2^n}|\{M_{\mathcal D}u>\lambda\}|
$$
である。
:::

:::proof
大きい立方体上では $\operatorname{avg}_Q|u|\le |Q|^{-1/r}\|u\|_r\to0$ なので、$\{M_{\mathcal D}u>\lambda\}$ を平均が $\lambda$ を超える極大立方体 $Q_i$ に分けられる。親の平均は $\lambda$ 以下であるから $\operatorname{avg}_{Q_i}|u|\le2^n\lambda$ である。$Q_i$ に $M^\#_{\mathcal D}u\le\gamma\lambda$ となる点が一つもなければ、その立方体は左辺へ寄与しない。そうした点が一つあれば $\int_{Q_i}|u-u_{Q_i}|\le\gamma\lambda |Q_i|$ となる。

$Q_i$ 内で $M_{\mathcal D}u>A\lambda$ となる点は、$Q_i$ の真の子孫に平均が $A\lambda$ を超える立方体を持つ。上の立方体や $Q_i$ 自身では親の評価があるためである。その子孫では $|u_{Q_i}|\le2^n\lambda$ を引いても $|u-u_{Q_i}|$ の平均が $(A-2^n)\lambda$ を超える。二進極大作用素の弱型評価は、該当する極大立方体で積分を足せば定数 $1$ で得られる。従って $Q_i$ 内の例外集合の測度は $(\gamma/(A-2^n))|Q_i|$ 以下である。互いに交わらない $Q_i$ について足すと結論が従う。
:::

:::corollary[鋭い極大関数の $L^p$ ノルム同値性]{#cor-four-sharp-maximal-lp}
$u\in L^p(\R^n)$、$1<p<\infty$ なら $\|u\|_p\asymp\|M^\#u\|_p$ である。
:::

:::proof
下界 $\|u\|_p\le C_{n,p}\|M^\#u\|_p$ は :ref[thm-four-dyadic-good-lambda] と $M^\#_{\mathcal D}u\le M^\#u$ による。上界は各立方体で $|u-u_Q|$ の平均を $|u|$ の平均の二倍で抑え、$x\in Q$ なら $Q$ を $x$ 中心の球に体積を一定倍だけ増やして含めれば従う。実際、この比較で $M^\#u\le C_nMu$ となるため、:ref[prop-four-maximal-strong] により $\|M^\#u\|_p\le C_{n,p}\|u\|_p$ である。
:::

:::example[小さい振動がある場所だけを数える]
$n=1$、$A=4$、$\gamma=1/4$ とする。平均が $\lambda$ を超える極大区間 $Q_i$ で、どこかに $M^\#_{\mathcal D}u\le\lambda/4$ の点があれば、$M_{\mathcal D}u>4\lambda$ の領域は高々 $|Q_i|/8$ である。親の平均が $\lambda$ 以下なら $|u_{Q_i}|\le2\lambda$ であり、残る $2\lambda$ は振動によって生じなければならない。平均振動が $\lambda/4$ 以下なので体積比は $(\lambda/4)/(2\lambda)=1/8$ 以下となる。
:::

### 入力を局所的な二乗平均で測る

BMO 端点の証明では、入力全体の本質上限を用いた。同じ近距離と遠距離の分割で、入力の局所平均を残しておくと点ごとの評価が得られる。$\mathfrak M_2 f(x)=\sup_{Q\ni x}(\operatorname{avg}_Q|f|^2)^{1/2}$ と置く。

:::proposition[特異積分の鋭い極大評価]{#prop-four-cz-sharp-maximal}
$T$ が :ref[thm-four-cz-bmo] の仮定を満たし、$f\in L^2(\R^n)$ とする。このとき $M^\#(Tf)(x)\le C\mathfrak M_2 f(x)$ a.e. である。
:::

:::proof
$x_0\in Q$ を固定し、前と同じ $Q^\dagger$ を取る。$T(f1_{Q^\dagger})$ の $Q$ 上の絶対平均は
$$
|Q|^{-1/2}\|T\|_{2\to2}\|f1_{Q^\dagger}\|_2
\le C\mathfrak M_2 f(x_0)
$$
である。遠い部分は、まず有界台に切って核表示を使う。$x\in Q$ とし、$2^{j+1}Q^\dagger\setminus2^jQ^\dagger$ の上で第 1 変数の差を積分すると
$$
\int_{2^{j+1}Q^\dagger\setminus2^jQ^\dagger}
|K(x,y)-K(c_Q,y)|\,|f(y)|\,dy
\le C2^{-j}\operatorname{avg}_{2^{j+1}Q^\dagger}|f|
\le C2^{-j}\mathfrak M_2 f(x_0).
$$
従って定数を引いた遠い部分は絶対収束する。切断した $f$ の $L^2$ 収束により $T$ の出力も $L^2$ 収束する。$Q$ 上の平均を引けば定数は消えるので、極限に移してよい。:ref[prop-four-bmo-centers] により $\operatorname{avg}_Q|Tf-(Tf)_Q|\le C\mathfrak M_2 f(x_0)$ を得る。$x_0$ を含む全立方体で上限を取れば結論が従う。
:::

:::example[遠い二つの尺度の寄与]
立方体 $Q$ に対し、入力が $f=1$ を $2Q^\dagger\setminus Q^\dagger$ と $8Q^\dagger\setminus4Q^\dagger$ 上で取り、他でゼロとする。$\mathfrak M_2 f\le1$ である。上の証明で非零なのは $j=0,2$ だけなので、定数を引いた遠い部分の評価は $C(1+1/4)$ となる。同じ高さの入力でも、中心から尺度で二段離れると、差分核による寄与は $1/4$ 倍に抑えられる。
:::

### Poisson 積分の勾配と箱の積分

BMO では任意の立方体上で定数を引ける。Poisson 積分を微分すると定数が消えるので、立方体の上の箱でエネルギーを積分する場合にも BMO の振動を使える。この積分では高さが小さいほど境界の局所的な変化を測る。

:::proposition[BMO の Poisson 拡張の箱評価]{#prop-four-bmo-poisson-carleson}
$b\in\mathrm{BMO}(\R^n)$ とし、$u(x,t)=P_t*b(x)$ を :ref[prop-four-upper-poisson-n] の核で定める。$u$ の積分は絶対収束し、全ての立方体 $Q$ について
$$
\int_0^{\ell(Q)}\int_Q t\,|\nabla_{x,t}u(x,t)|^2\,dx\,dt
\le C_n |Q|\,\|b\|_{\mathrm{BMO}}^2 .
$$
このように $Q\times(0,\ell(Q))$ 上の測度を $|Q|$ で一様に評価できる性質を :term[Carleson 測度]{en="Carleson measure" yomi="かーるそんそくど"} の条件という。
:::

:::proof
$t>0$ を固定すると $|P_t(x-y)|\le C_{x,t}(1+|y|)^{-n-1}$ であり、:ref[prop-four-bmo-weight-integrable] で絶対収束する。$t$ が正のコンパクト範囲にあり $x$ が有界なとき、各 1 階微分も同じ重みで抑えられるので積分内で微分できる。

$Q^\dagger$ を前と同じ十分大きい拡大とする。$b_{Q^\dagger}$ を引けば、$b_1=(b-b_{Q^\dagger})1_{Q^\dagger}$、$b_2=(b-b_{Q^\dagger})1_{(Q^\dagger)^c}$ の二つの拡張の勾配の和を評価すればよい。定数の Poisson 拡張は同じ定数となるためである。:ref[thm-four-john-nirenberg] から $\|b_1\|_2^2\le C|Q|\|b\|_{\mathrm{BMO}}^2$ である。Plancherel と $\widehat P_t(\xi)=e^{-2\pi t|\xi|}$ を使うと
$$
\begin{aligned}
\int_0^\infty\int_{\R^n}t|\nabla_{x,t}(P_t*b_1)|^2dx\,dt
&=\int_{\R^n}|\widehat b_1(\xi)|^2
   \int_0^\infty 8\pi^2 t|\xi|^2 e^{-4\pi t|\xi|}dt\,d\xi\\
&=\tfrac12\|b_1\|_2^2.
\end{aligned}
$$
$t$ 方向と空間方向の二乗がそれぞれ $4\pi^2|\xi|^2$ を与える。$\xi=0$ は積分で零測度である。

$0<t<\ell(Q)$、$x\in Q$、$y\notin Q^\dagger$ では核を直接微分して
$|\nabla_{x,t}P_t(x-y)|\le C_n|y-c_Q|^{-n-1}$ を得る。実際 $t(t^2+|x-y|^2)^{-(n+1)/2}$ の $t$ 微分の各項は $C(t^2+|x-y|^2)^{-(n+1)/2}$ 以下であり、空間微分も $t|x-y|\le(t^2+|x-y|^2)/2$ を使えば同じ評価になる。$2^{j+1}Q^\dagger$ 上の $|b-b_{Q^\dagger}|$ の平均は、隣接する平均の差を $j+1$ 回足して $C(j+2)\|b\|_{\mathrm{BMO}}$ 以下である。殻ごとの積分を足すと $|\nabla(P_t*b_2)(x)|\le C\ell(Q)^{-1}\|b\|_{\mathrm{BMO}}\sum_{j\ge0}(j+2)2^{-j}$ となる。二乗して $\int_0^{\ell(Q)}t\,dt=\ell(Q)^2/2$ を使えば遠い部分も評価できる。
:::

:::example[跳躍に対する箱のエネルギー]
$n=1$、$b=1_{(0,\infty)}$ の拡張は $u(x,t)=1/2+\pi^{-1}\arctan(x/t)$ である。微分を二乗して足すと $|\nabla u|^2=1/(\pi^2(x^2+t^2))$ となる。$Q=(-r,r)$ の箱で積分すれば
$$
\int_0^{2r}\int_{-r}^r t|\nabla u|^2\,dx\,dt
=\frac{2r}{\pi^2}\int_0^2\arctan(1/s)\,ds
=\frac{2r}{\pi^2}\left(2\arctan(1/2)+\tfrac12\log5\right).
$$
最後は $\int\arctan(1/s)ds=s\arctan(1/s)+\frac12\log(1+s^2)$ による。$r$ が変わっても $|Q|$ で割った値は同じである。
:::

## H¹ 空間と原子分解

$L^1$ では離れた場所に $|x|^{-n}$ の尾が残りやすい。平均ゼロの部品ではその先頭項が相殺される。そこで、局所化と打ち消しを同時に備えた部品を基本にする。

:::definition[実 Hardy 空間と原子]{#def-four-real-hardy}
$\mathcal A$ を $\operatorname{supp}\phi\subset B(0,1)$、$\phi\in C_c^\infty$、$\|\phi\|_\infty+\|\nabla\phi\|_\infty\le1$ を満たす試験関数の集合とし、$\phi_t(x)=t^{-n}\phi(x/t)$ とする。$f\in L^1$ について :term[大極大関数]{en="grand maximal function" yomi="だいきょくだいかんすう"} を $\mathcal M f(x)=\sup_{t>0,\phi\in\mathcal A}|\phi_t*f(x)|$ とする。$\mathcal M f\in L^1$ となる空間を実変数 Hardy 空間 $H^1(\R^n)$ といい $\|f\|_{H^1}=\|\mathcal M f\|_1$ とする。:term[原子]{en="atom" yomi="げんし"} は立方体 $Q$ に台を持ち、$\int a=0$、$\|a\|_\infty\le |Q|^{-1}$ の関数である。
:::

$\mathcal M f$ は連続な畳み込みの上限なので下半連続である。$\mathcal M f\le C_n Mf$ が絶対値による球平均から従う。また $\mathcal A$ の一つの非負関数を使って Lebesgue 点で近似すると $|f|\le C_n\mathcal M f$ である。

:::lemma[原子の極大ノルム]{#lem-four-hardy-atom-bound}
台が $Q$ で平均ゼロ、$\|a\|_2\le |Q|^{-1/2}$ を満たす関数は $\|\mathcal Ma\|_1\le C_n$ である。従って有界原子にもこの評価が成り立つ。
:::

:::proof
$Q$ を固定倍に拡大した内部では $\mathcal M a\le CM a$、$M$ の $L^2$ 有界性と Cauchy–Schwarz により $\int \mathcal M a\le C|Q|^{1/2}\|a\|_2$ である。外部では $c_Q$ を中心とし、試験関数の台が $Q$ と交わるには $t\ge c|x-c_Q|$ が必要である。平均ゼロなので $\phi_t(x-c_Q)$ を引くと、作用の絶対値は $C\ell(Q)t^{-n-1}\|a\|_1$ 以下である。$\|a\|_1\le1$ であり、$|x-c_Q|^{-n-1}$ を外部で積分すれば $C/\ell(Q)$ となる。
:::

:::theorem[原子分解]{#thm-four-hardy-atoms}
$f\in H^1$ なら $f=\sum_j\lambda_j a_j$ が $L^1$ で成り立つような有界原子と $\sum_j|\lambda_j|\le C_n\|f\|_{H^1}$ が存在する。逆にこの形の $L^1$ 和は $H^1$ に属し、ノルムは $C_n\sum_j|\lambda_j|$ 以下である。
:::

:::proof
逆方向は :ref[lem-four-hardy-atom-bound] を使う。$L^1$ の和を固定した $\phi_t$ に作用させると一様に収束するので、$\mathcal M f\le\sum_j|\lambda_j|\mathcal M a_j$ である。積分して評価を得る。

順方向では $\Omega_k=\{\mathcal Mf>2^k\}$ を考える。これは有限測度の開集合であり、空の場合は以下の和を空和とする。:ref[thm-four-whitney] の分割を取り、辺長を $\ell_i^k$ とする。そこで構成した $\eta_i^k$ の正規化による分割を $\psi_i^k$ と書く。$\int\psi_i^k\asymp |Q_i^k|$ である。実際、分母は有界重なりで上から抑えられ、分子は $Q_i^k$ 上で $1$ である。

$$
c_i^k=\frac{\int f\psi_i^k}{\int\psi_i^k},\qquad
b_i^k=(f-c_i^k)\psi_i^k,\qquad
g^k=f-\sum_i b_i^k
$$

と置く。各 $Q_i^k$ の中心から $C\ell_i^k$ 以内に $\Omega_k$ の外の点 $z$ がある。$\psi_i^k$ はその点から半径 $C'\ell_i^k$ 以内に台を持ち、勾配が $C/\ell_i^k$ である。従って $t=C'\ell_i^k$ に拡大し向きを反転した試験関数として、$\psi_i^k(w)/(C''t^n)=\phi_t(z-w)$、$\phi\in\mathcal A$ と書ける。$\mathcal M f(z)\le2^k$ より $|c_i^k|\le C2^k$ である。$\int b_i^k=0$ であり、$g^k$ は $\Omega_k$ の外で $f$、内部で $\sum_i c_i^k\psi_i^k$ である。従って $\|g^k\|_\infty\le C2^k$、$\|g^k-f\|_1\le\int_{\Omega_k}|f|+C2^k|\Omega_k|$ である。$\mathcal M f$ が可積分なので後者は $k\to\infty$ でゼロとなる。$k\to-\infty$ では $g^k\to0$ が超関数として成り立つ。これで差 $g^{k+1}-g^k$ の和が $f$ になる。

差の各部品を積分ゼロにするため

$$
d_{ij}=\frac{\int\psi_i^k b_j^{k+1}}{\int\psi_j^{k+1}},\qquad
h_i^k=b_i^k-\sum_j(\psi_i^k b_j^{k+1}-d_{ij}\psi_j^{k+1})
$$

と置く。交わらない台の項はゼロである。交わる場合、$\Omega_{k+1}\subset\Omega_k$ と Whitney の距離比較より $\ell_j^{k+1}\le C\ell_i^k$ である。従って $\psi_i^k\psi_j^{k+1}$ は第 $j$ 立方体の尺度で勾配が $C/\ell_j^{k+1}$ 以下である。第 $j$ 立方体の近くの $\Omega_{k+1}$ の外の点を使い、$c_i^k$ の評価と同じ手順をこの積に適用すると $|\int f\psi_i^k\psi_j^{k+1}|\le C2^{k+1}|Q_j^{k+1}|$ である。$c_j^{k+1}$ の評価も組み合わせれば $|d_{ij}|\le C2^{k+1}$ を得る。

$h_i^k$ は $Q_i^k$ の固定倍の立方体に台を持つ。積分は各項の定義によりゼロである。また $\Omega_{k+1}$ 上で $\sum_j\psi_j^{k+1}=1$ なので、$b_i^k-\sum_j\psi_i^k b_j^{k+1}$ の $f$ の項はその外側だけに残り、そこでは $|f|\le C2^{k+1}$ である。定数の項と補正項も、$\sum_j\psi_j^{k+1}\le1$ と上の評価により $C2^{k+1}$ 以下である。よって $h_i^k$ は $C2^k|Q_i^k|$ 倍の原子である。

$\sum_i d_{ij}=0$ は $\sum_i\psi_i^k=1$ が第 $j$ 立方体の拡大上で成り立ち、$\int b_j^{k+1}=0$ であることによる。和の交換は $\sum_i|\int\psi_i^k b_j^{k+1}|\le\|b_j^{k+1}\|_1$ と $\sum_j\|b_j^{k+1}\|_1<\infty$ で正当化できる。従って $\sum_i h_i^k=g^{k+1}-g^k$ である。係数の絶対値和は $C\sum_k2^k|\Omega_k|$ 以下であり、各点において $\sum_{2^k<\mathcal M f(x)}2^k\le2\mathcal M f(x)$ だから収束する。原子和の $L^1$ 極限と望遠和の超関数極限が一致し、目的の表示になる。
:::

:::example[左右で打ち消す入力]
$a=1_{[0,1/2)}-1_{[1/2,1)}$ は有界原子であり $\|a\|_1=1$ である。$x>2$ で Hilbert 変換を積分すると

$$
R_1a(x)=\pi^{-1}\log\frac{x(x-1)}{(x-1/2)^2}
$$

となる。分子と分母の差は $-1/4$ で、対数は $O(x^{-2})$ である。半直線の例で現れた対数の増大と異なり、ここでは打ち消しで尾が可積分になる。
:::

各 $f\in H^1$ は $\int f=0$ を満たす。原子和が $L^1$ 収束し、各項の積分がゼロだからである。例えば $1_{[0,1]}$ は $L^1$ だが $H^1$ に属さない。また、:ref[thm-four-cz] の仮定下では $T$ は $H^1\to L^1$ に有界に延長できる。原子の拡大立方体の内側は $L^2$ 評価、外側は Hörmander 条件で $\int|Ta|\le C$ となる。有限原子和での評価を全体へ移す際は、任意の原子表示の部分和に弱型評価を使って $T$ の測度極限との一致を得る。従って表示によらず、原子係数の評価から結論が従う。

:::example[原子を二つ足すときの係数]
$a=1_{[0,1/2)}-1_{[1/2,1)}$、$b(x)=a(x-2)$ とする。$f=3a-2b$ は二つの原子の和であり、$\|f\|_1=5$ である。台が交わらないためである。原子分解の逆向きから $\|f\|_{H^1}\le5C_n$ となる。$\|f\|_1$ だけでは $H^1$ への所属を判定できないが、各成分の平均ゼロと台を示せばこの評価まで得られる。
:::

### 完備性、原子の二乗条件、滑らかな近似

原子分解で使った有界性は、各原子を作るときに便利であった。すでに立方体と打ち消しが分かっている入力では、点ごとの上限の代わりに二乗平均を使うと、Plancherel や Hilbert 空間の議論を適用できる。

:::proposition[$H^1$ の完備性と二乗原子]{#prop-four-hardy-l2-atoms}
$H^1$ は Banach 空間である。また有界原子の代わりに、立方体 $Q$ に台を持ち $\int a=0$、$\|a\|_2\le |Q|^{-1/2}$ を満たす二乗原子の絶対可和な和を使っても同じ空間を得る。両方の表示で係数の絶対値和の下限は $\|f\|_{H^1}$ と同値である。
:::

:::proof
$\mathcal M$ の劣加法性と斉次性で三角不等式を得る。$|f|\le C\mathcal M f$ a.e. により、ノルムがゼロなら $f=0$ である。$H^1$ の Cauchy 列 $f_k$ は同じ評価で $L^1$ の Cauchy 列となり、$L^1$ 極限 $f$ を持つ。固定した $\phi,t$ の畳み込みは $L^1$ から $L^\infty$ へ連続なので、各 $x$ で
$$
\mathcal M(f_k-f)(x)\le\liminf_{j\to\infty}\mathcal M(f_k-f_j)(x).
$$
これは各 $\phi,t$ について左の値以下であることを確かめた後で上限を取った結果である。Fatou により、$k$ が十分大きければ $\|\mathcal M(f_k-f)\|_1$ は任意に小さくできる。これで $f\in H^1$ と分かり、完備性が従う。

有界原子は Cauchy–Schwarz の尺度で二乗条件を満たすので、:ref[thm-four-hardy-atoms] は二乗原子による表示も与える。逆向きでは全ての二乗原子が $\|a\|_1\le1$ を満たすため、係数が絶対可和なら和は $L^1$ で収束する。:ref[lem-four-hardy-atom-bound] と固定した試験関数での極限を用いると $\|f\|_{H^1}\le C\sum|\lambda_j|$ である。両方の原子の下限についてこの評価と原子分解の上界を合わせる。
:::

:::corollary[滑らかなコンパクト台関数による近似]{#cor-four-hardy-smooth-density}
$C_c^\infty(\R^n)$ に属して積分が $0$ の関数は $H^1$ に稠密である。
:::

:::proof
有限原子和は :ref[thm-four-hardy-atoms] の係数の尾の評価により稠密である。有界原子 $a$ を取り、質量 $1$ の非負の滑らかなコンパクト台近似単位 $\rho_\eps$ と畳み込む。$a_\eps=\rho_\eps*a$ は滑らかで積分が $0$ であり、$\eps$ を十分小さくすると $a$ と $a_\eps$ の台は同じ固定拡大立方体 $R$ に含まれる。移動の $L^2$ 連続性で $\|a-a_\eps\|_2\to0$ である。$\int(a-a_\eps)=0$ なので :ref[lem-four-hardy-atom-bound] を正規化して使えば $\|a-a_\eps\|_{H^1}\le C|R|^{1/2}\|a-a_\eps\|_2\to0$ となる。有限個の原子を同様に近似し、それらの和を取れば結論が従う。
:::

:::example[非有界な二乗原子]
$Q=(0,1)$ で $v=x^{-1/4}-4/3$ とし、外ではゼロとする。$\int_0^1x^{-1/4}dx=4/3$ なので平均ゼロである。$\|v\|_2^2=2-(4/3)^2=2/9$ だから $a=3v/\sqrt2$ は二乗原子である。$a$ は原点の近くで非有界だが、上の近似ではその特異性を $L^2$ で平滑化し、平均ゼロを保つことで $H^1$ の誤差まで小さくできる。
:::

### 遠くへ減衰する部品を原子に分ける

特異積分に原子を入力すると、出力は一般に原子の立方体から外へ広がる。ただし殻ごとに減衰するなら、各殻の平均を補正して再び原子に分けられる。全体の平均ゼロは、殻で引いた定数をまとめる際に使う。

:::lemma[殻ごとに減衰する部品の分解]{#lem-four-hardy-molecule}
立方体 $Q$ を固定し、$Q_j=2^jQ$、$E_0=Q$、$E_j=Q_j\setminus Q_{j-1}$（$j\ge1$）とする。$u\in L^1$ が $\int u=0$ かつ
$$
\|u1_{E_j}\|_2\le A\,2^{-\delta j}|Q_j|^{-1/2}\qquad (j\ge0)
$$
を $\delta>0$ について満たすなら、$u\in H^1$ で $\|u\|_{H^1}\le C_{n,\delta}A$ である。
:::

:::proof
$u_j=u1_{E_j}$、$c_j=\int u_j$、$e_j=1_{Q_j}/|Q_j|$ とする。Cauchy–Schwarz により $|c_j|\le A2^{-\delta j}$ である。$u_j-c_je_j$ は $Q_j$ に台を持つ平均ゼロの関数で、二乗ノルムは $2A2^{-\delta j}|Q_j|^{-1/2}$ 以下である。従ってこれらの和は係数可和な二乗原子和となる。

定数部分について $s_j=\sum_{k=0}^j c_k=-\sum_{k>j}c_k$ と置く。全平均ゼロによって $|s_j|\le C_\delta A2^{-\delta j}$ である。有限和について
$$
\sum_{j=0}^N c_j e_j=\sum_{j=0}^{N-1}s_j(e_j-e_{j+1})+s_N e_N
$$
と望遠和で計算できる。$\|s_N e_N\|_1=|s_N|\to0$ である。$e_j-e_{j+1}$ は $Q_{j+1}$ を台とする平均ゼロの関数で、二乗ノルムは $(2^{n/2}+1)|Q_{j+1}|^{-1/2}$ 以下である。従って定数部分も係数可和な二乗原子和になる。:ref[prop-four-hardy-l2-atoms] と二つの係数和の上界から結論が従う。
:::

:::corollary[滑らかな畳み込み型特異積分の Hardy 空間上の有界性]{#cor-four-hardy-cz-stability}
$T$ が :ref[thm-four-cz] の仮定を満たす $L^2$ 上の Fourier 乗算子であり、核が $K(x-y)$ の形で $|\nabla K(z)|\le C|z|^{-n-1}$ を満たすとする。このとき $T:H^1\to H^1$ が有界である。
:::

:::proof
有界原子 $a$ の台を $Q$ とする。十分大きい固定倍 $R$ の $Q$ 上では $\|Ta\|_2\le\|T\|_{2\to2}|Q|^{-1/2}$ である。$R$ の外では $a$ の平均がゼロなので、$K(x-y)-K(x-c_Q)$ を用い
$$
|Ta(x)|\le C\ell(Q)|x-c_Q|^{-n-1}\|a\|_1
$$
を得る。$\|a\|_1\le1$ であり、この尾は可積分である。また殻 $E_j=2^jR\setminus2^{j-1}R$ 上の二乗ノルムは $C2^{-j}|2^jR|^{-1/2}$ 以下である。内部 $R$ の二乗評価を合わせれば :ref[lem-four-hardy-molecule] の条件が、平均ゼロを除いて成立する。

$Ta\in L^1\cap L^2$ なので $\widehat{Ta}$ は連続な代表元を持つ。$L^2$ の恒等式 $\widehat{Ta}=m\widehat a$ と乗数の本質上限から $|\widehat{Ta}(\xi)|\le\|m\|_\infty|\widehat a(\xi)|$ a.e. である。両辺の連続性により全点で成り立つ。$\widehat a(0)=0$ だから $\widehat{Ta}(0)=0$ であり、$Ta$ は平均ゼロである。従って $\|Ta\|_{H^1}\le C$ が全ての原子で一様に成り立つ。

任意の原子表示の部分和に $T$ を掛けると、この評価と完備性で $H^1$ に収束する。一方、$T$ の弱型 $(1,1)$ 延長によれば同じ部分和は $Tf$ に測度収束する。$H^1$ 収束は $L^1$ 収束を含むので、二つの極限は一致し表示には依存しない。係数和の下限を取れば作用素の評価が従う。
:::

:::example[一段ずつ移した箱の差]
$\R$ で $Q=[-1,1]$、$e_j=1_{2^jQ}/|2^jQ|$ とする。$u=\sum_{j\ge0}2^{-j}(e_j-e_{j+1})$ では各差は平均ゼロであり、これを $(\sqrt2+1)$ で割ると二乗原子になる。従って $H^1$ ノルムは $C(\sqrt2+1)\sum_{j\ge0}2^{-j}=2C(\sqrt2+1)$ 以下である。例えば $x\in Q$ では $e_j(x)-e_{j+1}(x)=2^{-j-2}$ であるから $u(x)=\sum_j2^{-2j-2}=1/3$ となる。全積分がゼロなのは、$Q$ 外の裾がこの正の部分を打ち消すためである。
:::

### Poisson 平均による境界での測定

大極大関数では複数の試験関数を同時に使って原子を作った。調和拡張を調べる場合には Poisson 核一つを使いたい。以下では、原子から非接極大関数への向きの評価を直接示す。

:::proposition[Poisson 非接極大関数の Hardy 評価]{#prop-four-hardy-poisson-maximal}
$a>0$ を固定し、$N_a f(x)=\sup_{t>0,\ |z-x|<at}|P_t*f(z)|$ とする。$f\in H^1$ なら $\|N_a f\|_1\le C_{n,a}\|f\|_{H^1}$ である。
:::

:::proof
台が $Q$ の原子 $v$ を考える。$N_a v\le C_{n,a}Mv$ は :ref[thm-four-fatou-poisson] の殻の評価であり、$M$ の $L^2$ 有界性を使えば固定拡大立方体での $L^1$ 積分は $C_{n,a}|Q|^{1/2}\|v\|_2$ 以下である。外側では $r=|x-c_Q|$ とする。$|z-x|<at$ のもとで、$t\ge r/(2(a+1))$ と $t<r/(2(a+1))$ に分ける。前者では $\|\nabla P_t\|_\infty\le Ct^{-n-1}\le C_{n,a}r^{-n-1}$ である。後者では $|z-c_Q|\ge r/2$ となり、$Q$ を十分離れた外部では $c_Q$ と $y\in Q$ を結ぶ線分上でも距離は $r/4$ 以上である。核を直接微分すれば、この線分上で $|\nabla P_t(z-\cdot)|\le C_{n,a}r^{-n-1}$ となる。どちらの場合も $v$ の平均ゼロで $P_t(z-c_Q)$ を引けば、$|P_t*v(z)|\le C_{n,a}\ell(Q) r^{-n-1}\|v\|_1$ である。従って外部でも積分は定数以下である。

原子の絶対可和和は $L^1$ に収束する。固定した $t>0,z$ では $P_t*f(z)$ はこの和を通して計算できるため、$N_a f\le\sum|\lambda_j|N_a v_j$ となる。積分し、係数和を評価すれば結論を得る。
:::

:::theorem[Poisson 非接極大関数による特徴付け]{#thm-four-hardy-poisson-characterization}
$f\in L^1(\R^n)$ とする。ある $a>0$ について $N_a f\in L^1$ なら $f\in H^1$ であり、$\|f\|_{H^1}\le C_{n,a}\|N_a f\|_1$ である。
:::

:::proof
$U(y,t)=P_t*f(y)$、$N=N_af$ とし、$n/(n+1)<r<1$ を固定する。まず半径 $c_at$ の空間立方体 $I$ を考える。$y,x\in I$ なら $|y-x|\le C_nc_at$ である。$c_a$ を十分小さく取ると、$(y,t)$ の周りの半径 $c_at$ の $(n+1)$ 次元の球内の全ての点 $(z,s)$ について $|z-x|<as$ が成り立つ。従ってその球内の $|U|$ の上限は $\inf_I N$ 以下である。調和関数の平均値性から、中心での勾配は $C(c_at)^{-1}$ と球上の絶対値の上限の積以下になる。実際、各導関数に平均値性を使い、球上の微分の積分を発散定理で球面積分に換えれば、球面積と球の体積の比が現れる。以上から

$$
\sup_{y\in I}\bigl(|U(y,t)|+t|\nabla_{y,t}U(y,t)|\bigr)
\le C_{n,a}\left(\frac1{|I|}\int_I N^r\right)^{1/r}.
$$

$0<t\le\eps$ に対し、$B(x,\eps)$ を辺長が $c_at$ の互いに交わらない格子立方体で覆う。必要な立方体は全て $B(x,C_a\eps)$ に含まれる。上の評価を各立方体で積分し、非負数 $A_I$ に対する $\sum_I A_I^{1/r}\le(\sum_I A_I)^{1/r}$ を使う。$|I|=(c_at)^n$ なので

$$
\int_{B(x,\eps)}\bigl(|U(y,t)|+t|\nabla_{y,t}U(y,t)|\bigr)\,dy
\le C t^{-d}\left(\int_{B(x,C_a\eps)}N^r\right)^{1/r},
\qquad d=n(1/r-1)<1 .
$$

次に :ref[def-four-real-hardy] の任意の $\phi$ と $\eps>0$ に対し $F(t)=\int\phi_\eps(x-y)U(y,t)dy$ と置く。$F(0)$ は $\phi_\eps*f(x)$ として定義でき、近似単位性から $F(t)\to F(0)$ である。$0<\delta<\eps$ について積分形の Taylor 公式を上端 $\eps$ から下端 $\delta$ に使うと

$$
F(\delta)=F(\eps)-(\eps-\delta)F'(\eps)
+\int_\delta^\eps(t-\delta)F''(t)dt.
$$

調和性 $\partial_t^2U=-\Delta_y U$ により $F''(t)$ は $\int\nabla_y[\phi_\eps(x-y)]\cdot\nabla_y U(y,t)\,dy$ に等しい。従って $t|F''(t)|\le C\eps^{-n-1}t^{-d}(\int_{B(x,C_a\eps)}N^r)^{1/r}$ であり、$d<1$ なので $\delta\downarrow0$ を通せる。端の二項は $t=\eps$ の評価と $\|\phi_\eps\|_\infty\le\eps^{-n}$ で抑えられる。積分項では $\int_0^\eps t^{-d}dt=\eps^{1-d}/(1-d)$ を使う。三項を合わせると

$$
|\phi_\eps*f(x)|
\le C\eps^{-n-d}\left(\int_{B(x,C_a\eps)}N^r\right)^{1/r}
\le C[M(N^r)(x)]^{1/r},
$$

ここで $n+d=n/r$ を用いた。定数は試験関数と $\eps$ に依存しないから、同じ評価が $\mathcal M f$ に成り立つ。$1/r>1$ における $M$ の有界性を $N^r$ に使うと
$\|\mathcal M f\|_1\le C\|M(N^r)\|_{1/r}^{1/r}\le C\|N\|_1$ となる。
:::

:::corollary[鉛直極大関数による特徴付け]{#cor-four-hardy-radial-maximal}
$f\in L^1(\R^n)$ に対し $Vf(x)=\sup_{t>0}|P_t*f(x)|$ と置く。$Vf\in L^1$ と $f\in H^1$ は同値であり、$\|Vf\|_1\asymp\|f\|_{H^1}$ である。
:::

:::proof
$f\in H^1$ なら $Vf\le N_a f$ と :ref[prop-four-hardy-poisson-maximal] により評価できる。逆向きには、調和関数 $u$ の内部平均値評価を小さい冪まで拡張する。$0<q<1$、$u$ が閉球 $B(X,R)$ の近傍で調和的とする。半径 $R_j=R(1-2^{-j-1})$ の球上の上限を $S_j$、$B(X,R)$ 上の $|u|^q$ の平均を $L$、空間次元を $d$ とする。$B(X,R_j)$ 内の各点で半径 $R_{j+1}-R_j$ の平均値公式を使うと、
$$
S_j\le C_d\,2^{dj} S_{j+1}^{1-q}L
$$
となる。実際、平均に現れる $|u|$ を $|u|^q S_{j+1}^{1-q}$ で抑え、体積比を掛ければよい。任意の固定した $\eta>0$ について右辺は、積に対する Young の不等式で $\eta S_{j+1}+C_{d,q,\eta}2^{dj/q}L^{1/q}$ 以下となる。$\eta<2^{-d/q}$ を選んで繰り返せば、残る $\eta^jS_j$ はゼロに収束し、係数の和は有限である。従って $|u(X)|\le S_0\le C_{d,q}L^{1/q}$ である。

この評価を $u(z,t)=P_t*f(z)$ に対し、$(z,t)$ を中心とする半径 $t/2$ の $(n+1)$ 次元の球で使う。$|z-x|<at$ なら、球の水平射影は $B(x,(a+1)t)$ に含まれ、各水平点で $|u|$ は $Vf$ 以下である。高さ方向の長さが高々 $t$ なので、球平均は $C_{a,n}M((Vf)^q)(x)$ 以下になる。よって $N_a f\le C[M((Vf)^q)]^{1/q}$ である。$M$ の $L^{1/q}$ 有界性により $\|N_a f\|_1\le C\|Vf\|_1$ となり、:ref[thm-four-hardy-poisson-characterization] を使えばよい。
:::

### 共役調和系による Riesz 変換の判定

一つの調和関数の絶対値は劣調和だが、絶対値の $1$ 未満の冪に同じ性質が成り立つとは限らない。Poisson 拡張した $f$ と全ての Riesz 変換を成分として並べると、導関数の行列が対称でトレースがゼロになる。この条件によって、長さの $1$ 未満の冪にも劣調和性が残る。

:::lemma[共役調和系の冪の劣調和性]{#lem-four-hardy-riesz-subharmonic}
$d\ge2$ とし、$\R^d$ の開集合上で $U$ は $\R^d$ 値の調和関数とする。すなわち各成分が調和であるとする。Jacobian $A=DU$ が対称で $\operatorname{tr}A=0$ なら、$(d-2)/(d-1)<q<1$ に対し $|U|^q$ は劣調和である。
:::

:::proof
実対称行列 $A$ の固有値を $\lambda_1,\ldots,\lambda_d$ とする。$\sum_j\lambda_j=0$ なので、各 $j$ について
$$
|\lambda_j|^2=\left|\sum_{k\ne j}\lambda_k\right|^2
\le(d-1)\sum_{k\ne j}|\lambda_k|^2 .
$$
よって最大固有値の絶対値の二乗は $(d-1)\sum_k|\lambda_k|^2/d$ 以下である。$\|A\|_{\mathrm{HS}}^2=\sum_{j,k}|A_{jk}|^2=\sum_k|\lambda_k|^2$ と書くと、任意の $v$ について $|Av|^2\le(d-1)\|A\|_{\mathrm{HS}}^2|v|^2/d$ となる。

$W=|U|^2+\eps$、$\eps>0$ とする。調和性から $\Delta W=2\|A\|_{\mathrm{HS}}^2$、対称性から $\nabla W=2AU$ であり、連鎖律は
$$
\Delta(W^{q/2})
=qW^{q/2-1}\left(\|A\|_{\mathrm{HS}}^2+(q-2)\frac{|AU|^2}{W}\right)
$$
を与える。括弧は $(1-(2-q)(d-1)/d)\|A\|_{\mathrm{HS}}^2$ 以上であり、これは $q>(d-2)/(d-1)$ のため非負である。従って $W^{q/2}$ は劣調和である。閉球上の平均値不等式を $\eps\downarrow0$ の一様収束で移せば、$|U|^q$ も劣調和になる。
:::

:::theorem[Riesz 変換による実 Hardy 空間の特徴付け]{#thm-four-hardy-riesz-characterization}
$f\in L^1(\R^n)$ に対し、$\widehat{R_j f}(\xi)=-i\xi_j|\xi|^{-1}\widehat f(\xi)$ を有界関数と見て $R_j f$ を緩増加超関数として定める。$\xi=0$ での値は任意とする。全ての $R_j f$ が $L^1$ 関数で表されることと $f\in H^1$ は同値であり、
$$
\|f\|_{H^1}\asymp\|f\|_1+\sum_{j=1}^n\|R_j f\|_1
$$
である。
:::

:::proof
$f\in H^1$ のときは :ref[cor-four-hardy-cz-stability] から $L^1$ 評価が従う。有限原子和を $f$ に $H^1$ 近似すると、入力と出力の両方が $L^1$ で収束する。従って変換は一様に収束し、$L^2$ 上の乗数公式が a.e. の等式として極限へ移る。これで超関数として定めた $R_j f$ とこの延長が一致する。

逆向きではまず実数値の $f$ を取る。Riesz 乗数は $m(-\xi)=\overline{m(\xi)}$ を満たしており、各 $g_j=R_j f\in L^1$ は実数値である。$g_0=f$ とし $U_j(x,t)=P_t*g_j(x)$ と置く。座標の順を $(t,x_1,\ldots,x_n)$ と取ると、変換の式から $\partial_t U_j=\partial_{x_j}U_0$、$\partial_{x_k}U_j=\partial_{x_j}U_k$、$\partial_t U_0+\sum_j\partial_{x_j}U_j=0$ が従う。実際、$t$ 微分の乗数は $-2\pi|\xi|$、$x_j$ 微分の乗数は $2\pi i\xi_j$ であり、Riesz 乗数を掛け合わせれば各式が一致する。導関数の逆積分は $t>0$ で指数減衰により絶対収束する。各成分は調和であり、:ref[lem-four-hardy-riesz-subharmonic] により $n/(n+1)<q<1$ なら $V=|U|^q$ が劣調和である。この $q$ の下限は $(n-1)/n$ より大きい。

$\delta>0$ とし $h_\delta(x)=|U(x,\delta)|^q$ と置く。$V(x,t+\delta)\le P_t*h_\delta(x)$ を示す。各 $g_j\in L^1$ なので $U(x,s)\to0$ は $s\to\infty$ で $x$ に一様である。実際 $|P_s*g_j|\le C s^{-n}\|g_j\|_1$ である。$\delta\le s\le S$ を固定すれば、$g_j$ を大きな球の内外に分け、内側で核の一様な空間減衰、外側で $C\delta^{-n}$ と $L^1$ の尾を使うと、$|x|\to\infty$ でも一様にゼロになる。従って $V(x,t+\delta)$ は閉上半空間で無限遠に向かうとゼロになる。$h_\delta$ は連続で無限遠でゼロへ向かい、その Poisson 拡張は $t=0$ まで連続で同じ境界値を持つ。二つの差は劣調和で、底面ではゼロ、上半球の球面では非負な $P_t*h_\delta$ を引いているので上限は $\sup V(x,t+\delta)$ 以下である。半径を無限大へ移す最大値原理により、差は全点で非正である。

$P_\delta*g_j\to g_j$ は $L^1$ で成り立つ。ベクトル $v,w$ に対する $\big||v|^q-|w|^q\big|\le|v-w|^q$ と $0<q<1$ により、$h_\delta$ は $h=(\sum_{j=0}^n|g_j|^2)^{q/2}$ に $L^{1/q}$ で収束する。従って固定した $(x,t)$ で Hölder により $P_t*h_\delta(x)\to P_t*h(x)$ となり、$|U_0(x,t)|^q\le P_t*h(x)$ を得る。:ref[thm-four-fatou-poisson] の非接評価により $N_a f\le C_{n,a}[Mh]^{1/q}$ である。$1/q>1$ での $M$ の有界性から
$$
\|N_a f\|_1\le C_{n,a}\|h\|_{1/q}^{1/q}
\le C_{n,a}\sum_{j=0}^n\|g_j\|_1 .
$$
:ref[thm-four-hardy-poisson-characterization] により求める $H^1$ 評価を得る。複素値の $f$ は実部と虚部に分ける。Riesz 変換は実関数を実超関数へ移すので、仮定と評価は両部に適用でき、三角不等式で結論が従う。
:::

:::example[平均ゼロを欠く箱の Poisson 平均]
$\R$ で $f=1_{[0,1]}$ とする。$x\ge2$、$t=x$ なら $0\le x-y\le x$ が $y\in[0,1]$ で成り立つので
$$
(P_x*f)(x)=\frac1\pi\int_0^1\frac{x}{x^2+(x-y)^2}\,dy
\ge\frac1{2\pi x}.
$$
従って $N_a f$ は可積分でない。高さを固定すると尾は $O(x^{-2})$ だが、極大関数では $t=x$ も選べるため $x^{-1}$ の大きさが残る。原子で境界の平均ゼロが必要になる理由を、この一例から直接確認できる。
:::

### 周波数の原点における打ち消し

原子は平均ゼロなので Fourier 変換の原点に零点を持つ。その零点の強さは台の直径で決まる。高周波側は二乗ノルムを使い、両側を合わせると Hardy 空間で有効な積分不等式を得る。

:::proposition[Fourier 変換に対する Hardy の不等式]{#prop-four-hardy-fourier-inequality}
$f\in H^1(\R^n)$ なら
$$
\int_{\R^n}\frac{|\widehat f(\xi)|}{|\xi|^n}\,d\xi
\le C_n\|f\|_{H^1}.
$$
:::

:::proof
台が $Q$、辺長が $\ell$ の原子 $a$ を取る。$c_Q$ を中心とし、$\int a=0$ であることから
$$
|\widehat a(\xi)|
\le\int_Q |e^{-2\pi i x\cdot\xi}-e^{-2\pi i c_Q\cdot\xi}||a(x)|\,dx
\le C_n\ell|\xi|\|a\|_1.
$$
従って $|\xi|\le1/\ell$ では重み付き積分が $C_n\ell\int_0^{1/\ell}dr=C_n$ 以下になる。高周波側では Cauchy–Schwarz と Plancherel から
$$
\int_{|\xi|>1/\ell}|\widehat a(\xi)||\xi|^{-n}\,d\xi
\le\|a\|_2\left(|\mathbb S^{n-1}|\int_{1/\ell}^\infty r^{-n-1}dr\right)^{1/2}
\le C_n\ell^{-n/2}\ell^{n/2}.
$$
各原子に一様な評価を得た。原子和は $L^1$ に収束するので変換は一様収束し、三角不等式と Tonelli で原子ごとの上界を足せる。係数和の上界から結論が従う。
:::

:::example[左右を引いた箱の低周波]
$a=1_{[0,1/2)}-1_{[1/2,1)}$ とすると
$$
\widehat a(\xi)=\frac{(1-e^{-\pi i\xi})^2}{2\pi i\xi}.
$$
従って $\xi\to0$ で $\widehat a(\xi)=i\pi\xi/2+O(\xi^2)$ であり、$|\widehat a(\xi)|/|\xi|$ は原点近くで有界である。高周波では分子の絶対値は $4$ 以下なので重みを掛けた関数は $2/(\pi|\xi|^2)$ 以下となる。平均を引かなかった箱では $\widehat f(0)=1$ なので、同じ重み付き積分が原点で対数発散する。
:::

## H¹ と BMO の双対性

原子は平均がゼロなので、原子と $b$ の積分では $b$ の局所平均を引ける。逆に $H^1$ 上の線形汎関数は、固定された立方体の中で平均ゼロの二乗可積分関数に制限すると局所的な関数で表される。

:::theorem[$H^1$–BMO 双対性]{#thm-four-hardy-bmo-duality}
$H^1(\R^n)$ の連続な複素線形双対は $\mathrm{BMO}(\R^n)$ と同型である。対応は原子上で $\Lambda_b(a)=\int a b$ であり、$\|\Lambda_b\|\asymp\|b\|_{\mathrm{BMO}}$ となる。定数を加えても汎関数は同じである。
:::

:::proof
実数値 $b$ に対し $b_N=\max(-N,\min(N,b))$ と置く。切断写像が $1$-Lipschitz なので、中心を $b_Q$ の切断値に取って :ref[prop-four-bmo-centers] を使うと $\|b_N\|_{\mathrm{BMO}}\le2\|b\|_{\mathrm{BMO}}$ である。複素値には実部と虚部に同じ操作を行う。任意の $f\in H^1$ を原子展開し、有界関数 $b_N$ との積分を和に入れる。各原子で $|\int a b_N|\le C\|b\|_{\mathrm{BMO}}$ であり、局所優収束により $\int a b_N\to\int a b$ である。係数が絶対可和なので $\int f b_N$ の極限が存在し、その絶対値は $C\|b\|_{\mathrm{BMO}}\|f\|_{H^1}$ 以下である。この定義は表示によらず、有限和上で目的の組み合わせと一致する。

逆に有界汎関数 $\Lambda$ を取る。立方体 $Q$ に台を持つ平均ゼロの $L^2$ 関数 $v$ について :ref[lem-four-hardy-atom-bound] から $\|v\|_{H^1}\le C|Q|^{1/2}\|v\|_2$ である。この閉部分空間上で Riesz の表現定理を用い、共役を代表元側に吸収すれば平均ゼロの $b^{(Q)}\in L^2(Q)$ で $\Lambda(v)=\int_Q v b^{(Q)}$、$\|b^{(Q)}\|_2\le C|Q|^{1/2}\|\Lambda\|$ を得る。包含 $Q\subset R$ では $b^{(R)}|_Q-b^{(Q)}$ が平均ゼロの全関数に直交するので定数である。増加する立方体列で代表元の定数をそろえ、$b\in L^2_{\mathrm{loc}}$ を作る。任意の $Q$ 上でも同じ議論で $b-b_Q=b^{(Q)}$ となる。Cauchy–Schwarz より $\|b\|_{\mathrm{BMO}}\le C\|\Lambda\|$ である。

両汎関数は全原子で一致する。:ref[thm-four-hardy-atoms] とその逆方向のノルム評価で有限原子和は $H^1$ に稠密なので全体で一致する。同じ汎関数を表す二関数の差は、どの立方体でも平均ゼロの $L^2$ 関数の積分に寄与しない。:ref[thm-four-john-nirenberg] で差は局所 $L^2$ であり、各立方体で定数、従って全空間でも定数である。
:::

### 双対性を使う際の積分の範囲

:::example[原子の位置が変わる場合]
$n=1$ で $a_R=R^{-1}(1_{[0,R/2)}-1_{[R/2,R)})$ とする。$b(x)=\log|x|$ との積分は $x=Rt$ の置換で
$$
\int a_R b=\int_0^{1/2}\log t\,dt-\int_{1/2}^{1}\log t\,dt=-\log2
$$
である。$\log R$ は原子の積分がゼロなので消える。固定の BMO 関数と組み合わせた値が、原子の尺度によらない例である。
:::

:::example[絶対積分と連続延長]
互いに離れた区間 $Q_j=[e^{j^2},e^{j^2}+1]$（$j\ge1$）を取り、各区間の前半と後半でそれぞれ $1,-1$ の原子 $a_j$ を取る。$f=\sum_{j\ge1}j^{-2}a_j$ は :ref[thm-four-hardy-atoms] により $H^1$ に属する。$b=\log|x|$ とすると $Q_j$ 上で $b\ge j^2$、$\int|a_j|=1$ なので $\int|fb|\ge\sum_{j\ge1}1=\infty$ である。一方、各原子との組み合わせは正負の部分を合わせた後の値であり、半ノルムで一様に抑えられる。従って双対汎関数の値 $\sum_jj^{-2}\int a_j b$ は絶対収束する。一般の $H^1$ 関数に対し、双対作用を $\int fb$ の絶対積分で置き換えられないことが分かる。
:::

### Riesz 変換による BMO の表示

原子を使う双対性では BMO の要素を線形汎関数として扱った。Riesz 変換で Hardy ノルムを測ると、その汎関数を有限個の $L^1$ 関数の組の上へ延長できる。この延長から、有界関数とその Riesz 変換を足して BMO を表す方法が得られる。

:::corollary[有界関数と Riesz 変換による表示]{#cor-four-bmo-riesz-decomposition}
任意の $b\in\mathrm{BMO}(\R^n)$ は定数を除いて
$$
b=g_0+\sum_{j=1}^n R_jg_j,\qquad
\max_{0\le j\le n}\|g_j\|_\infty\le C_n\|b\|_{\mathrm{BMO}}
$$
と表せる。ここで $R_jg_j$ は :ref[thm-four-cz-bmo] の延長である。
:::

:::proof
$X=(L^1)^{n+1}$ に成分ノルムの和を入れ、$Jf=(f,R_1f,\ldots,R_nf)$ とする。:ref[thm-four-hardy-riesz-characterization] により $J$ は $H^1$ のノルムと同値なノルムを像に与える。従って $\Lambda_b\circ J^{-1}$ は像上の連続線形汎関数である。Hahn–Banach で全体へ同じノルムで延長し、各成分を $L^1$ の双対で表すと
$$
\Lambda_b(f)=\int f g_0-\sum_{j=1}^n\int (R_j f)g_j
$$
となる。符号は代表元の選択に吸収でき、各 $\|g_j\|_\infty$ は延長の作用素ノルム以下である。

各有界原子 $a$ について、$g_j$ を $B(0,R)$ に切断した $g_{j,R}$ を取る。これは $L^2$ に属する。Riesz 乗数は純虚数で奇関数なので、$L^2$ の随伴は $-R_j$ である。核も実数値のため、複素共役なしの双線形積分で $\int a R_j g_{j,R}=-\int(R_j a)g_{j,R}$ が成り立つ。台 $Q$ 上で、左辺は $R_j g_j$ の局所表示の積分へ収束する。:ref[thm-four-cz-bmo] の証明どおり切断出力から定数を引いて $L^1(Q)$ 収束させられ、$a$ の平均ゼロにより定数は寄与しない。右辺も $R_j a\in L^1$ なので優収束を使える。従って表示の右辺と $b$ は原子上で同じ汎関数を定める。:ref[thm-four-hardy-bmo-duality] の一意性から定数を除いて一致する。
:::

:::example[対数は有界入力の変換だけで表せる]
$\R$ で $g_0=0$、$g_1=\pi1_{(0,\infty)}$ とすれば、前節の切断計算により $R_1g_1=\log|x|$ が BMO の元として成り立つ。例えば $a=1_{[0,1/2)}-1_{[1/2,1)}$ に対して $\Lambda_{\log|\cdot|}(a)=-\log2$ である。随伴の表示を使うと同じ値が $-\pi\int_0^\infty R_1a(x)\,dx$ に等しい。$R_1 a$ は可積分なので、この積分は通常の絶対収束積分である。
:::

### 双対ノルムを原子で測る

:::proposition[原子での試験と BMO 半ノルム]{#prop-four-bmo-atom-testing}
$B=\|b\|_{\mathrm{BMO}}$、$A=\sup_{a}|\int ab|$ とし、上限は有界原子全体で取る。このとき $B/2\le A\le B$ である。同じ上限を取る際に、コンパクト台の滑らかな原子に限定しても値は変わらない。
:::

:::proof
台が $Q$ の原子に対し $\int ab=\int_Qa(b-b_Q)$ なので $|\int ab|\le B$ である。逆向きには、$b-b_Q\ne0$ の点で $v=\overline{b-b_Q}/|b-b_Q|$、他で $v=0$ とする。
$$
a=(v-v_Q)1_Q/(2|Q|)
$$
は平均ゼロで $\|a\|_\infty\le|Q|^{-1}$ である。従って $\int ab=(2|Q|)^{-1}\int_Q|b-b_Q|$ であり、$Q$ の上限で $B/2\le A$ が従う。

有界原子 $a$ と質量 $1$ の正の滑らかな近似単位 $\rho_\eps$ を取る。$a*\rho_\eps$ の台を $Q$ の各面を $\eps$ だけ広げた立方体 $Q_\eps$ に入れられるよう近似単位の台を選ぶ。$(|Q|/|Q_\eps|)a*\rho_\eps$ は滑らかな原子であり、$a$ に $L^2$ 収束する。固定拡大立方体 $R$ に全ての台が含まれ、$b\in L^2(R)$ は John–Nirenberg による。従ってこの原子との積分は $\int ab$ に近づく。任意の原子についてできるので上限の値が一致する。
:::

:::example[跳躍に対して上限が達成される場合]
$b=1_{(0,\infty)}$、$Q=(-1,1)$ とする。$a=(1_{(0,1)}-1_{(-1,0)})/2$ は原子で、$\int ab=1/2$ である。$B=1/2$ なので $A=B$ である。上の滑らかな近似は、同じ値を滑らかな原子で任意の精度まで近似できることを示す。
:::

:::example[有界関数による BMO 近似に対する下界]
$b(x)=\log|x|$ とし、$g\in L^\infty(\R)$ とする。$c_R=1_{(-R,R)}/(2R)$ と置くと、$c_R-c_{2R}$ は $(-2R,2R)$ を台とする有界原子である。実際、$|c_R-c_{2R}|=1/(4R)$ がこの区間内のほとんど全ての点で成り立つ。一方、$Q_R=(-R,R)$ における $b-g$ の平均は $\log R-1-g_{Q_R}$ である。$m$ 個の尺度で平均の差を足すと
$$
m\log2-|g_{Q_{2^m R}}-g_{Q_R}|
\le m\sup_{r>0}\left|\int (c_r-c_{2r})(b-g)\right|
\le m\|b-g\|_{\mathrm{BMO}}.
$$
左辺の第 2 項は $2\|g\|_\infty$ 以下である。$m\to\infty$ とすれば $\|b-g\|_{\mathrm{BMO}}\ge\log2$ を得る。BMO が $L^\infty$ より広いという事実よりさらに強く、この対数は BMO ノルムで有界関数から正の距離にある。
:::

### Poisson 拡張から BMO を判定する

箱の二乗評価と原子との組み合わせを結び付ければ、Poisson 拡張だけを使って境界関数の平均振動も評価できる。

:::proposition[箱評価からの平均振動の復元]{#prop-four-bmo-poisson-carleson-converse}
$b\in L^1_{\rm loc}(\R^n)$ は $\int |b(x)|(1+|x|)^{-n-1}dx<\infty$ を満たすとする。$u=P_t*b$ とし、
$$
A^2=\sup_Q |Q|^{-1}\int_0^{\ell(Q)}\int_Q t|\partial_t u(x,t)|^2\,dx\,dt
$$
が有限なら $\|b\|_{\mathrm{BMO}}\le C_n A$ である。
:::

:::proof
滑らかな有界原子 $a$ を取り、台を含む立方体を $Q$ とする。$v=P_t*a$ と置く。まず
$$
\int b(x)a(x)dx=4\int_0^\infty\int t\,\partial_t u(x,t)\partial_t v(x,t)\,dx\,dt
$$
を示す。$t$ を $\eps\le t\le R$ に制限すると、Fubini と $P_t*P_s=P_{t+s}$ により右辺は $\int b\,[(J_{2\eps}-J_{2R})*a]$ に等しい。ここで $J_s=P_s-s\partial_sP_s$ である。実際、二つの高さ微分の積は $\partial_s^2P_s|_{s=2t}$ の畳み込みに対応し、$\partial_s J_s=-s\partial_s^2P_s$ を積分すればよい。固定した高さ区間では核と導関数の絶対値は $C_{\eps,R}(1+|x|)^{-n-1}$ で抑えられる。この重みの自己畳み込みも同じ重みの定数倍以下である。積分を $|y|\ge|x|/2$ と $|x-y|\ge|x|/2$ に分ければ分かる。従って $b$ の仮定により Fubini が使える。

$J_s(x)=(n+1)c_n s^3(s^2+|x|^2)^{-(n+3)/2}$ は非負で質量 $1$ の尺度核である。$s\downarrow0$ で $J_s*a\to a$ 一様であり、$s\le1$ では台の十分遠方で $C_a(1+|x|)^{-n-1}$ 以下である。$s\ge1$ では $\int a=0$ を使い、核から $Q$ の中心での値を引く。微分を評価すると $|J_s*a(x)|\le C_a (s+|x|)^{-n-1}$ となる。ここで $x$ が台から離れたときは線分上の距離を使い、有界な範囲では $s\ge1$ と微分の上限 $Cs^{-n-1}$ を使う。従って $R\to\infty$ と $\eps\downarrow0$ の両方で $b$ との積分に優収束が使える。

右辺の絶対収束と評価には、$Q$ を十分大きい固定倍の立方体 $Q_0$ に拡大し、箱 $2^j Q_0\times(0,\ell(2^jQ_0))$ を使う。最初の箱では $u$ 側の二乗積分は $A^2|Q_0|$ 以下、$v$ 側は Plancherel で全上半空間の積分が $\|a\|_2^2/4\le1/(4|Q|)$ である。$j\ge1$ の箱と一つ前の箱の差では、$r_j=2^j\ell(Q_0)$ と置くと $t+|x-c_Q|\ge c r_j$ である。従って $\int a=0$ と混合微分の核評価から $|\partial_t v|\le C\ell(Q)r_j^{-n-2}$ となる。領域で $t\,dx\,dt$ の積分は $Cr_j^{n+2}$ 以下なので、$v$ 側の二乗積分は $C\ell(Q)^2r_j^{-n-2}$ 以下である。Cauchy–Schwarz で $u$ 側の箱評価と合わせると、第 $j$ 項は $CA\ell(Q)/r_j$ 以下となる。足し合わせて $|\int ba|\le C A$ を得る。

任意の有界原子は非負の近似単位で平滑化し、台の立方体の拡大に応じて係数を調整すれば、滑らかな有界原子で a.e. に近似できる。台は共通の有界集合に入り、近似列も一様有界なので $b$ との積分は優収束で移る。最後に各立方体 $Q$ で $w=\overline{b-b_Q}/|b-b_Q|$ を取り、分母がゼロの点では $w=0$ と置く。$(w-w_Q)1_Q/(2|Q|)$ は有界原子であり、$b$ との積分は $\operatorname{avg}_Q|b-b_Q|/2$ である。従って平均振動は一様に $C A$ 以下になる。
:::

## 演習問題

::::exercise{level=basic}
$n=1$ で $b=2\,1_{(0,\infty)}$ の BMO 半ノルムを求めよ。
:::solution
区間の正側比率を $\theta$ とすると平均は $2\theta$、平均振動は $4\theta(1-\theta)$ である。これは $\theta=1/2$ で最大 $1$ となる。
:::
::::

::::exercise{level=basic}
$Q=[0,2]$ 上の $a=(1_{[0,1)}-1_{[1,2)})/2$ が原子であることを確かめよ。
:::solution
積分は $1/2-1/2=0$、上限は $1/2=|Q|^{-1}$ である。台は $Q$ に含まれるので条件を満たす。
:::
::::

::::exercise{level=basic}
直前の原子と $b(x)=x$ の $Q$ 上の積分を求めよ。
:::solution
$\int_Q ab=\frac12(\int_0^1 x\,dx-\int_1^2x\,dx)=\frac12(1/2-3/2)=-1/2$ である。$b$ に定数を加えても $\int a=0$ なので同じ値となる。
:::
::::

::::exercise{level=basic}
$Q=(0,1)$ で $b=-\log x$ の平均振動を求めよ。
:::solution
$b_Q=1$ である。$t=-\log x$ と置くと積分は $\int_0^\infty|t-1|e^{-t}dt$ である。$\int_0^1(1-t)e^{-t}dt=1/e$、$\int_1^\infty(t-1)e^{-t}dt=1/e$ なので $2/e$ である。
:::
::::

::::exercise{level=advanced}
$b\in\mathrm{BMO}$ に対して $|b_{2^kQ}-b_Q|\le 2^n k\|b\|_{\mathrm{BMO}}$ を示せ。
:::solution
$|b_Q-b_{2Q}|\le |Q|^{-1}\int_Q|b-b_{2Q}|\le2^n\|b\|_{\mathrm{BMO}}$ である。$Q$ を $2^jQ$ に替え、連続する平均の差を $j=0,\ldots,k-1$ で足す。
:::
::::

::::exercise{level=advanced}
$H^1$ の任意の非負関数はゼロであることを示せ。
:::solution
原子分解は $L^1$ で収束し、原子の積分はすべてゼロである。従って $f$ の積分はゼロである。非負性により $f=0$ a.e. である。
:::
::::

::::exercise{level=advanced}
$u$ が $Q$ に台を持つ平均ゼロの $L^2$ 関数なら、$\left|\int b u\right|\le C_n\|b\|_{\mathrm{BMO}}|Q|^{1/2}\|u\|_2$ を示せ。
:::solution
平均ゼロにより $\int b u=\int_Q (b-b_Q)u$ である。Cauchy–Schwarz を使い、:ref[thm-four-john-nirenberg] の $q=2$ の評価を代入する。
:::
::::

::::exercise{level=advanced}
$\|b\|_{\mathrm{BMO}}$ を有界原子との積分から下から評価せよ。
:::solution
$Q$ を固定し、$g=\overline{b-b_Q}/|b-b_Q|$ を差がゼロの場所ではゼロと置く。$a=(g-g_Q)1_Q/(2|Q|)$ は原子である。$\int_Q(b-b_Q)=0$ より $\int a b=(2|Q|)^{-1}\int_Q|b-b_Q|$ である。$Q$ の上限を取ると $\sup_a|\int ab|\ge\|b\|_{\mathrm{BMO}}/2$ となる。
:::
::::

::::exercise{level=basic}
$b(x)=\log|x|$、$Q=(-4,4)$ とする。$b_Q-b_{(-1,1)}$ を求めよ。
:::solution
$\int_0^r\log x\,dx=r\log r-r$ である。対称性で $b_{(-r,r)}=\log r-1$ となる。従って差は $\log4=2\log2$ である。
:::
::::

::::exercise{level=basic}
$w(x)=x^{1/2}$（$x>0$）について、$Q=(0,2)$ で $p=2$ の $A_p$ の積を求めよ。
:::solution
$\operatorname{avg}_Q w=(2/3)2^{1/2}$、$\operatorname{avg}_Q w^{-1}=2\cdot2^{-1/2}$ である。積は $4/3$ になる。この計算は $Q$ 一つだけの評価であり、全区間での評価にはべき重みの例の比較が必要である。
:::
::::

::::exercise{level=advanced}
$n=1$ で区間 $Q$ に対し $u=1_Q/|Q|-1_{4Q}/|4Q|$ とする。$u$ を二つの有界原子の一次結合で表し、係数の絶対値和を $2$ 以下とせよ。
:::solution
$e_j=1_{2^jQ}/|2^jQ|$ と置くと $u=(e_0-e_1)+(e_1-e_2)$ である。各差は平均ゼロである。内側 $2^jQ$ で $e_j-e_{j+1}=(2^n-1)/|2^{j+1}Q|$、残りで $-1/|2^{j+1}Q|$ となるので、$a_j=(e_j-e_{j+1})/(2^n-1)$ は $2^{j+1}Q$ を台とする有界原子である。$n=1$ なので $2^n-1=1$ となり、係数はどちらも $1$ である。絶対値和は $2$ になる。
:::
::::

::::exercise{level=advanced}
$a=1_{[0,1/2)}-1_{[1/2,1)}$ に対し、$\int xa(x)\,dx$ と $\widehat a'(0)$ を求めよ。
:::solution
$\int xa(x)\,dx=\int_0^{1/2}x\,dx-\int_{1/2}^1x\,dx=1/8-3/8=-1/4$ である。コンパクト台なので変換を積分内で微分でき、$\widehat a'(0)=-2\pi i\int xa(x)\,dx=i\pi/2$ となる。これが低周波の一次の係数である。
:::
::::

## 確認テスト

::::quiz{answer="0" tol="1e-6"}
任意の実 Hardy 空間の原子の積分を求めよ。
:::explanation
定義で各原子は平均ゼロの条件を満たす。従って積分は $0$ であり、この性質によって積分核の定数項が打ち消される。
:::
::::

::::quiz
John–Nirenberg の不等式が与えるのはどの評価か。

- [x] 平均からの差が大きい集合の測度の指数減衰
- [ ] 全ての BMO 関数の本質上限の有限性
- [ ] 各立方体の平均値がゼロになる性質

:::explanation
平均を超える場所を子立方体に分け、各段でその体積を半分以下にする。閾値は段数に比例するため、その集合の測度は指数的に減る。$\log|x|$ は BMO だが非有界であり、平均は尺度に応じて変わる。
:::
::::

::::quiz{answer="0.5" tol="1e-6"}
$b=1_{(0,\infty)}$ の $Q=[-1,1]$ における平均振動を求めよ。
:::explanation
左右の長さが同じなので $b_Q=1/2$ である。両側で $|b-b_Q|=1/2$ になり、その平均は $1/2$ である。
:::
::::

::::quiz
$H^1$ と BMO の組み合わせに対してどれが成り立つか。

- [x] $b$ に定数を足しても汎関数は変わらない。
- [ ] $b$ が一定なら汎関数は常に非ゼロである。
- [ ] BMO の任意の代表元の平均は全立方体で共通である。

:::explanation
原子と定数の積分は平均ゼロにより消える。有限原子和が稠密なので同じことが $H^1$ 全体に移る。例えば半直線の指示関数の平均は立方体の位置で変わる。
:::
::::

::::quiz{answer="1.33333333333333" tol="1e-6"}
$w(x)=x^{1/2}$、$Q=(0,1)$ とする。この区間における $\bigl(\operatorname{avg}_Qw\bigr)\bigl(\operatorname{avg}_Qw^{-1}\bigr)$ を求めよ。
:::explanation
$\int_0^1 x^{1/2}dx=2/3$、$\int_0^1x^{-1/2}dx=2$ である。区間の長さは $1$ なので平均もそれぞれ同じ値であり、積は $4/3$ である。
:::
::::

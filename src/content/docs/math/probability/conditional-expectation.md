---
title: "条件付き期待値"
description: "既知の情報の下での平均を測度の導関数として定め、条件付き分布と正規回帰を導く。"
notation:
  - symbol: '\mathbb E[X\mid\mathcal G]'
    meaning: "情報 G に関する条件付き期待値"
  - symbol: 'K(\omega,B)'
    meaning: "正則条件付き分布の核"
  - symbol: '\Sigma_{XX}^{+}'
    meaning: "共分散ブロックの擬逆行列"
---

## この章の目標

- 正の確率をもつ事象への条件付けと、情報全体への条件付けを区別する。
- Radon–Nikodym の定理で条件付き期待値を構成し、計算法を証明する。
- 条件付き分布の核を構成し、条件付き独立性の判定に使う。
- 退化を含む正規ベクトルの条件付き平均と分散を計算する。

## 初等的な条件付き確率と Bayes の公式

事象 $B$ が起こったという情報を得たら、$B$ 内に残った確率を全体で割り直す。$\mathbb P(B)>0$ のとき :term[条件付き確率]{en="conditional probability" yomi="じょうけんつきかくりつ"} を $\mathbb P(A\mid B)=\mathbb P(A\cap B)/\mathbb P(B)$ と置く。この関数は $A$ について確率測度である。分割 $(B_i)$ が互いに素で全体を覆い、各確率が正なら、可算加法性より
$$
\mathbb P(A)=\sum_i\mathbb P(A\mid B_i)\mathbb P(B_i)
$$
である。$\mathbb P(A)>0$ なら両辺を分母に使い、$\mathbb P(B_j\mid A)=\mathbb P(A\mid B_j)\mathbb P(B_j)/\mathbb P(A)$ を得る。これが Bayes の公式である。

:::example[異なる母集団からの抽出]
箱 $B_1$ を確率 $1/3$、箱 $B_2$ を確率 $2/3$ で選ぶ。赤玉を引く確率は箱ごとに $3/4,1/4$ とする。赤の全確率は $1/4+1/6=5/12$ で、赤を見たときの箱 $B_1$ の確率は $(1/4)/(5/12)=3/5$ である。
:::

:::example[観測された個数]
二回の公平で独立な硬貨投げで表の回数を $Y$ とする。$Y=1$ と分かると、同確率の二通りが残るので一回目が表である確率は $1/2$ である。$Y=0,2$ では一回目の結果はそれぞれ裏、表と決まる。従って条件付き確率を観測値 $Y$ の関数として書ける。
:::

### 事前の比と観測の倍率

:::example[二種類の検査を順に使う]
状態 $D$ の事前確率を $1/5$ とする。最初の結果 $A$ は $D$ の下で確率 $3/4$、$D^c$ の下で $1/4$ である。従って $\mathbb P(A)=3/20+1/5=7/20$ で、$\mathbb P(D\mid A)=3/7$ である。次に $A\cap D$ の中で結果 $B$ の確率が $2/3$、$A\cap D^c$ の中では $1/3$ とする。このとき
$$
\mathbb P(D\mid A\cap B)=\frac{(3/20)(2/3)}{(3/20)(2/3)+(1/5)(1/3)}
=\frac35
$$
となる。第二段階で使ったのは $A$ を得た後の群ごとの確率である。$A$ と $B$ が無条件に独立だと置き換える必要はない。
:::

:::example[観測された数と順番の違い]
三回の公平な独立試行で、成功がちょうど二回だったという情報を得た。三つの成功位置の組はすべて元の確率が $1/8$ であり、条件付け後はそれぞれ確率 $1/3$ である。従って最初が成功である確率は $2/3$、さらに最初が成功だと分かった後に第二も成功である確率は $1/2$ になる。同じ報告の中でどの順番まで指定するかにより、分母に残る組が変わる。
:::

### 連続な母数の Bayes 計算

有限の事前分布では、各候補の重みに観測の確率を掛けて規格化する。候補が連続なら同じ操作を密度に対して行い、分母は和から積分に変わる。このとき、観測の分布に密度があるという条件は母数の事前分布が連続であることとは別である。

:::proposition[尤度による事前分布の更新]{#prop-prob-bayes-kernel}
母数 $T$ が確率測度 $\pi$ に従い、$T=t$ を固定したとき観測 $Y$ に共通の $\sigma$ 有限な基準測度 $\lambda$ に関する密度 $L(t,y)\ge0$ があるとする。$L$ は二変数で可測とし、各 $t$ で $\int L(t,y)\lambda(\dd y)=1$ とする。$m(y)=\int L(t,y)\pi(\dd t)$ と置くと、$0<m(y)<\infty$ の点で :term[事後分布]{en="posterior distribution" yomi="じごぶんぷ"} は
$$
\pi_y(\dd t)=\frac{L(t,y)}{m(y)}\pi(\dd t)
$$
である。
:::

:::proof
非負積分を使うと $m$ は可測で、$\int m\,\dd\lambda=1$ である。母数が $A$、観測が $B$ である同時確率は $\int_A\int_B L(t,y)\lambda(\dd y)\pi(\dd t)$ である。積分順序を交換すれば $\int_B \pi_y(A)m(y)\lambda(\dd y)$ となる。$m$ が零または無限の場所は観測の分布で零測度なので、そこで $\pi_y$ に固定の分布を補ってもこの等式は保たれる。核は各 $y$ で確率測度であり、各 $A$ で可測である。従って長方形上の同時確率の計算が、事後法則の定義式を満たす。
:::

:::example[一様事前から成功率を更新する]
成功率 $P$ の事前法則を $(0,1)$ 上の一様分布とする。条件付きで独立な三試行を行い、順序を指定して成功、失敗、成功だったとする。尤度は $p^2(1-p)$ で、周辺確率は $\int_0^1p^2(1-p)\dd p=1/12$ である。事後密度は $12p^2(1-p)$ だから Beta の形状は $(3,2)$ である。次回の成功の条件付き確率は $\int_0^1p\cdot12p^2(1-p)\dd p=3/5$ となる。もし順序を忘れて成功二回だけを条件とすれば、尤度は $3p^2(1-p)$、周辺確率は $1/4$ となる。三倍が規格化で消えるため事後分布は同じだが、観測した事象自体の確率は異なる。
:::

:::example[Poisson 個数で Gamma の母数を更新する]
強度 $T$ に形状 $2$、率 $3$ の Gamma 事前分布を置き、時間 $2$ に得る回数 $K$ は条件付きで $\operatorname{Pois}(2T)$ とする。$K=4$ を観測したとき、事前密度に比例する $t e^{-3t}$ と、尤度に比例する $t^4 e^{-2t}$ を掛けると $t^5e^{-5t}$ となる。従って事後は形状 $6$、率 $5$ の Gamma 分布で、平均強度は $6/5$ である。新たに時間 $1$ の回数 $J$ を $T$ の下で独立に取ると、$\mathbb E[J\mid K=4]=6/5$、$\mathbb P(J=0\mid K=4)=\mathbb E[e^{-T}\mid K=4]=(5/6)^6$ である。後者では密度を指数で重み付けして積分しており、事後平均強度をそのまま指数に入れた $e^{-6/5}$ とは異なる。
:::

## Radon–Nikodym の定理による定義と存在

:::definition[条件付き期待値]{#def-prob-conditional-expectation}
$X\in L^1(\Omega,\mathcal F,\mathbb P)$ とし、$\mathcal G\subset\mathcal F$ を部分 $\sigma$-加法族とする。$\mathcal G$ 可測な可積分変数 $Z$ が全 $A\in\mathcal G$ で
$$
\int_A Z\,\dd\mathbb P=\int_A X\,\dd\mathbb P \label{eq-prob-conditional-integrals}
$$
を満たすとき、$Z$ を :term[条件付き期待値]{en="conditional expectation" yomi="じょうけんつききたいち"} $\mathbb E[X\mid\mathcal G]$ の一つの版という。事象 $B$ に対し $\mathbb P(B\mid\mathcal G)=\mathbb E[\mathbf1_B\mid\mathcal G]$ とする。観測 $Y$ が生成する $\sigma(Y)$ を使う場合は $\mathbb E[X\mid Y]$、$\mathbb P(B\mid Y)$ とも書く。
:::

:::theorem[存在と概確実な一意性]{#thm-prob-conditional-existence}
:ref[def-prob-conditional-expectation] を満たす版は存在し、零事象上での違いを除いて一意である。
:::

:::proof
$\mathcal G$ に制限した符号付き測度 $\nu(A)=\int_A X\,\dd\mathbb P$ は有限で $\mathbb P|_{\mathcal G}$ に絶対連続である。従って :ref[thm-meas-rn] で導関数 $Z=\dd\nu/\dd(\mathbb P|_{\mathcal G})$ が存在し、:eqref[eq-prob-conditional-integrals] を満たす。二つの版の差が正である集合を $A$ に取ると、差の積分が零なのでその集合は零事象である。負の集合についても同様である。
:::

$\mathcal G$ が可算分割 $(B_i)$ で生成される場合、正の確率をもつ $B_i$ 上で定数 $\mathbb E[X\mathbf1_{B_i}]/\mathbb P(B_i)$ を使えば、各セルの積分が一致する。零確率の $B_i$ には有限な定数を置けばよい。特に自明な $\mathcal G=\{\varnothing,\Omega\}$ なら定数 $\mathbb EX$、$\mathcal G=\mathcal F$ なら $X$ 自身になる。非負 $X$ に拡張するときは $\mathbb E[X\wedge m\mid\mathcal G]$ の単調極限を取れば、単調収束で :eqref[eq-prob-conditional-integrals] が保たれる。

### 有限分割の上で積分条件を解く

:::proposition[可算分割に関する条件付き期待値]{#prop-prob-conditional-partition}
$C_i$ を可算可測分割とし $\mathcal G=\sigma(C_i:i\ge1)$、$X\in L^1$ とする。正の確率をもつ $C_i$ 上では、$\mathbb E[X\mid\mathcal G]$ は定数 $\mathbb E[X\mathbf1_{C_i}]/\mathbb P(C_i)$ である。零確率のセルには任意の有限定数を置ける。
:::

:::proof
候補 $Z$ は $\mathcal G$-可測である。$\int|Z|\,\dd\mathbb P=\sum_{i:\mathbb P(C_i)>0}|\int_{C_i}X\,\dd\mathbb P|\le \mathbb E|X|$ なので可積分である。各セルで $\int_{C_i}Z=\int_{C_i}X$ が成り立つ。$\mathcal G$ の事象はセルの可算和であり、両積分は可積分性により和に交換できる。従って $\mathcal G$ の全事象で一致し、定義を満たす。
:::

:::example[不等確率のセルで平均する]
$X=1,2,5$ の確率を $1/2,1/4,1/4$ とする。観測は $X\le2$ か否かだけとする。前者では積分が $1/2+2/4=1$、確率は $3/4$ なので予測値は $4/3$ である。後者では $5$ が確定する。$4/3$ と $5$ を元のセル確率で平均すると $9/4$ になり、全期待値と一致する。セル内で $1$ と $2$ を同じ重みで平均して $3/2$ とすると積分条件を満たさない。
:::

### テスト関数で条件付き期待値を同定する

定義は事象の指示関数との積分を比較する。実際の計算では観測と連続な関数を掛けた積分の方が求めやすいことがある。実数値観測では、そのような試験関数だけで候補を確認できる。

:::proposition[観測の関数との積分による同定]{#prop-prob-conditional-test-functions}
$X$ が実数値、$Y\in L^1$、$g$ が実 Borel 関数で $g(X)\in L^1$ とする。全ての有界連続関数 $h$ に対して $\mathbb E[h(X)Y]=\mathbb E[h(X)g(X)]$ なら $g(X)=\mathbb E[Y\mid X]$ 概確実である。逆も成り立つ。
:::

:::proof
二式の差 $Y-g(X)$ に $h(X)$ を掛ける。有界連続な $h_m(x)=\max(0,1-m(x-a)_+)$ は半直線 $(-\infty,a]$ の指示関数へ収束する。差の絶対値で支配して積分の極限を取れば、全半直線の逆像で差の積分は零になる。定数 $h=1$ で全体の積分も零である。差の積分が零になる観測値の Borel 集合の族は $\lambda$ 系なので、半直線から全 Borel 集合へ拡張できる。よって $g(X)$ は積分条件を満たす。逆向きは、試験関数を指示関数の線形結合で一様に近似して定義から従う。
:::

:::example[対称性を使うと二変数の積分が消える]
$X,Y$ が独立で同じ可積分分布に従い、$S=X+Y$ とする。任意の有界連続な $h$ に対して $\mathbb E[h(S)X]=\mathbb E[h(S)Y]$ であり、和を取ると各項は $\tfrac12\mathbb E[h(S)S]$ である。従って $\mathbb E[X\mid S]=S/2$ となる。例えば率が等しい独立指数変数なら、$s\mapsto s/2$ という連続な版では、合計時間 $s=4$ に対応する条件付き平均は $2$ である。和の値が固定で零確率の事象でも、観測に関して積分する等式ならそのまま使える。
:::

### 非負変数と符号付き変数の定義の境目

:::proposition[非負変数に対する拡張条件付き期待値]
$X$ を非負可測とし、$X$ の期待値が無限であることを許す。このとき $Z:\Omega\to[0,\infty]$ で $\mathcal G$-可測かつ $\int_A Z=\int_A X$ が全 $A\in\mathcal G$ で成り立つものが存在し、概確実に一意である。
:::

:::proof
$X_n=X\wedge n$ とする。非負な可積分変数の条件付き期待値が負になる集合で積分条件を使うと、その集合は零事象である。これを各 $X_n$ および $X_{n+1}-X_n$ に適用すると、$\mathbb E[X_n\mid\mathcal G]$ は $\mathcal G$-可測な共通の零集合の外で非負かつ非減少にできる。その集合で全ての版を零に補えば、極限 $Z$ は拡張実数値で可測である。任意の $A\in\mathcal G$ で二辺に単調収束を使うと積分等式を得る。

一意性のため別の $W$ を取る。$Z>W$ では $W<\infty$ であり、ある正整数 $m$ について $W\le m$、$Z\ge W+1/m$ となる。この二条件を満たす事象 $B_m$ は $\mathcal G$ に属し、$\int_{B_m} W\le m$ である。積分等式と $Z\ge W+1/m$ から $\mathbb P(B_m)=0$ が従う。可算和で $\{Z>W\}$ は零であり、反対の不等式にも同じ議論を使う。
:::

:::example[条件を分けると片方の期待値だけ無限になる]
公平なビット $A$ を先に選び、$A=0$ なら $X=1$、$A=1$ なら独立な $U\sim\operatorname{Unif}(0,1)$ から $X=1/U$ とする。$X$ は非負で全期待値は無限だが、$\mathbb E[X\mid A]$ は $A=0$ で $1$、$A=1$ で $\infty$ である。例えば $n\ge1$ では後者のセルで $\mathbb E[X\wedge n\mid A=1]=1+\log n$ であり、切断を外すと発散する。符号付き変数で、条件付けた正部分と負部分の両方が無限になると、差 $\infty-\infty$ は定義できない。例えば独立公平符号を $1/U$ に掛けて自明な情報で条件付けると、その問題が生じる。
:::

### 条件付き期待値を測度を変える前後で比較する

:::proposition[密度で測度を変えたときの条件付き期待値]{#prop-prob-conditional-change-of-measure}
$\mathbb Q\ll\mathbb P$ を確率測度、$L=\dd\mathbb Q/\dd\mathbb P$、$X\in L^1(\mathbb Q)$ とする。$\mathbb Q$ の下での条件付き期待値は、$\mathbb Q$-概確実に
$$
\mathbb E_{\mathbb Q}[X\mid\mathcal G]
=\frac{\mathbb E_{\mathbb P}[LX\mid\mathcal G]}
{\mathbb E_{\mathbb P}[L\mid\mathcal G]}
$$
である。分母が零の場所では任意の有限値を補える。
:::

:::proof
$D=\mathbb E_{\mathbb P}[L\mid\mathcal G]$、$N=\mathbb E_{\mathbb P}[LX\mid\mathcal G]$ と置く。$D$ は条件付き平均なので $\mathbb Q(D=0)=\mathbb E_{\mathbb P}[L\mathbf1_{\{D=0\}}]=\mathbb E_{\mathbb P}[D\mathbf1_{\{D=0\}}]=0$ である。さらに $\{D=0\}$ 上では $\mathbb P$ の下で $L=0$ 概確実だから、$N=0$ も同じ場所で概確実に成り立つ。$Z=(N/D)\mathbf1_{\{D>0\}}$（分母零では $Z=0$ と定める）とすると、$Z$ は $\mathcal G$-可測で、
$$
\mathbb E_{\mathbb Q}|Z|=\mathbb E_{\mathbb P}[L|Z|]
=\mathbb E_{\mathbb P}[D|Z|]\le\mathbb E_{\mathbb P}|LX|<\infty
$$
である。持ち出しには非負変数の単調近似を使う。同じ計算を任意の $A\in\mathcal G$ 上で符号付きに適用すると $\mathbb E_{\mathbb Q}[Z\mathbf1_A]=\mathbb E_{\mathbb P}[N\mathbf1_A]=\mathbb E_{\mathbb Q}[X\mathbf1_A]$ である。従って $\mathbb Q$ 下の定義を満たす。
:::

:::example[不公平なサイコロで偶奇の予測を更新する]
公平な測度 $\mathbb P$ で出目を $X$ とし、$L=2X/7$ で測度を変える。偶数セルでは公平な条件付き平均 $\mathbb E_{\mathbb P}[X\mid X\text{ は偶数}]=4$、条件付き二次積率は $(4+16+36)/3=56/3$ である。新しい条件付き平均は比 $(56/3)/4=14/3$ となる。奇数セルでは同じ計算で $(1+9+25)/(1+3+5)=35/9$ となる。公平な予測 $4,3$ に一定の補正を加えれば済むわけではなく、各セルの中で重みを規格化し直す。
:::

## 条件付き期待値の性質

:::proposition[計算規則]{#prop-prob-conditional-rules}
可積分変数について条件付き期待値は線形、単調で、$\mathbb E\mathbb E[X\mid\mathcal G]=\mathbb EX$ である。$\mathcal H\subset\mathcal G$ なら
$$
\mathbb E[\mathbb E[X\mid\mathcal G]\mid\mathcal H]=\mathbb E[X\mid\mathcal H].
$$
有界な $\mathcal G$ 可測 $Y$ について $\mathbb E[YX\mid\mathcal G]=Y\mathbb E[X\mid\mathcal G]$ である。$X$ が $\mathcal G$ と独立なら $\mathbb E[X\mid\mathcal G]=\mathbb EX$ である。
:::

:::proof
線形性と全期待値は :eqref[eq-prob-conditional-integrals] から従う。$X\ge0$ なら候補が負の集合に積分すると負になってはならないため、候補は非負である。この事実を差に適用すると単調性が従う。塔の公式では $\mathcal H$ の集合はすべて $\mathcal G$ にも属するため、積分を二度置き換えられる。既知の因子の公式は、$Y$ が指示関数のときに集合の交わりで積分すればよく、単関数と有界近似で一般化する。独立なら $\mathbb E[X\mathbf1_A]=\mathbb EX\,\mathbb P(A)$ で、定数が条件を満たす。一般の $\mathcal G$ 可測 $Y$ でも $YX$ が可積分なら、有界切断で近似して同じ式が成り立つ。実際、$|Y|\mathbb E[|X|\mid\mathcal G]$ の期待値は単調収束により $\mathbb E|YX|$ なので、右辺の切断にも優収束が使える。
:::

条件付き Jensen の不等式 $\phi(\mathbb E[X\mid\mathcal G])\le\mathbb E[\phi(X)\mid\mathcal G]$ は、有限凸関数 $\phi:\R\to\R$ を可算個の支えの直線の上限で表せば得られる。有理点ごとに一つ支えを取り、凸関数の連続性とコンパクトな区間上での勾配の局所有界性を使えば、任意の実点で有理点から選んだ直線の値が $\phi$ に近づく。各直線ごとに単調性と線形性を適用し、可算個の例外集合を一つにまとめて上限を取る。一次の支えがあるため、負部分は可積分である。特に $p\ge1$ では $\|\mathbb E[X\mid\mathcal G]\|_p\le\|X\|_p$ である。また $0\le X_n\uparrow X$ のときは、条件付き期待値の単調極限について各 $A$ 上で単調収束を用いると、条件付き単調収束が得られる。これを $\inf_{k\ge n}X_k$ に適用して条件付き Fatou を得る。$|X_n|\le H\in L^1$ かつ $X_n\to X$ 概確実なら、$H\pm X_n$ に Fatou を適用して条件付き優収束が従う。

:::example[情報を段階的に粗くする]
公平なサイコロの出目 $X$ について、偶奇を知ると偶数側の平均は $(2+4+6)/3=4$、奇数側は $3$ である。偶奇の情報も捨てると $(4+3)/2=7/2$ となり、元の平均と一致する。
:::

### 情報を変えても尾を一様に制御する

:::proposition[一つの変数を条件付けた族の一様可積分性]{#prop-prob-conditional-ui}
$X\in L^1$ を固定すると、$\mathcal G$ を任意の部分 $\sigma$-加法族に動かした $\{\mathbb E[X\mid\mathcal G]\}$ は一様可積分である。
:::

:::proof
$Z=\mathbb E[X\mid\mathcal G]$ とすれば $|Z|\le\mathbb E[|X|\mid\mathcal G]$ である。$A=\{|Z|>K\}\in\mathcal G$ 上で積分すると $\mathbb E[|Z|\mathbf1_A]\le\mathbb E[|X|\mathbf1_A]$ となる。また $\mathbb P(A)\le\mathbb E|X|/K$ である。右辺を $|X|\le M$ と他に分けると $M\mathbb E|X|/K+\mathbb E[|X|\mathbf1_{\{|X|>M\}}]$ 以下である。先に $M$ を大きく選び、次に $K$ を大きく選べば、$\mathcal G$ に依らず任意に小さくなる。
:::

:::example[塔の公式には包含が要る]
公平な独立ビット $U,V$ で $X=U+2V$ とする。$\mathbb E[X\mid U]=U+1$、$\mathbb E[X\mid V]=1/2+2V$ である。前者を $V$ で再度条件付けると $3/2$ となり、後者とは一般に等しくない。$\sigma(U)$ と $\sigma(V)$ のどちらも他方を含まない。一方 $\sigma(U,V)$ で条件付けてから $\sigma(V)$ に移る場合は、最初の値が $X$ そのものなので塔の等式が成立する。
:::

### 条件付きの Hölder 不等式

:::proposition[条件付き Hölder 不等式]{#prop-prob-conditional-holder}
$p,q>1$、$1/p+1/q=1$ とする。$X\in L^p$、$Y\in L^q$ に対して
$$
\mathbb E[|XY|\mid\mathcal G]
\le \mathbb E[|X|^p\mid\mathcal G]^{1/p}
       \mathbb E[|Y|^q\mid\mathcal G]^{1/q}
$$
が概確実に成り立つ。
:::

:::proof
$A=\mathbb E[|X|^p\mid\mathcal G]^{1/p}$、$B=\mathbb E[|Y|^q\mid\mathcal G]^{1/q}$ とする。$\eps>0$ を固定し、Young の点ごとの不等式に $|X|/(A+\eps)$、$|Y|/(B+\eps)$ を代入する。分母の逆数は有界な $\mathcal G$-可測変数なので条件付き期待値の外へ出せる。その結果
$$
\frac{\mathbb E[|XY|\mid\mathcal G]}{(A+\eps)(B+\eps)}
\le\frac{A^p}{p(A+\eps)^p}+\frac{B^q}{q(B+\eps)^q}\le1
$$
となる。通常の Hölder により $|XY|$ は可積分である。$\eps\downarrow0$ を可算列で取れば、$A=0$ または $B=0$ の事象も含めて結論が従う。
:::

:::example[観測ごとに異なる Cauchy–Schwarz 上界]
公平ビット $A$ を観測する。$A=0$ で $(X,Y)=(Z,2Z)$、$A=1$ で $(X,Y)=(Z,W)$ とし、$Z,W$ は互いに独立な公平符号であり、組として $A$ とも独立である。各セルで $\mathbb E[|XY|\mid A]$ は $2,1$、条件付き二次積率の平方根の積も $2,1$ である。しかし $\mathbb E[XY\mid A]$ は $2,0$ であり、後者のセルでは符号が打ち消し合う。絶対値を積分の前に取るか後に取るかで不等式の鋭さが変わる。
:::

### $L^p$ での縮小と測定前後の誤差

:::proposition[条件付き期待値の $L^p$ 縮小性]
$1\le p\le\infty$ とする。$X,Y\in L^p$ に対し
$$
\|\mathbb E[X\mid\mathcal G]-\mathbb E[Y\mid\mathcal G]\|_p
\le\|X-Y\|_p
$$
である。特に同じ $\mathcal G$ で条件を取るとき、$X_n\to X$ の $L^p$ 収束は条件付き期待値にも移る。
:::

:::proof
線形性で $Y=0$ の場合に帰着する。$p=1$ は条件付き絶対値の評価を全積分する。$1<p<\infty$ では $\phi(x)=|x|^p$ に対する条件付き Jensen で $|\mathbb E[X\mid\mathcal G]|^p\le\mathbb E[|X|^p\mid\mathcal G]$ であり、積分して $p$ 乗根を取る。$p=\infty$ では $-\|X\|_\infty\le X\le\|X\|_\infty$ から、単調性と定数の持ち出しで同じ有界区間に条件付き平均が入る。差に同じ不等式を使えば収束も得られる。
:::

:::example[条件付けてから切断を外す]
$X$ が平均 $1$ の指数変数とし、$\mathcal G=\sigma(\mathbf1_{\{X>1\}})$ とする。$X_K=X\wedge K$（$K>1$）と置く。平均の誤差は
$$
\mathbb E|X-X_K|=\int_K^\infty e^{-t}\,\dd t=e^{-K}
$$
である。従って条件付き平均の $L^1$ 誤差も高々 $e^{-K}$ である。実際、$X\le1$ のセルでは切断は影響せず、$X>1$ のセルでは $\mathbb P(X-1>t\mid X>1)=e^{-t}$（$t\ge0$）なので、切断による平均の差は $e^{-(K-1)}$ である。セルの確率 $e^{-1}$ を掛けると $e^{-K}$ となり、この例では $L^1$ の縮小評価に等号が成立する。差が非負であるため、条件付き平均で符号が打ち消されないことが理由である。
:::

:::example[$L^2$ で等号になる観測と情報を捨てる観測]
$X$ が $-1,0,1$ を等確率で取るとする。$\mathcal G=\sigma(|X|)$ では $|X|=1$ のセルで $-1,1$ が等確率、零のセルでは $0$ なので条件付き平均はいつも零である。従って $\|X\|_2^2=2/3$ が予測から全て失われる。一方 $\sigma(X)$ を使えば条件付き平均は $X$ 自身でありノルムは変わらない。一般に二乗平均の等号は、射影の誤差の二乗積分が零であることと同値で、観測で $X$ が概確実に復元できるときに起こる。
:::

### 過去から予測した確率の和で判定する

:::proposition[条件付き Borel–Cantelli の補題]{#prop-prob-conditional-borel-cantelli}
$\mathcal F_0\subset\mathcal F_1\subset\cdots$ を増加する $\sigma$-加法族、$A_k\in\mathcal F_k$ とし、$p_k=\mathbb P(A_k\mid\mathcal F_{k-1})$ と置く。確率 $1$ の集合上で、
$$
\sum_{k\ge1}\mathbf1_{A_k}=\infty
\quad\Longleftrightarrow\quad
\sum_{k\ge1}p_k=\infty
$$
が成り立つ。
:::

:::proof
まず $\sum_k p_k<\infty$ の部分を考える。$C_k=\sum_{j=1}^k p_j$ は $\mathcal F_{k-1}$-可測である。整数 $M\ge1$ ごとに
$$
\mathbb E\sum_k \mathbf1_{\{C_k\le M\}}\mathbf1_{A_k}
=\mathbb E\sum_k \mathbf1_{\{C_k\le M\}}p_k\le M
$$
が成り立つ。最初の等式には非負和と条件付き期待値を使い、最後は $C_k$ が非減少であることによる。従って各 $M$ で左の和は概確実に有限である。$C_\infty<\infty$ の標本では十分大きい整数 $M$ が全ての $C_k$ を上から抑えるので、実際の発生数も有限である。

次に固定した $m$ 以降に一回も発生しない事象を考える。$n\ge m$ で
$$
L_n=\mathbf1_{\cap_{k=m}^{n}A_k^c}\exp\left(\sum_{k=m}^{n}p_k\right)
$$
とし、$L_{m-1}=1$ とする。$\mathcal F_{n-1}$ で条件付けると $L_{n-1}e^{p_n}(1-p_n)\le L_{n-1}$ であり、従って各 $\mathbb E L_n\le1$ である。もし $m$ 以降発生せず、かつ $\sum p_k=\infty$ なら $L_n\to\infty$ となる。Fatou で $\mathbb E\liminf_nL_n\le1$ だから、この事象の確率は零である。発生回数が有限ならある整数 $m$ 以降に発生しないので、可算和を取ると第二の向きが従う。
:::

:::example[小さな確率を、過去の状態に応じて選ぶ]
最初に公平な $C$ を観測し、その後の独立一様数 $U_k$ で
$$
A_k=
\begin{cases}
\{U_k\le1/(k+1)\}& C=1,\\
\{U_k\le1/(k+1)^2\}& C=0
\end{cases}
$$
を判定する。過去には $C$ と使用済みの $U_j$ が入っており、$p_k$ は第一の枝で $1/(k+1)$、第二の枝で $1/(k+1)^2$ である。従って第一の枝で無限回、第二の枝で有限回となり、全体では無限回の確率が $1/2$ である。周辺確率の和は両枝を平均した値の和なので発散するが、条件付き確率の和は枝ごとに異なる。共通の観測によって強い依存がある場合に、その違いを表せる。
:::

### 条件付きの分散を足して正規極限を求める

:::proposition[有界なマルチンゲール差配列の中心極限定理]{#prop-prob-bounded-martingale-clt}
各行で $\mathcal F_{n,0}\subset\cdots\subset\mathcal F_{n,k_n}$ とする。$D_{n,k}$ は $\mathcal F_{n,k}$-可測で $\mathbb E[D_{n,k}\mid\mathcal F_{n,k-1}]=0$ を満たす。このような列をマルチンゲール差という。全 $k$ で $|D_{n,k}|\le\delta_n$、$\delta_n\to0$ とする。
$$
V_n=\sum_k\mathbb E[D_{n,k}^2\mid\mathcal F_{n,k-1}]
$$
が定数 $v\ge0$ に確率収束し、全行で $V_n\le C$ と一様有界なら、$\sum_kD_{n,k}\Rightarrow N(0,v)$ である。
:::

:::proof
一行で $S_k=\sum_{j\le k}D_{n,j}$、$V_k=\sum_{j\le k}v_j$、$v_j=\mathbb E[D_{n,j}^2\mid\mathcal F_{n,j-1}]$ と略記する。固定した $t$ について条件付き Taylor 展開を行えば、
$$
\mathbb E[e^{itD_{n,k}}\mid\mathcal F_{n,k-1}]
=1-\tfrac12t^2v_k+r_k,\qquad |r_k|\le |t|^3\delta_n v_k/6.
$$
条件付き一次項は零であり、$|D_{n,k}|^3\le\delta_n D_{n,k}^2$ を使った。$v_k$ は直前に可測なので、$Z_k=\exp(itS_k+\tfrac12t^2V_k)$ を用いると
$$
\mathbb E[Z_k\mid\mathcal F_{n,k-1}]
=Z_{k-1}e^{t^2v_k/2}(1-t^2v_k/2+r_k).
$$
括弧付きの因子と $1$ の差は $C_t(v_k^2+\delta_n v_k)$ 以下である。これは $0\le v_k\le\delta_n^2$ と指数の二次余りから得られる。$|Z_{k-1}|\le e^{t^2C/2}$、$\sum_k v_k\le C$ なので、条件付き期待値の差を足して終端で $|\mathbb E Z_{k_n}-1|\le C_t'(\delta_n^2+\delta_n)\to0$ となる。

さらに $V_n\to v$ は確率収束し一様有界であるから、$\mathbb E|\exp(-t^2V_n/2)-\exp(-t^2v/2)|\to0$ である。$|Z_{k_n}|$ の一様上界を使うと $\mathbb Ee^{itS_{k_n}}=\mathbb E[Z_{k_n}e^{-t^2V_n/2}]\to e^{-t^2v/2}$ となる。Lévy の連続性定理で結論が従う。
:::

:::example[先に観測した分散の型が途中で入れ替わる]
公平ビット $B$ を最初に観測し、これと独立な公平符号 $\xi_k$ を順に観測する。長さ $2n$ の行では、前半の $n$ 項を $a_B\xi_k/\sqrt n$、後半の $n$ 項を $b_B\xi_k/\sqrt n$ とする。$(a_0,b_0)=(1,2)$、$(a_1,b_1)=(2,1)$ とする。過去を固定した各符号の平均は零なので、これはマルチンゲール差の配列である。各項は $2/\sqrt n$ 以下の絶対値を持ち、条件付き分散和は $a_B^2+b_B^2=5$ である。従って極限は $N(0,5)$ となる。前半の二項の絶対値は共に $B$ によって決まり、無条件では独立でない。
:::

## L² 射影としての条件付き期待値

:::theorem[最小二乗予測]{#thm-prob-conditional-projection}
$X\in L^2$ に対し $Z=\mathbb E[X\mid\mathcal G]$ は閉部分空間 $L^2(\mathcal G)$ への直交射影であり、任意の $Y\in L^2(\mathcal G)$ で
$$
\mathbb E(X-Y)^2=\mathbb E(X-Z)^2+\mathbb E(Z-Y)^2
$$
が成り立つ。
:::

:::proof
条件付き Jensen で $Z\in L^2$ となる。$L^2(\mathcal G)$ の列が $L^2$ 収束すると概収束する部分列を取れるため、可測な極限の版が存在し、部分空間は閉じている。任意の $Y$ に対して既知の因子の公式が使える。実際、$XY$ と $ZY$ は Cauchy–Schwarz で可積分である。従って $\mathbb E[(X-Z)Y]=0$ である。二乗を展開すると交差項が零になり所望の分解を得る。
:::

条件付き分散を $\operatorname{Var}(X\mid\mathcal G)=\mathbb E[(X-Z)^2\mid\mathcal G]$ と置く。最小二乗式で $Y=\mathbb EX$ と取ると全分散の公式 $\operatorname{Var}(X)=\mathbb E\operatorname{Var}(X\mid\mathcal G)+\operatorname{Var}(Z)$ になる。例えば前節のサイコロでは $\operatorname{Var}(X)=35/12$、偶奇による予測値の分散は $1/4$ であり、残差分散の平均は $8/3$ である。

### 予測に残った誤差を共分散で調べる

:::proposition[全共分散の公式]
$X,Y\in L^2$、$M=\mathbb E[X\mid\mathcal G]$、$N=\mathbb E[Y\mid\mathcal G]$ とすると
$$
\operatorname{Cov}(X,Y)=
\mathbb E[(X-M)(Y-N)]+\operatorname{Cov}(M,N).
$$
:::

:::proof
$X-\mathbb EX=(X-M)+(M-\mathbb EX)$ とし、$Y$ も同様に分ける。残差と $\mathcal G$-可測な二乗可積分変数との積分は射影の直交性により零である。従って四つの積のうち、二つの残差の積と二つの予測の中心化した積だけが残る。すべての積は Cauchy–Schwarz で可積分である。
:::

:::example[粗い観測と最良の定数予測を比較する]
公平なサイコロの出目 $X$ に対し、$X\le3$ か否かを観測する。条件付き平均は $2,5$ であり、両セルの確率は $1/2$ なので予測の分散は $(3/2)^2=9/4$ である。元の分散 $35/12$ から引くと、最小二乗誤差は $2/3$ である。実際、どちらのセルでも残差は $-1,0,1$ を等確率で取る。一定値 $7/2$ で予測する場合の誤差 $35/12$ に比べ、セルの観測が区別した分だけ誤差が減る。
:::

### 追加の観測がどれだけ予測を改善するか

:::proposition[二段階の予測誤差の差]{#prop-prob-conditional-error-decomposition}
$X\in L^2$、$\mathcal H\subset\mathcal G$ とする。$H=\mathbb E[X\mid\mathcal H]$、$G=\mathbb E[X\mid\mathcal G]$ なら
$$
\mathbb E(X-H)^2=\mathbb E(X-G)^2+\mathbb E(G-H)^2.
$$
また条件付き分散についても
$$
\operatorname{Var}(X\mid\mathcal H)=\mathbb E[\operatorname{Var}(X\mid\mathcal G)\mid\mathcal H]+\mathbb E[(G-H)^2\mid\mathcal H]
$$
が成り立つ。
:::

:::proof
分解 $X-H=(X-G)+(G-H)$ の交差項を考える。$G-H$ は $\mathcal G$-可測で二乗可積分、$X-G$ は $\mathcal G$ で条件付き平均が零である。有界に切断した $G-H$ を用いて条件付き積分を計算し、Cauchy–Schwarz による $L^1$ 収束で切断を外すと、交差項の $\mathcal G$ 条件付き平均は零である。従って $\mathcal H$ で条件付けても塔の公式で零となる。二乗展開の残り二項に塔の公式を使い条件付き分散の式を得る。さらに全期待値を取れば最初の等式になる。
:::

:::example[観測を取得する順番で改善量を比べる]
公平な独立ビット $U,V$ を使い、$X=U+2V$ とする。観測なしの分散は $5/4$ である。$U$ を先に読むと予測は $U+1$、誤差の二乗平均は $1$ となり、改善量は $1/4$ である。さらに $V$ を読むと誤差は零となり、改善量は $1$ である。$V$ を先に読めば先に $1$、後に $1/4$ 改善する。どちらの順序でも合計が $5/4$ であるのは、最初と最後の情報が同じだからである。
:::

### 条件付き最小二乗を解いて別の損失と比べる

:::proposition[条件付き二乗損失の最小値]
$X\in L^2$ とし、$Z=\mathbb E[X\mid\mathcal G]$ とする。任意の $\mathcal G$-可測な $Y\in L^2$ について、概確実に
$$
\mathbb E[(X-Y)^2\mid\mathcal G]
=\operatorname{Var}(X\mid\mathcal G)+(Z-Y)^2
$$
である。
:::

:::proof
$X-Y=(X-Z)+(Z-Y)$ を二乗する。$Z-Y$ の有界切断は $\mathcal G$-可測なので、交差項の条件付き期待値では外へ出せる。$\mathbb E[X-Z\mid\mathcal G]=0$ だから、切断した交差項は条件付きで零である。切断を外したときの積の $L^1$ 誤差は Cauchy–Schwarz により零へ向かい、条件付き期待値も $L^1$ 縮小性で収束する。よって交差項は零である。残る二項のうち $(Z-Y)^2$ は $\mathcal G$-可測なのでそのまま取り出せる。
:::

:::example[絶対誤差では中央値が予測値になる]
$X$ は $0$ を確率 $3/4$、$4$ を確率 $1/4$ で取り、観測はないとする。二乗誤差の最小点は平均 $1$ で、最小誤差は $3$ である。絶対誤差では $a\in[0,4]$ に対して $\mathbb E|X-a|=(3/4)a+(1/4)(4-a)=1+a/2$ なので、最小点は $0$、最小誤差は $1$ になる。区間の外では端点へ寄せれば両方の距離が減るので最小にはならない。もし観測によって二点の条件付き確率が逆転するセルがあれば、そこでは絶対誤差の最小点は $4$ に移る。一方、条件付き二乗誤差の最小点はセルごとの平均であり、両端の値を滑らかに補間する。分布の条件付けと、どの損失を最小化するかは別の計算である。
:::

## 正則条件付き分布と分解定理

各 $B$ について条件付き確率の版を別々に取るだけでは、一つの標本 $\omega$ で全 $B$ の値が確率測度になるとは限らない。実直線では有理閾値だけを同時に選んでから、分布関数を復元できる。

:::definition[確率核]{#def-prob-conditional-kernel}
$K(\omega,B)$ が各 $\omega$ で $B$ について確率測度を定め、各 Borel 集合 $B$ で $\omega$ について $\mathcal G$ 可測であれば、$K$ を :term[確率核]{en="probability kernel" yomi="かくりつかく"} という。各 $B$ で $K(\cdot,B)$ が $\mathbb P(Y\in B\mid\mathcal G)$ の版でもあるとき、$Y$ の :term[正則条件付き分布]{en="regular conditional distribution" yomi="せいそくじょうけんつきぶんぷ"} という。
:::

:::theorem[実数値の正則条件付き分布と分解]{#thm-prob-regular-conditional}
実数値 $Y$ の正則条件付き分布が存在し、共通の零事象を除いて一意である。$\mathcal G=\sigma(X)$ で $X$ が実数値なら、$K(\omega,B)=k(X(\omega),B)$ と書ける Borel 確率核 $k$ を選べる。このとき $(X,Y)$ の同時分布は $\mu_X(\dd x)k(x,\dd y)$ である。
:::

:::proof
有理数 $q$ に対し $F_\omega(q)$ を $\mathbb P(Y\le q\mid\mathcal G)$ の版とする。単調性、整数閾値を両端に送った極限 $0,1$、有理 $q$ で $F_\omega(q)=\lim_{m\to\infty}F_\omega(q+1/m)$ はすべて概確実に成り立つ。これは前節の収束則を適用した結果であり、全有理点について可算個の例外しかない。よって共通の可測零集合を除いた上で、$F_\omega(x)=\inf_{q>x}F_\omega(q)$ は分布関数である。例外上では $\delta_0$ の分布関数を使う。:ref[thm-prob-cdf-characterization] により各 $\omega$ で確率測度を定められる。有理半直線では希望する性質が成立する。可測性と積分等式を満たす集合の族は $\lambda$ 系であり、:ref[thm-prob-pi-lambda] により全 Borel 集合に拡張できる。一意性は有理閾値で同時に一致することから従う。

:ref[lem-prob-doob-dynkin] により、各有理 $q$ について $F_\omega(q)=f_q(X(\omega))$ となる Borel 関数 $f_q$ を取れる。上の可算条件を満たさない $x$ の Borel 集合には $\mu_X$ が質量を置かない。その集合で $\delta_0$ を使い、それ以外で同じ構成を行うと $k$ が得られる。長方形 $A\times B$ について同時確率が $\int_A k(x,B)\mu_X(\dd x)$ になることは定義式から従う。生成族上の一意性で同時分布全体に拡張できる。
:::

一般に、被条件付け変数の値域が標準 Borel 空間、すなわち完備可分距離空間の Borel 部分集合と可測同型な空間である場合も扱える。この空間は実直線のある Borel 部分集合 $S$ と可測同型である。この特徴付けの証明は記述集合論に属するため :cite[kallenberg-fmp] を参照する。変数を $S$ に写し、実数値の核を作ると、$K(\omega,S)=1$ が概確実に成り立つ。これが成り立たない可測零集合では $S$ の一点に質量を置く。こうして全点で $S$ に集中する核を逆写像で戻せばよい。

例えば三角形 $0<y<x<1$ 上の同時密度 $2$ では $X$ の周辺密度は $2x$ であり、$0<x<1$ における $Y$ の条件付き密度は $x^{-1}\mathbf1_{(0,x)}(y)$ である。核はこの範囲の各 $x$ で正規化される。端点 $x=0$ での値は分解を変えない。

### 同時密度を核に変える

:::proposition[同時密度の商による分解]
$(X,Y)$ の同時密度を Borel 版 $f(x,y)$、$m(x)=\int f(x,y)\,\dd y$ とする。$0<m(x)<\infty$ の点では $k(x,\dd y)=f(x,y)\,\dd y/m(x)$ とし、他では一つの固定した点質量を置けば、$Y$ の $X=x$ での正則条件付き分布となる。
:::

:::proof
非負積分の可測性により $m$ と、任意の Borel 集合 $B$ に対する $\int_B f(x,y)\,\dd y$ は可測である。候補の核は各 $x$ で確率測度である。$\int m(x)\dd x=1$ なので $m=\infty$ は Lebesgue 零集合上でしか起こらず、$m=0$ の集合と合わせて $\mu_X$ に関して零である。積の長方形 $A\times B$ を考えると、$m(x)k(x,B)$ を $A$ 上で積分した値は同時密度を $A\times B$ 上で積分した値である。よって定義式を満たす。
:::

:::example[幅が変わる三角形で閾値を評価する]
同時密度 $2\mathbf1_{\{0<y<x<1\}}$ に対して、$\mathbb P(Y>1/2\mid X=x)$ は $0<x\le1/2$ で $0$、$1/2<x<1$ で $1-1/(2x)$ である。例えば $x=3/4$ では $1/3$ となる。この核を $X$ の密度 $2x$ で平均すると $\int_{1/2}^1(2x-1)\,\dd x=1/4$ であり、元の三角形から $Y>1/2$ の範囲を切り取った確率と一致する。
:::

:::example[零確率の点における核の違い]
$X=Y$ が $[0,1]$ 上の一様変数である。$k(x,\cdot)=\delta_x$ は正則条件付き分布である。一方 $x=1/2$ でだけ $\delta_0$ に変更しても、$\mu_X(\{1/2\})=0$ なので全ての分解積分は同じである。この場合、核の $x$ に関する弱連続性を要求すれば最初の版が選ばれる。弱連続な二つの版は $\mu_X$ の支え上で一致する。実際、その支えで一致する全測度集合は稠密であり、そこでの一致を連続性で延長できる。
:::

### 核で両方の変数に依存する量を積分する

:::proposition[正則条件付き分布による関数の積分]{#prop-prob-conditional-kernel-integration}
$(X,Y)$ の分解が $\mu_X(\dd x)k(x,\dd y)$ であるとする。$f$ が積空間上の非負可測関数なら $g(x)=\int f(x,y)k(x,\dd y)$ は可測で、
$$
\mathbb E[f(X,Y)\mid X]=g(X)
$$
が非負変数の拡張された意味で成り立つ。$f(X,Y)$ が可積分なら、正負部分の差を取り、例外集合で有限値を補って同じ式を使える。
:::

:::proof
まず $f=\mathbf1_{A\times B}$ では $g(x)=\mathbf1_A(x)k(x,B)$ であり、定義と持ち出しの性質から成立する。単調な集合列と補集合で可測性と各事象上の積分等式を保つので、公式の成立する集合の指示関数は $\lambda$ 系をなす。ここで補集合には、全体の積分が有限であることを使う。$\pi$–$\lambda$ 定理で長方形から積 Borel 集合へ拡張し、線形性と単調収束で非負 $f$ に拡張する。最後に可積分性は $g$ の正負部分に対応する二つの積分がともに $\mu_X$-概至るところ有限であることを保証する。
:::

:::example[観測に比例する閾値]
$0<y<x<1$ に同時密度 $2$ を置く。条件付きで $Y$ は $(0,x)$ 上一様なので $Y>X/3$ の条件付き確率は $2/3$ である。また $X=x$ において $\mathbb E[(Y-x/3)_+\mid X=x]=x^{-1}\int_{x/3}^x(y-x/3)\,\dd y=2x/9$ である。ここでは積分に使う関数自体が観測値 $x$ に依存する。最後に $X$ の密度が $2x$ であることを使えば、その無条件平均は $4/27$ となる。
:::

### 核の合成と周辺化

:::proposition[二段階の核の合成]{#prop-prob-conditional-kernel-composition}
実数値 $X,Y,Z$ に対し、$Y$ の $X=x$ における条件付き法則を $k(x,\dd y)$、$Z$ の $(X,Y)=(x,y)$ における条件付き法則を $l(x,y,\dd z)$ とする。すると
$$
m(x,B)=\int l(x,y,B)k(x,\dd y)
$$
は $Z$ の $X=x$ における条件付き法則である。
:::

:::proof
任意の $x$ で $m(x,\R)=1$ であり、可算加法性は $l(x,y,\cdot)$ の可算加法性と非負積分の順序交換による。$m$ の $x$ に関する可測性は、$k$ による二変数関数の積分の可測性から従う。次に任意の $A$ で $\{X\in A\}$ 上の $Z\in B$ の確率を先に $(X,Y)$ で条件付ける。条件付き確率は $l(X,Y,B)$ なので、さらに $X$ で条件付けると積分公式により $m(X,B)$ になる。塔の公式でこれが $Z$ の $X$ 条件付き確率である。
:::

:::example[観測から隠れた型を更新して次を予測する]
$X$ は公平な二値とする。$X=0$ のとき隠れたビット $Y=1$ の確率は $1/4$、$X=1$ では $3/4$ である。次の出力 $Z$ は $Y=0$ なら確率 $1/5$、$Y=1$ なら確率 $4/5$ で $1$ になる。$Z$ の条件付き法則は $Y$ だけに依存する。合成による成功確率は、$X=0$ で $(3/4)(1/5)+(1/4)(4/5)=7/20$、$X=1$ で $(1/4)(1/5)+(3/4)(4/5)=13/20$ である。$Z$ を直接 $X$ に関して計算したいとき、隠れた型を一点で代用せず、各型の予測を型の確率で平均する。
:::

:::example[離散成分と連続成分を同じ核で運ぶ]
$U$ は $(0,1)$ に一様、$D$ は独立な確率 $1/3$ の成功指示変数とし、$Y=DX$、$X=U$ とする。$x\in(0,1)$ での核は $(2/3)\delta_0+(1/3)\delta_x$ である。どの $x$ でも $Y$ の条件付き法則には Lebesgue 密度がない。ところが全体で $X$ を平均すると、零に質量 $2/3$ があり、残りは $(0,1)$ の密度 $1/3$ である。例えば $\mathbb P(1/2<Y\le1)=\int_{1/2}^1(1/3)\dd x=1/6$ である。正則条件付き分布を密度の商だけで作れない場合も、確率核として分解すれば同じ計算ができる。
:::

## 条件付き独立性

:::definition[条件付き独立性]{#def-prob-conditional-independence}
部分 $\sigma$-加法族 $\mathcal A,\mathcal B$ が $\mathcal G$ の下で :term[条件付き独立]{en="conditionally independent" yomi="じょうけんつきどくりつ"} であるとは、全 $A\in\mathcal A,B\in\mathcal B$ に対して $\mathbb P(A\cap B\mid\mathcal G)=\mathbb P(A\mid\mathcal G)\mathbb P(B\mid\mathcal G)$ が概確実に成り立つことをいう。
:::

指示関数から単関数と有界可測関数へ近似すると、$\mathcal A$ 可測な有界 $U$ と $\mathcal B$ 可測な有界 $V$ に対して $\mathbb E[UV\mid\mathcal G]=\mathbb E[U\mid\mathcal G]\mathbb E[V\mid\mathcal G]$ となる。逆は指示関数を代入すればよい。さらに条件付き独立性は、すべての有界 $\mathcal A$ 可測 $U$ について $\mathbb E[U\mid\mathcal G\vee\mathcal B]=\mathbb E[U\mid\mathcal G]$ となることと同値である。順方向では積分式を $G\cap B$（$G\in\mathcal G,B\in\mathcal B$）上で確かめ、$\pi$–$\lambda$ 定理で拡張する。逆方向では $B$ の指示関数を掛け、塔の公式を使う。

:::example[共通パラメータ]
成功率 $P$ が $1/4,3/4$ を等確率で取り、$P$ を固定すると $X,Y$ が独立な $\operatorname{Bern}(P)$ であるとする。無条件では $\mathbb EX=\mathbb EY=1/2$、$\mathbb EXY=\mathbb E P^2=5/16$ なので、共分散は $1/16$ で独立ではない。
:::

:::example[和を観測した後の依存]
公平な独立硬貨 $X,Y$ を零と一で表す。$X+Y=1$ を観測したときは各周辺の成功確率が $1/2$ だが、同時成功確率が零である。無条件に独立でも、新たな観測が二変数を結び付けることがある。
:::

### 条件付き独立性を外した後の相関

:::example[三つの観測で同時に生じる依存]
$P=0,1/2,1$ にそれぞれ確率 $1/4,1/2,1/4$ を置き、$P$ を固定して三つの Bernoulli 観測を独立に行う。観測の周辺成功率は $\mathbb EP=1/2$ である。最初の二つが成功する確率は $\mathbb EP^2=3/8$、三つとも成功する確率は $\mathbb EP^3=5/16$ である。従って二回の成功を見た後の第三の成功率は $(5/16)/(3/8)=5/6$ である。共通の $P$ を積分すると、成功したという情報によって第三の予測が上がる。
:::

:::example[条件付けで独立性が保たれる組]
公平な独立ビット $U,V,W$ に対し $T=U+V$ とする。$T$ は $(U,V)$ の関数であり、$W$ は $(U,V)$ と独立である。$T=t$ の確率が正なら、任意の $u,w$ に対し
$$
\mathbb P(U=u,W=w\mid T=t)
=\frac{\mathbb P(U=u,T=t)\mathbb P(W=w)}{\mathbb P(T=t)}.
$$
従って $U,W$ は $T$ の下で条件付き独立である。一方 $T=1$ の下では $U,V$ は $\mathbb P(U=1,V=1\mid T=1)=0$ なので独立ではない。どの観測を使って条件を作ったかにより、依存が誘導される組が変わる。
:::

### 指定した情報から予測を一段減らせる場合

:::proposition[独立な情報を加えた条件付き期待値]
$\mathcal A,\mathcal B$ が独立、$\mathcal G\subset\mathcal A$、$X$ が $\mathcal A$-可測で可積分とする。このとき
$$
\mathbb E[X\mid\mathcal G\vee\mathcal B]=\mathbb E[X\mid\mathcal G].
$$
:::

:::proof
$Z=\mathbb E[X\mid\mathcal G]$ は $\mathcal G\vee\mathcal B$ に関して可測である。$A\in\mathcal G,B\in\mathcal B$ に対し $\mathbb E[X\mathbf1_A\mathbf1_B]=\mathbb E[X\mathbf1_A]\mathbb P(B)$ である。これは $\mathcal A$-可測変数 $X\mathbf1_A$ と $\mathbf1_B$ の独立性を使った。同じ計算は $Z$ に対しても成立し、$\mathcal G$ での積分条件により両辺は等しい。$A\cap B$ の族は $\pi$ 系なので、差の積分が零となる集合の $\lambda$ 系は $\mathcal G\vee\mathcal B$ 全体を含む。
:::

:::example[実験の乱数と別の乱数を観測する]
公平な独立サイコロ $X,Y$ で、$X$ の偶奇を $\mathcal G$、$Y$ の出目を $\mathcal B$ とする。$X$ が偶数かだけを知ると条件付き平均は $4$、奇数なら $3$ である。さらに $Y=6$ を知っても、上の命題により $X$ の予測は変わらない。これに対し $X+Y=6$ を観測すれば新しい情報は $X$ も使っており、例えば $X$ が偶数のときは $2,4$ の二値に等確率で絞られるので予測は $3$ になる。
:::

### 条件付きの対独立と相互独立を分ける

:::example[各ペアの独立性は三つの同時確率を決めない]
公平な独立ビット $U,V$ と $W=U+V\bmod2$ を取り、これらと独立な別のビット $A$ を観測する。$A$ の各セルで $(U,V,W)$ の法則は元のままなので、どのペアも $A$ の下で条件付き独立である。ところが全て零の条件付き確率は $1/4$ であり、各周辺の積 $1/8$ と異なる。条件付き独立性を三つ以上に拡張するときは、任意の部分族で積の公式を要求する必要がある。条件を取る操作は、もともとの三変数の関係を自動的には消さない。
:::

:::proposition[粗い情報に移しても条件付き独立性が保たれる十分条件]
$\mathcal H\subset\mathcal G$ とし、$\mathcal A,\mathcal B$ が $\mathcal G$ の下で条件付き独立とする。全 $A\in\mathcal A$ について $\mathbb P(A\mid\mathcal G)$ が $\mathcal H$-可測な版を持つなら、$\mathcal A,\mathcal B$ は $\mathcal H$ の下でも条件付き独立である。
:::

:::proof
固定した $A\in\mathcal A,B\in\mathcal B$ について、塔の公式と $\mathcal G$ 条件付き独立性により
$$
\mathbb P(A\cap B\mid\mathcal H)
=\mathbb E[\mathbb P(A\mid\mathcal G)\mathbb P(B\mid\mathcal G)\mid\mathcal H].
$$
第一因子は $\mathcal H$-可測で有界なので取り出せる。各因子に塔の公式を使うと、第一因子は $\mathbb P(A\mid\mathcal H)$、残る条件付き期待値は $\mathbb P(B\mid\mathcal H)$ になる。これが主張の積公式である。
:::

:::example[共有した条件のどこまでが必要か]
独立な公平ビット $H,J$ を取り、両方を固定したとき $U,V$ を独立に引く。$U$ の成功率は $H=0$ で $1/4$、$H=1$ で $3/4$ とし、$V$ の成功率は $J=0$ で $1/5$、$J=1$ で $4/5$ とする。$\sigma(H,J)$ を $\sigma(H)$ に粗くしたとき、$U$ の条件付き分布は $J$ に依らず、上の命題が使える。例えば $H=1$ では両方成功の確率は $(3/4)(1/2)=3/8$ である。この場合さらに無条件でも独立である。一方、$V$ の成功率も $H$ に合わせて $1/5,4/5$ とすると、$H$ の条件を消した後の同時成功率は $(1/2)((1/4)(1/5)+(3/4)(4/5))=13/40$ であり、周辺の積 $1/4$ とは異なる。共有される母数を消すか、片側だけの母数を消すかで結果が変わる。
:::

### 二つの予測をどちらの順に行っても同じになる条件

:::proposition[共有情報の下で独立な二つの条件付け]
$\mathcal H\subset\mathcal A\cap\mathcal B$ とする。$\mathcal A,\mathcal B$ が $\mathcal H$ の下で条件付き独立なら、任意の可積分な $X$ について
$$
\mathbb E[\mathbb E[X\mid\mathcal A]\mid\mathcal B]
=\mathbb E[X\mid\mathcal H]
=\mathbb E[\mathbb E[X\mid\mathcal B]\mid\mathcal A]
$$
が成り立つ。逆に第一の等式が全ての有界 $X$ について成立すれば、条件付き独立である。
:::

:::proof
$U$ が有界な $\mathcal A$-可測変数なら、条件付き独立性により $\mathbb E[U\mathbf1_B\mid\mathcal H]=\mathbb E[U\mid\mathcal H]\mathbb E[\mathbf1_B\mid\mathcal H]$（$B\in\mathcal B$）である。全期待値を取り、右辺の第一因子は $\mathcal H$-可測なので $\mathbb E[\mathbb E[U\mid\mathcal H]\mathbf1_B]$ に等しい。$\mathcal H\subset\mathcal B$ なので、積分条件により $\mathbb E[U\mid\mathcal B]=\mathbb E[U\mid\mathcal H]$ である。有界切断と $L^1$ 縮小性でこの式は可積分な $U$ にも成立する。$U=\mathbb E[X\mid\mathcal A]$ を代入し、塔の公式を使う。逆順も対称性による。

逆に $X=\mathbf1_A$ を第一の等式に入れると $\mathbb E[\mathbf1_A\mid\mathcal B]=\mathbb P(A\mid\mathcal H)$ である。$B\in\mathcal B$ の指示関数を掛け、$\mathcal H$ で条件を取る。塔の公式と持ち出しにより $\mathbb P(A\cap B\mid\mathcal H)=\mathbb P(A\mid\mathcal H)\mathbb P(B\mid\mathcal H)$ が得られる。
:::

:::example[共有しているビットの情報だけが残る]
$H,U,V$ を独立公平ビットとし、$\mathcal A=\sigma(H,U)$、$\mathcal B=\sigma(H,V)$、$\mathcal H=\sigma(H)$ とする。$X=H+2U+4V$ では、まず $\mathcal A$ で条件を取ると $H+2U+2$、そこから $\mathcal B$ で条件を取ると $H+3$ である。逆に $\mathcal B$ では $H+1+4V$、そこから $\mathcal A$ に移すと同じ $H+3$ である。共有ビット $H$ は条件付けを二回行っても消えず、一方にしか入らない $U,V$ はどちらの順序でも平均値に置き換わる。
:::

## 正規分布の条件付き分布の計算

:::theorem[正規ベクトルの条件付き法則]{#thm-prob-gaussian-conditioning}
同時正規ベクトル $(X,Y)$ の平均を $m_X,m_Y$、共分散行列のブロックを $\Sigma_{XX},\Sigma_{XY},\Sigma_{YX},\Sigma_{YY}$ とし、$A=\Sigma_{YX}\Sigma_{XX}^{+}$ とする。ここで対称半正定値行列の擬逆 $\Sigma_{XX}^{+}$ は正の固有値を逆数に、零固有値を零に置き換えた行列である。$Y$ の $X=x$ における条件付き法則は、$\mu_X$ に関してほとんど至るところで
$$
N(m_Y+A(x-m_X),\,\Sigma_{YY}-\Sigma_{YX}\Sigma_{XX}^{+}\Sigma_{XY})
$$
で与えられる。
:::

:::proof
$\Sigma_{XX}v=0$ なら $v^\top(X-m_X)=0$ 概確実であり、$\Sigma_{YX}v=0$ である。従って $A\Sigma_{XX}=\Sigma_{YX}$ である。残差 $R=Y-m_Y-A(X-m_X)$ の $X$ との共分散は零であり、$(X,R)$ は同時正規なので特性関数が積に分かれる。:ref[thm-prob-vector-clt] に続く多次元特性関数の一意性により $X,R$ は独立である。残差の共分散を展開すると
$\Sigma_{YY}-A\Sigma_{XY}-\Sigma_{YX}A^\top+A\Sigma_{XX}A^\top=\Sigma_{YY}-A\Sigma_{XY}$ となる。独立残差に $x$ での予測値を足した分布を核とすれば、同時分布の長方形の確率が一致する。
:::

例えば平均零、共分散 $\begin{pmatrix}4&3\\3&9\end{pmatrix}$ の二変数では、$Y\mid X=x$ は $N(3x/4,27/4)$ である。一方 $X=(Z,Z)$、$Y=2Z+W$、$Z,W$ が独立な標準正規なら $\Sigma_{XX}$ は特異である。擬逆を使うと $A=(1,1)$ となり、$Y\mid X=(z,z)$ の法則は $N(2z,1)$ である。分布が質量を置く対角線上で式を読む。

### 精度で表す正規観測

:::example[雑音を含む測定から元の値を予測する]
独立な $Z\sim N(0,4)$、$\eta\sim N(0,1)$ に対して $X=Z+\eta$ を観測する。$\operatorname{Var}(X)=5$、$\operatorname{Cov}(Z,X)=4$ であり、$Z\mid X=x$ は $N(4x/5,4/5)$ である。例えば $x=5$ なら条件付き平均は $4$ であり、観測値そのものを予測に使うより平均が中央に寄る。共分散公式を使わず、密度の指数 $-z^2/8-(x-z)^2/2$ を $z$ で平方完成しても、係数 $-5/8$ と中心 $4x/5$ が得られる。
:::

:::example[独立な二つの測定を合わせる]
同じ $Z$ を、独立な $N(0,1)$ の雑音 $\eta_1,\eta_2$ を加えた $X_1,X_2$ で測る。事前の精度（分散の逆数）は $1/4$ であり、各測定の精度は $1$ である。指数 $-z^2/8-(x_1-z)^2/2-(x_2-z)^2/2$ の $z$ に関する二次係数は $-9/8$、一次項は $(x_1+x_2)z$ である。従って条件付き分散は $4/9$、条件付き平均は $4(x_1+x_2)/9$ である。例えば観測が $4,5$ なら平均は $4$ となる。一つの測定の分散 $4/5$ と比較すると、独立な二つ目の測定で残差が減る。
:::

### 時刻の間を正規条件付き分布で予測する

:::example[二つの端点の間での Brown 運動]
平均零の正規過程 $B$ で独立増分、$\operatorname{Var}(B_t-B_s)=t-s$ とする。$0<s<t<u$ とする。増分 $A=B_t-B_s$ と $C=B_u-B_t$ は互いに独立で、組としても $B_s$ と独立である。分散をそれぞれ $a=t-s$、$c=u-t$ と置く。$B_s=x$、$B_u=z$ を観測することは、$B_s=x$ と $A+C=z-x$ を観測することに等しいので、:ref[thm-prob-gaussian-conditioning] により
$$
B_t\mid(B_s=x,B_u=z)
\sim N\left(x+\frac{a}{a+c}(z-x),\frac{ac}{a+c}\right).
$$
例えば $s=1,t=2,u=4$、$x=1,z=4$ では平均は $2$、分散は $2/3$ である。始点だけを知る場合の分散は $1$ なので、終点の情報で残るゆらぎが小さくなる。
:::

:::example[正規な残差と単なる無相関を区別する]
$X\sim N(0,1)$ に対し $Y=X^2-1$ とする。$\mathbb EY=0$、$\operatorname{Cov}(X,Y)=0$ だが、$\mathbb E[Y\mid X]=Y$ であり一般に零ではない。例えば $X=2$ における正則な条件付き法則は $\delta_3$ である。最小二乗で直線を当てはめると傾きと切片はともに零になるが、$Y$ は $X$ の非線形関数として完全に予測できる。従って共分散の公式で得る直線が一般の条件付き期待値になるには、同時正規性のような条件が必要である。
:::

### 逆共分散行列から条件付き法則を読む

:::proposition[同時正規分布の精度行列による条件付け]{#prop-prob-gaussian-precision-conditioning}
$(X,Y)$ が正定値の共分散をもつ同時正規ベクトルで、その逆行列を精度行列 $K=\begin{pmatrix}K_{XX}&K_{XY}\\K_{YX}&K_{YY}\end{pmatrix}$ とする。平均を $(m_X,m_Y)$ とすると、$X=x$ の下で $Y$ の条件付き分布は
$$
N\left(m_Y-K_{YY}^{-1}K_{YX}(x-m_X),\,K_{YY}^{-1}\right)
$$
である。
:::

:::proof
$K$ は正定値であり、$Y$ の座標部分に制限すると $K_{YY}$ も正定値である。中心化した $u=x-m_X,v=y-m_Y$ と書けば、$y$ に依存する同時密度の指数は $-\tfrac12 v^\top K_{YY}v-v^\top K_{YX}u$ である。$a=K_{YY}^{-1}K_{YX}u$ と置くと、この部分は $-\tfrac12(v+a)^\top K_{YY}(v+a)+\tfrac12a^\top K_{YY}a$ である。最後の項と $x$ だけの残りは規格化で消える。$v+a$ を中心とした正規密度の積分で規格化すれば所望の平均と共分散を得る。
:::

:::example[条件付き独立性を精度行列の零成分から知る]
互いに独立な標準正規変数 $Z,E_1,E_2$ で $X=Z+E_1$、$Y=Z+E_2$ とする。$Z$ を固定すれば $X,Y$ は独立で、条件付き共分散は零、条件付き平均は両方 $Z$ である。一方で無条件の共分散は $\operatorname{Var}(Z)=1$ である。同時密度の指数は
$$
-\tfrac12(z^2+(x-z)^2+(y-z)^2)
$$
であり、$xy$ の項はない。従って $(X,Y)$ のブロックの精度は単位行列で、$Z$ の条件付き法則で二成分は独立である。逆に $Z$ を隠すと、それを通じて二成分が同じ向きに変わるので無条件での独立性は失われる。
:::

:::example[二変数の一方を観測した後に次を観測する]
$Z\sim N(0,2)$、$E_1,E_2\sim N(0,1)$ は全て独立とし、$X_i=Z+E_i$ とする。最初の $X_1=x_1$ での事後精度は $1/2+1=3/2$、平均は $(2/3)x_1$ である。次に $X_2=x_2$ を使うと精度は $5/2$、平均は $(2/5)(x_1+x_2)$ になる。実際、指数の一次係数は、最初の事後分布の精度と平均の積 $x_1$ に、次の測定の係数 $x_2$ を加えた値である。例えば $x_1=3,x_2=2$ なら最終平均は $2$、分散は $2/5$ である。二観測を最初からまとめて条件付けた場合と同じ結果になる。順に処理しても、前の事後分布に残った分散を捨てなければ計算は一致する。
:::

## 演習問題

::::exercise{level=basic}
事前確率 $1/4$ の箱から赤が出る確率は $2/3$、残りの箱からは $1/3$ とする。赤を得たときの最初の箱の確率を求めよ。

:::solution
赤の全確率は $(1/4)(2/3)+(3/4)(1/3)=5/12$ である。Bayes の公式で $(1/6)/(5/12)=2/5$ となる。
:::
::::

::::exercise{level=basic}
公平サイコロ $X$ の $X>4$ という事象上の条件付き平均を求めよ。

:::solution
生き残る出目 $5,6$ は等確率であり、平均は $11/2$ である。
:::
::::

::::exercise{level=basic}
$X$ が $\mathcal G$ 可測で可積分なら $\mathbb E[X\mid\mathcal G]$ は何か。

:::solution
$X$ 自身は $\mathcal G$ 可測で、すべての $A\in\mathcal G$ について定義の両辺が同じである。一意性により $X$ である。
:::
::::

::::exercise{level=basic}
$0<y<x<1$ に同時密度 $2$ があるとき、$\mathbb E[Y\mid X]$ を求めよ。

:::solution
$X=x$ で $Y$ の密度は $(0,x)$ 上で $1/x$ である。よって平均は $x^{-1}\int_0^x y\,\dd y=x/2$ であり、条件付き平均は $X/2$ である。
:::
::::

::::exercise{level=advanced}
$X$ と $Y$ は独立な標準正規である。$X+Y$ を知ったときの $X$ の法則を求めよ。

:::solution
$S=X+Y$ の分散は $2$、$\operatorname{Cov}(S,X)=1$ である。よって $X\mid S=s$ は $N(s/2,1-1/2)=N(s/2,1/2)$ となる。
:::
::::

::::exercise{level=advanced}
$X\in L^2$、$\mathcal H\subset\mathcal G$ とする。$L^2$ 誤差について $\mathbb E(X-\mathbb E[X\mid\mathcal H])^2\ge \mathbb E(X-\mathbb E[X\mid\mathcal G])^2$ を示せ。

:::solution
$\mathcal H$ に関する予測値は $L^2(\mathcal G)$ にも入る。最小二乗分解で左辺と右辺の差が二つの予測値の差の二乗平均になるため非負である。
:::
::::

::::exercise{level=advanced}
$U,V$ は $W$ の下で条件付き独立な有界変数とする。$\operatorname{Cov}(U,V)=\operatorname{Cov}(\mathbb E[U\mid W],\mathbb E[V\mid W])$ を示せ。

:::solution
条件付き独立性で条件付き積の平均が条件付き平均の積になる。全期待値を取り、$\mathbb EU=\mathbb E\mathbb E[U\mid W]$ と $V$ の同様の式を差し引けばよい。
:::
::::

::::exercise{level=advanced}
$X$ が $[0,1]$ 上で一様、$Y=X$ の場合、$k(x,B)=\mathbf1_B(x)$ が正則条件付き分布を与えることを示せ。

:::solution
各 $x$ で $\delta_x$ であり、$x$ について Borel 可測である。任意の Borel $A$ 上で積分すると $\mu_X(A\cap B)$ で、$\mathbb P(X\in A,Y\in B)$ と一致する。
:::
::::

::::exercise{level=advanced}
$X$ は $(0,1)$ 上で一様とし、$X=x$ の下で $Y$ は $(0,x)$ 上で一様とする。$\mathbb E[Y/X]$、$\operatorname{Cov}(X,Y)$ を求め、$X$ と $Y/X$ が独立であることを示せ。

:::solution
$X=x>0$ を固定すると、$Y/X=Y/x$ の条件付き法則は $(0,1)$ 上の一様分布である。任意の Borel 集合 $B$ に対するこの条件付き確率は $x$ に依らず、$B\cap(0,1)$ の Lebesgue 測度である。従って任意の事象 $X\in A$ との交叉は周辺確率の積になり、独立性と $\mathbb E[Y/X]=1/2$ が従う。一方 $\mathbb E[Y\mid X]=X/2$ であるから、$\mathbb EY=(1/2)(1/2)=1/4$、$\mathbb E[XY]=\mathbb E[X\cdot X/2]=1/6$ である。よって共分散は $1/6-(1/2)(1/4)=1/24$ となる。元の $Y$ には $X$ が決める尺度が残るが、比にすると尺度が除かれる。
:::
::::

::::exercise{level=advanced}
$X,Y$ は独立で、$X$ は率 $2$、$Y$ は率 $3$ の指数分布とする。$S=X+Y$ を観測したときの $X$ の条件付き密度と、$S=1$ での $\mathbb P(X>1/2\mid S=1)$ を求めよ。二つの率が同じ場合とも比較せよ。

:::solution
変数 $(x,y)$ を $(x,s=x+y)$ へ移すヤコビアンの絶対値は $1$ である。台は $0<x<s$ であり、同時密度は $6e^{-2x-3(s-x)}=6e^{-3s}e^x$ である。$s>0$ の周辺密度は $\int_0^s6e^{-3s}e^x\dd x=6e^{-3s}(e^s-1)$ だから、条件付き密度は
$$
k(s,x)=\frac{e^x}{e^s-1}\mathbf1_{\{0<x<s\}}
$$
である。従って $S=1$ での尾は $(e-e^{1/2})/(e-1)$ であり、分子分母から共通因子を消すと $e^{1/2}/(1+e^{1/2})$ である。この式は $1/2$ より大きい。率の小さい $X$ の方が、固定された合計の中で長い側の時間を占めやすい。

二つの率が同じ $\lambda$ なら、同時密度の $x$ に依存する因子が消え、$0<x<s$ 上で一定の $\lambda^2e^{-\lambda s}$ である。周辺密度で割ると $1/s$ となり、$S=1$ での尾は $1/2$ である。観測の値 $s$ が正確に等しい事象の確率は零だが、ここで求めた条件付き密度は分解積分を満たす版である。
:::
::::

## 確認テスト

::::quiz
条件付き期待値の一意性の意味はどれか。

- [x] 概確実に一致する。
- [ ] 全点で必ず一致する。
- [ ] 平均に関してだけ一致する。

:::explanation
定義の積分等式を二つの候補の差が正または負となる集合に適用すると、それらの集合の確率は零である。零集合上の変更は積分に現れない。
:::
::::

::::quiz{answer="1" tol="1e-6"}
$Y=3X+Z$、$X,Z$ が独立な標準正規とする。条件付き分散 $\operatorname{Var}(Y\mid X)$ を求めよ。

:::explanation
予測値は $3X$、残差は $Z$ である。$Z$ は $X$ と独立なので条件付き分散は元の分散 $1$ である。
:::
::::

::::quiz{answer="4" tol="1e-6"}
公平サイコロの出目について、偶数だと分かった場合の条件付き平均を求めよ。

:::explanation
偶数の三点 $2,4,6$ の条件付き確率はそれぞれ $1/3$ である。平均は $(2+4+6)/3=4$ となる。
:::
::::

::::quiz
塔の公式をそのまま適用できるのはどれか。

- [x] 大きい情報で条件付けた値を、そこに含まれる小さい情報で条件付ける場合。
- [ ] 二つの情報の包含関係を問わない場合。
- [ ] 条件付き期待値の代わりに任意の可測変数を置き換える場合。

:::explanation
$\mathcal H\subset\mathcal G$ なら $\mathcal H$ の各事象は $\mathcal G$ にも入る。従って $\mathbb E[\mathbb E[X\mid\mathcal G]\mid\mathcal H]=\mathbb E[X\mid\mathcal H]$ である。
:::
::::

::::quiz{answer="0.6" tol="1e-6"}
$\mathbb P(D)=1/5$、$\mathbb P(A\mid D)=3/4$、$\mathbb P(A\mid D^c)=1/4$ とする。さらに $\mathbb P(B\mid D\cap A)=2/3$、$\mathbb P(B\mid D^c\cap A)=1/3$ である。$\mathbb P(D\mid A\cap B)$ を答えよ。

:::explanation
$D\cap A\cap B$ の確率は $(1/5)(3/4)(2/3)=1/10$ である。$D^c\cap A\cap B$ の確率は $(4/5)(1/4)(1/3)=1/15$ である。両者を足した $1/6$ で前者を割るので、条件付き確率は $3/5=0.6$ である。
:::
::::

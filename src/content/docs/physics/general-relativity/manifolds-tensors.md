---
title: "時空の幾何学"
description: "時空の座標とテンソルの変換則を導き、計量・引き戻し・体積形式を計算する。"
notation:
  - symbol: "T_pM"
    meaning: "点 p の接空間"
  - symbol: "T_p^*M"
    meaning: "点 p の余接空間"
  - symbol: "g_{ab}"
    meaning: "Lorentz 計量"
  - symbol: "F^*"
    meaning: "写像 F による引き戻し"
  - symbol: "\\epsilon_g"
    meaning: "向きを選んだ計量体積形式"
---

## この章の目標

- 座標近傍と座標変換の意味を説明する。
- 接ベクトル、余接ベクトル、テンソルの成分を別の座標へ移す。
- Lorentz 計量で因果型と観測者の測定量を求める。
- 引き戻し、外微分、体積形式を積分に用いる。

特に断らない限り時空は四次元で、計量の符号は $(-,+,+,+)$、$c=1$ とする。接空間の線形代数だけを取り出した例では次元と署名を明記する。

## 多様体と座標変換

曲線に沿った測定値と位置の目盛りを分けるため、まず点の集合に座標を入れる。[多様体](/math/manifolds/) で扱う微分構造を、時空では四次元に取る。

:::definition[時空の座標近傍]{#def-gr-chart}
Hausdorff かつ第二可算な空間 $M$ を、$\R^4$ の開集合への同相写像 $x:U\to x(U)$ で覆う。重なりの写像 $y\circ x^{-1}$ と $x\circ y^{-1}$ が滑らかなとき、互いに両立するこの座標系の族は $M$ に滑らかな構造を与える。
:::

Jacobian を $J^\alpha{}_\mu=\partial y^\alpha/\partial x^\mu$、逆変換を $K^\mu{}_\alpha=\partial x^\mu/\partial y^\alpha$ とする。合成が恒等写像なので連鎖律から $KJ=I$ である。したがって局所的な変換則には可逆な行列が現れる。逆に、Jacobian の行列式が零になる点では、その写像をそこで座標変換に使えない。

:::example[円周の角度]
$S^1$ を切れ目の異なる二つの弧で覆う。重なりの各連結成分で角度は $\theta'=\theta$ または $\theta'=\theta+2\pi$ であり、どちらも微分は $1$ である。一方、区間 $[0,2\pi]$ は開集合ではなく、端の二つの角度は同じ点を指す。
:::

:::example[球面の極]
単位球の表示 $(\theta,\phi)\mapsto(\sin\theta\cos\phi,\sin\theta\sin\phi,\cos\theta)$ では $\theta=0$ で $\partial_\phi$ が零になる。極を含む北半球の座標には $(X,Y)=(\sin\theta\cos\phi,\sin\theta\sin\phi)$ を使い、$Z=\sqrt{1-X^2-Y^2}$ と表せる。極の問題は球面上の点の欠落でなく、角度表示の退化である。
:::

### 境界のない球面と座標が切れる位置

単位球面 $X^2+Y^2+Z^2=1$ の北極を除き、写像 $u=X/(1-Z),v=Y/(1-Z)$ を取る。$w=u^2+v^2$ とすると $X=2u/(1+w),Y=2v/(1+w),Z=(w-1)/(1+w)$ である。実際 $X^2+Y^2+Z^2=(4w+(w-1)^2)/(1+w)^2=1$ となり、元の射影比を戻せる。南極を除く射影 $\bar u=X/(1+Z),\bar v=Y/(1+Z)$ と重なる部分では $(\bar u,\bar v)=(u,v)/w$ である。この式は $w>0$ で滑らかで、同じ反転を二度行うと元に戻る。北極はこの射影の一つ目の座標に入らないが、二つ目では $(0,0)$ という通常の内点である。

:::example[座標変換を微分して向きを読む]
反転 $F(x)=x/(x^\top x)$ の微分は $DF=(I-2nn^\top)/w$ であり、$n=x/\sqrt w$ である。半径方向には固有値 $-1/w$、それと直交する方向には $1/w$ となる。例えば $(u,v)=(2,0)$ では行列は $\diag(-1/4,1/4)$ で、行列式は $-1/16$ である。球面上で整合する向きを選ぶには、一方の座標の片方の軸を反転する。このとき重なりでの変換の行列式は正になる。
:::

### 座標表示から多様体を作る

滑らかな構造に必要なのは、各座標の定義域が開集合で、重なりの変換が滑らかな微分同相となることである。両立する座標をすべて加えた極大アトラスは、はじめのアトラスから一意に定まる。実際、新しい座標二つの重なりの任意の点で元の座標を一つ選べば、二つの新座標の間の変換は元の座標を経由した滑らかな変換の合成で書ける。逆変換にも同じ議論が使える。このため、多様体の構造を指定するたびに全ての座標を列挙する必要はない。

:::example[三次元の空間断面と時空の座標]
$M=\R\times S^3$ とし、$S^3$ を $\R^4$ の中の単位球面とする。点 $X=(X^1,\ldots,X^4)$ で法線方向の成分 $X^A\ne0$ を選び、残りの三成分を座標に取る。元の $X^A$ は符号を固定して $X^A=\pm\sqrt{1-\sum_{i\ne A}(X^i)^2}$ から復元できる。球面を覆うには $A=1,\ldots,4$ と二通りの符号を用いる。重なりでは座標の入れ替えとこの平方根を合成する。時空の座標にはこの三成分と $\R$ の時刻の組を使えばよい。同じ時刻を表す条件は貼り合わせで保存するが、座標によらず定まる時刻関数が存在するという性質は、この直積で選んだ構造に由来する。
:::

### 制約で与えた曲面と座標数

$n+k$ 個の変数に対して滑らかな制約 $C^A(x)=0$ が $k$ 個あり、制約面上で微分 $\dd C^1,\ldots,\dd C^k$ が線形独立であるとする。Jacobian のある $k\times k$ 小行列式が非零なので、陰関数定理により対応する $k$ 変数を残る $n$ 変数の関数として解ける。これが制約面の局所座標になる。接ベクトルは全ての制約の微分を零にする。逆に、この線形条件を満たすベクトルは自由変数に関して任意の成分を持ち、従属成分は陰関数の微分に一致するため、面内の曲線で実現できる。

:::example[時間的な双曲面の接空間]
$\R^{1,4}$ で $-(X^0)^2+\sum_{i=1}^4(X^i)^2=L^2$ を考える。$L>0$ なら制約の勾配はどの点でも零にならず、この面は四次元である。点 $X=(0,0,0,0,L)$ では $\dd C=2L\dd X^4$ なので接空間は $V^4=0$ である。周囲の計量をこの接空間に制限すると $\diag(-1,1,1,1)$ が残る。他の点でも法線 $X$ は空間的であり、その直交補空間の署名は同じである。この構成では周囲の座標が五つあっても、時空内で動ける方向は四つである。
:::

### 座標の境界に近づく曲線を吟味する

座標の定義域が切れることと、その曲線が多様体内で終わることは別々の問題である。例えば Minkowski 平面で $t>|x|$、$t>0$ に $\tau=\sqrt{t^2-x^2}$、$\chi=\operatorname{artanh}(x/t)$ を取る。逆写像は $t=\tau\cosh\chi$、$x=\tau\sinh\chi$ であり、Jacobian の行列式は $\tau>0$ である。この座標が記述するのは原点の未来光円錐の内部である。例えば $(t,x)=(1,1)$ へ内部から近づくと $\tau\to0$、$\chi\to+\infty$ になる。元の Cartesian 座標ではその点も曲線も正則である。逆に各 $\chi$ を固定して $\tau\to0$ とするとどの曲線も同じ原点に近づく。座標の端でこの二種類の極限を同じとみなさないよう、元の事象への写像を調べる。

:::example[同じ平坦領域で三つの目盛りをつなぐ]
上の変換を微分すると
$$
\begin{pmatrix}\dd t\\\dd x\end{pmatrix}
=\begin{pmatrix}\cosh\chi&\tau\sinh\chi\\
\sinh\chi&\tau\cosh\chi\end{pmatrix}
\begin{pmatrix}\dd\tau\\\dd\chi\end{pmatrix}.
$$
従って計量は $-\dd\tau^2+\tau^2\dd\chi^2$ である。さらに $u=t-x=\tau e^{-\chi}$、$v=t+x=\tau e^\chi$ を取ると $-\dd u\,\dd v$ になる。時空の体積密度は Cartesian 座標で $1$、$(\tau,\chi)$ で $\tau$、$(u,v)$ で $1/2$ である。例えば $1<\tau<2$、$0<\chi<\log2$ の時空面積は $(3/2)\log2$ であり、$u,v$ では $\int\dd u\,\dd v/2$ で同じ値を得る。この場合、$u,v$ の積と比により積分領域は $1<uv<4$、$1<v/u<4$ となる。行列式の違いを比べるときには、同じ数値範囲ではなく、同じ事象の範囲を選ぶ。
:::

## 接ベクトルと余接ベクトル

曲線の接線を「関数をどれだけ変えるか」で測れば、座標を替えても同じ対象を比較できる。[接空間](/math/manifolds/tangent-spaces/) の定義を時空に適用する。

:::definition[接ベクトルと余接ベクトル]{#def-gr-tangent-cotangent}
滑らかな曲線 $\gamma$ が $\gamma(0)=p$ を満たすとき、$v(f)=\frac{\dd}{\dd s}(f\circ\gamma)(0)$ を接ベクトルとする。同じ作用を与える曲線を同値とし、接ベクトル全体を $T_pM$ と書く。:term[余接空間]{en="cotangent space" yomi="よせつくうかん"} $T_p^*M$ は $T_pM$ の線形汎関数の空間である。
:::

連鎖律から $v=v^\mu\partial_\mu$、$v^\mu=v(x^\mu)$ である。逆に任意の成分は座標で直線を取ると実現できるので、$T_pM$ は四次元である。双対基底を $\dd x^\mu$ と書けば $\dd x^\mu(\partial_\nu)=\delta^\mu_\nu$ となる。成分の変換は

$$
v'^\alpha=J^\alpha{}_\mu v^\mu,\qquad
\omega'_\alpha=K^\mu{}_\alpha\omega_\mu .
\label{eq-gr-dual-change}
$$

前者は $v(y^\alpha)$、後者は $\omega(\partial'_\alpha)$ の計算で得られる。縮約 $\omega(v)=\omega_\mu v^\mu$ では $KJ=I$ により変換が打ち消される。Lorentz 変換に限定した場合の具体例は :ref[def-sr-fourvector]{full} にある。

:::example[極座標での作用]
平面の $(x,y)=(3,4)$ で $\dd r=(3\dd x+4\dd y)/5$ である。$v=4\partial_x-3\partial_y$ なら $\dd r(v)=0$ なので、$v$ は円に接する。$\dd\phi=(-y\dd x+x\dd y)/r^2$ から $\dd\phi(v)=-1$ となる。
:::

:::example[関数の方向微分]
$f=x^2y$、$p=(2,3)$ なら $\dd f=12\dd x+4\dd y$ である。$v=\partial_x-2\partial_y$ による微分は $4$ であり、実際 $(2+s)^2(3-2s)=12+4s+\order{s^2}$ である。
:::

### 曲線の同値類から導分へ

同一点を通る二曲線 $\gamma_1,\gamma_2$ は、ある座標で同じ速度を持つとき同値とする。座標変換で速度を Jacobian が写すので、この条件は座標に依存しない。ここで点での導分は、その点の近傍で定まる関数を、関数がより小さい近傍で一致すれば同一視してから評価する。同値類に $D_\gamma f=(f\circ\gamma)'(0)$ を対応させると、通常の連鎖律から座標速度だけで決まり、積の微分則を満たす。逆に点での導分 $D$ を取る。座標を $x$、点を $x_0$ とし、座標の小球内で
$$
f(x)-f(x_0)=\sum_i(x^i-x_0^i)\int_0^1\partial_i f(x_0+s(x-x_0))\dd s
$$
と書く。$D(1)=0$ と積の法則を使うと $Df=\sum_i D(x^i)\partial_i f(x_0)$ となる。座標で直線 $x_0+sv$ を取れば $v^i=D(x^i)$ を実現できる。従って曲線と導分による定義は同じ接ベクトルを与える。導関数が一致すればよいため、加速度の異なる曲線も同じ接ベクトルを表し得る。詳細は :ref[thm-mfd-tangent-equivalent-models]{full} を参照せよ。

:::example[余接基底による速度の検算]
$x=r\cos\theta,y=r\sin\theta$、$r>0$ とする。$(r,\theta)=(2,0)$ で $(\dot r,\dot\theta)=(3,4)$ の曲線は直交座標で $(\dot x,\dot y)=(3,8)$ となる。同じ点で $\dd r=\dd x$、$\dd\theta=\dd y/2$ であるから、この速度への作用は $3,4$ を回復する。角方向の単位ベクトルは $\partial_\theta/r$ であるため、角速度 $4$ と単位枠での速度成分 $8$ は別の数になる。
:::

### 接束の貼り合わせとベクトル場の局所表示

各点の接空間を合わせた集合を $TM=\bigsqcup_{p\in M}T_pM$ とする。座標 $x=(x^1,\ldots,x^n)$ の領域では $(p,V)\mapsto(x(p),V^1,\ldots,V^n)$ が $\mathbb R^{2n}$ の開集合への全単射である。重なりで $y=y(x)$ とすれば $(x,V)\mapsto(y(x),(\partial y/\partial x)V)$ は滑らかで、逆写像も $x=x(y)$ の Jacobian から得られる。従ってこれらの座標で接束を滑らかな多様体にできる。射影 $\pi:TM\to M$ は局所的に最初の $n$ 成分を取る写像である。滑らかなベクトル場 $V$ は $\pi\circ V=\id_M$ を満たす切断であり、各座標で成分が滑らかという条件と同値である。同様に $T^*M$ では逆 Jacobian の転置に替えれば同じ構成ができる。一点を移動するときのベクトルの比較はこの貼り合わせだけからは定まらず、次章の接続を用いる。

### 座標基底に限らない局所枠

各点の接空間の基底を滑らかに選ぶと局所枠になる。座標 $x^\mu$ に対し $e_A=E^\mu{}_A\partial_\mu$ と書けば、双対枠は $\theta^A=(E^{-1})^A{}_\mu\dd x^\mu$ である。実際、$\theta^A(e_B)=(E^{-1})^A{}_\mu E^\mu{}_B=\delta^A_B$ である。座標を替えなくても、局所枠だけを変えれば観測者が使う時間軸と空間軸を指定できる。局所枠を座標基底と見なせるためには、少なくともどの二本も交換しなければならない。座標基底なら滑らかな関数への二階微分が交換するからである。

:::example[極座標の単位枠は座標基底にならない]
極座標の平面で $e_1=\partial_r,e_2=r^{-1}\partial_\phi$ とする。滑らかな関数 $F$ に作用させると
$$
e_1(e_2F)-e_2(e_1F)=-\frac1{r^2}\partial_\phi F=-\frac1r e_2F
$$
である。従ってこの組を同時に座標偏微分にすることはできない。双対枠は $\theta^1=\dd r,\theta^2=r\dd\phi$ であり、$r>0$ の局所領域でも $\theta^2$ は座標関数の微分ではない。角方向に同じだけ動くときの長さは半径に比例するので、この単位枠を座標網にするには一方の目盛りの間隔を保てない。
:::

### 高さ関数と等位面の法線

点 $p$ で $\dd f_p\ne0$ なら、関数 $f$ の等位面は余次元一の滑らかな面である。その接空間は $\ker\dd f_p$ である。この段階で法線として自然に与えられるのは余接ベクトル $\dd f$ である。方向を持つ法線ベクトルやその規格化には計量を使う。例えば平面の $f=x+2y$ で接線 $V=2\partial_x-\partial_y$ は $\dd f(V)=0$ であるが、Euclid 計量と $\diag(1,4)$ では法線ベクトルがそれぞれ $(1,2)$ と $(1,1/2)$ になる。どちらの法線も、対応する計量で接線と直交する。

### 二種類の微分を座標で照合する

前の例で $\partial_\tau=\cosh\chi\,\partial_t+\sinh\chi\,\partial_x$、$\partial_\chi=\tau(\sinh\chi\,\partial_t+\cosh\chi\,\partial_x)$ である。例えば $u=t-x$ に作用させると $\partial_\tau u=e^{-\chi}$、$\partial_\chi u=-\tau e^{-\chi}$ となり、$u=\tau e^{-\chi}$ を直接微分した結果を得る。一方 $\dd\tau$ の慣性枠での成分は $(t/\tau,-x/\tau)=(\cosh\chi,-\sinh\chi)$ である。二つを縮約すれば $\dd\tau(\partial_\tau)=1$、$\dd\tau(\partial_\chi)=0$ となる。ベクトルと一形式の成分は異なる変換で得られるため、行列のどちらを使うかは評価の不変性で点検できる。

:::example[計量で上げる操作と座標変更の順序]
$(\tau,\chi)$ で $\dd\tau$ の添字を上げると $-\partial_\tau$ である。先に慣性枠へ変換すると $(\cosh\chi,-\sinh\chi)$ になり、Minkowski の逆計量を掛ければ $(-\cosh\chi,-\sinh\chi)$ を得る。これは $-\partial_\tau$ を変換した成分である。一般に、一形式の成分が $K^{\mathsf T}\omega$ に移ると逆計量は $Jg^{-1}J^{\mathsf T}$ に移るので、積は $J g^{-1}\omega$ になる。添字を上げる演算が座標に依存しないためには計量と一形式の両方を変換する必要がある。
:::

## テンソル場と抽象添字記法

### 多重線形形式と縮約

$(r,s)$ 型の :term[テンソル]{en="tensor" yomi="てんそる"} は $r$ 個の余接ベクトルと $s$ 個の接ベクトルを受け取り、各入力について線形な実数値写像である。座標基底と双対基底を入力すると、すべての成分を得る。入力を基底展開して多重線形性を使うと、それらの成分から全入力での値が復元される。接ベクトルは余接ベクトルへの評価として $(1,0)$ 型である。

成分 $T^{\mu_1\cdots\mu_r}{}_{\nu_1\cdots\nu_s}$ は上添字ごとに $J$、下添字ごとに $K$ を受ける。添字の一つずつを同じ添字にして足す :term[縮約]{en="contraction" yomi="しゅくやく"} では、対応する $J^\alpha{}_\mu K^\nu{}_\alpha$ が $\delta^\nu_\mu$ となる。そのため結果は $(r-1,s-1)$ 型に変換する。各点でテンソルを選び、任意の座標で成分が滑らかに変わるものをテンソル場という。

### 抽象添字

$T^a{}_b$ のラテン小文字が :term[抽象添字]{en="abstract index notation" yomi="ちゅうしょうそえじ"} を表すとき、$a,b$ は成分番号でなく入力の種類を示す。$T^a{}_a$ は縮約されたスカラー、$T^a{}_bv^b$ はベクトルである。自由添字は等式の両辺と各加数で同じ型・並びにする。対称化・反対称化の例は $A_{(ab)}=(A_{ab}+A_{ba})/2$、$A_{[ab]}=(A_{ab}-A_{ba})/2$ である。この二つを足せば $A$ を復元する。

:::example[跡の座標不変性]
$A=\diag(2,5)$ を $J=\begin{pmatrix}1&1\\0&1\end{pmatrix}$ で変換すると $JAJ^{-1}=\begin{pmatrix}2&3\\0&5\end{pmatrix}$ である。成分 $A^1{}_2$ は変わるが縮約 $A^a{}_a=7$ は変わらない。
:::

:::example[反対称成分]
$A_{ab}$ の成分行列が $\begin{pmatrix}1&4\\2&3\end{pmatrix}$ なら、対称部分は非対角成分が $3$、反対称部分は $A_{[12]}=1$、$A_{[21]}=-1$ である。任意の $v$ に対し $A_{[ab]}v^av^b=0$ であり、添字を交換すると値が自分の負になることからも分かる。
:::

### 積と対称化は座標を選ばず定まる

$r$ 個の反変スロットと $s$ 個の共変スロットを持つテンソルは、余接と接ベクトルの引数に多重線形である。二つのテンソルの積を作るには、まず引数を二組に分けてそれぞれ評価し、得たスカラーを掛ける。その積を基底で評価すれば成分は各テンソルの成分の積となる。例えば二つの一形式 $\alpha,\beta$ に対し $A_{ab}=\alpha_a\beta_b$ とすると
$$
A_{(ab)}=\tfrac12(\alpha_a\beta_b+\alpha_b\beta_a),\qquad
A_{[ab]}=\tfrac12(\alpha_a\beta_b-\alpha_b\beta_a).
$$
これらを足せば元の $A$ になる。行列を転置したときの符号が異なるので、対称かつ反対称なテンソルは零であり、分解は一意である。

:::example[反対称部が消える評価と残る評価]
$\alpha=2\dd x+\dd y,\beta=\dd x-3\dd y$ とする。$A$ の成分は $\begin{pmatrix}2&-6\\1&-3\end{pmatrix}$、反対称部は $\begin{pmatrix}0&-7/2\\7/2&0\end{pmatrix}$ となる。この反対称部へ同じ $v$ を二度代入すると、交差する二項が相殺する。一方、$(\partial_x,\partial_y)$ へ代入すると $-7/2$ である。二つの引数の順序を独立に選べるから、計量の二乗ノルムだけを測る場合と異なり反対称部を復元できる。
:::

### テンソル作用素の跡と添字の上げ下げ

$(1,1)$ 型テンソル $A^a{}_b$ を線形写像と見る。基底を $S$ で変えると成分行列は $S^{-1}AS$ になり、$\tr(BC)=\tr(CB)$ は成分の総和の順序を交換して従う。従って $\tr(S^{-1}AS)=\tr A$ である。一方、共変二階テンソル $B_{ab}$ は $S^{\mathsf T}BS$ と変わるので行列の対角和は一般に保存しない。計量を使えば $g^{ab}B_{ab}=\tr(g^{-1}B)$ が定まり、$g^{-1}B$ の相似変換により不変である。例えば $g=\diag(-4,9)$、$B=\diag(8,9)$ では通常の対角和は $17$、計量との縮約は $-2+1=-1$ である。

### 成分を使ったテンソルの判定

同じ物理量について各座標で成分の候補を与えた場合、変換則を確認すればテンソル場として貼り合わせられる。例えば $A$ が成分 $A^\mu$ を持つとき、$\partial_\nu A^\mu$ には
$$
\partial'_\beta A'^\alpha
=K^\nu{}_\beta J^\alpha{}_\mu\partial_\nu A^\mu+
K^\nu{}_\beta(\partial_\nu J^\alpha{}_\mu)A^\mu
$$
という変換則がある。後半に座標変換の二階微分があるため、この成分列は一般には $(1,1)$ 型の変換則を満たさない。一方、スカラー $f$ の一階微分の成分は $K^\nu{}_\beta\partial_\nu f$ であり、二階微分には同じ種類の余分な項が生じる。共変微分を導入する際には、この余分な項を補う必要がある。

:::example[一次微分の対称・反対称部分の違い]
$x>0$ で $u=x^2$、$v=y$ と取り、$f=x+y$ とする。$x,y$ では $f$ のすべての二階偏微分は零だが、$u,v$ では $\partial_u^2 f=-1/(4u^{3/2})$ となる。従って二階偏微分をそのまま二階共変テンソルとみなせない。他方、任意の一形式 $\alpha=\alpha_\mu\dd x^\mu$ について $\partial'_\beta\alpha'_\gamma-\partial'_\gamma\alpha'_\beta$ を展開すると、$\partial'_\beta K^\nu{}_\gamma-\partial'_\gamma K^\nu{}_\beta=0$ により座標変換の二階微分が消える。これが一形式から二形式を作れる理由である。
:::

### 入力の交換と添字の添え方

四階の共変テンソル $A_{abcd}$ では、例えば $A_{[ab][cd]}$ は $a,b$ と $c,d$ を別々に反対称化する。$A_{[abcd]}$ は全四添字を反対称化する。前者と後者の条件は異なる。例えば $\alpha=\dd x\wedge\dd y$ に対して $A_{abcd}=\alpha_{ab}\alpha_{cd}$ は二組それぞれに反対称で、二組を交換しても不変である。一方、$A_{xyxy}=1$ だから全反対称ではない。全反対称化すると座標 $x,y$ しか現れず、その全ての成分が消える。曲率で使う四階テンソルも、二つの添字対を分けて扱う。

### 対称形式を測定から復元する

対称二階テンソル $S$ を測る際に、各方向の二次評価 $Q(v)=S(v,v)$ だけが与えられているとする。双線形性から $S(v,w)=[Q(v+w)-Q(v-w)]/4$ であり、全ての成分が復元できる。測れる方向が時間的な方向に限られる場合にも、一意性は残る。二つの候補の差の二次多項式が時間的な方向の開集合で零なら、そこでの二階微分が零なのでその係数が全て消えるためである。一方、反対称形式では同じ方向を二回入れると必ず零となる。観測可能な入力の組が、復元できるテンソルの部分を指定する。

:::example[観測方向を傾けて混合成分を分離する]
二次元の Minkowski 平面で $S_{00}=5,S_{01}=S_{10}=2,S_{11}=1$ とする。単位時間的な $u_\pm=\gamma(1,\pm v)$ に対し $S(u_\pm,u_\pm)=\gamma^2(5\pm4v+v^2)$ である。$v=3/5$ なら $\gamma=5/4$、二つの値は $97/8$ と $37/8$ である。一般に差は $4\gamma^2vS_{01}$、和は $2\gamma^2(S_{00}+v^2S_{11})$ である。静止枠の測定で $S_{00}$ を既知にすれば、二方向の測定から $S_{01},S_{11}$ を得る。観測者の速度の符号を交換する操作は、混合成分に対して奇、対角成分に対して偶になる。
:::

## 計量と Lorentz 署名

### 非退化な二次形式

:::definition[Lorentz 計量]{#def-gr-lorentz-metric}
:term[Lorentz 計量]{en="Lorentzian metric" yomi="ろーれんつけいりょう"} $g_{ab}$ は滑らかな対称二階共変テンソル場で、各点で非退化かつ負の固有値が一つ、正の固有値が三つである。基底に依存しないこの符号の組を署名という。
:::

対称行列を直交対角化し、各基底を固有値の絶対値の平方根で割ると $\eta=\diag(-1,1,1,1)$ になる。負定値となる部分空間の最大次元は $1$ である。二次元以上の部分空間は正値の三次元部分空間と交わるからである。正定値部分空間についても最大次元は $3$ となる。同じ議論は任意の基底で成り立つので署名は基底に依存しない。

非退化性は $g_{ac}g^{cb}=\delta_a{}^b$ で定まる逆があることを意味する。$v_a=g_{ab}v^b$ と $\omega^a=g^{ab}\omega_b$ は互いに逆の操作である。$v\ne0$ に対し $g(v,v)$ が負、零、正の順に時間的、光的、空間的と呼ぶ。未来向き単位時間ベクトル $u$ を観測者とし、$v=Eu+w$、$E=-g(u,v)$ とおけば $w=v-Eu$ は $g(u,w)=0$ を満たす。$u^\perp$ は正定値である。非退化性と直和 $\R u\oplus u^\perp$ によって $u^\perp$ は非退化であり、負方向があれば $u$ と合わせて負定値二次元平面を作ってしまう。

:::example[光のエネルギー]
二次元で光の運動量を $p=(6,6)$、観測者を $u=(5/4,3/4)$ とすると、$g(u,u)=-1$、$-g(u,p)=3$ である。$w=p-3u=(9/4,15/4)$ に対し $g(w,w)=9$ となるので、観測者の空間運動量の大きさはエネルギー $3$ と一致する。
:::

:::example[時空の交差成分]
$-\dd t^2+(\dd x+2\dd t)^2$ の行列は $\begin{pmatrix}3&2\\2&1\end{pmatrix}$ で、行列式が $-1$ である。$\partial_t$ は空間的だが $\partial_t-2\partial_x$ は時間的である。したがって $g_{tt}>0$ から時空計量の署名が変わったとは結論できない。
:::

### 非対角な光的基底

Minkowski 平面に $u=t-x,v=t+x$ を取ると $\dd s^2=-\dd u\,\dd v$ である。積の省略記法を展開すると $g_{uv}=g_{vu}=-1/2$ である。行列式は $-1/4$ なので非退化であり、$\partial_u,\partial_v$ はともに光的であっても二つの間の内積は非零である。逆行列は $g^{uv}=g^{vu}=-2$ である。従って $\dd u$ の添字を上げると $-2\partial_v$ になる。$u=\mathrm{const}$ に沿う方向が $\partial_v$ であるため、零超曲面では法線が接線にもなる。

:::example[署名を保つ斜交座標]
$t=\bar t,\ x=\bar x+2\bar t$ と置くと $g_{\bar t\bar t}=3,g_{\bar t\bar x}=2,g_{\bar x\bar x}=1$ である。行列式は $-1$、固有値は $2\pm\sqrt5$ なので、成分として正の対角項が並んでも時間方向は存在する。実際 $V=\partial_{\bar t}-2\partial_{\bar x}$ を代入すると $3-8+4=-1$ である。この方向は $x$ を固定した元の時間軸に対応する。
:::

### 時間の向きと観測者の射影

時間的ベクトルは一点で二つの連結成分をなす。正規直交枠では $|v^0|>\sqrt{\sum_i(v^i)^2}$ であり、$v^0>0$ と $v^0<0$ に分かれるからである。時空で連続に一方を選ぶことを :term[時間向き付け]{en="time orientation" yomi="じかんむきづけ"} という。大域的な滑らかな時間的ベクトル場 $X$ があれば、$g(X,v)<0$ の成分を未来向きと指定できる。逆に時間向きが指定されていれば、局所の未来時間的ベクトル場を分割の一で足し合わせれば大域的な場になる。未来時間的ベクトルの正係数和も未来時間的であることは、空間成分に対する三角不等式から従う。

観測者 $u$ に直交する方向への射影は $h^a{}_b=\delta^a_b+u^au_b$ である。$h(u)=0$、$h(V)=V$ が $V\in u^\perp$ で成り立つため $h^2=h$ である。四元運動量 $p^a$ でエネルギーを $E=-u\cdot p$、空間運動量を $P=h(p)$ とおく。質量殻 $g(p,p)=-m^2$ から $E^2=m^2+g(P,P)$ を得る。未来向きの $p$ に対して $E>0$ であり、$m>0$ なら観測速度は $P/E$ である。

:::example[同じ光線を二つの速度から測る]
局所慣性枠で $p=\omega(1,\cos\vartheta,\sin\vartheta,0)$ とする。別の観測者を $u=\gamma(1,v,0,0)$ とおき、その空間の第一軸を $e_1=\gamma(v,1,0,0)$ と取る。測るエネルギーは $\gamma\omega(1-v\cos\vartheta)$、第一空間成分は $\gamma\omega(\cos\vartheta-v)$、第二空間成分は $\omega\sin\vartheta$ である。従って第一軸に対する余弦は $(\cos\vartheta-v)/(1-v\cos\vartheta)$ となる。例えば $v=3/5$ で元の光が第一軸と直交する場合、エネルギーは $5\omega/4$、方向の成分は $(-3/5,4/5)$ になる。
:::

### 光的超曲面の退化と空間断面

$\Sigma$ を $\dd f\ne0$ を満たす $f$ の等位面とし、$n^a=g^{ab}\partial_b f$ とする。$\Sigma$ の接空間は $n^\perp$ である。$n$ が時間的なら誘導計量は正定値である。$n$ が空間的なら正規直交枠の空間軸を $n$ に沿わせると、残りの軸には時間軸が含まれ、誘導計量の署名は $(-,+,+)$ となる。$n$ が光的なら $n\in n^\perp$ である。そのため誘導計量は退化し、退化方向はちょうど $n$ の張る直線である。実際、非退化な周囲の計量では $(n^\perp)^\perp=\R n$ である。例えば $\dd s^2=-\dd u\,\dd v+\dd y^2+\dd z^2$ の $u=\mathrm{const}$ では、誘導形式は $\dd y^2+\dd z^2$ であり、$\partial_v$ に沿って零になる。一方、$u,v$ をともに固定した二次元断面では通常の正定値な面積を測れる。

### 正規直交枠を局所に選ぶ方法

点 $p$ で時間的なベクトルを取り、その座標成分を近傍へ滑らかに延長する。連続性により十分近くではノルムが負だから、平方根で割って単位軸 $e_0$ を得る。次に点で独立な三本の空間方向を選び、近傍へ延長して $h=\id+e_0\otimes e_0^\flat$ で射影する。独立性は開条件であり、$e_0^\perp$ は正定値なので通常の Gram–Schmidt 法を使える。こうして局所の正規直交枠を得る。二つの正規直交枠を結ぶ行列 $S$ は $S^{\mathsf T}\eta S=\eta$ を満たす。従って各点での変換は Lorentz 変換であり、$S$ が点ごとに違っても両方の枠は同じ計量を表す。

:::example[交差成分から観測者の正規枠を読む]
$\dd s^2=-N^2\dd t^2+A^2(\dd x+B\dd t)^2+\dd y^2+\dd z^2$ とし、$N,A>0$ とする。共枠を $\theta^0=N\dd t$、$\theta^1=A(\dd x+B\dd t)$ と取ると、双対枠は $e_0=N^{-1}(\partial_t-B\partial_x)$、$e_1=A^{-1}\partial_x$ である。$e_0$ は $t$ 一定の空間に直交する。$\partial_t$ に沿う時計も $A|B|<N$ なら時間的で、その速度は $w=\partial_t/\sqrt{N^2-A^2B^2}$ である。$w$ は $e_0,e_1$ の枠で速さ $|AB/N|$ を持つ。例えば $N=2,A=3,B=2/5$ なら $w$ は $e_0$ に対して $3/5$ の速度を持ち、同じ座標時間での固有時の率は $8/5$ である。$e_0$ に沿う時計の率は $2$ であり、両者は同じ点で比較しても別の世界線を持つ。
:::

## 写像の引き戻しと押し出し

### 微分と引き戻し

滑らかな $F:M\to N$ に対し、点での :term[押し出し]{en="pushforward" yomi="おしだし"} を $(F_*v)(f)=v(f\circ F)$ と定める。成分は $\partial_\mu F^\alpha v^\mu$ である。共変テンソル $S$ の :term[引き戻し]{en="pullback" yomi="ひきもどし"} は

$$
(F^*S)_p(v_1,\dots,v_s)=S_{F(p)}(F_*v_1,\dots,F_*v_s).
\label{eq-gr-pullback}
$$

連鎖律を二回用いると $(H\circ F)_*=H_*F_*$、$(H\circ F)^*=F^*H^*$ となる。$F$ が微分同相なら $F_*$ と $(F^{-1})_*$ を使ってどの型のテンソルも移せる。一般の写像では異なる原像のベクトルが同一点に来るので、点での微分から標的全体のベクトル場が決まるとは限らない。

:::example[曲面の誘導計量]
$F(u,v)=(u,v,u+2v)$ で Euclid 計量を引き戻すと $\dd u^2+\dd v^2+(\dd u+2\dd v)^2$ となる。行列は $\begin{pmatrix}2&2\\2&5\end{pmatrix}$ である。単位座標長方形の面積は $\sqrt{6}$ であり、$F_*\partial_u$ と $F_*\partial_v$ の Gram 行列の行列式からも同じ値を得る。
:::

:::example[退化する誘導形式]
Minkowski 平面への曲線 $F(s)=(s,s)$ では $F_*\partial_s=(1,1)\ne0$ であるが $F^*g=0$ である。光的な方向には時空の非退化な計量を制限しても零方向が残る。
:::

### 微分を合成する順序

$F:M\to N,G:N\to P$ とする。$f$ を $P$ 上の関数とし、$X\in T_pM$ に作用させると $d(G\circ F)X$ も $dG(dF X)$ も $X(f\circ G\circ F)$ を返す。従って微分は写像と同じ順で合成する。一形式 $\omega$ を引き戻すときは、これを任意の $X$ に評価して $\omega(dG(dF X))$ と書けるので、$(G\circ F)^*\omega=F^*(G^*\omega)$ である。同様に各スロットに微分を入れると、共変テンソルでも同じ順序が成り立つ。

:::example[曲線への引き戻しと周囲での微分]
曲線 $F(t)=(t,t^2)$ と一形式 $\omega=x\dd y-y\dd x$ を取る。$\dd y$ の引き戻しは $2t\dd t$ なので $F^*\omega=t^2\dd t$、$t=0$ から $1$ までの積分は $1/3$ である。周囲の二次元では $\dd\omega=2\dd x\wedge\dd y$ であるが、曲線に引き戻すと $2\dd t\wedge2t\dd t=0$ となる。次元が下がるので二形式の情報が失われても、引き戻しと外微分は交換する。
:::

### 写像に沿うベクトルと写像に関係するベクトル場

$F:M\to N$ に対し $F_*X$ は各 $p$ で $T_{F(p)}N$ のベクトルを与える。このように基点を $M$ の点として保持するものは、$F$ に沿うベクトル場である。$N$ 上の $Y$ と $F_*X(p)=Y(F(p))$ が成り立つとき、$X$ と $Y$ は $F$ で関係するという。この条件は任意の $f$ について $X(f\circ F)=(Yf)\circ F$ と書ける。もし $X_i$ と $Y_i$ がそれぞれ $F$ で関係すれば、作用の差を取って $[X_1,X_2](f\circ F)=([Y_1,Y_2]f)\circ F$ を得る。従って括弧も $F$ で関係する。

:::example[射影で同一の点へ送られる速度]
$F(x,y)=x$ とし、$X=y\partial_x+\partial_y$ とする。像の点 $x$ に送られる速度は $y\partial_x$ であり、同じ $x$ でも $y$ によって異なる。このため $F_*X$ は $F$ に沿うベクトル場としては定まるが、$\R$ 上の場にはならない。$X$ を $Z=(1+x^2)\partial_x+y\partial_y$ に替えると $F_*Z=(1+x^2)\partial_x$ で、標的の場を得る。縦の成分 $y\partial_y$ が消えることと、残る横成分が原像で一定であることの両方を確認している。
:::

### 直交座標への埋め込みで交差成分を調べる

滑らかなグラフ $F(t,x,y)=(t,x,y,s(t,x,y))$ を Minkowski 時空に入れると、誘導形式は
$$
F^*g=-\dd t^2+\dd x^2+\dd y^2
+(s_t\dd t+s_x\dd x+s_y\dd y)^2
$$
である。固定した $t$ の二次元切片では空間部分が $I+(\partial_i s)(\partial_j s)$ になり、接平面のどの方向にも正の長さを持つ。時間方向を含む誘導形式にはこの性質は自動ではない。例えば $s=bt+ax$ と置くと $t,x$ のブロックは
$$
\begin{pmatrix}b^2-1&ab\\ab&1+a^2\end{pmatrix}
$$
であり、行列式は $b^2-1-a^2$ である。$b^2<1+a^2$ なら時間的超曲面、等しいなら光的超曲面、$b^2>1+a^2$ なら空間的超曲面になる。$a=b=1$ なら $\partial_t$ は光的だが $y$ を含む誘導計量は Lorentz 署名を持つ。光的な接線の存在と誘導計量の退化を区別できる。

### 向きの異なる写像と保存する量

$F$ が微分同相なら、標的の共変テンソルを $F^*$ で移すと同時に、反変テンソルは逆写像の押し出しで移せる。これは $\omega(V)$ がスカラーであることから検算できる。$F^*\omega$ を $(F^{-1})_*V$ に評価すると、二つの微分は逆行列となって $\omega(V)\circ F$ になる。一般の $F$ について引き戻しがあるのは共変スロットであり、$\omega(V)$ を同じ形式で移すためには $V$ が $F$ で関係するベクトル場を持つなどの条件が要る。

:::example[同じ球面に異なる次数で写す]
円柱の帯 $(\theta,\phi)$ から単位球面への写像を
$$
F(\theta,\phi)=(\sin\theta\cos n\phi,\ \sin\theta\sin n\phi,\ \cos\theta)
$$
とする。$0<\theta<\pi$、$\phi$ は周期 $2\pi$、$n$ は零でない整数である。三次元の内積を引き戻すと $\dd\theta^2+n^2\sin^2\theta\,\dd\phi^2$、外向きの球面面積形式の引き戻しは $n\sin\theta\,\dd\theta\wedge\dd\phi$ となる。有向面積の積分は $4\pi n$、引き戻し計量の面積密度の積分は $4\pi|n|$ である。$n=-2$ でも計量は正定値であり、向きを反転して各点を二回覆うという違いが積分の符号と倍率に現れる。像の球面を一度だけ積分した面積 $4\pi$ とも区別する。
:::

## 微分形式と体積要素

### 外積と外微分

交代的な $k$ 階共変テンソル場を :term[微分形式]{en="differential form" yomi="びぶんけいしき"} という。座標では増加順の添字 $I=(i_1,\dots,i_k)$ に対し $\omega=\sum_I f_I\dd x^{i_1}\wedge\cdots\wedge\dd x^{i_k}$ と展開する。$\dd x^{i_1}\wedge\cdots\wedge\dd x^{i_k}$ を入力ベクトルに作用させると、各 $\dd x^{i_j}$ の評価を並べた行列の行列式になる。外積を作るときは基底を連結し、並べ替えの符号を付け、重複する添字を含む項を捨てる。外微分は $\dd\omega=\sum_I\dd f_I\wedge \dd x^I$ とする。$\partial_j\partial_l f_I$ が対称で基底の外積が反対称なので $\dd^2\omega=0$ である。

関数について $F^*\dd f=\dd(f\circ F)$ は連鎖律である。基底形式に適用すると $F^*\dd y^\alpha=\dd F^\alpha$ となり、外積への作用は行列式の多重線形性によって各因子への作用の外積になる。よって $F^*\omega=\sum_I(f_I\circ F)\dd F^{i_1}\wedge\cdots$ である。これを微分すると、$\dd^2F^\alpha=0$ により係数だけの微分が残る。したがって $\dd(F^*\omega)=F^*(\dd\omega)$ であり、外微分の定義は座標によらない。

:::example[回転と外微分]
$\omega=x\dd y-y\dd x$ なら $\dd\omega=2\dd x\wedge\dd y$ である。極座標へ引き戻すと $\omega=r^2\dd\theta$ となるので、微分は $2r\dd r\wedge\dd\theta$ である。半径 $2$ の円周で $\int\omega=8\pi$、円板で $\int\dd\omega=8\pi$ である。
:::

### 体積形式と発散

選んだ向きと整合する座標で $D=\det(g_{\mu\nu})<0$ とおく。:term[体積形式]{en="volume form" yomi="たいせきけいしき"} は

$$
\epsilon_g=\sqrt{-D}\,\dd x^0\wedge\dd x^1\wedge\dd x^2\wedge\dd x^3 .
\label{eq-gr-volume}
$$

計量行列が $K^\top gK$ に移ると平方根は $\abs{\det K}$ 倍になる。座標の外積は $\det J$ 倍なので、同じ向きの座標同士では両者が打ち消し合う。向きを選ばなくても $\sqrt{-D}\,\dd^4x$ を絶対値付きの密度として積分できる。ベクトル $V$ の収縮 $\iota_V\epsilon_g$ を微分する。座標で $V^\mu$ に対応する項は $\mu$ 番目の微分を取り除いた基底に $(-1)^\mu\sqrt{-D}V^\mu$ を掛けたものである。微分のうち $\partial_\mu$ だけが重複しないので、$\dd\iota_V\epsilon_g=\partial_\mu(\sqrt{-D}V^\mu)\dd^4x$ を得る。

:::example[球座標での体積密度]
Minkowski 時空で $g_{\mu\nu}=\diag(-1,1,r^2,r^2\sin^2\theta)$ なら $\sqrt{-D}=r^2\sin\theta$ である。$0<t<2$、$1<r<2$ の殻の時空体積は $2\cdot(8-1)/3\cdot4\pi=56\pi/3$ である。極の座標特異点は体積の積分に寄与しない。
:::

### Hodge 双対と時空の向き

電磁場のように二形式で表した量を方向ごとの流束へ変えるには、体積形式と逆計量を同時に使う。$k$ 形式の内積を $\langle\alpha,\beta\rangle_g=\alpha_{a_1\ldots a_k}\beta^{a_1\ldots a_k}/k!$ とすると、:term[Hodge 双対]{en="Hodge dual" yomi="ほっじそうつい"} $*\beta$ は $\alpha\wedge *\beta=\langle\alpha,\beta\rangle_g\epsilon_g$ を任意の $\alpha$ に対して満たす $(4-k)$ 形式である。この等式で $*\beta$ は一意に決まる。実際、増加順の基底外積を $\alpha$ として一つずつ入れると、補集合の添字を持つ係数が一つずつ判定できる。

向き付いた正規直交枠の基底形式を $\theta^A$ とする。$*\theta^I$ は補集合 $I^c$ の基底に、並べ替えの符号と $\theta^I$ の内積の符号を掛けたものになる。二度の補集合の交換による符号は $(-1)^{k(4-k)}$、内積の符号の積は時刻成分の負号一つである。よって Lorentz 署名で $*^2=(-1)^{k(4-k)+1}$ である。特に二形式に対して $*^2=-1$ であり、正定値の四次元の場合と符号が異なる。

:::example[二形式の双対]
$\epsilon_g=\theta^0\wedge\theta^1\wedge\theta^2\wedge\theta^3$ とすると、$\theta^0\wedge\theta^1$ のノルム二乗は $-1$ である。従って $*(\theta^0\wedge\theta^1)=-\theta^2\wedge\theta^3$ である。一方、$\theta^2\wedge\theta^3$ はノルム二乗 $1$ で、$(2,3,0,1)$ の置換は偶なので $*(\theta^2\wedge\theta^3)=\theta^0\wedge\theta^1$ である。二度適用したときの負号を直接読める。
:::

:::example[二球面での流束]
平坦時空の球座標で $\theta^0=\dd t,\theta^1=\dd r,\theta^2=r\dd\theta,\theta^3=r\sin\theta\dd\phi$ とする。電場に対応する $F=-Q_0\dd t\wedge\dd r/r^2$ には $*F=Q_0\sin\theta\,\dd\theta\wedge\dd\phi$ が対応する。球で積分すると $4\pi Q_0$ となり、電場の $1/r^2$ が面積の $r^2$ で打ち消される。
:::

### 外微分と積分の符号

有向長方形 $[0,2]\times[0,1]$ に $\dd x\wedge\dd y$ の向きを入れ、$\omega=x^2\dd y+y\dd x$ を取る。$\dd\omega=(2x-1)\dd x\wedge\dd y$ である。面積積分は $\int_0^1\int_0^2(2x-1)\dd x\dd y=2$ となる。境界を反時計回りにたどると、下辺は零、右辺は $4$、上辺は $-2$、左辺は零で、同じ $2$ を得る。一般の長方形でも $\omega=P\dd x+Q\dd y$ なら、左右の差は $\partial_x Q$ を面上で積分したもの、上下の差は $-\partial_y P$ を積分したものになる。これが一変数の微積分の基本定理から二次元の Stokes の定理が生じる仕組みである。

:::example[球面の面積形式を座標内で積分する]
単位球面の計量 $\dd\theta^2+\sin^2\theta\dd\phi^2$ に、$0<\theta<\pi$ で $\dd\theta\wedge\dd\phi$ の向きを取る。体積形式は $\sin\theta\dd\theta\wedge\dd\phi$ であり、北極から $\theta=\pi/3$ までの帽子の面積は $2\pi(1-\cos(\pi/3))=\pi$ である。極が座標に含まれない場合も、先に $\epsilon<\theta<\pi/3$ を積分し $\epsilon\to0$ とすれば、除いた部分の面積は $2\pi(1-\cos\epsilon)$ とともに零へ行く。
:::

### 向き・密度・体積形式の区別

計量行列の行列式は座標の Jacobian の二乗を掛けて変換する。従って $\sqrt{|\det g|}$ は Jacobian の絶対値を掛けて変換する密度である。一方、$\dd x^1\wedge\cdots\wedge\dd x^n$ の変換には Jacobian の符号が残る。向きの整合する座標だけを選ぶと正の Jacobian によって係数とくさび積が打ち消し合い、大域的な体積形式になる。向き付け不可能な場合にも密度は定まり、非負な体積の積分に使える。

:::example[向きを反転するときに何が変わるか]
二次元で $x'=x,y'=-y$ とし、元の体積形式を $\omega=\dd x\wedge\dd y$ とする。新しい座標では $\omega=-\dd x'\wedge\dd y'$ である。計量が直交なら両座標で $\sqrt{\det g}=1$ なので、体積形式の符号をこの数からだけは決められない。密度 $|\dd x\,\dd y|$ による長方形の面積は両方で同じ値を得る。
:::

:::example[Lorentz 署名の Hodge 双対を基底すべてで求める]
Cartesian 座標 $(t,x,y,z)$ の Minkowski 計量で、向きを $\dd t\wedge\dd x\wedge\dd y\wedge\dd z$ とする。$k$ 形式の基底の内積は、$k$ 個の一形式の内積を行列に並べたときの行列式で定める。順序付き基底形式 $\omega_I$ を一つ取ると、異なる集合 $I,J$ ではこの行列が零の行を持つので直交する。内積 $\langle\omega_I,\omega_I\rangle$ は $\dd t$ を含むと $-1$、含まないと $+1$ である。

$\omega_I$ と全く添字が重ならない補集合の形式を $\omega_{I^c}$ とする。$\omega_I\wedge\omega_{I^c}$ を体積形式の順序へ並べ替えたとき、転倒数が $N$ なら符号は $(-1)^N$ である。従って $\star\omega_I=(-1)^N\langle\omega_I,\omega_I\rangle\omega_{I^c}$ となる。他の基底形式を左から外積した場合は添字が重なって零となり、内積で定めた式にも一致する。以下の表は $k=0$ から $k=4$ の順に並べてある。

| 基底形式 | 転倒数 | 基底形式の内積 | Hodge 双対 |
|---|---|---|---|
| $1$ | 0 | $1$ | $\dd t\wedge \dd x\wedge \dd y\wedge \dd z$ |

| 基底形式 | 転倒数 | 基底形式の内積 | Hodge 双対 |
|---|---|---|---|
| $\dd t$ | 0 | $-1$ | $-\dd x\wedge \dd y\wedge \dd z$ |
| $\dd x$ | 1 | $1$ | $-\dd t\wedge \dd y\wedge \dd z$ |
| $\dd y$ | 2 | $1$ | $\dd t\wedge \dd x\wedge \dd z$ |
| $\dd z$ | 3 | $1$ | $-\dd t\wedge \dd x\wedge \dd y$ |

| 基底形式 | 転倒数 | 基底形式の内積 | Hodge 双対 |
|---|---|---|---|
| $\dd t\wedge \dd x$ | 0 | $-1$ | $-\dd y\wedge \dd z$ |
| $\dd t\wedge \dd y$ | 1 | $-1$ | $\dd x\wedge \dd z$ |
| $\dd t\wedge \dd z$ | 2 | $-1$ | $-\dd x\wedge \dd y$ |
| $\dd x\wedge \dd y$ | 2 | $1$ | $\dd t\wedge \dd z$ |
| $\dd x\wedge \dd z$ | 3 | $1$ | $-\dd t\wedge \dd y$ |
| $\dd y\wedge \dd z$ | 4 | $1$ | $\dd t\wedge \dd x$ |

| 基底形式 | 転倒数 | 基底形式の内積 | Hodge 双対 |
|---|---|---|---|
| $\dd t\wedge \dd x\wedge \dd y$ | 0 | $-1$ | $-\dd z$ |
| $\dd t\wedge \dd x\wedge \dd z$ | 1 | $-1$ | $\dd y$ |
| $\dd t\wedge \dd y\wedge \dd z$ | 2 | $-1$ | $-\dd x$ |
| $\dd x\wedge \dd y\wedge \dd z$ | 3 | $1$ | $-\dd t$ |

| 基底形式 | 転倒数 | 基底形式の内積 | Hodge 双対 |
|---|---|---|---|
| $\dd t\wedge \dd x\wedge \dd y\wedge \dd z$ | 0 | $-1$ | $-1$ |

例えば $\star(\dd t\wedge\dd x)=-\dd y\wedge\dd z$、$\star(\dd y\wedge\dd z)=\dd t\wedge\dd x$ である。従って二形式に二度作用させると符号が反転する。一般にも補集合をもう一度取ると、外積を入れ替える符号は $(-1)^{k(4-k)}$、全ての計量成分の符号の積は $-1$ なので $\star^2=(-1)^{k(4-k)+1}$ となる。向きを反転すると $\star$ は一回の作用ごとに全体が反転するが、二度作用する符号は変わらない。

電磁場の二形式を $F=E\dd t\wedge\dd x+B\dd y\wedge\dd z$ とすれば $\star F=B\dd t\wedge\dd x-E\dd y\wedge\dd z$ である。例えば $E=2,B=3$ なら $\langle F,F\rangle=-4+9=5$、$F\wedge F=2EB\,\mathrm{vol}=12\,\mathrm{vol}$ である。$F\wedge\star F=5\,\mathrm{vol}$ となるため、電場成分に付く負号はこの積からも検算できる。
:::

### 外積の規格化と次数付きの Leibniz 則

交代テンソルを成分で扱うときは
$$
\alpha=\frac{1}{k!}\alpha_{\mu_1\ldots\mu_k}\dd x^{\mu_1}\wedge\cdots\wedge\dd x^{\mu_k}
$$
と全添字を足すか、増加順の添字だけを足すかを選ぶ。全添字を足す場合は添字の置換の符号が成分と基底の両方に現れて打ち消し合い、各増加順の項が $k!$ 回数えられる。このため上の係数が必要である。$p$ 形式 $\alpha$ と $q$ 形式 $\beta$ の外積の成分は $(p+q)!/(p!q!)$ を掛けた $\alpha_{[\mu_1\ldots\mu_p}\beta_{\mu_{p+1}\ldots\mu_{p+q}]}$ である。実際、入力のうち $\alpha$ に渡す $p$ 個を選ぶ方法は $\binom{p+q}{p}$ 通りあり、選んだ組内での並べ替えは交代性によってすでに処理されている。

座標の基底形式だけの外微分は零なので、$\alpha=f\dd x^I,\beta=g\dd x^J$ として係数の積 $fg$ を微分する。$\dd (fg)=g\dd f+f\dd g$ の第二項で $\dd g$ を $p$ 個の基底の後ろに動かすと $(-1)^p$ が出る。従って
$$
\dd(\alpha\wedge\beta)=\dd\alpha\wedge\beta+(-1)^p\alpha\wedge\dd\beta
$$
である。任意の形式はこのような単項の和だから一般の場合も従う。なお $\alpha\wedge\beta=(-1)^{pq}\beta\wedge\alpha$ は全ての $p$ 個の基底を全ての $q$ 個の基底の向こうへ動かすと分かる。

:::example[時間と空間で分けた二形式]
三次元の空間座標 $x,y,z$ と時刻 $t$ を使い、時間依存する空間一形式 $E$ と二形式 $B$ を取る。$F=B-\dd t\wedge E$ とし、空間だけの外微分を $\dd_{\!s}$ とする。次数付き Leibniz 則から $\dd F=\dd_{\!s}B+\dd t\wedge(\partial_t B+\dd_{\!s}E)$ である。例えば $E=y\dd x$、$B=t\dd x\wedge\dd y$ では $\dd_{\!s}E=-\dd x\wedge\dd y$、$\partial_t B=\dd x\wedge\dd y$ で、和が消える。実際 $A=ty\dd x$ を取れば $-\dd A= -y\dd t\wedge\dd x+t\dd x\wedge\dd y$ となり、同じ $F$ を得る。全微分が零になる条件は、時間変化と空間微分の間の符号を決める。
:::

### 閉形式の局所ポテンシャルと大域的な障害

星形の座標領域で $H_t(x)=tx$ とする。$k>0$ の形式 $\omega$ に
$$
(h\omega)_x(v_1,\ldots,v_{k-1})
=\int_0^1 t^{k-1}\omega_{tx}(x,v_1,\ldots,v_{k-1})\dd t
$$
を対応させる。写像 $H(t,x)=tx$ の引き戻しを $\dd t\wedge A_t+B_t$ と分けると、$A_t$ が上の積分の被積分形式であり、$B_t=H_t^*\omega$ である。外微分を取った $\dd t$ の係数は $\partial_t B_t-\dd_{\!x}A_t$ だから、引き戻しとの交換を使うと $\dd h\omega+h\dd\omega=H_1^*\omega-H_0^*\omega=\omega$ である。従って $\dd\omega=0$ なら $\omega=\dd h\omega$ となる。この局所的な結論を :term[Poincaré の補題]{en="Poincaré lemma" yomi="ぽあんかれのほだい"} という。領域が原点以外の点に関して星形なら、原点をその点に移して同じ構成ができる。

:::example[穴の周りで閉形式の積分を測る]
$\R^2\setminus\{0\}$ 上の $\omega=(-y\dd x+x\dd y)/(x^2+y^2)$ は、極座標のパッチでは $\dd\phi$ である。従って外微分はパッチごとに零であり、全域で $\dd\omega=0$ である。半径 $R>0$ の円周を反時計回りに一周すると引き戻しは $\dd\phi$ なので積分は $2\pi$ になる。もし全域で $\omega=\dd f$ なら、円周の積分は $f$ の終値と初値の差で零になるから矛盾する。この障害は局所座標を改めても残る。電磁場をポテンシャルで表す場合も、局所で書けることと同じポテンシャルを全域に選べることを区別する。
:::

### 超曲面を貫く流束を計算する

時空体積形式 $\epsilon_g$ にベクトル場 $J$ を挿入した三形式を $j=\iota_J\epsilon_g$ とする。$\Sigma$ が空間的超曲面で未来単位法線を $n$ とすれば、誘導する向きを $\epsilon_\Sigma=\iota_n\epsilon_g$ と選べる。このとき $J=\rho n+J_\parallel$ の接線成分 $J_\parallel$ を入れると、残る三つの入力とともに三次元の接空間に収まるので外積が零になる。従って $\Sigma$ への $j$ の引き戻しは $\rho\epsilon_\Sigma$、ここで $\rho=-g(J,n)$ である。発散が零なら $\dd j=0$ である。:ref[thm-mfd-stokes]{full} により、二枚の空間切片に挟まれた領域で切片の荷の差は側面から出た流束の負になる。

:::example[運動する境界を通る流れ]
二次元時空 $(t,x)$ で $\epsilon=\dd t\wedge\dd x$、$J=\rho(\partial_t+v\partial_x)$ を定数 $\rho,v$ で取る。$j=\rho(\dd x-v\dd t)$ であり、境界 $x=R+wt$ に引き戻すと $\rho(w-v)\dd t$ となる。領域の右端では Stokes の外向きの境界向きは $t$ の減る向きだから、時刻 $0$ から $T$ の間に右端から外へ出る流束は $\rho(v-w)T$ である。例えば $w=v$ なら流れと境界が共に動くため零となる。右端だけを速さ $w$ で動かし左端を固定すると、空間切片で数えた荷は $\rho wT$ 増え、左端からの流入と右端からの流出の差に一致する。
:::

### 座標成分で Hodge 双対を計算する

向きを選び、体積形式の成分を $\epsilon_{a_1\ldots a_n}$ とする。$k$ 形式 $\alpha$ について
$$
(*\alpha)_{b_1\ldots b_{n-k}}
=\frac{1}{k!}\alpha^{a_1\ldots a_k}
\epsilon_{a_1\ldots a_k b_1\ldots b_{n-k}}
$$
である。点の正規直交枠で増加順の基底形式 $\alpha$ を入れると、上添字にする際の署名と、補集合を残す置換符号が現れる。$k!$ は同じ組の順序を全て数えた分を取り除くので、基底について上で得た定義に一致する。両辺はテンソルだから任意の座標でも同じ式で計算できる。例えば非対角な計量では $\alpha$ を先に逆計量で上げる操作が必要であり、下添字のまま置換の符号だけを付けても一般には正しくならない。

:::example[零座標で二形式を双対にする]
$\dd s^2=-2\dd u\,\dd v+\dd y^2+\dd z^2$ で向きを $\dd u\wedge\dd v\wedge\dd y\wedge\dd z$ とする。逆計量は $g^{uv}=g^{vu}=-1$ である。$\alpha=\dd u\wedge\dd y$ の上添字の非零成分は $\alpha^{vy}=-1$ とその逆順である。従って $*\alpha=-\dd u\wedge\dd z$ になる。$\beta=\dd u\wedge\dd z$ では $\beta^{vz}=-1$、置換 $(v,z,u,y)$ は奇なので $*\beta=\dd u\wedge\dd y$ である。二度作用すれば $*^2\alpha=-\alpha$ となる。また $\alpha\wedge*\alpha=0$ であり、$\alpha$ 自身は非零でも Lorentz 計量による二乗ノルムは零である。そのことは対になる $\dd v\wedge\dd y$ との内積が $-1$ であることと両立する。
:::

### 二形式が測る二つのスカラー

正規直交な共枠で $F=-\theta^0\wedge E+B$ とし、$E=E_i\theta^i$、$B=B_1\theta^2\wedge\theta^3+B_2\theta^3\wedge\theta^1+B_3\theta^1\wedge\theta^2$ とする。内積を取ると $\langle F,F\rangle=|\boldsymbol B|^2-|\boldsymbol E|^2$ である。二形式同士の外積は可換なので、
$$
F\wedge F=-2(\boldsymbol E\cdot\boldsymbol B)\epsilon_g,\qquad
F\wedge*F=(|\boldsymbol B|^2-|\boldsymbol E|^2)\epsilon_g
$$
となる。$B\wedge B$ は空間が三次元のため零であり、時間を含む項の二乗も $\theta^0$ の重複で零である。例えば $\boldsymbol E=(A,0,0)$、$\boldsymbol B=(0,A,0)$ では $F\ne0$ でも二つのスカラーはともに零となる。個々の観測者が分けて測る $\boldsymbol E,\boldsymbol B$ と、全ての向き付き正規直交枠で一致するスカラーを区別する。

### 体積の計算と簡単な境界積分

四次元 Minkowski 時空で $t_0<t<t_1$、$r_1<r<r_2$ の球殻を取り、$J=J^t(t,r)\partial_t+J^r(t,r)\partial_r$ とする。$\epsilon_g=r^2\sin\theta\,\dd t\wedge\dd r\wedge\dd\theta\wedge\dd\phi$ だから
$$
\iota_J\epsilon_g=r^2\sin\theta\,(J^t\dd r-J^r\dd t)\wedge\dd\theta\wedge\dd\phi .
$$
例えば $J^t=C,\ J^r=D/r^2$ を定数 $C,D$ で取ると、$\partial_t(r^2J^t)+\partial_r(r^2J^r)=0$ である。各切片の荷は $4\pi C(r_2^3-r_1^3)/3$、外面を出る流束は $4\pi D(t_1-t_0)$、内側から入る流束も同じ大きさである。従って時刻によらず荷は一定になる。$r_1\to0$ にしても $D\ne0$ の場は原点で正則にならない。原点を含まない殻で発散が零であることだけから、原点を含む全ての球面で流束が零とは結論できない。

:::example[向きの変化を積分で検算する]
Lorentz 計量の正規直交枠を、時間軸を保ったまま空間の一軸だけ反転させる。Hodge 双対の定義に用いる向きもその反転に合わせれば、空間断面を貫く流束の形式は全体の符号を変える。一方、有向断面の積分も向きを同時に反転すれば元と等しい。例えば向きを選んだ二球面での流束を $4\pi Q_0$ と規格化した場合、体積形式と二球面の向きをともに反転すると二つの負号が打ち消す。計量の密度だけで測る球面積 $4\pi r^2$ はどちらの向きでも等しい。
:::

## 演習問題

::::exercise{level=basic}
$u=2x+y$、$v=x+y$ について $\partial_u$、$\partial_v$ を求めよ。

:::solution
逆変換は $x=u-v$、$y=2v-u$ であるから $\partial_u=\partial_x-\partial_y$、$\partial_v=-\partial_x+2\partial_y$ である。
:::
::::

::::exercise{level=basic}
同じ変換で $3\dd x+\dd y$ の成分を求めよ。

:::solution
$\dd x=\dd u-\dd v$、$\dd y=2\dd v-\dd u$ より、$3\dd x+\dd y=2\dd u-\dd v$ である。
:::
::::

::::exercise{level=basic}
$g=\diag(-4,9)$ で $V=(3,2)$ の因果型と $V_a$ を求めよ。

:::solution
$g(V,V)=-4\cdot9+9\cdot4=0$ で $V\ne0$ なので光的である。$V_a=(-12,18)$ となる。
:::
::::

::::exercise{level=basic}
$F(u,v)=(u^2,v,u+v)$ の Euclid 計量の引き戻しの行列式を求めよ。

:::solution
接ベクトルは $(2u,0,1)$、$(0,1,1)$ なので、Gram 行列は $\begin{pmatrix}4u^2+1&1\\1&2\end{pmatrix}$ である。行列式は $8u^2+1$ となり、すべての点で正である。
:::
::::

::::exercise{level=advanced}
$\omega=xy\dd z+z\dd x$ に対して $\dd^2\omega=0$ を直接確認せよ。

:::solution
$\dd\omega=y\dd x\wedge\dd z+x\dd y\wedge\dd z+\dd z\wedge\dd x$ である。これを微分して得る項は $\dd y\wedge\dd x\wedge\dd z+\dd x\wedge\dd y\wedge\dd z=0$ である。
:::
::::

::::exercise{level=advanced}
$F(u,v)=(u,v,uv)$ に対し $\dd(F^*(x\dd z))=F^*\dd(x\dd z)$ を計算せよ。

:::solution
左辺の引き戻しは $uv\dd u+u^2\dd v$ で、その微分は $(-u+2u)\dd u\wedge\dd v=u\dd u\wedge\dd v$ である。右辺は $\dd u\wedge(v\dd u+u\dd v)=u\dd u\wedge\dd v$ となる。
:::
::::

::::exercise{level=advanced}
$(1,1)$ 型テンソルの跡が基底の選択によらないことを証明せよ。

:::solution
成分が $A'^\alpha{}_\beta=J^\alpha{}_\mu A^\mu{}_\nu K^\nu{}_\beta$ に変わる。$\alpha=\beta$ として足すと $K^\nu{}_\alpha J^\alpha{}_\mu=\delta^\nu_\mu$ となるので $A'^\alpha{}_\alpha=A^\mu{}_\mu$ である。
:::
::::

::::exercise{level=advanced}
球座標の空間で $V=r^2\partial_r$ の体積形式に関する発散を求めよ。

:::solution
体積密度は $r^2\sin\theta$ なので、$\dd\iota_V\epsilon=\partial_r(r^4\sin\theta)\dd r\wedge\dd\theta\wedge\dd\phi$ となる。体積形式で割ると $4r$ である。
:::
::::

::::exercise{level=advanced}
$F(u,v,w)=(u\cos v,u\sin v,w+\lambda v)$ を $u>0$ の局所座標変換とする。Euclid 計量の引き戻し、体積密度、および $\dd Z$ の引き戻しの平方ノルムを求めよ。$\lambda$ は実定数である。

:::solution
$\dd X^2+\dd Y^2=\dd u^2+u^2\dd v^2$、$\dd Z=\dd w+\lambda\dd v$ である。従って計量の $v,w$ ブロックは
$$
\begin{pmatrix}u^2+\lambda^2&\lambda\\ \lambda&1\end{pmatrix}
$$
で、$u$ 方向の成分は $1$ である。ブロックの行列式は $u^2$ なので体積密度は $u$ となる。逆ブロックは $u^{-2}\begin{pmatrix}1&-\lambda\\-\lambda&u^2+\lambda^2\end{pmatrix}$ である。これを共変成分 $(\lambda,1)$ に掛けると $(0,1)$ なので、平方ノルムは $1$ である。これは $\dd Z$ が元の座標では単位であったことと一致する。$u>0$ に制限したのは行列式が零になる軸を除くためで、$\lambda$ が非零でもその軸上では Jacobian は可逆にならない。
:::
::::

## 確認テスト

::::quiz
共変添字が一つの成分に作用する座標変換を選べ。

- [ ] $J^\alpha{}_\mu$
- [x] $K^\mu{}_\alpha$
- [ ] $\det J$

:::explanation
$\omega'_\alpha=\omega(\partial'_\alpha)=K^\mu{}_\alpha\omega_\mu$ である。ベクトルの変換は $J$ なので評価 $\omega(v)$ は保たれる。
:::
::::

::::quiz{answer="6" tol="1e-9"}
二次元の計量行列 $\diag(-4,9)$ の体積密度を答えよ。

:::explanation
行列式は $-36$ であるから $\sqrt{-\det g}=6$ である。
:::
::::

::::quiz
$F$ が任意の滑らかな写像であるとき、定義できるものをすべて選べ。

- [x] 点での接ベクトルの押し出し
- [x] 共変テンソル場の引き戻し
- [ ] 任意のベクトル場から標的全体のベクトル場への押し出し

:::explanation
点では微分 $F_*$ を適用でき、その作用によって共変テンソルのすべての入力を標的へ送れる。一方で $F$ が単射でない場合、同じ像に異なる接ベクトルが送られ得る。標的に $F$ の像に属さない点がある場合も、その点の値は得られない。
:::
::::

::::quiz{answer="7" tol="1e-9"}
固有値が $2,5$ の二次元の線形写像の跡を求めよ。

:::explanation
固有値が異なるので固有ベクトルが基底をなす。この基底で行列は $\diag(2,5)$ である。跡は縮約で基底によらず、$2+5=7$ となる。
:::
::::

::::quiz{answer="-1" tol="1e-7"}
計量行列 $\begin{pmatrix}3&2\\2&1\end{pmatrix}$ の行列式を答えよ。

:::explanation
対角成分の積 $3\cdot1$ から非対角の積 $2\cdot2$ を引くので $-1$ である。同じ符号の固有値では行列式が正になるため、この対称行列には正と負が一つずつある。
:::
::::

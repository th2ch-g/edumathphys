---
title: Abel–Jacobi の定理と Weierstrass 点
description: 正則微分の周期を用いて主因子を特徴付け、点ごとの極位数の違いを調べる
notation:
  - symbol: '\tau'
    meaning: a 周期を単位行列に正規化したときの b 周期行列
  - symbol: 'J(X)'
    meaning: 正則微分の周期格子で割った Jacobi 多様体
  - symbol: 'A(D)'
    meaning: 次数零の因子 D の Abel–Jacobi 像
---

## この章の目標

- 交点数で正規化した曲線に沿って正則微分の周期を取る。
- Riemann の双線形関係式を切り開いた曲面の境界から導く。
- Abel の定理と Jacobi の逆問題を証明する。
- Weierstrass 点を正則微分の消失位数で検出する。
- 楕円曲線で Abel–Jacobi 写像と群法則を結び付ける。

以下、$X$ は種数 $g$ のコンパクト連結 Riemann 面とする。

## 周期と正規化された正則微分

曲面の標準的な切断に対応する $1$ 次ホモロジーの基底を $a_1,b_1,\ldots,a_g,b_g$ とし、$a_i\cdot b_j=\delta_{ij}$、$a_i\cdot a_j=b_i\cdot b_j=0$ とする。基点を揃えて切り開くと境界語が $\prod_i a_i b_i a_i^{-1} b_i^{-1}$ となる多角形が得られる。[曲面の基本群とホモロジー](/math/algebraic-topology/) を参照する。形式 $\alpha$ の $a,b$ 周期を $A_i(\alpha),B_i(\alpha)$ と書く。

:::definition[正規化された基底]{#def-cg-normalized-periods}
$H^0(X,K_X)$ の基底 $(\omega_j)$ が $\int_{a_i}\omega_j=\delta_{ij}$ を満たすとき :term[正規化された正則微分]{en="normalized holomorphic differential" yomi="せいきかされたせいそくびぶん"} と呼ぶ。$b$ 周期行列を $\tau_{ij}=\int_{b_i}\omega_j$ とする。
:::

この基底の存在には、全 $a$ 周期が零なら微分が零であることが必要である。以下では微分の面積積分を周期から計算してこれを示す。

:::example[種数一の正規化]
$X=\C/(2\Z+(1+3i)\Z)$ では $2,1+3i$ を向きの合う順に取る。$\omega=\dd z/2$ により $a$ 周期は $1$、$b$ 周期は $(1+3i)/2$ である。
:::

### 周期行列を実線形写像として読む

$a$ 周期を単位行列にそろえた基底に対して、周期の実行列は
$$
(m,n)\longmapsto m+\tau n,\qquad
\begin{pmatrix}\Re(m+\tau n)\\\Im(m+\tau n)\end{pmatrix}
=\begin{pmatrix}I&\Re\tau\\0&\Im\tau\end{pmatrix}\begin{pmatrix}m\\n\end{pmatrix}.
$$
$\Im\tau$ の正定値性によりこの実行列は可逆である。従って実係数の周期の表現は一意である。特に $z\in\C^g$ が格子に入るかどうかは、$n=(\Im\tau)^{-1}\Im z$ と $m=\Re z-(\Re\tau)n$ がともに整数ベクトルかどうかで判定できる。

:::example[斜めの一次元格子で判定する]
$\tau=1/2+2i$ とする。$z=2+4i$ なら $n=2$、$m=1$ なので周期である。$z=2+i$ なら $n=1/2$ なので周期でないが、$2z=4+2i$ については $n=1,m=7/2$ で、依然として周期でない。$4z=8+4i$ で $n=2,m=7$ となり、$z$ の像は位数 $4$ を持つ。
:::

:::example[周期の基底を取り替える]
種数 $1$ で $(a',b')=(b,-a)$ とする。もとの $a$ 周期が $1$、$b$ 周期が $\tau$ なら、$\eta'=\eta/\tau$ は新しい $a'$ 周期が $1$ である。新しい $b'$ 周期は $\tau'=-1/\tau$ となる。$\tau=1+2i$ なら $\tau'=(-1+2i)/5$、その虚部は $2/5$ である。面積の倍率は $1/|\tau|^2=1/5$ で、正定値性が保たれる。
:::

### 非対角の周期行列の座標

:::example[周期格子で差の位数を求める]
次の周期行列を持つトーラスで、二つの代表ベクトル $p,q$ の差が周期かどうかを判定する。Abel 積分の計算が済んだ後には、この判定だけが必要になる。

$$
\begin{aligned}
\tau&=\left[\begin{matrix}2 i & i\\i & 3 i\end{matrix}\right] \\
\Im\tau&=\left[\begin{matrix}2 & 1\\1 & 3\end{matrix}\right] \\
p&=\left[\begin{matrix}2 + i\\\frac{1}{2} + 3 i\end{matrix}\right],\qquad q=\left[\begin{matrix}1\\2 i\end{matrix}\right] \\
p-q&=\left[\begin{matrix}1 + i\\\frac{1}{2} + i\end{matrix}\right]
\end{aligned}
$$

$p-q=m+\tau n$ と書く。虚部の方程式 $Bn=\Im(p-q)$ を、正定値の $B$ の逆行列で解く。その結果を実部の方程式に入れて $m$ を得る。

$$
\begin{aligned}
\det B&=5 \\
n&=B^{-1}\Im(p-q)=\left[\begin{matrix}\frac{2}{5}\\\frac{1}{5}\end{matrix}\right] \\
m&=\Re(p-q)-(\Re\tau)n=\left[\begin{matrix}1\\\frac{1}{2}\end{matrix}\right]
\end{aligned}
$$

格子に入る条件は全成分が整数となることに等しい。今回の全成分の分母の最小公倍数は $10$ なので、トーラス上の差の位数は $10$ である。

このデータが次数の等しい二因子から得られた場合には、この倍数だけ取った差が初めて主因子になる。
:::

:::example[周期格子で差の位数を求める]
次の周期行列を持つトーラスで、二つの代表ベクトル $p,q$ の差が周期かどうかを判定する。Abel 積分の計算が済んだ後には、この判定だけが必要になる。

$$
\begin{aligned}
\tau&=\left[\begin{matrix}1 + 3 i & -1 + i\\-1 + i & 2 i\end{matrix}\right] \\
\Im\tau&=\left[\begin{matrix}3 & 1\\1 & 2\end{matrix}\right] \\
p&=\left[\begin{matrix}3 + 2 i\\4 + i\end{matrix}\right],\qquad q=\left[\begin{matrix}0\\-1\end{matrix}\right] \\
p-q&=\left[\begin{matrix}3 + 2 i\\5 + i\end{matrix}\right]
\end{aligned}
$$

$p-q=m+\tau n$ と書く。虚部の方程式 $Bn=\Im(p-q)$ を、正定値の $B$ の逆行列で解く。その結果を実部の方程式に入れて $m$ を得る。

$$
\begin{aligned}
\det B&=5 \\
n&=B^{-1}\Im(p-q)=\left[\begin{matrix}\frac{3}{5}\\\frac{1}{5}\end{matrix}\right] \\
m&=\Re(p-q)-(\Re\tau)n=\left[\begin{matrix}\frac{13}{5}\\\frac{28}{5}\end{matrix}\right]
\end{aligned}
$$

格子に入る条件は全成分が整数となることに等しい。今回の全成分の分母の最小公倍数は $5$ なので、トーラス上の差の位数は $5$ である。

このデータが次数の等しい二因子から得られた場合には、この倍数だけ取った差が初めて主因子になる。
:::

## Riemann の双線形関係式

:::theorem[双線形積分公式]{#thm-cg-bilinear}
閉複素 $1$ 形式 $\alpha,\beta$ に対し
$$
\int_X\alpha\wedge\beta
=\sum_{i=1}^g \bigl(A_i(\alpha)B_i(\beta)-B_i(\alpha)A_i(\beta)\bigr).
\label{eq-cg-bilinear}
$$
$\alpha$ が正則で $\beta$ が有理型微分のとき、極を避ける切断多角形上で $F(p)=\int_{p_0}^p\alpha$ と取ると、右辺は $2\pi i\sum_p\operatorname{Res}_p(F\beta)$ になる。周期もこの多角形の辺に沿って取る。
:::

:::proof
多角形上で $\dd F=\alpha$ として $\int\dd(F\beta)=\int_{\partial} F\beta$ と書く。各四辺を進む順に $a,b,a^{-1},b^{-1}$ とする。向きを $a,b$ の向きに揃えて二つの対辺を比較する。$a^{-1}$ 上の $F$ は対応する $a$ 上より $B(\alpha)$ だけ大きいので、$a$ に平行な二辺は $-B(\alpha)A(\beta)$ を与える。$b$ 上では対応する $b^{-1}$ 上より $A(\alpha)$ だけ大きく、$A(\alpha)B(\beta)$ を得る。閉形式の一周期の総和は四辺で零だから、前の四辺の初期値は次の寄与に影響しない。

有理型の場合は極の小円を取り除く。$\alpha\wedge\beta=0$ なので多角形から穴を除いた面で $\dd(F\beta)=0$ である。円の境界は負向きとなり、対辺の寄与と釣り合う。従って右辺は小円の正向き積分の和であり、留数で表される。
:::

:::corollary[周期行列の性質]{#cor-cg-period-matrix}
正規化基底は一意に存在し、$\tau^\top=\tau$、$\Im\tau$ は正定値である。
:::

:::proof
正則 $\alpha$ に対し $(i/2)\int_X\alpha\wedge\bar\alpha$ は局所的に $|\alpha/\dd z|^2\dd x\dd y$ の積分なので、$\alpha\ne0$ なら正である。全 $a$ 周期が零なら :eqref[eq-cg-bilinear] で積分が零になる。従って $a$ 周期写像は単射で、:ref[prop-cg-genus-equality] から両辺の次元は $g$ だから同型である。$\omega_j\wedge\omega_k=0$ を公式に入れると $\tau_{jk}=\tau_{kj}$ となる。$\alpha=\sum_j v_j\omega_j$ と共役を入れると $(i/2)\int\alpha\wedge\bar\alpha=\bar v^\top(\Im\tau)v$ となり正定値性が従う。
:::

:::example[平行四辺形の面積]
$X=\C/(\Z+(2+i)\Z)$ では $\omega=\dd z$ で $\tau=2+i$ である。公式の右辺は $\bar\tau-\tau=-2i$ であり、$i/2$ を掛けると基本平行四辺形の面積 $1$ になる。
:::

### 面積行列の座標変換

正規化前の微分を行 $(\eta_1,\ldots,\eta_g)$ に並べ、$a$ 周期の行列を $A$、$b$ 周期の行列を $B$ とする。行が積分路、列が微分を表すなら、新しい微分は右から $A^{-1}$ を掛けて得られ、新しい $b$ 周期行列は $BA^{-1}$ である。$A$ の可逆性は :ref[cor-cg-period-matrix] の証明に含まれる。この演算ではどの添字を周期に割り当てたかを固定する必要がある。

:::example[二種類の基底を戻す]
$a$ 周期が $A=\left(\begin{smallmatrix}2&1\\0&1\end{smallmatrix}\right)$、$b$ 周期が $B=\left(\begin{smallmatrix}2i&i\\0&2i\end{smallmatrix}\right)$ であるとする。$A^{-1}=\left(\begin{smallmatrix}1/2&-1/2\\0&1\end{smallmatrix}\right)$ を掛けると $\tau=\diag(i,2i)$ であり、虚部は正定値である。これは周期正規化の行列演算の例であり、任意の対称正定値の周期行列が連結曲線から得られるという主張ではない。
:::

### 周期と積分の反対称性

微分 $\alpha,\beta$ が正則なら、点ごとに $\alpha\wedge\beta=0$ である。従って双線形公式の右辺も零になる。正規化した基底 $\omega_i,\omega_j$ を代入すれば $\tau_{ji}-\tau_{ij}=0$ となる。一方 $\alpha=\sum_j c_j\omega_j$ に対して $\frac i2\alpha\wedge\bar\alpha$ の積分は $c^\top(\Im\tau)\bar c$ である。これが $c\ne0$ で正なのは、$\alpha$ が零でない開集合で正の面積密度を与えるためである。従って対称性と正値性は別々の積から得られる。

:::example[候補行列の失敗を別々に診断する]
$\tau_1=i\begin{pmatrix}3&1\\1&2\end{pmatrix}$ の虚部の首座小行列式は $3,5$ なので正定値である。$\tau_2=i\begin{pmatrix}1&2\\2&1\end{pmatrix}$ は対称だが $(1,-1)$ に対する二次形式が $-2$ となる。$\tau_3=\begin{pmatrix}i&1\\0&i\end{pmatrix}$ の虚部は単位行列だが対称でない。後二者はそれぞれ異なる双線形関係に反する。
:::

:::example[対数微分の周期が整数を検出する]
球面の $f(z)=(z-a)^2/(z-b)^3$ について $\mathrm df/f=2\,\mathrm dz/(z-a)-3\,\mathrm dz/(z-b)$ である。二つの小円での積分は $4\pi i,-6\pi i$ であり、残りの無限遠の小円の積分は $2\pi i$ である。無限遠では $f$ が一位の零を持つため、この符号は有限の二留数の和の逆になる。
:::

### 実規格化した第三種微分

留数が実数で総和零の第三種微分 $\eta$ を取る。極の小円での周期は純虚数なので、閉曲線の周期の実部を考えるとき穴を埋めてよい。正則微分 $\omega$ の実部のコホモロジー類は全ての実一次類を尽くす。従って $\eta$ の周期の実部を消す $\omega$ が一意に存在する。別の規格化として $a$ 周期そのものを零にする条件があるが、両方の規格化は一般に同時には成立しない。

:::example[二つの規格化の補正係数]
周期行列が種数一で $\tau=1+2i$ とし、$a$ 周期を零とした第三種微分の $b$ 周期が $B=3+4i$ であったとする。$a$ 周期が $1$ の $\omega$ に対して $\eta+c\omega$ の $a$ 周期は $c$ である。実規格化のため $c=it$ と置くと、$b$ 周期の実部は $3-2t$ だから $t=3/2$ である。補正後の二つの周期は $3i/2$ と $11i/2$ で、両方とも純虚数である。$a$ 周期が零であるという最初の条件は補正で外れる。
:::

### 周期の正値性を数値で読む

:::example[周期格子で差の位数を求める]
次の周期行列を持つトーラスで、二つの代表ベクトル $p,q$ の差が周期かどうかを判定する。Abel 積分の計算が済んだ後には、この判定だけが必要になる。

$$
\begin{aligned}
\tau&=\left[\begin{matrix}2 i & -1 + i\\-1 + i & 4 i\end{matrix}\right] \\
\Im\tau&=\left[\begin{matrix}2 & 1\\1 & 4\end{matrix}\right] \\
p&=\left[\begin{matrix}\frac{1}{2} + 2 i\\3 + 3 i\end{matrix}\right],\qquad q=\left[\begin{matrix}1\\0\end{matrix}\right] \\
p-q&=\left[\begin{matrix}- \frac{1}{2} + 2 i\\3 + 3 i\end{matrix}\right]
\end{aligned}
$$

$p-q=m+\tau n$ と書く。虚部の方程式 $Bn=\Im(p-q)$ を、正定値の $B$ の逆行列で解く。その結果を実部の方程式に入れて $m$ を得る。

$$
\begin{aligned}
\det B&=7 \\
n&=B^{-1}\Im(p-q)=\left[\begin{matrix}\frac{5}{7}\\\frac{4}{7}\end{matrix}\right] \\
m&=\Re(p-q)-(\Re\tau)n=\left[\begin{matrix}\frac{1}{14}\\\frac{26}{7}\end{matrix}\right]
\end{aligned}
$$

格子に入る条件は全成分が整数となることに等しい。今回の全成分の分母の最小公倍数は $14$ なので、トーラス上の差の位数は $14$ である。

このデータが次数の等しい二因子から得られた場合には、この倍数だけ取った差が初めて主因子になる。
:::

## Abel の定理

### 道に依存する積分を格子で割る

$\boldsymbol\omega=(\omega_1,\ldots,\omega_g)^\top$ とし、$\Lambda=\Z^g+\tau\Z^g$ とする。$\Im\tau$ が可逆なので、$(x,y)\mapsto x+\tau y$ は $\R^{2g}$ から $\C^g$ への実同型である。従って $\Lambda$ は格子である。次数零の因子 $D=\sum n_p[p]$ について
$$
A(D)=\left[\sum_p n_p\int_{p_0}^p\boldsymbol\omega\right]\in\C^g/\Lambda
$$
と定める。積分路を変えると格子ベクトルが加わる。基点を変えると、共通の積分に $\sum n_p=0$ を掛けた分しか変わらない。

:::lemma[指定した留数をもつ微分]{#lem-cg-third-kind}
相異なる点 $p_1,\ldots,p_r$ と複素数 $n_1,\ldots,n_r$ が $\sum n_j=0$ を満たすとする。これらの点に限って単純極を許し、留数 $n_j$、全 $a$ 周期零である有理型微分 $\eta$ が一意に存在する。ただし周期を測るループは点を避けて固定する。
:::

:::proof
$r=0$ では $\eta=0$ とする。$r>0$ では $P=\sum[p_j]$ として、留数写像 $H^0(K_X(P))\to\C^r$ の核は $H^0(K_X)$ である。:eqref[eq-cg-classical-rr] と $h^0(\mathcal O_X(-P))=0$ により $h^0(K_X(P))=g+r-1$ である。従って像は次元 $r-1$ をもち、留数定理の総和零の超平面と一致する。微分を得たら各 $a$ 周期を係数として正規化済み基底を引けばよい。差が正則で $a$ 周期零なら零なので一意である。
:::

:::theorem[Abel の定理]{#thm-cg-abel}
次数零の因子 $D$ が主因子であることと、$A(D)=0$ が同値である。
:::

:::proof
$A(D)$ は基点に依存しないので、$p_0$ を $D$ の台の外に選ぶ。$D$ の係数を留数として :ref[lem-cg-third-kind] の $\eta$ を取る。切断多角形内で積分路を選び、:ref[thm-cg-bilinear] を $\omega_j,\eta$ に適用すると $B_j(\eta)=2\pi i\sum_p n_p\int_{p_0}^p\omega_j$ である。従って $A(D)=0$ なら $B(\eta)/(2\pi i)=m+\tau n$、$m,n\in\Z^g$ と書ける。$\eta'= \eta-2\pi i\sum_j n_j\omega_j$ は $a,b$ 周期が $2\pi i$ の整数倍である。また極を囲む小円の周期も留数が整数なので同じ形になる。これらの道が穴あき曲面の $H_1$ を生成するため、$f(p)=\exp\int_{p_0}^p\eta'$ は穴の外で単価である。極の近くで $\eta'=n_p\dd z/z+\dd h$ と書けば $f=Cz^{n_p}e^h$ となり、有理型に延長して $(f)=D$ を得る。

逆に $D=(f)$ なら $\dd f/f$ は要求された留数をもち、全周期が $2\pi i$ の整数倍である。これは $f$ の像の巻き数である。$a$ 周期 $2\pi i n$ を正則微分で引いて正規化すると、$b$ 周期を $2\pi i$ で割った値は整数ベクトルと $\tau$ の整数倍ベクトルの和になる。先の公式から $A(D)=0$ が従う。
:::

:::example[球面とトーラスでの違い]
球面では $g=0$ なので全ての次数零因子が主因子である。平方格子トーラスでは $(\tfrac12)+( \tfrac i2)-(0)-(\tfrac{1+i}2)$ の Abel 像は $0$ なので主因子である。一方 $(\tfrac12)-(0)$ は $1/2$ が格子に入らないため主因子でない。
:::

### 主因子を作った後の一意性と平行移動

同じ因子 $D$ を持つ二関数 $f,g$ は $f/g$ の零も極もないため、コンパクト性から定数倍だけ異なる。従って Abel の定理は存在を判定すると同時に、基点での非零値を一つ指定すれば関数が一意であることを保証する。因子を次数零に限る理由もここで分かる。基点を $p_0$ から $q_0$ へ変えると、全ての積分に同じ定数ベクトルが加わるため、$D$ の係数和が零である場合に限りその変化は必ず消える。

:::example[超楕円射影から主因子を直接作る]
$y^2=P(x)$ の正規化で、$a,b$ を分岐値でない相異なる有限値とする。$x=a$ 上の点を $p,p'$、$x=b$ 上を $q,q'$ と書けば
$$
\left(\frac{x-a}{x-b}\right)=[p]+[p']-[q]-[q'].
$$
従って Abel 像の和 各点の像を $A(p)=A([p]-[p_0])$ と書けば $A(p)+A(p')$ はファイバーを変えても同じである。$b$ が単純分岐値なら $[q]+[q']$ を $2[q]$ と理解する。$x-b$ の位数が分岐点で二位となるためである。
:::

:::example[格子トーラスの有効因子の比較]
$\Lambda=\Z+(1/2+i)\Z$ 上で $p=1/4,q=3/4,r=i/2,s=-i/2$ と置く。$p+q-r-s=1$ は格子に入るので $[p]+[q]-[r]-[s]$ は主因子である。$q$ を $q+1/3$ へ動かすと差は $4/3$ になり、虚部零の格子点は整数だけなので主因子でなくなる。
:::

### 周期判定と捩れ

:::example[周期格子で差の位数を求める]
次の周期行列を持つトーラスで、二つの代表ベクトル $p,q$ の差が周期かどうかを判定する。Abel 積分の計算が済んだ後には、この判定だけが必要になる。

$$
\begin{aligned}
\tau&=\left[\begin{matrix}2 i\end{matrix}\right] \\
\Im\tau&=\left[\begin{matrix}2\end{matrix}\right] \\
p&=\left[\begin{matrix}\frac{1}{3} + i\end{matrix}\right],\qquad q=\left[\begin{matrix}\frac{1}{2}\end{matrix}\right] \\
p-q&=\left[\begin{matrix}- \frac{1}{6} + i\end{matrix}\right]
\end{aligned}
$$

$p-q=m+\tau n$ と書く。虚部の方程式 $Bn=\Im(p-q)$ を、正定値の $B$ の逆行列で解く。その結果を実部の方程式に入れて $m$ を得る。

$$
\begin{aligned}
\det B&=2 \\
n&=B^{-1}\Im(p-q)=\left[\begin{matrix}\frac{1}{2}\end{matrix}\right] \\
m&=\Re(p-q)-(\Re\tau)n=\left[\begin{matrix}- \frac{1}{6}\end{matrix}\right]
\end{aligned}
$$

格子に入る条件は全成分が整数となることに等しい。今回の全成分の分母の最小公倍数は $6$ なので、トーラス上の差の位数は $6$ である。

このデータが次数の等しい二因子から得られた場合には、この倍数だけ取った差が初めて主因子になる。
:::

:::example[周期格子で差の位数を求める]
次の周期行列を持つトーラスで、二つの代表ベクトル $p,q$ の差が周期かどうかを判定する。Abel 積分の計算が済んだ後には、この判定だけが必要になる。

$$
\begin{aligned}
\tau&=\left[\begin{matrix}\frac{1}{2} + 3 i\end{matrix}\right] \\
\Im\tau&=\left[\begin{matrix}3\end{matrix}\right] \\
p&=\left[\begin{matrix}3 + i\end{matrix}\right],\qquad q=\left[\begin{matrix}\frac{1}{2}\end{matrix}\right] \\
p-q&=\left[\begin{matrix}\frac{5}{2} + i\end{matrix}\right]
\end{aligned}
$$

$p-q=m+\tau n$ と書く。虚部の方程式 $Bn=\Im(p-q)$ を、正定値の $B$ の逆行列で解く。その結果を実部の方程式に入れて $m$ を得る。

$$
\begin{aligned}
\det B&=3 \\
n&=B^{-1}\Im(p-q)=\left[\begin{matrix}\frac{1}{3}\end{matrix}\right] \\
m&=\Re(p-q)-(\Re\tau)n=\left[\begin{matrix}\frac{7}{3}\end{matrix}\right]
\end{aligned}
$$

格子に入る条件は全成分が整数となることに等しい。今回の全成分の分母の最小公倍数は $3$ なので、トーラス上の差の位数は $3$ である。

このデータが次数の等しい二因子から得られた場合には、この倍数だけ取った差が初めて主因子になる。
:::

## Jacobi の逆問題

:::theorem[Jacobi の逆問題]{#thm-cg-jacobi-inversion}
任意の $u\in\C^g/\Lambda$ は $A([p_1]+\cdots+[p_g]-g[p_0])$ と書ける。点は重複してよい。
:::

:::proof
$g=0$ は標的が一点なので自明である。$\Phi:X^g\to\C^g/\Lambda$ を各点の Abel 積分の和とする。点 $(p_1,\ldots,p_g)$ における微分の列は、それぞれの座標での $\boldsymbol\omega(p_j)$ の係数ベクトルである。これらのベクトルが張る空間が $\C^g$ より小さければ、線形汎関数で消せるので基底微分の非自明な線形結合が $X$ の全点で消える。従って $g$ 本の独立な列を選べて、$\Phi$ はどこかで局所同型になる。

標的の正の体積形式 $\nu$ の引き戻しは正則写像の実 Jacobian が $|\det\dd\Phi|^2$ なので非負で、局所同型の点の近くで正である。よって $\int_{X^g}\Phi^*\nu>0$ である。像が標的全体でなければ、像は閉なので補集合に台を持ち $\int\nu$ と同じ積分をもつ滑らかな最高次形式 $\nu'$ を選べる。連結コンパクト向き付け可能多様体の最高次 de Rham コホモロジーは積分で同定されるため、$\nu-\nu'$ は完全である。:ref[thm-mfd-top-cohomology] を参照する。Stokes の定理で引き戻しの積分は等しいが $\Phi^*\nu'=0$ なので矛盾する。
:::

代数的には次数 $g$ の任意の直線束に Riemann–Roch で少なくとも一つの非零切断がある。二つ以上あれば、その切断の零因子の取り方は一意でない。切断の比の因子が零なら比が定数なので、有効代表の選択は非零切断を定数倍で割った射影空間によって記述される。

### 一般の点の逆像を数える

順序付き和 $\Phi:X^g\to J(X)$ の位相的次数は $g!$ である。計算するため、$\nu$ を $\C^g$ の標準体積形式からトーラスに降ろす。$\Lambda$ の基本領域を $(x,y)\mapsto x+\tau y$ で表したとき、実 Jacobian の行列はブロック三角で行列式は $\det\Im\tau$ になる。従って $\int_J\nu=\det\Im\tau$ である。

$A_{jk}=(i/2)\int_X\omega_j\wedge\bar\omega_k=(\Im\tau)_{jk}$ とする。$\Phi$ の微分で第 $j$ 成分は $X$ の各因子への引き戻し $\sum_\ell\operatorname{pr}_\ell^*\omega_j$ である。標準体積の引き戻しを外積で展開し、各因子に $(1,0)$ と $(0,1)$ の形式が一つずつ現れる項だけを取る。$(1,0)$ 形式を各因子へ割り当てる置換は $g!$ 通りである。その配置を一つ固定すると、$(0,1)$ 形式の割り当ての置換が符号付きの積 $\prod_j A_{j,\sigma(j)}$ を与える。この和は $\det A$ であり、他の $(1,0)$ 配置でも同じ値になる。よって $\int_{X^g}\Phi^*\nu=g!\det A$ である。積分の比が写像の次数になることは [写像の次数](/math/differential-topology/degree/) を参照する。

Sard の定理で正則値 $u$ を取る。その逆像の各点で $\Phi$ の局所実次数は $+1$ であり、従って逆像は $g!$ 点である。二つの成分が一致する順序付き $g$ 組では微分の対応する二列が比例するため臨界点になる。従って正則値の逆像で重複はなく、$g!$ 個の並べ替えが一組の有効因子を与える。この意味で Jacobi の逆問題の解は一般に一意である。

:::example[種数二の逆像]
$g=2$ で正則値 $u$ を取ると、順序付き写像の逆像は $(p,q)$ と $(q,p)$ である。$p=q$ が起こる値では行列式が零になるので、この数え方の仮定を満たさない。特殊因子が含まれると、同じ Abel 像を持つ有効因子が正次元の族になる場合がある。
:::

### 対称積の座標と重複点での微分

$X^g$ の順序を忘れた商を $\operatorname{Sym}^gX$ と書く。一点の近くに $k$ 点が集まる場合には、その局所座標の値 $z_1,\ldots,z_k$ に対して基本対称式を使う。これらはモニック多項式 $\prod_j(T-z_j)$ の係数であり、根を重複込みで一意に復元できる。根が指定円板に入る条件は境界での Rouché の定理により係数に関して開である。このため係数は重複点でも座標になる。異なる点のクラスターには互いに素な円板を取り、これらの座標を組み合わせる。

:::example[二点が重なったときに残る微分]
局所原始関数が $F(z)=c_0+c_1z+c_2z^2+O(z^3)$ とする。対称積の座標 $e_1=z_1+z_2,e_2=z_1z_2$ では
$$
F(z_1)+F(z_2)=2c_0+c_1e_1+c_2(e_1^2-2e_2)+O((|e_1|+|e_2|)^2)
$$
と一次まで書ける。三次以上のべき和は漸化式 $p_k=e_1p_{k-1}-e_2p_{k-2}$ から一次項を持たない。従って $(e_1,e_2)=(0,0)$ での微分は $c_1\,\mathrm de_1-2c_2\,\mathrm de_2$ である。順序付きの積で二列が等しくなることから、対称積の写像の微分も特異だと結論してはならない。
:::

:::example[正の次元を持つ逆像]
種数 $2$ の超楕円曲線では $\operatorname{Sym}^2X$ 上のファイバー因子 $[p]+[\iota p]$ が $x$ の値ごとに変わるが、全て同じ Abel 像を持つ。実際、二つのファイバーの差は $(x-a)/(x-b)$ の主因子である。従って一般の逆像の一意性と、この一点の上での無限個の解は両立する。
:::

### 多重点を含む因子で Abel 写像を微分する

有効因子 $D=\sum_i k_i[p_i]$ において、各 $p_i$ に互いに素な座標円板を取る。各円板内の $k_i$ 点の基本対称式 $e_1,\ldots,e_{k_i}$ は対称積の局所座標になる。実際、モニック式の係数から重複度込みの根の集合が一意に決まり、根が円板内に留まることは境界上の Rouché の定理で開条件になる。局所的な Abel 積分の和は、係数に依存して正則である。原始関数 $F$ と、すべての根を囲む円上のモニック式 $P$ を用いると、その和は $(2\pi i)^{-1}\int F(\zeta)P'(\zeta)/P(\zeta)\,\dd\zeta$ だからである。

ある円板で $\omega=f(z)\dd z=(\sum_{m\ge0}a_m z^m)\dd z$ とする。原始関数の値の和は $\sum_{m\ge1}(a_{m-1}/m)s_m$ であり、$s_m=\sum_j z_j^m$ は根のべき和である。Newton の恒等式から、$e=0$ における $s_m$ の一次部分は $m\le k_i$ なら $(-1)^{m-1}m e_m$、$m>k_i$ なら零である。級数は根を囲む固定小円上で一様収束するため、積分表示を通して項ごとに微分できる。従って $\omega$ で評価した微分のこの円板での成分は $\sum_{m=1}^{k_i}(-1)^{m-1}a_{m-1}\dd e_m$ である。

よって Abel 写像の微分を双対に移したときの核は、全ての $p_i$ で $k_i$ 位以上消える微分、すなわち $H^0(X,K_X(-D))$ である。$\deg D=g$ なら Riemann–Roch によりこの次元は $\ell(D)-1$ となるため、対称積での微分が全階数となる条件は $\ell(D)=1$ である。この判定には、順序付きの積で同じ点が重なったことだけでは足りない。

### 特殊因子の逆像の次元

有効因子 $D$ に線束 $L=\mathcal O(D)$ を対応させる。同じ次数と Abel 像を持つ有効因子は Abel の定理からちょうど線形同値なものに等しい。任意の非零切断 $s$ は零因子を与え、同じ零因子を与える二切断の比は定数である。逆に $D'$ が $D$ と線形同値なら、$D'-D=(f)$ となる $f$ が $L$ の切断を与え、その零因子は $D'$ になる。従って集合としての逆像は $\mathbb P H^0(X,L)$ と同一視され、$\ell(D)-1$ 個の複素パラメータを持つ。

:::example[種数三で二次元の空間が逆像に現れる]
超楕円曲線 $y^2=P(x)$ を種数 $3$ とし、$D$ を $x$ のファイバーと一点 $p$ の和とする。$1,x$ を適切な共通分母で表せば、ファイバーを動かす二次元の切断空間が $L(D)$ に入る。従って次数 $3$ の Abel 写像の逆像には少なくとも一次元の射影直線が含まれる。固定した $p$ と、任意の $x=a$ のファイバーを加えた因子族がその直線を実現する。
:::

:::example[種数二で重複する二点の微分]
種数二の $y^2=P(x)$ の非分岐点 $p$ で $z=x-x(p)$ とする。微分 $h(z)\mathrm dz=\mathrm dx/y$、$(z+a)h(z)\mathrm dz=x\,\mathrm dx/y$ を使う。対称積の $2p$ での微分の二列は $(h(0),a h(0))^\top$ と $-(h'(0),h(0)+ah'(0))^\top$ である。行列式は $-h(0)^2$ で非零である。これに対し分岐点では局所座標が $x-a=t^2$ であり、第二微分の $a$ 倍を引いた差は位数 $2$ なので、接空間に対応する一次の係数には寄与しない。
:::

### 逆問題で指定するトーラス点の座標

:::example[周期格子で差の位数を求める]
次の周期行列を持つトーラスで、二つの代表ベクトル $p,q$ の差が周期かどうかを判定する。Abel 積分の計算が済んだ後には、この判定だけが必要になる。

$$
\begin{aligned}
\tau&=\left[\begin{matrix}1 + 4 i & 2 + i\\2 + i & 2 i\end{matrix}\right] \\
\Im\tau&=\left[\begin{matrix}4 & 1\\1 & 2\end{matrix}\right] \\
p&=\left[\begin{matrix}2 + 2 i\\3 + 4 i\end{matrix}\right],\qquad q=\left[\begin{matrix}1\\1\end{matrix}\right] \\
p-q&=\left[\begin{matrix}1 + 2 i\\2 + 4 i\end{matrix}\right]
\end{aligned}
$$

$p-q=m+\tau n$ と書く。虚部の方程式 $Bn=\Im(p-q)$ を、正定値の $B$ の逆行列で解く。その結果を実部の方程式に入れて $m$ を得る。

$$
\begin{aligned}
\det B&=7 \\
n&=B^{-1}\Im(p-q)=\left[\begin{matrix}0\\2\end{matrix}\right] \\
m&=\Re(p-q)-(\Re\tau)n=\left[\begin{matrix}-3\\2\end{matrix}\right]
\end{aligned}
$$

格子に入る条件は全成分が整数となることに等しい。今回の全成分の分母の最小公倍数は $1$ なので、トーラス上の差の位数は $1$ である。

このデータが次数の等しい二因子から得られた場合には、この倍数だけ取った差が初めて主因子になる。
:::

:::example[周期格子で差の位数を求める]
次の周期行列を持つトーラスで、二つの代表ベクトル $p,q$ の差が周期かどうかを判定する。Abel 積分の計算が済んだ後には、この判定だけが必要になる。

$$
\begin{aligned}
\tau&=\left[\begin{matrix}3 i & \frac{1}{2} + i\\\frac{1}{2} + i & 5 i\end{matrix}\right] \\
\Im\tau&=\left[\begin{matrix}3 & 1\\1 & 5\end{matrix}\right] \\
p&=\left[\begin{matrix}1 + i\\2 + i\end{matrix}\right],\qquad q=\left[\begin{matrix}1\\2\end{matrix}\right] \\
p-q&=\left[\begin{matrix}i\\i\end{matrix}\right]
\end{aligned}
$$

$p-q=m+\tau n$ と書く。虚部の方程式 $Bn=\Im(p-q)$ を、正定値の $B$ の逆行列で解く。その結果を実部の方程式に入れて $m$ を得る。

$$
\begin{aligned}
\det B&=14 \\
n&=B^{-1}\Im(p-q)=\left[\begin{matrix}\frac{2}{7}\\\frac{1}{7}\end{matrix}\right] \\
m&=\Re(p-q)-(\Re\tau)n=\left[\begin{matrix}- \frac{1}{14}\\- \frac{1}{7}\end{matrix}\right]
\end{aligned}
$$

格子に入る条件は全成分が整数となることに等しい。今回の全成分の分母の最小公倍数は $14$ なので、トーラス上の差の位数は $14$ である。

このデータが次数の等しい二因子から得られた場合には、この倍数だけ取った差が初めて主因子になる。
:::

## Jacobi 多様体

:::definition[Jacobi 多様体]{#def-cg-jacobian}
:term[Jacobi 多様体]{en="Jacobian variety" yomi="やこびたようたい"} は $J(X)=\C^g/(\Z^g+\tau\Z^g)$ である。座標に依存しない形では、$H^0(K_X)^*$ を $H_1(X,\Z)$ の周期汎関数の像で割ったものである。
:::

:ref[thm-cg-abel] と :ref[thm-cg-jacobi-inversion] により次数零因子類の群は $J(X)$ と同型である。:ref[thm-cg-curve-rr] で全直線束が因子に由来すると示したので、これは次数零部分の Picard 群でもある。格子上の交代形式
$$
E(m+\tau n,m'+\tau n')=m^\top n'-n^\top m'
$$
は交点形式に一致する。$Y=\Im\tau$ とし、第一変数で線形な Hermite 形式 $H(u,v)=u^\top Y^{-1}\bar v$ を取る。$u=m+\tau n$、$v=m'+\tau n'$ とすれば $-\Im H(u,v)=E(u,v)$ である。実部の交差項を展開すると $\Re\tau$ の対称性で余分な項が相殺する。この Hermite 形式は正定値で、交代形式は整数かつユニモジュラーである。この組を主偏極という。

:::example[基点付き曲線の Abel 写像]
$g\ge1$ では $p\mapsto A([p]-[p_0])$ が単射である。二点の像が一致すれば差は主因子になり、その非自明な場合は球面への次数一写像が生じるからである。また正則微分の全てが同じ点で消えることはない。これは $h^0(K_X(-p))=g-1$ を Riemann–Roch で確認すれば従う。そのため微分も単射であり、コンパクト性から埋め込みである。
:::

### テータ関数で Abel 像を調べる

正定値な $\Im\tau$ を持つ対称行列 $\tau$ について
$$
\vartheta(z,\tau)=\sum_{m\in\Z^g}
  \exp(\pi i m^\top\tau m+2\pi i m^\top z)
$$
を :term[Riemann のテータ関数]{en="Riemann theta function" yomi="りーまんのてーたかんすう"} と呼ぶ。$\Im\tau$ の最小固有値を $c>0$ とすると、コンパクトな $z$ 集合上で各項の絶対値は $\exp(-\pi c|m|^2+C|m|)$ 以下である。任意階の微分で生じる $m$ の多項式因子を掛けても総和は収束する。従って正則関数であり、格子平行移動に対して
$$
\begin{aligned}
\vartheta(z+a,\tau)&=\vartheta(z,\tau),\\
\vartheta(z+\tau b,\tau)&=
  e^{-\pi i b^\top\tau b-2\pi i b^\top z}\vartheta(z,\tau)
\end{aligned}
$$
を満たす。ここで $a,b\in\Z^g$ である。最初の式は指数に $2\pi i m^\top a$ を加え、次の式は $m+b$ で和を取り直せば得られる。指数の因子は非零なので零点集合はトーラスへ降りる。また実座標の $[0,1]^g$ で $\vartheta$ を積分すると $m=0$ の項だけが残って $1$ となるため、恒等零ではない。

Riemann の消滅定理は $g\ge1$ の曲線の Jacobi 多様体では、この零因子の台が $X^{g-1}$ の Abel 像の平行移動に等しいことをいう。対応する次数 $g-1$ の直線束 $L$ における零の重複度が $h^0(X,L)$ に等しいというのが Riemann の特異点定理である。これらには多角形上でのテータ関数の零点の計数と因子の変形の解析が必要である。この章では周期変換の導出までを行い、両定理の証明は :cite[griffiths-harris] を参照する。

:::example[種数一で見えるテータの零]
$g=1$ で $z=(1+\tau)/2$ とすると、第 $m$ 項は $(-1)^m\exp(\pi i\tau m(m+1))$ となる。$m$ と $-m-1$ の項は絶対収束和の中で相殺するため、ここは零点である。種数一では次数零の有効因子は零因子だけであり、Abel 像の平行移動で一つの点が得られるという一般の定理と両立する。
:::

### テータの因子の整合性を確かめる

$\lambda=a+\tau b$ に対して $j_\lambda(z)=\exp(-\pi i b^\top\tau b-2\pi i b^\top z)$ と書く。$\lambda'=a'+\tau b'$ について
$$
j_{\lambda'}(z+\lambda)j_\lambda(z)=j_{\lambda+\lambda'}(z)
$$
である。実際、二次項の交差項 $-2\pi i {b'}^\top\tau b$ がちょうど現れ、余分な項は $-2\pi i {b'}^\top a$ で、その指数は $1$ である。従って $(z,v)\mapsto(z+\lambda,j_\lambda(z)v)$ が格子の作用となり、商は Jacobi 多様体上の直線束である。テータ関数はその切断を与える。格子平行移動に対し値が変わることは、切断として降りることと矛盾しない。

:::example[種数一の基本平行四辺形で零を数える]
平行四辺形の境界に零が来ないように平行移動し、$h=\vartheta'/\vartheta$ を使う。変換則より $h(z+1)=h(z)$、$h(z+\tau)=h(z)-2\pi i$ である。横の辺の差は $2\pi i$、斜めの二辺は相殺する。従って $\oint h\,\mathrm dz/(2\pi i)=1$ であり、テータの零点は重複度込みで一つである。既に見つけた $(1+\tau)/2$ は単純零点である。
:::

### Jacobi トーラスの有限位数の点

:::example[周期格子で差の位数を求める]
次の周期行列を持つトーラスで、二つの代表ベクトル $p,q$ の差が周期かどうかを判定する。Abel 積分の計算が済んだ後には、この判定だけが必要になる。

$$
\begin{aligned}
\tau&=\left[\begin{matrix}3 i & 1 + i\\1 + i & 2 i\end{matrix}\right] \\
\Im\tau&=\left[\begin{matrix}3 & 1\\1 & 2\end{matrix}\right] \\
p&=\left[\begin{matrix}3 + 2 i\\2 + i\end{matrix}\right],\qquad q=\left[\begin{matrix}0\\0\end{matrix}\right] \\
p-q&=\left[\begin{matrix}3 + 2 i\\2 + i\end{matrix}\right]
\end{aligned}
$$

$p-q=m+\tau n$ と書く。虚部の方程式 $Bn=\Im(p-q)$ を、正定値の $B$ の逆行列で解く。その結果を実部の方程式に入れて $m$ を得る。

$$
\begin{aligned}
\det B&=5 \\
n&=B^{-1}\Im(p-q)=\left[\begin{matrix}\frac{3}{5}\\\frac{1}{5}\end{matrix}\right] \\
m&=\Re(p-q)-(\Re\tau)n=\left[\begin{matrix}\frac{14}{5}\\\frac{7}{5}\end{matrix}\right]
\end{aligned}
$$

格子に入る条件は全成分が整数となることに等しい。今回の全成分の分母の最小公倍数は $5$ なので、トーラス上の差の位数は $5$ である。

このデータが次数の等しい二因子から得られた場合には、この倍数だけ取った差が初めて主因子になる。
:::

:::example[周期格子で差の位数を求める]
次の周期行列を持つトーラスで、二つの代表ベクトル $p,q$ の差が周期かどうかを判定する。Abel 積分の計算が済んだ後には、この判定だけが必要になる。

$$
\begin{aligned}
\tau&=\left[\begin{matrix}1 + 2 i & i\\i & 3 i\end{matrix}\right] \\
\Im\tau&=\left[\begin{matrix}2 & 1\\1 & 3\end{matrix}\right] \\
p&=\left[\begin{matrix}1 + 3 i\\2 + 4 i\end{matrix}\right],\qquad q=\left[\begin{matrix}1\\2\end{matrix}\right] \\
p-q&=\left[\begin{matrix}3 i\\4 i\end{matrix}\right]
\end{aligned}
$$

$p-q=m+\tau n$ と書く。虚部の方程式 $Bn=\Im(p-q)$ を、正定値の $B$ の逆行列で解く。その結果を実部の方程式に入れて $m$ を得る。

$$
\begin{aligned}
\det B&=5 \\
n&=B^{-1}\Im(p-q)=\left[\begin{matrix}1\\1\end{matrix}\right] \\
m&=\Re(p-q)-(\Re\tau)n=\left[\begin{matrix}-1\\0\end{matrix}\right]
\end{aligned}
$$

格子に入る条件は全成分が整数となることに等しい。今回の全成分の分母の最小公倍数は $1$ なので、トーラス上の差の位数は $1$ である。

このデータが次数の等しい二因子から得られた場合には、この倍数だけ取った差が初めて主因子になる。
:::

### 因子の和とトーラスの群法則

Abel 写像は次数零の因子の加法からトーラスの加法への準同型である。基点を変えたときの補正は因子の次数を掛けた定数なので消える。よって次数零の線束のテンソル積は、対応する二つのトーラス点の和になる。例えば平方格子上で $p=1/3+i/2$ とすると、$2p=2/3+i$、$3p=1+3i/2$ はいずれも格子に入らない。$6p=2+3i$ は格子に入るため、$\mathcal O([p]-[0])$ の位数は $6$ である。次数零という条件だけでは、どの正のテンソルべきが自明になるかは決まらない。

## Weierstrass 点と超楕円曲線

### 極位数の欠番

点 $p$ を固定し $\ell(n)=h^0(\mathcal O_X(n[p]))$ と書く。増分が $0$ または $1$ であることは一点の完全列から従う。増分が $0$ なら、$p$ だけにちょうど $n$ 位の極を持つ関数がなく、$n$ を :term[欠番]{en="gap" yomi="けつばん"} と呼ぶ。$n\ge2g$ では $n-1>2g-2$ なので Riemann–Roch は増分 $1$ を与える。さらに $g\ge1$ なら $\ell(2g-1)=g$、$\ell(0)=1$ なので $1,\ldots,2g-1$ にちょうど $g$ 個の欠番がある。

Riemann–Roch の差を取ると
$$
\ell(n)-\ell(n-1)
=1+h^0(K_X(-np))-h^0(K_X(-(n-1)p)).
$$
したがって欠番は正則微分の可能な消失位数に $1$ を加えたものに等しい。多くの点で欠番は $1,\ldots,g$ である。この列から外れる点を :term[Weierstrass 点]{en="Weierstrass point" yomi="わいえすとらすてん"} と呼ぶ。

:::proposition[重みの総和]{#prop-cg-weierstrass-weight}
$g\ge2$ として、$p$ での欠番を $n_1<\cdots<n_g$ とする。重み $w(p)=\sum_j(n_j-j)$ は非負であり、$\sum_p w(p)=g^3-g$ である。
:::

:::proof
局所座標で基底を $\omega_j=f_j(z)\dd z$ と書き、Wronski 行列式 $W=\det(f_j^{(i-1)})_{i,j=1}^g$ を取る。別の座標 $z=z(t)$ では係数 $f_j(z(t))z'(t)$ を $i-1$ 回微分する。最高階の項は $(z')^i f_j^{(i-1)}$ で、残りは前の階の行の線形結合である。従って $W(\dd z)^{g(g+1)/2}$ は大域切断をなす。

$p$ での係数の位数が $a_j=n_j-1$ になる基底を、消失位数のフィルトレーションに合わせて取る。行列式の最初の候補次数は $\sum_j a_j-g(g-1)/2$ で、その係数は各 $f_j$ の初項係数と $\det((a_j)_{i-1})$ の積である。下降階乗 $(a)_k=a(a-1)\cdots(a-k+1)$ は次数 $k$ のモニック多項式だから、この行列式は Vandermonde 行列式で非零となる。従って位数は候補次数に等しく $w(p)$ である。切断の零因子の次数は $g(g+1)\deg K_X/2=g^3-g$ になる。
:::

:::example[超楕円曲線]
平方自由な次数 $2g+1$ または $2g+2$ の $P$ を用いて $y^2=P(x)$ とする。正則微分の基底は $x^j\dd x/y$、$0\le j<g$ である（:ref[prop-cpx-hyperelliptic-forms]）。有限分岐点 $x=a$ では $(x-a)^j\dd x/y$ の位数が $2j$ なので、欠番は $1,3,\ldots,2g-1$、重みは $g(g-1)/2$ である。無限遠が分岐する場合は $x$ が位数二の極であり、$\dd x/y$ の位数は $2g-2$ なので、基底の位数を逆順にすれば同じ列になる。全 $2g+2$ 分岐点の重みを足すと $(2g+2)g(g-1)/2=g^3-g$ である。よって $g\ge2$ ではこの曲線の Weierstrass 点は分岐点ちょうどである。
:::

### 標準微分の消失列と欠番

$p$ を中心とする座標で、正則微分の部分空間 $H^0(K-kp)$ は $p$ で $k$ 位以上消えるものからなる。各段階で高々一つ次元が落ち、落ちる段階の位数を $a_0<\cdots<a_{g-1}$ とする。Riemann–Roch を $kp$ と $(k-1)p$ について引き算すると
$$
\ell(kp)-\ell((k-1)p)
=1+h^0(K-kp)-h^0(K-(k-1)p)
$$
となる。従って欠番は $a_j+1$ であり、重みは $\sum_j(a_j-j)$ である。基底の選び方を変えても、消失で定めた部分空間の次元が同じなので欠番は変わらない。

:::example[分岐点とそうでない点を比較する]
$y^2=P(x)$、$\deg P=2g+2$ で平方自由とする。有限の非分岐点では $(x-a)^j\mathrm dx/y$ の位数は $0,1,\ldots,g-1$ である。一方、単純分岐点では $x-a=t^2$、$\mathrm dx/y$ は単元倍の $\mathrm dt$ なので位数は $0,2,\ldots,2g-2$ となる。$g=3$ では非分岐点の欠番は $1,2,3$ で重み零、分岐点の欠番は $1,3,5$ で重み $3$ である。八分岐点から $24$ となり、$g^3-g$ と一致する。
:::

### 超楕円曲線以外の欠番を求める

:::example[3 枚の根の被覆を完成する]
$y^3=x \left(x - 4\right) \left(x - 3\right) \left(x - 2\right) \left(x - 1\right)$ の根の被覆を取り、完成で加わる点と各点での局所次数を調べる。

根に $e^{2\pi i/3}$ を掛けて番号を付ける。指数 $m$ の点では、非消失因子の $3$ 乗根を選んで除く。一周は番号を $m$ だけ進め、$d=\gcd(3,m)$ 個の長さ $e=3/d$ の軌道に分ける。無限遠では $s=1/x$ として指数の和の負を使う。

$$
\begin{aligned}
0\colon\quad m=1,\quad d=\gcd(3,1)=1,\quad e=3,\quad R=2 \\
1\colon\quad m=1,\quad d=\gcd(3,1)=1,\quad e=3,\quad R=2 \\
2\colon\quad m=1,\quad d=\gcd(3,1)=1,\quad e=3,\quad R=2 \\
3\colon\quad m=1,\quad d=\gcd(3,1)=1,\quad e=3,\quad R=2 \\
4\colon\quad m=1,\quad d=\gcd(3,1)=1,\quad e=3,\quad R=2 \\
\infty\colon\quad m=-5,\quad d=\gcd(3,-5)=1,\quad e=3,\quad R=2
\end{aligned}
$$

完成では各軌道に一点を加え、底の中心座標を $t^e$ とする。全体は $1$ 成分であり、番号を一つ進める置換は各成分を入れ替える。従って各成分の次数は $3$、分岐寄与は全寄与を $1$ で割った $12$ である。分岐値ごとの補正は $d$ 個の逆像を合わせた $d(e-1)$ である。

$$
\begin{aligned}
2g-2&=-2\cdot 3+12=6 \\
g&=4
\end{aligned}
$$

無限遠の一点 $q$ で $x,y$ の極位数は $3,5$ であり、有限部分に極を持たない。次は $9$ 以下の極位数を実現する単項式である。

$$
\begin{aligned}
0&:\quad 1 \\
3&:\quad x \\
5&:\quad y \\
6&:\quad x^{2} \\
8&:\quad x y \\
9&:\quad x^{3}
\end{aligned}
$$

候補の欠番は $1,2,4,7$ であり、個数は種数 $4$ に等しい。これ以外の増分を全て実現したため、候補が実際の欠番だと分かる。
:::

:::example[4 枚の根の被覆を完成する]
$y^4=x \left(x - 2\right) \left(x - 1\right)$ の根の被覆を取り、完成で加わる点と各点での局所次数を調べる。

根に $e^{2\pi i/4}$ を掛けて番号を付ける。指数 $m$ の点では、非消失因子の $4$ 乗根を選んで除く。一周は番号を $m$ だけ進め、$d=\gcd(4,m)$ 個の長さ $e=4/d$ の軌道に分ける。無限遠では $s=1/x$ として指数の和の負を使う。

$$
\begin{aligned}
0\colon\quad m=1,\quad d=\gcd(4,1)=1,\quad e=4,\quad R=3 \\
1\colon\quad m=1,\quad d=\gcd(4,1)=1,\quad e=4,\quad R=3 \\
2\colon\quad m=1,\quad d=\gcd(4,1)=1,\quad e=4,\quad R=3 \\
\infty\colon\quad m=-3,\quad d=\gcd(4,-3)=1,\quad e=4,\quad R=3
\end{aligned}
$$

完成では各軌道に一点を加え、底の中心座標を $t^e$ とする。全体は $1$ 成分であり、番号を一つ進める置換は各成分を入れ替える。従って各成分の次数は $4$、分岐寄与は全寄与を $1$ で割った $12$ である。分岐値ごとの補正は $d$ 個の逆像を合わせた $d(e-1)$ である。

$$
\begin{aligned}
2g-2&=-2\cdot 4+12=4 \\
g&=3
\end{aligned}
$$

無限遠の一点 $q$ で $x,y$ の極位数は $4,3$ であり、有限部分に極を持たない。次は $7$ 以下の極位数を実現する単項式である。

$$
\begin{aligned}
0&:\quad 1 \\
3&:\quad y \\
4&:\quad x \\
6&:\quad y^{2} \\
7&:\quad x y
\end{aligned}
$$

候補の欠番は $1,2,5$ であり、個数は種数 $3$ に等しい。これ以外の増分を全て実現したため、候補が実際の欠番だと分かる。
:::

### 超楕円でない曲線の欠番と重み

:::example[3 枚の根の被覆を完成する]
$y^3=x \left(x - 6\right) \left(x - 5\right) \left(x - 4\right) \left(x - 3\right) \left(x - 2\right) \left(x - 1\right)$ の根の被覆を取り、完成で加わる点と各点での局所次数を調べる。

根に $e^{2\pi i/3}$ を掛けて番号を付ける。指数 $m$ の点では、非消失因子の $3$ 乗根を選んで除く。一周は番号を $m$ だけ進め、$d=\gcd(3,m)$ 個の長さ $e=3/d$ の軌道に分ける。無限遠では $s=1/x$ として指数の和の負を使う。

$$
\begin{aligned}
0\colon\quad m=1,\quad d=\gcd(3,1)=1,\quad e=3,\quad R=2 \\
1\colon\quad m=1,\quad d=\gcd(3,1)=1,\quad e=3,\quad R=2 \\
2\colon\quad m=1,\quad d=\gcd(3,1)=1,\quad e=3,\quad R=2 \\
3\colon\quad m=1,\quad d=\gcd(3,1)=1,\quad e=3,\quad R=2 \\
4\colon\quad m=1,\quad d=\gcd(3,1)=1,\quad e=3,\quad R=2 \\
5\colon\quad m=1,\quad d=\gcd(3,1)=1,\quad e=3,\quad R=2 \\
6\colon\quad m=1,\quad d=\gcd(3,1)=1,\quad e=3,\quad R=2 \\
\infty\colon\quad m=-7,\quad d=\gcd(3,-7)=1,\quad e=3,\quad R=2
\end{aligned}
$$

完成では各軌道に一点を加え、底の中心座標を $t^e$ とする。全体は $1$ 成分であり、番号を一つ進める置換は各成分を入れ替える。従って各成分の次数は $3$、分岐寄与は全寄与を $1$ で割った $16$ である。分岐値ごとの補正は $d$ 個の逆像を合わせた $d(e-1)$ である。

$$
\begin{aligned}
2g-2&=-2\cdot 3+16=10 \\
g&=6
\end{aligned}
$$

無限遠の一点 $q$ で $x,y$ の極位数は $3,7$ であり、有限部分に極を持たない。次は $13$ 以下の極位数を実現する単項式である。

$$
\begin{aligned}
0&:\quad 1 \\
3&:\quad x \\
6&:\quad x^{2} \\
7&:\quad y \\
9&:\quad x^{3} \\
10&:\quad x y \\
12&:\quad x^{4} \\
13&:\quad x^{2} y
\end{aligned}
$$

候補の欠番は $1,2,4,5,8,11$ であり、個数は種数 $6$ に等しい。これ以外の増分を全て実現したため、候補が実際の欠番だと分かる。
:::

:::example[5 枚の根の被覆を完成する]
$y^5=x \left(x - 2\right) \left(x - 1\right)$ の根の被覆を取り、完成で加わる点と各点での局所次数を調べる。

根に $e^{2\pi i/5}$ を掛けて番号を付ける。指数 $m$ の点では、非消失因子の $5$ 乗根を選んで除く。一周は番号を $m$ だけ進め、$d=\gcd(5,m)$ 個の長さ $e=5/d$ の軌道に分ける。無限遠では $s=1/x$ として指数の和の負を使う。

$$
\begin{aligned}
0\colon\quad m=1,\quad d=\gcd(5,1)=1,\quad e=5,\quad R=4 \\
1\colon\quad m=1,\quad d=\gcd(5,1)=1,\quad e=5,\quad R=4 \\
2\colon\quad m=1,\quad d=\gcd(5,1)=1,\quad e=5,\quad R=4 \\
\infty\colon\quad m=-3,\quad d=\gcd(5,-3)=1,\quad e=5,\quad R=4
\end{aligned}
$$

完成では各軌道に一点を加え、底の中心座標を $t^e$ とする。全体は $1$ 成分であり、番号を一つ進める置換は各成分を入れ替える。従って各成分の次数は $5$、分岐寄与は全寄与を $1$ で割った $16$ である。分岐値ごとの補正は $d$ 個の逆像を合わせた $d(e-1)$ である。

$$
\begin{aligned}
2g-2&=-2\cdot 5+16=6 \\
g&=4
\end{aligned}
$$

無限遠の一点 $q$ で $x,y$ の極位数は $5,3$ であり、有限部分に極を持たない。次は $9$ 以下の極位数を実現する単項式である。

$$
\begin{aligned}
0&:\quad 1 \\
3&:\quad y \\
5&:\quad x \\
6&:\quad y^{2} \\
8&:\quad x y \\
9&:\quad y^{3}
\end{aligned}
$$

候補の欠番は $1,2,4,7$ であり、個数は種数 $4$ に等しい。これ以外の増分を全て実現したため、候補が実際の欠番だと分かる。
:::

## 楕円曲線と Weierstrass の ℘ 関数

格子 $\Lambda$ に対し $\wp$ は二位の極を原点に持つ偶楕円関数であり、$\wp'$ は三位の極を持つ奇関数である。級数による構成は :ref[def-cpx-wp]、周期性は :ref[prop-cpx-wp-period] を参照する。$g_2=60\sum_{\lambda\ne0}\lambda^{-4}$、$g_3=140\sum_{\lambda\ne0}\lambda^{-6}$ と置く。Laurent 展開は $\wp=z^{-2}+(g_2/20)z^2+(g_3/28)z^4+\cdots$ である。従って $(\wp')^2-4\wp^3+g_2\wp+g_3$ の主要部と定数項は消え、極のない楕円関数はコンパクトなトーラス上で定数である。ここでは Laurent 展開の定数項も零だから恒等的に零となる。

:::example[原点を無限遠に写す]
$[\wp(z):\wp'(z):1]$ は原点で $z^3$ を全座標に掛けて $[0:1:0]$ に延長できる。$1,\wp,\wp'$ は極位数が異なり独立であり、次数三の $\mathcal O_X(3[0])$ の切断空間の基底である。従って :ref[thm-cg-curve-embedding] によりこの写像は埋め込みである。像のアフィン式は $y^2=4x^3-g_2x-g_3$ である。例えば $g_2=4,g_3=0$ の式は有限で $y=0$ となる三つの点の $x$ 座標が $-1,0,1$ であり、右辺の単根に対応する。
:::

直線を曲線に制限した切断の零は重複度込みで三点 $p,q,r$ を与える。無限遠の直線との交点因子は $3[0]$ である。二切断の比の因子に Abel の定理を適用して $p+q+r=0$ を得る。従って弦または接線をとった第三交点を $(x,y)\mapsto(x,-y)$ で反転すると和になる。接線では接点を二回数える。

### 微分方程式を Laurent 級数と突き合わせる

$\wp(z)=z^{-2}+Az^2+Bz^4+O(z^6)$ とする。偶性と定義の差し引きにより定数項はない。二乗と三乗の先頭項を比べると
$$
(\wp')^2-4\wp^3=-20Az^{-2}-28B+O(z^2).
$$
従って $g_2=20A,g_3=28B$ を足し戻すと極が消え、残りはトーラス上の正則関数である。原点での値が零なのでその関数は恒等的に零である。格子点以外に極がないため、一点の Laurent 係数だけで大域方程式を決められる。

:::example[格子の回転が定数を消す]
平方格子では $i\Lambda=\Lambda$ である。定義級数の変数を掛け直すと $\wp(iz)=i^{-2}\wp(z)=-\wp(z)$ である。$z^4$ の係数を比べると $B=-B$、従って $g_3=0$ である。同様に $\zeta=e^{\pi i/3}$ が保つ格子では $\wp(\zeta z)=\zeta^{-2}\wp(z)$ であり、$z^2$ の係数を比べて $(\zeta^4-1)A=0$ を得るため $g_2=0$ である。
:::

:::example[半周期から三次式の根を読む]
半周期 $p\ne0$ では $-p=p$ が格子を法として成り立つ。$\wp'$ は奇かつ周期的なので $\wp'(p)=0$ となる。三つの非零半周期での $\wp$ の値は三次式 $4x^3-g_2x-g_3$ の根である。二つが同じ値なら、偶関数 $\wp$ の次数 $2$ のファイバーに二つの局所次数 $2$ が含まれてしまうので相異なる。
:::

### 加法公式から接線法を導く

非零格子外の $v$ に対して $\wp'(v)\ne0$ とする。$\wp$ は次数二で偶関数なので $\wp(z)=\wp(v)$ となるのは $z=\pm v$ に限る。有理式
$$
H(z)=-\wp(z)-\wp(v)+\frac14\left(\frac{\wp'(z)-\wp'(v)}{\wp(z)-\wp(v)}\right)^2
$$
は $z=v$ では分子も消えるため正則である。$z=-v$ では差の分母が $-\wp'(v)(z+v)+O((z+v)^2)$ であり、分子は $-2\wp'(v)+O(z+v)$ である。さらに分子の一次係数 $\wp''(-v)$ と分母の二次係数 $\wp''(-v)/2$ を比べると、商の定数項は零である。従って $H$ の主要部は $(z+v)^{-2}$ だけである。

$z=0$ では $\wp(z)=z^{-2}+O(z^2)$、$\wp'(z)=-2z^{-3}+O(z)$ なので、商は $-2/z-2\wp(v)z+O(z^2)$ である。その二乗を $4$ で割ると $z^{-2}+2\wp(v)+O(z)$ となるから、$H(0)=\wp(v)$ として延長できる。よって $H(z)-\wp(z+v)$ はトーラス上で極のない正則関数となり、$z=0$ で零なので恒等的に零となる。$z=v$ かつ $\wp'(v)\ne0$ では、商の極限は $\wp''(v)/\wp'(v)$ であり、接線の傾きになる。$\wp'(v)=0$ では $v$ は非零半周期の一つである。実際、$\wp'$ の極の総位数が $3$ で、三つの非零半周期ですでに消えるためである。接線は垂直で、$2v=0$ がトーラス上で成り立つため倍点は無限遠に移る。加法公式は二変数の有理型恒等式として成り立ち、分母が零の値では必要に応じて元の $\wp(z+v)$ への延長で解釈する。既存の公式と導出は :ref[thm-cpx-wp-addition] を参照できる。

:::example[三次式の三つの交点を求める]
$Y^2=4X^3-4X$ 上の $(0,0)$ と $(2,\sqrt{24})$ を結ぶ直線の傾きは $\sqrt6$ である。$Y=\sqrt6 X$ を代入すると $4X^3-6X^2-4X=2X(X-2)(2X+1)$ であり、第三交点の $X$ は $-1/2$ である。加法公式の右辺は $-0-2+6/4=-1/2$ となる。和に対応するのは第三交点の $Y$ の符号を反転した点であり、$X$ 座標だけを求める段階では符号は見えない。
:::

### 楕円トーラスの自己写像を周期で計算する

:::example[格子の整数行列と核を照合する]
$E=\C/(\Z+i\Z)$ とし、$E^1$ に次の複素線形行列を作用させる。成分が Gaussian 整数なので、像格子は同じ格子に含まれる。

$$
\begin{aligned}
M&=\left[\begin{matrix}2 + 3 i\end{matrix}\right] \\
\det_\C M&=2 + 3 i
\end{aligned}
$$

実座標の順序を実部すべて、虚部すべてとする。成分を行列の列に並べると、実格子上の作用は次の整数行列になる。

$$
\begin{aligned}
A&=\left[\begin{matrix}2 & -3\\3 & 2\end{matrix}\right] \\
A^{-1}&=\left[\begin{matrix}\frac{2}{13} & \frac{3}{13}\\- \frac{3}{13} & \frac{2}{13}\end{matrix}\right] \\
\det A&=13
\end{aligned}
$$

核は $A^{-1}\Z^2/\Z^2$ であり、$A$ を作用させて $\Z^2/A\Z^2$ と同一視する。整数の行・列操作の Smith 標準形を取ると核の巡回因子の位数が分かる。

$$
\begin{aligned}
\operatorname{Smith}(A)&=\diag(1,13) \\
\#\Ker f&=13=\abs{\det_\C M}^2
\end{aligned}
$$

像の点の一つの逆像に核を加えれば全逆像を得られる。基本平行四辺形の面積は $|2+3i|^2=13$ 倍になり、局所同型のファイバーを数えた結果と一致する。
:::

## 演習問題

::::exercise{level=basic}
$\Lambda=3\Z+(1+2i)\Z$ のトーラスで $\tau$ と正規化微分の面積積分 $(i/2)\int\omega\wedge\bar\omega$ を求めよ。
:::solution
$\omega=\dd z/3$ なので $\tau=(1+2i)/3$ である。双線形公式から面積積分は $\Im\tau=2/3$ になる。基本平行四辺形の面積 $6$ に計量の倍率 $1/9$ を掛けても $2/3$ になる。
:::
::::

::::exercise{level=basic}
平方格子トーラスで $[1/2]-[0]$ は何倍したとき初めて主因子になるか。
:::solution
Abel 像は $1/2$ であり、元の格子には入らないが二倍の $1$ は入る。従って Abel の定理で最小の正整数は $2$ である。
:::
::::

::::exercise{level=basic}
種数 $3$ で欠番列が $1,2,5$ の点の重みを求めよ。
:::solution
重みは $(1-1)+(2-2)+(5-3)=2$ である。消失位数の列は $0,1,4$ である。
:::
::::

::::exercise{level=basic}
$y^2=x^7-1$ の無限遠で基底微分の位数と欠番列を求めよ。
:::solution
この面の種数は $3$ である。$x$ は二位、$y$ は七位の極なので、$x^j\dd x/y$ の位数は $4-2j$、$j=0,1,2$ となる。昇順にすると $0,2,4$、欠番は $1,3,5$ である。
:::
::::

::::exercise{level=advanced}
$x+\tau y$ が格子点になる条件が $x,y\in\Z^g$ であることを $x,y\in\R^g$ の範囲で示せ。
:::solution
格子点なら $x+\tau y=m+\tau n$ と書ける。虚部をとると $\Im\tau(y-n)=0$ となり、正定値性から $y=n$ である。次に実部で $x=m$ が従う。逆は格子の定義による。
:::
::::

::::exercise{level=advanced}
種数 $2$ で $p$ の欠番が $1,3$ であることと、$p$ だけに二位の極を持つ関数があることの同値を示せ。
:::solution
$1$ は次数一写像がないため常に欠番であり、欠番は $1,2,3$ のうち二つである。二位の極を持つ関数がある条件は $2$ が欠番でないことなので、残る欠番は $3$ になる。逆向きも同じ増分の定義に従う。
:::
::::

::::exercise{level=advanced}
種数 $1$ の面で $[p]+[q]-2[0]$ が主因子となる条件を求めよ。
:::solution
低種数の分類で原点 $0$ を持つトーラスに同一視すると、Abel 像は $p+q$ である。よって $p+q=0$ が必要十分である。$p\not=0$ で $q=-p$ なら $\wp(z)-\wp(p)$ はその因子をもつ。偶性から二点で消え、$p=-p$ では対称性から導関数も零となる。極の総位数が $2$ だからそれ以外の零点はない。
:::
::::

::::exercise{level=advanced}
$y^2=4x^3-4x$ で $p=(0,0)$ と $q=(1,0)$ の和を求めよ。
:::solution
二点の弦は $y=0$ で、交点の $x$ 座標は $x^3-x=0$ の三根である。第三点は $(-1,0)$。反転でも変わらないため $p+q=(-1,0)$ となる。
:::
::::

::::exercise{level=advanced}
$X:y^2=x^9-1$ の無限遠の Weierstrass 重みを求めよ。微分の消失列を使った計算と、極位数の半群を使った計算を照合せよ。

:::solution
有限の九根は単根であり、無限遠も二枚被覆の分岐値なので、分岐寄与は $10$ である。従って $2g-2=-4+10=6$、$g=4$ となる。無限遠では $x=t^{-2}$、$y=t^{-9}u$、$u^2=1-t^{18}$ と書ける。よって
$$
x^j\mathrm dx/y=-2t^{6-2j}u^{-1}\mathrm dt
$$
であり、$j=0,1,2,3$ の四形式は相異なる位数 $6,4,2,0$ を持つ。有限分岐でも正則であり、種数から四本で基底となる。昇順の消失列は $0,2,4,6$、欠番列は $1,3,5,7$、重みは $(0-0)+(2-1)+(4-2)+(6-3)=6$ である。

一方 $x$ の累乗から $2,4,6,8$ の極を実現できる。標準束の次数は $6$ なので $7$ 以下の正整数にちょうど四つ欠番がある。$2,4,6$ が欠番でないことから、残り $1,3,5,7$ が全欠番である。この結論は微分から得たものと一致する。
:::
::::

## 確認テスト

::::quiz{answer="24" tol="1e-6"}
種数 $3$ の面の Weierstrass 重みの総和を求めよ。
:::explanation
Wronski 切断は $K_X^{\otimes6}$ の切断であり $\deg K_X=4$ である。従って零位数の総和は $6\cdot4=24$ となる。
:::
::::

::::quiz
正規化した周期行列の正定値性を与える積分を選べ。

- [ ] $\int_X\omega\wedge\omega$
- [x] $(i/2)\int_X\omega\wedge\bar\omega$
- [ ] $\int_X\dd\omega$

:::explanation
局所表示 $\omega=f\dd z$ では $(i/2)\omega\wedge\bar\omega=|f|^2\dd x\wedge\dd y$ である。非零の $\omega$ に対して積分は正であり、双線形公式で $\Im\tau$ の二次形式になる。
:::
::::

::::quiz{answer="2" tol="1e-6"}
種数 $2$ の面の Jacobi 多様体の複素次元を求めよ。
:::explanation
正則微分の空間の次元は種数に等しい。Jacobi 多様体はその双対の $\C^2$ を離散な周期格子で割ったものなので、局所座標の複素次元は $2$ である。
:::
::::

::::quiz
次数零因子 $D$ に関して Abel の定理が主張する同値を選べ。

- [ ] $D$ が有効であることと $A(D)=0$。
- [x] $D$ が主因子であることと $A(D)=0$。
- [ ] $D=0$ と $A(D)=0$。

:::explanation
有理型関数の対数微分の周期は $2\pi i$ の整数倍なので、その因子は零に写る。逆は、留数を $D$ の係数とする微分の周期を正則微分で補正して指数積分を取ることで示される。
:::
::::

::::quiz{answer="6" tol="1e-6"}
平方格子で $1/3+i/2$ が表す Jacobi トーラス上の点の位数を求めよ。

:::explanation
$k(1/3+i/2)$ が格子に属するには $k/3,k/2$ がともに整数である必要がある。最小の正整数は $\operatorname{lcm}(3,2)=6$ である。
:::
::::

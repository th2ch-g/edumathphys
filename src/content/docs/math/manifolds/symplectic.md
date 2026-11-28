---
title: シンプレクティック多様体
description: 非退化な閉 2 形式の局所構造を調べ、Hamilton 力学と運動量写像による簡約を導く。
notation:
  - symbol: '\omega'
    meaning: シンプレクティック形式
  - symbol: 'X_f'
    meaning: 内部積が df となる Hamilton ベクトル場
  - symbol: '\{f,g\}'
    meaning: 本章の規約で ω(X_f,X_g) となる Poisson 括弧
  - symbol: 'J:M\to\mathfrak g^*'
    meaning: Lie 群作用の運動量写像
---

## この章の目標

- 非退化な交代形式の標準基底と直交補の性質を証明する。
- 余接束の標準形式から Hamilton 方程式と Poisson 括弧を導く。
- Moser の方法で Darboux の定理を証明する。
- Lagrange 部分多様体の例と、自由固有作用の正則簡約を計算する。

## シンプレクティック線形代数

面積に似た交代双線形形式が与えられたとする。どの非零ベクトルにも、それと組み合わせて面積が非零になる相手があれば、形式の値だけでベクトルを区別できる。この性質を基に、空間を二本ずつの基底に分けられる。

:::definition[シンプレクティックベクトル空間]{#def-mfd-symplectic-vector-space}
有限次元実ベクトル空間 $V$ と交代双線形形式 $\omega$ の組で、$v\mapsto\omega(v,-)$ が $V\to V^*$ の同型であるものを :term[シンプレクティックベクトル空間]{en="symplectic vector space" yomi="しんぷれくてぃっくべくとるくうかん"} という。この条件を非退化性という。部分空間 $W\subset V$ に対し $W^\omega=\{v:\omega(v,w)=0\ \forall w\in W\}$ をシンプレクティック直交補と呼ぶ。
:::

:::proposition[シンプレクティック直交補]{#prop-mfd-symplectic-orthogonal}
$\dim W^\omega=\dim V-\dim W$、$(W^\omega)^\omega=W$ である。制限 $\omega|_W$ が非退化なら $V=W\oplus W^\omega$ で、$\omega|_{W^\omega}$ も非退化である。
:::

:::proof
$W^\omega$ は同型 $V\to V^*$ で $W$ の消滅空間に写る。$W$ の基底を $V$ の基底に延ばせば、消滅空間の次元は $\dim V-\dim W$ である。交代性から $W\subset(W^\omega)^\omega$ であり、両者が同次元なので等号である。$\omega|_W$ の核は $W\cap W^\omega$ であり、非退化なら零である。次元公式から直和になる。$W^\omega$ 上の制限の核は $W^\omega\cap(W^\omega)^\omega=W^\omega\cap W=0$ である。
:::

:::theorem[シンプレクティック基底]{#thm-mfd-symplectic-basis}
シンプレクティックベクトル空間の次元は偶数 $2n$ であり、基底 $e_1,f_1,\ldots,e_n,f_n$ を選んで $\omega(e_i,f_j)=\delta_{ij}$、$\omega(e_i,e_j)=\omega(f_i,f_j)=0$ とできる。
:::

:::proof
非零 $e_1$ を選ぶと、非退化性から $\omega(e_1,f_1)\ne0$ となる $f_1$ がある。その値で $f_1$ を割って $1$ に合わせる。$W=\spn(e_1,f_1)$ の形式は行列 $\begin{pmatrix}0&1\\-1&0\end{pmatrix}$ を持ち非退化である。:ref[prop-mfd-symplectic-orthogonal] により直交補も非退化で、次元が $2$ 減る。次元零までこの操作を繰り返せる。もし一次元で止まれば交代形式が零となり非退化に反するため、次元は偶数である。異なる段階の二本は直交するので、全ての基底ベクトルが要求された値を持つ。
:::

双対基底では $\omega=\sum_i e^i\wedge f^i$ となり、その $n$ 回の外積は $n!\,e^1\wedge f^1\wedge\cdots\wedge e^n\wedge f^n$ である。次数 $2$ の因子は交換で符号が変わらず、重複因子は零なので $n!$ が現れる。逆に $2n$ 次元で $\omega^n\ne0$ なら、もし $v\ne0$ が核に入れば $v$ を含む基底への評価がすべて零になって矛盾する。従って最高次の外積の非消滅でも非退化性を判定できる。

### 交代形式を行列で検査する

順序付き基底で交代形式 $\omega$ の行列を $\Omega_{ij}=\omega(e_i,e_j)$ とすれば、$\Omega^{\mathsf T}=-\Omega$ である。ベクトル $v$ による評価は $v^{\mathsf T}\Omega w$ なので、核は行列の左核と一致する。従って非退化性は $\det\Omega\ne0$ と同値である。奇数次では $\det\Omega=\det\Omega^{\mathsf T}=\det(-\Omega)=-\det\Omega$ とも示せるので、必ず退化する。偶数次でも全ての交代行列が可逆とは限らない。例えば四次元で $\mathrm dq^1\wedge\mathrm dp_1$ だけを取ると、二番目の座標対の全方向が核に入る。

二つのベクトル空間でそれぞれの非退化な交代形式を取ったとき、直和で両因子の評価を足した形式も非退化である。全方向に対する評価が零なら、片側だけが非零の入力に限定することで両成分がそれぞれ零になるからである。例えば $\omega=2\,\mathrm dq^1\wedge\mathrm dp_1-3\,\mathrm dq^2\wedge\mathrm dp_2$ なら、$e_i=\partial_{q^i}$、$f_1=(1/2)\partial_{p_1}$、$f_2=-(1/3)\partial_{p_2}$ が標準基底になる。また $\omega^2/2=-6\,\mathrm dq^1\wedge\mathrm dp_1\wedge\mathrm dq^2\wedge\mathrm dp_2$ であり、係数の負号は、与えた座標の向きと $\omega$ が定める向きが反対であることを示す。

:::example[等方な方向と余等方な方向を数える]
四次元標準空間の $W=\spn(\partial_{q^1},\partial_{q^2},\partial_{p_1})$ に対し、直交補を求める。$v=a_1\partial_{q^1}+a_2\partial_{q^2}+b_1\partial_{p_1}+b_2\partial_{p_2}$ とすると、$v$ と $W$ の三つの基底との値は順に $-b_1,-b_2,a_1$ である。従って $W^\omega=\R\partial_{q^2}$ であり $W^\omega\subset W$ である。この部分空間に制限した形式は退化するが、退化する方向を商で除けば二次元の非退化な空間が残る。ここと同じ線形代数を、後で運動量の等位集合の接空間に使う。
:::

### 形式と両立する複素構造と内積

$J^2=-I$ を満たす実線形写像 $J$ は、$iv=Jv$ と置くことで複素数倍を定める。$\omega(Jv,Jw)=\omega(v,w)$ かつ $\omega(v,Jv)>0$（$v\ne0$）のとき、$J$ は $\omega$ と両立するという。このとき $g_J(v,w)=\omega(v,Jw)$ は実内積である。実際、$J$ の不変性から $\omega(w,Jv)=\omega(Jw,-v)=\omega(v,Jw)$ であり、対称性が従う。正定値性は仮定そのものである。

:::proposition[両立する複素構造の構成]{#prop-mfd-compatible-complex-structure}
シンプレクティックベクトル空間 $(V,\omega)$ と任意の実内積 $g$ を取る。$g(Av,w)=\omega(v,w)$ で定めた $A$ と $S=(-A^2)^{1/2}$ に対し、$J=S^{-1}A$ は $\omega$ と両立する。
:::

:::proof
$\omega$ が交代的なので $A^*=-A$ であり、非退化性から $A$ は可逆である。従って $-A^2=A^*A$ は正定値自己共役である。その直交固有空間ごとに固有値の正の平方根を取れば、正定値な $S$ を得る。$A$ は $-A^2$ と可換なのでその各固有空間を保ち、そこでは $S$ はスカラー倍だから、$A$ は $S$ とも可換である。従って $J^2=S^{-2}A^2=-I$、$J^*J=A^*S^{-2}A=I$ である。さらに $J$ は $A$ と可換なので、$\omega(Jv,Jw)=g(JAv,Jw)=g(Av,w)$ である。最後に $\omega(v,Jv)=g(Av,S^{-1}Av)>0$ が $v\ne0$ に対して成り立つ。
:::

:::example[平面で倍率が打ち消される]
Euclid 内積で $\omega=3\,\mathrm dx\wedge\mathrm dy$ とすると $A(x,y)=3(-y,x)$、$-A^2=9I$ である。従って $J(x,y)=(-y,x)$ となる。元の内積を使って $A$ を作ったが、新しい $g_J$ は $3$ 倍の Euclid 内積になる。例えば $(1,2)$ に対し $\omega((1,2),J(1,2))=3(1+4)=15$ である。$J$ と単なる作用素 $A$ を区別するのは、二回作用させたとき $-I$ にするためである。
:::

## シンプレクティック形式と例

:::definition[シンプレクティック多様体]{#def-mfd-symplectic-manifold}
多様体 $M$ 上の 2 形式 $\omega$ が点ごとに非退化で、$\mathrm d\omega=0$ なら、$\omega$ を :term[シンプレクティック形式]{en="symplectic form" yomi="しんぷれくてぃっくけいしき"}、$(M,\omega)$ を :term[シンプレクティック多様体]{en="symplectic manifold" yomi="しんぷれくてぃっくたようたい"} という。微分同相 $F$ が $F^*\omega_N=\omega_M$ を満たすとき、シンプレクティック同相という。
:::

$M$ は偶数次元で、$\omega^n/n!$ が非消滅なので向きを持つ。例えば $\R^{2n}$ の座標 $(q^i,p_i)$ では
$$
\omega_0=\sum_i\mathrm dq^i\wedge\mathrm dp_i
\label{eq-mfd-standard-symplectic}
$$
が標準形式である。係数が定数だから閉であり、座標基底がシンプレクティック基底となる。二つのシンプレクティック多様体の積では、各形式を射影で引き戻して足すと閉である。また片方の因子の速度だけをスロットに入れれば各成分の非退化性を検査できるので、和も非退化である。

以下では境界のないシンプレクティック多様体を扱う。境界点で局所標準形を考える場合は座標域が半空間となり、時間依存の流れで境界を保つには、場が境界に接するという条件も必要になる。

:::example[曲面・高次元での閉性]
向き付けられた曲面の非消滅 2 形式は、次元により自動的に閉であり、非退化である。例えば標準の向きを持つ $S^2$ の面積形式はシンプレクティック形式である。一方、$\R^4$ の $e^{q^2}\mathrm dq^1\wedge\mathrm dp_1+\mathrm dq^2\wedge\mathrm dp_2$ は非退化であるが、外微分が $e^{q^2}\mathrm dq^2\wedge\mathrm dq^1\wedge\mathrm dp_1$ なので閉でない。高次元では二条件を個別に検査する必要がある。
:::

非空で境界なしのコンパクトな $2n$ 次元シンプレクティック多様体（$n>0$）では、$\omega$ は完全形式にならない。実際 $\omega=\mathrm d\lambda$ なら $\omega^n=\mathrm d(\lambda\wedge\omega^{n-1})$ であり、Stokes の定理で積分は零である。一方、$\omega^n$ が定める向きでは積分は正なので矛盾する。これは局所の標準形式が完全であることと、大域の完全性を区別する例である。

### 同じ次元の局所座標における保存条件

二次元標準空間の写像 $F:(q,p)\mapsto(Q,P)$ は、微分同相であるときシンプレクティック同相となる条件が $\det(\partial(Q,P)/\partial(q,p))=1$ と同値である。実際、最高次形式の引き戻しは行列式倍なので、この二形式の保存が面積の保存に等しいからである。例えば $F(q,p)=(2q,p/2)$ は面積を保つ。交換 $(q,p)\mapsto(p,q)$ は行列式 $-1$ であり、面積の絶対値は同じでも $\omega$ の符号を反転する。代わりに $(q,p)\mapsto(p,-q)$ とすれば行列式は $1$ である。

四次元以上では、最高次形式の保存だけでは $\omega$ の保存は従わない。例えば $(q^1,p_1,q^2,p_2)$ を $(2q^1,p_1,q^2,p_2/2)$ に写すと、体積係数は $2\cdot1\cdot1\cdot(1/2)=1$ である。しかし $\omega$ の引き戻しは $2\,\mathrm dq^1\wedge\mathrm dp_1+\frac12\mathrm dq^2\wedge\mathrm dp_2$ となり、元と異なる。シンプレクティック同相は、最高次体積より細かく、どの二方向の組がどの値を持つかを保つ。線形写像 $A$ の場合には、この条件を行列で $A^{\mathsf T}\Omega A=\Omega$ と書ける。この等式が成り立てば外積を $n$ 回取って $\det A=1$ が従うが、逆はこの例の通り成り立たない。

:::example[局所標準形が同じ二つのトーラス]
周期 $1$ の角 $x,y$ を持つ $T^2$ で $\omega_a=a\,\mathrm dx\wedge\mathrm dy$（$a>0$）とする。閉性と非退化性からシンプレクティック形式であり、局所では $q=ax,p=y$ と置けば $\mathrm dq\wedge\mathrm dp$ になる。しかし $a=1$ と $a=2$ の場合はシンプレクティック同相でない。もし引き戻しで一致すれば同相は二つの形式の定める向きを保ち、変数変換で全積分も等しいはずだが、それぞれ $1,2$ である。$x\mapsto2x\bmod\Z$ は局所では正しい倍率を与えても二重被覆で、逆写像を持たない。
:::

:::example[同じ向きでも補間が退化する高次元の例]
$\R^4$ の標準形式を $\omega_0$ とすると、$-\omega_0$ もシンプレクティック形式であり、最高次の外積は $(-\omega_0)^2=\omega_0^2$ なので同じ向きを与える。しかし直線補間の中点は零形式になる。二次元では同じ向きの面積形式を補間すると係数は正のままだが、高次元の非退化性にはこの凸性がない。Moser の方法で補間を使うとき、この条件を独立に検査する。
:::

## 余接束の標準シンプレクティック構造

:::definition[余接束の標準 1 形式]{#def-mfd-tautological-one-form}
余接束の射影を $\pi:T^*Q\to Q$ とし、$\alpha\in T_q^*Q$、$W\in T_\alpha(T^*Q)$ に対して
$$
\lambda_\alpha(W)=\alpha(\mathrm d\pi_\alpha W)
$$
で :term[標準 1 形式]{en="tautological one-form" yomi="ひょうじゅんいちけいしき"} を定める。標準シンプレクティック形式を $\omega=-\mathrm d\lambda$ とする。
:::

$\alpha=\sum_i p_i\,\mathrm dq^i$ と書けば、束座標 $(q,p)$ で $\lambda=\sum_i p_i\,\mathrm dq^i$ である。実際、$\mathrm d\pi$ は速度の $q$ 成分だけを取り出す。この局所式により $\lambda$ は滑らかであり、$-\mathrm d\lambda=\sum_i\mathrm dq^i\wedge\mathrm dp_i$ となる。これは閉かつ非退化である。定義に底の座標がないため、どの束座標でも局所式が両立する。

:::example[微分同相の余接束への持ち上げ]
$F:Q\to Q'$ が微分同相なら $\widehat F:T^*Q\to T^*Q'$ を $\alpha\mapsto\alpha\circ\mathrm dF^{-1}$ とする。$\pi'\widehat F=F\pi$ であり、連鎖律により
$$
(\widehat F^*\lambda')_\alpha(W)=
\alpha\bigl(\mathrm dF^{-1}\mathrm d\pi'\mathrm d\widehat F(W)\bigr)
=\alpha(\mathrm d\pi(W)).
$$
従って $\widehat F^*\lambda'=\lambda$ で、外微分を取れば標準形式を保つ。例えば実直線で $Q'=2q$ なら、余接成分 $P$ は $p/2$ に変わり、$\mathrm dQ'\wedge\mathrm dP=\mathrm dq\wedge\mathrm dp$ となる。底を伸ばす倍率と余接成分を縮める倍率が逆になる。
:::

### 余接束での平行移動と一次形式

余接束 $T^*Q$ のファイバーに沿って、滑らかな 1 形式 $\alpha$ を足す写像 $F_\alpha$ を取る。$F_\alpha(q,p)=(q,p+\alpha_q)$ であり、逆は $F_{-\alpha}$ である。射影は変わらず、標準 1 形式の定義で新しい余接成分を評価すると
$$
F_\alpha^*\lambda=\lambda+\pi^*\alpha,\qquad
F_\alpha^*\omega=\omega-\pi^*\mathrm d\alpha
$$
となる。従って $F_\alpha$ は $\alpha$ が閉なら $\omega$ を保つ。逆に保つなら $\pi^*\mathrm d\alpha=0$ であり、零切断の引き戻しを取れば $\mathrm d\alpha=0$ が従う。平行移動が閉 1 形式で定まるとき、局所原始関数から作る場合と大域原始関数を持たない場合の両方が含まれる。

例えば $Q=\R^2$ で $\alpha=y\,\mathrm dx+x\,\mathrm dy$ とすれば、$p_x$ に $y$、$p_y$ に $x$ を足す。この引き戻しでは $\mathrm dx\wedge\mathrm dy$ と $\mathrm dy\wedge\mathrm dx$ の追加項が相殺する。一方 $\alpha=x\,\mathrm dy$ だけを使えば、追加される 2 形式は $-\mathrm dx\wedge\mathrm dy$ になり、形式を保たない。同じようにファイバーに何かを足す写像でも、底方向における追加成分の交差微分が整合する必要がある。

### 底点を動かさないシンプレクティック同相を分類する

:::proposition[底を固定する余接束のシンプレクティック同相]{#prop-mfd-cotangent-base-preserving}
$\pi F=\pi$ を満たす $F:T^*Q\to T^*Q$ が標準シンプレクティック形式を保つなら、$F$ はある閉 1 形式によるファイバー平行移動である。逆も成り立つ。
:::

:::proof
束座標では $F(q,p)=(q,P(q,p))$ と書ける。$\sum_i\mathrm dq^i\wedge\mathrm dP_i=\sum_i\mathrm dq^i\wedge\mathrm dp_i$ の $\mathrm dq^i\wedge\mathrm dp_j$ の係数を比較すると $\partial_{p_j}P_i=\delta_{ij}$ である。従って連結な各ファイバーで $P_i-p_i$ は $p$ に依存せず、$P(q,0)$ に等しい。この値は $F$ を零切断に合成した切断であり、大域的な 1 形式 $\alpha$ を与える。そこで $F=F_\alpha$ と分かる。前に求めた $F_\alpha^*\omega=\omega-\pi^*\mathrm d\alpha$ と、零切断による引き戻しから、保存条件は $\mathrm d\alpha=0$ と同値である。
:::

:::example[円の余接束で保存される 2 形式と変化する 1 形式]
$T^*S^1$ を周期 $1$ の角 $q$ とその余接成分 $p$ で書く。$F(q,p)=(q,p+2)$ は $\omega$ を保つが、$F^*\lambda-\lambda=2\,\mathrm dq$ である。この差の零切断上の周回積分は $2$ なので完全でない。従って余接束には $\omega$ を保ちつつ、標準 1 形式を大域原始関数の微分では表せないだけ変える写像がある。
:::

## Hamilton ベクトル場と Poisson 括弧

非退化性によって一形式から対応する場を一意に求められる。特に関数の微分から作る場は形式を保ち、もとの関数を一定に保つ。

:::definition[Hamilton ベクトル場と Poisson 括弧]{#def-mfd-hamiltonian-field}
滑らかな関数 $f$ に対し $\iota_{X_f}\omega=\mathrm df$ を満たす場を :term[Hamilton ベクトル場]{en="Hamiltonian vector field" yomi="はみるとんべくとるば"} $X_f$ と呼ぶ。:term[Poisson 括弧]{en="Poisson bracket" yomi="ぽあそんかっこ"} を
$$
\{f,g\}=\omega(X_f,X_g)=X_gf=-X_fg
$$
で定める。
:::

束の同型 $TM\to T^*M$ の逆の成分は局所で行列の逆なので滑らかであり、$X_f$ は存在し一意である。Cartan の公式から $\mathcal L_{X_f}\omega=\mathrm d^2 f=0$ である。よってその流れは定義される範囲で $\omega$ を保つ。また $X_ff=\omega(X_f,X_f)=0$ なので、$f$ は自分の流れに沿って保存される。標準座標では
$$
X_f=\sum_i \left(\partial_{p_i}f\,\partial_{q^i}-\partial_{q^i}f\,\partial_{p_i}\right),\qquad
\{f,g\}=\sum_i\left(\partial_{q^i}f\,\partial_{p_i}g-\partial_{p_i}f\,\partial_{q^i}g\right)
\label{eq-mfd-hamilton-coordinates}
$$
であり、$\dot q^i=\partial_{p_i}f$、$\dot p_i=-\partial_{q^i}f$ が Hamilton 方程式である。これは $\sum_i\mathrm dq^i\wedge\mathrm dp_i$ の内部積の係数を $\mathrm df$ と比較した結果である。

:::proposition[Poisson 括弧の法則]{#prop-mfd-poisson-laws}
Poisson 括弧は実双線形、反対称で、各引数で積の法則を満たす。さらに Jacobi 恒等式を満たし、
$$
[X_f,X_g]=-X_{\{f,g\}}
$$
である。
:::

:::proof
双線形性と反対称性は形式と対応 $f\mapsto X_f$ から従う。積の法則は、固定した一方について他方をベクトル場で微分する定義から従う。括弧の式を示すには、$X$ の流れで $\iota_Y\omega$ を引き戻すと $\iota_{\Phi_t^*Y}\Phi_t^*\omega$ に等しいことを使う。これは点での評価と連鎖律から直ちに従う。時刻零で微分すると $\mathcal L_X\iota_Y\omega=\iota_{[X,Y]}\omega+\iota_Y\mathcal L_X\omega$ となる。$X=X_f,Y=X_g$ と取ると $\mathcal L_{X_f}\omega=0$ なので、$\iota_{[X_f,X_g]}\omega=\mathcal L_{X_f}\mathrm dg=\mathrm d(X_fg)=-\mathrm d\{f,g\}$ である。非退化性から括弧の式を得る。

この式の両側を任意の関数 $h$ に作用させると、左辺は $\{f,\{g,h\}\}-\{g,\{f,h\}\}$、右辺は $\{\{f,g\},h\}$ となる。反対称性で並べ直すと Jacobi 恒等式になる。
:::

:::example[調和振動子]
$\omega=\mathrm dq\wedge\mathrm dp$、$H=(p^2+q^2)/2$ とする。場は $p\partial_q-q\partial_p$ である。初期値 $(q_0,p_0)$ に対し
$$
q(t)=q_0\cos t+p_0\sin t,\qquad p(t)=p_0\cos t-q_0\sin t
$$
と書け、微分すれば方程式と初期条件を確認できる。初期値 $(2,1)$ は時刻 $\pi/2$ で $(1,-2)$ へ移り、$H$ は $5/2$ のままである。一般の $g$ の時間変化は $X_Hg=\{g,H\}$ であり、$\{g,H\}=0$ かどうかが保存の判定になる。
:::

### 力学での保存量と Liouville の体積

$\mathcal L_{X_H}\omega=0$ なので、外積の積公式により $\mathcal L_{X_H}(\omega^n/n!)=0$ である。従って Hamilton 流はシンプレクティック体積を保つ。標準座標でこれを直接見ると、$X_H$ の通常の発散は
$$
\sum_i\partial_{q^i}\partial_{p_i}H-\sum_i\partial_{p_i}\partial_{q^i}H=0
$$
であり、混合微分が相殺する。例えば $\R^2$ で $H=qp$ の場は $q\partial_q-p\partial_p$ である。流れは $(q,p)\mapsto(e^t q,e^{-t}p)$ となり、片方の長さが増えてももう片方が逆倍率で減るため面積が保たれる。$H$ も $(e^t q)(e^{-t}p)=qp$ で保存される。従って体積を保つことは、軌道が円のように閉じていることや、Euclid 距離を保つことを意味しない。

また $\{f,H\}=\{g,H\}=0$ なら、Jacobi 恒等式により $\{\{f,g\},H\}=0$ である。実際、
$$
\{\{f,g\},H\}=\{f,\{g,H\}\}-\{g,\{f,H\}\}
$$
の両項が零になる。従って保存量から Poisson 括弧を取って保存量を作れる。その結果が非定数とは限らない。例えば $\R^4$ で $H=(q_1^2+p_1^2)/2$ とすると、第 2 座標対 $(q_2,p_2)$ はともに保存され、その括弧は定数 $1$ である。

:::example[自由粒子と運動量]
標準座標で $H=p^2/(2m)$（$m>0$）と置くと、$\dot q=p/m$、$\dot p=0$ となる。初期値 $(q_0,p_0)$ からの流れは $(q_0+t p_0/m,p_0)$ である。$H$ が $q$ に依存しないため、$p$ の時間変化が零になる。例えば $m=2,q_0=1,p_0=6,t=3$ では $q=10$ であり、$H=36/4=9$ のままである。流れのヤコビ行列は $\begin{pmatrix}1&t/m\\0&1\end{pmatrix}$ であり、行列式が $1$ なので面積は保存される。
:::

### シンプレクティックベクトル場と Hamilton ベクトル場の差

$\mathcal L_X\omega=0$ を満たす $X$ を :term[シンプレクティックベクトル場]{en="symplectic vector field" yomi="しんぷれくてぃっくべくとるば"} という。Cartan の公式は $\mathcal L_X\omega=\mathrm d(\iota_X\omega)$ を与えるので、この条件は $\iota_X\omega$ の閉性に等しい。非退化性は $X$ と 1 形式の間の一対一対応を与える。従ってシンプレクティックベクトル場の実ベクトル空間を Hamilton ベクトル場で割ったものは、$X\mapsto[\iota_X\omega]$ によって $H_{\mathrm{dR}}^1(M)$ と同型である。実際、核は完全な内部積を持つ場であり、任意の閉 1 形式に対応する滑らかな場を逆行列で作れる。

:::example[トーラスの平行移動は局所的には Hamilton である]
周期 $1$ の $x,y$ による $T^2$ で $\omega=\mathrm dx\wedge\mathrm dy$ とする。$X=2\partial_x+\partial_y$ は平行移動を生成し、$\iota_X\omega=2\,\mathrm dy-\mathrm dx$ は閉じている。角度の枝が取れる近傍なら $H=2y-x$ が Hamilton 関数である。一方、$x$ 方向の円周積分は $-1$、$y$ 方向は $2$ である。従ってこの 1 形式は完全でなく、同じ場を作る大域的な $H$ は存在しない。
:::

:::example[可換な流れを作る二つの保存量]
$\R^4$ の標準座標で $H_j=(q_j^2+p_j^2)/2$（$j=1,2$）とする。$\{H_1,H_2\}=0$ なのは各項で変数の組が異なるためである。$H_1+2H_2$ の流れは第 1 平面で角速度 $-1$、第 2 平面で角速度 $-2$ の回転になる。例えば $(q_1,p_1,q_2,p_2)=(1,0,0,2)$ は時刻 $\pi/2$ で $(0,-1,0,-2)$ に移る。$H_1,H_2>0$ の共通等位集合は二つの円の積で、両方の半径が保たれる。括弧が零であることを知れば、二つの平面の運動を別々に求めて合成できる。
:::

## Darboux の定理

閉性と非退化性がありさえすれば、局所では座標を変えて係数を標準形にできる。点での線形変換はすでに :ref[thm-mfd-symplectic-basis] で作った。次は近傍で残る係数の違いを、時間依存の写像で取り除く。

:::theorem[Darboux の定理]{#thm-mfd-darboux}
任意のシンプレクティック多様体の各点に、形式が式 :eqref[eq-mfd-standard-symplectic] となる局所座標が存在する。
:::

:::proof
座標の線形変換で点を原点に移し、$\omega$ が原点で定数形式 $\omega_0$ と同じ値を取るようにする。$\omega_t=\omega_0+t(\omega-\omega_0)$（$0\le t\le1$）と置く。原点では全ての $t$ で同じ非退化形式である。行列式は連続であり、時間区間がコンパクトなので、共通の小さい座標球で全ての $\omega_t$ が非退化になる。各 $\omega_t$ は閉である。

星形の球上で :ref[thm-mfd-poincare-lemma] の原始形式を使い、$\mathrm d\beta=\omega-\omega_0$ とする。ホモトピー作用素の被積分関数は原点で零なので $\beta(0)=0$ である。$\iota_{Y_t}\omega_t=-\beta$ で時間依存の場 $Y_t$ を定める。逆行列の滑らかさにより $Y_t$ は両変数で滑らかで、原点では零である。原点の軌道は定数なので、その近くの初期点からの解は全区間 $[0,1]$ で定義される。これには時間を追加変数にした場 $\partial_t+Y_t$ の局所解を使える。$\omega_t$ は時間の端の近傍にも非退化な範囲で延びるので、定数軌道のコンパクトな時間区間を有限の局所解で覆い、初期点を絞ればよい。初期時刻零の解写像を $\psi_t$ とする。時刻を逆に解くことにより各 $\psi_t$ は像への微分同相である。

$\psi_{t+h}\circ\psi_t^{-1}$ は $h=0$ でベクトル場 $Y_t$ の一次変化を持つ。局所座標で引き戻しの係数とヤコビ行列を微分すると、これは時間 $t$ の自律場 $Y_t$ による Lie 微分と同じ一次項を与える。時間依存の形式自体の微分も加えると
$$
\frac{\mathrm d}{\mathrm dt}\psi_t^*\omega_t
=\psi_t^*(\partial_t\omega_t+\mathcal L_{Y_t}\omega_t)
=\psi_t^*(\mathrm d\beta-\mathrm d\beta)=0
$$
となる。最後に Cartan の公式と $\mathrm d\omega_t=0$ を使った。従って $\psi_1^*\omega=\omega_0$ である。$\psi_1$ を座標の逆写像として取れば局所標準形を得る。
:::

例えば平面で $\omega=e^x\mathrm dx\wedge\mathrm dy$ なら $q=e^x-1$、$p=y$ で直接標準形を得る。一般の場合の証明も、非退化性を用いて形式の変化を打ち消す場を求める手順である。正定値計量なら点での線形標準化は可能だが、近傍全体で定数にできるとは限らない。Darboux の議論では、形式が閉で、その差を局所的に外微分として書けることを使った。

### 大域で形式を変形できる条件

Darboux の証明で使った時間依存の計算は、初期点の共通近傍が取れる場合に限らない。コンパクトで境界のない $M$ 上に $0\le t\le1$ で滑らかに変わるシンプレクティック形式 $\omega_t$ があり、滑らかに変わる 1 形式 $\beta_t$ で $\partial_t\omega_t=\mathrm d\beta_t$ が成り立つとする。$\iota_{Y_t}\omega_t=-\beta_t$ を解くと、$[0,1]\times M$ で滑らかな場を得る。時間を追加変数とする自律系ではその軌道は有限な時間区間中、コンパクトな集合に留まる。そのため時刻の端まで解を延長でき、全ての初期点から時刻 $t$ の微分同相 $\psi_t$ が得られる。逆は終端から逆に積分する操作である。Cartan の公式を用いると
$$
\frac{\mathrm d}{\mathrm dt}\psi_t^*\omega_t
=\psi_t^*\bigl(\mathrm d\beta_t+\mathrm d(\iota_{Y_t}\omega_t)\bigr)=0
$$
となり、$\psi_t^*\omega_t=\omega_0$ である。これが Moser の方法の大域的な形である。

例えばコンパクトで連結かつ境界のない向き付けられた曲面上の二つの正の面積形式 $\omega_0,\omega_1$ が同じ面積を持つとする。差は最高次形式で積分が零なので、:ref[thm-mfd-top-cohomology] により $\omega_1-\omega_0=\mathrm d\beta$ である。直線補間 $\omega_t=(1-t)\omega_0+t\omega_1$ は、各点で同じ向きの正の形式の和だから、全時刻で非退化かつ閉である。従って $\beta_t=\beta$ として上の方法を適用でき、$\omega_1$ は $\omega_0$ と微分同相で対応する。逆に面積が違えば、向きを保つ微分同相の変数変換により対応できない。この曲面の議論では正の最高次形式を補間したが、高次元では二つの非退化 2 形式の直線補間が非退化とは限らない。例えば $\omega_1=-\omega_0$ なら途中で零になる。そのため高次元でこの方法を使う場合には非退化性を別に確認する。

### Moser の方程式を直接解く

:::example[形式の変化を一次元の変数変換として求める]
平面で $\omega=(1+x^2)\,\mathrm dx\wedge\mathrm dy$ とし、$\omega_t=(1+tx^2)\,\mathrm dx\wedge\mathrm dy$、$\beta=(x^3/3)\,\mathrm dy$ と取る。$\mathrm d\beta=x^2\,\mathrm dx\wedge\mathrm dy$ なので、$\iota_{Y_t}\omega_t=-\beta$ は
$$
Y_t=-\frac{x^3}{3(1+tx^2)}\partial_x
$$
を与える。解の $x$ 成分を $u(t)$ とすると
$$
\frac{\mathrm d}{\mathrm dt}\big(u+t u^3/3\big)
=(1+tu^2)u'+u^3/3=0.
$$
従って初期点 $(a,b)$ からの解は $u+t u^3/3=a,\ y=b$ で定まる。各 $t\in[0,1]$ で左辺は $u$ の狭義単調で全射な関数であり、導関数も正なので、解は全域で滑らかである。終時刻では $q=x+x^3/3,\ p=y$ が Darboux 座標であり、$\mathrm dq=(1+x^2)\mathrm dx$ から直接にも検証できる。例えば $x=1$ は $q=4/3$ に移る。
:::

:::example[横断変数を残して一つの座標だけ積分する]
$\{(x,y):|x|<1/2,\ |y|<1/2\}$ で $\omega=(2+x+y)\,\mathrm dx\wedge\mathrm dy$ とする。係数は正であり、$p=y$ を固定して $q=\int_0^x(2+s+y)\,\mathrm ds=(2+y)x+x^2/2$ と置く。すると $\mathrm dq=(2+x+y)\,\mathrm dx+x\,\mathrm dy$ なので $\mathrm dq\wedge\mathrm dp=\omega$ である。ヤコビ行列式は $2+x+y>0$ で、各 $y$ で $q$ が狭義単調増加だから、この座標は領域全体で像への微分同相になる。例えば $(x,y)=(1/4,0)$ は $(q,p)=(17/32,0)$ である。固定した $y$ に関する積分なので、$\mathrm dq$ には $y$ の微分の項が残るが、$\mathrm dp$ との外積で消える。
:::

## Lagrange 部分多様体

:::definition[等方・余等方・Lagrange 部分空間]{#def-mfd-lagrangian}
シンプレクティックベクトル空間の部分空間 $W$ が $W\subset W^\omega$ を満たすとき等方、$W^\omega\subset W$ を満たすとき余等方という。$W=W^\omega$ のとき :term[Lagrange 部分空間]{en="Lagrangian subspace" yomi="らぐらんじゅぶぶんくうかん"} という。部分多様体のすべての接空間が Lagrange 部分空間であるとき、その部分多様体を :term[Lagrange 部分多様体]{en="Lagrangian submanifold" yomi="らぐらんじゅぶぶんたようたい"} という。
:::

次元 $2n$ の空間では、等方なら $\dim W\le n$ である。実際 $\dim W^\omega=2n-\dim W$ だからである。従って等方で次元 $n$ なら両空間の次元が同じなので Lagrange である。逆に Lagrange なら等方で次元 $n$ である。部分多様体 $i:L\hookrightarrow M$ の場合は、$\dim L=n$ と $i^*\omega=0$ が必要十分条件になる。

等方な $W$ が半次元より小さければ、$W^\omega$ に属するが $W$ に属さないベクトルを一つ足しても等方性が保たれる。この操作を繰り返せば、任意の等方部分空間を Lagrange 部分空間へ延ばせる。従って Lagrange 部分空間は包含に関して極大な等方部分空間である。

余等方部分空間 $C$ では $\omega|_C$ の核が $C^\omega$ に等しい。そこで $C/C^\omega$ に $\overline\omega([v],[w])=\omega(v,w)$ と定めると、代表を変えても核の成分しか加わらないので値は変わらない。商での値が全ての $[w]$ で零なら $v\in C^\omega$ であり、$[v]=0$ だから $\overline\omega$ は非退化である。その次元は $2\dim C-\dim V$ となる。この商が、後に扱う特性方向に沿った簡約の各点での形である。

:::example[グラフと余接ファイバー]
余接束で零切断と各ファイバーは Lagrange 部分多様体である。零切断では $p=0$、ファイバーでは $\mathrm dq^i=0$ であり、どちらも $\lambda$ の制限が零で、次元は $n$ である。一般の 1 形式 $\alpha$ のグラフを与える切断 $s_\alpha:Q\to T^*Q$ では $s_\alpha^*\lambda=\alpha$ である。実際、$\pi\circ s_\alpha=\id_Q$ を標準 1 形式の定義に代入すればよい。従って $s_\alpha^*\omega=-\mathrm d\alpha$ であり、グラフが Lagrange となる必要十分条件は $\alpha$ が閉であることである。例えば $\alpha=y\,\mathrm dx+x\,\mathrm dy$ は $\mathrm d(xy)$ なので、$\R^2$ の余接束でそのグラフは Lagrange である。$\alpha=x\,\mathrm dy$ のグラフでは引き戻した形式が $-\mathrm dx\wedge\mathrm dy$ となり、Lagrange ではない。
:::

### 余法束とグラフの例

埋め込まれた部分多様体 $S\subset Q$ の余法束を $T^*Q$ の部分集合と見なすと、これは Lagrange 部分多様体である。局所スライス座標 $(u,z)$ で $S$ は $z=0$、余法の条件は $p_u=0$ である。従って余法束の像は $(u,p_z)$ を自由変数とするスライスであり、次元は $\dim S+(\dim Q-\dim S)=\dim Q$ である。標準 1 形式 $\lambda=p_u\,\mathrm du+p_z\,\mathrm dz$ はこのスライスに制限すると零になる。なぜなら $p_u=0$ であり、$\mathrm dz$ も接方向に対して零だからである。従って $-\mathrm d\lambda$ の制限も零で、半次元という条件と併せて Lagrange である。例えば $S=\{y=0\}\subset\R^2$ では、$T^*\R^2$ の条件は $y=0,p_x=0$ であり、$x$ と $p_y$ が自由である。

別の例として、同じ次元のシンプレクティック多様体 $M,N$ に対して積に形式 $-\pi_M^*\omega_M+\pi_N^*\omega_N$ を入れる。滑らかな $F:M\to N$ のグラフの引き戻しは $-\omega_M+F^*\omega_N$ である。グラフは積の半次元なので、Lagrange となることは $F^*\omega_N=\omega_M$ と同値である。$F$ が微分同相ならシンプレクティック同相の条件そのものである。例えば $\R^2$ の $(q,p)\mapsto(2q,p/2)$ のグラフは Lagrange だが、$(q,p)\mapsto(2q,p)$ のグラフは引き戻しが $-\omega+2\omega=\omega$ になるので Lagrange ではない。積の第 1 因子の負号が、二つの形式を比較する条件を作っている。

### Lagrange 部分多様体に沿う Darboux 座標

:::proposition[Lagrange 部分多様体に適合する局所標準形]{#prop-mfd-lagrangian-darboux}
$L\subset M$ を埋め込まれた Lagrange 部分多様体とする。$L$ の各点に、$\omega=\sum_i\mathrm dq^i\wedge\mathrm dp_i$ かつ $L=\{p=0\}$ となる局所座標が存在する。
:::

:::proof
まずスライス座標 $(q,z)$ で $L=\{z=0\}$ とし、$e_i=\partial_{q^i}$ を $L$ 上の基底とする。$L$ の各接空間が $L^\omega$ に等しいため、線形写像 $v\mapsto(\omega(e_i,v))_i$ の核は $T L$ である。従って $z$ 方向から $\R^n$ への制限は同型であり、逆行列の成分を使えば $L$ 上に滑らかな横断枠 $f_j$ が取れて $\omega(e_i,f_j)=\delta_{ij}$ となる。そこで $C_{ij}=\omega(f_i,f_j)$ と置き、
$$
\widehat f_j=f_j+\tfrac12\sum_i C_{ij}e_i
$$
と修正する。$\omega(e_i,e_j)=0$ なので $e_i$ との値は変わらない。二つの新しい横断ベクトルの値は $C_{ij}+\tfrac12 C_{ji}-\tfrac12 C_{ij}=0$ である。従って $(e_i,\widehat f_j)$ は各点でシンプレクティック基底になる。

元の座標の中で、$L$ 上の点 $(q,0)$ に $\sum_j p_j\widehat f_j(q)$ を加える写像を作る。$p=0$ での微分の列は $e_i,\widehat f_j$ なので、逆写像定理で点の近くの微分同相になる。それを使って $\omega$ を新しい座標に引き戻すと、$\widetilde\omega$ は零切断の点で、法方向を含めて標準形式 $\omega_0$ と一致する。小さい直積の座標域に絞り、$H(t,q,p)=(q,tp)$ を使う。閉形式 $\widetilde\omega-\omega_0$ の $H_0$ での引き戻しは零なので、ホモトピー公式によりこの差は $\mathrm d\beta$ と書ける。時間方向の速度は $p=0$ で零であり、従って $K$ で作った $\beta$ は零切断上で余接ベクトルとして零である。

$\omega_t=(1-t)\omega_0+t\widetilde\omega$ は、零切断上で全時刻に非退化なので、指定点の近くに共通の小さい領域を取りそこで非退化にできる。Darboux の証明と同じく $\iota_{Y_t}\omega_t=-\beta$ を解くと、$Y_t$ は零切断で消える。指定点の軌道は定数なので、近くの点には時刻 $1$ までの解写像 $\psi_t$ がある。Cartan の公式により $\psi_t^*\omega_t=\omega_0$ である。さらに零切断上の軌道は定数のため、$\psi_t$ はそれを点ごとに固定する。従ってこの写像を新しい座標の逆と合わせれば、形式が標準でありながら $L$ が $p=0$ に留まるような座標を得る。
:::

例えば $\R^4$ の $q=(x,y)$ と $p$ の標準座標で、$L$ が $p_x=y,p_y=x$ のグラフなら、$P_x=p_x-y,P_y=p_y-x$ と取るだけで $L$ を零切断に移せる。実際、新しい座標の形式は二つの交差項 $\mathrm dx\wedge(-\mathrm dy)$ と $\mathrm dy\wedge(-\mathrm dx)$ が相殺して元と一致する。この特別な例ではグラフが $\mathrm d(xy)$ なので、閉 1 形式の平行移動で調整できる。上の証明は最初の座標でまだグラフの形が簡単でない場合にも、零切断に沿う微分をそろえてから同じ目的を達成する。

### コンパクトな Lagrange 部分多様体の近傍

:::theorem[コンパクトな場合の Weinstein の Lagrange 近傍定理]{#thm-mfd-lagrangian-neighborhood}
$L\subset(M,\omega)$ をコンパクトな埋め込まれた Lagrange 部分多様体とする。$T^*L$ の零切断の開近傍から $M$ の $L$ の開近傍へのシンプレクティック同相で、各 $0_q$ を $q$ に写すものが存在する。
:::

:::proof
束計量を使って $TM|_L$ における $TL$ の滑らかな補束を取る。$\omega(v,s_q\alpha)=\alpha(v)$（$v\in T_qL$）となるよう、補束への写像 $s_q:T_q^*L\to T_qM$ を選ぶ。この写像は、$\omega|_{T_qL}$ が零であり $T_qL=(T_qL)^\omega$ であることから、局所標準形の証明と同様に可逆な行列を解いて得られる。$C(\alpha,\beta)=\omega(s_q\alpha,s_q\beta)$ とし、$\alpha(A_q\beta)=C(\alpha,\beta)/2$ となる $A_q:T_q^*L\to T_qL$ を取る。双対性により $A_q$ は一意かつ滑らかである。$\widehat s_q=s_q+A_q$ と置くと、$TL$ との対は変わらず、$\widehat s_q$ 同士の評価は $C(\alpha,\beta)-\alpha(A_q\beta)+\beta(A_q\alpha)=0$ になる。

:ref[thm-mfd-proper-weak-whitney] で $M$ を Euclid 空間の閉部分多様体と見て、管状近傍から $M$ への滑らかな後退 $R$ を取る。$q+\widehat s_q\alpha$ がその近傍に入る範囲で $\chi(q,\alpha)=R(q+\widehat s_q\alpha)$ と置く。$R|_M=\id$ なので、$0_q$ での微分は $(v,\alpha)\mapsto v+\widehat s_q\alpha$ である。この写像は同型で、上の評価により標準形式を $\omega_q$ に対応させる。従って $\chi$ は零切断の近くで局所微分同相である。ファイバーのノルムで一様に小さく絞れば単射にもできる。そうでなければ、同じ像を持つ異なる二点をノルムが零へ近づくように選べる。$L$ のコンパクト性により両者は部分列で零切断へ収束し、像が等しいから極限は同じ $0_q$ である。その点の局所単射性に反する。

こうして得た微分同相で $\omega$ を引き戻し、標準形式との差を $\gamma$ とする。$\gamma$ は零切断の各点で全接方向に対して零である。ファイバーを $t$ 倍するホモトピーと :ref[thm-mfd-homotopy-operator] を使えば $\gamma=\mathrm d\beta$ で、$\beta$ 自体が零切断上で零となる。標準形式と引き戻した形式の直線補間は、零切断上で両端が一致するため、コンパクト性により零切断の共通近傍で非退化である。そこで Darboux の証明と同じ Moser の方程式を解く。場は零切断で消えるので、零切断の各点を固定する。解の定義域の開性と零切断のコンパクト性により、時刻 $1$ まで流せる共通近傍を取れる。時刻 $1$ の写像に $\chi$ を合成すると、標準形式を引き戻す所望のシンプレクティック同相となる。
:::

### 母関数から Lagrange なグラフを作る

標準座標の同じ次元の二空間を $(q,p)$、$(Q,P)$ とし、積には $-\omega+\omega'$ を入れる。$(q,P)$ の滑らかな関数 $S(q,P)$ から
$$
p_i=\partial_{q^i}S,\qquad Q^i=\partial_{P_i}S
$$
と置くと、$(q,P)$ をパラメータとする部分多様体は Lagrange である。実際、$q,P$ は出力の成分に含まれるためパラメータ写像は埋め込みであり、次元は積の半分である。その上で $\sum_i(p_i\,\mathrm dq^i+Q^i\,\mathrm dP_i)=\mathrm dS$ だから、外微分は零である。左辺の外微分は $-\sum_i\mathrm dq^i\wedge\mathrm dp_i+\sum_i\mathrm dQ^i\wedge\mathrm dP_i$ であり、これが積の形式の引き戻しである。

さらに混合二階微分行列 $(\partial_{q^i}\partial_{P_j} S)_{ij}$ が可逆なら、$(q,P)\mapsto(q,p)$ は局所微分同相になる。ヤコビ行列の対角ブロックは恒等行列とこの混合行列だからである。従ってこの部分多様体は局所的に $(q,p)\mapsto(Q,P)$ のグラフであり、グラフが Lagrange なので $F^*\omega'=\omega$ である。非退化性から微分は単射で、同じ次元なので局所微分同相を得る。こうして正準変換を与える $S$ を :term[母関数]{en="generating function" yomi="ぼかんすう"} と呼ぶ。

:::example[せん断を母関数から求める]
一次元で $S(q,P)=qP+q^3/3$ とすると、$p=P+q^2$、$Q=q$ である。よって $(Q,P)=(q,p-q^2)$ であり、$(q,p)=(2,5)$ は $(2,1)$ に移る。形式は $\mathrm dQ\wedge\mathrm dP=\mathrm dq\wedge(\mathrm dp-2q\,\mathrm dq)=\omega$ となる。混合二階微分は $1$ なので解法は全域で使える。一方 $S=q^2+P^2$ だと $p=2q$ が固定され、$(q,p)$ を自由に動かせない。得られる積の部分多様体は Lagrange だが、始域全体の局所変換のグラフにはならない。
:::

## 運動量写像とシンプレクティック簡約の概観

対称性を持つ系では、対称性が作る方向を運動から取り除いて少ない変数で記述できる。シンプレクティック形式を保つだけでは商に形式が降りないため、まず運動量の値を固定して、消す方向を形式の核にする。

:::definition[運動量写像]{#def-mfd-momentum-map}
シンプレクティック形式を保つ $G$ の左作用を考える。滑らかな $J:M\to\mathfrak g^*$ が、全 $\xi\in\mathfrak g$ に対して
$$
\iota_{\xi_M}\omega=\mathrm d J^\xi,\qquad J^\xi(p)=\langle J(p),\xi\rangle
\label{eq-mfd-momentum-equation}
$$
を満たすとき、:term[運動量写像]{en="momentum map" yomi="うんどうりょうしゃぞう"} という。余随伴作用 $g\cdot\mu$ を $\langle g\cdot\mu,\xi\rangle=\langle\mu,\operatorname{Ad}_{g^{-1}}\xi\rangle$ で定める。以下の簡約では $J(g\cdot p)=g\cdot J(p)$ も仮定する。
:::

作用が形式を保つなら $\iota_{\xi_M}\omega$ は、Cartan の公式で外微分が零となる。ただし完全とは限らない。例えば $T^2$ の定数形式 $\omega=\mathrm dx\wedge\mathrm dy$ に対し第 1 成分の平行移動の基本場は $\partial_x$ で、内部積は $\mathrm dy$ である。その第 2 円周に沿う積分は $1$ なので原始関数を持たない。この作用に運動量写像は存在しない。

$H$ が作用で不変なら $J^\xi$ は $H$ の流れで保存される。実際 $X_HJ^\xi=\{J^\xi,H\}=-\xi_M H=0$ となる。これが対称性から保存量を得る仕組みである。等変な場合、$t\mapsto\exp(t\eta)$ による等変性を微分すると $\eta_MJ^\xi=J^{[\xi,\eta]}$ となる。これは $\operatorname{Ad}_{\exp(-t\eta)}$ の微分が $-\operatorname{ad}_\eta$ だからである。従って $\{J^\xi,J^\eta\}=J^{[\xi,\eta]}$ である。

### 余随伴軌道のシンプレクティック形式

$\mu\in\mathfrak g^*$ の余随伴軌道 $\mathcal O$ に、:ref[thm-mfd-homogeneous-orbit] の $G/G_\mu$ の構造を入れる。$\nu\in\mathcal O$ で $\xi$ が作る接ベクトルを $\xi_{\mathcal O}(\nu)$ と書き、
$$
\omega_\nu(\xi_{\mathcal O},\eta_{\mathcal O})=\langle\nu,[\xi,\eta]\rangle
$$
と置く。$\xi_{\mathcal O}(\nu)=0$ である条件は、余随伴作用を微分すると、全 $\eta$ に対し $\langle\nu,[\xi,\eta]\rangle=0$ である。従って式は接ベクトルを作る Lie 代数の元の選択によらない。また基本場が軌道の接空間を張るので、同じ条件から形式は非退化である。ある点の接空間の基底を作る基本場を選ぶと近くでも基底になり、その枠で表示した式の係数は滑らかである。

閉性を示すには、三つの基本場で $\mathrm d\omega$ を評価すればよい。例えば係数を微分する項 $\xi_{\mathcal O}\big(\omega(\eta_{\mathcal O},\zeta_{\mathcal O})\big)$ は $\langle\nu,[[\eta,\zeta],\xi]\rangle$ になる。その巡回和は Jacobi 恒等式で零である。括弧を含む項では $[\xi_{\mathcal O},\eta_{\mathcal O}]=-[\xi,\eta]_{\mathcal O}$ を入れると、同じ三重括弧の巡回和が現れ、やはり零になる。従ってこれは :term[軌道の標準シンプレクティック形式]{en="Kirillov–Kostant–Souriau form" yomi="きどうのひょうじゅんしんぷれくてぃっくけいしき"} である。作用で $\xi_{\mathcal O}$ は $(\operatorname{Ad}_g\xi)_{\mathcal O}$ に移る。$\operatorname{Ad}_g$ が括弧を保ち、余随伴作用がその逆で双対成分を移すので、形式の値も保たれる。包含 $\mathcal O\to\mathfrak g^*$ は等変であり、$\eta_{\mathcal O}\langle\nu,\xi\rangle=\langle\nu,[\xi,\eta]\rangle$ により運動量写像である。

例えば $\SO(3)$ では、$\xi$ とベクトルの対応 $\xi(v)=u\times v$ と Euclid 内積を使って $\mathfrak g^*$ を $\R^3$ と同一視できる。$g\in\SO(3)$ は内積と行列式を保つので、外積も $g(u\times v)=(gu)\times(gv)$ と変換する。従って随伴作用は $u\mapsto gu$、余随伴作用も同じ回転である。半径 $r>0$ の軌道は球面であり、北極 $\nu=re_3$ では $e_1,e_2$ の基本場が $-re_2,re_1$ となる。形式の値は $\nu\cdot(e_1\times e_2)=r$ であり、同じ二本に対する標準面積形式の値は $r^2$ である。両形式は回転で保存されるため、軌道の形式は球面全体で標準面積形式の $1/r$ 倍であり、全積分は $4\pi r$ になる。

:::theorem[正則な値における簡約]{#thm-mfd-symplectic-reduction}
等変運動量写像 $J$ の正則値 $\mu$ に対し、$Z=J^{-1}(\mu)$ が非空とする。余随伴作用で $\mu$ を固定する閉部分群 $G_\mu$ が $Z$ に自由かつ固有に作用するとき、商 $Z/G_\mu$ に一意なシンプレクティック形式 $\omega_{\mathrm{red}}$ があり、
$$
\pi^*\omega_{\mathrm{red}}=i^*\omega
\label{eq-mfd-reduced-form}
$$
を満たす。ただし $i:Z\hookrightarrow M$、$\pi:Z\to Z/G_\mu$ である。
:::

:::proof
:ref[thm-mfd-regular-value] により $Z$ は部分多様体で、$T_pZ=\Ker\mathrm dJ_p$ である。等変性により $G_\mu$ の作用は $Z$ を保ち、その商は :ref[thm-mfd-quotient-manifold] によって多様体となる。

$p\in Z$ で $W=\{\xi_M(p):\xi\in\mathfrak g\}$ とすると、:eqref[eq-mfd-momentum-equation] により $T_p Z=W^\omega$ である。従って $(T_pZ)^\omega=W$ であり、$i^*\omega$ の核は $T_pZ\cap W$ である。$\xi_M(p)$ が $T_pZ$ に入る条件は、等変性によって $\xi$ が余随伴軌道の $\mu$ で零の速度を与えることである。:ref[prop-mfd-orbit-rank] を余随伴作用に適用すると、これは $\xi\in T_eG_\mu$ と同値である。従ってこの核は $G_\mu$ の軌道方向に一致する。

$\overline v,\overline w\in T_{\pi(p)}(Z/G_\mu)$ を $T_pZ$ のベクトルに持ち上げて $i^*\omega$ で評価する。持ち上げの差は $\Ker\mathrm d\pi_p$ に入り、上で求めた核なので値は変わらない。同じ軌道の別の点に移しても、作用が $\omega$ を保つので値は変わらない。これで $\omega_{\mathrm{red}}$ が定まり、局所切断を使うと滑らかである。定義によって式 :eqref[eq-mfd-reduced-form] を満たす。全射沈め込みの引き戻しは局所切断で元の形式を復元できるので単射である。従って $\mathrm d i^*\omega=0$ から $\mathrm d\omega_{\mathrm{red}}=0$ が従い、一意性も従う。最後に、商の核に属するベクトルの持ち上げは $i^*\omega$ の核に入るため商で零である。従って商の形式は非退化である。
:::

:::example[複素空間の回転と射影空間]
$\C^m$ を実座標 $z_j=x_j+iy_j$ で表し、$\omega=\sum_j\mathrm dx_j\wedge\mathrm dy_j$ とする。$e^{it}$ が全ての成分を $e^{-it}$ 倍する作用を取る。基本場は $\sum_j(y_j\partial_{x_j}-x_j\partial_{y_j})$ なので、内部積は $\sum_j(x_j\mathrm dx_j+y_j\mathrm dy_j)$ である。従って $J=\sum_j|z_j|^2/2$ である。$c>0$ では $J=c$ は半径 $\sqrt{2c}$ の球面であり、原点を含まないため正則である。円周作用はこの球面上で自由で、コンパクト群の作用だから固有である。商は複素直線の空間 $\C P^{m-1}$ になる。実際、非零の各直線と球面の交わりが一つの軌道である。各比の座標域で分母に選んだ成分を正実数にする局所切断を使えば、商と複素射影座標が滑らかに一致する。

$m=2$ で比の座標 $w=a+ib$ を使うと、局所切断は $z=\sqrt{2c}(1,w)/\sqrt{1+|w|^2}$ である。$\omega=\mathrm d(\sum_j(x_j\,\mathrm dy_j-y_j\,\mathrm dx_j)/2)$ なので、切断で括弧内を引き戻すと $c(a\,\mathrm db-b\,\mathrm da)/(1+a^2+b^2)$ を得る。分母の微分も含めて外微分を取ると
$$
\omega_{\mathrm{red}}=\frac{2c}{(1+a^2+b^2)^2}\,\mathrm da\wedge\mathrm db
$$
となる。この座標の外に残るのは一点であり、全積分は極座標で $2\pi\int_0^\infty 2cr/(1+r^2)^2\,\mathrm dr=2\pi c$ である。$c=0$ では前像が原点だけで作用は自由でなく、上の定理の仮定を満たさない。
:::

### 余接持ち上げの運動量と零での簡約

$G$ が $Q$ に滑らかに作用するとする。各 $g$ の微分同相を余接束へ持ち上げると、合成についての逆写像と連鎖律から $T^*Q$ 上の左作用が得られる。先に示したように、この作用は標準 1 形式 $\lambda$ を保つ。$T^*Q$ 上の基本場を $\xi_{T^*Q}$ と書くと、$\mathcal L_{\xi_{T^*Q}}\lambda=0$ である。従って $\omega=-\mathrm d\lambda$ と Cartan の公式から
$$
\iota_{\xi_{T^*Q}}\omega
=\mathrm d(\iota_{\xi_{T^*Q}}\lambda)
$$
となる。標準 1 形式の定義と、持ち上げが射影と可換であることにより、
$$
J^\xi(\alpha_q)=\alpha_q\bigl((\xi_Q)_q\bigr)
\label{eq-mfd-cotangent-momentum}
$$
が得られる。これは $\xi$ に関して線形なので $\mathfrak g^*$ に値を取る写像 $J$ を定める。さらに $g$ の作用は基点を $gq$ に移し余接成分に $\mathrm d(a_{g^{-1}})$ を合成するので、
$$
J^\xi(g\cdot\alpha_q)
=\alpha_q\bigl((\operatorname{Ad}_{g^{-1}}\xi)_Q(q)\bigr)
$$
となる。右辺の場の変換は $g^{-1}\exp(t\xi)g$ を微分した結果である。これで余随伴作用に対する等変性が分かる。

作用が自由かつ固有なら、$B=Q/G$ は多様体であり、商射影 $\pi_Q:Q\to B$ は全射沈め込みである。零等位集合 $Z=J^{-1}(0)$ は、$\alpha_q$ が軌道方向を消すという条件を表す。零が正則値であることは、$J$ をファイバー方向で微分すれば分かる。各 $q$ で $\xi\mapsto(\xi_Q)_q$ は自由性により単射であり、その双対 $T_q^*Q\to\mathfrak g^*$ は全射である。この双対写像が $J$ のファイバー方向の微分である。$Z$ の元の余接成分は $\pi_Q$ の微分の核を消すため、唯一の $\beta_{\pi_Q(q)}\in T_{\pi_Q(q)}^*B$ で $\alpha_q=(\mathrm d\pi_Q)_q^*\beta$ と書ける。ここで全射な線形写像の双対の像が核の消滅空間であることを使った。

対応 $Z\to T^*B$、$\alpha_q\mapsto\beta_{\pi_Q(q)}$ は局所モデル $Q|_V\cong V\times G$ で $Z$ が $(b,g;\beta,0)$ と表されることから、滑らかな全射沈め込みである。持ち上げた作用は $g$ 成分だけを動かすので、各ファイバーは軌道であり作用は自由である。固有性は、$Z$ の作用の対写像のコンパクトな値集合の逆像において、基点の対のコンパクト集合と $Q$ 上の固有性が群の変数をコンパクト集合に制限することから従う。さらに出発点の $Z$ 成分も値の第 2 射影のコンパクト集合に入り、逆像はそれらの積の閉部分集合である。従って正則簡約の条件を満たす。商と $T^*B$ の滑らかな構造は全射沈め込みの普遍性により一致する。

標準 1 形式を比較する。$Z$ 内で速度 $W$ を取り、$\alpha_q=(\mathrm d\pi_Q)_q^*\beta$ とする。$\lambda_Q(W)$ はまず $W$ を $T_qQ$ へ射影し、それを $\alpha_q$ で評価した値である。$\beta$ へ降ろしてから $\lambda_B$ を評価すると、その同じベクトルをさらに $\mathrm d\pi_Q$ で $T_{\pi_Q(q)}B$ へ移し、$\beta$ で評価する。$\alpha_q$ の式により両方の値が等しい。従って $Z$ に制限した $\lambda_Q$ は $\lambda_B$ の引き戻しである。負の外微分を取ると、簡約の形式が $T^*B$ の標準形式に一致する。例えば $Q=\R^2$ で $\R$ が $x$ を平行移動すると、$J=p_x$、零等位集合では $p_x=0$ である。$x$ を商で除いた後に $(y,p_y)$ が残り、形式は $\mathrm dy\wedge\mathrm dp_y$ である。

### 不変 Hamilton 関数の力学を降ろす

正則簡約の仮定のもとで $H$ が $G$ 不変なら、$X_H$ は運動量を保存するので $Z$ に接する。$H$ の $Z$ への制限は各軌道で一定だから、ある滑らかな $\overline H$ によって $\pi^*\overline H=i^*H$ と書ける。$\overline v\in T_{\pi(p)}(Z/G_\mu)$ を $v\in T_pZ$ に持ち上げると
$$
\omega_{\mathrm{red}}(\mathrm d\pi_p X_H(p),\overline v)
=\omega(X_H(p),v)=\mathrm dH_p(v)
=\mathrm d\overline H_{\pi(p)}(\overline v)
$$
となる。従って $\mathrm d\pi_p X_H(p)$ は $\overline H$ の Hamilton ベクトル場に等しく、代表点 $p$ の選択に依存しない。このため簡約空間で解いた軌道は元の軌道の射影である。例えば上の平行移動作用で $H=(p_x^2+p_y^2)/2+y^2/2$ とすれば、零等位集合上では $\overline H=(p_y^2+y^2)/2$ となる。元の力学で $p_x=0$ のまま $x$ が定数であり、商では残る調和振動子だけが記述される。

### 特性方向に沿う商

シンプレクティック多様体 $(M,\omega)$ の埋め込まれた部分多様体 $i:S\hookrightarrow M$ で $(T_pS)^\omega\subset T_pS$ が全点で成り立つとする。このとき $i^*\omega$ の核は $(T_pS)^\omega$ である。次元は余次元と一致するので階数一定であり、閉形式の核を扱った [分布と Frobenius の定理](/math/manifolds/distributions/) の議論により積分可能である。これを特性分布と呼ぶ。各点の近くで、その葉をファイバーとする沈め込み $\pi:S\to Q$ がある。ここで $Q$ は局所的な横断座標域である。

$S$ 上の局所座標 $(u,z)$ で $u$ を特性方向とすると、$i^*\omega$ の外積基底展開には $\mathrm du$ が含まれない。なぜなら $\partial_u$ を入れると零になるため、$\mathrm du$ を含む係数が全て零でなければならないからである。さらに $\mathrm d(i^*\omega)=0$ なので各係数の $u$ 微分は零になる。$\mathrm du$ が一つだけ含まれる外微分の係数を見れば、他の項から同じ成分が出ないためである。従ってこの二形式は $z$ の座標だけの形式の引き戻しである。横断方向で退化するなら、それに対応する速度も核に入り特性方向に限られるので、降ろした形式は非退化である。これは正則な簡約の局所線形代数であり、正則簡約ではさらに、これらの局所の商を群の自由固有作用によって大域にする。

例えば標準の $\R^4$ で $S=\{p_1=0\}$ とすると、制限は $\mathrm dq^2\wedge\mathrm dp_2$ だけであり、核は $\partial_{q^1}$ が張る。局所または大域で $Q$ を $(q^2,p_2)$ 平面とすれば、制限した形式はその標準形式の引き戻しである。より一般に正則なエネルギー面 $S=H^{-1}(c)$ では $T_pS=\Ker\mathrm dH_p$ であり、$X_H$ はその直交補を張る。実際 $\omega(X_H,v)=\mathrm dH(v)=0$ が $T_pS$ 上で成り立ち、直交補が一次元だからである。また $X_H$ は $H$ を保つので $T_pS$ に入り、この面は余等方である。従ってエネルギー面内で Hamilton 軌道に沿う方向を除いたものが、局所的にシンプレクティックな商になる。大域で軌道空間が多様体になるには別の仮定が必要である。

### 正則簡約で変わる次元

正則値の等位集合の余次元は $\dim G$ であり、さらに自由な $G_\mu$ の軌道を商に取ると $\dim G_\mu$ だけ減る。従って簡約空間の次元は
$$
\dim M-\dim G-\dim G_\mu
$$
である。$\mu=0$ では余随伴作用で常に固定されるので $G_\mu=G$ となり、群の次元の二倍が減る。これは一度目に運動量を固定し、二度目に同じ値を持つ軌道の点を同一視するためである。例えば $\C^3$ のスカラー回転を正の値で簡約すると、$6-1-1=4$ 次元の $\C P^2$ になる。非零の値に対して $G_\mu$ が $G$ の真の部分群なら、第二段階で消えるのは $\dim G_\mu$ だけである。また $G_\mu$ の次元を求めるには、余随伴軌道写像の $\mu$ での速度が零となる Lie 代数の元を数えればよい。

また運動量写像の条件は、単に値が軌道上で一定であるという条件ではない。余随伴軌道が非自明なときには、同じ $G$ 軌道に沿って運動量の値も変わる。そのため一つの値 $\mu$ の等位集合だけを扱う際には、その値を固定する $G_\mu$ の作用を使う。

### 各基本場に原始関数があっても等変性は別に要る

:::example[平面の並進の運動量写像]
$(q,p)$ 平面に加法群 $\R^2$ を並進で作用させる。二つの基本場は $\partial_q,\partial_p$ で、その内部積は $\mathrm dp,-\mathrm dq$ である。従って $J=(p,-q)$ は運動量写像である。ところが群は可換なので余随伴作用は恒等であり、等変性には $J$ が並進で不変なことが必要である。実際には $(a,b)$ の並進で $J$ は $(b,-a)$ だけ変わる。原始関数に定数を足してもこの差は消えない。また $\{p,-q\}=1$ だが二方向の Lie 括弧は零である。等変な場合の Poisson 括弧の式とは、この定数の分だけ異なる。
:::

一般の運動量写像でも $\{J^\xi,J^\eta\}-J^{[\xi,\eta]}$ の微分は零である。実際、:ref[prop-mfd-poisson-laws] と左作用の基本場の符号を合わせれば、二つの関数の Hamilton ベクトル場は一致する。$\omega$ との内部積を取ると二つの微分が等しい。従って連結な $M$ では差は定数である。等変性がこの定数を零にすることは、余随伴の微分からすでに確認した。上の並進では可換なため $J^{[\xi,\eta]}=0$ であり、定数の変更では差を除けない。

:::example[並進一方向だけを簡約する]
$(q_1,p_1,q_2,p_2)$ の標準空間で $\R$ が $q_1$ に平行移動する。運動量は $J=p_1$ であり、全点で $\mathrm dJ$ は全射である。等位集合 $p_1=c$ では $\mathrm dp_1$ の制限が零なので、制限した形式は $\mathrm dq_2\wedge\mathrm dp_2$ になる。その核は $\partial_{q_1}$ の張る直線で、ちょうど軌道方向である。作用は自由かつ固有である。実際、作用前後の点がコンパクト集合に入ると、第 1 成分の差が移動量を有界に制御する。商は $(q_2,p_2)$ の平面であり、制限した形式がその標準形式に降りる。例えば $H=p_1^2/2+(q_2^2+p_2^2)/2$ は商で $c^2/2+(q_2^2+p_2^2)/2$ となり、定数項は商の運動に影響しない。
:::

## 演習問題

::::exercise{level=basic}
$\omega=2\mathrm dq\wedge\mathrm dp$ とする。$e=\partial_q$ を固定したとき、$\omega(e,f)=1$ となる $p$ 方向の $f$ を求めよ。
:::solution
$f=a\partial_p$ と置けば $\omega(e,f)=2a$ である。従って $f=(1/2)\partial_p$ である。$e,f$ は線形独立であり、形式は交代的なのでこれでシンプレクティック基底になる。
:::
::::

::::exercise{level=basic}
標準形式 $\mathrm dq\wedge\mathrm dp$ と $H=p^2/2+q^4$ に対し Hamilton 方程式を求めよ。
:::solution
$\partial_p H=p$、$\partial_q H=4q^3$ であり、式 :eqref[eq-mfd-hamilton-coordinates] から $\dot q=p$、$\dot p=-4q^3$ である。これを $H$ の全微分に代入すると $4q^3p-4pq^3=0$ で、$H$ は保存される。
:::
::::

::::exercise{level=basic}
$\R^2$ の標準形式について $\{q^2,p^3\}$ を求め、$(q,p)=(2,1)$ で評価せよ。
:::solution
$\partial_q(q^2)=2q$、$\partial_p(p^3)=3p^2$ で、他の偏微分は零である。従って括弧は $6qp^2$ となり、指定点で $12$ である。
:::
::::

::::exercise{level=basic}
$T^*\R^2$ の $\alpha=x\,\mathrm dx+2y\,\mathrm dy$ のグラフは Lagrange か判定せよ。
:::solution
$\alpha=\mathrm d(x^2/2+y^2)$ なので閉である。そのグラフへの $\omega$ の引き戻しは $-\mathrm d\alpha=0$ である。グラフは二次元で周囲は四次元なので Lagrange である。
:::
::::

::::exercise{level=advanced}
$\omega=(1+x^2)\mathrm dx\wedge\mathrm dy$ を標準形にする大域座標を求めよ。
:::solution
$q=x+x^3/3$、$p=y$ と置けば $\mathrm dq=(1+x^2)\mathrm dx$ なので $\mathrm dq\wedge\mathrm dp=\omega$ である。$q$ は微分が正で、$x\to\pm\infty$ で両方向に発散する。従って $\R$ 上の全単射であり、局所逆は逆写像定理で滑らかである。よって $(q,p)$ は大域座標になる。
:::
::::

::::exercise{level=advanced}
$T^*\R^2$ 上の $\omega$ に対し、$F(q,p)=(q,p+\nabla S(q))$ がシンプレクティック同相であることを示せ。
:::solution
逆は $F^{-1}(q,p)=(q,p-\nabla S(q))$ なので微分同相である。標準 1 形式は $F^*\lambda=\lambda+\sum_i\partial_iS\,\mathrm dq^i=\lambda+\mathrm dS$ と変わる。負の外微分を取ると $-\mathrm d^2S=0$ なので $F^*\omega=\omega$ である。
:::
::::

::::exercise{level=advanced}
$\omega=\mathrm dq\wedge\mathrm dp$ で $\R$ が $q$ の平行移動で作用する。運動量写像を求め、各値での簡約空間を特定せよ。
:::solution
作用の基本場は $\partial_q$ で、内部積は $\mathrm dp$ なので $J=p$ とできる。各値 $c$ は正則で、等位集合は直線 $p=c$ である。平行移動作用は自由であり、始点と終点の差が移動量を決めるため固有である。一つの軌道が等位集合全体をなすので簡約空間は一点である。その接空間は零であり、零の 2 形式は零次元で非退化である。
:::
::::

::::exercise{level=advanced}
シンプレクティック多様体で、$f,g$ の Hamilton ベクトル場が可換であることと、$\{f,g\}$ が各連結成分で定数であることの同値性を示せ。
:::solution
:ref[prop-mfd-poisson-laws] から $[X_f,X_g]=0$ は $X_{\{f,g\}}=0$ と同値である。定義と非退化性によりこれは $\mathrm d\{f,g\}=0$ と同値である。微分が零の関数は各連結成分で定数であり、逆に成分上定数なら微分は零である。例えば $f=q$、$g=p$ では括弧が定数 $1$ で、二つの場は定数場なので可換である。
:::
::::

::::exercise{level=advanced}
$(\R^4,\omega=\mathrm dq^1\wedge\mathrm dp_1+\mathrm dq^2\wedge\mathrm dp_2)$ で $H=(p_1^2+(q^1)^2+p_2^2+(q^2)^2)/2$ とする。$c>0$ の等位面上の Hamilton 流を求め、特性方向の商がどの例になるか説明せよ。
:::solution
$\dot q^i=p_i$、$\dot p_i=-q^i$ である。複素成分 $z_i=q^i+ip_i$ にまとめると $z_i'=p_i-iq^i=-iz_i$ なので $z_i(t)=e^{-it}z_i(0)$ である。$c>0$ では少なくとも一つの成分が非零だから、全成分が元に戻る時刻は $2\pi$ の整数倍に限られる。従って流れはこの面上で自由な円周作用を与える。球面は半径 $\sqrt{2c}$ で、作用がコンパクト群のものなので固有である。

この作用の運動量が $H$ であるのは、基本場が $X_H$ で、$\iota_{X_H}\omega=\mathrm dH$ を満たすからである。よって正則簡約の例で求めた $\C P^1$ が商である。特性方向は $X_H$ の張る一次元の方向であり、この場合は大域で各葉が一つの円周軌道になる。比の座標 $w=a+ib$ で商の形式の係数は $2c/(1+a^2+b^2)^2$ なので、例えば $w=0$ で $2c$、$w=1$ で $c/2$ となる。最高次の全積分は $2\pi c$ であり、$c$ の変更で局所の商の次元は変わらないが面積は変わる。
:::
::::

## 確認テスト

::::quiz{answer="4" tol="1e-6"}
$T^*\R^4$ の Lagrange 部分多様体の次元はいくつか。
:::explanation
余接束の底とファイバーはともに $4$ 次元なので全体は $8$ 次元である。シンプレクティック直交補の次元公式により $W=W^\omega$ なら次元はその半分の $4$ である。
:::
::::

::::quiz
シンプレクティック形式の条件を選べ。

- [x] 閉 2 形式である。
- [x] 各点で非退化である。
- [ ] 全域で完全形式である。

:::explanation
最初の二条件が定義である。閉性は外微分が零、非退化性はベクトルからその内部積への写像が線形同型であることを意味する。コンパクトな境界なしの正次元シンプレクティック多様体では、完全であれば最高次の外積の積分が Stokes の定理で零になって矛盾する。
:::
::::

::::quiz{answer="6" tol="1e-6"}
標準形式で $f=q^2$、$g=p$ とする。点 $(3,0)$ での $\{f,g\}$ を求めよ。
:::explanation
局所公式の第 1 項は $(2q)\cdot1$、第 2 項は零である。従って $\{f,g\}=2q$ となり、指定点では $6$ である。
:::
::::

::::quiz
正則な運動量の等位集合 $Z$ で、$i^*\omega$ が商に降りるために使った性質を選べ。

- [x] 作用が形式を保つ。
- [x] 商射影の鉛直方向が $i^*\omega$ の核に入る。
- [ ] $i^*\omega$ が $Z$ 自体で非退化である。

:::explanation
同じ商の点の代表を軌道上で取り替えても値を一致させるため、作用による不変性を使う。商の接ベクトルの持ち上げを替えると鉛直成分だけが変わるため、その成分が核に入ることも必要である。正則簡約では核がちょうど鉛直方向に一致し、その方向が商で消えるため、降りた形式が非退化になる。
:::
::::

---
title: "非可換ゲージ理論"
description: "接続と曲率から Yang–Mills 作用、Wilson 線、物質の結合と自己相互作用を導く。"
notation:
  - symbol: 'T^a,f^{abc}'
    meaning: 生成子と構造定数
  - symbol: 'F_{\mu\nu}'
    meaning: ゲージ接続の曲率
  - symbol: 'W_\gamma'
    meaning: 経路 γ に沿う平行移動
---

## この章の目標

- 局所変換を保つ共変微分を作る。
- 曲率の共変性と Bianchi 恒等式を証明する。
- 表現の群不変量と Yang–Mills の三点・四点頂点を計算する。

自然単位系 $\hbar=c=1$、計量 $\eta=\diag(+1,-1,-1,-1)$ を用いる。特に断らない限りゲージ代数の成分は Hermite 生成子で表し、$D_\mu=\partial_\mu-ig A_\mu$ の符号を用いる。

## ゲージ不変性と共変微分

### 異なる点で内部ベクトルを比べる

行列群 $G$ のユニタリ表現で、物質場 $\psi$ を $\psi'(x)=U(x)\psi(x)$ に移す。通常微分には $(\partial_\mu U)\psi$ が加わるので、$D_\mu=\partial_\mu-igA_\mu$ に接続 $A_\mu$ を導入する。ここでは $g>0$ とし、Abel 極限では前の QED の符号つき電荷 $q=-g$ に対応する。$D_\mu'\psi'=UD_\mu\psi$ を任意の $\psi$ に課すと
$$
A_\mu'=UA_\mu U^{-1}-\frac{i}{g}(\partial_\mu U)U^{-1}
\label{eq-qft-ym-gauge}
$$
である。$U=1+ig\omega+\mathcal O(\omega^2)$ なら $\delta A_\mu=\partial_\mu\omega-ig[A_\mu,\omega]$。随伴場 $X'=UXU^{-1}$ の共変微分は $D_\mu X=\partial_\mu X-ig[A_\mu,X]$ であり、積の微分で $\partial U$ の二項が相殺する。

### 有限変換の合成と純ゲージ

:eqref[eq-qft-ym-gauge] で $U$ の後に $V$ を作用させると、微分項は $(\partial V)V^{-1}+V(\partial U)U^{-1}V^{-1}=(\partial(VU))(VU)^{-1}$ となる。従って合成は $VU$ であり、順序は一般に交換できない。$A_\mu=0$ を $U$ で移した場 $A_\mu=-i(\partial_\mu U)U^{-1}/g$ では $D_\mu=U\partial_\mu U^{-1}$ と作用する。$D_\mu D_\nu$ を合成すると中間の $U^{-1}U$ が消え、微分同士が可換だから曲率は零である。例えば $U=e^{ig\omega(x)T^3}$ なら $A_\mu=(\partial_\mu\omega)T^3$。一方、異なる色を回す二つの指数の積では、最初の微分項が後の変換で共役されるため単純な和にはならない。

### 有限変換を逆向きに確認する

任意の列ベクトル関数 $f$ に $U\partial_\mu U^{-1}$ を作用させると、$\partial_\mu f+U(\partial_\mu U^{-1})f$。$U U^{-1}=1$ の微分により $U(\partial_\mu U^{-1})=-(\partial_\mu U)U^{-1}$ である。従って $D_\mu'=U D_\mu U^{-1}$ から :eqref[eq-qft-ym-gauge] を得る。特に $U$ がユニタリなら $(\partial_\mu U)U^{-1}$ は反 Hermite なので、$-i/g$ を掛けた項は Hermite。$SU(N)$ では $\partial_\mu\det U=0$ からこの項のトレースも零になる。よって変換後の $A_\mu$ は同じ Lie 代数に留まる。

:::example[一定場を点ごとに回す]
$SU(2)$ で $A_2=aT^1$、他成分は零とし、$U(x)=e^{ig\omega(x^1)T^3}$ を作用させる。$\zeta=g\omega$ とすると $UT^1U^{-1}=T^1\cos\zeta-T^2\sin\zeta$。実際、$\zeta$ による微分は $i[T^3,\,\cdot\,]$ であり、その初期値と初期微分は $T^1,-T^2$。従って $A_2'=a(T^1\cos\zeta-T^2\sin\zeta)$、$A_1'=\omega' T^3$。$A_2'$ の普通微分は一般に零でないが、$D'_1 A'_2=\partial_1 A'_2-ig[A'_1,A'_2]=0$。この計算では接続の微分項を落とすと純ゲージな背景に見かけの曲率を作ってしまう。
:::

### 固定した点から半径方向の接続を消す

星形の座標領域で基準点を $x_0$ とし、各点 $x$ への直線上の移送行列を $W(x,x_0)$ とする。これは直線に沿う一次の方程式 $\partial_s W=ig(x-x_0)^\mu A_\mu(x_0+s(x-x_0))W$ の $s=1$ の値である。$U(x)=W(x,x_0)^{-1}$ と選ぶと $(x-x_0)^\mu A'_\mu(x)=0$。実際、直線に沿って $\partial_s U=-igU(x-x_0)^\mu A_\mu$ なので、:eqref[eq-qft-ym-gauge] の共役項と接線の微分項が相殺する。この条件を :term[動径ゲージ]{en="radial gauge" yomi="どうけいげーじ"} と呼ぶ。基準点で $U(x_0)=1$ とし、$A$ が滑らかなら同じ領域で局所的に構成できる。

このゲージで $y=x-x_0$ とする。条件 $y^\mu A'_\mu=0$ を $\nu$ で微分すると $y^\mu\partial_\nu A'_\mu=-A'_\nu$。従って $y^\mu F'_{\mu\nu}=y^\mu\partial_\mu A'_\nu+A'_\nu$ で、交換子も消える。直線上でこれは $\partial_s[s A'_\nu(x_0+sy)]=s y^\mu F'_{\mu\nu}(x_0+sy)$。従って
$$
A'_\nu(x_0+y)=\int_0^1 s\,\dd s\,y^\mu F'_{\mu\nu}(x_0+sy).
$$
下端では滑らかさにより $s A'_\nu\to0$。特に $F=0$ なら $A'=0$ になり、平坦な接続はこの領域では純ゲージである。閉じた空間や基準点から全ての点へ一意に直線を引けない領域では、この証明を大域的に適用できない。他の経路で移送した行列の情報は、後の Wilson ループで見る。

### 局所座標の接続を重なりで比較する

二つの開領域 $O_i,O_j$ の共通部分で、同じ物質場を $\psi_j=U_{ji}\psi_i$ と表す。$U_{ji}$ は点に依存するユニタリな遷移行列である。$D_j\psi_j=U_{ji}D_i\psi_i$ が成り立つなら、接続 $A_i,A_j$ の関係は :eqref[eq-qft-ym-gauge] と同じである。三つの領域の共通部分では、$\psi_i$ から $\psi_k$ へ直接変える場合と $j$ を経由する場合が同じなので $U_{ki}=U_{kj}U_{ji}$。各領域の場をさらに $V_i$ で変換すると $U_{ji}\mapsto V_jU_{ji}V_i^{-1}$ となる。従って曲率や共変微分を局所式で計算しても、重なりでは同じ物理量の成分になっている。例えば $\tr F_{\mu\nu}F^{\mu\nu}$ は共役で不変なので、領域の取り方によらない密度を定める。

経路が領域をまたぐ場合は、経路の断片ごとに移送し、重なりで $U_{ji}$ を挿入する。区間ごとの変換則を積に適用すると、経由点の $V_j^{-1}$ と $V_j$ が順に相殺される。最後に残るのは始点・終点の基底変換だけである。例えば同じ点で領域 $i$ から $j$ へ移り、すぐ $i$ へ戻れば $U_{ij}U_{ji}=1$。重なりで起こる成分の変更を、時空で実際に曲率に沿って移送した結果と混同しない。この局所記述を使うと、単一の行列関数 $A_\mu(x)$ で空間全体を覆う仮定を置かずにゲージ場を扱える。

:::example[二次元トーラス上の平坦な接続]
二つの周期方向の長さを $L_1,L_2$ とし、$A_i=a_i T^3$ の一定 $SU(2)$ 場を考える。接続は可換なので曲率は零。それでも各基本閉路の移送は $W_i=e^{igL_i a_iT^3}$ であり、周期な変換で消せるのは両方が $1$ となる場合だけである。例えば $gL_1a_1=2\pi$、$gL_2a_2=0$ なら第一の移送は中心の $-1$、第二は $1$。基本表現の規格化トレースは第一の閉路で $-1$、随伴のトレースは $1$ になる。局所の方程式が同じでも、どの表現を測定に使えるかで閉路が区別する情報が変わる。
:::

### 微小変換の閉じた交換子を確かめる

$\delta_\omega A_\mu=D_\mu\omega$ とし、二つのパラメータを場から独立に固定する。場の方向微分として $\delta_1(\delta_2 A_\mu)=-ig[D_\mu\omega_1,\omega_2]$。従って
$$
(\delta_1\delta_2-\delta_2\delta_1)A_\mu
=-ig D_\mu[\omega_1,\omega_2].
$$
共変微分が交換子の積の法則を満たすためである。ここで $\delta_1\delta_2$ は、$\delta_2 A$ の関数を $\delta_1 A$ の方向へ変化させるという順を採った。群の行列そのものの積で先にどちらを作用させるかという規約とは、符号を比べる前に作用の順を揃える。例えば $\omega_1,\omega_2$ が共に一定で一つの Cartan 生成子に比例すれば交換子は零。異なる非可換方向では、微小であっても閉じた変換は第三の色方向を生む。この性質が三つのゲージ脚の色係数を決める。

## Lie 代数と表現の基本事項

### 群因子の規格化

:ref[prop-lie-classical-algebras] に現れる Lie 代数のうち、ここでは $SU(N)$ を基準に Hermite 生成子 $T^a$ を取り、$[T^a,T^b]=if^{abc}T^c$ とする。基本表現で $\tr(T^aT^b)=\delta^{ab}/2$。$\tr([T^a,T^b]T^c)=\tr(T^a[T^b,T^c])$ なので $f^{abc}$ は全反対称である。Jacobi 恒等式は三重積を展開して各項が逆符号で二度出ることから従う。

:term[Dynkin 指数]{en="Dynkin index" yomi="でぃんきんしすう"} $T(R)$ と二次 Casimir $C_2(R)$ を $\tr_R(T^aT^b)=T(R)\delta^{ab}$、$\sum_a T^aT^a=C_2(R)1_R$ と定める。後者は全生成子と可換であり、実際、交換子の二項は $f$ の反対称性で相殺する。既約表現では Schur の補題を適用できる。トレースを取って $C_2(R)\dim R=T(R)\dim G$。基本表現なら $C_F=(N^2-1)/(2N)$ である。

完全性は正規直交基底に $1/\sqrt{2N}$ を加えた行列基底から
$$
\sum_a T^a_{ij}T^a_{kl}=\tfrac12(\delta_{il}\delta_{jk}-\delta_{ij}\delta_{kl}/N)
\label{eq-qft-su-completeness}
$$
となる。従って $\sum_a T^a X T^a=\frac12(\tr X\,1-X/N)$。二重交換子を展開して $X$ を無跡とすれば $\sum_a[T^a,[T^a,X]]=2C_F X+X/N=NX$。随伴 Casimir は $C_A=N$ である。例えば $SU(3)$ では $C_F=4/3$、$C_A=3$。

### 交換子と反交換子の色構造

基本表現の任意の積を単位行列と無跡成分に分けると
$$
T^aT^b=\frac{\delta^{ab}}{2N}1+\frac12(d^{abc}+if^{abc})T^c,
\qquad d^{abc}=2\tr(\{T^a,T^b\}T^c).
$$
最初の係数は両辺のトレース、後の係数は $T^c$ を掛けたトレースで決まる。$d$ は対称、$f$ は反対称である。例えば $SU(2)$ の $T^a=\sigma_a/2$ では $d=0$。$SU(3)$ の $T^8=\diag(1,1,-2)/(2\sqrt3)$ では $\tr[(T^8)^3]=-1/(4\sqrt3)$、従って $d^{888}=-1/\sqrt3$。この対称トレースは [アノマリー](/physics/quantum-field-theory/anomalies/) で使う。

### 表現を取り替えたときの不変量

単純なコンパクト Lie 代数を考え、随伴を $(T_A^a)_{bc}=-if^{abc}$ で表す。実際、基底 $T^c$ に $[T^a,\,\cdot\,]$ を作用させると $[T^a,T^c]=if^{acb}T^b=-if^{abc}T^b$ である。従って $\tr(T_A^a T_A^b)=f^{acd}f^{bcd}=C_A\delta^{ab}$。ここで $T(A)=C_A$ となるのは、Casimir と指数の関係で $\dim A=\dim G$ のためである。一般の既約表現の $\tr_R(T^aT^b)$ が基準の $\delta^{ab}$ に比例することは、トレースの巡回性を使ってこの双線形形式が随伴作用で不変であることを示し、基準内積で一方の添字を上げて Schur の補題を適用すれば分かる。直和を取ると $T(R)$ は加法的。Casimir は直和の各既約成分に対応する値を持ち、全体で一つの値を持つとは限らない。

### テンソル積の指数と二次 Casimir

テンソル積 $R\otimes S$ の生成子は $T_R^a\otimes1+1\otimes T_S^a$。これは $e^{itT_R^a}\otimes e^{itT_S^a}$ を $t=0$ で微分して得る。どの $T_R^a$ も無跡である。単純 Lie 代数は交換子で張られ、表示した交換子のトレースが零だからである。従って指数の積の交差項のトレースが消え、
$$
T(R\otimes S)=\dim S\,T(R)+\dim R\,T(S).
$$
例えば $SU(N)$ で $F\otimes\bar F$ の指数は $N$。行列空間に対応づけると、$U$ の作用は $X\mapsto UXU^{-1}$ であり、単位行列と無跡空間に分かれる。単位行列は一重項で指数零なので、無跡空間、すなわち随伴表現の指数は $N$ である。完全性から求めた $C_A$ と一致する。

### 共役表現と中心元を使って色を分類する

共役表現 $R^*$ では $U$ は $U^*$ に移るため、$U=e^{i\alpha^a T_R^a}$ の微分で $T_{R^*}^a=-(T_R^a)^{\mathsf T}$。指数と Casimir は符号が二回出るので $R$ と等しい。例えば二つの基本表現の反対称部分 $A_{ij}=-A_{ji}$ を $SU(3)$ で $v_k=\varepsilon_{kij}A_{ij}/2$ に写すと、$\det U=1$ から $v'_k=(U^{-1})_{jk}v_j$。従ってこの三成分は共役基本である。一方同じテンソル積の対称部分には六成分がある。次元は独立な対角 $3$ と異なる添字の組 $3$ の和である。

$SU(N)$ の中心元 $z=e^{2\pi i/N}1$ は、基本表現では $z$、共役では $z^{-1}$、随伴では一として作用する。中心元の作用が一である表現だけが、中心で割った群 $SU(N)/\mathbb Z_N$ の表現になる。実際、同じ剰余類の二つの元は中心元だけ異なるため、行列 $U_R$ が剰余類で一意になる必要がある。従って随伴だけを局所的に扱う理論と、基本表現の外部源を許す理論では、Lie 代数と局所の頂点が同じでも選べる大域的なゲージ群が違い得る。例えば $SU(2)$ の中心の $-1$ は基本二成分に負号、随伴の三成分には単位として作用する。

:::example[群の生成子の規格化を別の表現で検査する]
$SU(2)$ で角運動量 $j$ の既約表現を取る。$T^3$ の固有値は $m_j=-j,\ldots,j$ なので $\tr[(T^3)^2]=j(j+1)(2j+1)/3$。等差の二乗の和でこの式を得る。従って指数は $T(j)=j(j+1)(2j+1)/3$。Casimir は $3T(j)/(2j+1)=j(j+1)$ である。$j=1/2$ で指数 $1/2$、$j=1$ で指数 $2$。異なる表現を使って同じ曲率の作用をトレース表示する場合、$-\tr_R F^2/[4T(R)]$ とすれば成分表示の $-F^a_{\mu\nu}F^{a\mu\nu}/4$ を再現する。
:::

### $SU(N)$ の完全性を成分で導く

複素行列空間の Hilbert–Schmidt 内積を $(X,Y)=\tr(X^\dagger Y)$ とする。無跡 Hermite 行列を複素線形に拡張した空間は無跡行列空間になる。基本表現の $T^a$ に単位行列 $T^0=1/\sqrt{2N}$ を加えれば、$\sqrt2 T^A$ は全行列空間の正規直交基底。任意の $X$ に対し $X=2\sum_A T^A\tr(T^A X)$ と展開できる。$X=E_{lk}$ を取ると $(i,j)$ 成分は $\delta_{il}\delta_{jk}=2\sum_A T^A_{ij}T^A_{kl}$。$A=0$ の項 $\delta_{ij}\delta_{kl}/N$ を引けば :eqref[eq-qft-su-completeness] である。ここでは Hermite 行列の添字を数字で書いたため、複素共役が $\tr(T^A X)$ の中で添字の交換として現れる。この交換の順を替えると恒等写像ではなく転置の写像になってしまう。

完全性から構造定数の二つの縮約も得られる。例えば $\sum_a T^a T^b T^a=-T^b/(2N)$。従って $\sum_a[T^a,[T^a,T^b]]=N T^b$。これを成分で書けば $f^{acd}f^{bcd}=N\delta^{ab}$、全ての色を足した二乗は $N(N^2-1)$。さらに対称係数 $d^{abc}$ について
$$
\sum_a\{T^a,\{T^a,T^b\}\}
=2C_F T^b+2\sum_a T^aT^bT^a
=\frac{N^2-2}{N}T^b
$$
となる。一方 $\{T^a,T^b\}=\delta^{ab}1/N+d^{abc}T^c$ を二回使うと、単位行列の成分は $\sum_a d^{a b a}/N=0$。実際 $\sum_a d^{aba}=4\tr(T^b\sum_a T^aT^a)=0$。$T^b$ の項に $2/N$ が入り、残る縮約は
$$
d^{acd}d^{bcd}=\frac{N^2-4}{N}\delta^{ab}.
$$
例えば $N=2$ ではこの係数が零で、Pauli 行列の反交換子が単位行列になることと一致する。$N=3$ では $f$ の二乗の全色和が $24$、$d$ の方は $40/3$ である。

## Yang–Mills Lagrangian と場の強さ

二つの微分を異なる順に行った差は曲率になる。$[D_\mu,D_\nu]$ を展開すると $\partial$ が物質場に作用する交差項が相殺し、
$$
F_{\mu\nu}=\frac{i}{g}[D_\mu,D_\nu]
=\partial_\mu A_\nu-\partial_\nu A_\mu-ig[A_\mu,A_\nu]
\label{eq-qft-ym-curvature}
$$
である。従って $F'=UFU^{-1}$。$[D_\mu,[D_\nu,D_\rho]]+$ 巡回和が Jacobi 恒等式で零になるため $D_\mu F_{\nu\rho}+D_\nu F_{\rho\mu}+D_\rho F_{\mu\nu}=0$。これが :term[Bianchi 恒等式]{en="Bianchi identity" yomi="びあんきこうとうしき"} である。

$\mathcal L_{\rm YM}=-\frac12\tr F_{\mu\nu}F^{\mu\nu}$ は変換に対し不変である。$\delta F_{\mu\nu}=D_\mu\delta A_\nu-D_\nu\delta A_\mu$ なので $\delta S=-2\int\tr F^{\mu\nu}D_\mu\delta A_\nu$。トレースの巡回性と部分積分により運動方程式は $D_\mu F^{\mu\nu}=0$ となる。Abel の場合は交換子が消えて Maxwell 方程式になる。$SU(2)$ で一定な $A_1=aT^1,A_2=bT^2$ を非零成分とする場合は $F_{12}=gabT^3$ となり、場の微分が零でも曲率は残る。

### 場の強さが一定でも方程式を満たすとは限らない

既出の $A_1=aT^1,A_2=bT^2$ では $F_{12}=gabT^3$。空間添字が二つなので $F^{12}=F_{12}$ であり、$D_\mu F^{\mu2}=D_1 F^{12}=-g^2a^2bT^2$。同様に $D_\mu F^{\mu1}=-g^2ab^2T^1$。$ab\ne0$ なら外部源なしの停留解ではない。例えば $a=b=1\,\mathrm{GeV}$、$g=2$ では磁気エネルギー密度は $(gab)^2/2=2\,\mathrm{GeV}^4$、非零の電流成分の大きさは $4\,\mathrm{GeV}^3$ である。比較として $A_2=Bx^1 T^3$、他の成分が零なら $F_{12}=BT^3$ と $[A,F]=0$ から源なしの局所方程式を満たす。空間全体では一定磁場の全エネルギーは有限でない。

### Gauss 拘束が局所ゲージ変換を生成する

時間ゲージ $A_0=0$ を用い、空間成分を下添字で $A_i^a$、共役を $E^{ai}=F_{0i}^a$ と取る。$A_0$ を零にする前の作用の変分で $D_iE^i=0$ を得る。成分が同じ実随伴表現の内積で、コンパクトな台を持つ $\omega$ について
$$
\mathcal G[\omega]=-\int\dd^3x\,\omega^a(D_iE^i)^a
=\int\dd^3x\,E^{ai}(D_i\omega)^a
$$
と置く。部分積分では Lie 代数の不変内積が $[A,\,\cdot\,]$ を反対側へ移す符号を与える。$\{A_i^a,E^{bj}\}=\delta^{ab}\delta_i^j\delta^3$ から $\{A_i,\mathcal G[\omega]\}=D_i\omega$。また $\delta_A\mathcal G=-2ig\int\tr E[\delta A,\omega]$ を巡回させて $\{E^i,\mathcal G[\omega]\}=ig[\omega,E^i]$。従ってこれは電場を随伴場として同時に動かす。

Gauss の量 $D_iE^i$ 自身も随伴に変換する。従って $\{\mathcal G[\omega_1],\mathcal G[\omega_2]\}$ は同じ拘束の線形結合で、拘束面では零である。実際、$\delta_2\mathcal G[\omega_1]=-2\int\tr\omega_1(ig[\omega_2,D_iE^i])$ を巡回すれば $\mathcal G[ig[\omega_1,\omega_2]]$。この計算では変換パラメータを位相空間の変数として変分しない。局所拘束の消去は Maxwell 場と同様にゲージ選択と対にして行うが、ここでの Gauss 拘束の Poisson 括弧は、Lie 代数の交換子をパラメータとする同じ拘束に閉じる。

:::example[異なる色の残留ゲージ回転]
空間一定の $\omega_1=T^1,\omega_2=T^2$ を周期箱で取る。小さい変換は $A_i$ を $ig[\omega,A_i]$ だけ回す。括弧に現れるパラメータは $ig[T^1,T^2]=-gT^3$。例えば第三方向の場 $A_i\propto T^3$ にはこの合成の先頭は作用しないが、第一方向の場 $A_i\propto T^1$ は第二方向へ動く。この差は色の回転が可換でないことによる。
:::

### Yang–Mills 場のエネルギーと拘束

$E_i^a=F_{0i}^a$、$B_i^a=\varepsilon_{ijk}F_{jk}^a/2$ と置く。$\eta=\diag(+1,-1,-1,-1)$ では $\mathcal L=\sum_{a,i}[(E_i^a)^2-(B_i^a)^2]/2$。$E_i=\partial_0 A_i-D_i A_0$ なので、時間微分を含まない $A_0$ の共役は零、$A_i$ の共役は $E_i$。$\sum_{a,i}E_i^a\partial_0 A_i^a-\mathcal L$ に $\partial_0 A_i=E_i+D_i A_0$ を代入し、空間の境界項を捨てると
$$
H=\int\dd^3x\,\left[\frac12(E_i^aE_i^a+B_i^aB_i^a)-A_0^a(D_iE_i)^a\right].
$$
従って $A_0$ は Gauss 拘束を課し、その面では場のエネルギーは非負である。$A_0=0$ を最初に作用へ入れて変分するとこの拘束を取り損なう。境界を持つ系で部分積分した表面項を捨てるには、境界の場または許す変分を指定する必要がある。

$T^{\mu\nu}=-F^{a\mu\lambda} F^{a\nu}{}_\lambda+\eta^{\mu\nu}F^{a\rho\sigma}F^a_{\rho\sigma}/4$ と置く。$T^{00}=(E^2+B^2)/2$。色縮約の積を微分するときは普通微分を二つの共変微分に分けられる。接続に比例する交換子の項が不変内積で相殺するためである。従って $D_\mu F^{\mu\lambda}=0$ と Bianchi 恒等式を用いると $\partial_\mu T^{\mu\nu}=0$。具体的には $F^{\mu\lambda}D_\mu F_{\nu\lambda}=F^{\mu\lambda}D_\nu F_{\mu\lambda}/2$ であり、前者が運動項の発散、後者がトレース項の発散を相殺する。四次元の跡は $-F^2+4F^2/4=0$ である。これは古典の質量なしの作用についての等式であり、量子論での尺度の扱いは [アノマリー](/physics/quantum-field-theory/anomalies/) へつながる。

:::example[固定した色に沿った平面波]
$A_\mu^a T^a=\epsilon_\mu T\,e^{-ikx}$ の実部を取る。全ての成分が一つの生成子 $T$ に比例するので交換子は零。$k^2=0$、$k\cdot\epsilon=0$ なら Maxwell の平面波と同じ方程式を満たす。異なる生成子に比例する二つの波を単純に足すと、重なる領域で $-ig[A_\mu,A_\nu]$ に交差項が出る。従って各波が別々に解でも、その和が解となる保証はない。
:::

### 運動方程式と局所対称性の恒等式

任意の Lie 代数値のコンパクトな台を持つ $\omega$ に対し、$\delta A_\mu=D_\mu\omega$ である。作用の変分を既に求めた形で書くと
$$
0=2\int\dd^4x\,\tr[(D_\nu F^{\nu\mu})D_\mu\omega]
=-2\int\dd^4x\,\tr[\omega D_\mu D_\nu F^{\nu\mu}].
$$
ここで場の方程式は使っていない。$\omega$ を任意に取れるので $D_\mu D_\nu F^{\nu\mu}\equiv0$ が恒等式となる。曲率の反対称性から二つの共変微分を交換子にしても同じ結論を得る。一般の局所対称性は、異なる方程式が微分によって結びつく関係を与える。この例では、時間ゲージを選んだ後も $D_iE^i$ の初期値を指定する必要があるのは、$A_0$ の方程式が他の三成分の方程式とは違って拘束を表すためである。

時間ゲージで $D_\nu F^{\nu i}=0$ を課すと、上の恒等式は $\partial_0 D_\nu F^{\nu0}=0$ を要求する。$D_\nu F^{\nu0}=D_iE_i$ だから、初期時刻で零なら全時刻で零のままである。外部電流がある場合には $D_\mu j^\mu=0$ なら同じ議論を $D_\nu F^{\nu\mu}-j^\mu$ に適用できる。例えば静的な色電荷を指定して空間電流を零とする場合、$\partial_0 j^0=0$ だけでは一般の $A_0$ と両立せず、$[A_0,j^0]=0$ も確かめる必要がある。

### 曲率の電場と磁場から時間発展を組み立てる

時間ゲージで $E_i=F_{0i}=\dot A_i$、$B_i=\varepsilon_{ijk}F_{jk}/2$ とする。運動方程式 $D_\mu F^{\mu i}=0$ では $F^{0i}=-E_i$、$F^{ji}=F_{ji}$ なので $\dot E_i=D_j F_{ji}=-\varepsilon_{ijk}D_j B_k$。Bianchi 恒等式の時間空間成分を使うと $\dot B_i=\varepsilon_{ijk}D_jE_k$。空間三添字の成分は $D_iB_i=0$ を与える。この符号は $E_i=F_{0i}$ と下添字の空間ポテンシャルを使った定義による。成分の色内積でエネルギー密度を $\mathcal E=(E_i^aE_i^a+B_i^aB_i^a)/2$ と書けば、
$$
\partial_0\mathcal E
=-E_i^a\varepsilon_{ijk}(D_jB_k)^a
+B_i^a\varepsilon_{ijk}(D_jE_k)^a
=\partial_i(\varepsilon_{ijk}E_j^aB_k^a).
$$
最後では $\varepsilon$ の反対称性と、色内積の共変微分の積の法則を使った。従ってこの規約のエネルギー流は $-\vb E^a\times\vb B^a$。例えば $A_i$ の三つの成分が同じ色に比例すれば共変回転は普通の回転に戻る。異なる色の電場と磁場を取り、単独に Maxwell の時間発展を足しても、$D_jB_k$ にある非線形項を落とすため一般にはエネルギーの計算が合わない。

:::example[空間一様な $SU(2)$ の非線形振動]
時間ゲージで $A_i(t)=a(t)T^i$ とし、$i=1,2,3$。すると $F_{ij}=g a^2\varepsilon^{ijk}T^k$、$B_i=g a^2T^i$、$E_i=\dot a T^i$。Gauss 拘束では各 $[A_i,E_i]$ が零なので自動的に満たされる。$D_jF_{ji}$ は二つの $j\ne i$ がそれぞれ $-g^2a^3T^i$ を与える。従って $\ddot a=-2g^2a^3$。両辺に $\dot a$ を掛けると $\dot a^2+g^2a^4$ が保存する。時刻を移して $\dot a=0$、$a=a_0>0$ に始めると、$a=0$ で $|\dot a|=g a_0^2$。エネルギー密度は $3(\dot a^2+g^2a^4)/2$。周期箱の一様場として解釈でき、全空間での有限エネルギー解とするには別の境界を課さなければならない。
:::

### 双対曲率と全微分の密度

四次元で $\widetilde F^{\mu\nu}=\varepsilon^{\mu\nu\rho\sigma}F_{\rho\sigma}/2$、$\varepsilon^{0123}=+1$ と置く。Bianchi 恒等式をこの $\varepsilon$ で縮約すると $D_\mu\widetilde F^{\mu\nu}=0$。従って $\int\tr F_{\mu\nu}\widetilde F^{\mu\nu}$ の局所変分は境界項である。実際、双対を付けた縮約は二つの $F$ の間で対称なので
$$
\delta\tr(F_{\mu\nu}\widetilde F^{\mu\nu})
=4\tr(\widetilde F^{\mu\nu} D_\mu\delta A_\nu)
=4\partial_\mu\tr(\widetilde F^{\mu\nu}\delta A_\nu).
$$
最後で $D_\mu\widetilde F^{\mu\nu}=0$ とトレースの巡回性を用いた。従って境界を固定した古典の局所方程式はこの密度を加えても変わらない。$E_i=F_{0i}, B_i=\varepsilon_{ijk}F_{jk}/2$ の成分では $F_{\mu\nu}^a\widetilde F^{a\mu\nu}=4E_i^aB_i^a$。空間反転で電場はベクトル、磁場は軸性なのでこの密度は符号を変える。

局所の接続一形式を $\mathcal A=-ig A$、曲率二形式を $\mathcal F=\dd\mathcal A+\mathcal A\wedge\mathcal A=-ig F$ とする。外積と行列積を共に含めて書くと
$$
\tr(\mathcal F\wedge\mathcal F)=\dd\tr\left(\mathcal A\wedge\dd\mathcal A
+\frac23\mathcal A\wedge\mathcal A\wedge\mathcal A\right).
$$
左辺で $\tr\mathcal A^4=0$ なのは、一形式を三形式の前から後ろへ移すと負号が付くため。交差項は $2\tr(\dd\mathcal A\,\mathcal A^2)$、右の三次の微分は巡回性と外積の符号により $3\tr(\dd\mathcal A\,\mathcal A^2)$。残る $\tr(\dd\mathcal A\,\dd\mathcal A)$ も両辺で同じである。ただし右辺の三形式は接続の局所表示に依存する。境界や座標の張り合わせがある場合には、左辺の積分を零と結論してはならない。その大域的な役割は [非摂動的手法](/physics/quantum-field-theory/nonperturbative/) に続く。

### 線形化して伝播する自由度を数える

真空接続 $A=0$ の近くで曲率の一次だけを残すと、各色について独立な Maxwell 方程式が得られる。時間ゲージで空間の運動量 $\vb k\ne0$ の Fourier モードを見る。Gauss 拘束は $\vb k\cdot\dot{\vb A}^a=0$ なので、縦の共役運動量が零。時間に依存しない残留ゲージ $\delta A_i^a=\partial_i\omega^a$ は、縦の座標そのものを固定できる。従って三つの座標と三つの共役のうち一対が除かれ、各色に二つの横の自由度が残る。横の運動方程式は $\ddot{\vb A}_T^a+|\vb k|^2\vb A_T^a=0$。これは小さい場での線形化の結果であり、箱の一様モードや大域的な残留ゲージは別に扱う。古典の交換子は二次以上に現れるので、線形化で自由に見えることと、四点の自己相互作用が存在することは矛盾しない。

## Wilson 線と Wilson ループ

$x(s)$ に沿って $D_s\psi=0$ を解く。$W_\gamma(s)$ を初期点から $s$ までの移送行列とすると $\partial_s W=ig\dot x^\mu A_\mu W$、$W(0)=1$。反復積分で大きい $s$ の行列を左に置き
$$
W_\gamma=\mathcal P\exp\left(ig\int_\gamma A_\mu\,\dd x^\mu\right)
\label{eq-qft-wilson-line}
$$
を得る。:term[Wilson 線]{en="Wilson line" yomi="うぃるそんせん"} は $U(x_f)W_\gamma U(x_i)^{-1}$ に変わる。この候補を微分方程式と初期条件へ代入すれば確認できる。閉路では端点が一致するので $\tr_R W_\gamma/\dim R$ がゲージ不変な :term[Wilson ループ]{en="Wilson loop" yomi="うぃるそんるーぷ"} である。

向きつき微小矩形では辺を移すと $\partial_\mu A_\nu-\partial_\nu A_\mu$ が、四つの指数の積の交差項では $-ig[A_\mu,A_\nu]$ が現れる。従って $\mu$ 方向から $\nu$ 方向へ進む微小閉路の $W$ は $1+igF_{\mu\nu}\delta x^\mu\delta x^\nu$ までで曲率を測る。例えば Cartan 方向に一定の $F_{12}=BT^3$ を持つ $SU(2)$ の可換な背景では面積 $S$ のループが $\cos(gBS/2)$ を与える。

### 経路の合成と向きの反転

$\gamma_1$ の後に $\gamma_2$ を進む場合、常微分方程式の初期値を第二区間の始点に移すと $W_{\gamma_2\circ\gamma_1}=W_{\gamma_2}W_{\gamma_1}$。逆経路をつなげると逐次的に逆写像となるため $W_{\gamma^{-1}}=W_\gamma^{-1}$ である。また $A$ が Hermite なので $\partial_s(W^\dagger W)=0$、初期値から $W^{-1}=W^\dagger$。例えば独立な線分の移送が $W_1=e^{i\pi\sigma_1/4},W_2=e^{i\pi\sigma_2/4}$ の場合、順序の違いは Pauli 積の $\sigma_3$ 項を反転する。閉じた交換子 $W_2^{-1}W_1^{-1}W_2W_1$ の規格化トレースは $1/2$ であり、各線分が可換なら同じ交換子は $1$ になる。

### 微小矩形の積を二次まで展開する

$A_\mu$ と $A_\nu$ を始点 $x$ で評価し、辺の長さを $a,b$ とする。まず $+\mu$、次に $+\nu$、次に $-\mu$、最後に $-\nu$ に進む。二辺目の評価値は $A_\nu+a\partial_\mu A_\nu$、逆向きの $\mu$ 辺では $A_\mu+b\partial_\nu A_\mu$ となる。順序は移送行列の左から後の辺を掛ける。微分の二次は $igab(\partial_\mu A_\nu-\partial_\nu A_\mu)$ を与える。一方、各辺で接続を一定に置いた積は
$$
e^{-ig b A_\nu}e^{-ig a A_\mu}e^{ig b A_\nu}e^{ig a A_\mu}
=1+g^2ab[A_\mu,A_\nu]+O(a^2b+ab^2).
$$
一次は打ち消し合い、残る交差項の符号は因子の順序で決まる。二つを足せば $igab F_{\mu\nu}$。閉路の向きを逆にすると移送行列が逆になり、この二次項は反転する。

### 平坦な場でも閉路に情報が残る例

長さ $L$ の空間円で、周期なゲージ変換だけを考える。$SU(2)$ の一定接続 $A_1=aT^3$ は曲率零である。一周の移送は $W=e^{ig a L T^3}$、規格化トレースは $\cos(g a L/2)$。$A_1$ を零へ移す局所変換は $U(x)=e^{-igaxT^3}$ であり、$U(L)=U(0)$ の条件は $g a L\in4\pi\Z$。例えば $g a L=\pi$ ではトレースは零であり、$A_1=0$ のトレース一とはゲージ同値でない。曲率零から局所的な純ゲージを作れても、空間の境界条件は大域的な変換の可否を決める。この例では基本表現にも作用する周期変換を採ったため、中心元まで一致する変換を許す別のゲージ群を選ぶ場合には境界条件を再検討する。

### 二つの端点を物質場で閉じる

同じ表現の物質場 $\psi$ を用い、$\bar\psi(x_f)W_\gamma\psi(x_i)$ を作る。始点の $U(x_i)^{-1}$ は $\psi(x_i)$ の変換と、終点の $U(x_f)$ は $\bar\psi(x_f)$ の変換と相殺する。従って始点と終点が異なってもこの双局所挿入はゲージ不変である。一方、開いた $W_\gamma$ だけのトレースは左右の変換が異なるので一般に不変でない。例えば $W_\gamma=1$ の接続に $U(x_i)=1$、$U(x_f)=e^{i\pi T^3}$ の変換を行うと、基本表現の規格化トレースは一から零へ変わる。端点に物質を置いた場合は同じ位相が物質側で消える。

### 経路を少し動かして曲率を挿入する

滑らかな経路 $x(s)$、$0\le s\le1$ と、変形 $\delta x^\mu(s)$ を考える。$A_s=\dot x^\mu A_\mu(x)$ とし、部分移送 $W(s_2,s_1)$ を用いる。運動方程式を一回変分して $W(1,s)$ を左から掛け、$s$ を積分すると
$$
\delta W(1,0)=ig\int_0^1\dd s\, W(1,s)\delta A_s W(s,0).
$$
ここで $\delta A_s=\dd(\delta x^\nu A_\nu)/\dd s+\delta x^\nu \dot x^\mu(\partial_\nu A_\mu-\partial_\mu A_\nu)$。全微分の項を部分積分すると、$W(1,s)$ と $W(s,0)$ の微分がそれぞれ $-igW(1,s)A_s$ と $ig A_s W(s,0)$ なので、交換子の項が生じる。従って
$$
\begin{aligned}
\delta W(1,0)
={}&ig\,\delta x_f^\mu A_\mu(x_f) W(1,0)
-ig\, W(1,0)\delta x_i^\mu A_\mu(x_i)\\
&+ig\int_0^1\dd s\,W(1,s) F_{\nu\mu}(x(s)) W(s,0)\,
   \delta x^\nu\dot x^\mu .
\end{aligned}
$$
挿入点の曲率を両端へ移送する行列が必要なのは、異なる位置で色の基底が独立に回るためである。固定端点で平坦な接続なら微小変形に対する $\delta W$ は零。従ってホモトピーで連続に動かした経路の移送は同じである。ただしつながらないホモトピー類は異なる値を持ち得る。閉路で基準点も動く場合、端の二項はトレースで消える。

:::example[可換な平行色の二領域]
長方形の中で $F_{12}=B(x^1,x^2)T^3$、移送も同じ $T^3$ で可換な $SU(2)$ の背景を考える。曲率挿入は包んだ面積の変分に積分され、$W=\exp[igT^3\int_S B\,\dd x^1\dd x^2]$。例えば二つの等面積 $S/2$ の領域で $B_1=-B_2$ なら、一周の規格化トレースは一である。それぞれの曲率が零という意味ではなく、この経路での全磁束が打ち消し合った結果である。色の方向が異なる場合には、面の二箇所の値を持ち運んだ行列が可換と限らず、同じ単純な面積積分を使えない。
:::

### 長い時間の閉路で短距離のエネルギーを読む

先頭の摂動論では Yang–Mills の二次作用は、一色ごとに前の Abel 自由場と同じである。Euclid の Feynman ゲージで $\langle A_\mu^a(x)A_\nu^b(y)\rangle=\delta^{ab}\delta_{\mu\nu}D(x-y)$、$D(x)=1/(4\pi^2|x|^2)$。この $D$ は $1/k^2=\int_0^\infty e^{-u k^2}\dd u$ と四次元 Fourier Gaussian から得る。規格化した閉路 $W_R=\tr_R W_\gamma/\dim R$ を $g^2$ まで展開する。順序つき二点の積はこの次数ではトレースの二生成子を対称にできるので、全経路積分の半分を取ればよい。従って
$$
\log\langle W_R\rangle=-\frac{g^2 C_2(R)}2
\oint_\gamma\dd x_\mu\oint_\gamma\dd y_\mu\,D(x-y)+O(g^4).
$$
同じ線上の近接点の積分は正則化し、外部源一個の自己エネルギーとして分ける。ここでは距離 $r$ に依存する先頭だけを求める。時間幅 $T$、空間間隔 $r>0$ の矩形では、二本の時間線の向きが逆である。従って交差部分は $-2I(T,r)$、ただし
$$
I(T,r)=\frac{2}{4\pi^2}\int_0^T
\frac{T-u}{u^2+r^2}\,\dd u
=\frac{T}{2\pi^2 r}\arctan\frac T r
-\frac{1}{4\pi^2}\log(1+T^2/r^2).
$$
時間線と空間線はこのゲージでは計量の縮約が零。二本の空間線の積は $T\to\infty$ で線形の $T$ を作らない。従って $-T^{-1}\log\langle W_R\rangle$ の $r$ に依存する部分は $-g^2 C_2(R)/(4\pi r)$。反対向きの二つの源を閉じて色一重項にした短距離ポテンシャルであり、同じ係数を交換振幅からも得る。閉路の接合部に必要な正則化や、強結合での長距離の面積則は、この先頭の計算からは決まらない。

:::example[同じ形の閉路で表現を替える]
$SU(3)$ で基本表現と随伴表現の矩形を同じ $r,T$ で比べる。$T$ に比例し距離に依存する先頭の係数は $C_F=4/3$、$C_A=3$ なので比は $9/4$。これは同じ二次の伝播関数に二生成子を掛けた結果である。高次の色因子や大きな距離でこの比が保存されることまでは、この二次展開は主張しない。
:::

### 表現を変えた閉路のトレース

$SU(2)$ の移送を基本表現で $W=\cos\zeta\,1+i\sin\zeta\,\widehat{\vb n}\cdot\boldsymbol\sigma$ と書く。その固有値は $e^{\pm i\zeta}$。従って角運動量 $j$ の表現の固有値は $e^{2im\zeta}$、$m=-j,\ldots,j$ である。$W$ は $e^{2i\zeta T^3}$ と共役であり、$T^3$ の固有値を各表現で用いた。規格化したトレースは幾何級数から
$$
\frac{\tr_j W}{2j+1}
=\frac{\sin[(2j+1)\zeta]}{(2j+1)\sin\zeta} .
$$
分母が零の端点は左辺の有限な和で延長する。例えば $\zeta=\pi/2$ では基本の値は零、随伴 $j=1$ は $-1/3$。一方 $\zeta=\pi$ は中心元であり、整数 $j$ では $+1$、半整数 $j$ では $-1$。随伴の閉路観測量だけでは基本表現の二つの中心を区別できない。局所の $F$ に対する同じ作用から出発しても、閉路の表現が持つ大域的な情報は異なる。

## 物質場との結合とゲージ群の表現

物質の表現 $R$ では $A_\mu^a T_R^a$ を用いる。$\bar\psi(i\not D-m)\psi$ は $U$ と $U^{-1}$ が相殺して不変である。複素スカラーでは $(D_\mu\phi)^\dagger D^\mu\phi-V(\phi)$ を用いる。四次元では $[g]=0$ である。Dirac の結合から $j^{a\mu}=-g\bar\psi\gamma^\mu T_R^a\psi$ と定義すれば $D_\nu F^{\nu\mu}=j^\mu$。$D_\mu D_\nu F^{\nu\mu}$ は微分の交換子の半分であり、同じ曲率同士の交換子なので零になる。従って $D_\mu j^\mu=0$ である。例えば随伴表現の場は $X'=UXU^{-1}$ で変換し、ゲージ一重項なら $T_R^a=0$ なのでこの直接結合を持たない。$SU(3)$ の基本表現の色を平均した単一生成子の二乗和は $(1/3)\sum_a\tr(T^aT^a)=4/3$ である。

### 静的な粒子対の色電流

基本表現の重い粒子と共役表現の反粒子を摂動的な外部源として扱う。単独の色荷のゲージ不変な漸近状態は仮定しない。木の一ゲージ粒子交換で非相対論的な時間電流だけを残すと $V(r)=g^2(\sum_a T_1^a T_2^a)/(4\pi r)$ となる。時間電流と外脚規格化は Abel の静的極限と同じである。$T_{\rm tot}^a=T_1^a+T_2^a$ の二乗から $\sum_a T_1^aT_2^a=(C_{\rm tot}-2C_F)/2$。従って一重項では $V=-g^2 C_F/(4\pi r)$、随伴チャネルでは $V=g^2/(8\pi Nr)$。例えば $SU(3)$ の二つの係数は共通の $g^2/(4\pi r)$ を取って $-4/3$ と $1/6$ である。この式は結合が小さい短距離の近似である。長距離の Wilson ループは [非摂動的手法](/physics/quantum-field-theory/nonperturbative/) で扱う。

### 二つの同じ基本電荷の色チャネル

二基本表現のテンソル積で :eqref[eq-qft-su-completeness] を作用させる。$\sum_a T^a\otimes T^a=(P-1/N)/2$ であり、$P$ は二つの色添字を交換する。従って対称チャネルでは $(1-1/N)/2$、反対称では $(-1-1/N)/2$。$N=3$ では順に $1/3,-2/3$ である。短距離の静的ポテンシャルではこれに $g^2/(4\pi r)$ を掛けるので、対称色の交換は斥力、反対称色の交換は引力となる。この色の交換はフェルミ統計に課される全波動関数の交換とは区別して計算する。

### スカラー物質の電流と二光子型の接触項

$D_\mu\phi=(\partial_\mu-igA_\mu^aT_R^a)\phi$ とすると、密度の $A$ に一次の部分は
$$
ig A_\mu^a[\phi^\dagger T_R^a\partial^\mu\phi
-(\partial^\mu\phi)^\dagger T_R^a\phi],
$$
二次は $g^2 A_\mu^aA^{b\mu}\phi^\dagger T_R^a T_R^b\phi$ である。従って二つのゲージ脚について微分すると $ig^2\eta^{\mu\nu}\{T_R^a,T_R^b\}$ を得る。Abel の電荷に帰着すれば反交換子は数の二倍であり、スカラー QED の局所項と一致する。非可換では異なる色の二つの順序が行列として別なので、先に一色に固定した式を二色へ無条件に延長できない。

電流を $j^{a\mu}=-\partial\mathcal L_{\rm matter}/\partial A_\mu^a$ と定めると、上のスカラーの場合は $-ig[\phi^\dagger T_R^aD^\mu\phi-(D^\mu\phi)^\dagger T_R^a\phi]$。この式の $D$ を展開すれば二次密度の変分も含む。有限変換では二つの物質場と $D\phi$ が同じ $U_R$ で変わるため、$j$ は随伴に変換する。随伴の共変発散を取ると $D_\mu\phi$ を二つ含む項が相殺し、$\phi^\dagger T_R^aD^2\phi-(D^2\phi)^\dagger T_R^a\phi$ だけが残る。$V(\phi^\dagger\phi)$ の方程式では $D^2\phi=-V'(\phi^\dagger\phi)\phi$ で $V'$ は実数なので、二項も相殺する。これは Yang–Mills 方程式から従う整合条件 $D_\mu j^\mu=0$ を物質だけから示したものである。

:::example[基本二成分で二色を選ぶ]
$SU(2)$ の基本表現では $\{T^a,T^b\}=\delta^{ab}1/2$。従って色 $1$ と色 $2$ のスカラー接触頂点は零、一方で順に入る二つの三点頂点の積は $T^1T^2=iT^3/2$ で非零である。二色の順序を逆にすると積の符号が反転する。接触項が零であることだけから、散乱の各逐次挿入まで零になるとは言えない。
:::

### 軌道に沿って色荷が回る古典模型

時空の曲線 $x^\mu(\tau)$ に沿い、Lie 代数値の色荷 $Q(\tau)=Q^a T^a$ を随伴の平行移動に従わせる。従って
$$
\dot Q-ig[\dot x^\mu A_\mu(x(\tau)),Q]=0 .
$$
解は $Q(\tau)=W(\tau,\tau_0)Q(\tau_0)W(\tau,\tau_0)^{-1}$。この規約では $j^\mu(x)=g\int\dd\tau\,\dot x^\mu Q(\tau)\delta^4(x-x(\tau))$ は共変保存する。実際、任意の試験関数へ積分した $\partial_\mu j^\mu$ では $\dot x^\mu\partial_\mu$ が曲線に沿う微分を与え、部分積分で $g\int\dot Q\,\delta^4$。交換子の項を加えると上の方程式が出る。端点の寄与は曲線が試験関数の台を通って始終する場合には残るので、保存を主張するのは閉じた曲線か、端点が台の外にある場合である。

$\dot x^\mu A_\mu=aT^3$ が一定の $SU(2)$ の背景なら、$Q(0)=T^1$ から $Q(\tau)=T^1\cos(ga\tau)-T^2\sin(ga\tau)$。$\tr Q^2$ は一定である。成分の変化は局所的な色基底に依存し、軌道上では電荷は共変に一定。一定なベクトル電荷をそのまま全てのゲージで固定する模型では $D_\mu j^\mu=0$ を満たさないことがある。Wilson 線が外部の重い色源を記述するのは、この移送と電流の条件を一緒に含むためである。

### 色対を射影して交換を対角化する

基本表現と共役の積の状態を行列 $X_{ij}$ で表す。第一の生成子は左から $T^a$、第二は右から $-T^a$ を掛ける。従って交換色演算子は
$$
QX=-\sum_a T^a X T^a
=-\frac12\tr(X)1+\frac{X}{2N}.
$$
一重項への射影は $P_1X=\tr(X)1/N$、随伴の射影は $P_{\rm adj}X=X-P_1X$。実際、トレース内積で $1$ と無跡行列は直交し、それぞれ $U X U^{-1}$ で保たれる。従って $Q=-C_F P_1+(1/(2N))P_{\rm adj}$。例えば $X=E_{ij}$、$i\ne j$ なら随伴に属するので交換係数は $1/(2N)$。同じ添字の $X=E_{ii}$ は $\tr X=1$ なので二つの射影を両方持つ。色と反色のラベルが同じことだけでは、一重項であると判断できない。

:::example[色を選ばずに二乗を平均する]
任意の $F\otimes\bar F$ の入射について色を均等に平均した交換演算子の二乗は $\tr Q^2/N^2$。射影で評価すると
$$
\frac1{N^2}\left[C_F^2+\frac{N^2-1}{4N^2}\right]
=\frac{N^2-1}{4N^2}.
$$
例えば $N=3$ なら $2/9$。一重項だけを選んだ場合の $C_F^2=16/9$ より小さい。直接トレースを計算しても $\sum_{a,b}\tr(T^aT^b)\tr(T^{a*}T^{b*})=(N^2-1)/4$。この平均は全ての色の状態を等しい確率で準備したという仮定を含むので、閉路で選んだ一重項の物理量にはそのまま使わない。
:::

### 物質へ移るエネルギーと力

物質電流で $D_\mu F^{\mu\lambda}=j^\lambda$ とする。場だけの対称な $T^{\mu\nu}$ の発散を計算したとき、Bianchi 恒等式で相殺する二項は源があっても同じ。一方、場の方程式を使う項は $-j^{a\lambda}F^{a\nu}{}_\lambda$ となる。従って全系の平行移動対称性で $\partial_\mu T_{\rm tot}^{\mu\nu}=0$ なら、物質の側は $+j^{a\lambda}F^{a\nu}{}_\lambda$ を受ける。これは色内積を取った量であり、両因子を同じ随伴変換で回しても変わらない。例えば色荷の軌道の模型 $j^\lambda=g\int Q^a\dot x^\lambda\delta^4$ と点粒子のテンソルを入れると、質量 $m$ の固有時で $m\ddot x^\nu=gQ^aF^{a\nu}{}_\lambda\dot x^\lambda$ を得る。速度と縮約すると $F$ の反対称性により $m\dot x^2$ の微分は零。色荷の長さの保存は別に、軌道の共変平行移動から従う。

例えば一定な電場が一色 $T^3$ だけを持ち、色荷も同じ方向に揃えば、軌道の色はその中で回らない。力の式は電荷を $gQ^3$ とした Lorentz 力に戻る。一方 $Q$ が直交する色方向だけを持ち、接続も常に $T^3$ に沿う場合は、色荷は回っても $Q^3$ は零のまま。従ってこの外部場による力は零である。単に個々の色成分が時間発展することと、色を縮約した力の存在は区別できる。

## 古典的 Yang–Mills 理論の相互作用頂点

$F^a_{\mu\nu}=f^a_{\mu\nu}+g f^{abc}A_\mu^b A_\nu^c$、$f^a_{\mu\nu}=\partial_\mu A_\nu^a-\partial_\nu A_\mu^a$ と分けると
$$
\begin{aligned}
\mathcal L_3&=-g f^{abc}(\partial_\mu A_\nu^a)A^{b\mu}A^{c\nu},\\
\mathcal L_4&=-\tfrac14 g^2 f^{abc}f^{ade} A_\mu^b A_\nu^c A^{d\mu}A^{e\nu}.
\end{aligned}
\label{eq-qft-ym-vertices}
$$
第一式は曲率二乗の二つの交差項と $f^a_{\mu\nu}$ の反対称性から得る。第二式は交換子部分の二乗である。運動量頂点は各外部脚の $(a,\mu,p)$ を全ての場に割り当てる。全入射を $e^{-ipx}$ とすれば三点の各割付けは $i(-g)f^{abc}(-ip_\mu)$ と計量で決まり、六通りを集めると
$$
g f^{abc}\bigl[\eta^{\mu\nu}(p-q)^\rho+\eta^{\nu\rho}(q-r)^\mu+\eta^{\rho\mu}(r-p)^\nu\bigr]
$$
を得る。四点では四つの場へ $4!$ 通りを代入し、同じ色の交換子対をまとめる。例えば全脚が同じ Cartan 方向なら $f^{a b b}=0$ なので三点も四点も消える。非可換の自己相互作用は、物質がない理論でもゲージ量子同士を結ぶ。

### 四つの脚の色を固定する

全て入射の色と Lorentz 添字を $(a,\mu),(b,\nu),(c,\rho),(d,\sigma)$ とする。:eqref[eq-qft-ym-vertices] の四点項では、第一と第二の交換子のどちらに脚対 $(a,b)$ を入れるかで二通り、各交換子内の順序で二通りずつある。前者は同じ項を与え、後者は構造定数と計量の反対称な組を作る。三つの脚対の分割を集めた頂点は
$$
\begin{aligned}
-ig^2\sum_e\big[&
 f^{eab}f^{ecd}(\eta^{\mu\rho}\eta^{\nu\sigma}-\eta^{\mu\sigma}\eta^{\nu\rho})\\
&+f^{eac}f^{ebd}(\eta^{\mu\nu}\eta^{\rho\sigma}-\eta^{\mu\sigma}\eta^{\rho\nu})\\
&+f^{ead}f^{ebc}(\eta^{\mu\nu}\eta^{\sigma\rho}-\eta^{\mu\rho}\eta^{\sigma\nu})\big].
\end{aligned}
$$
例えば $SU(2)$ で $A_1^1=a$、$A_2^2=b$ だけが非零の一様場なら $F_{12}^3=gab$、$F_{21}^3=-gab$。従って作用密度の四次部分は $-g^2a^2b^2/2$ であり、四点係数は $i\partial_a^2\partial_b^2\mathcal L_4=-2ig^2$。上の色・計量の式に $(1,1),(1,1),(2,2),(2,2)$ を入れた値と一致する。

### 三点頂点の運動量縮約

色因子 $g f^{abc}$ を外した三点頂点を $V^{\mu\nu\rho}(p,q,r)$ とする。$p+q+r=0$ を使って $p_\mu V^{\mu\nu\rho}$ を展開すると、計量の項は $(r^2-q^2)\eta^{\nu\rho}$、残りは $q^\nu q^\rho-r^\nu r^\rho$ になる。従って
$$
p_\mu V^{\mu\nu\rho}
=(r^2\eta^{\nu\rho}-r^\nu r^\rho)-(q^2\eta^{\nu\rho}-q^\nu q^\rho).
$$
各括弧は運動項の横の逆核である。例えば $q,r$ をそれぞれ質量ゼロの外脚に置き、それらに横偏極をつなぐと両括弧の縮約が消える。一方、内部線につなぐと逆核が伝播関数を消すため、四点接触項や他チャネルとの関係が必要になる。

### 頂点の対称性と非可換な色因子

三点の運動量テンソル $V^{\mu\nu\rho}(p,q,r)$ は二つの脚を入れ替えると符号が反転する。例えば $(\mu,p)$ と $(\nu,q)$ を入れ替えると、第一項は $q-p$、残り二項は互いに入れ替わり符号が変わる。従って同時に $f^{abc}$ の二つの色を入れ替えた全頂点は不変である。同じ結論は、古典作用を三つの可換な場の成分で微分する順序に依存しないことからも従う。四点では例えば第一の色積 $f^{eab}f^{ecd}$ と第一の計量差は、$a,b$ および対応する Lorentz 添字の交換でそれぞれ反転する。残り二分割はその交換で互いに移り、全体がやはり不変となる。

三つの交換の色テンソルは Jacobi 恒等式
$$
f^{abe}f^{ecd}+f^{bce}f^{ead}+f^{cae}f^{ebd}=0
$$
を満たす。実際、$[[T^a,T^b],T^c]+[[T^b,T^c],T^a]+[[T^c,T^a],T^b]=0$ の $T^d$ の係数である。例えば $SU(2)$ の $a,b,c,d=1,2,1,2$ では三項は $1,-1,0$。従って三つの色分割を独立な三種類の結合として調整することはできない。四次の頂点の全体が、二次と三次と共通の $g$ を持つ局所作用から出る。

### 四外脚の木で縦の成分を相殺する

四つのゲージ波を全て入射の規約で $p_i,\epsilon_i,a_i$ とし、$p_i^2=0$、$\epsilon_i\cdot p_i=0$、$\sum_i p_i=0$。二つの三点を結んだ三図と四点接触があり、基準として自由な Feynman ゲージの伝播関数を用いる。色の略記を $C_s=f^{a_1a_2e}f^{ea_3a_4}$、$C_t=f^{a_1a_3e}f^{ea_2a_4}$、$C_u=f^{a_1a_4e}f^{ea_2a_3}$ とする。Jacobi より $C_s-C_t+C_u=0$ と分かる。共通の規約で全入射の運動量を使うので、これらの符号は三点の脚順と対になっている。

第一の偏極を $p_1$ に替える。$s$ 図の内部運動量を $P_s=p_1+p_2$ とすると、第一の三点の縮約は $P_s^2\epsilon_2^\alpha-(P_s\cdot\epsilon_2)P_s^\alpha$。第二の三点では $p_3,p_4$ が横かつ殻上なので $P_s$ の縮約が零。従って伝播関数の $P_s^2$ が消え、共通の $-ig^2C_s$ を除いた残りは
$$
\begin{aligned}
W_s={}&-2(\epsilon_2\cdot\epsilon_3)(p_3\cdot\epsilon_4)
+(\epsilon_3\cdot\epsilon_4)((p_3-p_4)\cdot\epsilon_2)\\
&+2(\epsilon_2\cdot\epsilon_4)(p_4\cdot\epsilon_3).
\end{aligned}
$$
これは第二の三点 $V(P_s,p_3,p_4)$ の第一脚に $\epsilon_2$ をつないで得る。四点の同じ色の項は $K_s=(p_1\cdot\epsilon_3)(\epsilon_2\cdot\epsilon_4)-(p_1\cdot\epsilon_4)(\epsilon_2\cdot\epsilon_3)$。従って $p_1=-(p_2+p_3+p_4)$ を使うと、
$$
\begin{aligned}
W_s+K_s={}&(\epsilon_2\cdot\epsilon_3)((p_2-p_3)\cdot\epsilon_4)
+(\epsilon_3\cdot\epsilon_4)((p_3-p_4)\cdot\epsilon_2)\\
&+(\epsilon_4\cdot\epsilon_2)((p_4-p_2)\cdot\epsilon_3).
\end{aligned}
$$
これは $2,3,4$ の巡回で同じ、交換で反転する。従って $t$ 図と対応する接触の組はこれの負、$u$ 図と対応する組は同じ。全色の和は $C_s-C_t+C_u$ を含むので消える。内部運動量が零の極では振幅の処方を別に取るが、一般のチャネルでこの縮約は恒等的に零である。四点接触と三つの交換を全て含めることで、古典作用から得た頂点のゲージの関係が木の振幅に現れる。量子のゲージ条件とループの恒等式は [非可換ゲージ理論の量子化](/physics/quantum-field-theory/gauge-quantization/) に続く。

:::example[全偏極が散乱面に垂直な四外脚]
全ての外部運動量を同じ空間平面に置き、$\epsilon_i=n$ をその法線の純空間ベクトルで $n^2=-1$ とする。四点の計量差は同じ二つの積の差なので零。三点では二つの外部偏極の内積だけが残り、$s$ 図の二つの電流は $-(p_1-p_2)$ と $-(p_3-p_4)$。従って $s$ 図は $i\mathcal M_s=-ig^2 C_s(t-u)/s$、ここで $s=(p_1+p_2)^2$、$t=(p_1+p_3)^2$、$u=(p_1+p_4)^2$。$SU(2)$ で色 $(1,2,1,2)$ を選ぶと $C_s=1,C_t=0,C_u=-1$。従って $\mathcal M=-g^2[(t-u)/s-(s-t)/u]$。例えば $t=u=-s/2$ では $\mathcal M=-3g^2$。Cartan の同じ色に全てを替えると三つの色係数も零になり、この木の散乱は消える。
:::

## 演習問題

::::exercise{#exr-qft-gauge-theories-1 level=basic}
$SU(2)$ で $T^a=\sigma^a/2$ とする。$[T^1,T^2]$ を求めよ。

:::solution
$\sigma^1\sigma^2=i\sigma^3$、逆順は $-i\sigma^3$ である。四で割ると $[T^1,T^2]=iT^3$。
:::
::::

::::exercise{#exr-qft-gauge-theories-2 level=basic}
$SU(4)$ の基本表現で $C_F$ を求めよ。

:::solution
$C_F=(4^2-1)/(2\cdot4)=15/8$ である。
:::
::::

::::exercise{#exr-qft-gauge-theories-3 level=basic}
一定の $A_1=aT^1,A_2=bT^2$ の $F_{12}$ を求めよ。

:::solution
微分が消え、$-igab[T^1,T^2]=gabT^3$ が残る。
:::
::::

::::exercise{#exr-qft-gauge-theories-4 level=basic}
Abel の閉路で経路順序が不要な理由を述べよ。

:::solution
接続の値がどの位置でも可換なので、順序領域の積分が全領域の $1/n!$ 倍である。従って通常の指数が出る。
:::
::::

::::exercise{#exr-qft-gauge-theories-5 level=advanced}
$U(x)$ が局所的な場合の $W_\gamma$ の変換則を示せ。

:::solution
$W'(s)=U(x(s))W(s)U(x_i)^{-1}$ を微分する。$\partial_s W'=[(\partial_s U)U^{-1}+ig U A_s U^{-1}]W'=igA_s'W'$ で初期値は一なので一意解である。
:::
::::

::::exercise{#exr-qft-gauge-theories-6 level=advanced}
:eqref[eq-qft-su-completeness] を使って無跡な $X$ への二重交換子を評価せよ。

:::solution
$\sum[T^a,[T^a,X]]=2C_FX-2\sum T^a X T^a=(2C_F+1/N)X=NX$ である。
:::
::::

::::exercise{#exr-qft-gauge-theories-7 level=advanced}
曲率の二重共変発散が消えることを示せ。

:::solution
$D_\mu D_\nu F^{\nu\mu}=\frac12[D_\mu,D_\nu]F^{\nu\mu}=-ig[F_{\mu\nu},F^{\nu\mu}]/2$ である。対角計量では各組について $F^{\nu\mu}$ が $F_{\mu\nu}$ の定数倍なので零である。
:::
::::

::::exercise{#exr-qft-gauge-theories-8 level=advanced}
$SU(2)$ の一定可換背景で $gBS=\pi$ の Wilson ループを求めよ。

:::solution
$T^3$ の固有値は $\pm1/2$。二つの位相の平均は $\cos(\pi/2)=0$ である。
:::
::::

::::exercise{#exr-qft-gauge-theories-9 level=advanced}
同じ二つの基本表現の積の対称成分と反対称成分について、$SU(N)$ の二次 Casimir を求めよ。$N=3$ では反対称成分が共役基本表現と同じ Casimir を持つことを確かめよ。

:::solution
全生成子は $T^a\otimes1+1\otimes T^a$。二乗を和にすると $2C_F+2\sum_aT^a\otimes T^a$ である。:eqref[eq-qft-su-completeness] から後半は $P-1/N$。従って対称と反対称の値は
$$
C_{\rm S}=\frac{N^2+N-2}{N},\qquad
C_{\rm A2}=\frac{N^2-N-2}{N}.
$$
$N=3$ で順に $10/3,4/3$。反対称成分は $\epsilon_{ijk}$ で二つの上添字を一つの下添字へ写せる。$U\in SU(3)$ では $\epsilon$ が $\det U=1$ により不変なので、これは共役基本への同変写像である。その Casimir $4/3$ と一致する。
:::
::::

::::exercise{#exr-qft-gauge-theories-10 level=advanced}
$SU(2)$ の Wilson 閉路の基本表現の規格化したトレースが $1/2$ とする。随伴の規格化したトレースを求めよ。

:::solution
基本の固有値は $e^{\pm i\zeta}$ で、規格化したトレースは $\cos\zeta$。随伴の固有値は $1,e^{\pm2i\zeta}$。従って随伴の規格化したトレースは $(1+2\cos2\zeta)/3=(4\cos^2\zeta-1)/3$。$\cos\zeta=1/2$ を入れると零である。基本のトレースを $-1/2$ に替えても随伴は零であり、中心の $-1$ を掛ける変換が随伴では単位であることと一致する。
:::
::::

::::exercise{#exr-qft-gauge-theories-11 level=advanced}
空間一様な時間ゲージの $A_i(t)=a(t)T^i$ を $a(0)=a_0>0,\dot a(0)=0$ から出発させる。最初に $a=0$ へ達する時間の $a_0,g$ に対する依存を求めよ。

:::solution
$\ddot a=-2g^2a^3$ の初期加速度は負なので、その区間では $\dot a=-g\sqrt{a_0^4-a^4}$。従って到達時間は
$$
T_0=\frac1g\int_0^{a_0}\frac{\dd a}{\sqrt{a_0^4-a^4}}
=\frac1{g a_0}\int_0^1\frac{\dd u}{\sqrt{1-u^4}}
=\frac{B(1/4,1/2)}{4ga_0}.
$$
最後は $z=u^4$ と置いた Beta 積分。初期の折り返し点の特異性は $(1-u)^{-1/2}$ なので積分できる。初期振幅を二倍にすると時間は半分。四次のポテンシャルによる振動であり、振幅に依存しない調和振動子の周期と異なる。
:::
::::

::::exercise{#exr-qft-gauge-theories-12 level=advanced}
四つの $SU(2)$ のゲージ波の色を $(1,2,1,2)$、偏極を全て散乱面の法線に取る。重心で $\theta=\pi/3$ の木の振幅を求めよ。

:::solution
$p_3,p_4$ を全入射にする符号を含め、$t=-s(1-\cos\theta)/2=-s/4$、$u=-s(1+\cos\theta)/2=-3s/4$。同じ偏極なので四点接触は消え、色は $C_s=1,C_t=0,C_u=-1$。本文の $\mathcal M=-g^2[(t-u)/s-(s-t)/u]$ に入れると第一項が $1/2$、第二項の比は $(5s/4)/(-3s/4)=-5/3$。従って $\mathcal M=-13g^2/6$。この値は偏極と色を指定した振幅であり、断面積にするときは初期状態の平均を行う測定かどうかを別に指定する。
:::
::::

::::exercise{#exr-qft-gauge-theories-13 level=advanced}
二実ベクトル $u,v$ を用い、$SU(2)$ の同じ空間領域に $A_1=u^aT^a$、$A_2=v^aT^a$ の一定成分だけがあるとする。曲率の色不変なノルムと、場が平坦になる条件を求めよ。

:::solution
普通の微分が零で $F_{12}$ の交換子部分だけが残る。この章の Hermite 成分の記法で $F_{12}=g (u\times v)^aT^a$。従って $2\tr F_{12}^2=g^2|u\times v|^2=g^2(|u|^2|v|^2-(u\cdot v)^2)$。$g>0$ で平坦になるのは二つが線形従属のときで、どちらかが零の場合を含む。例えば $|u|=|v|=a$ で角度が $\pi/3$ ならノルム二乗は $3g^2a^4/4$。物理的なエネルギー密度には $F_{21}=-F_{12}$ も数えた作用の係数が入り、磁気の寄与は $3g^2a^4/8$ となる。局所で測る曲率が零でも、周期境界での Wilson 移送が自動的に単位とは限らない。
:::
::::

## 確認テスト

::::quiz
曲率の局所ゲージ変換則を選べ。

- [x] $F'=UFU^{-1}$
- [ ] $F'=F+\partial U$
- [ ] $F'=UF$

:::explanation
$D'=UDU^{-1}$ を交換子の積に代入すると中間の $U^{-1}U$ が恒等写像になる。従って曲率は両側から共役される。
:::
::::

::::quiz{answer="3" tol="1e-6"}
$SU(3)$ の随伴 Casimir $C_A$ を答えよ。

:::explanation
基本生成子の完全性で二重交換子が $NX$ となるため、$N=3$ では三である。
:::
::::

::::quiz
閉路のトレースが端点のゲージ変換に依存しない理由を選べ。

- [x] 同じ端点の $U$ と $U^{-1}$ がトレースの巡回性で消える
- [ ] 開路でも両端の変換は常に同一
- [ ] 接続行列が位置によらない

:::explanation
閉路では端点が同じなので $W'=U WU^{-1}$ である。トレースは共役で不変である。
:::
::::

::::quiz{answer="0.75" tol="1e-6"}
$SU(2)$ の基本表現の $C_F$ を答えよ。

:::explanation
$\sum_a(\sigma_a/2)^2=3\,1/4$ より $3/4$ である。
:::
::::

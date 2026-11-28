---
title: "非可換ゲージ理論の量子化"
description: "非可換ゲージの固定と BRST 恒等式を構成し、背景場の一ループ発散からベータ関数を得る。"
notation:
  - symbol: 'c,\bar c,B'
    meaning: ゴースト、反ゴースト、補助場
  - symbol: 's'
    meaning: BRST 微分
  - symbol: 'b_0'
    meaning: ゲージ結合の一ループ係数
---

## この章の目標

- ゲージ軌道の Jacobian をゴースト積分で表す。
- BRST 微分と物理状態の条件を構成する。
- 背景場の二次演算子から漸近的自由性を導く。

## Faddeev–Popov の方法

### 軌道の交点と Jacobian

ゲージ場を積分するとき、ゲージ変換で結ばれた場は同じ物理配置を表す。場の積分からこの重複を除くため、各軌道が $\partial^\mu A^a_\mu=f^a$ と一度だけ交わる摂動的な近傍を取る。:eqref[eq-qft-ym-gauge] を微分すると $\delta A_\mu^a=(D_\mu\omega)^a$ だから、交点でのゲージ条件の軌道微分は $M^{ab}=\partial^\mu D_\mu^{ab}$。この微分の行列式を用いて正しい切片測度を得るのが :term[Faddeev–Popov 法]{en="Faddeev–Popov method" yomi="ふぁでーえふぽぽふほう"} である。

変数を軌道に沿って移し、軌道の体積を分母で除く。$f$ を Gaussian 平均すると $-(\partial\cdot A)^2/(2\xi)$ が加わる。Abel では $M=\Box$ なので行列式が場に依存しない。非可換では $D_\mu$ に $A_\mu$ が含まれるため、行列式が量子補正を作る。ゲージ条件の複数交点が現れる大域的な問題は、この摂動的な構成の仮定の外にある。

### 有限次元の軌道で Jacobian を調べる

平面の回転に対して $\chi(x,y)=y=0$ を選ぶと、半径 $r>0$ の軌道は $(r,0)$ と $(-r,0)$ の二箇所で交わる。回転角に関する $\chi$ の微分は各点で $+r,-r$。$\delta(\chi)$ の変数変換で絶対値 $r$ を掛けると全円周で積分値は $2$、符号つきで掛けると $0$ になる。これに対し $x>0$ の半軸を選ぶ近傍では交点は一つで積分値は $1$。局所的な Faddeev–Popov の式ではこの一交点条件を用いる。場の空間でも、軌道の複数交点や $M$ の零固有値を含む領域には同じ局所式をそのまま使えない。

### ゲージ体積を除いた積分の導出

経路積分の変数変換を使う前に、測度を有限個のモードに制限した模型として考える。境界で恒等変換となるゲージ変換を用いるか、残留する定数変換を別に除く。そうすれば自由な $\Box$（Euclid では $\partial^2$）の定数零モードを行列式へ誤って含めずに済む。Minkowski での自由逆演算子は真空境界条件を指定して定める。以下ではまず Euclid で軌道方向の Jacobian を考え、摂動則を求める段階で Minkowski に戻る。

:::definition[ゲージ切片と軌道方向の微分]
ゲージ条件を $\chi^a[A]=f^a$ とする。この方程式を満たす場の集合を局所的な :term[ゲージ切片]{en="gauge slice" yomi="げーじせっぺん"} と呼ぶ。軌道の接方向 $D_\mu\omega$ に沿って微分した
$$
(M[A]\omega)^a(x)=
\int\dd^dy\,\frac{\delta\chi^a[A](x)}{\delta A_\mu^b(y)}
(D_\mu\omega)^b(y)
$$
が :term[Faddeev–Popov 演算子]{en="Faddeev–Popov operator" yomi="ふぁでーえふぽぽふえんざんし"} である。ゲージ固定された積分で評価する $A$ は切片上にある。
:::

同じ軌道を記述する二つの場で行列式の値を無条件に等しいと置くと、軌道体積を正しく除けない。軌道上の各点を始点とし、
$$
\Delta_f[A]^{-1}=\int\mathcal DU\,\delta(\chi[A^U]-f)
$$
と定義する。変換を合成しても Haar 測度が不変なので $\Delta_f[A^V]=\Delta_f[A]$。切片上の $A$ では、積分を恒等変換の近くの座標 $\omega$ で評価できる。このとき $y=\chi[A^{e^{ig\omega}}]-f$ の Jacobian が $M[A]$ であり、デルタ関数が $y=0$ でその逆絶対値を抜き出す。ゆえに切片で $\Delta_f[A]=|\det M[A]|$。切片外で同じ式を使うなら $M$ は交点で評価する。

ゲージ不変な挿入 $O[A]$ の非規格化積分に $1=\Delta_f\int\mathcal DU\,\delta(\chi[A^U]-f)$ を入れる。$A\mapsto A^{U^{-1}}$ により $O,S,\Delta_f,\mathcal DA$ が不変なら、$U$ 積分は共通な体積因子になる。分子と分母で約して
$$
\langle O\rangle_f=
\frac{\int\mathcal DA\,e^{iS}\,O\,\delta(\chi[A]-f)\,|\det M[A]|}
{\int\mathcal DA\,e^{iS}\,\delta(\chi[A]-f)\,|\det M[A]|}
$$
を得る。同じ変数変換を規格化前に用いれば、分母のゲージ固定積分は $f$ を各切片が存在する範囲で変えても一定である。従って分子と分母を同じ Gaussian $\exp[-i\int f^2/(2\xi)]$ で積分してから比を取れる。通常の摂動的な線形ゲージではこの操作を自由 Gaussian の処方として定義する。Euclid では $\exp[-\int f^2/(2\xi)]$ を使えば $\xi>0$ で通常の積分になる。

:::example[非線形な切片では場のない行列式にならない]
一成分の Abel 接続に $\chi[A]=\partial\cdot A+\kappa A_\mu A^\mu$ を試す。$\delta A_\mu=\partial_\mu\omega$ により $M\omega=\Box\omega+2\kappa A^\mu\partial_\mu\omega$。従って Abel でゴーストが離れるのは線形条件の性質である。$\kappa=0$ なら対応するゴースト頂点が消える。正準次元は $[A]=1$ なので $d=4$ では $[\kappa]=0$ となる。ここではこの切片が摂動真空の近くで使える範囲を想定した。
:::

### Gribov コピーと軌道の零モード

:term[Gribov コピー]{en="Gribov copy" yomi="ぐりぼふこぴー"} は、同じ条件を満たすゲージ同値な異なる場である。コピー同士が離れているときは、片方の $M$ に零固有値があるとは限らない。これに対して $M\omega=0$ なら軌道方向で切片に接し、逆関数定理で交点を一意に構成できない。例えば $\chi(\theta)=r\sin\theta$ で $f$ を動かすと、$|f|<r$ には二つの交点があり $|f|=r$ で一つに合流する。そこで $\chi'(\theta)=0$ となる。有限の Euclid 領域で Landau 条件を選び、定数変換を除いた上で $-\partial_\mu D_\mu$ の最小固有値が零になる境界を考えると、この接触の場の例が得られる。全てのコピーを除いた非摂動的な測度の構成は、局所式からは決まらない。

### 軌道を割った測度を極座標で求める

先の回転の例を、正則な点での積分まで計算する。群の Haar 測度を $\dd\theta$ とし、観測量を半径だけの $f(r)$ とする。Euclid の重みを $e^{-ar^2}$、$a>0$ として、回転の体積 $2\pi$ で割ると
$$
\frac1{2\pi}\int_{\R^2}\dd x\,\dd y\,e^{-a(x^2+y^2)} f(r)
=\int_0^\infty r\,\dd r\,e^{-ar^2}f(r).
$$
同じものを切片 $y=0,x>0$ で積分すると、$\delta(y)$ と正の Jacobian $x$ により右辺の $r\,\dd r$ が出る。例えば $f=1$ の積分は $1/(2a)$、$f=r^2$ では $1/(2a^2)$、規格化した期待値は $1/a$。切片を選んだだけで $x$ を掛け忘れると、後者の期待値は $1/(2a)$ になってしまう。$x=0$ は軌道の次元が下がる固定点で、そこの Jacobian は零である。その一点はこの積分で測度零だが、場の理論で固定点や零モードが常に無視できるとは限らない。

### 残留変換と非共変ゲージ

軸型の条件 $n^\mu A_\mu=0$ を取り、一定の $n^\mu$ を用いる。切片では $n\cdot D\omega=n\cdot\partial\omega$ なのでゴースト行列式は形式的に場に依存しない。ただし $n\cdot\partial\omega=0$ となる変換はこの条件で残る。例えば $n$ が時間方向なら、時間によらない全ての空間関数が零モードである。線の端で変換を固定するなど、零モードの除去と逆演算子の境界条件が必要になる。運動量空間で $1/(n\cdot k)$ に遭遇したとき、どの境界条件で逆を取るかは切片の定義の一部である。共変な Gaussian ゲージを選ぶと離れたゴーストの相互作用は復活するが、ローレンツ共変な伝播核と次元正則化を同時に使いやすい。

## ゴースト場とその Feynman 則

Grassmann 積分で $\det M$ を表し、Minkowski で $\mathcal L_{\rm gh}=-\bar c^a\partial^\mu(D_\mu c)^a=(\partial^\mu\bar c^a)(D_\mu c)^a$ と取る。自由項の逆は $\langle T c^a\bar c^b\rangle=i\delta^{ab}/(p^2+i0)$。相互作用は $g f^{abc}(\partial^\mu\bar c^a)A_\mu^b c^c$ である。全ての運動量を入射として $e^{-ipx}$ を使えば、反ゴーストの入射運動量 $p$ の頂点係数は $g f^{abc}p^\mu$ となる。Grassmann 積の順序を反ゴースト、ゴーストの順で固定する。閉ゴースト線にはフェルミオンの縮約と同じ負号が付く。ゴーストは Lorentz スカラーだが、外部の物理状態には入れない。

### 行列式から閉ループを復元する

$M=M_0+V[A]$ とし、自由行列式を割る。有限モードで $\log(\det M/\det M_0)=\Tr\log(1+M_0^{-1}V)=\Tr X-\Tr X^2/2+\cdots$、$X=M_0^{-1}V$。$\Tr X^n$ では $n$ 個の挿入を伝播核が巡回してつなぐ。例えば二色だけを残した模型で $X=\begin{pmatrix}0&x\\y&0\end{pmatrix}$ とすると一次項は零、二次は $-xy$。実際、$\det(1+X)=1-xy$ の対数の先頭と一致する。ボースの複素場では逆行列式になるため、この対数全体の符号が逆である。

### Grassmann 行列式の規格化とゴーストの向き

有限個の対 $(c_i,\bar c_i)$ について積分順と符号を $\int\prod_i\dd\bar c_i\dd c_i\,e^{-\bar c c}=1$ となるよう決める。$e^{-\bar c Kc}$ を展開すると、積分で残る次数は各 $\bar c_i,c_i$ を一度ずつ含む項である。$c$ の順を並べ替えると置換の符号が付くので、残る係数は $\det K$。さらに独立な奇の源を入れて平方完成すれば
$$
\int\mathcal D\bar c\mathcal Dc\,
e^{-\bar c Kc+\bar\eta c+\bar c\eta}
=(\det K)e^{\bar\eta K^{-1}\eta}.
$$
源は書いた順で保つ。これにより Euclid の $\langle c_i\bar c_j\rangle=(K^{-1})_{ij}$ となる。反対順の相関では符号が反転する。図の矢印はこの二つの端の区別に使う。実スカラーの積分では二次形式の逆平方根、複素スカラーでは逆行列式が出るので、ゴーストでは複素スカラーとは逆符号の閉ループを得る。これは伝播粒子の統計を物理的な外部状態に割り当てる議論とは独立である。

:::example[一次および二次の背景頂点]
運動量空間で $K=K_0+V$ とし、$V_{p',p}$ は背景の運動量 $p'-p$ を受け取るとする。$K^{-1}=G_0-G_0VG_0+\cdots$ なので、ゴースト二端の先頭補正は $-G_0(p')V_{p',p}G_0(p)$。これに対して背景だけの有効作用は $-\Tr\log K$ であり、二次は
$$
\tfrac12\sum_{p',p}\tr_{\rm col}\bigl[G_0(p')V_{p',p}G_0(p)V_{p,p'}\bigr]
$$
になる。色の一頂点は随伴生成子一つに比例し $\tr_{\rm col}T^a_{\rm adj}=0$ なので、真空の背景一脚には寄与しない。二次では $\tr_{\rm col}T^a_{\rm adj}T^b_{\rm adj}=C_A\delta^{ab}$ が残る。二頂点を外場で微分すると二つの配置が生じ、表示の $1/2$ が消える。
:::

### 非零の $\xi$ で使う自由ゲージ核

二次作用を部分積分し、ゲージ固定を含む積分核を $K_{\mu\nu}(k)=-k^2P^T_{\mu\nu}-(k^2/\xi)P^L_{\mu\nu}$ とする。ただし $P^L_{\mu\nu}=k_\mu k_\nu/k^2$、$P^T_{\mu\nu}=\eta_{\mu\nu}-P^L_{\mu\nu}$ である。直交した射影なので $iK^{-1}$ は
$$
D_{\mu\nu}^{ab}(k)=\frac{-i\delta^{ab}}{k^2+i0}
\left[\eta_{\mu\nu}-(1-\xi)\frac{k_\mu k_\nu}{k^2+i0}\right]
$$
となる。式は質量ゼロ極への同じ真空処方を前提とする。$\xi=1$ では計量だけ、$\xi\to0$ では横射影だけが残る。極の縦方向もループに現れるため、内部線を二つの物理偏極の和だけに置き換えることは、共変ゲージのこの Feynman 則とは両立しない。物理状態の選別は、ゲージ場とゴーストを共に扱う次節で行う。

### ループの色因子と外部運動量

ゴースト線に外部ゲージ場を二本結合する例を取る。入射外運動量を $p,-p$、内部線の運動量を $\ell,\ell+p$ とすれば、それぞれの反ゴースト端は $\ell$ と $\ell+p$ の成分を頂点へ渡す。随伴の二つの行列トレースは $C_A\delta^{ab}$。よって、全二点関数の符号を除いたテンソル積分は
$$
I_{\mu\nu}=\int\frac{\dd^d\ell}{(2\pi)^d}
\frac{\ell_\mu(\ell+p)_\nu}{\ell^2(\ell+p)^2}.
$$
Euclid の $p^2>0$ とする。Feynman パラメータ $x$ により $r=\ell+xp$ を取れば、偶部分は $r_\mu r_\nu-x(1-x)p_\mu p_\nu$。積分では $r_\mu r_\nu\mapsto\delta_{\mu\nu}r^2/d$ を使う。$\Delta=x(1-x)p^2$ とし、$I_2$ の極を $1/(16\pi^2\bar\epsilon)$ とすると $\int r^2/(r^2+\Delta)^2$ の極は $-2\Delta I_2|_{\rm pole}$。これは $r^2=(r^2+\Delta)-\Delta$ と :eqref[eq-qft-loop-master] の一次分母の極を使った。したがって
$$
I_{\mu\nu}\big|_{\rm pole}
=-\frac{p^2\delta_{\mu\nu}/12+p_\mu p_\nu/6}{16\pi^2\bar\epsilon}.
$$
$\int_0^1x(1-x)\dd x=1/6$ を用いた。このテンソルだけを $p^\mu$ で縮約すると零にならない。横の恒等式は、残りのゲージ図と同じ規約で合算した二点に対して検査する。

:::example[対称因子を行列式から読む]
複素のゴーストは $c$ から $\bar c$ へ向かう線を持つ。二背景挿入の $-\Tr\log K$ には $1/2$ があるが、二つの外場で微分すると二配置があるので係数は一。これに対し実ゲージ場の $\tfrac12\Tr\log \Delta$ では、実 Gaussian の $1/2$ がその後も残る。同じ二本の内部線でも、数える場の積分が違うため、図の係数を同じにすることはできない。
:::

## BRST 対称性と物理的状態

### 奇の対称変換を作る

:term[BRST 対称性]{en="BRST symmetry" yomi="びーあーるえすてぃーたいしょうせい"} の微分 $s$ は奇の導分で、$s(XY)=(sX)Y+(-1)^{|X|}XsY$ である。行列記法で
$$
sA_\mu=D_\mu c,\quad sc=ig c^2,\quad s\bar c=B,\quad sB=0,\quad s\psi=ig c\psi .
\label{eq-qft-brst}
$$
$sc^a=-gf^{abc}c^bc^c/2$ は奇の積から得る。$s^2c=ig[(sc)c-c(sc)]=0$、$s^2\psi=ig(sc)\psi-igc(igc\psi)=0$。また $s^2A=D(sc)-ig[(Dc)c+c(Dc)]=0$ である。ゲージ固定とゴースト項は $s[\bar c^a(\partial\cdot A^a+\xi B^a/2)]$ であり、展開すると $B\partial\cdot A+\xi B^2/2-\bar c\partial\cdot Dc$。$B$ の方程式を代入すれば元のゲージ固定を与える。$s^2=0$ なのでこの項は BRST 不変である。

量子論で $s$ を生成する保存電荷 $Q_B$ が冪零であれば、摂動論の物理状態は $\ker Q_B/\operatorname{im}Q_B$ のゴースト数零の部分で定義される。閉状態 $\ket\alpha$ に対し $Q_B$ が適切な随伴性を持てば $\langle\alpha|Q_B|\beta\rangle=0$ なので、BRST 完全な状態を足しても物理行列要素は変わらない。自由理論では $sA_\mu=\partial_\mu c$、$s\bar c=B$ が縦方向とゴーストを組にする。組 $(u,v)$ で $su=v,sv=0$ なら $h(v)=u,h(u)=0$ と拡張した $h$ に対し $\{s,h\}$ は組の変数の総次数である。次数 $n>0$ の閉多項式 $X$ は $X=s(hX/n)$ なのでコホモロジーから消える。残る自由粒子は横偏極である。

### ゲージ固定に依存しない行列要素の条件

ゲージ固定とゴースト作用を $s\Psi_\xi$ と書くと、$\partial_\xi S=s\,\partial_\xi\Psi_\xi$ である。:eqref[eq-qft-brst] では $\partial_\xi\Psi_\xi=\frac12\int\bar c^a B^a$。物理挿入 $O$ が $sO=0$ で、測度と積分領域が BRST を保つなら、変数変換で $\langle s X\rangle=0$ である。従って規格化した期待値で
$$
\partial_\xi\langle O\rangle
=i\langle O\,s(\partial_\xi\Psi_\xi)\rangle_c=0
$$
となる。実際、分母由来の真空挿入は $\langle s(\partial_\xi\Psi_\xi)\rangle=0$、分子は $O$ が偶なら $\langle s(O\,\partial_\xi\Psi_\xi)\rangle=0$ である。例えば閉じた Wilson ループは端点で共役を受けるだけでトレースが不変なので閉挿入になる。外部時空点のゲージポテンシャル単独は $sA=Dc$ を持つので、同じ手順からその相関の $\xi$ 独立性は出ない。

### ゴースト数による閉状態の区別

$A,B$ のゴースト数を $0$、$c,\bar c$ を $+1,-1$ とする。:eqref[eq-qft-brst] で $s$ の作用は常にゴースト数を一つ上げ、$\Psi$ のゴースト数が $-1$ なら $s\Psi$ のゴースト数は零である。自由理論の二重項 $su=v,sv=0$ で例えば $u$ を偶、$v$ を奇とする。$s(u^n)=n u^{n-1}v$、$s(u^{n-1}v)=0$ なので、$u,v$ を含む閉多項式は次数零を除いて全て完全である。例えば $uv=s(u^2/2)$。この消去では奇変数の $v^2=0$ が使われ、通常の数で $v$ を置き換えると同じ複体にはならない。

### 自由一粒子空間の二重項

運動方程式を満たす一つの光的運動量 $k$ を取る。Feynman ゲージの自由な解では $sA_\mu\propto k_\mu c$。$k\cdot n=1$ を満たす補助ベクトル $n$ で縮約すると $s(n\cdot A)\propto c$ となる。補助場を消した自由系では $B=-\partial\cdot A$ なので $s\bar c\propto k\cdot A$、そして $s(k\cdot A)\propto k^2c=0$。二つのモードの係数を非零定数で規格化して $(n\cdot A,c)$ と $(\bar c,k\cdot A)$ の二重項を得る。$k,n$ に直交する偏極を二つ取ればその $s$ は零である。例えば $k=(E,0,0,E)$、$n=(1,0,0,-1)/(2E)$ に対し、$x,y$ 方向の空間偏極が残る。

このモード計算を自由 Fock 空間に広げるには、生成演算子から作る有限粒子ベクトルに二重項を数える作用素 $N_{\rm dbl}$ を作用させる。$h$ を逆向きの奇写像としてすべての生成変数に導分で延長すると $\{s,h\}=N_{\rm dbl}$。核にあるベクトルのうち $N_{\rm dbl}$ の固有値 $n>0$ の成分は $s(h/n)$ で書ける。$N_{\rm dbl}=0$ の成分は横偏極と物質の自由粒子から作る。共変状態空間には不定内積があるが、この自由な余剰モードを商で除いた内積は正である。相互作用のある摂動論での主張には保存された BRST 電荷と漸近状態の構成を使う。色のある自由な漸近状態が閉じ込め相の厳密な粒子であることは、この自由計算からは決まらない。

### 行列積で物質と曲率の BRST 変換を追う

$\bar\psi$ も奇変数として扱う。$\psi\mapsto U\psi$ と $\bar\psi\mapsto\bar\psi U^{-1}$ のゲージ変換に奇の定数パラメータ $\eta$ を左に置き $U=1+ig\eta c$ とすると $s\bar\psi=ig\bar\psi c$。この符号では $s(\bar\psi\psi)=(s\bar\psi)\psi-\bar\psi s\psi=0$。また $s(D_\mu\psi)=igcD_\mu\psi$ は、$s$ を $D_\mu$ と $\psi$ の両方に当てたときの $D_\mu c$ が打ち消すことから出る。従って Dirac の運動項も $s$ で閉じる。曲率には $sF_{\mu\nu}=D_\mu D_\nu c-D_\nu D_\mu c=-ig[F_{\mu\nu},c]$ を得るので、$\tr F^2$ は巡回性によって不変である。このように古典作用の不変性を示した後、固定項の不変性は $s^2=0$ から従う。

:::example[ゲージ条件を変えたときの挿入]
$\Psi_t=\int\bar c^a(\chi_t^a+\xi B^a/2)$ とし $\chi_t^a=\partial\cdot A^a+t\,n\cdot A^a$ とする。境界条件と正則化が全ての $t$ で BRST を保つ範囲で $\partial_tS=s\int\bar c^a n\cdot A^a$。二つの分離した $\tr F^2$ を挿入した積は閉じるので、同じ規格化でくりこんだ相関の $t$ 微分は零である。これは $t$ と独立に定義した演算子で、接触点を除いた主張である。$t$ を変えて零モードや特異な境界条件を通過するときは、この連続変形だけで同じ積分になるとは言えない。
:::

### 自由モードの行列要素と商に入る状態

$v=Q_B u$ と書ける状態 $v$ の内積を考える。$Q_B$ が不定内積に関する自己随伴（あるいは規約上の位相付き自己随伴）なら、$w\in\ker Q_B$ に対して $(w,v)=(Q_Bw,u)$ に比例し零である。$v$ 自身も核にあるので $v$ のノルムは零。従って商空間の上に閉状態の内積を誘導できる。単にノルム零の任意のベクトルを取り除けばよいという意味ではない。例えば時刻偏極と縦偏極の和と差はそれぞれ零ノルムにできるが、両者の混合内積は一般に非零。$s(n\cdot A)$ はゴーストに写るため、その組の片方はそもそも閉状態にならない。

相互作用描像で散乱作用素が漸近 BRST 電荷と可換であれば $S(\ker Q_B)\subset\ker Q_B$、$S(\operatorname{im}Q_B)\subset\operatorname{im}Q_B$。後者は $S Q_Bu=Q_B S u$ を使った。従って $S$ は物理状態の商へ降り、元の不定内積での $S$ の保存性から商での内積も保つ。漸近場の物理的な空間が自由横偏極で構成される摂動論では、その商でユニタリな条件を課せる。この議論の前提は保存する BRST 対称性である。ゲージアノマリーで測度が変化する場合は、同じ変数変換を使った恒等式が破れる。一般の異常の分類は [アノマリー](/physics/quantum-field-theory/anomalies/) で扱う。

### 非線形ゲージ条件にも使う冪零性

固定関数 $\chi$ に場依存性があっても、$\Psi=\int\bar c^a(\chi^a+\xi B^a/2)$ はゴースト数 $-1$。$s\Psi=B^a\chi^a+\xi B^2/2-\bar c^a\,s\chi^a$。$s^2\chi=0$ により $s(B\chi)=B\,s\chi$ と $s(-\bar c\,s\chi)=-B\,s\chi$ は相殺する。例えば一成分 Abel の $\chi=\partial\cdot A+\kappa A^2$ では $sc=0$ でも $s\chi=\Box c+2\kappa A^\mu\partial_\mu c$。$s$ をもう一度作用させると $\partial^\mu c\,\partial_\mu c$ が出るが、奇変数の交換に対して計量が対称なので零。非線形ゲージのゴースト相互作用と、冪零な変換とが矛盾しないことを直接確認できる。

## Slavnov–Taylor 恒等式

場ごとに通常の源 $J_i$ と、非線形変換 $s\Phi_i$ の源 $K_i$ を入れる。規格化した積分の BRST 変数変換で作用と測度が不変とすれば $\sum_i\langle s\Phi_i\rangle J_i=0$（源項の順を $\Phi_iJ_i$ とする） となる。$K$ の微分で $\langle s\Phi_i\rangle$ を表し、$J$ について Legendre 変換すると、純ゲージ理論では左・右微分を整合させて
$$
\int\dd^4x\,\left(
\frac{\delta\Gamma}{\delta A_\mu^a}\frac{\delta\Gamma}{\delta K_A^{a\mu}}
+\frac{\delta_l\Gamma}{\delta c^a}\frac{\delta_l\Gamma}{\delta K_c^a}
+B^a\frac{\delta_l\Gamma}{\delta\bar c^a}
\right)=0
\label{eq-qft-slavnov-taylor}
$$
となる。これが :term[Slavnov–Taylor 恒等式]{en="Slavnov–Taylor identity" yomi="すらぶのふていらーこうとうしき"} である。物質があれば同様の対を加える。通常のゲージ変換を頂点ごとに課すだけではゴーストの項を失う。例えばくりこみ後も一つの $g$ を使うと、三ゲージ頂点は $Z_gZ_A^{3/2}$、ゴースト頂点は $Z_g Z_c Z_A^{1/2}$ なので両者を独立には調整できない。

### 補助場の恒等式と縦伝播

$B$ を積分する前にその変数の平行移動を行うと、線形ゲージ条件では有効作用にも $\delta\Gamma/\delta B^a=\partial\cdot A^a+\xi B^a$ が成り立つ。右辺が場に一次であり、その期待値を古典場の値へ置けるためである。また $B=s\bar c$ から、BRST 不変真空で $\langle BB\rangle=\langle s(\bar c B)\rangle=0$。Gaussian 積分で $B=B'-\partial\cdot A/\xi$ と分けると $B'$ は局所核 $i/\xi$ を持つ。従って Fourier 空間の完全伝播関数の二重縮約は $k_\mu k_\nu D^{\mu\nu}=-i\xi$ である。相互作用で変わるスカラー関数を横部分に置く理由がここにある。この議論は接触項 $i/\xi$ を含めた時間順序の規約を用いている。

### Legendre 変換と線形化

右微分と左微分の使い分けには、変分の位置を定めると便利である。ここで $\delta f=(\delta_r f/\delta z)\delta z=\delta z(\delta_l f/\delta z)$ とする。左に置いた奇定数 $\eta$ で $\delta z=\eta sz$ と書けば、偶の $f$ の変分は $\delta f=\eta\sum_i(sz_i)(\delta_l f/\delta z_i)$。また $sA$ は奇、$sc$ は偶である。非線形変換の源項を $K_A sA+K_c sc$ と書くと、$\delta_l\Gamma/\delta K_A=\langle sA\rangle_J$、$\delta_l\Gamma/\delta K_c=\langle sc\rangle_J$ が古典場を固定したまま成立する。$\Gamma=W-J\Phi$ の微分では、$J$ の変分から来る $W$ と積 $J\Phi$ の寄与が等しく消える。

運動方程式に関係する $J$ の因子を $\Gamma$ の場微分に置き換えると、$c$ の項には $\delta_l\Gamma/\delta c$ が現れる。$c$ の BRST 変換は偶なので、:eqref[eq-qft-slavnov-taylor] の二つの因子はこの項では交換できる。$\bar c$ の変換は $B$ と線形なので、非線形源を入れなくても $\langle B\rangle_J$ はそのまま古典場となる。$B$ の恒等式はゲージを固定した二点の縦成分を制約する一方、ST 恒等式の非線形部分は異なる頂点を結び付ける。

$S$ を非線形源込みの古典作用とし、$\Gamma=S+\hbar\Gamma_1+\cdots$ とループ次数を形式的に記録する。Slavnov 汎関数を $\mathcal S$ と書けば一次は
$$
\mathcal B_S\Gamma_1=0,\qquad
\mathcal B_S X=\left.\frac{\dd}{\dd t}\mathcal S(S+tX)\right|_{t=0}
$$
である。この式には $S$ の場微分と $X$ の源微分の積と、$S$ の源微分と $X$ の場微分の積の両方がある。線形ゲージの $B\,\delta/\delta\bar c$ も残る。例えば場と源を逆の因子で再規格化すると、各場・源対の積の規格化が打ち消し合う。ゲージ結合は三ゲージ、四ゲージ、ゴースト、物質のいずれにも同じ値を割り当てるため、この再規格化から頂点ごとの因子の関係が得られる。

:::example[一つの結合と四つの頂点]
$Z_3,Z_4,\widetilde Z_1,Z_{1\psi}$ をそれぞれ三ゲージ、四ゲージ、ゴースト–ゲージ、Dirac–ゲージ頂点に掛ける係数とする。場の因子を $Z_A,Z_c,Z_\psi$ としたとき
$$
\frac{Z_3}{Z_A^{3/2}}=
\frac{\sqrt{Z_4}}{Z_A}=
\frac{\widetilde Z_1}{Z_c Z_A^{1/2}}=
\frac{Z_{1\psi}}{Z_\psi Z_A^{1/2}}=Z_g .
$$
平方根は木の値が $1$ の枝である。例えば $Z_i=1+\delta Z_i$ とした一次では $\delta Z_4=2\delta Z_3-\delta Z_A$。非可換理論の $\widetilde Z_1=1$ を全てのゲージ・全ての有限運動量で仮定する必要はない。
:::

### 反ゴースト恒等式

線形共変ゲージでは、$\bar c$ は $\partial_\mu\bar c\,sA^\mu$、非線形源は $K_A^\mu sA_\mu$ として入る。この二項は $K_A^\mu+\partial^\mu\bar c$ によってのみ依存する。通常の場の源を Legendre 変換した後でも、両者を逆向きに変える恒等式は保たれる。例えば $\delta\bar c=\varepsilon(x)$、$\delta K_A^\mu=-\partial^\mu\varepsilon$ と置く。古典源付き積分の $\bar c$ の平行移動で普通の反ゴースト源の線形項だけが出て、これが Legendre 変換で移る。境界が消えると
$$
\frac{\delta_l\Gamma}{\delta\bar c^a}
+\partial^\mu\frac{\delta_l\Gamma}{\delta K_A^{a\mu}}=0 .
$$
ここでは $\delta\bar c$ を左に置いて部分積分した。この式は線形ゲージの作用に依存し、任意のゲージ固定で同じ形を仮定できない。

$b^a=\delta\Gamma/\delta B^a=\partial\cdot A^a+\xi B^a$ と組み合わせると、補助場とゴーストを含む縦の二点を切り分けられる。例えば接触項を含めた純ゲージ二点について $D_{\rm full}=D_T P^T-i\xi P^L/(k^2+i0)$ と書ける。この形で $\xi$ を裸量で書くかくりこんだ量で書くかは、$A$ の規格化と同時に決める。$A_0=\sqrt{Z_A}A$ とすると $D_0=Z_A D$、縦の縮約から $\xi_0=Z_A\xi$。同じ関係は $B_0=B/\sqrt{Z_A}$ を選び $B_0\partial\cdot A_0$ を保つと $\xi_0 B_0^2/2=\xi B^2/2$ からも得られる。従って $\xi$ の走りは、一般の共変ゲージの $\beta_g$ と同じ独立な結合の走りではない。

## 1 ループ発散とくりこみ

### 背景曲率に比例する発散を選ぶ

背景を $\bar A$、量子揺らぎを $a$ とし Euclid へ移る。背景ゲージの $\xi=1$ では二次演算子が
$$
(\Delta_1)_{\mu\nu}=-\bar D^2\delta_{\mu\nu}-2\Omega_{\mu\nu},\qquad
\Delta_0=-\bar D^2,\qquad
\Omega_{\mu\nu}=[\bar D_\mu,\bar D_\nu]
$$
となる。曲率を $\bar F+\bar D_\mu a_\nu-\bar D_\nu a_\mu-ig[a_\mu,a_\nu]$ と展開する。一次曲率の二乗の非対角二微分はゲージ固定で消え、微分を交換した項と背景曲率の交差項がそれぞれ $-\Omega_{\mu\nu}$ を与える。従って一ループ作用は $\frac12\Tr\log\Delta_1-\Tr\log\Delta_0$。Dirac 物質は位相を無視できるベクトル表現で $-\frac12\Tr\log[-(\not {\bar D})^2+m^2]$ を加える。伝播関数ごとに非共変な頂点を集める代わりに、背景の局所な $\bar F^2$ を計算すれば結合の極が求まる。

### 背景曲率の二次項を展開する

$\bar D$ について部分積分し、Euclid の二次作用を取る。線形化した曲率の二乗は $\frac12\int[(\bar D_\mu a_\nu)^2-(\bar D_\mu a_\nu)(\bar D_\nu a_\mu)]$。ここに $\frac12\int(\bar D\cdot a)^2$ を足すと、異なる順序の微分の差が $\frac12\int a_\nu[\bar D_\mu,\bar D_\nu]a_\mu=-\frac12\int a_\mu\Omega_{\mu\nu}a_\nu$ となる。さらに $\bar F$ と $-ig[a_\mu,a_\nu]$ の交差項も随伴内積の不変性で $-\frac12\int a_\mu\Omega_{\mu\nu}a_\nu$。二つを合わせ、$S^{(2)}=\frac12\int a_\mu(\Delta_1)_{\mu\nu}a_\nu$ と書くと $-2\Omega_{\mu\nu}$ が得られる。例えば可換ゲージ場では随伴作用が自明なので $\Omega=0$。背景との交差項はなく、光子自身にはこの一ループ補正が生じない。

### 一ループ対項の局所性と背景の不変性

連続する各頂点で大きいループ運動量 $\ell$ に外部運動量 $p$ が足される。例えば $(\ell+p)^{-2}=\ell^{-2}-(2\ell\cdot p+p^2)\ell^{-4}+\cdots$。高い運動量を切り分けた領域で外部運動量と質量を展開すると、紫外発散するのは表面的次数までの多項式である。低運動量の特異性をこの Taylor 展開に混ぜないため、外部 Euclid 運動量を非例外点に取るか、熱核では赤外側に減衰因子を入れる。一ループ 1PI では真のループの部分発散がないので、この操作から局所な対項を取れる。

ゲージ場とゴーストに次元一、Dirac 場に $3/2$、$B$ に二を与える。純ゲージなら三ゲージ頂点とゴースト頂点の微分は一、四ゲージは零である。内部線を $I_A,I_c,I_\psi$、微分付き頂点の総数を $V_\partial$ とすると
$$
\omega=4L-2I_A-2I_c-I_\psi+V_\partial .
$$
各頂点の次元和が四であるため、内部線の両端で数えた次元を消すと $\omega=4-E_A-E_c-\tfrac32 E_\psi$ になる。$E$ は該当する外脚を数える。例えばゲージ五点の全体の $\omega=-1$ なので、一ループで独立なゲージ五点の発散はない。ゴースト二点の $\omega=2$ に対し $\bar c$ に微分が付くことは運動量の依存をさらに制約する。

背景法では $\bar A$ の変換を接続として、$a,c,\bar c,B$ は随伴表現で斉次に変換する。$\bar D^\mu a_\mu$ は斉次に変換し、その二乗は不変である。背景測度も保たれる場合、背景だけの有効作用は通常のゲージ不変性を持つ。質量次元四以下で平坦空間の純ゲージ背景を使い、CP 偶の局所発散を取ると、場に依存する候補は $\tr\bar F^2$。例えば $\bar A_\mu^2$ は接続の非斉次項が消えないため候補から外れる。CP 奇の $\tr\bar F\widetilde{\bar F}$ は全微分であり、CP 偶の作用を保つ正則化ではこの計算に混ざらない。

### 背景と量子場のくりこみ

背景の二点の極は背景三点・四点の極も決める。背景共変微分を $\partial-ig\bar A$ に保つ規格化を採用すると、ゲージ変換パラメータを $U$ として書いた変換の非斉次項の係数は $1/g$ である。裸の場 $\bar A_0=Z_{\bar A}^{1/2}\bar A$ と裸結合 $g_0=\mu^\epsilon Z_g g$ を使った式を比較すると $Z_g Z_{\bar A}^{1/2}=1$。一方、$a$ は斉次変換なので同じ理由から $Z_g Z_a^{1/2}=1$ とは言えない。例えばゲージ固定の $\xi$ が現れる普通の量子二点と、$\beta_g$ の極を読む背景二点は規格化の対象が異なる。背景をゼロとした散乱のくりこみには通常の ST の関係を使う。

:::example[二点の係数から三点・四点を読む]
$\bar F=\bar F^{(1)}+g\bar F^{(2)}$ とし、背景を正準規格化する前の局所対項を $\delta Z_{\bar A}\bar F^2/4$ とする。その成分は $\delta Z_{\bar A}(\bar F^{(1)})^2/4$、$\delta Z_{\bar A}g\bar F^{(1)}\bar F^{(2)}/2$、$\delta Z_{\bar A}g^2(\bar F^{(2)})^2/4$。例えば二次だけを引くと、$\delta\bar A_\mu$ の交換子部分による二次項の変分を打ち消す三次項が残らない。背景ゲージ不変な極を消すには、背景の全曲率の平方を引く。
:::

### 普通の Feynman ゲージで二点の極を求める

背景を使わない二点の極には縦・ゴーストの寄与が個別に現れる。純ゲージ、$\xi=1$、外部 Euclid 運動量 $p^2>0$ とする。三ゲージ頂点の三本の入射を $p,\ell,-p-\ell$、外部添字を $\mu$、内部を $\alpha,\beta$ として、$g, f$ を除いたテンソルを
$$
V_{\mu\alpha\beta}
=\delta_{\mu\alpha}(p-\ell)_\beta+
 \delta_{\alpha\beta}(2\ell+p)_\mu+
 \delta_{\beta\mu}(-2p-\ell)_\alpha
$$
とする。Euclid の逆二点の自己エネルギーは、ゲージ泡とゴースト泡の和で
$$
-\delta^{ab}g^2 C_A
\int\frac{\dd^d\ell}{(2\pi)^d}
\frac{\tfrac12 V_{\mu\alpha\beta}V_{\nu\alpha\beta}-
\ell_\mu(\ell+p)_\nu}{\ell^2(\ell+p)^2}
$$
を得る。例えば三ゲージ頂点を二次で有効作用へ入れる符号は、実の Gaussian の $\tfrac12\Tr\log(K+V)$ の二次が $-\tfrac14\Tr K^{-1}V K^{-1}V$ であることから読める。頂点の二つの微分の $i$ と反対側での運動量反転を同時に取ると表示の二次テンソルになる。ゴーストは $-\Tr\log$ から相対的に逆符号。四点タドポールはこの質量なしの次元正則化では尺度を持たないので零である。

$r=\ell+xp$ として偶部分を取る。分子の組を四次元で計量縮約し、$r_\mu r_\nu\mapsto r^2\delta_{\mu\nu}/4$ とすると
$$
\begin{aligned}
N_{\mu\nu}
={}&\delta_{\mu\nu}\left[2r^2+\left(x^2-x+\frac52\right)p^2\right]\\
&+\left(4x^2-4x-1\right)p_\mu p_\nu .
\end{aligned}
$$
この計算では極に掛かる次数零の部分だけが要る。例えば $p$ を第四方向に取ると横の対角は第一行だけ、第四成分には二行目が加わる。$r^2$ の極積分は $-2x(1-x)p^2$ 倍の二分母の極なので、$\delta_{\mu\nu}p^2$ の係数は $\int_0^1(5x^2-5x+5/2)\dd x=5/3$。$p_\mu p_\nu$ の係数は $\int_0^1(4x^2-4x-1)\dd x=-5/3$。従って自己エネルギーの極は横であり、$Z_A$ に $+(5/3)C_A g^2/(16\pi^2\bar\epsilon)$ が必要になる。この $5/3$ は量子場の二点をくりこむ係数であって、背景の $11/3$ ではない。

同じゲージのゴースト逆二点の補正は $-g^2C_A\int p\cdot(\ell+p)/[\ell^2(\ell+p)^2]$。線形分子の奇な部分を除くと $(1-x)p^2$ となり、極は $-g^2C_Ap^2/(32\pi^2\bar\epsilon)$。従って $Z_c-1=+(1/2)C_A g^2/(16\pi^2\bar\epsilon)$。純ゲージの $Z_g-1=-(11/6)C_A g^2/(16\pi^2\bar\epsilon)$ を背景から入れると、ST 関係からゴースト頂点の極は $\widetilde Z_1-1=-(1/2)C_A g^2/(16\pi^2\bar\epsilon)$。ここではゲージ二点だけから $Z_g$ を推定せず、別の頂点の規格化と整合させた。

## 背景場法による β 関数の計算

### 短時間展開で係数を数える

平坦空間の $\Delta=-D^2+E$ の熱核で、発散に必要な局所係数は全微分を除き
$$
a_4=\tr\bigl(\Omega_{\mu\nu}\Omega_{\mu\nu}/12+E^2/2\bigr).
$$
$E^2/2$ は定数 $E$ を持つ自由熱核の $e^{-tE}$ の二次項で決まる。$\Omega^2$ は局所性と共変性からこの次数で許される残りの項であり、係数は平面の定数可換磁場で固定できる。荷電スカラーの横振動子の固有値は $(2n+1)|qB|$、単位面積あたりの縮退は $|qB|/(2\pi)$。これらは $A_y=Bx$ で $p_y$ ごとに中心 $p_y/(qB)$ の振動子となることから従う。自由熱核で割ると $t|qB|/\sinh(t|qB|)=1-(tqB)^2/6+\cdots$。$\Omega_{\mu\nu}\Omega_{\mu\nu}=-2q^2B^2$ なので係数は $1/12$ である。

ベクトルの $E_{\mu\nu}=-2\Omega_{\mu\nu}$ では Lorentz 添字のトレースで $E^2$ が $-4\Omega_{\mu\nu}\Omega_{\mu\nu}$ になる。従ってベクトルとゴーストの重みを含む係数は $\frac12(4/12-2)-1/12=-11/12$。随伴 $\tr\Omega^2=-g^2 C_A\bar F^2$ により $11g^2C_A\bar F^2/12$ である。Dirac の $E=-[\gamma^\mu,\gamma^\nu]\Omega_{\mu\nu}/4$ では四スピン成分と四行列トレースにより $4/12-1=-2/3$。重み $-1/2$ と $\tr_R\Omega^2=-g^2T(R)\bar F^2$ を掛けると $-g^2T(R)\bar F^2/3$ になる。

$\Tr\log\Delta=-\int_0^\infty(\dd t/t)\Tr e^{-t\Delta}+$ 定数とし、紫外の $\int_0\dd t\,t^{-1+\epsilon}$ から極を取る。背景二点の対項は $\delta Z_{\bar A}=b_0g^2/(16\pi^2\bar\epsilon)$、$b_0=11C_A/3-4\sum_{\rm Dirac}T(R)/3$ である。背景共変微分の $g\bar A$ を全体でくりこむので $Z_g Z_{\bar A}^{1/2}=1$。$g_0=\mu^\epsilon gZ_g$ を固定して微分すると
$$
\beta_g=-\epsilon g-\frac{b_0}{16\pi^2}g^3+\mathcal O(g^5).
\label{eq-qft-ym-beta}
$$
球対称四次元の熱核係数 $(4\pi t)^{-2}$ と、作用の $\bar F^2/4$ が $16\pi^2$ と係数の四を与える。

### 複素スカラーの寄与と多重度

複素スカラー $X$ を表現 $R$ に置くと、Euclid の二次積分は $(\det[-\bar D^2+m^2])^{-1}$ である。一ループ作用の重みは $+\Tr\log$ で、スピン曲率 $E$ は零。従って $a_4$ の重み付き寄与は $\tr_R\Omega^2/12=-g^2T(R)\bar F^2/12$。前節と同じ紫外極の換算で $b_0$ に $-T(R)/3$ を足す。実スカラーなら重みは $1/2$ なので $-T(R)/6$。例えば $SU(2)$ で基本複素スカラー一つは $-1/6$、随伴実スカラー一つは $-1/3$ である。色の次元だけでなく、曲率二乗のトレースの規格化と実成分の数を区別する。

### 背景二点から結合の極を決める微分

:eqref[eq-qft-ym-beta] の前の関係を一般の一ループ係数 $b_0$ について展開する。$Z_{\bar A}=1+b_0g^2/(16\pi^2\bar\epsilon)$ なら $Z_g=1-b_0g^2/(32\pi^2\bar\epsilon)$。固定した $g_0=\mu^\epsilon g Z_g$ で微分して
$$
0=\epsilon g\left(1-\frac{b_0g^2}{32\pi^2\bar\epsilon}\right)+
\beta_g\left(1-\frac{3b_0g^2}{32\pi^2\bar\epsilon}\right)+\cdots .
$$
$\beta_g=-\epsilon g+Bg^3$ としたとき、$g^3$ の係数には二つの有限寄与がある。明示的な $\epsilon$ の項は $-b_0/(32\pi^2)$、$\beta$ の先頭と括弧の積は $+3b_0/(32\pi^2)$。これらと $B$ の和を零にして $B=-b_0/(16\pi^2)$。従って $d=4$ に置くのは極と $\epsilon$ を掛けた後である。

:::example[Abel 極限を独立に確認する]
ゲージ場の随伴曲率を消し、単位電荷 Dirac 場一つだけを残すと $C_A=0,T(R)=1$。従って $b_0=-4/3$、$\beta_g=g^3/(12\pi^2)$。$g=e$、$\alpha=e^2/(4\pi)$ で $\beta_\alpha=2\alpha^2/(3\pi)$ となる。:eqref[eq-qft-vacuum-polarization] の大きい $Q$ の $\log Q^2$ の係数に $\alpha$ を掛け、$\log Q$ で微分した値と一致する。
:::

### 質量と調整尺度を入れた紫外極

熱核の局所係数は適切な Euclid 背景で短い固有時間の展開を行って求める。$\Tr$ は時空と内部自由度の両方、$\tr$ は内部自由度を表す。赤外を除く補助質量 $m_{\rm IR}>0$ を掛けると、曲率二乗に掛かる積分は
$$
\frac{\mu^{2\epsilon}}{(4\pi)^{2-\epsilon}}
\int_0^\infty\dd t\, t^{-1+\epsilon}e^{-m_{\rm IR}^2t}
=\frac{1}{16\pi^2}\left(\frac1{\bar\epsilon}+
\log\frac{\mu^2}{m_{\rm IR}^2}\right)+\mathcal O(\epsilon),
\quad
\frac1{\bar\epsilon}=\frac1\epsilon-\gamma_{\rm E}+\log4\pi .
$$
ここでは背景の大きい $t$ の積分をこの式で近似することは要らない。極の係数は $e^{-m_{\rm IR}^2t}=1+\mathcal O(t)$ で確定し、質量が変えるのは有限部と、質量を含む別の局所演算子の係数である。従って定数磁場におけるベクトル演算子の低い固有値の性質から、紫外極の符号を判断する必要はない。

Dirac の計算では $\{\gamma_\mu,\gamma_\nu\}=2\delta_{\mu\nu}$ と $\tr_{\rm spin}1=4$ を使う。反対称な $\Omega_{\mu\nu}$ と縮約すると $E=-\gamma_\mu\gamma_\nu\Omega_{\mu\nu}/2$。四行列のトレースは
$$
4(\delta_{\mu\nu}\delta_{\rho\sigma}
-\delta_{\mu\rho}\delta_{\nu\sigma}+\delta_{\mu\sigma}\delta_{\nu\rho}) .
$$
第一項は反対称性で零、後二項はともに $-\Omega_{\mu\nu}\Omega_{\mu\nu}$ を与えるので $\tr_{\rm spin}E^2=-2\Omega_{\mu\nu}\Omega_{\mu\nu}$。従ってスピンを取った $E^2/2$ は $-\Omega^2$、接続の項は $4\Omega^2/12$。ベクトルの場合は $\tr_{\rm vec}E^2=4\sum_{\mu,\nu}\Omega_{\mu\nu}\Omega_{\nu\mu}=-4\Omega^2$ であり、Grassmann のゴーストを含めた重み付き係数と比較できる。次元続行した内部成分数の $\mathcal O(\epsilon)$ の違いは、この一ループの極の係数を変えない。

:::example[作用から極の符号を照合する]
熱核の重み付き係数を一ループ作用へ入れると $\Tr\log$ の積分表示に負号が一つある。純ゲージなら極は $-(11C_A/12)(g^2/16\pi^2\bar\epsilon)\int\bar F^2$。これを消す対項には反対符号が必要なので、$S_E=\int\bar F^2/4$ の係数の極は正の $11C_Ag^2/(48\pi^2\bar\epsilon)$ となる。Dirac のみが荷電する Abel 理論では $a_4$ の符号が逆で、$\delta Z_{\bar A}<0$。この二つの例で、行列式の符号と対項の符号を混同すると $\beta$ が両方とも逆向きになることが分かる。
:::

### 背景熱核における縮退と磁気結合

定数な可換背景で、接続の色行列を対角化して一つの荷電成分 $q$ に限る。二次元の箱 $L_xL_y$ で $p_y=2\pi n/L_y$ の振動子の中心が箱内にある範囲は長さ $|qB|L_x$ だから、状態数は $|qB|L_xL_y/(2\pi)$。横の熱核の面積密度は
$$
\frac{|qB|}{2\pi}\sum_{n=0}^\infty e^{-(2n+1)t|qB|}
=\frac{|qB|}{4\pi\sinh(t|qB|)} .
$$
もう二つの方向の自由運動の熱核を掛けると $(4\pi t)^{-2}(t|qB|/\sinh(t|qB|))$ となる。この計算で $1/12$ の符号は、$[D_\mu,D_\nu]$ を反 Hermite とする規約で $\Omega_{\mu\nu}\Omega_{\mu\nu}$ が負であることに対応する。仮に $F_{\mu\nu}$ を直接 Hermite 行列で使うなら、$\Omega=-igF$ の $i^2$ を別に入れなければならない。

ベクトルでは $E_{\mu\nu}$ に背景の磁気結合が入る。上の二平面では二成分の固有値が $\pm 2qB$、残りの二成分が零で、$\tr_{\rm vec}e^{-t E}=2+2\cosh(2tqB)$。$z=tqB$ とおくと $4+4z^2+\cdots$。軌道の熱核を掛けると二次は $4-4/6=10/3$ である。ベクトルの $1/2$ とゴーストの $-1$ を加えると二次は $5/3+1/6=11/6$。$\Omega^2=-2q^2B^2$ と比較して $-11/12$ を得る。Dirac のスピン因子は $\pm qB$ が各二つなので $4\cosh z$。軌道と掛けた二次は $2-4/6=4/3$、$\Omega^2$ に対する係数は $-2/3$ である。磁気結合の違いと Gaussian の統計の違いの両方が、遮蔽と反遮蔽の符号を決める。

:::example[基本表現と随伴表現の規格化を区別する]
$SU(2)$ の可換背景を $T^3$ に選ぶ。基本の固有値は $+1/2,-1/2$ なので二乗和は $1/2$。随伴の固有値は $+1,0,-1$ で二乗和は $2$。従って同じ背景の $F^2$ から読んでも荷電成分を単に数えるのではなく、電荷の二乗を足す必要がある。例えば一つの基本 Dirac と一つの複素随伴スカラーを入れると $b_0=22/3-2/3-2/3=6$。随伴スカラーを実にすると最後の項は $1/3$、同じ模型の $b_0$ は $19/3$ になる。
:::

## 漸近的自由性

四次元で $b_0>0$ なら $\dd(1/g^2)/\dd\log\mu=b_0/(8\pi^2)$。従って高い尺度で
$$
g^2(\mu)\simeq\frac{8\pi^2}{b_0\log(\mu/\Lambda)}
$$
となる。結合が高エネルギーで零へ近づくことを :term[漸近的自由性]{en="asymptotic freedom" yomi="ぜんきんてきじゆうせい"} という。$SU(3)$ の基本 Dirac が $n_f$ 種なら $b_0=11-2n_f/3$。例えば $n_f=6$ なら $b_0=7$ である。$\Lambda$ は積分定数であり、先頭摂動式の分母が零になる尺度の近くにはこの式を適用できない。結合が赤外で増えることだけから、色を持つ漸近粒子が存在しないことまでは証明できない。

### 物質を加えたときの境界

$N$ 色の基本 Dirac 場を $n_f$ 種、質量を無視できる範囲で加えると $C_A=N$、$T(F)=1/2$ から $b_0=(11N-2n_f)/3$。従って先頭で紫外自由になる条件は $n_f<11N/2$ である。$N=3$ なら $n_f=16$ で $b_0=1/3$、$n_f=17$ で $b_0=-1/3$ になる。$b_0$ が小さい場合は低尺度で高次係数の相対寄与が大きくなり、一ループだけで赤外固定点の有無を判断できない。

四次元で $\alpha_g=g^2/(4\pi)$ を使うと、:eqref[eq-qft-ym-beta] の変数変換で $\mu\,\dd\alpha_g/\dd\mu=-b_0\alpha_g^2/(2\pi)+\cdots$ である。活性粒子を固定できる二つの高尺度の間では $\alpha_g(\mu_2)^{-1}-\alpha_g(\mu_1)^{-1}=b_0\log(\mu_2/\mu_1)/(2\pi)$。途中に重い質量がある場合は両側の $b_0$ を変え、閾値で同じ物理振幅へ照合する。

### 重い色荷が閾値をまたぐ例

$SU(3)$ で尺度 $\mu_0=10M$ において六つの基本 Dirac 場を活性とし、一つだけ質量 $M$ で残りを無質量とする。$\mu_1=M/10$ がなお摂動領域なら、閾値で先頭照合したとき
$$
\alpha_g(\mu_1)^{-1}-\alpha_g(\mu_0)^{-1}
=-\frac{7+23/3}{2\pi}\log 10 .
$$
上側の $b_0$ は $7$、下側は $23/3$ であり、二つの対数区間は同じ長さである。$M$ を切替点としたこの近似の補正には、有限な照合と高次の $\beta$ が入る。低い尺度でこの式の逆結合が零に近づけば、続けて同じ近似を使えるという仮定が崩れる。

### 閾値を動かして照合の対数を確かめる

先頭で上側と下側の結合を等しいと置くのは、照合する尺度を $\mu=M$ と選んだときである。別の尺度 $\mu_{\rm m}$ で逆結合の差 $\Delta(\mu_{\rm m})=\alpha_L^{-1}(\mu_{\rm m})-\alpha_H^{-1}(\mu_{\rm m})$ を取ると、一ループで $\dd\Delta/\dd\log\mu_{\rm m}=(b_L-b_H)/(2\pi)$。従って $\Delta(M)=0$ から
$$
\Delta(\mu_{\rm m})=\frac{b_L-b_H}{2\pi}\log\frac{\mu_{\rm m}}M .
$$
上側から $\mu_{\rm m}$ まで走り、この差を足し、そこから下側の観測尺度まで走れば、$\mu_{\rm m}$ の対数は両区間の端点と照合の間で相殺する。例えば基本 Dirac を一種除くと $b_L-b_H=2/3$。$\mu_{\rm m}=2M$ で照合の差を落とした計算には、逆結合で $\log 2/(3\pi)$ のずれが残る。ここでの $\mu_{\rm m}$ は有効理論を交換する計算上の尺度である。完全な摂動次数で照合した物理量は、この任意の尺度には依存しない。

### 高エネルギーの対数と次の次数

境界値 $\alpha_\star=\alpha_g(\mu_\star)>0$ を所与とすると、一ループの解は
$$
\alpha_g(\mu)=
\frac{\alpha_\star}{1+\dfrac{b_0\alpha_\star}{2\pi}\log(\mu/\mu_\star)}.
$$
$\Lambda=\mu_\star\exp[-2\pi/(b_0\alpha_\star)]$ と置けば $\alpha_g$ は無次元の境界値を使わずに $\mu/\Lambda$ で表せる。この置換を :term[次元転換]{en="dimensional transmutation" yomi="じげんてんかん"} と呼ぶ。$\Lambda$ は質量の次元を持ち、くりこみ処方によって値が変わる。例えば $b_0=7$、$\alpha_\star=0.1$ で $\log(\mu/\mu_\star)=2\pi/(7\alpha_\star)$ なら $\alpha_g(\mu)=0.05$。結合の半減には尺度の大きな比が必要で、これは $1/\log\mu$ で減ることの具体例である。

四次元で一つのゲージ結合だけの質量独立な処方を考え、$a=g^2/(16\pi^2)$ と置く。$d a/d\log\mu=-2b_0a^2-2b_1a^3+\cdots$ とする。$a'=a+c a^2+\cdots$ なら $d a'/d\log\mu=(1+2ca+\cdots)d a/d\log\mu$。$a=a'-c a'^2+\cdots$ を代入すると $b_0,b_1$ は同じであり、$a'^3$ の余分な二項が相殺する。三次以降の再定義の係数は一般により高い $\beta$ 係数を変える。二ループ係数の図式計算には重なった発散の差し引きが要るため :cite[weinberg-qft2] を参照する。基本 Dirac 場だけを持つ単純群の結果は
$$
b_1=\frac{34}{3}C_A^2-
\left(\frac{20}{3}C_A+4C_F\right)T(F)n_f .
$$
$C_F$ は基本表現で $\sum_aT^aT^a=C_F1$ を満たす二次 Casimir である。

:::example[小さい赤外固定点の条件]
二ループ式で $b_0>0,b_1<0$ なら $a_*=-b_0/b_1$ に正の零がある。$\beta_a'(a_*)=-2b_1a_*^2>0$ なので $\log\mu$ を減らす向きに近づく。例えば $SU(3)$、$n_f=16$ で $C_A=3,C_F=4/3,T(F)=1/2$ を入れると $b_0=1/3,b_1=-302/3$、従って $a_*=1/302$。クォークは質量ゼロと仮定した。小さい $a_*$ では高次の寄与を摂動的に見積もれるが、$b_0$ が大きくなって固定点が強結合へ動くと二ループまでの判断は使えない。
:::

### 結合の走りと物理的な対数

単一の大きい運動量 $Q$ で決まる無次元の物理量に、摂動展開のある寄与が $K a(\mu)+a(\mu)^2[c+d\log(Q^2/\mu^2)]$ を持つとする。これに $\mu$ 微分を作用させると、二次までに $-2K b_0 a^2-2da^2$ が残る。尺度独立性から $d=-K b_0$。$\mu=Q$ に取れば対数を零にでき、別の尺度で現れたこの対数列は走る結合の展開に移る。例えば $Q$ が閾値をまたぐ場合は、単一の $b_0$ でこの対数を延長せず、各有効理論で一定の活性場を使う。摂動的な漸近的自由性からはこの短距離での計算法が従う。長距離の面積則は [非摂動的手法](/physics/quantum-field-theory/nonperturbative/) で取り上げる。

### 二ループの積分定数と近似の範囲

$a=g^2/(16\pi^2)$、$t=\log\mu$ として、二ループの小さい $a$ では
$$
\frac{\dd t}{\dd a}
=-\frac1{2b_0a^2}+\frac{b_1}{2b_0^2a}+\mathcal O(1).
$$
積分により、適切な定数 $\Lambda$ を用いて
$$
\log\frac{\mu}{\Lambda}
=\frac1{2b_0a}+\frac{b_1}{2b_0^2}\log(b_0a)+\mathcal O(a).
$$
$L=\log(\mu^2/\Lambda^2)$ とし、$a=1/(b_0 L)+c\,\log L/L^2+\cdots$ を右辺に戻す。$a^{-1}=b_0 L-b_0^2 c\,\log L+\cdots$、$\log(b_0a)=-\log L+\cdots$ から、$c=-b_1/b_0^3$。例えば $b_1>0$ なら、同じ定義の $\Lambda$ に対する大きい $L$ での二ループの対数補正は負になる。この式は $L\gg1$ を要し、$\Lambda$ に近い領域の極の構造を与えるものではない。

結合を変える再定義 $a'=a+c a^2+\cdots$ では $(a')^{-1}=a^{-1}-c+\cdots$ だから、上の規約の $\Lambda'$ は $\Lambda e^{c/(2b_0)}$。二つの処方が同じ $b_0,b_1$ を持つことと、積分定数の数値が同じであることは別である。物理量を比較するには同じ励起や同じ短距離の振幅へ照合する。例えば一ループの閾値で $\alpha^{-1}$ が連続で、上側・下側の係数を $b_H,b_L>0$ とすれば $b_H\log(M/\Lambda_H)=b_L\log(M/\Lambda_L)$。従って $\Lambda_L=M^{1-b_H/b_L}\Lambda_H^{b_H/b_L}$。$M$ を使った式の右辺は全体として次元一であり、閾値より下の結合を同一の $\Lambda_H$ で走らせることはできない。

## 演習問題

::::exercise{#exr-qft-gauge-quantization-1 level=basic}
Abel 理論のゴーストが光子と結合しない理由を示せ。

:::solution
随伴共変微分の交換子が零なので $M=\Box$ は $A$ に依存しない。Grassmann 積分の行列式は規格化に吸収される。
:::
::::

::::exercise{#exr-qft-gauge-quantization-2 level=basic}
$s^2\bar c$ を求めよ。

:::solution
:eqref[eq-qft-brst]より $s\bar c=B$、$sB=0$ なので零である。
:::
::::

::::exercise{#exr-qft-gauge-quantization-3 level=basic}
$SU(2)$ の純ゲージ理論で $b_0$ を求めよ。

:::solution
$C_A=2$、Dirac 場なしより $b_0=22/3$ である。
:::
::::

::::exercise{#exr-qft-gauge-quantization-4 level=basic}
$SU(3)$、基本 Dirac $n_f=3$ の $b_0$ を求めよ。

:::solution
$C_A=3$、$T(F)=1/2$ を代入し $11-2n_f/3=9$ である。
:::
::::

::::exercise{#exr-qft-gauge-quantization-5 level=advanced}
$s^2\psi=0$ を奇微分の Leibniz 則で導け。

:::solution
$s^2\psi=ig[(igc^2)\psi-c(igc\psi)]$ となる。二項が逆符号なので零である。
:::
::::

::::exercise{#exr-qft-gauge-quantization-6 level=advanced}
補助場 $B$ を消してゲージ固定項を得よ。

:::solution
$B\,F+\xi B^2/2=\xi(B+F/\xi)^2/2-F^2/(2\xi)$。$\xi\ne0$ なら運動方程式で平方が零になる。
:::
::::

::::exercise{#exr-qft-gauge-quantization-7 level=advanced}
$g(\mu_0)=g_\star$ として一ループの $g(\mu)$ を求めよ。

:::solution
:eqref[eq-qft-ym-beta] の四次元式を積分して $1/g^2(\mu)=1/g_\star^2+(b_0/(8\pi^2))\log(\mu/\mu_0)$。正の枝の平方根の逆が求める結合である。
:::
::::

::::exercise{#exr-qft-gauge-quantization-8 level=advanced}
四次元熱核のベクトルとゴーストの群トレースに掛かる係数を再現せよ。

:::solution
ベクトルは $4/12-2=-5/3$、Gaussian の係数は $1/2$。ゴーストは $1/12$ に $-1$ を掛ける。和は $-5/6-1/12=-11/12$ である。
:::
::::

::::exercise{#exr-qft-gauge-quantization-9 level=advanced}
純ゲージのベクトルに加え、随伴の Weyl 場四個と実スカラー六個を含む理論を考える。これらが全て活性なときの一ループ $b_0$ を求めよ。一個の実スカラーだけが重く、その直下では除去されたとすると、閾値の下での係数はいくらか。

:::solution
随伴の $T(R)=C_A$。質量ゼロの Weyl は Dirac の半分のスピン成分を与えるので、一個あたり $-2C_A/3$。実スカラーは一個あたり $-C_A/6$ である。従って全て活性なら
$$
b_0=\left(\frac{11}{3}-4\frac23-6\frac16\right)C_A=0.
$$
実スカラーを一個取り除くと負の項が一つ減るので $b_0=C_A/6$。ここでは閾値の両側で活性な場だけを使う先頭の有効理論を考えている。上側で一ループ係数が零であるという計算だけでは、任意の相互作用を持つ同じ粒子内容の理論が全次数で尺度不変とは結論できない。
:::
::::

::::exercise{#exr-qft-gauge-quantization-10 level=advanced}
$SU(3)$ の Feynman ゲージで、純ゲージの一ループの $Z_A,Z_c,Z_g$ を用いる。三ゲージ、四ゲージ、ゴースト–ゲージ頂点に必要な極の係数を $g^2/(16\pi^2\bar\epsilon)$ の単位で求めよ。量子場の二点だけを背景の二点として使った場合、$Z_g$ の係数は何になるか。

:::solution
この単位で $\delta Z_A=(5/3)C_A$、$\delta Z_c=C_A/2$、$\delta Z_g=-(11/6)C_A$ である。従って三ゲージは $\delta Z_g+(3/2)\delta Z_A=(-11/6+5/2)C_A=(2/3)C_A=2$。四ゲージは $2\delta Z_g+2\delta Z_A=(-11/3+10/3)C_A=-C_A/3=-1$。ゴースト頂点は $\delta Z_g+\delta Z_c+\delta Z_A/2=(-11/6+1/2+5/6)C_A=-C_A/2=-3/2$。ここでは各式の $Z$ を一次まで展開している。誤って量子場の因子を $Z_g Z_A^{1/2}=1$ に使うと $\delta Z_g=-(5/6)C_A=-5/2$ になる。背景場の非斉次ゲージ変換に対して成り立つ関係を量子場の規格化へ使ったため、正しい値 $-11/2$ と一致しない。
:::
::::

::::exercise{#exr-qft-gauge-quantization-11 level=advanced}
純ゲージで、$SU(4)$ の可換背景に $H=\diag(1,-1,0,0)/2$ を選ぶ。基本と随伴における $H$ の二乗のトレースを求めよ。さらに四つの基本 Dirac、多重度一の随伴実スカラーが全て活性なら $b_0$ はいくらか。

:::solution
$H$ は $\tr_F H^2=1/2$ なので規格化済みである。随伴の作用は交換子 $[H,\cdot]$。行列単位 $E_{ij}$ への固有値は $h_i-h_j$、ただし $h_i$ は $H$ の対角成分。対角の三方向はすべて中性。非対角で $i,j=1,2$ の二順序は $+1,-1$、$i,j$ の一方だけが $3,4$ なら八順序で絶対値 $1/2$、双方が $3,4$ なら零。従って随伴トレースは $2+8/4=4$。この $H$ では $C_A=4$ と $T(F)=1/2$ を直接得た。実随伴スカラーの寄与を $-4/6$ として
$$
b_0=\frac{44}{3}-\frac43\cdot4\cdot\frac12-\frac46
=\frac{34}{3}.
$$
随伴を複素場と誤って数えるとスカラーの寄与を二倍にしてしまう。背景の荷電成分を対にまとめても、場の実性を決めるのは元の Gaussian の規格化である。
:::
::::

::::exercise{#exr-qft-gauge-quantization-12 level=advanced}
$\mu$ が閾値 $M$ を下向きに通る。上側・下側で一ループ係数をそれぞれ $7,23/3$ とし、$\alpha_g(M)=\pi/14$ で先頭照合する。$\Lambda_H,\Lambda_L$ を $M$ の倍数で表し、$\mu=M/e$ での下側の逆結合を求めよ。この尺度で、上側の結合をそのまま走らせた結果と比較せよ。

:::solution
一ループの $\alpha_g^{-1}=b_0\log(\mu/\Lambda)/(2\pi)$ に $\alpha_g(M)^{-1}=14/\pi$ を代入する。従って $\log(M/\Lambda_H)=28/7=4$、$\log(M/\Lambda_L)=28/(23/3)=84/23$ であり、$\Lambda_H=M e^{-4}$、$\Lambda_L=M e^{-84/23}$。$\mu=M/e$ で対数区間は $-1$ なので
$$
\alpha_g(M/e)^{-1}=\frac{14}{\pi}-\frac{23}{6\pi}=\frac{61}{6\pi}.
$$
上側を延長した結果は $14/\pi-7/(2\pi)=63/(6\pi)$。下側で重いフェルミオンの遮蔽がなくなると $b_0$ が増えるため、同じ下向きの区間で逆結合がより減る。この計算は閾値で結合の有限な高次補正を落としている。得た $\mu=M/e$ の比較は二つの走りの差を示すもので、重い粒子の運動学的な生成率を近似したものではない。
:::
::::

## 確認テスト

::::quiz
ゴースト頂点に現れる運動量はどの微分に由来するか。

- [x] $\partial_\mu\bar c$
- [ ] $\partial_\mu A^\mu$ だけ
- [ ] $B$ の時間微分

:::explanation
$\bar c\,\partial^\mu(D_\mu c)$ を部分積分すると微分は反ゴーストに作用する。$D$ のうち接続を含む部分は他の微分を含まない。
:::
::::

::::quiz{answer="7" tol="1e-6"}
$SU(3)$ で基本 Dirac が六種の場合の $b_0$ を答えよ。

:::explanation
$11C_A/3=11$、物質は $(4/3)(6/2)=4$ を引くので七になる。
:::
::::

::::quiz
BRST の物理状態空間を選べ。

- [x] $\ker Q_B/\operatorname{im}Q_B$ のゴースト数零の部分
- [ ] ゴーストを含む全空間
- [ ] $\operatorname{im}Q_B$ だけ

:::explanation
冪零性から像は核に含まれる。ゴースト数零で閉条件を満たし、BRST 完全な状態の差による同値関係で割った空間を用いる。
:::
::::

::::quiz{answer="11" tol="1e-6"}
$SU(3)$ の純ゲージ理論の $b_0$ を答えよ。

:::explanation
$C_A=3$ で物質がないため $11\cdot3/3=11$ である。
:::
::::

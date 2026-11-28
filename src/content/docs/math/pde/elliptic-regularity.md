---
title: 楕円型方程式の正則性と最大値原理
description: 差分商による二階弱微分の評価を証明し、最大値原理、Harnack 不等式、Hölder 評価の関係を調べる。
notation:
  - symbol: 'D_h^j u'
    meaning: 第 j 方向の差分商
  - symbol: '\nu'
    meaning: 領域の外向き単位法線
  - symbol: 'u_+'
    meaning: u の正部分
---

## この章の目標

- 差分商で弱解の二階弱微分を評価する。
- 境界を平坦化し、接方向と法線方向の微分を回復する。
- 弱最大値原理、Hopf の補題、強最大値原理を使い分ける。
- Harnack、Schauder、De Giorgi–Nash–Moser の各評価に必要な仮定を比較する。

## 差分商と内部正則性

正則性を示すため方程式を微分したくても、弱解ではその微分が存在するとは限らない。そこで微分の代わりに平行移動の有限差を使い、刻み幅を零にするときの一様評価を得る。

:::definition[差分商]{#def-pde-difference-quotient}
$D_h^jv(x)=(v(x+he_j)-v(x))/h$ を :term[差分商]{en="difference quotient" yomi="さぶんしょう"} という。積分領域が平行移動を許すか、試験関数が内部に台を持てば $\int D_h^jv\,\phi=-\int v\,D_{-h}^j\phi$ である。
:::

最後の式は変数変換で従う。滑らかな $v$ では $D_h^j v(x)=\int_0^1\partial_jv(x+\theta he_j)\,\dd\theta$ なので、Jensen の不等式から差分商の $L^2$ ノルムは少し広い領域での $\partial_jv$ の $L^2$ ノルム以下である。平滑化で同じ評価を $H^1$ に拡張できる。逆に差分商が $L^2$ で有界なら、$h\to0$ の部分列で弱収束する。上の積分関係から極限は $\partial_jv$ の分布微分であり、そのノルムは一様上界以下である。

:::theorem[内部 $H^2$ 評価]{#thm-pde-elliptic-interior-h2}
$A$ は一様楕円で局所 Lipschitz、$f\in L^2_{\rm loc}(\Omega)$ とする。$u\in H^1_{\rm loc}$ が $-\divg(A Du)=f$ の弱解であれば $u\in H^2_{\rm loc}(\Omega)$ であり、$U\Subset V\Subset\Omega$ なら

$$
\norm{u}_{H^2(U)}\le C\big(\norm{f}_{L^2(V)}+\norm u_{H^1(V)}\big)
\label{eq-pde-interior-h2}
$$

となる。定数は楕円性、$A$ の Lipschitz ノルム、切断の幅に依存する。
:::

:::proof
$\zeta=1$ を $U$ 上で満たす $C_c^\infty(V)$ を選ぶ。方向を固定し $w=D_hu$、試験関数を $-D_{-h}(\zeta^2w)$ とする。離散部分積分と積の法則により

$$
\int A(x+he_j)Dw\cdot D(\zeta^2w)
=-\int (D_hA)Du\cdot D(\zeta^2w)-\int fD_{-h}(\zeta^2w).
$$

左辺は主項と交差項の和である。主項は $\lambda\norm{\zeta Dw}_2^2$ 以上であり、交差項の絶対値は $2\norm A_\infty\norm{\zeta Dw}_2\norm{wD\zeta}_2$ 以下である。右辺では $\norm{D_hA}_\infty$ が有界であり、$\norm{D_{-h}(\zeta^2w)}_2\le\norm{D(\zeta^2w)}_2$ を用いられる。$D(\zeta^2w)=\zeta^2Dw+2\zeta wD\zeta$ なので、左辺の交差項を右辺へ移すと、右側にある $\norm{\zeta Dw}_2$ を含む積の総和は Young の不等式で $\lambda\norm{\zeta Dw}_2^2/2$ と $C(\norm f_2^2+\norm{Du}_2^2)$ に分けられる。ここで $\norm w_2\le\norm{Du}_2$ を使った。従って $D_h Du$ は $U$ 上で一様に有界である。差分商の逆の評価を各方向に適用すると結論を得る。
:::

有界な下位項は $f$ から $b\cdot Du+cu$ を引けば同じ評価に含められる。高階正則性では、例えば $A$ が滑らかで $f\in H^1_{\rm loc}$ なら $-\divg(A D\partial_j u)=\partial_jf+\divg((\partial_jA)Du)$ を得る。$u\in H^2$ のため右辺は $L^2_{\rm loc}$ であり、再度定理を使い $u\in H^3_{\rm loc}$ を得る。同じ微分の手順を積の法則で帰納すれば、滑らかな係数と右辺には内部の任意の次数の Sobolev 正則性が従う。Sobolev 埋め込みで古典的にも滑らかである。

:::example[可測係数で二階正則性が失われる]
$(-1,1)$ で $a=1$（左）、$a=2$（右）とし、$u=x$（左）、$u=x/2$（右）とする。$a u'=1$ なので弱い斉次方程式を満たすが、$u''=-\delta_0/2$ は $L^2$ に属さない。Lipschitz 係数という仮定を有界可測性に落とすことはできない。
:::

:::example[平行移動が許す方向]
平面境界 $x_n=0$ でトレースが $0$ の関数 $u(x)=x_n$ を考える。接方向の差分商は $0$ でトレースも $0$ だが、法線方向では $D_h^n u=1$ で境界でも $1$ である。境界付近で内部証明の試験関数を法線方向にずらすと、零境界値を失い得る。
:::

### Caccioppoli 評価から内部ノルムを減らす

$u$ が $-\divg(A Du)=f$ の弱解とし、試験関数を $\eta^2 u$ とする。左辺は $\int \eta^2 A Du\cdot Du+2\int\eta u A Du\cdot D\eta$ である。第二項を Young の不等式で $\lambda\norm{\eta Du}_2^2/2+C\norm{uD\eta}_2^2$ に分ける。半径 $r<R$ の同心球で $\eta=1$ を $B_r$ に取り、$|D\eta|\le C/(R-r)$ とすると
$$
\norm{Du}_{L^2(B_r)}^2\le C\left((R-r)^{-2}\norm u_{L^2(B_R)}^2
+\norm f_{L^2(B_R)}\norm u_{L^2(B_R)}\right).
$$
これを内部 $H^2$ 評価で用いれば、少し広い球での $H^1$ ノルムの代わりに $L^2$ ノルムと源項を使える。定数解の勾配に関心がある場合、$u$ を $u-m$ に置き換えても発散形の式は変わらない。そのため $m$ を球の平均と取って、測定する量を振動に絞れる。

:::example[係数の微分と曲率の寄与]
$(0,1)$ で $A=1+2x$、零境界を考える。$u=x(1-x)$ を候補にして、必要な外力を決める。

まず流束 $q=Au_x$ を未知量にする。$q_x=-f$ を積分すると、任意の定数 $q(0)$ について $q(x)=q(0)-\int_0^x f$ である。さらに $q/A$ を積分し、指定された両端の値または平均に合わせて積分定数を決める。以下では計算後の $u$ について関係を逆に確かめる。

$$
\begin{aligned}
u(x)&=x \left(1 - x\right)\\
u_x&=\frac{d}{d x} x \left(1 - x\right)=1 - 2 x\\
q&=\left(1 - 2 x\right) \left(2 x + 1\right)=- \left(2 x - 1\right) \left(2 x + 1\right)\\
-f&=\frac{d}{d x} \left(- \left(2 x - 1\right) \left(2 x + 1\right)\right)=- 8 x
\end{aligned}
$$

次に値の条件と流出量を両端で区別して評価する。区間の外法線は左端で $-1$、右端で $1$ である。このため源項の積分は左側の流束から右側の流束を引いたものになる。流束を係数で割った積分は両端の値の差に等しい。

$$
\begin{aligned}
u(0)&=0,\qquad x=0\\
u(L)&=0,\qquad x=1\\
q(0)&=1,\qquad x=0\\
q(L)&=-3,\qquad x=1\\
\int_0^L f\,\mathrm{d}x&=\int\limits_{0}^{1} 8 x\, dx=4\\
\int_0^L (q/A)\,\mathrm{d}x&=\int\limits_{0}^{1} \left(1 - 2 x\right)\, dx=0
\end{aligned}
$$

弱形式も具体的に比べる。試験関数 $\phi=x(L-x)$ は両端で零であり、積分の間に境界項は出ない。以下の値が一致するのは、任意の零トレース試験関数に対しても部分積分できることの一例である。

$$
\begin{aligned}
\phi&=x \left(1 - x\right)\\
\phi^{\prime}&=\frac{d}{d x} x \left(1 - x\right)=1 - 2 x\\
\int_0^L q\phi^{\prime}\,\mathrm{d}x&=\int\limits_{0}^{1} \left(- \left(1 - 2 x\right) \left(2 x - 1\right) \left(2 x + 1\right)\right)\, dx=\frac{2}{3}\\
\int_0^L f\phi\,\mathrm{d}x&=\int\limits_{0}^{1} 8 x^{2} \left(1 - x\right)\, dx=\frac{2}{3}
\end{aligned}
$$

最後に内部での値、勾配、流束を同じ点で測る。区間全体での勾配の寄与には $A$ の重みが付くので、点での値だけでは全体のエネルギーは分からない。

$$
\begin{aligned}
u&=\frac{3}{16},\qquad x=\frac{1}{4}\\
u_x&=\frac{1}{2},\qquad x=\frac{1}{4}\\
q&=\frac{3}{4},\qquad x=\frac{1}{4}\\
\frac12\int_0^L A(u_x)^2\,\mathrm{d}x&=\int\limits_{0}^{1} \frac{\left(1 - 2 x\right)^{2} \left(2 x + 1\right)}{2}\, dx=\frac{1}{3}
\end{aligned}
$$

流束の微分では $A u_{xx}=-2-4x$ と $A_xu_x=2-4x$ が加わり、定数項は相殺する。$u_{xx}=-2$ なので $x=0$ では外力が $0$ でも曲率は零ではない。差分商評価で係数の Lipschitz 定数が現れる理由はこの積の法則にある。
:::

## 境界までの正則性

:::theorem[Dirichlet 境界での $H^2$ 評価]{#thm-pde-elliptic-boundary-h2}
有界 $C^{1,1}$ 領域で $A$ は一様楕円かつ Lipschitz、$f\in L^2$ とする。$u\in H^1_0$ が $-\divg(A Du)=f$ の弱解なら $u\in H^2$ であり $\norm u_{H^2}\le C(\norm f_2+\norm u_{H^1})$ である。
:::

:::proof
境界を $x_n=\psi(x')$ と表し、$y_n=x_n-\psi(x')$ の座標で半領域に変える。この写像の行列式は $1$ で、変分形式の係数は $B=D\Phi A(D\Phi)^\top$ になる。$D\Phi$ とその逆は有界で Lipschitz なので、一様楕円性と Lipschitz 性を保つ。$u$ の境界値は平面 $y_n=0$ 上で $0$ になる。接方向の差分商は零トレースを保つ。従って半球に支えられた切断を使い :ref[thm-pde-elliptic-interior-h2] と同じ評価を行うと、接方向に少なくとも一回微分した全ての二階弱微分が $L^2$ に入る。

分布の式を積の法則で展開すると

$$
-B_{nn}\partial_{nn}u
=f+\sum_{(i,j)\ne(n,n)} B_{ij}\partial_{ij}u+\sum_{i,j}(\partial_iB_{ij})\partial_j u
$$

である。右辺の各項は $L^2$ に入る。より厳密には、まず $\partial_n(B_{nn}\partial_n u)$ を左に置き、$\partial_n$ を含む他の項では接方向の二階微分が既知であることを用いる。次に $1/B_{nn}$ は有界 Lipschitz なので積の法則で $\partial_{nn}u$ を回復できる。有限個の境界パッチと内部の評価を合わせると全領域での結論が従う。
:::

:::example[角のある領域]
平面の角 $0<\theta<3\pi/2$ で $u=r^{2/3}\sin(2\theta/3)$ は、極座標の Laplace 作用素に代入すると調和で、両側で $0$ になる。勾配の大きさは $Cr^{-1/3}$ で局所二乗可積分である。一方、Hessian の大きさは一定係数を除き $r^{-4/3}$ であり、二乗の半径積分は $\int_0^1r^{-5/3}\,\dd r=\infty$ となる。外側で滑らかに切断すれば、右辺 $L^2$ で $H_0^1$ の弱解だが境界近くで $H^2$ でない例が得られる。
:::

:::example[接方向の微分から法線微分を求める]
上半平面で $A=\begin{pmatrix}2&1\\1&3\end{pmatrix}$、$u=x_1x_2+x_2^2$ とする。境界で $u=0$、$u_{11}=0$、$u_{12}=1$、$f=-\divg(A Du)=-8$ である。接微分から分かる項を移すと $-3u_{22}=-8+2=-6$ になり、法線二階微分 $2$ を回復できる。
:::

### 高階の Dirichlet 正則性

$m$ を非負整数とする。有界 $C^{m+1,1}$ 領域で $A\in W^{m+1,\infty}$ が一様楕円、$f\in H^m$ なら、零トレースの弱解は $H^{m+2}$ に属する。評価は $\norm u_{H^{m+2}}\le C(\norm f_{H^m}+\norm u_{H^1})$ である。$m=0$ は上の定理である。

帰納段階では $u\in H^{m+1}$ を既知とする。平坦化した境界パッチで $m$ 回の純接微分 $D_{y'}^\gamma u$ を取ると、平面上の跡は零である。これはトレースが接方向の微分と可換であり、元の跡が零であるためである。係数を微分した項を右辺へ移すと、例えば一つの項は $\divg((D_{y'}^\beta B)D D_{y'}^{\gamma-\beta}u)$ となる。$|\beta|\ge1$ なので、発散を展開しても解に掛かる微分は高々 $m+1$ 階である。従って右辺は $L^2$ であり、上の接方向差分商と法線回復の議論を局所パッチで再び使える。必要な切断を掛けると生じる項にも解の $m+1$ 階までの微分しか現れない。これで全階数 $m+2$、法線微分の回数が高々二回の項が求まる。

三回以上の法線微分は、$u_{nn}$ を他の二階微分と $f$ で表した式を使う。欲しい法線微分が $r$ 回なら、その式を法線方向に $r-2$ 回、接方向に $m+2-r$ 回微分する。右辺の同じ全階数の微分には、もとの二階項が $nn$ を含まないため、法線微分は高々 $r-1$ 回しか現れない。係数に微分が当たった項は $u\in H^{m+1}$ で評価できる。従って $r$ の小さい順に全て回復できる。これらの微分は弱微分として求める。平坦化した係数と $1/B_{nn}$ の $W^{m+1,\infty}$ 有界性は座標写像の $C^{m+1,1}$ 正則性と楕円性から従う。内点の評価と有限個のパッチの評価を合わせ、帰納仮定で中間の $H^{m+1}$ ノルムを抑えると上の評価になる。

### 非零の Dirichlet データを引く

境界データが $G\in H^2(\Omega)$ のトレースであるとする。$v=u-G\in H_0^1$ は $-\divg(A Dv)=f+\divg(A DG)$ を満たす。$A\in W^{1,\infty}$ の場合、積の微分は $A:D^2 G+(\divg A)\cdot DG$ であり、絶対値を $L^2$ で評価できる。従って :ref[thm-pde-elliptic-boundary-h2] により $\norm u_{H^2}\le C(\norm f_2+\norm u_{H^1}+\norm G_{H^2})$ となる。延長 $G$ の取り方はこの上界の大きさに影響するが、解自身は変わらない。

:::example[角に生じる指数を分離で求める]
平面の角領域 $0<\theta<\omega$ で $u=r^\beta\sin(\beta\theta)$ は原点の外で調和であり、$\beta=\pi/\omega$ なら二辺で零である。$|Du|$ は $r^{\beta-1}$ の定数倍なので $\int_0^1r^{2\beta-1}\dd r$ が有限である。二階微分は $\beta\ne1$ で $r^{\beta-2}$ の成分を持つため $\int_0^1 r^{2\beta-3}\dd r$ が有限となるには $\beta>1$ が要る。例えば $\omega=3\pi/2$ なら $\beta=2/3$ で $H^1$ だが $H^2$ ではない。境界の曲率を有界にする仮定はこの角を排除する。
:::

:::example[アフィン延長が消さない外力]
$(0,2)$ で $A=2+x$、$f=-2$、$u(0)=1,u(2)=5$ とする。アフィン延長の流束が一定になるか調べる。

まず流束 $q=Au_x$ を未知量にする。$q_x=-f$ を積分すると、任意の定数 $q(0)$ について $q(x)=q(0)-\int_0^x f$ である。さらに $q/A$ を積分し、指定された両端の値または平均に合わせて積分定数を決める。以下では計算後の $u$ について関係を逆に確かめる。

$$
\begin{aligned}
u(x)&=2 x + 1\\
u_x&=\frac{d}{d x} \left(2 x + 1\right)=2\\
q&=2 x + 4=2 \left(x + 2\right)\\
-f&=\frac{d}{d x} 2 \left(x + 2\right)=2
\end{aligned}
$$

次に値の条件と流出量を両端で区別して評価する。区間の外法線は左端で $-1$、右端で $1$ である。このため源項の積分は左側の流束から右側の流束を引いたものになる。流束を係数で割った積分は両端の値の差に等しい。

$$
\begin{aligned}
u(0)&=1,\qquad x=0\\
u(L)&=5,\qquad x=2\\
q(0)&=4,\qquad x=0\\
q(L)&=8,\qquad x=2\\
\int_0^L f\,\mathrm{d}x&=\int\limits_{0}^{2} \left(-2\right)\, dx=-4\\
\int_0^L (q/A)\,\mathrm{d}x&=\int\limits_{0}^{2} 2\, dx=4
\end{aligned}
$$

弱形式も具体的に比べる。試験関数 $\phi=x(L-x)$ は両端で零であり、積分の間に境界項は出ない。以下の値が一致するのは、任意の零トレース試験関数に対しても部分積分できることの一例である。

$$
\begin{aligned}
\phi&=x \left(2 - x\right)\\
\phi^{\prime}&=\frac{d}{d x} x \left(2 - x\right)=- 2 \left(x - 1\right)\\
\int_0^L q\phi^{\prime}\,\mathrm{d}x&=\int\limits_{0}^{2} 2 \left(2 - 2 x\right) \left(x + 2\right)\, dx=- \frac{8}{3}\\
\int_0^L f\phi\,\mathrm{d}x&=\int\limits_{0}^{2} \left(- 2 x \left(2 - x\right)\right)\, dx=- \frac{8}{3}
\end{aligned}
$$

最後に内部での値、勾配、流束を同じ点で測る。区間全体での勾配の寄与には $A$ の重みが付くので、点での値だけでは全体のエネルギーは分からない。

$$
\begin{aligned}
u&=3,\qquad x=1\\
u_x&=2,\qquad x=1\\
q&=6,\qquad x=1\\
\frac12\int_0^L A(u_x)^2\,\mathrm{d}x&=\int\limits_{0}^{2} \left(2 x + 4\right)\, dx=12
\end{aligned}
$$

アフィン延長 $G=1+2x$ の曲率は零だが $-(AG_x)_x=-2$ である。従って $u-G$ の源項は零であり、両端でも零なので $u=G$ である。変数係数の発散形作用素はアフィン関数を必ずしも零にしない。
:::

## 弱最大値原理

:::proposition[正部分を用いる比較]{#prop-pde-weak-elliptic-maximum}
$\Omega$ は有界開集合とする。$L=-\divg(A D\cdot)+c$ で、$A$ は有界可測かつ一様楕円、$c$ は有界可測で非負とする。$u\in H^1$、$u_+\in H^1_0$ であり、全ての非負 $v\in H^1_0$ について $a(u,v)\le0$ なら $u\le0$ a.e. である。
:::

:::proof
$u_+=\max(u,0)$ を試験関数に選ぶ。Sobolev の連鎖律により $Du_+=\mathbf1_{\{u>0\}}Du$ a.e. である。これは正部分の Lipschitz 近似を用いて得られる。従って $a(u,u_+)\ge\lambda\norm{Du_+}_2^2+\int c u_+^2\ge0$ である。不等式と合わせると $Du_+=0$ であり、Poincaré の不等式から $u_+=0$ となる。
:::

古典的な非発散形の符号規約として $P=a_{ij}\partial_{ij}+b_i\partial_i+c$、$c\le0$ を用いる。$\Omega$ 有界、係数有界、一様楕円とし、$u\in C^2(\Omega)\cap C(\overline\Omega)$ が内部で $Pu\ge0$、境界で $u\le0$ を満たすとする。$\beta$ を大きく取ると $P(e^{\beta x_1})=(\beta^2 a_{11}+\beta b_1+c)e^{\beta x_1}>0$ である。$u+\eps e^{\beta x_1}$ は正の内部最大点を持てない。そこで最大点では勾配が零、Hessian は非正で、$c$ の項も非正だからである。従って境界値の最大と $0$ で上から抑え、$\eps\downarrow0$ とすれば $u\le0$ である。

:::example[符号が逆の場合]
$(0,\pi)$ で $P=\partial_{xx}+1$、$u=\sin x$ とすると $Pu=0$ かつ境界で $0$ だが、内部では正である。零階項の条件がなければ、内部最大点で $Pu\le0$ の符号は保証されない。
:::

:::example[弱解の順序を源項から求める]
単位球で $-\Delta u=f$、$u$ のトレースが $0$、$0\le f\le 2n$ とする。$-u$ と $u-(1-|x|^2)$ に正部分の原理を適用すると $0\le u\le1-|x|^2$ である。上界は $f=2n$ で達成される。境界上の点ごとの微分がなくても、この比較は弱形式で成立する。
:::

### 下位項を含む比較の条件

実係数の発散形 $Lu=-\divg(A Du)+b\cdot Du+cu$ を考える。$b$ が Lipschitz、$c-\divg b/2\ge0$ であり、$Lu\le0$、境界で $u\le0$ とする。$u_+\in H_0^1$ を試験すると $Du=Du_+$ が $u>0$ 上で成り立ち、移流項は $\int b\cdot D(u_+^2)/2=-\int (\divg b)u_+^2/2$ となる。従って $\lambda\norm{Du_+}_2^2+\int(c-\divg b/2)u_+^2\le0$ となり、Poincaré から $u_+=0$。係数に関するこの条件は、正部分を使ったエネルギー証明の十分条件である。

:::example[負の零次項が比較を破る]
$(0,\pi)$ で $L=-\partial_{xx}-2$、$u=\sin x$ とすると、境界値は零で $Lu=-\sin x\le0$ なのに内部は正である。$-2$ の項が正の勾配エネルギーを上回るためである。一方 $L=-\partial_{xx}-1/2$ では $\int|v'|^2\ge\int|v|^2$ なので $\int(|v'|^2-|v|^2/2)\ge\int|v'|^2/2$。後者では零次項が負でも同じ正部分の証明が使える。
:::

有界な $b$ に対しても、$c\ge0$ なら $Lu\le0$、$u_+\in H_0^1$ から $u\le0$ が従う。ここでは $A$ は有界かつ一様楕円、$c$ は有界とする。$k\ge0$ として $w=(u-k)_+$ を入れる。零次項は $\int cuw\ge0$ であり、ドリフト項は $E_k=\{u>k,\ Du\ne0\}$ にだけ寄与する。そこで $H_0^1(\Omega)\hookrightarrow L^q(\Omega)$ となる $q>2$ を固定し、Hölder と Sobolev の不等式を用いると

$$
\lambda\norm{Dw}_2^2
\le \norm b_\infty\norm{Dw}_2\norm{w}_{L^2(E_k)}
\le C\norm b_\infty |E_k|^{1/2-1/q}\norm{Dw}_2^2 .
$$

$M=\operatorname*{ess\,sup}u$ が正と仮定する。$M=+\infty$ なら $k\to+\infty$ で $|E_k|\to0$ である。$M<+\infty$ でも $k\uparrow M$ で $|E_k|\to0$ となる。実際、Sobolev 関数の勾配はレベル集合 $\{u=M\}$ 上で a.e. に零である。従って $k<M$ のまま右辺の係数を $\lambda$ 未満にできる。すると $Dw=0$ で $w=0$ となり、$k<M$ に矛盾する。

### 上界を作る障壁の使い方

$-\divg(A Du)\le f$ で、境界で $u\le0$ とする。非負な関数 $W$ が $-\divg(A DW)\ge f$ を満たせば、差の正部分を試験して $u\le W$ である。例えば区間 $(0,1)$ で $A(x)=1+x$、$f=1$ とする。$W(x)=x(1-x)$ なら $-((1+x)W')'=1+4x\ge1$ であるから $0\le u$ が分かっている場合には $0\le u\le1/4$ の上界を得る。比較定理の適用では $W$ が厳密な方程式を解く必要はなく、境界条件と外力の不等式がそろえばよい。逆に $f=-6$ の下界を作るなら、単純に $W$ の符号を変えるだけでなく、比較する源の大きさに合わせて係数を選ぶ。

## Hopf の補題と強最大値原理

:::proposition[Hopf の境界点補題]{#prop-pde-hopf-boundary}
$P$ は上の非発散形の仮定を満たし、$u\in C^2$ は $Pu\ge0$ とする。内部球 $B_R(y)$ が $\xi$ で境界に接し、$u$ は閉球で連続であり、$u<M=u(\xi)$ が球内で成り立つとする。$M\ge0$ とし、球の外向き半径方向 $\nu=(\xi-y)/R$ の内側からの微分が存在すれば $\partial_\nu u(\xi)>0$ である。
:::

:::proof
環 $R/2<r=|x-y|<R$ で $w=e^{-\beta r^2}-e^{-\beta R^2}$ を取る。$w$ の二階項の $\beta^2$ の部分は $4\beta^2 e^{-\beta r^2}(x-y)^\top A(x-y)$ である。二階項の残りと一階項は $C\beta e^{-\beta r^2}$ 以下、零階項は $C e^{-\beta r^2}$ 以下なので、大きい $\beta$ で $Pw>0$ になる。$P(u-M)=Pu-cM\ge0$ である。内周で $M-u$ の正の下限があるので、小さい $\delta>0$ に対して $u-M+\delta w\le0$ が内外周で成り立つ。最大値原理により環全体でも同じである。接点での半径差商を取ると $\partial_\nu u(\xi)\ge-\delta\partial_r w(R)=2\delta\beta R e^{-\beta R^2}>0$ である。
:::

:::theorem[強最大値原理]{#thm-pde-strong-elliptic-maximum}
連結領域で上の $P$ と $Pu\ge0$ を満たす $u$ が非負の最大値 $M$ を内点で取れば定数である。
:::

:::proof
$u<M$ の集合が空でなければ、連結性により領域内部に $u<M$ と最大点集合の境界点がある。その境界点に十分近い $u<M$ の点 $y$ を選び、最大点集合までの距離を $R$ とすれば、閉球は領域内にあり、球は最大点集合にある $\xi$ に接する。:ref[prop-pde-hopf-boundary] により $\xi$ で外向き微分は正である。一方 $\xi$ は領域の内部最大点なので $Du(\xi)=0$ であり矛盾する。
:::

:::example[球内の符号]
単位球の $u=|x|^2-1$ は $\Delta u=2n>0$ で、境界最大値は $0$ である。法線微分は $2$ なので Hopf の結論と符号が一致する。$-\Delta$ を正の作用素として書くときは劣解の不等号も反転する。
:::

:::example[内接球の条件と角]
第一象限の円盤部分で $u=-x_1x_2$ は調和で、内部では負、原点では $0$ であるが $Du(0)=0$ である。原点に接して象限内に含まれる正半径の球は存在しない。実際、中心 $(a,b)$ の距離が半径 $R=\sqrt{a^2+b^2}$ であるなら、象限に入る条件 $a,b\ge R$ と矛盾する。従って Hopf の結論をこの角に適用できない。
:::

### 境界点の非零な勾配を定量化するモデル

定数係数の環 $r_0<|x|<R$（$n\ge3$）で $u=0$ を内周、$u=1$ を外周に課す。放射状の調和解は
$$
u(r)=\frac{r_0^{2-n}-r^{2-n}}{r_0^{2-n}-R^{2-n}},\qquad
\partial_r u(R)=\frac{(n-2)R^{1-n}}{r_0^{2-n}-R^{2-n}}>0 .
$$
実際 $r^{n-1}u'$ は定数であり、二つの境界値で決まる。三次元で $r_0=1,R=2$ なら $u=2(1-1/r)$、外向き微分は $1/2$ である。内接球の環で作る一般の指数障壁は、このような放射状の補助関数で係数の誤差を上回るように調整したものである。

:::example[内外の境界で Hopf の符号を照合する]
四次元の環 $1<r<2$ で内側を $0$、外側を $1$ とする。外法線で測ると同じ $u_r>0$ でも端ごとに法線微分の符号が変わる。

放射方向だけに依存する解を考える。半径 $r$ の球面積が $\sigma_{n-1}r^{n-1}$ であることから、発散を計算すると $-(r^{n-1}u_r)^{\prime}=r^{n-1}f$ となる。これを一度積分し、中心を含む場合には中心に余分な点源がないよう流束の極限を零にする。環状領域では内側の境界条件で積分定数を決める。

$$
\begin{aligned}
u&=\frac{4}{3} - \frac{4}{3 r^{2}}\\
u_r&=\frac{d}{d r} \left(\frac{4}{3} - \frac{4}{3 r^{2}}\right)=\frac{8}{3 r^{3}}\\
r^{n-1}u_r&=\frac{8}{3}\\
-r^{n-1}f&=\frac{d}{d r} \frac{8}{3}=0\\
-\Delta u&=0
\end{aligned}
$$

次に値と流束を境界で読む。外側では法線は半径が増える向きである。内側にも境界がある場合は法線が反対向きになる。源項の体積積分を角度方向に積分すると面積係数が現れ、流束との差を確かめられる。

$$
\begin{aligned}
u&=0,\qquad r=1\\
u&=1,\qquad r=2\\
u_r&=\frac{1}{3},\qquad r=2\\
\sigma_{n-1}&=2 \pi^{2}\\
\text{全源量}&=\int\limits_{1}^{2} 0\, dr=0\\
\text{負の全流束}&=0
\end{aligned}
$$

観測半径で値と半径微分を計算する。勾配の二乗を用いてエネルギーも計算できる。観測点で勾配が小さい場合でも、球面積の因子があるため、そこを通る全流束が小さいとは限らない。

$$
\begin{aligned}
u&=\frac{20}{27},\qquad r=\frac{3}{2}\\
u_r&=\frac{64}{81},\qquad r=\frac{3}{2}\\
\sigma_{n-1}r^{n-1}u_r&=\frac{16 \pi^{2}}{3},\qquad r=\frac{3}{2}\\
\frac12\int |Du|^2&=\int\limits_{1}^{2} \frac{64 \pi^{2}}{9 r^{3}}\, dr=\frac{8 \pi^{2}}{3}
\end{aligned}
$$

外側では最大値を取るため外法線微分は正である。内側では最小値を取るため、同じ向きの半径微分を外法線に変えると負である。$-u$ に Hopf の境界点補題の向きを合わせれば、内側でも同じ定理でこの符号が得られる。
:::

## Harnack の不等式

正の解では最大値原理が等号の可否を決めるだけでなく、内球における二つの値の比を評価できる。Laplace 作用素では :ref[thm-pde-harmonic-harnack] により、非負性と球体平均だけで証明できる。定数対称行列 $A>0$ では $x=A^{1/2}y$ により $a_{ij}\partial_{ij}u=\Delta_yv$ となる。従って同じ評価が対応する楕円体上で成り立つ。

:::example[異方性の向き]
$A=\diag(4,1)$ では $y$ の単位円は $x$ の半軸 $2,1$ の楕円に写る。$v(y)=3+y_1$ を引き戻すと $u(x)=3+x_1/2$ である。この楕円内では $2<u<4$ で、値の比は $2$ 以下である。
:::

### 有界可測な係数での Moser の方法

:::theorem[Moser の Harnack 不等式]{#thm-pde-moser-harnack}
有界可測な一様楕円行列 $A$ に対し、非負の弱解 $u\in H^1(B_{2r})$ が $-\divg(A Du)=0$ を満たすとする。このとき

$$
\operatorname*{ess\,sup}_{B_r}u\le C\operatorname*{ess\,inf}_{B_r}u
\label{eq-pde-moser-harnack}
$$

であり、$C$ は次元と係数の上界と楕円定数の比だけに依存する。
:::

:::proof[正のべきと局所有界性]
平行移動と拡大で固定した大きさの球を考えればよい。$v=u+\delta$、$\delta>0$ とおく。$v$ も解である。立方体 $Q_R$ と同じ中心の小立方体 $Q_s$ を考え、$\zeta=1$ を $Q_s$ で満たし $Q_R$ に台を持つ切断で $|D\zeta|\le C/(R-s)$ とする。$p\ge2$ について $\psi(v)=v\min(v,N)^{p-2}$ を用いて $\zeta^2\psi(v)$ で試験する。$\psi$ は固定した $N$ で Lipschitz なので試験が許される。$\Phi(v)=v\min(v,N)^{(p-2)/2}$ と置くと各側で $(\Phi')^2\le p^2\psi'$、$\psi^2/\psi'\le\Phi^2$ が成り立つ。例えば $v\le N$ では $\psi'=(p-1)v^{p-2}$、$\Phi'=p v^{p/2-1}/2$ であり、$v>N$ ではどちらもべきを固定した線形式になる。楕円性と Young の不等式で $\int\zeta^2\psi'|Dv|^2\le C\int\Phi^2|D\zeta|^2$ である。よって $\norm{D(\zeta\Phi)}_2\le Cp\norm{\Phi D\zeta}_2$ となる。

$n>2$ なら $\chi=n/(n-2)$、$n\le2$ なら例えば $\chi=2$ とし、Sobolev 埋め込みを適用する。領域を固定した単位大の立方体に制限すれば $\norm{\zeta\Phi}_{2\chi}\le C\norm{D(\zeta\Phi)}_2$ である。$v\in L^p(Q_R)$ なら $N\uparrow\infty$ として

$$
\norm{v}_{L^{p\chi}(Q_s)}
\le\left(\frac{Cp}{R-s}\right)^{2/p}\norm{v}_{L^p(Q_R)}
\label{eq-pde-moser-step}
$$

を得る。$p_j=2\chi^j$、$R_j=s+(R-s)2^{-j}$ を使って反復すると、定数の積は有限である。実際、その対数は $\sum_j(2/p_j)(C+\log p_j+j\log2-\log(R-s))$ 以下である。従ってある $d>0$ について $\sup_{Q_s}v\le C(R-s)^{-d}\norm v_{L^2(Q_R)}$ となり、$v$ は局所有界である。

同じ上界は、定数とべき $d$ を変えれば任意の固定した $0<p_0<2$ にも使える。$M(R)=\operatorname*{ess\,sup}_{Q_R}v$ とし、$\norm v_2\le M(R)^{1-p_0/2}\norm v_{p_0}^{p_0/2}$ を使う。Young の不等式で任意の $\theta>0$ について $M(s)\le\theta M(R)+C_\theta (R-s)^{-2d/p_0}\norm v_{p_0}$ を得る。半径を $s$ から少し大きい半径まで増分 $(R-s)2^{-j-1}$ で増やし、$\theta\,2^{2d/p_0}<1$ と選ぶと幾何級数が収束する。最後の $M$ は局所有界なので $\theta^jM\to0$ である。従って固定した内外の立方体に対し $\sup_{\rm inner}v\le C\norm v_{p_0,\rm outer}$ が従う。
:::

:::proof[負のべきと対数の評価]
$v$ は $\delta$ 以上で局所有界である。従って $p<0$ では $\zeta^2v^{p-1}$ で直接試験できる。$v^{p-1}$ の微分の符号を反転して、絶対値を取った交差項に Young を使うと

$$
\int\zeta^2|D(v^{p/2})|^2
\le C\frac{p^2}{(p-1)^2}\int v^p|D\zeta|^2
$$

を得る。$v^{-1}$ に対して指数を $p_0\chi^j$ と増やし同じ Sobolev 反復を行えば、$\sup_{\rm inner}v^{-1}\le C\norm{v^{-1}}_{p_0,\rm outer}$ となる。$p-1$ は負の指数では零に近づかないため、定数の積の収束も正のべきと同じである。

残る問題は、$v$ と $v^{-1}$ の小さいべきの積である。$\zeta^2/v$ を試験すると $\int\zeta^2|D\log v|^2\le C\int |D\zeta|^2$ となる。固定した立方体 $Q$ の任意の部分立方体 $Q'$ で、同心で一回り大きい立方体が方程式の領域内に入るよう最初の $Q$ を選ぶ。この切断を用いると $\int_{Q'}|D\log v|^2\le C(\operatorname{side} Q')^{n-2}$ であり、Poincaré と Cauchy–Schwarz で $|Q'|^{-1}\int_{Q'}|\log v-(\log v)_{Q'}|\le C$ である。 :ref[thm-four-john-nirenberg]{full} の立方体上の証明は $Q$ 内の部分立方体だけを用いる。従ってこの局所的な平均振動の上界から、ある $p_0\in(0,2)$ に対し $|Q|^{-1}\int_Q\exp(p_0|\log v-(\log v)_Q|)\le C$ となる。従って $\big(\int_Qv^{p_0}\big)\big(\int_Qv^{-p_0}\big)\le |Q|^2C^2$ であり、平均の指数因子は相殺する。

内側の立方体で、正べきの上界と逆数の上界を掛けると $\sup v\le C\inf v$ を得る。$B_1$ を $B_2$ 内に十分な余白を持つ同じ大きさの小立方体で覆い、隣の内側立方体が重なるよう結べば、任意の二つの値を有限回の評価で比較できる。必要な個数は次元だけに依存する。$\delta\downarrow0$ とし、尺度を元に戻すと結論を得る。
:::

:::example[不連続係数でも正の値を比較する]
$(-1,1)$ で $a(x)=1$（左）、$2$（右）とし、$u(x)=2+x$（左）、$2+x/2$（右）とする。$au'=1$ なので弱い斉次方程式である。$(-1/2,1/2)$ の本質的上限と下限は $9/4,3/2$ で、その比は $3/2$ である。係数が跳んでも、二階微分を使わずに Harnack の評価を適用できる。
:::

### 有界な外力がある場合

同じ係数条件で $-\divg(A Du)=f$、$f\in L^\infty(B_{2r})$、$u\ge0$ とする。このとき

$$
\operatorname*{ess\,sup}_{B_r}u
\le C\left(\operatorname*{ess\,inf}_{B_r}u+\frac{r^2}{\lambda}\norm f_\infty\right).
$$

実際、半径を $1$ に規格化すると右辺は $r^2 f$ になる。$v=u+\delta+r^2\norm f_\infty/\lambda$ と置くと $r^2|f|\le\lambda v$ である。正べきの試験 $\zeta^2\psi(v)$ では、外力からの項は $\lambda\int\zeta^2v\psi(v)=\lambda\int\zeta^2\Phi(v)^2$ 以下である。従って先の反復で $|D\zeta|^2$ を $|D\zeta|^2+\zeta^2$ に替えれば足りる。負べき $p<0$ でも外力項は $\lambda\int\zeta^2v^p$ 以下である。$p-1$ は零から離れているため、反復定数が増えるのは高々 $p$ の多項式の範囲であり、その対数の級数は収束する。対数の試験では新しい項は $\lambda\int\zeta^2$ 以下なので、小立方体での平均振動の一様上界も保たれる。こうして $\sup v\le C\inf v$ を得て、$\delta\downarrow0$ とすれば示した式になる。例えば $u=x^2$ は $-\partial_{xx}u=-2$ を満たし、原点を含む区間で下限が零である。このとき上限の尺度 $r^2$ は外力の項が捉えている。

### 球ごとの比較を道に沿ってつなぐ

$u>0$ が方程式を満たすとし、弱解なら後節の振動縮小で得る連続代表を使う。$B_{2r}(x_j)\subset\Omega$ に同じ Harnack 定数 $C_H$ が使えるとする。$x_{j+1}\in B_r(x_j)$ なら $u(x_{j+1})\le C_Hu(x_j)$ であり、$N$ 段を掛けると $u(x_N)\le C_H^N u(x_0)$ となる。逆向きにも同じ被覆が使える場合は比の逆も評価できる。例えば境界から距離 $d$ 以上で長さ $L$ の曲線なら、$r=d/4$ として弧長を $r/2$ 以下に分ければ $N\le\lceil 8L/d\rceil+1$ とできる。従って領域に細い通路があると、この方法の定数は通路の幅に依存する。連結だけから一様な値の比が出るわけではない。

一次元では別の計算もできる。$(a u')'=0$、$0<\lambda\le a\le\Lambda$ なら $u(x)=c+q\int_0^x a^{-1}$ である。尺度変数 $s(x)=\int_0^x a^{-1}$ は増加し、$u$ はその変数で線形になる。$u\ge0$ を両端まで課せば、内側の点では線分補間係数の両方が正の下限を持つ。例えば $(-2,2)$ から $[-1,1]$ を見ると各係数は少なくとも $(1/\Lambda)/(4/\lambda)=\lambda/(4\Lambda)$ である。従って内側での最大・最小の比は $4\Lambda/\lambda$ 以下である。

:::example[Harnack 球を楕円体として戻す]
二変数の一定の拡散行列について、変数変換後の半径を比較する。

拡散行列を $A$ とする。対角成分は座標方向での値を与える。固有方程式から直交する方向を取り出すため、特性多項式 $\det(\lambda I-A)$ を計算する。

$$
\begin{aligned}
A&=\left[\begin{matrix}\frac{53}{5} & - \frac{36}{5}\\- \frac{36}{5} & \frac{32}{5}\end{matrix}\right]\\
\det(\lambda I-A)&=\lambda^{2} - 17 \lambda + 16=\left(\lambda - 16\right) \left(\lambda - 1\right)\\
d_1&=1,\quad v_1=\left[\begin{matrix}\frac{3}{5}\\\frac{4}{5}\end{matrix}\right]\\
d_2&=16,\quad v_2=\left[\begin{matrix}- \frac{4}{5}\\\frac{3}{5}\end{matrix}\right]
\end{aligned}
$$

$V$ の列に固有ベクトルを並べる。$V^\top V=I$ と $V^\top AV=D$ を使えば、元の座標 $X$ と固有座標 $Y=V^\top X$ の双方で勾配を記述できる。二次形式は固有座標にすると交差項を持たない。

$$
\begin{aligned}
V&=\left[\begin{matrix}\frac{3}{5} & - \frac{4}{5}\\\frac{4}{5} & \frac{3}{5}\end{matrix}\right]\\
V^\top AV&=\left[\begin{matrix}1 & 0\\0 & 16\end{matrix}\right]\\
X^\top AX&=\frac{53 x_{1}^{2} - 72 x_{1} x_{2} + 32 x_{2}^{2}}{5}\\
\sum_j d_jY_j^2&=Y_{1}^{2} + 16 Y_{2}^{2}
\end{aligned}
$$

一つの勾配で固有係数を計算し、勾配の二次形式を求める。直交変換なので二乗ノルムは変わらない。二次形式を二乗ノルムで割ると、固有値を非負の重みで平均した値になる。

$$
\begin{aligned}
X&=\left[\begin{matrix}1\\1\end{matrix}\right]\\
Y&=\left[\begin{matrix}\frac{7}{5}\\- \frac{1}{5}\end{matrix}\right]\\
|X|^2&=2\\
X^\top AX&=\frac{13}{5}\\
\frac{X^\top AX}{|X|^2}&=\frac{13}{10}\\
\min_j d_j&=1\\
\max_j d_j&=16
\end{aligned}
$$

$x=A^{1/2}y$ で $A:D^2$ が Laplacian になる。$y$ の単位球は $x$ 座標で半軸 $1,4$ の楕円であり、半径 $1/2$ の球の像では半軸も半分になる。元の Euclid 球で評価する場合、この軸比は被覆の数に影響する。
:::

### Harnack の定数が境界まで保たれない例

上半平面で $u(x,y)=y$ は非負調和で、$y=0$ で消える。中心 $(0,1)$ の半径 $0<r<1$ の球では最大が $1+r$、最小が $1-r$（開球では上限と下限）であり、比は $(1+r)/(1-r)$ である。例えば $r=1/2$ では $3$ だが、$r\uparrow1$ では比は発散する。一方、中心を $(0,h)$、半径を $h/2$ とすれば同じ比 $3$ である。Harnack 不等式の一様性は、比較する球の大きさと境界からの距離を同じ割合で保てるときに使える。

## Schauder 評価

:term[Schauder 評価]{en="Schauder estimate" yomi="しゃうだーひょうか"} は、右辺や係数の Hölder 連続性を用いて、解の二階微分の Hölder 半ノルムまで抑える評価である。まず $-\Delta u=f$ の内部評価を得る。$B_2$ での方程式を考え、$B_{3/2}$ で $1$ の切断 $\chi\in C_c^\infty(B_2)$ を選ぶ。$v=N(\chi f)$ とすれば :ref[prop-pde-newton-holder] により $v$ の $C^{2,\alpha}$ ノルムは $C\norm f_{C^{0,\alpha}(B_2)}$ 以下である。$u-v$ は $B_{3/2}$ で調和であり、:ref[prop-pde-harmonic-derivative] の $m\le3$ の評価により $B_1$ の $C^{2,\alpha}$ ノルムが $C\norm{u-v}_\infty$ 以下である。従って

$$
\norm u_{C^{2,\alpha}(B_1)}
\le C\big(\norm u_{C^0(B_2)}+\norm f_{C^{0,\alpha}(B_2)}\big).
$$

:::theorem[変数係数の内部 Schauder 評価]{#thm-pde-interior-schauder}
$P=A(x):D^2+b(x)\cdot D+c(x)$ は一様楕円、係数は $B_2$ 上で $C^{0,\alpha}$ とする。$u\in C^{2,\alpha}_{\rm loc}(B_2)$、$Pu=f$ のとき、

$$
\norm{u}_{C^{2,\alpha}(B_1)}\le C(\norm u_{C^0(B_2)}+\norm f_{C^{0,\alpha}(B_2)})
\label{eq-pde-schauder-interior}
$$

である。定数は次元、指数、楕円性、係数の Hölder ノルムに依存する。
:::

:::proof[境界までの距離を重みにする]
$D=B_{3/2}$ を使い、$d_x=\operatorname{dist}(x,\partial D)$、$d_{xy}=\min(d_x,d_y)$ とする。次の量は $u$ が $\overline D$ の近傍で $C^{2,\alpha}$ なので有限である。

$$
S=\sup_D d_x^2|D^2u(x)|,\qquad G=\sup_D d_x|Du(x)|,\qquad
K=\sup_{x\ne y\in D}d_{xy}^{2+\alpha}\frac{|D^2u(x)-D^2u(y)|}{|x-y|^\alpha}.
$$

単位方向 $e$ に沿い、刻み $h=\delta d_x$ の二次中心差分を取る。$0<\delta<1/4$ なら両側の線分で $d_y\ge d_x/2$ であり、二次の Taylor 余りは $C K h^{2+\alpha}d_x^{-2-\alpha}$ 以下である。二次差分そのものは $4\norm u_\infty$ 以下なので、$h^2$ で割り $d_x^2$ を掛けると

$$
S\le C\delta^{-2}\norm u_\infty+C\delta^\alpha K .
$$

全方向の二階微分の評価から Hessian のノルムを抑えた。一階差分を同じように用いると $G\le C\delta^{-1}\norm u_\infty+C\delta S$ である。特に任意の $\varepsilon>0$ に対し $S+G\le\varepsilon K+C_\varepsilon\norm u_\infty$ である。
:::

:::proof[係数を固定し最高階の項を吸収する]
二点 $x,y$ で $d_x=d_{xy}$ とする。固定した小さい $\theta>0$ に対し $|x-y|\ge\theta d_x/2$ なら $K$ のこの組の商は $C\theta^{-\alpha}S$ 以下である。残る場合は $r=\theta d_x$ として $B_r(x)$ を取り、$A(x)$ を固定する。この球では全ての距離の重みが $d_x$ と同程度である。定数係数の評価は楕円行列の平方根で Laplace 方程式に変換し、同程度の球を選べば得られる。尺度を元に戻すと、$A(x):D^2u=h$ に対して

$$
[D^2u]_{\alpha,B_{r/2}(x)}
\le C(r^{-2-\alpha}\norm u_{\infty,B_r}
+r^{-\alpha}\norm h_{\infty,B_r}+[h]_{\alpha,B_r})
$$

となる。ここで $h=f+(A(x)-A(\cdot)):D^2u-b\cdot Du-cu$ である。

$h$ をこの評価に代入し $d_x^{2+\alpha}$ を掛ける。$A$ の差の大きさは $[A]_\alpha r^\alpha$ 以下なので、$[D^2u]_{\alpha,B_r}$ を含む唯一の項は $C[A]_\alpha r^\alpha K$ 以下である。積の半ノルムの他方の項は $[A]_\alpha\sup|D^2u|$ であり、重み付きでは $CS$ 以下になる。$r^{-\alpha}\norm h_\infty$ からも主部の差は $CS$ 以下である。下位項では $[Du]_{\alpha,B_r}\le(2r)^{1-\alpha}\norm{D^2u}_{\infty,B_r}$ と $[u]_{\alpha,B_r}\le(2r)^{1-\alpha}\norm{Du}_{\infty,B_r}$ を使えば、$S+G+\norm u_\infty$ に定数を掛けたもの以下になる。従って遠い組も合わせて

$$
K\le C[A]_\alpha(2\theta)^\alpha K+C_\theta(S+G+\norm u_\infty+\norm f_{C^{0,\alpha}(B_2)})
$$

を得る。最初に $\theta$ を小さくして第 1 項の係数を $1/4$ 以下とし、その後で補間の $\varepsilon$ を選んで $C_\theta(S+G)$ からの $K$ の係数も $1/4$ 以下にする。残る評価で $B_1$ では $d_x\ge1/2$ なので全ての重みを外せる。
:::

### Dirichlet 境界での Schauder 評価

まず平面 $x_n=0$ の上側で $\Delta u=f$、$u=0$ とする。右辺の平面上の値が零なら、奇延長で内部評価に帰着する。実際、接方向の二階微分は境界値から零、法線二階微分も $\Delta u=f=0$ から零であり、奇延長は二階まで接続する。混合微分は偶に延長されるので接続を妨げない。右辺の奇延長は平面の両側の点に対しても Hölder 連続である。二点を結ぶ線分と平面との交点で値が消えるため、両側からその点までの距離の $\alpha$ 乗を加えればよい。

$f$ が平面上で消えない場合は、境界値を壊さずに右辺だけを修正する。平面上の $f_0(x')=f(x',0)$ を切断し、平滑化核 $\rho_t(z)=t^{-(n-1)}\rho(z/t)$ を用いて $F(x',t)=(\rho_t*f_0)(x')$、$t=x_n>0$ と置く。$\int\rho=1$ とする。$\psi=t^2 F/2$ とすれば、境界で $\psi=0$ かつ $\Delta\psi=f_0$ を満たす。以下で正則性を確認する。$F-f_0$ の一様ノルムは $C t^\alpha[f_0]_\alpha$ 以下であり、全ての $k\ge1$ の $x',t$ に関する導関数は $Ct^{\alpha-k}[f_0]_\alpha$ 以下である。実際、微分した核の積分は零であり、差 $f_0(x'-z)-f_0(x')$ を使った積分の範囲は $|z|\le Ct$ である。従って $D^2\psi$ は、法線法線成分が $f_0$、他成分が零の境界値との差を $Ct^\alpha[f_0]_\alpha$ で抑えられる。また $|D^3\psi|\le Ct^{\alpha-1}[f_0]_\alpha$ である。二点の距離が小さい側の高さの半分未満なら、線分上で最後の勾配評価を使い $C|x-y|^\alpha$ で抑える。それ以外なら両点の高さは距離の定数倍以下なので、各点とその平面への射影の差、そして平面上での差を合わせて同じ Hölder 上界を得る。よって $\norm{\psi}_{C^{2,\alpha}}\le C\norm{f_0}_{C^{0,\alpha}}$ である。$u-\psi$ に奇延長を使えば、単位半球の中央の半球で内部と同じ形の評価を得る。切断は中央のさらに広い範囲で $1$ にしておけばよい。

凍結した正定値係数 $A_0$ の場合も、平面を保つ線形変換で Laplace 作用素に直せる。$A_0=\begin{pmatrix}A'&b\\b^\top&a\end{pmatrix}$ と書けば、$y'=x'-b x_n/a$、$y_n=x_n$ で混合項が消え、対角ブロックは $A'-bb^\top/a$ と $a$ になる。各ブロックの平方根で規格化すればよい。楕円性と行列の上界から、この変換と逆変換のノルムは一様に抑えられる。従って適切な内側の半球で、尺度 $r$ に対する定数係数評価が得られる。

:::proposition[境界までの Schauder の事前評価]{#prop-pde-boundary-schauder}
$\Omega$ は有界 $C^{2,\alpha}$ 領域、$P$ は :ref[thm-pde-interior-schauder] と同じ仮定を閉包上で満たすとする。$u\in C^{2,\alpha}(\overline\Omega)$ が $Pu=f$、$u|_{\partial\Omega}=g$ を満たすとき

$$
\norm u_{C^{2,\alpha}(\overline\Omega)}
\le C\big(\norm u_{C^0(\overline\Omega)}+\norm f_{C^{0,\alpha}(\overline\Omega)}
+\norm g_{C^{2,\alpha}(\partial\Omega)}\big).
$$

:::

:::proof
境界を有限個の $C^{2,\alpha}$ 座標で平坦化する。偏微分を二度変換すると、主部の係数の Hölder 性は保たれ、座標の二階微分は一階項に入る。各平坦な座標で、境界データ $g(x')$ を法線方向に一定に延ばし切断を掛ける。この延長の $C^{2,\alpha}$ ノルムは元の境界ノルム以下の定数倍である。これを引けば零境界の局所問題になり、右辺にも同じ定数倍の変化しか生じない。

平坦な半球状パッチでは、$d_x$ を曲面側の端までの距離とし、前の重み $S,G,K$ を用いる。補間は二次中心差分の代わりに、$e_n$ 成分が非負の単位方向 $e$ に沿う前進二次差分 $u(x+2he)-2u(x+he)+u(x)$ を使う。どの線分も半空間の内部に留まり、曲面側の距離は $h$ を十分小さく選べば半分以上である。二次形式の値は $e$ と $-e$ で同じなので、この半球の方向だけで Hessian を評価できる。一階微分の各成分は $e=e_j$ を使った前進差分で評価できる。従って補間式はそのまま成立する。

近い二点の組で、平面までの距離が $d_x$ の小さい固定倍より大きければ内部球で係数を凍結する。平面に近い場合は $x$ の平面への射影を中心とする半球で凍結する。具体的に $|x-y|<\theta d_x$、$x_n\le2\theta d_x$ なら半径 $8\theta d_x$ の平坦な半球を選ぶ。内側半球に両点が入り、外側半球の中でも曲面端までの距離は $d_x$ と同程度である。$x_n>2\theta d_x$ なら半径 $2\theta d_x$ の内部球を使う。いずれも係数の差は $C[A]_\alpha(\theta d_x)^\alpha$ 以下なので、内部評価と同じ積の評価で $C[A]_\alpha\theta^\alpha K$ を吸収できる。遠い二点は $S$ で抑えられ、残る $S,G$ は補間で消せる。これでパッチの中央で評価が得られる。

有限個の境界パッチと内球を合わせる。二点が十分近ければ同じ拡大パッチの内側に入り、その局所評価を使える。十分離れた二点の半ノルムは Hessian の一様ノルムで抑えられる。従って閉包全体で所望の評価が成り立つ。
:::

この命題で得た事前評価を、後で連続法による古典解の存在証明に使う。一次元の場合、境界の平坦化に代えて $u''=(f-bu'-cu)/a$ と書けば、$u,u'$ の一様ノルムから同じ Hölder 評価が得られる。実際、短い線分での一階差分から $\norm{u'}_\infty\le\varepsilon\norm{u''}_\infty+C_\varepsilon\norm u_\infty$ として低次項を吸収し、$u'$ の半ノルムを $\norm{u''}_\infty$ で抑える。

:::example[Hölder 指数を上げられない例]
$u(x)=|x|^{2+\alpha}$ を一次元で取ると $u''=(2+\alpha)(1+\alpha)|x|^\alpha$ である。$u$ は $C^{2,\alpha}$ だが $\beta>\alpha$ では $u''(x)-u''(0)$ を $|x|^\beta$ で割った比が発散する。
:::

:::example[球の縮小と Hessian の尺度]
$-\Delta u=f$ を $B_{2r}(x_0)$ で考え、$v(y)=u(x_0+ry)$ とする。$-\Delta_yv=r^2 f(x_0+ry)$ であり、$[D_y^2v]_\alpha=r^{2+\alpha}[D_x^2u]_\alpha$ である。従って単位球の評価から $r^{2+\alpha}[D^2u]_{\alpha,B_r}\le C(\norm u_{\infty,B_{2r}}+r^2\norm f_{\infty,B_{2r}}+r^{2+\alpha}[f]_{\alpha,B_{2r}})$ を得る。右辺のノルムを同じ尺度に直してから比べる必要がある。
:::

### 連続法で事前評価を可解性に使う

Banach 空間 $X,Y$ と有界な作用素 $L_s=(1-s)L_0+sL_1:X\to Y$、$0\le s\le1$ を考える。$L_0$ は同型、全ての $s$ に共通の $\norm u_X\le C\norm{L_su}_Y$ があるとする。$L_s$ が同型なら $L_{s+h}=L_s(I+hL_s^{-1}(L_1-L_0))$ である。$|h|C\norm{L_1-L_0}<1$ のとき括弧の逆は Neumann 級数で得られるので、共通の幅で可逆性を前へ進められる。この幅は $s$ に依存しないため、有限回で $s=1$ に達する。従って連続法に必要なものは連続な係数変化、共通の事前評価、始点の可解性である。

古典的 Dirichlet 問題では $X=\{u\in C^{2,\alpha}(\overline\Omega):u|_{\partial\Omega}=0\}$、$Y=C^{0,\alpha}(\overline\Omega)$ を使う。Schauder 評価に残る $\norm u_\infty$ を比較原理から $\norm{L_su}_\infty$ で押さえれば、上の形になる。例えば有界領域の座標 $x_1$ 方向に $a_{11}\ge\lambda$、$|b_1|\le B$、$c\le0$ なら $w=C_0-e^{kx_1}\ge0$ を取れる。$P w\le-(\lambda k^2-Bk)e^{kx_1}$ なので $k>B/\lambda$ で厳密な負の上界が得られ、$Pw$ の符号で $u$ と $-u$ を挟める。下位項の符号と境界条件を押さえることが、事前評価を同型の評価へ変える段階である。

:::theorem[古典的 Dirichlet 問題の可解性]{#thm-pde-schauder-dirichlet}
:ref[prop-pde-boundary-schauder] と同じ領域と係数を仮定し、$c\le0$ とする。任意の $f\in C^{0,\alpha}(\overline\Omega)$ と $g\in C^{2,\alpha}(\partial\Omega)$ に対して、$Pu=f$、$u|_{\partial\Omega}=g$ は一意な $C^{2,\alpha}(\overline\Omega)$ 解を持つ。
:::

:::proof
境界の座標上で $g$ を法線方向に一定に延長し、切断と分割の一でつなぐ。延長を引いて零境界値に帰着する。$P_s=(1-s)\Delta+sP$ では楕円定数の下界、係数の Hölder ノルムの上界、$c_s\le0$ が一様である。上の障壁を $P_s$ に共通に使い、$P_s w\le-1$ と規格化する。比較を $u-\norm{P_su}_\infty w$ と $-u-\norm{P_su}_\infty w$ に用いれば $\norm u_\infty\le C\norm{P_su}_\infty$ である。境界 Schauder 評価で $\norm u_X\le C\norm{P_su}_Y$ となるので、始点の Laplace 作用素の可解性を示せばよい。

境界が滑らかで $f$ も滑らかなら、まず :ref[thm-pde-elliptic-coercive] で弱解を得る。上で証明した高階 Dirichlet 正則性と Sobolev 埋め込みにより、この解は閉包まで $C^{2,\alpha}$ である。一般の領域では、境界を内側から滑らかな領域 $\Omega_j$ で近似する。具体的には境界の近くで、各座標のグラフの上下差を内部が負になる符号で取り、分割の一でつないで定義関数 $\rho$ を作る。境界上では各局所関数は零であり、各勾配は外法線方向に正の倍数なので $D\rho\ne0$ である。境界から離れた側で正負の定数に滑らかにつなげられる。$\rho$ を平滑化した $\rho_j$ に対し、一様近似誤差より大きい $\delta_j\downarrow0$ を使って $\Omega_j=\{\rho_j<-\delta_j\}$ とすれば $\Omega_j\subset\Omega$ で、$\Omega$ の任意のコンパクト部分は十分大きい $j$ で $\Omega_j$ に入る。

$\rho_j$ の $C^{2,\alpha}$ ノルムは一様に有界で、$\rho_j\to\rho$ が $C^2$ で成り立つ。境界近くで一つの座標方向の微分を零から離しておけば、陰関数の微分公式により各グラフの $C^{2,\alpha}$ ノルムも一様に有界である。従って固定した有限個のパッチと共通の縮小幅で上の境界評価が使える。$f$ を閉包の外へ延長するときは、例えば $f_{\rm ext}(x)=\inf_{y\in\overline\Omega}(f(y)+[f]_\alpha |x-y|^\alpha)$ とすればよい。$|x-y|^\alpha$ の三角不等式で元の値と半ノルムを保つ。切断して平滑化した $f_j$ は一様に $C^{0,\alpha}$ 有界で $f$ に一様収束する。

$\Delta u_j=f_j$ を $\Omega_j$ 上の零境界値で解く。共通の障壁と境界評価から $u_j$ の $C^{2,\alpha}(\overline{\Omega_j})$ ノルムは一様有界である。コンパクト集合上で Ascoli の定理と対角化を使うと、二階まで局所一様に収束する部分列を得る。極限は $\Delta u=f$ を満たす。また勾配の一様上界により $|u_j(x)|\le C\operatorname{dist}(x,\partial\Omega_j)\le C\operatorname{dist}(x,\partial\Omega)$ なので、極限の境界値は零である。二階微分の Hölder 評価は任意の内点の組で極限に移せる。境界の近くではグラフ座標で二点を距離の定数倍の長さの曲線で結べるため、一階微分も連続に閉包へ延長される。これで $u\in X$ を得た。以上で $P_0$ は同型となり、連続法で $P_1=P$ の結論が従う。
:::

### $L^p$ データと二階微分

Hölder 連続性を仮定できない右辺には、二階微分を $L^p$ で測る方法がある。まず $\R^n$ で $v\in C_c^\infty$、$1<p<\infty$ とする。Fourier 乗数を比較すると $\partial_i\partial_j v=R_iR_j(-\Delta v)$ である。実際、$R_iR_j$ の乗数は $-\xi_i\xi_j/|\xi|^2$ で、$-\Delta$ の乗数は $4\pi^2|\xi|^2$ である。:ref[thm-four-riesz-transform] により

$$
\norm{D^2 v}_{L^p(\R^n)}\le C_{n,p}\norm{\Delta v}_{L^p(\R^n)} .
$$

定数正定値対称行列に対しては行列の平方根で座標を変えればよい。この評価を局所化するには、$v=\zeta u$ に適用し、積の法則で生じる項を分ける。例えば $P=A(x):D^2$ で $A$ が一様楕円かつ連続なら、小球上で $A_0=A(x_0)$ と凍結する。$v$ の台がこの球に入る場合、

$$
\norm{D^2v}_p\le C\big(\norm{Pv}_p+\norm{A-A_0}_{\infty}\norm{D^2v}_p\big)
$$

となる。半径を小さくして最後の係数を $1/2$ 以下にすると吸収できる。従って $U\Subset V\Subset\Omega$ で、有限個の小球と切断を使い $\norm{D^2u}_{L^p(U)}\le C(\norm{Pu}_{L^p(V)}+\norm u_{W^{1,p}(V)})$ を得る。これはまず $W^{2,p}$ の候補に対する事前評価である。係数が有界可測なだけでは、小球上の振幅を一様に小さくする操作は使えない。

有界 $C^{1,1}$ 領域で $A$ が閉包上で連続、下位係数が有界とする。零 Dirichlet 条件まで評価を広げるには、境界を平面に直した後に同じ凍結を行う。凍結した対称行列を、Schauder 評価で使った平面を保つ変換で Laplace 作用素にする。平面上で零の $W^{2,p}$ 関数は、奇反射すると全空間で $W^{2,p}$ に属する。実際、接微分の跡は零、法線一階微分の反射は偶であり、一階微分の跡が両側で一致するので二階微分に面上の測度は生じない。右辺も奇反射される。従って切断が平面以外の境界に触れない関数について、全空間の評価で半空間の Hessian を抑えられる。平坦化で生じる下位係数は有界なので、有限個の小パッチで凍結と吸収を行うと

$$
\norm{D^2u}_{L^p(\Omega)}
\le C\big(\norm{Pu}_{L^p(\Omega)}+\norm{u}_{W^{1,p}(\Omega)}\big)
$$

を得る。$\norm{Du}_p\le\varepsilon\norm{D^2u}_p+C_\varepsilon\norm u_p$ を使って最後の一階項を吸収できる。この補間不等式が破れるなら、$\norm{Du_j}_p=1$、$\norm{D^2u_j}_p\le1/\varepsilon$、$\norm{u_j}_p\to0$ となる列がある。Rellich の定理で $W^{1,p}$ の強極限を取ると極限は零であり、勾配ノルムが $1$ であることに反する。従って $W^{2,p}$ ノルムは $\norm{Pu}_p+\norm u_p$ で抑えられる。一般の $L^p$ データでの可解性には、候補の正則性を仮定せずに解を構成する段階が別に必要である。係数の条件と可解性の一般論は :cite[gilbarg-trudinger] を参照する。

:::example[係数を凍結した行列の楕円性]
ある点で係数を凍結し、方向別の値と上限・下限を算出する。

拡散行列を $A$ とする。対角成分は座標方向での値を与える。固有方程式から直交する方向を取り出すため、特性多項式 $\det(\lambda I-A)$ を計算する。

$$
\begin{aligned}
A&=\left[\begin{matrix}\frac{31}{5} & - \frac{12}{5}\\- \frac{12}{5} & \frac{24}{5}\end{matrix}\right]\\
\det(\lambda I-A)&=\lambda^{2} - 11 \lambda + 24=\left(\lambda - 8\right) \left(\lambda - 3\right)\\
d_1&=3,\quad v_1=\left[\begin{matrix}\frac{3}{5}\\\frac{4}{5}\end{matrix}\right]\\
d_2&=8,\quad v_2=\left[\begin{matrix}- \frac{4}{5}\\\frac{3}{5}\end{matrix}\right]
\end{aligned}
$$

$V$ の列に固有ベクトルを並べる。$V^\top V=I$ と $V^\top AV=D$ を使えば、元の座標 $X$ と固有座標 $Y=V^\top X$ の双方で勾配を記述できる。二次形式は固有座標にすると交差項を持たない。

$$
\begin{aligned}
V&=\left[\begin{matrix}\frac{3}{5} & - \frac{4}{5}\\\frac{4}{5} & \frac{3}{5}\end{matrix}\right]\\
V^\top AV&=\left[\begin{matrix}3 & 0\\0 & 8\end{matrix}\right]\\
X^\top AX&=\frac{31 x_{1}^{2} - 24 x_{1} x_{2} + 24 x_{2}^{2}}{5}\\
\sum_j d_jY_j^2&=3 Y_{1}^{2} + 8 Y_{2}^{2}
\end{aligned}
$$

一つの勾配で固有係数を計算し、勾配の二次形式を求める。直交変換なので二乗ノルムは変わらない。二次形式を二乗ノルムで割ると、固有値を非負の重みで平均した値になる。

$$
\begin{aligned}
X&=\left[\begin{matrix}2\\3\end{matrix}\right]\\
Y&=\left[\begin{matrix}\frac{18}{5}\\\frac{1}{5}\end{matrix}\right]\\
|X|^2&=13\\
X^\top AX&=\frac{196}{5}\\
\frac{X^\top AX}{|X|^2}&=\frac{196}{65}\\
\min_j d_j&=3\\
\max_j d_j&=8
\end{aligned}
$$

作用素ノルムでの摂動が $1$ 以下なら、新しい行列の二次形式は単位ベクトル上で $2$ 以上 $9$ 以下になる。係数の変動が十分小さい球を選び、凍結した係数との差を最高階微分の小さい係数として左辺へ吸収する。
:::

### 補間で低い微分の項を吸収する

一次元で $u$ は $[x-h,x+h]$ 上で $C^2$ とする。Taylor 公式から
$$
u(x+h)-u(x-h)=2h u'(x)
+\int_0^h(h-s)(u''(x+s)-u''(x-s))\dd s .
$$
従って $|u'(x)|\le h^{-1}\norm u_\infty+(h/2)\norm{u''}_\infty$ である。内部評価で右辺に $C\norm{u'}_\infty$ が出た場合、$h$ を十分小さくすれば $C h\norm{u''}_\infty/2$ を左へ吸収し、残りは $C/h$ と $\norm u_\infty$ で測れる。境界に近い点では幅 $h$ の対称差を取れないため、内部だけに適用するか片側差と境界の条件を用いる。例えば $u(x)=\cos(kx)$ では $\norm{u'}_\infty=k$、$\norm{u''}_\infty=k^2$ であり、$h=1/k$ とすれば右辺は $3k/2$。二つの項を同じ尺度にそろえることで、高周波での正しい次数が得られる。区間を半分に縮めるときに $h$ も半分にすれば、この評価での二階微分の係数も半分になる。従って縮小する内球では、吸収する項と残る定数を半径とともに記録する必要がある。

## De Giorgi–Nash–Moser の理論（概観）

一様楕円な有界可測係数 $A$ に対する $-\divg(A Du)=0$ の $H^1_{\rm loc}$ 弱解は局所 Hölder 連続な代表を持つ。指数と内部評価の定数は次元と楕円性の比に依存する。この :term[De Giorgi–Nash–Moser 理論]{en="De Giorgi–Nash–Moser theory" yomi="でじょるじなっしゅもーざーりろん"} は係数の微分を用いない。一方で、先の界面例が示すように二階弱微分まで一般には得られない。

:::proposition[切断レベルのエネルギー評価]{#prop-pde-caccioppoli}
$w=(u-k)_+$、$\zeta\in C_c^\infty(\Omega)$ とすると

$$
\int \zeta^2|Dw|^2\le C\int w^2|D\zeta|^2
\label{eq-pde-caccioppoli}
$$

である。
:::

:::proof
$\zeta^2 w$ で試験する。$w>0$ の範囲で $Du=Dw$ であり、左辺の主部は $\lambda\int\zeta^2|Dw|^2$ 以上である。交差項は $2\norm A_\infty\int|\zeta Dw|\,|wD\zeta|$ 以下で、Young の不等式により主部の半分を吸収して結論を得る。
:::

半径 $r$ の球で $\zeta=1$、半径 $2r$ で $\zeta=0$ とすれば右辺は $Cr^{-2}\int_{B_{2r}}(u-k)_+^2$ である。小さい球と高い切断レベルを順に選び、Sobolev 埋め込みで高いレベル集合の測度を抑えるのが De Giorgi の反復の出発点である。別の帰結として正の解に Harnack を適用すると振動が減る。$M,m$ を $B_{2r}$ 上の本質的上下限、$M_r,m_r$ を $B_r$ 上のものとする。非負解 $M-u,u-m$ についての Harnack は $M-m_r\le C(M-M_r)$ と $M_r-m\le C(m_r-m)$ を与える。足して整理すると $M_r-m_r\le(C-1)(M-m)/(C+1)$ である。これを半径ごとに反復すれば指数 $\alpha=-\log_2((C-1)/(C+1))$ の振動評価を得る。$C>1$ に取り、必要なら指数を $1$ より小さく選べる。符号を持つ $u$ でも、$u_+$ と $(-u)_+$ に正べきの評価を使うと局所有界性が従う。実際、正部分の凸な滑らかな近似 $\psi$ に対し、非負の $\phi$ を用いて $\psi'(u)\phi$ で元の方程式を試験すると $\int A D\psi(u)\cdot D\phi=-\int\phi\psi''(u) A Du\cdot Du\le0$ になる。極限で両正部分が非負な弱劣解であり、正べきの試験関数は非負なので同じ上界を得られる。従ってこの振動評価は一般の弱解にも適用できる。半径が零に向かう球の本質的上下限が同じ極限を持つため、その値を各点での代表とすれば局所 Hölder 連続である。例えば共通の大球の内側で $x,y$ の距離を $d$ とすると、$x$ を中心とした半径 $2d$ の球上の振動が $C d^\alpha$ 以下となり、$x,y$ での値の差も同じ評価を満たす。

:::example[切断するレベルを超える集合]
$u=x_1$ を二次元の単位球で考え、$k=1/2$ とする。$w=(x_1-1/2)_+$ に対し $|Dw|^2$ は $x_1>1/2$ の部分で $1$、残りで $0$ である。従ってエネルギーは全円盤ではなく円弧で囲まれた部分の面積 $\int_{1/2}^1 2\sqrt{1-s^2}\,\dd s=\pi/3-\sqrt3/4$ である。切断レベルを上げると評価する領域が小さくなる。
:::

:::example[振動の幾何級数]
もし $B_2$ の振動が $8$ で、一回の縮小での Harnack 定数が $C=3$ なら係数は $1/2$ である。半径を三回半分にした $B_{1/4}$ での振動は $1$ 以下となる。解の値がどこで最大かを求めなくても、局所振動に上界を与えられる。
:::

### 振動の縮小から Hölder 代表へ

ある $\theta\in(0,1)$ について、解の本質的振動が $\operatorname{osc}_{B_{r/2}}u\le\theta\operatorname{osc}_{B_r}u$ を満たすとする。反復すると半径 $2^{-k}R$ で $\theta^k$ 倍以下である。$2^{-k-1}R<r\le2^{-k}R$ と選べば
$$
\operatorname{osc}_{B_r}u\le 2^\alpha(r/R)^\alpha\operatorname{osc}_{B_R}u,qquad
\alpha=-\log\theta/\log2.
$$
この式を中心を動かして使う。任意の二つの Lebesgue 点 $x,y$ を距離 $2|x-y|$ の球に入れれば、差は右の振動で抑えられる。Lebesgue 点は稠密なので一意な連続拡張が得られ、それが Hölder 代表である。

斉次方程式なら Harnack からも縮小が得られる。$m,M$ を外球での本質的下限・上限、$m',M'$ を内球での値とする。$u-m$ と $M-u$ に Harnack を使って $M'-m\le C(m'-m)$、$M-m'\le C(M-M')$ を足す。すると $(C+1)(M'-m')\le(C-1)(M-m)$ なので $\theta=(C-1)/(C+1)$ である。非斉次の場合は外力に応じた項を残して反復する必要がある。

:::example[中心での振動を半径とともに追う]
単位円板の調和関数 $u=r^2\cos2\theta$ を使い、振動縮小の指数がこの例でどう現れるか調べる。

円周データを定数と三角関数に分ける。角振動数 $j$ について半径の方程式は $R^{\prime\prime}+r^{-1}R^{\prime}-j^2r^{-2}R=0$ である。元の条件に合う根を選び、境界半径で係数を規格化する。この例での関数は次のとおりである。

$$
\begin{aligned}
u(r,\theta)&=r^{2} \cos{\left(2 \theta \right)}\\
g&=\cos{\left(2 \theta \right)},\qquad r=1\\
u_r&=\frac{\partial}{\partial r} r^{2} \cos{\left(2 \theta \right)}=2 r \cos{\left(2 \theta \right)}\\
u_\theta&=\frac{\partial}{\partial \theta} r^{2} \cos{\left(2 \theta \right)}=- 2 r^{2} \sin{\left(2 \theta \right)}
\end{aligned}
$$

二階微分を別々に計算する。半径方向の二階微分だけを角度方向の二階微分と足しても相殺しない。極座標の Laplacian の一階項も入れると各振動数で相殺する。

$$
\begin{aligned}
u_{rr}&=\frac{\partial^{2}}{\partial r^{2}} r^{2} \cos{\left(2 \theta \right)}=2 \cos{\left(2 \theta \right)}\\
u_{\theta\theta}&=\frac{\partial^{2}}{\partial \theta^{2}} r^{2} \cos{\left(2 \theta \right)}=- 4 r^{2} \cos{\left(2 \theta \right)}\\
\Delta u&=0
\end{aligned}
$$

境界における値の平均と法線微分を読む。$j\ne0$ の角方向の平均は零であり、全体の平均は定数項である。またエネルギーは体積積分の方法と Green の境界積分の方法の両方で計算できる。

$$
\begin{aligned}
\text{円周の平均}&=\int\limits_{- \pi}^{\pi} \frac{\cos{\left(2 \theta \right)}}{2 \pi}\, d\theta=0\\
\partial_\nu u&=2 \cos{\left(2 \theta \right)}\\
\int_{-\pi}^{\pi}|Du|^2\,\mathrm{d}\theta&=8 \pi r^{2}\\
\int |Du|^2&=\int\limits_{0}^{1} 8 \pi r^{3}\, dr=2 \pi\\
\int_{\partial\Omega} u\partial_\nu u&=\int\limits_{- \pi}^{\pi} 2 \cos^{2}{\left(2 \theta \right)}\, d\theta=2 \pi
\end{aligned}
$$

次に同じ振動数の半径因子を比較し、角度を固定した一点で値と半径方向の傾きを調べる。角度の取り方によって、振幅のあるモードでも点での値が零になることがある。

$$
\begin{aligned}
u&=- \frac{1}{8},\qquad r=\frac{1}{2},\theta=\frac{\pi}{3}\\
u_r&=- \frac{1}{2},\qquad r=\frac{1}{2},\theta=\frac{\pi}{3}\\
u_\theta&=- \frac{\sqrt{3}}{4},\qquad r=\frac{1}{2},\theta=\frac{\pi}{3}
\end{aligned}
$$

各 $r$ の円板で振動は $2r^2$ なので、中心を固定して半径を半分にすると振動は $1/4$ になる。ただし一般の点での一次微分は非零なので、領域のどこでも振動が $r^2$ で減るわけではない。一般定理では全ての内点で使える共通の指数を求める。
:::

### レベル集合を減らして局所有界性を得る

$n\ge3$、$u\in H^1(B_1)$ を一様楕円な $-\divg(A Du)=0$ の弱解とする。$\norm{u_+}_2=0$ のときは既に $u\le0$ なので、以下では $\norm{u_+}_2>0$ とする。半径を $r_j=1/2+2^{-j-1}$、レベルを $k_j=K(1-2^{-j})$、$w_j=(u-k_j)_+$、$Y_j=\int_{B_{r_j}}w_j^2$ と置く。$B_{r_{j+1}}$ で $1$、$B_{r_j}$ の外で $0$ の切断 $\zeta_j$ は $|D\zeta_j|\le C2^j$ とできる。切断レベルの評価から $\int|D(\zeta_j w_j)|^2\le C4^jY_j$。Sobolev の不等式で $\norm{\zeta_j w_j}_{2n/(n-2)}^2\le C4^jY_j$ を得る。一方 $B_{r_{j+1}}$ 上で $w_{j+1}>0$ なら $w_j\ge k_{j+1}-k_j=K2^{-j-1}$ なので、その集合の測度は $4^{j+1}Y_j/K^2$ 以下である。従って Hölder の不等式から
$$
Y_{j+1}\le C4^jY_j \left(4^{j+1}Y_j/K^2\right)^{2/n}.
$$
$Z_j=Y_j/K^2$ とすると $Z_{j+1}\le C_0 b^j Z_j^{1+\delta}$、$\delta=2/n$、$b=4^{1+2/n}>1$ である。$Z_0\le C_0^{-1/\delta}b^{-1/\delta^2}$ なら帰納法で $Z_j\le C_0^{-1/\delta}b^{-j/\delta-1/\delta^2}$ となる。実際右の上界を漸化式に代入すると指数は $-(j+1)/\delta-1/\delta^2$ になる。従って $K=C\norm{u_+}_{L^2(B_1)}$ と十分大きく選べば $u\le K$ が $B_{1/2}$ で a.e. に成立する。$-u$ にも同じ評価を使えば局所有界性を得る。

この反復では初期の $u$ が有界であると仮定していない。他方、上下に有界になっただけでは振動の縮小は出ない。Harnack を非負な二つの差 $u-m$ と $M-u$ へ使うか、レベル集合の測度の減少をさらに評価して、半径を縮めたときの上下界の差を比べる。二次元と一次元では、勾配から有限な $L^p$ ノルムを抑える Sobolev の不等式を使う。$p>2$ を一つ固定すれば、上の $2/n$ を $1-2/p$ に置き換えて同じ漸化関係を得る。

## 演習問題

::::exercise{level=basic}
$u(x)=x^2$ に対する $D_hu$ とその極限を求めよ。

:::solution
$(x+h)^2-x^2=2xh+h^2$ なので $D_hu=2x+h$、$h\to0$ で $2x$ に収束する。任意の有界区間で $L^2$ 収束でもある。
:::
::::

::::exercise{level=basic}
$-\Delta u=0$ の単位球の弱解で $u_+\in H^1_0$ とするとき $u$ の符号を求めよ。

:::solution
$u_+$ で試験すると $\int|Du_+|^2=0$ である。Poincaré の不等式から $u_+=0$ なので $u\le0$ である。
:::
::::

::::exercise{level=basic}
$u=|x|^2-1$ の単位球境界での外向き法線微分を求めよ。

:::solution
外側を向く半径微分は $2r$ であり、$r=1$ なので値は $2$ である。
:::
::::

::::exercise{level=basic}
区間 $(0,\pi)$ の $u=\sin 2x$ について $\norm{u''}_2^2$ を求めよ。

:::solution
$u''=-4\sin2x$ なので二乗積分は $16(\pi/2)=8\pi$ である。
:::
::::

::::exercise{level=advanced}
楕円型の一様 Lipschitz 係数と滑らかな試験関数に対し差分商の積の法則を導け。

:::solution
$A(x+h)Du(x+h)-A(x)Du(x)$ に $-A(x+h)Du(x)+A(x+h)Du(x)$ を加える。$h$ で割ると $D_h(ADu)=A(x+h)D_hDu+(D_hA)Du(x)$ となる。$\norm{D_hA}_\infty$ は Lipschitz 定数以下である。
:::
::::

::::exercise{level=advanced}
$u=r^\gamma\sin(\gamma\theta)$ の勾配と Hessian の局所二乗可積分性の閾値を、$0<\gamma<1$ で求めよ。

:::solution
大きさはそれぞれ $C_\gamma r^{\gamma-1}$ と $C'_\gamma r^{\gamma-2}$ である。半径積分の指数は $2\gamma-1$ と $2\gamma-3$ である。積分可能性には指数 $>-1$ が必要なので、勾配は $\gamma>0$ で入るが Hessian には $\gamma>1$ が必要である。
:::
::::

::::exercise{level=advanced}
Harnack 定数 $C=3$ で半径を半分にする振動評価を 4 回使うと、元の振動の何倍以下か。

:::solution
一回の係数は $(C-1)/(C+1)=1/2$ であるから、4 回で $(1/2)^4=1/16$ 以下である。このとき半径も元の $1/16$ である。
:::
::::

::::exercise{level=advanced}
$A=\diag(9,4)$ の定数係数作用素を局所座標で Laplace 作用素に変えよ。

:::solution
$y_1=x_1/3,y_2=x_2/2$ とすれば $9\partial_{x_1x_1}+4\partial_{x_2x_2}=\partial_{y_1y_1}+\partial_{y_2y_2}$ である。$y$ の単位球に対応する $x$ の図形は半軸 $3,2$ の楕円である。
:::
::::

## 確認テスト

::::quiz
Lipschitz 係数の差分商を用いる直接の理由はどれか。

- [x] その係数の差分商を刻み幅によらず有界にできる
- [ ] 初めから解の二階古典微分が存在する
- [ ] 試験関数の境界値が任意になる

:::explanation
差分商で方程式を移すと $(D_hA)Du$ が出る。Lipschitz 性は $D_hA$ を一様に抑え、エネルギー不等式で吸収できる形にする。
:::
::::

::::quiz{answer="2" tol="1e-6"}
$u=|x|^2-1$ の単位球境界での外向き法線微分を求めよ。

:::explanation
半径方向の微分は $2r$ であり、境界の $r=1$ に代入して $2$ である。
:::
::::

::::quiz
一般の有界可測係数で成り立つものはどれか。

- [x] 斉次発散形の弱解には局所 Hölder 連続な代表がある
- [ ] 全ての弱解は $C^2$ である
- [ ] 全ての弱解で Hessian は有界である

:::explanation
一様楕円なら De Giorgi–Nash–Moser の結論が適用できる。区分一定の係数で勾配の跳びを持つ解は Hölder 連続でも二階微分が測度になる。
:::
::::

::::quiz{answer="0.5" tol="1e-6"}
Harnack 定数 $C=3$ で半径を半分にしたときの振動の上界係数 $(C-1)/(C+1)$ を求めよ。

:::explanation
$(3-1)/(3+1)=2/4=1/2$ なので係数は $0.5$ である。
:::
::::

::::quiz{answer="0.5" tol="1e-6"}
三次元の環 $1<r<2$ の放射状調和解で $u(1)=0,u(2)=1$ とする。外側の境界での法線微分を求めよ。

:::explanation
放射状調和性から $r^2u^{\prime}$ は定数で、$u=A+B/r$ である。境界値から $A=2,B=-2$ なので $u^{\prime}=2/r^2$。外側の法線は半径が増える向きであり、$r=2$ での値は $1/2$ である。
:::
::::


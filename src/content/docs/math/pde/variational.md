---
title: 変分法
description: エネルギーの第一変分と第二変分を計算し、直接法、拘束、峠の補題、対称性から偏微分方程式を導く。
notation:
  - symbol: 'J[u]'
    meaning: 関数 u に割り当てる汎関数の値
  - symbol: 'F_p'
    meaning: 被積分関数の勾配変数に関する偏微分
  - symbol: '\delta^2J[u](v,v)'
    meaning: u における方向 v の第二変分
---

## この章の目標

- 第一変分から Euler–Lagrange 方程式と自然境界条件を導く。
- 第二変分と凸性で停留点や最小化元の性質を調べる。
- 直接法と峠の補題を使って弱解を構成する。
- 拘束から乗数を求め、対称性から保存する流束を導く。

## 第一変分と Euler–Lagrange 方程式

体積や長さを測る式を関数に対して最小化すると、未知量は有限個の数から関数に変わる。そこで一点の偏微分の代わりに、許される関数の方向に微小に動かしたときの変化を考える。

:::definition[第一変分]{#def-pde-first-variation}
$J[u]=\int_\Omega F(x,u,Du)\,\dd x$ に対し、許される方向 $v$ の微分 $\frac{\dd}{\dd s}J[u+sv]|_{s=0}$ を :term[第一変分]{en="first variation" yomi="だいいちへんぶん"} という。すべての方向について $0$ なら $u$ は停留点である。
:::

:::proposition[変分方程式と自然境界条件]{#prop-pde-euler-lagrange}
$F,u$ が十分滑らかで、$\Omega$ は滑らかな有界領域とする。内部の滑らかな変分に対する停留条件は

$$
F_z(x,u,Du)-\divg F_p(x,u,Du)=0
\label{eq-pde-euler-lagrange}
$$

である。境界値も自由であれば、停留条件は内部方程式と境界条件 $F_p\cdot\nu=0$ の両方が成り立つことと同値である。
:::

:::proof
積分の下で微分して $\delta J[u](v)=\int(F_z v+F_p\cdot Dv)$ を得る。部分積分で体積項は :eqref[eq-pde-euler-lagrange] と $v$ の積、境界項は $\int_{\partial\Omega}(F_p\cdot\nu)v$ である。内部試験関数を選べば変分法の基本補題により PDE が従う。実際、連続な残差が一点で正なら同符号の小球内に試験関数を置くと積分が正になる。境界値を固定すれば境界項は消える。自由な場合、滑らかな境界関数を局所の法線座標で延長して $v$ に使えるため、境界項が全ての $v$ で消えることと $F_p\cdot\nu=0$ が同値である。逆も部分積分で従う。
:::

:::example[膜のエネルギー]
$F=|p|^2/2-f(x)z$ では $F_z=-f$、$F_p=p$ なので $-\Delta u=f$ となる。境界を自由にすると自然境界条件は $\partial_\nu u=0$ である。この場合定数変分 $v=1$ も可能なので $\int f=0$ が必要になる。
:::

:::example[グラフの面積]
$F=\sqrt{1+|p|^2}$ の変分から $\divg(Du/\sqrt{1+|Du|^2})=0$ を得る。一次元では $u'/\sqrt{1+(u')^2}$ が一定であり、左辺は $u'$ の狭義単調関数なので $u'$ は一定である。$(0,2)$ で端を $1,5$ に固定すれば $u=1+2x$ で、グラフの長さは $2\sqrt5$ である。
:::

### 自由端での Robin 条件と境界汎関数

$J[u]=\int_\Omega(|Du|^2/2-fu)+\int_{\partial\Omega}(\beta u^2/2-gu)$ とする。境界値を固定しないとき、任意の滑らかな $v$ について
$$
J'[u]v=\int_\Omega(-\Delta u-f)v+\int_{\partial\Omega}(\partial_\nu u+\beta u-g)v.
$$
内部に台を持つ $v$ により $-\Delta u=f$ を得てから、境界付近で任意の値を取る $v$ により $\partial_\nu u+\beta u=g$ を得る。この順番は、内部と境界で独立に試験できるためである。境界の一部で値を固定するなら、その部分で $v=0$ だから自然条件は残りの部分にだけ出る。$\beta\ge0$ のとき差のエネルギーは $\int|Dv|^2/2+\int\beta v^2/2$ であり、定数方向まで正にできれば停留点は唯一の最小化元である。

:::example[棒の一端を固定しもう一端をばねで支える]
$(0,1)$ で $u(0)=0$ を固定し、$J=\int_0^1((u')^2/2-2u)+u(1)^2$ とする。方程式は $-u''=2$、ばね条件は $u'(1)+2u(1)=0$ である。$u=-x^2+Ax$ を代入すると $A-2+2(A-1)=0$、従って $A=4/3$ である。端点の値 $u(1)=1/3$ に対し流束 $u'(1)=-2/3$ がばね力と打ち消し合う。$u(1)$ を固定して試験するとこの関係は導けない。
:::

:::example[多項式の変分方向に自然条件を入れる]
区間 $(0,1)$ で $u(0)=0$ とし、一次関数 $u(x)=ax$ を用いる。汎関数を $\int_0^1 (u_x^2/2-3u)+u(1)^2/2$ とする。

まず候補の選択を一変数 $a$ の問題にする。停留条件は $J_a=0$ であり、境界や無限遠も許される場合は停留点と分けて扱う。次の式では一階微分を取って零になる点を求め、二階微分を残しておく。

$$
\begin{aligned}
J(a)&=a^{2} - \frac{3 a}{2}\\
J_a&=\frac{d}{d a} \left(a^{2} - \frac{3 a}{2}\right)=\frac{4 a - 3}{2}\\
J_{aa}&=\frac{d^{2}}{d a^{2}} \left(a^{2} - \frac{3 a}{2}\right)=2\\
a_*&=\frac{3}{4}\\
J_a(a_*)&=0,\qquad a=\frac{3}{4}
\end{aligned}
$$

停留値と、その点から $h$ だけ動かした差を計算する。一次項が相殺しているため、展開の先頭から局所の振る舞いが分かる。一方、停留点が全体で最小か最大かを結論するには、許容範囲全体でこの差の符号を調べる必要がある。

$$
\begin{aligned}
J(a_*)&=- \frac{9}{16},\qquad a=\frac{3}{4}\\
J_{aa}(a_*)&=2\\
J(a_*+h)-J(a_*)&=- \frac{3 h}{2} + \left(h + \frac{3}{4}\right)^{2} - \frac{9}{16}=h^{2}\\
\text{二次係数}&=1
\end{aligned}
$$

指定した値で目的関数を比較する。$a$ 自体の値と、そこでの目的関数の値は別の量である。

$$
\begin{aligned}
a_*\text{ の具体値}&=\frac{3}{4}\\
J_*\text{ の具体値}&=- \frac{9}{16}\\
J&=0,\qquad a=0\\
J&=- \frac{1}{2},\qquad a=1\\
\text{二つの値の差}&=- \frac{1}{2}\\
\text{停留点まわりの差}&=h^{2}
\end{aligned}
$$

この $a_*$ はアフィン関数の族に限定した最小化元であり、一般の $H^1$ 変分で停留するとは限らない。例えば $u=ax$ では $-u_{xx}=0$ なので内部方程式 $-u_{xx}=3$ を満たさない。有限次元族での停留と全ての試験方向での停留を区別するための例である。
:::

### 積分条件と弱い Euler–Lagrange 方程式

$F\in C^1(\R^n)$、$|DF(p)|\le C(1+|p|)$、$f\in L^2$ とし、有界領域で $J(u)=\int (F(Du)-fu)$ を $H^1$ のアフィン空間で考える。$v\in H_0^1$ と $|s|\le1$ に対し $|DF(Du+sDv)\cdot Dv|\le C(1+|Du|+|Dv|)|Dv|\in L^1$ である。よって優収束で微分を積分に入れられるので、最小化元は $\int F_p(Du)\cdot Dv=\int fv$ を満たす。実際、$s$ の正負で $J(u+sv)-J(u)\ge0$ なので、差商の極限は零である。逆に $F$ が凸なら
$$
J(u+v)-J(u)\ge\int F_p(Du)\cdot Dv-\int fv=0
$$
であり、弱式だけで最小性を確認できる。さらに $D^2F\ge\lambda I$ なら、線分上の Taylor 公式から右辺に $\lambda\norm{Dv}_2^2/2$ が加わり、固定トレースのもとで最小化元は一意である。古典解として $\divg F_p$ を計算できることは、この弱い最小性の確認には必要ない。

## 第二変分

:::definition[第二変分]{#def-pde-second-variation}
$J[u+sv]$ の $s=0$ における二階導関数を :term[第二変分]{en="second variation" yomi="だいにへんぶん"} といい、滑らかな場合は

$$
\delta^2J[u](v,v)=\int_\Omega\big(F_{zz}v^2+2vF_{zp}\cdot Dv+(Dv)^\top F_{pp}Dv\big)\,\dd x
\label{eq-pde-second-variation}
$$

である。
:::

最小化元では二階導関数が非負でなければならない。一方、各方向で正という条件だけでは一様な局所最小は従わない。Hilbert 空間で $J$ が $C^2$、$J'(u)=0$、$J''(u)[v,v]\ge c\norm v^2$ なら局所最小になる。実際、$J''$ の作用素ノルム連続性により近傍では下界を $c/2$ にでき、$J[u+v]-J[u]=\int_0^1(1-s)J''(u+sv)[v,v]\dd s\ge c\norm v^2/4$ になる。

:::proposition[勾配変数についての必要条件]{#prop-pde-legendre-necessary}
滑らかな最小化元の内部では $F_{pp}(x,u,Du)$ は半正定値である。
:::

:::proof
点 $x_0$ と方向 $\xi\ne0$ に対し、近傍に台を持つ $\zeta$ と $v_N=N^{-1}\zeta(x)\sin(N\xi\cdot x)$ を使う。$v_N\to0$ 一様で、$Dv_N=\zeta\xi\cos(N\xi\cdot x)+O(N^{-1})$ である。:eqref[eq-pde-second-variation] の非負性で極限をとると $\frac12\int\zeta^2 \xi^\top F_{pp}\xi\ge0$ となる。ここで $\cos^2$ を定数項 $1/2$ と振動項に分ける。振動項の積分は Riemann–Lebesgue の補題で消える。もし $\xi^\top F_{pp}(x_0)\xi<0$ なら連続性により同符号の小球に $\zeta$ を取れ、矛盾する。
:::

未知関数が $\R^m$ に値を取る場合、同じ試験を $v_N=N^{-1}a\,\zeta\sin(N\xi\cdot x)$、$a\in\R^m$ に替える。主項の勾配は行列 $a\otimes\xi$ にスカラーを掛けたものである。従って最小化元に必要な条件は全ての $a,\xi$ に対し $F_{pp}[a\otimes\xi,a\otimes\xi]\ge0$ であり、これを :term[Legendre–Hadamard 条件]{en="Legendre–Hadamard condition" yomi="るじゃんどるあだまーるじょうけん"} という。スカラー未知関数では全ての勾配方向がこの形になる。一方、行列では階数二以上の方向もある。例えば二次元の二成分変分で $F(p)=\det p$ とすると、$D^2F(p_0)[q,q]=2\det q$ は階数一の方向で零だが、$q=\diag(1,-1)$ では負である。実際、$\det Du=\partial_1(u_1\partial_2u_2)-\partial_2(u_1\partial_1u_2)$ であり、内部に台を持つ変分でその積分は変わらない。従ってこの例では第二変分に点ごとの負の方向があっても、積分汎関数では全ての内部変分に対して第二変分が零である。

:::example[低いモードから不安定になる]
$(0,\pi)$ の Dirichlet 条件で $J[u]=\frac12\int((u')^2-\mu u^2)$ とする。第二変分は各モードで $(j^2-\mu)\norm v_2^2$ となる。$\mu=2$ なら $v=\sin x$ で $-\pi/2$、$v=\sin2x$ で $\pi$ となる。従って $u=0$ は停留点だが鞍点である。勾配変数の Hessian は $1$ で正でも、低いモードでは零階項の影響で負になり得る。
:::

:::example[零方向がある第二変分]
$J(a,b)=a^2-b^4$ の原点では第一変分が零、第二変分は $2a^2$ で半正定値である。ところが $b$ 軸上では任意に近くで $J<0$ となる。二階の項が消える方向では高次項が極値を決める。これに対して $a^2+b^4$ は同じ第一・第二変分を持つが原点が最小である。
:::

### 一様な二次下界から局所最小へ

$J$ は Banach 空間で $C^2$ とし、$J'(u)=0$、$J''(u)[h,h]\ge c\norm h^2$（$c>0$）を仮定する。二階微分が作用素ノルムで連続なら、十分小さい $h$ に対し線分上の二次形式は $J''(u+sh)[h,h]\ge(c/2)\norm h^2$ である。一次元の Taylor 公式を線分に適用すると
$$
J[u+h]-J[u]=\int_0^1(1-s)J''(u+sh)[h,h]\dd s\ge(c/4)\norm h^2.
$$
従って強い意味での厳密局所最小である。$J''(u)[h,h]>0$ が各 $h\ne0$ で成立しても、無限次元では一様な $c$ は出ない場合がある。

:::example[各方向には正でも一様性がない]
$\ell^2$ で $J(x)=\sum_j x_j^2/j-\sum_j x_j^4$ とする。四乗和は $\sum |x_j|^4\le\norm x_2^4$ であり多重線形形式から滑らかである。$J''(0)[h,h]=2\sum h_j^2/j>0$ だが、$x^{(j)}=\sqrt{2/j}\,e_j$ に対し $J(x^{(j)})=2/j^2-4/j^2=-2/j^2<0$ である。しかもノルムは零に近づく。有限次元における単位球のコンパクト性に相当する一様性が欠けている。
:::

:::example[方向によって符号が変わる第二変分]
二変数で臨界点の第二変分が対称であるとし、鞍点か判定する。

第二変分の行列を $H$ とする。対角成分は座標方向での値を与える。固有方程式から直交する方向を取り出すため、特性多項式 $\det(\lambda I-H)$ を計算する。

$$
\begin{aligned}
H&=\left[\begin{matrix}\frac{71}{25} & - \frac{72}{25}\\- \frac{72}{25} & \frac{29}{25}\end{matrix}\right]\\
\det(\lambda I-H)&=\lambda^{2} - 4 \lambda - 5=\left(\lambda - 5\right) \left(\lambda + 1\right)\\
d_1&=-1,\quad v_1=\left[\begin{matrix}\frac{3}{5}\\\frac{4}{5}\end{matrix}\right]\\
d_2&=5,\quad v_2=\left[\begin{matrix}- \frac{4}{5}\\\frac{3}{5}\end{matrix}\right]
\end{aligned}
$$

$V$ の列に固有ベクトルを並べる。$V^\top V=I$ と $V^\top HV=D$ を使えば、元の座標 $X$ と固有座標 $Y=V^\top X$ の双方で変分の係数を記述できる。二次形式は固有座標にすると交差項を持たない。

$$
\begin{aligned}
V&=\left[\begin{matrix}\frac{3}{5} & - \frac{4}{5}\\\frac{4}{5} & \frac{3}{5}\end{matrix}\right]\\
V^\top HV&=\left[\begin{matrix}-1 & 0\\0 & 5\end{matrix}\right]\\
X^\top HX&=\frac{71 x_{1}^{2} - 144 x_{1} x_{2} + 29 x_{2}^{2}}{25}\\
\sum_j d_jY_j^2&=- Y_{1}^{2} + 5 Y_{2}^{2}
\end{aligned}
$$

一つの変分の係数で固有係数を計算し、二次のエネルギー差を求める。直交変換なので二乗ノルムは変わらない。二次形式を二乗ノルムで割ると、固有値を非負の重みで平均した値になる。

$$
\begin{aligned}
X&=\left[\begin{matrix}1\\1\end{matrix}\right]\\
Y&=\left[\begin{matrix}\frac{7}{5}\\- \frac{1}{5}\end{matrix}\right]\\
|X|^2&=2\\
X^\top HX&=- \frac{44}{25}\\
\frac{X^\top HX}{|X|^2}&=- \frac{22}{25}\\
\min_j d_j&=-1\\
\max_j d_j&=5
\end{aligned}
$$

一方向で計算した第二変分が正でも、$v_1$ では $-1$、$v_2$ では $5$ になる。$J^{\prime}(0)=0$ なら Taylor 展開により一方の方向で値が下がり他方で上がる。従って臨界点は鞍点である。
:::

:::example[第二変分が零だと第四次まで必要になる]
実変数のエネルギー $J(a)=3+a^4$ を考える。臨界点零で二次の判定が使えないとき、差の最初の非零項を求める。

まず候補の選択を一変数 $a$ の問題にする。停留条件は $J_a=0$ であり、境界や無限遠も許される場合は停留点と分けて扱う。次の式では一階微分を取って零になる点を求め、二階微分を残しておく。

$$
\begin{aligned}
J(a)&=a^{4} + 3\\
J_a&=\frac{d}{d a} \left(a^{4} + 3\right)=4 a^{3}\\
J_{aa}&=\frac{d^{2}}{d a^{2}} \left(a^{4} + 3\right)=12 a^{2}\\
a_*&=0\\
J_a(a_*)&=0,\qquad a=0
\end{aligned}
$$

停留値と、その点から $h$ だけ動かした差を計算する。一次項が相殺しているため、展開の先頭から局所の振る舞いが分かる。一方、停留点が全体で最小か最大かを結論するには、許容範囲全体でこの差の符号を調べる必要がある。

$$
\begin{aligned}
J(a_*)&=3,\qquad a=0\\
J_{aa}(a_*)&=0\\
J(a_*+h)-J(a_*)&=h^{4}\\
\text{二次係数}&=0
\end{aligned}
$$

指定した値で目的関数を比較する。$a$ 自体の値と、そこでの目的関数の値は別の量である。

$$
\begin{aligned}
a_*\text{ の具体値}&=0\\
J_*\text{ の具体値}&=3\\
J&=4,\qquad a=-1\\
J&=\frac{49}{16},\qquad a=\frac{1}{2}\\
\text{二つの値の差}&=- \frac{15}{16}\\
\text{停留点まわりの差}&=h^{4}
\end{aligned}
$$

$J^{\prime\prime}(0)=0$ でも $J(h)-J(0)=h^4>0$ なので厳密な最小点である。反対に $3-a^4$ は同じ一階・二階微分を持つが最大点になる。第二変分の退化を知っただけでは、最小か最大かは決まらない。
:::

:::example[境界の重みが二次の安定性を変える]
$(0,1)$ で左端を零に固定し、$J(u)=\int_0^1 (u')^2/2-\beta u(1)^2/2$ とする。任意の変分で $v(1)=\int_0^1 v'$ だから $v(1)^2\le\int_0^1|v'|^2$ である。従って $0\le\beta<1$ なら第二変分は $(1-\beta)\norm{v'}_2^2$ 以上である。$\beta=1$ では $v(x)=x$ に対して等号で零となり、この方向は二次では安定性を与えない。$\beta>1$ なら同じ変分で $1-\beta<0$ となる。内部の密度は変化していないが、自然境界の二次項が符号を決める。
:::

## 凸性と直接法による最小化元の存在

:::theorem[直接法の基本形]{#thm-pde-direct-method}
Banach 空間 $X$ は反射的、許容集合 $\mathcal A$ は空でない弱閉集合とする。$J:\mathcal A\to\R\cup\{+\infty\}$ は下から有界で有限値を持つ点があり、弱逐次下半連続、かつ有界な劣位集合を持つとする。このとき最小化元が存在する。$\mathcal A$ が凸で $J$ が狭義凸なら、最小化元は一意である。
:::

:::proof
$m=\inf_\mathcal A J$ とし、$J(u_j)\downarrow m$ に近づく列を選ぶ。劣位集合の有界性により反射空間で弱収束部分列 $u_{j_k}\rightharpoonup u$ があり、弱閉性で $u\in\mathcal A$ である。下半連続性から $J(u)\le\liminf J(u_{j_k})=m$ である。異なる最小化元が二つあれば中点の値が狭義に小さくなるので一意性が従う。
:::

これを :term[直接法]{en="direct method" yomi="ちょくせつほう"} という。典型例では $X=W^{1,p}$、$1<p<\infty$ とする。$F(p)=|p|^p/p$ と線形の右辺を考えると、弱収束で勾配のノルムは下半連続である。これはノルムが双対単位球上の線形汎関数の上限で表されるためである。斉次 Dirichlet では Poincaré の不等式で $J[u]\ge \norm{Du}_p^p/p-C\norm{Du}_p$ となり劣位集合は有界である。境界を固定した集合は閉部分空間の平行移動なので弱閉である。$|p|^p$ の狭義凸性と零トレースから最小化元が一意になる。第一変分は $|Du|^{p-2}Du\in L^{p/(p-1)}$ と Hölder の不等式により意味を持ち、弱形式は $\int |Du|^{p-2}Du\cdot Dv=\ell(v)$ となる。

:::example[凸性で停留点が最小になる]
$F(x,z,p)$ が $(z,p)$ について微分可能かつ凸なら、接平面の不等式で $J[w]\ge J[u]+\delta J[u](w-u)$ となる。従って凸な許容集合の内部で全方向の第一変分が $0$ なら最小化元である。例えば固定端の $J=\int (u')^4/4$ ではアフィンな $u$ の第一変分が $\int(\mathrm{const})v'=0$ なので最小を取る。区間 $(0,1)$ で端点差 $2$ なら最小値は $4$ である。
:::

:::example[非凸性により弱極限の値が上がる]
$(0,1)$ の各小区間で導関数が半分ずつ $1,-1$ となる三角波 $u_n$ を、幅 $1/n$、両端 $0$ で作る。$u_n\to0$ は一様で、導関数は各区間で平均 $0$ なので連続な試験関数との積分が $0$ に収束し、密度により弱収束する。$J[u]=\int((u')^2-1)^2$ では $J[u_n]=0$ だが $J[0]=1$ である。従ってこの汎関数は弱下半連続でない。
:::

### 凸積分汎関数の弱下半連続性

:::proposition[積分汎関数の下半連続性]{#prop-pde-convex-integral-lsc}
$1<p<\infty$、$F(x,z,q)$ は $x$ に可測で a.e. の $x$ について $(z,q)$ に関して凸かつ連続とする。$F(x,z,q)\ge-a(x)$ となる $a\in L^1(\Omega)$ があれば $J[u]=\int_\Omega F(x,u,Du)$ は $W^{1,p}$ の弱収束に関して逐次下半連続である。
:::

:::proof
$u_j\rightharpoonup u$ とする。下極限が $+\infty$ なら自明なので、有限な下極限 $\ell$ へ $J[u_j]$ が収束する部分列を取る。弱収束列の各尾の凸包のノルム閉包に $u$ が入る。入らなければ Hahn–Banach の分離定理でその閉凸集合と $u$ を連続線形汎関数で厳密に分離でき、弱収束に反するからである。従って各尾からの凸結合 $v_k$ で $v_k\to u$ を $W^{1,p}$ の強収束にできる。強収束部分列をとると $v_k,Dv_k$ が a.e. で $u,Du$ に収束する。凸性で $J[v_k]\le\sup_{j\ge k}J[u_j]$ であり、Fatou の補題を $F+a$ に用いると $J[u]\le\liminf J[v_k]\le\ell$ である。
:::

この証明では強収束させるために凸結合を取った。例えば勾配変数だけが凸で $u$ への依存は非凸でも、$u_j$ が強い $L^p$ ノルムで収束し、その差によって積分値が小さく変わる条件があれば同じ結論が出る。$F(x,z,q)=|q|^p+W(z)$ で $|W'(z)|\le C(1+|z|^{p-1})$ なら、$W(u_j)-W(u)$ の積分を平均値の定理と Hölder で $C(1+\norm{u_j}_p^{p-1}+\norm u_p^{p-1})\norm{u_j-u}_p$ 以下にできる。勾配項は弱下半連続なので、この型では $W$ が凸でなくても下半連続性を証明できる。

### 勾配についてのみ凸な積分汎関数

$\Omega$ は有界開集合、$1<p<\infty$ とする。$F(x,z,q)$ は $x$ に可測、$(z,q)$ について連続かつ $q$ について凸であり、連続な $D_qF$ を持つとする。次の成長条件を仮定する。

$$
|F(x,z,q)|\le C(1+|z|^p+|q|^p),\qquad
|D_qF(x,z,q)|\le C(1+|z|^{p-1}+|q|^{p-1}) .
$$

$u_j\to u$ が $L^p$ で強収束し、$Du_j\rightharpoonup Du$ が $L^p$ で弱収束するとき、$J(u)\le\liminf J(u_j)$ である。実際、勾配変数の接平面で

$$
F(x,u_j,Du_j)\ge F(x,u_j,Du)+D_qF(x,u_j,Du)\cdot(Du_j-Du)
$$

である。下極限を実現する部分列から、さらに $\sum_j\norm{u_j-u}_p<\infty$ となるよう選ぶ。このとき $|u_j|\le |u|+\sum_i|u_i-u|$ は一つの $L^p$ 関数で抑えられ、点ごとの収束も成り立つ。成長条件と優収束で第一項は $L^1$ で極限に移り、第二項の係数は $L^{p/(p-1)}$ で強収束する。勾配の差は弱収束なので第二項の積分は零に近づく。これで下半連続性が従う。

例えば $F(x,z,q)=a(z)|q|^p/p+W(z)$ で $a$ は正の有界連続関数、$W$ は連続で $p$ 次以下の成長とする。勾配についての凸性は $a>0$ から従うので、$z$ について非凸でも上の議論が使える。固定トレースの場合、$W_0^{1,p}$ の有界列を零延長して Rellich の定理を使うと、必要な $L^p$ 強収束が得られる。ただし最小化元の存在には、これに加えて劣位集合の有界性が必要である。

### 障害物制約と変分不等式

床を表す関数 $\psi$ に対し $\mathcal A=\{u\in H^1_0(\Omega):u\ge\psi\ {\rm a.e.}\}$ とし、$\Omega$ は有界開集合、$f\in L^2(\Omega)$、$\psi$ は可測で $\mathcal A$ は空でないとする。許容集合は凸で強閉である。実際、$H^1_0$ での強収束から a.e. 収束部分列があり不等式を保つ。凸性により強閉性から弱閉性が従う。Dirichlet エネルギー $J[u]=\int |Du|^2/2-\int fu$ は強圧的で弱下半連続なので、:ref[thm-pde-direct-method] により一意な最小化元を持つ。任意の $v\in\mathcal A$ に対し線分 $u+s(v-u)$ の $s=0$ の右微分は非負でなければならず、

$$
\int Du\cdot D(v-u)\ge\int f(v-u)
\label{eq-pde-obstacle-inequality}
$$

を得る。逆にこの不等式があれば二次のエネルギー差を展開して最小性が従う。非負の変分は全て許されるため $-\Delta u\ge f$ が分布で成り立つ。連続な $u,\psi$ を持つ場合、$u>\psi$ の開集合に台を持つ変分には両符号が許されるので、その集合では $-\Delta u=f$ となる。

:::example[接触点に残る力]
$(0,1)$ で $\psi(x)=1-2|x-1/2|$、$f=0$ とする。$u=\psi$ では弱微分が左 $2$、右 $-2$ であり、$-u''=4\delta_{1/2}$ である。任意の許容関数 $v$ は一次元の Sobolev 連続代表を持ち、$v(1/2)\ge1$ である。従って $\int u'(v-u)'=4(v(1/2)-1)\ge0$ であり、:eqref[eq-pde-obstacle-inequality] により $u$ が最小化元である。ここで点に集中した力は床の制約が許す反作用である。
:::

### 弱収束から強収束を得る条件

$J$ が二点の中点で
$$
J\big((u+v)/2\big)\le\tfrac12J(u)+\tfrac12J(v)-\kappa\norm{u-v}^p
$$
を満たし、許容集合が凸で、$u$ がその最小化元とする。最小化列 $u_j$ との中点を入れると $\kappa\norm{u-u_j}^p\le (J(u_j)-J(u))/2$ となる。従って列は強収束する。Dirichlet エネルギーは勾配の平行四辺形恒等式により $p=2$、$\kappa=1/8$ の形を満たす。固定トレースのアフィン集合では $J(u_j)-J(u)=\norm{D(u_j-u)}_2^2/2$ という等式が第一変分から成立する。障害物のような凸制約では、これに非負の一次項 $J^{\prime}(u)(u_j-u)$ が加わる。安定性に必要なのは単なる凸性だけでなく、二点が離れたときのエネルギーの落ち方である。

:::example[床に接触する区間を求める]
$(-1,1)$ で $J(u)=\int((u')^2/2+2u)$、$u(\pm1)=1/4$、制約 $u\ge0$ とする。$u(x)=(|x|-1/2)_+^2$ は境界値を満たし、接触の端点で $u'$ が連続なので $u''=2$ が非接触部分、$0$ が接触部分に成立する。変分不等式の左辺は $\int(-u''+2)(v-u)=2\int_{-1/2}^{1/2}v\ge0$ である。従って凸性から最小化元である。非接触部分だけで通常の Euler–Lagrange 方程式が成立し、接触部分では制約の反力 $2$ が残る。
:::

:::example[勾配に重みを付けたエネルギーの最小化]
$(0,1)$ で $J(v)=\int_0^1((2+x)v_x^2/2-fv)$ を考える。零境界を課し、$u=x(1-x)$ となるよう $f$ を決める。

まず流束 $q=Au_x$ を未知量にする。$q_x=-f$ を積分すると、任意の定数 $q(0)$ について $q(x)=q(0)-\int_0^x f$ である。さらに $q/A$ を積分し、指定された両端の値または平均に合わせて積分定数を決める。以下では計算後の $u$ について関係を逆に確かめる。

$$
\begin{aligned}
u(x)&=x \left(1 - x\right)\\
u_x&=\frac{d}{d x} x \left(1 - x\right)=1 - 2 x\\
q&=\left(1 - 2 x\right) \left(x + 2\right)=- \left(x + 2\right) \left(2 x - 1\right)\\
-f&=\frac{d}{d x} \left(- \left(x + 2\right) \left(2 x - 1\right)\right)=- 4 x - 3
\end{aligned}
$$

次に値の条件と流出量を両端で区別して評価する。区間の外法線は左端で $-1$、右端で $1$ である。このため源項の積分は左側の流束から右側の流束を引いたものになる。流束を係数で割った積分は両端の値の差に等しい。

$$
\begin{aligned}
u(0)&=0,\qquad x=0\\
u(L)&=0,\qquad x=1\\
q(0)&=2,\qquad x=0\\
q(L)&=-3,\qquad x=1\\
\int_0^L f\,\mathrm{d}x&=\int\limits_{0}^{1} \left(4 x + 3\right)\, dx=5\\
\int_0^L (q/A)\,\mathrm{d}x&=\int\limits_{0}^{1} \left(1 - 2 x\right)\, dx=0
\end{aligned}
$$

弱形式も具体的に比べる。試験関数 $\phi=x(L-x)$ は両端で零であり、積分の間に境界項は出ない。以下の値が一致するのは、任意の零トレース試験関数に対しても部分積分できることの一例である。

$$
\begin{aligned}
\phi&=x \left(1 - x\right)\\
\phi^{\prime}&=\frac{d}{d x} x \left(1 - x\right)=1 - 2 x\\
\int_0^L q\phi^{\prime}\,\mathrm{d}x&=\int\limits_{0}^{1} \left(- \left(1 - 2 x\right) \left(x + 2\right) \left(2 x - 1\right)\right)\, dx=\frac{5}{6}\\
\int_0^L f\phi\,\mathrm{d}x&=\int\limits_{0}^{1} x \left(1 - x\right) \left(4 x + 3\right)\, dx=\frac{5}{6}
\end{aligned}
$$

最後に内部での値、勾配、流束を同じ点で測る。区間全体での勾配の寄与には $A$ の重みが付くので、点での値だけでは全体のエネルギーは分からない。

$$
\begin{aligned}
u&=\frac{1}{4},\qquad x=\frac{1}{2}\\
u_x&=0,\qquad x=\frac{1}{2}\\
q&=0,\qquad x=\frac{1}{2}\\
\frac12\int_0^L A(u_x)^2\,\mathrm{d}x&=\int\limits_{0}^{1} \frac{\left(1 - 2 x\right)^{2} \left(x + 2\right)}{2}\, dx=\frac{5}{12}
\end{aligned}
$$

$J(u+h)-J(u)=\int_0^1(2+x)h_x^2/2$ は弱形式を任意の $h\in H_0^1$ で使って一次項を消すと従う。$2\le A\le3$ なので誤差エネルギーは $\norm{h_x}_2^2$ 以上、$3\norm{h_x}_2^2/2$ 以下である。この両側の評価は強圧性と連続性を同時に反映する。
:::

:::example[直接法と有限次元近似のエネルギー差]
$(0,1)$ で $A=1$、$c=1$、厳密解 $u=x^3(1-x)$ とする。二つの多項式の族に制限した変分解を求める。

剛性行列を作る前に、方程式と試験空間を決める。ここでは全ての基底が両端で消える。厳密解を含まない試験空間に制限すると、近似は割り当てた方向の弱形式を満たすように決まる。厳密解は誤差を測るときに使う。

$$
\begin{aligned}
A&=1\\
u&=x^{3} \left(1 - x\right)\\
f&=- x \left(x^{3} - x^{2} - 12 x + 6\right)\\
\psi_1&=x \left(1 - x\right)\\
\psi_2&=x^{2} \left(1 - x\right)^{2}
\end{aligned}
$$

剛性の成分は $K_{ij}=\int A\psi^{\prime}_j\psi^{\prime}_i+\int c\psi_j\psi_i$、右辺は $F_i=\int f\psi_i$ である。基底を正規直交には取っていないので、内積の質量行列 $M_{ij}=\int\psi_j\psi_i$ も別に求める。

$$
\begin{aligned}
K&=\left[\begin{matrix}\frac{11}{30} & \frac{31}{420}\\\frac{31}{420} & \frac{13}{630}\end{matrix}\right]\\
M&=\left[\begin{matrix}\frac{1}{30} & \frac{1}{140}\\\frac{1}{140} & \frac{1}{630}\end{matrix}\right]\\
F&=\left[\begin{matrix}\frac{23}{210}\\\frac{41}{2520}\end{matrix}\right]\\
\det(K-\lambda M)&=\frac{\lambda^{2}}{529200} - \frac{19 \lambda}{88200} + \frac{1121}{529200}=\frac{\lambda^{2} - 114 \lambda + 1121}{529200}
\end{aligned}
$$

連立方程式 $Kc=F$ を解き、元の関数に戻す。残差 $r=f-Lu_N$ を計算すると一般に零ではないが、選んだ試験関数に対する積分は全て零になる。PDE を点で満たすことと有限次元の弱形式を満たすことを区別する。

$$
\begin{aligned}
c&=\left[\begin{matrix}\frac{1}{2}\\-1\end{matrix}\right]\\
u_N&=- \frac{x \left(x - 1\right) \left(2 x^{2} - 2 x + 1\right)}{2}\\
e&=u-u_N=- \frac{x \left(x - 1\right) \left(2 x - 1\right)}{2}\\
r&=- \frac{\left(x - 3\right) \left(x + 2\right) \left(2 x - 1\right)}{2}\\
(\int r\psi_i)_i&=\left[\begin{matrix}0\\0\end{matrix}\right]
\end{aligned}
$$

最後に誤差の二乗ノルムを二通りに計算する。近似空間へのエネルギー直交性により、誤差の勾配エネルギーは全体のエネルギーから近似の分を引いた値になる。ここでは二次形式そのものを用いており、エネルギー汎関数の二分の一は付けていない。

$$
\begin{aligned}
a(u,u)&=\frac{113}{1260}\\
a(u_N,u_N)&=\frac{97}{2520}\\
a(e,e)&=\frac{43}{840}\\
\norm e_2^2&=\int\limits_{0}^{1} \frac{x^{2} \left(x - 1\right)^{2} \left(2 x - 1\right)^{2}}{4}\, dx=\frac{1}{840}
\end{aligned}
$$

二つの基底関数 $\psi_1,\psi_2$ は $x=1/2$ に関して対称である。従って厳密解の反対称部分は試験空間で表せず、誤差に残る。変分法は許容集合内で最小値を達成するが、空間を制限した場合には省いた方向まで誤差が消えるわけではない。
:::

### 強制項が弱極限に移る積分の組合せ

$J(u)=\int F(Du)+\int W(u)-\int fu$ を $H_0^1$ で最小にする。$F$ は凸で二次の強圧下界を持ち、$W$ は下に有界な連続関数で $|W(s)|\le C(1+s^2)$、$f\in L^2$ とする。最小化列は強圧性で $H_0^1$ に有界である。弱収束部分列を選び、Rellich から $L^2$ 強収束と a.e. 収束を得る。従って $\int fu_j\to\int fu$。一方 $W$ の項では $W$ の下界を引いて Fatou の補題を使うと $\int W(u)\le\liminf\int W(u_j)$ である。$F$ の項は勾配の弱収束に関する凸積分の下半連続性から同じ評価を得る。従って $W$ が非凸でも、非凸性が値の変数にだけ入り、値の強収束が得られるなら直接法が使える。この結論を勾配変数に非凸な $F$ へそのまま移すことはできない。

## 最小化元の正則性

最小化元の存在はある Sobolev 空間の元を得ることを意味する。古典的微分が使えるかは、第一変分の方程式に帰って調べる。例えば $J[u]=\int(|Du|^2/2-fu)$ の最小化元には :ref[thm-pde-elliptic-interior-h2] が適用できる。$f$ が局所で滑らかなら繰り返し微分して古典解になる。

:::proposition[一様凸な勾配エネルギーの内部二階正則性]{#prop-pde-variational-h2}
$F\in C^2(\R^n)$ が $\lambda I\le F_{pp}\le\Lambda I$ を満たし、$u\in H^1_{\rm loc}$ が $-\divg F_p(Du)=f$ の弱解、$f\in L^2_{\rm loc}$ なら $u\in H^2_{\rm loc}$ である。
:::

:::proof
差分商 $w=D_hu$ に対して

$$
D_h(F_p(Du))=A_h Dw,\qquad
A_h(x)=\int_0^1F_{pp}((1-s)Du(x)+sDu(x+he_j))\,\dd s
$$

であり、$A_h$ の一様楕円性と上界は $h$ に依存しない。:ref[thm-pde-elliptic-interior-h2] と同じ試験 $-D_{-h}(\zeta^2w)$ を用いると左辺は $\int A_h Dw\cdot D(\zeta^2w)$ になる。ここでは係数 $A_h$ を微分しない。右辺は $\norm f_2\norm{D(\zeta^2w)}_2$ 以下であり、Young の不等式と $\norm w_2\le\norm{Du}_2$ で $\int\zeta^2|Dw|^2$ を一様に抑えられる。差分商の弱収束により二階弱微分を得る。
:::

:::example[退化により二階微分が発散する]
$J=\int_{-1}^1((u')^4/4-u)\,\dd x$ の零境界最小化元では $-((u')^3)'=1$ である。$u'(x)=\sqrt[3]{-x}$、$u(x)=3(1-|x|^{4/3})/4$ と置くと端点値は $0$、流束は $-x$ で、弱形式が従う。凸性からこの停留点が最小化元である。$u''$ は原点外で $-|x|^{-2/3}/3$ なので $L^2$ には入らない。この例では $F_{pp}=3p^2$ に正の一様下界がない。
:::

一般の $p$-Laplace 方程式の正則性やベクトル値最小化元の部分正則性には、勾配が小さい領域と大きい領域に応じた評価が必要である。完全な正則性理論は :cite[evans-pde] と :cite[gilbarg-trudinger] を参照する。

:::example[滑らかな強圧エネルギーの非線形方程式]
一次元で $F(p)=p^2/2+1-\cos p$ では $F_{pp}=1+\cos p$ が零になることがある。$F$ に $p^2/2$ を足せば $2+\cos p\in[1,3]$ となり一様凸である。後者で $f=0$ なら $\partial_x(2u'+\sin u')=0$ であり、$2u'+\sin u'$ は一定である。その写像の導関数が $1$ 以上なので $u'$ が一定となる。両端を $0,2$ とすると区間 $(0,1)$ の解は $u=2x$ である。
:::

### 勾配の方程式へ引き継ぐ楕円性

$J(u)=\int F(Du)-fu$ で $F\in C^2$、$\lambda I\le D^2F\le\Lambda I$ とする。内部 $H^2$ が得られた後、$w_j=\partial_j u$ は
$$
-\divg(D^2F(Du)Dw_j)=\partial_j f
$$
を満たす。実際 $DF$ は大域的 Lipschitz なので Sobolev の連鎖律で $D_j(DF(Du))=D^2F(Du)Dw_j$ と書ける。滑らかな $u_m\to u$ を $H^2$ で近似すれば、$D^2F$ の有界性と連続性によりこの連鎖律を確認できる。右辺が零の場合、係数 $D^2F(Du)$ が当初は有界可測でも一様楕円なので :ref[thm-pde-moser-harnack] と振動縮小を使い、各 $w_j$ は局所 Hölder 連続になる。ここでは未知量がスカラーで、成分ごとに同じスカラー方程式が現れることを使う。ベクトル値の $u$ にそのままスカラーの Harnack 評価を適用することはできない。

:::example[二次の近似で残る高次項]
$F(p)=p^2/2+\eps(1-\cos p)$、$0<\eps<1$ とすると $F''=1+\eps\cos p$ で一様凸である。一次元で $f=0$ なら $(DF(u'))'=0$ なので $u'+\eps\sin u'$ が定数である。この写像の微分は $1-\eps$ 以上だから $u'$ は定数であり、端点値だけで決まる。$u(0)=0,u(2)=3$ では $u'=3/2$ で、流束は $3/2+\eps\sin(3/2)$ である。同じ解の勾配でも、流束は二次エネルギーの場合と異なる。
:::

### 一次元での正則性と保存流束

$J(u)=\int F(x,u')-f(x)u$ とし、$F_{pp}\ge\lambda>0$、$F_p$ は滑らかとする。弱い第一変分は $(F_p(x,u'))'=-f$ を意味するので、$f$ が連続なら流束は $C^1$ の代表 $q$ を持つ。各 $x$ で $F_p(x,p)$ は $p$ に関して傾きが $\lambda$ 以上なので全実数へ増加全単射である。$F_p(x,p)=q(x)$ を陰関数定理で解いて $p$ が $C^1$ となり、$u$ は $C^2$ の代表を持つ。微分すると
$$
u''=\frac{-f(x)-F_{xp}(x,u')}{F_{pp}(x,u')}.
$$
従って二階微分の上界を得るために必要なのは、右辺の滑らかさだけでなく、$F_{pp}$ が零から離れることである。多次元では $F_p$ はベクトルであり発散だけが指定されるので、同じ一回の積分で各成分を決められない。

:::example[位置によって重みが変わる棒]
$F(x,p)=(1+x)p^2/2$、$f=0$、$(0,1)$ で $u(0)=0,u(1)=1$ とする。$F_p=(1+x)u'$ は定数なので $u'=C/(1+x)$、$u=C\log(1+x)$。端点から $C=1/\log2$ である。従って $u''=-1/((1+x)^2\log2)$ であり、力がない場合でも係数の変化に応じて勾配は減る。一方、流束は空間全体で一定である。弾性エネルギーは $C^2\log 2/2=1/(2\log2)$ となる。
:::

## 拘束条件と Lagrange 乗数

:::proposition[一つの等式拘束の乗数]{#prop-pde-lagrange-variational}
Banach 空間の $C^1$ 汎関数 $J,G$ に対し、$G(u)=c$ のもとで $u$ が局所極値で $G'(u)\ne0$ なら、ある $\lambda\in\R$ に対し $J'(u)=\lambda G'(u)$ である。
:::

:::proof
$G'(u)w=1$ となる $w$ を選ぶ。任意の $v\in\ker G'(u)$ に対して $G(u+sv+r w)=c$ を $r$ について陰関数定理で解くと $r(0)=0$、$r'(0)=0$ である。拘束内の曲線に沿った $J$ の微分が消えることから $J'(u)v=0$ である。任意の $v$ は $v-G'(u)v\,w$ と $G'(u)v\,w$ の和なので、$J'(u)v=J'(u)w\cdot G'(u)v$ である。
:::

:::example[固有値が乗数になる]
$J[u]=\frac12\int|Du|^2$ を $G[u]=\frac12\int u^2=1/2$ のもとで最小化する。$G'(u)u=1$ であり、第一変分は $\int Du\cdot Dv=\lambda\int uv$ になる。$v=u$ で $\lambda=\int|Du|^2$ を得る。コンパクト埋め込みにより $\norm u_2=1$ が最小化列の弱 $H^1_0$ 極限で保たれるので、直接法の存在証明と :eqref[eq-pde-elliptic-minmax] の $k=1$ が一致する。
:::

:::example[積分値を固定した区間問題]
$(0,1)$ の両端を零に固定し、$\int_0^1u=1$ の下で $\frac12\int (u')^2$ を最小化する。乗数方程式は $-u''=\lambda$ なので $u=\lambda x(1-x)/2$ である。拘束から $\lambda/12=1$、従って $u=6x(1-x)$ である。同じ拘束を満たす摂動では一次項が $\lambda\int v=0$ となり、エネルギー差は $\int(v')^2/2$ であるため、これが最小である。
:::

### 複数の独立な拘束

$G=(G_1,\ldots,G_m):X\to\R^m$ が $C^1$、$DG(u)$ が全射であるとする。$DG(u)e_j$ が標準基底になる $e_j\in X$ を選べる。$h\in\ker DG(u)$ に対し $G(u+th+\sum z_j(t)e_j)=G(u)$ を陰関数定理で解くと $z(0)=0$、$z'(0)=0$ である。従って拘束上の停留点では $J'(u)h=0$。任意の $v$ を $v-\sum DG_j(u)v\,e_j$ と残りに分けると
$$
J'(u)v=\sum_j J'(u)e_j\,DG_j(u)v .
$$
従って乗数が存在し、全射性から一意である。全射性がなければ、拘束の接空間が $\ker DG(u)$ に一致するとは限らない。例えば実変数で拘束 $G(x)=x^2=0$ の許容集合は $\{0\}$ なので $J(x)=x$ は拘束上の最小値を持つが、$J'(0)=1$ は $DG(0)=0$ の倍数にできない。

:::example[拘束の接方向に限った第二変分]
線形拘束 $v_1^\top X=0$ を課し、三変数の二次形式をこの平面に制限する。

第二変分の行列を $H$ とする。対角成分は座標方向での値を与える。固有方程式から直交する方向を取り出すため、特性多項式 $\det(\lambda I-H)$ を計算する。

$$
\begin{aligned}
H&=\left[\begin{matrix}\frac{2473}{845} & - \frac{3756}{845} & \frac{240}{169}\\- \frac{3756}{845} & \frac{282}{845} & - \frac{180}{169}\\\frac{240}{169} & - \frac{180}{169} & \frac{463}{169}\end{matrix}\right]\\
\det(\lambda I-H)&=\lambda^{3} - 6 \lambda^{2} - 13 \lambda + 42=\left(\lambda - 7\right) \left(\lambda - 2\right) \left(\lambda + 3\right)\\
d_1&=-3,\quad v_1=\left[\begin{matrix}\frac{3}{5}\\\frac{4}{5}\\0\end{matrix}\right]\\
d_2&=2,\quad v_2=\left[\begin{matrix}- \frac{4}{13}\\\frac{3}{13}\\\frac{12}{13}\end{matrix}\right]\\
d_3&=7,\quad v_3=\left[\begin{matrix}\frac{48}{65}\\- \frac{36}{65}\\\frac{5}{13}\end{matrix}\right]
\end{aligned}
$$

$V$ の列に固有ベクトルを並べる。$V^\top V=I$ と $V^\top HV=D$ を使えば、元の座標 $X$ と固有座標 $Y=V^\top X$ の双方で変分の係数を記述できる。二次形式は固有座標にすると交差項を持たない。

$$
\begin{aligned}
V&=\left[\begin{matrix}\frac{3}{5} & - \frac{4}{13} & \frac{48}{65}\\\frac{4}{5} & \frac{3}{13} & - \frac{36}{65}\\0 & \frac{12}{13} & \frac{5}{13}\end{matrix}\right]\\
V^\top HV&=\left[\begin{matrix}-3 & 0 & 0\\0 & 2 & 0\\0 & 0 & 7\end{matrix}\right]\\
X^\top HX&=\frac{2473 x_{1}^{2} - 7512 x_{1} x_{2} + 2400 x_{1} x_{3} + 282 x_{2}^{2} - 1800 x_{2} x_{3} + 2315 x_{3}^{2}}{845}\\
\sum_j d_jY_j^2&=- 3 Y_{1}^{2} + 2 Y_{2}^{2} + 7 Y_{3}^{2}
\end{aligned}
$$

一つの変分の係数で固有係数を計算し、二次のエネルギー差を求める。直交変換なので二乗ノルムは変わらない。二次形式を二乗ノルムで割ると、固有値を非負の重みで平均した値になる。

$$
\begin{aligned}
X&=\left[\begin{matrix}1\\1\\2\end{matrix}\right]\\
Y&=\left[\begin{matrix}\frac{7}{5}\\\frac{23}{13}\\\frac{62}{65}\end{matrix}\right]\\
|X|^2&=6\\
X^\top HX&=\frac{5703}{845}\\
\frac{X^\top HX}{|X|^2}&=\frac{1901}{1690}\\
\min_j d_j&=-3\\
\max_j d_j&=7
\end{aligned}
$$

試験ベクトルの第1固有係数は非零であり、そのままでは拘束を満たさない。射影して第1係数を零にすれば二次形式は正であり、単位ベクトルでの下限は $2$ になる。拘束のない空間での負の方向は許容集合から除かれる。
:::

### 拘束付き第二変分の式

$J'(u)=\sum_i\lambda_iG_i'(u)$ とし、$G(u)=0$ を満たす $C^2$ 曲線 $w(s)$ で $w(0)=u,w'(0)=h$ とする。二階微分で
$$
G_i'(u)w''(0)+G_i''(u)[h,h]=0
$$
である。これを $J(w(s))$ の二階微分に代入すると
$$
(J(w(s)))''|_{s=0}=J''(u)[h,h]-\sum_i\lambda_iG_i''(u)[h,h]
$$
となる。従って拘束の接方向 $DG(u)h=0$ で局所最小なら右辺が非負でなければならない。$DG(u)$ が全射なら、前節の陰関数の構成により各接方向をこのような曲線の速度として実現できる。曲がった拘束面では $J''$ だけを見ると必要な符号が変わり得る。

:::example[線形汎関数を球面で最小にする]
$J(x,y)=x$、$G=x^2+y^2-1=0$ とする。$u=(-1,0)$ で $J'=(1,0)$、$G'=(-2,0)$ なので $\lambda=-1/2$ である。接方向は $h=(0,a)$ に限られ、$J''=0$ だが $-\lambda G''[h,h]=a^2$ となる。実際 $w(s)=(-\cos(as),\sin(as))$ で上の二階微分は $a^2$ である。逆に $(1,0)$ では $\lambda=1/2$ となり、同じ接方向で第二変分は負になる。
:::

## 峠の補題と臨界点理論

:::definition[Palais–Smale 条件]{#def-pde-palais-smale}
$J(u_j)$ が有界で $\norm{J'(u_j)}\to0$ となる列が強収束部分列を持つとき、$J$ は :term[Palais–Smale 条件]{en="Palais–Smale condition" yomi="ぱれすめーるじょうけん"} を満たすという。
:::

:::theorem[Hilbert 空間における峠の補題]{#thm-pde-mountain-pass}
実 Hilbert 空間上の $J\in C^2$ が Palais–Smale 条件を満たすとする。$J(0)=0$、$\norm u=r$ 上で $J(u)\ge a>0$、ある $e$ に対し $\norm e>r$、$J(e)<0$ とする。$0$ から $e$ への連続な道全体 $\mathcal P$ に対して

$$
c=\inf_{\gamma\in\mathcal P}\max_{s\in[0,1]}J(\gamma(s))
\label{eq-pde-mountain-pass}
$$

は臨界値であり $c\ge a$ である。
:::

:::proof
任意の道は半径 $r$ の球面を横切るので $c\ge a$ であり、直線の道で $c$ は有限である。$c$ に臨界点がないと仮定する。もし任意の小さい帯 $|J-c|\le2\eps$ で勾配ノルムの下限が $0$ に近づくなら、Palais–Smale 列を選べてその極限が $c$ にある臨界点になり矛盾する。従ってある帯で $\norm{\nabla J}\ge d>0$ である。$\eps$ を小さくして端点の値を帯の外に置き、$|J-c|\le\eps$ で $1$、$|J-c|\ge2\eps$ で $0$ になる滑らかな切断 $\chi(J)$ を使う。ベクトル場 $-\chi(J)\nabla J/\norm{\nabla J}^2$ は帯の外で $0$ と延長でき、局所 Lipschitz でノルムは $1/d$ 以下である。有限の時間に軌道は Cauchy 上界を持ち、完全性と局所存在により流れを延長できる。流れに沿って $J'_{\rm time}=-\chi(J)$ なので値は増えず、中央の帯では速度 $1$ で減る。最大値が $c+\eps/2$ 以下の道を時間 $2\eps$ 流すと、端点は動かず最大値は $c-\eps/2$ 以下になる。これは $c$ の定義に反する。
:::

:::example[劣臨界の半線形方程式]
有界領域で $X=H^1_0$、$2<p<2n/(n-2)$（$n\le2$ では有限の $p>2$）とし、$J[u]=\frac12\norm{Du}_2^2-\frac1p\norm u_p^p$ とする。埋め込みにより半径 $r$ が小さいと $J[u]\ge r^2/2-Cr^p>0$ であり、任意の非零の $\phi$ を十分大きく倍すると $J[t\phi]<0$ である。$p>2$ のため $J$ は $C^2$ であり、二階変分の下位項の連続性は Hölder の不等式と $L^p$ 埋め込みから従う。Palais–Smale 列には

$$
J(u_j)-p^{-1}\langle J'(u_j),u_j\rangle=(1/2-1/p)\norm{Du_j}_2^2
$$

より有界性がある。Rellich により部分列は $X$ で弱収束、$L^p$ で強収束する。弱極限 $u$ に対し $\langle J'(u_j),u_j-u\rangle\to0$ であり、非線形項は $L^p$ 強収束と $\norm{u_j}_p^{p-1}$ の有界性で消える。従って $\int Du_j\cdot D(u_j-u)\to0$ となり、弱収束を引いて $\norm{D(u_j-u)}_2^2\to0$ を得る。峠の臨界点は $J>0$ なので非零で、$-\Delta u=|u|^{p-2}u$ の弱解である。
:::

:::example[臨界点は全体の最小値でなくても得られる]
$J(s)=s^2-s^4$ とすると $J(0)=0$、$r=1/2$ で $J=3/16>0$、$e=2$ で $J=-12$ である。この汎関数は下に非有界だが、$|J'(s)|=|2s-4s^3|\to\infty$ が $|s|\to\infty$ で成り立つため Palais–Smale 列は有界である。$0$ と $2$ を結ぶ道は $s=1/\sqrt2$ を通り、その値は $1/4$ である。直線上の最大もそこで達成されるため、これが峠の値である。
:::

### 非零解を Nehari 拘束で捉える

$p>2$ とし、$H_0^1(\Omega)\cap L^p(\Omega)$ に二つのノルムの和を入れる。$J(u)=\tfrac12\norm{Du}_2^2-\tfrac1p\norm u_p^p$ はこの空間上で $C^1$ である。$u\ne0$ の方向に $t>0$ を動かすと $J(tu)=t^2 A/2-t^pB/p$（$A=\norm{Du}_2^2,B=\norm u_p^p$）である。導関数が零になる点は $t=(A/B)^{1/(p-2)}$ で、その前後で導関数の符号が正から負に変わる。従って各半直線と $G(u)=A-B=0$ の非零集合は一点で交わる。これを :term[Nehari 多様体]{en="Nehari manifold" yomi="ねはりたようたい"} という。そこで $DG(u)u=2A-pB=(2-p)A\ne0$ である。拘束上の停留点は $J'(u)=\lambda DG(u)$ を満たすが、$u$ に作用させれば $0=\lambda(2-p)A$、従って $\lambda=0$ となる。このためこの拘束で得た停留点は元の方程式の解である。

:::example[一つの方向での峠の高さ]
$A=6,B=2,p=4$ なら $t=\sqrt3$ で最大値 $J(tu)=9/2$ である。$t>\sqrt6$ では $J(tu)<0$ になる。零からその点まで直線をたどる経路の最大値は $9/2$ なので、多数の経路の最小を取る峠の値に対しては上界を与える。終点を同じ半直線上の負の点に選んだ場合の上界である。別の方向を選ぶと $A^2/(4B)$ が変わるため、一方向の最大値を全経路の最小値と同一視できない。
:::

:::example[一つの方向での峠と負の領域]
変分方向に沿ったエネルギーを $J(a)=4a^2-a^4$ とする。非負の $a$ だけを考え、零から負のエネルギーまで到達する。

まず候補の選択を一変数 $a$ の問題にする。停留条件は $J_a=0$ であり、境界や無限遠も許される場合は停留点と分けて扱う。次の式では一階微分を取って零になる点を求め、二階微分を残しておく。

$$
\begin{aligned}
J(a)&=- a^{4} + 4 a^{2}\\
J_a&=\frac{d}{d a} \left(- a^{4} + 4 a^{2}\right)=- 4 a \left(a^{2} - 2\right)\\
J_{aa}&=\frac{d^{2}}{d a^{2}} \left(- a^{4} + 4 a^{2}\right)=- 4 \left(3 a^{2} - 2\right)\\
a_*&=\sqrt{2}\\
J_a(a_*)&=0,\qquad a=\sqrt{2}
\end{aligned}
$$

停留値と、その点から $h$ だけ動かした差を計算する。一次項が相殺しているため、展開の先頭から局所の振る舞いが分かる。一方、停留点が全体で最小か最大かを結論するには、許容範囲全体でこの差の符号を調べる必要がある。

$$
\begin{aligned}
J(a_*)&=4,\qquad a=\sqrt{2}\\
J_{aa}(a_*)&=-16\\
J(a_*+h)-J(a_*)&=- \left(h + \sqrt{2}\right)^{4} + 4 \left(h + \sqrt{2}\right)^{2} - 4=- h^{2} \left(h^{2} + 4 \sqrt{2} h + 8\right)\\
\text{二次係数}&=-8
\end{aligned}
$$

指定した値で目的関数を比較する。$a$ 自体の値と、そこでの目的関数の値は別の量である。

$$
\begin{aligned}
a_*\text{ の具体値}&=\sqrt{2}\\
J_*\text{ の具体値}&=4\\
J&=3,\qquad a=1\\
J&=-45,\qquad a=3\\
\text{二つの値の差}&=-48\\
\text{停留点まわりの差}&=- h^{2} \left(h^{2} + 4 \sqrt{2} h + 8\right)
\end{aligned}
$$

$J(a)-J(\sqrt2)=-(a^2-2)^2$ であり、非負の領域で最大点は $\sqrt2$ である。$a=3$ での値は負なので、零から $3$ への連続経路は中間値の定理から必ず $\sqrt2$ を通る。従ってこの一次元経路の峠の値は $4$ である。元の無限次元問題では経路がこの方向を離れられるため、これは一般には上界だけを与える。
:::

## Noether の定理

エネルギーが連続的な変換の下で変わらないなら、その微小変換を変分に使うことができる。内部での積分の相殺を点ごとの発散式に戻すと保存量が得られる。

:::proposition[変分対称性の発散形]{#prop-pde-noether}
独立変数を速度 $X(x)$、従属変数を速度 $\eta(x,u)$ で変える微小変換を考える。グラフの変化を同じ座標で比べた変分を $Q=\eta-X\cdot Du$ とする。任意の $u$ に対しエネルギー密度の微小変化が $\divg B$ であるとき、滑らかな停留点では

$$
\divg(F_pQ+FX-B)=0
\label{eq-pde-noether}
$$

となる。
:::

:::proof
固定座標での変化から $F_zQ+F_p\cdot DQ$ が生じる。さらに移動する積分の座標変換から $DF\cdot X+F\divg X=\divg(FX)$ が生じる。例えば Jacobian の一階項は $\divg X$ であるため、積の法則でこの形になる。従って $\divg B=(F_z-\divg F_p)Q+\divg(F_pQ+FX)$ である。停留点で最初の項が :eqref[eq-pde-euler-lagrange] により消えるので、この :term[Noether の定理]{en="Noether theorem" yomi="ねーたーのていり"} の形を得る。
:::

:::example[値の平行移動と空間の平行移動]
$F=|Du|^2/2$ に対する $u\mapsto u+s$ では $X=0$、$Q=1$、$B=0$ なので保存流束は $Du$ であり $\Delta u=0$ になる。空間の第 $k$ 方向の平行移動では $X=e_k,\ Q=-u_{x_k}$ であり、流束は $F e_k-u_{x_k}Du$ になる。$u=x_1^2-x_2^2$、$k=1$ では流束が $(2(x_2^2-x_1^2),4x_1x_2)$ で、二成分の対応する偏微分は $-4x_1$ と $4x_1$ だから発散は $0$ である。
:::

:::example[一次元の第一積分]
$F(u,u')=(u')^2/2+u^2/2$ は $x$ に陽に依存しないので、平行移動から $F-u'F_{u'}=(u^2-(u')^2)/2$ が一定になる。Euler–Lagrange 方程式は $u''=u$ である。$u=e^x$ ではこの定数は $0$、$u=\cosh x$ では $\cosh^2x-\sinh^2x=1$ より $1/2$ である。
:::

### 回転と応力の流束

$F=|Du|^2/2$ では、平行移動に対応する流束を列に並べると $T_{ik}=F\delta_{ik}-u_i u_k$ となる。積を微分すれば $\sum_i\partial_iT_{ik}=-(\Delta u)u_k$ が成り立つ。回転の速度場 $X(x)=Kx$ は反対称行列 $K$ で与えられ、$Q=-X\cdot Du$ なので Noether の流束は $TX$ である。実際 $\divg(TX)=\sum_{i,k}(\partial_iT_{ik})X_k+\sum_{i,k}T_{ik}K_{ki}$ であり、調和解では前者が零、後者も対称行列と反対称行列の縮約なので零である。

:::example[回転に対応する保存流束]
二次元で $u=x^2-y^2$、$X=(-y,x)$ とすると $Q=4xy$ である。$F=2(x^2+y^2)$ を使うと流束は $(6x^2y-2y^3,\ 2x^3-6xy^2)$ となる。その発散は $12xy-12xy=0$ である。従って円周を横切る全流束はどの半径でも零であり、角度ごとの流束は正負に変わり得る。
:::

### 時空間の平行移動とエネルギー

波の作用密度を $L(u_t,Du)=u_t^2/2-c^2|Du|^2/2$ とする。Euler–Lagrange 方程式は $u_{tt}-c^2\Delta u=0$ である。時間の平行移動に対して $u_t$ を掛けると
$$
\partial_t\big(u_t^2/2+c^2|Du|^2/2\big)
-\divg(c^2u_tDu)=0 .
$$
この保存流は $L$ の時間方向の平行移動から得られる。空間の第 $i$ 方向への移動では $u_i$ を掛ける。積の微分で
$$
\partial_t(u_tu_i)-\divg(c^2 u_i Du)
-\partial_i\big(u_t^2/2-c^2|Du|^2/2\big)=0
$$
となる。最後に密度 $L$ の導関数が残るため、$u_t$ を $u_i$ へ置き換えるだけでは流束を正しく読めない。周期領域で積分すると空間微分は消え、時間平行移動の積分からエネルギー、空間平行移動から各運動量が一定になる。

:::example[調和エネルギーの境界表現と回転]
円板で $g=\cos\theta+2\sin2\theta$ を固定し Dirichlet エネルギーの停留点を考える。

円周データを定数と三角関数に分ける。角振動数 $j$ について半径の方程式は $R^{\prime\prime}+r^{-1}R^{\prime}-j^2r^{-2}R=0$ である。元の条件に合う根を選び、境界半径で係数を規格化する。この例での関数は次のとおりである。

$$
\begin{aligned}
u(r,\theta)&=2 r^{2} \sin{\left(2 \theta \right)} + r \cos{\left(\theta \right)}\\
g&=2 \sin{\left(2 \theta \right)} + \cos{\left(\theta \right)}=\left(4 \sin{\left(\theta \right)} + 1\right) \cos{\left(\theta \right)},\qquad r=1\\
u_r&=\frac{\partial}{\partial r} \left(2 r^{2} \sin{\left(2 \theta \right)} + r \cos{\left(\theta \right)}\right)=\left(8 r \sin{\left(\theta \right)} + 1\right) \cos{\left(\theta \right)}\\
u_\theta&=\frac{\partial}{\partial \theta} \left(2 r^{2} \sin{\left(2 \theta \right)} + r \cos{\left(\theta \right)}\right)=r \left(4 r \cos{\left(2 \theta \right)} - \sin{\left(\theta \right)}\right)
\end{aligned}
$$

二階微分を別々に計算する。半径方向の二階微分だけを角度方向の二階微分と足しても相殺しない。極座標の Laplacian の一階項も入れると各振動数で相殺する。

$$
\begin{aligned}
u_{rr}&=\frac{\partial^{2}}{\partial r^{2}} \left(2 r^{2} \sin{\left(2 \theta \right)} + r \cos{\left(\theta \right)}\right)=4 \sin{\left(2 \theta \right)}\\
u_{\theta\theta}&=\frac{\partial^{2}}{\partial \theta^{2}} \left(2 r^{2} \sin{\left(2 \theta \right)} + r \cos{\left(\theta \right)}\right)=- r \left(16 r \sin{\left(\theta \right)} + 1\right) \cos{\left(\theta \right)}\\
\Delta u&=4 \sin{\left(2 \theta \right)} + \frac{4 r \sin{\left(2 \theta \right)} + \cos{\left(\theta \right)}}{r} - \frac{8 r \sin{\left(2 \theta \right)} + \cos{\left(\theta \right)}}{r}=0
\end{aligned}
$$

境界における値の平均と法線微分を読む。$j\ne0$ の角方向の平均は零であり、全体の平均は定数項である。またエネルギーは体積積分の方法と Green の境界積分の方法の両方で計算できる。

$$
\begin{aligned}
\text{円周の平均}&=\int\limits_{- \pi}^{\pi} \frac{2 \sin{\left(2 \theta \right)} + \cos{\left(\theta \right)}}{2 \pi}\, d\theta=0\\
\partial_\nu u&=4 \sin{\left(2 \theta \right)} + \cos{\left(\theta \right)}=\left(8 \sin{\left(\theta \right)} + 1\right) \cos{\left(\theta \right)}\\
\int_{-\pi}^{\pi}|Du|^2\,\mathrm{d}\theta&=32 \pi r^{2} + 2 \pi=2 \pi \left(16 r^{2} + 1\right)\\
\int |Du|^2&=\int\limits_{0}^{1} r \left(32 \pi r^{2} + 2 \pi\right)\, dr=9 \pi\\
\int_{\partial\Omega} u\partial_\nu u&=\int\limits_{- \pi}^{\pi} \left(2 \sin{\left(2 \theta \right)} + \cos{\left(\theta \right)}\right) \left(4 \sin{\left(2 \theta \right)} + \cos{\left(\theta \right)}\right)\, d\theta=9 \pi
\end{aligned}
$$

次に同じ振動数の半径因子を比較し、角度を固定した一点で値と半径方向の傾きを調べる。角度の取り方によって、振幅のあるモードでも点での値が零になることがある。

$$
\begin{aligned}
u&=0,\qquad r=\frac{1}{2},\theta=\frac{\pi}{2}\\
u_r&=0,\qquad r=\frac{1}{2},\theta=\frac{\pi}{2}\\
u_\theta&=- \frac{3}{2},\qquad r=\frac{1}{2},\theta=\frac{\pi}{2}
\end{aligned}
$$

円板と境界値を同じ角度だけ回転すると各振動数の正弦・余弦の係数が直交変換され、その二乗和は変わらない。従って境界値が動いても、同じ回転を領域と解に施すとエネルギーは一定である。固定した境界値の問題における内部の変分では、境界による仕事がない試験方向を選ぶ必要がある。
:::

:::example[振動子の時間並進対称性を直接確かめる]
$J(q)=\int(q_t^2-9q^2)/2\dd t$ に対し、端点で零になる変分から $q_{tt}+9q=0$ を得る。初期値 $(q(0),q_t(0))=(1,6)$ を置く。

この例は一自由度なので、質量係数と剛性係数はともにスカラーである。時間二階微分と剛性項の和を計算して、候補が Euler–Lagrange 方程式を満たすことを確認する。

$$
\begin{aligned}
M&=\left[\begin{matrix}1\end{matrix}\right]\\
K&=\left[\begin{matrix}9\end{matrix}\right]\\
U&=\left[\begin{matrix}2 \sin{\left(3 t \right)} + \cos{\left(3 t \right)}\end{matrix}\right]\\
U_{tt}&=\left[\begin{matrix}- 9 \left(2 \sin{\left(3 t \right)} + \cos{\left(3 t \right)}\right)\end{matrix}\right]\\
F&=\left[\begin{matrix}0\end{matrix}\right]
\end{aligned}
$$

初期時刻へ戻して、先に指定したデータと比較する。次に正の時刻における成分の値を読む。成分が混ざっていても、初期データは方程式を積分する前の元の座標で指定される。

$$
\begin{aligned}
U(0)&=\left[\begin{matrix}1\end{matrix}\right],\qquad t=0\\
U_t(0)&=\left[\begin{matrix}6\end{matrix}\right],\qquad t=0\\
U&=\left[\begin{matrix}2\end{matrix}\right],\qquad t=\frac{\pi}{6}\\
U_t&=\left[\begin{matrix}-3\end{matrix}\right],\qquad t=\frac{\pi}{6}
\end{aligned}
$$

二次量を具体的に求め、微分した値と流束または外力の寄与を照合する。質量行列が単位行列でない場合は係数の普通の二乗和がエネルギーになるとは限らない。

$$
\begin{aligned}
e&=\frac{45}{2}\\
e_t&=\frac{d}{d t} \frac{45}{2}=0\\
e_t&=0
\end{aligned}
$$

時間をずらしても Lagrangian の形は同じなので、対応する保存量は $q_t^2/2+9q^2/2$ である。計算した時刻では変位は $2$、速度は $-3$ になり、変位と速度の二項の入れ替わりを確認できる。時間依存のばね定数なら $K_t q^2/2$ の項が残り、保存性は失われる。
:::

### 発散だけ変わる対称性

場の変分 $u\mapsto u+sQ(x,u)+o(s)$ に対し、作用密度の変分が厳密に零である代わりに発散 $\divg B$ に等しい場合を考える。連鎖律と積の微分で
$$
\delta L=L_u Q+L_p\cdot DQ=(L_u-\divg L_p)Q+\divg(QL_p)
$$
である。従って Euler–Lagrange 方程式の解に沿って $\divg(QL_p-B)=0$ となる。密度が点ごとに不変という強い仮定を外しても、境界項を差し引けば保存量は成り立つ。ただし境界からの作用を固定する問題では、変分がその境界条件を保つかどうかも別に確認する。

:::example[一定の外力を持つ粒子の位置並進]
$L(q,q')=(q')^2/2-Fq$、$F$ を定数とする。$q\mapsto q+s$ で $Q=1$、$\delta L=-F$ は $B(t)=-Ft$ の微分である。よって保存量は $q'+Ft$ である。実際 Euler–Lagrange 方程式は $q''=-F$ なのでその微分は零になる。例えば $q(0)=1,q'(0)=4,F=3$ なら $q=1+4t-3t^2/2$ であり、$t=2$ の速度は $-2$、保存量は $4$ のままである。$q'$ 自体は保存しないが、作用密度の変化分と合わせると一定になる。
:::

## 演習問題

::::exercise{level=basic}
$J[u]=\int_0^1((u')^2+3u)\,\dd x$ の零境界での停留点を求めよ。

:::solution
第一変分は $\int(2u'v'+3v)$ で、内部方程式は $-2u''+3=0$ である。両端を $0$ にして $u=3x(x-1)/4$ を得る。二次の勾配項が狭義凸なのでこの停留点は最小である。
:::
::::

::::exercise{level=basic}
自由端で $F(p,z)=p^2/2+z^2/2$ の自然境界条件を求めよ。

:::solution
一階の偏微分は $F_p=p$ である。境界項は $[u'v]_0^1$ で、$v$ の端点値を独立に選べる。従って $u'(0)=u'(1)=0$ が必要十分である。
:::
::::

::::exercise{level=basic}
$(0,\pi)$ で $J[u]=\tfrac12\int((u')^2-2u^2)$ の $v=\sin3x$ 方向の第二変分を求めよ。

:::solution
このモードの固有値は $9$ なので第二変分は $(9-2)\int\sin^2(3x)=7\pi/2$ である。
:::
::::

::::exercise{level=basic}
両端の値が $0,3$ の区間 $(0,1)$ で $\int (u')^4/4$ の最小値を求めよ。

:::solution
平均勾配は $3$ で、Jensen の不等式から値は $3^4/4$ 以上である。$u=3x$ で等号なので最小値は $81/4$ である。
:::
::::

::::exercise{level=advanced}
$\int_0^1 u=2$、両端 $0$ のもとで $\frac12\int_0^1(u')^2$ を最小化せよ。

:::solution
弱方程式による候補 $u=\lambda x(1-x)/2$ の積分は $\lambda/12$ である。$\lambda=24$ なので $u=12x(1-x)$ である。摂動 $v$ の平均が零なら差は $\int(v')^2/2$ なので候補が最小で、値は $\frac12\,144\int(1-2x)^2=24$ である。
:::
::::

::::exercise{level=advanced}
$J(s)=s^2/2-s^4/4$ を $\R$ 上で考え、$0$ と $2$ を結ぶ峠の値を求めよ。

:::solution
連続な道は中間値の定理で $s=1$ を通り、$J(1)=1/4$ なので下界は $1/4$ である。直線の道が通る $[0,2]$ で $J'=s-s^3$ は $0<s<1$ で正、$s>1$ で負なので最大は $1/4$ である。従って峠の値は $1/4$ で、$s=1$ は臨界点である。
:::
::::

::::exercise{level=advanced}
$F(p)=|p|^4/4$ に対する一次元の空間平行移動の保存量を求めよ。

:::solution
$F_p=p^3$、$Q=-u'$ なので $F+F_pQ=-3(u')^4/4$ である。停留解では $((u')^3)'=0$ なので $u'$ は一定となり、流束の微分が $0$ である。
:::
::::

::::exercise{level=advanced}
$u(x)=3(1-|x|^{4/3})/4$ の二階弱微分が $L^q(-1,1)$ に入る $q\ge1$ の範囲を求めよ。

:::solution
一階導関数の連続性により、原点の左右で部分積分した端点項は相殺する。原点外の二階導関数は局所可積分なので、その関数が二階弱微分である。値は $-|x|^{-2/3}/3$ であるから $\int_0^1 x^{-2q/3}$ が有限となる $q<3/2$ が必要十分である。
:::
::::

## 確認テスト

::::quiz
境界値を固定した変分では試験関数の境界値はどうなるか。

- [x] 零である
- [ ] 常に $1$ である
- [ ] 任意である

:::explanation
$u+sv$ が全て同じ境界値を保つには $v$ の境界値が $0$ でなければならない。このとき部分積分の境界項は消える。
:::
::::

::::quiz{answer="12" tol="1e-6"}
$(0,1)$ の両端 $0$、$\int_0^1u=1$ の最小化で $-u''=\lambda$ となるとき $\lambda$ を求めよ。

:::explanation
端点条件から $u=\lambda x(1-x)/2$ で、その積分は $\lambda/12$ である。積分値 $1$ より $\lambda=12$ である。
:::
::::

::::quiz
直接法で弱収束から最小値を得るために必要な性質はどれか。

- [x] 弱逐次下半連続性
- [ ] 常に二階微分可能であること
- [ ] 最小化列の勾配が一様に連続すること

:::explanation
弱極限 $u$ で $J(u)\le\liminf J(u_n)$ が成り立てば、最小化列の値の極限以上に上がらない。弱閉性で $u$ が許容されることと組み合わせる。
:::
::::

::::quiz{answer="0.25" tol="1e-6"}
$J(s)=s^2/2-s^4/4$ で $0$ から $2$ に至る峠の値を求めよ。

:::explanation
すべての道が $1$ を通るので $J(1)=1/4$ が下界である。区間 $[0,2]$ の最大が $1$ で生じるため、直線の道がこの下界を達成する。
:::
::::

::::quiz{answer="4" tol="1e-6"}
$J(a)=4a^2-a^4$ の $a\ge0$ での最大値を求めよ。

:::explanation
$J(a)=4-(a^2-2)^2$ と平方完成できる。最大値は $a=\sqrt2$ で達成され $4$ である。$a=0$ も臨界点だが、その点の値は $0$ である。
:::
::::


---
title: Lagrange 形式
description: d'Alembert の原理と Hamilton の原理から一般化座標の運動方程式を導き、未定乗数による拘束力の計算、作用の最小性、Lagrangian の任意性と Galilei 不変性、電磁場中の荷電粒子、散逸を扱う。
notation:
  - symbol: 'q^i,\ \dot q^i'
    meaning: 第 i 一般化座標とその時間微分
  - symbol: 'T,\ V'
    meaning: 運動エネルギーと位置エネルギー
  - symbol: 'L'
    meaning: 運動を記述する Lagrangian
  - symbol: 'p_i'
    meaning: 第 i 一般化座標に共役な運動量
  - symbol: 'E_i(L)'
    meaning: Lagrangian L に対する Euler–Lagrange 式の左辺
  - symbol: 'S'
    meaning: 時間区間にわたる作用
  - symbol: 'U(q,\dot q,t)'
    meaning: 速度に依存する一般化ポテンシャル
  - symbol: '\phi,\ \vb{A}'
    meaning: 電磁場のスカラーポテンシャルとベクトルポテンシャル
  - symbol: '\mathcal R'
    meaning: 速度に依存する Rayleigh 散逸関数
---

## この章の目標

- d'Alembert の原理から理想拘束の反力を消去して Lagrange 方程式を導き、それが一般化座標の初期値から運動を一意に決めることを示す。
- 拘束力が必要な場合や非ホロノミック拘束がある場合に、未定乗数を含む Lagrange 方程式を立てて解く。
- Hamilton の原理と Lagrange 方程式の同値性を示し、作用が実際に最小となる条件を調べる。
- Lagrangian の全時間微分の不定性と点変換に対する共変性を証明し、Galilei 不変性から自由粒子の Lagrangian を定め、回転座標系の慣性力を導く。
- 一般化ポテンシャルとして電磁場中の荷電粒子を扱い、散逸関数で粘性抵抗とエネルギー収支を記述する。

二重振り子の軌道を予測するには二本の棒の角度がわかればよい。しかし直交座標の運動方程式には、両方の棒から受ける力が現れる。
この章では、位置を初めから角度で表し、運動エネルギーと位置エネルギーを微分して二つの角度の方程式を求める。
以下では質量一定の非相対論的な質点系を慣性系で扱い、和はすべて明示する。一般化座標そのものは時間に依存する座標系から選んでもよい。
標準的な扱いは :cite[goldstein-cm]、:cite[landau-mechanics]、:cite[onuki-am] にある。

## d'Alembert の原理から Lagrange 方程式へ

### 実際の加速度を仮想変位に射影する

質点 $a$ の質量を $m_a$、位置を $\vb{r}_a$、速度を $\vb{v}_a$ とする。
印加力 $\vb{F}_a$ と拘束力 $\vb{R}_a$ を分けると、Newton 方程式は $m_a\dot{\vb{v}}_a=\vb{F}_a+\vb{R}_a$ である。
[Newton 力学と拘束・一般化座標](/physics/analytical-mechanics/constraints/) の :ref[def-am-ideal-constraint] では、許されたどの仮想変位に対しても拘束力の仮想仕事の和が零であった。
そのため動いている系でも、時刻を固定して運動方程式と仮想変位の内積をとれば、反力を消去できる。

:::principle[d'Alembert の原理]{#prin-am-dalembert}
理想拘束の下を運動する $N$ 質点について、各時刻に許される任意の仮想変位に対し

$$
\sum_{a=1}^N(\vb{F}_a-m_a\dot{\vb{v}}_a)\cdot\delta\vb{r}_a=0.
\label{eq-am-dalembert}
$$

これを **d'Alembert の原理**という。
:::

実際、括弧内は $-\vb{R}_a$ に等しいのでこの式が従う。静止して加速度も零ならば :ref[prin-am-virtual-work] のつり合い条件に戻る。
$-m_a\dot{\vb{v}}_a$ を「慣性力」とみなせば、動力学を静力学のつり合いの形に書いたものと読める。
以下の導出では、さらに拘束はホロノミックで、独立な座標 $q^1,\ldots,q^n$ によって局所的に解けていると仮定する。
すなわち $\vb{r}_a=\vb{r}_a(q,t)$ は $C^2$ 級で、座標に関する微分の階数は $n$ である。非ホロノミック拘束は後の小節で扱う。

### 速度の二つの微分公式

式 :eqref[eq-am-coordinate-velocity] によれば、直交座標で測った各質点の速度は

$$
\vb{v}_a(q,\dot q,t)=
\sum_{j=1}^n\pdv{\vb{r}_a}{q^j}\dot q^j+\pdv{\vb{r}_a}{t}
$$

となる。軌道上では $\dot q^j$ は $q^j(t)$ の導関数であるが、
速度を $(q,\dot q,t)$ の関数として偏微分するときには、これらの引数を独立に扱う。

:::lemma[位置と速度の微分公式]{#lem-am-velocity-identities}
上の仮定の下で、どの $a,i$ に対しても

$$
\pdv{\vb{v}_a}{\dot q^i}=\pdv{\vb{r}_a}{q^i},\qquad
\pdv{\vb{v}_a}{q^i}=\dv{}{t}\pdv{\vb{r}_a}{q^i}.
\label{eq-am-velocity-identities}
$$

ここで右の式の右辺は、軌道 $q(t)$ に沿った全時間微分である。
:::

:::proof
第一式は、$\vb{v}_a$ のうち $\dot q^i$ にかかる係数を読むと得られる。
第二式では $q^i$ による偏微分と軌道に沿った全時間微分をそれぞれ展開する。

$$
\begin{aligned}
\pdv{\vb{v}_a}{q^i}
&=\sum_{j=1}^n
\frac{\partial^2\vb{r}_a}{\partial q^i\partial q^j}\dot q^j
+\frac{\partial^2\vb{r}_a}{\partial q^i\partial t},\\
\dv{}{t}\pdv{\vb{r}_a}{q^i}
&=\sum_{j=1}^n
\frac{\partial^2\vb{r}_a}{\partial q^j\partial q^i}\dot q^j
+\frac{\partial^2\vb{r}_a}{\partial t\partial q^i}.
\end{aligned}
$$

$C^2$ 級なので混合偏導関数を交換でき、二つの右辺は等しい。最後の項を含めたため、支点などが時間とともに動く場合にも成り立つ。
:::

### 慣性項を運動エネルギーで表す

全運動エネルギーを $T=\tfrac12\sum_{a=1}^N m_a\vb{v}_a\cdot\vb{v}_a$ とする。
第 $i$ 座標の方向に射影した慣性項に積の微分則を適用し、:ref[lem-am-velocity-identities] を使うと

$$
\begin{aligned}
\sum_{a=1}^N m_a\dot{\vb{v}}_a\cdot\pdv{\vb{r}_a}{q^i}
&=\dv{}{t}\sum_{a=1}^N m_a\vb{v}_a\cdot\pdv{\vb{r}_a}{q^i}
-\sum_{a=1}^N m_a\vb{v}_a\cdot\dv{}{t}\pdv{\vb{r}_a}{q^i}\\
&=\dv{}{t}\sum_{a=1}^N m_a\vb{v}_a\cdot\pdv{\vb{v}_a}{\dot q^i}
-\sum_{a=1}^N m_a\vb{v}_a\cdot\pdv{\vb{v}_a}{q^i}\\
&=\dv{}{t}\pdv{T}{\dot q^i}-\pdv{T}{q^i}.
\end{aligned}
$$

一方、$\delta\vb{r}_a=\sum_{i=1}^n(\partial\vb{r}_a/\partial q^i)\delta q^i$ を :eqref[eq-am-dalembert] に代入すると、
印加力の項は :eqref[eq-am-generalized-force] の一般化力 $Q_i$ で表せる。したがって

$$
\sum_{i=1}^n
\left(Q_i-\dv{}{t}\pdv{T}{\dot q^i}+\pdv{T}{q^i}\right)\delta q^i=0.
$$

:ref[prop-am-virtual-space] により各 $\delta q^i$ を独立に選べるため、すべての $i$ について括弧内が零である。
この段階では印加力が保存力であることを仮定していない。

位置に依存するポテンシャル $V(\vb{r}_1,\ldots,\vb{r}_N,t)$ が与えられ、その力が $-\nabla_a V$ であるとする。
$V$ を座標 $q$ で表し直すと、連鎖律より
$-\sum_{a=1}^N\nabla_a V\cdot(\partial\vb{r}_a/\partial q^i)=-\partial V/\partial q^i$ となる。
ポテンシャルに含めなかった力の一般化力を $Q_i^{\mathrm{nc}}$ と書けば、$Q_i=-\partial V/\partial q^i+Q_i^{\mathrm{nc}}$ である。

:::definition[Lagrangian と一般化運動量]{#def-am-lagrangian}
系の運動方程式を定める関数 $L(q,\dot q,t)$ を :term[Lagrangian]{en="Lagrangian" yomi="らぐらんじあん"}（Lagrange 関数）という。
上の力学系では

$$
L=T-V,\qquad p_i=\pdv{L}{\dot q^i}
$$

とおく。$p_i$ を $q^i$ に共役な :term[一般化運動量]{en="generalized momentum" yomi="いっぱんかうんどうりょう"} という。
:::

$V$ が速度に依存しないことから、上で導いた $T$ の方程式は次の :term[Lagrange 方程式]{en="Lagrange equations" yomi="らぐらんじゅほうていしき"}（第二種）になる。

$$
\dv{}{t}\pdv{L}{\dot q^i}-\pdv{L}{q^i}=Q_i^{\mathrm{nc}},
\qquad i=1,\ldots,n.
\label{eq-am-lagrange-equations}
$$

左辺を $E_i(L)$ と書く。これは前章の Euler–Lagrange 方程式 :eqref[eq-am-vector-variation] の左辺と同じ形であり、$x$ を $t$ に、$y^i$ を $q^i$ に読み替えたものである。
力をすべて $L$ に含めた場合には右辺は零である。力を $L$ に含める場合と $Q_i^{\mathrm{nc}}$ に入れる場合のどちらかを選び、二重に数えない。
反力が必要なら、解 $q(t)$ から各 $\ddot{\vb{r}}_a$ を計算し、Newton 方程式に戻って求められる。後の小節では反力を直接求める方法を述べる。

Lagrange 方程式が運動を決めることを確かめよう。

:::proposition[Lagrange 方程式の初期値問題]{#prop-am-lagrange-ode}
$L=T-V$ とし、一般化座標の定義式 $\vb{r}_a(q,t)$ は $C^3$ 級、$V$ は $C^2$ 級、$Q^{\mathrm{nc}}_i(q,\dot q,t)$ は $C^1$ 級とする。
このとき :eqref[eq-am-lagrange-equations] は

$$
\sum_{j=1}^nM_{ij}(q,t)\ddot q^j=\Phi_i(q,\dot q,t)
$$

の形に書ける。ここで $M_{ij}$ は :eqref[eq-am-kinetic-generalized] の質量行列、$\Phi_i$ は $C^1$ 級の関数である。
したがって正則な座標の範囲では $\ddot q=M^{-1}\Phi$ と解け、初期値 $(q(t_0),\dot q(t_0))$ に対して解は存在する限り一意である。
:::

:::proof
:eqref[eq-am-kinetic-generalized] から $\partial T/\partial\dot q^i=\sum_jM_{ij}\dot q^j+b_i$ である。これを軌道に沿って微分すると

$$
\dv{}{t}\pdv{T}{\dot q^i}=\sum_{j=1}^nM_{ij}\ddot q^j+\sum_{j,k}\pdv{M_{ij}}{q^k}\dot q^k\dot q^j+\sum_j\pdv{M_{ij}}{t}\dot q^j+\sum_k\pdv{b_i}{q^k}\dot q^k+\pdv{b_i}{t}
$$

となり、$\ddot q$ は第一項にだけ現れる。残りの項と $\partial T/\partial q^i$、$\partial V/\partial q^i$、$Q^{\mathrm{nc}}_i$ をまとめたものが $\Phi_i$ である。
$M,b,T_0$ は $\vb{r}_a$ の一階偏導関数の積なので $C^2$ 級であり、$\Phi$ は $C^1$ 級となる。
:ref[prop-am-mass-matrix] により $M$ は正則で逆行列も $C^1$ 級なので、$\ddot q=M^{-1}\Phi$ は一階の連立方程式 $\dot q=v$、$\dot v=M^{-1}\Phi$ に書ける。
右辺は $C^1$ 級なので、常微分方程式の解の一意性定理（[解の存在と一意性](/math/ode/existence-uniqueness/)）が適用できる。
:::

これで、理想拘束の下の運動が $(q,\dot q)$ の初期値で一意に決まることが示された。:ref[def-am-equilibrium] の後の注意で予告した事実である。
一般に、$\det(\partial^2L/\partial\dot q^i\partial\dot q^j)\ne0$ となる Lagrangian を :term[正則な Lagrangian]{en="regular Lagrangian" yomi="せいそくならぐらんじあん"} という。正則な Lagrangian では同じ議論で方程式を $\ddot q$ について解ける。この条件は [Hamilton 形式](/physics/analytical-mechanics/hamiltonian/) への移行でも必要になる。

### 単振り子と二重振り子

固定支点から長さ $\ell$ の棒で質量 $m$ の質点をつるす。棒は軽く、支点での摩擦はなく、重力は一様とする。
鉛直下向きからの角 $\theta$ で位置を $(\ell\sin\theta,-\ell\cos\theta)$ と表すと、速さの二乗は $\ell^2\dot\theta^2$ である。
したがって

$$
L=\tfrac12m\ell^2\dot\theta^2+mg\ell\cos\theta,\qquad
p_\theta=m\ell^2\dot\theta,\qquad
\ddot\theta+\frac{g}{\ell}\sin\theta=0 .
$$

最後の式で小振幅近似は行っていない。棒の張力が消えたのは、許される変位が円の接線方向だからである。

次に、一つ目のおもりの先にもう一本の軽い棒と質点を取り付ける。
長さを $\ell_1,\ell_2$、質量を $m_1,m_2$ とし、各棒の角 $\theta_1,\theta_2$ はそれぞれ**鉛直下向きから**同じ向きに測る。
$\vb{e}(\theta)=(\sin\theta,-\cos\theta)$ と書くと、

$$
\vb{r}_1=\ell_1\vb{e}(\theta_1),\qquad
\vb{r}_2=\vb{r}_1+\ell_2\vb{e}(\theta_2).
$$

```tikz caption="二本の棒の角度は、それぞれの支点を通る鉛直線から測る。図では第一の角は正、第二の角は負である。"
\begin{tikzpicture}[>=Stealth]
\coordinate (O) at (0,0);
\coordinate (P) at (-52:2.2);
\coordinate (Q) at ($(P)+(-125:1.8)$);
\draw[dashed] (O) -- ++(0,-2);
\draw[dashed] (P) -- ++(0,-1.9);
\draw[thick] (O) -- node[above right] {$\ell_1$} (P);
\draw[thick] (P) -- node[left] {$\ell_2$} (Q);
\draw[->] (0,-0.65) arc (-90:-52:0.65);
\node at (-72:0.98) {$\theta_1$};
\draw[->] ($(P)+(0,-0.65)$) arc (-90:-125:0.65);
\node at ($(P)+(-106:1.04)$) {$\theta_2$};
\fill (O) circle (1.6pt) node[above] {$O$};
\fill (P) circle (2.5pt) node[right] {$m_1$};
\fill (Q) circle (2.5pt) node[below left] {$m_2$};
\end{tikzpicture}
```

$\vb{e}'(\theta_1)\cdot\vb{e}'(\theta_2)=\cos(\theta_1-\theta_2)$ なので、二つ目の質点の速さの二乗には交差項がある。
$M=m_1+m_2$、$K=m_2\ell_1\ell_2$、$\Delta=\theta_1-\theta_2$ とおけば、

$$
\begin{aligned}
T&=\tfrac12 M\ell_1^2\dot\theta_1^2+
\tfrac12 m_2\ell_2^2\dot\theta_2^2+
K\dot\theta_1\dot\theta_2\cos\Delta,\\
V&=-Mg\ell_1\cos\theta_1-m_2g\ell_2\cos\theta_2 .
\end{aligned}
$$

ここで結合の符号を確かめよう。$L=T-V$ の速度微分と角度微分は

$$
\begin{aligned}
p_1&=M\ell_1^2\dot\theta_1+K\dot\theta_2\cos\Delta,&
\pdv{L}{\theta_1}&=-K\dot\theta_1\dot\theta_2\sin\Delta-Mg\ell_1\sin\theta_1,\\
p_2&=m_2\ell_2^2\dot\theta_2+K\dot\theta_1\cos\Delta,&
\pdv{L}{\theta_2}&= K\dot\theta_1\dot\theta_2\sin\Delta-m_2g\ell_2\sin\theta_2 .
\end{aligned}
$$

$\dd(\cos\Delta)/\dd t=-(\dot\theta_1-\dot\theta_2)\sin\Delta$ を使うと、$\dot p_i-\partial L/\partial\theta_i$ で二つの角速度の積が相殺する。よって

$$
\begin{aligned}
M\ell_1^2\ddot\theta_1+K\cos\Delta\,\ddot\theta_2+
K\sin\Delta\,\dot\theta_2^2+Mg\ell_1\sin\theta_1&=0,\\
m_2\ell_2^2\ddot\theta_2+K\cos\Delta\,\ddot\theta_1-
K\sin\Delta\,\dot\theta_1^2+m_2g\ell_2\sin\theta_2&=0.
\end{aligned}
\label{eq-am-double-pendulum}
$$

正の質量と長さに対して、加速度の係数行列（質量行列）の行列式は
$m_2\ell_1^2\ell_2^2(m_1+m_2\sin^2\Delta)>0$ である。したがって同時刻の角度と角速度から二つの角加速度が一意に決まり、:ref[prop-am-lagrange-ode] の具体例になっている。
小振幅でなければこの方程式は初等関数で解けず、数値計算によれば、エネルギーが大きいと初期値にきわめて敏感に依存する運動を示す。このような非可積分系の振る舞いの一般論は [可積分系とカオス入門](/physics/analytical-mechanics/integrability-chaos/) で扱う。

### 球面振り子・回転する円環・Atwood の装置

:::example[球面振り子]{#exm-am-spherical-pendulum}
長さ $\ell$ の軽い棒の先の質量 $m$ の質点が、固定支点のまわりの球面上を動く。鉛直下向きから測った極角 $\theta$ と方位角 $\varphi$ を使うと、

$$
L=\tfrac12m\ell^2(\dot\theta^2+\sin^2\theta\,\dot\varphi^2)+mg\ell\cos\theta .
$$

Lagrange 方程式は

$$
m\ell^2\ddot\theta-m\ell^2\sin\theta\cos\theta\,\dot\varphi^2+mg\ell\sin\theta=0,\qquad
\dv{}{t}\big(m\ell^2\sin^2\theta\,\dot\varphi\big)=0
$$

である。$L$ は $\varphi$ を含まないので $p_\varphi=m\ell^2\sin^2\theta\,\dot\varphi$ が一定であり、これは鉛直軸まわりの角運動量である（[対称性と保存則](/physics/analytical-mechanics/symmetries/)）。
$\theta=\theta_0$（$0<\theta_0<\pi/2$）一定で回る円錐振り子の解は、第一式から $\dot\varphi^2=g/(\ell\cos\theta_0)$ を満たす。
:::

:::example[回転する円環上の小環]{#exm-am-rotating-hoop}
半径 $a$ の滑らかな鉛直円環を、その鉛直直径のまわりに一定角速度 $\Omega$ で回す。
円環に通した質量 $m$ の小環の角 $\theta$ を最下点から測る。慣性系での位置は

$$
\vb{r}=a(\sin\theta\cos\Omega t,\ \sin\theta\sin\Omega t,\ -\cos\theta).
$$

円に沿う接線方向と、円を回す方位方向とは直交する。そのため
$\vb{v}^2=a^2(\dot\theta^2+\Omega^2\sin^2\theta)$ となり、

$$
L=\tfrac12 ma^2(\dot\theta^2+\Omega^2\sin^2\theta)+mga\cos\theta .
$$

速度微分は $ma^2\dot\theta$、角度微分は $ma^2\Omega^2\sin\theta\cos\theta-mga\sin\theta$ である。したがって

$$
\ddot\theta=(\Omega^2\cos\theta-g/a)\sin\theta.
\label{eq-am-rotating-hoop}
$$

$\Omega=0$ では固定した円上の振り子に戻る。$\Omega^2>g/a$ になると最下点、最上点のほかに $\cos\theta=g/(a\Omega^2)$ を満たす二つの相対平衡が現れる。
ここでの相対平衡とは、$\theta$ が一定で円環と一緒に回る運動である。
拘束は非定常であり、円環を回す装置が小環に仕事をする。
:::

:::example[Atwood の装置]
伸びない軽い糸を、軸の摩擦と回転慣性が無視できる固定滑車に掛ける。両端に質量 $m_1,m_2>0$ のおもりをつけ、糸は張っているとする。
下向きを正とし、片方の高さ座標を $x$ とすると、もう一方は糸の長さから $c-x$ である（$c$ は定数）。
よって速さは等しく、重力の位置エネルギーのうち $x$ に依存する部分は $-(m_1-m_2)gx$ である。したがって

$$
L=\tfrac12(m_1+m_2)\dot x^2+(m_1-m_2)gx,\qquad
\ddot x=\frac{m_1-m_2}{m_1+m_2}g .
$$

重いほうが下がる符号になっている。$m_1=m_2$ の場合は、初速度が零なら静止し、初速度を与えれば等速で動く。
:::

### 未定乗数法：拘束力と非ホロノミック拘束

これまでは拘束をすべて解いた独立な座標を使った。しかし (i) 拘束力そのものを知りたい場合、(ii) 拘束を解くのが面倒な場合、(iii) 非ホロノミック拘束があって独立な座標で解けない場合には、座標の一部の関係を拘束として残したほうがよい。

座標 $q^1,\ldots,q^n$ は、一部のホロノミック拘束を解いて得られた正則な一般化座標とする。これに加えて、$m$ 個の独立な Pfaff 型拘束

$$
\sum_{i=1}^na_{\alpha i}(q,t)\dot q^i+a_{\alpha0}(q,t)=0,\qquad\alpha=1,\ldots,m
\label{eq-am-pfaffian-generalized}
$$

が課されているとする。ホロノミック拘束 $f_\alpha(q,t)=0$ は $a_{\alpha i}=\partial f_\alpha/\partial q^i$、$a_{\alpha0}=\partial f_\alpha/\partial t$ の場合である。
:eqref[eq-am-pfaffian-virtual] と同様に、仮想変位は $\sum_ia_{\alpha i}\delta q^i=0$ を満たす $\delta q$ とし、これらの拘束も理想的とする。

:::theorem[未定乗数を含む Lagrange 方程式]{#thm-am-lagrange-multipliers}
上の仮定の下で、実際の運動に対して関数 $\lambda_1(t),\ldots,\lambda_m(t)$ がただ一組存在して

$$
\dv{}{t}\pdv{L}{\dot q^i}-\pdv{L}{q^i}=Q_i^{\mathrm{nc}}+\sum_{\alpha=1}^m\lambda_\alpha a_{\alpha i},\qquad i=1,\ldots,n
\label{eq-am-lagrange-multiplier-equations}
$$

が成り立つ。$\sum_\alpha\lambda_\alpha a_{\alpha i}$ は、追加した拘束の拘束力の一般化力である。
:eqref[eq-am-lagrange-multiplier-equations] の $n$ 本と拘束 :eqref[eq-am-pfaffian-generalized] の $m$ 本が、$q$ の $n$ 成分と $\lambda$ の $m$ 成分を決める。
:::

:::proof
$q$ を座標とする段階の拘束（すでに解いたもの）の拘束力と、追加した拘束の拘束力は、どちらも理想的である。したがって d'Alembert の原理 :eqref[eq-am-dalembert] は、追加した拘束を満たすすべての仮想変位 $\delta q$ について成り立つ。
上の計算と同じく、それは

$$
\sum_{i=1}^n\big(E_i(L)-Q_i^{\mathrm{nc}}\big)\delta q^i=0\qquad\Big(\textstyle\sum_ia_{\alpha i}\delta q^i=0\ \ (\alpha=1,\ldots,m)\Big)
$$

と書ける。:ref[lem-am-annihilator] を行列 $A=(a_{\alpha i})$ に適用すると、ベクトル $(E_i(L)-Q_i^{\mathrm{nc}})_i$ は $\transp{A}\lambda$ の形にただ一通りに書ける。
後半の主張は、追加した拘束の拘束力の一般化力を $Q'_i$ とすると、その座標の段階での d'Alembert の原理が $E_i(L)=Q_i^{\mathrm{nc}}+Q'_i$ を与えることからわかる。
:::

直交座標 $q=X$ をそのまま使い、ホロノミック拘束 $f_\alpha(\vb{r}_1,\ldots,\vb{r}_N,t)=0$ をすべて乗数で扱えば、$L=\sum_a\tfrac12m_a\vb{v}_a^2-V$ に対して

$$
m_a\ddot{\vb{r}}_a=-\nabla_aV+\sum_{\alpha=1}^k\lambda_\alpha\nabla_af_\alpha,\qquad a=1,\ldots,N
\label{eq-am-lagrange-first-kind}
$$

となる。これを **Lagrange 方程式（第一種）** という。右辺第二項は :eqref[eq-am-constraint-force-multiplier] の拘束力そのものであり、:ref[prop-am-ideal-multiplier] と同じ符号の規約である。

:::proposition[乗数の決定]{#prop-am-multiplier-determined}
:eqref[eq-am-lagrange-first-kind] で $V$ と $f_\alpha$ は $C^2$ 級、拘束は独立とする。拘束を満たす位置と、:eqref[eq-am-admissible-velocity] を満たす速度に対し、$\lambda$ は $(\vb{r},\vb{v},t)$ の連続関数としてただ一通りに決まる。
その $\lambda$ を代入した方程式の解は、初めに拘束とその時間微分を満たせば以後も拘束を満たす。
:::

:::proof
拘束を二回時間微分すると

$$
\sum_{a=1}^N\nabla_af_\alpha\cdot\ddot{\vb{r}}_a+w_\alpha(\vb{r},\vb{v},t)=0
$$

である。ここで $w_\alpha$ は $\ddot{\vb{r}}$ を含まない項（$f_\alpha$ の二階偏導関数と速度の二次式など）をまとめたものである。
:eqref[eq-am-lagrange-first-kind] を代入すると $\sum_\beta G_{\alpha\beta}\lambda_\beta=\sum_a\nabla_af_\alpha\cdot\nabla_aV/m_a-w_\alpha$ を得る。ここで

$$
G_{\alpha\beta}=\sum_{a=1}^N\frac{1}{m_a}\nabla_af_\alpha\cdot\nabla_af_\beta
$$

である。$\sum_{\alpha,\beta}G_{\alpha\beta}\xi_\alpha\xi_\beta=\sum_a m_a^{-1}\abs{\sum_\alpha\xi_\alpha\nabla_af_\alpha}^2$ は、拘束の独立性により $\xi\ne0$ なら正である。よって $G$ は正則で、$\lambda$ は一意に決まる。
この $\lambda$ を代入した方程式では、構成から $\ddot f_\alpha=0$ が恒等的に成り立つ。初期に $f_\alpha=\dot f_\alpha=0$ なら $f_\alpha\equiv0$ である。
:::

:::example[振り子の張力]{#exm-am-pendulum-tension}
平面振り子を極座標 $(r,\theta)$（$\theta$ は鉛直下向きから）で表し、拘束 $f=r-\ell=0$ を乗数で扱う。
$L=\tfrac12m(\dot r^2+r^2\dot\theta^2)+mgr\cos\theta$、$a_r=1$、$a_\theta=0$ なので、

$$
m\ddot r-mr\dot\theta^2-mg\cos\theta=\lambda,\qquad
\dv{}{t}(mr^2\dot\theta)+mgr\sin\theta=0 .
$$

$r=\ell$ を代入すると第二式は単振り子の方程式になり、第一式から $\lambda=-m\ell\dot\theta^2-mg\cos\theta$ を得る。
拘束力は $\lambda\nabla f=\lambda\vb{e}_r$ であり、$\lambda<0$ は支点向きの張力を意味する。張力の大きさは $S=-\lambda=m\ell\dot\theta^2+mg\cos\theta$ である。
これは向心力 $m\ell\dot\theta^2$ と重力の動径成分の和として理解できる。
:::

:::example[斜面を転がる輪]{#exm-am-rolling-hoop}
質量 $M$、半径 $r$ の薄い輪が、傾斜角 $\alpha$ の斜面を滑らずに転がり下りる。斜面に沿って下向きに測った中心の位置を $x$、輪の回転角を $\vartheta$ とする。輪の中心まわりの慣性モーメントは $Mr^2$ である。
転がりの条件 $r\dot\vartheta-\dot x=0$ は積分可能だが、摩擦力を求めるために乗数で扱う（$a_x=-1$、$a_\vartheta=r$）。

$$
L=\tfrac12M\dot x^2+\tfrac12Mr^2\dot\vartheta^2+Mgx\sin\alpha
$$

から

$$
M\ddot x-Mg\sin\alpha=-\lambda,\qquad Mr^2\ddot\vartheta=\lambda r .
$$

拘束から $\ddot x=r\ddot\vartheta$ なので、第二式は $M\ddot x=\lambda$ となり、

$$
\ddot x=\tfrac12g\sin\alpha,\qquad \lambda=\tfrac12Mg\sin\alpha
$$

を得る。$x$ 方向の一般化力 $-\lambda$ は斜面上向きの静止摩擦力、$\vartheta$ 方向の $\lambda r$ はその力のモーメントである。
摩擦がなく滑り落ちる場合の加速度 $g\sin\alpha$ の半分になるのは、位置エネルギーの減少が並進と回転に等分されるためである。
:::

:::example[鉛直円板の転がり運動]{#exm-am-rolling-disk-motion}
:ref[exm-am-rolling-disk] の円板（質量 $M$）について運動方程式を立てる。円板は支持機構によって鉛直に保たれ、その機構の力は理想的とする。
自転軸まわりと鉛直な直径まわりの慣性モーメントを $I_1,I_2$ とすれば、重力の位置エネルギーは一定なので

$$
L=\tfrac12M(\dot x^2+\dot y^2)+\tfrac12I_1\dot\phi^2+\tfrac12I_2\dot\psi^2 .
$$

二つの拘束 $\dot x-a\cos\psi\,\dot\phi=0$、$\dot y-a\sin\psi\,\dot\phi=0$ の係数は $(1,0,-a\cos\psi,0)$、$(0,1,-a\sin\psi,0)$（変数の順は $x,y,\phi,\psi$）である。:eqref[eq-am-lagrange-multiplier-equations] は

$$
M\ddot x=\lambda_1,\qquad M\ddot y=\lambda_2,\qquad I_1\ddot\phi=-a(\lambda_1\cos\psi+\lambda_2\sin\psi),\qquad I_2\ddot\psi=0
$$

となる。拘束を微分した $\ddot x=a\cos\psi\,\ddot\phi-a\sin\psi\,\dot\psi\dot\phi$、$\ddot y=a\sin\psi\,\ddot\phi+a\cos\psi\,\dot\psi\dot\phi$ を代入すると、$\lambda_1\cos\psi+\lambda_2\sin\psi=Ma\ddot\phi$ より $(I_1+Ma^2)\ddot\phi=0$ を得る。
したがって $\dot\phi$ と $\dot\psi=\Omega$ は一定であり、$\Omega\ne0$ なら中心は半径 $a\abs{\dot\phi/\Omega}$ の円を一定の速さで回る。
接地点の摩擦力は $(\lambda_1,\lambda_2)=Ma\dot\phi\,\Omega(-\sin\psi,\cos\psi)$ で、この円運動の向心力を与える。
:::

## Hamilton の原理

### 経路を変えて作用を微分する

:::definition[作用]{#def-am-action}
時刻 $t_0<t_1$ と、両端の配位 $q(t_0),q(t_1)$ を固定する。経路 $q(t)$ に対し

$$
S[q]=\int_{t_0}^{t_1}L(q(t),\dot q(t),t)\,\dd t
\label{eq-am-action}
$$

を :term[作用]{en="action" yomi="さよう"} という。$L$ の単位は $\mathrm J$、$S$ の単位は $\mathrm{J\,s}$ である。
:::

:::principle[Hamilton の原理]{#prin-am-hamilton}
理想的なホロノミック拘束を独立な座標で解き、印加力をすべて $L$ に含められるとき、実際の運動の経路は、
両端の時刻と配位を固定した変形に対して作用を停留させる。すなわち $\delta S=0$ である。
:::

$L$ と位置表示が $C^2$ 級で実軌道も $C^2$ 級とする。端点で消える任意の $C^1$ 級関数 $\eta^i(t)$ により、
近傍の経路を $q_\eps^i=q^i+\eps\eta^i$ とおく。すべての近傍の経路でホロノミック拘束を満たすように一般化座標をとっている。
$\delta\dot q^i=\dot\eta^i$ を使って第一変分をとると、

$$
\begin{aligned}
\delta S
&=\int_{t_0}^{t_1}\sum_{i=1}^n
\left(\pdv{L}{q^i}\eta^i+\pdv{L}{\dot q^i}\dot\eta^i\right)\dd t\\
&=\left[\sum_{i=1}^np_i\eta^i\right]_{t_0}^{t_1}
+\int_{t_0}^{t_1}\sum_{i=1}^n
\left(\pdv{L}{q^i}-\dv{}{t}\pdv{L}{\dot q^i}\right)\eta^i\,\dd t.
\end{aligned}
$$

境界項は零である。よって :ref[lem-am-fundamental-variation] を一成分ずつ使えば、$\delta S=0$ は
:eqref[eq-am-lagrange-equations] の右辺を零とした式と同値である。これが同じ運動を経路の性質として述べられる理由である。
「最小作用の原理」という呼び方もあるが、ここで用いた条件が与えるのは停留であり、最小性には別の判定が要る（次小節）。

運動方程式を初期値問題として解く場合には、通常 $q(t_0),\dot q(t_0)$ を与える。
これに対し、上の証明では実軌道上の二点を選び、変形の間だけ両端をその位置に留めた。
したがって二つの端点を結ぶ運動が一意であるという仮定は使っていない。

ポテンシャルに含められない力 $Q_i^{\mathrm{nc}}$ がある場合も、上の計算から

$$
\delta S+\int_{t_0}^{t_1}\sum_{i=1}^nQ_i^{\mathrm{nc}}\eta^i\,\dd t=0
$$

が、端点を固定したすべての $\eta$ について成り立つことが :eqref[eq-am-lagrange-equations] と同値である。ただしこれは一つの汎関数の停留条件ではない。

Hamilton の原理の利点は、座標の選び方によらない形で運動法則を述べられることにある。作用は座標を変えても同じ値をもつので、停留条件も座標によらない（:ref[prop-am-point-covariance] で方程式の側から確かめる）。
さらに場の理論や量子力学の経路積分など、Newton 方程式から出発しにくい理論も同じ原理で定式化できる（[連続体と場の解析力学](/physics/analytical-mechanics/field-mechanics/)）。

:::remark[転がり拘束での変分]
非ホロノミック拘束がある場合、Hamilton の原理を素朴に適用すると誤った方程式を得る。
:ref[exm-am-rolling-disk] の拘束 $\dot x=a\cos\psi\,\dot\phi$ を例にとる。

1. 拘束を満たす経路の族 $(x_\eps,\ldots)$ を作って $\eps$ で微分すると、変分は $\delta\dot x=a\cos\psi\,\delta\dot\phi-a\sin\psi\,\dot\phi\,\delta\psi$ を満たす。
2. d'Alembert の原理で使う仮想変位は、各時刻で $\delta x=a\cos\psi\,\delta\phi$ を満たす。これを時間微分すると $\dv{}{t}\delta x=a\cos\psi\,\delta\dot\phi-a\sin\psi\,\dot\psi\,\delta\phi$ である。

最後の項が $\dot\phi\,\delta\psi$ と $\dot\psi\,\delta\phi$ で異なるので、二つの条件は一般に両立しない。
実際の転がり運動は d'Alembert の原理、すなわち 2 に従い、その結果が :ref[thm-am-lagrange-multipliers] の方程式である。
1 のように拘束を満たす経路だけで作用を停留させると、一般には転がる円板の運動とは異なる方程式が得られる。
ホロノミック拘束では二つの条件は一致し、このような食い違いは生じない。
:::

### 作用は最小か

前章 :ref[exm-am-jacobi-example] の汎関数は、調和振動子の作用に比例している。この例が示すように、作用は時間区間が長いと最小にならない。逆に、区間が短ければ最小になることを証明しよう。

:::proposition[短い時間区間では作用は最小]{#prop-am-short-time-minimum}
$q\in\R^n$ に対し $L=\tfrac12m\sum_{i=1}^n(\dot q^i)^2-V(q,t)$ とし、$V$ は $C^2$ 級で、定数 $K$ があってすべての $(q,t)$ と $\xi\in\R^n$ で

$$
\sum_{i,j}\frac{\partial^2V}{\partial q^i\partial q^j}(q,t)\,\xi^i\xi^j\le K\abs{\xi}^2
$$

を満たすとする。$q(t)$ を $[t_0,t_1]$ 上の Lagrange 方程式の $C^2$ 級の解とし、$\tau=t_1-t_0$ とおく。
$K\le0$、または $K>0$ かつ $\tau<\pi\sqrt{m/K}$ ならば、端点で零になる任意の $C^1$ 級関数 $\eta\ne0$ について $S[q+\eta]>S[q]$ である。
:::

:::proof
差を厳密に書くと

$$
S[q+\eta]-S[q]=\int_{t_0}^{t_1}\big(m\dot q\cdot\dot\eta-\grad V\cdot\eta\big)\dd t
+\int_{t_0}^{t_1}\Big(\tfrac12m\abs{\dot\eta}^2-\big[V(q+\eta)-V(q)-\grad V(q)\cdot\eta\big]\Big)\dd t
$$

である。第一項は部分積分により $\int(-m\ddot q-\grad V)\cdot\eta\,\dd t=0$ となる。
Taylor の定理の積分形 $V(q+\eta)-V(q)-\grad V(q)\cdot\eta=\int_0^1(1-s)\sum_{i,j}V_{ij}(q+s\eta,t)\eta^i\eta^j\,\dd s$ と仮定から、括弧内は $\tfrac12K\abs\eta^2$ 以下である。よって

$$
S[q+\eta]-S[q]\ge\frac12\int_{t_0}^{t_1}\big(m\abs{\dot\eta}^2-K\abs{\eta}^2\big)\dd t .
$$

$K\le0$ なら右辺は $\tfrac12m\int\abs{\dot\eta}^2\,\dd t$ 以上で、$\eta\ne0$ なら $\dot\eta\ne0$（端点で零なので）より正である。
$K>0$ なら、:ref[lem-am-wirtinger] を区間 $[t_0,t_1]$ に平行移動して各成分に適用し、$\int\abs{\dot\eta}^2\ge(\pi/\tau)^2\int\abs\eta^2$ を得る。右辺は $\tfrac12(m\pi^2/\tau^2-K)\int\abs\eta^2\,\dd t>0$ である。
:::

:::example[調和振動子の作用]{#exm-am-oscillator-action}
$L=\tfrac12m\dot x^2-\tfrac12m\omega^2x^2$ では $K=m\omega^2$ であり、:ref[prop-am-short-time-minimum] は $\tau<\pi/\omega$、すなわち半周期より短い区間で実軌道が作用の唯一の最小を与えることを示す。
この評価は最良である。差は二次式なので展開に余りがなく、

$$
S[x+\eta]-S[x]=\frac m2\int_{t_0}^{t_1}\big(\dot\eta^2-\omega^2\eta^2\big)\,\dd t
$$

が厳密に成り立つ。$\tau>\pi/\omega$ なら $\eta=\sin(\pi(t-t_0)/\tau)$ で右辺は $\tfrac{m\tau}{4}(\pi^2/\tau^2-\omega^2)<0$ となり、実軌道は作用の極小を与えない（鞍点である）。
$\tau=\pi/\omega$ では、$x=0$ から出発する解はすべて半周期後に $x=0$ に戻るので、同じ端点を結ぶ実軌道が無数にある。これが共役点の力学的な意味である。
:::

したがって Hamilton の原理を「最小作用の原理」と呼ぶのは、十分短い時間区間についてだけ正しい。全区間については停留としか言えない。

### 中心力を極座標で表す

平面上を動く質量 $m$ の質点について、位置エネルギーが原点からの距離 $r$ だけの関数 $V(r)$ であるとする。
$r>0$ の領域で極座標 $(r,\varphi)$ を使うと、位置は $(r\cos\varphi,r\sin\varphi)$ であり、
半径方向と円周方向の速度成分は $\dot r$ と $r\dot\varphi$ である。作用の被積分関数と二つの方程式は

$$
\begin{aligned}
L&=\tfrac12 m(\dot r^2+r^2\dot\varphi^2)-V(r),\\
m\ddot r&=mr\dot\varphi^2-V'(r),\qquad
\dv{}{t}(mr^2\dot\varphi)=0.
\end{aligned}
\label{eq-am-polar-central}
$$

第一式の運動エネルギーを $r$ で偏微分したときの $mr\dot\varphi^2$ が、半径方向に必要な向心加速度を与える。
実際、動径成分をまとめると $m(\ddot r-r\dot\varphi^2)=-V'(r)$ である。
また $\partial L/\partial\varphi=0$ なので、$p_\varphi=mr^2\dot\varphi$ は一定である。直交座標に戻せば $p_\varphi=m(x\dot y-y\dot x)$ であり、原点まわりの角運動量に等しい。
この系の解析は [中心力問題と Kepler 問題](/physics/analytical-mechanics/central-force/) で行う。

## Lagrangian の任意性

### 全時間微分を加える

Lagrangian を一つ見つけても、その式だけが同じ運動を記述するわけではない。
まず、零でない定数 $c$ を掛けた $cL$ は同じ方程式を与える。これは作用の単位を変えることに当たる。
次に、$C^2$ 級の一価関数 $F(q,t)$ から、次の $L'$ を作る。

$$
L'=L+\dv{F}{t}
=L+\sum_{i=1}^n\pdv{F}{q^i}\dot q^i+\pdv{F}{t}.
\label{eq-am-total-derivative}
$$

このとき作用の差は $F(q(t_1),t_1)-F(q(t_0),t_0)$ だけであり、固定端ではその変分が零になる。
したがって :ref[prin-am-hamilton] は同じ軌道を定める。方程式の左辺でも差が消えることを確かめられる。実際、

$$
\pdv{}{\dot q^i}\left(\dv{F}{t}\right)=\pdv{F}{q^i},
\qquad
\pdv{}{q^i}\left(\dv{F}{t}\right)=\dv{}{t}\pdv{F}{q^i}
$$

である。後者は混合偏微分を交換すると従う。よって $E_i(\dd F/\dd t)=0$ が恒等的に成り立つ。
一方、一般化運動量の値は $p'_i=p_i+\partial F/\partial q^i$ に移る。どの $L$ を使ったかで運動量の表示が変わる例を、後の電磁場中の荷電粒子の節で扱う。

逆に、Euler–Lagrange 式を恒等的に零にする Lagrangian は全時間微分に限られる。

:::proposition[恒等的に零の Euler–Lagrange 式]{#prop-am-null-lagrangian}
$D\subset\R^n\times\R$ を原点に関して星形の開集合（$(q,t)\in D$ なら $0\le s\le1$ で $(sq,st)\in D$）とし、$\Lambda(q,\dot q,t)$ を $D\times\R^n$ 上の $C^2$ 級関数とする。
すべての $(q,\dot q,\ddot q,t)$ で $E_i(\Lambda)=0$ $(i=1,\ldots,n)$ が成り立つならば、$C^2$ 級関数 $F(q,t)$ があって $\Lambda=\dd F/\dd t$ である。
:::

:::proof
$E_i(\Lambda)$ を展開すると

$$
E_i(\Lambda)=\sum_j\frac{\partial^2\Lambda}{\partial\dot q^i\partial\dot q^j}\ddot q^j+\sum_j\frac{\partial^2\Lambda}{\partial\dot q^i\partial q^j}\dot q^j+\frac{\partial^2\Lambda}{\partial\dot q^i\partial t}-\pdv{\Lambda}{q^i}
$$

である。$\ddot q$ は任意なので、その係数 $\partial^2\Lambda/\partial\dot q^i\partial\dot q^j$ は恒等的に零である。
よって $\partial\Lambda/\partial\dot q^i$ は $\dot q$ によらず、$\Lambda=\sum_jA_j(q,t)\dot q^j+B(q,t)$ の形をもつ（$A_j$ と $B=\Lambda(q,0,t)$ は $C^1$ 級）。これを代入すると

$$
E_i(\Lambda)=\sum_j\left(\pdv{A_i}{q^j}-\pdv{A_j}{q^i}\right)\dot q^j+\pdv{A_i}{t}-\pdv{B}{q^i}
$$

で、これが任意の $\dot q$ で零なので

$$
\pdv{A_i}{q^j}=\pdv{A_j}{q^i},\qquad \pdv{A_i}{t}=\pdv{B}{q^i}
$$

である。そこで

$$
F(q,t)=\int_0^1\Big(\sum_iA_i(sq,st)\,q^i+B(sq,st)\,t\Big)\dd s
$$

とおく。積分記号の下で微分し、上の関係を使うと

$$
\pdv{F}{q^k}=\int_0^1\Big(A_k+s\sum_i\pdv{A_k}{q^i}q^i+s\pdv{A_k}{t}t\Big)(sq,st)\,\dd s=\int_0^1\dv{}{s}\big(sA_k(sq,st)\big)\dd s=A_k(q,t)
$$

となる。同様に $\partial F/\partial t=\int_0^1\dd(sB(sq,st))/\dd s\,\dd s=B$ である。したがって $\dd F/\dd t=\sum_kA_k\dot q^k+B=\Lambda$ を得る。
:::

この命題から、二つの Lagrangian が恒等的に同じ Euler–Lagrange 式を与えるのは、（星形の領域で）全時間微分だけ異なる場合に限ることがわかる。
ただし、Euler–Lagrange 式が**異なる式**でも**同じ解**をもつことはある。例えば一次元の自由粒子では、$f''\ne0$ を満たす任意の $f$ について $L=f(\dot x)$ の方程式は $f''(\dot x)\ddot x=0$ であり、解は等速運動である。
したがって「同じ運動を与える Lagrangian」は全時間微分の差よりも広い。物理的な Lagrangian を選ぶには、次の小節のような追加の要請が必要である。

### Galilei 不変性と自由粒子

慣性系 $K$ に対して一定速度 $\vb{V}$ で動く座標系 $K'$ を考え、$\vb{r}=\vb{r}'+\vb{V}t$ とする。これを :term[Galilei 変換]{en="Galilean transformation" yomi="がりれいへんかん"} という。
Newton 方程式は Galilei 変換で形を変えない（Galilei の相対性原理）。これを Lagrangian で表すと次のようになる。
自由粒子の $L=\tfrac12m\vb{v}^2$ に $\vb{v}=\vb{v}'+\vb{V}$ を代入すると

$$
\tfrac12m\vb{v}^2=\tfrac12m\vb{v}'^2+\dv{}{t}\Big(m\vb{V}\cdot\vb{r}'+\tfrac12m\vb{V}^2t\Big)
$$

であり、$K'$ での Lagrangian は同じ形 $\tfrac12m\vb{v}'^2$ と全時間微分の差しかない。

逆に、この性質から自由粒子の Lagrangian の形が決まる。空間と時間の一様性から自由粒子の $L$ は $\vb{r}$ と $t$ によらず、空間の等方性から速度の向きにもよらないとすれば、$L=f(\abs{\vb{v}}^2)$ と書ける :cite[landau-mechanics]。

:::proposition[自由粒子の Lagrangian]{#prop-am-free-particle-lagrangian}
$f$ を $[0,\infty)$ 上の $C^2$ 級関数とし、$L(\vb{v})=f(\abs{\vb{v}}^2)$ とする。すべての $\vb{V}\in\R^3$ について、$\Lambda_{\vb V}(\vb{v})=L(\vb{v}+\vb{V})-L(\vb{v})$ の Euler–Lagrange 式が恒等的に零であるとする。このとき定数 $c,d$ があって $f(s)=cs+d$ である。
:::

:::proof
$\Lambda_{\vb V}$ は位置と時刻によらないので、$E_i(\Lambda_{\vb V})=\sum_j(\partial^2\Lambda_{\vb V}/\partial v_i\partial v_j)\dot v_j$ である。これが任意の $\dot{\vb v}$ で零なので、$L$ の Hesse 行列 $H(\vb{v})$ は $H(\vb{v}+\vb{V})=H(\vb{v})$ を満たし、定数行列である。
$g(\vb v)=f(\abs{\vb v}^2)$ の Hesse 行列は $H_{ij}=2f'(\abs{\vb v}^2)\delta_{ij}+4f''(\abs{\vb v}^2)v_iv_j$ である。$\vb v\ne0$ と、それに垂直な単位ベクトル $\vb e$ をとれば $\sum_{ij}H_{ij}e_ie_j=2f'(\abs{\vb v}^2)$ であり、これが $\vb v=0$ での値 $2f'(0)$ に等しい。
よって $f'$ は定数であり、$f(s)=cs+d$ となる。
:::

定数 $d$ は運動に影響しない。$c$ を $m/2$ と書いたものが粒子の質量の定義を与える。
$m>0$ でなければならないのは、$m<0$ なら速く往復する経路で作用がいくらでも小さくなり、:ref[prop-am-short-time-minimum] のような最小性が失われるからである。
相互作用する質点系では $L=\sum_a\tfrac12m_a\vb{v}_a^2-U(\vb{r}_1,\ldots,\vb{r}_N)$ をとり、$U$ が相互作用を表す。$U$ が相対位置だけに依存するときに成り立つ保存則は [対称性と保存則](/physics/analytical-mechanics/symmetries/) で扱う。

### 点変換に対する共変性

同じ配位を別の一般化座標 $u^1,\ldots,u^n$ で指定し、$q^i=q^i(u,t)$ と書く。
これを :term[点変換]{en="point transformation" yomi="てんへんかん"} という。以下では時刻は変えず、変換は $C^2$ 級で、Jacobi 行列 $(\partial q^i/\partial u^j)$ は正則とする。
新しい Lagrangian は、速度も変換して

$$
\widetilde L(u,\dot u,t)=L(q(u,t),\dot q(u,\dot u,t),t),\qquad
\dot q^i=\sum_{j=1}^n\pdv{q^i}{u^j}\dot u^j+\pdv{q^i}{t}
$$

と定める。時間に依存する変換では、速度の最後の項を残すことが大切である。

:::proposition[Lagrange 方程式の共変性]{#prop-am-point-covariance}
運動方程式の左辺を $E_i(L)=\dd(\partial L/\partial\dot q^i)/\dd t-\partial L/\partial q^i$ と書くと、

$$
E_j(\widetilde L)=\sum_{i=1}^n E_i(L)\pdv{q^i}{u^j}.
\label{eq-am-point-covariance}
$$

一般化力の右辺も同じ変換に従うので、新旧の座標で書いた方程式は同値である。
:::

:::proof
:ref[lem-am-velocity-identities] の位置を $q^i(u,t)$ に置き換えると、
$\partial\dot q^i/\partial\dot u^j=\partial q^i/\partial u^j$、
$\partial\dot q^i/\partial u^j=\dd(\partial q^i/\partial u^j)/\dd t$ である。ゆえに

$$
\begin{aligned}
\pdv{\widetilde L}{\dot u^j}
&=\sum_{i=1}^np_i\pdv{q^i}{u^j},\\
\pdv{\widetilde L}{u^j}
&=\sum_{i=1}^n\left(
\pdv{L}{q^i}\pdv{q^i}{u^j}
+p_i\dv{}{t}\pdv{q^i}{u^j}\right).
\end{aligned}
$$

上段を全時間微分して下段を引くと、$p_i$ と Jacobi 行列の微分との積が相殺して主張を得る。
また、仮想仕事に $\delta q^i=\sum_{j=1}^n(\partial q^i/\partial u^j)\delta u^j$ を代入すると
$\widetilde Q_j^{\mathrm{nc}}=\sum_{i=1}^nQ_i^{\mathrm{nc}}(\partial q^i/\partial u^j)$ となる。
行列が正則なので逆方向にも方程式を移せる。
:::

例えば自由粒子の直交座標 $x$ から、原点が $b(t)$ にある座標 $u=x-b(t)$ に移ると
$\widetilde L=m(\dot u+\dot b)^2/2$ である。その方程式は $m(\ddot u+\ddot b)=0$ となる。
移動原点から見た加速度 $-\ddot b$ は、速度の変換を含めるだけで現れる。
$\widetilde L=\tfrac12m\dot u^2+m\dot u\dot b+\tfrac12m\dot b^2$ の第二項は $\dd(mu\dot b)/\dd t-mu\ddot b$ と書けるので、全時間微分を除けば $\widetilde L\simeq\tfrac12m\dot u^2-m\ddot b\,u$ であり、慣性力 $-m\ddot b$ のポテンシャル $m\ddot b\,u$ が現れる。

:::example[回転座標系と慣性力]{#exm-am-rotating-frame}
慣性系の座標 $\vb{x}$ と、$z$ 軸のまわりに一定の角速度 $\Omega$ で回る座標系の座標 $\vb{x}'$ を $\vb{x}=\mathsf R(t)\vb{x}'$ で結ぶ。$\mathsf R(t)$ は角 $\Omega t$ の回転行列である。
$\boldsymbol\Omega=(0,0,\Omega)$ とすると、$\dot{\vb{x}}=\mathsf R(t)\big(\dot{\vb{x}}'+\boldsymbol\Omega\times\vb{x}'\big)$ であり、回転は長さを保つので

$$
\widetilde L=\tfrac12m\abs{\dot{\vb{x}}'}^2+m\dot{\vb{x}}'\cdot(\boldsymbol\Omega\times\vb{x}')+\tfrac12m\abs{\boldsymbol\Omega\times\vb{x}'}^2-V(\vb{x}')
$$

となる。ここで $V$ は $\vb{x}'$ で表した位置エネルギーである。
$\dot{\vb{x}}'\cdot(\boldsymbol\Omega\times\vb{x}')=\vb{x}'\cdot(\dot{\vb{x}}'\times\boldsymbol\Omega)$ の $\vb{x}'$ 勾配は $\dot{\vb{x}}'\times\boldsymbol\Omega$、
$\tfrac12\abs{\boldsymbol\Omega\times\vb{x}'}^2=\tfrac12(\Omega^2\abs{\vb{x}'}^2-(\boldsymbol\Omega\cdot\vb{x}')^2)$ の勾配は $-\boldsymbol\Omega\times(\boldsymbol\Omega\times\vb{x}')$ である。
一般化運動量は $m\dot{\vb{x}}'+m\boldsymbol\Omega\times\vb{x}'$ なので、Lagrange 方程式は

$$
m\ddot{\vb{x}}'=-\grad'V-2m\boldsymbol\Omega\times\dot{\vb{x}}'-m\boldsymbol\Omega\times(\boldsymbol\Omega\times\vb{x}')
$$

となる。右辺第二項が Coriolis 力、第三項が遠心力である。
:eqref[eq-am-kinetic-generalized] でいえば、Coriolis 力は一般化速度の一次の項 $b_i\dot q^i$ から、遠心力は零次の項 $T_0$ から生じる。
:::

:::example[鉛直に振動する支点をもつ振り子]{#exm-am-vertical-support}
支点の高さが既知の関数 $h(t)$ に従って上下する振り子を考える。位置は $(\ell\sin\theta,\ h(t)-\ell\cos\theta)$ なので

$$
L=\tfrac12m\big(\ell^2\dot\theta^2+2\ell\dot h\dot\theta\sin\theta+\dot h^2\big)-mg(h-\ell\cos\theta).
$$

$m\ell\dot h\dot\theta\sin\theta=\dv{}{t}(-m\ell\dot h\cos\theta)+m\ell\ddot h\cos\theta$ であり、$\tfrac12m\dot h^2-mgh$ は時刻だけの関数なのでその原始関数の全時間微分である。これらの全時間微分を除くと

$$
L'=\tfrac12m\ell^2\dot\theta^2+m(g+\ddot h)\ell\cos\theta,\qquad
\ddot\theta+\frac{g+\ddot h(t)}{\ell}\sin\theta=0
$$

を得る。支点の加速度は重力加速度を $g+\ddot h$ に置き換える。例えば支点が自由落下する $\ddot h=-g$ では、振り子は重さを感じない。
:::

## 電磁場中の荷電粒子

### 一般化ポテンシャル

これまでの Lagrangian $T-V$ では、力は位置だけのポテンシャルから導かれた。速度に依存する力でも同じ形に取り込める場合がある。

:::definition[一般化ポテンシャル]{#def-am-generalized-potential}
関数 $U(q,\dot q,t)$ があって、一般化力が

$$
Q_i=\dv{}{t}\pdv{U}{\dot q^i}-\pdv{U}{q^i},\qquad i=1,\ldots,n
$$

と書けるとき、$U$ を :term[一般化ポテンシャル]{en="generalized potential" yomi="いっぱんかぽてんしゃる"}（速度に依存するポテンシャル）という。
:::

このとき慣性項の式 $E_i(T)=Q_i$ は $E_i(T-U)=0$ と同値なので、$L=T-U$ が運動方程式を与える。
右辺の時間微分は一般に $\ddot q$ を含む。加速度によらない力を与えるには、$\ddot q$ の係数 $\partial^2U/\partial\dot q^i\partial\dot q^j$ が零、すなわち $U$ が速度について一次式 $U=\phi(q,t)-\sum_iA_i(q,t)\dot q^i$ でなければならない。電磁場はちょうどこの形である。

### ポテンシャルと一般化運動量

非相対論的な質量 $m$、電荷 $q$ の粒子が、外から与えた滑らかな電磁場の中を動くとする。
この節の $q$ は一定の電荷であり、位置座標には $x_1,x_2,x_3$ を用いる。粒子自身が作る場と放射は扱わない。
SI 系で電場 $\vb{E}$ と磁束密度 $\vb{B}$ を、スカラーポテンシャル $\phi$ とベクトルポテンシャル $\vb{A}$ から

$$
\vb{E}=-\grad\phi-\pdv{\vb{A}}{t},\qquad
\vb{B}=\nabla\times\vb{A}
$$

と表す。ここでの時間偏微分は位置を固定して行う（[電磁誘導](/physics/electromagnetism/induction/)）。

:::example[荷電粒子の Lagrangian]{#exm-am-electromagnetic-lagrangian}
$v_i=\dot x_i$ として、一般化ポテンシャル $U=q\phi-q\vb{v}\cdot\vb{A}$ を用いた Lagrangian は

$$
L=\tfrac12m\vb{v}^2-q\phi+q\vb{v}\cdot\vb{A}.
\label{eq-am-em-lagrangian}
$$

運動エネルギー以外の二項はどちらもエネルギーの単位をもつ。
実際、$\phi$ の単位は $\mathrm{J/C}$、$\vb{A}$ の単位は $\mathrm{T\,m}=\mathrm{kg\,m/(C\,s)}$ である。
この式に速度に比例する項があるため、一般化運動量は

$$
p_i=mv_i+qA_i
$$

となる。これを :term[正準運動量]{en="canonical momentum" yomi="せいじゅんうんどうりょう"} と呼び、質点の運動から測る :term[力学的運動量]{en="kinetic momentum" yomi="りきがくてきうんどうりょう"} $mv_i$ と区別する。
:::

### Lorentz 力を成分から組み立てる

位置微分では速度を固定し、運動量の時間微分では粒子の位置も時間とともに動かす。この違いに注意すると

$$
\begin{aligned}
\pdv{L}{x_i}
&=-q\pdv{\phi}{x_i}+q\sum_{j=1}^3v_j\pdv{A_j}{x_i},\\
\dv{p_i}{t}
&=m\dot v_i+q\pdv{A_i}{t}
+q\sum_{j=1}^3v_j\pdv{A_i}{x_j}.
\end{aligned}
$$

両者を等置して、$\phi$ と時間偏微分の項を $E_i$ にまとめると

$$
m\dot v_i=qE_i+q\sum_{j=1}^3 v_j
\left(\pdv{A_j}{x_i}-\pdv{A_i}{x_j}\right).
$$

例えば $i=1$ の和では $j=1$ の項が消える。残りは
$v_2(\partial A_2/\partial x_1-\partial A_1/\partial x_2)
+v_3(\partial A_3/\partial x_1-\partial A_1/\partial x_3)=v_2B_3-v_3B_2$ である。
添字を巡回させれば他の二成分も得られ、運動方程式は

$$
m\dot{\vb{v}}=q(\vb{E}+\vb{v}\times\vb{B})
\label{eq-am-lorentz-force}
$$

となる。右辺を :term[Lorentz 力]{en="Lorentz force" yomi="ろーれんつりょく"} という。
両辺と速度の内積をとると $\dot T=q\vb{v}\cdot\vb{E}$ であり、磁場からの力は速度の向きだけを変える。
Lorentz 力は速度に依存し、位置だけのポテンシャルからは導けないが、一般化ポテンシャルにより Lagrange 形式に含まれる。
相対論的な粒子では運動エネルギーの項を $-mc^2\sqrt{1-\vb{v}^2/c^2}$ に置き換える（[特殊相対論](/physics/special-relativity/)）。

### ゲージ変換と境界項

任意の滑らかな関数 $\chi(\vb{r},t)$ に対し、
$\vb{A}'=\vb{A}+\grad\chi$、$\phi'=\phi-\partial\chi/\partial t$ と変えることを
:term[ゲージ変換]{en="gauge transformation" yomi="げーじへんかん"} という。
回転 $\nabla\times\grad\chi$ は零であり、$\vb{E}$ の二つの項に現れる $\grad(\partial\chi/\partial t)$ は相殺するので、電磁場は同じである。
さらに :eqref[eq-am-em-lagrangian] の差を計算すれば

$$
L'-L=q\pdv{\chi}{t}+q\vb{v}\cdot\grad\chi
=\dv{}{t}(q\chi).
$$

したがって :eqref[eq-am-total-derivative] と同じ仕組みで運動が保たれる。
実際 $p'_i=p_i+q\,\partial\chi/\partial x_i$ だが、$p'_i-qA'_i=p_i-qA_i=mv_i$ である。
正準運動量はゲージに依存し、力学的運動量はゲージに依存しない。量子力学ではこの区別が本質的になる。

:::example[一様な磁場と直交する電場]{#exm-am-crossed-fields}
$\vb{B}=(0,0,B)$、$\vb{E}=(E,0,0)$（$B>0$、$E$ は定数）を $\vb{A}=(0,Bx,0)$、$\phi=-Ex$ で表す。

$$
L=\tfrac12m(\dot x^2+\dot y^2+\dot z^2)+qEx+qBx\dot y .
$$

$L$ は $y,z$ を含まないので $p_y=m\dot y+qBx$ と $p_z=m\dot z$ が一定であり、$x$ の方程式は $m\ddot x=qE+qB\dot y$ である。
$\omega=qB/m$、$u=\dot y+E/B$ とおくと、$\ddot y=-\omega\dot x$ と合わせて

$$
\ddot x=\omega u,\qquad \dot u=-\omega\dot x
$$

となる。$(\dot x,u)$ は角振動数 $\omega$ で回転し、$\dot x=w\sin(\omega t+\gamma)$、$\dot y=-E/B+w\cos(\omega t+\gamma)$（$w,\gamma$ は定数）である。
$E=0$ では、粒子は :term[サイクロトロン角振動数]{en="cyclotron frequency" yomi="さいくろとろんかくしんどうすう"} $\omega=qB/m$ で円運動する。
$E\ne0$ では円運動に一定速度 $-E/B$ の $y$ 方向の移動が重なる。これは $\vb{E}\times\vb{B}/B^2=(0,-E/B,0)$ に等しく、電荷の符号や質量によらない（$E\times B$ ドリフト）。
:::

## 散逸と Rayleigh 散逸関数

### 速度に比例する抵抗

抵抗力がある系では、座標と速度だけの通常の $T-V$ に加え、右辺の一般化力を評価する。
例えば静止した媒質中で直線上を動く粒子に $-c\dot x$（$c>0$）が働くとき、速度による偏微分でこの力を得るには
$\mathcal R=c\dot x^2/2$ を用意し、$-\partial\mathcal R/\partial\dot x$ をとればよい。

:::definition[Rayleigh 散逸関数]{#def-am-rayleigh}
粘性抵抗の一般化力が $Q_i^{\mathrm{diss}}=-\sum_{j=1}^n c_{ij}(q,t)\dot q^j$ と表され、
行列 $(c_{ij})$ が対称かつ半正定値であるとする。このとき

$$
\mathcal R(q,\dot q,t)=\tfrac12\sum_{i=1}^n\sum_{j=1}^n c_{ij}(q,t)\dot q^i\dot q^j
$$

を :term[Rayleigh 散逸関数]{en="Rayleigh dissipation function" yomi="れいりーさんいつかんすう"} という。単位は $\mathrm W$ である。
:::

$\dot q^i$ で微分すると二次式の二つの因子から等しい寄与が出るので、
$\partial\mathcal R/\partial\dot q^i=\sum_{j=1}^n c_{ij}\dot q^j$ となる。
他に駆動力 $Q_i^{\mathrm{drive}}$ があれば、方程式は

$$
\dv{}{t}\pdv{L}{\dot q^i}-\pdv{L}{q^i}+\pdv{\mathcal R}{\dot q^i}
=Q_i^{\mathrm{drive}} .
\label{eq-am-rayleigh}
$$

各質点に $\vb{f}_a=-k_a\vb{v}_a$（$k_a\ge0$）の抵抗が働く場合には、直交座標で $\mathcal R=\tfrac12\sum_ak_a\vb{v}_a^2$ とし、一般化座標の速度を代入すればよい。
実際、:ref[lem-am-velocity-identities] から $-\partial\mathcal R/\partial\dot q^i=-\sum_ak_a\vb{v}_a\cdot\partial\vb{r}_a/\partial q^i$ であり、これは :eqref[eq-am-generalized-force] で定めた抵抗力の一般化力に等しい。

静止した座標系でばね定数 $k>0$ のばねに質量 $m$ を取り付け、駆動力を $f(t)$ とする。
$L=m\dot x^2/2-kx^2/2$ と $\mathcal R=c\dot x^2/2$ を代入すると
$m\ddot x+c\dot x+kx=f(t)$ となる。その解の性質は [振動と波動](/physics/vibrations-waves/) で詳しく扱う。

:::example[電気回路との対応]{#exm-am-lcr}
自己インダクタンス $L_{\mathrm{ind}}$、抵抗 $r$、静電容量 $C$ を直列につなぎ、起電力 $\mathcal E(t)$ で駆動する。コンデンサーの電荷 $Q_{\mathrm c}$ を一般化座標とし、

$$
L=\tfrac12L_{\mathrm{ind}}\dot Q_{\mathrm c}^2-\frac{Q_{\mathrm c}^2}{2C},\qquad
\mathcal R=\tfrac12r\dot Q_{\mathrm c}^2,\qquad Q^{\mathrm{drive}}=\mathcal E(t)
$$

とおくと、:eqref[eq-am-rayleigh] は $L_{\mathrm{ind}}\ddot Q_{\mathrm c}+r\dot Q_{\mathrm c}+Q_{\mathrm c}/C=\mathcal E(t)$ となり、Kirchhoff の電圧則に一致する。
インダクタンスの磁場エネルギーが運動エネルギー、コンデンサーの静電エネルギーが位置エネルギー、Joule 熱の発生率の半分が $\mathcal R$ に対応する。
:::

散逸関数の考え方は線形抵抗以外にも使える。速度の $C^1$ 級関数 $\mathcal R$ が速度について $k$ 次の同次関数なら、Euler の同次関数定理 $\sum_i\dot q^i\partial\mathcal R/\partial\dot q^i=k\mathcal R$ により、抵抗力の仕事率は $-k\mathcal R$ である。
例えば速さの二乗に比例する空気抵抗 $-c\abs{\vb{v}}\vb{v}$ は $\mathcal R=\tfrac c3\abs{\vb v}^3$ から得られ、仕事率は $-3\mathcal R$ である。

### 仕事率と力学的エネルギー

抵抗力の一般化力による仕事率は
$\sum_{i=1}^nQ_i^{\mathrm{diss}}\dot q^i=-\sum_{i=1}^n\sum_{j=1}^n c_{ij}\dot q^i\dot q^j=-2\mathcal R\le0$ である。
係数 $2$ は $\mathcal R$ が速度の二次式であることから生じる。

移動する拘束を使う場合には、一般化力と一般化速度の積を、実験室での仕事率と区別する必要がある。
実際、任意の力 $\vb{f}_a$ とそれに対応する一般化力 $Q_i$ について、:eqref[eq-am-coordinate-velocity] より

$$
\sum_{a=1}^N\vb{f}_a\cdot\vb{v}_a
=\sum_{i=1}^n Q_i\dot q^i+\sum_{a=1}^N\vb{f}_a\cdot\pdv{\vb{r}_a}{t}.
$$

実験室で測る仕事率には、座標値を固定したままでも生じる速度の寄与が加わる。

より一般に $E_L=\sum_{i=1}^n p_i\dot q^i-L$ とおいて軌道に沿って微分すると、

$$
\begin{aligned}
\dv{E_L}{t}
&=\sum_{i=1}^n\left(\dot p_i-\pdv{L}{q^i}\right)\dot q^i-\pdv{L}{t}\\
&=\sum_{i=1}^n Q_i^{\mathrm{drive}}\dot q^i-2\mathcal R-\pdv{L}{t}.
\end{aligned}
$$

第一行では $p_i\ddot q^i$ が $\dot L$ の項と相殺した。
座標の位置表示が時刻に陽に依存せず、$T$ が一般化速度の同次二次式で $L=T-V(q)$ なら、
$\sum_{i=1}^n(\partial T/\partial\dot q^i)\dot q^i=2T$ なので $E_L=T+V$ である。
この条件で駆動力がなければ、力学的エネルギーの減少率は $2\mathcal R$ である。
一方、:ref[exm-am-rotating-hoop] のように外部から拘束を動かす場合には、$E_L$ と $T+V$ は一般には一致しない。
$E_L$ の保存の条件と意味は [対称性と保存則](/physics/analytical-mechanics/symmetries/) で詳しく調べる。散逸力を一般化力として扱う方法は :cite[goldstein-cm] も参照。

:::remark[時間に依存する Lagrangian による減衰振動]
散逸関数を使わずに減衰振動を Lagrangian だけで表すこともできる。$\gamma>0$ として

$$
L=e^{\gamma t}\left(\tfrac12m\dot x^2-\tfrac12kx^2\right)
$$

とおくと、$E(L)=e^{\gamma t}(m\ddot x+\gamma m\dot x+kx)$ なので、方程式は $m\ddot x+\gamma m\dot x+kx=0$ となる。
ただしこの $L$ は $T-V$ の形ではなく、一般化運動量 $e^{\gamma t}m\dot x$ も力学的運動量と一致しない。散逸を Lagrangian に含める方法は一意ではなく、物理的な解釈には注意を要する。
:::

## 演習問題

::::exercise{level=basic}
Atwood の装置で、$m_1$ の下向き変位を $x$ とする。糸の張力を求め、$m_1=m_2$ の場合を確かめよ。
:::solution
理想的な糸と滑車では両側の張力の大きさは等しい。これを $H$ とすると第 1 質点の Newton 方程式は
$m_1\ddot x=m_1g-H$ である。上で求めた $\ddot x$ を代入すれば

$$
H=m_1g\left(1-\frac{m_1-m_2}{m_1+m_2}\right)=\frac{2m_1m_2}{m_1+m_2}g .
$$

第 2 質点の下向き加速度は $-\ddot x$ なので、$-m_2\ddot x=m_2g-H$ からも同じ値を得る。等質量 $m$ なら $H=mg$ である。
:::
::::

::::exercise{level=basic}
自由粒子の $L=m\dot x^2/2$ に、$F(x,t)=k(t)x$ の全時間微分を加える。
新しい運動量と運動方程式を求めよ。
:::solution
$L'=m\dot x^2/2+k\dot x+\dot kx$ である。よって $p'=m\dot x+k$、$\partial L'/\partial x=\dot k$ である。
運動量を時間微分すると $m\ddot x+\dot k$ となるから、方程式は $m\ddot x+\dot k-\dot k=0$ である。
作用は端点での $kx$ の差だけ変わり、軌道は自由粒子の直線運動である。
:::
::::

::::exercise{level=basic}
一様な磁場 $\vb{B}=(0,0,B)$ を $\vb{A}=(0,Bx,0)$、$\phi=0$ で表す。
$y$ に共役な運動量を求め、$x,y$ に関する運動方程式と整合することを示せ。
:::solution
$L=m(\dot x^2+\dot y^2+\dot z^2)/2+qBx\dot y$ は $y$ に依存しない。
したがって $p_y=m\dot y+qBx$ は一定であり、微分すると $m\ddot y=-qB\dot x$ となる。
一方、$p_x=m\dot x$ と $\partial L/\partial x=qB\dot y$ により $m\ddot x=qB\dot y$ である。
これは $q\vb{v}\times\vb{B}=q(B\dot y,-B\dot x,0)$ の各成分と一致する。
:::
::::

::::exercise{level=basic}
:ref[exm-am-pendulum-tension] の振り子を角 $\theta_0$（$0<\theta_0<\pi/2$）から静かに放す。張力を $\theta$ の関数として求め、その最大値を求めよ。
:::solution
エネルギー保存 $\tfrac12m\ell^2\dot\theta^2-mg\ell\cos\theta=-mg\ell\cos\theta_0$ から $\ell\dot\theta^2=2g(\cos\theta-\cos\theta_0)$ である。
:ref[exm-am-pendulum-tension] の結果に代入すると

$$
S=m\ell\dot\theta^2+mg\cos\theta=mg(3\cos\theta-2\cos\theta_0).
$$

$\abs\theta\le\theta_0<\pi/2$ の範囲で $\cos\theta$ は $\theta=0$ で最大なので、最大値は最下点での $mg(3-2\cos\theta_0)$ である。$\theta_0\to0$ で $mg$ に近づく。
:::
::::

::::exercise{level=basic}
:ref[exm-am-vertical-support] で支点が $h(t)=h_0\cos\gamma t$ と振動するとき、運動方程式を書け。$\theta$ が小さい範囲では、これがどのような線形方程式になるか。
:::solution
$\ddot h=-h_0\gamma^2\cos\gamma t$ を代入して

$$
\ddot\theta+\frac{g-h_0\gamma^2\cos\gamma t}{\ell}\sin\theta=0
$$

である。$\sin\theta\simeq\theta$ とすると $\ddot\theta+\big(g/\ell-(h_0\gamma^2/\ell)\cos\gamma t\big)\theta=0$ となり、係数が周期的に変わる線形方程式（Mathieu 方程式）になる。
係数の変調によって振幅が成長するパラメータ共振や、速い振動による倒立振り子の安定化はこの方程式から説明される。
:::
::::

::::exercise{level=advanced}
:ref[exm-am-rotating-hoop] で $\theta=0$ の近くと、$\Omega^2>g/a$ のときの最下点以外の下半円の相対平衡 $\theta_*$ の近くで運動を線形化せよ。安定な場合の微小振動の角振動数を求めよ。
:::solution
:eqref[eq-am-rotating-hoop] の右辺を $f(\theta)$ とおけば、$\theta=\theta_*+\eta$ の一次では $\ddot\eta=f'(\theta_*)\eta$ である。$f'(\theta)=-\Omega^2\sin^2\theta+(\Omega^2\cos\theta-g/a)\cos\theta$ である。
最下点では $f'(0)=\Omega^2-g/a$ であり、$\Omega^2<g/a$ なら安定で角振動数は $\sqrt{g/a-\Omega^2}$、$\Omega^2>g/a$ なら不安定である。
等号の場合は線形項が消える。$\Omega^2=g/a$ では $f(\theta)=(g/a)(\cos\theta-1)\sin\theta\simeq-(g/2a)\theta^3$ で、最下点からの小さな変位に対する力は復元方向だが、線形の角振動数は零である。
一方、$\cos\theta_*=g/(a\Omega^2)$ を満たす二点では第二項が消え

$$
f'(\theta_*)=-\Omega^2\sin^2\theta_*<0.
$$

したがって二点とも安定で、角振動数は $\abs{\Omega\sin\theta_*}=\Omega\sqrt{1-g^2/(a^2\Omega^4)}$ となる。
:::
::::

::::exercise{level=advanced}
中心力の問題で角運動量 $p_\varphi=h\ne0$ を固定する。半径 $r_0>0$ の円軌道の条件を求め、$V(r)=-\kappa/r$（$\kappa>0$）で半径をわずかに変えたときの角振動数を計算せよ。
:::solution
$\dot\varphi=h/(mr^2)$ を :eqref[eq-am-polar-central] の動径方程式に代入すると
$m\ddot r=h^2/(mr^3)-V'(r)$ である。半径一定の条件は $V'(r_0)=h^2/(mr_0^3)$ となる。
引力ポテンシャルでは $V'(r)=\kappa/r^2$ なので $h^2=m\kappa r_0$ を得る。右辺 $h^2/(mr^3)-\kappa/r^2$ の $r$ 微分は $-3h^2/(mr^4)+2\kappa/r^3$ で、$r_0$ では $-\kappa/r_0^3$ となる。$r=r_0+\eta$ の一次では

$$
m\ddot\eta=-\frac{\kappa}{r_0^3}\eta,\qquad
\omega=\sqrt{\frac{\kappa}{mr_0^3}} .
$$

$\omega^2>0$ であるから、半径の小さなずれは復元力を受ける。
円軌道の公転角速度も $\dot\varphi=h/(mr_0^2)=\sqrt{\kappa/(mr_0^3)}$ で $\omega$ に等しい。動径の振動と公転の周期が一致するので、わずかにずれた軌道も閉じる（楕円軌道）。
:::
::::

::::exercise{level=advanced}
直線上の二質点を、相対速度に比例する抵抗をもつダンパーで結ぶ。$\mathcal R=c(\dot x_2-\dot x_1)^2/2$（$c>0$）とし、他に水平方向の力はない。各質点の力と、重心の速度、運動エネルギーの変化率を求めよ。
:::solution
速度で偏微分して負号を付けると $Q_1=c(\dot x_2-\dot x_1)$、$Q_2=-c(\dot x_2-\dot x_1)$ である。
それぞれ $m_i\ddot x_i=Q_i$ に従い、両辺を足せば全運動量が一定とわかる。したがって重心速度は初期値 $(m_1\dot x_1(0)+m_2\dot x_2(0))/(m_1+m_2)$ を保つ。
一方、力と速度の積を加えると

$$
\dot T=Q_1\dot x_1+Q_2\dot x_2=-c(\dot x_2-\dot x_1)^2=-2\mathcal R.
$$

二質点が同じ速さで同じ向きに動けば抵抗は消え、相対運動があれば運動エネルギーが減少する。
:::
::::

::::exercise{level=advanced}
水平面内で一定角速度 $\omega$ で回る滑らかな細い棒に質量 $m$ の小環を通す。極座標 $(r,\varphi)$ と拘束 $\varphi-\omega t=0$ を用い、未定乗数法で小環の運動と棒からの力を求めよ。棒を回す装置がする仕事率も求めよ。
:::solution
$L=\tfrac12m(\dot r^2+r^2\dot\varphi^2)$、拘束の係数は $a_r=0$、$a_\varphi=1$ である。:eqref[eq-am-lagrange-multiplier-equations] は

$$
m\ddot r-mr\dot\varphi^2=0,\qquad \dv{}{t}(mr^2\dot\varphi)=\lambda .
$$

$\dot\varphi=\omega$ を代入すると $\ddot r=\omega^2r$、すなわち $r=A\cosh\omega t+B\sinh\omega t$ であり、$\lambda=2mr\dot r\omega$ を得る。
$\lambda$ は $\varphi$ 方向の一般化力、すなわち棒が小環に及ぼす力の原点まわりのモーメントであり、棒に垂直な力の大きさは $\lambda/r=2m\omega\dot r$ である。
拘束は非定常なので、この力は実際の運動で仕事をする。仕事率は $\lambda\dot\varphi=2m\omega^2r\dot r=\dv{}{t}(m\omega^2r^2)$ であり、運動エネルギー $\tfrac12m(\dot r^2+\omega^2r^2)$ の増加率 $m\dot r\ddot r+m\omega^2r\dot r=2m\omega^2r\dot r$ に等しい。
:::
::::

::::exercise{level=advanced}
調和振動子 $L=\tfrac12m\dot x^2-\tfrac12m\omega^2x^2$ について、$x(0)=0$、$x(\tau)=X$ を満たす実軌道の作用を求めよ（$0<\omega\tau<\pi$）。$\omega\tau\to\pi$ で何が起こるか。
:::solution
$\sin\omega\tau\ne0$ なので、端点条件を満たす解は $x(t)=X\sin\omega t/\sin\omega\tau$ だけである。部分積分と運動方程式 $\ddot x=-\omega^2x$ を使うと

$$
S=\frac m2\int_0^\tau(\dot x^2-\omega^2x^2)\,\dd t=\frac m2\big[x\dot x\big]_0^\tau-\frac m2\int_0^\tau x(\ddot x+\omega^2x)\,\dd t=\frac m2X\dot x(\tau)
$$

であり、$\dot x(\tau)=X\omega\cos\omega\tau/\sin\omega\tau$ から

$$
S=\frac{m\omega X^2}{2}\cot\omega\tau
$$

を得る。$\omega\tau\to\pi$ では $X\ne0$ なら $\abs S\to\infty$ となる。$\omega\tau=\pi$ では $x(\tau)=0$ 以外の終点に達する解がなく、$X=0$ なら $c\sin\omega t$ がすべて解になる。これが :ref[exm-am-oscillator-action] の共役点である。
:::
::::

::::exercise{level=advanced}
平面上の粒子の Lagrangian $L=\tfrac12m(\dot x^2+\dot y^2)-V(x,y)$ に $\Lambda=\tfrac k2(x\dot y-y\dot x)$（$k$ は定数）を加える。$\Lambda$ は全時間微分か。運動方程式はどう変わるか。
:::solution
$\Lambda=A_x\dot x+A_y\dot y$、$A=(-ky/2,kx/2)$ の形である。:ref[prop-am-null-lagrangian] の証明から、全時間微分であるためには $\partial A_x/\partial y=\partial A_y/\partial x$ が必要だが、左辺は $-k/2$、右辺は $k/2$ で、$k\ne0$ なら成り立たない。よって全時間微分ではない。
実際、$E_x(\Lambda)=\dv{}{t}(-ky/2)-k\dot y/2=-k\dot y$、$E_y(\Lambda)=k\dot x$ であり、方程式は

$$
m\ddot x=-\pdv{V}{x}+k\dot y,\qquad m\ddot y=-\pdv{V}{y}-k\dot x
$$

となる。これは $z$ 方向の一様な磁場 $B=k/q$ の中の Lorentz 力 $q\vb{v}\times\vb{B}$ と同じであり、$\Lambda$ は対称ゲージ $\vb{A}=\tfrac12\vb{B}\times\vb{r}$ の $q\vb{v}\cdot\vb{A}$ に当たる。
:::
::::

## 確認テスト

::::quiz
$L$ に $\dd F(q,t)/\dd t$ を加えた。固定端での作用と一般化運動量について正しいものをすべて選べ。

- [x] 作用の第一変分は同じである
- [ ] 作用の値はどの端点でも必ず同じである
- [x] $p_i$ には $\partial F/\partial q^i$ が加わる

:::explanation
作用の差は端点での $F$ の差 $F(q(t_1),t_1)-F(q(t_0),t_0)$ であり、固定端ならその変分は零である。しかし値そのものは一般に変わる。速度で偏微分すると $\dd F/\dd t$ からは $\partial F/\partial q^i$ が残る。
:::
::::

::::quiz{answer="2.45" tol="0.001"}
Atwood の装置で $m_1=5\,\mathrm{kg}$、$m_2=3\,\mathrm{kg}$、$g=9.8\,\mathrm{m\,s^{-2}}$ とする。$m_1$ の下向き加速度を $\mathrm{m\,s^{-2}}$ で答えよ。

:::explanation
$L=\tfrac12(m_1+m_2)\dot x^2+(m_1-m_2)gx$ の方程式から $\ddot x=(m_1-m_2)g/(m_1+m_2)=2\times9.8/8=2.45$ である。
:::
::::

::::quiz
SI 系の電磁場中の荷電粒子に、$\phi,\vb{A}$ で表した Lagrangian を用いる。一般化運動量はどれか。

- [ ] $m\vb{v}-q\vb{A}$
- [x] $m\vb{v}+q\vb{A}$
- [ ] $m\vb{v}+q\vb{E}$

:::explanation
$L=\tfrac12m\vb{v}^2-q\phi+q\vb{v}\cdot\vb{A}$ の $q\vb{v}\cdot\vb{A}$ を速度の各成分で偏微分すると $qA_i$ が得られる。運動エネルギーの微分 $mv_i$ と合わせればよい。これは正準運動量であり、ゲージ変換で変わる。
:::
::::

::::quiz{answer="4.9" tol="0.001"}
質量 $2\,\mathrm{kg}$ の薄い輪が傾斜角 $30^\circ$ の斜面を滑らずに転がり下りる。$g=9.8\,\mathrm{m\,s^{-2}}$ として、斜面からの静止摩擦力の大きさを $\mathrm N$ で答えよ。

:::explanation
:ref[exm-am-rolling-hoop] の未定乗数法から、摩擦力は $\lambda=\tfrac12Mg\sin\alpha=\tfrac12\times2\times9.8\times0.5=4.9\,\mathrm N$ である。加速度は $\tfrac12g\sin\alpha=2.45\,\mathrm{m\,s^{-2}}$ である。
:::
::::

::::quiz
一定角速度 $\boldsymbol\Omega$ で回る座標系で Lagrangian を書いた。一般化速度について一次の項から生じる力はどれか。

- [ ] 遠心力 $-m\boldsymbol\Omega\times(\boldsymbol\Omega\times\vb{x}')$
- [x] Coriolis 力 $-2m\boldsymbol\Omega\times\dot{\vb{x}}'$
- [ ] 回転座標系では慣性力は Lagrangian に現れない

:::explanation
:ref[exm-am-rotating-frame] の $\widetilde L$ のうち一次の項 $m\dot{\vb{x}}'\cdot(\boldsymbol\Omega\times\vb{x}')$ が Coriolis 力を、速度によらない $\tfrac12m\abs{\boldsymbol\Omega\times\vb{x}'}^2$ が遠心力を与える。慣性力は時間に依存する点変換で速度を変換すると自動的に現れる。
:::
::::

::::quiz{answer="12"}
ばねにつないだ質点に抵抗力 $-c\dot x$ が働き、他に駆動力はない。$c=3\,\mathrm{kg/s}$、$\dot x=2\,\mathrm{m/s}$ の瞬間に、力学的エネルギーが減る速さを正の値で $\mathrm W$ を単位として答えよ。

:::explanation
$\mathcal R=\tfrac12c\dot x^2$ で、定常な系では力学的エネルギーの減少率は $2\mathcal R=c\dot x^2=3\times2^2=12\,\mathrm W$ である。
:::
::::

---
title: Newton 力学と拘束・一般化座標
description: 質点系の Newton 方程式と保存則から出発し、拘束の分類と積分可能性、一般化座標と配位空間、仮想仕事の原理と Lagrange 乗数による静力学を扱う。
notation:
  - symbol: '\vb{r}_a'
    meaning: 第 a 質点の位置ベクトル
  - symbol: '\vb{F}_a,\ \vb{R}_a'
    meaning: 第 a 質点に働く印加力と拘束力
  - symbol: '\vb{P},\ \vb{L}'
    meaning: 全運動量と原点まわりの全角運動量
  - symbol: 'q^i'
    meaning: 第 i 一般化座標
  - symbol: '\mathcal Q_t'
    meaning: 時刻 t に許される配位の集合
  - symbol: 'M_{ij}(q,t)'
    meaning: 一般化座標で表した質量行列
  - symbol: '\delta\vb{r}_a'
    meaning: 時刻を固定した第 a 質点の仮想変位
  - symbol: 'Q_i'
    meaning: 第 i 一般化座標に対応する一般化力
  - symbol: '\lambda_\alpha'
    meaning: 第 α 拘束に対する Lagrange 乗数
---

## この章の目標

- 質点系の Newton 方程式から運動量・角運動量・エネルギーの関係式を導き、重心の運動と内部の運動を分離する。
- 拘束をホロノミックと非ホロノミック、定常と非定常、両側と片側に分類し、速度に関する拘束が位置の拘束に積分できるかを判定する。
- 正則なホロノミック拘束の下で自由度を数え、一般化座標で速度と運動エネルギーを表す。
- 仮想変位と理想拘束を定義し、仮想仕事の原理と Lagrange 乗数を使って静力学のつり合いと反力を求める。

固定した棒の端に取り付けたおもりは、空間内の任意の方向に動けるわけではない。
[力学](/physics/mechanics/) の運動方程式でこの事実を表すには、棒がおもりに及ぼす力を未知量に加える。
一方、おもりが動ける方向に沿って式を立てると、その力を求める前に運動が決まる。
本章ではこの二つの記述を結びつける準備をし、次の二章で変分原理と Lagrange 形式に進む。
拘束を力学系の記述に取り入れる標準的な扱いは :cite[goldstein-cm] や :cite[yamamoto-nakamura-am] にある。

## 質点系の Newton 方程式

### 質点ごとに力を数える

慣性系に直交座標をとり、質量 $m_a>0$ が一定の $N$ 個の質点を考える。
番号 $a$ は $1$ から $N$ までを走る。位置 $\vb{r}_a(t)$ は時刻 $t$ について二度微分できるとする。
速度を $\vb{v}_a=\dot{\vb{r}}_a$ と書く。SI では位置、時間、質量の単位はそれぞれ $\mathrm{m},\mathrm{s},\mathrm{kg}$ である。

:::law[質点の運動方程式]{#law-am-newton}
第 $a$ 質点に働く合力を $\vb{f}_a$ とすると、

$$
m_a\ddot{\vb{r}}_a=\vb{f}_a,\qquad a=1,\ldots,N.
\label{eq-am-newton-system}
$$

力の単位は $\mathrm{N}=\mathrm{kg\,m\,s^{-2}}$ である。
:::

:eqref[eq-am-newton-system] は $3N$ 個の未知関数に対する二階の常微分方程式系である。
合力が位置・速度・時刻の関数として $\vb{f}_a=\vb{f}_a(\vb{r}_1,\ldots,\vb{r}_N,\vb{v}_1,\ldots,\vb{v}_N,t)$ と与えられ、それが連続微分可能であれば、初期位置と初速度を与えたとき解は存在する限り一意である（[解の存在と一意性](/math/ode/existence-uniqueness/)）。
これが古典力学の決定論の数学的な内容である。

運動が指定された面や棒から外れないように働く力を :term[拘束力]{en="constraint force" yomi="こうそくりょく"} という。
拘束力を $\vb{R}_a$、それ以外の力（印加力）を $\vb{F}_a$ とすれば、$\vb{f}_a=\vb{F}_a+\vb{R}_a$ である。
例えば滑らかな面上の粒子には、重力と面からの垂直抗力が働く。重力は既知でも、抗力の大きさは粒子の運動を調べるまでわからない。
拘束力は前もって位置と速度の関数として与えられないので、上の一意性定理をそのままでは適用できない。
拘束力を未知量として残すか、それとも方程式から消去するかが、本章と [Lagrange 形式](/physics/analytical-mechanics/lagrangian/) の主題である。

これとは別に、力を**誰が及ぼすか**によって分類できる。
系の外から及ぼす力を外力、系の中の別の質点から受ける力を内力と呼ぶ。
例えば二つのおもりを棒で結び、棒を質量のない部材として扱うと、棒が伝える力はおもり間の内力であり、同時に拘束力でもある。

### 運動量・角運動量・エネルギー

各質点の運動量を $m_a\vb{v}_a$ とし、系の全運動量、原点まわりの全角運動量、全運動エネルギーを

$$
\vb{P}=\sum_{a=1}^N m_a\vb{v}_a,\qquad
\vb{L}=\sum_{a=1}^N \vb{r}_a\times m_a\vb{v}_a,\qquad
T=\sum_{a=1}^N\tfrac12 m_a\vb{v}_a\cdot\vb{v}_a
$$

と定める。
$\vb{f}_{ab}$ を質点 $b$ が $a$ に及ぼす力、$\vb{f}^{\mathrm{ext}}_a$ を第 $a$ 質点に働く外力の合計とし、$\vb{f}_a=\vb{f}^{\mathrm{ext}}_a+\sum_{b\ne a}\vb{f}_{ab}$ と書く。

:::definition[作用反作用の法則]{#def-am-action-reaction}
すべての $a\ne b$ について $\vb{f}_{ab}=-\vb{f}_{ba}$ が成り立つとき、内力は :term[作用反作用の法則]{en="law of action and reaction" yomi="さようはんさようのほうそく"} の**弱い形**に従うという。
さらに各 $\vb{f}_{ab}$ が $\vb{r}_a-\vb{r}_b$ に平行であるとき、**強い形**に従うという。
:::

万有引力や Coulomb 力のように二点を結ぶ直線に沿って働く力は強い形に従う。
これに対し、運動する荷電粒子どうしの磁気的な力は一般に弱い形すら満たさない。そのときは電磁場自身が運動量を運ぶ（[電磁気学](/physics/electromagnetism/)）。

:::proposition[運動量・角運動量・運動エネルギーの時間変化]{#prop-am-momentum-theorems}
:eqref[eq-am-newton-system] に従う運動について、次が成り立つ。

1. 内力が弱い形に従えば、$\dot{\vb{P}}=\sum_{a=1}^N\vb{f}^{\mathrm{ext}}_a$ である。
2. 内力が強い形に従えば、$\dot{\vb{L}}=\sum_{a=1}^N\vb{r}_a\times\vb{f}^{\mathrm{ext}}_a$ である。
3. 力の性質によらず、$\dot T=\sum_{a=1}^N\vb{f}_a\cdot\vb{v}_a$ である。
:::

:::proof
1. :eqref[eq-am-newton-system] を全質点について足すと、内力の和は対 $\{a,b\}$ ごとにまとめられ、

$$
\dot{\vb{P}}
=\sum_{a=1}^N\vb{f}^{\mathrm{ext}}_a+
\sum_{1\le a<b\le N}(\vb{f}_{ab}+\vb{f}_{ba})
=\sum_{a=1}^N\vb{f}^{\mathrm{ext}}_a.
$$

2. $\dv{}{t}(\vb{r}_a\times m_a\vb{v}_a)=\vb{v}_a\times m_a\vb{v}_a+\vb{r}_a\times m_a\dot{\vb{v}}_a=\vb{r}_a\times\vb{f}_a$ である。和をとると内力の寄与は対ごとに

$$
\vb{r}_a\times\vb{f}_{ab}+\vb{r}_b\times\vb{f}_{ba}
=(\vb{r}_a-\vb{r}_b)\times\vb{f}_{ab}
$$

となる。最初の等号で弱い形を使い、強い形によりこれは零である。

3. $\dv{}{t}\big(\tfrac12m_a\vb{v}_a\cdot\vb{v}_a\big)=m_a\dot{\vb{v}}_a\cdot\vb{v}_a=\vb{f}_a\cdot\vb{v}_a$ を足せばよい。
:::

外力がなければ全運動量と全角運動量は保存する。
この二つの保存則が空間の一様性と等方性から導かれることは [対称性と保存則](/physics/analytical-mechanics/symmetries/) で示す。

### 重心の運動と内部の運動

全質量を $M=\sum_{a=1}^N m_a$、重心の位置を $\vb{r}_{\mathrm G}=M^{-1}\sum_{a=1}^N m_a\vb{r}_a$ とする。
重心に対する相対位置を $\vb{r}'_a=\vb{r}_a-\vb{r}_{\mathrm G}$、相対速度を $\vb{v}'_a=\dot{\vb{r}}'_a$ と書く。

:::proposition[重心運動と相対運動への分解]{#prop-am-com-decomposition}
重心の速度を $\vb{v}_{\mathrm G}=\dot{\vb{r}}_{\mathrm G}$ とすると

$$
\begin{aligned}
\vb{P}&=M\vb{v}_{\mathrm G},\\
\vb{L}&=\vb{r}_{\mathrm G}\times M\vb{v}_{\mathrm G}+\sum_{a=1}^N\vb{r}'_a\times m_a\vb{v}'_a,\\
T&=\tfrac12M\vb{v}_{\mathrm G}\cdot\vb{v}_{\mathrm G}+\sum_{a=1}^N\tfrac12m_a\vb{v}'_a\cdot\vb{v}'_a .
\end{aligned}
$$

特に内力が弱い形に従えば $M\ddot{\vb{r}}_{\mathrm G}=\sum_{a=1}^N\vb{f}^{\mathrm{ext}}_a$ である。
:::

:::proof
重心の定義から $\sum_a m_a\vb{r}'_a=\sum_a m_a\vb{r}_a-M\vb{r}_{\mathrm G}=0$ であり、時間微分して $\sum_a m_a\vb{v}'_a=0$ でもある。
$\vb{P}=\sum_a m_a(\vb{v}_{\mathrm G}+\vb{v}'_a)=M\vb{v}_{\mathrm G}$ である。
角運動量では

$$
\vb{L}=\sum_{a=1}^N(\vb{r}_{\mathrm G}+\vb{r}'_a)\times m_a(\vb{v}_{\mathrm G}+\vb{v}'_a)
=\vb{r}_{\mathrm G}\times M\vb{v}_{\mathrm G}
+\vb{r}_{\mathrm G}\times\sum_a m_a\vb{v}'_a
+\Big(\sum_a m_a\vb{r}'_a\Big)\times\vb{v}_{\mathrm G}
+\sum_a\vb{r}'_a\times m_a\vb{v}'_a
$$

の中央の二項が消える。運動エネルギーでも同様に交差項 $\vb{v}_{\mathrm G}\cdot\sum_a m_a\vb{v}'_a$ が消える。
最後の主張は :ref[prop-am-momentum-theorems] の 1 に $\vb{P}=M\vb{v}_{\mathrm G}$ を代入したものである。
:::

運動エネルギーの分解は König の定理とも呼ばれる。
重心だけを追うと内力は式から消える。しかし、この式だけでは質点間の距離がどう変わるかは決まらない。

:::example[同じ重心運動をする二つの質点系]
水平面上の質量が等しい二つの質点を考え、面内の外力を零とする。
一方の系では両者を軽いばねで結び、もう一方では結ばずに放す。どちらも重心の速度は一定である。
違いを調べるには相対位置 $\vb{s}=\vb{r}_2-\vb{r}_1$ を観測する。
各質点の質量を $m$ とすれば、二つの運動方程式の差は

$$
m\ddot{\vb{s}}=\vb{f}_{21}-\vb{f}_{12}=2\vb{f}_{21}
$$

である。ばねを伸ばして放した直後には右辺は零でなく、結ばない系では零である。
質量が異なる場合は、換算質量 $\mu=m_1m_2/(m_1+m_2)$ を使って $\mu\ddot{\vb{s}}=\vb{f}_{21}$ となる。これは [中心力問題と Kepler 問題](/physics/analytical-mechanics/central-force/) の出発点である。
:::

### 保存力と力学的エネルギー

:::definition[保存力とポテンシャル]{#def-am-conservative-force}
領域 $D\subset\R^3$ 上の力の場 $\vb{f}(\vb{r})$ に対して $C^1$ 級の関数 $V$ があり、$\vb{f}=-\grad V$ となるとき、$\vb{f}$ を :term[保存力]{en="conservative force" yomi="ほぞんりょく"}、$V$ をその :term[ポテンシャルエネルギー]{en="potential energy" yomi="ぽてんしゃるえねるぎー"}（位置エネルギー）という。
:::

保存力が曲線 $\gamma$ に沿ってする仕事は、線積分の基本定理により $\int_\gamma\vb{f}\cdot\dd\vb{r}=V(\text{始点})-V(\text{終点})$ であり、経路によらない。
星形の領域では $\nabla\times\vb{f}=0$ が保存力であることと同値である（[線積分・面積分とベクトル解析](/math/analysis/vector-calculus/)）。

内力については、各対 $\{a,b\}$ に距離だけの関数 $V_{ab}(s)$ が与えられ、

$$
\vb{f}_{ab}=-\nabla_a V_{ab}(\abs{\vb{r}_a-\vb{r}_b})
=-V'_{ab}(s)\frac{\vb{r}_a-\vb{r}_b}{s},\qquad s=\abs{\vb{r}_a-\vb{r}_b}
$$

となる場合が典型的である。ここで $\nabla_a$ は $\vb{r}_a$ の三成分に関する勾配である。
$\nabla_b$ で微分すると符号だけが変わるので $\vb{f}_{ba}=-\vb{f}_{ab}$ であり、しかも力は $\vb{r}_a-\vb{r}_b$ に平行である。
したがってこの形の内力は作用反作用の法則の強い形に従う。

:::proposition[質点系の力学的エネルギー保存]{#prop-am-system-energy}
外力が時刻によらないポテンシャル $V_a(\vb{r}_a)$ から $\vb{f}^{\mathrm{ext}}_a=-\nabla_a V_a$ と導かれ、内力が上の対ポテンシャル $V_{ab}$ から導かれ、ほかに力が働かないとする。このとき

$$
E=T+\sum_{a=1}^NV_a(\vb{r}_a)+\sum_{1\le a<b\le N}V_{ab}(\abs{\vb{r}_a-\vb{r}_b})
$$

は運動の間一定である。
:::

:::proof
ポテンシャルの部分を $V$ と書いて軌道に沿って微分すると、

$$
\dot V=\sum_{a=1}^N\nabla_aV_a\cdot\vb{v}_a
+\sum_{a<b}\big(\nabla_aV_{ab}\cdot\vb{v}_a+\nabla_bV_{ab}\cdot\vb{v}_b\big)
=-\sum_{a=1}^N\vb{f}^{\mathrm{ext}}_a\cdot\vb{v}_a-\sum_{a<b}\big(\vb{f}_{ab}\cdot\vb{v}_a+\vb{f}_{ba}\cdot\vb{v}_b\big)
=-\sum_{a=1}^N\vb{f}_a\cdot\vb{v}_a .
$$

:ref[prop-am-momentum-theorems] の 3 と合わせれば $\dot E=0$ である。
:::

拘束力が働く場合には、拘束力が実際の運動でする仕事も勘定に入る。
定常な理想拘束ではこの仕事が零になることを :ref[prop-am-ideal-no-work] で示す。

### 既知の運動から反力を求める

滑らかな水平板を鉛直方向に動かし、高さを既知の関数 $h(t)$ とする。上向きを正とする。
質点が板に乗り続ける間、その鉛直位置と加速度は $z=h(t)$、$\ddot z=\ddot h(t)$ である。
垂直抗力の上向き成分を $R$、重力加速度の大きさを $g>0$ とすると、

$$
m\ddot z=R-mg,\qquad R=m(g+\ddot h).
$$

この結果から、板で押すだけで指定の運動が実現するかを判定できる。
$\ddot h<-g$ では上式は下向きの抗力を要求する。接着していない板は引けないため、この場合は接触を仮定した方程式を使えない。
このように、拘束を保つのに必要な拘束力は運動方程式から**事後的に**決まる。

## 拘束条件の分類

### 位置の条件と時間の役割

:::definition[ホロノミック拘束]{#def-am-holonomic}
時刻と位置の等式

$$
f_\alpha(\vb{r}_1,\ldots,\vb{r}_N,t)=0,\qquad \alpha=1,\ldots,k
\label{eq-am-holonomic}
$$

により表される拘束を :term[ホロノミック拘束]{en="holonomic constraint" yomi="ほろのみっくこうそく"} という。
速度に対する条件が、積分定数を初期位置で決めればこの形に帰着する場合も含む。
本章では特に断らない限り、$f_\alpha$ は $C^2$ 級で、等式が運動中ずっと成り立つ拘束を扱う。
:::

時間への依存は別の分類を与える。拘束自体が時間とともに変わらない場合を :term[定常拘束]{en="scleronomic constraint" yomi="ていじょうこうそく"}、
時刻に陽に依存して変わる場合を :term[非定常拘束]{en="rheonomic constraint" yomi="ひていじょうこうそく"} という。
ここでは実験室に固定した座標で判定する。例えば $z=h(t)$ は、$h$ が定数なら定常、上下する板なら非定常であり、どちらもホロノミックである。

:::example[平面振り子]{#exm-am-pendulum-constraint}
鉛直面内で、質点を長さ $\ell>0$ の軽い剛体棒で原点に結ぶ。鉛直上向きを $y$ 軸にとる。
棒の長さが一定である条件は

$$
f(x,y)=\tfrac12(x^2+y^2-\ell^2)=0
$$

である。したがってこれは定常なホロノミック拘束である。
支点の水平位置を $b(t)$ と指定する場合は

$$
f(x,y,t)=\tfrac12\big((x-b(t))^2+y^2-\ell^2\big)=0
$$

となる。許される位置は依然として円をなすが、時刻とともに円の中心が動く。
:::

:::example[剛体]{#exm-am-rigid-constraint}
剛体を、互いの距離が変わらない質点の集まりとみなす。拘束は

$$
\abs{\vb{r}_a-\vb{r}_b}^2-c_{ab}^2=0,\qquad 1\le a<b\le N
$$

であり、定常なホロノミック拘束である。$c_{ab}>0$ は定数である。
:::

### 両側拘束と片側拘束

:::definition[両側拘束と片側拘束]{#def-am-unilateral}
等式で表される拘束を :term[両側拘束]{en="bilateral constraint" yomi="りょうがわこうそく"}、
$g(\vb{r}_1,\ldots,\vb{r}_N,t)\ge0$ の形の不等式で表される拘束を :term[片側拘束]{en="unilateral constraint" yomi="かたがわこうそく"} という。
:::

容器に閉じ込めた気体分子や、球面の外側に置いた粒子は片側拘束を受ける。
片側拘束では、$g>0$ の間は拘束力が働かず、$g=0$ の間は両側拘束と同じように運動する。
ただし後者では、拘束力が許される側を向いていること、すなわち滑らかな面なら $\vb{R}=\lambda\grad g$ で $\lambda\ge0$ となることを確かめなければならない。
$\lambda$ が負になる時刻が来れば、粒子は面を離れる（演習問題を参照）。

糸で質点をつるす場合も、支点からの距離は糸の長さ**以下**であればよいので片側拘束である。
したがって :ref[exm-am-pendulum-constraint] の等式で運動を記述できるのは、糸が張り続ける区間に限られる。
以下で「拘束」と言えば、特に断らない限り両側拘束を指す。

### 速度に関する式は積分できるか

位置の条件 :eqref[eq-am-holonomic] を運動に沿って微分する。
$\nabla_a f_\alpha$ を質点 $a$ の三つの位置成分についての勾配とすれば、連鎖律から

$$
\sum_{a=1}^N \nabla_a f_\alpha\cdot\vb{v}_a
+\pdv{f_\alpha}{t}=0.
\label{eq-am-admissible-velocity}
$$

逆に、この形の全微分が零なら、軌道に沿って $f_\alpha$ は一定である。
初期位置で $f_\alpha=0$ ならば、以後も位置の条件を満たす。

一般に、位置成分を順に並べた $X=(X_1,\ldots,X_{3N})$ を使って、速度に関する拘束を次の形で考える。

:::definition[Pfaff 型拘束]{#def-am-pfaffian}
係数 $a_{\alpha j}(X,t)$、$a_{\alpha0}(X,t)$ を $C^1$ 級の関数として

$$
\sum_{j=1}^{3N} a_{\alpha j}(X,t)\dot X_j+a_{\alpha0}(X,t)=0,\qquad\alpha=1,\ldots,k'
\label{eq-am-pfaffian}
$$

と書ける拘束を、速度について線形な拘束、または :term[Pfaff 型拘束]{en="Pfaffian constraint" yomi="ぷふぁふがたこうそく"} という。
係数行列 $(a_{\alpha j})$ の階数が $k'$ のとき、拘束は独立であるという。
:::

ホロノミック拘束の微分 :eqref[eq-am-admissible-velocity] は $a_{\alpha j}=\partial f_\alpha/\partial X_j$、$a_{\alpha0}=\partial f_\alpha/\partial t$ とした Pfaff 型拘束である。
逆に、ある Pfaff 型拘束について、各 $(X,t)$ で同じ速度の集合を定めるホロノミック拘束の微分が存在するとき、その拘束は :term[積分可能な拘束]{en="integrable constraint" yomi="せきぶんかのうなこうそく"} であるという。
このとき Pfaff 型拘束を満たす運動は、$f_\alpha$ を一定に保つ運動にほかならない。

:::example[積分因子]
平面上の粒子に $y\dot x-x\dot y=0$ を課す。$x>0$ の領域でこの式に $-1/x^2$ を掛けると

$$
-\frac{y\dot x-x\dot y}{x^2}=\dv{}{t}\left(\frac{y}{x}\right)=0
$$

となる。したがって $y/x$ は一定であり、粒子は原点を通る直線上を動く。
同じ速度の集合を定める式に掛ける関数 $-1/x^2$ を :term[積分因子]{en="integrating factor" yomi="せきぶんいんし"} という。
元の式そのものは全微分の形をしていなくても、積分可能なことがある。
:::

例えば車輪が一本の直線上だけを転がるとき、中心の速さに符号を付けたものを $\dot x$、回転角速度を $\dot\phi$ とすれば、半径 $a$ に対して $\dot x=a\dot\phi$ である。
よって角度を連続的に数えた区間では $x-a\phi$ が一定になり、この転がり条件は積分可能である。

これに対し、許される速度を指定しても、位置の等式に積分できない拘束を :term[非ホロノミック拘束]{en="nonholonomic constraint" yomi="ひほろのみっくこうそく"} という。
その違いを確かめるために、次に車輪の向きも変えられるようにする。

:::example[鉛直に保った円板の転がり]{#exm-am-rolling-disk}
半径 $a$ の薄い円板が水平面上を滑らずに転がる。
円板は常に鉛直とし、中心の水平位置を $(x,y)$、進行方向の方位角を $\psi$、円板の自転角を $\phi$ とする。
角度の正の向きは、自転によって中心が $(\cos\psi,\sin\psi)$ の向きに進むようにとる。
接地点の速度が零になるには、進行方向の中心速度が $a\dot\phi$、横方向の中心速度が零でなければならない。したがって

$$
\dot x=a\cos\psi\,\dot\phi,\qquad
\dot y=a\sin\psi\,\dot\phi.
\label{eq-am-rolling-disk}
$$

中心の高さはすでに $a$ に固定した。円板を鉛直に保つこと自体も一つの拘束であり、それには車輪を支える枠のような支持機構からの力のモーメントが必要である。
残る四つの座標のうち、瞬間に独立に選べる速度は $\dot\psi,\dot\phi$ の二つである。
:::

なぜこの二式から、位置の座標を二つ消去できないのだろうか。
方位を $\psi=0$ として距離 $s$ を転がし、その場で方位を $\beta$ だけ変える。
次に同じ距離だけ逆向きに転がし、その場で方位を元に戻す。理想化した点接触では、その場での方位変更は :eqref[eq-am-rolling-disk] を満たす。
$\phi,\psi$ は両方とも初期値に戻るが、中心は

$$
(\Delta x,\Delta y)=s(1-\cos\beta,-\sin\beta)
$$

だけ移る。特に $\beta$ を小さくしても横方向のずれが生じる。
もし二つの独立な位置の条件 $f_1=f_2=0$ に積分できるなら、出発点と同じ $(\psi,\phi)$ をもつ配位はその条件から（局所的に）一点に決まるはずであり、上の結果と矛盾する。

積分できないことを局所的にも確かめよう。位置に関する保存条件 $H(x,y,\psi,\phi)$ があるとすれば、すべての許される速度で $\dot H=0$ である。
$\dot\psi,\dot\phi$ の係数を別々に零とおくと、

$$
H_\psi=0,\qquad aH_x\cos\psi+aH_y\sin\psi+H_\phi=0.
$$

最初の式により、二つ目の式の $H_x,H_y,H_\phi$ は $\psi$ に依存しない。
区間上で $1,\cos\psi,\sin\psi$ は一次独立なので、三つの偏導関数はすべて零である。
したがって局所的にも、転がり条件を与える非定数の位置の積分は一つも存在しない。

拘束が一つで変数が三つの場合には、次の簡単な判定法がある。

:::proposition[積分可能性の必要条件]{#prop-am-pfaff-test}
$\R^3$ の領域 $D$ 上の $C^1$ 級ベクトル場 $\vb{a}(\vb{x})\ne0$ に対し、拘束 $\vb{a}(\vb{x})\cdot\dot{\vb{x}}=0$ を考える。
点 $\vb{x}_0\in D$ の近傍で $\grad H\ne0$ を満たす $C^2$ 級関数 $H$ があり、そこで $\vb{a}\cdot\dot{\vb{x}}=0$ と $\grad H\cdot\dot{\vb{x}}=0$ が同じ速度の集合を定めるならば、その近傍で

$$
\vb{a}\cdot(\nabla\times\vb{a})=0
\label{eq-am-pfaff-test}
$$

である。
:::

:::proof
二つの式はどちらも $\dot{\vb{x}}$ に一つの一次条件を課し、同じ平面を定める。法線ベクトルは平行なので、零にならない関数 $\mu$ によって $\grad H=\mu\vb{a}$ と書ける。
$\mu=\grad H\cdot\vb{a}/\abs{\vb{a}}^2$ であるから $\mu$ は $C^1$ 級である。勾配の回転は零なので

$$
0=\nabla\times(\mu\vb{a})=\grad\mu\times\vb{a}+\mu\,\nabla\times\vb{a}.
$$

両辺と $\vb{a}$ の内積をとると、$\vb{a}\cdot(\grad\mu\times\vb{a})=0$ であるから $\mu\,\vb{a}\cdot(\nabla\times\vb{a})=0$ となる。$\mu\ne0$ より主張を得る。
:::

逆に、:eqref[eq-am-pfaff-test] が成り立てば局所的に積分因子が存在する。これは Frobenius の定理の特別な場合であり、証明には多様体上の分布の理論を用いるので本章では扱わない（[分布と Frobenius の定理](/math/manifolds/distributions/)、:cite[lee-smooth]）。

:::example[スケートの刃]{#exm-am-skate}
氷上のスケートの刃を、平面上の点 $(x,y)$ と刃の向き $\psi$ で表す。刃は横滑りしないので、刃に垂直な方向の速度は零である。

$$
-\sin\psi\,\dot x+\cos\psi\,\dot y=0.
$$

$(x,y,\psi)$ を $\R^3$ の座標とみなすと $\vb{a}=(-\sin\psi,\cos\psi,0)$ であり、

$$
\nabla\times\vb{a}=\left(\pdv{a_3}{y}-\pdv{a_2}{\psi},\ \pdv{a_1}{\psi}-\pdv{a_3}{x},\ \pdv{a_2}{x}-\pdv{a_1}{y}\right)
=(\sin\psi,-\cos\psi,0),
$$

$\vb{a}\cdot(\nabla\times\vb{a})=-1\ne0$ となる。よってこの拘束は非ホロノミックである。
実際、刃を前後に滑らせることと向きを変えることを組み合わせれば、任意の $(x,y,\psi)$ に到達できる。自動車の縦列駐車はこの性質を利用している。
:::

## 一般化座標と配位空間

### 正則な拘束と自由度

固定した時刻 $t$ における**配位**とは系全体の位置 $X=(X_1,\ldots,X_{3N})$ の指定である。その時刻に許される配位全体を :term[配位空間]{en="configuration space" yomi="はいいくうかん"} $\mathcal Q_t$ と呼ぶ。
正則なホロノミック拘束では、次のように局所的な座標の数が決まる。

:::proposition[独立な位置拘束と自由度]{#prop-am-constraint-count}
式 :eqref[eq-am-holonomic] の $f_\alpha$ を $C^1$ 級とする。許された点で $k\times 3N$ 行列
$D_X f=(\partial f_\alpha/\partial X_j)$ の階数が $k$ であれば、その近傍の配位は

$$
n=3N-k
$$

個の独立な実数で記述できる。この局所的な座標の数を、ホロノミック系の :term[自由度]{en="degree of freedom" yomi="じゆうど"} という。
:::

:::proof
階数が $k$ なので、$k$ 列を選んで作る正方行列の行列式を零でないようにできる。
変数を並べ替え、その列に対応する位置成分を $u$、残る $3N-k$ 個を $w$ とする。
陰関数定理（[逆関数定理と陰関数定理](/math/analysis/implicit-function/)）により、$f(u,w,t)=0$ はこの点の近くでただ一つの $C^1$ 級関数 $u=u(w,t)$ を定める。
したがって時刻を固定したとき、$w$ を近傍で自由に選べば配位が一意に決まり、すべての近傍内の配位がこの形で得られる。
:::

この数え方で重要なのは拘束の独立性である。
例えば $f_1=z,\ f_2=2z$ とすれば、方程式は二つあっても零にすべき成分は一つである。
さらに $f=z^2$ は同じ平面を表すのに、平面上で勾配が消える。
このような表示では :ref[prop-am-constraint-count] の仮定を確認できないため、$f=z$ と書き直してから数える。

:::example[剛体の自由度]
:ref[exm-am-rigid-constraint] の拘束の数 $N(N-1)/2$ は、$N\ge 8$ で $3N$ を超える。したがってこれらの拘束は独立ではなく、自由度を $3N-N(N-1)/2$ と数えてはならない。
同一直線上にない三点 $\vb{r}_1,\vb{r}_2,\vb{r}_3$ を選ぶと、それらの九成分には三つの独立な距離の条件が課されるので自由度は $6$ である。
他の点は三点からの距離で鏡像の二点のどちらかに決まり、連続な運動では一方に固定される。よって剛体の自由度は $6$（重心の並進 $3$ と回転 $3$）である。
詳しくは [剛体の運動](/physics/analytical-mechanics/rigid-body/) で扱う。
:::

### 運動を位置座標で表す

:::definition[一般化座標]{#def-am-generalized-coordinates}
許される配位を局所的に一意に指定し、位置を $C^2$ 級の関数で

$$
\vb{r}_a=\vb{r}_a(q^1,\ldots,q^n,t),\qquad a=1,\ldots,N
$$

と表す独立なパラメータを :term[一般化座標]{en="generalized coordinates" yomi="いっぱんかざひょう"} という。
正則な座標では、$q^i$ に関する位置の偏導関数を全質点分並べた $3N$ 成分の $n$ 本の列 $\partial X/\partial q^i$ が一次独立である。
上付きの $i$ は座標の番号を表す。
:::

:ref[prop-am-constraint-count] の証明で得た $w$ は一般化座標の一例である。
しかし実際には、拘束の幾何に合わせて角度などを選ぶほうが便利である。
例えば :ref[exm-am-pendulum-constraint] の固定支点の場合、下向きの鉛直線から測った角 $\theta$ を使うと

$$
(x,y)=(\ell\sin\theta,-\ell\cos\theta),\qquad
\pdv{\vb{r}}{\theta}=\ell(\cos\theta,\sin\theta).
$$

角 $\theta$ と $\theta+2\pi$ は同じ配位を表すので、振り子の配位空間は円 $S^1$ である。
軌道に沿って角を連続して数えると何周したかも記録できるが、その周回数は瞬間の配位を区別しない。

```tikz caption="棒の長さを固定した平面振り子。許される変位は円に接する。"
\begin{tikzpicture}[>=Stealth]
\coordinate (O) at (0,0);
\coordinate (P) at (-55:2.2);
\draw[dashed] (O) -- (0,-2.4);
\draw[gray] (-115:2.2) arc (-115:-35:2.2);
\draw[thick] (O) -- node[right] {$\ell$} (P);
\draw[->] (0,-0.72) arc (-90:-55:0.72);
\node at (-72:1) {$\theta$};
\draw[->,thick] (P) -- ++(35:1.15) node[right] {$\delta\mathbf r$};
\fill (O) circle (1.6pt) node[above] {$O$};
\fill (P) circle (2.5pt);
\end{tikzpicture}
```

平面内を自由に動く二重振り子なら、二本の棒が鉛直下向きとなす角を別々に指定できる。
棒同士の衝突を考えなければ、配位空間は二つの円の直積 $S^1\times S^1$ である。

一本の棒を空間の任意の向きに振れるようにした :term[球面振り子]{en="spherical pendulum" yomi="きゅうめんふりこ"} では、質点の配位空間は球面 $S^2$ となる。
鉛直下向きから測った極角 $\theta$ と方位角 $\varphi$ を使えば

$$
\vb{r}=\ell(\sin\theta\cos\varphi,\ \sin\theta\sin\varphi,\ -\cos\theta),\qquad
\pdv{\vb{r}}{\theta}\times\pdv{\vb{r}}{\varphi}=\ell^2\sin\theta\,\vb{e}
$$

である。$\vb{e}$ は単位ベクトルである。したがって二本の列が一次独立なのは $\sin\theta\ne0$ のときに限る。
最下点 $\theta=0$ では方位を変えても点は動かず、この座標は正則でない。
この例が示すように、配位空間全体を一組の正則な座標で覆えるとは限らず、一般には複数の座標系を貼り合わせる必要がある。この見方を整理したものが [多様体](/math/manifolds/) の概念である。

### 許される速度の集合

一般化座標を時刻の関数として位置の式を微分すると、

$$
\vb{v}_a=\sum_{i=1}^n\pdv{\vb{r}_a}{q^i}\dot q^i+\pdv{\vb{r}_a}{t}.
\label{eq-am-coordinate-velocity}
$$

$\dot q^i$ を :term[一般化速度]{en="generalized velocity" yomi="いっぱんかそくど"} という。
最後の項は、座標値を一定に保っても生じる速度である。
例えば水平移動する支点に対して $\theta$ を測ると、

$$
\vb{r}=(b(t)+\ell\sin\theta,-\ell\cos\theta),\qquad
\vb{v}=(\dot b,0)+\ell(\cos\theta,\sin\theta)\dot\theta.
$$

したがって非定常拘束では、許される速度の集合は、固定時刻の接ベクトル全体を支点の速度などの分だけ平行移動したものになる。
線形方程式でいえば、:eqref[eq-am-admissible-velocity] の一般解は一つの特解と斉次方程式の一般解の和である（:ref[prop-linalg-solution-structure]）。

ある時刻の系の**状態**は、配位 $q$ だけでなく一般化速度 $\dot q$ も合わせて指定して初めて決まる。
Newton 方程式が二階であるため、$(q,\dot q)$ を与えれば以後の運動が決まると期待される。これは [Lagrange 形式](/physics/analytical-mechanics/lagrangian/) で確かめる。

### 一般化座標で表した運動エネルギー

:eqref[eq-am-coordinate-velocity] を $T=\sum_a\tfrac12m_a\vb{v}_a\cdot\vb{v}_a$ に代入して展開すると、

$$
T=\frac12\sum_{i=1}^n\sum_{j=1}^nM_{ij}(q,t)\dot q^i\dot q^j+\sum_{i=1}^nb_i(q,t)\dot q^i+T_0(q,t)
\label{eq-am-kinetic-generalized}
$$

となる。ここで係数は

$$
M_{ij}=\sum_{a=1}^Nm_a\pdv{\vb{r}_a}{q^i}\cdot\pdv{\vb{r}_a}{q^j},\qquad
b_i=\sum_{a=1}^Nm_a\pdv{\vb{r}_a}{q^i}\cdot\pdv{\vb{r}_a}{t},\qquad
T_0=\frac12\sum_{a=1}^Nm_a\pdv{\vb{r}_a}{t}\cdot\pdv{\vb{r}_a}{t}
$$

である。行列 $(M_{ij})$ を :term[質量行列]{en="mass matrix" yomi="しつりょうぎょうれつ"} という。
:eqref[eq-am-kinetic-generalized] の三つの項は一般化速度についてそれぞれ二次・一次・零次の同次式である。

:::proposition[質量行列の正定値性]{#prop-am-mass-matrix}
一般化座標が正則な点では、質量行列 $(M_{ij}(q,t))$ は実対称な正定値行列である。
特に逆行列をもつ。
:::

:::proof
対称性は定義から明らかである。任意の $\xi\in\R^n$ に対して

$$
\sum_{i,j}M_{ij}\xi^i\xi^j=\sum_{a=1}^Nm_a\Big|\sum_{i=1}^n\pdv{\vb{r}_a}{q^i}\xi^i\Big|^2\ge0
$$

であり、等号が成り立てば、$m_a>0$ よりすべての $a$ で $\sum_i(\partial\vb{r}_a/\partial q^i)\xi^i=0$ である。
これは $3N$ 成分の列 $\partial X/\partial q^i$ の一次結合が零であることを意味し、正則性により $\xi=0$ となる。
また $M\xi=0$ なら $\sum_{i,j}M_{ij}\xi^i\xi^j=0$ なので $\xi=0$ であり、$M$ は単射、したがって正則である。
:::

拘束が定常で、一般化座標の定義式が時刻に陽に依存しなければ、$b_i=0$、$T_0=0$ であり、運動エネルギーは一般化速度の正定値な二次形式になる。

:::example[運動エネルギーの例]
1. 平面の極座標 $(r,\varphi)$ では $\vb{r}=(r\cos\varphi,r\sin\varphi)$ から $M=\diag(m,mr^2)$、$b=0$、$T_0=0$ である。$r=0$ で正則性が失われる。
2. 水平に動く支点をもつ振り子では、$M_{\theta\theta}=m\ell^2$、$b_\theta=m\ell\dot b\cos\theta$、$T_0=\tfrac12m\dot b^2$ である。
3. 球面振り子では $M=\diag(m\ell^2,m\ell^2\sin^2\theta)$ であり、最下点 $\sin\theta=0$ で正定値性が失われる。これは座標の特異性による。
:::

## 仮想変位と仮想仕事の原理

### 時刻を固定して変位をとる

:::definition[仮想変位]{#def-am-virtual-displacement}
ホロノミック系で、時刻 $t$ を固定し、許される配位を通る曲線 $q(\eps)$ をとる。
$\eps=0$ における接ベクトルを :term[仮想変位]{en="virtual displacement" yomi="かそうへんい"} と呼ぶ。
$\delta q^i=(\dd q^i/\dd\eps)_{\eps=0}$ と書けば

$$
\delta\vb{r}_a=\sum_{i=1}^n\pdv{\vb{r}_a}{q^i}\delta q^i,\qquad
\sum_{a=1}^N\nabla_a f_\alpha\cdot\delta\vb{r}_a=0.
\label{eq-am-virtual-displacement}
$$

右の式は $f_\alpha(\vb{r}_1(\eps),\ldots,\vb{r}_N(\eps),t)=0$ を $\eps$ で微分したものである。
:::

時刻を固定しているので、:eqref[eq-am-coordinate-velocity] の $\partial\vb{r}_a/\partial t$ に当たる項は現れない。

:::proposition[仮想変位の全体]{#prop-am-virtual-space}
独立なホロノミック拘束の下で、正則な一般化座標をとる。各 $(q,t)$ で次の二つの $3N$ 次元ベクトルの集合は一致し、ともに $n$ 次元の部分空間である。

1. 一般化座標による仮想変位 $\delta X=\sum_i(\partial X/\partial q^i)\delta q^i$ 全体
2. 一次方程式 $D_Xf\,\delta X=0$ の解全体
:::

:::proof
1 は列 $\partial X/\partial q^i$ の張る空間で、列の一次独立性から次元 $n$ をもつ。
恒等式 $f(X(q,t),t)=0$ を $q^i$ で微分すると $D_Xf\,(\partial X/\partial q^i)=0$ なので、1 は 2 に含まれる。
2 は階数 $k$ の行列の核であり、次元公式（:ref[thm-linalg-rank-nullity]）により次元は $3N-k=n$ である。
次元の等しい部分空間の包含関係なので両者は一致する。
:::

位置が同じ二つの許容運動を考えると、速度の差では :eqref[eq-am-coordinate-velocity] の $\partial\vb{r}_a/\partial t$ が相殺する。
そのため仮想変位の向きは、許される速度どうしの差の向きからも読み取れる。
水平な板が上昇する例では、実際の変位には鉛直成分があるが、時刻を固定した板上の仮想変位は水平である。

速度について線形な非ホロノミック拘束にも、この「速度どうしの差」を用いる。
すなわち Pfaff 型拘束 :eqref[eq-am-pfaffian] の仮想変位を、$a_{\alpha0}$ を落とした

$$
\sum_{j=1}^{3N}a_{\alpha j}(X,t)\,\delta X_j=0,\qquad\alpha=1,\ldots,k'
\label{eq-am-pfaffian-virtual}
$$

の解と定める。
例えば円板の :eqref[eq-am-rolling-disk] で許される仮想変位は

$$
\delta x=a\cos\psi\,\delta\phi,\qquad
\delta y=a\sin\psi\,\delta\phi
$$

に従う。これはその瞬間の転がりに関する規則であって、四つの位置座標を二つに減らす規則ではない。

### 理想拘束と一般化力

力が許される変位の方向にどれだけ働くかを測るため、力と仮想変位の内積をとる。
各質点でとった内積の和を、その力の :term[仮想仕事]{en="virtual work" yomi="かそうしごと"} という。

:::definition[理想拘束]{#def-am-ideal-constraint}
すべての許される仮想変位に対して、拘束力の仮想仕事が

$$
\sum_{a=1}^N\vb{R}_a\cdot\delta\vb{r}_a=0
\label{eq-am-ideal-work}
$$

となる拘束を :term[理想拘束]{en="ideal constraint" yomi="りそうこうそく"} という。
:::

滑らかな面からの抗力は面に垂直で、仮想変位は面に接するため、両者の内積は零である。
さらに、質点を結ぶ棒では、両端に働く力の仮想仕事を**合計**して確かめる必要がある。
$\vb{s}=\vb{r}_2-\vb{r}_1$ とし、棒からの力を $\vb{R}_2=-\vb{R}_1=\lambda\vb{s}$ とすると、

$$
\vb{R}_1\cdot\delta\vb{r}_1+\vb{R}_2\cdot\delta\vb{r}_2
=\lambda\vb{s}\cdot(\delta\vb{r}_2-\delta\vb{r}_1)
=\tfrac{\lambda}{2}\delta(\vb{s}\cdot\vb{s})=0
$$

となる。最後の等号で棒の長さが一定であることを使った。
滑らずに転がる円板の接地点に働く静止摩擦力も、接地点の仮想変位が零なので仮想仕事をしない。
一方、粗い面上を滑る粒子では、動摩擦力は面に沿う変位に仕事をするので、印加力の側に含める。

理想拘束の拘束力がどのような形をとるかは、次の線形代数の補題で決まる。

:::lemma[核に直交するベクトル]{#lem-am-annihilator}
$A$ を $k\times m$ 実行列とする。$w\in\R^m$ が、$A\xi=0$ を満たすすべての $\xi\in\R^m$ に対して $w\cdot\xi=0$ を満たすならば、$w=\transp{A}\lambda$ となる $\lambda\in\R^k$ が存在する。
$A$ の階数が $k$ ならば、この $\lambda$ は一意である。
:::

:::proof
$A$ の行ベクトルが張る部分空間を $W$、核を $K=\{\xi:A\xi=0\}$、$K$ に直交するベクトル全体を $K^\perp$ とする。
$A\xi=0$ は各行と $\xi$ の内積が零であることなので、$W\subset K^\perp$ である。
次元を比べる。$r=\rank A$ とすると、行階数と列階数の一致（:ref[thm-linalg-row-column-rank]）から $\dim W=r$、次元公式から $\dim K=m-r$ である。
$K$ の基底 $\xi_1,\ldots,\xi_{m-r}$ を行とする行列 $B$ を作れば $K^\perp=\Ker B$ であり、$\rank B=m-r$ なので再び次元公式から $\dim K^\perp=m-(m-r)=r$ となる。
したがって $W=K^\perp$ であり、仮定 $w\in K^\perp$ から $w$ は行の一次結合、すなわち $w=\transp{A}\lambda$ と書ける。
$r=k$ なら行は一次独立なので、$\transp{A}\lambda=\transp{A}\lambda'$ から $\lambda=\lambda'$ が従う。
:::

:::proposition[理想拘束力と Lagrange 乗数]{#prop-am-ideal-multiplier}
独立なホロノミック拘束 :eqref[eq-am-holonomic] の下で、拘束力 $\vb{R}_a$ が理想的であるための必要十分条件は、各時刻に実数 $\lambda_1,\ldots,\lambda_k$ があって

$$
\vb{R}_a=\sum_{\alpha=1}^k\lambda_\alpha\nabla_af_\alpha,\qquad a=1,\ldots,N
\label{eq-am-constraint-force-multiplier}
$$

となることである。$\lambda_\alpha$ は $\vb{R}_a$ から一意に定まり、:term[Lagrange 乗数]{en="Lagrange multiplier" yomi="らぐらんじゅじょうすう"} と呼ばれる。
独立な Pfaff 型拘束の場合も、仮想変位を :eqref[eq-am-pfaffian-virtual] で定めれば、$\nabla_af_\alpha$ を係数 $a_{\alpha j}$ のうち第 $a$ 質点の三成分に置き換えて同じことが成り立つ。
:::

:::proof
:eqref[eq-am-constraint-force-multiplier] が成り立てば、:eqref[eq-am-virtual-displacement] により
$\sum_a\vb{R}_a\cdot\delta\vb{r}_a=\sum_\alpha\lambda_\alpha\sum_a\nabla_af_\alpha\cdot\delta\vb{r}_a=0$ である。
逆に理想拘束なら、拘束力を並べた $3N$ 成分のベクトルは $D_Xf$ の核のすべてのベクトルと直交する（:ref[prop-am-virtual-space]）。
:ref[lem-am-annihilator] を $A=D_Xf$ に適用すれば、それは行 $\nabla f_\alpha$ の一次結合として一意に書ける。
Pfaff 型拘束では $A=(a_{\alpha j})$ とすればよい。
:::

本サイトでは、拘束力を :eqref[eq-am-constraint-force-multiplier] のように $+\lambda_\alpha\nabla_af_\alpha$ の形に書く規約を一貫して用いる。
例えば振り子の $f=\tfrac12(x^2+y^2-\ell^2)$ では $\nabla f=\vb{r}$ なので $\vb{R}=\lambda\vb{r}$ である。棒が質点を支点の側へ引く張力のとき $\lambda<0$ で、張力の大きさは $\abs{\lambda}\ell$ である。

理想拘束が実際の運動でする仕事については、次のことがいえる。

:::proposition[定常な理想拘束は仕事をしない]{#prop-am-ideal-no-work}
定常なホロノミック拘束の下で、一般化座標の定義式 $\vb{r}_a(q)$ が時刻に陽に依存しないとする。拘束が理想的ならば、実際の運動に沿って $\sum_{a=1}^N\vb{R}_a\cdot\vb{v}_a=0$ である。
したがって、印加力が :ref[prop-am-system-energy] の仮定を満たせば、拘束があっても力学的エネルギーは保存する。
:::

:::proof
$\partial\vb{r}_a/\partial t=0$ なので、:eqref[eq-am-coordinate-velocity] から $\vb{v}_a=\sum_i(\partial\vb{r}_a/\partial q^i)\dot q^i$ である。
これは $\delta q^i=\dot q^i$ とした仮想変位に等しいから、:eqref[eq-am-ideal-work] により拘束力の仕事率は零である。
エネルギーについては、:ref[prop-am-system-energy] の証明で $\dot T=\sum_a(\vb{F}_a+\vb{R}_a)\cdot\vb{v}_a$ の第二項が消えることを使えばよい。
:::

非定常拘束ではこの結論は成り立たない。上下する板の抗力は実際の運動に対して仕事率 $R\dot h$ をもつ。

印加力の仮想仕事を一般化座標で表してみよう。式 :eqref[eq-am-virtual-displacement] を代入して和を並べ替えると、

$$
\delta W=\sum_{a=1}^N\vb{F}_a\cdot\delta\vb{r}_a
=\sum_{i=1}^n Q_i\,\delta q^i,\qquad
Q_i=\sum_{a=1}^N\vb{F}_a\cdot\pdv{\vb{r}_a}{q^i}.
\label{eq-am-generalized-force}
$$

係数 $Q_i$ を $q^i$ に対応する :term[一般化力]{en="generalized force" yomi="いっぱんかりょく"} という。
$Q_i\delta q^i$ の単位は仕事の単位 $\mathrm{J}$ である。$q^i$ が長さなら $Q_i$ は力、回転角なら力のモーメントの単位をもつ。
同じ定義を拘束力に適用すると、$\partial\vb{r}_a/\partial q^i$ 自身が仮想変位なので、理想拘束力の一般化力はすべて零である。
これが一般化座標を使うと拘束力が方程式から消える理由である。

印加力がポテンシャル $V(\vb{r}_1,\ldots,\vb{r}_N,t)$ から $\vb{F}_a=-\nabla_aV$ と導かれる場合、連鎖律により

$$
Q_i=-\sum_{a=1}^N\nabla_aV\cdot\pdv{\vb{r}_a}{q^i}=-\pdv{}{q^i}V(\vb{r}_1(q,t),\ldots,\vb{r}_N(q,t),t)
$$

である。

### 静力学のつり合い

以下この小節では、拘束は定常で独立なホロノミック拘束とし、印加力は配位だけの関数 $\vb{F}_a(X)$ とする。

:::definition[平衡配位]{#def-am-equilibrium}
許される配位 $X_0$ について、静止状態 $X(t)\equiv X_0$ が、ある理想拘束力 $\vb{R}_a$ とともに Newton 方程式 $m_a\ddot{\vb{r}}_a=\vb{F}_a+\vb{R}_a$ と拘束条件を満たすとき、$X_0$ を :term[平衡配位]{en="equilibrium configuration" yomi="へいこうはいい"} といい、系は $X_0$ でつり合うという。
:::

静止状態では加速度が零なので、この条件は $\vb{F}_a(X_0)+\vb{R}_a=0$ となる理想拘束力が存在することと同じである。

:::principle[仮想仕事の原理]{#prin-am-virtual-work}
上の仮定の下で、許される配位 $X_0$ について次の三条件は同値である。

1. $X_0$ は平衡配位である。
2. $X_0$ におけるすべての仮想変位に対して $\sum_{a=1}^N\vb{F}_a(X_0)\cdot\delta\vb{r}_a=0$ である。
3. 実数 $\lambda_1,\ldots,\lambda_k$ があって

$$
\vb{F}_a(X_0)+\sum_{\alpha=1}^k\lambda_\alpha\nabla_af_\alpha(X_0)=0,\qquad a=1,\ldots,N.
\label{eq-am-static-multipliers}
$$

正則な一般化座標で書けば、2 は一般化力がすべて零、$Q_i(q_0)=0$ $(i=1,\ldots,n)$ であることと同値である。
:::

:::proof
1 ならば 2：$\vb{F}_a=-\vb{R}_a$ と :eqref[eq-am-ideal-work] による。

2 ならば 3：印加力を並べた $3N$ 成分のベクトルは $D_Xf$ の核と直交する（:ref[prop-am-virtual-space]）。:ref[lem-am-annihilator] によりそれは $\transp{(D_Xf)}\mu$ と書けるので、$\lambda=-\mu$ とすればよい。

3 ならば 1：$\vb{R}_a=\sum_\alpha\lambda_\alpha\nabla_af_\alpha$ とおけば、:ref[prop-am-ideal-multiplier] によりこれは理想拘束力であり、$\vb{F}_a(X_0)+\vb{R}_a=0$ を満たす。したがって $X(t)\equiv X_0$ は Newton 方程式と拘束条件を満たす。

最後の主張は、:ref[prop-am-virtual-space] により仮想変位が $\delta q\in\R^n$ を任意に動かして得られ、$\delta W=\sum_iQ_i\delta q^i$ がすべての $\delta q$ で零であることが $Q=0$ と同値であることによる。
:::

条件 3 は、未知数 $X_0$ の $3N$ 成分と $\lambda$ の $k$ 成分に対し、:eqref[eq-am-static-multipliers] の $3N$ 本と拘束 $f_\alpha=0$ の $k$ 本を連立させるものである。
解けば平衡配位と同時に、反力 $\vb{R}_a=\sum_\alpha\lambda_\alpha\nabla_af_\alpha$ も求まる。
一方、条件 2 を一般化座標で使えば、反力を一切計算せずに $n$ 本の方程式 $Q_i=0$ で平衡配位が決まる。

:::remark[静止を保つこと]
定義では、静止状態が運動方程式の**一つの**解であることを要求した。初めに静止していた系が実際にその状態にとどまるには、さらに拘束付きの運動の一意性が必要である。
理想拘束の下では運動が $(q,\dot q)$ の初期値から一意に決まることを [Lagrange 形式](/physics/analytical-mechanics/lagrangian/) で示す。これによって、平衡配位から初速度零で出発した系は静止し続ける。
摩擦のように理想的でない拘束では仮想仕事の原理はこの形では成り立たない。粗い斜面上に静止した物体には、印加力の仮想仕事が零でなくても静止摩擦力がつり合いを保つ。
:::

印加力が位置エネルギー $V(X)$ から $\vb{F}_a=-\nabla_a V$ と書ける場合は、
$\delta W=-\delta V$ であり、平衡配位は拘束面上に制限した $V$ の停留点である。
:eqref[eq-am-static-multipliers] は $V-\sum_\alpha\lambda_\alpha f_\alpha$ の停留条件であり、多変数関数の条件付き極値に対する Lagrange の乗数法そのものである。
停留点には不安定なつり合いもある。平衡の安定性は [微小振動](/physics/analytical-mechanics/small-oscillations/) で扱う。

:::example[水平力を受ける振り子]{#exm-am-pendulum-horizontal}
:ref[exm-am-pendulum-constraint] のおもりに、重力 $(0,-mg)$ と一定の水平力 $(P,0)$ を加える。
棒の反力は接線に直交するので、印加力の仮想仕事だけを計算すると

$$
\delta W=\ell(P\cos\theta-mg\sin\theta)\delta\theta.
$$

したがって平衡角は $P\cos\theta=mg\sin\theta$ を満たす。
例えば $P=0$ なら下向きと上向きの二つの平衡がある。$V=-mg\ell\cos\theta$ は前者で極小、後者で極大なので、両者の安定性は異なる。

棒の張力を求めるには乗数を使う。$f=\tfrac12(x^2+y^2-\ell^2)$ に対し、:eqref[eq-am-static-multipliers] は

$$
P+\lambda x=0,\qquad -mg+\lambda y=0,\qquad x^2+y^2=\ell^2
$$

となる。$y<0$ の解は $\lambda=-\sqrt{P^2+m^2g^2}/\ell$、$(x,y)=\ell(P,-mg)/\sqrt{P^2+m^2g^2}$ である。
棒の力は $\vb{R}=\lambda\vb{r}$ で支点を向き、その大きさは $\sqrt{P^2+m^2g^2}$ である。これは印加力の合力の大きさに等しい。
:::

:::example[てこ]
支点のまわりに回転できる軽い棒の両側、支点からの距離 $\ell_1,\ell_2$ の点に質量 $m_1,m_2$ のおもりを置く。
棒が水平から角 $\theta$ だけ傾くと、二つのおもりの高さは $\ell_1\sin\theta$ と $-\ell_2\sin\theta$ である。$\theta=0$ における仮想仕事は

$$
\delta W=-m_1g\ell_1\delta\theta+m_2g\ell_2\delta\theta
$$

なので、水平につり合う条件は $m_1\ell_1=m_2\ell_2$ である。これは Archimedes のてこの原理である。
支点の反力は仮想仕事に現れない。
:::

特定の拘束力だけを知りたい場合には、**その拘束を外し、拘束力を未知の印加力とみなす**方法が有効である。
外した拘束の分だけ仮想変位が増え、その新しい方向の仮想仕事の条件から、外した拘束力が求まる。残った拘束は依然として理想的なので、その反力は現れない。

```tikz caption="滑らかな壁と床に立てかけ、下端を糸で壁につないだはしご。"
\begin{tikzpicture}[>=Stealth]
\draw[thick] (0,0) -- (0,3);
\draw[thick] (0,0) -- (4,0);
\coordinate (A) at (2.6,0);
\coordinate (B) at (0,2.2);
\draw[very thick] (A) -- (B);
\draw[dashed] (A) -- (0,0);
\draw[->] (2.6,0.0) ++(180:0.6) arc (180:140:0.6);
\node at (1.75,0.25) {$\alpha$};
\draw[->] (1.3,1.1) -- (1.3,0.3) node[right] {$mg$};
\draw[->,thick] (A) -- ++(-0.8,0) node[below] {$S$};
\end{tikzpicture}
```

:::example[糸でつないだはしご]{#exm-am-ladder}
長さ $\ell$、質量 $m$ の一様なはしごを、滑らかな鉛直の壁と滑らかな水平の床に立てかけ、下端を水平な糸で壁の根元につなぐ。床とはしごのなす角を $\alpha$ とする。糸の張力 $S$ を求めよう。
糸を外すと、はしごは壁と床に接したまま $\alpha$ を変えられる。壁を $x=0$ とすれば、下端の位置は $x=\ell\cos\alpha$、重心の高さは $(\ell/2)\sin\alpha$ である。
糸の力は下端に $-S$ の $x$ 成分として働くので、$\delta\alpha$ に対する印加力の仮想仕事は

$$
\delta W=(-S)\,\delta(\ell\cos\alpha)-mg\,\delta\!\left(\tfrac{\ell}{2}\sin\alpha\right)
=\ell\left(S\sin\alpha-\tfrac12mg\cos\alpha\right)\delta\alpha .
$$

つり合いでは係数が零なので $S=\tfrac12mg\cot\alpha$ である。
壁と床の抗力は理想拘束力なので計算に現れない。Newton 方程式で確かめると、水平方向のつり合いから壁の抗力は $S$ に等しく、下端まわりのモーメントのつり合い $S\ell\sin\alpha=mg(\ell/2)\cos\alpha$ から同じ結果を得る。
:::

## 演習問題

::::exercise{level=basic}
水平面内の二質点を長さ $\ell>0$ の軽い剛体棒で結ぶ。自由度を数え、中点 $(X,Y)$ を使って両端の位置を表せ。
:::solution
四つの位置成分に対し、二点間の距離を指定する拘束は一つである。
拘束 $f=\tfrac12(\abs{\vb{r}_1-\vb{r}_2}^2-\ell^2)$ の勾配は $(\vb{r}_1-\vb{r}_2,\vb{r}_2-\vb{r}_1)$ で、距離が零でないので消えない。よって :ref[prop-am-constraint-count] により自由度は $4-1=3$ となる。
残る一変数を棒の方位角 $\theta$ とし、$\vb{e}=(\cos\theta,\sin\theta)$ とおくと、

$$
\vb{r}_1=(X,Y)-\tfrac{\ell}{2}\vb{e},\qquad
\vb{r}_2=(X,Y)+\tfrac{\ell}{2}\vb{e}.
$$

差の長さは $\ell$ であり、逆に両端の位置から中点と棒の向きが一意に決まる（角は $2\pi$ を法とする）。
:::
::::

::::exercise{level=basic}
角度 $\alpha$ だけ右上がりの滑らかな固定斜面上に質量 $m$ の物体を置く。
$0<\alpha<\pi/2$ として、水平右向きにいくらの力を加えれば静止するか。抗力の大きさも求めよ。
:::solution
斜面上向きの単位ベクトルは $(\cos\alpha,\sin\alpha)$ である。斜面に沿う任意の変位を $\delta s$ とすると、

$$
\delta W=(P\cos\alpha-mg\sin\alpha)\delta s.
$$

:ref[prin-am-virtual-work] により係数が零であるから $P=mg\tan\alpha$ を得る。
抗力 $\vb{R}$ は印加力 $(P,-mg)$ とつり合うので $\vb{R}=(-P,mg)$ である。面の外向き法線 $(-\sin\alpha,\cos\alpha)$ に射影すると、

$$
R=P\sin\alpha+mg\cos\alpha=mg/\cos\alpha>0
$$

であり、接触を保てる。接線成分は $-P\cos\alpha+mg\sin\alpha=0$ で、抗力は確かに面に垂直である。
:::
::::

::::exercise{level=basic}
支点が固定された振り子の棒の長さを、装置で $\ell(t)>0$ に変える。$\vb{r}=\ell(t)(\sin\theta,-\cos\theta)$ として、
棒に沿う力の仮想仕事と、実際の運動で単位時間にする仕事を求めよ。
:::solution
$\vb{e}_r=(\sin\theta,-\cos\theta)$、$\vb{e}_\theta=(\cos\theta,\sin\theta)$ とおく。
時刻を固定すると $\ell$ は変わらないから $\delta\vb{r}=\ell\vb{e}_\theta\delta\theta$ である。
一方で実速度は

$$
\vb{v}=\dot\ell\vb{e}_r+\ell\dot\theta\vb{e}_\theta
$$

となる。力を $R\vb{e}_r$ と書くと、$\vb{e}_r\cdot\vb{e}_\theta=0$ よりその仮想仕事は零だが、実際の仕事率は $R\dot\ell$ である。
拘束は理想的だが非定常なので、:ref[prop-am-ideal-no-work] は適用できない。棒を伸縮する装置と質点の間でエネルギーを受け渡せることがわかる。
:::
::::

::::exercise{level=basic}
次の速度拘束が積分可能かどうかを判定せよ。積分可能なら位置の条件を求めよ。

1. $y\,\dot x+x\,\dot y+2z\,\dot z=0$
2. $\dot z-y\,\dot x=0$
:::solution
1. 左辺は $\dd(xy+z^2)/\dd t$ に等しい。したがって $xy+z^2$ が一定という位置の条件に積分でき、ホロノミックである。
実際 $\vb{a}=(y,x,2z)$ は勾配なので $\nabla\times\vb{a}=0$ であり、:eqref[eq-am-pfaff-test] も満たす。

2. $\vb{a}=(-y,0,1)$ とすると

$$
\nabla\times\vb{a}=\left(\pdv{a_3}{y}-\pdv{a_2}{z},\ \pdv{a_1}{z}-\pdv{a_3}{x},\ \pdv{a_2}{x}-\pdv{a_1}{y}\right)=(0,0,1)
$$

であり、$\vb{a}\cdot(\nabla\times\vb{a})=1\ne0$ である。:ref[prop-am-pfaff-test] により、どの点の近傍でも積分可能でなく、非ホロノミックである。
:::
::::

::::exercise{level=basic}
水平面内で一定の角速度 $\omega$ で回転する細い直線の棒に、質量 $m$ の小環を通す。回転軸からの距離 $r$ を一般化座標とし、:eqref[eq-am-kinetic-generalized] の $M$、$b$、$T_0$ を求めよ。
:::solution
小環の位置は $\vb{r}=r(\cos\omega t,\sin\omega t)$ である。

$$
\pdv{\vb{r}}{r}=(\cos\omega t,\sin\omega t),\qquad
\pdv{\vb{r}}{t}=r\omega(-\sin\omega t,\cos\omega t).
$$

二つは直交するので $b=0$ であり、$M=m$、$T_0=\tfrac12m\omega^2r^2$ である。よって $T=\tfrac12m(\dot r^2+\omega^2r^2)$ となる。
一次の項が消えても、拘束が非定常であることは $T_0\ne0$ に現れる。
:::
::::

::::exercise{level=advanced}
内力が作用反作用の法則の弱い形に従うが強い形には従わない二質点系を考える。$\vb{f}_{12}=-\vb{f}_{21}=\vb{F}$ とし、外力はないとする。全角運動量の時間変化を求め、$\vb{F}$ が $\vb{r}_1-\vb{r}_2$ に垂直な場合に $\vb{L}$ が保存しないことを示せ。
:::solution
:ref[prop-am-momentum-theorems] の証明と同じ計算で

$$
\dot{\vb{L}}=\vb{r}_1\times\vb{f}_{12}+\vb{r}_2\times\vb{f}_{21}=(\vb{r}_1-\vb{r}_2)\times\vb{F}
$$

である。$\vb{F}\perp(\vb{r}_1-\vb{r}_2)$ で $\vb{F}\ne0$、$\vb{r}_1\ne\vb{r}_2$ なら、外積の大きさは $\abs{\vb{r}_1-\vb{r}_2}\abs{\vb{F}}>0$ である。よって $\dot{\vb{L}}\ne0$ である。
一方、$\dot{\vb{P}}=\vb{f}_{12}+\vb{f}_{21}=0$ なので全運動量は保存する。角運動量の保存には強い形が必要である。
:::
::::

::::exercise{level=advanced}
半径 $a$ の鉛直円板を、初めは $\psi=0$ として距離 $s>0$ だけ転がす。その場で左に直角を向き、同じ距離を後退し、最後に方位を元に戻す。
中心の変位と二つの角の変化を求めよ。また、方位を常に零に固定した場合に積分できる条件を二つ挙げよ。
:::solution
最初の転がりで $(\Delta x,\Delta y)=(s,0)$、自転角の増加は $s/a$ である。
方位変更では中心は動かない。後退の方位は $\pi/2$ なので変位は $(0,-s)$、自転角の増加は $-s/a$ である。
よって中心の正味の変位は $(s,-s)$ で、$\psi,\phi$ の正味の変化はともに零となる。本文の式で $\beta=\pi/2$ とした場合に当たる。
これに対して方位を固定すると :eqref[eq-am-rolling-disk] は $\dot y=0,\ \dot x=a\dot\phi$ となる。
したがって $y$ および $x-a\phi$ が一定である。
:::
::::

::::exercise{level=advanced}
滑らかな放物線 $y=x^2/(2a)$（$a>0$、上向きが $y$）に質量 $m$ の小環を通す。
重力と一定の水平力 $P$ の下で平衡位置を求め、安定性と針金からの反力を調べよ。反力は Lagrange 乗数を用いて求めよ。
:::solution
独立な座標として $x$ を選ぶと、$\delta y=(x/a)\delta x$ である。
したがって $\delta W=(P-mgx/a)\delta x$ がすべての $\delta x$ で零となる位置は

$$
x_*=\frac{aP}{mg},\qquad y_*=\frac{aP^2}{2m^2g^2}.
$$

拘束上の位置エネルギーは $V(x)=mgx^2/(2a)-Px$ であり、$V''=mg/a>0$ なので極小である。よって平衡は安定である。
反力は、拘束 $f=y-x^2/(2a)$ に対して :eqref[eq-am-static-multipliers] を書くと

$$
P-\lambda\frac{x}{a}=0,\qquad -mg+\lambda=0
$$

となる。よって $\lambda=mg$、$x=aP/(mg)$ であり、反力は $\vb{R}=\lambda\nabla f=mg(-x_*/a,1)=(-P,mg)$ である。
接ベクトル $(1,x_*/a)$ と内積をとれば零になり、反力は確かに針金に垂直である。
:::
::::

::::exercise{level=advanced}
長さ $\ell$、質量 $m$ の一様な二本の棒の上端を滑らかな蝶番でつなぎ、滑らかな床に A 字形に立てる。二本の下端を水平な糸で結び、各棒が鉛直となす角を $\alpha$ とする。糸の張力を仮想仕事の原理で求めよ。
:::solution
糸を外して張力 $S$ を印加力とみなす。対称な変形 $\delta\alpha$ を考えると、下端の間隔は $2\ell\sin\alpha$、各棒の重心の高さは $(\ell/2)\cos\alpha$ である。
張力は各下端を内側に引くので、間隔の増加 $\delta(2\ell\sin\alpha)$ に対して仕事 $-S\,\delta(2\ell\sin\alpha)$ をする。したがって

$$
\delta W=-S\cdot2\ell\cos\alpha\,\delta\alpha-2mg\cdot\left(-\tfrac{\ell}{2}\sin\alpha\right)\delta\alpha
=\ell(mg\sin\alpha-2S\cos\alpha)\delta\alpha .
$$

つり合いでは係数が零なので $S=\tfrac12mg\tan\alpha$ である。
床の抗力と蝶番の力は理想拘束力なので現れない。
Newton 方程式による確認：右の棒について、床の抗力は対称性と鉛直方向のつり合いから $mg$、蝶番の力は対称性から水平である。上端まわりのモーメントのつり合いは $mg\ell\sin\alpha-S\ell\cos\alpha-mg(\ell/2)\sin\alpha=0$ となり、同じ結果を与える。
:::
::::

::::exercise{level=advanced}
半径 $b$ の滑らかな固定球の頂上から、質量 $m$ の小物体が初速度ほぼ零で滑り始める。鉛直上向きから測った角を $\theta$ とし、物体が球面を離れる角を求めよ。片側拘束であることに注意せよ。
:::solution
拘束は $g(\vb{r})=\abs{\vb{r}}-b\ge0$ である。接触している間は球面上の運動で、抗力 $\vb{R}=\lambda\grad g=\lambda\vb{r}/b$ は面に垂直なので仕事をしない（:ref[prop-am-ideal-no-work]）。
エネルギー保存より $\tfrac12mv^2=mgb(1-\cos\theta)$ である。動径方向の運動方程式は

$$
-m\frac{v^2}{b}=-mg\cos\theta+\lambda
$$

であり、$\lambda=mg\cos\theta-mv^2/b=mg(3\cos\theta-2)$ を得る。
片側拘束の抗力は外向き、すなわち $\lambda\ge0$ でなければならない。$\lambda$ は $\theta$ とともに減少し、$\cos\theta=2/3$（$\theta\simeq48.2^\circ$）で零になる。
これより先では等式拘束の解は $\lambda<0$ を要求するので、物体はこの角で球面を離れ、以後は放物運動をする。
:::
::::

::::exercise{level=advanced}
空間の一質点に二つの拘束 $f_1=\tfrac12(x^2+y^2+z^2-1)=0$、$f_2=z-c=0$ を課す（$c$ は定数）。:ref[prop-am-constraint-count] の仮定が成り立つ $c$ の範囲と、そのときの自由度を求めよ。仮定が破れる $c$ では何が起こるか。
:::solution
勾配は $\nabla f_1=(x,y,z)$、$\nabla f_2=(0,0,1)$ である。両者が一次従属になるのは $x=y=0$ のときに限る。
拘束を満たす点で $x=y=0$ なら $z=c$、$z^2=1$ なので $c=\pm1$ である。
したがって $\abs{c}<1$ では階数は $2$ で、自由度は $3-2=1$（配位空間は円 $x^2+y^2=1-c^2$）である。$\abs{c}>1$ では許される配位がない。
$c=\pm1$ では平面が球に接し、許される配位は一点 $(0,0,\pm1)$ だけで自由度は零である。この点で階数は $1$ に落ちており、公式 $n=3N-k$ は $1$ を与えて実際と合わない。
:::
::::

## 確認テスト

::::quiz
平面振り子の支点が水平方向に指定された運動をしている。実験室の座標でみた棒の長さの条件はどれか。

- [ ] 定常な非ホロノミック拘束
- [x] 非定常なホロノミック拘束
- [ ] 非定常な非ホロノミック拘束

:::explanation
$(x-b(t))^2+y^2=\ell^2$ と位置と時間の等式で書けるのでホロノミックである。支点の位置 $b(t)$ に応じて許される円が移動し、条件が時刻に陽に依存するので非定常である。
:::
::::

::::quiz{answer="5"}
三次元空間の二質点を長さが正で一定の棒で結ぶ。他に拘束がなければ自由度はいくつか。

:::explanation
六つの位置成分に対し、距離の条件は一つである。距離が正なのでその勾配は零でなく、:ref[prop-am-constraint-count] より自由度は $6-1=5$ である（重心の並進 $3$ と棒の向き $2$）。
:::
::::

::::quiz
理想拘束について常に成り立つものを選べ。

- [ ] 拘束力は実際の運動で仕事をしない
- [ ] 各質点ごとに拘束力の仮想仕事が零である
- [x] 全質点について足した拘束力の仮想仕事が、許されるどの仮想変位でも零である

:::explanation
:ref[def-am-ideal-constraint] では和をとる。二質点を結ぶ棒では各端の仮想仕事は零でないが、両端の和が相殺する。
また、実際の運動で仕事をしないのは定常な場合に限られ（:ref[prop-am-ideal-no-work]）、上下する板では抗力が実際の運動に対して仕事をする。
:::
::::

::::quiz
三変数の速度拘束 $\vb{a}\cdot\dot{\vb{x}}=0$ について、$\vb{a}\cdot(\nabla\times\vb{a})$ がある点で零でない。正しいものを選べ。

- [x] その点の近傍で、この拘束は位置の条件に積分できない
- [ ] 適当な積分因子を掛ければ必ず位置の条件に積分できる
- [ ] この判定からは何もいえない

:::explanation
:ref[prop-am-pfaff-test] により、積分可能ならばその近傍で $\vb{a}\cdot(\nabla\times\vb{a})=0$ でなければならない。対偶をとれば、零でない点の近傍では積分可能でない。スケートの刃の拘束がその例である。
:::
::::

::::quiz{answer="0.5" tol="1e-6"}
固定支点の振り子に右向き $4.9\,\mathrm N$ の力を加えた。$m=1.0\,\mathrm{kg}$、$g=9.8\,\mathrm{m\,s^{-2}}$ として、平衡での $\tan\theta$ はいくらか。角は鉛直下向きから右へ測る。

:::explanation
仮想仕事 $\delta W=\ell(P\cos\theta-mg\sin\theta)\delta\theta$ が零となる条件 $P\cos\theta=mg\sin\theta$ から、$\tan\theta=P/(mg)=4.9/9.8=0.5$ となる。
:::
::::

::::quiz{answer="0.5" tol="1e-6"}
:ref[exm-am-ladder] で $\alpha=45^\circ$ のとき、糸の張力は $mg$ の何倍か。

:::explanation
$S=\tfrac12mg\cot\alpha$ で $\cot45^\circ=1$ なので $S/(mg)=0.5$ である。
:::
::::

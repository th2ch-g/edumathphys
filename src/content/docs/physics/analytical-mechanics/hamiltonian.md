---
title: Hamilton 形式
description: 速度と運動量を交換する条件を調べ、位相空間での変分原理と体積保存を導く。
notation:
  - symbol: 'W'
    meaning: 速度に関する Lagrangian の Hessian
  - symbol: 'H'
    meaning: Hamiltonian
  - symbol: 'z'
    meaning: 位相空間の点
  - symbol: '\Phi_{t,s}'
    meaning: 時刻 s の初期値を時刻 t に送る写像
---

## この章の目標

- 速度から運動量への変換が逆に解ける条件を調べる。
- Hamiltonian を構成し、その偏微分から時間発展を求める。
- 位相空間で独立に変分する変数と端点条件を区別する。
- Liouville の定理を写像の行列式と密度の式で示す。

## Legendre 変換

位置と速度を与えれば運動が決まる。
速度の代わりに、運動の方向に作用がどれだけ変わるかを表す共役運動量を使うと、運動方程式が一対の一階方程式になる。
単に速度の記号を書き換えるだけでは、この対称な形は得られない。

### 接線の傾きを独立変数にする

開区間上の $C^2$ 級の関数で $f''(v)>0$ とする。
傾き $p$ を $f'$ の値域から選ぶ。$f'(v)=p$ の解は一意であり、$pv-f(v)$ の微分 $p-f'(v)$ はその点の左で正、右で負となる。従ってその点で最大になる。
接線 $y=px-f^*(p)$ の切片を用いて $f^*(p)=pv-f(v)$ とおくと、$\dd f^*/\dd p=v+(p-f')\dd v/\dd p=v$。
元の変数が傾きとして戻る。
この操作を :term[Legendre 変換]{en="Legendre transform" yomi="るじゃんどるへんかん"} という。
例えば $f=mv^2/2$ なら $v=p/m$ で $f^*=p^2/(2m)$ である。

### 力学に必要な正則性

$q,t$ を固定し、:ref[def-am-lagrangian] の定義で

$$
p_i=L_{v^i}(q,v,t),\qquad W_{ij}=\frac{\partial^2 L}{\partial v^i\partial v^j}
$$

とおく。
$C^2$ 級の $L$ で $\det W\ne0$ なら :ref[thm-anl-inverse-function] により近傍内で $v=v(q,p,t)$ が定まる。
このような $L$ を正則という。
力学の局所的な変換に必要なのは正則性である。
大域的な最大値という解釈には凸性と最大点の存在が別に必要になる。

:::definition[Hamiltonian]{#def-am-hamiltonian}
速度を運動量から求められる領域で

$$
H(q,p,t)=\sum_{i=1}^np_i v^i(q,p,t)-L(q,v(q,p,t),t)
\label{eq-am-legendre}
$$

を :term[Hamiltonian]{en="Hamiltonian" yomi="はみるとにあん"}（Hamilton 関数）という。
:ref[def-am-energy-function] を $(q,p)$ で表示したもので、単位は $\mathrm J$ である。
:::

:::example[速度に一次の項]
$L=m v^2/2+A(x,t)v-U(x,t)$ で $A$ は運動量と同じ単位をもつ。
$p=mv+A$、従って

$$
H=(p-A)^2/(2m)+U .
$$

$A=2\,\mathrm{kg\,m/s}$、$p=5\,\mathrm{kg\,m/s}$、$m=3\,\mathrm{kg}$ なら速さは $1\,\mathrm{m/s}$、
$H-U=1.5\,\mathrm J$。
共役運動量 $p$ と力学的運動量 $mv$ の区別は :ref[exm-am-electromagnetic-lagrangian] の荷電粒子にも現れる。
:::

:::example[逆が存在しないとき]
$L=m(v^1-v^2)^2/2$ なら $p_1=m(v^1-v^2)=-p_2$。
和 $v^1+v^2$ はどの運動量からも復元できず、許される運動量は $p_1+p_2=0$ の面にある。
ここで速度を自由に二つとも解く手順は使えない。
このような特異な系で拘束を保ちながら時間発展を構成する一般論は :cite[yamamoto-nakamura-am] に譲る。
:::

### 傾きから曲率まで復元する

速度と運動量の対応を $p=L_v$ として、$q,t$ を固定したまま $p_j$ で微分する。
連鎖律から

$$
\delta_{ij}=\sum_{k=1}^n W_{ik}\frac{\partial v^k}{\partial p_j}.
$$

従って $\partial v/\partial p=W^{-1}$。
他方、$H_p=v$ なので $H$ の運動量 Hessian は $W^{-1}$ となる。
元の $L$ の速度 Hessian が正定値なら、$H$ の運動量 Hessian も正定値である。厳密な凸性だけでは速度 Hessian の非退化性は従わない。例えば $L=v^4/4$ では $v=0$ で Hessian が消える。
逆に $H$ が与えられた場合は $v=H_p$ を解き、$L=\sum_i p_iv^i-H$ とすればよい。
往復すれば元に戻るのは、どちらの変換でも鎖を通じた微分が相殺するからである。

:::example[非対角な運動エネルギー]
二つの変位を m で、行列を kg の数値で表し、

$$
M=\begin{pmatrix}2&1\\1&2\end{pmatrix},\qquad L=\tfrac12 v^\top Mv-V(q)
$$

とする。
$p=Mv$ より $v=M^{-1}p$、$M^{-1}=\frac13\begin{pmatrix}2&-1\\-1&2\end{pmatrix}$ となる。
従って

$$
H=\tfrac12p^\top M^{-1}p+V=\frac{p_1^2-p_1p_2+p_2^2}{3}+V
$$

である（係数は上記 SI 表示）。
$p=(3,0)$ なら $v=(2,-1)$、運動エネルギーは $3\,\mathrm J$。
$p_2=0$ でも $\dot q^2=0$ とは限らない。
:::

:::example[正則な枝が二つある場合]
一つの速度を基準単位で無次元化し $L/E_0=v^3/3$ とする。
共役な無次元運動量は $p=v^2$ で、$v\ne0$ なら Hessian $2v$ は非零である。
しかし正の $p$ は $v=\pm\sqrt p$ の両方を持ち、$H/E_0=\pm2p^{3/2}/3$ に二つの枝がある。
どちらでも微分は選んだ $v$ に戻る。
局所的に正準方程式が書けることと、全速度域で一つの Hamiltonian にできることの違いがここに現れる。
:::

### 相対論的な自由粒子での速度の値域

速度の上限があっても運動量の上限があるとは限らない。質量 $m>0$ の粒子に対し、真空での光速を $c$ として、慣性系の時刻を用いた $L=-mc^2\sqrt{1-|\vb v|^2/c^2}$ を考える。$|\vb v|<c$ で

$$
\vb p=\frac{m\vb v}{\sqrt{1-|\vb v|^2/c^2}}.
$$

$\vb p$ と $\vb v$ の向きは同じである。大きさの式を二乗して解けば $\vb v=c^2\vb p/\sqrt{m^2c^4+c^2|\vb p|^2}$ となり、任意の有限な $\vb p$ が一つの速度を指定する。$\vb p\cdot\vb v-L=mc^2/\sqrt{1-|\vb v|^2/c^2}$ だから

$$
H=\sqrt{m^2c^4+c^2|\vb p|^2}.
$$

実際に $H_{\vb p}=c^2\vb p/H$ は元の速度である。$|\vb p|\ll mc$ では $H=mc^2+|\vb p|^2/(2m)+O(|\vb p|^4/(m^3c^2))$ となり、静止エネルギーを引けば Newton 力学の運動エネルギーが先頭に現れる。

### 一部の速度だけを変換する

二つの自由度で $\theta$ が循環し、$q$ が循環しない場合を考える。
$p_\theta=L_{\dot\theta}=c$ を $\dot\theta$ について解いて $R(q,\dot q;c)=L-c\dot\theta$ を作ったのが Routh の縮約である。
残りの速度をさらに Legendre 変換すると

$$
p_q=\frac{\partial R}{\partial\dot q}
=L_{\dot q}+(L_{\dot\theta}-c)\frac{\partial\dot\theta}{\partial\dot q}=L_{\dot q}.
$$

従って全 Hamiltonian を直接作って $p_\theta=c$ と置いた結果は、$R$ から $p_q\dot q-R$ を作ったものに等しい。
符号の違う二つの部分変換を混同しなくて済む。

:::example[回転と伸びの結合をもつ模型]
$L=A\dot q^2/2+B\dot q\dot\theta+C\dot\theta^2/2-V(q)$ とする。
係数は定数で $A,C>0,\ AC>B^2$。
$c=B\dot q+C\dot\theta$ を解いて代入すると

$$
R=\frac12\left(A-\frac{B^2}{C}\right)\dot q^2
+\frac{Bc}{C}\dot q-V(q)-\frac{c^2}{2C}.
$$

よって $p_q=(A-B^2/C)\dot q+Bc/C$ であり、

$$
H_c=\frac{(p_q-Bc/C)^2}{2(A-B^2/C)}+V(q)+\frac{c^2}{2C}.
$$

全質量行列 $\begin{pmatrix}A&B\\B&C\end{pmatrix}$ の逆行列から直接 $H$ を作っても同じになる。
$B=0$ なら単純に回転のエネルギー $c^2/(2C)$ を足すだけだが、結合があれば残る運動量の原点もずれる。
:::

### 速度を消せない系で何が失われるか

正則性は変数の個数についての条件である。速度から運動量への写像の階数が落ちると、相空間に任意の初期点を選べない。例えば補助変数 $y$ とその速度を含む表示で $L=m\dot x^2/2-U(x,y)$、$U_{yy}\ne0$ とする。$p_y=0$ はどんな速度でも満たす :term[一次拘束]{en="primary constraint" yomi="いちじこうそく"} である。$p_x=m\dot x$ だけを消去し、未定速度 $w$ を残した関数 $H_T=p_x^2/(2m)+U+w p_y$ を作る。$\dot y=w$ となる。$p_y=0$ の時間微分も零であるためには $\dot p_y=-U_y=0$ が必要で、$U_y=0$ は二次拘束となる。さらに微分して $0=U_{yx}p_x/m+U_{yy}w$ とすれば未定速度も決まる。拘束を最初から関数に代入して全偏微分を消してしまうと、この順序が見えなくなる。

:::example[支点に質量のない中間ばね]
座標を $x,y$、その数値を m で取る。式中の係数は $U$ が J となる SI 値で、$m=2\,\mathrm{kg}$、$U=2 x^{2} - 4 x y + 3 y^{2}$ とする。$y$ に慣性を付けない模型では、それに付随する力は随時釣り合わなければならない。$U_y=2 \left(- 2 x + 3 y\right)$ を零とおくと $y=\frac{2 x}{3}$。この条件を保つ速度と、独立な変数に戻した Hamiltonian は

$$
\dot y=\frac{p}{3},\qquad H_{\rm red}=\frac{3 p^{2} + 8 x^{2}}{12}.
$$

ここで $p=p_x$ である。縮約後は $\dot p=-\dd U_{\rm eff}/\dd x=- \frac{4 x}{3}$。元の式で $-U_x$ を先に計算してから枝の値を代入しても同じになる。$x=1\,\mathrm m,p=3\,\mathrm{kg\,m/s}$ の瞬間で $y=\frac{2}{3}\,\mathrm m,\dot y=1\,\mathrm{m/s}$、残る自由度は一つである。二本のばねが直列のとき、中間点の力の釣り合いを静的に取ってから運動する端だけを追う操作である。
:::

:::example[曲がった平衡面をもつ補助座標]
座標を $x,y$、その数値を m で取る。式中の係数は $U$ が J となる SI 値で、$m=1\,\mathrm{kg}$、$U=\frac{3 x^{4} - 4 x^{2} y + 2 y^{2}}{4}$ とする。$y$ に慣性を付けない模型では、それに付随する力は随時釣り合わなければならない。$U_y=- x^{2} + y$ を零とおくと $y=x^{2}$。この条件を保つ速度と、独立な変数に戻した Hamiltonian は

$$
\dot y=2 p x,\qquad H_{\rm red}=\frac{2 p^{2} + x^{4}}{4}.
$$

ここで $p=p_x$ である。縮約後は $\dot p=-\dd U_{\rm eff}/\dd x=- x^{3}$。元の式で $-U_x$ を先に計算してから枝の値を代入しても同じになる。$x=1\,\mathrm m,p=3\,\mathrm{kg\,m/s}$ の瞬間で $y=1\,\mathrm m,\dot y=6\,\mathrm{m/s}$、残る自由度は一つである。補助座標の速度は位置に比例して変わる。運動量零はその変数の速度零を意味しないことが特に見やすい。
:::

:::example[ばねに線形荷重を受ける中間点]
座標を $x,y$、その数値を m で取る。式中の係数は $U$ が J となる SI 値で、$m=3\,\mathrm{kg}$、$U=x^{2} - 2 x y + 4 y^{2} - 2 y$ とする。$y$ に慣性を付けない模型では、それに付随する力は随時釣り合わなければならない。$U_y=2 \left(- x + 4 y - 1\right)$ を零とおくと $y=\frac{x + 1}{4}$。この条件を保つ速度と、独立な変数に戻した Hamiltonian は

$$
\dot y=\frac{p}{12},\qquad H_{\rm red}=\frac{2 p^{2} + 9 x^{2} - 6 x - 3}{12}.
$$

ここで $p=p_x$ である。縮約後は $\dot p=-\dd U_{\rm eff}/\dd x=- \frac{3 x - 1}{2}$。元の式で $-U_x$ を先に計算してから枝の値を代入しても同じになる。$x=1\,\mathrm m,p=3\,\mathrm{kg\,m/s}$ の瞬間で $y=\frac{1}{2}\,\mathrm m,\dot y=\frac{1}{4}\,\mathrm{m/s}$、残る自由度は一つである。一定荷重によってばねの自然状態からのずれが生じる。定数のエネルギー原点と平衡位置の移動は異なる。
:::

### 同じ特異性でも未定のまま残る速度

二つの長さ $x,y$ を用いて $L=m(\dot x-\dot y)^2/2-V(x-y)$ と記した場合、$p_x+p_y=0$ であり、$x,y$ の同時移動は任意の時刻関数 $\epsilon(t)$ に対して作用を変えない。$r=x-y,s=y$ と変えれば $p_r=m\dot r,p_s=0$ となり、$H_T=p_r^2/(2m)+V(r)+w p_s$ である。今度は $\dot p_s=0$ が自動的に成り立つので $w(t)$ は決まらない。$\dot s=w$ を変えても $r,p_r$ は変化しない。先の補助座標の例では未定速度が拘束を保つ条件で決まったが、ここでは同じ物理状態に複数の座標表示を許している。$s=0$ を時間全体で選ぶと初めて $w=0$ が決まる。拘束の一般的な分類には [Poisson 括弧](/physics/analytical-mechanics/canonical-transformations/) が必要になる。

## Hamilton の正準方程式

### 何を固定して偏微分するか

:eqref[eq-am-legendre] の全微分をとると

$$
\dd H=\sum_{i=1}^n\left(v^i\dd p_i-L_{q^i}\dd q^i+
(p_i-L_{v^i})\dd v^i\right)-L_t\dd t .
$$

最後の速度の変化は係数が零である。
Lagrange 方程式の解では $\dot p_i=L_{q^i}$ なので、

$$
\dot q^i=\frac{\partial H}{\partial p_i},\qquad
\dot p_i=-\frac{\partial H}{\partial q^i},\qquad H_t=-L_t
\label{eq-am-hamilton-equations}
$$

を得る。
最初の二式を :term[正準方程式]{en="canonical equations" yomi="せいじゅんほうていしき"} という。
逆にこの二式を満たすと、第一式から速度と運動量の定義が戻り、第二式から Lagrange 方程式が戻る。
両形式は正則な領域で同値である。
運動に沿う $H$ の微分では $\sum_i(H_{q^i}H_{p_i}-H_{p_i}H_{q^i})=0$ だから $\dot H=H_t$ となる。

:::example[平面極座標の遠心力]
$H=p_r^2/(2m)+p_\phi^2/(2mr^2)+V(r)$ から $\dot r=p_r/m$、$\dot\phi=p_\phi/(mr^2)$、$\dot p_\phi=0$、
$\dot p_r=p_\phi^2/(mr^3)-V'$。
$r$ 微分で固定するのは $p_\phi$ である。
$\dot\phi$ を固定して微分すれば遠心項の符号を誤る。
:::

:::example[位置に依存する質量係数]
$L=M(q)\dot q^2/2-V(q)$、$M>0$ とすると $H=p^2/(2M)+V$。
従って $\dot p=p^2M'/(2M^2)-V'$。
$p=M\dot q$ の微分を代入して $M\ddot q+M'\dot q^2/2=-V'$ を得る。
運動エネルギーが位置に依存する場合に $\dot p=-V'$ と置くことはできない。
:::

### 磁場中で正準方程式を最後まで戻す

電荷 $e$ の粒子に対して、成分表示で

$$
H=\frac1{2m}\sum_{j=1}^3(p_j-e A_j)^2+e\Phi,\qquad
v_j=(p_j-e A_j)/m .
$$

位置の式は $\dot x_j=v_j$。
第二の正準方程式では $p$ を固定するので

$$
\dot p_i=e\sum_{j=1}^3v_j\partial_i A_j-e\partial_i\Phi .
$$

実際の加速度に直すには $\dot p_i=m\dot v_i+e\partial_tA_i+e\sum_jv_j\partial_jA_i$ を使う。
従って

$$
m\dot v_i=e(-\partial_i\Phi-\partial_tA_i)+
e\sum_{j=1}^3v_j(\partial_i A_j-\partial_j A_i)
$$

となり :eqref[eq-am-lorentz-force] を得る。
例えば $A_y=Bx$、他成分と $\Phi$ を零とすると $\dot p_y=0$ だが $m\dot v_y=-eBv_x$。
ここで $\dot p_y$ を $m\dot v_y$ と読み替えると磁場の効果を一成分失う。

### 保存量を使って先に一次の式を積分する

$q^a$ が循環する場合、$\dot p_a=-H_{q^a}=0$ だから $p_a=c$ と置き、残る組 $(q^i,p_i)$ だけを $H|_{p_a=c}$ で発展させられる。
除いた $q^a$ は最後に $\dot q^a=H_{p_a}|_{p_a=c}$ から復元する。
例えば円筒面上の自由粒子では $H=p_z^2/(2m)+p_\phi^2/(2ma^2)$ で両運動量が保存する。
軸方向は等速、方位は等角速度だから空間ではらせんを描く。
$p_\phi\ne0$ なら一回転で進む高さは $2\pi a^2p_z/|p_\phi|$ である。

### 位置に依存する質量と運動量の偏微分

線上に制限した粒子でも、長さ以外を座標に選べば有効質量は変化する。$L=M(q)v^2/2+A(q,t)v-V(q,t)$、$M>0$ とすると、$H=(p-A)^2/(2M)+V$。正準方程式の第二式では $p$ を固定して微分し、最後に $v=(p-A)/M$ を代入する。この順序から
$$
\dot p=\tfrac12 M'v^2+A_qv-V_q,\qquad
M\dot v=-\tfrac12M'v^2-V_q-A_t
\label{eq-am-coordinate-mass}
$$
となる。後式は $p=Mv+A$ を時間微分して代入したものである。一次元では空間だけに依存する $A$ は力に寄与しない。一方、磁場のある二次元では $A$ の二つの交差微分の差が残る。$M'v^2/2$ は座標の選択から生じた慣性項であり、位置エネルギーの微分に置き換えられない。

:::example[放物線の横座標で記述した小環]
SI の基準量で座標と係数を数値化し、$M=q^{2} + 1,A=0,V=5 q^{2}$ とする。初めに $p=Mv+A$ を解いて、$H$ の引数から速度を除く。この模型の Hamiltonian と正準方程式は

$$
H=\frac{p^{2} + 10 q^{4} + 10 q^{2}}{2 \left(q^{2} + 1\right)},\qquad \dot q=\frac{p}{q^{2} + 1},\qquad \dot p=\frac{q \left(p^{2} - 10 q^{4} - 20 q^{2} - 10\right)}{\left(q^{2} + 1\right)^{2}}.
$$

次に一階の式を合成して座標加速度を求める。$\dot v=(\partial_qv)\dot q+(\partial_pv)\dot p+\partial_tv$ には、$v$ に含まれるすべての引数の変化が入る。整理した結果は $\dot v=- \frac{q \left(p^{2} + 10 q^{4} + 20 q^{2} + 10\right)}{\left(q^{2} + 1\right)^{3}}$。$(q,p,t)=(1,4,0)$ の瞬間で $(\dot q,\dot p,\ddot q)=(2,-6,-7)$ となる。形状を $y=x^2/2$ とした場合、弧長速度は横速度に $\sqrt{1+x^2}$ を掛けたものになる。曲率による質量の変化はこの幾何に由来する。
:::

:::example[時間変化する速度一次の結合]
SI の基準量で座標と係数を数値化し、$M=2,A=q t,V=0$ とする。初めに $p=Mv+A$ を解いて、$H$ の引数から速度を除く。この模型の Hamiltonian と正準方程式は

$$
H=\frac{\left(p - q t\right)^{2}}{4},\qquad \dot q=\frac{p - q t}{2},\qquad \dot p=\frac{t \left(p - q t\right)}{2}.
$$

次に一階の式を合成して座標加速度を求める。$\dot v=(\partial_qv)\dot q+(\partial_pv)\dot p+\partial_tv$ には、$v$ に含まれるすべての引数の変化が入る。整理した結果は $\dot v=- \frac{q}{2}$。$(q,p,t)=(2,6,1)$ の瞬間で $(\dot q,\dot p,\ddot q)=(2,2,-1)$ となる。ここでは $A$ の空間微分が消えても時間微分が力として残る。$A$ を時間依存のゲージ変換で消すなら、同時にスカラーポテンシャルを変える必要がある。
:::

:::example[指数座標で表した自由粒子]
SI の基準量で座標と係数を数値化し、$M=e^{2 q},A=0,V=0$ とする。初めに $p=Mv+A$ を解いて、$H$ の引数から速度を除く。この模型の Hamiltonian と正準方程式は

$$
H=\frac{p^{2} e^{- 2 q}}{2},\qquad \dot q=p e^{- 2 q},\qquad \dot p=p^{2} e^{- 2 q}.
$$

次に一階の式を合成して座標加速度を求める。$\dot v=(\partial_qv)\dot q+(\partial_pv)\dot p+\partial_tv$ には、$v$ に含まれるすべての引数の変化が入る。整理した結果は $\dot v=- p^{2} e^{- 4 q}$。$(q,p,t)=(0,2,0)$ の瞬間で $(\dot q,\dot p,\ddot q)=(2,4,-4)$ となる。基準長 $1\,\mathrm m$ に対し実際の長さを $x=e^q$ 倍とすると $q$ は無次元である。$q$ の加速度は零でなくても $\ddot x=e^q(\ddot q+\dot q^2)$ は零になる。
:::

:::example[流れるベルトに沿う速度のずれ]
SI の基準量で座標と係数を数値化し、$M=q^{2} + 1,A=q^{2} + 1,V=- \frac{q^{2} + 1}{2}$ とする。初めに $p=Mv+A$ を解いて、$H$ の引数から速度を除く。この模型の Hamiltonian と正準方程式は

$$
H=\frac{p \left(p - 2 q^{2} - 2\right)}{2 \left(q^{2} + 1\right)},\qquad \dot q=\frac{p - q^{2} - 1}{q^{2} + 1},\qquad \dot p=\frac{p^{2} q}{\left(q^{2} + 1\right)^{2}}.
$$

次に一階の式を合成して座標加速度を求める。$\dot v=(\partial_qv)\dot q+(\partial_pv)\dot p+\partial_tv$ には、$v$ に含まれるすべての引数の変化が入る。整理した結果は $\dot v=- \frac{p q \left(p - 2 q^{2} - 2\right)}{\left(q^{2} + 1\right)^{3}}$。$(q,p,t)=(1,0,0)$ の瞬間で $(\dot q,\dot p,\ddot q)=(-1,0,0)$ となる。これは $L=M(q)(v+1)^2/2$ の形である。$p=0$ では瞬間の一般化速度は負であり、座標とともに動く観測系の速度が運動量に含まれる。
:::

## 相空間と相流

### 点で指定する初期値

$(q^1,\ldots,q^n,p_1,\ldots,p_n)$ の空間を :term[相空間]{en="phase space" yomi="そうくうかん"} という。
局所的には $2n$ 次元である。
座標を $q=q(Q)$ と取り直すとき、作用の境界項 $\sum_i p_i\dd q^i=\sum_j P_j\dd Q^j$ を保つには $P_j=\sum_i p_i\partial q^i/\partial Q^j$ とする。
従って運動量は配位空間上の余接ベクトルであり、この変換則に従う全共役変数の空間は余接束 $T^*\mathcal Q$ である。与えた $L$ から実際に得られる Hamiltonian の定義域は、その開部分や局所的な枝に限られることがある。
例えば平面振り子の角度は $2\pi$ 周期なので、その相空間は円筒である。

$z=(q,p)^\top$、$J=\begin{pmatrix}0&I_n\\-I_n&0\end{pmatrix}$ とすれば $\dot z=J\nabla_z H$。
局所的に解が一意に存在する時間域で、時刻 $s$ の初期点 $z$ を時刻 $t$ の点へ送る写像 $\Phi_{t,s}(z)$ を :term[相流]{en="phase flow" yomi="そうりゅう"} という。
一意性により $\Phi_{t,s}\circ\Phi_{s,r}=\Phi_{t,r}$ かつ逆は $\Phi_{s,t}$ である。
自律系なら $\Phi_{t,s}=\Phi_{t-s,0}$ と一つの時間差で書ける。
同じ自律系の相平面上で、異なる軌道が交わることはない。交点を時刻の原点に取り直せば、二つの解が同じ初期値を持ち、一意性により一致するからである。時間依存系の相平面への射影では、異なる時刻に同じ点を通ることがある。

### 位相図を読む

振動子 $H=p^2/(2m)+m\omega^2 q^2/2$ の流れは

$$
\begin{pmatrix}q(t)\\p(t)\end{pmatrix}=
\begin{pmatrix}\cos\omega t&\sin\omega t/(m\omega)\\
-m\omega\sin\omega t&\cos\omega t\end{pmatrix}
\begin{pmatrix}q(0)\\p(0)\end{pmatrix}.
$$

微分すれば正準方程式を満たす。
エネルギー等高線は楕円であり、$(q,p)=(A,0)$ からは負の $p$ へ進む。
エネルギーを $E$ とすれば半軸の積は $2E/\omega$、囲む面積は $2\pi E/\omega$ である。

振り子 $H=p^2/(2ml^2)+mgl(1-\cos q)$ では $E<2mgl$ の正のエネルギーは往復振動、$E>2mgl$ は回転に対応する。
境界の $E=2mgl$ は :term[セパラトリクス]{en="separatrix" yomi="せぱらとりくす"} である。
上端付近で $q=\pi-x$ と置くと $|\dot x|\simeq\sqrt{g/l}|x|$ となり、$\int\dd x/|x|$ は発散する。
従ってこの軌道上で上端に到着するには無限時間を要する。

### 固定点の周囲を線形化する

自律系で $J$ は可逆だから固定点は $\nabla H=0$ の点である。
一自由度の固定点 $z_*$ で

$$
B=H''(z_*)=\begin{pmatrix}a&b\\b&c\end{pmatrix}
$$

とすると、変位の発展行列は $JB=\begin{pmatrix}b&c\\-a&-b\end{pmatrix}$。
跡は零で行列式は $ac-b^2=\det B$ であり、固有値は $\lambda^2=-\det B$ を満たす。
$\det B<0$ なら伸長と収縮の実方向を持つ双曲型、$\det B>0$ なら線形化は楕円型である。
一自由度で後者なら $H$ は厳密な局所極値で、近くの閉じた等高線が運動を囲む。
$\det B=0$ では線形化だけで判定できない。

:::example[点を固定する流れと面積]
反転した振動子 $H=p^2/(2m)-m\nu^2q^2/2$ を考える。
$p=\pm m\nu q$ は $H=0$ の直線で、この上では $\dot q=\pm\nu q$。
$u=p+m\nu q,\ v=p-m\nu q$ と置けば $\dot u=\nu u,\ \dot v=-\nu v$ になる。
$u$ 方向に $e^{\nu t}$ 倍伸び、$v$ 方向に同じ割合で縮む。
近くの初期値が指数的に離れても、相体積の保存とは矛盾しない。
:::

### 二次 Hamiltonian の流れと分布の形

固定点に集めた初期値の誤差は、相面積が一定でも時間とともに形を変える。次元をそろえるため $q=q_0 Q,p=p_0 P,t=t_0\tau$、$H=(p_0q_0/t_0)h$ とする。無次元の列 $x=(Q,P)^\top$ について $h=x^\top Bx/2$、$B$ を実対称行列とすれば $x'=Ax,\ A=JB$。二次行列の恒等式から $A^2=-\det B\,1$ である。従って行列指数の偶数次と奇数次を分ければ、$\det B=\omega^2>0$ なら $S=1\cos\omega\tau+A\sin(\omega\tau)/\omega$、$\det B=-\nu^2<0$ なら三角関数を双曲線関数に替え、$\det B=0$ なら $S=1+\tau A$ となる。いずれも $S(0)=1,S'=AS$ である。

初期値を多数測り、平均を $\bar x$、共分散を $\Sigma=\langle(x-\bar x)(x-\bar x)^\top\rangle$ とする。線形な流れでは各偏差に $S$ を掛けるので $\Sigma(\tau)=S\Sigma(0)S^\top$。$\det S=1$ なら共分散行列の行列式も変わらない。これは線形な流れでの主張であり、非線形な流れで雲が曲がると体積保存だけでは共分散を決められない。

:::example[寄与が混ざる楕円型]
基準時刻を $t_0=1\,\mathrm s$ として無次元表示を使う。$B=\left(\begin{matrix}2 & 1\\1 & 1\end{matrix}\right)$ の Hamiltonian で $\Sigma(0)=\left(\begin{matrix}1 & 0\\0 & 4\end{matrix}\right)$ とする。運動方程式の行列と指定時刻の写像は

$$
A=\left(\begin{matrix}1 & 1\\-2 & -1\end{matrix}\right),\quad A^2=\left(\begin{matrix}-1 & 0\\0 & -1\end{matrix}\right),\quad S(\frac{\pi}{2})=\left(\begin{matrix}1 & 1\\-2 & -1\end{matrix}\right).
$$

$\det S=1$ と $S^\top B S=B$ を行列積で確認できる。前者は面積、後者は各初期値のエネルギーの保存であり、別の等式である。中心が初め $(1,0)$ にあれば $\left(\begin{matrix}1 & -2\end{matrix}\right)^\top$ へ移る。中心の移動を引いた後の共分散は

$$
S\Sigma(0)S^\top=\left(\begin{matrix}5 & -6\\-6 & 8\end{matrix}\right),\qquad \det\Sigma(\tau)=4.
$$

固有振動数が一でも、測定した二軸で見た楕円の主軸は傾く。初期の分布がエネルギーだけの関数でないなら、各粒子のエネルギーが保存しても分布は静止しない。
:::

:::example[鞍点の近くの細長い束]
基準時刻を $t_0=1\,\mathrm s$ として無次元表示を使う。$B=\left(\begin{matrix}0 & 1\\1 & 0\end{matrix}\right)$ の Hamiltonian で $\Sigma(0)=\left(\begin{matrix}2 & 1\\1 & 1\end{matrix}\right)$ とする。運動方程式の行列と指定時刻の写像は

$$
A=\left(\begin{matrix}1 & 0\\0 & -1\end{matrix}\right),\quad A^2=\left(\begin{matrix}1 & 0\\0 & 1\end{matrix}\right),\quad S(\ln{\left(2 \right)})=\left(\begin{matrix}2 & 0\\0 & \frac{1}{2}\end{matrix}\right).
$$

$\det S=1$ と $S^\top B S=B$ を行列積で確認できる。中心が初め $(1,0)$ にあれば $\left(\begin{matrix}2 & 0\end{matrix}\right)^\top$ へ移る。中心の移動を引いた後の共分散は

$$
S\Sigma(0)S^\top=\left(\begin{matrix}8 & 1\\1 & \frac{1}{4}\end{matrix}\right),\qquad \det\Sigma(\tau)=1.
$$

正の時間で一方の幅が広がり他方が狭くなる。位置だけを測ると粒子の密度は薄くなっても、同時に運動量を測れば曲線束の面積は保たれる。
:::

:::example[ばねがない自由な流れ]
基準時刻を $t_0=1\,\mathrm s$ として無次元表示を使う。$B=\left(\begin{matrix}0 & 0\\0 & 2\end{matrix}\right)$ の Hamiltonian で $\Sigma(0)=\left(\begin{matrix}1 & -1\\-1 & 2\end{matrix}\right)$ とする。運動方程式の行列と指定時刻の写像は

$$
A=\left(\begin{matrix}0 & 2\\0 & 0\end{matrix}\right),\quad A^2=\left(\begin{matrix}0 & 0\\0 & 0\end{matrix}\right),\quad S(2)=\left(\begin{matrix}1 & 4\\0 & 1\end{matrix}\right).
$$

$\det S=1$ と $S^\top B S=B$ を行列積で確認できる。中心が初め $(1,0)$ にあれば $\left(\begin{matrix}1 & 0\end{matrix}\right)^\top$ へ移る。中心の移動を引いた後の共分散は

$$
S\Sigma(0)S^\top=\left(\begin{matrix}25 & 7\\7 & 2\end{matrix}\right),\qquad \det\Sigma(\tau)=1.
$$

固定点は孤立せず、行列の二乗が零になる。時間に線形なずれを指数的な不安定性と取り違えてはいけない。
:::

:::example[逆向きに進む閉軌道]
基準時刻を $t_0=1\,\mathrm s$ として無次元表示を使う。$B=\left(\begin{matrix}-1 & 0\\0 & -1\end{matrix}\right)$ の Hamiltonian で $\Sigma(0)=\left(\begin{matrix}3 & 1\\1 & 1\end{matrix}\right)$ とする。運動方程式の行列と指定時刻の写像は

$$
A=\left(\begin{matrix}0 & -1\\1 & 0\end{matrix}\right),\quad A^2=\left(\begin{matrix}-1 & 0\\0 & -1\end{matrix}\right),\quad S(\frac{\pi}{2})=\left(\begin{matrix}0 & -1\\1 & 0\end{matrix}\right).
$$

$\det S=1$ と $S^\top B S=B$ を行列積で確認できる。中心が初め $(1,0)$ にあれば $\left(\begin{matrix}0 & 1\end{matrix}\right)^\top$ へ移る。中心の移動を引いた後の共分散は

$$
S\Sigma(0)S^\top=\left(\begin{matrix}1 & -1\\-1 & 3\end{matrix}\right),\qquad \det\Sigma(\tau)=2.
$$

エネルギーが極大でも孤立した一自由度の系では近くの等高線が閉じる。散逸が加わった系の安定性とは区別する。
:::

## 修正された Hamilton の原理

### 位相空間で経路を独立に動かす

従来の :ref[prin-am-hamilton] では位置経路だけを変分し、速度の変化はその時間微分で決まった。
今回は $q,p$ を独立な曲線として

$$
\mathcal S[q,p]=\int_{t_1}^{t_2}\left(\sum_{i=1}^n p_i\dot q^i-H(q,p,t)\right)\dd t
\label{eq-am-phase-action}
$$

を用いる。
時刻と端点の位置を固定して $\delta q^i(t_1)=\delta q^i(t_2)=0$ とし、運動量端点は自由にする。
積分を一回部分積分すれば

$$
\delta\mathcal S=\left[\sum_{i=1}^n p_i\delta q^i\right]_{t_1}^{t_2}
+\int_{t_1}^{t_2}\sum_{i=1}^n
\left[(\dot q^i-H_{p_i})\delta p_i-(\dot p_i+H_{q^i})\delta q^i\right]\dd t.
$$

端点項が消え、:ref[lem-am-fundamental-variation] を各変分に用いると停留条件は :eqref[eq-am-hamilton-equations] である。
第一式を解いて運動量を消去すれば、正則な場合に元の $\int L\,\dd t$ が戻る。

### 端点項とエネルギー固定の原理

両端で運動量を指定するなら、$\mathcal S$ から $[\sum_i p_i q^i]_{t_1}^{t_2}$ を引く。
新しい端点項は $-[\sum_i q^i\delta p_i]_{t_1}^{t_2}$ となり、内部の運動方程式は変わらない。
例えば自由粒子の軌道で時間差を $\Delta t$、変位を $\Delta q$ とすると $S=m(\Delta q)^2/(2\Delta t)$。
終点位置の微分は $p=m\Delta q/\Delta t$、終点時刻の微分は $-p^2/(2m)$ である。

自律系でエネルギー $E$ を指定して経路の形だけを求めるとき、短縮作用 $W=\int\sum_i p_i\,\dd q^i$ を使える。
任意のパラメータ $s$ と乗数 $\lambda(s)$ を用い $\int[\sum_i p_i(q^i)'-\lambda(H-E)]\dd s$ を変分すると $H=E,\ (q^i)'=\lambda H_{p_i},\ p_i'=-\lambda H_{q^i}$。
$\lambda>0$ の枝をとり、$\dd t=\lambda\,\dd s$ と置けば、元の時間向きの正準方程式になる。端点の位置は固定し、運動量はエネルギー条件の範囲で変分する。
この原理では移動時間があらかじめ固定されていない。
$H=|\vb p|^2/(2m)+V$ なら運動方向に $\vb p$ が平行なので $W=\int\sqrt{2m(E-V)}\,\dd l$。
一定の $V$ なら停留曲線は直線となる。

### 両端を固定した解を代入してよい順序

変分をとる前に $p=m\dot q$ を使えば Lagrange 形式に戻れる。
ただしこの消去は $p$ の変分による代数方程式を解いたものである。
$q$ の軌道まで先に具体的な解に固定すると、変分すべき方向が消えてしまう。
例として $H=p^2/(2m)+m\omega^2q^2/2$ で $p$ を平方完成すると、

$$
p\dot q-H=-\frac{(p-m\dot q)^2}{2m}
+\frac m2\dot q^2-\frac{m\omega^2}{2}q^2 .
$$

各時刻での $p$ の停留値は明らかに $m\dot q$ である。
その後に残る $q$ の変分から $\ddot q+\omega^2q=0$ を得る。
$p$ の方向には二次変分が負なので、停留作用を全変数に関する最小値と呼ぶこともできない。

:::example[固定時刻と固定エネルギーの違い]
距離 $d$ の直線を時間 $\Delta t$ で進む自由粒子なら $S_{\rm cl}=md^2/(2\Delta t)$。
これに $E\Delta t$ を加えて $\Delta t$ を変分すると

$$
-\frac{md^2}{2\Delta t^2}+E=0,\qquad
\Delta t=d\sqrt{m/(2E)} .
$$

停留値は $d\sqrt{2mE}=W$。
時刻とエネルギーも、端点データとしては Legendre 変換で交換できる。
$m=2\,\mathrm{kg},d=3\,\mathrm m,E=9\,\mathrm J$ では $\Delta t=1\,\mathrm s$、
$S_{\rm cl}=9\,\mathrm{J\,s}$、$W=18\,\mathrm{J\,s}$ である。
:::

### 時刻も座標に取る拡張相空間

時間依存系を自律系として扱うため、$q^0=t$ と共役な $p_0$ を追加し、任意のパラメータ $s$ に対して

$$
\mathcal K(q^0,q,p_0,p)=p_0+H(q,p,q^0)
$$

とする。
正準方程式は

$$
\begin{aligned}
\frac{\dd q^0}{\dd s}&=1,&
\frac{\dd p_0}{\dd s}&=-\partial_tH,\\
\frac{\dd q^i}{\dd s}&=H_{p_i},&
\frac{\dd p_i}{\dd s}&=-H_{q^i}.
\end{aligned}
$$

従って $s$ と $t$ は定数差を除いて等しい。
$\mathcal K$ は自律的なのでその値を零と選べば $p_0=-H$ が保たれる。
拡張空間の正準一形式 $\sum_{i=0}^n p_i\dd q^i$ をこの面に制限したものが $\sum_{i=1}^n p_i\dd q^i-H\dd t$ である。
時間とエネルギーを一組として扱う理由が、端点の作用からも正準方程式からも得られる。

例えば $H=p^2/(2m)+f(t)q$ では $\dot p_0=-f'(t)q$。
実際の全エネルギーが $f'(t)q$ だけ変わる分を $p_0$ が相殺する。
拡張変数を加えたことだけで実験室のエネルギーが保存するわけではない。

### 一つの軌道から軌道の束へ

初期の閉曲線をパラメータ $u\in[0,1]$ で表し、各点を同じ Hamiltonian で時間発展させる。
各時刻で初期点 $u=0$ と $u=1$ は同じ点に移る。
曲線に沿う積分

$$
C(t)=\int_0^1\sum_i p_i(u,t)\,\partial_u q^i(u,t)\,\dd u
$$

を考え、その時間微分を取る。
積の第二項を $u$ で部分積分すると端点の差は相殺され、

$$
\begin{aligned}
\dot C(t)
&=\int_0^1\sum_i(\dot p_i\,\partial_u q^i-\partial_u p_i\,\dot q^i)\dd u\\
&=-\int_0^1\sum_i(H_{q^i}\partial_u q^i+H_{p_i}\partial_u p_i)\dd u\\
&=-\int_0^1\partial_u H\,\dd u=0 .
\end{aligned}
$$

ここでは閉曲線上を動く間は $t$ が同じなので $H_t$ の項は生じない。
この :term[Poincaré の積分不変量]{en="Poincaré integral invariant" yomi="ぽあんかれのせきぶんふへんりょう"} は周期軌道一つに限らず、
任意の閉曲線を相流に乗せた場合に成り立つ。
一自由度で閉曲線が単純なら $\oint p\,\dd q$ の絶対値が囲む面積に等しく、位相面積保存の別証明になる。

### 駆動系で初期値を戻す

時間に依存する Hamiltonian でも相流は可逆である。
ただし、ある時刻で速度の符号を変えてそのまま同じ駆動を続けることと、過去へ積分して初期値を戻すことは一般に一致しない。
例として

$$
H=\frac{p^2}{2m}-f(t)q,\qquad f(t)=f_0+\gamma t
$$

を考える。
時刻 $0$ に $(q_0,p_0)$ を与えて積分すると

$$
\begin{aligned}
p(t)&=p_0+f_0t+\gamma t^2/2,\\
q(t)&=q_0+p_0t/m+f_0t^2/(2m)+\gamma t^3/(6m).
\end{aligned}
$$

この流れは自由なせん断に時間依存の並進を合成したものだから、初期値に関する行列式は $1$ である。
終点 $(q_f,p_f)$ と駆動の履歴を知っていれば、第一式で $p_0$ を戻し、第二式で $q_0$ を戻せる。

終点で運動量の符号だけを反転してその後も $f(t)$ を増やし続ける実験では、振り返るべき古い力の履歴を再現していない。
逆向きに同じ位置列をたどらせるには、時間反転で奇な変数を反転させ、駆動スケジュールも逆順にする必要がある。
相流に逆写像が存在する条件と、装置が時間反転に対称である条件を区別できる例である。

### 経路が交わる時刻と端点作用の枝

二次の流れを物理的な正準変数で $(q_2,p_2)^\top=\left(\begin{smallmatrix}a&b\\c&d\end{smallmatrix}\right)(q_1,p_1)^\top$ と表す。$ad-bc=1$ とする。時刻を固定した実軌道の作用を端点で微分すれば $\dd S=p_2\dd q_2-p_1\dd q_1$ である。$b\ne0$ なら $p_1=(q_2-aq_1)/b,\ p_2=(d q_2-q_1)/b$。交差微分がともに $-1/b$ となるので積分でき、
$$
S(q_2,q_1)=\frac{d q_2^2-2q_2q_1+a q_1^2}{2b}+\text{定数}
$$
を得る。定数は端点には依存せず、Hamiltonian に加えた時刻だけの関数や作用の基準で決まる。$b=0$ のときは初期位置を固定すると、すべての初期運動量が同じ終点位置に集まる。位置間の作用を一価に作れないが、相空間全体での写像はなお可逆である。

二つの区間の作用 $S_{21},S_{32}$ をつなぐには $S_{32}(q_3,q_2)+S_{21}(q_2,q_1)$ を $q_2$ について停留させる。微分が零という条件は中間時刻の左右で運動量が等しいことを表す。自由粒子なら二項は $m(q_3-q_2)^2/(2T_2)+m(q_2-q_1)^2/(2T_1)$。停留点は $q_2=(T_1q_3+T_2q_1)/(T_1+T_2)$、代入した作用は $m(q_3-q_1)^2/[2(T_1+T_2)]$。時間の和に対応する流れの合成則が、作用では中間点の停留条件として現れる。

### 固定エネルギーの変分を Jacobi 計量で書く

拘束面に沿う自然系で $T=\sum_{i,j}M_{ij}(q)\dot q^i\dot q^j/2$、$M$ は正定値とする。$E>V(q)$ の領域では速度の大きさが位置で決まる。短縮作用の要素は
$$
\sum_i p_i\dd q^i
 =\sqrt{2(E-V)\sum_{i,j}M_{ij}\dd q^i\dd q^j}.
$$
実際 $\dd q=\dot q\,\dd t$ を入れると両辺は $2(E-V)\dd t$ になる。従って $g^J_{ij}=2(E-V)M_{ij}$ を :term[Jacobi 計量]{en="Jacobi metric" yomi="やこびけいりょう"} とすれば、物理軌道はこの計量で長さを停留させる曲線である。距離パラメータ $s_J$ は作用の単位をもち、物理時刻の復元には $\dd t=\dd s_J/[2(E-V)]$ を使う。$E=V$ の転回点では計量が退化する。そこを通常の正則な距離の点として通過させず、物理時刻で二つの枝を接続する必要がある。

例えば $M_{ij}=m\delta_{ij}$ の平面で $V=mgy$ とし、$y$ 方向の傾き $y'=\dd y/\dd x$ が有限な区間を取る。長さの被積分関数は $\sqrt{2m(E-mgy)}\sqrt{1+y'^2}$。$x$ が陽に現れないので :ref[prop-am-beltrami] より $\sqrt{2m(E-mgy)}/\sqrt{1+y'^2}=p_x$ が一定となる。$x$ が増える枝を取った。平方して $y'^2=2m(E-mgy)/p_x^2-1$、傾きが非零の区間で微分すると $y''=-m^2g/p_x^2$。滑らかさで頂点にも続けられる。初期の傾き $s_0$ と高さ $y_0$ を用いれば $y=y_0+s_0x-m^2gx^2/(2p_x^2)$。等時刻の運動方程式を一度も解かずに放物線の形が出た。$p_x=0$ の鉛直軌道は $x$ を独立変数にした段階で除いたので、別に扱う。

## Liouville の定理

数多くの初期値を小さな領域に置くと、異なる点が異なる速度で動くため領域の形は変わる。
その体積がどう変わるかを、微小な直方体の辺の時間変化から求める。

:::theorem[Liouville の体積保存定理]{#thm-am-liouville-volume}
$H$ が相空間変数について $C^2$ 級で、相空間変数の二階までの偏微分が時刻についても連続であるとする。相流が存在する区間では、正準座標での体積 $\dd\Gamma=\dd q^1\cdots\dd q^n\,\dd p_1\cdots\dd p_n$ は保存する。
$H$ は時間に陽に依存してもよい。
:::

:::proof
$X=J\nabla H$ の発散は混合微分の交換により

$$
\sum_{i=1}^n(\partial_{q^i}H_{p_i}-\partial_{p_i}H_{q^i})=0 .
$$

$D=\partial\Phi_{t,s}/\partial z$ を初期値で微分した行列とすれば $\dot D=(\partial_z X)D,\ D(s)=I$。
従って微小時間 $h$ に $\det D(t+h)=\det(I+h\partial_z X+o(h))\det D(t)$。
行列式の一次項は対角成分の和だから $\dd(\det D)/\dd t=(\divg X)\det D=0$。
初期値が $1$ なので常に $\det D=1$。
積分の変数変換により任意の運ばれる領域の体積が等しい。
:::

例えば自由粒子の束を $0<q<d,\ p_-<p<p_+$ に一様に置く。時刻 $t$ の四角形の左右の辺は $q=tp/m$ と $q=d+tp/m$ であり、上下の辺は初めの高さにある。各 $p$ で水平に切った幅は常に $d$ だから、面積は $d(p_+-p_-)$ のままである。一方、位置の最小値から最大値までの幅は $t>0$ で $d+t(p_+-p_-)/m$ に広がる。相空間の面積と、位置軸に射影して測る粒子束の幅との違いを図から直接読み取れる。

```tikz caption="自由粒子の位相流。初期に同じ位置にある二点も、運動量が異なれば異なる距離を進む。"
\draw[->] (-.4,0)--(4.8,0) node[right] {$q$};
\draw[->] (0,-.3)--(0,2.7) node[above] {$p$};
\draw[dashed] (.55,.35) rectangle (1.8,1.7);
\draw[thick] (.9,.35)--(2.15,.35)--(3.5,1.7)--(2.25,1.7)--cycle;
\draw[->] (1.8,1.73)--(3.45,1.73);
\draw[->] (1.8,.25)--(2.14,.25);
```

### 密度の輸送

体積要素に含まれる系の個数を $\rho(z,t)\dd\Gamma$ とする。
流れとともに動く領域から出入りがなければ、固定領域で $\partial_t\rho+\divg(\rho X)=0$。
発散零を使うと

$$
\partial_t\rho+\sum_{i=1}^n\left(H_{p_i}\partial_{q^i}\rho-H_{q^i}\partial_{p_i}\rho\right)=0 .
$$

これは軌道に沿って $\dd\rho/\dd t=0$ という意味である。
空間の同じ点で密度が不変とは限らない。
自由粒子なら $\rho(q,p,t)=\rho_0(q-pt/m,p)$ であり、長方形は平行四辺形に変わる。
自律系では任意の $\rho=f(H)$ が定常解になる。

:::example[抵抗のある振動子]
$\dot q=p/m,\ \dot p=-kq-\gamma p$ の発散は $-\gamma$。
面積は $e^{-\gamma t}$ 倍に縮み、$\gamma=1\,\mathrm{s^{-1}}$ なら $t=\ln2\,\mathrm s$ で半分になる。
この $(q,p=m\dot q)$ の流れは閉じた正準 Hamilton 系の流れと異なる。
:::

```tikz caption="自由粒子の相流。破線の初期領域は運動量ごとに異なる距離だけ右へ移り、同じ面積の平行四辺形になる。"
\draw[->] (-.3,0) -- (5.1,0) node[right] {$q$};
\draw[->] (0,-.2) -- (0,2.5) node[above] {$p$};
\draw[dashed,thick] (.35,.55) rectangle (1.75,1.85);
\draw[thick] (.9,.55) -- (2.3,.55) -- (3.6,1.85) -- (2.2,1.85) -- cycle;
\draw[->] (.35,2.02) -- (2.2,2.02);
\draw[->] (.35,.4) -- (.9,.4);
```

### 正準座標以外で体積を数える

同じ物理的な領域でも座標の選び方によって体積要素に係数が現れる。
例えば平面極座標で速度 $(\dot r,\dot\phi)$ を使うと $p_r=m\dot r,\ p_\phi=mr^2\dot\phi$。
$r,\phi$ を固定したときの運動量変換の行列式が $m^2r^2$ だから、

$$
\dd\Gamma=m^2r^2\,\dd r\,\dd\phi\,\dd\dot r\,\dd\dot\phi .
$$

位置によって係数が変わる。
速度座標の小箱を単に数えたときに体積が変わっても Liouville の定理には反しない。
同じ極座標で、直交速度成分の動径・接線成分 $v_r,v_\phi$ を使えば $v_\phi=r\dot\phi$ で、係数は $m^2r$ になる。

### 分布の定常性と混合

自律系で $e^{-\beta H}$ の相空間積分が有限となる $\beta$ を選び、初期密度を $\rho_0(q,p)=C e^{-\beta H(q,p)}$ とすると、$\beta$ の単位は $\mathrm{J^{-1}}$ である。
$H$ の値が軌道上で変わらないため $\rho_0$ 自体が定常である。
一方、同じ Hamiltonian に非対称な分布を載せれば、密度の値を保ったまま模様が移動する。
例えば振動子で細い初期集団を $q>0$ に集めると、半周期後には $q<0$ に移る。
相流は初期の模様を引き戻せるので、体積保存だけから一様分布への緩和時間は決まらない。

:::example[流体粒子のように運ばれる確率]
自由粒子で $0<q_0<L,\ 0<p_0<P$ に一様な確率密度を置く。
規格化定数は $1/(LP)$。
時刻 $t$ には $0<q-pt/m<L,\ 0<p<P$ が支持で、密度値は依然 $1/(LP)$ である。
位置だけの分布は、各 $q$ に対してこの条件を満たす $p$ の区間長を $LP$ で割ったものになり、一般には初期値と異なる。
全相空間の保存則と部分変数の分布は別々に扱う必要がある。
:::

### エネルギー殻の体積要素

Liouville の定理で保存される $\dd\Gamma$ は $2n$ 次元の体積である。
自律系でエネルギーを固定すると軌道は $2n-1$ 次元の面にあるため、面上の分布を考えるには厚さの取り方を決める必要がある。
相変数の各軸に一定の基準量を選んで無次元座標 $z$ にそろえる。その座標で $\dd\Gamma=C\,\dd z$ とし、$C$ は各軸の基準量の積である。正則なエネルギー $E$ で幅 $\Delta E$ の層を切り出す。以下の展開は連続な $g$ がコンパクトな台をもつ場合、または同じ極限を積分の下で取れる場合に使う。
法線方向の距離は一次まで $\Delta E/|\nabla H|$ なので

$$
\int_{E<H<E+\Delta E} g(z)\dd\Gamma
=C\Delta E\int_{H=E}\frac{g(z)}{|\nabla_z H|}\dd\Sigma+o(\Delta E)
$$

となる。
これは局所的に $H$ を法線方向の座標の代わりに使う変数変換である。
相流は $\dd\Gamma$ とエネルギーの各層を保つ。層の積分を幅で割って極限を取れば、殻上の測度 $C\,\dd\Sigma/|\nabla_z H|$ も保存する。測度の重み $1/|\nabla_z H|$ は同じエネルギー幅が座標ごとに異なる厚さに対応することから生じる。
臨界点では分母が零になるため正則面の議論をそのまま使えない。

:::example[振動子の状態数に対応する面積]
振動子で $\Gamma(E)=\int_{H<E}\dd q\dd p$ を求める。
楕円の半軸の長さは

$$
q_{\max}=\sqrt{2E/(m\omega^2)},\qquad p_{\max}=\sqrt{2mE}.
$$

従って $\Gamma(E)=2\pi E/\omega$、薄い層の単位エネルギーあたりの面積は $\Gamma'(E)=2\pi/\omega$。
これは振動の周期にも等しい。
一般の一次元束縛運動についても上半分と下半分の運動量の差を積分すれば

$$
\Gamma(E)=2\int_{q_-}^{q_+}\sqrt{2m(E-V(q))}\,\dd q .
$$

転回点では被積分関数が零なので端点微分は寄与せず、$E$ で微分すると往復時間が得られる。
面積の変化率が実際に観測する時間を決める一例である。
:::

### 位相図の密度と位置の出現確率

一つの振動子を同じエネルギーで長時間観測する。
軌道が一つの周期曲線を往復するだけでも、位置を等間隔に通るわけではない。
幅 $\dd q$ の位置区間を通るのに、片道で $\dd t=\dd q/|\dot q|$ かかるから、一周期における位置の確率密度は

$$
P(q)=\frac{2}{T|\dot q|}=\frac{1}{\pi\sqrt{A^2-q^2}},\qquad |q|<A .
$$

$q=A\sin\alpha$ と置けば積分は $1$ となる。
転回点で発散するのはそこで遅く動くためで、有限な幅にいる確率は有限である。
例えば $|q|<A/2$ にいる時間の割合は $2\arcsin(1/2)/\pi=1/3$ であり、残りの $2/3$ を外側で過ごす。
相空間での流れの非圧縮性と、位置だけで見た滞在時間の偏りは両立する。

### 有限体積での再帰とその限界

自律系の不変領域 $D$ が有限の正準体積を持ち、相流が $D$ 上で全時間定義されるとする。
一定時間 $\tau>0$ ごとの写像を $T$ とする。
$A\subset D$ を可測集合として、その中で将来一度も $A$ に戻らない点の集合を $B$ とする。
$B,TB,T^2B,\ldots$ は対ごとに交わらない。
実際、$T^jb_1=T^kb_2$、$k>j$ なら $T^{k-j}b_2=b_1\in A$ となり $b_2$ の定義に反する。
すべてが $D$ に入り、体積は等しいから $B$ の体積は零である。

有限回だけ戻る点も、最後の帰還点は $B$ に入るので $\bigcup_{j\ge0}T^{-j}B$ に含まれ、その体積は零である。
従って 体積に関して $A$ のほとんどすべての点が $A$ に何度でも戻る。
これが Poincaré の再帰定理の集合に関する形である。
戻るまでの時間の有用な上限や、点が領域全体を一様に訪れることはこの証明からは出ない。
無限空間の自由粒子は有限の不変領域という仮定から外れる。

### 中心力で位相体積の要素を具体化する

平面の直交座標 $(x,y)$ から極座標 $(r,\phi)$ に替えると位置面積は $r\,\dd r\,\dd\phi$ となる。
一方、角度運動量 $p_\phi$ は接線方向の通常の運動量 $p_{\rm tan}$ の $r$ 倍である。
そこで位置を固定して比較すると

$$
\dd p_x\dd p_y=\dd p_r\dd p_{\rm tan}
=\frac1r\,\dd p_r\dd p_\phi .
$$

二つを掛けると位置の $r$ と運動量の $1/r$ が相殺される。
正準座標での体積が $\dd r\,\dd\phi\,\dd p_r\,\dd p_\phi$ という同じ形に保たれることを、正準変換の一般論を使わずに確かめられた。

例えば半径を $r$ から $2r$ へ移しても同じ接線速度の幅だけを選ぶと、対応する $p_\phi$ の幅は二倍になる。
各半径で**同じ角運動量の幅**を選ぶ集団と、**同じ接線速度の幅**を選ぶ集団は、物理的に異なる初期条件を数えている。
位相体積を比較するには、どの変数をどの範囲で固定したかを明示する。

:::example[回転する基準で測った時刻と正準運動量]
半径 $a$ の滑らかな水平円周を質量 $m$ の粒子が動く。
角度の目盛りだけを角速度 $\Omega$ で回し、新しい角を $q=\phi-\Omega t$ とする。
$C=ma^2$ と置けば

$$
L=\tfrac C2(\dot q+\Omega)^2,\qquad p=C(\dot q+\Omega).
$$

従って $H=p^2/(2C)-\Omega p$。
$H_p=p/C-\Omega$ は目盛りから見た角速度であり、$p$ は実験室での角運動量である。
$C=2\,\mathrm{kg\,m^2}$、$\Omega=1\,\mathrm{s^{-1}}$、$p=4\,\mathrm{J\,s}$ なら実験室の角速度は $2\,\mathrm{s^{-1}}$、
目盛りからは $1\,\mathrm{s^{-1}}$ となる。
正準エネルギーはこの数値例で $0\,\mathrm J$ だが、実験室の運動エネルギーは $4\,\mathrm J$ である。
両者の差 $\Omega p$ は、時間を進めながら角の原点も回したことに由来する。
:::

### エネルギー殻の密度と平均の重み

定常な位相密度として $H$ の関数を選べる。そのうち薄いエネルギー層で一様な測度を取り、層の厚さで割った極限を考える。エネルギー $E$ 以下の正準体積を $\Gamma(E)=\int_{H<E}\dd q\,\dd p$、その微分を $g(E)$ とする。エネルギー区間 $[E,E+\dd E]$ の体積は一次まで $g(E)\dd E$ なので、殻上の平均は
$$
\langle A\rangle_E=\frac{1}{g(E)}\frac{\dd}{\dd E}\int_{H<E} A(q,p)\,\dd q\,\dd p
$$
である。$A$ は $E$ に陽に依存しないとした。相空間に追加の積分があれば一本の軌道は殻全体を訪れない場合がある。測度を定義しただけで、この平均を任意の初期点からの時間平均と同一視する根拠は得られない。

独立な調和振動子 $n$ 個では各 $q_j,p_j$ の平面を $x_j=\sqrt{m_j}\omega_j q_j,\ y_j=p_j/\sqrt{m_j}$ に変えると、$H=\sum_j(x_j^2+y_j^2)/2$、面積要素は $\dd q_j\dd p_j=\dd x_j\dd y_j/\omega_j$。各平面で角を積分すれば $2\pi\,\dd E_j/\omega_j$ が残る。非負の $E_j$ の和が $E$ 以下の単体の体積は $E^n/n!$。これは $n-1$ 個の場合を $E_n=0$ から $E$ まで積分して帰納的に示せる。従って
$$
\Gamma(E)=\frac{(2\pi)^n E^n}{n!\prod_{j=1}^n\omega_j},\qquad
g(E)=\frac{(2\pi)^n E^{n-1}}{(n-1)!\prod_{j=1}^n\omega_j}.
\label{eq-am-oscillator-shell}
$$

:::example[1 自由度の振動子のエネルギーを数える]
角振動数を $(2)\,\mathrm{s^{-1}}$、全エネルギーを $4\,\mathrm J$ とする。質量の違いは $q,p$ 軸を逆方向に伸縮するので体積には現れない。角積分のあとに残る単体を積分すると、SI の数値で

$$
\Gamma(E)=\pi E,\qquad g(E)=\pi,\qquad g(4)=\pi.
$$

$\Gamma$ の単位は $(\mathrm{J\,s})^{1}$、$g$ はそれを J で割ったものである。殻ではエネルギー変数に関して対称なので平均は $\langle E_1\rangle=E/n=4\,\mathrm J$。殻の密度は周期に等しく、エネルギーに依存しない。
:::

:::example[2 自由度の振動子のエネルギーを数える]
角振動数を $(2,3)\,\mathrm{s^{-1}}$、全エネルギーを $6\,\mathrm J$ とする。角積分のあとに残る単体を積分すると、SI の数値で

$$
\Gamma(E)=\frac{\pi^{2} E^{2}}{3},\qquad g(E)=\frac{2 \pi^{2} E}{3},\qquad g(6)=4 \pi^{2}.
$$

$\Gamma$ の単位は $(\mathrm{J\,s})^{2}$、$g$ はそれを J で割ったものである。殻ではエネルギー変数に関して対称なので平均は $\langle E_1\rangle=E/n=3\,\mathrm J$。第一モードのエネルギーは殻上で一様に分布する。初期点を一つ選んだあとでは各モードのエネルギーが常に一定であることと両立する。
:::

:::example[3 自由度の振動子のエネルギーを数える]
角振動数を $(1,2,4)\,\mathrm{s^{-1}}$、全エネルギーを $8\,\mathrm J$ とする。角積分のあとに残る単体を積分すると、SI の数値で

$$
\Gamma(E)=\frac{\pi^{3} E^{3}}{6},\qquad g(E)=\frac{\pi^{3} E^{2}}{2},\qquad g(8)=32 \pi^{3}.
$$

$\Gamma$ の単位は $(\mathrm{J\,s})^{3}$、$g$ はそれを J で割ったものである。殻ではエネルギー変数に関して対称なので平均は $\langle E_1\rangle=E/n=\frac{8}{3}\,\mathrm J$。第一モードに大きなエネルギーを渡すほど、残りのモードがとれる状態の体積が小さくなる。
:::

:::example[4 自由度の振動子のエネルギーを数える]
角振動数を $(1,2,3,5)\,\mathrm{s^{-1}}$、全エネルギーを $10\,\mathrm J$ とする。角積分のあとに残る単体を積分すると、SI の数値で

$$
\Gamma(E)=\frac{\pi^{4} E^{4}}{45},\qquad g(E)=\frac{4 \pi^{4} E^{3}}{45},\qquad g(10)=\frac{800 \pi^{4}}{9}.
$$

$\Gamma$ の単位は $(\mathrm{J\,s})^{4}$、$g$ はそれを J で割ったものである。殻ではエネルギー変数に関して対称なので平均は $\langle E_1\rangle=E/n=\frac{5}{2}\,\mathrm J$。殻の幅を変えると状態数はほぼ比例して増えるが、平均した各モードの取り分は正規化によって幅に依存しなくなる。
:::

### 殻の一部だけを観測した分布

殻の平均値だけでなく、観測した一つのモードのエネルギー分布も求められる。$n\ge2$ とし、第一モードのエネルギーを $s$ から $s+\dd s$ に限る。第一平面が占める面積は $(2\pi/\omega_1)\dd s$。残りの $n-1$ 個は全エネルギー $E-s$ の殻にいる。従って :eqref[eq-am-oscillator-shell] を使って規格化すると
$$
P_E(s)=\frac{(2\pi/\omega_1) g_{n-1}(E-s)}{g_n(E)}
=\frac{n-1}{E}\left(1-\frac{s}{E}\right)^{n-2},\qquad 0<s<E .
$$
外では密度を零とする。$u=s/E$ と置いて積分すれば一に戻る。$\langle s\rangle=E/n$ だけでなく、$\langle s^2\rangle=2E^2/[n(n+1)]$ も二回の部分積分で得られる。例えば三自由度で $E=8\,\mathrm J$ なら $\Pr(s>4\,\mathrm J)=\int_{1/2}^1 2(1-u)\dd u=1/4$。一つの振動子を見続けて、この確率でエネルギーが勝手に移るという意味ではない。調和モードが独立なら個々の $s$ は保存し、分布は初期状態の集合に関するものである。

$n,E$ をともに大きくして $E/n=\varepsilon_0$ を一定に取ると、固定した $s$ に対して $P_E(s)\to e^{-s/\varepsilon_0}/\varepsilon_0$。対数を展開すると $(n-2)\log[1-s/(n\varepsilon_0)]\to-s/\varepsilon_0$ である。個々のモードの統計が全体系の殻の体積と結び付くことを、相流の混合を仮定せず計算できる。

### 断面を単位時間に通る相体積

密度を薄い層に制限した場合、断面上で点を一様に数える操作と、断面付近に滞在する点を一様に数える操作を区別する必要がある。$q^1=0$ を正向きに横切る系で、時間 $\dd t$ に通る初期点は断面手前の厚さ $\dot q^1\dd t$ にある。エネルギー殻の測度をデルタ関数で略記すると、正向きの流束は
$$
\int_{\dot q^1>0}\delta(E-H)\,\dot q^1\,
\dd p_1\prod_{j=2}^{n}\dd q^j\dd p_j .
$$
積分には $q^1=0$ を代入している。単純根で $p_1$ を $H$ に替えると $\dd H=H_{p_1}\dd p_1=\dot q^1\dd p_1$ であり、二つの因子は相殺する。従って正の速度の根が一つだけなら、残った自由度の正準体積を許容域で積分すればよい。接するだけの点では変数変換が破れるので、横断を仮定する理由も見える。

:::example[鞍点を越える調和モードの流束]
鞍点の近くを二自由度の $H=E_b+p_x^2/(2m)-k x^2/2+p_y^2/(2M)+M\omega^2 y^2/2$ で近似する。$x=0,p_x>0$ の断面では、$y$ 振動のエネルギー $E_y<E-E_b$ のときに限り $p_x$ が決まる。よって上の流束は楕円の面積 $2\pi(E-E_b)/\omega$。$E_b=4\,\mathrm J,E=5\,\mathrm J,\omega=2\,\mathrm{s^{-1}}$ なら $\pi\,\mathrm{J\,s}$ であり、$E$ がしきい値を下回れば零である。これはエネルギー幅あたりの相体積の流れである。実際の粒子数毎秒に換算するには入射する分布の正規化とエネルギー幅が別途要る。近傍を出た後に粒子が戻る現象まで、この局所的な二次模型だけで排除することもできない。
:::

## 演習問題

::::exercise{level=basic}
$L=m\dot x^2/2-U$ で $m=2\,\mathrm{kg}, p=6\,\mathrm{kg\,m/s},U=4\,\mathrm J$ のとき $H$ を求めよ。

:::solution
$\dot x=p/m=3\,\mathrm{m/s}$ で $H=p^2/(2m)+U=9+4=13\,\mathrm J$ となる。
:::
::::

::::exercise{level=basic}
$L=I(\dot q-\Omega)^2/2$（$I,\Omega$ 一定）の $H$ を求め、角速度を復元せよ。

:::solution
$p=I(\dot q-\Omega)$ より $H=p^2/(2I)+\Omega p$。
$H_p=p/I+\Omega=\dot q$ である。
:::
::::

::::exercise{level=basic}
$H=p^2/(2m)+f(t)q$ について $\dot H$ と力を求めよ。

:::solution
正準方程式は $\dot p=-f(t),\dot q=p/m$。
$\dot H=p\dot p/m+f\dot q+f'q=f'q$ であり、力は $-f(t)$ となる。
:::
::::

::::exercise{level=basic}
自由粒子の流れが面積を保つことを直接示せ。

:::solution
$q=q_0+tp_0/m,\ p=p_0$ を初期値で微分すると行列は $\begin{pmatrix}1&t/m\\0&1\end{pmatrix}$。
行列式は $1$ である。
:::
::::

::::exercise{level=advanced}
$L=mv^2/2+\alpha v^4/4$、$m,\alpha>0$ で速度と運動量の対応は大域的に一対一であることを示し、$H$ を $v$ で媒介表示せよ。

:::solution
$p=mv+\alpha v^3$ は導関数 $m+3\alpha v^2>0$ で、両端で $\pm\infty$ に発散する。
従って実数全体への全単射であり、$H=pv-L=mv^2/2+3\alpha v^4/4$ と表せる。
:::
::::

::::exercise{level=advanced}
振動子の相流で初期長方形の各辺が直交し続ける必要がない理由を行列から示せ。

:::solution
$q,p$ を一定の基準単位で無次元表示し、変換行列の二列の内積をとると $\sin\omega t\cos\omega t[(m\omega)^{-1}-m\omega]$（$m\omega$ も同じ基準で表示）となる。
一般に零でない。
一方、行列式は $\cos^2\omega t+\sin^2\omega t=1$ であり面積は一定である。
:::
::::

::::exercise{level=advanced}
質量 $m$ の粒子に $L'=L+DF(q,t)$ を用いたとき、$H'$ を元の $H$ で表せ。

:::solution
$p'=p+F_q$ を Legendre 変換に入れると $\dot qF_q$ が消え、$H'(q,p',t)=H(q,p'-F_q,t)-F_t$ となる。
運動量の引数も同時にずらす必要がある。
:::
::::

::::exercise{level=advanced}
$\dot q=aq,\ \dot p=-ap$（定数 $a$）の Hamiltonian と初期密度の時間発展を求めよ。

:::solution
$H_p=aq$ から $H=aqp+f(q,t)$、第二式から $f_q=0$ なので $f=0$ を選べる。
流れは $(e^{at}q_0,e^{-at}p_0)$。
行列式が $1$ なので $\rho(q,p,t)=\rho_0(e^{-at}q,e^{at}p)$ である。
:::
::::

::::exercise{level=advanced}
球面振り子で $H=p_\theta^2/(2ml^2)+p_\phi^2/(2ml^2\sin^2\theta)-mgl\cos\theta$ とする。
$\theta$ の正準方程式を立て、速度で書いた二階式と比べよ。

:::solution
正準方程式から $\dot\theta=p_\theta/(ml^2),\ \dot\phi=p_\phi/(ml^2\sin^2\theta)$ と

$$
\dot p_\theta=\frac{p_\phi^2\cos\theta}{ml^2\sin^3\theta}-mgl\sin\theta
$$

を得る。
二本の速度の式を代入すると $\ddot\theta=\dot\phi^2\sin\theta\cos\theta-(g/l)\sin\theta$。
$\theta$ と $\phi$ を独立に変分した Lagrange 方程式に等しい。
極で分母が零になるのは方位角を選んだためであり、そこでの運動には別の座標を使う。
:::
::::

::::exercise{level=advanced}
位相平面上で、原点を中心に $q,p$ をともに $e^{-\gamma t}$ 倍にする流れを考える。
$\gamma>0$ とする。
この変数を正準対とした Hamiltonian が存在しないことを二通りで示せ。

:::solution
流れの発散は $-2\gamma$、時間発展の行列式は $e^{-2\gamma t}$ であり、Liouville の定理に矛盾する。
別に、$\dot q=-\gamma q=H_p$ から $H=-\gamma qp+f(q,t)$ と積分すると、$-H_q=\gamma p-f_q$ が $\dot p=-\gamma p$ に等しくなるには $f_q=2\gamma p$ が全 $p$ で必要である。
$f$ は $p$ に依存しないので不可能である。
:::
::::

::::exercise{level=advanced}
ばね定数の代わりに振動子の質量を外部から $m(t)>0$ と指定し、$H=p^2/[2m(t)]+kq^2/2$ とする。全エネルギーの変化と相体積の変化を求めよ。$m=2\,\mathrm{kg}$、$\dot m=0.10\,\mathrm{kg/s}$、$p=4\,\mathrm{kg\,m/s}$ の瞬間を評価せよ。

:::solution
正準方程式は $\dot q=p/m,\dot p=-kq$。全時間微分では $kq\dot q+(p/m)\dot p$ が消え、$\dot H=-p^2\dot m/(2m^2)$ となる。指定の時点では $-0.20\,\mathrm W$。一方、固定時刻で $\partial_q\dot q+\partial_p\dot p=0$ なので、初期領域からの相体積倍率は一である。この模型は与えた時間依存 Lagrangian に対応する。物質が実際に出入りする台車に適用するときは、流入・流出する物質の運動量も別に指定する必要がある。
:::
::::

## 確認テスト

::::quiz{answer="4" tol="1e-6"}
$m=2\,\mathrm{kg},p=8\,\mathrm{kg\,m/s}$ の自由粒子の $\dot q$ を m/s で答えよ。

:::explanation
$H_p=p/m=4\,\mathrm{m/s}$ である。
:::
::::

::::quiz
通常の位相空間作用で両端に固定する量はどれか。

- [x] 時刻と位置
- [ ] 時刻と運動量だけ
- [ ] 速度だけ

:::explanation
部分積分で生じる端点項は $\sum_i p_i\delta q^i$ だから、時刻を固定した変分では位置を固定すれば消える。
:::
::::

::::quiz{answer="1" tol="1e-6"}
滑らかな Hamilton 相流を初期値で微分した行列の行列式はいくらか。

:::explanation
発散が零であるため初期の恒等写像の行列式 $1$ を保つ。
:::
::::

::::quiz
時間に依存する $H$ で常に成り立つものを選べ。

- [ ] $H$ が一定
- [x] 正準相体積が一定
- [ ] 固定した点で密度が一定

:::explanation
体積保存には相変数の混合偏微分の交換だけを使う。
$\dot H=H_t$ は一般には零でなく、密度は軌道に沿って保存する。
:::
::::

::::quiz
移動する閉曲線に沿う積分 $\oint\sum_i p_i\dd q^i$ が保存される条件について正しいものを選べ。

- [ ] 閉曲線は単一の周期軌道でなければならない
- [x] 閉曲線上の各初期点を同じ Hamilton 相流で運べばよい
- [ ] 閉曲線に沿って $H=0$ でなければならない

:::explanation
積分の時間微分は曲線のパラメータ $u$ に関する $-\int\partial_uH\,\dd u$ に等しい。
閉曲線なのでその両端で $H$ の値が一致し、軌道が周期的かどうかには依存しない。
:::
::::

::::quiz{answer="0.5" tol="1e-6"}
単位をそろえた正準座標で、初期値の $q$ 方向の幅を二倍にし、$p$ 方向の幅をある倍率で伸縮する。
一自由度 Hamilton 相流による軸に沿う伸縮なら、$p$ 方向の倍率はいくらか。

:::explanation
行列式は二つの倍率の積で $1$。
よって $2a=1$ を解いて $a=0.5$ を得る。
:::
::::

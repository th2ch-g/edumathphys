---
title: 対称性と保存則
description: 作用の境界項と端点の移動を追い、時間・空間の対称性から保存量を求める。
notation:
  - symbol: 'p_i'
    meaning: 座標 qⁱ に共役な運動量
  - symbol: 'E_L'
    meaning: エネルギー関数
  - symbol: '\tau,\xi^i'
    meaning: 時刻と座標の無限小変化の係数
  - symbol: 'B'
    meaning: 作用の変化を表す境界関数
---

## この章の目標

- 循環座標を選び、保存量を使って自由度を減らす。
- 時間に依存する拘束でエネルギー関数の意味を判定する。
- 時刻と座標を同時に変換する Noether の定理を導く。
- 並進、回転、等速系の変更に対応する保存量を計算する。

## 循環座標と一般化運動量

球の表面を滑る小球を、鉛直軸のまわりで少し回してから放す。
重力も球面もこの回転で変化しない。
したがって運動の式には、経度の原点を選んだ痕跡が現れない。
このように座標の選び方によって保存則を式から読み取れる。

### 座標を一つ消す

:::definition[循環座標]{#def-am-cyclic}
独立な一般化座標の一つ $s=q^a$ が $L(q,\dot q,t)$ に陽に現れず、$\partial L/\partial s=0$ が恒等的に成り立つとき、$s$ を :term[循環座標]{en="cyclic coordinate" yomi="じゅんかんざひょう"} という。
速度 $\dot s$ は現れてよい。
:::

:ref[def-am-lagrangian] の運動量を $p_i=\partial L/\partial\dot q^i$ と書く。非保存力がないとき、:eqref[eq-am-lagrange-equations] より

$$
\dot p_a=\frac{\partial L}{\partial s}=0.
$$

逆に、座標の定数移動 $s\mapsto s+\eps$ がすべての状態で $L$ を不変にするなら、$\eps=0$ での微分は $\partial_s L=0$ である。
一つの特殊解の上で $\dot p_a=0$ を観測しただけでは、この恒等式は得られない。

:::example[球面上の振り子]
糸の長さ $l$、質量 $m$ とし、下向き鉛直からの角を $\theta$、方位角を $\phi$ とする。
水平面への射影の半径は $l\sin\theta$ だから、

$$
L=\tfrac12 ml^2(\dot\theta^2+\sin^2\theta\,\dot\phi^2)+mgl\cos\theta.
$$

$\phi$ が循環座標なので $p_\phi=ml^2\sin^2\theta\,\dot\phi$ は一定である。
$m=1\,\mathrm{kg}$、$l=1\,\mathrm m$、$\theta=\pi/6$、$\dot\phi=4\,\mathrm{s^{-1}}$ なら $p_\phi=1\,\mathrm{kg\,m^2/s}$ である。
同じ解が赤道を通るときには $\dot\phi=1\,\mathrm{s^{-1}}$ になる。
:::

```tikz caption="球面振り子を鉛直軸を含む平面で切った図。方位方向に動く円の半径は糸の長さより小さい。"
\coordinate (O) at (0,2.7);
\coordinate (P) at (1.35,0.36);
\draw[dashed] (0,-.15) -- (0,2.95);
\draw[thick] (O) -- node[right] {$l$} (P);
\draw[<->] (0,.36) -- node[below] {$l\sin\theta$} (P);
\draw (0,2.06) arc[start angle=-90,end angle=-60,radius=.64];
\node at (.26,1.88) {$\theta$};
\fill (O) circle (1.5pt);
\fill (P) circle (2.1pt);
```

### 消去後の変分で何を固定するか

循環座標の運動量を $c$ に固定し、$L_{\dot s\dot s}\ne0$ とする。
$p_s=c$ を解いて $\dot s=v(q^\alpha,\dot q^\alpha,c,t)$ と書ける。
ここで $\alpha$ は残った座標の番号である。
縮約した :term[Routh 関数]{en="Routhian" yomi="らうすかんすう"} を

$$
R(q^\alpha,\dot q^\alpha,c,t)=(L-c\dot s)_{\dot s=v}
$$

ととる。
$c$ を固定して $q^\alpha$ で微分すると、$v$ を通じた寄与は $(L_{\dot s}-c)\partial_\alpha v=0$ である。
$\dot q^\alpha$ で微分しても同じ相殺が起きる。
よって $R$ に通常の Lagrange 方程式を適用できる。
最後に $\dot s=v$ を積分すれば、消去した座標も復元できる。

:::example[回転の遠心項]
平面極座標で $L=m(\dot r^2+r^2\dot\phi^2)/2-V(r)$ とすると、$c=mr^2\dot\phi$ である。
$R=m\dot r^2/2-V-c^2/(2mr^2)$ より $m\ddot r=-V'+c^2/(mr^3)$ を得る。
$L$ に $\dot\phi=c/(mr^2)$ を直接代入すると最後の項の符号が逆になる。
固定した $c$ のもとで変分するには $-c\dot\phi$ が必要である。
:::

### 対称な方向を座標に選ぶ

並進に対してポテンシャルが一定でも、もとの座標のどれかが見た目で循環するとは限らない。
例えば斜面上の座標では、等高線に沿った向きを座標軸にとると並進対称性が初めて見える。

:::example[斜めの溝に沿う運動量]
水平面に質量 $m$ の小球を置き、

$$
V(x,y)=\frac{k}{2}(x-2y)^2
$$

という溝を作る。
$V_x=k(x-2y)$、$V_y=-2k(x-2y)$ だから、どちらの偏微分も恒等的には零でない。
ただし $(x,y)\mapsto(x+2s,y+s)$ は $x-2y$ を保つ。
そこで直交座標

$$
u=\frac{2x+y}{\sqrt5},\qquad v=\frac{x-2y}{\sqrt5}
$$

を選ぶ。
運動エネルギーは $m(\dot u^2+\dot v^2)/2$、位置エネルギーは $5kv^2/2$ になる。
従って

$$
p_u=m\dot u=(2p_x+p_y)/\sqrt5,\qquad \ddot v=-5kv/m .
$$

溝に沿う並進と、溝を横切る振動が同時に起こる。
例えば $m=1\,\mathrm{kg}, k=5\,\mathrm{N/m}$ なら振動の角振動数は $5\,\mathrm{s^{-1}}$ である。
$p_x,p_y$ を別々に一定としてしまうとこの横振動を失う。
:::

:::example[ねじに沿う対称性]
円筒面上で角と高さを $\phi,z$、半径を $a$ とする。
位置エネルギー $V=U(z-h\phi)$ を角をほどいた領域で考える。
円筒全体で一価にするには、$h\ne0$ のとき $U$ に周期 $2\pi|h|$ を課す。$h=0$ なら $U$ は高さだけの関数でよい。
角を $\eps$ だけ回し高さを $h\eps$ 上げると作用は変わらない。
$w=z-h\phi$ に替えると

$$
L=\tfrac m2[a^2\dot\phi^2+(\dot w+h\dot\phi)^2]-U(w).
$$

従って新しい $\phi$ の共役運動量は $ma^2\dot\phi+hm\dot z$、すなわち古い $p_\phi+h p_z$ である。
$h$ は角度あたりの上昇量で長さの単位をもち、二つの項の単位はともに角運動量になる。
:::

### 縮約で得る初期条件

循環座標を一つ除くと、二階方程式の本数が一つ減る。
除いた座標の二つの初期条件は、運動量の定数 $c$ と、最後の積分に現れる初期位置 $s(0)$ として残る。
従って縮約した解を比較するときには、同じ $c$ を固定した比較なのかも指定する。

例えば同じ半径で静かに動径運動を始めても、$c$ が違えば $c^2/(2mr^2)$ の高さが異なる。
$m=1\,\mathrm{kg},r=2\,\mathrm m$ で $c=2$ と $4\,\mathrm{kg\,m^2/s}$ を比較すると、遠心項のエネルギーは $0.5$ と $2\,\mathrm J$ である。
この差は、縮約前には円周方向の運動エネルギーの差であった。

### いくつかの運動量を同時に固定する

二つの車輪が回転できる台車を考える。回転角を測る目盛の原点はどちらの車輪にも任意に置ける。二つの循環速度を消すとき、対角成分だけを拾ってよいのは慣性の交差項が消える場合に限る。一般に残す座標を列 $x$、循環座標を列 $s$ とし、速度について二次の系を

$$
L=\tfrac12\dot x^\top A\dot x+\dot x^\top B\dot s+
\tfrac12\dot s^\top C\dot s+a^\top\dot x+b^\top\dot s-V
$$

と書こう。行列積の各成分は行と列の積の和である。係数は $x,t$ に依存してよい。質量行列が正定値なら主小行列 $C$ は正定値である。したがって $c=B^\top\dot x+C\dot s+b$ から $\dot s=C^{-1}(c-b-B^\top\dot x)$ が一意に決まる。$d=c-b$ とおいて平方を完成すると

$$
R=\tfrac12\dot x^\top(A-BC^{-1}B^\top)\dot x+
(a+BC^{-1}d)^\top\dot x-V-\tfrac12 d^\top C^{-1}d .
\label{eq-am-block-routh}
$$

残る運動の慣性は $A$ より小さくなり得る。これを確かめるには、任意の $\dot x$ を固定して運動エネルギーの二次部分を $\dot s$ について最小化すればよい。最小値は上式の最初の項に等しい。元の二次形式が正定値なので、$\dot x\ne0$ なら最小値も正である。速度に一次の項は一般には全微分にならない。二つの残存座標があれば、その係数の交差微分の差が力の速度比例項として残る。

消した運動は $-\partial R/\partial c=\dot s$ から戻る。従って $x$ が一周して出発点に戻っても $s$ は元の値に戻るとは限らない。例えば角運動量を零とし、$b=0$ で $B,C$ が時刻に依存しなければ $\Delta s=-\oint C^{-1}B^\top\dd x$ である。被積分量が一価の関数の全微分ならどの閉路でも零である。被積分量が定義された領域内で閉路が曲面を囲む場合には、Stokes の定理によりその外微分の面積分でずれを求められる。この条件では閉路をたどる速さにずれは依存しない。

:::example[台とばねの相対変位]
水平な床を台が滑り、その上を質点が滑る。台と質点はばねで結ばれている。両者の質量を $3,2\,\mathrm{kg}$、台の位置を $s$、ばねの自然長からの伸びを $q$、ばね定数を $6\,\mathrm{N/m}$ とする。質点の実験室速度は $w+u$ である。
ここでは式中の長さ・時刻・質量をそれぞれ m、s、kg で測った数値で書く。$u=\dot q,\ w=\dot s$ とすれば、運動エネルギーを展開した結果と共役運動量は

$$
L=- \frac{6 q^{2} - 2 u^{2} - 4 u w - 5 w^{2}}{2},\qquad p_s=2 u + 5 w.
$$

速度を消す際には $p_s=c$ を固定する。残る座標で変分する前に $cw$ を引き、項を速度の次数ごとに整理すると

$$
R=\frac12\left(\frac{6}{5}\right)u^2+\left(\frac{2 c}{5}\right)u-\left(\frac{c^{2} + 30 q^{2}}{10}\right).
$$

$q$ 一つだけの問題で、$u$ に一次の係数は $q$ だけの関数である。従って局所的にはその原始関数の時間微分となり、運動方程式には寄与しない。$M(q)=R_{uu}$ とすれば $M\ddot q=-\partial_q U_{\rm eff}-M'u^2/2$ を得る。この装置では

$$
\ddot q=- 5 q,\qquad \dot s=- \frac{- c + 2 u}{5}.
$$

初期値を $(q,u,w)=(1,1,-1)$ とすると、保存される運動量の数値は $c=-3$、初期の加速度は $\ddot q=-5$ となる。エネルギーは縮約前の $T+V$ でも縮約後の $uR_u-R$ でも $\frac{9}{2}\,\mathrm J$ である。このとき負の全運動量は重心が左へ進むことを示す。ばねが縮み始める時刻を知るには、重心の初期速度を相対運動の初期速度と取り違えないようにする。
:::

:::example[円錐面上の小環]
鉛直上向きに開いた滑らかな円錐を $z=q,\ r=q>0$ とし、方位角を $s$ とする。質量 $1\,\mathrm{kg}$、重力加速度 $10\,\mathrm{m/s^2}$ とする。半径と高さが同時に変わるため、子午線方向の速さの二乗は $2u^2$ となる。
ここでは式中の長さ・時刻・質量をそれぞれ m、s、kg で測った数値で書く。$u=\dot q,\ w=\dot s$ とすれば、運動エネルギーを展開した結果と共役運動量は

$$
L=\frac{q^{2} w^{2} - 20 q + 2 u^{2}}{2},\qquad p_s=q^{2} w.
$$

速度を消す際には $p_s=c$ を固定する。残る座標で変分する前に $cw$ を引き、項を速度の次数ごとに整理すると

$$
R=\frac12\left(2\right)u^2+\left(0\right)u-\left(\frac{c^{2} + 20 q^{3}}{2 q^{2}}\right).
$$

$q$ 一つだけの問題で、$u$ に一次の係数は $q$ だけの関数である。$M(q)=R_{uu}$ とすれば $M\ddot q=-\partial_q U_{\rm eff}-M'u^2/2$ を得る。この装置では

$$
\ddot q=- \frac{- c^{2} + 10 q^{3}}{2 q^{3}},\qquad \dot s=\frac{c}{q^{2}}.
$$

初期値を $(q,u,w)=(2,0,2)$ とすると、保存される運動量の数値は $c=8$、初期の加速度は $\ddot q=-1$ となる。エネルギーは縮約前の $T+V$ でも縮約後の $uR_u-R$ でも $28\,\mathrm J$ である。初期には重力の斜面方向成分が勝つ。円軌道には $c^2=10q^3$ を要し、頂点に近づくと同じ非零の $c$ を保つなら方位角速度の絶対値は増す。ただし有限エネルギーの軌道は有効ポテンシャルの発散により頂点まで届かない。
:::

:::example[回転台に沿って動く質点]
水平な回転台の慣性モーメントを $2\,\mathrm{kg\,m^2}$、半径方向の溝を滑る質点の質量を $1\,\mathrm{kg}$ とする。$q$ を中心からの距離、$s$ を台の角とする。質点を原点に向けてばね定数 $4\,\mathrm{N/m}$ で引く。軸受はトルクを伝えない。
ここでは式中の長さ・時刻・質量をそれぞれ m、s、kg で測った数値で書く。$u=\dot q,\ w=\dot s$ とすれば、運動エネルギーを展開した結果と共役運動量は

$$
L=\frac{q^{2} w^{2} - 4 q^{2} + u^{2} + 2 w^{2}}{2},\qquad p_s=w \left(q^{2} + 2\right).
$$

速度を消す際には $p_s=c$ を固定する。残る座標で変分する前に $cw$ を引き、項を速度の次数ごとに整理すると

$$
R=\frac12\left(1\right)u^2+\left(0\right)u-\left(\frac{c^{2} + 4 q^{4} + 8 q^{2}}{2 \left(q^{2} + 2\right)}\right).
$$

$q$ 一つだけの問題で、$u$ に一次の係数は $q$ だけの関数である。$M(q)=R_{uu}$ とすれば $M\ddot q=-\partial_q U_{\rm eff}-M'u^2/2$ を得る。この装置では

$$
\ddot q=- \frac{q \left(- c + 2 q^{2} + 4\right) \left(c + 2 q^{2} + 4\right)}{\left(q^{2} + 2\right)^{2}},\qquad \dot s=\frac{c}{q^{2} + 2}.
$$

初期値を $(q,u,w)=(1,0,2)$ とすると、保存される運動量の数値は $c=6$、初期の加速度は $\ddot q=0$ となる。エネルギーは縮約前の $T+V$ でも縮約後の $uR_u-R$ でも $8\,\mathrm J$ である。この初期状態では動径加速度が零である。台をモーターで一定速に保つ別の実験と同じ初期加速度になっても、半径を乱した後の角速度の変化は異なる。
:::


### 縮約後の速度一次の項が力になる例

一つの循環角 $s$ をもち、残りの座標が二つあるときに :eqref[eq-am-block-routh] の一次項を実際に計算する。局所模型
$$
L=\tfrac m2(\dot x^2+\dot y^2)+\tfrac C2(\dot s+\beta x\dot y)^2
$$
で $x,y$ は長さ、$C$ は慣性モーメント、$\beta$ の単位は $\mathrm{m^{-2}}$ とする。角運動量 $c=C(\dot s+\beta x\dot y)$ を固定すると $\dot s=c/C-\beta x\dot y$。$R=L-c\dot s$ に戻せば $R=m(\dot x^2+\dot y^2)/2+c\beta x\dot y-c^2/(2C)$ である。一次形式 $c\beta x\,\dd y$ の外微分は $c\beta\,\dd x\wedge\dd y$ であり、$c\beta\ne0$ では全微分として消せない。実際に二つの Lagrange 方程式は $m\ddot x=c\beta\dot y,\ m\ddot y=-c\beta\dot x$ となる。縮約で残った力は速度に垂直で、縮約前のエネルギー $m(\dot x^2+\dot y^2)/2+c^2/(2C)$ を保つ。

例えば $m=2\,\mathrm{kg},c=3\,\mathrm{kg\,m^2/s},\beta=2\,\mathrm{m^{-2}}$ の場合、初速度が $(1,0)\,\mathrm{m/s}$ なら初期加速度は $(0,-3)\,\mathrm{m/s^2}$。実験室平面での軌跡は角速度の大きさ $3\,\mathrm{s^{-1}}$、半径 $1/3\,\mathrm m$ の円になる。右向きに出た点は下へ曲がる。一次項を直接捨てると等速直線運動になり、エネルギーの値の検算だけではこの誤りを検出できない。

### 角運動量が零でも回転できる機構

軸受が鉛直軸まわりのトルクを及ぼさない台の上で、内部の小球を動かす。左右対称に二球を置けば重心は軸上にとどまる。台の慣性モーメントを $I_0$、二球合わせた質量を $m$、軸からの距離を $r$ とする。台が実験室に対して回った角度を $s$、台から見た小球の方位角を $\alpha$ とする。内部の機構は $r,\alpha$ を指定できるが、外からの軸トルクはない。このとき全角運動量は

$$
p_s=I_0\dot s+mr^2(\dot s+\dot\alpha).
$$

初め静止していれば $p_s=0$ である。従って形を変える間の台の角度は $\dd s=-h(r)\dd\alpha$、$h(r)=mr^2/(I_0+mr^2)$ から求まる。$r$ を固定して内部角を $\alpha_0$ 増やし、半径を $r_1$ から $r_2$ に変え、角を $\alpha_0$ 戻し、半径も戻す。この順序なら台の最終角は $[h(r_2)-h(r_1)]\alpha_0$ である。

例えば $I_0=2\,\mathrm{kg\,m^2},m=2\,\mathrm{kg},r_1=1\,\mathrm m,r_2=2\,\mathrm m,\alpha_0=\pi/2$ なら $h(r_1)=1/2,h(r_2)=4/5$ であり、台は $3\pi/20$、すなわち $27^\circ$ 回る。二球を半径方向だけに移す区間では角度は変わらない。周回を逆に行うとずれの符号が反転する。一周後に形と速さを元に戻しても姿勢まで元に戻さなくてよいのは、角運動量が角度そのものを拘束していないためである。

図の縦辺をゆっくりたどっても速くたどっても、そこでは $\dd\alpha=0$ のため台は回らない。横辺を同じ量だけ往復して回転を取り消すには、往路と復路で同じ半径を保つ必要がある。小さな長方形なら角度のずれは $h'(r)\,\Delta r\,\Delta\alpha$ であり、$h'=2mI_0r/(I_0+mr^2)^2$。したがって球を軸のごく近くに置いた場合には、半径を少し変えても回転の効果は小さい。

```tikz caption="内部角と半径の閉路。台の回転は、横の二辺で異なる慣性モーメントを使うことから生じる。"
\draw[->] (0,0)--(4.5,0) node[right] {$\alpha$};
\draw[->] (0,0)--(0,3.1) node[above] {$r$};
\draw[dashed] (0,.8)--(.65,.8) (0,2.35)--(.65,2.35);
\node[left] at (0,.8) {$r_1$};
\node[left] at (0,2.35) {$r_2$};
\draw[thick,->] (.7,.8)--(3.55,.8);
\draw[thick,->] (3.55,.8)--(3.55,2.35);
\draw[thick,->] (3.55,2.35)--(.7,2.35);
\draw[thick,->] (.7,2.35)--(.7,.8);
\node[below] at (2.1,.8) {$+\alpha_0$};
\node[above] at (2.1,2.35) {$-\alpha_0$};
```

## エネルギー関数とその保存

時計の原点をずらしても同じ実験ができるなら、運動に沿って一定となる時間方向の量がある。
一般化座標ではその量と実験室での力学的エネルギーの関係を調べる必要がある。

:::definition[エネルギー関数]{#def-am-energy-function}
速度の関数として

$$
E_L(q,\dot q,t)=\sum_{i=1}^n p_i\dot q^i-L
\label{eq-am-energy-function}
$$

を :term[エネルギー関数]{en="energy function" yomi="えねるぎーかんすう"} という。
各項の単位は $\mathrm J$ である。
:::

### 源のあるエネルギー収支

外から加わる一般化力を $Q_i$ として $\dot p_i=L_{q^i}+Q_i$ と書く。
積の微分を展開すると

$$
\begin{aligned}
\dot E_L
&=\sum_{i=1}^n(\dot p_i\dot q^i+p_i\ddot q^i)
 -\sum_{i=1}^n(L_{q^i}\dot q^i+L_{\dot q^i}\ddot q^i)-L_t\\
&=\sum_{i=1}^nQ_i\dot q^i-L_t .
\end{aligned}
\label{eq-am-energy-balance}
$$

特に $Q_i=0$ なら、自律系（$L_t=0$）で $E_L$ が保存する。
これは :ref[prop-am-beltrami] を時間に適用した結果でもある。
$Q_i$ が速度に線形な抵抗なら :ref[def-am-rayleigh] の二次形式 $\mathcal F$ により仕事率は $-2\mathcal F$ となる。

### 力学的エネルギーとの比較

$T$ が速度について二次の斉次式なら、$T(q,\lambda\dot q)=\lambda^2T(q,\dot q)$ を $\lambda=1$ で微分して $\sum_i\dot q^i T_{\dot q^i}=2T$ を得る。
したがって $L=T-V(q)$ では $E_L=T+V$ である。
より一般に $L=T_2+T_1+T_0-V$ と次数別に書くと $E_L=T_2-T_0+V$ である。

:::example[回転する水平な棒]
角速度 $\Omega$ で駆動する水平な棒に小球を通す。
棒に沿う座標 $r$ では $L=m(\dot r^2+\Omega^2r^2)/2$ である。
一定の $\Omega$ では $E_L=m(\dot r^2-\Omega^2r^2)/2$ が保存する。
一方、実験室の運動エネルギーを微分し $\ddot r=\Omega^2r$ を使うと $\dot T=2m\Omega^2r\dot r$ を得る。
駆動装置が棒の反力を通してこの仕事をする。
:::

:::example[ばね定数を変える]
$L=m\dot x^2/2-k(t)x^2/2$ では $\dot E_L=\dot k x^2/2$ である。
$x=0.10\,\mathrm m$、$\dot k=20\,\mathrm{N/(m\,s)}$ の瞬間、装置は $0.10\,\mathrm W$ を供給する。
ばねの自然長を通る瞬間には、この仕事率は零になる。
:::

### 観測者と一緒に移動するポテンシャル

一定速度 $u$ で右へ動くばねの支持台を考える。
実験室の座標 $x$ を使えば

$$
L=\tfrac12m\dot x^2-\tfrac12 k(x-ut)^2.
$$

ポテンシャルの谷が移動するから、実験室のエネルギーの時間微分は $-ku(x-ut)$ である。
では支持台に乗って観測すると何が一定か。
$s=x-ut$ と置き換えてみると、

$$
L'=\tfrac12m(\dot s+u)^2-\tfrac12ks^2,\qquad
E_{L'}=\tfrac12m\dot s^2+\tfrac12ks^2-\tfrac12mu^2 .
$$

新しい座標では $L'_t=0$ なのでこの量が一定になる。
元の変数で展開すれば $E_{L'}=E_L-up_x$。
実際、$\dot E_L=-ku(x-ut)$ と $\dot p_x=-k(x-ut)$ の右辺の比は $u$ である。
時刻を $\eps$ だけ遅らせると同時に場所を $u\eps$ だけ進めて観測すれば装置は元通りに見える。
この結合した対称性を次節で系統的に取り出す。

### 速度依存ポテンシャルからエネルギーを読み取る

電荷 $e$ の粒子で $L=m|\vb v|^2/2+e\vb A\cdot\vb v-e\Phi$ とすると、

$$
\vb p=m\vb v+e\vb A,\qquad E_L=\tfrac12m|\vb v|^2+e\Phi .
$$

磁場に関する項は速度に一次なので打ち消される。
静電場では $\vb E=-\nabla\Phi$、磁気力は速度に垂直であるから

$$
\frac{\dd}{\dd t}\left(\tfrac12m|\vb v|^2\right)=e\vb E\cdot\vb v=-e\frac{\dd\Phi}{\dd t}
$$

となり、力学からも保存を確認できる。
時間依存ゲージを使うと $E_L$ の式は変わる。
その場合の保存量は作用の境界項を含めて求める。

### 動く拘束が渡すエネルギーと運動量

動く壁に球を当てると、衝突が弾性的でも実験室での運動エネルギーは変わる。滑らかな拘束についても、このエネルギー交換を反力から読み取れる。Cartesian 座標 $\vb r$ に $f(\vb r,t)=0$ を課し、反力を :ref[prop-am-ideal-multiplier] の規約 $\vb R=+\lambda\nabla f$ で表す。実際の速度では $\nabla f\cdot\dot{\vb r}=-f_t$ であるから、反力の仕事率は $-\lambda f_t$ である。仮想変位に対しては時計を止めるので、この仕事率が非零でも仮想仕事は零でよい。複数の独立な拘束でも各反力の仕事率を加えればよい。

壁が法線方向に一定速度 $U$ で動く場合、滑らかな急峻なポテンシャル $V(x-Ut)$ でまず考える。力を $F=-V'$ とおけば $\dot p=F$、$\dot E=UF$ である。壁の反力が大きく働く時間を積分すると、入射と反射の自由区間を結ぶ $\Delta E=U\Delta p$ が得られる。硬い壁の極限で壁の系の相対速度は反転するので $v_{\rm out}=2U-v_{\rm in}$。これは $E-Up$ が等しい二状態のうち、相対速度が非零の入射から反射へ移る方の解である。

:::example[後退するピストンへの衝突]
質量 $0.20\,\mathrm{kg}$ の滑車用のおもりが右向きに $5\,\mathrm{m/s}$ で進み、右へ $1\,\mathrm{m/s}$ で動く壁に当たる一次元模型を考える。反射後は $-3\,\mathrm{m/s}$、力積は $0.20(-3-5)=-1.6\,\mathrm{N\,s}$ となる。粒子は運動エネルギーを $0.20(25-9)/2=1.6\,\mathrm J$ 失う。壁を一定速に保つ装置が引き受けるエネルギーまで含めれば収支が閉じる。壁の質量も有限で自由に動くなら $U$ は衝突中に変わるので、ここで用いた時空並進は粒子単独の対称性にならない。
:::

### 電動機を系に含めた場合と角速度を指定した場合

回転軸に摩擦がなくても、角速度を指定するためにトルクが要ることがある。その大きさを計算して系の境界を定めよう。$r$ と $\phi$ を独立とした段階で $L=m\dot r^2/2+I(r)\dot\phi^2/2-U(r)$ とする。回転を $\phi=\Phi(t)$ と指定する拘束 $f=\phi-\Phi(t)=0$ に対して、:ref[prop-am-ideal-multiplier] の符号では $\lambda=\dot p_\phi=\dd(I\dot\Phi)/\dd t$ である。$\lambda$ の単位はトルクであり、装置の仕事率は $-\lambda f_t=\lambda\dot\Phi$。$\Omega=\dot\Phi$ と書けば

$$
\begin{aligned}
m\ddot r&=\tfrac12 I'(r)\Omega^2-U'(r),\\
\frac{\dd E_{\rm lab}}{\dd t}
&=\left[I'(r)\dot r\Omega+I(r)\dot\Omega\right]\Omega,\qquad
E_{\rm lab}=\tfrac12m\dot r^2+\tfrac12 I\Omega^2+U .
\end{aligned}
$$

一方、$\phi=\Phi(t)$ を代入してから作ったエネルギー関数は $E_{\rm red}=E_{\rm lab}-I\Omega^2$。差の微分を取ると $\dot E_{\rm red}=-I\Omega\dot\Omega$ であり、縮約後の $-L_t$ と等しい。角速度が一定なら $E_{\rm lab}-\Omega p_\phi$ が保存する。これは回転と時間並進を組み合わせた保存量である。

:::example[モーターのトルクを二つの測定値から分ける]
$1\,\mathrm{kg}$ の小球を慣性モーメント $2\,\mathrm{kg\,m^2}$ の台に乗せる。$r=1\,\mathrm m,\dot r=0.50\,\mathrm{m/s}$ の瞬間、$I=3\,\mathrm{kg\,m^2},\dot I=1\,\mathrm{kg\,m^2/s}$ である。$\Omega=2\,\mathrm{s^{-1}},\dot\Omega=-1\,\mathrm{s^{-2}}$ ならモーターのトルクは $1\times2+3(-1)=-1\,\mathrm{N\,m}$ であり、仕事率は $-2\,\mathrm W$。縮約したエネルギー関数は同時に $6\,\mathrm W$ の割合で増える。この相違は両者の差 $I\Omega^2$ が $4-12=-8\,\mathrm W$ の割合で減ることから説明できる。仕事の向きを測りたいときは軸トルクと回転速度を測る。
:::

## Noether の定理

### 同じ時刻で比べる変化と端点の移動

変換前の経路上の点 $(t,q)$ が、変換後に $(t',q')$ になるとする。
微小な群パラメータ $\eps$ に対して

$$
t'=t+\eps\tau(q,t),\qquad
q'^i(t')=q^i(t)+\eps\xi^i(q,t)
$$

と一次まで書く。
以下、$D$ は元の経路に沿う全時間微分である。
同じ時計の読みで経路を比べた変化は $\eta^i=\xi^i-\dot q^i\tau$ である。
分母の時間間隔も変わるため、変換した点での速度は

$$
\frac{\dd q'^i}{\dd t'}=\dot q^i+\eps(D\xi^i-\dot q^iD\tau)+O(\eps^2)
$$

となる。
:ref[def-am-virtual-displacement] で固定していた時刻が、ここでは変換の対象に含まれる。

:::theorem[境界項を含む Noether の定理]{#thm-am-noether}
$C^2$ 級の $L$ と滑らかな生成係数 $\tau,\xi^i$ を考える。上の点変換が任意の経路に沿って一次まで

$$
L(q',\dot q',t')\,\dd t'= (L+\eps DB)\,\dd t+O(\eps^2)\,\dd t
$$

を満たすとする。
$B(q,t)$ は一つの境界関数である。
すると外力項のない運動方程式の解に沿って

$$
J=\sum_{i=1}^np_i\xi^i-E_L\tau-B
\label{eq-am-noether-charge}
$$

が保存する。
この $J$ を :term[Noether 量]{en="Noether charge" yomi="ねーたーりょう"} という。
:::

:::proof
時間間隔の変化を含めて $\eps$ の係数を $\Delta$ と書くと、

$$
\Delta=L_t\tau+\sum_{i=1}^n\bigl(L_{q^i}\xi^i+p_iD\xi^i-p_i\dot q^iD\tau\bigr)+LD\tau.
$$

$\xi^i=\eta^i+\dot q^i\tau$ を入れ、$DL=L_t+\sum_i(L_{q^i}\dot q^i+p_i\ddot q^i)$ を用いる。
積の微分を集めれば解の上に限らず

$$
\Delta-DB=
\sum_{i=1}^n(L_{q^i}-Dp_i)\eta^i
+D\left(\sum_{i=1}^np_i\eta^i+L\tau-B\right).
$$

左辺は仮定により零であり、解に沿えば右辺の和も零である。
括弧内は :eqref[eq-am-energy-function] により $J$ に等しいので $DJ=0$ である。
:::

### 境界項が観測量にどう効くか

$L$ を $L+DF$ に取り替えると :eqref[eq-am-total-derivative] より同じ軌道が得られる。
このとき $p_i$ は $F_{q^i}$ だけ増え、$E_L$ は $F_t$ だけ減る。
変換による端点での $F$ の変化は $\tau F_t+\sum_i\xi^iF_{q^i}$ だから、$B$ にもこれを加える。
式 :eqref[eq-am-noether-charge] で三つの変化は打ち消し合う。

:::example[等速で移動する慣性系]
自由粒子の能動的な変換 $x'=x+\eps t$ を考える。
$\Delta=m\dot x=D(mx)$ なので $B=mx$ であり、$J=pt-mx$ となる。
$m=2\,\mathrm{kg}$、$x(0)=3\,\mathrm m$ ならどの速度でも $J=-6\,\mathrm{kg\,m}$ である。
厳密な $L$ の不変性だけを調べると、この保存則を見落とす。
:::

:::example[時間と長さを同時に変える]
$x>0$ で $L=m\dot x^2/2-g/x^2$ とする。
$x'=e^\eps x,\ t'=e^{2\eps}t$ では速度は $e^{-\eps}$ 倍になり、$L\,\dd t$ は不変である。
$\xi=x,\ \tau=2t$ より $J=px-2tE_L$ である。
直接微分しても $x\dot p+p\dot x=2g/x^2+m\dot x^2=2E_L$ となる。
:::

### 運動の一部だけを取り出した場合

解に沿った変化と任意の試験経路での変化を区別する。
例えば $J=p_x$ がある特別な運動で一定であっても、すべての試験経路に対して並進が作用の境界項だけを変えるとは限らない。
調和振動子が原点に静止する解では $p_x=0$ が続くが、原点からずらした試験経路で $L$ の並進変化は $-kx$ となる。

保存系に外力 $Q_i$ を加えたときは、もとの $L$ の対称性を用いて収支を記述できる。
証明の最後で $Dp_i=L_{q^i}+Q_i$ を残すと

$$
DJ=\sum_{i=1}^n Q_i(\xi^i-\dot q^i\tau).
$$

空間並進なら右辺はその方向の外力であり、時間並進なら $J=-E_L$ に対応して右辺は仕事率の負になる。
例えば平面回転で $Q$ が粘性抵抗 $-b\dot{\vb r}$ なら $\dot L_z=-(b/m)L_z$、従って $L_z(t)=L_z(0)e^{-bt/m}$ と減る。
対称な形をした抵抗でも、作用から導く閉じた保存系の仮定は外れている。

:::example[一様磁場で保存する運動量]
$\vb B=B\vb e_z$、電場零とし、平面で $\vb A=(0,Bx,0)$ を選ぶ。
粒子の Lagrangian は

$$
L=\tfrac m2(\dot x^2+\dot y^2)+eBx\dot y .
$$

$y$ は循環し $p_y=m\dot y+eBx$ が保存する。
$x$ の並進では $\Delta=eB\dot y=D(eBy)$ であるから、もう一つの保存量は $m\dot x-eBy$ である。
$eB\ne0$ とすると

$$
X_c=x+\frac{m\dot y}{eB},\qquad
Y_c=y-\frac{m\dot x}{eB}
$$

が一定である。
これらは円運動の中心の座標であり、Lorentz 力 $m\ddot x=eB\dot y,\ m\ddot y=-eB\dot x$ を代入して確認できる。
ゲージを変えてもこの中心は動かない。
:::

:::example[時間依存の座標変換による積分]
調和振動子で $\xi=f(t)$、$\tau=0$ とする。
$f$ 自身が $\ddot f+\omega^2f=0$ を満たせば、

$$
\Delta=m\dot q\dot f-m\omega^2 qf=D(mq\dot f).
$$

得られる $J=m(\dot qf-q\dot f)$ は二つの振動解の Wronskian である。
$f=\cos\omega t$ と $\sin\omega t$ を選ぶとそれぞれ $m\dot q(0)$ と $-m\omega q(0)$ を与える。
時間に陽に依存する式でも初期条件を保存できる。
:::

### 候補となる変換を方程式から探す

対称性が図から明らかでない場合、任意の試験経路で成り立つという仮定を使って $\xi,\tau,B$ を逆算できる。
一自由度の自然系

$$
L=\frac m2 v^2-V(q,t),\qquad v=\dot q
$$

で試す。
$\tau,\xi,B$ を $q,t$ の関数として、$D=\partial_t+v\partial_q$ を Noether の条件に入れる。
速度の変化と積分要素の変化を一緒に展開すると

$$
\begin{aligned}
\Delta={}&-\tfrac m2\tau_q v^3
+m\bigl(\xi_q-\tfrac12\tau_t\bigr)v^2\\
&+(m\xi_t-V\tau_q)v
-\xi V_q-\tau V_t-V\tau_t .
\end{aligned}
$$

右辺は $DB=B_t+vB_q$ と等しくなければならない。
$q,t$ を固定しても試験速度 $v$ は任意だから、各次数の係数が等しい。
順に

$$
\begin{aligned}
\tau_q&=0,\\
\xi_q&=\tfrac12\tau_t,\\
B_q&=m\xi_t,\\
B_t&=-\xi V_q-\tau V_t-V\tau_t .
\end{aligned}
$$

を得る。
例えば $\xi^i,\tau$ に最初から定数だけを代入する方法では、等速変換や時間に依存する振動解の加算は見つからない。

最初の三条件は $V$ に依存しない。
$\tau=\tau(t)$ として積分すると、ある関数 $a(t),b(t)$ により

$$
\begin{aligned}
\xi(q,t)&=\tfrac12\dot\tau(t)q+a(t),\\
B(q,t)&=\tfrac m4\ddot\tau(t)q^2+m\dot a(t)q+b(t).
\end{aligned}
$$

残る最後の式に具体的な $V$ を入れると、許される時間依存が決まる。

:::example[調和振動子の変換を系統的に求める]
$V=m\omega^2q^2/2$ とする。
最後の条件を $q$ の二次、一次、定数に分ければ

$$
\dddot\tau+4\omega^2\dot\tau=0,\qquad
\ddot a+\omega^2a=0,\qquad \dot b=0
$$

となる。
$\tau=1$ は時間並進、$\tau=0$ では $a=\cos\omega t$ または $\sin\omega t$ が先の例を再現する。
さらに $\tau=\sin(2\omega t)/\omega$ と選ぶと $\xi=q\cos(2\omega t)$、

$$
B=-m\omega q^2\sin(2\omega t)
$$

である。
この変換の Noether 量を初期時刻で評価すると $J(0)=p(0)q(0)$。
運動の途中では

$$
J=pq\cos(2\omega t)+\left(m\omega q^2-\frac{E_L}{\omega}\right)\sin(2\omega t)
$$

という時間依存の式で同じ積を読み取れる。
:::

### 境界項を座標の選択で吸収する

時刻を動かさない一パラメータ対称性があり、その生成方向が零でない近傍を考える。
変換に沿う曲線を積分して、そのパラメータを $s$、曲線を区別する横座標を $u$ とする。
この座標では変換は $s\mapsto s+\epsilon$ となる。
境界関数が $B(u,s,t)$ なら

$$
F(u,s,t)=-\int_{s_0}^s B(u,\sigma,t)\,\dd\sigma
$$

を選べば、$L'=L+DF$ の境界関数は $B+\partial_sF=0$。
従って $s$ は $L'$ の循環座標である。
新しい運動量は $p'_s=p_s+\partial_sF=p_s-B$ であり、古い表示の Noether 量と一致する。

この構成は局所的である。
$s$ が周期座標なら、一周で $\int B\,\dd s$ が零でない場合、選んだ $B$ に対して上の $F$ は一価に戻らない。ただし厳密に周期的な $L$ が $\partial_sL=DB$ を満たすなら、周期積分を取って $D\int B\,\dd s=0$ を得る。任意の試験速度に対して成り立つので、その周期積分は横座標と時刻に依存しない。$B$ に定数を足して平均を零にすれば、角を一周したときの障害は消える。
保存量があるという局所的な事実と、一つの座標表現でどこでも循環座標にできるという主張にはこの差がある。

### 質量行列を保つ変換を見つける

一般化座標に移ると運動エネルギーには座標に依存する係数 $g_{ij}$ が現れる。ポテンシャルの等高線に沿って点を動かすだけでは、速度の長さまで保てるとは限らない。この二つの条件を分けて調べる。時間を変えない変換 $\delta q^i=\epsilon X^i(q)$ と自然系 $L=\sum_{i,j}g_{ij}\dot q^i\dot q^j/2-V$ を考える。任意の試験速度に対して一次の変化を展開すると

$$
\delta L/\epsilon=\frac12\sum_{i,j}D_{ij}\dot q^i\dot q^j-\sum_k X^k\partial_k V,\qquad
D_{ij}=\sum_k(X^k\partial_k g_{ij}+g_{kj}\partial_iX^k+g_{ik}\partial_jX^k).
$$

定数項と二次形式は独立である。$D$ は対称なので二次形式が恒等的に零であれば、単位ベクトルと二つの単位ベクトルの和を代入してすべての成分が零と分かる。従って厳密な対称性の条件は $X(V)=0$ および $D_{ij}=0$ である。後者を満たす方向場を :term[Killing ベクトル場]{en="Killing vector field" yomi="きりんぐべくとるば"} と呼ぶ。対応する保存量は $J_X=\sum_i p_iX^i$ である。一定質量の平面でも極座標なら $g_{\phi\phi}=mr^2$ が変化する。並進が存在することと行列の全成分が定数であることを同一視してはならない。

:::example[極座標で横向きの並進を検出する]
位置エネルギーは高さにだけ依存する。直交座標で横へ一歩移したベクトルを極座標の基底に分解すると、半径方向と方位方向の両方に成分が現れる。質量を $1\,\mathrm{kg}$ として SI の数値を用いる。$\delta r=\epsilon X^r,\ \delta\phi=\epsilon X^\phi$ と書いたとき、座標上で比較すべき関数は

$$
V=r^{2} \sin^{2}{\left(\phi \right)},\qquad X^r=\cos{\left(\phi \right)},\quad X^\phi=- \frac{\sin{\left(\phi \right)}}{r}.
$$

まず力の変化を調べる。$\partial_rV$ と $\partial_\phi V$ はそれぞれ $2 r \sin^{2}{\left(\phi \right)}$ と $2 r^{2} \sin{\left(\phi \right)} \cos{\left(\phi \right)}$ である。これらに $X^r,X^\phi$ を掛けて足すと零になる。次に線素 $\dd r^2+r^2\dd\phi^2$ を調べる。方向場の微分は

$$
\frac{\partial(X^r,X^\phi)}{\partial(r,\phi)}=\left(\begin{matrix}0 & - \sin{\left(\phi \right)}\\\frac{\sin{\left(\phi \right)}}{r^{2}} & - \frac{\cos{\left(\phi \right)}}{r}\end{matrix}\right).
$$

従って $D_{rr}=2\partial_r X^r$、$D_{r\phi}=\partial_\phi X^r+r^2\partial_rX^\phi$、$D_{\phi\phi}=2rX^r+2r^2\partial_\phi X^\phi$ はすべて消える。保存量を速度 $u=\dot r,w=\dot\phi$ で表示すると

$$
J_X=X^ru+r^2X^\phi w=- r w \sin{\left(\phi \right)} + u \cos{\left(\phi \right)}.
$$

例えば $(r,\phi,u,w)=(2,\frac{\pi}{3},1,2)$ でその数値は $- \frac{-1 + 4 \sqrt{3}}{2}$ である。原点 $r=0$ ではこの座標表示が使えないので、そこまで解を続ける場合には直交座標へ戻る。
:::

:::example[軸を原点からずらした回転]
回転の中心を直交座標の $(-2,0)$ に置く。中心からの相対座標は $(r\cos\phi+2,r\sin\phi)$ であり、ポテンシャルをこの長さの二乗に比例させる。質量を $1\,\mathrm{kg}$ として SI の数値を用いる。$\delta r=\epsilon X^r,\ \delta\phi=\epsilon X^\phi$ と書いたとき、座標上で比較すべき関数は

$$
V=r^{2} + 4 r \cos{\left(\phi \right)} + 4,\qquad X^r=2 \sin{\left(\phi \right)},\quad X^\phi=\frac{r + 2 \cos{\left(\phi \right)}}{r}.
$$

まず力の変化を調べる。$\partial_rV$ と $\partial_\phi V$ はそれぞれ $2 \left(r + 2 \cos{\left(\phi \right)}\right)$ と $- 4 r \sin{\left(\phi \right)}$ である。これらに $X^r,X^\phi$ を掛けて足すと零になる。次に線素 $\dd r^2+r^2\dd\phi^2$ を調べる。方向場の微分は

$$
\frac{\partial(X^r,X^\phi)}{\partial(r,\phi)}=\left(\begin{matrix}0 & 2 \cos{\left(\phi \right)}\\- \frac{2 \cos{\left(\phi \right)}}{r^{2}} & - \frac{2 \sin{\left(\phi \right)}}{r}\end{matrix}\right).
$$

従って $D_{rr}=2\partial_r X^r$、$D_{r\phi}=\partial_\phi X^r+r^2\partial_rX^\phi$、$D_{\phi\phi}=2rX^r+2r^2\partial_\phi X^\phi$ はすべて消える。保存量を速度 $u=\dot r,w=\dot\phi$ で表示すると

$$
J_X=X^ru+r^2X^\phi w=r^{2} w + 2 r w \cos{\left(\phi \right)} + 2 u \sin{\left(\phi \right)}.
$$

例えば $(r,\phi,u,w)=(2,\frac{\pi}{2},1,2)$ でその数値は $10$ である。
:::

:::example[鉛直方向へ自由な溝]
井戸の断面を四次にしても、その軸に沿って断面が同じであれば保存量がある。ここでは平面上の $y$ 方向を軸とし、$x$ にだけ復元力が働く。質量を $1\,\mathrm{kg}$ として SI の数値を用いる。$\delta r=\epsilon X^r,\ \delta\phi=\epsilon X^\phi$ と書いたとき、座標上で比較すべき関数は

$$
V=r^{4} \cos^{4}{\left(\phi \right)},\qquad X^r=\sin{\left(\phi \right)},\quad X^\phi=\frac{\cos{\left(\phi \right)}}{r}.
$$

まず力の変化を調べる。$\partial_rV$ と $\partial_\phi V$ はそれぞれ $4 r^{3} \cos^{4}{\left(\phi \right)}$ と $- 4 r^{4} \sin{\left(\phi \right)} \cos^{3}{\left(\phi \right)}$ である。これらに $X^r,X^\phi$ を掛けて足すと零になる。次に線素 $\dd r^2+r^2\dd\phi^2$ を調べる。方向場の微分は

$$
\frac{\partial(X^r,X^\phi)}{\partial(r,\phi)}=\left(\begin{matrix}0 & \cos{\left(\phi \right)}\\- \frac{\cos{\left(\phi \right)}}{r^{2}} & - \frac{\sin{\left(\phi \right)}}{r}\end{matrix}\right).
$$

従って $D_{rr}=2\partial_r X^r$、$D_{r\phi}=\partial_\phi X^r+r^2\partial_rX^\phi$、$D_{\phi\phi}=2rX^r+2r^2\partial_\phi X^\phi$ はすべて消える。保存量を速度 $u=\dot r,w=\dot\phi$ で表示すると

$$
J_X=X^ru+r^2X^\phi w=r w \cos{\left(\phi \right)} + u \sin{\left(\phi \right)}.
$$

例えば $(r,\phi,u,w)=(1,\frac{\pi}{4},2,1)$ でその数値は $\frac{3 \sqrt{2}}{2}$ である。
:::

### 時間に依存する振動子にも残る二次の積分

ばねの硬さを変えると $E_L$ は一般に保存しない。それでも線形方程式では、伸びを測る物差しと時計の両方を変えると一定の振動数に帰着できる場合がある。$\ddot q+\omega(t)^2q=0$ とし、正の無次元関数 $\rho(t)$ と定数 $\omega_*>0$ を導入する。$Q=q/\rho$ と $\dd\sigma=\dd t/\rho^2$ を定義すれば

$$
\frac{\dd Q}{\dd\sigma}=\rho\dot q-\dot\rho q,\qquad
\frac{\dd^2 Q}{\dd\sigma^2}
=\rho^2(\rho\ddot q-\ddot\rho q)
=-\rho^3(\ddot\rho+\omega^2\rho)Q .
$$

従って $\ddot\rho+\omega^2\rho=\omega_*^2/\rho^3$ を満たす物差しを選べば、$\sigma$ で見た振動子のエネルギーが保存する。$\rho>0$ の初期値を与えれば右辺は滑らかなので少なくとも局所解がある。保存量は

$$
I=\frac m2\left[(\rho\dot q-\dot\rho q)^2+
\omega_*^2 q^2/\rho^2\right].
\label{eq-am-ermakov-energy}
$$

これを :term[Ermakov–Lewis 不変量]{en="Ermakov–Lewis invariant" yomi="えるまこふるいすふへんりょう"} という。ここでの規格化では単位は J である。Noether の形式との対応を見るには先の決定方程式に $\tau=\rho^2,\xi=\rho\dot\rho q$ を代入する。$B=m(\dot\rho^2+\rho\ddot\rho)q^2/2$ として $p\xi-E_L\tau-B$ を整理すれば $-I$ になる。$\omega$ がゆっくり変わるという条件はこの導出にはない。

:::example[物差しを先に決めてばねを設計する]
$\rho=1+t/t_0$、$t\ge0$ と決めれば $\ddot\rho=0$ なので $\omega=\omega_*/\rho^2$ とすればよい。$t_0=1\,\mathrm s$、$\omega_*=2\,\mathrm{s^{-1}}$、$m=1\,\mathrm{kg}$、$q(0)=0.10\,\mathrm m,\dot q(0)=0.10\,\mathrm{m/s}$ とする。変換後は $Q(0)=0.10\,\mathrm m$、$\dd Q/\dd\sigma(0)=0$。従って $\sigma=t_0(1-1/\rho)$、$q=0.10\rho\cos(\omega_*\sigma)\,\mathrm m$ であり、$I=0.020\,\mathrm J$ が全時刻で保存する。通常の初期エネルギーは $0.025\,\mathrm J$。どちらもエネルギーと同じ単位を持つが、ばねを通じて外部へ渡った仕事を数えるには通常の $E_L$ を使う。
:::

### 第一積分から変分対称性を逆に作る

保存則を得る向きの逆も局所的にはたどれるが、許す変換の範囲が重要である。正則な $L$ で速度の Hessian を $W$ とし、任意の初期値からの解に沿って $J(q,\dot q,t)$ が保存するとする。運動方程式を解いた加速度を $a_*(q,\dot q,t)$ とすれば、任意の試験経路の加速度 $a$ について
$$
D J=\sum_i J_{\dot q^i}(a^i-a_*^i),\qquad
D p_i-L_{q^i}=\sum_j W_{ij}(a^j-a_*^j).
$$
第一式では加速度を $a_*$ にしたときに全微分が零になることを使った。従って $\eta=W^{-1}\nabla_{\dot q}J$ と選べば $DJ=\sum_i(Dp_i-L_{q^i})\eta^i$ が試験経路上でも成立する。時刻を動かさず $\delta q^i=\epsilon\eta^i$ とし、速度の変化を $\epsilon D\eta^i$ とすれば
$$
\delta L/\epsilon=\sum_i(L_{q^i}\eta^i+p_iD\eta^i)
=D\!\left(\sum_i p_i\eta^i-J\right).
$$
得られた境界項から Noether 量を作れば初めの $J$ が戻る。この変換は一般に速度を使い、境界項も速度に依存し得る。任意の第一積分が配位空間の点変換に由来する、という結論を出すことはできない。

例えば二次元等方振動子で $J=m\dot x\dot y+m\omega^2xy$ とする。二つの方程式を代入すると $\dot J=0$。質量行列は $m1$ だから $\eta=(\dot y,\dot x)$ となり、境界関数は $B=m\dot x\dot y-m\omega^2xy$。実際に $DB=m(\ddot x\dot y+\dot x\ddot y)-m\omega^2(\dot x y+x\dot y)$ が試験経路の $\delta L/\epsilon$ に一致する。この例の積分は二軸の振動の位相関係を記憶する。時間並進だけから得る全エネルギーではその関係を知ることはできない。

## 空間の一様性・等方性と保存則

### 原点を移す

孤立した $N$ 質点の $L=\sum_{a=1}^Nm_a|\dot{\vb r}_a|^2/2-V$ で、$V$ が相対位置だけに依存すれば、全粒子を同じ $\eps\vb e$ だけ移しても $L$ は変わらない。
任意の定ベクトル $\vb e$ に対する保存量が $\vb e\cdot\sum_a m_a\dot{\vb r}_a$ なので、独立な三方向を選ぶと :term[全運動量]{en="total momentum" yomi="ぜんうんどうりょう"} $\vb P=\sum_a m_a\dot{\vb r}_a$ が保存する。
全質量を $M$、重心を $\vb R$ とすれば $M\dot{\vb R}=\vb P$ であり、重心は等速直線運動をする。

例えば質量 $2,3\,\mathrm{kg}$ の台車が同一直線上を $3,-1\,\mathrm{m/s}$ で走るとき、全運動量は $3\,\mathrm{kg\,m/s}$ である。
台車間でばねが伸縮しても重心速度 $0.6\,\mathrm{m/s}$ は保たれる。

### 軸の向きを変える

定単位ベクトル $\vb n$ のまわりの回転で $\xi_a=\vb n\times\vb r_a$ である。
全体の回転により相対距離も速さも変わらない。
従って距離に依存する内部ポテンシャルに対し、保存量は

$$
\sum_{a=1}^N\vb p_a\cdot(\vb n\times\vb r_a)
=\vb n\cdot\sum_{a=1}^N\vb r_a\times\vb p_a
$$

となる。
任意の $\vb n$ が使えるとき、:term[角運動量]{en="angular momentum" yomi="かくうんどうりょう"} の三成分が保存する。
鉛直な一様重力のように軸対称性だけが残るときには、その軸への射影だけが保証される。
固定された原点を $\vb a$ だけ移せば角運動量は $-\vb a\times\vb P$ だけ変わる。

:::example[どの対称性が残るか]
$V=k(x^2+y^2)/2+mgz$ は鉛直軸まわりで不変である。
したがって $m(x\dot y-y\dot x)$ が保存する。
鉛直並進 $z\mapsto z+\eps$ では $\Delta=-mg=D(-mgt)$ だから、厳密な不変性は失われても $p_z+mgt$ は保存する。
固定した鉛直平面での鏡映にも $V$ は不変だが、離散変換には恒等変換からの無限小パラメータがなく、上の定理を直接適用できない。
:::

### 時計の原点を移す

$t'=t+\eps$、$\xi^i=0$ とすれば、$L_t=0$ のとき $B=0$ で $J=-E_L$ となる。
粒子の配置に関する一様性と時刻に関する一様性は、同じ作用の微分の異なる成分を調べている。
空間についても時間についても、変換する対象には外部ポテンシャルと拘束の条件を含める必要がある。

### 重心の積分と角運動量の分解

自由粒子で見た等速変換を孤立多体系に施す。
$x$ 方向に各位置を $\epsilon t$ だけ移すと相互距離は不変で、$\Delta=\sum_a m_a\dot x_a=D(MX)$ である。
Noether 量は $tP_x-MX$。
三方向を合わせると

$$
M\vb R(t)-t\vb P=M\vb R(0)
$$

が成り立つ。
並進対称性が一定の重心速度を与え、等速変換の対称性が重心位置の積分定数を与える。

重心からの変位を $\rho_a=\vb r_a-\vb R$ とし、角運動量を展開すると

$$
\sum_a\vb r_a\times m_a\dot{\vb r}_a
=\vb R\times\vb P+\sum_a\rho_a\times m_a\dot\rho_a .
$$

$\sum_am_a\rho_a=\sum_am_a\dot\rho_a=0$ なので交差項は消える。
左辺が原点まわりの全角運動量、第一項が重心の公転、第二項が重心まわりの内部運動である。
外力の合力が零なら $\vb R\times\vb P$ の微分も零である。さらに全角運動量が保存するなら、差を取ると内部の角運動量も保存する。

### 運動方程式が不変な変換と作用の変化

時間を $t'=\lambda^\alpha t$、位置を $\vb r'=\lambda\vb r$ に変える。
斉次ポテンシャル $V(\lambda\vb r)=\lambda^s V(\vb r)$ のもとで Newton 方程式の両辺はそれぞれ $\lambda^{1-2\alpha}$ と $\lambda^{s-1}$ 倍になる。
従って

$$
\alpha=1-s/2
$$

なら元の解から別の解が作れる。
Kepler ポテンシャルの $s=-1$ では $t'=\lambda^{3/2}t$ となり、軌道の尺度と周期の関係が直ちに分かる。

同じ変換で $T,V$ はともに $\lambda^s$ 倍だが、時間要素は $\lambda^\alpha$ 倍なので、作用全体は $\lambda^{1+s/2}$ 倍になる。
停留条件は保たれても、一般には作用の差は境界項にならない。
例外の $s=-2$ では倍率が $1$ となり、上で求めた膨張の保存量 $px-2tE_L$ が現れる。
この比較により、解を別の解に写すだけの変換と、保存量を生む変換を区別できる。

同じ斉次性から時間平均についても恒等式が得られる。
$G=\sum_a\vb r_a\cdot\vb p_a$ と置けば

$$
\dot G=2T-\sum_a\vb r_a\cdot\nabla_aV=2T-sV .
$$

周期軌道では一周期で $G$ が戻る。
より一般に $G$ が有界なら、時間 $t$ で割った端点差は長時間で零になる。
従って平均が存在する場合 $2\langle T\rangle=s\langle V\rangle$ である。
例えば重力束縛系なら $\langle T\rangle=-\langle V\rangle/2$。
この :term[ビリアル定理]{en="virial theorem" yomi="びりあるていり"} は時間平均の関係であり、各瞬間に $2T=sV$ を課すものではない。

### 保存を破る外部装置を同定する

重力中で長さ $l(t)$ をモーターで指定した平面振り子を考える。
支点固定、糸の方向を下向きから $\theta$ とすると

$$
L=\frac m2(\dot l^2+l^2\dot\theta^2)+mgl\cos\theta .
$$

$l$ は指定された時刻の関数であり、一般化座標は $\theta$ だけである。
共役運動量とエネルギー関数は

$$
p_\theta=ml^2\dot\theta,\qquad
E_L=\frac m2l^2\dot\theta^2-\frac m2\dot l^2-mgl\cos\theta .
$$

Lagrange 方程式から

$$
\dot p_\theta=-mgl\sin\theta,\qquad
\ddot\theta+2\frac{\dot l}{l}\dot\theta+\frac gl\sin\theta=0 .
$$

伸長時に現れる速度比例項は、摩擦を仮定した結果ではない。
慣性モーメント $ml^2$ 自体が増加するためである。

長さの時間変化があると $E_L$ は一般に保存せず、実験室の力学的エネルギーは $E_{\rm lab}=E_L+m\dot l^2$ である。
上式を用いて微分を整理すると

$$
\dot E_{\rm lab}=m(\ddot l-l\dot\theta^2-g\cos\theta)\dot l .
$$

括弧に $m$ を掛けたものは外向きを正に取った拘束力である。
例えば $\theta=\dot\theta=\ddot l=0$ で糸を一定速度で延ばせば $\dot E_{\rm lab}=-mg\dot l$。
糸が下向きに伸びる間、上向きの張力が負の仕事をする。
系と装置の仕事の境界を決めれば、どの量を保存と呼べるかを実験の条件に結び付けられる。

### 力学的相似を解の初期条件まで追う

同じ装置を大きく作ったとき、元の軌道を定規と時計の目盛の変更だけで再現できる場合がある。先の斉次性を使うには力の係数を変えず、すべての初期条件にも同じ変換を施す必要がある。長さを $a$ 倍にして $t'=a^{1-s/2}t$ とすれば $\dot q'=a^{s/2}\dot q$、$\ddot q'=a^{s-1}\ddot q$ である。自然系では運動量は速度と同じ倍率、角運動量と一周期の作用積分は長さと速度の積の倍率となる。時間平均は新しい周期で規格化するから、エネルギーと同じ倍率になる。以下では軌道上の量を比較し、周期が存在しない運動については対応する二点間の飛行時間を比較する。

:::example[どの測定量を何倍にするか]
質量と力の係数を固定し、元の軌道の二点間の時間を $2\,\mathrm s$、第一の点での速さを $3\,\mathrm{m/s}$ とする。次の表で、$a$ は初期位置を含む全長さの倍率である。対応する二点で $t'/t=a^{1-s/2},\,v'/v=a^{s/2}$ とし、作用の倍率は $a^{1+s/2}$ を用いた。

| ポテンシャルと装置 | $s$ | $a$ | 飛行時間（s） | 初速（m/s） | 作用の倍率 |
| --- | --- | --- | --- | --- | --- |
| 線形ばね | $2$ | $2$ | $2$ | $6$ | $4$ |
| 円錐面上の重力 | $1$ | $4$ | $4$ | $6$ | $8$ |
| 二体間の重力 | $-1$ | $4$ | $16$ | $3/2$ | $2$ |
| 逆二乗ポテンシャル | $-2$ | $3$ | $18$ | $1$ | $1$ |

第一行では長さを倍にすると復元加速度も倍になり、初速度を倍に選んだ質点は同じ時間で対応点に着く。第二行では円錐の頂点を原点に取る必要がある。円錐の開き角が同じなら斜面に沿う重力加速度が同じなので、飛行時間は距離の平方根に比例する。第三行では $r\mapsto4r$ で引力が $1/16$ になり、距離を速さで割ると時間の倍率は $4/(1/2)=8$。第四行では $p$ と $\dd q$ の倍率が打ち消し合うため一周に限らず任意の対応区間で短縮作用が等しい。散乱を含む比較では初期位置だけでなく、入射速度と衝突径数にも相似を適用する。
:::

相似から平均エネルギーにも独立な検算が得られる。代表長さ $r_*$ の軌道で $V$ の尺度を $V_*\propto r_*^s$ とすると、運動量の尺度は $\sqrt{m|V_*|}$。従って作用積分の単位をそろえた値は $r_*\sqrt{m|V_*|}$ に比例する。Kepler 楕円を相似に拡大すれば、束縛エネルギーの絶対値は小さくなるのに作用積分は増す。この結果は、定数倍を除いて $|E|\propto I^{-2}$ という関係に相当し、[作用角変数による計算](/physics/analytical-mechanics/hamilton-jacobi/)でも確かめられる。

### 離散対称性が与える初期値の制限

鏡映と時間反転には恒等変換から連続に近づく群パラメータがない。従って変換をパラメータで微分して Noether 量を作る手順は適用できない。一方、解の一意性と組み合わせれば、運動をある面に閉じ込めたり、折返し点の前後を結び付けたりできる。

例えば平面上の $V(x,y)=V(x,-y)$ は $y=0$ で $V_y=0$ を満たす。$y(0)=\dot y(0)=0$ の初期条件を与えると、鏡映した解と元の解は初期状態を共有する。力が局所 Lipschitz なら一意性により両者は一致し、従って $y(t)=0$ である。しかし $y(0)$ をずらした粒子の $p_y$ はふつう一定にならない。例として $V=k(x^2+4y^2)/2$ では $\dot p_y=-4ky$ である。二つの座標の方程式が独立なことも、鏡映対称性そのものからは従わない。$x^2y^2$ の相互作用を加えても鏡映は保たれる。

時間反転については、$L(q,-\dot q)=L(q,\dot q)$ の自律系なら $q(-t)$ も解である。全速度が零となる時刻 $t_*$ があると、$q(t_*+\delta)$ と $q(t_*-\delta)$ は同じ初期位置と零の速度から出発するので一致する。一つの座標の速さだけが零になる時刻にはこの議論を使えない。磁場中では $m\ddot{\vb r}=e\dot{\vb r}\times\vb B$ の左辺が時間反転で変わらず、右辺の速度だけが反転する。磁場を固定した同一装置で単に同じ空間曲線を逆走させることは一般にできない。磁場まで反転した装置なら変換後の運動が得られる。

## 演習問題

::::exercise{level=basic}
半径 $0.5\,\mathrm m$ の円上を質量 $2\,\mathrm{kg}$ の小球が角速度 $3\,\mathrm{s^{-1}}$ で進む。
共役運動量を求めよ。

:::solution
$L=mr^2\dot\phi^2/2$ を $\dot\phi$ で微分すると $p_\phi=mr^2\dot\phi=1.5\,\mathrm{kg\,m^2/s}$ である。
:::
::::

::::exercise{level=basic}
$L=m\dot x^2/2+ax\dot x-kx^2/2$ で $p$ と $E_L$ を求めよ。

:::solution
$p=m\dot x+ax$ となる。
$p\dot x-L=m\dot x^2/2+kx^2/2$ であり、$ax\dot x=D(ax^2/2)$ はエネルギーに寄与しない。
:::
::::

::::exercise{level=basic}
$V=mgz$ で保存する通常の運動量成分と角運動量成分を、原点を固定して選べ。

:::solution
$x,y$ 並進と $z$ 軸回転は $V$ を変えないので $p_x,p_y,L_z$ である。
$\dot{\vb p}=-mg\vb e_z$、$\dot{\vb L}=\vb r\times(-mg\vb e_z)$ を計算しても同じ三成分の微分が零になる。
:::
::::

::::exercise{level=basic}
$\mathcal F=b\dot x^2/2$ の抵抗下で $L=m\dot x^2/2-kx^2/2$ とする。
単位時間の損失を求めよ。

:::solution
$Q=-\partial_{\dot x}\mathcal F=-b\dot x$ をエネルギー収支に入れると $\dot E_L=-b\dot x^2$ である。
$b=0.2\,\mathrm{kg/s}$、$\dot x=3\,\mathrm{m/s}$ なら損失は $1.8\,\mathrm W$ である。
:::
::::

::::exercise{level=advanced}
$L=m(\dot x^2+\dot y^2)/2-k(x-y)^2/2$ に対し、循環座標を選んで二つの運動を分けよ。

:::solution
$X=(x+y)/2,\ s=x-y$ を用いる。
$L=m\dot X^2+m\dot s^2/4-ks^2/2$ より $2m\dot X$ は一定で、$(m/2)\ddot s=-ks$ である。
重心の等速運動と角振動数 $\sqrt{2k/m}$ の相対振動が得られる。
:::
::::

::::exercise{level=advanced}
一様重力下で $z'=z+\eps t$ とする。
境界関数と保存量を求めよ。

:::solution
$\Delta=m\dot z-mgt=D(mz-mgt^2/2)$ である。
$\xi=t$ より $J=m\dot zt-mz+mgt^2/2$ で、微分は $mt(\ddot z+g)=0$ となる。
:::
::::

::::exercise{level=advanced}
球面振り子について $p_\phi=c$ として $\theta$ の有効ポテンシャルを求め、定傾角の条件を示せ。

:::solution
$R=ml^2\dot\theta^2/2-U$、$U=c^2/(2ml^2\sin^2\theta)-mgl\cos\theta$ となる。
$U'=0$ は $mgl\sin\theta=c^2\cos\theta/(ml^2\sin^3\theta)$ である。
軸以外で $c=ml^2\sin^2\theta\,\dot\phi$ を戻すと $\dot\phi^2\cos\theta=g/l$ を得る。
:::
::::

::::exercise{level=advanced}
一自由度で $L\mapsto L+DF$ により、循環座標が式に陽に現れる例を作り、元の保存量を復元せよ。

:::solution
自由粒子で $F=ax^2t/2$ とする。
$L'=m\dot x^2/2+axt\dot x+ax^2/2$ は $x$ に依存する。
並進の変化 $\partial_x L'=at\dot x+ax=D(axt)$ より $B=axt$ であり、$p'-B=m\dot x$ が保存する。
:::
::::

::::exercise{level=advanced}
二質点の高さを固定した水平面上で、二つの位置の差だけに依存するばねを考える。
回転対称性の保存量を重心運動と相対運動に分け、それぞれが保存する条件を明らかにせよ。

:::solution
相対位置を $\vb r=\vb r_1-\vb r_2$、重心を $\vb R$ と書く。
変換は $\delta\vb r_a=\epsilon\vb e_z\times\vb r_a$ なので、保存量の候補は

$$
J=\vb e_z\cdot\left(\vb R\times M\dot{\vb R}+\vb r\times\mu\dot{\vb r}\right).
$$

ここで $\mu=m_1m_2/M$、相対項の係数は二つの質点の速度を展開して得られる。
外力がなければ第一項は自動的に保存する。
第二項の微分は $-\vb e_z\cdot(\vb r\times\nabla_{\vb r}V)$。
従って二体間距離だけに依存するばねなら零になる。
差ベクトルだけに依存しても、例えば $V=k(x_1-x_2)^2/2$ なら相対位置と力が一般に平行でなく、運動量は保存しても回転の対称性は保証されない。
:::
::::

::::exercise{level=advanced}
自由粒子で時刻の変換係数を $\tau=t^2$、座標の変換係数を $\xi=tq$ とする。
境界項を見つけ、保存量を初期条件で表せ。

:::solution
$L=m\dot q^2/2$ で $D\tau=2t,\ D\xi=q+t\dot q$ だから、

$$
\Delta=m\dot q(q-t\dot q)+2tL=mq\dot q=D(mq^2/2).
$$

従って $B=mq^2/2$。
Noether 量は

$$
J=ptq-E_Lt^2-mq^2/2=-\frac{m}{2}(q-t\dot q)^2.
$$

自由粒子の解 $q=q_0+v_0t$ を入れると $J=-mq_0^2/2$。
これは等速変換の保存量から作れるので、連続な対称性をもう一つ見つけても独立な初期定数が増えるとは限らない。
:::
::::

::::exercise{level=advanced}
一様な重力中の鉛直座標を $z$、上向きを正とする。二つの変換 $z\mapsto z+a$ と $z\mapsto z+bt$ の境界関数を求め、得られた保存量から運動を復元せよ。$m=2\,\mathrm{kg}$、$g=9.8\,\mathrm{m/s^2}$、$z(0)=1\,\mathrm m$、$\dot z(0)=3\,\mathrm{m/s}$ とする。

:::solution
並進に対する単位パラメータあたりの変化は $-mg=D(-mgt)$ なので $J_a=m\dot z+mgt$。等速変換では $m\dot z-mgt=D(mz-mgt^2/2)$ なので $J_b=m\dot zt-mz+mgt^2/2$ である。初期値から $J_a=6\,\mathrm{kg\,m/s},J_b=-2\,\mathrm{kg\,m}$。第一式で速度を消すと $z=(J_at-J_b)/m-gt^2/2$、すなわち SI の数値で $z=1+3t-4.9t^2$ となる。保存量の値だけで軌道が決まるのは、二つの式が初期位置と初期速度を独立に指定するためである。
:::
::::

## 確認テスト

::::quiz{answer="2" tol="1e-6"}
$m=1\,\mathrm{kg},r=2\,\mathrm m,\dot\phi=0.5\,\mathrm{s^{-1}}$ で $p_\phi$ は SI 単位でいくらか。

:::explanation
角度に共役な運動量は $mr^2\dot\phi=1\cdot4\cdot0.5=2$ である。
:::
::::

::::quiz
時間によらない $L$ に全時間微分 $DF(q,t)$ を加えると $E_L$ の変化はどれか。

- [ ] $+DF$
- [x] $-F_t$
- [ ] $+F_t$

:::explanation
運動量の増分と $\dot q$ の積を全時間微分から差し引くと、$-F_t$ だけが残る。
:::
::::

::::quiz{answer="3" tol="1e-6"}
$k(t)$ が毎秒 $6\,\mathrm{N/m}$ 増えるばねが $1\,\mathrm m$ 伸びているとき、供給される仕事率を W で答えよ。

:::explanation
$\dot E_L=-L_t=\dot k x^2/2=3\,\mathrm W$ である。
:::
::::

::::quiz
回転対称な運動で、無限小回転の係数はどれか。

- [ ] $\vb n\cdot\vb r$
- [ ] $\dot{\vb r}\times\vb n$
- [x] $\vb n\times\vb r$

:::explanation
回転角を微分すると半径に垂直な接線方向の変位が得られる。
右手系の正の回転ではその係数が $\vb n\times\vb r$ である。
:::
::::

::::quiz
空間並進の対称性で境界関数 $B$ が零でないとき、保存量として選べるものはどれか。
変換は $x'=x+\epsilon$ とする。

- [ ] 常に $m\dot x$ だけである
- [x] $p_x-B$ である
- [ ] $p_x+\dot B$ である

:::explanation
$\xi=1,\tau=0$ を Noether の式に代入すれば $J=p_x-B$ が得られる。
運動方程式を使うと $\dot p_x=\partial_x L=DB$ であり、直接微分しても同じ結論になる。
:::
::::

::::quiz{answer="0.75" tol="1e-6"}
質量 $1\,\mathrm{kg}$ の二質点が $x$ 軸上で $x_1=0,\ x_2=1.5\,\mathrm m$ にある。
孤立系の重心位置の初期値を m で答えよ。

:::explanation
重心位置は $(m_1x_1+m_2x_2)/(m_1+m_2)=1.5/2=0.75\,\mathrm m$ である。
運動開始後も $X(t)-P_xt/M$ にこの初期値が保存される。
:::
::::

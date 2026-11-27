---
title: 中心力問題と Kepler 問題
description: 二体運動を一つの動径に還元し、束縛軌道の閉鎖性と散乱角を同じ運動方程式から求める。
notation:
  - symbol: 'M,\mu'
    meaning: 全質量と換算質量
  - symbol: '\ell'
    meaning: 相対運動の角運動量の大きさ
  - symbol: 'u'
    meaning: 距離の逆数
  - symbol: 'b,\chi'
    meaning: 衝突径数と散乱角
---

## この章の目標

- 重心運動を分離し、転回点と円軌道の安定性を判定する。
- Binet の式を導き、閉じる束縛軌道の条件を調べる。
- Kepler 軌道の形、周期、近日点の向きを保存量から求める。
- 入射粒子の面積分布を立体角に変換して散乱断面積を求める。

## 二体問題の換算

### どちらを中心に置くか

互いに引き合う二つの質点が動くとき、一方を原点に固定した座標系も加速する。
両者の質量が同程度ならその効果は無視できない。
そこでまず慣性系で位置 $\vb r_1,\vb r_2$ を測り、固定質量 $m_1,m_2>0$、外力なし、相互作用 $V(|\vb r_1-\vb r_2|)$ を仮定する。

:::definition[重心と換算質量]{#def-am-reduced-mass}
$M=m_1+m_2$、$\vb R=(m_1\vb r_1+m_2\vb r_2)/M$、$\vb r=\vb r_1-\vb r_2$ とおく。
相対運動の慣性を表す $\mu=m_1m_2/M$ を :term[換算質量]{en="reduced mass" yomi="かんさんしつりょう"} という。
:::

逆変換は $\vb r_1=\vb R+(m_2/M)\vb r$、$\vb r_2=\vb R-(m_1/M)\vb r$ である。
速度の二乗を展開すると交差項の係数は $m_1m_2/M-m_2m_1/M=0$ となり、相対速度の二乗の係数は $\mu/2$ になる。
従って

$$
L=\tfrac12 M|\dot{\vb R}|^2+\tfrac12\mu|\dot{\vb r}|^2-V(r),\qquad r=|\vb r|.
$$

重心の三座標は :ref[def-am-cyclic] であり、$M\ddot{\vb R}=0$。
一方、相対位置の方程式は $\mu\ddot{\vb r}=F(r)\vb r/r$、$F=-V'$ である。
$\vb R$ を消去しても、残った問題には両質量が入る。

:::example[質量比の効果]
$m_1=2\,\mathrm{kg},m_2=6\,\mathrm{kg}$ なら $\mu=1.5\,\mathrm{kg}$ である。
重心から軽い方までの距離は $3r/4$、重い方までは $r/4$ である。
逆に $m_2/m_1=100$ なら $\mu/m_1=100/101$ となる。
中心を固定して $\mu=m_1$ と近似する誤差の尺度は $m_1/m_2$ である。
:::

### 面積速度

$\vb\ell=\vb r\times\mu\dot{\vb r}$ の微分は $\vb r\times F\vb r/r=0$ である。
$\ell=|\vb\ell|>0$ なら $\vb r$ は一定の $\vb\ell$ に垂直な平面内にあり、$\vb\ell$ 側から見て $\phi$ が増える向きを選べる。
単位時間に掃く扇形の面積は

$$
\dot A=r^2\dot\phi/2=\ell/(2\mu).
$$

例えば $\ell=6\,\mathrm{kg\,m^2/s}$、$\mu=2\,\mathrm{kg}$ なら $\dot A=1.5\,\mathrm{m^2/s}$ である。
$\ell=0$ の運動は直線上の動径運動として別に解く。

### 初期値を二つの系の間で移す

実験室で片方が静止する衝突を考える。
入射粒子の速度が $\vb v$、標的の速度が零なら、

$$
\dot{\vb R}=\frac{m_1}{M}\vb v,\qquad
\dot{\vb r}=\vb v,\qquad E_{\rm rel}=\tfrac12\mu|\vb v|^2
$$

である。
入射運動エネルギー $E_{\rm lab}=m_1|\vb v|^2/2$ のうち、相対運動に使える割合は $m_2/M$ である。
例えば等質量なら半分だけで、残りは重心の並進として反応の前後で同じ値をもつ。

:::example[等質量の重力二体]
質量 $m$ の二粒子が間隔 $d$ を保って回る。
各粒子の回転半径は $d/2$、角速度は共通の $\Omega$ である。
各粒子の式は $m(d/2)\Omega^2=Gm^2/d^2$ であるから

$$
\Omega^2=2Gm/d^3 .
$$

相対座標を使うと $\mu=m/2$、回転半径は $d$ であり、$(m/2)d\Omega^2=Gm^2/d^2$ と同じ答えになる。
一粒子に換算した後でも、実際の各質点が重心のまわりに描く円の半径は $d/2$ のままである。
:::

### 実験室の観測から相対軌道へ移る

片方だけが見える二体系では、観測した位置をそのまま相対位置に代入すると力定数を誤って推定する。重心を原点とし、質点 1 の位置を $\vb x$ とすれば $\vb x=\alpha\vb r$、$\alpha=m_2/(m_1+m_2)$ である。逆二乗引力を $V=-k/r$ とすると $\ddot{\vb r}=-(k/\mu)\vb r/r^3$ なので、測定した軌道では

$$
\ddot{\vb x}=-\frac{\alpha^3 k}{\mu}\frac{\vb x}{|\vb x|^3}.
$$

重力なら $k=Gm_1m_2$ であり、係数は $Gm_2^3/(m_1+m_2)^2$。従って片方の軌道の長半径と周期だけで全質量を決めることはできない。両方が観測できれば、重心からの二つの長半径の和が相対運動の長半径であり、その比が逆質量比である。さらに視線速度だけを測った場合、速度の振幅には軌道面の傾き $i$ による $\sin i$ が掛かる。円軌道で速度振幅を $K_1$、周期を $T$ とすると $a_1\sin i=K_1T/(2\pi)$。これを上の係数に代入して

$$
\frac{T K_1^3}{2\pi G}=\frac{m_2^3\sin^3 i}{(m_1+m_2)^2}
$$

を得る。$0<\sin i\le1$ なので、左辺の測定値は $m_2^3/(m_1+m_2)^2$ の下限になる。$m_1$ が別に分かれば、右辺が $m_2>0$ で単調増加することから伴星質量の下限を決められる。

:::example[軽い方だけを測った場合]
位置は m、速度は $\mathrm{m/s}$、質量は kg で表す。重心系で $m_1=1,m_2=3$、質点 1 の観測値を $\vb x=\left(\begin{matrix}3 & 0 & 0\end{matrix}\right)^\top,\ \dot{\vb x}=\left(\begin{matrix}0 & 3 & 0\end{matrix}\right)^\top$ とする。質点 2 の位置は重心条件 $m_1\vb x+m_2\vb x_2=0$ から決まり、速度についてもその時間微分を使える。従って二体間の位置と速度は

$$
\vb r=\left(\begin{matrix}4 & 0 & 0\end{matrix}\right)^\top,\qquad \dot{\vb r}=\left(\begin{matrix}0 & 4 & 0\end{matrix}\right)^\top.
$$

換算質量は $\mu=\frac{3}{4}\,\mathrm{kg}$ である。外積をとると $\boldsymbol\ell=\mu\vb r\times\dot{\vb r}=\left(\begin{matrix}0 & 0 & 12\end{matrix}\right)^\top\,\mathrm{kg\,m^2/s}$、相対運動の面積速度は $|\boldsymbol\ell|/(2\mu)=8\,\mathrm{m^2/s}$。一方、測定した軌道の面積速度は $|\vb x\times\dot{\vb x}|/2=\frac{9}{2}\,\mathrm{m^2/s}$ であり、両者の比は長さの比の二乗になる。半径をそのまま相対運動の半径にすると、小さな軌道を見ている分だけ力定数を過小評価する。
:::

:::example[質量が等しい連星の斜めの速度]
位置は m、速度は $\mathrm{m/s}$、質量は kg で表す。重心系で $m_1=2,m_2=2$、質点 1 の観測値を $\vb x=\left(\begin{matrix}0 & 2 & 0\end{matrix}\right)^\top,\ \dot{\vb x}=\left(\begin{matrix}1 & 1 & 0\end{matrix}\right)^\top$ とする。従って二体間の位置と速度は

$$
\vb r=\left(\begin{matrix}0 & 4 & 0\end{matrix}\right)^\top,\qquad \dot{\vb r}=\left(\begin{matrix}2 & 2 & 0\end{matrix}\right)^\top.
$$

換算質量は $\mu=1\,\mathrm{kg}$ である。外積をとると $\boldsymbol\ell=\mu\vb r\times\dot{\vb r}=\left(\begin{matrix}0 & 0 & -8\end{matrix}\right)^\top\,\mathrm{kg\,m^2/s}$、相対運動の面積速度は $|\boldsymbol\ell|/(2\mu)=4\,\mathrm{m^2/s}$。一方、測定した軌道の面積速度は $|\vb x\times\dot{\vb x}|/2=1\,\mathrm{m^2/s}$ であり、両者の比は長さの比の二乗になる。動径速度もあるため、観測した速さをすべて接線速度に算入してはいけない。
:::

:::example[重心からの位置と速度が直交しない場合]
位置は m、速度は $\mathrm{m/s}$、質量は kg で表す。重心系で $m_1=4,m_2=2$、質点 1 の観測値を $\vb x=\left(\begin{matrix}1 & 2 & 2\end{matrix}\right)^\top,\ \dot{\vb x}=\left(\begin{matrix}2 & 0 & 1\end{matrix}\right)^\top$ とする。従って二体間の位置と速度は

$$
\vb r=\left(\begin{matrix}3 & 6 & 6\end{matrix}\right)^\top,\qquad \dot{\vb r}=\left(\begin{matrix}6 & 0 & 3\end{matrix}\right)^\top.
$$

換算質量は $\mu=\frac{4}{3}\,\mathrm{kg}$ である。外積をとると $\boldsymbol\ell=\mu\vb r\times\dot{\vb r}=\left(\begin{matrix}24 & 36 & -48\end{matrix}\right)^\top\,\mathrm{kg\,m^2/s}$、相対運動の面積速度は $|\boldsymbol\ell|/(2\mu)=\frac{9 \sqrt{29}}{2}\,\mathrm{m^2/s}$。一方、測定した軌道の面積速度は $|\vb x\times\dot{\vb x}|/2=\frac{\sqrt{29}}{2}\,\mathrm{m^2/s}$ であり、両者の比は長さの比の二乗になる。観測点での面積速度に寄与するのは速度の横成分である。軌道面を先に水平と仮定すると、傾いた面の射影面積しか測れない。
:::

## 有効ポテンシャルと軌道の分類

### 許容域、円運動、微小なずれ

:eqref[eq-am-polar-central] の Lagrangian で質量を $\mu$ とし、エネルギー $E=\mu(\dot r^2+r^2\dot\phi^2)/2+V(r)$ を作る。この式に $\dot\phi=\ell/(\mu r^2)$ を代入すると

$$
E=\tfrac12\mu\dot r^2+U_\ell(r),\qquad U_\ell=V(r)+\frac{\ell^2}{2\mu r^2}.
\label{eq-am-radial-energy}
$$

$U_\ell$ を :term[有効ポテンシャル]{en="effective potential" yomi="ゆうこうぽてんしゃる"} という。
運動可能なのは $E\ge U_\ell$ の連結区間である。
単純根 $E=U_\ell$ では $\dot r=0$ となり、$\mu\ddot r=-U_\ell'$ に従って戻る。
許容区間が両側で閉じていれば動径は束縛されるが、軌道が平面内で閉じるかどうかは別に調べる。

半径 $r_0$ の円では $U_\ell'(r_0)=0$、すなわち $\ell^2=\mu r_0^3V'(r_0)$ である。
$r=r_0+x$ として一次までとると $\mu\ddot x=-U_\ell''(r_0)x$ である。
正の曲率なら安定な動径振動の角振動数は $\kappa=\sqrt{U_\ell''/\mu}$、負ならずれが指数的に成長する。
曲率零では高次を調べる。
円運動自体の角速度 $\Omega=\ell/(\mu r_0^2)$ と比較すると

$$
\frac{\kappa^2}{\Omega^2}=3+\frac{r_0 V''(r_0)}{V'(r_0)}.
$$

:::example[異なる二つの周期]
$V=kr^2/2$ では $\kappa/\Omega=2$。
楕円の中心が力の中心であるため、一周の間に中心までの距離は二度極小になる。
$V=-k/r$ では比は $1$ で、一周に一度だけ近点を通る（円を除く）。
:::

### Binet の軌道方程式

時間を消し、$u=1/r$、プライムをこの段落では $\phi$ 微分とする。
$\dot\phi=\ell u^2/\mu$ だから

$$
\dot r=-\ell u'/\mu,\qquad
\ddot r=-\ell^2u^2u''/\mu^2,\qquad
r\dot\phi^2=\ell^2u^3/\mu^2 .
$$

動径方向の Newton 方程式 $\mu(\ddot r-r\dot\phi^2)=F$ に代入して

$$
u''+u=-\frac{\mu F(1/u)}{\ell^2u^2}
\label{eq-am-binet}
$$

を得る。
これが :term[Binet の公式]{en="Binet equation" yomi="びねのこうしき"} である。
例えば無力なら $u=C\cos(\phi-\phi_0)$、すなわち $r\cos(\phi-\phi_0)=1/C$ は直線である。

### 等方ばねで転回点と楕円を照合する

$V=kr^2/2$、$k>0$ として $\omega=\sqrt{k/\mu}$ とおく。直交座標で両成分が同じ角振動数をもつので、$\vb r(t)=\vb a\cos\omega t+\vb b\sin\omega t$ となる。時刻の原点を移すと二つの振幅ベクトルは同じ二次元の回転で混ざる。その Gram 行列を直交対角化する角を選べば両ベクトルを直交させられる。$\ell>0$ のとき独立な二軸の長さを $A\ge B>0$ とすると、力の中心を楕円の中心に置ける。$\ell=\mu\omega AB$、$E=k(A^2+B^2)/2$ である。動径の式も $r^2=A^2\cos^2\omega t+B^2\sin^2\omega t$ となり、非円の場合の動径周期は公転周期の半分である。

この結果は有効ポテンシャルでも検算できる。転回点を $s=r^2$ で表すと $ks^2-2Es+\ell^2/\mu=0$ なので、二根の和は $2E/k$、積は $\ell^2/(\mu k)$ である。上で得た $A^2,B^2$ はこの和と積を満たす。例えば $\mu=1\,\mathrm{kg}$、$k=4\,\mathrm{N/m}$、$A=2\,\mathrm m$、$B=1\,\mathrm m$ なら $\ell=4\,\mathrm{kg\,m^2/s}$、$E=10\,\mathrm J$。動径周期は $\pi/2\,\mathrm s$ で、公転周期は $\pi\,\mathrm s$ である。

### なぜすべての束縛軌道が閉じることは特別なのか

:::theorem[Bertrand の定理]{#thm-am-bertrand}
$r>0$ で滑らかな引力ポテンシャルを考え、各半径に非退化な安定円軌道があるとする。
半径が正のコンパクト区間にとどまる軌道がすべて閉じるための必要十分条件は、加法定数を除いて $V=-k/r$ または $V=kr^2/2$（$k>0$）である。
:::

:::proof
まず、近点から次の近点までに進む角を $\Phi$ とする。
閉軌道なら $\Phi/(2\pi)$ は有理数である。
円近傍では転回点を単純根として動かせば、動径振動の周期は初期振幅と円半径に連続に依存する。
実際、角度の積分で二つの転回点の間を $r=r_-+(r_+-r_-)\sin^2 s$ と置くと、端点の平方根特異性が相殺される。
連結区間から有理数への連続関数は一定だから、円軌道極限でも $\kappa/\Omega=\beta>0$ は一定である。
よって $rV''/V'=\beta^2-3$ を積分すると $V'=Cr^{\beta^2-3}$。
ここには $V$ が対数になる場合も含まれる。

円運動の一次近似だけでは、まだ無数の指数が残る。
有限振幅でも周期が変わらない条件を求める。
$y=u/u_0$ とすると :eqref[eq-am-binet] は $y''+y=y^{1-\beta^2}$ となる。
$x=y-1$ として展開すれば

$$
x''+\beta^2x=A x^2+B x^3+O(x^4),\quad
A=\frac{\beta^2(\beta^2-1)}2,\quad B=-\frac{\beta^2(\beta^2-1)(\beta^2+1)}6.
$$

角振動数を $\nu=\beta(1+d a^2+\cdots)$、位相を $s=\nu\phi$ とし、$x=a\cos s+a^2h(s)+O(a^3)$ とおく。
二次の方程式 $\beta^2(h''+h)=A(1+\cos 2s)/2$ の周期解は、基本調波を振幅に吸収すれば $h=A(1-\cos 2s/3)/(2\beta^2)$。
三次の右辺の基本調波は $2A\cos s\,h+B\cos^3s+2d\beta^2\cos s$ から生じる。
その係数が非零なら解に $s\sin s$ が現れ、周期解にならない。
従って

$$
d=-\frac{5A^2}{12\beta^4}-\frac{3B}{8\beta^2}
=\frac{(\beta^2-1)(4-\beta^2)}{24}.
$$

一定周期には $d=0$ が必要で、$\beta=1,2$ だけが残る。
この展開で使った周期と振幅の滑らかさは、変数 $x$ の一次元エネルギーから確かめられる。復元力のポテンシャルを $U(0)=0$ として $U(x)=\beta^2x^2/2+O(x^3)$ と書き、$y=x\sqrt{2U(x)/(\beta^2x^2)}$ と置く。$y=x+O(x^2)$ なので原点近傍で逆関数 $x=f(y)$ がある。動径エネルギーを $\beta^2a_0^2/2$ と書けば、転回点は $y=\pm a_0$ であり、$\phi$ で測った周期は
$$
\frac{2}{\beta}\int_{-\pi/2}^{\pi/2} f'(a_0\sin s)\,\dd s
$$
となる。積分範囲が固定されているので周期は滑らかで、$a_0$ の一次項は積分で消える。正の転回点から始めた解も初期振幅と周期に滑らかに依存する。基本調波の振幅は $a=a_0+O(a_0^2)$ だから、二次の周波数係数を求める上の計算が正当化される。

$\beta=1$ なら $V=-k/r$ で Binet 方程式の右辺は定数、従って正の $u$ が全角度で続く束縛解は $2\pi$ 周期である。
$\beta=2$ なら $V=kr^2/2$ で、直交二軸の運動はいずれも $\ddot x=-(k/\mu)x$ に従うから同じ周期で閉じる。
これで十分性も示した。
より一般の仮定における中心力の分類は :cite[goldstein-cm] を参照されたい。
:::

### 定理の仮定を弱められる範囲

上の必要性の議論は、ある半径 $r_*$ の安定円と十分近い振幅の閉軌道があれば、その円の近くの半径区間だけで実行できる。
実際、$V'(r_*)>0$ と $3V'(r_*)+r_*V''(r_*)>0$ は厳密な不等式なので近くでも続く。
各半径に $\ell^2=\mu r^3V'$ を選べばそこに安定円があり、その近傍の二つの単純転回点の間に束縛軌道がある。
したがって「すべての束縛軌道が閉じる」の仮定から、まずこの区間上で二つのポテンシャルのどちらかが得られる。

$V$ が連結した $r>0$ の領域で実解析的なら、同じ解析関数が開区間で一致するため領域全体で一致する。
解析性を用いない場合でも、滑らかな $V$ について非退化安定円の半径の集合を考えられる。
上の区間を含むその連結成分上では $V'=k/r^2$ または $V'=kr$ であり、$k>0$ は同じ定数である。
この成分が有限な正の半径で終わると仮定すると、端点での $V',V''$ は同じ式の極限をもち、二つの厳密な不等式は端点でもなお成り立つ。
すると成分を延長できて矛盾する。
従って $V$ が $r>0$ 全域で滑らかなら、少なくとも一つ非退化な安定円があるという通常の仮定からも同じ二種類に限定される。
束縛運動が一つもない斥力の系を、空虚な「すべて閉じる」の例に含めないために、安定円の存在は大切である。

### 同じエネルギーで異なる動径運動

中心力なら常に遠心障壁が落下を防ぐとは限らない。
原点に近づいたときに $V$ と $\ell^2/(2\mu r^2)$ のどちらが優勢かを比較する。

:::example[逆二乗ポテンシャルへの落下]
$V=-g/(2r^2)$、$g>0$ とする。
$U_\ell=(\ell^2/\mu-g)/(2r^2)$ である。
$\ell^2>\mu g$ なら原点近くが禁止領域となり、$E>0$ の粒子は

$$
r_{\min}=\sqrt{(\ell^2/\mu-g)/(2E)}
$$

で戻る。
$\ell^2<\mu g$ では $E=0$ の内向き解は $\dot r=-a/r$、$a=\sqrt{g/\mu-\ell^2/\mu^2}$ に従い、
$r(t)^2=r_0^2-2at$ で有限時刻に原点へ至る。
境界の $\ell^2=\mu g$ では有効ポテンシャルは平坦で $\ddot r=0$。
Bertrand の定理で非退化な安定円を仮定した理由がこの例に現れる。
:::

:::example[転回点から軌道を復元する積分]
二つの単純転回点 $r_-<r_+$ があれば、一往復の時間と角度は別々の積分

$$
T_r=2\int_{r_-}^{r_+}\frac{\dd r}{\sqrt{2(E-U_\ell)/\mu}},\qquad
\Phi_r=2\int_{r_-}^{r_+}
\frac{\ell\,\dd r}{\mu r^2\sqrt{2(E-U_\ell)/\mu}}
$$

で得られる。
両者は同じ積分要素に角速度を掛けるかどうかで異なる。
例えば $\Phi_r/(2\pi)=3/2$ なら動径二往復の間に方位が三周して軌道は閉じる。
単純根の近くの被積分関数は $|r-r_\pm|^{-1/2}$ の桁で積分可能である。
これに対し不安定円の二重根に近づく軌道では $|r-r_0|^{-1}$ が現れ、到達時間は対数的に発散する。
:::

### 角運動量を選んで有効ポテンシャルの分岐を調べる

極小をもつ曲線を図に一つ描いただけでは、隣の角運動量でも同じ型の運動があるか分からない。半径を入力して円運動を保つ角運動量を出力する関数 $C(r)=\mu r^3V'(r)$ を作ると、複数の円の出現と消滅を一つの図で判定できる。実際に

$$
U_\ell'=\frac{C(r)-\ell^2}{\mu r^3},\qquad
U_\ell''(r_*)=\frac{C'(r_*)}{\mu r_*^3}\quad(C(r_*)=\ell^2)
\label{eq-am-circle-branch}
$$

である。正の水平線 $\ell^2$ が $C$ の増加する枝と交わると極小、減少する枝なら極大である。接するとき二次近似による判定は失敗する。この方法では $C<0$ の領域には円が存在しないことも同時に分かる。以下では $\mu=1\,\mathrm{kg}$ とし、$V$ は J、$r$ は m の数値で記す。

:::example[対数型の井戸]
$V(r)=\ln{\left(r \right)}$ とする。$r_*=2\,\mathrm m$ に円があるかを調べる。力は $F=-V'=- \frac{1}{r}$ だから、この半径で必要な角運動量の二乗と有効ポテンシャルは

$$
C(r)=r^{2},\quad \ell_*^2=4,\quad U_{\ell_*}(r)=\frac{r^{2} \ln{\left(r \right)} + 2}{r^{2}}.
$$

従って接線速度は $v_\phi=\ell_*/(\mu r_*)=1\,\mathrm{m/s}$、円のエネルギーは $U_{\ell_*}(r_*)=\frac{1 + 2 \ln{\left(2 \right)}}{2}\,\mathrm J$。半径をずらすときは今求めた $\ell_*$ を固定する。微分した後に初期半径を代入すれば

$$
U'_{\ell_*}=\frac{\left(r - 2\right) \left(r + 2\right)}{r^{3}},\quad U''_{\ell_*}(r_*)=\frac{1}{2}\,\mathrm{N/m}.
$$

従って安定で、近点から次の近点までの角は $2\pi\Omega/\kappa=\sqrt{2} \pi$ である。対数の引数は基準長 $1\,\mathrm m$ に対する比である。円速度はどの半径でも一定になるが、動径周期と公転周期は一致しない。
:::

:::example[斥力性の芯をもつ井戸]
$V(r)=- \frac{4 r - 1}{r^{2}}$ とする。$r_*=1\,\mathrm m$ に円があるかを調べる。力は $F=-V'=- \frac{2 \left(2 r - 1\right)}{r^{3}}$ だから、この半径で必要な角運動量の二乗と有効ポテンシャルは

$$
C(r)=2 \left(2 r - 1\right),\quad \ell_*^2=2,\quad U_{\ell_*}(r)=- \frac{2 \left(2 r - 1\right)}{r^{2}}.
$$

従って接線速度は $v_\phi=\ell_*/(\mu r_*)=\sqrt{2}\,\mathrm{m/s}$、円のエネルギーは $U_{\ell_*}(r_*)=-2\,\mathrm J$。微分した後に初期半径を代入すれば

$$
U'_{\ell_*}=\frac{4 \left(r - 1\right)}{r^{3}},\quad U''_{\ell_*}(r_*)=4\,\mathrm{N/m}.
$$

従って安定で、近点から次の近点までの角は $2\pi\Omega/\kappa=\sqrt{2} \pi$ である。芯の斥力が勝つ小半径には引力の円はない。外側の安定な枝では、同じ動径振動の間に進む角が Kepler 楕円の場合からずれる。
:::

:::example[急に減衰する引力が作る不安定円]
$V(r)=- \frac{1}{r^{4}}$ とする。$r_*=2\,\mathrm m$ に円があるかを調べる。力は $F=-V'=- \frac{4}{r^{5}}$ だから、この半径で必要な角運動量の二乗と有効ポテンシャルは

$$
C(r)=\frac{4}{r^{2}},\quad \ell_*^2=1,\quad U_{\ell_*}(r)=\frac{r^{2} - 2}{2 r^{4}}.
$$

従って接線速度は $v_\phi=\ell_*/(\mu r_*)=\frac{1}{2}\,\mathrm{m/s}$、円のエネルギーは $U_{\ell_*}(r_*)=\frac{1}{16}\,\mathrm J$。微分した後に初期半径を代入すれば

$$
U'_{\ell_*}=- \frac{\left(r - 2\right) \left(r + 2\right)}{r^{5}},\quad U''_{\ell_*}(r_*)=- \frac{1}{8}\,\mathrm{N/m}.
$$

従って不安定で、変位の増大解は $\exp(\gamma t)$ に比例し、$\gamma=\frac{\sqrt{2}}{4}\,\mathrm{s^{-1}}$ である。この円を内側にずらした運動と外側にずらした運動は分かれる。負の曲率をもつ円に実数の動径周期を割り当ててはいけない。
:::

:::example[ばねと一定の外向き力]
$V(r)=\frac{r \left(r - 2\right)}{2}$ とする。$r_*=2\,\mathrm m$ に円があるかを調べる。力は $F=-V'=1 - r$ だから、この半径で必要な角運動量の二乗と有効ポテンシャルは

$$
C(r)=r^{3} \left(r - 1\right),\quad \ell_*^2=8,\quad U_{\ell_*}(r)=\frac{r^{4} - 2 r^{3} + 8}{2 r^{2}}.
$$

従って接線速度は $v_\phi=\ell_*/(\mu r_*)=\sqrt{2}\,\mathrm{m/s}$、円のエネルギーは $U_{\ell_*}(r_*)=1\,\mathrm J$。微分した後に初期半径を代入すれば

$$
U'_{\ell_*}=\frac{\left(r - 2\right) \left(r^{3} + r^{2} + 2 r + 4\right)}{r^{3}},\quad U''_{\ell_*}(r_*)=\frac{5}{2}\,\mathrm{N/m}.
$$

従って安定で、近点から次の近点までの角は $2\pi\Omega/\kappa=\frac{2 \sqrt{5} \pi}{5}$ である。力が零になる半径そのものは、正の速さの円軌道ではない。遠心項と実際の力を別々に評価して向きも確かめられる。
:::

## Kepler 問題の厳密解

### 軌道の幾何

$V=-k/r$ とする。
重力なら $k=Gm_1m_2$ で $k/\mu=GM$ である。
Binet 方程式を積分して

$$
r=\frac{p}{1+e\cos(\phi-\phi_0)},\qquad
p=\frac{\ell^2}{\mu k},\qquad e^2=1+\frac{2E\ell^2}{\mu k^2}
\label{eq-am-kepler-conic}
$$

を得る。
最後の式は $E=\ell^2[(u')^2+u^2]/(2\mu)-ku$ に解を代入すれば得られ、$e\ge0$ を :term[離心率]{en="eccentricity" yomi="りしんりつ"} という。
$e>0$ なら角の原点を近点にとると $r+ex=p$ であり、焦点からの距離と直線 $x=p/e$ までの距離の比が $e$ である。$e=0$ では軌道式は直接 $r=p$ を与える。

| 条件 | 軌道 | 近点と遠方の挙動 |
| --- | --- | --- |
| $e=0$ | 円 | 半径 $p$ |
| $0<e<1$ | 楕円 | $r_-=p/(1+e),\ r_+=p/(1-e)$ |
| $e=1$ | 放物線 | $E=0$、無限遠で速さが零に近づく |
| $e>1$ | 双曲線の一枝 | $E>0$、無限遠で速さは $\sqrt{2E/\mu}$ |

楕円の長半径は $a=(r_++r_-)/2=p/(1-e^2)=-k/(2E)$。
$r=p-ex$ を二乗して平方完成すると短半径は $a\sqrt{1-e^2}$ である。
楕円面積を面積速度で割れば

$$
T=\frac{2\pi\mu a^2\sqrt{1-e^2}}{\ell}
=2\pi\sqrt{\frac{\mu a^3}{k}} .
$$

同じ $a$ なら離心率によらず同じ周期をもつ。

```tikz caption="楕円の右焦点を力の中心にとった例。近点までの距離と遠点までの距離の和は長軸長に等しい。"
\draw[thick] (0,0) ellipse (2.5 and 2);
\draw[dashed] (-2.65,0) -- (2.7,0);
\fill (1.5,0) circle (1.8pt) node[above left] {$F$};
\fill (0,0) circle (1pt) node[above] {$C$};
\draw[<->] (1.5,-.35) -- node[below] {$r_-$} (2.5,-.35);
\draw[<->] (-2.5,-.35) -- node[below] {$r_+$} (1.5,-.35);
\draw[->,thick] (1.5,0) -- (1.25,1.732) node[midway,right] {$\mathbf r$};
\fill (1.25,1.732) circle (1.5pt);
```

### 時刻から位置を求める

離心近点角 $\mathcal E$ を $x=a(\cos\mathcal E-e)$、$y=a\sqrt{1-e^2}\sin\mathcal E$ で定義する。
$r=a(1-e\cos\mathcal E)$ である。
面積微分 $2\,\dd A=x\,\dd y-y\,\dd x$ に代入し、近点時刻 $t_p$ から積分すると

$$
n(t-t_p)=\mathcal E-e\sin\mathcal E,\qquad n=2\pi/T .
$$

右辺の微分は $1-e\cos\mathcal E>0$ だから $\mathcal E$ は一意に決まる。
例えば $a=4\,\mathrm m,\ e=1/2,\ k/\mu=4\,\mathrm{m^3/s^2}$ なら $n=1/4\,\mathrm{s^{-1}}$。
$\mathcal E=\pi/2$ で $r=4\,\mathrm m$、近点からの時間は $(2\pi-2)\,\mathrm s\simeq4.28\,\mathrm s$である。

### 正のエネルギーと零のエネルギーでの飛行時間

楕円では遠点を通過して同じ周期を繰り返すが、開いた軌道に周期はない。
その代わり最近接時刻を原点にして、近点からどれだけ時間が経ったかを求める。

双曲線で $a_h=k/(2E)>0$、$e>1$ とし、

$$
x=a_h(e-\cosh H_h),\qquad
y=a_h\sqrt{e^2-1}\sinh H_h
$$

と置く。
$r=a_h(e\cosh H_h-1)$ であり、近点は $H_h=0$。
面積速度を使うと

$$
x\,\dd y-y\,\dd x
=a_h^2\sqrt{e^2-1}(e\cosh H_h-1)\dd H_h .
$$

ここで $\ell^2=\mu k a_h(e^2-1)$ を用いて積分すると

$$
t-t_p=\sqrt{\frac{\mu a_h^3}{k}}\,(e\sinh H_h-H_h).
$$

$t$ は負から正へ単調に増え、$|H_h|$ が大きいと $r/|t|\to\sqrt{2E/\mu}$。
従って遠方の速さはエネルギーの漸近値とも一致する。

放物線では $e=1$、$D=\tan[(\phi-\phi_0)/2]$ と置けば $r=p(1+D^2)/2$ となる。
$\dd\phi=2\,\dd D/(1+D^2)$ と $\dd t=\mu r^2\,\dd\phi/\ell$ より

$$
t-t_p=\frac12\sqrt{\frac{\mu p^3}{k}}\left(D+\frac{D^3}{3}\right).
$$

例えば $p=2\,\mathrm m,\ k/\mu=2\,\mathrm{m^3/s^2}$ の放物線で $\phi-\phi_0=\pi/2$ に達する時間は $D=1$ より $4/3\,\mathrm s$ である。

### 観測位置での速さと円軌道への接続

エネルギー式で速度を解くと、楕円では

$$
v^2=\frac{k}{\mu}\left(\frac2r-\frac1a\right).
$$

同じ位置 $r$ でも長半径が違えば速さが違う。
例えば一つの楕円の近点と遠点で、速度はともに動径に垂直である。
従って $\ell=\mu r_-v_-=\mu r_+v_+$ から $v_-/v_+=(1+e)/(1-e)$。
$e=1/2$ なら距離が三倍になる間に速さは三分の一になる。
$e\to0$ とすればどの位置でも $v^2=k/(\mu a)$ となり、円運動の向心力条件に連続につながる。

### 近点時刻の数値計算と精度

観測時刻から位置を戻すには、平均近点角 $M_a=n(t-t_p)$ が与えられている。
$f(z)=z-e\sin z-M_a$ の根は

$$
M_a-e\le z\le M_a+e
$$

にある。
これは $z-M_a=e\sin z$ から従う。
さらに $f'(z)=1-e\cos z\ge1-e$ より、近似値 $\hat z$ に対する角の誤差は

$$
|\hat z-\mathcal E|\le |f(\hat z)|/(1-e)
$$

で抑えられる。
Newton 法の更新は

$$
z_{\rm next}=z-\frac{z-e\sin z-M_a}{1-e\cos z}.
$$

高離心率で分母が小さい場合にも、上の区間で符号を保った二分法を併用すれば一意な根を追える。

:::example[時刻から位置・速度を一緒に決める]
$e=1/2$、$M_a=\pi/2-1/2$ とすれば厳密解は $\mathcal E=\pi/2$ である。
速度を得るには時刻方程式を微分し、

$$
\dot{\mathcal E}=\frac{n}{1-e\cos\mathcal E}
$$

を $x,y$ の媒介表示に入れる。
従って

$$
\begin{aligned}
\dot x&=-\frac{an\sin\mathcal E}{1-e\cos\mathcal E},\\
\dot y&=\frac{an\sqrt{1-e^2}\cos\mathcal E}{1-e\cos\mathcal E}.
\end{aligned}
$$

$a=4\,\mathrm m,\ n=1/4\,\mathrm{s^{-1}}$ なら、上の時刻で位置は $(-2,2\sqrt3)\,\mathrm m$、速度は $(-1,0)\,\mathrm{m/s}$。
$x\dot y-y\dot x=2\sqrt3\,\mathrm{m^2/s}$ から、面積速度も直接確かめられる。
:::

### 角運動量零の退化した極限

Binet 方程式では $\phi$ を独立変数に替えるため、$\ell=0$ の軌道を扱えない。
半径 $R$ で静止した粒子が引力 $-k/r^2$ によって落下する場合を、エネルギーだけで解いてみよう。
$E=-k/R$ より

$$
\dot r=-\sqrt{\frac{2k}{\mu}\left(\frac1r-\frac1R\right)} .
$$

$r=R\cos^2\eta$ と置けば $0\le\eta\le\pi/2$ に対して $\dd r=-2R\sin\eta\cos\eta\,\dd\eta$。
積分すると

$$
t=\sqrt{\frac{\mu R^3}{2k}}
  \left(\eta+\tfrac12\sin2\eta\right),\qquad
t_{\rm fall}=\frac\pi2\sqrt{\frac{\mu R^3}{2k}} .
$$

これは初点を遠点とする $a=R/2$ の退化楕円の半周期に一致する。
ただし $r=0$ は点粒子ポテンシャルの特異点であり、衝突の後にどの軌道を接続するかはこの方程式だけでは指定されない。

### 衝突近くの時計を変える Levi-Civita の正則化

角運動量が小さい束縛軌道では近点付近だけ速さが大きくなり、等しい時間刻みで軌道を追いにくい。平面の位置を複素数 $z=x+iy$、$\gamma=k/\mu$ とすれば $\ddot z=-\gamma z/|z|^3$。原点を除く軌道の区間で $z=w^2$ の平方根を連続に一枚選び、時計を $\dd t=|w|^2\dd s$ と変える。$w$ の単位は $\mathrm{m^{1/2}}$、$s$ は $\mathrm{s/m}$ である。プライムを $s$ 微分とすると
$$
\dot z=\frac{2w'}{\bar w},\qquad
\ddot z=\frac{2(w''\bar w-|w'|^2)}{w\bar w^3}.
$$
力の式の右辺は $-\gamma/(w\bar w^3)$ だから $2w''\bar w=2|w'|^2-\gamma$。一方、$\varepsilon=E/\mu$ に対するエネルギー式に $|w|^2$ を掛けると $\varepsilon |w|^2=2|w'|^2-\gamma$。従って $w''=\varepsilon w/2$ となる。負のエネルギーでは二つの実成分が同じ振動数で振動し、正なら双曲線関数、零なら直線で進む。この線形方程式に初期位置と任意の初期速度を入れてよいわけではない。元の $k$ を保つには $2|w'|^2-\varepsilon |w|^2=\gamma$ が必要で、線形方程式を掛けて微分すればこの条件が保たれる。

例えば $\gamma=2\,\mathrm{m^3/s^2}$ とし、$z(0)=1\,\mathrm m$ から静かに落とす。SI の数値で $\varepsilon=-2$、$w(0)=1,\ w'(0)=0$ なので $w=\cos s$。従って $r=\cos^2s,\ t=s/2+\sin(2s)/4$。$s=\pi/2$ で $r=0$、物理的な到達時間は $\pi/4\,\mathrm s$ になる。$w'$ は有限だが $\dd t/\dd s$ は零となるので $\dot r$ の発散とは矛盾しない。正則化した変数で衝突後へ解を続けることと、実物が接触したあとの反発や合体を指定することは異なる。後者には有限半径と短距離相互作用の模型が要る。

### 一時刻の状態から軌道を空間に配置する

円錐曲線の式は軌道面をすでに知っているときの表示である。望遠鏡や計算からは空間の位置と速度が一組だけ与えられることも多い。基準面を $xy$ 面とし、単位質量あたりの引力定数を $\gamma=k/\mu$ と書く。$\vb h=\vb r\times\vb v$ が零でなければ平面の正の法線が決まる。近点方向を表すベクトルは $\vb e=\vb v\times\vb h/\gamma-\vb r/r$ である。これが一定であることは次節で運動方程式から証明する。ここでは $\vb h\cdot(\vb r\times\vb v)=h^2$ より $\vb e\cdot\vb r=h^2/\gamma-r$ を使う。この式を $r$ で割れば、$\vb e$ から測った角 $\nu$ に対して $r=(h^2/\gamma)/(1+e\cos\nu)$ が得られる。

面の傾きは $\cos i=h_z/h$、昇交点方向は $\vb n=\vb e_z\times\vb h$ である。この向きで速度の $z$ 成分は正になる。$\vb n\ne0$ なら $\operatorname{atan2}(n_y,n_x)$ で昇交点の方位角を求める。面内の二つの方向 $\vb a,\vb b$ の有向角も、$\operatorname{atan2}(\hat{\vb h}\cdot(\vb a\times\vb b),\vb a\cdot\vb b)$ から求まる。$\vb n$ から $\vb e$ までが近点引数、$\vb e$ から $\vb r$ までが真近点角である。円では近点が定まらず、基準面内の軌道では昇交点も定まらない。その場合でも位置・速度を使う表現は正則なままである。

:::example[昇交点を近点で通る楕円]
相対位置 $\vb r=\left(\begin{matrix}2 & 0 & 0\end{matrix}\right)^\top\,\mathrm m$、相対速度 $\vb v=\left(\begin{matrix}0 & 1 & 1\end{matrix}\right)^\top\,\mathrm{m/s}$、$\gamma=3\,\mathrm{m^3/s^2}$ とする。まず空間の情報を使わず、単位質量あたりのエネルギー $\varepsilon=v^2/2-\gamma/r=- \frac{1}{2}\,\mathrm{m^2/s^2}$ を求める。次に軌道の向きを求めると

$$
\vb h=\left(\begin{matrix}0 & -2 & 2\end{matrix}\right)^\top,\qquad \vb n=\left(\begin{matrix}2 & 0 & 0\end{matrix}\right)^\top,\qquad \vb e=\left(\begin{matrix}\frac{1}{3} & 0 & 0\end{matrix}\right)^\top.
$$

この結果から $p=h^2/\gamma=\frac{8}{3}\,\mathrm m$、$e=\frac{1}{3}$。不変量の関係 $1+2\varepsilon h^2/\gamma^2=\frac{1}{9}$ は直接計算した $e^2$ と一致する。近点距離は $p/(1+e)=2\,\mathrm m$、面法線と鉛直との角は $\cos i=\frac{\sqrt{2}}{2}$ で定まる。傾きは $45^\circ$ であり、この瞬間には高さが増えている。近点引数を選ぶとき、降交点を出発方向に選ぶと半周の誤りとなる。
:::

:::example[近点を過ぎた放物線]
相対位置 $\vb r=\left(\begin{matrix}1 & 0 & 0\end{matrix}\right)^\top\,\mathrm m$、相対速度 $\vb v=\left(\begin{matrix}1 & 1 & 0\end{matrix}\right)^\top\,\mathrm{m/s}$、$\gamma=1\,\mathrm{m^3/s^2}$ とする。まず空間の情報を使わず、単位質量あたりのエネルギー $\varepsilon=v^2/2-\gamma/r=0\,\mathrm{m^2/s^2}$ を求める。次に軌道の向きを求めると

$$
\vb h=\left(\begin{matrix}0 & 0 & 1\end{matrix}\right)^\top,\qquad \vb n=\left(\begin{matrix}0 & 0 & 0\end{matrix}\right)^\top,\qquad \vb e=\left(\begin{matrix}0 & -1 & 0\end{matrix}\right)^\top.
$$

この結果から $p=h^2/\gamma=1\,\mathrm m$、$e=1$。不変量の関係 $1+2\varepsilon h^2/\gamma^2=1$ は直接計算した $e^2$ と一致する。近点距離は $p/(1+e)=\frac{1}{2}\,\mathrm m$、面法線と鉛直との角は $\cos i=1$ で定まる。この軌道は基準面内にあるため昇交点の方位は任意である。直交成分からは近点と観測点が同じ方向でないことも分かる。
:::

:::example[基準面へ沈む双曲線]
相対位置 $\vb r=\left(\begin{matrix}2 & 0 & 0\end{matrix}\right)^\top\,\mathrm m$、相対速度 $\vb v=\left(\begin{matrix}0 & 1 & -2\end{matrix}\right)^\top\,\mathrm{m/s}$、$\gamma=2\,\mathrm{m^3/s^2}$ とする。まず空間の情報を使わず、単位質量あたりのエネルギー $\varepsilon=v^2/2-\gamma/r=\frac{3}{2}\,\mathrm{m^2/s^2}$ を求める。次に軌道の向きを求めると

$$
\vb h=\left(\begin{matrix}0 & 4 & 2\end{matrix}\right)^\top,\qquad \vb n=\left(\begin{matrix}-4 & 0 & 0\end{matrix}\right)^\top,\qquad \vb e=\left(\begin{matrix}4 & 0 & 0\end{matrix}\right)^\top.
$$

この結果から $p=h^2/\gamma=10\,\mathrm m$、$e=4$。不変量の関係 $1+2\varepsilon h^2/\gamma^2=16$ は直接計算した $e^2$ と一致する。近点距離は $p/(1+e)=2\,\mathrm m$、面法線と鉛直との角は $\cos i=\frac{\sqrt{5}}{5}$ で定まる。観測位置は正の $x$ 軸上でも速度の鉛直成分は負であり、ここは降交点である。近点距離は正で、遠点は存在しない。
:::

:::example[近点を定められない傾斜円]
相対位置 $\vb r=\left(\begin{matrix}1 & 0 & 0\end{matrix}\right)^\top\,\mathrm m$、相対速度 $\vb v=\left(\begin{matrix}0 & 1 & 1\end{matrix}\right)^\top\,\mathrm{m/s}$、$\gamma=2\,\mathrm{m^3/s^2}$ とする。まず空間の情報を使わず、単位質量あたりのエネルギー $\varepsilon=v^2/2-\gamma/r=-1\,\mathrm{m^2/s^2}$ を求める。次に軌道の向きを求めると

$$
\vb h=\left(\begin{matrix}0 & -1 & 1\end{matrix}\right)^\top,\qquad \vb n=\left(\begin{matrix}1 & 0 & 0\end{matrix}\right)^\top,\qquad \vb e=\left(\begin{matrix}0 & 0 & 0\end{matrix}\right)^\top.
$$

この結果から $p=h^2/\gamma=1\,\mathrm m$、$e=0$。不変量の関係 $1+2\varepsilon h^2/\gamma^2=0$ は直接計算した $e^2$ と一致する。近点距離は $p/(1+e)=1\,\mathrm m$、面法線と鉛直との角は $\cos i=\frac{\sqrt{2}}{2}$ で定まる。離心率ベクトルが零なので、その方向の正規化はできない。位置は昇交点からの角で指定できる。
:::

### 楕円上の二点を位置と速度で結ぶ

軌道が決まっていても、座標軸を毎回近点の方向に取り直すのは不便である。そこで初期の二ベクトルを基底の代わりに使う。楕円の長半径を $a$、平均運動を $n=\sqrt{k/(\mu a^3)}$、初期離心近点角を $u_0$、その増分を $d$ とする。近点軸とその垂線に沿った単位ベクトルを $\vb i,\vb j$ とすれば

$$
\vb r(u)=a[(\cos u-e)\vb i+\sqrt{1-e^2}\sin u\,\vb j],\qquad
\vb v(u)=\frac{na^2}{r(u)}[-\sin u\,\vb i+\sqrt{1-e^2}\cos u\,\vb j].
$$

非退化な楕円では $\vb r_0,\vb v_0$ が独立なので $\vb r=f\vb r_0+g\vb v_0$ と書ける。両辺と $\vb v_0$ の外積をとれば $f=(\vb r\times\vb v_0)\cdot\hat{\vb h}/h$。同様に $g=(\vb r_0\times\vb r)\cdot\hat{\vb h}/h$ となる。上の媒介表示で $h=na^2\sqrt{1-e^2}$ を使い、$\cos(u_0+d)\cos u_0+\sin(u_0+d)\sin u_0=\cos d$ を代入すると

$$
f=1-\frac a{r_0}(1-\cos d),\qquad
g=\frac{\sin d-e[\sin(u_0+d)-\sin u_0]}{n}
  =\Delta t-\frac{d-\sin d}{n}.
$$

最後の等号は Kepler 方程式を二時刻で引いて得る。$f$ は無次元、$g$ は時間である。終点の速度も初期ベクトルで書ける。初期時刻を固定し、$\dot d=na/r$ を掛けて上の式を微分すれば

$$
\dot f=-\frac{na^2}{rr_0}\sin d,\qquad
\dot g=1-\frac a{r}(1-\cos d),\qquad
\vb v=\dot f\,\vb r_0+\dot g\,\vb v_0 .
$$

面積速度の保存から $f\dot g-\dot f g=1$ である。これは四つの係数を別々に丸めたときの検査にも使える。ただし $d$ 自体は初期条件と飛行時間に依存するため、これらの係数を任意の異なる軌道に使い回すことはできない。

:::example[近点から四分の一だけ離心近点角を進める]
$a=2\,\mathrm m,e=1/2,n=1\,\mathrm{s^{-1}}$ とする。近点から $d=\pi/2$ だけ進めば $r_0=1\,\mathrm m,r=2\,\mathrm m,\Delta t=(\pi-1)/2\,\mathrm s$。一周期の四分の一より短いのは、近点で速く動くためである。係数は $f=-1,g=1/2\,\mathrm s,\dot f=-2\,\mathrm{s^{-1}},\dot g=0$。従って $\vb r_0=(1,0)\,\mathrm m,\vb v_0=(0,2\sqrt3)\,\mathrm{m/s}$ から $\vb r=(-1,\sqrt3)\,\mathrm m,\vb v=(-2,0)\,\mathrm{m/s}$ が出る。初期と終点の外積はともに $2\sqrt3\,\mathrm{m^2/s}$ である。$\dot g=0$ は終点速度に初期接線方向の成分がないことを意味し、速度が零であることを意味しない。
:::

この手順は時間を逆にしても使える。終点から出発点へ戻る係数は行列式が一の行列の逆であるから $f_{\rm back}=\dot g,\ g_{\rm back}=-g$ である。上の例で逆向きに代入すると $\vb r_0=-(1/2)\vb v$ となり、向きも単位も元の位置に戻る。一般の観測データでは外積とエネルギーを併せて確かめると、時間の符号の取り違えと長さ尺度の取り違えを区別しやすい。

補助円に対応する点を $P'$、楕円上の点を $P$ とする。離心近点角は中心 $C$ から $P'$ を見た角であり、真近点角は焦点 $F$ から $P$ を見た角である。図では $e=3/5,\cos u=3/5$ を選んだので $P$ は焦点の真上に来る。このとき真近点角は直角だが、$u$ は直角より小さい。円の中心から焦点までの水平距離 $ae$ を引いてから $\operatorname{atan2}$ を使うと、この違いを正しく再現できる。

```tikz caption="楕円と補助円。同じ離心近点角に対し、楕円上の高さだけを短軸と長軸の比で縮める。"
\draw[->] (-2.6,0)--(2.8,0) node[right] {$x$};
\draw[dashed] (0,0) circle (2);
\draw[thick] (0,0) ellipse (2 and 1.6);
\draw[dashed] (1.2,0)--(1.2,1.6);
\draw (0,0)--(1.2,1.6);
\draw[thick] (1.2,0)--(1.2,1.28);
\fill (1.2,0) circle (1.6pt) node[below] {$F$};
\fill (0,0) circle (1.3pt) node[below] {$C$};
\fill (1.2,1.28) circle (1.6pt) node[below left] {$P$};
\fill (1.2,1.6) circle (1.3pt) node[above right] {$P'$};
\draw (.46,0) arc[start angle=0,end angle=53.13,radius=.46];
\node at (.66,.29) {$u$};
```

## Runge–Lenz ベクトル

### 保存量から近点の方向を読む

エネルギーと角運動量だけでは、軌道を平面内で回したときに近点がどちらを向くかは決まらない。
$\vb p=\mu\dot{\vb r}$ とし、方向まで指定する :term[Runge–Lenz ベクトル]{en="Runge–Lenz vector" yomi="るんげれんつべくとる"} を

$$
\vb A=\vb p\times\vb\ell-\mu k\hat{\vb r}
\label{eq-am-runge-lenz}
$$

と定める。
$\dot{\vb p}=-k\hat{\vb r}/r^2$ と三重積の公式から

$$
\dot{\vb p}\times\vb\ell=\frac{k}{r}(\vb p-\hat{\vb r}(\hat{\vb r}\cdot\vb p))
=\mu k\,\dot{\hat{\vb r}}
$$

となり、$\dot{\vb A}=0$ である。
また $\vb A\cdot\vb\ell=0$、$\vb A\cdot\vb r=\ell^2-\mu kr$ なので、$\vb A$ 方向からの角で軌道式を解けば :eqref[eq-am-kepler-conic] と一致し、
$|\vb A|=\mu k e$ となる。
円では $\vb A=0$ で近点方向は定まらない。
さらに二乗を展開すると $A^2=\mu^2 k^2+2\mu E\ell^2$ となるので、これらの保存量はすべてが独立というわけではない。

:::example[接線方向に放した位置]
$\vb r=r_0\vb e_x,\ \vb p=\mu v\vb e_y$ なら $\vb A=\mu(\mu v^2r_0-k)\vb e_x$ である。
非零の接線速度が円速度より大きければ出発点は近点であり、円速度より小さければ遠点である。脱出速度以上なら遠点のない開軌道になる。
例えば $v^2=3k/(2\mu r_0)$ なら $e=1/2$ で、遠点は $3r_0$ になる。
:::

### 運動量の先端が描く円

固定ベクトルがあると、位置だけでなく速度の描く曲線も決まる。
$\hat{\vb n}=\vb\ell/\ell$ として :eqref[eq-am-runge-lenz] の両辺に左から $\hat{\vb n}\times$ を作用させる。
$\vb p\perp\hat{\vb n}$ なので

$$
\vb p=\frac{\hat{\vb n}\times\vb A}{\ell}
+\frac{\mu k}{\ell}\hat{\vb n}\times\hat{\vb r}.
$$

第一項は一定、第二項は長さが一定で向きだけが回る。
従って運動量の先端が描く曲線（ホドグラフ）は、半径 $\mu k/\ell$ の円またはその一部である。
中心までの距離と半径の比が $e$ なので、楕円運動なら運動量の原点は円の内部、双曲線なら外部にある。

:::example[一つの初期状態から向きを決める]
数値を SI 単位で $\mu=k=1$、$\vb r=(1,0,0)$、$\vb p=(1/2,1,0)$ とする。
$\vb\ell=(0,0,1)$、$\vb p\times\vb\ell=(1,-1/2,0)$ だから $\vb A=(0,-1/2,0)$。
従って近点は負の $y$ 方向、$e=1/2$ である。
エネルギーは $E=(1/4+1)/2-1=-3/8$ で、$a=4/3$、$p=1$、$r_-=2/3,\ r_+=2$ と得られる。
出発点では動径速度が正で、近点を過ぎて遠点へ向かっていることまで読み取れる。
:::

### 逆二乗力を乱したときの近点の移動

近点がじわじわ回る軌道を一周ずつ測ると、理想的な逆二乗則からのずれを調べられる。$\vb p=\mu\dot{\vb r}$、$\vb A=\vb p\times\boldsymbol\ell-\mu k\hat{\vb r}$ とし、逆二乗引力に小さな中心力 $f(r)\hat{\vb r}$ を加える。中心力なので $\boldsymbol\ell$ は一定で、逆二乗力だけによる項は相殺する。従って $\dot{\vb A}=f(r)\hat{\vb r}\times\boldsymbol\ell=-\ell f(r)\hat{\boldsymbol\phi}$。非円の軌道で、無摂動の近点を $x$ 軸にとると $\vb A=\mu ke\,\vb e_x$。有向角の変化は一次まで $\dot\varpi=\dot A_y/(\mu ke)=-\ell f(r)\cos\nu/(\mu ke)$ である。$\dd t=\mu r^2\dd\nu/\ell$ を用いて一周積分すれば

$$
\Delta\varpi=-\frac1{ke}\int_0^{2\pi}r^2 f(r)\cos\nu\,\dd\nu,\qquad
r=\frac{p}{1+e\cos\nu}
\label{eq-am-apsis-perturbation}
$$

を得る。右辺がすでに摂動の一次なので、積分には無摂動楕円を使える。大きさの変化には $\sin\nu$ が掛かるため、一周積分は零になる。瞬間の離心率まで一定という結論ではない。円に近い極限では近点方向の測定が悪条件になるので、微小動径振動の周波数と照合する。

:::example[逆二乗力の係数だけを変える]
$f(r)=-h/r^{2}$、$h>0$ とする。$h$ の単位は $\mathrm{N\,m^{2}}$ である。基準の逆二乗力に対する比が近点でも小さい範囲を扱う。半通径を $p$ とすると、歳差の積分に現れる角度部分は

$$
(1+e\cos\nu)^{0}\cos\nu=\cos{\left(\nu \right)}.
$$

奇数乗の余弦は半周ずつの寄与が相殺する。偶数乗は $\cos^2\nu=(1+\cos2\nu)/2$ を繰り返し用いて定数項だけを積分する。従って一周の角積分は $0$、求める一次の角は $\Delta\varpi=(h/k)\left(0\right)$。$e=1/2$、無次元の強さ $h/(kp^{0})=10^{-3}$ とした例では $\Delta\varpi=0$ ラジアンとなる。これは有効な $k$ の取り替えであり、一周の歳差は生じない。初期状態から推定した離心率はずれても、楕円が回り続けることはない。
:::

:::example[逆三乗の付加引力]
$f(r)=-h/r^{3}$、$h>0$ とする。$h$ の単位は $\mathrm{N\,m^{3}}$ である。基準の逆二乗力に対する比が近点でも小さい範囲を扱う。半通径を $p$ とすると、歳差の積分に現れる角度部分は

$$
(1+e\cos\nu)^{1}\cos\nu=\left(e \cos{\left(\nu \right)} + 1\right) \cos{\left(\nu \right)}.
$$

従って一周の角積分は $\pi e$、求める一次の角は $\Delta\varpi=(h/k)\left(\frac{\pi}{p}\right)$。$e=1/2$、無次元の強さ $h/(kp^{1})=10^{-3}$ とした例では $\Delta\varpi=\frac{\pi}{1000}$ ラジアンとなる。逆三乗の力は Binet 方程式の線形項の係数を変える。従ってこの場合には厳密な軌道の周期を摂動展開して比較できる。
:::

:::example[逆四乗の付加引力]
$f(r)=-h/r^{4}$、$h>0$ とする。$h$ の単位は $\mathrm{N\,m^{4}}$ である。基準の逆二乗力に対する比が近点でも小さい範囲を扱う。半通径を $p$ とすると、歳差の積分に現れる角度部分は

$$
(1+e\cos\nu)^{2}\cos\nu=\left(e \cos{\left(\nu \right)} + 1\right)^{2} \cos{\left(\nu \right)}.
$$

従って一周の角積分は $2 \pi e$、求める一次の角は $\Delta\varpi=(h/k)\left(\frac{2 \pi}{p^{2}}\right)$。$e=1/2$、無次元の強さ $h/(kp^{2})=10^{-3}$ とした例では $\Delta\varpi=\frac{\pi}{500}$ ラジアンとなる。近点を通る短い時間に強く働く補正である。同じ長半径でも離心率を変えると半通径も変わり、歳差の大きさが変わる。
:::

:::example[逆五乗の付加引力]
$f(r)=-h/r^{5}$、$h>0$ とする。$h$ の単位は $\mathrm{N\,m^{5}}$ である。基準の逆二乗力に対する比が近点でも小さい範囲を扱う。半通径を $p$ とすると、歳差の積分に現れる角度部分は

$$
(1+e\cos\nu)^{3}\cos\nu=\left(e \cos{\left(\nu \right)} + 1\right)^{3} \cos{\left(\nu \right)}.
$$

従って一周の角積分は $\frac{3 \pi e \left(e^{2} + 4\right)}{4}$、求める一次の角は $\Delta\varpi=(h/k)\left(\frac{3 \pi \left(e^{2} + 4\right)}{4 p^{3}}\right)$。$e=1/2$、無次元の強さ $h/(kp^{3})=10^{-3}$ とした例では $\Delta\varpi=\frac{51 \pi}{16000}$ ラジアンとなる。三次の余弦まで出るため平均に離心率の二次補正が残る。円近傍だけを測った周波数比からこの依存性は得られない。
:::

## 散乱断面積

### 入射する環と出射する円錐

$V(r)\to0$ とし、相対速度の大きさが遠方で $v_\infty>0$ に近づく弾性散乱を重心系で観測する。以下では入射漸近線が存在し、外側の根 $r_m>0$ が単純で、そこから無限遠へ戻る軌道を扱う。
入射漸近線と力の中心との距離を :term[衝突径数]{en="impact parameter" yomi="しょうとつけいすう"} $b$ と呼ぶ。
$\ell=\mu b v_\infty$、$E=\mu v_\infty^2/2$ である。
最近接距離 $r_m$ は $E=U_\ell(r_m)$ の外側の根である。
最近接から無限遠までの方位角は

$$
\Phi=\int_{r_m}^\infty\frac{b\,\dd r}{r^2\sqrt{1-V(r)/E-b^2/r^2}} .
$$

これは $\dot\phi/\dot r$ の積分である。
近点を方位角の原点にすると、入射と出射の漸近位置はそれぞれ $-\Phi,+\Phi$ の方向にある。入射速度は入射側の漸近位置と逆向きだから、主値の散乱角は $\chi=\arccos[-\cos(2\Phi)]$ である。通常の斥力散乱で $0\le\Phi\le\pi/2$ なら $\chi=\pi-2\Phi$ となる。引力で何周も回り込む場合には方位角の増分をこの主値へ戻す必要がある。
無力なら $b/r=\sin s$ と置いて $\Phi=\pi/2$ となる。

一様な入射束では環の面積 $2\pi b|\dd b|$ に比例する粒子が、立体角 $2\pi\sin\chi|\dd\chi|$ に入る。
よって :term[微分散乱断面積]{en="differential cross section" yomi="びぶんさんらんだんめんせき"} は

$$
\frac{\dd\sigma}{\dd\Omega}=\sum_j\frac{b_j}{\sin\chi}\left|\frac{\dd b_j}{\dd\chi}\right|
$$

となる。
和は同じ $\chi$ に至る枝についてとる。
単位は $\mathrm{m^2/sr}$ である。

### Rutherford 散乱

斥力 $V=+k/r$、$k>0$ では $u=(\mu k/\ell^2)(e\cos\phi-1)$ と書ける。
最近接を $\phi=0$ とすると、漸近線で $\cos\Phi=1/e$。
従って $\sin(\chi/2)=1/e$ である。
エネルギーから $e^2-1=2E\ell^2/(\mu k^2)=(2Eb/k)^2$ なので、

$$
b=\frac{k}{2E}\cot\frac{\chi}{2},\qquad
\frac{\dd\sigma}{\dd\Omega}=
\left(\frac{k}{4E}\right)^2\frac{1}{\sin^4(\chi/2)} .
\label{eq-am-rutherford}
$$

これが :term[Rutherford の散乱公式]{en="Rutherford formula" yomi="らざふぉーどのさんらんこうしき"} である。
SI で同符号の点電荷なら $k=q_1q_2/(4\pi\epsilon_0)$。
引力 $V=-k/r$ でも $b>0$ なら $\cos\Phi=-1/e$ なので $\chi=2\Phi-\pi$、従って $\sin(\chi/2)=1/e$ となる。同じ $|k|$ なら偏向する側は反対だが微分断面積は等しい。小角への発散は、任意に大きな $b$ でも力が及ぶという仮定に由来する。
検出限界 $\chi_{\min}>0$ を置けば、$\chi\ge\chi_{\min}$ の断面積は $\pi b(\chi_{\min})^2$ となる。

:::example[検出器の角度を変える]
$k/(2E)=2\,\mathrm{nm}$ のとき $\chi=\pi/2$ に対応する $b$ は $2\,\mathrm{nm}$、微分断面積は $4\,\mathrm{nm^2/sr}$ である。
$\chi=\pi/3$ に移すと $\sin(\chi/2)=1/2$ より $16\,\mathrm{nm^2/sr}$ になる。
:::

### 重心系と実験室系の角度

実験室で $m_2$ が初め静止し、$m_1$ が速さ $v$ で $x$ 方向に入射するとする。
散乱面を $xy$ 平面とし、相対速度が重心系で $\chi$ 曲がったとする。
弾性散乱だから相対速度の大きさは $v$ のままであり、実験室での出射速度は

$$
\vb v_1'=\frac{v}{M}\bigl(m_1+m_2\cos\chi,\ m_2\sin\chi\bigr).
$$

従って実験室の角 $\theta$ は、この二成分から象限も含めて定める。
$\alpha=m_1/m_2$ と置くと

$$
\cos\theta=\frac{\alpha+\cos\chi}{\sqrt{1+2\alpha\cos\chi+\alpha^2}},\quad
\frac{\dd\cos\theta}{\dd\cos\chi}
=\frac{1+\alpha\cos\chi}{(1+2\alpha\cos\chi+\alpha^2)^{3/2}} .
$$

断面積は同じ事象の個数を表すので、実験室の $\dd\sigma/\dd\Omega$ は重心系の値をこの導関数の絶対値で割り、同じ $\theta$ に対応する枝を足せばよい。
$\alpha\to0$ なら二つの角と立体角は一致する。
$\alpha=1$ なら $\theta=\chi/2$（出射粒子が止まる端点を除く）である。
例えば重心系の $90^\circ$ 散乱を、等質量では実験室の $45^\circ$ で観測する。

### 断面積から数え上げへ

粒子束 $j$ を単位面積・単位時間あたりの粒子数とし、互いに影響しない標的が $N_s$ 個ある薄い試料を用いる。
微小検出立体角 $\Delta\Omega$ 内で断面積がほぼ一定なら

$$
\dot N=jN_s\frac{\dd\sigma}{\dd\Omega}\Delta\Omega.
$$

例えば $j=10^{10}\,\mathrm{m^{-2}s^{-1}}$、$N_s=10^{12}$、
$\dd\sigma/\dd\Omega=10^{-20}\,\mathrm{m^2/sr}$、$\Delta\Omega=0.02\,\mathrm{sr}$ なら毎秒 $2$ 個となる。
試料が厚く多重散乱がある場合や角度範囲が広い場合にはこの単一角度での積は使えない。

:::example[引力による捕捉半径の増大]
半径 $R$ の球に接触した粒子をすべて数えるとする。
外では $V=-k/r$、入射速度は $v_\infty$ とする。
ちょうど $R$ で動径速度が零となる境界の軌道にエネルギー式を使うと

$$
\tfrac12\mu v_\infty^2=\frac{\mu v_\infty^2 b_{\max}^2}{2R^2}-\frac{k}{R},
\qquad
\sigma_{\rm hit}=\pi R^2\left(1+\frac{2k}{\mu Rv_\infty^2}\right).
$$

従って引力が遠方から軌道を曲げる分だけ幾何学的な面積より大きくなる。
脱出速度と入射速度が等しい場合は $\sigma_{\rm hit}=2\pi R^2$。
ここで数えたのは球への接触であり、散逸のない点粒子二体が無限遠から束縛状態に移る確率ではない。
:::

### 散乱角の小さい極限を力積で照合する

Rutherford の式を小角の極限で別の近似から確かめる。
$b$ が十分大きいとき、入射する相対位置を一次近似で $(v_\infty t,b)$ に置ける。
斥力 $+k/r^2$ の横成分をこの直線上で積分すれば

$$
\Delta p_y=\int_{-\infty}^{\infty}
\frac{kb\,\dd t}{(b^2+v_\infty^2t^2)^{3/2}}
=\frac{2k}{bv_\infty}.
$$

最後の積分は $s=v_\infty t/b$ とし、$s/\sqrt{1+s^2}$ の微分を使った。
縦力は時間の奇関数なのでこの次数で正味の力積を与えない。
従って

$$
\chi\simeq\frac{\Delta p_y}{\mu v_\infty}=\frac{k}{Eb}.
$$

厳密式 $\tan(\chi/2)=k/(2Eb)$ の一次項に一致する。
直線軌道での近似を使える条件は、この式で予測した $\chi$ 自体が十分小さいことである。

:::example[全断面積と運動量輸送で重みが違う]
前方へのごく小さな偏向は散乱事象としては一個でも、入射方向の運動量損失はわずかである。
$p_\infty=\mu v_\infty$ とすると一回の前方運動量減少は $p_\infty(1-\cos\chi)$。
そこで重み $1-\cos\chi$ を付けた輸送断面積を定義する。
$d=k/(2E)$、$\chi\ge\chi_{\min}>0$ について

$$
\begin{aligned}
\sigma_{\rm tr}
&=\int (1-\cos\chi)\frac{\dd\sigma}{\dd\Omega}\,\dd\Omega\\
&=4\pi d^2\int_{\sin(\chi_{\min}/2)}^1\frac{\dd u}{u}
=4\pi d^2\log\frac1{\sin(\chi_{\min}/2)} .
\end{aligned}
$$

ここでは $u=\sin(\chi/2)$ と置いた。
無遮蔽の逆二乗力では、小角発散は重みをかけても対数として残る。
例えば $\chi_{\min}=\pi/3$ なら $\sigma_{\rm tr}=4\pi d^2\log2$ であり、同じ範囲の全断面積 $3\pi d^2$ と違う。
:::

### 多価な偏向関数と虹散乱

偏向角 $\chi(b)$ が単調でない場合、一つの検出角へ複数の入射環が粒子を送る。異なる環の粒子数は加算されるので、各枝で計算した $b|\dd b/\dd\chi|/\sin\chi$ の和が微分断面積になる。$b_*>0$、$0<\chi_*<\pi$ で $\chi'(b_*)=0$、$\chi''(b_*)\ne0$ なら、$\chi-\chi_*\simeq\chi''(b_*)(b-b_*)^2/2$。寄与のある側での二枝の和は
$$
\frac{\dd\sigma}{\dd\Omega}\simeq
\frac{\sqrt{2}\,b_*}{\sin\chi_*\,\sqrt{|\chi''(b_*)|}}\,
|\chi-\chi_*|^{-1/2}.
$$

これを :term[虹散乱]{en="rainbow scattering" yomi="にじさんらん"} という。角幅 $\delta$ の検出器に入る粒子数は $\int_0^\delta x^{-1/2}\dd x=2\sqrt\delta$ に比例して有限である。ここでの特異性は入射面積から角度への写像の折返しによる。前方での Rutherford 断面積の発散とは、原因も角度積分の振る舞いも異なる。

### 実験室角に上限がある衝突

静止標的より入射粒子が重いとき、重心系で後方へ散乱されても実験室で前へ進むことがある。入射方向の単位ベクトルを $\vb i$ とする。散乱後の入射粒子の実験室速度の先端は、重心速度 $U\vb i$ を中心に、半径 $v_*$ の円を描く。入射質量 $m$、標的質量 $M$ なら $U/v_*=m/M$。$U>v_*$ の場合に原点からこの円へ二本の接線が引け、その角度から $\sin\theta_{\max}=v_*/U=M/m$ を得る。$0<\theta<\theta_{\max}$ の同じ観測方向に速い粒子と遅い粒子がある。円と半直線の交点を速度 $w$ について解けば
$$
w_\pm=U\cos\theta\pm\sqrt{v_*^2-U^2\sin^2\theta}.
$$
二枝のエネルギーまで測れば、角度だけでは区別できなかった二つの重心系散乱角を復元できる。例えば $m=2M$、入射速度 $3\,\mathrm{m/s}$ なら $U=2,\ v_*=1$（$\mathrm{m/s}$）で $\theta_{\max}=30^\circ$。前方では $w_+=3,\ w_-=1$、上限角では $w_+=w_-=\sqrt3$ となる。速度が大きい枝に限って集計した微分断面積と、角度だけを選んだ検出器の計数とは一般に異なる。

### 斥力による反跳粒子の実験室角

静止標的の反跳も、相対運動が分かれば求められる。二体の速度の逆変換から

$$
\vb v_2'=\frac{m_1v}{M}(1-\cos\chi,-\sin\chi).
$$

従って $0<\chi\le\pi$ では反跳速さが $2m_1v\sin(\chi/2)/M$、入射方向との角が $\vartheta=\pi/2-\chi/2$ である。反跳の散乱面内の横成分は入射粒子の横成分と逆向きである。円錐の立体角について $|\dd\cos\chi/\dd\cos\vartheta|=4\cos\vartheta$ だから、反跳を数える検出器では重心系の微分断面積をこの係数で掛ける。例えば $m_1=m_2$ で重心系の散乱角が $60^\circ$ なら、入射粒子は $30^\circ$、標的は反対側に $60^\circ$ で進む。両者の速度ベクトルは直交し、運動エネルギーの比は $3:1$ である。$\chi=0$ では標的が静止したままなので反跳方向を定義できない。

## 演習問題

::::exercise{level=basic}
$m_1=3\,\mathrm{kg},m_2=6\,\mathrm{kg}$ の二体が距離 $3\,\mathrm m$ にある。
換算質量と軽い方の重心距離を求めよ。

:::solution
$\mu=18/9=2\,\mathrm{kg}$。
軽い方は相対距離の $6/9$ を担うので $2\,\mathrm m$ である。
:::
::::

::::exercise{level=basic}
近点 $2\,\mathrm m$、遠点 $6\,\mathrm m$ の Kepler 楕円の $a,e,p$ を求めよ。

:::solution
和と差から $a=4\,\mathrm m$、$e=(6-2)/(6+2)=1/2$。
$p=a(1-e^2)=3\,\mathrm m$ となる。
:::
::::

::::exercise{level=basic}
$V=kr^4/4$ の円軌道の近くで動径振動の間に進む角を求めよ。

:::solution
$rV''/V'=3$ なので $\kappa/\Omega=\sqrt6$。
動径の周期は $2\pi/\kappa$ だから進む角は $2\pi/\sqrt6$ である。
:::
::::

::::exercise{level=basic}
Kepler 力で半径 $R$ における脱出速度と円速度の比を求めよ。

:::solution
無限遠で速さ零なら $\mu v_e^2/2=k/R$。
円では $\mu v_c^2/R=k/R^2$ だから $v_e/v_c=\sqrt2$ となる。
:::
::::

::::exercise{level=advanced}
半径 $R$ の弾性硬球による点粒子の微分断面積を求めよ。

:::solution
接触点の法線と入射線の鋭角を $\alpha$ とすれば $b=R\sin\alpha$、鏡面反射により $\chi=\pi-2\alpha$。
従って $b=R\cos(\chi/2)$ を断面積の式に入れると $R^2/4$。
全立体角で積分して $\pi R^2$ である。
:::
::::

::::exercise{level=advanced}
$V=-k/r+h/(2r^2)$、$k>0$ とし、$\ell>0$、$\ell^2+\mu h>0$ を満たす非円軌道を考える。
束縛軌道の近点移動を一動径周期につき求めよ。

:::solution
$F=-ku^2+hu^3$ を Binet の式に入れると $u''+(1+\mu h/\ell^2)u=\mu k/\ell^2$。
角方向の周期は $2\pi/\sqrt{1+\mu h/\ell^2}$。
一周との差は $2\pi[(1+\mu h/\ell^2)^{-1/2}-1]$ である。
:::
::::

::::exercise{level=advanced}
Rutherford 散乱で $\chi\ge\pi/3$ の断面積を積分と入射面積の両方で求めよ。

:::solution
$d=k/(2E)$ とおくと、$b=d\cot(\chi/2)$ より入射円板は $\pi d^2\cot^2(\pi/6)=3\pi d^2$。
角積分では $s=\sin(\chi/2)$ により $\sin\chi\,\dd\chi=4s\,\dd s$、よって $2\pi d^2\int_{1/2}^1s^{-3}\dd s=3\pi d^2$ で一致する。
:::
::::

::::exercise{level=advanced}
楕円軌道の時間平均 $\langle1/r\rangle$ を求め、平均運動エネルギーを計算せよ。

:::solution
$\dd t=(1-e\cos\mathcal E)\dd\mathcal E/n$ なので $\dd t/r=\dd\mathcal E/(na)$。
一周期で積分して $T$ で割ると $1/a$ である。
$\langle T_{\rm kin}\rangle=E+k\langle1/r\rangle=k/(2a)$ となる。
:::
::::

::::exercise{level=advanced}
引力と小さな調和力を合わせ $V=-k/r+hr^2/2$、$k,h>0$ とする。
半径 $r_0$ のほぼ円軌道で、近点移動を $hr_0^3/k$ の一次まで求めよ。

:::solution
$V'=k/r^2+hr$、$V''=-2k/r^3+h$ なので $\delta=hr_0^3/k$ として

$$
\frac{\kappa^2}{\Omega^2}=3+\frac{-2+\delta}{1+\delta}=\frac{1+4\delta}{1+\delta}.
$$

近点から次の近点までの角は $2\pi\Omega/\kappa=2\pi\sqrt{(1+\delta)/(1+4\delta)}$。
$\delta\ll1$ では $2\pi(1-3\delta/2)+O(\delta^2)$。
従って近点は一動径周期につき約 $3\pi\delta$ だけ逆行する。
この評価には小離心率と小摂動という二つの仮定がある。
:::
::::

::::exercise{level=advanced}
入射粒子と静止標的の質量が等しい。
重心系で一定の微分断面積 $C$ を実験室の入射粒子の角度で表示し、全断面積が一致することを示せ。

:::solution
$\theta=\chi/2$ より $0\le\theta<\pi/2$。
円錐の立体角を比較すると

$$
\frac{\dd\Omega_{\rm cm}}{\dd\Omega_{\rm lab}}
=\frac{\sin(2\theta)\,2\dd\theta}{\sin\theta\,\dd\theta}=4\cos\theta .
$$

従って実験室では前方で $4C\cos\theta$、後方では零である。
全断面積は

$$
2\pi\int_0^{\pi/2}4C\cos\theta\sin\theta\,\dd\theta=4\pi C
$$

であり重心系の全立体角積分に一致する。
$\chi=\pi$ では入射粒子が静止して方向を定義できないが、連続分布の積分値には影響しない。
:::
::::

::::exercise{level=advanced}
吸収中心の外側で $V=-C/r^4$、$C>0$ とする。無限遠から速度 $v$ で入る粒子について、遠心障壁を越えて内側へ進む断面積を求めよ。

:::solution
$\ell>0$ で $U=\ell^2/(2\mu r^2)-C/r^4$ の極値は $r_*^2=4\mu C/\ell^2$ にあり、$U''(r_*)=-8C/r_*^6<0$。障壁の高さは $\ell^4/(16\mu^2 C)$ である。$E=\mu v^2/2,\ell=\mu vb$ を代入すると $b^4<4C/E$ が通過条件になる。中心が障壁の十分内側にあり、到達した粒子を吸収するとすれば $\sigma=2\pi\sqrt{C/E}$。境界の軌道は不安定円に漸近するが、その一つの入射環の面積は零である。速度との積 $v\sigma=2\pi\sqrt{2C/\mu}$ が一定になるのは、ポテンシャルの指数が四だからである。
:::
::::

## 確認テスト

::::quiz{answer="1.5" tol="1e-6"}
ともに $3\,\mathrm{kg}$ の二粒子の換算質量を kg で答えよ。

:::explanation
$\mu=3\cdot3/(3+3)=1.5\,\mathrm{kg}$ である。
:::
::::

::::quiz
円軌道を固定した $\ell$ で動径方向に少し乱したとき、安定な微小振動を保証するのはどれか。

- [x] $U_\ell''>0$
- [ ] $V''<0$
- [ ] $U_\ell=0$

:::explanation
円では $U_\ell'=0$ であり、線形化した力は $-U_\ell''\,\delta r$ であり、曲率が正なら復元力が生じる。曲率が零なら二次近似から有限の振動数は得られない。
:::
::::

::::quiz{answer="8" tol="1e-6"}
同じ $k,\mu$ で長半径を4倍にすると周期は何倍か。

:::explanation
$T^2$ は $a^3$ に比例するため、$4^{3/2}=8$ 倍である。
:::
::::

::::quiz
円を除く Kepler 楕円で $\vb A$ はどこを向くか。

- [ ] 軌道面の法線
- [x] 力の中心から近点への方向
- [ ] 遠点での速度方向

:::explanation
$\vb A\cdot\vb r=\ell^2-\mu kr$ から、$\vb A$ と $\vb r$ が平行なときに半径が最小となる。
円では $\vb A=0$ である。
:::
::::

::::quiz
同じ逆二乗引力で同じエネルギーをもつ二つの非円楕円について、必ず等しい量をすべて選べ。

- [x] 長半径
- [x] 周期
- [ ] 近点の向き
- [ ] 角運動量の大きさ

:::explanation
$\mu,k$ が同じなら $a=-k/(2E)$ と $T=2\pi\sqrt{\mu a^3/k}$ が定まる。
離心率は同じとは限らず、$\ell^2=\mu ka(1-e^2)$ と近点方向はエネルギーだけから決まらない。
:::
::::

::::quiz{answer="0.25" tol="1e-6"}
同じ $k$ で入射の相対エネルギーを二倍にする。
一定の散乱角で Rutherford の微分断面積は何倍か。

:::explanation
$k/(4E)$ の二乗に比例するため $1/2^2=0.25$ 倍になる。
ここで角度を固定して比較したので、同じ衝突径数の軌道を比較した場合とは条件が異なる。
:::
::::

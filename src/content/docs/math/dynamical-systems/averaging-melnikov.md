---
title: 平均化法と Melnikov の方法
description: 周期摂動の累積と、分離軌道の安定・不安定枝の分裂を積分で評価する。
notation:
  - symbol: '\bar f'
    meaning: 周期にわたる f の時刻平均
  - symbol: '\varepsilon'
    meaning: 無次元の摂動パラメータ
  - symbol: 'M(\theta)'
    meaning: 位相 θ での分裂の一次係数
---

## この章の目標

- 周期摂動で速い位相と遅い変数を分ける。
- 平均化定理を近恒等変換と誤差評価から証明する。
- エネルギー差を用いて Melnikov 関数を導く。
- 強制 Duffing 振動子について横断的交差が起こる条件を計算する。

## 周期摂動系

時間に周期的に依存する小さな摂動でも、各周期で同じ向きに変位が蓄積すれば長い時間で有限の変化を起こす。そこで一周期内で往復する成分と、周期をまたいで蓄積する成分を分ける。

:::definition[周期摂動と平均]{#def-dyn-averaging}
$T>0$ とし、$f,g$ は時刻について $T$ 周期とする。状態の方程式と平均場を

$$
\dot x=\eps f(t,x)+\eps^2 g(t,x,\eps),\qquad
\bar f(x)=\frac1T\int_0^T f(t,x)\dd t
\label{eq-dyn-averaging-system}
$$

と書く。$\bar f$ を時刻平均、$\dot{\bar x}=\eps\bar f(\bar x)$ を一次の :term[平均化方程式]{en="averaged equation" yomi="へいきんかほうていしき"} と呼ぶ。$\eps>0$ に対して遅い時刻 $\tau=\eps t$ を使うと $\dd\bar x/\dd\tau=\bar f(\bar x)$ となる。
:::

周期が $T$ の外力は、位相が $\R/T\Z$ 上を速度 $1$ で進む自律系として扱える（:ref[prop-dyn-ode-flow] の後の拡張）。元の変数を時刻 $T$ ごとに観測すると、その位相が同じところへ戻る。平均は一周期の変位の一次係数に対応する。

:::example[周期内の変位と平均の変位]
$\dot x=\eps(2+\cos t)$ なら $x(t)=x_0+2\eps t+\eps\sin t$ である。平均化方程式との差は $\eps\sin t$ だけで、$t=2\pi n$ では零である。$t=1/\eps$ では平均による変位が $2$ である一方、周期内の変位は $\eps$ 以下である。
:::

:::example[状態に依存する平均]
$\dot x=\eps(-1+2\sin t)x$ なら $x(t)=x_0\exp[-\eps t+2\eps(1-\cos t)]$ である。$\bar f(x)=-x$ であり、一時的に係数が正になっても、一周期の倍率 $e^{-2\pi\eps}$ は $1$ 未満である。
:::

### 弱い外力を受ける振動子の座標

速く振動する変数自身が小さくなくても、振幅と位相がゆっくり変わる座標を選べば平均化が使える。例として、時刻と変位を無次元化した方程式

$$
\ddot q+\omega_0^2q
=\eps[-2\gamma\dot q-\alpha q^3+F\cos(\omega_d t)],
\qquad
\omega_d=\omega_0+\eps\Delta>0
$$

を考える。$\omega_0>0,\gamma>0$ とし、$\gamma,\alpha,F,\Delta$ を固定して $\eps\to0$ とする。駆動角 $\theta=\omega_d t$ に対し

$$
q=u\cos\theta+v\sin\theta,\qquad
\dot q=\omega_d(-u\sin\theta+v\cos\theta)
$$

で $(u,v)$ を状態座標として定義する。これは各時刻で可逆な回転とスケーリングであり、$\dot u,\dot v$ を零とする仮定ではない。二式を逆に解いて微分すると

$$
\dot u=-\frac{\sin\theta}{\omega_d}(\ddot q+\omega_d^2q),\qquad
\dot v=\frac{\cos\theta}{\omega_d}(\ddot q+\omega_d^2q)
$$

である。振動の主部が消え、$\ddot q+\omega_d^2q=(2\eps\omega_0\Delta+\eps^2\Delta^2)q+\eps[-2\gamma\dot q-\alpha q^3+F\cos\theta]$ なので右辺は $\eps$ のオーダーである。周期が $\eps$ に少し依存することを避けるには、先に $\theta$ を独立変数とし $\dd/\dd t=\omega_d\dd/\dd\theta$ とおけばよい。その場合は周期が常に $2\pi$ となる。

一周期で平方の正弦・余弦は $1/2$、その積は零に平均される。また $r^2=u^2+v^2$ とすると

$$
\langle q^3\sin\theta\rangle=\frac38r^2v,\qquad
\langle q^3\cos\theta\rangle=\frac38r^2u.
$$

例えば前者では $3u^2v\cos^2\theta\sin^2\theta$ と $v^3\sin^4\theta$ の積分だけが残る。平均はそれぞれ $3u^2v/8$ と $3v^3/8$ である。遅い時刻 $\tau=\eps t$ に戻して、一次平均は

$$
\begin{aligned}
\dv{u}{\tau}&=-\gamma u-(\Delta-\kappa r^2)v,\\
\dv{v}{\tau}&=(\Delta-\kappa r^2)u-\gamma v+D,
\qquad \kappa=\frac{3\alpha}{8\omega_0},\quad D=\frac{F}{2\omega_0}
\end{aligned}
\label{eq-dyn-oscillator-average}
$$

となる。$\kappa=0$ では線形で、平衡振幅は $r=|D|/\sqrt{\gamma^2+\Delta^2}$ である。元の線形方程式の定常振幅は

$$
\frac{\eps|F|}{\sqrt{(\omega_0^2-\omega_d^2)^2+(2\eps\gamma\omega_d)^2}}
$$

であり、$\omega_d=\omega_0+\eps\Delta$ を代入して展開すると同じ値に $O(\eps)$ で近づく。外力の振幅は $\eps|F|$ でも、共鳴付近では分母も $\eps$ のオーダーになるため、応答は一般に小さいとは限らない。

$\kappa\ne0$ では、平衡条件の両式から $r^2(\gamma^2+(\Delta-\kappa r^2)^2)=D^2$ を得る。実際、$a=\Delta-\kappa r^2$ と置き、$-\gamma u-av=0$、$au-\gamma v=-D$ を解けば、$u=-aD/(\gamma^2+a^2)$、$v=\gamma D/(\gamma^2+a^2)$ である。従って非負の根 $s=r^2$ ごとに平衡が一つである。$D=0$ なら $s=0$ だけが該当する。右辺を定数として左辺を $G(s)$ と書くと $G'(s)=\gamma^2+\Delta^2-4\Delta\kappa s+3\kappa^2s^2$ である。$\Delta\kappa>0$ かつ $\Delta^2>3\gamma^2$ のとき正の極大・極小が生じ、$D^2$ の値によっては三根をもつ。この応答の多価性は、非線形で周波数が振幅に依存することから出る。

:::example[線形共鳴の幅]
$\omega_0=2,\gamma=1/2,F=2$ なら $D=1/2$ である。$\Delta=0$ で $r=1$、$\Delta=1/2$ で $r=1/\sqrt2$ となる。実際の駆動周波数の差は $\eps/2$ なので、この応答の幅は $\eps$ とともに狭くなる。
:::

:::example[非線形応答に現れる二つの折り畳み]
$\gamma=1,\Delta=2,\kappa=1$ なら $G(s)=s[1+(2-s)^2]$ で、$G'(s)=3s^2-8s+5$ の根は $1,5/3$ である。対応する $G$ の値は $2,50/27$ であり、その間の $D^2$ を選ぶと三つの振幅が得られる。これは一次平均の平衡点の数である。元の周期解へ移すには、各平衡点が非退化であることを調べる。
:::

### 多価の振幅応答を正確な根で調べる

非線形の周波数ずれを含む振幅条件は三次式である。したがって加えた力の大きさが一つでも、定常状態が複数に分かれる場合がある。以下では同じ減衰のもとで、周波数差や力を変えたときの正の根を取り、元の二成分を復元してから安定性を判定する。

:::example[同じ外力に対する三つの振幅]
平均系の平衡点で $s=u^2+v^2$ と置き、振幅条件と力の大きさの差を因数分解する。外力は $v$ 成分の定数項として加える。

$$
\begin{aligned}
\gamma&=\sqrt2\\
\Delta&=3\\
\kappa&=1\\
D&=\sqrt{6}\\
G(s)-D^2&=s[2+(3-s)^2]-6=(s-1)(s-2)(s-3)
\end{aligned}
$$

正の実根だけ残し、対応する位置 $(u,v)$ を取り直す。$\beta=\Delta-\kappa s$ とおくと二元一次方程式の逆行列により次の式となる。

$$
u=-\frac{\beta D}{\gamma^2+\beta^2},\qquad v=\frac{\gamma D}{\gamma^2+\beta^2}
$$

得られた二成分から二乗和を計算すると $s$ に一致する。安定性は、跡が $-2\gamma$ で一定であることと、行列式 $G^\prime(s)$ の符号から決まる。

| $s$ | $u$ | $v$ | 行列式 | 平均系での型 |
| --- | --- | --- | --- | --- |
| $1$ | $-D/3$ | $\sqrt2D/6$ | $2$ | 吸引 |
| $2$ | $-D/3$ | $\sqrt2D/3$ | $-1$ | 鞍点 |
| $3$ | $0$ | $D/\sqrt2$ | $2$ | 吸引 |

中央の根でだけ行列式が負である。従って同じ外力に対して、振幅の小さい応答と大きい応答はともに吸引的で、その間に鞍型の応答がある。三根はいずれも単根なので、後述する周期解の存在判定により、元の方程式でも対応する周期解が持続する。
:::

:::example[駆動周波数を固有周波数より下げたときの単一応答]
平均系の平衡点で $s=u^2+v^2$ と置き、振幅条件と力の大きさの差を因数分解する。外力は $v$ 成分の定数項として加える。

$$
\begin{aligned}
\gamma&=1\\
\Delta&=-1\\
\kappa&=1\\
D&=\sqrt{5}\\
G(s)-D^2&=\left(s - 1\right) \left(s^{2} + 3 s + 5\right)
\end{aligned}
$$

正の実根だけ残し、対応する位置 $(u,v)$ を取り直す。$\beta=\Delta-\kappa s$ とおくと二元一次方程式の逆行列により次の式となる。

$$
u=-\frac{\beta D}{\gamma^2+\beta^2},\qquad v=\frac{\gamma D}{\gamma^2+\beta^2}
$$

得られた二成分から二乗和を計算すると $s$ に一致する。安定性は、跡が $-2\gamma$ で一定であることと、行列式 $G^\prime(s)$ の符号から決まる。

| $s$ | $u$ | $v$ | 行列式 | 平均系での型 |
| --- | --- | --- | --- | --- |
| $1$ | $\frac{2 \sqrt{5}}{5}$ | $\frac{\sqrt{5}}{5}$ | $9$ | 吸引 |

この場合は $G^\prime(s)=3s^2+4s+2$ で、正の $s$ で厳密に正である。従ってどの正の外力でも振幅は一意である。根の一意性を、三次式に現れる複素根の値まで調べずに示せる。
:::

## 平均化定理

平均を引いた関数は各周期で積分が零になる。従って時刻で積分しても原始関数は周期的で有界である。この原始関数を座標変換に使うと、速い周期変動を二次の剰余に移せる。

:::theorem[一次平均化定理]{#thm-dyn-first-averaging}
式 :eqref[eq-dyn-averaging-system] で $0<\eps\ll1$ とする。$f,g$ は時刻で連続、状態について $C^1$ とし、必要な近傍で関数と状態微分は一様に有界とする。$0\le t\le L/\eps$ で平均解が定義域の内側のコンパクト集合にとどまるとする。同じ初期値からの元の解はその期間に存在し、$\sup_{0\le t\le L/\eps}\norm{x(t)-\bar x(t)}=O(\eps)$ となる。定数は $L$ とコンパクト集合の近傍の評価だけで選べる。
:::

:::proof
$u(t,y)=\int_0^t(f(s,y)-\bar f(y))\dd s$ と置く。周期積分が零なので $u$ と $D_yu$ は時間一様に有界である。$x=y+\eps u(t,y)$ は $D_yx=I+\eps D_yu$ が一様に恒等写像に近く、必要な近傍で局所的に反転できる。関係式を微分すると

$$
(I+\eps D_yu)\dot y=\eps f(t,y+\eps u)-\eps(f(t,y)-\bar f(y))+\eps^2 g(t,y+\eps u,\eps).
$$

状態に関する Lipschitz 評価と $I+\eps D_yu$ の逆の評価から $\dot y=\eps\bar f(y)+\eps^2 R(t,y,\eps)$、$\norm R\le C$ を得る。$u(0,y)=0$ なので $y(0)=\bar x(0)$ である。$K$ を平均場の Lipschitz 定数として

$$
\norm{y(t)-\bar x(t)}
\le \eps K\int_0^t\norm{y(s)-\bar x(s)}\dd s+C\eps^2 t
$$

となる。Gronwall の不等式により右辺から $C\eps^2 t e^{\eps Kt}$ の上界を得る。$t\le L/\eps$ では $CLe^{KL}\eps$ 以下である。平均解のコンパクト集合から定義域の外まで正の距離をとれば、$\eps$ を小さくすることで退出時刻が起こらない。$x-y=\eps u$ の上界を加えればよい。
:::

平均場に双曲型平衡点があると、周期摂動系では近くに周期軌道が現れる。実際に一周期写像について定量的に計算すると、この主張を確認できる。

:::proposition[平均場の平衡点と周期解]{#prop-dyn-averaged-periodic}
$f,g$ が $C^2$ で時刻について $T$ 周期とする。$\bar f(p)=0$、$D\bar f(p)$ が可逆なら、小さな非零の $\eps$ に対し $p$ 近くに $T$ 周期解が一つある。$D\bar f(p)$ の全固有値の実部が非零なら、十分小さい正の $\eps$ でこの周期解は双曲型である。安定方向の次元は実部が負の固有値の代数的重複度の和に等しい。
:::

:::proof
初期値 $x$ から出る解を $X_\eps(t,x)$ とし、その一周期写像を $P_\eps$ とする。$\eps\ne0$ なら

$$
Q(\eps,x)=\frac{P_\eps(x)-x}{\eps T}
=\frac1T\int_0^T
\bigl[f(t,X_\eps(t,x))+\eps g(t,X_\eps(t,x),\eps)\bigr]\dd t
$$

である。右辺は有限時間の解の滑らかな依存性により $\eps=0$ でも $C^1$ であり、$X_0(t,x)=x$ なので $Q(0,x)=\bar f(x)$ である。陰関数定理により $Q(\eps,p_\eps)=0$ を満たす $p_\eps$ が $p$ 近くに一意に存在する。これは $\eps\ne0$ で $P_\eps$ の不動点である。時刻周期性と解の一意性から、対応する解は $T$ 周期である。

$DP_\eps(p_\eps)=I+\eps T D_xQ(\eps,p_\eps)$ であり、$D_xQ(\eps,p_\eps)\to D\bar f(p)$ である。従って $D_xQ$ の固有値 $\lambda$ の実部の符号とその重複度は、極限行列の固有値が虚軸を避ければ変わらない。しかも $|1+\eps T\lambda|^2-1=2\eps T\Re\lambda+\eps^2T^2|\lambda|^2$ なので、$\eps>0$ が十分小さいとき、これは $\Re\lambda$ と同じ符号である。従って :ref[thm-dyn-stable-manifold] により周期解の安定方向の次元が得られる。全固有値の実部が負の場合は局所吸引的である。
:::

:::example[線形の周期解を直接求める]
$\dot x=\eps(-2x+\cos t)$ の周期解を $x=a\cos t+b\sin t$ と置くと $b=\eps(1-2a)$、$a=2\eps b$ である。よって $a=2\eps^2/(1+4\eps^2)$、$b=\eps/(1+4\eps^2)$ である。平均場の零点 $0$ からの差は $O(\eps)$、一周期の変分倍率は $e^{-4\pi\eps}$ である。
:::

:::example[平均零なら二次が必要になる]
$\bar f=0$ とし、$g=0$ とする。上の変換で $\partial_tu=f$ だから、$C^2$ の仮定の下では $\dot y=\eps^2 D_y f(t,y)u(t,y)+O(\eps^3)$ である。例えば二変数で $f=(\cos t,x\sin t)$ とすると $u=(\sin t,x(1-\cos t))$、$D_y f\cdot u=(0,\sin^2 t)$ なので二次の平均は $(0,1/2)$ である。実際 $x(t)=x_0+\eps\sin t$ を第二式に代入すると、一周期ごとに $y$ が $\pi\eps^2$ 増える。一次平均が零でも、より長い時間に変位が蓄積する。
:::

### 積分だけで得る一次平均の誤差

:ref[thm-dyn-first-averaging] の近恒等変換を作らず、振動成分を積分して評価することもできる。$u_t=f-\bar f$ の原始関数を使い、実際の軌道に沿って連鎖律を適用すると

$$
\int_0^t(f(s,x(s))-\bar f(x(s)))\dd s
=u(t,x(t))-u(0,x(0))-\int_0^t D_xu(s,x(s))\dot x(s)\dd s.
$$

$u,D_xu$ は時間一様に有界であり、評価領域で $\dot x=O(\eps)$ である。従って平均を引いた積分の大きさは $C(1+\eps t)$ 以下である。元の積分方程式を平均解と比較すれば

$$
\norm{x(t)-\bar x(t)}\le C(\eps+\eps^2t)
  +\eps L_f\int_0^t\norm{x(s)-\bar x(s)}\dd s
$$

となる。最初に固定幅の管を出る時刻まで適用し、Gronwall の上界が $O(\eps)$ で管幅より小さいことから退出を排除できる。この証明では座標変換の反転が不要であり、平均化の誤差が時刻積分の端の二点と、状態変動に由来する一つの積分から来ることが分かる。

### 吸引平衡点の近くで全時間の誤差を抑える

有限の遅い時間での評価では、Gronwall の指数因子により上界が時間とともに大きくなる。平均場が局所的に軌道間の距離を縮める場合は、摂動で生じた差も減衰するので全時間の評価ができる。

:::proposition[吸引型の平均平衡点付近での一様近似]{#prop-dyn-uniform-averaging-sink}
:ref[thm-dyn-first-averaging] の局所的な滑らかさと有界性を仮定する。$\bar f(p)=0$ で、$A=D\bar f(p)$ の全固有値の実部が負であるとする。このとき $p$ の十分小さい固定近傍から共通の初期値で出る元の解と平均解の差は、十分小さい正の $\eps$ について全 $t\ge0$ で一様に $O(\eps)$ である。
:::

:::proof
$P=\int_0^\infty e^{tA^\top}e^{tA}\dd t$ と置く。$A$ の指数減衰により積分は収束し、$v\ne0$ に対する被積分関数 $\|e^{tA}v\|^2$ は $t=0$ で正なので $P$ は正定値である。被積分関数を微分して積分すると $A^\top P+PA=-I$ となる。よって $p$ の近くの小さい凸な楕円体では、$D\bar f$ の連続性により、ある $\mu>0$ について

$$
\langle y-z,P(\bar f(y)-\bar f(z))\rangle_{\R^d}
\le-\mu\|y-z\|_P^2,\qquad
\|v\|_P=\sqrt{v^\top Pv}
$$

が成り立つ。実際、微分行列の対称部の評価を $z$ と $y$ の線分上で積分すればよい。特に $z=p$ として、平均解が初めに内側の小さい楕円体にいれば、そこから出ずに $p$ に近づく。

近恒等変換後の軌道 $y$ は $\dot y=\eps\bar f(y)+\eps^2R(t,y,\eps)$、$\|R\|_P\le C$ を満たす。差 $d=y-\bar x$ のノルムの上右微分は $-\eps\mu\|d\|_P+\eps^2 C$ 以下である。$d\ne0$ ではノルムの微分公式、$d=0$ では右辺のベクトルのノルムの評価による。従って

$$
\|y(t)-\bar x(t)\|_P
\le e^{-\eps\mu t}\|y(0)-\bar x(0)\|_P
   +\frac{\eps C}{\mu}(1-e^{-\eps\mu t})
$$

である。初期差は $O(\eps)$ であり、$\eps$ を十分小さくすればこの上界は内側の楕円体から評価領域の境界までの距離より小さい。最初の退出時刻まで評価して退出を排除できる。元の座標へ戻しても近恒等変換で増える差は一様に $O(\eps)$ である。
:::

### 二次平均の一般形

$f,g$ が状態について十分滑らかで、時刻 $T$ 周期とする。$u_t=f-\bar f$ かつ $u(0,y)=0$ となる $u$ を使う。$x=y+\eps u(t,y)$ の式を二次まで展開すると

$$
\dot y=\eps\bar f(y)+\eps^2 R_2(t,y)+O(\eps^3),\qquad
R_2=Df(t,y)u+g(t,y,0)-D_yu\,\bar f.
$$

右辺の最初の項は $f(t,y+\eps u)$ の Taylor の一次係数から、最後の項は $(I+\eps D_yu)^{-1}$ の一次係数から来る。切断せず同じコンパクト近傍で剰余の上界を取るには、$g$ の $\eps$ 微分にも上界を置く。次に $v_t=R_2-\bar R_2$ を満たす周期的な原始関数を選び、$y=w+\eps^2 v(t,w)$ に変える。この変換の状態微分と $\dot w$ の積は三次であるため、二次までで

$$
\dot w=\eps\bar f(w)+\eps^2\bar R_2(w)+O(\eps^3)
$$

となる。$\bar f=0$ ならこの式が $1/\eps^2$ の時間尺度での平均方程式になる。実際、剰余は三次、二次平均の Lipschitz 定数は $\eps^2$ のオーダーなので、Gronwall の評価でその尺度まで $O(\eps)$ の誤差に抑えられる。一方、$\bar f$ が非零なら $\eps$ の線形増幅が残り、このまま $1/\eps^2$ まで評価を伸ばせるとは限らない。

:::example[符号の異なる順序の効果]
$\dot x=\eps\cos t$、$\dot y=\eps x\sin t$ では $y$ の一周期変位は $\pi\eps^2$ である。二つの外力の位相を交換して $\dot x=\eps\sin t$、$\dot y=\eps x\cos t$ とすれば、一周期変位は $-\pi\eps^2$ になる。いずれも一次は零である。第二式の外力が作用するときの $x$ は第一式の過去の積分を記憶しているので、平均する順序は二次で効く。
:::

### 振幅応答から周期解と安定性を読む

式 :eqref[eq-dyn-oscillator-average] の右辺を $(-\gamma u-av,au-\gamma v+D)$、$a=\Delta-\kappa s$ と書く。その Jacobi 行列は

$$
\begin{pmatrix}
-\gamma+2\kappa uv&-a+2\kappa v^2\\
a-2\kappa u^2&-\gamma-2\kappa uv
\end{pmatrix}.
$$

跡は $-2\gamma$、行列式は $\gamma^2+a^2-2\kappa as=G'(s)$ である。従って $G'(s)>0$ の枝は線形吸引型、$G'(s)<0$ の枝は鞍点である。$G'(s)\ne0$ なら :ref[prop-dyn-averaged-periodic] の陰関数の議論を、角を独立変数にした $2\pi$ 周期系へ適用できる。吸引枝は元の強制振動子にも吸引周期解を与え、鞍点枝は変分乗数の一方が単位円の外になる。折り畳みに近づいて $G'(s)$ が零に近いと、一次平均の誤差で実際の分岐位置も移りうる。一時刻の符号に加え、非退化性を確かめるのはこのためである。

### 原始関数の積分定数と二次平均

二次項を $W=R_2$ と書く。原始関数を $u$ から $u+b(y)$ に取り替えると、$W$ の時刻平均は $D\bar f\,b-Db\,\bar f$ だけ変わる。実際、二次項の差は $D_yf\,b-Db\,\bar f$ であり、$b$ は時刻に依らないので平均の外に出せる。この差は一次平均の座標を $y+\eps b(y)$ に取り替えたときに生じる項である。特に $\bar f=0$ なら二次平均そのものが積分定数に依らない。以下では $u$ の時刻平均を零に選ぶ。

この選択では $u(0,y)$ は一般に零でない。そこで二次平均の解を元の初期値から出すなら、変換した座標との初期差 $e(0)=O(\eps)$ を含める。$\bar f=0$ のとき、変換後の軌道と二次平均の差 $e(t)$ は、両方が評価領域内にある間、

$$
e(t)\le \bigl(e(0)+C\eps^3t\bigr)e^{L\eps^2t}
$$

を満たす。任意の固定した $S>0$ に対し、二次平均解が定義域内部のコンパクト集合に留まるなら、$t\le S/\eps^2$ で右辺は $O(\eps)$ である。最初の退出時刻までこの不等式を使えば、一次平均と同じように退出を排除できる。

:::example[零平均の線形時変係数]
時刻平均をとり、時刻で微分した値が振動成分になる周期関数 $u$ を作る。積分定数は $u$ の平均が零になるように選ぶ。

$$
\begin{aligned}
f&=x \sin{\left(t \right)} + \cos{\left(t \right)}\\
g&=0\\
\bar f&=0\\
u&=- x \cos{\left(t \right)} + \sin{\left(t \right)}
\end{aligned}
$$

座標変換 $x=y+\eps u(t,y)$ を微分してから、$f(t,x)$ を $y$ で Taylor 展開する。二次の残りは、$f$ の展開項、座標変換の微分による項、および元の二次項の和になる。

$$
\begin{aligned}
W&=f_xu-u_x\bar f+g\\
W&=- \frac{x \sin{\left(2 t \right)}}{2} + \sin^{2}{\left(t \right)}\\
\bar W&=\frac{1}{2}
\end{aligned}
$$

ここで平均を取るのは $f_x$ と $u$ のそれぞれでなく、その積である。積の平均が零でない場合には、一次の平均が消えても長時間で変位が見える。

二次平均は定数 $1/2$ である。従って二次の遅い時刻 $\tau=\eps^2 t$ では、平均解は初期値に $\tau/2$ を加えた値である。例えば $\tau=2$ では速度の一次平均が零にもかかわらず、変位は $1+O(\eps)$ になる。
:::

:::example[振動による二次増分と散逸]
時刻平均をとり、時刻で微分した値が振動成分になる周期関数 $u$ を作る。積分定数は $u$ の平均が零になるように選ぶ。

$$
\begin{aligned}
f&=x^{2} \sin{\left(2 t \right)} + x \cos{\left(2 t \right)}\\
g&=- x\\
\bar f&=0\\
u&=\frac{x \left(- x \cos{\left(2 t \right)} + \sin{\left(2 t \right)}\right)}{2}
\end{aligned}
$$

座標変換 $x=y+\eps u(t,y)$ を微分してから、$f(t,x)$ を $y$ で Taylor 展開する。二次の残りは、$f$ の展開項、座標変換の微分による項、および元の二次項の和になる。

$$
\begin{aligned}
W&=f_xu-u_x\bar f+g\\
W&=- \frac{x^{3} \sin{\left(4 t \right)}}{2} + x^{2} \sin^{2}{\left(2 t \right)} - \frac{x^{2} \cos^{2}{\left(2 t \right)}}{2} + \frac{x \sin{\left(4 t \right)}}{4} - x\\
\bar W&=\frac{x \left(x - 4\right)}{4}
\end{aligned}
$$

ここで平均を取るのは $f_x$ と $u$ のそれぞれでなく、その積である。積の平均が零でない場合には、一次の平均が消えても長時間で変位が見える。

二次平均は $x^2/4-x$ である。$0<x<4$ では平均速度が負で、$x>4$ では正になる。これは $\tau$ 時刻での方程式の安定性である。元の時変系で原点が吸引的であることは線形方程式 $\dot x=(\eps\cos(2t)-\eps^2)x$ の一周期の倍率からも確認できる。
:::

## 摂動されたハミルトン系

保存系の分離曲線は、同じエネルギー値をもつ安定枝と不安定枝が一致してできる。摂動があるとエネルギーは保存しない。そこで二つの枝を同じ断面まで追い、エネルギーの差を測る。

:::definition[周期摂動付き平面 Hamilton 系]{#def-dyn-perturbed-hamiltonian}
標準行列 $J=\left(\begin{smallmatrix}0&1\\-1&0\end{smallmatrix}\right)$ を取り、

$$
\dot z=J\grad H(z)+\eps G(z,t),\qquad G(z,t+T)=G(z,t)
\label{eq-dyn-perturbed-hamiltonian}
$$

とする。非摂動系の双曲型鞍点 $z_*$ へ両方向で近づくループを $\gamma(t)$ とする。$H$ の軌道微分は $\eps\grad H(z)\cdot G(z,t)$ である。
:::

$q=\gamma(0)$ において $\grad H(q)\ne0$ である。$\grad H(q)$ の向きの小線分を横断面として使うと、その上で $H$ は狭義単調で局所座標になる。時間位相 $\theta$ を固定して、この線分で不安定枝の $H$ 値から安定枝の値を引けば分裂量が得られる。

:::example[散逸と強制の分け方]
$H(q,p)=p^2/2+q^2/2$、$G=(0,-\delta p+F\cos\omega t)$ なら $\dot H=\eps(-\delta p^2+Fp\cos\omega t)$ である。散逸は常に非正、強制は位相によって正負が変わる。$q=\cos t$、$p=-\sin t$ の一周期では $\omega=1$ の余弦強制の仕事が零で、散逸の仕事は $-\eps\delta\pi$ となる。
:::

:::example[Hamilton 型摂動は常に分裂するわけではない]
$G=J\grad K(z)$ が時刻に依らないなら全系は $H+\eps K$ を保存する。非摂動ループに沿う $\grad H\cdot J\grad K=-\dv{}{t}K(\gamma(t))$ を積分すると、両端で同じ鞍点 $z_*$ に近づくので零になる。保存量の変更だけで起こる一次のエネルギー差は零である。
:::

### 時間依存のポテンシャルと保存量の区別

$K(q,p,t)$ を摂動 Hamiltonian とすれば、全方程式は $X_{H+\eps K}$ である。発散は混合偏微分の差 $\partial_q\partial_p(H+\eps K)-\partial_p\partial_q(H+\eps K)$ なので零である。従って固定時刻の写像は Liouville の行列式公式で面積を保存する。一方、$H+\eps K$ の全微分では状態の微分との内積が相殺しても明示的な時刻微分が残る。これにより面積保存とエネルギー保存を区別できる。

:::example[剛性が周期変化する Duffing 系]
エネルギーと摂動の関数を選び、標準座標の方程式により速度を組み立てる。非摂動の速度と摂動の速度は別々に求める。

$$
\begin{aligned}
H&=\frac{p^{2}}{2} + \frac{q^{4}}{4} - \frac{q^{2}}{2}\\
K&=\frac{q^{2} \cos{\left(2 t \right)}}{2}\\
X_H&=\left[\begin{matrix}p\\- q^{3} + q\end{matrix}\right]\\
X_K&=\left[\begin{matrix}0\\- q \cos{\left(2 t \right)}\end{matrix}\right]
\end{aligned}
$$

時間変化の一次係数は勾配と摂動速度の内積である。もとの $X_H$ との内積は反対称な行列のため零である。

$$
\begin{aligned}
\nabla H\cdot X_K&=- p q \cos{\left(2 t \right)}\\
\frac{\dd}{\dd t}(H+\eps K)&=\eps\partial_t K
\end{aligned}
$$

最後の式は摂動系の実際の軌道に沿う全微分である。時刻依存の Hamiltonian でも写像は面積を保存するが、関数値は外部の時刻変化によって動く。

この摂動は原点で零になり、鞍点は留まる。ただし $q=0$ 以外では仕事率は一般に零でない。$K$ に $p$ が現れないので第一式の速度は変わらず、第二式の力だけに $-\eps q\cos(2t)$ が加わる。
:::

:::example[位置と運動量の結合を変える]
エネルギーと摂動の関数を選び、標準座標の方程式により速度を組み立てる。非摂動の速度と摂動の速度は別々に求める。

$$
\begin{aligned}
H&=\frac{p^{2}}{2} + \frac{q^{4}}{4} - \frac{q^{2}}{2}\\
K&=p q \sin{\left(t \right)}\\
X_H&=\left[\begin{matrix}p\\- q^{3} + q\end{matrix}\right]\\
X_K&=\left[\begin{matrix}q \sin{\left(t \right)}\\- p \sin{\left(t \right)}\end{matrix}\right]
\end{aligned}
$$

時間変化の一次係数は勾配と摂動速度の内積である。もとの $X_H$ との内積は反対称な行列のため零である。

$$
\begin{aligned}
\nabla H\cdot X_K&=- \left(p^{2} - q^{4} + q^{2}\right) \sin{\left(t \right)}\\
\frac{\dd}{\dd t}(H+\eps K)&=\eps\partial_t K
\end{aligned}
$$

最後の式は摂動系の実際の軌道に沿う全微分である。時刻依存の Hamiltonian でも写像は面積を保存するが、関数値は外部の時刻変化によって動く。

位置に掛かる速度にも運動量に掛かる力にも摂動が入るが、互いの発散は相殺する。仕事率には $p^2$ と $q^2,q^4$ の項が同時に入る。ループに沿って積分するときは、計算を短くするため $\partial_t K$ と端点の項へ書き直せる。
:::

:::example[時刻に依らない運動エネルギー修正]
エネルギーと摂動の関数を選び、標準座標の方程式により速度を組み立てる。非摂動の速度と摂動の速度は別々に求める。

$$
\begin{aligned}
H&=\frac{p^{2}}{2} - \cos{\left(q \right)} + 1\\
K&=\frac{p^{2}}{2}\\
X_H&=\left[\begin{matrix}p\\- \sin{\left(q \right)}\end{matrix}\right]\\
X_K&=\left[\begin{matrix}p\\0\end{matrix}\right]
\end{aligned}
$$

時間変化の一次係数は勾配と摂動速度の内積である。もとの $X_H$ との内積は反対称な行列のため零である。

$$
\begin{aligned}
\nabla H\cdot X_K&=p \sin{\left(q \right)}\\
\frac{\dd}{\dd t}(H+\eps K)&=\eps\partial_t K
\end{aligned}
$$

最後の式は摂動系の実際の軌道に沿う全微分である。時刻依存の Hamiltonian でも写像は面積を保存するが、関数値は外部の時刻変化によって動く。

今度は全エネルギー $H+\eps K$ が一定である。元の $H$ の値は一般に変化するが、ループの両端で $K$ が等しければ、その変化率の一次積分は零になる。別の Hamiltonian を採用しただけの場合、時間依存の仕事を仮定しない。
:::

## Melnikov 関数

枝のエネルギー差の一次係数は、摂動による仕事率を非摂動ループ全体で積分したものになる。鞍点が周期的に動く場合は、移動した周期解の仕事率を引いてから積分する。

:::theorem[Melnikov の分裂公式]{#thm-dyn-melnikov}
式 :eqref[eq-dyn-perturbed-hamiltonian] の $H$ は $C^4$、$G$ は時刻と状態について $C^3$ とする。非摂動系に双曲型鞍点 $z_*$ へのホモクリニック軌道 $\gamma$ があるとする。摂動後も鞍点は時刻周期解として持続する。時刻位相 $\theta$ で前節の断面を通る両枝のエネルギー差は

$$
\Delta_\eps(\theta)=\eps M(\theta)+O(\eps^2),\qquad
M(\theta)=\int_{-\infty}^{\infty}\grad H(\gamma(t))\cdot G(\gamma(t),t+\theta)\dd t
\label{eq-dyn-melnikov}
$$

である。$M$ を :term[Melnikov 関数]{en="Melnikov function" yomi="めるにこふかんすう"} と呼ぶ。この展開の剰余は $\theta$ について $C^1$ で一様である。
:::

:::proof
鞍点の近傍には一周期写像の不動点の持続性（:ref[prop-dyn-hyperbolic-continuation]）により、双曲型周期解 $z_{*,\eps}(t)$ がある。$G(z_*,t)=0$ ならこれは定数解 $z_*$ である。時刻位相 $\theta$ の一周期写像を $P_{\eps,\theta}$ とし、$z_{*,\eps}(\theta)$ における局所安定・不安定多様体を取る。非摂動ループは十分遠い正負の整数周期でその局所多様体に入る。従って各局所曲線を有限回の $P_{\eps,\theta}$ またはその逆写像で断面の近傍まで延ばすと、$\eps=0$ で $\gamma$ と同じ接線をもつ曲線が得られる。固定した位相で、曲線上のパラメータを陰関数定理で解いて横断面との交点を求める。この点から延びる半軌道を $z^u_\eps(t;\theta)$（$t\le0$）、$z^s_\eps(t;\theta)$（$t\ge0$）と書く。

尾部の評価には :ref[thm-dyn-stable-manifold] の固定点方程式を使う。周期解の断面点を原点に平行移動し、$\eps=0$ の安定・不安定線形部を固定する。$\eps,\theta$ に依存する線形部の差を非線形剰余と合わせて扱う。この剰余は原点で常に零で、座標近傍を縮めれば状態に関する Lipschitz 定数を一様に小さくできる。従って重み付き列空間の固定点方程式は一様な縮小となる。剰余をパラメータだけで微分した項も原点で零なので、軌道値の減衰と同じ重みで有界である。混合微分で現れる項は、軌道の一次微分を掛けた有界な係数と、固定点方程式の線形化からなる。線形化の逆は縮小率の幾何級数で抑えられるため、周期点から引いた尾部とその $\eps,\theta,\eps\theta$ 微分はいずれも一様に指数減衰する。整数周期の間では有限時間の流れの滑らかな依存性を使えば、同じ評価を定数倍の範囲で保てる。

$A(z,\phi)=\grad H(z)\cdot G(z,\phi)$ とおく。半軌道に沿って $H(z^\ast_\eps(t;\theta))-H(z_{*,\eps}(t+\theta))$ を微分すると、$\eps[A(z^\ast_\eps(t;\theta),t+\theta)-A(z_{*,\eps}(t+\theta),t+\theta)]$ となる。このエネルギー差は各半直線の無限遠で零へ行く。過去の積分は断面での不安定枝の値から周期解の値を引いたものを、未来の積分はその周期解の値から安定枝の値を引いたものを与える。二つを足すと周期解での値が相殺し、$\Delta_\eps/\eps$ の積分表示を得る。

被積分関数とその $\eps,\theta,\eps\theta$ 微分は、尾部の評価と $A$ の滑らかさにより時間について可積分な同じ指数関数で抑えられる。そこで $\eps=0$ を代入すると、$\grad H(z_*)=0$ なので基準の項が消え、式 :eqref[eq-dyn-melnikov] の $M$ になる。積分の下で $\eps$ 微分と混合微分を取って一様に評価し、平均値の定理を使えば $\|\Delta_\eps/\eps-M\|_{C^1}=O(|\eps|)$ を得る。
:::

:::example[位相に依らない摂動]
$G=(0,-\delta p)$ とし、ループに沿う第一式が $\dot q=p$ なら $M=-\delta\int_\R p(t)^2\dd t$ である。$\delta>0$ なら厳密に負で、同じループの近くで一次の分裂はどの位相でも消えない。
:::

:::example[偶奇性による簡約]
ループの $q(t)$ が偶関数、$p(t)=q'(t)$ が奇関数とする。$G=(0,F\cos(\omega(t+\theta)))$ なら $p(t)\cos\omega t$ は奇関数、$p(t)\sin\omega t$ は偶関数である。従って $M(\theta)=-F\sin(\omega\theta)\int_\R p(t)\sin\omega t\dd t$ となる。位相を $\pi/\omega$ だけ進めると符号が逆になる。位相をどこから数えても、正弦項の振幅は同じである。
:::

### 周期 Hamiltonian 摂動と Melnikov ポテンシャル

$G(z,t)=J\grad_z K(z,t)$ とする。$K$ が滑らかで時刻について $T$ 周期なら、:term[Melnikov ポテンシャル]{en="Melnikov potential" yomi="めるにこふぽてんしゃる"} を

$$
L(\theta)=\int_\R[K(\gamma(t),t+\theta)-K(z_*,t+\theta)]\dd t
$$

と定められる。$\gamma(t)-z_*$ は両端で指数的に減衰するため、$K$ の差とその時刻偏微分の差はどちらも絶対可積分である。非摂動方程式 $\dot\gamma=J\grad H(\gamma)$ と $J^\top=-J$ により

$$
\begin{aligned}
\frac{\dd}{\dd t}\bigl[K(\gamma(t),t+\theta)-K(z_*,t+\theta)\bigr]
&=-\grad H(\gamma(t))\cdot J\grad_z K(\gamma(t),t+\theta)\\
&\quad+\partial_tK(\gamma(t),t+\theta)-\partial_tK(z_*,t+\theta)
\end{aligned}
$$

となる。両端で括弧内は零に向かうので、積分して $M=L'$ を得る。従って $M$ の周期平均は零である。$L$ に非退化な極値があれば $M$ はそこに単純零点をもつ。一方、$L$ が一定でなければ $M$ は符号を変えるが、極値で $L''$ が零になる可能性があるため、横断性には別途 $L''$ の確認が要る。

:::example[周期外力のポテンシャルから仕事を求める]
$K(q,p,t)=-Fq\cos(\omega t)$ では $J\grad_z K=(0,F\cos\omega t)$ である。正側の Duffing ループ $q_h(t)=\sqrt2\operatorname{sech}t$ と鞍点 $z_*=(0,0)$ に対し、$q_h$ は偶関数なので

$$
L(\theta)=-F\cos(\omega\theta)\int_\R q_h(t)\cos(\omega t)\dd t.
$$

従って $M(\theta)=F\omega\sin(\omega\theta)\int_\R q_h(t)\cos(\omega t)\dd t$ となる。この積分は後の Duffing 振動子の節で求める。定数外力では $\omega=0$ なので $L$ は一定であり、分裂の一次項は零である。
:::

### 非保存の平面系での重み付き分裂

Hamilton 系でエネルギーを差の座標に使えた理由は、非摂動場に沿って保存されるからである。一般の平面の非摂動場 $F$ にはその保存量がなくても、軌道方向のベクトルと一次変分が張る面積を追えば分裂の一次係数を求められる。

:::proposition[平面の重み付き Melnikov 積分]{#prop-dyn-weighted-melnikov}
平面の $C^3$ 場 $F$ が双曲型鞍点 $p$ へのホモクリニック軌道 $\gamma$ をもち、$G$ は滑らかで時刻周期的な摂動とする。$q=\gamma(0)$ に横断面をとり、不安定枝と安定枝の交点の $\eps$ 微分を $w^u(0),w^s(0)$ とする。二次元行列式を $u\wedge v$ と書くと

$$
F(q)\wedge(w^u(0)-w^s(0))
=\int_\R e^{-\int_0^t\divg F(\gamma(r))\dd r}\,
  F(\gamma(t))\wedge G(\gamma(t),t+\theta)\dd t .
\label{eq-dyn-weighted-melnikov}
$$

:::

:::proof
鞍点の周期的な続きから二つの枝を作り、その位相を $\theta$ に固定する。:ref[thm-dyn-melnikov] の証明と同じ滑らかな依存性により、各時間半直線で枝を微分できる。その変分 $w$ は $w'=DF(\gamma)w+G(\gamma,t+\theta)$ を満たす。$a(t)=F(\gamma(t))$ は同次の変分方程式を満たすため、$D(t)=a(t)\wedge w(t)$ と置くと

$$
D'=(\divg F)(\gamma(t))D+a(t)\wedge G(\gamma(t),t+\theta)
$$

である。これは任意の二次行列 $A$ について $(Aa)\wedge w+a\wedge(Aw)=(\tr A)(a\wedge w)$ となる成分計算による。両辺に跡から作った積分因子を掛け、過去半直線と未来半直線で別々に積分する。

端の項が消えることを確認する。鞍点の固有値を $\lambda>0,-\sigma<0$ とすれば、正方向では $\gamma-p$ と $a$ は安定多様体上で $e^{(-\sigma+\delta)t}$ 以下、負方向では $e^{(\lambda-\delta)t}$ 以下である。任意の小さい $\delta>0$ に対し、安定・不安定曲線の一次元方程式を積分すればこの評価を得る。発散は $\lambda-\sigma+O(\norm{\gamma-p})$ であり、剰余は可積分である。従って積分因子を掛けた $a(t)$ は正の尾で $e^{(-\lambda+\delta)t}$、負の尾で $e^{(\sigma-\delta)t}$ の定数倍以下となる。$w^s,w^u$ はそれぞれの尾で、周期軌道の有界な変分と指数減衰する変分との差として評価できる。従って端の項は零であり、積分も絶対収束する。二つの半直線の結果を合わせ、時刻零では積分因子が一であることを使えば式を得る。
:::

$F=J\grad H$ の場合は $\divg F=0$ で、$F\wedge G=\grad H\cdot G$ だから :eqref[eq-dyn-melnikov] に戻る。非保存系では発散が分裂への寄与を時間ごとに重み付けする。従って Hamilton 系の積分式をそのまま非保存系へ移すと、この因子を落としてしまう。

### 仕事の積分が進行方向を区別する

二つの領域を隔てる障壁の高さが等しければ、二つの鞍点を結ぶ保存軌道が得られる。この場合は速度の積分が零とは限らない。位置の両端が異なるため、同じエネルギー値を取ることと仕事の距離が零になることを区別する。

:::example[左の鞍点から右の鞍点へ進む]
Hamiltonian は二つの鞍点 $q=\pm1$、$p=0$ に同じ高さの障壁を持つ。エネルギーをその高さに置き、$|q|<a$ の部分で $p=\dot q$ を解く。

$$
\begin{aligned}
H&=\frac{p^{2}}{2} - \frac{q^{4}}{4} + \frac{q^{2}}{2}\\
H_{\rm sad}&=\frac{1}{4}\\
\dot q&=\frac{\sqrt{2} \left(1 - q^{2}\right)}{2}
\end{aligned}
$$

符号は接続の進む向きに対応する。固定した時刻 $t=0$ で $q=0$ となるよう積分定数を選ぶ。

解を双曲線正接で表し、その導関数を使って外力と摩擦の仕事を積分する。外力が定数なら速度の積分は両端の位置の差である。

$$
\begin{aligned}
q(t)&=\tanh{\left(\frac{\sqrt{2} t}{2} \right)}\\
p(t)&=\frac{\sqrt{2} \left(1 - \tanh^{2}{\left(\frac{\sqrt{2} t}{2} \right)}\right)}{2}\\
\int_{\R}p(t)\dd t&=2\\
\int_{\R}p(t)^2\dd t&=\frac{2 \sqrt{2}}{3}
\end{aligned}
$$

二乗の積分は $s=\tanh(t/\sqrt2)$ で $-1$ から $1$ までの多項式の積分に直せる。この値は速度を反転しても同じである。

定数の力と速度に比例する減衰を足した場合の一次の仕事を、二つの寄与に分けて足す。

$$
\begin{aligned}
F&=1\\
\delta&=\frac{1}{2}\\
M&=2 - \frac{\sqrt{2}}{3}
\end{aligned}
$$

同一点へ戻るループなら定数の力の寄与は零になるが、ここでは両端が異なるため零にならない。二つの鞍点で元の $H$ の勾配が零であるため、摂動で鞍点が動いても元の $H$ の両端の値の変化は一次で現れない。
:::

:::example[逆の向きに進む]
Hamiltonian は二つの鞍点 $q=\pm1$、$p=0$ に同じ高さの障壁を持つ。エネルギーをその高さに置き、$|q|<a$ の部分で $p=\dot q$ を解く。

$$
\begin{aligned}
H&=\frac{p^{2}}{2} - \frac{q^{4}}{4} + \frac{q^{2}}{2}\\
H_{\rm sad}&=\frac{1}{4}\\
\dot q&=\frac{\sqrt{2} \left(q^{2} - 1\right)}{2}
\end{aligned}
$$

符号は接続の進む向きに対応する。固定した時刻 $t=0$ で $q=0$ となるよう積分定数を選ぶ。

解を双曲線正接で表し、その導関数を使って外力と摩擦の仕事を積分する。外力が定数なら速度の積分は両端の位置の差である。

$$
\begin{aligned}
q(t)&=- \tanh{\left(\frac{\sqrt{2} t}{2} \right)}\\
p(t)&=- \frac{\sqrt{2} \left(1 - \tanh^{2}{\left(\frac{\sqrt{2} t}{2} \right)}\right)}{2}\\
\int_{\R}p(t)\dd t&=-2\\
\int_{\R}p(t)^2\dd t&=\frac{2 \sqrt{2}}{3}
\end{aligned}
$$

二乗の積分は $s=\tanh(t/\sqrt2)$ で $-1$ から $1$ までの多項式の積分に直せる。この値は速度を反転しても同じである。

定数の力と速度に比例する減衰を足した場合の一次の仕事を、二つの寄与に分けて足す。

$$
\begin{aligned}
F&=1\\
\delta&=\frac{1}{2}\\
M&=-2 - \frac{\sqrt{2}}{3}
\end{aligned}
$$

逆向きの計算では仕事がさらに負になり、同じ摩擦と定数の外力では一次で釣り合わない。一つの横断面で符号を比較するときには、どちらの接続を基準にしたかを指定する。
:::

### 周期的な剛性変化の積分

力を加える代わりに、ポテンシャルの曲率を周期的に変えることもできる。非摂動 Duffing ループ $q_h=\sqrt2\operatorname{sech}t$、$p_h=\dot q_h$ に対し、$K=q^2\cos(\omega t)/2$ を加える。このときエネルギー変化の一次係数は $-p_hq_h\cos(\omega(t+\theta))$ である。$p_h q_h=(q_h^2/2)'$ と書いて部分積分すると、境界項は零であり
$$
M_K(\theta)
=-\omega\sin(\omega\theta)\int_{\R}\operatorname{sech}^2t\cos(\omega t)\,\dd t
$$
を得る。正弦成分のうち $\sin(\omega t)$ に比例する積分は、$\operatorname{sech}^2t$ の偶性により零である。

:::lemma[双曲線正割の二乗の Fourier 積分]{#lem-dyn-sech-square-transform}
$\omega>0$ で
$$
\int_{\R}\frac{e^{i\omega t}}{\cosh^2t}\,\dd t=\frac{\pi\omega}{\sinh(\pi\omega/2)}
$$
である。
:::

:::proof
:ref[thm-cpx-residue-theorem]{text="留数定理"} を使い、下辺が $[-R,R]$、上辺の虚部が $\pi$ の長方形を反時計回りに回る。$\cosh(z+i\pi)=-\cosh z$ なので、上下の積分の和は実軸上の積分の $1-e^{-\pi\omega}$ 倍である。縦の辺では $|\cosh(\pm R+iy)|^2=\sinh^2R+\cos^2y$ であり、$|e^{i\omega(\pm R+iy)}|\le1$ だから、縦の積分は $R\to\infty$ で零に向かう。

長方形内の極は $z=i\pi/2$ の二位の極だけである。そこからのずれを $w$ とすれば、$\cosh(i\pi/2+w)=i\sinh w$ であり、$\cosh^{-2}(i\pi/2+w)=-w^{-2}+O(1)$ である。従って留数は $-i\omega e^{-\pi\omega/2}$ であり、周回積分は $2\pi\omega e^{-\pi\omega/2}$ となる。これを $1-e^{-\pi\omega}$ で割ると所望の式になる。
:::

:::example[位相で正負が交換される剛性変化]
$\omega=2$ とすれば $M_K(\theta)=-4\pi\sin(2\theta)/\sinh\pi$ である。例えば $\theta=\pi/4$ では最小値、$\theta=3\pi/4$ では最大値になる。正弦の零点では導関数は零でないため、単独でも一周期内に二つの横断的な零を与える。摂動が時刻に依存せず $\omega=0$ なら積分は恒等的に零であり、単純零点の議論は使えない。
:::

:::example[散逸があると必要なポテンシャル振幅が変わる]
$K=Bq^2\cos(2t)/2$、減衰力を $-\gamma p$ とすれば、$M=-4\gamma/3-(4\pi B/\sinh\pi)\sin(2\theta)$ である。$\gamma=1/2$ なら定数損失は $2/3$ であり、$|B|>\sinh\pi/(6\pi)$ の場合に単純零点が二つある。加法外力で計算した閾値の双曲線余弦をそのまま流用せず、ここでは剛性による仕事に対応する双曲線正弦で評価する。
:::

## 横断的ホモクリニック点の存在判定

横断面と位相を合わせた二次元の断面では、二つの枝は位相を横座標とする曲線になる。高さの差が零で、そこで傾きが違えば、交点は横断的である。

:::proposition[単純零点の判定]{#prop-dyn-melnikov-simple}
$M(\theta_0)=0$、$M'(\theta_0)\ne0$ なら、小さな非零 $\eps$ に対して近くに横断的なホモクリニック点がある。$M$ が全位相で零から正の距離をもつなら、この断面で追ったループの枝は小さな非零 $\eps$ に対して交わらない。
:::

:::proof
:ref[thm-dyn-melnikov] により $\Delta_\eps/\eps$ は $C^1$ で $M$ に近づく。単純零点近くで値が両側で異符号、微分が零から離れるので、近くに一意の零点があり $\Delta_\eps'\ne0$ となる。$H$ が断面の局所座標なので、二曲線の高さが等しく傾きが異なる。流れ方向を断面で除いているため、断面の横断性は拡大系での多様体の横断性に対応する。後半は一様剰余の上界と $M$ の絶対値の正の下界を比較すれば従う。
:::

交差から双曲的な軌道の族を作る議論は [Smale の馬蹄とカオス](/math/dynamical-systems/chaos/) で扱う。この判定はループ近くの特定の枝に対するものである。$M=0$ が重根なら二次以上の差が交差を消したり二つに分けたりしうる。

:::example[強制と散逸の閾値]
$M(\theta)=3\sin\theta-2$ では $\sin\theta=2/3$ の二点があり、$\cos\theta\ne0$ なので両方単純である。$M=2\sin\theta-2$ の $\theta=\pi/2$ では微分が零で、例えば剰余の定数符号を変えると零点数が零個と二個に分かれる。
:::

:::example[位相シフトの効果]
$M(\theta)=\sin\theta+\cos\theta$ は $\sqrt2\sin(\theta+\pi/4)$ なので、零点は $-\pi/4$ と $3\pi/4$ である。ループ上の断面を時間だけずらすと位相の名前が変わり、積分に入る $M$ は平行移動するが、単純零点の有無は変わらない。
:::

### 二つの調和波がある場合の零点

一つの周波数では振幅と定数項の大小だけで判定できる。複数の周波数があると、その区間内でどの極値が定数項を越えるかを調べる必要がある。Duffing の加法外力を $\sum_{k=1}^2F_k\cos(kt+\varphi_k)$ とすれば、線形性により加法外力の一次 Melnikov 関数は $\sum_{k=1}^2 A_k\sin(k\theta+\varphi_k)$ であり、$A_k=\sqrt2\pi k F_k/\cosh(\pi k/2)$ である。各 $F_k$ と $\varphi_k$ を選んで係数と位相を独立に設定できる。

:::example[減衰を下げると四つの単純零点]
位相を $\pi/2$ にそろえ、$A_1=1/2$、$A_2=1$、定数損失 $1/2$ を選ぶ。すなわち減衰係数を $3/8$ に置く。$\cos2\theta=2\cos^2\theta-1$ を使うと
$$
M(\theta)=\tfrac12\cos\theta+\cos2\theta-\tfrac12
         = 2u^2+\tfrac12u-\tfrac32,\qquad u=\cos\theta.
$$
根は $u=3/4,-1$ である。このままでは $u=-1$ が区間の端にあり、$\theta=\pi$ で $\dd u/\dd\theta=0$ のため単純零ではない。そこで定数損失を $1/4$ に減らす。今度の多項式 $P(u)=2u^2+u/2-5/4$ は $P(-1)=1/4>0$、$P(0)=-5/4<0$、$P(1)=5/4>0$ である。従って $(-1,0)$、$(0,1)$ に異なる二つの単根があり、それぞれで $\cos\theta=u$ は一周期に二つの解をもつ。$\sin\theta\ne0$ なので四つとも単純零点である。
:::

:::example[同じ外力で二零点に減る]
定数損失を $1$ にすると多項式は $2u^2+u/2-2$ になる。$P(-1)=-1/2<0$、$P(0)=-2<0$、$P(1)=1/2>0$ である。負根は $-1$ より小さく、正根だけが $(0,1)$ に入る。従って一周期で二つだけ単純零点をもつ。負根が区間から抜けたのは、減衰の増加により極値が正から負へ移動したためである。係数の絶対値の和と減衰の比較だけではこの個数は決まらない。
:::

:::proposition[単純零の個数の局所持続]{#prop-dyn-melnikov-zero-count}
$M$ が周期的で、零点が有限個あり全て単純とする。$M+R_\eps$、$\|R_\eps\|_{C^1}\to0$ であれば、十分小さい $\eps$ で各零の近くに一つずつ単純零があり、それ以外にはない。
:::

:::proof
各零に対し互いに素な短い弧を取り、その閉弧で $M'$ の符号が一定かつ絶対値が正の下限をもつようにする。端点では $M$ の符号が互いに反対になる。$C^1$ 誤差が小さければ端点での符号と弧の内部の狭義単調性は変わらないから、中間値の定理と単調性で各弧にちょうど一つの単純零がある。補集合はコンパクトで $|M|$ が正の最小値をもつため、$C^0$ 誤差をこれより小さくすると零は追加されない。
:::

:::example[重根は小さな同符号の補正で消える]
$M(\theta)=1-\cos\theta$ では $\theta=0\bmod2\pi$ が二重零である。定数の正の補正を足すと零は全くなくなる。負の補正 $-\epsilon$、$0<\epsilon<1$ を足すなら $\cos\theta=1-\epsilon$ に二解があり、それぞれは単純になる。この二つの場合はどちらも $C^1$ で元の関数へ収束する。従って重根で交差を主張するには、一次の式だけでは足りない。
:::

## 強制 Duffing 振動子の解析

時間と変位を無次元化した振動子

$$
\ddot q=q-q^3+\eps(-\delta\dot q+F\cos\omega t),\qquad
\delta>0,\ \omega>0
$$

を考える。非摂動の Hamilton 関数は $H=p^2/2-q^2/2+q^4/4$ である。原点は固有値 $\pm1$ の鞍点で、正側ループは

$$
q(t)=\sqrt2\,\operatorname{sech} t,\qquad
p(t)=-\sqrt2\,\operatorname{sech}t\tanh t
$$

である。$\operatorname{sech}''t=\operatorname{sech}t-2\operatorname{sech}^3t$ から方程式を満たし、両時刻端で原点に収束する。負側のループは両成分の符号を反転したものである。

散逸の積分は $v=\tanh t$ で $\int_\R p^2\dd t=2\int_{-1}^1 v^2\dd v=4/3$ である。強制の積分は部分積分で

$$
\int_\R p(t)\sin\omega t\dd t=-\sqrt2\,\omega\int_\R\frac{\cos\omega t}{\cosh t}\dd t
=-\sqrt2\,\pi\omega\,\operatorname{sech}(\pi\omega/2).
$$

最後の積分を求めるため、$e^{i\omega z}/\cosh z$ を頂点 $\pm R,\pm R+i\pi$ の長方形で積分する。縦辺は $\cosh z$ が指数的に大きくなるため $R\to\infty$ で消える。上辺では $\cosh(z+i\pi)=-\cosh z$ で、向きも逆になるため、上下の和は実軸積分の $(1+e^{-\pi\omega})$ 倍である。唯一の極 $z=i\pi/2$ の留数は $e^{-\pi\omega/2}/i$ なので、留数定理で積分は $2\pi e^{-\pi\omega/2}/(1+e^{-\pi\omega})$ となる。これは $\pi\operatorname{sech}(\pi\omega/2)$ に等しい。

従って正側のループの判定式は

$$
M(\theta)=-\frac43\delta+\sqrt2\,\pi\omega F\,\operatorname{sech}(\pi\omega/2)\sin(\omega\theta).
\label{eq-dyn-duffing-melnikov}
$$

負側では正弦項の符号が変わる。$|F|>(4\delta/3)\cosh(\pi\omega/2)/(\sqrt2\pi\omega)$ なら各周期に二つの単純零点がある。逆の厳密不等式なら零点がなく、等号では二次以上を調べる。例えば $\delta=1,\ \omega=1$ の閾値は約 $0.753$ である。$F=1$ としても、この結論は実際の摂動 $\eps F$、$\eps\delta$ が十分小さいときに適用できる。

:::example[高周波と低周波の極限]
$\omega\to\infty$ では $\omega\operatorname{sech}(\pi\omega/2)\sim2\omega e^{-\pi\omega/2}$ となり、固定した $F$ に対する強制の仕事は小さい。$\omega\downarrow0$ でも同じ係数は $\omega$ とともに零へ行く。定数外力なら $\int_\R p(t)\dd t=q(+\infty)-q(-\infty)=0$ となることと一致する。
:::

:::example[ループの側を変える]
$F>0$ で正側の零点が $\theta_0$ なら、負側の対応する零点は $\theta_0+\pi/\omega$ である。散逸項は $p^2$ なのでループを反転しても同じで、外力の仕事だけが符号を変える。このため必要な強制振幅の閾値は両側で同じになる。
:::

### 閾値を周波数の関数として比較する

式 :eqref[eq-dyn-duffing-melnikov] で散逸 $\delta$ を固定し、横断交差に必要な外力の最小の振幅を $F_{\rm c}(\omega)$ とする。$\omega>0$ で $F_{\rm c}$ は $\cosh(\pi\omega/2)/\omega$ に比例する。その対数微分は $(\pi/2)\tanh(\pi\omega/2)-1/\omega$ なので、$u=\pi\omega/2$ とおけば極値条件は $u\tanh u=1$ である。左辺の微分 $\tanh u+u\operatorname{sech}^2u$ は正で、値は $0$ から無限大へ増える。従って正根はただ一つであり、低周波と高周波での発散を合わせると大域最小である。数値的には $u\approx1.19968$、$\omega\approx0.76374$ となる。これは小摂動の一次閾値の最適点であり、有限の $\eps$ における交差条件の極小値まで同じとは限らない。

### 鞍点の曲率も変化させた場合

$\ddot q=a q-q^3$ のエネルギーは $p^2/2-a q^2/2+q^4/4$ であり、曲率 $a>0$ の鞍点へのループのエネルギーは零である。提示した $q_h$ を二回微分して方程式を確認し、$t\to\pm\infty$ で $q_h,p_h\to0$ となることから接続を確認できる。$\sqrt a\,t$ への置換により余弦変換では周波数が $\omega/\sqrt a$ に変わるが、$q_h$ の振幅と $\dd t$ の係数が相殺する。その結果、損失と入力は曲率に対して異なる仕方で変わる。

:::example[障壁の曲率が 1 の Duffing 系]
非摂動の式を $\ddot q=a q-q^3$ とし、右側ループの $q_h(t)=\sqrt{2a}\,\operatorname{sech}(\sqrt a\,t)$ を使う。減衰と外力をそれぞれ $-\gamma p$、$F\cos(\omega(t+\theta)+\varphi)$ とする。

$$
\begin{aligned}
a&=1\\
\gamma&=\frac{1}{4}\\
F&=1\\
\omega&=1\\
\varphi&=0
\end{aligned}
$$

損失積分では $v=\tanh(\sqrt a\,t)$ と置く。すると $p_h^2\dd t=2a^{3/2}v^2\dd v$ である。外力は $p_h=q_h^\prime$ で部分積分し、$q_h$ の偶性により余弦の Fourier 変換だけを残す。

$$
\begin{aligned}
\int_{\R}p_h^2\dd t&=\frac{4}{3}\\
\lvert M_{\rm osc}\rvert_{\max}&=\frac{\sqrt{2} \pi}{\cosh{\left(\frac{\pi}{2} \right)}}\\
M(\theta)&=\frac{\sqrt{2} \pi \sin{\left(\theta \right)}}{\cosh{\left(\frac{\pi}{2} \right)}} - \frac{1}{3}
\end{aligned}
$$

正の振幅で定数損失を割る。比が $1$ 未満なら二つの単純零点をもち、$1$ より大きければどの位相も損失を打ち消せない。

$$
\frac{\gamma\int p_h^2\dd t}{\lvert M_{\rm osc}\rvert_{\max}}=\frac{\sqrt{2} \cosh{\left(\frac{\pi}{2} \right)}}{6 \pi}
$$

この比は $0.18825$ 程度である。従って一周期で二位相の横断的交差を予測する。

外力による振幅を $Y=|M_{\rm osc}|_{\max}$ と書くと、この例での位相微分は $Y\omega\cos(\omega\theta+\varphi)$ である。比が $1$ 未満なので余弦はどちらの零点でも零にならず、二つの交差が得られる。
:::

:::example[障壁の曲率が 4 の Duffing 系]
非摂動の式を $\ddot q=a q-q^3$ とし、右側ループの $q_h(t)=\sqrt{2a}\,\operatorname{sech}(\sqrt a\,t)$ を使う。減衰と外力をそれぞれ $-\gamma p$、$F\cos(\omega(t+\theta)+\varphi)$ とする。

$$
\begin{aligned}
a&=4\\
\gamma&=1\\
F&=1\\
\omega&=1\\
\varphi&=\frac{\pi}{4}
\end{aligned}
$$

損失積分では $v=\tanh(\sqrt a\,t)$ と置く。すると $p_h^2\dd t=2a^{3/2}v^2\dd v$ である。外力は $p_h=q_h^\prime$ で部分積分し、$q_h$ の偶性により余弦の Fourier 変換だけを残す。

$$
\begin{aligned}
\int_{\R}p_h^2\dd t&=\frac{32}{3}\\
\lvert M_{\rm osc}\rvert_{\max}&=\frac{\sqrt{2} \pi}{\cosh{\left(\frac{\pi}{4} \right)}}\\
M(\theta)&=\frac{\sqrt{2} \pi \sin{\left(\theta + \frac{\pi}{4} \right)}}{\cosh{\left(\frac{\pi}{4} \right)}} - \frac{32}{3}
\end{aligned}
$$

正の振幅で定数損失を割る。比が $1$ 未満なら二つの単純零点をもち、$1$ より大きければどの位相も損失を打ち消せない。

$$
\frac{\gamma\int p_h^2\dd t}{\lvert M_{\rm osc}\rvert_{\max}}=\frac{16 \sqrt{2} \cosh{\left(\frac{\pi}{4} \right)}}{3 \pi}
$$

この比は $3.18018$ 程度である。従ってこのループの近傍に対し一次の分裂は符号を変えない。

位相を $\pi/4$ だけずらしても振幅は変わらない。今回比が大きくなった原因は、損失の係数とともに $a^{3/2}$ が増えたことにある。
:::

:::example[障壁の曲率が 9 の Duffing 系]
非摂動の式を $\ddot q=a q-q^3$ とし、右側ループの $q_h(t)=\sqrt{2a}\,\operatorname{sech}(\sqrt a\,t)$ を使う。減衰と外力をそれぞれ $-\gamma p$、$F\cos(\omega(t+\theta)+\varphi)$ とする。

$$
\begin{aligned}
a&=9\\
\gamma&=\frac{1}{16}\\
F&=1\\
\omega&=6\\
\varphi&=\frac{\pi}{3}
\end{aligned}
$$

損失積分では $v=\tanh(\sqrt a\,t)$ と置く。すると $p_h^2\dd t=2a^{3/2}v^2\dd v$ である。外力は $p_h=q_h^\prime$ で部分積分し、$q_h$ の偶性により余弦の Fourier 変換だけを残す。

$$
\begin{aligned}
\int_{\R}p_h^2\dd t&=36\\
\lvert M_{\rm osc}\rvert_{\max}&=\frac{6 \sqrt{2} \pi}{\cosh{\left(\pi \right)}}\\
M(\theta)&=\frac{6 \sqrt{2} \pi \sin{\left(6 \theta + \frac{\pi}{3} \right)}}{\cosh{\left(\pi \right)}} - \frac{9}{4}
\end{aligned}
$$

正の振幅で定数損失を割る。比が $1$ 未満なら二つの単純零点をもち、$1$ より大きければどの位相も損失を打ち消せない。

$$
\frac{\gamma\int p_h^2\dd t}{\lvert M_{\rm osc}\rvert_{\max}}=\frac{3 \sqrt{2} \cosh{\left(\pi \right)}}{16 \pi}
$$

この比は $0.97841$ 程度である。従って一周期で二位相の横断的交差を予測する。

この例で周波数 $6$ は鞍点の指数 $3$ の二倍である。$\cosh$ の引数は $\pi$ になり、減衰係数を小さくしても入力の効果を計算するときにはこの指数的な抑制を含める必要がある。
:::

## 演習問題

::::exercise{level=basic}
$\dot x=\eps(3+\sin t)$ の平均解と元の解の差の上限を求めよ。
:::solution
平均は $3$ なので $\bar x=x_0+3\eps t$ である。直接積分すると $x-\bar x=\eps(1-\cos t)$ となる。$\cos t$ の範囲から絶対値の最大は $2|\eps|$ で、$t=\pi$ で達する。
:::
::::

::::exercise{level=basic}
$\dot x=\eps(-x+\sin t)$ の周期解を求めよ。ただし $\eps>0$ とする。
:::solution
$x=a\sin t+b\cos t$ とすると $-b=\eps(1-a)$、$a=-\eps b$ である。従って $a=\eps^2/(1+\eps^2)$、$b=-\eps/(1+\eps^2)$ となる。任意の別の解との差は $d'=-\eps d$ を満たすため、この周期解は一意で吸引的である。
:::
::::

::::exercise{level=basic}
$H=(q^2+p^2)/2$、$G=(0,-2p)$ の摂動で、非摂動の単位円を一周したときの一次のエネルギー変化を求めよ。
:::solution
$\dot H=-2\eps p^2$ であり、非摂動で $p=-\sin t$ と選べる。従って一次の変化は $-2\eps\int_0^{2\pi}\sin^2 t\dd t=-2\pi\eps$ である。実際の軌道との有限時間の差は $O(\eps)$ なので、この積分に代入して生じる誤差は $O(\eps^2)$ である。
:::
::::

::::exercise{level=basic}
$M(\theta)=2\cos\theta-1$ の単純零点を一周期で求めよ。
:::solution
$\cos\theta=1/2$ から $\theta=\pi/3,5\pi/3$ である。微分は $-2\sin\theta$ なので両点でそれぞれ $-\sqrt3,\sqrt3$、従って単純零点である。
:::
::::

::::exercise{level=advanced}
$\dot x=\eps\sin t$、$\dot y=\eps x\cos t$ で一周期の変位を求め、一次平均との違いを述べよ。
:::solution
$x=x_0+\eps(1-\cos t)$ であり一周期の $\Delta x$ は零である。第二式を積分すると $\int_0^{2\pi}(\eps x_0+\eps^2)\cos t\dd t-\eps^2\int_0^{2\pi}\cos^2t\dd t=-\pi\eps^2$ である。一次平均場は両方零だが、二次で $y$ に一定の変位が残る。
:::
::::

::::exercise{level=advanced}
$G=J\grad K$、$K(q,p)=q^2/2$ とする。Duffing 正側ループに沿う $M$ を求めよ。
:::solution
$\grad H\cdot J\grad K=-pq=-\dv{}{t}(q^2/2)$ である。$q(t)\to0$ が両方向で成り立つので積分は $0$ である。横断的な分裂について :ref[prop-dyn-melnikov-simple] の仮定を満たさない。
:::
::::

::::exercise{level=advanced}
Duffing 正側ループで粘性項を $-\delta p^3$ に替える。散逸の寄与を求めよ。
:::solution
仕事率は $-\delta p^4$ である。$v=\tanh t$ により $p^4\dd t=4(1-v^2)v^4\dd v$ となる。積分は $4(2/5-2/7)=16/35$ なので寄与は $-16\delta/35$ である。
:::
::::

::::exercise{level=advanced}
$\omega=1$ の強制 Duffing で $F$ を選んで正弦項の振幅を $8\delta/3$ にする。単純零点とそこでの位相微分を求めよ。
:::solution
式 :eqref[eq-dyn-duffing-melnikov] で $F=(8\delta/3)\cosh(\pi/2)/(\sqrt2\pi)$ とすればよい。零点は $\sin\theta=1/2$ より $\pi/6,5\pi/6$ である。微分は $(8\delta/3)\cos\theta$ なので両点で $\pm4\sqrt3\delta/3$ となり、$\delta>0$ だから横断性が従う。
:::
::::

::::exercise{level=advanced}
非摂動 Duffing ループ $q_h=\sqrt2\operatorname{sech}t$ に対し、摂動 Hamiltonian を $K=Bq^2\cos(2t)/2$ とする。これに減衰力 $-\gamma p$ を加え、$\gamma>0$、$B>0$ とする。一次の分裂に零点が現れる閾値と、その閾値で零点が単純かを求めよ。さらに閾値の二倍の $B$ で零点の位相を一周期分求めよ。
:::solution
剛性変化の項は $-p_hq_h\cos(2(t+\theta))$ に $B$ を掛けたものを積分する。$q_h^2/2=\operatorname{sech}^2t$ と :ref[lem-dyn-sech-square-transform] によってその積分は $-(4\pi B/\sinh\pi)\sin2\theta$ となる。損失は $-\gamma\int p_h^2\dd t=-4\gamma/3$ である。従って
$$
M=-\frac{4\gamma}{3}-\frac{4\pi B}{\sinh\pi}\sin2\theta
$$
であり、零が存在するための条件は $B\ge B_*=\gamma\sinh\pi/(3\pi)$ である。境界の等号では $\sin2\theta=-1$ でしか零にならず、そのとき $\cos2\theta=0$ なので $M'=0$ である。この場合は一次の符号変化がなく、摂動後の横断交差は一次計算だけでは結論できない。

$B=2B_*$ なら $\sin2\theta=-1/2$ なので $0\le\theta<\pi$ における解は $\theta=7\pi/12,11\pi/12$ である。両方で $|\cos2\theta|=\sqrt3/2$ となるから単純零であり、摂動が十分小さい場合の横断的な分裂の零が各々に続く。
:::
::::

::::exercise{level=basic}
$\dot x=\eps(\sin(2t)+3)$、$x(0)=1$ とする。一周期後の値と、平均解との差の最大値を求めよ。
:::solution
積分して $x(t)=1+3\eps t+\eps(1-\cos2t)/2$ であり、平均解は $\bar x(t)=1+3\eps t$ である。一周期は $\pi$ だから $x(\pi)=1+3\pi\eps$ となる。差の絶対値は $\cos2t=-1$ で最大となり、その値は $|\eps|$ である。
:::
::::

## 確認テスト

::::quiz{answer="2" tol="1e-6"}
$f(t,x)=2+3\cos t$ の $2\pi$ 周期平均を求めよ。
:::explanation
$\int_0^{2\pi}\cos t\dd t=0$ であり、定数項の積分は $4\pi$ である。周期 $2\pi$ で割って $2$ となる。
:::
::::

::::quiz
一次平均化定理の結論を選べ。

- [x] 適切なコンパクト性の下で、時刻 $L/\eps$ までの誤差が $O(\eps)$ である。
- [ ] 平均が零なら元の系の解は定数である。
- [ ] 有限時間の評価だけで全時間での誤差が零になる。

:::explanation
振動成分を除いた残りは $O(\eps^2)$ であり、$L/\eps$ まで積分すると $O(\eps)$ になる。平均零でも周期内に移動し、二次以上で累積しうる。全時間の結論には別の安定性などが必要である。
:::
::::

::::quiz{answer="1.3333333333333333" tol="1e-8"}
Duffing の非摂動ループに沿う $\int_\R p(t)^2\dd t$ を数値で答えよ。
:::explanation
$p^2=2\operatorname{sech}^2(t)\tanh^2(t)$ であり、$v=\tanh t$ と置くと積分は $2\int_{-1}^1 v^2\dd v=4/3$ になる。
:::
::::

::::quiz
$M(\theta)=1-\cos\theta$ で $\theta=0$ の零点について言えることを選べ。

- [x] 単純零点の判定では交差を結論できない。
- [ ] 任意の小さな摂動で必ず横断交差する。
- [ ] 零点が存在しない。

:::explanation
$M(0)=0$ だが $M'(0)=0$ なので、剰余の符号によって零点の個数が変わる余地がある。単純零点の仮定を使う判定はこのとき適用できない。
:::
::::

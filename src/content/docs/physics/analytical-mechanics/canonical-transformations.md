---
title: Poisson 括弧と正準変換
description: Poisson 括弧、作用の境界項、シンプレクティック形式を結び、対称性の生成子を構成する。
notation:
  - symbol: '\{F,G\}'
    meaning: 正準座標に関する Poisson 括弧
  - symbol: 'J'
    meaning: 正準方程式の反対称行列
  - symbol: 'S'
    meaning: 変換の Jacobi 行列
  - symbol: 'F_1,F_2,F_3,F_4'
    meaning: 独立変数の選び方に対応する母関数
---

## この章の目標

- Poisson 括弧の計算規則と Jacobi 恒等式を証明する。
- 四種類の母関数を作り、時間に依存する変換で Hamiltonian を求める。
- 正準性を行列条件と基本括弧の両方で検査する。
- 保存量を無限小変換の生成子として読む。

## Poisson 括弧と保存量

位置と運動量が混ざる新しい座標を選ぶと、複雑な運動が回転や並進として見えることがある。
そこで、同じ正準構造を保つ座標の間で使える時間微分の公式を用意する。

:::definition[Poisson 括弧]{#def-am-poisson}
正準座標 $(q^1,\ldots,q^n,p_1,\ldots,p_n)$ の滑らかな関数 $F,G$ に対して :term[Poisson 括弧]{en="Poisson bracket" yomi="ぽあそんかっこ"} を

$$
\{F,G\}=\sum_{i=1}^n(F_{q^i}G_{p_i}-F_{p_i}G_{q^i})
=(\nabla F)^\top J\nabla G,\qquad J=\begin{pmatrix}0&I\\-I&0\end{pmatrix}
\label{eq-am-poisson}
$$

で定める。
微分は時刻を固定してとる。
:::

:eqref[eq-am-hamilton-equations] を連鎖律に代入すると $\dot F=F_t+\{F,H\}$。従って解の上で一定となる第一積分の条件は $F_t+\{F,H\}=0$ である。座標自身の括弧は $\{q^i,p_j\}=\delta_{ij}$、$\{q^i,q^j\}=\{p_i,p_j\}=0$。括弧は各引数について線形で、交換すると符号が変わり、積に関して $\{FG,K\}=F\{G,K\}+G\{F,K\}$ となる。いずれも各偏微分についての同じ性質から従う。

:::example[面積を掃く運動]
平面の $L_z=xp_y-yp_x$ を第二引数に入れると $\{x,L_z\}=-y,\{y,L_z\}=x$、$\{p_x,L_z\}=-p_y,\{p_y,L_z\}=p_x$。
この四式で $p_x^2+p_y^2$ と $x^2+y^2$ の変化を計算すると両方零となる。
従って中心力 Hamiltonian で $\{L_z,H\}=0$ である。
:::

### Jacobi 恒等式

:::proposition[Jacobi 恒等式と保存量の閉性]{#prop-am-jacobi}
二階偏微分が連続な関数の括弧には

$$
\{F,\{G,K\}\}+\{G,\{K,F\}\}+\{K,\{F,G\}\}=0
$$

が成り立つ。
特に、二つの第一積分の括弧も第一積分である。
:::

:::proof
$z=(q,p)$ とし $F_a=\partial F/\partial z^a$ と略記する。
第一項を展開すると

$$
\sum_{a,b,c,d=1}^{2n}F_a J_{ab}J_{cd}(G_{cb}K_d+G_cK_{db})
$$

である。
巡回和で $G$ の二階微分を含む項を集めると一つは上の $F_aJ_{ab}J_{cd}G_{cb}K_d$ の和、他方は $\{K,\{F,G\}\}$ から来る $\sum_{a,b,c,d}K_dJ_{dc}F_aJ_{ab}G_{bc}$ である。
$G_{cb}=G_{bc}$ と $J_{dc}=-J_{cd}$ により相殺する。
$F,K$ の二階微分も文字を巡回すれば消える。

この恒等式を並べ替えると $\{\{F,G\},H\}=\{\{F,H\},G\}+\{F,\{G,H\}\}$ である。
明示的時間微分の積則も加えると $\dd\{F,G\}/\dd t=\{\dot F,G\}+\{F,\dot G\}$。
両者が保存すれば右辺は零である。
:::

例えば自由粒子の $p$ と $mx-pt$ は保存し、その括弧 $\{p,mx-pt\}=-m$ も一定である。
括弧によって得た量が既知の第一積分から独立とは限らない。

### 三つの回転と括弧の順序

$L_x=yp_z-zp_y$、$L_y=zp_x-xp_z$ とし、$\{L_x,L_y\}$ を計算すると $xp_y-yp_x=L_z$ となる。
添字を巡回すると

$$
\{L_y,L_z\}=L_x,\qquad \{L_z,L_x\}=L_y .
$$

従って $L^2=L_x^2+L_y^2+L_z^2$ は各 $L_i$ と可換である。
例えば $\{L^2,L_x\}=2L_y(-L_z)+2L_z L_y=0$。
大きさが回転で変わらないという幾何が、括弧の積則に現れている。

関数 $F$ に対し演算子 $D_G F=\{F,G\}$ と置くと、二つの変化の差は

$$
(D_GD_K-D_KD_G)F
=\{\{F,K\},G\}-\{\{F,G\},K\}
=\{F,\{K,G\}\}.
$$

最後の引数が $\{K,G\}$ となるのは、この章では生成子を第二引数に置いたためである。
有限回転の順序交換で差が出ることの無限小版になっている。

:::example[外部関数を含む保存量]
一定重力の下で $H=p^2/(2m)+mgq$ とし、

$$
F=p+mgt,\qquad G=q-pt/m-gt^2/2
$$

を考える。
$F_t=mg,\ \{F,H\}=-mg$ だから $F$ が保存する。
$G_t=-p/m-gt$、$\{G,H\}=p/m+gt$ なので $G$ も保存する。
$\{G,F\}=1$ であり、$(G,F)$ は初期の位置と運動量をそのまま与える正準な組である。
:::

### 逆二乗力で保存量の括弧を閉じる

Kepler 問題のベクトル $\vb A$ は回転のもとで通常のベクトルとして変換する。
さらに $\{A_i,A_j\}$ まで計算すると、余分な保存量の代数がエネルギーで変わることが分かる。
:eqref[eq-am-runge-lenz] の定義を使い、三次元で $H=p^2/(2\mu)-k/r$ とし、$p^2=\vb p\cdot\vb p$、$s=\vb r\cdot\vb p$、$f=p^2-\mu k/r$ と略記する。
三重積を展開すれば

$$
A_i=f r_i-s p_i .
$$

基本括弧から $\{r_i,f\}=2p_i,\ \{p_i,f\}=-\mu kr_i/r^3$、$\{r_i,s\}=r_i,\ \{p_i,s\}=-p_i$ である。
さらに $\{f,s\}=\mu k/r-2p^2$。
積則を使い $L_{ij}=r_i p_j-r_jp_i$ とすると

$$
\begin{aligned}
\{fr_i,fr_j\}&=-2f L_{ij},\\
\{sp_i,sp_j\}&=0,\\
-\{fr_i,sp_j\}-\{sp_i,fr_j\}&=-(f+\{f,s\})L_{ij}.
\end{aligned}
$$

最後の係数は $f+\{f,s\}=-p^2$ だから

$$
\{A_i,A_j\}=(-p^2+2\mu k/r)L_{ij}=-2\mu H L_{ij}.
$$

例えば $i=x,\ j=y$ なら $L_{xy}=L_z$。
負エネルギー面上では $\vb D=\vb A/\sqrt{-2\mu H}$ とすると $\{D_x,D_y\}=L_z$ で、回転の関係と同じ符号が現れる。
これに $\{L_x,D_y\}=D_z$ を加えると、$(\vb L+\vb D)/2$ と $(\vb L-\vb D)/2$ はそれぞれ回転の括弧をもち、二組の間の括弧は零となる。
正エネルギー面では $\vb D=\vb A/\sqrt{2\mu H}$ と置くと $\{D_x,D_y\}=-L_z$ となる。$\{L_x,D_y\}=D_z$ は変わらない。この符号は、空間回転に三つの Lorentz ブーストを加えた $\mathfrak{so}(3,1)$ の括弧に対応する。負エネルギー側の二組の回転は $\mathfrak{so}(4)$ に対応する。いずれも $H$ は $\vb L,\vb A$ と括弧が零だから、規格化因子の微分による追加項は出ない。

$H=0$ では規格化せず、$\{A_i,A_j\}=0$ を使う。$\vb A$ の三成分は互いに可換で、$\vb L$ の括弧に対してベクトルの成分として回る。これは回転と並進の Lie 代数である。これらは保存量の括弧の関係であり、各生成子を相空間の全点で任意の時間だけ積分できるかは別の条件である。

### 拘束面へ括弧を制限する

二つの関数が拘束面上で同じ値をもっていても、面外での延長の仕方が違うと通常の Poisson 括弧は一致しない。例えば $p=0$ の直線上で $p$ と零は等しいが、$q$ との括弧はそれぞれ $1$ と零である。拘束 $\chi_i=0$ の相互の括弧 $C_{ij}=\{\chi_i,\chi_j\}$ が正則なとき、この違いを消す操作がある。独立な拘束を取り、$C^{-1}$ の定義域をその近傍に限って
$$
\{F,G\}_D=\{F,G\}
-\sum_{i,j}\{F,\chi_i\}(C^{-1})_{ij}\{\chi_j,G\}
\label{eq-am-dirac-bracket}
$$
と置く。これを :term[Dirac 括弧]{en="Dirac bracket" yomi="でぃらっくかっこ"} という。第一引数に $\chi_k$ を入れれば $\sum_i C_{ki}(C^{-1})_{ij}=\delta_{kj}$ で全体が零になる。従って $F$ に $\sum_i f_i\chi_i$ を加えても、積の微分則により拘束面上での結果は変わらない。

この式は力学からも導ける。全 Hamiltonian を $H_T=H+\sum_i v_i\chi_i$ とし、時間に依存しない拘束を保つ条件を課すと、面上で $\{\chi_j,H\}+\sum_i C_{ji}v_i=0$。これを $v_i$ について解き $\dot F=\{F,H_T\}+F_t$ に代入すれば $\dot F=\{F,H\}_D+F_t$ を得る。面上での二形式は非退化である。実際、面の接空間にシンプレクティック直交する方向は各拘束の Hamilton ベクトル場の張る空間であり、その組合せが再び接するなら $C$ を掛けた係数列が零になるため係数はすべて零である。引き戻した二形式と括弧の関係は、場の射影から確かめられる。$X_G=J\nabla G$ に $\sum_i v_i X_{\chi_i}$ を加え、$v_i=-\sum_j(C^{-1})_{ij}\{\chi_j,G\}$ とすれば接ベクトルになる。接ベクトル $Y$ に対し $\Omega(X_{\chi_i},Y)=\dd\chi_i(Y)=0$ だから、射影した場は拘束面上の $G$ の Hamilton ベクトル場である。この場で $F$ を微分した値が Dirac 括弧となる。後で使う符号規約では、二形式の成分行列を $B$ とすると Poisson 行列は $-B^{-1}$ である。局所成分で逆行列を微分する公式 $\partial B^{-1}=-B^{-1}(\partial B)B^{-1}$ を使えば、Jacobi 恒等式の三項は二形式の外微分の三項に一致して消える。

### 球面上の位置と運動量の成分

半径 $a$ の球面上の粒子を三組の Cartesian 正準座標で表し、$\chi_1=(\vb r^2-a^2)/2,\ \chi_2=\vb r\cdot\vb p$ とする。第二の条件は速度が接平面内にあることに相当する。$\{\chi_1,\chi_2\}=\vb r^2$ だから逆行列は球面上で存在する。$\{r_i,\chi_1\}=0,\{r_i,\chi_2\}=r_i,\{p_i,\chi_1\}=-r_i,\{p_i,\chi_2\}=-p_i$ を代入すると
$$
\{r_i,r_j\}_D=0,\quad \{r_i,p_j\}_D=\delta_{ij}-r_ir_j/a^2,\quad
\{p_i,p_j\}_D=(p_i r_j-r_i p_j)/a^2 .
$$
混合括弧は接平面への射影行列である。$\{r_i,p_j\}_D$ が単位行列でないのは、三成分すべてを自由には変えられないためである。$H=\vb p^2/(2m)+V$ に適用すると $\dot{\vb r}=\vb p/m,\ \dot{\vb p}=-(1-\vb r\vb r^\top/a^2)\nabla V-\vb r\,\vb p^2/(ma^2)$。法線加速度は瞬間の運動エネルギーを変えず、速度が変化する接平面にとどまるための曲率を与える。:ref[prop-am-ideal-multiplier] の規約で法線力を $+\lambda\vb r$ と置いて $\ddot\chi_1=0$ を解くことでも同じ式が出る。

:::example[北極を通過する自由粒子]
位置を $\vb r=\left(\begin{matrix}0 & 0 & 2\end{matrix}\right)^\top\,\mathrm m$、接線運動量を $\vb p=\left(\begin{matrix}3 & 4 & 0\end{matrix}\right)^\top\,\mathrm{kg\,m/s}$ とする。$m=2\,\mathrm{kg}$ で $\nabla V=\left(\begin{matrix}0 & 0 & 0\end{matrix}\right)^\top\,\mathrm N$。$\vb r\cdot\vb p=0$ を確かめてから、混合括弧の行列と運動量同士の括弧を計算すると

$$
(\{r_i,p_j\}_D)=\left(\begin{matrix}1 & 0 & 0\\0 & 1 & 0\\0 & 0 & 0\end{matrix}\right),\qquad(\{p_i,p_j\}_D)=\left(\begin{matrix}0 & 0 & \frac{3}{2}\\0 & 0 & 2\\- \frac{3}{2} & -2 & 0\end{matrix}\right).
$$

ここから得る $\dot{\vb p}$ は $\left(\begin{matrix}0 & 0 & - \frac{25}{4}\end{matrix}\right)^\top\,\mathrm N$、元の外力を引いた拘束力は $\left(\begin{matrix}0 & 0 & - \frac{25}{4}\end{matrix}\right)^\top\,\mathrm N$。運動量の法線成分の時間微分には法線自身の変化も入る。実際に $\vb r\cdot\dot{\vb p}=- \frac{25}{2}$ と $\dot{\vb r}\cdot\vb p=\frac{25}{2}$ が相殺し、接線条件が保たれる。球座標はこの点で方位の目盛が定まらないが、射影行列と拘束力は有限である。
:::

:::example[赤道で鉛直荷重を受ける質点]
位置を $\vb r=\left(\begin{matrix}2 & 0 & 0\end{matrix}\right)^\top\,\mathrm m$、接線運動量を $\vb p=\left(\begin{matrix}0 & 0 & 2\end{matrix}\right)^\top\,\mathrm{kg\,m/s}$ とする。$m=1\,\mathrm{kg}$ で $\nabla V=\left(\begin{matrix}0 & 0 & 10\end{matrix}\right)^\top\,\mathrm N$。$\vb r\cdot\vb p=0$ を確かめてから、混合括弧の行列と運動量同士の括弧を計算すると

$$
(\{r_i,p_j\}_D)=\left(\begin{matrix}0 & 0 & 0\\0 & 1 & 0\\0 & 0 & 1\end{matrix}\right),\qquad(\{p_i,p_j\}_D)=\left(\begin{matrix}0 & 0 & -1\\0 & 0 & 0\\1 & 0 & 0\end{matrix}\right).
$$

ここから得る $\dot{\vb p}$ は $\left(\begin{matrix}-2 & 0 & -10\end{matrix}\right)^\top\,\mathrm N$、元の外力を引いた拘束力は $\left(\begin{matrix}-2 & 0 & 0\end{matrix}\right)^\top\,\mathrm N$。運動量の法線成分の時間微分には法線自身の変化も入る。実際に $\vb r\cdot\dot{\vb p}=-4$ と $\dot{\vb r}\cdot\vb p=4$ が相殺し、接線条件が保たれる。鉛直力はこの点ではすべて接線方向である。それでも球面を外れないための法線力が速度から生じる。
:::

:::example[斜めの法線に対する拘束反力]
位置を $\vb r=\left(\begin{matrix}1 & 2 & 2\end{matrix}\right)^\top\,\mathrm m$、接線運動量を $\vb p=\left(\begin{matrix}2 & -1 & 0\end{matrix}\right)^\top\,\mathrm{kg\,m/s}$ とする。$m=2\,\mathrm{kg}$ で $\nabla V=\left(\begin{matrix}0 & 0 & 6\end{matrix}\right)^\top\,\mathrm N$。$\vb r\cdot\vb p=0$ を確かめてから、混合括弧の行列と運動量同士の括弧を計算すると

$$
(\{r_i,p_j\}_D)=\left(\begin{matrix}\frac{8}{9} & - \frac{2}{9} & - \frac{2}{9}\\- \frac{2}{9} & \frac{5}{9} & - \frac{4}{9}\\- \frac{2}{9} & - \frac{4}{9} & \frac{5}{9}\end{matrix}\right),\qquad(\{p_i,p_j\}_D)=\left(\begin{matrix}0 & \frac{5}{9} & \frac{4}{9}\\- \frac{5}{9} & 0 & - \frac{2}{9}\\- \frac{4}{9} & \frac{2}{9} & 0\end{matrix}\right).
$$

ここから得る $\dot{\vb p}$ は $\left(\begin{matrix}\frac{19}{18} & \frac{19}{9} & - \frac{35}{9}\end{matrix}\right)^\top\,\mathrm N$、元の外力を引いた拘束力は $\left(\begin{matrix}\frac{19}{18} & \frac{19}{9} & \frac{19}{9}\end{matrix}\right)^\top\,\mathrm N$。運動量の法線成分の時間微分には法線自身の変化も入る。実際に $\vb r\cdot\dot{\vb p}=- \frac{5}{2}$ と $\dot{\vb r}\cdot\vb p=\frac{5}{2}$ が相殺し、接線条件が保たれる。荷重が法線にも接線にも射影をもつので、外力の大きさをそのまま接線方向の加速に使えない。
:::

## 正準変換と母関数

### 作用の積分が境界だけ変わる場合

各時刻で $\sum_i\dd q^i\wedge\dd p_i=\sum_i\dd Q^i\wedge\dd P_i$ を保つ可逆で滑らかな変換を :term[正準変換]{en="canonical transformation" yomi="せいじゅんへんかん"} という。
この二形式は、各位置・運動量の組で測る有向面積を足したものである。これを保つと正準方程式の形も保たれることを以下で示す。
:eqref[eq-am-phase-action] の被積分一形式について、局所的に

$$
\sum_i p_i\dd q^i-H\dd t
=\sum_i P_i\dd Q^i-K\dd t+\dd F_1
\label{eq-am-generator-differential}
$$

が成り立つなら、内部で台をもつ変分に対して二つの作用の変化は同じであり、新しい運動は $K$ の正準方程式に従う。
$F_1$ を $(q,Q,t)$ で表せる領域では係数比較により $p_i=\partial_{q^i}F_1,\ P_i=-\partial_{Q^i}F_1,\ K=H+\partial_t F_1$ である。
最後の式の $H$ は新変数で表し直す。

### 四つの独立変数の選び方

必要に応じて $F_1$ に $\sum_iQ^iP_i$ を加え、または $\sum_iq^ip_i$ を引く。
微分を展開すると $P_i\dd Q^i$ と $p_i\dd q^i$ が相殺され、選んだ変数だけの微分が残る。
この関数を :term[母関数]{en="generating function" yomi="ぼかんすう"} という。
下表で下付き文字の偏微分は、同じ行の他の独立変数を固定する。

| 型 | 基本微分から作る母関数 | 旧変数の式 | 新変数の式 |
| --- | --- | --- | --- |
| 1 | $F_1(q,Q,t)$ | $p_i=F_{1,q^i}$ | $P_i=-F_{1,Q^i}$ |
| 2 | $F_2(q,P,t)=F_1+\sum_i Q^iP_i$ | $p_i=F_{2,q^i}$ | $Q^i=F_{2,P_i}$ |
| 3 | $F_3(p,Q,t)=F_1-\sum_iq^ip_i$ | $q^i=-F_{3,p_i}$ | $P_i=-F_{3,Q^i}$ |
| 4 | $F_4(p,P,t)=F_1+\sum_i(Q^iP_i-q^ip_i)$ | $q^i=-F_{4,p_i}$ | $Q^i=F_{4,P_i}$ |

各行で $K=H+F_{j,t}$ となる。
例えば第二型の微分は $\dd F_2=\sum_i(p_i\dd q^i+Q^i\dd P_i)+(K-H)\dd t$ である。
逆に任意に選んだ $F_2$ が局所変換を定めるには、$p=F_{2,q}$ を $P$ について解ける必要がある。
従って混合 Hessian $(F_{2,q^iP_j})$ の行列式が非零ならよい。
他の型についても表の対応する混合微分が同じ役割を果たす。
恒等変換では $Q=q$ のため $(q,Q)$ は独立にならないが、第二型 $F_2=\sum_iq^iP_i$ は使える。

:::example[尺度の変更と位置・運動量の交換]
一自由度で $F_2=aqP$、$a\ne0$ とすれば $Q=aq,\ P=p/a$。
$a=2$ なら位置方向の区間は二倍、運動量方向は半分で面積が同じになる。
一方 $F_1=cqQ$ とすれば $Q=p/c,\ P=-cq$。
$Q$ を長さとして使うには $c$ に質量／時間の単位を持たせる。
振動子では $c=m\omega$ とするとこの変換は周期の四分の一の時間発展になる。
:::

:::example[動く原点]
$F_2=(q-a(t))P$ から $p=P,\ Q=q-a$。
$H=p^2/(2m)+V(q)$ なら

$$
K=P^2/(2m)+V(Q+a)-\dot a P .
$$

$\dot Q=P/m-\dot a$ は動く原点から見た速度である。
例えば原点が $2\,\mathrm{m/s}$ で走り粒子が $5\,\mathrm{m/s}$ なら相対速度は $3\,\mathrm{m/s}$。
$F_{2,t}$ を落とすとこの差を再現できない。
:::

### 自由度ごとに独立変数を変える

二自由度で第一の正準対をそのままにし、第二の対を交換して $Q^1=q^1,\ P_1=p_1,\ Q^2=p_2,\ P_2=-q^2$ とする。第一型の独立変数には $q^1,Q^1$ が重複し、第二型には $q^2,P_2$ が重複する。同様に、第三型では $p_2,Q^2$、第四型では $p_1,P_1$ が独立でない。そこで第一の対には第二型、第二の対には第一型を使う。混合した母関数は
$$
F(q^1,P_1,q^2,Q^2)=q^1P_1+q^2Q^2
$$
である。微分して $\dd F=p_1\dd q^1+Q^1\dd P_1+p_2\dd q^2-P_2\dd Q^2$ の係数と比較すると、四つの変換式がすべて戻る。多自由度では四種類のうち一つを全自由度に適用する方法だけでは足りず、各対の独立性を確かめて選ぶ必要がある。

### 同じ回転を四種類の母関数で表す

位置・運動量を一定の基準値で規格化し、作用の単位もその積にとる。
無次元の回転

$$
Q=q\cos\alpha+p\sin\alpha,\qquad
P=-q\sin\alpha+p\cos\alpha
$$

を考える。
$\sin\alpha\ne0$ なら $p=(Q-q\cos\alpha)/\sin\alpha$ を $q$ で積分し、$P=-F_{1,Q}$ で $Q$ だけの積分定数を決める。
結果は

$$
F_1=\frac{qQ-\frac12(q^2+Q^2)\cos\alpha}{\sin\alpha}.
$$

同じ手順で、それぞれの独立変数が使える範囲において

$$
\begin{aligned}
F_2&=\frac{qP+\frac12(q^2+P^2)\sin\alpha}{\cos\alpha},\\
F_3&=\frac{-pQ+\frac12(p^2+Q^2)\sin\alpha}{\cos\alpha},\\
F_4&=\frac{pP-\frac12(p^2+P^2)\cos\alpha}{\sin\alpha}.
\end{aligned}
$$

例えば第四型では $q=(p\cos\alpha-P)/\sin\alpha$、$Q=(p-P\cos\alpha)/\sin\alpha$ となり、元の式を $p,P$ について解いたものと一致する。
$\alpha\to0$ では第二型は $qP$、第三型は $-pQ$ に近づくが、第一・第四型の分母は消える。
逆に $\alpha=\pi/2$ の交換には第一・第四型が使える。
母関数の特異性は、必ずしも変換の特異性を意味しない。

### 点変換で運動量を決める

$q\mapsto Q=f(q,t)$ を与えたとき、第二型を

$$
F_2(q,P,t)=\sum_{j=1}^n f^j(q,t)P_j
$$

に選べば、$p_i=\sum_jP_j\partial_i f^j$ となり、$K=H+\sum_jP_j\partial_t f^j$ が得られる。
位置だけを新座標に取り替え、運動量に同じ変換則を機械的に用いると面積を保てない。

:::example[直交座標から極座標へ]
$r=\sqrt{x^2+y^2}$、$\phi$ は極角の滑らかな枝とし $F_2=rP_r+\phi P_\phi$ と置く。
$\partial_x r=\cos\phi,\ \partial_x\phi=-\sin\phi/r$ などを用い、

$$
p_x=P_r\cos\phi-\frac{P_\phi}{r}\sin\phi,\quad
p_y=P_r\sin\phi+\frac{P_\phi}{r}\cos\phi
$$

を得る。
従って $P_r=p_x\cos\phi+p_y\sin\phi,\ P_\phi=xp_y-yp_x$。
$p_x^2+p_y^2=P_r^2+P_\phi^2/r^2$ なので Hamiltonian に遠心項が現れる。
$\phi$ の枝を替えても微分は一致するが、$r=0$ では極座標そのものが使えない。
:::

### 点変換の余接持ち上げと計量

曲線座標での運動量を作るには速度のように変換してはいけない。新座標を $Q$、旧位置を $q=f(Q)$ と書き、$D=\partial f/\partial Q$ を正則とする。作用に現れる一次形式の等式 $p^\top\dd q=P^\top\dd Q$ をすべての $\dd Q$ に課すと $P=D^\top p$。これは $p=D^{-\top}P$ と解ける。この変換を :term[余接持ち上げ]{en="cotangent lift" yomi="よせつもちあげ"} という。旧運動量と新位置の母関数を使うなら $F_3(p,Q)=-p^\top f(Q)$。微分すれば $q=-\partial_p F_3,\ P=-\partial_Q F_3$ である。自由粒子の $H=p^\top p/(2m)$ は $P^\top(D^\top D)^{-1}P/(2m)$ に移る。変換後に交差項が現れることは、元の粒子に新しい相互作用が生じたことを意味しない。

:::example[せん断した直交格子]
格子線は直線だが直交しないため、計量に交差項が出る。物理量は基準の長さ $1\,\mathrm m$、質量 $1\,\mathrm{kg}$、時間 $1\,\mathrm s$ で無次元化する。新座標を $(a,b)$ と書く。位置と Jacobi 行列は

$$
\binom{x}{y}=\left(\begin{matrix}a + 2 b\\b\end{matrix}\right),\qquad D=\left(\begin{matrix}1 & 2\\0 & 1\end{matrix}\right),\quad \det D=1.
$$

従って $D^\top D=\left(\begin{matrix}1 & 2\\2 & 5\end{matrix}\right)$。新運動量を $(P_a,P_b)$ と書けば、運動エネルギーは

$$
H=\tfrac12 (P_a,P_b)(D^\top D)^{-1}(P_a,P_b)^\top=\frac{5 P_a^{2} - 4 P_a P_b + P_b^{2}}{2}.
$$

$(a,b)=(1,2)$ で $(p_x,p_y)=(2,-1)$ のとき、$(P_a,P_b)^\top=D^\top p=\left(\begin{matrix}2\\3\end{matrix}\right)$。元の座標でも新座標でも運動エネルギーは $\frac{5}{2}\,\mathrm J$ である。逆行列の転置を位置変換と同じ行列に置き換えるとこの一致が失われる。
:::

:::example[放物線座標の一枚の座標域]
座標線は直交する。原点を除いても二つの符号を同時に変えると同じ点に移るため、逆写像の枝を一枚に選ぶ。物理量は基準の長さ $1\,\mathrm m$、質量 $1\,\mathrm{kg}$、時間 $1\,\mathrm s$ で無次元化する。新座標を $(a,b)$ と書く。位置と Jacobi 行列は

$$
\binom{x}{y}=\left(\begin{matrix}\frac{a^{2}}{2} - \frac{b^{2}}{2}\\a b\end{matrix}\right),\qquad D=\left(\begin{matrix}a & - b\\b & a\end{matrix}\right),\quad \det D=a^{2} + b^{2}.
$$

従って $D^\top D=\left(\begin{matrix}a^{2} + b^{2} & 0\\0 & a^{2} + b^{2}\end{matrix}\right)$。新運動量を $(P_a,P_b)$ と書けば、運動エネルギーは

$$
H=\tfrac12 (P_a,P_b)(D^\top D)^{-1}(P_a,P_b)^\top=\frac{P_a^{2} + P_b^{2}}{2 \left(a^{2} + b^{2}\right)}.
$$

$(a,b)=(1,2)$ で $(p_x,p_y)=(2,1)$ のとき、$(P_a,P_b)^\top=D^\top p=\left(\begin{matrix}4\\-3\end{matrix}\right)$。元の座標でも新座標でも運動エネルギーは $\frac{5}{2}\,\mathrm J$ である。
:::

:::example[指数で半径を表した極座標]
二方向の尺度因子は同じで、半径に比例する。角度の周期をまたいで一つの実数の逆関数を定めることはできない。物理量は基準の長さ $1\,\mathrm m$、質量 $1\,\mathrm{kg}$、時間 $1\,\mathrm s$ で無次元化する。新座標を $(a,b)$ と書く。位置と Jacobi 行列は

$$
\binom{x}{y}=\left(\begin{matrix}e^{a} \cos{\left(b \right)}\\e^{a} \sin{\left(b \right)}\end{matrix}\right),\qquad D=\left(\begin{matrix}e^{a} \cos{\left(b \right)} & - e^{a} \sin{\left(b \right)}\\e^{a} \sin{\left(b \right)} & e^{a} \cos{\left(b \right)}\end{matrix}\right),\quad \det D=e^{2 a}.
$$

従って $D^\top D=\left(\begin{matrix}e^{2 a} & 0\\0 & e^{2 a}\end{matrix}\right)$。新運動量を $(P_a,P_b)$ と書けば、運動エネルギーは

$$
H=\tfrac12 (P_a,P_b)(D^\top D)^{-1}(P_a,P_b)^\top=\frac{\left(P_a^{2} + P_b^{2}\right) e^{- 2 a}}{2}.
$$

$(a,b)=(0,\frac{\pi}{2})$ で $(p_x,p_y)=(3,-1)$ のとき、$(P_a,P_b)^\top=D^\top p=\left(\begin{matrix}-1\\-3\end{matrix}\right)$。元の座標でも新座標でも運動エネルギーは $5\,\mathrm J$ である。
:::

:::example[曲がった溝に沿う座標]
元の位置は全平面で一対一である。測定した運動量の混合は、溝の接線の傾きに応じて位置とともに変わる。物理量は基準の長さ $1\,\mathrm m$、質量 $1\,\mathrm{kg}$、時間 $1\,\mathrm s$ で無次元化する。新座標を $(a,b)$ と書く。位置と Jacobi 行列は

$$
\binom{x}{y}=\left(\begin{matrix}a\\\frac{a^{2}}{2} + b\end{matrix}\right),\qquad D=\left(\begin{matrix}1 & 0\\a & 1\end{matrix}\right),\quad \det D=1.
$$

従って $D^\top D=\left(\begin{matrix}a^{2} + 1 & a\\a & 1\end{matrix}\right)$。新運動量を $(P_a,P_b)$ と書けば、運動エネルギーは

$$
H=\tfrac12 (P_a,P_b)(D^\top D)^{-1}(P_a,P_b)^\top=\frac{P_a^{2} - 2 P_a P_b a + P_b^{2} a^{2} + P_b^{2}}{2}.
$$

$(a,b)=(2,1)$ で $(p_x,p_y)=(1,2)$ のとき、$(P_a,P_b)^\top=D^\top p=\left(\begin{matrix}5\\2\end{matrix}\right)$。元の座標でも新座標でも運動エネルギーは $\frac{5}{2}\,\mathrm J$ である。
:::

## シンプレクティック条件

### 面積を測る二形式

相空間の二つの微小変位 $u,v$ に対する有向面積の和を

$$
\Omega(u,v)=u^\top Jv=\sum_i(u_{q^i}v_{p_i}-u_{p_i}v_{q^i})
$$

とする。
これを二形式で $\Omega=\sum_i\dd q^i\wedge\dd p_i$ と書く。
固定時刻で :eqref[eq-am-generator-differential] の外微分をとると $\dd^2F=0$ より旧新の二形式が等しい。
逆に二形式が一致すれば、一形式の差は閉じている。
小さい星形の座標近傍では閉一形式 $\sum_i a_i(z)\dd z^i$ は、原点を近傍の中心に選んで $F(z)=\int_0^1\sum_i z^i a_i(sz)\dd s$ の微分になる。
実際、$\partial_j a_i=\partial_i a_j$ より積分を $z^j$ で微分した被積分関数は $\dd[sa_j(sz)]/\dd s$ である。
大域的には閉曲線に沿う積分値も確かめる必要がある。

:::proposition[シンプレクティック条件と括弧の不変性]{#prop-am-symplectic}
$S=\partial(Q,P)/\partial(q,p)$ が正則なら、固定時刻で次は同値である。

1. あらゆる $u,v$ に対して $\Omega(Su,Sv)=\Omega(u,v)$。
2. $S^\top JS=J$。
3. 新変数の基本括弧が旧変数で計算しても正準形である。
4. 任意の関数の Poisson 括弧を旧新どちらの座標で計算しても一致する。
:::

:::proof
第一条件の左辺は $u^\top S^\top JSv$ なので、標準基底の全組を代入すれば第二条件と同値である。
$J^{-1}=-J$ を使って逆行列をとり、左から $S$、右から $S^\top$ を掛けると第二条件は $SJS^\top=J$ と同値になる。
この左辺の $(a,b)$ 成分は新変数 $Z^a,Z^b$ の旧変数での括弧であり第三条件を与える。
さらに $\nabla_z F=S^\top\nabla_Z F$ より括弧の変換で両勾配の間に現れる行列は $SJS^\top$。
従って第三条件から第四条件、座標関数に適用して逆も従う。
:::

このように二形式を保つことを :term[シンプレクティック]{en="symplectic" yomi="しんぷれくてぃっく"} という。
$\Omega$ を $n$ 回くさび積した形は $n!$ と符号を除いて正準体積要素であるため、正準変換は向きを含めて体積を保つ。
自由度が一つなら $\Omega$ は面積そのもので、$2\times2$ 行列では $S^\top JS=(\det S)J$ となる。

:::example[体積は保つが正準でない変換]
二自由度で $(Q^1,Q^2,P_1,P_2)=(2q^1,q^2,p_1,p_2/2)$。
行列式は $1$ だが $\{Q^1,P_1\}=2$。
一対の面積を拡げ、他方の面積を縮めたため全体積だけでは検出できない。
:::

```tikz caption="同じ正準面積を保つ尺度変換。横幅を二倍にすると共役運動量の幅は半分になる。"
\draw[->] (0,0) -- (2.5,0) node[right] {$q$};
\draw[->] (0,0) -- (0,2) node[above] {$p$};
\draw[thick] (.3,.5) rectangle (1.3,1.5);
\draw[->,thick] (2.9,1) -- (3.8,1);
\draw[->] (4,0) -- (7.2,0) node[right] {$Q$};
\draw[->] (4,0) -- (4,2) node[above] {$P$};
\draw[thick] (4.6,.25) rectangle (6.6,.75);
```

### ブロック行列で検査する

$S=\begin{pmatrix}A&B\\C&D\end{pmatrix}$ と分け、$S^\top JS$ を実際に掛けると、必要十分条件は

$$
A^\top C=C^\top A,\qquad B^\top D=D^\top B,\qquad
A^\top D-C^\top B=I .
$$

残る一つのブロックは第三式の転置である。
例えば $Q=q,\ P=p+Cq$ では $A=D=I,\ B=0$ だから、$C=C^\top$ の場合に限り正準である。
$C$ が対称なら $F_2=q^\top P-q^\top Cq/2$ が実際にこの変換を与える。

一方、$Q=q,\ P=p-eA(q)$ として力学的運動量を座標にしたいとする。
各成分の括弧は

$$
\{P_i,P_j\}=e(\partial_i A_j-\partial_j A_i).
$$

磁場があると右辺は一般に零でなく、$(q,P)$ の基本括弧は正準形でない。
電磁場の計算でどの運動量を選んだかを明記すべき理由である。

### Hamilton 相流が二形式を保つ直接証明

:ref[thm-am-liouville-volume] で使った相流の微分を $S(t)$ とする。変分方程式は $\dot S=JB S$、$B=H''(z(t),t)=B^\top$。従って

$$
\frac{\dd}{\dd t}(S^\top JS)
=S^\top[(JB)^\top J+J(JB)]S
=S^\top(B-B)S=0 .
$$

初期時刻の $S=I$ より $S^\top JS=J$。
$H$ が時間依存でも成立する。
すなわち時間発展はすべての有向面積の組を保存し、体積保存はその帰結である。

:::example[一歩だけ時間を進める近似]
$H=p^2/(2m)+kq^2/2$ を時間幅 $h$ だけ進める。
両方の右辺を古い時刻で評価する陽的 Euler 法の行列は

$$
S_{\rm E}=\begin{pmatrix}1&h/m\\-hk&1\end{pmatrix},\qquad
\det S_{\rm E}=1+h^2k/m .
$$

初期点の小領域は一歩ごとに拡大する。
これに対し $P=p-hkq$ を先に求め、続いて $Q=q+hP/m$ とすれば

$$
S_{\rm SE}=\begin{pmatrix}1-h^2k/m&h/m\\-hk&1\end{pmatrix}
$$

となり行列式は $1$。
二つのせん断を合成したこの方法は一自由度で正準である。
$q=1,p=0,m=k=1,h=0.1$ を基準単位で入れると後者のエネルギーは初期 $0.5$ から $0.49505$ へ変わる。
正準性とエネルギーの厳密保存は独立に検査すべき性質である。
:::

### 変換が時間に依存する場合にも同じ行列条件を使えるか

時刻ごとに $Z=f(z,t)$ がシンプレクティックであるとする。
時間発展の速度は $\dot Z=S\dot z+f_t$ であり、Jacobi 行列の条件だけを使って第二項を忘れると、新しい運動を誤る。
時刻を固定して旧 Hamiltonian を新座標で $\tilde H(Z,t)=H(f^{-1}(Z,t),t)$ と表す。
勾配の変換と $SJS^\top=J$ により

$$
S\dot z=SJ\nabla_z H=J\nabla_Z\tilde H .
$$

残る $Y(Z,t)=f_t(f^{-1}(Z,t),t)$ は、動く座標自体の速度である。
各時刻に正準な $f$ を微分して $\partial_tS=(D_ZY)S$。
$S^\top JS=J$ の時間微分を取ると

$$
(D_ZY)^\top J+J D_ZY=0 .
$$

従って行列 $-J D_ZY$ は対称であり、列 $-JY$ の各成分は混合微分の整合条件を満たす。
星形の座標近傍では局所関数 $g$ が存在して $\nabla_Z g=-JY$。
ゆえに $Y=J\nabla_Zg$、$K=\tilde H+g$ なら新変数でも正準方程式が成立する。
$g$ に時刻だけの関数を足す自由度は残る。
母関数が使える領域ではこの $g$ が $F_{j,t}$ を新座標で表したものになる。

:::example[動く座標で運動量だけをずらす]
$Q=q,\ P=p-a(t)$ とすると $Y=(0,-\dot a)$。
従って $g$ は $g_P=0,\ -g_Q=-\dot a$ を満たし $g=\dot a Q$ と選べる。
$K=H(Q,P+a,t)+\dot a Q$ の右辺から計算すれば $\dot P=-H_q-\dot a$ となる。
実際 $F_2=q(P+a)$ を使えば $F_{2,t}=q\dot a$ で、同じ項が戻る。
:::

### 周期座標での Hamiltonian の存在

角 $q\pmod{2\pi}$ をもつ円筒で $Q=q,\ P=p-a(t)$ とする。各時刻でこの変換は正準である。ところが $\dot a\ne0$ のとき、動く座標に由来する部分は局所的に $g=\dot a Q$ となり、一周すると $2\pi\dot a$ だけ増える。実際、対応する閉一形式 $\dot a\,\dd Q$ の一周積分が零でないため、大域的に一価な $g$ は存在しない。元の $H$ が一価でも $K=\tilde H+g$ は角をほどいた座標域ごとに作る必要がある。一方、$a$ が定数なら $Y=0$ なので $K=\tilde H$ は大域的に使える。その場合でも旧新の正準一形式の差 $a\,\dd Q$ は $a\ne0$ なら大域的な全微分にならない。動く座標での Hamiltonian の存在と、固定時刻での大域的な母関数の存在は、別々の閉一形式の周期で決まる。

### 局所的に母関数を積分して作る

一自由度で $Q=q+p,\ P=p$ を、特定の母関数を先に仮定せずに調べる。
$q,Q$ が独立なため $p=P=Q-q$ と書ける。
旧新の一形式の差は

$$
p\,\dd q-P\,\dd Q=(Q-q)(\dd q-\dd Q).
$$

従って $F_1=-(Q-q)^2/2$ である。
第二型を使うなら $F_2=F_1+QP=qP+P^2/2$ となる。
これらは見た目は別の関数だが、異なる端点データから同じせん断写像を生成する。
$q,p$ に異なる SI 単位を保つ場合には、$Q=q+cp$ とし $c$ に長さ／運動量の単位を付ければよい。

また $\int(p\,\dd q-P\,\dd Q)$ を積分して母関数を作る際には、経路を変えても積分が変わらないことが必要である。
例えば $Q=q+p,\ P=2p$ では $S$ の行列式が $2$。
$q,Q$ で表した一形式は $(Q-q)\dd q-2(Q-q)\dd Q$ であり、$\dd q$ の係数の $Q$ 微分は $1$、$\dd Q$ の係数の $q$ 微分は $2$。
混合微分が一致しないのでどの局所関数の全微分にもならない。

### Lagrange 括弧と積分不変量

正準座標をパラメータ $u^a$ の関数として与えた場合、逆写像を解かずに正準性を調べたい。$z=(q,p)$ とし $S=\partial z/\partial u$ を作る。二形式 $\sum_i\dd q^i\wedge\dd p_i$ の $u$ 表示の係数
$$
[u^a,u^b]_L=\sum_i\left(\frac{\partial q^i}{\partial u^a}
\frac{\partial p_i}{\partial u^b}-
\frac{\partial q^i}{\partial u^b}\frac{\partial p_i}{\partial u^a}\right)
$$
を :term[Lagrange 括弧]{en="Lagrange bracket" yomi="らぐらんじゅかっこ"} と呼ぶ。行列では $B=S^\top JS$、同じ $u$ を関数として取る Poisson 括弧の行列は $C=S^{-1}JS^{-\top}$。$J^2=-1$ なので $BC=-1$、すなわち $C=-B^{-1}$。一自由度だけなら $B_{12}=\det S,\ C_{12}=1/\det S$ である。多自由度では成分ごとの逆数にはならない。

時間発展で運ぶ小平行四辺形の二辺を接ベクトル $v,w$ とすると、符号付き面積の和は $v^\top Jw$。相流の微分 $S$ が正準ならこの数が保たれる。滑らかな曲面を微小平行四辺形に分けて積分すれば $\int\sum_i\dd q^i\wedge\dd p_i$ が保存する。さらに Stokes の定理から境界をなす閉曲線に沿う $\oint\sum_i p_i\dd q^i$ が保存される。元の閉曲線が境界とは限らない場合にも、Hamilton 流に沿って時間微分すれば被積分量の変化は $\dd(\sum_i p_iH_{p_i}-H)$ で、閉路積分は零である。これが Poincaré の相対積分不変量である。二形式を $k$ 回外積して得る $2k$ 次形式も保たれ、最高次 $k=n$ で相体積の保存が出る。逆方向には体積保存だけから二形式の保存を結論できない。

:::example[二つの運動量に等しく交差項を足す]
基準作用をそろえた無次元座標を、位置二つ、運動量二つの順に並べる。$z=Su$ に対して

$$
S=\left(\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 2 & 1 & 0\\2 & 0 & 0 & 1\end{matrix}\right),\qquad \det S=1.
$$

微分を取る向きに注意して、座標曲線の接ベクトルから $B$ を、逆行列の行から $C$ を作ると

$$
B=S^\top JS=\left(\begin{matrix}0 & 0 & 1 & 0\\0 & 0 & 0 & 1\\-1 & 0 & 0 & 0\\0 & -1 & 0 & 0\end{matrix}\right),\qquad C=S^{-1}JS^{-\top}=\left(\begin{matrix}0 & 0 & 1 & 0\\0 & 0 & 0 & 1\\-1 & 0 & 0 & 0\\0 & -1 & 0 & 0\end{matrix}\right).
$$

$BC=-1$ は可逆性から成り立ち、この式単独では正準性を検査できない。この場合には $B=J$ なので正準である。付加した一次形式は $2\dd(q^1q^2)$ であり、微分が零になる。生成関数の対称な二階微分が行列条件に現れている。
:::

:::example[等体積でも共役平面を違う倍率で伸ばす]
基準作用をそろえた無次元座標を、位置二つ、運動量二つの順に並べる。$z=Su$ に対して

$$
S=\left(\begin{matrix}2 & 0 & 0 & 0\\0 & \frac{1}{2} & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right),\qquad \det S=1.
$$

微分を取る向きに注意して、座標曲線の接ベクトルから $B$ を、逆行列の行から $C$ を作ると

$$
B=S^\top JS=\left(\begin{matrix}0 & 0 & 2 & 0\\0 & 0 & 0 & \frac{1}{2}\\-2 & 0 & 0 & 0\\0 & - \frac{1}{2} & 0 & 0\end{matrix}\right),\qquad C=S^{-1}JS^{-\top}=\left(\begin{matrix}0 & 0 & \frac{1}{2} & 0\\0 & 0 & 0 & 2\\- \frac{1}{2} & 0 & 0 & 0\\0 & -2 & 0 & 0\end{matrix}\right).
$$

$BC=-1$ は可逆性から成り立ち、この式単独では正準性を検査できない。この場合には $B\ne J$ なので正準条件を満たさない。第一平面と第二平面の面積倍率が逆数になり、四次元体積の検査だけでは変換の失敗を見抜けない。
:::

:::example[位置の基底をせん断し、運動量に逆転置を掛ける]
基準作用をそろえた無次元座標を、位置二つ、運動量二つの順に並べる。$z=Su$ に対して

$$
S=\left(\begin{matrix}1 & 1 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & -1 & 1\end{matrix}\right),\qquad \det S=1.
$$

微分を取る向きに注意して、座標曲線の接ベクトルから $B$ を、逆行列の行から $C$ を作ると

$$
B=S^\top JS=\left(\begin{matrix}0 & 0 & 1 & 0\\0 & 0 & 0 & 1\\-1 & 0 & 0 & 0\\0 & -1 & 0 & 0\end{matrix}\right),\qquad C=S^{-1}JS^{-\top}=\left(\begin{matrix}0 & 0 & 1 & 0\\0 & 0 & 0 & 1\\-1 & 0 & 0 & 0\\0 & -1 & 0 & 0\end{matrix}\right).
$$

$BC=-1$ は可逆性から成り立ち、この式単独では正準性を検査できない。この場合には $B=J$ なので正準である。位置のせん断による交差面積が運動量の変換でちょうど取り消される。
:::

:::example[運動量の片側だけをずらす]
基準作用をそろえた無次元座標を、位置二つ、運動量二つの順に並べる。$z=Su$ に対して

$$
S=\left(\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 1 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right),\qquad \det S=1.
$$

微分を取る向きに注意して、座標曲線の接ベクトルから $B$ を、逆行列の行から $C$ を作ると

$$
B=S^\top JS=\left(\begin{matrix}0 & 1 & 1 & 0\\-1 & 0 & 0 & 1\\-1 & 0 & 0 & 0\\0 & -1 & 0 & 0\end{matrix}\right),\qquad C=S^{-1}JS^{-\top}=\left(\begin{matrix}0 & 0 & 1 & 0\\0 & 0 & 0 & 1\\-1 & 0 & 0 & -1\\0 & -1 & 1 & 0\end{matrix}\right).
$$

$BC=-1$ は可逆性から成り立ち、この式単独では正準性を検査できない。この場合には $B\ne J$ なので正準条件を満たさない。変換後に二つの旧位置方向だけで張る平行四辺形にも符号付き面積が生じる。共役平面それぞれの倍率だけを測っても十分でない。
:::

## 無限小正準変換と対称性

### 生成子から流れを作る

関数 $G(q,p,t)$ を選び、時刻を止めてパラメータ $\eps$ に沿う方程式を $\partial_\eps z=J\nabla G$ とする。
このとき任意の観測量は $\delta F=\eps\{F,G\}$ だけ変わる。
Jacobi 恒等式から $\delta\{F,K\}=\{\delta F,K\}+\{F,\delta K\}$。
従って基本括弧の値は有限変換を作る各段階でも保存される。
$G$ をこの変換の :term[生成子]{en="generator" yomi="せいせいし"} という。
逆に $F_2=\sum_i q^iP_i+\eps G(q,P,t)$ とすれば一次まで

$$
Q^i=q^i+\eps G_{p_i},\qquad P_i=p_i-\eps G_{q^i}
$$

となる。
例えば $G=p_j$ は $q^j$ の並進、$G=L_z$ は面内回転、$G=H$ は自律系の時間発展を生成する。

### Hamilton 形式での Noether の定理

$F_2$ による新 Hamiltonian を同じ引数 $Z$ で比較する。
$z=Z-\eps J\nabla G+O(\eps^2)$ と戻すため、

$$
K(Z,t)-H(Z,t)=\eps(G_t-\{H,G\})=\eps(G_t+\{G,H\})
$$

が一次で成り立つ。
従って保存量 $G$ は $H$ の形を保つ正準変換を生成する。
右辺が時刻だけの関数なら $G$ からその時間積分を引いて保存量にできる。
これが :ref[thm-am-noether] の Hamilton 表示である。
点変換の Noether 量を Legendre 変換して $G=\sum_i p_i\xi^i-H\tau-B$ とすれば、$G_{p_i}=\xi^i-\dot q^i\tau$ となり、
同時刻の変位が正しく戻る。

:::example[自由粒子の Galilei 変換]
$G=tp-mq$ とすると $\delta q=\eps t,\ \delta p=\eps m$。
明示的時間微分は $p$ だが $\{G,p^2/(2m)\}=-p$ だから $G$ は保存する。
この変換は時刻を固定したまま、その時刻に比例する位置の変位と、一定の運動量の変化を与える。
:::

### 回転座標で保存する Hamiltonian

平面で一定角速度 $\Omega$ の回転軸に乗る。
空間の位置を $q$、回転座標の位置を $Q=R(-\Omega t)q$ とし、$F_2=P^\top R(-\Omega t)q$ と置く。
$p=R(\Omega t)P$ だから速度二乗に対応する $p^\top p$ は $P^\top P$ と等しい。
平面で $\vb e_z\times Q=(-Q_y,Q_x)$ を使えば明示的時間微分は

$$
F_{2,t}=-\Omega\,P\cdot(\vb e_z\times Q)=-\Omega(Q_xP_y-Q_yP_x).
$$

実験室系のポテンシャルが回転座標で静止する $V(Q)$ なら、新 Hamiltonian は

$$
K=\frac{P^\top P}{2m}+V(Q)-\Omega L_z .
$$

従って実験室でエネルギーが変化していても $H-\Omega L_z$ は保存する。
これは時間並進と空間回転を組み合わせた Noether 量の符号を選び直したものである。

運動量を消去すると $\dot Q=P/m-\Omega\vb e_z\times Q$ であり、$\dot P=-\nabla V-\Omega\vb e_z\times P$。
従って

$$
m\ddot Q=-\nabla V-2m\Omega\vb e_z\times\dot Q
-m\Omega^2\vb e_z\times(\vb e_z\times Q).
$$

第二項は Coriolis 力、第三項は遠心力である。
例えば回転座標で静止している小球を半径 $a$ に保つには、内向きの力 $m\Omega^2a$ が必要になる。

### 有限な対称変換の読み方

$G$ が時間によらず、$D_G F=\{F,G\}$ の流れが存在するとする。
観測量を流れに沿って何度も微分して Taylor 展開すると

$$
F(\eps)=F(0)+\eps D_GF(0)+\tfrac12\eps^2D_G^2F(0)+\cdots .
$$

解析的な場合は収束域内で有限変換を再現する。
例えば $G=L_z$ なら $D_Gx=-y,\ D_G^2x=-x,\ D_G^3x=y$ なので $x(\eps)=x\cos\eps-y\sin\eps$、同様に $y(\eps)=x\sin\eps+y\cos\eps$。
$G=p_x$ なら $D_Gx=1,\ D_G^2x=0$ で $x(\eps)=x+\eps$。
同じ括弧から、一方は回転、他方は並進が得られる。

### 保存される量と生成される変換を両方向にたどる

保存量が与えられたとき、それが動かすのは配置だけとは限らない。
例えば調和振動子のエネルギーを $G=H$ として、時刻と独立にパラメータ $s$ を導入すると

$$
\frac{\partial Q}{\partial s}=P/m,\qquad
\frac{\partial P}{\partial s}=-m\omega^2Q .
$$

これは運動そのものに沿う正準変換である。
一方、位相を指定する保存量

$$
G(q,p,t)=p\cos\omega t+m\omega q\sin\omega t
$$

を使うと、$G_t=-\omega p\sin\omega t+m\omega^2q\cos\omega t$ と $\{G,H\}=\omega p\sin\omega t-m\omega^2 q\cos\omega t$ が打ち消し合う。
固定時刻での生成方向は

$$
\delta q=\epsilon\cos\omega t,\qquad
\delta p=-\epsilon m\omega\sin\omega t
$$

であり、既存の軌道に初期変位 $\epsilon$ の振動解を足す。
運動を時間に沿って送る変換と、異なる軌道に移す変換をこの二つで比較できる。

:::example[生成子の次元を合わせる]
角運動量による回転で $\delta x=\epsilon\{x,L_z\}=-\epsilon y$ となるから $\epsilon$ は角度を表す無次元数である。
一方、$G=p_x$ なら $\{x,p_x\}=1$ で $\epsilon$ は長さ、$G=H$ なら $\{x,H\}=v_x$ で $\epsilon$ は時間を持つ。
一般に括弧は二つの引数の単位の積を作用の単位で割った次元をもつ。
従って $\epsilon G$ が作用の単位をもつようにパラメータを選べば $\delta F$ と $F$ の単位が一致する。
:::

### 幾何学的な座標変換と系の対称性

正準変換を行えば同じ現象を新しい変数で記述できる。
一方、物理的な対称性として利用するには、新 Hamiltonian が同じ変数の関数としても元の形を保つ必要がある。
具体的に $G=qp$ の尺度変換を $H=p^2/(2m)+kq^2/2$ で調べると

$$
\{G,H\}=\frac{p^2}{m}-kq^2 .
$$

この量は軌道の大部分で零でない。
従って $Q=e^\epsilon q,\ P=e^{-\epsilon}p$ は正準だがこの振動子の対称性ではない。
実際、新変数では

$$
K=e^{2\epsilon}P^2/(2m)+e^{-2\epsilon}kQ^2/2 .
$$

位相楕円の面積は保たれ、軸の長さは変わる。
周期平均に限れば $\langle\{G,H\}\rangle=0$ が成り立ち、運動エネルギーとばねエネルギーの平均が等しいことが分かる。

生成子の流れでどちらの方向に変換するかによって、関数の式に出る符号が変わる場合もある。
新しい点で同じ関数を評価する能動変換では $F(z+\epsilon J\nabla G)=F(z)+\epsilon\{F,G\}$。
同じ物理量を新変数で書き直す受動変換では逆写像を代入して $F(z)=F(Z)-\epsilon\{F,G\}(Z)$。
この区別を明示すれば、時計回りと反時計回りを別々の括弧規約と誤解せずに済む。

### 微小なせん断を合成すると何が残るか

二つの簡単な生成子の流れを順番に作用させると、それぞれの変換の一次項を打ち消しても二次に正味の動きが残る。
一自由度の無次元正準対で $G=p^2/2,\ K=q^2/2$ を選ぶ。
流れの行列は

$$
S_G(s)=\begin{pmatrix}1&s\\0&1\end{pmatrix},\qquad
S_K(s)=\begin{pmatrix}1&0\\-s&1\end{pmatrix}.
$$

$G$、$K$、$G$ の逆、$K$ の逆の順に点を動かすと、全行列は右から順に掛けて

$$
S_K(-s)S_G(-s)S_K(s)S_G(s)
=\begin{pmatrix}1+s^2&s^3\\s^3&1-s^2+s^4\end{pmatrix}
$$

となる。
従って二次まででは $Q=(1+s^2)q,\ P=(1-s^2)p$。
$\{K,G\}=qp$ の生成する尺度変換をパラメータ $s^2$ だけ行ったものに一致する。
順序を逆にすれば二次の符号も逆になる。
交換子と Poisson 括弧の順序を自分で確かめるには、このように点を実際に一周させるとよい。

### 無限小の条件を有限の流れに戻す

初期値 $z$ でのベクトル場 $Y(z)$ が与えられ、その一ステップ $z\mapsto z+sY(z)$ の正準性を調べる。
$S=I+sDY$ を展開すると必要条件は

$$
(DY)^\top J+J(DY)=0 .
$$

これは一次の条件であり、そのまま有限の $s$ を代入した行列が正準になるとは限らない。
$Y=J\nabla G$ なら上式を満たすが、一般に $s^2(DY)^\top J(DY)$ の項が残る。
そこで真の流れを積分する必要がある。
例えば $G=(q^2+p^2)/2$ で一次近似を使うと $Q=q+sp,\ P=p-sq$ で行列式は $1+s^2$。
微分方程式を最後まで解けば $Q=q\cos s+p\sin s,\ P=p\cos s-q\sin s$ で行列式が厳密に $1$ となる。

:::example[生成子に定数を足す自由度]
$G'=G+c(t)$ としても、位置と運動量についての微分は同じなので $\{F,G'\}=\{F,G\}$。
従って固定時刻で同じ変換を生成する。
一方、保存条件の左辺には $\dot c(t)$ が加わる。
例えば自由粒子に対し $G=p+t f_0$ なら同じ並進を生成するが、その値は一定の割合 $f_0$ で動く。
$t f_0$ を引けば保存する生成子 $p$ が戻る。
運動方程式の対称性から生成子の数値まで一意には決まらないことが分かる。
:::

### Galilei 変換の生成子と質量

自由粒子の並進の生成子 $P_i=p_i$ と等速変換の生成子 $B_i=tp_i-mq^i$ を比べる。能動変換に対し $\{q^j,B_i\}=t\delta_{ij},\ \{p_j,B_i\}=m\delta_{ij}$。従って $B_i$ のパラメータは速度である。時間に依存するにもかかわらず $\partial_t B_i+\{B_i,H\}=p_i-p_i=0$ なので保存される。二つの異なる生成子の括弧は
$$
\{P_i,B_j\}=m\delta_{ij},\qquad \{B_i,B_j\}=0 .
$$
右辺の質量は定数であり、どの関数との括弧も零である。従って生成子の括弧は零でなくても、位相点に働く並進と等速変換の順序を交換できる。質点系なら $P_i=\sum_a p_{ai},\ B_i=tP_i-\sum_a m_a q_a^i$ として右辺に全質量が出る。これは:ref[thm-am-noether] による等速変換の境界項 $m\vb v\cdot\vb q$ の同じ情報を正準形式で読み直した結果である。

例えば $m=2\,\mathrm{kg}$ で、位置を $a=3\,\mathrm m$ 移し、速度 $b=4\,\mathrm{m/s}$ の変換をかける。一方の順序では境界関数 $mbq$ にすでに $q+a$ が入るので、他方との境界関数の差は $mab=24\,\mathrm{J\,s}$。これは位置と時刻に依存しない定数である。境界関数自体にずれがあっても、作用が受ける終端と始端での差は相殺し、観測した位置・運動量には順序の違いが現れない。

### Lie 変換で弱い非調和項を平均する

無限小正準変換の用途は既知の対称性を追うことに限らない。小さな摂動をもつ Hamiltonian の速く振動する項を新座標から取り除ける。無次元の振動子で $q=\sqrt{2I}\sin\theta,\ p=\sqrt{2I}\cos\theta$ と置くと $\dd q\wedge\dd p=\dd\theta\wedge\dd I$、$H_0=\omega I$ とする。$H=H_0+\epsilon H_1$ の旧点を新点から生成子 $W$ のパラメータ $\epsilon$ の流れで得る規約を採ると、新表示は $H+\epsilon\{H,W\}+O(\epsilon^2)$。一次項は $H_1-\omega W_\theta$ なので
$$
W(\theta,I)=\frac1{\omega}\int_0^\theta[H_1(s,I)-\overline H_1(I)]\dd s,\qquad
\overline H_1(I)=\frac1{2\pi}\int_0^{2\pi}H_1(s,I)\dd s .
\label{eq-am-lie-average}
$$
角平均を引いたため $W$ は一価な周期関数である。近似変換の後は $I$ の変動は二次の桁で、周波数は $\omega+\epsilon\,\partial_I\overline H_1$ となる。旧座標で測る振動波形には変換自体による高調波が残る。多自由度では一つの $\omega$ を組合せ周波数に置き換えるため、零に近い分母でこの手順が破れる。

:::example[ばねの四次の硬化]
$\omega=1$ とし、無次元の摂動項 $H_1=I^{2} \sin^{4}{\left(\theta \right)}$ を取る。角を一周積分して $2\pi$ で割ると $\overline H_1=\frac{3 I^{2}}{8}$。周期関数の原始関数を一つ選べば

$$
W=- \frac{I^{2} \left(8 \sin{\left(2 \theta \right)} - \sin{\left(4 \theta \right)}\right)}{32},\qquad H_1+\{H_0,W\}=\frac{3 I^{2}}{8}.
$$

右の等式は $\{H_0,W\}=-W_\theta$ を微分で直接確かめられる。よって $\partial_I\overline H_1=\frac{3 I}{4}$。作用 $I=2$、$\epsilon=0.01$ では一次までの無次元周波数は $\frac{203}{200}$ である。値を毎秒に直すには時間単位の逆数を掛ける。ばねの復元力が硬化する場合、高いエネルギーほど周波数が増す。
:::

:::example[運動量の四次項]
$\omega=1$ とし、無次元の摂動項 $H_1=- I^{2} \cos^{4}{\left(\theta \right)}$ を取る。角を一周積分して $2\pi$ で割ると $\overline H_1=- \frac{3 I^{2}}{8}$。周期関数の原始関数を一つ選べば

$$
W=- \frac{I^{2} \left(8 \sin{\left(2 \theta \right)} + \sin{\left(4 \theta \right)}\right)}{32},\qquad H_1+\{H_0,W\}=- \frac{3 I^{2}}{8}.
$$

よって $\partial_I\overline H_1=- \frac{3 I}{4}$。作用 $I=2$、$\epsilon=0.01$ では一次までの無次元周波数は $\frac{197}{200}$ である。値を毎秒に直すには時間単位の逆数を掛ける。大きな運動量をもつ時間帯でエネルギーの補正が負になり、同じ作用の周波数は低くなる。
:::

:::example[位置と運動量の混合]
$\omega=1$ とし、無次元の摂動項 $H_1=4 I^{2} \sin^{2}{\left(\theta \right)} \cos^{2}{\left(\theta \right)}$ を取る。角を一周積分して $2\pi$ で割ると $\overline H_1=\frac{I^{2}}{2}$。周期関数の原始関数を一つ選べば

$$
W=- \frac{I^{2} \sin{\left(4 \theta \right)}}{8},\qquad H_1+\{H_0,W\}=\frac{I^{2}}{2}.
$$

よって $\partial_I\overline H_1=I$。作用 $I=2$、$\epsilon=0.01$ では一次までの無次元周波数は $\frac{51}{50}$ である。値を毎秒に直すには時間単位の逆数を掛ける。軌道のどの四分円でも同符号の寄与が出る。二倍角の項が消え、角依存には四倍角だけが残る。
:::

:::example[三次の左右非対称な項]
$\omega=1$ とし、無次元の摂動項 $H_1=I^{\frac{3}{2}} \sin^{3}{\left(\theta \right)}$ を取る。角を一周積分して $2\pi$ で割ると $\overline H_1=0$。周期関数の原始関数を一つ選べば

$$
W=\frac{I^{\frac{3}{2}} \left(\cos^{2}{\left(\theta \right)} - 3\right) \cos{\left(\theta \right)}}{3},\qquad H_1+\{H_0,W\}=0.
$$

よって $\partial_I\overline H_1=0$。作用 $I=2$、$\epsilon=0.01$ では一次までの無次元周波数は $1$ である。値を毎秒に直すには時間単位の逆数を掛ける。平均は零でも波形の左右非対称は一次から現れる。周波数補正を知るには二次まで変換する必要がある。
:::

### 交互に厳密な正準変換を行う時間積分

多くの力学系では $H=T(p)+V(q)$ の全体は解けなくても、二つの部分流は厳密に解ける。$T$ だけなら運動量一定で $q\mapsto q+h\,\partial_pT$。$V$ だけなら位置一定で $p\mapsto p-h\,\partial_qV$。それぞれ有限時間の Hamilton 流なので正準である。後者を半刻み、前者を一刻み、後者を半刻みの順に作用させると
$$
p_{1/2}=p_0-\tfrac h2 V'(q_0),\quad
q_1=q_0+hT'(p_{1/2}),\quad p_1=p_{1/2}-\tfrac h2V'(q_1)
$$
となる。これは :term[Störmer–Verlet 法]{en="Störmer–Verlet method" yomi="すてるめるゔぇるれほう"} の一つの表示である。刻みの符号を逆にすれば同じ三段階が逆順の逆写像になる。従って逆写像は負の刻みの同じ方法で得られる。運動量の符号反転で物理的な時間反転も表せるのは、$T(-p)=T(p)$ の場合である。正準性と対称性だけでは精度の次数は定まらないので、展開して確認する。

$T=p^2/(2m)$ とし、$V$ が三階まで滑らかな有界領域を考える。位置は $q_1=q_0+hp_0/m-h^2 V'(q_0)/(2m)$。運動量については $V'(q_1)=V'(q_0)+hV''(q_0)p_0/m+O(h^2)$ より $p_1=p_0-hV'(q_0)-h^2 V''(q_0)p_0/(2m)+O(h^3)$。厳密解の二階までの Taylor 係数と一致するため、一段の誤差は $O(h^3)$。有限の時間区間で一段の写像の Lipschitz 定数が $1+Ch$ 以下なら、$N$ 段の誤差は各局所誤差に $(1+Ch)^j$ を掛けた和で抑えられる。$Nh$ を固定して $O(h^2)$ の大域誤差になる。

調和振動子では刻みを $\tau=\omega h$、座標を $Q=q/q_*,P=p/(m\omega q_*)$ で無次元化する。一段の行列は
$$
S=\begin{pmatrix}
1-\tau^2/2&\tau\\
-\tau(1-\tau^2/4)&1-\tau^2/2
\end{pmatrix}.
$$
$\det S=1$ は任意の刻みで成り立つ。しかし非零の刻みで、すべての初期値からの反復が有界であるには $0<|\tau|<2$ が必要である。$\tau=0$ の写像は恒等写像である。実際に $S^\top B S=B$ を満たす $B=\diag(1-\tau^2/4,1)$ がこの範囲で正定値になり、各反復点はその楕円上にある。保存される二次式は真のエネルギー $(P^2+Q^2)/2$ と $O(\tau^2)$ 違う。数値的な位相の一段の進みは $\arccos(1-\tau^2/2)=2\arcsin(|\tau|/2)$。エネルギー誤差が小さく有界でも、周波数のわずかな違いは長時間の位相に蓄積する。

:::example[小刻みで進めた振動子]
時刻刻みを $\tau=\frac{1}{2}$ とし、自然長を通過する瞬間 $(Q,P)=(0,1)$ から始める。最初の半刻みでは復元力は零、次に位置が $\tau$ 進み、最後にその位置での復元力を半刻み掛ける。行列と最初の二状態は

$$
S=\left(\begin{matrix}\frac{7}{8} & \frac{1}{2}\\- \frac{15}{32} & \frac{7}{8}\end{matrix}\right),\qquad x_1=\left(\begin{matrix}\frac{1}{2}\\\frac{7}{8}\end{matrix}\right),\qquad x_2=\left(\begin{matrix}\frac{7}{8}\\\frac{17}{32}\end{matrix}\right).
$$

第一段の後、真の無次元エネルギーは $\frac{65}{128}$、保存二次式の半分は $\frac{1}{2}$ である。一段の固有方程式は $\lambda^2-(\frac{7}{4})\lambda+1=0$。保存する楕円の二つの半径は近く、真のエネルギーは初期値に近い範囲を振動する。
:::

:::example[安定限界ちょうどの反復]
時刻刻みを $\tau=2$ とし、自然長を通過する瞬間 $(Q,P)=(0,1)$ から始める。最初の半刻みでは復元力は零、次に位置が $\tau$ 進み、最後にその位置での復元力を半刻み掛ける。行列と最初の二状態は

$$
S=\left(\begin{matrix}-1 & 2\\0 & -1\end{matrix}\right),\qquad x_1=\left(\begin{matrix}2\\-1\end{matrix}\right),\qquad x_2=\left(\begin{matrix}-4\\1\end{matrix}\right).
$$

第一段の後、真の無次元エネルギーは $\frac{5}{2}$、保存二次式の半分は $\frac{1}{2}$ である。一段の固有方程式は $\lambda^2-(-2)\lambda+1=0$。行列式は一でも固有ベクトルが一本に合流する。初期値によっては反復数に比例する成長がある。
:::

:::example[安定限界を越す反復]
時刻刻みを $\tau=3$ とし、自然長を通過する瞬間 $(Q,P)=(0,1)$ から始める。最初の半刻みでは復元力は零、次に位置が $\tau$ 進み、最後にその位置での復元力を半刻み掛ける。行列と最初の二状態は

$$
S=\left(\begin{matrix}- \frac{7}{2} & 3\\\frac{15}{4} & - \frac{7}{2}\end{matrix}\right),\qquad x_1=\left(\begin{matrix}3\\- \frac{7}{2}\end{matrix}\right),\qquad x_2=\left(\begin{matrix}-21\\\frac{47}{2}\end{matrix}\right).
$$

第一段の後、真の無次元エネルギーは $\frac{85}{8}$、保存二次式の半分は $\frac{1}{2}$ である。一段の固有方程式は $\lambda^2-(-7)\lambda+1=0$。固有値の一方の絶対値が一を超え、正準性を厳密に満たしても数値解は真の振動から離れる。
:::

図では初期位置と運動量をともに正に取った。第一の半刻みで負の復元力を受けても運動量はまだ正なので、中央の一刻みでは位置が増える。最後の半刻みでは、ばねが初めより伸びており、第一の半刻みより大きい力積を受ける。二つの半刻みをまとめて最初の位置で力を一回だけ評価すると、図の終点は変わる。

```tikz caption="分割した準ステップの順序。鉛直線分は位置を固定した力積、水平線分は運動量を固定した移動を表す。"
\draw[->] (-.3,0)--(3.8,0) node[right] {$Q$};
\draw[->] (0,-.25)--(0,3.15) node[above] {$P$};
\fill (.8,2.6) circle (1.6pt) node[above] {$z_0$};
\draw[thick,->] (.8,2.6)--(.8,2.18);
\draw[thick,->] (.8,2.18)--(2.85,2.18);
\draw[thick,->] (2.85,2.18)--(2.85,.78);
\fill (2.85,.78) circle (1.6pt) node[right] {$z_1$};
```

## 演習問題

::::exercise{level=basic}
一自由度で $\{q^2,p^3\}$ と $\{p^3,q^2\}$ を求めよ。

:::solution
前者は $(2q)(3p^2)=6qp^2$、後者は反対称性より $-6qp^2$ である。
:::
::::

::::exercise{level=basic}
$Q=aq,\ P=bp$ が正準となる条件を求めよ。

:::solution
$\{Q,P\}=ab$。
一自由度では他の自己括弧は零だから $ab=1$ が必要十分である。
:::
::::

::::exercise{level=basic}
第三型 $F_3=-pQ/a$（$a\ne0$）から変換を求めよ。

:::solution
$q=-F_{3,p}=Q/a$、$P=-F_{3,Q}=p/a$。
従って $Q=aq,\ P=p/a$ となる。
:::
::::

::::exercise{level=basic}
第四型 $F_4=pP/c$（$c\ne0$）から変換を求め、基本括弧を確認せよ。

:::solution
$q=-P/c,\ Q=p/c$ より $P=-cq$。
$\{p/c,-cq\}=1$ である。
:::
::::

::::exercise{level=advanced}
$Q=\ln(q/q_0)$（$q,q_0>0$）の点変換で、$P$ を $p$ に比例させたい。
母関数を作れ。

:::solution
$F_2=P\ln(q/q_0)$ とすれば $Q=F_{2,P}$ であり、$p=P/q$ より $P=qp$。
$\{Q,P\}=(1/q)q=1$ とも確認できる。
:::
::::

::::exercise{level=advanced}
$G=qp$ の無限小変換を有限変換まで積分せよ。

:::solution
$\partial_\eps Q=Q,\ \partial_\eps P=-P$ を初期値 $q,p$ で解くと $Q=e^\eps q,\ P=e^{-\eps}p$。
この二つの尺度の積は $1$ なので任意の $\eps$ で正準である。
:::
::::

::::exercise{level=advanced}
$L_x=yp_z-zp_y,\ L_y=zp_x-xp_z$ の括弧を計算し、保存量の閉性と結び付けよ。

:::solution
$x$ に対応する寄与と $y$ に対応する寄与は零である。
$z$ の項は $(-p_y)(-x)-y p_x=xp_y-yp_x=L_z$。
従って $L_x,L_y$ が保存する系では Jacobi 恒等式により $L_z$ も保存する。
:::
::::

::::exercise{level=advanced}
$F_2=qP+f(t)q$ を $H=p^2/(2m)-Fq$（一定の力 $F$）に適用し、$Q$ が循環する $f$ を選べ。

:::solution
$p=P+f,\ Q=q$ より $K=(P+f)^2/(2m)+(f'-F)Q$。
$f=Ft$ とすれば $K$ は $Q$ に依存せず $P=p-Ft$ が保存する。
$\dot Q=(P+Ft)/m$ を積分すれば一様加速度の解となる。
:::
::::

::::exercise{level=advanced}
周期座標 $q\pmod{2\pi}$ の円筒上で $Q=q,\ P=p+c$（$c$ 定数）を考える。
正準性を確かめ、$c\ne0$ で旧新の正準一形式の差が大域的に一価な関数の微分にならないことを示せ。

:::solution
$\{Q,P\}=1$ なので円筒の各点で正準である。
一形式の差は $p\,\dd q-P\,\dd Q=-c\,\dd q$。
$q$ をほどいた近傍では $-cq$ の微分だが、一周に沿う積分は $-2\pi c$ である。$c=0$ は恒等変換なので、以下 $c\ne0$ とする。
一価な関数の全微分なら閉路積分は零であるから、大域的な一価の原始関数は存在しない。
作用変数は基本閉路について $c$ だけずれる。
:::
::::

::::exercise{level=advanced}
一次元の正準変換 $F_2=qP+hH(q,P)$ を小さな時間幅 $h$ の近似として用いる。
$H=p^2/(2m)+V(q)$ のときの更新式とその行列式を求めよ。

:::solution
母関数の二つの偏微分から $P=p-hV'(q)$、$Q=q+hP/m$。
従って先に力による運動量更新、次に新しい運動量による位置更新を行う。
Jacobi 行列は

$$
S=\begin{pmatrix}1-h^2V''/m&h/m\\-hV''&1\end{pmatrix}
$$

だから行列式は $1$ である。
$h$ の一次までで正準方程式と一致し、非線形 $V$ でも一歩の写像は正準になる。
:::
::::

::::exercise{level=advanced}
陽的 Euler 法を無次元振動子 $\dot Q=P,\dot P=-Q$ に一段適用し、初期面積とエネルギーが何倍になるかを調べよ。刻みが小さい極限と、一定刻みで長く進めた極限を区別せよ。

:::solution
一段の行列は $\left(\begin{smallmatrix}1&h\\-h&1\end{smallmatrix}\right)$。面積倍率は行列式の $1+h^2$、新しい二乗和も $(Q+hP)^2+(P-hQ)^2=(1+h^2)(Q^2+P^2)$ となる。$N$ 段ではどちらも $(1+h^2)^N$ 倍。$Nh=T$ を固定して $h\to0$ なら対数は $N\log(1+h^2)\sim Th\to0$ で倍率は一に戻る。固定した $h\ne0$ で $N\to\infty$ とすれば無限に増え、真の保存運動を長時間たどれない。
:::
::::

## 確認テスト

::::quiz{answer="12" tol="1e-6"}
$q=2,p=1$ で $\{q^2,p^3\}$ の数値はいくらか。

:::explanation
各引数を一回ずつ微分し $6qp^2=12$ となる。
:::
::::

::::quiz
$F_2(q,P,t)$ での $K$ はどれか。

- [ ] $H-F_{2,t}$
- [x] $H+F_{2,t}$
- [ ] 常に $H$

:::explanation
$\dd F_2$ の $\dd t$ の係数が $K-H$ だから、新しい変数で表示した $H$ に $F_{2,t}$ を加える。
:::
::::

::::quiz{answer="0.25" tol="1e-6"}
$Q=4q$、$P=bp$ の線形正準変換の $b$ を求めよ。

:::explanation
$\{Q,P\}=4b=1$ から $b=1/4$ である。
:::
::::

::::quiz
時間に陽に依存する生成子 $G$ が保存する条件はどれか。

- [ ] $\{G,H\}=0$ だけ
- [x] $G_t+\{G,H\}=0$
- [ ] $G_t=\{G,H\}$

:::explanation
軌道に沿う微分は明示的な時間依存と正準方程式による変化を加えたものである。
:::
::::

::::quiz
一自由度の円筒で $Q=q,\ P=p+c$ とする。
$c\ne0$ のとき、正しいものをすべて選べ。

- [x] 基本括弧を保存する
- [x] $\oint(p\,\dd q-P\,\dd Q)$ は円筒を一周すると零でない
- [ ] どの閉曲線上でも差の積分が消える

:::explanation
$\{Q,P\}=1$ だが一形式の差は $-c\,\dd q$。
角が一周増えると積分は $-2\pi c$ である。
固定時刻の局所的な正準性と、大域的な一価の母関数の存在は同じ条件ではない。
:::
::::

::::quiz{answer="6" tol="1e-6"}
$\{q,p\}=1$ とし、$G=qp$ の無限小変換を考える。
$q=3$ のとき、$F=q^2/3$ の変化率 $\partial_\epsilon F$ はいくらか。

:::explanation
$\partial_\epsilon F=\{q^2/3,qp\}=(2q/3)q=2q^2/3$ である。
$q=3$ を代入して $6$ を得る。
:::
::::

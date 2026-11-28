---
title: "対称性と Killing ベクトル"
description: "Lie 微分から対称性を表し、測地線の保存量、最大対称性、静的性を調べる。"
notation:
  - symbol: "\\mathcal L_X"
    meaning: "流れ X による Lie 微分"
  - symbol: "\\xi^a"
    meaning: "Killing ベクトル"
  - symbol: "K_{ab}"
    meaning: "Killing テンソル"
---

## この章の目標

- 流れから Lie 微分の成分表示を得る。
- 計量の対称性を Killing 方程式で判定する。
- ベクトルやテンソルから測地線の保存量を作る。
- 最大対称時空と静的時空の条件を説明する。

## Lie 微分

ベクトル場 $X$ の局所流 $\varphi_s$ を $\partial_s\varphi_s=X\circ\varphi_s$、$\varphi_0=\id$ で定める。異なる点でのテンソルの値を流れによる引き戻しで比べる。:term[Lie 微分]{en="Lie derivative" yomi="りーびぶん"} は $\mathcal L_X T=(\dd/\dd s)|_0\varphi_s^*T$ である。ベクトルの引き戻しでは逆流の微分を使う。$\varphi_s^\mu=x^\mu+sX^\mu+\order{s^2}$ を使うと、評価点の移動は $X^\nu\partial_\nu$ を与える。さらに各下添字から $+\partial_\mu X^\nu$、各上添字から $-\partial_\mu X^\nu$ を得る。従って $\mathcal L_XY=[X,Y]$、$(\mathcal L_X\omega)_\mu=X^\nu\partial_\nu\omega_\mu+\omega_\nu\partial_\mu X^\nu$ となる。引き戻しが積と縮約、外微分を保つため、Lie 微分は積の微分則を満たし、縮約および外微分と交換する。

例えば $X=x\partial_y$、$Y=y\partial_x$ なら $[X,Y]=x\partial_x-y\partial_y$ である。また $\mathcal L_X\dd y=\dd X(y)=\dd x$ となる。

### Lie 微分と共変微分の関係

Lie 微分は流れで引き戻したテンソルの変化率である。一形式 $\alpha$ について、引き戻しの定義を任意のベクトル $Y$ に適用し $\mathcal L_X(\alpha(Y))=X(\alpha(Y))$ を使うと
$$
(\mathcal L_X\alpha)(Y)=X(\alpha(Y))-\alpha([X,Y])
$$
である。ねじれのない接続で $[X,Y]=\nabla_X Y-\nabla_Y X$ を代入すると、右辺は $(\nabla_X\alpha)(Y)+\alpha(\nabla_Y X)$ となる。多くの共変スロットについて同じ操作を行えば各スロットに $\nabla X$ が加わる。反変スロットでは符号が反対になる。例えば $Y^a$ に対して $\mathcal L_X Y^a=X^b\nabla_bY^a-Y^b\nabla_b X^a$ である。どのねじれのない接続を使って記述しても、差が相殺して同じ Lie 微分になる。

:::example[面積を変える流れと保存する流れ]
平面で $X=x\partial_x+2y\partial_y$ とすると流れは $(x,y)\mapsto(e^s x,e^{2s}y)$ である。$\dd x\wedge\dd y$ を引き戻すと $e^{3s}$ 倍なので Lie 微分は $3\dd x\wedge\dd y$ になる。これに対し $Y=-y\partial_x+x\partial_y$ の流れは回転で、Jacobian の行列式が $1$ なので面積を保存する。面積の Lie 微分は各成分の偏微分の和 $\partial_x Y^x+\partial_yY^y=0$ からも確認できる。
:::

### 微分形式の Lie 微分と流束

ベクトル $X$ を $p$ 形式 $\omega$ の第一引数に入れた $(p-1)$ 形式を $\iota_X\omega$ とする。:ref[thm-mfd-cartan-formula]{full} は $\mathcal L_X=\dd\iota_X+\iota_X\dd$ と述べる。この式を使うと、閉形式の流れによる局所的な変化は完全微分である。時間 $s$ の微分は $\varphi_s^*(\mathcal L_X\omega)$ だから、閉形式に対して
$$
\varphi_s^*\omega-\omega
=\dd\int_0^s\varphi_\sigma^*(\iota_X\omega)\dd\sigma
$$
である。従って境界のないコンパクトな有向部分多様体への積分は、その面を流れで移しても同じである。ここでベクトル場 $X$ は必ずしも計量の対称性を生成する必要はない。

:::example[発散から体積の変化を読む]
時空の体積形式 $\epsilon_g$ は最高次なので $\dd\epsilon_g=0$ である。従って $\mathcal L_X\epsilon_g=(\nabla_aX^a)\epsilon_g$ となる。平坦な時空で $X=t\partial_t+2x\partial_x-y\partial_y$ とすれば、流れの Jacobian は $e^se^{2s}e^{-s}=e^{2s}$、体積の Lie 微分は $2\epsilon_g$ である。一般の計量でも、ある領域を流れで運ぶ間の体積を $V(s)=\int_{\varphi_s(D)}\epsilon_g$ と置くと、$V'(s)$ は運んだ領域での発散の積分になる。発散が全域で一定 $a$ であれば $V(s)=e^{as}V(0)$ である。
:::

### 定数倍と関数倍の違い

Lie 微分は生成する場に関して実線形だが、その場を関数倍すると $\mathcal L_{fX}g=f\mathcal L_Xg+\dd f\otimes X^\flat+X^\flat\otimes\dd f$ となる。これは座標表示で $\partial_a(fX^b)=f\partial_aX^b+(\partial_a f)X^b$ を使えば得られる。例えば Minkowski 計量で $X=\partial_t$ は Killing だが、$f=1+x$ とすると $tx$ 成分は $-1$ になり、$fX$ は Killing ではない。曲線ごとに固有時で規格化した観測者の速度場は、元の時間並進の Killing 場と同じ軌道を持っていても等長変換を生むとは限らない。

### 流れを合成して交換子を読む

$X$ の流れの引き戻しは初めの二項が $I+s\mathcal L_X$ となる。$Y$ の流れも $I+t\mathcal L_Y$ である。共通の定義域で二つの流れと逆流を合成し、引き戻しの順序を $X,Y,-X,-Y$ とすれば、$st$ の係数は $\mathcal L_X\mathcal L_Y-\mathcal L_Y\mathcal L_X$ になる。積を順に展開すると一次項がすべて相殺し、二つの順序の積の差が残るためである。各引き戻しの順序は点を動かす写像の合成順と逆である。従って括弧の符号を確認するには、点を辿る順と作用素を掛ける順を最初に指定する。

:::example[剪断を二つの順序で行う]
平面で $X=y\partial_x$、$Y=x\partial_y$ とする。点への流れは $\varphi_s(x,y)=(x+sy,y)$、$\psi_t(x,y)=(x,y+tx)$ である。$\psi_{-t}\circ\varphi_{-s}\circ\psi_t\circ\varphi_s$ を点へ作用させると、$st$ の項は $(-x,y)st$ である。例えば第一段階の $x+sy$ を第二段階の $Y$ に入れると $y+tx+sty$ が生じ、逆の二段階で一次項だけが差し引かれる。この係数は $[X,Y]=-x\partial_x+y\partial_y$ の成分に一致する。各流れは面積を保ち、括弧の発散も $-1+1=0$ になる。
:::

## Killing ベクトルと等長変換

:::definition[Killing ベクトル]{#def-gr-killing}
$\mathcal L_\xi g=0$ を満たすベクトル場 $\xi$ を :term[Killing ベクトル]{en="Killing vector" yomi="きりんぐべくとる"} という。計量を保存する微分同相を等長変換という。
:::

$\nabla g=0$ と接続の下添字の対称性により Lie 微分の偏微分を共変微分に置き換えられる。成分で整理すると

$$
(\mathcal L_\xi g)_{ab}=\nabla_a\xi_b+\nabla_b\xi_a .
\label{eq-gr-killing}
$$

流れの合成則より $\frac{\dd}{\dd s}\varphi_s^*g=\varphi_s^*(\mathcal L_\xi g)$ なので、流れが等長変換であることとこの右辺が零であることは同値である。例えば Minkowski 平面の $x\partial_t+t\partial_x$ は $\partial_x\xi_t=-1$、$\partial_t\xi_x=1$ により Killing で、boost を生成する。Euclid 平面の $x\partial_y-y\partial_x$ も回転を生成する Killing 場だが、$x\partial_x+y\partial_y$ では $\mathcal L_Xg=2g$ となる。

### Killing 場の初期値を一点に置く

$\xi$ を Killing 場とし、その流れで測地線を写す。等長変換は Levi-Civita 接続を保つ。実際、引き戻した接続もねじれが零で同じ計量と整合するため、一意性によって元の接続に等しい。従って流れで得る測地線の族に対し、変分場 $\xi$ は Jacobi 方程式を満たす。初期値は $\xi_p$ と $\nabla_v\xi$ であり、後者は一点での反対称テンソル $\nabla_a\xi_b$ から決まる。両方が零なら正規近傍のすべての放射測地線上で $\xi=0$ である。さらに $\{\xi=0,\nabla\xi=0\}$ は閉集合で、この議論により開集合でもある。連結なら全体が零になる。従って局所自由度は最大 $n+n(n-1)/2$ 個である。

:::example[平坦時空で上限を満たす]
Cartesian 座標で $\xi^a=A^a+B^a{}_b x^b$ とし、$B_{ab}+B_{ba}=0$ とする。偏微分は定数で、その対称化が零なので Killing 方程式を満たす。四次元では四並進、三空間回転、三 boost が独立で、原点での値と反対称な微分を任意に指定できる。例えば $x\partial_t+t\partial_x$ は原点で零だが、微分が非零なので時空全体で零にはならない。
:::

### Killing 場の括弧が再び Killing になること

流れによる引き戻しをテンソルに作用させると、すべてのスカラーに対して $[\mathcal L_X,\mathcal L_Y]f=[X,Y]f$ である。ベクトルに対しては $[X,[Y,Z]]-[Y,[X,Z]]=[[X,Y],Z]$ で、作用素の交換子の Jacobi 恒等式から従う。一形式の任意のベクトルとの評価に積の法則を使えば同じ恒等式が従う。一般のテンソルは座標上でスカラー成分と基底ベクトル・一形式の積の和に展開できるため、結局 $[\mathcal L_X,\mathcal L_Y]=\mathcal L_{[X,Y]}$ である。従って $\mathcal L_Xg=\mathcal L_Yg=0$ なら $\mathcal L_{[X,Y]}g=0$ で、Killing 場全体は Lie 代数になる。

:::example[boost の括弧から回転を得る]
Minkowski 時空で $B_x=x\partial_t+t\partial_x$、$B_y=y\partial_t+t\partial_y$ を使う。$[B_x,B_y]$ は時間成分で $B_x(y)-B_y(x)=0$、$x$ 成分で $0-B_y(t)=-y$、$y$ 成分で $B_x(t)=x$ となる。従って $[B_x,B_y]=x\partial_y-y\partial_x$ である。異なる軸の boost は有限変換としても一般に可換でなく、微小な非可換性は空間回転を生む。一方、同じ軸ではベクトルが比例するので括弧は零である。
:::

### 等長性から座標と曲率を読み取る

Killing 場 $\xi$ が点 $p$ で零でなければ、その場の横断面を選び、流れの時刻と横断面上の座標を組み合わせられる。実際、横断面の接ベクトルと $\xi_p$ を並べた微分は可逆である。逆関数定理により流れの地図が局所座標を与え、$\xi=\partial_s$ となる。Killing 方程式はこの座標で $\partial_s g_{\mu\nu}=0$ と同値である。零点ではこの座標は作れないが、Killing 場自体は滑らかであり得る。球面の回転生成子が極で零になるのが一例である。

また等長変換は $g$ と $\nabla$ を保つので、:ref[def-gr-riemann]{full} の定義から $R^a{}_{bcd}$ も保つ。その結果 $\mathcal L_\xi R^a{}_{bcd}=0$ であり、曲率を縮約して作ったスカラー $I$ について $\xi(I)=0$ となる。この必要条件は Killing 方程式の候補を絞る際に使える。二次元で二つの曲率スカラー $I,J$ の勾配が点で独立なら、その点の近傍で両者を零にする Killing 場は零になる。連結領域全体の Killing 場なら、一点で値と微分がともに零であることから全域で零になる。

:::example[曲面の計量から回転を求める]
$g=\dd x^2+F(x)^2\dd y^2$、$F>0$ とする。$\xi=A(x,y)\partial_x+B(x,y)\partial_y$ と置くと、三つの Killing 方程式は $\partial_x A=0$、$\partial_y A+F^2\partial_x B=0$、$FF'A+F^2\partial_y B=0$ である。曲率 $-F''/F$ の $x$ 微分が非零の点では、曲率を保つため $A=0$ となる。その点の近傍では残りの二式から $B$ は定数になる。例えば $F=e^{x^2}$ で $x>0$ なら $K'=-8x\ne0$ なので、連結領域の Killing 場は $\partial_y$ の定数倍だけである。一方 $F$ が定数ならこの消去はできず、局所では平面の三つの場を持つ。
:::

### 等長変換の値と微分から局所の作用を決める

等長変換 $F$ が共変微分を保つとき、測地線の初期値問題の一意性から $F(\exp_p v)=\exp_{F(p)}((dF)_p v)$ である。従って同じ一点の像と同じ微分を持つ二つの等長変換は、共通の正規近傍で一致する。大域で定義されて連結な領域を考えているなら、二つの値と微分が一致する点の集合は閉であり、この関係で開でもあるから全域で一致する。これは連続な Killing 流だけに限らない。例えば空間の反射はある超曲面を点ごとに固定できるが、法線方向の微分は $-1$ であり、恒等変換とは区別できる。

:::example[球面の北極を固定する等長変換]
単位二球面で北極を固定する等長変換の微分は、その接平面の直交変換である。接平面で角 $\alpha$ の回転なら、指数曲線は初期方位角が $\alpha$ 増えるので、球面上では $\phi\mapsto\phi+\alpha$ になる。一方、接平面で一軸を反射すれば、その軸を含む子午線を固定して反対側を交換する。例えば $\phi\mapsto-\phi$ は北極と南極を固定するが、非零の角速度を持つ連続な回転の時刻 $s$ としては得られない。このように Killing 場の次元は、等長変換の恒等成分の連続自由度に対応する。
:::

### Killing 場を一点から運ぶ方程式を展開する

$A_{ab}=\nabla_a\xi_b$ とすると $A$ は反対称である。$T_{abc}=\nabla_a A_{bc}$ に対して $T_{abc}+T_{bca}=-R^d{}_{cab}\xi_d$ である。これと巡回した二式を順に加え引きすると
$$
2T_{abc}=\big[-R^d{}_{cab}+R^d{}_{abc}-R^d{}_{bca}\big]\xi_d=2R^d{}_{abc}\xi_d
$$
となる。最後は第一 Bianchi 恒等式を使った。任意の曲線の接線を $u$ とすれば $\nabla_u\xi_b=u^a A_{ab}$、$\nabla_u A_{bc}=u^a R^d{}_{abc}\xi_d$ である。これは Killing 場の候補を一点の初期データから輸送する方程式になる。ただし一般の計量では、どの経路で運んでも同じになるとは限らない。このため一点での値と反対称微分が任意に指定できるという主張には、上限に達する場合など追加の条件を要する。

## Killing ベクトルと保存量

アフィン測地線の接線 $p$ に沿う微分は

$$
\dv{(\xi_ap^a)}{\lambda}=p^ap^b\nabla_{(a}\xi_{b)}+\xi_a\nabla_p p^a=0 .
\label{eq-gr-killing-charge}
$$

逆にすべてのアフィン測地線で一定なら、初速度 $p$ は任意なので偏極によって $\nabla_{(a}\xi_{b)}=0$ を得る。定常座標の $\xi=\partial_t$ では $E=-p_t$、軸対称な $\zeta=\partial_\phi$ では $L=p_\phi$ を保存する。$g_{tt}=-1/4$ の地点で静止観測者の四元速度は $2\partial_t$ なので、$E=3$ の光の観測エネルギーは $6$ である。加速度が零でない場合は $\xi\cdot\nabla_p p$ が残る。

### 粒子の保存量から流体の保存流へ

$\nabla_aT^{ab}=0$ の対称テンソルと Killing 場 $\xi$ に対し $J^a=T^{ab}\xi_b$ とする。積を微分すると
$$
\nabla_aJ^a=(\nabla_aT^{ab})\xi_b+T^{ab}\nabla_{(a}\xi_{b)}=0
$$
である。二つの空間切片と側面に囲まれた領域で発散定理を使うと、切片上の荷の差は側面の流束で決まる。時間並進 $\xi=\partial_t$ に対するエネルギー流を $-J^a$ と取る。$\xi$ が未来時間的で、エネルギー流 $-J$ が未来因果的または零であれば、未来単位法線 $n$ を持つ切片で測る $-g(n,-J)$ は非負になる。この仮定は [エネルギー条件](/physics/general-relativity/einstein-equations/) で応力テンソルの条件として扱う。この保存には物質流の発散零と、背景で対応する対称性があるという二つの条件が必要である。

:::example[塵のエネルギーと粒子数]
静的背景の点粒子のエネルギーは $-p_a\xi^a$ であり、その粒子の静止質量が $m_s$、静止枠に対する Lorentz 因子が $\gamma$、$\sqrt{-\xi^2}=N$ なら $m_sN\gamma$ である。例えば $N=1/2$、$\gamma=3/2$ なら $3m_s/4$ である。粒子数流でその粒子を数えると寄与は $1$ であるから、荷によって重みが異なる。粒子数にエネルギーを掛けて集団の流束を得るには、この重みを各粒子の軌道に沿って用いる。
:::

:::example[時空の対称性で定めた荷と局所観測の比較]
$X$ を未来時間的な Killing 場、$N=\sqrt{-g(X,X)}>0$ とする。光の運動量を $p^a$ とすれば、保存量 $E=-p\cdot X$ と局所エネルギー $\mathcal E=-p\cdot X/N$ の違いは、観測者の固有時に合わせるための係数である。さらに観測者が静止枠に対して光の進む方向へ速さ $v$ で動く場合、Doppler 因子 $\gamma(1-v)$ を掛ける。

ここでは$N$ は時間並進のノルム、$v$ は静止枠に対する無次元速度であるとする。計算に用いる式を先に並べる。

$$
\begin{aligned}
\mathcal E_{\rm stat}/E&=\frac{1}{N}\\
\mathcal E_{\rm mov}/E&=\frac{\sqrt{\frac{1 - v}{1 + v}}}{N}
\end{aligned}
$$

次の条件を比較する。入力を区別してから計算し、定数の規格化をそろえた値を表に示す。

| $N$ | $v$ | $\mathcal E_{\rm stat}/E$ | $\mathcal E_{\rm mov}/E$ |
|---|---|---|---|
| $1$ | $0$ | $1$ | $1$ |
| $\frac{1}{2}$ | $0$ | $2$ | $2$ |
| $\frac{1}{2}$ | $\frac{3}{5}$ | $2$ | $1$ |
| $\frac{1}{2}$ | $- \frac{3}{5}$ | $2$ | $4$ |
| $\frac{1}{3}$ | $\frac{4}{5}$ | $3$ | $1$ |

例えば $N=1/2$ で進行方向に $v=3/5$ の観測者には静止枠との Doppler 因子が $1/2$ となり、二つの因子は相殺する。一方、逆向きの観測では因子が $2$ になる。差が出るのは測っている四元速度が異なるためで、同じ光軌道で Killing 荷は一定である。
:::

### 共形 Killing 場は光の荷を与える

$\mathcal L_X g=2\psi g$ とする。アフィン測地線の速度 $p$ に対して
$$
\dv{}{s} g(X,p)=\nabla_{(a}X_{b)}p^ap^b=\psi\,g(p,p)
$$
である。従って光的測地線ではこの右辺が零であり、$X$ が等長変換を生成しなくても $g(X,p)$ は保存する。時間的で単位の測地線なら右辺は $-\psi$ である。例えば Minkowski 時空で $X=x^a\partial_a$ は拡大を生成し、$g$ の Lie 微分は $2g$ である。測地線 $x=x_0+sp$ に入れると $g(X,p)=g(x_0,p)+s g(p,p)$ となり、一般式を直接確認できる。曲がった宇宙で共形時刻の並進から光の周波数の尺度を求める際にもこの関係を使える。

### 等長変換の保存量と力の項

粒子を測る計量に加え、外部の反対称場 $F=\dd A$ があり、質量 $m_s$、電荷 $q$ の粒子が $m_s\nabla_u u^a=qF^a{}_bu^b$ に従うとする。Killing 場 $\xi$ に対し力学的運動量の荷を微分すると
$$
\dv{}{\tau}(m_su_a\xi^a)=q(\iota_\xi F)_a u^a
$$
となる。もし $\mathcal L_\xi A=\dd\chi$ なら Cartan の公式から $\iota_\xi F=\dd(\chi-\iota_\xi A)$ である。従って $m_su_a\xi^a+q(\xi^a A_a-\chi)$ が一定になる。物質の力が非零でも、計量と外部場に共通する対称性から保存量を作れる。$A\to A+\dd f$ では $\chi\to\chi+\xi(f)$ と選べば括弧の値は保たれる。

:::example[静電ポテンシャル中で保存するエネルギー]
$A=-\varphi(x)\dd t$、$\xi=\partial_t$ とする。$\chi=0$ と選べ、負号を付けた保存量は $-m_su_t+q\varphi$ である。Minkowski 計量で $\varphi=-\mathcal E x$ とすれば $F_{tx}=-\mathcal E$ で、粒子が $\Delta x$ 進むと $m_su^t$ は $q\mathcal E\Delta x$ 増える。例えば初めに静止し $q\mathcal E\Delta x=m_s$ なら終点の Lorentz 因子は $2$ となる。粒子の機械的エネルギーの変化を、ゲージで記述するポテンシャルの変化が補う。
:::

### 共形 Killing の逆を光円錐で示す

滑らかな場 $X$ に対し、全ての光的測地線のアフィンな接線 $p$ で $g(X,p)$ が保存するとする。その一点での初速度には任意の光的ベクトルを選べるので、$S_{ab}=\nabla_{(a}X_{b)}$ は $S(k,k)=0$ を全ての光的 $k$ で満たす。正規直交枠で $k=e_0\pm e_i$ を入れると $S_{0i}=0$、$S_{ii}=-S_{00}$ である。さらに $e_0+(e_i+e_j)/\sqrt2$ を入れれば $S_{ij}=0$ となる。従って $S=\psi g$ で、$\psi=\nabla_aX^a/n$ である。光の保存量は計量の共形対称性を特徴付け、質量を持つ測地線での保存まで要求すれば $\psi=0$ が必要になる。

### Hamilton 形式で保存量の括弧を計算する

余接束の座標を $(x^a,p_a)$ とする。$H=g^{ab}p_ap_b/2$ に対して保存量候補を $Q_X=X^a p_a$ とし、Poisson 括弧の規約を $\{F,G\}=\partial_{x^a}F\,\partial_{p_a}G-\partial_{p_a}F\,\partial_{x^a}G$ とする。微分を代入すれば
$$
\{Q_X,H\}=-\tfrac12(\mathcal L_X g^{ab})p_ap_b,\qquad
\{Q_X,Q_Y\}=-Q_{[X,Y]}
$$
を得る。第一式の右辺は $g^{ab}$ の座標微分と、二つの上添字に対する $-\partial X$ の作用をまとめたものである。第二式の二項は $Y^a\partial_a X^b p_b-X^a\partial_a Y^b p_b$ である。従って Killing 場の線形保存量は Poisson 括弧で再び線形保存量を与える。この符号は $Q_X=X\cdot p$ と選んだ規約による。

:::example[回転の保存量のうち同時に指定できるもの]
Euclid 三次元で $L_i=(\boldsymbol x\times\boldsymbol p)_i$ とする。$L_z=xp_y-yp_x$、$L_x=yp_z-zp_y$、$L_y=zp_x-xp_z$ である。例えば $\{L_x,L_y\}=L_z$ は成分を直接微分して得られる。$\{\sum_i L_i^2,L_j\}=2\sum_i L_i\{L_i,L_j\}=0$ は対称な二次積と反対称記号の縮約なので零である。球対称な測地線の問題では、$L^2$ と一軸の $L_z$ を用いて角運動を整理できる。三つの成分は保存するが、二つずつの括弧は一般には零ではない。
:::

## Killing テンソルと隠れた対称性

対称テンソル $K_{a_1\cdots a_n}$ に $\nabla_{(b}K_{a_1\cdots a_n)}=0$ を課すと :term[Killing テンソル]{en="Killing tensor" yomi="きりんぐてんそる"} になる。測地線に沿って $K_{a_1\cdots a_n}p^{a_1}\cdots p^{a_n}$ を微分すると、$p$ の微分は零、残りは完全対称化との縮約なので零である。逆も任意の初速度について同次多項式の係数を比較すれば従う。Killing 場の対称積は線形保存量の積を生む。例えば球対称計量で $p_\theta=r^2\dot\theta$、$p_\phi=r^2\sin^2\theta\dot\phi$ とすると

$$
L^2=p_\theta^2+\frac{p_\phi^2}{\sin^2\theta}
$$

が保存する。$\dot p_\theta=p_\phi^2\cos\theta/(r^2\sin^3\theta)$ と $\dot p_\phi=0$ を使うと、$L^2$ の微分中の二つの非零項が打ち消し合う。赤道で $p_\theta=0,p_\phi=4$ なら $L^2=16$ である。一方、$p_\phi=4$ の軌道が $\theta=\pi/6$ を通るには $L^2\ge64$ が要る。線形保存量から分解できない場合には別の対称性を示す。例は [Kerr 時空](/physics/general-relativity/charged-rotating-bh/) の Carter 定数である。

### 積で得られる Killing テンソル

二つの Killing 場 $X,Y$ に対して $K_{ab}=X_{(a}Y_{b)}$ を取る。任意の点の任意の速度 $v$ を初速に測地線を取れる。測地線上で $K_{ab}v^av^b=(X_av^a)(Y_bv^b)$ は保存する。その微分は $\nabla_{(c}K_{ab)}v^cv^av^b$ であり、すべての $v$ で零である。対称多重線形形式は対角上の三次多項式で決まるので $\nabla_{(c}K_{ab)}=0$ となる。計量自身も計量整合性により Killing テンソルである。これらから作れない場合には、二次保存量が一次保存量から独立に現れる。

:::example[球面の三つの回転から角運動量の二乗を作る]
三次元 Euclid 空間の単位球面で、位置を $\vb n$、速度を $\dot{\vb n}$ とする。軸 $\vb e_i$ の回転の保存量は $L_i=(\vb e_i\times\vb n)\cdot\dot{\vb n}$ である。従って $\sum_i L_i^2=\abs{\vb n\times\dot{\vb n}}^2=\abs{\dot{\vb n}}^2$ となる。最後は $\vb n\cdot\dot{\vb n}=0,\abs{\vb n}=1$ を使った。例えば赤道で $\dot\phi=3$ なら $L_z=3$、他は零で、二乗の和は $9$ である。従ってこの球面の二次保存量は回転の一次保存量から構成できる。
:::

### 計量の係数を分離して二次保存量を作る

正定値二次元計量 $g=S(\dd x^2+\dd y^2)$、$S=A(x)+B(y)>0$ を取る。測地線の Hamilton 関数は $H=(p_x^2+p_y^2)/(2S)$ である。Hamilton 方程式から $\dot p_x=H A'/S$、$\dot x=p_x/S$ が従う。よって
$$
Q=p_x^2-2AH=\frac{B p_x^2-Ap_y^2}{S}
$$
の微分は $2p_x H A'/S-2HA'\dot x=0$ である。従って $Q$ の二次係数は反変 Killing テンソルの成分である。対応する共変テンソルは、二つの添字を計量で下ろして得る。この計算は $x,y$ のどちらも循環座標でない場合にも使える。例えば $A=1+x^2,B=2y^2$ では $Q$ は一般には $H$ の定数倍にはならない。$p_y=0$ なら $Q/H=2B$、$p_x=0$ なら $Q/H=-2A$ であり、同じ点でも値が異なるからである。

:::example[一次の対称性がない近傍で二次保存量が残る]
上の $A=1+x^2,B=2y^2$ を用いる。二次元の共形計量に対して $K=-(\partial_x^2+\partial_y^2)(\log S)/(2S)$ であり、$K=-(3+x^2-2y^2)/S^3$ となる。例えばこの公式は :eqref[eq-gr-curvature-components] に $g_{ij}=e^{2f}\delta_{ij}$、$f=\log S/2$ を代入すると、$R^x{}_{yxy}=-(f_{xx}+f_{yy})$ となることから得られる。$J=(K_x^2+K_y^2)/S$ も曲率から定まるスカラーである。点 $(1,1)$ では
$$
(\partial_x K,\partial_y K)=(1/64,5/32),\qquad
(\partial_x J,\partial_y J)=(-577/32768,-885/16384).
$$
二組の微分の行列式は $125/65536$ で非零である。従って両方を保つ Killing ベクトルはこの近傍では零となる。一方で $Q$ は上で示したように保存する。二次保存量が見つかっても、それを生成する等長流があるとは限らない。
:::

### Killing–Yano 形式を二乗する

反対称二階テンソル $Y_{ab}=-Y_{ba}$ が $\nabla_{(a}Y_{b)c}=0$ を満たすとき、これを :term[Killing–Yano テンソル]{en="Killing–Yano tensor" yomi="きりんぐやのてんそる"} と呼ぶ。アフィン測地線の速度 $p$ に対し $w_c=p^bY_{bc}$ と置くと
$$
\nabla_p w_c=p^ap^b\nabla_aY_{bc}=0
$$
である。$w$ のノルムが保存するため、$K_{ab}=Y_{ac}Y_b{}^c$ は二次の保存量を与える。どの初速度でも成立するので、Killing テンソルの逆の判定により $K$ は Killing テンソルになる。例えば向き付き二次元の計量体積形式は平行なので Killing–Yano 形式である。その二乗は正定値署名なら $g$、Lorentz 署名なら $-g$ になる。これは正規直交枠で体積形式の唯一の独立成分を $1$ として直接縮約すれば分かる。

## 最大対称空間

Killing 方程式を微分した $T_{abc}=\nabla_a\nabla_b\xi_c$ は $T_{abc}=-T_{acb}$ を満たす。交換子から $T_{abc}+T_{bca}=-R^d{}_{cab}\xi_d$ を得る。巡回した式を順に加減し、第一 Bianchi 恒等式を用いると $T_{abc}=R^d{}_{abc}\xi_d$ となる。従って曲線上で $(\xi,\nabla\xi)$ は閉じた線形常微分方程式を満たす。一点での $n$ 個の値と $n(n-1)/2$ 個の反対称微分で解が決まるので、連結領域での上限は $n(n+1)/2$ である。局所的にこの上限に達することを :term[最大対称]{en="maximal symmetry" yomi="さいだいたいしょう"} という。

上限に達すると初期値への評価写像は全射である。従って点を任意の方向に動かす局所等長流と、点を固定し任意の無限小 Lorentz 変換を与える流が存在する。等方群の恒等成分は正規直交な空間的二平面に推移的に作用するので、空間次元が二以上なら、その断面曲率はすべて等しい。二次元時空では曲率の独立成分がもともと一つなので、同じ定曲率形の結論を直接得る。差 $S_{abcd}=R_{abcd}-K(g_{ac}g_{bd}-g_{ad}g_{bc})$ の $S(v,w,v,w)$ は空間的二平面を張る開集合で零である。これは多項式だから全域で零である。まず $v$ について偏極し $S(v,w,z,w)=0$、次に $w$ について偏極し $S(v,w,z,u)+S(v,u,z,w)=0$ を得る。対の反対称性により $S(v,w,z,u)=S(v,u,w,z)$ となり、巡回すると第一 Bianchi 恒等式の三項が等しく、各項が零となる。よって曲率は定曲率形であり、点を移す等長変換によって $K$ は一定である。Minkowski 時空の並進 $4$ 個と Lorentz 変換 $6$ 個は上限を実現する。

### 定曲率から等長変換を作る

逆向きの局所命題も、指数座標内で証明できる。一定曲率 $K$ なら :eqref[eq-gr-jacobi] は、速度 $v$ に直交する分離 $J$ に $J''=-K g(v,v)J$ を与える。平行枠で $J(0)=0,J'(0)=w$ を解くと、一定係数の振動子の解で終点の成分が決まる。$w$ が $v$ に平行なら右辺が零なので $J=sw$ である。非光的な $v$ における指数写像での計量は、こうして $K$ と $g_p$ だけで決まる。光的な $v$ での値も滑らかさで従う。従って同じ署名と曲率を持つ二点の接空間を線形等長写像で結ぶと、二つの指数写像を介して局所等長変換が得られる。

この構成で初期点 $q$ と接空間の正規枠を滑らかに動かせる。点を固定して枠だけを動かす場合、元の点で零になり一次微分が任意の反対称写像になる Killing 場を得る。一方、点を動かせば任意の初期値を持てる。従って局所次元は上限に達する。ここで用いる等長写像は指数座標の範囲で定義したものである。例えば一つの空間方向を周期同一視した平坦時空は局所では Minkowski 時空だが、その空間方向と時間を混ぜる boost は一般に同一視と両立しない。

:::example[静的パッチで一定曲率を確認する]
$H>0$、$f=1-H^2r^2$ とし、$0<r<H^{-1}$ で計量 $-f\dd t^2+f^{-1}\dd r^2+r^2(\dd	heta^2+\sin^2	heta\,\dd\phi^2)$ を考える。計量成分は $r$ で変わるが、曲率成分を計量の積と比較すると一定になる。 座標の順序を $(t,r,\theta,\phi)$ と固定する。成分行列と逆行列は次の組になる。

$$
\begin{aligned}
(g_{ab})&=\left[\begin{matrix}H^{2} r^{2} - 1 & 0 & 0 & 0\\0 & \frac{1}{- H^{2} r^{2} + 1} & 0 & 0\\0 & 0 & r^{2} & 0\\0 & 0 & 0 & r^{2} \sin^{2}{\left(\theta \right)}\end{matrix}\right]\\
(g^{ab})&=\left[\begin{matrix}\frac{1}{\left(H r - 1\right) \left(H r + 1\right)} & 0 & 0 & 0\\0 & - \left(H r - 1\right) \left(H r + 1\right) & 0 & 0\\0 & 0 & \frac{1}{r^{2}} & 0\\0 & 0 & 0 & \frac{1}{r^{2} \sin^{2}{\left(\theta \right)}}\end{matrix}\right]\\
\det g&=- r^{4} \sin^{2}{\left(\theta \right)}
\end{aligned}
$$

逆行列を掛けて単位行列を得ることを確かめてから、:eqref[eq-gr-christoffel] を使う。同じ下添字の順序を二度計算せず、$b\le c$ の組を選ぶ。中央の列には、座標に依存する計量成分の微分だけを残した。最後の列は微分を実行した値であり、表に出ない成分は添字交換で得られるものを除いて零である。

| 成分 | 微分への代入 | 結果 |
|---|---|---|
| $\Gamma^{t}_{t r}$ | $\frac{1}{2 \left(H r - 1\right) \left(H r + 1\right)}\left(\partial_{r}\left(H^{2} r^{2} - 1\right)\right)$ | $\frac{H^{2} r}{\left(H r - 1\right) \left(H r + 1\right)}$ |
| $\Gamma^{r}_{t t}$ | $- \frac{\left(H r - 1\right) \left(H r + 1\right)}{2}\left(-\partial_{r}\left(H^{2} r^{2} - 1\right)\right)$ | $H^{2} r \left(H r - 1\right) \left(H r + 1\right)$ |
| $\Gamma^{r}_{r r}$ | $- \frac{\left(H r - 1\right) \left(H r + 1\right)}{2}\left(\partial_{r}\left(\frac{1}{- H^{2} r^{2} + 1}\right)+\partial_{r}\left(\frac{1}{- H^{2} r^{2} + 1}\right)-\partial_{r}\left(\frac{1}{- H^{2} r^{2} + 1}\right)\right)$ | $- \frac{H^{2} r}{\left(H r - 1\right) \left(H r + 1\right)}$ |
| $\Gamma^{r}_{\theta \theta}$ | $- \frac{\left(H r - 1\right) \left(H r + 1\right)}{2}\left(-\partial_{r}\left(r^{2}\right)\right)$ | $r \left(H r - 1\right) \left(H r + 1\right)$ |
| $\Gamma^{r}_{\phi \phi}$ | $- \frac{\left(H r - 1\right) \left(H r + 1\right)}{2}\left(-\partial_{r}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $r \left(H r - 1\right) \left(H r + 1\right) \sin^{2}{\left(\theta \right)}$ |
| $\Gamma^{\theta}_{r \theta}$ | $\frac{1}{2 r^{2}}\left(\partial_{r}\left(r^{2}\right)\right)$ | $\frac{1}{r}$ |
| $\Gamma^{\theta}_{\phi \phi}$ | $\frac{1}{2 r^{2}}\left(-\partial_{\theta}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $- \frac{\sin{\left(2 \theta \right)}}{2}$ |
| $\Gamma^{\phi}_{r \phi}$ | $\frac{1}{2 r^{2} \sin^{2}{\left(\theta \right)}}\left(\partial_{r}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{1}{r}$ |
| $\Gamma^{\phi}_{\theta \phi}$ | $\frac{1}{2 r^{2} \sin^{2}{\left(\theta \right)}}\left(\partial_{\theta}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{1}{\tan{\left(\theta \right)}}$ |

曲率の計算では :eqref[eq-gr-curvature-components] の微分差と積差を分ける。第二列は $\partial_c\Gamma^a{}_{db}$、第三列は $-\partial_d\Gamma^a{}_{cb}$ であり、第四列は接続の積を符号付きで足したものである。三列を足して最後の列を得る。

| 成分 | 第一の微分 | 第二の微分 | 接続の積 | 結果 |
|---|---|---|---|---|
| $R^{t}_{r t r}$ | $0$ | $\frac{H^{3} r}{\left(H r - 1\right) \left(H r + 1\right)^{2}} + \frac{H^{3} r}{\left(H r - 1\right)^{2} \left(H r + 1\right)} - \frac{H^{2}}{\left(H r - 1\right) \left(H r + 1\right)}$ | $- \frac{H^{4} r^{2}}{\left(H r - 1\right)^{2} \left(H r + 1\right)^{2}} - \frac{H^{4} r^{2}}{\left(H r - 1\right)^{2} \left(H r + 1\right)^{2}}$ | $- \frac{H^{2}}{\left(H r - 1\right) \left(H r + 1\right)}$ |
| $R^{t}_{\theta t \theta}$ | $0$ | $0$ | $H^{2} r^{2}$ | $H^{2} r^{2}$ |
| $R^{t}_{\phi t \phi}$ | $0$ | $0$ | $H^{2} r^{2} \sin^{2}{\left(\theta \right)}$ | $H^{2} r^{2} \sin^{2}{\left(\theta \right)}$ |
| $R^{r}_{\theta r \theta}$ | $H r \left(H r - 1\right) + H r \left(H r + 1\right) + \left(H r - 1\right) \left(H r + 1\right)$ | $0$ | $- H^{2} r^{2} - \left(H r - 1\right) \left(H r + 1\right)$ | $H^{2} r^{2}$ |
| $R^{r}_{\phi r \phi}$ | $H r \left(H r - 1\right) \sin^{2}{\left(\theta \right)} + H r \left(H r + 1\right) \sin^{2}{\left(\theta \right)} + \left(H r - 1\right) \left(H r + 1\right) \sin^{2}{\left(\theta \right)}$ | $0$ | $- H^{2} r^{2} \sin^{2}{\left(\theta \right)} - \left(H r - 1\right) \left(H r + 1\right) \sin^{2}{\left(\theta \right)}$ | $H^{2} r^{2} \sin^{2}{\left(\theta \right)}$ |
| $R^{\theta}_{\phi \theta \phi}$ | $- \cos{\left(2 \theta \right)}$ | $0$ | $\left(H r - 1\right) \left(H r + 1\right) \sin^{2}{\left(\theta \right)} + \frac{\sin{\left(2 \theta \right)}}{2 \tan{\left(\theta \right)}}$ | $H^{2} r^{2} \sin^{2}{\left(\theta \right)}$ |

第一添字を計量で下ろすと $R_{abcd}=g_{ak}R^k{}_{bcd}$ となる。対の反対称性と交換対称性を用いて、非零の独立な組を下表に示す。

| 成分 | 結果 |
|---|---|
| $R_{t r t r}$ | $- H^{2}$ |
| $R_{t \theta t \theta}$ | $H^{2} r^{2} \left(H r - 1\right) \left(H r + 1\right)$ |
| $R_{t \phi t \phi}$ | $H^{2} r^{2} \left(H r - 1\right) \left(H r + 1\right) \sin^{2}{\left(\theta \right)}$ |
| $R_{r \theta r \theta}$ | $- \frac{H^{2} r^{2}}{\left(H r - 1\right) \left(H r + 1\right)}$ |
| $R_{r \phi r \phi}$ | $- \frac{H^{2} r^{2} \sin^{2}{\left(\theta \right)}}{\left(H r - 1\right) \left(H r + 1\right)}$ |
| $R_{\theta \phi \theta \phi}$ | $H^{2} r^{4} \sin^{2}{\left(\theta \right)}$ |

Ricci テンソルは第一・第三添字を縮約して得る。曲率が残る場合も相殺が起こるので、各項を落とさずに加える。

| 成分 | 縮約の項 | 代入 | 結果 |
|---|---|---|---|
| $R_{t t}$ | $R^{r}_{t r t}+R^{\theta}_{t \theta t}+R^{\phi}_{t \phi t}$ | $H^{2} \left(H r - 1\right) \left(H r + 1\right) + H^{2} \left(H r - 1\right) \left(H r + 1\right) + H^{2} \left(H r - 1\right) \left(H r + 1\right)$ | $3 H^{2} \left(H r - 1\right) \left(H r + 1\right)$ |
| $R_{r r}$ | $R^{t}_{r t r}+R^{\theta}_{r \theta r}+R^{\phi}_{r \phi r}$ | $- \frac{H^{2}}{\left(H r - 1\right) \left(H r + 1\right)} - \frac{H^{2}}{\left(H r - 1\right) \left(H r + 1\right)} - \frac{H^{2}}{\left(H r - 1\right) \left(H r + 1\right)}$ | $- \frac{3 H^{2}}{\left(H r - 1\right) \left(H r + 1\right)}$ |
| $R_{\theta \theta}$ | $R^{t}_{\theta t \theta}+R^{r}_{\theta r \theta}+R^{\phi}_{\theta \phi \theta}$ | $H^{2} r^{2} + H^{2} r^{2} + H^{2} r^{2}$ | $3 H^{2} r^{2}$ |
| $R_{\phi \phi}$ | $R^{t}_{\phi t \phi}+R^{r}_{\phi r \phi}+R^{\theta}_{\phi \theta \phi}$ | $H^{2} r^{2} \sin^{2}{\left(\theta \right)} + H^{2} r^{2} \sin^{2}{\left(\theta \right)} + H^{2} r^{2} \sin^{2}{\left(\theta \right)}$ | $3 H^{2} r^{2} \sin^{2}{\left(\theta \right)}$ |

逆計量と縮約すると Ricci スカラーとなる。さらに $G_{ab}=R_{ab}-g_{ab}R/2$ を成分ごとに計算する。

$$
\begin{aligned}
R&=\sum_{a,b}g^{ab}R_{ab}\\
&=3 H^{2} + 3 H^{2} + 3 H^{2} + 3 H^{2}\\
&=12 H^{2}
\end{aligned}
$$

非零の Einstein テンソル成分は次のとおりである。

| 成分 | Ricci 項 | 跡の項 | 差 |
|---|---|---|---|
| $G_{t t}$ | $3 H^{2} \left(H r - 1\right) \left(H r + 1\right)$ | $6 H^{2} \left(H^{2} r^{2} - 1\right)$ | $- 3 H^{2} \left(H r - 1\right) \left(H r + 1\right)$ |
| $G_{r r}$ | $- \frac{3 H^{2}}{\left(H r - 1\right) \left(H r + 1\right)}$ | $\frac{6 H^{2}}{- H^{2} r^{2} + 1}$ | $\frac{3 H^{2}}{\left(H r - 1\right) \left(H r + 1\right)}$ |
| $G_{\theta \theta}$ | $3 H^{2} r^{2}$ | $6 H^{2} r^{2}$ | $- 3 H^{2} r^{2}$ |
| $G_{\phi \phi}$ | $3 H^{2} r^{2} \sin^{2}{\left(\theta \right)}$ | $6 H^{2} r^{2} \sin^{2}{\left(\theta \right)}$ | $- 3 H^{2} r^{2} \sin^{2}{\left(\theta \right)}$ |

例えば $R_{trtr}=H^2g_{tt}g_{rr}=-H^2$ である。角方向の $R_{\theta\phi\theta\phi}$ を $g_{\theta\theta}g_{\phi\phi}$ で割っても $H^2$ を得る。表の独立な全成分で $R_{abcd}=H^2(g_{ac}g_{bd}-g_{ad}g_{bc})$ が成立するので、このパッチは局所的に一定曲率である。
:::

### 高次元の平坦空間に埋め込んだ定曲率の模型

周囲が平坦な $n+1$ 次元空間で、位置 $X$ が $\langle X,X\rangle=\sigma L^2$（$\sigma=\pm1$）を満たす面を取る。$L>0$ とし、周囲の署名を面の誘導計量が Lorentz 署名になるように選ぶ。単位法線は $N=X/L$ で $\langle N,N\rangle=\sigma$ である。接ベクトル場に対する周囲の微分を $D$ とすると $D_U N=U/L$ だから、$\langle D_U V,N\rangle=-g(U,V)/L$ となる。従って
$$
D_U V=\nabla_U V-\sigma L^{-1}g(U,V)N.
$$
ここで $\nabla$ が Levi-Civita 接続なのは、接射影が計量整合性とねじれ零を保つことから従う。周囲の曲率が零である式の接成分を取ると
$$
R(U,V)W=\frac{\sigma}{L^2}\big[g(V,W)U-g(U,W)V\big].
$$
実際、$D_U D_V W$ の法線補正を微分した接成分は $-\sigma g(V,W)U/L^2$ で、$U,V$ を入れ替えた項を差し引けば上式になる。従って曲率は $\sigma/L^2$ である。

周囲の内積に関して反対称な線形写像 $A$ から $X\mapsto AX$ を作る。$\langle X,AX\rangle=0$ なので面に接し、$e^{sA}$ は周囲の等長変換で面を保つため Killing 場となる。こうして作る空間の次元は $n(n+1)/2$ である。これらが独立であることは、もし $AX$ が面の開パッチで消えればその一点の $X$ と、そこでの全ての接ベクトルを $A$ が零にすることから分かる。$X$ は法線で非光的だから、これらが周囲空間全体を張る。$\sigma=1$ では周囲の署名は $(-,+,\ldots,+)$、$\sigma=-1$ では時間方向を二つ取る。いずれも局所の主張であり、同一視や被覆を変えた場合には大域の等長変換がどう作用するかを別に確認する。

:::example[四次元 de Sitter 時空の静的座標]
$\sigma=1$、周囲の座標を $(X^0,X^1,X^2,X^3,X^4)$ とし、$0<r<L$ で
$$
X^0=\sqrt{L^2-r^2}\sinh(t/L),\quad
X^4=\sqrt{L^2-r^2}\cosh(t/L),\quad
(X^1,X^2,X^3)=r\boldsymbol n
$$
と置く。$\boldsymbol n$ は単位二球面の点である。最初と最後の座標の時間微分が作る計量は $-(1-r^2/L^2)\dd t^2$ で、半径の微分は $r^2/(L^2-r^2)\dd r^2$ を加える。残る三成分は $\dd r^2+r^2\dd\Omega^2$ を加えるので、合計して先の静的パッチの計量を得る。時間並進は周囲の $X^0$ と $X^4$ の boost である。例えば $r=L/2$ の静止時計は座標時刻の $\sqrt3/2$ 倍の固有時を測る。
:::

### 等方性だけと最大対称性の差

空間の各切片が定曲率でも、時空全体が最大対称になるとは限らない。例えば $g=-\dd t^2+L^2\dd\Omega_3^2$ とする。これは Minkowski の時間直線と半径 $L$ の三次元球面の直積である。時間方向には平坦なので $R_{0i0j}=0$、空間方向には $R_{ijkl}=L^{-2}(\delta_{ik}\delta_{jl}-\delta_{il}\delta_{jk})$ がある。従って断面をどちらに取るかで曲率が変わる。Ricci 作用素の固有値は時間軸で零、空間で $2/L^2$ である。等長変換はこれを保つので時間軸の直線分布を保存し、恒等成分なら未来単位時間軸 $\partial_t$ を固定する。従って任意の Killing 場は $\partial_t$ と可換でその成分は時間に依存しない。残る混合 Killing 方程式は時間成分の空間微分が零であることを要求する。このため局所の Killing 場は時間並進一つと三次元球の六つであり、四次元時空で可能な十には達しない。

:::example[周期同一視で大域の等長性が減る]
四次元の平坦時空で一つの空間座標 $x$ を長さ $\ell$ で同一視する。被覆上の任意の Killing 場は、上限の議論と $A^\mu+B^\mu{}_\nu x^\nu$ の解から、このアフィン形で尽くされる。商へ降りるには $B^\mu{}_x\ell=0$ が必要十分である。反対称性から添字 $x$ を含む全ての回転・boost が消え、時間と残り二空間の Lorentz 代数は三次元である。四つの並進は残るので、大域的な Killing 場は七次元になる。一つの座標だけ周期的にすることは局所計量を変更しないが、取り得る大域的対称性には制約を加える。
:::

## 定常時空と静的時空

時間的 Killing 場が存在する領域は :term[定常]{en="stationary" yomi="ていじょう"} である。そのような場が超曲面に直交するとき :term[静的]{en="static" yomi="せいてき"} であるという。流れに適合する座標では $\xi=\partial_t$、$\partial_tg_{\mu\nu}=0$ である。$\alpha=g(\xi,\cdot)$ を $A(\dd t+\beta_i\dd x^i)$ と書くと、$\alpha\wedge\dd\alpha=A^2(\dd t+\beta)\wedge\dd\beta$ である。係数が $t$ に依存しないので、この量が零であることと $\dd\beta=0$ は同値である。局所的に $\beta=\dd f$ と書けば $T=t+f$ 一定の面が直交する。逆に $\alpha=A\dd T$ から $\alpha\wedge\dd\alpha=0$ が従う。従ってこの三形式は、選んだ場に適合する時刻の取り方によらない局所的な判定を与える。回転座標の平坦時空に $g_{t\phi}=\Omega r^2$ があっても、$\partial_t-\Omega\partial_\phi$ は慣性系の時間並進であり、直交面を持つ。

### 一成分の交差項から直交性を判定する

計量 $g=-(\dd t+b(x)\dd y)^2+\dd x^2+\dd y^2+\dd z^2$ では $\partial_t$ はノルム $-1$ の Killing 場である。この場の一形式は $-\dd t-b(x)\dd y$ だから、その外微分との外積は $b'(x)\dd t\wedge\dd x\wedge\dd y$ である。例えば $b=2$ なら $T=t+2y$ という座標で直交する Minkowski 形に戻せる。$b(x)=x$ では同じ Killing 場に直交する面を作れない。判定の対象は選んだ Killing 場なので、時空全体に別の時間的な直交 Killing 場が存在するかどうかはさらに調べる。

一方、回転座標で $g=-\dd t^2+\dd r^2+r^2(\dd\phi+\Omega\dd t)^2+\dd z^2$ と書くと、$\partial_t$ のノルムは $-1+\Omega^2r^2$ である。$r=|\Omega|^{-1}$ で光的になるが、この面の法線 $\dd r$ のノルムは $1$ であり、零超曲面にはならない。Killing 場が光的になる面と Killing 地平線は、このように面の法線を調べて区別する。

### 直交する切片を持つかを調べる

時間的 Killing 場 $\xi=\partial_t$ に適合する座標で、計量の時間部分を $-N^2(\dd t+A_i\dd x^i)^2$ と書き、残りを空間二次形式にまとめる。全係数は $t$ に依存しない。$\xi^\flat=-N^2(\dd t+A)$ であるため
$$
\xi^\flat\wedge\dd\xi^\flat=N^4(\dd t+A)\wedge\dd A
$$
となる。ここで $\dd A$ に $\dd t$ はない。従って $\partial_t$ で縮約して $\dd A=0$ が必要であり、逆に $\dd A=0$ なら三形式は零である。可縮なパッチでは $A=\dd F$ と書けるため、時刻を $t+F$ に替えると交差成分が消える。大域的には閉形式が完全微分であるかを別に調べる。

:::example[回転座標の時間並進]
回転する Minkowski 座標で $g_{tt}=-(1-\Omega^2r^2)$、$g_{t\phi}=\Omega r^2$ とする。$|\Omega|r<1$ では $A=-\Omega r^2\dd\phi/(1-\Omega^2r^2)$ になる。微分は $-2\Omega r(1-\Omega^2r^2)^{-2}\dd r\wedge\dd\phi$ であり、$\Omega\ne0,r>0$ では $\xi=\partial_t$ に直交する切片がない。それでも Minkowski 時空自体には慣性時刻の別の静的 Killing 場がある。静的かどうかは、指定された場を調べる場合と、適切な場の存在を問う場合とで区別する。
:::

:::example[交差項のある座標で静的 Killing 場を探す]
Minkowski 時空を回転座標に取り、$\Omega$ を一定とする。$t$ に依存しないので座標場 $\partial_t$ は Killing であるが、そこから静的と判定するには、その場の直交面も調べる必要がある。 座標の順序を $(t,r,\phi,z)$ と固定する。成分行列と逆行列は次の組になる。

$$
\begin{aligned}
(g_{ab})&=\left[\begin{matrix}\Omega^{2} r^{2} - 1 & 0 & \Omega r^{2} & 0\\0 & 1 & 0 & 0\\\Omega r^{2} & 0 & r^{2} & 0\\0 & 0 & 0 & 1\end{matrix}\right]\\
(g^{ab})&=\left[\begin{matrix}-1 & 0 & \Omega & 0\\0 & 1 & 0 & 0\\\Omega & 0 & - \frac{\left(\Omega r - 1\right) \left(\Omega r + 1\right)}{r^{2}} & 0\\0 & 0 & 0 & 1\end{matrix}\right]\\
\det g&=- r^{2}
\end{aligned}
$$

逆行列を掛けて単位行列を得ることを確かめてから、:eqref[eq-gr-christoffel] を使う。同じ下添字の順序を二度計算せず、$b\le c$ の組を選ぶ。中央の列には、座標に依存する計量成分の微分だけを残した。最後の列は微分を実行した値であり、表に出ない成分は添字交換で得られるものを除いて零である。

| 成分 | 微分への代入 | 結果 |
|---|---|---|
| $\Gamma^{r}_{t t}$ | $\frac{1}{2}\left(-\partial_{r}\left(\Omega^{2} r^{2} - 1\right)\right)$ | $- \Omega^{2} r$ |
| $\Gamma^{r}_{t \phi}$ | $\frac{1}{2}\left(-\partial_{r}\left(\Omega r^{2}\right)\right)$ | $- \Omega r$ |
| $\Gamma^{r}_{\phi \phi}$ | $\frac{1}{2}\left(-\partial_{r}\left(r^{2}\right)\right)$ | $- r$ |
| $\Gamma^{\phi}_{t r}$ | $\frac{\Omega}{2}\left(\partial_{r}\left(\Omega^{2} r^{2} - 1\right)\right)- \frac{\left(\Omega r - 1\right) \left(\Omega r + 1\right)}{2 r^{2}}\left(\partial_{r}\left(\Omega r^{2}\right)\right)$ | $\frac{\Omega}{r}$ |
| $\Gamma^{\phi}_{r \phi}$ | $\frac{\Omega}{2}\left(\partial_{r}\left(\Omega r^{2}\right)\right)- \frac{\left(\Omega r - 1\right) \left(\Omega r + 1\right)}{2 r^{2}}\left(\partial_{r}\left(r^{2}\right)\right)$ | $\frac{1}{r}$ |

アフィン径数の微分を点で表し、$w_a$ を座標 $a$ の速度成分とする。ここでは添字を下げる操作を行っていない。測地線の右辺 $-\Gamma^a{}_{bc}w_bw_c$ では $b\ne c$ の項を二度数える。

$$
\begin{aligned}
\ddot{t}&=0\\
\ddot{r}&=r \left(\Omega w_{t} + w_{\phi}\right)^{2}\\
\ddot{\phi}&=- \frac{2 w_{r} \left(\Omega w_{t} + w_{\phi}\right)}{r}\\
\ddot{z}&=0
\end{aligned}
$$

また計量ノルムは

$$
g(\dot x,\dot x)=w_{\phi} \left(\Omega r^{2} w_{t} + r^{2} w_{\phi}\right) + w_{r}^{2} + w_{t} \left(\Omega r^{2} w_{\phi} + w_{t} \left(\Omega^{2} r^{2} - 1\right)\right) + w_{z}^{2}
$$
である。これを径数で微分して上の加速度を代入すると、係数の微分と速度の微分の和が零になる。従って初期のノルムが保存される。

$Y=\partial_t-\Omega\partial_\phi$ は二つの Killing 場の定数係数の線形結合であり、成分を掛けると $Y^\flat=-\dd t$、$g(Y,Y)=-1$ になる。従って $Y$ に直交する $t=\mathrm{const}$ は空間的切片である。他方 $\partial_t$ のノルムは $-1+\Omega^2r^2$ で、$r=1/|\Omega|$ で光的になり、その場の直交面は一般には存在しない。同じ交差項を持つ座標でも、どの定常流が静的性を持つかはベクトルを指定して計算する。
:::

### 定常な観測者と静的な切片のラプス

$g=-N^2(\dd t+A)^2+h_{ij}\dd x^i\dd x^j$ と書いたとき、座標に静止する観測者の速度は $u=N^{-1}\partial_t$ である。この $N=\sqrt{-g_{tt}}$ は時刻切片の単位法線が測るラプスとは一般には一致しない。時刻切片が空間的で、その切片を含む分解を
$$
g=-\mathcal N^2\dd t^2+\gamma_{ij}(\dd x^i+\beta^i\dd t)(\dd x^j+\beta^j\dd t)
$$
とする。ここで $\gamma_{ij}=g_{ij}$、$\beta^i=\gamma^{ij}g_{tj}$ である。二次形式の平方完成から $g_{tt}=-\mathcal N^2+\gamma_{ij}\beta^i\beta^j$、逆計量から $\mathcal N=(-g^{tt})^{-1/2}$ となる。$\partial_t$ が時間的であっても切片の法線とは異なる向きを向き得るので、二種類の時間の率が現れる。静的 Killing 場に直交する切片を選べば $\beta^i=0$ となり、両者は一致する。

:::example[回転する座標で二つの時計を比べる]
回転する Minkowski 計量の $t,\phi$ ブロックで $g^{tt}=-1$ だから、$t=\mathrm{const}$ 面のラプスは $\mathcal N=1$ である。シフトは $\beta^\phi=\Omega$ で、未来単位法線は $\partial_t-\Omega\partial_\phi$ である。円周上の座標静止時計は $|\Omega|r<1$ の範囲で $N=\sqrt{1-\Omega^2r^2}$ を測る。例えば $\Omega r=3/5$ なら $N=4/5$ であり、法線に沿う時計と同じ時刻間隔で比べると短く進む。この違いは Killing 場のノルムを赤方偏移の式に入れるときにも重要である。
:::

### 直交性の局所判定と大域の時刻

時間並進の一形式 $\alpha=g(\xi,\cdot)$ に $\alpha\wedge\dd\alpha=0$ を課すことは、各点を通る直交切片を局所に作れることを意味する。適合座標での $\dd A=0$ はこの条件と同じである。しかし局所の静的な時刻 $t+F$ を大域で一つに選べるかは、閉じた経路上の $\int A$ で調べる必要がある。

:::example[角変数の一周で静的な時刻がずれる]
円周 $\phi\sim\phi+2\pi$ を持つ計量
$$
g=-(\dd t+a\dd\phi)^2+R^2\dd\phi^2+\dd z^2+\dd x^2,\qquad |a|<R,
$$
を考える。$\partial_t$ は時間的 Killing 場で、$A=a\dd\phi$ は閉じている。従って局所では $T=t+a\phi$ が直交時刻である。一方、$\phi$ を一周すると $\int A=2\pi a$ であり、$a\ne0$ ならこの形の一価な大域時刻を作れない。$\xi$ に直交する曲線は $\dd t+a\dd\phi=0$ に従うので、一周の後に同じ空間ラベルへ戻っても $t$ が $-2\pi a$ だけ変わる。これはこの Killing 場とその時刻の組の大域的な性質である。
:::

### 静止観測者の加速度と Killing 場のねじれ

$\xi$ が未来時間的な Killing 場、$N=\sqrt{-\xi^2}$ とする。Killing 方程式を $\xi$ で縮約すると $\xi^a\nabla_a\xi_b=-\tfrac12\nabla_b(\xi^2)$ であり、もう一本の $\xi$ を掛ければ $\xi(N)=0$ である。従って $u=\xi/N$ の固有加速度は
$$
a_b=u^a\nabla_a u_b=\partial_b\log N
$$
である。この加速度は直交切片の存在を仮定しなくても成り立つ。直交性の方は $\alpha=\xi^\flat$ に対して $\alpha\wedge\dd\alpha$ で調べる。$\xi$ のノルムによる規格化は、$\alpha\to\alpha/N$ の微分に出る $\dd N\wedge\alpha$ が最初の $\alpha$ との外積で零になるため、この三形式の消滅条件を変えない。

観測者の空間射影を $h^a{}_b$ とし、$h_a{}^c h_b{}^d\nabla_cu_d$ を対称部と反対称部に分ける。$u=\xi/N$ を代入すると、$N$ の微分を伴う項は空間射影で消え、残りの対称部は Killing 方程式で零になる。従ってこの観測者の合同には剪断と膨張がない。ただし空間内の反対称部は一般に残る。$u^\flat\wedge\dd u^\flat$ は正規直交枠で $u^\flat=-\theta^0$ と置けば空間・空間成分だけを選ぶので、反対称部が零であることと直交切片を持つことが同値である。

:::example[ノルムが一定でも直交しない定常観測者]
$g=-(\dd t+b(x)\dd y)^2+\dd x^2+\dd y^2+\dd z^2$ で $\xi=\partial_t$ とすると、$N=1$ なので $\nabla_u u=0$ である。共枠 $\theta^0=\dd t+b\dd y$、$\theta^1=\dd x$、$\theta^2=\dd y$ を用いれば $\dd u^\flat=-b'\theta^1\wedge\theta^2$ である。例えば $b=2x$ なら空間射影した $\nabla_{[a}u_{b]}$ の $12$ 成分は $-1$ である。各観測者が自由落下していても、全ての観測者に直交する局所的な同時面があるとは限らない。
:::

### 同じ閉路で逆向きに光を送る

定常な領域で $g=-N^2(\dd t+A_i\dd x^i)^2+h_{ij}\dd x^i\dd x^j$ とする。案内路が時間に依存せず、静止する装置上の同じ空間閉曲線 $C$ を二方向に光に通らせる。未来向きの光は $-\xi\cdot k=N^2(\dot t+A_i\dot x^i)>0$ なので
$$
\dd t=\frac{\sqrt{h_{ij}\dd x^i\dd x^j}}{N}-A_i\dd x^i
$$
を選ぶ。二つの向きで第一項は同じ、第二項は符号が反転する。したがって座標時間の差 $t_C-t_{-C}=-2\oint_C A$、発光点の静止時計ではこれに $N_{\rm e}$ を掛けた値になる。時間ラベルを $T=t+F(x)$ に替えると $A$ は $A-\dd F$ に移るから、閉路での積分は同じである。閉路が小さい曲面を囲むときは $\oint A=\int\dd A$ であり、同時面の判定に使った $\dd A$ が測定にも現れる。

:::example[光の移動時刻から回転を検出する]
回転する Minkowski 座標の円周 $r=R$ では $N^2=1-\Omega^2R^2$、$A=-\Omega R^2\dd\phi/N^2$ である。$\phi$ が増える向きを正とすると $\Delta t=4\pi\Omega R^2/N^2$、発光点の静止時計の差は $\Delta\tau=4\pi\Omega R^2/N$ になる。例えば $\Omega R=3/5$ なら $N=4/5$ なので $\Delta\tau=3\pi R$ となる。この計算は静止装置が円周に光を沿わせるという仮定を使う。任意に放った光がこの円周上の零測地線になるという主張ではない。
:::

### 時間並進と軸回転を組み合わせた観測者

可換な時間並進 $\xi=\partial_t$ と軸回転 $\eta=\partial_\phi$ がある領域を取る。一定の角速度 $\Omega$ に対する $X=\xi+\Omega\eta$ は Killing であり、そのノルムは $g_{tt}+2\Omega g_{t\phi}+\Omega^2g_{\phi\phi}$ である。$X$ が未来時間的な地点では、光に対し $E=-p_t,L=p_\phi$ とすればこの観測者が測るエネルギーは
$$
\mathcal E_X=\frac{E-\Omega L}{\sqrt{-(g_{tt}+2\Omega g_{t\phi}+\Omega^2g_{\phi\phi})}}
$$
である。分子は一点での観測のために便宜的に作った組ではなく、$X$ 自身の保存量である。$\Omega=\Omega(x)$ を位置ごとに変えると観測者の速度は定義できても、この組が光路の全体で保存するとは限らない。Killing 方程式に関数の微分が加わるためである。

:::example[慣性枠の荷から回る時計の測定を得る]
Minkowski 空間の円柱座標で $g_{t\phi}=0$、$g_{\phi\phi}=r^2$ とする。$r=R$ で接線方向へ進む光に対し $L=RE$ なら、$\Omega R=3/5$ の観測者は $(1-3/5)E/(4/5)=E/2$ を測る。同じ点で逆向きの光では $L=-RE$ であり、測定値は $2E$ である。これは同じ角速度の観測者について、光が進む向きの符号だけを変えた結果である。もし $|\Omega|R=1$ まで大きくすると、$X$ が光的になりこの方法で単位四元速度を作れない。
:::

## 演習問題

::::exercise{level=basic}
$X=y\partial_x$、$Y=x\partial_y$ の括弧を求めよ。

:::solution
$X(Y^y)=y$、$Y(X^x)=x$ なので $[X,Y]=-x\partial_x+y\partial_y$ である。
:::
::::

::::exercise{level=basic}
$X=2x\partial_x+y\partial_y$ で $\mathcal L_X\dd x$ を求めよ。

:::solution
外微分との交換から $\mathcal L_X\dd x=\dd(Xx)=2\dd x$ である。
:::
::::

::::exercise{level=basic}
$\partial_t$ が Killing、$g_{tt}=-1/9$ の地点で $-p_t=2$ の光のエネルギーを静止観測者で求めよ。

:::solution
$u=3\partial_t$ なので $-u\cdot p=3(-p_t)=6$ である。
:::
::::

::::exercise{level=basic}
$g=\dd r^2+r^2\dd\phi^2$ の回転ベクトルが Killing であることを示せ。

:::solution
$\partial_\phi$ の成分は一定で、Lie 微分は $\partial_\phi g_{\mu\nu}=0$ となる。
:::
::::

::::exercise{level=advanced}
二つの Killing 場 $\xi,\zeta$ から $\xi_{(a}\zeta_{b)}$ が Killing テンソルになることを示せ。

:::solution
共変微分を積に適用する。第一因子を微分すると完全対称化で $\nabla_{(c}\xi_{a)}=0$ が現れ、第二因子も同様である。
:::
::::

::::exercise{level=advanced}
$\alpha=-(\dd t+b(x)\dd y)$ が直交面を持つ条件を求めよ。

:::solution
$\alpha\wedge\dd\alpha=b'\dd t\wedge\dd x\wedge\dd y$ なので $b'=0$ が必要十分である。このとき $t+by$ 一定の面を取れる。
:::
::::

::::exercise{level=advanced}
Killing 場の発散が零であることを示せ。

:::solution
:eqref[eq-gr-killing] の逆計量との縮約から $2\nabla_a\xi^a=0$ を得る。
:::
::::

::::exercise{level=advanced}
粒子の $(p_t,p_\phi)=(-5,2)$ とする。$X=\partial_t+\partial_\phi/2$ が $g(X,X)=-1/4$ を満たす点で、$X$ の方向の観測者が測るエネルギーを求めよ。

:::solution
$u=2X$ で規格化できるので $-u\cdot p=2(5-2/2)=8$ である。
:::
::::

::::exercise{level=advanced}
Minkowski 平面の Killing ベクトル $X=x\partial_t+t\partial_x$ の軌道を求めよ。$x>|t|$ の領域で $X$ が時間的であることを確かめ、この流れに沿った観測者が測る光のエネルギーを、保存量 $-p\cdot X$ で表せ。

:::solution
流れの方程式は $t'=x,x'=t$ なので、$t(s)=t_0\cosh s+x_0\sinh s$、$x(s)=x_0\cosh s+t_0\sinh s$ である。$\partial_t X_x=1,\partial_xX_t=-1$ で、対角微分は零だから Killing 方程式を満たす。ノルムは $-x^2+t^2$ で負、$X^t=x>0$ なので未来時間的である。観測者の単位速度は $X/\sqrt{x^2-t^2}$ であるため、$E_{\rm obs}=(-p\cdot X)/\sqrt{x^2-t^2}$ である。分子は光の測地線上で保存し、分母は観測者の軌道に沿って一定だが、別の観測者の軌道では変わる。
:::
::::

## 確認テスト

::::quiz
Killing 場について従う命題を選べ。

- [x] 局所流が計量を保つ。
- [ ] ノルムは常に一定である。
- [ ] 常に時間的である。

:::explanation
流れで引き戻した計量の微分が零である。例えば平面の回転生成子は時空へ拡張しても空間的で、ノルムは半径に依存する。
:::
::::

::::quiz{answer="10" tol="1e-9"}
四次元時空の局所 Killing 場の最大数を求めよ。

:::explanation
一点での値が $4$、一次微分の反対称部分が $6$ 成分であり、これで解が決まる。Minkowski 時空が上限 $10$ を実現する。
:::
::::

::::quiz
静的時空の直交座標で成立する条件をすべて選べ。

- [x] $\partial_tg_{\mu\nu}=0$
- [x] $g_{ti}=0$
- [ ] $\partial_ig_{tt}=0$

:::explanation
Killing 流に適合する座標では計量が時間に依存せず、直交面の空間座標との混合成分が消える。一方、時間成分の空間変化は静的条件で制約されない。
:::
::::

::::quiz{answer="16" tol="1e-9"}
赤道で $p_\theta=0,p_\phi=4$ の粒子の $L^2$ を求めよ。

:::explanation
$\sin(\pi/2)=1$ であり、$L^2=p_\theta^2+p_\phi^2/\sin^2\theta=16$ となる。
:::
::::

::::quiz{answer="12" tol="1e-7"}
$N=1/4$ の静止枠で、Killing エネルギー $E=3$ の光の局所エネルギーを求めよ。

:::explanation
静止枠の四元速度は $X/N$ である。従って $-p\cdot X/N=E/N=12$ である。
:::
::::

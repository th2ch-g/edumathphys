---
title: "曲率"
description: "平行移動の交換子から曲率を定め、縮約、共形変換、潮汐と移動枠を計算する。"
notation:
  - symbol: "R^a{}_{bcd}"
    meaning: "MTW 規約の Riemann 曲率"
  - symbol: "R_{ab}"
    meaning: "Ricci 曲率"
  - symbol: "G_{ab}"
    meaning: "Einstein テンソル"
  - symbol: "C_{abcd}"
    meaning: "Weyl テンソル"
  - symbol: "\\Omega^A{}_B"
    meaning: "曲率二形式"
---

## この章の目標

- 接続の交換子を計算し、曲率の符号を固定する。
- 曲率の対称性と Bianchi 恒等式から成分数と保存則を導く。
- 共形変換で残る Weyl 曲率と測地線偏差を求める。
- 正規直交枠で Cartan の式を使う。

## Riemann テンソルの定義と幾何学的意味

二方向に順に平行移動した結果の差を取り、微小な矩形の面積で割ると、接続係数そのものに含まれる座標依存性を取り除ける。

:::definition[Riemann 曲率]{#def-gr-riemann}
:term[Riemann 曲率テンソル]{en="Riemann curvature tensor" yomi="りーまんきょくりつてんそる"} を $R(X,Y)Z=\nabla_X\nabla_Y Z-\nabla_Y\nabla_XZ-\nabla_{[X,Y]}Z$ と定め、$R(\partial_c,\partial_d)\partial_b=R^a{}_{bcd}\partial_a$ と書く。
:::

$\nabla_X(fZ)=X(f)Z+f\nabla_XZ$ を二回用いると、$Z$ の係数の微分は $[X,Y]f$ だけ残り最後の項で消える。$X\mapsto fX$ では $-Y(f)\nabla_XZ$ が括弧の補正で消え、反対称性から $Y$ についても関数線形である。座標で展開して得る式は

$$
R^a{}_{bcd}=\partial_c\Gamma^a{}_{db}-\partial_d\Gamma^a{}_{cb}
+\Gamma^a{}_{ce}\Gamma^e{}_{db}-\Gamma^a{}_{de}\Gamma^e{}_{cb}.
\label{eq-gr-curvature-components}
$$

矩形を $+\partial_c,+\partial_d,-\partial_c,-\partial_d$ の順に回る輸送では、各辺の輸送行列は先頭で $I-\Gamma_i\Delta x^i$ である。向かい合う辺の差から $-(\partial_c\Gamma_d-\partial_d\Gamma_c)$、積の順序から $-[\Gamma_c,\Gamma_d]$ が残る。よって戻ったベクトルの変化は $-R(\partial_c,\partial_d)Z$ と向き付き面積の積を先頭項とする。平面極座標では $R^r{}_{\theta r\theta}=-1-(-r)(1/r)=0$、単位球では $R^\theta{}_{\phi\theta\phi}=(\sin^2\theta-\cos^2\theta)+\cos^2\theta=\sin^2\theta$ である。

### 一形式とテンソルに曲率を作用させる

座標基底で接続のねじれが零なら、二階共変微分の差を $[\nabla_c,\nabla_d]$ と略記できる。二階微分で下添字 $d$ に作用する接続項は、順序交換で相殺する。従ってベクトルに対してこの交換子は $R^a{}_{bcd}V^b$ である。一形式 $\omega$ についてはスカラー $\omega(V)$ に交換子を作用させると零になる。積の法則を二度用いると一階微分同士の交差項が消えるので、
$$
([\nabla_c,\nabla_d]\omega)_a=-R^b{}_{acd}\omega_b
$$
である。各添字に同じ議論を適用すると、一般のテンソルでは上添字に正、下添字に負の作用が現れる。例えば計量に適用すると $-R_{bacd}-R_{abcd}=0$ となり、最初の二添字の反対称性を別の方法でも確認できる。

:::example[Hessian の三階微分を交換する]
スカラー $f$ の Hessian は対称だが、その微分には $[\nabla_c,\nabla_d]\nabla_a f=-R^b{}_{acd}\nabla_b f$ が残る。単位球面で $f=\cos\theta$ を取る。接続を使うと $H_{\theta\theta}=-\cos\theta$、$H_{\phi\phi}=-\cos\theta\sin^2\theta$、$H_{\theta\phi}=0$ なので $H=-fg$ である。従って $\nabla_\theta H_{\phi\phi}-\nabla_\phi H_{\theta\phi}=\sin^3\theta$ となる。一形式への曲率作用では $-R^\theta{}_{\phi\theta\phi}\partial_\theta f=\sin^3\theta$ であり、下添字への作用の負号を検算できる。
:::

### 平坦な接続と局所慣性座標

開集合の全点で Levi-Civita 曲率が零なら、十分小さい座標箱に平行な枠がある。これを行列で見るため、$\nabla_i e_A=0$ を $\partial_i U=-\Gamma_i U$ と書く。まず一つの軸に沿って常微分方程式を解き、次に残る軸を一つずつ加えて解を延長する。最後に加えた軸を $x^n$ とする。以前の軸 $i<n$ に関する残差 $W_i=\partial_i U+\Gamma_iU$ は
$$
\partial_n W_i+\Gamma_nW_i
=(\partial_n\Gamma_i-\partial_i\Gamma_n+[\Gamma_n,\Gamma_i])U
$$
を満たす。右辺は零で、初期切片では $W_i=0$ だから、線形方程式の一意性により引き続き零である。枠は初期値を正規直交に選べば計量整合性から全域で正規直交になる。さらにねじれが零なので $[e_A,e_B]=0$ である。この可換な枠の局所流を順に合成すると、各流の径数を座標にできる。流が可換であるため各座標微分は対応する $e_A$ になり、この座標では計量が定数である。逆に計量が局所的に定数なら接続と曲率は零となる。

:::example[平坦性と大域的な輸送]
円錐の頂点を除いて $g=\dd r^2+\beta^2r^2\dd\phi^2$、$\phi\sim\phi+2\pi$、$\beta>0$ とする。局所変数 $\varphi=\beta\phi$ では平面極座標の形であり、曲率は零である。円周に沿って単位枠の成分を $(A,B)$ とすれば、:eqref[eq-gr-christoffel] から $A'=\beta B$、$B'=-\beta A$ となる。従って一周で $-2\pi\beta$ 回り、$2\pi$ を法として角度は $2\pi(1-\beta)$ である。例えば $\beta=3/4$ では $\pi/2$ 回転する。頂点を除いた領域の閉曲線の中には領域内で一点に縮められないものがあり、曲率が零でもその周回で輸送が変わり得る。
:::

## Riemann テンソルの対称性と独立成分の数

$R_{abcd}=g_{ae}R^e{}_{bcd}$ とする。点の正規座標では接続の積が消えるので、:eqref[eq-gr-curvature-components] に :eqref[eq-gr-christoffel] を代入して

$$
R_{abcd}=\tfrac12(g_{ad,bc}+g_{bc,ad}-g_{bd,ac}-g_{ac,bd})
$$

を得る。計量と二階微分の対称性により $R_{abcd}=-R_{bacd}=-R_{abdc}=R_{cdab}$ である。$b,c,d$ を巡回させると各二階微分が反対符号の項と対になり、$R_{abcd}+R_{acdb}+R_{adbc}=0$ となる。テンソルの等式なので任意の座標で成り立つ。

四次元では反対称な添字対が六つあり、対交換の対称性から曲率は対空間で対称な $6\times6$ 行列になる。ここから $21$ 成分を得る。巡回和は対の対称性と反対称性により四添字の全反対称化の三倍である。実際、24 項の全反対称化では三つの対の分け方が各 8 回ずつ現れる。したがって四添字がすべて異なる場合だけ独立な条件が残り、成分数は $21-1=20$ である。例えば $R_{0123}$ と $R_{0231}$ を決めると $R_{0312}=-R_{0123}-R_{0231}$ が決まる。

### 添字対の行列と四形式

$n$ 次元で曲率成分を二つの反対称な添字対にまとめると、各対は $N=n(n-1)/2$ 通りあり、対の交換対称性は $N$ 行 $N$ 列の対称行列を与える。残る第一 Bianchi 条件がどの程度独立かを判定するため、この対称行列に対応するテンソル $A_{abcd}$ について
$$
B(A)_{abcd}=\frac{A_{abcd}+A_{acdb}+A_{adbc}}3
$$
とする。対の反対称性と交換対称性からこれは完全反対称化 $A_{[abcd]}$ に等しい。実際、全 $24$ 個の並べ替えを、二対への分け方 $ab|cd,ac|db,ad|bc$ の三組に分けると、それぞれ八項が同じ符号付き成分になる。四形式は二対の反対称性と対交換対称性を持ち、$B$ で自分自身に写るので $B$ は全射である。従って曲率の独立成分数は $N(N+1)/2-\binom n4=n^2(n^2-1)/12$ となる。

:::example[四次元の対行列で一つだけ引く]
四次元では六組から成る対行列が $21$ 成分を持つ。$0123$ の係数で $B(A)=0$ を課すと $A_{0123}+A_{0231}+A_{0312}=0$ である。四つの添字のいずれかが等しい場合には完全反対称化が既に零なので、追加条件は生じない。従ってこの一つの関係で $20$ 成分になる。
:::

### 断面曲率からテンソルが決まる理由

$u,v$ が張る平面の計量が非退化なら、:term[断面曲率]{en="sectional curvature" yomi="だんめんきょくりつ"} を
$$
K(u,v)=\frac{R(u,v,u,v)}{g(u,u)g(v,v)-g(u,v)^2}
$$
とする。分子の $R(u,v,u,v)$ は下添字四つへ順に入力した値である。基底を可逆な二行二列の行列で取り替えると、反対称性から分子は行列式の二乗倍になり、分母の Gram 行列式も同じ倍率で変わる。従って値は平面だけに依存する。退化平面では分母が零であるからこの比による定義を適用できない。

曲率の対称性を持つ二つのテンソルが全ての非退化平面で同じ断面曲率を与えるとする。その差を $A$ とすると、$A(u,v,u,v)$ は非退化な組 $(u,v)$ の開集合で零となる。多項式なので全ての組で零である。まず $u$ を偏極すると $A(u,v,w,v)=0$ であり、次に $v$ を偏極すると $A(u,v,w,z)=-A(u,z,w,v)=A(u,z,v,w)$ となる。これを続けると最後の三つの入力を巡回させた各項が等しい。第一 Bianchi 恒等式でその和が零なので $A=0$ である。従って全断面の曲率を測れば全テンソルを決められる。

:::example[一つの枠の座標平面だけでは足りない]
四次元では正規直交枠の六つの座標平面から対行列の対角成分だけを得る。例えば $A_{0102}=a$ とその対称性で決まる成分だけを持つテンソルを取る。全反対称化は零なので第一 Bianchi を満たす。各座標平面の断面曲率は零だが、$e_0$ と $(e_1+e_2)/\sqrt2$ を張る平面で分子は $a$、分母は $-1$ であり、断面曲率は $-a$ となる。斜めの平面に回すと交差成分を測定できる。
:::

### 代数的な曲率を計量で実現する

曲率対称性を持つ任意の実テンソル $A_{abcd}$ は、ある点での滑らかな計量の曲率として実現できる。このことを確かめるため、原点の近くで
$$
g_{ab}(x)=\eta_{ab}-\tfrac13 A_{acbd}x^cx^d
$$
と置く。対交換で計量は対称であり、原点で $\eta$ なので近傍で非退化かつ同じ署名になる。原点での計量の二階微分は $g_{ab,cd}=-(A_{acbd}+A_{adbc})/3$ である。これを前節の四項の式に代入すると、係数 $1/6$ を除いた残りは
$$
\begin{aligned}
&-A_{abdc}-A_{acdb}-A_{bacd}-A_{bdca}\\
&\quad+A_{badc}+A_{bcda}+A_{abcd}+A_{adcb}.
\end{aligned}
$$
対の反対称性と交換対称性で整理すると $4A_{abcd}-2A_{acdb}-2A_{adbc}=6A_{abcd}$ となる。従って曲率は原点で $A$ である。点ごとの独立成分数を数えたとき、代数的な条件以外の制約が成分数をさらに減らさないこともこれで確認できる。近傍で任意の曲率を指定する問題では、さらに微分 Bianchi 恒等式が関わる。

:::example[指定した一つの潮汐成分を計量の二階項に入れる]
$A_{0101}=a$ と対称性から従う成分だけを残す。原点の時間・空間二次元のブロックは $g_{00}=-1-a(x^1)^2/3$、$g_{11}=1-a(x^0)^2/3$、$g_{01}=ax^0x^1/3$ である。原点で $R_{0101}=(2g_{01,01}-g_{11,00}-g_{00,11})/2=a$ となる。この計量の他の点まで曲率が一定と考えるには別の確認が必要である。ここでの計量は、原点の曲率が所定の値となるように作った局所模型である。
:::

## Bianchi 恒等式

前節の巡回和を第一 :term[Bianchi 恒等式]{en="Bianchi identity" yomi="びあんきこうとうしき"} という。第二恒等式は

$$
\nabla_e R^a{}_{bcd}+\nabla_c R^a{}_{bde}+\nabla_d R^a{}_{bec}=0
\label{eq-gr-bianchi}
$$

である。正規座標の原点で $\Gamma=0$ なので各共変微分は偏微分となる。曲率の接続二次項を微分しても原点で零であり、残る項は $\partial_e\partial_c\Gamma^a{}_{db}$ などの六項である。各項が逆の微分順と反対の符号で現れるので和は零である。原点を任意に選べるのでこの等式が全域で成り立つ。

例えば二つの曲面の直積 $M_1\times M_2$ の計量で各因子の曲率が一定なら、接続係数は同じ因子内だけに存在し、曲率はその因子の計量の積となる。したがって $\nabla R=0$ であり、第二恒等式の各項が個別に零となる。一般の計量では巡回和の打ち消しが必要である。

### 接続を消した点で微分 Bianchi を確認する

正規座標で考えている点では $\Gamma=0$ である。そこで $\nabla_e R^a{}_{bcd}$ は $\partial_e R^a{}_{bcd}$ に等しく、$\partial_e(\Gamma\Gamma)$ も各項に微分されない接続が残るため零である。残る項は $\partial_e\partial_c\Gamma^a{}_{db}-\partial_e\partial_d\Gamma^a{}_{cb}$ である。$(e,c,d)$ を循環させ、六項を並べると、すべて二階偏微分の交換により一対ずつ消える。共変微分した曲率の循環和はテンソルなので、点で消えれば全座標で消える。ここで $\partial R$ と $\nabla R$ を等置したのは一点だけであり、近くの点で曲率を接続の微分差だけで表したわけではない。

## Ricci テンソル・スカラー曲率・Einstein テンソル

:term[Ricci テンソル]{en="Ricci tensor" yomi="りっちてんそる"} は $R_{bd}=R^a{}_{bad}$、スカラー曲率は $R=g^{bd}R_{bd}$ とする。対を交換して縮約添字を入れ替えると $R_{bd}=R_{db}$ を得る。第二恒等式の $a,c$ を縮約すると $\nabla_e R_{bd}+\nabla_a R^a{}_{bde}-\nabla_dR_{be}=0$ となる。さらに $g^{bd}$ を掛けると、$g^{bd}R^a{}_{bde}=-R^a{}_e$ なので $\nabla_e R-2\nabla^aR_{ae}=0$ である。従って

$$
G_{ab}=R_{ab}-\tfrac12Rg_{ab},\qquad \nabla^aG_{ab}=0 .
\label{eq-gr-einstein-div}
$$

この :term[Einstein テンソル]{en="Einstein tensor" yomi="あいんしゅたいんてんそる"} の発散が零であることは、後で物質の保存則と場の方程式を両立させる。例えば $R_{abcd}=K(g_{ac}g_{bd}-g_{ad}g_{bc})$ なら四次元で $R_{ab}=3Kg_{ab}$、$R=12K$、$G_{ab}=-3Kg_{ab}$ となる。$\nabla g=0$ よりこの形では $K$ が一定である。

:::example[単位球面の成分曲率を取り出す]
単位二球面で $\theta$ を極角、$\phi$ を方位角とする。球面を含む四次元計量を扱う前に、角部分だけを計算する。 座標の順序を $(\theta,\phi)$ と固定する。成分行列と逆行列は次の組になる。

$$
\begin{aligned}
(g_{ab})&=\left[\begin{matrix}1 & 0\\0 & \sin^{2}{\left(\theta \right)}\end{matrix}\right]\\
(g^{ab})&=\left[\begin{matrix}1 & 0\\0 & \frac{1}{\sin^{2}{\left(\theta \right)}}\end{matrix}\right]\\
\det g&=\sin^{2}{\left(\theta \right)}
\end{aligned}
$$

逆行列を掛けて単位行列を得ることを確かめてから、:eqref[eq-gr-christoffel] を使う。同じ下添字の順序を二度計算せず、$b\le c$ の組を選ぶ。中央の列には、座標に依存する計量成分の微分だけを残した。最後の列は微分を実行した値であり、表に出ない成分は添字交換で得られるものを除いて零である。

| 成分 | 微分への代入 | 結果 |
|---|---|---|
| $\Gamma^{\theta}_{\phi \phi}$ | $\frac{1}{2}\left(-\partial_{\theta}\left(\sin^{2}{\left(\theta \right)}\right)\right)$ | $- \frac{\sin{\left(2 \theta \right)}}{2}$ |
| $\Gamma^{\phi}_{\theta \phi}$ | $\frac{1}{2 \sin^{2}{\left(\theta \right)}}\left(\partial_{\theta}\left(\sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{1}{\tan{\left(\theta \right)}}$ |

曲率の計算では :eqref[eq-gr-curvature-components] の微分差と積差を分ける。第二列は $\partial_c\Gamma^a{}_{db}$、第三列は $-\partial_d\Gamma^a{}_{cb}$ であり、第四列は接続の積を符号付きで足したものである。三列を足して最後の列を得る。

| 成分 | 第一の微分 | 第二の微分 | 接続の積 | 結果 |
|---|---|---|---|---|
| $R^{\theta}_{\phi \theta \phi}$ | $- \cos{\left(2 \theta \right)}$ | $0$ | $\frac{\sin{\left(2 \theta \right)}}{2 \tan{\left(\theta \right)}}$ | $\sin^{2}{\left(\theta \right)}$ |

第一添字を計量で下ろすと $R_{abcd}=g_{ak}R^k{}_{bcd}$ となる。対の反対称性と交換対称性を用いて、非零の独立な組を下表に示す。

| 成分 | 結果 |
|---|---|
| $R_{\theta \phi \theta \phi}$ | $\sin^{2}{\left(\theta \right)}$ |

Ricci テンソルは第一・第三添字を縮約して得る。曲率が残る場合も相殺が起こるので、各項を落とさずに加える。

| 成分 | 縮約の項 | 代入 | 結果 |
|---|---|---|---|
| $R_{\theta \theta}$ | $R^{\phi}_{\theta \phi \theta}$ | $1$ | $1$ |
| $R_{\phi \phi}$ | $R^{\theta}_{\phi \theta \phi}$ | $\sin^{2}{\left(\theta \right)}$ | $\sin^{2}{\left(\theta \right)}$ |

逆計量と縮約すると Ricci スカラーとなる。さらに $G_{ab}=R_{ab}-g_{ab}R/2$ を成分ごとに計算する。

$$
\begin{aligned}
R&=\sum_{a,b}g^{ab}R_{ab}\\
&=1 + 1\\
&=2
\end{aligned}
$$

Einstein テンソルはすべての成分で零となる。

$R_{\theta\phi\theta\phi}/(g_{\theta\theta}g_{\phi\phi}-g_{\theta\phi}^2)=1$ である。半径 $L$ の球面では計量が $L^2$ 倍となり、上添字を持つ曲率は同じである。従って分子は $L^2$ 倍、分母は $L^4$ 倍なので断面曲率は $1/L^2$ となる。
:::

:::example[Ricci が零で Weyl が残る光的曲率]
$\dd s^2=-2\dd u\,\dd v+a(x^2-y^2)\dd u^2+\dd x^2+\dd y^2$ を用いる。ここで $a$ は定数である。$u,v$ 部分の行列式が $-1$ であり、空間二方向が正なので Lorentz 署名となる。 座標の順序を $(u,v,x,y)$ と固定する。成分行列と逆行列は次の組になる。

$$
\begin{aligned}
(g_{ab})&=\left[\begin{matrix}a \left(x^{2} - y^{2}\right) & -1 & 0 & 0\\-1 & 0 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]\\
(g^{ab})&=\left[\begin{matrix}0 & -1 & 0 & 0\\-1 & - a \left(x - y\right) \left(x + y\right) & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]\\
\det g&=-1
\end{aligned}
$$

逆行列を掛けて単位行列を得ることを確かめてから、:eqref[eq-gr-christoffel] を使う。同じ下添字の順序を二度計算せず、$b\le c$ の組を選ぶ。中央の列には、座標に依存する計量成分の微分だけを残した。最後の列は微分を実行した値であり、表に出ない成分は添字交換で得られるものを除いて零である。

| 成分 | 微分への代入 | 結果 |
|---|---|---|
| $\Gamma^{v}_{u x}$ | $- \frac{1}{2}\left(\partial_{x}\left(a \left(x^{2} - y^{2}\right)\right)\right)$ | $- a x$ |
| $\Gamma^{v}_{u y}$ | $- \frac{1}{2}\left(\partial_{y}\left(a \left(x^{2} - y^{2}\right)\right)\right)$ | $a y$ |
| $\Gamma^{x}_{u u}$ | $\frac{1}{2}\left(-\partial_{x}\left(a \left(x^{2} - y^{2}\right)\right)\right)$ | $- a x$ |
| $\Gamma^{y}_{u u}$ | $\frac{1}{2}\left(-\partial_{y}\left(a \left(x^{2} - y^{2}\right)\right)\right)$ | $a y$ |

曲率の計算では :eqref[eq-gr-curvature-components] の微分差と積差を分ける。第二列は $\partial_c\Gamma^a{}_{db}$、第三列は $-\partial_d\Gamma^a{}_{cb}$ であり、第四列は接続の積を符号付きで足したものである。三列を足して最後の列を得る。

| 成分 | 第一の微分 | 第二の微分 | 接続の積 | 結果 |
|---|---|---|---|---|
| $R^{v}_{x u x}$ | $0$ | $a$ | $0$ | $a$ |
| $R^{v}_{y u y}$ | $0$ | $- a$ | $0$ | $- a$ |

第一添字を計量で下ろすと $R_{abcd}=g_{ak}R^k{}_{bcd}$ となる。対の反対称性と交換対称性を用いて、非零の独立な組を下表に示す。

| 成分 | 結果 |
|---|---|
| $R_{u x u x}$ | $- a$ |
| $R_{u y u y}$ | $a$ |

Ricci テンソルは第一・第三添字を縮約して得る。曲率が残る場合も相殺が起こるので、各項を落とさずに加える。

| 成分 | 縮約の項 | 代入 | 結果 |
|---|---|---|---|
| $R_{u u}$ | $R^{x}_{u x u}+R^{y}_{u y u}$ | $- a + a$ | $0$ |

逆計量と縮約すると Ricci スカラーとなる。さらに $G_{ab}=R_{ab}-g_{ab}R/2$ を成分ごとに計算する。

$$
\begin{aligned}
R&=\sum_{a,b}g^{ab}R_{ab}\\
&=0\\
&=0
\end{aligned}
$$

Einstein テンソルはすべての成分で零となる。

Ricci が零なので、後で定義する Weyl テンソルはこの例では Riemann テンソルに一致する。しかも曲率二乗のスカラーは零である。各非零の下添字成分には $u$ が二つ現れ、$g^{uu}=0$ だから、縮約でそれらを相手に結ぶには相手の $v$ 成分が必要である。そのような曲率成分は表にない。従ってスカラーが零でも曲率テンソルが零とは限らない。
:::

### 平面の曲率を Ricci の和として読む

計量の正規直交基底 $e_a$ の符号を $\epsilon_a=g(e_a,e_a)$ とする。$R_{ab}=\sum_c\epsilon_c R_{cacb}$ だから、$a$ 方向の Ricci 成分は $a,c$ 平面の断面曲率 $K_{ca}=R_{caca}/(\epsilon_c\epsilon_a)$ を用いて $R_{aa}=\epsilon_a\sum_{c\ne a}K_{ca}$ と書ける。跡はさらに各平面を二度数えるため $R=2\sum_{a<c}K_{ac}$ である。Lorentz 時空でこの表示を用いるときは、平面が退化する零方向を基底軸として選ばない。例えば空間が二球面の半径 $L$ と直線の積で、時間方向も平坦な直積なら、球面内の一平面だけが $1/L^2$ なので $R=2/L^2$ である。その一平面に属さない時間・直線方向の Ricci 成分は零であり、Einstein テンソルでは跡の項によってそれらの方向にも成分が生じる。

:::example[曲がる二平面と平坦な二平面の直積]
$\dd s^2=-N(z)^2\dd t^2+\dd z^2+\dd x^2+\dd y^2$ を用いる。$N>0$ で、時間・鉛直平面だけに曲がりを持たせる。ほかの二方向に Riemann 成分がなくても、Einstein テンソルの成分がどこに現れるかを調べる。 座標の順序を $(t,z,x,y)$ と固定する。成分行列と逆行列は次の組になる。

$$
\begin{aligned}
(g_{ab})&=\left[\begin{matrix}- N^{2}{\left(z \right)} & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]\\
(g^{ab})&=\left[\begin{matrix}- \frac{1}{N^{2}{\left(z \right)}} & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]\\
\det g&=- N^{2}{\left(z \right)}
\end{aligned}
$$

逆行列を掛けて単位行列を得ることを確かめてから、:eqref[eq-gr-christoffel] を使う。同じ下添字の順序を二度計算せず、$b\le c$ の組を選ぶ。中央の列には、座標に依存する計量成分の微分だけを残した。最後の列は微分を実行した値であり、表に出ない成分は添字交換で得られるものを除いて零である。

| 成分 | 微分への代入 | 結果 |
|---|---|---|
| $\Gamma^{t}_{t z}$ | $- \frac{1}{2 N^{2}{\left(z \right)}}\left(\partial_{z}\left(- N^{2}{\left(z \right)}\right)\right)$ | $\frac{N'{\left(z \right)}}{N{\left(z \right)}}$ |
| $\Gamma^{z}_{t t}$ | $\frac{1}{2}\left(-\partial_{z}\left(- N^{2}{\left(z \right)}\right)\right)$ | $N{\left(z \right)} N'{\left(z \right)}$ |

曲率の計算では :eqref[eq-gr-curvature-components] の微分差と積差を分ける。第二列は $\partial_c\Gamma^a{}_{db}$、第三列は $-\partial_d\Gamma^a{}_{cb}$ であり、第四列は接続の積を符号付きで足したものである。三列を足して最後の列を得る。

| 成分 | 第一の微分 | 第二の微分 | 接続の積 | 結果 |
|---|---|---|---|---|
| $R^{t}_{z t z}$ | $0$ | $- \frac{N''{\left(z \right)}}{N{\left(z \right)}} + \frac{\left(N'{\left(z \right)}\right)^{2}}{N^{2}{\left(z \right)}}$ | $- \frac{\left(N'{\left(z \right)}\right)^{2}}{N^{2}{\left(z \right)}}$ | $- \frac{N''{\left(z \right)}}{N{\left(z \right)}}$ |

第一添字を計量で下ろすと $R_{abcd}=g_{ak}R^k{}_{bcd}$ となる。対の反対称性と交換対称性を用いて、非零の独立な組を下表に示す。

| 成分 | 結果 |
|---|---|
| $R_{t z t z}$ | $N{\left(z \right)} N''{\left(z \right)}$ |

Ricci テンソルは第一・第三添字を縮約して得る。曲率が残る場合も相殺が起こるので、各項を落とさずに加える。

| 成分 | 縮約の項 | 代入 | 結果 |
|---|---|---|---|
| $R_{t t}$ | $R^{z}_{t z t}$ | $N{\left(z \right)} N''{\left(z \right)}$ | $N{\left(z \right)} N''{\left(z \right)}$ |
| $R_{z z}$ | $R^{t}_{z t z}$ | $- \frac{N''{\left(z \right)}}{N{\left(z \right)}}$ | $- \frac{N''{\left(z \right)}}{N{\left(z \right)}}$ |

逆計量と縮約すると Ricci スカラーとなる。さらに $G_{ab}=R_{ab}-g_{ab}R/2$ を成分ごとに計算する。

$$
\begin{aligned}
R&=\sum_{a,b}g^{ab}R_{ab}\\
&=- \frac{N''{\left(z \right)}}{N{\left(z \right)}} - \frac{N''{\left(z \right)}}{N{\left(z \right)}}\\
&=- \frac{2 N''{\left(z \right)}}{N{\left(z \right)}}
\end{aligned}
$$

非零の Einstein テンソル成分は次のとおりである。

| 成分 | Ricci 項 | 跡の項 | 差 |
|---|---|---|---|
| $G_{x x}$ | $0$ | $- \frac{N''{\left(z \right)}}{N{\left(z \right)}}$ | $\frac{N''{\left(z \right)}}{N{\left(z \right)}}$ |
| $G_{y y}$ | $0$ | $- \frac{N''{\left(z \right)}}{N{\left(z \right)}}$ | $\frac{N''{\left(z \right)}}{N{\left(z \right)}}$ |

$N=1+az$ では二階微分が零になり平坦である。$N=e^{az}$ では $R=-2a^2$ であり、$G_{xx}=G_{yy}=a^2$ となる。時間・鉛直方向の Einstein 成分が相殺するのは二次元の Ricci が計量に比例するからで、横方向には同じ Ricci 項がないため跡の補正が残る。
:::

:::example[一定曲率の係数を決める]
$n>2$ で $R_{abcd}=K(x)(g_{ac}g_{bd}-g_{ad}g_{bc})$ と仮定する。縮約により $R_{ab}=(n-1)Kg_{ab},R=n(n-1)K$ となる。従って $\nabla^aR_{ab}-(1/2)\partial_bR=(n-1)(1-n/2)\partial_bK$ である。係数が非零なので $K$ は連結成分上で一定である。一方 $n=2$ ではこの関係は常に零で、例えば曲面の Gauss 曲率には位置依存が許される。
:::

### 曲率の指標を選ぶ

スカラー曲率 $R$ は Ricci の跡であり、真空のように Ricci が消える場合に全曲率の大きさを測れない。この場合に使われる縮約の一つが $R_{abcd}R^{abcd}$ である。ただし Lorentz 署名ではこの縮約にも負の寄与がある。正規直交枠で縮約すると各成分の二乗に添字ごとの署名の積が掛かる。時間添字が奇数個の成分には負号が付く。従って曲率を消す判定には、必要なら全テンソルの成分または測地線偏差を調べる。

:::example[二つの曲面が作る同じスカラーと異なる曲率]
時間・空間二次元面 $Q$ と空間二次元面 $S$ の直積を取る。各因子の断面曲率を $K_Q,K_S$ とし、正規直交枠を二つの因子に合わせる。混合成分は零で、独立な成分は $R_{0101}=-K_Q$、$R_{2323}=K_S$ である。スカラーは $2(K_Q+K_S)$、曲率二乗は $4(K_Q^2+K_S^2)$ となる。後者の係数 $4$ は各添字対で二つの順序を数えたものである。例えば $K_Q=-1/L^2,K_S=1/L^2$ なら $R=0$ だが曲率二乗は $8/L^4$ である。跡の相殺と点ごとの平坦性を別々に検査できる。
:::

## Weyl テンソルと共形変換

### 跡を除く分解

四次元では :term[Schouten テンソル]{en="Schouten tensor" yomi="しゃうてんてんそる"} を $P_{ab}=(R_{ab}-Rg_{ab}/6)/2$ と置く。$B(P)_{abcd}=g_{ac}P_{bd}-g_{ad}P_{bc}-g_{bc}P_{ad}+g_{bd}P_{ac}$ とすれば、$g^{ac}B(P)_{abcd}=2P_{bd}+g_{bd}\tr_gP=R_{bd}$ である。したがって :term[Weyl テンソル]{en="Weyl tensor" yomi="わいるてんそる"} $C=Riem-B(P)$ は跡が零である。任意の対称二階テンソルを $P$ に選べば $B(P)$ で対応する Ricci 部分が実現されるので、四次元での独立成分は $20-10=10$ である。

### 計量の局所倍率

$\widetilde g=e^{2f}g$、$f_a=\nabla_af$ とする。:eqref[eq-gr-christoffel] より接続差は $D^a{}_{bc}=\delta^a_b f_c+\delta^a_c f_b-g_{bc}f^a$ である。曲率の差を展開すると $\nabla_cD^a{}_{db}-\nabla_dD^a{}_{cb}+D^a{}_{ce}D^e{}_{db}-D^a{}_{de}D^e{}_{cb}$ となる。第一添字を元の計量で下ろすと、微分部分は $B(-\nabla\nabla f)$ である。積部分では $\delta$ 同士から $g_{ac}f_bf_d-g_{ad}f_bf_c$、$g_{db}f^e$ との交差から $-g_{bc}f_af_d+g_{bd}f_af_c$ が生じ、二つの $f$ が縮約された項は $-(g_{ac}g_{bd}-g_{ad}g_{bc})f_ef^e$ となる。このため $\widetilde R_{abcd}=e^{2f}(R_{abcd}+B(S)_{abcd})$、$S=-\nabla\nabla f+\dd f\otimes\dd f-\tfrac12(\nabla f)^2 g$ である。縮約により $\widetilde P=P+S$ を得るので $\widetilde C_{abcd}=e^{2f}C_{abcd}$、$\widetilde C^a{}_{bcd}=C^a{}_{bcd}$ となる。

例えば $\widetilde g=a(\eta)^2\eta_{ab}$ では $C=0$ だが、$a=e^{H\eta}$ では $P_{\eta\eta}=H^2/2$、$P_{ij}=H^2\delta_{ij}/2$ となる。Weyl 曲率が零でも Ricci 曲率は残り得る。

### 低次元の Weyl テンソルと曲率の復元

三次元では対行列の成分数は $6$ であり、対称な Ricci テンソルの成分数も $6$ である。この数の一致だけでは復元できるとは限らないため、任意の対称テンソル $S$ に対し $B(S)_{abcd}=g_{ac}S_{bd}-g_{ad}S_{bc}-g_{bc}S_{ad}+g_{bd}S_{ac}$ の縮約を取る。結果は $(n-2)S_{bd}+(\tr S)g_{bd}$ なので、$n\ge3$ では任意の Ricci 候補に対して $S$ を一意に解ける。特に三次元では Ricci への縮約が全射かつ同次元であり、核は零である。このため Weyl テンソルは恒等的に零になる。四次元では同じ縮約が全射であっても核の次元が $20-10=10$ である。

:::example[共形に平坦でも Ricci は残る]
四次元で $g=a(\eta)^2\eta_{ab}$ とすると、上の共形変換則から Weyl は零である。$a(\eta)=-1/(H\eta)$、$H>0$、$\eta<0$ とする。$f=\log a$ の微分は $f'=-1/\eta,f''=1/\eta^2$ である。平坦計量の Schouten は零なので、変換則に代入すると $\widetilde P_{\eta\eta}=-1/(2\eta^2)$、空間対角は $1/(2\eta^2)$ となる。従って $\widetilde P=H^2g/2$ であり、四次元の縮約で $R_{ab}=3H^2g_{ab}$ となる。共形に平坦という条件は、等長に平坦という条件と異なる。
:::

### Weyl の発散と Cotton テンソル

$n\ge3$ とし、$P_{ab}=(R_{ab}-Rg_{ab}/[2(n-1)])/(n-2)$ とする。先の跡の計算から $\nabla^aP_{ab}=\partial_b\tr_g P$ である。縮約した微分 Bianchi 恒等式を添字順を保って書くと $\nabla^a R_{abcd}=\nabla_c R_{bd}-\nabla_dR_{bc}$ である。分解の $B(P)$ を直接微分すると
$$
\nabla^a B(P)_{abcd}=\nabla_c P_{bd}-\nabla_dP_{bc}
-g_{bc}\partial_d\tr_gP+g_{bd}\partial_c\tr_gP.
$$
一方 $R_{bd}=(n-2)P_{bd}+(\tr_gP)g_{bd}$ なので、差を取って
$$
\nabla^a C_{abcd}=(n-3)(\nabla_cP_{bd}-\nabla_dP_{bc})
$$
を得る。括弧の量は :term[Cotton テンソル]{en="Cotton tensor" yomi="こっとんてんそる"} の一つの添字規約である。四次元の Ricci が計量の定数倍であればこの微分が消える。一方、三次元では Weyl が恒等的に零でも、右辺の係数が零なので Cotton まで消えるとは限らない。

:::example[三次元で点ごとの曲率とその微分を比べる]
計量を直線と二次元曲面の直積 $\dd z^2+\dd x^2+F(x)^2\dd y^2$ とし、$F>0$ とする。二次元因子の曲率を $K=-F''/F$ とすると $R_{zz}=0$、曲面上の Ricci は $K$ 倍の計量、$R=2K$ である。三次元の Schouten は $P_{zz}=-K/2$、曲面成分は $Kg_{ij}/2$ となる。直積の接続には混合項がないので、$\nabla_xP_{zz}-\nabla_zP_{zx}=-K'/2$ である。例えば $F=e^{x^2}$ では $K=-2-4x^2$、この Cotton 成分は $4x$ となる。Weyl の成分数だけを調べてもこの変化は拾えない。
:::

### Weyl の電気型と磁気型

四次元の正規直交枠 $(u,e_i)$ を向き付きで選ぶ。$E_{ij}=C_{i0j0}$ とし、空間の完全反対称記号 $\varepsilon_{123}=1$ で $B_{ij}=\frac12\varepsilon_i{}^{kl}C_{klj0}$ とする。$E$ の対称性は対交換から、無跡性は Weyl の縮約零から従う。$B$ の跡 $\varepsilon^{ikl}C_{kli0}/2$ は第一 Bianchi の巡回和のため零になる。また反対称部を $\varepsilon^{ijm}$ で縮約すると
$$
\varepsilon^{ijm}B_{ij}
=\tfrac12(\delta^{jk}\delta^{ml}-\delta^{jl}\delta^{mk})C_{klj0}
=\sum_j C_{jmj0}=0 .
$$
従って $B$ も対称無跡である。各々五成分を持つ。混合ブロックは $C_{klj0}=\varepsilon_{kl}{}^iB_{ij}$ で復元できる。残る空間ブロックを $C_{ijkl}=\varepsilon_{ijm}\varepsilon_{kln}D_{mn}$ と書くと、対交換から $D$ は対称である。跡の条件 $\sum_k C_{kikj}=E_{ij}$ は $\delta_{ij}\tr D-D_{ij}=E_{ij}$ だから $\tr D=0$、$D=-E$ となる。従ってこの二つの行列で Weyl 全体が決まる。

:::example[真空潮汐を二つの行列から再構成する]
ある観測者に対して $E=\diag(-2b,b,b)$、$B=0$ なら、$C_{0101}=-2b$、$C_{0202}=C_{0303}=b$ である。空間成分は $C_{2323}=2b$、$C_{3131}=C_{1212}=-b$ となり、時間添字が一つだけの成分は零である。例えば $R_{11}$ の縮約は $-C_{0101}+C_{2121}+C_{3131}=2b-b-b=0$ になる。曲率二乗は $4$ 倍の六成分の二乗和で $48b^2$ である。真空でも局所の潮汐測定で曲率は検出できる。
:::

### 共形平坦性の十分条件を局所に示す

$n\ge3$ とする。共形に平坦なら Weyl は零である。逆に $C=0$ とし、さらに三次元の場合は Cotton が零と仮定する。四次元以上では Weyl の発散の式から自動的に Cotton が零になる。これらの条件の下で元の計量を $e^{2f}$ 倍し、曲率を零にする $f$ を局所で求める。$u_a=\partial_a f$ と書くと、$\widetilde P=0$ の条件は
$$
\nabla_a u_b=P_{ab}+u_a u_b-\tfrac12(u_cu^c)g_{ab}=:S_{ab}(x,u)
$$
である。右辺は対称だから、この一階系の解 $u$ があれば $\dd u=0$ となる。星形の座標領域で閉一形式を積分すれば $f$ を得る。

この一階系の可積分条件を調べる。$\nabla_a S_{bc}-\nabla_b S_{ac}$ を展開し、現れる $\nabla_a u_b$ を $S_{ab}$ で置き換える。$S_{ab}u_c$ の項は反対称化で消え、残りは
$$
\nabla_aP_{bc}-\nabla_bP_{ac}+u_b S_{ac}-u_a S_{bc}
-u^d S_{ad}g_{bc}+u^d S_{bd}g_{ac}
$$
である。$u$ の三次項は $S$ の二次項を入れると消える。第一・第二項は Cotton 零により消え、残りは $-u^d B(P)_{dcab}$ になる。$C=0$ ならこれは $-u_d R^d{}_{cab}$ であり、一形式への曲率作用とちょうど一致する。

この一致は、座標と $u$ の成分を独立変数とする空間で
$$
D_a=\partial_{x^a}+\big(\Gamma^c{}_{ab}u_c+S_{ab}\big)\partial_{u_b}
$$
と書いたとき $[D_a,D_b]=0$ となることを意味する。実際、括弧の座標成分は零で、縦の成分は二通りの二階偏微分の差であり、上の曲率作用を右辺の微分から引いたものである。従って一点の初期値から $D_a$ の局所流を順に合成できる。可換性により順序の違いで結果が変わらず、座標 $x^a$ への射影は各流でその座標だけを平行移動する。よって合成でできる積分面は $u=u(x)$ のグラフとなり、所要の一階系を解く。$\tilde g=e^{2f}g$ では $\widetilde P=0$ と $\widetilde C=0$ が両方成立するので全曲率が零となる。逆向きに三次元で共形に平坦なら、上の一階系を満たす $u$ が既に存在し、同じ可積分条件で $C=0$ と合わせて Cotton が零だと分かる。従って、四次元以上では Weyl の消滅、三次元では Cotton の消滅が局所共形平坦性を特徴付ける。

:::example[平坦計量に戻せる曲がった計量]
四次元の $g=e^{2H t}\eta$ では $f=-Ht$ を選べば $e^{2f}g=\eta$ になる。$g$ の Schouten は $P_{tt}=H^2/2$、$P_{ij}=H^2\delta_{ij}/2$ である。例えば $u_t=-H$、$u_i=0$ を一階系に入れると、時間成分は $S_{tt}=H^2$、空間成分は $S_{ij}=H^2\delta_{ij}$ になる。$\nabla_a u_b=-\Gamma^t{}_{ab}u_t$ と $\Gamma^t{}_{tt}=H$、$\Gamma^t{}_{ij}=H\delta_{ij}$ からも同じ値を得る。共形平坦性を示す一階系の解を直接書ける例である。
:::

## 測地線偏差と潮汐力

測地線族 $x(\lambda,s)$ の接線を $u=\partial_\lambda x$、分離を $\xi=\partial_s x$ とする。各測地線の $\lambda$ はアフィン径数とする。二つの変分方向はパラメータ面上で可換であり、ねじれ零から $D_\lambda\xi=D_su$ である。パラメータ面上の交換子 $D_\lambda D_s-D_s D_\lambda$ を $u$ に作用させ、$D_\lambda u=0$ を使えば、次の式になる。曲線が交差する点でも、微分はパラメータ面に引き戻した接続を使うため意味を持つ。

$$
\nabla_u^2\xi=R(u,\xi)u,\qquad
(\nabla_u^2\xi)^a=-R^a{}_{bcd}u^b\xi^c u^d .
\label{eq-gr-jacobi}
$$

この :term[測地線偏差]{en="geodesic deviation" yomi="そくちせんへんさ"} は有限の距離の力ではなく、隣接した軌道の一次差である。基準となる時間的測地線に沿い、正規直交枠を平行移動して時間軸を接線に取る。この自由落下枠では $\dd^2\xi^i/\dd\tau^2=-R^i{}_{0j0}\xi^j$ となる。例えば弱く静的な場で $g_{00}=-(1+2\Phi)$ とおくと、一次では $\Gamma^i{}_{00}=\partial_i\Phi$、$R^i{}_{0j0}=\partial_j\partial_i\Phi$ である。$\Phi=-GM/r$ なら半径方向に $2GM/r^3$、接線方向に $-GM/r^3$ の係数を得る。逆に定曲率 $K>0$ の時空では $-R^i{}_{0j0}=K\delta^i_j$ なので、一定方向の初期分離 $L$、初相対速度零に対して $\xi=L\cosh(\sqrt K\tau)$ となる。

### 正規座標の計量を二次まで求める

正規座標の原点から初速度 $v$ を少し変えると、終点の変化は Jacobi 場で得られる。そこで放射測地線 $\gamma(s)=\exp_p(sv)$ に沿って平行な枠を取り、$J(0)=0$, $D_sJ(0)=w$ の変分を考える。偏差式 :eqref[eq-gr-jacobi] の $D_s^2J=R(v,J)v$ から $D_s^2J(0)=0$、$D_s^3J(0)=R(v,w)v$ である。従って平行枠で $J(s)=s w+s^3 R(v,w)v/6+\order{s^4}$ になる。ここで正規座標の座標ベクトル $\partial_a$ は初速度を変える $J$ の $s$ 倍を除いたものである。二本の座標ベクトルの内積を取ると、曲率の対の交換と反対称性により

$$
g_{ab}(x)=\eta_{ab}-\tfrac13 R_{acbd}(p)x^cx^d+\order{\abs{x}^3}.
\label{eq-gr-normal-metric}
$$

従って計量の二階偏微分は任意ではなく、測地線を座標上の動径直線に固定すると曲率で決まる。

:::example[定曲率での二次項]
式の曲率に $R_{abcd}=K(\eta_{ac}\eta_{bd}-\eta_{ad}\eta_{bc})$ を入れると $g_{ab}=\eta_{ab}-K(\eta_{ab}\eta(x,x)-x_ax_b)/3+\order{\abs{x}^3}$ である。単位球の正定値計量に同じ導出を適用すると動径に垂直な座標長が $\sin r/r$ 倍になるので、その二乗は $1-r^2/3+\order{r^4}$ となる。動径方向には補正がなく、測地線のパラメータがそのまま半径の長さである。
:::

### 偏差方程式の初期値と有限基線

測地線の族を $x(\tau,s)$ とし、分離場 $J=\partial_s x$ を作る。隣の曲線との座標差は $sJ+\order{s^2}$ であり、偏差式はこの一次の分離に作用する。例えば平行枠で一定の電気型曲率 $R_{i0j0}=\diag(-2b,b,b)$ が与えられる模型では、$J^i(0)=J^i_0$、$\dot J^i(0)=0$ の解は
$$
J^1=J^1_0\cosh(\sqrt{2b}\,\tau),\qquad
J^{2,3}=J^{2,3}_0\cos(\sqrt b\,\tau)
$$
である。これは $b>0$ とした場合で、各固有方向に分けて二階微分すれば検算できる。実際に $b$ が変わるときは上の解を全区間には使えない。短い区間なら $J^i=J^i_0-\tau^2R^i{}_{0j0}(0)J_0^j/2+\order{\tau^3}$ となり、その点の潮汐を抽出できる。

:::example[曲率の跡と体積の初期変化]
三方向の初期分離を単位とする。上の解が囲む体積は $\cosh(\sqrt{2b}\,\tau)\cos^2(\sqrt b\,\tau)$ 倍である。二次では前の因子の $+b\tau^2$ と後ろの $-b\tau^2$ が消えるが、四次まで展開すると $1-b^2\tau^4/2+\order{\tau^6}$ となる。初期の跡が零でも個々の長さは変わるため、後の体積に変化が現れる。
:::

### 偏差方程式の解を測地線族から作る

アフィンな測地線 $\gamma$ に初期位置と初速度の小さな変化を与える。位置の変化を $\xi_0$、初速度の共変変化を $\zeta_0$ とする。初期位置を曲線 $p(s)$ に沿って動かし、元の速度と $s\zeta_0$ を加えたベクトルをその曲線に沿って平行移動して初速度に用いれば、二つの値を独立に指定できる。常微分方程式の初期値依存性により変分場は $\xi(0)=\xi_0$、$D_\tau\xi(0)=\zeta_0$ となる。偏差式は線形二階方程式で初期値が一意だから、任意のその解はこの方法で測地線族の一次変分として得られる。

:::example[平坦時空の焦点は初期速度だけでも生じる]
平坦な時空で直線の族 $\gamma_s(\lambda)=(\lambda,s(1-\lambda/T),0,0)$ を取る。$|s|<T$ なら時間的であり、分離場は $\xi^1=1-\lambda/T$ である。$R=0$ の偏差方程式を満たしながら、$\lambda=T$ ですべてが交わる。これは初期位置が異なる測地線を、内向きの速度で放したためである。これに対し同じ始点から異なる方向へ出る族の分離は $\xi=\lambda \zeta_0$ であり、非零 $\zeta_0$ では二度目の零点を持たない。交差の意味は初期条件と合わせて読む必要がある。
:::

### 正規座標で体積と発散を展開する

正規座標の二次の計量補正を $h_{ab}=-R_{acbd}x^cx^d/3$ とする。行列式を展開する際に二つの列から同時に $h$ を取ると四次になるため、二次精度では $\det(\eta+h)=\det\eta(1+\eta^{ab}h_{ab})$ である。平方根を取ると
$$
\sqrt{-\det g(x)}=1-\frac16 R_{cd}(p)x^cx^d+\order{\norm{x}^3}
$$
である。正規座標で $X=x^a\partial_a$ を動径ベクトル場とする。体積での発散は $\nabla_a X^a=4+X(\log\sqrt{-\det g})$ だから、$4-R_{cd}x^cx^d/3+\order{\norm{x}^3}$ となる。この式の補正では Riemann の全成分ではなく Ricci が現れる。

:::example[空間の小球の面積と体積を照合する]
半径 $L$ の正定値三次元球面で曲率が $1/L^2$ の場合、測地半径 $s$ の球面の面積は $4\pi L^2\sin^2(s/L)$ である。従って面積は $4\pi s^2(1-s^2/(3L^2)+\order{s^4/L^4})$、体積は $4\pi s^3(1-s^2/(5L^2)+\order{s^4/L^4})/3$ となる。正規座標では $R_{ij}=2\delta_{ij}/L^2$ だから密度は $1-\norm{x}^2/(3L^2)$ から始まり、Euclid 球でこの密度を積分しても同じ係数を得る。ここでは正の計量行列式を使う。
:::

### 世界線に沿う Fermi 正規座標の二次項

単位時間的測地線 $z(\tau)$ に平行な空間枠 $e_i$ を取り、$x=\exp_{z(\tau)}(X^i e_i)$ と置く。この局所座標で計量の二次項を得るには、空間放射測地線の初期点と初速度を別々に変分する。$V=X^i e_i$ とし、放射測地線を $s\in[0,1]$ で表す。始点を $\tau$ に沿って動かす変分場 $J_0$ の初期値は $(u,0)$ である。初期速度の変化が零になるのは、枠を平行移動しているからである。空間座標を変えるときの場 $J_i$ の初期値は $(0,e_i)$ である。双方に偏差式を適用し平行枠で積分すると、終点で
$$
J_0=u+\tfrac12 R(V,u)V+\order{\norm X^3},\qquad
J_i=e_i+\tfrac16 R(V,e_i)V+\order{\norm X^3}
$$
となる。実際、$J=J(0)+sJ'(0)+\int_0^s(s-s')R(V,J(s'))V\,\dd s'$ であり、曲率の項は $X$ の二次である。曲率を始点で評価し、$J(s')$ を先頭の初期値で置き換えると、$J_0$ では $\int_0^1(1-s')\dd s'=1/2$、$J_i$ では $\int_0^1(1-s')s'\dd s'=1/6$ が現れる。曲率の評価点を動かした差は $X$ の三次以上になる。

終点で $J_A$ は座標基底である。内積を取り、対の対称性を使うと
$$
\begin{aligned}
g_{00}&=-1-R_{0i0j}X^i X^j+\order{\norm X^3},\\
g_{0i}&=-\tfrac23 R_{0jik}X^j X^k+\order{\norm X^3},\\
g_{ij}&=\delta_{ij}-\tfrac13 R_{ikjl}X^k X^l+\order{\norm X^3}.
\end{aligned}
$$
曲率成分は $\tau$ における枠で評価する。例えば零交差成分の二次の二つの寄与は $1/2$ と $1/6$ が同じ $-R_{0jik}X^jX^k$ を掛けている。この座標の中心線では各時刻の計量が $\eta$ であり、中心から少し離れた時計では $g_{00}$ の二次項が時刻の率を変える。中心線に沿う曲率が変わる場合でも、横方向の展開は各時刻で適用できる。

:::example[一定曲率で閉じた形の時間成分を求める]
$K$ が一定の時空で、$X^i$ の大きさを $\ell$ とする。横に出る放射線の速度 $V$ は空間的で $g(V,V)=\ell^2$ であり、$R(V,J)V=-K\ell^2 J$ が $V$ に直交する場で成立する。$J_0$ はそのような場で、初期の導関数は零である。従って $K>0$ なら $J_0(1)=\cos(\sqrt K\ell)u$、$K<0$ なら $\cosh(\sqrt{-K}\ell)u$ となる。計量の時間成分はそれぞれ $-\cos^2(\sqrt K\ell)$、$-\cosh^2(\sqrt{-K}\ell)$ である。どちらも $-1+K\ell^2+\order{\ell^4}$ となり、$R_{0i0j}=-K\delta_{ij}$ を入れた式と一致する。$K>0$ で余弦が零となる距離ではこの座標の微分が退化し、局所展開をそのまま座標の大域的な有効範囲とはみなせない。
:::

### 固有時の二次変分に現れる潮汐項

時間的測地線 $z(\tau)$ の端点を固定し、変分の速度 $J$ が $u=\dot z$ と直交するとする。変分曲面の接線を $U=\partial_\tau z_s$ と書けば、固有時の被積分関数は $\sqrt{-g(U,U)}$ である。$s=0$ で二回微分すると
$$
-g(D_sU,D_sU)-g(u,D_s^2U)-g(u,D_sU)^2
$$
となる。ねじれ零より $D_sU=D_\tau J$ で、$g(u,D_\tau J)=0$ は直交条件と平行な $u$ から従う。また $D_s^2U=D_\tau D_sJ+R(J,u)J$ である。最初の項の $u$ との積は全微分で、固定端点では積分して零となる。従って平行枠での第二変分は
$$
\delta^2\tau=\int\big[-\delta_{ij}\dot J^i\dot J^j+R_{i0j0}J^iJ^j\big]\dd\tau .
$$
潮汐項は負にも正にもなり得る。例えば $K=-b<0$ が一定なら $R_{i0j0}=b\delta_{ij}$ であり、長さ $T$ の区間の $J=\sin(\pi\tau/T)e_1$ について $\delta^2\tau=T[b-(\pi/T)^2]/2$ である。$T>\pi/\sqrt b$ ではこの変分で固有時が増える。第一変分は零であっても長区間で最大性を失う例であり、$T=\pi/\sqrt b$ で零になる変分方向は偏差式を満たす。[共役点と焦点](/physics/general-relativity/causal-structure/) ではこの判定を一般の曲率へ拡張する。

### 初期の三本の分離から体積の変化を測る

平行な空間枠で、三本の分離ベクトルを列に持つ行列を $J$ とする。偏差式は $J''=-EJ$ であり、$E^i{}_j=R^i{}_{0j0}$ である。$\det J\ne0$ の区間では $B=J'J^{-1}$ を定義できる。逆行列の微分から $B'=-E-B^2$ を得る。また行列式の微分から $(\log|\det J|)'=\tr B$ なので
$$
(\log|\det J|)''=-R_{00}-\tr(B^2)
$$
である。ここで $R_{00}=\tr E$ を使った。$B^2$ の跡は $B$ が対称であれば非負だが、一般の実行列では非負にならない。対称性を使える一例は初期に $J=I$、$J'=J'^{\mathsf T}$ の場合である。実際、$B^{\mathsf T}$ は対称な $E$ の下で同じ行列方程式を満たし、同じ初期値からの解の一意性で $B^{\mathsf T}=B$ となる。分離の体積を扱う際に必要な条件は [測地線の束](/physics/general-relativity/causal-structure/) でさらに調べる。

:::example[初期の膨張がある平坦な束]
$E=0$ で $J(0)=I,J'(0)=bI$ とすると $J=(1+b\tau)I$ である。従って $\log|\det J|=3\log|1+b\tau|$、その二階微分は $-3b^2/(1+b\tau)^2$ である。これは $B=bI/(1+b\tau)$ を代入した式と一致する。$b<0$ なら $\tau=-1/b$ で体積が零になり、それ以後に同じ式で $B$ を滑らかに保つことはできない。平坦な場合でも初期の相対速度によって、このような焦点が生じ得る。
:::

## 正規直交枠と Cartan の構造方程式

:term[正規直交枠]{en="orthonormal frame" yomi="せいきちょっこうわく"} $e_A$ の双対 $\theta^A$ に対し $\nabla_X e_B=\omega^A{}_B(X)e_A$ と置く。計量整合性は $\omega_{AB}=-\omega_{BA}$ を与える。$\dd\theta^A(X,Y)=X\theta^A(Y)-Y\theta^A(X)-\theta^A([X,Y])$ を使うと $\theta^A(T(X,Y))=(\dd\theta^A+\omega^A{}_B\wedge\theta^B)(X,Y)$ となる。同様に $R(X,Y)e_B$ に積の微分則を使うと

$$
\dd\theta^A+\omega^A{}_B\wedge\theta^B=0,\qquad
\Omega^A{}_B=\dd\omega^A{}_B+\omega^A{}_C\wedge\omega^C{}_B
=\tfrac12 R^A{}_{BCD}\theta^C\wedge\theta^D .
\label{eq-gr-cartan}
$$

これが :term[Cartan の構造方程式]{en="Cartan structure equations" yomi="かるたんのこうぞうほうていしき"} である。球面で $\theta^1=\dd\theta,\theta^2=\sin\theta\dd\phi$ と取れば $\omega^1{}_2=-\cos\theta\dd\phi$、$\Omega^1{}_2=\sin\theta\dd\theta\wedge\dd\phi$ となる。平面で $\theta^2=r\dd\phi$ と取ると $\omega^1{}_2=-\dd\phi$ であり、曲率は零である。接続が零かどうかと曲率が零かどうかを区別できる。

### 基底が交換しないときの接続形式

正規直交共枠を $\theta^a$、双対の枠を $e_a$ とする。添字付きの一次形式 $\omega^a{}_b$ を $\nabla_X e_b=\omega^a{}_b(X)e_a$ で定める。ねじれが零なら
$$
\dd\theta^a(e_b,e_c)=-\theta^a([e_b,e_c])
=-\omega^a{}_c(e_b)+\omega^a{}_b(e_c).
$$
従って $\dd\theta^a+\omega^a{}_b\wedge\theta^b=0$ である。計量整合性から $\omega_{ab}+\omega_{ba}=0$ である。回転する枠で $\dd\theta^a$ が非零でも、ねじれはこの接続項と合計して判定する。曲率を $R(X,Y)e_b$ の係数として展開すると、$X[\omega^a{}_b(Y)]-Y[\omega^a{}_b(X)]-\omega^a{}_b([X,Y])$ と接続の積差になる。前者は外微分、後者はくさび積なので第二構造方程式が得られる。

:::example[双曲平面の共枠]
上半平面 $y>0$ の計量 $(\dd x^2+\dd y^2)/y^2$ で $\theta^1=\dd x/y,\theta^2=\dd y/y$ とすると $\dd\theta^1=\theta^1\wedge\theta^2,\dd\theta^2=0$ である。第一式から $\omega^1{}_2=-\theta^1$ を得て、反対称性から $\omega^2{}_1=\theta^1$ となる。二次元では積項がなく、$\Omega^1{}_2=-\theta^1\wedge\theta^2$ なので断面曲率は $-1$ である。単位球面の $+1$ と符号が変わるのは、測地線の分離の仕方も反対になることと対応する。
:::

### 曲率形式の変換と Bianchi 恒等式

枠を $e'_b=e_c S^c{}_b$ と変え、接続形式を行列 $\omega$ と見なす。ベクトルを微分する積の法則から $\omega'=S^{-1}\omega S+S^{-1}\dd S$ になる。点ごとに違う $S$ を使うと第二項があるため、接続形式は同じ型のベクトル値一次形式の単純な基底変換ではない。一方、$\Omega=\dd\omega+\omega\wedge\omega$ に $\omega'$ を代入する。$\dd S^{-1}=-S^{-1}(\dd S)S^{-1}$ と $\dd^2S=0$ を使うと、$\dd S$ を含む項は $\dd\omega'$ と積の間で消え、$\Omega'=S^{-1}\Omega S$ である。これは曲率がテンソルであることの枠での表示である。

さらに $\dd(\omega\wedge\omega)=\dd\omega\wedge\omega-\omega\wedge\dd\omega$ なので
$$
\dd\Omega+\omega\wedge\Omega-\Omega\wedge\omega=0
$$
となる。$\Omega$ を展開すれば、$\dd\omega$ を含む二組が逆符号、三重の積も逆符号で消える。この式を三本の接ベクトルで評価したものが接続形式に関する微分 Bianchi 恒等式である。

:::example[回転する直交枠の接続は非零でも曲率は零]
平坦な二次元で直線座標の枠を $S=\begin{pmatrix}\cos u&-\sin u\\\sin u&\cos u\end{pmatrix}$ で回す。元の接続は零なので $\omega'=S^{-1}\dd S=\begin{pmatrix}0&-1\\1&0\end{pmatrix}\dd u$ である。$\dd\omega'=0$、積も $\dd u\wedge\dd u=0$ であり、曲率は零である。一方、$u$ に沿って変わる枠成分の平行輸送には $\omega'$ の項が必要である。接続の非零さは幾何の曲率だけを表す量ではない。
:::

### 構造係数から接続形式を一意に決める

枠の括弧を $[e_B,e_C]=c^A{}_{BC}e_A$ と書き、$c_{ABC}=\eta_{AD}c^D{}_{BC}$ とする。$\omega_{AB}=\omega_{AB C}\theta^C$ の未知係数は第一・第二添字について反対称である。第一構造方程式を $(e_B,e_C)$ で評価すると $\omega_{AC B}-\omega_{AB C}=c_{ABC}$ になる。三つの添字を巡回させて足し引きすると
$$
2\omega_{AB C}=c_{ACB}-c_{CBA}+c_{BAC}
$$
である。これは正規直交枠に Koszul の式を適用しても得られる。例えば $\theta^0=N(x)\dd t$、$\theta^1=\dd x$ なら $[e_1,e_0]=-(N'/N)e_0$ であり、$\omega^0{}_1=(N'/N)\theta^0$、$\omega^1{}_0=\omega^0{}_1$ を得る。両者の符号が同じなのは時間添字を下げる際に負号が出るからである。

:::example[静的ラプスと回転対称な空間の組合せ]
$\theta^0=N(x)\dd t,\theta^1=\dd x,\theta^2=F(x)\dd y,\theta^3=F(x)\dd z$ とし、$N,F>0$ とする。$a=N'/N,b=F'/F$ と書く。外微分は $\dd\theta^0=a\theta^1\wedge\theta^0$、$\dd\theta^I=b\theta^1\wedge\theta^I$（$I=2,3$）である。第一構造方程式から
$$
\omega^0{}_1=a\theta^0,\quad \omega^1{}_0=a\theta^0,\quad
\omega^I{}_1=b\theta^I,\quad \omega^1{}_I=-b\theta^I
$$
を得る。第二構造方程式で横方向との混合項を追うと、
$$
\begin{aligned}
\Omega^0{}_1&=-(a'+a^2)\theta^0\wedge\theta^1,\\
\Omega^0{}_I&=-ab\,\theta^0\wedge\theta^I,\\
\Omega^1{}_I&=-(b'+b^2)\theta^1\wedge\theta^I,\\
\Omega^2{}_3&=-b^2\theta^2\wedge\theta^3 .
\end{aligned}
$$
例えば $\Omega^0{}_I$ には $\dd\omega^0{}_I$ がなく、$\omega^0{}_1\wedge\omega^1{}_I$ だけが寄与する。$N=F=e^{kx}$ なら $a=b=k$ は一定で、六平面の断面曲率はいずれも $-k^2$ である。時間添字のある $\Omega$ から下添字曲率を読むときは第一添字を下げる負号を含める。その結果、この例では $R_{0101}=k^2$ であっても断面曲率は $-k^2$ になる。
:::

### 回転しながら運んだ枠の曲率を再確認する

平坦領域で選んだ正規直交枠を点ごとに変えると、接続の値を零とする枠と非零とする枠がある。例えば Minkowski 平面で Cartesian 枠に boost $S=e^{\chi J}$ を作用させる。$J^0{}_1=J^1{}_0=1$、その他の成分を零とする。変換則から $\omega'=J\dd\chi$ である。外微分は $\dd^2\chi=0$、積項も同じ一形式の外積なので零となる。従って boost の速さが場所によって変わっても曲率は零である。直交枠の移動に伴う回転や boost と、計量そのものが持つ曲率を区別できる。

:::example[正規直交枠で一次と二次の構造方程式を連続して使う]
曲がった時空の二次元部分を $-\dd t^2+a(t)^2\dd x^2$、残り二方向を平坦な直積とする。$\theta^0=\dd t,\theta^1=a\dd x$ を使うと $\dd\theta^1=(\dot a/a)\theta^0\wedge\theta^1$ である。第一式と $\omega^0{}_1=\omega^1{}_0$ から $\omega^1{}_0=(\dot a/a)\theta^1$ になる。従って
$$
\Omega^0{}_1=\left[\dv{}{t}(\dot a/a)+(\dot a/a)^2\right]\theta^0\wedge\theta^1
=\frac{\ddot a}{a}\theta^0\wedge\theta^1.
$$
時間添字を下げると $R_{0101}=-\ddot a/a$ なので、この二平面の断面曲率は $+\ddot a/a$ である。例えば $a=t$、$t>0$ なら接続形式は $\dd x$ と非零だが外微分は零になる。一方、$a=e^{Ht}$ では $\omega^1{}_0=H\theta^1$ であり、係数 $H$ が一定でも $\theta^1$ 自体の外微分により $\Omega^1{}_0\ne0$ となる。枠の係数だけを微分するとこの項を見落とす。
:::

## 演習問題

::::exercise{level=basic}
半径 $b$ の球面のスカラー曲率を求めよ。

:::solution
枠を $\theta^1=b\dd\theta,\theta^2=b\sin\theta\dd\phi$ とすると接続は単位球と同じで、曲率形式は $b^{-2}\theta^1\wedge\theta^2$ である。二次元では $R=2R_{1212}=2/b^2$ となる。
:::
::::

::::exercise{level=basic}
四次元で $K=-2$ の定曲率計量の $R$ と $G_{ab}$ を求めよ。

:::solution
縮約は $R=12K=-24$、$G_{ab}=-3Kg_{ab}=6g_{ab}$ である。
:::
::::

::::exercise{level=basic}
$R_{0123}=2,R_{0231}=5$ のとき $R_{0312}$ を求めよ。

:::solution
第一 Bianchi 恒等式から和が零なので $R_{0312}=-7$ となる。
:::
::::

::::exercise{level=basic}
真空弱重力で $\Phi=3(x^2-y^2)$ の測地線偏差行列を求めよ。

:::solution
Hessian は $\diag(6,-6,0)$ である。偏差行列はその負の $\diag(-6,6,0)$ である。
:::
::::

::::exercise{level=advanced}
二次元の静的計量 $-\!A(z)^2\dd t^2+\dd z^2$ の $R$ を求めよ。

:::solution
$\Gamma^z{}_{tt}=AA'$、$\Gamma^t{}_{zt}=A'/A$ より $R^z{}_{tzt}=(AA')'-AA'(A'/A)=AA''$ である。署名を含めて $R=2g^{zz}g^{tt}R_{ztzt}=-2A''/A$ となる。
:::
::::

::::exercise{level=advanced}
定曲率の時間的測地線束の偏差で $K=-4$、$L\ne0$、$\xi(0)=L$、$\xi'(0)=0$ の解と最初の零点を求めよ。

:::solution
$\xi''=-4\xi$ の一般解は $a\cos(2\tau)+b\sin(2\tau)$ である。初期値から $\xi=L\cos(2\tau)$、$L\ne0$ のとき最初の正の零点は $\pi/4$ となる。
:::
::::

::::exercise{level=advanced}
四次元の定曲率計量で Weyl 曲率が零であることを示せ。

:::solution
$P=(3K-2K)g/2=Kg/2$ であり、$B(P)=K(g_{ac}g_{bd}-g_{ad}g_{bc})$ は元の曲率に一致する。このため差が零である。
:::
::::

::::exercise{level=advanced}
$\omega^1{}_2=-f'(r)\dd\phi$ で表される $\dd r^2+f(r)^2\dd\phi^2$ の断面曲率を求めよ。

:::solution
$f>0$ の領域で曲率形式は $-f''\dd r\wedge\dd\phi$ である。枠に直して $\Omega^1{}_2=-(f''/f)\theta^1\wedge\theta^2$ となるので $K=-f''/f$ である。
:::
::::

::::exercise{level=advanced}
$g=-N(z)^2\dd t^2+\dd z^2+\dd x^2+\dd y^2$ とする。正規直交枠での Ricci 成分と Einstein 成分を求め、真空になる条件を示せ。$N$ は正の滑らかな関数とする。

:::solution
$e_0=N^{-1}\partial_t,e_1=\partial_z$ とする。接続係数は $\Gamma^t{}_{tz}=N'/N$、$\Gamma^z{}_{tt}=NN'$ であり、$R^t{}_{ztz}=-(N'/N)'-(N'/N)^2=-N''/N$ である。第一添字を下げ枠へ直すと $R_{0101}=N''/N=B$ となる。残りの独立な曲率成分は零である。従って $R_{00}=B$、$R_{11}=-B$、$R_{22}=R_{33}=0$、$R=-2B$ である。ここで枠の $g_{00}=-1$ を使うため、$G_{00}=R_{00}+R/2=0$ となる。同様に $G_{11}=0$、$G_{22}=G_{33}=B$ である。真空は $N''=0$ の場合に限り、このとき全曲率も零である。例えば $N=1+az$ は加速座標で、$N=e^{az}$ では $B=a^2$ の曲率が残る。
:::
::::

## 確認テスト

::::quiz
曲率について任意の座標で成立する等式を選べ。

- [x] $R_{abcd}=R_{cdab}$
- [ ] $R_{abcd}=R_{bacd}$
- [ ] $R_{abcd}=0$ が接続係数の非零性から従う。

:::explanation
正規座標の計量二階微分の表示で対を交換しても値は同じである。各対の内部は反対称である。接続は曲率が零でも極座標のように非零になり得る。
:::
::::

::::quiz{answer="20" tol="1e-9"}
四次元の Riemann 曲率の独立成分数を求めよ。

:::explanation
六つの添字対での対称行列に $21$ 成分があり、全反対称四形式に対応する第一 Bianchi 条件が一つある。したがって $20$ 成分である。
:::
::::

::::quiz
真空で $R_{ab}=0$ とする。曲率のうち必ず零になるものをすべて選べ。

- [x] スカラー曲率
- [x] Einstein テンソル
- [ ] Weyl テンソル

:::explanation
$R$ と $G_{ab}$ は Ricci の縮約と線形結合なので零である。このとき $P=0$ であり $C=Riem$ なので、独立な跡なし曲率は残り得る。
:::
::::

::::quiz{answer="2" tol="1e-9"}
半径 $1$ の球面のスカラー曲率を求めよ。

:::explanation
正規直交枠で $R_{1212}=1$ であり、二次元のスカラー縮約は同じ値を二方向から数えるので $R=2$ である。
:::
::::

::::quiz{answer="50" tol="1e-7"}
五次元の Riemann テンソルの独立成分数を求めよ。

:::explanation
添字対は $\binom52=10$ 通りなので対称行列の成分は $55$ 個である。完全反対称化の像の次元 $\binom54=5$ を引いて $50$ 個になる。
:::
::::

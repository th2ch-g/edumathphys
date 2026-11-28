---
title: "漸近的平坦性と質量"
description: "無限遠の幾何を定め、静的質量、空間無限遠のエネルギー、放射によって減る Bondi 質量を計算する。"
notation:
  - symbol: "\\Omega"
    meaning: "共形因子"
  - symbol: "\\mathscr I^\\pm"
    meaning: "未来・過去の零無限遠"
  - symbol: "E_{\\rm ADM}"
    meaning: "空間無限遠のエネルギー"
  - symbol: "m_B,N_{AB}"
    meaning: "Bondi の質量面密度と news"
---

## この章の目標

- 零無限遠と空間無限遠を共形図で区別する。
- Komar 積分と ADM エネルギーを Schwarzschild 解で計算する。
- Bondi 質量が減る理由を角方向の積分から説明する。
- 正エネルギー定理と BMS 対称性の仮定を把握する。

## 共形コンパクト化と無限遠

無限遠で光がどこへ向かうかを調べるには、円錐を保ちながら距離の目盛りを縮める。:term[共形コンパクト化]{en="conformal compactification" yomi="きょうけいこんぱくとか"} では物理計量 $g$ に対し $\tilde g=\Omega^2g$ を用いる。$\Omega>0$ の範囲では零測地線の像が同じであり、$\Omega=0$ で $\tilde g$ が滑らかに延長する境界を調べる。

任意の基準長で単位を取った Minkowski 時空で $u=t-r$、$v=t+r$、$U=\arctan u$、$V=\arctan v$ とする。$\Omega=2\cos U\cos V$ とすると

$$
\tilde g=-4\dd U\dd V+\sin^2(V-U)\dd\Omega_2^2 .
$$

動径部は $\dd u=\sec^2U\dd U$、$\dd v=\sec^2V\dd V$ から得られ、角部は $2r=v-u$ を使う。$V=\pi/2$、$-\pi/2<U<\pi/2$ が未来零無限遠 $\mathscr I^+$、$U=-\pi/2$ の対応する面が過去零無限遠 $\mathscr I^-$ である。$i^+$ は $U=V=\pi/2$、$i^-$ は両方 $-\pi/2$、$i^0$ は $(U,V)=(-\pi/2,\pi/2)$ である。例えば $t$ 一定で $r\to\infty$ なら $i^0$ へ向かい、$u$ 一定で $r\to\infty$ なら $\mathscr I^+$ へ向かう。$\mathscr I^+$ は零面であり、そこで法線 $\dd V$ の逆計量ノルムが零である。

### アフィン径数は共形変換で変わる

$\tilde g=\Omega^2 g$ とし、$k=\dd x/\dd\lambda$ を物理計量のアフィン光的接線とする。二つの接続の差に $k$ を二度入れると $\tilde\nabla_k k=2k(\log\Omega)k$ である。$\tilde k=\Omega^{-2} k$ とすれば、係数を微分する項とこの項が相殺して $\tilde\nabla_{\tilde k}\tilde k=0$ となる。従って $\dd\tilde\lambda/\dd\lambda=\Omega^2$ である。例えば Minkowski 時空の外向き光で $\lambda=r$、$\Omega=1/r$ とすると、$\tilde\lambda$ の残りは $\int_r^\infty r'^{-2}\dd r'=1/r$ で有限である。光的無限遠を有限の座標で表せることと、物理的アフィン径数が無限であることが両立する。

### 光的無限遠の面が光的になること

:::example[零方向の共形計量を直接拡張する]
Minkowski 時空の遅れ時刻 $u=t-r$ で $\Omega=1/r$ とすると $\dd r=-\Omega^{-2}\dd\Omega$ である。従って $\Omega^2g=-\Omega^2\dd u^2+2\dd u\,\dd\Omega+\dd\Omega_2^2$ となる。ここで $\dd\Omega_2^2$ は単位球面の計量である。$u,\Omega$ ブロックの行列式は $-1$ で $\Omega=0$ でも非退化、逆計量の $\Omega\Omega$ 成分は $\Omega^2$ である。従って $\Omega=0$ の法線 $\dd\Omega$ は光的で、接方向にもなっている。この座標ではそのベクトルは $\partial_u$ である。
:::

### 共形ゲージと端点

Minkowski の変換で $T=U+V$、$\chi=V-U$ とすると、計量は $-\dd T^2+\dd\chi^2+\sin^2\chi\,\dd\Omega_2^2$ となる。条件 $U,V\in(-\pi/2,\pi/2)$ と $r\ge0$ は $0\le\chi<\pi$、$|T|+\chi<\pi$ に等しい。$\chi=0$ では球面が正則に縮み、物理時空の中心になる。一方、$\chi=\pi,T=0$ では球面が無限遠の点 $i^0$ に縮む。$\Omega=\cos T+\cos\chi$ なので $i^0$ と $i^\pm$ では $\dd\Omega=0$ である。従って $\dd\Omega\ne0$ を課す滑らかな光的境界の議論はこれらの端点を含まない。

同じ境界を表す $\Omega'=\omega\Omega$ を選ぶ。ただし $\omega$ は滑らかで正である。境界上では $\dd\Omega'=\omega\dd\Omega$、$\tilde g'{}^{-1}=\omega^{-2}\tilde g^{-1}$ なので、法線ベクトルは $n'=\omega^{-1}n$ である。境界への引き戻し計量は $q'=\omega^2q$ になる。従って法線の零性とその張る方向は共形因子の取り方によらない。例えば出射座標で $\Omega=1/r$ を $L/r$ に替えると、球の誘導計量は $L^2q$、生成子は元の $1/L$ 倍である。

### Einstein 方程式と境界の署名

四次元で $g=\Omega^{-2}\tilde g$ とし、$R[g]$ が $\Omega\to0$ で $4\Lambda$ に近づくと仮定する。共形変換のスカラー曲率式に $f=-\log\Omega$ を入れる。$\tilde\Box f=-\Omega^{-1}\tilde\Box\Omega+\Omega^{-2}|\tilde\nabla\Omega|^2$ だから
$$
R[g]=\Omega^2\tilde R+6\Omega\tilde\Box\Omega-12\tilde g^{-1}(\dd\Omega,\dd\Omega).
$$
$\tilde g,\Omega$ が滑らかなら先頭の二項は境界で零である。従って $\tilde g^{-1}(\dd\Omega,\dd\Omega)|_{\Omega=0}=-\Lambda/3$ となる。$\dd\Omega\ne0$ の部分で $\Lambda=0$ なら光的、$\Lambda>0$ なら空間的、$\Lambda<0$ なら時間的な境界を得る。真空なら $R[g]=4\Lambda$ が成り立つ。物質があれば、その応力の跡がどう減衰するかも調べる。

:::example[反 de Sitter の境界]
$f=1+r^2/L^2$ の静的計量で $x=1/r$ とする。共形計量の $tt$ 成分は $-(x^2+L^{-2})$、$xx$ 成分は $(x^2+L^{-2})^{-1}$、角成分は単位球面計量である。従って $x=0$ でも非退化で $\tilde g^{-1}(\dd x,\dd x)=L^{-2}>0$ である。これは $\Lambda=-3/L^2$ とした一般式に一致する。これに対し Schwarzschild で静的時刻を固定して同じ変換をすると、動径係数は $x^{-2}(1-2mx)^{-1}$ であり発散する。光的境界で非退化な計量を得るためには、時刻も遅延時刻へ変える必要がある。
:::

## 漸近的平坦な時空

:term[漸近的平坦性]{en="asymptotic flatness" yomi="ぜんきんてきへいたんせい"} はどの無限遠を見るかで条件が異なる。空間無限遠では、ある端の漸近 Cartesian 座標で $h_{ij}=\delta_{ij}+\order{r^{-1}}$、$\partial h=\order{r^{-2}}$、$\partial^2h=\order{r^{-3}}$、$K_{ij}=\order{r^{-2}}$ などを課す。質量積分の極限には曲率の可積分性、座標不変性には許容する座標の変化と落ち方も指定する。零無限遠では、滑らかな $\tilde g$、$\Omega=0$、$\dd\Omega\ne0$ の零境界を求め、放射を含む各 $u$ 一定の断面で量を定める。条件を満たす無限遠があっても内部が平坦になるとは限らない。

例えば Schwarzschild の外向き座標で $g=-f\dd u^2-2\dd u\dd r+r^2\dd\Omega_2^2$ とする。$x=1/r$ として $x^2g$ を計算すると $-x^2(1-2mx)\dd u^2+2\dd u\dd x+\dd\Omega_2^2$ である。$x=0$ で非退化、法線のノルムは $x^2(1-2mx)=0$ である。この展開で $m$ の寄与は $\tilde g_{uu}$ の $x^3$ 項に現れる。

### 空間方向と光方向を取り違えない

外向きの座標 $u=t-r$ を取る。固定した $t$ で $r\to\infty$ に進むと $u\to-\infty$ であり、一方、固定した $u$ で進むと $t=r+u\to\infty$ になる。例えば外向きの波形を $h(t,r)=F(t-r)/r$、$F$ は有界な区間に台を持つとする。固定時刻では十分遠方で $h$ が零になるが、$F(u)\ne0$ の固定した $u$ では $rh\to F(u)$ となる。従って同じ $1/r$ の展開でも、固定している変数を指定しなければ境界で測るデータが変わる。漸近無限遠で質量や放射を定義するときも、空間切片上の極限と零方向の極限を別々に取る。

### 波に伴う質量の境界データ

:::example[放射の前の切片と後の切片を比較する]
球対称の外向き零殻を理想化して、殻の前方を $m(u)=m_i$（$u<u_0$）、後方を $m(u)=m_f<m_i$（$u>u_0$）とする。空間端では、前方の Schwarzschild 領域の静的時刻 $t_i$ を固定して $r\to\infty$ を取る。この端では $u=t_i-r_{*,i}$ が $-\infty$ へ向かうため、漸近計量の質量は $m_i/G$ である。一方、未来光的無限遠の断面を殻の通過後の $u>u_0$ に取ると $m_f/G$ になる。二つの値の差は殻が運んだエネルギーであり、以下の滑らかな Vaidya 模型で流束積分からも得られる。
:::

### 漸近端と光線の終点

初期切片のコンパクト部分を除いて残る連結成分が $\R^3$ の球の外部と微分同相であるとする。その座標を $x^i$、$r=(\delta_{ij}x^ix^j)^{1/2}$ とし、無限大へ向かう部分を一つの :term[漸近端]{en="asymptotic end" yomi="ぜんきんたん"} と呼ぶ。拘束を満たすデータで $\rho$ の絶対積分が有限なら、$R^{(3)}=16\pi G\rho+K_{ij}K^{ij}-K^2$ も絶対可積分である。後二項は $r^{-4}$ 程度なので $r^2\dd r$ を掛けても積分できる。同様に $j_i$ の可積分性が運動量積分の極限を制御する。

滑らかな $\mathscr I^+$ が存在することを、全ての光線がそこへ逃げる条件と区別する。:term[漸近単純性]{en="asymptotic simplicity" yomi="ぜんきんたんじゅんせい"} は全ての延長不能な光的測地線が共形境界に未来端と過去端を持つことを含む。ブラックホールへ落ちる光線を含む時空全体は、滑らかな外部の光的無限遠を持ってもこの条件を満たさない。共形計量について境界を横断するアフィン光線では、横断性により $\Omega=a(\tilde\lambda_0-\tilde\lambda)+\order{(\tilde\lambda_0-\tilde\lambda)^2}$、$a>0$ である。従って物理径数 $\int\Omega^{-2}\dd\tilde\lambda$ は発散する。共形計量で有限の径数に端点を持ちながら、物理計量では無限の径数に及ぶ。

:::example[遅延時刻に含まれる対数]
Schwarzschild で $r_*=r+2m\log(r/(2m)-1)$、$u=t-r_*$ を固定する。このとき平坦時空の式から得た時刻 $t-r$ は $u+2m\log(r/(2m)-1)$ となり、$r\to\infty$ で一定値に向かわない。同じ計量が $g_{tt}=-1+\order{r^{-1}}$ を満たしても、光線に沿って時刻を作る際に $1/r$ を積分すると対数が生じる。空間側で使う近似座標から、光的境界で同じ断面を読むことはできない。
:::

## Komar 積分

時間並進の Killing 場 $\xi$ を無限遠で単位時間的に規格化する。:term[Komar 積分]{en="Komar integral" yomi="こまーせきぶん"} による質量を $M_K=-(8\pi G)^{-1}\int_S *\dd\xi^\flat$ とする。向きは Minkowski の時間と外向き球に合わせ、静的計量 $-N^2\dd t^2+h_{ij}\dd x^i\dd x^j$ なら $M_K=(4\pi G)^{-1}\int_S D_iN\,s^i\dd A$ となる。実際 $\xi^\flat=-N^2\dd t$ なので $\dd\xi^\flat=-2N\dd N\wedge\dd t$ であり、単位時間の一形式 $N\dd t$ を用いて Hodge 双対を取るとこの流束を得る。

Killing 場では $\nabla_b\nabla^b\xi^a=-R^a{}_b\xi^b$ が成り立つ。Killing 方程式で第一微分を移し、交換子の縮約と $\nabla_b\xi^b=0$ を使うからである。従って $\nabla_b(\dd\xi^\flat)^{ba}=-2R^a{}_b\xi^b$ で、真空では $*\dd\xi^\flat$ が閉形式となる。Stokes により真空領域で相同な面の積分は等しい。物質がある場合は曲率が流束の源になる。Schwarzschild では $N=\sqrt f$、$s^r=\sqrt f$ なので $s^iD_iN=f'/2=m/r^2$ となり、任意の外部球で $M_K=m/G=M$ である。Killing 場を定数倍すると積分も定数倍になるため、規格化が必要である。

### 静的球で Komar 積分の係数を検算する

静的計量 $-N^2\dd t^2+h_{ij}\dd x^i\dd x^j$ で、$\xi=\partial_t$ を採用する。法線 $n=N^{-1}\xi$ と切片内の外向き法線 $s$ に対し、$\nabla_a\xi_b=(\dd\xi^\flat)_{ab}/2$ から Komar 質量は $(4\pi G)^{-1}\int s^iD_iN\dd A$ と書ける。Schwarzschild 外部で $N=\sqrt f,s^r=\sqrt f$ を掛けると $s(N)=f'/2=m/r^2$ なので積分は $m/G$ で、どの球でも同じである。球を別の半径へ動かしても同じなのは $R_{ab}=0$ による。例えば荷電外部で $f=1-2m/r+q^2/r^2$ にすると同じ計算は $(m-q^2/r)/G$ を与える。二球で値が異なる分は間の電磁場の Ricci に由来する。

### 電磁場のエネルギーと Komar 質量の変化

:::example[電場の密度から二つの球の差を積分する]
静的な荷電外部で $8\pi G\rho=q^2/r^4$、電磁場の跡は零である。Komar の体積式には $2N\rho$ が入り、固有体積は $4\pi r^2 \dd r/\sqrt f$、$N=\sqrt f$ である。半径 $R_1,R_2$ の差は $\int_{R_1}^{R_2}8\pi r^2\rho\,\dd r=(q^2/G)(1/R_1-1/R_2)$ となる。これは面積積分で求めた $(m-q^2/R_2)/G-(m-q^2/R_1)/G$ と一致する。ここに係数 $2$ があるのは Komar の体積式が $T_{ab}$ だけでなく跡反転を通じて圧力にも依存するためである。
:::

### 物質がある場合の体積式

$\Lambda=0$ とし、空間的切片の未来単位法線を $n$ とする。Killing 恒等式と Stokes を用いると、外側と内側の球の差は
$$
M_K(S_2)-M_K(S_1)
=\frac{1}{4\pi G}\int_{\Sigma_{12}}R_{ab}n^a\xi^b\,\dd V_h
=2\int_{\Sigma_{12}}(T_{ab}-\tfrac12 T g_{ab})n^a\xi^b\,\dd V_h .
$$
係数は $\nabla_b(\dd\xi^\flat)^{ba}=-2R^a{}_b\xi^b$ の $2$ と未来時間的法線の符号から得る。静的な場合は $\xi=Nn$、$S=h^{ij}T_{ij}$、$T=S-\rho$ だから被積分は $N(\rho+S)$ となる。静止完全流体では $S=3p$ であり、圧力と赤方偏移の効果が入る。

同じ式を静的計量から確かめる。$\Gamma^t{}_{ti}=D_i\log N$、$\Gamma^i{}_{tt}=ND^iN$ である。$R_{tt}$ の空間微分と空間接続は $D_i(ND^iN)$ となり、残る積は $-|DN|^2$ なので $R_{tt}=N\Delta_hN$ である。従って $\Delta_hN=4\pi G N(\rho+S)$ となり、固有体積で積分すれば同じ面積分を得る。例えば弱い静的場で $N=1+\Phi+\order{\Phi^2}$、$p\ll\rho$ なら $\Delta\Phi=4\pi G\rho$ に近づく。

### 軸対称 Killing 場の角運動量

軸対称場 $\varphi$ の軌道を $2\pi$ 周期に規格化し、同じ向きの球で $J_K=(16\pi G)^{-1}\int *\dd\varphi^\flat$ と定義する。軸対称な切片に $\varphi$ が接していれば、$\nabla_s\varphi\cdot n=-K_{ij}s^i\varphi^j$ である。法線二平面で Hodge 双対を評価すると $J_K=-(8\pi G)^{-1}\int K_{ij}s^i\varphi^j\dd A$ を得る。二球の差を取る体積式では跡項が $n\cdot\varphi=0$ で消え、$J_K(S_2)-J_K(S_1)=-\int T_{ab}n^a\varphi^b\dd V_h=\int j_i\varphi^i\dd V_h$ となる。従って $K=0$ の静的切片では、Komar 質量が非零でも角運動量は零である。

時間並進と軸回転の両方があれば $\xi+\omega\varphi$ に対する質量の規格化の積分は $M_K-2\omega J_K$ である。$\omega$ は定数でなければならない。$\omega$ の勾配があると Killing 方程式の対称化に $\varphi_{(a}\partial_{b)}\omega$ が残るからである。特にこの組合せが地平線を生成する場合は、定常真空領域を通して積分を地平線へ移し、面積と結び付けられる。[ブラックホール熱力学](/physics/general-relativity/bh-thermodynamics/)

## ADM 質量と Bondi 質量

### 空間無限遠

:term[ADM エネルギー]{en="ADM energy" yomi="えーでぃーえむえねるぎー"} は漸近 Cartesian 座標で

$$
E_{\rm ADM}=(16\pi G)^{-1}\lim_{R\to\infty}\int_{S_R}(\partial_jh_{ij}-\partial_i h_{jj})\,s^i\dd A .
\label{eq-gr-adm-energy}
$$

線形化した切片のスカラー曲率は、接続の偏微分の差から $\partial_i\partial_j(h_{ij}-\delta_{ij})-\Delta(h_{jj}-3)$ である。先頭の体積積分を発散定理で境界へ移した式が上の組を与える。さらに ADM の Hamilton 拘束を $N\to1$ で変分すると同じ法線微分の境界項が出るため、これを相殺するようエネルギーを Hamiltonian に加える。外的曲率の符号は :eqref[eq-gr-adm-constraints] とそろえる。空間運動量は $P_i=(8\pi G)^{-1}\lim\int(Kh_{ij}-K_{ij})s^j\dd A$ である。これは shift の運動量拘束の境界項を、空間並進のエネルギー運動量の符号で表したものである。以下では座標変換の減衰を明示して、この条件でのエネルギーの不変性を示す。

Schwarzschild の等方半径 $r$ では $h_{ij}=(1+m/(2r))^4\delta_{ij}$ である。$h_{ij}=(1+2m/r)\delta_{ij}+\order{r^{-2}}$ を入れれば積分の法線成分は $4m/r^2$ となり、$E_{\rm ADM}=m/G$ である。静的切片では $K=0$ なので $\vb P=0$ である。

### 極限が存在する理由と座標を変えたときの一次項

$p_{ij}=h_{ij}-\delta_{ij}$、$V_i=\partial_jp_{ij}-\partial_ip_{jj}$ とする。曲率の式で $h^{-1}=\delta^{-1}+\order{r^{-1}}$ を使うと
$$
R^{(3)}=\partial_i V_i+\order{|p||\partial^2p|+|\partial p|^2}.
$$
誤差は $\order{r^{-4}}$ で積分可能である。従って $R^{(3)}$ が可積分なら、半径 $R_1,R_2$ の表面積分の差は、球殻内の $\partial_iV_i$ の積分として $R_1\to\infty$ で一様に零へ近づく。従って表面積分は Cauchy 条件を満たす。積分球の法線と面積は背景 Euclid 計量で評価できる。物理計量に替えても相対的な差が $\order{r^{-1}}$ で、$V_i=\order{r^{-2}}$ だから極限で消える。

例えば漸近座標変換を、固定の回転と並進を除いて $x^i=y^i+\zeta^i(y)$、$\zeta=\order{1}$、$\partial^k\zeta=\order{r^{-k}}$（$1\le k\le3$）に限定する。引き戻し計量の $1/r$ の変化は $\partial_i\zeta_j+\partial_j\zeta_i$ であり、$p$ の引数をずらした差や Jacobian の積は $\order{r^{-2}}$ である。積分核の変化は
$$
\Delta V_i=\Delta\zeta_i-\partial_i\partial_j\zeta_j+\order{r^{-3}} .
$$
先頭は反対称テンソル $\partial_j\zeta_i-\partial_i\zeta_j$ の発散である。$\zeta$ を積分球より内側で滑らかに零へ延ばせば、先頭の流束は $\partial_i\Delta\zeta_i-\Delta\partial_j\zeta_j=0$ から厳密に零となる。残りの流束は $\order{r^{-1}}$ に減るのでこの座標変換でエネルギーは不変である。回転では積分核と球面が一緒に回り、固定並進による球面の変更は可積分な発散の積分として消える。この証明では指定した落ち方を超えて変換を許していない。

### 切片の時間方向を変えるときのエネルギー

ADM の四元運動量について、未来向き漸近単位ベクトル $U^a$ に対応する測定エネルギーは $-P_aU^a$ である。例えば静止 Schwarzschild の漸近慣性枠で $P^a=(m/G,0,0,0)$ なら、速度 $v$ の枠で $E'=\gamma m/G$、$P'^x=-\gamma vm/G$ である。従って不変な ADM 質量は $\sqrt{E'^2-|\vb P'|^2}=m/G$ である。計量を同じ切片内で引き戻す座標変換ではこのような boost は生じない。漸近 Lorentz 変換で四元運動量として振る舞うという主張には、比較する切片を含む時空の無限遠の仮定が要る。静止切片の表面積分だけから任意の時間切片でのエネルギー保存を結論しない。

:::example[物質の運動量から無限遠の符号を読む]
$A^j{}_i=K^j{}_i-K\delta^j_i$ と置くと拘束は $D_j A^j{}_i=-8\pi G j_i$ である。外向き法線への $A^j{}_i$ の流束に負号を付けた量が $P_i$ なので、弱い場の最初の次数では $P_i=\int j_i\,\dd^3x$ になる。実際、$K$ と接続係数の積を落とすと $D_j$ が $\partial_j$ になり、発散定理を直接適用できる。$j_i=-T_{ni}$ は、未来向き $n=\partial_t$ の局所 Minkowski 枠で $T^{0i}$ に等しい。従って正の $x$ 方向に動く弱い物質源に対して、無限遠で定義した $P_x$ も正となる。
:::

### 零無限遠

:term[Bondi 質量]{en="Bondi mass" yomi="ぼんでぃしつりょう"} は遅延時刻 $u$ の断面ごとに定まる。面積半径 $r$ の Bondi 座標で $g_{uu}=-1+2m_B/r+\order{r^{-2}}$、$g_{AB}=r^2q_{AB}+rC_{AB}+\order{1}$、$q^{AB}C_{AB}=0$ とし、$N_{AB}=\partial_uC_{AB}$ を :term[news]{en="Bondi news" yomi="ぼんでぃにゅーす"} と呼ぶ。$q$ は単位球計量である。質量は $M_B=(4\pi G)^{-1}\int m_B\,\dd\Omega_2$ である。

滑らかな Bondi 展開を仮定し、$r^2T_{rr},r^2T_{rA}$ と物質の跡の $r^2T$ が零へ向かう場合を考える。以下で導く $uu$ 方程式は $\partial_um_B=\tfrac14D_AD_BN^{AB}-\tfrac18N_{AB}N^{AB}-4\pi G T_{uu}^{(2)}$ を与える。ここで $T_{uu}^{(2)}=\lim_{r\to\infty}r^2T_{uu}$ とする。球面で全発散は積分すると零なので

$$
\dv{M_B}{u}=-\frac1{32\pi G}\int N_{AB}N^{AB}\dd\Omega_2-\int T_{uu}^{(2)}\dd\Omega_2
\label{eq-gr-bondi-loss}
$$

である。放射だけなら正定値な球面計量で news の二乗が非負なので $M_B$ は非増加である。遠方の弱い波の $rh^{\rm TT}$ を角枠で $C$ と同定すれば、第一項の係数は :eqref[eq-gr-isaacson] の全流束と一致する。

### 面積座標の決定と $rr$ 方程式

Bondi 座標で $\gamma_{AB}=g_{AB}/r^2$、$g_{ur}=-e^{2b}$ とし、$\det\gamma=\det q$ を固定する。$k=\partial_r$ は出射光線に接し、$\nabla_k k=2(\partial_rb)k$ である。実際、$g_{rr}=g_{rA}=0$ のとき $\Gamma^a{}_{rr}=g^{au}\partial_rg_{ru}$ で、$g^{au}$ は $a=r$ でのみ非零である。面積密度が $r^2\sqrt{\det q}$ なので膨張率は $2/r$、せん断の混合成分は $\gamma^{AC}\partial_r\gamma_{CB}/2$ である。跡が零になるのは $\partial_r\det\gamma=0$ から従う。非アフィン光的 Raychaudhuri 方程式から
$$
R_{rr}=\frac{4}{r}\partial_r b
-\frac14\operatorname{tr}\left[(\gamma^{-1}\partial_r\gamma)^2\right].
$$
従って真空、または $T_{rr}=o(r^{-4})$ なら、$\gamma=q+C/r+\order{r^{-2}}$ を使い $\partial_r b=C_{AB}C^{AB}/(16r^3)+o(r^{-3})$ を得る。無限遠で $b\to0$ と規格化すれば $b=-C_{AB}C^{AB}/(32r^2)+o(r^{-2})$ である。従って $g_{ur}=-1+C_{AB}C^{AB}/(16r^2)+o(r^{-2})$ となる。$T_{rr}$ がこの次数を変える場合も、上で仮定した $g_{ur}$ の形が保たれれば、以下の $uu$ 方程式において係数 $B$ は相殺する。

行列式条件から角計量の次の跡も読める。$\gamma=q+C/r+D/r^2+\cdots$ で対数行列式を展開すると、$\tr(q^{-1}C)/r+[\tr(q^{-1}D)-\tr((q^{-1}C)^2)/2]/r^2$ が出る。従って $\tr(q^{-1}D)=C_{AB}C^{AB}/2$ である。この関係は角方向に出る見かけの跡の時間微分を、独立な放射自由度として数えないために使う。

### 質量面密度の変化を与える成分

Bondi 座標で $g_{rr}=g_{rA}=0$ とし、$\det(g_{AB}/r^2)=\det q$ で面積半径を固定する。$g_{ur}=-1+B/r^2+\order{r^{-3}}$、$g_{uA}=J_A+\order{r^{-1}}$ と置く。この落ち方を仮定して $rA$ の方程式を調べる。$\Gamma^B{}_{rA}=\delta^B_A/r-C^B{}_A/(2r^2)+\order{r^{-3}}$、$\Gamma^r{}_{rA}=J_A/r+\order{r^{-2}}$ である。$\partial_r\Gamma^r{}_{rA}$ と接続の積を合わせると $J_A/r^2$、$\Gamma^B{}_{rA}$ の角微分と角接続を合わせると $-D_B C^B{}_A/(2r^2)$ が残る。従って $r^{-2}$ で物質源がない場合に $J_A=D^BC_{AB}/2$ となる。

次に $\Gamma^a{}_{ua}=\partial_u\log\sqrt{-g}=-\partial_u B/r^2+\order{r^{-3}}$ を使う。$R_{uu}$ の微分部分で $\partial_u\Gamma^u{}_{uu}-\partial_u\Gamma^a{}_{ua}=-\partial_um_B/r^2+\order{r^{-3}}$、$\Gamma^r{}_{uu}=-(\partial_u m_B)/r+\order{r^{-2}}$ なのでその動径微分が同じ値を逆符号で与える。$\Gamma^a{}_{ar}=2/r+\order{r^{-3}}$ との積から $-2\partial_um_B/r^2$ が残る。角方向では $\Gamma^A{}_{uu}=(\partial_u J^A)/r^2+\order{r^{-3}}$ の発散と、$\Gamma^A{}_{uB}=N^A{}_B/(2r)+\order{r^{-2}}$ の積がそれぞれ $D_A\partial_u J^A/r^2$ と $-N_{AB}N^{AB}/(4r^2)$ を与える。従って跡項が $r^{-2}$ の $uu$ 成分に寄与しない場合、:eqref[eq-gr-trace-einstein] を適用すると上の質量面密度の式になる。$B$ は質量流束の式では相殺する。完全な漸近展開の存在や、他の座標で同じ量を定義するには別の解析が必要である。

:::example[球対称な放射模型で Bondi 質量を検算する]
外向きの Vaidya 計量を使い、$u$ を遅れ時刻とする。角微分が消える模型で $m(u)$ の非定常な項を追う。 座標の順序を $(u,r,\theta,\phi)$ と固定する。成分行列と逆行列は次の組になる。

$$
\begin{aligned}
(g_{ab})&=\left[\begin{matrix}-1 + \frac{2 m{\left(u \right)}}{r} & -1 & 0 & 0\\-1 & 0 & 0 & 0\\0 & 0 & r^{2} & 0\\0 & 0 & 0 & r^{2} \sin^{2}{\left(\theta \right)}\end{matrix}\right]\\
(g^{ab})&=\left[\begin{matrix}0 & -1 & 0 & 0\\-1 & \frac{r - 2 m{\left(u \right)}}{r} & 0 & 0\\0 & 0 & \frac{1}{r^{2}} & 0\\0 & 0 & 0 & \frac{1}{r^{2} \sin^{2}{\left(\theta \right)}}\end{matrix}\right]\\
\det g&=- r^{4} \sin^{2}{\left(\theta \right)}
\end{aligned}
$$

逆行列を掛けて単位行列を得ることを確かめてから、:eqref[eq-gr-christoffel] を使う。同じ下添字の順序を二度計算せず、$b\le c$ の組を選ぶ。中央の列には、座標に依存する計量成分の微分だけを残した。最後の列は微分を実行した値であり、表に出ない成分は添字交換で得られるものを除いて零である。

| 成分 | 微分への代入 | 結果 |
|---|---|---|
| $\Gamma^{u}_{u u}$ | $- \frac{1}{2}\left(-\partial_{r}\left(-1 + \frac{2 m{\left(u \right)}}{r}\right)\right)$ | $- \frac{m{\left(u \right)}}{r^{2}}$ |
| $\Gamma^{u}_{\theta \theta}$ | $- \frac{1}{2}\left(-\partial_{r}\left(r^{2}\right)\right)$ | $r$ |
| $\Gamma^{u}_{\phi \phi}$ | $- \frac{1}{2}\left(-\partial_{r}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $r \sin^{2}{\left(\theta \right)}$ |
| $\Gamma^{r}_{u u}$ | $- \frac{1}{2}\left(\partial_{u}\left(-1 + \frac{2 m{\left(u \right)}}{r}\right)+\partial_{u}\left(-1 + \frac{2 m{\left(u \right)}}{r}\right)-\partial_{u}\left(-1 + \frac{2 m{\left(u \right)}}{r}\right)\right)+\frac{r - 2 m{\left(u \right)}}{2 r}\left(-\partial_{r}\left(-1 + \frac{2 m{\left(u \right)}}{r}\right)\right)$ | $- \frac{r^{2} m'{\left(u \right)} - r m{\left(u \right)} + 2 m^{2}{\left(u \right)}}{r^{3}}$ |
| $\Gamma^{r}_{u r}$ | $- \frac{1}{2}\left(\partial_{r}\left(-1 + \frac{2 m{\left(u \right)}}{r}\right)\right)$ | $\frac{m{\left(u \right)}}{r^{2}}$ |
| $\Gamma^{r}_{\theta \theta}$ | $\frac{r - 2 m{\left(u \right)}}{2 r}\left(-\partial_{r}\left(r^{2}\right)\right)$ | $- r + 2 m{\left(u \right)}$ |
| $\Gamma^{r}_{\phi \phi}$ | $\frac{r - 2 m{\left(u \right)}}{2 r}\left(-\partial_{r}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $- \left(r - 2 m{\left(u \right)}\right) \sin^{2}{\left(\theta \right)}$ |
| $\Gamma^{\theta}_{r \theta}$ | $\frac{1}{2 r^{2}}\left(\partial_{r}\left(r^{2}\right)\right)$ | $\frac{1}{r}$ |
| $\Gamma^{\theta}_{\phi \phi}$ | $\frac{1}{2 r^{2}}\left(-\partial_{\theta}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $- \frac{\sin{\left(2 \theta \right)}}{2}$ |
| $\Gamma^{\phi}_{r \phi}$ | $\frac{1}{2 r^{2} \sin^{2}{\left(\theta \right)}}\left(\partial_{r}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{1}{r}$ |
| $\Gamma^{\phi}_{\theta \phi}$ | $\frac{1}{2 r^{2} \sin^{2}{\left(\theta \right)}}\left(\partial_{\theta}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{1}{\tan{\left(\theta \right)}}$ |

曲率の計算では :eqref[eq-gr-curvature-components] の微分差と積差を分ける。第二列は $\partial_c\Gamma^a{}_{db}$、第三列は $-\partial_d\Gamma^a{}_{cb}$ であり、第四列は接続の積を符号付きで足したものである。三列を足して最後の列を得る。

| 成分 | 第一の微分 | 第二の微分 | 接続の積 | 結果 |
|---|---|---|---|---|
| $R^{u}_{\theta u \theta}$ | $0$ | $0$ | $- \frac{m{\left(u \right)}}{r}$ | $- \frac{m{\left(u \right)}}{r}$ |
| $R^{u}_{\phi u \phi}$ | $0$ | $0$ | $- \frac{m{\left(u \right)} \sin^{2}{\left(\theta \right)}}{r}$ | $- \frac{m{\left(u \right)} \sin^{2}{\left(\theta \right)}}{r}$ |
| $R^{r}_{r u r}$ | $0$ | $\frac{2 m{\left(u \right)}}{r^{3}}$ | $0$ | $\frac{2 m{\left(u \right)}}{r^{3}}$ |
| $R^{r}_{\theta u \theta}$ | $2 m'{\left(u \right)}$ | $0$ | $\frac{\left(- r + 2 m{\left(u \right)}\right) m{\left(u \right)}}{r^{2}} - \frac{r^{2} m'{\left(u \right)} - r m{\left(u \right)} + 2 m^{2}{\left(u \right)}}{r^{2}}$ | $m'{\left(u \right)}$ |
| $R^{r}_{\theta r \theta}$ | $-1$ | $0$ | $\frac{r - 2 m{\left(u \right)}}{r} + \frac{m{\left(u \right)}}{r}$ | $- \frac{m{\left(u \right)}}{r}$ |
| $R^{r}_{\phi u \phi}$ | $2 \sin^{2}{\left(\theta \right)} m'{\left(u \right)}$ | $0$ | $- \frac{\left(r - 2 m{\left(u \right)}\right) m{\left(u \right)} \sin^{2}{\left(\theta \right)}}{r^{2}} - \frac{\left(r^{2} m'{\left(u \right)} - r m{\left(u \right)} + 2 m^{2}{\left(u \right)}\right) \sin^{2}{\left(\theta \right)}}{r^{2}}$ | $\sin^{2}{\left(\theta \right)} m'{\left(u \right)}$ |
| $R^{r}_{\phi r \phi}$ | $- \sin^{2}{\left(\theta \right)}$ | $0$ | $\frac{\left(r - 2 m{\left(u \right)}\right) \sin^{2}{\left(\theta \right)}}{r} + \frac{m{\left(u \right)} \sin^{2}{\left(\theta \right)}}{r}$ | $- \frac{m{\left(u \right)} \sin^{2}{\left(\theta \right)}}{r}$ |
| $R^{\theta}_{\phi \theta \phi}$ | $- \cos{\left(2 \theta \right)}$ | $0$ | $\frac{\sin{\left(2 \theta \right)}}{2 \tan{\left(\theta \right)}} - \frac{\left(r - 2 m{\left(u \right)}\right) \sin^{2}{\left(\theta \right)}}{r}$ | $\frac{2 m{\left(u \right)} \sin^{2}{\left(\theta \right)}}{r}$ |

第一添字を計量で下ろすと $R_{abcd}=g_{ak}R^k{}_{bcd}$ となる。対の反対称性と交換対称性を用いて、非零の独立な組を下表に示す。

| 成分 | 結果 |
|---|---|
| $R_{u r u r}$ | $- \frac{2 m{\left(u \right)}}{r^{3}}$ |
| $R_{u \theta u \theta}$ | $- \frac{r^{2} m'{\left(u \right)} - r m{\left(u \right)} + 2 m^{2}{\left(u \right)}}{r^{2}}$ |
| $R_{u \theta r \theta}$ | $\frac{m{\left(u \right)}}{r}$ |
| $R_{u \phi u \phi}$ | $- \frac{\left(r^{2} m'{\left(u \right)} - r m{\left(u \right)} + 2 m^{2}{\left(u \right)}\right) \sin^{2}{\left(\theta \right)}}{r^{2}}$ |
| $R_{u \phi r \phi}$ | $\frac{m{\left(u \right)} \sin^{2}{\left(\theta \right)}}{r}$ |
| $R_{\theta \phi \theta \phi}$ | $2 r m{\left(u \right)} \sin^{2}{\left(\theta \right)}$ |

Ricci テンソルは第一・第三添字を縮約して得る。曲率が残る場合も相殺が起こるので、各項を落とさずに加える。

| 成分 | 縮約の項 | 代入 | 結果 |
|---|---|---|---|
| $R_{u u}$ | $R^{r}_{u r u}+R^{\theta}_{u \theta u}+R^{\phi}_{u \phi u}$ | $- \frac{2 \left(r - 2 m{\left(u \right)}\right) m{\left(u \right)}}{r^{4}} - \frac{r^{2} m'{\left(u \right)} - r m{\left(u \right)} + 2 m^{2}{\left(u \right)}}{r^{4}} - \frac{r^{2} m'{\left(u \right)} - r m{\left(u \right)} + 2 m^{2}{\left(u \right)}}{r^{4}}$ | $- \frac{2 m'{\left(u \right)}}{r^{2}}$ |
| $R_{u r}$ | $R^{u}_{u u r}+R^{\theta}_{u \theta r}+R^{\phi}_{u \phi r}$ | $- \frac{2 m{\left(u \right)}}{r^{3}} + \frac{m{\left(u \right)}}{r^{3}} + \frac{m{\left(u \right)}}{r^{3}}$ | $0$ |
| $R_{\theta \theta}$ | $R^{u}_{\theta u \theta}+R^{r}_{\theta r \theta}+R^{\phi}_{\theta \phi \theta}$ | $- \frac{m{\left(u \right)}}{r} - \frac{m{\left(u \right)}}{r} + \frac{2 m{\left(u \right)}}{r}$ | $0$ |
| $R_{\phi \phi}$ | $R^{u}_{\phi u \phi}+R^{r}_{\phi r \phi}+R^{\theta}_{\phi \theta \phi}$ | $- \frac{m{\left(u \right)} \sin^{2}{\left(\theta \right)}}{r} - \frac{m{\left(u \right)} \sin^{2}{\left(\theta \right)}}{r} + \frac{2 m{\left(u \right)} \sin^{2}{\left(\theta \right)}}{r}$ | $0$ |

逆計量と縮約すると Ricci スカラーとなる。さらに $G_{ab}=R_{ab}-g_{ab}R/2$ を成分ごとに計算する。

$$
\begin{aligned}
R&=\sum_{a,b}g^{ab}R_{ab}\\
&=0\\
&=0
\end{aligned}
$$

非零の Einstein テンソル成分は次のとおりである。

| 成分 | Ricci 項 | 跡の項 | 差 |
|---|---|---|---|
| $G_{u u}$ | $- \frac{2 m'{\left(u \right)}}{r^{2}}$ | $0$ | $- \frac{2 m'{\left(u \right)}}{r^{2}}$ |

$G_{uu}=-2m'(u)/r^2$ なので、外向きの正のエネルギーで $m$ が減る。そこで $M_B=m/G$ と $T^{(2)}_{uu}=r^2T_{uu}=-m'(u)/(4\pi G)$ を :eqref[eq-gr-bondi-loss] に入れると、球面積分が $m'/G$ を与える。球面の全発散と news が零であるため、質量変化は物質の流束だけに由来する。
:::

:::example[等方座標で空間無限遠の係数を出す]
ここでの $r$ は面積半径ではなく等方座標の半径で、$r>m/2$ とする。$\psi=1+m/(2r)$、$N=(1-m/(2r))/\psi$ を使う。面積半径は $r\psi^2$ である。 座標の順序を $(t,r,\theta,\phi)$ と固定する。成分行列と逆行列は次の組になる。

$$
\begin{aligned}
(g_{ab})&=\left[\begin{matrix}- \frac{\left(- \frac{m}{2 r} + 1\right)^{2}}{\left(\frac{m}{2 r} + 1\right)^{2}} & 0 & 0 & 0\\0 & \left(\frac{m}{2 r} + 1\right)^{4} & 0 & 0\\0 & 0 & r^{2} \left(\frac{m}{2 r} + 1\right)^{4} & 0\\0 & 0 & 0 & r^{2} \left(\frac{m}{2 r} + 1\right)^{4} \sin^{2}{\left(\theta \right)}\end{matrix}\right]\\
(g^{ab})&=\left[\begin{matrix}- \frac{\left(m + 2 r\right)^{2}}{\left(- m + 2 r\right)^{2}} & 0 & 0 & 0\\0 & \frac{16 r^{4}}{\left(m + 2 r\right)^{4}} & 0 & 0\\0 & 0 & \frac{16 r^{2}}{\left(m + 2 r\right)^{4}} & 0\\0 & 0 & 0 & \frac{16 r^{2}}{\left(m + 2 r\right)^{4} \sin^{2}{\left(\theta \right)}}\end{matrix}\right]\\
\det g&=- \frac{\left(- m + 2 r\right)^{2} \left(m + 2 r\right)^{10} \sin^{2}{\left(\theta \right)}}{4096 r^{8}}
\end{aligned}
$$

逆行列を掛けて単位行列を得ることを確かめてから、:eqref[eq-gr-christoffel] を使う。同じ下添字の順序を二度計算せず、$b\le c$ の組を選ぶ。中央の列には、座標に依存する計量成分の微分だけを残した。最後の列は微分を実行した値であり、表に出ない成分は添字交換で得られるものを除いて零である。

| 成分 | 微分への代入 | 結果 |
|---|---|---|
| $\Gamma^{t}_{t r}$ | $- \frac{\left(m + 2 r\right)^{2}}{2 \left(- m + 2 r\right)^{2}}\left(\partial_{r}\left(- \frac{\left(- \frac{m}{2 r} + 1\right)^{2}}{\left(\frac{m}{2 r} + 1\right)^{2}}\right)\right)$ | $\frac{4 m}{\left(- m + 2 r\right) \left(m + 2 r\right)}$ |
| $\Gamma^{r}_{t t}$ | $\frac{8 r^{4}}{\left(m + 2 r\right)^{4}}\left(-\partial_{r}\left(- \frac{\left(- \frac{m}{2 r} + 1\right)^{2}}{\left(\frac{m}{2 r} + 1\right)^{2}}\right)\right)$ | $\frac{64 m r^{4} \left(- m + 2 r\right)}{\left(m + 2 r\right)^{7}}$ |
| $\Gamma^{r}_{r r}$ | $\frac{8 r^{4}}{\left(m + 2 r\right)^{4}}\left(\partial_{r}\left(\left(\frac{m}{2 r} + 1\right)^{4}\right)+\partial_{r}\left(\left(\frac{m}{2 r} + 1\right)^{4}\right)-\partial_{r}\left(\left(\frac{m}{2 r} + 1\right)^{4}\right)\right)$ | $- \frac{2 m}{r \left(m + 2 r\right)}$ |
| $\Gamma^{r}_{\theta \theta}$ | $\frac{8 r^{4}}{\left(m + 2 r\right)^{4}}\left(-\partial_{r}\left(r^{2} \left(\frac{m}{2 r} + 1\right)^{4}\right)\right)$ | $- \frac{r \left(- m + 2 r\right)}{m + 2 r}$ |
| $\Gamma^{r}_{\phi \phi}$ | $\frac{8 r^{4}}{\left(m + 2 r\right)^{4}}\left(-\partial_{r}\left(r^{2} \left(\frac{m}{2 r} + 1\right)^{4} \sin^{2}{\left(\theta \right)}\right)\right)$ | $- \frac{r \left(- m + 2 r\right) \sin^{2}{\left(\theta \right)}}{m + 2 r}$ |
| $\Gamma^{\theta}_{r \theta}$ | $\frac{8 r^{2}}{\left(m + 2 r\right)^{4}}\left(\partial_{r}\left(r^{2} \left(\frac{m}{2 r} + 1\right)^{4}\right)\right)$ | $\frac{- m + 2 r}{r \left(m + 2 r\right)}$ |
| $\Gamma^{\theta}_{\phi \phi}$ | $\frac{8 r^{2}}{\left(m + 2 r\right)^{4}}\left(-\partial_{\theta}\left(r^{2} \left(\frac{m}{2 r} + 1\right)^{4} \sin^{2}{\left(\theta \right)}\right)\right)$ | $- \frac{\sin{\left(2 \theta \right)}}{2}$ |
| $\Gamma^{\phi}_{r \phi}$ | $\frac{8 r^{2}}{\left(m + 2 r\right)^{4} \sin^{2}{\left(\theta \right)}}\left(\partial_{r}\left(r^{2} \left(\frac{m}{2 r} + 1\right)^{4} \sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{- m + 2 r}{r \left(m + 2 r\right)}$ |
| $\Gamma^{\phi}_{\theta \phi}$ | $\frac{8 r^{2}}{\left(m + 2 r\right)^{4} \sin^{2}{\left(\theta \right)}}\left(\partial_{\theta}\left(r^{2} \left(\frac{m}{2 r} + 1\right)^{4} \sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{1}{\tan{\left(\theta \right)}}$ |

空間部分を直交座標に移すと $h_{ij}=\psi^4\delta_{ij}=(1+2m/r+\order{r^{-2}})\delta_{ij}$ である。積分核の各項は $\partial_jh_{ij}=\partial_i\psi^4$、$\partial_i h_{jj}=3\partial_i\psi^4$ なので、球の法線と縮約すると $-2\partial_r\psi^4=4m/r^2+\order{r^{-3}}$ になる。積分面積の先頭を $4\pi r^2$ と取れば ADM エネルギーは $m/G$ である。法線と面積の補正はともに相対的な $O(r^{-1})$ なので、極限には残らない。
:::

### 適切なゲージでの ADM 積分の計算例

球対称な時間対称切片で面積半径 $r$ を使い、$h=f(r)^{-1}\dd r^2+r^2\dd\Omega^2$ とする。漸近 Cartesian 座標では $h_{ij}=\delta_{ij}+b(r)n_i n_j$、$b=f^{-1}-1$、$n_i=x_i/r$ である。$\partial_j n_i=(\delta_{ij}-n_i n_j)/r$ を用いると
$$
\partial_j h_{ij}=b'n_i+\frac{2b}{r}n_i,\qquad
\partial_i h_{jj}=b'n_i .
$$
従って ADM 積分の結果は $\lim_{r\to\infty}rb(r)/(2G)$ である。$f=1-2m/r$ なら $rb\to2m$ なので $m/G$ になる。等方座標の $2m\,\delta_{ij}/r$ とは見た目が異なるが、同じ積分値を与える。面積半径表示で有限の $r$ において $rb/(2G)$ を計算しても、その値は一般には有限半径の物理的質量ではない。幾何的な球面の質量には $r(1-f)/(2G)$ を使い、二つの値が一致するのは上の減衰条件での極限である。

### Bondi エネルギーと運動量の流束

球面上の $n_i=x_i/r$ に対し $D_AD_B n_i=-q_{AB}n_i$ である。Bondi の空間運動量を $P_{B,i}=(4\pi G)^{-1}\int m_B n_i\,\dd\Omega_2$ と定める。$uu$ で得た式に $n_i$ を掛け、二回部分積分を行うと、全発散の寄与は $\int (D_AD_B n_i)N^{AB}=0$ である。従って一つの式にまとめて
$$
\frac{\dd P_B^\mu}{\dd u}
=-\int (1,\vb n)^\mu\left(\frac{N_{AB}N^{AB}}{32\pi G}
+T^{(2)}_{uu}\right)\dd\Omega_2 .
$$
右辺の括弧が非負なら、失う四元運動量は未来因果的または零である。実際、各方向が光的な $(1,\vb n)$ を与え、空間成分の大きさを三角不等式で評価すれば時間成分以下となる。例えば南北が対称な流束では $n_z$ の積分が相殺し、エネルギーだけが減る。非対称な場合には反作用として系の線形運動量も変わる。

ここで $M_B$ と書いた全角平均は、選んだ Bondi 枠でのエネルギー成分 $P_B^0$ である。空間成分が非零なら静止質量は $({P_B^0}^2-|\vb P_B|^2)^{1/2}$ であり、枠を変える場合には区別する。系が初期に放射を持たない定常外部に近く、空間無限遠と $\mathscr I^+$ の過去端で座標と電荷を対応させられるなら、初期の Bondi エネルギーを ADM エネルギーと照合できる。一般の時空でこの同一視をするには、$i^0$ 近傍の適合条件が要る。

:::example[運動量を持ち去る方向を積分する]
ある断面の全流束密度を $F(\theta)=F_0(1+\alpha\cos\theta)$、$F_0>0,|\alpha|\le1$ とする。角積分で $\dot P_B^0=-4\pi F_0$、$\dot P_B^z=-4\pi F_0\alpha/3$ となる。実際 $\int\cos\theta\dd\Omega_2=0$、$\int\cos^2\theta\dd\Omega_2=4\pi/3$ である。$\alpha>0$ では北向きの放射が多いので、残る系の運動量は負方向へ変化する。運び出される運動量の大きさは、運び出されるエネルギーの $|\alpha|/3$ 倍である。
:::

## 正エネルギー定理

:term[正エネルギー定理]{en="positive energy theorem" yomi="せいえねるぎーていり"} は、完全で向き付け可能な境界のない漸近平坦な三次元初期データが Einstein の拘束と $\rho\ge\abs{\vb j}_h$ を満たすとき、適切な減衰と可積分性の下で各端の ADM 四元運動量が未来因果的または零、すなわち $E_{\rm ADM}\ge\abs{\vb P}$ であると述べる。零エネルギーの場合の剛性はデータが Minkowski 時空の切片として実現されることを意味する。スピノル法では漸近定数の境界値を持つ Dirac 型方程式を解き、境界の四元運動量を非負なバルク積分として表す。非コンパクト多様体上の楕円型方程式の可解性と剛性の証明はこの章の範囲を超えるため :cite[wald-gr] を参照する。例えば形式的に $m<0$ とした Schwarzschild 真空外部で $E_{\rm ADM}<0$ となっても、中心の特異点で空間の完全性が破れており、反例ではない。

### 正値性を積分で確認できる共形データ

時間対称データ $h_{ij}=\psi^4\delta_{ij}$ を全 $\R^3$ で滑らかとし、$\psi=1+c/r+\order{r^{-2}}$、一階微分も対応して減衰するとする。拘束は $-8\psi^{-5}\Delta\psi=16\pi G\rho$、従って $\Delta\psi=-2\pi G\rho\psi^5$ である。球で積分してから $r\to\infty$ とすると $-4\pi c=-2\pi G\int\rho\psi^5\dd^3x$ である。ADM 式で $h_{ij}\simeq(1+4c/r)\delta_{ij}$ を使うと $E_{\rm ADM}=2c/G$ となり、従って $\int\rho\psi^5\dd^3x$ と等しい。$\rho\ge0$ ならこの特殊な場合の正値性を直接得る。$\rho=0$ なら $\psi-1$ は無限遠で零の調和関数なので、球での最大値原理から $\psi=1$ である。この証明には共形平坦性・時間対称性と内部境界がないことを用いた。

### 内部境界が正値性の仮定を変える例

:::example[負の質量を持つ形式解の端]
$\psi=1-b/r$、$b>0$ とし、$r>b$ だけを持つ時間対称なデータを考える。この領域では $\Delta\psi=0$、先頭の ADM エネルギーは $-2b/G$ である。だが $r\downarrow b$ に近づく固有距離は $\int_b^R (1-b/r)^2\dd r$ で、被積分は零に近づくため有限である。$r=b$ は計量が退化して正定値の内点としては含まれない。従ってこのデータは完全性を仮定する正エネルギー定理から外れる。負の符号を生むのは、原点まで滑らかな全空間の発散定理を使えないことに現れる。
:::

### 時間対称な場合のスピノルによる境界評価

時間対称性だけを仮定し、共形平坦性は課さない。向き付けられた三次元空間のスピン構造を取り、Clifford 乗法を $c_i c_j+c_j c_i=-2\delta_{ij}$、$c_i^\dagger=-c_i$ に規格化する。ただし $c_i=c(e_i)$ である。$\omega_{iab}=g(e_a,D_{e_i}e_b)$ とすると、スピン接続はこの枠で $\nabla^s_{e_i}=e_i-\omega_{iab}c_a c_b/4$ である。符号は $[-\omega_{iab}c_a c_b/4,c_j]=\omega_{iaj}c_a$ として枠の接続と照合できる。その曲率は $[\nabla^s_{e_i},\nabla^s_{e_j}]-\nabla^s_{[e_i,e_j]}=-R_{abij}c_a c_b/4$ である。$D_s=\sum_i c_i\nabla^s_{e_i}$ の二乗を正規直交枠で展開すると、同じ方向の項の和は $\nabla^{s*}\nabla^s$ になる。異なる方向の項は $-\sum_{ijab}R_{abij}c_i c_j c_a c_b/8$ になる。この和で四添字が相異なる部分は Bianchi 恒等式で消え、一組だけが一致する部分は対称な Ricci テンソルと反対称な Clifford 積の縮約で消える。残る $i=a,j=b$ と $i=b,j=a$ ではそれぞれ $c_i c_j c_i c_j=-1$、$c_i c_j c_j c_i=1$ なので $R^{(3)}/4$ が得られる。従って
$$
D_s^2=\nabla^{s*}\nabla^s+R^{(3)}/4 .
$$
有限領域で部分積分すると、$\int(|\nabla^s\psi|^2+R^{(3)}|\psi|^2/4-|D_s\psi|^2)$ は境界の $\operatorname{Re}\int\langle\psi,(\nabla_s^s+c(s)D_s)\psi\rangle$ と等しい。ここで $s$ は切片内の外向き単位法線である。$\nabla^s$ と $D_s$ それぞれに Green の公式を用いればこの組を得る。

$D_s\psi=0$ を解き、対象の端で $\psi\to\psi_\infty$、他の端で $\psi\to0$ を課す。漸近端では $h_{ij}=\delta_{ij}+p_{ij}$、$e_i=(\delta_i^j-p_i{}^j/2)\partial_j+\order{r^{-2}}$ の枠で定数スピノルを表す。境界積分に $\nabla_s^s+c(s)D_s$ を用いると、$\psi-\psi_\infty$ の一次微分は Clifford 二乗の対角成分が打ち消され、反対称な成分の発散の流束となるため積分で消える。一方、接続の一次項は
$$
\omega_{iab}=g(e_a,\nabla_{e_i}e_b)
=\tfrac12(\partial_b p_{ia}-\partial_a p_{ib})+\order{r^{-3}}
$$
であり、実部を取った Clifford 積のスカラー成分は $(\partial_jp_{ij}-\partial_ip_{jj})s^i|\psi_\infty|^2/4$ である。従って極限の境界積分は $4\pi G E_{\rm ADM}|\psi_\infty|^2$ となる。Hamilton 拘束が $R^{(3)}\ge0$ を与えるので、バルクの積分は非負である。$E_{\rm ADM}=0$ なら $\nabla^s\psi=0$ であり、非零な平行スピノルが存在する。曲率交換子を Clifford 乗法で縮約すると $\operatorname{Ric}_{ij}c(e_j)\psi=0$ を得る。実ベクトル $v$ の Clifford 乗法の二乗は $-|v|^2$ なので、$\operatorname{Ric}=0$ となる。三次元では Ricci が Riemann を決めるため局所的に平坦である。

この導出で未証明なのは、指定した端の境界条件を満たすスピノルの存在と減衰、一般の $K$ を持つ場合の時空接続による評価、および大域的な剛性である。これらには漸近端での楕円型作用素と拘束方程式を組み合わせた評価が必要となる。定理の証明は :cite[wald-gr] を参照する。

## 漸近対称性と BMS 群の概観

:term[BMS 群]{en="Bondi–Metzner–Sachs group" yomi="ぼんでぃめつなーさっくすぐん"} は漸近平坦な零無限遠の普遍構造を保つ変換である。$\mathscr I^+\simeq\R_u\times S^2$ の退化計量 $q$ と零生成子 $n=\partial_u$ を、共形同値 $(q,n)\sim(\omega^2 q,\omega^{-1}n)$ の組として扱う。ただし $\omega$ は角だけに依存する。角写像 $x'\!=F(x)$ は $F^*q=K^2q$ の球の共形変換であり、$u'=K(x)(u+\alpha(x))$ とすれば $q$ と $n$ が相応に変換する。逆に生成子の向きと共形構造を保つには角度が $u$ に依存せず、$u'$ の傾きはこの $K$ になるので、積分定数の $\alpha(x)$ が残る。滑らかな大域共形変換の向き保存成分は Lorentz 群に対応し、任意の滑らかな $\alpha$ が :term[超並進]{en="supertranslation" yomi="ちょうへいしん"} を与える。従って群は超並進と Lorentz 群の半直積である。通常の並進は $\alpha$ の $\ell=0,1$ 成分である。例えば平坦時空で $t\mapsto t+\delta t$、$\vb r\mapsto\vb r+\vb d$ による $u$ の先頭変化は $\delta t-\hat{\vb r}\cdot\vb d$ となる。$\ell\ge2$ の $\alpha$ は方向ごとに遅延時刻の原点を別々に変える。

### 超並進と波の記憶

弱い波で $\ddot h$ が有限期間だけ非零であり、前後で $\dot h=0$ とする。距離変化の式を波の前後で比較すると $\Delta L/L=\tfrac12 n^i n^j\Delta h^{\rm TT}_{ij}$ となる。遠方では $\Delta h^{\rm TT}_{AB}\simeq\Delta C_{AB}/r$ である。従って $N=\partial_u C$ を時刻で積分した $\Delta C_{AB}$ が、波の通過後に残る相対変位を表す。これを :term[重力波メモリー]{en="gravitational wave memory" yomi="じゅうりょくはめもりー"} という。

同じ $C_{AB}$ は零無限遠の座標の取り方でも変わる。例えば Minkowski の外向き座標で微小な変換ベクトルを $\xi^u=\alpha$、$\xi^A=-D^A\alpha/r$、$\xi^r=D^2\alpha/2$ とする。$rr$ と $rA$ の変化は零であり、最後の成分は角計量の行列式を変えないよう選んだ。$g_{AB}=r^2q_{AB}$ の Lie 微分を取ると、$r^{-1}$ の角変位から $-2rD_AD_B\alpha$、動径の変位から $r q_{AB}D^2\alpha$ が現れる。従って同じ $\alpha$ の変換は二つの定常断面の $C$ の双方に同じ量を足す。二断面の差はこの原点の変更で消えず、光を用いた実験の永続変位と対応する。

### 球面調和関数で並進の四成分を取り出す

超並進 $f$ が剪断に加える角テンソル $-2D_AD_Bf+q_{AB}D^2 f$ に着目する。定数は微分が零なので時間並進に対応する。単位球面の座標関数 $n_i$ は $D_AD_B n_i=-q_{AB}n_i$ を満たす。実際、埋込み位置の二階微分の接成分を球面接続で引くと、法線方向の成分が $-q_{AB}\vb n$ になる。従って $D^2n_i=-2n_i$ で無跡 Hessian も零であり、三空間並進に対応する。例えば $f=P_2(\cos\theta)$ では $D_\theta D_\theta f-\frac12D^2f=3\sin^2\theta/2$ となり、無跡部が残る。球面で任意の関数を選べることで、四つの通常の並進以外の変換も定義される。

### 超並進を残した Lie 括弧

未来光的無限遠で $u$ を生成子の座標、$x^A$ を角座標とする。球面の共形 Killing 場を $Y^A$、$\psi_Y=D_AY^A$ とし、無限遠での生成子を
$$
X_{(f,Y)}=(f+u\psi_Y/2)\partial_u+Y^A\partial_A
$$
と書く。二つの括弧の角成分は $[Y,Z]$ である。$u$ 成分では $u\psi_Y\psi_Z$ の積が交換して消え、$u$ に依存しない部分は
$$
h=Y^A\partial_Ag-Z^A\partial_Af+\tfrac12(f\psi_Z-g\psi_Y)
$$
になる。$u$ に比例する部分は $u(Y\psi_Z-Z\psi_Y)/2=uD_A[Y,Z]^A/2$ である。最後は発散を体積形式の Lie 微分で表し、$[\mathcal L_Y,\mathcal L_Z]=\mathcal L_{[Y,Z]}$ を使えば従う。$Y=Z=0$ なら括弧は零なので超並進は可換である。一方 $Y$ が回転なら $\psi_Y=0$ なので、角度依存の超並進は $f$ の角微分によって回転の作用を受ける。

### Bondi ゲージを保つベクトル場を探す

平坦背景を出射座標で $g=-\dd u^2-2\dd u\dd r+r^2q_{AB}\dd x^A\dd x^B$ と書く。$\xi$ による変換で $g_{rr}=g_{rA}=0$ を保つと、順に $\partial_r\xi^u=0$、$r^2q_{AB}\partial_r\xi^B=\partial_A\xi^u$ である。従って $\xi^u=f(u,x)$、$\xi^A=Y^A(u,x)-r^{-1}D^Af$ となる。面積座標の条件を保つため $(\mathcal L_\xi g_{AB})q^{AB}=0$ を課すと
$$
\xi^r=-\tfrac r2D_AY^A+\tfrac12 D^2 f
$$
である。角成分の $r^2$ の係数の無跡部が零であるためには $D_{(A}Y_{B)}=\tfrac12 q_{AB}D_CY^C$ を要する。$uA$ 成分の $r^2$ の項を消すには $\partial_uY^A=0$、$ur$ 成分の先頭を保つには $-\partial_uf-\partial_r\xi^r=0$ であり、$f=f_0(x)+u D_AY^A/2$ を得る。$Y=0$ の部分では角成分の $r$ の係数が $q_{AB}D^2 f_0-2D_AD_B f_0$ である。従って背景の剪断が零でも、角依存の超並進は新しい剪断の成分を作る。時間に依存しないこの項は news の定義 $\partial_u C_{AB}$ では消える。

### 角依存の電荷とメモリーの制約

任意の滑らかな実関数 $f(x^A)$ について $Q_f(u)=(4\pi G)^{-1}\int f m_B\,\dd\Omega_2$ と定める。ここでは固定した Bondi 枠の一連の断面で、この定義による変化だけを求める。$uu$ 方程式に $f$ を掛けて部分積分すると
$$
\frac{\dd Q_f}{\dd u}
=\frac{1}{16\pi G}\int (D_AD_Bf)N^{AB}\dd\Omega_2
-\int f\left(\frac{N_{AB}N^{AB}}{32\pi G}+T_{uu}^{(2)}\right)\dd\Omega_2.
$$
並進に対応する $f=1,n_i$ では第一項が消えるが、それより高い角次数では一般に残る。また $f$ が負の部分を持てば最後の積分にも一定の符号はない。従って全エネルギーの単調性だけから、これらの角度で重みづけした量の単調性を結論できない。

メモリーのうち、球面上のスカラー $F$ を使って $\Delta C_{AB}=2D_AD_B F-q_{AB}D^2F$ と表せる部分を考える。単位球面の Ricci テンソルは $q_{AB}$ なので、共変微分の交換により $D^B D_A D_B F=D_A D^2F+D_A F$ である。従って $D^B\Delta C_{AB}=D_A(D^2+2)F$、さらに二重発散は $D^2(D^2+2)F$ となる。$uu$ 方程式を $u_1$ から $u_2$ まで積分して
$$
D^2(D^2+2)F
=4\Delta m_B
+\int_{u_1}^{u_2}\left(\tfrac12 N_{AB}N^{AB}+16\pi G T_{uu}^{(2)}\right)\dd u
$$
を得る。ただし $\Delta C$ のうち上の形以外の成分がある場合、この式で決まるのは二重発散に寄与する成分である。例えば $2\epsilon_{C(A}D_{B)}D^CG$ の二重発散は零である。第一発散は $D^2G+2G$ の回転勾配に比例し、第二発散は Hessian の対称性と面積形式の反対称性で消える。各次数 $\ell\ge2$ に対し、逆演算の固有値は $[\ell(\ell+1)(\ell(\ell+1)-2)]^{-1}$ である。$\ell=0,1$ では左辺が零なので、右辺は先に求めたエネルギー運動量の収支を満たさなければならない。ここに $N^2$ が現れるため、重力波そのものが運ぶエネルギーによってもメモリーが生じ得る。

:::example[四重極型の news を一つ積分する]
$f=P_2(\cos\theta)$、$C_{AB}=b(u)(-2D_AD_B f+q_{AB}D^2 f)$ とする。$\theta$ 方向のテンソル成分は $-3b\sin^2\theta$、正規化した $\phi$ 方向では $+3b\sin^2\theta$ であり、混合成分は零である。従って $N_{AB}N^{AB}=18(b')^2\sin^4\theta$ である。$x=\cos\theta$ を使って $\int \sin^4\theta\,\dd\Omega_2=2\pi\int_{-1}^1(1-x^2)^2\dd x=32\pi/15$ と計算すると、真空の瞬間の損失は $-\dot M_B=6(b')^2/(5G)$ になる。$b$ が前後で異なる定数ならメモリーはその差に比例する。一方、$b$ が元の値に戻っても途中で $b'\ne0$ なら総エネルギーは減る。これはこの特定の news の形での結論であり、定常な端点を持つ真空解として実現するには角度依存の質量変化も収支式に合わせる必要がある。
:::

### Lorentz 変換が通常の並進を保つ理由

未来向きの光的ベクトル $l^a=(1,\vb n)$ を取り、正時的な Lorentz 変換 $\Lambda$ により $\Lambda l^a=D(\vb n)(1,\vb n')$ と書く。例えば $z$ 方向の boost では $D=\gamma(1-v n_z)>0$、$n'_z=(n_z-v)/(1-v n_z)$ である。球面上の接ベクトルを $w^a=(0,\delta\vb n)$ とすると、変換後の $\delta l'$ は $D^{-1}\Lambda w$ と光的な $l'$ に比例する補正の和である。$\Lambda w$ は $l'$ に直交するため、接ベクトルのノルムは $D^{-1}$ 倍になる。従って角写像の共形因子は $K=D^{-1}$ である。

また $x^a=r l^a+u(1,0,0,0)$ と書き、変換後の $t'-|\vb x'|$ を $r$ の大きいところで展開する。$rD(1,\vb n')$ の寄与は消え、残りは $-u\,\Lambda(1,0,0,0)\cdot(1,\vb n')=u/D$ である。従って同じ $K$ が角計量と遅延時刻の両方に現れる。並進 $a^a$ の超並進関数は $f_a(\vb n)=-a\cdot l$ である。Lorentz 変換後には $Kf_a(\vb n)=-(\Lambda a)\cdot(1,\vb n')$ となるから、$\ell=0,1$ の空間は保持される。例えばもとの時間並進 $a=(A,0,0,0)$ は、上の boost で $\Lambda a=(\gamma A,0,0,-\gamma v A)$ になる。新枠の関数は $\gamma A(1+v n'_z)$ であり、時間と空間の通常の並進だけの組合せである。

## 演習問題

::::exercise{level=basic}
$u=0$ として未来向きに $r\to\infty$ とすると Minkowski の共形座標でどこへ行くか。

:::solution
$U=\arctan0=0$、$v=2r\to\infty$ から $V\to\pi/2$ なので、$\mathscr I^+$ の $U=0$ の断面である。
:::
::::

::::exercise{level=basic}
Schwarzschild で $m=2G$ の静止切片の ADM 質量を答えよ。

:::solution
:eqref[eq-gr-adm-energy] の面積積分は $4m \cdot4\pi$ なので $m/G=2$ である。
:::
::::

::::exercise{level=basic}
$N_{AB}=0$ かつ $T_{uu}^{(2)}=0$ の間、Bondi 質量はどう変化するか。

:::solution
:eqref[eq-gr-bondi-loss] の右辺が両方とも零なので、その区間で一定である。
:::
::::

::::exercise{level=basic}
$\alpha(\theta,\phi)=3-2\cos\theta$ の超並進の調和次数を求めよ。

:::solution
定数は $\ell=0$、$\cos\theta$ は $\ell=1$ の軸対称成分である。よって通常の時間並進と軸方向の空間並進の組である。
:::
::::

::::exercise{level=advanced}
等方的な $h_{ij}=(1+B/r)\delta_{ij}+\order{r^{-2}}$ の $E_{\rm ADM}$ を求めよ。

:::solution
$\partial_jh_{ij}-\partial_i h_{jj}=-2B\partial_i(1/r)=2Bs_i/r^2$ であり、面積積分から $E_{\rm ADM}=B/(2G)$ となる。
:::
::::

::::exercise{level=advanced}
$\xi$ を $3\xi$ に置き換えたときの Komar 積分の変化を求めよ。

:::solution
添字を下げる操作、外微分、Hodge 双対、積分はいずれも線形なので三倍になる。通常の質量では無限遠での単位ノルム規格化を行う。
:::
::::

::::exercise{level=advanced}
平坦な三次元切片の線形化された曲率で、$h_{ij}=\delta_{ij}+2\Phi\delta_{ij}$ の摂動を計算せよ。

:::solution
二重発散は $2\Delta\Phi$、跡の Laplacian は $6\Delta\Phi$ である。差は $-4\Delta\Phi$ となる。
:::
::::

::::exercise{level=advanced}
$\int N_{AB}N^{AB}\dd\Omega_2=32\pi G C$ が一定で物質流が零の間、$u_0$ からの質量を求めよ。

:::solution
:eqref[eq-gr-bondi-loss] により $\dot M_B=-C$ であり、積分すると $M_B(u)=M_B(u_0)-C(u-u_0)$ である。
:::
::::

::::exercise{level=advanced}
空間無限遠の計量に純粋な座標摂動 $h_{ij}=\delta_{ij}+2\partial_i\partial_j\varphi$ を線形で入れる。適切な減衰のある $\varphi$ に対して ADM エネルギーの一次寄与を求めよ。

:::solution
一次の積分核は $\partial_jh_{ij}-\partial_i h_{jj}$ である。摂動の最初の項は $2\partial_j\partial_i\partial_j\varphi=2\partial_i\Delta\varphi$、次も $2\partial_i\Delta\varphi$ なので差は各点で零である。従ってこの形の摂動の一次 ADM 寄与は零である。例えば $\varphi=c\log r$ を遠方だけで用いると Hessian は $O(r^{-2})$ であり、同じ相殺が成立する。この証明は純粋な Hessian の一次項に関するもので、有限な非線形座標変換の一般の不変性を示すには、落ち方と境界項の制御も必要になる。
:::
::::

## 確認テスト

::::quiz
ADM エネルギーの定義で積分面を向かわせる先を選べ。

- [x] 空間無限遠
- [ ] 遅延時刻一定の未来零無限遠
- [ ] 中心の時間的世界線

:::explanation
一つの空間切片の漸近球面で $R\to\infty$ を取る。一定遅延時刻の零無限遠では Bondi 質量を定義する。
:::
::::

::::quiz{answer="2" tol="1e-9"}
$G=1$、$h_{ij}=(1+4/r)\delta_{ij}$ の ADM エネルギーを求めよ。

:::explanation
積分の法線成分は $8/r^2$ なので $E_{\rm ADM}=8\cdot4\pi/(16\pi)=2$ である。
:::
::::

::::quiz
真空の Bondi 質量が減少し得る理由を選べ。

- [x] 重力波の news が非零である。
- [ ] 物質密度の非零が必須である。
- [ ] 空間無限遠への座標変換で質量を変更できる。

:::explanation
物質流が零でも :eqref[eq-gr-bondi-loss] の news 二乗項は非正である。ある断面で news が非零なら、その断面で重力波による質量の変化率は負になる。
:::
::::

::::quiz{answer="4" tol="1e-9"}
BMS 超並進のうち通常の時空並進に対応する独立成分数を答えよ。

:::explanation
球面上の定数が一成分、$\ell=1$ の三つの方向余弦が三成分で、合計 $4$ である。
:::
::::

::::quiz{answer="6" tol="1e-7"}
$G=1$ で、共形データの $\psi=1+3/r+\order{r^{-2}}$ が与えられた。ADM エネルギーを求めよ。

:::explanation
$h_{ij}=\psi^4\delta_{ij}$ の $1/r$ の係数は $12$ である。等方な係数 $B$ は $E_{\rm ADM}=B/(2G)$ を与えるため、$E_{\rm ADM}=6$ である。
:::
::::

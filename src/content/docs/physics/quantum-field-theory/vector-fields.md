---
title: "ベクトル場と電磁場の量子化"
description: "拘束を含むベクトル場を量子化し、物理偏極と共変伝播関数を求める。"
notation:
  - symbol: 'F_{\mu\nu}'
    meaning: 電磁場の強さ
  - symbol: '\xi'
    meaning: 共変ゲージのパラメータ
  - symbol: '\epsilon^\mu_\lambda'
    meaning: 偏極ベクトル
---

## この章の目標

- Coulomb ゲージの二つの正準自由度を取り出す。
- Gupta–Bleuler 条件と零ノルム商を説明する。
- 質量の有無に応じて偏極和を構成する。

## Maxwell 場の正準量子化と Coulomb ゲージ

### Gauss 則を先に解く

$\mathcal L=-F_{\mu\nu}F^{\mu\nu}/4-j_\mu A^\mu$ とし、$j$ は保存する外部源とする。$\vb E=-\dot{\vb A}-\nabla A^0$、$\vb B=\nabla\times\vb A$ によって運動項は $(|\vb E|^2-|\vb B|^2)/2$ になる。ここでは電荷を有理化した自然単位系を使う。:ref[def-em-gauge] のゲージ自由度を使い、無限遠で消える場では $\nabla^2\chi=\nabla\cdot\vb A$ を解いて $\vb A\mapsto\vb A-\nabla\chi$ とすれば $\nabla\cdot\vb A=0$ にできる。これを :term[Coulomb ゲージ]{en="Coulomb gauge" yomi="くーろんげーじ"} という。$A^0$ の時間微分が無いので共役運動量は零であり、変分すると $-\nabla^2 A^0=\rho$ になる。無限遠での境界条件では解は $A^0(\vb x)=\int\dd^3y\,\rho(\vb y)/(4\pi|\vb x-\vb y|)$ である。

ゲージを固定する前の $\vb A$ の共役運動量は $-\vb E$ である。Gauss 則と Coulomb 条件を課した後、独立変数 $\vb A_T$ に共役なのはその横成分 $\boldsymbol\pi_T=\dot{\vb A}_T$ である。縦と横の交差積分は発散が零のため消え、源つき Hamiltonian は
$$
H=\int\dd^3x\,\left(\tfrac12\boldsymbol\pi_T^2+\tfrac12\vb B^2-\vb j\cdot\vb A_T\right)
 +\tfrac12\int\dd^3x\,\rho A^0 .
$$
$A^0$ の式を使う前は静的な部分に $-\frac12|\nabla A^0|^2+\rho A^0$ が現れ、部分積分すると右辺の Coulomb 項になる。二つの横偏極を正準量子化すると
$$
[A^i_T(\vb x),\pi^j_T(\vb y)]
=i\int\frac{\dd^3p}{(2\pi)^3}\left(\delta_{ij}-\frac{p_i p_j}{|\vb p|^2}\right)e^{i\vb p\cdot(\vb x-\vb y)}
\label{eq-qft-transverse-bracket}
$$
である。行列は横方向への直交射影であり、:eqref[eq-qft-scalar-modes] を偏極ごとに適用すれば得られる。$A^0$ は時刻ごとに $\rho(t,\vb x)$ で決まるため、その応答には瞬時項がある。一方、局所場 $F_{\mu\nu}$ は Maxwell 方程式により $\Box F_{\mu\nu}=\partial_\mu j_\nu-\partial_\nu j_\mu$ を満たす。従って入射する自由波がなければ、縦成分と横成分を足した解は遅延解になる。

### 拘束と横射影

$A^0$ を含めて独立な座標とすると $\pi^0=0$ が一次拘束であり、その時間保存で Gauss 則が出る。ゲージ条件 $\partial_i A^i=0$ と縦運動量の拘束を課すことは、場と運動量を横部分に制限することに対応する。各 $\vb p\ne0$ で二つの単位横ベクトル $e_\lambda$ を選ぶと $A^i=e_\lambda^i q_\lambda$、$\pi^j_T=e_\lambda^j p_\lambda$。正準一形式 $\sum_i\pi^i_T\delta A^i$ は $\sum_\lambda p_\lambda\delta q_\lambda$ である。従って横変数を正準量子化し、その交換関係を元の三成分へ戻すと :eqref[eq-qft-transverse-bracket] になる。例えば $\vb p$ が $(1,1,0)$ に平行なら
$$
P_T=\begin{pmatrix}1/2&-1/2&0\\-1/2&1/2&0\\0&0&1\end{pmatrix}
$$
であるので、元の添字で異なる成分の正準括弧も零とは限らない。

### Dirac 括弧で横の自由度を取り出す

源のない空間 Fourier モード $\vb k\ne0$ を実基底で選び、縦座標と運動量を $q_L,p_L$ とする。Gauss 則と Coulomb 条件によって拘束は $\chi_1=q_L=0,\chi_2=p_L=0$ と書ける。外部源があれば $p_L$ は $\rho$ が決める既知の値に固定される。この値を引いて拘束を定義しても、正準括弧は変わらない。$C_{ab}=\{\chi_a,\chi_b\}$ は $\begin{pmatrix}0&1\\-1&0\end{pmatrix}$。:term[Dirac 括弧]{en="Dirac bracket" yomi="でぃらっくかっこ"} を $\{F,G\}_D=\{F,G\}-\{F,\chi_a\}(C^{-1})_{ab}\{\chi_b,G\}$ とする。$\{\chi_a,G\}_D=0$ なので拘束を厳密な等式として代入できる。$q_i=q_{Ti}+\hat k_i q_L$、$p_j=p_{Tj}+\hat k_j p_L$ に代入すれば $\{q_i,p_j\}_D=\delta_{ij}-\hat k_i\hat k_j$。例えば $\hat k=\vb e_3$ では $z$ 成分の括弧は零、$x,y$ の二組は正準形のままである。

### 横射影が局所デルタでない理由

任意の滑らかな急減少ベクトル $\vb f$ を $\vb f=\vb f_T+\nabla u$ に分けるには、$\nabla^2u=\nabla\cdot\vb f$ を無限遠で消える解で解く。$\vb f_T=\vb f-\nabla u$ は発散が零。Fourier 空間で $u$ の逆 Laplacian は $-1/|\vb k|^2$ を掛けるので $\vb f_T=P_T\vb f$。従って等時刻の横デルタの核は
$$
\delta_{ij}^{T}(\vb r)
=\delta_{ij}\delta^3(\vb r)+\partial_i\partial_j\frac1{4\pi r}.
$$
微分が分布に作用する式であり、$\nabla^2(1/(4\pi r))=-\delta^3(\vb r)$ から $\partial_i\delta_{ij}^{T}=0$ が従う。$r>0$ では第二項は $(3r_ir_j-r^2\delta_{ij})/(4\pi r^5)$。従って Coulomb ゲージのポテンシャルの座標と正準運動量の括弧は空間的に離れた点でも零とは限らない。横成分を選ぶ前に、全空間の分布を使って Gauss 則を解いたためである。

:::example[場の強さの括弧で縦の核が消える]
源がないとき $E_i=-\pi_{T,i}$、$B_j=\epsilon_{jkl}\partial_k A_{T,l}$ とする。等時刻で $[E_i(\vb x),B_j(\vb y)]=i\epsilon_{jkl}\partial_{y^k}\delta^T_{li}(\vb y-\vb x)$。横デルタの非局所項には $\epsilon_{jkl}\partial_k\partial_l$ が作用するため零となる。従って括弧は $i\epsilon_{jki}\partial_{y^k}\delta^3(\vb y-\vb x)$ に等しい。例えば離れた二点に台をもつ試験関数で $E,B$ を塗りならせば、この微分デルタに台が重ならず括弧は零。観測可能な電磁場の局所性と、ゲージポテンシャルの横射影の非局所性は両立する。
:::

### 一次拘束から Gauss 則が出る過程

空間ベクトルの成分を $A_i^{\rm sp}=A^i$、その共役を $\pi_i=-E^i$ と書き、$A^0$ の共役を $\pi^0$ とする。上付きの $\mathrm{sp}$ は Euclid の三成分を表し、計量で下げた $A_i=-A_i^{\rm sp}$ と区別する。Legendre 変換で $\dot{\vb A}=\boldsymbol\pi-\nabla A^0$ を用い、境界項が消えるとする。このとき
$$
H_{\rm tot}=\int\dd^3x\left[
 \tfrac12\boldsymbol\pi^2+\tfrac12\vb B^2-\vb j\cdot\vb A
 +A^0(\nabla\cdot\boldsymbol\pi+\rho)+u\pi^0\right].
$$
$\dot A^0$ の係数は Lagrangian で零なので $\pi^0=0$ が一次拘束であり、$u$ は乗数である。:eqref[eq-qft-field-poisson] を使うと $\dot\pi^0=-\nabla\cdot\boldsymbol\pi-\rho$。従って $G=\nabla\cdot\boldsymbol\pi+\rho=0$ が二次拘束となる。さらに $\dot{\boldsymbol\pi}=-\nabla\times\vb B+\vb j$ なので $\dot G=\nabla\cdot\vb j+\dot\rho=0$。源の保存がここで必要である。$\pi^0$ と $G$ は互いに正準括弧が零で、ゲージ条件を課す前は第一級拘束である。

ゲージを $A^\mu\mapsto A^\mu+\partial^\mu\chi$ とすると、$\vb A\mapsto\vb A-\nabla\chi$、$A^0\mapsto A^0+\dot\chi$。これを時刻ごとの正準変換で実現する生成子は
$$
\mathcal G[\chi]=\int\dd^3x\,(\chi G+\dot\chi\,\pi^0).
$$
確かめるには $\{\vb A,\int\chi\,\nabla\cdot\boldsymbol\pi\}=-\nabla\chi$ と、$\{A^0,\int\dot\chi\,\pi^0\}=\dot\chi$ を用いればよい。$\pi_i$ は不変なので電場も保たれる。Coulomb 条件の時間微分は $\partial_t(\nabla\cdot\vb A)=\nabla\cdot\boldsymbol\pi-\nabla^2 A^0$ であり、$G=0$ と合わせると $-\nabla^2 A^0=\rho$ が得られる。境界で消える残留変換は $\nabla^2\chi=0$ を満たす調和関数なので、全空間で滑らかなら零である。箱の周期条件では空間的に一定なモードが残り、Poisson 方程式の可解条件は全電荷が零であることになる。以下では伝播する $\vb p\ne0$ のモードを量子化する。周期箱で空間的に一定なポテンシャルを含める場合は、その境界条件とゲージ群に応じて別に扱う。

:::example[源があるときの縦運動量]
時刻を固定し、周期境界に $\rho(\vb x)=\rho_0\cos(kz)$、$k\ne0$ を取る。Coulomb の解は $A^0=\rho_0\cos(kz)/k^2$。従って $\pi_L^z=\partial_z A^0=-\rho_0\sin(kz)/k$ であり、$\partial_z\pi_L^z=-\rho$ が成り立つ。縦運動量の値は零でないが、新たな正準量子を作る振動の自由度でもない。源が $\rho_0(t)$ に依存する場合は $j^z=-\dot\rho_0(t)\sin(kz)/k$ を伴えば、同じ時刻ごとの Gauss 則と連続の式を満たす。
:::

### 横場のモード規格化とエネルギー

源を零にし、$\vb p\ne0$ に対し横偏極 $e_\lambda(\vb p)$ を二つ選ぶ。$e_\lambda^*\cdot e_\sigma=\delta_{\lambda\sigma}$、$\vb p\cdot e_\lambda=0$。自由方程式は $(\partial_t^2-\nabla^2)\vb A_T=0$ なので
$$
\vb A_T(x)=\sum_{\lambda=1}^2\int_{\vb p\ne0}
 \frac{\dd^3p}{(2\pi)^3 2|\vb p|}
 \left[e_\lambda(\vb p)a_\lambda(\vb p)e^{-ipx}
       +e_\lambda^*(\vb p)a_\lambda^\dagger(\vb p)e^{ipx}\right],
\qquad p^0=|\vb p|.
$$
各 $a_\lambda$ の交換子はスカラーと同じ $(2\pi)^3 2|\vb p|\delta^3$ を持ち、異なる偏極では零とする。$\boldsymbol\pi_T=\dot{\vb A}_T$ と合わせて等時刻交換子を取ると、消滅・生成の二項の一つずつに $\sum_\lambda e^i_\lambda e^{*j}_\lambda=P_T^{ij}(\vb p)$ が入る。これは実で $\vb p$ の反転に対して偶だから、二項を同じ Fourier 指数へ直せる。従って :eqref[eq-qft-transverse-bracket] の規格化を再現する。

エネルギーには $\pi_T^2$ と $|\nabla\times\vb A_T|^2$ がある。各平面波で $\nabla\times$ は $i\vb p\times$ を掛け、横偏極で $(\vb p\times e_\lambda)^*\cdot(\vb p\times e_\sigma)=|\vb p|^2\delta_{\lambda\sigma}$。箱で実の横基底ごとにモードを取れば、各独立座標は明示的に $\frac12(\pi_\lambda^2+|\vb p|^2 q_\lambda^2)$ となる。進行波の規格化へ移せば
$$
:P^\mu:
=\sum_{\lambda=1}^2\int\frac{\dd^3p}{(2\pi)^3 2|\vb p|}
       p^\mu a_\lambda^\dagger(\vb p)a_\lambda(\vb p).
$$
空間成分は Poynting 流 $\int\vb E\times\vb B$ からも得られる。実際、空間成分で $(\vb E\times\vb B)^i=E_j^{\rm sp}\partial_i A_j^{\rm sp}-E_j^{\rm sp}\partial_j A_i^{\rm sp}$、$\vb E=-\boldsymbol\pi_T$。第二項は $\nabla\cdot\boldsymbol\pi_T=0$ と境界条件で積分が消え、残りは $\int\pi_j\partial^i A_j^{\rm sp}\dd^3x$。スカラーの並進電荷と同様に一量子へ $\vb p$ を与える。例えばエネルギー $4\,\mathrm{GeV}$ の同じ向きに二つの直交偏極を一個ずつ置くと、全エネルギーは $8\,\mathrm{GeV}$ であり、時間成分の振動子を別に足す必要はない。

### 外部電流が変位を作る

外部源 $j^\mu$ を滑らかで時空にコンパクトな台を持つ関数とする。Coulomb の Hamiltonian のうち場に線形なのは $-\int \vb j\cdot\vb A_T$。従ってスカラー源の場合と同様、時間順序を外すときの交換子は数である。$\rho A^0/2$ の項も源だけで決まる数なので、場の状態には位相を与える。源を止めた後の真空からの状態は、各横モードの変位
$$
\alpha_\lambda(\vb p)=i\int\dd^4x\, e^{ipx}
             e_\lambda(\vb p)^*\cdot\vb j(x),
\qquad
\langle N\rangle=\sum_\lambda\int\dd\Pi_p\,|\alpha_\lambda(\vb p)|^2
$$
で決まる。これは古典電流で励起した自由光子に対する式であり、電流を量子演算子にすると交換子は一般に数ではなくなる。$p^0=|\vb p|$ を満たす成分が現れるので、$\vb j$ の非零成分が全て放射に寄与するわけではない。

例えば $p$ が $z$ 軸にあり、$\widetilde{\vb j}(p)=(j_x,ij_x,j_z)$ なら横線偏極で数を足した重みは $2|j_x|^2$。円偏極 $e_\pm=(\vb e_x\pm i\vb e_y)/\sqrt2$ で見ると $e_+^*\cdot \widetilde{\vb j}=\sqrt2 j_x$、$e_-^*\cdot\widetilde{\vb j}=0$ で、一方の偏極だけが励起される。$j_z$ は電流保存による $j^0$ との関係を満たしていても、この放射の数には現れない。これは源の非放射部分が Coulomb 場を作ることと整合する。

### 横射影を用いた遅延解から局所な場を復元する

空間 Fourier で $\nabla=i\vb k$、Coulomb の条件を $\vb k\cdot\vb A_T=0$ とする。$\vb A_T$ は
$$
(\partial_t^2+|\vb k|^2)\vb A_T=P_T\vb j,\qquad
A^0=\rho/|\vb k|^2
$$
を満たす。$\vb E=-\dot{\vb A}_T-i\vb k\rho/|\vb k|^2$ に同じ波動演算子を作用させる。第一項は $-P_T\dot{\vb j}$。第二項では連続の式 $\dot\rho=-i\vb k\cdot\vb j$ をもう一回微分し、
$$
-\frac{i\vb k}{|\vb k|^2}(\ddot\rho+|\vb k|^2\rho)
=-\frac{\vb k\vb k^\mathsf T}{|\vb k|^2}\dot{\vb j}
 -i\vb k\rho
$$
となる。二つを足すと $-\dot{\vb j}-i\vb k\rho$、すなわち $\Box\vb E=-\dot{\vb j}-\nabla\rho$。$\vb B=i\vb k\times\vb A_T$ に作用させると横射影の縦成分は外積で消え、$\Box\vb B=\nabla\times\vb j$。時刻以前に源と場が零なら、これらの局所的な源に対する遅延解である。瞬時の $A^0$ と非局所な射影で選んだ $\vb A_T$ の両方を入れると、観測される電磁場の応答に局所方程式が戻る。

## 共変量子化と Gupta–Bleuler 法

### 共変性のために四成分を残す

$-(\partial\cdot A)^2/2$ を Lagrangian に加え、部分積分すると自由密度は $-\partial_\mu A_\nu\partial^\mu A^\nu/2$ である。各成分は波動方程式に従うが時間成分の運動項の符号が負なので、モード交換子は $[a^\mu(\vb p),a^{\nu\dagger}(\vb q)]=-\eta^{\mu\nu}(2\pi)^3 2E\delta^3(\vb p-\vb q)$ となる。時間方向の一量子のノルムは負である。

:term[Gupta–Bleuler 条件]{en="Gupta–Bleuler condition" yomi="ぐぷたぶろいらーじょうけん"} は $(\partial_\mu A^\mu)^{(+)}\ket{\mathrm{phys}}=0$ である。ここで正周波数部分は消滅演算子だけを含む。条件により物理状態間で $\partial\cdot A$ の行列要素は消える。$\vb p$ が $z$ 軸のとき条件は $(a^0-a^3)\ket{\mathrm{phys}}=0$。この消滅演算子は $(a^{0\dagger}-a^{3\dagger})$ と可換であり、後者が作る状態は零ノルムである。有限励起の状態を生成演算子の多項式で書くと、条件は時間・縦方向を合わせた一変数に関する微分が零という条件になる。残る変数は零ノルム方向だけであり、それを割った空間には横の二偏極が残る。

### 零ノルム方向の商を取る

固定した波数に対し生成変数 $x=a^{0\dagger}$、$y=a^{3\dagger}$ を使うと $a^0-a^3$ は多項式上で $-\partial_x-\partial_y$ である。条件を解くには $u=x-y,v=x+y$ と変数変換すればよく、$-2\partial_v$ が零を与える多項式は $u$ だけに依存する。さらに $[u^\dagger,u]=0$、$u^\dagger\ket0=0$ なので $u^n$ を含む物理状態と任意の物理状態の内積は $n\ge1$ で零である。従って $u$ が生成する成分を全て同一視した商に横二振動子の正定値内積が入る。例えば $(a^{1\dagger})^2\ket0+u^2\ket0$ の物理ノルム二乗は $2$ であり、$u^2\ket0$ の係数を変えても同じである。

### 共変ゲージでの正準交換関係を求める

Feynman ゲージで表面項を除いた自由作用は $-\partial_\mu A_\nu\partial^\mu A^\nu/2$。座標を $A^\nu$ に取ると共役は $\Pi_\nu=-\dot A_\nu$ である。従って $[A^\mu,\Pi_\nu]=i\delta^\mu_\nu\delta^3$ から
$$
[A^\mu(t,\vb x),\dot A^\nu(t,\vb y)]
=-i\eta^{\mu\nu}\delta^3(\vb x-\vb y).
$$
各成分を :eqref[eq-qft-scalar-expansion] の質量零版で展開し、等時刻の逆モード積分を使うと、右辺の $-\eta^{\mu\nu}$ がそのままモードの交換関係に移る。時間成分について $[a^0,a^{0\dagger}]=-1$ とした一モードでは、真空で $\langle0|a^0 a^{0\dagger}|0\rangle=-1$。従ってここでの $\dagger$ は不定内積についての随伴を表す。四成分をそのまま正定値 Hilbert 空間に置いてこの交換関係と消滅真空を両立させることはできない。

元の源なし Maxwell 方程式は $\partial_\mu F^{\mu\nu}=\Box A^\nu-\partial^\nu(\partial\cdot A)=0$ である。ゲージ固定作用では $\Box A^\nu=0$ だが、全空間で $\partial\cdot A=0$ を演算子恒等式とすると、正準交換子に矛盾する。実際、$t$ を一致させた $[\partial\cdot A(t,\vb x),A^0(t,\vb y)]=i\delta^3(\vb x-\vb y)$ は時間微分から来る。そこで正周波数だけを物理状態に課し、負周波数に対しては物理ブラに随伴条件を課す。これで物理状態の間の $\partial\cdot A$ とその微分の行列要素が消え、Maxwell 方程式が物理的な行列要素で成り立つ。

### 商の上で局所演算子が作用すること

任意の光的 $p$ に対し消滅側の拘束を $p_\mu a^\mu(p)$ と書く。場の強さの生成側は $p_\rho a_\sigma^\dagger-p_\sigma a_\rho^\dagger$ である。前者と後者の交換子は $-p_\rho p_\sigma+p_\sigma p_\rho=0$ に比例する。このため $F_{\rho\sigma}$ は Gupta–Bleuler 条件の核を保つ。随伴も同じ核を保つので、任意の物理状態 $v$ と零ノルムの方向 $n$ に対し $\langle v,F_{\rho\sigma}n\rangle=\langle F_{\rho\sigma}v,n\rangle=0$。従って $F_{\rho\sigma}n$ は再び物理内積の根基に入り、$F_{\rho\sigma}$ は零ノルムを同一視した商へ降りる。この議論は試験関数で塗りならした演算子を有限粒子の共通領域で用いたものである。

同じ結論を偏極で見る。$p$ に比例する零方向の一粒子偏極は $p\wedge p=0$ により $F$ を生成しない。例えば $p=(\omega,0,0,\omega)$ に対し $A$ の時間・縦成分が等しい平面波では $F^{03}=-i\omega(A^3-A^0)=0$。一方、横の $A^1$ が一なら $F^{01}=-i\omega$、$F^{31}=-i\omega$ が残る。商の正の二偏極が場の強さの実際の行列要素を作る。

:::example[規格化した物理状態を商で比べる]
一箱モードで $u=a^{0\dagger}-a^{3\dagger}$、$e^\dagger=a^{1\dagger}$ とする。$|\Psi\rangle=(e^\dagger+c u)|0\rangle$ は任意の複素係数 $c$ で物理条件を満たす。$\langle\Psi|\Psi\rangle=1$ かつ、二つの異なる $c$ の差は $u|0\rangle$ に比例する零方向である。例えば $F_{01}$ を空間・時間の同じモードに射影すれば、真空との行列要素に寄与するのは $e^\dagger|0\rangle$ だけである。不定内積で自己内積が零の任意のベクトルを割るのではなく、ここでは物理条件の核の中で全物理状態と直交する方向だけを割っている。
:::

### 保存源の発展が物理条件を保つ

場に線形な源 $-j_\mu A^\mu$ を Lagrangian に加える。Feynman ゲージでは $\Box A^\mu=j^\mu$ なので、発散を取ると $\Box(\partial\cdot A)=\partial\cdot j=0$。従って保存源は物理条件に入る自由な発散場を直接励起しない。モード変位で見るには、源が入る前後を自由とし $\int\dd^4x\,j_\mu A^\mu$ を用いた発展を作る。この指数と $p_\mu a^\mu(p)$ の交換子は $p_\mu\widetilde j^\mu(p)$ に比例する。コンパクトな台の保存源では Fourier 変換した連続の式によりこれは零。よって相互作用による変位は物理部分空間を保ち、随伴も保つから零方向の同値類も保つ。例えば $p$ に平行な $\widetilde j^\mu$ は光的な殻で $p\cdot\widetilde j=0$ を満たすが、物理商では零方向しか作らない。

### 一般の共変ゲージと発散場

$\xi$ を有限の非零の実数とし、
$$
\mathcal L_\xi=-\tfrac14F_{\mu\nu}F^{\mu\nu}
        -\frac{(\partial\cdot A)^2}{2\xi}-j_\mu A^\mu
$$
を変分する。第二項の変分は部分積分で $+(\partial^\nu\partial\cdot A)\delta A_\nu/\xi$。従って
$$
\Box A^\nu-(1-1/\xi)\partial^\nu(\partial\cdot A)=j^\nu,\qquad
\Box(\partial\cdot A)=\xi\,\partial\cdot j.
$$
保存源なら全ての $\xi\ne0$ で発散場は自由である。ただし $\xi\ne1$ では各 $A^\nu$ は独立な単純波動方程式にならず、縦の部分に二重極が現れる。四個の独立な振動子に $-\eta^{\mu\nu}$ の交換子を付けた先の明示構成は $\xi=1$ のものである。Landau ゲージは逆核の $\xi\to0$ の極限として使うか、補助場の変分で $\partial\cdot A=0$ を課す。後者では $b\,\partial\cdot A+\xi b^2/2$ と置き、$\xi\ne0$ で $b=-\partial\cdot A/\xi$ を消去すれば同じ密度になる。$\xi=0$ では $b$ の変分が直接 $\partial\cdot A=0$ を課すため、$\xi\ne0$ のように $b$ の値を解くことはできない。

## 光子の伝播関数とゲージ依存性

$-(\partial\cdot A)^2/(2\xi)$ を加えたとき、運動量での二次演算子は $K_{\mu\nu}=-k^2\eta_{\mu\nu}+(1-1/\xi)k_\mu k_\nu$。この節では三次元の横射影とは別に、四次元の射影を使う。$k^2\ne0$ で $P_L^{\mu\nu}=k^\mu k^\nu/k^2$、$P_T=\eta-P_L$ と置けば横成分の固有値は $-k^2$、縦成分は $-k^2/\xi$ である。従って Feynman 処方をとった伝播関数は
$$
D_F^{\mu\nu}(k)=\frac{-i}{k^2+i0}
 \left(\eta^{\mu\nu}-(1-\xi)\frac{k^\mu k^\nu}{k^2+i0}\right)
\label{eq-qft-photon-propagator}
$$
となる。$\xi=1$ が Feynman ゲージ、$\xi\to0$ が Landau ゲージである。保存電流との交換振幅では縦項の両端に $k\cdot j$ が現れて消える。また自由な $F_{\mu\nu}$ の二点関数で $A$ の各端点を反対称微分すると $k_\mu k_\nu$ による違いが消える。従って場 $A^\mu$ の相関と、局所可観測な $F_{\mu\nu}$ の相関を区別する。

### 保存源で二つのゲージを照合する

Coulomb ゲージで $D_C^{00}(k)=i/|\vb k|^2$、$D_C^{0i}=0$、$D_C^{ij}=iP_T^{ij}/(k^2+i0)$ である。時間成分は Gauss 則 $-\nabla^2 A^0=\rho$ を積分に戻した瞬時項、空間成分は横振動子の伝播関数から出る。保存する二つの電流で挟むと、Feynman ゲージでは :eqref[eq-qft-photon-propagator] により $-i(j^0j'^0-\vb j\cdot\vb j')/(k^2+i0)$。保存条件 $j_L=k^0j^0/|\vb k|$ を入れると
$$
-\frac{i}{k^2+i0}(j^0j'^0-\vb j\cdot\vb j')
=\frac{i}{|\vb k|^2}j^0j'^0+\frac{i}{k^2+i0}\vb j_T\cdot\vb j'_T .
$$
$1-(k^0)^2/|\vb k|^2=-k^2/|\vb k|^2$ を使った。右辺は Coulomb ゲージの電流の縮約そのものである。例えば $k^0=0$ で空間電流の無い静電配置では、どちらも $i j^0j'^0/|\vb k|^2$ だけを与える。

### 場の強さの二点関数

Fourier 変換 $\partial_\mu\mapsto-ik_\mu$ の下で二点関数の両端に $F$ を置くと、もう一方の点には逆の符号で微分が作用する。従って Feynman ゲージで得た縮約は
$$
\frac{-i}{k^2+i0}\bigl(k_\mu k_\rho\eta_{\nu\sigma}
-k_\mu k_\sigma\eta_{\nu\rho}
-k_\nu k_\rho\eta_{\mu\sigma}
+k_\nu k_\sigma\eta_{\mu\rho}\bigr).
$$
ここでは時間順序積に微分を含めた共変な分布を計算した。演算子の時刻で時間順序を取る定義と比較すると、時間微分が階段関数に当たる接触項を区別する。各端の反対称化により :eqref[eq-qft-photon-propagator] の $k_\alpha k_\beta$ は全て消える。例えば $\vb k$ が $z$ 軸のとき $F_{01}$ の二点関数の分子は $-(k^0)^2$、$F_{12}$ の分子は零である。

### 伝播関数の二重極と分布の境界値

共変ゲージの $D_F^{\mu\nu}$ で、$\xi\ne1$ の項には $(k^2+i0)^{-2}$ がある。これは二つの独立な光子の質量殻の積ではなく、同じ二次逆核の縦成分に由来する。実変数 $z$ について分布を微分すれば
$$
\frac1{(z+i0)^2}=\operatorname{Pf}\frac1{z^2}+i\pi\delta'(z)
$$
であり、$\operatorname{Pf}1/z^2=-\partial_z\operatorname{PV}1/z$ と定めた。$z\delta'(z)=-\delta(z)$、また $z\operatorname{Pf}(1/z^2)=\operatorname{PV}(1/z)$ なので、$z/(z+i0)^2=1/(z+i0)$。従って $k_\mu D_F^{\mu\nu}=-i\xi k^\nu/(k^2+i0)$ が分布として成り立つ。ここでは零運動量の扱いを別にした、通常の Feynman 境界値で縮約している。$\delta(k^2)^2$ を使って二重極を表すと、この微分関係を失う。

:::example[ゲージパラメータに依らない線の両端]
二つの保存電流で $D^{\mu\nu}$ を挟めば、$\xi$ に比例する差は $(j\cdot k)(j'\cdot k)/(k^2+i0)^2$。各縮約が恒等的に零なので、この縦項は境界値を取る前に消える。例えば $k=(5,0,0,3)$、$j=(3,1,0,5)$ とすると $j\cdot k=15-15=0$。ここでは $k^2=16$ で質量殻外だが、同じ電流保存の計算は光的な極の近傍にも使える。
:::

### 静電の項は時間成分の自由振動ではない

Coulomb 条件で $A^0$ の二次作用だけを分ける。空間逆核を $K=-\nabla^2$ とし、与えた境界条件で正の逆 $K^{-1}$ がある部分空間に制限する。この作用は時刻ごとに $\int \frac12 A^0 K A^0-\rho A^0$ である。$A^0=K^{-1}\rho+a$ とずらすと交差項が消え、
$$
S_{A^0}=\frac12\int\dd t\,\dd^3x\,aKa
-\frac12\int\dd t\,\dd^3x\,\rho K^{-1}\rho
$$
となる。後者が Hamiltonian の Coulomb エネルギーと逆符号である。源を二つ変分した交換核は $i K^{-1}\delta(t-t')$ になり、運動量では $i/|\vb k|^2$ と表される。自由の Coulomb ゲージで制約を解いた Hilbert 空間では、$\rho=0$ のとき $A^0=0$ である。従ってここでの $D_C^{00}$ は、その空間で二個の自由な時間成分を時間順序で掛けた期待値という意味ではなく、源どうしの瞬時の結合を表す規則である。

:::example[静止した二電荷の交差項]
異なる場所 $\vb x_1,\vb x_2$ に点電荷 $q_1,q_2$ を固定する。$\frac12\int\rho K^{-1}\rho$ に $\rho=q_1\delta_1+q_2\delta_2$ を入れると、対角には点源の自己エネルギーが生じる。離れた二点の交差項は二つあり、前の $1/2$ が消えて $q_1q_2/(4\pi|\vb x_1-\vb x_2|)$ になる。$q_1q_2>0$ なら、距離を拡げるとこのエネルギーは下がる。静的な源の間の力は横の光子の放射を要せず、拘束を消去した項だけからも得られる。
:::

### 共変ゲージ間の差を時空で読む

:eqref[eq-qft-photon-propagator] の二つの $\xi$ の差は $\delta D^{\mu\nu}=-i(\xi_2-\xi_1)k^\mu k^\nu/(k^2+i0)^2$。時空に戻すと $\partial^\mu\partial^\nu$ が掛かる純縦の核である。試験関数として保存する電流 $j,j'$ を入れ、境界積分が消えるとする。例えば一端の微分を部分積分すれば $\partial_\mu j^\mu=0$ が出て消える。従って一端が厳密に保存すれば差は零になる。不連続に切った電流では、時間の端に $\partial\cdot j$ の分布が出る場合がある。計算の切替えを時刻だけで滑らかにしても、$\rho\ne0$ なら $\partial_t g(t)\rho$ が生じるので、自動的に保存則は保たれない。電流保存を使うゲージ比較では、電流を作る物理過程も含めて発散を取る必要がある。

:::example[自由な二局所点での因果性]
Feynman ゲージで $[A^\mu(x),A^\nu(y)]=-\eta^{\mu\nu}C_0(x-y)$。$C_0$ は :eqref[eq-qft-scalar-commutator] の $m=0$ の値である。$F$ 同士の交換子はこれを両端で反対称微分した有限和になる。例えば $F_{01}$ 同士では $(\partial_{x^0}\partial_{y^0}-\partial_{x^1}\partial_{y^1})C_0$ を得る。分布の微分は台を広げないため、離れた二領域が互いに空間的なら塗りならした $F$ の交換子は零。Gupta–Bleuler 商に移した後も、$F$ が条件を保つので同じ式が作用する。真空の時間順序相関は空間的分離で非零でも、この交換子の局所性とは矛盾しない。
:::

### 時間微分の接触項を横振動子で求める

源を零とする。横空間の自由二点を $D_T^{ij}(x,y)=\langle0|T A_T^i(x)A_T^j(y)|0\rangle$ と置く。第一時刻の微分では等時刻の $[A_T^i,A_T^j]=0$ により境界項が消える。第二時刻でも微分すると、階段関数の境界から $[A_T^j(y),\dot A_T^i(x)]$ が出る。:eqref[eq-qft-transverse-bracket] に従って
$$
\partial_{x^0}\partial_{y^0} D_T^{ij}
=\langle0|T\dot A_T^i(x)\dot A_T^j(y)|0\rangle
+i\delta(x^0-y^0)\delta_{ij}^T(\vb x-\vb y).
$$
従って $k^0$ を両端の微分係数とする共変な計算では、この接触項も含む。例えば一横振動子で $D(t)=e^{-i\omega|t|}/(2\omega)$ とすれば、$\partial_x\partial_y D=-D''=\omega^2D+i\delta(t)$。演算子の時間順序で $\dot A$ 同士を取ると $\omega^2D$ だけであり、二つの方法の違いがこの式で検算できる。共変ゲージ間の縦項が $F$ の反対称化で消えるという主張と、異なる時間微分の定義の間に接触項があるという主張は両立する。

### ゲージ変換を演算子の二次逆核で追う

$k^2\ne0$ で共変ゲージの運動量空間の行列を $K=-k^2 P_T-(k^2/\xi)P_L$ とする。添字の片方を下げれば $P_T,P_L$ は冪等で積は零。逆 $K^{-1}=-(P_T+\xi P_L)/k^2$ に $i$ を掛けて境界値を取ると :eqref[eq-qft-photon-propagator] に戻る。$k$ を二度縮約すると $k_\mu D^{\mu\nu}k_\nu=-i\xi$。位置ではこの式は局所の接触分布である。各 $k^2$ を極へ置いて零とすると、この接触項を落とす。$\xi=1$ では $D^{\mu\nu}=-\eta^{\mu\nu}\widetilde\Delta_F$ で、四つの共変振動子の交換子から得られる。$\xi\to0$ では一端の $k$ 縮約が分布として零になるが、三次元の Coulomb 条件 $\vb k\cdot\vb A=0$ を課した結果とは成分が一致しない。

:::example[静的核と共変縦成分の違い]
$k^0=0$、$\vb k=(0,0,k_z)$、$k_z\ne0$ とする。Landau と Feynman のいずれでも $D^{00}=i/k_z^2$。一方 Feynman ゲージでは $D^{33}=-i/k_z^2$、Landau では縦を射影して零になる。この $D^{33}$ は同じ静的な保存源の交換には現れない。実際 $k_\mu j^\mu=-k_z j^3=0$ なので、$j^3$ を一端に入れることはできない。この例では二つのゲージで成分を比較すると違いがあり、電流保存まで入れると観測に対応する縮約が一致する。
:::

## 質量のあるベクトル場（Proca 場）

質量項 $m^2 A_\mu A^\mu/2$ を加えた :term[Proca 場]{en="Proca field" yomi="ぷろかば"} は $m>0$ でゲージ不変ではない。変分で $\partial_\mu F^{\mu\nu}+m^2A^\nu=j^\nu$、発散で保存源に対し $\partial\cdot A=0$。残りは $(\Box+m^2)A^\nu=j^\nu$ である。$A^0$ は依然として共役運動量を持たず、保存源の下でも独立に初期値を設定できない。静的な点源では $(-\nabla^2+m^2)A^0=q\delta^3$ であり、球対称で原点を除く解は $c\,e^{-mr}/r$。微小球で勾配の流束を積分すると $4\pi c=q$ である。

自由な二次演算子は $-(k^2-m^2)P_T+m^2P_L$ なので、逆に $i$ を掛けて
$$
D_{\rm Proca}^{\mu\nu}(k)=\frac{-i}{k^2-m^2+i0}
\left(\eta^{\mu\nu}-\frac{k^\mu k^\nu}{m^2}\right)
\label{eq-qft-proca-propagator}
$$
を得る。縦項は $m\to0$ で大きくなるが保存源間の振幅で消える。他方、運動量に直交する物理的縦偏極は有限質量で存在する。

### Proca 場の正の Hamiltonian

源を零にし、$\boldsymbol\pi=\dot{\vb A}+\nabla A^0$ を用いる。部分積分後の Hamiltonian 密度は
$$
\mathcal H=\tfrac12\boldsymbol\pi^2+\tfrac12\vb B^2+\tfrac12m^2\vb A^2
  +A^0\nabla\cdot\boldsymbol\pi-\tfrac12m^2(A^0)^2
$$
となる。$A^0$ の変分から $m^2A^0=\nabla\cdot\boldsymbol\pi$。従って消去した Hamiltonian に残る追加項は $(\nabla\cdot\boldsymbol\pi)^2/(2m^2)$ であり、全体が非負である。縦モードで $\vb B=0$、波数の大きさを $p$ とすると、一モードは $\frac12(1+p^2/m^2)\pi_L^2+\frac12m^2 A_L^2$。正準方程式を二回使うと $\ddot A_L+(p^2+m^2)A_L=0$ である。従って縦方向もエネルギーが正の物理モードであり、質量ゼロの拘束系とは自由度の数が異なる。

### Stueckelberg 場を用いた記述

$A_\mu-\partial_\mu\sigma/m$ を Proca の質量項に入れると、$A_\mu\mapsto A_\mu+\partial_\mu\alpha,\sigma\mapsto\sigma+m\alpha$ で不変になる。$\sigma=0$ は $m>0$ で選べるゲージである。一方、質量項を展開すると混合項 $-m A^\mu\partial_\mu\sigma$ は部分積分で $+m\sigma\partial_\mu A^\mu$ となる。$-(\partial\cdot A+\xi m\sigma)^2/(2\xi)$ を加えればこれを相殺し、$\sigma$ の二次項は $(\partial\sigma)^2/2-\xi m^2\sigma^2/2$ となる。例えば $m=3\,\mathrm{GeV},\xi=4$ ならこの補助記述のスカラーの極は $p^2=36\,\mathrm{GeV}^2$、Proca の物理極は $9\,\mathrm{GeV}^2$ である。

### Proca の拘束行列と正準変数の数

一次拘束を $\chi_1=\pi^0$、二次拘束を $\chi_2=\nabla\cdot\boldsymbol\pi-m^2 A^0$ と取る。等時刻の正準括弧から $\{\chi_1(\vb x),\chi_2(\vb y)\}=m^2\delta^3(\vb x-\vb y)$、逆順はその負号である。対応する拘束行列が可逆なものを第二級拘束と呼ぶ。従って位相空間の八成分から二区分の条件を除き六次元、すなわち三組の正準変数が残る。$A^i,\pi^j$ の Dirac 括弧を計算すると、第一拘束との括弧は両方とも零である。拘束行列の逆の対角が零なので補正全体が消え、三組の空間成分には標準の正準括弧が残る。

:::example[質量ゼロで同じ手順が使えない理由]
$m=0$ では二つの拘束が互いに可換で、その行列は反転できない。この場合は一次・二次拘束に加えて、ゲージ条件を選んでゲージ軌道の自由度を取り除く。例えば $\vb p\ne0$ で時間方向の拘束の組と縦方向の拘束の組を一つずつ取り除くと、四つの位相空間成分が残る。これらは横二偏極の二組であり、有限質量の三組の単なる数値極限ではない。
:::

### 質量が小さいときの縦成分の放射

縦偏極 $\epsilon_L=(p,0,0,E)/m$ に対して自由平面波の場の強さは $-i(k^\mu\epsilon_L^\nu-k^\nu\epsilon_L^\mu)$ に比例する。$0,3$ 成分では $E(E/m)-p(p/m)=m$ なので、$A^0,A^3$ がそれぞれ $1/m$ で増えても $F^{03}$ は $m$ に比例して小さくなる。横偏極 $\epsilon_x=(0,1,0,0)$ なら $F^{01}$ の係数は $E$、$F^{31}$ の係数は $p$。従って $E$ を固定した $m\to0$ で、縦の場の強さは横の電場成分に比べて $m/E$ となる。

:::example[偏極の規格化と場の強さの比較]
$m=3\,\mathrm{GeV},E=5\,\mathrm{GeV},p=4\,\mathrm{GeV}$ では縦の偏極の時間成分は $4/3$、方向成分は $5/3$。反対称化した振幅は $25/3-16/3=3\,\mathrm{GeV}$ となり、横の対応する電場振幅 $5\,\mathrm{GeV}$ との比は $3/5$ である。同じ $p$ で $m$ だけを小さくすると $\epsilon_L$ の各成分の大きさは増えるが、$F$ のこの比は逆に小さくなる。縦モードの数え上げと保存源への結合の強さを別々に扱う理由がこれで分かる。
:::

### Proca のモードを正準括弧と一致させる

自由場で $E=\sqrt{p^2+m^2}$ を使い、$\vb p$ を $z$ 方向とする。横二成分は $\pi_T^2/2+E^2 A_T^2/2$ の振動子、縦成分は $(E/m)^2\pi_L^2/2+m^2A_L^2/2$。$Q_L=mA_L/E$、$P_L=E\pi_L/m$ とすると縦成分も $P_L^2/2+E^2 Q_L^2/2$ になる。従って $A_L$ の平面波の係数には横成分と比べて $E/m$ が必要である。拘束 $m^2A^0=\nabla\cdot\boldsymbol\pi$ と正周波数の $\pi_L=-i(m^2/E^2)E A_L$ を用いると、$A^0=(p/E)A_L$。結果として縦偏極は $(p,0,0,E)/m$ となる。

任意の運動量へ回して、三偏極を $\epsilon_\lambda^\mu(p)$ とする。自由演算子
$$
A^\mu(x)=\sum_{\lambda=1}^3\int\dd\Pi_p\,
  [\epsilon_\lambda^\mu(p)a_\lambda(\vb p)e^{-ipx}
   +\epsilon_\lambda^\mu(p)^*a_\lambda^\dagger(\vb p)e^{ipx}]
$$
に対し、三組の $a_\lambda$ に正ノルムを与える。すると等時刻の $\vb A,\boldsymbol\pi$ で、縦の $E/m$ とその共役に出る $m/E$ は相殺し、横二成分に縦一成分を加えた単位行列を得る。これは Proca の拘束を消去した Dirac 括弧と一致する。各正規化した振動子のエネルギーは $E$ であり、正規順序化した $H$ は三偏極にわたる $E\,a^\dagger a$ の和になる。

同じモード式から $[A^\mu(x),A^\nu(y)]$ を取れば
$$
[A^\mu(x),A^\nu(y)]
=\left(-\eta^{\mu\nu}-m^{-2}\partial^\mu\partial^\nu\right)C_m(x-y)
$$
になる。微分は $x$ に作用する。$C_m$ は自由実スカラーの交換子、$\partial^\mu\partial^\nu$ は質量殻の二つの指数へともに $-p^\mu p^\nu$ を掛ける。従って偏極和を再現し、空間的分離で交換子は零。Coulomb ゲージの横だけを使った交換子とは異なり、Proca の三振動子で構成した四元ポテンシャル自体に局所な交換子がある。

:::example[縦の共役運動量が発散を持つモード]
$m=3,p=4,E=5\,\mathrm{GeV}$ とする。縦の $A_L$ の一波の係数を $5/3$、波の時間依存を $e^{-iEt}$ とすると、$\pi_L$ の係数は $-i m^2 A_L/E=-3i\,\mathrm{GeV}$。空間の微分を掛けると $\nabla\cdot\boldsymbol\pi$ の係数は $12\,\mathrm{GeV}^2$。$m^2$ で割った $A^0$ の係数は $4/3$ であり、ブーストした偏極の時間成分に一致する。この時間成分に別の生成演算子を置くと、同じ拘束を二度数える。
:::

### Proca の時間順序積に現れる瞬時項

運動量空間の逆核 :eqref[eq-qft-proca-propagator] の時間成分は
$$
D_{\rm Proca}^{00}(k)
=\frac{i|\vb k|^2/m^2}{(k^0)^2-|\vb k|^2-m^2+i0}+\frac{i}{m^2}.
$$
この分解は $(k^0)^2-m^2=((k^0)^2-E^2)+|\vb k|^2$ を使った。第一項は三つの自由振動子の時間順序積から得られる。実際、時間成分を持つ縦偏極の二乗は $|\vb k|^2/m^2$ である。第二項は $i\delta^4(x-y)/m^2$ に対応する局所項で、源の逆核を使うときに必要である。これを Hamiltonian で確認するには、外部源のある $A^0$ に関する部分
$$
A^0(\nabla\cdot\boldsymbol\pi+\rho)-\tfrac12m^2(A^0)^2
$$
を消去する。解は $A^0=(\nabla\cdot\boldsymbol\pi+\rho)/m^2$ であり、代入した密度に $\rho^2/(2m^2)$ が加わる。これが時間成分の源どうしの瞬時項を作る。混合 $0i$ 成分では正・負の時間枝で $p^i$ の符号が入れ替わるので、時刻積分の分子に $k^0 k^i$ が出る。空間 $ij$ 成分の分子は $k^0$ に依らず、自由振動子の時間順序積と同じである。

:::example[静的な源で二つの項を足す]
$k^0=0$ と取ると動的な部分は $-i|\vb k|^2/[m^2(|\vb k|^2+m^2)]$。瞬時の $i/m^2$ を足すと $i/(|\vb k|^2+m^2)$ になり、静的 Proca 方程式の Yukawa 核を再現する。特に $\vb k\to0$ では動的部分は零へ向かうが、静的核は $i/m^2$ に近づく。ポテンシャルの時間順序積だけを源の交換に使うと、一様に近い源でこの値を失う。
:::

### 保存しない源を入れたときの発散

Proca 方程式の発散は $m^2\partial\cdot A=\partial\cdot j$。従って一般の源には
$$
(\Box+m^2)A^\nu=j^\nu+m^{-2}\partial^\nu(\partial\cdot j)
$$
が成り立つ。遅延解を選び、両端で同じ境界条件を使えば、発散の右辺が零でなくても元の拘束を満たす。これは右辺の発散が $m^{-2}(\Box+m^2)\partial\cdot j$ であることから分かる。任意の二源の伝播関数の縮約は
$$
j_\mu D_{\rm Proca}^{\mu\nu}j'_\nu
=\frac{-i}{k^2-m^2+i0}
 \left(j\cdot j'-\frac{(k\cdot j)(k\cdot j')}{m^2}\right).
$$
従って $m\to0$ の振る舞いは、各源がどのように保存条件に近づくかに依存する。例えば $j^\mu=k^\mu f,j'^\mu=k^\mu g$ の純縦の源では、分子に $(k^2-m^2)$ が掛かって極が消え、縮約は $i k^2 f g/m^2$ になる。殻上の三偏極には $k\cdot\epsilon_\lambda=0$ により結合しないが、源の間には局所な項が残る。逆に保存源では $m^{-2}$ の全項が消え、$k^2\ne0$ で滑らかに共変 Maxwell の核に近づく。

### Stueckelberg の逆核と補助極の相殺

ゲージ固定後の $A_\mu$ の二次演算子を、$k^2\ne0$ で $P_T,P_L$ に分ける。固有値は $-(k^2-m^2)$ と $m^2-k^2/\xi$。従って
$$
\begin{aligned}
D_\xi^{\mu\nu}(k)
&=-i\left[\frac{P_T^{\mu\nu}}{k^2-m^2+i0}
      +\frac{\xi P_L^{\mu\nu}}{k^2-\xi m^2+i0}\right]\\
&=\frac{-i}{k^2-m^2+i0}
 \left[\eta^{\mu\nu}-(1-\xi)\frac{k^\mu k^\nu}{k^2-\xi m^2+i0}\right].
\end{aligned}
$$
$\sigma$ の核は $i/(k^2-\xi m^2+i0)$ であり、両場の混合はこのゲージで消える。ゲージ不変な $\widehat A^\mu=A^\mu-\partial^\mu\sigma/m$ の二点の共変核では、両端の微分が $k^\mu k^\nu/m^2$ を掛ける。従って
$$
D_{\widehat A}^{\mu\nu}
=D_\xi^{\mu\nu}
  +\frac{i k^\mu k^\nu}{m^2(k^2-\xi m^2+i0)}
=D_{\rm Proca}^{\mu\nu}.
$$
縦の射影で見ると、足した二項は $i P_L(k^2/m^2-\xi)/(k^2-\xi m^2+i0)=iP_L/m^2$。補助的な極が消え、元の Proca の縦逆核になる。この計算では時間微分の接触項を含む共変核をそろえている。質量付きの三物理偏極とは別に $\sigma$ の極を足すと、同じゲージを違う変数で記述した自由度を二重に数える。

:::example[二つのゲージの極を比較する]
$m=3\,\mathrm{GeV}$、$\xi=4$ では $\widehat A$ を作る前の $D_\xi$ の縦の極は $36\,\mathrm{GeV}^2$、横の極は $9\,\mathrm{GeV}^2$。$\xi=1$ なら両極の位置は同じ $9\,\mathrm{GeV}^2$ になって成分表示は単純だが、$\widehat A$ にはどちらの場合も物理偏極の留数だけが残る。例えば $k_\mu D_{\widehat A}^{\mu\nu}=i k^\nu/m^2$ は極を持たない。時空で自由場の発散が零であることと、源を含んだ共変核の接触項が非零であることを、この縮約で区別できる。
:::

## 偏極ベクトルと偏極和

### 完全性による和

$m>0$ の質量殻で $p\cdot\epsilon_\lambda=0$、$\epsilon_\lambda\cdot\epsilon_\sigma^*=-\delta_{\lambda\sigma}$ の三偏極を選ぶ。静止系で時間成分は零、空間成分は通常の正規直交基底である。従って静止系で和は $\diag(0,1,1,1)$ になり、Lorentz 変換すれば
$$
\sum_{\lambda=1}^3\epsilon_\lambda^\mu\epsilon_\lambda^{*\nu}
=-\eta^{\mu\nu}+p^\mu p^\nu/m^2
\label{eq-qft-massive-polarizations}
$$
を得る。例えば $p=(5,0,0,4)\,\mathrm{GeV}$、$m=3\,\mathrm{GeV}$ では $\epsilon_L=(4/3,0,0,5/3)$ が規格化された縦偏極である。

光子については $p^2=0$ であるため $p$ 自身が $p^\perp$ に入る。$p\cdot n\ne0$ の時間的な $n$ を選び、$p,n$ の双方と直交する二偏極を使う。和は
$$
-\eta^{\mu\nu}+\frac{p^\mu n^\nu+n^\mu p^\nu}{p\cdot n}
-\frac{n^2p^\mu p^\nu}{(p\cdot n)^2}
$$
である。両辺は $p,n$ を縮約すれば零で、二次元の残りの空間上で負の計量の逆を与えるため一致する。$n=(1,\vb0)$ を選ぶと空間ブロックは :eqref[eq-qft-transverse-bracket] の射影になる。保存電流で両側を挟む場合は $p$ を含む項が消える。

### 縦偏極と保存電流

$p=(E,0,0,k)$ で $\epsilon_L=(k,0,0,E)/m$ とする。保存電流なら $EJ^0-kJ^3=0$ だから $k\ne0$ で
$$
J\cdot\epsilon_L=\frac{k J^0-EJ^3}{m}=-\frac{m}{k}J^0.
$$
従って $J^0$ と $k$ を固定した $m\to0$ 極限では縦偏極への結合が消える。例えば $E=5,k=4,m=3\,\mathrm{GeV}$、$J=(4,0,0,5)$ なら縮約は $-3$ である。これに対し $J=(1,0,0,0)$ は保存条件を満たさず、縮約は $k/m$。質量ゼロ極限での違いは、偏極の規格化だけでは判断できない。

### 偏極和を射影演算子として構成する

質量付きの $p^\perp$ に任意の四元ベクトル $v$ を射影するには、$p$ に沿う成分を引き $v_\perp=v-p(p\cdot v)/m^2$ とする。三偏極の計量内積は $-\delta_{\lambda\sigma}$ なので、その基底での係数は $-\epsilon_\lambda^*\cdot v$。従って
$$
v^\mu-\frac{p^\mu(p\cdot v)}{m^2}
=-\sum_{\lambda=1}^3\epsilon_\lambda^\mu(\epsilon_\lambda^*\cdot v).
$$
任意の $v$ で成り立つため、下添字で縮約する逆計量を表に出すと :eqref[eq-qft-massive-polarizations] が得られる。偏極和そのものを $R^{\mu\nu}$ と書くと、混合添字 $-R^\mu{}_\nu$ が冪等である。二つの反変添字のままユークリッドの行列積を取ると、この冪等性の判定を誤る。

質量ゼロで $n$ を時間的とする。$v_\perp=v-a p-b n$ として $p\cdot v_\perp=n\cdot v_\perp=0$ を解く。$p^2=0$ より $b=(p\cdot v)/(p\cdot n)$、$a=(n\cdot v)/(p\cdot n)-n^2(p\cdot v)/(p\cdot n)^2$。従って $-v_\perp$ の係数から、前の光子偏極和の三項が得られる。ここでは $p,n$ の張る二次元平面の Gram 行列の行列式が $-(p\cdot n)^2$ で、$p$ が零でなければ可逆である。この射影の像は負定値な二次元空間であり、$n$ の静止系で普通の三空間の横面として見える。

:::example[時間軸を取り直しても同じ光子を表す]
$p=(E,0,0,E)$ とし、最初は $n=(1,0,0,0)$、$\epsilon=(0,1,0,0)$ を選ぶ。$n'=(\gamma,\gamma v,0,0)$ に変えると $\epsilon\cdot n'=-\gamma v$、$p\cdot n'=\gamma E$。$\epsilon'=\epsilon+ (v/E)p=(v,1,0,v)$ とすれば $n'\cdot\epsilon'=0$ となる。$p\cdot\epsilon'=0$、$\epsilon'^2=-1$ は $p^2=0$ により保たれる。保存電流 $j=(a,b,c,a)$ との縮約は両方とも $-b$。個々の時間成分は変わるが、局所な平面波の場の強さ $p\wedge\epsilon$ と放射振幅は同じである。
:::

### 二つの実偏極と円偏極の密度行列

$p$ が $+z$ 方向なら横面に $e_x,e_y$ を取れる。$e_\pm=(e_x\pm i e_y)/\sqrt2$ は互いに Hermite 内積で直交する。$z$ まわりの回転は $e_\pm\mapsto e^{\mp i\theta}e_\pm$ なので、ヘリシティはそれぞれ $\pm1$。電子や他の荷電源が作る横振幅を $\vb J_T$ とすると、振幅 $M_\pm=e_\pm^*\cdot\vb J_T$ から偏極の密度行列 $\rho_{\lambda\sigma}=M_\lambda M_\sigma^*$ を定められる。純粋な一光子の偏極状態なら、全確率で割った $\rho$ は半正定値で階数一、トレース一である。未観測の源の状態も和に取ると、異なる $M$ の外積の和になり、階数が増えることがある。偏極を測定しないときに使う和は $\sum_\lambda|M_\lambda|^2=\vb J_T^*\cdot\vb J_T$ である。

:::example[線偏極の位相情報と円偏極の混合を比較する]
$x$ 向きの線偏極では $e_x=(e_++e_-)/\sqrt2$ だから円偏極の測定では両方の確率が $1/2$。一方、正負の円偏極を古典的に半分ずつ混ぜても同じ二確率になる。両者を区別するには、再び $x,y$ 基底で測る。線偏極なら $x$ が確率一、混合では $x,y$ がともに $1/2$。密度行列の対角成分が一致しても、線偏極の円基底では非対角成分 $\rho_{+-}=1/2$ が残るので、別の基底で同じ確率にはならない。
:::

### 光子の偏極和を $-\eta$ に置く条件

外部光子に接続する振幅を $M_\lambda=\epsilon_\lambda^{*\mu}\mathcal J_\mu$ とする。$p^\mu\mathcal J_\mu=0$ が成り立つなら、
$$
\sum_\lambda|M_\lambda|^2
=-\eta^{\mu\nu}\mathcal J_\mu\mathcal J_\nu^*
$$
となる。証明は光子の偏極和の $p$ を含む各項が両端の縮約で消えることである。この右辺は、条件を課した空間で正である。例えば $p=E(1,0,0,1)$ では $\mathcal J^0=\mathcal J^3$ となり、右辺は $|\mathcal J^1|^2+|\mathcal J^2|^2$ である。従って時間成分の負ノルムを物理的な偏極として足したわけではない。電流保存が一つ一つの寄与で成り立たず、足した振幅でだけ成り立つこともある。この場合、各項の絶対値を先に取ると保存条件を使えず、干渉を含めた和に対して偏極和を取る必要がある。

:::example[保存しない源には共変な置換を使えない]
同じ光的 $p$ で $\mathcal J^\mu=(1,0,0,0)$ とする。$n=(1,\vb0)$ の物理光子への縮約は二つとも零だが、$-\eta$ だけを使うと $-1$ が出る。ここで $\mathcal J\cdot p=E\ne0$ だから、$n$ を変えるゲージ成分を加えると振幅そのものが変わる。負の値は光子の観測確率という意味を持たず、縮約の条件を落とした結果である。
:::

### 質量付きの偏極分解と保存流の光的極限

$m>0$ で $p=(E,0,0,k)$、$E^2=k^2+m^2$ とする。$p\cdot J=0$ なら $J^3=E J^0/k$。三偏極の全和は
$$
-\eta_{\mu\nu}J^\mu J^{*\nu}
=|J^1|^2+|J^2|^2+\frac{m^2}{k^2}|J^0|^2
$$
である。最後の項は $\epsilon_L$ との縮約 $|J\cdot\epsilon_L|^2$ に等しい。$|k|>0$ と有限の $J$ を固定して $m\to0$ とすれば、右辺が光子の横二成分に近づく。極限の途中で $k$ も零にすると $m^2/k^2$ の見積りは使えない。例えば静止系の保存条件は $J^0=0$ だけであり、空間三成分は対等に偏極と結合する。空間運動量を固定した極限と静止系を毎回選んだ極限は、調べている外部源と運動量の族が異なる。

### 場の強さが選ぶ円偏極と双対

実波を最後に作る前に、複素な正周波数平面波 $A^\mu=(0,e_\lambda)e^{-ipx}$ を考える。$\vb E=i\omega e_\lambda e^{-ipx}$、$\vb B=i\,\vb p\times e_\lambda e^{-ipx}$。$\vb p=\omega\vb e_z$、$\vb e_z\times e_\lambda=-i\lambda e_\lambda$ なので $\vb B=-i\lambda\vb E$ となる。複素共役の負周波数側ではこの関係の $i$ の符号が逆になる。実電磁場は両側を足したものであり、実数の $\vb E,\vb B$ を複素係数で互いに比例させるという意味ではない。各複素波の反対称テンソルは双対作用の固有空間の一方に属し、空間反転は二つの円偏極を交換する。これは [Lorentz 群と Poincaré 群の表現](/physics/quantum-field-theory/lorentz-poincare/) の双対分解を、自由運動方程式の解に適用したものである。

:::example[重ね合わせのビームの時間平均]
同じ周波数の正周波数振幅 $\vb E_0=E_0 e_+$ に対し $\vb B_0=-iE_0 e_+$。実場を $\Re(\vb E_0 e^{-ipx})$、$\Re(\vb B_0 e^{-ipx})$ と取ると、$|\vb E|^2=|\vb B|^2=|E_0|^2/2$、$\vb E\cdot\vb B=0$。例えば $E_0$ が実で時空位相 $p x=\pi/2$ なら $\vb E=E_0\vb e_y/\sqrt2$、$\vb B=-E_0\vb e_x/\sqrt2$ である。Poynting 流は正の $z$ 方向で、その大きさ $|E_0|^2/2$ はエネルギー密度と一致する。偏極は回転するが、伝播の向きは波の運動量で決まる。
:::

### 偏極を数えて自由光子の熱エネルギーを求める

自由光子の周期箱で零点を引き、温度 $T>0$、$k_{\rm B}=1$ とする。角周波数 $p=|\vb p|>0$ ごとに規格化した一振動子の占有数は $n\ge0$ で、分配関数の幾何級数は $1/(1-e^{-p/T})$。平均は $n_B(p)=1/(e^{p/T}-1)$。横二偏極なので単位体積当たりの熱エネルギーは
$$
u(T)=2\int\frac{\dd^3p}{(2\pi)^3}\frac{p}{e^{p/T}-1}
     =\frac{T^4}{\pi^2}\int_0^\infty\dd x\,\frac{x^3}{e^x-1}.
$$
正項の幾何級数 $\sum_{r\ge1}e^{-rx}$ を入れて積分すると $6\sum_{r\ge1}r^{-4}=\pi^4/15$ となる。従って $u=\pi^2 T^4/15$。方向平均では光子の運動量流に $p_i p_j/p$ が入り、$\langle p_i p_j\rangle=p^2\delta_{ij}/3$ なので圧力は $u/3$。例えば $T$ を二倍にすると $u$ は $16$ 倍になる。温度のある状態は静止系を選ぶため、真空テンソルの Lorentz 不変性を熱平均へそのまま課すことはできない。共変ゲージの四振動子をそのまま熱的な物理状態として数えると、この係数を二倍に誤る。

### 箱の運動量と偏極の独立性

二つの異なる運動量を持つ光子の一粒子状態は箱の Fourier 直交性で分かれる。同じ $\vb p$ でも直交偏極は $[a_\lambda(\vb p),a_\sigma^\dagger(\vb p)]=\delta_{\lambda\sigma}$ により区別される。例えば一つの $\vb p$ に光子二個を入れると横二モードに対して $(2,0),(1,1),(0,2)$ の三状態がある。$e_+$ の生成演算子を $a_+^\dagger=(a_x^\dagger+i a_y^\dagger)/\sqrt2$ とすれば、
$$
\frac{(a_+^\dagger)^2}{\sqrt2}|0\rangle
=\frac12 |2,0\rangle+\frac{i}{\sqrt2}|1,1\rangle-\frac12 |0,2\rangle .
$$
従って線偏極で二個を数える測定では、各確率は $1/4,1/2,1/4$。どの測定結果でも全エネルギーは $2|\vb p|$、全運動量は $2\vb p$ である。偏極を選び直すユニタリ変換が自由 Hamiltonian と可換な理由は、同じ波数の二偏極が同じ周波数を持つことにある。

### 横偏極の位相を運動量ごとに選び直す

$\epsilon_\lambda(p)$ の選び方は固定した $p$ の横面の正規直交基底を選ぶことに相当する。例えば $\epsilon'_\lambda=\sum_\sigma\epsilon_\sigma U_{\sigma\lambda}(p)$ で $U$ がユニタリなら、$a'_\lambda=\sum_\rho U_{\rho\lambda}^*a_\rho$ と変えれば $\sum_\lambda\epsilon'_\lambda a'_\lambda=\sum_\lambda\epsilon_\lambda a_\lambda$ である。新しいモードの交換関係では運動量デルタが $p=q$ を選ぶため、係数に $U^\dagger(p)U(p)=1$ が残り、規格化が保たれる。偏極和も $\sum_\lambda U_{\sigma\lambda}U_{\rho\lambda}^*=\delta_{\sigma\rho}$ により不変である。従って波束で複数の運動量を重ねても、基底の位相を吸収したモードで同じ局所場を構成できる。

例えば基準光線を回して線偏極を定義する方法では、その回転に光線の周りの回転を足せる。角 $\alpha(p)$ の回転で円偏極は $\epsilon_\pm'(p)=e^{\mp i\alpha(p)}\epsilon_\pm(p)$。このとき $a_\pm'(p)=e^{\pm i\alpha(p)}a_\pm(p)$ なので、消滅側の積は保たれる。個々の偏極ベクトルの座標表示がある方向で不連続でも、そこで別の基底を使って場と偏極和をつなげられる。基底を選び直すこの自由度と、四元偏極に $p^\mu$ に比例する項を足すゲージの自由度は操作が異なる。前者は横面の中でのユニタリ変換、後者は物理商で零と同一視する方向への移動である。

## 演習問題

::::exercise{#exr-qft-vector-fields-1 level=basic}
$\vb p$ が $z$ 軸なら横射影の行列は何か。

:::solution
$\widehat{\vb p}=(0,0,1)$ なので $1-\widehat{\vb p}\widehat{\vb p}^{\mathsf T}=\diag(1,1,0)$ である。
:::
::::

::::exercise{#exr-qft-vector-fields-2 level=basic}
自由質量ゼロ場の偏極数を求めよ。

:::solution
三次元空間の横部分は $\vb p\ne0$ で二次元なので二つである。
:::
::::

::::exercise{#exr-qft-vector-fields-3 level=basic}
$m=3,p_z=4,E=5$ で $\epsilon_L^2$ を求めよ。

:::solution
$\epsilon_L=(4/3,0,0,5/3)$ から $(16-25)/9=-1$ である。
:::
::::

::::exercise{#exr-qft-vector-fields-4 level=basic}
Coulomb ゲージで静止点電荷の $A^0$ を求めよ。

:::solution
球対称で原点を除く Laplace 方程式を満たし、無限遠で消える解は $c/r$ である。$-\nabla A^0$ の微小球を通る流束が $q$ なので $c=q/(4\pi)$ となる。
:::
::::

::::exercise{#exr-qft-vector-fields-5 level=advanced}
横射影が冪等であることを示せ。

:::solution
$u^{\mathsf T}u=1$ として $(1-uu^{\mathsf T})^2=1-2uu^{\mathsf T}+uu^{\mathsf T}=1-uu^{\mathsf T}$ である。
:::
::::

::::exercise{#exr-qft-vector-fields-6 level=advanced}
Gupta–Bleuler の縦・時間の零ノルムを計算せよ。

:::solution
規格化因子を一とすればノルム二乗は $\mel{0}{(a^0-a^3)(a^{0\dagger}-a^{3\dagger})}{0}=-1+1=0$ である。交差項は対角計量で零である。
:::
::::

::::exercise{#exr-qft-vector-fields-7 level=advanced}
保存源に対し Proca 方程式が $\partial\cdot A=0$ を含むことを示せ。

:::solution
発散の第一項は二階微分が対称で $F$ が反対称なので零になる。$\partial\cdot j=0$ なら $m^2\partial\cdot A=0$ であり、$m>0$ から従う。
:::
::::

::::exercise{#exr-qft-vector-fields-8 level=advanced}
:eqref[eq-qft-massive-polarizations] を計量で縮約せよ。

:::solution
$\eta_{\mu\nu}(-\eta^{\mu\nu}+p^\mu p^\nu/m^2)=-4+1=-3$。これは各偏極のノルム $-1$ の三項の和である。
:::
::::

::::exercise{#exr-qft-vector-fields-9 level=advanced}
Proca 場の自由な縦モードの座標と運動量 $A_L,\pi_L$ に対し、Hamiltonian が $\frac12(1+p^2/m^2)\pi_L^2+\frac12m^2A_L^2$ である。標準振動子の正準変数を作り、基底の $\langle A_L^2\rangle$ を求めよ。$m=3,p=4\,\mathrm{GeV}$ の値を答えよ。

:::solution
$E=\sqrt{p^2+m^2}$、$C=E/m$ と置くと運動量の二次係数は $C^2$。$Q=A_L/C$、$P=C\pi_L$ とすれば $P\delta Q=\pi_L\delta A_L$ なので正準である。代入すると $H=P^2/2+E^2Q^2/2$。従って $\langle Q^2\rangle=1/(2E)$、$\langle A_L^2\rangle=C^2/(2E)=E/(2m^2)$ となる。指定の値では $5/18\,\mathrm{GeV}^{-1}$ である。単位は一つの規格化済み箱モードの座標に対するもので、元の位置空間の場へ戻すときは箱のモード関数を掛ける。

固定した $p>0$ で $m\to0$ とすると $A_L$ の分散は発散する。しかし $m>0$ での正規化した座標 $Q$ の分散は有限である。また保存電流との縮約には $m$ の抑制がある。従ってポテンシャルの成分の大きさからだけでは物理的な結合を判断できない。
:::
::::

::::exercise{#exr-qft-vector-fields-contact-poles level=advanced}
$m>0$ の Proca 場に対する空間的に一様な時間成分の外部源を考える。この源は一般には保存しない。源の周波数 $k^0\ne0$ に対し、$D_{\rm Proca}^{00}(k^0,\vb0)$ と $F_{\mu\nu}$ の応答を求めよ。静止モードの物理偏極の時間成分と矛盾しない理由を示せ。

:::solution
$\vb k=0$ では :eqref[eq-qft-proca-propagator] の $00$ 成分は
$$
\frac{-i}{(k^0)^2-m^2+i0}
\left(1-\frac{(k^0)^2}{m^2}\right)=\frac{i}{m^2}.
$$
分子の零が分母の極を消すため、$k^0=\pm m$ も含めて分布は定数である。方程式でも時間成分は $m^2 A^0=\rho$、空間成分は無入射の条件で $\vb A=0$ となる。従って $F_{0i}=-\partial_i A_0=0$、$F_{ij}=0$。発散条件は $m^2\partial_t A^0=\partial_t\rho$ である。正の周波数の自由場で $A^0=0$ なのは、対応する源なしの拘束が $m^2 A^0=0$ だからである。この源つきの応答では自由振動の極がなく、時間成分は拘束から直接決まる。
:::
::::

::::exercise{#exr-qft-vector-fields-linear-polarizations level=advanced}
光子が $p=(E,0,0,E)$ を持ち、二つの保存電流から同じ終状態へ放射されるとする。複素振幅としての電流を $J=(2,1,0,2)$、$J'=(1,0,i,1)$ とし、相対位相 $\beta$ を付けて $J+e^{i\beta}J'$ とする。円偏極ごとの強度と、その和を求めよ。

:::solution
どちらの電流も $p\cdot J=p\cdot J'=0$。横成分は $(1,i e^{i\beta})$ だから、振幅の全体の負号を除き $e_\pm^*\cdot\vb J_T=(1\pm e^{i\beta})/\sqrt2$。従って強度は $1\pm\cos\beta$、和は $2$。例えば $\beta=0$ では正側が $2$、負側が零、$\beta=\pi$ では逆になる。時間成分と縦成分は等しいため $-\eta$ を用いた総和では互いに消える。一方、二つの電流の全過程を非干渉に足すと、どの $\beta$ でも各円偏極に一ずつとなり、偏極ごとの位相依存が失われる。
:::
::::

## 確認テスト

::::quiz
Coulomb ゲージで $A^0$ の方程式を選べ。

- [x] $-\nabla^2A^0=\rho$
- [ ] $\partial_t^2 A^0=0$
- [ ] $A^0=\rho$

:::explanation
$\nabla\cdot\vb E=\rho$ に $\vb E=-\nabla A^0-\dot{\vb A}$ と $\nabla\cdot\vb A=0$ を代入すると得られる。
:::
::::

::::quiz{answer="3" tol="1e-6"}
Proca 場の偏極数を答えよ。

:::explanation
静止系では $p\cdot\epsilon=m\epsilon^0=0$ である。残る空間成分が三つある。
:::
::::

::::quiz
保存電流で光子のゲージ依存項が消える理由を選べ。

- [x] 分子が $k^\mu k^\nu$ に比例する
- [ ] 伝播関数の極が存在しない
- [ ] 電流の全成分が零である

:::explanation
電流保存は $k_\mu j^\mu=0$ である。従って端点のどちらで縮約しても縦項は消える。
:::
::::

::::quiz{answer="2" tol="1e-6"}
三次元横射影のトレースを答えよ。

:::explanation
$\tr(1-uu^{\mathsf T})=3-u^{\mathsf T}u=2$ である。
:::
::::

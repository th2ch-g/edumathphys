---
title: "曲がった時空の場の量子論入門"
description: "場の正準量子化から観測者による粒子の違いを導き、地平線、膨張、くりこみの計算を扱う。"
notation:
  - symbol: "(f,g)_{\\rm KG}"
    meaning: "Klein–Gordon 内積"
  - symbol: "\\alpha_{ij},\\beta_{ij}"
    meaning: "Bogoliubov 係数"
  - symbol: "\\xi"
    meaning: "スカラー場の曲率結合"
  - symbol: "\\eta"
    meaning: "宇宙の共形時刻"
---

## この章の目標

- スカラー場のモードを保存内積で規格化する。
- Bogoliubov 変換によって真空と粒子数の関係を求める。
- Unruh 効果と Hawking 放射の温度を零座標から導く。
- 膨張宇宙の粒子生成と局所演算子のくりこみを区別する。

導出では $c=1$ とし、この章の場の規格化には $\hbar=1$ も用いる。温度の式では $\hbar,k_B$ を戻す。背景計量を古典的なものとして与え、場の量子状態をその上で解く。

## 曲がった時空のスカラー場とモード展開

大域的双曲な背景 $g$ 上の実スカラー場に

$$
S=-\tfrac12\int\sqrt{-g}\bigl(g^{ab}\partial_a\phi\partial_b\phi+(m_s^2+\xi R)\phi^2\bigr)\dd^4x
$$

を用いる。ここで $m_s$ は逆長さで表した場の質量である。部分積分により $(\Box_g-m_s^2-\xi R)\phi=0$ を得る。複素化した二解 $f,g_1$ に対して $j^a=i(f^*\nabla^ag_1-g_1\nabla^af^*)$ は、同じ実ポテンシャルを満たすため $\nabla_aj^a=0$ である。未来法線 $n$ を使い :term[Klein–Gordon 内積]{en="Klein–Gordon inner product" yomi="くらいんごるどんないせき"} を

$$
(f,g_1)_{\rm KG}=i\int_\Sigma(f^*\nabla_n g_1-g_1\nabla_n f^*)\dd\Sigma
\label{eq-gr-kg-product}
$$

と取る。Cauchy 面間の発散定理と空間の境界流束がないことにより、これは面に依存しない。複素解全体で正定値ではなく、平坦時空の $e^{-i\omega t}$ には正、共役には負のノルムを与える。

正ノルムのモード $u_i$ に $(u_i,u_j)=\delta_{ij}$、$(u_i,u_j^*)=0$ を課し、$u_i,u_i^*$ が完全とする。場の展開 $\widehat\phi=\sum_i(a_i u_i+a_i^\dagger u_i^*)$ では $a_i=(u_i,\widehat\phi)_{\rm KG}$ である。等時刻の場と正準運動量の交換関係をこの積に入れると $[a_i,a_j^\dagger]=(u_i,u_j)_{\rm KG}=\delta_{ij}$ になる。有限体積の周期境界で Minkowski のモードを $u_{\vb k}=e^{i\vb k\cdot\vb x-i\omega t}/\sqrt{2\omega V}$ と取れば、時間微分の差は規格化が一致するモードで $2\omega$ を、空間積分は $V\delta_{\vb k,\vb k'}$ を与える。静的背景では、適切な境界条件で空間作用素が正のときに時間並進で正周波数を選べる。定常場が空間的になるエルゴ領域や、空間作用素に負の固有値がある場合には、この選び方で正ノルムの基底を得られるとは限らない。一般の計量にはそもそも共通の時間並進がない。:ref[def-aqm-fock]{full} の表現を組み立てる一粒子空間はモードの選択に依存する。

### 静的時空でモード問題を自己共役形にする

静的計量 $-\!N^2\dd t^2+h_{ij}\dd x^i\dd x^j$ で $N>0$ とする。波動方程式に $\phi=e^{-i\omega t}f(\vb x)$ を入れると
$$
\omega^2 f=\mathcal A f,\qquad
\mathcal A=-\frac{N}{\sqrt h}\partial_i(N\sqrt h\,h^{ij}\partial_j)
+N^2(m_s^2+\xi R)
$$
である。重み $\sqrt h/N$ で $\overline f\mathcal A g$ を積分し、境界で流束が消える条件を選べば、微分項は $\int N\sqrt h h^{ij}\partial_i\overline f\,\partial_jg$ になる。従ってこの作用素は対称であり、ポテンシャルが非負ならその二次形式も非負である。自己共役作用素として使うには関数空間と境界条件を指定する必要がある。正の固有値のモードでは KG ノルムは $2\omega$ とこの重みでの空間ノルムの積である。$\omega=0$ のモードではこの規格化が使えないので、例えば質量ゼロの有限周期箱では定数の自由度を別に量子化する。

:::example[反射端を持つ平坦な区間]
$0<x<L$ の両端で場を零とする。空間の作用素は $-\partial_x^2+m_s^2$ であり、固有関数 $f_n=\sqrt{2/L}\sin(n\pi x/L)$、$\omega_n^2=m_s^2+(n\pi/L)^2$ となる。境界項は両端で $f_n=0$ のため消える。場のモードを $f_n e^{-i\omega_n t}/\sqrt{2\omega_n}$ とすると、時間微分の差が $\omega_n+\omega_{n'}$ を与え、空間積分の $\delta_{nn'}$ によってノルムは $1$ となる。
:::

### Cauchy データから正準交換関係を導く

一般の時刻切片について、時間発展を $\partial_t=Nn+\beta^i\partial_i$ とする。スカラー場の運動項では $\partial_t\phi$ は $n\phi=N^{-1}(\partial_t\phi-\beta^i\partial_i\phi)$ の中だけに現れる。体積密度 $\sqrt{-g}=N\sqrt h$ を使うと正準運動量密度は $\pi_\phi=\sqrt h\,n\phi$ である。従って等時刻の交換関係は
$$
[\widehat\phi(\vb x),\widehat\pi_\phi(\vb y)]=i\delta^3(\vb x-\vb y),\qquad
[\widehat\phi,\widehat\phi]=[\widehat\pi_\phi,\widehat\pi_\phi]=0
$$
である。ここで $\delta^3$ は座標体積 $\dd^3x$ に対する分布とする。:eqref[eq-gr-kg-product] は
$$
a_i=i\int_\Sigma\bigl[u_i^*\widehat\pi_\phi-\sqrt h(nu_i^*)\widehat\phi\bigr]\dd^3x
$$
になる。共役で $i$ の符号を反転し、$\phi$ と $\pi_\phi$ の交差項だけを交換に残すと、$[a_i,a_j^\dagger]=i\int\sqrt h(u_i^*nu_j-u_jnu_i^*)=(u_i,u_j)_{\rm KG}$ となる。$[a_i,a_j]$ では $u_j$ が共役になり、正負の直交性から零である。

### 交換子は真空の選択に依存しない

場を実の古典解 $f$ に対して線形な演算子 $\widehat A_f=\int_\Sigma(f\widehat\pi_\phi-\sqrt h(nf)\widehat\phi)\dd^3x$ として読む。等時刻交換関係から
$$
[\widehat A_f,\widehat A_g]=i\int_\Sigma (f\,ng-g\,nf)\dd\Sigma
$$
である。右辺の反対称形式は波動方程式と発散定理で保存する。一般に $\widehat\phi(x)$ は演算子値分布であり、点での方程式や交換子は試験関数で積分した意味で使う。場の交換子 $[\widehat\phi(x),\widehat\phi(x')]=i\Delta(x,x')$ は、$x'$ を固定した波動方程式と、等時刻の値が零、第二点の未来法線微分が空間デルタ関数というデータで決まる。波動方程式の有限伝播速度により、空間的に離れた二点ではこの交換子は零である。モード基底で作る二点関数 $W(x,x')$ は、$\Delta$ による反対称部分を共有するが、対称部分は状態に依存する。

:::example[箱の零モードは振動子にならない]
質量ゼロの実場を体積 $V$ の周期箱に入れ、空間一定の成分を $\phi=q(t)/\sqrt V$ と書く。作用のこの部分は $\int\dot q^2\dd t/2$ なので $p=\dot q$、$H=p^2/2$ となる。解は $q(t)=q(0)+tp$ であり、交換子は $[q(t),q(t')]=i(t'-t)$ である。有限の正周波数を持つモードの定常真空をこの自由粒子の零モードに代入できない。非零の波数だけを使うか、零モードに独立な波束状態を指定すれば、この退化を明示して扱える。
:::

### 不安定な固有モードと正周波数

静的な空間作用素 $\mathcal A$ に負の固有値 $-\gamma^2$ が一つある場合、そのモードの時間発展は $q=c_+e^{\gamma t}+c_-e^{-\gamma t}$ である。実の指数関数それぞれの KG 内積でのノルムは零だが、両者の反対称形式は非零である。例えば空間固有関数を重み付きで規格化し、$q_\pm=e^{\pm\gamma t}$ とすれば $q_+\dot q_--q_-\dot q_+=-2\gamma$ である。この組を複素化して正ノルムの基底を作っても、時間並進はそのモードと共役を一般に混ぜる。実際 $\partial_t$ の固有モードは実指数であり、正のエネルギーで振動するモードではない。量子化の時間零のデータを与えることと、時間並進不変な基底状態を要求することの違いがここに現れる。

## Bogoliubov 変換と粒子の定義

別の正ノルムモードを $v_i=\sum_j(\alpha_{ij}u_j+\beta_{ij}u_j^*)$ とする。:term[Bogoliubov 変換]{en="Bogoliubov transformation" yomi="ぼごりゅーぼふへんかん"} の係数は $\alpha_{ij}=(u_j,v_i)$、$\beta_{ij}=-(u_j^*,v_i)$ である。正・負ノルムと直交性を代入すると $\alpha\alpha^\dagger-\beta\beta^\dagger=I$、$\alpha\beta^\top-\beta\alpha^\top=0$ を得る。新消滅演算子は $b_i=\sum_j(\alpha_{ij}^*a_j-\beta_{ij}^*a_j^\dagger)$ なので、旧真空 $a_i\ket{0_u}=0$ に対して

$$
\bra{0_u}b_i^\dagger b_i\ket{0_u}=\sum_j|\beta_{ij}|^2
\label{eq-gr-bogoliubov-number}
$$

である。例えば一モードで $\alpha=5/4,\beta=3/4$ なら $\abs{\alpha}^2-\abs{\beta}^2=1$ で、$v$ 粒子の平均数は $9/16$ である。連続モードの $\delta$ 規格化では波束を作ってから期待値を取る。モードが無限にある場合は全粒子数の和が発散することがあり、二つの Fock 表示が単一のユニタリで結ばれると仮定できない。

### 一モードの変換を真空の係数から解く

$[a,a^\dagger]=1$ とし、$b=\alpha a-\beta a^\dagger$、$\alpha,\beta$ は実で $\alpha^2-\beta^2=1,\alpha>0$ とする。$\ket{0_b}=\sum c_n\ket n_a$ を $b$ で消すと
$$
\alpha\sqrt{n+1}\,c_{n+1}=\beta\sqrt n\,c_{n-1},\qquad c_{-1}=0.
$$
従って奇数係数は零で、$c_{2k}=c_0(\beta/\alpha)^k\sqrt{(2k)!}/(2^k k!)$ となる。$z=(\beta/\alpha)^2<1$ とすると、二項展開 $\sum_{k\ge0}\binom{2k}{k}(z/4)^k=(1-z)^{-1/2}$ から $|c_0|^2=\sqrt{1-z}=1/\alpha$ である。粒子数の平均は級数を $z$ で微分して $2z$ を掛け、規格化すると $z/(1-z)=\beta^2$ となる。偶数だけ現れるのは一モード振動子での対生成を反映する。

:::example[真空の重なりと粒子数の比較]
$\alpha=5/4,\beta=3/4$ とすると、旧真空が新真空に含まれる確率は $|c_0|^2=4/5$ である。ちょうど二粒子の確率は $(4/5)(9/25)/2=18/125$ である。平均粒子数は $9/16$ で、二粒子への確率の二倍だけでは足りない。四粒子以上の重みも平均に寄与するためである。
:::

### 有限モードの変換を逆に解く

正・負周波数モードをまとめ、変換の中で $\alpha,\beta$ がどの添字に作用するかを固定する。列ベクトルに並べた旧演算子 $a,a^\dagger$ を使えば
$$
\begin{pmatrix}b\\ b^\dagger\end{pmatrix}
=\begin{pmatrix}\bar\alpha&-\bar\beta\\-\beta&\alpha\end{pmatrix}
\begin{pmatrix}a\\a^\dagger\end{pmatrix}
$$
である。ここで上線は成分ごとの複素共役である。交換多重行列を $J=\diag(I,-I)$ とすると、係数行列 $S$ は $SJS^\dagger=J$ を満たす。実際、上左ブロックは $\bar\alpha\alpha^{\mathsf T}-\bar\beta\beta^{\mathsf T}=I$、上右は $-\bar\alpha\bar\beta^{\mathsf T}+\bar\beta\bar\alpha^{\mathsf T}=0$ を使う。有限次元で $S^{-1}=JS^\dagger J$ である。従って旧演算子で新演算子を作る式と、新演算子で旧演算子を作る式の符号を交換関係から検査できる。無限モードではこの行列形式だけで全体の Fock 変換が実装できるとは限らず、係数の二乗総和を調べる必要がある。

### 多モードの真空を一つの式で作る

有限個の旧モードからなる列 $a$ を使い、$Z=(\bar\alpha)^{-1}\bar\beta$ とする。正準条件から $\alpha\alpha^\dagger=I+\beta\beta^\dagger$ で $\alpha$ は可逆である。また $\alpha\beta^{\mathsf T}=\beta\alpha^{\mathsf T}$ を左右から $\alpha^{-1}$ とその転置で掛けると $Z^{\mathsf T}=Z$ が従う。$I-ZZ^\dagger=(\bar\alpha)^{-1}[(\bar\alpha)^{-1}]^\dagger$ なので $Z$ の特異値は $1$ 未満である。従って $b$ に消される真空は
$$
\ket{0_b}=\mathcal N\exp\left(\tfrac12\sum_{jk}Z_{jk}a_j^\dagger a_k^\dagger\right)\ket{0_a}
$$
で表せる。実際、$a_i$ と指数の交換子が $\sum_jZ_{ij}a_j^\dagger$ であり、さらにこれが指数と可換なので $(a-Za^\dagger)\ket{0_b}=0$ となる。

複素対称行列 $Z$ はユニタリ基底を選んで $Z=U\operatorname{diag}(z_i)U^{\mathsf T}$、$0\le z_i<1$ と分解できる。これは対称行列の特異値分解である。その基底では各モードの指数が別々になり、一モードの級数規格化 $|c_0|^2=\sqrt{1-z_i^2}$ を積にして $|\mathcal N|^2=\det(I-ZZ^\dagger)^{1/2}$ を得る。例えば $Z_{12}=Z_{21}=z$ だけが非零の二モードの場合、指数は $z a_1^\dagger a_2^\dagger$ であり、$\ket{0_b}=\sqrt{1-|z|^2}\sum_{n\ge0}z^n\ket{n,n}$ となる。実際 $(a_1^\dagger a_2^\dagger)^n\ket{0,0}=n!\ket{n,n}$ である。単一モードの場合は同じモードに二個ずつ入り、ここでは異なる二モードに一個ずつ入る。

### 無限モードで有限個の直感が破れる例

互いに独立な振動子を可算個用意し、全てに同じ実係数 $\alpha>1,\beta\ne0$ の変換をかける。最初の $n$ 個を変換した真空と、最初の $n+m$ 個を変換した真空の内積の絶対値は $\alpha^{-m/2}$ である。共通の $n$ 個の因子は内積 $1$、加わった $m$ 個で旧真空との重なりが $\alpha^{-1/2}$ ずつだからである。従ってこの列は旧 Fock 空間のノルムで Cauchy 列にならない。独立な振動子で各 $\alpha_j>0$ の場合には $\sum_j\log\alpha_j=\frac12\sum_j\log(1+|\beta_j|^2)$ が有限であることと、真空の積が収束することが上の内積から分かる。この条件は $\sum_j|\beta_j|^2<\infty$ と同値である。

例えば箱の体積 $V$ を無限にし、スペクトルが $\sum_{\vb k}|\beta_{\vb k}|^2\simeq V\int\dd^3k\,|\beta_{\vb k}|^2/(2\pi)^3$ になる場合、右の積分が有限で正なら粒子密度は有限でも総数は無限になる。この極限では、箱の規格化で得た一個の波数の占有と、無限空間で正規化された波束に対する期待値を区別する。

### 一般の有界な変換を実装できる条件

可分な一粒子空間で、有界かつ有界な逆を持つ正準変換 $S$ を考える。$SJS^\dagger=S^\dagger JS=J$ の両式により、$\bar\alpha$ とその随伴はともに下から $1$ によってノルムを抑えられる。前者から閉値域、後者から稠密値域なので $\bar\alpha$ は可逆である。従って $Z=(\bar\alpha)^{-1}\bar\beta$ は有界で対称であり、有限の場合と同じ恒等式から $I-ZZ^\dagger\ge \|\alpha\|^{-2}I$ である。

変換が Fock 空間のユニタリで実装されるなら、新真空 $\Psi=\bigoplus_{n\ge0}\Psi_n$ が存在し、$(a-Za^\dagger)\Psi=0$ である。これを有限粒子試験ベクトルと内積を取って粒子数ごとに読む。もし $\Psi_0=0$ であれば、非零の最低粒子数 $n$ で全ての $a_i\Psi_n$ が零となる。$\sum_i\|a_i\Psi_n\|^2=n\|\Psi_n\|^2$ なので矛盾する。次に一粒子成分を読むと $a_i\Psi_2=\Psi_0\sum_j Z_{ij}\ket{1_j}$ である。従って
$$
|\Psi_0|^2\sum_{ij}|Z_{ij}|^2=2\|\Psi_2\|^2<\infty .
$$
この二乗可和性を :term[Hilbert–Schmidt 条件]{en="Hilbert–Schmidt condition" yomi="ひるべるとしゅみっとじょうけん"} という。有界な可逆係数 $\bar\alpha$ を左から掛けても条件が変わらないため、$\beta$ についても同値である。

逆にこの条件があれば、成分共役を $C$ として反線形作用素 $ZC$ の二乗 $ZZ^\dagger$ はコンパクトな正作用素である。正の固有値 $z_i^2$ の各固有空間で $ZC/z_i$ は等長で二乗が恒等になる。この共役で固定される実部分の正規直交基底を取ると $ZC$ は各基底ベクトルで $z_i$ 倍になる。核も基底に加えれば、有限モードで使った対称行列の分解を得る。同じ理由で有限の場合の分解も成り立つ。

基底変換と、その後に残る複素線形なユニタリの回転は、各粒子数のテンソル積に同じユニタリを作用させれば実装できる。実際、$Z$ が対角なら $a-Za^\dagger$ に $(I-ZZ^\dagger)^{-1/2}$ を掛けた列は正規化された独立な一モード変換であり、$b=\bar\alpha(a-Za^\dagger)$ とは左からのユニタリだけ異なる。一モードについて $a=(X+iP)/\sqrt2$ と置くと、$\alpha=\cosh s$、$\beta=\sinh s$ の変換は $X\mapsto e^{-s}X$、$P\mapsto e^sP$ である。位置表示で $U_s\psi(x)=e^{-s/2}\psi(e^{-s}x)$ とすればノルムが保存し、$U_s XU_s^{-1}=e^{-s}X$、$U_s P U_s^{-1}=e^sP$ となるので、これが一モードの実装である。

対角の一モードの係数は $\alpha_j=(1-z_j^2)^{-1/2}$ で、$\sup z_j<1$ と $\sum_j z_j^2<\infty$ から $\sum_j\log\alpha_j<\infty$ である。有限個ずつの一モードの実装を掛けた積は、有限個のモードにのみ励起がある状態上では、残りの因子が上で示した収束する真空の積になる。従ってこの稠密部分で強収束し、ノルムが一で一様有界なので全空間に拡張する。逆変換も同じ真空の重なりを持って強収束するため、極限同士は逆となる。こうして条件が十分であることも分かる。

## Unruh 効果

Minkowski 時空の右 Rindler 領域で $T=\rho\sinh(a\tau)$、$X=\rho\cosh(a\tau)$、$\rho>0$ とする。$\rho=1/a$ の曲線では $\tau$ が固有時で固有加速度は $a$ である。右領域の外向き零座標 $u=\tau-a^{-1}\log(a\rho)$ と Minkowski の $U=T-X$ の間に $U=-a^{-1}e^{-au}$ がある。無質量場の右モードは $e^{-i\omega u}=(-aU)^{i\omega/a}$ を $U<0$ に持つ。Minkowski の正周波数 $e^{-ikU}$、$k>0$ は $\Im U<0$ で減衰する。正周波数波束の境界値として右モードを下半平面から $U>0$ に回すと、$-U$ の偏角が $+\pi$ に変わるので係数 $e^{-\pi\omega/a}$ を伴い左モードの共役になる。Fourier 積分で負周波数がないことは、波束を下半平面へ閉じる積分で確かめられる。

従って Minkowski の正周波数の組合せには $\abs{\beta/\alpha}=e^{-\pi\omega/a}$ があり、$\abs\alpha^2-\abs\beta^2=1$ を用いると $\abs\beta^2=(e^{2\pi\omega/a}-1)^{-1}$ である。左右のモードの真空条件は $(b_R-e^{-\pi\omega/a}b_L^\dagger)\ket0=0$ と左右を入れ替えた式である。占有数基底で係数の漸化式を解くと同じ数 $n$ の左右成分が $e^{-\pi n\omega/a}$ の比で並ぶ。左の部分跡で $n$ の確率が $e^{-2\pi n\omega/a}$ に比例するため、右での温度は

$$
T_U=\frac{\hbar a}{2\pi k_B}.
$$

これが :term[Unruh 効果]{en="Unruh effect" yomi="うんるーこうか"} である。SI では $a$ を固有加速度 $\mathrm{m\,s^{-2}}$ として分母に $c$ が加わる。$a=9.8\,\mathrm{m\,s^{-2}}$ なら $T_U\simeq4.0\times10^{-20}\,\mathrm K$ である。Minkowski 真空の慣性観測者には粒子がないことと、一様加速する検出器の励起率が熱分布を示すことは、異なる時間発展を用いるため両立する。

### 左右のモードを同じ時間向きで規格化する

右領域の boost $\partial_\tau$ は Minkowski の未来を向くが、左領域では同じ boost ベクトルは過去を向く。右を $U<0<V$、左を $V<0<U$ とし、左で未来向きの座標を $T=\rho\sinh(a\tau_L)$、$X=-\rho\cosh(a\tau_L)$ と取る。左の $U=\rho e^{a\tau_L}$ に依存する未来正周波数の枝は $(aU)^{-i\omega/a}$ に比例する。従って右の $(-aU)^{i\omega/a}$ を $U$ の下半平面で延長した枝は、左の正周波数の共役を $z_\omega=e^{-\pi\omega/a}$ の係数で含む。

右と左の正ノルムモードを $w_R,w_L$ とすれば、この解析接続の組合せは $p_R=(w_R+z_\omega w_L^*)/\sqrt{1-z_\omega^2}$ である。右と左の支持が異なり、共役のノルムが負であるため分母は $\sqrt{1-z_\omega^2}$ になる。左から同じ操作をした $p_L=(w_L+z_\omega w_R^*)/\sqrt{1-z_\omega^2}$ も同様に規格化され、$(p_R,p_L)$ および $(p_R,p_L^*)$ は零である。後者では右と左の二つの交差項が逆符号で相殺する。連続モードは波束で扱う。二つの $p$ は慣性正周波数の組であり、その消滅条件は前節の有限二モード計算で $Z_{RL}=Z_{LR}=z_\omega$ を選ぶ場合に当たる。厳密な無限領域で全周波数の真空を左と右の通常の密度行列の積空間で書けると仮定する必要はない。有限モードの模型で得る Boltzmann 比は、局所検出器の相関関数でも確認できる。

### 検出器の時間相関

粒子の定義を局所測定に結び付けるため、二準位検出器を軌道 $x(\tau)$ 上で弱く場と結合する。準位差を $\omega>0$ とし、一次摂動の遷移振幅に検出器の内部時間変化 $e^{i\omega\tau}$ を入れる。場の終状態をすべて足して振幅の二乗を作ると、長時間の励起率は $\int_{-\infty}^{\infty}e^{-i\omega s}W(s)\dd s$ に比例する。ここで $W(s)=\langle\phi(x(\tau))\phi(x(\tau-s))\rangle$ であり、定常な軌道では $\tau$ に依存しない。

四次元の質量ゼロ Minkowski 真空では

$$
W(x,x')=-\frac1{4\pi^2[(T-T'-i0)^2-\abs{\vb X-\vb X'}^2]}.
$$

平面波積分で角度を積分し、$k$ の正周波数積分に減衰因子を入れるとこの境界値になる。双曲線軌道では時間差と空間差から $(T-T')^2-(X-X')^2=4a^{-2}\sinh^2(as/2)$ を得る。従って $W(s)=-a^2/[16\pi^2\sinh^2(a(s-i0)/2)]$ である。$\omega>0$ のとき積分を下半平面へ閉じると、極 $s=-2\pi i n/a+i0$、$n=1,2,\ldots$ が寄与する。各極では $W(s)\sim-1/(4\pi^2(s-s_n)^2)$ なので、時計回りに留数を足すと応答が $\omega/(2\pi)\sum_{n\ge1}e^{-2\pi n\omega/a}=\omega/[2\pi(e^{2\pi\omega/a}-1)]$ となる。なお厳密な積分では先に滑らかな切り替えで波束を作り、長時間の割合を取る。慣性軌道の真空では下半平面に同じ極列がなく、正の準位差の定常励起率は零である。

### 有限の観測時間と真空中の遷移率

検出器に滑らかで実の切り替え関数 $\chi(\tau)$ を掛け、結合を $\lambda\chi(\tau)\mu(\tau)\widehat\phi(x(\tau))$ とする。$\mu$ は検出器の演算子で、二準位の間の行列要素を $\mu_{eg}$ とする。場の終状態を完全系で足すと、励起確率の $\lambda^2|\mu_{eg}|^2$ で割った量は
$$
\mathcal F_\chi(\omega)
=\int\dd\tau\,\dd\tau'\,\chi(\tau)\chi(\tau')
e^{-i\omega(\tau-\tau')} W(x(\tau),x(\tau')) .
$$
符号は、励起の振幅に $e^{i\omega\tau}$ が入り、その二乗で場の二点関数の順序を上の形にそろえると時間の添字が入れ替わることから来る。慣性静止の真空で無質量なら $W(s)=\int_0^\infty k e^{-iks}\dd k/(4\pi^2)$ である。$\hat\chi(\nu)=\int e^{-i\nu\tau}\chi(\tau)\dd\tau$ と定めれば
$$
\mathcal F_\chi(\omega)=\frac1{4\pi^2}\int_0^\infty k\,|\hat\chi(\omega+k)|^2\dd k .
$$
従って有限の観測で慣性真空の励起が起こり得る。切り替えをゆっくりにして $\chi_T(\tau)=\chi(\tau/T)$ とすると、$\hat\chi_T(\nu)=T\hat\chi(T\nu)$ である。$\omega>0$ なら引数は常に $T\omega$ より大きく、滑らかなコンパクト台の Fourier 変換の急減衰により $\mathcal F_{\chi_T}(\omega)/T\to0$ である。この極限が真空中の定常励起率の零に対応する。

例えば $\omega>0$ の逆遷移で $\omega\mapsto-\omega$ とすると、積分核は $k=\omega$ の近くで幅 $\order{1/T}$ に集中する。$z=T(k-\omega)$ と置けば極限は $\omega\int|\hat\chi(z)|^2\dd z/(4\pi^2)=\omega\int|\chi(\tau)|^2\dd\tau/(2\pi)$ である。従って切り替え時間の規格化で割った自発遷移率は $\omega/(2\pi)$ になり、加速した軌道の詳細釣り合いにおける同じ項と一致する。

### 温度と詳細釣り合い

定常な検出器の遷移率を $\mathcal R(\omega)=\int e^{-i\omega s}W(s)\dd s$ とする。KMS 条件 $W(s-i\beta)=W(-s)$ の成り立つ帯の中で、適切に波束化した積分の輪郭を $-i\beta$ だけ移せるとする。すると指数から $e^{-\beta\omega}$ が出て、残る積分の $s$ を反転すれば
$$
\mathcal R(\omega)=e^{-\beta\omega}\mathcal R(-\omega)
$$
となる。Unruh の場合は $\beta=2\pi/a$ である。応答式の $\omega$ を負に続けると、$\omega>0$ について下向き率は $\omega/[2\pi(1-e^{-\beta\omega})]$ となり、励起率との差は $\omega/(2\pi)$ である。この真空でも残る項は、自発的な下向き遷移に対応する。例えば $\beta\omega=\log 3$ なら励起と下向き率の比は $1/3$、長時間で二準位の占有率は上準位 $1/4$、下準位 $3/4$ となる。

### 空間の運動量を残した Rindler モード

四次元の右領域で $\dd s^2=-a^2\rho^2\dd\tau^2+\dd\rho^2+\dd y^2+\dd z^2$ とし、$\phi=e^{-i\omega\tau+i\vb k_\perp\cdot\vb x_\perp}F(\rho)$ を入れる。体積密度 $a\rho$ から動径式は
$$
\rho^2 F''+\rho F'+\left(\frac{\omega^2}{a^2}-(m_s^2+k_\perp^2)\rho^2\right)F=0 .
$$
$\kappa_\perp=\sqrt{m_s^2+k_\perp^2}>0$ と置けば、大きい $\rho$ では $e^{\pm\kappa_\perp\rho}$ の指数が先頭になるため、減衰側を選ぶ。$\rho$ が小さい領域では質量と横運動量の項は $\rho^2$ に比例して消え、$F\sim\rho^{\pm i\omega/a}$ の二枝になる。これは $u=\tau-a^{-1}\log(a\rho)$ とその反対向きの零座標の位相である。従って質量や横運動量は領域全体でのモードの形を変えるが、地平線に十分近い対数位相の係数には現れない。ここでの $\omega$ は $\tau$ に関する周波数なので、別の $\rho$ の静止観測者の固有時間で測ると $\omega/(a\rho)$ となる。

:::example[同じ boost 時間を使う二観測者の温度]
$\rho_1=1/a,\rho_2=2/a$ に静止する二観測者を取る。boost 時間の周期が $2\pi/a$ なら、それぞれの固有時間での周期は $2\pi\rho_1,2\pi\rho_2$ である。従って第二の温度は第一の半分である。また加速度計での固有加速度は $1/\rho$ で、局所の温度 $1/(2\pi\rho)$ と比例する。異なる $\rho$ にある観測で同一の座標周波数を比較する際は、局所周波数への換算も同じ比で効く。
:::

## Hawking 放射の導出

崩壊で非退化な地平線が形成されると、未来無限遠へ届く遅い光を過去へ逆追跡する際に赤方偏移が大きくなる。外部で $u=t-r_*$ とし、地平線を横切る正則な零変数を $U=-C e^{-\kappa u}$ と取る。Schwarzschild では :eqref[eq-gr-kruskal] によりこの形で $\kappa=1/(4GM)$ である。滑らかな崩壊を通して過去無限遠のアフィン時刻 $v$ と結ぶと、最後に逃げる光線の $v_0$ の近傍で $v_0-v=Ae^{-\kappa u}+\order{e^{-2\kappa u}}$ となる。これは正則座標間に零でない一次の Jacobian がある場合の展開である。遅い時刻の波束を追えば $\omega$ の出射モードは先頭で $\Theta(v_0-v)[(v_0-v)/A]^{i\omega/\kappa}$ に比例する。

係数の大きさを比較するため $x=v_0-v>0$、$b=\omega/\kappa$ と置く。正負の過去周波数との Fourier 積分の核は

$$
\int_0^\infty x^{ib}e^{-(\eps\pm i\omega')x}\dd x
=\Gamma(1+ib)(\eps\pm i\omega')^{-(1+ib)},\qquad \omega'>0 .
$$

式は Gamma 積分の複素解析接続で、まず $\eps>0$ と取る。$\eps\to0^+$ では偏角 $\pm\pi/2$ による絶対値が $e^{\pm\pi b/2}$ であり、正周波数係数は $+i\omega'$ の核に対応する。従って負対正の絶対値の比は $e^{-\pi b}$ である。波束で正準規格化すると期待粒子数は $(e^{2\pi\omega/\kappa}-1)^{-1}$ となり、:eqref[eq-gr-hawking-temperature] に一致する。無限遠の各角モードには外側ポテンシャルの透過確率 $\Gamma_\ell(\omega)$ が掛かる。例えば質量ゼロの Schwarzschild スカラー場では $\phi=e^{-i\omega t}Y_{\ell m}\psi(r)/r$ を $\Box_g\phi=0$ に入れると $\psi''_{r_*}+[\omega^2-f(\ell(\ell+1)/r^2+f'/r)]\psi=0$ となり、ポテンシャルによる一部反射を計算できる。この :term[灰色体因子]{en="greybody factor" yomi="はいいろたいんし"} により遠方の流束は完全な黒体のものとは異なる。

### 薄い崩壊殻で零座標を接ぐ

入射する光的な殻を $v=v_s$ に置き、内側は平坦、外側は $m=GM>0$ の Schwarzschild とする。外の遅延時刻を $u=v-2r_*$、内の遅延時刻を $U_{\rm in}=T-r$ とする。平坦領域の入射時刻 $T+r$ を殻上で $v_s$ に合わせる。半径 $r_c$ で殻を横切る出射光線では $U_{\rm in}=v_s-2r_c$、$u=v_s-2r_*(r_c)$ である。最後に逃げ得る光線の内の値は $U_0=v_s-4m$ である。実際、殻上の $r_c=2m$ を未来にたどると外の地平線の生成子になる。

$x=U_0-U_{\rm in}=2(r_c-2m)$ として外の $r_*=r+2m\log(r/(2m)-1)$ に代入すれば
$$
u=U_0-x-4m\log\left(\frac{x}{4m}\right).
$$
従って $u\to+\infty$ では $x=4m e^{(U_0-u)/(4m)}+\order{e^{-u/(2m)}}$ となる。実際、指数を取れば $x e^{x/(4m)}=4m e^{(U_0-u)/(4m)}$ であり、左辺は $x=0$ の近傍で微分が $1$ である。平坦な中心で反射する動径光線には、中心で入射時刻と出射時刻が等しい。従って幾何光学のこの模型では $U_{\rm in}$ は過去無限遠から来るアフィン時刻へ結び付く。この計算で指数の係数は $\kappa=1/(4m)$ と直接分かる。殻と中心の扱いを滑らかにする場合にも、正則な零座標の間の一次の微分が非零なら先頭の指数は保たれる。

### 波束の粒子数を遅延時刻ごとに読む

正の出射周波数を $\omega$ とし、有限幅の区間 $I_j=[j\varepsilon,(j+1)\varepsilon]$ で波束を作る。$\varepsilon>0$ とし、
$$
p_{jn}=\varepsilon^{-1/2}\int_{I_j} e^{2\pi i n\omega/\varepsilon}p_\omega\,\dd\omega .
$$
$p_\omega\sim e^{-i\omega u}$ なら波束の中心は $u=2\pi n/\varepsilon$ 付近である。元の KG 内積に $\delta(\omega-\omega')$ を入れると同じ $j$ では $\varepsilon^{-1}\int_{I_j}e^{2\pi i(n'-n)\omega/\varepsilon}\dd\omega=\delta_{nn'}$ になり、異なる $j$ では零である。従って個々の波束の粒子数には規格化に由来する $\delta(0)$ が入らない。以下では区間が $\omega=0$ から離れた場合を考える。

散乱の過渡期を過ぎた定常近似では、周波数の粒子数相関は $\langle b^\dagger_\omega b_{\omega'}\rangle=n_\omega\delta(\omega-\omega')$ である。この形は遅延時刻のずれで二点の周波数成分の位相が同時に変わることから従う。そのとき波束の粒子数は $\varepsilon^{-1}\int_{I_j}n_\omega\dd\omega$ であり、狭い区間では代表値に近づく。一方、連続な時刻の全ての波束について足すと、一定の放射率が無限時間続く理想化では総数は発散する。従って平面波の規格化に出る無限大を、有限時間に観測した粒子数と解釈してはならない。遅い波束の近似を使うときは、幾何光学の誤差、初期状態、計量の時間変化の時間尺度を指定する。

### 動径散乱の流束で灰色体因子を定義する

実のポテンシャル $V(r_*)$ に対する $\psi''+(\omega^2-V)\psi=0$ を考える。$\mathcal J=(\overline\psi\psi'-\psi\overline\psi')/(2i)$ を微分すると、方程式と共役を代入して $\mathcal J'=0$ となる。無限遠で $\psi=e^{-i\omega r_*}+\mathcal R e^{i\omega r_*}$、未来地平線で $\psi=\mathcal T e^{-i\omega r_*}$ とする。前者で $\mathcal J=\omega(|\mathcal R|^2-1)$、後者で $\mathcal J=-\omega|\mathcal T|^2$ である。干渉項は共役との差で消える。従って $|\mathcal T|^2=1-|\mathcal R|^2$、これがこの静的な模型での透過率である。地平線への吸収を計算して放射の透過率に使えるのは、実ポテンシャルによる左右の透過率が等しいためである。実際、左入射解と右入射解の複素共役を取らない Wronskian は一定で、両端で評価するとそれぞれの透過振幅の等しさを与える。

:::example[障壁で出射スペクトルが変わる周波数]
Schwarzschild の $\ell=0$ 無質量スカラーでは $V=2m(1-2m/r)/r^3$ である。$V'=2m(8m-3r)/r^5$ から極大は $r=8m/3$ で、その値は $27/(1024m^2)$ である。従って角運動量が零でも $m\omega\sim\sqrt{27}/32$ の尺度が散乱に入る。波の透過は障壁の高さだけでなく幅にも依存するので、この値で透過率が鋭く切り替わるとは仮定しない。
:::

### Bogoliubov 比からエネルギー流束へ

外部に散乱障壁のない二次元で、一つの出射する無質量スカラーだけをとる。正の周波数 $\omega$ のモードの熱占有は $n_\omega=(e^{\beta\omega}-1)^{-1}$ である。周期箱の長さを $L$ とすると右向きのモードの間隔は $2\pi/L$ である。単位長さ当たりのエネルギーは $L^{-1}\sum_{\omega>0}\omega n_\omega$、伝播速度が $1$ なので同じ値がエネルギー流になる。連続極限で
$$
F=\frac{1}{2\pi}\int_0^\infty\frac{\omega\,\dd\omega}{e^{\beta\omega}-1}
=\frac{1}{2\pi\beta^2}\sum_{n\ge1}\int_0^\infty x e^{-nx}\dd x
=\frac{\pi}{12\beta^2}.
$$
級数は正項なので積分と和を交換でき、最後は $\sum n^{-2}=\pi^2/6$ を用いた。$\beta=2\pi/\kappa$ に対して $F=\kappa^2/(48\pi)$ である。四次元で使う場合は角モードの個数と灰色体因子が入るため、この係数をそのまま四次元の出力へ移さない。

:::example[粒子数流とエネルギー流の赤外極限]
同じ二次元の連続モードで粒子数流だけを積分すると、低周波で $n_\omega\sim1/(\beta\omega)$ のため対数発散する。一方、エネルギー流の被積分関数には $\omega$ が掛かり、低周波で定数 $1/\beta$ へ近づく。従って低周波のカットオフを $\omega_0$ へ入れた場合、取り除いたエネルギー流は先頭で $\omega_0/(2\pi\beta)$ にすぎない。総エネルギーが有限であることだけから、連続基底で数えた総粒子数も有限とは結論できない。
:::

### 永遠の外部領域で異なる量子状態を選ぶ

定常な Schwarzschild 外部では、無限遠から来る入射モードと、地平線側から来る出射モードを独立な散乱データに取れる。両方を Killing 時刻の正周波数で定めた空の状態を :term[Boulware 真空]{en="Boulware vacuum" yomi="ぶーるうぇあしんくう"} という。遠方では静止観測者の真空だが、非退化地平線を横切る正則な状態にはならない。出射部分の $u$ 真空と、横切る座標 $U=-C e^{-\kappa u}$ での真空が一致しないことは、先の正規化比からも分かる。

将来の地平線を通る崩壊の遅い外部に対応させるには、過去無限遠からの入射モードを空に保ち、出射部分には地平線の正則な零座標の正周波数を使う。永遠の外部上でこの選択をした状態は :term[Unruh 真空]{en="Unruh vacuum" yomi="うんるーしんくう"} である。遠方に熱因子を持つ出射がある一方、過去無限遠から補う熱浴はない。これに対し、両方の地平線で正則になるよう出射と入射に平衡の熱占有を置くのが :term[Hartle–Hawking 状態]{en="Hartle–Hawking state" yomi="はーとるほーきんぐじょうたい"} である。時間を反転すると入射と出射を交換するので、平衡では無限遠での正味のエネルギー流は零である。これらは永遠の静的ブラックホールにおける区別である。崩壊初期から全時刻の状態が Unruh 真空で記述されると仮定するものではない。正則性の分布論的な証明と大域的な状態の構成は :cite[wald-gr] を参照する。

:::example[二次元で状態ごとの入射・出射の流束を比較する]
散乱のない二次元の外部で、無限遠の慣性座標 $(u,v)$ の真空を基準に流束を測る。Boulware 型では両方向で零、Unruh 型では出射が $\kappa^2/(48\pi)$、入射が零である。Hartle–Hawking 型では両方向の流束の大きさが $\kappa^2/(48\pi)$ で、正味の動径流は零である。従って二方向とも空の状態と、両方向に熱占有がある状態は、正味の流束だけでは区別できない。二点相関やエネルギー密度を比較すれば区別できる。この例では無限遠での流束を比べているため、これらの数値を地平線上の座標成分へそのまま代入しない。
:::

### 回転があると周波数の閾値が変わる

軸対称な定常背景で、モードの時間と方位角の依存を $e^{-i\omega t+i m_\phi\phi}$ とする。地平線の生成子は $\partial_t+\Omega_H\partial_\phi$ なので、地平線に沿う位相の微分で $\tilde\omega=\omega-m_\phi\Omega_H$ が出る。吸収を調べる動径式を実係数の散乱形に直し、両端で波が平面波になる場合を考える。単位の入射を無限遠に置けば、動径 Wronskian 流は無限遠で $\omega(|\mathcal R|^2-1)$、地平線で $-\tilde\omega|\mathcal T|^2$ である。従って
$$
1-|\mathcal R|^2=\frac{\tilde\omega}{\omega}|\mathcal T|^2 .
$$
$0<\omega<m_\phi\Omega_H$ では右辺が負で、反射が入射より強い。これが古典波の :term[超放射]{en="superradiance" yomi="ちょうほうしゃ"} である。従って回転系の灰色体係数を吸収率として定義すると、この周波数域では負になり得る。Hawking 放射のスカラーの熱因子に入るのも $\tilde\omega$ であり、占有の式は吸収率を $e^{2\pi\tilde\omega/\kappa}-1$ で割る。超放射域では分子と分母が両方負なので、出射占有は正である。非回転のように各モードの吸収率を $|\mathcal T|^2$ と同一視してはいけない。ここでは非退化の定常背景で散乱と地平線の近傍計算を分けた。回転外部の全体で静的な熱浴を作れるという結論は従わない。

## 膨張宇宙での粒子生成

平坦 Robertson–Walker 計量を $a(\eta)^2(-\dd\eta^2+\dd\vb x^2)$ と書き、$\phi=a^{-1}\chi$ と置く。$\Box_g\phi=a^{-2}(-\phi''-2a'\phi'/a+\Delta\phi)$、$R=6a''/a^3$ である。Fourier モードは

$$
\chi_{\vb k}''+\left[k^2+a^2m_s^2+(6\xi-1)\frac{a''}a\right]\chi_{\vb k}=0 .
\label{eq-gr-flrw-modes}
$$

従って質量ゼロかつ $\xi=1/6$ では一定周波数 $k$ となり、共形時刻で定義した正周波数は混合しない。これは共形真空の粒子数についての結論である。一般には時間依存する周波数 $\Omega$ の振動子となる。$w=(2\Omega)^{-1/2}e^{-i\int\Omega\dd\eta}$ を代入すると残差は $(-\Omega''/(2\Omega)+3\Omega'^2/(4\Omega^2))w$ であり、$\abs{\Omega'}/\Omega^2$ や $\abs{\Omega''}/\Omega^3$ が小さいと断熱近似が有効になる。

比較のため振動子の周波数が $\eta=0$ で $\omega_i$ から $\omega_f$ に急変するとする。ここでは初期の正周波数解を後期の正・負周波数解で展開した係数を $\alpha,\beta$ とする。$w,w'$ の連続性から $\alpha+\beta=\sqrt{\omega_f/\omega_i}$、$\alpha-\beta=\sqrt{\omega_i/\omega_f}$ なので $\abs\beta^2=(\omega_f-\omega_i)^2/(4\omega_f\omega_i)$ である。この向きの変換を逆にすると負周波数係数の符号が反転するが、占有数は同じである。例えば $\omega_f=2\omega_i$ なら $1/8$ である。この模型は有界な周波数の急変を理想化したもので、計量の微分にデルタ関数がある模型では一致条件から求め直す必要がある。

### 断熱振幅の連立方程式

$\chi''+\Omega^2\chi=0$ で $\Omega>0$ とし、$\Theta'=\Omega$ とする。二つの係数 $c_\pm$ を
$$
\begin{aligned}
\chi &= \frac{c_+e^{-i\Theta}+c_-e^{i\Theta}}{\sqrt{2\Omega}},\\
\chi'&=-i\sqrt{\frac{\Omega}{2}}(c_+e^{-i\Theta}-c_-e^{i\Theta})
\end{aligned}
$$
で定める。一行目を微分して二行目から引くと $c_+'e^{-i\Theta}+c_-'e^{i\Theta}=(\Omega'/2\Omega)(c_+e^{-i\Theta}+c_-e^{i\Theta})$ となる。二行目を微分して波動方程式を使うと、差の式 $c_+'e^{-i\Theta}-c_-'e^{i\Theta}=-(\Omega'/2\Omega)(c_+e^{-i\Theta}-c_-e^{i\Theta})$ を得る。足し引きして
$$
c_+'=\frac{\Omega'}{2\Omega}e^{2i\Theta}c_-,
\qquad c_-'=\frac{\Omega'}{2\Omega}e^{-2i\Theta}c_+
$$
となる。従って $|c_+|^2-|c_-|^2$ は一定で、初期に $c_+=1,c_-=0$ なら $1$ である。混合が小さいとき $c_+\simeq1$ を代入した $c_-$ の積分で粒子生成の先頭を得る。$\Omega$ の変化が緩やかなら位相の速い振動が寄与を相殺する。$\Omega$ が零に近づくとこの展開の小パラメータが失われるので、元の方程式をそこで解く必要がある。

:::example[周波数の急変で規格化・平均・真空確率を照合する]
$s=\omega_f/\omega_i>0$ の急変モデルで、$\chi,\chi'$ の一致条件を使う。実の $\alpha=(\sqrt s+1/\sqrt s)/2,\beta=(\sqrt s-1/\sqrt s)/2$ となり、$\alpha^2-\beta^2=1$ が保たれる。単一実振動子で測る平均粒子数は $\beta^2$、真空の重なりの確率は $1/|\alpha|$ である。

ここでは$s$ は後と前の角周波数の比であるとする。計算に用いる式を先に並べる。

$$
\begin{aligned}
\alpha&=\frac{\sqrt{s}}{2} + \frac{1}{2 \sqrt{s}}\\
\beta&=\frac{\sqrt{s}}{2} - \frac{1}{2 \sqrt{s}}\\
\langle N\rangle&=\frac{\left(s - 1\right)^{2}}{4 s}\\
P_0&=\frac{2 \sqrt{s}}{s + 1}
\end{aligned}
$$

次の条件を比較する。入力を区別してから計算し、定数の規格化をそろえた値を表に示す。

| $s$ | $\alpha$ | $\beta$ | $\langle N\rangle$ | $P_0$ |
|---|---|---|---|---|
| $\frac{1}{4}$ | $\frac{5}{4}$ | $- \frac{3}{4}$ | $\frac{9}{16}$ | $\frac{4}{5}$ |
| $\frac{1}{2}$ | $\frac{3 \sqrt{2}}{4}$ | $- \frac{\sqrt{2}}{4}$ | $\frac{1}{8}$ | $\frac{2 \sqrt{2}}{3}$ |
| $1$ | $1$ | $0$ | $0$ | $1$ |
| $2$ | $\frac{3 \sqrt{2}}{4}$ | $\frac{\sqrt{2}}{4}$ | $\frac{1}{8}$ | $\frac{2 \sqrt{2}}{3}$ |
| $4$ | $\frac{5}{4}$ | $\frac{3}{4}$ | $\frac{9}{16}$ | $\frac{4}{5}$ |
| $9$ | $\frac{5}{3}$ | $\frac{4}{3}$ | $\frac{16}{9}$ | $\frac{3}{5}$ |

$s$ を逆数に変えると $\beta$ は符号を変えるが、平均数と真空の確率は変わらない。急変がない $s=1$ では $\beta=0$ で、旧真空と新真空が一致する。対になる複素 Fourier モードを数える場合は、どの実自由度を一振動子とみなすかを指定して重複計数を避ける。
:::

### 共形時刻での二つの厳密解

平坦な Robertson–Walker で $a(\eta)=A\eta^p$、場は質量ゼロ、$\eta>0$ とする。:eqref[eq-gr-flrw-modes] の周波数二乗は $k^2+(6\xi-1)p(p-1)/\eta^2$ である。放射の模型で $p=1$ なら、最小結合でもこの追加項が零で、$\chi=e^{-ik\eta}/\sqrt{2k}$ を同じ規格化で使える。物質の模型で $p=2,\xi=0$ なら $\chi''+(k^2-2/\eta^2)\chi=0$ となる。$\chi=e^{-ik\eta}F(\eta)$ と置けば $F''-2ik F'-2F/\eta^2=0$ で、$F=1+C/\eta$ を入れると $2ikC-2=0$ である。従って
$$
\chi=\frac{1-i/(k\eta)}{\sqrt{2k}}e^{-ik\eta}.
$$
Wronskian $i(\chi^*\chi'-\chi\chi'^*)$ は動径方程式で一定であり、$\eta\to\infty$ の平面波で $1$ なので全域で $1$ である。この選択は遅い時刻での平面波に合う。一方、$\eta\to0^+$ では $\Omega^2=k^2-2/\eta^2$ が負であり、その時刻で正の実周波数を用いた初期真空はこの式からは選べない。

### 滑らかな膨張に対する最初の摂動

共形結合 $\xi=1/6$ とし、$a^2=a_i^2+\delta a^2\,F(\eta)$、$F$ が過去で $0$、未来で $1$ へ滑らかに近づくとする。$m_s>0$ のとき $\omega_i^2=k^2+m_s^2a_i^2$、$\Omega^2=\omega_i^2+m_s^2\delta a^2 F$ である。$|m_s^2\delta a^2|\ll\omega_i^2$ で、$F'$ が有限時間の近くに集中する場合、断熱振幅の式を $\delta a^2$ の一次まで展開して
$$
c_-^{(1)}(+\infty)=\frac{m_s^2\delta a^2}{4\omega_i^2}
\int_{-\infty}^{\infty}F'(\eta)e^{-2i\omega_i\eta}\dd\eta
$$
を得る。係数は $\Omega'/(2\Omega)=(\Omega^2)'/(4\Omega^2)$ から来る。位相の補正は、それ自体が一次で微分係数も一次なので、変化のある有限の時間尺度では二次に入る。無限に続く定常部分の位相は粒子数の絶対値に影響しない。

例えば $F'(\eta)=e^{-\eta^2/T^2}/(\sqrt\pi T)$ と取る。Gauss 積分の Fourier 変換で $c_-^{(1)}=(m_s^2\delta a^2/(4\omega_i^2))e^{-\omega_i^2T^2}$ である。$T\to0$ では急変の式を $\omega_f^2-\omega_i^2=m_s^2\delta a^2$ の一次まで展開した係数に一致する。固定の $T>0$ では大きな $\omega_i$ の振幅が抑えられる。一定の小さい摂動で $T$ を非常に大きく取ると位相の補正が長時間に蓄積するため、この一次式の絶対誤差と指数の精度は同じではない。滑らかで遅い変化の一般の抑制は、振動積分や高次の断熱近似で評価する。

### 時間変化があっても粒子を残さない厳密例

振動子の周波数を $\Omega^2(\eta)=\omega^2+2\mu^2\operatorname{sech}^2(\mu\eta)$、$\omega,\mu>0$ とする。例えば $\xi=1/6$ の質量を持つ場で、漸近値 $a_0^2$ に $2\mu^2\operatorname{sech}^2(\mu\eta)/m_s^2$ を足した共形計量がこの形を与える。この背景は一度膨らんで戻り、過去と未来の $\omega=\sqrt{k^2+m_s^2a_0^2}$ は同じである。$W(\eta)=\mu\tanh(\mu\eta)$ とすると $W'=\mu^2-W^2$ であり、
$$
(\partial_\eta+W)(\partial_\eta-W)=\partial_\eta^2-\mu^2,\qquad
(\partial_\eta-W)(\partial_\eta+W)=\partial_\eta^2-\mu^2+2W'
$$
である。従って初期平面波に $\partial_\eta-W$ を作用させた
$$
\chi(\eta)=\frac{-i\omega-W(\eta)}{-i\omega+\mu}
\frac{e^{-i\omega\eta}}{\sqrt{2\omega}}
$$
は振動子方程式を満たし、過去で正しく規格化される。未来では係数が $(-i\omega-\mu)/(-i\omega+\mu)$ になり、絶対値は $1$ である。負周波数は現れない。この例では周波数は途中で変わり、単純な断熱条件の成否にかかわらず、混合が全体で打ち消される。「時間依存の周波数なら必ず粒子が生じる」という命題は成り立たない。

### 空間の並進対称性と対の生成

実スカラーの $\vb k$ と $-\vb k$ のモードでは、複素共役を取ると空間位相が反転する。時間だけが変わる背景では空間の全運動量は保存するので、波数 $\vb k$ の出射消滅演算子が生成演算子と混ざるとき、相手は $-\vb k$ になる。有限周期箱で両者が異なるとき、例えば $b_{\vb k}=\alpha^*a_{\vb k}-\beta^*a_{-\vb k}^\dagger$ とその反転の組に対して、真空は前節の二モードの形である。各側の平均占有数は $|\beta|^2$、全運動量は各項で零である。各側の粒子数を足すと $2|\beta|^2$ だが、同じ二モードの組を数える回数は一回である。例えば $\omega_f=2\omega_i$ の急変なら各側が $1/8$、組全体で $1/4$ である。エネルギーは時間並進がないため保存を要求できず、後期の定常領域での真空差のエネルギーは、この組で $2\omega_f|\beta|^2$ になる。

## 発散の処理とエネルギー運動量テンソルのくりこみの概観

点 $\widehat\phi(x)^2$ や $\widehat T_{ab}(x)$ は分布の同一点積を含む。例えば平坦な質量ゼロ真空の $\langle\widehat\phi^2\rangle$ は運動量の上限を $K$ とすると $\int^{K}\dd^3k/[(2\pi)^3\,2k]=K^2/(8\pi^2)$ である。エネルギー密度では一振動子当たり $k/2$ が入り、$K^4/(16\pi^2)$ になる。このカットオフは特定の慣性枠を使い、一般の時空でそのまま共変な応力を与えない。

:term[点分離くりこみ]{en="point-splitting renormalization" yomi="てんぶんりくりこみ"} では二点関数 $W(x,x')$ から、局所幾何で決まる Hadamard 特異部 $H(x,x')$ を引く。その差に応力を定める微分演算子を作用させて一致極限を取る。:term[Hadamard 状態]{en="Hadamard state" yomi="あだまーるじょうたい"} は、局所的にこの差が滑らかになる短距離構造を持つ状態である。$\sigma$ を二点間の符号付き測地距離二乗の半分とすると、四次元の特異部の形は $(8\pi^2)^{-1}[U/\sigma_++V\log(\mu^2\sigma_+)]$ である。光円錐側の境界値を $\sigma_+$ として指定し、$U,V$ は波動演算子を作用させて生じる各特異次数の係数が消えるよう運搬方程式で決める。二つの Hadamard 状態の応力差は同じ $H$ が相殺するので、二点関数の滑らかな差から求められる。絶対値には保存性を保つ局所補正が必要で、宇宙項、$R$、曲率二乗作用の変分に対応する有限な自由度が残る。局所的な減算を共変な演算子値分布として構成する全般の証明には分布の特異性の伝播を制御する必要があるため、この章では具体的な最低次数と状態差を計算する。

### 状態差なら局所項を消せる

同じ背景に二状態を置き $W_1-W_2$ が滑らかであるとする。最小結合の実スカラーでは古典的応力は $\partial_a\phi\partial_b\phi-\tfrac12g_{ab}[(\nabla\phi)^2+m_s^2\phi^2]$ である。そこで二点関数の滑らかな差に、添字を一致点へ平行移動してから対称化した $\nabla_a\nabla_{b'}$ とその跡、さらに質量項を作用させ、一致極限を取る。状態と独立な Hadamard 減算と局所項は双方から引くと消えるため、絶対値の規約を指定しなくてもこの差を得られる。例えば有限周期箱の真空と単一の正規化モードに $n$ 個の量子を置いた数状態の差は $n[u(x)u^*(x')+u^*(x)u(x')]$ である。$u=e^{i\vb k\cdot\vb x-i\omega t}/\sqrt{2\omega V}$ に空間・時間微分をかけるとエネルギー密度の差は $n(\omega^2+k^2+m_s^2)/(2\omega V)=n\omega/V$ となる。従って空間積分で量子 $n$ 個のエネルギーが得られる。

### 場の曲率結合が応力にも現れる

作用の曲率項 $-\xi\int\sqrt{-g}R\phi^2/2$ を逆計量で変分する。$\delta(\sqrt{-g}R)$ のバルク項は $\sqrt{-g}G_{ab}\delta g^{ab}$ である。一方、通常は発散となる $g^{ab}\delta R_{ab}$ に $\phi^2$ を掛けて二度部分積分すると
$$
\delta S_{\rm curvature}
=-\frac{\xi}{2}\int\sqrt{-g}
[\phi^2 G_{ab}+g_{ab}\Box(\phi^2)-\nabla_a\nabla_b(\phi^2)]\delta g^{ab}.
$$
従って角括弧に $\xi$ を掛けた量が局所応力に加わる。平坦背景でも $\phi^2$ の微分の項は一般に残る。四次元で $m_s=0,\xi=1/6$ とすると、古典応力の跡は $-(\nabla\phi)^2+[-R\phi^2+3\Box\phi^2]/6$ である。$\Box\phi^2=2(\nabla\phi)^2+2\phi\Box\phi$ と運動方程式 $\Box\phi=R\phi/6$ を入れると零になる。量子論でくりこんだ応力では、差し引く局所特異項にも微分が作用するため、この古典の代入で跡が零になるとは限らない。

### 世界関数から先頭の Hadamard 係数を求める

凸正規近傍で世界関数 $\sigma$ を測地距離の符号付き二乗の半分とする。放射測地線の Gauss の補題から $\nabla^a\sigma\nabla_a\sigma=2\sigma$ である。$P=\Box_g-m_s^2-\xi R$ とし、点が一致せず零面から外れた部分で積の微分を行うと
$$
P(U/\sigma)=\frac{PU}{\sigma}
+\frac{(4-\Box\sigma)U-2(\nabla^a\sigma)\nabla_a U}{\sigma^2}.
$$
実際、$\nabla_a(\sigma^{-1})=-\sigma^{-2}\nabla_a\sigma$、$\Box(\sigma^{-1})=(4-\Box\sigma)/\sigma^2$ である。従って最も強い特異項を消す運搬方程式は
$$
2\nabla^a\sigma\,\nabla_aU+(\Box\sigma-4)U=0,\qquad [U]=1
$$
である。$[U]$ は一致極限を表し、局所 Minkowski 規格化で $1$ とする。平坦な四次元では $\Box\sigma=4$ で、$\nabla\sigma$ は基準点からの放射方向なので $U=1$ となる。次の次数では $PU/\sigma$ と $P(V\log\sigma)$ の $1/\sigma$ 部分が打ち消し合う条件がある。$P(V\log\sigma)=(PV)\log\sigma+[2\nabla\sigma\cdot\nabla V+(\Box\sigma-2)V]/\sigma$ なので、平坦時空の一致点では $-m_s^2+2[V]=0$、従って $[V]=m_s^2/2$ である。これは境界値にして分布を拡張する前に行える局所計算であり、点分離処方で保存される応力全体の構成には追加の局所項が必要である。

### 曲率が対数発散の係数に入る理由

基準点 $x'$ からの正規座標を $y^a$ とする。Gauss の補題により $\nabla^a\sigma=y^a$ である。正規座標の二次の計量を $g_{ab}=\eta_{ab}-R_{acbd}(x')y^cy^d/3+\order{|y|^3}$ とすると、行列式の一次変化は跡なので $\sqrt{-g}=1-R_{cd}y^cy^d/6+\order{|y|^3}$ である。従って
$$
\Box\sigma=\nabla_a y^a
=4+y^a\partial_a\log\sqrt{-g}
=4-\tfrac13R_{cd}y^cy^d+\order{|y|^3}.
$$
$U$ の運搬式にこの結果を代入する。一次係数は零であり、二次項は放射方向の微分で二倍されるので $U=1+R_{ab}y^ay^b/12+\order{|y|^3}$ である。従って $[\Box U]=R/6$ となる。$1/\sigma$ の運搬式を一致点で取ると $[PU]+2[V]=0$ だから、
$$
[V]=\tfrac12\bigl(m_s^2+(\xi-1/6)R\bigr)
$$
である。これは膨張模型で紫外の対数項を作る組合せと一致する。実際、共形時刻の $R=6a''/a^3$ により $\mathcal V=a^2(\xi-1/6)R$ である。質量ゼロで共形結合ならこの係数は消える。ただし $U$ の曲率補正と $V$ の高次係数は別に残り得るので、局所応力の跡まで零とは結論しない。

### 点分離で定数の差が有限に残る例

:::example[曲率項に比例する有限な変更]
局所作用に $\delta S=C\int\sqrt{-g}R$ を足すと、コンパクトな台の逆計量変分の係数は $CG_{ab}$ である。従ってこれを物質側の応力の変更と数える規約では $\Delta T_{ab}=-2CG_{ab}$ で、縮約した Bianchi 恒等式により保存する。Minkowski では $G_{ab}=0$ で寄与がなく、四次元の de Sitter で $R_{ab}=3H^2g_{ab}$ なら $G_{ab}=-3H^2g_{ab}$ なので $\Delta T_{ab}=6CH^2g_{ab}$ となる。この差は背景が変われば変化するが、同じ背景での二状態の応力の差では両方に同じ量が足されて相殺する。
:::

### 膨張模型で断熱減算の最初の係数を導く

:eqref[eq-gr-flrw-modes] の周波数を $\omega^2+\mathcal V$、$\omega^2=k^2+a^2m_s^2$、$\mathcal V=(6\xi-1)a''/a$ と書く。$\mathcal V$ は尺度因子の微分を二回含むと数える。モードを $\chi=(2W_\mathrm{ad})^{-1/2}e^{-i\int W_\mathrm{ad}\dd\eta}$ とすれば、方程式を満たすための条件は
$$
W_\mathrm{ad}^2=\omega^2+\mathcal V
-\frac{W_\mathrm{ad}''}{2W_\mathrm{ad}}+\frac{3(W_\mathrm{ad}')^2}{4W_\mathrm{ad}^2}
$$
である。振幅の一階微分と位相の二階微分による虚部が相殺するため、実の式だけが残る。$W_\mathrm{ad}=\omega+\delta W$ とし二階まで残すと
$$
\delta W=\frac{\mathcal V}{2\omega}
-\frac{\omega''}{4\omega^2}+\frac{3(\omega')^2}{8\omega^3}.
$$
従って $|\chi|^2=(2W_\mathrm{ad})^{-1}$ の二階までの減算項は $1/(2\omega)-\mathcal V/(4\omega^3)+\omega''/(8\omega^4)-3(\omega')^2/(16\omega^5)$ である。$k\to\infty$ では $1/(2\omega)=1/(2k)-a^2m_s^2/(4k^3)+\cdots$、また $\omega',\omega''=\order{k^{-1}}$ だから、$|\chi|^2$ の二次と対数の紫外発散を作るのは $1/(2k)$ と $-(a^2m_s^2+\mathcal V)/(4k^3)$ である。$\langle\phi^2\rangle$ には $a^{-2}\int\dd^3 k/(2\pi)^3$ が掛かるので、この二項は運動量測度を入れて $k$ と $1/k$ に比例する。高波数でこの展開に従う状態に対して、これらを減算すれば紫外側の発散を除ける。

応力では時間・空間微分が加わり、運動量についてさらに二つの次数が必要になる。四階まで再帰的に $W_\mathrm{ad}$ を求めて応力全体を同じ規約で減算する。二階までの展開は局所的な高波数の式であり、$k=0$ までこの漸近展開が一様に正しいとは限らない。特に $m_s=0$ で $\mathcal V\ne0$ の場合、単に $1/k^3$ を全域で引くと逆に赤外側の対数発散を入れる。減算の局所性を保つため、質量尺度と有限局所項を含めた規約を指定する。例えば $a$ が定数で Minkowski の真空を選ぶと $\mathcal V=\omega'=\omega''=0$ であり、正確な $|\chi|^2=1/(2\omega)$ と零階の減算項が一致する。

### 点分離から二次元の放射流を再計算する

散乱のない二次元の無質量スカラーの出射部分を考える。漸近的に平坦な領域の時刻を $u$ とし、状態を選ぶための零時刻を $p=p(u)$、$\dd p/\dd u>0$ とする。出射微分場の真空二点関数は、二点の時刻差を $\Delta p$ と書けば $\int_0^\infty \omega e^{-i\omega(\Delta p-i0)}\dd\omega/(4\pi)=-[4\pi(\Delta p-i0)^2]^{-1}$ である。場の微分を $u$ に関して取れば、分子に各点の $\dd p/\dd u$ が掛かる。同じ平坦な漸近領域で $u$ 真空の二点関数を引けば、局所減算は共通なので状態差として流束を読める。

点を $u\pm s/2$ とし、$p_\pm=p(u\pm s/2)$、$\dot p=\dd p/\dd u$ と書く。Taylor 展開により
$$
\begin{aligned}
p_+-p_- &= s\dot p+s^3\dddot p/24+\order{s^5},\\
\dot p_+\dot p_- &= \dot p^{\,2}
+\tfrac14 s^2(\dot p\,\dddot p-\ddot p^{\,2})+\order{s^4}.
\end{aligned}
$$
従って比 $\dot p_+\dot p_-/(p_+-p_-)^2$ と $1/s^2$ の差は
$$
\frac16\left[\frac{\dddot p}{\dot p}-\frac32\left(\frac{\ddot p}{\dot p}\right)^2\right]
$$
に収束する。角括弧は :term[Schwarz 微分]{en="Schwarzian derivative" yomi="しゅわるつびぶん"} である。二点関数の係数を掛けると、出射のくりこまれた $\Delta T_{uu}$ はこの微分の $-1/(24\pi)$ 倍である。$p=p_0-C e^{-\kappa u}$ では $\ddot p/\dot p=-\kappa$、$\dddot p/\dot p=\kappa^2$ となり、$\Delta T_{uu}=\kappa^2/(48\pi)$ である。これは熱占有のエネルギー積分と一致する。例えば $p$ が $u$ のアフィン変換なら $\ddot p=0$ なので流束差はない。この導出は微分場を用いているため、二次元の無質量場の定数モードの曖昧さが流束に入らない。

### 曲率二乗作用で跡の規約を動かす

$b\int\sqrt{-g}R^2\dd^4x$ を有限な局所項として選ぶ。その逆計量変分では体積の変分から $-g_{ab}R^2/2$、$2R\delta R$ の Ricci 項から $2RR_{ab}$、二回部分積分する項から $2(g_{ab}\Box R-\nabla_a\nabla_bR)$ が出る。従って応力側にこの項を置く規約では
$$
\Delta T_{ab}=-2b\left[2RR_{ab}-\tfrac12g_{ab}R^2
+2(g_{ab}\Box R-\nabla_a\nabla_bR)\right].
$$
括弧の発散は零である。実際、最初の二項の発散は縮約 Bianchi 恒等式を使うと $2R_{ab}\nabla^aR$ であり、最後の項から交換子によりその負が出る。四次元の跡は $\Delta T^a{}_a=-12b\Box R$ である。このため共形場の :term[跡異常]{en="trace anomaly" yomi="せきいじょう"} のうち $\Box R$ の係数はこの有限項で変えられる。古典の共形方程式だけから量子の跡が零とは言えず、さらに跡のうちどの局所項が規約に依存するかを分ける必要がある。

例えば平坦な Robertson–Walker で宇宙時刻の尺度因子 $a=t^p$、$t>0$ とする。$R=6p(2p-1)/t^2$、$\Box R=-\ddot R-3H\dot R$ である。従って $\Box R=36p(2p-1)(p-1)/t^4$ である。$p=1/2$ では $R=0$ なのでこの変更は零だが、時空の Ricci テンソルや Riemann テンソルまで零とは限らない。例えば同じ放射模型で $R_{00}=3/(4t^2)$ である。従ってこの一項が消えることから、くりこんだ応力や他の曲率による寄与が消えるとは結論しない。

## 演習問題

::::exercise{level=basic}
Minkowski モード $e^{-i\omega t+i\vb k\cdot\vb x}$ のオンシェル条件を求めよ。

:::solution
$(\Box-m_s^2)$ を作用させると $\omega^2-k^2-m_s^2$ が掛かるので、$\omega^2=k^2+m_s^2$ が必要である。
:::
::::

::::exercise{level=basic}
一モードで $\beta=3/4$ の新粒子数を旧真空で求めよ。

:::solution
:eqref[eq-gr-bogoliubov-number] により期待数は $\abs\beta^2=9/16$ である。
:::
::::

::::exercise{level=basic}
周波数変化 $\omega_f=3\omega_i$ の振動子の粒子数を求めよ。

:::solution
$\abs\beta^2=(3-1)^2/(4\cdot3)=1/3$ となる。
:::
::::

::::exercise{level=basic}
Unruh 検出器で $a$ が二倍なら温度は何倍か。

:::solution
$T_U=\hbar a/(2\pi k_B)$ なので二倍である。
:::
::::

::::exercise{level=advanced}
同じ状態の二つの異なる Cauchy 面で KG 内積が一致することを示せ。

:::solution
$j^a$ の発散が零なので、面間の体積積分に発散定理を用いると未来面の流束から過去面の流束を引いたものが零である。空間無限遠に境界流束がないと仮定すれば、両面で :eqref[eq-gr-kg-product] が一致する。
:::
::::

::::exercise{level=advanced}
$\alpha=\cosh s,\beta=\sinh s$ の一モードで $s=\log 2$ の粒子数を求めよ。

:::solution
$\alpha^2-\beta^2=1$ で正準規格化される。$\sinh(\log2)=(2-1/2)/2=3/4$ なので粒子数は $9/16$ である。
:::
::::

::::exercise{level=advanced}
$a(\eta)=A\eta^2$、$\eta>0$ の質量ゼロ最小結合でモード周波数の二乗を求めよ。

:::solution
$\xi=0$ なので :eqref[eq-gr-flrw-modes] の補正は $-a''/a=-2/\eta^2$ である。従って $\Omega^2=k^2-2/\eta^2$ となり、$k\eta\gg1$ で平坦な分散関係に近づく。
:::
::::

::::exercise{level=advanced}
Schwarzschild の Hawking 放射で $\Gamma_\ell(\omega)=1/5$、$e^{2\pi\omega/\kappa}=3$ のとき無限遠の当該モードの占有数を求めよ。

:::solution
熱分布の占有数は $(3-1)^{-1}=1/2$ であり、透過で $1/5$ 倍される。従って $1/10$ になる。
:::
::::

::::exercise{level=advanced}
周波数が $\omega$ から $s^2\omega$ へ瞬時に変わった後、時間 $\Delta$ 経過して元の $\omega$ に瞬時に戻る実振動子を考える。$s>0$ とする。最終粒子数を $\Delta$ の関数で求めよ。

:::solution
最初の跳びの係数は $\alpha=(s+s^{-1})/2$、$\beta=(s-s^{-1})/2$ である。中間周波数 $\omega_2=s^2\omega$ での位相は $\theta=\omega_2\Delta$ であり、第二の跳びの直前の係数は $(\alpha e^{-i\theta},\beta e^{i\theta})$ になる。逆の跳びは同じ $\alpha$ と反対の $\beta$ を持つため、最終負周波数係数は $-\beta\alpha e^{-i\theta}+\alpha\beta e^{i\theta}=2i\alpha\beta\sin\theta$ である。従って
$$
N_{\rm final}=4\alpha^2\beta^2\sin^2(s^2\omega\Delta)
=\tfrac14(s^2-s^{-2})^2\sin^2(s^2\omega\Delta).
$$
最終周波数が初期と等しくても、中間で位相が進むので一般に粒子が残る。$\theta$ が $\pi$ の整数倍なら二度の跳びが相殺し、$\theta=\pi/2$ ならこの $s$ での最大値を得る。
:::
::::

## 確認テスト

::::quiz
Bogoliubov 変換で旧真空に新粒子を作る係数を選べ。

- [x] $\beta$
- [ ] $\alpha$ だけ
- [ ] どの周波数でも $\alpha+\beta$ のみ

:::explanation
新消滅演算子には旧生成演算子が $-\beta^*$ の係数で入る。旧真空では旧消滅演算子の項が消え、期待数に $\sum_j|\beta_{ij}|^2$ が残る。
:::
::::

::::quiz{answer="0.125" tol="1e-9"}
振動子の周波数が瞬時に $2\omega_i$ になった例の生成粒子数を求めよ。

:::explanation
モードとその微分を連続に接ぐと $\abs\beta^2=(2\omega_i-\omega_i)^2/(4(2\omega_i)\omega_i)=1/8$ である。
:::
::::

::::quiz
平坦な Robertson–Walker 計量で共形時刻のモードが混合しない条件をすべて選べ。ただしこの組合せで任意の尺度因子を許すものとする。

- [x] 質量が零
- [x] $\xi=1/6$
- [ ] $\xi=0$ が必要

:::explanation
:eqref[eq-gr-flrw-modes] の $a$ 依存項は質量項と $(6\xi-1)a''/a$ である。二条件を同時に課すと周波数が常に $k$ で、正周波数のモードの時間発展に負周波数は混ざらない。
:::
::::

::::quiz{answer="0.1" tol="1e-9"}
$e^{2\pi\omega/\kappa}=6$、透過確率 $1/2$ の Hawking モードの無限遠での平均占有数を求めよ。

:::explanation
地平線の熱因子による占有数は $1/(6-1)=1/5$ である。透過確率を掛けると $1/10=0.1$ になる。
:::
::::

::::quiz{answer="0.25" tol="1e-7"}
定常な Unruh 検出器で $\beta\omega=\log 4$ とする。励起と下向き率の比を求めよ。

:::explanation
詳細釣り合いの比は $e^{-\beta\omega}=e^{-\log4}=1/4$ である。温度を変えると励起率と下向き率は同じ倍率では変わらない。
:::
::::

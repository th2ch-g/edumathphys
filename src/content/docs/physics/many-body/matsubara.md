---
title: "有限温度 Green 関数と松原形式"
description: "熱平衡の虚時間相関を松原周波数で計算し、スペクトルと実時間応答へ変換する。"
notation:
  - symbol: '\beta'
    meaning: 逆温度
  - symbol: '\omega_n,\nu_m'
    meaning: Fermi および Bose の松原周波数
  - symbol: 'f,n_B'
    meaning: Fermi および Bose 分布
---

## この章の目標

- 熱平均から虚時間 Green 関数の境界条件を導く。
- 留数と部分分数を使って松原和を計算する。
- スペクトル表示を解析接続し、高周波の総和則を導く。
- 二次 Hamiltonian の Gibbs 状態で Wick の定理を証明する。

## 虚時間発展と温度 Green 関数

実時間の発展は状態に位相を付ける。一方、$t=-i\tau$ とすると位相因子が指数重みになるため、熱平衡の密度演算子と同じ指数にまとまる。:ref[def-stm-quantum-statistics] の熱平均を $\langle X\rangle_\beta=Z^{-1}\Tr(e^{-\beta K}X)$、$Z=\Tr e^{-\beta K}$ と書く。$K=H-\mu\hat N$、$\beta=1/T$ であり、温度を K で指定するときは $\beta=1/(k_BT)$ である。以下、まず有限体積でトレースが収束し、必要な演算子積の期待値も有限であるとする。

:::definition[温度 Green 関数]{#def-mb-matsubara-green}
$X(\tau)=e^{\tau K}Xe^{-\tau K}$ を虚時間発展とする。$T_\tau$ は虚時間が大きい演算子を左へ置き、Fermi 場どうしの交換には負号を掛ける。統計の符号を Bose では $\eta=1$、Fermi では $\eta=-1$ とし、:term[温度 Green 関数]{en="thermal Green function" yomi="おんどぐりーんかんすう"} を
$$
\mathcal G_{ij}(\tau)=-\langle T_\tau a_i(\tau)a_j^\dagger\rangle_\beta
$$
と定義する。$-\beta<\tau<0$ では $-\eta\langle a_j^\dagger a_i(\tau)\rangle_\beta$ である。
:::

$\partial_\tau X=[K,X]$ であり、虚時間発展では一般に $X(\tau)^\dagger=X^\dagger(-\tau)$ となる。実時間の $G^T$ と比べ、定義には $i$ が付かない。したがって跳びは $\mathcal G(0^+)-\mathcal G(0^-)=-I$、密度は $\gamma=-\eta\mathcal G(0^-)$ である。

:::example[一つの自由軌道]
$K=\xi a^\dagger a$ なら $\mathcal G(\tau)=-(1+\eta n)e^{-\xi\tau}$（$0<\tau<\beta$）、$n=\langle a^\dagger a\rangle_\beta$ である。Fermi の二状態の重み $1,e^{-\beta\xi}$ から $n=f(\xi)=1/(e^{\beta\xi}+1)$ を得る。Bose では $\xi>0$ を仮定し、等比級数の重みから $n=n_B(\xi)=1/(e^{\beta\xi}-1)$ となる。例えば $\beta\xi=\log 3$ ならそれぞれ $1/4$ と $1/2$ である。
:::

:::example[逆温度の尺度]
$k_BT=25\,\mathrm{meV}$ なら $\beta=40\,\mathrm{eV^{-1}}$ である。$\xi=25\log 3\,\mathrm{meV}$ の Fermi 軌道は確率 $1/4$ で占められる。温度を固定し、軌道が化学ポテンシャルを越える向きを反転させると占有率は $3/4$ になる。
:::

### 虚時間の端で測るエネルギー

平衡で $\tau>0$ とし、$\partial_\tau\mathcal G_{ij}=-\langle [K,a_i(\tau)]a_j^\dagger\rangle_\beta$ である。左極限でも同じ微分を負時間の順序で行うと、
$$
\partial_\tau\mathcal G(0^+)-\partial_\tau\mathcal G(0^-)
=\langle[[a_i,K],a_j^\dagger]_\eta\rangle_\beta
$$
を得る。すなわち値の跳びは正準代数を、傾きの跳びは相互作用を含む一次モーメントを検査する。Hubbard 原子で $K=\xi(n_\uparrow+n_\downarrow)+Un_\uparrow n_\downarrow$ なら、傾きの跳びは $\xi+U\langle n_\downarrow\rangle_\beta$ である。これは $[a_\uparrow,K]=(\xi+Un_\downarrow)a_\uparrow$ と、$n_\downarrow$ が上スピン場と可換であることから従う。

:::example[粒子・正孔対称点の傾き]
$\xi=-2\,\mathrm{eV}$、$U=4\,\mathrm{eV}$ とする。対称点では $\langle n_\downarrow\rangle_\beta=1/2$ なので傾きの跳びは零である。一方、二つの枝は $\pm2\,\mathrm{eV}$ にあり、二階微分の跳びは $\mathcal G$ の負符号を含めて $-4\,\mathrm{eV^2}$ である。一次モーメントが零でも、熱状態の相関が $\tau$ に依存しないわけではない。
:::

### Gibbs 状態と固定粒子数の違い

$[H,\hat N]=0$ で大正準の重みは $e^{-\beta H}e^{\beta\mu\hat N}$ に分かれる。各粒子数部門の射影を $P_N$、その分配関数を $Z_N=\Tr(P_N e^{-\beta H})$ とすれば
$$
Z=\sum_N e^{\beta\mu N}Z_N,\qquad
p_N=e^{\beta\mu N}Z_N/Z .
$$
粒子数を保つ演算子の平均は、各固定粒子数平均を $p_N$ で足したものである。一粒子相関の積全体は粒子数を保つが、二つの場の間の発展は $N\pm1$ 部門に入る。従って固定粒子数の相関では、熱重みのトレースの中の $P_N$ が各場と可換でない点に注意する。例えば $P_N a=aP_{N+1}$ であるから、トレースを循環して $a$ を熱因子の反対側へ移すと $P_{N+1}$ が現れる。規格化も $Z_N$ から $Z_{N+1}$ へ変わり得る。通常の松原反周期性は、$P_N$ を入れない大正準 Gibbs トレースの式として成立する。

:::example[固定した空状態で反周期性を試す]
一つの自由 Fermi 軌道について、固定粒子数 $N=0$ で平均する。$\tau>0$ では $\mathcal G=-e^{-\xi\tau}$、$\tau<0$ では $\mathcal G=0$ である。正時間の値を反周期に延長してもこの負時間の期待値にならない。大正準 Gibbs 状態では $p=\langle n\rangle=f(\xi)$ が有限であり、$0<\tau<\beta$ では $-(1-p)e^{-\xi\tau}$、負時間では $pe^{-\xi\tau}$ になる。$(1-p)e^{-\beta\xi}=p$ が成り立つため、正負の値は一つの反周期関数へつながる。
:::

### 虚時間発展を行列要素から安定に評価する

有限次元で $K\ket r=E_r\ket r$ とすると、正時間の熱トレースは
$$
-Z^{-1}\sum_{r,s}e^{-(\beta-\tau)E_r-\tau E_s}
 \mel r{a_i}s\mel s{a_j^\dagger}r .
$$
各エネルギーを一律に $E_{\min}$ だけ引くと、分子と分母は同じ $e^{-\beta E_{\min}}$ を持つ。従って $0\le\tau\le\beta$ の指数は全て非正にできる。これに対し $a_i(\tau)=e^{\tau K}a_i e^{-\tau K}$ を行列として先に作り、後で熱因子を掛ける手順では、相殺される大きな指数を別々に計算する。このためスペクトル和を使う有限模型では、区間の距離 $\tau,\beta-\tau$ の指数を組にする方が安定である。

:::example[一つの負エネルギーと端の占有]
$\xi=-E$、$E>0$ の一軌道では、空状態のエネルギーは選んだ基底から $E$ 高い。正時間の $\mathcal G$ は $-e^{-(\beta-\tau)E}/(1+e^{-\beta E})$ である。$\beta E=10$、$\tau=\beta/2$ なら絶対値は $e^{-5}/(1+e^{-10})$ である。一方、$\tau\to\beta^-$ では絶対値は $1/(1+e^{-10})$ へ近づく。区間の最後からの距離が短いときは、負エネルギーの軌道の占有が直接現れる。
:::

### 有限系で零温度極限を取る

大正準 Hamiltonian の基底エネルギーを零に選び、縮退度を $d_0$、次のエネルギーを $\Delta>0$ とする。有限次元 $d$ の場合、Gibbs 重みを取り分ければ $Z=d_0+Z_{\rm ex}$、$0\le Z_{\rm ex}\le(d-d_0)e^{-\beta\Delta}$ である。従って基底への射影の部分は $P_0/(d_0+Z_{\rm ex})$ であり、励起を含む確率は $Z_{\rm ex}/(d_0+Z_{\rm ex})$ 以下である。任意の有界な観測量 $O$ の平均と基底等重み平均との差は、両部門の平均の差の絶対値が $2\|O\|$ 以下であることから
$$
\left|\langle O\rangle_\beta-\frac{\Tr(P_0O)}{d_0}\right|
\le 2\|O\|\frac{Z_{\rm ex}}{d_0+Z_{\rm ex}}
$$
となる。二時刻の虚時間相関へこの評価を適用するには時刻の差を固定して、その場の積のノルムが温度とともに増えないようにする。例えば $\tau$ を $\beta/2$ に置いたまま $\beta$ を伸ばす操作では、二時刻の距離も増えるため、固定した $\tau$ の極限と同じ議論にはならない。

:::example[縮退を温度で分解する]
基底候補が二つあってエネルギーが $0,\delta>0$ であり、他の状態を無視できるとする。二番目の確率は $1/(1+e^{\beta\delta})$ である。$\beta\delta=\log3$ なら $1/4$、$\delta\to0$ を先に取れば $1/2$、$\beta\to\infty$ を先に取れば零である。このモデルでは弱い場が縮退を選ぶ場合と、場を零とした Gibbs 極限で基底を混合する場合で、熱平均の出発点が変わる。
:::

## 周期性・反周期性と松原周波数

### トレースの循環性

$0<\tau<\beta$ に対して
$$
\begin{aligned}
\Tr(e^{-\beta K}a_i(\tau)a_j^\dagger)
&=\Tr(e^{-(\beta-\tau)K}a_i e^{-\tau K}a_j^\dagger)\\
&=\Tr(e^{-\beta K}a_j^\dagger a_i(\tau-\beta)).
\end{aligned}
$$
ここで指数因子をまとめるには $K$ の関数どうしの可換性を用いた。したがって $\mathcal G(\tau+\beta)=\eta\mathcal G(\tau)$ が境界値を含めて成立する。トレースの順序の移動そのものには符号が付かない。$\eta$ が出るのは負の時間差で時間順序を定義するときである。

:::definition[松原周波数]{#def-mb-matsubara-frequency}
境界条件 $e^{-i\omega\beta}=\eta$ を満たす離散周波数を :term[松原周波数]{en="Matsubara frequency" yomi="まつばらしゅうはすう"} と呼ぶ。$\eta=-1$ では $\omega_n=(2n+1)\pi/\beta$、$\eta=+1$ では $\nu_m=2m\pi/\beta$、$n,m\in\Z$ である。
:::

同じ境界条件の指数関数の比は周期的であり、$[0,\beta]$ 上の内積は $\beta\delta_{mn}$ である。Fermi の系は周期的な Fourier 基底全体を $e^{-i\pi\tau/\beta}$ 倍して得るため、同じ区間で完全である。それゆえ
$$
\mathcal G(i\omega_n)=\int_0^\beta e^{i\omega_n\tau}\mathcal G(\tau)\dd\tau,\qquad
\mathcal G(\tau)=\frac1\beta\sum_n e^{-i\omega_n\tau}\mathcal G(i\omega_n).
\label{eq-mb-matsubara-fourier}
$$
右辺は不連続点の片側を求める場合、先に $\tau\ne0$ の値を再構成してから極限を取る。境界条件から $\gamma=-\eta\mathcal G(0^-)=-\mathcal G(\beta^-)$ である。周波数和で左極限を表すときは $e^{i\omega_n0^+}$ の収束因子を残す。

:::example[自由関数の Fourier 係数]
$0<\tau<\beta$ で $-(1+\eta n)e^{-\xi\tau}$ を積分する。$e^{i\omega_n\beta}=\eta$ と $(1+\eta n)(1-\eta e^{-\beta\xi})=1$ から $\mathcal G_0(i\omega_n)=1/(i\omega_n-\xi)$ になる。Fermi の $\xi=0$、$\beta=2$ なら最小正周波数の係数は $-2i/\pi$ である。
:::

:::example[密度の伝える周波数]
二本の Fermi 線の周波数差 $(2n+1)\pi/\beta-(2l+1)\pi/\beta$ は $2(n-l)\pi/\beta$ である。密度は Fermi 演算子二つの積で、虚時間境界条件は周期的である。このため密度応答には Bose 周波数を割り当てる。
:::

### 境界つきデルタ関数と自由逆行列

$\delta_\eta(\tau)=\beta^{-1}\sum_n e^{-i\omega_n\tau}$ を境界条件を満たす関数上の分布とする。試験関数を Fourier 展開すると、各係数に恒等的に作用するため、この分布を畳み込むと試験関数が戻る。自由 Green 関数は区間内部で $\partial_\tau\mathcal G=-\xi\mathcal G$ を満たすが、場の跳びが $-1$ なので
$$
(-\partial_\tau-\xi)\mathcal G(\tau)=\delta_\eta(\tau)
$$
である。従って $i\omega_n-\xi$ が各モードで非零なら、逆は一意である。$\eta=1,\xi=0$ では定数モードの逆が存在せず、周期箱で自由 Bose 零モードを積分するには正のギャップ、固定粒子数、または相互作用を考慮する必要がある。

:::example[虚時間上の平均を引く逆]
Bose の $\xi=0$ で定数モードを除いて $-\partial_\tau$ を逆にする。$0<\tau<\beta$ で $g(\tau)=\tau/\beta-1/2$ を周期延長すると、境界で $-1$ 跳び、区間の平均は零である。従って $-\partial_\tau g=\delta_+(\tau)-1/\beta$ であり、非零の松原係数は $1/(i\nu_m)$ となる。零点で関数値を選んでも、跳びや非零係数は変わらない。
:::

### Kubo–Martin–Schwinger 境界条件を演算子の順序で書く

時間に依存しない $A,B$ に対し、有限次元で $F_{AB}(z)=\langle A(z)B\rangle_\beta$、$A(z)=e^{iKz}Ae^{-iKz}$ とする。$z=t-i\beta$ を代入すると $A(t-i\beta)=e^{\beta K}A(t)e^{-\beta K}$ であり、
$$
F_{AB}(t-i\beta)=Z^{-1}\Tr[A(t)e^{-\beta K}B]
=\langle BA(t)\rangle_\beta .
$$
これが :term[Kubo–Martin–Schwinger 条件]{en="Kubo–Martin–Schwinger condition" yomi="くぼまるてぃんしゅうぃんがーじょうけん"} の有限次元での形である。この段階では通常のトレースを移しただけなので、Fermi 場にも符号を付けない。$F_{AB}(t-i\tau)$ は $0<\tau<\beta$ の帯内で熱指数 $e^{-(\beta-\tau)E_r-\tau E_s}$ を持つ。作用させる演算子の順序はこの帯の二つの境界で交換され、時間順序 Green 関数へ組み込むときに統計に応じた符号が付く。

:::example[自由 Fermi の大小の境界値]
$A=c,B=c^\dagger$ では $F_{AB}(t)=(1-f)e^{-i\xi t}$、$F_{AB}(t-i\beta)=f e^{-i\xi t}$ である。例えば $\beta\xi=\log 4$ では正時間の交換前の係数は $4/5$、交換後は $1/5$ である。実時間の $G^>=-iF_{AB}$、$G^<=i\langle B A(t)\rangle$ に移すと $G^<(t)=-e^{-\beta\xi}G^>(t)$ になる。熱の比は正であり、この式の負号は二つの実時間成分の定義に由来する。
:::

### 化学ポテンシャルを時間発展から外した場合

$H$ で虚時間を発展させる場を $a_H(\tau)$、$K$ で発展させる場を $a_K(\tau)$ とする。$[\hat N,a]=-a$ から $a_K(\tau)=e^{\mu\tau}a_H(\tau)$ である。従って同じ Gibbs 重みで期待値を取っても $H$ 発展の温度 Green 関数 $\mathcal G_H=e^{-\mu\tau}\mathcal G_K$ は
$$
\mathcal G_H(\tau+\beta)=\eta e^{-\beta\mu}\mathcal G_H(\tau)
$$
を満たす。反周期または周期の Fourier 基底をそのまま使うのは $\mathcal G_K$ の方である。$H$ 発展の指数 $e^{-(i\omega_n+\mu)\tau}$ はこの捻れた境界条件に合う。エネルギーの原点を変えるとき、分母の $\xi=\epsilon-\mu$ と虚時間の境界条件は同じ選択に合わせる必要がある。

:::example[同じ軌道の二種類の指数]
$\beta=1\,\mathrm{meV^{-1}}$、$\epsilon=\log 6\,\mathrm{meV}$、$\mu=\log 2\,\mathrm{meV}$ とする。$K$ の軌道エネルギーでは $e^{-\beta\xi}=1/3$、従って Fermi 占有は $1/4$ である。$H$ の正時間相関の指数は $\epsilon$ を用い、$\beta^-$ での値は $-1/8$ になる。$K$ のものは $\beta^-$ で $-1/4$ である。$H$ での末端値から密度を得るなら、境界の $e^{-\beta\mu}=1/2$ を戻してから読む。
:::

### Fourier 級数が跳びの中点を返す理由

反周期の $\mathcal G_F(\tau)$ に $e^{i\pi\tau/\beta}$ を掛けると周期関数になり、Matsubara 展開は通常の Fourier 展開へ移る。区間内で区分的に微分可能で、跳びがある場合、対称な Fourier 和は跳びの中点へ収束する。実際、Dirichlet 核は偶関数であり、左右の極限を引いた残りを核で積分した寄与は、原点近傍での差が $O(|\tau|)$ であれば振動積分として消える。従って一粒子で $\mathcal G(0^+)-\mathcal G(0^-)=-1$ のとき、対称和での値は $\mathcal G(0^-)-1/2$ である。特に $e^{i\omega_n\delta}$ を付けず単純に対称打ち切りを進める場合は、期待値の $0^-$ に対応する $1/2$ を戻す。

:::example[対称な軌道の Matsubara 項を組にする]
$\xi=0$ の Fermi では $\mathcal G(i\omega_n)=1/(i\omega_n)$ である。$n$ と $-n-1$ の二項は和が零なので、対称切断した値は常に零である。実空間の左極限は $+1/2$、右極限は $-1/2$ なので中点の零と一致する。ここから占有を零と結論すると片側と中点を取り違える。
:::

### 密度の時間順序と化学ポテンシャルの打ち消し

粒子数保存の二次場 $B=\sum_{ij}B_{ij}a_i^\dagger a_j$ は $[\hat N,B]=0$ を満たす。従って $H$ 発展と $K$ 発展のどちらでも $B(\tau)$ は同じである。これは各生成と消滅の因子 $e^{-\mu\tau},e^{+\mu\tau}$ が打ち消し合うことからも分かる。二つの密度の時間順序では軌道の間で四回の Fermi 場の交換になるため、時刻を入れ替えるときに負符号は付かない。従って外部の周波数は周期的な Bose 周波数であり、内部に現れる二本の線は Fermi の周波数を持つ。外部周波数の統計は構成粒子の統計の名前から直接決めるのではなく、交換で動かした場のパリティを数えて決める。

:::example[スピン反転の外部周波数]
$S^+=a_\uparrow^\dagger a_\downarrow$ とし、分散 $\xi_\uparrow,\xi_\downarrow$ の自由な二軌道を取る。$S^+(\tau)=e^{(\xi_\uparrow-\xi_\downarrow)\tau}S^+$ であり、$\mu$ は差から消える。正時間の相関 $\langle S^+(\tau)S^-\rangle$ は $e^{(\xi_\uparrow-\xi_\downarrow)\tau}f_\uparrow(1-f_\downarrow)$ である。$\tau=\beta^-$ では $e^{\beta(\xi_\uparrow-\xi_\downarrow)}f_\uparrow(1-f_\downarrow)=f_\downarrow(1-f_\uparrow)$ となり、負時間で二つの演算子を反対に置いた等時刻と一致する。例えば $e^{\beta\xi_\uparrow}=3$、$e^{\beta\xi_\downarrow}=1$ なら正時間の初めの係数は $1/8$、終わりは $3/8$ である。
:::

## 松原和の計算法

### 留数と収束因子

$f(z)=1/(e^{\beta z}+1)$ の $i\omega_n$ での留数は $-1/\beta$ である。$n_B(z)$ の $i\nu_m$ での留数は $+1/\beta$ である。$F$ が有理関数でこれらの点に極を持たず、無限遠で $O(z^{-2})$ とする。虚軸上の熱極の中間に水平辺を持つ長方形を広げると、熱関数は辺上で有界であり、積の輪郭積分は零に収束する。そこで Fermi の和は $\sum_{z_p}\operatorname{Res}_{z_p}(fF)$、Bose の和は $-\sum_{z_p}\operatorname{Res}_{z_p}(n_B F)$ となる。ここで $z_p$ は $F$ の極である。繰り返し極には留数の微分公式を使う。

$F(z)=1/(z-E)$ では絶対収束しない。自由関数の虚時間表示を $\tau=0^-$ に近づけると、因子 $e^{i\omega_n0^+}$ 付きの和は Fermi で $f(E)$、Bose で $-n_B(E)$ となる。収束因子を逆に選ぶと $0^+$ に対応し、どちらも値が $1$ だけ減る。

:::example[二本の分母]
$a\ne b$ に対して
$$
\frac1\beta\sum_n\frac1{(i\omega_n-a)(i\omega_n-b)}=\frac{f(a)-f(b)}{a-b}.
\label{eq-mb-matsubara-two-poles}
$$
$F$ の留数は $1/(a-b)$ と $1/(b-a)$ である。$b\to a$ のとき $f'(a)=-\beta f(a)(1-f(a))$ になる。例えば $\beta=2$、$a=b=0$ では $-1/2$ である。
:::

:::example[対を作る二本の線]
$a=E,b=-E$、$E>0$ とすると、分母は $-\omega_n^2-E^2$ である。よって
$$
\frac1\beta\sum_n\frac1{\omega_n^2+E^2}=\frac{\tanh(\beta E/2)}{2E}.
\label{eq-mb-pair-sum}
$$
Bose 周波数については $n_B(-E)=-1-n_B(E)$ を用いて同じ和が $\coth(\beta E/2)/(2E)$ となる。$\beta E=\log 3$ なら Fermi の双曲線正接は $1/2$、Bose の双曲線余接は $2$ である。
:::

外から Bose 周波数 $i\nu_m$ が入る場合、$F(z)=[(z-a)(z+i\nu_m-b)]^{-1}$ の後者の極に $f(b-i\nu_m)=f(b)$ を用いる。和は $(f(a)-f(b))/(a-b+i\nu_m)$ になる。和を先に評価した上で外部周波数を連続変数へ移す。途中で $i\nu_m$ を任意の複素数に変えると、ここで使った熱関数の周期性が失われる。

### 電荷の静的ゆらぎと自由エネルギー

自由軌道に静的な準位変化 $\lambda n$ を加えると熱自由エネルギーは $\Omega(\lambda)=-\beta^{-1}\log(1+e^{-\beta(\xi+\lambda)})$ である。微分は $\partial_\lambda\Omega=f(\xi+\lambda)$、二階微分は $-\beta f(1-f)$ になる。一方、二本の自由線を同じ周波数に取った松原和は式 :eqref[eq-mb-matsubara-two-poles] で $a=b=\xi$ の極限なので、同じ値を与える。零の外部 Matsubara 周波数を基準に取ることがここで重要である。非零の外部周波数で同じ $a=b$ を代入すれば、分子が零になる。

:::example[静的な一点の応答]
$\beta=2$、$\xi=0$ の一軌道の熱揺らぎは $\operatorname{Var}(n)=1/4$、静的な応答は $-1/2$ である。この数を非零の Bose 周波数でのバブルへ加える理由はない。保存量の動的応答にはこのゆらぎが現れない。
:::

### 近接する極の処理と高周波尾部

$a_1,\ldots,a_r$ が互いに異なる実数で $r\ge2$ なら、Fermi の和は
$$
\frac1\beta\sum_n\prod_{j=1}^r(i\omega_n-a_j)^{-1}
=\sum_{j=1}^r \frac{f(a_j)}{\prod_{l\ne j}(a_j-a_l)} .
$$
各項の係数は $z=a_j$ で残る他の分母の逆である。部分分数の和を元の分母で通分すると、分子は各 $a_j$ で $1$ となる。分子と $1$ の差は次数 $r-1$ 以下で根が $r$ 個あるので恒等的に零である。極が重なるときは $f$ の微分へ合流させる。全て $a$ に重なる場合は $f^{(r-1)}(a)/(r-1)!$ である。分母を $a$ で微分した結果が $(r-1)!(z-a)^{-r}$ であり、$r\ge2$ の導関数の和はコンパクトな実軸上で絶対一様収束する。数値計算では近接した $a_j$ の個別の留数が大きくなるため、合流の式が桁落ちを避ける。

:::example[三重極が零になる位置]
$r=3$ では $f''=\beta^2 f(1-f)(1-2f)$ である。半充填の $a=0$ で $f=1/2$ なので、三重の和は零となる。各 Matsubara 項でも $n$ と $-n-1$ を組にすると $1/(i\omega_n)^3$ が相殺する。$\beta=2$、$a=(\log3)/2$ では $f=1/4$ なので $f''/2=3/16$ になる。
:::

有限のカットオフで密度を計算するには、虚軸上で $G(z)=I/z+M_1/z^2+O(z^{-3})$ と展開し、最初の二項を引く。例えば $1/z$ の Fermi 等時刻和は $1/2$、$1/z^2$ の和は $f'(0)=-\beta/4$ である。従って引いた二項の寄与は $1/2-\beta M_1/4$ と解析的に足せる。引いた後の尾部を絶対値で評価すると $\sum_{|n|>N}|n|^{-3}=O(N^{-2})$ となる。$M_1$ の符号を誤ると $1/z^2$ が残るため、打ち切りの誤差は $O(N^{-1})$ に戻る。

### Bose の二分母と零モードの扱い

$E>0$ として $\nu_m=2\pi m T$ を用いる。$1/(\nu_m^2+E^2)$ を $z=i\nu_m$ による部分分数へ変え、Bose の留数式を使うと
$$
T\sum_m\frac1{\nu_m^2+E^2}=\frac{1+2n_B(E)}{2E}.
$$
零モードは $T/E^2$ なので、$E/T\ll1$ で発散する部分が全てこの一項にある。非零モードの和は $E=0$ を項ごとに取れる。実際、各項は $T/\nu_m^2$ で支配され、その和は $2T(2\pi T)^{-2}\sum_{m\ge1}m^{-2}=1/(12T)$ である。従って $T$ を固定して $E\to0$ とすると、非零モードの先頭は有限で、零モードの扱いを分けることで数値和の特異性を取り除ける。Fermi の対の和には零モードがなく、その小さな $E$ の極限は $1/(4T)$ である。

:::example[古典的な軟モードと Fermi 対を比較する]
$T=1\,\mathrm{meV}$、$E=0.1\,\mathrm{meV}$ なら Bose の零モードは $100\,\mathrm{meV^{-1}}$、非零モードは先頭で約 $0.0833\,\mathrm{meV^{-1}}$ である。同じ $E,T$ の Fermi の対の和は $\tanh(0.05)/(0.2\,\mathrm{meV})\simeq0.250\,\mathrm{meV^{-1}}$ であり、軟化で $E^{-2}$ は現れない。有限温度で低い周波数だけを先に扱うとき、この差が有効理論に残す自由度を決める。
:::

### 混合する周波数の和を一つ評価する

外部が Fermi 周波数 $i\omega_n$、内部が Bose 周波数 $i\nu_m$ とする。電子・Bose の一線ずつについて
$$
S(i\omega_n)=T\sum_m
\frac{1}{i\nu_m-b}\frac{1}{i\omega_n-i\nu_m-a}
$$
を考える。$b$ は零でない実数、$a$ は実数である。内部変数 $z$ の二極は $b$ と $i\omega_n-a$ であり、留数はそれぞれ $(i\omega_n-a-b)^{-1}$ とその負である。Bose の式では各留数に $-n_B(z_p)$ を掛ける。$e^{\beta(i\omega_n-a)}=-e^{-\beta a}$ なので、$n_B(i\omega_n-a)=-(1-f(a))$ である。従って
$$
S(i\omega_n)=-\frac{n_B(b)+1-f(a)}{i\omega_n-a-b}.
$$
変位の自由関数では $b>0$ の線と $-b$ の線を差し引く。後者に $n_B(-b)=-1-n_B(b)$ を入れると、二つの分子が $n_B(b)+1-f(a)$ と $n_B(b)+f(a)$ に分かれる。全体の符号は、図で変位の温度関数を $-\langle T_\tau XX\rangle$ とするか、正の Gaussian 共分散を使うかに合わせる。

:::example[低温で放出と吸収を区別する]
$b=2\,\mathrm{meV}>0$、$a=1\,\mathrm{meV}$ とし、$T\to0$ では $n_B(b)\to0$、$f(a)\to0$ である。従って $b$ の項の分子は $1$、$-b$ の項の分子は $0$ となる。$a=-1\,\mathrm{meV}$ に変えると $f(a)\to1$ なので二つの役割が交換する。$T>0$ ではどちらの分子も非零になり、熱により振動を吸収する過程も残る。
:::

### 微分で高次の松原和を得る

$J(E)=T\sum_n(\omega_n^2+E^2)^{-1}$ とする。$E>0$ で各項を微分すると $-2E(\omega_n^2+E^2)^{-2}$ であり、コンパクトな $E$ の範囲で一様収束する。従って
$$
T\sum_n\frac{1}{(\omega_n^2+E^2)^2}
=-\frac{J'(E)}{2E}
=\frac{\tanh(\beta E/2)}{4E^3}
-\frac{\beta\operatorname{sech}^2(\beta E/2)}{8E^2}.
$$
この式は二つの項が共に $E^{-2}$ のように見えるが、$E\to0$ で相殺する。$x=\beta E/2$ で $\tanh x=x-x^3/3+O(x^5)$ と $\operatorname{sech}^2x=1-x^2+O(x^4)$ を代入すると、極限は $\beta^3/48$ である。元の和では Fermi の周波数に零がないため、極限が有限となることも直接に分かる。大きな項の差を数値で評価する場合は、小さな $E$ のところでこの Taylor 展開を使う。

:::example[四線の和を臨界点で使う]
$\beta=2\,\mathrm{meV^{-1}}$ なら $E\to0$ の四線の和は $1/6\,\mathrm{meV^{-3}}$ である。$T=0$ で $E>0$ を保持すると二項目が消え、値は $1/(4E^3)$ になる。温度を先に零へ近づけた計算を、有限温度の超伝導の臨界点の係数へ代入すると赤外の依存性を取り違える。
:::

### Bose 零モードを含む和を時間分割で求める

自由な実変位の二次作用の逆核を考える。$g_\beta(\tau)=T\sum_m e^{-i\nu_m\tau}/(\nu_m^2+E^2)$ と置くと、$(-\partial_\tau^2+E^2)g_\beta=\delta_+(\tau)$ である。区間内部では斉次方程式なので、$\tau$ と $\beta-\tau$ の対称性と周期境界を使えば $g_\beta=C(e^{-E\tau}+e^{-E(\beta-\tau)})$ と書ける。境界で微分の跳びは $g_\beta'(0^+)-g_\beta'(\beta^-)=-1$ である。従って $C=1/[2E(1-e^{-\beta E})]$ となり、
$$
g_\beta(0)=\frac{1+e^{-\beta E}}{2E(1-e^{-\beta E})}
=\frac{\coth(\beta E/2)}{2E}.
$$
この方法では、零モードを含む周波数和と同じ答えを、虚時間の端の微分条件から求められる。$\tau=\beta/2$ では分子の二つの指数が同じであり、相関が端と中央でどれだけ異なるかも分かる。零モードを除く場合は $T/E^2$ を全時刻で引く。

:::example[熱の両端が近いとき]
$\beta E=\log4$ なら $e^{-\beta E}=1/4$ であり、端は $g_\beta(0)=5/(6E)$、中央は $g_\beta(\beta/2)=2/(3E)$ である。比は $4/5$ であり、長い虚時間でも一方向だけの指数 $e^{-E\beta/2}=1/2$ より大きい。周期境界で逆向きに伝わった項が同時に寄与するためである。
:::

## 解析接続と遅延 Green 関数

### スペクトルによる接続

$K\ket r=E_r\ket r$、$p_r=e^{-\beta E_r}/Z$ とする。正の虚時間の相関に固有状態を挿入すると
$$
\mathcal G_{ij}(\tau)=-\sum_{rs}p_r e^{-\tau(E_s-E_r)}
\mel r{a_i}s\mel s{a_j^\dagger}r .
$$
積分では $1-\eta e^{-\beta(E_s-E_r)}$ が現れる。そのため周波数の分母は $i\omega_n-(E_s-E_r)$、分子の熱重みは $p_r-\eta p_s$ である。Fermi 系では全ての松原点で、Bose 系では少なくとも非零の松原点で、これが :eqref[eq-mb-green-resolvent] と同じ形を与える。Bose の零松原点で同じ置換を使うには、等エネルギー部門の行列要素を処理する必要がある。$G(z)$ の上半平面境界値 $z\to\omega+i0$ を取れば遅延関数になる。有限温度でも、正の実時間と実軸上のスペクトルを使った Laplace 積分が同じ表示を与える。正の松原周波数は遅延側、負の松原周波数は先進側の半平面に属し、$G(-i\omega)=G(i\omega)^\dagger$ で結ばれる。実軸に近づく際は、実周波数の正負にかかわらず遅延側では $+i0$ を取る。

:::remark[静的な成分]
Bose 型相関に $E_r=E_s$ の行列要素があると、$0<\tau<\beta$ の定数項は $i\nu_m=0$ だけに寄与する。交換子のスペクトルには $p_r-p_s=0$ によりこの成分が現れない。従って保存量の平衡ゆらぎとその遅延応答は区別する。
:::

:::example[極を持つ関数と補間]
$\mathcal G(i\omega_n)=1/(i\omega_n-2)$ が自由軌道の表示から出たと分かっているなら、遅延関数は $1/(\omega-2+i0)$ である。離散点に値を合わせるだけでは、例えば $\cosh(\beta z/2)$ を加えた関数も全ての Fermi 点で同じ値を取る。後者はスペクトル表示の $1/z$ 漸近性を壊すため、物理的な補間から外れる。
:::

:::example[二つの準位の遷移]
占有差が $f(a)-f(b)=1/2$、準位差 $b-a=2$ なら、二本の Fermi 線の和は $1/[2(i\nu_m-2)]$ である。上から接続するとその虚部は $-\pi\delta(\omega-2)/2$ になる。逆向き遷移を持つ Hermite な演算子では、負の周波数の重みも加わる。
:::

### 有限個の極を持つ場合の一意性

有限体積かつ有限次元で、$K$ の固有値を全て含む表示から $\mathcal G(z)$ は有限個の実軸の極を持つ有理関数である。このクラスで二つの候補がすべての Matsubara 点で同じ値なら、その差の分母を共通にした分子は無限個の相異なる根を持つ。有限次の多項式ではこれは零の場合だけである。従ってこのクラスの中では接続が一意である。Bose の零周波数に独立な熱成分がある場合は、非零の周波数列で同じ議論を適用し、零モードを別に扱う。連続スペクトルを持つ系では有理関数性を仮定できず、数値的には有限個の誤差つきの値から接続するため、さらに情報が要る。

:::example[一つの点が残す自由度]
$z=i$ でだけ自由軌道 $1/(z-1)$ と合う有理関数なら、例えば $\lambda(z^2+1)/[(z-2)(z-3)(z-4)(z-5)]$ を足せる。共役な点 $z=-i$ でもこの候補は一致し、$1/z$ の先頭係数も保つ。新しい線の重みが非負になるかは別の制約である。有限個の点だけでは、候補クラスの自由度を先に固定する必要がある。
:::

### 正の有限スペクトルによる補間の非一意性

$G_1(z)=1/z$ と比較して、規格化した対称二線 $A_a=(\delta_{-a}+\delta_a)/2$ の関数は $G_a(z)=z/(z^2-a^2)$ である。$G_a(i\omega)=-i\omega/(\omega^2+a^2)$ なので、$|G_a(i\omega)-G_1(i\omega)|=a^2/[|\omega|(\omega^2+a^2)]$ である。$\omega\ge\pi/\beta$ のデータに有限精度しかないと、十分小さな $a$ では片方を他方から区別できない。一方、厳密に同じ一点のデータをもつ二つの正スペクトルも作れる。例えば $\omega=1$ に対して $A_1$ と $(3A_0+5A_2)/8$ はともに $G(i)=-i/2$ である。後者では $-i$ と $-i/5$ をそれぞれ $3/8,5/8$ で平均するからである。どちらも総重みは $1$、一次モーメントは零である。

### 有限極の接続と実軸の極処方

$\mathcal G(i\omega_n)=\sum_{\alpha=1}^M w_\alpha/(i\omega_n-E_\alpha)$、$E_\alpha\in\R$ とする。正の松原点は上半平面の点列なので、その半平面で同じ分数和を使い実軸へ近づける。このときどの $E_\alpha$ にも $+i0$ を付ける。負の松原点は下半平面から得られ、実軸への接近は $-i0$ である。これは正または負の $\omega$ の値を選ぶ話ではなく、複素平面のどちらの側で $G$ が定義されているかの違いである。連続帯に対しても積分の分母で同じ区別をする。

:::example[一つの正周波数と一つの負周波数の線]
規格化した Fermi スペクトルを $\epsilon=-1,2$ の二線とし、重みを $1/4,3/4$ とする。$G(i)$ の第一項は $(1-i)/8$、第二項は $-3(2+i)/20$ なので合計 $-7/40-11i/40$ である。$G(-i)$ はその複素共役である。実時間の遅延関数では $-1$ と $2$ のどちらも $-\Im G^R/\pi$ に正の線を作る。零温度の時間順序関数に移す場合だけ、除去線に $-i0$ が入り虚部の符号が逆になる。
:::

### 虚時間での積分核と復元の感度

Fermi の正時間表示は
$$
-\mathcal G(\tau)=\int_\R
\frac{e^{-\tau\epsilon}}{1+e^{-\beta\epsilon}}A(\epsilon)\dd\epsilon,\qquad 0<\tau<\beta
$$
である。有限温度の Lehmann 和で分子の $p_r$ を $(p_r+p_s)/(1+e^{-\beta\epsilon})$ に直せば得られる。正の高エネルギー線では $e^{-\tau\epsilon}$、負の高エネルギー線では $e^{-(\beta-\tau)|\epsilon|}$ が重みを抑える。従って区間の中央だけを精密に測定しても、高いエネルギーの線を見分けにくい。一方、両端に近い点は広いエネルギーを測り、粒子の付加と除去の総量を制約する。端の位置を動かしたときに同じ雑音でどれほどピークが見えるかは、この積分核から比較できる。

:::example[中心は正負の線を取り違えやすい]
規格化した一本の線が $+E$ にある場合と $-E$ にある場合を比べる。$\tau=\beta/2$ の核はどちらも $1/[2\cosh(\beta E/2)]$ なので完全に一致する。時間の微分は $-E$ または $+E$ が掛かり、符号を区別する。例えば $\beta E=\log9$ なら中心での $-\mathcal G$ は $3/10$ である。$\tau=0^+$ では $+E$ の線が $9/10$、$-E$ の線が $1/10$ となり、中心だけと端も測る場合で情報量が異なる。
:::

### 行列のスペクトルを復元するときの制約

Fermi の行列値スペクトルの重みは各エネルギー領域で正半定値なので、二軌道の重みを $A_{ij}$ としたとき $|A_{12}|^2\le A_{11}A_{22}$ が成り立つ。実際、$A$ を正半定値の平方根で分解すれば $A_{12}$ は二つの列の内積であり、Cauchy–Schwarz が使える。従って $A_{11}=0$ の領域では非対角の線も零である。さらに $\int A_{12}=0$ は異なる軌道の等時刻反交換子から来る。二つの対角成分だけを独立に合わせてから非対角を任意に加えると、ある軌道の重ね合わせに負のスペクトルが出ることがある。

:::example[同じ対角線で悪い非対角を検出する]
二つの周波数の各々で $A_{11}=A_{22}=1/2$、$A_{12}=\pm3/4$ と置くと、総和は $I$ になる。しかし各行列の固有値は $5/4,-1/4$ なので、負の重みを持つ重ね合わせがある。$3/4$ を $1/2$ へ弱めると固有値は $1,0$ になり、各線は直交する一軌道の投影を表せる。総和則だけでは両者を区別できず、線ごとの正値性が必要である。
:::

### 全 Matsubara 点が連続スペクトルを一意に決める条件

正規化した有限の Fermi スペクトル測度二つが同じ全ての Matsubara 係数を持つとする。その差を有限符号付き測度 $D$ とする。虚時間区間の中の関数
$$
F(\tau)=\int_\R\frac{e^{-\tau E}}{1+e^{-\beta E}}\dd D(E)
$$
の全ての反周期 Fourier 係数が零になる。実際、各 $E$ で $\tau$ の積分を先に取るとスペクトルの分数和が得られ、核は区間内で有界なので積分を交換できる。基底の完全性から $F$ は区間内で零に等しく、連続だから全ての点で零である。次に $\tau$ を複素変数とする。帯 $0<\Re\tau<\beta$ のコンパクト集合では、$E$ の正負に応じて $e^{-\delta|E|}$ の抑制があり、$\tau$ 微分も $|E|e^{-\delta|E|}$ で有界である。従って帯内で解析的であり、実区間で零なので恒等定理により帯全体で零である。

特に $\tau=\beta/2+it$ で考えると、$F$ は有限測度 $\dd M(E)=[2\cosh(\beta E/2)]^{-1}\dd D(E)$ の Fourier 変換である。この変換が全ての $t$ で零なら $M=0$ である。実際、標準偏差 $\sigma>0$ の規格化 Gaussian $\rho_\sigma$ と $M$ の畳み込みの Fourier 表現は、変換 $F$ に $e^{-\sigma^2t^2/2}$ を掛けて逆変換したものであり、積分順序を交換できるので各 $\sigma$ で零である。任意の連続でコンパクト台の試験関数と積分し、$\sigma\to0$ で Gaussian がその関数を一様に近似することを使うと、元の $M$ も全ての試験関数に対し零である。重み $1/[2\cosh(\beta E/2)]$ はどの有限 $E$ でも正であるから、各コンパクト区間で逆重みを掛けて $D=0$ を得る。行列値のスペクトルにも各成分の有限測度に同じ証明が使える。

この一意性は、全ての Matsubara 点を誤差なく与え、スペクトル積分を持つクラスへ絞った主張である。有限点の補間例では Fourier 係数の残りを零と置けないので、帯内で $F$ が恒等的に消える段階に進めない。増えた情報を何点にするかだけでなく、端の和則や測定誤差の大きさを指定して復元法を選ぶ。

## スペクトル表示と総和則

### 熱重みと総和則

有限温度の一粒子密度を明示すると
$$
A_{ij}(\epsilon)=\sum_{rs}(p_r-\eta p_s)
\mel r{a_i}s\mel s{a_j^\dagger}r\,\delta(\epsilon-E_s+E_r).
\label{eq-mb-thermal-spectral}
$$
$r,s$ の完全性を使うと $\int A(\epsilon)\dd\epsilon=\langle[a_i,a_j^\dagger]_\eta\rangle_\beta=\delta_{ij}$ である。Fermi の場合は $\epsilon=0$ も含めて $f(\epsilon)(p_r+p_s)=p_s$ なので、$\gamma_{ij}=\int f(\epsilon)A_{ij}(\epsilon)\dd\epsilon$ となる。Bose の場合、$n_B(\epsilon)(p_r-p_s)=p_s$ を使えるのは $\epsilon\ne0$ である。よって $\int_{\epsilon\ne0}n_B(\epsilon)A_{ij}(\epsilon)\dd\epsilon$ に、$E_s=E_r$ に限った和 $\gamma^{(0)}_{ij}=\sum_{rs:E_s=E_r}p_s\,\mel r{a_i}s\mel s{a_j^\dagger}r$ を足して密度を得る。対応する虚時間相関は定数 $-\gamma^{(0)}$ で、松原関数には $-\beta\gamma^{(0)}\delta_{m0}$ を与える。

$1/(z-\epsilon)=1/z+\epsilon/z^2+\epsilon^2/[z^2(z-\epsilon)]$ なので、有限の二次絶対モーメントがあれば虚軸上で
$$
\mathcal G_{ij}(z)=\frac{\delta_{ij}}z+
\frac{\langle[[a_i,K],a_j^\dagger]_\eta\rangle_\beta}{z^2}
+O(|z|^{-3}).
$$
分子の二重積は $(E_s-E_r)\mel r{a_i}s$ が $\mel r{[a_i,K]}s$ に等しいことから得た。

:::example[有限温度の Hubbard 原子]
$K=-\mu(n_\uparrow+n_\downarrow)+Un_\uparrow n_\downarrow$ で $x=e^{\beta\mu}$、$y=e^{-\beta U}$ とすると $Z=1+2x+x^2y$、$p=\langle n_\downarrow\rangle_\beta=(x+x^2y)/Z$ である。上スピンの $A$ は $\epsilon=-\mu,U-\mu$ に $1-p,p$ の重みを持つ。$\mu=U/2$ では $x^2y=1$ なので、どの温度でも両者は $1/2$ である。
:::

:::example[安定な Bose 準位の零周波数遷移]
一モードで $K=U n(n-1)/2$、$U>0$ とする。$n=0,1$ はともにエネルギー零だが、$n\ge2$ ではエネルギーが増えるため有限温度のトレースは収束する。$n=l$ と $l+1$ の間の消滅演算子の行列要素の二乗は $l+1$、エネルギー差は $lU$ である。従って交換子スペクトルの重みは $(l+1)(p_l-p_{l+1})$ であり、$l=0$ は $p_0=p_1$ により消える。一方、虚時間の順序相関にはこの遷移から $-p_0$ の定数が残る。分配関数 $Z=\sum_l e^{-\beta U l(l-1)/2}$ に対し、零松原点の追加項は $-\beta/Z$、占有に足す重みは $1/Z$ である。自由な Bose の零エネルギー準位と違い、この零モードは相互作用で熱トレースが収束した状態の内部にある。
:::

:::example[占有と正負の線]
$\beta=\log 3$、Fermi の二線 $\epsilon=\pm1$、重みをそれぞれ $1/2$ とする。占有数は $(f(1)+f(-1))/2=1/2$ である。一方で正の線も $f(1)=1/4$ の係数で占有に寄与するので、負周波数だけを積分する零温度の公式とは各線の寄与が異なる。
:::

### 詳細釣り合いと占有の情報

Fermi の $C^>(\epsilon)=\sum_{rs}p_r|\mel r a s|^2\delta(\epsilon-E_s+E_r)$、$C^<(\epsilon)$ は同じ和の $p_r$ を $p_s$ に置き換えたものとする。デルタの支持で $p_s=e^{-\beta\epsilon}p_r$ なので、$C^<=e^{-\beta\epsilon}C^>$、$A=C^>+C^<$ となる。従って $C^>= (1-f)A$、$C^<=fA$ である。周波数の分解は正の測度について成立し、線の重なりがあっても各線を別々に同定する必要はない。

:::example[同じピークの出入りを比べる]
$\beta\epsilon=\log 7$ の位置に積分重み $3/4$ のピークがある。$f=1/8$ より空き状態の寄与は $21/32$、占有の寄与は $3/32$ である。両者の比は $7$ となる。同じ $\epsilon$ で温度を二倍にすると $\beta\epsilon$ は半分になり、重みの比は $\sqrt7$ になる。相互作用系では $A$ 自体にも温度依存があり、この比は同じ温度で二つの測定を比べるときに使う。
:::

### 収束因子と熱エネルギーの一次モーメント

二体 Fermi Hamiltonian の高周波係数を求める。$[a_i,K]=\sum_j h_{ij}a_j+\frac12\sum_{jkl}\bar v_{ij;kl}a_j^\dagger a_l a_k$ である。ここで $h$ は化学ポテンシャルを含む。反交換子を $a_m^\dagger$ と取ると、三因子の項では消滅二つのそれぞれからデルタを拾う。実際、
$$
\{a_j^\dagger a_l a_k,a_m^\dagger\}
=\delta_{km} a_j^\dagger a_l-\delta_{lm}a_j^\dagger a_k .
$$
$\bar v_{ij;kl}=-\bar v_{ij;lk}$ を用いると両者は同じになり、$1/2$ を消す。従って $M_{1,im}=h_{im}+\sum_{jl}\bar v_{ij;ml}\gamma_{lj}$ である。この式は相互作用を自由平均で置き換えた近似ではない。$\gamma$ が相互作用を含む熱平均である限り、一次モーメントの正確な演算子関係である。他方、Green 関数の全周波数を定数の $M_1-h$ だけで補正するのは、二次以上のモーメントを一致させるとは限らない。

:::example[原子で一次と二次の検査を分ける]
原子の上スピンの枝が $\xi,\xi+U$、重み $1-p,p$ なら、$M_1=\xi+pU$、$M_2=\xi^2+2pU\xi+pU^2$ である。平均場の一本の準位を $\xi+pU$ に置くと、一次は一致するが二次の差は $p(1-p)U^2$ である。$\xi=-1\,\mathrm{meV}$、$U=4\,\mathrm{meV}$、$p=1/2$ では一次が $1\,\mathrm{meV}$、厳密の二次が $5\,\mathrm{meV^2}$、一本の近似が $1\,\mathrm{meV^2}$ である。熱により占有の確率 $p$ が変わると、差も $p(1-p)$ に応じて変わる。
:::

### 熱微分から二体相関を取り出す

$H(\lambda)=H_0+\lambda V$、粒子数を保ち、$V$ と $H_0$ が可換とは限らないとする。指数の微分は
$$
\partial_\lambda e^{-\beta K}
=-\int_0^\beta e^{-(\beta-\tau)K}V e^{-\tau K}\dd\tau
$$
である。有限行列では常微分方程式 $\partial_\beta e^{-\beta K}=-K e^{-\beta K}$ を $\lambda$ で微分し、初期値零の方程式を積分するとこの式を得る。トレース内では循環して全ての $\tau$ で同じ $-\Tr(e^{-\beta K}V)$ になる。従って $\partial_\lambda\Omega=\langle V\rangle$、$\Omega=-T\log Z$ である。もう一度微分すると
$$
\partial_\lambda^2\Omega=-\int_0^\beta
\langle \delta V(\tau)\delta V(0)\rangle_\beta\dd\tau.
$$
可換であれば $V(\tau)=V$ なので $-\beta\operatorname{Var} V$ である。可換でないとき、同じ時間の分散だけを掛けると遷移エネルギーの分母を落とす。

:::example[原子の二重占有から自由エネルギーを微分する]
$\mu$ を固定して Hubbard 原子の $U$ を動かす。$x=e^{\beta\mu}$、$y=e^{-\beta U}$ とすると $\partial_U\Omega=x^2y/(1+2x+x^2y)=d$ である。$d$ は二重占有の確率であり、二階微分は $-\beta d(1-d)$ である。例えば $x=2,y=1/4$ では $Z=6,d=1/6$ なので二階微分は $-5\beta/36$ となる。ここで計算したのは $\mu$ を固定して $U$ だけを変える微分であり、粒子・正孔対称の $\mu=U/2$ に沿う微分には粒子数項が加わる。
:::

### Bose の負周波数の符号を一振動子で確認する

実変位 $X=b+b^\dagger$ の自由 Hamiltonian $K=E b^\dagger b$、$E>0$ を取る。$X$ は粒子数を変える成分を両方含むため、交換子スペクトルは
$$
A_X(\epsilon)=\delta(\epsilon-E)-\delta(\epsilon+E)
$$
である。実際に $[X(t),X(0)]=e^{-iEt}-e^{iEt}$ であり、等時刻の交換子が零なので総重みも零である。正負の線を両方含めると一次モーメントは $2E$ である。負周波数の符号を正に取り替えると、総和則と高周波の次数の両方が変わってしまう。

等時刻の二点は $\langle X^2\rangle_\beta=\int n_B(\epsilon)A_X(\epsilon)\dd\epsilon$ である。$n_B(-E)=-[1+n_B(E)]$ なので、これは $1+2n_B(E)$ となる。負周波数の Bose 分布だけを確率として扱うことはできず、スペクトルの負号と積にして交換子から順序を戻す。例えば $\beta E=\log3$ なら $n_B(E)=1/2$、$n_B(-E)=-3/2$ であり、二つの寄与の和は $2$ である。$T\to0$ でも負周波数側が $1$ を与え、変位の零点ゆらぎが残る。

### 変位の高周波尾部と線形な場の違い

同じ振動子の消滅場 $b$ の交換子スペクトルは $\delta(\epsilon-E)$ だけで総重み $1$、従って $\mathcal G_b(z)\sim1/z$ である。一方、変位の関数は $1/(z-E)-1/(z+E)=2E/(z^2-E^2)$ なので $2E/z^2$ から始まる。これは $b$ と $b^\dagger$ の両方の符号を含む変位の等時刻交換から従う。従って松原の高周波側を差し引くとき、同じ振動子であっても場の選び方によって差し引く次数が異なる。例えば $E=2\,\mathrm{meV}$ の $D_X(i\nu)=-4/(\nu^2+4)$ は先頭から負であり、$1/(i\nu)$ の奇な尾を足して総和をとる必要はない。

## 有限温度の Wick の定理

### 熱平均の帰納法

$K_0=\sum_j\xi_j a_j^\dagger a_j+C$ とし、Bose では $\xi_j>0$ とする。$q_j=e^{-\beta\xi_j}$ と書く。$a_j e^{-\beta K_0}=q_j e^{-\beta K_0}a_j$ なので
$\langle X a_j\rangle_\beta=q_j\langle a_j X\rangle_\beta$ である。生成演算子の場合は $q_j^{-1}$ が対応する。

:::theorem[有限温度の Wick の定理]{#thm-mb-thermal-wick}
上の二次 Hamiltonian の熱状態で、生成・消滅演算子の線形結合の時間順序積の期待値は、全ての完全対の選び方について二因子の時間順序平均を掛けて和を取ったものに等しい。Fermi 系では対を隣接させる置換の符号を含める。奇数個の積は零である。
:::

:::proof
時間を先に並べ、最初の因子を固有モードの消滅演算子 $X_1=a_j$ とする。残り $m-1$ 因子の右に運ぶ際に数 $[X_1,X_l]_\eta$ と交換符号が出る。最後まで運んだ項の熱平均は、循環性で $\eta^{m-1}q_j$ 倍の元の積に戻る。$m$ 偶数では左辺に移して $1-\eta q_j$ を得る。二因子だけに同じ手順を使えば $\langle X_1X_l\rangle_\beta=[X_1,X_l]_\eta/(1-\eta q_j)$ である。それゆえ
$$
\langle X_1\cdots X_m\rangle_\beta
=\sum_{l=2}^m\eta^{l-2}\langle X_1X_l\rangle_\beta
\langle X_2\cdots X_{l-1}X_{l+1}\cdots X_m\rangle_\beta .
$$
生成では $q_j^{-1}$ に置き換えれば分母は同様に非零である。奇数個の単項式は少なくとも一モードの正味の生成数が零でないため、占有数基底での対角要素が消える。線形性と偶数の再帰で任意の時間順序積に拡張できる。
:::

:::example[熱揺らぎの違い]
一モード Bose 状態で $\langle a^\dagger a^\dagger aa\rangle_\beta=2n^2$ である。$n^2=a^\dagger a^\dagger aa+n$ の演算子関係から占有数の分散は $n(1+n)$ となる。Fermi では占有数演算子が射影なので分散は $n(1-n)$ である。平均 $n=1/2$ なら分散はそれぞれ $3/4,1/4$ である。
:::

:::example[相互作用のある Gibbs 状態]
半充填の Hubbard 原子で $U>0$ とする。二重占有確率は $1/[2(1+e^{\beta U/2})]$ であり、$U=0$ で $1/4$、$\beta U/2=\log 3$ で $1/8$ となる。一体密度は各スピンで常に $1/2$ なので、この差は二次の熱状態を仮定した Wick 分解が相互作用状態では成立しないことを示す。
:::

### 熱 Gaussian 状態の数の生成関数

一体密度を $\gamma$ とする数保存の自由 Fermi の大正準 Gibbs 状態で、観測する軌道部分空間への射影を $P$、個数を $N_P$ とする。この部分空間の直交基底で $e^{sN_P}=\prod_j[1+(e^s-1)n_j]$ である。各 $n_j$ は射影で互いに可換だからである。積を展開すると一組の軌道を選んだ項は、Wick の定理により $P\gamma P$ の対応する主小行列式になる。行列式 $\det[I+(e^s-1)P\gamma P]$ を行ごとの多重線形性で展開しても同じ主小行列式が現れる。従って両者は等しい。制限した行列の固有値を $\lambda_j$ とすると、分布は成功確率 $\lambda_j$ の独立な Bernoulli 変数の和と同じである。観測軌道の $n_j$ が独立とは限らないが、個数全体の分布にはこの形が使える。

:::example[非対角成分が分散を変える]
二軌道で $\gamma=\begin{pmatrix}1/2&1/4\\1/4&1/2\end{pmatrix}$ とする。固有値は $1/4,3/4$ なので $P(N_P=0)=P(N_P=2)=3/16$、$P(N_P=1)=5/8$ である。平均は $1$、分散は $3/8$ となる。元の軌道が独立に確率 $1/2$ で占有されるとすると分散は $1/2$ になるので、非対角の相関を落とすと違いが出る。
:::

### Wick の行列式と温度の符号

数を保存する Gaussian の熱状態で、$m$ 個の生成と $m$ 個の消滅を逆順に並べる。完全対の全てについて、一番左の生成 $a_1^\dagger$ が縮約する消滅を選ぶ。生成列と消滅列を揃えて書いた行列 $C$ の中で第 $l$ 列と結ぶと、対応する置換の符号は $(-1)^{1+l}$ である。残った対の集合は $C$ の第 $1$ 行と第 $l$ 列を除いた完全対に一致する。従って全体の和は Laplace 展開によって $\det C$ になる。各成分は $\langle a_i^\dagger a_{j_l}\rangle$ である。時間順序で別の並びにした場合は、元の場の並びの置換で全体の符号が変わる。Bose の同じ並びでは符号がないので、行列式と同じ置換和から符号だけを外した permanent が現れる。

:::example[Fermi と Bose で二つの同じ縮約を使う]
二モードの数の相関を考え、対角の占有の平均を $n_1,n_2$、非対角の熱期待値を $c$ とする。モードが異なる場合、完全対の二つを取ると $\langle n_1n_2\rangle=n_1n_2-|c|^2$ が Fermi、$n_1n_2+|c|^2$ が Bose である。例えば $n_1=n_2=1/2$、$|c|=1/4$ ならそれぞれ $3/16,5/16$ である。これはよく混ざった二軌道で同時占有を計算するとき、同じ一体密度の値を入れても統計により相関が反対向きに補正される例である。
:::

### 異常成分を含む二次状態の四点平均

二次 Hamiltonian に対生成や対消滅がある場合でも、正準 Bogoliubov 変換で自由な準粒子へ対角化できる有限模型では、元の場は準粒子の線形結合である。従って熱 Wick の定理を線形性で移せる。このとき $\langle a_i a_j\rangle$ と $\langle a_i^\dagger a_j^\dagger\rangle$ の縮約も非零であり、単に生成と消滅の行列式へまとめることはできない。例えば異なる二モードの $n_1n_2=a_1^\dagger a_1 a_2^\dagger a_2$ では
$$
\langle n_1 n_2\rangle=\langle n_1\rangle\langle n_2\rangle
-\langle a_1^\dagger a_2^\dagger\rangle\langle a_1 a_2\rangle
+\langle a_1^\dagger a_2\rangle\langle a_1 a_2^\dagger\rangle .
$$
四因子の二、三番目を交換して一対を作る項には負が付き、二回交換する最後の対には正が付く。異常平均の二つは随伴と反交換を使うと積が $-|\langle a_1a_2\rangle|^2$ になり、最初の負号と合わせて同時占有を増やす。通常の非対角の積は反対の符号に寄与する。

:::example[零温度の対状態を直接検査する]
$\ket\Psi=(\ket0+a_1^\dagger a_2^\dagger\ket0)/\sqrt2$ は二モードの対 Hamiltonian の Gaussian 基底である。両方の占有は $1/2$、通常の非対角は零、$|\langle a_1a_2\rangle|=1/2$ である。上の Wick 式は $\langle n_1n_2\rangle=1/4+1/4=1/2$ を与える。直接に占有部門を数えても、二電子部門の重みは $1/2$ である。通常の一体密度だけを使った行列式では $1/4$ となり、対の振幅が消えてしまう。
:::

### 有限温度の密度相関から自由な応答へつなぐ

数を保存する二次 Hamiltonian $K=\sum_r \xi_r c_r^\dagger c_r$ で、$B=\sum_{rs}B_{rs}c_r^\dagger c_s$ を Hermite とする。連結相関では二つの頂点内の縮約の積を引くので、$0<\tau<\beta$ で残るのは
$$
\langle B(\tau)B(0)\rangle_c
=\sum_{rs}|B_{rs}|^2 e^{-(\xi_s-\xi_r)\tau}f_r(1-f_s).
$$
最初の生成と最後の消滅を縮約し、残りの消滅と生成の縮約を取る。Kronecker デルタで添字は反転した二組に揃い、符号はこの順序では正である。$\xi_s-\xi_r\ne0$ なら Fourier 積分で $f_r(1-f_s)(1-e^{-\beta(\xi_s-\xi_r)})=f_r-f_s$ を得る。従って $-\langle T_\tau BB\rangle_c$ の $i\nu_m$ での係数は $(f_r-f_s)/(i\nu_m+\xi_r-\xi_s)$ を含む。等エネルギーの組では先の指数が定数なので、$-\beta f_r(1-f_s)\delta_{m0}$ が残る。これは Wick 計算でも、保存される部分の零モードを明示的に保持する方法である。

:::example[同じ場に対角と非対角の成分がある]
二準位 $\xi_1=0,\xi_2=\Delta>0$ で $B=n_1+c_1^\dagger c_2+c_2^\dagger c_1$ とする。$\beta\Delta=\log3$ では $f_1=1/2,f_2=1/4$ であり、非対角の寄与は
$\frac14[(i\nu_m-\Delta)^{-1}-(i\nu_m+\Delta)^{-1}]$ である。$n_1$ の部分は零モードだけの $-\beta/4$ となる。従って $i\nu_m=0$ での全静的な値は $-1/(2\Delta)-\beta/4$ である。非零点を遅延応答へ接続してから零周波数を取った値は $-1/(2\Delta)$ である。二種類の測定が違うのは、再平衡した占有分布を含めるかどうかによる。
:::

## 演習問題

::::exercise{level=basic}
逆温度 $\beta=4$ の最低正 Fermi 周波数と最低非零正 Bose 周波数を求めよ。

:::solution
境界条件から順に $\pi/4$ と $2\pi/4=\pi/2$ である。Bose の $\nu_0=0$ は「非零」の条件から除く。
:::
::::

::::exercise{level=basic}
$\beta\xi=\log 2$ の自由モードの Fermi 占有と Bose 占有を求めよ。

:::solution
熱分布への代入により $f(\xi)=1/(2+1)=1/3$、$n_B(\xi)=1/(2-1)=1$ である。
:::
::::

::::exercise{level=basic}
$\beta=4$ における $\beta^{-1}\sum_n(i\omega_n)^{-2}$ を求めよ。

:::solution
:eqref[eq-mb-matsubara-two-poles] で二つの極をともに $0$ に近づける。結果は $f'(0)=-\beta f(0)(1-f(0))=-4/4=-1$ である。
:::
::::

::::exercise{level=basic}
$E=2$、$\beta E=\log 3$ における $\beta^{-1}\sum_n(\omega_n^2+E^2)^{-1}$ を求めよ。

:::solution
:eqref[eq-mb-pair-sum] の $\tanh(\beta E/2)=(e^{\beta E}-1)/(e^{\beta E}+1)=1/2$ を使う。従って $(1/2)/(2E)=1/8$ である。
:::
::::

::::exercise{level=advanced}
有限次元の系で $[K,Q]=0$、$Q^\dagger=Q$ とする。連結した虚時間相関 $C(\tau)=\langle T_\tau Q(\tau)Q(0)\rangle-\langle Q\rangle^2$ と遅延交換子の違いを求めよ。

:::solution
可換性により $Q(\tau)=Q$ なので $C(\tau)=\operatorname{Var}(Q)$ である。Fourier 変換では $\beta\operatorname{Var}(Q)\delta_{m0}$ になる。遅延関数では同じ $Q$ どうしの交換子が零のため、全時刻で零である。
:::
::::

::::exercise{level=advanced}
二準位の熱分布 $f(a)=3/4$、$f(b)=1/4$、$b-a=2$ に対する Bose 外部周波数付きの二本の線の和を計算せよ。

:::solution
$z=a$ と $z=b-i\nu_m$ で留数を取る。後者の熱重みは周期性から $f(b)$ である。二つを足せば $(3/4-1/4)/(i\nu_m-2)=1/[2(i\nu_m-2)]$ となる。
:::
::::

::::exercise{level=advanced}
独立な二つの熱 Fermi モードの占有が $p_1=1/4,p_2=3/4$ のとき、$N=n_1+n_2$ の分散を Wick の定理で求めよ。

:::solution
異なる軌道の一体縮約が零なので $\langle n_1n_2\rangle=p_1p_2$ である。各 $n_i^2=n_i$ であるため $\operatorname{Var}(N)=p_1(1-p_1)+p_2(1-p_2)=3/8$ になる。
:::
::::

::::exercise{level=advanced}
半充填 Hubbard 原子で $\beta U/2=\log 7$ とする。二重占有と上スピンのスペクトルの二次モーメントを求めよ。

:::solution
このとき $x=7$ と $x^2y=1$ より $Z=16$、二重占有は $1/16$ である。二つのスペクトル線 $\pm U/2$ の重みはそれぞれ半分であるので、二次モーメントは $U^2/4$ である。
:::
::::

::::exercise{level=advanced}
一つの自由 Fermi 軌道について、$\beta^{-1}\sum_n(i\omega_n-a)^{-4}$ を占有 $f(a)$ の式で表し、$a=0$ での値を求めよ。

:::solution
分母を $a$ で三回微分すると $6(i\omega_n-a)^{-4}$ である。二因子以上の和は絶対収束するので、二重極の式をさらに二回微分して $f'''(a)/6$ となる。$f'=-\beta f(1-f)$ と $f''=\beta^2 f(1-f)(1-2f)$ を使うと
$$
f'''=-\beta^3 f(1-f)(1-6f+6f^2).
$$
$a=0$ では $f=1/2$ なので $f'''(0)/6=\beta^3/48$ である。同じ点で $(i\omega_n)^{-4}=\omega_n^{-4}$ であり、全ての項が正になることとも一致する。
:::
::::

::::exercise{level=advanced}
$K=h\sigma_z+\lambda\sigma_x$、$h>0$ の二準位で $\lambda=0$ における $\partial_\lambda^2\Omega$ を求め、虚時間相関の積分からも確かめよ。

:::solution
固有エネルギーは $\pm\sqrt{h^2+\lambda^2}$ であり、$\Omega=-\beta^{-1}\log[2\cosh(\beta\sqrt{h^2+\lambda^2})]$ である。最初の微分で $\lambda/\sqrt{h^2+\lambda^2}$ が出るため、二階微分は $\lambda=0$ で $-\tanh(\beta h)/h$ である。一方、$B=\sigma_x$ は二つの状態を結び、対角要素がないので $\langle B\rangle=0$ である。熱重みを $p_-=e^{\beta h}/Z$、$p_+=e^{-\beta h}/Z$ とすると正時間の相関は $p_- e^{-2h\tau}+p_+e^{2h\tau}$ である。積分は
$$
\int_0^\beta \langle B(\tau)B(0)\rangle\dd\tau
=\frac{p_-(1-e^{-2\beta h})+p_+(e^{2\beta h}-1)}{2h}
=\frac{p_--p_+}{h}.
$$
従ってその負は固有エネルギーからの結果と一致する。等時刻の $B^2=I$ だけから $-\beta$ とすると、量子遷移の時間依存を落とす。
:::
::::

## 確認テスト

::::quiz
Fermi の温度 Green 関数の境界条件はどれか。

- [ ] 周期的である。
- [x] 反周期的である。
- [ ] $\tau=\beta$ で必ず零になる。

:::explanation
熱トレースを循環させると正の時間差の積が順序を入れ替えた負の時間差の積になる。Fermi の時間順序の符号を戻すと $\mathcal G(\tau+\beta)=-\mathcal G(\tau)$ となる。
:::
::::

::::quiz{answer="0.3333333333" tol="1e-6"}
自由 Fermi 軌道で $\beta\xi=\log 2$ の場合、占有率を答えよ。

:::explanation
二状態の重みは $1$ と $1/2$ である。占有率は $(1/2)/(1+1/2)=1/3$ となる。
:::
::::

::::quiz{answer="0.125" tol="1e-6"}
$E=2$、$\beta E=\log 3$ のとき、$\beta^{-1}\sum_n(\omega_n^2+E^2)^{-1}$ を求めよ。

:::explanation
二極の留数和は $\tanh(\beta E/2)/(2E)$ である。分子が $(3-1)/(3+1)=1/2$、分母が $4$ になるため、結果は $1/8$ である。
:::
::::

::::quiz
有限温度について正しい記述を選べ。

- [x] 密度の相関には Bose 周波数を使う。
- [x] 遅延関数はスペクトル表示の上半平面からの境界値である。
- [ ] 任意の相互作用を含む熱状態で Wick 分解が成り立つ。

:::explanation
偶数個の Fermi 演算子からなる密度は周期的なので Bose 周波数を取る。遅延関数の正時間の Laplace 積分は上半平面で収束する。有限温度 Wick 分解の証明には二次 Hamiltonian の熱重みを使い、Hubbard 原子の二重占有が一般の Gibbs 状態での反例になる。
:::
::::

::::quiz{answer="0.1666666666667" tol="1e-6"}
$\beta=2\,\mathrm{meV^{-1}}$ とする。$E\to0$ の Fermi 和 $\beta^{-1}\sum_n(\omega_n^2+E^2)^{-2}$ を $\mathrm{meV^{-3}}$ で答えよ。

:::explanation
二線の和 $J(E)=\tanh(\beta E/2)/(2E)$ は $J(E)=\beta/4-\beta^3 E^2/48+O(E^4)$ である。求める和は $-J'(E)/(2E)$ なので、極限は $\beta^3/48$ である。$\beta=2\,\mathrm{meV^{-1}}$ では $8/48=1/6\,\mathrm{meV^{-3}}$ となる。最小 Fermi 周波数が非零なので、項ごとに $E=0$ を代入しても和は収束する。
:::
::::

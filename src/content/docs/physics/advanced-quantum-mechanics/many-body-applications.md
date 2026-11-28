---
title: 第二量子化の多体系への応用
description: 平均場と集団励起の近似を導出し、電子ガス、Bose 気体、磁性体に適用する。
notation:
  - symbol: 'k_F'
    meaning: Fermi 波数
  - symbol: 'E_k'
    meaning: Bogoliubov 準粒子エネルギー
  - symbol: 'S(\vb q,\omega)'
    meaning: 動的構造因子
  - symbol: 'S'
    meaning: 格子点当たりのスピン量子数
---

## この章の目標

- Slater 行列式の変分から Hartree–Fock 方程式を導く。
- 一様電子ガスの交換エネルギーを Fermi 球の重なりから求める。
- 希薄 Bose 気体と強磁性体の集団励起を二次 Hamiltonian で計算する。
- 励起の分散と密度応答から、散逸の閾値や散乱スペクトルを解釈する。

## Hartree–Fock 近似

相互作用するフェルミオンを一つの Slater 行列式で近似すると、反対称性を厳密に課したまま軌道を変分できる。正規直交する占有軌道を $\phi_1,\ldots,\phi_N$、一体演算子を $h$、粒子対の相互作用を $W(x,y)$ とする。スピンも $x$ に含める。

占有基底の状態で $\langle a_i^\dagger a_j^\dagger a_l a_k\rangle=n_i n_j(\delta_{ik}\delta_{jl}-\delta_{il}\delta_{jk})$ である。消滅した軌道を同じ組で戻す必要があり、二通りの順序の符号が逆になるからである。:eqref[eq-aqm-two-body] に代入して

$$
E=\sum_i\mel{i}{h}{i}+\frac12\sum_{ij}(W_{ij;ij}-W_{ij;ji})
\label{eq-aqm-hf-energy}
$$

となる。第一の二体項を直接項、第二を :term[交換項]{en="exchange term" yomi="こうかんこう"} と呼ぶ。$i=j$ の項は相殺し、一粒子が自分自身と相互作用する項を除く。

直交条件を Hermite 行列 $\lambda$ で課し、$E-\sum_{ij}\lambda_{ji}(\braket{\phi_i}{\phi_j}-\delta_{ij})$ を各 $\overline{\phi_i(x)}$ で変分する。二体項では変分する軌道が第 1 因子または第 2 因子に出る二寄与が等しく、$1/2$ が消える。すると

$$
\begin{aligned}
F\phi_i&=\sum_j\phi_j\lambda_{ji},\qquad F=h+\sum_j(J_j-K_j),\\
(J_j\phi)(x)&=\phi(x)\int W(x,y)|\phi_j(y)|^2\dd y,\\
(K_j\phi)(x)&=\phi_j(x)\int W(x,y)\overline{\phi_j(y)}\phi(y)\dd y .
\end{aligned}
\label{eq-aqm-hf}
$$

この自己無撞着方程式を :term[Hartree–Fock 方程式]{en="Hartree–Fock equation" yomi="はーとりーふぉっくほうていしき"} と呼ぶ。右辺の $\lambda$ を占有軌道のユニタリ変換で対角化できる。その変換は Slater 行列式に行列式の位相を掛けるだけで、占有射影 $P$ と $F$ を保つ。解を得た後、軌道エネルギーを足すと二体エネルギーを二回含むので $E=\sum_i\epsilon_i-\tfrac12\sum_{ij}(W_{ij;ij}-W_{ij;ji})$ と戻す。各 Slater 行列式は規格化した試行状態なので、その集合でのエネルギーの下限は厳密な基底エネルギー以上である。自己無撞着解の中には極小以外の停留点もある。

:::example[同じ空間軌道で反対スピン]
空間軌道 $\phi$ に上向きと下向きを一個ずつ入れると、スピン非依存相互作用の交換項は $\braket{\uparrow}{\downarrow}=0$ で消える。従って二体エネルギーは空間座標だけの直接積分 $\int |\phi(\vb r)|^2W(\vb r,\vb r')|\phi(\vb r')|^2\dd^3r\dd^3r'$ の一回分である。
:::

:::example[自己項の相殺]
$N=1$ なら $J_1\phi_1=K_1\phi_1$ なので占有軌道の式は $h\phi_1=\epsilon_1\phi_1$ に戻る。仮に $W_{11;11}=3\,\mathrm{eV}$ でも二体寄与は零である。
:::

### 占有射影での停留条件と一粒子エネルギー

有限モードで、正規直交する占有軌道の射影 $P=\sum_{i\in{\rm occ}}\ketbra{i}{i}$ を使えば、Hartree–Fock のエネルギーは $P$ の線形項と二次項である。交換対称性で二次項を微分すると二つの寄与が等しく、$\delta E=\tr(F\,\delta P)$ となる。$P(\lambda)=e^{\lambda X}Pe^{-\lambda X}$、$X^\dagger=-X$ とすれば $\delta P=[X,P]$ である。巡回性から $\delta E=\tr([P,F]X)$ となる。全ての $X$ で零なら特に $X=[P,F]^\dagger$ を選んでノルム二乗が零になり、$[P,F]=0$ を得る。したがって占有・非占有を結ぶ Fock 行列要素は零である。逆も同じ式で従う。この条件は停留点を特徴づけるので、安定性を判定するには二次変分または他の停留解とのエネルギー比較が要る。

:::example[凍結軌道で粒子を取り除く]
占有 $i$ を除き、残った軌道を変えない。エネルギーの変化は $-h_{ii}-\sum_j(W_{ij;ij}-W_{ij;ji})=-F_{ii}$ である。実際、二体和では $i$ が先頭または末尾に来る二組があり、自己項は零なので $1/2$ と相殺する。$h_{ii}=-7\,\mathrm{eV}$、対の差の和が $2\,\mathrm{eV}$ なら除去コストは $5\,\mathrm{eV}$ である。除去後の Slater 行列を同じ残数で最適化すると、変分原理でそのエネルギーは凍結値以下になる。従って同じ初期状態からの緩和後の除去コストは $5\,\mathrm{eV}$ 以下になる。
:::

### 一つの占有軌道を回したときの一次と二次

直交軌道一組 $i\in{\rm occ},a\notin{\rm occ}$ だけを動かし、占有を $\cos\theta\ket i+e^{i\phi}\sin\theta\ket a$ に替える。他の占有との正規直交性は保たれる。射影の一次変化は $\theta(e^{i\phi}\ket a\bra i+e^{-i\phi}\ket i\bra a)$ だから、全エネルギーの一次は $2\theta\Re(e^{i\phi}F_{ia})$ である。位相を任意に選んでも零になるのは $F_{ia}=0$ の場合である。この静止性は :term[Brillouin の定理]{en="Brillouin theorem" yomi="ぶりるあんのていり"} の一電子励起の形に対応する。実際、行列状態を微分すると $\ket{\Phi_i^a}$ に比例するので、同じエネルギーの一次は $2\Re(e^{i\phi}\mel{\Phi}{H}{\Phi_i^a})$ である。軌道の順序と励起状態の位相を揃えれば両者が一致する。

この特定の一つの軌道だけを動かす場合、エネルギーには、その軌道の一体項と動かない占有との対だけが含まれる。従って $\theta$ 依存部分を $2\times2$ の Hermite 行列 $B^{(i)}$ の期待値にできる。$B^{(i)}_{ab} = h_{ab}+\sum_{j\ne i}(W_{aj;bj}-W_{aj;jb})$ とし、その $i,a$ の小ブロックを取る。自己項を除いたため、$B_{ii}^{(i)}=F_{ii}$、$B_{ia}^{(i)}=F_{ia}$ だが $B_{aa}^{(i)}$ は $F_{aa}$ と一般に異なる。停留した $\theta=0$ での二次係数は $B_{aa}^{(i)}-B_{ii}^{(i)}$ である。全ての単独回転でこの係数が正でも、複数を同時に動かす変分には交差項があり、安定の十分条件にはならない。

:::example[自己無撞着な停留解が最低でない場合]
二つの占有のうち動かす $i$ に対し $B_{ii}^{(i)}=-2\,\mathrm{eV}$、候補 $a$ に対し $B_{aa}^{(i)}=-3\,\mathrm{eV}$、$B_{ia}^{(i)}=0$ とする。$\theta$ を小さくした変化は $-\sin^2\theta\,\mathrm{eV}$ である。$\theta=\pi/6$ では $-1/4\,\mathrm{eV}$ になる。固定した $P$ から作った $F$ と $P$ が可換というだけでは、変分の最小を保証できない。$\phi$ によらず負の曲率なので、占有の条件を固定し直す必要がある。
:::

### 平均場から外れた二粒子の励起

占有空間と非占有空間をそれぞれ対角化して $F$ の正準軌道を取る。$H_0$ の Slater 固有値を軌道エネルギー $\epsilon_i$ の和で与える。$H-H_0$ の二体部分は占有二個を取り除いて非占有二個を作ることができる。一体励起は上の静止条件で零で、三個以上の励起は二体演算子一回では到達できない。二次摂動の公式をこの直交基底に挿入すると、非縮退で分母が零でない場合の補正は
$$
E^{(2)}=
-\sum_{i<j\in{\rm occ}}\sum_{a<b\notin{\rm occ}}
\frac{|W_{ab;ij}-W_{ab;ji}|^2}{\epsilon_a+\epsilon_b-\epsilon_i-\epsilon_j}.
$$
反対称化された行列要素は、二体項の作る $ij$ と $ji$ の二順序の差である。順序条件を外して全添字を和にする場合は各対の入れ替えで四回数えるので $1/4$ が掛かる。占有が低い軌道に揃えば分母が正で、この補正はエネルギーを下げる。ただし分母が小さいと混合が強く、選んだ一つの Slater 行列のまわりの摂動の適用を見直す。

:::example[二つの配置で二次の尺度を検査する]
基準のエネルギーを零、励起の差を $\Delta=5\,\mathrm{eV}$、反対称化行列要素を $W=1\,\mathrm{eV}$ とする。二次の式は $-W^2/\Delta=-0.2\,\mathrm{eV}$ である。二配置だけなら厳密な低い根は $(\Delta-\sqrt{\Delta^2+4W^2})/2\simeq-0.193\,\mathrm{eV}$ になる。$W$ を二倍にする場合、摂動の補正は四倍と予測するが、厳密根は約 $-0.702\,\mathrm{eV}$ である。$W/\Delta$ を評価することで、補正値だけから精度を判断する誤りを避けられる。
:::

### 二サイト一重項と単一 Slater 行列を比較する

二サイトの Hubbard 模型で二電子、$t>0,U\ge0$ を取る。スピンが上と下の一電子ずつを同じ空間軌道 $\sqrt x\,\ket1+\sqrt{1-x}\ket2$ に置く。位相を揃えると移動エネルギーは $-4t\sqrt{x(1-x)}$、同じサイトに二つがいる確率は $x^2+(1-x)^2$ なので相互作用は $U[x^2+(1-x)^2]$ である。$0\le x\le1$ で $x(1-x)$ が最大なのは $x=1/2$ であり、両項がその点で最小となる。従って同じ空間軌道を占有するよう制限した Hartree–Fock の値は $E_{\rm RHF}=-2t+U/2$ である。

別々の軌道を許し、上に $(\cos\theta,\sin\theta)$、下に $(\sin\theta,\cos\theta)$ を用いる。各移動エネルギーは $-t\sin2\theta$、同じサイトの確率の和は $2\sin^2\theta\cos^2\theta$ である。一般に各スピンの軌道を実かつ非負に取ると移動エネルギーが最小になる。両軌道の角を $\theta_\uparrow,\theta_\downarrow$ として $s_\sigma=\sin2\theta_\sigma$ と置き、サイト占有の偏り $\cos 2\theta_\sigma$ を逆符号に選べば、相互作用エネルギーは $U[1-\sqrt{(1-s_\uparrow^2)(1-s_\downarrow^2)}]/2$ になる。$s=(s_\uparrow+s_\downarrow)/2$ を固定すると、相加・相乗平均の不等式と二乗和から $\sqrt{(1-s_\uparrow^2)(1-s_\downarrow^2)}\le 1-s^2$ である。等号は同じ $s_\sigma=s$ を選ぶ上の軌道対で達成される。このため $s=\sin2\theta\in[0,1]$ としたエネルギーの最小値が、上・下各一粒子の共線な Slater 行列式の範囲での最小値になる。
$$
E_{\rm UHF}(s)=-2ts+\frac U2 s^2
$$
になる。$U\le2t$ では最小は $s=1$、$U>2t$ なら $s=2t/U$ に移って $-2t^2/U$ になる。後者は一サイトに上、他方に下を置く傾向を示すが、その Slater 行列は一般に全スピンの固有状態ではない。

:::example[相関の欠落が定量的に見える強さ]
$U=4t$ なら制限した値は零、別軌道の最適値は $-t/2$ である。厳密な最低値は二重占有と一重項を対角化した $(U-\sqrt{U^2+16t^2})/2=(2-2\sqrt2)t\simeq-0.828t$ である。二つの Slater 行列の範囲の拡大でエネルギーが下がっても、厳密な一重項の相関までは再現できない。$U/t\to\infty$ で別軌道の先頭 $-2t^2/U$ を厳密な $-4t^2/U$ と比較しても係数 2 の差が残る。
:::

## 一様電子ガスと交換エネルギー

一様な正電荷背景と体積 $\mathcal V$ の周期箱に電子を入れる。SI 単位で Coulomb 相互作用の非零 Fourier 成分は $\widetilde W(q)=e^2/(\epsilon_0 q^2)$ である。原点の成分は背景と打ち消し合わせる。均一な無偏極状態を仮定し、各スピンで $|\vb k|\le k_F$ を占有する。箱の和を $\mathcal V\int\dd^3k/(2\pi)^3$ に替えると密度 $n=k_F^3/(3\pi^2)$、運動エネルギーは一粒子当たり $3\hbar^2k_F^2/(10m)$ である。後者の平均は $\int_0^{k_F}k^4\dd k/\int_0^{k_F}k^2\dd k=3k_F^2/5$ から得る。

直接項は $q=0$ にだけ寄与し、交換項は同じスピンの $\vb k-\vb k'$ を含む。二スピンの和が $1/2$ と相殺するので

$$
E_x/\mathcal V
=-\frac{e^2}{\epsilon_0(2\pi)^6}
 \int_{|\vb k|,|\vb k'|\le k_F}\frac{\dd^3k\dd^3k'}{|\vb k-\vb k'|^2}.
$$

$\vb q=\vb k-\vb k'$ を固定した積分は半径 $K=k_F$ の二球の重なり体積 $V_{\cap}(q)$ になる。これは高さ $h=K-q/2$ の球冠二個であり、$2\int_{K-h}^K\pi(K^2-z^2)\dd z=\pi(4K+q)(2K-q)^2/12$ である。$q^2$ は角積分の測度と打ち消し合う。残る $\int_0^{2K}V_\cap(q)\dd q=\pi K^4$ から

$$
\frac{E_x}{N}=-\frac{3}{4\pi}\frac{e^2 k_F}{4\pi\epsilon_0}.
\label{eq-aqm-jellium-exchange}
$$

負の符号は同一スピン電子が近くに来る確率が交換により減ることに対応する。Hartree–Fock は密度の変動に伴う相関エネルギーを含まない。

:::example[密度を変える]
密度を 8 倍にすると $k_F$ が 2 倍、運動エネルギーは 4 倍、交換の絶対値は 2 倍である。従ってこの二項の比では高密度ほど運動項が相対的に大きい。
:::

:::example[完全スピン偏極]
全密度を保って一種のスピンだけを使うと $K=2^{1/3}k_F$ になる。交換密度はスピン種当たり $-K^4$ に比例するため、無偏極の場合との絶対値の比は $2^{4/3}/2=2^{1/3}$ であり、運動項も $2^{2/3}$ 倍になる。交換だけを見て基底状態の偏極を決めることはできない。
:::

### スピン分極による交換と運動の競合

実空間密度 $n$ を固定し、$n_\uparrow=n(1+\zeta)/2$、$n_\downarrow=n(1-\zeta)/2$ とする。各 Fermi 球の半径は非偏極の半径 $k_F$ に $(1\pm\zeta)^{1/3}$ を掛けたものである。交換は同じスピンの球だけを使うので、体積当たりの交換エネルギーは半径四乗の和に比例する。運動エネルギーは $\int_{|\vb k|\le K}k^2\dd^3k\propto K^5$ である。従って非偏極を基準にした倍率はそれぞれ
$$
X(\zeta)=\frac{(1+\zeta)^{4/3}+(1-\zeta)^{4/3}}2,\qquad
T(\zeta)=\frac{(1+\zeta)^{5/3}+(1-\zeta)^{5/3}}2 .
$$
両方増加するが、交換のエネルギーは負である。二項展開すると $X=1+2\zeta^2/9+\order{\zeta^4}$、$T=1+5\zeta^2/9+\order{\zeta^4}$ なので、非偏極での運動項 $t$ と交換項 $e_x<0$ に対して局所の二次係数は $(5t+2e_x)/9$ になる。この判定は一様 Slater 行列の中のもので、密度変調や相関の自由度を含まない。

:::example[完全分極と非偏極を二つの項で比べる]
$\zeta=1$ で $X=2^{1/3}$、$T=2^{2/3}$ である。$e_x=-2t$ なら完全分極との差は $t[2^{2/3}-1-2(2^{1/3}-1)]\simeq0.0676t$ で、非偏極の方が低い。$e_x=-3t$ なら差は約 $-0.192t$ となる。この二端点の比較だけから、途中の分極も含む大域的な最小値は求まらない。小分極の局所安定性は上の二次係数で別に判定する。
:::

### 交換穴を座標から求める

スピン $\sigma$ の占有射影の核を $\gamma_\sigma(\vb r)=\int_{|\vb k|<K_\sigma}e^{i\vb k\cdot\vb r}\dd^3k/(2\pi)^3$ とする。球の角積分 $4\pi\sin kr/(kr)$ を用い、$\int_0^K k\sin(kr)\dd k$ を部分積分すると
$$
\gamma_\sigma(r)=n_\sigma \frac{3 j_1(K_\sigma r)}{K_\sigma r},\qquad
n_\sigma=\frac{K_\sigma^3}{6\pi^2}.
$$
Slater 状態の二体縮約を使うと、同じスピンの二点の対相関を一様密度の積で割った値は $g_{\sigma\sigma}(r)=1-|\gamma_\sigma(r)|^2/n_\sigma^2$ である。逆スピンでは交換項が零なので $g_{\uparrow\downarrow}=1$ である。$\gamma_\sigma$ は射影の核なので $\int|\gamma_\sigma(\vb r)|^2\dd^3r=n_\sigma$ である。これは Fourier の Parseval 等式で占有の指示関数を二乗しても同じことからも分かる。従って固定した一粒子の周りの交換だけによる密度の不足 $n_\sigma[g_{\sigma\sigma}(r)-1]$ は空間積分で $-1$ になる。

:::example[最近接点の抑制を波数で比べる]
$z=K_\sigma r$ が小さいとき $3j_1(z)/z=1-z^2/10+\order{z^4}$ であり、同スピンの $g$ は $z^2/5+\order{z^4}$ である。例えば $z=1/2$ の先頭は $1/20$ であり、逆スピンの値 1 より小さい。この低下は反対称性だけで得られ、Coulomb の動的な相関を入れていない。一粒子当たりの交換エネルギーは、この不足を $e^2/(4\pi\epsilon_0r)$ で重みづけて積分し、対の重複を除く $1/2$ を掛けた値である。スピン分極がある場合は、さらに参照粒子のスピンを $n_\sigma/n$ の重みで平均する。
:::

### Fock の自己エネルギーを球の積分で得る

一個の同スピンの波数 $k$ に対する交換項は、Coulomb Fourier を $e^2/(\epsilon_0 q^2)$ としたとき
$$
\Sigma_x(k)=-\frac{e^2}{\epsilon_0}\int_{p<K}
\frac{\dd^3p}{(2\pi)^3}\frac1{k^2+p^2-2kp\cos\theta}.
$$
$x=\cos\theta$ の積分で $(kp)^{-1}\log|(k+p)/(k-p)|$ を得る。残る $\int_0^K p\log|(k+p)/(k-p)|\dd p$ を部分積分すると
$$
kK+\frac{K^2-k^2}{2}\log\left|\frac{k+K}{k-K}\right|.
$$
例えば微分すれば上限 $K$ での被積分関数となり、$K=0$ で零であるから確認できる。$p=k$ の対数は可積分であり、部分積分は左右を同じ切断で合わせて行える。従って
$$
\Sigma_x(k)=-\frac{e^2}{4\pi^2\epsilon_0}
\left[K+\frac{K^2-k^2}{2k}\log\left|\frac{k+K}{k-K}\right|\right]
$$
である。$k=0$ では角積分に戻って $-e^2 K/(2\pi^2\epsilon_0)$、$k=K$ では対数の係数が零になりその半分になる。Fermi 面で値は有限だが傾きが対数的に発散する。したがって裸の Coulomb 交換だけで電子ガスの励起速度を精密に読むには限界があり、遮蔽や相関を含める動機になる。

:::example[Fermi 球の中心と境界の交換を区別する]
$k=K/2$ では括弧が $K(1+(3/4)\log3)$ なので約 $1.824K$ である。中心 $2K$ と境界 $K$ の中間だが、単純な線形補間の $1.5K$ ではない。この関数は各軌道のエネルギーへの交換項である。多体の全交換は占有軌道に足して半分にする。同じ対が互いの自己エネルギーに一度ずつ現れるためである。
:::

### Fermi 球の全エネルギーから圧力を取る

無偏極で一様な電子ガスについて、Hartree 項は中和背景で相殺したとする。運動の一粒子当たりの平均は $(3/5)E_F$ である。波数積分の比 $\int_0^{k_F} k^4\dd k/\int_0^{k_F} k^2\dd k=3k_F^2/5$ による。交換を $-C_x k_F$ と書くと、固定粒子数でのエネルギーは $E=N(C_t n^{2/3}-C_e n^{1/3})$、$n=N/\mathcal V$、$C_t,C_e>0$ である。従って
$$
P=-\pdv E{\mathcal V}=n^2\dv{(E/N)}n
=n\left(\frac23e_t+\frac13e_x\right),
\quad
\mu=\pdv EN=\frac53e_t+\frac43e_x.
$$
最後は体積固定で微分した。中和背景を一緒に密度変化させる模型の式である。圧縮率の逆は $n\partial_nP=n[10e_t/9+4e_x/9]$ になる。この近似の値が負なら、仮定した一様変分の範囲と相関の補正を再検討する必要がある。

:::example[密度を変えたときに符号が変わる量]
ある密度で $e_t=3\,\mathrm{eV}$、$e_x=-5\,\mathrm{eV}$ と仮定する。このとき $P/n=1/3\,\mathrm{eV}$、$\mu=-5/3\,\mathrm{eV}$、$n^{-1}\partial_{\log n}P$ は $10/9\,\mathrm{eV}$ である。密度を八分の一にすると $e_t$ は $1/4$、$e_x$ は $1/2$ になり、$P/n=1/2-5/6=-1/3\,\mathrm{eV}$ に替わる。圧力や圧縮率の符号は、交換の絶対値が大きいことだけでなく密度の冪の違いを反映する。実際の低密度の基底状態はこの二項だけから決定しない。
:::

## 弱く相互作用する Bose 気体と Bogoliubov 理論

零温度の一様な三次元 Bose 気体を考える。短距離の斥力の低エネルギー結合を $g=4\pi\hbar^2 a_s/m>0$ とする。二体の換算質量 $m/2$ で :eqref[eq-aqm-first-born] の接触頂点を $f(0)=-a_s$ に合わせるとこの値になる。ここで $g$ は低エネルギーの有効結合である。紫外まで裸の定数として積分する場合には二体の散乱長に合わせる再正規化が必要である。

密度 $n=N/\mathcal V$、$n a_s^3\ll1$ を仮定する。凝縮モードとの粒子の移動の行列要素を $\sqrt{N_0}$ の次数で評価し、非零モードについて二次まで残す。粒子数 $N$ を固定する展開では $N_0=N-\sum_{\vb k\ne0}N_{\vb k}$ として凝縮エネルギー $gN_0^2/(2\mathcal V)$ の一次変化 $-gn\sum N_{\vb k}$ を含める。他方、相互作用の一個ずつの非零生成・消滅の位置は四通りで $2gn\sum N_{\vb k}$ となる。二個生成または二個消滅では波数保存により $\pm\vb k$ の組だけが残る。固定 $N$ での対生成は $a_{\vb k}^\dagger a_{-\vb k}^\dagger a_0^2$ であり、凝縮から二粒子を取り除く。非零モードの占有列を固定 $N$ の基底に同一視し、残る凝縮占有が十分大きいとき $\sqrt{N_0(N_0-1)}\simeq N$ とする。以下の $a_{\vb k}$ はこの二次近似における励起の占有を変える演算子である。固定 $N$ の表示では、その作用には凝縮との粒子の移し替えが含まれる。厳密には励起の全占有は $N$ 以下だが、枯渇が小さい領域でこの境界を外し、ボソンの正準交換関係を使う。よって

$$
H_2=\frac{gn^2\mathcal V}{2}+\sum_{\vb k\ne0}
 \left[(\epsilon_k+gn)a_{\vb k}^\dagger a_{\vb k}
 +\frac{gn}{2}(a_{\vb k}^\dagger a_{-\vb k}^\dagger+a_{\vb k}a_{-\vb k})\right],
\quad\epsilon_k=\frac{\hbar^2 k^2}{2m}.
$$

$\vb k,-\vb k$ ごとに $a_{\vb k}=u_k b_{\vb k}-v_k b_{-\vb k}^\dagger$ と置き、$u=\cosh r$、$v=\sinh r$ を用いる。異常項が消える条件は $\tanh(2r)=gn/(\epsilon_k+gn)$ で、:eqref[eq-aqm-bogoliubov-one] と同じ二次式から

$$
E_k=\sqrt{\epsilon_k(\epsilon_k+2gn)},\qquad
v_k^2=\frac12\left(\frac{\epsilon_k+gn}{E_k}-1\right).
\label{eq-aqm-bose-dispersion}
$$

$k\xi\ll1$（$\xi=\hbar/\sqrt{2mgn}$）で $E_k\simeq\hbar k c_s$、$c_s=\sqrt{gn/m}$ である。有効結合が使える範囲で $\epsilon_k\gg gn$ なら $E_k\simeq\epsilon_k+gn$ である。$c_s$ は音速であり、光速 $c$ と区別する。

### 凝縮から外れる粒子

準粒子真空でも $\langle a_k^\dagger a_k\rangle=v_k^2$ である。$x=k\hbar/\sqrt{2mgn}$ と置くと

$$
\frac{N-N_0}{\mathcal V}
=\frac{(2mgn)^{3/2}}{4\pi^2\hbar^3}
\int_0^\infty\left[\frac{x(x^2+1)}{\sqrt{x^2+2}}-x^2\right]\dd x .
$$

第一項の原始関数は $(x^2+2)^{3/2}/3-\sqrt{x^2+2}$、第二項は $-x^3/3$ で、その和は無限遠で零に近づく。原点との差から積分は $\sqrt2/3$ になる。従って枯渇割合は $8\sqrt{n a_s^3}/(3\sqrt\pi)$ である。

:::example[枯渇の大きさ]
$n a_s^3=10^{-6}$ なら枯渇割合は約 $1.50\times10^{-3}$ である。この小ささは非零モードの高次項を落とした近似の自己整合性を支える。
:::

:::example[運動項と相互作用の境界]
$\epsilon_k=gn$ の点で $E_k=\sqrt3 gn$、$v_k^2=1/\sqrt3-1/2$ となる。この波数の逆数 $\xi=\hbar/\sqrt{2mgn}$ を :term[回復長]{en="healing length" yomi="かいふくちょう"} と呼ぶ。実空間で凝縮密度が変化する際に、勾配の運動エネルギーと相互作用を比較する尺度でもある。
:::

### 基底エネルギーの最初の補正

二次 Hamiltonian の対角化で出る定数は $\tfrac12\sum_{\vb k\ne0}(E_k-\epsilon_k-gn)$ である。接触結合を全波数で一定にすると、$E_k=\epsilon_k+gn-(gn)^2/(2\epsilon_k)+\cdots$ の最後の項の和が三次元で線形発散する。裸の結合 $g_0$ を上限 $\Lambda$ 付きで使うと、真空の零エネルギーの二体散乱は接触頂点の級数で $g=g_0-g_0^2\int_{k<\Lambda}\dd^3k/[(2\pi)^3 2\epsilon_k]+\cdots$ になる。中間二体状態のエネルギーは $2\epsilon_k$ であり、分母は入射零エネルギーからそれを引いた値である。従って同じ摂動次数で $g_0=g+g^2\int\dd^3k/[(2\pi)^3 2\epsilon_k]+\cdots$ を凝縮エネルギーに入れると、発散を相殺する項が出る。この置換は固定した cutoff で次数をそろえてから、有限になった結果で cutoff を外す。

$\mu_0=gn$ とし、$x=\hbar k/\sqrt{2m\mu_0}$ を使うと

$$
\begin{aligned}
\frac{\Delta E}{\mathcal V}
&=\frac12\int\frac{\dd^3k}{(2\pi)^3}
\left[E_k-\epsilon_k-\mu_0+\frac{\mu_0^2}{2\epsilon_k}\right]\\
&=\frac{(2m\mu_0)^{3/2}\mu_0}{4\pi^2\hbar^3}
 \int_0^\infty\left[x^3\sqrt{x^2+2}-x^4-x^2+\frac12\right]\dd x .
\end{aligned}
$$

原始関数は $(x^2+2)^{5/2}/5-2(x^2+2)^{3/2}/3-x^5/5-x^3/3+x/2$ である。無限遠では各正冪が相殺して零に近づき、原点で $-8\sqrt2/15$ となる。従って

$$
\frac{E}{\mathcal V}
=\frac{gn^2}{2}\left[1+\frac{128}{15\sqrt\pi}\sqrt{n a_s^3}+\cdots\right].
$$

この $\sqrt{n a_s^3}$ の項を :term[Lee–Huang–Yang 補正]{en="Lee–Huang–Yang correction" yomi="りーほあんやんほせい"} という。例えば $n a_s^3=10^{-6}$ では平均場エネルギーから約 $0.48\%$ 増える。同じ密度の組でも有効距離や高次の多体相互作用の効果を比較するには、ここで残した次数より先を評価する必要がある。

### 音速、回復長、圧縮率

$g>0$ の凝縮体で、全粒子が規格化した同じ軌道 $\phi$ を占有する積状態を試行状態に取る。一体項は $N$ 個、接触相互作用は $N(N-1)/2$ 対を数えるので、$\psi=\sqrt N\phi$ を使って先頭のエネルギーは $\int[\hbar^2|\grad\psi|^2/(2m)+g|\psi|^4/2]\dd^3r$ である。粒子数の制約 $\int|\psi|^2\dd^3r=N$ を $\mu$ で課し、$\bar\psi$ で変分すると
$$
[-\hbar^2\nabla^2/(2m)+g|\psi|^2-\mu]\psi=0
$$
となる。これが静的な :term[Gross–Pitaevskii 方程式]{en="Gross–Pitaevskii equation" yomi="ぐろすぴたえふすきーほうていしき"} である。時間依存の式は $\int\dd t\,[\tfrac{i\hbar}{2}\int(\bar\psi\,\partial_t\psi-\psi\,\partial_t\bar\psi)\dd^3r-E[\psi]]$ を変分して得られ、$\mu\psi$ の位置に $i\hbar\partial_t\psi$ が入る。周期境界または無限遠での十分な減衰を仮定し、部分積分の境界項を消した。一様解 $\psi=\sqrt n$ から $\mu=gn$ である。運動エネルギーの尺度 $\hbar^2/(2m\xi^2)$ と $gn$ を等しくする長さ $\xi=\hbar/\sqrt{2mgn}$ を :term[回復長]{en="healing length" yomi="かいふくちょう"} と呼ぶ。Bogoliubov の分散は $E_k=\hbar c_s k\sqrt{1+(k\xi)^2/2}$、$c_s^2=gn/m$ になる。長波長の流体方程式でも、連続式 $\partial_t\delta n=-n\grad\cdot\vb v$ と $m\partial_t\vb v=-\grad\delta\mu=-g\grad\delta n$ を組み合わせると $\partial_t^2\delta n=(gn/m)\nabla^2\delta n$ を得る。圧縮で化学ポテンシャルが増えることが音の復元力になる。

:::example[回復長と同じ尺度の波を作る]
$k\xi=1$ では $\epsilon_k=gn$、$E_k=\sqrt3\,gn$ である。一方、音波の直線を外挿すると $\hbar c_sk=\sqrt2\,gn$ なので、厳密な二次理論の値は $\sqrt{3/2}$ 倍になる。粒子と空孔の係数は $u_k^2+v_k^2=(\epsilon_k+gn)/E_k=2/\sqrt3$ であり、$u_k^2-v_k^2=1$ と合わせると $v_k^2=(2/\sqrt3-1)/2\simeq0.0774$ となる。低波数での音波と高波数での自由粒子の中間にある。
:::

### 基底エネルギーの補正を音速に反映する

一様な流体では化学ポテンシャルは $\mu=\partial(E/\mathcal V)/\partial n$ である。零温度の連続式と運動方程式を使った音速は $c_s^2=(n/m)\partial_n\mu$ となる。$g,a_s$ を固定し、先に求めた基底エネルギー密度の $n^2$ と $n^{5/2}$ の項を微分すると
$$
\begin{aligned}
\mu&=gn\left[1+\frac{32}{3\sqrt\pi}\sqrt{na_s^3}+\cdots\right],\\
c_s^2&=\frac{gn}{m}\left[1+\frac{16}{\sqrt\pi}\sqrt{na_s^3}+\cdots\right].
\end{aligned}
$$
二次の分散を導いたときの音速は、この先頭の項である。例えば $na_s^3=10^{-6}$ では $c_s^2$ の増分は先頭の値の約 $0.90\%$ であり、音速そのものの増分は展開の同じ次数で約 $0.45\%$ になる。

### 凝縮体の位相と密度の線形方程式

時間依存の平均場で $\psi=\sqrt{n+\delta n}\,e^{i(\theta-gn t/\hbar)}$ と置く。$\psi$ の Laplacian は振幅の Laplacian と位相の勾配、両者の積の項を持つ。一次では $i\hbar \partial_t\psi$ の虚部から $\partial_t\delta n=-(\hbar n/m)\nabla^2\theta$、実部から $\hbar\partial_t\theta=-g\delta n+\hbar^2\nabla^2\delta n/(4mn)$ である。後者の勾配補正は $\sqrt{n+\delta n}=\sqrt n+\delta n/(2\sqrt n)+\cdots$ から出る。二式を組み合わせれば
$$
\partial_t^2\delta n=\frac{gn}{m}\nabla^2\delta n-\frac{\hbar^2}{4m^2}\nabla^4\delta n.
$$
平面波の代入は $\omega^2=c_s^2k^2+\hbar^2k^4/(4m^2)$ を与え、二次の Bogoliubov 分散と一致する。この線形方程式を用いるには、密度と速度の変動が二次以上の項を無視できるほど小さいことが必要である。

:::example[位相だけの初期変調が密度へ変換する]
実の一つの Fourier 成分で $\delta n(0)=0,\theta(0)=\theta_0$ とする。$\partial_t\delta n(0)=(\hbar n/m)k^2\theta_0$ だから $\delta n(t)=(\hbar n k^2\theta_0/(m\omega_k))\sin\omega_k t$ である。$k\xi=1$ の場合 $\omega_k=\sqrt3\,gn/\hbar$ で、振幅を $n$ で割ると $2\theta_0/\sqrt3$ になる。$\theta_0=0.03$ rad では約 $0.0346$ であり、一次の仮定に合う。この応答で密度が初期に変わっていなくても、位相の勾配による速度が密度の変化を起こす。
:::

### 密度を回復する静的な界面

剛体壁 $x=0$ の右側で一様凝縮体が遠方に密度 $n$ を持ち、壁で $\psi(0)=0$ とする。無次元な $f=\psi/\sqrt n$ を実に選ぶと $-\xi^2 f''+(f^2-1)f=0$、$f(\infty)=1$ である。$f'$ を掛けて積分すると $\xi^2(f')^2=(1-f^2)^2/2$ となる。単調増加する解では分離して $\operatorname{artanh}f=x/(\sqrt2\,\xi)$、従って $f=\tanh[x/(\sqrt2\,\xi)]$ になる。回復長の定義で $\sqrt2$ の置き場所が異なる場合は導出に使った $\xi$ を比較する。この解は一様な壁を横に無限とし、境界に沿った速度と変化を無視した静的模型である。

:::example[密度が半分に達する距離]
$|\psi|^2/n=1/2$ は $f=1/\sqrt2$ の点なので、$x=\sqrt2\,\xi\operatorname{artanh}(1/\sqrt2)\simeq1.246\,\xi$ である。振幅が半分になる位置 $x=\sqrt2\,\xi\operatorname{artanh}(1/2)$ とは違う。壁の近くの $x\ll\xi$ では密度は $n x^2/(2\xi^2)$、遠方の不足 $1-f^2$ は $4e^{-\sqrt2 x/\xi}$ に近づく。どの位置で近似が使えるかを測った密度の範囲から決められる。
:::

### 接触結合と二体の散乱長を同じ cutoff で結ぶ

体積規格化の箱で $\widetilde V(\vb q)=g_0$ を低波数の台内で用いる。相対運動の二体 $T$ を定数として LS の反復をまとめると
$$
T(E)^{-1}=g_0^{-1}-\int_{|\vb q|<\Lambda}
\frac{\dd^3q}{(2\pi)^3}\frac1{E+i0-2\epsilon_q}
$$
である。これは各中間波数で同じ接触係数を左と右から掛けるため、積分方程式がスカラーになることによる。$E=0$ で[散乱理論](/physics/advanced-quantum-mechanics/scattering/) の振幅係数に換算すると、同質量で換算質量 $m/2$ より $T(0)=4\pi\hbar^2a_s/m\defeq g$ である。従って
$$
\frac1{g_0}=\frac1g-\int_{|\vb q|<\Lambda}\frac{\dd^3q}{(2\pi)^3}\frac1{2\epsilon_q}
=\frac1g-\frac{m\Lambda}{2\pi^2\hbar^2}.
$$
ここで $g_0$ と同じ cutoff を Bogoliubov 展開にも使う。平均場 $g_0 n^2/2$ を $g$ の二次まで展開すると $+(gn)^2\int(2\epsilon_q)^{-1}/2$ が現れる。一方、零点エネルギーの $\tfrac12[E_q-\epsilon_q-gn]$ の大波数項は $-\tfrac12(gn)^2/(2\epsilon_q)$ である。両者を足すと発散項が消える。これは有限 cutoff での結合展開の同じ次数を組み合わせてから $\Lambda$ を外す操作であり、$g_0$ を固定したまま無制限な波数を足すこととは異なる。

:::example[波数 cutoff を変えると裸の係数が動く]
無次元な組 $2a_s\Lambda/\pi$ を $1/5$ とすると $g_0=g/(1-1/5)=5g/4$ である。cutoff を半分にすれば $g_0=10g/9$ であるが、両方の $T(0)$ は $g$ に揃う。有限な cutoff の計算で $g_0$ を変えずに比較すると、実際に指定した散乱長まで変わることになる。多体展開ではさらに $na_s^3\ll1$ が必要で、二体の係数が揃うだけで高次の多体効果が無視できるわけではない。
:::

## 超流動と Landau の判定条件

速度 $\vb v$ で流れる流体が静止した障害物に運動量を渡すことを考える。Galilei 不変な流体を考え、流体の静止系で等方的な励起のエネルギーと運動量を $\varepsilon(p),\vb p$ とすると、Galilei 変換の運動エネルギーの交差項により障害物の系での変化は $\varepsilon(p)+\vb v\cdot\vb p$ である。最小は $\vb p$ が流れと反対のときで、負のエネルギー変化が可能になる境界は

$$
v_{\mathrm L}=\inf_{p>0}\frac{\varepsilon(p)}{p}.
\label{eq-aqm-landau-velocity}
$$

これは :term[Landau の判定条件]{en="Landau criterion" yomi="らんだうのはんていじょうけん"} である。励起を作る行列要素と実際の境界が到達する過程も必要である。渦や壁の粗さが重要な系では、ここで用いた単一のバルク励起分散だけでは測定される臨界速度を決められない。

:::example[Bogoliubov 音波]
$p=\hbar k$ として :eqref[eq-aqm-bose-dispersion] を割ると $E_k/p=\sqrt{c_s^2+p^2/(4m^2)}$ である。単調増加なので下限は $c_s$ であり、$gn$ を 4 倍にすると $v_{\mathrm L}$ は 2 倍になる。
:::

:::example[自由 Bose 気体]
$\varepsilon=p^2/(2m)$ の比は $p/(2m)$ であり下限が零である。凝縮は存在し得るが、運動量を十分細かく変えられる熱力学極限でこの判定条件の正の閾値はない。
:::

### 有限の質量を持つ不純物の反跳

流体に対して速度 $\vb v$ を持つ質量 $M$ の不純物が運動量 $\vb p$ を励起に渡すと、失う運動エネルギーは $\vb v\cdot\vb p-p^2/(2M)$ である。これは初期 $\vb P=M\vb v$ と終状態 $\vb P-\vb p$ の二乗の差を $2M$ で割った値である。従って必要な速度は
$$
v\ge \inf_{p>0}\left[\frac{\varepsilon(p)}p+\frac p{2M}\right].
$$
エネルギー保存を最も満たしやすいのは $\vb p$ と $\vb v$ が平行な方向である。無限に重い外部の撹乱では反跳項が消える。どちらの場合も結合する行列要素が零ならその励起は作れないため、運動学的な許容条件と実際の散逸率を区別する。

:::example[ギャップのある枝での反跳]
$\varepsilon(p)=\Delta+p^2/(2m_*)$ とすると、最小化する式は $\Delta/p+p/(2\mu_*)$、$\mu_*^{-1}=m_*^{-1}+M^{-1}$ である。導関数から $p_*=\sqrt{2\mu_*\Delta}$、臨界速度は $\sqrt{2\Delta/\mu_*}$ となる。$M=m_*$ では $\mu_*=m_*/2$ なので、重い撹乱の結果より速度は $\sqrt2$ 倍で、選ばれる運動量は $1/\sqrt2$ 倍になる。
:::

### ギャップを持つ最小点の近くでの Landau 速度

分散が有限運動量 $p_0>0$ に最小値を持ち、その近くで $\varepsilon(p)=\Delta+(p-p_0)^2/(2m_*)$ と近似できるとする。$\varepsilon(p)/p$ の微分を零にする条件は $p\varepsilon'(p)=\varepsilon(p)$ である。これは原点から分散へ引いた接線の傾きがそこでの比と一致することに対応する。この模型を代入すると $p^2=p_0^2+2m_*\Delta$ で、対応する比は
$$
v_*=\frac{\sqrt{p_0^2+2m_*\Delta}-p_0}{m_*}.
$$
もし $\Delta m_*/p_0^2\ll1$ なら $v_*\simeq\Delta/p_0$ だが、厳密な最小比は分散の最低点そのものではなく、その少し右にある。この停留点が二次模型の範囲内にない場合や、他に音波などの枝がある場合は、全ての枝に対して比を取り直す。

:::example[低速音波と有限運動量の枝を比べる]
$m_*\Delta/p_0^2=1/8$ なら $p_*/p_0=\sqrt5/2$、$m_*v_*/p_0=(\sqrt5-2)/2\simeq0.118$ である。最小点での近似 $\Delta/p_0$ の無次元値は $1/8=0.125$ である。別に音速 $m_*c_s/p_0=0.15$ の枝があり、その比の下限が音速なら、有限運動量の枝の方が先にエネルギー・運動量保存を満たす。音速が 0.1 に下がる場合は音波の方が先である。
:::

### 音速を超えた動く弱い外場が選ぶモード

流体静止系で弱い外場 $U(\vb r-\vb v t)$ を Fourier に分ける。一成分は $U_{\vb q}e^{i\vb q\cdot\vb r-i\vb q\cdot\vb v t}$ で、密度に結合する。一次の時間積分の長時間極限のデルタは $E_{\vb q}-\hbar\vb q\cdot\vb v=0$ を課す。すなわち外場が供給する周波数は $\omega=\vb q\cdot\vb v$ であり、この周波数で密度に結合できるモードだけが励起される。理想的な無限時間と均一流体を使い、多準粒子や有限時間の幅を無視する。$\vb v$ の方向を軸にして $\theta$ を励起波数との角度とすれば、許される条件は $\cos\theta=\omega(q)/(qv)$ であり、$|\cos\theta|\le1$ が必要である。

:::example[Bogoliubov 分散で許される波数の上端]
$\omega(q)/q=c_s\sqrt{1+(q\xi)^2/2}$ とし、$v=\sqrt2c_s$ を取る。$\omega/(qv)\le1$ より $q\xi\le\sqrt2$ である。例えば $q\xi=1$ なら $\cos\theta=\sqrt3/2$ なので $\theta=\pi/6$ となる。$q\to0$ の角度は $\arccos(c_s/v)=\pi/4$、上端では $\theta=0$ となる。$v<c_s$ ならどの非零波数にも解がない。実際の発生強度は、ここで許された波数に対して $U_{\vb q}$ の大きさと密度の行列要素も掛けたものである。
:::

## 密度相関関数と動的構造因子

$\rho_{\vb q}=\sum_j e^{i\vb q\cdot\vb r_j}$ を密度 Fourier 成分とする。零温度で基底 $\ket0$ と固有状態 $\ket n$ を用いて

$$
S(\vb q,\omega)=\frac1N\sum_n |\mel n{\rho_{\vb q}}0|^2
   \delta(\omega-(E_n-E_0)/\hbar)
\label{eq-aqm-structure-factor}
$$

と定め、これを :term[動的構造因子]{en="dynamic structure factor" yomi="どうてきこうぞういんし"} と呼ぶ。Heisenberg 時間で $\langle\rho_{-\vb q}(t)\rho_{\vb q}(0)\rangle$ に固有状態を挿入すると各項に $e^{-i(E_n-E_0)t/\hbar}$ が付く。従って $S$ はこの相関を $\int e^{i\omega t}\dd t/(2\pi N)$ で変換した値である。$\omega$ 積分は静的構造因子 $\langle\rho_{-\vb q}\rho_{\vb q}\rangle/N$ を与える。

座標だけの相互作用に対し $\comm{\rho_{-\vb q}}{[H,\rho_{\vb q}]}=N\hbar^2q^2/m$ である。確かめるには $e^{-i\vb q\cdot\vb r}\vb p\,e^{i\vb q\cdot\vb r}=\vb p+\hbar\vb q$ を各粒子で用いる。同じ粒子の $\vb p^2/(2m)$ が出す二つの運動量シフトの和が $\hbar^2 q^2/m$ になる。Hamiltonian と基底の射線が空間反転で不変なら、励起にも同じエネルギーの反転像があるため $S(-\vb q,\omega)=S(\vb q,\omega)$ である。二重交換子に完全性を挿入すると二方向の一次モーメントの和になるため、反転対称な場合は

$$
\int_0^\infty \hbar\omega S(\vb q,\omega)\dd\omega=\frac{\hbar^2 q^2}{2m}
$$

となる。これは密度の $f$ 和則である。

:::example[一準粒子のピーク]
凝縮近似で $\rho_{\vb q}\simeq\sqrt N(a_{\vb q}^\dagger+a_{-\vb q})$ である。真空に作用する項の係数は $\sqrt N(u_q-v_q)$。$(u-v)^2=e^{-2r}=\epsilon_q/E_q$ なので $S\simeq(\epsilon_q/E_q)\delta(\omega-E_q/\hbar)$ となり、一次モーメントは $\epsilon_q$ と和則を飽和する。
:::

:::example[点 $\vb q=0$]
$\rho_0=\hat N$ であり、固定粒子数では $\rho_0\ket0=N\ket0$ である。従って $S(0,\omega)=N\delta(\omega)$ であり、弾性寄与だけが残る。一次モーメントは零である。
:::

### 熱平衡での詳細釣り合いと線形応答

熱平衡で各初期固有状態 $a$ の重みを $p_a=e^{-\beta E_a}/Z$ とする。:eqref[eq-aqm-structure-factor] を初期と終状態についての和に拡張すると
$$
S_\beta(\vb q,\omega)=\frac1N
\sum_{a,b}p_a|\mel b{\rho_{\vb q}}a|^2
 \delta(\omega-(E_b-E_a)/\hbar)
$$
となる。負周波数側で $\vb q$ の向きを反転し、$a,b$ を交換すると、$\rho_{-\vb q}=\rho_{\vb q}^\dagger$ から同じ行列要素が出る。ただし重みは $p_b=e^{-\beta\hbar\omega}p_a$ になる。従って $S_\beta(-\vb q,-\omega)=e^{-\beta\hbar\omega}S_\beta(\vb q,\omega)$ である。これが熱平衡の詳細釣り合いであり、負周波数は物質が励起を失う散乱に対応する。

弱い外部摂動として複素成分 $U(t)\rho_{\vb q}$ を Hamiltonian に足し、共役成分を必要に応じて加える。密度 $\rho_{-\vb q}/N$ の一次変化を求めるには、相互作用表示の密度演算子を $\delta\rho=-(i/\hbar)\int[U(t')\rho_{\vb q}(t'),\rho_\beta]\dd t'$ とする。測定演算子を掛けてトレースし、巡回性で交換子の位置を移せば応答核
$$
\chi_{\rm R}(\vb q,t)=-\frac{i}{N\hbar}\mathbf1_{t>0}
  \langle[\rho_{-\vb q}(t),\rho_{\vb q}(0)]\rangle_\beta
$$
を得る。二項を固有基底で展開し、片側 Fourier 変換 $\int_0^\infty e^{i(\omega+i0)t}\dd t$ を使うと、分母の虚部が $-\pi\delta$ となる。重みの差は $p_a-p_b=p_a(1-e^{-\beta\hbar\omega})$ なので
$$
\Im\chi_{\rm R}(\vb q,\omega)
=-\frac{\pi}{\hbar}(1-e^{-\beta\hbar\omega})S_\beta(\vb q,\omega)
$$
となる。熱揺らぎとして測る構造因子と、外場による散逸を表す応答の虚部を結ぶ式である。符号はここで外場を $+U\rho_{\vb q}$ として加えた規約に従う。例えば $\beta\hbar\omega=\log4$ では反対符号のピーク比は $1/4$ で、応答には両ピークの差に対応する $3/4$ の因子が現れる。

### 静的応答を逆エネルギーのモーメントで求める

基底状態に密度に結合する静的外場を加える。励起エネルギーを $\varepsilon_\nu=E_\nu-E_0>0$ とすれば、一次の状態補正に $\mel{\nu}{\rho_{\vb k}}0/(-\varepsilon_\nu)$ が掛かる。$\vb k$ と $-\vb k$ に対称性があり、$S_E(k,E)=S(k,E/\hbar)/\hbar$ を単位エネルギー当たりの重みとすると、密度の静的な応答は粒子一個当たりで
$$
\chi(k,0)=-2\int_0^\infty\frac{S_E(k,E)}{E}\dd E
$$
になる。二つの項は状態補正を bra と ket のどちらに入れるかに対応する。Bogoliubov の一極近似で $S_E=(\epsilon_k/E_k)\delta(E-E_k)$ を使うと $\chi(k,0)=-2\epsilon_k/E_k^2=-2/(\epsilon_k+2gn)$ になる。単位体積当たりはこれに $n$ を掛けるので、$k\to0$ で $-1/g$ に近づき、静的な $\partial n/\partial\mu=1/g$ に対応する。

:::example[同じ積分強度でも静的応答が異なる場合]
二つのピークをエネルギー $E_0,3E_0$ に置き、重みを $1/2$ ずつとする。全重みは 1、エネルギー一次モーメントは $2E_0$ であるが、逆エネルギーモーメントは $2/(3E_0)$ である。同じ全重みと一次モーメントを持つ単一ピークを $2E_0$ に置くと逆モーメントは $1/(2E_0)$ となる。従って前者の静的応答の絶対値は $4/3$ 倍である。低エネルギーの重みを解像しないと、全重みと平均エネルギーだけでは感受率を復元できない。
:::

### 理想 Fermi 気体の静的構造因子

零温度で $\rho_{\vb q}=\sum_{\vb k,\sigma}a^\dagger_{\vb k+\vb q,\sigma}a_{\vb k,\sigma}$ を作用させる。$\vb q\ne0$ では各 $\vb k$ で占有スピン軌道を一つ除き、空いた $\vb k+\vb q$ を埋めた場合だけ寄与する。異なる $\vb k,\sigma$ の行列状態は直交なので、静的な全重みは $\sum_{\vb k,\sigma} n_{\vb k\sigma}(1-n_{\vb k+\vb q,\sigma})/N$ である。無偏極で同じ大きさの Fermi 球なら、球全体のうち球を $\vb q$ だけずらしたとき重ならない体積の割合である。交換エネルギーで求めた重なり体積を使うと
$$
S(q)=\begin{cases}
3q/(4k_F)-q^3/(16k_F^3),&0<q<2k_F,\\
1,&q\ge2k_F.
\end{cases}
$$
前式は $1-V_\cap(q)/(4\pi k_F^3/3)$ を展開したもので、$q=2k_F$ で 1 に達する。$q=0$ そのものには固定粒子数の弾性デルタがあるので、この極限の式と区別する。自由な Fermi 気体で低波数の重みが線形なのは、Fermi 球の表面近くの一部だけが空の状態へ移れるためである。

:::example[同じ波数で粒子・空孔の割合を求める]
$q=k_F$ では静的重みは $3/4-1/16=11/16$ である。残り $5/16$ は終波数が既占有で行列要素が零になる割合である。二つの運動量の間の励起エネルギーは $\hbar^2(q^2+2\vb k\cdot\vb q)/(2m)$ であり、許容された $\vb k$ 全体でばらつく。Bogoliubov の単一モード近似で同じ静的重みを作ったとしても、エネルギースペクトルまで同じにはならない。
:::

### 自由粒子・空孔の連続帯を動的に計算する

理想 Fermi 気体で $q\ge2k_F$ を取れば、占有球上の各初期波数 $\vb k$ から $\vb k+\vb q$ は空である。ただし球面上の一点の例外は積分に寄与しない。軸を $\vb q$ に合わせ $k_\parallel$ とする。励起エネルギーは $E=E_q+\hbar^2q k_\parallel/m$、$E_q=\hbar^2q^2/(2m)$ である。固定 $k_\parallel$ での初期の切片の面積は $\pi(k_F^2-k_\parallel^2)$ である。球の体積 $4\pi k_F^3/3$ で割り、$\dd k_\parallel/\dd E=m/(\hbar^2q)$ を掛ける。$W=\hbar^2q k_F/m$ と置けば単位エネルギー当たりの重みは
$$
S_E(q,E)=\frac3{4W}\left[1-\left(\frac{E-E_q}{W}\right)^2\right]\mathbf1_{|E-E_q|<W}.
$$
積分は $3/4\int_{-1}^1(1-u^2)\dd u=1$ であり、一次モーメントは奇な補正が消えて $E_q$ になる。従って静的な構造因子と $f$ 和則を同じ連続分布で満たす。$0<q<2k_F$ では終状態が Fermi 球の内側にある部分を各切片から引く。$u=(E-E_q)/W$、$\eta=q/k_F$ とすると、初状態の縦成分は $k_F u$、終状態の縦成分は $k_F(u+\eta)$ である。$E>0$ なら後者の絶対値の方が大きいので、除くべき切片は初状態の円板の中に同心円板として収まる。したがって全波数範囲で
$$
S_E(q,E)=\frac3{4W}\big[1-u^2-\max(0,1-(u+\eta)^2)\big]\mathbf1_{|u|<1}\mathbf1_{E>0}
$$
である。$-\eta/2<u<1-\eta$ では括弧が $2\eta u+\eta^2$ になり、$1-\eta<u<1$ では $1-u^2$ に戻る。積分を二つの区間に分けると、全重みは $3\eta/4-\eta^3/16$、エネルギー一次モーメントは $E_q$ となる。これは球の重なりと $f$ 和則から得た値である。例えば $q=k_F$ のとき $E=E_F/2$ と $E_F$ での密度はそれぞれ $3/(16E_F)$、$3/(8E_F)$ である。

:::example[一極に近いかを相対幅で評価する]
$q=3k_F$ のとき $E_q=9E_F$、$W=6E_F$ であり、帯の下と上は $3E_F,15E_F$ になる。標準偏差は $W/\sqrt5$ である。実際、$\langle u^2\rangle=(3/4)\int_{-1}^1(u^2-u^4)\dd u=1/5$ となる。平均に対する標準偏差は $2/(3\sqrt5)\simeq0.298$ である。$q=10k_F$ にするとこの比は $1/(5\sqrt5)$ に下がり、相対幅は狭くなる一方、エネルギーの絶対幅 $W$ は増える。
:::

### 自由電子の静的な応答

体積当たりの応答 $\chi_0$ を求める。密度演算子の粒子・空孔展開を、前の遅延応答の交換子に入れる。二項では逆の遷移の占有因子が $n_{\vb k}(1-n_{\vb k+\vb q})$ と $n_{\vb k+\vb q}(1-n_{\vb k})$ になる。差をとれば無偏極な理想 Fermi 気体について
$$
\chi_0(q,\omega)=2\int\frac{\dd^3k}{(2\pi)^3}
\frac{n_{\vb k}-n_{\vb k+\vb q}}
{\hbar(\omega+i0)+\epsilon_k-\epsilon_{|\vb k+\vb q|}},
\qquad n_{\vb k}=\mathbf1_{|\vb k|<k_F}
$$
である。これを :term[Lindhard 関数]{en="Lindhard function" yomi="りんどはるとかんすう"} と呼ぶ。係数 2 はスピンの和である。$\omega=0$ で第 2 項の積分変数を $\vb k+\vb q$ に変えると、両項を同じ Fermi 球に制限できる。軸を $\vb q$ にそろえ、対称な角度積分を主値で行うと
$$
\chi_0(q,0)=-\frac{m}{\pi^2\hbar^2 q}\int_0^{k_F} k\log\left|\frac{q+2k}{q-2k}\right|\dd k
$$
となる。たとえば第 1 項の分母は $-\hbar^2q(q+2k\cos\theta)/(2m)$ であり、第 2 項も $\cos\theta$ の反転で同じ実部を与える。虚部は零周波数なので消える。交換自己エネルギーで用いた対数の積分に $q/2$ を入れると、$x=q/(2k_F)$ として
$$
\chi_0(q,0)=-\frac{mk_F}{\pi^2\hbar^2}
\left[\frac12+\frac{1-x^2}{4x}\log\left|\frac{1+x}{1-x}\right|\right]
$$
を得る。$x\to0$ の括弧は 1 であり、$n=k_F^3/(3\pi^2)$ と $\mu=\hbar^2k_F^2/(2m)$ から求めた $-\partial n/\partial\mu$ に一致する。$x=1$ では括弧が連続に $1/2$ となるが、対数により傾きは発散する。

### 自己無撞着な密度応答を閉じる

外場 $U_{\rm ext}$ による密度の変動 $\delta n$ に対し、Hartree の変動を Fourier 空間で $v(q)\delta n$ とする。自由な一次応答なら $\delta n=\chi_0(q,\omega)[U_{\rm ext}+v(q)\delta n]$ である。従って $\chi=\chi_0/[1-v(q)\chi_0]$ が得られる。これは揺らぐ密度を平均場として反復する :term[ランダム位相近似]{en="random phase approximation" yomi="らんだむいそうきんじ"} であり、交換の頂点や多体相関の補正を全て含むわけではない。高周波で $\chi_0\simeq nq^2/(m\omega^2)$ となるのは、応答のスペクトルの分母を展開し定数項の交換子が零、次の項が $f$ 和則で決まるためである。Coulomb の $v(q)=e^2/(\epsilon_0 q^2)$ を使うと分母の零は $\omega_p^2=ne^2/(m\epsilon_0)$ に近づく。個別の粒子・空孔の連続帯とは異なる集団的な密度振動である。小さい $q$ により粒子・空孔の典型周波数より $\omega_p$ が十分大きい領域でこの展開を使う。

:::example[プラズマの尺度を自由粒子の帯から切り分ける]
Coulomb の自己無撞着な分母で $q\to0$ の根を考える。電子密度を四倍にすると $\omega_p=\sqrt{ne^2/(m\epsilon_0)}$ は二倍になるが、同じ波数での自由な粒子・空孔の帯の幅の尺度 $qv_F=\hbar q k_F/m$ は $4^{1/3}$ 倍である。従って $\omega_p/(qv_F)$ は $2^{1/3}$ 倍になる。二つの周波数の密度の冪が違うのは、集団振動では電荷密度が復元力を作る一方、粒子・空孔の励起幅は占有球の半径から決まるためである。
:::

## スピン波と Holstein–Primakoff 変換

強磁性の整列状態からの小さいずれを、各サイトのスピン反転数で表す。無次元スピン $\vb s_i=\vb J_i/\hbar$ を使い、長さ $S$ の表現を $0\le n_i\le2S$ のボソン部分空間と同一視する。

:::definition[Holstein–Primakoff 変換]{#def-aqm-holstein-primakoff}
$N_i=a_i^\dagger a_i$ として

$$
s_i^z=S-N_i,\qquad
s_i^+=\sqrt{2S-N_i}\,a_i,\qquad
s_i^-=a_i^\dagger\sqrt{2S-N_i}
$$

と表す。この :term[Holstein–Primakoff 変換]{en="Holstein–Primakoff transformation" yomi="ほるすたいんぷりまこふへんかん"} は物理部分空間で厳密である。
:::

$m=S-n$ なので $s^+\ket n=\sqrt{n(2S-n+1)}\ket{n-1}$ となり、梯子係数と一致する。$n=2S$ で $s^-$ が零になるから、物理部分空間はすべてのスピン作用で閉じる。低密度 $n\ll2S$ なら平方根を $\sqrt{2S}$ で近似する。

最近接結合の $H=-J\sum_{\langle ij\rangle}\vb s_i\cdot\vb s_j$、$J>0$ とする。各結合は一回だけ数える。縦成分は $S^2-S(N_i+N_j)$、横成分 $\tfrac12(s_i^+s_j^-+s_i^-s_j^+)$ は $S(a_i^\dagger a_j+a_j^\dagger a_i)$ となる。従って

$$
H_2=E_{\mathrm{pol}}+JS\sum_{\langle ij\rangle}
 (a_i^\dagger-a_j^\dagger)(a_i-a_j).
$$

超立方格子で格子間隔 $a$、次元 $d$ なら正方向の隣接だけを和に入れる。Fourier 変換で一結合方向の因子は $|1-e^{ik_\alpha a}|^2$ だから、:term[スピン波]{en="spin wave" yomi="すぴんは"} の分散は $E_{\vb k}=2JS\sum_{\alpha=1}^d(1-\cos(k_\alpha a))$ である。低波数では $JSa^2 k^2$ となり、周波数は $E_{\vb k}/\hbar$ である。一個の反転がある部分空間では、二個以上の反転を必要とする相互作用項が作用しないため、この分散は厳密になる。

:::example[一次元の帯域幅]
$k=0$ で $E=0$、$k=\pi/a$ で $E=4JS$ である。例えば $J=2\,\mathrm{meV}$、$S=1/2$ の帯域幅は $4\,\mathrm{meV}$ になる。
:::

:::example[一様磁場とギャップ]
$H_B=-g_m\mu_B B\sum_i s_i^z$ を足せば、一準粒子のエネルギーに $g_m\mu_B B$ が加わる。磁場の向きに整列させた場合には $k=0$ のエネルギーが正となり、一様なスピン歳差の周波数は $g_m\mu_B B/\hbar$ になる。
:::

### 二副格子での反強磁性スピン波

最近接反強磁性 $H=J\sum_{\langle i,j\rangle}\vb s_i\cdot\vb s_j$、$J>0$ とし、格子を二つの副格子 $A,B$ に分ける。古典 Néel 配置で $A$ は上向き、$B$ は下向きである。各局所軸に Holstein–Primakoff を適用すると実験室系では $s_A^z=S-a^\dagger a$、$s_B^z=-S+b^\dagger b$、$s_A^+\simeq\sqrt{2S}a$、$s_B^+\simeq\sqrt{2S}b^\dagger$ である。一結合の二次項は $JS(a^\dagger a+b^\dagger b+ab+a^\dagger b^\dagger)$ になる。近接数を $z$、$\gamma_{\vb k}=z^{-1}\sum_{\delta}e^{i\vb k\cdot\delta}$ とすれば Fourier 変換後の二モード対の係数は $A=zJS$、$B=zJS\gamma_{\vb k}$ である。格子が反転対称で $\gamma$ が実の場合に $b_{\vb k}$ と $a_{-\vb k}^\dagger$ の混合を行い、
$$
E_{\vb k}=zJS\sqrt{1-\gamma_{\vb k}^2}
$$
となる。正則な波数で $\tanh2r=\gamma_{\vb k}$ とすれば対項が消え、零周波数点は分散の極限として扱う。波数は二副格子を一組にした単位胞の逆格子で一度ずつ数え、各波数に二つの準粒子モードを置く。磁化を議論するときは、微弱な交代磁場などで先に整列軸を定め、熱力学極限でその場を外す。二次理論の真空でも各副格子の元の反転数は $v_{\vb k}^2=(1/\sqrt{1-\gamma_{\vb k}^2}-1)/2$ である。従って有限系の零モードを除き、副格子上の磁化の減少は $2\sum_{\vb k}v_{\vb k}^2/N_{\rm site}$ となる。無ギャップ点の近くでは $1-\gamma_{\vb k}^2$ が $|\vb k|^2$ に比例するので、零温度の小波数積分は $\int_0 k^{d-2}\dd k$ の形になる。一次元では対数発散し、単純な Néel 配置のまわりで反転密度が小さいという仮定が破れる。

:::example[平方格子の二つの波数]
格子定数 $a$ の平方格子で $\gamma=(\cos k_xa+\cos k_ya)/2$ である。$(\pi/(2a),0)$ で $\gamma=1/2$ だから $E=2\sqrt3\,JS$、$(\pi/a,0)$ で $E=4JS$ になる。原点近傍では $\gamma\simeq1-a^2k^2/4$ なので $E\simeq2\sqrt2\,JSa|\vb k|$ である。強磁性の二次分散と異なるのは、局所軸が交互に逆で対生成項がすでに二次式に現れるためである。
:::

### 強磁性の一マグノンでは二次式を厳密に照合できる

最近接強磁性 $H=-J\sum_{\langle i,j\rangle}\vb s_i\cdot\vb s_j$、$J>0$ を対象とする。整列状態 $\ket F$ は全結合で $S^2$ を与える。規格化した一反転 $\ket i=s_i^-\ket F/\sqrt{2S}$ を作る。サイト $i$ と隣の結合では $s_i^z$ が $S-1$ となるため、対角は一結合当たり $JS$ 上がる。一方 $-J s_i^+s_j^-/2$ の係数は二つの $\sqrt{2S}$ の積から $-JS$ であり、反転を隣へ移す。その他の結合の動作は整列と同じである。従って $H-E_F$ は一反転空間で格子 Laplacian $JS(zI-A)$ と等しく、$A$ は隣接行列である。二次 Holstein–Primakoff の予測がこの空間では線形化の誤差なしに成り立つ。二つ以上の反転が隣接すると縦成分の積や平方根が変わるため、独立なマグノンだけの和は一般に厳密でない。

:::example[四サイトの周期鎖で固有値を計算する]
$N=4$、周期境界で波数 $ka=0,\pi/2,\pi,3\pi/2$ を取る。隣接行列の各固有値は $2\cos ka=2,0,-2,0$ なので励起は $0,2JS,4JS,2JS$ となる。初期反転が一サイトにあるなら全四波数の重みは各 $1/4$ で、平均励起は $2JS$ である。これは局所の反転に接する二結合が各 $JS$ 上がるという計算と一致する。$ka=0$ の零は全スピンを下げた状態で、回転対称性による同じ多重項に入る。
:::

### マグノンの熱占有が磁化に与える次数

三次元の強磁性で低波数の分散を $E_{\vb k}=Dk^2$、$D>0$ と近似する。各マグノンは $s^z$ を一つ減らすので、単位体積当たりの減少は
$$
n_{\mathrm m}=\int\frac{\dd^3k}{(2\pi)^3}\frac1{e^{\beta Dk^2}-1}
=\frac{\zeta(3/2)}{(4\pi\beta D)^{3/2}} .
$$
第二式は $1/(e^x-1)=\sum_{n\ge1} e^{-nx}$ を入れ、各三次元 Gaussian の積分 $(4\pi n\beta D)^{-3/2}$ を足して得る。正の項の和なのでこの交換は許される。実際の格子では第一 Brillouin 域に限るが、温度が十分低ければ典型波数 $(\beta D)^{-1/2}$ が格子の逆長より小さく、上端の延長は先頭次数を変えない。二次元では小波数で被積分関数が $k^{-2}$、動径測度が $k\dd k$ なので対数発散する。これは異方性ギャップのない二次近似が有限温度で有限な減少を与えないことを意味する。これだけで相互作用系の一般の定理が証明できるわけではない。

:::example[温度とギャップを変えたときの磁化]
ギャップのない三次元で温度を四倍にすると、低温のマグノン密度は $(4)^{3/2}=8$ 倍になる。分散を $\Delta+D k^2$ に替えれば、同じ積分は $\sum_{n\ge1}e^{-n\beta\Delta}n^{-3/2}$ になり、$\beta\Delta\gg1$ では $e^{-\beta\Delta}$ が先頭になる。無ギャップの冪則とギャップありの抑制は同じ分散の積分で区別できる。どちらも、磁化の減少が格子の全スピン密度より十分小さい領域での計算である。
:::

## 演習問題

::::exercise{level=basic}
占有軌道二個について直接積分 $5\,\mathrm{eV}$、交換積分 $2\,\mathrm{eV}$ とする。全二体エネルギーを求めよ。
:::solution
自己項は相殺する。$(1,2),(2,1)$ の二組で同じ値 $5-2$ が出るので、$1/2$ を掛けると $3\,\mathrm{eV}$ である。
:::
::::

::::exercise{level=basic}
無偏極電子ガスの密度を 27 倍にする。Fermi 波数と交換エネルギーの絶対値の倍率を求めよ。
:::solution
$n=k_F^3/(3\pi^2)$ から $k_F$ は 3 倍である。交換の一粒子当たりの絶対値は $k_F$ に比例するので 3 倍、体積当たりの絶対値は $n$ の増分も含み 81 倍になる。
:::
::::

::::exercise{level=basic}
$S=1$ のサイトで $n=1$ の状態に $s^-$ を作用させよ。
:::solution
$\sqrt{2S-N}$ は 1、$a^\dagger$ は $\sqrt2$ を掛けるので $\sqrt2\ket{n=2}$ である。さらに $s^-$ を作用させると平方根が零になる。
:::
::::

::::exercise{level=basic}
$\epsilon_k=2gn$ の Bose 気体の励起エネルギーを $gn$ の単位で求めよ。
:::solution
$E_k^2=\epsilon_k(\epsilon_k+2gn)=8(gn)^2$ なので $E_k/(gn)=2\sqrt2$ である。
:::
::::

::::exercise{level=advanced}
励起にギャップ $\Delta>0$ があり、$\varepsilon(p)=\Delta+p^2/(2M)$ とする。Landau 速度を計算せよ。
:::solution
比の微分は $-\Delta/p^2+1/(2M)$ である。$p=\sqrt{2M\Delta}$ で零、そこまで減少しその後増加する。従って最小値は $\sqrt{2\Delta/M}$ である。
:::
::::

::::exercise{level=advanced}
一様電子ガスの重なり体積を積分し、交換の尺度 $E_x/\mathcal V\propto-k_F^4$ を求めよ。
:::solution
$V_\cap(q)=\pi(16K^3-12K^2q+q^3)/12$ を $0$ から $2K$ まで積分すると $\pi(32-24+4)K^4/12=\pi K^4$ である。分母 $q^2$ が測度の $q^2$ と相殺するためこれが波数依存性の全てである。
:::
::::

::::exercise{level=advanced}
Bogoliubov ピークの静的重みを $k\to0$ で展開せよ。
:::solution
$E_k\simeq\hbar c_s k$ で $\epsilon_k=\hbar^2 k^2/(2m)$ だから $\epsilon_k/E_k\simeq\hbar k/(2m c_s)$ である。零に向かって線形に減る。重みと励起エネルギーの積は $\epsilon_k$ を保つ。
:::
::::

::::exercise{level=advanced}
一次元強磁性鎖のスピン波の群速度の最大絶対値を求めよ。
:::solution
$\omega=2JS(1-\cos ka)/\hbar$ を $k$ で微分すると $v_g=2JSa\sin(ka)/\hbar$ である。$|ka|=\pi/2$ で絶対値が最大 $2JSa/\hbar$ となる。
:::
::::

::::exercise{level=advanced}
中和背景を持つ一様な無偏極電子ガスの Hartree–Fock で、ある密度 $n_0$ において $e_x(n_0)=-2e_t(n_0)$ とする。圧力と圧縮率の逆の符号を求めよ。さらに同じ状態族で密度を下げると圧縮率の逆が初めて零になる密度を求めよ。

:::solution
$e_t=t_0(n/n_0)^{2/3}$、$e_x=-2t_0(n/n_0)^{1/3}$ と置く。$t_0>0$ であり、運動と交換の密度の冪は Fermi 球の半径から出る。固定粒子数で体積を微分すれば $P/n=(2e_t+e_x)/3$、従って $n=n_0$ では零である。一方 $n^{-1} n\partial_nP=\partial_nP=(10e_t+4e_x)/9$ なので、その点では $2t_0/9>0$ である。圧力が零でも傾きは正であり得る。

密度比の三乗根 $u=(n/n_0)^{1/3}$ を用いると $10e_t+4e_x=t_0(10u^2-8u)$ である。正の根は $u=4/5$、密度は $(4/5)^3n_0=64n_0/125$ である。この点の圧力は $P/n=(2t_0/3)(u^2-u)=-8t_0/75$ である。低密度側で圧縮率の逆が負になる結果は、この一様な平均場の式を続けた結果であり、実際の基底状態の相や安定性をこの継続だけで決めない。
:::
::::

::::exercise{level=advanced}
分散 $E_{\vb q}=\hbar c_s q\sqrt{1+q^2\xi^2/2}$ の Bose 気体で、一次の弱い動く摂動が波数 $q$ の一モードを作ることを考える。摂動の速度 $v=3c_s/2$ のとき許される最大波数と $q\xi=1$ での余弦を求めよ。$q$ が十分小さいときと最大波数での方向も求めよ。

:::solution
運動量の方向とのなす角を $\theta$ とすると外場から渡るエネルギー $\hbar\vb q\cdot\vb v$ が励起エネルギーに等しい必要がある。従って
$$
\cos\theta=\frac{2}{3}\sqrt{1+\frac{q^2\xi^2}{2}} .
$$
右辺は正で波数とともに増えるので、許されるには右辺が 1 以下でなければならない。二乗して $1+(q\xi)^2/2\le9/4$、よって $q_{\max}\xi=\sqrt{5/2}$ である。$q\xi=1$ のとき余弦は $\tfrac23\sqrt{3/2}=\sqrt{2/3}$ である。$q\downarrow0$ の極限は $\cos\theta=2/3$、最大波数では $\theta=0$ になる。波数がこれを超えると、同じ速度では平行方向でも十分なエネルギーを渡せない。

密度の行列要素の二乗は一粒子当たりで $\epsilon_q/E_q$ であり、$q\downarrow0$ で $\hbar q/(2mc_s)$ に比例して小さくなる。したがって許される角が有限であることから、零波数で有限な励起の率があるとは言えない。実際の率には摂動の Fourier 成分と位相空間の重みも掛ける。
:::
::::

::::exercise{level=advanced}
理想 Fermi 気体で $q=3k_F$ とする。動的構造因子の重みを、エネルギーが平均より大きく、平均から半幅 $W/2$ までの範囲で積分せよ。その範囲と反対の範囲との比、および両範囲を足した重みを求めよ。

:::solution
$q>2k_F$ なので終状態の Pauli 禁止は出ず、$S_E\,\dd E=(3/4)(1-u^2)\dd u$、$u=(E-E_q)/W$ である。指定範囲は $0<u<1/2$ なので
$$
\int_0^{1/2}\frac34(1-u^2)\dd u
=\frac34\left(\frac12-\frac1{24}\right)=\frac{11}{32}.
$$
反対側 $-1/2<u<0$ の重みも被積分関数の偶性から同じであり、比は 1、和は $11/16$ である。平均 $E_q=9E_F$、$W=6E_F$ なので、二範囲を合わせると $6E_F<E<12E_F$ に当たる。なお、どちらか片側の範囲だけで計算した平均エネルギーは $E_q$ と異なる。実際、上側での $\langle u\rangle$ は $\int_0^{1/2}u(1-u^2)\dd u/\int_0^{1/2}(1-u^2)\dd u=(7/64)/(11/24)=21/88$ である。両側を同じ重みで足すと $E_q$ に戻る。全域積分の $f$ 和則と窓を選んだ測定を区別する。
:::
::::

## 確認テスト

::::quiz
Hartree–Fock の自己項が消える理由を選べ。

- [ ] Coulomb 力は同一点で零である
- [x] 同じ軌道に対する直接積分と交換積分が等しい
- [ ] 軌道エネルギーが零になる

:::explanation
$i=j$ の場合、$W_{ij;ij}$ と $W_{ij;ji}$ は同じ積分である。差を取ると零になり、同じフェルミオンを一対として数えない。
:::
::::

::::quiz{answer="4" tol="1e-9"}
一次元強磁性鎖で $E_{\pi/a}/(JS)$ を求めよ。
:::explanation
分散 $2JS(1-\cos ka)$ に $ka=\pi$ を代入すると $2JS(1-(-1))=4JS$ になる。
:::
::::

::::quiz
長波長の Bogoliubov 分散と、最近接強磁性の分散の波数依存性を選べ。

- [x] 前者は一次、後者は二次
- [ ] 両者とも一次
- [ ] 前者は二次、後者は一次

:::explanation
Bose 気体では $\epsilon_k\propto k^2$ が平方根の中にあり、$E_k$ は $k$ に比例する。強磁性では $1-\cos ka$ を展開して $k^2$ が先頭となる。
:::
::::

::::quiz{answer="2" tol="1e-9"}
相互作用パラメータ $gn$ が 4 倍になり質量が変わらない。音速の倍率を求めよ。
:::explanation
$c_s$ は $\sqrt{gn/m}$ であるため、比は $\sqrt4=2$ になる。
:::
::::

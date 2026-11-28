---
title: 宇宙論的摂動論
description: 一様宇宙の一次摂動を分解し、座標に依存しない重力と流体の進化式を導く。
notation:
  - symbol: '\eta'
    meaning: 共形時刻
  - symbol: '\mathcal H'
    meaning: 共形 Hubble パラメータ
  - symbol: '\Psi'
    meaning: 時間計量のポテンシャル
  - symbol: '\Phi'
    meaning: 空間曲率のポテンシャル
  - symbol: '\zeta'
    meaning: 一様密度面の曲率摂動
---

## この章の目標

- 計量の自由度を回転対称性によって分ける。
- ゲージ変換の影響を追い、Bardeen 変数を構成する。
- Einstein 制約と保存則を組み合わせて曲率摂動の保存条件を導く。

## 計量摂動のスカラー・ベクトル・テンソル分解

$c=\hbar=k_B=1$、計量符号 $(-,+,+,+)$ とし、空間的に平坦な背景の非零 Fourier モードを扱う。Fourier 核は $e^{i\vb k\cdot\vb x}$ とし、$k=|\vb k|$ とする。$\dd\eta=\dd t/a$、$\mathcal H=a'/a$ とし、プライムは $\eta$ 微分である。$k_i$ に沿う射影 $k_i k_j/k^2$ とその余射影でベクトルを縦と横に分ける。対称テンソルは跡、縦の跡なし成分、横ベクトルの微分、横・跡なし成分に分ける。これらは $k$ 周りの回転でそれぞれヘリシティ $0,\pm1,\pm2$ を持つので、背景が等方なら線形方程式の中で互いに混ざらない。計量を

$$
\dd s^2=a^2[-(1+2A)\dd\eta^2+2(\partial_i B+S_i)\dd\eta\,\dd x^i+((1-2C)\delta_{ij}+2\partial_i\partial_j E+\partial_iF_j+\partial_j F_i+h_{ij})\dd x^i\dd x^j] \label{eq-cosmo-svt}
$$

と書く。$\partial_i S_i=\partial_i F_i=\partial_i h_{ij}=h_{ii}=0$ である。独立な成分数はスカラー $4$、ベクトル $4$、テンソル $2$ である。$k=0$ ではこの分解が特異になり背景の再定義と区別する。

:::example[波数を $z$ 方向に取る場合]
テンソルは $h_{xx}=-h_{yy}$ と $h_{xy}=h_{yx}$ だけ独立である。
:::


:::example[縦速度]
フーリエ振幅が $(0,0,v)$ の速度で $k_i=(0,0,k)$ なら横成分は零である。
:::

### 波数ごとの直交射影

非零波数 $\vb k$ に直交する平面への射影を $P_{ij}=\delta_{ij}-\hat k_i\hat k_j$ とする。任意の対称行列 $X_{ij}$ の横・跡なし部分は
$$
X_{ij}^{\rm TT}=
\left(P_i{}^mP_j{}^n-\tfrac12P_{ij}P^{mn}\right)X_{mn}.
$$
第一項は両添字を横平面へ移し、第二項はその二次元の跡を引く。$P^2=P$、$P_i{}^i=2$ だから結果は横・跡なしであり、もう一度作用させても変わらない。残りの成分は平面内の跡 $P^{ij}X_{ij}$、法線成分 $\hat k^i\hat k^jX_{ij}$、混合成分 $P_i{}^m\hat k^nX_{mn}$ で決まる。これで対称行列の六自由度が $2+2+2$ に分解される。

:::example[波数が第三軸を向く場合]
$X=\diag(4,2,9)$ では横ブロックは $\diag(4,2)$ である。半分の跡 $3$ を引くと $X^{\rm TT}=\diag(1,-1,0)$ となる。第三成分の値を変えても横・跡なし成分は変わらない。
:::

### 射影による成分の復元

:term[スカラー・ベクトル・テンソル分解]{en="scalar-vector-tensor decomposition" yomi="すからーべくとるてんそるぶんかい"} を計量から復元する。$M_i=a^{-2}\delta g_{0i}$、$X_{ij}=a^{-2}\delta g_{ij}$ とおく。Fourier 空間で $M_i=ik_i B+S_i$ を $k^i$ と $P_i{}^j$ にそれぞれ縮約すると $B=-ik^jM_j/k^2$、$S_i=P_i{}^jM_j$ となる。空間成分について、横平面の跡は $P^{ij}X_{ij}=-4C$、縦成分は $\hat k^i\hat k^jX_{ij}=-2C-2k^2E$、混合成分は $P_i{}^j\hat k^lX_{jl}=ikF_i$ である。従って
$$
C=-\frac14P^{ij}X_{ij},\qquad
E=\frac{-\hat k^i\hat k^jX_{ij}-2C}{2k^2},\qquad
F_i=-\frac{i}{k}P_i{}^j\hat k^lX_{jl},\qquad h_{ij}=X_{ij}^{\rm TT}.
$$
最後に $A=-a^{-2}\delta g_{00}/2$ である。全変数を復元できるので分解は一意である。実空間では $P_{ij}=\delta_{ij}-\partial_i\partial_j\nabla^{-2}$ と読める。周期境界なら空間平均を除外し、無限空間なら無限遠での減衰を仮定して、逆 Laplace 作用素の零モードを指定する。

$\vb e_1,\vb e_2$ を波数軸に垂直な正規直交基底とし、$\vb e_\pm=(\vb e_1\pm i\vb e_2)/\sqrt2$ とする。回転角 $\alpha$ で $\vb e_\pm$ は逆向きの位相 $e^{\mp i\alpha}$ を持ち、そのテンソル積 $e_{\pm i}e_{\pm j}$ は $e^{\mp2i\alpha}$ を持つ。回転と可換な線形作用素の行列要素 $L_{mn}$ は任意の $\alpha$ で $(e^{im\alpha}-e^{in\alpha})L_{mn}=0$ を満たす。従って異なるヘリシティを結ぶ要素は零である。二次では二つの波数を持つ積が現れるため、一次の分離は一般に成り立たない。

:::example[横成分と縦成分が同時にある計量]
$\vb k=(0,0,k)$、$X_{xx}=5$、$X_{yy}=1$、$X_{zz}=-1$、$X_{xz}=X_{zx}=2$、他の成分は零とする。横平面の跡は $6$ なので $C=-3/2$ である。縦成分の式 $-1=3-2k^2E$ から $E=2/k^2$ を得る。また $ikF_x=2$ なので $F_x=-2i/k$ である。横方向の跡を引くと $h_{xx}=2$、$h_{yy}=-2$ となる。実計量を作るときには $-\vb k$ の係数も加え、二つを複素共役にする。
:::

:::example[偏光振幅の正規化]
波数を第三軸に取ると $e^+_{ij}=(e_{xi}e_{xj}-e_{yi}e_{yj})/\sqrt2$、$e^\times_{ij}=(e_{xi}e_{yj}+e_{yi}e_{xj})/\sqrt2$ は直交し、各自の縮約ノルムは $1$ である。$h_{xx}=b$、$h_{xy}=c$ なら $h_+=\sqrt2b$、$h_\times=\sqrt2c$ となる。$\sum_\lambda |h_\lambda|^2=2(|b|^2+|c|^2)=h^*_{ij}h_{ij}$ により、両偏光を足すときの規格化を確かめられる。
:::

### 応力の分解と幾何学的な曲率

物質の対称な空間応力の一次変化を $\delta T^i{}_j=\delta p\,\delta^i_j+\Pi^i{}_j$ と分ける。$\Pi^i{}_i=0$ であり、非零波数のスカラー部分は $\Pi^i{}_j=B_s Q^i{}_j$、$Q^i{}_j=\hat k^i\hat k_j-\delta^i_j/3$ である。$Q^i{}_j Q^j{}_i=2/3$ なので、異方性応力を $\sigma_f=-Q^i{}_j\Pi^j{}_i/(\bar\rho+\bar p)$ と規格化した場合は $B_s=-3(\bar\rho+\bar p)\sigma_f/2$ となる。従ってスカラー応力の縦の成分は $\Pi^i{}_j\hat k_i\hat k^j=-(\bar\rho+\bar p)\sigma_f$、各横成分はその反対符号の半分である。例えば波数が第三軸方向で $(\Pi^x{}_x,\Pi^y{}_y,\Pi^z{}_z)=(q,q,-2q)$ なら、$\sigma_f=2q/(\bar\rho+\bar p)$ である。

同じ射影で計量から空間の内在曲率を計算できる。平坦な三次元背景で小変化 $X_{ij}$ があるとき、Ricci スカラーの一次は $a^{-2}(\partial^i\partial^jX_{ij}-\nabla^2 X^i{}_i)$ である。これは線形接続の $\partial_k\Gamma^k{}_{ij}-\partial_j\Gamma^k{}_{ik}$ に $\delta^{ij}/a^2$ を掛けると得る。$X_{ij}=-2C\delta_{ij}+2\partial_i\partial_jE+\partial_iF_j+\partial_jF_i+h_{ij}$ を代入すると、$E$ の二項は $2\nabla^4 E$ 同士で相殺し、$F,h$ の項は横条件と跡を使えば零である。従って $\delta R^{(3)}=4a^{-2}\nabla^2 C$ だけが残る。$-C$ を曲率摂動と呼んでいても、それは Fourier 振幅を扱うための変数で、空間の Ricci スカラー自体には $k^2$ が掛かる。どの面で曲率を測ったかは別途指定する。

## ゲージ自由度とゲージ不変変数

受動変換 $\tilde\eta=\eta+T$、$\tilde x^i=x^i+\partial^iL$ で計量摂動は背景計量の Lie 微分だけ差し引かれる。各成分から $\tilde A=A-\mathcal H T-T'$、$\tilde C=C+\mathcal HT$、$\tilde B=B+T-L'$、$\tilde E=E-L$ となる。背景だけに依存する物理量の一次変化は $\widetilde{\delta\rho}=\delta\rho-\bar\rho'T$ となる。$\sigma=B-E'$ は $T$ だけ増えるので

$$
\Psi=A+\mathcal H\sigma+\sigma',\qquad \Phi=C-\mathcal H\sigma,\qquad \zeta=-C-\mathcal H\,\delta\rho/\bar\rho' . \label{eq-cosmo-bardeen}
$$

初めの二つは :term[Bardeen ポテンシャル]{en="Bardeen potentials" yomi="ばーでぃーんぽてんしゃる"} である。$L=E$、$T=-\sigma$ を選ぶと $B=E=0$ となり、:term[Newton ゲージ]{en="Newtonian gauge" yomi="にゅーとんげーじ"} の計量に $\Psi,\Phi$ が直接現れる。$\zeta$ は :term[一様密度面の曲率摂動]{en="curvature perturbation on uniform-density slices" yomi="いちようみつどめんのきょくりつせつどう"} であり、$\bar\rho'\ne0$ の場合に定義する。$T=\delta\rho/\bar\rho'$ を選んだ一様密度面では $-C$ になる。

:::example[純粋な座標変更]
背景で全摂動零から $T$ だけ動かすと $A=-\mathcal HT-T'$, $\sigma=T$ になり $\Psi=0$ のままである。
:::


:::example[時間のずれと密度]
$\bar\rho'<0$ なら $T>0$ の純粋な時間変更で座標上の $\delta\rho$ は正になる。
:::

### 共動面と一様密度面の比較

Newton ゲージから時刻を $T$ だけ変える。混合成分 $T^0{}_i$ に背景の Lie 微分を作用させると速度の発散は $\tilde\theta=\theta-k^2T$ となる。従って共動面は $T=\theta/k^2$ で得られる。その面の :term[共動曲率摂動]{en="comoving curvature perturbation" yomi="きょうどうきょくりつせつどう"} を $\mathcal R=-\tilde C$、密度を $\Delta\rho=\widetilde{\delta\rho}$ とおくと
$$
\mathcal R=-\Phi-\mathcal H\theta/k^2,\qquad
\Delta\rho=\delta\rho+3\mathcal H(\bar\rho+\bar p)\theta/k^2.
$$
:eqref[eq-cosmo-einstein-scalar] の最初の二式で時間項を消すと $k^2\Phi=-4\pi Ga^2\Delta\rho$ となる。背景保存式から $\zeta-\mathcal R=\Delta\rho/[3(\bar\rho+\bar p)]$ である。従って正則な $\Phi$ と非零のエンタルピーのもとで、十分長い波長では両者の差は勾配の二次に抑えられる。

:::example[物質期の大尺度]
平坦な物質優勢なら $4\pi Ga^2\bar\rho=3\mathcal H^2/2$ である。従って $\zeta-\mathcal R=-(2/9)(k/\mathcal H)^2\Phi$ である。$k/\mathcal H=0.1$ なら差は $-0.00222\Phi$ であり、二つの曲率は近い。
:::

### Lie 微分による変換の導出

背景と摂動を比較する座標の同定の変更を :term[ゲージ変換]{en="gauge transformation" yomi="げーじへんかん"} と呼ぶ。$\xi^\mu=(T,\partial^iL+L_T^i)$、$\partial_iL_T^i=0$ とする。変換後の座標点で旧座標を $\tilde x-\xi$ と評価し、二つの添字に $\partial x/\partial\tilde x$ を掛ければ $\widetilde{\delta g}_{\mu\nu}=\delta g_{\mu\nu}-(\mathcal L_\xi\bar g)_{\mu\nu}$ を得る。直接計算すると
$$
\begin{aligned}
a^{-2}(\mathcal L_\xi\bar g)_{00}&=-2(\mathcal HT+T'),\\
a^{-2}(\mathcal L_\xi\bar g)_{0i}&=\partial_i(L'-T)+(L_T)_i',\\
a^{-2}(\mathcal L_\xi\bar g)_{ij}&=2\mathcal HT\delta_{ij}
+2\partial_i\partial_jL+\partial_i(L_T)_j+\partial_j(L_T)_i .
\end{aligned}
$$
$0i$ では $\bar g_{00}\partial_i\xi^0$ と $\bar g_{ij}\partial_0\xi^j$ が逆符号で入る。これにより $\tilde S_i=S_i-(L_T)_i'$、$\tilde F_i=F_i-(L_T)_i$ であり、$V_i=S_i-F_i'$ は不変である。また空間成分の右辺は TT 射影で消えるため、$h_{ij}$ も不変である。

物質の背景は $\bar T^0{}_0=-\bar\rho$、$\bar T^j{}_i=\bar p\delta^j_i$ であり、$(\mathcal L_\xi\bar T)^0{}_i=-(\bar\rho+\bar p)\partial_iT$ を得る。従って $\tilde v_i=v_i+\partial_iT$ である。ここで $v_i$ は $\delta T^0{}_i=(\bar\rho+\bar p)v_i$ で定義する。一般スカラーゲージで $v_i=\partial_i v$ なら Newton ゲージでのポテンシャルは $v_N=v-\sigma$、密度は $\delta\rho_N=\delta\rho+\bar\rho'\sigma$ である。これらは $T=-\sigma$ を代入して得られ、任意の $T$ に対して不変である。

### 同期ゲージと残留変換

:term[同期ゲージ]{en="synchronous gauge" yomi="どうきげーじ"} は $A=B=0$ を課す。任意のゲージから移るには $T'+\mathcal H T=A$、$L'=B+T$ を解く。前者に積分因子 $a$ を掛けると $T=a^{-1}[\int aA\dd\eta+T_0(\vb x)]$ である。後者にも空間だけの積分定数 $L_0(\vb x)$ が残る。同期ゲージ内では残留自由度は $T=T_0/a$、$L=L_0+T_0\int\dd\eta/a$ となる。従って同期ゲージ内で異なる見かけの密度を得ても、物理的な成長が変わったと判断する前に残留変換を確かめる。

圧力なしの流体は同期座標では測地線に沿い、線形速度が $a^{-1}$ に減る。従って $T_0$ の選択により $\tilde\theta=\theta-k^2T_0/a=0$ と固定できる。これはその流体と共動する同期座標である。空間だけのずれ $L_0$ を指定するには、空間計量の初期条件を用いる。

:::example[座標だけから現れる物質密度]
無摂動の物質背景に残留変換 $T=T_0/a$ を掛けると $\delta=3\mathcal H T$ となる。$a\propto\eta^2$ なら見かけの密度は $\eta^{-3}$ に比例する。しかし $\sigma=T$ のため $\delta\rho_N=\delta\rho+\bar\rho'\sigma=0$ である。座標に依存しない密度は元の背景と一致する。
:::

### 空間平坦面と一様膨張面

$\mathcal H\ne0$ とする。$T=-C/\mathcal H$ で $\tilde C=0$ にでき、$L=E$ も指定すれば空間計量のスカラー摂動が全部消える。この面での密度は $\delta\rho_{\rm flat}=\delta\rho+\bar\rho'C/\mathcal H$ である。従って $\zeta=-\mathcal H\delta\rho_{\rm flat}/\bar\rho'$ となる。平坦面で測った密度が正で、膨張により $\bar\rho'<0$ なら、一様密度面における曲率振幅 $\zeta$ は正である。例えば圧力なしの物質で $\delta_{\rm flat}=3\times10^{-5}$ なら $\zeta=10^{-5}$ である。この規格化はどの時間座標から平坦面へ移ったかに依らない。

一様膨張面を選ぶには、どの膨張を指すかを先に決める。時刻面の法線の膨張は面そのものの取り方に依存する。一方、指定した流体の四元速度の $\mathcal K=\nabla_\mu u^\mu$ は時空のスカラーであり、一次変換は $\delta\mathcal K\mapsto\delta\mathcal K-\bar{\mathcal K}'T$ である。$\bar{\mathcal K}'\ne0$ なら $T=\delta\mathcal K/\bar{\mathcal K}'$ でこの流体の一様膨張面へ移れる。Newton ゲージでは $\bar{\mathcal K}=3\mathcal H/a$、$\delta\mathcal K=(\theta-3\Phi'-3\mathcal H\Psi)/a$ なので、同じ時刻面の法線と流体の膨張の差は $\theta/a$ である。この差が小さい大尺度では面の違いは勾配で抑えられる。定数の膨張率を持つ真空背景では $\bar{\mathcal K}'=0$ となるので、同じ方法で時刻を指定できない。

## Einstein 方程式の摂動

Newton ゲージで Fourier の規約を $e^{i\vb k\cdot\vb x}$ とする。$\theta=i k_i v^i$ は速度の発散であり、$\sigma_f$ は $k^i k_j/k^2-\delta^i_j/3$ をかけた跡なし応力の負を $\bar\rho+\bar p$ で割ったものとする。計量の空間曲率の摂動 $4a^{-2}\nabla^2\Phi$ と体積膨張の摂動 $-3(\Phi'+\mathcal H\Psi)/a$ を Hamiltonian 制約に入れると $4\nabla^2\Phi-12\mathcal H(\Phi'+\mathcal H\Psi)=16\pi Ga^2\delta\rho$ である。運動量制約は膨張摂動の空間勾配と $T^0{}_i=(\bar\rho+\bar p)v_i$ を結ぶ。空間の跡なし成分の二階微分を射影すると $k^2(\Phi-\Psi)$ が応力を測る。まとめると

$$
\begin{aligned}k^2\Phi+3\mathcal H(\Phi'+\mathcal H\Psi)&=-4\pi Ga^2\delta\rho,\\ k^2(\Phi'+\mathcal H\Psi)&=4\pi Ga^2(\bar\rho+\bar p)\theta,\\ k^2(\Phi-\Psi)&=12\pi Ga^2(\bar\rho+\bar p)\sigma_f .\end{aligned} \label{eq-cosmo-einstein-scalar}
$$

テンソルの波動は空間の横・跡なし成分から得る。

$$
h_{ij}''+2\mathcal Hh_{ij}'+k^2h_{ij}=16\pi Ga^2\Pi_{ij}^{\rm TT}. \label{eq-cosmo-tensor-wave}
$$

$\Pi^{\rm TT}$ は応力の横・跡なし成分である。この形の幾何学的な計算には :eqref[eq-gr-einstein-field] を一次に展開する。

:::example[亜地平線の物質]
時間項が小さいと過密 $\delta\rho>0$ は負の $\Phi$ を作る。
:::


:::example[応力のない流体]
非零の $k$ で第三式の右辺が零なら $\Phi=\Psi$ になる。
:::

### 空間の跡とポテンシャルの伝播

空間の跡を取ると、時間微分を含む Einstein 方程式が残る。Newton ゲージでは混合成分の跡の摂動に $a^2/6$ を掛けると、接続の時間部分が $\Phi''+\mathcal H(\Psi'+2\Phi')+(2\mathcal H'+\mathcal H^2)\Psi$ を、空間の微分が $\nabla^2(\Psi-\Phi)/3$ を与える。物質側は $\delta T^i{}_i=3\delta p$ なので
$$
\Phi''+\mathcal H(\Psi'+2\Phi')+(2\mathcal H'+\mathcal H^2)\Psi+\frac{k^2}{3}(\Phi-\Psi)=4\pi Ga^2\delta p .
\label{eq-cosmo-einstein-trace}
$$
異方性応力を持たない一成分で、$w=\bar p/\bar\rho$ が一定かつ $\delta p=w\delta\rho$ とする。:eqref[eq-cosmo-einstein-scalar] の第一式を代入して
$$
\Phi''+3(1+w)\mathcal H\Phi'+[2\mathcal H'+(1+3w)\mathcal H^2+wk^2]\Phi=0 .
$$
平坦で $w$ 一定の単成分なら、背景保存式で $2\mathcal H'=-(1+3w)\mathcal H^2$ となる。従って音速と膨張摩擦だけがポテンシャルの変化を決める。

:::example[物質期の一定モード]
$w=0$、$\mathcal H=2/\eta$ を用いると $(\eta^6\Phi')'=0$ である。その解は $\Phi=A+B\eta^{-5}$ であり、成長する密度摂動に対応するポテンシャルは時間一定となる。
:::

:::example[放射期のポテンシャルの減衰]
$w=1/3$、$\mathcal H=1/\eta$ とし $y=k\eta/\sqrt3$ とすると、方程式は $\Phi_{yy}+4\Phi_y/y+\Phi=0$ になる。$\Phi_i$ を $y=0$ での値とする正則解は $3\Phi_i(\sin y-y\cos y)/y^3$ である。分子の三次の項が $y^3/3$ なので初期値を再現し、$y\gg1$ で包絡が $y^{-2}$ に減る。
:::

### 自由重力波が地平線に入るとき

応力のないテンソル波で放射優勢なら $a\propto\eta$、$\mathcal H=1/\eta$ である。$x=k\eta$ を使い :eqref[eq-cosmo-tensor-wave] の左辺を $\partial_x^2h+2\partial_xh/x+h=0$ とする。$u=xh$ に置き換えると $u''+u=0$ であり、中心で有界かつ初期値 $h_i$ の解は $h/h_i=\sin x/x$ である。物質優勢では摩擦が $4/x$ となり、正則解 $3(\sin x-x\cos x)/x^3$ を二回微分して確認できる。大きな $x$ で包絡はそれぞれ $x^{-1},x^{-2}$ であり、いずれも $a^{-1}$ となる。

:::example[最初の零点の違い]
放射期では $x=\pi$ が最初の零点である。物質期では $x=\pi$ において $h/h_i=3/\pi^2\simeq0.304$ で、零点は $\tan x=x$ の最初の正の非零根へ移る。従って地平線に入った後の振幅だけでなく位相にも背景が現れる。
:::

### 制約から求める物質優勢期の解析解

物質のみの平坦な背景で $\mathcal H=2/\eta$、$\sigma_f=0$ とする。減衰モードを除けば $\Phi=\Psi=\Phi_i$ は一定である。速度の式 $\theta'+2\theta/\eta=k^2\Phi_i$ に $\eta^2$ を掛けて積分し、中心で正則な解を選ぶと $\theta=k^2\eta\Phi_i/3$ である。連続式 $\delta'=-\theta$ を積分すると $\delta=\delta_i-k^2\eta^2\Phi_i/6$ である。時間制約では $4\pi Ga^2\bar\rho=6/\eta^2$ だから、$k\eta\ll1$ の定数は $\delta_i=-2\Phi_i$ となる。この値は任意の積分定数ではなく、Einstein 制約を満たすために決まる。

共動ゲージの密度差は $\Delta\rho/\bar\rho=\delta+3\mathcal H\theta/k^2$ である。上の解を入れると定数項が相殺され、$\Delta\rho/\bar\rho=-(k\eta)^2\Phi_i/6$ である。従って $\eta^2\propto a$ の成長則は共動密度で全波数に見えるが、Newton ゲージの密度では長波長に定数部分がある。$\mathcal R=-\Phi_i-\mathcal H\theta/k^2=-5\Phi_i/3$ も一定であり、通常の Sachs–Wolfe 項の正規化に使える。

:::example[二つのゲージで異なる密度を比較する]
$k\eta=1$、$\Phi_i=-6\times10^{-6}$ なら $\Delta\rho/\bar\rho=10^{-6}$、Newton ゲージの $\delta=1.3\times10^{-5}$ である。差の $1.2\times10^{-5}$ は密度を評価する時間面の違いによる。Einstein 制約を解かずに一方をもう一方の Poisson 方程式へ入れると、この尺度では大きな誤差になる。
:::

### Hamiltonian 制約と運動量制約の成分計算

空間計量 $\gamma_{ij}=a^2(1-2\Phi)\delta_{ij}$、lapse $N=a(1+\Psi)$、シフト零とする。法線を $n^\mu=(N^{-1},0,0,0)$ とし、外部曲率は $K_{ij}=\gamma'_{ij}/(2N)$ の符号で定義する。空間接続の線形項と空間 Ricci テンソルは順に
$$
\begin{aligned}
{}^{(3)}\Gamma^i{}_{jk}&=-\delta^i_j\partial_k\Phi-\delta^i_k\partial_j\Phi
+\delta_{jk}\partial^i\Phi,\\
{}^{(3)}R_{ij}&=\partial_l{}^{(3)}\Gamma^l{}_{ij}
-\partial_j{}^{(3)}\Gamma^l{}_{il}
=\partial_i\partial_j\Phi+\delta_{ij}\nabla^2\Phi .
\end{aligned}
$$
跡は $4a^{-2}\nabla^2\Phi$ である。$K^i{}_j=a^{-1}[\mathcal H-\Phi'-\mathcal H\Psi]\delta^i_j$ だから、$K^2-K^i{}_jK^j{}_i$ は中括弧の係数の二乗の六倍に $a^{-2}$ を掛けたものである。一次項は $-12a^{-2}\mathcal H(\Phi'+\mathcal H\Psi)$ となる。Gauss の関係 $2G_{nn}=R^{(3)}+K^2-K^i{}_jK^j{}_i$ に入れ、背景式を引くと第一制約を得る。

Codazzi の関係では $G_{ni}=D_jK^j{}_i-D_iK=2a^{-1}\partial_i(\Phi'+\mathcal H\Psi)$ となる。一方、$T_{ni}=-a(\bar\rho+\bar p)v_i$ である。従って $ik^i$ を掛けると、$-k^2(\Phi'+\mathcal H\Psi)=-4\pi Ga^2(\bar\rho+\bar p)\theta$ となる。法線射影には一つの下付き時間添字があるため、混合成分 $T^0{}_i$ と符号が異なる。

空間の跡なし成分では、空間 Ricci テンソルと $-N^{-1}D^iD_jN$ の Hessian の差だけが残り、$(\delta G^i{}_j)_{\rm TF}=a^{-2}(\partial^i\partial_j-\delta^i_j\nabla^2/3)(\Phi-\Psi)$ になる。$Q^i{}_j=\hat k^i\hat k_j-\delta^i_j/3$ とおくと $Q^i{}_jQ^j{}_i=2/3$ である。Einstein 方程式を $Q^j{}_i$ と縮約すれば、左辺は $-(2/3)k^2(\Phi-\Psi)/a^2$、右辺は $-8\pi G(\bar\rho+\bar p)\sigma_f$ となる。これで :eqref[eq-cosmo-einstein-scalar] の第三式の係数が定まる。

### テンソルの Ricci テンソルと波の摩擦

$\delta g_{ij}=a^2h_{ij}$、他の計量摂動を零とする。$\tilde g_{\mu\nu}=g_{\mu\nu}/a^2$ の線形 Ricci テンソルは $\delta\tilde R_{ij}=(h_{ij}''-\nabla^2h_{ij})/2$ である。実際、$\delta R_{\mu\nu}=(\partial_\alpha\partial_\mu h^\alpha{}_\nu+\partial_\alpha\partial_\nu h^\alpha{}_\mu-\Box h_{\mu\nu}-\partial_\mu\partial_\nu h^\alpha{}_\alpha)/2$ のうち横・跡なし条件で $-\Box$ だけが残る。共形変換した $R_{ij}$ では $-2\tilde\nabla_i\tilde\nabla_j\log a=2\mathcal H\tilde\Gamma^0{}_{ij}=\mathcal H h_{ij}'$ である。他の共形項は背景の係数を $\tilde g_{ij}$ に掛けた形を持つ。混合成分 $R^i{}_j$ では逆計量がその部分の $h_{ij}$ を打ち消し、跡を引く補正は TT 射影で消える。従って :eqref[eq-cosmo-tensor-wave] が得られる。

応力なしで $u_{ij}=ah_{ij}$ を導入すると $u_{ij}''+(k^2-a''/a)u_{ij}=0$ になる。$k^2\gg|a''/a|$ で係数が緩やかなら短時間の振動は平面波に近く、$h_{ij}$ の振幅は $a^{-1}$ に比例する。一方、勾配を捨てた式 $(a^2h_{ij}')'=0$ は定数解と $\int\dd\eta/a^2$ の解を持つ。$a\propto\eta^p$ の場合、第二解の時間依存部は $p\ne1/2$ で $\eta^{1-2p}$、$p=1/2$ では $\log\eta$ である。従って膨張だけから全ての背景で第二解が減衰するとは結論できない。

### ベクトル重力と渦度

$F_i=0$ を選び、シフトを $S_i=V_i$ とする。横の運動量密度を $q_i^{\rm V}=(\delta T^0{}_i)_{\rm T}$、物質の横ベクトル応力を $(\delta T^i{}_j)_{\rm V}=\partial^iU_j+\partial_j U^i$ と定義する。外部曲率のベクトル部分は $-(\partial^iV_j+\partial_jV^i)/(2a)$ である。横の運動量制約では $D_jK^j{}_i-D_iK=-\nabla^2V_i/(2a)$ なので $T_{ni}=-a q_i^{\rm V}$ と合わせて $\nabla^2 V_i=16\pi Ga^2 q_i^{\rm V}$ となる。空間の混合成分のベクトル部分は $(K^i{}_j)'/a$ と $3\mathcal HK^i{}_j/a$ の摂動から $-[\partial^i(V_j'+2\mathcal HV_j)+\partial_j(V^{i\prime}+2\mathcal HV^i)]/(2a^2)$ である。対称化した勾配は非零波数の横ベクトル上で単射なので
$$
-k^2V_i=16\pi Ga^2q_i^{\rm V},\qquad
V_i'+2\mathcal HV_i=-16\pi Ga^2U_i
$$
を得る。第一式を微分して第二式を代入すると $q_i^{{\rm V}\prime}+4\mathcal Hq_i^{\rm V}=k^2U_i$ となる。従って完全流体のベクトル重力は $a^{-2}$、運動量密度は $a^{-4}$ に減る。定数 $w$ の完全流体なら速度の横成分は $a^{3w-1}$ に比例する。放射の横速度が一定でも、そのベクトル重力は膨張で減る。

:::example[外部応力によるベクトル摂動]
$a^4U_i=U_{0i}$ が一定なら $(a^2V_i)'=-16\pi G U_{0i}$ である。初期値が零の場合、時間 $\Delta\eta$ の間に $V_i=-16\pi G U_{0i}\Delta\eta/a^2$ が生じる。運動量の制約から、それを支える $q_i^{\rm V}$ も決まる。
:::

### 空間の跡を接続から計算する

計量を $g_{\mu\nu}=a^2\tilde g_{\mu\nu}$ とし、$\tilde g_{\mu\nu}$ を Minkowski 計量の摂動として扱う。添字 $\tilde{}$ の接続で、一次に現れる独立な成分は
$$
\begin{aligned}
\tilde\Gamma^0{}_{00}&=\Psi',&
\tilde\Gamma^0{}_{0i}&=\partial_i\Psi,&
\tilde\Gamma^0{}_{ij}&=-\Phi'\delta_{ij},\\
\tilde\Gamma^i{}_{00}&=\partial^i\Psi,&
\tilde\Gamma^i{}_{0j}&=-\Phi'\delta^i_j,&
\tilde\Gamma^i{}_{jk}&=-\delta^i_j\partial_k\Phi-\delta^i_k\partial_j\Phi
+\delta_{jk}\partial^i\Phi .
\end{aligned}
$$
$\tilde R_{\mu\nu}=\partial_\alpha\tilde\Gamma^\alpha{}_{\mu\nu}-\partial_\nu\tilde\Gamma^\alpha{}_{\mu\alpha}$ を使うと
$$
\begin{aligned}
\tilde R_{00}&=\nabla^2\Psi+3\Phi'',\\
\tilde R_{0i}&=2\partial_i\Phi',\\
\tilde R_{ij}&=\partial_i\partial_j(\Phi-\Psi)
+\delta_{ij}(\nabla^2\Phi-\Phi'').
\end{aligned}
$$
特に $\tilde R_{00}$ には $\Psi''$ が残らず、lapse が独立な波動変数にならないことが分かる。次に $\lambda=\log a$ と書き、四次元の共形変換
$$
R_{\mu\nu}=\tilde R_{\mu\nu}-2\tilde\nabla_\mu\partial_\nu\lambda
-\tilde g_{\mu\nu}\tilde\Box\lambda
+2\partial_\mu\lambda\partial_\nu\lambda
-2\tilde g_{\mu\nu}(\tilde\nabla\lambda)^2
$$
を使う。$\tilde\Box\lambda=-\mathcal H'+2\Psi\mathcal H'+\mathcal H(\Psi'+3\Phi')$、$(\tilde\nabla\lambda)^2=-(1-2\Psi)\mathcal H^2$ である。空間の混合成分を取ると、上の $\tilde R_{ij}$ に加えて対角方向に
$$
a^{-2}\big[\mathcal H'+2\mathcal H^2
-2(\mathcal H'+2\mathcal H^2)\Psi-\mathcal H(\Psi'+5\Phi')\big]\delta^i_j
$$
が入る。時間成分は $R_{00}=-3\mathcal H'+\nabla^2\Psi+3\Phi''+3\mathcal H(\Psi'+\Phi')$ であり、$g^{00}=-a^{-2}(1-2\Psi)$ を掛けて $R^0{}_0$ にする。最後に $G^i{}_i=-(R^i{}_i+3R^0{}_0)/2$ と整理すれば、一次の跡は
$$
(\delta G^i{}_i)=\frac6{a^2}\left[
\Phi''+\mathcal H(\Psi'+2\Phi')+(2\mathcal H'+\mathcal H^2)\Psi
+\frac{\nabla^2(\Psi-\Phi)}3\right]
$$
になる。この結果に $\delta T^i{}_i=3\delta p$ を入れ、Fourier 変換すれば :eqref[eq-cosmo-einstein-trace] である。混合成分を使う前に背景の $R_{ij}$ を引くと、逆計量の一次摂動を落とすことがある。ここでは混合成分を作ってから背景を引いた。

### 同期ゲージで制約を確認する

同期ゲージのスカラー空間計量を $X_{ij}=\hat k_i\hat k_j h_s+6(\hat k_i\hat k_j-\delta_{ij}/3)\eta_s$ と書く。ここで $h_s$ は空間計量の跡、$\eta_s$ は曲率の係数であり、共形時刻 $\eta$ と区別する。:eqref[eq-cosmo-svt] に合わせれば $C=\eta_s$、$E=-(h_s+6\eta_s)/(2k^2)$ である。$B=0$ なので $\alpha_s=-E'=(h_s'+6\eta_s')/(2k^2)$ が Newton ゲージへ移る時間変換の負である。従って
$$
\Psi=\alpha_s'+\mathcal H\alpha_s,\quad
\Phi=\eta_s-\mathcal H\alpha_s,\quad
\delta\rho_N=\delta\rho_s+\bar\rho'\alpha_s,\quad
\theta_N=\theta_s+k^2\alpha_s .
$$
平坦背景では $\mathcal H^2-\mathcal H'=4\pi Ga^2(\bar\rho+\bar p)$ である。Newton ゲージの運動量制約に両側の変換を入れると、$\alpha_s$ に比例する項がこの背景式で打ち消し合う。時間制約でも $\bar\rho'=-3\mathcal H(\bar\rho+\bar p)$ で同じ相殺が起きる。残りは
$$
k^2\eta_s-\frac12\mathcal H h_s'=-4\pi Ga^2\delta\rho_s,\qquad
k^2\eta_s'=4\pi Ga^2(\bar\rho+\bar p)\theta_s
$$
である。同期ゲージを使う計算と Newton ゲージを使う計算の出力を比べる際は、同じ物理変数へ変換してから差を取る。

:::example[冷たい物質に共動する同期計量]
純物質期で一定の $\Phi=\Psi$ を持つ正則解を取る。$\theta_N=k^2\eta\Phi/3$ であるから、$\theta_s=0$ の同期ゲージへの変換で $\alpha_s=\eta\Phi/3$ となる。$\mathcal H=2/\eta$ を代入すれば $\eta_s=5\Phi/3$ で一定である。$\delta_s=\delta_N-\bar\rho'\alpha_s/\bar\rho=-(k\eta)^2\Phi/6$ であり、同期ゲージで成長する密度は共動密度と一致する。$h_s'=2k^2\alpha_s$ を積分すると $h_s=(k\eta)^2\Phi/3+$ 空間ゲージの残余である。この残余は密度成長に寄与しない。
:::

### 制約は時間発展でどう伝わるか

制約二本を各時刻で独立な発展式として積分すると、同じ自由度に矛盾する条件を課すことがある。これを調べるため $A_G=4\pi Ga^2$、$h=\bar\rho+\bar p$、$J=\Phi'+\mathcal H\Psi$ と置き、制約残差を
$$
F_0=k^2\Phi+3\mathcal H J+A_G\delta\rho,\quad
F_1=k^2J-A_Gh\theta,\quad
F_2=k^2(\Phi-\Psi)-3A_Gh\sigma_f
$$
と定義する。背景で $\mathcal H^2-\mathcal H'=A_Gh$ である。空間の跡から
$$
J'=A_G\delta p-2\mathcal HJ+A_Gh\Psi-\frac{k^2}3(\Phi-\Psi)
$$
を得る。エネルギー保存 $\delta\rho'=-3\mathcal H(\delta\rho+\delta p)+h(3\Phi'-\theta)$ を $F_0'$ へ代入する。$A_G'=2\mathcal HA_G$ のため $\delta\rho$ に掛かる係数は $-\mathcal HA_G$、$\delta p$ の項は相殺する。背景式で $J$ と $\Phi'$ の係数も整理すると $F_0'=-\mathcal HF_0+F_1$ となる。また $(h\theta)'+4\mathcal Hh\theta=k^2[\delta p+h(\Psi-\sigma_f)]$ を $F_1'$ へ代入すると、重力と圧力の項が打ち消し合い
$$
F_0'=-\mathcal HF_0+F_1,\qquad
F_1'=-2\mathcal HF_1-\frac{k^2}3F_2
$$
となる。従って跡なしの制約を保ち、初期の $F_0,F_1$ が零なら、連立一次方程式の一意性から両者は全時刻で零である。これは Bianchi 恒等式と物質保存の関係を成分で見たものである。数値解で残差が増える場合、物質の保存、背景式、跡なし応力の規格化が一貫しているかを検査できる。

:::example[初期の制約誤差が残る場合]
$F_2=0$ を保つが、初期に $F_{1i}\ne0$ とする。$F_1=F_{1i}(a_i/a)^2$ なので、$(aF_0)'=aF_1=F_{1i}a_i^2/a$ である。初期に $F_{0i}=0$ でも後の $F_0$ は零でない。運動の時間積分が安定でも、初期データが制約に反していれば本来の解を得ない。
:::

### 空間曲率が小さい背景で何が変わるか

一定曲率 $K$ の三次元背景計量 $\bar\gamma_{ij}$ を用い、ここだけ $g_{ij}=a^2(1-2\Phi)\bar\gamma_{ij}$ とする。$R[\bar\gamma]=6K$ である。三次元の共形変換を一次まで使うと
$$
R^{(3)}=\frac1{a^2}(6K+12K\Phi+4\bar\nabla^2\Phi)
$$
となる。空間依存する共形係数 $e^{-2\Phi}$ により背景の $6K$ 自体が $12K\Phi$ の摂動を持つ点が、平坦の場合との違いである。運動量制約はスカラー勾配の形を保つので、時間制約と組み合わせると $(\bar\nabla^2+3K)\Phi=4\pi Ga^2\Delta\rho$ になる。Laplace 固有値 $\bar\nabla^2=-k^2$ で読めば、平坦の $k^2$ は $k^2-3K$ に替わる。これは Fourier 空間の全ての式を一律に置換するという意味ではない。テンソルやベクトルでは共変微分の順序の交換にも曲率が現れ、空間の境界条件によって固有モードが異なる。平坦の式との局所比較には少なくとも $|K|/k^2\ll1$ が必要である。

## 流体成分の摂動方程式

エネルギー保存の一次式は $\delta\rho'+3\mathcal H(\delta\rho+\delta p)=(\bar\rho+\bar p)(3\Phi'-\theta)$ である。これは局所膨張の一次変化 $\theta-3\Phi'$ を背景のエンタルピー密度に掛けた式である。運動量保存は $(\bar\rho+\bar p)v_i$ の共形時刻微分と $4\mathcal H$ の希釈、圧力勾配、重力、跡なし応力からなる。$\delta=\delta\rho/\bar\rho$、定数 $w=\bar p/\bar\rho$、$\delta p=w\delta\rho$ の場合は

$$
\delta'=-(1+w)(\theta-3\Phi'),\qquad \theta'=-\mathcal H(1-3w)\theta+k^2\left(\Psi+\frac{w}{1+w}\delta-\sigma_f\right). \label{eq-cosmo-fluid-perturbations}
$$

$w$ が変化する成分では背景の $\bar p'$ や静止系の音速を別に管理する。真空のように $\bar\rho+\bar p=0$ の成分に速度方程式を適用できない。

:::example[冷たい物質]
自身の重力を無視した試験成分では、外部の $\Psi=0$ で初期速度があっても $\theta\propto a^{-1}$ に減衰する。
:::


:::example[理想放射]
重力と異方性を無視すると $\delta''+(k^2/3)\delta=0$ になり、音速は $1/\sqrt3$ である。
:::

### 時間変化する状態方程式

交換のない成分の運動量密度に対し、$T^0{}_i=(\bar\rho+\bar p)v_i$ の微分の項をそのまま残す。$\bar\rho'/\bar\rho=-3\mathcal H(1+w)$ を使って $(\bar\rho+\bar p)'/(\bar\rho+\bar p)=-3\mathcal H(1+w)+w'/(1+w)$ を得る。ここから
$$
\begin{aligned}
\delta'&=-(1+w)(\theta-3\Phi')-3\mathcal H(\delta p/\bar\rho-w\delta),\\
\theta'&=-\left[\mathcal H(1-3w)+\frac{w'}{1+w}\right]\theta+
k^2\left[\Psi+\frac{\delta p}{\bar\rho(1+w)}-\sigma_f\right]
\end{aligned}
$$
となる。$w=-1$ で表式が特異なのは速度を定義する分母が零になるためである。背景に $w(a)$ を指定するだけでは $\delta p$ と $\sigma_f$ は分からない。構成する場や粒子の性質が必要である。

:::example[一定の背景比でも圧力摂動が別にある場合]
$w=0$ と $w'=0$ を設定しても、$\delta p/\bar\rho>0$ なら速度の式に $k^2\delta p/\bar\rho$ が残る。正の圧力摂動で流れが外向きに発散して、連続式の $-\theta$ を通じて過密が弱まる。
:::

### 静止系の音速から Newton ゲージへ戻す

成分の静止系で圧力と密度の揺らぎを $\delta p_{\rm rest}=c_s^2\delta\rho_{\rm rest}$ と閉じる。時刻変換 $T=\theta/k^2$ で $\delta p_{\rm rest}=\delta p-\bar p'T$、$\delta\rho_{\rm rest}=\delta\rho-\bar\rho'T$ なので
$$
\delta p=c_s^2\delta\rho+3\mathcal H(\bar\rho+\bar p)(c_s^2-c_a^2)\frac{\theta}{k^2},
\qquad c_a^2=\bar p'/\bar\rho'.
$$
従って非断熱圧力は $(c_s^2-c_a^2)\Delta\rho$ である。$w=\mathrm{const.}$ でも $c_s^2\ne w$ は可能であり、$w$ の値だけでは伝播速度を指定できない。

:::example[一定圧力比の場]
$w=-1/2$ を背景で保つ正準場の静止系では $c_s^2=1$ である。$c_a^2=w$ なので $\delta p_{\rm nad}=(3/2)\Delta\rho$ である。$\delta p=w\delta\rho$ と置くとこの圧力応答を失う。
:::

### 共変保存則の直接展開

四元速度を $u^\mu=a^{-1}(1-\Psi,v^i)$ とする。$\sqrt{-g}=a^4(1+\Psi-3\Phi)$ なので、$\nabla_\mu u^\mu=(\sqrt{-g})^{-1}\partial_\mu(\sqrt{-g}u^\mu)$ に代入すると
$$
\nabla_\mu u^\mu=a^{-1}[3\mathcal H(1-\Psi)-3\Phi'+\theta].
$$
エネルギー方程式は $u^\mu\partial_\mu\rho+(\rho+p)\nabla_\mu u^\mu=0$ である。跡なし応力の仕事は速度勾配との積なので等方背景では一次に現れない。$u^\mu\partial_\mu\rho=a^{-1}(\bar\rho'+\delta\rho'-\Psi\bar\rho')$ を用いると、$\Psi$ の係数は $-\bar\rho'-3\mathcal H(\bar\rho+\bar p)=0$ になる。残りを整理するとエネルギーの一次式を得る。

運動量式では $u^\mu\nabla_\mu u_i=v_i'+\mathcal H v_i+\partial_i\Psi$ である。圧力の空間射影を計算すると $\partial_i\delta p+\bar p'v_i$ が現れる。応力の発散を加えると
$$
(\bar\rho+\bar p)(v_i'+\mathcal H v_i+\partial_i\Psi)
+\bar p'v_i+\partial_i\delta p+\partial_j\Pi^j{}_i=0 .
$$
$i k^i$ を掛けると圧力とポテンシャルの勾配には $-k^2$ が付き、$\partial_i\partial_j\Pi^{ij}=k^2(\bar\rho+\bar p)\sigma_f$ となる。さらに $\mathcal H+\bar p'/(\bar\rho+\bar p)=\mathcal H(1-3w)+w'/(1+w)$ なので、先ほどの速度式が得られる。

### 重力を与えた場合の密度応答

一定 $w=c_s^2>0$ の完全流体の自己重力を無視し、$\Phi,\Psi$ を外から与える。連続式を微分して速度式を代入すると
$$
\delta''+\mathcal H(1-3w)\delta'+wk^2\delta
=(1+w)[3\Phi''+3\mathcal H(1-3w)\Phi'-k^2\Psi].
$$
静的ポテンシャルで $k\ne0$ なら、定常特解は $\delta=-(1+w)\Psi/w$ である。これは静水圧の釣合い $\partial_i\delta p=-(\bar\rho+\bar p)\partial_i\Psi$ に一致する。一般解には初期値で決まる同次解も加わる。

:::example[放射の静水圧温度]
$w=1/3$ では $\delta=-4\Psi$、温度摂動は $\Theta_0=\delta/4=-\Psi$ である。従って $\Theta_0+\Psi$ は零となる。井戸の中の放射温度は高いが、そこから出た光は赤方偏移する。
:::

### 成分間のエネルギー交換

全体で運動量・エネルギーを交換する複数成分を取る。成分 $i$ の四元ベクトル $Q_i^\nu$ を $\nabla_\mu T_i^{\mu\nu}=Q_i^\nu$、$\sum_i Q_i^\nu=0$ で定義する。成分の静止系における固有時当たりのエネルギー受取を $Q_i=-u_{i\nu}Q_i^\nu$ とする。背景では全成分が同じ四元速度を持つので $\bar\rho_i'+3\mathcal H h_i=a\bar Q_i$、$h_i=\bar\rho_i+\bar p_i$ である。

エネルギー保存の先ほどの展開を繰り返すと、lapse の項は $-\Psi(\bar\rho_i'+3\mathcal Hh_i)$ となる。背景の右辺が零でないため、これを移項して
$$
\delta\rho_i'+3\mathcal H(\delta\rho_i+\delta p_i)-h_i(3\Phi'-\theta_i)
=a(\delta Q_i+\Psi\bar Q_i)
$$
を得る。$Q_i$ の測定が固有時当たりだから、共形時へ換算すると $\Psi\bar Q_i$ の補正が付く。エネルギー受取の合計は背景で零であり、成分間の速度差と運動量交換の積は二次なので、一次でも $\sum_i\delta Q_i=0$ である。全式を足すと孤立した系の保存式が戻る一方、各成分の $\zeta_i$ は個別に保存するとは限らない。この場合にも $\zeta_i=-C-\mathcal H\delta\rho_i/\bar\rho_i'$ は不変であるが、分母を $-3\mathcal Hh_i$ と置くことはできない。

:::example[一定寿命で崩壊する物質]
圧力なしの成分 $m$ が固有時当たり一定の率 $\Gamma$ で放射へ崩壊するとする。$Q_m=-\Gamma\rho_m$ より、Newton ゲージで $\delta Q_m=-\Gamma\bar\rho_m\delta_m$ である。数密度の背景減衰を用いて分数摂動の式に直すと $\delta_m'=-\theta_m+3\Phi'-a\Gamma\Psi$ が残る。$\delta_m$ に直接 $-a\Gamma\delta_m$ を足すと、背景と摂動の同じ割合の減少を二度数える。この式には崩壊が成分の静止系で等方になるという条件も用いた。
:::

### 一定音速の厳密解と波の領域

応力零の一成分で $p=w\rho$、$w>0$ 一定、平坦背景を考える。$\mathcal H=2/[(1+3w)\eta]$ と $x=\sqrt w\,k\eta$ をポテンシャルの発展式に入れると
$$
\Phi_{xx}+\frac{6(1+w)}{(1+3w)x}\Phi_x+\Phi=0 .
$$
$m=[6(1+w)/(1+3w)-1]/2=(5+3w)/[2(1+3w)]$、$\Phi=x^{-m}F(x)$ とおく。微分項を整理すると $F_{xx}+F_x/x+(1-m^2/x^2)F=0$ となる。正則性を課すと $Y_m$ の枝は除かれ、$J_m(x)\sim(x/2)^m/\Gamma(m+1)$ から
$$
\Phi(x)=\Phi_i\,2^m\Gamma(m+1)x^{-m}J_m(x)
$$
を得る。$w=1/3$ で $m=3/2$ となり、球 Bessel の放射解へ戻る。大きい $x$ では $J_m$ の包絡が $x^{-1/2}$ なので、ポテンシャルは $x^{-(m+1/2)}$ の包絡で振動する。$w\to0$ と $x\gg1$ の極限は同じ固定時刻で交換できない。$x$ 自体が $\sqrt w$ を含むため、圧力が無視できる物質の解は $w=0$ の方程式から求める。

## 曲率摂動の保存

$c_a^2=\bar p'/\bar\rho'$、$\delta p_{\rm nad}=\delta p-c_a^2\delta\rho$ とする。$\delta p_{\rm nad}$ は時間変換で追加される項が相殺するので不変である。背景の保存則により Newton ゲージで $\zeta=-\Phi+\delta\rho/[3(\bar\rho+\bar p)]$ と書ける。時間微分して上のエネルギー式を入れると

$$
\zeta'=-\frac{\mathcal H}{\bar\rho+\bar p}\delta p_{\rm nad}-\frac{\theta}{3}. \label{eq-cosmo-zeta-conservation}
$$

$\theta=k^2$ に比例する正則な速度ポテンシャルなら長波長で最後の項が消える。断熱なら $\delta p_{\rm nad}=0$ なので、その極限で $\zeta$ は保存する。成分ごとの圧力が密度だけの関数でも、成分間の相対密度のずれが総系の非断熱圧力を作り得る。

:::example[断熱な物質]
$\delta p=0=c_a^2$ なので $\theta\to0$ で曲率が保存する。
:::


:::example[非断熱源]
一定区間に $\delta p_{\rm nad}/(\bar\rho+\bar p)=s_0$、$\theta=0$ なら $\Delta\zeta=-s_0\log(a_f/a_i)$ である。
:::

### 相対エントロピーが曲率を変える

交換のない成分 $i$ に対して $h_i=\bar\rho_i+\bar p_i$、$\zeta_i=-\Phi+\delta\rho_i/(3h_i)$ とする。全成分で $h=\sum_i h_i$ を使うと $\zeta=\sum_i h_i\zeta_i/h$ である。成分間の差 $S_{ij}=3(\zeta_i-\zeta_j)$ は時間の共通ずれが消えるゲージ不変量である。$c_i^2=\bar p_i'/\bar\rho_i'$ とすると $c_a^2=\sum_i h_i c_i^2/h$ だから
$$
\delta p_{\rm nad}=\sum_i\delta p_{{\rm nad},i}
+\frac1{2h}\sum_{ij}h_i h_j(c_i^2-c_j^2)S_{ij}.
$$
第二項は $\sum_i(c_i^2-c_a^2)\delta\rho_i$ に $c_a^2$ を代入し、$i,j$ を交換した式と平均して得る。各成分の固有の非断熱圧力が零でも、二成分の圧力応答が違えば総系の曲率は変わりうる。

:::example[物質と放射での転換]
長波長で $\zeta_m=10^{-5}$、$\zeta_r=0$ が個別に保存するとする。$\rho_m=\rho_r$ では $\zeta=3\zeta_m/7$ であり、物質が優勢になると $\zeta$ は $\zeta_m$ に近づく。全系の曲率が変わるのは、それぞれの保存値のエンタルピー重みが変化するためである。
:::

### 長波長極限をとる条件

有限区間 $[\eta_i,\eta_f]$ で $h=\bar\rho+\bar p$ が零から離れ、$|\delta p_{\rm nad}/h|\le k^2b(\eta)$、$|\theta|\le k^2c(\eta)$ とする。$|\mathcal H|b,c$ が積分可能なら :eqref[eq-cosmo-zeta-conservation] から
$$
|\zeta(\eta_f)-\zeta(\eta_i)|\le k^2\int_{\eta_i}^{\eta_f}(|\mathcal H|b+c/3)\dd\eta .
$$
従って固定した区間では $k\to0$ で保存する。$h$ が小さくなる背景や積分区間も伸びる極限では、この評価を再確認する。非零波数のゲージ分解をそのまま背景の零モードに適用してはいけない。

一様密度面で $E=0$ なら空間の体積要素は $a^3(1-3C)$ である。従って二つの一様密度面の間の局所的な対数膨張の差は、勾配を無視すれば $\zeta(\eta_f)-\zeta(\eta_i)$ である。圧力を密度だけの一価関数とし、局所流出入も無視すれば $\dd\rho/\dd N=-3(\rho+p(\rho))$ となる。ここで $N$ は局所体積の対数を三で割った量である。両面で初期密度と終点密度が等しければ、積分 $\int-\dd\rho/[3(\rho+p)]$ は全領域で同じであり、曲率は保存する。組成が独立なら、同じ密度で圧力が違うため、この推論は成立しない。

### 放射期と物質期の規格化

平坦な単成分背景、応力零、$w>-1$ の一定な長波長モードでは、第一制約から $\delta=-2\Phi$ となる。$\zeta=-\Phi+\delta/[3(1+w)]$ と組み合わせれば $\Phi=-3(1+w)\zeta/(5+3w)$ である。放射期に $\Phi=-2\zeta/3$、物質期に $\Phi=-3\zeta/5$ なので、両時代の定数モードの振幅比は $9/10$ である。断熱性で保存される量は $\zeta$ であり、ポテンシャルには状態方程式の違いが現れる。

:::example[成分間の相対摂動だけがある初期条件]
$\zeta_r=0$、$S_{mr}=S_0$、$y=\bar\rho_m/\bar\rho_r$ とすると $\zeta_m=S_0/3$ であり、$\zeta=yS_0/(3y+4)$ である。$y\propto a$ なので $\dd\zeta/\dd\log a=4yS_0/(3y+4)^2$ を得る。$c_m^2-c_r^2=-1/3$ を相対非断熱圧力の式に入れ、曲率保存式を使っても同じ結果が出る。
:::

### 共動曲率の時間微分

総系について $h=\bar\rho+\bar p\ne0$ とし、Newton ゲージの $\mathcal R=-\Phi-\mathcal H\theta/k^2$ を微分する。運動量保存式における $\theta$ の摩擦を $\mathcal H+\bar p'/h$ と書けば
$$
\begin{aligned}
\mathcal R'&=-(\Phi'+\mathcal H\Psi)
+(\mathcal H^2-\mathcal H')\frac{\theta}{k^2}
-\frac{\mathcal H}{h}\left(\delta p-\bar p'\frac{\theta}{k^2}\right)
+\mathcal H\sigma_f\\
&=-\frac{\mathcal H}{h}\delta p_{\rm com}+\mathcal H\sigma_f .
\end{aligned}
$$
第二行では有限 $k$ の運動量制約と $\mathcal H^2-\mathcal H'=4\pi Ga^2h$ を使った。$\delta p_{\rm com}$ は共動面での圧力である。単成分の閉じ方が $\delta p_{\rm com}=c_s^2\Delta\rho$ なら、共動 Poisson 制約から $\mathcal R'=\mathcal H c_s^2 k^2\Phi/(4\pi Ga^2h)+\mathcal H\sigma_f$ である。圧力と応力がともに零なら $\mathcal R'$ は全波数で零になる。$\zeta$ の保存式には速度発散があるので、有限 $k$ では二者の時間依存が異なり得る。

:::example[物質期の有限波数における曲率の差]
物質だけで $\Phi$ が一定なら、$\mathcal R=-5\Phi/3$ は一定である。一方 $\zeta=\mathcal R-(k\eta)^2\Phi/18$ となる。従って $\zeta'=-k^2\eta\Phi/9=-\theta/3$ である。断熱であっても勾配を保持すると一様密度面の曲率が変わるという点を、速度式と制約式で別々に確認できる。
:::

### 等密度期を通る長波長の重力ポテンシャル

交換のない完全放射と物質、共通の一定曲率 $\zeta$、異方性応力零を仮定する。$y=\bar\rho_m/\bar\rho_r\propto a$ とする。$k/\mathcal H$ の二次を捨てた時間制約から $\delta=-2(\Phi+y\partial_y\Phi)$ であり、$\zeta=-\Phi+\delta/[3(1+w)]$ に $w=1/[3(1+y)]$ を代入すると
$$
y\partial_y\Phi+\frac{5y+6}{2(1+y)}\Phi
=-\frac{3y+4}{2(1+y)}\zeta .
$$
積分因子は $y^3/\sqrt{1+y}$ である。従って正則解は
$$
\Phi(y)=-\zeta\frac{\sqrt{1+y}}{y^3}
\int_0^y\frac{u^2(3u+4)}{2(1+u)^{3/2}}\dd u .
$$
下端に定数を足すと $\sqrt{1+y}/y^3$ の発散解を加えるので、正則な枝では積分定数は零である。積分を $s=\sqrt{1+u}$ に変換すると被積分部分は $(s^2-1)^2(3s^2+1)/s^2=3s^4-5s^2+1+s^{-2}$ である。原始関数 $3s^5/5-5s^3/3+s-s^{-1}$ の両端の差を取り、
$$
\Phi(y)=-\frac{\zeta}{15y^3}
  [\,16\sqrt{1+y}+9y^3+2y^2-8y-16\,]
$$
を得る。小さい $y$ では角括弧が $10y^3+\order{y^4}$、大きい $y$ では $9y^3+\order{y^2}$ である。従って $\Phi$ は $-2\zeta/3$ から $-3\zeta/5$ へ連続的に移る。この解で $k/\mathcal H\ll1$ を課しているため、等密度付近までに Hubble 半径内へ入る波数には適用できない。

:::example[等密度におけるポテンシャル]
$y=1$ では $\Phi/\zeta=-(16\sqrt2-13)/15\simeq-0.642$ となる。等密度で圧力比を急に変えた定数解に切り替えると、この中間値を再現できない。初期値と終値の比だけが必要なら両極限を使い、時間幅が物理に効く場合は連続解を使う。
:::

### エントロピーを曲率に変える局所的な崩壊

長波長で圧力なしの成分 $m$ が放射 $r$ に一斉に崩壊する模型を取る。崩壊は全密度が共通の面で起き、全密度と空間計量はその面で連続すると仮定する。崩壊前は個別の $\zeta_m,\zeta_r$ が保存している。崩壊面では $\delta\rho_m=3\bar\rho_m(\zeta_m-\zeta)$、$\delta\rho_r=4\bar\rho_r(\zeta_r-\zeta)$ で、全密度摂動が零なので
$$
\zeta_{\rm out}=\zeta_r+
\frac{3\bar\rho_m}{3\bar\rho_m+4\bar\rho_r}(\zeta_m-\zeta_r)
$$
となる。崩壊後は単一の放射なので、勾配を無視できる間はこの曲率が保存する。例えば崩壊前のエネルギー比が $\bar\rho_m:\bar\rho_r=4:3$ なら係数は $1/2$ である。実際に崩壊率が空間で変化して崩壊面が一様密度面からずれる場合や、崩壊が長く続く場合には、成分交換の時間発展を解く必要がある。

### 一般の平坦背景における断熱ポテンシャル

平坦背景に総エンタルピー $h\ne0$ があり、長波長で異方性応力を無視でき、$\zeta$ が一定の解を取る。宇宙時刻で $H=\dot a/a$ とする。時間制約で勾配を落とすと $\delta\rho=-3H(\dot\Phi+H\Phi)/(4\pi G)$ である。$\dot H=-4\pi G h$ と $\zeta=-\Phi+\delta\rho/(3h)$ を合わせれば
$$
\dot\Phi+\left(H-\frac{\dot H}{H}\right)\Phi
=\frac{\dot H}{H}\zeta .
$$
$\dd(a/H)/\dd t=(a/H)(H-\dot H/H)$ なので積分因子は $a/H$ である。$\zeta$ が一定なら $\int a \dot H/H^2\,\dd t=-a/H+\int a\,\dd t$ を使って
$$
\Phi(t)=-\zeta\left[1-\frac{H(t)}{a(t)}\int^{t}a(t')\dd t'\right]
+C_d\frac{H(t)}{a(t)}
$$
と書ける。積分の下端の変更は同次解の係数 $C_d$ に吸収される。$a\propto t^p$、$p>0$ の場合は下端を $0$ に選ぶと $\Phi=-\zeta/(p+1)$ である。例えば $p=1/2$ では $-2\zeta/3$、$p=2/3$ では $-3\zeta/5$ へ戻る。小さくても有限の $k$ に適用するには、その期間に $k/\mathcal H$ が小さいことも必要である。

:::example[純真空に近い背景の曲率]
$a\propto t^p$ で大きい $p$ を選べば、一定の $\zeta$ に対する $\Phi$ は $1/(p+1)$ に抑えられる。$H$ が定数へ近づく背景でも、残りの物質が曲率を定義する間は同様に小さくなり得る。しかし純真空の $\dot H=0$ では、導出で割った $h$ が零なので、一様密度面の曲率をこの方法で定義できない。真空に近い背景の極限と、純真空の摂動変数の定義は分けて扱う。
:::

### 長波長の断熱解と局所的な座標変更

厳密な零波数で、背景に時間移動 $\tilde\eta=\eta+T(\eta)$ と空間の小さい一様拡大 $\tilde x^i=(1+\lambda)x^i$ を施す。この変換は空間勾配を持たないため、Newton ゲージのシフト零を保つ。計量の時間と空間の一次変化は $\Psi=-\mathcal HT-T'$、$\Phi=\mathcal HT+\lambda$ である。全成分の密度は $\delta\rho_i=-\bar\rho_i'T$ なので、同じ時刻だけずれた断熱条件を満たす。背景にない波を作ったのではなく、零波数では座標の再定義である。この表示の曲率振幅は $-\Phi-\mathcal H\delta\rho/\bar\rho'=-\lambda$ となる。

小さい非零波数へ延長した正則解で $\Phi=\Psi$ を課せる場合を考える。零波数の式にもその等式を続けると $T'+2\mathcal HT=-\lambda$ である。従って $T=-\lambda a^{-2}\int a^2\dd\eta$ に $a^{-2}$ の同次項を足した解となる。$\int a^2\dd\eta=\int a\dd t$ と $\mathcal H/a^2=H/a$ なので、前のポテンシャル積分を再現する。有限波数で異方性応力が $k^2$ 程度あれば、跡なし制約で $k^2$ を割るため $\Phi-\Psi$ が零次に残ることがある。その場合はこの特定の延長条件を外して応力も解く。零波数の背景の座標変更だけから、粒子の応力込みの初期値を決めることはできない。

### 二成分の確率的な揺らぎを移す

長波長の交換なし系について $\zeta=(1-r)\zeta_1+r\zeta_2$、$r=h_2/(h_1+h_2)$ とする。初期に各曲率のパワーが $P_1,P_2$、交差パワーが $P_{12}$ なら、各成分が個々に保存する間の総パワーは
$$
P_\zeta=(1-r)^2P_1+r^2 P_2+2r(1-r)\operatorname{Re}P_{12}.
$$
これは Fourier 係数の二乗を平均する定義を展開したものなので、$P_1,P_2$ がそれぞれ非負でも、交差の符号により全パワーが小さくなる場合がある。共通の完全な断熱モードでは $\zeta_1=\zeta_2$ なので、重みが変わっても同じパワーになる。一方、無相関かつ等しいパワーなら $r=1/2$ の時点で全パワーは個別の $1/2$ になる。曲率の加重平均を平方せずにパワーを加重平均すると、この差を再現できない。転換を観測に結び付けるときは、成分間の平均振幅だけでなく相関も指定する。

## 演習問題

::::exercise{#exr-cosmo-perturbations-1 level=basic}
非零波数のテンソル摂動の独立数を求めよ。
:::solution
対称行列の六成分から横条件三つと跡条件一つを除いて二つである。
:::
::::

::::exercise{#exr-cosmo-perturbations-2 level=advanced}
$\Phi=C-\mathcal H\sigma$ の不変性を示せ。
:::solution
$\tilde C=C+\mathcal HT$ と $\tilde\sigma=\sigma+T$ を代入すると $T$ の二項が相殺する。
:::
::::

::::exercise{#exr-cosmo-perturbations-3 level=basic}
$\bar\rho'=-4\mathcal H\bar\rho$ のとき $\zeta$ を Newton ゲージで表せ。
:::solution
:eqref[eq-cosmo-bardeen] に入れて $-\Phi+\delta/4$ を得る。
:::
::::

::::exercise{#exr-cosmo-perturbations-4 level=advanced}
二つの Einstein 制約から速度を含む Poisson 形を導け。
:::solution
第二式で $\Phi'+\mathcal H\Psi$ を消去すると $k^2\Phi=-4\pi Ga^2[\delta\rho+3\mathcal H(\bar\rho+\bar p)\theta/k^2]$ となる。
:::
::::

::::exercise{#exr-cosmo-perturbations-5 level=basic}
重力を無視した放射の振動角周波数を求めよ。
:::solution
二つの式を消去すると係数は $k^2/3$ になるので $\omega=k/\sqrt3$ である。
:::
::::

::::exercise{#exr-cosmo-perturbations-6 level=advanced}
$\delta p_{\rm nad}$ が不変であることを示せ。
:::solution
二つの摂動の変換項は $-\bar p'T$ と $-\bar\rho'T$ であり、その差は $-\bar p'T+c_a^2\bar\rho'T=0$ である。
:::
::::

::::exercise{#exr-cosmo-perturbations-7 level=basic}
$\theta=0$、$\delta p_{\rm nad}/(\bar\rho+\bar p)=1/100$ で $a$ が $e^2$ 倍になる間の $\zeta$ の変化を求めよ。
:::solution
:eqref[eq-cosmo-zeta-conservation] を $\log a$ で積分して $-2/100=-0.02$ となる。
:::
::::

::::exercise{#exr-cosmo-perturbations-8 level=advanced}
自身の重力を無視した試験流体を、外部から与えた静的な Newton ゲージの計量に置く。圧力と異方性応力は零である。速度が零で、$\Psi$ の Fourier 振幅が正の初期時刻における、密度コントラスト $\delta$ の二階微分を求めよ。
:::solution
一定の $\Phi$ と応力零なら $\delta'=-\theta$、初期には $\theta'=k^2\Psi$ なので $\delta''=-k^2\Psi<0$ である。
:::
::::

::::exercise{#exr-cosmo-perturbations-9 level=advanced}
放射優勢下に自己重力を無視できる圧力なしの試験流体を置く。外部ポテンシャルを一定の $\Psi_0$ とし、$\theta(\eta_i)=0$ とする。速度発散と、その $\eta_i\to0$ の極限を求めよ。
:::solution
$\mathcal H=1/\eta$ なので速度式は $(\eta\theta)'=k^2\eta\Psi_0$ となる。$\eta_i$ から積分し $\eta$ で割ると $\theta=k^2\Psi_0(\eta-\eta_i^2/\eta)/2$ である。固定した $\eta>0$ で $\eta_i\to0$ とすると正則な $\theta=k^2\eta\Psi_0/2$ を得る。
:::
::::

::::exercise{#exr-cosmo-perturbations-10 level=advanced}
平坦な物質・放射宇宙の応力を零とする。長波長で両成分の曲率が共通の $\zeta$ であり、放射期の発散解はない。$y=\bar\rho_m/\bar\rho_r$ で $\Phi(y)$ が満たす一次方程式を導き、$\Phi(y\to\infty)$ を得よ。
:::solution
$k/\mathcal H$ の二次を無視した時間制約から $\delta=-2(\Phi+\Phi'/\mathcal H)=-2(\Phi+y\Phi_y)$ である。$1+w=(3y+4)/[3(1+y)]$ と $\zeta=-\Phi+\delta/[3(1+w)]$ を入れると $y\Phi_y+(5y+6)\Phi/[2(1+y)]=-(3y+4)\zeta/[2(1+y)]$ となる。積分因子 $y^3/\sqrt{1+y}$ を掛け、正則性から積分を $0$ から取る。被積分関数は $u^2(3u+4)/[2(1+u)^{3/2}]$ であり、$u\gg1$ では $(3/2)u^{3/2}$ である。積分の先頭は $(3/5)y^{5/2}$、外の因子は $y^{-5/2}$ なので、$\Phi\to-3\zeta/5$ となる。固定の非零波数で無限の未来まで長波長とするとは限らないため、この極限は勾配を捨てた解に対して取っている。
:::
::::

::::exercise{#exr-cosmo-perturbations-11 level=advanced}
ある有限の時間区間で、制約残差の初期値は $F_0=F_1=0$ だが $F_2=F$ を一定とする。膨張が無視できる短時間の近似で、終了時の $F_0,F_1$ を求めよ。特に跡なし制約の誤差が時間制約に入る次数を示せ。
:::solution
$\mathcal H$ の項を落とすと $F_1'=-k^2F/3$ なので、時間差 $u$ で $F_1=-k^2Fu/3$ である。$F_0'=F_1$ に入れて $F_0=-k^2Fu^2/6$ を得る。初期に時間制約と運動量制約が成り立っていても、跡なし制約が破れていれば、前者に時間の二次、後者に一次で誤差が入る。背景が変わる場合も、初期の Taylor 展開のこれらの先頭次数は同じである。
:::
::::

::::exercise{#exr-cosmo-perturbations-12 level=advanced}
定数 $w=0$ で $\delta p=0$、$\sigma_f=0$ とする。一方、ポテンシャルには成長モードと減衰モードが両方存在する。Newton ゲージで $\Phi=A+B\eta^{-5}$ から $\theta$ および $\mathcal R$ を求め、曲率に減衰成分が残るか確かめよ。
:::solution
$\mathcal H=2/\eta$、$4\pi Ga^2\bar\rho=6/\eta^2$、$\Psi=\Phi$ である。運動量制約から $\theta=k^2\eta^2[\Phi'+(2/\eta)\Phi]/6=k^2[\eta A/3-B\eta^{-4}/2]$ である。従って $\mathcal R=-\Phi-2\theta/(k^2\eta)=-5A/3$ となる。圧力がないので $\mathcal R'=0$ であり、減衰モードはこの変数の係数を持たない。ただし $\Phi$ や密度の減衰モードは残るので、$\mathcal R$ 一つだけで任意の有限波数の初期条件を指定したことにはならない。
:::
::::

::::exercise{#exr-cosmo-perturbations-13 level=advanced}
二成分の一様密度面の曲率が $\zeta_1=g$、$\zeta_2=-g$ であり、長波長で個別に保存するとする。$r=h_2/(h_1+h_2)$ を用いて全曲率とパワーを求め、$r=1/2$ の状況を説明せよ。
:::solution
$\zeta=(1-r)g-rg=(1-2r)g$ である。従ってパワーは $(1-2r)^2P_g$ となる。$r=1/2$ の面で全曲率は零だが、$S_{12}=3(\zeta_1-\zeta_2)=6g$ は非零である。二成分の圧力応答が違い、時間とともにエンタルピー比が変わるなら、その後に全曲率を再び生じる。この例では同じ場 $g$ から両成分が作られるため、二つのパワーを独立として足してはいけない。
:::
::::

::::exercise{#exr-cosmo-perturbations-14 level=advanced}
一成分の静止系音速が $c_s^2=1$、断熱音速が $c_a^2=0$ で、異方性応力はないとする。Newton ゲージのある時刻で $\delta\rho=0$、$\theta\ne0$ である。$h=\bar\rho+\bar p>0$ の下で $\delta p_{\rm nad}$ と $\mathcal R'$ を求めよ。
:::solution
共動面への移動は $T=\theta/k^2$ であり、$\Delta\rho=\delta\rho-\bar\rho'T=3\mathcal Hh\theta/k^2$ になる。従って $\delta p_{\rm nad}=(c_s^2-c_a^2)\Delta\rho=3\mathcal Hh\theta/k^2$ である。一方、共動圧力は $c_s^2\Delta\rho=\Delta\rho$ であるから、共動曲率の式で $\mathcal R'=-3\mathcal H^2\theta/k^2$ である。Newton ゲージの密度が零の瞬間でも、共動面の圧力は零と限らない。Einstein 制約を保つには、この時刻のポテンシャルと時間微分を指定した速度から合わせる。
:::
::::

::::exercise{#exr-cosmo-perturbations-15 level=basic}
Newton ゲージで背景の $\bar\rho'=-3\mathcal Hh$、$h>0$、$\delta\rho=6h A$、$\Phi=A$ とする。平坦面の密度と一様密度面の曲率を求めよ。
:::solution
Newton ゲージでは $C=\Phi=A$ なので、平坦面への時刻変換は $T=-A/\mathcal H$ である。従って $\delta\rho_{\rm flat}=\delta\rho+\bar\rho'A/\mathcal H=3h A$ となる。これを $\zeta=-\mathcal H\delta\rho_{\rm flat}/\bar\rho'$ に入れると $\zeta=A$ である。Newton ゲージの密度を直接 $\zeta=-\Phi+\delta\rho/(3h)$ に入れても $-A+2A=A$ となる。
:::
::::

## 確認テスト

::::quiz
$\Phi$ と $\Psi$ の差を作るものを選べ。

- [x] 異方性応力
- [ ] 一様な背景圧力だけ
- [ ] 宇宙時刻の単位

:::explanation
空間 Einstein 方程式の跡なし成分に入り、非零波数で差を与えるのは異方性応力である。
:::
::::

::::quiz{answer="0.333333" tol="1e-5"}
理想放射の音速の二乗を求めよ。
:::explanation
$\bar p=\bar\rho/3$ なので断熱な圧力微分は $1/3$ である。
:::
::::

::::quiz{answer="-0.01" tol="1e-6"}
$\theta=0$、$\delta p_{\rm nad}/(\bar\rho+\bar p)=0.01$、$\Delta\log a=1$ の $\Delta\zeta$ を求めよ。
:::explanation
:eqref[eq-cosmo-zeta-conservation] を積分して $-0.01$ になる。
:::
::::

::::quiz
一様密度ゲージで $\zeta$ に等しい量は何か。

- [x] $-C$
- [ ] $C$
- [ ] $\Psi$

:::explanation
$\delta\rho=0$ を :eqref[eq-cosmo-bardeen] に入れると $-C$ になる。
:::
::::

---
title: Dirac 方程式
description: 相対論的な一次波動方程式を構成し、スピノル、磁気モーメント、低速展開を導く。
notation:
  - symbol: '\gamma^\mu'
    meaning: 計量 (+,−,−,−) に対応する Dirac 行列
  - symbol: '\bar\psi'
    meaning: Dirac 随伴スピノル
  - symbol: '\boldsymbol\alpha,\beta'
    meaning: Hamiltonian 形式の Dirac 行列
  - symbol: '\boldsymbol\pi'
    meaning: 運動学的運動量 p−qA
---

## この章の目標

- 自由 Hamiltonian の線形化から Clifford 代数を得る。
- Lorentz 変換とスピノルの双一次形式を計算する。
- 自由平面波を構成し、スピンとヘリシティを区別する。
- 低速極限と Foldy–Wouthuysen 展開から磁場と電場の補正を求める。

## Dirac 方程式の導出

時間について一次の自己共役な Hamiltonian なら、初期波動関数一個から時間発展を決められる。空間の微分も一次にするため $H=c\sum_i\alpha_i p_i+\beta mc^2$ と置く。自由の $p_i$ は可換なので、$H^2=c^2\vb p^2+m^2c^4$ が $m>0$ を固定して任意の $\vb p$ で成り立つ条件は、係数ごとに

$$
\acomm{\alpha_i}{\alpha_j}=2\delta_{ij}I,\qquad
\acomm{\alpha_i}{\beta}=0,\qquad \beta^2=I
\label{eq-aqm-dirac-alpha}
$$

となる。$\alpha_i,\beta$ を Hermite とする。どれも他の一つとの共役で符号が反転するからトレースは零で、内部次元は偶数である。二次元の traceless Hermite 行列は三次元実ベクトルと Pauli 行列の積で表されるため、四つの互いに反交換する Hermite 行列で、二乗がすべて $I$ となるものには四本の直交ベクトルが必要となり不可能である。四次元なら

$$
\beta=\begin{pmatrix}I_2&0\\0&-I_2\end{pmatrix},\qquad
\alpha_i=\begin{pmatrix}0&\sigma_i\\\sigma_i&0\end{pmatrix}
$$

が条件を満たす。従って最小の複素表示は四成分である。

$(i\hbar\partial_t-H)\psi=0$ に $\beta/c$ を掛け、$\gamma^0=\beta$、$\gamma^i=\beta\alpha_i$、$x^0=ct$ と書く。すると :term[Dirac 方程式]{en="Dirac equation" yomi="でぃらっくほうていしき"} は

$$
(i\hbar\gamma^\mu\partial_\mu-mc)\psi=0
\label{eq-aqm-dirac}
$$

となる。$\bar\psi=\psi^\dagger\gamma^0$ を :term[Dirac 随伴]{en="Dirac adjoint" yomi="でぃらっくずいはん"} という。Hermite 共役に $\gamma^0$ を掛けて $i\hbar(\partial_\mu\bar\psi)\gamma^\mu+mc\bar\psi=0$ を得る。二式から $\partial_\mu(c\bar\psi\gamma^\mu\psi)=0$ である。時間成分は $j^0=c\psi^\dagger\psi$ なので、密度は両エネルギー枝で非負である。

:::example[静止 Hamiltonian]
$\vb p=0$ なら $H=\beta mc^2$ であり、上二成分は $+mc^2$、下二成分は $-mc^2$ を持つ。どちらの枝でも $\psi^\dagger\psi$ はノルム二乗である。
:::

:::example[一次元への制限]
$p_x=p_y=0$ なら $\sigma_z$ の各固有値 $s=\pm1$ が独立である。各二成分の行列は $\begin{pmatrix}mc^2&scp_z\\scp_z&-mc^2\end{pmatrix}$ で、特性式は $E^2-m^2c^4-c^2p_z^2$ になる。
:::

### 成分数が四つ必要になる理由

$\beta^2=I$ で $\alpha_i$ が $\beta$ と反交換するため、$\alpha_i$ は $\beta$ の正負部分空間を逆に移す。$\alpha_i^2=I$ で逆写像もあるから両次元は等しい。もし全次元が 2 なら $\beta=\sigma_z$ と対角化でき、これと反交換する Hermite 行列は $\alpha=a\sigma_x+b\sigma_y$ に限る。$\alpha^2=I$ は $(a,b)$ が単位ベクトルである条件であり、異なる $\alpha_i,\alpha_j$ の反交換は対応する二次元ベクトルの内積が零である条件である。二次元には互いに垂直な単位ベクトルが三つない。従って全次元は少なくとも 4 であり、本節の明示形で 4 を達成する。

:::example[二つの空間方向での二成分表示]
空間が二次元だけなら $\alpha_x=\sigma_x,\alpha_y=\sigma_y,\beta=\sigma_z$ で $H^2=(c^2p_x^2+c^2p_y^2+m^2c^4)I$ である。$p_x=mc,p_y=0$ では $H/(mc^2)=\sigma_x+\sigma_z$、固有値は $\pm\sqrt2$ で各一重である。$m>0$ の三次元問題の各枝の二重度を、任意の運動に対してこの二次元の行列で表すことはできない。
:::

### 確率流の時間成分を一階方程式から得る

$H=-i\hbar c\sum_i\alpha_i\partial_i+\beta mc^2+V$ とし、$V$ は局所な Hermite 行列とする。時間方程式から $\partial_t\psi=-c\alpha_i\partial_i\psi-i(\beta mc^2+V)\psi/\hbar$ である。随伴の式を加えると、$V$ と質量の項は互いに符号が逆で相殺する。従って
$$
\partial_t(\psi^\dagger\psi)=-c\sum_i[(\partial_i\psi)^\dagger\alpha_i\psi+\psi^\dagger\alpha_i\partial_i\psi]
=-\grad\cdot(c\psi^\dagger\boldsymbol\alpha\psi).
$$
二階時間の Klein–Gordon と違い、初期条件は $\psi$ 一つであり、時間成分は同じ時刻の $\psi^\dagger\psi$ から定まる。境界を流束が通らないなら合計は一定である。さらに正エネルギーに射影するかどうかはスペクトルの条件であり、この局所保存式を得るための仮定ではない。

:::example[周期箱で相殺する境界流束]
一辺 $L$ の箱の対向面に同じ値を課すと、$j_i$ は向かい側で同じで外向き法線は逆である。例えば $p_z=3mc/4$ の正枝で箱平均が $1/L^3$ の波なら、$\partial_t\rho=0$、$z$ 面の単位面積当たりの流束は $(3c/5)L^{-3}$ である。反対の面を合わせると総流出は零になる。局所の電流が零でなくても全確率を保てる。
:::

### 初期条件を保つ一階の Green 関数

自由運動量ごとに $H(\vb p)=c\boldsymbol\alpha\cdot\vb p+\beta M$、$M=mc^2$ とする。$H^2=E_{\vb p}^2 I$ のため、その指数の偶数項と奇数項は
$$
U(\vb p,t)=I\cos(E_{\vb p}t/\hbar)
-i\,\frac{H(\vb p)}{E_{\vb p}}\sin(E_{\vb p}t/\hbar)
$$
となる。初期 $\widetilde\psi(\vb p,0)$ にこれを掛けて空間 Fourier を戻せば、滑らかな自由の初期データからの解が得られる。各 $\vb p$ で $H$ は Hermite なので $U^\dagger U=I$、よって全確率も Fourier の等長性で保存される。遅延の逆核は $G_R(t,\vb p)=-(i/\hbar)\mathbf1_{t>0}U(\vb p,t)$ であり、微分の境界項が $(i\hbar)\cdot(-i/\hbar)\delta(t)I$ になる。Klein–Gordon の核では導関数の跳びが二階方程式の源を作るが、Dirac の一階方程式では核そのものの跳びが源になる。

空間でコンパクトな台を持つ初期データに対し、因果的な伝播はさらに局所保存式から証明できる。$|\vb j|\le c\rho$ なので、速度 $c$ で縮む球の中の $\int\rho$ の時間微分は、流束による $-\int\vb j\cdot\hat{\vb n}$ と境界の移動による $-c\int\rho$ の和で零以下である。初期に球内で $\psi=0$ なら、その未来の縮む球でも $\rho=0$ である。静的な正負射影を空間で先に行う場合の非局所性とは別の主張である。

:::example[自由な上成分の時間発展を一つの波数で見る]
$p_z=3mc/4$ で $z$ 上向きスピンを選ぶと、上下の二成分に $H/M=\begin{pmatrix}1&3/4\\3/4&-1\end{pmatrix}$ が作用する。$E=5M/4$、初期 $(1,0)^\mathsf T$ なので、下の振幅は $-i(3/5)\sin(Et/\hbar)$ である。下を測る最大確率は $9/25$ である。初期の上成分が正エネルギーだけであれば、定常な正スピノルの上・下比は変わらない。ここでの変化は、初期の上成分だけを指定するとその運動量では両エネルギー枝を含むことによる。
:::

## γ 行列と Clifford 代数

:::definition[Clifford 代数]{#def-aqm-clifford}
計量 $\eta=\diag(1,-1,-1,-1)$ に対する :term[Clifford 代数]{en="Clifford algebra" yomi="くりふぉーどだいすう"} の Dirac 表現では

$$
\acomm{\gamma^\mu}{\gamma^\nu}=2\eta^{\mu\nu}I,\qquad
(\gamma^\mu)^\dagger=\gamma^0\gamma^\mu\gamma^0 .
\label{eq-aqm-gamma}
$$

ベクトルの収縮を $\not{a}=\gamma^\mu a_\mu$ と書く。
:::

対称な積 $a_\mu a_\nu$ に反対称成分が寄与しないため、$\not a^2=a^2I$ である。:eqref[eq-aqm-dirac] に $(i\hbar\gamma^\mu\partial_\mu+mc)$ を作用させると、各成分は :eqref[eq-aqm-kg] を満たす。逆向きは一般に成り立たず、Dirac 方程式は成分の間に追加の関係を課す。

$\gamma^5=i\gamma^0\gamma^1\gamma^2\gamma^3$ と置く。異なる行列同士を三回交換すれば $\{\gamma^5,\gamma^\mu\}=0$ となり、二乗の並べ替えに六回の交換を使うと $(\gamma^5)^2=I$ である。随伴を取ると三つの空間成分の負号と $i$ の共役が相殺し、四因子の反転は偶置換なので $\gamma^5$ は Hermite である。よって $P_\pm=(I\pm\gamma^5)/2$ は互いに直交する射影であり、これをカイラリティ射影という。

:::example[二行列のトレース]
$\tr(\gamma^\mu\gamma^\nu)=4\eta^{\mu\nu}$ である。巡回性により左右の順序のトレースが等しく、反交換子をトレースして得る。異なる二成分では零である。
:::

:::example[スラッシュ行列の逆]
$a^2\ne m^2c^2$ なら $(\not a-mc)^{-1}=(\not a+mc)/(a^2-m^2c^2)$ である。両順序の積を計算してスカラーになる。例えば $a^\mu=(2mc,0,0,0)$ では分母は $3m^2c^2$ となる。
:::

### 四行列のトレース

$T=\tr(\gamma^\mu\gamma^\nu\gamma^\rho\gamma^\sigma)$ とする。最初の $\gamma^\mu$ を右へ順に移すたびに反交換子が一項ずつ出る。最後に残る積は負号が三回掛かり、巡回性で $-T$ になる。従って
$$
T=4(\eta^{\mu\nu}\eta^{\rho\sigma}
-\eta^{\mu\rho}\eta^{\nu\sigma}
+\eta^{\mu\sigma}\eta^{\nu\rho}).
$$
例えば最初に現れる項は $2\eta^{\mu\nu}\tr(\gamma^\rho\gamma^\sigma)=8\eta^{\mu\nu}\eta^{\rho\sigma}$ であり、最後に $2T$ を割る。奇数個の $\gamma^\mu$ のトレースは $\gamma^5$ で共役すると符号が替わるため零である。

:::example[内積だけでスラッシュ行列の積を調べる]
$a^\mu=(2,1,0,0)$、$b^\mu=(3,0,1,0)$ を同じ運動量単位で取る。$a^2=3,b^2=8,a\cdot b=6$ で、$\tr(\not a\not b)=24$ である。順序を交互にした四因子のトレースは $4[2(6)^2-3\cdot8]=192$、$\not a^2\not b^2$ では $4\cdot3\cdot8=96$ である。トレースの巡回性では隣り合わない因子の交換はできない。
:::

### 四成分の別表示を等長変換で結ぶ

Hermite な $\beta$ と $\alpha_i$ で代数を満たす四次元の表示を取る。$\beta$ の正負部分空間は各二次元で、$\alpha_3$ により相互に等長に対応する。$\beta=\diag(I,-I)$ とし、上の基底から下の基底を $\alpha_3$ で移して、$\alpha_3=\begin{pmatrix}0&I\\I&0\end{pmatrix}$ とできる。$i=1,2$ の $\alpha_i$ はオフ対角で、上右ブロックを $A_i$ とすれば $\alpha_i=\begin{pmatrix}0&A_i\\A_i^\dagger&0\end{pmatrix}$ である。$\alpha_3$ との反交換から $A_i^\dagger=-A_i$、二乗が $I$ から $A_i^\dagger A_i=I$ である。従って $A_i=iB_i$ で $B_i$ は二次元の Hermite な反転であり、$B_1B_2=-B_2B_1$ である。$B_1$ の固有値は $\pm1$ が一つずつで、$B_2$ は両固有ベクトルを逆に移す。上・下の基底を同じユニタリで調整すると $B_1=\sigma_z,B_2=\sigma_x$ に固定できる。この規約への構成は任意の四成分表示に使えるため、二つの表示はユニタリ同値である。

:::example[Dirac 表現から Weyl 表現へ移す]
$W=2^{-1/2}\begin{pmatrix}I&-I\\I&I\end{pmatrix}$ を Dirac 表現の列に作用させる。$W^\dagger W=I$ であり、$\beta$ はオフ対角の $\begin{pmatrix}0&I\\I&0\end{pmatrix}$、$\alpha_i$ は $\diag(-\sigma_i,\sigma_i)$ に移る。さらに $\gamma^5=\begin{pmatrix}0&I\\I&0\end{pmatrix}$ は $\diag(-I,I)$ になる。例えば静止の $(\chi,0)$ は $(\chi,\chi)/\sqrt2$ へ移るので、上下のノルムが大・小を表すのは Dirac 表現の性質である。エネルギー固有値と四元電流は、演算子も同時に共役して計算すれば同じになる。
:::

### 行列基底の展開をトレースで行う

独立な十六個を $I,\gamma^\mu,\sigma^{\mu\nu}\ (\mu<\nu),\gamma^\mu\gamma^5,\gamma^5$ とする。それぞれの積のトレースには異なる部分集合を持つ二因子の間で零になる性質がある。実際、積は対称差に当たる $\gamma$ の積に比例し、偶数個なら積に含む一因子、奇数個なら積に含まない一因子で共役して符号を反転できる。同じ部分集合の積は非零のスカラーなので、これら十六個は線形独立である。ここで $\sigma^{\mu\nu}=(i/2)[\gamma^\mu,\gamma^\nu]$ とした。このため二つの添字を上げ下げして双対を作れば、任意の $4\times4$ 行列 $A$ を
$$
\frac14\left[
I\tr A+\gamma^\mu\tr(\gamma_\mu A)
+\frac12\sigma^{\mu\nu}\tr(\sigma_{\mu\nu} A)
-\gamma^\mu\gamma^5\tr(\gamma_\mu\gamma^5 A)
+\gamma^5\tr(\gamma^5 A)
\right]
$$
に展開できる。各軸性の積では $\gamma^5$ を一回通すと $-$ になるので、双対に負号が必要である。反対称テンソルの $\mu,\nu$ を制限なしで足すと二回数えるため $1/2$ が掛かる。例えば $\tr(\sigma^{\rho\sigma}\sigma_{\mu\nu})=4(\delta^\rho_\mu\delta^\sigma_\nu-\delta^\rho_\nu\delta^\sigma_\mu)$ であり、これと $1/2$ がちょうど反対称成分を取り出す。式を各基底行列と掛けて跡を取ると $A$ と一致し、基底の完全性から行列として等しい。

:::example[静止の一つのスピン射影のトレース展開]
$u=(1,0,0,0)^\mathsf T$ とする。$u\bar u$ のスカラー係数の跡は 1、ベクトルの時間成分も 1 で、空間電流は零である。軸性の空間 $z$ 成分と $\sigma^{12}$ は 1 である。従って
$$
u\bar u=\frac14(I+\gamma^0+\gamma^3\gamma^5+\sigma^{12})
$$
となる。$\gamma_3=-\gamma^3$ と軸性の双対の負号が相殺する。右辺は $\frac14(I+\beta)(I+\Sigma_z)$ であり、上成分の $z$ スピン正だけを残す。二つの制約のどちらかを足すだけでは一成分の射影にはならない。
:::

## Lorentz 共変性とスピノルの変換

$x'=\Lambda x$ とする。$\psi'(x')=S(\Lambda)\psi(x)$ が同じ方程式を満たすには

$$
S^{-1}\gamma^\mu S=\Lambda^\mu{}_\nu\gamma^\nu
\label{eq-aqm-spin-covariance}
$$

が成り立てばよい。連鎖律による $\partial'_\mu=(\Lambda^{-1})^\nu{}_\mu\partial_\nu$ とこの関係を組み合わせると、変換後の方程式は元の式に左から $S$ を掛けたものになる。

$\sigma^{\mu\nu}=(i/2)[\gamma^\mu,\gamma^\nu]$ と定める。$\gamma^\rho$ を二因子の積の反対側へ動かし、各反交換子を使うと

$$
[\sigma^{\mu\nu},\gamma^\rho]
=2i(\eta^{\nu\rho}\gamma^\mu-\eta^{\mu\rho}\gamma^\nu)
$$

である。これが回転と boost の生成子を作る理由である。具体的に $\Sigma_i=\diag(\sigma_i,\sigma_i)$ とすれば $z$ 回転は $S_R=e^{-i\theta\Sigma_z/2}$ である。$x'^0=x^0\cosh\zeta+x^1\sinh\zeta$、$x'^1=x^0\sinh\zeta+x^1\cosh\zeta$ の boost には $S_B=e^{\zeta\alpha_x/2}$ を取る。例えば $\gamma^0\alpha_x=\gamma^1$、$\alpha_x\gamma^0=-\gamma^1$ なので $S_B^{-1}\gamma^0S_B=\gamma^0\cosh\zeta+\gamma^1\sinh\zeta$ になる。回転と boost の積で固有順時 Lorentz 変換を作れるため、これらの積で共変性を実現できる。

回転は $\beta$ と可換でユニタリである。boost では $\beta$ と $\alpha_x$ が反交換なので $S_B^\dagger\beta S_B=\beta$ となり、同じ関係は任意の積にも成り立つ。従って $\bar\psi'(x')=\bar\psi(x)S^{-1}$ である。boost は通常の四成分内積でユニタリとは限らない。保存される空間積分は平行時刻面だけに固定せず、カレントを空間的超曲面の法線に収縮した $\int_\Sigma j^\mu\dd\Sigma_\mu/c$ である。面の間の境界流束が消えるなら発散定理で二つの面の積分が一致する。

:::example[$2\pi$ の回転]
$\Sigma_z$ の固有値は $\pm1$ なので $S_R(2\pi)=-I$、$S_R(4\pi)=I$ である。一方 $\gamma^\mu$ の共役作用はどちらも恒等である。
:::

:::example[静止スピノルの boost]
$x$ 方向 boost を $(\chi,0)^\mathsf T$ に作用させると $(\cosh(\zeta/2)\chi,\sinh(\zeta/2)\sigma_x\chi)^\mathsf T$ となる。$\|\chi\|=1$ なら $\bar uu=\cosh^2(\zeta/2)-\sinh^2(\zeta/2)=1$、$u^\dagger u=\cosh\zeta$ である。
:::

### 異方向の boost を交換する

$\zeta_x,\zeta_y$ が小さいとき、$e^{\zeta_x\alpha_x/2}e^{\zeta_y\alpha_y/2}e^{-\zeta_x\alpha_x/2}e^{-\zeta_y\alpha_y/2}$ を展開する。一次項は相殺し、二次で異なる角を含む項は $\zeta_x\zeta_y[\alpha_x,\alpha_y]/4=i\zeta_x\zeta_y\Sigma_z/2$ である。同じ角を二回含む項は各逆元と相殺する。従って一周の作用はこの次数で $z$ 回転の角 $-\zeta_x\zeta_y$ に対応する。行列の共役作用を使えば座標変換にも同じ回転が現れる。単一方向の boost に限ると生成子が同じなので、この回転は現れない。

:::example[向きを保った boost の連続合成]
$x$ 方向だけなら $S(\zeta)S(\eta)=S(\zeta+\eta)$ である。$\zeta=\eta=\log2$ では各 $\cosh\zeta=5/4,\sinh\zeta=3/4$ で、単一 boost の速さは $3c/5$ である。組み合わせると $\cosh(2\zeta)=17/8,\sinh(2\zeta)=15/8$、速さは $15c/17$ になる。静止スピノルを移した $u^\dagger u$ は $17/8$ であり、Dirac 随伴の $\bar uu=1$ は保たれる。
:::

### 微小な変換の群を復元する

$S=e^X$ とし $X$ を二つの $\gamma$ の反対称積から作る。微小な変換で $e^{-X}\gamma^\mu e^X=\gamma^\mu+[\gamma^\mu,X]+\cdots$ である。共役で生じる一次の変化は $\gamma^\mu$ の線形結合であり、同じ $X$ で繰り返しても四次元の部分空間 $\operatorname{span}(\gamma^\mu)$ から出ない。従って制御変数を連続に増やした有限変換は、この部分空間で一次行列を指数にした Lorentz 変換になる。$\gamma$ の反交換子を共役するとスカラーの $2\eta^{\mu\nu}I$ を保つので、この共役作用は必ず計量を保つ。随伴で $\beta$ を保つ条件を別に確認すれば、双一次形式の変換が得られる。

同じ Lorentz 行列を実現する二つの $S$ の商は、すべての $\gamma^\mu$ と可換である。$\gamma$ の積が $4\times4$ 行列を張ることを使うと商はスカラーである。実際、すべての行列単位と可換な行列は対角で、対角成分がすべて等しい。このスカラーは $\beta$ を保つ条件で絶対値 1 になる。通常のスピン表現では $S$ を反対称積の指数で作り、Weyl 表現では左右の二成分に分かれ、それぞれの行列式が 1 になる。実際、回転と boost の生成子の各二次ブロックは Pauli 行列に比例し、トレースが零である。従って Lorentz 行列が同じ二つの積の商を $zI$ とすると、二成分ごとの行列式から $z^2=1$ となる。$-I$ も $2\pi$ 回転として実現するので、固有順時 Lorentz 群の各元に二つのスピノル変換が対応する。このとき連続な曲線を持ち上げる際は初期の $S$ を指定する必要があり、閉じた $2\pi$ 回転が $-I$ に終わる。密度などの双一次形式には共役作用だけが出るので、この符号自体はそこから消える。

:::example[boost の前後でスカラーと密度を照合する]
静止スピノル $(\chi,0)$、$\chi=(1,0)$ を $\zeta=\log 2$ の $z$ boost で移す。半角は $\cosh(\zeta/2)=3/(2\sqrt2)$、$\sinh(\zeta/2)=1/(2\sqrt2)$ である。従って $u^\dagger u=9/8+1/8=5/4$、$\bar uu=9/8-1/8=1$ となる。$u^\dagger\alpha_z u$ は交差の二倍で $3/4$ である。時間成分と空間成分が $\gamma(1,v/c)=(5/4,3/4)$ となり、正の局所密度がフレームで数値を変えることと、カレントの共変性が両立する。
:::

### 固有時向きの変換で保存積分を運ぶ

座標変換 $x'=\Lambda x$ とスピノル変換 $\psi'(x')=S\psi(x)$ を同時に行う。双一次ベクトルの共変性から $j'^\mu(x')=\Lambda^\mu{}_\nu j^\nu(x)$ であり、空間的超曲面 $\Sigma$ の向き付き面要素 $\dd\Sigma_\mu$ は共変に変換する。従って $\int_\Sigma j^\mu\dd\Sigma_\mu$ は同じ幾何学的な面を別の座標で表示しても不変である。別の面と比較する場合には保存則を使う。例えば $t'=0$ 面は $t=0$ 面と一般に一致しないので、二つの面の間の四体積に発散定理を使い、無限遠の境界の流束が消えることを仮定する。このため Lorentz の boost が四成分ベクトルの通常の等長変換でなくても、解空間の内積は保存される。$\rho=\psi^\dagger\psi\ge0$ と $|\vb j|\le c\rho$ により、未来向きの法線との縮約も非負である。

:::example[局所密度と平行移動した薄い体積の積]
静止近傍の電流が $(c\rho,0,0,0)$ と近似できる波束の一部分を考える。$v=3c/5$ の $z$ boost で時間成分は $5c\rho/4$ になる。一方、元の波束の $z$ の幅を boost 後の同時刻面で比較するには、静止波形の $z$ 依存に $z=\gamma z'$ を入れる。この面上で幅は $4/5$ 倍なので、定常な密度を仮定した計算では積が保たれる。時間変化のある波束の場合にこの単純な縮みで置き換えず、超曲面のカレントを使う。
:::

## 双一次形式

スピノルを一つずつ変換しても、その成分は直接には四元ベクトルでない。bra に Dirac 随伴を使うと、共役作用が $\gamma$ 行列だけにかかり、Lorentz 変換の型を分類できる。

固有順時変換で $\bar\psi\psi$ はスカラー、$\bar\psi\gamma^\mu\psi$ はベクトル、$\bar\psi\sigma^{\mu\nu}\psi$ は反対称二階テンソルである。第一は $S^{-1}S=I$、第二は :eqref[eq-aqm-spin-covariance]、第三は同式を二因子に適用して従う。$\gamma^5$ は $\gamma^\mu\gamma^\nu$ の偶積と可換なので、生成子の指数とも可換である。従って $i\bar\psi\gamma^5\psi$ は固有順時変換でスカラー、$\bar\psi\gamma^\mu\gamma^5\psi$ はベクトルとして変換する。空間反転を $S_P=\beta$ とすると $\beta\gamma^5\beta=-\gamma^5$ になり、これらは :term[擬スカラー]{en="pseudoscalar" yomi="ぎすからー"} と軸性ベクトルである。反転で通常のベクトルは時間成分が偶・空間成分が奇、軸性ベクトルは逆となる。

これら五種類の双一次形式の成分は合計 $1+4+6+4+1=16$ 個である。対応する行列は前節で示した基底をなし、スピノルの外積 $\psi\bar\psi$ をこれらの成分から再構成できる。

:::example[静止電流]
$\psi=(\chi,0)^\mathsf T$ なら $j^0=c\chi^\dagger\chi$、$j^i=c\psi^\dagger\alpha_i\psi=0$ である。静止状態の空間電流は消えるが、軸性ベクトルの空間成分 $\psi^\dagger\alpha_i\gamma^5\psi$ は $\chi^\dagger\sigma_i\chi$ であり、スピンを記述する。
:::

:::example[電流の大きさの上限]
任意の単位方向 $\hat{\vb n}$ に対して $\hat{\vb n}\cdot\boldsymbol\alpha$ の二乗は $I$ なので固有値は $\pm1$ である。従って $|\hat{\vb n}\cdot\vb j|\le c\psi^\dagger\psi$ であり、$\hat{\vb n}$ を $\vb j$ に沿わせると $|\vb j|\le j^0$ を得る。
:::

### 双一次形式の実性と軸性カレント

$\bar\psi M\psi$ の複素共役は $\bar\psi(\beta M^\dagger\beta)\psi$ である。従って $M=I,\gamma^\mu,\sigma^{\mu\nu},i\gamma^5,\gamma^\mu\gamma^5$ の場合は $\beta M^\dagger\beta=M$ を使って実数となる。$M=\gamma^5$ では $\beta\gamma^5\beta=-\gamma^5$ なので、$i$ が必要になる。最小結合の共変微分を使い、二つの Dirac 方程式から
$$
(\partial_\mu\bar\psi)\gamma^\mu=+\frac{imc}{\hbar}\bar\psi+\text{電磁項},
\qquad
\gamma^\mu(\partial_\mu\psi)=-\frac{imc}{\hbar}\psi+\text{電磁項}
$$
を得る。積の微分では電磁項が反対の電荷で相殺する。普通のベクトルでは質量項も相殺する。一方、軸性カレントに $\gamma^5$ を挟むと $\gamma^\mu$ との反交換で第二項の符号が逆になり、
$$
\partial_\mu(\bar\psi\gamma^\mu\gamma^5\psi)=\frac{2imc}{\hbar}\bar\psi\gamma^5\psi
$$
になる。これは古典スピノル、または一粒子方程式の恒等式である。

:::example[静止スピンの軸性カレントを boost する]
静止で $\langle\sigma_z\rangle=1$ とすると軸性ベクトル $A^\mu$ は $(0,0,0,1)$ である。$z$ 方向の boost の $\cosh\zeta=5/4,\sinh\zeta=3/4$ を使うと $(3/4,0,0,5/4)$ になる。同じ規格化のベクトル電流 $j^\mu/c$ は $(5/4,0,0,3/4)$ である。内積は $(15-15)/16=0$、軸性ベクトルの二乗は $-1$ である。スピンに結び付く軸性成分は運動する系で時間成分も持つ。
:::

### 反対称成分の双対をカイラリティで表す

$\epsilon^{0123}=+1$ の記号を取る。Dirac 表現の $\gamma^5=i\gamma^0\gamma^1\gamma^2\gamma^3$ から
$$
\sigma^{\mu\nu}\gamma^5=\frac i2\epsilon^{\mu\nu\rho\sigma}\sigma_{\rho\sigma}
$$
となる。例えば $(\mu,\nu)=(0,1)$ では $\sigma^{01}=i\gamma^0\gamma^1$、積は $-\gamma^0\gamma^1\gamma^0\gamma^1\gamma^2\gamma^3=-\gamma^2\gamma^3$ である。右辺は $i\sigma_{23}=-\gamma^2\gamma^3$ で一致する。残りは添字の置換で、そのつど計量の符号と $\gamma$ の交換符号を追う。例えば $(2,3)$ では $\gamma^2\gamma^3$ の二乗が $-I$ なので左辺は $+\gamma^0\gamma^1$、右辺は $i\sigma_{01}=+\gamma^0\gamma^1$ になる。従って $\gamma^5$ を付けた二階テンソルを双一次形式の独立な分類にもう六個加える必要はない。

:::example[磁気テンソルを静止状態で測る]
$\psi=(\ket{+z},0)^\mathsf T$ なら $\bar\psi\sigma^{12}\psi=1$ であり、$\sigma^{12}=\diag(\sigma_z,\sigma_z)$ を使った。$\bar\psi\sigma^{03}\psi=0$ は $\sigma^{03}=i\alpha_z$ が上下を結ぶからである。しかし $\bar\psi\sigma^{03}\gamma^5\psi=i$ であり、双対式の右辺 $i\bar\psi\sigma_{12}\psi$ と一致する。双一次形式の実性を保つ規約で双対を使う場合には、この $i$ を含む規約を扱う。
:::

### スピノルの質量項を作用の変分で確認する

Lagrange 密度を
$$
\mathcal L=\frac{i\hbar c}{2}
[\bar\psi\gamma^\mu\partial_\mu\psi-(\partial_\mu\bar\psi)\gamma^\mu\psi]
-mc^2\bar\psi\psi
$$
と取る。この対称化で、微分項は共役しても同じ実値を持つ。$\bar\psi$ を独立に変分し第二の項を部分積分すると、境界を除いて $i\hbar c\gamma^\mu\partial_\mu\psi-mc^2\psi$ が係数になる。これを零にすると :eqref[eq-aqm-dirac] が戻る。$\psi\to e^{i\alpha(x)}\psi$ を変分のため局所の関数へ広げると、一次の変化は $-\hbar c(\partial_\mu\alpha)\bar\psi\gamma^\mu\psi$ である。質量項は位相が打ち消す。方程式を満たすと作用の一次変分が零なので、コンパクト台の任意の $\alpha$ について部分積分して $\partial_\mu(\bar\psi\gamma^\mu\psi)=0$ を得る。この導出は時間並進からの保存量と、連続な位相からの保存量を区別する。

:::example[純粋な Weyl 成分と質量のスカラー]
Weyl 表現で $\psi=(\psi_L,0)$ とすれば $\bar\psi\psi=\psi_R^\dagger\psi_L+\psi_L^\dagger\psi_R=0$ である。しかし $\psi^\dagger\psi=\|\psi_L\|^2$ は正である。$\bar\psi\psi$ 単独にこの代入をして零になっても、先に $\psi_R=0$ に制約して変分すれば質量による左右の結合を落とす。このため純粋な一方の Weyl 成分が独立な運動方程式を満たすのは、質量を零にした場合である。
:::

## 平面波解とスピン・ヘリシティ

$E=\sqrt{c^2p^2+m^2c^4}>0$ とし、正枝では$\psi=u(\vb p)e^{-i(Et-\vb p\cdot\vb r)/\hbar}$ とする。$(c\boldsymbol\alpha\cdot\vb p+\beta mc^2)u=Eu$ の下二成分から、上二成分 $\chi$ に対し下は $c\boldsymbol\sigma\cdot\vb p\,\chi/(E+mc^2)$ になる。$\chi^\dagger\chi=1$、$(\boldsymbol\sigma\cdot\vb p)^2=p^2I$ として

$$
u(\vb p,\chi)=\sqrt{\frac{E+mc^2}{2E}}
\begin{pmatrix}\chi\\ c\boldsymbol\sigma\cdot\vb p\,\chi/(E+mc^2)\end{pmatrix}
\label{eq-aqm-dirac-u}
$$

とすれば $u^\dagger u=1$ である。同じ空間波数で $-E$ の固有ベクトルは同じ前因子に $(-c\boldsymbol\sigma\cdot\vb p\,\xi/(E+mc^2),\xi)^\mathsf T$ を掛ける。内積で上下の寄与が相殺するので正枝と負枝は直交する。各枝二状態で全四次元を張り、射影は $(I\pm H(\vb p)/E)/2$ である。これは二乗と積を $H^2=E^2I$ で確認できる。

$\hat{\vb p}\cdot\boldsymbol\Sigma/2$ の固有値を :term[ヘリシティ]{en="helicity" yomi="へりしてぃ"} という。$\boldsymbol\sigma\cdot\vb p$ と同じ軸の Pauli 行列であるから、固定運動量で $H$ と可換である。質量があると粒子を追い越す boost によって $\vb p$ の向きが逆となり、縦スピンを保ったままヘリシティが反転し得る。$m=0$ では $\gamma^5$ が $H$ と可換で、正枝において $\gamma^5=\boldsymbol\Sigma\cdot\hat{\vb p}$ が成り立つ。実際 Dirac 表現では $\gamma^5\boldsymbol\Sigma=\boldsymbol\alpha$ なので、正枝の $H/(cp)=1$ から従う。

:::example[速度の期待値]
速度演算子は $(i/\hbar)[H,\vb r]=c\boldsymbol\alpha$ である。:eqref[eq-aqm-dirac-u] に入れると、上下の交差項に $\sigma_i(\boldsymbol\sigma\cdot\vb p)+(\boldsymbol\sigma\cdot\vb p)\sigma_i=2p_i I$ が出て、$u^\dagger c\alpha_i u=c^2p_i/E$ になる。
:::

:::example[小成分の重み]
$E=2mc^2$ の正エネルギー状態では、下二成分のノルム二乗は $(E-mc^2)/(2E)=1/4$ である。静止に近づくとこの重みは零に近づく。
:::

### 運動量とスピン電流の分離

同じ質量の二つの正枝スピノル $u(p),u(p')$ を使う。$\not p\,u(p)=mc\,u(p)$ と $\bar u(p')\not p'=mc\,\bar u(p')$ を左右から入れると

$$
2mc\,\bar u(p')\gamma^\mu u(p)
=\bar u(p')(\not p'\gamma^\mu+\gamma^\mu\not p)u(p).
$$

ここに $\gamma^\mu\gamma^\nu=\eta^{\mu\nu}I-i\sigma^{\mu\nu}$ を使えば

$$
\bar u(p')\gamma^\mu u(p)=
\frac1{2mc}\bar u(p')\big[(p'+p)^\mu+i\sigma^{\mu\nu}(p'-p)_\nu\big]u(p).
$$

これを Gordon 分解という。最初の項は平均運動量に比例し、二番目は移行運動量を含むスピンの項である。後者にさらに $(p'-p)_\mu$ を掛けると反対称性で零である。前者にも同じ収縮をすると $p'^2-p^2=0$ なので、行列要素での電流保存が確認できる。$p'=p$ の場合は最初の項だけが残る。:eqref[eq-aqm-dirac-u] の規格化なら $\bar uu=mc^2/E$ なので、$j^\mu=c^2p^\mu/E$ となり、時間成分は $c$ である。

### スピン和と二種類の規格化

:eqref[eq-aqm-dirac-u] と同じ $u_s^\dagger u_{s'}=\delta_{ss'}$ の基底なら $\sum_s u_su_s^\dagger=P_+$ である。右から $\beta$ を掛けると
$$
\sum_s u_s(\vb p)\bar u_s(\vb p)
=\frac12(I+H/E)\beta=\frac{c}{2E}(\not p+mc).
$$
最後は $(I+H/E)\beta=\gamma^0+(mc^2-c\gamma^i p_i^{\rm sp})/E$ を使った。ここで $p_i^{\rm sp}$ は三次元の Cartesian 成分である。Lorentz 不変なスカラー規格化 $\bar uu=1$ を採るには $u\mapsto\sqrt{E/(mc^2)}\,u$ とする。するとスピン和は $(\not p+mc)/(2mc)$ となり、エネルギーに依存する前の規格化とは係数が異なる。負枝を $e^{+ip\cdot x/\hbar}$ に付けて $p^0>0$ とするなら $\sum v\bar v=c(\not p-mc)/(2E)$ である。空間位相の符号が正枝と逆なので $H(-\vb p)$ の負射影を使って求める。

:::example[規格化で見かけの密度が変わる場合]
$E=2mc^2$ では $u^\dagger u=1$ のスピノルの上・下の重みは $3/4,1/4$、差は $1/2$ である。スカラー規格化にすると両方を二倍し、重みは $3/2,1/2$、差は 1 になる。後者で成分を確率として解釈したいときは、全ノルム 2 でもう一度割る。両方の表記で運動方向の $j^z/j^0$ は $cp_z/E$ に揃う。
:::

### スピン偏極の射影を共変に組み立てる

静止のスピン軸 $\vb n$ を単位とし、$s^\mu=(0,\vb n)$ を運動量 $p^\mu$ と同じ boost で動かす。$s^2=-1$、$s\cdot p=0$ が保たれる。$G=\gamma^5\not s$ とすると $\{\gamma^5,\not s\}=0$、$\not s^2=-I$ から $G^2=I$ である。$\not p$ は $\gamma^5$ と $\not s$ の双方に反交換するので $G$ と可換である。静止では $G=\diag(\boldsymbol\sigma\cdot\vb n,-\boldsymbol\sigma\cdot\vb n)$ であり、正枝の上成分の偏極射影は $(I+G)/2$ になる。共変性で動く系でも同じ条件が成り立つ。正規化を $\bar uu=1$ とすれば、一つの偏極の外積は
$$
u(p,s)\bar u(p,s)=\frac{\not p+mc}{2mc}\frac{I+\gamma^5\not s}{2}.
$$
静止で両辺が同じ階数一の射影であることを見れば、同じ boost を共役することによって全運動量に成り立つ。局所の確率密度で固有ベクトルを規格化する場合は、前節の $mc^2/E$ を掛ける。

:::example[進行方向と垂直な偏極]
$\vb p=p\hat{\vb z}$、静止の $\vb n=\hat{\vb x}$ なら同方向の boost は $x$ 成分を変えず、$s=(0,1,0,0)$ のままである。これはヘリシティの一固有状態ではなく、正枝では $\hat z$ の二偏極の等振幅の和である。したがって $\pm1/2$ の確率を各 $1/2$ と得る。$\vb n=\hat{\vb z}$ にした場合には $s=(p/(mc),0,0,E/(mc^2))$ となり、$E=5mc^2/4,p=3mc/4$ なら $s=(3/4,0,0,5/4)$ である。この値は $s^2=-1$ と $s\cdot p=0$ を直接満たす。
:::

### 正負の極で射影が異なることを確認する

固定 $\vb p$ のエネルギー変数 $z$ に対する逆行列は
$$
(zI-H)^{-1}=\frac{zI+H}{z^2-E^2}
=\frac{P_+}{z-E}+\frac{P_-}{z+E}.
$$
最後は $I=P_++P_-$ と $H=E(P_+-P_-)$ で確認できる。これに右から $\beta$ を掛ければ共変なスラッシュ行列の分子と対応する。可換でない行列に $\beta$ を掛ける側を変えると、空間の $\gamma^i$ の項に符号差が出るので、Hamiltonian の逆と共変な方程式の逆を組み替えるときは順序を確認する。両方の極を下に置く遅延条件では、時間に戻すと $P_+ e^{-iEt/\hbar}+P_- e^{+iEt/\hbar}$ の両方を同じステップ関数で使う。量子場で時間順序を指定する Feynman の境界値は、ここでの初期値からの遅延と同じではない。

:::example[あるスピノルに対する二つの留数]
$p_z=\sqrt3 mc$、上向きスピンの上成分を入力にする。この方向で $\beta$ の期待値は 1、$\alpha_z$ は零であり、$E=2M$ なので対角の resolvent の留数は $(1+M/E)/2=3/4$ と $(1-M/E)/2=1/4$ である。和が 1 になるのは大きい $z$ での先頭が $1/z$ になるためである。静止 $p_z=0$ なら下の極の重みは零で、上成分だけの入力が正エネルギーの固有状態になる。別の運動量で同じ列を使っても枝の重みは一定ではない。
:::

### ヘリシティと横方向の偏極

固定した非零 $\vb p$ の正枝のスピン自由度を二成分 $\chi$ で表す。$\vb n$ が $\hat{\vb p}$ なら $\boldsymbol\Sigma\cdot\vb n$ は $H$ と可換だが、別の軸では一般に可換でない。実際、$[\boldsymbol\Sigma\cdot\vb n,\boldsymbol\alpha\cdot\vb p]=2i\boldsymbol\alpha\cdot(\vb n\times\vb p)$ である。正枝内を保つ偏極を測りたいときには、静止での Pauli 軸を同じ boost で運んだ射影を使う。非相対論のスピン軸を $4\times4$ の $\boldsymbol\Sigma$ に単に挿入したものとは一般に異なる。

:::example[垂直軸の期待値が運動量に依存する]
$\vb p=p\hat z$、$\chi=\ket{+x}$ とすると、正枝の $u$ の上部分では $\sigma_x$ の期待値は $+1$、下部分は $\sigma_z\chi=\ket{-x}$ なので $-1$ である。従って $u^\dagger\Sigma_x u=(E+M)/(2E)-(E-M)/(2E)=M/E$ である。一方、静止の $\ket{+x}$ を共変に運んだ偏極の射影にはこの状態が確率 1 で入る。$E=2M$ なら単純な $\Sigma_x$ の期待値は $1/2$ で、同じ「スピン」という語でも指定する演算子が違う。
:::

## 非相対論的極限と g 因子

電磁結合は :eqref[eq-aqm-kg-minimal] と同じ規約で $\vb\pi=\vb p-q\vb A$、$H=\beta mc^2+c\boldsymbol\alpha\cdot\vb\pi+q\Phi$ とする。$\psi=e^{-imc^2t/\hbar}(\varphi,\chi)^\mathsf T$ と置くと

$$
(i\hbar\partial_t-q\Phi)\varphi=c\boldsymbol\sigma\cdot\vb\pi\,\chi,\qquad
(2mc^2+i\hbar\partial_t-q\Phi)\chi=c\boldsymbol\sigma\cdot\vb\pi\,\varphi .
$$

エネルギー変化が $mc^2$ より小さければ $\chi\simeq(\boldsymbol\sigma\cdot\vb\pi)\varphi/(2mc)$ となる。$\comm{\pi_i}{\pi_j}=iq\hbar\epsilon_{ijk}B_k$ であり、Pauli 積の反対称部にはこの交換子だけが寄与する。よって

$$
(\boldsymbol\sigma\cdot\vb\pi)^2=\vb\pi^2-q\hbar\boldsymbol\sigma\cdot\vb B.
$$

上成分の方程式は

$$
i\hbar\partial_t\varphi=
\left(\frac{\vb\pi^2}{2m}+q\Phi-\frac q m\vb S\cdot\vb B\right)\varphi,
\qquad \vb S=\hbar\boldsymbol\sigma/2.
\label{eq-aqm-pauli-limit}
$$

$q\ne0$ として磁気モーメントの規約 $\vb\mu=gq\vb S/(2m)$ と比べると $g=2$ である。この値には輻射補正を含めていない。

:::example[電子の Zeeman 分裂]
$q=-e$、$\vb B=B\hat{\vb z}$ ならスピン項は $+\mu_B B\sigma_z$、$\mu_B=e\hbar/(2m)$ である。二準位差は $2\mu_B B$ であり、$\mu_B\simeq5.79\times10^{-5}\,\mathrm{eV\,T^{-1}}$ を用いると $B=1\,\mathrm T$ で約 $1.16\times10^{-4}\,\mathrm{eV}$ になる。
:::

:::example[正電荷の場合]
$q>0$、$B>0$ なら $s_z=+\hbar/2$ のスピン項は $-q\hbar B/(2m)$ である。軌道項は電荷の符号に応じて $\vb\pi$ を変えるが、スピン項の符号は電荷だけで直接反転する。
:::

### 磁場だけの場合の厳密な分散

$\Phi=0$、一様な $\vb B=B\hat{\vb z}$ とし、$H^2=M^2+c^2(\boldsymbol\alpha\cdot\vb\pi)^2$ を計算する。Pauli 積を上下に適用すると $(\boldsymbol\alpha\cdot\vb\pi)^2=\vb\pi^2-q\hbar B\Sigma_z$ である。:eqref[eq-aqm-kg-landau] の Landau 振動子の固有値を用いれば
$$
E_{n,s}^2=M^2+c^2p_z^2+\hbar c^2|qB|\big[2n+1-s\,\sgn(qB)\big],
\quad s=\pm1.
$$
ここで $s$ は平方化した二成分のラベルである。$H$ 自身では横運動が上下の成分を結ぶ。正枝は、上成分 $\varphi$ を $(\boldsymbol\sigma\cdot\vb\pi)^2$ の固有状態に取り、下を $c\boldsymbol\sigma\cdot\vb\pi\varphi/(E+M)$ として構成できる。

:::example[最低 Landau 準位のスピン]
$qB>0$、$p_z=0$ とする。$n=0,s=+1$ では横の係数が零なので $E=M$ である。$n=0,s=-1$ と $n=1,s=+1$ はどちらも係数が 2 で同じ平方根を持つ。非相対論的には、前者の最低準位で軌道の $\hbar\omega_c/2$ とスピンの $-\hbar\omega_c/2$ が相殺する。軌道とスピンの両方を入れて初めてこの零点が出る。
:::

### 弱い磁場のエネルギーを接線で読む

非相対論的な自由な二成分で運動量の交換子は $[\pi_i,\pi_j]=iq\hbar\epsilon_{ijk}B_k$ である。$\sigma_i\sigma_j=\delta_{ij}I+i\epsilon_{ijk}\sigma_k$ から、反対称な後半は $\tfrac i2\sum\epsilon_{ijk}\sigma_k[\pi_i,\pi_j]$ となる。$\sum_{ij}\epsilon_{ijk}\epsilon_{ij\ell}=2\delta_{k\ell}$ なので $-q\hbar\boldsymbol\sigma\cdot\vb B$ を得る。$\vb S=\hbar\boldsymbol\sigma/2$ と書けば磁気項は $-(q/m)\vb S\cdot\vb B$、磁気能率を $g(q/(2m))\vb S$ と定義したとき $g=2$ である。軌道の磁気能率を同じ定義で求めるには対称ゲージ $\vb A=\vb B\times\vb r/2$ を使い、$\vb B$ に一次の $(\vb p-q\vb A)^2/(2m)$ は $-q\vb L\cdot\vb B/(2m)$ になる。$\vb A\cdot\vb p=\vb B\cdot(\vb r\times\vb p)/2$ のためであり、軌道側は係数 1 である。

:::example[軌道とスピンの補正が逆の寄与を持つ準位]
一様な弱場の $z$ 方向に対し、電子の $\ell=1,m_\ell=1$、$s_z=-\hbar/2$ の積状態を使う。スピン軌道相互作用を無視できる模型で、磁気エネルギーの一次は $\mu_BB(m_\ell+2m_s)=0$ となる。$m_s=+1/2$ に替えれば $2\mu_BB$ になる。これは積基底での期待値であり、スピン軌道相互作用が大きく固定 $J$ が良い量子数のときは結合基底に磁気演算子を射影して求める。
:::

### Dirac 方程式を二乗すると磁場と電場の項が残る

$D_\mu=\partial_\mu+iq A_\mu/\hbar$、$\pi_\mu=i\hbar D_\mu$ とする。$[\pi_\mu,\pi_\nu]=-iq\hbar(\partial_\mu A_\nu-\partial_\nu A_\mu)=-iq\hbar F_{\mu\nu}$ である。$\gamma^\mu\gamma^\nu=\eta^{\mu\nu}I-i\sigma^{\mu\nu}$ を使い、後者が反対称なので
$$
(\gamma^\mu\pi_\mu)^2=\pi_\mu\pi^\mu-\frac{q\hbar}{2}\sigma^{\mu\nu}F_{\mu\nu}
$$
となる。これは最小結合を入れた後に二乗した結果であり、共変運動量同士がもう可換ではないことが効く。$A^\mu=(\Phi/c,\vb A)$ の規約で純磁場 $\vb B=B\hat z$ なら $F_{12}=-B$ なので、二項目は $+q\hbar\Sigma_z B$ となり、時間側の $E^2$ へ移せば先の Landau 分散の $-q\hbar c^2\Sigma_z B$ と一致する。電場があると $\sigma^{0i}=i\alpha_i$ に結合する項も出る。二乗した方程式は一階の Dirac 条件を含意しないので、スペクトルを求めた後で成分間の式を満たすスピノルを構成する。

:::example[同じ運動量の二成分で磁気項の符号を追う]
$q>0$、$\vb B=B\hat z$、同じ軌道振動子 $n$ を固定する。上向き $s=+1$ と下向き $s=-1$ で平方エネルギーの差は $2q\hbar c^2B$ で下向きが上である。低場なら二平方根の差は $2q\hbar c^2B/(2M)=q\hbar B/m$ となる。Pauli のスピンエネルギーの差と一致する。平方の段階の差にはエネルギー二乗の単位があり、そのまま分光線のエネルギー差とはならない。
:::

### 外場が一つの保存電流と結びつく

最小結合した Dirac の作用で $\partial_\mu\to D_\mu$ とすると、作用の相互作用部分は $-qc\,\bar\psi\gamma^\mu A_\mu\psi$ である。$D_\mu$ を微分とゲージ場に分ければ、もとの微分の二項で同符号の電位が出る。$j^\mu=c\bar\psi\gamma^\mu\psi$ を使って $-qj^\mu A_\mu$ になり、Hamiltonian の相互作用は $q\Phi-qc\boldsymbol\alpha\cdot\vb A$ と一致する。$A_\mu\mapsto A_\mu-\partial_\mu\Lambda$、$\psi\mapsto e^{iq\Lambda/\hbar}\psi$ で $D_\mu\psi$ は同じ位相を得るので全作用が不変である。場の変換をしないで外部ポテンシャルだけをゲージ変換した場合、電流との積の変化は作用の境界と $q\Lambda\partial_\mu j^\mu$ になる。保存電流を使う限り、境界で消えるゲージの変更は一次の遷移振幅に寄与しない。

:::example[平面波間の縦なゲージ成分を除く]
同じ質量の二つの正枝で四元運動量の移行を $Q_\mu=p'_\mu-p_\mu$ とする。移動する波の電流は $\bar u'\gamma^\mu u$ に比例する。$Q_\mu\bar u'\gamma^\mu u=\bar u'(\not p'-\not p)u=mc\,\bar u'u-mc\,\bar u'u=0$ である。従って外場の Fourier 成分に $Q_\mu$ に比例する項を加えても、行列要素は変わらない。例えば静的な純ゲージの縦ベクトルポテンシャルは、両枝が同じエネルギーを持つ条件で散乱に寄与しない。積の中で $\not p'$ は左の bra、$\not p$ は右の ket に作用させる必要がある。
:::

:::example[最低準位との差を磁場の二次まで比べる]
$qB>0,p_z=0$ の分散で、横エネルギーが零の枝と横の係数 2 の枝を比較する。差は $\sqrt{M^2+2q\hbar c^2B}-M$ である。$\hbar\omega_c=q\hbar B/m$ とすれば展開は $\hbar\omega_c-(\hbar\omega_c)^2/(2M)+\cdots$ となる。$\hbar\omega_c/M=1/100$ のとき、非相対論的な先頭に対する二次の相対補正は $-1/200$ である。軌道とスピンを一度合わせた平方エネルギーから平方根を取るので、各項の非相対論的な結果に独立な百分率を掛けて足す方法ではこの補正を得られない。
:::

## Foldy–Wouthuysen 変換とスピン軌道相互作用

大成分・小成分の混合をユニタリに取り除けば、正枝の内積を保ちながら高次の有効 Hamiltonian を求められる。静的場で $H=\beta M+\mathcal E+\mathcal O$、$M=mc^2$ とし、偶演算子は $[\beta,\mathcal E]=0$、奇演算子は $\{\beta,\mathcal O\}=0$ とする。$A=\beta\mathcal O/(2M)$ は反 Hermite であり、$e^A H e^{-A}$ の一次交換子は $[A,\beta M]=-\mathcal O$ なので元の奇項を消す。このブロック対角化を :term[Foldy–Wouthuysen 変換]{en="Foldy–Wouthuysen transformation" yomi="ふぉるでぃーうーすはいぜんへんかん"} と呼ぶ。

交換子展開の偶項は $[A,\mathcal O]/2=\beta\mathcal O^2/(2M)$ となる。$[A,\mathcal E]=\beta[\mathcal O,\mathcal E]/(2M)$ は奇で、次の小回転で消す。同じ展開の $\frac12[A,[A,\mathcal E]]$ は $-[\mathcal O,[\mathcal O,\mathcal E]]/(8M^2)$ となる。新しい奇項の二乗から出る偶補正は $M^{-3}$ 以降である。静電場 $V=q\Phi$ と $\vb A=0$ の場合は

$$
[\mathcal O,[\mathcal O,V]]
=-\hbar^2c^2\nabla^2 V-2\hbar c^2\boldsymbol\Sigma\cdot(\grad V\times\vb p).
$$

最初の交換子は $-i\hbar c\boldsymbol\alpha\cdot\grad V$ である。次は $\alpha_i\alpha_j=\delta_{ij}I+i\epsilon_{ijk}\Sigma_k$ を用い、二階微分の対称成分と運動量を残す反対称成分に分ける。

自由奇項の高次補正は正枝の平方根 $\sqrt{M^2+\mathcal O^2}$ を展開して求まる。したがって低速 $p/(mc)\ll1$ で、静電場の勾配の補正を含めた正枝は

$$
H_+=mc^2+\frac{p^2}{2m}+V-\frac{p^4}{8m^3c^2}
 +\frac{\hbar^2}{8m^2c^2}\nabla^2 V
 +\frac{\hbar}{4m^2c^2}\boldsymbol\sigma\cdot(\grad V\times\vb p)
\label{eq-aqm-fw}
$$

となる。ここでは $c^{-2}$ までを残した。第五項を :term[Darwin 項]{en="Darwin term" yomi="だーうぃんこう"}、第六項を :term[スピン軌道相互作用]{en="spin–orbit interaction" yomi="すぴんきどうそうごさよう"} という。中心力では $\grad V\times\vb p=(V'(r)/r)\vb L$ なので係数は $V'(r)\vb L\cdot\vb S/(2m^2c^2r)$ となる。

:::example[静電調和井戸]
$V=m\omega^2 r^2/2$ なら $\nabla^2V=3m\omega^2$ である。Darwin 項は $3\hbar^2\omega^2/(8mc^2)$ の一様シフトを与える。スピン軌道の係数 $V'/r=m\omega^2$ も一定で、角運動量の結合で対角化できる。
:::

:::example[Coulomb 場の符号]
$V=-g/r$、$g>0$ では $V'/r=g/r^3>0$ なので、同じ $\ell>0$ なら $\vb L\cdot\vb S$ が正の $j=\ell+1/2$ でこの項は上方にずらす。Darwin 項は $\nabla^2V=4\pi g\delta^{(3)}(\vb r)$ により $s$ 波に作用する。
:::

### 自由 Dirac Hamiltonian の厳密な分離

自由の各運動量で $\mathcal O=c\boldsymbol\alpha\cdot\vb p$、$E=\sqrt{M^2+c^2p^2}$ を用いる。$\beta\mathcal O$ は反 Hermite で $(\beta\mathcal O)^2=-c^2p^2 I$ なので、

$$
U=\frac{E+M+\beta\mathcal O}{\sqrt{2E(E+M)}}
$$

はユニタリである。分子と随伴の積は $(E+M)^2+c^2p^2=2E(E+M)$ になるからである。さらに $(E+M+\beta\mathcal O)(\beta M+\mathcal O)=E[(E+M)\beta+\mathcal O]$ で、$UHU^\dagger=\beta E$ を得る。計算は $\mathcal O\beta=-\beta\mathcal O$ と $\mathcal O^2=(E^2-M^2)I$ を使う。$|\vb p|\ll mc$ なら $U=I+\beta\mathcal O/(2M)+\order{p^2}$ なので、静的場で用いた最初の小変換と一致する。$E=2M$ なら $U$ の対角の成分は $\sqrt3/2$、上下混合の絶対値は $1/2$ で、元のスピノルの成分比を規格化した回転に当たる。

### 偶なポテンシャルと厳密な自由回転

自由の $U$ は運動量に依存するので、位置に依存する $V$ があると $UVU^\dagger=V$ とは限らない。低速で $U=I+A+A^2/2+\cdots$、$A=\beta\mathcal O/(2M)$ として展開すると、最初の補正 $[A,V]$ は奇であり、次の $\tfrac12[A,[A,V]]$ が偶である。係数を順に動かすと $[A,[A,V]]=-[\mathcal O,[\mathcal O,V]]/(4M^2)$ である。最初の奇補正を除く次の発生子は $[A,V]/M$ の次数なので、この除去の偶補正は $c^{-2}$ の静電勾配項より高い。このように運動とポテンシャルの交換子が、古典分散の展開には見えない二項を生む。

:::example[一方向で二重交換子の符号を確認する]
$V(x)=m\omega^2x^2/2$、$\mathcal O=c\alpha_x p_x$ とする。$[\mathcal O,V]=-i\hbar c\alpha_x m\omega^2x$ である。次の交換子では $\alpha_x^2=I$、$[p_x,x]=-i\hbar$ により $-\hbar^2c^2m\omega^2I$ となる。$-1/(8M^2)$ を掛けると Darwin シフト $+\hbar^2\omega^2/(8mc^2)$ である。$\grad V$ と運動が同じ向きなのでこの一次元計算にスピン軌道項はない。三方向に同じ井戸を置く場合は Laplacian の寄与が三倍になる。
:::

### 静電項のスペクトルでブロック変換を検証する

座標の静電項を運動量空間に移すと $\mel{\vb p'}{V}{\vb p}=\widetilde V(\vb p'-\vb p)$ である。ただし Fourier 規格化の係数をこの記号に含める。自由に正枝へ写す規格化スピノル $u$ を左右につけると、正枝に一次の行列要素は
$$
\widetilde V(\vb p'-\vb p)N' N
\left[I+\frac{c^2(\boldsymbol\sigma\cdot\vb p')(\boldsymbol\sigma\cdot\vb p)}{(E'+M)(E+M)}\right]
$$
となる。$N=\sqrt{(E+M)/(2E)}=1-p^2/(8m^2c^2)+\cdots$ と分母 $E+M\simeq2M$ を使い、$c^{-2}$ で展開すると
$$
\widetilde V(\vb p'-\vb p)\left[I-\frac{|\vb p'-\vb p|^2}{8m^2c^2}
+\frac{i\boldsymbol\sigma\cdot(\vb p'\times\vb p)}{4m^2c^2}\right]
$$
となる。実空間では $\nabla^2V$ が $-|\vb p'-\vb p|^2/\hbar^2$、$\grad V\times\vb p$ が $i(\vb p'-\vb p)\times\vb p/\hbar$ を与えるので :eqref[eq-aqm-fw] と一致する。これは $V$ に一次の静電行列要素の検算であり、強い場で自由射影を課せば全ての束縛スペクトルが得られるという仮定は置いていない。

:::example[向かい合う二運動量で Darwin 項だけを残す]
$\vb p'=-\vb p$ とすると交差積が零で、係数は $1-p^2/(2m^2c^2)$ になる。同じ $\vb p'=\vb p$ なら補正は零であり、均一なポテンシャルが内部の運動量依存の結合を作らないことと合う。$\vb p'=p\hat x,\vb p=p\hat y$ では交差積が $p^2\hat z$ になるため、$z$ 偏極の行列要素に $\pm i p^2/(4m^2c^2)$ の差が現れる。運動量交換の方向を逆にすれば虚部が共役になり、Hermite 性が保たれる。
:::

### 時変な変換で電場の項の一部を残す

$\psi_{\rm new}=U(t)\psi$ とすると、新しい Hamiltonian は $UHU^\dagger+i\hbar\dot U U^\dagger$ である。$U=e^{A(t)}$、$A=\beta\mathcal O/(2M)$ とする。指数を二次まで展開すると
$$
\dot U U^\dagger=\dot A+\tfrac12[A,\dot A]+\cdots
$$
であり、$\dot A$ は奇、$[A,\dot A]=-[\mathcal O,\dot{\mathcal O}]/(4M^2)$ は偶である。場の時間変化を速過ぎないものとして $c^{-2}$ の次数を保つと、この偶寄与は
$$
-\frac1{8M^2}[\mathcal O,[\mathcal O,\mathcal E]+i\hbar\dot{\mathcal O}]
$$
に静的な二重交換子と合わせられる。残る奇項を次の小変換で除くとき、その二乗はこの次数より高いと仮定する。時間変化が速く $\hbar\partial_t$ の尺度が $mc^2$ に近い場合はこの整理を使わない。

例えば $\mathcal E=0$、空間で一定の $\vb A(t)$ では $\dot{\mathcal O}=qc\boldsymbol\alpha\cdot\vb E$ である。$[\alpha_i,\alpha_j]=2i\epsilon_{ijk}\Sigma_k$ を用いると $[\mathcal O,\dot{\mathcal O}]=2iqc^2\boldsymbol\Sigma\cdot(\vb\pi\times\vb E)$ である。従って正枝の偶補正は $-q\hbar\boldsymbol\sigma\cdot(\vb E\times\vb\pi)/(4m^2c^2)$ であり、静電場で $\grad V=-q\vb E$ を入れた符号と同じになる。これは空間一様な場での照合である。空間変化があっても、$\mathcal E=q\Phi$ と $\mathcal O=c\boldsymbol\alpha\cdot\vb\pi$ なら内側の和は $[\mathcal O,\mathcal E]+i\hbar\dot{\mathcal O}=i\hbar qc\boldsymbol\alpha\cdot\vb E$ となる。そこで $[\pi_i,E_j]=-i\hbar\partial_i E_j$ と Pauli 積を使うと、外側の交換子は
$$
q\hbar^2c^2\grad\cdot\vb E
+q\hbar c^2\boldsymbol\Sigma\cdot(\vb E\times\vb\pi-\vb\pi\times\vb E)
$$
である。第一項は $\alpha_i\alpha_j$ の対称成分から、第二項は反対称成分から生じる。従って正枝の電場補正は
$$
-\frac{q\hbar^2}{8m^2c^2}\grad\cdot\vb E
-\frac{q\hbar}{8m^2c^2}\boldsymbol\sigma\cdot(\vb E\times\vb\pi-\vb\pi\times\vb E).
$$
順序を保つと $(\vb E\times\vb\pi)^\dagger=-\vb\pi\times\vb E$ であり、第二項も Hermite になる。磁場を含む運動項は $\mathcal O^2$ と $\mathcal O^4$ から別に得られる。

:::example[同じ電場を二つのゲージで見る]
空間一定の $\vb E=E_x\hat x$ を短い時間で考える。$\Phi=-E_x x,\vb A=0$ のゲージでは静電の二重交換子が $\sigma_z$ に比例する項を生む。$\Phi=0,\vb A=-E_x t\hat x$ で同じ $\vb E$ を表すなら、同じ項は $\dot{\mathcal O}$ を含む時間依存の変換から出る。瞬間の $\pi_y=p_y$ を固定すると $\delta H=-q\hbar E_xp_y\sigma_z/(4m^2c^2)$ である。電位の勾配だけを残した式を後者のゲージに移すと、この電場の効果を落とすことになる。
:::

:::example[自由な平方根の補正の規模を見積もる]
$p=mc/2$ の正枝では $E-M=(\sqrt5/2-1)M\simeq0.1180M$ である。非相対論的な $p^2/(2m)$ は $0.125M$、$-p^4/(8m^3c^2)$ を加えると $(1/8-1/128)M=0.1171875M$ となる。補正を入れた後の差は約 $8.46\times10^{-4}M$ である。自由な場合は $p/(mc)$ がこの展開の精度を決める。静電場がある場合はその運動量の尺度に加えて勾配による二項の大きさも確かめる。
:::

## 演習問題

::::exercise{level=basic}
$\gamma^i$ の二乗が $-I$ であることを $\beta,\alpha_i$ から示せ。
:::solution
$\gamma^i=\beta\alpha_i$ で $\alpha_i\beta=-\beta\alpha_i$ なので $(\gamma^i)^2=-\beta^2\alpha_i^2=-I$ である。空間計量の負符号と一致する。
:::
::::

::::exercise{level=basic}
静止した正枝で $\chi=(1,0)^\mathsf T$ の $j^\mu$ と $\bar uu$ を求めよ。
:::solution
$u=(1,0,0,0)^\mathsf T$ なので $\bar uu=u^\dagger\beta u=1$ である。$\alpha_i$ は上下を交換するので空間電流は零であり、時間成分は $c$ となる。
:::
::::

::::exercise{level=basic}
$p_z>0$ で $\chi=(0,1)^\mathsf T$ の正枝のヘリシティを求めよ。
:::solution
$\sigma_z\chi=-\chi$ で、下成分も $\sigma_z\chi$ に比例して同じ固有値を持つ。従って $\Sigma_z/2$ の固有値は $-1/2$ である。
:::
::::

::::exercise{level=basic}
自由の正枝と負枝の射影 $P_\pm=(I\pm H/E)/2$ が直交することを示せ。
:::solution
積は $(I-H^2/E^2)/4=0$ である。和は $I$、各二乗は $H^2=E^2I$ により元に等しい。$H^\dagger=H$ なのでどちらも直交射影である。
:::
::::

::::exercise{level=advanced}
$u$ を :eqref[eq-aqm-dirac-u] の規格化とする。$\bar uu$ を求めよ。
:::solution
上成分のノルム二乗は $(E+mc^2)/(2E)$、下成分は $(E-mc^2)/(2E)$ である。$\beta$ は下の符号を反転するので差 $mc^2/E$ が得られる。Lorentz スカラーを 1 に揃える別の規格化では $u$ に $\sqrt{E/(mc^2)}$ を掛ける。
:::
::::

::::exercise{level=advanced}
$\mathcal O$ が奇なら $[A,\beta M]=-\mathcal O$、$[A,\mathcal O]=\beta\mathcal O^2/M$ を証明せよ。ただし $A=\beta\mathcal O/(2M)$ とする。
:::solution
$\beta\mathcal O\beta=-\mathcal O$ なので $A\beta M-\beta M A=(-\mathcal O-\mathcal O)/2=-\mathcal O$ である。次の交換子では $\mathcal O\beta\mathcal O=-\beta\mathcal O^2$ を使うと $2\beta\mathcal O^2/(2M)$ となる。
:::
::::

::::exercise{level=advanced}
$V=m\omega^2r^2/2$、$\ell=1$ のスピン軌道項だけによる $j=3/2$ と $j=1/2$ のエネルギー差を求めよ。
:::solution
係数は $\omega^2/(2mc^2)$ である。$\vb L\cdot\vb S=(\vb J^2-\vb L^2-\vb S^2)/2$ の固有値は $\hbar^2/2,-\hbar^2$ なので、差は $3\hbar^2\omega^2/(4mc^2)$ である。
:::
::::

::::exercise{level=advanced}
$\gamma^5$ と Dirac 方程式について、質量が零なら二つのカイラリティを独立に解けることを示せ。
:::solution
$P_\pm\gamma^\mu=\gamma^\mu P_\mp$ である。$m=0$ の式に $P_\pm$ を左から掛けると $\gamma^\mu\partial_\mu(P_\mp\psi)=0$ を得る。$m\ne0$ なら右辺に $mcP_\pm\psi$ が残るため、両成分は結合する。
:::
::::

::::exercise{level=advanced}
自由 Dirac 粒子の $p_z=3mc/4$ で、初期を Dirac 表現の $\psi_0=(1,0,0,0)^\mathsf T$ とする。正負の枝への射影の重み、下成分の時間依存、平均エネルギーを求め、三つの結果の整合性を確認せよ。

:::solution
$M=mc^2$ とすれば $E=5M/4$ であり、スピン上向きは質量項と $p_z$ の運動項によって他のスピンへ移らない。上向きの上下二成分の $H$ を $M\sigma_z+(3M/4)\sigma_x$ と書ける。$P_\pm=(I\pm H/E)/2$ を初期の上成分に挟むと、重みは $(1\pm4/5)/2$ なので $9/10,1/10$ である。従って平均エネルギーは $E(9/10-1/10)=M$ である。最初の期待値を $H$ で直接求めても、$\sigma_z$ は $+1$、$\sigma_x$ は零なので同じになる。

時間発展を $\cos(Et/\hbar)I-i(H/E)\sin(Et/\hbar)$ と書けば、下成分は $-(3i/5)\sin(Et/\hbar)$、上成分は $\cos(Et/\hbar)-(4i/5)\sin(Et/\hbar)$ である。確率の和は $\cos^2+(16/25+9/25)\sin^2=1$ である。下成分の最大確率 $9/25$ は最初の負枝の重み $1/10$ と異なる。枝は $H$ の射影であり、上下成分は $\beta$ の射影である。前者は自由の発展で保存され、後者は一般に二枝の干渉を含む。
:::
::::

::::exercise{level=advanced}
静電項の一つの Fourier 成分で $\vb p=p\hat{\vb y}$ から $\vb p'=p\hat{\vb x}$ への正枝の行列要素を求める。$z$ のスピンがそれぞれ上・下の場合の :eqref[eq-aqm-fw] の二次までの係数を求め、等重みのスピン平均で虚の項がどうなるか調べよ。

:::solution
運動量の移行は $\vb Q=\vb p'-\vb p$ であり、$\abs{\vb Q}^2=2p^2$ である。空間 Fourier の規約により $\nabla^2V$ の係数は $-\abs{\vb Q}^2/\hbar^2$ なので、Darwin の項は $\widetilde V$ に $-p^2/(4m^2c^2)$ を掛ける。$(\grad V)\times\vb p$ の行列要素は $i\widetilde V(\vb Q)(\vb p'\times\vb p)/\hbar$ であり、$\vb p'\times\vb p=p^2\hat{\vb z}$ である。よって静電項を含めた係数は
$$
1-\frac{p^2}{4m^2c^2}\ \pm\ i\frac{p^2}{4m^2c^2}.
$$
ここでは自由な運動エネルギーは移行 $\vb p'\ne\vb p$ の行列要素に寄与しない。スピンの等重みの平均を取ると虚の寄与は打ち消し、実の補正だけが二次までに残る。個別の絶対値二乗を平均しても、虚の係数は二乗で四次にしか寄与しない。相対スピン位相を観測する場合には、強度の平均とは違いこの虚の項を残す。$p=mc/5$ なら二次の各補正の大きさは $1/100$ である。
:::
::::

## 確認テスト

::::quiz
計量 $(+,-,-,-)$ の Dirac 密度 $j^0/c$ を選べ。

- [ ] $\psi^\dagger\beta\psi$
- [x] $\psi^\dagger\psi$
- [ ] $\psi^\dagger\gamma^5\psi$

:::explanation
$j^0=c\bar\psi\gamma^0\psi=c\psi^\dagger(\gamma^0)^2\psi$ であり $(\gamma^0)^2=I$ である。従ってノルム二乗となる。
:::
::::

::::quiz{answer="2" tol="1e-9"}
最小電磁結合した Dirac 粒子の非相対論的なスピンの $g$ 因子を求めよ。
:::explanation
Pauli Hamiltonian のスピン結合は $-(q/m)\vb S\cdot\vb B$ である。定義 $-gq\vb S\cdot\vb B/(2m)$ と係数を比べると $g=2$ になる。
:::
::::

::::quiz
質量が非零の自由粒子で、固定運動量の Hamiltonian と可換なものを選べ。

- [x] $\hat{\vb p}\cdot\boldsymbol\Sigma$
- [ ] 一般の軸に沿う $\boldsymbol\Sigma$ の成分
- [ ] $\gamma^5$

:::explanation
$\hat{\vb p}\cdot\boldsymbol\Sigma$ は運動項の Pauli 行列と同じ軸なので可換で、$\beta$ とも可換である。一般の軸では Pauli 行列の交換子が残り、$\gamma^5$ は質量の $\beta$ と反交換する。
:::
::::

::::quiz{answer="0.25" tol="1e-9"}
$E=2mc^2$ の正枝の規格化スピノルで、下二成分のノルム二乗を求めよ。
:::explanation
下二成分の重みは $(E-mc^2)/(2E)$ である。$E=2mc^2$ を入れると $mc^2/(4mc^2)=1/4$ となる。
:::
::::

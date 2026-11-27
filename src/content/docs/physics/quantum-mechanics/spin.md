---
title: スピンと角運動量の合成
description: 二成分のスピン状態を磁気測定に結びつけ、角運動量を合成して多重項とスピン軌道相互作用を求める。
notation:
  - symbol: '\boldsymbol\sigma'
    meaning: Pauli 行列
  - symbol: '\gamma'
    meaning: 符号を含む磁気回転比
  - symbol: 'C_{m_1m_2}^{JM}'
    meaning: Clebsch–Gordan 係数
---

## この章の目標

- 任意方向のスピン測定と磁場中の時間発展を計算する。
- 直積基底から全角運動量の固有基底を構成する。
- スピン軌道相互作用の符号と準位間隔を判定する。

## スピン 1/2 と Pauli 行列

銀原子線の二分裂を記述するには、位置座標だけでは指定できない二状態が要る。
電子の :term[スピン]{en="spin" yomi="すぴん"} は内部角運動量である。
$s=1/2$ とすると :ref[thm-qm-ladder] から $S_z$ の固有値は $\pm\hbar/2$。
対応する $\ket\uparrow,\ket\downarrow$ をこの順に基底とすれば $S_+=\hbar\ketbra\uparrow\downarrow$ なので $S_i=\hbar\sigma_i/2$ で

$$
\sigma_x=\begin{pmatrix}0&1\\1&0\end{pmatrix},\quad \sigma_y=\begin{pmatrix}0&-i\\i&0\end{pmatrix},\quad \sigma_z=\begin{pmatrix}1&0\\0&-1\end{pmatrix}.
$$

これらを :term[Pauli 行列]{en="Pauli matrices" yomi="ぱうりぎょうれつ"} と呼び、以下では $X,Y,Z$ とも書く。積を取ると $\sigma_i^2=I$、$\sigma_x\sigma_y=i\sigma_z$ とその巡回置換、逆順では負になる。
よって $(\vb a\cdot\boldsymbol\sigma)(\vb b\cdot\boldsymbol\sigma)=(\vb a\cdot\vb b)I+i(\vb a\times\vb b)\cdot\boldsymbol\sigma$。
特に単位方向 $\vb n$ の $\vb n\cdot\boldsymbol\sigma$ の二乗は一であり、射影は $P_\pm=(I\pm\vb n\cdot\boldsymbol\sigma)/2$ である。

$\vb n=(\sin\theta\cos\phi,\sin\theta\sin\phi,\cos\theta)$ の上向きケットは $\cos(\theta/2)\ket\uparrow+e^{i\phi}\sin(\theta/2)\ket\downarrow$ である。
行列の対角項 $\pm\cos\theta$ と非対角項 $e^{\mp i\phi}\sin\theta$ を掛けて確かめられる。
確率 $p_k$ で規格化されたケット $\ket{\psi_k}$ を準備した集団は、$\rho=\sum_k p_k\ketbra{\psi_k}{\psi_k}$ という :term[密度行列]{en="density matrix" yomi="みつどぎょうれつ"} で表す。各準備の期待値を平均すれば $\langle A\rangle=\tr(\rho A)$ になる。$I,\sigma_x,\sigma_y,\sigma_z$ はエルミートな $2\times 2$ 行列の実基底なので、跡が一の $\rho$ は $\rho=(I+\vb s\cdot\boldsymbol\sigma)/2$ と一意に書ける。$\vb s$ は :term[Bloch ベクトル]{en="Bloch vector" yomi="ぶろっほべくとる"} である。
固有値 $(1\pm|\vb s|)/2$ が非負である条件は $|\vb s|\le1$ であり、階数一の純粋状態で等号となる。

:::example[分析器を傾ける]
$z$ 上向きを準備し $\theta=60^\circ$ の軸を測ると、上の確率は $\cos^2(30^\circ)=3/4$ である。
$\theta=90^\circ$ では半々となる。
一般の状態での上の確率は $\tr(\rho P_+)=(1+\vb s\cdot\vb n)/2$ である。
:::

:::example[一周の符号]
回転演算子のべき級数を偶奇に分けると $U=\cos(\theta/2)I-i\sin(\theta/2)\vb n\cdot\boldsymbol\sigma$。
$\theta=2\pi$ なら $-I$、$4\pi$ なら $I$ である。
:::

### スピノルの回転と三次元の回転

スピンの回転 $U$ を作用させると、$U(\vb v\cdot\boldsymbol\sigma)U^\dagger=(R_U\vb v)\cdot\boldsymbol\sigma$ により実ベクトルの回転 $R_U$ が定まる。Pauli 行列の積の公式から、$U=e^{-i\theta\vb n\cdot\boldsymbol\sigma/2}$ に対し

$$
R_U\vb v=\vb v\cos\theta+(\vb n\times\vb v)\sin\theta+
\vb n(\vb n\cdot\vb v)(1-\cos\theta)
$$

となる。第一項と第三項を合わせると軸に平行な成分はそのままであり、軸に垂直な成分は $\vb v_\perp\cos\theta+(\vb n\times\vb v_\perp)\sin\theta$ となる。これは [角運動量の章](/physics/quantum-mechanics/angular-momentum/) で構成した $\SU(2)$ から $\SO(3)$ への二重被覆の具体式である。

$U=-I$ の全体位相は一経路の測定確率からは読めない。例えば経路を $(\ket a+\ket b)/\sqrt2$ に準備して $b$ にだけ $2\pi$ の回転を施せば、再合成の直前には経路が $(\ket a-\ket b)/\sqrt2$ になる。この二つの経路状態は直交するため、同じ分析器で左右の出力が入れ替わる。

### 三軸の測定による状態の再構成

同じ状態を多数回準備して $\sigma_x,\sigma_y,\sigma_z$ を別々の試料で測る。
上向きの確率 $p_i$ から $s_i=2p_i-1$ が得られる。
行列で書くと

$$
\rho=\tfrac12\begin{pmatrix}1+s_z&s_x-is_y\\s_x+is_y&1-s_z\end{pmatrix}.
$$

$z$ 測定だけでは対角分布は決まるが、非対角成分の位相は決まらない。
例えば $p_x=1/2,p_y=1,p_z=1/2$ なら $\vb s=(0,1,0)$ であり、状態は $\ket{+y}$ である。
各測定が別の試料で行われるのは、最初の測定自身がその試料の後続統計を変えるためである。

### スピノルの成分から測定方向を組み立てる

規格化した二成分列を $(a,b)^\mathsf T$ とする。$a\ne0$ なら全体位相を選んで $a\ge0$ とし、$a=\cos(\theta/2)$、$b=e^{i\phi}\sin(\theta/2)$ と書ける。$a=0$ では南極に対応して $\phi$ は不要である。外積の対角は $\cos^2(\theta/2),\sin^2(\theta/2)$、非対角は $e^{\mp i\phi}\sin\theta/2$。従って密度行列は $(I+\vb n\cdot\boldsymbol\sigma)/2$、$\vb n=(\sin\theta\cos\phi,\sin\theta\sin\phi,\cos\theta)$ になる。球面上の一点と純粋な射線が対応し、上下二成分の複素数を三次元の古典的な角運動量の成分と見なす必要はない。

二つの純粋状態の重なりは、対応する射影の積の跡で求められる。$\tr\sigma_i=0,\ \tr(\sigma_i\sigma_j)=2\delta_{ij}$ より
$$
|\braket{+\vb m}{+\vb n}|^2=\tr(P_{\vb m}P_{\vb n})
=\frac{1+\vb m\cdot\vb n}{2}=\cos^2(\gamma/2).
$$
角 $\gamma$ が直角な二軸の上向きは直交しない。直交する二つの状態は球上で正反対に位置する。同じ計算は $\rho=(I+\vb s\cdot\boldsymbol\sigma)/2$ の混合試料にも使え、読みの確率は $(1+\vb s\cdot\vb m)/2$ である。混合状態は球の内部、純粋状態はその表面を占める。


### 三方向の計数からスピノルの位相を決める

同じ方法で準備した試料を三群に分け、各群で $x,y,z$ のどれか一軸を測る。正の出口の確率を $p_x,p_y,p_z$ とすると $r_i=2p_i-1$ は Bloch ベクトルの成分である。例えば $p_x=0.74,\ p_y=0.82,\ p_z=0.80$ なら $\vb r=(0.48,0.64,0.60)$ で $|\vb r|^2=1$。この場合は純粋ケットを復元できる。上の成分 $a=\sqrt{p_z}$ を実正に選ぶと、密度行列の左下成分 $ba^*=(r_x+ir_y)/2$ より $b=(0.48+0.64i)/(2\sqrt{0.8})$。従って下の成分の確率は $0.20$、相対位相の正接は $4/3$ となる。

$|\vb r|<1$ なら同じ操作で純粋スピノルを決めることはできない。測定で決まるのは固有値 $(1\pm|\vb r|)/2$ をもつ混合密度行列である。実測計数をそのまま上式に入れて $|\vb r|>1$ になった場合は、負の固有値を物理的に解釈する前に標本数と誤差を調べる。三軸は一粒子で順に測った値ではなく同じ準備の三群の統計であり、途中に状態更新を挟む実験と比較するにはその測定順を明示する。


### 位相のあるスピノルに対する異なる二軸の相関

二つの分析器を互いに直交しない方向に向け、単独で測った分散と、強い不確定性関係の共分散を区別する。

Bloch 表示での各軸の期待値は、状態ベクトルから作った Bloch ベクトルと測定方向との内積に $\hbar/2$ を掛けた量である。

各物理量を $S_z=a_0 A$、$S_n=b_0 B$ と分け、$a_0=5.2729\times10^{-35}\,\mathrm{J\,s}$、$b_0=5.2729\times10^{-35}\,\mathrm{J\,s}$ とする。計算に必要な無次元の列と行列は

$$
\ket\psi=\begin{pmatrix}\frac{\sqrt{3}}{3}\\\frac{\sqrt{3} \left(1 + i\right)}{3}\end{pmatrix},\qquad A=\begin{pmatrix}1 & 0\\0 & -1\end{pmatrix},\qquad B=\begin{pmatrix}\frac{4}{5} & \frac{3}{5}\\\frac{3}{5} & - \frac{4}{5}\end{pmatrix}.
$$

である。まず作用後の列を作り、元の状態に沿う成分を測る。

$$
\begin{aligned}A\ket\psi&=\begin{pmatrix}\frac{\sqrt{3}}{3}\\\frac{\sqrt{3} \left(-1 - i\right)}{3}\end{pmatrix},&\quad \mel\psi A\psi&=- \frac{1}{3},\\ B\ket\psi&=\begin{pmatrix}\frac{\sqrt{3} \left(7 + 3 i\right)}{15}\\\frac{\sqrt{3} \left(-1 - 4 i\right)}{15}\end{pmatrix},& \mel\psi B\psi&=\frac{2}{15}.\end{aligned}
$$

平均をスカラー倍した元の列を引けば、どちらも $\psi$ と直交する。これらの残りを $u,w$ とし、対角に分散、非対角に複素相関を含む小さなグラム行列を作る。

$$
u=\begin{pmatrix}\frac{4 \sqrt{3}}{9}\\\frac{2 \sqrt{3} \left(-1 - i\right)}{9}\end{pmatrix},\quad w=\begin{pmatrix}\frac{\sqrt{3} \left(19 + 9 i\right)}{45}\\\frac{\sqrt{3} \left(-5 - 14 i\right)}{45}\end{pmatrix},\quad \begin{pmatrix}\braket u u&\braket u w\\\braket w u&\braket w w\end{pmatrix}=\begin{pmatrix}\frac{8}{9} & \frac{38}{45} + \frac{2 i}{5}\\\frac{38}{45} - \frac{2 i}{5} & \frac{221}{225}\end{pmatrix}.
$$

実部は対称化した共分散、虚部は交換子の平均に対応する。積の順序を入れ替えれば相関の虚部だけが反転するため、次の二つを個別に照合できる。

$$
\frac{\langle AB+BA\rangle}{2}-\langle A\rangle\langle B\rangle=\frac{38}{45},\qquad \frac{\langle AB-BA\rangle}{2i}=\frac{2}{5}.
$$



| 読み | 単位目盛による表示 | SI での値 |
|---|---|---|
| 第一の標準偏差 | $\frac{2 \sqrt{2}}{3}\,a_0$ | $4.9713\times10^{-35}\,\mathrm{J\,s}$ |
| 第二の標準偏差 | $\frac{\sqrt{221}}{15}\,b_0$ | $5.2258\times10^{-35}\,\mathrm{J\,s}$ |

二つの幅の積が取り得る下限を比較するには、無次元の行列式を計算すればよい。

$$
\|u\|^2\|w\|^2-(\Re\braket u w)^2-(\Im\braket u w)^2=\frac{1768}{2025}-(\frac{38}{45})^2-(\frac{2}{5})^2=0.
$$

これは $\|\alpha u+\beta w\|^2\ge0$ に対応する条件である。ここでは二軸が非直交なので、反交換子の平均には二軸の内積に由来する項も入る。


## 磁場中のスピンと歳差運動

磁気モーメントを $\boldsymbol\mu=\gamma\vb S$ と書く。符号を含む $\gamma$ は :term[磁気回転比]{en="gyromagnetic ratio" yomi="じきかいてんひ"} である。
電子では電荷が負なので $\gamma=-g_se/(2m_e)<0$、$g_s\simeq2.00232$。
静的で一様な $\vb B=B\vb e_z$ について運動との結合を無視すると $H=-\gamma BS_z$、準位差は $\hbar|\gamma|B$ である。
$z$ 上向きの電子は磁気モーメントが磁場と逆向きで高エネルギー側にある。

$$
\begin{aligned}
U(t)&=\diag(e^{i\gamma Bt/2},e^{-i\gamma Bt/2}), \\
\frac{\dd\expval{\vb S}}{\dd t}&=\gamma\expval{\vb S}\times\vb B .
\end{aligned}
$$

運動方程式は $[S_i,S_j]=i\hbar\epsilon_{ijk}S_k$ を Heisenberg 方程式に入れたもの。
初期スピンが $+x$ なら、期待値の方位角は $-\gamma Bt$ となる。
この回転を :term[Larmor 歳差運動]{en="Larmor precession" yomi="らーもあさいさうんどう"} という。

:::example[一テスラの電子]
$\mu_{\rm B}=e\hbar/(2m_e)=5.7884\times10^{-5}\,\mathrm{eV\,T^{-1}}$ とすると、$B=1.00\,\mathrm T$ で準位差 $g_s\mu_{\rm B}B=1.16\times10^{-4}\,\mathrm{eV}$、通常の周波数は約 $28.0\,\mathrm{GHz}$ である。
:::

:::example[横磁場パルス]
$z$ 上向きを $B_x$ に時間 $t$ だけさらす。
回転の式で $z$ 下向きへの確率は $\sin^2(\gamma B_xt/2)$。
$|\gamma|B_xt=\pi$ にすれば完全に反転する。
:::

### 磁場の向きを急に変える

磁場の大きさを保って $z$ 軸から $x$ 軸へ短時間で向け直す。
切替直後にケットは元のままである。
切替後の固有基底では $\ket{\uparrow_z}=(\ket{+x}+\ket{-x})/\sqrt2$。
二成分に逆向きの位相が掛かり、再び $z$ 軸で測る下向き確率は $\sin^2(\gamma Bt/2)$ である。
一方、向きをゆっくり変える場合には瞬間固有状態への追従を検討する必要があり、この突然切替の初期条件を用いることはできない。

:::example[途中で磁場を反転させる]
一定の $z$ 磁場で時間 $\tau$ 歳差させた後、同じ大きさの逆向き磁場で時間 $\tau$ 歳差させる。
二つのユニタリは $U_{-B}(\tau)=U_B(\tau)^\dagger$ なので積は $I$。
この相殺にはずれを含む全磁場の反転が必要である。反転しない背景磁場や途中で変動した磁場に対しては、二つの発展は一般には逆にならない。
:::

### 磁場反転と時間反転

古典軌道を時間の逆順にたどれば、位置は同じで運動量と角運動量の向きが反転する。
スピン半分に対する時間反転は、$z$ 基底での複素共役を $K$ として $\Theta=-i\sigma_yK$ と表せる。
$\Theta$ は係数を複素共役にする反線形写像で、$\braket{\Theta u}{\Theta v}=\braket v u$ を満たす。
$\sigma_x,\sigma_z$ は実行列、$\sigma_y$ は純虚行列なので、直接計算すると $\Theta S_i\Theta^{-1}=-S_i$、$\Theta^2=-I$ となる。
この性質は全体位相の選び方によらない。
磁場も時間反転で逆向きになり、固定された外部磁場は時間反転対称性を破る。

任意の自己共役な $H$ がこの $\Theta$ と交換すれば、$H\ket u=E\ket u$ から $H\Theta\ket u=E\Theta\ket u$ が従う。
二つの状態の直交性も示せる。
反ユニタリ性を $u,\Theta u$ に用いると $\braket{\Theta u}{\Theta^2u}=\braket{\Theta u}u$ だが、左辺は $\Theta^2=-I$ により右辺の負である。
従って内積は零であり、離散準位は少なくとも二重に縮退する。
これが Kramers 縮退である。
整数スピンの場合はこの証明に使う負号がない。

:::example[磁場の実験での判定]
孤立した一つの二準位スピンで時間反転を保つエルミート行列を $H=aI+\vb b\cdot\boldsymbol\sigma$ と書く。
実係数 $\vb b$ の項は反転で符号を変えるため $\vb b=0$ が必要で、二準位は等しい。
外から $B_z$ を加えると $b_z$ が生じ、歳差の周波数として分裂を読み取れる。
:::

### 磁場反転と横方向パルスで起きる再収束

試料の各位置で静磁場がわずかに違い、角周波数のずれが $\delta$ になるとする。横成分を同じ方向に準備しても、時間 $\tau$ の後には異なる位置のスピンが別々の位相を蓄え、集団の横磁化が見かけ上減る。単に外部磁場の設定値の符号を反転しても、同時に反転しない背景ずれは相殺されない。これに対し途中で $x$ 軸に沿う理想的な $\pi$ パルス $R_x=-i\sigma_x$ を加えると $R_x\sigma_z R_x^\dagger=-\sigma_z$ である。従って $Z(\delta\tau)=e^{-i\delta\tau\sigma_z/2}$ として
$$
Z(\delta\tau)R_x Z(\delta\tau)=R_x Z(-\delta\tau)Z(\delta\tau)=R_x .
$$
$\delta$ が試料内で異なっていても、同じスピンで前後の $\delta$ が同じなら最終式に残らない。これがスピンエコーの再収束である。逆に二つの区間で $\delta_1,\delta_2$ が異なれば残差は $R_xZ[(\delta_1-\delta_2)\tau]$ となる。本来の不可逆な緩和と、静的なばらつきによる可逆な位相の散らばりを区別する測定である。

### 回転座標系で横磁場を静止させる

Hamilton 演算子を $H(t)=\hbar[\omega_0 Z+\Omega(X\cos\omega t+Y\sin\omega t)]/2$ とする。これは二つの横成分を最初から回転させる円偏光の磁場である。$W=e^{-i\omega tZ/2}$ とおくと $WXW^\dagger=X\cos\omega t+Y\sin\omega t$。$\ket\psi=W\ket{\psi'}$ を代入して $-i\hbar W^\dagger\dot W=-\hbar\omega Z/2$ を足せば
$$
H'=\frac{\hbar}{2}[(\omega_0-\omega)Z+\Omega X].
$$
$\delta=\omega_0-\omega$、$\Omega_{\rm eff}=\sqrt{\Omega^2+\delta^2}$ とすれば $(H')^2=(\hbar\Omega_{\rm eff}/2)^2 I$。初期上向きから下向きへの振幅の絶対値は $(|\Omega|/\Omega_{\rm eff})|\sin(\Omega_{\rm eff}t/2)|$ である。離調があると、回転座標系での有効軸が傾き、パルス幅をいくら延ばしても完全に下へは移せない。線形に振動する横磁場をこの式に用いるには、逆向きに回転する成分を無視できる条件を別に確かめる。


### 二つの横パルスの位相を変えた Ramsey 測定

二つの短いパルスの間に、共鳴からの角周波数のずれに比例した相対位相が蓄積する。離調は中央の区間だけで加える。両端のパルスは共鳴させ、回転座標系でそれぞれ純粋な横軸の回転とする。

各区間の作用を二次元空間で解く。無次元の行列 $G=\boldsymbol n\cdot\boldsymbol\sigma$ に対して $G^2=I$ なら、指数級数の偶数次と奇数次をまとめて

$$
e^{-i\theta G/2}=I\cos(\theta/2)-iG\sin(\theta/2),\qquad \theta=\Omega\,\delta t.
$$

となる。ここで $H=\hbar\Omega G/2$、一定の角周波数の大きさは $\Omega=6.2832\times10^{5}\,\mathrm{s^{-1}}$ とする。符号反転は軸の反転に含める。初期列は $\begin{pmatrix}1\\0\end{pmatrix}$ とし、各段の直後に読みを取る場合と、全段を連続して加える場合を区別する。

**コヒーレンスの準備。** 角度 $\frac{\pi}{2}$ に要する時間は $2.5\times10^{-6}\,\mathrm s$。作用を左から掛けて次の列を得る。

$$
G_1=\begin{pmatrix}0 & 1\\1 & 0\end{pmatrix},\quad U_1=\begin{pmatrix}\frac{\sqrt{2}}{2} & - \frac{\sqrt{2} i}{2}\\- \frac{\sqrt{2} i}{2} & \frac{\sqrt{2}}{2}\end{pmatrix},\quad \ket{\psi_1}=\begin{pmatrix}\frac{\sqrt{2}}{2}\\- \frac{\sqrt{2} i}{2}\end{pmatrix}.
$$



**自由な位相蓄積。** 角度 $\frac{\pi}{3}$ に要する時間は $1.6667\times10^{-6}\,\mathrm s$。作用を左から掛けて次の列を得る。

$$
G_2=\begin{pmatrix}1 & 0\\0 & -1\end{pmatrix},\quad U_2=\begin{pmatrix}\frac{\sqrt{3}}{2} - \frac{i}{2} & 0\\0 & \frac{\sqrt{3}}{2} + \frac{i}{2}\end{pmatrix},\quad \ket{\psi_2}=\begin{pmatrix}\frac{\sqrt{2} \left(\sqrt{3} - i\right)}{4}\\\frac{\sqrt{2}}{4} - \frac{\sqrt{6} i}{4}\end{pmatrix}.
$$



**逆向きの読み出し。** 角度 $\frac{\pi}{2}$ に要する時間は $2.5\times10^{-6}\,\mathrm s$。作用を左から掛けて次の列を得る。

$$
G_3=\begin{pmatrix}0 & -1\\-1 & 0\end{pmatrix},\quad U_3=\begin{pmatrix}\frac{\sqrt{2}}{2} & \frac{\sqrt{2} i}{2}\\\frac{\sqrt{2} i}{2} & \frac{\sqrt{2}}{2}\end{pmatrix},\quad \ket{\psi_3}=\begin{pmatrix}\frac{\sqrt{3}}{2}\\\frac{1}{2}\end{pmatrix}.
$$

各区間の位相を検査するため、三つの軸の平均と、検出列 $f=\begin{pmatrix}0\\1\end{pmatrix}$ への射影確率を併記する。各行はその段まで加えた別の試料で検査できる。

| 区間の終端 | $\vert \langle f\vert \psi_k\rangle\vert ^2$ | $\langle X\rangle$ | $\langle Y\rangle$ | $\langle Z\rangle$ |
|---|---|---|---|---|
| 1 | $\frac{1}{2}$ | $0$ | $-1$ | $0$ |
| 2 | $\frac{1}{2}$ | $\frac{\sqrt{3}}{2}$ | $- \frac{1}{2}$ | $0$ |
| 3 | $\frac{1}{4}$ | $\frac{\sqrt{3}}{2}$ | $0$ | $\frac{1}{2}$ |

同じ区間をすべて逆順に加えた場合も、作用時間の総和は等しい。積の順だけを入れ替えて

$$
U_{\rm forward}=\begin{pmatrix}\frac{\sqrt{3}}{2} & - \frac{1}{2}\\\frac{1}{2} & \frac{\sqrt{3}}{2}\end{pmatrix},\quad U_{\rm reversed}=\begin{pmatrix}\frac{\sqrt{3}}{2} & \frac{1}{2}\\- \frac{1}{2} & \frac{\sqrt{3}}{2}\end{pmatrix},\quad |\langle f|U_{\rm reversed}|\psi_0\rangle|^2=\frac{1}{4}.
$$

後段と前段が互いに逆回転なので、中央の位相が零なら初期状態へ戻る。この戻りからのずれは自由区間の角度を読む信号である。


### 磁場をスピン反転パルスの長さで校正する

局所の実効横磁場の大きさを、共鳴でスピンが初めて反転する時刻から求める。場は回転座標系で一定、緩和時間はパルスより十分長いとする。

磁気回転比の絶対値を $g=|\gamma|$ と略記すると回転角は $gB_1t$。初めの完全反転には角度 $\pi$ が必要である。

**順方向の予測。** 以下では丸める前に同じ定数を使い、最後に表示する値だけを丸める。

| 入力 | 記号 | 設定値 |
|---|---|---|
| 核スピンの磁気回転比の絶対値 | $g$ | $2.6752\times10^{8}\,\mathrm{s^{-1}\,T^{-1}}$ |
| 実効横磁場 | $B_1$ | $1\times10^{-4}\,\mathrm{T}$ |



$$
t_\pi=\frac{\pi}{B_1 g}\simeq 1.1743393591\times10^{-4}\,\mathrm{s}.
$$

**測定値からの逆算。** 次に実効横磁場だけを未知とする。反転時間を測り、同じ関係を未知量について解く。ここで物理的な範囲にある解を選ぶ。

$$
B_1=\frac{\pi}{g t_\pi}\simeq 1\times10^{-4}\,\mathrm{T}.
$$

この二つの計算は同じ理想関係を逆にたどる。実際の記録では、観測値や別に校正した入力にも誤差がある。各入力の相対標準偏差を下表のように仮定し、誤差は小さく互いに独立とする。対数を微分した係数 $c_i=(u_i/v)(\partial v/\partial u_i)$ は、微小な相対変化の比を符号つきで表す。近い測定値の差が分母にある場合には、その差に比べても小さな変化に限ってこの線形近似を使う。

| 揺れる量 | 対数微分 | 係数の値 | 入力の幅（%） | 推定値への寄与（%） |
|---|---|---|---|---|
| 反転時間 | $-1$ | $-1$ | $0.8$ | $0.8$ |
| 核スピンの磁気回転比の絶対値 | $-1$ | $-1$ | $0.01$ | $0.01$ |

分散を加算してから平方根を取ると、独立な誤差について

$$
\frac{\sigma_v}{|v|}\simeq\sqrt{\sum_i c_i^2\frac{\sigma_i^2}{u_i^2}}\simeq 0.0080006,\qquad \sigma_v\simeq 8.0006\times10^{-7}\,\mathrm{T}.
$$

線形振動場の実験室での最大振幅をそのまま $B_1$ とすると、回転方向を選んだ円成分の振幅との間で係数二が生じる。


### 共鳴から外れた駆動で届く最大反転率

パルス時間を各離調で最適に選ぶと、共鳴から外れたときの反転の限界を有効軸の傾きだけで調べられる。

離調を横角周波数で割った $x=\delta/\Omega$ を用いる。厳密な最大値の Taylor 展開は $1-x^2+O(x^4)$ である。

入力の目盛を $b=6.2832\times10^{5}\,\mathrm{s^{-1}}$、出力の目盛を $d=1$ と定める。入力を $b$ で割った値を $x$ として、最大反転確率の値を $d f(x)$ と書くと

$$
f(x)=\frac{1}{x^{2} + 1},\qquad f_{\rm app}(x)=- \left(x - 1\right) \left(x + 1\right).
$$

となる。ここで $f_{\rm app}$ は比較に使う近似式である。例えば入力が $3.1416\times10^{5}\,\mathrm{s^{-1}}$ のとき、$x=\frac{1}{2}$ なので

$$
\begin{aligned}P&=d f(\frac{1}{2})=d\left(\frac{4}{5}\right)\simeq 0.8,\\ P_{\rm app}&=d f_{\rm app}(\frac{1}{2})=d\left(\frac{3}{4}\right)\simeq 0.75.\end{aligned}
$$

同じ近似の精度をほかの目盛でも調べる。下表の差は近似値から元の式の値を引いているので、符号から過大評価か過小評価かを判定できる。

| $x$ | 入力 $[\mathrm{s^{-1}}]$ | 元の式 $[1]$ | 近似 $[1]$ | 差 $[1]$ |
|---|---|---|---|---|
| $\frac{1}{10}$ | $62832$ | $0.9901$ | $0.99$ | $-9.901\times10^{-5}$ |
| $\frac{3}{10}$ | $1.885\times10^{5}$ | $0.91743$ | $0.91$ | $-0.0074312$ |
| $\frac{1}{2}$ | $3.1416\times10^{5}$ | $0.8$ | $0.75$ | $-0.05$ |
| $1$ | $6.2832\times10^{5}$ | $0.5$ | $0$ | $-0.5$ |
| $\frac{6}{5}$ | $7.5398\times10^{5}$ | $0.40984$ | $-0.44$ | $-0.84984$ |

入力に小さな幅があるときには、表の行を選ぶだけでは出力の幅を決められない。連鎖律により、実際の入力 $u=bx$ に対する傾きは

$$
\frac{\dd P}{\dd u}=\frac{d}{b} f^\prime(x),\qquad f^\prime(x)=- \frac{2 x}{\left(x^{2} + 1\right)^{2}}.
$$

である。入力が基準値から $1\%$ 上にずれた場合を直接計算した差と、接線を使った推定を照合する。

$$
\begin{aligned}\delta P_{\rm lin}&=d f^\prime(\frac{1}{2})\,\frac{1}{200}\simeq -0.0032,\\ \delta P_{\rm exact}&=d\{f(\frac{101}{200})-f(\frac{1}{2})\}\simeq -0.0032031.\end{aligned}
$$

大離調で二次式が負になるのは反転率が負になることを意味しない。級数展開の範囲を越えたことを示す。


### 円錐を一周する磁場で残る幾何学的位相

磁場の大きさを保ちつつ向きだけをゆっくり一周させる。最初に磁場に沿った固有スピンを選べば、最後に同じ方向の測定をしても反転はほぼ検出されない。それでも基準経路と再合成すれば、固有エネルギーの時間積分からは説明できない位相を測れる。

$H=-\hbar\Omega\vb n\cdot\vb\sigma/2$、$\Omega>0$ として、$\vb n$ の極座標を $\theta,\phi$ とする。北極で滑らかな下側の固有ケットを $u=(\cos(\theta/2),e^{i\phi}\sin(\theta/2))^\mathsf T$ と選ぶ。$\psi=e^{i\gamma(t)}e^{-i\int E\,\dd t/\hbar}u$ を Schrödinger 方程式に代入し、左から $u^\dagger$ を掛けると $\dot\gamma=i\langle u|\dot u\rangle$。極角を一定にすれば $\langle u|\dot u\rangle=i\dot\phi\sin^2(\theta/2)$ である。固有空間の外への漏れを無視できる断熱極限で


$$
\gamma=-\int_0^{2\pi}\sin^2(\theta/2)\,\dd\phi=-\pi(1-\cos\theta)=-\frac{\Omega_{\rm solid}}2\pmod{2\pi}.\label{eq-qm-spin-geometric-phase}
$$

ここで $\Omega_{\rm solid}=2\pi(1-\cos\theta)$ はこの向きに囲った立体角で、角周波数 $\Omega$ と区別する。$u\mapsto e^{i\chi}u$ とすると積分は $-\int\dot\chi\,\dd t$ だけ変わるが、初末のケットの位相も同じだけ取り直すので干渉の予測は変わらない。逆回りでは幾何学的位相だけが符号を変え、同じ時間と磁場強度なら力学的位相は等しい。

一定の $\dot\phi$ で回す場合には、回転座標系で $H'=-\hbar[\Omega\sin\theta\,X+(\Omega\cos\theta+\dot\phi)Z]/2$ を得る。有効軸 $\Omega\vb n(0)+\dot\phi\vb e_z$ が非零とし、初期スピン方向との角度を $\beta$ とする。初期スピンを基準軸に取ると、有効軸の垂直成分だけが反転振幅を作るため、その確率の上限は $\sin^2\beta$ である。軸の外積から $\sin\beta=|\dot\phi|\sin\theta/|\Omega\vb n(0)+\dot\phi\vb e_z|$。従って $|\dot\phi|\ll\Omega$ がこの周回での断熱近似の十分条件になる。

```tikz caption="一定の極角を保つ磁場方向。楕円は単位球上の緯線を投影したものである。"
\begin{tikzpicture}[>=Stealth,x=1cm,y=1cm]
\draw[->] (0,-0.2)--(0,3.0) node[above] {$z$};
\draw[dashed] (0,1.8) ellipse (1.65 and 0.48);
\draw[thick,->] (0,0)--(1.62,1.68) node[right] {$\hat n$};
\draw (0,0.72) arc[start angle=90,end angle=48,radius=0.72];
\node at (0.34,0.97) {$\theta$}; \draw[->] (0.35,1.34)--(0.94,1.4);
\end{tikzpicture}
```


### 磁場の静的なばらつきから自由誘導の包絡線を求める

横向きに用意した独立なスピンの集団を観測する。基準角周波数とともに回る座標を使い、個々の局所磁場のずれを零平均の Gaussian 分布とする。観測時間内ではずれは各スピンについて一定とする。



| 設定する量 | 記号 | 計算に用いる値 |
|---|---|---|
| 磁気回転比の絶対値 | $\gamma$ | $2.6752\times10^{8}\,\mathrm{s^{-1}\,T^{-1}}$ |
| 縦磁場の標準偏差 | $\sigma_B$ | $2\times10^{-5}\,\mathrm{T}$ |
| 観測時間 | $t$ | $2.5\times10^{-4}\,\mathrm{s}$ |
| 独立な横緩和の時間 | $T_2$ | $0.008\,\mathrm{s}$ |

**角周波数の標準偏差。** 磁場のずれに実効磁気回転比を掛けて、異なるスピンの位相の進み方に換算する。


$$
\begin{aligned}\sigma_\omega&=\sigma_B \gamma,\\ \frac{\sigma_\omega}{\mathrm{s^{-1}}}&\simeq (2\times10^{-5}) (2.6752\times10^{8})\simeq 5350.4.\end{aligned}
$$

**位相の標準偏差。** 各スピンでずれが固定されているため位相は時間に比例し、分散は時間の二乗に比例する。


$$
\begin{aligned}\sigma_\phi&=t \sigma_\omega,\\ \sigma_\phi&\simeq (2.5\times10^{-4}) 5350.4\simeq 1.3376.\end{aligned}
$$

**ばらつきに由来する横成分。** 正規化した Gaussian の Fourier 変換を時刻 t で評価する。中心の位相はすでに回転座標系で引いた。


$$
\begin{aligned}D&=e^{- \frac{\sigma_\phi^{2}}{2}},\\ D&\simeq e^{- \frac{1.3376^{2}}{2}}\simeq 0.40878.\end{aligned}
$$

**横緩和だけの残存率。** 静的な磁場分布とは独立な指数減衰を仮定する。理想的な中央反転で戻らない成分である。


$$
\begin{aligned}A&=e^{- \frac{t}{T_2}},\\ A&\simeq e^{- \frac{(2.5\times10^{-4})}{0.008}}\simeq 0.96923.\end{aligned}
$$

**自由誘導の振幅。** 独立な二つの減衰要因を順に掛ける。エコーの頂点では静的な要因だけが一に戻る。


$$
\begin{aligned}S&=A D,\\ S&\simeq 0.96923323 \cdot 0.40877644\simeq 0.3962.\end{aligned}
$$

各時間において自由誘導と理想的な中央反転を入れたエコーを比較する。最後の列がエコー頂点に相当する。


| $t\;[\mathrm{s}]$ | $D$ | $S$ | $A$ |
|---|---|---|---|
| $1\times10^{-4}$ | $0.86664$ | $0.85587$ | $0.98758$ |
| $2\times10^{-4}$ | $0.56409$ | $0.55017$ | $0.97531$ |
| $3\times10^{-4}$ | $0.27577$ | $0.26562$ | $0.96319$ |
| $5\times10^{-4}$ | $0.027922$ | $0.02623$ | $0.93941$ |

Gaussian 因子では指数が時間の二乗に、不可逆な横緩和のモデルでは時間の一次に比例する。同じ一時点で横成分が小さいというだけでは二つを区別できない。


## 角運動量の合成と Clebsch–Gordan 係数

二つの角運動量 $\vb J_1,\vb J_2$ は別々の因子に作用するので交差する成分は可換である。
従って $\vb J=\vb J_1+\vb J_2$ も同じ角運動量代数を満たす。
$J_1^2,J_2^2$ は固定し、直積基底 $\ket{m_1,m_2}$ と合成基底 $\ket{J,M}$ の間を

$$
\ket{J,M}=\sum_{m_1+m_2=M} C_{m_1m_2}^{JM}\ket{m_1,m_2}
\label{eq-qm-clebsch-gordan}
$$

で変換する。
係数 $C_{m_1m_2}^{JM}=\braket{m_1,m_2}{J,M}$ を :term[Clebsch–Gordan 係数]{en="Clebsch–Gordan coefficient" yomi="くれぶしゅごるだんけいすう"} という。
許された $M$ に限る理由は $J_z=J_{1z}+J_{2z}$ である。
二基底とも正規直交なので係数行列はユニタリである。

### 許される全角運動量

$j_1\ge j_2$ としてよい。
最大の $M=j_1+j_2$ をもつ積状態はただ一つで $J_+$ に消される。
それを下げて一つの既約多重項を作る。
その直交補空間は $J_\pm^\dagger=J_\mp$ から不変であり、残りの最大 $M$ を次の最高状態としてよい。
具体的に $M\ge0$ の積状態の数 $d_M$ は、$m_2$ の区間 $[-j_2,j_2]\cap[M-j_1,M+j_1]$ の格子点数である。
$j_1-j_2\le M\le j_1+j_2$ では $d_M=j_1+j_2-M+1$、その下では $2j_2+1$ で一定となる。
各多重項は $|M|\le J$ を一度ずつ含むから、$J$ の多重度は $d_J-d_{J+1}$。
従って

$$
J=|j_1-j_2|,|j_1-j_2|+1,\ldots,j_1+j_2
$$

が一度ずつ現れる。
実際、各次元の和は $(2j_1+1)(2j_2+1)$ である。

### 係数の計算法

最高状態を $J_+\ket{J,J}=0$ と規格化から求める。
全体位相は最高状態で最大の $m_1$ をもつ係数を正に選ぶ。
次の行は $J_-=J_{1-}+J_{2-}$ によって求まる。
$M-1=m_1+m_2$ について

$$
\begin{aligned}
\sqrt{(J+M)(J-M+1)}\,C_{m_1m_2}^{J,M-1}&=\sqrt{(j_1-m_1)(j_1+m_1+1)}\,C_{m_1+1,m_2}^{JM}\\
&\quad+\sqrt{(j_2-m_2)(j_2+m_2+1)}\,C_{m_1,m_2+1}^{JM}.
\end{aligned}
$$

範囲外の係数は零である。
これはどの組合せにも使える有限手順で、既知の表を必要としない。

:::example[二つのスピン半分]
$\ket{1,1}=\ket{\uparrow\uparrow}$ に $J_-$ を作用させ $\sqrt2\hbar$ で割ると $\ket{1,0}=(\ket{\uparrow\downarrow}+\ket{\downarrow\uparrow})/\sqrt2$。
もう一度下げると $\ket{1,-1}=\ket{\downarrow\downarrow}$。
直交する $\ket{0,0}=(\ket{\uparrow\downarrow}-\ket{\downarrow\uparrow})/\sqrt2$ は $J_\pm$ で消える。
前三者が三重項、最後が一重項である。
:::

:::example[軌道とスピン半分の一般係数]
$a=\ket{M-1/2,\uparrow}$、$b=\ket{M+1/2,\downarrow}$ と置く。
$A=\sqrt{(\ell+M+1/2)/(2\ell+1)}$、$B=\sqrt{(\ell-M+1/2)/(2\ell+1)}$ とすると、$J=\ell+1/2$ は $Aa+Bb$、$J=\ell-1/2$ は $-Ba+Ab$ である。
前者は最高状態で $A=1$ となり、降下の漸化式に代入すると各段で上式に戻る。
後者は同じ $M$ の直交補で、最高状態の $b$ の符号を正に選んだ。
端で存在しない積状態は省く。
例えば $\ell=1,M=1/2$ では上の多重項の係数は $\sqrt{2/3},\sqrt{1/3}$ である。
:::

### 最高状態だけを直接求める方法

$J=j_1+j_2-k$ とし、$M=J$ で可能な直積を $\ket{j_1-r,j_2-k+r}$ と書く。
三角条件より $0\le k\le2\min(j_1,j_2)$ なので $r=0,\ldots,k$ がすべて許される。
その係数を $b_r$ とすると $J_+$ が消す条件は、$0\le r<k$ で

$$
b_{r+1}\sqrt{(r+1)(2j_1-r)}+b_r\sqrt{(k-r)(2j_2-k+r+1)}=0.
$$

従って交代する符号とすべての大きさの比が一意に決まり、$\sum|b_r|^2=1$、$b_0>0$ で解が完成する。
$M<J$ の係数は本文の降下式で求まる。
各最高状態を既に求めた多重項から Gram–Schmidt で取り出すことと等価だが、この手順なら連立方程式は不要である。

:::example[$1$ と $1$ の合成を上から始める]
$k=0$ なら $\ket{2,2}=\ket{1,1}$（右は二つの磁気量子数）。
$k=1$ なら係数比 $b_1/b_0=-1$ で $\ket{1,1}_{\rm total}=(\ket{1,0}-\ket{0,1})/\sqrt2$。
$k=2$ では $b_1/b_0=b_2/b_1=-1$ だから一重項の三係数は $(1,-1,1)/\sqrt3$。
三つの梯子で $5+3+1=9$ 状態を覆う。
:::

### 一重項の測定で条件を付ける

一重項は全角運動量が零なので、両スピンを同じ回転で動かしても変わらない。
従ってどの軸で第一スピンを測っても上下は半々で、同じ軸で第二スピンを測ると逆向きが確実に出る。
第一を方向 $\vb a$ で上、第二を $\vb b$ で上と記録する同時確率は $(1-\vb a\cdot\vb b)/4$。
初めの結果を捨てると第二の結果はやはり半々である。
条件を付けた確率と、離れた検出器だけで読める周辺確率を区別する例である。

:::example[四組の測定軸で一重項相関を比べる]
二つの検出器で $\pm1$ を返すスピン測定を行い、第一の軸を $\vb a,\vb a'$、第二の軸を $\vb b,\vb b'$ から選ぶ。四つの値 $A,A',B,B'\in\{-1,1\}$ が共通の確率変数として定まり、遠い検出器の軸の選択に依存しない模型では、$A(B+B')+A'(B-B')$ は各試行で $\pm2$ になる。二つの括弧のうち一方が零、他方が $\pm2$ だからである。軸の選択が変数の分布を変えないなら、四組の相関を足し引きした値の絶対値は高々二である。これが :term[CHSH 不等式]{en="CHSH inequality" yomi="しーえいちえすえいちふとうしき"} である。

一重項では相関が $-\vb a\cdot\vb b$ になる。$\vb a=\vb e_z,\ \vb a'=\vb e_x,\ \vb b=(\vb e_z+\vb e_x)/\sqrt2,\ \vb b'=(\vb e_z-\vb e_x)/\sqrt2$ を選ぶと、同じ組合せは $-2\sqrt2$ になる。四組は同じ準備をした別々の試料で測る。各検出器だけの確率は軸の選び方によらず上下半々であり、この違いが現れるのは二つの記録を照合した相関である。
:::

### 行列要素を角度部分と内部構造に分ける

先の合成係数を、[角運動量の章](/physics/quantum-mechanics/angular-momentum/) で定義したベクトル演算子の球成分 $V_q$ に適用する。
同じ $j$ が複数回現れる場合は、残りの量子数を $\alpha$ として状態を $\ket{\alpha,j,m}$ と区別する。
写像 $F$ を直積空間上で

$$
F(\ket{j,m}\otimes\ket{1,q})=V_q\ket{\alpha,j,m}
$$

と定める。
$J_iV_q\ket{\alpha,j,m}=V_qJ_i\ket{\alpha,j,m}+[J_i,V_q]\ket{\alpha,j,m}$ の右辺で、第一項は最初の因子、第二項はスピン一の行列として第二の因子に作用する。
従って $F$ は全角運動量の作用と交換する写像である。

合成基底の最高状態 $\ket{J,J}$ を $F$ で写すと、結果はやはり磁気量子数 $J$ をもち $J_+$ で消える。
指定した終多重項 $\alpha',j'$ への射影は $j'=J$ の場合にだけ残る。
いったん最高状態で比例定数が定まれば、$J_-$ と $F$ が交換するので、その多重項の下の全状態で同じ定数になる。
直積基底へ戻して

$$
\mel{\alpha',j',m'}{V_q}{\alpha,j,m}=C_{m q}^{j'm'}\,v_{\alpha'j',\alpha j}.
\label{eq-qm-wigner-eckart}
$$

ここでは右辺の $C$ は $j$ と $1$ の合成係数、$v$ は $m,m',q$ に依存しない定数である。
同じ証明は、第二因子を階数 $k$ の成分 $T_q^{(k)}$ に置き換えても働く。ここで階数 $k$ の :term[球面テンソル演算子]{en="spherical tensor operator" yomi="きゅうめんてんそるえんざんし"} は $[J_z,T_q^{(k)}]=\hbar q T_q^{(k)}$、$[J_\pm,T_q^{(k)}]=\hbar\sqrt{(k\mp q)(k\pm q+1)}\,T_{q\pm1}^{(k)}$ を満たす $2k+1$ 個の成分である。この交換関係が、第二因子の作用が角運動量 $k$ になることを保証する。従って行列要素には $j$ と $k$ の合成係数が現れ、$m'=m+q$ と $|j-k|\le j'\le j+k$ が必要である。これが :term[Wigner–Eckart の定理]{en="Wigner–Eckart theorem" yomi="うぃぐなーえっかるとのていり"} である。
規約によっては $v$ に $\sqrt{2j'+1}$ を掛けたものを換算行列要素と呼ぶ。
上の定義に従えば行列要素の比にはその違いが入らない。

:::example[強度の比を角度だけで求める]
下準位 $j=1/2,m=1/2$ から上準位 $j'=3/2$ へ、$q=+1$ なら伸び切った積状態なので係数の大きさは一。
$q=0$ なら、その上状態に $J_-$ を一度作用させた式から係数は $\sqrt{2/3}$。
同じ球成分の大きさの外力で比べると、遷移強度の比は $1:(2/3)=3:2$ になる。
:::

### 三つの角運動量を合成するときの多重度

三つの角運動量では、最初に二つを合成して得た $J_{12}$ ごとに、残りとの合成を行う。$J_{12}^2$ は全 $J^2,J_z$ と可換なので、中間の量子数を直交基底の区別に使える。例えば三つのスピン半分では、初めの二個を $0$ にしても $1$ にしても、最後に全 $J=1/2$ を作れる。

:::example[三スピンの二つの二重項]
$M=1/2$ で $u=\ket{\uparrow\uparrow\downarrow}$、$v=\ket{\uparrow\downarrow\uparrow}$、$w=\ket{\downarrow\uparrow\uparrow}$ と書く。全対称な $(u+v+w)/\sqrt3$ は $J=3/2$ の最高状態を一回下げて得る。$(v-w)/\sqrt2$ は初めの二スピンが一重項である。残る $(2u-v-w)/\sqrt6$ は初めの二スピンが対称で、全対称な列に直交する。後二者はどちらも $J=1/2$。$J^2,J_z$ の二つだけを読んでも、この二重項を見分けられない。
:::


### 軌道二とスピン半分を固定した磁気量子数で結ぶ

軌道成分が大きくても、スピン半分を結ぶ一つの $M$ ブロックの次元は高々二である。状態の符号と強度を一つずつ求める。

固定した全磁気量子数 $M=\frac{1}{2}$ の中だけを扱う。$j_1=2,j_2=\frac{1}{2}$ で、積基底の順を $\ket{1,- \frac{1}{2}}$、$\ket{0,\frac{1}{2}}$ と決める。他の $M$ への行列要素は、上げ下げが第一・第二の添字で必ず逆向きに起こるので零である。

$$
D=\frac{\boldsymbol J_1\cdot\boldsymbol J_2}{\hbar^2}=\frac{J_{1z}J_{2z}}{\hbar^2}+\frac{J_{1+}J_{2-}+J_{1-}J_{2+}}{2\hbar^2}=\begin{pmatrix}- \frac{1}{2} & \frac{\sqrt{6}}{2}\\\frac{\sqrt{6}}{2} & 0\end{pmatrix}.
$$

対角成分は $m_1m_2$。例えば第一の添字を下げる成分は $\sqrt{(j_1+m_1)(j_1-m_1+1)(j_2-m_2)(j_2+m_2+1)}/2$ である。逆向きはその共役であり、状態の位相をこの順に取るとどちらも実である。$J^2=J_1^2+J_2^2+2\hbar^2D$ を使えば $D$ の根を全角運動量のラベルと対応させられる。

$$
\det(zI-D)=\frac{\left(z - 1\right) \left(2 z + 3\right)}{2},\qquad d_J=\tfrac12[J(J+1)-j_1(j_1+1)-j_2(j_2+1)].
$$

根ごとに $(D-d_J)c=0$ を解く。全体位相は標準の Clebsch–Gordan 係数に合わせた。独立な行を選んで成分比を求め、最後にノルム一へ直した列と、その列を零に写す行列は次の通りである。

| $J$ | $d_J$ | 零化する行列 $D-d_J I$ | 規格化した列 $c_J$ |
|---|---|---|---|
| $\frac{5}{2}$ | $1$ | $\begin{pmatrix}- \frac{3}{2} & \frac{\sqrt{6}}{2}\\\frac{\sqrt{6}}{2} & -1\end{pmatrix}$ | $\begin{pmatrix}\frac{\sqrt{10}}{5}\\\frac{\sqrt{15}}{5}\end{pmatrix}$ |
| $\frac{3}{2}$ | $- \frac{3}{2}$ | $\begin{pmatrix}1 & \frac{\sqrt{6}}{2}\\\frac{\sqrt{6}}{2} & \frac{3}{2}\end{pmatrix}$ | $\begin{pmatrix}\frac{\sqrt{15}}{5}\\- \frac{\sqrt{10}}{5}\end{pmatrix}$ |

従って積基底から見た結合基底の行列は $U=\begin{pmatrix}\frac{\sqrt{10}}{5} & \frac{\sqrt{15}}{5}\\\frac{\sqrt{15}}{5} & - \frac{\sqrt{10}}{5}\end{pmatrix}$。各列が直交して $U^\dagger U=I$ である。特に表の最初の積状態を準備したときは、最初の行を読み、絶対値の二乗を取ればよい。ここでは分裂の有効係数を $1.0\,\mathrm{meV}$ に選んだ。

| 測定された $J$ | 確率 | 相互作用 $bD$ の値 |
|---|---|---|
| $\frac{5}{2}$ | $\frac{2}{5}$ | $1.6022\times10^{-22}\,\mathrm J$ |
| $\frac{3}{2}$ | $\frac{3}{5}$ | $-2.4033\times10^{-22}\,\mathrm J$ |

ここで相互作用の目盛は $b=1.6022\times10^{-22}\,\mathrm J$ とした。同じ $J$ の内部で $M$ を変えても相互作用のエネルギーは一致する。外部磁場で方向を選べばこの縮退を分離できる。


### 一方のスピンだけで相対位相を測れるか

二粒子を $\ket{\uparrow\downarrow}$ に準備して交換結合を一定時間働かせる。一重項と三重項に異なる位相が付く。両項の相対位相が $\pi/2$ になったときの状態は局所位相の規約を選べば $(\ket{\uparrow\downarrow}+i\ket{\downarrow\uparrow})/\sqrt2$ と書ける。

二つの自由度を順に $A,B$ と呼ぶ。この例で用いる係数は $\ket\Psi=\sum_{ab}C_{ab}\ket a\ket b$ の形で一枚の行列に並べられる。第二の添字に沿って列を読むと、$B$ に対応する読みごとに $A$ の振幅が得られる。

$$
C=\begin{pmatrix}0 & \frac{\sqrt{2}}{2}\\\frac{\sqrt{2} i}{2} & 0\end{pmatrix},\qquad \sum_{ab}|C_{ab}|^2=\tr(CC^\dagger)=1.
$$

**片側だけを読む。** $A$ の量 $M$ の期待値を成分で書くと $\sum_{a,a^{\prime},b} C_{ab}^* M_{a a^{\prime}} C_{a^{\prime}b}$ である（$a'$ は第二の和の添字）。$b$ を先に足すと、$\rho_A=CC^\dagger$ として $\tr(\rho_A M)$ が得られる。逆に $B$ を調べる場合は共役を取る側が入れ替わる。

$$
\rho_A=\begin{pmatrix}\frac{1}{2} & 0\\0 & \frac{1}{2}\end{pmatrix},\qquad \rho_B=C^{\mathsf T}C^*=\begin{pmatrix}\frac{1}{2} & 0\\0 & \frac{1}{2}\end{pmatrix}.
$$

行列の対角成分はこの基底での局所確率である。両者の非零固有値が一致する理由は $CC^\dagger u=\lambda u$ に左から $C^\dagger$ を掛ければ分かる。$\lambda>0$ なら $C^\dagger u\ne0$ で、$C^\dagger C$ の同じ固有値に属する。逆向きも同様であり、転置で特性多項式は変わらない。今回の純粋性は

$$
\tr\rho_A^2=\frac{1}{2}=\tr\rho_B^2.
$$

である。**条件を付けた測定。** $B$ をこの基底で読み、値を選別できる場合は、同じ実験でも $A$ の統計をさらに分けられる。各列 $c_b$ に対して $p_b=c_b^\dagger c_b$、$p_b>0$ で $\ket{v_b}=c_b/\sqrt{p_b}$ となる。第一のスピンで $z$ 成分を読む。 $b_0=5.2729\times10^{-35}\,\mathrm{J\,s}$ を実際の目盛、$M=\begin{pmatrix}1 & 0\\0 & -1\end{pmatrix}$ を無次元行列とする。

| $B$ の読み（列番号は零から） | 割合 | 条件付きの $A$ の列 | 条件付き平均 |
|---|---|---|---|
| 0 | $\frac{1}{2}$ | $\begin{pmatrix}0\\i\end{pmatrix}$ | $-5.2729\times10^{-35}\,\mathrm{J\,s}$ |
| 1 | $\frac{1}{2}$ | $\begin{pmatrix}1\\0\end{pmatrix}$ | $5.2729\times10^{-35}\,\mathrm{J\,s}$ |

すべての列で平均すると $\sum_b p_b\langle v_b|M|v_b\rangle=0$ であり、縮約した行列から計算した $\tr(\rho_A M)$ と一致する。さらに $B$ の測定基底をユニタリ行列 $U$ で変えると、列を作る行列は $CU^*$ になる。$U^*U^{\mathsf T}=I$ より局所の $\rho_A$ は変わらない。条件付きの列そのものは基底に依存する。

相対位相が違う Bell 状態の局所密度が同じになることがある。その位相を復元するには、第二のスピンも測る。例えば掲載した位相では $\langle X\otimes Y\rangle=-1$、$i$ を $-i$ に替えると $+1$ となる。


### 対称な三スピンの一反転状態を回転させる

同じ向きの回転を三つのスピン半分に与える。最大全スピンの部分空間は保たれるので、$j=3/2$、$m=1/2$ の対称な一反転状態を一つのスピンとして追う。

回転軸を $y$ とし、正方向の能動回転を $U=\exp(-i\beta J_y/\hbar)$ で定める。角度は $c=\cos(\beta/2)=\frac{\sqrt{2}}{2}$、$s=\sin(\beta/2)=\frac{\sqrt{2}}{2}$ とする。まずスピン半分の二列に作用させると $\ket+\mapsto c\ket+ +s\ket-$、$\ket-\mapsto-s\ket++c\ket-$。後者のマイナスは二列の直交を保つために必要である。任意の $j$ については [角運動量の章](/physics/quantum-mechanics/angular-momentum/) と同じ多項式表示を使う。$p=j+m,q=j-m$ として $u^pv^q/\sqrt{p!q!}$ を $\ket{j,m}$ に対応させる。変数を掛ける操作と微分する操作は、それぞれ二つの独立な調和振動子の生成・消滅演算子と同じ係数をもつ。従って $J_+=\hbar a_u^\dagger a_v$、$J_-=\hbar a_v^\dagger a_u$、$J_z=\hbar(a_u^\dagger a_u-a_v^\dagger a_v)/2$ とも表せる。全励起数が $2j$ の部分空間で角運動量一組を実現するこの見方を、Schwinger の振動子表示という。

$$
U\frac{u^pv^q}{\sqrt{p!q!}}=\frac{(cu+sv)^p(-su+cv)^q}{\sqrt{p!q!}}.
$$

全次数が固定されるので、$v=1$ と置いても各単項式の二つの次数は一意に戻せる。今回の $j=\frac{3}{2},m=\frac{1}{2}$ を入れ、展開前後の係数を示す。次式の多項式は初期規格化係数をまだ掛けていない。

$$
(cu+s)^p(-su+c)^q=- \frac{\sqrt{2} u^{3}}{4} - \frac{\sqrt{2} u^{2}}{4} + \frac{\sqrt{2} u}{4} + \frac{\sqrt{2}}{4}.
$$



| 最後の $m$ | 単項式の係数 | 規格化した振幅 | 確率 | 確率の概数 |
|---|---|---|---|---|
| $\frac{3}{2}$ | $- \frac{\sqrt{2}}{4}$ | $- \frac{\sqrt{6}}{4}$ | $\frac{3}{8}$ | $0.375$ |
| $\frac{1}{2}$ | $- \frac{\sqrt{2}}{4}$ | $- \frac{\sqrt{2}}{4}$ | $\frac{1}{8}$ | $0.125$ |
| $- \frac{1}{2}$ | $\frac{\sqrt{2}}{4}$ | $\frac{\sqrt{2}}{4}$ | $\frac{1}{8}$ | $0.125$ |
| $- \frac{3}{2}$ | $\frac{\sqrt{2}}{4}$ | $\frac{\sqrt{6}}{4}$ | $\frac{3}{8}$ | $0.375$ |

表を直接足す以外にも、軸を逆向きに回した演算子で一次、二次モーメントを確かめられる。$U^\dagger J_z U=J_z\cos\beta-J_x\sin\beta$ だから初状態で $\langle J_x\rangle=0$ を使い、共分散も零であることを梯子から読む。

$$
\frac{\langle J_z\rangle_f}{\hbar}=m_i\cos\beta=0,\quad \frac{(\Delta J_z)^2_f}{\hbar^2}=\frac{j(j+1)-m_i^2}{2}\sin^2\beta=\frac{7}{4}.
$$

角速度が $\Omega=6.2832\times10^{5}\,\mathrm{s^{-1}}$ の一様な場なら必要な時間は $2\arccos(c)/\Omega=2.5\times10^{-6}\,\mathrm s$。同じ大きさの逆回転を直後に加えれば、途中で測らない限り逆行列により初状態に戻る。ここで途中の $J_z$ を測って結果を捨てると、最後の復帰確率は各経路の確率を足した $\sum_{m_f}|d_{m_fm_i}|^4=\frac{5}{16}$ に下がる。粒子ごとに独立な向きの場を加えると全スピンの大きさも変わり、この四次元部分空間だけでは発展を記述できない。


## スピン軌道相互作用

電子から見ると核の電場は磁場の成分を生じ、スピンの磁気エネルギーが電子の運動方向に依存する。
低速の電子で $g_s=2$、核固定、中心ポテンシャル $V$ を仮定する。
電子静止系の誘導磁場 $-\vb v\times\vb E/c^2$、$e\vb E=\nabla V$ と $\boldsymbol\mu=-e\vb S/m_e$ を使うと、磁気相互作用は $(m_e^2c^2r)^{-1}V'(r)\vb L\cdot\vb S$ となる。
ところが電子の瞬間静止系は運動に伴って変わる。磁気トルクをいったん除き、静止系でトルクなしに運ばれるスピンを考える。計量の符号を $(-,+,+,+)$ とし、スピン四元ベクトル $s^\mu$ と四元速度 $u^\mu$ は $u\cdot s=0$ を満たす。瞬間静止系で空間成分の固有時微分を零にするには、$\dd s^\mu/\dd\tau$ が $u^\mu$ に比例しなければならない。$u\cdot s=0$ の微分と $u\cdot u=-c^2$ から、その係数は $a\cdot s/c^2$ となる。ここで $a^\mu=\dd u^\mu/\dd\tau$ である。

静止系の空間スピンを $\vb S$ とすれば、Lorentz 変換で実験室系の空間成分は低速で $\vb s=\vb S+\vb v(\vb v\cdot\vb S)/(2c^2)$ になる。上の輸送式の空間成分 $\dot{\vb s}=\vb v(\dot{\vb v}\cdot\vb S)/c^2$ から補正項の時間微分を引くと

$$
\dot{\vb S}=\frac{\dot{\vb v}\times\vb v}{2c^2}\times\vb S
$$

が得られる。これが最低次の :term[Thomas 歳差]{en="Thomas precession" yomi="とーますさいさ"} である。$\dot{\vb v}=-\nabla V/m_e$ を代入すると歳差の角速度は $-V'(r)\vb L/(2m_e^2c^2r)$。歳差を生成する項はこの角速度と $\vb S$ の内積だから、先の磁気相互作用の半分を打ち消す。従って有効 Hamilton 演算子は

$$
\begin{aligned}
H_{\rm so}&=\frac{1}{2m_e^2c^2r}\frac{\dd V}{\dd r}\,\vb L\cdot\vb S, \\
\vb L\cdot\vb S&=\frac{J^2-L^2-S^2}{2}.\label{eq-qm-spin-orbit}
\end{aligned}
$$

角因子は $j=\ell+1/2$ で $\hbar^2\ell/2$、$j=\ell-1/2$ で $-\hbar^2(\ell+1)/2$。
$\ell=0$ では前者だけがあり角因子は零である。固定した動径状態に射影すると、中心力の詳細は動径の係数の平均に入る。[時間に依存しない摂動論](/physics/quantum-mechanics/stationary-perturbation/) では、同じ相互作用を Dirac 方程式から導く。

:::example[$p$ 準位の分裂]
$\ell=1$ で $H_{\rm so}=\xi\,\vb L\cdot\vb S/\hbar^2$ とする。
$j=3/2$ の四状態は $+\xi/2$、$j=1/2$ の二状態は $-\xi$。
縮退度を掛けて足すと零となり、重心は変わらない。
:::

### 相互作用が保存する角運動量

$F(r)$ は軌道角運動量の三成分と可換であるから、$H_{\rm so}=F(r)\vb L\cdot\vb S$ は $L^2,S^2,J^2,J_z$ と可換となる。
ただし一般には $L_z$ と $S_z$ に個別の保存則はない。
実際、$L_-S_+$ は $m_\ell$ を一つ下げ $m_s$ を一つ上げるため $m_j$ は保たれる。
外部磁場も加えたときの基底の選び方は [Zeeman 効果](/physics/quantum-mechanics/stationary-perturbation/#zeeman-効果) で扱う。

:::example[$p$ 軌道での対角化を二通りで行う]
$m_j=1/2$ の直積 $\ket{0,\uparrow},\ket{1,\downarrow}$ の順に並べる。
$L_zS_z$ の対角は $0,-\hbar^2/2$、$L_+S_-/2$ による非対角は $\hbar^2/\sqrt2$ である。
$\vb L\cdot\vb S/\hbar^2$ の固有方程式は $\lambda^2+\lambda/2-1/2=0$、固有値は $1/2,-1$ となり $j=3/2,1/2$ の値と一致する。
:::

### 微細な分裂の重心と間隔則

ある一つの動径配置で有効係数 $\xi$ を一定とできる範囲で $H_{\rm so}=\xi\,\vb L\cdot\vb S/\hbar^2$ と書く。$E_J=\xi[J(J+1)-L(L+1)-S(S+1)]/2$ だから、固定した $L,S$ に属する隣接準位の差は $E_J-E_{J-1}=\xi J$ である。これを Landé の間隔則という。スペクトルの線強度は状態の準備や遷移行列要素にも依存するので、この式だけで各線の高さは決まらない。

多重項全体での重心は分裂だけでは動かない。直積基底で $\tr(L_i\otimes S_i)=(\tr L_i)(\tr S_i)=0$ である。$z$ 成分は正負対で跡が零、$x,y$ 成分も梯子を一段動かすので対角成分が零だからである。一方、結合基底で同じ跡を取ると $\sum_J(2J+1)E_J=0$。従って線の位置を比較する際の基準は、各準位の個数 $2J+1$ で重み付けした平均である。


### 三つの J に分かれる多重項の中心を調べる

軌道 $L=1$ と全スピン $S=1$ をもつ配置では、角運動量の大きさの組を固定しても三つの分裂が起こる。$M=0$ で三本とも現れるので同じ行列から検証できる。

固定した全磁気量子数 $M=0$ の中だけを扱う。$j_1=1,j_2=1$ で、積基底の順を $\ket{1,-1}$、$\ket{0,0}$、$\ket{-1,1}$ と決める。他の $M$ への行列要素は、上げ下げが第一・第二の添字で必ず逆向きに起こるので零である。

$$
D=\frac{\boldsymbol J_1\cdot\boldsymbol J_2}{\hbar^2}=\frac{J_{1z}J_{2z}}{\hbar^2}+\frac{J_{1+}J_{2-}+J_{1-}J_{2+}}{2\hbar^2}=\begin{pmatrix}-1 & 1 & 0\\1 & 0 & 1\\0 & 1 & -1\end{pmatrix}.
$$

対角成分は $m_1m_2$。例えば第一の添字を下げる成分は $\sqrt{(j_1+m_1)(j_1-m_1+1)(j_2-m_2)(j_2+m_2+1)}/2$ である。逆向きはその共役であり、状態の位相をこの順に取るとどちらも実である。$J^2=J_1^2+J_2^2+2\hbar^2D$ を使えば $D$ の根を全角運動量のラベルと対応させられる。

$$
\det(zI-D)=\left(z - 1\right) \left(z + 1\right) \left(z + 2\right),\qquad d_J=\tfrac12[J(J+1)-j_1(j_1+1)-j_2(j_2+1)].
$$

根ごとに $(D-d_J)c=0$ を解く。全体位相は標準の Clebsch–Gordan 係数に合わせた。独立な行を選んで成分比を求め、最後にノルム一へ直した列と、その列を零に写す行列は次の通りである。

| $J$ | $d_J$ | 零化する行列 $D-d_J I$ | 規格化した列 $c_J$ |
|---|---|---|---|
| $2$ | $1$ | $\begin{pmatrix}-2 & 1 & 0\\1 & -1 & 1\\0 & 1 & -2\end{pmatrix}$ | $\begin{pmatrix}\frac{\sqrt{6}}{6}\\\frac{\sqrt{6}}{3}\\\frac{\sqrt{6}}{6}\end{pmatrix}$ |
| $1$ | $-1$ | $\begin{pmatrix}0 & 1 & 0\\1 & 1 & 1\\0 & 1 & 0\end{pmatrix}$ | $\begin{pmatrix}\frac{\sqrt{2}}{2}\\0\\- \frac{\sqrt{2}}{2}\end{pmatrix}$ |
| $0$ | $-2$ | $\begin{pmatrix}1 & 1 & 0\\1 & 2 & 1\\0 & 1 & 1\end{pmatrix}$ | $\begin{pmatrix}\frac{\sqrt{3}}{3}\\- \frac{\sqrt{3}}{3}\\\frac{\sqrt{3}}{3}\end{pmatrix}$ |

従って積基底から見た結合基底の行列は $U=\begin{pmatrix}\frac{\sqrt{6}}{6} & \frac{\sqrt{2}}{2} & \frac{\sqrt{3}}{3}\\\frac{\sqrt{6}}{3} & 0 & - \frac{\sqrt{3}}{3}\\\frac{\sqrt{6}}{6} & - \frac{\sqrt{2}}{2} & \frac{\sqrt{3}}{3}\end{pmatrix}$。各列が直交して $U^\dagger U=I$ である。特に表の最初の積状態を準備したときは、最初の行を読み、絶対値の二乗を取ればよい。軌道とスピンの有効結合を正の $2\,\mathrm{meV}$ とする。

| 測定された $J$ | 確率 | 相互作用 $bD$ の値 |
|---|---|---|
| $2$ | $\frac{1}{6}$ | $3.2044\times10^{-22}\,\mathrm J$ |
| $1$ | $\frac{1}{2}$ | $-3.2044\times10^{-22}\,\mathrm J$ |
| $0$ | $\frac{1}{3}$ | $-6.4087\times10^{-22}\,\mathrm J$ |

ここで相互作用の目盛は $b=3.2044\times10^{-22}\,\mathrm J$ とした。それぞれ $J=0,1,2$ に対する重み $1,3,5$ を用いて三つのエネルギーを足すと零であり、隣接ギャップの比は一対二になる。


### 微細な二本の線からスピン軌道定数を抽出する

同じ動径配置にある $d$ 電子の二つの $j$ の差を測る。外場は十分弱く、二本は同じ下準位に落ちるとする。

スピン半分に対して $j=\ell+1/2,\ell-1/2$ の間隔は $(2\ell+1)|\xi|/2$。$\ell=2$ を代入して周波数差へ直す。

**順方向の予測。** 以下では丸める前に同じ定数を使い、最後に表示する値だけを丸める。

| 入力 | 記号 | 設定値 |
|---|---|---|
| 正の有効結合定数 | $\xi$ | $3.2044\times10^{-23}\,\mathrm{J}$ |



$$
\delta\nu=\frac{5 \xi}{2 h}\simeq 1.208994621\times10^{11}\,\mathrm{Hz}.
$$

**測定値からの逆算。** 次に正の有効結合定数だけを未知とする。線の周波数差を測り、同じ関係を未知量について解く。ここで物理的な範囲にある解を選ぶ。

$$
\xi=\frac{2 h \delta\nu}{5}\simeq 3.2044\times10^{-23}\,\mathrm{J}.
$$

この二つの計算は同じ理想関係を逆にたどる。実際の記録では、観測値や別に校正した入力にも誤差がある。各入力の相対標準偏差を下表のように仮定し、誤差は小さく互いに独立とする。対数を微分した係数 $c_i=(u_i/v)(\partial v/\partial u_i)$ は、微小な相対変化の比を符号つきで表す。近い測定値の差が分母にある場合には、その差に比べても小さな変化に限ってこの線形近似を使う。

| 揺れる量 | 対数微分 | 係数の値 | 入力の幅（%） | 推定値への寄与（%） |
|---|---|---|---|---|
| 線の周波数差 | $1$ | $1$ | $0.5$ | $0.5$ |

分散を加算してから平方根を取ると、独立な誤差について

$$
\frac{\sigma_v}{|v|}\simeq\sqrt{\sum_i c_i^2\frac{\sigma_i^2}{u_i^2}}\simeq 0.005,\qquad \sigma_v\simeq 1.6022\times10^{-25}\,\mathrm{J}.
$$

測った差の大きさだけでは相互作用の符号は一意でない。どちらの $j$ が上だったかという割り当てが符号を決める。


### d 軌道でスピン軌道結合から磁場優勢へ移る過程

一定の動径準位について、全磁気量子数だけを保存したまま磁場を強くする。同じ準位の呼び名がどの範囲で有効かを積状態の占有から判断する。

一定の $\ell=2$ と $M=m_\ell+m_s=- \frac{1}{2}$ を選ぶ。磁場に沿う磁気量子数は保存するので、積基底 $\ket{-1,\uparrow},\ket{0,\downarrow}$ 内で対角化できる。$g_s=2$ とし、スピン軌道係数を $\xi=1.6022\times10^{-23}\,\mathrm J$ に選ぶ。

$$
H/\xi=D+xZ,\qquad D=\begin{pmatrix}- \frac{1}{2} & \frac{\sqrt{6}}{2}\\\frac{\sqrt{6}}{2} & 0\end{pmatrix},\quad Z=\begin{pmatrix}0 & 0\\0 & -1\end{pmatrix},\quad x=\mu_{\rm B}B/\xi.
$$

対角要素は $m_\ell m_s+x(m_\ell+2m_s)$。非対角要素は $L_+S_-/2$ と $L_-S_+/2$ の間で結ぶ係数である。$A=(D+xZ)_{11}, C=(D+xZ)_{22},v=D_{12}$ と書けば、二次方程式から

$$
\epsilon_\pm=\frac{A+C}{2}\pm\frac12\sqrt{(A-C)^2+4v^2},\quad P_\pm=\frac{D+xZ-\epsilon_\mp I}{\epsilon_\pm-\epsilon_\mp}.
$$

となる。上の状態で最初の積ケットを見つける確率は射影の $(1,1)$ 成分である。磁場を正の方向に変え、エネルギーの順と列の成分を同時に追う。

| 磁場 $[\mathrm T]$ | $x$ | 下枝 $[\mathrm J]$ | 上枝 $[\mathrm J]$ | 上枝での第一成分の確率 |
|---|---|---|---|---|
| $0$ | $0$ | $-2.4033\times10^{-23}$ | $1.6022\times10^{-23}$ | $0.4$ |
| $0.05$ | $0.028942$ | $-2.4219\times10^{-23}$ | $1.5745\times10^{-23}$ | $0.40558$ |
| $0.5$ | $0.28942$ | $-2.6019\times10^{-23}$ | $1.3371\times10^{-23}$ | $0.45717$ |
| $1$ | $0.57884$ | $-2.8275\times10^{-23}$ | $1.099\times10^{-23}$ | $0.51608$ |
| $5$ | $2.8942$ | $-5.463\times10^{-23}$ | $2.4859\times10^{-25}$ | $0.84949$ |
| $10$ | $5.7884$ | $-9.7064\times10^{-23}$ | $-3.6871\times10^{-24}$ | $0.9537$ |

弱磁場では $D$ の固有列で $Z$ の対角平均を読む。強磁場では $Z$ の固有列、すなわち積基底で $D$ の対角平均を読む。この二つの操作は同じ行列を異なる比 $x$ で近似したものである。異なる $M$ 同士を比較するときには、それぞれの行列の平均エネルギーも保つ必要がある。交換される二列の成分を直接読むことで、弱磁場の $J$ と強磁場の積状態という二つの記述を連続的に結べる。


## 演習問題

::::exercise{level=basic}
$\ket{+y}$ を $z$ 基底で書き、$\sigma_y$ を作用させよ。

:::solution
列 $(1,i)^\top/\sqrt2$ に $\sigma_y$ を掛けると $(-i\,i,i)^\top/\sqrt2=(1,i)^\top/\sqrt2$ となる。
固有値は一である。
:::
::::

::::exercise{level=basic}
$z$ 上向きの状態で $\Delta S_x$ を求めよ。

:::solution
$\expval{S_x}=0$、$S_x^2=\hbar^2 I/4$ より標準偏差は $\hbar/2$ である。
:::
::::

::::exercise{level=basic}
$j_1=2,j_2=1$ の直積を分解して次元を照合せよ。

:::solution
$J=3,2,1$ の次元は $7,5,3$。
和は $15=(2\cdot2+1)(2\cdot1+1)$ となる。
:::
::::

::::exercise{level=basic}
$\ket{\uparrow\downarrow}$ で $J^2$ を測る結果と確率を求めよ。

:::solution
三重項と一重項の $M=0$ の等重ね合わせなので、$2\hbar^2$ と零を確率 $1/2$ ずつ得る。
:::
::::

::::exercise{level=advanced}
一重項で任意方向 $\vb a,\vb b$ のスピン相関 $\expval{(\vb a\cdot\boldsymbol\sigma_1)(\vb b\cdot\boldsymbol\sigma_2)}$ を求めよ。

:::solution
一重項は全回転で変わらず、$\expval{\sigma_{1i}\sigma_{2j}}$ は $\delta_{ij}$ の定数倍となる。
$z$ 測定は必ず逆向きなので定数は $-1$。
従って相関は $-\vb a\cdot\vb b$ である。
:::
::::

::::exercise{level=advanced}
$j_1=j_2=1$ の $J=0$ の状態を構成せよ。

:::solution
$M=0$ の状態を $a\ket{1,-1}+b\ket{0,0}+c\ket{-1,1}$ と直積で書く。
$J_+$ による二成分が零となる条件は $a+b=b+c=0$。
規格化と位相の規約から $(a,b,c)=(1,-1,1)/\sqrt3$ である。
:::
::::

::::exercise{level=advanced}
$\ell=2$ にスピン半分を加えたときのスピン軌道角因子を求めよ。

:::solution
$j=5/2,3/2$ で $\vb L\cdot\vb S/\hbar^2=[j(j+1)-6-3/4]/2$。
順に $1,-3/2$ となる。
次元は $6,4$ で加重和は零である。
:::
::::

::::exercise{level=advanced}
$\rho=(I+s\sigma_z)/2$ を $y$ 軸周りに $\pi/2$ 回転して $S_x$ 上向きの確率を求めよ。

:::solution
能動回転で $+\vb e_z$ は $+\vb e_x$ へ移る。
従って $\rho'=(I+s\sigma_x)/2$ であり、$\tr(\rho'P_{+x})=(1+s)/2$ となる。
:::
::::

::::exercise{level=advanced}
$\ell=1$ の軌道とスピン半分を初めに $\ket{m_\ell=0,\uparrow}$ に用意し、$H=\xi\vb L\cdot\vb S/\hbar^2$ で時間発展させる。
$\ket{m_\ell=1,\downarrow}$ に変わる確率を求めよ。

:::solution
同じ $M=1/2$ の二状態を $a,b$ とし、$\ket{3/2,1/2}=\sqrt{2/3}a+\sqrt{1/3}b$、$\ket{1/2,1/2}=-\sqrt{1/3}a+\sqrt{2/3}b$ と書く。
二つのエネルギーは $\xi/2,-\xi$ である。
初状態をこの基底に展開し、時刻 $t$ に $b$ に射影すると

$$
\begin{aligned}
c_b&=\frac{\sqrt2}{3}\left(e^{-i\xi t/(2\hbar)}-e^{i\xi t/\hbar}\right), \\
P_b&=\frac89\sin^2\frac{3\xi t}{4\hbar}.
\end{aligned}
$$

全角運動量の $z$ 成分は一定だが、各成分の間で角運動量をやり取りする。
確率の最大が一に届かないのは、直積基底での二状態の対角エネルギーも異なるためである。
:::
::::

::::exercise{level=advanced}
任意方向のスピン半分の測定を行って結果を捨てる。
初期 Bloch ベクトル $\vb s$ と測定方向 $\vb n$ を用い、測定後の Bloch ベクトルを求めよ。

:::solution
$N=\vb n\cdot\boldsymbol\sigma$ と置くと $P_\pm=(I\pm N)/2$。
二つの結果を足して $\rho'=(\rho+N\rho N)/2$ を得る。
Pauli 行列の積の公式から

$$
N(\vb s\cdot\boldsymbol\sigma)N=[2(\vb n\cdot\vb s)\vb n-\vb s]\cdot\boldsymbol\sigma .
$$

従って $\vb s'=(\vb n\cdot\vb s)\vb n$。
測定軸に平行な偏りだけが残り、垂直成分は消える。
例えば初期 $\vb s=\vb e_z$、$\vb n=\vb e_x$ では $\rho'=I/2$ となる。
:::
::::

::::exercise{level=advanced}
次の計算には :eqref[eq-qm-spin-geometric-phase] を使える。$\theta=\pi/3$、$\Omega=2\pi\times1.0\,\mathrm{MHz}$、磁場の向きが一定の角速度で一周する時間を $1.0\,\mathrm{ms}$ とする。逆向きに回る二経路の位相差と、断熱近似を管理する行列要素の比を求めよ。

:::solution
一方向では $\gamma=-\pi/2$、逆向きは $+\pi/2$ なので位相差は $\pi$。同位相なら明るい出力が理想的には暗くなる。直交固有状態に漏れる係数は $|\langle u_-|\dot u_+\rangle|=|\dot\phi|\sin\theta/2$。これを準位差の角周波数 $\Omega$ で割ると $(10^3/10^6)\sqrt3/4=4.33\times10^{-4}$ である。強度の大小だけでなく回転速度が十分小さいことを確かめられる。
:::
::::


## 確認テスト

::::quiz{answer="0.75" tol="1e-6"}
スピン半分で $S^2/\hbar^2$ はいくつか。

:::explanation
$s(s+1)=(1/2)(3/2)=3/4$ である。
:::
::::

::::quiz
電子を $+z$ 磁場に置くと、どちらのスピンが低エネルギーか。

- [ ] $+z$ 向き
- [x] $-z$ 向き
- [ ] 同じ

:::explanation
電子の磁気モーメントはスピンと逆向きである。
$-z$ スピンの磁気モーメントが磁場とそろい、エネルギーが下がる。
:::
::::

::::quiz{answer="0.5" tol="1e-6"}
三重項の $M=0$ で第一スピンが上を向く確率はいくつか。

:::explanation
二つの直積状態の係数は $1/\sqrt2$ なので、その一方の確率は $1/2$ である。
:::
::::

::::quiz
中心力のスピン軌道項の角度因子を対角にする基底はどれか。

- [ ] 一般の位置固有基底
- [x] $L^2,S^2,J^2,J_z$ の同時固有基底
- [ ] $L_x,S_z$ の任意の同時固有基底

:::explanation
$\vb L\cdot\vb S=(J^2-L^2-S^2)/2$ と書けるため、合成基底で角因子が対角になる。
:::
::::

::::quiz{answer="2" tol="0.0001"}
三つのスピン半分を合成すると $J=1/2$ の既約表現は何組現れるか。

:::explanation
初めの二つを $J_{12}=0$ にした場合に一組、$J_{12}=1$ にした場合にも一組である。残る $J=3/2$ の四状態と合わせて全八状態を尽くす。
:::
::::

---
title: "電子-フォノン相互作用"
description: "格子変位から電子–フォノン結合を導き、電子の質量、対形成、ポーラロンに及ぼす効果を計算する。"
notation:
  - symbol: 'g_{\lambda}(\vb k,\vb q)'
    meaning: 電子–フォノン頂点
  - symbol: '\Omega_{\lambda\vb q}'
    meaning: フォノンのエネルギー
  - symbol: 'D^R'
    meaning: フォノンの遅延 Green 関数
  - symbol: '\alpha^2F,\lambda'
    meaning: 結合のスペクトルと積分強度
---

## この章の目標

- 格子変位を量子化して電子の散乱振幅を得る。
- フォノンの Green 関数と媒介相互作用を導く。
- Migdal 近似の小さい量と質量増大の条件を説明する。
- 独立振動子の変位変換からポーラロンを計算する。

## 電子-フォノン結合の導出と Fröhlich 模型

### イオンポテンシャルの微分

[格子振動](/physics/solid-state/phonons/) の固有モードを量子化し、その変位を電子の Hamiltonian に入れる。イオンの平衡位置を $\vb R_j$、変位を $\vb u_j$ とする。電子が見るポテンシャルを $\sum_j v(\vb r-\vb R_j-\vb u_j)$ と書くと、一次の摂動は $-\sum_j\vb u_j\cdot\nabla v(\vb r-\vb R_j)$ である。電子の Bloch 軌道に射影した行列要素を $d_{j\alpha}=\mel{\vb k'}{-\partial_\alpha v}{\vb k}$ とする。平行移動で $d_{j\alpha}=e^{i(\vb k-\vb k')\cdot\vb R_j}d_{0\alpha}$ となる。:eqref[eq-mb-ep-lattice] と同様、フォノンの変位に $e^{i\vb q\cdot\vb R_j}$ が掛かるためサイト和が $\vb k'=\vb k+\vb q$ を逆格子ベクトルの差まで選ぶ。多原子基底では各原子の質量 $M_\kappa$ と偏極ベクトル $e_{\kappa\alpha,\lambda}(\vb q)$ を用い、$g_\lambda=\sum_{\kappa\alpha}d_{\kappa\alpha}e_{\kappa\alpha,\lambda}/\sqrt{2M_\kappa\Omega_{\lambda\vb q}}$ とする。$d$ はエネルギーを長さで割った次元、零点変位は長さであり、$g$ はエネルギーである。

### 極性縦光学フォノン

極性媒質の縦光学振動では、偏極 $\vb P$ の発散が電荷密度 $-i\vb q\cdot\vb P$ を作る。Poisson 方程式の分母 $q^2$ と合わせるとポテンシャル振幅は $1/q$ に比例する。連続な三次元媒質で一つの分散しない縦光学分枝 $\Omega_{\rm LO}$ を使う :term[Fröhlich 模型]{en="Fröhlich model" yomi="ふれーりっひもけい"} では、箱の規格化を $1/\sqrt V$ とし、
$$
|g_q|^2=\frac{e^2\Omega_{\rm LO}}{2\epsilon_{\rm vac}q^2}
  (\epsilon_\infty^{-1}-\epsilon_s^{-1})
$$
と置く。この連続規格化でエネルギーの次元を持つのは $g_q/\sqrt V$ である。誘電定数 $\epsilon_s,\epsilon_\infty$ は相対値である。体積因子を分離した静的な媒介核 $-2|g_q|^2/\Omega_{\rm LO}$ を高周波の Coulomb 項へ足すと $e^2/(\epsilon_{\rm vac}\epsilon_s q^2)$ に戻るため、この規格化が定まる。SI では $\Omega_{\rm LO}$ は角周波数に $\hbar$ を掛けたエネルギーである。

:::example[誘電遮蔽の差]
$\epsilon_s=20,\epsilon_\infty=5$ なら逆数の差は $3/20$ である。同じ光学分枝上で $q$ が二倍になると $|g_q|^2$ は四分の一になる。音響分枝のひずみによる結合が小波数で弱くなる場合とは、偏極の長距離場の有無が異なる。
:::

:::example[音響モードの変位と頂点]
小波数で $\Omega_q=c|q|$ とし、ひずみに結合する頂点を用いる。波数を四倍にすると零点変位係数は半分、ひずみの微分は四倍になるので $|g_q|$ は二倍、$|g_q|^2$ は四倍になる。光学モードの $1/q$ とは異なる波数依存である。
:::

### 長波長結合と音響のスペクトル重み

三次元の縦音響分枝 $\Omega_q=cq$ と変形ポテンシャルの $|g_q|^2\propto q$ を取り、移る運動量の全方向を同じ重みで積分した裸の環境スペクトルを $J(\Omega)=\int\dd^3q\,(2\pi)^{-3}|g_q|^2\delta(\Omega-cq)$ とする。半径の測度は $4\pi q^2\dd q$ であり、デルタ関数の Jacobian は $1/c$ である。従って $J\propto\Omega^3$ となる。電子の自己エネルギーに入る $\alpha^2F$ では、さらに電子の初期・終状態を Fermi 面へ制限する。そのため $J$ のべきをそのまま自己エネルギーのべきとして使えない。

:::example[電子状態を指定する前の強度の比]
音速の線形範囲で $\Omega$ を二倍にすると、頂点の二乗は二倍、フォノン状態密度は四倍なので $J$ は八倍になる。電子の散乱ではエネルギー保存を別に課すため、同じ周波数変更で散乱率自体が八倍とは限らない。
:::

### Bloch 状態と偏極の規格化

単位胞の数を $N_c$ とし、その中のイオンを $\kappa$ で区別する。イオンの平衡位置を $\vb R_l+\vb d_\kappa$ とする。力定数行列を質量で両側から割って得る動力学行列の固有ベクトルを $\vb e_{\kappa\lambda}(\vb q)$ とし、$\sum_\kappa\vb e_{\kappa\lambda}^*\cdot\vb e_{\kappa\lambda'}=\delta_{\lambda\lambda'}$ と規格化する。時間反転対称で実の力定数を取り、分枝の基底を $\vb e_{\kappa\lambda}(-\vb q)=\vb e_{\kappa\lambda}(\vb q)^*$ となるよう選べる場合、
$$
\vb u_{l\kappa}=\frac1{\sqrt{N_c}}\sum_{\lambda\vb q}
\frac{\vb e_{\kappa\lambda}(\vb q)e^{i\vb q\cdot\vb R_l}}
{\sqrt{2M_\kappa\Omega_{\lambda\vb q}}}
(b_{\lambda\vb q}+b_{\lambda,-\vb q}^\dagger)
$$
である。セル内位相を $\vb d_\kappa$ にも掛ける Fourier 規約なら、偏極を逆の位相で変えれば変位は同じである。Bloch 状態の規格化には全結晶の $1/\sqrt{N_c}$ を含めるので、局所微分の行列要素から $d_{\kappa\alpha}^{ba}(\vb k,\vb q)=N_c\mel{b,\vb k+\vb q}{-\partial_\alpha v_\kappa}{a,\vb k}$ を作る。和 $\sum_l e^{i(\vb k-\vb k'+\vb q)\cdot\vb R_l}$ が $N_c$ を与え、最終ハミルトニアンには全体の $1/\sqrt{N_c}$ が残る。ここで $a,b$ はバンド番号であり、$\vb k+\vb q$ は Brillouin ゾーンへ折り返す。従って
$$
H_{\rm ep}=\frac1{\sqrt{N_c}}\sum_{ab\lambda\vb k\vb q\sigma}
g_\lambda^{ba}(\vb k,\vb q)c_{b,\vb k+\vb q,\sigma}^\dagger
c_{a,\vb k,\sigma}X_{\lambda\vb q}.
$$
$X_{\lambda\vb q}^\dagger=X_{\lambda,-\vb q}$ なので、Hermite 条件は $g_\lambda^{ab}(\vb k+\vb q,-\vb q)=g_\lambda^{ba}(\vb k,\vb q)^*$ となる。フォノンを吸収すると電子は結晶運動量 $\vb q$ を得る。放出の項ではフォノンの運動量が $-\vb q$ なので、電子の増分はやはり $\vb q$ である。局所ポテンシャルは自己無撞着な遮蔽を既に含むことがあり、その場合には同じ電子の静的遮蔽を頂点へもう一度掛けない。

### 一様並進と変形ポテンシャル

全イオンを同じ変位 $\vb u$ だけ動かすと、電子の一体 Hamiltonian は空間並進で元のものとユニタリ同値である。$T(\vb u)=e^{-i\vb u\cdot\vb p}$ と書けば、一次の変化は $-i[\vb u\cdot\vb p,H_0]$ である。対角の行列要素は固有エネルギーの差を含むので零になる。従って、絶縁体や適切に遮蔽された金属では、波長の長い音響モードに対する通常のバンド内結合は一様変位の項から始まらない。単純な等方媒質でバンド端のエネルギーを圧縮によって動かすなら、最初のスカラーは $D_{\rm ac}\nabla\cdot\vb u$ である。$D_{\rm ac}$ は :term[変形ポテンシャル]{en="deformation potential" yomi="へんけいぽてんしゃる"} で、ひずみの変化あたりのエネルギーを表す。SI の質量密度を $\rho_M$、縦音速を $c_s$ とすると、箱で一量子に対応する結合は
$$
M_{\vb q}=i D_{\rm ac}\,\vb q\cdot\hat e_{\vb q}
\sqrt{\frac{\hbar}{2\rho_M V c_s q}} .
$$
純横波では $\vb q\cdot\hat e_{\vb q}=0$ である。結晶の対称性が許すせん断の結合や圧電場まで無いことを意味しない。圧電や極性の長距離場では局所ひずみだけの展開自体を見直す。

:::example[電子散乱を持たない変位の成分]
$\vb q=q(1,0,0)$ の等方媒質で $\hat e=(1,0,0)$ なら、$|M|^2=D_{\rm ac}^2\hbar q/(2\rho_M V c_s)$ である。$\hat e=(0,1,0)$ ならこの項は零になる。縦と横の二つの音速の差を考慮する前に、結合するひずみの成分を選ぶ必要がある。
:::

### 極性結晶の誘電応答から結合を固定する

単一の極性光学分枝を持ち、低い周波数でイオンの相対変位 $u$ が電場 $E$ を受けるとする。単位体積あたり $n_I$ 個の双極子、相対質量 $M_r$、有効電荷 $Q$ に対して
$$
M_r(\omega_{\rm TO}^2-\omega^2)u=QE,\qquad
P_{\rm ion}=n_I Q u .
$$
電子だけの相対誘電率を $\epsilon_\infty$ とすれば $\epsilon(\omega)=\epsilon_\infty+n_IQ^2/[\epsilon_{\rm vac}M_r(\omega_{\rm TO}^2-\omega^2)]$ である。自由電荷のない純縦波で電束密度が零のとき、非零の縦電場は $\epsilon(\omega_{\rm LO})=0$ を要求する。従って
$$
\epsilon(\omega)=\epsilon_\infty\frac{\omega_{\rm LO}^2-\omega^2}{\omega_{\rm TO}^2-\omega^2},
\qquad
\frac{\epsilon_s}{\epsilon_\infty}=\frac{\omega_{\rm LO}^2}{\omega_{\rm TO}^2}.
$$
これを逆数にして縦 Coulomb 核 $e^2/[\epsilon_{\rm vac}q^2\epsilon(\omega)]$ に代入する。周波数に依らない項は高周波の Coulomb 核である。残りの項は $\Omega_{\rm LO}^2(\epsilon_\infty^{-1}-\epsilon_s^{-1})/[(\hbar\omega)^2-\Omega_{\rm LO}^2]$ を共通因子 $e^2/(\epsilon_{\rm vac}q^2)$ に掛けた形になる。これを $2\Omega_{\rm LO}|g_q|^2/[(\hbar\omega)^2-\Omega_{\rm LO}^2]$ と照合すれば、本文の Fröhlich 頂点が得られる。この導出では $\epsilon_\infty$ をフォノンより高く電子の共鳴より低い周波数での背景とみなした。

:::example[縦横分裂の確認]
$\epsilon_s/\epsilon_\infty=4$ の一分枝模型で $\hbar\omega_{\rm TO}=10\,\mathrm{meV}$ なら $\Omega_{\rm LO}=20\,\mathrm{meV}$ である。$\epsilon_s=\epsilon_\infty$ なら縦横の分裂は消え、媒介核の留数も零になる。
:::

### 二次の格子変位と一様並進の相殺

電子 Hamiltonian の Taylor 展開には一次の変位だけでなく、$\frac12\sum_j u_{j\alpha}u_{j\beta}\partial_\alpha\partial_\beta v(\vb r-\vb R_j)$ がある。調和フォノンの熱平均で一次の変位の平均は零でも、二次の変位は零点揺らぎと熱揺らぎを持つ。従って電子の準位の変化を変位の二次まで求めるときには、一次頂点を二度使う項と二次頂点を一度使う項の両方が出る。前者は電子とフォノンを含む中間状態の伝播を持ち、後者は同じ時刻の変位相関を持つ。しばしばそれぞれを Fan 項と :term[Debye–Waller 項]{en="Debye–Waller term" yomi="でばいわらーこう"} と呼ぶ。

一様な実の変位 $u$ による非縮退の固有エネルギーの二次微分を考え、$H'= \partial_u H$、$H''=\partial_u^2 H$ とする。固有方程式の一次微分を他の固有状態へ射影すると $\langle m|\partial_u n\rangle=\langle m|H'|n\rangle/(\epsilon_n-\epsilon_m)$ である。これを $\partial_u\epsilon_n=\langle n|H'|n\rangle$ の微分へ代入すると
$$
\partial_u^2\epsilon_n=\langle n|H''|n\rangle+
2\sum_{m\ne n}\frac{|\langle m|H'|n\rangle|^2}{\epsilon_n-\epsilon_m}.
$$
全イオンの一様変位では $H(u)=T(u)H(0)T(u)^\dagger$ であるため左辺は零になる。従って各項を単独で残すと、この並進に対する相殺を壊すことがある。縮退がある場合は固有空間を保った射影で処理する。静的な一様変位の証明であり、有限波数・有限周波数の準位の変化まで全て消えるということではない。

:::example[二準位のユニタリ変形]
$H_0=\diag(0,E)$ を実の角度 $u$ の回転行列で両側から変換する。$u=0$ の基底に対する $H''$ の期待値は $2E$、$H'$ の非対角要素の絶対値は $E$ である。そのため中間状態の項は $2E^2/(0-E)=-2E$ となり、二次微分は零になる。実空間での並進と同じユニタリ同値性の検査である。
:::

## フォノン Green 関数

### 振動子の交換子

無次元の変位 $X_{\vb q}=b_{\vb q}+b^\dagger_{-\vb q}$ に対し $D^R=-i\theta(t)\langle[X_{\vb q}(t),X_{-\vb q}]\rangle$ とする。自由発展と交換関係から交換子は $e^{-i\Omega t}-e^{i\Omega t}$ である。従って
$$
D_0(z)=\frac1{z-\Omega}-\frac1{z+\Omega}
=\frac{2\Omega}{z^2-\Omega^2},\qquad
D_0(i\nu_m)=-\frac{2\Omega}{\nu_m^2+\Omega^2}.
\label{eq-mb-phonon-propagator}
$$
座標そのものの Green 関数ではこれを $2M\Omega$ で割る。虚時間では正時間の平均が $-(n_B+1)e^{-\Omega\tau}-n_B e^{\Omega\tau}$ であり、 Fourier 変換は同じ $D_0(i\nu_m)$ を与える。電子の密度に $gX$ が結合すると、電子密度の変化は変位に対する力を変える。RPA 型の閉じ方では $D^{-1}=D_0^{-1}-|g|^2\chi_0$ になる。特に静的な $\chi_0<0$ は振動のばね定数を下げる。

:::example[フォノンの軟化]
$\Omega=2$、$|g|^2\chi_0=-1/4$ を周波数に依存しないと近似すると、極の二乗は $\Omega^2+2\Omega|g|^2\chi_0=3$ である。静的な二乗が負になった場合は、調和近似で置いた平衡位置が不安定である。
:::

:::example[虚周波数が分枝と等しいとき]
$\nu=\Omega=4\,\mathrm{meV}$ とすると $D_0(i\nu)=-1/\Omega=-0.25\,\mathrm{meV^{-1}}$ である。静的な $-2/\Omega=-0.5\,\mathrm{meV^{-1}}$ に対して大きさは半分である。
:::

### 電子の遅延応答が作る周波数と幅

弱い電子補正 $\Pi^R(z)=|g|^2\chi^R(z)$ を入れると、変位の逆伝播は $(z^2-\Omega^2)/(2\Omega)-\Pi^R(z)$ である。正周波数で $\Omega$ の近くに狭い共鳴があり、$\Pi$ の周波数依存が幅にわたり小さいと仮定する。$z=\Omega+\delta z$ を入れると最初の項は $\delta z+O(\delta z^2/\Omega)$ なので $\delta\Omega\simeq\Re\Pi^R(\Omega)$、半幅 $\gamma\simeq-\Im\Pi^R(\Omega)$ である。電子の粒子–正孔連続体に共鳴が入ると $\Im\chi^R<0$ になり、フォノンがそこでエネルギーを失う。静的極限 $\Pi^R(0)$ は軟化への手掛かりだが、動的幅を与えない。

:::example[狭い共鳴の見積もり]
$\Omega=20\,\mathrm{meV}$、$\Pi^R(\Omega)=(-1-i/4)\,\mathrm{meV}$ なら、ピークは約 $19\,\mathrm{meV}$、全幅は $1/2\,\mathrm{meV}$ である。ピーク変化が裸のエネルギーの五分の一を超える場合などは、線形近似の二次項と $\Pi$ の周波数依存を別に検査する。
:::

### 変位の熱揺らぎとスペクトルの規格化

スペクトルを $\rho_X(\omega)=-\Im D_0^R(\omega)/\pi$ と定義する。式 :eqref[eq-mb-phonon-propagator] は $\rho_X(\omega)=\delta(\omega-\Omega)-\delta(\omega+\Omega)$ を与える。積分は零であり、これは同時刻の $[X_{\vb q},X_{-\vb q}]=0$ に対応する。一方、一次モーメントは $2\Omega$ である。直接に $[\dot X_{\vb q},X_{-\vb q}]=-2i\Omega$ を計算するか、$D^R(z)=2\Omega/z^2+O(z^{-4})$ と展開して確かめられる。物理座標 $u=X/\sqrt{2M\Omega}$ ならその一次モーメントは $1/M$ になる。自己エネルギーが高周波でこの係数を変えない限り、相互作用した系で鋭い留数が減った分は他の周波数に移る。

調和振動子で $\langle X_{\vb q}X_{-\vb q}\rangle=1+2n_B(\Omega)=\coth(\beta\Omega/2)$ である。一般の平衡で変位の非順序相関の Fourier 重みは、正と負の周波数を別々に数えると $\rho_X(\omega)/(1-e^{-\beta\omega})$ になる。状態 $n$ から $m$ への遷移の重み $e^{-\beta E_n}$ と、逆の $e^{-\beta E_m}$ の差が交換子に入るためである。従って線の位置と幅だけを与える $D^R$ から、詳細釣り合いの係数を掛けて熱散乱の強度を求められる。

:::example[零点揺らぎと熱揺らぎ]
$\beta\Omega=\log3$ なら $n_B=1/2$ であり、$\langle X X^\dagger\rangle=2$ である。正周波数での非順序相関の重みは $3/2$、負周波数では $1/2$ となる。二本の差は $1$ なので、熱分布が変わっても自由振動子の交換子は同じである。
:::

### Kohn 異常を静的バブルから調べる

電子の密度応答はフォノンのばね定数を変える。その波数依存を見るため、一つの放物線帯の三次元の自由電子を零温度で取る。両スピンの状態密度を $N(0)$、$x=q/(2k_F)$ とすれば、[電子ガス](/physics/many-body/electron-gas/) の Lindhard 関数の静的実部は
$$
\chi_0(q,0)=-N(0) L(x),\qquad
L(x)=\frac12+\frac{1-x^2}{4x}\log\left|\frac{1+x}{1-x}\right| .
$$
分散が球対称なので、$q$ の方向を極軸に選んだ角度積分に二本の分母の対数が現れる。$x=1+\eta$ では $(1-x^2)/(4x)=-\eta/2+O(\eta^2)$、対数は $\log(2/|\eta|)+O(\eta)$ なので
$$
L(1+\eta)=\frac12-\frac\eta2\log\frac{2}{|\eta|}
+O(\eta^2|\log|\eta||).
$$
値は有限であるが傾きには対数がある。$g(q)$ と裸の分散がこの波数で滑らかなら、$2\Omega_q |g_q|^2\chi_0(q,0)$ を加えたフォノンの二乗周波数も同じ非解析性を持つ。これを :term[Kohn 異常]{en="Kohn anomaly" yomi="こーんいじょう"} と呼ぶ。有限の周波数、温度や寿命を入れると粒子–正孔の鋭い境界が丸くなる。一次元の自由放物線では二つの Fermi 点を結ぶため静的関数そのものが対数で発散する。そこでは調和振動の補正が大きくなり、低温で新しい平衡構造を調べる必要がある。

:::example[静的軟化と分散の傾き]
$|g_q|^2 N(0)$ が $q=2k_F$ で $\Omega_q/4$ に等しいとする。その点で $L(1)=1/2$ なので静的補正を入れた周波数の二乗は $3\Omega_q^2/4$ である。有限の軟化であっても、静的極限で計算した傾きには $\log|q-2k_F|$ が残る。軟化の大きさと傾きの非解析性は同じ条件で消えるとは限らない。
:::

### 粒子–正孔対へ崩壊できる周波数

伝導電子の粒子–正孔スペクトルをフォノンが取り込む例を、小 $q$ の三次元放物線帯で調べる。$\epsilon_{\vb k+\vb q}-\epsilon_{\vb k}\simeq v_Fq x$ と線形化し、$x=\hat k\cdot\hat q$ とする。零温度で占有数の差は $v_F qx\,\delta(\xi_{\vb k})$ であり、遅延バブルは両スピンの状態密度 $N(0)$ を用いて
$$
\chi_0^R(q,\omega)\simeq \frac{N(0)}2
\int_{-1}^1\dd x\,\frac{v_F q x}{\omega+i0-v_F qx}.
$$
遅延の虚部のデルタ関数を取ると
$$
-\Im\chi_0^R(q,\omega)
=\frac{\pi N(0)}2\frac{\omega}{v_Fq}\theta(v_Fq-|\omega|),\qquad \omega>0 .
$$
正の周波数の範囲では虚部は負であり、フォノンの半幅は正になる。縦音響の $\Omega=c_s q$ では $c_s<v_F$ なら連続体に入る。変形ポテンシャルの $|g_q|^2\propto q$ を使うと半幅は $q$ に比例し、$\gamma/\Omega$ は先頭で一定である。これに対し、固定した非零エネルギーを持つ光学分枝では、小さな $q$ のこの粒子–正孔連続体の外に位置する。有限温度、複数バンド、電子の有限寿命、非調和なフォノンの過程を入れた場合には、ここで零とした領域でも幅が現れ得る。

### 浴との結合が変える極の留数

$D^{-1}(z)=(z^2-\Omega^2)/(2\Omega)-\Pi(z)$ に対し、幅を無視できる単純な正エネルギーの根 $\widetilde\Omega$ を取る。その留数は
$$
R_X=\left[\frac{\widetilde\Omega}{\Omega}
-\left.\frac{\partial\Re\Pi}{\partial\omega}\right|_{\widetilde\Omega}\right]^{-1}
$$
である。$\Pi$ の虚部が小さければ $\gamma\simeq-R_X\Im\Pi^R(\widetilde\Omega)$ になる。正の留数を持つ安定な座標の枝でこの式を適用する。例えば定数の負の $\Pi$ は周波数を軟化させ、$R_X=\Omega/\widetilde\Omega$ によって変位のピークを強める。一方、ピークの面積とエネルギーを掛けた一次の総和則の寄与は $\widetilde\Omega R_X=\Omega$ である。周波数に依る $\Pi$ では鋭いピークから連続体へ重みが移るため、この一枝の積は必ずしも裸の値と等しくない。

:::example[周波数と留数が逆の倍率で動く場合]
無次元座標 $X$ の $D_0$ で定数の補正により $\widetilde\Omega=\Omega/2$ となったとする。正の線の留数は $2$、負の線は $-2$ である。一次モーメントは $2\Omega$ を保つ。元の座標 $X$ の定義を固定したまま評価した結果であり、新しい振動数で改めて $X$ を規格化した場合は表示する留数自体が変わる。
:::

## Migdal の定理

### 頂点補正の尺度

一ループの電子自己エネルギーはフォノン線一本と電子線一本を持つ。一方、最低次の頂点補正ではフォノン線に加え、運動量差 $\vb q$ を持つ二本の電子線が入る。:term[Migdal の定理]{en="Migdal theorem" yomi="みぐだるのていり"} の通常の金属での結論は、広い滑らかな Fermi 面、典型的な散乱運動量 $q\sim k_F$、フォノンエネルギー $\Omega\ll E_F$ のとき、先頭の相対頂点補正が典型的に $O(\lambda\Omega/E_F)$ であることだ。

小さい量の出方を見るため、角度 $u$ を速度の $\vb q$ への射影とするパッチで電子二線を $(i\omega-\xi)^{-1}$ と $[i(\omega+\nu)-\xi-v_Fqu]^{-1}$ に近似する。$\xi$ 積分が非零なのは二つの極が反対側にある幅 $|\nu|$ の周波数区間であり、留数の分母は $i\nu-v_Fqu$ である。対称な角度積分の大きさは $2\arctan(v_Fq/|\nu|)/(v_Fq)$ となる。そのため $\nu\sim\Omega$ の補正を自己エネルギーの無次元強度と比較すると $\Omega/(v_Fq)$ が残る。この評価は典型的な遊動金属のパッチでの評価であり、全角度と全次数の図の評価には :cite[mahan-many] を参照する。小さい運動量に結合が集中する場合、帯域が狭い場合、フォノンが不安定になる場合には同じ見積もりを適用できない。

:::example[断熱比]
$E_F=5\,\mathrm{eV}$、$\Omega=25\,\mathrm{meV}$、$\lambda=1$ なら尺度は $0.005$ である。$E_F=30\,\mathrm{meV}$ の小さいポケットに同じ $\Omega$ を入れると比は約 $0.83$ になり、小さい展開変数を得られない。
:::

:::example[小角の補正を点検する]
同じ $\Omega$ で $v_Fq=10\Omega$ のとき $\Omega/(v_Fq)=0.1$ である。$q$ が百倍小さくなると、この比は $10$ になり、典型角度に限った抑制で補正を評価できない。角度積分の正接の逆関数も、小波数ではその引数に比例する範囲へ移る。
:::

### 同じ近似で得る小波数の交差

二本の電子線を積分した後の角度平均を $I=\frac12\int_{-1}^1\dd u/(i\nu-v_Fqu)$ とする。$\nu>0$ では奇部分が消え、
$$
I=-i\nu\int_0^1\frac{\dd u}{\nu^2+(v_Fq)^2u^2}
=-\frac i{v_Fq}\arctan\frac{v_Fq}{\nu}.
$$
周波数窓の幅 $\nu$ を掛けると、無次元の因子は $\arctan x/x$、$x=v_Fq/\nu$ である。$x\gg1$ では $\pi/(2x)+O(x^{-2})$、$x\ll1$ では $1-x^2/3+O(x^4)$ となる。従って同じ留数計算でも前方散乱に集中すると抑制が消える。これは角度平均の交差であり、運動量依存した頂点全体の非摂動的評価ではない。

:::example[典型散乱と前方散乱の比較]
$x=1$ で因子は $\pi/4$、$x=10$ で約 $0.147$ である。$x=1/10$ では約 $0.997$ になる。フォノンのエネルギーが同じでも、電子の二本の殻がどれほど離れているかが結果を左右する。
:::

### 頂点の一ループ積分の周波数窓

スカラーの外部頂点 $g_{\rm ext}$ を挟んだ一ループの相対補正は、符号と滑らかな内部頂点因子を除くと
$$
T\sum_m\int\frac{\dd^d p}{(2\pi)^d}
|g_{\rm int}|^2D_0(i\omega_n-i\omega_m,\vb k-\vb p)
G_0(i\omega_m,\vb p)G_0(i\omega_m+i\nu,\vb p+\vb q)
$$
を持つ。広い帯の殻で状態密度を $\nu_0$ と見なし、二本の電子のずれを $a=v_Fqu$ とする。$i\omega_m-\xi$ と $i(\omega_m+\nu)-\xi-a$ の差で部分分数にすると、$\int_\R\dd\xi/(i\omega-\xi)=-i\pi\operatorname{sgn}\omega$ という対称切断の値から
$$
\int_\R\frac{\dd\xi}{(i\omega_m-\xi)(i\omega_m+i\nu-\xi-a)}
=\frac{i\pi[\operatorname{sgn}(\omega_m+\nu)-\operatorname{sgn}\omega_m]}{i\nu-a}.
$$
差は $\nu>0$ で $-\nu<\omega_m<0$ のみ非零になる。$\nu\sim\Omega$ の場合、内部フォノン線の大きさは $1/\Omega$ の桁であり、周波数積分の幅は $\Omega$ である。角平均で分母を $v_Fq$ の桁に置くと、相対補正は $\nu_0 |g_{\rm int}|^2/(v_Fq)$ の桁になる。$\lambda$ は同じ分枝で $\nu_0 |g_{\rm int}|^2/\Omega$ の桁なので、両者の比に $\Omega/(v_Fq)$ が残る。実際の定数は転送する $\vb q$、外部周波数、偏極と頂点の角度分布に左右される。この積分では曲率の必要なパッチや外部運動量が $2k_F$ に近い特異領域を除いている。

一様な密度に外場を掛ける動的極限では、準粒子の自己エネルギーと頂点は Ward 恒等式 :eqref[eq-mb-ward] で結ばれる。例えば運動量に依らない $\Sigma$ の傾きが $-\lambda$ なら、$\vb q=0$ で $\nu\to0$ とした密度頂点は $1+\lambda$ である。電子の密度に結合するフォノン頂点をこの極限で調べると、自己エネルギーから必要となる頂点は典型運動量の小さい補正とは異なる。格子の変形で密度以外の演算子に結合する頂点には、そのまま同じ Ward 恒等式を適用できない。正確な小波数の応答では、外場が何に結合するかを最初に定める。

:::example[断熱極限の質量依存]
ばね定数と電子帯域を固定し、同位体の質量 $M$ を四倍にする。調和近似で $\Omega\propto M^{-1/2}$ なのでエネルギーは半分である。局所ポテンシャルの変位微分が同じなら $|g|^2\propto 1/(M\Omega)$ で半分になり、$|g|^2/\Omega$ と $\lambda$ は同じである。通常の断熱比による頂点補正の尺度だけが半分になる。
:::

### どの電子帯域を断熱比の分母に置くか

広い帯で $q$ が典型的でも、実際の二線のエネルギー差が $v_Fq$ に比べて小さくなる面の部分がある。近接した二つの Fermi 面、角度の広い範囲で速度が小さい面、バンド端に近い小ポケットなどが例である。頂点を落とす前に、角度平均の前の電子の分母とフォノンのスペクトルの積を見て、重みの大きな領域のエネルギー差が大きいかを確かめる。球対称のパッチ模型で $\nu>0$ の角度を $u=0$ の周りの $|u|<u_c$ だけに制限し、その範囲で規格化すると、
$$
\frac{\nu}{2u_c}\int_{-u_c}^{u_c}
  \frac{\dd u}{i\nu-v_F q u}
=-i\,\frac{\arctan(v_F q u_c/\nu)}{v_F q u_c/\nu}
$$
となる。従って $q\sim k_F$ でも $u_c\ll\nu/(v_F q)$ なら、規格化した頂点の補正に抑制がない。この例は $\lambda$ の規格化でも同じ角領域に重みを置くとした比較である。角領域が単に狭くなって全結合まで弱くなる場合には、その弱さも同時に係数に入る。

頂点を無視した自己無撞着な近似では、一個のフォノン線を渡る電子の自己エネルギーを、幅や質量を持った電子線で計算し直す。各段階で対角自己エネルギーの大きさが $\lambda$ を生むので、$\lambda$ が $1$ の桁でも、典型的な相対頂点補正 $\lambda\Omega/E_{\rm el}$ が小さいなら両者を別の近似として扱える。$E_{\rm el}$ は実際に電子二線の分母が離れるエネルギーである。断熱比が小さいことはペアリングを弱く置く必要性を意味しない。一方、フォノン自身の軟化や多体の不安定性で基準にした正常相が変わると、同じ近似をそのまま反復できない。

:::example[同じ大きな波数で角度だけが変わる]
$v_F q=100\Omega$、$\nu=\Omega$ とする。全角度 $u_c=1$ での因子は約 $0.0156$ だが、$u_c=1/100$ なら $\pi/4$ になる。波数の大きさだけを見て頂点補正を捨てると、後者の狭い角度への集中を見落とす。
:::

## 電子の質量くりこみ

### フェルミ面で平均した自己エネルギー

状態密度を一定、温度を零とし、粒子側と正孔側にフォノン帯域より十分広い共通の切断 $W$ を取る。非対称な切断を用いる場合は定数のシフト $\Re\Sigma^R(0)$ を化学ポテンシャルに吸収する。Fermi 面で平均した $|g|^2$ とフォノンの状態密度の積を $\alpha^2F(\Omega)\ge0$ と書く。規格化は
$$
\lambda=2\int_0^\infty\frac{\alpha^2F(\Omega)}{\Omega}\dd\Omega
\label{eq-mb-ep-lambda}
$$
とする。一ループでは空いた終状態の積分 $\int_0^W\dd\xi/(\omega-\Omega-\xi+i0)$ と占有状態の積分 $\int_{-W}^0\dd\xi/(\omega+\Omega-\xi+i0)$ が入る。実部を取って $W\to\infty$ とすると、二つの端点の対数が相殺し $\log|(\omega-\Omega)/(\omega+\Omega)|$ が残る。虚部は $\xi$ の積分範囲に極が入ると $-\pi$ になる。従って
$$
\Re\Sigma^R(\omega)=\int_0^\infty\alpha^2F(\Omega)
\log\left|\frac{\omega-\Omega}{\omega+\Omega}\right|\dd\Omega,\quad
\Im\Sigma^R(\omega)=-\pi\int_0^{|\omega|}\alpha^2F(\Omega)\dd\Omega .
$$
零周波数の傾きは $-\lambda$ である。運動量依存を無視できる放物線帯では $Z=(1+\lambda)^{-1}$、$m^*/m=1+\lambda$ となる。

:::example[Einstein 分枝]
$\alpha^2F=A\delta(\Omega-\Omega_0)$、$A=\Omega_0/4$ なら $\lambda=1/2$ であり、$Z=2/3$ になる。$|\omega|<\Omega_0$ では零温度の虚部が零である。励起の減衰にはフォノン量子を作れるエネルギーが要る。
:::

### 熱フォノンによる Fermi 面の幅

幅広い帯域に対し、有限温度で一ループの二つの分母の分子は $n_B(\Omega)+1-f(\xi)$ と $n_B(\Omega)+f(\xi)$ になる。前者は電子がフォノンを出した後の空状態と熱フォノンの誘導過程、後者はフォノンの吸収と占有側からの逆過程を足したものである。虚部の二本のデルタ関数で電子のエネルギーは $\xi=\omega-\Omega,\omega+\Omega$ に固定される。従って
$$
-\Im\Sigma^R(0,T)=2\pi\int_0^\infty\dd\Omega\,\alpha^2F(\Omega)
[n_B(\Omega)+f(\Omega)] .
$$
$\omega=0$ で $1-f(-\Omega)=f(\Omega)$ を使った。フォノン帯域に対し高温でも電子帯域に対して低温なら $n_B(\Omega)\simeq T/\Omega$ が支配し、:eqref[eq-mb-ep-lambda] より $-\Im\Sigma^R(0,T)\simeq\pi\lambda T$ である。伝導度の輸送幅にはさらに速度変化の角度重みが入る。

:::example[Einstein モードを熱で励起する]
$\beta\Omega_0=\log3$ なら $n_B=1/2,f=1/4$ である。$\alpha^2F=A\delta(\Omega-\Omega_0)$ の場合は $-\Im\Sigma^R(0,T)=3\pi A/2$ となる。温度を零に下げれば両占有が消え、Fermi 面の幅も零に近づく。
:::

### 複数分枝の寄与を逆周波数で比較する

$\alpha^2F=\sum_j A_j\delta(\Omega-\Omega_j)$ なら $\lambda=2\sum_j A_j/\Omega_j$ である。一方、零温度の $-\Im\Sigma^R$ は分枝を横切るたびに $\pi A_j$ だけ増える。従って同じ幅の階段を作る高周波モードと低周波モードは傾きに異なる寄与をもつ。

:::example[二本から傾きと幅を復元する]
$\Omega_1=10\,\mathrm{meV}$、$\Omega_2=40\,\mathrm{meV}$、$A_1=A_2=2\,\mathrm{meV}$ とする。低い分枝の $\lambda_1=2/5$、高い分枝の $\lambda_2=1/10$ なので $Z=2/3$ である。$|\omega|=20\,\mathrm{meV}$ の $-\Im\Sigma$ は $2\pi\,\mathrm{meV}$、$|\omega|=50\,\mathrm{meV}$ では $4\pi\,\mathrm{meV}$ になる。傾きでは低い分枝が四倍効くが、自己エネルギー虚部の段差は同じである。
:::

### Fermi 面平均と Bloch–Grüneisen のべき

初期軌道と終軌道を Fermi 面に束縛したスペクトルを明示する。体積の規格化をすべて積分測度に入れ、一スピンの状態密度を $\nu_0=\int_{\vb k}\delta(\xi_{\vb k})$ と書く規約では
$$
\alpha^2F(\Omega)=\frac1{\nu_0}
\sum_\lambda\int_{\vb k}\int_{\vb k'}
\delta(\xi_{\vb k})\delta(\xi_{\vb k'})
|g_\lambda(\vb k,\vb k'-\vb k)|^2
\delta(\Omega-\Omega_{\lambda,\vb k'-\vb k}).
$$
取り得る波数はバンドを含めて数える。ここでは一つの等方的な面と縦分枝 $\Omega_q=c_s q$ を取る。$\vb k$ を面上の一点に固定し、$\vb q$ をその方向から角度 $u$ で測ると、放物線帯では $\xi_{\vb k+\vb q}=(q^2+2k_F q u)/(2m)$ である。従って $\int_{-1}^1\dd u\,\delta(\xi_{\vb k+\vb q})=(v_F q)^{-1}$ は $0<q<2k_F$ でだけ非零になる。測度 $q^2\dd q$、変形ポテンシャルの $|g|^2\propto q$、終状態の $1/q$ を合わせると、$\Omega$ が小さい領域で $\alpha^2F(\Omega)\propto\Omega^2$ である。従って本文の虚部の積分で、零温度の幅は $\omega^3$、熱による Fermi 面の幅は $T^3$ になる。温度の式では $\Omega=Tx$ と置き、$\int_0^\infty x^2[n_B(x)+f(x)]\dd x=(7/2)\zeta(3)$ が有限であることを使う。

電流の緩和では低角度の散乱が弱く効く。[不純物散乱と量子輸送](/physics/many-body/impurity-transport/) で弾性散乱の角度重みを導く。ここでは音速が $v_F$ より小さく、散乱前後の速度を面上の値で評価する。弾性近似の角度重みは $1-\hat k\cdot\hat k'=q^2/(2k_F^2)$ である。これをスペクトルに掛けた $\alpha_{\rm tr}^2F$ は小波数で $\Omega^4$ に比例する。電気抵抗を求める際は、フォノンが熱浴に緩和して電子の運動量を受け取れると仮定する。清浄な閉じた系で全運動量が保存される場合、平衡のフォノン分布を固定した輸送式だけでは直流抵抗を決められない。

線形化した衝突率の二次形式に、分布変形を速度の方向に比例するとして入れる。吸収の平衡重みを電子エネルギーで積分すると
$$
n_B(\Omega)\int_\R\dd\xi\,f(\xi)[1-f(\xi+\Omega)]
=\frac{\Omega}{4\sinh^2(\Omega/(2T))}
$$
である。実際、角括弧の積は $[f(\xi)-f(\xi+\Omega)]/(1-e^{-\Omega/T})$ で、その分子の積分は $\Omega$ になる。放出の平衡重みも詳細釣り合いで同じ値である。分布変分の $1/T$ と速度差の二乗を合わせると、輸送率の温度依存は
$$
\frac1T\int_0^{\Omega_{\max}}\dd\Omega\,
\frac{\Omega\,\alpha_{\rm tr}^2F(\Omega)}{\sinh^2(\Omega/(2T))}
$$
に比例する。周波数上端は三次元の単純な模型で $\Omega_{\max}=\min(2c_s k_F,\Omega_D)$ である。$\alpha_{\rm tr}^2F\propto\Omega^4$ を使えば $T\ll\Omega_{\max}$ で $T^5$ の :term[Bloch–Grüneisen 則]{en="Bloch–Grüneisen law" yomi="ぶろっほぐりゅーないぜんそく"} を得る。$T\gg\Omega_{\max}$ でも電子の縮退が保たれるなら $\sinh(\Omega/(2T))\simeq\Omega/(2T)$ なので $T$ に比例する。低温の無次元積分は $\int_0^\infty x^5/\sinh^2(x/2)\dd x=4(5!)\zeta(5)$ である。$1/[4\sinh^2(x/2)]=\sum_{j\ge1}j e^{-jx}$ を用いた。この計算は一つの球状の面での低角度の変形結合を前提としており、他の次元や散乱機構に指数を移すには測度と角度重みを計算し直す。

:::example[線幅と抵抗を同じ試料で比較する]
温度を低温域で二倍にすると、上の模型の Fermi 面の熱フォノン幅は八倍、熱フォノンによる抵抗は三十二倍である。高温域で両方が $T$ に比例しても、質量の $\lambda$ と輸送の重みを含めた $\lambda_{\rm tr}$ が一致するとは限らない。
:::

### 運動量微分を含む速度

電子の極は $\omega-\xi_{\vb k}-\Re\Sigma(\omega,\vb k)=0$ を満たす。面の法線で暗黙微分を取ると
$$
v_F^*=\frac{v_F+\partial_{k_\perp}\Re\Sigma(0,\vb k_F)}
{1-\partial_\omega\Re\Sigma(0,\vb k_F)} .
$$
ここで $\hbar=1$ であり、SI で波数微分を速度へ換算するには $1/\hbar$ を掛ける。状態密度の球面平均で求める $\lambda$ は分母の補正を近似するが、分子は小波数散乱などで同じ精度では決まらない。非等方的な電子帯では $\lambda$ も局所的な $\lambda(\vb k_F)$ として扱う。実験で測る変化として比較できるのは、同じ面の法線速度である。

:::example[自己エネルギーの二つの勾配]
$-\partial_\omega\Re\Sigma=1/2$、$\partial_{k_\perp}\Re\Sigma=v_F/5$ なら $v_F^*/v_F=(6/5)/(3/2)=4/5$ である。球状の放物線を基準にした質量比は $5/4$ になる。周波数勾配だけを使った $3/2$ とは、分子の変化だけ違う。
:::

### 有限帯域での傾きと非対称な実部

Einstein フォノンの Fermi 面平均として $\alpha^2F=A\delta(\Omega-\Omega_0)$ を取り、電子帯を $-W_-<\xi<W_+$ とする。正孔と粒子の寄与を別々に積分した零温度の実部は
$$
\Re\Sigma(\omega)=A\left[
 \log\left|\frac{\omega-\Omega_0}{\omega-\Omega_0-W_+}\right|
 +\log\left|\frac{\omega+\Omega_0+W_-}{\omega+\Omega_0}\right|\right].
$$
従って $\Re\Sigma(0)=A\log[(\Omega_0+W_-)/(\Omega_0+W_+)]$ である。一定項を化学ポテンシャルに移し、零周波数の傾きを取ると
$$
-\partial_\omega\Re\Sigma(0)= A\left[\frac2{\Omega_0}
-\frac1{\Omega_0+W_+}-\frac1{\Omega_0+W_-}\right].
$$
この式の $W_\pm$ を無限へ延ばして初めて $\lambda=2A/\Omega_0$ に戻る。帯の端はフォノンを出した後の電子の終状態にも効く。例えば正周波数で $\Omega_0<\omega<\Omega_0+W_+$ の間だけ $-\Im\Sigma=\pi A$ が残り、その上の周波数でこの分枝の寄与は消える。従って自己エネルギーの実部の切断と虚部の終状態を別の帯域で評価すると、解析性を保たない。

:::example[狭い対称帯での質量の傾き]
$W_+=W_-=\Omega_0$ なら傾きの絶対値は $A/\Omega_0$ であり、広帯域の半分である。$A=\Omega_0/2$ とすると、運動量微分を落としたとき $Z=2/3$ となる。広帯域にした同じ結合では $Z=1/2$ である。
:::

## フォノンを介した電子間引力

### 交換したエネルギーによる引力

虚時間の電子密度に $gX$ を線形に結合する。正定値の二次核 $A_X$ に対し $S_X=\tfrac12 X A_X X+gXn$ と書き、$X$ を平方完成で積分すると $-\tfrac12 g^2n A_X^{-1}n$ が作用に加わる。$A_X^{-1}=-D_0$ なので
$$
V_{\rm ph}(\vb q,i\nu_m)=-\frac{2|g_q|^2\Omega_{\vb q}}{\nu_m^2+\Omega_{\vb q}^2}.
\label{eq-mb-phonon-attraction}
$$
周波数差が小さい散乱を引き付け、大きい周波数差には効きにくい。静的な反発 $\mu>0$ を電子の高エネルギー殻から消去すると、一定状態密度の Cooper チャネルで $d\mu/d\ell=-\mu^2$ となる。実際、幅 $\dd\ell=\dd\xi/\xi$ の仮想対で二次の減少が $-\mu^2\dd\ell$ である。$\ell=\log(E_F/\Omega)$ まで積分すれば $\mu^*=\mu/[1+\mu\log(E_F/\Omega)]$ を得る。

:::example[遅延と反発]
$\mu=1/2$、$\log(E_F/\Omega)=4$ なら $\mu^*=1/6$ である。引力チャネルが $\lambda=1/3$ なら、単純な弱結合の差は $1/6$ になる。同じ周波数に両方を集中させるとこの高エネルギーの減少は得られない。
:::

:::example[二つの虚周波数の引力]
$\Omega=2\,\mathrm{meV}$、$|g|=1\,\mathrm{meV}$ の格子規格化なら $V_{\rm ph}(0)=-1\,\mathrm{meV}$ である。周波数差が $4\,\mathrm{meV}$ なら分母は $4+16$ なので $V_{\rm ph}=-0.2\,\mathrm{meV}$ となる。
:::

### 虚時間で見た遅延の長さ

零温度の :eqref[eq-mb-phonon-attraction] を $\nu$ で Fourier 変換する。$\int_\R(\dd\nu/2\pi)e^{-i\nu\tau}/(\nu^2+\Omega^2)=e^{-\Omega|\tau|}/(2\Omega)$ である。これは $\tau>0$ で下半平面の $-i\Omega$、$\tau<0$ で上半平面の $i\Omega$ を拾うことで得られる。従って相互作用は $-|g|^2e^{-\Omega|\tau|}$ となり、等時刻での頂点と、全時間を積分した静的値 $-2|g|^2/\Omega$ とを区別できる。周期 $\beta$ の場合は像を足し、$0<\tau<\beta$ で $-|g|^2[e^{-\Omega\tau}+e^{-\Omega(\beta-\tau)}]/(1-e^{-\beta\Omega})$ となる。二つの指数は周期像の正負を等比級数で和にしたものである。

:::example[遅延のある相互作用の平均]
零温度で $\Omega=20\,\mathrm{meV}$ の減衰時間は $\hbar/\Omega\simeq33\,\mathrm{fs}$ である。虚時間の距離をこの値の二倍へ伸ばすと核の絶対値は $e^{-2}$ 倍になる。静的に同じ強さを持つよう $|g|^2/\Omega$ を一定にしながら $\Omega$ を二倍にすると、時刻零の高さは二倍、時間幅は半分になる。
:::

### 反発の高エネルギー殻を消去する

:term[Coulomb 擬ポテンシャル]{en="Coulomb pseudopotential" yomi="くーろんぎぽてんしゃる"} の分母を、線形対方程式の二つの殻から求める。引力が $|\xi|<\Omega$ だけにあり、反発 $\mu$ はそこから $W\gg\Omega$ まで一定とする。小さい温度では低殻の Cooper 積分を $L_l$、高殻の積分を $L_h\simeq\log(W/\Omega)$ と書ける。殻内の対の応答をそれぞれ $\Delta_l,\Delta_h$ とすれば、定数の核に対して
$$
\Delta_l=(\lambda_{\rm pair}-\mu)L_l\Delta_l-\mu L_h\Delta_h,\qquad
\Delta_h=-\mu L_l\Delta_l-\mu L_h\Delta_h .
$$
この二式は、対の伝播の積分が殻ごとに $L$ となり、高い外部エネルギーには反発だけが残るという近似である。第二式を解いて第一式へ代入すると
$$
\Delta_l=\left[\lambda_{\rm pair}-\frac{\mu}{1+\mu L_h}\right] L_l\Delta_l .
$$
反発による低殻の振幅と符号を変えた高殻の振幅が、低エネルギーで反発の一部を打ち消す。$\mu>0,L_h>0$ で、$\mu^*$ は正で $\mu$ より小さい。$\lambda_{\rm pair}$ は選んだ角度チャネルの対の結合である。質量の式の $\lambda$ と同じ値を使うには等方的な核などの前提が必要となる。

### 高振動数の局所模型と低振動数の遅延

格子の一サイトの作用に対して $\nu$ の小さい範囲で展開すると、核は
$$
U+V_{\rm ph}(i\nu)=U-\frac{2g^2}{\Omega}+
\frac{2g^2\nu^2}{\Omega^3}+O(\nu^4)
$$
となる。$U$ は裸の局所反発である。$\Omega$ を電子の関係する周波数より大きく取る反断熱極限では、最初の二項の定数を有効な Hubbard 反発として扱える。このとき $2g^2/\Omega$ を固定して $\Omega\to\infty$ とすれば $g^2/\Omega^2$ は零へ行き、後で導くホッピングの指数補正も消える。一方、$\Omega$ が電子帯域より小さい場合、静的な負の $U+V_{\rm ph}(0)$ を全周波数に使うと高エネルギーで実際より強い引力を与える。低周波数の符号だけから転移温度を計算することもできない。高周波数側の反発や質量の周波数依存を対の方程式に残す必要がある。

:::example[静的相殺が有限周波数で外れる]
$\Omega=4$、$g=2$、$U=2$ の格子規格化では $U+V_{\rm ph}(0)=0$ である。$\nu=\Omega$ では引力が静的値の半分の $-1$ となるため、全核は $+1$ である。高エネルギーの殻の反発が静的相殺で消えないことが分かる。
:::

:::example[遅延を利用できる帯域幅]
$\mu=1/2$ とすると $L_h=0$ では $\mu^*=1/2$、$L_h=2$ では $\mu^*=1/4$ である。$\lambda_{\rm pair}=1/3$ なら二殻模型で前者にはこの対の不安定性がなく、後者の低殻には $1/12$ の正の係数が残る。高殻を消した後に初めてこの符号が分かる。
:::

### 一次元の軟化と二量化のエネルギー

スピン縮退のある一次元の最近接帯を半充填し、ホップを $t+\delta t,t-\delta t$ と交互に変える。この変形は $2k_F=\pi/a$ の波数を持つ。二サイトを単位胞とし、$t>0$、$|\delta t|\ll t$ とすると、一スピンの二バンドは
$$
E_\pm(k)=\pm|t+\delta t+(t-\delta t)e^{2ika}|
=\pm2\sqrt{t^2\cos^2(ka)+\delta t^2\sin^2(ka)} .
$$
折り畳んだゾーンの端でエネルギーの差は $4|\delta t|$ である。端の左右を同じ一次の座標 $\xi=v_F p$ で表せば、新しい占有枝は $-\sqrt{\xi^2+m_d^2}$、元の枝は $-|\xi|$ であり、$m_d=2|\delta t|$ である。狭い $\xi$ の殻の寄与を、正側の測度に両端と二スピンを含めた正の係数 $C$ で書くと
$$
\delta\mathcal E_{\rm el}=-C\int_0^W
(\sqrt{\xi^2+m_d^2}-\xi)\dd\xi
=-\frac C2 m_d^2\log\frac{2W}{m_d}+O(m_d^2).
$$
原始関数は $[\xi\sqrt{\xi^2+m_d^2}+m_d^2\operatorname{arsinh}(\xi/m_d)-\xi^2]/2$ であり、対数の係数はそこから出る。静的な格子の弾性費用は $K_d m_d^2/2$ と書けるので、零温度かつ無限鎖の断熱近似では、小さい変形で電子の負の対数が正の弾性費用より大きくなる。これが :term[Peierls 不安定性]{en="Peierls instability" yomi="ぱいえるすふあんていせい"} の一つの導出である。この一次元の模型では規則的な二量化を考える。有限温度、有限サイズ、量子的なイオンの揺らぎ、電子間反発を入れた実際の長距離秩序は別に判定する必要がある。

:::example[折り畳んだ帯の端にできる分離]
$t=100\,\mathrm{meV}$、$\delta t=2\,\mathrm{meV}$ なら、端の新しい二本は $\pm4\,\mathrm{meV}$ であり、差は $8\,\mathrm{meV}$ である。裸の帯域 $4t$ に比べると小さいが、元は同じエネルギーにあった二つの Fermi 点をつなぐので、低エネルギーの積分に影響する。
:::

## ポーラロン

### 局所振動子を変位させる

電子一個が格子を歪ませ、歪みとともに動く励起を :term[ポーラロン]{en="polaron" yomi="ぽーらろん"} と呼ぶ。局所 Holstein 模型 $H=H_t+\Omega\sum_j b_j^\dagger b_j+g\sum_jn_j(b_j+b_j^\dagger)$ で $\mathcal U=e^{(g/\Omega)\sum_j n_j(b_j^\dagger-b_j)}$ とすると、$[b^\dagger-b,b]=-1$ から $\mathcal U b_j\mathcal U^\dagger=b_j-gn_j/\Omega$ である。さらに $[n_j,c_j]=-c_j$ を順に使うと $c_j$ に $X_j=e^{-(g/\Omega)(b_j^\dagger-b_j)}$ が掛かる。結果は
$$
\mathcal U H\mathcal U^\dagger=
-t\sum_{\langle ij\rangle,\sigma}(c_i^\dagger c_j X_i^\dagger X_j+\mathrm{h.c.})
+\Omega\sum_j b_j^\dagger b_j-\frac{g^2}{\Omega}\sum_j n_j^2 .
$$
$t\ll\Omega$ で真空への射影を第一近似に取る。交換子が数であるため $X_j=e^{-u^2/2}e^{-u b_j^\dagger}e^{u b_j}$、$u=g/\Omega$ であり、別サイト二つの真空平均は $e^{-u^2}$ になる。そのため $t_{\rm eff}=te^{-u^2}$ である。$n_j^2=n_j+2n_{j\uparrow}n_{j\downarrow}$ から局所の二体引力は $-2g^2/\Omega$ になる。

:::example[移動度の減少]
$g/\Omega=2$ では $t_{\rm eff}/t=e^{-4}\simeq0.0183$ である。$\Omega=20\,\mathrm{meV}$ なら一電子の局所エネルギーの低下は $80\,\mathrm{meV}$ である。長距離 Fröhlich 模型では変位が多数のサイトに分布するため、この局所模型の指数をそのまま用いることはできない。
:::

:::example[二電子の局所束縛]
$\Omega=4,g=2$ の一サイトのエネルギーは、一電子で $-1$、二電子で $-4$ だけずれる。離れた二サイトを一粒子ずつ占めたシフトは $-2$ であり、同じサイトの裸の $U$ を足した後の差は $U-2$ である。従ってこの原子極限では $U=1$ なら同サイトが低く、$U=3$ なら別サイトが低い。
:::

### 原子極限の振動サテライト

電子が零個の振動子真空から一電子を加える。変位の大きさを $u=g/\Omega$ とすれば、電子のある固有状態は変位した振動子の数状態である。元の真空を逆向きに変位すると $e^{-u^2/2}\sum_{l\ge0}(\pm u)^l\ket l/\sqrt{l!}$ なので、付加のスペクトルに
$$
A_+(\omega)=e^{-u^2}\sum_{l\ge0}\frac{u^{2l}}{l!}\delta(\omega-\epsilon_d+g^2/\Omega-l\Omega)
$$
が現れる。総重みは指数級数で $1$ となり、平均の振動数は $u^2$ である。従って一次モーメントは $\epsilon_d-g^2/\Omega+u^2\Omega=\epsilon_d$ となる。最低ピークのシフトとスペクトルの平均のシフトは異なる。

:::example[零フォノンと一フォノンの重み]
$u=1$ では基底と一フォノンの重みはともに $e^{-1}$、二フォノンは $e^{-1}/2$ である。$\Omega=30\,\mathrm{meV}$ のとき最低線は裸の準位より $30\,\mathrm{meV}$ 低く、一フォノン線で裸の準位に戻る。幅を与える連続体を入れていないので、各線は離散的である。
:::

### 弱結合 Fröhlich ポーラロンのエネルギーと質量

一電子を真空中の放物線帯へ入れた極性結晶を考える。多電子の Fermi 面がないので Migdal の金属での議論の代わりに、真空のフォノンを一量子作る二次摂動を使う。$\hbar=1$ とし、$C=e^2(\epsilon_\infty^{-1}-\epsilon_s^{-1})/(4\pi\epsilon_{\rm vac})$ とおけば、体積で規格化した Fröhlich 頂点は $|g_q|^2=2\pi C\Omega/q^2$ である。裸の運動量 $\vb p$ の中間状態は電子の運動量 $\vb p-\vb q$ とフォノン $\vb q$ を持つ。そのエネルギー差は $\Omega+q^2/(2m)-\vb p\cdot\vb q/m$ である。従って二次のシフトは
$$
\delta E_{\vb p}=-\int\frac{\dd^3q}{(2\pi)^3}
\frac{2\pi C\Omega}{q^2[\Omega+q^2/(2m)-\vb p\cdot\vb q/m]} .
$$
$p^2/(2m)<\Omega$ なら分母は正である。さらに $p/\sqrt{2m\Omega}\ll1$ なら角度平均で一次項が消え、二次まで展開できる。$q_0=\sqrt{2m\Omega}$、$x=q/q_0$ として、定数項は $-Cq_0/\pi$ に $\int_0^\infty\dd x/(1+x^2)=\pi/2$ を掛けた値である。二次項は同じ係数に $4(p/q_0)^2/3$ と $\int_0^\infty x^2\dd x/(1+x^2)^3=\pi/16$ を掛ける。後者は $x=\tan\theta$ と置いて $\sin^2\theta\cos^2\theta$ を積分した。従って $\alpha=C\sqrt{m/(2\Omega)}$ とすると
$$
E_{\vb p}=-\alpha\Omega+
\frac{p^2}{2m}\left(1-\frac{\alpha}{6}\right)+O(\alpha^2,\alpha p^4),
\qquad \frac{m_{\rm pol}}m=1+\frac{\alpha}{6}+O(\alpha^2).
$$
SI で $\alpha$ には $1/\hbar$ が掛かる。バンドの幅が有限または $\Omega$ が波数によって変わる場合、連続模型で全波数にわたるこの積分の定数も変わる。強結合の連続ポーラロンでは電子が自分で作る偏極の空間形状を変分する必要があり、高次の展開と変分の比較は :cite[mahan-many] を参照する。

:::example[弱結合の二つの補正]
$\alpha=0.3$、$\Omega=30\,\mathrm{meV}$ なら二次摂動のエネルギーシフトは $-9\,\mathrm{meV}$、先頭の質量増加は $5\%$ である。局所模型の指数的な質量増加は、一サイトの大きな変位に伴う全く別の展開から出る。
:::

### 有限温度の局所重なりと仮想フォノン

変位演算子 $X=e^{-u(b^\dagger-b)}$ の調和振動子の熱平均を求める。正規順序に直すと $e^{-u^2/2}e^{-ub^\dagger}e^{ub}$ であり、熱状態で $\langle (b^\dagger)^l b^m\rangle=\delta_{lm} l!\,n_B^l$ である。占有の幾何分布を $l$ 回微分すればこの階乗モーメントを得る。従って $\langle X\rangle=e^{-u^2(n_B+1/2)}$ である。別サイトの独立な振動子の積を取ると $\langle X_i^\dagger X_j\rangle=e^{-u^2(2n_B+1)}$ になる。これは熱平均でのコヒーレントな一次のホップであり、残りの揺らぎが作る温度依存した非弾性ホップは別に計算する必要がある。

零温度で一電子が隣接へ移って戻る仮想過程では、変位 $X_i^\dagger X_j$ が二サイトに振動量子を作る。総数 $l$ の確率は二つの独立な Poisson 分布の和から $e^{-2u^2}(2u^2)^l/l!$ である。実際、総数を固定した和 $\sum_{j=0}^l u^{2j}u^{2(l-j)}/[j!(l-j)!]$ に二項定理を使える。真空へ射影した一次のホップを再度数えないために $l=0$ を除く。ループ一結合のエネルギーシフトは
$$
-\frac{t^2}{\Omega}\,e^{-2u^2}\sum_{l=1}^\infty
\frac{(2u^2)^l}{l!\,l}.
$$
分母の $l\Omega$ は実際の中間フォノンエネルギーを使った。$u$ が大きければ確率は $l\simeq 2u^2$ にあり、戻りのシフトの尺度は $-t^2/(2u^2\Omega)$ である。これに対して真空から真空へ直接動く $te^{-u^2}$ は指数的に小さくなる。二次の全ての経路が戻りだけになるわけではないため、分散の補正を得るには別のサイトへ出る経路も保持する。

:::example[熱平均と一ループの大きさ]
$u=1$、$\beta\Omega=\log3$ なら $\langle X_i^\dagger X_j\rangle=e^{-2}$ であり、零温度の $e^{-1}$ より $e^{-1}$ 倍になる。同じ $u=1$ の零温度の戻りでは、量子を一つだけ作る項は $-2e^{-2}t^2/\Omega$ である。符号は励起したフォノンが中間状態を高くすることから出る。
:::

### 格子の変位を変分して強結合の尺度を作る

連続 Fröhlich 模型の強い結合を考え、一電子波動関数を $\psi(\vb r)$ とし、フォノンを偏極のコヒーレント状態として変分する。包絡関数を規格化し、密度の Fourier 成分 $\rho_\psi(\vb q)=\int\dd^3r\,|\psi(\vb r)|^2e^{i\vb q\cdot\vb r}$ を取る。各分枝で $\langle b_{\vb q}\rangle=\beta_{\vb q}$ とすればエネルギーの二次と一次の項は $\Omega|\beta_{\vb q}|^2+2\Re[g_q\rho_\psi(\vb q)\beta_{\vb q}]/\sqrt V$ である。各 $\beta_{\vb q}$ を平方完成で最小化して
$$
\mathcal E[\psi]=
\int\frac{|\nabla\psi|^2}{2m}\dd^3r
-\int\frac{\dd^3q}{(2\pi)^3}\frac{|g_q|^2}{\Omega}|\rho_\psi(\vb q)|^2
$$
を得る。$|g_q|^2/\Omega=2\pi C/q^2$ を使うと第二項は $-\frac C2\int\dd^3r\,\dd^3r'\,|\psi(\vb r)|^2|\psi(\vb r')|^2/|\vb r-\vb r'|$ となる。ここで $C=e^2(\epsilon_\infty^{-1}-\epsilon_s^{-1})/(4\pi\epsilon_{\rm vac})>0$ は先の連続結合である。規格化した同じ形を長さ $R$ に伸縮すれば、運動項は $c_1/(mR^2)$、偏極項は $-c_2 C/R$ であり、$c_1,c_2>0$ は形だけで決まる。微分すると $R=2c_1/(m c_2 C)$、最小は $-m c_2^2 C^2/(4c_1)$ である。従ってこの変分上界の尺度は $-\alpha^2\Omega$ である。波動関数と格子を積で取った上界であり、移動質量や厳密な基底エネルギーをそのまま与えない。並進対称な模型では変分解の中心を一つに選んだため、運動量固有状態を作るにはその位置の重ね合わせも考える。

例えば $\psi(\vb r)=(\pi R^2)^{-3/4}e^{-r^2/(2R^2)}$ では $\langle|\nabla\psi|^2\rangle=3/(2R^2)$、$\rho_\psi(q)=e^{-q^2R^2/4}$ である。偏極項は $-(C/\pi)\int_0^\infty e^{-q^2R^2/2}\dd q=-C/(\sqrt{2\pi}R)$ である。従って $R=3\sqrt{2\pi}/(2mC)$、エネルギーは $-m C^2/(6\pi)=-\alpha^2\Omega/(3\pi)$ となる。長さが格子間隔に近づく場合、連続の有効質量と $1/q$ の小波数の近似をすべての $q$ に使えなくなる。

### 閉じ込めのエネルギーと局所模型の束縛を比べる

Holstein 模型では一サイトの変位を消した後の密度項が $-E_p n_j^2$、$E_p=g^2/\Omega$ である。二電子が同じサイトのとき、別サイトに一個ずつの状態に比べて $2E_p$ 低くなる。裸の局所反発 $U$ を足すと、原子極限の束縛の条件は $2E_p>U$ である。ホップを戻すと、空サイトを通って二電子が別々に動く利得もあるため、原子の境界をそのまま有限帯幅の境界とは呼べない。反断熱でフォノンの励起を落とし、強い負の有効 $U_{\rm eff}=U-2E_p$ と $|U_{\rm eff}|\gg t_{\rm eff}$ を仮定すれば、同じサイトの対が隣へ動くには二回の一電子ホップが必要である。途中の別サイトの配置の費用は $|U_{\rm eff}|$ で、二つのスピンの順序があるので対のホップの絶対値は $2t_{\rm eff}^2/|U_{\rm eff}|$ となる。フォノンも中間に励起するパラメータ域では、分母に $l\Omega$ の異なる経路が入る。

:::example[電子一個と局所対の帯幅]
$|U_{\rm eff}|=10t_{\rm eff}$ なら上の射影模型で対のホップの絶対値は $t_{\rm eff}/5$ である。同じ格子で一電子のホップより小さい。束縛エネルギーがさらに大きくなれば、二次の対移動の尺度はさらに小さくなる。
:::

## 演習問題

::::exercise{level=basic}
$\Omega=2$ の $D_0(0)$ を求めよ。

:::solution
式 :eqref[eq-mb-phonon-propagator] で分母は $-4$、分子は $4$ なので $-1$ である。
:::
::::

::::exercise{level=basic}
$|g|=2$、$\Omega=4$ の静的な媒介相互作用を求めよ。

:::solution
$-2|g|^2/\Omega=-2$ である。単位をエネルギーにそろえると二乗した結合をフォノンエネルギーで割った量になる。
:::
::::

::::exercise{level=basic}
$\lambda=2$、運動量に依存しない自己エネルギーと放物線帯で $m^*/m$ と $Z$ を求めよ。

:::solution
周波数微分が $-\lambda$ なので $Z=1/(1+2)=1/3$ である。速度が $Z$ 倍だから質量比は $3$ になる。
:::
::::

::::exercise{level=basic}
$\mu=1/3$、$\log(E_F/\Omega)=3$ の $\mu^*$ を求めよ。

:::solution
分母は $1+(1/3)\cdot3=2$ なので $\mu^*=1/6$ である。
:::
::::

::::exercise{level=advanced}
Fröhlich 結合が静的な誘電率を再現することを示せ。

:::solution
$e^2/(\epsilon_{\rm vac}\epsilon_\infty q^2)$ に $-2|g_q|^2/\Omega$ を足す。逆数差の $\epsilon_\infty^{-1}$ が消え、$e^2/(\epsilon_{\rm vac}\epsilon_s q^2)$ が残る。
:::
::::

::::exercise{level=advanced}
角度積分 $\int_{-1}^1\dd u/(i\nu-au)$ を $\nu,a>0$ で求めよ。

:::solution
$(-au-i\nu)/(a^2u^2+\nu^2)$ と書くと実部は奇関数で消える。虚部を積分して $-2i\arctan(a/\nu)/a$ を得る。
:::
::::

::::exercise{level=advanced}
Holstein 一サイトの $\Omega=4$、$g=2$ で一電子シフトと二体引力を求めよ。

:::solution
変位後のエネルギーは $-g^2n^2/\Omega=-n^2$ である。$n=1$ のシフトは $-1$、$n=2$ のシフト $-4$ から一電子シフトを二回引くと二体分は $-2$ になる。
:::
::::

::::exercise{level=advanced}
Einstein 分枝で $\Omega_0=3$、$\lambda=1$ の $\alpha^2F$ のデルタ重みと、$\omega=4$ の $\Im\Sigma^R$ を求めよ。

:::solution
正規化 $\lambda=2A/\Omega_0$ から $A=3/2$ を得る。積分の上端 $|\omega|=4$ は極 $3$ を含むため虚部は $-3\pi/2$ である。
:::
::::

::::exercise{level=advanced}
静的なフォノン補正 $\Pi=-3\Omega/8$ が一定で、$D^{-1}(z)=(z^2-\Omega^2)/(2\Omega)-\Pi$ とする。正の極の位置と留数を求め、一次モーメントを確認せよ。

:::solution
分母の零は $z^2=\Omega^2+2\Omega\Pi=\Omega^2/4$ なので、正の極は $\Omega/2$ である。逆関数の微分は $z/\Omega$ であり、正の点ではその逆数が $2$ となる。負の点では留数は $-2$ である。二本のエネルギー重みの和は $(\Omega/2)2+(-\Omega/2)(-2)=2\Omega$ で、裸の交換子の一次モーメントを保つ。
:::
::::

::::exercise{level=advanced}
三次元の等方的な音響分枝に対し、$|g_q|^2\propto q^3$ と仮定する。他の条件は変形ポテンシャルを用いた本文の模型と同じとし、零温度、低エネルギーの自己エネルギーの幅と低温のフォノン抵抗のべきを求めよ。

:::solution
Fermi 面の測度では $q^2$ に終状態の $1/q$ が掛かり、頂点と合わせて $\alpha^2F\propto q^4$ である。$\Omega=c_s q$ なので $\int_0^{|\omega|}\alpha^2F\dd\Omega$ は $|\omega|^5$ となる。輸送では角度の $q^2$ が付き $\alpha^2_{\rm tr}F\propto\Omega^6$ である。$T^{-1}\int\dd\Omega\,\Omega\alpha^2_{\rm tr}F/\sinh^2(\Omega/(2T))$ で $\Omega=Tx$ と置けば $T^7$ になる。熱浴としての格子の緩和を仮定した場合の抵抗である。
:::
::::

## 確認テスト

::::quiz
Migdal の見積もりに必要な条件を選べ。

- [x] 典型的な $v_F q$ がフォノンエネルギーより大きい。
- [ ] 任意の $q$ で頂点補正が厳密に零になる。
- [ ] 電子の状態密度が必ず零になる。

:::explanation
二電子線の位相空間を比較すると $\Omega/(v_Fq)$ が小さい領域で補正が抑えられる。典型的な散乱が小波数だけに集中するとこの前提を使えない。
:::
::::

::::quiz{answer="-2" tol="1e-6"}
$|g|=2$、$\Omega=4$ の静的なフォノン媒介相互作用を求めよ。

:::explanation
虚周波数を零とすると $-2|g|^2/\Omega=-2$ である。
:::
::::

::::quiz{answer="3" tol="1e-6"}
$\lambda=2$ の放物線帯で自己エネルギーの運動量微分を無視したときの質量比を答えよ。

:::explanation
周波数微分が $-2$ なので留数は $1/3$ であり、速度比の逆数は $3$ である。
:::
::::

::::quiz
局所 Holstein 模型を真空へ射影したときに得る量を選べ。

- [x] $t_{\rm eff}=te^{-(g/\Omega)^2}$
- [x] $-2g^2/\Omega$ の二体引力
- [ ] 速度が必ず自由電子の二倍になる

:::explanation
別サイトの変位演算子の真空重なりはそれぞれ $e^{-(g/\Omega)^2/2}$ である。エネルギーシフトの $n^2$ を一体と二体に分ければ二体引力を得る。
:::
::::

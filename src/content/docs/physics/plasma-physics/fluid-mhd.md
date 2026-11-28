---
title: 流体記述と二流体理論
description: 速度モーメントを取り、二流体の運動から一般化 Ohm の法則を導く。
notation:
  - symbol: 'f_s'
    meaning: 粒子種 s の分布関数
  - symbol: '\mathsf P_s'
    meaning: 圧力テンソル
  - symbol: '\vb j'
    meaning: 電流密度
  - symbol: '\vb U'
    meaning: 質量平均速度
---

## この章の目標

- 速度積分から保存則を導く。
- 圧力による流体ドリフトを解釈する。
- 二流体から一流体への近似を判定する。

## 速度モーメントと流体方程式の導出

粒子が同じ位置を異なる速度で通るとき、密度と平均速度で捉えられない運動量輸送が残る。この残りが圧力である。応力による運動量輸送の位置づけは [応力と保存則](/physics/continuum-mechanics/stress/) にある。

### 分布と平均

$f_s(\vb x,\vb v,t)\dd^3x\dd^3v$ を粒子数とする。$n_s=\int f_s\dd^3v$、$n_s\vb u_s=\int\vb v f_s\dd^3v$、$\vb c=\vb v-\vb u_s$ と置き、
$$
P_{s,ij}=m_s\int c_i c_j f_s\dd^3v,\qquad
\vb h_s=\frac{m_s}{2}\int c^2\vb c f_s\dd^3v
$$
を :term[圧力テンソル]{en="pressure tensor" yomi="あつりょくてんそる"} と :term[熱流束]{en="heat flux" yomi="ねつりゅうそく"} とする。$P_{s,ij}$ の単位は Pa、$\vb h_s$ の単位は $\mathrm{W\,m^{-2}}$ である。

粒子輸送式 $\partial_t f_s+\vb v\cdot\nabla f_s+\vb a_s\cdot\nabla_v f_s=C_s$ で、$\vb a_s=q_s(\vb E+\vb v\times\vb B)/m_s$、$C_s$ は種の粒子数を保存する衝突項とする。$\nabla_v\cdot\vb a_s=0$ であり、$f_s$ と各積の境界項は速度無限遠で消えるとする。重み $1$ で積分すると $\partial_t n_s+\nabla\cdot(n_s\vb u_s)=0$ である。重み $m_s v_i$ では、速度微分の部分積分が $-m_s\int a_{s,i} f_s\dd^3v$ を与える。二次の流束は $m_s n_s u_{s,i}u_{s,j}+P_{s,ij}$ であり、混合項は $\int\vb c f_s=0$ で消える。よって
$$
m_s n_s(\partial_t+\vb u_s\cdot\nabla)\vb u_s
=q_s n_s(\vb E+\vb u_s\times\vb B)-\nabla\cdot\mathsf P_s+\vb R_s,\quad
\vb R_s=\int m_s\vb v C_s\dd^3v .
\label{eq-plasma-species-momentum}
$$
重み $m_sv^2/2$ での流束は $(K_s+U_s)\vb u_s+\mathsf P_s\vb u_s+\vb h_s$ である。$K_s=m_s n_su_s^2/2$、$U_s=\tr\mathsf P_s/2$ である。力による供給は $q_s n_s\vb u_s\cdot\vb E$ となる。運動量式と $\vb u_s$ の内積を引けば
$$
\partial_t U_s+\nabla\cdot(U_s\vb u_s+\vb h_s)
=-\mathsf P_s:\nabla\vb u_s+\int (m_s c^2/2)C_s\dd^3v.
\label{eq-plasma-internal-energy}
$$

:::example[対向ビーム]
粒子の半数が $\vb v=u\hat{\vb x}$、残りが $-u\hat{\vb x}$ なら $\vb u_s=0$ だが $P_{xx}=m_s n_su^2$ である。横成分はゼロなのでスカラー圧力への置換で透過する運動量が変わる。
:::

### 二次中心モーメントの発展

温度の方向依存を追うには $Q_{ijk}=m_s\int c_i c_j c_k f_s\dd^3v$ が必要である。運動論式に $m_s c_i c_j$ を掛けて部分積分する。$\partial_t c_i=-\partial_t u_{s,i}$、$\partial_{x_k} c_i=-\partial_{x_k}u_{s,i}$ であり、$\int c_i f_s=0$ から一様な外力と平均流の加速度は中心モーメントに残らない。速度に一次の磁気力だけが二次モーメントを回転させる。よって種添字を省いて
$$
D_t P_{ij}+P_{ij}\partial_k u_k+P_{ik}\partial_k u_j+P_{jk}\partial_k u_i+\partial_k Q_{ijk}
=\frac qm(\epsilon_{i\ell n}P_{\ell j}+\epsilon_{j\ell n}P_{i\ell})B_n+\int mc_i c_j C\dd^3v
$$
となる。$\epsilon_{i\ell n}$ は Levi-Civita 記号である。跡をとると磁気項は反対称性により消え、式 :eqref[eq-plasma-internal-energy] の内部エネルギー式が戻る。

:::example[無衝突な一方向圧縮]
$\partial_x u_x=-a$、他の速度勾配ゼロ、熱流束・磁場・衝突なしで $D_t P_{xx}=3aP_{xx}$、$D_t P_{yy}=aP_{yy}$ である。密度は $D_t n=an$ なので $P_{xx}\propto n^3$、$P_{yy}\propto n$ となる。等方閉包にはこの差を消す機構が必要である。
:::

### 分布の圧力と熱流束を別々に測る

二つの狭いビームが速度 $v_+$ と $v_-$ で同じ軸上を進む。正側の粒子数割合を $a$、速度差を $d=v_+-v_-$ とすると、系の速度 $u=v_-+ad$ からの差は $(1-a)d$ と $-ad$ である。ゆえに二次中心モーメントは $a(1-a)d^2$、三次中心モーメントは $a(1-a)(1-2a)d^3$ である。横の速度分散はゼロとし、$P_{xx}/(mn)$ は二次モーメントに等しく、$h_x/(mn)$ は三次モーメントの半分に等しい。この軸方向の速度分布を、三次元の等方分布に置き換えると同じ内部エネルギーでも別の圧力テンソルになる。Cauchy の応力の規約 :eqref[eq-cont-traction] では運動論的な圧力テンソルは応力に負符号で入る。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $a$ | 正側ビームの粒子数割合 | 1 |
| $v_+$ | 正側ビームの速度 | $\mathrm{m\,s^{-1}}$ |
| $v_-$ | 負側ビームの速度 | $\mathrm{m\,s^{-1}}$ |

平均流速を $u$ とする。粒子数で重み付けする。

$$
u=a \left(v_{+} - v_{-}\right) + v_{-}
$$

ビーム間の速度差を $d$ とする。平均流とは独立な相対速度である。

$$
d=v_{+} - v_{-}
$$

単位質量当たりの軸圧力を $P_{xx}/(mn)$ とする。二本の偏差の二乗を平均する。

$$
P_{xx}/(mn)=a d^{2} \left(1 - a\right)
$$

単位質量当たりの内部エネルギーを $U/(mn)$ とする。この模型では横寄与がゼロである。

$$
U/(mn)=\frac{P_{xx}/(mn)}{2}
$$

単位密度当たりの熱流束を $h_x/(mn)$ とする。中心速度に対する運動エネルギーを中心速度で運ぶ。

$$
h_x/(mn)=\frac{P_{xx}/(mn) d \left(1 - 2 a\right)}{2}
$$

単位質量当たりの全運動エネルギーを $(K+U)/(mn)$ とする。交差項は平均してゼロになる。

$$
(K+U)/(mn)=U/(mn) + \frac{u^{2}}{2}
$$

:::example[対称で低速のビーム]
平均速度がゼロのときに圧力を求める。

入力値を次のように置く。

$$
\begin{aligned}
a&=\frac{1}{2} \\v_+&=20\,\mathrm{m\,s^{-1}} \\v_-&=-20\,\mathrm{m\,s^{-1}}
\end{aligned}
$$

これらを順に代入する。

平均流速は

$$
u=\left(-20\right) + \left(\frac{1}{2}\right) \left(- \left(-20\right) + \left(20\right)\right)=0\,\mathrm{m\,s^{-1}}
$$

である。

ビーム間の速度差は

$$
d=- \left(-20\right) + \left(20\right)=40\,\mathrm{m\,s^{-1}}
$$

である。

単位質量当たりの軸圧力は

$$
P_{xx}/(mn)=\left(40\right)^{2} \left(\frac{1}{2}\right) \left(1 - \left(\frac{1}{2}\right)\right)=400\,\mathrm{m^2\,s^{-2}}
$$

である。

単位質量当たりの内部エネルギーは

$$
U/(mn)=\frac{\left(400\right)}{2}=200\,\mathrm{m^2\,s^{-2}}
$$

である。

単位密度当たりの熱流束は

$$
h_x/(mn)=\frac{\left(400\right) \left(40\right) \left(1 - 2 \left(\frac{1}{2}\right)\right)}{2}=0\,\mathrm{m^3\,s^{-3}}
$$

である。

単位質量当たりの全運動エネルギーは

$$
(K+U)/(mn)=\frac{\left(0\right)^{2}}{2} + \left(200\right)=200\,\mathrm{m^2\,s^{-2}}
$$

である。

圧力は有限であるが熱流束は対称性で消える。
:::

:::example[希薄で高速のビーム]
熱流束がどちら向きか調べる。

入力値を次のように置く。

$$
\begin{aligned}
a&=\frac{1}{10} \\v_+&=100\,\mathrm{m\,s^{-1}} \\v_-&=-10\,\mathrm{m\,s^{-1}}
\end{aligned}
$$

これらを順に代入する。

平均流速は

$$
u=\left(-10\right) + \left(\frac{1}{10}\right) \left(- \left(-10\right) + \left(100\right)\right)=1\,\mathrm{m\,s^{-1}}
$$

である。

ビーム間の速度差は

$$
d=- \left(-10\right) + \left(100\right)=110\,\mathrm{m\,s^{-1}}
$$

である。

単位質量当たりの軸圧力は

$$
P_{xx}/(mn)=\left(110\right)^{2} \left(\frac{1}{10}\right) \left(1 - \left(\frac{1}{10}\right)\right)=1089\,\mathrm{m^2\,s^{-2}}
$$

である。

単位質量当たりの内部エネルギーは

$$
U/(mn)=\frac{\left(1089\right)}{2}=\frac{1089}{2}\,\mathrm{m^2\,s^{-2}}
$$

である。

単位密度当たりの熱流束は

$$
h_x/(mn)=\frac{\left(1089\right) \left(110\right) \left(1 - 2 \left(\frac{1}{10}\right)\right)}{2}=47916\,\mathrm{m^3\,s^{-3}}
$$

である。

単位質量当たりの全運動エネルギーは

$$
(K+U)/(mn)=\frac{\left(1\right)^{2}}{2} + \left(\frac{1089}{2}\right)=545\,\mathrm{m^2\,s^{-2}}
$$

である。

平均速度は小さい一方、中心速度の三乗では高速側が優勢である。
:::

:::example[多数派が正向きの場合]
同じ二速度でビームの割合だけを交換する。

入力値を次のように置く。

$$
\begin{aligned}
a&=\frac{9}{10} \\v_+&=100\,\mathrm{m\,s^{-1}} \\v_-&=-10\,\mathrm{m\,s^{-1}}
\end{aligned}
$$

これらを順に代入する。

平均流速は

$$
u=\left(-10\right) + \left(\frac{9}{10}\right) \left(- \left(-10\right) + \left(100\right)\right)=89\,\mathrm{m\,s^{-1}}
$$

である。

ビーム間の速度差は

$$
d=- \left(-10\right) + \left(100\right)=110\,\mathrm{m\,s^{-1}}
$$

である。

単位質量当たりの軸圧力は

$$
P_{xx}/(mn)=\left(110\right)^{2} \left(\frac{9}{10}\right) \left(1 - \left(\frac{9}{10}\right)\right)=1089\,\mathrm{m^2\,s^{-2}}
$$

である。

単位質量当たりの内部エネルギーは

$$
U/(mn)=\frac{\left(1089\right)}{2}=\frac{1089}{2}\,\mathrm{m^2\,s^{-2}}
$$

である。

単位密度当たりの熱流束は

$$
h_x/(mn)=\frac{\left(1089\right) \left(110\right) \left(1 - 2 \left(\frac{9}{10}\right)\right)}{2}=-47916\,\mathrm{m^3\,s^{-3}}
$$

である。

単位質量当たりの全運動エネルギーは

$$
(K+U)/(mn)=\frac{\left(89\right)^{2}}{2} + \left(\frac{1089}{2}\right)=4505\,\mathrm{m^2\,s^{-2}}
$$

である。

二次の値は前例と等しいが、熱流束は反転する。
:::

:::example[並進する対称ビーム]
ビームを平行移動した場合に内部量が変化するか求める。

入力値を次のように置く。

$$
\begin{aligned}
a&=\frac{1}{2} \\v_+&=120\,\mathrm{m\,s^{-1}} \\v_-&=80\,\mathrm{m\,s^{-1}}
\end{aligned}
$$

これらを順に代入する。

平均流速は

$$
u=\left(80\right) + \left(\frac{1}{2}\right) \left(\left(120\right) - \left(80\right)\right)=100\,\mathrm{m\,s^{-1}}
$$

である。

ビーム間の速度差は

$$
d=\left(120\right) - \left(80\right)=40\,\mathrm{m\,s^{-1}}
$$

である。

単位質量当たりの軸圧力は

$$
P_{xx}/(mn)=\left(40\right)^{2} \left(\frac{1}{2}\right) \left(1 - \left(\frac{1}{2}\right)\right)=400\,\mathrm{m^2\,s^{-2}}
$$

である。

単位質量当たりの内部エネルギーは

$$
U/(mn)=\frac{\left(400\right)}{2}=200\,\mathrm{m^2\,s^{-2}}
$$

である。

単位密度当たりの熱流束は

$$
h_x/(mn)=\frac{\left(400\right) \left(40\right) \left(1 - 2 \left(\frac{1}{2}\right)\right)}{2}=0\,\mathrm{m^3\,s^{-3}}
$$

である。

単位質量当たりの全運動エネルギーは

$$
(K+U)/(mn)=\frac{\left(100\right)^{2}}{2} + \left(200\right)=5200\,\mathrm{m^2\,s^{-2}}
$$

である。

平均流と全運動エネルギーが変わっても中心モーメントは変わらない。
:::


## 二流体方程式

電子と一価イオンを別々に追う。電子は軽いため、同じ力密度での速度変化が大きい。$q_i=e$、$q_e=-e$ を :eqref[eq-plasma-species-momentum] に入れて二本の運動方程式を得る。これらを二本の連続式と :ref[law-em-maxwell] に結合する。電荷密度は $e(n_i-n_e)$、電流は $\vb j=e(n_i\vb u_i-n_e\vb u_e)$ である。衝突は内部の運動量交換である場合 $\vb R_i=-\vb R_e$ を満たす。

電荷保存は二本の連続式に電荷を掛けた和である。Maxwell 方程式の Ampère–Maxwell 則の発散と Gauss 則からも同じ式が得られる。初期時刻に Gauss 則を課せば、この発展で保たれる。

:::example[共通流と相対流]
$n_i=n_e=n$、$\vb u_i=4\hat{\vb x}$、$\vb u_e=\hat{\vb x}$ を $\mathrm{m\,s^{-1}}$ で与えると $\vb j=3en\hat{\vb x}$ である。両種に同じ速度を足してもこの電流は変わらない。
:::

:::example[種の間の運動量交換]
同じ数密度の二種について力 $\vb R_e=m_e n\nu(\vb u_i-\vb u_e)$ のみを残す。差の時間微分は $-(1+m_e/m_i)\nu(\vb u_i-\vb u_e)$ であり、共通の質量流は一定である。$\nu$ は電子がほぼ固定したイオンに対し緩和する率であり、相対流の率には小さな反跳補正がある。
:::

### 種ごとの渦度と磁場

無衝突で圧力が $p_s=p_s(n_s)$ のとき、種の運動方程式を $n_s$ で割り、両辺の回転を取る。$\vb\omega_s=\nabla\times\vb u_s$ とし、$(\vb u_s\cdot\nabla)\vb u_s=\nabla(u_s^2/2)-\vb u_s\times\vb\omega_s$、$\nabla p_s/n_s=\nabla\int^{n_s}p_s'(n')\dd n'/n'$ を使う。Faraday 則により電場の回転を磁場の時間微分で置き換えると
$$
\partial_t(m_s\vb\omega_s+q_s\vb B)=
\nabla\times[\vb u_s\times(m_s\vb\omega_s+q_s\vb B)].
$$
これは :term[正準渦度]{en="canonical vorticity" yomi="せいじゅんうずど"} の種の流速による輸送式である。衝突や非順圧の圧力は右辺に新しい生成項を加える。場が流体に凍結するかを考えるとき、どの種の慣性を残したかで保存される組み合わせが変わる。

:::example[電子慣性が小さい場合]
電子について $m_e\abs{\vb\omega_e}\ll eB$ なら、その正準渦度は $-e\vb B$ で近似できる。電子圧力が順圧で衝突がなければ磁場が電子流で輸送される。イオンの渦度が無視できない尺度では、同じ議論から磁場がイオン流で輸送されるとは言えない。
:::

### 二種間摩擦の運動量と加熱

一様な二種に外力がなく、一定の摩擦係数 $A>0$ に対して $\rho_e\dot u_e=A(u_i-u_e)$、$\rho_i\dot u_i=-A(u_i-u_e)$ とする。差 $w=u_i-u_e$ は $\dot w=-A(\rho_e^{-1}+\rho_i^{-1})w$、質量平均 $u$ は定数である。相対運動のエネルギー密度 $K_{\rm rel}=\rho_e\rho_iw^2/[2(\rho_e+\rho_i)]$ の時間微分は $-Aw^2$ となる。全衝突項がエネルギーを保存するなら、この減少が両種の内部エネルギーを増やす。摩擦力の指定だけから各種への熱の分配は決まらない。熱平衡へ向かう温度交換は相対流の緩和より遅い場合がある。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $r$ | 質量密度比 $\rho_i/\rho_e$ | 1 |
| $u_{i0}$ | イオン初速 | $\mathrm{m\,s^{-1}}$ |
| $u_{e0}$ | 電子初速 | $\mathrm{m\,s^{-1}}$ |
| $s$ | 時間 $At/\rho_e$ | 1 |

保存する質量平均流を $u$ とする。掛け合わせて二種の摩擦項を相殺する。

$$
u=\frac{r u_{i0} + u_{e0}}{r + 1}
$$

相対流を $w$ とする。摩擦の両種の加速度が加わる。

$$
w=\left(- u_{e0} + u_{i0}\right) e^{- s \left(1 + \frac{1}{r}\right)}
$$

イオン速度を $u_i$ とする。質量平均と相対流を連立する。

$$
u_i=u + \frac{w}{r + 1}
$$

電子速度を $u_e$ とする。質量中心からの偏差には電子の軽さが現れる。

$$
u_e=- \frac{r w}{r + 1} + u
$$

残る相対流エネルギーを $K_{\rm rel}/\rho_e$ とする。密度の換算質量を使う。

$$
K_{\rm rel}/\rho_e=\frac{r w^{2}}{2 r + 2}
$$

内部エネルギーへの増分を $\Delta U/\rho_e$ とする。摩擦の仕事を初期から積分する。

$$
\Delta U/\rho_e=- K_{\rm rel}/\rho_e + \frac{r \left(- u_{e0} + u_{i0}\right)^{2}}{2 r + 2}
$$

:::example[電子流だけを初期に与える]
流速差の緩和で重い種が動く量を求める。

入力値を次のように置く。

$$
\begin{aligned}
r&=1836 \\u_{i0}&=0\,\mathrm{m\,s^{-1}} \\u_{e0}&=100\,\mathrm{m\,s^{-1}} \\s&=1
\end{aligned}
$$

これらを順に代入する。

保存する質量平均流は

$$
u=\frac{\left(0\right) \left(1836\right) + \left(100\right)}{\left(1836\right) + 1}=\frac{100}{1837}\,\mathrm{m\,s^{-1}}
$$

である。

相対流は

$$
w=\left(\left(0\right) - \left(100\right)\right) e^{- \left(1\right) \left(1 + \frac{1}{\left(1836\right)}\right)}\simeq -36.768\,\mathrm{m\,s^{-1}}
$$

である。

イオン速度は

$$
u_i=\frac{\left(- \frac{100}{e^{\frac{1837}{1836}}}\right)}{\left(1836\right) + 1} + \left(\frac{100}{1837}\right)\simeq 0.034421\,\mathrm{m\,s^{-1}}
$$

である。

電子速度は

$$
u_e=- \frac{\left(- \frac{100}{e^{\frac{1837}{1836}}}\right) \left(1836\right)}{\left(1836\right) + 1} + \left(\frac{100}{1837}\right)\simeq 36.802\,\mathrm{m\,s^{-1}}
$$

である。

残る相対流エネルギーは

$$
K_{\rm rel}/\rho_e=\frac{\left(\left(- \frac{100}{e^{\frac{1837}{1836}}}\right)\right)^{2} \left(1836\right)}{2 \left(1836\right) + 2}\simeq 675.57\,\mathrm{m^2\,s^{-2}}
$$

である。

内部エネルギーへの増分は

$$
\Delta U/\rho_e=\frac{\left(1836\right) \left(\left(0\right) - \left(100\right)\right)^{2}}{2 \left(1836\right) + 2} - \left(\frac{9180000}{1837 e^{\frac{1837}{918}}}\right)\simeq 4321.7\,\mathrm{m^2\,s^{-2}}
$$

である。

電子の初期運動量が小さいためイオンの流れは小さい。
:::

:::example[等質量の対向流]
熱へ移った割合を求める。

入力値を次のように置く。

$$
\begin{aligned}
r&=1 \\u_{i0}&=10\,\mathrm{m\,s^{-1}} \\u_{e0}&=-10\,\mathrm{m\,s^{-1}} \\s&=1
\end{aligned}
$$

これらを順に代入する。

保存する質量平均流は

$$
u=\frac{\left(-10\right) + \left(10\right) \left(1\right)}{\left(1\right) + 1}=0\,\mathrm{m\,s^{-1}}
$$

である。

相対流は

$$
w=\left(- \left(-10\right) + \left(10\right)\right) e^{- \left(1\right) \left(1 + \frac{1}{\left(1\right)}\right)}\simeq 2.7067\,\mathrm{m\,s^{-1}}
$$

である。

イオン速度は

$$
u_i=\left(0\right) + \frac{\left(\frac{20}{e^{2}}\right)}{\left(1\right) + 1}\simeq 1.3534\,\mathrm{m\,s^{-1}}
$$

である。

電子速度は

$$
u_e=\left(0\right) - \frac{\left(1\right) \left(\frac{20}{e^{2}}\right)}{\left(1\right) + 1}\simeq -1.3534\,\mathrm{m\,s^{-1}}
$$

である。

残る相対流エネルギーは

$$
K_{\rm rel}/\rho_e=\frac{\left(1\right) \left(\left(\frac{20}{e^{2}}\right)\right)^{2}}{2 \left(1\right) + 2}\simeq 1.8316\,\mathrm{m^2\,s^{-2}}
$$

である。

内部エネルギーへの増分は

$$
\Delta U/\rho_e=\frac{\left(1\right) \left(- \left(-10\right) + \left(10\right)\right)^{2}}{2 \left(1\right) + 2} - \left(\frac{100}{e^{4}}\right)\simeq 98.168\,\mathrm{m^2\,s^{-2}}
$$

である。

質量平均はゼロなので残った運動エネルギーは相対流だけである。
:::

:::example[共通流の中で緩和]
摩擦を加えても共通の並進が残るか調べる。

入力値を次のように置く。

$$
\begin{aligned}
r&=10 \\u_{i0}&=120\,\mathrm{m\,s^{-1}} \\u_{e0}&=80\,\mathrm{m\,s^{-1}} \\s&=2
\end{aligned}
$$

これらを順に代入する。

保存する質量平均流は

$$
u=\frac{\left(10\right) \left(120\right) + \left(80\right)}{\left(10\right) + 1}=\frac{1280}{11}\,\mathrm{m\,s^{-1}}
$$

である。

相対流は

$$
w=\left(\left(120\right) - \left(80\right)\right) e^{- \left(2\right) \left(1 + \frac{1}{\left(10\right)}\right)}\simeq 4.4321\,\mathrm{m\,s^{-1}}
$$

である。

イオン速度は

$$
u_i=\left(\frac{1280}{11}\right) + \frac{\left(\frac{40}{e^{\frac{11}{5}}}\right)}{\left(10\right) + 1}\simeq 116.77\,\mathrm{m\,s^{-1}}
$$

である。

電子速度は

$$
u_e=- \frac{\left(10\right) \left(\frac{40}{e^{\frac{11}{5}}}\right)}{\left(10\right) + 1} + \left(\frac{1280}{11}\right)\simeq 112.33\,\mathrm{m\,s^{-1}}
$$

である。

残る相対流エネルギーは

$$
K_{\rm rel}/\rho_e=\frac{\left(10\right) \left(\left(\frac{40}{e^{\frac{11}{5}}}\right)\right)^{2}}{2 \left(10\right) + 2}\simeq 8.929\,\mathrm{m^2\,s^{-2}}
$$

である。

内部エネルギーへの増分は

$$
\Delta U/\rho_e=\frac{\left(10\right) \left(\left(120\right) - \left(80\right)\right)^{2}}{2 \left(10\right) + 2} - \left(\frac{8000}{11 e^{\frac{22}{5}}}\right)\simeq 718.34\,\mathrm{m^2\,s^{-2}}
$$

である。

相対流が減っても質量平均は保存する。
:::

:::example[短時間のイオン流]
観測の初めに生じる加熱を求める。

入力値を次のように置く。

$$
\begin{aligned}
r&=100 \\u_{i0}&=10\,\mathrm{m\,s^{-1}} \\u_{e0}&=0\,\mathrm{m\,s^{-1}} \\s&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

保存する質量平均流は

$$
u=\frac{\left(0\right) + \left(100\right) \left(10\right)}{\left(100\right) + 1}=\frac{1000}{101}\,\mathrm{m\,s^{-1}}
$$

である。

相対流は

$$
w=\left(- \left(0\right) + \left(10\right)\right) e^{- \left(\frac{1}{100}\right) \left(1 + \frac{1}{\left(100\right)}\right)}\simeq 9.8995\,\mathrm{m\,s^{-1}}
$$

である。

イオン速度は

$$
u_i=\left(\frac{1000}{101}\right) + \frac{\left(\frac{10}{e^{\frac{101}{10000}}}\right)}{\left(100\right) + 1}\simeq 9.999\,\mathrm{m\,s^{-1}}
$$

である。

電子速度は

$$
u_e=- \frac{\left(100\right) \left(\frac{10}{e^{\frac{101}{10000}}}\right)}{\left(100\right) + 1} + \left(\frac{1000}{101}\right)\simeq 0.099497\,\mathrm{m\,s^{-1}}
$$

である。

残る相対流エネルギーは

$$
K_{\rm rel}/\rho_e=\frac{\left(100\right) \left(\left(\frac{10}{e^{\frac{101}{10000}}}\right)\right)^{2}}{2 \left(100\right) + 2}\simeq 48.515\,\mathrm{m^2\,s^{-2}}
$$

である。

内部エネルギーへの増分は

$$
\Delta U/\rho_e=\frac{\left(100\right) \left(- \left(0\right) + \left(10\right)\right)^{2}}{2 \left(100\right) + 2} - \left(\frac{5000}{101 e^{\frac{101}{5000}}}\right)\simeq 0.98997\,\mathrm{m^2\,s^{-2}}
$$

である。

相対エネルギーは振幅の二乗なので、その割合は流速差より速く減る。
:::


## 反磁性ドリフトと流体ドリフト

低周波で横慣性と横摩擦を無視し、$\mathsf P_s=p_s I$ とする。$\vb E+\vb u_s\times\vb B=\nabla p_s/(q_sn_s)$ の横成分を解くと
$$
\vb u_{s\perp}=\vb v_E+\vb u_{*s},\qquad
\vb u_{*s}=\frac{\vb b\times\nabla p_s}{q_sn_sB}.
\label{eq-plasma-diamagnetic}
$$
後者は :term[反磁性ドリフト]{en="diamagnetic drift" yomi="はんじせいどりふと"} である。一様磁場中で外力のない一粒子は閉じた円を描く。それでも密度勾配があれば、ある位置を通過する旋回の上半分と下半分で粒子数が違い、局所平均速度が残る。

等温で案内中心の密度が $n(\vb R)$ の場合、位置 $\vb r$ での分布を $f_0(v)n(\vb r+\vb v\times\vb b/\Omega)$ と書く。正規化した等方分布で $\int v_i v_j f_0\dd^3v=(T/m)\delta_{ij}$ を使い、一次の速度モーメントをとると $\vb\Gamma_*=(T/(m\Omega))\vb b\times\nabla n$ である。これは :eqref[eq-plasma-diamagnetic] の $n\vb u_*$ である。磁場が一様なら $\nabla\cdot(\vb b\times\nabla p_s)=0$ である。

:::example[圧力勾配が同じ向きの場合]
$\vb B=B\hat{\vb z}$、$\partial_x(p_i+p_e)>0$ なら $\vb j_*=\vb b\times\nabla(p_i+p_e)/B$ は $+\hat{\vb y}$ を向く。イオンと電子の速度は逆向きで、電流は足し合わさる。
:::

:::example[磁化電流としての反磁性]
一様 $\vb B$ で等方な分布を考える。旋回の磁気双極子の密度は $\vb M=-p\vb b/B$ である。$\nabla\times\vb M=-(\nabla p)\times\vb b/B=\vb b\times\nabla p/B$ は式 :eqref[eq-plasma-diamagnetic] に電荷密度 $q_sn_s$ を掛けた電流に一致する。粒子が案内中心を横切って移流しなくても電流が現れる。
:::

### 反磁性流の発散と磁化電流

一様場 $\vb B=B\hat{\vb z}$ で $p_s=p_s(x,y)$ とする。最低次の力のつりあいから粒子流束は $n_s\vb u_{*s}=\hat{\vb z}\times\grad p_s/(q_sB)$ である。発散は混合偏微分の差だけなのでゼロになる。他方、旋回する粒子の磁気モーメントは磁場と逆向きであり、等方圧で磁化は $\vb M=-p_s\hat{\vb z}/B$ となる。$\rot\vb M=\hat{\vb z}\times\grad p_s/B$ は反磁性電流に等しい。一様場を外れれば曲率と強度勾配の案内中心流も加える必要がある。以下は局所勾配を使う例であり、境界の外へ流れる粒子数をこの流速だけから評価しない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $G_x$ | 種の圧力の $x$ 勾配 | $\mathrm{Pa\,m^{-1}}$ |
| $G_y$ | 種の圧力の $y$ 勾配 | $\mathrm{Pa\,m^{-1}}$ |
| $B$ | 一様な磁束密度 | $\mathrm{T}$ |
| $Q$ | 符号付きの種の電荷密度 | $\mathrm{C\,m^{-3}}$ |

反磁性流の $x$ 成分を $u_{*x}$ とする。外積で $y$ 勾配に負符号が付く。

$$
u_{*x}=- \frac{G_{y}}{B Q}
$$

反磁性流の $y$ 成分を $u_{*y}$ とする。電荷符号で流れが反転する。

$$
u_{*y}=\frac{G_{x}}{B Q}
$$

反磁性電流の $x$ 成分を $j_{*x}$ とする。流れに電荷密度を掛ける。

$$
j_{*x}=Q u_{*x}
$$

反磁性電流の $y$ 成分を $j_{*y}$ とする。同じ圧力なら種の電荷符号に依存しない。

$$
j_{*y}=Q u_{*y}
$$

反磁性流に沿った圧力変化を $\vb u_*\cdot\grad p_s$ とする。流れは等圧線に接する。

$$
\vb u_*\cdot\grad p_s=G_{x} u_{*x} + G_{y} u_{*y}
$$

電流密度の大きさを $\abs{\vb j_*}$ とする。二つの直交成分から求める。

$$
\abs{\vb j_*}=\sqrt{j_{*x}^{2} + j_{*y}^{2}}
$$

:::example[一方向のイオン圧]
圧力勾配に垂直な流れを求める。

入力値を次のように置く。

$$
\begin{aligned}
G_x&=4\,\mathrm{Pa\,m^{-1}} \\G_y&=0\,\mathrm{Pa\,m^{-1}} \\B&=2\,\mathrm{T} \\Q&=1\,\mathrm{C\,m^{-3}}
\end{aligned}
$$

これらを順に代入する。

反磁性流の $x$ 成分は

$$
u_{*x}=- \frac{\left(0\right)}{\left(1\right) \left(2\right)}=0\,\mathrm{m\,s^{-1}}
$$

である。

反磁性流の $y$ 成分は

$$
u_{*y}=\frac{\left(4\right)}{\left(1\right) \left(2\right)}=2\,\mathrm{m\,s^{-1}}
$$

である。

反磁性電流の $x$ 成分は

$$
j_{*x}=\left(0\right) \left(1\right)=0\,\mathrm{A\,m^{-2}}
$$

である。

反磁性電流の $y$ 成分は

$$
j_{*y}=\left(1\right) \left(2\right)=2\,\mathrm{A\,m^{-2}}
$$

である。

反磁性流に沿った圧力変化は

$$
\vb u_*\cdot\grad p_s=\left(0\right) \left(2\right) + \left(0\right) \left(4\right)=0\,\mathrm{Pa\,s^{-1}}
$$

である。

電流密度の大きさは

$$
\abs{\vb j_*}=\sqrt{\left(0\right)^{2} + \left(2\right)^{2}}=2\,\mathrm{A\,m^{-2}}
$$

である。

$y$ 方向の電流にのみ寄与する。
:::

:::example[同じ電子圧]
電子圧に替えたときに電流が打ち消すか調べる。

入力値を次のように置く。

$$
\begin{aligned}
G_x&=4\,\mathrm{Pa\,m^{-1}} \\G_y&=0\,\mathrm{Pa\,m^{-1}} \\B&=2\,\mathrm{T} \\Q&=-1\,\mathrm{C\,m^{-3}}
\end{aligned}
$$

これらを順に代入する。

反磁性流の $x$ 成分は

$$
u_{*x}=- \frac{\left(0\right)}{\left(-1\right) \left(2\right)}=0\,\mathrm{m\,s^{-1}}
$$

である。

反磁性流の $y$ 成分は

$$
u_{*y}=\frac{\left(4\right)}{\left(-1\right) \left(2\right)}=-2\,\mathrm{m\,s^{-1}}
$$

である。

反磁性電流の $x$ 成分は

$$
j_{*x}=\left(-1\right) \left(0\right)=0\,\mathrm{A\,m^{-2}}
$$

である。

反磁性電流の $y$ 成分は

$$
j_{*y}=\left(-1\right) \left(-2\right)=2\,\mathrm{A\,m^{-2}}
$$

である。

反磁性流に沿った圧力変化は

$$
\vb u_*\cdot\grad p_s=\left(-2\right) \left(0\right) + \left(0\right) \left(4\right)=0\,\mathrm{Pa\,s^{-1}}
$$

である。

電流密度の大きさは

$$
\abs{\vb j_*}=\sqrt{\left(0\right)^{2} + \left(2\right)^{2}}=2\,\mathrm{A\,m^{-2}}
$$

である。

流れは逆だが電流は同符号になる。
:::

:::example[斜めの勾配]
圧力に沿う変化を内積で検査する。

入力値を次のように置く。

$$
\begin{aligned}
G_x&=3\,\mathrm{Pa\,m^{-1}} \\G_y&=4\,\mathrm{Pa\,m^{-1}} \\B&=1\,\mathrm{T} \\Q&=1\,\mathrm{C\,m^{-3}}
\end{aligned}
$$

これらを順に代入する。

反磁性流の $x$ 成分は

$$
u_{*x}=- \frac{\left(4\right)}{\left(1\right)^{2}}=-4\,\mathrm{m\,s^{-1}}
$$

である。

反磁性流の $y$ 成分は

$$
u_{*y}=\frac{\left(3\right)}{\left(1\right)^{2}}=3\,\mathrm{m\,s^{-1}}
$$

である。

反磁性電流の $x$ 成分は

$$
j_{*x}=\left(-4\right) \left(1\right)=-4\,\mathrm{A\,m^{-2}}
$$

である。

反磁性電流の $y$ 成分は

$$
j_{*y}=\left(1\right) \left(3\right)=3\,\mathrm{A\,m^{-2}}
$$

である。

反磁性流に沿った圧力変化は

$$
\vb u_*\cdot\grad p_s=\left(-4\right) \left(3\right) + \left(3\right) \left(4\right)=0\,\mathrm{Pa\,s^{-1}}
$$

である。

電流密度の大きさは

$$
\abs{\vb j_*}=\sqrt{\left(-4\right)^{2} + \left(3\right)^{2}}=5\,\mathrm{A\,m^{-2}}
$$

である。

二成分の内積がちょうど相殺する。
:::

:::example[磁場を倍にする]
同じ勾配で磁化電流の強さを比較する。

入力値を次のように置く。

$$
\begin{aligned}
G_x&=3\,\mathrm{Pa\,m^{-1}} \\G_y&=4\,\mathrm{Pa\,m^{-1}} \\B&=2\,\mathrm{T} \\Q&=1\,\mathrm{C\,m^{-3}}
\end{aligned}
$$

これらを順に代入する。

反磁性流の $x$ 成分は

$$
u_{*x}=- \frac{\left(4\right)}{\left(1\right) \left(2\right)}=-2\,\mathrm{m\,s^{-1}}
$$

である。

反磁性流の $y$ 成分は

$$
u_{*y}=\frac{\left(3\right)}{\left(1\right) \left(2\right)}=\frac{3}{2}\,\mathrm{m\,s^{-1}}
$$

である。

反磁性電流の $x$ 成分は

$$
j_{*x}=\left(-2\right) \left(1\right)=-2\,\mathrm{A\,m^{-2}}
$$

である。

反磁性電流の $y$ 成分は

$$
j_{*y}=\left(1\right) \left(\frac{3}{2}\right)=\frac{3}{2}\,\mathrm{A\,m^{-2}}
$$

である。

反磁性流に沿った圧力変化は

$$
\vb u_*\cdot\grad p_s=\left(-2\right) \left(3\right) + \left(4\right) \left(\frac{3}{2}\right)=0\,\mathrm{Pa\,s^{-1}}
$$

である。

電流密度の大きさは

$$
\abs{\vb j_*}=\sqrt{\left(-2\right)^{2} + \left(\frac{3}{2}\right)^{2}}=\frac{5}{2}\,\mathrm{A\,m^{-2}}
$$

である。

流れと電流は半分になる。
:::


## 閉包問題と状態方程式

未知量の数を数えると、運動量方程式には二次中心モーメントが、エネルギー式には三次中心モーメントが現れる。指定した次数で上位モーメントを既知量から与える操作を :term[閉包]{en="closure" yomi="へいほう"} という。

### 等温閉包と断熱閉包

熱交換が速く $T$ が指定される模型では $p=nT$ である。等方な三自由度の模型で、熱流束と内部加熱を無視すると $U_s=3p_s/2$ である。式 :eqref[eq-plasma-internal-energy] を展開すると $D_s p_s/Dt=-(5/3)p_s\nabla\cdot\vb u_s$ が出る。連続式で $D_s n_s/Dt=-n_s\nabla\cdot\vb u_s$ なので $D_s(p_s n_s^{-5/3})/Dt=0$ である。

圧力が磁場に対し非等方なら、$p_\parallel,p_\perp$ を別に追う必要がある。衝突と熱流束を無視し、低周波で粒子の磁気モーメントが保存される模型では $p_\perp/(nB)$ が流体要素で保存される。流体要素は断面積 $A\propto B^{-1}$、長さ $\ell\propto B/n$ の磁束管として変形する。平行運動に作用保存を適用できる均一な往復模型では $v_\parallel\ell$ が一定なので $p_\parallel B^2/n^3$ が一定となる。一般の無衝突な系では熱流束や共鳴を無視できるか別に確認する。

:::example[同じ密度変化]
密度を八倍にしたとき圧力は等温で八倍、三自由度断熱で $8^{5/3}=32$ 倍である。熱を逃がせない模型ほど圧縮への復元力が強い。
:::

:::example[平行圧縮と横圧縮]
二重断熱の模型で長さだけが半分の磁束管では $B$ は一定、$n$ は二倍になる。この場合 $p_\parallel$ は八倍、$p_\perp$ は二倍である。長さ一定で断面積が半分のときは $B$ と $n$ が二倍なので、$p_\parallel$ は二倍、$p_\perp$ は四倍である。
:::

### 圧力テンソルからの二重断熱式

衝突を無視し、第三中心モーメントの発散を落とす。場の変化が旋回より遅く、圧力テンソルは磁場の周りの回転に対称であると仮定する。このとき $\mathsf P=p_\perp(I-\vb b\vb b)+p_\parallel\vb b\vb b$ である。ここで $\Theta=\nabla\cdot\vb u$、$K=b_i b_j\partial_j u_i$ とする。二次モーメントの方程式を $b_i b_j$ と縮約すれば磁気項は消え、$D_t p_\parallel=-p_\parallel(\Theta+2K)$ となる。$D_t(\vb b\cdot\vb b)=0$ なので、$\vb b$ の時間微分はこの縮約に寄与しない。

同じ方程式の跡を取ると $D_t(2p_\perp+p_\parallel)=-(2p_\perp+p_\parallel)\Theta-2[p_\perp\Theta+(p_\parallel-p_\perp)K]$ である。平行の式を引いて二で割れば $D_t p_\perp=p_\perp(K-2\Theta)$ を得る。流速 $\vb u$ に磁束が凍結する誘導式を併用すると $D_t\ln B=K-\Theta$、連続式では $D_t\ln n=-\Theta$ となる。よって
$$
D_t\left(\frac{p_\perp}{nB}\right)=0,\qquad
D_t\left(\frac{p_\parallel B^2}{n^3}\right)=0.
$$
この組を :term[二重断熱閉包]{en="double adiabatic closure" yomi="にじゅうだんねつへいほう"} と呼ぶ。粒子種の圧力を個別に扱う場合には、磁束を運ぶ速度と、その圧力を移流する速度を同じと近似できるか確認する。

### 等方断熱と二重断熱の圧力比

一様な物質要素が磁束を保って変形し、熱流束・散乱を無視する。旋回磁気モーメントが一定なら $T_\perp\propto B$ なので $p_\perp/(nB)$ が一定である。平行方向に周期的に反射する粒子で長さ $\ell$ の変化が往復より遅ければ、平行作用が $m v_\parallel\ell$ に比例して保存する。磁束の保存から横断面は $B^{-1}$ に比例し、粒子数の保存から $\ell\propto B/n$ となるため $p_\parallel B^2/n^3$ が一定である。この要素模型は二重断熱閉包の二つの不変式を再現する。粒子のピッチ角を再分配する衝突が速い場合には、単一圧力に対して三次元理想気体の断熱式 $p\propto n^{5/3}$ を使う。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $N$ | 密度の増加倍率 | 1 |
| $M$ | 磁場の増加倍率 | 1 |

垂直圧力比を $p_\perp/p_0$ とする。磁気モーメントと粒子数から求める。

$$
p_\perp/p_0=M N
$$

平行圧力比を $p_\parallel/p_0$ とする。平行作用と磁束を使う。

$$
p_\parallel/p_0=\frac{N^{3}}{M^{2}}
$$

圧力異方性を $p_\perp/p_\parallel$ とする。初期圧力は両方向で $p_0$ とする。

$$
p_\perp/p_\parallel=\frac{p_\perp/p_0}{p_\parallel/p_0}
$$

二重断熱の内部エネルギーを $U/p_0$ とする。垂直の自由度は二つである。

$$
U/p_0=\frac{p_\parallel}{2p_0}+\frac{p_\perp}{p_0}
$$

等方断熱の圧力を $p_{\rm iso}/p_0$ とする。同じ密度変化に対する比較値である。

$$
p_{\rm iso}/p_0=N^{\frac{5}{3}}
$$

等方断熱の内部エネルギーを $U_{\rm iso}/p_0$ とする。三自由度の和を取る。

$$
U_{\rm iso}/p_0=\frac{3 p_{\rm iso}/p_0}{2}
$$

:::example[等方な圧縮]
二つの閉包が同じ圧力を与える変形を確かめる。

入力値を次のように置く。

$$
\begin{aligned}
N&=8 \\M&=4
\end{aligned}
$$

これらを順に代入する。

垂直圧力比は

$$
p_\perp/p_0=\left(4\right) \left(8\right)=32
$$

である。

平行圧力比は

$$
p_\parallel/p_0=\frac{\left(8\right)^{3}}{\left(4\right)^{2}}=32
$$

である。

圧力異方性は

$$
p_\perp/p_\parallel=1
$$

である。

二重断熱の内部エネルギーは

$$
U/p_0=\frac{3 \left(32\right)}{2}=48
$$

である。

等方断熱の圧力は

$$
p_{\rm iso}/p_0=\left(8\right)^{\frac{5}{3}}=32
$$

である。

等方断熱の内部エネルギーは

$$
U_{\rm iso}/p_0=\frac{3 \left(32\right)}{2}=48
$$

である。

横断面と平行長が同じ長さ倍率で縮むと異方性を作らない。
:::

:::example[磁場平行だけの圧縮]
磁場が不変の圧縮でどの成分が増えるか求める。

入力値を次のように置く。

$$
\begin{aligned}
N&=2 \\M&=1
\end{aligned}
$$

これらを順に代入する。

垂直圧力比は

$$
p_\perp/p_0=\left(1\right) \left(2\right)=2
$$

である。

平行圧力比は

$$
p_\parallel/p_0=\frac{\left(2\right)^{3}}{\left(1\right)^{2}}=8
$$

である。

圧力異方性は

$$
p_\perp/p_\parallel=\frac{\left(2\right)}{\left(8\right)}=\frac{1}{4}
$$

である。

二重断熱の内部エネルギーは

$$
U/p_0=\left(2\right) + \frac{\left(8\right)}{2}=6
$$

である。

等方断熱の圧力は

$$
p_{\rm iso}/p_0=\left(2\right)^{\frac{5}{3}}\simeq 3.1748
$$

である。

等方断熱の内部エネルギーは

$$
U_{\rm iso}/p_0=\frac{3 \left(2 \cdot 2^{\frac{2}{3}}\right)}{2}\simeq 4.7622
$$

である。

平行の圧力が大きく増加するので等方の閉包とは異なる。
:::

:::example[横断面だけの圧縮]
平行長を変えないときの圧力比を求める。

入力値を次のように置く。

$$
\begin{aligned}
N&=2 \\M&=2
\end{aligned}
$$

これらを順に代入する。

垂直圧力比は

$$
p_\perp/p_0=\left(2\right)^{2}=4
$$

である。

平行圧力比は

$$
p_\parallel/p_0=\left(2\right)=2
$$

である。

圧力異方性は

$$
p_\perp/p_\parallel=\frac{\left(4\right)}{\left(2\right)}=2
$$

である。

二重断熱の内部エネルギーは

$$
U/p_0=\frac{\left(2\right)}{2} + \left(4\right)=5
$$

である。

等方断熱の圧力は

$$
p_{\rm iso}/p_0=\left(2\right)^{\frac{5}{3}}\simeq 3.1748
$$

である。

等方断熱の内部エネルギーは

$$
U_{\rm iso}/p_0=\frac{3 \left(2 \cdot 2^{\frac{2}{3}}\right)}{2}\simeq 4.7622
$$

である。

平行温度は一定、横温度は磁場と同じ倍率で上がる。
:::

:::example[密度一定の磁場増加]
体積一定でも内部エネルギーが変わるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
N&=1 \\M&=2
\end{aligned}
$$

これらを順に代入する。

垂直圧力比は

$$
p_\perp/p_0=\left(1\right) \left(2\right)=2
$$

である。

平行圧力比は

$$
p_\parallel/p_0=\frac{\left(1\right)^{3}}{\left(2\right)^{2}}=\frac{1}{4}
$$

である。

圧力異方性は

$$
p_\perp/p_\parallel=\frac{\left(2\right)}{\left(\frac{1}{4}\right)}=8
$$

である。

二重断熱の内部エネルギーは

$$
U/p_0=\left(2\right) + \frac{\left(\frac{1}{4}\right)}{2}=\frac{17}{8}
$$

である。

等方断熱の圧力は

$$
p_{\rm iso}/p_0=\left(1\right)^{\frac{5}{3}}=1
$$

である。

等方断熱の内部エネルギーは

$$
U_{\rm iso}/p_0=\frac{3 \left(1\right)}{2}=\frac{3}{2}
$$

である。

垂直圧力と平行圧力が逆に変わり、等方圧力の体積仕事の式だけでは説明できない。
:::


## 一般化 Ohm の法則

電子にかかる力を電場について解くと、磁束を運ぶ速度と散逸を起こす項を分けられる。一価で準中性を使い $n_i=n_e=n$、$\vb j=en(\vb u_i-\vb u_e)$ とする。$\vb R_e=m_e n\nu_{ei}(\vb u_i-\vb u_e)$ を仮定して
$$
\vb E+\vb u_i\times\vb B
=\eta\vb j+\frac{\vb j\times\vb B}{en}
-\frac{\nabla\cdot\mathsf P_e}{en}
-\frac{m_e}{e}(\partial_t+\vb u_e\cdot\nabla)\vb u_e,\qquad
\eta=\frac{m_e\nu_{ei}}{ne^2}.
\label{eq-plasma-ohm}
$$
これを :term[一般化 Ohm の法則]{en="generalized Ohm law" yomi="いっぱんかおーむのほうそく"} という。Hall 項 $\vb j\times\vb B/(en)$ と電流の内積はゼロである。従って Hall 項は電場を変えるが直接の Joule 熱は作らない。イオン速度を質量平均速度 $\vb U$ に変える場合、$\vb u_i-\vb U=(\rho_e/\rho)\vb j/(en)$ である。従って左辺を $\vb E+\vb U\times\vb B$ にすると、Hall 項には係数 $\rho_i/\rho$ が付く。電子質量をイオン質量より十分小さくすると、この係数は一に近づく。

:::example[電子慣性の空間尺度]
一様 $n$、ほぼ静止イオンの線形問題では電子慣性項が $(m_e/(ne^2))\partial_t\vb j$ になる。Ampère 則と Faraday 則を組み合わせると、その誘導効果は $d_e^2 k^2$ が目安になる。ここで $d_e^2=m_e/(\mu_0ne^2)$ である。
:::

### 電子の磁束と圧力の回転

電子慣性と抵抗を落とすと、式 :eqref[eq-plasma-ohm] は $\vb E+\vb u_e\times\vb B=-\nabla p_e/(en)$ である。$p_e=p_e(n)$ のとき右辺は $-\nabla[\int^n p_e'(s)\dd s/(es)]$ であり、curl が消える。従って Hall 項を保っても磁束は電子流に凍結する。この関係に対し、$p_e=nT_e$ で $\nabla T_e$ と $\nabla n$ が平行でない場合は
$$
\partial_t\vb B-\nabla\times(\vb u_e\times\vb B)
=\nabla(1/(en))\times\nabla p_e
=\frac1e\nabla T_e\times\nabla\ln n.
$$
これは :term[Biermann 電池]{en="Biermann battery" yomi="びあまんでんち"} の磁場生成項である。

:::example[直交勾配]
$\nabla T_e=T'\hat{\vb x}$、$\nabla\ln n=N'\hat{\vb y}$ なら、磁場がゼロの時点でも $\partial_t B_z=T'N'/e$ である。両勾配が平行ならこの生成項はゼロになる。
:::

### Hall 項と電子慣性の波数による選別

一価で準中性な電子–イオン系で $d_s^2=m_s/(\mu_0ne^2)$ とおく。波数を $k$、時間の尺度を $1/(kv_A)$ とし、$\delta B/B_0\sim u/v_A$ を仮定する。Ampère の法則から $J/(ne)\sim k d_i u$ である。従って Hall 電場と対流電場の比は $kd_i$ の尺度になる。一方、$m_e\partial_t\vb J/(ne^2)$ と $\vb u\times\vb B$ の比は $k^2 d_e^2$ の尺度である。抵抗項の比は $\eta_m k/v_A$ となる。この順序評価は Alfvén 的な振幅比を基準にしたものであり、対流電場が消える偏波や周波数が異なる分枝では項を再評価する。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $K$ | イオン慣性長で規格化した波数 $kd_i$ | 1 |
| $R$ | 質量比 $m_i/m_e$ | 1 |
| $D$ | 磁気拡散の比 $\eta_m/(v_A d_i)$ | 1 |

Hall 電場の相対尺度を $H$ とする。イオンと電子の流速差を評価する。

$$
H=K
$$

電子慣性長で規格化した波数を $kd_e$ とする。慣性長は質量の平方根に比例する。

$$
kd_e=\frac{K}{\sqrt{R}}
$$

電流の時間微分による電場尺度を $M$ とする。Ampère の法則と Alfvén 時間を使う。

$$
M=(kd_e)^2
$$

抵抗性電場の相対尺度を $P$ とする。磁場の拡散時定数と伝播時間の比である。

$$
P=D K
$$

電子慣性と Hall の比を $M/H$ とする。対流電場による規格化を約分する。

$$
M/H=\frac{M}{H}
$$

抵抗と Hall の比を $P/H$ とする。この振幅順序では波数因子が消える。

$$
P/H=\frac{P}{H}
$$

:::example[MHD の長波長]
三つの補正が同時に小さいか求める。

入力値を次のように置く。

$$
\begin{aligned}
K&=\frac{1}{100} \\R&=1836 \\D&=\frac{1}{10000}
\end{aligned}
$$

これらを順に代入する。

Hall 電場の相対尺度は

$$
H=\left(\frac{1}{100}\right)=\frac{1}{100}
$$

である。

電子慣性長で規格化した波数は

$$
kd_e=\frac{\left(\frac{1}{100}\right)}{\sqrt{\left(1836\right)}}\simeq 0.00023338
$$

である。

電流の時間微分による電場尺度は

$$
M=\left(\frac{\sqrt{51}}{30600}\right)^{2}=\frac{1}{18360000}
$$

である。

抵抗性電場の相対尺度は

$$
P=\left(\frac{1}{10000}\right) \left(\frac{1}{100}\right)=\frac{1}{1000000}
$$

である。

電子慣性と Hall の比は

$$
M/H=\frac{\left(\frac{1}{18360000}\right)}{\left(\frac{1}{100}\right)}=\frac{1}{183600}
$$

である。

抵抗と Hall の比は

$$
P/H=\frac{\left(\frac{1}{1000000}\right)}{\left(\frac{1}{100}\right)}=\frac{1}{10000}
$$

である。

この比較で最大の Hall 項も一パーセントの尺度である。
:::

:::example[イオン慣性長]
Hall だけを残す順序を調べる。

入力値を次のように置く。

$$
\begin{aligned}
K&=1 \\R&=1836 \\D&=\frac{1}{10000}
\end{aligned}
$$

これらを順に代入する。

Hall 電場の相対尺度は

$$
H=\left(1\right)=1
$$

である。

電子慣性長で規格化した波数は

$$
kd_e=\frac{\left(1\right)}{\sqrt{\left(1836\right)}}\simeq 0.023338
$$

である。

電流の時間微分による電場尺度は

$$
M=\left(\frac{\sqrt{51}}{306}\right)^{2}=\frac{1}{1836}
$$

である。

抵抗性電場の相対尺度は

$$
P=\left(1\right) \left(\frac{1}{10000}\right)=\frac{1}{10000}
$$

である。

電子慣性と Hall の比は

$$
M/H=\frac{\left(\frac{1}{1836}\right)}{\left(1\right)}=\frac{1}{1836}
$$

である。

抵抗と Hall の比は

$$
P/H=\frac{\left(\frac{1}{10000}\right)}{\left(1\right)}=\frac{1}{10000}
$$

である。

電子慣性と抵抗の比はまだ小さい。
:::

:::example[電子慣性長]
$kd_e=1$ の場合の相対項を求める。

入力値を次のように置く。

$$
\begin{aligned}
K&=6 \sqrt{51} \\R&=1836 \\D&=\frac{1}{10000}
\end{aligned}
$$

これらを順に代入する。

Hall 電場の相対尺度は

$$
H=\left(6 \sqrt{51}\right)\simeq 42.849
$$

である。

電子慣性長で規格化した波数は

$$
kd_e=\frac{\left(6 \sqrt{51}\right)}{\sqrt{\left(1836\right)}}=1
$$

である。

電流の時間微分による電場尺度は

$$
M=\left(1\right)^{2}=1
$$

である。

抵抗性電場の相対尺度は

$$
P=\left(6 \sqrt{51}\right) \left(\frac{1}{10000}\right)\simeq 0.0042849
$$

である。

電子慣性と Hall の比は

$$
M/H=\frac{\left(1\right)}{\left(6 \sqrt{51}\right)}\simeq 0.023338
$$

である。

抵抗と Hall の比は

$$
P/H=\frac{\left(\frac{3 \sqrt{51}}{5000}\right)}{\left(6 \sqrt{51}\right)}=\frac{1}{10000}
$$

である。

この波数では Alfvén 的な基準から外れる分枝があり、実際には分散と振幅を解く必要がある。
:::

:::example[高抵抗の流体]
長波長でも抵抗が無視できない例を計算する。

入力値を次のように置く。

$$
\begin{aligned}
K&=\frac{1}{10} \\R&=1836 \\D&=5
\end{aligned}
$$

これらを順に代入する。

Hall 電場の相対尺度は

$$
H=\left(\frac{1}{10}\right)=\frac{1}{10}
$$

である。

電子慣性長で規格化した波数は

$$
kd_e=\frac{\left(\frac{1}{10}\right)}{\sqrt{\left(1836\right)}}\simeq 0.0023338
$$

である。

電流の時間微分による電場尺度は

$$
M=\left(\frac{\sqrt{51}}{3060}\right)^{2}=\frac{1}{183600}
$$

である。

抵抗性電場の相対尺度は

$$
P=\left(5\right) \left(\frac{1}{10}\right)=\frac{1}{2}
$$

である。

電子慣性と Hall の比は

$$
M/H=\frac{\left(\frac{1}{183600}\right)}{\left(\frac{1}{10}\right)}=\frac{1}{18360}
$$

である。

抵抗と Hall の比は

$$
P/H=\frac{\left(\frac{1}{2}\right)}{\left(\frac{1}{10}\right)}=5
$$

である。

質量比の小ささは抵抗項を抑えず、電子慣性だけが特に小さい。
:::


### Biermann 電池による磁場生成の符号と向き

電子慣性と摩擦を省き、準中性で電子圧を等方とする。一般化 Ohm 則を電子流で書くと $\vb E+\vb u_e\times\vb B=-\grad p_e/(en)$ である。Faraday の法則に代入すると
$$
\partial_t\vb B=\rot(\vb u_e\times\vb B)
+\grad\left(\frac1{en}\right)\times\grad p_e
=\rot(\vb u_e\times\vb B)+\frac1e\grad T_e\times\grad\ln n.
$$
最後の等式では $p_e=nT_e$ を使い、$\grad n$ 同士の外積がゼロになることを用いた。これを :term[Biermann バッテリー]{en="Biermann battery" yomi="びあまんばってりー"} という初期磁場源と呼ぶ。符号は Faraday の法則の負符号と $1/n$ の微分に由来する。等温なら源が消えるが、等温でなくても $\grad T_e\parallel\grad n$ なら消える。微小時間に限って背景勾配を一定とし $\vb B(0)=0$ と置くと $\vb B(t)=t\grad(T_e/e)\times\grad\ln n+O(t^2)$ である。単位は $(\mathrm{V\,m^{-1}})(\mathrm{m^{-1}})=\mathrm{T\,s^{-1}}$ となる。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $G_{T,x}$ | $T_e/e$ の $x$ 方向勾配 | $\mathrm{V\,m^{-1}}$ |
| $G_{T,y}$ | $T_e/e$ の $y$ 方向勾配 | $\mathrm{V\,m^{-1}}$ |
| $G_{n,x}$ | $\ln n$ の $x$ 方向勾配 | $\mathrm{m^{-1}}$ |
| $G_{n,y}$ | $\ln n$ の $y$ 方向勾配 | $\mathrm{m^{-1}}$ |
| $t$ | 経過時間 | $\mathrm s$ |

温度の $x$ 勾配と密度の $y$ 勾配の積を $S_{xy}$ とする。右手系で正の $z$ 寄与である。

$$
S_{xy}=G_{T,x} G_{n,y}
$$

温度の $y$ 勾配と密度の $x$ 勾配の符号付き積を $S_{yx}$ とする。外積の順を保つ。

$$
S_{yx}=- G_{T,y} G_{n,x}
$$

初期誘導率を $\dot B_z(0)$ とする。二つの寄与を足す。

$$
\dot B_z(0)=S_{xy} + S_{yx}
$$

短時間で生じる磁場を $B_z^{(1)}$ とする。移流と背景の変化を高次とする。

$$
B_z^{(1)}=\dot B_z(0) t
$$

温度勾配の大きさを $\abs{\grad(T_e/e)}$ とする。外積の方向とは別に評価する。

$$
\abs{\grad(T_e/e)}=\sqrt{G_{T,x}^{2} + G_{T,y}^{2}}
$$

勾配がなす有向角の正弦を $\sin\theta$ とする。二勾配がゼロでない場合に求める。

$$
\sin\theta=\frac{\dot B_z(0)}{\abs{\grad(T_e/e)}\abs{\grad\ln n}}
$$

:::example[直交する勾配]
誘導の初期方向と大きさを求める。

入力値を次のように置く。

$$
\begin{aligned}
G_{T,x}&=1\,\mathrm{V\,m^{-1}} \\G_{T,y}&=0\,\mathrm{V\,m^{-1}} \\G_{n,x}&=0\,\mathrm{m^{-1}} \\G_{n,y}&=1\,\mathrm{m^{-1}} \\t&=\frac{1}{1000}\,\mathrm s
\end{aligned}
$$

これらを順に代入する。

温度の $x$ 勾配と密度の $y$ 勾配の積は

$$
S_{xy}=\left(1\right)^{2}=1\,\mathrm{T\,s^{-1}}
$$

である。

温度の $y$ 勾配と密度の $x$ 勾配の符号付き積は

$$
S_{yx}=- \left(0\right)^{2}=0\,\mathrm{T\,s^{-1}}
$$

である。

初期誘導率は

$$
\dot B_z(0)=\left(0\right) + \left(1\right)=1\,\mathrm{T\,s^{-1}}
$$

である。

短時間で生じる磁場は

$$
B_z^{(1)}=\left(1\right) \left(\frac{1}{1000}\right)=\frac{1}{1000}\,\mathrm T
$$

である。

温度勾配の大きさは

$$
\abs{\grad(T_e/e)}=\sqrt{\left(0\right)^{2} + \left(1\right)^{2}}=1\,\mathrm{V\,m^{-1}}
$$

である。

勾配がなす有向角の正弦は

$$
\sin\theta=\frac{1}{\sqrt{\left(0\right)^{2} + \left(1\right)^{2}}}=1
$$

である。

正の $z$ 方向に磁場が生じる。
:::

:::example[直交の交換]
勾配の順を交換したときに符号を確認する。

入力値を次のように置く。

$$
\begin{aligned}
G_{T,x}&=0\,\mathrm{V\,m^{-1}} \\G_{T,y}&=1\,\mathrm{V\,m^{-1}} \\G_{n,x}&=1\,\mathrm{m^{-1}} \\G_{n,y}&=0\,\mathrm{m^{-1}} \\t&=\frac{1}{1000}\,\mathrm s
\end{aligned}
$$

これらを順に代入する。

温度の $x$ 勾配と密度の $y$ 勾配の積は

$$
S_{xy}=\left(0\right)^{2}=0\,\mathrm{T\,s^{-1}}
$$

である。

温度の $y$ 勾配と密度の $x$ 勾配の符号付き積は

$$
S_{yx}=- \left(1\right)^{2}=-1\,\mathrm{T\,s^{-1}}
$$

である。

初期誘導率は

$$
\dot B_z(0)=\left(-1\right) + \left(0\right)=-1\,\mathrm{T\,s^{-1}}
$$

である。

短時間で生じる磁場は

$$
B_z^{(1)}=\left(-1\right) \left(\frac{1}{1000}\right)=- \frac{1}{1000}\,\mathrm T
$$

である。

温度勾配の大きさは

$$
\abs{\grad(T_e/e)}=\sqrt{\left(0\right)^{2} + \left(1\right)^{2}}=1\,\mathrm{V\,m^{-1}}
$$

である。

勾配がなす有向角の正弦は

$$
\sin\theta=\frac{\left(-1\right)}{\left(1\right) \sqrt{\left(0\right)^{2} + \left(1\right)^{2}}}=-1
$$

である。

二勾配の大きさは同じでも磁場は逆向きになる。
:::

:::example[同一直線上の勾配]
両方の勾配があっても初期磁場が生じるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
G_{T,x}&=1\,\mathrm{V\,m^{-1}} \\G_{T,y}&=2\,\mathrm{V\,m^{-1}} \\G_{n,x}&=2\,\mathrm{m^{-1}} \\G_{n,y}&=4\,\mathrm{m^{-1}} \\t&=\frac{1}{1000}\,\mathrm s
\end{aligned}
$$

これらを順に代入する。

温度の $x$ 勾配と密度の $y$ 勾配の積は

$$
S_{xy}=\left(1\right) \left(4\right)=4\,\mathrm{T\,s^{-1}}
$$

である。

温度の $y$ 勾配と密度の $x$ 勾配の符号付き積は

$$
S_{yx}=- \left(2\right)^{2}=-4\,\mathrm{T\,s^{-1}}
$$

である。

初期誘導率は

$$
\dot B_z(0)=\left(-4\right) + \left(4\right)=0\,\mathrm{T\,s^{-1}}
$$

である。

短時間で生じる磁場は

$$
B_z^{(1)}=\left(0\right) \left(\frac{1}{1000}\right)=0\,\mathrm T
$$

である。

温度勾配の大きさは

$$
\abs{\grad(T_e/e)}=\sqrt{\left(1\right)^{2} + \left(2\right)^{2}}\simeq 2.2361\,\mathrm{V\,m^{-1}}
$$

である。

勾配がなす有向角の正弦は

$$
\sin\theta=\frac{\left(0\right)}{\left(\sqrt{5}\right) \sqrt{\left(2\right)^{2} + \left(4\right)^{2}}}=0
$$

である。

外積がゼロで、圧力電場の回転は消える。
:::

:::example[斜めの弱い源]
勾配の大きさと源の小ささが両立するか調べる。

入力値を次のように置く。

$$
\begin{aligned}
G_{T,x}&=1\,\mathrm{V\,m^{-1}} \\G_{T,y}&=1\,\mathrm{V\,m^{-1}} \\G_{n,x}&=1\,\mathrm{m^{-1}} \\G_{n,y}&=\frac{11}{10}\,\mathrm{m^{-1}} \\t&=\frac{1}{1000}\,\mathrm s
\end{aligned}
$$

これらを順に代入する。

温度の $x$ 勾配と密度の $y$ 勾配の積は

$$
S_{xy}=\left(1\right) \left(\frac{11}{10}\right)=\frac{11}{10}\,\mathrm{T\,s^{-1}}
$$

である。

温度の $y$ 勾配と密度の $x$ 勾配の符号付き積は

$$
S_{yx}=- \left(1\right)^{2}=-1\,\mathrm{T\,s^{-1}}
$$

である。

初期誘導率は

$$
\dot B_z(0)=\left(-1\right) + \left(\frac{11}{10}\right)=\frac{1}{10}\,\mathrm{T\,s^{-1}}
$$

である。

短時間で生じる磁場は

$$
B_z^{(1)}=\left(\frac{1}{1000}\right) \left(\frac{1}{10}\right)=\frac{1}{10000}\,\mathrm T
$$

である。

温度勾配の大きさは

$$
\abs{\grad(T_e/e)}=\sqrt{2} \sqrt{\left(1\right)^{2}}\simeq 1.4142\,\mathrm{V\,m^{-1}}
$$

である。

勾配がなす有向角の正弦は

$$
\sin\theta=\frac{\left(\frac{1}{10}\right)}{\left(\sqrt{2}\right) \sqrt{\left(1\right)^{2} + \left(\frac{11}{10}\right)^{2}}}\simeq 0.047565
$$

である。

ほぼ平行なので、各寄与に対して合計は小さい。
:::


## 一流体方程式への縮約

両流体の運動量を足すときは、質量平均を先にとる。$\rho=\rho_i+\rho_e$、$\vb U=(\rho_i\vb u_i+\rho_e\vb u_e)/\rho$、$\rho_s=m_s n_s$ とする。質量流束は正確に $\rho\vb U$ である。速度差 $\vb w=\vb u_i-\vb u_e$ で書けば
$$
\sum_s\rho_s\vb u_s\vb u_s
=\rho\vb U\vb U+
\frac{\rho_i\rho_e}{\rho}\vb w\vb w
$$
である。実際 $\vb u_i-\vb U=(\rho_e/\rho)\vb w$、$\vb u_e-\vb U=-(\rho_i/\rho)\vb w$ を代入すると一次の交差項が消える。各種の連続式を足すと $\partial_t\rho+\nabla\cdot(\rho\vb U)=0$ になる。運動量式の移流項を保存形にしてから足し、内部摩擦を相殺すると
$$
\partial_t(\rho\vb U)+\nabla\cdot\left(
\rho\vb U\vb U+\frac{\rho_i\rho_e}{\rho}\vb w\vb w+\mathsf P_i+\mathsf P_e\right)
=\rho_q\vb E+\vb j\times\vb B .
$$
相対流による運動量流束を省く場合は、他の流束や力に比べてその発散が小さいことを確認する。

### MHD の尺度

準中性、等方圧、相対流の応力が小さい条件では $\rho D\vb U/Dt=\vb j\times\vb B-\nabla p$ となる。低周波で変位電流を無視し、$\vb j\sim B/(\mu_0L)$ と見積もる。$U$ が Alfvén 速度 $v_A=B/\sqrt{\mu_0\rho}$ 程度なら Hall 項と $U B$ の比は $d_i/L$、$d_i=\sqrt{m_i/(\mu_0 ne^2)}$ である。$d_i/L\ll1$ は Hall 項を省く判断の一つであり、電子圧力や抵抗は別に見積もる。

:::example[質量平均と電流]
$m_i/m_e=2000$、$\vb u_i=0$、$\vb u_e=2001\,\hat{\vb x}\,\mathrm{m\,s^{-1}}$ なら $\vb U=1\,\hat{\vb x}\,\mathrm{m\,s^{-1}}$ だが $\vb j=-2001en\hat{\vb x}$ である。質量流の小ささだけでは磁気力の小ささは決まらない。
:::

:::example[二流体の残る領域]
代表速度が $v_A$ の一価イオン系で $L=d_i/2$ なら Hall 項は移流電場の二倍程度と見積もられる。一方 $d_e/L=2\sqrt{m_e/m_i}$ は陽子で約 $0.047$ である。準中性が成り立つことを別に確認すれば、イオンと電子の流れを分け、電子慣性だけを最初に無視する根拠になる。
:::

### 質量流と電流のエネルギーを分離する

一価イオンを仮定し $\rho_s=m_s n$ とする。$\vb u=(\rho_i\vb u_i+\rho_e\vb u_e)/\rho$、$\vb w=\vb j/(en)$ と書けば $\vb u_i=\vb u+(\rho_e/\rho)\vb w$、$\vb u_e=\vb u-(\rho_i/\rho)\vb w$ である。速度積の和を展開すると混合項は消え、$\sum_s\rho_s u_{si}u_{sj}=\rho u_i u_j+(\rho_e\rho_i/\rho)w_iw_j$ となる。最後の項は圧力テンソルと同じ次元の運動量流束である。ただし生じる理由は種の相対流であり、各種の中心化された圧力へ含めない。エネルギーでは跡を半分にした値を使う。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $r$ | 質量比 $m_i/m_e$ | 1 |
| $u$ | 一方向の質量平均速度 | $\mathrm{m\,s^{-1}}$ |
| $w$ | 同じ方向にとった電流の速度 | $\mathrm{m\,s^{-1}}$ |

イオンの速度を $u_i$ とする。イオンの偏差は電子の質量分率で決まる。

$$
u_i=u + \frac{w}{r + 1}
$$

電子の速度を $u_e$ とする。イオンの偏差と符号が逆である。

$$
u_e=- \frac{r w}{r + 1} + u
$$

質量流の単位質量当たりエネルギーを $K_{\rm bulk}/\rho$ とする。一流体の運動に対応する。

$$
K_{\rm bulk}/\rho=\frac{u^{2}}{2}
$$

電流の単位質量当たりエネルギーを $K_{\rm rel}/\rho$ とする。電流の速度差の二乗を使う。

$$
K_{\rm rel}/\rho=\frac{r w^{2}}{2 \left(r + 1\right)^{2}}
$$

種の平均運動の総和を $K_{\rm tot}/\rho$ とする。熱運動を含める場合には圧力の跡を別に足す。

$$
K_{\rm tot}/\rho=K_{\rm bulk}/\rho + K_{\rm rel}/\rho
$$

相対流が担う割合を $K_{\rm rel}/K_{\rm tot}$ とする。電子に小さな質量しかなくても速度差が大きい場合には無視できない。

$$
K_{\rm rel}/K_{\rm tot}=\frac{K_{\rm rel}/\rho}{K_{\rm tot}/\rho}
$$

:::example[電流のない質量流]
二種が同じ速度のときの縮約を確かめる。

入力値を次のように置く。

$$
\begin{aligned}
r&=1836 \\u&=100\,\mathrm{m\,s^{-1}} \\w&=0\,\mathrm{m\,s^{-1}}
\end{aligned}
$$

これらを順に代入する。

イオンの速度は

$$
u_i=\frac{\left(0\right)}{\left(1836\right) + 1} + \left(100\right)=100\,\mathrm{m\,s^{-1}}
$$

である。

電子の速度は

$$
u_e=- \frac{\left(0\right) \left(1836\right)}{\left(1836\right) + 1} + \left(100\right)=100\,\mathrm{m\,s^{-1}}
$$

である。

質量流の単位質量当たりエネルギーは

$$
K_{\rm bulk}/\rho=\frac{\left(100\right)^{2}}{2}=5000\,\mathrm{m^2\,s^{-2}}
$$

である。

電流の単位質量当たりエネルギーは

$$
K_{\rm rel}/\rho=\frac{\left(0\right)^{2} \left(1836\right)}{2 \left(\left(1836\right) + 1\right)^{2}}=0\,\mathrm{m^2\,s^{-2}}
$$

である。

種の平均運動の総和は

$$
K_{\rm tot}/\rho=\left(0\right) + \left(5000\right)=5000\,\mathrm{m^2\,s^{-2}}
$$

である。

相対流が担う割合は

$$
K_{\rm rel}/K_{\rm tot}=\frac{\left(0\right)}{\left(5000\right)}=0
$$

である。

相対項が厳密に消える。
:::

:::example[電流だけが流れる配置]
質量流がゼロでも二種の運動エネルギーが残るか調べる。

入力値を次のように置く。

$$
\begin{aligned}
r&=1836 \\u&=0\,\mathrm{m\,s^{-1}} \\w&=100\,\mathrm{m\,s^{-1}}
\end{aligned}
$$

これらを順に代入する。

イオンの速度は

$$
u_i=\left(0\right) + \frac{\left(100\right)}{\left(1836\right) + 1}=\frac{100}{1837}\,\mathrm{m\,s^{-1}}
$$

である。

電子の速度は

$$
u_e=\left(0\right) - \frac{\left(100\right) \left(1836\right)}{\left(1836\right) + 1}=- \frac{183600}{1837}\,\mathrm{m\,s^{-1}}
$$

である。

質量流の単位質量当たりエネルギーは

$$
K_{\rm bulk}/\rho=\frac{\left(0\right)^{2}}{2}=0\,\mathrm{m^2\,s^{-2}}
$$

である。

電流の単位質量当たりエネルギーは

$$
K_{\rm rel}/\rho=\frac{\left(100\right)^{2} \left(1836\right)}{2 \left(\left(1836\right) + 1\right)^{2}}\simeq 2.7203\,\mathrm{m^2\,s^{-2}}
$$

である。

種の平均運動の総和は

$$
K_{\rm tot}/\rho=\left(0\right) + \left(\frac{9180000}{3374569}\right)\simeq 2.7203\,\mathrm{m^2\,s^{-2}}
$$

である。

相対流が担う割合は

$$
K_{\rm rel}/K_{\rm tot}=1
$$

である。

残るエネルギーは荷電粒子の相対運動だけである。
:::

:::example[電子が大きく流れる配置]
電子質量が小さくても相対項が重要な例を求める。

入力値を次のように置く。

$$
\begin{aligned}
r&=1836 \\u&=10\,\mathrm{m\,s^{-1}} \\w&=1000\,\mathrm{m\,s^{-1}}
\end{aligned}
$$

これらを順に代入する。

イオンの速度は

$$
u_i=\frac{\left(1000\right)}{\left(1836\right) + 1} + \left(10\right)=\frac{19370}{1837}\,\mathrm{m\,s^{-1}}
$$

である。

電子の速度は

$$
u_e=- \frac{\left(1000\right) \left(1836\right)}{\left(1836\right) + 1} + \left(10\right)=- \frac{1817630}{1837}\,\mathrm{m\,s^{-1}}
$$

である。

質量流の単位質量当たりエネルギーは

$$
K_{\rm bulk}/\rho=\frac{\left(10\right)^{2}}{2}=50\,\mathrm{m^2\,s^{-2}}
$$

である。

電流の単位質量当たりエネルギーは

$$
K_{\rm rel}/\rho=\frac{\left(1000\right)^{2} \left(1836\right)}{2 \left(\left(1836\right) + 1\right)^{2}}\simeq 272.03\,\mathrm{m^2\,s^{-2}}
$$

である。

種の平均運動の総和は

$$
K_{\rm tot}/\rho=\left(50\right) + \left(\frac{918000000}{3374569}\right)\simeq 322.03\,\mathrm{m^2\,s^{-2}}
$$

である。

相対流が担う割合は

$$
K_{\rm rel}/K_{\rm tot}=\frac{\left(\frac{918000000}{3374569}\right)}{\left(\frac{1086728450}{3374569}\right)}\simeq 0.84474
$$

である。

速度差が質量流の百倍なので、相対項の係数が小さいことだけでは落とせない。
:::

:::example[等質量の場合]
対プラズマの速度とエネルギーを求める。

入力値を次のように置く。

$$
\begin{aligned}
r&=1 \\u&=10\,\mathrm{m\,s^{-1}} \\w&=20\,\mathrm{m\,s^{-1}}
\end{aligned}
$$

これらを順に代入する。

イオンの速度は

$$
u_i=\left(10\right) + \frac{\left(20\right)}{\left(1\right) + 1}=20\,\mathrm{m\,s^{-1}}
$$

である。

電子の速度は

$$
u_e=\left(10\right) - \frac{\left(1\right) \left(20\right)}{\left(1\right) + 1}=0\,\mathrm{m\,s^{-1}}
$$

である。

質量流の単位質量当たりエネルギーは

$$
K_{\rm bulk}/\rho=\frac{\left(10\right)^{2}}{2}=50\,\mathrm{m^2\,s^{-2}}
$$

である。

電流の単位質量当たりエネルギーは

$$
K_{\rm rel}/\rho=\frac{\left(1\right) \left(20\right)^{2}}{2 \left(\left(1\right) + 1\right)^{2}}=50\,\mathrm{m^2\,s^{-2}}
$$

である。

種の平均運動の総和は

$$
K_{\rm tot}/\rho=2 \left(50\right)=100\,\mathrm{m^2\,s^{-2}}
$$

である。

相対流が担う割合は

$$
K_{\rm rel}/K_{\rm tot}=\frac{\left(50\right)}{\left(100\right)}=\frac{1}{2}
$$

である。

電子だけを質量なしとした縮約が適用できない。
:::


## 演習問題

::::exercise{level=basic}
等方 Maxwell 分布で $P_{xx}$ と $P_{xy}$ を求めよ。

:::solution
ガウス積分で $\langle c_x^2\rangle=T/m$、奇関数積分で $\langle c_xc_y\rangle=0$ なので $P_{xx}=nT$、$P_{xy}=0$ である。
:::
::::

::::exercise{level=basic}
$n_i=n_e=n$、$u_i=1$、$u_e=-2$ の一次元電流を求めよ。

:::solution
$j=en(u_i-u_e)=3en$ である。速度は同じ単位でそろえる。
:::
::::

::::exercise{level=basic}
断熱で密度が二倍なら $T$ の倍率を求めよ。

:::solution
$p\propto n^{5/3}$ と $T=p/n$ から $2^{2/3}$ 倍である。
:::
::::

::::exercise{level=basic}
一様磁場で $\nabla p_e=0$ で、横慣性と摩擦を無視できるとき電子の横速度は何か。

:::solution
式 :eqref[eq-plasma-diamagnetic] の反磁性項がゼロなので $\vb E\times\vb B/B^2$ である。
:::
::::

::::exercise{level=advanced}
電子が静止し、衝突を無視できるとき、等温圧力勾配を支える電場を求めよ。

:::solution
電子運動量式で $0=-en\vb E-T_e\nabla n$ だから $\vb E=-(T_e/e)\nabla\ln n$ である。
:::
::::

::::exercise{level=advanced}
種の摩擦で失われる流れのエネルギー率を求めよ。

:::solution
$\vb u_e\cdot\vb R_e+\vb u_i\cdot(-\vb R_e)=-(\vb u_i-\vb u_e)\cdot\vb R_e=-\eta j^2$ である。内部衝突ならこの分が熱に移る。
:::
::::

::::exercise{level=advanced}
$B$ と $n$ がともに二倍の要素について二重断熱の圧力倍率を求めよ。

:::solution
$p_\perp\propto nB$ より四倍、$p_\parallel\propto n^3/B^2$ より二倍である。
:::
::::

::::exercise{level=advanced}
等質量、等密度で $\vb u_i=-\vb u_e=u\hat{\vb x}$ の追加応力を求めよ。

:::solution
$\vb U=0$ であり、$\rho_i\rho_e/\rho=\rho/4$ と $\vb w=2u\hat{\vb x}$ から $\rho u^2\hat{\vb x}\hat{\vb x}$ である。
:::
::::

## 確認テスト

::::quiz
Hall 項の Joule 加熱への直接寄与はどれか。

- [x] ゼロ
- [ ] 常に正
- [ ] 常に負

:::explanation
$\vb j\cdot(\vb j\times\vb B)=0$ は外積が両因子に垂直であることから従う。
:::
::::

::::quiz{answer="32" tol="1e-6"}
断熱で密度が八倍のとき圧力の倍率を答えよ。

:::explanation
$p/n^{5/3}$ が一定なので $8^{5/3}=32$ である。
:::
::::

::::quiz
等温スカラー圧力の式を選べ。

- [x] $p=nT$
- [ ] $p=T/n$
- [ ] $p=n/T$

:::explanation
速度分布の一成分の分散は $T/m$ であり $P_{xx}=mn(T/m)=nT$ となる。
:::
::::

::::quiz{answer="0.01" tol="1e-6"}
$d_i=0.02\,\mathrm m$、$L=2\,\mathrm m$、$U=v_A$ とすると Hall 項の比の目安はいくらか。

:::explanation
$d_i/L=0.02/2=0.01$ である。
:::
::::

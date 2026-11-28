---
title: プラズマの基本的性質
description: 遮蔽と集団振動から、プラズマの空間尺度と時間尺度を導く。
notation:
  - symbol: '\lambda_D'
    meaning: Debye 長
  - symbol: '\omega_{ps}'
    meaning: 粒子種 s のプラズマ角周波数
  - symbol: 'T_s'
    meaning: 粒子種 s の温度を表すエネルギー
  - symbol: 'N_D'
    meaning: Debye 球内の電子数
---

## この章の目標

- 電離度と集団性を区別する。
- 静電遮蔽と集団振動を導く。
- 準中性近似を使える尺度を見積もる。

## プラズマの定義と存在領域

中性気体の一部を電離すると、電子とイオンは外から加えた電場で移動する。その移動で電荷が偏ると、粒子が作る電場も変わる。個々の衝突に加えて、この場を通じた応答を追う必要がある。

:::definition[集団応答]{#def-plasma-collective}
自由な荷電粒子を含み、粒子が作る電磁場によって複数の粒子が協調して応答する媒質を :term[プラズマ]{en="plasma" yomi="ぷらずま"} という。電子数密度を $n_e$、電荷数 $Z_s$ のイオンの数密度を $n_s$ とすると、電荷密度は $\rho_q=e(\sum_s Z_sn_s-n_e)$ である。ここで $e>0$ は電気素量である。
:::

この章は古典的で非相対論的な粒子を扱う。数密度の単位は $\mathrm{m^{-3}}$、温度 $T_s$ は $k_BT_{s,\mathrm K}$ を表すエネルギーであり、$1\,\mathrm{eV}=1.602176634\times10^{-19}\,\mathrm J$ である。種ごとの温度を使うときは、各種の速度分布が等方的であるかも確かめる。非等方分布には方向別の温度が必要になる。

### 電離度と熱平衡

一価イオンと中性原子だけがあるとき、電離度は $x=n_i/(n_i+n_n)$ である。熱平衡での電離を概算するには [大正準分布](/physics/statistical-mechanics/grand-canonical/) を使う。希薄理想気体では一粒子分配関数が $Vg_s e^{-E_{s0}/T}/\ell_s^3$、$\ell_s=h/\sqrt{2\pi m_sT}$ である。基底エネルギーを中性原子で $-I$、イオンと電子でゼロにとれば、化学ポテンシャルは $\mu_s=E_{s0}+T\ln(n_s\ell_s^3/g_s)$ となる。電離を一回進めたときの自由エネルギー変化 $\mu_e+\mu_i-\mu_n$ がゼロなので、
$$
\frac{n_en_i}{n_n}=\frac{g_eg_i}{g_n}\frac{\ell_n^3}{\ell_e^3\ell_i^3}e^{-I/T}.
$$
これは :term[Saha の電離式]{en="Saha ionization equation" yomi="さはのでんりしき"} である。照射や放電で電離が維持される場合、反応の収支で電離度を決める。

:::example[同じ電離度をもつ二つの気体]
$n_i=n_e=10^{16}\,\mathrm{m^{-3}}$、$n_n=9\times10^{16}\,\mathrm{m^{-3}}$ なら $x=0.1$ である。全密度を百倍にしても $x$ は $0.1$ だが、同じ温度なら後で導く Debye 長は十分の一となる。
:::

:::example[電離式から電離度を求める]
核の全密度を $N=n_n+n_i$、電離式の右辺を $S$ とすると、$n_e=n_i$ より $x^2/(1-x)=S/N$ である。$S/N=1/2$ なら $2x^2+x-1=(2x-1)(x+1)=0$ なので物理的な解は $x=1/2$ である。
:::

### 熱平衡にない希薄気体の電離

発生した光が再吸収されずに逃げると、物質の衝突と放射をまとめた熱平衡を仮定できない。例えば希薄で光学的に薄い水素気体を考え、電子衝突による電離と放射再結合だけを残す。電子と原子の電離の速度係数を $K_{\rm ion}$、電子とイオンの再結合の速度係数を $K_{\rm rec}$ とする。両係数の単位は $\mathrm{m^3\,s^{-1}}$ である。核密度 $N=n_i+n_n$ が一定で、$n_e=n_i$ なら
$$
\dv{x}{t}=N x\left[K_{\rm ion}(1-x)-K_{\rm rec}x\right],\qquad x=n_i/N .
$$
これは体積当たりの収支 $\partial_t n_i=n_e n_n K_{\rm ion}-n_e n_i K_{\rm rec}$ を $N$ で割った式である。少量の電子が存在し、係数が正で一定なら、$x=K_{\rm ion}/(K_{\rm ion}+K_{\rm rec})$ に近づく。実際、右辺はこの値より小さい正の $x$ で正、大きい $x$ では負である。この近似を :term[コロナ平衡]{en="coronal equilibrium" yomi="ころなへいこう"} と呼ぶ。たとえば $K_{\rm ion}/K_{\rm rec}=2$ なら定常電離度は $2/3$ である。この模型では核密度を増やすと平衡値は同じまま緩和が速くなる。三体再結合や放射の再吸収まで効く条件では、反応収支を追加して電離度を求める。

### 一様な交流電場での集団性

粒子と中性気体との運動量緩和率を $\nu_{en}$、交流の角周波数を $\omega$ とする。後で導く電子プラズマ角周波数を $\omega_{pe}$ とする。イオン電流を無視し、電子の一様平均流が $m_e(\partial_t+\nu_{en})\vb u_e=-e\vb E$ を満たすとき、$e^{-i\omega t}$ の振幅に対して $\vb j=\sigma\vb E$、$\sigma=n_e e^2/[m_e(\nu_{en}-i\omega)]$ である。Maxwell 方程式でこの電流と変位電流の大きさを比べると
$$
\frac{\abs{\sigma}}{\varepsilon_0\omega}
=\frac{\omega_{pe}^2}{\omega\sqrt{\omega^2+\nu_{en}^2}}.
$$
小さい電離度でも、この比が大きければ電子は場の応答に大きく寄与する。逆に周波数が高くなると、同じ媒質でも電子の速度変化が間に合わない。この比較は一様な横応答の模型であり、磁場による異方性や電離・再結合による密度の変化は省いた。幾何学的に境界の影響を受けるかは、後で求める Debye 長と試料の寸法で別に判定する。

:::example[電離度と交流応答]
$\omega_{pe}=10^8\,\mathrm{s^{-1}}$、$\nu_{en}=10^7\,\mathrm{s^{-1}}$、$\omega=10^7\,\mathrm{s^{-1}}$ なら上の比は $100/\sqrt2$ である。$\omega=10^9\,\mathrm{s^{-1}}$ に上げると比は $10^{-2}/\sqrt{1+10^{-4}}$ となる。同じ密度でも低い周波数で電子電流が強く寄与する。
:::

### 熱平衡の電離比と密度の測定

低密度の水素気体で、原子の基底状態と完全電離した状態だけを考える。全水素核密度は $n_{\rm H}=n_a+n_i$、電離度は $x=n_i/n_{\rm H}$ である。非縮退かつ熱平衡で $n_e=n_i$ とする。イオン化エネルギーを $\chi$、内部状態の統計重率を $g_s$、熱的 de Broglie 波長を $\lambda_s=h/\sqrt{2\pi m_s T}$ とすると、理想混合気体の化学ポテンシャルは $\mu_s=E_{0s}+T\ln(n_s\lambda_s^3/g_s)$ である。この式は一粒子分配関数 $Vg_s e^{-E_{0s}/T}/\lambda_s^3$ を用い、$N_s!$ を含む分配関数の対数を $N_s$ で微分すれば得られる。

反応 $a\rightleftharpoons i+e$ で全粒子数が一つ増える。そのため温度が一定でも、希釈すると電離側の並進エントロピーの寄与が大きくなる。実際、$E_{0i}+E_{0e}-E_{0a}=\chi$ と $\mu_a=\mu_i+\mu_e$ を代入すると
$$
\frac{n_i n_e}{n_a}
=\frac{g_i g_e}{g_a}\frac{\lambda_a^3}{\lambda_i^3\lambda_e^3}e^{-\chi/T}
=:S(T),\qquad
\frac{x^2}{1-x}=\frac{S(T)}{n_{\rm H}}.
$$
右辺を $s$ とすると $x^2+s x-s=0$ である。負の根は密度の条件を満たさない。正の根は $2/(1+\sqrt{1+4/s})$ と書くと、大きな $s$ においても、同程度の二つの数を差し引かずに計算できる。電離の密度感度は対数微分で求められる。温度を固定すると $(2/x+1/(1-x))\,\dd x=-\dd\ln n_{\rm H}$ だから、感度の絶対値は電離が進むほど小さくなる。同じ温度で密度を変えてこの感度を測れば、Saha 平衡の予測を追加で検証できる。外部の放射で電離を維持している気体にこの平衡定数を適用するには、放射と物質の熱平衡を別に確認する。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $n_{\rm H}$ | 水素核の密度 | $\mathrm{m^{-3}}$ |
| $s$ | Saha 比 $S(T)/n_{\rm H}$ | 1 |

電離度を $x$ とする。密度を非負にする根を選ぶ。

$$
x=\frac{2}{\sqrt{1 + \frac{4}{s}} + 1}
$$

電子密度を $n_e$ とする。各イオンは電子を一個放出する。

$$
n_e=n_{\rm H} x
$$

中性原子密度を $n_a$ とする。反応前後で水素核の総数は一定である。

$$
n_a=n_{\rm H} \left(1 - x\right)
$$

電子と中性原子の密度比を $n_e/n_a$ とする。中性粒子が衝突相手になる条件を読む第一歩である。衝突頻度には断面積も要る。

$$
n_e/n_a=\frac{x}{1 - x}
$$

電離度の密度感度を $\left.\partial\ln x/\partial\ln n_{\rm H}\right\rvert_T$ とする。$x^2/(1-x)=s$ の対数微分を用いる。

$$
\left.\partial\ln x/\partial\ln n_{\rm H}\right\rvert_T=\frac{x - 1}{2 - x}
$$

電子密度の感度を $\left.\partial\ln n_e/\partial\ln n_{\rm H}\right\rvert_T$ とする。$n_e=n_{\rm H}x$ なので水素核数の変化も加わる。

$$
\left.\partial\ln n_e/\partial\ln n_{\rm H}\right\rvert_T=\left.\partial\ln x/\partial\ln n_{\rm H}\right\rvert_T + 1
$$

:::example[ほぼ中性な平衡]
電離した粒子の割合が小さいとき、核密度を増やすと電子密度がどの程度増えるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
n_{\rm H}&=100000000000000000000\,\mathrm{m^{-3}} \\s&=\frac{1}{10000}
\end{aligned}
$$

これらを順に代入する。

電離度は

$$
x=\frac{2}{\sqrt{1 + \frac{4}{\left(\frac{1}{10000}\right)}} + 1}\simeq 0.0099501
$$

である。

電子密度は

$$
n_e=\left(100000000000000000000\right) \left(\frac{2}{1 + \sqrt{40001}}\right)\simeq 9.9501\times 10^{17}\,\mathrm{m^{-3}}
$$

である。

中性原子密度は

$$
n_a=\left(100000000000000000000\right) \left(1 - \left(\frac{2}{1 + \sqrt{40001}}\right)\right)\simeq 9.9005\times 10^{19}\,\mathrm{m^{-3}}
$$

である。

電子と中性原子の密度比は

$$
n_e/n_a=\frac{\left(\frac{2}{1 + \sqrt{40001}}\right)}{1 - \left(\frac{2}{1 + \sqrt{40001}}\right)}\simeq 0.01005
$$

である。

電離度の密度感度は

$$
\left.\partial\ln x/\partial\ln n_{\rm H}\right\rvert_T=\frac{\left(\frac{2}{1 + \sqrt{40001}}\right) - 1}{2 - \left(\frac{2}{1 + \sqrt{40001}}\right)}\simeq -0.4975
$$

である。

電子密度の感度は

$$
\left.\partial\ln n_e/\partial\ln n_{\rm H}\right\rvert_T=\left(\frac{-1 + \frac{2}{1 + \sqrt{40001}}}{2 - \frac{2}{1 + \sqrt{40001}}}\right) + 1\simeq 0.5025
$$

である。

電子密度の増加率は核密度の増加率のおよそ半分になる。この極限では $x\simeq\sqrt{s}$ であり、衝突を調べる場合には中性粒子の寄与も評価する。
:::

:::example[中性原子とイオンが同数になる条件]
$S(T)$ が核密度の半分に相当するとき、電離度が半分になるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
n_{\rm H}&=100000000000000000000\,\mathrm{m^{-3}} \\s&=\frac{1}{2}
\end{aligned}
$$

これらを順に代入する。

電離度は

$$
x=\frac{2}{\sqrt{1 + \frac{4}{\left(\frac{1}{2}\right)}} + 1}=\frac{1}{2}
$$

である。

電子密度は

$$
n_e=\left(100000000000000000000\right) \left(\frac{1}{2}\right)=50000000000000000000\,\mathrm{m^{-3}}
$$

である。

中性原子密度は

$$
n_a=\left(100000000000000000000\right) \left(1 - \left(\frac{1}{2}\right)\right)=50000000000000000000\,\mathrm{m^{-3}}
$$

である。

電子と中性原子の密度比は

$$
n_e/n_a=\frac{\left(\frac{1}{2}\right)}{1 - \left(\frac{1}{2}\right)}=1
$$

である。

電離度の密度感度は

$$
\left.\partial\ln x/\partial\ln n_{\rm H}\right\rvert_T=\frac{\left(\frac{1}{2}\right) - 1}{2 - \left(\frac{1}{2}\right)}=- \frac{1}{3}
$$

である。

電子密度の感度は

$$
\left.\partial\ln n_e/\partial\ln n_{\rm H}\right\rvert_T=\left(- \frac{1}{3}\right) + 1=\frac{2}{3}
$$

である。

$x=1/2$ であり、電子密度は核密度と単純な比例関係にはない。密度を一割増やす微小変更では、先頭次数で電子密度は約 $6.7$ パーセント増える。
:::

:::example[ほぼ完全電離した平衡]
中性原子が残る量と、密度を増やしたときの感度を見積もる。

入力値を次のように置く。

$$
\begin{aligned}
n_{\rm H}&=100000000000000000000\,\mathrm{m^{-3}} \\s&=100
\end{aligned}
$$

これらを順に代入する。

電離度は

$$
x=\frac{2}{\sqrt{1 + \frac{4}{\left(100\right)}} + 1}\simeq 0.9902
$$

である。

電子密度は

$$
n_e=\left(100000000000000000000\right) \left(\frac{2}{1 + \frac{\sqrt{26}}{5}}\right)\simeq 9.902\times 10^{19}\,\mathrm{m^{-3}}
$$

である。

中性原子密度は

$$
n_a=\left(100000000000000000000\right) \left(1 - \left(\frac{2}{1 + \frac{\sqrt{26}}{5}}\right)\right)\simeq 9.8049\times 10^{17}\,\mathrm{m^{-3}}
$$

である。

電子と中性原子の密度比は

$$
n_e/n_a=\frac{\left(\frac{2}{1 + \frac{\sqrt{26}}{5}}\right)}{1 - \left(\frac{2}{1 + \frac{\sqrt{26}}{5}}\right)}\simeq 100.99
$$

である。

電離度の密度感度は

$$
\left.\partial\ln x/\partial\ln n_{\rm H}\right\rvert_T=\frac{\left(\frac{2}{1 + \frac{\sqrt{26}}{5}}\right) - 1}{2 - \left(\frac{2}{1 + \frac{\sqrt{26}}{5}}\right)}\simeq -0.0097097
$$

である。

電子密度の感度は

$$
\left.\partial\ln n_e/\partial\ln n_{\rm H}\right\rvert_T=\left(\frac{-1 + \frac{2}{1 + \frac{\sqrt{26}}{5}}}{2 - \frac{2}{1 + \frac{\sqrt{26}}{5}}}\right) + 1\simeq 0.99029
$$

である。

中性原子は少数であるが、核密度が大きいので絶対密度は有限である。電離度の感度が小さいため、電子密度は核密度にほぼ比例する。
:::

:::example[同一温度での希釈]
$S(T)=5\times10^{19}\,\mathrm{m^{-3}}$ の気体を希釈したとき、電子密度と電離度を比較する。

入力値を次のように置く。

$$
\begin{aligned}
n_{\rm H}&=1000000000000000000\,\mathrm{m^{-3}} \\s&=50
\end{aligned}
$$

これらを順に代入する。

電離度は

$$
x=\frac{2}{\sqrt{1 + \frac{4}{\left(50\right)}} + 1}\simeq 0.98076
$$

である。

電子密度は

$$
n_e=\left(1000000000000000000\right) \left(\frac{2}{1 + \frac{3 \sqrt{3}}{5}}\right)\simeq 9.8076\times 10^{17}\,\mathrm{m^{-3}}
$$

である。

中性原子密度は

$$
n_a=\left(1000000000000000000\right) \left(1 - \left(\frac{2}{1 + \frac{3 \sqrt{3}}{5}}\right)\right)\simeq 1.9238\times 10^{16}\,\mathrm{m^{-3}}
$$

である。

電子と中性原子の密度比は

$$
n_e/n_a=\frac{\left(\frac{2}{1 + \frac{3 \sqrt{3}}{5}}\right)}{1 - \left(\frac{2}{1 + \frac{3 \sqrt{3}}{5}}\right)}\simeq 50.981
$$

である。

電離度の密度感度は

$$
\left.\partial\ln x/\partial\ln n_{\rm H}\right\rvert_T=\frac{\left(\frac{2}{1 + \frac{3 \sqrt{3}}{5}}\right) - 1}{2 - \left(\frac{2}{1 + \frac{3 \sqrt{3}}{5}}\right)}\simeq -0.018875
$$

である。

電子密度の感度は

$$
\left.\partial\ln n_e/\partial\ln n_{\rm H}\right\rvert_T=\left(\frac{-1 + \frac{2}{1 + \frac{3 \sqrt{3}}{5}}}{2 - \frac{2}{1 + \frac{3 \sqrt{3}}{5}}}\right) + 1\simeq 0.98113
$$

である。

二番目の例と $n_{\rm H}s$ が等しいので同一温度の比較になっている。電離度は増えるが、初めの核密度が小さいため、電子密度は減る。そのため遮蔽長や電子振動数は電離度だけでは比較できない。
:::


## Debye 遮蔽

正の試験電荷を置くと、電子が集まり、可動な正イオンは押し出される。どちらの変化も試験電荷の電場を弱める。遮蔽距離を求めるには、密度変化を [Poisson 方程式](/physics/electromagnetism/electrostatics/) に戻す。

### 静電平衡の密度

温度 $T_s$ が一定で平均流のない粒子種 $s$ では、圧力勾配と電気力密度のつりあいが $T_s\grad n_s=-q_s n_s\grad\phi$ である。連結領域で積分し、遠方に $n_s=n_{s0},\ \phi=0$ を課すと $n_s=n_{s0}e^{-q_s\phi/T_s}$ を得る。衝突がなくても、全粒子エネルギー $m_s v^2/2+q_s\phi$ の指数関数を平衡分布に選べば同じ密度になる。定常性は粒子軌道上でエネルギーが保存されることから従う。

:::definition[Debye 長]{#def-plasma-debye}
背景が中性で $\abs{q_s\phi}\ll T_s$ なら $q_s(n_s-n_{s0})=-n_{s0}q_s^2\phi/T_s+\order{\phi^2}$ である。応答できる種について和をとり、
$$
\lambda_D^{-2}=\sum_s\lambda_{Ds}^{-2},\qquad
\lambda_{Ds}^2=\frac{\varepsilon_0 T_s}{n_{s0}q_s^2}
\label{eq-plasma-debye-length}
$$
を定める。$\lambda_D$ を :term[Debye 長]{en="Debye length" yomi="でばいちょう"} という。
:::

### 点電荷と平面境界

:::proposition[線形遮蔽ポテンシャル]{#prop-plasma-yukawa}
無限一様媒質に外部点電荷 $Q$ を置くと、線形化した静電方程式の解は
$$
\phi(r)=\frac{Q}{4\pi\varepsilon_0r}e^{-r/\lambda_D}.
\label{eq-plasma-yukawa}
$$
:::

:::proof
:ref[law-em-gauss] により $(\nabla^2-\lambda_D^{-2})\phi=-Q\delta^{(3)}(\vb r)/\varepsilon_0$ である。球対称で $r>0$ なら $u=r\phi$ は $u''=\lambda_D^{-2}u$ を満たす。遠方で減衰する解は $u=C e^{-r/\lambda_D}$ である。半径 $\epsilon$ の球で方程式を積分すると $4\pi\epsilon^2\phi'(\epsilon)$ が残り、$\int_0^\epsilon r^2\phi\,\dd r=\order{\epsilon^2}$ は消える。極限から $-4\pi C=-Q/\varepsilon_0$ を得る。

解の差は原点で正則な同次解である。その差を $w$ とし、方程式に $w$ を掛けて空間で積分すると、表面項は減衰条件で消え、$\int(\abs{\grad w}^2+\lambda_D^{-2}w^2)\,\dd^3r=0$ となる。よって解は一意である。
:::

この解を :term[Debye 遮蔽]{en="Debye shielding" yomi="でばいしゃへい"} と呼ぶ。点電荷のすぐ近くでは線形化条件が破れる。その領域の微視的構造までこの解で記述することはできない。

:::example[外部電荷を打ち消す電荷]
誘導電荷は $\rho_{\mathrm{ind}}=-\varepsilon_0\phi/\lambda_D^2$ である。全空間で積分すると $-Q\lambda_D^{-2}\int_0^\infty r e^{-r/\lambda_D}\,\dd r=-Q$ となる。有限の球では遮蔽雲の外側が残るため、包む総電荷は一般にゼロにならない。
:::

:::example[壁による微小な電位変化]
$x>0$ に媒質があり $\phi(0)=\phi_w$ なら、$\phi''=\phi/\lambda_D^2$ の減衰解は $\phi_w e^{-x/\lambda_D}$ である。壁から $3\lambda_D$ の点では電位は壁の約 $5\,\%$ まで減る。大振幅の壁シースにはこの指数解を直接使えない。
:::

### 分布した外部電荷への応答

$\rho_{\rm ext}$ が点電荷に限らないとき、一定の $\lambda_D$ を持つ無限媒質では Fourier 変換で遮蔽を求められる。波数 $\vb k$ の成分は
$$
\phi_{\vb k}=\frac{\rho_{{\rm ext},\vb k}}{\varepsilon_0(k^2+\lambda_D^{-2})},\qquad
\rho_{{\rm ind},\vb k}=-\frac{\rho_{{\rm ext},\vb k}}{1+k^2\lambda_D^2}.
$$
前者は Poisson 方程式の $\nabla^2\mapsto-k^2$ と密度の一次応答から得られる。後者を外部電荷に足すと、残る正味電荷の割合は $k^2\lambda_D^2/(1+k^2\lambda_D^2)$ である。長波長ではほとんど全てが打ち消される。短波長では遮蔽雲を作る電荷の再分布が小さく、同じ外部電荷に対する電位は真空解に近づく。

:::example[長短二成分の違い]
外部電荷が $\rho_1\cos(x/\lambda_D)+\rho_1\cos(3x/\lambda_D)$ なら、残留成分はそれぞれ $1/2$ と $9/10$ 倍となる。電位は残留成分の大きさをさらに $k^2$ で割って決まるので、第 2 成分の電位は第 1 成分の $1/5$ である。
:::

### 一成分の非線形平面解

固定された一価イオンと Boltzmann 電子を半空間 $x>0$ に置く。$y=e\phi/T_e$、$X=x/\lambda_{De}$ とすれば $y''=e^y-1$ である。遠方で $y,y'\to0$ とし、$y'$ を掛けて積分すると
$$
\frac{(y')^2}{2}=e^y-1-y .
$$
右辺は指数関数の凸性により $y\ne0$ で正である。減衰枝は $y'=-\sgn(y)\sqrt{2(e^y-1-y)}$ であり、$y'$ の符号は壁面の電位と反対である。$e^y-1-y=y^2/2+\order{y^3}$ のため遠方では指数遮蔽に戻る。これは固定イオンという仮定の解であり、壁へ流れるイオンの密度変化は別のシース模型を用いて求める。

:::example[表面電場の比較]
壁で $y=-1$ なら $y'=\sqrt{2/e}\simeq0.858$ である。壁で $y=+1$ なら $y'=-\sqrt{2(e-2)}\simeq-1.20$ となる。線形解が同じ大きさを与える二つの境界値でも、電子がたまる側で電場が強い。
:::

### 二つの可動種による遮蔽と平面 Green 関数

点電荷から見た電子とイオンの再配置を分けて計算する。温度はエネルギー単位で測り、一次応答 $\delta n_s=-q_sn_{s0}\phi/T_s$ を使う。これは $T_s$ が一定で任意の種が位置を再分布でき、電位エネルギーが $T_s$ より小さいときのモデルである。Gaussian 単位系への換算は不要で、各種の逆遮蔽長の二乗は $\kappa_s^2=n_{s0}q_s^2/(\varepsilon_0 T_s)$ となる。Poisson 方程式に代入すると
$$
(-\nabla^2+\kappa^2)\phi=\rho_{\rm ext}/\varepsilon_0,
\qquad \kappa^2=\sum_s\kappa_s^2.
$$

面電荷 $\rho_{\rm ext}=\sigma\delta(x)$ では、$x\ne0$ の減衰解の連続性により $\phi=C e^{-\kappa|x|}$ である。電位を連続とするのは、隙間を挟む二重層に相当する $\delta'$ の電荷源を与えていないためである。原点をまたいで微分方程式を積分すると $-\phi'(0^+)+\phi'(0^-)=\sigma/\varepsilon_0$ となる。従って $C=\sigma/(2\varepsilon_0\kappa)$ である。種 $s$ の誘起面電荷は $\int_{-\infty}^{\infty}-\varepsilon_0\kappa_s^2\phi\,\dd x=-\sigma\kappa_s^2/\kappa^2$ である。このように、全遮蔽長には各種の逆二乗が加わり、打ち消す電荷の配分も同じ逆二乗の比で決まる。

空間的に滑らかな電荷源であれば、平面解を点 $x'$ に平行移動して重ねる。その結果、電位は $\phi(x)=\int e^{-\kappa|x-x'|}\rho_{\rm ext}(x')\,\dd x'/(2\varepsilon_0\kappa)$ となる。境界条件は無限遠で電位が消えることである。差を取って同次方程式を得ると、全実軸で滑らかで両側で減衰する解はゼロだから、この条件のもとで解は一意である。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $r$ | 温度比 $T_i/T_e$。一価かつ等しい背景密度とする。 | 1 |
| $X$ | 観測距離 $\abs{x}/\lambda_{De}$ | 1 |
| $A$ | 電子だけを可動とした場合の面上での振幅 $e\sigma\lambda_{De}/(2\varepsilon_0 T_e)$ | 1 |

逆遮蔽長の倍率を $K$ とする。$K=\kappa\lambda_{De}$ とおけば電子寄与は $1$、イオン寄与は $1/r$ となる。

$$
K=\sqrt{1 + \frac{1}{r}}
$$

面上の規格化電位を $U_0$ とする。微分の飛びを固定したので、遮蔽長が短くなると面上の電位も小さくなる。

$$
U_0=\frac{A}{K}
$$

観測点の電子規格化電位を $U_e(X)$ とする。$U_e=e\phi/T_e$ であり、線形化は面上の $U_0$ の絶対値で判定できる。

$$
U_e(X)=U_{0} e^{- K X}
$$

観測点のイオン規格化電位を $U_i(X)$ とする。一価なので $e\phi/T_i$ である。こちらの絶対値も小さい必要がある。

$$
U_i(X)=\frac{U_e(X)}{r}
$$

誘起電荷の電子分率を $f_e$ とする。電子種の逆長さの二乗を合計で割る。

$$
f_e=\frac{1}{K^{2}}
$$

誘起電荷のイオン分率を $f_i$ とする。二つの誘起電荷を足せば $-\sigma$ である。

$$
f_i=1 - f_{e}
$$

:::example[同温度の弱い面電荷]
面電荷から電子 Debye 長だけ離れた電位と電子・イオンの分担を求める。

入力値を次のように置く。

$$
\begin{aligned}
r&=1 \\X&=1 \\A&=\frac{1}{10}
\end{aligned}
$$

これらを順に代入する。

逆遮蔽長の倍率は

$$
K=\sqrt{1 + \frac{1}{\left(1\right)}}\simeq 1.4142
$$

である。

面上の規格化電位は

$$
U_0=\frac{\left(\frac{1}{10}\right)}{\left(\sqrt{2}\right)}\simeq 0.070711
$$

である。

観測点の電子規格化電位は

$$
U_e(X)=\left(\frac{\sqrt{2}}{20}\right) e^{- \left(1\right) \left(\sqrt{2}\right)}\simeq 0.017191
$$

である。

観測点のイオン規格化電位は

$$
U_i(X)=\frac{\left(\frac{\sqrt{2}}{20 e^{\sqrt{2}}}\right)}{\left(1\right)}\simeq 0.017191
$$

である。

誘起電荷の電子分率は

$$
f_e=\frac{1}{\left(\sqrt{2}\right)^{2}}=\frac{1}{2}
$$

である。

誘起電荷のイオン分率は

$$
f_i=1 - \left(\frac{1}{2}\right)=\frac{1}{2}
$$

である。

同じ温度でも電荷の符号が逆なので電子は増え、イオンは減る。電荷密度の変化はどちらも面電荷を打ち消す向きである。
:::

:::example[低温イオンの静的な応答]
イオン温度を低くして、電子規格化電位が小さいだけで十分かを確かめる。

入力値を次のように置く。

$$
\begin{aligned}
r&=\frac{1}{100} \\X&=\frac{1}{10} \\A&=\frac{1}{1000}
\end{aligned}
$$

これらを順に代入する。

逆遮蔽長の倍率は

$$
K=\sqrt{1 + \frac{1}{\left(\frac{1}{100}\right)}}\simeq 10.05
$$

である。

面上の規格化電位は

$$
U_0=\frac{\left(\frac{1}{1000}\right)}{\left(\sqrt{101}\right)}\simeq 9.9504\times 10^{-5}
$$

である。

観測点の電子規格化電位は

$$
U_e(X)=\left(\frac{\sqrt{101}}{101000}\right) e^{- \left(\frac{1}{10}\right) \left(\sqrt{101}\right)}\simeq 3.6423\times 10^{-5}
$$

である。

観測点のイオン規格化電位は

$$
U_i(X)=\frac{\left(\frac{\sqrt{101}}{101000 e^{\frac{\sqrt{101}}{10}}}\right)}{\left(\frac{1}{100}\right)}\simeq 0.0036423
$$

である。

誘起電荷の電子分率は

$$
f_e=\frac{1}{\left(\sqrt{101}\right)^{2}}=\frac{1}{101}
$$

である。

誘起電荷のイオン分率は

$$
f_i=1 - \left(\frac{1}{101}\right)=\frac{100}{101}
$$

である。

面におけるイオンの規格化振幅は $U_0/r$ で、約 $0.01$ である。従ってこの入力ではイオンについても線形化できる。有限時間でイオンが動けるかは別の時間尺度の比較で決まる。
:::

:::example[高温イオンを含む場合]
高温のイオンが可動であっても電子のみの Debye 長に近づくか調べる。

入力値を次のように置く。

$$
\begin{aligned}
r&=100 \\X&=2 \\A&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

逆遮蔽長の倍率は

$$
K=\sqrt{1 + \frac{1}{\left(100\right)}}\simeq 1.005
$$

である。

面上の規格化電位は

$$
U_0=\frac{\left(\frac{1}{100}\right)}{\left(\frac{\sqrt{101}}{10}\right)}\simeq 0.0099504
$$

である。

観測点の電子規格化電位は

$$
U_e(X)=\left(\frac{\sqrt{101}}{1010}\right) e^{- \left(2\right) \left(\frac{\sqrt{101}}{10}\right)}\simeq 0.0013333
$$

である。

観測点のイオン規格化電位は

$$
U_i(X)=\frac{\left(\frac{\sqrt{101}}{1010 e^{\frac{\sqrt{101}}{5}}}\right)}{\left(100\right)}\simeq 1.3333\times 10^{-5}
$$

である。

誘起電荷の電子分率は

$$
f_e=\frac{1}{\left(\frac{\sqrt{101}}{10}\right)^{2}}=\frac{100}{101}
$$

である。

誘起電荷のイオン分率は

$$
f_i=1 - \left(\frac{100}{101}\right)=\frac{1}{101}
$$

である。

打ち消しのほとんどを電子が担い、遮蔽長は電子だけの値より約 $0.5$ パーセント短い。この極限では静的な電位からイオン温度を推定する感度が低い。
:::

:::example[外部電荷が負の場合]
負の面電荷で電位と各種の分担がどう変わるかを調べる。

入力値を次のように置く。

$$
\begin{aligned}
r&=4 \\X&=\frac{1}{2} \\A&=- \frac{1}{10}
\end{aligned}
$$

これらを順に代入する。

逆遮蔽長の倍率は

$$
K=\sqrt{1 + \frac{1}{\left(4\right)}}\simeq 1.118
$$

である。

面上の規格化電位は

$$
U_0=\frac{\left(- \frac{1}{10}\right)}{\left(\frac{\sqrt{5}}{2}\right)}\simeq -0.089443
$$

である。

観測点の電子規格化電位は

$$
U_e(X)=\left(- \frac{\sqrt{5}}{25}\right) e^{- \left(\frac{1}{2}\right) \left(\frac{\sqrt{5}}{2}\right)}\simeq -0.051141
$$

である。

観測点のイオン規格化電位は

$$
U_i(X)=\frac{\left(- \frac{\sqrt{5}}{25 e^{\frac{\sqrt{5}}{4}}}\right)}{\left(4\right)}\simeq -0.012785
$$

である。

誘起電荷の電子分率は

$$
f_e=\frac{1}{\left(\frac{\sqrt{5}}{2}\right)^{2}}=\frac{4}{5}
$$

である。

誘起電荷のイオン分率は

$$
f_i=1 - \left(\frac{4}{5}\right)=\frac{1}{5}
$$

である。

電位は負なので電子は減り、イオンは増える。正の誘起電荷の分担は面電荷が正のときと同じで、振幅だけが符号反転する。これは線形方程式の重ね合わせと整合する。
:::


## プラズマ振動とプラズマ周波数

電子がイオンに対してずれると、遮蔽のときと同じ電気力が働く。ただし、変位が時間に依存する場合には電子の慣性があり、元の位置を通り過ぎる。

### 冷たい電子の変位

一様なイオン背景を固定し、電子流体の微小変位を $\xi(x,t)$、平衡密度を $n_0$ とする。粒子数保存 $n_e(x+\xi)(1+\partial_x\xi)=n_0$ により $\delta n_e=-n_0\partial_x\xi$ である。波数 $k\ne0$ の摂動に限れば $\partial_x E=e n_0\partial_x\xi/\varepsilon_0$ を積分した定数成分はゼロとなる。圧力と衝突を無視して $m_e\ddot\xi=-eE$ を代入すると
$$
\ddot\xi+\omega_{pe}^2\xi=0,\qquad
\omega_{ps}\defeq\sqrt{\frac{n_{s0}q_s^2}{\varepsilon_0m_s}}.
\label{eq-plasma-frequency}
$$
$\omega_{ps}$ を種 $s$ の :term[プラズマ角周波数]{en="plasma angular frequency" yomi="ぷらずまかくしゅうはすう"} という。Hz で測る周波数は $\omega_{ps}/(2\pi)$ である。

振動のエネルギー密度は $n_0m_e\dot\xi^2/2+\varepsilon_0 E^2/2$ である。時間微分すると第 1 項は $-e n_0 \dot\xi E$、第 2 項は $e n_0\dot\xi E$ なので和は保存される。電気エネルギーが支配的なのは最大変位の時点であり、運動エネルギーが支配的なのは変位ゼロの時点である。

### 可動イオンと衝突

イオンも冷たく一価なら $E=e n_0(\xi_e-\xi_i)/\varepsilon_0$ である。相対変位 $\xi_e-\xi_i$ の時間二階微分は $-(\omega_{pe}^2+\omega_{pi}^2)(\xi_e-\xi_i)$ となり、質量中心には電気力が残らない。電子に中性粒子との摩擦 $-m_e\nu\dot\xi$ が加わる場合、時間指数は $-\nu/2\pm\sqrt{\nu^2/4-\omega_{pe}^2}$ となる。

:::example[電子と陽子の応答時間]
電子と陽子の密度が同じなら $\omega_{pi}/\omega_{pe}=\sqrt{m_e/m_p}\simeq0.0233$ である。したがって電子周期だけを見る実験ではイオン固定は良い近似である。
:::

:::example[初期位相]
$\xi(0)=0$、$\dot\xi(0)=u_0$ なら $\xi(t)=(u_0/\omega_{pe})\sin(\omega_{pe}t)$ である。変位の振幅は $|u_0|/\omega_{pe}$ であり、その時点の電場エネルギー密度は初期の $n_0m_e u_0^2/2$ に等しい。
:::

### 駆動と自由振動

固定イオンの冷たい模型で、外部電場 $E_a\cos\omega t$ が電子に作用する場合を考える。電子の変位で生じる自己場とは別にこの力を入れると
$$
\ddot\xi+\nu\dot\xi+\omega_{pe}^2\xi=-(e/m_e)E_a\cos\omega t .
$$
複素振幅 $\xi_\omega$ は $-(e/m_e)E_a/(\omega_{pe}^2-\omega^2-i\nu\omega)$ である。電子が受ける外部電場の平均仕事率は $\frac12\Re[(-en_0)(-i\omega\xi_\omega)E_a^*]$ であり、
$$
\frac{n_0e^2\nu \omega^2\abs{E_a}^2}{2m_e[(\omega_{pe}^2-\omega^2)^2+\nu^2\omega^2]} .
$$
これは摩擦で失う平均運動エネルギーに等しい。自由振動の周波数 $\sqrt{\omega_{pe}^2-\nu^2/4}$ と定常駆動のエネルギー吸収の最大点 $\omega=\omega_{pe}$ は別の測り方である。後者は分母を $\omega^2$ で割り、$(\omega_{pe}^2/\omega-\omega)^2$ がゼロになるとき最大と分かる。

:::example[臨界減衰]
$\nu=2\omega_{pe}$ では自由方程式の指数が重根 $-\omega_{pe}$ である。$\xi(0)=a,\ \dot\xi(0)=0$ なら $\xi=a(1+\omega_{pe}t)e^{-\omega_{pe}t}$ である。単に指数 $ae^{-\omega_{pe}t}$ だけにすると、初期速度が指定に合わない。
:::

### 電子・イオン変位の初期値問題

境界の寄与を一様な層の相対変位で表す。単一電荷の一様な電子・イオンがともに動けるとし、初期に共通密度が $n_0$ である。各層の変位を $\xi_e,\xi_i$ とすると、内部の電場は $E=en_0(\xi_e-\xi_i)/\varepsilon_0$ である。これは層の境界に蓄積する電荷を Gauss の法則で評価した結果である。外部からの力と衝突を無視して、$m_e\ddot\xi_e=-eE$、$m_i\ddot\xi_i=eE$ を差し引く。総力がゼロであるため質量中心は等速で、相対変位は調和振動する。
$$
u=\xi_e-\xi_i,\qquad
X=\frac{m_e\xi_e+m_i\xi_i}{m_e+m_i},\qquad
\ddot u+(\omega_{pe}^2+\omega_{pi}^2)u=0,\quad \ddot X=0.
$$
$\omega_*^2=\omega_{pe}^2+\omega_{pi}^2$、初期条件を $u(0)=u_0,\dot u(0)=v_0$ とすると $u(t)=u_0\cos\omega_*t+(v_0/\omega_*)\sin\omega_*t$ である。運動エネルギーの和は質量中心成分と換算質量 $m_r=m_em_i/(m_e+m_i)$ の相対成分に分かれる。電場のエネルギー密度は $\varepsilon_0 E^2/2=n_0 m_r\omega_*^2 u^2/2$ なので、相対運動は粒子の運動エネルギーと場のエネルギーを交換する。$X$ だけを励起した初期条件は $\xi_e=\xi_i$ を保ち、電場を生じない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $r$ | 質量比 $m_e/m_i$ | 1 |
| $a$ | 初期相対変位 $u_0/L$ | 1 |
| $b$ | 初期相対速度 $v_0/(\omega_*L)$ | 1 |
| $\theta$ | 位相 $\omega_*t$。$L$ は任意の変位単位である。 | 1 |

電子周波数に対する角周波数を $\omega_*/\omega_{pe}$ とする。可動イオンの逆質量が復元力を受ける相対加速度に加わる。

$$
\omega_*/\omega_{pe}=\sqrt{r + 1}
$$

相対変位 $u/L$を $U$ とする。初期変位・初期速度から出る二つの項を位相に合わせて足す。

$$
U=a \cos{\left(\theta \right)} + b \sin{\left(\theta \right)}
$$

規格化相対速度を $V$ とする。時間微分して $\omega_* L$ で割る。

$$
V=- a \sin{\left(\theta \right)} + b \cos{\left(\theta \right)}
$$

質量中心から見た電子の変位を $(\xi_e-X)/L$ とする。質量中心の定義と相対変位を連立する。

$$
(\xi_e-X)/L=\frac{U}{r + 1}
$$

質量中心から見たイオンの変位を $(\xi_i-X)/L$ とする。電子とイオンの質量重み付き変位は互いに打ち消す。

$$
(\xi_i-X)/L=- \frac{U r}{r + 1}
$$

相対エネルギーに占める場の分率を $f_{\rm field}$ とする。共通因子 $n_0m_r\omega_*^2 L^2/2$ を除き $U^2/(U^2+V^2)$ を取る。

$$
f_{\rm field}=\frac{U^{2}}{a^{2} + b^{2}}
$$

:::example[陽子を伴う電子の初期変位]
電子だけをずらしてから四分の一周期後の、質量中心に対する状態を求める。

入力値を次のように置く。

$$
\begin{aligned}
r&=\frac{1}{1836} \\a&=1 \\b&=0 \\\theta&=\frac{\pi}{2}
\end{aligned}
$$

これらを順に代入する。

電子周波数に対する角周波数は

$$
\omega_*/\omega_{pe}=\sqrt{\left(\frac{1}{1836}\right) + 1}\simeq 1.0003
$$

である。

相対変位 $u/L$は

$$
U=\left(0\right) \sin{\left(\left(\frac{\pi}{2}\right) \right)} + \left(1\right) \cos{\left(\left(\frac{\pi}{2}\right) \right)}=0
$$

である。

規格化相対速度は

$$
V=\left(0\right) \cos{\left(\left(\frac{\pi}{2}\right) \right)} - \left(1\right) \sin{\left(\left(\frac{\pi}{2}\right) \right)}=-1
$$

である。

質量中心から見た電子の変位は

$$
(\xi_e-X)/L=\frac{\left(0\right)}{\left(\frac{1}{1836}\right) + 1}=0
$$

である。

質量中心から見たイオンの変位は

$$
(\xi_i-X)/L=- \frac{\left(0\right) \left(\frac{1}{1836}\right)}{\left(\frac{1}{1836}\right) + 1}=0
$$

である。

相対エネルギーに占める場の分率は

$$
f_{\rm field}=\frac{\left(0\right)^{2}}{\left(0\right)^{2} + \left(1\right)^{2}}=0
$$

である。

この時刻は電場がゼロになる時刻で、相対運動エネルギーが最大である。角周波数のイオンによる増加は約 $0.027$ パーセントである。
:::

:::example[電子・陽電子の対プラズマ]
対プラズマで半周期後の変位と場のエネルギーを求める。

入力値を次のように置く。

$$
\begin{aligned}
r&=1 \\a&=1 \\b&=0 \\\theta&=\pi
\end{aligned}
$$

これらを順に代入する。

電子周波数に対する角周波数は

$$
\omega_*/\omega_{pe}=\sqrt{\left(1\right) + 1}\simeq 1.4142
$$

である。

相対変位 $u/L$は

$$
U=\left(0\right) \sin{\left(\left(\pi\right) \right)} + \left(1\right) \cos{\left(\left(\pi\right) \right)}=-1
$$

である。

規格化相対速度は

$$
V=\left(0\right) \cos{\left(\left(\pi\right) \right)} - \left(1\right) \sin{\left(\left(\pi\right) \right)}=0
$$

である。

質量中心から見た電子の変位は

$$
(\xi_e-X)/L=\frac{\left(-1\right)}{\left(1\right) + 1}=- \frac{1}{2}
$$

である。

質量中心から見たイオンの変位は

$$
(\xi_i-X)/L=- \frac{\left(-1\right) \left(1\right)}{\left(1\right) + 1}=\frac{1}{2}
$$

である。

相対エネルギーに占める場の分率は

$$
f_{\rm field}=\frac{\left(-1\right)^{2}}{\left(0\right)^{2} + \left(1\right)^{2}}=1
$$

である。

両種は質量中心に対して同振幅で逆向きに動く。重い正電荷を固定した模型では、この周波数を再現できない。
:::

:::example[初期速度による励起]
電荷分離を初期に与えず、速度差を与えたときのエネルギー分配を求める。

入力値を次のように置く。

$$
\begin{aligned}
r&=\frac{1}{1836} \\a&=0 \\b&=2 \\\theta&=\frac{\pi}{6}
\end{aligned}
$$

これらを順に代入する。

電子周波数に対する角周波数は

$$
\omega_*/\omega_{pe}=\sqrt{\left(\frac{1}{1836}\right) + 1}\simeq 1.0003
$$

である。

相対変位 $u/L$は

$$
U=\left(0\right) \cos{\left(\left(\frac{\pi}{6}\right) \right)} + \left(2\right) \sin{\left(\left(\frac{\pi}{6}\right) \right)}=1
$$

である。

規格化相対速度は

$$
V=- \left(0\right) \sin{\left(\left(\frac{\pi}{6}\right) \right)} + \left(2\right) \cos{\left(\left(\frac{\pi}{6}\right) \right)}\simeq 1.7321
$$

である。

質量中心から見た電子の変位は

$$
(\xi_e-X)/L=\frac{\left(1\right)}{\left(\frac{1}{1836}\right) + 1}=\frac{1836}{1837}
$$

である。

質量中心から見たイオンの変位は

$$
(\xi_i-X)/L=- \frac{\left(1\right) \left(\frac{1}{1836}\right)}{\left(\frac{1}{1836}\right) + 1}=- \frac{1}{1837}
$$

である。

相対エネルギーに占める場の分率は

$$
f_{\rm field}=\frac{\left(1\right)^{2}}{\left(0\right)^{2} + \left(2\right)^{2}}=\frac{1}{4}
$$

である。

電場のエネルギーは全相対エネルギーの四分の一である。質量中心の並進が同時にあれば、装置から見た全エネルギーにはその運動も含める。
:::

:::example[位相の相殺]
初期変位と初期速度がともにある系で、電場が消える位相を調べる。

入力値を次のように置く。

$$
\begin{aligned}
r&=\frac{1}{4} \\a&=1 \\b&=1 \\\theta&=\frac{3 \pi}{4}
\end{aligned}
$$

これらを順に代入する。

電子周波数に対する角周波数は

$$
\omega_*/\omega_{pe}=\sqrt{\left(\frac{1}{4}\right) + 1}\simeq 1.118
$$

である。

相対変位 $u/L$は

$$
U=\left(1\right) \sin{\left(\left(\frac{3 \pi}{4}\right) \right)} + \left(1\right) \cos{\left(\left(\frac{3 \pi}{4}\right) \right)}=0
$$

である。

規格化相対速度は

$$
V=- \left(1\right) \sin{\left(\left(\frac{3 \pi}{4}\right) \right)} + \left(1\right) \cos{\left(\left(\frac{3 \pi}{4}\right) \right)}\simeq -1.4142
$$

である。

質量中心から見た電子の変位は

$$
(\xi_e-X)/L=\frac{\left(0\right)}{\left(\frac{1}{4}\right) + 1}=0
$$

である。

質量中心から見たイオンの変位は

$$
(\xi_i-X)/L=- \frac{\left(0\right) \left(\frac{1}{4}\right)}{\left(\frac{1}{4}\right) + 1}=0
$$

である。

相対エネルギーに占める場の分率は

$$
f_{\rm field}=\frac{\left(0\right)^{2}}{2 \left(1\right)^{2}}=0
$$

である。

変位を生む余弦と速度を生む正弦の寄与が打ち消す。その一方、速度では両方の項が同符号になるので、振動そのものが消えたという解釈はできない。
:::


## プラズマパラメータと結合の強さ

平均場で遮蔽を記述するには、遮蔽雲を作る粒子が十分に多いことが望ましい。同じ条件は、隣接粒子の相互作用が熱運動に比べて弱いことにも対応する。

:::definition[結合定数と Debye 球]{#def-plasma-coupling}
電荷の大きさ $q$、密度 $n$、温度 $T$ の同一種について、半径 $a=(3/(4\pi n))^{1/3}$ の球に平均一個の粒子が入る。結合定数を $\Gamma=q^2/(4\pi\varepsilon_0 aT)$、Debye 球の粒子数を $N_D=4\pi n\lambda_D^3/3$ とする。$a$ は :term[Wigner–Seitz 半径]{en="Wigner–Seitz radius" yomi="うぃぐなーさいつはんけい"} である。
:::

単一可動種と固定中和背景では、$n=3/(4\pi a^3)$ より $\lambda_D^2=a^2/(3\Gamma)$ である。したがって
$$
N_D=(3\Gamma)^{-3/2}.
\label{eq-plasma-coupling}
$$
別の慣例ではプラズマパラメータを $n\lambda_D^3$ やその逆数で表す。定数係数や大小の向きは定義を確認してから比較する。複数種の遮蔽を使った $N_D$ には上の単一種の関係をそのまま適用できない。

独立な粒子数が Poisson 分布に従う場合、相対標準偏差は $1/\sqrt{N_D}$ である。これは粗視化の目安である。遮蔽のために正負の粒子数には相関が生じ、正味電荷の揺らぎを独立な粒子数だけで決めることはできない。さらに古典条件として $n\ell_s^3/g_s\ll1$ を確認する。これは $\Gamma\ll1$ とは別の条件である。

:::example[弱結合の値]
$\Gamma=1/300$ なら $\lambda_D=10a$、$N_D=1000$ である。逆に $\Gamma=1/3$ では $N_D=1$ となり、平均場の線形遮蔽を微視的な相関の計算に使う前提が崩れる。
:::

:::example[圧縮と加熱]
密度が八倍になると $a$ は半分になる。同時に温度も二倍なら $aT$ は一定なので $\Gamma$ と $N_D$ は保存されるが、Debye 長は半分になる。
:::

### 距離の階層とクーロン衝突

同一種の電位エネルギーと温度が等しくなる距離を $b_T=q^2/(4\pi\varepsilon_0T)$ とする。$\lambda_D^2=\varepsilon_0 T/(nq^2)$ から
$$
b_T=\Gamma a,\qquad
\frac{\lambda_D}{b_T}=\frac{1}{\sqrt3\Gamma^{3/2}}=3N_D
$$
である。弱結合では $b_T\ll a\ll\lambda_D$ の順になる。熱速度での二体散乱は最初の距離で強く曲がり、それより大きい距離では小さな偏向が累積する。遮蔽で上限が切られると、その累積は二つの距離の比の対数になる。[Coulomb 衝突と輸送](/physics/plasma-physics/collisions/) では換算質量と相対速度を入れ、正確な小角の係数を求める。

:::example[相互作用の尺度]
$\Gamma=1/300$ なら $b_T=a/300$、$\lambda_D=10a$ なので $\lambda_D/b_T=3000$ である。この比の対数は約 $8.01$ であり、距離の比そのものを衝突率に掛けるわけではない。
:::

### 結合、古典性、遮蔽粒子数の独立な判定

平均粒子間隔の指標に Wigner–Seitz 半径 $a=(3/(4\pi n))^{1/3}$ を採用し、電子と価数一の静止背景を考える。$b_T=e^2/(4\pi\varepsilon_0 T)$ は Coulomb エネルギーが $T$ に達する長さである。$\Gamma=b_T/a$ という定義と $\lambda_{De}^2=\varepsilon_0 T/(n e^2)$ から、因子 $3$ を含む変換が得られる。
$$
\lambda_{De}^2=\frac{a^3}{3b_T},
\qquad
\frac{\lambda_{De}}a=\frac1{\sqrt{3\Gamma}},
\qquad
N_D=\frac{4\pi n\lambda_{De}^3}{3}=(3\Gamma)^{-3/2}.
$$
球内粒子数 $N_D$ を $n\lambda_{De}^3$ と定義する流儀は $4\pi/3$ だけ異なる。弱結合と古典性も分けて判定する。熱的波長 $\lambda_{\rm th}=h/\sqrt{2\pi m_e T}$ を用いると、非縮退の目安は $n\lambda_{\rm th}^3\ll1$ である。古典統計を使うためには状態当たりの占有数が小さくなければならない。相互作用が弱いかは別に $\Gamma$ によって調べる。

基準状態から密度を $A$ 倍、温度を $C$ 倍に変えると $a$ は $A^{-1/3}$ 倍、$b_T$ は $C^{-1}$ 倍である。従って結合は $A^{1/3}/C$ 倍、縮退の指標は $A/C^{3/2}$ 倍になる。この二つの感度は異なる。さらに短距離衝突の回折を考えるときは、占有数以外に相対運動量 $p$ に対する $\hbar/p$ と $b_{90}$ の比も必要である。前者が後者より長ければ、非縮退であっても散乱体が古典的な点軌道として見える最小衝突径数が変わる。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $A$ | 基準面に対する粒子密度の倍率 | 1 |
| $C$ | 基準状態に対する温度倍率 | 1 |
| $\Gamma_0$ | 基準状態の Coulomb 結合 | 1 |
| $D_0$ | 基準状態の縮退指標 $n_0\lambda_{{\rm th},0}^3$ | 1 |

変更後の結合を $\Gamma$ とする。間隔の短縮と熱運動の活発化を同時に計上する。

$$
\Gamma=\frac{\sqrt[3]{A} \Gamma_{0}}{C}
$$

電子 Debye 球の粒子数を $N_D$ とする。古典的線形遮蔽で計算した形式的な値である。

$$
N_D=\frac{\sqrt{3}}{9 \Gamma^{\frac{3}{2}}}
$$

間隔を単位とした遮蔽長を $\lambda_{De}/a$ とする。粒子間隔と同程度になれば連続な密度を使う遮蔽の前提を見直す。

$$
\lambda_{De}/a=\frac{\sqrt{3}}{3 \sqrt{\Gamma}}
$$

変更後の縮退指標を $D$ とする。波長が温度の平方根の逆数なので三乗は $C^{-3/2}$ となる。

$$
D=\frac{A D_{0}}{C^{\frac{3}{2}}}
$$

短距離と長距離の比を $b_T/\lambda_{De}$ とする。古典的な近接衝突の尺度と遮蔽の尺度を比べる。

$$
b_T/\lambda_{De}=\frac{\Gamma}{\lambda_{De}/a}
$$

距離比の対数を $L_T$ とする。$b_T$ を下限としたときの値である。十分な分離がなければ対数近似を使えない。

$$
L_T=\ln{\left(\frac{\lambda_{De}/a}{\Gamma} \right)}
$$

:::example[密度だけを上げた弱結合気体]
千倍の圧縮で、結合と量子統計の指標がどのように変わるかを判定する。

入力値を次のように置く。

$$
\begin{aligned}
A&=1000 \\C&=1 \\\Gamma_0&=\frac{1}{10000} \\D_0&=\frac{1}{100000000}
\end{aligned}
$$

これらを順に代入する。

変更後の結合は

$$
\Gamma=\frac{\sqrt[3]{\left(1000\right)} \left(\frac{1}{10000}\right)}{\left(1\right)}=\frac{1}{1000}
$$

である。

電子 Debye 球の粒子数は

$$
N_D=\frac{\sqrt{3}}{9 \left(\frac{1}{1000}\right)^{\frac{3}{2}}}\simeq 6085.8
$$

である。

間隔を単位とした遮蔽長は

$$
\lambda_{De}/a=\frac{\sqrt{3}}{3 \sqrt{\left(\frac{1}{1000}\right)}}\simeq 18.257
$$

である。

変更後の縮退指標は

$$
D=\frac{\left(1000\right) \left(\frac{1}{100000000}\right)}{\left(1\right)^{\frac{3}{2}}}=\frac{1}{100000}
$$

である。

短距離と長距離の比は

$$
b_T/\lambda_{De}=\frac{\left(\frac{1}{1000}\right)}{\left(\frac{10 \sqrt{30}}{3}\right)}\simeq 5.4772\times 10^{-5}
$$

である。

距離比の対数は

$$
L_T=\ln{\left(\frac{\left(\frac{10 \sqrt{30}}{3}\right)}{\left(\frac{1}{1000}\right)} \right)}\simeq 9.8123
$$

である。

結合は十倍、縮退指標は千倍になる。この入力では両方とも小さいので、古典弱結合模型の出発点を保てる。
:::

:::example[結合を一定に保つ圧縮]
加熱によって結合を保ったまま圧縮し、縮退指標を比較する。

入力値を次のように置く。

$$
\begin{aligned}
A&=1000 \\C&=10 \\\Gamma_0&=\frac{1}{1000} \\D_0&=\frac{1}{10000}
\end{aligned}
$$

これらを順に代入する。

変更後の結合は

$$
\Gamma=\frac{\sqrt[3]{\left(1000\right)} \left(\frac{1}{1000}\right)}{\left(10\right)}=\frac{1}{1000}
$$

である。

電子 Debye 球の粒子数は

$$
N_D=\frac{\sqrt{3}}{9 \left(\frac{1}{1000}\right)^{\frac{3}{2}}}\simeq 6085.8
$$

である。

間隔を単位とした遮蔽長は

$$
\lambda_{De}/a=\frac{\sqrt{3}}{3 \sqrt{\left(\frac{1}{1000}\right)}}\simeq 18.257
$$

である。

変更後の縮退指標は

$$
D=\frac{\left(1000\right) \left(\frac{1}{10000}\right)}{\left(10\right)^{\frac{3}{2}}}\simeq 0.0031623
$$

である。

短距離と長距離の比は

$$
b_T/\lambda_{De}=\frac{\left(\frac{1}{1000}\right)}{\left(\frac{10 \sqrt{30}}{3}\right)}\simeq 5.4772\times 10^{-5}
$$

である。

距離比の対数は

$$
L_T=\ln{\left(\frac{\left(\frac{10 \sqrt{30}}{3}\right)}{\left(\frac{1}{1000}\right)} \right)}\simeq 9.8123
$$

である。

遮蔽粒子数は不変である。しかし縮退指標は約 $31.6$ 倍に増える。結合を一定に保つ条件と量子統計の条件は一致しない。
:::

:::example[冷却で遮蔽模型が限界に近づく場合]
量子縮退の指標が小さくても古典的な遮蔽の計算が破れるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
A&=1 \\C&=\frac{1}{10} \\\Gamma_0&=\frac{1}{20} \\D_0&=\frac{1}{1000000}
\end{aligned}
$$

これらを順に代入する。

変更後の結合は

$$
\Gamma=\frac{\sqrt[3]{\left(1\right)} \left(\frac{1}{20}\right)}{\left(\frac{1}{10}\right)}=\frac{1}{2}
$$

である。

電子 Debye 球の粒子数は

$$
N_D=\frac{\sqrt{3}}{9 \left(\frac{1}{2}\right)^{\frac{3}{2}}}\simeq 0.54433
$$

である。

間隔を単位とした遮蔽長は

$$
\lambda_{De}/a=\frac{\sqrt{3}}{3 \sqrt{\left(\frac{1}{2}\right)}}\simeq 0.8165
$$

である。

変更後の縮退指標は

$$
D=\frac{\left(1\right) \left(\frac{1}{1000000}\right)}{\left(\frac{1}{10}\right)^{\frac{3}{2}}}\simeq 3.1623\times 10^{-5}
$$

である。

短距離と長距離の比は

$$
b_T/\lambda_{De}=\frac{\left(\frac{1}{2}\right)}{\left(\frac{\sqrt{6}}{3}\right)}\simeq 0.61237
$$

である。

距離比の対数は

$$
L_T=\ln{\left(\frac{\left(\frac{\sqrt{6}}{3}\right)}{\left(\frac{1}{2}\right)} \right)}\simeq 0.49041
$$

である。

$N_D$ は一個を下回り、$L_T\simeq0.49$ である。遠距離の切断と強散乱の距離が十分に離れていないため、この値を弱結合の衝突周波数へ代入できない。
:::

:::example[加熱を伴う圧縮]
密度と温度の両方を上げ、縮退指標を保ったときの結合を求める。

入力値を次のように置く。

$$
\begin{aligned}
A&=1000 \\C&=100 \\\Gamma_0&=\frac{1}{1000} \\D_0&=\frac{1}{10000}
\end{aligned}
$$

これらを順に代入する。

変更後の結合は

$$
\Gamma=\frac{\sqrt[3]{\left(1000\right)} \left(\frac{1}{1000}\right)}{\left(100\right)}=\frac{1}{10000}
$$

である。

電子 Debye 球の粒子数は

$$
N_D=\frac{\sqrt{3}}{9 \left(\frac{1}{10000}\right)^{\frac{3}{2}}}\simeq 1.9245\times 10^{5}
$$

である。

間隔を単位とした遮蔽長は

$$
\lambda_{De}/a=\frac{\sqrt{3}}{3 \sqrt{\left(\frac{1}{10000}\right)}}\simeq 57.735
$$

である。

変更後の縮退指標は

$$
D=\frac{\left(1000\right) \left(\frac{1}{10000}\right)}{\left(100\right)^{\frac{3}{2}}}=\frac{1}{10000}
$$

である。

短距離と長距離の比は

$$
b_T/\lambda_{De}=\frac{\left(\frac{1}{10000}\right)}{\left(\frac{100 \sqrt{3}}{3}\right)}\simeq 1.7321\times 10^{-6}
$$

である。

距離比の対数は

$$
L_T=\ln{\left(\frac{\left(\frac{100 \sqrt{3}}{3}\right)}{\left(\frac{1}{10000}\right)} \right)}\simeq 13.266
$$

である。

熱的波長の減少が高密度化を打ち消して $D$ は一定になる。その間に結合は一桁小さくなり、遮蔽球には約 $31.6$ 倍の粒子が入る。
:::


## 準中性条件

各種の密度がほぼ等しくても、圧力勾配を支える電場は残る。電荷差を捨てる操作と、電位をゼロにする操作を区別する。

:::definition[準中性]{#def-plasma-quasineutral}
電子と複数の正イオン種からなる系について、$\abs{\sum_s Z_sn_s-n_e}/n_e\ll1$ の状態を :term[準中性]{en="quasineutrality" yomi="じゅんちゅうせい"} と呼ぶ。
:::

電位のスケールを $T_e/e$、長さを $L$ にとり、$\phi=T_e\widehat\phi/e$、$\vb x=L\widehat{\vb x}$、基準密度を $n_0$ とする。Poisson 方程式は
$$
\frac{\sum_sZ_sn_s-n_e}{n_0}
=-\frac{\lambda_{De}^2}{L^2}\widehat\nabla^2\widehat\phi.
\label{eq-plasma-neutral-order}
$$
したがって無次元電位の曲率が有界なら密度差は $\order{(\lambda_{De}/L)^2}$ である。この近似は電位を消去するという意味ではない。例えば等温電子が速く平行方向に応答すれば、$E_\parallel=-(T_e/e)\partial_\parallel\ln n_e$ が残る。

### 壁シース

電子がイオンより速く吸収される壁では、壁が負に帯電し、壁近くに正味の正電荷が残ることがある。この層を :term[シース]{en="sheath" yomi="しーす"} と呼ぶ。冷たい一価イオンが速度 $u_0$ で周辺から入り、電子が Boltzmann 応答をすると仮定する。シース入口を $\phi=0$ とすれば、イオンの粒子束とエネルギーから
$$
n_i=\frac{n_0u_0}{\sqrt{u_0^2-2e\phi/m_i}},\qquad
n_e=n_0 e^{e\phi/T_e}.
$$
$\phi<0$ の入口近傍では $n_i-n_e\simeq n_0e\phi(1/(m_i u_0^2)-1/T_e)$ である。正の空間電荷を生むには $u_0\ge\sqrt{T_e/m_i}$ が必要である。等号では一次項が相殺されるため、高次項と外部領域との接続が重要になる。これはこのモデルの :term[Bohm 条件]{en="Bohm criterion" yomi="ぼーむじょうけん"} である。

:::example[体積内の誤差]
$L/\lambda_{De}=100$、$\widehat\nabla^2\widehat\phi=-2$ なら、基準密度に対する正の密度差は $2\times10^{-4}$ である。
:::

:::example[入口での空間電荷]
$u_0=2\sqrt{T_e/m_i}$、$e\phi/T_e=-0.01$ を使うと、一次近似で $(n_i-n_e)/n_0=0.0075$ である。入口速度が音速の半分なら同じ近似で符号が逆になり、単調なイオンシースへの接続条件を満たさない。
:::

### 冷イオンシースの第一積分

$u_0=M\sqrt{T_e/m_i}$、$y=e\phi/T_e$、$X=x/\lambda_{De}$ としてイオンが壁に流れる模型を規格化する。Poisson 方程式は $y''=e^y-(1-2y/M^2)^{-1/2}$ となる。入口から区間全体にわたる解を調べるため、$y'=0$ となる $y=0$ を基準に積分すると
$$
\frac12(y')^2=e^y-1+M^2(\sqrt{1-2y/M^2}-1).
$$
$y\to0$ の二次係数は $(1-M^{-2})/2$ である。$M<1$ では負になり、ゼロから負の電位へ単調につながる実解はない。$M=1$ では二次が消え、三次係数 $-1/3$ が $y<0$ で正になる。入口近くで $M>1$ なら $y'$ の大きさは $\sqrt{1-M^{-2}}|y|$ に比例し、入口へ向かう極限は指数的である。$M=1$ では $|y'|$ が $\sqrt{2/3}|y|^{3/2}$ に比例するため、漸近的に $|y|$ は入口側への距離の二乗に反比例する。

$M\ge1$ で負の電位全体への接続を確認するには、$w=-y>0$ とおく。第一積分の右辺の $w$ 微分は $(1+2w/M^2)^{-1/2}-e^{-w}$ である。$e^{2w}>1+2w\ge1+2w/M^2$ よりこれは正である。従って右辺は $w>0$ で正となり、壁に向かう向きの平方根を取れば、任意の有限な負の壁電位まで単調に接続できる。この結論は背景から電子を Boltzmann 分布で供給する冷イオン模型に限る。

:::example[壁の電位と負電荷の束]
壁へのイオン粒子束を $n_0u_0$ とし、入口で等方 Maxwell 電子の、壁へ向かう速度成分だけを数える。平面での障壁 $-e\phi_w>0$ を越えるには $v_x>\sqrt{-2e\phi_w/m_e}$ が要る。一次元 Maxwell 分布を $v_x$ で重み付けしてこの範囲で積分すると、電子束は $n_0\sqrt{T_e/(2\pi m_e)}e^{e\phi_w/T_e}$ である。放出も二次電子もない浮遊壁の電荷保存は両束の等しさを要求する。$u_0=\sqrt{T_e/m_i}$ なら $e\phi_w/T_e=\frac12\ln(2\pi m_e/m_i)$ である。陽子なら約 $-2.84$ であり、束を抑える負電位が自発的に生まれる。これは Maxwell 分布を供給する貯蔵領域と、途中に電位の極値がないシースを仮定した見積もりである。
:::

### 準中性近似で解く両極性の静電場

一次元で電子と一価イオンがともに静止し、種ごとの温度が空間的に一定であるとする。電子にだけ外力が働かず、正のイオンには単位粒子当たりの力 $F_i=-\partial_x U$ が働く。磁場と流れを無視すると運動量収支は $0=-T_e\partial_x n_e-en_e E$ および $0=-T_i\partial_x n_i+en_i E-n_i U'$ である。$n_i=n_e=n$ を最低次で課して二式を足し、さらに電子側へ戻すと
$$
\partial_x\ln n=-\frac{U'}{T_i+T_e},
\qquad
eE=\frac{T_e}{T_i+T_e}U',
\qquad
n(x)=n(0)\exp\!\left[-\frac{U(x)-U(0)}{T_i+T_e}\right].
$$
電場は電子に加わる圧力勾配を支え、イオンの外力の一部を打ち消す。この最低次で得た電場を Gauss の法則に代入すれば、次の次数で必要な電荷差の見積もりが得られる。
$$
\frac{n_i-n_e}{n}\simeq
\lambda_{De}(x)^2\,\frac{U''}{T_i+T_e},
\qquad
\lambda_{De}(x)^2=\frac{\varepsilon_0 T_e}{n(x)e^2}.
$$
従って一様な外力では電場が定数になり、この体積項による電荷密度差はゼロである。ただし、有限区間の両端で実際の電位を定めるには境界条件を課す。例えば壁の電位が体積解と整合しない場合には Debye 層が必要になる。逆に外力の勾配があり密度が低下する部分では、背景で小さく見積もった Debye 補正が大きくなり得る。最も希薄な領域でも誤差が小さいかを確認する。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $r$ | 温度比 $T_i/T_e$ | 1 |
| $u$ | 位置 $x$ の外力ポテンシャル差 $[U(x)-U(0)]/T_e$ | 1 |
| $f$ | 規格化外力勾配 $L U'(x)/T_e$ | 1 |
| $h$ | 規格化外力曲率 $L^2 U''(x)/T_e$ | 1 |
| $d$ | 基準位置での長さ比 $\lambda_{De}(0)/L$ | 1 |

密度比を $n/n(0)$ とする。全圧力がイオンの外力を支えるので指数の分母は温度の和である。

$$
n/n(0)=e^{- \frac{u}{r + 1}}
$$

規格化電場を $eEL/T_e$ とする。電子の運動量収支から求める。

$$
eEL/T_e=\frac{f}{r + 1}
$$

規格化電位差を $e[\phi(x)-\phi(0)]/T_e$ とする。$E=-\partial_x\phi$ を積分する。

$$
e[\phi(x)-\phi(0)]/T_e=- \frac{u}{r + 1}
$$

局所の Debye 長を $\lambda_{De}(x)/L$ とする。希薄な位置では同じ温度でも長くなる。

$$
\lambda_{De}(x)/L=\frac{d}{\sqrt{n/n(0)}}
$$

体積電荷の相対量を $(n_i-n_e)/n$ とする。Gauss の法則で得る最低次の推定値である。

$$
(n_i-n_e)/n=\frac{[\lambda_{De}(x)/L]^2 h}{1+r}
$$

:::example[一様な外力による成層]
温度が等しい二種が高度とともに希薄になるとき、電場と体積電荷を求める。

入力値を次のように置く。

$$
\begin{aligned}
r&=1 \\u&=2 \\f&=2 \\h&=0 \\d&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

密度比は

$$
n/n(0)=e^{- \frac{\left(2\right)}{\left(1\right) + 1}}\simeq 0.36788
$$

である。

規格化電場は

$$
eEL/T_e=\frac{\left(2\right)}{\left(1\right) + 1}=1
$$

である。

規格化電位差は

$$
e[\phi(x)-\phi(0)]/T_e=- \frac{\left(2\right)}{\left(1\right) + 1}=-1
$$

である。

局所の Debye 長は

$$
\lambda_{De}(x)/L=\frac{\left(\frac{1}{100}\right)}{\sqrt{\left(e^{-1}\right)}}\simeq 0.016487
$$

である。

体積電荷の相対量は

$$
(n_i-n_e)/n=\frac{\left(0\right) \left(\left(\frac{e^{\frac{1}{2}}}{100}\right)\right)^{2}}{\left(1\right) + 1}=0
$$

である。

電場が有限でもその空間微分はゼロなので、体積内の電荷はゼロである。一定の電場をどの境界条件で支えるかは、体積解と合わせて考える。
:::

:::example[冷イオンを支える電場]
曲率のある外力で、冷たいイオンにどの程度の電場が必要か求める。

入力値を次のように置く。

$$
\begin{aligned}
r&=\frac{1}{10} \\u&=1 \\f&=2 \\h&=2 \\d&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

密度比は

$$
n/n(0)=e^{- \frac{\left(1\right)}{\left(\frac{1}{10}\right) + 1}}\simeq 0.40289
$$

である。

規格化電場は

$$
eEL/T_e=\frac{\left(2\right)}{\left(\frac{1}{10}\right) + 1}=\frac{20}{11}
$$

である。

規格化電位差は

$$
e[\phi(x)-\phi(0)]/T_e=- \frac{\left(1\right)}{\left(\frac{1}{10}\right) + 1}=- \frac{10}{11}
$$

である。

局所の Debye 長は

$$
\lambda_{De}(x)/L=\frac{\left(\frac{1}{100}\right)}{\sqrt{\left(e^{- \frac{10}{11}}\right)}}\simeq 0.015755
$$

である。

体積電荷の相対量は

$$
(n_i-n_e)/n=\frac{\left(2\right) \left(\left(\frac{e^{\frac{5}{11}}}{100}\right)\right)^{2}}{\left(\frac{1}{10}\right) + 1}\simeq 0.00045128
$$

である。

$eE/U'=1/(1+r)$ であるため、この入力では外力の約九割を電気力が支える。曲率は正なので、推定される空間電荷も正である。
:::

:::example[希薄側で近似が悪くなる場合]
基準点の Debye 長だけでは誤差判定が不十分であることを確かめる。

入力値を次のように置く。

$$
\begin{aligned}
r&=1 \\u&=10 \\f&=4 \\h&=1 \\d&=\frac{1}{20}
\end{aligned}
$$

これらを順に代入する。

密度比は

$$
n/n(0)=e^{- \frac{\left(10\right)}{\left(1\right) + 1}}\simeq 0.0067379
$$

である。

規格化電場は

$$
eEL/T_e=\frac{\left(4\right)}{\left(1\right) + 1}=2
$$

である。

規格化電位差は

$$
e[\phi(x)-\phi(0)]/T_e=- \frac{\left(10\right)}{\left(1\right) + 1}=-5
$$

である。

局所の Debye 長は

$$
\lambda_{De}(x)/L=\frac{\left(\frac{1}{20}\right)}{\sqrt{\left(e^{-5}\right)}}\simeq 0.60912
$$

である。

体積電荷の相対量は

$$
(n_i-n_e)/n=\frac{\left(1\right) \left(\left(\frac{e^{\frac{5}{2}}}{20}\right)\right)^{2}}{\left(1\right) + 1}\simeq 0.18552
$$

である。

局所 Debye 長は $L$ の約六割になる。体積電荷の推定値も無視できないので、得られた数値を精密解とみなさず、Poisson 方程式を残して解き直す。
:::

:::example[曲率が負の区間]
外力の強度勾配が減少する区間で、空間電荷の符号を求める。

入力値を次のように置く。

$$
\begin{aligned}
r&=1 \\u&=1 \\f&=1 \\h&=-2 \\d&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

密度比は

$$
n/n(0)=e^{- \frac{\left(1\right)}{\left(1\right) + 1}}\simeq 0.60653
$$

である。

規格化電場は

$$
eEL/T_e=\frac{\left(1\right)}{\left(1\right) + 1}=\frac{1}{2}
$$

である。

規格化電位差は

$$
e[\phi(x)-\phi(0)]/T_e=- \frac{\left(1\right)}{\left(1\right) + 1}=- \frac{1}{2}
$$

である。

局所の Debye 長は

$$
\lambda_{De}(x)/L=\frac{\left(\frac{1}{100}\right)}{\sqrt{\left(e^{- \frac{1}{2}}\right)}}\simeq 0.01284
$$

である。

体積電荷の相対量は

$$
(n_i-n_e)/n=\frac{\left(-2\right) \left(\left(\frac{e^{\frac{1}{4}}}{100}\right)\right)^{2}}{\left(1\right) + 1}\simeq -0.00016487
$$

である。

電場は正の向きだが $E'$ は負になる。従って電子がイオンより多くなる。電場と空間電荷の符号が独立であることが分かる。
:::


## 宇宙・実験室プラズマの特徴的スケール

比較する現象によって代表尺度を選ぶ。静電的な密度応答には $\lambda_{De}$、非磁化電磁波の電流応答には $d_e=c/\omega_{pe}$ を使う。$c=1/\sqrt{\mu_0\varepsilon_0}$ は真空の光速である。磁化を判定するときは $\Omega_s=q_s B/m_s$ と熱速度 $v_{ts}=\sqrt{T_s/m_s}$ から、回転半径 $\rho_s=v_{ts}/\abs{\Omega_s}$ を得る。これらの定義から $\lambda_{De}=v_{te}/\omega_{pe}$、$d_e/\lambda_{De}=c/v_{te}$ である。

### 代表値の読み方

以下は比較のために選んだ一様模型の入力値である。陽子と電子の密度は等しい。表の $\lambda_{De}$ はイオンの静電応答を除いた値である。

| 模型 | $n_e$（$\mathrm{m^{-3}}$） | $T_e$（eV） | $B$（T） | $\lambda_{De}$（m） | $\omega_{pe}$（$\mathrm{s^{-1}}$） |
| --- | ---: | ---: | ---: | ---: | ---: |
| 希薄な太陽風 | $5\times10^6$ | $10$ | $5\times10^{-9}$ | $10.5$ | $1.26\times10^5$ |
| 低温放電 | $10^{16}$ | $2$ | $10^{-2}$ | $1.05\times10^{-4}$ | $5.64\times10^9$ |
| 磁気閉じ込めの高温領域 | $10^{20}$ | $10^4$ | $5$ | $7.43\times10^{-5}$ | $5.64\times10^{11}$ |

:::example[太陽風の電子尺度]
表の太陽風模型で $d_e\simeq 2.38\,\mathrm{km}$、$\rho_e\simeq1.51\,\mathrm{km}$ である。電荷の遮蔽は約 $10\,\mathrm m$ で起きるが、その事実から km スケールで電子慣性が小さいとは言えない。
:::

:::example[低温放電の時間尺度]
表の放電模型では $\abs{\Omega_e}=1.76\times10^9\,\mathrm{s^{-1}}$ である。電子と中性粒子の有効衝突周波数を $10^7\,\mathrm{s^{-1}}$ とすれば、その比は約 $176$ となる。電子は衝突間に何度も回転できる。波動に対して磁場の影響を判定するときは、衝突周波数の代わりに波の角周波数とも比較する。
:::

### 磁気圧と境界を含む尺度

電子・イオンの圧力を $p=sum_s n_sT_s$ とする。電子一種の半径が小さくても、磁場が流体の運動を強く拘束するかは圧力の比
$$
\beta=\frac{p}{B^2/(2\mu_0)}
$$
で別に調べる。例えば表の太陽風模型に $T_i=T_e$ を与えると $p\simeq1.60\times10^{-11}\,\mathrm{Pa}$、$B^2/(2\mu_0)\simeq9.95\times10^{-12}\,\mathrm{Pa}$ なので $\beta\simeq1.61$ である。放電模型で $T_i\ll T_e$ とすれば $p\simeq3.20\times10^{-3}\,\mathrm{Pa}$、磁気圧は約 $39.8\,\mathrm{Pa}$ なので $\beta\simeq8.05\times10^{-5}$ である。磁場を含む力のつりあいは [MHD 平衡と磁場閉じ込め](/physics/plasma-physics/equilibrium/) で扱う。

有限のプラズマでは、遮蔽の重なりと境界への粒子損失も別の問題である。例えば放電模型で平行な壁の間隔が $1\,\mathrm{cm}$ なら $L/\lambda_{De}\simeq95$ であり、線形遮蔽だけを考える限り小さな壁電位の摂動は中央にほとんど届かない。それでも壁への粒子束があれば、密度を維持するには体積内の電離や外部からの供給が必要である。

### 実験の尺度選択

低周波の流体模型では、波長と時間尺度に加え、衝突が局所 Maxwell 化を起こすかを確認する。仮に熱速度 $v_t$ と緩和率 $\nu$ が与えられれば、平均自由行程 $\ell_{\rm mfp}=v_t/\nu$ と巨視的長さ $L$ の比が局所の衝突閉包を使う目安になる。同じ媒質でも $L$ が異なると適用できる閉包が変わる。磁気的な低周波の電流応答では、$d_i=c/\omega_{pi}$ が $L$ と比べて小さいかを確認する。$d_i$ はイオン慣性長であり、$\lambda_{De}$ と異なる。

:::example[低周波でも二流体になる領域]
表の太陽風模型と陽子質量 $m_p$ を使うと $d_i=d_e\sqrt{m_p/m_e}\simeq102\,\mathrm{km}$ である。長さ $10\,\mathrm{km}$ では $\lambda_{De}/L\simeq10^{-3}$ だが $d_i/L$ は約 $10$ である。電荷差を小さく見積もれても、電子とイオンの速度差は無視できない。
:::

### 周波数・長さ・磁化を同時に比較する

電子の一方向熱速さは $v_{te}=\sqrt{T_e/m_e}$ である。温度を eV で示したときはジュールへ直して式に入れる。磁場の大きさ $B$ に対して電子の旋回角周波数は $\abs{\Omega_e}=eB/m_e$ であり、熱旋回半径は $\rho_e=v_{te}/\abs{\Omega_e}$ である。Debye 長と慣性長の定義を使うと
$$
\frac{\rho_e}{\lambda_{De}}=\frac{\omega_{pe}}{\abs{\Omega_e}},
\qquad
\frac{\lambda_{De}}{d_e}=\frac{v_{te}}{c},
\qquad
d_e=\frac{c}{\omega_{pe}}.
$$
第一の等式は両方の長さを同じ熱速さで定義したために成り立つ。第二の等式では電流の応答を伝える $c$ が現れる。これらの長さを測定器の分解能と比較すると、電荷分離と電流による慣性をどこまで落とせるか判断できる。

長さ $L$ と実験の角周波数 $\omega$ が与えられたとき、$L/\lambda_{De}$ と $\omega/\omega_{pe}$ を用いて電子が電荷分離を小さく保てるか調べる。それに加え $\omega/\abs{\Omega_e}$ と $L/\rho_e$ が軌道平均の条件を与える。これらは電子だけの判定であり、イオンに同じ模型を使うには $m_i$ を入れ直す。例えば一価イオンでは $\omega_{pi}/\omega_{pe}=\sqrt{m_e/m_i}$、$\abs{\Omega_i}/\abs{\Omega_e}=m_e/m_i$ である。衝突がある場合、場の時間変化が遅くても $\nu/\abs{\Omega_s}$ が大きければ一周の軌道を保てないため、さらに独立な比較が要る。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $A$ | 長さ比 $L/\lambda_{De}$ | 1 |
| $H$ | 周波数比 $\omega_{pe}/\abs{\Omega_e}$ | 1 |
| $V$ | 電子の熱速さ比 $v_{te}/c$ | 1 |
| $W$ | 駆動周波数比 $\omega/\omega_{pe}$ | 1 |

電荷分離を評価する長さ比を $\lambda_{De}/L$ とする。長波長の展開にはこの数が小さいことが必要である。

$$
\lambda_{De}/L=\frac{1}{A}
$$

電子旋回の空間比を $\rho_e/L$ とする。熱速さを両方の長さで共通に取っている。

$$
\rho_e/L=\frac{H}{A}
$$

電子慣性の長さ比を $d_e/L$ とする。電流層を扱うときは、この比と $\lambda_{De}/L$ を区別する。

$$
d_e/L=\frac{1}{A V}
$$

旋回に対する駆動の周波数比を $\omega/\abs{\Omega_e}$ とする。旋回の時間平均には、この値も小さいことが必要である。

$$
\omega/\abs{\Omega_e}=H W
$$

陽子のプラズマ振動に対する比を $\omega/\omega_{pi}$ とする。陽子と電子の質量比を $1836$ とした近似値である。

$$
\omega/\omega_{pi}=6 \sqrt{51} W
$$

陽子の旋回に対する比を $\omega/\abs{\Omega_i}$ とする。イオンの旋回は電子より遅いため、実験によってはイオンを別の模型で扱う。

$$
\omega/\abs{\Omega_i}=1836 \omega/\abs{\Omega_e}
$$

:::example[電子だけが磁化される駆動]
電子の旋回平均と陽子の旋回平均がともに可能かを調べる。

入力値を次のように置く。

$$
\begin{aligned}
A&=1000 \\H&=10 \\V&=\frac{1}{100} \\W&=\frac{1}{1000}
\end{aligned}
$$

これらを順に代入する。

電荷分離を評価する長さ比は

$$
\lambda_{De}/L=\frac{1}{\left(1000\right)}=\frac{1}{1000}
$$

である。

電子旋回の空間比は

$$
\rho_e/L=\frac{\left(10\right)}{\left(1000\right)}=\frac{1}{100}
$$

である。

電子慣性の長さ比は

$$
d_e/L=\frac{1}{\left(1000\right) \left(\frac{1}{100}\right)}=\frac{1}{10}
$$

である。

旋回に対する駆動の周波数比は

$$
\omega/\abs{\Omega_e}=\left(10\right) \left(\frac{1}{1000}\right)=\frac{1}{100}
$$

である。

陽子のプラズマ振動に対する比は

$$
\omega/\omega_{pi}=6 \sqrt{51} \left(\frac{1}{1000}\right)\simeq 0.042849
$$

である。

陽子の旋回に対する比は

$$
\omega/\abs{\Omega_i}=1836 \left(\frac{1}{100}\right)=\frac{459}{25}
$$

である。

電子については空間比と時間比が小さい。一方、陽子の駆動比は一を越えるので、電子で成り立つ平均を陽子に移すことはできない。
:::

:::example[準中性でも電子慣性が必要な薄層]
長さが多数の Debye 長に及ぶ場合に、電子慣性を落とせるか判定する。

入力値を次のように置く。

$$
\begin{aligned}
A&=100 \\H&=3 \\V&=\frac{1}{100} \\W&=\frac{1}{1000}
\end{aligned}
$$

これらを順に代入する。

電荷分離を評価する長さ比は

$$
\lambda_{De}/L=\frac{1}{\left(100\right)}=\frac{1}{100}
$$

である。

電子旋回の空間比は

$$
\rho_e/L=\frac{\left(3\right)}{\left(100\right)}=\frac{3}{100}
$$

である。

電子慣性の長さ比は

$$
d_e/L=\frac{1}{\left(100\right) \left(\frac{1}{100}\right)}=1
$$

である。

旋回に対する駆動の周波数比は

$$
\omega/\abs{\Omega_e}=\left(3\right) \left(\frac{1}{1000}\right)=\frac{3}{1000}
$$

である。

陽子のプラズマ振動に対する比は

$$
\omega/\omega_{pi}=6 \sqrt{51} \left(\frac{1}{1000}\right)\simeq 0.042849
$$

である。

陽子の旋回に対する比は

$$
\omega/\abs{\Omega_i}=1836 \left(\frac{3}{1000}\right)=\frac{1377}{250}
$$

である。

$d_e/L=1$ である。電荷分離を抑えた模型でも電流の慣性は競合し得るので、長尺度の理想 MHD への縮約は別途検討する。
:::

:::example[低密度の強磁化を想定した入力]
電子と陽子が旋回を完了する条件を周波数から比べる。

入力値を次のように置く。

$$
\begin{aligned}
A&=10000 \\H&=\frac{1}{100} \\V&=\frac{1}{100} \\W&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

電荷分離を評価する長さ比は

$$
\lambda_{De}/L=\frac{1}{\left(10000\right)}=\frac{1}{10000}
$$

である。

電子旋回の空間比は

$$
\rho_e/L=\frac{\left(\frac{1}{100}\right)}{\left(10000\right)}=\frac{1}{1000000}
$$

である。

電子慣性の長さ比は

$$
d_e/L=\frac{1}{\left(10000\right) \left(\frac{1}{100}\right)}=\frac{1}{100}
$$

である。

旋回に対する駆動の周波数比は

$$
\omega/\abs{\Omega_e}=\left(\frac{1}{100}\right)^{2}=\frac{1}{10000}
$$

である。

陽子のプラズマ振動に対する比は

$$
\omega/\omega_{pi}=6 \sqrt{51} \left(\frac{1}{100}\right)\simeq 0.42849
$$

である。

陽子の旋回に対する比は

$$
\omega/\abs{\Omega_i}=1836 \left(\frac{1}{10000}\right)=\frac{459}{2500}
$$

である。

陽子の比は約 $0.18$ であり、先頭近似の精度を確かめる必要がある。また、陽子の旋回半径を判定するには陽子温度を別に与えなければならない。
:::

:::example[横変位が旋回半径に達する観測]
周波数が低くても旋回円を無視できるかを判定する。

入力値を次のように置く。

$$
\begin{aligned}
A&=30 \\H&=30 \\V&=\frac{1}{50} \\W&=\frac{1}{10000}
\end{aligned}
$$

これらを順に代入する。

電荷分離を評価する長さ比は

$$
\lambda_{De}/L=\frac{1}{\left(30\right)}=\frac{1}{30}
$$

である。

電子旋回の空間比は

$$
\rho_e/L=1
$$

である。

電子慣性の長さ比は

$$
d_e/L=\frac{1}{\left(30\right) \left(\frac{1}{50}\right)}=\frac{5}{3}
$$

である。

旋回に対する駆動の周波数比は

$$
\omega/\abs{\Omega_e}=\left(30\right) \left(\frac{1}{10000}\right)=\frac{3}{1000}
$$

である。

陽子のプラズマ振動に対する比は

$$
\omega/\omega_{pi}=6 \sqrt{51} \left(\frac{1}{10000}\right)\simeq 0.0042849
$$

である。

陽子の旋回に対する比は

$$
\omega/\abs{\Omega_i}=1836 \left(\frac{3}{1000}\right)=\frac{1377}{250}
$$

である。

電子旋回の時間平均は可能だが、場の空間変動を旋回円の一点で代表する近似は使えない。有限旋回半径を保持するときは [磁化プラズマの運動論](/physics/plasma-physics/magnetized-kinetics/) の平均を使う。
:::


## 演習問題

::::exercise{level=basic}
一価イオン密度 $3\times10^{17}\,\mathrm{m^{-3}}$、中性原子密度 $7\times10^{17}\,\mathrm{m^{-3}}$ の電離度を求めよ。

:::solution
核の総数密度は $10^{18}\,\mathrm{m^{-3}}$ なので $x=3/10$ である。
:::
::::

::::exercise{level=basic}
$\lambda_{Di}=2\lambda_{De}$ のとき全 Debye 長を求めよ。

:::solution
逆二乗を足すと $\lambda_D^{-2}=(5/4)\lambda_{De}^{-2}$ だから、正の根をとり $\lambda_D=2\lambda_{De}/\sqrt5$ を得る。
:::
::::

::::exercise{level=basic}
点電荷の電位と真空中の電位との比が $1/2$ となる距離を求めよ。

:::solution
共通の $Q/(4\pi\varepsilon_0 r)$ が消え、$e^{-r/\lambda_D}=1/2$ なので $r=\lambda_D\ln2$ である。
:::
::::

::::exercise{level=basic}
密度を九倍にしたとき電子プラズマ周期がどう変わるか求めよ。

:::solution
$\omega_{pe}^2\propto n$ なので $\omega_{pe}$ は三倍、周期 $2\pi/\omega_{pe}$ は三分の一になる。
:::
::::

::::exercise{level=advanced}
遮蔽された点電荷を中心に半径 $R$ の球をとる。球内の総電荷を求めよ。

:::solution
Gauss 則から $Q_{\rm tot}=-4\pi\varepsilon_0R^2\phi'(R)$ である。$r^{-1}e^{-r/\lambda_D}$ を微分すると $Q_{\rm tot}=Q(1+R/\lambda_D)e^{-R/\lambda_D}$ となる。$R$ がゼロに近づくと $Q$ に戻る。
:::
::::

::::exercise{level=advanced}
固定イオンの模型で初期条件 $\xi(0)=a,\dot\xi(0)=0$ とし、時間平均した電場エネルギーを求めよ。

:::solution
解は $a\cos\omega_{pe}t$ であり $E=e n_0a\cos\omega_{pe}t/\varepsilon_0$ である。$\cos^2$ の周期平均は $1/2$ なので電場エネルギー密度の平均は $n_0m_e\omega_{pe}^2a^2/4$ である。
:::
::::

::::exercise{level=advanced}
単一種で温度を四倍、密度を二倍にすると $N_D$ は何倍になるか。

:::solution
$\lambda_D^2\propto T/n$ なので長さは $\sqrt2$ 倍である。密度の倍率も掛けると $N_D$ は $2(\sqrt2)^3=4\sqrt2$ 倍である。
:::
::::

::::exercise{level=advanced}
$e\phi/T_e=-\delta\cos(kx)$ とする。$\delta>0$ のときどこで正の電荷密度が最大か。

:::solution
Poisson 方程式から $\rho_q=-\varepsilon_0(T_e/e)\delta k^2\cos(kx)$ である。$k\ne0$ なら $kx=(2j+1)\pi$ の点で最大となり、最大値は $\varepsilon_0(T_e/e)\delta k^2$ である。
:::
::::

## 確認テスト

::::quiz
両方の粒子種が静電遮蔽に応答すると全 Debye 長はどうなるか。

- [x] 各種の長さより短い。
- [ ] 各種の長さの和に等しい。
- [ ] 最大の長さに等しい。

:::explanation
有限の正の $\lambda_{De},\lambda_{Di}$ について逆二乗の和は各逆二乗より大きい。逆数の平方根を取ると各長さより小さい。
:::
::::

::::quiz{answer="4" tol="1e-6"}
密度を十六倍にしたときの冷たい電子プラズマ角周波数の倍率を答えよ。

:::explanation
$\omega_{pe}=\sqrt{n_e e^2/(\varepsilon_0 m_e)}$ なので倍率は $\sqrt{16}=4$ である。
:::
::::

::::quiz
準中性近似について正しいものを選べ。

- [x] 圧力勾配とつりあう電場を許す。
- [ ] 任意の境界で正味電荷を厳密にゼロにする。
- [ ] 電位を至る所で一定にする。

:::explanation
等温電子の力のつりあいは $E_\parallel=-(T_e/e)\partial_\parallel\ln n_e$ であり、密度勾配があれば電場が残る。壁シースでは電荷分離も残る。
:::
::::

::::quiz{answer="3" tol="1e-6"}
$N_D=27$ の単一種模型で $\lambda_D/a$ を答えよ。

:::explanation
定義より $N_D=(\lambda_D/a)^3$ である。正の立方根をとり $\lambda_D/a=3$ を得る。
:::
::::

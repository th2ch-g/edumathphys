---
title: 核融合プラズマ
description: 核反応率とプラズマのエネルギー収支から、閉じ込めと点火の条件を求める。
notation:
  - symbol: '\langle\sigma v\rangle'
    meaning: 熱核反応の反応度
  - symbol: 'Q'
    meaning: プラズマの核融合出力と外部加熱の比
  - symbol: '\tau_E'
    meaning: エネルギー閉じ込め時間
  - symbol: '\rho R'
    meaning: 面密度
---

## この章の目標

- 熱速度分布から反応率を導く。
- 点火と有限利得の条件を区別する。
- 磁気閉じ込めと慣性閉じ込めをエネルギー収支で比較する。

## 核融合反応と反応率

軽い核の融合では、核同士が近づく確率と、その後の反応確率を断面積で表す。反応のエネルギーは [原子核反応](/physics/nuclear-physics/reactions/) の質量差から得られる。D と T の反応は $\mathrm D+\mathrm T\to\alpha+n$ であり、発熱量は約 $17.6\,\mathrm{MeV}$ である。初期の質量中心運動を無視すれば反応後の二粒子の運動量の大きさが等しい。非相対論近似ではエネルギーが質量の逆比で分かれ、$\alpha$ に約 $3.5\,\mathrm{MeV}$、中性子に約 $14.1\,\mathrm{MeV}$ が渡る。

代表的な反応には、$\mathrm D+\mathrm D\to{}^3\mathrm{He}+n$ と $\mathrm D+\mathrm D\to \mathrm T+p$ の二枝や、$\mathrm D+{}^3\mathrm{He}\to\alpha+p$ もある。二つの DD 枝では同じ反応対を数え、生成物ごとに断面積と放出エネルギーを区別する。D と $\,{}^3\mathrm{He}$ の反応は中性子を直接出さないが、D を含む燃料では DD 枝も起こり得る。従って炉のエネルギー収支には、主反応だけでなく副反応による生成物も入れる。 :cite[chen-plasma]

### 相対速度の平均

独立な二粒子の速度から相対速度 $\vb u$ を取り、断面積 $\sigma(u)$ を相対速度と掛けて平均する。体積・時間当たりの反応数は
$$
R_{ab}=\frac{n_a n_b}{1+\delta_{ab}}\langle\sigma u\rangle .
\label{eq-plasma-fusion-rate}
$$
同種なら順序付きの対が二度数えられるため分母は二である。両種が温度 $T_i$ の Maxwell 分布なら、$\vb u$ の各成分の分散は $T_i(1/m_a+1/m_b)=T_i/m_r$ である。換算質量 $m_r$ で $E=m_ru^2/2$ と変数変換すると :term[反応度]{en="reactivity" yomi="はんのうど"} は
$$
\langle\sigma u\rangle=\sqrt{\frac{8}{\pi m_r}}\frac1{T_i^{3/2}}
\int_0^\infty \sigma(E)E e^{-E/T_i}\dd E .
$$
単位は $\mathrm{m^3\,s^{-1}}$ である。トンネル透過が支配的で $E\sigma(E)\simeq S(E)e^{-\sqrt{E_G/E}}$ と書ける範囲では、$S(E)$ が緩やかなら指数の最大寄与は $E/T_i+\sqrt{E_G/E}$ の極小にある。微分から $E_0=(E_G T_i^2/4)^{1/3}$ である。DT の具体的な反応度は核の共鳴も含む断面積を積分して求める。

:::example[燃料比]
全イオン密度 $n=n_D+n_T$、混合比 $x=n_D/n$ を固定すると $R_{DT}=n^2 x(1-x)\langle\sigma u\rangle$ である。積 $x(1-x)=1/4-(x-1/2)^2$ は等量で最大となる。
:::

:::example[反応出力の単位]
$n_D=n_T=5\times10^{19}\,\mathrm{m^{-3}}$、例として $\langle\sigma u\rangle=10^{-22}\,\mathrm{m^3\,s^{-1}}$ を与えると反応率は $2.5\times10^{17}\,\mathrm{m^{-3}s^{-1}}$ である。$17.6\,\mathrm{MeV}$ を J に直すと体積当たりの出力は約 $0.705\,\mathrm{MW\,m^{-3}}$ となる。
:::

### 二つの温度での相対速度

イオン種の温度が $T_a,T_b$ と異なっても、両方が等方 Gaussian で平均流が同じなら相対速度の分散は $T_a/m_a+T_b/m_b$ である。従って反応度の式に入る有効温度は $T_{\rm rel}=m_r(T_a/m_a+T_b/m_b)$ である。平均流が異なる場合には非中心 Gaussian になるので、有効温度への置換だけでは反応率を計算できない。

:::example[等量混合と温度差]
質量を $m_D=2m_p$、$m_T=3m_p$ と近似し、$T_D=10\,\mathrm{keV}$、$T_T=5\,\mathrm{keV}$ とすれば $T_{\rm rel}=8\,\mathrm{keV}$ である。算術平均 $7.5\,\mathrm{keV}$ と異なるのは軽い側の速度の分散が大きいためである。
:::

### 反応に寄与するエネルギー幅

温度を下げると Maxwell 分布の高速成分が減り、逆に極めて低速の衝突ではクーロン障壁の透過率が小さくなる。両者を合わせた指数を $h(E)=E/T+\sqrt{E_G/E}$ とする。$E_0$ で $h(E_0)=3E_0/T$、$h''(E_0)=3/(2TE_0)$ である。$E_0/T\gg1$、$S(E)$ がこの幅で緩やかなら、二次近似したガウス積分により
$$
\langle\sigma u\rangle\simeq
\sqrt{\frac8{\pi m_r}}T^{-3/2}S(E_0)
e^{-3E_0/T}\sqrt{\frac{4\pi TE_0}{3}}
$$
となる。標準偏差は $\sqrt{2TE_0/3}$ であり、温度に変化があると反応率だけでなく、積分を支配するエネルギーも移る。この推定は核の共鳴構造が幅内にない場合を想定する。

:::example[幅と平均エネルギー]
$E_0/T=10$ なら相対的な標準偏差は $\sqrt{2/(3\cdot10)}=1/\sqrt{15}$ である。積分への寄与には有限の幅があるため、試験粒子のエネルギーを $E_0$ の一点に固定した値と、熱分布での反応度はそのまま比較できない。
:::

### 相対速度に平均流がある場合の反応度

二つの等方 Maxwell 分布の相対速度が平均 $\vb U$、各成分の分散 $v_r^2$ を持つとする。速度空間の極軸を $\vb U$ に選ぶと $F\propto e^{-(u^2+U^2)/(2v_r^2)}e^{uU\cos\theta/v_r^2}$ となる。角積分は $2\sinh y/y$、$y=uU/v_r^2$ であり、相対速度の大きさの確率密度は
$$
p(u)=\sqrt{\frac2\pi}\frac{u^2}{v_r^3}e^{-(u^2+U^2)/(2v_r^2)}
\frac{\sinh(uU/v_r^2)}{uU/v_r^2}.
$$
$U\to0$ の極限は通常の Maxwell の速さ分布である。反応度は $\int_0^\infty \sigma(u)u p(u)\dd u$ で評価する。平均エネルギーを合わせた温度だけで置き換えられないことを多項式の断面積で確かめる。$\sigma(u)=C u$ なら $\langle\sigma u\rangle=C(U^2+3v_r^2)$ である。$\sigma=C' u^3$ では $\langle u^4\rangle=U^4+10U^2v_r^2+15v_r^4$ である。後式は $\vb u=\vb U+\vb w$ の四乗を展開し、$\langle|\vb w|^4\rangle=15v_r^4$ と $\langle(\vb U\cdot\vb w)^2\rangle=U^2v_r^2$ を入れて求める。これらは平均の仕方を確認する断面積模型である。実際の DT 反応の断面積には核の反応データを用いる。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $r$ | 相対平均流と熱速さの比 $U/v_r$ | 1 |
| $x$ | 観測する相対速さの比 $u/v_r$ | 1 |

平均エネルギーを合わせた等方温度比を $H$ とする。平均流の運動エネルギーも熱へ入れた仮想値である。

$$
H=\frac{r^{2}}{3} + 1
$$

規格化した速さ分布を $v_r p(u)$ とする。$xr\ne0$ で角積分を評価する。

$$
v_r p(u)=\frac{\sqrt{2} x \left(e^{r x} - e^{- r x}\right) e^{- \frac{r^{2}}{2} - \frac{x^{2}}{2}}}{2 \sqrt{\pi} r}
$$

速さの二乗平均を $\langle u^2\rangle/v_r^2$ とする。熱と流れの両方を数える。

$$
\langle u^2\rangle/v_r^2=r^{2} + 3
$$

速さの四乗平均を $\langle u^4\rangle/v_r^4$ とする。方向成分の Gaussian モーメントを展開する。

$$
\langle u^4\rangle/v_r^4=r^{4} + 10 r^{2} + 15
$$

等方温度で置き換えた四乗平均を $M_{\rm iso}$ とする。二乗平均は厳密値と同じに取る。

$$
M_{\rm iso}=15 H^{2}
$$

四乗に対する置換の相対誤差を $\delta_4$ とする。強い速度依存の反応で差が増え得る。

$$
\delta_4=\frac{M_{\rm iso}}{\langle u^4\rangle/v_r^4} - 1
$$

:::example[弱い相対ドリフト]
温度への換算誤差が小さい場合を求める。

入力値を次のように置く。

$$
\begin{aligned}
r&=\frac{1}{10} \\x&=1
\end{aligned}
$$

これらを順に代入する。

平均エネルギーを合わせた等方温度比は

$$
H=\frac{\left(\frac{1}{10}\right)^{2}}{3} + 1=\frac{301}{300}
$$

である。

規格化した速さ分布は

$$
v_r p(u)=\frac{\sqrt{2} \left(1\right) \left(e^{\left(1\right) \left(\frac{1}{10}\right)} - e^{- \left(1\right) \left(\frac{1}{10}\right)}\right) e^{- \frac{\left(1\right)^{2}}{2} - \frac{\left(\frac{1}{10}\right)^{2}}{2}}}{2 \sqrt{\pi} \left(\frac{1}{10}\right)}\simeq 0.48233
$$

である。

速さの二乗平均は

$$
\langle u^2\rangle/v_r^2=\left(\frac{1}{10}\right)^{2} + 3=\frac{301}{100}
$$

である。

速さの四乗平均は

$$
\langle u^4\rangle/v_r^4=\left(\frac{1}{10}\right)^{4} + 10 \left(\frac{1}{10}\right)^{2} + 15=\frac{151001}{10000}
$$

である。

等方温度で置き換えた四乗平均は

$$
M_{\rm iso}=15 \left(\frac{301}{300}\right)^{2}=\frac{90601}{6000}
$$

である。

四乗に対する置換の相対誤差は

$$
\delta_4=-1 + \frac{\left(\frac{90601}{6000}\right)}{\left(\frac{151001}{10000}\right)}=\frac{2}{453003}
$$

である。

平均流補正は弱く、二乗と四乗とも熱運動が主である。
:::

:::example[熱速さほどの流れ]
平均エネルギーが一致しても四乗が一致するか調べる。

入力値を次のように置く。

$$
\begin{aligned}
r&=1 \\x&=2
\end{aligned}
$$

これらを順に代入する。

平均エネルギーを合わせた等方温度比は

$$
H=\frac{\left(1\right)^{2}}{3} + 1=\frac{4}{3}
$$

である。

規格化した速さ分布は

$$
v_r p(u)=\frac{\sqrt{2} \left(2\right) \left(e^{\left(1\right) \left(2\right)} - e^{- \left(1\right) \left(2\right)}\right) e^{- \frac{\left(1\right)^{2}}{2} - \frac{\left(2\right)^{2}}{2}}}{2 \sqrt{\pi} \left(1\right)}\simeq 0.47508
$$

である。

速さの二乗平均は

$$
\langle u^2\rangle/v_r^2=\left(1\right)^{2} + 3=4
$$

である。

速さの四乗平均は

$$
\langle u^4\rangle/v_r^4=\left(1\right)^{4} + 10 \left(1\right)^{2} + 15=26
$$

である。

等方温度で置き換えた四乗平均は

$$
M_{\rm iso}=15 \left(\frac{4}{3}\right)^{2}=\frac{80}{3}
$$

である。

四乗に対する置換の相対誤差は

$$
\delta_4=-1 + \frac{\left(\frac{80}{3}\right)}{\left(26\right)}=\frac{1}{39}
$$

である。

温度置換の値の方が大きい。
:::

:::example[大きい平均流]
流速付近の確率と四乗平均を求める。

入力値を次のように置く。

$$
\begin{aligned}
r&=4 \\x&=4
\end{aligned}
$$

これらを順に代入する。

平均エネルギーを合わせた等方温度比は

$$
H=\frac{\left(4\right)^{2}}{3} + 1=\frac{19}{3}
$$

である。

規格化した速さ分布は

$$
v_r p(u)=\frac{\sqrt{2} \left(e^{\left(4\right)^{2}} - e^{- \left(4\right)^{2}}\right) e^{- \left(4\right)^{2}}}{2 \sqrt{\pi}}\simeq 0.39894
$$

である。

速さの二乗平均は

$$
\langle u^2\rangle/v_r^2=\left(4\right)^{2} + 3=19
$$

である。

速さの四乗平均は

$$
\langle u^4\rangle/v_r^4=\left(4\right)^{4} + 10 \left(4\right)^{2} + 15=431
$$

である。

等方温度で置き換えた四乗平均は

$$
M_{\rm iso}=15 \left(\frac{19}{3}\right)^{2}=\frac{1805}{3}
$$

である。

四乗に対する置換の相対誤差は

$$
\delta_4=-1 + \frac{\left(\frac{1805}{3}\right)}{\left(431\right)}=\frac{512}{1293}
$$

である。

速度分布は中心ゼロの Maxwell 分布と異なる。
:::

:::example[流れのある系の低速側]
平均流を同じにし、分布の低速側を比較する。

入力値を次のように置く。

$$
\begin{aligned}
r&=4 \\x&=1
\end{aligned}
$$

これらを順に代入する。

平均エネルギーを合わせた等方温度比は

$$
H=\frac{\left(4\right)^{2}}{3} + 1=\frac{19}{3}
$$

である。

規格化した速さ分布は

$$
v_r p(u)=\frac{\sqrt{2} \left(1\right) \left(e^{\left(1\right) \left(4\right)} - e^{- \left(1\right) \left(4\right)}\right) e^{- \frac{\left(1\right)^{2}}{2} - \frac{\left(4\right)^{2}}{2}}}{2 \sqrt{\pi} \left(4\right)}\simeq 0.0011076
$$

である。

速さの二乗平均は

$$
\langle u^2\rangle/v_r^2=\left(4\right)^{2} + 3=19
$$

である。

速さの四乗平均は

$$
\langle u^4\rangle/v_r^4=\left(4\right)^{4} + 10 \left(4\right)^{2} + 15=431
$$

である。

等方温度で置き換えた四乗平均は

$$
M_{\rm iso}=15 \left(\frac{19}{3}\right)^{2}=\frac{1805}{3}
$$

である。

四乗に対する置換の相対誤差は

$$
\delta_4=-1 + \frac{\left(\frac{1805}{3}\right)}{\left(431\right)}=\frac{512}{1293}
$$

である。

速さ分布の場所によっては数値が大幅に変わる。反応断面積の重みが要る。
:::


## Lawson 条件と点火

荷電した $\alpha$ 粒子はプラズマ中で減速して加熱に寄与できる。中性子は磁場で閉じ込められず、主に炉の周囲へ出る。ここでは一様な DT 燃料、$n=n_e=n_D+n_T$、$T_i=T_e=T$、等量混合、完全な $\alpha$ エネルギー堆積を仮定する。

蓄積熱エネルギー密度 $W=3nT$ に対し、放射を含む総損失を $W/\tau_E$ と定義する。外部加熱密度 $P_{\rm ext}$、$P_\alpha=n^2\langle\sigma u\rangle E_\alpha/4$ なら $\dot W=P_{\rm ext}+P_\alpha-W/\tau_E$ である。$P_{\rm ext}=0$ で平衡を保つ :term[点火]{en="ignition" yomi="てんか"} の境界は
$$
n\tau_E=\frac{12 T}{\langle\sigma u\rangle E_\alpha}.
\label{eq-plasma-lawson}
$$
実際に持続するには燃料と不純物の収支、平衡の温度安定性も必要である。式 :eqref[eq-plasma-lawson] はこの模型の :term[Lawson 条件]{en="Lawson criterion" yomi="ろーそんじょうけん"} であり、$nT\tau_E$ を比較する表現も使われる。

### 三重積を温度の関数として比較する

点火に必要な三重積 $nT\tau_E$ は、この一様な模型では $12T^2/[K(T)E_\alpha]$、$K(T)=\langle\sigma u\rangle$ である。従って極値の条件は $T K'(T)/K(T)=2$ となる。実際、三重積の対数を $\ln T$ で微分すると $2-TK'/K$ である。反応度の対数勾配が温度の上昇とともに二を上から下へ横切るなら、その交点で必要な三重積は局所的に最小になる。$n\tau_E$ を最小にするときの条件は、分子が $T$ なので対数勾配が一である。例えば狭い温度範囲で $K\propto T^3$ なら三重積は $T^{-1}$ に比例するが、$K\propto T$ の範囲では $T$ に比例する。具体的な温度を決めるには核反応断面積から得た $K(T)$ と、採用した損失模型を使う。

### 有限利得と放射

$P_f=n^2\langle\sigma u\rangle E_f/4$、$Q=P_f/P_{\rm ext}$ なら定常式から $n\tau_E=12T/[\langle\sigma u\rangle(E_\alpha+E_f/Q)]$ である。$Q=1$ は外部加熱と融合出力の等しさを表し、点火は $Q\to\infty$ の極限である。加熱器の電力や熱から電気への変換はこの $Q$ の定義に入っていない。

輸送の損失だけに $\tau_{\rm tr}$ を用いる慣例では、放射 $P_{\rm rad}$ を別に引いて $3nT/\tau_{\rm tr}=P_{\rm ext}+P_\alpha-P_{\rm rad}$ とする。非相対論的で光学的に薄い制動放射では $P_{\rm rad}\propto n_e\sum_i Z_i^2n_i\sqrt{T_e}$ が目安である。$P_\alpha<P_{\rm rad}$ なら輸送損失がゼロでも自己加熱だけでは釣り合わない。等号の場合、釣り合うには輸送損失もゼロでなければならない。

:::example[点火条件の尺度]
$T=10\,\mathrm{keV}$、反応度を例として $10^{-22}\,\mathrm{m^3\,s^{-1}}$、$E_\alpha=3500\,\mathrm{keV}$ とすれば $n\tau_E\simeq3.43\times10^{20}\,\mathrm{m^{-3}s}$ である。エネルギーの比には同じ keV 単位を使える。
:::

### 温度の小さな摂動

密度を固定した一様模型で $C_V=3n$、$L(T)=P_{\rm loss}(T)-P_\alpha(T)-P_{\rm ext}(T)$ とする。定常温度 $T_0$ の周りでは $C_V\partial_t\delta T=-L'(T_0)\delta T$ である。従って $L'>0$ はこの温度モードの減衰、$L'<0$ は成長を意味する。$\tau_E$ が温度や加熱で変わる場合も、損失の微分に含める。定常収支が等しいことだけでは安定性は分からない。

:::example[局所的な加熱の傾き]
$P_{\rm ext}=0$、$P_{\rm loss}\propto T$、近傍で $P_\alpha\propto T^2$ と仮定する。両者が $P_0$ で交わっても $L'=-P_0/T_0$ なので温度摂動は成長する。別の温度依存の損失や加熱制御があればこの符号は変わる。
:::

### 有限の自己加熱と密度の希釈

全イオン密度を基準にした等量混合の模型で、電子密度が燃料イオン密度と異なるとする。温度が両者で等しく、全てのイオン数密度を $n_i$、燃料密度を $n_f$、$\alpha$ の堆積割合を $f_\alpha$ とすると、熱エネルギーは $W=3(n_e+n_i)T/2$、加熱は $f_\alpha n_f^2\langle\sigma u\rangle E_\alpha/4$ である。したがって点火に必要な時間は
$$
\tau_E=\frac{6(n_e+n_i)T}{f_\alpha\,n_f^2\langle\sigma u\rangle E_\alpha}
$$
になる。堆積割合と燃料比は異なる要因であり、比べる炉ごとに何を一定にしたかを指定する。純 DT で $n_e=n_i=n_f$、$f_\alpha=1$ なら式 :eqref[eq-plasma-lawson] に戻る。

:::example[部分堆積]
同一の純 DT 密度と温度で $f_\alpha=0.8$ とすれば、必要な $\tau_E$ は完全堆積の $1/0.8=1.25$ 倍である。密度を変えて補う場合、同じ $f_\alpha$ が保たれるとは限らず、粒子の減速も確認する。
:::

### 自己加熱と制動放射を区別した利得

一様な等量 DT で融合電力密度を $P_f$、堆積する自己加熱の割合を $a=P_{\alpha,\rm dep}/P_f$、放射を $b=P_{\rm rad}/P_f$ とする。輸送にのみ $\tau_{\rm tr}$ を使うと定常収支は $3nT/\tau_{\rm tr}=P_f(a-b+1/Q)$ となる。$P_f=n^2 K E_f/4$ なので $n\tau_{\rm tr}KE_f/(12T)=1/(a-b+1/Q)$ である。この正の右辺が存在するためには括弧が正でなければならない。自己加熱で運転する極限では $a>b$ が必要である。$Q$ が有限なら $1/Q$ の外部入力を足せる。これらの比較は密度・温度を指定しているため、平衡点の立ち上げや安定性を保証しない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $a$ | 堆積する自己加熱の融合電力に対する比 | 1 |
| $b$ | 放射損失の融合電力に対する比 | 1 |
| $u$ | 外部加熱の融合電力に対する比 $1/Q$ | 1 |

放射前の規格化加熱を $H_0$ とする。自己加熱と外部加熱を足す。

$$
H_0=a + u
$$

自己加熱から放射を引いた値を $H_\alpha$ とする。点火極限の正値条件になる。

$$
H_\alpha=a - b
$$

規格化した輸送への残りを $H_{\rm tr}$ とする。定常時に輸送が運び去る分である。

$$
H_{\rm tr}=H_{0} - b
$$

規格化した必要閉じ込めを $C_{\rm tr}$ とする。$H_{\rm tr}>0$ の場合に意味がある。

$$
C_{\rm tr}=\frac{1}{H_{\rm tr}}
$$

放射前の熱源のうち自己加熱の比を $f_{\rm self}$ とする。融合電力全体に対する割合とは異なる。

$$
f_{\rm self}=\frac{a}{H_{0}}
$$

加熱に占める放射損失を $P_{\rm rad}/P_{\rm heat}$ とする。残りは輸送へ渡る。

$$
P_{\rm rad}/P_{\rm heat}=\frac{b}{H_{0}}
$$

:::example[放射の小さい点火模型]
外部加熱ゼロで輸送へ渡る電力を求める。

入力値を次のように置く。

$$
\begin{aligned}
a&=\frac{1}{5} \\b&=\frac{1}{50} \\u&=0
\end{aligned}
$$

これらを順に代入する。

放射前の規格化加熱は

$$
H_0=\left(0\right) + \left(\frac{1}{5}\right)=\frac{1}{5}
$$

である。

自己加熱から放射を引いた値は

$$
H_\alpha=- \left(\frac{1}{50}\right) + \left(\frac{1}{5}\right)=\frac{9}{50}
$$

である。

規格化した輸送への残りは

$$
H_{\rm tr}=- \left(\frac{1}{50}\right) + \left(\frac{1}{5}\right)=\frac{9}{50}
$$

である。

規格化した必要閉じ込めは

$$
C_{\rm tr}=\frac{1}{\left(\frac{9}{50}\right)}=\frac{50}{9}
$$

である。

放射前の熱源のうち自己加熱の比は

$$
f_{\rm self}=1
$$

である。

加熱に占める放射損失は

$$
P_{\rm rad}/P_{\rm heat}=\frac{\left(\frac{1}{50}\right)}{\left(\frac{1}{5}\right)}=\frac{1}{10}
$$

である。

放射がなければ必要な閉じ込めは五だが、この条件ではより長い。
:::

:::example[助熱のある定常点]
外部加熱を加えて必要な閉じ込めを比べる。

入力値を次のように置く。

$$
\begin{aligned}
a&=\frac{1}{5} \\b&=\frac{1}{50} \\u&=\frac{1}{10}
\end{aligned}
$$

これらを順に代入する。

放射前の規格化加熱は

$$
H_0=\left(\frac{1}{10}\right) + \left(\frac{1}{5}\right)=\frac{3}{10}
$$

である。

自己加熱から放射を引いた値は

$$
H_\alpha=- \left(\frac{1}{50}\right) + \left(\frac{1}{5}\right)=\frac{9}{50}
$$

である。

規格化した輸送への残りは

$$
H_{\rm tr}=- \left(\frac{1}{50}\right) + \left(\frac{3}{10}\right)=\frac{7}{25}
$$

である。

規格化した必要閉じ込めは

$$
C_{\rm tr}=\frac{1}{\left(\frac{7}{25}\right)}=\frac{25}{7}
$$

である。

放射前の熱源のうち自己加熱の比は

$$
f_{\rm self}=\frac{\left(\frac{1}{5}\right)}{\left(\frac{3}{10}\right)}=\frac{2}{3}
$$

である。

加熱に占める放射損失は

$$
P_{\rm rad}/P_{\rm heat}=\frac{\left(\frac{1}{50}\right)}{\left(\frac{3}{10}\right)}=\frac{1}{15}
$$

である。

与えた温度で必要な閉じ込めの下限が緩和される。
:::

:::example[放射が自己加熱を超える場合]
定常状態と点火の条件を区別する。

入力値を次のように置く。

$$
\begin{aligned}
a&=\frac{1}{10} \\b&=\frac{1}{5} \\u&=\frac{1}{4}
\end{aligned}
$$

これらを順に代入する。

放射前の規格化加熱は

$$
H_0=\left(\frac{1}{10}\right) + \left(\frac{1}{4}\right)=\frac{7}{20}
$$

である。

自己加熱から放射を引いた値は

$$
H_\alpha=\left(\frac{1}{10}\right) - \left(\frac{1}{5}\right)=- \frac{1}{10}
$$

である。

規格化した輸送への残りは

$$
H_{\rm tr}=- \left(\frac{1}{5}\right) + \left(\frac{7}{20}\right)=\frac{3}{20}
$$

である。

規格化した必要閉じ込めは

$$
C_{\rm tr}=\frac{1}{\left(\frac{3}{20}\right)}=\frac{20}{3}
$$

である。

放射前の熱源のうち自己加熱の比は

$$
f_{\rm self}=\frac{\left(\frac{1}{10}\right)}{\left(\frac{7}{20}\right)}=\frac{2}{7}
$$

である。

加熱に占める放射損失は

$$
P_{\rm rad}/P_{\rm heat}=\frac{\left(\frac{1}{5}\right)}{\left(\frac{7}{20}\right)}=\frac{4}{7}
$$

である。

定常値は正でも、入力をゼロへ下げると点火の条件を満たさない。
:::

:::example[放射で余裕をほとんど失う場合]
放射差に対する閉じ込めの感度を求める。

入力値を次のように置く。

$$
\begin{aligned}
a&=\frac{1}{5} \\b&=\frac{19}{100} \\u&=0
\end{aligned}
$$

これらを順に代入する。

放射前の規格化加熱は

$$
H_0=\left(0\right) + \left(\frac{1}{5}\right)=\frac{1}{5}
$$

である。

自己加熱から放射を引いた値は

$$
H_\alpha=- \left(\frac{19}{100}\right) + \left(\frac{1}{5}\right)=\frac{1}{100}
$$

である。

規格化した輸送への残りは

$$
H_{\rm tr}=- \left(\frac{19}{100}\right) + \left(\frac{1}{5}\right)=\frac{1}{100}
$$

である。

規格化した必要閉じ込めは

$$
C_{\rm tr}=\frac{1}{\left(\frac{1}{100}\right)}=100
$$

である。

放射前の熱源のうち自己加熱の比は

$$
f_{\rm self}=1
$$

である。

加熱に占める放射損失は

$$
P_{\rm rad}/P_{\rm heat}=\frac{\left(\frac{19}{100}\right)}{\left(\frac{1}{5}\right)}=\frac{19}{20}
$$

である。

残る電力が小さいので必要な閉じ込めが大きくなる。
:::


## トカマクの閉じ込めと閉じ込め時間

磁場閉じ込めの時間には、粒子数を保持する時間と熱エネルギーを保持する時間がある。後者を :term[エネルギー閉じ込め時間]{en="energy confinement time" yomi="えねるぎーとじこめじかん"} と呼ぶ。体積積分した熱エネルギー $W_p$ と総加熱 $P_{\rm heat}$ を測り、$P_{\rm loss}=P_{\rm heat}-\dot W_p>0$ のとき $\tau_E=W_p/P_{\rm loss}$ と評価する。定常でのみ $P_{\rm loss}=P_{\rm heat}$ である。

円筒近似で熱拡散率 $\chi$ が一定、辺で温度摂動がゼロとする。$\partial_t T=\chi r^{-1}\partial_r(r\partial_r T)$ の分離解は $T=e^{-\chi k^2t}J_0(kr)$ である。軸で発散するもう一つの Bessel 解を捨て、$J_0(ka)=0$ を課す。最低固有値で $ka\simeq2.405$ なので減衰時間は $a^2/(\chi\,2.405^2)$ である。実際の装置では輸送係数は一定でなく、加熱分布、周辺の勾配、乱流状態で変わる。経験的スケーリングは有限の運転領域での近似であり、拡散模型の指数と同一視できない。

:::example[過渡の閉じ込め評価]
$W_p=10\,\mathrm{MJ}$、$P_{\rm heat}=8\,\mathrm{MW}$、$\dot W_p=3\,\mathrm{MW}$ なら損失は $5\,\mathrm{MW}$、$\tau_E=2\,\mathrm s$ である。
:::

:::example[拡散時間との比較]
$\chi=1\,\mathrm{m^2s^{-1}}$、$a=1\,\mathrm m$ なら、円筒の最低モードの時間は $1/2.405^2\simeq0.173\,\mathrm s$ である。同じ $\chi$ で半径を二倍にすると時間は四倍になる。実機で $\chi$ が同じ値を保つとは限らないので、この倍率は一定係数の模型の結果である。
:::

### 閉じ込めの経験則と入力電力

トカマクでは外周の輸送が変わると、同じ磁場と密度でも熱の蓄積が変わる。低閉じ込めの状態を L モード、周辺に輸送障壁が形成される高閉じ込めの状態を H モードという。周辺で温度や密度の勾配を支えられるようになれば、内部の圧力分布全体も変化する。従って異なる運転状態を一つの拡散率の変更だけで比較するのは難しい。閉じ込めの輸送と経験則の背景は :cite[chen-plasma] を参照する。

同じ運転状態のデータで、入力のほかの変数を固定した熱閉じ込めが $\tau_E=\tau_0(P_{\rm loss}/P_0)^{-\beta}$ と近似できるとする。定常なら蓄積熱は
$$
W_p=P_{\rm loss}\tau_E
=P_0\tau_0(P_{\rm loss}/P_0)^{1-\beta}.
$$
$0<\beta<1$ のとき、加熱を増して蓄積熱が増えても閉じ込め時間は短くなる。この二つの量を分ければ、装置の加熱を変えたときに何を観測すべきかが分かる。例えば $\beta=2/3$ の仮定で損失と釣り合う加熱を八倍にすると、$\tau_E$ は四分の一、$W_p$ は二倍である。指数は経験則に固有であり、別の運転状態に外挿する場合はデータの範囲や放射の扱いを確認する。

### 蓄積熱の過渡から閉じ込めを同定する

単一の熱容量に対して加熱を $P$ に変え、一定の閉じ込め $\tau$ を仮定すると $\dot W=P-W/\tau$ である。積分因子を用いて $W(t)=P\tau+(W_0-P\tau)e^{-t/\tau}$ となる。$\dot W$ を測ったときの瞬時比 $W/(P-\dot W)$ はすべての時刻で $\tau$ である。$W/P$ を用いると過渡の間はこの値を誤る。装置で推定される値は放射を総損失に含めるか、熱以外の蓄積を含めるかでも変わる。ここでは一定熱容量で等温の一領域模型のみを計算する。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $W_0$ | 初期の蓄積熱 | $\mathrm{MJ}$ |
| $P$ | 変更後一定の加熱 | $\mathrm{MW}$ |
| $\tau$ | 仮定した閉じ込め時間 | $\mathrm s$ |
| $s$ | 加熱変更後の時間 $t/\tau$ | 1 |

新しい定常熱エネルギーを $W_\infty$ とする。損失が加熱と一致する。

$$
W_\infty=P \tau
$$

過渡時の蓄積熱を $W$ とする。初期と定常との差が指数で減る。

$$
W=W_{\infty} + \left(W_{0} - W_{\infty}\right) e^{- s}
$$

過渡時の損失を $P_{\rm loss}$ とする。熱エネルギーを閉じ込め時間で割る。

$$
P_{\rm loss}=\frac{W}{\tau}
$$

熱エネルギーの増加率を $\dot W$ とする。加熱と損失を別々に計上する。

$$
\dot W=P - P_{\rm loss}
$$

加熱を損失としたときの形式時間を $W/P$ とする。過渡の増減を無視した値である。

$$
W/P=\frac{W}{P}
$$

収支から求めた閉じ込めを $W/(P-\dot W)$ とする。蓄積の変化も測定に入れる。

$$
W/(P-\dot W)=\frac{W}{P - \dot W}
$$

:::example[加熱を立ち上げる過渡]
放出されずに蓄積する電力を求める。

入力値を次のように置く。

$$
\begin{aligned}
W_0&=2\,\mathrm{MJ} \\P&=5\,\mathrm{MW} \\\tau&=2\,\mathrm s \\s&=\frac{1}{2}
\end{aligned}
$$

これらを順に代入する。

新しい定常熱エネルギーは

$$
W_\infty=\left(2\right) \left(5\right)=10\,\mathrm{MJ}
$$

である。

過渡時の蓄積熱は

$$
W=\left(10\right) + \left(- \left(10\right) + \left(2\right)\right) e^{- \left(\frac{1}{2}\right)}\simeq 5.1478\,\mathrm{MJ}
$$

である。

過渡時の損失は

$$
P_{\rm loss}=\frac{\left(10 - \frac{8}{e^{\frac{1}{2}}}\right)}{\left(2\right)}\simeq 2.5739\,\mathrm{MW}
$$

である。

熱エネルギーの増加率は

$$
\dot W=- \left(5 - \frac{4}{e^{\frac{1}{2}}}\right) + \left(5\right)\simeq 2.4261\,\mathrm{MW}
$$

である。

加熱を損失としたときの形式時間は

$$
W/P=\frac{\left(10 - \frac{8}{e^{\frac{1}{2}}}\right)}{\left(5\right)}\simeq 1.0296\,\mathrm s
$$

である。

収支から求めた閉じ込めは

$$
W/(P-\dot W)=\frac{\left(10 - \frac{8}{e^{\frac{1}{2}}}\right)}{\left(5\right) - \left(\frac{4}{e^{\frac{1}{2}}}\right)}\simeq 2\,\mathrm s
$$

である。

$W/P$ は閉じ込め時間を過小評価する。
:::

:::example[加熱を下げた過渡]
失う電力が現在の加熱を上回るか調べる。

入力値を次のように置く。

$$
\begin{aligned}
W_0&=10\,\mathrm{MJ} \\P&=2\,\mathrm{MW} \\\tau&=2\,\mathrm s \\s&=\frac{1}{2}
\end{aligned}
$$

これらを順に代入する。

新しい定常熱エネルギーは

$$
W_\infty=\left(2\right)^{2}=4\,\mathrm{MJ}
$$

である。

過渡時の蓄積熱は

$$
W=\left(4\right) + \left(\left(10\right) - \left(4\right)\right) e^{- \left(\frac{1}{2}\right)}\simeq 7.6392\,\mathrm{MJ}
$$

である。

過渡時の損失は

$$
P_{\rm loss}=\frac{\left(\frac{6}{e^{\frac{1}{2}}} + 4\right)}{\left(2\right)}\simeq 3.8196\,\mathrm{MW}
$$

である。

熱エネルギーの増加率は

$$
\dot W=\left(2\right) - \left(\frac{3}{e^{\frac{1}{2}}} + 2\right)\simeq -1.8196\,\mathrm{MW}
$$

である。

加熱を損失としたときの形式時間は

$$
W/P=\frac{\left(\frac{6}{e^{\frac{1}{2}}} + 4\right)}{\left(2\right)}\simeq 3.8196\,\mathrm s
$$

である。

収支から求めた閉じ込めは

$$
W/(P-\dot W)=\frac{\left(\frac{6}{e^{\frac{1}{2}}} + 4\right)}{- \left(- \frac{3}{e^{\frac{1}{2}}}\right) + \left(2\right)}\simeq 2\,\mathrm s
$$

である。

蓄積熱も放出されるため、$W/P$ は閉じ込めを過大に見積もる。
:::

:::example[初めから定常の場合]
二つの評価法が一致する条件を確認する。

入力値を次のように置く。

$$
\begin{aligned}
W_0&=10\,\mathrm{MJ} \\P&=5\,\mathrm{MW} \\\tau&=2\,\mathrm s \\s&=1
\end{aligned}
$$

これらを順に代入する。

新しい定常熱エネルギーは

$$
W_\infty=\left(2\right) \left(5\right)=10\,\mathrm{MJ}
$$

である。

過渡時の蓄積熱は

$$
W=\left(10\right)=10\,\mathrm{MJ}
$$

である。

過渡時の損失は

$$
P_{\rm loss}=\frac{\left(10\right)}{\left(2\right)}=5\,\mathrm{MW}
$$

である。

熱エネルギーの増加率は

$$
\dot W=0\,\mathrm{MW}
$$

である。

加熱を損失としたときの形式時間は

$$
W/P=\frac{\left(10\right)}{\left(5\right)}=2\,\mathrm s
$$

である。

収支から求めた閉じ込めは

$$
W/(P-\dot W)=\frac{\left(10\right)}{- \left(0\right) + \left(5\right)}=2\,\mathrm s
$$

である。

$\dot W$ はゼロなので両者が一致する。
:::

:::example[十分長く待った場合]
定常の近似に残る誤差を求める。

入力値を次のように置く。

$$
\begin{aligned}
W_0&=2\,\mathrm{MJ} \\P&=5\,\mathrm{MW} \\\tau&=2\,\mathrm s \\s&=5
\end{aligned}
$$

これらを順に代入する。

新しい定常熱エネルギーは

$$
W_\infty=\left(2\right) \left(5\right)=10\,\mathrm{MJ}
$$

である。

過渡時の蓄積熱は

$$
W=\left(10\right) + \left(- \left(10\right) + \left(2\right)\right) e^{- \left(5\right)}\simeq 9.9461\,\mathrm{MJ}
$$

である。

過渡時の損失は

$$
P_{\rm loss}=\frac{\left(10 - \frac{8}{e^{5}}\right)}{\left(2\right)}\simeq 4.973\,\mathrm{MW}
$$

である。

熱エネルギーの増加率は

$$
\dot W=- \left(5 - \frac{4}{e^{5}}\right) + \left(5\right)\simeq 0.026952\,\mathrm{MW}
$$

である。

加熱を損失としたときの形式時間は

$$
W/P=\frac{\left(10 - \frac{8}{e^{5}}\right)}{\left(5\right)}\simeq 1.9892\,\mathrm s
$$

である。

収支から求めた閉じ込めは

$$
W/(P-\dot W)=\frac{\left(10 - \frac{8}{e^{5}}\right)}{\left(5\right) - \left(\frac{4}{e^{5}}\right)}\simeq 2\,\mathrm s
$$

である。

加熱変更から数倍の閉じ込め時間がたつと過渡の寄与は小さくなる。
:::


## プラズマの加熱法

外部加熱では、注入したエネルギーのうち、熱粒子へ移った量を収支に使う。Ohm 加熱は :eqref[eq-plasma-ohm] の抵抗項から $P_{\rm Ohm}=\eta j^2$ である。高温で抵抗率が下がるため、同じ電流密度で加熱を増やすことは難しくなる。

中性粒子ビームは電荷がないときに磁気面を横切れる。プラズマ内で電離すると速いイオンの軌道になり、衝突で熱粒子にエネルギーを渡す。均一な標的で有効電離断面積 $\sigma_{\rm ion}$ を一定とすれば、距離 $x$ での生存粒子束は $\Gamma(x)=\Gamma(0)e^{-n\sigma_{\rm ion}x}$ である。微小区間で失う割合が $n\sigma_{\rm ion}\dd x$ であることを積分した。外に抜ける成分や壁に当たる軌道は熱に数えない。

高周波加熱では電子またはイオンのサイクロトロン共鳴、混成波や平行共鳴を選ぶ。局所的な堆積には偏波、アクセスできる伝播分枝、吸収の勾配が関わる。周波数が裸の旋回周波数に等しくても、波が領域へ伝播できるとは限らない。波の伝播と吸収は [誘電テンソルと波動伝播](/physics/plasma-physics/dielectric-tensor/) の手順で分けて評価する。

:::example[未吸収の比率]
単純なビーム模型で $n\sigma_{\rm ion}L=2$ なら、$L$ を通っても電離していない割合は $e^{-2}\simeq0.135$ である。
:::

:::example[高周波で選ぶ共鳴]
電子サイクロトロンの基本周波数は $f_{ce}=eB/(2\pi m_e)$ である。$B=2\,\mathrm T$ なら約 $56\,\mathrm{GHz}$ である。平行速度がある電子ではここから Doppler ずれがあり、磁場が空間で変われば周波数が共鳴する場所も変わる。
:::

### 中性粒子ビーム源の空間的なエネルギー収支

幅一定のビームを一様に吸収する標的を考える。$P'(x)=-\kappa P(x)$、$\kappa=n\sigma$ とすると $P(x)=P_0e^{-\kappa x}$ である。区間 $[x_1,x_2]$ でイオン化する粒子が持ち込む電力は $P_0(e^{-\kappa x_1}-e^{-\kappa x_2})$ になる。その位置での熱入力を求めるには、イオン化した粒子が局所で熱化するという近似も要る。捕捉された高速イオンのうち加熱を担う割合を $f_{\rm th}$ とし、残りが軌道損失で逃げる単純模型では、熱電力はこの差に $f_{\rm th}$ を掛ける。以下では長さ $L$ の層の中央で二分し、入口側と奥側の差を計算する。実際の分布が不均一なら $\kappa x$ を光路上の積分で置き換える。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $P_0$ | 入射ビーム電力 | $\mathrm{MW}$ |
| $\tau$ | 全層の光学厚さ $\kappa L$ | 1 |
| $f_{\rm th}$ | イオン化後に熱化する割合 | 1 |

通り抜ける中性ビームの電力を $P_{\rm pass}$ とする。イオン化していない成分である。

$$
P_{\rm pass}=P_{0} e^{- \tau}
$$

前半で失われたビーム電力を $P_{\rm front}$ とする。入口と中点の差である。

$$
P_{\rm front}=P_{0} \left(1 - e^{- \frac{\tau}{2}}\right)
$$

後半で失われたビーム電力を $P_{\rm back}$ とする。中点と出口の差である。

$$
P_{\rm back}=P_{0} e^{- \frac{\tau}{2}} - P_{\rm pass}
$$

全層で熱化する電力を $P_{\rm heat}$ とする。速度分布への遅延を省いた値である。

$$
P_{\rm heat}=f_{\rm th} \left(P_{\rm back} + P_{\rm front}\right)
$$

イオン化後の軌道損失を $P_{\rm orbit}$ とする。中性のままの透過とは異なる。

$$
P_{\rm orbit}=\left(1 - f_{\rm th}\right) \left(P_{\rm back} + P_{\rm front}\right)
$$

前半に対する後半の堆積比を $P_{\rm back}/P_{\rm front}$ とする。十分薄ければ一に近づく。

$$
P_{\rm back}/P_{\rm front}=\frac{P_{\rm back}}{P_{\rm front}}
$$

:::example[吸収の弱い標的]
奥まで届くビームで加熱分布を求める。

入力値を次のように置く。

$$
\begin{aligned}
P_0&=10\,\mathrm{MW} \\\tau&=\frac{1}{10} \\f_{\rm th}&=1
\end{aligned}
$$

これらを順に代入する。

通り抜ける中性ビームの電力は

$$
P_{\rm pass}=\left(10\right) e^{- \left(\frac{1}{10}\right)}\simeq 9.0484\,\mathrm{MW}
$$

である。

前半で失われたビーム電力は

$$
P_{\rm front}=\left(10\right) \left(1 - e^{- \frac{\left(\frac{1}{10}\right)}{2}}\right)\simeq 0.48771\,\mathrm{MW}
$$

である。

後半で失われたビーム電力は

$$
P_{\rm back}=\left(10\right) e^{- \frac{\left(\frac{1}{10}\right)}{2}} - \left(\frac{10}{e^{\frac{1}{10}}}\right)\simeq 0.46392\,\mathrm{MW}
$$

である。

全層で熱化する電力は

$$
P_{\rm heat}=\left(1\right) \left(\left(- \frac{10}{e^{\frac{1}{10}}} + \frac{10}{e^{\frac{1}{20}}}\right) + \left(10 - \frac{10}{e^{\frac{1}{20}}}\right)\right)\simeq 0.95163\,\mathrm{MW}
$$

である。

イオン化後の軌道損失は

$$
P_{\rm orbit}=\left(1 - \left(1\right)\right) \left(\left(- \frac{10}{e^{\frac{1}{10}}} + \frac{10}{e^{\frac{1}{20}}}\right) + \left(10 - \frac{10}{e^{\frac{1}{20}}}\right)\right)=0\,\mathrm{MW}
$$

である。

前半に対する後半の堆積比は

$$
P_{\rm back}/P_{\rm front}=\frac{\left(- \frac{10}{e^{\frac{1}{10}}} + \frac{10}{e^{\frac{1}{20}}}\right)}{\left(10 - \frac{10}{e^{\frac{1}{20}}}\right)}\simeq 0.95123
$$

である。

前後の差は小さいが、多くの電力が通り抜ける。
:::

:::example[中程度の吸収]
深部と入口の発熱を比べる。

入力値を次のように置く。

$$
\begin{aligned}
P_0&=10\,\mathrm{MW} \\\tau&=2 \\f_{\rm th}&=1
\end{aligned}
$$

これらを順に代入する。

通り抜ける中性ビームの電力は

$$
P_{\rm pass}=\left(10\right) e^{- \left(2\right)}\simeq 1.3534\,\mathrm{MW}
$$

である。

前半で失われたビーム電力は

$$
P_{\rm front}=\left(10\right) \left(1 - e^{- \frac{\left(2\right)}{2}}\right)\simeq 6.3212\,\mathrm{MW}
$$

である。

後半で失われたビーム電力は

$$
P_{\rm back}=\left(10\right) e^{- \frac{\left(2\right)}{2}} - \left(\frac{10}{e^{2}}\right)\simeq 2.3254\,\mathrm{MW}
$$

である。

全層で熱化する電力は

$$
P_{\rm heat}=\left(1\right) \left(\left(- \frac{10}{e^{2}} + \frac{10}{e}\right) + \left(10 - \frac{10}{e}\right)\right)\simeq 8.6466\,\mathrm{MW}
$$

である。

イオン化後の軌道損失は

$$
P_{\rm orbit}=\left(1 - \left(1\right)\right) \left(\left(- \frac{10}{e^{2}} + \frac{10}{e}\right) + \left(10 - \frac{10}{e}\right)\right)=0\,\mathrm{MW}
$$

である。

前半に対する後半の堆積比は

$$
P_{\rm back}/P_{\rm front}=\frac{\left(- \frac{10}{e^{2}} + \frac{10}{e}\right)}{\left(10 - \frac{10}{e}\right)}\simeq 0.36788
$$

である。

ビームの減衰によって前半の加熱が強い。
:::

:::example[試験粒子の軌道損失を加える]
イオン化量を保ったまま熱化を減らす。

入力値を次のように置く。

$$
\begin{aligned}
P_0&=10\,\mathrm{MW} \\\tau&=2 \\f_{\rm th}&=\frac{3}{4}
\end{aligned}
$$

これらを順に代入する。

通り抜ける中性ビームの電力は

$$
P_{\rm pass}=\left(10\right) e^{- \left(2\right)}\simeq 1.3534\,\mathrm{MW}
$$

である。

前半で失われたビーム電力は

$$
P_{\rm front}=\left(10\right) \left(1 - e^{- \frac{\left(2\right)}{2}}\right)\simeq 6.3212\,\mathrm{MW}
$$

である。

後半で失われたビーム電力は

$$
P_{\rm back}=\left(10\right) e^{- \frac{\left(2\right)}{2}} - \left(\frac{10}{e^{2}}\right)\simeq 2.3254\,\mathrm{MW}
$$

である。

全層で熱化する電力は

$$
P_{\rm heat}=\left(\frac{3}{4}\right) \left(\left(- \frac{10}{e^{2}} + \frac{10}{e}\right) + \left(10 - \frac{10}{e}\right)\right)\simeq 6.485\,\mathrm{MW}
$$

である。

イオン化後の軌道損失は

$$
P_{\rm orbit}=\left(1 - \left(\frac{3}{4}\right)\right) \left(\left(- \frac{10}{e^{2}} + \frac{10}{e}\right) + \left(10 - \frac{10}{e}\right)\right)\simeq 2.1617\,\mathrm{MW}
$$

である。

前半に対する後半の堆積比は

$$
P_{\rm back}/P_{\rm front}=\frac{\left(- \frac{10}{e^{2}} + \frac{10}{e}\right)}{\left(10 - \frac{10}{e}\right)}\simeq 0.36788
$$

である。

透過電力は同じで、失われる内訳が変わる。
:::

:::example[吸収の大きい標的]
全体の吸収と奥への到達を区別する。

入力値を次のように置く。

$$
\begin{aligned}
P_0&=10\,\mathrm{MW} \\\tau&=8 \\f_{\rm th}&=\frac{3}{4}
\end{aligned}
$$

これらを順に代入する。

通り抜ける中性ビームの電力は

$$
P_{\rm pass}=\left(10\right) e^{- \left(8\right)}\simeq 0.0033546\,\mathrm{MW}
$$

である。

前半で失われたビーム電力は

$$
P_{\rm front}=\left(10\right) \left(1 - e^{- \frac{\left(8\right)}{2}}\right)\simeq 9.8168\,\mathrm{MW}
$$

である。

後半で失われたビーム電力は

$$
P_{\rm back}=\left(10\right) e^{- \frac{\left(8\right)}{2}} - \left(\frac{10}{e^{8}}\right)\simeq 0.1798\,\mathrm{MW}
$$

である。

全層で熱化する電力は

$$
P_{\rm heat}=\left(\frac{3}{4}\right) \left(\left(- \frac{10}{e^{8}} + \frac{10}{e^{4}}\right) + \left(10 - \frac{10}{e^{4}}\right)\right)\simeq 7.4975\,\mathrm{MW}
$$

である。

イオン化後の軌道損失は

$$
P_{\rm orbit}=\left(1 - \left(\frac{3}{4}\right)\right) \left(\left(- \frac{10}{e^{8}} + \frac{10}{e^{4}}\right) + \left(10 - \frac{10}{e^{4}}\right)\right)\simeq 2.4992\,\mathrm{MW}
$$

である。

前半に対する後半の堆積比は

$$
P_{\rm back}/P_{\rm front}=\frac{\left(- \frac{10}{e^{8}} + \frac{10}{e^{4}}\right)}{\left(10 - \frac{10}{e^{4}}\right)}\simeq 0.018316
$$

である。

ほぼ全てがイオン化しても均一には加熱されない。
:::


## 慣性閉じ込め核融合

慣性閉じ込めでは閉じ込め時間が流体の膨張時間で決まる。半径 $R$ の一様な球で音速 $c_s$ を代表速度とすると、$\tau_{\rm dis}\sim R/c_s$ であり、粒子密度 $n=\rho/\overline m$ に対して $n\tau_{\rm dis}\sim\rho R/(\overline m c_s)$ となる。ここで $\overline m$ は燃料イオンの平均質量である。従って :term[面密度]{en="areal density" yomi="めんみつど"} $\rho R$ は、密度と膨張までの時間の積を表す尺度になる。荷電生成物の停止は経路に沿う物質量にも依存する。面密度が同じでも温度や密度が異なれば、停止の詳細は変わる。

慣性閉じ込めの実験では、短時間で供給されたエネルギーで燃料を圧縮し、その慣性が膨張を遅らせる間に反応が進む。直接駆動では光が標的に直接作用し、間接駆動では一度放射に変換したエネルギーで駆動する。投入エネルギーの定義はこの方式と計測面に依存する。標的へ送ったエネルギーを $E_{\rm drive}$、燃料が受け取ったエネルギーを $E_{\rm abs}$ とし、$E_{\rm abs}/E_{\rm drive}=\eta_c$ と書けば、融合エネルギー $E_f^{\rm pulse}$ に対し
$$
\frac{E_f^{\rm pulse}}{E_{\rm drive}}=\eta_c\frac{E_f^{\rm pulse}}{E_{\rm abs}}.
$$
たとえば $\eta_c=1/4$ と仮定すれば、燃料への入力量を分母にした利得が $4$ のときに標的へ送ったエネルギーとの比が $1$ になる。装置に供給する電力量で評価する場合には、この比にさらに駆動源の効率が掛かる。

### 単純な燃焼率

等量の燃料が体積一定で混合し、反応度 $K$ が一定、燃料補給がない模型では $\dot n=-2R=-K n^2/2$ である。$n=n_0/(1+n_0Kt/2)$ を得る。時刻 $t$ の消費割合は $n_0Kt/(2+n_0Kt)$ である。実際の膨張では体積と温度が変わり、$\alpha$ 粒子の加熱と放射も連成する。総合的な燃焼と流体の非一様性を扱うには放射流体のモデルが必要であり、閉じ込めの概説は :cite[chen-plasma] を参照する。

:::example[低燃焼率極限]
$n_0Kt/2=0.1$ なら消費割合は $1/11$ である。反応率を初期値に固定すると $0.1$ となるが、燃料の減少を追うと小さくなる。
:::

:::example[膨張時間の相似比較]
温度を保ち、二つの一様膨張模型で $\overline m,c_s$ が同じとする。片方の密度が十倍で半径が十分の一なら、面密度と $n\tau_{\rm dis}$ は同じだが、膨張時間は十分の一である。面密度の一致から時間や総反応数の一致を結論することはできない。
:::

### 膨張する球の燃料密度を追う

半径が $R(t)=R_0(1+t/t_0)$ と変わる一様球を考える。温度を一様かつ一定とし、したがって熱反応度も一定とする理想化を用い、現実の断熱冷却は別に扱う。全燃料数 $N$ の収支は等量混合なら $\dot N=-K N^2/[2V(t)]$ である。$K$ は反応度、$V=V_0(1+t/t_0)^3$ である。$\partial_t(1/N)=K/(2V)$ を積分して
$$
\frac{N(t)}{N_0}=
\left[1+\frac{K n_0 t_0}{4}
\left(1-\frac{1}{(1+t/t_0)^2}\right)\right]^{-1}
$$
を得る。密度比にはさらに体積比の逆数が掛かる。無限時間でも積分された反応機会が有限になるのは、体積の増加が三乗だからである。定温を仮定した運動学的な模型であり、自己加熱や温度・圧縮の設計を与えない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $A$ | 反応の無次元係数 $K n_0 t_0/4$ | 1 |
| $x$ | 膨張時間 $t/t_0$ | 1 |

体積比を $V/V_0$ とする。球の半径比の三乗である。

$$
V/V_0=\left(x + 1\right)^{3}
$$

累積反応機会を $B$ とする。体積の逆数の積分を表す。

$$
B=A \left(1 - \frac{1}{\left(x + 1\right)^{2}}\right)
$$

全燃料数の残存率を $N/N_0$ とする。定積の反応式と同様に逆数を解く。

$$
N/N_0=\frac{1}{B + 1}
$$

密度の残存率を $n/n_0$ とする。反応による減少と希釈を両方含める。

$$
n/n_0=\frac{N/N_{0}}{V/V_{0}}
$$

初期燃料数のうち反応した割合を $f_{\rm burn}$ とする。密度比だけからは得られない。

$$
f_{\rm burn}=1 - N/N_{0}
$$

この定温模型の最終割合を $f_{\rm burn,\infty}$ とする。無限時間の幾何的積分から求める。

$$
f_{\rm burn,\infty}=\frac{A}{A + 1}
$$

:::example[膨張の早い段階]
密度低下と燃料の消費を比較する。

入力値を次のように置く。

$$
\begin{aligned}
A&=\frac{1}{10} \\x&=\frac{1}{10}
\end{aligned}
$$

これらを順に代入する。

体積比は

$$
V/V_0=\left(\left(\frac{1}{10}\right) + 1\right)^{3}=\frac{1331}{1000}
$$

である。

累積反応機会は

$$
B=\left(\frac{1}{10}\right) \left(1 - \frac{1}{\left(\left(\frac{1}{10}\right) + 1\right)^{2}}\right)=\frac{21}{1210}
$$

である。

全燃料数の残存率は

$$
N/N_0=\frac{1}{\left(\frac{21}{1210}\right) + 1}=\frac{1210}{1231}
$$

である。

密度の残存率は

$$
n/n_0=\frac{\left(\frac{1210}{1231}\right)}{\left(\frac{1331}{1000}\right)}=\frac{10000}{13541}
$$

である。

初期燃料数のうち反応した割合は

$$
f_{\rm burn}=1 - \left(\frac{1210}{1231}\right)=\frac{21}{1231}
$$

である。

この定温模型の最終割合は

$$
f_{\rm burn,\infty}=\frac{\left(\frac{1}{10}\right)}{\left(\frac{1}{10}\right) + 1}=\frac{1}{11}
$$

である。

密度の減少の大半は体積の拡大による。
:::

:::example[一膨張時間後]
半径が二倍になった場合を計算する。

入力値を次のように置く。

$$
\begin{aligned}
A&=\frac{1}{10} \\x&=1
\end{aligned}
$$

これらを順に代入する。

体積比は

$$
V/V_0=\left(\left(1\right) + 1\right)^{3}=8
$$

である。

累積反応機会は

$$
B=\left(\frac{1}{10}\right) \left(1 - \frac{1}{\left(\left(1\right) + 1\right)^{2}}\right)=\frac{3}{40}
$$

である。

全燃料数の残存率は

$$
N/N_0=\frac{1}{\left(\frac{3}{40}\right) + 1}=\frac{40}{43}
$$

である。

密度の残存率は

$$
n/n_0=\frac{\left(\frac{40}{43}\right)}{\left(8\right)}=\frac{5}{43}
$$

である。

初期燃料数のうち反応した割合は

$$
f_{\rm burn}=1 - \left(\frac{40}{43}\right)=\frac{3}{43}
$$

である。

この定温模型の最終割合は

$$
f_{\rm burn,\infty}=\frac{\left(\frac{1}{10}\right)}{\left(\frac{1}{10}\right) + 1}=\frac{1}{11}
$$

である。

希釈で反応が進みにくくなり、多くの燃料が残る。
:::

:::example[時間を長くとる]
最終割合への接近を調べる。

入力値を次のように置く。

$$
\begin{aligned}
A&=\frac{1}{10} \\x&=10
\end{aligned}
$$

これらを順に代入する。

体積比は

$$
V/V_0=\left(\left(10\right) + 1\right)^{3}=1331
$$

である。

累積反応機会は

$$
B=\left(\frac{1}{10}\right) \left(1 - \frac{1}{\left(\left(10\right) + 1\right)^{2}}\right)=\frac{12}{121}
$$

である。

全燃料数の残存率は

$$
N/N_0=\frac{1}{\left(\frac{12}{121}\right) + 1}=\frac{121}{133}
$$

である。

密度の残存率は

$$
n/n_0=\frac{\left(\frac{121}{133}\right)}{\left(1331\right)}=\frac{1}{1463}
$$

である。

初期燃料数のうち反応した割合は

$$
f_{\rm burn}=1 - \left(\frac{121}{133}\right)=\frac{12}{133}
$$

である。

この定温模型の最終割合は

$$
f_{\rm burn,\infty}=\frac{\left(\frac{1}{10}\right)}{\left(\frac{1}{10}\right) + 1}=\frac{1}{11}
$$

である。

密度は大きく落ちるが、全燃料数は有限の値に近づく。
:::

:::example[反応度を変えた形式的な比較]
同じ体積変化で反応項を大きくする。

入力値を次のように置く。

$$
\begin{aligned}
A&=1 \\x&=1
\end{aligned}
$$

これらを順に代入する。

体積比は

$$
V/V_0=\left(\left(1\right) + 1\right)^{3}=8
$$

である。

累積反応機会は

$$
B=\left(1\right) \left(1 - \frac{1}{\left(\left(1\right) + 1\right)^{2}}\right)=\frac{3}{4}
$$

である。

全燃料数の残存率は

$$
N/N_0=\frac{1}{\left(\frac{3}{4}\right) + 1}=\frac{4}{7}
$$

である。

密度の残存率は

$$
n/n_0=\frac{\left(\frac{4}{7}\right)}{\left(8\right)}=\frac{1}{14}
$$

である。

初期燃料数のうち反応した割合は

$$
f_{\rm burn}=1 - \left(\frac{4}{7}\right)=\frac{3}{7}
$$

である。

この定温模型の最終割合は

$$
f_{\rm burn,\infty}=\frac{\left(1\right)}{\left(1\right) + 1}=\frac{1}{2}
$$

である。

同じ希釈でも燃料消費は増える。現実の温度変化は反応度を変えるため別に解く。
:::


## 核融合炉の課題

発生したエネルギーを連続して取り出すには、プラズマ内部の収支と周囲の部材の収支を同時に考える。壁への熱流束を単純化して面積 $A$ に一様に分ければ $P/A$ になる。粒子と熱が細い磁束管に集まると接触面積は小さくなり、局所の熱負荷が平均を大きく上回りうる。中性子による発熱と部材の損傷は、表面の熱流とは別にエネルギーと粒子の輸送を追う。

### 発電と燃料の収支

熱交換から電力への効率を $\eta_{\rm th}$、外部加熱器への電力からプラズマに入る加熱への効率を $\eta_{\rm h}$、補機電力を $P_{\rm aux}$ とする。熱出力を融合出力だけで近似した模型で、正味電力は $\eta_{\rm th}P_f-P_{\rm ext}/\eta_{\rm h}-P_{\rm aux}$ である。したがって補機を無視しても $Q>1/(\eta_{\rm th}\eta_{\rm h})$ が必要になる。実際の熱収支では外部加熱の回収とブランケットでの発熱を含める。

DT 系では燃焼一回につき T 一個が消費される。ブランケットで得た T のうち利用できる割合を $f_{\rm ret}$、生成個数の燃焼数に対する比を $b$ とすれば、最も単純な定常補給条件は $b f_{\rm ret}\ge1$ である。回収の時間遅れ、放射壊変、在庫と立ち上げの消費を入れると余裕が必要になる。例えば燃焼数の時間率 $R_T$ が一定で、回収の遅れを無視し、利用可能な在庫 $N_T$ が消滅率 $\lambda_T N_T$ で減るとする。外部補給がないときの在庫式は
$$
\dot N_T=(b f_{\rm ret}-1)R_T-\lambda_T N_T
$$
である。有限の正の在庫を定常に保つには、単に燃焼量と増殖量を等しくするだけでは足りず、右辺の最初の項が壊変を埋める必要がある。燃料系が運転に必要とする最小在庫が $N_{\min}$ なら、$R_T>0$ のもとで $b f_{\rm ret}\ge 1+\lambda_T N_{\min}/R_T$ を得る。

:::example[二つの利得]
$\eta_{\rm th}=0.4$、$\eta_{\rm h}=0.5$ なら補機を省いた模型でさえ $Q>5$ が必要である。プラズマ上の $Q=1$ ではこの収支の正味電力は負である。
:::

:::example[燃料回収]
$f_{\rm ret}=0.9$ なら、在庫変化を無視しても $b\ge10/9$ が必要である。
:::

:::example[非燃料イオンの希釈]
完全電離した $\alpha$ 粒子の密度を $n_\alpha=0.1n_e$ とし、電子密度 $n_e$ を固定する。準中性より燃料イオン密度は $n_f=n_e-2n_\alpha=0.8n_e$ である。他の条件を保った等量 DT の出力は純燃料の $0.8^2=0.64$ 倍になる。実際には灰の排出が放射と閉じ込めにも影響するため、反応率の比だけでは運転点を比較できない。
:::

### 発電所の補機と正味電力の条件

融合電力を $P_f$、プラズマ利得を $Q=P_f/P_{\rm ext}$ とし、熱変換効率は $\eta_{\rm th}$、加熱器の効率は $\eta_{\rm inj}$ と記す。融合に由来する熱だけを回収すると $P_{\rm net}/P_f=\eta_{\rm th}-1/(\eta_{\rm inj}Q)-p_a$、$p_a=P_{\rm aux}/P_f$ である。この模型の境界は $Q_c=1/[\eta_{\rm inj}(\eta_{\rm th}-p_a)]$ で、$\eta_{\rm th}>p_a$ が必要である。外部加熱も炉内で熱として回収するなら熱源と効率を別々に書き、分母にすでに数えた電力を二重に加えない。以下の入力は効率と規格化した運転点の仮定値である。

:::remark[科学利得の分母]
$P_{\rm ext}$ はプラズマへ実際に渡った加熱である。同じ利得を名乗る実験でも入射電力や放電全体のエネルギーで規格化している場合があるため、値を比べるときは分母と時間窓を確認する。 例えば放電の一部だけで利得が高い場合、立ち上げや停止の間も含む平均では補機の消費が残る。発電量を求めるには運転周期全体の電力を積分する。
:::

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $\eta_{\rm th}$ | 融合熱から電力への変換効率 | 1 |
| $\eta_{\rm inj}$ | 購入電力からプラズマ加熱への効率 | 1 |
| $p_a$ | 融合電力に対する補機電力比 | 1 |
| $Q$ | プラズマ利得 | 1 |

プラズマに渡る加熱比を $P_{\rm ext}/P_f$ とする。電源の損失をまだ含まない。

$$
P_{\rm ext}/P_f=\frac{1}{Q}
$$

電源から消費する加熱電力比を $P_{\rm inj}/P_f$ とする。効率で割る。

$$
P_{\rm inj}/P_f=\frac{P_{\rm ext}/P_f}{\eta_{\rm inj}}
$$

正味電力がゼロの利得を $Q_c$ とする。熱効率から補機比を引く。

$$
Q_c=\frac{1}{\eta_{\rm inj} \left(\eta_{\rm th} - p_{a}\right)}
$$

正味発電比を $P_{\rm net}/P_f$ とする。熱から得る電力から二つの消費を引く。

$$
P_{\rm net}/P_f=- P_{\rm inj}/P_f + \eta_{\rm th} - p_{a}
$$

発生した電力のうち構内で使う比を $P_{\rm recirc}/P_{\rm gross}$ とする。一を越えると外部の電力供給が必要になる。

$$
P_{\rm recirc}/P_{\rm gross}=\frac{P_{\rm inj}/P_f + p_{a}}{\eta_{\rm th}}
$$

正味電力境界に対する利得比を $Q/Q_c$ とする。一より大きければ正である。

$$
Q/Q_c=\frac{Q}{Q_{c}}
$$

:::example[科学的な損益分岐]
プラズマ利得が一のときの電力収支を求める。

入力値を次のように置く。

$$
\begin{aligned}
\eta_{\rm th}&=\frac{2}{5} \\\eta_{\rm inj}&=\frac{1}{2} \\p_a&=\frac{1}{20} \\Q&=1
\end{aligned}
$$

これらを順に代入する。

プラズマに渡る加熱比は

$$
P_{\rm ext}/P_f=\frac{1}{\left(1\right)}=1
$$

である。

電源から消費する加熱電力比は

$$
P_{\rm inj}/P_f=\frac{\left(1\right)}{\left(\frac{1}{2}\right)}=2
$$

である。

正味電力がゼロの利得は

$$
Q_c=\frac{1}{\left(\frac{1}{2}\right) \left(- \left(\frac{1}{20}\right) + \left(\frac{2}{5}\right)\right)}=\frac{40}{7}
$$

である。

正味発電比は

$$
P_{\rm net}/P_f=- \left(2\right) - \left(\frac{1}{20}\right) + \left(\frac{2}{5}\right)=- \frac{33}{20}
$$

である。

発生した電力のうち構内で使う比は

$$
P_{\rm recirc}/P_{\rm gross}=\frac{\left(2\right) + \left(\frac{1}{20}\right)}{\left(\frac{2}{5}\right)}=\frac{41}{8}
$$

である。

正味電力境界に対する利得比は

$$
Q/Q_c=\frac{\left(1\right)}{\left(\frac{40}{7}\right)}=\frac{7}{40}
$$

である。

正味電力は負であり、融合電力と購入電力の収支は別である。
:::

:::example[補機なしで境界]
単純模型の正味発電境界を確認する。

入力値を次のように置く。

$$
\begin{aligned}
\eta_{\rm th}&=\frac{2}{5} \\\eta_{\rm inj}&=\frac{1}{2} \\p_a&=0 \\Q&=5
\end{aligned}
$$

これらを順に代入する。

プラズマに渡る加熱比は

$$
P_{\rm ext}/P_f=\frac{1}{\left(5\right)}=\frac{1}{5}
$$

である。

電源から消費する加熱電力比は

$$
P_{\rm inj}/P_f=\frac{\left(\frac{1}{5}\right)}{\left(\frac{1}{2}\right)}=\frac{2}{5}
$$

である。

正味電力がゼロの利得は

$$
Q_c=\frac{1}{\left(\frac{1}{2}\right) \left(- \left(0\right) + \left(\frac{2}{5}\right)\right)}=5
$$

である。

正味発電比は

$$
P_{\rm net}/P_f=- \left(0\right)=0
$$

である。

発生した電力のうち構内で使う比は

$$
P_{\rm recirc}/P_{\rm gross}=\frac{\left(0\right) + \left(\frac{2}{5}\right)}{\left(\frac{2}{5}\right)}=1
$$

である。

正味電力境界に対する利得比は

$$
Q/Q_c=1
$$

である。

加熱用電力だけで発生電力を使う。
:::

:::example[補機を考慮した利得]
利得が高い場合の正味発電比を求める。

入力値を次のように置く。

$$
\begin{aligned}
\eta_{\rm th}&=\frac{2}{5} \\\eta_{\rm inj}&=\frac{1}{2} \\p_a&=\frac{1}{20} \\Q&=10
\end{aligned}
$$

これらを順に代入する。

プラズマに渡る加熱比は

$$
P_{\rm ext}/P_f=\frac{1}{\left(10\right)}=\frac{1}{10}
$$

である。

電源から消費する加熱電力比は

$$
P_{\rm inj}/P_f=\frac{\left(\frac{1}{10}\right)}{\left(\frac{1}{2}\right)}=\frac{1}{5}
$$

である。

正味電力がゼロの利得は

$$
Q_c=\frac{1}{\left(\frac{1}{2}\right) \left(- \left(\frac{1}{20}\right) + \left(\frac{2}{5}\right)\right)}=\frac{40}{7}
$$

である。

正味発電比は

$$
P_{\rm net}/P_f=- \left(\frac{1}{20}\right) - \left(\frac{1}{5}\right) + \left(\frac{2}{5}\right)=\frac{3}{20}
$$

である。

発生した電力のうち構内で使う比は

$$
P_{\rm recirc}/P_{\rm gross}=\frac{\left(\frac{1}{20}\right) + \left(\frac{1}{5}\right)}{\left(\frac{2}{5}\right)}=\frac{5}{8}
$$

である。

正味電力境界に対する利得比は

$$
Q/Q_c=\frac{\left(10\right)}{\left(\frac{40}{7}\right)}=\frac{7}{4}
$$

である。

残る割合は利得だけでなく補機と効率で決まる。
:::

:::example[加熱器の効率低下]
同じプラズマを維持する電力の影響を調べる。

入力値を次のように置く。

$$
\begin{aligned}
\eta_{\rm th}&=\frac{2}{5} \\\eta_{\rm inj}&=\frac{1}{4} \\p_a&=\frac{1}{20} \\Q&=10
\end{aligned}
$$

これらを順に代入する。

プラズマに渡る加熱比は

$$
P_{\rm ext}/P_f=\frac{1}{\left(10\right)}=\frac{1}{10}
$$

である。

電源から消費する加熱電力比は

$$
P_{\rm inj}/P_f=\frac{\left(\frac{1}{10}\right)}{\left(\frac{1}{4}\right)}=\frac{2}{5}
$$

である。

正味電力がゼロの利得は

$$
Q_c=\frac{1}{\left(\frac{1}{4}\right) \left(- \left(\frac{1}{20}\right) + \left(\frac{2}{5}\right)\right)}=\frac{80}{7}
$$

である。

正味発電比は

$$
P_{\rm net}/P_f=- \left(\frac{1}{20}\right)=- \frac{1}{20}
$$

である。

発生した電力のうち構内で使う比は

$$
P_{\rm recirc}/P_{\rm gross}=\frac{\left(\frac{1}{20}\right) + \left(\frac{2}{5}\right)}{\left(\frac{2}{5}\right)}=\frac{9}{8}
$$

である。

正味電力境界に対する利得比は

$$
Q/Q_c=\frac{\left(10\right)}{\left(\frac{80}{7}\right)}=\frac{7}{8}
$$

である。

プラズマ側の $Q$ が同じでも正味電力の符号が変わる。
:::


## 演習問題

::::exercise{level=basic}
等量 DT で全イオン密度を三倍にすると、温度固定で反応率は何倍か。

:::solution
$n_D n_T=n^2/4$ なので九倍である。
:::
::::

::::exercise{level=basic}
T の割合が $1/4$ のとき、等量と同じ全密度での反応率比を求めよ。

:::solution
積は $(1/4)(3/4)$ であり、等量の $1/4$ で割って $3/4$ である。
:::
::::

::::exercise{level=basic}
$W_p=24\,\mathrm{MJ}$、損失 $12\,\mathrm{MW}$ の $\tau_E$ を求めよ。

:::solution
$W_p/P_{\rm loss}=24/12=2\,\mathrm s$ である。
:::
::::

::::exercise{level=basic}
光学的厚さ $n\sigma_{\rm ion}L=\ln 4$ の生存割合を求めよ。

:::solution
指数減衰に代入し $e^{-\ln 4}=1/4$ である。
:::
::::

::::exercise{level=advanced}
温度を八倍にすると Gamow の $E_0$ は何倍か。

:::solution
指数の最小化で $E_0\propto T^{2/3}$ と得たので、$8^{2/3}=4$ 倍である。
:::
::::

::::exercise{level=advanced}
完全堆積模型の定常で $Q=E_f/E_\alpha$ のとき、$n\tau_E$ と点火値の比を求めよ。

:::solution
$E_f/Q=E_\alpha$ なので有効加熱エネルギーは $2E_\alpha$ である。点火値の分母 $E_\alpha$ と比較すると、必要な $n\tau_E$ は半分である。
:::
::::

::::exercise{level=advanced}
等量の固定体積燃焼で燃料を半分消費する時間を求めよ。

:::solution
$n(t)/n_0=1/(1+n_0Kt/2)$ を $1/2$ に等置して $n_0Kt/2=1$、従って $t=2/(n_0K)$ を得る。
:::
::::

::::exercise{level=advanced}
回収効率 $4/5$、増殖比 $6/5$ で定常維持できるか。

:::solution
燃焼一回当たりの回収数は $(4/5)(6/5)=24/25$ で一個を下回る。外から補給しなければ在庫が減る。
:::
::::

## 確認テスト

::::quiz
DT の反応生成物のうち主に磁場閉じ込めの加熱に寄与するものはどれか。

- [x] 荷電した $\alpha$ 粒子
- [ ] 中性子だけ
- [ ] いずれの生成物も必ず同じ割合

:::explanation
$\alpha$ は荷電しており閉じ込められた軌道上で減速できる。中性子は磁力を受けず、エネルギーの主な堆積先はプラズマの外である。
:::
::::

::::quiz{answer="3" tol="1e-6"}
$W_p=9\,\mathrm{MJ}$、総損失 $3\,\mathrm{MW}$ の閉じ込め時間を秒で答えよ。

:::explanation
定義より $W_p/P_{\rm loss}=9/3=3\,\mathrm s$ である。
:::
::::

::::quiz
プラズマの点火を表す模型条件を選べ。

- [x] 外部加熱なしで $\alpha$ 加熱と損失がつりあう。
- [ ] 外部加熱電力と発電電力が等しい。
- [ ] 反応率がゼロである。

:::explanation
熱収支が $\dot W=P_{\rm ext}+P_\alpha-P_{\rm loss}$ であるため、$P_{\rm ext}=0$ での定常条件は $P_\alpha=P_{\rm loss}$ となる。
:::
::::

::::quiz{answer="4" tol="1e-6"}
$\eta_{\rm th}=0.5$、$\eta_{\rm h}=0.5$、補機無視の単純な発電模型における損益ゼロの $Q$ を答えよ。

:::explanation
$\eta_{\rm th}P_f=P_{\rm ext}/\eta_{\rm h}$ だから $Q=1/(\eta_{\rm th}\eta_{\rm h})=4$ である。
:::
::::

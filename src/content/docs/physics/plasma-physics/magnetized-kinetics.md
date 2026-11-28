---
title: 磁化プラズマの運動論
description: 旋回軌道に沿って熱いプラズマの応答を求め、共鳴と低周波近似を導く。
notation:
  - symbol: '\Omega_s'
    meaning: 符号付き旋回角周波数
  - symbol: 'b_s'
    meaning: 熱 Larmor 半径と垂直波数の積の二乗
  - symbol: '\mu'
    meaning: 磁気モーメント
  - symbol: '\Gamma_n'
    meaning: 旋回平均に現れる Bessel 因子
---

## この章の目標

- Vlasov 方程式を非摂動軌道上で積分する。
- 旋回共鳴と有限旋回半径を区別する。
- ドリフト運動論とジャイロ運動論の順序を説明する。

## 磁化プラズマの線形化 Vlasov 方程式

一様な背景 $\vb B_0=B_0\hat{\vb z}$ に対し、磁場に垂直な速度を $v_\perp(\cos\vartheta,\sin\vartheta)$ とする。符号つき角周波数を $\Omega_s=q_s B_0/m_s$、熱速度を $v_{ts}=\sqrt{T_s/m_s}$ とする。$\dot\vartheta=-\Omega_s$ なので、空間 Fourier 化した :eqref[eq-plasma-vlasov] の線形化は
$$
[-i\omega+i\vb k\cdot\vb v-\Omega_s\partial_\vartheta]f_{s1}
=-\frac{q_s}{m_s}(\vb E_1+\vb v\times\vb B_1)\cdot\nabla_v f_{s0}.
$$
平衡分布は旋回位相に依存しないとする。$\vb r_{-\tau},\vb v_{-\tau}$ を現在の位相点から時間 $\tau$ だけ遡った軌道とし、$\Im\omega>0$ で過去からの応答を取る。微分方程式を軌道に沿い積分すると
$$
f_{s1}(\vb v)=-\frac{q_s}{m_s}\int_0^\infty
e^{i\omega\tau+i\vb k\cdot(\vb r_{-\tau}-\vb r)}
[\vb E_1+\vb v_{-\tau}\times\vb B_1]\cdot
\nabla_{v_{-\tau}} f_{s0}\,\dd\tau.
\label{eq-plasma-hot-orbit}
$$
摂動に初期値がある場合には軌道を遡った初期分布も加わる。式 :eqref[eq-plasma-hot-orbit] は外場に対する因果的な応答である。

:::example[磁場に平行な静電波]
$\vb k=k_\parallel\hat{\vb z}$、$\vb E_1\parallel\vb B_0$ なら $\vb k\cdot(\vb r_{-\tau}-\vb r)=-k_\parallel v_\parallel\tau$ である。力の項も $\vartheta$ に依存しない。時間積分は $i/(\omega-k_\parallel v_\parallel)$ となり、縦応答は無磁化の :eqref[eq-plasma-kinetic-epsilon] に戻る。
:::

:::example[位相のずれの上限]
$k_\parallel=0$ で、軌道が波数方向に移動する距離は $2\rho$ を越えない。従って $k_\perp\rho=0.01$ なら位相差の絶対値は常に $0.02$ 以下である。時間を長く遡ってもこの横位相は直線運動のように発散しない。
:::

### 因果応答の高調波の幅と位相

軌道に沿った一高調波の応答には $\int_0^\infty e^{i(\omega-k_\parallel v_\parallel-n\Omega)\tau}\dd\tau$ が現れる。$\omega=\omega_r+i\gamma$、$\gamma>0$ ならこの積分は $1/(\gamma-i\Delta)$、$\Delta=\omega_r-k_\parallel v_\parallel-n\Omega$ である。実部は $\gamma/(\gamma^2+\Delta^2)$ であり、$\Delta=\gamma u$ と置換して滑らかな試験関数と積分すれば $\gamma\downarrow0$ で $\pi\delta(\Delta)$ に収束する。虚部は $\Delta/(\gamma^2+\Delta^2)$ である。原点近くで試験関数からその原点値を差し引くと可積分になり、残った奇関数部分は対称積分で消える。この結果が主値 $1/\Delta$ である。有限の $\gamma$ を物理的な衝突率と同一視しない。ここでは過去の初期条件を選ぶために上半平面から極限を取っている。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $\Delta$ | 軌道から見た離調 | $\mathrm{s^{-1}}$ |
| $\gamma$ | 上半平面の収束率 | $\mathrm{s^{-1}}$ |
| $k_\parallel v_t$ | 熱運動の Doppler 幅 | $\mathrm{s^{-1}}$ |

収束幅単位の離調を $X$ とする。静止粒子の極からの距離である。

$$
X=\frac{\Delta}{\gamma}
$$

軌道積分の実部を $I_R$ とする。離調に対して偶関数である。

$$
I_R=\frac{1}{\gamma \left(X^{2} + 1\right)}
$$

軌道積分の虚部を $I_I$ とする。離調に対して奇関数である。

$$
I_I=I_{R} X
$$

軌道積分の絶対値を $\abs{I}$ とする。二成分を二乗して加える。

$$
\abs{I}=\sqrt{I_{I}^{2} + I_{R}^{2}}
$$

収束幅と Doppler 幅の比を $\gamma/(k_\parallel v_t)$ とする。速度分布を積分する前に収束幅を選ぶ。

$$
\gamma/(k_\parallel v_t)=\frac{\gamma}{k_\parallel v_t}
$$

$[-を $F_{\rm band}$ とする。\Delta|,|\Delta|]$ 内の実部面積の割合 | Lorentz 形を積分する。

$$
F_{\rm band}=\frac{2 \operatorname{atan}{\left(\left|{X}\right| \right)}}{\pi}
$$

:::example[離調ゼロ]
極の実部を有限幅で計算する。

入力値を次のように置く。

$$
\begin{aligned}
\Delta&=0\,\mathrm{s^{-1}} \\\gamma&=\frac{1}{10}\,\mathrm{s^{-1}} \\k_\parallel v_t&=1\,\mathrm{s^{-1}}
\end{aligned}
$$

これらを順に代入する。

収束幅単位の離調は

$$
X=\frac{\left(0\right)}{\left(\frac{1}{10}\right)}=0
$$

である。

軌道積分の実部は

$$
I_R=\frac{1}{\left(\frac{1}{10}\right) \left(\left(0\right)^{2} + 1\right)}=10\,\mathrm s
$$

である。

軌道積分の虚部は

$$
I_I=\left(0\right) \left(10\right)=0\,\mathrm s
$$

である。

軌道積分の絶対値は

$$
\abs{I}=\sqrt{\left(0\right)^{2} + \left(10\right)^{2}}=10\,\mathrm s
$$

である。

収束幅と Doppler 幅の比は

$$
\gamma/(k_\parallel v_t)=\frac{\left(\frac{1}{10}\right)}{\left(1\right)}=\frac{1}{10}
$$

である。

$[-は

$$
F_{\rm band}=\frac{2 \operatorname{atan}{\left(\left|{\left(0\right)}\right| \right)}}{\pi}=0
$$

である。

位相の虚成分が消え、積分は過去の有効時間を表す。
:::

:::example[収束幅に等しい離調]
実部と虚部を比べる。

入力値を次のように置く。

$$
\begin{aligned}
\Delta&=\frac{1}{10}\,\mathrm{s^{-1}} \\\gamma&=\frac{1}{10}\,\mathrm{s^{-1}} \\k_\parallel v_t&=1\,\mathrm{s^{-1}}
\end{aligned}
$$

これらを順に代入する。

収束幅単位の離調は

$$
X=1
$$

である。

軌道積分の実部は

$$
I_R=\frac{1}{\left(\frac{1}{10}\right) \left(\left(1\right)^{2} + 1\right)}=5\,\mathrm s
$$

である。

軌道積分の虚部は

$$
I_I=\left(1\right) \left(5\right)=5\,\mathrm s
$$

である。

軌道積分の絶対値は

$$
\abs{I}=\sqrt{2} \sqrt{\left(5\right)^{2}}\simeq 7.0711\,\mathrm s
$$

である。

収束幅と Doppler 幅の比は

$$
\gamma/(k_\parallel v_t)=\frac{\left(\frac{1}{10}\right)}{\left(1\right)}=\frac{1}{10}
$$

である。

$[-は

$$
F_{\rm band}=\frac{2 \operatorname{atan}{\left(\left|{\left(1\right)}\right| \right)}}{\pi}=\frac{1}{2}
$$

である。

実部は半値、虚部も同じ大きさになる。
:::

:::example[遠く負側の粒子]
因果応答の虚部の向きを求める。

入力値を次のように置く。

$$
\begin{aligned}
\Delta&=-1\,\mathrm{s^{-1}} \\\gamma&=\frac{1}{10}\,\mathrm{s^{-1}} \\k_\parallel v_t&=1\,\mathrm{s^{-1}}
\end{aligned}
$$

これらを順に代入する。

収束幅単位の離調は

$$
X=\frac{\left(-1\right)}{\left(\frac{1}{10}\right)}=-10
$$

である。

軌道積分の実部は

$$
I_R=\frac{1}{\left(\frac{1}{10}\right) \left(\left(-10\right)^{2} + 1\right)}=\frac{10}{101}\,\mathrm s
$$

である。

軌道積分の虚部は

$$
I_I=\left(-10\right) \left(\frac{10}{101}\right)=- \frac{100}{101}\,\mathrm s
$$

である。

軌道積分の絶対値は

$$
\abs{I}=\sqrt{\left(- \frac{100}{101}\right)^{2} + \left(\frac{10}{101}\right)^{2}}\simeq 0.99504\,\mathrm s
$$

である。

収束幅と Doppler 幅の比は

$$
\gamma/(k_\parallel v_t)=\frac{\left(\frac{1}{10}\right)}{\left(1\right)}=\frac{1}{10}
$$

である。

$[-は

$$
F_{\rm band}=\frac{2 \operatorname{atan}{\left(\left|{\left(-10\right)}\right| \right)}}{\pi}\simeq 0.93655
$$

である。

離調を反転すると実部は同じで虚部は反転する。
:::

:::example[Doppler 幅に対して鋭い極]
狭い収束幅がどれほど積分面積を局在化するか調べる。

入力値を次のように置く。

$$
\begin{aligned}
\Delta&=\frac{1}{10}\,\mathrm{s^{-1}} \\\gamma&=\frac{1}{1000}\,\mathrm{s^{-1}} \\k_\parallel v_t&=1\,\mathrm{s^{-1}}
\end{aligned}
$$

これらを順に代入する。

収束幅単位の離調は

$$
X=\frac{\left(\frac{1}{10}\right)}{\left(\frac{1}{1000}\right)}=100
$$

である。

軌道積分の実部は

$$
I_R=\frac{1}{\left(\frac{1}{1000}\right) \left(\left(100\right)^{2} + 1\right)}=\frac{1000}{10001}\,\mathrm s
$$

である。

軌道積分の虚部は

$$
I_I=\left(100\right) \left(\frac{1000}{10001}\right)=\frac{100000}{10001}\,\mathrm s
$$

である。

軌道積分の絶対値は

$$
\abs{I}=\sqrt{\left(\frac{100000}{10001}\right)^{2} + \left(\frac{1000}{10001}\right)^{2}}\simeq 9.9995\,\mathrm s
$$

である。

収束幅と Doppler 幅の比は

$$
\gamma/(k_\parallel v_t)=\frac{\left(\frac{1}{1000}\right)}{\left(1\right)}=\frac{1}{1000}
$$

である。

$[-は

$$
F_{\rm band}=\frac{2 \operatorname{atan}{\left(\left|{\left(100\right)}\right| \right)}}{\pi}\simeq 0.99363
$$

である。

広い積分区間はほとんど全吸収面積を含むが、ピーク高さは収束幅に依存する。
:::


## 熱いプラズマの誘電テンソル

熱い媒質では電流が場だけでなく波数にも依存する。磁場摂動を $\vb B_1=\vb k\times\vb E_1/\omega$ で除く。三重外積から $[\vb E_1+\vb v'\times\vb B_1]_i=A_{ij}(\vb v')E_{1j}$、$A_{ij}=(1-\vb k\cdot\vb v'/\omega)\delta_{ij}+k_i v'_j/\omega$ である。式 :eqref[eq-plasma-hot-orbit] の電流モーメントは
$$
\sigma_{ij}=-\sum_s\frac{q_s^2}{m_s}
\int\dd^3v\,v_i\int_0^\infty\dd\tau\,
e^{i\omega\tau+i\vb k\cdot(\vb r_{-\tau}-\vb r)}
A_{\ell j}(\vb v_{-\tau})\partial_{v'_\ell}f_{s0}(\vb v_{-\tau}).
\label{eq-plasma-hot-conductivity}
$$
繰り返す $\ell$ は和をとる。$\mathsf K=I+i\mathsf\sigma/(\varepsilon_0\omega)$ が熱い誘電テンソルの積分表示である。任意の旋回対称分布を直接入れられる。等方分布では $\nabla_v f_0$ が $\vb v$ と平行なので、磁気力摂動とその内積はゼロである。温度異方性があるとこの項を捨てられない。

### 旋回高調波

$\vb k_\perp=k_\perp\hat{\vb x}$ とすれば $r_x=R_x-v_\perp\sin\vartheta/\Omega_s$ である。軌道上の位相因子は $\exp[-i(k_\perp v_\perp/\Omega_s)\{\sin(\vartheta+\Omega_s\tau)-\sin\vartheta\}]$ である。$e^{ia\sin\vartheta}=\sum_{n\in\Z}J_n(a)e^{in\vartheta}$ を二つの因子に用いると、時間積分の分母は $\omega-k_\parallel v_\parallel-n\Omega_s$ の和に整理できる。$\vb v'$ の成分に含まれる $\sin\vartheta'$ や $\cos\vartheta'$ は指数を一つずらす。従って直交成分を求めるには $J_n,J_{n\pm1}$ の積を積分する。式 :eqref[eq-plasma-hot-conductivity] を成分展開した標準表については、全成分の規約の照合が必要な場合に :cite[stix-waves] を参照する。

:::example[小さい垂直波数]
$a=k_\perp v_\perp/\Omega_s$ が小さいと $J_0(a)=1-a^2/4+\cdots$、$J_{\pm1}(a)=\pm a/2+\cdots$ である。平行電場の密度応答では $n=0$ が優勢だが、横電場は $v_x,v_y$ の因子を通じて $k_\perp=0$ でも $n=\pm1$ を駆動する。
:::

### 平行電磁応答の Maxwell 積分

等方 Maxwell 分布で $\vb k=k_\parallel\hat{\vb z}$、$k_\parallel>0$ とする。R の円偏波に対し $\vb E_R\cdot\vb v_{-\tau}=E_R(v_x+iv_y)e^{i\Omega_s\tau}/\sqrt2$ である。横速度の Maxwell 平均から $T_s/m_s$ がかかり、$\tau$ 積分は $i/(\omega-k_\parallel v_\parallel+\Omega_s)$ を与える。平行の一変数積分は $-Z(\zeta_s^+)/(\sqrt2 k_\parallel v_{ts})$ である。従って
$$
K_R=1+\sum_s\frac{\omega_{ps}^2}{\omega\sqrt2 k_\parallel v_{ts}}\,
Z\left(\frac{\omega+\Omega_s}{\sqrt2 k_\parallel v_{ts}}\right).
$$
L の式には $\omega-\Omega_s$ が入る。$k_\parallel v_{ts}\to0$ で極から離れれば $Z(\zeta)\sim-1/\zeta$ だから :eqref[eq-plasma-rl] に戻る。有限温度でのことであり、$\omega+\Omega_s=0$ を保ったまま冷たい極限をとると展開は使えない。

:::example[流体近似と共鳴幅]
$|\omega+\Omega_s|=10 k_\parallel v_{ts}$ なら $\zeta_s^+=10/\sqrt2$ である。実部の先頭項に対する次項の比は $1/(2\zeta^2)=1/100$、共鳴粒子数の代表的な因子は $e^{-50}$ である。
:::

### 等方 Maxwell 平衡の熱いテンソルの高調波表示

等方平衡 $F_0$ では $\grad_vF_0=-(m/T)\vb v F_0$ なので摂動磁気力は分布の微分と直交する。$\vb k_\perp=k_\perp\hat{\vb x}$ として $a=k_\perp v_\perp/\Omega$ と定義する。旋回成分を
$$
\vb v(\vartheta)e^{-ia\sin\vartheta}=\sum_{n\in\Z}\vb U_n e^{-in\vartheta},
\qquad
\vb U_n=
\left(\frac{n\Omega}{k_\perp}J_n(a),\,i v_\perp J_n'(a),\,v_\parallel J_n(a)\right)
$$
と書ける。実際、$v_x=v_\perp\cos\vartheta$ は $J_{n-1}+J_{n+1}=2nJ_n/a$ を、$v_y=v_\perp\sin\vartheta$ は $J_{n-1}-J_{n+1}=2J_n'$ を与える。軌道積分で現在の速度因子には複素共役が付き、過去の因子には $e^{-in\Omega\tau}$ が付く。旋回平均は異なる高調波の積を消す。従って :eqref[eq-plasma-hot-conductivity] から、各成分を
$$
K_{ij}=\delta_{ij}-\sum_s\frac{q_s^2}{\varepsilon_0\omega T_s}
\sum_{n\in\Z}\int_0^\infty 2\pi v_\perp\dd v_\perp
\int_{\R}\dd v_\parallel
\frac{F_{s0}\overline{U_{ni}}U_{nj}}{\omega-k_\parallel v_\parallel-n\Omega_s}
$$
と書ける。$k_\perp=0$ は $\vb U_n$ を三角関数の Fourier 係数の極限から求める。縦横の結合を含む式であり、$\omega$ の上半平面では因果軌道積分と同じである。共鳴の残差を求めると、対角成分では $|U_{ni}|^2$、非対角成分では偏波間の位相差が重みに入る。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $a$ | 一粒子の垂直波数と符号付き半径の積 | 1 |
| $u$ | 平行速度と横速さの比 $v_\parallel/v_\perp$ | 1 |
| $n$ | 旋回高調波番号 | 1 |

Fourier 位相の高調波を $J_n$ とする。Fourier 積分の係数である。

$$
J_n=J_{n}\left(a\right)
$$

Bessel 関数の微分を $J_n'$ とする。正弦変調された横成分に入る。

$$
J_n'=\frac{J_{n - 1}\left(a\right)}{2} - \frac{J_{n + 1}\left(a\right)}{2}
$$

横の $x$ 成分の対角重みを $\abs{U_{nx}}^2/v_\perp^2$ とする。$a\ne0$ で直接計算する。

$$
\abs{U_{nx}}^2/v_\perp^2=\frac{J_{n}^{2} n^{2}}{a^{2}}
$$

横の $y$ 成分の対角重みを $\abs{U_{ny}}^2/v_\perp^2$ とする。虚数因子は絶対値の二乗で消える。

$$
\abs{U_{ny}}^2/v_\perp^2=(J_n')^2
$$

平行成分の対角重みを $\abs{U_{nz}}^2/v_\perp^2$ とする。平行速度の二乗を掛ける。

$$
\abs{U_{nz}}^2/v_\perp^2=J_{n}^{2} u^{2}
$$

横の二成分間の偏波重みを $\Im(\overline U_{nx}U_{ny})/v_\perp^2$ とする。$i$ の係数を取る。

$$
\Im(\overline U_{nx}U_{ny})/v_\perp^2=\frac{J_{n} J_n' n}{a}
$$

:::example[平行応答の零次調波]
零次調波でどの成分が残るか調べる。

入力値を次のように置く。

$$
\begin{aligned}
a&=\frac{1}{10} \\u&=1 \\n&=0
\end{aligned}
$$

これらを順に代入する。

Fourier 位相の高調波は

$$
J_n=J_{\left(0\right)}\left(\left(\frac{1}{10}\right)\right)\simeq 0.9975
$$

である。

Bessel 関数の微分は

$$
J_n'=\frac{J_{\left(0\right) - 1}\left(\left(\frac{1}{10}\right)\right)}{2} - \frac{J_{\left(0\right) + 1}\left(\left(\frac{1}{10}\right)\right)}{2}\simeq -0.049938
$$

である。

横の $x$ 成分の対角重みは

$$
\abs{U_{nx}}^2/v_\perp^2=\frac{\left(0\right)^{2} \left(J_{0}\left(\frac{1}{10}\right)\right)^{2}}{\left(\frac{1}{10}\right)^{2}}=0
$$

である。

横の $y$ 成分の対角重みは

$$
\abs{U_{ny}}^2/v_\perp^2=\left(- J_{1}\left(\frac{1}{10}\right)\right)^{2}\simeq 0.0024938
$$

である。

平行成分の対角重みは

$$
\abs{U_{nz}}^2/v_\perp^2=\left(1\right)^{2} \left(J_{0}\left(\frac{1}{10}\right)\right)^{2}\simeq 0.99501
$$

である。

横の二成分間の偏波重みは

$$
\Im(\overline U_{nx}U_{ny})/v_\perp^2=\frac{\left(- J_{1}\left(\frac{1}{10}\right)\right) \left(0\right) \left(J_{0}\left(\frac{1}{10}\right)\right)}{\left(\frac{1}{10}\right)}=0
$$

である。

$x$ 成分はゼロだが平行成分は大きい。
:::

:::example[一次の旋回共鳴]
横の二成分の重みを比べる。

入力値を次のように置く。

$$
\begin{aligned}
a&=\frac{1}{10} \\u&=1 \\n&=1
\end{aligned}
$$

これらを順に代入する。

Fourier 位相の高調波は

$$
J_n=J_{\left(1\right)}\left(\left(\frac{1}{10}\right)\right)\simeq 0.049938
$$

である。

Bessel 関数の微分は

$$
J_n'=\frac{J_{\left(1\right) - 1}\left(\left(\frac{1}{10}\right)\right)}{2} - \frac{J_{\left(1\right) + 1}\left(\left(\frac{1}{10}\right)\right)}{2}\simeq 0.49813
$$

である。

横の $x$ 成分の対角重みは

$$
\abs{U_{nx}}^2/v_\perp^2=\frac{\left(1\right)^{2} \left(J_{1}\left(\frac{1}{10}\right)\right)^{2}}{\left(\frac{1}{10}\right)^{2}}\simeq 0.24938
$$

である。

横の $y$ 成分の対角重みは

$$
\abs{U_{ny}}^2/v_\perp^2\simeq 0.24813
$$

である。

平行成分の対角重みは

$$
\abs{U_{nz}}^2/v_\perp^2=\left(1\right)^{2} \left(J_{1}\left(\frac{1}{10}\right)\right)^{2}\simeq 0.0024938
$$

である。

横の二成分間の偏波重みは

$$
\Im(\overline U_{nx}U_{ny})/v_\perp^2=\frac{\left(- \frac{J_{2}\left(\frac{1}{10}\right)}{2} + \frac{J_{0}\left(\frac{1}{10}\right)}{2}\right) \left(1\right) \left(J_{1}\left(\frac{1}{10}\right)\right)}{\left(\frac{1}{10}\right)}\simeq 0.24875
$$

である。

両者は約四分の一で、平行成分より大きい。
:::

:::example[二次の旋回共鳴]
垂直位相が一になると二次調波が現れるか計算する。

入力値を次のように置く。

$$
\begin{aligned}
a&=1 \\u&=\frac{1}{2} \\n&=2
\end{aligned}
$$

これらを順に代入する。

Fourier 位相の高調波は

$$
J_n=J_{\left(2\right)}\left(\left(1\right)\right)\simeq 0.1149
$$

である。

Bessel 関数の微分は

$$
J_n'=\frac{J_{\left(2\right) - 1}\left(\left(1\right)\right)}{2} - \frac{J_{\left(2\right) + 1}\left(\left(1\right)\right)}{2}\simeq 0.21024
$$

である。

横の $x$ 成分の対角重みは

$$
\abs{U_{nx}}^2/v_\perp^2=\frac{\left(2\right)^{2} \left(J_{2}\left(1\right)\right)^{2}}{\left(1\right)^{2}}\simeq 0.052811
$$

である。

横の $y$ 成分の対角重みは

$$
\abs{U_{ny}}^2/v_\perp^2\simeq 0.044202
$$

である。

平行成分の対角重みは

$$
\abs{U_{nz}}^2/v_\perp^2=\left(J_{2}\left(1\right)\right)^{2} \left(\frac{1}{2}\right)^{2}\simeq 0.0033007
$$

である。

横の二成分間の偏波重みは

$$
\Im(\overline U_{nx}U_{ny})/v_\perp^2=\frac{\left(- \frac{J_{3}\left(1\right)}{2} + \frac{J_{1}\left(1\right)}{2}\right) \left(2\right) \left(J_{2}\left(1\right)\right)}{\left(1\right)}\simeq 0.048315
$$

である。

高調波の重みは有限になるが、波の分母も共鳴の評価に要る。
:::

:::example[回転の向きを固定し高調波を反転]
対角と非対角の重みの対称性を比較する。

入力値を次のように置く。

$$
\begin{aligned}
a&=1 \\u&=\frac{1}{2} \\n&=-2
\end{aligned}
$$

これらを順に代入する。

Fourier 位相の高調波は

$$
J_n=J_{\left(-2\right)}\left(\left(1\right)\right)\simeq 0.1149
$$

である。

Bessel 関数の微分は

$$
J_n'=\frac{J_{\left(-2\right) - 1}\left(\left(1\right)\right)}{2} - \frac{J_{\left(-2\right) + 1}\left(\left(1\right)\right)}{2}\simeq 0.21024
$$

である。

横の $x$ 成分の対角重みは

$$
\abs{U_{nx}}^2/v_\perp^2=\frac{\left(-2\right)^{2} \left(J_{2}\left(1\right)\right)^{2}}{\left(1\right)^{2}}\simeq 0.052811
$$

である。

横の $y$ 成分の対角重みは

$$
\abs{U_{ny}}^2/v_\perp^2\simeq 0.044202
$$

である。

平行成分の対角重みは

$$
\abs{U_{nz}}^2/v_\perp^2=\left(J_{2}\left(1\right)\right)^{2} \left(\frac{1}{2}\right)^{2}\simeq 0.0033007
$$

である。

横の二成分間の偏波重みは

$$
\Im(\overline U_{nx}U_{ny})/v_\perp^2=\frac{\left(- \frac{J_{3}\left(1\right)}{2} + \frac{J_{1}\left(1\right)}{2}\right) \left(-2\right) \left(J_{2}\left(1\right)\right)}{\left(1\right)}\simeq -0.048315
$$

である。

対角は同じになり、偏波の虚結合の符号は反転する。
:::


### 垂直速度まで積分した Maxwell 応答

上の高調波表示では等方性を使ったため、垂直速度積分を解析的にできる。一種について $w=v_\perp/v_t$、$z=v_\parallel/v_t$ とし、$\langle Q\rangle_\perp=\int_0^\infty Q(w)w e^{-w^2/2}\dd w$ と書く。$b=k_\perp^2v_t^2/\Omega^2>0$、$s_\Omega=\sgn\Omega$ とすると $a=s_\Omega\sqrt b\,w$ である。まず垂直速度の平均を求める。二時刻の旋回角の差を $\lambda$ とすると、二つの Jacobi–Anger 展開を掛けて初期角を平均した位相因子は $\sum_n J_n(a)^2 e^{-in\lambda}$ である。同じ位相差は、垂直速度の独立な二つの Gaussian 成分に係数 $(k_\perp/\Omega)(\sin\lambda,\cos\lambda-1)$ を掛けた和である。従って速度平均した位相因子は $\exp[-b(1-\cos\lambda)]$ になる。$e^{b\cos\lambda}=\sum_n I_n(b)e^{-in\lambda}$ と係数を比較すれば
$$
\langle J_n(a)^2\rangle_\perp=\Gamma_n(b),\qquad\Gamma_n(b)=e^{-b}I_n(b)
\label{eq-plasma-thermal-bessel}
$$
である。$b$ で微分し、また $w$ で部分積分すれば
$$
\langle aJ_n J_n'\rangle_\perp=b\Gamma_n',
\qquad
\langle w^2J_n^2\rangle_\perp=2\Gamma_n+2b\Gamma_n'
$$
となる。第一式の微分は $J_n(a)$ の引数に、第二式の部分積分は重み $w e^{-w^2/2}$ に作用する。

$J_n$ の方程式 $J_n''+J_n'/a+(1-n^2/a^2)J_n=0$ を使って $\Gamma_n$ の二階微分を整理する。すると
$$
2b\Gamma_n''=\langle w^2(J_n')^2\rangle_\perp
-\langle w^2J_n^2\rangle_\perp+(n^2/b)\Gamma_n-2\Gamma_n'
$$
を得る。一方、$\Gamma_n=e^{-b}I_n(b)$ と変形 Bessel 方程式から $b\Gamma_n''+(1+2b)\Gamma_n'+(1-n^2/b)\Gamma_n=0$ である。これで $yy$ 成分の平均を求めると、$\langle w^2(J_n')^2\rangle_\perp=(n^2/b)\Gamma_n-2b\Gamma_n'$ が得られる。

$\overline{U_{ni}}U_{nj}/v_t^2$ の垂直平均を $\mathsf P_n(z)$ とすれば、これらの結果と $xz,yz$ 成分の積から
$$
\mathsf P_n(z)=
\begin{pmatrix}
n^2\Gamma_n/b & i n\Gamma_n' & s_\Omega n z\Gamma_n/\sqrt b\\
-i n\Gamma_n' & n^2\Gamma_n/b-2b\Gamma_n' & -i s_\Omega \sqrt b\,z\Gamma_n'\\
s_\Omega n z\Gamma_n/\sqrt b & i s_\Omega \sqrt b\,z\Gamma_n' & z^2\Gamma_n
\end{pmatrix}.
$$
各成分を平行 Maxwell 分布で積分し、$\mathsf K=I-\sum_s (\omega_{ps}^2/\omega)\sum_n\int \mathsf P_n(z)e^{-z^2/2}\dd z/[\sqrt{2\pi}(\omega-k_\parallel v_{ts}z-n\Omega_s)]$ とする。$k_\parallel>0$ なら、$z^0,z^1,z^2$ の三つの積分は、共通因子 $-1/(\sqrt2 k_\parallel v_{ts})$ を除いて
$$
Z(\zeta_n),\qquad
\sqrt2[1+\zeta_n Z(\zeta_n)],\qquad
2\zeta_n[1+\zeta_n Z(\zeta_n)]
$$
である。$\zeta_n=(\omega-n\Omega_s)/(\sqrt2 k_\parallel v_{ts})$ とする。これらは $x/(x-\zeta)=1+\zeta/(x-\zeta)$ と $x^2/(x-\zeta)=x+\zeta+\zeta^2/(x-\zeta)$ から得られる。したがって数値計算には $I_n$ と $Z$ の評価だけが残る。

:::example[垂直波数ゼロでの偏波]
$b\to0$ では $n=\pm1$ の $xx,yy$ 成分が $1/2$、$xy$ 成分が $\pm i/2$ に近づき、$n=0$ の $zz$ 成分が $z^2$ に近づく。$x+iy$ 方向の円偏波には $n=-1$ の行列だけが作用し、前に求めた $K_R$ の分母が得られる。この計算では $\Omega$ の符号を保持する。
:::

## サイクロトロン減衰

粒子から見た波の周波数 $\omega-k_\parallel v_\parallel$ が旋回高調波 $n\Omega_s$ と一致すると、毎回の作用が同位相で足される。この過程による吸収を :term[サイクロトロン減衰]{en="cyclotron damping" yomi="さいくろとろんげんすい"} という。$k_\parallel>0$ なら共鳴速度は $v_{\parallel r}=(\omega_r-n\Omega_s)/k_\parallel$ である。実周波数へ上から近づけると $1/(\omega-k_\parallel v_\parallel-n\Omega_s+i0)$ の虚部は $-\pi\delta(\omega-k_\parallel v_\parallel-n\Omega_s)$ である。速度積分で $\delta(v_\parallel-v_{\parallel r})/\abs{k_\parallel}$ が付く。

波と粒子の共鳴では、力に沿う分布勾配は $k_\parallel\partial_{v_\parallel}f_0+(n\Omega_s/v_\perp)\partial_{v_\perp}f_0$ に比例する。これは波の位相に沿うエネルギー・平行運動量交換からも導ける。$\delta p_\parallel=(k_\parallel/\omega)\delta W$ とし、共鳴上で $\delta W=m(v_\parallel\delta v_\parallel+v_\perp\delta v_\perp)$ を使うと速度変化の方向は上の微分と一致する。等方 Maxwell 分布なら方向微分は $-(m/T)(k_\parallel v_\parallel+n\Omega_s)f_0=-(m/T)\omega_r f_0$ である。正エネルギー波を減衰させる。異方性分布では符号が逆になることがある。

:::example[Doppler 補正]
$n\Omega_s=8\,\mathrm{s^{-1}}$、$\omega_r=10\,\mathrm{s^{-1}}$、$k_\parallel=2\,\mathrm{m^{-1}}$ なら共鳴速度は $1\,\mathrm{m\,s^{-1}}$ である。
:::

:::example[温度による幅]
共鳴速度が平均流から $3v_{ts}$ 離れていると、Maxwell 分布のそこでの値は中心の $e^{-9/2}$ 倍である。$k_\parallel$ を二倍にして周波数差を保つと共鳴速度は半分になる。ただし吸収には分布の値だけでなく、方向微分と $1/\abs{k_\parallel}$ も入る。
:::

### 吸収の偏波依存性

等方 Maxwell 平衡では、導電率を高調波で表したとき各共鳴項の実電力は非負である。複素振幅の実部を物理場とし、$\omega_r>0$、$k_\parallel\ne0$ とする。$P_s=\Re(\vb E^\dagger\mathsf\sigma_s\vb E)/2$ に軌道積分を入れ、前節の高調波平均と $\Re[i/(x+i0)]=\pi\delta(x)$ を使うと
$$
P_s=\frac{\pi q_s^2}{2T_s}\sum_n
\int_0^\infty 2\pi v_\perp\dd v_\perp\int_\R\dd v_\parallel\,
F_{s0}|\vb U_n\cdot\vb E|^2\,\delta(\omega_r-k_\parallel v_\parallel-n\Omega_s).
$$
この式は、共鳴粒子がいても $\vb U_n\cdot\vb E=0$ ならその高調波による吸収がないことを示す。正の波エネルギー密度を $W$ とし弱減衰で評価すれば、$2\gamma W=-\sum_s P_s$ である。係数二は波のエネルギーが振幅の二乗に比例するためである。

### 異方性 Maxwell 分布の共鳴方向微分

$F_0\propto\exp[-m(v_\parallel-U)^2/(2T_\parallel)-mv_\perp^2/(2T_\perp)]$ とする。共鳴方向 $\mathcal G=k_\parallel\partial_{v_\parallel}+(n\Omega/v_\perp)\partial_{v_\perp}$ に作用させると $\mathcal G F_0=-(m/T_\parallel)F_0[(\omega-n\Omega)-k_\parallel U+(T_\parallel/T_\perp)n\Omega]$ である。角括弧の符号は単なる速度勾配の符号と一致しない場合がある。背景に静止した等方な種では正周波数の波に対し角括弧は $\omega$ となる。正エネルギー波でこの種が放出側に回る条件は角括弧が負になることである。全成長率の符号を決めるには他種、他高調波、波のエネルギーを含める。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $N$ | $n\Omega/\omega$ の符号付き比 | 1 |
| $U$ | 平均流の Doppler 比 $k_\parallel U/\omega$ | 1 |
| $A$ | 異方性 $T_\perp/T_\parallel$ | 1 |
| $W$ | $\omega/(k_\parallel\sqrt{T_\parallel/m})>0$ | 1 |

熱速さで規格化した共鳴速度のずれを $x_r$ とする。平均流の速度を差し引く。

$$
x_r=W \left(- N - U + 1\right)
$$

同じ横速さでの共鳴粒子の比を $F_r/F_c$ とする。Maxwell の平行速度因子である。

$$
F_r/F_c=e^{- \frac{x_{r}^{2}}{2}}
$$

平行微分の符号付き寄与を $G_\parallel$ とする。分布を微分した項から共通の負係数を除く。

$$
G_\parallel=- N - U + 1
$$

横微分の符号付き寄与を $G_\perp$ とする。異方性によって重みが変わる。

$$
G_\perp=\frac{N}{A}
$$

角括弧の周波数単位での値を $G$ とする。正なら正エネルギー波への吸収に寄与する。

$$
G=G_{\parallel} + G_{\perp}
$$

Maxwell の粒子数を含めた重みを $D_r$ とする。偏波因子と波エネルギーは含まない。

$$
D_r=\frac{F_r}{F_c} G
$$

:::example[等方・静止の共鳴]
平行微分が負の寄与でも全体を評価する。

入力値を次のように置く。

$$
\begin{aligned}
N&=2 \\U&=0 \\A&=1 \\W&=1
\end{aligned}
$$

これらを順に代入する。

熱速さで規格化した共鳴速度のずれは

$$
x_r=\left(1\right) \left(- \left(0\right) - \left(2\right) + 1\right)=-1
$$

である。

同じ横速さでの共鳴粒子の比は

$$
F_r/F_c=e^{- \frac{\left(-1\right)^{2}}{2}}\simeq 0.60653
$$

である。

平行微分の符号付き寄与は

$$
G_\parallel=- \left(0\right) - \left(2\right) + 1=-1
$$

である。

横微分の符号付き寄与は

$$
G_\perp=\frac{\left(2\right)}{\left(1\right)}=2
$$

である。

角括弧の周波数単位での値は

$$
G=\left(-1\right) + \left(2\right)=1
$$

である。

Maxwell の粒子数を含めた重みは

$$
D_r=\left(1\right) \left(e^{- \frac{1}{2}}\right)\simeq 0.60653
$$

である。

横運動からの項が加わるため吸収側になる。
:::

:::example[高い垂直温度]
横の微分が弱くなったときの符号を求める。

入力値を次のように置く。

$$
\begin{aligned}
N&=2 \\U&=0 \\A&=4 \\W&=1
\end{aligned}
$$

これらを順に代入する。

熱速さで規格化した共鳴速度のずれは

$$
x_r=\left(1\right) \left(- \left(0\right) - \left(2\right) + 1\right)=-1
$$

である。

同じ横速さでの共鳴粒子の比は

$$
F_r/F_c=e^{- \frac{\left(-1\right)^{2}}{2}}\simeq 0.60653
$$

である。

平行微分の符号付き寄与は

$$
G_\parallel=- \left(0\right) - \left(2\right) + 1=-1
$$

である。

横微分の符号付き寄与は

$$
G_\perp=\frac{\left(2\right)}{\left(4\right)}=\frac{1}{2}
$$

である。

角括弧の周波数単位での値は

$$
G=\left(-1\right) + \left(\frac{1}{2}\right)=- \frac{1}{2}
$$

である。

Maxwell の粒子数を含めた重みは

$$
D_r=\left(- \frac{1}{2}\right) \left(e^{- \frac{1}{2}}\right)\simeq -0.30327
$$

である。

符号が反転し、この調波に限れば波へエネルギーを与える。
:::

:::example[境界の垂直温度]
二つの微分が打ち消す条件を確かめる。

入力値を次のように置く。

$$
\begin{aligned}
N&=2 \\U&=0 \\A&=2 \\W&=1
\end{aligned}
$$

これらを順に代入する。

熱速さで規格化した共鳴速度のずれは

$$
x_r=\left(1\right) \left(- \left(0\right) - \left(2\right) + 1\right)=-1
$$

である。

同じ横速さでの共鳴粒子の比は

$$
F_r/F_c=e^{- \frac{\left(-1\right)^{2}}{2}}\simeq 0.60653
$$

である。

平行微分の符号付き寄与は

$$
G_\parallel=- \left(0\right) - \left(2\right) + 1=-1
$$

である。

横微分の符号付き寄与は

$$
G_\perp=1
$$

である。

角括弧の周波数単位での値は

$$
G=\left(-1\right) + \left(1\right)=0
$$

である。

Maxwell の粒子数を含めた重みは

$$
D_r=\left(0\right) \left(e^{- \frac{1}{2}}\right)=0
$$

である。

共鳴粒子が存在してもこの調波の符号付き寄与はゼロである。
:::

:::example[大きい平行流]
等方でも流れのために符号が変わるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
N&=\frac{1}{2} \\U&=2 \\A&=1 \\W&=2
\end{aligned}
$$

これらを順に代入する。

熱速さで規格化した共鳴速度のずれは

$$
x_r=\left(2\right) \left(- \left(2\right) - \left(\frac{1}{2}\right) + 1\right)=-3
$$

である。

同じ横速さでの共鳴粒子の比は

$$
F_r/F_c=e^{- \frac{\left(-3\right)^{2}}{2}}\simeq 0.011109
$$

である。

平行微分の符号付き寄与は

$$
G_\parallel=- \left(2\right) - \left(\frac{1}{2}\right) + 1=- \frac{3}{2}
$$

である。

横微分の符号付き寄与は

$$
G_\perp=\frac{\left(\frac{1}{2}\right)}{\left(1\right)}=\frac{1}{2}
$$

である。

角括弧の周波数単位での値は

$$
G=\left(- \frac{3}{2}\right) + \left(\frac{1}{2}\right)=-1
$$

である。

Maxwell の粒子数を含めた重みは

$$
D_r=\left(-1\right) \left(e^{- \frac{9}{2}}\right)\simeq -0.011109
$$

である。

波を追い越す平均流は供給源になり得るが、速度分布の中心から共鳴が離れると寄与は弱まる。
:::


## Bernstein 波

垂直静電波 $\vb k=k_\perp\hat{\vb x}$ を考える。非相対論的な Maxwell 分布なら、共鳴周波数の Doppler 幅がこの方向で消えるため離散的な高調波の間に固有波が現れる。これを :term[Bernstein 波]{en="Bernstein wave" yomi="ばーんすたいんは"} という。

### 縦誘電関数の導出

$\vb E=-i\vb k\phi$、$\nabla_v f_{s0}=-(m_s/T_s)\vb v f_{s0}$ とする。軌道位相 $G(\tau)=e^{i\vb k\cdot(\vb r_{-\tau}-\vb r)}$ は $G'=-i\vb k\cdot\vb v_{-\tau}G$ を満たす。式 :eqref[eq-plasma-hot-orbit] を部分積分すると $f_{s1}=-(q_s\phi/T_s)f_{s0}[1+i\omega\int_0^\infty e^{i\omega\tau}G\dd\tau]$ となる。旋回平均した $G$ は $\sum_n J_n(a)^2 e^{-in\Omega_s\tau}$ である。さらに :eqref[eq-plasma-thermal-bessel] で速度を平均すると
$$
\langle G(\tau)\rangle=\sum_{n\in\Z}\Gamma_n(b_s)e^{-in\Omega_s\tau},\qquad
b_s=\frac{k_\perp^2 T_s}{m_s\Omega_s^2}
$$
となる。$\sum_n\Gamma_n=1$ と $\Gamma_{-n}=\Gamma_n$ を用いて Poisson 方程式へ入れると
$$
\varepsilon_\perp=1-\sum_s\frac{2\omega_{ps}^2}{b_s}
\sum_{n=1}^\infty \frac{n^2\Gamma_n(b_s)}{\omega^2-n^2\Omega_s^2}=0.
\label{eq-plasma-bernstein}
$$
一つの可動種で $b_s>0$ なら、隣り合う二つの極の間で $\partial_{\omega^2}\varepsilon_\perp>0$ であり、両端で $-\infty,+\infty$ となる。従って各間隔に一つの根がある。$b_s\to0$ では $\Gamma_1\simeq b_s/2$、$\Gamma_n=\order{b_s^n}$ なので極から離れて冷たい :eqref[eq-plasma-hybrid] に一致する。

:::example[二次高調波の重み]
$\Gamma_2\simeq b^2/8$ なので $2(2^2)\Gamma_2/b\simeq b$ である。小さな $b$ では極から離れた応答は弱いが、共鳴分母が小さい場所では捨てられない。
:::

:::example[垂直波の非共鳴区間]
イオンを固定し、$b_e>0$ とする。$\omega^2<\Omega_e^2$ なら式 :eqref[eq-plasma-bernstein] の全ての分母が負であり $\varepsilon_\perp>1$ となる。従って最初の旋回高調波より下にはこの模型の垂直静電固有根がない。低域混成波を求めるときにイオンを固定できない理由でもある。
:::

### 任意の角度での静電応答

同じ等方 Maxwell 平衡で $k_\parallel>0$ を残す。上で部分積分した $G$ の旋回・速度平均を行うと、$v_\parallel$ ごとに $\sum_n\Gamma_n(b)e^{-i(k_\parallel v_\parallel+n\Omega_s)\tau}$ が得られる。時間積分では $\int_0^\infty e^{i(\omega-k_\parallel v_\parallel-n\Omega_s)\tau}\dd\tau=i/(\omega-k_\parallel v_\parallel-n\Omega_s)$ である。平行速度をさらに積分して Poisson 方程式に入れると
$$
\varepsilon_L(\omega,\vb k)
=1+\sum_s\frac{1}{k^2\lambda_{Ds}^2}
\left[1+\frac{\omega}{\sqrt2 k_\parallel v_{ts}}
\sum_{n\in\Z}\Gamma_n(b_s)
Z\!\left(\frac{\omega-n\Omega_s}{\sqrt2 k_\parallel v_{ts}}\right)\right].
$$
ここで $\varepsilon_L=\hat{\vb k}\cdot\mathsf K\hat{\vb k}$ であり、横電場を捨てた静電近似で零点を探す。$b_s=0$ なら $n=0$ のみが残り、無磁化縦波と同じ応答である。$k_\parallel\to0$ を極から離れてとると $Z(\zeta)\sim-1/\zeta$ なので括弧は $1-\sum_n\Gamma_n \omega/(\omega-n\Omega_s)$ となる。$\sum_n\Gamma_n=1$ を用い、正負の $n$ を合わせれば式 :eqref[eq-plasma-bernstein] に戻る。$k_\parallel$ が小さくても $\abs{\omega-n\Omega_s}$ が $k_\parallel v_{ts}$ まで小さくなると、この極限と速度積分を交換できない。

:::example[位相速度の比較]
$n=0$ の共鳴速度は $\omega/k_\parallel$、$n=1$ では $(\omega-\Omega_s)/k_\parallel$ である。例えば $\Omega_s/(k_\parallel v_{ts})=20$、$\omega-\Omega_s=k_\parallel v_{ts}$ なら前者は $21v_{ts}$、後者は $v_{ts}$ である。有限の $b_s$ では同じ波の虚部を二つの共鳴が大きく異なる重みで担う。
:::

### 混成共鳴付近で電磁波と静電波が接近する条件

電子だけが運動する垂直伝播で、$\omega$ は旋回高調波の極から離れているとする。式 :eqref[eq-plasma-bernstein] の小 $b$ 展開に $\Gamma_1=b/2-b^2/2+\order{b^3}$、$\Gamma_2=b^2/8+\order{b^3}$ を入れると
$$
K_{xx}=S+b A+\order{b^2},\qquad
A=\omega_p^2\left(\frac{1}{\omega^2-\Omega^2}-\frac{1}{\omega^2-4\Omega^2}\right)
$$
である。$S$ は冷たい縦成分である。$|\Omega|<\omega<2|\Omega|$ では $A>0$ となる。静電近似は $K_{xx}=0$ であるが、高域混成の近くでは小さな $S$ を有限の電磁効果とも比較する必要がある。

垂直 X 波では $x,y$ の Maxwell 行列の行列式が $K_{xx}(K_{yy}-N_\perp^2)-K_{xy}K_{yx}$、$N_\perp=ck_\perp/\omega$ である。$N_\perp^2\gg|K_{yy}|$、$b\ll1$ で、非対角成分に冷たい値 $K_{xy}=-iD$、$K_{yx}=iD$ を用いると
$$
S+b A+\frac{D^2}{N_\perp^2}\simeq0 .
$$
$b=k_\perp^2v_t^2/\Omega^2$ を使えば $A b^2+Sb+C\simeq0$、$C=D^2(\omega v_t/(c\Omega))^2$ となる。$A,C>0$、$S<-2\sqrt{AC}$ では正の二根があり、$|S|\gg\sqrt{AC}$ なら一方は $C/|S|$、他方は $|S|/A$ に近い。前者の根では電磁項、後者では有限旋回半径の項が $S$ を打ち消す。例えば $S=-3\sqrt{AC}$ なら $b/\sqrt{C/A}=(3\pm\sqrt5)/2$ であり、$S=-2\sqrt{AC}$ で二根が重なる。この計算は二根がともに $b\ll1$、$N_\perp^2\gg|K_{yy}|$ の領域にある場合に使う。密度がゆっくり変化する媒質で根が重なる近くでは、分枝ごとの伝播近似が破れる。変換率を求めるにはその領域の波動方程式を境界条件とともに解く。 :cite[bellan-plasma]

### 有限の高調波和で得た根の誤差

一つの可動種について $X=\omega^2/\Omega^2$ と書く。第 $N$ 高調波で和を切り、$0<X<(N+1)^2$ の非共鳴区間を調べる。捨てた高調波による補正 $R_N(X)=\sum_{n>N} A_n/(n^2-X)$ は正である。したがって切断後の根 $X_N$ と同じ極の間に真の根 $X_\infty$ があるとき、誘電関数の単調性から $X_\infty<X_N$ である。両根を含む閉区間で切断した誘電関数の導関数の下限を $d>0$ とすれば
$$
0<X_N-X_\infty\le \frac{R_N(X_N)}{d}.
$$
実際、切断した関数は二根の間で $R_N(X_\infty)$ だけ増え、$R_N$ 自体も増加する。次の極に近いと $R_N(X_N)$ が大きくなるので、高調波の重み $\Gamma_n$ だけで切断の妥当性を決められない。

### 二極を保持した Bernstein 波の零点

一つの可動種で垂直静電の式 :eqref[eq-plasma-bernstein] を二調波まで切る。$X=\omega^2/\Omega^2$ とし $A_n=2(\omega_p^2/\Omega^2)n^2\Gamma_n(b)/b$ とすると $\varepsilon=1-A_1/(X-1)-A_2/(X-4)$ である。分母を払うと $X^2-(5+A_1+A_2)X+4+4A_1+A_2=0$ となる。$A_1,A_2>0$ なら一根は $1<X<4$、他根は $X>4$ にある。実際、これらの区間では $\varepsilon'>0$ で区間の両端の符号が逆である。実際の上側の枝には三次以降の極がある。ここで求める上側根は二極近似であり、次の極に接近すると誤差を再評価する。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $b$ | 波数と熱旋回半径の積の二乗 | 1 |
| $P$ | $\omega_p^2/\Omega^2$ | 1 |

一次の留数を $A_1$ とする。垂直熱平均を含む。

$$
A_1=\frac{2 P e^{- b} I_{1}\left(b\right)}{b}
$$

二次の留数を $A_2$ とする。$n^2=4$ が掛かる。

$$
A_2=\frac{8 P e^{- b} I_{2}\left(b\right)}{b}
$$

分母を払った二根の和を $S$ とする。二つの極の位置の和が含まれる。

$$
S=A_{1} + A_{2} + 5
$$

分母を払った二根の積を $P_X$ とする。定数項をまとめる。

$$
P_X=4 A_{1} + A_{2} + 4
$$

第二極より上の形式根を $X_+$ とする。三極目は含まない。

$$
X_+=\frac{S}{2} + \frac{\sqrt{- 4 P_{X} + S^{2}}}{2}
$$

第一極と第二極の間の根を $X_-$ とする。小さい根は積から求める。

$$
X_-=\frac{P_{X}}{X_{+}}
$$

:::example[小さい半径補正]
第二調波を弱くした場合の二根を求める。

入力値を次のように置く。

$$
\begin{aligned}
b&=\frac{1}{100} \\P&=1
\end{aligned}
$$

これらを順に代入する。

一次の留数は

$$
A_1=\frac{2 \left(1\right) e^{- \left(\frac{1}{100}\right)} I_{1}\left(\left(\frac{1}{100}\right)\right)}{\left(\frac{1}{100}\right)}\simeq 0.99006
$$

である。

二次の留数は

$$
A_2=\frac{8 \left(1\right) e^{- \left(\frac{1}{100}\right)} I_{2}\left(\left(\frac{1}{100}\right)\right)}{\left(\frac{1}{100}\right)}\simeq 0.0099006
$$

である。

分母を払った二根の和は

$$
S=\left(\frac{200 I_{1}\left(\frac{1}{100}\right)}{e^{\frac{1}{100}}}\right) + \left(\frac{800 I_{2}\left(\frac{1}{100}\right)}{e^{\frac{1}{100}}}\right) + 5\simeq 6
$$

である。

分母を払った二根の積は

$$
P_X=4 \left(\frac{200 I_{1}\left(\frac{1}{100}\right)}{e^{\frac{1}{100}}}\right) + \left(\frac{800 I_{2}\left(\frac{1}{100}\right)}{e^{\frac{1}{100}}}\right) + 4\simeq 7.9701
$$

である。

第二極より上の形式根は

$$
X_+=\frac{\left(\frac{800 I_{2}\left(\frac{1}{100}\right)}{e^{\frac{1}{100}}} + \frac{200 I_{1}\left(\frac{1}{100}\right)}{e^{\frac{1}{100}}} + 5\right)}{2} + \frac{\sqrt{\left(\left(\frac{800 I_{2}\left(\frac{1}{100}\right)}{e^{\frac{1}{100}}} + \frac{200 I_{1}\left(\frac{1}{100}\right)}{e^{\frac{1}{100}}} + 5\right)\right)^{2} - 4 \left(\frac{800 I_{2}\left(\frac{1}{100}\right)}{e^{\frac{1}{100}}} + \frac{800 I_{1}\left(\frac{1}{100}\right)}{e^{\frac{1}{100}}} + 4\right)}}{2}\simeq 4.0147
$$

である。

第一極と第二極の間の根は

$$
X_-=\frac{\left(\frac{800 I_{2}\left(\frac{1}{100}\right)}{e^{\frac{1}{100}}} + \frac{800 I_{1}\left(\frac{1}{100}\right)}{e^{\frac{1}{100}}} + 4\right)}{\left(\frac{400 I_{2}\left(\frac{1}{100}\right)}{e^{\frac{1}{100}}} + \frac{100 I_{1}\left(\frac{1}{100}\right)}{e^{\frac{1}{100}}} + \frac{\sqrt{-16 - \frac{3200 I_{1}\left(\frac{1}{100}\right)}{e^{\frac{1}{100}}} - \frac{3200 I_{2}\left(\frac{1}{100}\right)}{e^{\frac{1}{100}}} + \left(\frac{800 I_{2}\left(\frac{1}{100}\right)}{e^{\frac{1}{100}}} + \frac{200 I_{1}\left(\frac{1}{100}\right)}{e^{\frac{1}{100}}} + 5\right)^{2}}}{2} + \frac{5}{2}\right)}\simeq 1.9852
$$

である。

上の根は第二極に寄る。留数が正の有限値である間だけその枝が存在する。
:::

:::example[垂直波数が増加]
一次と二次の重みの比の変化を調べる。

入力値を次のように置く。

$$
\begin{aligned}
b&=\frac{1}{2} \\P&=1
\end{aligned}
$$

これらを順に代入する。

一次の留数は

$$
A_1=\frac{2 \left(1\right) e^{- \left(\frac{1}{2}\right)} I_{1}\left(\left(\frac{1}{2}\right)\right)}{\left(\frac{1}{2}\right)}\simeq 0.62568
$$

である。

二次の留数は

$$
A_2=\frac{8 \left(1\right) e^{- \left(\frac{1}{2}\right)} I_{2}\left(\left(\frac{1}{2}\right)\right)}{\left(\frac{1}{2}\right)}\simeq 0.30963
$$

である。

分母を払った二根の和は

$$
S=\left(\frac{16 I_{2}\left(\frac{1}{2}\right)}{e^{\frac{1}{2}}}\right) + \left(\frac{4 I_{1}\left(\frac{1}{2}\right)}{e^{\frac{1}{2}}}\right) + 5\simeq 5.9353
$$

である。

分母を払った二根の積は

$$
P_X=\left(\frac{16 I_{2}\left(\frac{1}{2}\right)}{e^{\frac{1}{2}}}\right) + 4 \left(\frac{4 I_{1}\left(\frac{1}{2}\right)}{e^{\frac{1}{2}}}\right) + 4\simeq 6.8124
$$

である。

第二極より上の形式根は

$$
X_+=\frac{\left(\frac{16 I_{2}\left(\frac{1}{2}\right)}{e^{\frac{1}{2}}} + \frac{4 I_{1}\left(\frac{1}{2}\right)}{e^{\frac{1}{2}}} + 5\right)}{2} + \frac{\sqrt{- 4 \left(\frac{16 I_{2}\left(\frac{1}{2}\right)}{e^{\frac{1}{2}}} + \frac{16 I_{1}\left(\frac{1}{2}\right)}{e^{\frac{1}{2}}} + 4\right) + \left(\left(\frac{16 I_{2}\left(\frac{1}{2}\right)}{e^{\frac{1}{2}}} + \frac{4 I_{1}\left(\frac{1}{2}\right)}{e^{\frac{1}{2}}} + 5\right)\right)^{2}}}{2}\simeq 4.38
$$

である。

第一極と第二極の間の根は

$$
X_-=\frac{\left(\frac{16 I_{2}\left(\frac{1}{2}\right)}{e^{\frac{1}{2}}} + \frac{16 I_{1}\left(\frac{1}{2}\right)}{e^{\frac{1}{2}}} + 4\right)}{\left(\frac{8 I_{2}\left(\frac{1}{2}\right)}{e^{\frac{1}{2}}} + \frac{2 I_{1}\left(\frac{1}{2}\right)}{e^{\frac{1}{2}}} + \frac{\sqrt{-16 - \frac{64 I_{1}\left(\frac{1}{2}\right)}{e^{\frac{1}{2}}} - \frac{64 I_{2}\left(\frac{1}{2}\right)}{e^{\frac{1}{2}}} + \left(\frac{16 I_{2}\left(\frac{1}{2}\right)}{e^{\frac{1}{2}}} + \frac{4 I_{1}\left(\frac{1}{2}\right)}{e^{\frac{1}{2}}} + 5\right)^{2}}}{2} + \frac{5}{2}\right)}\simeq 1.5553
$$

である。

二次の留数が増え、上の根は第二極から離れる。
:::

:::example[密度が低い場合]
両方の留数を一様に弱める。

入力値を次のように置く。

$$
\begin{aligned}
b&=\frac{1}{2} \\P&=\frac{1}{10}
\end{aligned}
$$

これらを順に代入する。

一次の留数は

$$
A_1=\frac{2 \left(\frac{1}{10}\right) e^{- \left(\frac{1}{2}\right)} I_{1}\left(\left(\frac{1}{2}\right)\right)}{\left(\frac{1}{2}\right)}\simeq 0.062568
$$

である。

二次の留数は

$$
A_2=\frac{8 \left(\frac{1}{10}\right) e^{- \left(\frac{1}{2}\right)} I_{2}\left(\left(\frac{1}{2}\right)\right)}{\left(\frac{1}{2}\right)}\simeq 0.030963
$$

である。

分母を払った二根の和は

$$
S=\left(\frac{2 I_{1}\left(\frac{1}{2}\right)}{5 e^{\frac{1}{2}}}\right) + \left(\frac{8 I_{2}\left(\frac{1}{2}\right)}{5 e^{\frac{1}{2}}}\right) + 5\simeq 5.0935
$$

である。

分母を払った二根の積は

$$
P_X=4 \left(\frac{2 I_{1}\left(\frac{1}{2}\right)}{5 e^{\frac{1}{2}}}\right) + \left(\frac{8 I_{2}\left(\frac{1}{2}\right)}{5 e^{\frac{1}{2}}}\right) + 4\simeq 4.2812
$$

である。

第二極より上の形式根は

$$
X_+=\frac{\left(\frac{8 I_{2}\left(\frac{1}{2}\right)}{5 e^{\frac{1}{2}}} + \frac{2 I_{1}\left(\frac{1}{2}\right)}{5 e^{\frac{1}{2}}} + 5\right)}{2} + \frac{\sqrt{\left(\left(\frac{8 I_{2}\left(\frac{1}{2}\right)}{5 e^{\frac{1}{2}}} + \frac{2 I_{1}\left(\frac{1}{2}\right)}{5 e^{\frac{1}{2}}} + 5\right)\right)^{2} - 4 \left(\frac{8 I_{2}\left(\frac{1}{2}\right)}{5 e^{\frac{1}{2}}} + \frac{8 I_{1}\left(\frac{1}{2}\right)}{5 e^{\frac{1}{2}}} + 4\right)}}{2}\simeq 4.0316
$$

である。

第一極と第二極の間の根は

$$
X_-=\frac{\left(\frac{8 I_{2}\left(\frac{1}{2}\right)}{5 e^{\frac{1}{2}}} + \frac{8 I_{1}\left(\frac{1}{2}\right)}{5 e^{\frac{1}{2}}} + 4\right)}{\left(\frac{4 I_{2}\left(\frac{1}{2}\right)}{5 e^{\frac{1}{2}}} + \frac{I_{1}\left(\frac{1}{2}\right)}{5 e^{\frac{1}{2}}} + \frac{\sqrt{-16 - \frac{32 I_{1}\left(\frac{1}{2}\right)}{5 e^{\frac{1}{2}}} - \frac{32 I_{2}\left(\frac{1}{2}\right)}{5 e^{\frac{1}{2}}} + \left(\frac{8 I_{2}\left(\frac{1}{2}\right)}{5 e^{\frac{1}{2}}} + \frac{2 I_{1}\left(\frac{1}{2}\right)}{5 e^{\frac{1}{2}}} + 5\right)^{2}}}{2} + \frac{5}{2}\right)}\simeq 1.0619
$$

である。

二根はそれぞれの極へ近づく。
:::

:::example[密度が高い場合]
上側根が次の高調波をどれほど離れるか見る。

入力値を次のように置く。

$$
\begin{aligned}
b&=1 \\P&=10
\end{aligned}
$$

これらを順に代入する。

一次の留数は

$$
A_1=\frac{2 \left(10\right) e^{- \left(1\right)} I_{1}\left(\left(1\right)\right)}{\left(1\right)}\simeq 4.1582
$$

である。

二次の留数は

$$
A_2=\frac{8 \left(10\right) e^{- \left(1\right)} I_{2}\left(\left(1\right)\right)}{\left(1\right)}\simeq 3.9951
$$

である。

分母を払った二根の和は

$$
S=\left(\frac{20 I_{1}\left(1\right)}{e}\right) + \left(\frac{80 I_{2}\left(1\right)}{e}\right) + 5\simeq 13.153
$$

である。

分母を払った二根の積は

$$
P_X=4 \left(\frac{20 I_{1}\left(1\right)}{e}\right) + \left(\frac{80 I_{2}\left(1\right)}{e}\right) + 4\simeq 24.628
$$

である。

第二極より上の形式根は

$$
X_+=\frac{\left(\frac{80 I_{2}\left(1\right)}{e} + \frac{20 I_{1}\left(1\right)}{e} + 5\right)}{2} + \frac{\sqrt{- 4 \left(\frac{80 I_{2}\left(1\right)}{e} + 4 + \frac{80 I_{1}\left(1\right)}{e}\right) + \left(\frac{80 I_{2}\left(1\right)}{e} + \frac{20 I_{1}\left(1\right)}{e} + 5\right)^{2}}}{2}\simeq 10.892
$$

である。

第一極と第二極の間の根は

$$
X_-=\frac{\left(\frac{80 I_{2}\left(1\right)}{e} + 4 + \frac{80 I_{1}\left(1\right)}{e}\right)}{\left(\frac{40 I_{2}\left(1\right)}{e} + \frac{10 I_{1}\left(1\right)}{e} + \frac{5}{2} + \frac{\sqrt{- \frac{320 I_{1}\left(1\right)}{e} - 16 - \frac{320 I_{2}\left(1\right)}{e} + \left(\frac{80 I_{2}\left(1\right)}{e} + \frac{20 I_{1}\left(1\right)}{e} + 5\right)^{2}}}{2}\right)}\simeq 2.261
$$

である。

二極近似で第三調波付近かそれ以上に移る値を実際の枝と同一視できない。
:::


## ドリフト運動論

$\omega/\abs\Omega\ll1$、$k_\perp\rho\ll1$ なら、旋回位相の分解能は不要である。位相を平均した分布 $F(\vb R,v_\parallel,\mu,t)$ は、:eqref[eq-plasma-gc-first] の特性で
$$
\partial_tF+\dot{\vb R}\cdot\nabla F+
\frac{qE_\parallel-\mu\nabla_\parallel B}{m}\partial_{v_\parallel}F=C[F]
$$
を満たす。これを :term[ドリフト運動論]{en="drift kinetics" yomi="どりふとうんどうろん"} と呼ぶ。ここで $\dot\mu=0$ を使った。速度の積分測度は $2\pi v_\perp\dd v_\perp\dd v_\parallel=(2\pi B/m)\dd\mu\dd v_\parallel$ である。曲がった場の完全な Hamilton 形式では、これに案内中心座標の Jacobian 補正が入る。

:::example[等温の鏡内平衡]
$\phi=0$ で $F\propto\exp[-(mv_\parallel^2/2+\mu B)/T]$ とする。平行移流の寄与 $-Fv_\parallel\mu \partial_\parallel B/T$ と鏡力の寄与は相殺する。
:::

:::example[速度積分の変数変換]
一様場で $F=A\exp[-(m v_\parallel^2/2+\mu B)/T]$ を積分する。$\mu$ の積分は $T/B$、平行のガウス積分は $\sqrt{2\pi T/m}$ であり、測度の $2\pi B/m$ を掛けて $n=A(2\pi T/m)^{3/2}$ となる。$B$ の因子が相殺されるので一様温度のこの平衡に磁場の強度だけで密度勾配は生じない。
:::

### 鏡力を含むドリフト運動論の特性

磁力線座標 $s$ に沿い $B(s)=B_0(1+s^2/L^2)$、$E_\parallel=0$ とする。一定の $\mu$ に対し $m\ddot s=-2\mu B_0s/L^2$ だから往復率は $\omega_b=\sqrt{2\mu B_0/(mL^2)}$ である。分布 $F(s,v,\mu,t)$ はこの振動子の逆時間写像で初期値を評価すればよい。例えば $F_0$ が $s,v$ の中心を持つ Gaussian なら、共移動の楕円が保たれる。速度積分の測度は物理位置で $B(s)$ に比例するので、この調和振動子の位相図だけから物理密度を読まない。磁束管断面が $B^{-1}$ に比例する場合、断面積と速度測度の磁場因子は相殺する。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $X_0$ | 初期位置 $s_0/L$ | 1 |
| $V_0$ | 初期規格化平行速度 $v_0/(\omega_b L)$ | 1 |
| $u$ | 往復位相 $\omega_b t$ | 1 |

現在の位置を $X$ とする。振動子の特性解である。

$$
X=V_{0} \sin{\left(u \right)} + X_{0} \cos{\left(u \right)}
$$

現在の平行速度を $V$ とする。位置を位相で微分する。

$$
V=V_{0} \cos{\left(u \right)} - X_{0} \sin{\left(u \right)}
$$

特性が通る磁場を $B/B_0$ とする。放物形の鏡分布を評価する。

$$
B/B_0=X^{2} + 1
$$

特性上保存するエネルギーを $H/(\mu B_0)$ とする。$\omega_b$ の定義から平行運動の係数が一になる。

$$
H/(\mu B_0)=V^{2} + X^{2} + 1
$$

反射点の絶対距離を $s_{\max}/L$ とする。反射点では平行速度がゼロである。

$$
s_{\max}/L=\sqrt{V_{0}^{2} + X_{0}^{2}}
$$

平行運動のエネルギー分率を $K_\parallel/H$ とする。残りは $\mu B$ である。

$$
K_\parallel/H=\frac{V^{2}}{H/(\mu B_{0)}}
$$

:::example[磁気面の中心から出発]
最初の反射点とエネルギーを求める。

入力値を次のように置く。

$$
\begin{aligned}
X_0&=0 \\V_0&=1 \\u&=\frac{\pi}{2}
\end{aligned}
$$

これらを順に代入する。

現在の位置は

$$
X=\left(0\right) \cos{\left(\left(\frac{\pi}{2}\right) \right)} + \left(1\right) \sin{\left(\left(\frac{\pi}{2}\right) \right)}=1
$$

である。

現在の平行速度は

$$
V=- \left(0\right) \sin{\left(\left(\frac{\pi}{2}\right) \right)} + \left(1\right) \cos{\left(\left(\frac{\pi}{2}\right) \right)}=0
$$

である。

特性が通る磁場は

$$
B/B_0=\left(1\right)^{2} + 1=2
$$

である。

特性上保存するエネルギーは

$$
H/(\mu B_0)=\left(0\right)^{2} + \left(1\right)^{2} + 1=2
$$

である。

反射点の絶対距離は

$$
s_{\max}/L=\sqrt{\left(0\right)^{2} + \left(1\right)^{2}}=1
$$

である。

平行運動のエネルギー分率は

$$
K_\parallel/H=\frac{\left(0\right)^{2}}{\left(2\right)}=0
$$

である。

反射点では平行運動エネルギーが横の $\mu B$ に替わる。
:::

:::example[端から中心へ]
反射点の初期条件で特性を追う。

入力値を次のように置く。

$$
\begin{aligned}
X_0&=1 \\V_0&=0 \\u&=\frac{\pi}{2}
\end{aligned}
$$

これらを順に代入する。

現在の位置は

$$
X=\left(0\right) \sin{\left(\left(\frac{\pi}{2}\right) \right)} + \left(1\right) \cos{\left(\left(\frac{\pi}{2}\right) \right)}=0
$$

である。

現在の平行速度は

$$
V=\left(0\right) \cos{\left(\left(\frac{\pi}{2}\right) \right)} - \left(1\right) \sin{\left(\left(\frac{\pi}{2}\right) \right)}=-1
$$

である。

特性が通る磁場は

$$
B/B_0=\left(0\right)^{2} + 1=1
$$

である。

特性上保存するエネルギーは

$$
H/(\mu B_0)=\left(-1\right)^{2} + \left(0\right)^{2} + 1=2
$$

である。

反射点の絶対距離は

$$
s_{\max}/L=\sqrt{\left(0\right)^{2} + \left(1\right)^{2}}=1
$$

である。

平行運動のエネルギー分率は

$$
K_\parallel/H=\frac{\left(-1\right)^{2}}{\left(2\right)}=\frac{1}{2}
$$

である。

速度の符号が負となり、赤道相当の中心を横切る。
:::

:::example[小さな振幅]
鏡内の位置と磁場差を求める。

入力値を次のように置く。

$$
\begin{aligned}
X_0&=\frac{1}{10} \\V_0&=0 \\u&=\frac{\pi}{4}
\end{aligned}
$$

これらを順に代入する。

現在の位置は

$$
X=\left(0\right) \sin{\left(\left(\frac{\pi}{4}\right) \right)} + \left(\frac{1}{10}\right) \cos{\left(\left(\frac{\pi}{4}\right) \right)}\simeq 0.070711
$$

である。

現在の平行速度は

$$
V=\left(0\right) \cos{\left(\left(\frac{\pi}{4}\right) \right)} - \left(\frac{1}{10}\right) \sin{\left(\left(\frac{\pi}{4}\right) \right)}\simeq -0.070711
$$

である。

特性が通る磁場は

$$
B/B_0=\left(\frac{\sqrt{2}}{20}\right)^{2} + 1=\frac{201}{200}
$$

である。

特性上保存するエネルギーは

$$
H/(\mu B_0)=\left(- \frac{\sqrt{2}}{20}\right)^{2} + \left(\frac{\sqrt{2}}{20}\right)^{2} + 1=\frac{101}{100}
$$

である。

反射点の絶対距離は

$$
s_{\max}/L=\sqrt{\left(0\right)^{2} + \left(\frac{1}{10}\right)^{2}}=\frac{1}{10}
$$

である。

平行運動のエネルギー分率は

$$
K_\parallel/H=\frac{\left(- \frac{\sqrt{2}}{20}\right)^{2}}{\left(\frac{101}{100}\right)}=\frac{1}{202}
$$

である。

放物形近似が実際の滑らかな鏡に適用しやすい領域である。
:::

:::example[位置と速度の両方が有限]
エネルギー交換で反射点へ到達するか調べる。

入力値を次のように置く。

$$
\begin{aligned}
X_0&=1 \\V_0&=1 \\u&=\frac{\pi}{4}
\end{aligned}
$$

これらを順に代入する。

現在の位置は

$$
X=\left(1\right) \sin{\left(\left(\frac{\pi}{4}\right) \right)} + \left(1\right) \cos{\left(\left(\frac{\pi}{4}\right) \right)}\simeq 1.4142
$$

である。

現在の平行速度は

$$
V=- \left(1\right) \sin{\left(\left(\frac{\pi}{4}\right) \right)} + \left(1\right) \cos{\left(\left(\frac{\pi}{4}\right) \right)}=0
$$

である。

特性が通る磁場は

$$
B/B_0=\left(\sqrt{2}\right)^{2} + 1=3
$$

である。

特性上保存するエネルギーは

$$
H/(\mu B_0)=\left(0\right)^{2} + \left(\sqrt{2}\right)^{2} + 1=3
$$

である。

反射点の絶対距離は

$$
s_{\max}/L=\sqrt{2} \sqrt{\left(1\right)^{2}}\simeq 1.4142
$$

である。

平行運動のエネルギー分率は

$$
K_\parallel/H=\frac{\left(0\right)^{2}}{\left(3\right)}=0
$$

である。

位置と速度の初期位相を両方使うと反射点が得られる。
:::


## ジャイロ運動論の概観

乱流の垂直尺度がイオンの旋回半径ほどでも、周波数が旋回周波数より低ければ高速位相を消去できる。$k_\perp\rho_i=\order{1}$、$k_\parallel/k_\perp\sim\omega/\abs{\Omega_i}\sim\rho_i/L\ll1$ を典型的な順序とする。:term[ジャイロ運動論]{en="gyrokinetics" yomi="じゃいろうんどうろん"} では、同じ電位でも旋回円上の値を平均して使う。

平面波 $\phi_{\vb k}$ の平均は
$$
\langle\phi\rangle_{\vb R}
=\frac1{2\pi}\int_0^{2\pi}\phi_{\vb k}e^{i\vb k\cdot(\vb R+\vb\varrho(\vartheta))}\dd\vartheta
=J_0(k_\perp\rho)\phi_{\vb k}e^{i\vb k\cdot\vb R}.
$$
均一な磁場、静電、Maxwell 平衡で、断熱応答を引いた分布を $h$ とすると、線形の式は $(\partial_t+v_\parallel\partial_\parallel)h=(q F_0/T)\partial_t\langle\phi\rangle_{\vb R}$ である。これは静電 Vlasov 方程式で粒子応答を $-(q\phi/T)F_0+h$ と書き、最低次の旋回平均を取って得られる。粒子密度への写し戻しにも $J_0$ が出るため、$\omega/( \omega-k_\parallel v_\parallel)$ が速度に依存しない極限では積分に $\Gamma_0$ が現れる。残りの $1-\Gamma_0$ が有限旋回半径の分極である。非線形・電磁的な展開では、分布と場を変分原理から求め、保存則を同じ次数にそろえる必要がある。その導出は :cite[bellan-plasma] の案内中心の方法を出発点とし、本章では一様場の静電応答に限る。

:::example[長波長極限]
$\Gamma_0=e^{-b}I_0(b)=1-b+\order{b^2}$ なので $1-\Gamma_0$ は $k_\perp^2$ に比例する。$k_\perp=0$ ではこの分極は消える。
:::

:::example[分極の数値評価]
$b=1$ なら $\Gamma_0=e^{-1}I_0(1)\simeq0.466$、$1-\Gamma_0\simeq0.534$ である。長波長近似 $1-\Gamma_0\simeq b$ はこの値をほぼ二倍に見積もる。垂直波長が旋回半径と同程度では因子を残して計算する。
:::

### ジャイロ運動論で分極をどこまで残すか

Maxwell の旋回平均で $\Gamma_0(b)=e^{-b}I_0(b)$ である。密度への写し戻しで $J_0$ がもう一つ掛かるため、粒子と案内中心を区別した分極は $1-\Gamma_0$ である。$I_0=1+b^2/4+\order{b^4}$ を指数と掛けると $1-\Gamma_0=b-3b^2/4+\order{b^3}$ となる。これに対し、与えられた一粒子の電位平均の低波数展開は $J_0(a)=1-a^2/4+\order{a^4}$ である。$a$ は粒子の速度を含み、$b$ は温度を含むため、$J_0$ の二次項を直接 $\Gamma_0$ と同一視しない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $b$ | 垂直波数と熱旋回半径の積の二乗 | 1 |
| $K$ | $k_\parallel/k_\perp$ の大きさ | 1 |
| $W$ | $\omega/\abs{\Omega}$ の大きさ | 1 |

Maxwell の二重旋回平均を $\Gamma_0$ とする。一粒子因子を二乗して熱平均する。

$$
\Gamma_0=e^{- b} I_{0}\left(b\right)
$$

正確な分極係数を $P$ とする。一から平均を引く。

$$
P=1 - \Gamma_{0}
$$

長波長での先頭項を $P_1$ とする。有限半径の二乗に比例する。

$$
P_1=b
$$

長波長での第二項込みの値を $P_2$ とする。$b$ が一を越えて大きい場合には延長できない。

$$
P_2=- \frac{3 b^{2}}{4} + b
$$

先頭項の相対誤差を $\delta_1$ とする。分極が小さくても相対誤差を比較する。

$$
\delta_1=\left|{1 - \frac{P_{1}}{P}}\right|
$$

時間と空間の異方性の尺度比を $W/K$ とする。低周波と長い平行尺度を別々に確認する。

$$
W/K=\frac{W}{K}
$$

:::example[ドリフト運動論と共通の領域]
分極の先頭項の誤差を求める。

入力値を次のように置く。

$$
\begin{aligned}
b&=\frac{1}{100} \\K&=\frac{1}{100} \\W&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

Maxwell の二重旋回平均は

$$
\Gamma_0=e^{- \left(\frac{1}{100}\right)} I_{0}\left(\left(\frac{1}{100}\right)\right)\simeq 0.99007
$$

である。

正確な分極係数は

$$
P=1 - \left(\frac{I_{0}\left(\frac{1}{100}\right)}{e^{\frac{1}{100}}}\right)\simeq 0.0099254
$$

である。

長波長での先頭項は

$$
P_1=\left(\frac{1}{100}\right)=\frac{1}{100}
$$

である。

長波長での第二項込みの値は

$$
P_2=- \frac{3 \left(\frac{1}{100}\right)^{2}}{4} + \left(\frac{1}{100}\right)=\frac{397}{40000}
$$

である。

先頭項の相対誤差は

$$
\delta_1=\left|{1 - \frac{\left(\frac{1}{100}\right)}{\left(- \frac{I_{0}\left(\frac{1}{100}\right)}{e^{\frac{1}{100}}} + 1\right)}}\right|\simeq 0.0075146
$$

である。

時間と空間の異方性の尺度比は

$$
W/K=1
$$

である。

両方の順序で軌道を平均でき、有限半径補正も小さい。
:::

:::example[有限旋回半径]
ジャイロ平均が必要な領域を評価する。

入力値を次のように置く。

$$
\begin{aligned}
b&=1 \\K&=\frac{1}{100} \\W&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

Maxwell の二重旋回平均は

$$
\Gamma_0=e^{- \left(1\right)} I_{0}\left(\left(1\right)\right)\simeq 0.46576
$$

である。

正確な分極係数は

$$
P=1 - \left(\frac{I_{0}\left(1\right)}{e}\right)\simeq 0.53424
$$

である。

長波長での先頭項は

$$
P_1=\left(1\right)=1
$$

である。

長波長での第二項込みの値は

$$
P_2=- \frac{3 \left(1\right)^{2}}{4} + \left(1\right)=\frac{1}{4}
$$

である。

先頭項の相対誤差は

$$
\delta_1=\left|{1 - \frac{\left(1\right)}{\left(- \frac{I_{0}\left(1\right)}{e} + 1\right)}}\right|\simeq 0.87182
$$

である。

時間と空間の異方性の尺度比は

$$
W/K=1
$$

である。

先頭の長波長式だけでは分極を大幅に過大評価する。
:::

:::example[並行順序を破る例]
時間が遅くても典型的な異方性にあるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
b&=1 \\K&=1 \\W&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

Maxwell の二重旋回平均は

$$
\Gamma_0=e^{- \left(1\right)} I_{0}\left(\left(1\right)\right)\simeq 0.46576
$$

である。

正確な分極係数は

$$
P=1 - \left(\frac{I_{0}\left(1\right)}{e}\right)\simeq 0.53424
$$

である。

長波長での先頭項は

$$
P_1=\left(1\right)=1
$$

である。

長波長での第二項込みの値は

$$
P_2=- \frac{3 \left(1\right)^{2}}{4} + \left(1\right)=\frac{1}{4}
$$

である。

先頭項の相対誤差は

$$
\delta_1=\left|{1 - \frac{\left(1\right)}{\left(- \frac{I_{0}\left(1\right)}{e} + 1\right)}}\right|\simeq 0.87182
$$

である。

時間と空間の異方性の尺度比は

$$
W/K=\frac{\left(\frac{1}{100}\right)}{\left(1\right)}=\frac{1}{100}
$$

である。

平行方向が十分長くないため典型的な低周波の順序とは異なる。
:::

:::example[高周波を含む例]
分極因子を計算できることとジャイロ運動論の適用を区別する。

入力値を次のように置く。

$$
\begin{aligned}
b&=1 \\K&=\frac{1}{100} \\W&=1
\end{aligned}
$$

これらを順に代入する。

Maxwell の二重旋回平均は

$$
\Gamma_0=e^{- \left(1\right)} I_{0}\left(\left(1\right)\right)\simeq 0.46576
$$

である。

正確な分極係数は

$$
P=1 - \left(\frac{I_{0}\left(1\right)}{e}\right)\simeq 0.53424
$$

である。

長波長での先頭項は

$$
P_1=\left(1\right)=1
$$

である。

長波長での第二項込みの値は

$$
P_2=- \frac{3 \left(1\right)^{2}}{4} + \left(1\right)=\frac{1}{4}
$$

である。

先頭項の相対誤差は

$$
\delta_1=\left|{1 - \frac{\left(1\right)}{\left(- \frac{I_{0}\left(1\right)}{e} + 1\right)}}\right|\simeq 0.87182
$$

である。

時間と空間の異方性の尺度比は

$$
W/K=\frac{\left(1\right)}{\left(\frac{1}{100}\right)}=100
$$

である。

同じ円平均でも、時間分離が破れているので旋回高調波を除去できない。
:::


### 一様場の非線形静電系と自由エネルギー

周期箱に一様な背景磁場 $B\hat{\vb z}$、位置によらない Maxwell 平衡 $F_{s0}$ を置く。場の微小さを $q_s\phi/T_s=\order{\epsilon}$、異方性を $k_\parallel/k_\perp=\order{\epsilon}$、対象とする横波数を $k_\perp\rho_s=\order{1}$ とする。ここでは電磁摂動、背景勾配、衝突を省く。$\Phi_s=\langle\phi\rangle_{\vb R}$ とし、案内中心の粒子分布の摂動を $h_s-(q_sF_{s0}/T_s)\Phi_s$ と表す。横の移流は $\vb v_E(\Phi_s)\cdot\grad h_s=\{\Phi_s,h_s\}/B$ である。ただし $\{a,b\}=\partial_Xa\,\partial_Y b-\partial_Ya\,\partial_Xb$ と定義する。

平行の電気力の最低次は平衡分布だけに作用し、$(q_s/m_s)\partial_z\Phi_s\,\partial_{v_\parallel}F_{s0}=-(q_sF_{s0}/T_s)v_\parallel\partial_z\Phi_s$ を与える。これは断熱部分の平行移流と相殺する。また、横の移流が断熱部分に作用すると $\{\Phi_s,\Phi_s\}=0$ を得る。平行加速度が摂動分布にも作用する項はさらに $\order{\epsilon}$ 小さいので、ここで保持する次数では
$$
\partial_t h_s+v_\parallel\partial_z h_s
+\frac{1}{B}\{\Phi_s,h_s\}
=\frac{q_sF_{s0}}{T_s}\partial_t\Phi_s
\label{eq-plasma-slab-gyrokinetic}
$$
となる。横移流は、横の勾配に $k_\perp\rho_s=\order{1}$ が入るため、時間発展と同じ次数で残る。

粒子密度に変換するには、位置 $\vb r$ を固定して $h_s$ を旋回平均する。その平均を $\langle h_s\rangle_{\vb r}$ と書く。Debye 長を無視できる準中性近似では
$$
\sum_s q_s\int \langle h_s\rangle_{\vb r}\,\dd^3v
=\sum_s\frac{q_s^2 n_s}{T_s}\phi .
$$
周期境界のもとで粒子位置と案内中心位置を平行移動して積分変数を変えると、二種類の旋回平均は互いに随伴である。従って :eqref[eq-plasma-slab-gyrokinetic] を $T_s h_s/F_{s0}$ で重み付けし、空間と速度で積分した右辺は $\int\partial_t\phi \sum_s q_s\int\langle h_s\rangle_{\vb r}\dd^3v\,\dd^3r$ になる。左辺の平行・垂直移流項は、周期性と $\divg\vb v_E=0$ から積分して消える。種を足して準中性式を代入すれば
$$
\mathcal W=
\sum_s\int\frac{T_s h_s^2}{2F_{s0}}\,\dd^3R\,\dd^3v
-\frac12\sum_s\frac{q_s^2n_s}{T_s}\int\phi^2\,\dd^3r
\label{eq-plasma-gyrokinetic-free-energy}
$$
が時間に依存しない。

式 :eqref[eq-plasma-gyrokinetic-free-energy] の非負性は、同じ準中性式によって示せる。$g_s=h_s-(q_sF_{s0}/T_s)\Phi_s$ と置き、二乗を展開する。$\Phi_s$ の二乗を速度で平均すると各 Fourier モードの $\Gamma_0(b_s)|\phi_{\vb k}|^2$ が得られる。準中性式に現れる $g_s$ の速度積分は $\sum_s(q_s^2n_s/T_s)(1-\Gamma_0)\phi_{\vb k}$ である。展開の交差項にこの積分を入れると、規格化された Fourier 係数で
$$
\mathcal W=
\sum_s\int\frac{T_s g_s^2}{2F_{s0}}\,\dd^3R\,\dd^3v
+\frac12\sum_{\vb k,s}\frac{q_s^2 n_s}{T_s}
[1-\Gamma_0(b_s)]|\phi_{\vb k}|^2
$$
となる。$\Gamma_0$ は $J_0^2$ の Maxwell 平均であり、円周平均の絶対値が一以下なので $0\le\Gamma_0\le1$ である。従って両項は非負である。この保存則は、電位が減っても速度分布の細かな構造にエネルギーが残ることを表せる。

:::example[磁場に平行な一様 Fourier 成分]
横一様な成分に限る低周波・長波長極限では全種で $\Phi_s=\phi$ かつ $\Gamma_0=1$ となり、横分極の項はゼロである。有限の $k_\parallel$ を保つと、平行移流によって速度方向の位相は混合する。この極限では静電の場の決定に準中性条件を用いる。分極項がないことは $g_s$ の二乗積分がゼロであることを意味しない。
:::

:::example[一つの磁化種の分極項の比較]
同じ電位振幅で $b=1/10$ と $b=1$ を比較する。$\Gamma_0=e^{-b}I_0(b)$ より $1-\Gamma_0$ はそれぞれ約 $0.09290$、$0.53424$ である。その比は約 $5.75$ で、先頭の長波長式から得る十倍とは異なる。短い横尺度では電位のエネルギー規格化にも有限旋回半径が効く。
:::

:::example[自由エネルギーが転送される場合]
あるモードで速度分布二乗項が $3$、分極項が $2$ と同じエネルギー単位で評価されるとする。全モードでの分極項がその後 $1/2$ に減り、他のモードに初期エネルギーがなかったなら、全分布二乗項は $5-1/2=9/2$ である。式が保存するのは全モードでの和であり、非線形の横移流は各モードの値を保存しない。
:::

## 演習問題

::::exercise{level=basic}
$k_\parallel=2$、$\omega_r=10$、$n\Omega_s=6$ を SI で与える。共鳴速度を求めよ。

:::solution
$(\omega_r-n\Omega_s)/k_\parallel=2\,\mathrm{m\,s^{-1}}$ である。
:::
::::

::::exercise{level=basic}
$k_\perp=0$ の旋回平均因子を求めよ。

:::solution
周回位相 $e^{i\vb k\cdot\vb\varrho}=1$ なので $J_0(0)=1$ である。
:::
::::

::::exercise{level=basic}
$T_s,m_s,B_0$ を保ち $k_\perp$ を三倍にする。$b_s$ は何倍か。

:::solution
定義で $b_s\propto k_\perp^2$ なので九倍である。
:::
::::

::::exercise{level=basic}
磁気モーメントを固定して磁場が二倍になると $v_\perp$ は何倍か。

:::solution
$v_\perp^2=2\mu B/m$ なので $\sqrt2$ 倍である。
:::
::::

::::exercise{level=advanced}
一種の垂直静電誘電関数で $b\to0$ としたとき、旋回高調波の極に近づく枝を除いて根を求めよ。$\omega_p^2+\Omega^2$ は高調波の角周波数の二乗と一致しないとする。

:::solution
式 :eqref[eq-plasma-bernstein] で極から離れて $n=1$ のみ残る。$1-\omega_p^2/(\omega^2-\Omega^2)=0$ だから $\omega^2=\Omega^2+\omega_p^2$ である。
:::
::::

::::exercise{level=advanced}
$\mathsf A$ の磁気力部分が速度との仕事をしないことを示せ。

:::solution
成分で $v_i A_{ij}=(1-\vb k\cdot\vb v/\omega)v_j+(\vb v\cdot\vb k)v_j/\omega=v_j$ である。従って $\vb v\cdot(\mathsf A-I)\vb E=0$ である。
:::
::::

::::exercise{level=advanced}
$T_\parallel,T_\perp$ の二温度 Maxwell 分布に共鳴方向微分を作用させよ。

:::solution
$\partial_{v_\parallel}f=-(mv_\parallel/T_\parallel)f$、$\partial_{v_\perp}f=-(mv_\perp/T_\perp)f$ なので $-mf[(\omega_r-n\Omega_s)/T_\parallel+n\Omega_s/T_\perp]$ となる。
:::
::::

::::exercise{level=advanced}
$\Gamma_0$ の二次係数を求めよ。

:::solution
$I_0(b)=1+b^2/4+\order{b^4}$、$e^{-b}=1-b+b^2/2+\order{b^3}$ なので積の二次係数は $3/4$ である。
:::
::::

## 確認テスト

::::quiz
垂直 Bernstein 波の有限温度効果に直接関わる量を選べ。

- [x] $k_\perp\rho$
- [ ] $k_\parallel\rho$ だけ
- [ ] 熱速度と無関係な $k_\perp$ の値だけ

:::explanation
位相 $e^{i\vb k\cdot\vb\varrho}$ は旋回円で変化するので、平均に $J_n(k_\perp\rho)$ が現れる。
:::
::::

::::quiz{answer="4" tol="1e-6"}
$\omega_r=17\,\mathrm{s^{-1}},\,n\Omega=5\,\mathrm{s^{-1}},\,k_\parallel=3\,\mathrm{m^{-1}}$ の共鳴速度を SI で答えよ。

:::explanation
速度は周波数差を波数で割った $(17-5)/3=4\,\mathrm{m\,s^{-1}}$ である。
:::
::::

::::quiz
ジャイロ平均を取る代表的な順序を選べ。

- [x] $k_\perp\rho\sim 1$ かつ $\omega/\abs{\Omega}\ll 1$
- [ ] $\omega/\abs{\Omega}\gg 1$
- [ ] $\omega=\abs{\Omega}$ の共鳴

:::explanation
高速旋回を時間的に分離する。垂直波長は旋回半径に近くてよいため円周平均の Bessel 因子が残る。
:::
::::

::::quiz{answer="0.25" tol="1e-6"}
温度を変えず磁場を二倍にしたとき $b_s$ の倍率を答えよ。

:::explanation
$\Omega_s\propto B$、$b_s\propto\Omega_s^{-2}$ より $1/4$ である。
:::
::::

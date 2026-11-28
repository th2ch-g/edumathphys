---
title: 非線形プラズマ現象
description: 波と分布が互いに変化する過程を、捕捉、孤立波、乱流の模型で調べる。
notation:
  - symbol: 'D_v'
    meaning: 準線形速度拡散係数
  - symbol: '\omega_b'
    meaning: 捕捉粒子の往復角周波数
  - symbol: '\phi_1'
    meaning: 規格化電位の先頭摂動
---

## この章の目標

- 平均分布に対する波の速度拡散を導く。
- 捕捉と流体の弱非線形波を区別する。
- パラメトリック成長と無衝突衝撃波の収支を評価する。

## 準線形理論

線形理論では背景分布を固定する。多数の弱い波が長時間作用する場合、共鳴粒子の蓄積した変化で背景も変わる。空間または位相平均を $\langle\cdot\rangle$ とし、$f=F+f_1$、$\langle f_1\rangle=0$、$\langle E\rangle=0$ とする。一次元 Vlasov 式を平均すると $\partial_tF=-(q/m)\partial_v\langle E f_1\rangle$ である。波ごとの位相の相関を無視し、$F$ を相関時間中一定と置く。自由な軌道に沿う一次式から
$$
f_1(x,v,t)=-\frac qm\,\partial_vF
\int_0^\infty E(x-v\tau,t-\tau)\dd\tau
$$
を得る。ただし積分には小さい収束因子を入れる。従って :term[準線形理論]{en="quasilinear theory" yomi="じゅんせんけいりろん"} の分布式は
$$
\partial_tF=\partial_v(D_v\partial_vF),\qquad
D_v=\frac{q^2}{m^2}\int_0^\infty\langle E(x,t) E(x-v\tau,t-\tau)\rangle\dd\tau .
\label{eq-plasma-ql}
$$
実振幅 $E_j$ の余弦波をランダム位相で重ねれば、相関は $\sum_j |E_j|^2\cos[(\omega_j-k_jv)\tau]/2$ である。$\int_0^\infty e^{-\epsilon\tau}\cos(a\tau)\dd\tau=\epsilon/(\epsilon^2+a^2)\to\pi\delta(a)$ から $D_v=\pi q^2\sum_j |E_j|^2\delta(\omega_j-k_jv)/(2m^2)$ となる。幅を持つ滑らかな拡散係数として使うには、連続スペクトルや有限相関時間を入れる。単一の強い波では捕捉が起こり、直線軌道の積分は使えない。

:::example[平坦化と粒子数]
閉じた速度区間で無流束なら $\partial_t\int F\dd v=0$ である。$\dv{}{t}\int F^2\dd v=-2\int D_v(\partial_vF)^2\dd v\le0$ なので、$D_v>0$ の連結区間で定常分布は平坦になる。
:::

:::example[有限の相関時間]
単一周波数の相関に $e^{-\nu_c\tau}$ を掛けると、$\nu_c>0$ で $D_v=(q^2 |E|^2/(2m^2))\,\nu_c/[\nu_c^2+(\omega-kv)^2]$ である。半値になる速度差は $\nu_c/\abs k$ である。相関時間が短くなると共鳴幅は広がり、中心での拡散係数は小さくなる。
:::

### 波と粒子のエネルギー収支

振幅がゆっくり変わる正周波数の縦波を取り、$\omega_r$ を誘電関数の実部の単純零点とする。実振幅 $E_k$ の余弦波による種 $s$ の速度拡散を運動エネルギー密度に代入する。速度境界で流束をゼロとすれば
$$
\dot W_{{\rm p},s}
=-\int m_s vD_{v,s}\partial_v F_s\dd v
=-\frac{\pi q_s^2 |E_k|^2}{2m_s |k|}\frac{\omega_r}{k}
\left.\partial_v F_s\right|_{v=\omega_r/k}
$$
となる。ここで $F_s$ の速度積分は種の密度である。$k>0$ なら、この式は :eqref[eq-plasma-landau-rate] の種別の虚部 $\varepsilon_{I,s}$ を使って $\varepsilon_0\omega_r\varepsilon_{I,s}|E_k|^2/2$ と書ける。従って波のエネルギー密度を $W_{\rm w}=(\varepsilon_0/4)\omega_r(\partial_\omega\varepsilon_R)|E_k|^2$ と評価する次数では、$\partial_t|E_k|^2=2\gamma |E_k|^2$ と Landau の一次補正から $\dot W_{\rm w}=-\sum_s\dot W_{{\rm p},s}$ である。分布の変化に伴う $\omega_r$ と規格化の変化はここでは高次として落とした。外から固定したスペクトルによる拡散では、その外部入力もエネルギー収支に加える。

### 準線形拡散を一定係数で解く

狭い共鳴帯の内部で $D(v)=D_0>0$ と近似する。速度幅 $2L$ の区間 $v\in[-L,L]$ に反射端を置けば、速度流束 $-D_0\partial_vF$ は端でゼロである。$x=(v+L)/(2L)$ に対し $\cos(n\pi x)$ の端微分はゼロであるため、摂動 $F-\bar F$ の各係数は $\dot F_n=-D_0(n\pi/(2L))^2F_n$ を満たす。定数成分は保存し、細かな構造から消える。さらに二乗偏差の微分は部分積分により $-2D_0\int|\partial_vF|^2\dd v$ である。これは波スペクトルを固定して平坦化を見る模型である。実際には波エネルギーも分布の傾きに応じて変化する。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $n$ | ゼロでないモード番号 | 1 |
| $s$ | 時間 $D_0t/(2L)^2$ | 1 |
| $a$ | 初期モード振幅を平均分布で割った値 | 1 |
| $x$ | 評価する規格化速度 | 1 |

規格化速度空間波数を $K$ とする。Neumann 条件で離散化する。

$$
K=\pi n
$$

減衰指数を $Z$ とする。分布の二階微分から出る。

$$
Z=K^{2} s
$$

残るモード振幅を $A(t)$ とする。指数関数で減る。

$$
A(t)=a e^{- Z}
$$

局所速度分布を $F/\bar F$ とする。平均密度は保存する。

$$
F/\bar F=A(t) \cos{\left(K x \right)} + 1
$$

局所の速度勾配を $\partial_x(F/\bar F)$ とする。共鳴による波との交換の尺度である。

$$
\partial_x(F/\bar F)=- A(t) K \sin{\left(K x \right)}
$$

二乗偏差の残存率を $W/W_0$ とする。拡散はこの量を単調に減らす。

$$
W/W_0=e^{- 2 Z}
$$

:::example[大域の傾きの平坦化]
局所の値と勾配を別々に求める。

入力値を次のように置く。

$$
\begin{aligned}
n&=1 \\s&=\frac{1}{10} \\a&=\frac{1}{5} \\x&=\frac{1}{2}
\end{aligned}
$$

これらを順に代入する。

規格化速度空間波数は

$$
K=\pi \left(1\right)\simeq 3.1416
$$

である。

減衰指数は

$$
Z=\left(\frac{1}{10}\right) \left(\pi\right)^{2}\simeq 0.98696
$$

である。

残るモード振幅は

$$
A(t)=\left(\frac{1}{5}\right) e^{- \left(\frac{\pi^{2}}{10}\right)}\simeq 0.074542
$$

である。

局所速度分布は

$$
F/\bar F=\left(\frac{1}{5 e^{\frac{\pi^{2}}{10}}}\right) \cos{\left(\left(\frac{1}{2}\right) \left(\pi\right) \right)} + 1=1
$$

である。

局所の速度勾配は

$$
\partial_x(F/\bar F)=- \left(\frac{1}{5 e^{\frac{\pi^{2}}{10}}}\right) \left(\pi\right) \sin{\left(\left(\frac{1}{2}\right) \left(\pi\right) \right)}\simeq -0.23418
$$

である。

二乗偏差の残存率は

$$
W/W_0=e^{- 2 \left(\frac{\pi^{2}}{10}\right)}\simeq 0.13891
$$

である。

区間中央では値が平均に等しくても傾きが残る。
:::

:::example[同じ時間の細かい成分]
高次モードの減衰を比較する。

入力値を次のように置く。

$$
\begin{aligned}
n&=3 \\s&=\frac{1}{10} \\a&=\frac{1}{5} \\x&=\frac{1}{2}
\end{aligned}
$$

これらを順に代入する。

規格化速度空間波数は

$$
K=\pi \left(3\right)\simeq 9.4248
$$

である。

減衰指数は

$$
Z=\left(3 \pi\right)^{2} \left(\frac{1}{10}\right)\simeq 8.8826
$$

である。

残るモード振幅は

$$
A(t)=\left(\frac{1}{5}\right) e^{- \left(\frac{9 \pi^{2}}{10}\right)}\simeq 2.7755\times 10^{-5}
$$

である。

局所速度分布は

$$
F/\bar F=\left(\frac{1}{5 e^{\frac{9 \pi^{2}}{10}}}\right) \cos{\left(\left(3 \pi\right) \left(\frac{1}{2}\right) \right)} + 1=1
$$

である。

局所の速度勾配は

$$
\partial_x(F/\bar F)=- \left(3 \pi\right) \left(\frac{1}{5 e^{\frac{9 \pi^{2}}{10}}}\right) \sin{\left(\left(3 \pi\right) \left(\frac{1}{2}\right) \right)}\simeq 0.00026159
$$

である。

二乗偏差の残存率は

$$
W/W_0=e^{- 2 \left(\frac{9 \pi^{2}}{10}\right)}\simeq 1.9259\times 10^{-8}
$$

である。

速度の細かな構造は早く平坦になる。
:::

:::example[反射端の速度流束]
端に密度偏差がありながら流束がゼロになるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
n&=1 \\s&=\frac{1}{10} \\a&=\frac{1}{5} \\x&=0
\end{aligned}
$$

これらを順に代入する。

規格化速度空間波数は

$$
K=\pi \left(1\right)\simeq 3.1416
$$

である。

減衰指数は

$$
Z=\left(\frac{1}{10}\right) \left(\pi\right)^{2}\simeq 0.98696
$$

である。

残るモード振幅は

$$
A(t)=\left(\frac{1}{5}\right) e^{- \left(\frac{\pi^{2}}{10}\right)}\simeq 0.074542
$$

である。

局所速度分布は

$$
F/\bar F=\left(\frac{1}{5 e^{\frac{\pi^{2}}{10}}}\right) \cos{\left(\left(0\right) \left(\pi\right) \right)} + 1\simeq 1.0745
$$

である。

局所の速度勾配は

$$
\partial_x(F/\bar F)=- \left(\frac{1}{5 e^{\frac{\pi^{2}}{10}}}\right) \left(\pi\right) \sin{\left(\left(0\right) \left(\pi\right) \right)}=0
$$

である。

二乗偏差の残存率は

$$
W/W_0=e^{- 2 \left(\frac{\pi^{2}}{10}\right)}\simeq 0.13891
$$

である。

無流束条件は分布の値でなく微分を指定する。
:::

:::example[短い時間]
高次構造も拡散の初期に残るか判定する。

入力値を次のように置く。

$$
\begin{aligned}
n&=3 \\s&=\frac{1}{1000} \\a&=\frac{1}{5} \\x&=\frac{1}{2}
\end{aligned}
$$

これらを順に代入する。

規格化速度空間波数は

$$
K=\pi \left(3\right)\simeq 9.4248
$$

である。

減衰指数は

$$
Z=\left(3 \pi\right)^{2} \left(\frac{1}{1000}\right)\simeq 0.088826
$$

である。

残るモード振幅は

$$
A(t)=\left(\frac{1}{5}\right) e^{- \left(\frac{9 \pi^{2}}{1000}\right)}\simeq 0.183
$$

である。

局所速度分布は

$$
F/\bar F=\left(\frac{1}{5 e^{\frac{9 \pi^{2}}{1000}}}\right) \cos{\left(\left(3 \pi\right) \left(\frac{1}{2}\right) \right)} + 1=1
$$

である。

局所の速度勾配は

$$
\partial_x(F/\bar F)=- \left(3 \pi\right) \left(\frac{1}{5 e^{\frac{9 \pi^{2}}{1000}}}\right) \sin{\left(\left(3 \pi\right) \left(\frac{1}{2}\right) \right)}\simeq 1.7247
$$

である。

二乗偏差の残存率は

$$
W/W_0=e^{- 2 \left(\frac{9 \pi^{2}}{1000}\right)}\simeq 0.83723
$$

である。

指数は一より小さく、同じ波数でも時刻によって残存量が異なる。
:::


### 位相混合の反転とプラズマエコー

初期に電場を短く印加し、波数 $-k_1$ の分布摂動 $g_{-k_1}(v,0)$ を作る。電場の自己応答をいったん省いた自由飛行の段階では、時間 $\tau$ 後に $g_{-k_1}(v,0)e^{ik_1v\tau}$ となる。そこで別の短い電場を波数 $k_2>k_1>0$ で印加する。Vlasov 方程式の $E\partial_v g$ はこの二成分を掛け合わせるため、$k_3=k_2-k_1$ で $\partial_v[g_{-k_1}(v,0)e^{ik_1v\tau}]$ に比例する分布を作る。その後の自由飛行では $e^{-ik_3v(t-\tau)}$ が加わる。密度の速度積分に含まれる速度依存の位相は
$$
\exp\{iv[k_1\tau-k_3(t-\tau)]\}
$$
であり、$t_{\rm e}=k_2\tau/(k_2-k_1)$ の周辺で位相混合が戻る。これが :term[プラズマエコー]{en="plasma echo" yomi="ぷらずまえこー"} の時間を与える。速度微分や自己無撞着な電場も波形に効くため、その時刻が電場の絶対値の最大点になるとは限らない。例えば $k_2=3k_1$ なら $t_{\rm e}=3\tau/2$ であり、二度目の印加から $\tau/2$ 遅れる。自己無撞着な波形を求める場合は、差波数の二次の源を Landau の初期値応答に入れる。 :cite[bellan-plasma]

## 粒子捕捉と非線形 Landau 減衰

単色静電波 $\phi=\Phi\cos(kx-\omega t)$、$q\Phi$ の符号を取り込んだ位相でポテンシャルを $-\abs{q\Phi}\cos\theta$ とする。波の座標で $\theta=k(x-\omega t/k)$ の時間発展は
$$
\ddot\theta=-\omega_b^2\sin\theta,\qquad
\omega_b=|k|\sqrt{\abs{q\Phi}/m}.
\label{eq-plasma-bounce}
$$
:term[粒子捕捉]{en="particle trapping" yomi="りゅうしほそく"} は共鳴付近で波の井戸から脱出できない運動である。分離軌道のエネルギーは $+\abs{q\Phi}$、井戸底は $-\abs{q\Phi}$ なので底での速度半幅は $2\sqrt{\abs{q\Phi}/m}$ である。小振幅の往復周波数は $\omega_b$ だが、分離軌道付近の周期は長くなる。

線形減衰が粒子の分布を変えない近似は、振幅減衰時間内に捕捉運動が小さい場合に使える。$\omega_b/\abs{\gamma_L}$ が大きいと、捕捉粒子が速度勾配を変えて線形の減衰率が維持できない。実際の時間発展は位相混合、波振幅の変化、初期粒子分布に依存する。非線形 Landau 減衰の長時間定理には、分布の解析性や共鳴の連鎖の評価が必要である。この章では振り子模型によるしきい値の目安に限る。

:::example[振幅変化と往復]
一定 $k$ で $\Phi$ を九倍にすれば小振幅の往復周波数は三倍、分離軌道の半幅も三倍になる。
:::

### 分離軌道の周期

振り子の最大角を $\theta_m<\pi$、$\ddot\theta=-\omega_b^2\sin\theta$ とする。エネルギー保存から $\dot\theta^2=2\omega_b^2(\cos\theta-\cos\theta_m)$ なので周期は四分の一を積分して
$$
\mathcal T=\frac4{\omega_b}\int_0^{\pi/2}\frac{\dd\alpha}{\sqrt{1-\sin^2(\theta_m/2)\sin^2\alpha}}.
$$
ここで $\sin(\theta/2)=\sin(\theta_m/2)\sin\alpha$ を用いた。$\theta_m\to0$ で被積分関数は一に近づき、周期は $2\pi/\omega_b$ に戻る。$\theta_m\to\pi$ では $\alpha\to\pi/2$ で分母が $\cos\alpha$ に近づき、積分が対数的に発散する。このように捕捉域内でも位相混合を起こす周期差がある。

:::example[弱い振り子の有限振幅補正]
$k_b=\sin(\theta_m/2)$ が小さいと $(1-k_b^2\sin^2\alpha)^{-1/2}=1+k_b^2\sin^2\alpha/2+\cdots$ である。積分すると小振幅周期に $(1+k_b^2/4)$ が掛かる。大きい軌道は小さい軌道より少し遅れて巡る。
:::

### 捕捉条件と波の減衰の時間比較

電子が見る波の座標を $z=x-v_\phi t$ とし、位相を選んで $H=m\dot z^2/2-U\cos(kz)$ とする。$U=e|\phi|$ はエネルギーである。鞍点の $H=U$ が捕捉域の境界となり、井戸の中心での捕捉半幅は $2\sqrt{U/m}$ である。微小振動の角周波数は $\omega_b=k\sqrt{U/m}$ である。線形振幅が $e^{-g t}$ で減る初期段階では、初期 $\omega_{b0}$ を用いても実際の回転率は $\omega_{b0}e^{-gt/2}$ となる。これを時間積分すると位相の形式上限は $2\omega_{b0}/g$ である。既に捕捉で波が変化すればこの線形減衰の外挿は破れるため、この値を捕捉の厳密なしきい値にしない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $R$ | 初期往復角周波数と線形減衰率の比 $\omega_{b0}/g$ | 1 |
| $s$ | 時間 $gt$ | 1 |
| $X$ | 初期の速度を $\sqrt{U_0/m}$ で割った値 | 1 |
| $\alpha$ | 初期波の位相 $kz$ | 1 |

初期の規格化 Hamiltonianを $H_0/U_0$ とする。一未満なら固定初期振幅で捕捉される。

$$
H_0/U_0=\frac{X^{2}}{2} - \cos{\left(\alpha \right)}
$$

線形外挿の波振幅を $U/U_0$ とする。粒子の非線形な反作用を入れない。

$$
U/U_0=e^{- s}
$$

時刻での捕捉半幅を $\Delta v_{\rm tr}/\sqrt{U_0/m}$ とする。瞬間的な分離軌道を使う。

$$
\Delta v_{\rm tr}/\sqrt{U_0/m}=2 \sqrt{U/U_{0}}
$$

瞬間の往復率と減衰率の比を $\omega_b/g$ とする。振幅の平方根で変わる。

$$
\omega_b/g=R \sqrt{U/U_{0}}
$$

微小往復の線形外挿位相を $\Theta$ とする。時間変化する周波数を積分する。

$$
\Theta=2 R \left(1 - e^{- \frac{s}{2}}\right)
$$

外挿の往復回数を $\Theta/(2\pi)$ とする。分離軌道近傍の実周期は長い。

$$
\Theta/(2\pi)=\frac{\Theta}{2 \pi}
$$

:::example[弱い波の井戸中心近く]
初期のエネルギー条件と一往復の可能性を比較する。

入力値を次のように置く。

$$
\begin{aligned}
R&=\frac{1}{10} \\s&=1 \\X&=0 \\\alpha&=0
\end{aligned}
$$

これらを順に代入する。

初期の規格化 Hamiltonianは

$$
H_0/U_0=\frac{\left(0\right)^{2}}{2} - \cos{\left(\left(0\right) \right)}=-1
$$

である。

線形外挿の波振幅は

$$
U/U_0=e^{- \left(1\right)}\simeq 0.36788
$$

である。

時刻での捕捉半幅は

$$
\Delta v_{\rm tr}/\sqrt{U_0/m}=2 \sqrt{\left(e^{-1}\right)}\simeq 1.2131
$$

である。

瞬間の往復率と減衰率の比は

$$
\omega_b/g=\left(\frac{1}{10}\right) \sqrt{\left(e^{-1}\right)}\simeq 0.060653
$$

である。

微小往復の線形外挿位相は

$$
\Theta=2 \left(\frac{1}{10}\right) \left(1 - e^{- \frac{\left(1\right)}{2}}\right)\simeq 0.078694
$$

である。

外挿の往復回数は

$$
\Theta/(2\pi)=\frac{\left(\frac{1}{5} - \frac{1}{5 e^{\frac{1}{2}}}\right)}{2 \pi}\simeq 0.012525
$$

である。

固定場では捕捉点でも、波が速く消える場合は捕捉による緩和を期待しにくい。
:::

:::example[強い波と初期捕捉]
初期のエネルギーと往復の時間を求める。

入力値を次のように置く。

$$
\begin{aligned}
R&=10 \\s&=1 \\X&=\frac{1}{2} \\\alpha&=0
\end{aligned}
$$

これらを順に代入する。

初期の規格化 Hamiltonianは

$$
H_0/U_0=\frac{\left(\frac{1}{2}\right)^{2}}{2} - \cos{\left(\left(0\right) \right)}=- \frac{7}{8}
$$

である。

線形外挿の波振幅は

$$
U/U_0=e^{- \left(1\right)}\simeq 0.36788
$$

である。

時刻での捕捉半幅は

$$
\Delta v_{\rm tr}/\sqrt{U_0/m}=2 \sqrt{\left(e^{-1}\right)}\simeq 1.2131
$$

である。

瞬間の往復率と減衰率の比は

$$
\omega_b/g=\left(10\right) \sqrt{\left(e^{-1}\right)}\simeq 6.0653
$$

である。

微小往復の線形外挿位相は

$$
\Theta=2 \left(10\right) \left(1 - e^{- \frac{\left(1\right)}{2}}\right)\simeq 7.8694
$$

である。

外挿の往復回数は

$$
\Theta/(2\pi)=\frac{\left(20 - \frac{20}{e^{\frac{1}{2}}}\right)}{2 \pi}\simeq 1.2525
$$

である。

外挿では往復位相が大きくなり、線形減衰の前提を再検討する。
:::

:::example[井戸の外側]
往復率が大きくても粒子が捕捉されるか確かめる。

入力値を次のように置く。

$$
\begin{aligned}
R&=10 \\s&=1 \\X&=3 \\\alpha&=0
\end{aligned}
$$

これらを順に代入する。

初期の規格化 Hamiltonianは

$$
H_0/U_0=\frac{\left(3\right)^{2}}{2} - \cos{\left(\left(0\right) \right)}=\frac{7}{2}
$$

である。

線形外挿の波振幅は

$$
U/U_0=e^{- \left(1\right)}\simeq 0.36788
$$

である。

時刻での捕捉半幅は

$$
\Delta v_{\rm tr}/\sqrt{U_0/m}=2 \sqrt{\left(e^{-1}\right)}\simeq 1.2131
$$

である。

瞬間の往復率と減衰率の比は

$$
\omega_b/g=\left(10\right) \sqrt{\left(e^{-1}\right)}\simeq 6.0653
$$

である。

微小往復の線形外挿位相は

$$
\Theta=2 \left(10\right) \left(1 - e^{- \frac{\left(1\right)}{2}}\right)\simeq 7.8694
$$

である。

外挿の往復回数は

$$
\Theta/(2\pi)=\frac{\left(20 - \frac{20}{e^{\frac{1}{2}}}\right)}{2 \pi}\simeq 1.2525
$$

である。

初期のエネルギーがしきい値以上である。
:::

:::example[鞍点付近]
微小振動の尺度を分離軌道へ使えるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
R&=10 \\s&=1 \\X&=0 \\\alpha&=\pi
\end{aligned}
$$

これらを順に代入する。

初期の規格化 Hamiltonianは

$$
H_0/U_0=\frac{\left(0\right)^{2}}{2} - \cos{\left(\left(\pi\right) \right)}=1
$$

である。

線形外挿の波振幅は

$$
U/U_0=e^{- \left(1\right)}\simeq 0.36788
$$

である。

時刻での捕捉半幅は

$$
\Delta v_{\rm tr}/\sqrt{U_0/m}=2 \sqrt{\left(e^{-1}\right)}\simeq 1.2131
$$

である。

瞬間の往復率と減衰率の比は

$$
\omega_b/g=\left(10\right) \sqrt{\left(e^{-1}\right)}\simeq 6.0653
$$

である。

微小往復の線形外挿位相は

$$
\Theta=2 \left(10\right) \left(1 - e^{- \frac{\left(1\right)}{2}}\right)\simeq 7.8694
$$

である。

外挿の往復回数は

$$
\Theta/(2\pi)=\frac{\left(20 - \frac{20}{e^{\frac{1}{2}}}\right)}{2 \pi}\simeq 1.2525
$$

である。

この点の固定場の運動は井戸中心の調和振動ではない。
:::


## イオン音波ソリトンと KdV 方程式

低周波のイオン音波では有限振幅によって圧縮部分が速く進む。一方、有限 Debye 長による分散は短波長成分を遅らせる。冷たいイオンと Boltzmann 電子を仮定し、長さを $\lambda_{De}$、時間を $\omega_{pi}^{-1}$、速度を $c_s$、電位を $T_e/e$、密度を $n_0$ で規格化する。式は $n_t+(nu)_x=0$、$u_t+u u_x=-\phi_x$、$\phi_{xx}=e^\phi-n$ である。

### 多尺度展開

$\xi=\epsilon^{1/2}(x-t)$、$\tau=\epsilon^{3/2}t$、$n=1+\epsilon n_1+\epsilon^2n_2$、$u=\epsilon u_1+\epsilon^2u_2$、$\phi=\epsilon\phi_1+\epsilon^2\phi_2$ とする。一次は $n_1=u_1=\phi_1$ を与える。次の次数では
$$
\begin{aligned}
-\partial_\xi n_2+\partial_\xi u_2&=-\partial_\tau\phi_1-\partial_\xi(\phi_1^2),\\
-\partial_\xi u_2+\partial_\xi\phi_2&=-\partial_\tau\phi_1-\phi_1\partial_\xi\phi_1,\\
n_2&=\phi_2+\phi_1^2/2-\partial_\xi^2\phi_1.
\end{aligned}
$$
初めの二式を足し、最後の式を代入して
$$
\partial_\tau\phi_1+\phi_1\partial_\xi\phi_1+\tfrac12\partial_\xi^3\phi_1=0
\label{eq-plasma-kdv}
$$
を得る。:term[Korteweg–de Vries 方程式]{en="Korteweg–de Vries equation" yomi="こるてうぇぐどふりーすほうていしき"} である。イオン温度、共鳴減衰、高次の非線形性は落とした。

### 孤立波

$z=\xi-U\tau$、$\phi_1=F(z)$ とし、$F,F',F''\to0$ を課す。一回積分で $F''=2UF-F^2$、さらに $F'$ を掛けて $(F')^2=2UF^2-(2/3)F^3$ となる。山で $F'=0$ だから振幅は $3U$、正の $U$ で減衰尾部が可能である。$F=3U\,\mathrm{sech}^2(\sqrt{U/2}\,z)$ は $(\mathrm{sech}^2 y)''=4\,\mathrm{sech}^2 y-6\,\mathrm{sech}^4 y$ により微分方程式を満たす。

:::example[幅と振幅]
$U=2$ では振幅 $6$、$\mathrm{sech}^2$ の引数の逆係数は $1$ である。物理電位の振幅には $\epsilon T_e/e$ が掛かるので $\epsilon$ の小ささが弱非線形性を支える。
:::

:::example[分散項の符号の確認]
式 :eqref[eq-plasma-kdv] を線形化し $e^{i(K\xi-\Omega\tau)}$ とすると $\Omega=-K^3/2$ である。元の長波長イオン音波でも無次元で $\omega=k/\sqrt{1+k^2}=k-k^3/2+\order{k^5}$ であり、音速とともに動く座標での先頭の周波数差に一致する。
:::

### 有限振幅の定常波との接続

同じ冷イオンと等温電子の無次元方程式に、速さ $M>0$ で静かに進む解 $\phi(x-Mt)$ を求める。遠方で $n=1,u=0,\phi=0$ とし、波の座標で積分すると $n(u-M)=-M$、$(u-M)^2/2+\phi=M^2/2$ である。遠方の符号を保つ枝では $n=(1-2\phi/M^2)^{-1/2}$ であり、$\phi=M^2/2$ で冷イオン密度が発散する。Poisson 方程式の第一積分は
$$
\frac12(\phi')^2=W_M(\phi),\qquad
W_M(\phi)=e^\phi-1+M^2(\sqrt{1-2\phi/M^2}-1).
$$
小さい正の $\phi$ に滑らかに減衰する尾部を持つには、二次係数 $(1-M^{-2})/2>0$ が必要である。$M>1$ で $W_M'=e^\phi-(1-2\phi/M^2)^{-1/2}$ を調べる。符号は $\phi+\ln(1-2\phi/M^2)/2$ の符号に等しい。この関数は厳密に凹で、原点で正の微分を持ち、$\phi\to M^2/2$ で $-\infty$ に向かう。従って $W_M$ は最初増加し、一度極大を経て減る。境界値 $e^{M^2/2}-1-M^2$ が負なら、内部に正の単純な根 $\phi_m$ が一つある。$\phi_m$ で $\phi'=0$ を選び、左右に減少する枝を積分すれば滑らかな孤立波が得られる。単純根の近傍では $\phi_m-\phi$ が距離の二乗に比例し、原点付近の二次根では遠方尾部が指数減衰する。境界値が非負なら、この冷イオン枝に同型の滑らかな正の転回点はない。

$M=1+\epsilon U$、$\phi=O(\epsilon)$ で展開すると $W_M(\phi)=\epsilon U\phi^2-\phi^3/3$ が先頭の非零項である。根の高さは $3\epsilon U$、尾部の指数は $\sqrt{2\epsilon U}$ であり、KdV 解を元の無次元座標に戻した結果と一致する。この有限振幅の第一積分の利用は、流体閉包と単一の通過枝が維持される範囲に限る。

:::example[定常波の存在の判定]
$M=11/10$ なら端点での値は $e^{121/200}-1-121/100\simeq-0.379$ である。$M>1$ で二次係数が正であり、先の単調性を合わせると正の転回点がある。$M=2$ では端点が $e^2-5>0$ なので、同じ冷イオン枝はそのような転回点を持たない。
:::

### KdV 孤立波の積分不変量

$u_t+Auu_x+Bu_{xxx}=0$、$A,B>0$ とする。進行座標 $\xi=x-Vt$ で全微分を一度積分すると $Bu''=Vu-Au^2/2$ である。さらに $2u'$ を掛けて積分し、遠方で $u,u'\to0$ とすれば $B(u')^2=Vu^2-Au^3/3$ となる。$u=U\operatorname{sech}^2(\xi/\ell)$ を代入すると $V=AU/3$、$\ell=\sqrt{4B/V}$ である。置換 $y=\tanh(\xi/\ell)$ から $\int u\dd\xi=2U\ell$ と $\int u^2\dd\xi=4U^2\ell/3$ を得る。方程式を直接積分し、また $u$ を掛けて部分積分するとこの二量は時間一定である。従って早い孤立波と遅い孤立波の幅・積分量を別々に比較できる。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $A$ | KdV の非線形係数 | 1 |
| $B$ | KdV の分散係数 | 1 |
| $U$ | 孤立波の最大振幅 | 1 |
| $X$ | 境界からの距離を幅で割った値 $\xi/\ell$ | 1 |

孤立波の速度補正を $V$ とする。基準の線形速度からの相対速度である。

$$
V=\frac{A U}{3}
$$

孤立波の幅を $\ell$ とする。振幅が増すと狭くなる。

$$
\ell=2 \sqrt{\frac{B}{V}}
$$

規格化位置での振幅を $u(X)$ とする。双曲線正割の二乗を使う。

$$
u(X)=\frac{4 U}{\left(e^{X} + e^{- X}\right)^{2}}
$$

一次積分を $I_1$ とする。保存する面積である。

$$
I_1=2 U \ell
$$

二乗積分を $I_2$ とする。一次積分とは異なる振幅依存性を持つ。

$$
I_2=\frac{4 U^{2} \ell}{3}
$$

面積で規格化した二乗積分を $I_2/I_1$ とする。形を保つ解では振幅に比例する。

$$
I_2/I_1=\frac{I_{2}}{I_{1}}
$$

:::example[小振幅の孤立波]
速度と幅を計算する。

入力値を次のように置く。

$$
\begin{aligned}
A&=1 \\B&=1 \\U&=\frac{3}{100} \\X&=0
\end{aligned}
$$

これらを順に代入する。

孤立波の速度補正は

$$
V=\frac{\left(1\right) \left(\frac{3}{100}\right)}{3}=\frac{1}{100}
$$

である。

孤立波の幅は

$$
\ell=2 \sqrt{\frac{\left(1\right)}{\left(\frac{1}{100}\right)}}=20
$$

である。

規格化位置での振幅は

$$
u(X)=\frac{4 \left(\frac{3}{100}\right)}{\left(e^{\left(0\right)} + e^{- \left(0\right)}\right)^{2}}=\frac{3}{100}
$$

である。

一次積分は

$$
I_1=2 \left(20\right) \left(\frac{3}{100}\right)=\frac{6}{5}
$$

である。

二乗積分は

$$
I_2=\frac{4 \left(20\right) \left(\frac{3}{100}\right)^{2}}{3}=\frac{3}{125}
$$

である。

面積で規格化した二乗積分は

$$
I_2/I_1=\frac{\left(\frac{3}{125}\right)}{\left(\frac{6}{5}\right)}=\frac{1}{50}
$$

である。

速度補正は小さく、幅は長波長側にある。
:::

:::example[振幅を四倍にする]
速度と二種類の積分の倍率を調べる。

入力値を次のように置く。

$$
\begin{aligned}
A&=1 \\B&=1 \\U&=\frac{3}{25} \\X&=0
\end{aligned}
$$

これらを順に代入する。

孤立波の速度補正は

$$
V=\frac{\left(1\right) \left(\frac{3}{25}\right)}{3}=\frac{1}{25}
$$

である。

孤立波の幅は

$$
\ell=2 \sqrt{\frac{\left(1\right)}{\left(\frac{1}{25}\right)}}=10
$$

である。

規格化位置での振幅は

$$
u(X)=\frac{4 \left(\frac{3}{25}\right)}{\left(e^{\left(0\right)} + e^{- \left(0\right)}\right)^{2}}=\frac{3}{25}
$$

である。

一次積分は

$$
I_1=2 \left(10\right) \left(\frac{3}{25}\right)=\frac{12}{5}
$$

である。

二乗積分は

$$
I_2=\frac{4 \left(10\right) \left(\frac{3}{25}\right)^{2}}{3}=\frac{24}{125}
$$

である。

面積で規格化した二乗積分は

$$
I_2/I_1=\frac{\left(\frac{24}{125}\right)}{\left(\frac{12}{5}\right)}=\frac{2}{25}
$$

である。

幅は半分でも二乗積分は八倍になる。
:::

:::example[強い分散]
振幅を保って分散が幅に与える効果を求める。

入力値を次のように置く。

$$
\begin{aligned}
A&=1 \\B&=4 \\U&=\frac{3}{100} \\X&=1
\end{aligned}
$$

これらを順に代入する。

孤立波の速度補正は

$$
V=\frac{\left(1\right) \left(\frac{3}{100}\right)}{3}=\frac{1}{100}
$$

である。

孤立波の幅は

$$
\ell=2 \sqrt{\frac{\left(4\right)}{\left(\frac{1}{100}\right)}}=40
$$

である。

規格化位置での振幅は

$$
u(X)=\frac{4 \left(\frac{3}{100}\right)}{\left(e^{\left(1\right)} + e^{- \left(1\right)}\right)^{2}}\simeq 0.012599
$$

である。

一次積分は

$$
I_1=2 \left(40\right) \left(\frac{3}{100}\right)=\frac{12}{5}
$$

である。

二乗積分は

$$
I_2=\frac{4 \left(40\right) \left(\frac{3}{100}\right)^{2}}{3}=\frac{6}{125}
$$

である。

面積で規格化した二乗積分は

$$
I_2/I_1=\frac{\left(\frac{6}{125}\right)}{\left(\frac{12}{5}\right)}=\frac{1}{50}
$$

である。

同じ規格化位置の振幅は変わらないが物理的な幅は二倍になる。
:::

:::example[波の裾の評価]
中心から数幅離れた位置を求める。

入力値を次のように置く。

$$
\begin{aligned}
A&=1 \\B&=1 \\U&=\frac{3}{100} \\X&=3
\end{aligned}
$$

これらを順に代入する。

孤立波の速度補正は

$$
V=\frac{\left(1\right) \left(\frac{3}{100}\right)}{3}=\frac{1}{100}
$$

である。

孤立波の幅は

$$
\ell=2 \sqrt{\frac{\left(1\right)}{\left(\frac{1}{100}\right)}}=20
$$

である。

規格化位置での振幅は

$$
u(X)=\frac{4 \left(\frac{3}{100}\right)}{\left(e^{\left(3\right)} + e^{- \left(3\right)}\right)^{2}}\simeq 0.00029598
$$

である。

一次積分は

$$
I_1=2 \left(20\right) \left(\frac{3}{100}\right)=\frac{6}{5}
$$

である。

二乗積分は

$$
I_2=\frac{4 \left(20\right) \left(\frac{3}{100}\right)^{2}}{3}=\frac{3}{125}
$$

である。

面積で規格化した二乗積分は

$$
I_2/I_1=\frac{\left(\frac{3}{125}\right)}{\left(\frac{6}{5}\right)}=\frac{1}{50}
$$

である。

広い孤立波でも規格化距離が大きければ振幅は指数的に小さい。
:::


## パラメトリック不安定性

振動する場が別の波の復元力を周期的に変えると、外部のポンプからエネルギーを受け取れる。:term[パラメトリック不安定性]{en="parametric instability" yomi="ぱらめとりっくふあんていせい"} の最小模型は $\ddot x+\omega_0^2[1+h\cos(2\omega_0t)]x=0$ である。$x=A(t)e^{i\omega_0t}+\conj A(t)e^{-i\omega_0t}$ と置き、$\ddot A$ と非共鳴な $3\omega_0$ の項を捨てると $2i\omega_0\dot A+(h\omega_0^2/2)\conj A=0$。共役の式を微分に入れると $\ddot A=(h^2\omega_0^2/16)A$ なので成長率は $\abs h\omega_0/4$ である。

三波ならポンプの $(\omega_0,\vb k_0)$ と娘波の和が一致すると速い位相が相殺される。共鳴に近く、ポンプ振幅が一定の模型で $\dot a_1=g a_2^* e^{-i\Delta t}$、$\dot a_2^*=g^* a_1 e^{i\Delta t}$ とすれば、$a_1=b_1 e^{-i\Delta t/2}$、$a_2^*=b_2 e^{i\Delta t/2}$ で定係数になる。行列の対角が $\pm i\Delta/2$、非対角が $g,g^*$ なので $\gamma^2=\abs g^2-\Delta^2/4$ である。実際の散乱では $g$ は波の偏波と非線形電流から求める。

:::example[離調のしきい値]
$\abs g=3\,\mathrm{s^{-1}}$、$\Delta=4\,\mathrm{s^{-1}}$ なら成長率は $\sqrt5\,\mathrm{s^{-1}}$ である。$\abs\Delta=6\,\mathrm{s^{-1}}$ に達するとこの指数成長率はゼロとなる。
:::

:::example[娘波の減衰を含める]
二つの振幅に等しい減衰率 $\nu_d$ を加えると、定係数行列に $-\nu_d I$ が加わる。実数の正の成長には $\sqrt{\abs g^2-\Delta^2/4}>\nu_d$ が必要である。例えば $\abs g=3,\Delta=0,\nu_d=1$ を $\mathrm{s^{-1}}$ で与えれば最大成長率は $2\,\mathrm{s^{-1}}$ である。
:::

### 娘波の減衰率が異なる場合

三波結合に正の減衰率 $d_1,d_2$ を入れると、特性方程式は
$$
(\lambda+d_1-i\Delta/2)(\lambda+d_2+i\Delta/2)=|g|^2
$$
となる。中立点 $\lambda=i\nu$ では、虚部の条件から $\nu=(d_2-d_1)\Delta/[2(d_1+d_2)]$ である。この値を実部に戻せば、臨界結合は
$$
|g|_{\rm c}^2=d_1d_2\left(1+\frac{\Delta^2}{(d_1+d_2)^2}\right)
$$
となる。$|g|=0$ では両根の実部が負である。大きな $|g|$ では一根が正の実部を持ち、実部がゼロになる点は上の一つだけなので、この値を超えた側が不安定である。例えば $d_1=1\,\mathrm{s^{-1}}$、$d_2=3\,\mathrm{s^{-1}}$、$\Delta=4\,\mathrm{s^{-1}}$ なら $|g|_{\rm c}=\sqrt6\,\mathrm{s^{-1}}$ である。等減衰の場合と異なり中立点は一般に $\lambda=0$ からずれる。

### ポンプも変化する三波結合

離調と散逸を省き、正エネルギーの三波の作用を $N_j=|a_j|^2$ と規格化する。結合定数を実数 $g$ に取れる位相を選ぶと、共鳴三波の閉じた方程式は
$$
\dot a_0=-g a_1a_2,\qquad
\dot a_1=g a_0a_2^*,\qquad
\dot a_2=g a_0a_1^*
$$
となる。このとき $\dot N_0=-\dot N_1=-\dot N_2=-2g\Re(a_0^*a_1a_2)$ である。従って $N_0+N_1$ と $N_0+N_2$ は保存される。これらの :term[Manley–Rowe 関係]{en="Manley–Rowe relations" yomi="まんりーろーかんけい"} は、$\omega_0=\omega_1+\omega_2$ と合わせてエネルギー $\sum_j\omega_jN_j$ の保存を与える。一定ポンプの指数成長は、$N_0$ の減少が小さい初期区間で近づく。例えば娘波の作用が各々 $2$ 増えれば、ポンプの作用は $2$ 減る。三波の周波数が異なるので、エネルギー変化に直した各波の値は一般に等しくない。

### 高周波の包絡が作る力

一次元で $E=E_0(x)\cos\omega_0t$ とし、場の包絡が粒子の一周期の変位に比べて緩やかに変わるとする。中心 $X$ のまわりの高速変位は $\xi=-qE_0(X)\cos\omega_0t/(m\omega_0^2)$ である。$qE_0(X+\xi)\cos\omega_0t$ のうち、$qE_0(X)\cos\omega_0t$ は平均ゼロである。残る一次の勾配項を平均すると
$$
\langle F\rangle=q E_0'(X)\langle\xi\cos\omega_0t\rangle
=-\partial_X\left(\frac{q^2E_0(X)^2}{4m\omega_0^2}\right).
$$
この :term[ポンデロモーティブ力]{en="ponderomotive force" yomi="ぽんでろもーてぃぶりょく"} は、電荷の符号によらず強い場から弱い場へ向く。高周波の包絡と低周波の密度摂動の相互作用を考える際の一つの結合機構である。磁化媒質や粒子共鳴の近くでは、高速変位の計算をその応答に合わせて行う。 :cite[bellan-plasma]

### 娘波の減衰が等しい三波結合

ポンプを一定とし娘波振幅 $a,b$ の位相を回して、行列 $\begin{pmatrix}-d+i\Delta/2&G\\G&-d-i\Delta/2\end{pmatrix}$ を得る。ここで $b$ には一方の娘波の複素共役を選んでいる。行列式から $(\lambda+d)^2=G^2-\Delta^2/4$ となる。成長の条件は $G^2>d^2+\Delta^2/4$ である。離調の時間因子は除去しても対角項の周波数差として残る。根が複素の場合、両娘波は同じ減衰 $d$ を受ける。実の正根が得られてもポンプの減衰と娘波への移行を無視した初期成長の結果である。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $G$ | 結合率 | $\mathrm{s^{-1}}$ |
| $d$ | 娘波に共通の振幅減衰率 | $\mathrm{s^{-1}}$ |
| $\Delta$ | 周波数の離調 | $\mathrm{s^{-1}}$ |

減衰を除く二乗成長率を $Q$ とする。負の場合は固有値に虚部がある。

$$
Q=G^{2} - \frac{\Delta^{2}}{4}
$$

成長のしきい値を $G_c$ とする。最大実部がゼロになる値である。

$$
G_c=\sqrt{\frac{\Delta^{2}}{4} + d^{2}}
$$

しきい値との比を $G/G_c$ とする。一より大きいと活性である。

$$
G/G_c=\frac{G}{G_{c}}
$$

最大の実成長率を $\Re\lambda_+$ とする。負の二乗根は実部に寄与しない。

$$
\Re\lambda_+=- d + \sqrt{\frac{Q}{2} + \frac{\left|{Q}\right|}{2}}
$$

回転座標での振動率を $\abs{\Im\lambda_+}$ とする。離調が大きいと有限になる。

$$
\abs{\Im\lambda_+}=\sqrt{- \frac{Q}{2} + \frac{\left|{Q}\right|}{2}}
$$

固有値の積を $\lambda_+\lambda_-$ とする。境界でゼロになる。

$$
\lambda_+\lambda_-=- Q + d^{2}
$$

:::example[共鳴上の成長]
離調を消した場合の正根を求める。

入力値を次のように置く。

$$
\begin{aligned}
G&=3\,\mathrm{s^{-1}} \\d&=1\,\mathrm{s^{-1}} \\\Delta&=0\,\mathrm{s^{-1}}
\end{aligned}
$$

これらを順に代入する。

減衰を除く二乗成長率は

$$
Q=- \frac{\left(0\right)^{2}}{4} + \left(3\right)^{2}=9\,\mathrm{s^{-2}}
$$

である。

成長のしきい値は

$$
G_c=\sqrt{\frac{\left(0\right)^{2}}{4} + \left(1\right)^{2}}=1\,\mathrm{s^{-1}}
$$

である。

しきい値との比は

$$
G/G_c=\frac{\left(3\right)}{\left(1\right)}=3
$$

である。

最大の実成長率は

$$
\Re\lambda_+=- \left(1\right) + \sqrt{\frac{\left(9\right)}{2} + \frac{\left|{\left(9\right)}\right|}{2}}=2\,\mathrm{s^{-1}}
$$

である。

回転座標での振動率は

$$
\abs{\Im\lambda_+}=\sqrt{- \frac{\left(9\right)}{2} + \frac{\left|{\left(9\right)}\right|}{2}}=0\,\mathrm{s^{-1}}
$$

である。

固有値の積は

$$
\lambda_+\lambda_-=\left(1\right)^{2} - \left(9\right)=-8\,\mathrm{s^{-2}}
$$

である。

共通減衰が結合率より小さい。
:::

:::example[離調による抑制]
同じポンプで離調だけを増やす。

入力値を次のように置く。

$$
\begin{aligned}
G&=3\,\mathrm{s^{-1}} \\d&=1\,\mathrm{s^{-1}} \\\Delta&=8\,\mathrm{s^{-1}}
\end{aligned}
$$

これらを順に代入する。

減衰を除く二乗成長率は

$$
Q=\left(3\right)^{2} - \frac{\left(8\right)^{2}}{4}=-7\,\mathrm{s^{-2}}
$$

である。

成長のしきい値は

$$
G_c=\sqrt{\left(1\right)^{2} + \frac{\left(8\right)^{2}}{4}}\simeq 4.1231\,\mathrm{s^{-1}}
$$

である。

しきい値との比は

$$
G/G_c=\frac{\left(3\right)}{\left(\sqrt{17}\right)}\simeq 0.72761
$$

である。

最大の実成長率は

$$
\Re\lambda_+=- \left(1\right) + \sqrt{\frac{\left(-7\right)}{2} + \frac{\left|{\left(-7\right)}\right|}{2}}=-1\,\mathrm{s^{-1}}
$$

である。

回転座標での振動率は

$$
\abs{\Im\lambda_+}=\sqrt{- \frac{\left(-7\right)}{2} + \frac{\left|{\left(-7\right)}\right|}{2}}\simeq 2.6458\,\mathrm{s^{-1}}
$$

である。

固有値の積は

$$
\lambda_+\lambda_-=- \left(-7\right) + \left(1\right)^{2}=8\,\mathrm{s^{-2}}
$$

である。

実部は負で固有値には振動成分が出る。
:::

:::example[離調しても成長]
結合を増して回復するか調べる。

入力値を次のように置く。

$$
\begin{aligned}
G&=5\,\mathrm{s^{-1}} \\d&=1\,\mathrm{s^{-1}} \\\Delta&=6\,\mathrm{s^{-1}}
\end{aligned}
$$

これらを順に代入する。

減衰を除く二乗成長率は

$$
Q=\left(5\right)^{2} - \frac{\left(6\right)^{2}}{4}=16\,\mathrm{s^{-2}}
$$

である。

成長のしきい値は

$$
G_c=\sqrt{\left(1\right)^{2} + \frac{\left(6\right)^{2}}{4}}\simeq 3.1623\,\mathrm{s^{-1}}
$$

である。

しきい値との比は

$$
G/G_c=\frac{\left(5\right)}{\left(\sqrt{10}\right)}\simeq 1.5811
$$

である。

最大の実成長率は

$$
\Re\lambda_+=- \left(1\right) + \sqrt{\frac{\left(16\right)}{2} + \frac{\left|{\left(16\right)}\right|}{2}}=3\,\mathrm{s^{-1}}
$$

である。

回転座標での振動率は

$$
\abs{\Im\lambda_+}=\sqrt{- \frac{\left(16\right)}{2} + \frac{\left|{\left(16\right)}\right|}{2}}=0\,\mathrm{s^{-1}}
$$

である。

固有値の積は

$$
\lambda_+\lambda_-=- \left(16\right) + \left(1\right)^{2}=-15\,\mathrm{s^{-2}}
$$

である。

離調と減衰の両方を越えて実の正根が現れる。
:::

:::example[共鳴上で中立]
結合と減衰が同じ境界を確認する。

入力値を次のように置く。

$$
\begin{aligned}
G&=1\,\mathrm{s^{-1}} \\d&=1\,\mathrm{s^{-1}} \\\Delta&=0\,\mathrm{s^{-1}}
\end{aligned}
$$

これらを順に代入する。

減衰を除く二乗成長率は

$$
Q=- \frac{\left(0\right)^{2}}{4} + \left(1\right)^{2}=1\,\mathrm{s^{-2}}
$$

である。

成長のしきい値は

$$
G_c=\sqrt{\frac{\left(0\right)^{2}}{4} + \left(1\right)^{2}}=1\,\mathrm{s^{-1}}
$$

である。

しきい値との比は

$$
G/G_c=1
$$

である。

最大の実成長率は

$$
\Re\lambda_+=- \left(1\right) + \sqrt{\frac{\left(1\right)}{2} + \frac{\left|{\left(1\right)}\right|}{2}}=0\,\mathrm{s^{-1}}
$$

である。

回転座標での振動率は

$$
\abs{\Im\lambda_+}=\sqrt{- \frac{\left(1\right)}{2} + \frac{\left|{\left(1\right)}\right|}{2}}=0\,\mathrm{s^{-1}}
$$

である。

固有値の積は

$$
\lambda_+\lambda_-=\left(1\right)^{2} - \left(1\right)=0\,\mathrm{s^{-2}}
$$

である。

一方の根はゼロ、他方は減衰する。
:::


## 無衝突衝撃波

平均自由行程が層の厚さより長くても、電磁場で粒子が反射・散乱されると衝撃波ができる。これを :term[無衝突衝撃波]{en="collisionless shock" yomi="むしょうとつしょうげきは"} と呼ぶ。層の細部の機構と、上下流の保存則は分けて計算できる。[圧縮性流体と衝撃波](/physics/continuum-mechanics/compressible-flow/) と同じ保存形の議論を用いる。

### 垂直 MHD 衝撃波の収支

衝撃波静止系で法線速度 $u$、接線磁場 $B_t$、上流を 1、下流を 2 とする。定常の質量、誘導、法線運動量、総エネルギー式を積分して
$$
[\rho u]=[uB_t]=[\rho u^2+p+B_t^2/(2\mu_0)]
=\left[u\left(\rho u^2/2+\gamma_{\rm ad}p/(\gamma_{\rm ad}-1)+B_t^2/\mu_0\right)\right]=0.
$$
括弧 $[a]$ は $a_2-a_1$ である。最初の二式から圧縮比 $r=\rho_2/\rho_1$ に対し $u_2=u_1/r$、$B_{t2}=rB_{t1}$ となる。磁気圧は圧縮を妨げる。上流磁場と熱圧を極限としてゼロにすると、運動量式は $p_2=\rho_1u_1^2(1-1/r)$。エネルギー式に入れて $(r-1)[r+1-2\gamma_{\rm ad}/(\gamma_{\rm ad}-1)]=0$ から、非自明な根は $r=(\gamma_{\rm ad}+1)/(\gamma_{\rm ad}-1)$ である。

保存式は層の幅や粒子反射を決めない。電位差、磁気ミラー、分布の不安定性が反射粒子を作り、散乱と位相混合で粗視化した分布を変化させる。非熱成分や軌道の非局所性が無視できない場合、等方 MHD で上下流を閉じる仮定も誤差を持つ。

:::example[強い非磁化極限]
$\gamma_{\rm ad}=5/3$ なら $r=4$ である。磁場を無視できない垂直衝撃波にこの値をそのまま代入すると、磁気圧の仕事が抜ける。
:::

:::example[運動量流束の分担]
垂直衝撃波で $r=2$ なら下流の動圧は上流の半分になる。磁気圧は四倍になるので、運動量式から $p_2-p_1=\rho_1u_1^2/2-3B_{t1}^2/(2\mu_0)$ である。この式だけでは $r$ を選べず、エネルギー流束の条件も必要である。
:::

### 垂直衝撃波の下流圧力を保存則で検査する

衝撃波静止系で $\rho u$ と $uB$ が一定とし、法線磁場をゼロとする。圧縮比 $r=\rho_2/\rho_1$ なら $u_2=u_1/r$、$B_2=rB_1$ である。規格化上流熱圧 $a=p_1/(\rho_1u_1^2)$、磁気圧 $b=B_1^2/(2\mu_0\rho_1u_1^2)$ に対して、法線運動量から $p_2/(\rho_1u_1^2)=a+1-1/r+b(1-r^2)$ となる。検査にはさらにエネルギーが必要で、単位質量当たりの流束は $u^2/2+\gamma p/[(\gamma-1)\rho]+B^2/(\mu_0\rho)$ である。保存則だけでは層厚と粒子加速は分からない。非放射の垂直 MHD 跳躍としては下流が正圧で、エネルギー残差がゼロかつエントロピーが増えなければならない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $r$ | 仮定する圧縮比 | 1 |
| $a$ | 上流の規格化熱圧 | 1 |
| $b$ | 上流の規格化磁気圧 | 1 |

下流速度の比を $u_2/u_1$ とする。質量流束を保つ。

$$
u_2/u_1=\frac{1}{r}
$$

下流の磁場は $B_2=M_B B_1$ と書ける。磁場がゼロの場合にもこの式を使える。

$$
M_B=r
$$

下流熱圧を $p_2/(\rho_1u_1^2)$ とする。運動量収支で残る圧力である。

$$
p_2/(\rho_1u_1^2)=a + b \left(1 - r^{2}\right) - (u_2/u_1) + 1
$$

上流の質量流束で割ったエネルギー流束を $E_1/u_1^2$ とする。$\gamma=5/3$ とする。

$$
E_1/u_1^2=\frac{5 a}{2} + 2 b + \frac{1}{2}
$$

下流の質量流束で割ったエネルギー流束を $E_2/u_1^2$ とする。横磁場の輸送に係数二が付く。

$$
E_2/u_1^2=2 b r + \frac{5 [p_2/(\rho_1u_1^2)]}{2 r} + \frac{(u_2/u_1)^{2}}{2}
$$

エネルギー保存の残差を $(E_2-E_1)/u_1^2$ とする。保存する跳躍に必要な条件はゼロである。

$$
(E_2-E_1)/u_1^2=\frac{E_2}{u_1^2}-\frac{E_1}{u_1^2}
$$

:::example[強い非磁化衝撃波]
冷たい上流の圧縮比四を代入する。

入力値を次のように置く。

$$
\begin{aligned}
r&=4 \\a&=0 \\b&=0
\end{aligned}
$$

これらを順に代入する。

下流速度の比は

$$
u_2/u_1=\frac{1}{\left(4\right)}=\frac{1}{4}
$$

である。

下流磁場の倍率は

$$
M_B=\left(4\right)=4
$$

である。

下流熱圧は

$$
p_2/(\rho_1u_1^2)=\left(0\right) \left(1 - \left(4\right)^{2}\right) + \left(0\right) - \left(\frac{1}{4}\right) + 1=\frac{3}{4}
$$

である。

上流の比エネルギー流束は

$$
E_1/u_1^2=\frac{9 \left(0\right)}{2} + \frac{1}{2}=\frac{1}{2}
$$

である。

下流の比エネルギー流束は

$$
E_2/u_1^2=2 \left(0\right) \left(4\right) + \frac{\left(\frac{1}{4}\right)^{2}}{2} + \frac{5 \left(\frac{3}{4}\right)}{2 \left(4\right)}=\frac{1}{2}
$$

である。

エネルギー保存の残差は

$$
(E_2-E_1)/u_1^2=0
$$

である。

運動量とエネルギーの両方を満たす強衝撃波の極限である。
:::

:::example[圧縮比の不適切な外挿]
質量と運動量だけで生じる誤りを確認する。

入力値を次のように置く。

$$
\begin{aligned}
r&=5 \\a&=0 \\b&=0
\end{aligned}
$$

これらを順に代入する。

下流速度の比は

$$
u_2/u_1=\frac{1}{\left(5\right)}=\frac{1}{5}
$$

である。

下流磁場の倍率は

$$
M_B=\left(5\right)=5
$$

である。

下流熱圧は

$$
p_2/(\rho_1u_1^2)=\left(0\right) \left(1 - \left(5\right)^{2}\right) + \left(0\right) - \left(\frac{1}{5}\right) + 1=\frac{4}{5}
$$

である。

上流の比エネルギー流束は

$$
E_1/u_1^2=\frac{9 \left(0\right)}{2} + \frac{1}{2}=\frac{1}{2}
$$

である。

下流の比エネルギー流束は

$$
E_2/u_1^2=2 \left(0\right) \left(5\right) + \frac{\left(\frac{1}{5}\right)^{2}}{2} + \frac{5 \left(\frac{4}{5}\right)}{2 \left(5\right)}=\frac{21}{50}
$$

である。

エネルギー保存の残差は

$$
(E_2-E_1)/u_1^2=- \left(\frac{1}{2}\right) + \left(\frac{21}{50}\right)=- \frac{2}{25}
$$

である。

正圧であってもエネルギー残差が消えない。
:::

:::example[圧縮なし]
保存式に現れる自明根を求める。

入力値を次のように置く。

$$
\begin{aligned}
r&=1 \\a&=\frac{1}{10} \\b&=\frac{1}{10}
\end{aligned}
$$

これらを順に代入する。

下流速度の比は

$$
u_2/u_1=\frac{1}{\left(1\right)}=1
$$

である。

下流磁場の倍率は

$$
M_B=\left(1\right)=1
$$

である。

下流熱圧は

$$
p_2/(\rho_1u_1^2)=- \left(1\right) + \left(\frac{1}{10}\right) \left(1 - \left(1\right)^{2}\right) + \left(\frac{1}{10}\right) + 1=\frac{1}{10}
$$

である。

上流の比エネルギー流束は

$$
E_1/u_1^2=\frac{9 \left(\frac{1}{10}\right)}{2} + \frac{1}{2}=\frac{19}{20}
$$

である。

下流の比エネルギー流束は

$$
E_2/u_1^2=\frac{\left(1\right)^{2}}{2} + 2 \left(1\right) \left(\frac{1}{10}\right) + \frac{5 \left(\frac{1}{10}\right)}{2 \left(1\right)}=\frac{19}{20}
$$

である。

エネルギー保存の残差は

$$
(E_2-E_1)/u_1^2=0
$$

である。

全て連続なので衝撃波とは区別する。
:::

:::example[磁気圧が大きすぎる圧縮]
下流圧力の正値を確かめる。

入力値を次のように置く。

$$
\begin{aligned}
r&=4 \\a&=\frac{1}{10} \\b&=\frac{1}{10}
\end{aligned}
$$

これらを順に代入する。

下流速度の比は

$$
u_2/u_1=\frac{1}{\left(4\right)}=\frac{1}{4}
$$

である。

下流磁場の倍率は

$$
M_B=\left(4\right)=4
$$

である。

下流熱圧は

$$
p_2/(\rho_1u_1^2)=\left(\frac{1}{10}\right) \left(1 - \left(4\right)^{2}\right) + \left(\frac{1}{10}\right) - \left(\frac{1}{4}\right) + 1=- \frac{13}{20}
$$

である。

上流の比エネルギー流束は

$$
E_1/u_1^2=\frac{9 \left(\frac{1}{10}\right)}{2} + \frac{1}{2}=\frac{19}{20}
$$

である。

下流の比エネルギー流束は

$$
E_2/u_1^2=\frac{5 \left(- \frac{13}{20}\right)}{2 \left(4\right)} + 2 \left(4\right) \left(\frac{1}{10}\right) + \frac{\left(\frac{1}{4}\right)^{2}}{2}=\frac{17}{40}
$$

である。

エネルギー保存の残差は

$$
(E_2-E_1)/u_1^2=\left(\frac{17}{40}\right) - \left(\frac{19}{20}\right)=- \frac{21}{40}
$$

である。

圧縮された磁場に運動量が使われ、指定の圧縮比は負圧を与える。
:::


### 垂直速波衝撃波の圧縮比を解く

前の保存式で $\gamma=5/3$ とし、上流熱圧比を $a$、磁気圧比を $b$ とおく。下流の圧力を運動量の式で消去してエネルギー式に代入すると
$$
2r^2(E_2-E_1)/u_1^2
=-(r-1)\bigl[b r^2+(1+5a+5b)r-4\bigr].
\label{eq-plasma-perp-shock-jump}
$$
$r=1$ は前後が同じ根である。$a,b\ge0$ なら角括弧の多項式は正の半直線で単調増加し、原点で負である。唯一の正根が一より大きい条件は $5a+6b<3$ である。上流の垂直速波速度は $c_{f1}^2/u_1^2=5a/3+2b$ だから、同じ条件は上流が法線方向に超速波であることを表す。$b\to0$ で有限の根を保つため二次式を有理化すると $8/[1+5a+5b+\sqrt{(1+5a+5b)^2+16b}]$ となる。熱圧は先の運動量式で求める。十分条件を明示するため、$a>0$ として二次式から $5a=4/r-1-b(r+5)$ を入れると
$$
\begin{aligned}
5r\left(\frac{p_2}{\rho_1u_1^2}-ar^{5/3}\right)
&=4r-1-(4-r)r^{5/3}\\
&\quad+br^2\big[(r+5)r^{2/3}-(1+5r)\big].
\end{aligned}
$$
第一行と角括弧の関数は $r=1$ で値と一階微分がともにゼロである。二階微分は順に $(40/9)(r-1)r^{-1/3}$、$(10/9)(r-1)r^{-4/3}$ となり、$r>1$ で正である。従って圧縮根は $p_2>p_1r^{5/3}>0$ を満たし、気体のエントロピーは増える。:eqref[eq-plasma-perp-shock-jump] は巨視的な跳躍だけを与え、無衝突層での散逸の機構を指定しない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $a$ | 上流熱圧 $p_1/(\rho_1u_1^2)$ | 1 |
| $b$ | 上流磁気圧 $B_1^2/(2\mu_0\rho_1u_1^2)$ | 1 |

上流速波の二乗と流速の比を $c_{f1}^2/u_1^2$ とする。圧力と磁場がともに復元力を与える。

$$
c_{f1}^2/u_1^2=\frac{5 a}{3} + 2 b
$$

上流の速波 Mach 数を $M_{f1}$ とする。一より大きい圧縮枝を考える。

$$
M_{f1}=\frac{1}{\sqrt{c_{f1}^2/u_1^2}}
$$

圧縮比の二次式の一次係数を $C$ とする。上流の状態だけで決まる。

$$
C=5 a + 5 b + 1
$$

非自明な圧縮比を $r$ とする。$b=0$ でも有限となる形を使う。

$$
r=\frac{8}{C + \sqrt{C^{2} + 16 b}}
$$

下流の熱圧を $p_2/(\rho_1u_1^2)$ とする。磁場圧縮分を全圧から引く。

$$
p_2/(\rho_1u_1^2)=a + b \left(1 - r^{2}\right) + 1 - \frac{1}{r}
$$

理想気体の規格化エントロピー差を $\Delta s/c_v$ とする。正なら不可逆の圧縮と両立する。

$$
\Delta s/c_v=\ln{\left(\frac{[p_2/(\rho_1u_1^2)]}{a r^{\frac{5}{3}}} \right)}
$$

:::example[弱磁場の超速波流]
非自明な圧縮比と下流圧を求める。

入力値を次のように置く。

$$
\begin{aligned}
a&=\frac{1}{10} \\b&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

上流速波の二乗と流速の比は

$$
c_{f1}^2/u_1^2=2 \left(\frac{1}{100}\right) + \frac{5 \left(\frac{1}{10}\right)}{3}=\frac{14}{75}
$$

である。

上流の速波 Mach 数は

$$
M_{f1}=\frac{1}{\sqrt{\left(\frac{14}{75}\right)}}\simeq 2.3146
$$

である。

圧縮比の二次式の一次係数は

$$
C=5 \left(\frac{1}{100}\right) + 5 \left(\frac{1}{10}\right) + 1=\frac{31}{20}
$$

である。

非自明な圧縮比は

$$
r=\frac{8}{\left(\frac{31}{20}\right) + \sqrt{16 \left(\frac{1}{100}\right) + \left(\frac{31}{20}\right)^{2}}}\simeq 2.5391
$$

である。

下流の熱圧は

$$
p_2/(\rho_1u_1^2)=\left(\frac{1}{100}\right) \left(1 - \left(\frac{8}{\frac{31}{20} + \frac{\sqrt{41}}{4}}\right)^{2}\right) + \left(\frac{1}{10}\right) + 1 - \frac{1}{\left(\frac{8}{\frac{31}{20} + \frac{\sqrt{41}}{4}}\right)}\simeq 0.65168
$$

である。

理想気体の規格化エントロピー差は

$$
\Delta s/c_v=\ln{\left(\frac{\left(- \frac{\sqrt{41}}{32} - \frac{16}{25 \left(\frac{31}{20} + \frac{\sqrt{41}}{4}\right)^{2}} + \frac{733}{800}\right)}{\left(\frac{1}{10}\right) \left(\frac{8}{\frac{31}{20} + \frac{\sqrt{41}}{4}}\right)^{\frac{5}{3}}} \right)}\simeq 0.32141
$$

である。

圧縮枝で熱的なエントロピーが増える。
:::

:::example[強い横磁場]
同じ上流熱圧で圧縮比の変化を調べる。

入力値を次のように置く。

$$
\begin{aligned}
a&=\frac{1}{10} \\b&=\frac{1}{5}
\end{aligned}
$$

これらを順に代入する。

上流速波の二乗と流速の比は

$$
c_{f1}^2/u_1^2=\frac{5 \left(\frac{1}{10}\right)}{3} + 2 \left(\frac{1}{5}\right)=\frac{17}{30}
$$

である。

上流の速波 Mach 数は

$$
M_{f1}=\frac{1}{\sqrt{\left(\frac{17}{30}\right)}}\simeq 1.3284
$$

である。

圧縮比の二次式の一次係数は

$$
C=5 \left(\frac{1}{10}\right) + 5 \left(\frac{1}{5}\right) + 1=\frac{5}{2}
$$

である。

非自明な圧縮比は

$$
r=\frac{8}{\left(\frac{5}{2}\right) + \sqrt{16 \left(\frac{1}{5}\right) + \left(\frac{5}{2}\right)^{2}}}\simeq 1.4352
$$

である。

下流の熱圧は

$$
p_2/(\rho_1u_1^2)=\left(\frac{1}{10}\right) + \left(\frac{1}{5}\right) \left(1 - \left(\frac{8}{\frac{5}{2} + \frac{3 \sqrt{105}}{10}}\right)^{2}\right) + 1 - \frac{1}{\left(\frac{8}{\frac{5}{2} + \frac{3 \sqrt{105}}{10}}\right)}\simeq 0.19127
$$

である。

理想気体の規格化エントロピー差は

$$
\Delta s/c_v=\ln{\left(\frac{\left(- \frac{64}{5 \left(\frac{5}{2} + \frac{3 \sqrt{105}}{10}\right)^{2}} - \frac{3 \sqrt{105}}{80} + \frac{79}{80}\right)}{\left(\frac{1}{10}\right) \left(\frac{8}{\frac{5}{2} + \frac{3 \sqrt{105}}{10}}\right)^{\frac{5}{3}}} \right)}\simeq 0.046338
$$

である。

上流の磁場が強いほど速波 Mach 数が下がり、圧縮率が小さくなる。
:::

:::example[Mach 数が一]
しきい値で別の圧縮根があるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
a&=\frac{3}{10} \\b&=\frac{1}{4}
\end{aligned}
$$

これらを順に代入する。

上流速波の二乗と流速の比は

$$
c_{f1}^2/u_1^2=2 \left(\frac{1}{4}\right) + \frac{5 \left(\frac{3}{10}\right)}{3}=1
$$

である。

上流の速波 Mach 数は

$$
M_{f1}=\frac{1}{\sqrt{\left(1\right)}}=1
$$

である。

圧縮比の二次式の一次係数は

$$
C=5 \left(\frac{1}{4}\right) + 5 \left(\frac{3}{10}\right) + 1=\frac{15}{4}
$$

である。

非自明な圧縮比は

$$
r=\frac{8}{\left(\frac{15}{4}\right) + \sqrt{\left(\frac{15}{4}\right)^{2} + 16 \left(\frac{1}{4}\right)}}=1
$$

である。

下流の熱圧は

$$
p_2/(\rho_1u_1^2)=\left(\frac{1}{4}\right) \left(1 - \left(1\right)^{2}\right) + \left(\frac{3}{10}\right) + 1 - \frac{1}{\left(1\right)}=\frac{3}{10}
$$

である。

理想気体の規格化エントロピー差は

$$
\Delta s/c_v=\ln{\left(\frac{1}{\left(1\right)^{\frac{5}{3}}} \right)}=0
$$

である。

非自明根は一に合流し、エントロピー差はゼロになる。
:::

:::example[磁場のない極限]
$b=0$ で有理化した式を確認する。

入力値を次のように置く。

$$
\begin{aligned}
a&=\frac{1}{10} \\b&=0
\end{aligned}
$$

これらを順に代入する。

上流速波の二乗と流速の比は

$$
c_{f1}^2/u_1^2=2 \left(0\right) + \frac{5 \left(\frac{1}{10}\right)}{3}=\frac{1}{6}
$$

である。

上流の速波 Mach 数は

$$
M_{f1}=\frac{1}{\sqrt{\left(\frac{1}{6}\right)}}\simeq 2.4495
$$

である。

圧縮比の二次式の一次係数は

$$
C=5 \left(0\right) + 5 \left(\frac{1}{10}\right) + 1=\frac{3}{2}
$$

である。

非自明な圧縮比は

$$
r=\frac{8}{\left(\frac{3}{2}\right) + \sqrt{16 \left(0\right) + \left(\frac{3}{2}\right)^{2}}}=\frac{8}{3}
$$

である。

下流の熱圧は

$$
p_2/(\rho_1u_1^2)=\left(0\right) \left(1 - \left(\frac{8}{3}\right)^{2}\right) + \left(\frac{1}{10}\right) + 1 - \frac{1}{\left(\frac{8}{3}\right)}=\frac{29}{40}
$$

である。

理想気体の規格化エントロピー差は

$$
\Delta s/c_v=\ln{\left(\frac{\left(\frac{29}{40}\right)}{\left(\frac{1}{10}\right) \left(\frac{8}{3}\right)^{\frac{5}{3}}} \right)}\simeq 0.34629
$$

である。

$4/(1+5a)$ となり、冷たい強衝撃波の極限で四に近づく。
:::


## プラズマ乱流の概観

単独の波より多くのモードが同時に成長すると、波数の間でエネルギーをやり取りする。例として不可圧縮 MHD で $\vb b=\vb B/\sqrt{\mu_0\rho}$、$\vb z^\pm=\vb u\pm\vb b$ とする。理想運動式と誘導式を足し引きすると
$$
\partial_t\vb z^\pm+(\vb z^\mp\cdot\nabla)\vb z^\pm=-\nabla\Pi,\qquad
\nabla\cdot\vb z^\pm=0
$$
である。$\Pi$ は全圧を質量密度で割ったものである。周期境界で $\vb z^\pm$ と内積を取り体積積分すると、移流項は $\nabla\cdot(\vb z^\mp|\vb z^\pm|^2/2)$、圧力項も発散となるので $\int|\vb z^\pm|^2\dd V$ はそれぞれ保存される。従って非線形項は総量を消さず、尺度の間を移す。背景一様磁場を分離すれば、非線形相互作用は反対向きに伝播する二群の揺らぎの積になる。

:term[臨界釣り合い]{en="critical balance" yomi="りんかいつりあい"} の見積もりでは $k_\parallel v_A\sim k_\perp\delta u_\ell$ と置く。局所的で一定のエネルギー流束 $\varepsilon\sim\delta u_\ell^3/\ell_\perp$ も仮定すれば $\delta u_\ell\propto\ell_\perp^{1/3}$、$k_\parallel\propto k_\perp^{2/3}$ となる。この推定は強い乱流と局所的な転送を仮定する。運動論的なスケールでは電子・イオンの加熱へ至る各過程を別に追う必要がある。乱流の非局所的相互作用とジャイロ運動論の保存系はさらに場の統計を要するため、本章では上の保存則と尺度評価までを扱う。 :cite[bellan-plasma]

:::example[片側だけの波]
$\vb z^-=0$ で背景磁場を分離した系では、$\vb z^+$ はそのまま線形移流できる。非線形項に必要な反対向きの場がないためである。
:::

:::example[尺度比を変える]
局所カスケード $\delta u_\ell\propto\ell_\perp^{1/3}$ が成り立つと仮定する。$\ell_\perp$ を八分の一にすると速度差は半分、時間 $\ell_\perp/\delta u_\ell$ は四分の一になる。臨界釣り合いを用いれば平行の相関長も四分の一になるので、垂直ほど急には細かくならない。
:::

### 弱い Alfvén 波のカスケード

反対向きの波束との一回の遭遇を $\tau_A=\ell_\parallel/v_A$、非線形時間を $\tau_{\rm nl}=\ell_\perp/z_\ell$ とし、$\tau_A\ll\tau_{\rm nl}$ とする。一回の振幅変化の比は $\tau_A/\tau_{\rm nl}$ である。遭遇ごとの位相が独立と仮定すれば、$N$ 回後に比の二乗平均は $N(\tau_A/\tau_{\rm nl})^2$ であり、変化が一程度になる時刻は $\tau_{\rm cas}\sim\tau_{\rm nl}^2/\tau_A$ である。

伝播の向きが反対の波が Alfvén の分散式 $\omega^\pm=\mp k_\parallel v_A$ で三波共鳴するとき、波数と周波数の和を課すと反対向きの一成分の $k_\parallel$ がゼロになる。例えば $+$ の波が $+$ と $-$ の波の和でできる場合、$k_\parallel=k_{\parallel 1}+k_{\parallel 2}$ と $-k_\parallel=-k_{\parallel 1}+k_{\parallel 2}$ から $k_{\parallel 2}=0$ である。従って先頭の共鳴では伝播成分の平行波数が変わらない。そこで $\tau_A$ を一定とし、横への一定のエネルギー流束 $\varepsilon\sim z_\ell^2/\tau_{\rm cas}$ を仮定すると、$z_\ell\propto\ell_\perp^{1/2}$、$E(k_\perp)\propto k_\perp^{-2}$ となる。横尺度を四分の一にすればこの近似で $z_\ell$ は半分、$\tau_{\rm nl}$ も半分になる。従って小さい横尺度では弱い相互作用の仮定が破れうる。

### 局所カスケード模型の臨界釣り合い

磁場に垂直な尺度 $\ell_\perp$ で反対向き Elsässer 波の振幅が等しく $z_\ell$ であるとする。局所相互作用の時間を $\ell_\perp/z_\ell$、線形時間を $\ell_\parallel/v_A$ とし、両者が等しいと仮定する。エネルギー流束 $\epsilon\sim z_\ell^3/\ell_\perp$ が尺度に依存しなければ $z_\ell\propto\ell_\perp^{1/3}$、$\ell_\parallel\propto\ell_\perp^{2/3}$ である。一次元の横スペクトルで $k_\perp E(k_\perp)\sim z_\ell^2$ と対応させれば $E(k_\perp)\propto k_\perp^{-5/3}$ となる。これはバランスした強い乱流の仮定からの帰結であり、減衰、動的整列、間欠性がある場合に指数を固定する証明ではない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $\ell_\perp$ | 横尺度 | $\mathrm m$ |
| $\epsilon$ | 単位質量のエネルギー伝達率 | $\mathrm{m^2\,s^{-3}}$ |
| $v_A$ | 背景の Alfvén 速度 | $\mathrm{m\,s^{-1}}$ |

波の振幅の尺度を $z_\ell$ とする。定数流束を仮定する。

$$
z_\ell=\sqrt[3]{\ell_{\perp} \epsilon}
$$

非線形時間を $t_{\rm nl}$ とする。反対向き波の振幅で決まる。

$$
t_{\rm nl}=\frac{\ell_{\perp}}{z_{\ell}}
$$

平行尺度を $\ell_\parallel$ とする。伝播時間を非線形時間と等しくする。

$$
\ell_\parallel=t_{\rm nl} v_{A}
$$

異方性の長さ比を $\ell_\parallel/\ell_\perp$ とする。横尺が短いほど大きくなる。

$$
\ell_\parallel/\ell_\perp=\frac{\ell_{\parallel}}{\ell_{\perp}}
$$

横波数の尺度を $k_\perp$ とする。$2\pi$ の定数を吸収する。

$$
k_\perp=\frac{1}{\ell_{\perp}}
$$

横波数スペクトルの尺度を $E(k_\perp)$ とする。対数波数幅のエネルギーを波数で割る。

$$
E(k_\perp)=\frac{z_{\ell}^{2}}{k_{\perp}}
$$

:::example[基準の乱流尺度]
一つの尺度で釣り合う時間を求める。

入力値を次のように置く。

$$
\begin{aligned}
\ell_\perp&=1\,\mathrm m \\\epsilon&=1\,\mathrm{m^2\,s^{-3}} \\v_A&=10\,\mathrm{m\,s^{-1}}
\end{aligned}
$$

これらを順に代入する。

波の振幅の尺度は

$$
z_\ell=\sqrt[3]{\left(1\right)^{2}}=1\,\mathrm{m\,s^{-1}}
$$

である。

非線形時間は

$$
t_{\rm nl}=1\,\mathrm s
$$

である。

平行尺度は

$$
\ell_\parallel=\left(10\right) \left(1\right)=10\,\mathrm m
$$

である。

異方性の長さ比は

$$
\ell_\parallel/\ell_\perp=\frac{\left(10\right)}{\left(1\right)}=10
$$

である。

横波数の尺度は

$$
k_\perp=\frac{1}{\left(1\right)}=1\,\mathrm{m^{-1}}
$$

である。

横波数スペクトルの尺度は

$$
E(k_\perp)=\left(1\right)=1\,\mathrm{m^3\,s^{-2}}
$$

である。

背景磁場の方向に長い渦となる。
:::

:::example[小さな横尺度]
さらに小さい尺度での平行長と異方性を計算する。

入力値を次のように置く。

$$
\begin{aligned}
\ell_\perp&=\frac{1}{8}\,\mathrm m \\\epsilon&=1\,\mathrm{m^2\,s^{-3}} \\v_A&=10\,\mathrm{m\,s^{-1}}
\end{aligned}
$$

これらを順に代入する。

波の振幅の尺度は

$$
z_\ell=\sqrt[3]{\left(1\right) \left(\frac{1}{8}\right)}=\frac{1}{2}\,\mathrm{m\,s^{-1}}
$$

である。

非線形時間は

$$
t_{\rm nl}=\frac{\left(\frac{1}{8}\right)}{\left(\frac{1}{2}\right)}=\frac{1}{4}\,\mathrm s
$$

である。

平行尺度は

$$
\ell_\parallel=\left(10\right) \left(\frac{1}{4}\right)=\frac{5}{2}\,\mathrm m
$$

である。

異方性の長さ比は

$$
\ell_\parallel/\ell_\perp=\frac{\left(\frac{5}{2}\right)}{\left(\frac{1}{8}\right)}=20
$$

である。

横波数の尺度は

$$
k_\perp=\frac{1}{\left(\frac{1}{8}\right)}=8\,\mathrm{m^{-1}}
$$

である。

横波数スペクトルの尺度は

$$
E(k_\perp)=\frac{\left(\frac{1}{2}\right)^{2}}{\left(8\right)}=\frac{1}{32}\,\mathrm{m^3\,s^{-2}}
$$

である。

平行長も短くなるが、横長さとの比は増す。
:::

:::example[同じ横尺度で小さなエネルギー流束]
振幅と平行長の変化を求める。

入力値を次のように置く。

$$
\begin{aligned}
\ell_\perp&=1\,\mathrm m \\\epsilon&=\frac{1}{8}\,\mathrm{m^2\,s^{-3}} \\v_A&=10\,\mathrm{m\,s^{-1}}
\end{aligned}
$$

これらを順に代入する。

波の振幅の尺度は

$$
z_\ell=\sqrt[3]{\left(1\right) \left(\frac{1}{8}\right)}=\frac{1}{2}\,\mathrm{m\,s^{-1}}
$$

である。

非線形時間は

$$
t_{\rm nl}=\frac{\left(1\right)}{\left(\frac{1}{2}\right)}=2\,\mathrm s
$$

である。

平行尺度は

$$
\ell_\parallel=\left(10\right) \left(2\right)=20\,\mathrm m
$$

である。

異方性の長さ比は

$$
\ell_\parallel/\ell_\perp=\frac{\left(20\right)}{\left(1\right)}=20
$$

である。

横波数の尺度は

$$
k_\perp=\frac{1}{\left(1\right)}=1\,\mathrm{m^{-1}}
$$

である。

横波数スペクトルの尺度は

$$
E(k_\perp)=\frac{\left(\frac{1}{2}\right)^{2}}{\left(1\right)}=\frac{1}{4}\,\mathrm{m^3\,s^{-2}}
$$

である。

振幅の低下で非線形時間が長くなり、釣り合う平行長も増す。
:::

:::example[強い案内磁場]
同じ伝達率のもとで必要な平行長を調べる。

入力値を次のように置く。

$$
\begin{aligned}
\ell_\perp&=1\,\mathrm m \\\epsilon&=1\,\mathrm{m^2\,s^{-3}} \\v_A&=100\,\mathrm{m\,s^{-1}}
\end{aligned}
$$

これらを順に代入する。

波の振幅の尺度は

$$
z_\ell=\sqrt[3]{\left(1\right)^{2}}=1\,\mathrm{m\,s^{-1}}
$$

である。

非線形時間は

$$
t_{\rm nl}=1\,\mathrm s
$$

である。

平行尺度は

$$
\ell_\parallel=\left(100\right) \left(1\right)=100\,\mathrm m
$$

である。

異方性の長さ比は

$$
\ell_\parallel/\ell_\perp=\frac{\left(100\right)}{\left(1\right)}=100
$$

である。

横波数の尺度は

$$
k_\perp=\frac{1}{\left(1\right)}=1\,\mathrm{m^{-1}}
$$

である。

横波数スペクトルの尺度は

$$
E(k_\perp)=\left(1\right)=1\,\mathrm{m^3\,s^{-2}}
$$

である。

横のスペクトルの形式値は同じで、線形伝播が速いため平行長が伸びる。
:::


## 演習問題

::::exercise{level=basic}
$D_v$ が一定で $\partial_vF=0$ の領域では準線形粒子束はどうなるか。

:::solution
流束は $-D_v\partial_vF$ でありゼロである。場の振幅が非ゼロでも勾配がなければこの流束は生じない。
:::
::::

::::exercise{level=basic}
波の振幅 $\abs\Phi$ が四倍になったとき往復角周波数を比較せよ。

:::solution
式 :eqref[eq-plasma-bounce] は $\abs\Phi^{1/2}$ に比例するので二倍である。
:::
::::

::::exercise{level=basic}
$U=1/2$ の KdV 孤立波の振幅を求めよ。

:::solution
孤立波の最大値は $3U$ なので $3/2$ である。
:::
::::

::::exercise{level=basic}
三波模型で $\abs g=2,\Delta=5$ に指数成長はあるか。

:::solution
$\gamma^2=4-25/4=-9/4<0$ なので、固有指数は純虚数になり指数成長しない。
:::
::::

::::exercise{level=advanced}
KdV の $\int \phi_1^2\dd\xi$ の保存を示せ。

:::solution
式 :eqref[eq-plasma-kdv] に $2\phi_1$ を掛ける。非線形項は $(2/3)\partial_\xi\phi_1^3$ となり、分散項は $\phi_1\partial_\xi^3\phi_1=\partial_\xi(\phi_1\partial_\xi^2\phi_1-(\partial_\xi\phi_1)^2/2)$ となる。周期または減衰境界で発散の積分はゼロである。
:::
::::

::::exercise{level=advanced}
垂直衝撃波で $\rho_2=3\rho_1$ のとき $u_2,B_{t2}$ を求めよ。

:::solution
質量流束一定より $u_2=u_1/3$、接線電場一定より $B_{t2}=3B_{t1}$ である。
:::
::::

::::exercise{level=advanced}
$U=2$ の孤立波が積分後の二階方程式を満たすことを確認せよ。

:::solution
$F=6\,\mathrm{sech}^2z$ なので $F''=24\,\mathrm{sech}^2z-36\,\mathrm{sech}^4z=4F-F^2$ となる。両端で $F$ と導関数が消える。
:::
::::

::::exercise{level=advanced}
準線形拡散による粒子運動エネルギーの変化を求めよ。

:::solution
$\int mv^2/2$ を重みとして式 :eqref[eq-plasma-ql] を積分し、境界での流束をゼロとすると $\dot W=-\int mv D_v \partial_vF\dd v$ である。正速度の単調減少部分だけが共鳴すれば粒子はエネルギーを得る。
:::
::::

## 確認テスト

::::quiz
有限の捕捉が線形減衰を修正しうる条件を選べ。

- [x] 粒子の往復が線形減衰と同程度以上に速い。
- [ ] 波振幅が厳密にゼロである。
- [ ] 全粒子の速度が波と大きく異なる。

:::explanation
捕捉の時間尺度 $1/\omega_b$ が減衰時間以下なら、分布は減衰の途中で変化できる。
:::
::::

::::quiz{answer="1" tol="1e-6"}
$\omega_0=100\,\mathrm{s^{-1}}$、$h=0.04$ の共鳴振動子の成長率を SI で答えよ。

:::explanation
包絡線の成長率は $\abs h\omega_0/4=0.04\cdot100/4=1\,\mathrm{s^{-1}}$ である。
:::
::::

::::quiz
上流の接線磁場が非零の垂直衝撃波で、理想上下流をつなぐ圧縮は何を変えるか。

- [x] 接線磁場の大きさが増える。
- [ ] 接線磁場が必ず消える。
- [ ] 質量流束が変わる。

:::explanation
$uB_t$ と $\rho u$ がともに保存されるので $B_t/\rho$ は一定である。従って圧縮で磁場の大きさが増える。
:::
::::

::::quiz{answer="9" tol="1e-6"}
KdV 孤立波の $U=3$ での無次元振幅を答えよ。

:::explanation
積分式の正の転回点は $F=3U$ なので $9$ である。
:::
::::

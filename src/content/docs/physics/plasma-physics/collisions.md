---
title: Coulomb 衝突と輸送
description: 二体散乱から輸送断面積を導き、抵抗と拡散を比較する。
notation:
  - symbol: 'b'
    meaning: 衝突径数
  - symbol: '\ln\Lambda'
    meaning: Coulomb 対数
  - symbol: '\nu'
    meaning: 運動量緩和周波数
---

## この章の目標

- 小角散乱の累積を計算する。
- 衝突時間と抵抗率の関係を導く。
- 磁場の有無で粒子束の向きを比較する。

## Coulomb 衝突と Coulomb 対数

荷電粒子は遠くを通過しても運動量を少し交換する。輸送を求めるときは、この小さな変化を断面積で重み付けして足す。古典的で弱結合の二種を考え、相対速度を $u$、換算質量を $m_r=m_a m_b/(m_a+m_b)$、相互作用を $\kappa/r$、$\kappa=q_aq_b/(4\pi\varepsilon_0)$ とする。

### 軌道と輸送断面積

相対軌道には角運動量 $\ell=m_ru b$ が保存される。$y=1/r$ を方位角の関数とすると、動径方程式に $\dot\theta=\ell y^2/m_r$ を代入して $y''+y=-m_r\kappa/\ell^2$ を得る。解は $y=A\cos(\theta-\theta_0)-m_r\kappa/\ell^2$ である。エネルギー $m_ru^2/2=\ell^2((y')^2+y^2)/(2m_r)+\kappa y$ より $A^2=b^{-2}+(m_r\kappa/\ell^2)^2$ となる。二つの漸近線 $y=0$ のなす角から、偏向角 $\chi$ の大きさは $\tan(\chi/2)=b_0/b$、$b_0=\abs\kappa/(m_ru^2)$ である。環の面積 $2\pi b\dd b$ を立体角 $2\pi\sin\chi\dd\chi$ で割ると Rutherford の微分断面積 $b_0^2/(4\sin^4(\chi/2))$ が得られる。

入射方向の速度の減少は $u(1-\cos\chi)$ である。したがって上限 $b_{\max}$ までの運動量輸送断面積は
$$
\sigma_m=\int_0^{b_{\max}}(1-\cos\chi)2\pi b\,\dd b
=2\pi b_0^2\ln(1+b_{\max}^2/b_0^2).
\label{eq-plasma-coulomb-cross}
$$
弱結合では $b_{\max}$ の目安を遮蔽長、下限の目安を $b_{\min}\simeq\max(b_0,\hbar/(m_ru))$ とする。後者の量子スケール以下では古典軌道で偏向を追えない。対数精度で $\sigma_m\simeq 4\pi b_0^2\ln\Lambda$、$\ln\Lambda=\ln(b_{\max}/b_{\min})$ である。これを :term[Coulomb 対数]{en="Coulomb logarithm" yomi="くーろんたいすう"} という。遮蔽が動的な場合の厳密な上限は応答関数に依存する。

:::example[対数の寄与]
$b_0\ll b_1<b_2\le b_{\max}$ で積分すると $4\pi b_0^2\ln(b_2/b_1)$ になる。$10b_0$ から $100b_0$ と、$100b_0$ から $1000b_0$ の寄与は小角近似で等しい。
:::

### 偏向の大小で積分を分ける

$b<b_0$ の寄与を式 :eqref[eq-plasma-coulomb-cross] で切り出すと $2\pi b_0^2\ln2$ であり、$b>b_0$ の対数が大きい場合に比べて小さい。このため、少数の大角散乱を落としても対数精度で輸送を再現できる。$\ln\Lambda$ が数程度しかないときは、強散乱と遮蔽の移行部を別々に見積もって足す近似の精度が下がる。

:::example[90 度偏向の領域]
$b_{\max}=100 b_0$ と仮定すると $b<b_0$ の寄与割合は $\ln2/\ln10001\simeq0.0753$ である。強い一回の散乱より、小さな偏向の総和が輸送を主に担う。
:::

### Rutherford 散乱の厳密積分と切断

相対速度 $u$、換算質量 $m_r$ の斥力軌道では $\tan(\theta/2)=b_0/b$、$b_0=q_1q_2/(4\pi\varepsilon_0 m_r u^2)>0$ である。輸送断面積には軌道の向きの変化 $1-\cos\theta=2b_0^2/(b^2+b_0^2)$ を掛ける。衝突径数 $b$ と $b+\dd b$ の入射面積は $2\pi b\,\dd b$ なので、上限 $b_m$ まで積分すると $2\pi b_0^2\ln[1+(b_m/b_0)^2]$ となる。全散乱断面積 $\pi b_m^2$ とは異なる。$b\le b_0$ だけの寄与は $2\pi b_0^2\ln2$ であり、上限を増やしても一定である。この計算は Debye 遮蔽を鋭い切断で置き換えた模型である。量子回折が $b_0$ より長い場合にはこの軌道の積分をそのまま使えない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $b_m/b_0$ | 最大衝突径数と九十度偏向の径数の比。$b_m\ge b_0$ とする。 | 1 |
| $b/b_0$ | 一つの軌道の衝突径数 | 1 |

散乱角を $\theta$ とする。半角を逆正接で計算する。

$$
\theta=\frac{360 \arctan{\left(\frac{1}{b/b_{0}} \right)}}{\pi}
$$

軌道の運動量輸送の重みを $w$ とする。角度の余弦差を用いる。

$$
w=\frac{2}{(b/b_0)^2 + 1}
$$

規格化輸送断面積を $\sigma_T/(2\pi b_0^2)$ とする。強散乱も含めて積分する。

$$
\sigma_T/(2\pi b_0^2)=\ln{\left((b_m/b_0)^{2} + 1 \right)}
$$

対数近似での断面積を $\sigma_{T,\rm log}/(2\pi b_0^2)$ とする。$b_m\gg b_0$ の先頭項である。

$$
\sigma_{T,\rm log}/(2\pi b_0^2)=2 \ln{\left((b_m/b_0) \right)}
$$

強散乱の寄与の割合を $f_{90}$ とする。積分上限を $b_0$ とした結果を全体で割る。

$$
f_{90}=\frac{\ln{\left(2 \right)}}{[\sigma_T/(2\pi b_0^2)]}
$$

入射一軌道当たりの平均輸送重みを $\sigma_T/(\pi b_m^2)$ とする。弱散乱が多くても平均重みは小さくなり得る。

$$
\sigma_T/(\pi b_m^2)=\frac{2 [\sigma_T/(2\pi b_0^2)]}{(b_m/b_0)^{2}}
$$

:::example[上限が九十度散乱の径数]
対数だけを残すと何を落とすか計算する。

入力値を次のように置く。

$$
\begin{aligned}
b_m/b_0&=1 \\b/b_0&=1
\end{aligned}
$$

これらを順に代入する。

散乱角は

$$
\theta=\frac{360 \arctan{\left(\frac{1}{\left(1\right)} \right)}}{\pi}=90\,\mathrm{deg}
$$

である。

軌道の運動量輸送の重みは

$$
w=\frac{2}{\left(1\right)^{2} + 1}=1
$$

である。

規格化輸送断面積は

$$
\sigma_T/(2\pi b_0^2)=\ln{\left(\left(1\right)^{2} + 1 \right)}\simeq 0.69315
$$

である。

対数近似での断面積は

$$
\sigma_{T,\rm log}/(2\pi b_0^2)=2 \ln{\left(\left(1\right) \right)}=0
$$

である。

強散乱の寄与の割合は

$$
f_{90}=\frac{\ln{\left(2 \right)}}{\left(\ln{\left(2 \right)}\right)}=1
$$

である。

入射一軌道当たりの平均輸送重みは

$$
\sigma_T/(\pi b_m^2)=\frac{2 \left(\ln{\left(2 \right)}\right)}{\left(1\right)^{2}}\simeq 1.3863
$$

である。

対数近似がゼロでも大角散乱は有限である。弱結合の衝突論をこの結果から適用できない。
:::

:::example[一桁広い範囲]
最大径数の軌道の角度と全体の輸送を比較する。

入力値を次のように置く。

$$
\begin{aligned}
b_m/b_0&=10 \\b/b_0&=10
\end{aligned}
$$

これらを順に代入する。

散乱角は

$$
\theta=\frac{360 \arctan{\left(\frac{1}{\left(10\right)} \right)}}{\pi}\simeq 11.421\,\mathrm{deg}
$$

である。

軌道の運動量輸送の重みは

$$
w=\frac{2}{\left(10\right)^{2} + 1}=\frac{2}{101}
$$

である。

規格化輸送断面積は

$$
\sigma_T/(2\pi b_0^2)=\ln{\left(\left(10\right)^{2} + 1 \right)}\simeq 4.6151
$$

である。

対数近似での断面積は

$$
\sigma_{T,\rm log}/(2\pi b_0^2)=2 \ln{\left(\left(10\right) \right)}\simeq 4.6052
$$

である。

強散乱の寄与の割合は

$$
f_{90}=\frac{\ln{\left(2 \right)}}{\left(\ln{\left(101 \right)}\right)}\simeq 0.15019
$$

である。

入射一軌道当たりの平均輸送重みは

$$
\sigma_T/(\pi b_m^2)=\frac{2 \left(\ln{\left(101 \right)}\right)}{\left(10\right)^{2}}\simeq 0.092302
$$

である。

個々の最外軌道は弱い偏向だが、多数の遠距離軌道が輸送に寄与する。
:::

:::example[広い範囲にわたる積分]
一つの強散乱と積分での強散乱の分担を比べる。

入力値を次のように置く。

$$
\begin{aligned}
b_m/b_0&=1000 \\b/b_0&=\frac{1}{2}
\end{aligned}
$$

これらを順に代入する。

散乱角は

$$
\theta=\frac{360 \arctan{\left(\frac{1}{\left(\frac{1}{2}\right)} \right)}}{\pi}\simeq 126.87\,\mathrm{deg}
$$

である。

軌道の運動量輸送の重みは

$$
w=\frac{2}{\left(\frac{1}{2}\right)^{2} + 1}=\frac{8}{5}
$$

である。

規格化輸送断面積は

$$
\sigma_T/(2\pi b_0^2)=\ln{\left(\left(1000\right)^{2} + 1 \right)}\simeq 13.816
$$

である。

対数近似での断面積は

$$
\sigma_{T,\rm log}/(2\pi b_0^2)=2 \ln{\left(\left(1000\right) \right)}\simeq 13.816
$$

である。

強散乱の寄与の割合は

$$
f_{90}=\frac{\ln{\left(2 \right)}}{\left(\ln{\left(1000001 \right)}\right)}\simeq 0.050172
$$

である。

入射一軌道当たりの平均輸送重みは

$$
\sigma_T/(\pi b_m^2)=\frac{2 \left(\ln{\left(1000001 \right)}\right)}{\left(1000\right)^{2}}\simeq 2.7631\times 10^{-5}
$$

である。

強散乱は向きを大きく変える一方、積分への寄与は一割より小さい。
:::

:::example[中間の径数での弱散乱]
断面積を全衝突の面積で割ったときの尺度を求める。

入力値を次のように置く。

$$
\begin{aligned}
b_m/b_0&=100 \\b/b_0&=10
\end{aligned}
$$

これらを順に代入する。

散乱角は

$$
\theta=\frac{360 \arctan{\left(\frac{1}{\left(10\right)} \right)}}{\pi}\simeq 11.421\,\mathrm{deg}
$$

である。

軌道の運動量輸送の重みは

$$
w=\frac{2}{\left(10\right)^{2} + 1}=\frac{2}{101}
$$

である。

規格化輸送断面積は

$$
\sigma_T/(2\pi b_0^2)=\ln{\left(\left(100\right)^{2} + 1 \right)}\simeq 9.2104
$$

である。

対数近似での断面積は

$$
\sigma_{T,\rm log}/(2\pi b_0^2)=2 \ln{\left(\left(100\right) \right)}\simeq 9.2103
$$

である。

強散乱の寄与の割合は

$$
f_{90}=\frac{\ln{\left(2 \right)}}{\left(\ln{\left(10001 \right)}\right)}\simeq 0.075257
$$

である。

入射一軌道当たりの平均輸送重みは

$$
\sigma_T/(\pi b_m^2)=\frac{2 \left(\ln{\left(10001 \right)}\right)}{\left(100\right)^{2}}\simeq 0.0018421
$$

である。

平均重みは小さいため、遭遇回数をすべて九十度衝突と数えると緩和が過大になる。
:::


## 衝突周波数と緩和時間

一つの緩和時間ですべての衝突効果は表せない。重い背景に軽い試験粒子を入れる模型では、試験粒子の速さを固定して向きだけの緩和を測れる。時間 $\dd t$ の間に $n_b u\,\dd\sigma\,\dd t$ 回の散乱が見込まれるため、速度の方向平均の減衰率は $\nu_m(u)=n_bu\sigma_m$ である。:eqref[eq-plasma-coulomb-cross] から $\nu_m\propto u^{-3}\ln\Lambda$ を得る。

### エネルギー交換の時間

静止標的に渡る運動エネルギーは反跳運動量 $2m_ru\sin(\chi/2)$ の二乗を $2m_b$ で割った値である。入射エネルギー $m_a u^2/2$ に対する割合は $2m_r^2(1-\cos\chi)/(m_am_b)$ である。$m_a\ll m_b$ なら、エネルギー緩和は方向緩和よりおよそ $m_b/(2m_a)$ 倍遅い。熱い背景については相対速度を積分し、背景からのエネルギー流入も含める。

:::example[高速粒子]
同じ $\ln\Lambda$ を使えば、速さが二倍の試験粒子の $\nu_m$ は八分の一になる。平均自由行程 $u/\nu_m$ は十六倍になる。
:::

:::example[異なる緩和]
静止陽子に電子が衝突する模型では $2m_e/m_p\simeq 1.09\times10^{-3}$ である。ほぼ同じエネルギーのまま進行方向を何度も変えられる。
:::

### 移流と衝突の比較

熱速度を $v_t=\sqrt{T/m}$ とすると、平均自由行程を $\ell_{\rm mfp}\sim v_t/\nu_m(v_t)$ で見積もれる。局所閉包の条件の一つは $L\gg\ell_{\rm mfp}$ である。一方、緩和率と観測角周波数 $\omega$ の比較は一周期中に衝突が起こるかを判定する。両方を同じ判定とみなせるのは $\omega\sim v_t/L$ の場合である。

:::example[速度分布の高温尾部]
熱速度で $\ell_{\rm mfp}/L=1/100$ の場合、速度 $4v_t$ では対数を固定した同じ模型で $\ell_{\rm mfp}/L=4^4/100=2.56$ となる。バルクが局所的に振る舞っても、尾部は境界まで届きうる。
:::

### 熱平均と試験粒子の緩和を区別する

重く静止した散乱体の中を通る試験粒子に対し、輸送断面積を $\sigma_T=\sigma_m$ と書き、$\nu_T=n_b u\sigma_T$ を定義する。対数を一定とおくと $b_0\propto u^{-2}$ なので $\nu_T\propto u^{-3}$ である。飛行時間の逆数に対応する $\nu_T$ と、距離を測る平均自由行程 $\ell_T=u/\nu_T$ の速度依存性は異なる。散乱で角度だけが変わる極限では試験粒子の運動エネルギーは保たれる。有限質量 $M$ で標的の反跳を含むと、標的へ渡るエネルギーは運動量移行 $2m_r u\sin(\theta/2)$ の二乗を $2M$ で割った値であり、試験粒子の初期エネルギーに対する比は $2mM(1-\cos\theta)/(m+M)^2$ となる。従って重い標的への角度緩和とエネルギー緩和を同じ時間定数で置けない。動く熱背景で低速の粒子を扱う場合には、相対速度を背景分布で平均しなければならない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $x$ | 基準相対速度に対する速度比 $u/u_0$ | 1 |
| $\tau$ | 基準輸送周波数を掛けた観測時間 $\nu_{T0}t$ | 1 |
| $r$ | 標的と試験粒子の質量比 $M/m$ | 1 |

角度輸送の周波数比を $\nu_T/\nu_{T0}$ とする。Coulomb 対数の速度変化を無視する。

$$
\nu_T/\nu_{T0}=\frac{1}{x^{3}}
$$

平均自由行程の倍率を $\ell_T/\ell_{T0}$ とする。$\ell_T=u/\nu_T$ を使う。

$$
\ell_T/\ell_{T0}=x^{4}
$$

観測時間内の角度緩和数を $N_T$ とする。時間と周波数の積である。

$$
N_T=\nu_T/\nu_{T0} \tau
$$

輸送の重みに対するエネルギー移行の係数を $C_E$ とする。標的の反跳を含む二体力学の値である。

$$
C_E=\frac{2 r}{\left(r + 1\right)^{2}}
$$

初期のエネルギー緩和率と角度緩和率の比を $\nu_{E0}/\nu_T$ とする。減速に伴う周波数変化を含める前の比である。

$$
\nu_{E0}/\nu_T=C_{E}
$$

方向の一次相関を $\exp(-N_T)$ とする。速さの変化を無視できる Lorentz 極限の値を比較用に用いる。

$$
\exp(-N_T)=e^{- N_{T}}
$$

:::example[重いイオンに当たる電子]
方向が緩和する時間でのエネルギー交換の尺度を求める。

入力値を次のように置く。

$$
\begin{aligned}
x&=1 \\\tau&=1 \\r&=1836
\end{aligned}
$$

これらを順に代入する。

角度輸送の周波数比は

$$
\nu_T/\nu_{T0}=\frac{1}{\left(1\right)^{3}}=1
$$

である。

平均自由行程の倍率は

$$
\ell_T/\ell_{T0}=\left(1\right)^{4}=1
$$

である。

観測時間内の角度緩和数は

$$
N_T=\left(1\right)^{2}=1
$$

である。

輸送の重みに対するエネルギー移行の係数は

$$
C_E=\frac{2 \left(1836\right)}{\left(\left(1836\right) + 1\right)^{2}}=\frac{3672}{3374569}
$$

である。

初期のエネルギー緩和率と角度緩和率の比は

$$
\nu_{E0}/\nu_T=\left(\frac{3672}{3374569}\right)=\frac{3672}{3374569}
$$

である。

方向の一次相関は

$$
\exp(-N_T)=e^{- \left(1\right)}\simeq 0.36788
$$

である。

一次相関は大きく減るが、標的へ渡る比は千分の一程度である。
:::

:::example[速い尾部]
同じ観測時間で高速電子の角度緩和を調べる。

入力値を次のように置く。

$$
\begin{aligned}
x&=4 \\\tau&=1 \\r&=1836
\end{aligned}
$$

これらを順に代入する。

角度輸送の周波数比は

$$
\nu_T/\nu_{T0}=\frac{1}{\left(4\right)^{3}}=\frac{1}{64}
$$

である。

平均自由行程の倍率は

$$
\ell_T/\ell_{T0}=\left(4\right)^{4}=256
$$

である。

観測時間内の角度緩和数は

$$
N_T=\left(1\right) \left(\frac{1}{64}\right)=\frac{1}{64}
$$

である。

輸送の重みに対するエネルギー移行の係数は

$$
C_E=\frac{2 \left(1836\right)}{\left(\left(1836\right) + 1\right)^{2}}=\frac{3672}{3374569}
$$

である。

初期のエネルギー緩和率と角度緩和率の比は

$$
\nu_{E0}/\nu_T=\left(\frac{3672}{3374569}\right)=\frac{3672}{3374569}
$$

である。

方向の一次相関は

$$
\exp(-N_T)=e^{- \left(\frac{1}{64}\right)}\simeq 0.9845
$$

である。

平均自由行程は $256$ 倍であり、熱速度で局所性が成り立っても高速成分には再確認が要る。
:::

:::example[同程度の質量]
反跳が無視できない散乱で角度とエネルギーの尺度を比べる。

入力値を次のように置く。

$$
\begin{aligned}
x&=1 \\\tau&=1 \\r&=1
\end{aligned}
$$

これらを順に代入する。

角度輸送の周波数比は

$$
\nu_T/\nu_{T0}=\frac{1}{\left(1\right)^{3}}=1
$$

である。

平均自由行程の倍率は

$$
\ell_T/\ell_{T0}=\left(1\right)^{4}=1
$$

である。

観測時間内の角度緩和数は

$$
N_T=\left(1\right)^{2}=1
$$

である。

輸送の重みに対するエネルギー移行の係数は

$$
C_E=\frac{2 \left(1\right)}{\left(\left(1\right) + 1\right)^{2}}=\frac{1}{2}
$$

である。

初期のエネルギー緩和率と角度緩和率の比は

$$
\nu_{E0}/\nu_T=\left(\frac{1}{2}\right)=\frac{1}{2}
$$

である。

方向の一次相関は

$$
\exp(-N_T)=e^{- \left(1\right)}\simeq 0.36788
$$

である。

エネルギー変化は小さな補正で済まないため、固定速さの方向相関は比較値に留まる。
:::

:::example[方向緩和が先行する時間帯]
速い電子を長く観測した場合の緩和数を求める。

入力値を次のように置く。

$$
\begin{aligned}
x&=2 \\\tau&=16 \\r&=1836
\end{aligned}
$$

これらを順に代入する。

角度輸送の周波数比は

$$
\nu_T/\nu_{T0}=\frac{1}{\left(2\right)^{3}}=\frac{1}{8}
$$

である。

平均自由行程の倍率は

$$
\ell_T/\ell_{T0}=\left(2\right)^{4}=16
$$

である。

観測時間内の角度緩和数は

$$
N_T=\left(16\right) \left(\frac{1}{8}\right)=2
$$

である。

輸送の重みに対するエネルギー移行の係数は

$$
C_E=\frac{2 \left(1836\right)}{\left(\left(1836\right) + 1\right)^{2}}=\frac{3672}{3374569}
$$

である。

初期のエネルギー緩和率と角度緩和率の比は

$$
\nu_{E0}/\nu_T=\left(\frac{3672}{3374569}\right)=\frac{3672}{3374569}
$$

である。

方向の一次相関は

$$
\exp(-N_T)=e^{- \left(2\right)}\simeq 0.13534
$$

である。

方向相関は小さい。一方、初期のエネルギー緩和率を使った累積移行割合の見積もり $C_E N_T$ は一パーセント未満である。
:::


## Spitzer 抵抗率

電流を担う電子は電場で加速される一方、イオンへ運動量を渡す。電子の相対流速を $\vb u_e-\vb u_i$ とし、摩擦を $-m_e n_e\nu(\vb u_e-\vb u_i)$ と近似する。定常かつ一様なら $-en_e\vb E-m_e n_e\nu(\vb u_e-\vb u_i)=0$ であり、$\vb j=-en_e(\vb u_e-\vb u_i)$ から $\eta=m_e\nu/(n_e e^2)$ を得る。

熱速度での散乱から $\nu\propto n_i Z^2T_e^{-3/2}\ln\Lambda$ が出る。$n_e=Zn_i$ を使うと $\eta\propto ZT_e^{-3/2}\ln\Lambda$ である。係数を求めるには速度依存の衝突率を使い、電子同士が分布を変える効果も含める。平行方向の :term[Spitzer 抵抗率]{en="Spitzer resistivity" yomi="すぴっつぁーていこうりつ"} は、熱平均した電子イオン衝突時間を $\tau_e$ とする慣例で水素に対し約 $0.51m_e/(n_e e^2\tau_e)$ である。定義は [衝突の運動論](/physics/plasma-physics/collisional-kinetics/) に示す。$0.51$ は速度分布の非平衡補正を解いた値であり、ここでは温度依存の導出に限る。係数を与える線形衝突演算子の解法は :cite[goldston-rutherford] を参照する。

:::example[加熱後の抵抗]
$\ln\Lambda$ を一定とすると温度が四倍で抵抗率は八分の一になる。同じ電流密度を流すときの加熱率 $\eta j^2$ も八分の一である。
:::

### 伝導の時間応答

摩擦を一定とする一様場で $\dot{\vb j}+\nu\vb j=(n_e e^2/m_e)\vb E$ である。したがって $\sigma(\omega)=\sigma_0/(1-i\omega/\nu)$、$\sigma_0=n_e e^2/(m_e\nu)$ となる。平均加熱に寄与する $\Re\sigma=\sigma_0/[1+(\omega/\nu)^2]$ は、位相を無視した電流振幅の比から直接は得られない。

:::example[衝突が少ない交流伝導]
$\omega=3\nu$ では同じ電場振幅に対する平均散逸は定常伝導度を使った値の $1/10$ である。電流振幅の減少率 $1/\sqrt{10}$ とは異なる。
:::

### Drude 応答から抵抗性と誘導性を分ける

運動量輸送の単一周波数 $\nu>0$ で電子を記述し、電場を $\Re(E_0e^{-i\omega t})$ とする。$m_e(\partial_t+\nu)u_e=-eE$ を解いて $\sigma=n_e e^2/[m_e(\nu-i\omega)]$ を得る。実部は電場と同相、虚部は四分の一周期ずれた応答である。一周期平均の加熱は $\Re(\sigma)|E_0|^2/2$ で、粒子慣性があるので定常抵抗率の逆数だけでは決まらない。抵抗率の $T_e^{-3/2}$ という依存性は、$\nu\propto n_e T_e^{-3/2}$ を $\eta=m_e\nu/(n_e e^2)$ に代入すれば分かる。これはスケーリングの説明であり、Spitzer の数係数を得るには速度依存衝突と電子同士の衝突を含めて線形化する。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $a$ | 電子温度の基準値に対する倍率 | 1 |
| $N$ | 電子密度の基準値に対する倍率 | 1 |
| $W$ | 基準衝突周波数に対する駆動周波数 $\omega/\nu_0$ | 1 |

衝突周波数の倍率を $\nu/\nu_0$ とする。価数と Coulomb 対数を固定する。

$$
\nu/\nu_0=\frac{N}{a^{\frac{3}{2}}}
$$

定常抵抗率の倍率を $\eta/\eta_0$ とする。密度因子が分子と分母で打ち消す。

$$
\eta/\eta_0=\frac{\nu/\nu_{0}}{N}
$$

駆動と衝突の比を $r$ とする。変更後の衝突周波数を分母にする。

$$
r=\frac{W}{\nu/\nu_{0}}
$$

基準直流伝導度で規格化した実部を $\Re(\sigma)/\sigma_0$ とする。同相電流が平均の加熱を担う。

$$
\Re(\sigma)/\sigma_0=\frac{1}{\eta/\eta_{0} \left(r^{2} + 1\right)}
$$

規格化した虚部を $\Im(\sigma)/\sigma_0$ とする。$e^{-i\omega t}$ の位相規約を使う。

$$
\Im(\sigma)/\sigma_0=\Re(\sigma)/\sigma_{0} r
$$

同じ媒質の直流伝導度から見積もった電力に対する割合を $P/P_{\rm dc}$ とする。比較にはどちらも同じ実効電場を用いる。

$$
P/P_{\rm dc}=\frac{1}{r^{2} + 1}
$$

:::example[直流で四倍に加熱]
電子温度が上がったときの抵抗率と電流を比べる。

入力値を次のように置く。

$$
\begin{aligned}
a&=4 \\N&=1 \\W&=0
\end{aligned}
$$

これらを順に代入する。

衝突周波数の倍率は

$$
\nu/\nu_0=\frac{\left(1\right)}{\left(4\right)^{\frac{3}{2}}}=\frac{1}{8}
$$

である。

定常抵抗率の倍率は

$$
\eta/\eta_0=\frac{\left(\frac{1}{8}\right)}{\left(1\right)}=\frac{1}{8}
$$

である。

駆動と衝突の比は

$$
r=\frac{\left(0\right)}{\left(\frac{1}{8}\right)}=0
$$

である。

基準直流伝導度で規格化した実部は

$$
\Re(\sigma)/\sigma_0=\frac{1}{\left(\frac{1}{8}\right) \left(\left(0\right)^{2} + 1\right)}=8
$$

である。

規格化した虚部は

$$
\Im(\sigma)/\sigma_0=\left(0\right) \left(8\right)=0
$$

である。

同じ媒質の直流伝導度から見積もった電力に対する割合は

$$
P/P_{\rm dc}=\frac{1}{\left(0\right)^{2} + 1}=1
$$

である。

基準の八分の一の定常抵抗率となる。
:::

:::example[密度を増やして同じ温度を保つ]
直流抵抗率と交流の加熱率が同じ応答を示すか調べる。

入力値を次のように置く。

$$
\begin{aligned}
a&=1 \\N&=4 \\W&=1
\end{aligned}
$$

これらを順に代入する。

衝突周波数の倍率は

$$
\nu/\nu_0=\frac{\left(4\right)}{\left(1\right)^{\frac{3}{2}}}=4
$$

である。

定常抵抗率の倍率は

$$
\eta/\eta_0=1
$$

である。

駆動と衝突の比は

$$
r=\frac{\left(1\right)}{\left(4\right)}=\frac{1}{4}
$$

である。

基準直流伝導度で規格化した実部は

$$
\Re(\sigma)/\sigma_0=\frac{1}{\left(1\right) \left(\left(\frac{1}{4}\right)^{2} + 1\right)}=\frac{16}{17}
$$

である。

規格化した虚部は

$$
\Im(\sigma)/\sigma_0=\left(\frac{16}{17}\right) \left(\frac{1}{4}\right)=\frac{4}{17}
$$

である。

同じ媒質の直流伝導度から見積もった電力に対する割合は

$$
P/P_{\rm dc}=\frac{1}{\left(\frac{1}{4}\right)^{2} + 1}=\frac{16}{17}
$$

である。

定常抵抗率は一定だが、相対的に駆動が遅くなり電流実部は直流値に近づく。
:::

:::example[加熱によって駆動が高速側へ移る場合]
直流伝導度が向上しても同相電流が増えるか比較する。

入力値を次のように置く。

$$
\begin{aligned}
a&=4 \\N&=1 \\W&=1
\end{aligned}
$$

これらを順に代入する。

衝突周波数の倍率は

$$
\nu/\nu_0=\frac{\left(1\right)}{\left(4\right)^{\frac{3}{2}}}=\frac{1}{8}
$$

である。

定常抵抗率の倍率は

$$
\eta/\eta_0=\frac{\left(\frac{1}{8}\right)}{\left(1\right)}=\frac{1}{8}
$$

である。

駆動と衝突の比は

$$
r=\frac{\left(1\right)}{\left(\frac{1}{8}\right)}=8
$$

である。

基準直流伝導度で規格化した実部は

$$
\Re(\sigma)/\sigma_0=\frac{1}{\left(\frac{1}{8}\right) \left(\left(8\right)^{2} + 1\right)}=\frac{8}{65}
$$

である。

規格化した虚部は

$$
\Im(\sigma)/\sigma_0=\left(8\right) \left(\frac{8}{65}\right)=\frac{64}{65}
$$

である。

同じ媒質の直流伝導度から見積もった電力に対する割合は

$$
P/P_{\rm dc}=\frac{1}{\left(8\right)^{2} + 1}=\frac{1}{65}
$$

である。

衝突が少なすぎると電流は主に虚部であり、定常伝導度の増加だけでは加熱を説明できない。
:::

:::example[四分の一密度の交流応答]
密度を下げたときの虚部と実部を求める。

入力値を次のように置く。

$$
\begin{aligned}
a&=1 \\N&=\frac{1}{4} \\W&=1
\end{aligned}
$$

これらを順に代入する。

衝突周波数の倍率は

$$
\nu/\nu_0=\frac{\left(\frac{1}{4}\right)}{\left(1\right)^{\frac{3}{2}}}=\frac{1}{4}
$$

である。

定常抵抗率の倍率は

$$
\eta/\eta_0=1
$$

である。

駆動と衝突の比は

$$
r=\frac{\left(1\right)}{\left(\frac{1}{4}\right)}=4
$$

である。

基準直流伝導度で規格化した実部は

$$
\Re(\sigma)/\sigma_0=\frac{1}{\left(1\right) \left(\left(4\right)^{2} + 1\right)}=\frac{1}{17}
$$

である。

規格化した虚部は

$$
\Im(\sigma)/\sigma_0=\left(4\right) \left(\frac{1}{17}\right)=\frac{4}{17}
$$

である。

同じ媒質の直流伝導度から見積もった電力に対する割合は

$$
P/P_{\rm dc}=\frac{1}{\left(4\right)^{2} + 1}=\frac{1}{17}
$$

である。

虚部が実部の四倍であり、平均加熱を電流の絶対値から計算すると過大になる。
:::


### 衝突時間の規約をそろえた伝導度

電子温度をエネルギー単位で表し、$v_t=\sqrt{T_e/m_e}$ とする。$q_i=Ze$ のイオンを静止散乱体とした極限では、速度 $v$ での角度拡散率は $\nu_D(v)=4\pi n_i\kappa^2\ln\Lambda/(m_e^2v^3)$ である。ここで $\kappa=Ze^2/(4\pi\varepsilon_0)$ の大きさを使った。熱速さでの値 $\nu_*=\nu_D(v_t)$ と、電子イオン衝突時間の慣例
$$
\tau_e=\frac{3 m_e^{1/2}T_e^{3/2}}{4\sqrt{2\pi}n_i\kappa^2\ln\Lambda}
$$
の積は $3\sqrt{\pi/2}$ である。伝導度を求めるには、平行電場 $E$ による小さな分布変化 $f_1$ を考える。等方な背景 Maxwell 分布 $f_0(v)$ に対し、電場が駆動する成分は方向角の余弦 $v_z/v$ に比例する。この成分は等方な散乱で率 $\nu_D(v)$ だけ減衰する。散乱後の軸方向速度の平均が、方位角の平均により散乱前の値に $\cos\chi$ を掛けたものになるためである。したがって定常線形方程式は $-(eE/m_e)\partial_{v_z}f_0=-\nu_D f_1$ であり、
$$
f_1=-\frac{eE v_z}{T_e\nu_D(v)} f_0,qquad
\sigma_{\rm L}=\frac{4\pi e^2}{3T_e}\int_0^\infty\frac{v^4 f_0(v)}{\nu_D(v)}\,\dd v .
$$
角度積分には $\int\cos^2\vartheta\,\dd\Omega=4\pi/3$ を用いた。$f_0=n_e(2\pi v_t^2)^{-3/2}e^{-v^2/(2v_t^2)}$ と $\nu_D(v)=\nu_*(v_t/v)^3$ を入れると、残る積分は $\int_0^\infty v^7 e^{-v^2/(2v_t^2)}\dd v=48v_t^8$ である。従って $\sigma_{\rm L}=16\sqrt{2/\pi}\,n_e e^2/(m_e\nu_*)$ となる。同じ $\tau_e$ で伝導度と抵抗率を書くと
$$
\sigma_{\rm L}=\frac{32}{3\pi}\frac{n_e e^2\tau_e}{m_e},
\qquad
\eta_{\rm L}=\frac{3\pi}{32}\frac{m_e}{n_e e^2\tau_e}.
$$
これは伝導を担う速度が一つに固定されず、$\nu_D$ の小さい高速側を強く重み付けすることを示す。単一速度の Drude 模型で $\nu_D(v_t)$ を代入して求めた伝導度とはおよそ $12.77$ 倍異なる。

水素の $0.51$ という Spitzer 抵抗率の係数との違いは、電子同士の小角散乱が電流の速度分布を変えることから生じる。電子同士の衝突は電子全体の運動量を保存するが、電子からイオンへの運動量移行率は速度に依存する。電子種の全運動量が衝突で保存されるという事実だけから、電子同士の衝突を線形な定常分布の計算で省くことはできない。

:::example[同じ衝突時間での係数比較]
$\eta_{\rm L}/[m_e/(n_e e^2\tau_e)]=3\pi/32\simeq0.2945$ である。水素の Spitzer 値 $0.51$ をこれで割ると約 $1.73$ となる。速さの代表値を変えた比較ではなく、同じ $\tau_e$ で衝突演算子を変えた比較である。
:::

:::example[高い価数での近似の見通し]
準中性 $n_iZ=n_e$ により電子イオン衝突の密度因子 $n_i Z^2$ は $Z n_e$ となる。電子間衝突の密度因子は $n_e$ なので、同程度の Coulomb 対数では両率の比は $\order{1/Z}$ である。例えば $Z=10$ なら密度・電荷因子での比は $0.1$ であり、Lorentz 極限に近づく。これだけで厳密な抵抗率の誤差を一割と評価することはできない。
:::

## 弱電離プラズマの拡散と両極性拡散

中性粒子が多い放電では、イオンと電子の摩擦を中性気体との衝突で決める。一定温度、静止中性粒子、無磁場で慣性を無視し、種 $s$ の数束を
$$
\vb\Gamma_s=\sgn(q_s)\mu_s n_s\vb E-D_s\nabla n_s,\qquad
\mu_s=\frac{\abs{q_s}}{m_s\nu_s},\quad D_s=\frac{T_s}{m_s\nu_s}
$$
とする。これは力のつりあい $q_s n_s\vb E-T_s\nabla n_s-m_s\nu_s\vb\Gamma_s=0$ である。

電荷保存で電荷密度の変化を無視できる一次元準中性系では、電流の法線成分は位置によらない。境界でこの電流がゼロなら $\Gamma_e=\Gamma_i$ である。差をゼロと置くと
$$
E=\frac{D_i-D_e}{\mu_i+\mu_e}\partial_x\ln n,\qquad
\Gamma_i=\Gamma_e=-D_A\partial_x n,\qquad
D_A=\frac{\mu_iD_e+\mu_eD_i}{\mu_i+\mu_e}.
\label{eq-plasma-ambipolar}
$$
これを :term[両極性拡散]{en="ambipolar diffusion" yomi="りょうきょくせいかくさん"} と呼ぶ。$\mu_e\gg\mu_i$ の場合は $D_A\simeq D_i(1+T_e/T_i)$ となる。電子が速く広がろうとして生じる電場がイオンを加速する。

:::example[同温度]
$\mu_e=9\mu_i$、$T_i=T_e$ なら $D_e=9D_i$ なので $D_A=1.8D_i$ である。
:::

:::example[壁への減衰]
区間 $0<x<L$ で $n(0,t)=n(L,t)=0$、$D_A$ 一定なら $n=n_1\sin(\pi x/L)e^{-D_A(\pi/L)^2t}$ が $\partial_t n=D_A\partial_x^2n$ を満たす。壁シースの外側から見た有効境界の模型である。
:::

### 横磁場中の両極性電場

一次元の密度勾配が $\vb B$ に垂直で、中性粒子が静止する場合を考える。電場も勾配方向だけにあるとする。摩擦つき横運動の解から、勾配方向の移動度は $\mu_{s\perp}=\abs{q_s}\nu_s/[m_s(\nu_s^2+\Omega_s^2)]$、拡散係数は $D_{s\perp}=(T_s/\abs{q_s})\mu_{s\perp}$ である。勾配に直交する束は別に生じるが、全ての係数が一次元に依存するならその成分の発散はゼロである。よって法線電流がゼロという境界条件では、式 :eqref[eq-plasma-ambipolar] の $\mu_s,D_s$ を各横方向係数に置き換えればよい。これは両種の全ベクトルの粒子束が一致するという主張ではない。

無磁場で $\mu_e\gg\mu_i$ になりやすい場合でも、磁場をかければその順が逆転できる。電子だけが強磁化されるなら電子の法線移動度は $B^{-2}$ で減り、イオンの移動度は弱磁化の範囲でほぼ一定である。$\mu_{e\perp}\ll\mu_{i\perp}$ では $D_{A\perp}\simeq D_{e\perp}(1+T_i/T_e)$ となる。イオンが速く広がろうとする側での電場の向きには $D_{i\perp}-D_{e\perp}$ の符号も要り、移動度比だけでは決まらない。

:::example[電場が逆向きになる条件]
$T_i=T_e$ で無磁場の比が $\mu_e/\mu_i=100$ とする。イオンの磁化を無視し、$\abs{\Omega_e}/\nu_e=20$ とすれば横の比は $100/401<1$ となる。Einstein 関係により拡散係数の大小も逆転し、同じ密度勾配に対する両極性電場は無磁場の場合と反対になる。横の比を $x$ とすると $D_{A\perp}/D_{i\perp}=2x/(1+x)=200/501$ である。
:::

### 両極性拡散の電位と損失固有値

電子と一価イオンの流束を $\Gamma_i=\mu_i nE-D_i n'$、$\Gamma_e=-\mu_e nE-D_e n'$ と書く。両方が同じになる条件は $(\mu_e+\mu_i)nE=(D_i-D_e)n'$ である。これを戻すと $D_a=(D_i\mu_e+D_e\mu_i)/(\mu_e+\mu_i)$ が得られる。$D_s=\mu_s T_s/e$ のとき、$\mu_e\gg\mu_i$ では $D_a\simeq\mu_i(T_i+T_e)/e$ である。吸収端の一次元区間 $0<x<L$ で、一定の $D_a$ に対する $\partial_t n=D_a n''$ の固有関数は $\sin(n\pi x/L)$、減衰率は $D_a(n\pi/L)^2$ である。壁での密度ゼロは拡散模型の有効境界条件である。そこの局所電場を $\ln n$ の勾配から求めると発散するため、壁直近にはこの体積模型を延長しない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $r$ | 移動度比 $\mu_e/\mu_i$ | 1 |
| $t$ | 温度比 $T_e/T_i$ | 1 |
| $j$ | 閉じ込め区間の空間モード番号 | 1 |
| $a$ | 観測時間をイオン拡散時間で割った値 $D_it_{\rm obs}/L^2$ | 1 |

電子とイオンの拡散比を $D_e/D_i$ とする。Einstein 関係を両種に使う。

$$
D_e/D_i=r t
$$

両極性拡散係数の倍率を $D_a/D_i$ とする。相対流を電場で消した係数である。

$$
D_a/D_i=\frac{D_e/D_i + r}{r + 1}
$$

密度勾配に対する電場係数を $eE/(T_i\partial_x\ln n)$ とする。電子が速く拡散する場合には勾配と逆向きになる。

$$
eE/(T_i\partial_x\ln n)=\frac{1 - D_e/D_i}{r + 1}
$$

規格化した密度減衰率を $\gamma_j L^2/D_i$ とする。二階微分の固有値を使う。

$$
\gamma_j L^2/D_i=\pi^{2} (D_a/D_i) j^{2}
$$

規格化減衰時間を $\tau_j D_i/L^2$ とする。振幅が $1/e$ になる時間である。

$$
\tau_j D_i/L^2=\frac{1}{(\gamma_j L^2/D_i)}
$$

モード振幅の残存割合を $n_j(t_{\rm obs})/n_j(0)$ とする。線形方程式なのでモードごとに指数で減衰する。

$$
n_j(t_{\rm obs})/n_j(0)=e^{- (\gamma_j L^2/D_i) a}
$$

:::example[電子移動度が大きい放電]
イオン単独の拡散より損失が速いか求める。

入力値を次のように置く。

$$
\begin{aligned}
r&=100 \\t&=10 \\j&=1 \\a&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

電子とイオンの拡散比は

$$
D_e/D_i=\left(100\right) \left(10\right)=1000
$$

である。

両極性拡散係数の倍率は

$$
D_a/D_i=\frac{\left(1000\right) + \left(100\right)}{\left(100\right) + 1}=\frac{1100}{101}
$$

である。

密度勾配に対する電場係数は

$$
eE/(T_i\partial_x\ln n)=\frac{1 - \left(1000\right)}{\left(100\right) + 1}=- \frac{999}{101}
$$

である。

規格化した密度減衰率は

$$
\gamma_j L^2/D_i=\pi^{2} \left(1\right)^{2} \left(\frac{1100}{101}\right)\simeq 107.49
$$

である。

規格化減衰時間は

$$
\tau_j D_i/L^2=\frac{1}{\left(\frac{1100 \pi^{2}}{101}\right)}\simeq 0.0093031
$$

である。

モード振幅の残存割合は

$$
n_j(t_{\rm obs})/n_j(0)=e^{- \left(\frac{1100 \pi^{2}}{101}\right) \left(\frac{1}{100}\right)}\simeq 0.34133
$$

である。

電子圧の効果で両極性拡散は速くなるが、単独の電子拡散には及ばない。
:::

:::example[個別の拡散係数が等しい場合]
両極性電場が要るか計算する。

入力値を次のように置く。

$$
\begin{aligned}
r&=10 \\t&=\frac{1}{10} \\j&=1 \\a&=\frac{1}{10}
\end{aligned}
$$

これらを順に代入する。

電子とイオンの拡散比は

$$
D_e/D_i=\left(10\right) \left(\frac{1}{10}\right)=1
$$

である。

両極性拡散係数の倍率は

$$
D_a/D_i=\frac{\left(10\right) + \left(1\right)}{\left(10\right) + 1}=1
$$

である。

密度勾配に対する電場係数は

$$
eE/(T_i\partial_x\ln n)=\frac{1 - \left(1\right)}{\left(10\right) + 1}=0
$$

である。

規格化した密度減衰率は

$$
\gamma_j L^2/D_i=\pi^{2} \left(1\right)^{3}\simeq 9.8696
$$

である。

規格化減衰時間は

$$
\tau_j D_i/L^2=\frac{1}{\left(\pi^{2}\right)}\simeq 0.10132
$$

である。

モード振幅の残存割合は

$$
n_j(t_{\rm obs})/n_j(0)=e^{- \left(\frac{1}{10}\right) \left(\pi^{2}\right)}\simeq 0.37271
$$

である。

拡散だけで流束が一致するため体積に補償電場は生じない。
:::

:::example[高次モードの平滑化]
最初の例と同じ媒質で細かな空間構造を追う。

入力値を次のように置く。

$$
\begin{aligned}
r&=100 \\t&=10 \\j&=3 \\a&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

電子とイオンの拡散比は

$$
D_e/D_i=\left(100\right) \left(10\right)=1000
$$

である。

両極性拡散係数の倍率は

$$
D_a/D_i=\frac{\left(1000\right) + \left(100\right)}{\left(100\right) + 1}=\frac{1100}{101}
$$

である。

密度勾配に対する電場係数は

$$
eE/(T_i\partial_x\ln n)=\frac{1 - \left(1000\right)}{\left(100\right) + 1}=- \frac{999}{101}
$$

である。

規格化した密度減衰率は

$$
\gamma_j L^2/D_i=\pi^{2} \left(3\right)^{2} \left(\frac{1100}{101}\right)\simeq 967.42
$$

である。

規格化減衰時間は

$$
\tau_j D_i/L^2=\frac{1}{\left(\frac{9900 \pi^{2}}{101}\right)}\simeq 0.0010337
$$

である。

モード振幅の残存割合は

$$
n_j(t_{\rm obs})/n_j(0)=e^{- \left(\frac{1}{100}\right) \left(\frac{9900 \pi^{2}}{101}\right)}\simeq 6.2887\times 10^{-5}
$$

である。

三倍の波数に対して減衰率は九倍になるので、遅い時刻には長い波長が残る。
:::

:::example[イオンの拡散が速い側]
イオンを止める向きの電場になるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
r&=2 \\t&=\frac{1}{10} \\j&=1 \\a&=\frac{1}{20}
\end{aligned}
$$

これらを順に代入する。

電子とイオンの拡散比は

$$
D_e/D_i=\left(2\right) \left(\frac{1}{10}\right)=\frac{1}{5}
$$

である。

両極性拡散係数の倍率は

$$
D_a/D_i=\frac{\left(2\right) + \left(\frac{1}{5}\right)}{\left(2\right) + 1}=\frac{11}{15}
$$

である。

密度勾配に対する電場係数は

$$
eE/(T_i\partial_x\ln n)=\frac{1 - \left(\frac{1}{5}\right)}{\left(2\right) + 1}=\frac{4}{15}
$$

である。

規格化した密度減衰率は

$$
\gamma_j L^2/D_i=\pi^{2} \left(1\right)^{2} \left(\frac{11}{15}\right)\simeq 7.2377
$$

である。

規格化減衰時間は

$$
\tau_j D_i/L^2=\frac{1}{\left(\frac{11 \pi^{2}}{15}\right)}\simeq 0.13817
$$

である。

モード振幅の残存割合は

$$
n_j(t_{\rm obs})/n_j(0)=e^{- \left(\frac{11 \pi^{2}}{15}\right) \left(\frac{1}{20}\right)}\simeq 0.69636
$$

である。

電場係数は正で、密度が外側へ減る場合に電場は内側を向く。
:::


## 磁場を横切る古典拡散

一様磁場の横方向に一定の力 $\vb F$ があるとき、衝突摩擦つきの定常式は $(\nu I-\Omega C)\vb u=\vb F/m$ である。$C\vb u=\vb u\times\vb b$ と定義すれば、横平面で $C^2=-I$ なので逆作用素は $(\nu I+\Omega C)/(\nu^2+\Omega^2)$ である。$\vb F=q\vb E-T\nabla\ln n$ とすれば、勾配方向の拡散係数は
$$
D_\perp=\frac{T}{m}\frac{\nu}{\nu^2+\Omega^2},\qquad
D_\parallel=\frac{T}{m\nu}.
\label{eq-plasma-diffusion-magnetic}
$$
直交する粒子束も残る。$\nu\ll\abs\Omega$ なら $D_\perp\simeq\rho_t^2\nu$ である。これは衝突で案内中心が典型的に $\rho_t$ だけずれるという描像と合う。

:::example[抑制率]
$\Omega/\nu=10$ なら $D_\perp/D_\parallel=1/101$ である。衝突を増やしたとき、強磁場側では横拡散が増え、平行拡散は減る。
:::

### 完全電離系で拡散を担う摩擦

一価で準中性な電子とイオンを考え、両種の密度と全圧が $x$ だけに依存し、$\vb B=B\hat{\vb z}$ はこの局所計算で一様とする。慣性と粘性を省き、弱い摩擦による横流を考える。全圧を $p$ とすれば $x$ 方向の力の和から最低次で $j_y=p'/B$ である。電子の $y$ 方向の式で $E_y=0$ とし、電子への摩擦を $R_{e,y}=m_e\nu_{ei}j_y/e$ とすれば
$$
enB u_{e,x}+R_{e,y}=0,qquad
\Gamma_{e,x}=-\frac{m_e\nu_{ei}}{e^2 B^2}p' .
$$
イオンの式では電荷と摩擦の符号がともに反転し、同じ $x$ 方向の束を得る。$p=n(T_e+T_i)$ で温度が一定なら $D_\perp=m_e\nu_{ei}(T_e+T_i)/(e^2B^2)$ である。これは $E_y=0$ で、膨張が緩やかな局所の二流体模型である。電流は Ampère 則により $B_z'=-\mu_0 j_y$ の勾配を伴う。例えば幅にわたる圧力変化が磁気圧に比べて小さければ、$B_z$ の変化も小さく、係数では代表値 $B$ を使える。磁場の時間発展まで求める場合には誘導式も要る。

同種粒子同士だけの衝突は、その種の全運動量を保存するため、この局所摩擦には寄与しない。一様磁場の中で同じ $q/m$ を持つ二粒子が同じ場所近くで衝突する場合、案内中心の変化の和は $(\Delta\vb v_1+\Delta\vb v_2)\times\vb b/\Omega=0$ である。このため、標識した一粒子の拡散と、種全体の最低次の拡散束は区別する。

### 速度相関による導出

一様磁場中の速度が、率 $\nu$ の Poisson 過程で独立な等方熱速度に取り直される模型を考える。$t$ まで衝突がない確率は $e^{-\nu t}$ であり、この間は旋回するから $\langle v_x(t)v_x(0)\rangle=(T/m)e^{-\nu t}\cos\Omega t$ となる。衝突がある場合には独立性で相関がゼロになる。変位の分散を二重時間積分し長時間でその $2t$ 当たりの極限をとると
$$
D_\perp=\int_0^\infty \langle v_x(t)v_x(0)\rangle\dd t
=\frac{T}{m}\frac{\nu}{\nu^2+\Omega^2}
$$
である。実際、定常な相関 $C(\tau)$ なら変位分散は $2\int_0^t(t-\tau)C(\tau)\dd\tau$ であり、指数減衰で積分可能なので極限はこの値になる。この模型では流束の計算と確率過程の計算が一致する。完全電離系の種間・種内衝突には総運動量の制約もあるため、独立に速度を取り直す模型の係数をそのまま用いない。

:::example[低磁化の極限]
$\abs\Omega\ll\nu$ では $\cos\Omega t$ が相関時間内にほぼ $1$ なので $D_\perp\simeq T/(m\nu)$ である。逆に強磁化では正負の速度相関が打ち消し合い、$D_\perp$ が小さくなる。
:::

### 磁化ブラウン運動の相関から有限時間拡散へ

等方な散乱を一定の $\nu$ で近似し、各直交成分の熱速度の二乗平均を $v_t^2=T/m$ とする。Langevin 方程式で速度の条件付き平均は $e^{-\nu t}$ を掛けた旋回であるため、$C_{xx}(t)=v_t^2 e^{-\nu t}\cos\Omega t$ となる。変位分散は速度の二重積分であり、定常性から $\langle(\Delta x)^2\rangle=2\int_0^t(t-s)C_{xx}(s)\dd s$ である。よって瞬時の半分散勾配 $D(t)=\int_0^t C_{xx}(s)\dd s$ は初め $v_t^2t$ で増え、長時間で $v_t^2\nu/(\nu^2+\Omega^2)$ へ近づく。旋回によって短い時間の相関は負になることもあるので、有限時間の測定を一律に拡散近似へ置き換えない。平行成分には旋回がないため $D_\parallel=v_t^2/\nu$ である。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $r$ | 磁化比 $\Omega/\nu$ | 1 |
| $s$ | 規格化時間 $\nu t$ | 1 |
| $v_t$ | 一方向熱速さ | $\mathrm{m\,s^{-1}}$ |
| $\nu$ | 速度相関の緩和率 | $\mathrm{s^{-1}}$ |

平行拡散係数を $D_\parallel$ とする。磁場平行の相関は指数だけである。

$$
D_\parallel=\frac{v_{t}^{2}}{\nu}
$$

横の長時間拡散係数を $D_\perp$ とする。余弦を伴う Laplace 積分から求める。

$$
D_\perp=\frac{D_{\parallel}}{r^{2} + 1}
$$

規格化した速度相関を $C_{xx}(t)/v_t^2$ とする。符号は旋回位相にも依存する。

$$
C_{xx}(t)/v_t^2=e^{- s} \cos{\left(r s \right)}
$$

横変位の瞬時半分散勾配を $D(t)$ とする。有限時間まで相関を積分する。

$$
D(t)=D_{\perp} \left(- \left(- r \sin{\left(r s \right)} + \cos{\left(r s \right)}\right) e^{- s} + 1\right)
$$

長時間極限に対する瞬時値を $D(t)/D_\perp$ とする。短い時間には一より大きい場合も負の場合もある。

$$
D(t)/D_\perp=\frac{D(t)}{D_{\perp}}
$$

平行の相関長の尺度を $\ell_\parallel$ とする。緩和時間に熱速さを掛ける。

$$
\ell_\parallel=\frac{v_{t}}{\nu}
$$

:::example[無磁場の短時間観測]
長時間拡散を短い窓に適用できるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
r&=0 \\s&=\frac{1}{100} \\v_t&=100\,\mathrm{m\,s^{-1}} \\\nu&=10\,\mathrm{s^{-1}}
\end{aligned}
$$

これらを順に代入する。

平行拡散係数は

$$
D_\parallel=\frac{\left(100\right)^{2}}{\left(10\right)}=1000\,\mathrm{m^2\,s^{-1}}
$$

である。

横の長時間拡散係数は

$$
D_\perp=\frac{\left(1000\right)}{\left(0\right)^{2} + 1}=1000\,\mathrm{m^2\,s^{-1}}
$$

である。

規格化した速度相関は

$$
C_{xx}(t)/v_t^2=e^{- \left(\frac{1}{100}\right)} \cos{\left(\left(0\right) \left(\frac{1}{100}\right) \right)}\simeq 0.99005
$$

である。

横変位の瞬時半分散勾配は

$$
D(t)=\left(1000\right) \left(- \left(- \left(0\right) \sin{\left(\left(0\right) \left(\frac{1}{100}\right) \right)} + \cos{\left(\left(0\right) \left(\frac{1}{100}\right) \right)}\right) e^{- \left(\frac{1}{100}\right)} + 1\right)\simeq 9.9502\,\mathrm{m^2\,s^{-1}}
$$

である。

長時間極限に対する瞬時値は

$$
D(t)/D_\perp=\frac{\left(1000 - \frac{1000}{e^{\frac{1}{100}}}\right)}{\left(1000\right)}\simeq 0.0099502
$$

である。

平行の相関長の尺度は

$$
\ell_\parallel=\frac{\left(100\right)}{\left(10\right)}=10\,\mathrm{m}
$$

である。

瞬時勾配は長時間値の約一パーセントしかなく、初期は弾道的である。
:::

:::example[磁化が始まる領域]
旋回時間と衝突時間が同程度のときの抑制を求める。

入力値を次のように置く。

$$
\begin{aligned}
r&=1 \\s&=1 \\v_t&=100\,\mathrm{m\,s^{-1}} \\\nu&=10\,\mathrm{s^{-1}}
\end{aligned}
$$

これらを順に代入する。

平行拡散係数は

$$
D_\parallel=\frac{\left(100\right)^{2}}{\left(10\right)}=1000\,\mathrm{m^2\,s^{-1}}
$$

である。

横の長時間拡散係数は

$$
D_\perp=\frac{\left(1000\right)}{\left(1\right)^{2} + 1}=500\,\mathrm{m^2\,s^{-1}}
$$

である。

規格化した速度相関は

$$
C_{xx}(t)/v_t^2=e^{- \left(1\right)} \cos{\left(\left(1\right)^{2} \right)}\simeq 0.19877
$$

である。

横変位の瞬時半分散勾配は

$$
D(t)=\left(500\right) \left(- \left(- \left(1\right) \sin{\left(\left(1\right)^{2} \right)} + \cos{\left(\left(1\right)^{2} \right)}\right) e^{- \left(1\right)} + 1\right)\simeq 555.4\,\mathrm{m^2\,s^{-1}}
$$

である。

長時間極限に対する瞬時値は

$$
D(t)/D_\perp=\frac{\left(- \frac{500 \left(- \sin{\left(1 \right)} + \cos{\left(1 \right)}\right)}{e} + 500\right)}{\left(500\right)}\simeq 1.1108
$$

である。

平行の相関長の尺度は

$$
\ell_\parallel=\frac{\left(100\right)}{\left(10\right)}=10\,\mathrm{m}
$$

である。

横拡散は平行の半分となるが、有限時間値がそのまま半分になるわけではない。
:::

:::example[強磁化で四分の一旋回]
短い時間で大きな瞬時勾配が現れるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
r&=10 \\s&=\frac{\pi}{20} \\v_t&=100\,\mathrm{m\,s^{-1}} \\\nu&=10\,\mathrm{s^{-1}}
\end{aligned}
$$

これらを順に代入する。

平行拡散係数は

$$
D_\parallel=\frac{\left(100\right)^{2}}{\left(10\right)}=1000\,\mathrm{m^2\,s^{-1}}
$$

である。

横の長時間拡散係数は

$$
D_\perp=\frac{\left(1000\right)}{\left(10\right)^{2} + 1}=\frac{1000}{101}\,\mathrm{m^2\,s^{-1}}
$$

である。

規格化した速度相関は

$$
C_{xx}(t)/v_t^2=e^{- \left(\frac{\pi}{20}\right)} \cos{\left(\left(10\right) \left(\frac{\pi}{20}\right) \right)}=0
$$

である。

横変位の瞬時半分散勾配は

$$
D(t)=\left(\frac{1000}{101}\right) \left(- \left(- \left(10\right) \sin{\left(\left(10\right) \left(\frac{\pi}{20}\right) \right)} + \cos{\left(\left(10\right) \left(\frac{\pi}{20}\right) \right)}\right) e^{- \left(\frac{\pi}{20}\right)} + 1\right)\simeq 94.518\,\mathrm{m^2\,s^{-1}}
$$

である。

長時間極限に対する瞬時値は

$$
D(t)/D_\perp=\frac{\left(\frac{1000}{101} + \frac{10000}{101 e^{\frac{\pi}{20}}}\right)}{\left(\frac{1000}{101}\right)}\simeq 9.5464
$$

である。

平行の相関長の尺度は

$$
\ell_\parallel=\frac{\left(100\right)}{\left(10\right)}=10\,\mathrm{m}
$$

である。

長時間値に比べて一時的に勾配が大きい。時間平均を長くすると正負の寄与が相殺する。
:::

:::example[強磁化で四分の三旋回]
変位分散の瞬時傾きの符号を求める。

入力値を次のように置く。

$$
\begin{aligned}
r&=10 \\s&=\frac{3 \pi}{20} \\v_t&=100\,\mathrm{m\,s^{-1}} \\\nu&=10\,\mathrm{s^{-1}}
\end{aligned}
$$

これらを順に代入する。

平行拡散係数は

$$
D_\parallel=\frac{\left(100\right)^{2}}{\left(10\right)}=1000\,\mathrm{m^2\,s^{-1}}
$$

である。

横の長時間拡散係数は

$$
D_\perp=\frac{\left(1000\right)}{\left(10\right)^{2} + 1}=\frac{1000}{101}\,\mathrm{m^2\,s^{-1}}
$$

である。

規格化した速度相関は

$$
C_{xx}(t)/v_t^2=e^{- \left(\frac{3 \pi}{20}\right)} \cos{\left(\left(10\right) \left(\frac{3 \pi}{20}\right) \right)}=0
$$

である。

横変位の瞬時半分散勾配は

$$
D(t)=\left(\frac{1000}{101}\right) \left(- \left(- \left(10\right) \sin{\left(\left(10\right) \left(\frac{3 \pi}{20}\right) \right)} + \cos{\left(\left(10\right) \left(\frac{3 \pi}{20}\right) \right)}\right) e^{- \left(\frac{3 \pi}{20}\right)} + 1\right)\simeq -51.904\,\mathrm{m^2\,s^{-1}}
$$

である。

長時間極限に対する瞬時値は

$$
D(t)/D_\perp=\frac{\left(\frac{1000}{101} - \frac{10000}{101 e^{\frac{3 \pi}{20}}}\right)}{\left(\frac{1000}{101}\right)}\simeq -5.2423
$$

である。

平行の相関長の尺度は

$$
\ell_\parallel=\frac{\left(100\right)}{\left(10\right)}=10\,\mathrm{m}
$$

である。

軌道が出発点へ近づく期間には変位分散が減る。長時間拡散係数は依然として正である。
:::


## Bohm 拡散と新古典輸送の概観

磁場を横切る粒子束の一部は、古典的な衝突で説明できない。例として乱れた電位が $E_\perp\sim T_e/(e\ell)$ を作ると、電場ドリフトは $v_E\sim T_e/(eB\ell)$ となる。相関時間を $\ell/v_E$ と見積もれば $D\sim v_E^2\ell/v_E=T_e/(eB)$ である。経験的な :term[Bohm 拡散]{en="Bohm diffusion" yomi="ぼーむかくさん"} は $D_B=T_e/(16 eB)$ である。相関時間の推定から係数 $1/16$ は出ない。

### トーラスでの軌道幅

軸対称トーラスでは、捕捉粒子が往復する間に径方向へ移動する。この幅を $\Delta$ とする。大アスペクト比 $\epsilon=r/R\ll1$ で磁場の相対変化が $\order{\epsilon}$ なら、ミラー条件から捕捉領域は $\abs{v_\parallel}/v\lesssim\sqrt\epsilon$ である。安全係数を $q_{\mathrm{sf}}$ とすると、往復時間は $q_{\mathrm{sf}}R/(v\sqrt\epsilon)$ 程度である。径方向ドリフトの代表速度 $v\rho/R$ を掛けて $\Delta\sim q_{\mathrm{sf}}\rho/\sqrt\epsilon$ となる。ピッチ角拡散が幅 $\sqrt\epsilon$ に達する時間は $\epsilon/\nu$ なので、捕捉粒子割合 $\sqrt\epsilon$ を掛けて $D_{\rm banana}\sim q_{\mathrm{sf}}^2\rho^2\nu/\epsilon^{3/2}$ となる。これは衝突間に往復できる領域のオーダー評価である。このように、トーラス内の軌道の効果を衝突輸送に含める理論を :term[新古典輸送]{en="neoclassical transport" yomi="しんこてんゆそう"} という。対照となる局所的な Coulomb 輸送の運動論は :cite[goldston-rutherford] を参照する。幾何学と衝突領域をまたぐ係数を得るには、トーラスのドリフト運動論を解く必要がある。

:::example[磁場の倍率で区別する]
他の条件が固定なら $B$ を二倍にすると Bohm 係数は半分、強磁場極限の古典係数は四分の一となる。新古典係数を比較する場合は磁気面の形と安全係数も保つ必要がある。
:::

:::example[往復可能かの判定]
$\epsilon=1/9$ として、ピッチ角の散乱時間 $\epsilon/\nu$ と往復時間 $q_{\rm sf}R/(v\sqrt\epsilon)$ の比をとると $v/(27\nu q_{\rm sf}R)$ となる。$\nu q_{\rm sf}R/v=10^{-3}$ なら比は約 $37$ であり、一度の離脱までに数十回の往復が可能という見積もりである。
:::

### 粒子の軌道と衝突率の領域分け

バナナの軌道幅を使った結果は、衝突がなければ同じ磁気面の近傍を往復する粒子に、衝突が小さな不可逆変位を与える、という近似である。捕捉域内の散乱時間が往復時間より短い場合は、粒子がバナナ全体をたどれない。強い衝突極限では局所の平行流が磁気面上の粒子収支を支える。この間の過渡領域を含めた輸送率は、衝突演算子と平行移流の境界条件をそろえて計算する。例えば、衝突率を一度に変えたときに横拡散が単調増加すると、バナナの低衝突側の比例式だけから結論することはできない。

:::example[捕捉割合を含める理由]
同じ粒子速度で捕捉粒子の径方向係数を $D_{\rm tr}\sim \Delta^2/(\epsilon/\nu)$ と見積もる。捕捉割合 $f_{\rm tr}\sim\sqrt\epsilon$ を掛けた集団の係数は $f_{\rm tr}D_{\rm tr}\sim q_{\rm sf}^2\rho^2\nu\epsilon^{-3/2}$ である。割合を掛けない計算では $\epsilon^{-2}$ になり、低アスペクト比への依存性が異なる。
:::

### 通過粒子が支配する衝突領域

バナナ領域から衝突を増やすと、捕捉粒子が一往復する前にピッチ角を変える。大アスペクト比の軸対称トーラスでは、$\nu qR/v$ が $\epsilon^{3/2}$ と $1$ の間にある遷移域を :term[プラトー領域]{en="plateau regime" yomi="ぷらとーりょういき"}、$\nu qR/v\gg1$ を :term[Pfirsch–Schlüter 領域]{en="Pfirsch–Schlüter regime" yomi="ぷふぃるしゅしゅるーたーりょういき"} と呼ぶ。前者の先頭の輸送係数は衝突率にほぼ依存せず、後者では再び衝突率に比例する。

後者の比例を拡散型の模型で見積もる。衝突間に磁気面を一周できないため、通過粒子は磁力線方向に $D_\parallel\sim v^2/\nu$ で拡散するとみなす。半径方向ドリフトの符号を決めるポロイダル位置が相関を失うまでには $(qR)^2/D_\parallel$ ほどかかる。ドリフトの大きさ $v\rho/R$ の二乗にその時間を掛けると $D_{\rm PS}\sim q^2\rho^2\nu$ になる。これは平均値がゼロのポロイダル高調波で偏りを記述する見積もりであり、輸送係数の数因子には速度平均と幾何学が要る。

プラトー領域に衝突率が現れない理由は、通過粒子の角度平均を簡単な模型で見ると分かる。径方向ドリフトをポロイダル角の一つの高調波で表し、$u=v_\parallel$ の粒子の相関を $e^{-\nu t}\cos(ut/(qR))$ とする。$u$ の規格化分布を $g(u)$ とすれば、相関の時間積分の速度平均は
$$
\int_{\R} g(u)\frac{\nu}{\nu^2+(u/(qR))^2}\,\dd u
=qR\int_{\R}\frac{g(qR\nu y)}{1+y^2}\,\dd y .
$$
$g$ が $0$ で連続かつ全域で有界であれば $\nu qR/v\ll1$ で $\pi qR g(0)$ に近づく。$g(0)\sim 1/v$ とドリフト速度の二乗を使うと、$D\sim qv\rho^2/R$ という衝突率によらない尺度になる。この相関模型は捕捉粒子の離脱や境界層の構造を持たない。従ってこの計算からプラトー領域の下限や係数の数因子は決められない。

### 捕捉粒子の幅と輸送領域を見分ける

大アスペクト比 $\epsilon=r/R\ll1$ の円形磁気面で、磁場の相対変動は $\order{\epsilon}$ である。捕捉境界は $v_\parallel^2/v^2=\order{\epsilon}$ であり、等方分布の捕捉割合は $\order{\sqrt\epsilon}$ となる。ピッチ角方向の拡散でこの幅を横切る時間は $\tau_{\rm det}\sim\epsilon/\nu$、一往復の時間は $\tau_b\sim qR/(\sqrt\epsilon v)$ である。両者が分離する条件は $\nu_*=\nu qR/(\epsilon^{3/2}v)\ll1$ である。捕捉粒子の半径方向の幅を $\Delta_b\sim q\rho/\sqrt\epsilon$ とし、衝突で軌道が更新されると考えれば、半径方向の係数は捕捉割合と $\Delta_b^2/\tau_{\rm det}$ の積、すなわち $D_b\sim q^2\rho^2\nu/\epsilon^{3/2}$ となる。係数と磁気面上の平均は装置の形状に依存する。以下はオーダーの比較であり、$\nu_*\gtrsim1$ では $D_b$ を輸送係数に採用しない。

:::remark[新古典輸送の比較条件]
$\rho$ を比較の単位に取る場合は温度と粒子種を指定する。電子とイオンでは同じ磁場でも $\rho$ が異なるため、無次元の倍率が同じでも物理単位の拡散係数は変わる。
:::

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $\epsilon$ | 逆アスペクト比 | 1 |
| $q$ | 安全係数 | 1 |
| $C$ | 衝突の規格化値 $\nu R/v$ | 1 |
| $r_\rho$ | 軌道半径の規格化値 $\rho/R$ | 1 |

捕捉割合の尺度を $f_t$ とする。$v_\parallel/v$ が $\order{\sqrt\epsilon}$ の区間を考える。

$$
f_t=\sqrt{\epsilon}
$$

往復時間の尺度を $\tau_b v/R$ とする。磁力線方向の距離は $qR$ 程度である。

$$
\tau_b v/R=\frac{q}{\sqrt{\epsilon}}
$$

捕捉域を出る時間の尺度を $\tau_{\rm det}v/R$ とする。拡散で角幅の二乗だけ時間が掛かる。

$$
\tau_{\rm det}v/R=\frac{\epsilon}{C}
$$

一往復に対する衝突の尺度を $\nu_*$ とする。一より十分小さい領域でバナナ軌道が成立する。

$$
\nu_*=\frac{\tau_b v/R}{\tau_{\rm det}v/R}
$$

バナナ幅の規格化値を $\Delta_b/R$ とする。$\Delta_b/r$ も小さい必要がある。

$$
\Delta_b/R=\frac{q r_{\rho}}{\sqrt{\epsilon}}
$$

古典拡散に対する単純模型の倍率を $D_b/(\rho^2\nu)$ とする。捕捉割合と一回の半径方向変位を含める。

$$
D_b/(\rho^2\nu)=\frac{(\Delta_b/R)^2 f_{t}}{\epsilon r_{\rho}^{2}}
$$

:::example[十分に低い衝突率]
往復の後に散乱されるか調べ、拡散の倍率を求める。

入力値を次のように置く。

$$
\begin{aligned}
\epsilon&=\frac{1}{10} \\q&=2 \\C&=\frac{1}{10000} \\r_\rho&=\frac{1}{10000}
\end{aligned}
$$

これらを順に代入する。

捕捉割合の尺度は

$$
f_t=\sqrt{\left(\frac{1}{10}\right)}\simeq 0.31623
$$

である。

往復時間の尺度は

$$
\tau_b v/R=\frac{\left(2\right)}{\sqrt{\left(\frac{1}{10}\right)}}\simeq 6.3246
$$

である。

捕捉域を出る時間の尺度は

$$
\tau_{\rm det}v/R=\frac{\left(\frac{1}{10}\right)}{\left(\frac{1}{10000}\right)}=1000
$$

である。

一往復に対する衝突の尺度は

$$
\nu_*=\frac{\left(2 \sqrt{10}\right)}{\left(1000\right)}\simeq 0.0063246
$$

である。

バナナ幅の規格化値は

$$
\Delta_b/R=\frac{\left(2\right) \left(\frac{1}{10000}\right)}{\sqrt{\left(\frac{1}{10}\right)}}\simeq 0.00063246
$$

である。

古典拡散に対する単純模型の倍率は

$$
D_b/(\rho^2\nu)=\frac{\left(\frac{\sqrt{10}}{10}\right) \left(\frac{\sqrt{10}}{5000}\right)^{2}}{\left(\frac{1}{10000}\right)^{2} \left(\frac{1}{10}\right)}\simeq 126.49
$$

である。

$\nu_*$ は小さく、$\Delta_b/r$ にも尺度の分離があるので、低衝突模型の前提がそろう。
:::

:::example[同じ磁気面で衝突を増やす]
バナナ係数の形式値が変わらなくても条件が破れるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
\epsilon&=\frac{1}{10} \\q&=2 \\C&=\frac{1}{10} \\r_\rho&=\frac{1}{10000}
\end{aligned}
$$

これらを順に代入する。

捕捉割合の尺度は

$$
f_t=\sqrt{\left(\frac{1}{10}\right)}\simeq 0.31623
$$

である。

往復時間の尺度は

$$
\tau_b v/R=\frac{\left(2\right)}{\sqrt{\left(\frac{1}{10}\right)}}\simeq 6.3246
$$

である。

捕捉域を出る時間の尺度は

$$
\tau_{\rm det}v/R=1
$$

である。

一往復に対する衝突の尺度は

$$
\nu_*=\frac{\left(2 \sqrt{10}\right)}{\left(1\right)}\simeq 6.3246
$$

である。

バナナ幅の規格化値は

$$
\Delta_b/R=\frac{\left(2\right) \left(\frac{1}{10000}\right)}{\sqrt{\left(\frac{1}{10}\right)}}\simeq 0.00063246
$$

である。

古典拡散に対する単純模型の倍率は

$$
D_b/(\rho^2\nu)=\frac{\left(\frac{\sqrt{10}}{10}\right) \left(\frac{\sqrt{10}}{5000}\right)^{2}}{\left(\frac{1}{10000}\right)^{2} \left(\frac{1}{10}\right)}\simeq 126.49
$$

である。

$\nu_*$ は一より大きく、完全な往復を前提とした低衝突の係数は使えない。
:::

:::example[細い軌道をさらに細くする]
旋回半径とバナナ幅を比較する。

入力値を次のように置く。

$$
\begin{aligned}
\epsilon&=\frac{1}{10} \\q&=2 \\C&=\frac{1}{10000} \\r_\rho&=\frac{1}{100000}
\end{aligned}
$$

これらを順に代入する。

捕捉割合の尺度は

$$
f_t=\sqrt{\left(\frac{1}{10}\right)}\simeq 0.31623
$$

である。

往復時間の尺度は

$$
\tau_b v/R=\frac{\left(2\right)}{\sqrt{\left(\frac{1}{10}\right)}}\simeq 6.3246
$$

である。

捕捉域を出る時間の尺度は

$$
\tau_{\rm det}v/R=\frac{\left(\frac{1}{10}\right)}{\left(\frac{1}{10000}\right)}=1000
$$

である。

一往復に対する衝突の尺度は

$$
\nu_*=\frac{\left(2 \sqrt{10}\right)}{\left(1000\right)}\simeq 0.0063246
$$

である。

バナナ幅の規格化値は

$$
\Delta_b/R=\frac{\left(2\right) \left(\frac{1}{100000}\right)}{\sqrt{\left(\frac{1}{10}\right)}}\simeq 6.3246\times 10^{-5}
$$

である。

古典拡散に対する単純模型の倍率は

$$
D_b/(\rho^2\nu)=\frac{\left(\frac{\sqrt{10}}{10}\right) \left(\frac{\sqrt{10}}{50000}\right)^{2}}{\left(\frac{1}{100000}\right)^{2} \left(\frac{1}{10}\right)}\simeq 126.49
$$

である。

規格化した拡散の倍率は一定でも、実際の拡散は $\rho^2$ の減少で二桁小さくなる。
:::

:::example[軸付近の局所近似]
一往復と半径方向の局所性を別々に調べる。

入力値を次のように置く。

$$
\begin{aligned}
\epsilon&=\frac{1}{100} \\q&=3 \\C&=\frac{1}{100000} \\r_\rho&=\frac{1}{1000}
\end{aligned}
$$

これらを順に代入する。

捕捉割合の尺度は

$$
f_t=\sqrt{\left(\frac{1}{100}\right)}=\frac{1}{10}
$$

である。

往復時間の尺度は

$$
\tau_b v/R=\frac{\left(3\right)}{\sqrt{\left(\frac{1}{100}\right)}}=30
$$

である。

捕捉域を出る時間の尺度は

$$
\tau_{\rm det}v/R=\frac{\left(\frac{1}{100}\right)}{\left(\frac{1}{100000}\right)}=1000
$$

である。

一往復に対する衝突の尺度は

$$
\nu_*=\frac{\left(30\right)}{\left(1000\right)}=\frac{3}{100}
$$

である。

バナナ幅の規格化値は

$$
\Delta_b/R=\frac{\left(3\right) \left(\frac{1}{1000}\right)}{\sqrt{\left(\frac{1}{100}\right)}}=\frac{3}{100}
$$

である。

古典拡散に対する単純模型の倍率は

$$
D_b/(\rho^2\nu)=\frac{\left(\frac{1}{10}\right) \left(\frac{3}{100}\right)^{2}}{\left(\frac{1}{1000}\right)^{2} \left(\frac{1}{100}\right)}=9000
$$

である。

$\nu_*$ は小さいが幅は $r$ より大きい尺度になるので、局所的な磁気面上の輸送模型は破れる。
:::


## 演習問題

::::exercise{level=basic}
$b=b_0$ の偏向角を求めよ。

:::solution
$\tan(\chi/2)=1$ より $\chi=\pi/2$ である。
:::
::::

::::exercise{level=basic}
対数を固定して $u$ が三倍のとき $\nu_m$ を比較せよ。

:::solution
$u^{-3}$ のため $1/27$ 倍である。
:::
::::

::::exercise{level=basic}
$\Omega=3\nu$ のとき $D_\perp/D_\parallel$ を求めよ。

:::solution
$\nu^2/(\nu^2+\Omega^2)=1/(1+9)=1/10$ である。
:::
::::

::::exercise{level=basic}
電子温度が九倍、対数が一定なら抵抗率は何倍か。

:::solution
依存性は $T_e^{-3/2}$ なので $9^{-3/2}=1/27$ である。
:::
::::

::::exercise{level=advanced}
二種類の移動度が等しいとき $D_A$ を導け。

:::solution
両種の数束を足すと電場項は打ち消し合う。その半分が共通の数束であり $D_A=(D_i+D_e)/2$ となる。
:::
::::

::::exercise{level=advanced}
$T,m$ と非零の $\Omega$ を固定し、$\nu>0$ を変えたとき $D_\perp$ はどこで最大か。

:::solution
$\partial_\nu [\nu/(\nu^2+\Omega^2)]=(\Omega^2-\nu^2)/(\nu^2+\Omega^2)^2$ である。正から負へ変わる $\nu=\abs\Omega$ で最大、値は $T/(2m\abs\Omega)$ である。
:::
::::

::::exercise{level=advanced}
二体の質量が等しいとき、一回の偏向で失う入射エネルギー割合を求めよ。

:::solution
$m_r=m_a/2$ なので割合は $(1-\cos\chi)/2=\sin^2(\chi/2)$ である。正面の後方散乱では全エネルギーが移る。
:::
::::

::::exercise{level=advanced}
遮蔽上限が $b_0$ の十倍なら、厳密な断面積と対数近似との差を求めよ。

:::solution
厳密値から近似値を引くと $2\pi b_0^2(\ln101-2\ln10)=2\pi b_0^2\ln(101/100)$ である。対数近似はこの模型ではわずかに小さい。
:::
::::

## 確認テスト

::::quiz
遮蔽のない小角散乱を距離で積分するとどうなるか。

- [x] 遠方で対数発散する。
- [ ] 遠方で距離の二乗に比例する。
- [ ] 有限の値に収束する。

:::explanation
重み $1-\cos\chi$ が $2b_0^2/b^2$ に近づき、環の面積が $2\pi b\,\dd b$ だから、積分は $\int\dd b/b$ を含む。
:::
::::

::::quiz{answer="0.2" tol="1e-6"}
$\Omega/\nu=2$ のときの $D_\perp/D_\parallel$ を答えよ。

:::explanation
両係数の比は $1/(1+\Omega^2/\nu^2)=1/5$ である。
:::
::::

::::quiz
無電流の一次元両極性拡散で等しいものはどれか。

- [x] イオンと電子の粒子束
- [ ] イオンと電子の移動度
- [ ] イオンと電子の温度

:::explanation
$\vb j=e(\vb\Gamma_i-\vb\Gamma_e)=0$ としたため粒子束が等しい。他の係数は種によって異なってよい。
:::
::::

::::quiz{answer="0.5" tol="1e-6"}
$T_e/e=16\,\mathrm V$、$B=2\,\mathrm T$ の Bohm 係数を $\mathrm{m^2\,s^{-1}}$ 単位で答えよ。

:::explanation
経験式 $T_e/(16eB)$ に代入すると $16/(16\cdot2)=0.5$ である。
:::
::::

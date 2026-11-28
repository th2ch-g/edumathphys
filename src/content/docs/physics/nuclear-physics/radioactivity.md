---
title: 放射性崩壊と放射線計測
description: 崩壊系列の時間発展と、放射線が物質に残す信号から強度や年代を求める方法を扱う。
notation:
  - symbol: '\lambda'
    meaning: 崩壊定数
  - symbol: '\mathcal A'
    meaning: 放射能
  - symbol: 'S=-\dd E/\dd x'
    meaning: 阻止能
  - symbol: '\mu'
    meaning: 線減弱係数
---

## この章の目標

- 崩壊確率から放射能と半減期を導く。
- 親核と娘核の崩壊を連立方程式で解く。
- 物質内での荷電粒子、光子、中性子の相互作用を比較する。
- 計数補正、加速、放射年代測定の計算を行う。

## 崩壊法則と半減期

### 一定ハザードの模型

短時間 $\dd t$ で生き残った核一つが崩壊する条件付き確率を $\lambda\dd t+o(\dd t)$ とする。$\lambda>0$ は :term[崩壊定数]{en="decay constant" yomi="ほうかいていすう"} である。核が互いに影響せず、崩壊率が時間に依存しないなら、生存確率 $P(t)$ は $P'=-\lambda P$、$P(0)=1$ を満たす。従って核数の期待値と :term[放射能]{en="activity" yomi="ほうしゃのう"} は
$$
N(t)=N_0e^{-\lambda t},\qquad
\mathcal A(t)=\lambda N(t).
\label{eq-nuc-decay-law}
$$
単位 $\mathrm{Bq}$ は $\mathrm{s}^{-1}$ である。崩壊時刻の密度 $\lambda e^{-\lambda t}$ を使って部分積分すると平均寿命は $\tau=\lambda^{-1}$ となる。一方、:term[半減期]{en="half-life" yomi="はんげんき"} は $P(T_{1/2})=1/2$ の解 $\ln2/\lambda$ である。崩壊確率が指数則に従うとき、生存振幅を $\theta(t)e^{-iE_0t/\hbar-t/(2\tau)}$ と置く模型では、Fourier 積分が $[\frac1{2\tau}-i(E-E_0)/\hbar]^{-1}$ に比例する。その絶対値二乗は Lorentz 型であり、半値点は $|E-E_0|=\hbar/(2\tau)$ なので自然幅の半値全幅は $\Gamma=\hbar/\tau$ となる。

競合チャネルの崩壊定数を $\lambda_i$ とすると和 $\sum_i\lambda_i$ が全崩壊定数になる。チャネル $i$ が選ばれる確率は $\int_0^\infty\lambda_i e^{-\lambda t}\dd t=\lambda_i/\lambda$ である。時刻を固定した生存数は二項分布であり、平均のほかに分散 $N_0P(1-P)$ をもつ。

:::example[半減期と分岐]
半減期 $1\,\mathrm h$ で初期放射能 $80\,\mathrm{Bq}$ の線源は 3 時間後に $10\,\mathrm{Bq}$ となる。観測対象の分岐が $0.2$ なら、そのチャネルの瞬間の崩壊は毎秒平均 2 回である。
:::

:::example[生存数のばらつき]
半減期が一回過ぎると、独立な 100 核の生存数は平均 50、分散 25 である。標準偏差 5 があるため、個々の実験でちょうど 50 核が残るとは限らない。
:::

### 分岐比と待ち時間

独立に開いた崩壊路 $i$ の短時間の発生確率が $\lambda_i\dd t+o(\dd t)$ であるとする。同じ核で二つの崩壊は重ならないから、生存確率 $S$ は $\dot S=-\sum_i\lambda_i S$ を満たす。従って経路 $i$ が実現する確率は $\int_0^\infty\lambda_i e^{-\lambda t}\dd t=\lambda_i/\lambda$、ただし $\lambda=\sum_i\lambda_i$ である。部分半減期の逆数は加算できるが、各経路に割り当てられる核の数は分岐比で重み付けする。

:::example[分岐ごとの放射能]
全放射能 $A$ のうち一つの線に対応する分岐比を $b$ とする。

量の定義に従い、共通の単位で次の値を用いる。
$$
A=2400\,\mathrm{Bq},\qquad b=\frac{3}{25}
$$

この条件で求める量を ${A_i}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{A_i}&=A \cdot b\\&=\left(2400\right) \cdot \left(\frac{3}{25}\right)\\&=288\,\mathrm{Bq}.\end{aligned}
$$

毎秒の放出数であり、検出器の計数率にするには効率も掛ける。

分岐比の違う線を比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
b=\frac{3}{10}
$$

このときの値を ${A_i}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{A_i}_{\mathrm{new}}&=\left(2400\right) \cdot \left(\frac{3}{10}\right)\\&=720\,\mathrm{Bq},\\{A_i}_{\mathrm{new}}-{A_i}&=432\,\mathrm{Bq}.\end{aligned}
$$
:::

:::example[二経路の全半減期]
部分半減期を $t_1,t_2$ とおく。独立な率を足してから $\ln2$ を割る。

量の定義に従い、共通の単位で次の値を用いる。
$$
t_{1}=6\,\mathrm{h},\qquad t_{2}=12\,\mathrm{h}
$$

この条件で求める量を ${T}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{T}&=\frac{1}{1 \cdot \frac{1}{t_{2}} + 1 \cdot \frac{1}{t_{1}}}\\&=\frac{1}{\frac{1}{\left(6\right)} + \frac{1}{\left(12\right)}}\\&=4\,\mathrm{h}.\end{aligned}
$$

全半減期はどの部分半減期よりも短い。

第二の経路の率が増す場合を考える。そこで次の値に変更し、ほかの入力を固定する。
$$
t_{2}=3\,\mathrm{h}
$$

このときの値を ${T}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{T}_{\mathrm{new}}&=\frac{1}{\frac{1}{\left(6\right)} + \frac{1}{\left(3\right)}}\\&=2\,\mathrm{h},\\{T}_{\mathrm{new}}-{T}&=-2\,\mathrm{h}.\end{aligned}
$$
:::

:::example[観測時間内の減少]
観測窓 $t$ の前後で核数を比べる。半減期を $T$ とすると崩壊済みの割合は生存確率の余事象である。

量の定義に従い、共通の単位で次の値を用いる。
$$
t=1\,\mathrm{h},\qquad T=8\,\mathrm{h}
$$

この条件で求める量を ${p}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{p}&=1 - 2^{\frac{\left(-1\right) \cdot t}{T}}\\&=1 - 2^{- \frac{\left(1\right)}{\left(8\right)}}\\&\simeq 0.082996.\end{aligned}
$$

窓が半減期より短い場合でも、分率を時間に比例とみなせる範囲はさらに狭い。

窓を延長して比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
t=4\,\mathrm{h}
$$

このときの値を ${p}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{p}_{\mathrm{new}}&=1 - 2^{- \frac{\left(4\right)}{\left(8\right)}}\\&\simeq 0.292893,\\{p}_{\mathrm{new}}-{p}&\simeq 0.209897.\end{aligned}
$$
:::

:::example[初期核数が有限な場合]
独立な $N$ 個から始め、時刻 $t$ の一核の生存確率を $p$ とする。指示変数の分散 $p(1-p)$ を足す。

量の定義に従い、共通の単位で次の値を用いる。
$$
N=100,\qquad p=\frac{1}{2}
$$

この条件で求める量を ${\sigma}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\sigma}&=\sqrt{N \cdot p \cdot \left(1 - p\right)}\\&=\sqrt{\left(100\right) \cdot \left(\frac{1}{2}\right) \cdot \left(1 - \left(\frac{1}{2}\right)\right)}\\&=5.\end{aligned}
$$

母数 $N$ を増やすと絶対揺らぎは増えるが平均に対する相対揺らぎは減る。

試料量を増やす。そこで次の値に変更し、ほかの入力を固定する。
$$
N=400
$$

このときの値を ${\sigma}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\sigma}_{\mathrm{new}}&=\sqrt{\left(400\right) \cdot \left(\frac{1}{2}\right) \cdot \left(1 - \left(\frac{1}{2}\right)\right)}\\&=10,\\{\sigma}_{\mathrm{new}}-{\sigma}&=5.\end{aligned}
$$
:::

## 逐次崩壊と放射平衡

親核 $1$ が全崩壊の割合 $b$ で娘核 $2$ を作るとする。初期娘核数を $N_{20}$ とすると $\dot N_2+\lambda_2N_2=b\lambda_1N_{10}e^{-\lambda_1t}$ である。積分因子 $e^{\lambda_2t}$ を掛けて時間積分すると
$$
N_2(t)=N_{20}e^{-\lambda_2t}+
b\lambda_1N_{10}\frac{e^{-\lambda_1t}-e^{-\lambda_2t}}{\lambda_2-\lambda_1}.
\label{eq-nuc-two-member-decay}
$$
$\lambda_1=\lambda_2=\lambda$ では積分因子を掛けた生成率が一定になり、第 2 項は $b\lambda N_{10}t e^{-\lambda t}$ となる。長い系列では各段で $N_j=e^{-\lambda_jt}[N_j(0)+\int_0^t e^{\lambda_js}\sum_i b_{ij}\lambda_iN_i(s)\dd s]$ を適用する。安定な末端は $\lambda_j=0$ とすればよい。

$\lambda_2>\lambda_1$ で初期娘核の寄与が消えると $\mathcal A_2/\mathcal A_1\to b\lambda_2/(\lambda_2-\lambda_1)$ となる。寿命が十分に離れていれば比は $b$ に近づく。この状況を :term[永続平衡]{en="secular equilibrium" yomi="えいぞくへいこう"} という。比はほぼ一定だが、両方の放射能は親核の半減期で減る。比が一定で $b$ から有意に離れる場合は過渡平衡である。

:::example[娘核のピーク]
$\lambda_2=2\lambda_1$、$N_{20}=0$、$b=1$ なら、$N_2=N_{10}(e^{-\lambda_1 t}-e^{-2\lambda_1t})$ である。微分をゼロとすると $t=\ln2/\lambda_1$、そのとき $N_2=N_{10}/4$ である。そこでは生成と消滅が等しく $\mathcal A_2=\mathcal A_1$ となる。
:::

:::example[同じ寿命]
$\lambda_2=\lambda_1$ で初期娘核がなければ、$t=\tau$ で娘核数は $bN_{10}/e$ となる。混合物の全活動を一つの指数関数に当てはめると、この $t$ に比例する係数が残差を作る。
:::

### 一直線の系列の一般解

初期には親核だけがあり、核種 $i$ から $i+1$ への分岐が $b_i$、崩壊定数が $\lambda_i$ とする。$N_i$ の Laplace 変換 $\widetilde N_i(s)=\int_0^\infty e^{-st}N_i(t)\dd t$ を取る。部分積分で $\mathcal L(\dot N_i)=s\widetilde N_i-N_i(0)$ となるので
$$
\widetilde N_j(s)=\frac{N_1(0)\prod_{i=1}^{j-1}b_i\lambda_i}
{\prod_{i=1}^{j}(s+\lambda_i)}
$$
である。$j=1$ では空積を 1 とする。$\lambda_i$ が相異なるなら、分母の逆数の部分分数係数は $c_i=[\prod_{h\ne i}(\lambda_h-\lambda_i)]^{-1}$ である。実際、両辺に分母を掛けて $s=-\lambda_i$ を代入すると、その項だけが残る。次数 $j-1$ 以下の差が $j$ 点でゼロになるので等式が成立する。逆変換すると
$$
N_j(t)=N_1(0)\left(\prod_{i=1}^{j-1}b_i\lambda_i\right)
\sum_{i=1}^j\frac{e^{-\lambda_i t}}{\prod_{h\ne i}(\lambda_h-\lambda_i)}
$$
を得る。これが :term[Bateman の解]{en="Bateman solution" yomi="べいとまんのかい"} である。重根では $(s+\lambda)^{-m}$ の逆変換 $t^{m-1}e^{-\lambda t}/(m-1)!$ が現れる。近い二つの $\lambda_i$ に一般式を直接使うと大きな係数同士を引くことになるため、数値計算では極限形や連立方程式の積分が役立つ。

:::example[安定末端]
$\lambda_1=\lambda$、$\lambda_2=2\lambda$、$\lambda_3=0$、すべての分岐が 1 なら、親と娘の合計は $N_1(0)(2e^{-\lambda t}-e^{-2\lambda t})$ である。核数保存により $N_3=N_1(0)(1-e^{-\lambda t})^2$ となる。$t=\ln2/\lambda$ で末端には初期数の $1/4$ が達する。
:::

### 一定供給と放射平衡への接近

生成率 $P$ が一定なら $\dot N+\lambda N=P$ である。この式に $e^{\lambda t}$ を掛け、$0$ から $t$ まで積分すると $N=N(0)e^{-\lambda t}+(P/\lambda)(1-e^{-\lambda t})$ となる。したがって照射直後の値を知るには照射時間と輸送時間を別に数える。親核の供給が指数的に減る場合は $P(t)$ を積分の内側に入れれば逐次崩壊の式を回復する。

:::example[飽和放射能への接近]
初期の放射性核種はゼロで、毎秒 $P$ 個を一定に作る。半減期は $T$ である。

量の定義に従い、共通の単位で次の値を用いる。
$$
P=500\,\mathrm{s}^{-1},\qquad t=1\,\mathrm{h},\qquad T=1\,\mathrm{h}
$$

この条件で求める量を ${A}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{A}&=P \cdot \left(1 - 2^{\frac{\left(-1\right) \cdot t}{T}}\right)\\&=\frac{\left(500\right)}{2}\\&=250\,\mathrm{Bq}.\end{aligned}
$$

生成率と崩壊率が等しくなるので飽和時には放射能が $P$ に近づく。

照射時間を延ばす。そこで次の値に変更し、ほかの入力を固定する。
$$
t=4\,\mathrm{h}
$$

このときの値を ${A}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{A}_{\mathrm{new}}&=\left(500\right) \cdot \left(1 - 2^{- \frac{\left(4\right)}{\left(1\right)}}\right)\\&\simeq 468.75\,\mathrm{Bq},\\{A}_{\mathrm{new}}-{A}&\simeq 218.75\,\mathrm{Bq}.\end{aligned}
$$
:::

:::example[照射後の測定]
照射直後に $A_0$ で、輸送に $t$ がかかる。輸送中に新規生成はない。

量の定義に従い、共通の単位で次の値を用いる。
$$
A_{0}=1600\,\mathrm{Bq},\qquad t=3\,\mathrm{min},\qquad T=3\,\mathrm{min}
$$

この条件で求める量を ${A}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{A}&=2^{\frac{\left(-1\right) \cdot t}{T}} \cdot A_{0}\\&=\frac{\left(1600\right)}{2}\\&=800\,\mathrm{Bq}.\end{aligned}
$$

半減期の短い核では試料採取時刻の定義が結果に強く効く。

測定開始を遅らせる。そこで次の値に変更し、ほかの入力を固定する。
$$
t=6\,\mathrm{min}
$$

このときの値を ${A}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{A}_{\mathrm{new}}&=2^{- \frac{\left(6\right)}{\left(3\right)}} \cdot \left(1600\right)\\&=400\,\mathrm{Bq},\\{A}_{\mathrm{new}}-{A}&=-400\,\mathrm{Bq}.\end{aligned}
$$
:::

:::example[過渡平衡での比]
分岐比 $b$ の親から生じる娘について、十分後の比 $A_2/A_1$ を求める。$\lambda_2>\lambda_1$ とする。

量の定義に従い、共通の単位で次の値を用いる。
$$
b=\frac{4}{5},\qquad l_{1}=1\,\mathrm{h}^{-1},\qquad l_{2}=6\,\mathrm{h}^{-1}
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=\frac{b \cdot l_{2}}{- l_{1} + l_{2}}\\&=\frac{\left(6\right) \cdot \left(\frac{4}{5}\right)}{- \left(1\right) + \left(6\right)}\\&\simeq 0.96.\end{aligned}
$$

娘の崩壊がはるかに速い場合はこの比が $b$ に近づく。差が小さいほど平衡に達するまで長い。

親の寿命を短くして比を調べる。そこで次の値に変更し、ほかの入力を固定する。
$$
l_{1}=2\,\mathrm{h}^{-1}
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=\frac{\left(6\right) \cdot \left(\frac{4}{5}\right)}{- \left(2\right) + \left(6\right)}\\&\simeq 1.2,\\{R}_{\mathrm{new}}-{R}&\simeq 0.24.\end{aligned}
$$
:::

:::example[同じ崩壊定数の二段階]
初めに娘はなく、率は両段階で $\lambda$ とする。娘の密度に積分因子を使うと $N_2/N_0=\lambda t e^{-\lambda t}$ になる。

量の定義に従い、共通の単位で次の値を用いる。
$$
u=1
$$

この条件で求める量を ${f}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{f}&=u \cdot e^{- u}\\&=\left(1\right) \cdot e^{- \left(1\right)}\\&\simeq 0.367879.\end{aligned}
$$

ここで $u=\lambda t$ である。この分率は初期には増え、$u=1$ で $\dd f/\dd u=0$ となる。

最大を過ぎた時刻と比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
u=2
$$

このときの値を ${f}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{f}_{\mathrm{new}}&=\left(2\right) \cdot e^{- \left(2\right)}\\&\simeq 0.270671,\\{f}_{\mathrm{new}}-{f}&\simeq -0.0972089.\end{aligned}
$$
:::

## 荷電粒子と物質の相互作用（Bethe–Bloch の式）

### 衝突損失の起源

速度 $v$、電荷 $ze$ の重い入射粒子が電子の近くを通ると、Coulomb 力が電子に運動量を渡す。電子を最初静止とし、入射粒子の軌道を直線、衝突径数を $b$ とすると横方向の力の時間積分は
$$
\Delta p_\perp=\frac{ze^2}{4\pi\varepsilon_0}\int_{-\infty}^\infty \frac{b\,\dd t}{(b^2+v^2t^2)^{3/2}}
=\frac{2ze^2}{4\pi\varepsilon_0bv}.
$$
電子のエネルギー増分 $\Delta p_\perp^2/(2m_e)$ に、単位長さの通過で出会う電子数 $2\pi b n_e\dd b$ を掛けると、$z^2/v^2$ と $\int\dd b/b$ が現れる。下限には大きなエネルギー移行、上限には電子の束縛や物質の分極が関与する。

物質の電子密度 $n_e$、古典電子半径 $r_e=e^2/(4\pi\varepsilon_0m_ec^2)$、平均励起エネルギー $I$ を用いると、重い粒子の衝突 :term[阻止能]{en="stopping power" yomi="そしのう"} の代表式は
$$
S=4\pi n_e r_e^2m_ec^2\frac{z^2}{\beta^2}
\left[\frac12\ln\frac{2m_ec^2\beta^2\gamma^2 T_{\max}}{I^2}-\beta^2-\frac{\delta}{2}\right],
\quad
T_{\max}=\frac{2m_ec^2\beta^2\gamma^2}{1+2\gamma m_e/M+(m_e/M)^2}.
$$
ここで $\beta=v/c$、$\gamma=(1-\beta^2)^{-1/2}$、$M$ は入射粒子の質量、$\delta$ は密度効果補正である。$T_{\max}$ は重心系で運動量が反転する二体弾性衝突で、標的電子が受け取る最大の運動エネルギーである。式は速い重荷電粒子に対する衝突損失の平均である。低速では殻補正や有効電荷が重要になり、電子入射の場合は同質量の交換や制動放射を別に扱う。

古典的な衝突積分に $\Delta p_\perp$ を代入すると、$4\pi n_e r_e^2m_ec^2 z^2/\beta^2$ に $\ln(b_{\max}/b_{\min})$ を掛けた形になる。束縛の効果を入れるには、衝突時間 $b/v$ と電子の励起角振動数 $\omega_j$ を比べる。$b\omega_j/v\gg1$ の摂動では電子が断熱的に追従するため励起が抑えられる。速い重粒子の Born 近似では、各励起エネルギー $\hbar\omega_j$ が低運動量移行側の対数に入り、その重みは双極子振動子強度 $f_j$ である。$I$ は $\ln I=\sum_j f_j\ln(\hbar\omega_j)/\sum_j f_j$ で定義され、連続状態は和を積分に替える。媒質の分極は遠距離での電場を遮蔽し、高速側で $\delta$ の補正を生む。相対論的な全補正の導出には束縛電子の応答と自由電子散乱の接続が必要であるため、ここでは適用条件と衝突積分を示した。阻止能と飛程の標準的な扱いは :cite[krane-nuclear] を参照する。

### Born 近似での対数項

直線軌道で得た対数は、速い重粒子の非相対論的な Born 近似からも導ける。密度 $n_{\rm at}$ の原子に対し、$q$ を移行波数、$A_j(\vb q)=\langle j|\sum_{\ell}e^{i\vb q\cdot\vb r_\ell}|0\rangle$ を電子の遷移振幅とする。$g_e=e^2/(4\pi\varepsilon_0)$ と書く。Coulomb Fourier 成分 $4\pi g_e/q^2$ から Born 振幅は $2Mz g_e A_j/(\hbar^2 q^2)$ に比例する。始状態と終状態の波数を $k_i,k_f$ とすれば、断面積には流束の比 $k_f/k_i$ が入る。また $q^2=k_i^2+k_f^2-2k_ik_f\cos\theta$ なので、方位の積分を含めた立体角は $2\pi q\dd q/(k_i k_f)$ である。従って角度平均した $|A_j|^2$ に対し
$$
\dd\sigma_j=\frac{8\pi z^2g_e^2}{\hbar^2v^2}|A_j|^2\frac{\dd q}{q^3}
$$
を得る。小さい $q$ では遷移エネルギー $E_j=\hbar\omega_j>0$ に対し $|A_j|^2=q^2\hbar^2 f_j/(2m_eE_j)$ となる。これは双極子行列要素の方向平均の係数 $1/3$ と、$f_j=2m_eE_j|\langle j|\sum_\ell\vb r_\ell|0\rangle|^2/(3\hbar^2)$ を用いた結果である。最小波数は $k_i-k_f\simeq E_j/(\hbar v)$ なので、低波数側の積分は $E_j$ の対数を与える。

電子の全振動子強度の和は電子数 $Z_{\rm at}$ に等しい。座標だけに依存する束縛ポテンシャルでは、双極子の一成分 $D=\sum_\ell x_\ell$ の二重交換子 $[D,[H,D]]=Z_{\rm at}\hbar^2/m_e$ を計算し、固有状態を挿入すれば確かめられる。高波数では電子をほぼ自由とみなし、上限の波数を $2m_ev/\hbar$ とする。これらを対数の精度で合わせると、$n_e=Z_{\rm at}n_{\rm at}$ を用いて $S\simeq[4\pi n_e z^2 g_e^2/(m_ev^2)]\ln(2m_ev^2/I)$ となる。これは $\beta^2\ll1$、$M\gg m_e$、かつ電子の束縛運動より入射粒子が速い領域で、Bethe–Bloch 式の主項と一致する。

### 飛程

損失を平均値で置き換えると $-\dd E/\dd x=S(E)$ である。分離積分した $R(E_0)=\int_0^{E_0}\dd E/S(E)$ を連続減速近似の飛程と呼ぶ。実際の停止位置は衝突のばらつきと横散乱で分布する。遅くなる重粒子で単位長さ当たりの損失が増える領域は Bragg ピークにつながる。

:::example[電荷による比較]
速度と物質を固定し、$M\gg m_e\gamma$、完全電離を仮定すると $z=2$ と $z=1$ の阻止能の比は 4 である。運動エネルギーを固定すれば質量の違いで $\beta$ も変わるので、この比を直接適用できない。
:::

:::example[単純な飛程]
区間 $0<E<E_0$ の模型で $S=a/E$ とすれば $R=E_0^2/(2a)$ である。この式では入射エネルギーを 2 倍にすると飛程は 4 倍になる。模型の低速での発散は、実際の阻止能を停止点まで表す性質ではない。
:::

### 一回の運動学から最大移行を求める

電子の質量を $m_e$、入射粒子の質量を $M$ とする。電子が最初静止しているとき、二体の四元運動量の不変量をエネルギー二乗で $s=(M^2+m_e^2+2\gamma Mm_e)c^4$ と書く。重心系の運動量 $p_*$ は、二粒子の全エネルギーの和を $\sqrt s$ に等置して
$$
p_*^2c^2=\frac{[s-(M+m_e)^2c^4][s-(M-m_e)^2c^4]}{4s}
$$
と得られる。入射電子から出射電子への移行の四元運動量の不変量を $t$ と書く。重心系ではエネルギー移行がゼロなので $-t=2p_*^2c^2(1-\cos\theta_*)$、実験室では電子の質量殻条件から $-t=2m_ec^2 T_e$ である。最大は $\cos\theta_*=-1$ で $T_e=2p_*^2/m_e$ となる。$s$ を代入すると Bethe–Bloch 式の $T_{\max}$ を得る。

:::example[質量比の意味]
$M=m_e$ と置くと $T_{\max}=(\gamma-1)m_ec^2$ となり、全運動エネルギーを標的へ渡せる。ただし入射電子と標的電子は同一粒子であり、重粒子向けの阻止能の式をこのまま適用することはできない。
:::

### 阻止能を運動量移行で表示する

衝突パラメータ $b$ の直線軌道近似で、電子に加わる横向きの力は定数倍の $b/(b^2+v^2t^2)^{3/2}$ である。$u=vt/b$ と置いた時間積分は $2/(bv)$ を与える。移行エネルギーは運動量の二乗を $2m_e$ で割って求まり、環の面積 $2\pi b\dd b$ を掛けると $\dd b/b$ が残る。衝突の集団平均に出る対数はこの積分に由来する。遠い衝突で電子の束縛を無視すると積分が発散するので、束縛電子の応答を入れて上下限に物理的な意味を与える。

:::example[二体衝突の最大移行]
電子より十分重い入射粒子で質量比の補正を無視する。$\gamma^2\beta^2=\beta^2/(1-\beta^2)$ とおいて上限を得る。

量の定義に従い、共通の単位で次の値を用いる。
$$
m_{ec2}=\frac{511}{1000}\,\mathrm{MeV},\qquad \beta=\frac{2}{5}
$$

この条件で求める量を ${T_{\max}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{T_{\max}}&=\frac{2 \cdot \beta^{2} \cdot m_{ec2}}{1 - \beta^{2}}\\&=\frac{2 \cdot \left(\frac{2}{5}\right)^{2} \cdot \left(\frac{511}{1000}\right)}{1 - \left(\frac{2}{5}\right)^{2}}\\&\simeq 0.194667\,\mathrm{MeV}.\end{aligned}
$$

速度を上げると移行できる上限が増える。これが各衝突の平均移行そのものではない。

相対論的な速度域で比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
\beta=\frac{4}{5}
$$

このときの値を ${T_{\max}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{T_{\max}}_{\mathrm{new}}&=\frac{2 \cdot \left(\frac{4}{5}\right)^{2} \cdot \left(\frac{511}{1000}\right)}{1 - \left(\frac{4}{5}\right)^{2}}\\&\simeq 1.81689\,\mathrm{MeV},\\{T_{\max}}_{\mathrm{new}}-{T_{\max}}&\simeq 1.62222\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[密度を用いた阻止能の換算]
質量阻止能を $S_m$、物質の密度を $\rho$ とする。短い区間 $d$ では阻止能を一定に近似する。

量の定義に従い、共通の単位で次の値を用いる。
$$
\rho=\frac{6}{5}\,\mathrm{g\,cm}^{-3},\qquad d=\frac{1}{20}\,\mathrm{cm},\qquad S_{m}=4\,\mathrm{MeV\,cm^2/g}
$$

この条件で求める量を ${\Delta E}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\Delta E}&=S_{m} \cdot d \cdot \rho\\&=\left(4\right) \cdot \left(\frac{1}{20}\right) \cdot \left(\frac{6}{5}\right)\\&\simeq 0.24\,\mathrm{MeV}.\end{aligned}
$$

阻止能を質量当たりで表す単位では、厚みは面密度に換算して掛ける。

試料の厚みを二倍にする。そこで次の値に変更し、ほかの入力を固定する。
$$
d=\frac{1}{10}\,\mathrm{cm}
$$

このときの値を ${\Delta E}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\Delta E}_{\mathrm{new}}&=\left(4\right) \cdot \left(\frac{1}{10}\right) \cdot \left(\frac{6}{5}\right)\\&\simeq 0.48\,\mathrm{MeV},\\{\Delta E}_{\mathrm{new}}-{\Delta E}&\simeq 0.24\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[低速近似での飛程]
ある狭いエネルギー帯で $S(E)=a/E$ と近似する。$-\dd E/\dd x=S$ から二点間を積分する。

量の定義に従い、共通の単位で次の値を用いる。
$$
E_{0}=5\,\mathrm{MeV},\qquad E_{1}=4\,\mathrm{MeV},\qquad a=10\,\mathrm{MeV^2/cm}
$$

この条件で求める量を ${\Delta x}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\Delta x}&=\frac{E_{0}^{2} - E_{1}^{2}}{2 \cdot a}\\&=\frac{- \left(4\right)^{2} + \left(5\right)^{2}}{2 \cdot \left(10\right)}\\&\simeq 0.45\,\mathrm{cm}.\end{aligned}
$$

ここでの結果は $E_1$ 以下に外挿しない。狭い帯の外では阻止能の別の近似が必要となる。

出口のエネルギーを下げる。そこで次の値に変更し、ほかの入力を固定する。
$$
E_{1}=3\,\mathrm{MeV}
$$

このときの値を ${\Delta x}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\Delta x}_{\mathrm{new}}&=\frac{- \left(3\right)^{2} + \left(5\right)^{2}}{2 \cdot \left(10\right)}\\&\simeq 0.8\,\mathrm{cm},\\{\Delta x}_{\mathrm{new}}-{\Delta x}&\simeq 0.35\,\mathrm{cm}.\end{aligned}
$$
:::

:::example[移行エネルギー分布の高い側]
微分率が $T^{-2}$ に比例する模型で、下限 $T_0$ と上限 $T_m$ の間を積分する。$T_0\le T_c\le T_m$ として、$T\ge T_c$ の衝突の割合を求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
T_{0}=1\,\mathrm{keV},\qquad T_{c}=4\,\mathrm{keV},\qquad T_{m}=16\,\mathrm{keV}
$$

この条件で求める量を ${p}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{p}&=\frac{- \frac{1}{T_{m}} + 1 \cdot \frac{1}{T_{c}}}{- \frac{1}{T_{m}} + 1 \cdot \frac{1}{T_{0}}}\\&=\frac{\frac{1}{\left(4\right)} - \frac{1}{\left(16\right)}}{\frac{1}{\left(1\right)} - \frac{1}{\left(16\right)}}\\&\simeq 0.2.\end{aligned}
$$

大きな移行は個数では少なくても、エネルギー損失には個数より大きい重みをもつ。

しきい値を上げて比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
T_{c}=8\,\mathrm{keV}
$$

このときの値を ${p}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{p}_{\mathrm{new}}&=\frac{\frac{1}{\left(8\right)} - \frac{1}{\left(16\right)}}{\frac{1}{\left(1\right)} - \frac{1}{\left(16\right)}}\\&\simeq 0.0666667,\\{p}_{\mathrm{new}}-{p}&\simeq -0.133333.\end{aligned}
$$
:::

## γ 線と中性子の物質との相互作用

細い単色ビームの未衝突強度 $I_\gamma(x)$ を区別する。一様物質で各過程の微視的断面積を $\sigma_i$、標的数密度を $n$ とすると厚さ $\dd x$ での衝突確率は $\mu\dd x$、$\mu=n\sum_i\sigma_i$ である。従って
$$
I_\gamma(x)=I_\gamma(0)e^{-\mu x}.
\label{eq-nuc-attenuation}
$$
混合物では成分ごとの $n_j\sigma_{ji}$ を足す。散乱光子も検出される広い配置では測定強度に散乱成分が加わる。

光電効果では束縛電子が放出され、空孔の緩和が後続する。Compton 散乱では光子がエネルギーの一部を電子に渡す。自由静止電子近似で、散乱角 $\theta$、入射と出射光子エネルギー $E,E'$ に対し、反跳電子の四元運動量の二乗を $m_e^2c^2$ に固定すると $m_ec^2(E-E')=EE'(1-\cos\theta)$ となる。電子陽電子対生成では標的核が運動量を受け取る。質量 $M$ の静止核が基底状態に残る場合、初期系の不変質量二乗はエネルギー単位で $M^2c^4+2Mc^2E$ である。しきい値では重心系で終粒子が相対運動しないので、これを $(M+2m_e)^2c^4$ に等置すると $E_{\rm th}=2m_ec^2(1+m_e/M)$ を得る。重い核では約 $1.022\,\mathrm{MeV}$ である。

中性子は電荷をもたず、主に弾性散乱、非弾性散乱、捕獲を起こす。捕獲の励起エネルギーはしばしば $\gamma$ 線に移る。弾性衝突では反跳核が物質に電離を残す。中性子でも未衝突ビームに対して :eqref[eq-nuc-attenuation] を使える。散乱後も系内に残る中性子については、位置、方向、エネルギーを変える輸送として扱う。

### 弾性散乱での中性子のエネルギー移行

非相対論的な中性子の質量を $m$、静止標的の質量を $Am$、入射速度を $v$ とする。重心の速度は $v/(A+1)$、重心系での中性子の速さは $Av/(A+1)$ である。重心系の散乱角を $\theta_*$ とすると、散乱後はこの二つの速度のなす角が $\theta_*$ になる。実験室系のエネルギー比は
$$
\frac{E'}{E}=\frac{A^2+1+2A\cos\theta_*}{(A+1)^2},
\qquad
\eta E\le E'\le E,\quad \eta=\left(\frac{A-1}{A+1}\right)^2.
$$
重心系で等方散乱なら $\cos\theta_*$ は区間 $[-1,1]$ で一様なので、$E'/E$ も $[\eta,1]$ で一様になる。平均のエネルギー移行割合は $(1-\eta)/2$ である。水素ではその値が $1/2$、質量比 $A=12$ では $24/169\simeq0.142$ となる。等方性の近似が破れるエネルギーでは、実際の微分断面積で平均を取る。

:::example[変換距離]
$\mu=0.2\,\mathrm{cm}^{-1}$ の物質で平均自由行程は $1/\mu=5\,\mathrm{cm}$ である。半価層は $\ln2/\mu\simeq3.47\,\mathrm{cm}$ となる。
:::

:::example[Compton 端]
$E=m_ec^2$ の光子が $\theta=\pi$ に散乱すると $E'=E/3$ である。電子には $2E/3\simeq0.341\,\mathrm{MeV}$ が渡る。出射光子が検出器外に逃げれば、この運動エネルギーが連続スペクトルの端を与える。
:::

### γ 線のエネルギー分布と除去断面積

自由な静止電子に散乱する光子では、電子の四元運動量の二乗が散乱前後で同じである。光子の運動量を $k,k'$ とすると $(p_e+k-k')^2=p_e^2$、$k^2=k'^2=0$ なので $m_ec^2(E-E')=EE'(1-\cos\theta)$ である。従って散乱角から $E'$ が決まり、光子が逃げた検出事象で電子に残るエネルギーは $E-E'$ である。

:::example[Compton 端]
散乱角 $\theta=\pi$ で光子のエネルギーが最小となる。入射エネルギー $E$ に対し電子の最大運動エネルギーを求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
E=\frac{331}{500}\,\mathrm{MeV},\qquad m_{ec2}=\frac{511}{1000}\,\mathrm{MeV}
$$

この条件で求める量を ${T_{\rm edge}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{T_{\rm edge}}&=E - \frac{E}{\frac{2 \cdot E}{m_{ec2}} + 1}\\&=\left(\frac{331}{500}\right) - \frac{\left(\frac{331}{500}\right)}{\frac{2 \cdot \left(\frac{331}{500}\right)}{\left(\frac{511}{1000}\right)} + 1}\\&\simeq 0.47765\,\mathrm{MeV}.\end{aligned}
$$

最大移行の位置は、全吸収ピークより低い。両者の差は逃げる光子のエネルギーである。

入射する光子のエネルギーを変更する。そこで次の値に変更し、ほかの入力を固定する。
$$
E=\frac{5}{4}\,\mathrm{MeV}
$$

このときの値を ${T_{\rm edge}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{T_{\rm edge}}_{\mathrm{new}}&=\left(\frac{5}{4}\right) - \frac{\left(\frac{5}{4}\right)}{1 + \frac{2 \cdot \left(\frac{5}{4}\right)}{\left(\frac{511}{1000}\right)}}\\&\simeq 1.03786\,\mathrm{MeV},\\{T_{\rm edge}}_{\mathrm{new}}-{T_{\rm edge}}&\simeq 0.560211\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[透過から減弱係数を推定]
散乱光子が受光器に戻らない細いビームを仮定する。厚さ $d$ と透過比 $u$ から減弱係数を求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
d=2\,\mathrm{cm},\qquad u=\frac{3}{5}
$$

この条件で求める量を ${\mu}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\mu}&=\frac{\left(-1\right) \cdot \log{\left(u \right)}}{d}\\&=- \frac{\log{\left(\left(\frac{3}{5}\right) \right)}}{\left(2\right)}\\&\simeq 0.255413\,\mathrm{cm}^{-1}.\end{aligned}
$$

大きな受光面では散乱成分を含めるので、同じ測定式で求める係数が経路の除去率と一致するとは限らない。

透過比が低い材料に変える。そこで次の値に変更し、ほかの入力を固定する。
$$
u=\frac{1}{2}
$$

このときの値を ${\mu}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\mu}_{\mathrm{new}}&=- \frac{\log{\left(\left(\frac{1}{2}\right) \right)}}{\left(2\right)}\\&\simeq 0.346574\,\mathrm{cm}^{-1},\\{\mu}_{\mathrm{new}}-{\mu}&\simeq 0.0911608\,\mathrm{cm}^{-1}.\end{aligned}
$$
:::

:::example[エネルギー別の透過]
単一エネルギーの光子について、平均自由行程 $\ell$ と厚さ $d$ を用いる。元のビームのうち一度も相互作用せずに透過する光子を数える。

量の定義に従い、共通の単位で次の値を用いる。
$$
d=2\,\mathrm{cm},\qquad l=1\,\mathrm{cm}
$$

この条件で求める量を ${u}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{u}&=e^{\frac{\left(-1\right) \cdot d}{l}}\\&=e^{- \frac{\left(2\right)}{\left(1\right)}}\\&\simeq 0.135335.\end{aligned}
$$

多エネルギーの入射では透過は各エネルギー成分の和になる。単一の指数に当てはめると見かけの係数が厚みに依存しうる。

より長い行程をもつ成分と比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
l=4\,\mathrm{cm}
$$

このときの値を ${u}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{u}_{\mathrm{new}}&=e^{- \frac{\left(2\right)}{\left(4\right)}}\\&\simeq 0.606531,\\{u}_{\mathrm{new}}-{u}&\simeq 0.471195.\end{aligned}
$$
:::

:::example[競合する最初の衝突]
一様媒質で吸収の巨視的断面積を $\Sigma_a$、弾性散乱を $\Sigma_s$ とする。厚さ $d$ までに最初の事象が吸収となる確率を積分する。

量の定義に従い、共通の単位で次の値を用いる。
$$
S_{a}=\frac{2}{5}\,\mathrm{cm}^{-1},\qquad S_{s}=\frac{3}{5}\,\mathrm{cm}^{-1},\qquad d=1\,\mathrm{cm}
$$

この条件で求める量を ${p_a}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{p_a}&=\frac{S_{a} \cdot \left(1 - e^{d \cdot \left(- S_{a} - S_{s}\right)}\right)}{S_{a} + S_{s}}\\&=\frac{\left(\frac{2}{5}\right) \cdot \left(1 - e^{\left(1\right) \cdot \left(- \left(\frac{2}{5}\right) - \left(\frac{3}{5}\right)\right)}\right)}{\left(\frac{2}{5}\right) + \left(\frac{3}{5}\right)}\\&\simeq 0.252848.\end{aligned}
$$

散乱後に改めて吸収される事象は別に扱う。厚い極限でも最初の事象の種別で選ぶ割合には上限がある。

物質を厚くして極限への近づきを見る。そこで次の値に変更し、ほかの入力を固定する。
$$
d=3\,\mathrm{cm}
$$

このときの値を ${p_a}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{p_a}_{\mathrm{new}}&=\frac{\left(\frac{2}{5}\right) \cdot \left(1 - e^{\left(3\right) \cdot \left(- \left(\frac{2}{5}\right) - \left(\frac{3}{5}\right)\right)}\right)}{\left(\frac{2}{5}\right) + \left(\frac{3}{5}\right)}\\&\simeq 0.380085,\\{p_a}_{\mathrm{new}}-{p_a}&\simeq 0.127237.\end{aligned}
$$
:::

## 放射線検出器

気体検出器では一次電離の電子とイオンを電場で回収する。電離箱では一次電荷を主に測り、比例計数管では電子雪崩による増倍を使う。Geiger–Müller 管では放電が広がるため信号は粒子が渡したエネルギーへの比例性を失う。シンチレーション検出器は励起からの光を測る。半導体検出器では電子正孔対を作り、少ない生成エネルギーで多数のキャリアを得る。中性子検出では事象を荷電粒子や光子に変換する。

全エネルギー $E$ が検出器に沈着し、キャリア一対当たりの平均生成エネルギーを $w$ とすると、生成した対の個数 $n$ の平均は $E/w$ である。分散を $\operatorname{Var}(n)=F_{\rm Fano}E/w$ と書くときの $F_{\rm Fano}$ を :term[Fano 因子]{en="Fano factor" yomi="ふぁのいんし"} という。エネルギーの推定値は $wn$ なので、その標準偏差は $\sqrt{F_{\rm Fano}wE}$ である。Gaussian 近似で半値全幅は $2\sqrt{2\ln2}$ 倍である。電子回路の雑音と電荷回収のばらつきも測定幅に寄与する。

真の事象率を $r$、記録率を $m$ とし、事象を記録した直後の不感時間 $\tau_d$ 内では新たな記録も不感時間の延長もないとする。長時間で記録可能な時間割合は $1-m\tau_d$ である。到着が定常 Poisson 過程なら $m=r(1-m\tau_d)$ となり、$r=m/(1-m\tau_d)$ を得る。不感時間中の事象も不感時間を伸ばす模型では、到着前の $\tau_d$ が空である確率 $e^{-r\tau_d}$ を使い $m=re^{-r\tau_d}$ となる。

:::example[放射能の逆算]
バックグラウンドを引いた光電ピークが毎秒 20 計数で、光子の放出割合が $0.5$、ピーク効率が $0.1$ とする。不感時間を補正した計数なら、放射能は $20/(0.5\cdot0.1)=400\,\mathrm{Bq}$ である。
:::

:::example[不感時間]
非延長型で $m=800\,\mathrm{s}^{-1}$、$\tau_d=10^{-4}\,\mathrm s$ なら $m\tau_d=0.08$、$r\simeq870\,\mathrm{s}^{-1}$ である。
:::

### 計数の統計と効率

十分短い測定で放射能や背景が一定なら、独立な事象のカウントを Poisson 分布で近似できる。線源を含む時間 $T_s$ の窓で $C_s$、別の背景測定時間 $T_b$ で $C_b$ を得たとする。背景の平均率が同じなら、線源由来の率の推定値は $\widehat r=C_s/T_s-C_b/T_b$ である。各 Poisson カウントでは分散が平均に等しく、測定が独立なら分散の推定値は $C_s/T_s^2+C_b/T_b^2$ となる。効率 $\epsilon$ と分岐 $b$ が既知なら、率と標準偏差をそれぞれ $b\epsilon$ で割る。計数が少ないときは対称な Gaussian 誤差より、カウントの尤度を直接使う方が幅を記述しやすい。

:::example[背景を差し引いた率]
$T_s=T_b=100\,\mathrm s$、$C_s=500,C_b=100$ なら、差し引きの率は $4\,\mathrm{s}^{-1}$、統計誤差は $\sqrt{600}/100\simeq0.245\,\mathrm{s}^{-1}$ である。背景を引くと平均は減るが、測定誤差の二乗は加わる。
:::

### エネルギー応答を折り込む

入射エネルギーが $E$ のとき、検出された信号エネルギーの確率密度を $R(E_{\rm obs}|E)$ とする。未検出を除いた規格化で $\int R\dd E_{\rm obs}=1$ とし、検出確率を $\epsilon(E)$ と分けると、観測スペクトルは
$$
s_{\rm obs}(E_{\rm obs})=b(E_{\rm obs})+\int R(E_{\rm obs}|E)\epsilon(E)s(E)\dd E
$$
となる。$s(E)$ は線源の入射スペクトル、$b(E_{\rm obs})$ は背景である。モノエネルギー線 $N\delta(E-E_0)$ は $N\epsilon(E_0)R(E_{\rm obs}|E_0)$ を与える。核の自然幅が検出器の幅より十分に狭ければ、観測したピーク幅は主に $R$ による。単一エネルギーを全て吸収しない事象は連続成分や逃避ピークを作り、単一 Gaussian で表せない。

:::example[二光子の逃避]
対生成で電子・陽電子の運動エネルギーが全て測られ、陽電子消滅で生じた $0.511\,\mathrm{MeV}$ の光子が一つ逃げたとする。観測値は入射光子のエネルギーより $0.511\,\mathrm{MeV}$ 小さい。両方逃げれば差は $1.022\,\mathrm{MeV}$ になる。全吸収とこの二つの逃避が区別できる装置では、既知の差がエネルギー目盛りの検査に使える。
:::

### 不感時間とパルスの分布

非延長型の不感時間 $\tau$ では、記録した一事象の後の $\tau$ の間だけ次の事象を取りこぼす。長時間の測定で記録率が $r$ なら時間の割合 $r\tau$ は受け付け不能である。一様 Poisson 入射の真の率 $R$ は受け付け時間に対しても一定なので $r=R(1-r\tau)$ となる。これは各事象で不感時間を再開する延長型とはモデルが異なる。後者では各入射の直前の区間 $\tau$ に事象がない確率 $e^{-R\tau}$ が検出率に掛かる。

:::example[非延長型の補正]
記録率 $r$ と非延長型不感時間 $\tau$ を使って真の率を復元する。

量の定義に従い、共通の単位で次の値を用いる。
$$
r=1000\,\mathrm{s}^{-1},\qquad \tau=\frac{1}{10000}\,\mathrm{s}
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=\frac{r}{1 - r \cdot \tau}\\&=\frac{\left(1000\right)}{- \left(1000\right) \cdot \left(\frac{1}{10000}\right) + 1}\\&\simeq 1111.11\,\mathrm{s}^{-1}.\end{aligned}
$$

記録率が $1/\tau$ に近づくほど逆補正は不安定になる。

入射強度が増している測定と比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
r=4000\,\mathrm{s}^{-1}
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=\frac{\left(4000\right)}{- \left(4000\right) \cdot \left(\frac{1}{10000}\right) + 1}\\&\simeq 6666.67\,\mathrm{s}^{-1},\\{R}_{\mathrm{new}}-{R}&\simeq 5555.56\,\mathrm{s}^{-1}.\end{aligned}
$$
:::

:::example[独立した背景測定]
試料ありと背景のみのカウント数を $C,B$ とする。それぞれの時間を $t,t_b$ とした純率の標準偏差を Poisson 分散から求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
C=900,\qquad B=100,\qquad t=30\,\mathrm{s},\qquad t_{b}=30\,\mathrm{s}
$$

この条件で求める量を ${\sigma_r}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\sigma_r}&=\sqrt{\frac{B}{t_{b}^{2}} + \frac{C}{t^{2}}}\\&=\sqrt{\frac{\left(100\right)}{\left(30\right)^{2}} + \frac{\left(900\right)}{\left(30\right)^{2}}}\\&\simeq 1.05409\,\mathrm{s}^{-1}.\end{aligned}
$$

後の設定では同じ背景率で背景の時間を延ばしている。純率を引き算しても分散は足される。

背景の測定を長くする利得を比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
B=200,\qquad t_{b}=60\,\mathrm{s}
$$

このときの値を ${\sigma_r}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\sigma_r}_{\mathrm{new}}&=\sqrt{\frac{\left(200\right)}{\left(60\right)^{2}} + \frac{\left(900\right)}{\left(30\right)^{2}}}\\&\simeq 1.0274\,\mathrm{s}^{-1},\\{\sigma_r}_{\mathrm{new}}-{\sigma_r}&\simeq -0.0266902\,\mathrm{s}^{-1}.\end{aligned}
$$
:::

:::example[Fano 因子と分解能]
一対当たりの平均生成エネルギーを $w$、生成した対の個数の分散を $F E/w$ とする。Gaussian ピークの半値全幅は標準偏差の $\sqrt{8\ln2}$ 倍である。

量の定義に従い、共通の単位で次の値を用いる。
$$
F=\frac{1}{10},\qquad w=3\,\mathrm{eV},\qquad E=100000\,\mathrm{eV}
$$

この条件で求める量を ${\Delta E}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\Delta E}&=\sqrt{8 \cdot E \cdot F \cdot w \cdot \log{\left(2 \right)}}\\&=2 \cdot \sqrt{2} \cdot \sqrt{\left(100000\right) \cdot \left(3\right) \cdot \left(\frac{1}{10}\right)} \cdot \sqrt{\log{\left(2 \right)}}\\&\simeq 407.867\,\mathrm{eV}.\end{aligned}
$$

絶対幅は増えるが、ピークのエネルギーで割った相対幅は改善する。

同じ検出器で高いエネルギーを測る。そこで次の値に変更し、ほかの入力を固定する。
$$
E=400000\,\mathrm{eV}
$$

このときの値を ${\Delta E}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\Delta E}_{\mathrm{new}}&=2 \cdot \sqrt{2} \cdot \sqrt{\left(3\right) \cdot \left(400000\right) \cdot \left(\frac{1}{10}\right)} \cdot \sqrt{\log{\left(2 \right)}}\\&\simeq 815.734\,\mathrm{eV},\\{\Delta E}_{\mathrm{new}}-{\Delta E}&\simeq 407.867\,\mathrm{eV}.\end{aligned}
$$
:::

:::example[幾何効率の厳密な球冠積分]
小さな等方点源が半径 $a$ の円板を正面に見込む。軸上距離 $d$ では半頂角の余弦が $d/\sqrt{a^2+d^2}$ である。

量の定義に従い、共通の単位で次の値を用いる。
$$
a=2\,\mathrm{cm},\qquad d=4\,\mathrm{cm}
$$

この条件で求める量を ${\epsilon_\Omega}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\epsilon_\Omega}&=\frac{1 - \frac{d}{\sqrt{a^{2} + d^{2}}}}{2}\\&=- \frac{\left(4\right)}{2 \cdot \sqrt{\left(2\right)^{2} + \left(4\right)^{2}}} + \frac{1}{2}\\&\simeq 0.0527864.\end{aligned}
$$

立体角は $2\pi(1-\cos\theta)$ なので全立体角で割った結果である。検出器の吸収効率はこの値に別途掛ける。

検出器を遠ざけた設定を比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
d=8\,\mathrm{cm}
$$

このときの値を ${\epsilon_\Omega}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\epsilon_\Omega}_{\mathrm{new}}&=- \frac{\left(8\right)}{2 \cdot \sqrt{\left(2\right)^{2} + \left(8\right)^{2}}} + \frac{1}{2}\\&\simeq 0.0149287,\\{\epsilon_\Omega}_{\mathrm{new}}-{\epsilon_\Omega}&\simeq -0.0378577.\end{aligned}
$$
:::

## 加速器の原理

電場は粒子に仕事をし、静磁場は軌道を曲げる。静電場で電位 $\phi_i$ から $\phi_f$ まで進んだ粒子の運動エネルギーの増分は $q(\phi_i-\phi_f)$ である。加速方向を選んだとき、電位差の大きさ $V$ に対して利得は $|q|V$ となる。直線加速器では高周波電場が加速向きになる位相で粒子を通過させる。同期粒子の速さを $v$、周波数を $f$ とし、隣接する空隙の位相が反転する構成で粒子が位相反転一回分で到達するなら、空隙の間隔は $v/(2f)$ となる。

一様な磁場 $B$ に垂直な運動では Lorentz 力から $pv/R=|q|vB$、従って $p=|q|BR$ である。非相対論的サイクロトロンでは $f=|q|B/(2\pi m)$ が運動量によらない。相対論的には $\gamma$ が分母に入り、一定の周波数では位相がずれる。シンクロトロンは運動量に合わせて磁場と高周波系を同期させ、ほぼ一定の周回軌道を保つ。

:::example[静電加速]
一価イオンの静電加速で $V=300\,\mathrm{kV}$ なら $300\,\mathrm{keV}$ を得る。二価イオンは同じ電位差で $600\,\mathrm{keV}$ を得る。
:::

:::example[磁気剛性]
電荷の大きさが同じで運動量が 2 倍なら、同じ磁場で曲率半径は 2 倍である。半径を固定して加速するには、この式から磁場も 2 倍になる。
:::

### 粒子の束と位相安定性

高周波電圧の一回当たりのエネルギー増分を $qV_0\sin\phi$ とし、予定された軌道を通る同期粒子の位相を $\phi_s$ とする。粒子が少し早くまたは遅く着いたときの増分の差は $qV_0\cos\phi_s\,\Delta\phi$ である。相対的な運動量差が次の到着位相を変えるような系では、増分の位相勾配と周回時間のエネルギー勾配の符号で位相の復元が決まる。

連続時間の簡単な模型を $\dot{\Delta\phi}=a\Delta E,\ \dot{\Delta E}=b\Delta\phi$ と書く。$a$ は到着位相のエネルギー応答、$b$ は電圧による加速応答である。定数 $a,b$ で $\ddot{\Delta\phi}=ab\Delta\phi$ となるため、$ab<0$ なら振動、$ab>0$ なら指数的なずれである。$\sin\phi$ の勾配が最大の位相と加速電圧が最大の位相は異なるので、加速率と位相安定性を区別する。この模型は同期粒子の軌道付近の線形化である。

:::example[共鳴周波数のずれ]
相対論的な質量の効果で $\gamma$ が 1 から $1.01$ に変わると、一定磁場での周波数は $1/1.01\simeq0.9901$ 倍になる。一周では約 1 パーセントでも、位相を補正せずに何周も加速すれば空隙を通る位相がずれていく。
:::

### 磁気剛性と周回周波数

磁場に垂直な運動では $\dot{\vb p}=q\vb v\times \vb B$ である。速さが一定なら運動量の方向だけが毎秒 $v/R$ 回転し、力の大きさは $pv/R$ となる。よって $p=|q|BR$ であり、これは相対論的な場合にも成立する。周回時間は $2\pi R/v=2\pi\gamma m/(|q|B)$ となり、固定磁場での加速に伴う周波数のずれを与える。

:::example[電荷の違いと剛性]
同じ $p$ の二種を曲げる。GeV/$c$、T、m の単位で $p=0.299792458\,zBR$ とする。

量の定義に従い、共通の単位で次の値を用いる。
$$
p=\frac{3}{10}\,\mathrm{GeV}/c,\qquad c_{p}=\frac{149896229}{500000000},\qquad z=1,\qquad B=1\,\mathrm{T}
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=\frac{p}{B \cdot c_{p} \cdot z}\\&=\frac{\left(\frac{3}{10}\right)}{\left(1\right)^{2} \cdot \left(\frac{149896229}{500000000}\right)}\\&\simeq 1.00069\,\mathrm{m}.\end{aligned}
$$

電荷の絶対値を二倍にすると、同じ磁場で曲率は二倍である。

二価のイオンと比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
z=2
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=\frac{\left(\frac{3}{10}\right)}{\left(1\right) \cdot \left(2\right) \cdot \left(\frac{149896229}{500000000}\right)}\\&\simeq 0.500346\,\mathrm{m},\\{R}_{\mathrm{new}}-{R}&\simeq -0.500346\,\mathrm{m}.\end{aligned}
$$
:::

:::example[周回周波数のずれ]
固定磁場の非相対論的予測を $f_0$、Lorentz 因子を $\gamma$ とする。周回周波数を計算する。

量の定義に従い、共通の単位で次の値を用いる。
$$
f_{0}=12\,\mathrm{MHz},\qquad \gamma=\frac{101}{100}
$$

この条件で求める量を ${f}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{f}&=\frac{f_{0}}{\gamma}\\&=\frac{\left(12\right)}{\left(\frac{101}{100}\right)}\\&\simeq 11.8812\,\mathrm{MHz}.\end{aligned}
$$

運動量が大きくなると速度だけでなく曲率半径も変化する。固定周波数とのずれが周回ごとに累積する。

加速後の値を比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
\gamma=\frac{11}{10}
$$

このときの値を ${f}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{f}_{\mathrm{new}}&=\frac{\left(12\right)}{\left(\frac{11}{10}\right)}\\&\simeq 10.9091\,\mathrm{MHz},\\{f}_{\mathrm{new}}-{f}&\simeq -0.972097\,\mathrm{MHz}.\end{aligned}
$$
:::

:::example[一空隙の位相利得]
峰電圧 $V$ の高周波を、エネルギー利得が最大となる位相から $\phi$ だけずれて通過する。電荷数 $z$ の粒子を考える。

量の定義に従い、共通の単位で次の値を用いる。
$$
z=1,\qquad V=80\,\mathrm{kV},\qquad \phi=\frac{\pi}{6}\,\mathrm{rad}
$$

この条件で求める量を ${\Delta E}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\Delta E}&=V \cdot z \cdot \cos{\left(\phi \right)}\\&=\left(1\right) \cdot \left(80\right) \cdot \cos{\left(\left(\frac{\pi}{6}\right) \right)}\\&\simeq 69.282\,\mathrm{keV}.\end{aligned}
$$

通過時間が周期に比べて短いとして電圧を評価した。長い空隙では時間平均も必要である。

位相が遅れた粒子と比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
\phi=\frac{\pi}{3}\,\mathrm{rad}
$$

このときの値を ${\Delta E}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\Delta E}_{\mathrm{new}}&=\left(1\right) \cdot \left(80\right) \cdot \cos{\left(\left(\frac{\pi}{3}\right) \right)}\\&=40\,\mathrm{keV},\\{\Delta E}_{\mathrm{new}}-{\Delta E}&\simeq -29.282\,\mathrm{keV}.\end{aligned}
$$
:::

:::example[粒子流とビーム電流]
粒子率を $R$、電荷数を $z$ とする。電荷の大きさは $e_0$ で与える。

量の定義に従い、共通の単位で次の値を用いる。
$$
R=1000000000\,\mathrm{s}^{-1},\qquad z=1,\qquad e_{0}=\frac{801088317}{5000000000000000000000000000}\,\mathrm{C}
$$

この条件で求める量を ${I}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{I}&=R \cdot e_{0} \cdot z\\&=\left(1000000000\right) \cdot \left(1\right) \cdot \left(\frac{801088317}{5000000000000000000000000000}\right)\\&\simeq 1.60218\times10^{-10}\,\mathrm{A}.\end{aligned}
$$

パルスビームの時間平均を求めるときは $R$ も同じ時間間隔で平均する。

照射率を増やした値を比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
R=4000000000\,\mathrm{s}^{-1}
$$

このときの値を ${I}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{I}_{\mathrm{new}}&=\left(1\right) \cdot \left(4000000000\right) \cdot \left(\frac{801088317}{5000000000000000000000000000}\right)\\&\simeq 6.40871\times10^{-10}\,\mathrm{A},\\{I}_{\mathrm{new}}-{I}&\simeq 4.80653\times10^{-10}\,\mathrm{A}.\end{aligned}
$$
:::

## 放射年代測定

放射性の親核数 $P$ を年代の指標にするには、初期量または娘核の増分を推定する。閉じた系で娘核が安定、分岐が 1 なら、$\dot D=\lambda P$ と $\dot P=-\lambda P$ を加えて $D-D_0=P_0-P=P(e^{\lambda t}-1)$ を得る。保存された安定同位体数 $R$ で割ると
$$
\frac{D}{R}=\left(\frac{D}{R}\right)_0+\frac{P}{R}(e^{\lambda t}-1).
$$
生成時刻と初期比が同じ鉱物群でこの直線を当てはめると、傾きから年代が求まる。娘核の移動や混合で閉鎖系の仮定が崩れると、傾きの解釈を検討する必要がある。

炭素年代測定では生物の交換が止まった後に ${}^{14}\mathrm C$ が減る。初期の同位体比を $R_0$、測定値を $R$ とすると単純模型で $t=\lambda^{-1}\ln(R_0/R)$ である。実年代への換算には大気の初期比の変動と試料由来のずれを扱う校正が必要である。測る量が親核の残存比か、娘核の増分かで入力が異なる。

:::example[残存比]
$T_{1/2}=5730$ 年とし、初期比の $1/4$ が残っていれば単純模型の年代は $11460$ 年となる。
:::

:::example[等時線]
娘核比の回帰直線の傾きが 3 なら $e^{\lambda t}=4$ である。従って $t=2T_{1/2}$ となる。切片が非ゼロでもこの傾きの読み方は変わらない。
:::

### 年代の誤差と混入の影響

残存比 $f=R/R_0$ が分かる場合、$t=-\ln f/\lambda$ なので、微小誤差は $\delta t\simeq-\delta f/(\lambda f)$ である。一定の絶対誤差でも、残存比が小さい測定ほど年代の誤差が大きくなる。比を測るときはゼロとの差を統計的に区別できるかも必要である。

元の分率を $f$ とし、混入した現代炭素が参照同位体 ${}^{12}\mathrm C$ の総数に占める割合を $x$ とする。両成分の ${}^{14}\mathrm C/{}^{12}\mathrm C$ を共通の初期比で割り、現代炭素の値を 1 とすれば、混合後の分率は $f_{\rm mix}=(1-x)f+x$ である。$f<1$ なら $x>0$ により残存比が増え、混入を省いた年代は小さくなる。これが試料の前処理が重要になる理由である。

:::example[少量の混入]
真の比が $f=0.02$、現代炭素が $x=0.01$ 混ざると、観測比は $0.0298$ である。年代のずれは $-\lambda^{-1}\ln(0.0298/0.02)$、半減期を単位にすると約 $-0.58$ である。古い試料では同じ割合の混入でも、元の残存量に対して相対的な増加が大きい。
:::

### 閉鎖系と等時線

娘核が安定で親核の漏出もない場合、現在の親数 $P$ と放射起源の娘数 $D^*$ の和は初期の親数である。$P_0/P=e^{\lambda t}$ から $D^*/P=e^{\lambda t}-1$ を得る。初期の娘数が未知の場合は安定な同元素同位体 $S$ を参照として、共通の初期比 $D_0/S$ をもつ複数の試料について $D/S=D_0/S+(e^{\lambda t}-1)P/S$ を用いる。傾きから年代、切片から初期比を求める。各試料が共通の初期比と年代をもつかは、直線からの残差も使って検討する。

:::example[親と放射起源の娘]
生成後に閉鎖した試料で、放射起源の娘核数 $D^*$ と現在の親核数 $P$ の比を $r=D^*/P$ とする。親核は分岐比 1 で安定娘核を作るものとし、半減期 $T$ から年代を求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
T=100\,\mathrm{Ma},\qquad r=1
$$

この条件で求める量を ${t}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{t}&=\frac{T \cdot \log{\left(r + 1 \right)}}{\log{\left(2 \right)}}\\&=\frac{\left(100\right) \cdot \log{\left(\left(1\right) + 1 \right)}}{\log{\left(2 \right)}}\\&=100\,\mathrm{Ma}.\end{aligned}
$$

初めに含まれた娘がある場合は、その個数を引いた放射起源の分だけを $r$ に使う。

娘の多い試料と比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
r=3
$$

このときの値を ${t}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{t}_{\mathrm{new}}&=\frac{\left(100\right) \cdot \log{\left(\left(3\right) + 1 \right)}}{\log{\left(2 \right)}}\\&=200\,\mathrm{Ma},\\{t}_{\mathrm{new}}-{t}&=100\,\mathrm{Ma}.\end{aligned}
$$
:::

:::example[等時線の勾配]
半減期 $T$ の系で、試料間の比 $P/S,D/S$ を作る。勾配 $m$ から閉鎖後の時間を求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
T=4\,\mathrm{Ga},\qquad m=\frac{1}{10}
$$

この条件で求める量を ${t}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{t}&=\frac{T \cdot \log{\left(m + 1 \right)}}{\log{\left(2 \right)}}\\&=\frac{\left(4\right) \cdot \log{\left(\left(\frac{1}{10}\right) + 1 \right)}}{\log{\left(2 \right)}}\\&\simeq 0.550014\,\mathrm{Ga}.\end{aligned}
$$

傾きが小さい場合は $t\simeq Tm/\ln2$ だが、勾配の誤差が相対的に大きくなりやすい。

勾配が大きい等時線と比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
m=\frac{1}{5}
$$

このときの値を ${t}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{t}_{\mathrm{new}}&=\frac{\left(4\right) \cdot \log{\left(\left(\frac{1}{5}\right) + 1 \right)}}{\log{\left(2 \right)}}\\&\simeq 1.05214\,\mathrm{Ga},\\{t}_{\mathrm{new}}-{t}&\simeq 0.502124\,\mathrm{Ga}.\end{aligned}
$$
:::

:::example[炭素年代の一次誤差]
残存比 $f$ の標準偏差を $\sigma_f$ とすると、$t=-T\ln f/\ln2$ の微分の絶対値を掛ける。

量の定義に従い、共通の単位で次の値を用いる。
$$
T=5730\,\mathrm{yr},\qquad s=\frac{1}{200},\qquad f=\frac{1}{2}
$$

この条件で求める量を ${\sigma_t}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\sigma_t}&=\frac{T \cdot s}{f \cdot \log{\left(2 \right)}}\\&=\frac{\left(5730\right) \cdot \left(\frac{1}{200}\right)}{\left(\frac{1}{2}\right) \cdot \log{\left(2 \right)}}\\&\simeq 82.6664\,\mathrm{yr}.\end{aligned}
$$

同じ絶対誤差では古い試料の方が年代の誤差が大きい。

残存比の小さい試料に変える。そこで次の値に変更し、ほかの入力を固定する。
$$
f=\frac{1}{10}
$$

このときの値を ${\sigma_t}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\sigma_t}_{\mathrm{new}}&=\frac{\left(5730\right) \cdot \left(\frac{1}{200}\right)}{\left(\frac{1}{10}\right) \cdot \log{\left(2 \right)}}\\&\simeq 413.332\,\mathrm{yr},\\{\sigma_t}_{\mathrm{new}}-{\sigma_t}&\simeq 330.666\,\mathrm{yr}.\end{aligned}
$$
:::

:::example[二成分の混合年代]
残存比の異なる成分 $f_1,f_2$ を混ぜ、第二成分の炭素割合を $w$ とする。混合比を先に平均し、それを閉鎖系の式へ入れる。

量の定義に従い、共通の単位で次の値を用いる。
$$
T=5730\,\mathrm{yr},\qquad w=0,\qquad f_{1}=\frac{1}{8},\qquad f_{2}=1
$$

この条件で求める量を ${t_{\rm app}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{t_{\rm app}}&=\frac{- T \cdot \log{\left(f_{1} \cdot \left(1 - w\right) + f_{2} \cdot w \right)}}{\log{\left(2 \right)}}\\&=- \frac{\left(5730\right) \cdot \log{\left(\left(0\right) \cdot \left(1\right) + \left(\frac{1}{8}\right) \cdot \left(1 - \left(0\right)\right) \right)}}{\log{\left(2 \right)}}\\&=17190\,\mathrm{yr}.\end{aligned}
$$

対数を取る前に残存比を混ぜるので、二つの年代の重み付き算術平均にはならない。

若い炭素がわずかに入る場合を考える。そこで次の値に変更し、ほかの入力を固定する。
$$
w=\frac{1}{20}
$$

このときの値を ${t_{\rm app}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{t_{\rm app}}_{\mathrm{new}}&=- \frac{\left(5730\right) \cdot \log{\left(\left(1\right) \cdot \left(\frac{1}{20}\right) + \left(\frac{1}{8}\right) \cdot \left(1 - \left(\frac{1}{20}\right)\right) \right)}}{\log{\left(2 \right)}}\\&\simeq 14709.1\,\mathrm{yr},\\{t_{\rm app}}_{\mathrm{new}}-{t_{\rm app}}&\simeq -2480.86\,\mathrm{yr}.\end{aligned}
$$
:::

## 演習問題

::::exercise{level=basic}
半減期 4 時間の線源の放射能が $160\,\mathrm{Bq}$ である。12 時間後の放射能を求めよ。
:::solution
:eqref[eq-nuc-decay-law] で半減期が 3 回経過するので $160\cdot2^{-3}=20\,\mathrm{Bq}$ である。
:::
::::

::::exercise{level=basic}
$\lambda_2=10\lambda_1$、$b=1$ の系列の長時間における放射能比を求めよ。
:::solution
:eqref[eq-nuc-two-member-decay] で $e^{-\lambda_2t}/e^{-\lambda_1 t}\to0$ だから $\mathcal A_2/\mathcal A_1=10/9$ に近づく。
:::
::::

::::exercise{level=basic}
$\mu=0.5\,\mathrm{cm}^{-1}$ の物質を $4\,\mathrm{cm}$ 通過した未衝突光子の割合を求めよ。
:::solution
厚さと減弱係数の積は 2 なので、残る割合は $e^{-2}\simeq0.1353$ である。
:::
::::

::::exercise{level=basic}
非延長型不感時間が $2\times10^{-4}\,\mathrm s$、計数率が $1000\,\mathrm{s}^{-1}$ のとき真の率を求めよ。
:::solution
記録可能時間の割合は $1-m\tau_d=0.8$ である。従って真の率は $1000/0.8=1250\,\mathrm{s}^{-1}$ である。
:::
::::

::::exercise{level=advanced}
$\lambda_2>\lambda_1>0$ かつ $b>0,N_{20}=0$ の系列で、娘核数最大の時刻を求めよ。
:::solution
:eqref[eq-nuc-two-member-decay] の微分がゼロとなる条件は $\lambda_2 e^{-\lambda_2t}=\lambda_1 e^{-\lambda_1t}$ である。解は $\ln(\lambda_2/\lambda_1)/(\lambda_2-\lambda_1)$ である。微分は最初正で、比 $(\lambda_2/\lambda_1)e^{-(\lambda_2-\lambda_1)t}$ が単調に 0 へ減るので、その後は負となる。
:::
::::

::::exercise{level=advanced}
Compton 散乱で入射エネルギー $E=2m_ec^2$ のとき反跳電子の最大エネルギーを求めよ。
:::solution
$E'=E/[1+(E/m_ec^2)(1-\cos\theta)]$ の最小値は $\theta=\pi$ で $E/5$ となる。従って反跳電子の最大運動エネルギーは $4E/5=8m_ec^2/5$ である。
:::
::::

::::exercise{level=advanced}
安定娘核を生む分岐が $b$ のとき等時線の傾きを導け。
:::solution
$\dot D=b\lambda P=-b\dot P$ なので $D-D_0=b(P_0-P)=bP(e^{\lambda t}-1)$ である。$R$ で割った回帰の傾きを $s$ とすれば $s=b(e^{\lambda t}-1)$、$b>0$ で $t=\lambda^{-1}\ln(1+s/b)$ となる。
:::
::::

::::exercise{level=advanced}
時間 $T$ の間の期待崩壊数を積分で求め、初期放射能に $T$ を掛けた近似の相対的なずれの先頭項を求めよ。
:::solution
$\int_0^T \lambda N_0e^{-\lambda t}\dd t=N_0(1-e^{-\lambda T})$ である。一定率近似 $N_0\lambda T$ に対して真の値の比は $1-\lambda T/2+O((\lambda T)^2)$ となる。従って近似は先頭で $\lambda T/2$ だけ過大である。
:::
::::

## 確認テスト

::::quiz{answer="8" tol="1e-6"}
放射能が $64\,\mathrm{Bq}$ の線源について、半減期が 3 回過ぎたときの放射能を Bq で求めよ。
:::explanation
半減期ごとに $1/2$ 倍なので $64/8=8\,\mathrm{Bq}$ である。
:::
::::

::::quiz
永続平衡で全分岐が娘核を生むとき、娘核の放射能は何に近いか。

- [x] 親核の放射能
- [ ] 親核の核数
- [ ] 時刻によらない定数

:::explanation
$\lambda_2\gg\lambda_1$ では娘核の消滅が親核からの生成にほぼ追いつくため、$\lambda_2N_2\simeq\lambda_1N_1$ となる。両者は時間とともに減る。
:::
::::

::::quiz{answer="4" tol="1e-6"}
$S=a/E$ の連続減速模型で、エネルギーを 2 倍にすると飛程は何倍になるか。
:::explanation
積分により $R=E_0^2/(2a)$ なので $2^2=4$ 倍である。
:::
::::

::::quiz
磁場だけが粒子に力を及ぼす場合、曲げた粒子の運動エネルギーはどうなるか。

- [ ] 曲率が大きいほど増える。
- [x] 一定である。
- [ ] 符号が反対の粒子では減る。

:::explanation
Lorentz 力 $q\vb v\times\vb B$ は速度に直交しており、仕事率 $(q\vb v\times\vb B)\cdot\vb v=0$ となる。電荷の符号は曲がる向きに関係する。
:::
::::

::::quiz{answer="2" tol="1e-6"}
全分岐が安定娘核を生む閉鎖系で等時線の傾きが 3 である。年代は半減期の何倍か。
:::explanation
$e^{\lambda t}-1=3$ なので $\lambda t=\ln4=2\ln2$ となり、半減期の 2 倍である。
:::
::::

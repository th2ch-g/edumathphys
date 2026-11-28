---
title: 原子核天体物理と元素合成
description: 熱平均反応率、恒星の核燃焼、中性子捕獲、初期宇宙と高密度天体を扱う。
notation:
  - symbol: 'S(E)'
    meaning: Coulomb 透過と運動学因子を除いた天体物理的 S 因子
  - symbol: 'E_G,E_0'
    meaning: Gamow エネルギーと Gamow ピークの近似位置
  - symbol: 'Y_i'
    meaning: 核種 i のバリオン当たり個数
  - symbol: 'H'
    meaning: 宇宙の膨張率
---

## この章の目標

- Maxwell 分布から熱核反応率と Gamow ピークを導く。
- 水素、ヘリウムと重い核の燃焼を反応網で追う。
- 中性子捕獲と $\beta$ 崩壊の競合を使って元素合成を分類する。
- 宇宙の膨張と星の崩壊が反応の平衡を破る条件を考える。

## 熱核反応率と天体物理的 S 因子

### 相対速度の平均

異なる質量 $m_1,m_2$ の古典的な非相対論粒子が共通の温度 $T$ にあるとする。運動エネルギーの和を重心と相対運動に分けると分布の指数も分離し、相対速度 $\vb v$ の密度は $(\mu/(2\pi k_{\rm B}T))^{3/2}\exp[-\mu v^2/(2k_{\rm B}T)]$ である。換算質量は $\mu=m_1m_2/(m_1+m_2)$ である。速さの測度 $4\pi v^2\dd v$ に $\sigma(v)v$ を掛け、$E=\mu v^2/2$ に置き換えると
$$
\langle\sigma v\rangle=
\left(\frac8{\pi\mu}\right)^{1/2}(k_{\rm B}T)^{-3/2}
\int_0^\infty E\sigma(E)e^{-E/(k_{\rm B}T)}\dd E.
\label{eq-nuc-thermal-rate}
$$
次元は体積を時間で割ったものである。体積当たりの反応数は $r_{12}=n_1n_2\langle\sigma v\rangle/(1+\delta_{12})$ となる。$\delta_{12}=1$ は同種対の二重計数を除く。同じ粒子を二つ消費する反応では、粒子数の減少率は反応数の 2 倍である。

### Gamow ピーク

:term[天体物理的 S 因子]{en="astrophysical S factor" yomi="てんたいぶつりてきえすいんし"} の定義は :eqref[eq-nuc-fusion-s-factor] と同じである。$\alpha=e^2/(4\pi\varepsilon_0\hbar c)$ として $E_G=2\mu c^2(\pi\alpha Z_1Z_2)^2$ と置けば $2\pi\eta=\sqrt{E_G/E}$ となる。非共鳴で $S$ がゆっくり変化する場合、積分の主な寄与を
$$
\Phi(E)=E/(k_{\rm B}T)+\sqrt{E_G/E}
$$
の最小の近傍から求める。$\Phi'=1/(k_{\rm B}T)-\sqrt{E_G}/(2E^{3/2})$ をゼロとする解は $E_0=(E_G(k_{\rm B}T)^2/4)^{1/3}$ である。二階微分に代入すると $\Phi''(E_0)=3/(2E_0 k_{\rm B}T)$、$\Phi(E_0)=3E_0/(k_{\rm B}T)$ となる。従って Gaussian の標準偏差は $\sqrt{2E_0k_{\rm B}T/3}$ であり、最大値の $1/e$ までの全幅は $4\sqrt{E_0 k_{\rm B}T/3}$ である。温度を上げると熱分布と透過の重なるエネルギー域が上がる。

孤立した狭い共鳴では別の近似が有効である。:eqref[eq-nuc-breit-wigner] を積分に入れ、$E/k_a^2=\hbar^2/(2\mu)$ とし、幅 $\Gamma$ 内で熱因子と部分幅を一定とみなす。Lorentz 型を全実数上で積分すると $2\pi/\Gamma$ なので
$$
\langle\sigma v\rangle_R
=\left(\frac{2\pi}{\mu k_{\rm B}T}\right)^{3/2}\hbar^2
(\omega\gamma)e^{-E_R/(k_{\rm B}T)},\qquad
\omega\gamma=g_J\Gamma_a\Gamma_b/\Gamma .
$$
ここでは共鳴中心付近の寄与を求めるため、$E_R\gg\Gamma$、$\Gamma\ll k_{\rm B}T$ とした。広い共鳴や、熱重みが低エネルギー側の裾を強く拾う場合は、実際のエネルギー依存を保って積分する。

:::example[熱分布の典型値との比較]
$E_G=500\,\mathrm{keV}$、$k_{\rm B}T=1\,\mathrm{keV}$ で緩やかな $S$ なら、$E_0=125^{1/3}=5\,\mathrm{keV}$ である。反応への寄与は単に温度のエネルギー $1\,\mathrm{keV}$ に集中しない。
:::

:::example[共鳴の重み]
温度と $\omega\gamma$ が同じ二つの狭い共鳴で、後者の $E_R$ が $k_{\rm B}T\ln2$ だけ高い場合、後者の熱平均への寄与は半分である。
:::

### 温度指数を見積もる

Gamow ピークの幅が $E_0$ より十分小さく、その幅の中で $S$ を $S(E_0)$ に固定できる場合、指数の二次展開により積分を全実軸の Gaussian 積分で近似できる。:eqref[eq-nuc-thermal-rate] の積分に $S(E_0)e^{-\Phi(E_0)}\sqrt{2\pi/\Phi''(E_0)}$ を入れると

$$
\langle\sigma v\rangle\simeq S(E_0)
\sqrt{\frac{32 E_0}{3\mu}}\frac{e^{-3E_0/(k_{\rm B}T)}}{k_{\rm B}T}
$$

となる。$S(E_0)$ の温度依存を省くと、$E_0\propto T^{2/3}$ だから前因子は $T^{-2/3}$ に比例する。局所的な温度指数 $\partial\ln\langle\sigma v\rangle/\partial\ln T$ は $-2/3+E_0/(k_{\rm B}T)$ である。$S$ の傾きを残す場合には、さらに $(2/3)\,\partial\ln S/\partial\ln E$ を $E_0$ で加える。換算質量を固定して $Z_1Z_2$ を増やすと $E_G$ が増え、この指数の障壁の寄与も増える。燃焼網全体の温度指数を求めるには、共鳴や核種の存在量が温度でどう変わるかも考える。

:::example[同じ温度での比較]
$E_0/(k_{\rm B}T)=5$ と 15 の二つの非共鳴反応なら、局所指数は $13/3$ と $43/3$ である。微小な温度上昇は後者の率により大きい相対変化を与える。一方、率の絶対値には $S$ と存在量も必要である。
:::

### 狭い共鳴と反応ネットワークの単位

幅 $\Gamma$ の孤立共鳴が $E_r$ にあり、非相対論二体の断面積を $\sigma(E)=(\pi/k^2)\omega\Gamma_i\Gamma_f/[(E-E_r)^2+\Gamma^2/4]$ とする。$k^2=2\mu E/\hbar^2$ なので熱平均積分の因子 $E$ が相殺される。$\Gamma\ll E_r,k_{\rm B}T$ で幅と Boltzmann 因子を中心値に固定し、積分を全実軸に延ばすと Lorentz 形の面積は $2\pi/\Gamma$ である。従って
$$
\langle\sigma v\rangle_r=\left(\frac{2\pi}{\mu k_{\rm B}T}\right)^{3/2}\hbar^2(\omega\gamma)e^{-E_r/(k_{\rm B}T)},\qquad
\omega\gamma=\omega\,\frac{\Gamma_i\Gamma_f}{\Gamma}.
$$
これには速度と面積の積の単位がある。狭い共鳴でない場合やしきいで線が切れる場合は実際の断面積で積分する。非共鳴の反応に使う $S(E)$ の滑らかさと、狭い共鳴の近似は別の条件を課す。

組成を粒子の $Y_i=n_i/n_b$ で表し、体積が一様に変わっても $n_b$ と各 $n_i$ の希釈率が等しいとする。二体反応 $i+j$ の規格化率は $n_bY_iY_j\langle\sigma v\rangle/(1+\delta_{ij})$ であり、差の方程式には一反応で生成・消費される個数を掛ける。これにより各核種の内部状態を平均した場合でも、ネットワークの各項でバリオン数保存を確認できる。

:::example[Gamow ピークの温度依存]
$S(E)$ が一定に近い場合、指数を最小にするエネルギーは $E_0=(E_G(k_{\rm B}T)^2/4)^{1/3}$ である。

計算には次の値を用いる。
$$
E_{G}=500\,\mathrm{keV},\qquad t=1\,\mathrm{keV}
$$

この条件で求める量を ${E_0}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{E_0}&=\left(\frac{E_{G} \cdot t^{2}}{4}\right)^{1 \cdot \frac{1}{3}}\\&=\left(\frac{\left(1\right)^{2} \cdot \left(500\right)}{4}\right)^{1 \cdot \frac{1}{3}}\\&=5\,\mathrm{keV}.\end{aligned}
$$

熱平均を支配する相対エネルギーは $k_{\rm B}T$ 自身と異なる。

温度を二倍にする。次の値を変え、ほかは固定する。
$$
t=2\,\mathrm{keV}
$$

このときの値を ${E_0}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{E_0}_{\mathrm{new}}&=\left(\frac{\left(2\right)^{2} \cdot \left(500\right)}{4}\right)^{1 \cdot \frac{1}{3}}\\&\simeq 7.93701\,\mathrm{keV},\\{E_0}_{\mathrm{new}}-{E_0}&\simeq 2.93701\,\mathrm{keV}.\end{aligned}
$$
:::

:::example[狭い共鳴率の相対値]
同じ系で $\omega\gamma$ と $\mu$ を固定し、温度のエネルギー $t_1,t_2$ で比を求める。

計算には次の値を用いる。
$$
E_{r}=30\,\mathrm{keV},\qquad t_{1}=5\,\mathrm{keV},\qquad t_{2}=6\,\mathrm{keV}
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=\left(\frac{t_{1}}{t_{2}}\right)^{\frac{3}{2}} \cdot e^{E_{r} \cdot \left(- \frac{1}{t_{2}} + 1 \cdot \frac{1}{t_{1}}\right)}\\&=\left(\frac{\left(5\right)}{\left(6\right)}\right)^{\frac{3}{2}} \cdot e^{\left(30\right) \cdot \left(- \frac{1}{\left(6\right)} + \frac{1}{\left(5\right)}\right)}\\&\simeq 2.06787.\end{aligned}
$$

指数の増加には前因子の温度低下が一部対抗する。

高温の値を比べる。次の値を変え、ほかは固定する。
$$
t_{2}=8\,\mathrm{keV}
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=\left(\frac{\left(5\right)}{\left(8\right)}\right)^{\frac{3}{2}} \cdot e^{\left(30\right) \cdot \left(- \frac{1}{\left(8\right)} + \frac{1}{\left(5\right)}\right)}\\&\simeq 4.68795,\\{R}_{\mathrm{new}}-{R}&\simeq 2.62008.\end{aligned}
$$

**一次近似の検証。** 一つの温度の近傍で求めた反応率の傾きが、温度差を広げた場合にも使えるかを調べる。

${t_{2}}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${R}$ の一階微分を求めると、
$$
\frac{\mathrm d{R}}{\mathrm dt_{2}}=\frac{\left(\frac{t_{1}}{t_{2}}\right)^{\frac{3}{2}} \cdot \left(E_{r} - \frac{3 \cdot t_{2}}{2}\right) \cdot e^{- \frac{E_{r} \cdot \left(t_{1} - t_{2}\right)}{t_{1} \cdot t_{2}}}}{t_{2}^{2}}.
$$

この導関数を変更前の点で評価して $h=2$ を掛ける。ここで $t_{2}$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq 2.41251,\\\Delta-\Delta_{\mathrm{lin}}&\simeq 0.207568.\end{aligned}
$$

残差が示すのは、指定した二点間で生じる線形化の誤差である。基準点からさらに離れた値について同じ精度を保証するものではない。
:::

:::example[同種粒子による燃焼]
密度 $n$ の同種粒子の二体反応率係数を $K$ とする。体積当たりの反応数を求める。

計算には次の値を用いる。
$$
n=100000000000000000000\,\mathrm{cm}^{-3},\qquad K=\frac{1}{1000000000000000000000000}\,\mathrm{cm^3/s}
$$

この条件で求める量を ${r}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{r}&=\frac{K \cdot n^{2}}{2}\\&=\frac{\left(100000000000000000000\right)^{2} \cdot \left(\frac{1}{1000000000000000000000000}\right)}{2}\\&=5000000000000000\,\mathrm{cm}^{-3}\,\mathrm{s}^{-1}.\end{aligned}
$$

消費率はこれの二倍であり、反応数と異なる。

同じ温度で密度を二倍にする。次の値を変え、ほかは固定する。
$$
n=200000000000000000000\,\mathrm{cm}^{-3}
$$

このときの値を ${r}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{r}_{\mathrm{new}}&=\frac{\left(200000000000000000000\right)^{2} \cdot \left(\frac{1}{1000000000000000000000000}\right)}{2}\\&=20000000000000000\,\mathrm{cm}^{-3}\,\mathrm{s}^{-1},\\{r}_{\mathrm{new}}-{r}&=15000000000000000\,\mathrm{cm}^{-3}\,\mathrm{s}^{-1}.\end{aligned}
$$
:::

:::example[二体寿命の見積もり]
希薄な種が、ほぼ変わらない相手の密度 $n$ で捕獲される。係数 $K$ と時間 $t$ から残存比を得る。

計算には次の値を用いる。
$$
n=100000000000000000000\,\mathrm{cm}^{-3},\qquad K=\frac{1}{10000000000000000000000}\,\mathrm{cm^3/s},\qquad t=30\,\mathrm{s}
$$

この条件で求める量を ${f}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{f}&=e^{K \cdot \left(- n\right) \cdot t}\\&=e^{- \left(100000000000000000000\right) \cdot \left(30\right) \cdot \left(\frac{1}{10000000000000000000000}\right)}\\&\simeq 0.740818.\end{aligned}
$$

両成分が同程度に消費されると、相手を固定した指数にはならない。

燃焼時間を延ばす。次の値を変え、ほかは固定する。
$$
t=120\,\mathrm{s}
$$

このときの値を ${f}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{f}_{\mathrm{new}}&=e^{- \left(100000000000000000000\right) \cdot \left(120\right) \cdot \left(\frac{1}{10000000000000000000000}\right)}\\&\simeq 0.301194,\\{f}_{\mathrm{new}}-{f}&\simeq -0.439624.\end{aligned}
$$
:::

## 水素燃焼（pp 連鎖と CNO サイクル）

核種 $i$ の数密度を $n_i$、バリオン数密度を $n_{\rm b}=\sum_i A_i n_i$ とし、$Y_i=n_i/n_{\rm b}$ と置く。各反応 $r$ による核種 $i$ の増減数を $\nu_{ir}$、反応の体積率を $R_r$ とすると、同じ流体要素で $\dot Y_i=\sum_r \nu_{ir}R_r/n_{\rm b}$ である。反応ごとに $\sum_i A_i\nu_{ir}=0$ なら $\sum_i A_iY_i$ が保たれる。膨張で個々の密度が下がるときも、この比では共通の体積変化が消える。

### pp 連鎖

pp I の主な経路は
$$
p+p\to d+e^++\nu_e,\qquad
d+p\to{}^3\mathrm{He}+\gamma,\qquad
{}^3\mathrm{He}+{}^3\mathrm{He}\to{}^4\mathrm{He}+2p
$$
である。最初の二段を二回ずつ、最後の段を一回足すと中間の $d$ と ${}^3\mathrm{He}$ が消え、陽子が正味 4 個減る。入口の弱い相互作用が反応を遅くする。同じ重陽子を作る別の入口に $p+e^-+p\to d+\nu_e$（pep 反応）がある。

pp II と pp III は ${}^3\mathrm{He}(\alpha,\gamma){}^7\mathrm{Be}$ までが共通である。pp II では ${}^7\mathrm{Be}+e^-\to{}^7\mathrm{Li}+\nu_e$ に続き、${}^7\mathrm{Li}(p,\alpha){}^4\mathrm{He}$ により二つの $\alpha$ が生じる。pp III では ${}^7\mathrm{Be}(p,\gamma){}^8\mathrm B$ の後に、${}^8\mathrm B\to{}^8\mathrm{Be}^*+e^++\nu_e$、${}^8\mathrm{Be}^*\to2\alpha$ と進む。いずれも途中で使った $\alpha$ を差し戻すと一つが正味で増える。電子捕獲の線と $\beta^+$ の連続スペクトルの違いがニュートリノによる経路の識別に使える。

### CNO サイクル

代表的な閉じた経路は、${}^{12}\mathrm C(p,\gamma){}^{13}\mathrm N$、その $\beta^+$ で ${}^{13}\mathrm C$、次の陽子捕獲で ${}^{14}\mathrm N$、さらに陽子捕獲で ${}^{15}\mathrm O$、その $\beta^+$ で ${}^{15}\mathrm N$、最後の $(p,\alpha)$ で ${}^{12}\mathrm C$ に戻る。正味の燃料と生成物は pp と同じだが、C、N、O の核が触媒になる。重い電荷の相手との反応を含み、温度に対する感度は太陽に近い条件で pp 連鎖より大きい。どの段が律速するかは各反応率の比較で決まり、触媒の存在量も全体の率に関わる。

:::example[正味エネルギー]
原子質量 $m_{\rm H}=1.007825\,\mathrm u$、$M_{\rm at}({}^4\mathrm{He})=4.002603\,\mathrm u$ から、陽電子消滅まで含む質量差は $(4m_{\rm H}-M_{\rm at})c^2\simeq26.73\,\mathrm{MeV}$ となる。ニュートリノへ渡る分は恒星に熱として残るエネルギーから引く。
:::

:::example[最も遅い段と存在量]
二段の閉じた触媒模型で、各触媒粒子の反応率が $\lambda_1,\lambda_2$ なら定常時に $\lambda_1Y_1=\lambda_2Y_2$ である。$\lambda_2=\lambda_1/10$ なら $Y_2=10Y_1$ となり、遅い段の前に触媒が多くたまる。
:::

### pp 分岐の収支と CNO の定常流

pp 連鎖の各枝では、陽子四つから正味で $\alpha$ 核一つができる。pp 反応を入口とする場合、pp I と pp III では陽電子が二つ生じる。pp II では一つの陽電子放出と一つの電子捕獲があり、pep 反応を通る場合も陽電子放出の一つが電子捕獲に置き換わる。陽電子の消滅まで数えると、いずれも $4p+2e^-\to\alpha+2\nu_e$ が正味の粒子収支である。中性原子質量で求めた $Q$ はこの収支に対応する。枝が変わるとニュートリノが持ち去るエネルギーも変わるため、同じ $\alpha$ の生成率でも局所発熱は異なる。

CNO の主経路では ${}^{14}\mathrm N(p,\gamma){}^{15}\mathrm O$ が遅い。この段が他より十分遅い定常模型では、触媒の多くが ${}^{14}\mathrm N$ として存在する。さらに ${}^{15}\mathrm N$ の $(p,\gamma)$ 分岐は ${}^{16}\mathrm O$ を作る。この核は陽子捕獲で ${}^{17}\mathrm F$ となり、その $\beta^+$ 崩壊で ${}^{17}\mathrm O$ ができる。続く $(p,\alpha)$ で ${}^{14}\mathrm N$ に戻るので、主経路の一部を置き換える閉じた枝になる。複数の枝を同時に扱う場合は各接続点で流入と流出を足し合わせる。

閉じた理想的な CNO 循環で触媒種 $i$ が次へ移る率を $\lambda_i$ とする。定常では $\lambda_{i-1}n_{i-1}-\lambda_i n_i=0$ なので全リンクの流量は同じ $J$ である。総触媒密度 $n_{\rm cat}=\sum_i n_i$ を固定すると $J=n_{\rm cat}/\sum_i\lambda_i^{-1}$ となる。この証明は各リンクが一方向で、枝分かれや外への漏れが無視でき、触媒組成が定常に達した場合に限る。律速が一つならそのリンクの前に触媒が溜まる。

:::example[水素とヘリウムの質量差]
中性の水素とヘリウムの質量を原子質量単位で扱い、その単位を $u_c$ MeV とする。

計算には次の値を用いる。
$$
m_{H}=\frac{40313}{40000},\qquad m_{He}=\frac{4002603}{1000000},\qquad u_{c}=\frac{465747}{500}
$$

この条件で求める量を ${Q}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{Q}&=u_{c} \cdot \left(4 \cdot m_{H} - m_{He}\right)\\&=\left(\frac{465747}{500}\right) \cdot \left(- \left(\frac{4002603}{1000000}\right) + 4 \cdot \left(\frac{40313}{40000}\right)\right)\\&\simeq 26.7311\,\mathrm{MeV}.\end{aligned}
$$

第二の値は換算を概算にした結果である。

換算係数の丸めによる影響を比べる。次の値を変え、ほかは固定する。
$$
u_{c}=932
$$

このときの値を ${Q}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{Q}_{\mathrm{new}}&=\left(932\right) \cdot \left(- \left(\frac{4002603}{1000000}\right) + 4 \cdot \left(\frac{40313}{40000}\right)\right)\\&\simeq 26.7456\,\mathrm{MeV},\\{Q}_{\mathrm{new}}-{Q}&\simeq 0.0145207\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[ニュートリノが持ち去る分]
一回の完成鎖での中性原子の窓を $Q$、ニュートリノの合計エネルギーを $E_\nu$ とする。

計算には次の値を用いる。
$$
Q=\frac{2673}{100}\,\mathrm{MeV},\qquad E_{\nu}=\frac{3}{5}\,\mathrm{MeV}
$$

この条件で求める量を ${Q_h}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{Q_h}&=- E_{\nu} + Q\\&=\left(\frac{2673}{100}\right) - \left(\frac{3}{5}\right)\\&\simeq 26.13\,\mathrm{MeV}.\end{aligned}
$$

光子と荷電粒子のエネルギーが局所的に熱化すると仮定する。

大きなニュートリノ損失と比較する。次の値を変え、ほかは固定する。
$$
E_{\nu}=2\,\mathrm{MeV}
$$

このときの値を ${Q_h}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{Q_h}_{\mathrm{new}}&=- \left(2\right) + \left(\frac{2673}{100}\right)\\&\simeq 24.73\,\mathrm{MeV},\\{Q_h}_{\mathrm{new}}-{Q_h}&\simeq -1.4\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[二経路の局所競合]
ある中間核の二つの消失経路を $l_1,l_2$ とする。待ち時間を積分して第二への割合を求める。

計算には次の値を用いる。
$$
l_{1}=4\,\mathrm{s}^{-1},\qquad l_{2}=1\,\mathrm{s}^{-1}
$$

この条件で求める量を ${b_2}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{b_2}&=\frac{l_{2}}{l_{1} + l_{2}}\\&=\frac{\left(1\right)}{\left(1\right) + \left(4\right)}\\&\simeq 0.2.\end{aligned}
$$

これは一つの中間核の分岐であり、全連鎖の絶対フラックスは別の式で決まる。

第二の相手密度が増す場合を比較する。次の値を変え、ほかは固定する。
$$
l_{2}=2\,\mathrm{s}^{-1}
$$

このときの値を ${b_2}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{b_2}_{\mathrm{new}}&=\frac{\left(2\right)}{\left(2\right) + \left(4\right)}\\&\simeq 0.333333,\\{b_2}_{\mathrm{new}}-{b_2}&\simeq 0.133333.\end{aligned}
$$
:::

:::example[触媒の定常分配]
二つに縮約したサイクルの個々の滞在時間を $t_1,t_2$ とする。第一の占有割合を求める。

計算には次の値を用いる。
$$
t_{1}=8\,\mathrm{s},\qquad t_{2}=2\,\mathrm{s}
$$

この条件で求める量を ${f_1}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{f_1}&=\frac{t_{1}}{t_{1} + t_{2}}\\&=\frac{\left(8\right)}{\left(2\right) + \left(8\right)}\\&\simeq 0.8.\end{aligned}
$$

両時間が等しければ一方への蓄積はない。

速い側を遅くする。次の値を変え、ほかは固定する。
$$
t_{2}=8\,\mathrm{s}
$$

このときの値を ${f_1}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{f_1}_{\mathrm{new}}&=\frac{1}{2}\\&\simeq 0.5,\\{f_1}_{\mathrm{new}}-{f_1}&\simeq -0.3.\end{aligned}
$$
:::

## ヘリウム燃焼とそれ以降の核燃焼

四つの核子を持つ $\alpha$ 粒子二つからできる ${}^8\mathrm{Be}$ は束縛されず短寿命である。しかし熱平衡に近い系では一部が一時的に存在し、さらに $\alpha$ を捕獲できる。炭素の励起状態のうち、三 $\alpha$ のしきい値近傍にある Hoyle 状態が率を増す。ここから光子遷移で ${}^{12}\mathrm C$ の安定な準位に移る過程を :term[三重アルファ反応]{en="triple-alpha process" yomi="さんじゅうあるふぁはんのう"} という。

固定温度の簡単な逐次模型で $\dot n_8=k_1n_\alpha^2/2-\lambda_8n_8-k_2n_8n_\alpha$ とする。$n_8$ が速く緩和し、$\lambda_8\gg k_2n_\alpha$ なら $n_8\simeq k_1n_\alpha^2/(2\lambda_8)$ である。炭素への反応率 $k_2n_8n_\alpha$ はしたがって $n_\alpha^3$ に比例する。$k_1,k_2$ は熱平均二体反応率係数、$\lambda_8$ は崩壊定数である。共鳴準位の量子的な寄与の評価には三体と核準位の扱いが必要であるため :cite[krane-nuclear] を参照する。

生成された炭素は ${}^{12}\mathrm C(\alpha,\gamma){}^{16}\mathrm O$ で酸素を作る。ヘリウムが尽きる時点の炭素と酸素の比は、この二つの経路の競合で変わる。十分重い星ではさらに炭素、ネオン、酸素、ケイ素燃焼が進む。高温のケイ素燃焼では捕獲と光分解が互いに進み、反応群が平衡に近づく。鉄族付近では結合エネルギーの一核子当たりの増分が小さく、より重い核へ融合して長く熱を供給することは難しくなる。

燃焼段階の呼び名は、最初に消費する核と主要な流れを示す。炭素燃焼での代表的な入口は ${}^{12}\mathrm C+{}^{12}\mathrm C$ であり、${}^{20}\mathrm{Ne}+\alpha$ と ${}^{23}\mathrm{Na}+p$ の出口などがある。ネオン燃焼では、高エネルギー光子による ${}^{20}\mathrm{Ne}(\gamma,\alpha){}^{16}\mathrm O$ が $\alpha$ を出し、それを別の ${}^{20}\mathrm{Ne}$ が捕獲して ${}^{24}\mathrm{Mg}$ を作る。二段を足すと $\alpha$ は相殺する。酸素燃焼では ${}^{16}\mathrm O+{}^{16}\mathrm O\to{}^{28}\mathrm{Si}+\alpha$ などでケイ素族を作る。次のケイ素燃焼では、光分解で出た軽粒子の捕獲と逆反応が速く、反応網の一部が準平衡に達する。生成する鉄族核の内訳には全バリオンに占める陽子数の割合が効く。

:::example[密度依存の比較]
温度とヘリウム分率を固定し、$n_\alpha$ を 2 倍にすると上の三重反応の体積率は 8 倍になる。単純な二体反応の体積率は 4 倍である。質量当たりの値に直すと密度で割るため、倍率はそれぞれ 4 と 2 になる。
:::

:::example[中間体の補正]
$k_2n_\alpha=\lambda_8$ なら、定常解の分母は $2\lambda_8$ となる。従って $k_2n_\alpha$ を無視した場合の $n_8$ の半分である。逆に $k_2n_\alpha\gg\lambda_8$ なら、炭素の生成率は ${}^8\mathrm{Be}$ の生成率 $k_1n_\alpha^2/2$ に近づく。
:::

### 三重 α の熱平衡近似と燃焼の競合

短寿命の $^8\mathrm{Be}$ が $\alpha+\alpha$ と速い準平衡にあるとする。希薄な Maxwell–Boltzmann 気体の密度を $n_i=g_i(m_i k_{\rm B}T/2\pi\hbar^2)^{3/2}\exp[(\mu_i-m_i c^2)/k_{\rm B}T]$ と書くと、$\mu_8=2\mu_\alpha$ により $n_8/n_\alpha^2$ の依存は $T^{-3/2}\exp[-(m_8-2m_\alpha)c^2/(k_{\rm B}T)]$ となる。次の共鳴の狭幅近似で、もう一つの $T^{-3/2}$ と $\exp(-E_r/k_{\rm B}T)$ が掛かる。指数を合わせると三つの $\alpha$ のしきいから Hoyle 状態までの差が残り、全生成率は $n_\alpha^3 T^{-3}\exp[-E_H/(k_{\rm B}T)]$ に比例する。ここでは共鳴幅の温度依存と別経路を無視している。最終的な炭素と酸素の比には炭素を消費する $\alpha$ 捕獲も関わるので、この生成率だけでは決まらない。

:::example[定常な中間核の損失]
作られた中間体が率 $\lambda$ で戻るか、擬一次率 $u=k_2n_\alpha$ で炭素へ進む。

計算には次の値を用いる。
$$
lam=100\,\mathrm{s}^{-1},\qquad u=1\,\mathrm{s}^{-1}
$$

この条件で求める量を ${b}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{b}&=\frac{u}{lam + u}\\&=\frac{\left(1\right)}{\left(100\right) + \left(1\right)}\\&\simeq 0.00990099.\end{aligned}
$$

中間核を経由した後の条件付き確率である。

捕獲に使う相手密度の効果を比べる。次の値を変え、ほかは固定する。
$$
u=5\,\mathrm{s}^{-1}
$$

このときの値を ${b}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{b}_{\mathrm{new}}&=\frac{\left(5\right)}{\left(100\right) + \left(5\right)}\\&\simeq 0.047619,\\{b}_{\mathrm{new}}-{b}&\simeq 0.0377181.\end{aligned}
$$
:::

:::example[温度に対する局所指数]
$E_H$ を一定にして $T^{-3}e^{-E_H/(k_{\rm B}T)}$ を $\ln T$ で微分する。

計算には次の値を用いる。
$$
E_{H}=379\,\mathrm{keV},\qquad t=8\,\mathrm{keV}
$$

この条件で求める量を ${\nu}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\nu}&=\frac{E_{H}}{t} - 3\\&=\frac{\left(379\right)}{\left(8\right)} - 3\\&\simeq 44.375.\end{aligned}
$$

温度を変えると局所のべき指数自体が変わる。

より高い温度を比較する。次の値を変え、ほかは固定する。
$$
t=15\,\mathrm{keV}
$$

このときの値を ${\nu}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\nu}_{\mathrm{new}}&=-3 + \frac{\left(379\right)}{\left(15\right)}\\&\simeq 22.2667,\\{\nu}_{\mathrm{new}}-{\nu}&\simeq -22.1083.\end{aligned}
$$
:::

:::example[密度とヘリウム組成]
体積当たりの三体の生成数を考える。二つの状態での $\alpha$ 数密度比を $r$ とする。

計算には次の値を用いる。
$$
r=\frac{3}{2}
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=r^{3}\\&=\left(\frac{3}{2}\right)^{3}\\&\simeq 3.375.\end{aligned}
$$

質量当たりの生成率を求めるなら各結果を質量密度で割る。

密度比を大きくする。次の値を変え、ほかは固定する。
$$
r=2
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=\left(2\right)^{3}\\&=8,\\{R}_{\mathrm{new}}-{R}&\simeq 4.625.\end{aligned}
$$
:::

:::example[捕獲による炭素の残存]
酸素へ進む捕獲の一次率を $\lambda$ とし、新たな炭素生成を止めた後を考える。

計算には次の値を用いる。
$$
lam=\frac{1}{500}\,\mathrm{s}^{-1},\qquad t=100\,\mathrm{s}
$$

この条件で求める量を ${f}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{f}&=e^{- lam \cdot t}\\&=e^{- \left(100\right) \cdot \left(\frac{1}{500}\right)}\\&\simeq 0.818731.\end{aligned}
$$

進行中のヘリウム燃焼では供給項も同じ式に入る。

時間を延長する。次の値を変え、ほかは固定する。
$$
t=500\,\mathrm{s}
$$

このときの値を ${f}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{f}_{\mathrm{new}}&=e^{- \left(500\right) \cdot \left(\frac{1}{500}\right)}\\&\simeq 0.367879,\\{f}_{\mathrm{new}}-{f}&\simeq -0.450851.\end{aligned}
$$
:::

## s 過程と r 過程

中性子は Coulomb 障壁を持たず、捕獲で質量数 $A$ を一つ上げる。中性子過剰になった核は $\beta^-$ で $Z$ を上げる。捕獲の単位時間確率を $\lambda_n=n_n\langle\sigma v\rangle$、崩壊定数を $\lambda_\beta$ とすると、競合する二指数過程で捕獲が先に起こる確率は $\lambda_n/(\lambda_n+\lambda_\beta)$ である。これは :eqref[eq-nuc-decay-law] の分岐の考え方と同じである。

:term[s 過程]{en="s-process" yomi="えすかてい"} では不安定核の多くが次の捕獲までに崩壊でき、経路は安定線に近い。分岐のない定常な系列では $\dot Y_A=\lambda_{A-1}Y_{A-1}-\lambda_A Y_A=0$ となる。従って $\lambda_A Y_A$ が一定で、捕獲率の小さい閉殻付近では存在量が溜まりやすい。分岐する核種ではその量の比から中性子密度や温度への制約を得られる。

:term[r 過程]{en="r-process" yomi="あーるかてい"} では捕獲が $\beta$ 崩壊より速く、経路は中性子過剰側へ進む。高温では中性子捕獲と光分解が同位体の列で釣り合うことがある。化学ポテンシャルの条件 $\mu_{Z,A+1}=\mu_{Z,A}+\mu_n$ に古典気体の $n_i=g_i(m_i k_{\rm B}T/2\pi\hbar^2)^{3/2}e^{(\mu_i-m_ic^2)/(k_{\rm B}T)}$ を入れると、隣の同位体の比は $n_n$ と $e^{S_n/(k_{\rm B}T)}$ に比例する。$S_n=(m_{Z,A}+m_n-m_{Z,A+1})c^2$ であるから、核質量が経路に効く。中性子照射が終わると、$\beta$ 崩壊や遅発中性子放出で安定側に近づく。中性子星合体はこの過程を考える重要な環境であるが、全元素の起源を断定するには生成物の質量や放出頻度も必要である。

:::example[分岐比]
$\lambda_n=2\,\mathrm{s}^{-1}$、$\lambda_\beta=8\,\mathrm{s}^{-1}$ なら捕獲が先の確率は $0.2$ である。密度だけ 4 倍になれば $\lambda_n=8\,\mathrm{s}^{-1}$ となり、確率は $1/2$ になる。
:::

:::example[定常流と閉殻]
反応流が同じで核種 $A$ の捕獲率が隣の 1/5 なら $Y_A$ は隣の 5 倍になる。燃料が尽きて定常でなくなると、この関係をそのまま最終存在量へ使えない。
:::

### 捕獲と光分解の比を求める

一般の $i+n\rightleftharpoons j+\gamma$ で、核種 $i$ と $j$ のスピン縮退度と内部励起をまとめた無次元の分配関数を $G_i(T),G_j(T)$ とする。密度を古典気体として $\mu_j=\mu_i+\mu_n$ を使えば
$$
\frac{n_j}{n_i n_n}
=\frac{G_j}{G_i G_n}
 \left(\frac{m_j}{m_i m_n}\right)^{3/2}
 \left(\frac{2\pi\hbar^2}{k_{\rm B}T}\right)^{3/2}
 e^{S_n/(k_{\rm B}T)}
$$
である。これは $j$ の静止質量の差が指数 $S_n$ を与え、化学ポテンシャルが消えることによる。熱平衡で $n_in_n\langle\sigma v\rangle_{in}=n_j\lambda_{\gamma j}$ なので、逆反応率は右辺の逆数に $\langle\sigma v\rangle_{in}$ を掛けたものになる。非平衡の反応網でも、各反応率が同じ温度の熱分布で計算されていればこの詳細釣合いによる率の関係を使える。実際に存在量の比が平衡になるかは反応の時定数との比較が要る。

:::example[分離エネルギーの差の感度]
$k_{\rm B}T=0.1\,\mathrm{MeV}$ で、統計重みと質量の前因子を共通とする。$S_n$ が $0.2\,\mathrm{MeV}$ 大きい場合、同じ中性子密度で平衡の隣接比は $e^2\simeq7.39$ 倍になる。質量測定が存在量の予測に重要になるしくみである。
:::

### 中性子捕獲流と分岐点

$s$ 過程で $\beta$ 崩壊が捕獲より十分速く、各質量数を一つの同位体で代表できる区間を考える。$n_n v$ を共通の流束 $\Phi$ とし、有効断面積を $\sigma_A$ とする。生成消滅式 $\dot Y_A=\Phi(\sigma_{A-1}Y_{A-1}-\sigma_AY_A)$ が定常なら $\sigma_A Y_A$ は一定になる。小さな断面積の場所には核が溜まる。ある核で捕獲と $\beta$ が競合すると、流量が二つの道に分かれてこの単一路の関係が崩れる。

$r$ 過程のある同位体列内で捕獲と光分解がともに外部進化や $\beta$ より速いとき、化学平衡の比が使える。中性子分離エネルギーが指数に入るので、未知の核質量が平衡経路に効く。捕獲と光分解がともに遅くなった後は、平衡比を固定して進化を代用できず、残存中性子と $\beta$ 等を時間発展で数える。

恒星の $s$ 過程では ${}^{13}\mathrm C(\alpha,n){}^{16}\mathrm O$ と ${}^{22}\mathrm{Ne}(\alpha,n){}^{25}\mathrm{Mg}$ が代表的な中性子源である。中性子を放出する割合と照射時間が異なると、分岐点を通る割合も変わる。閉殻の中性子数 $N=50,82,126$ 付近では捕獲の断面積が小さくなる核があり、流れの中で存在量が増えやすい。一方、$r$ 過程では中性子過剰側の閉殻に近い同位体で次の元素への $\beta$ 崩壊を待つことがある。その後に $\beta$ 崩壊で安定側へ移るので、同じ中性子閉殻でも $s$ 過程と $r$ 過程でピークの質量数は一般に異なる。遅発中性子放出や凍結後の捕獲は、この単純な対応からピークをずらす。

:::example[分岐点の捕獲確率]
捕獲率を $\lambda_n$、半減期 $T_\beta$ から $\lambda_\beta=\ln2/T_\beta$ とする。

計算には次の値を用いる。
$$
l_{n}=\frac{1}{100}\,\mathrm{s}^{-1},\qquad T_{b}=60\,\mathrm{s}
$$

この条件で求める量を ${b_n}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{b_n}&=\frac{l_{n}}{l_{n} + \frac{\log{\left(2 \right)}}{T_{b}}}\\&=\frac{\left(\frac{1}{100}\right)}{\left(\frac{1}{100}\right) + \frac{\log{\left(2 \right)}}{\left(60\right)}}\\&\simeq 0.463984.\end{aligned}
$$

率が同程度のときに分岐が環境の変化に敏感になる。

捕獲の強い条件と比較する。次の値を変え、ほかは固定する。
$$
l_{n}=\frac{3}{100}\,\mathrm{s}^{-1}
$$

このときの値を ${b_n}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{b_n}_{\mathrm{new}}&=\frac{\left(\frac{3}{100}\right)}{\left(\frac{3}{100}\right) + \frac{\log{\left(2 \right)}}{\left(60\right)}}\\&\simeq 0.721979,\\{b_n}_{\mathrm{new}}-{b_n}&\simeq 0.257995.\end{aligned}
$$
:::

:::example[定常流での濃度比]
隣り合う二核の断面積 $\sigma_1,\sigma_2$ を使い $Y_2/Y_1$ を計算する。

計算には次の値を用いる。
$$
s_{1}=80\,\mathrm{mb},\qquad s_{2}=20\,\mathrm{mb}
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=\frac{s_{1}}{s_{2}}\\&=\frac{\left(80\right)}{\left(20\right)}\\&=4.\end{aligned}
$$

小さい断面積は同じ流れにより大きい存在量を必要とする。

閉殻で断面積がさらに小さい模型を比べる。次の値を変え、ほかは固定する。
$$
s_{2}=5\,\mathrm{mb}
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=\frac{\left(80\right)}{\left(5\right)}\\&=16,\\{R}_{\mathrm{new}}-{R}&=12.\end{aligned}
$$
:::

:::example[平衡比の質量感度]
他の因子を固定し、分離エネルギーが $\Delta S$ 増えた場合の平衡比の倍率を求める。

計算には次の値を用いる。
$$
dS=\frac{1}{5}\,\mathrm{MeV},\qquad t=\frac{1}{10}\,\mathrm{MeV}
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=e^{\frac{dS}{t}}\\&=e^{\frac{\left(\frac{1}{5}\right)}{\left(\frac{1}{10}\right)}}\\&\simeq 7.38906.\end{aligned}
$$

核質量の差が温度のエネルギーより大きいと比が強く変わる。

質量評価の違いが大きい場合を比較する。次の値を変え、ほかは固定する。
$$
dS=\frac{1}{2}\,\mathrm{MeV}
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=e^{\frac{\left(\frac{1}{2}\right)}{\left(\frac{1}{10}\right)}}\\&\simeq 148.413,\\{R}_{\mathrm{new}}-{R}&\simeq 141.024.\end{aligned}
$$
:::

:::example[中性子照射量の規格化]
積算流束 $F$ と一定断面積 $\sigma$ を用い、捕獲前に残る割合を求める。

計算には次の値を用いる。
$$
F=1000000000000000000000000\,\mathrm{cm}^{-2},\qquad \sigma=\frac{1}{2000000000000000000000000}\,\mathrm{cm^2}
$$

この条件で求める量を ${f}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{f}&=e^{- F \cdot \sigma}\\&=e^{- \left(1000000000000000000000000\right) \cdot \left(\frac{1}{2000000000000000000000000}\right)}\\&\simeq 0.606531.\end{aligned}
$$

捕獲後の核種の追加捕獲は含めていない。

積算流束を増やす。次の値を変え、ほかは固定する。
$$
F=3000000000000000000000000\,\mathrm{cm}^{-2}
$$

このときの値を ${f}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{f}_{\mathrm{new}}&=e^{- \left(3000000000000000000000000\right) \cdot \left(\frac{1}{2000000000000000000000000}\right)}\\&\simeq 0.22313,\\{f}_{\mathrm{new}}-{f}&\simeq -0.3834.\end{aligned}
$$
:::

### 分離エネルギーによる待機点の目安

同位体列内の平衡比を $Y_{A+1}/Y_A=C(T)n_n\exp[S_n(A+1)/(k_{\rm B}T)]$ とまとめる。$C(T)$ は熱波長比と内部状態の重みを含む。比が 1 となる条件は $S_n(A+1)=-k_{\rm B}T\ln[C(T)n_n]$ である。同じ温度で $n_n$ を増やすと右辺は小さくなるため、通常は分離エネルギーの小さい中性子過剰側の核が占有されやすい。ただし奇偶で $S_n$ の変化が跳ぶ場合は、最も多い核を一つの連続式で選べない。その核から次の元素へ進む時間は $\beta$ 崩壊の率で決まるため、平衡比を決める速い反応と元素を変える遅い反応を区別する。

## ビッグバン元素合成

宇宙初期の温度が高い間は中性子と陽子が弱い過程で変換する。レプトン化学ポテンシャルを無視した熱平衡では、質量差 $\Delta=(m_n-m_p)c^2\simeq1.293\,\mathrm{MeV}$ に対して $n_n/n_p\simeq e^{-\Delta/(k_{\rm B}T)}$ である。両者のスピン縮退度は同じで、質量の $3/2$ 乗の比は 1 に近い。放射優勢な宇宙では膨張率 $H$ は温度の二乗にほぼ比例する。これは Friedmann 方程式 $H^2\propto$ エネルギー密度と熱放射の $T^4$ から従う。弱い変換の率が $H$ より小さくなると平衡比を追えなくなる。その後も自由中性子は崩壊する。

重陽子の結合エネルギーは約 $2.22\,\mathrm{MeV}$ だが、バリオン数に比べて光子数が多いので、温度がそれより低くても光分解する高エネルギー光子が残る。光分解が十分に減ると重陽子が蓄積し、ヘリウムなどの軽核ができる。この遅れを重陽子のボトルネックという。核種ごとの存在量の予測には膨張中の逆反応と分岐を数値積分する必要がある。冷却が速いことに加え、質量数 5 と 8 に安定核がないことが重い元素の効率的な生成を妨げる。

重陽子が残るようになると、$d(d,p){}^3\mathrm H$ と $d(d,n){}^3\mathrm{He}$ が進み、続く重陽子捕獲で ${}^4\mathrm{He}$ ができる。少量の ${}^7\mathrm{Li}$ と、その後に電子捕獲で ${}^7\mathrm{Li}$ になる ${}^7\mathrm{Be}$ も生じる。標準的な計算の比較では、残存重水素はバリオン密度を上げると後続の反応でより消費される。一方、${}^4\mathrm{He}$ の分率は主に合成開始時の中性子数を反映し、その数には放射の自由度を通じて膨張率が効く。このため重水素量とヘリウム量は異なる入力を検査できる。

:::example[ヘリウムの質量分率]
合成に入る時点で $r=n_n/n_p=1/7$ とし、残存中性子をすべて ${}^4\mathrm{He}$ に入れる。ヘリウム個数は $n_n/2$ なので、バリオン数で測る分率は $4(n_n/2)/(n_n+n_p)=2r/(1+r)=1/4$ である。
:::

:::example[平衡比と温度]
$k_{\rm B}T=1\,\mathrm{MeV}$ では平衡近似で中性子対陽子の比は $e^{-1.293}\simeq0.274$ である。$k_{\rm B}T$ を半分にすると指数が 2 倍になり、比は約 $0.075$ となる。凍結後は、下がる温度をこの式に代入して比を予測することはできない。
:::

### 重陽子が残る条件の平衡近似

光子の化学ポテンシャルをゼロとし、$n+p\rightleftharpoons d+\gamma$ の平衡条件 $\mu_d=\mu_n+\mu_p$ を使う。核種 $i$ の熱波長を $\lambda_i=(2\pi\hbar^2/(m_i k_{\rm B}T))^{1/2}$、スピン縮退度を $g_i$ とすると、古典気体の密度は $n_i=(g_i/\lambda_i^3)e^{(\mu_i-m_ic^2)/(k_{\rm B}T)}$ である。代入して
$$
\frac{n_d}{n_n n_p}
=\frac{g_d}{g_n g_p}
 \left(\frac{m_d}{m_n m_p}\right)^{3/2}
 \left(\frac{2\pi\hbar^2}{k_{\rm B}T}\right)^{3/2}
 e^{B_d/(k_{\rm B}T)}
$$
を得る。$B_d=(m_n+m_p-m_d)c^2$ であり、指数の符号は正である。重陽子の分率を求めるには両辺に $n_n n_p/n_{\rm b}$ を掛けるので、固定温度ではバリオン密度が低いほど平衡での分率も低い。逆に温度が下がると指数により束縛状態が有利になる。平衡を追えるかどうかは膨張率に対する個別の反応率で検査する。

:::example[縮退度の因子]
中性子と陽子はともにスピン $1/2$、重陽子は $J=1$ である。内部励起を無視すると $g_d/(g_n g_p)=3/4$ である。中性子と陽子は区別できるため、平衡条件から得た式に同種二体反応の $1/2$ は掛からない。
:::

### 膨張と弱反応の時間尺度

放射優勢で相対論的な自由度 $g_*$ を一定とすると、放射のエネルギー密度は $g_*T^4$ に比例し、Friedmann 方程式から $H\propto\sqrt{g_*}T^2$ となる。低エネルギーの四点弱相互作用で相対論的レプトンの密度が $T^3$、断面積が $T^2$ の尺度をもつと $\Gamma_{\rm weak}\propto T^5$ である。従って $\Gamma_{\rm weak}/H\propto T^3/\sqrt{g_*}$ が冷却とともに低下する。弱反応率を固定した比較で $g_*$ を増やすと同じ温度で比が下がるので、弱反応はより高い温度から平衡を追いにくくなる。凍結は一時刻で切り替わる現象ではなく、その後の自由中性子崩壊や電子・陽電子対の消滅に伴う $g_*$ の変化も含めて存在量を求める。

すべての残存中性子がほぼ $^4\mathrm{He}$ に入る近似では、$n_{\rm He}=n_n/2$ である。全バリオンに占める比は $Y_p^{({\rm mass})}\simeq4n_{\rm He}/(n_n+n_p)=2r/(1+r)$、ただし $r=n_n/n_p$ である。ここで $Y_p^{({\rm mass})}$ の添字は primordial の意味で、陽子の個数比ではない。この近似は残存重水素・軽元素と質量欠損の寄与を省いている。

:::example[凍結温度の前後]
基準温度 $T_f$ で弱反応と膨張が同速とし、温度比を $u$ とする。

計算には次の値を用いる。
$$
u=\frac{3}{2}
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=u^{3}\\&=\left(\frac{3}{2}\right)^{3}\\&\simeq 3.375.\end{aligned}
$$

低い温度では膨張の時間より反応の時間が長くなる。

凍結の前後を比べる。次の値を変え、ほかは固定する。
$$
u=\frac{7}{10}
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=\left(\frac{7}{10}\right)^{3}\\&\simeq 0.343,\\{R}_{\mathrm{new}}-{R}&\simeq -3.032.\end{aligned}
$$
:::

:::example[中性子の自由崩壊と個数比]
初期比を $r_0=n_n/n_p$ とし、中性子の平均寿命を $\tau$ とする。崩壊した中性子は陽子に加わるとし、時刻 $t$ の比を求める。

計算には次の値を用いる。
$$
r_{0}=\frac{1}{5},\qquad t=100\,\mathrm{s},\qquad \tau=880\,\mathrm{s}
$$

この条件で求める量を ${r}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{r}&=\frac{r_{0} \cdot e^{\frac{\left(-1\right) \cdot t}{\tau}}}{r_{0} \cdot \left(1 - e^{\frac{\left(-1\right) \cdot t}{\tau}}\right) + 1}\\&=\frac{\left(\frac{1}{5}\right) \cdot e^{- \frac{\left(100\right)}{\left(880\right)}}}{\left(\frac{1}{5}\right) \cdot \left(1 - e^{- \frac{\left(100\right)}{\left(880\right)}}\right) + 1}\\&\simeq 0.174762.\end{aligned}
$$

分母の陽子の増加も保存則から入れている。

核反応までの待ち時間を延ばす。次の値を変え、ほかは固定する。
$$
t=200\,\mathrm{s}
$$

このときの値を ${r}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{r}_{\mathrm{new}}&=\frac{\left(\frac{1}{5}\right) \cdot e^{- \frac{\left(200\right)}{\left(880\right)}}}{\left(\frac{1}{5}\right) \cdot \left(1 - e^{- \frac{\left(200\right)}{\left(880\right)}}\right) + 1}\\&\simeq 0.153115,\\{r}_{\mathrm{new}}-{r}&\simeq -0.0216469.\end{aligned}
$$
:::

:::example[予想されるヘリウムの質量比]
合成開始時の中性子と陽子の比を $r$ とする。

計算には次の値を用いる。
$$
r=\frac{1}{7}
$$

この条件で求める量を ${Y}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{Y}&=\frac{2 \cdot r}{r + 1}\\&=\frac{2 \cdot \left(\frac{1}{7}\right)}{\left(\frac{1}{7}\right) + 1}\\&\simeq 0.25.\end{aligned}
$$

中性子がより多く残れば、近似でのヘリウム質量比も増える。

別の残存比を比較する。次の値を変え、ほかは固定する。
$$
r=\frac{1}{6}
$$

このときの値を ${Y}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{Y}_{\mathrm{new}}&=\frac{2 \cdot \left(\frac{1}{6}\right)}{\left(\frac{1}{6}\right) + 1}\\&\simeq 0.285714,\\{Y}_{\mathrm{new}}-{Y}&\simeq 0.0357143.\end{aligned}
$$
:::

:::example[光子の尾部に対する尺度]
温度をエネルギー $t$ とし、束縛 $B$ のしきいを越える光子の指数因子を $e^{-B/t}$ で見積もる。

計算には次の値を用いる。
$$
B=\frac{11123}{5000}\,\mathrm{MeV},\qquad t=\frac{1}{10}\,\mathrm{MeV}
$$

この条件で求める量を ${f}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{f}&=e^{\frac{\left(-1\right) \cdot B}{t}}\\&=e^{- \frac{\left(\frac{11123}{5000}\right)}{\left(\frac{1}{10}\right)}}\\&\simeq 2.18115\times10^{-10}.\end{aligned}
$$

これは光子の高い側の指数部分だけであり、実際の分解率には光子数と断面積も入る。

核分解が減る向きの冷却を考える。次の値を変え、ほかは固定する。
$$
t=\frac{2}{25}\,\mathrm{MeV}
$$

このときの値を ${f}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{f}_{\mathrm{new}}&=e^{- \frac{\left(\frac{11123}{5000}\right)}{\left(\frac{2}{25}\right)}}\\&\simeq 8.38216\times10^{-13},\\{f}_{\mathrm{new}}-{f}&\simeq -2.17276\times10^{-10}.\end{aligned}
$$
:::

### 光子の高エネルギー尾とバリオンの希薄さ

温度のエネルギーを $\vartheta$、重陽子の光分解しきいを $B$ とする。黒体光子のうち $E\ge B$ となる割合は、エネルギーを $x=E/\vartheta$ で無次元化して
$$
f_\gamma=\frac{1}{2\zeta(3)}\int_b^\infty\frac{x^2}{e^x-1}\dd x,\qquad b=B/\vartheta
$$
と書ける。全積分は $1/(e^x-1)=\sum_{n\ge1}e^{-nx}$ を項別積分して $2\sum_{n\ge1}n^{-3}=2\zeta(3)$ である。$b\gg1$ では分子の積分関数を $x^2e^{-x}$ で置き換え、部分積分を二回すると
$$
f_\gamma\simeq\frac{(b^2+2b+2)e^{-b}}{2\zeta(3)}
$$
となる。$x\ge b$ で近似した分母の補正は $1/(1-e^{-x})$ なので、相対誤差の上限には $e^{-b}/(1-e^{-b})$ を使える。しきいを越える光子数を一バリオン当たりに換算するには、この割合をバリオン・光子数比 $\eta$ で割る。指数因子が小さくても、$\eta$ 自体が小さければ一バリオン当たりの光子数は大きい。これは重陽子を遅い時刻まで光分解しうる理由を示す。実際の凍結の温度を得るには、しきい以上の光子の数だけでなく光分解の断面積、逆向きの捕獲、膨張時間を比較する。

## 超新星と中性子星物質

重い星の中心で鉄族の芯が成長すると、電子の縮退圧が重力に抗する。電子捕獲で電子の割合が下がり、光分解で熱の一部が結合の解除に費やされる。芯が崩壊して高い密度に達すると核物質の応答が重要になる。衝撃波とニュートリノ輸送の組合せが、外層の放出と元素合成の条件を変える。爆発的な加熱では捕獲と光分解が速くなり、冷えると反応の一部が膨張に追いつけず、平衡から外れる。超新星の定量的な計算には重力、流体、輸送と反応網を同時に解く枠組みが必要である。

冷えてニュートリノが逃げられる中性子星物質を考える。電子だけで電荷を補う単純化の下では、荷電中性から $\rho_p=n_e$、$\beta$ 平衡から $\mu_n=\mu_p+\mu_e$ となる。両方向の反応で物質の自由エネルギーの一次変化がゼロになる条件である。ニュートリノを捕捉する場合は $\mu_n+\mu_{\nu_e}=\mu_p+\mu_e$ となる。極端相対論的な縮退電子では $\mu_e\simeq\hbar c(3\pi^2 n_e)^{1/3}$ である。核子質量差を省き :eqref[eq-nuc-symmetry-energy] の二次近似を使う。数密度 $\rho=\rho_p+\rho_n$ を固定して陽子分率 $Y_p=\rho_p/\rho$ を増やすと、核子を中性子から陽子に変えるので、核物質の一核子当たりのエネルギーの微分は $\partial(E/A)/\partial Y_p=\mu_p-\mu_n$ である。二次項 $S(\rho)(1-2Y_p)^2$ を微分すると $-4S(\rho)(1-2Y_p)$ となる。これを $\beta$ 平衡と荷電中性の条件に入れると
$$
4S(\rho)(1-2Y_p)=\hbar c(3\pi^2\rho Y_p)^{1/3}
\label{eq-nuc-beta-equilibrated-matter}
$$
を満たす。左辺は $0\le Y_p\le1/2$ で減少し、右辺は増加する。$S>0$ なら両端で大小関係が逆転するので、この範囲に解が一つある。実際にはミューオンが現れると中性条件が変わり、高密度での $S$ の外挿も不確かである。圧力の状態方程式は星の質量と半径を決める入力となる。外層の核格子、自由中性子、核子超流動は、均質冷物質の模型とは別に扱う。

:::example[小さな陽子分率]
$\rho=0.16\,\mathrm{fm}^{-3}$、$S=32\,\mathrm{MeV}$ とし、左辺で一旦 $2Y_p$ を無視する。$Y_p\simeq(128/197.327)^3/(3\pi^2\cdot0.16)\simeq0.058$ である。この値では $2Y_p$ が約 0.12 あるため、元の式を解くと約 $0.044$ に下がる。
:::

:::example[捕捉の効果]
同じ $\rho,Y_p$ で $\mu_{\nu_e}>0$ とする。平衡条件では $\mu_n-\mu_p=\mu_e-\mu_{\nu_e}$ であり、自由に逃げる場合の $\mu_e$ より右辺が小さい。$\mu_{\nu_e}$ を固定した比較では、左辺を下げ右辺を上げる方向である $Y_p$ の増加が必要になる。
:::

### 縮退した物質での直接 Urca 過程

冷たい中性子星の $n,p,e$ の模型では、$\beta$ 平衡の付近で $n\to p+e^-+\bar\nu_e$ と $p+e^-\to n+\nu_e$ がニュートリノを運び出せる。Fermi 面近くの励起だけが寄与すると、ニュートリノの運動量は熱エネルギーの尺度で小さい。主要な運動量保存は $\vb p_{Fn}\simeq\vb p_{Fp}+\vb p_{Fe}$ である。荷電中性と :eqref[eq-nuc-fermi-gas] から $p_{Fp}=p_{Fe}$ なので、ベクトルの三角条件は $p_{Fn}\le2p_{Fp}$ になる。これを三乗し、$p_{Fn}^3/p_{Fp}^3=(1-Y_p)/Y_p$ を代入すると $Y_p\ge1/9$ が得られる。これは電子だけで電荷を補い、温度を低くした極限の条件である。満たさないときは、別の核子が運動量を受け持つような過程を考える。追加の粒子穴状態が位相空間を狭める。

:::example[平衡分率と冷却チャネル]
:eqref[eq-nuc-beta-equilibrated-matter] の飽和密度での例では $Y_p\simeq0.044<1/9$ である。この電子だけの模型では、低温の Fermi 面の運動量だけを使う直接 Urca は成立しない。高密度へ進んだときの分率は $S(\rho)$ の形に依存する。
:::

### 核統計平衡と凍結

強い相互作用と電磁反応が速く、$\beta$ 過程が遅い高温の系を考える。陽子と中性子の総数を固定し、その範囲で組成が緩和するなら、核種 $(A,Z)$ の化学ポテンシャルは $Z\mu_p+(A-Z)\mu_n$ である。これは全陽子数と全中性子数の制約で自由エネルギーを最小にする Lagrange 乗数の式でもある。古典的な核種の数密度に代入すると
$$
n_{A,Z}=\frac{G_{A,Z}}{\lambda_{A,Z}^3}
\bigg(\frac{n_p\lambda_p^3}{G_p}\bigg)^Z
\bigg(\frac{n_n\lambda_n^3}{G_n}\bigg)^{A-Z}e^{B(A,Z)/(k_{\rm B}T)}
$$
を得る。$\lambda_i$ は前節と同じ熱波長である。指数の結合は核種の静止質量を自由核子の質量和から差し引いて得た。$n_p,n_n$ は自由核子だけの密度であり、核に入った陽子や中性子も含めた全数の条件から決める。これが :term[核統計平衡]{en="nuclear statistical equilibrium" yomi="かくとうけいへいこう"} の理想気体近似である。高密度では核間相互作用や電子の遮蔽、不均一相が効くため、この形では不足する。

:::example[平衡を維持する条件]
体積が変化する時間尺度を $\tau_{\rm exp}$ とする。小さな組成のずれ $u$ が $\dot u=-\lambda u$ で緩和する模型なら、その時間でずれは $e^{-\lambda\tau_{\rm exp}}$ 倍になる。$\lambda\tau_{\rm exp}=10$ ならほぼ緩和するが、$0.1$ なら約 $0.90$ が残る。温度降下によって反応が後者の範囲に入ると、平衡式に現在温度を入れただけでは組成を求められない。
:::

### 崩壊する芯と冷たい高密度物質

ニュートリノが逃げられる低温・高密度の物質では $n\leftrightarrow p+e$ の化学平衡条件は $\mu_n=\mu_p+\mu_e$ である。ニュートリノが捕捉される場合は $\mu_{\nu_e}$ が残り、条件は $\mu_n+\mu_{\nu_e}=\mu_p+\mu_e$ になる。超新星の時間発展ではニュートリノの輸送と組成を合わせて求める。中性子星では量としては圧力とエネルギー密度の関係を使って重力平衡を解く。その一般相対論的な構造方程式への接続は [一般相対論](/physics/general-relativity/) を参照する。

ミューオンを無視し、電荷中性で低温の自由な Fermi 面を使うと、直接 Urca の運動量条件は $k_{Fn}\le k_{Fp}+k_{Fe}=2k_{Fp}$ である。$k_F$ は密度の三乗根に比例するので $(1-x_p)^{1/3}\le2x_p^{1/3}$、すなわち $x_p\ge1/9$ となる。電子以外の負電荷が中性条件を分担すると $k_{Fe}$ が小さくなり、このしきいは変わる。これは反応が運動学的に可能になる条件であり、冷却率の大きさは弱行列要素や対ギャップにも依存する。

:::example[Fermi 波数の三角条件]
電子と陽子の波数が同じ理想系で、右辺に対する中性子の波数比を求める。

計算には次の値を用いる。
$$
x=\frac{2}{25}
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=\frac{\left(\frac{1 - x}{x}\right)^{1 \cdot \frac{1}{3}}}{2}\\&=\frac{\left(\frac{1 - \left(\frac{2}{25}\right)}{\left(\frac{2}{25}\right)}\right)^{\frac{1}{3}}}{2}\\&\simeq 1.12859.\end{aligned}
$$

比が 1 以下で三角形を作れる。

陽子分率を増やす。次の値を変え、ほかは固定する。
$$
x=\frac{3}{20}
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=\frac{\left(\frac{1 - \left(\frac{3}{20}\right)}{\left(\frac{3}{20}\right)}\right)^{\frac{1}{3}}}{2}\\&\simeq 0.891414,\\{R}_{\mathrm{new}}-{R}&\simeq -0.237176.\end{aligned}
$$
:::

:::example[電子圧力の近似]
超相対論の縮退電子で運動エネルギー密度が $\epsilon_e=3 n_e\mu_e/4$、圧力がその三分の一となる。

計算には次の値を用いる。
$$
n_{e}=\frac{1}{50}\,\mathrm{fm}^{-3},\qquad \mu_{e}=150\,\mathrm{MeV}
$$

この条件で求める量を ${P_e}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{P_e}&=\frac{\mu_{e} \cdot n_{e}}{4}\\&=\frac{\left(150\right) \cdot \left(\frac{1}{50}\right)}{4}\\&\simeq 0.75\,\mathrm{MeV\,fm}^{-3}.\end{aligned}
$$

ここでは比較のため $\mu_e$ を固定した。自由気体の自己無撞着な密度変化では $\mu_e$ も変わる。

化学ポテンシャルを固定したまま圧力の依存を見る。次の値を変え、ほかは固定する。
$$
n_{e}=\frac{1}{25}\,\mathrm{fm}^{-3}
$$

このときの値を ${P_e}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{P_e}_{\mathrm{new}}&=\frac{\left(150\right) \cdot \left(\frac{1}{25}\right)}{4}\\&\simeq 1.5\,\mathrm{MeV\,fm}^{-3},\\{P_e}_{\mathrm{new}}-{P_e}&\simeq 0.75\,\mathrm{MeV\,fm}^{-3}.\end{aligned}
$$
:::

:::example[ニュートリノを捕捉した化学差]
荷電反応の平衡から $\mu_n-\mu_p=\mu_e-\mu_\nu$ を求める。

計算には次の値を用いる。
$$
\mu_{e}=180\,\mathrm{MeV},\qquad \mu_{\nu}=0\,\mathrm{MeV}
$$

この条件で求める量を ${\Delta\mu}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\Delta\mu}&=\mu_{e} - \mu_{\nu}\\&=- \left(0\right) + \left(180\right)\\&=180\,\mathrm{MeV}.\end{aligned}
$$

ニュートリノの化学ポテンシャルが正なら同じ電子化学ポテンシャルで必要な核子の差は小さい。

輸送で逃げられない場合を比較する。次の値を変え、ほかは固定する。
$$
\mu_{\nu}=40\,\mathrm{MeV}
$$

このときの値を ${\Delta\mu}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\Delta\mu}_{\mathrm{new}}&=\left(180\right) - \left(40\right)\\&=140\,\mathrm{MeV},\\{\Delta\mu}_{\mathrm{new}}-{\Delta\mu}&=-40\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[膨張による凍結の指標]
熱核のネットワークの反応時間を $\tau_r$、膨張時間を $\tau_h$ として、一膨張時間で進める反応の尺度を求める。

計算には次の値を用いる。
$$
\tau_{h}=\frac{1}{10}\,\mathrm{s},\qquad \tau_{r}=\frac{1}{100}\,\mathrm{s}
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=\frac{\tau_{h}}{\tau_{r}}\\&=\frac{\left(\frac{1}{10}\right)}{\left(\frac{1}{100}\right)}\\&=10.\end{aligned}
$$

比が 1 より十分小さければ、その反応だけでは組成を平衡に追随させられない。

温度低下で反応が遅い場合を比較する。次の値を変え、ほかは固定する。
$$
\tau_{r}=\frac{1}{2}\,\mathrm{s}
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=\frac{\left(\frac{1}{10}\right)}{\left(\frac{1}{2}\right)}\\&\simeq 0.2,\\{R}_{\mathrm{new}}-{R}&\simeq -9.8.\end{aligned}
$$
:::

## 演習問題

::::exercise{level=basic}
$k_{\rm B}T$ を 8 倍にすると $E_G$ が一定の $E_0$ は何倍か。
:::solution
$E_0=(E_G(k_{\rm B}T)^2/4)^{1/3}$ なので、倍率は $8^{2/3}=4$ である。
:::
::::

::::exercise{level=basic}
$n_\alpha$ が 3 倍、温度が同じで中間体の消滅が崩壊に支配されるとき、炭素を作る体積反応率は何倍か。
:::solution
中間の ${}^8\mathrm{Be}$ は $n_\alpha^2$、さらに一つの捕獲で $n_\alpha$ が掛かる。従って $3^3=27$ 倍である。
:::
::::

::::exercise{level=basic}
中性子捕獲率が $\beta$ 崩壊率の 3 倍のとき、先に捕獲する確率はいくらか。
:::solution
競合する二過程の確率は率の比で、$\lambda_n/(\lambda_n+\lambda_\beta)=3/(3+1)=3/4$ である。
:::
::::

::::exercise{level=basic}
合成時の中性子対陽子比が $1/9$ のとき、全中性子が ${}^4\mathrm{He}$ に入る模型でヘリウムのバリオン分率を求めよ。
:::solution
$r=1/9$ とすると分率は $2r/(1+r)=2/10=0.2$ である。他の中性子を含む核を無視した模型での値である。
:::
::::

::::exercise{level=advanced}
非共鳴で $S(E)=S_0$ として :eqref[eq-nuc-thermal-rate] に Gaussian 近似を適用し、点火前の温度依存の主項を求めよ。
:::solution
積分の $E\sigma$ は $S_0 e^{-\sqrt{E_G/E}}$ である。$\Phi$ を $E_0$ で二次まで展開し、幅が $E_0$ に比べて小さいとき積分範囲を実数全体に延長する。積分は $S_0 e^{-3E_0/(k_{\rm B}T)}\sqrt{4\pi E_0 k_{\rm B}T/3}$ となる。前因子 $(k_{\rm B}T)^{-3/2}$ と合わせ、$E_0\propto T^{2/3}$ なので $\langle\sigma v\rangle\propto T^{-2/3}\exp[-3(E_G/(4k_{\rm B}T))^{1/3}]$ である。
:::
::::

::::exercise{level=advanced}
同種核子 $a$ が二つ反応して質量数が二倍の核 $b$ を作る。$\dot n_a$ と $\dot n_b$ を $\langle\sigma v\rangle$ で表し、バリオン数が保存することを示せ。
:::solution
率は $R=n_a^2\langle\sigma v\rangle/2$ である。固定体積で $\dot n_a=-2R$、$\dot n_b=R$ となる。$A_b=2A_a$ なので $A_a\dot n_a+A_b\dot n_b=-2A_aR+2A_aR=0$ である。
:::
::::

::::exercise{level=advanced}
中性子捕獲だけの三核種の直列模型で、最初の核種への一定供給を $s$ とする。各核種の捕獲率を $\lambda_1,\lambda_2,\lambda_3>0$ として、定常存在量を求めよ。
:::solution
定常では最初の式が $s-\lambda_1Y_1=0$、残りが $\lambda_{i-1}Y_{i-1}-\lambda_iY_i=0$ である。順に代入すると $\lambda_1Y_1=\lambda_2Y_2=\lambda_3Y_3=s$ なので $Y_i=s/\lambda_i$ となる。最後の捕獲による流出がなければ有限の定常解は成立しない。
:::
::::

::::exercise{level=advanced}
:eqref[eq-nuc-beta-equilibrated-matter] で $\rho$ を固定して $S$ を増やした場合、$Y_p$ がどう変化するか示せ。
:::solution
元の解で $0<Y_p<1/2$ だから $S$ を増やすと左辺が増える。固定した新しい $S$ で左辺は $Y_p$ の減少関数、右辺は増加関数である。差を再びゼロにするには $Y_p$ を増やす必要がある。これは電子だけの平衡と二次近似の下での帰結である。
:::
::::

## 確認テスト

::::quiz{answer="5" tol="1e-6"}
$E_G=500\,\mathrm{keV}$、$k_{\rm B}T=1\,\mathrm{keV}$ の Gamow ピークの位置を keV で求めよ。
:::explanation
指数の一次微分がゼロの式から $E_0^3=500(1)^2/4=125\,\mathrm{keV}^3$ なので $E_0=5\,\mathrm{keV}$ である。
:::
::::

::::quiz
s 過程の分岐しない定常流で各核種間で一定となる量を選べ。

- [ ] 存在量 $Y_A$ だけ。
- [ ] 捕獲率 $\lambda_A$ だけ。
- [x] 捕獲率と存在量の積 $\lambda_A Y_A$。

:::explanation
各核種で上流からの生成と自身の捕獲が釣り合うため、$\lambda_{A-1}Y_{A-1}=\lambda_A Y_A$ である。率が小さい核種ほど定常量は増える。
:::
::::

::::quiz{answer="0.25" tol="1e-6"}
中性子対陽子比が $1/7$ で全中性子がヘリウムに入る場合、ヘリウムに属するバリオンの分率を答えよ。
:::explanation
ヘリウム一つには中性子が二つ、バリオンが四つ入る。従って分率は $2n_n/(n_n+n_p)=2/(1+7)=0.25$ である。
:::
::::

::::quiz
冷たくニュートリノが逃げる中性子星物質での $\beta$ 平衡の式を選べ。

- [x] $\mu_n=\mu_p+\mu_e$
- [ ] $\mu_n+\mu_e=\mu_p$
- [ ] $\mu_n=\mu_p=0$

:::explanation
$n\leftrightarrow p+e^-+\bar\nu_e$ を進めても自由エネルギーの一次変化がゼロであり、逃げるニュートリノの化学ポテンシャルをゼロとする。よって核子間の差は電子の化学ポテンシャルに等しい。
:::
::::

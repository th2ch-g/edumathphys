---
title: β 崩壊と弱い相互作用
description: β 崩壊の質量収支とスペクトルから、弱い遷移演算子、パリティ非保存、ニュートリノを扱う。
notation:
  - symbol: 'E_0'
    meaning: 電子の全エネルギーで表した端点
  - symbol: 'F(Z,E)'
    meaning: 娘核の電場による Fermi 関数
  - symbol: 'ft'
    meaning: 位相空間因子と部分半減期の積
  - symbol: 'M_F,M_{\rm GT}'
    meaning: Fermi 型と Gamow–Teller 型の核行列要素
---

## この章の目標

- 原子質量から $\beta^\pm$ 崩壊と電子捕獲の $Q$ 値を求める。
- 三体の位相空間から許容 $\beta$ スペクトルを導く。
- $ft$ 値と角運動量を用いて弱い遷移を分類する。
- パリティ非保存とニュートリノ質量が観測量に現れる仕方を説明する。

## β 崩壊のエネルギー論と電子捕獲

核子レベルでは $\beta^-$ が $n\to p+e^-+\bar\nu_e$、$\beta^+$ が $p\to n+e^++\nu_e$、:term[電子捕獲]{en="electron capture" yomi="でんしほかく"} が $p+e^-\to n+\nu_e$ である。自由陽子では質量収支が違うため、真空での $\beta^+$ は起こらない。核の結合エネルギーを含めると一部の核種で開く。いずれも核の $A$ を保ち、:ref[def-nuc-nuclide] の核図表上で $Z$ が 1 変わる。

### 質量規約

核質量で書くと、放出型では $Q_{\beta^\pm}=(M_{\rm nuc,i}-M_{\rm nuc,f}-m_e)c^2$ であり、ニュートリノ質量はここでは省く。中性原子質量を用い、電子束縛の差やイオン化のエネルギーを無視する近似では
$$
Q_{\beta^-}=(M_{{\rm at},i}-M_{{\rm at},f})c^2,\qquad
Q_{\beta^+}=(M_{{\rm at},i}-M_{{\rm at},f}-2m_e)c^2,\qquad
Q_{\rm EC}=(M_{{\rm at},i}-M_{{\rm at},f})c^2 .
\label{eq-nuc-beta-q}
$$
$\beta^-$ では中性娘原子の電子数が親より一つ多く、その電子の静止質量が放出電子の質量と相殺する。実際の終状態が一価の娘陽イオンと自由電子なら、中性原子質量差から娘の第一イオン化エネルギーを引くと、その終状態の $Q$ になる。$\beta^+$ では中性娘原子の電子数が一つ少ない。中性娘原子と余った自由電子および陽電子を終状態に選ぶと、原子質量差から二電子質量を引く。娘の負イオンができる場合や電子の励起・放出を伴う場合は、その終状態のエネルギーで補正する。電子捕獲では初期原子の電子が一つ失われ、娘原子の電子数は中性の場合と一致する。捕獲直後の空孔を含む電子励起と、娘核の励起 $E^*$ は、中性原子の基底状態間の質量差から差し引く。

放出型は核の反跳を省いても電子とニュートリノがエネルギーを分けるので連続スペクトルである。捕獲では、特定の空孔と娘核準位を指定すればほぼ二体終状態であり、ニュートリノは線を作る。空孔の緩和による X 線や Auger 電子を測って捕獲を同定できる。

娘核の励起エネルギー $E^*$ が中性子分離エネルギー $S_n$ より大きければ、$\beta$ 崩壊で作った娘核から中性子が放出され得る。この二段階の過程を :term[β 遅発中性子放出]{en="beta-delayed neutron emission" yomi="べーたちはつちゅうせいしほうしゅつ"} という。親の $\beta$ 崩壊で到達できる領域と、中性子放出が開く領域の両方に入る娘準位は、反跳とニュートリノ質量を省けば $S_n<E^*<Q_{\beta^-}$ を満たす。実際の分岐はその領域へ入る $\beta$ 強度と、娘準位から中性子や $\gamma$ 線を出す部分幅で決まる。この過程が分裂生成物で現れる例は [核分裂と核融合](/physics/nuclear-physics/fission-fusion/) で扱う。

許容の電子捕獲では、核半径内で束縛電子の波動関数がほぼ一定として振幅を評価できる。非相対論近似で軌道角運動量が $l$ の動径波動関数は原点で $r^l$ に比例するので、先頭項には $s$ 軌道が寄与する。軌道 $a$ の占有数を $n_a$、電子波動関数を $\psi_a$ とすると、一つの娘核準位への部分率は $n_a|\psi_a(0)|^2 E_{\nu,a}^2$ に比例する。最初の二因子が捕獲される電子の密度であり、最後の因子は黄金律で $\int q^2\delta(E_{\nu,a}-cq)\dd q$ を積分して得られる。異なる殻の率を比べるときは、核行列要素が共通であることを使う。内殻の電子密度が大きくても、対応する空孔を作るエネルギーが足りなければ、その殻からの捕獲は閉じる。

:::example[殻による捕獲率の違い]
核電荷 $Z$ の非相対論的な水素様軌道なら、規格化した $ns$ 軌道の核位置での密度は $|\psi_{ns}(0)|^2=Z^3/(\pi a_0^3 n^3)$ である。$a_0$ は Bohr 半径である。$1s$ と $2s$ の占有数が等しく、両者のニュートリノ線エネルギーの差を省ける模型では、$1s$ の捕獲率は $2s$ の 8 倍になる。多電子原子では遮蔽と電子間の相関がこの密度比を変える。
:::

:::example[二つのしきい値]
原子質量差のエネルギーが $0.8\,\mathrm{MeV}$ なら、娘核基底状態への電子捕獲はエネルギー的に開き得る。$\beta^+$ では $0.8-1.022=-0.222\,\mathrm{MeV}$ なので閉じる。
:::

:::example[空孔の効果]
$Q_{\rm EC}=800\,\mathrm{keV}$、娘核の励起が $500\,\mathrm{keV}$、空孔励起が $10\,\mathrm{keV}$ なら、反跳を省いたニュートリノのエネルギーは $290\,\mathrm{keV}$ である。
:::

### 原子質量での三種のしきい値

裸の核質量を $M_N$ とし、電子の全結合エネルギーを捨てた近似で $M_{\rm at}=M_N+Zm_e$ とする。$\beta^-$ の核質量差からは電子の質量を引くが、親と娘の $Z$ の差がこれを打ち消す。$\beta^+$ では娘の $Z$ が一つ小さいので電荷差による補正と放出する陽電子の質量が同符号に足され、原子質量差から $2m_ec^2$ を引く。電子捕獲では初期に入る電子と原子質量の差の補正が打ち消し合う。励起した終核や空孔を指定するときは、その励起分をさらに引く。電子捕獲で光子線を測っても、空孔緩和の X 線をニュートリノの線と混同しない。

:::example[陽電子が開くか]
親と娘の原子質量差をエネルギーで $D$ とし、励起はゼロとする。

計算には次の値を用いる。
$$
D=\frac{4}{5}\,\mathrm{MeV},\qquad m_{ec2}=\frac{511}{1000}\,\mathrm{MeV}
$$

この条件で求める量を ${Q_+}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{Q_+}&=D - 2 \cdot m_{ec2}\\&=\left(\frac{4}{5}\right) - 2 \cdot \left(\frac{511}{1000}\right)\\&\simeq -0.222\,\mathrm{MeV}.\end{aligned}
$$

最初の設定では電子捕獲の窓が正でも陽電子は閉じる。

差を大きくした場合を比べる。次の値を変え、ほかは固定する。
$$
D=\frac{7}{5}\,\mathrm{MeV}
$$

このときの値を ${Q_+}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{Q_+}_{\mathrm{new}}&=- 2 \cdot \left(\frac{511}{1000}\right) + \left(\frac{7}{5}\right)\\&\simeq 0.378\,\mathrm{MeV},\\{Q_+}_{\mathrm{new}}-{Q_+}&\simeq 0.6\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[捕獲に残るエネルギー]
原子質量差 $Q$ から、終核の励起 $E_f$ と最終原子の空孔の励起 $B_h$ を引く。

計算には次の値を用いる。
$$
Q=500\,\mathrm{keV},\qquad E_{f}=200\,\mathrm{keV},\qquad B_{h}=30\,\mathrm{keV}
$$

この条件で求める量を ${E_\nu}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{E_\nu}&=- B_{h} - E_{f} + Q\\&=- \left(200\right) - \left(30\right) + \left(500\right)\\&=270\,\mathrm{keV}.\end{aligned}
$$

これは反跳を無視したニュートリノのエネルギーである。

外側軌道の空孔と比較する。次の値を変え、ほかは固定する。
$$
B_{h}=5\,\mathrm{keV}
$$

このときの値を ${E_\nu}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{E_\nu}_{\mathrm{new}}&=- \left(200\right) + \left(500\right) - \left(5\right)\\&=295\,\mathrm{keV},\\{E_\nu}_{\mathrm{new}}-{E_\nu}&=25\,\mathrm{keV}.\end{aligned}
$$
:::

:::example[捕獲の反跳]
質量のないニュートリノがほぼ全エネルギー $E$ を持つ。娘核の静止エネルギーを $m_c$ とすると、反跳は最低次で $E^2/(2m_c)$ である。

計算には次の値を用いる。
$$
E=100000\,\mathrm{eV},\qquad m_{c}=50000000000\,\mathrm{eV}
$$

この条件で求める量を ${E_R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{E_R}&=\frac{E^{2}}{2 \cdot m_{c}}\\&=\frac{\left(100000\right)^{2}}{2 \cdot \left(50000000000\right)}\\&\simeq 0.1\,\mathrm{eV}.\end{aligned}
$$

分光で小さな差を比較するときは、窓に比べて小さいことだけでは省略の根拠にならない。

高いニュートリノエネルギーと比べる。次の値を変え、ほかは固定する。
$$
E=500000\,\mathrm{eV}
$$

このときの値を ${E_R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{E_R}_{\mathrm{new}}&=\frac{\left(500000\right)^{2}}{2 \cdot \left(50000000000\right)}\\&\simeq 2.5\,\mathrm{eV},\\{E_R}_{\mathrm{new}}-{E_R}&\simeq 2.4\,\mathrm{eV}.\end{aligned}
$$

**一次近似の検証。** 反跳エネルギーはニュートリノの線エネルギーの二乗に比例する。線の位置を変えたとき、一次近似で反跳の変化を見積もる。

${E}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${E_R}$ の一階微分を求めると、
$$
\frac{\mathrm d{E_R}}{\mathrm dE}=\frac{E}{m_{c}}.
$$

この導関数を変更前の点で評価して $h=400000$ を掛ける。ここで $E$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq 0.8\,\mathrm{eV},\\\Delta-\Delta_{\mathrm{lin}}&\simeq 1.6\,\mathrm{eV}.\end{aligned}
$$

この変化では一次近似の誤差が大きいため、増分の計算に元の式を使う。
:::

:::example[β の励起分岐の端点]
基底間の $Q$ で終状態の励起を $E_f$ とする。ニュートリノ質量と核反跳を省く。

計算には次の値を用いる。
$$
Q=\frac{5}{2}\,\mathrm{MeV},\qquad E_{f}=\frac{3}{10}\,\mathrm{MeV}
$$

この条件で求める量を ${T_0}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{T_0}&=- E_{f} + Q\\&=- \left(\frac{3}{10}\right) + \left(\frac{5}{2}\right)\\&\simeq 2.2\,\mathrm{MeV}.\end{aligned}
$$

線ではなく連続スペクトルの上端である。

より高い終状態と比較する。次の値を変え、ほかは固定する。
$$
E_{f}=\frac{3}{2}\,\mathrm{MeV}
$$

このときの値を ${T_0}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{T_0}_{\mathrm{new}}&=- \left(\frac{3}{2}\right) + \left(\frac{5}{2}\right)\\&=1\,\mathrm{MeV},\\{T_0}_{\mathrm{new}}-{T_0}&\simeq -1.2\,\mathrm{MeV}.\end{aligned}
$$
:::

## Fermi の理論と Kurie プロット

移行エネルギーと移行運動量をエネルギー単位に直した値がともに $m_Wc^2$ より十分小さい場合、核の電流とレプトン電流を局所的に結合できる。これは [電弱相互作用](/physics/standard-model/electroweak/) の $W$ 交換振幅を低エネルギーで展開した近似である。さらに、電子とニュートリノの波長が核半径より十分長い許容近似では、レプトン運動量によらない核行列要素を先頭項に取る。以下では初期核を無偏極とし、レプトンの角度を積分したスペクトルを :ref[law-qm-golden] で求める。

電子の運動量を $p$、全エネルギーを $E$ とすると $p^2\dd p=(pE/c^2)\dd E$ である。質量ゼロのニュートリノでは $q^2\dd q=E_\nu^2\dd E_\nu/c^3$ であり、エネルギー保存で $E_\nu=E_0-E$ と置く。電子波が電場によって核付近で増減する効果を $F(Z_f,E)$ として掛けると
$$
\frac{\dd\lambda}{\dd E}=C F(Z_f,E)pE(E_0-E)^2,\qquad
K(E)=\sqrt{\frac{\dd\lambda/\dd E}{F(Z_f,E)pE}}=\sqrt C(E_0-E).
\label{eq-nuc-kurie}
$$
$C$ は核の行列要素と定数を含む。$K(E)$ を描く :term[Kurie プロット]{en="Kurie plot" yomi="きゅりーぷろっと"} の直線は、横軸の切片で端点を与える。$\beta^-$ の電子は娘核に引かれ、$\beta^+$ の陽電子は反発されるため、特に低い速度で $F$ が大きく変わる。点電荷の非相対論近似では $F=2\pi\eta/(1-e^{-2\pi\eta})$ と表せる。$\alpha$ は微細構造定数である。電子なら $\eta=Z_f\alpha c/v>0$、陽電子なら負にする。厳密な波の有限核サイズと遮蔽の補正には原子と核の波動関数が必要なので :cite[krane-nuclear] を参照する。

質量 $m_\nu$ をもつ一成分では、$q^2\dd q=E_\nu \sqrt{E_\nu^2-m_\nu^2c^4}\dd E_\nu/c^3$ になる。このため置き換えは $(E_0-E)^2\mapsto(E_0-E)\sqrt{(E_0-E)^2-m_\nu^2c^4}$ であり、端点は $E_0-m_\nu c^2$ に移る。反跳や分子状態を含む端点精密測定では、これらの寄与を分離する。

:::example[直線の比]
$E_0=2\,\mathrm{MeV}$ とする。$E=1.0$ と $1.5\,\mathrm{MeV}$ での $K$ の比は $1/(0.5)=2$ である。未補正の計数の比には $pE$ と電場の効果がさらに掛かる。
:::

:::example[相対論的な理想極限]
$F=1$、$m_e=0$、$c=1$ とするとスペクトルは $E^2(E_0-E)^2$ であり、積分は $E_0^5\int_0^1x^2(1-x)^2\dd x=E_0^5/30$ である。端点の違いだけでも率が大きく変わる理由が分かる。
:::

### 点電荷の Fermi 関数

非相対論の Coulomb 波で、単位入射振幅に規格化した原点の波動関数の絶対値二乗を考える。引力で $\eta>0$ とする規約では、波動関数の原点での規格化は $e^{\pi\eta/2}\Gamma(1-i\eta)$ である。この解の構成には Coulomb 散乱方程式の合流型超幾何関数解を使うので [散乱理論](/physics/advanced-quantum-mechanics/scattering/) の範囲に属する。原点での因子から先は Gamma 関数で計算できる。$\Gamma(1+i\eta)=i\eta\Gamma(i\eta)$ と反射公式 $\Gamma(i\eta)\Gamma(1-i\eta)=\pi/\sin(i\pi\eta)$ から
$$
e^{\pi\eta}|\Gamma(1+i\eta)|^2
=e^{\pi\eta}\frac{\pi\eta}{\sinh\pi\eta}
=\frac{2\pi\eta}{1-e^{-2\pi\eta}}
$$
を得る。$\eta\to0$ の極限は 1 である。引力の低速度極限では $F\sim2\pi\eta$ であり、スペクトルに含まれる $p$ の減少を補う。反発で $\eta\to-\infty$ なら $F\sim2\pi|\eta|e^{-2\pi|\eta|}$ でさらに抑えられる。

:::example[電子と陽電子での強度比]
速度と $|Z_f|$ が同じなら、上の近似で $F(\eta)/F(-\eta)=e^{2\pi\eta}$ である。この比は電場の因子だけを比較した値であり、異なる核種の二つの物理スペクトルの比には $Q$ と核行列要素も入る。
:::

### 位相空間の積分と Kurie 切片

電子の球対称な状態密度は $p^2\dd p$ であり、$E^2=p^2c^2+m_e^2c^4$ により $p^2\dd p=pE\dd E/c^2$ となる。質量ゼロのニュートリノには $p_\nu^2\dd p_\nu=E_\nu^2\dd E_\nu/c^3$ を使う。重い核の反跳を無視すると黄金律の $\delta(E+E_\nu-E_0)$ を後者に積分して $(E_0-E)^2$ が残る。核行列要素を一定、電子の Coulomb 歪みを $F(Z,E)$ とした場合、$\sqrt{(\dd N/\dd E)/(pEF)}$ が $E_0-E$ に比例する。非平坦な検出効率や複数終状態の重なりがある場合には、単一の切片に回帰する前に測定応答を扱う。

:::example[二点の Kurie 切片]
直線 $K=c(E_0-E)$ の二点 $(E_1,K_1),(E_2,K_2)$ から傾きを消去する。

計算には次の値を用いる。
$$
E_{1}=1\,\mathrm{MeV},\qquad E_{2}=2\,\mathrm{MeV},\qquad K_{1}=4,\qquad K_{2}=2
$$

この条件で求める量を ${E_0}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{E_0}&=\frac{- E_{1} \cdot K_{2} + E_{2} \cdot K_{1}}{K_{1} - K_{2}}\\&=\frac{- \left(1\right) \cdot \left(2\right) + \left(2\right) \cdot \left(4\right)}{- \left(2\right) + \left(4\right)}\\&=3\,\mathrm{MeV}.\end{aligned}
$$

$K$ の縦軸の単位を共通にすれば絶対規格化は消える。

高エネルギー点の補正が変わる場合を比べる。次の値を変え、ほかは固定する。
$$
K_{2}=1
$$

このときの値を ${E_0}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{E_0}_{\mathrm{new}}&=\frac{- \left(1\right)^{2} + \left(2\right) \cdot \left(4\right)}{- \left(1\right) + \left(4\right)}\\&\simeq 2.33333\,\mathrm{MeV},\\{E_0}_{\mathrm{new}}-{E_0}&\simeq -0.666667\,\mathrm{MeV}.\end{aligned}
$$

**一次近似の検証。** 縦軸の不確かさが切片にどう伝わるかを、第一点を固定して調べる。

${K_{2}}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${E_0}$ の一階微分を求めると、
$$
\frac{\mathrm d{E_0}}{\mathrm dK_{2}}=\frac{- E_{1} \cdot K_{2} - E_{1} \cdot \left(K_{1} - K_{2}\right) + E_{2} \cdot K_{1}}{\left(K_{1} - K_{2}\right)^{2}}.
$$

この導関数を変更前の点で評価して $h=-1$ を掛ける。ここで $K_{2}$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&=-1\,\mathrm{MeV},\\\Delta-\Delta_{\mathrm{lin}}&\simeq 0.333333\,\mathrm{MeV}.\end{aligned}
$$

この変化では一次近似の誤差が大きいため、増分の計算に元の式を使う。
:::

:::example[無質量極限のピーク]
電子質量と Coulomb 力を無視して $w(E)=E^2(E_0-E)^2$ とする。規格化積分は $E_0^5/30$、最大は $E=E_0/2$ である。

計算には次の値を用いる。
$$
E_{0}=3\,\mathrm{MeV}
$$

この条件で求める量を ${f_{\rm max}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{f_{\rm max}}&=\frac{30}{16 \cdot E_{0}}\\&=\frac{15}{8 \cdot \left(3\right)}\\&\simeq 0.625\,\mathrm{MeV}^{-1}.\end{aligned}
$$

確率密度の最大値は端点が広がるほど下がる。

二倍の端点と比べる。次の値を変え、ほかは固定する。
$$
E_{0}=6\,\mathrm{MeV}
$$

このときの値を ${f_{\rm max}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{f_{\rm max}}_{\mathrm{new}}&=\frac{15}{8 \cdot \left(6\right)}\\&\simeq 0.3125\,\mathrm{MeV}^{-1},\\{f_{\rm max}}_{\mathrm{new}}-{f_{\rm max}}&\simeq -0.3125\,\mathrm{MeV}^{-1}.\end{aligned}
$$
:::

:::example[位相空間の比の尺度]
無質量で形状因子の同じ二分岐を取り、端点 $E_1,E_2$ から積分の比を求める。

計算には次の値を用いる。
$$
E_{1}=3\,\mathrm{MeV},\qquad E_{2}=\frac{33}{10}\,\mathrm{MeV}
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=\left(\frac{E_{2}}{E_{1}}\right)^{5}\\&=\frac{\left(\frac{33}{10}\right)^{5}}{\left(3\right)^{5}}\\&\simeq 1.61051.\end{aligned}
$$

強度比に直すには核行列要素の比も必要である。

窓がより広い分岐と比べる。次の値を変え、ほかは固定する。
$$
E_{2}=4\,\mathrm{MeV}
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=\frac{\left(4\right)^{5}}{\left(3\right)^{5}}\\&\simeq 4.21399,\\{R}_{\mathrm{new}}-{R}&\simeq 2.60348.\end{aligned}
$$
:::

:::example[点核の Coulomb 補正]
非相対論の引力の Sommerfeld 因子を用いる。$u=2\pi Z\alpha/v_{\!c}$、$v_c=v/c$ とする。

計算には次の値を用いる。
$$
u=\frac{1}{5}
$$

この条件で求める量を ${F}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{F}&=\frac{u}{1 - e^{- u}}\\&=\frac{\left(\frac{1}{5}\right)}{1 - e^{- \left(\frac{1}{5}\right)}}\\&\simeq 1.10333.\end{aligned}
$$

この式は点核の $s$ 波の原点での増強を表し、有限サイズと相対論の全補正を省いている。

引力の効果が強い条件と比べる。次の値を変え、ほかは固定する。
$$
u=1
$$

このときの値を ${F}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{F}_{\mathrm{new}}&=\frac{\left(1\right)}{1 - e^{- \left(1\right)}}\\&\simeq 1.58198,\\{F}_{\mathrm{new}}-{F}&\simeq 0.478646.\end{aligned}
$$
:::

## ft 値と許容遷移・禁止遷移

$W=E/(m_ec^2)$、$p_W=\sqrt{W^2-1}$、$W_0=E_0/(m_ec^2)$ とする。許容スペクトルの無次元積分を
$$
f=\int_1^{W_0} F(Z_f,W)p_W W(W_0-W)^2\dd W
$$
と置く。全半減期を $t_{1/2}$、対象分岐の割合を $b$ とすると、部分半減期は $t=t_{1/2}/b$ である。:term[比較半減期]{en="comparative half-life" yomi="ひかくはんげんき"} $ft$ は、$\lambda_b=\ln2/t$ の位相空間の違いを取り除く。$\log ft$ は $ft$ の秒単位を 10 を底として対数にした値である。逆比例の関係は許容近似で同じ種類の遷移を比べるときに使える。電子捕獲の位相空間には束縛電子の波動関数とニュートリノの線エネルギーが入る。

レプトン波を核内で $e^{i\vb k\cdot\vb r}=1+i\vb k\cdot\vb r+\cdots$ と展開する。定数項で進むものを :term[許容遷移]{en="allowed transition" yomi="きょようせんい"} と呼ぶ。一次の項はパリティを変える座標ベクトルを含み、Fermi 型のスカラーまたは Gamow–Teller 型のスピンと結合して階数 $0,1,2$ の候補を作る。これが :term[第一禁止遷移]{en="first-forbidden transition" yomi="だいいちきんしせんい"} の基本である。$\Delta J=2$ の第一禁止では最高階数が必要となる。下位の階数には複数の核電流項があり、スペクトルの形状因子を単なる定数として扱えない。反跳次数の電流も必要な系統的展開は電流の多重極分解を要するため :cite[bohr-mottelson] を参照する。

レプトンの軌道階数 $L$ の先頭項は $(pR/\hbar)^L$ の次数であり、軌道パリティは $(-1)^L$ である。核子の主要なベクトル電荷は階数 0、軸性スピン電流は階数 1 なので、最高の合成階数は $L+1$ である。$|\Delta J|=L+1$ かつ $\pi_i\pi_f=(-1)^L$ の遷移は、この最高の階数を最初に必要とし、$L$ 次の一意禁止遷移と呼ばれる（$L\ge1$）。例えば $3^+\to0^+$ では許容項も第一禁止項も角運動量を担えず、同パリティの二次の階数 3 が先頭の候補になる。一方、同じ $\Delta J$ でもパリティの積が違えば、必要な軌道階数も変わる。

:::example[分岐の規格化]
全半減期が $10\,\mathrm s$、分岐が $0.1$ なら $t=100\,\mathrm s$ である。$f=50$ なら $ft=5000\,\mathrm s$、$\log ft\simeq3.699$ となる。
:::

:::example[許容でないパリティ]
$2^-\to0^+$ ではパリティが変わり $\Delta J=2$ である。定数項の遷移は消え、第一禁止の最高階数が候補になる。小さい率は禁制の目安だが、エネルギーが小さい場合にも位相空間によって率は小さくなる。
:::

### 小さな端点での位相空間

運動エネルギーを $T_e=E-m_ec^2$、利用可能な全運動エネルギーを $Q$ とする。$F\simeq1$、$Q\ll m_ec^2$、ニュートリノ質量なしとすると、$p=\sqrt{2m_eT_e}[1+O(T_e/(m_ec^2))]$、$E\simeq m_ec^2$ である。従って位相空間積分の主な $Q$ 依存は
$$
\int_0^Q T_e^{1/2}(Q-T_e)^2\dd T_e
=Q^{7/2}\left(\frac23-\frac45+\frac27\right)
=\frac{16}{105}Q^{7/2}
$$
となる。電子の引力で $F\propto1/p$ が支配的なら $p$ が相殺し、積分は $\int_0^Q(Q-T_e)^2\dd T_e=Q^3/3$ になる。したがって一律に $Q^5$ を使えない。比較する領域で電子が相対論的か、Coulomb の補正が何に比例するかを確認する。

:::example[低エネルギーでの率差]
上の $F=1$ の非相対論極限で $Q$ を 4 倍にすると、核行列要素を固定した率は $4^{7/2}=128$ 倍になる。相対論的な理想極限での $Q^5$ なら $4^5=1024$ である。どちらも同じ領域に同時に適用できるわけではない。
:::

### 部分半減期と禁止度の形状

同じ親から複数の娘準位に分岐がある。親の全半減期 $t_{1/2}$ と分岐 $b$ に対し、部分率 $b\ln2/t_{1/2}$ に対応する半減期は $t_{1/2}/b$ となる。位相空間積分 $f$ は終点だけでなく $F(Z,E)$ と形状因子の規約にも依存するため、比較には同一の規約が必要である。

核内でのレプトンの波を $e^{i\vb q\cdot\vb r}=1+i\vb q\cdot\vb r+\cdots$ と展開すると、定数項は軌道階数 0 で偶、一次は階数 1 で奇である。一次の項とスピン階数 1 を合成すると階数 $0,1,2$ が得られる。従って第一次禁止ではパリティが反転し、最大階数の $\Delta J=2$ では先頭の核演算子が限られる。その他の一次禁止では複数の振幅と干渉を含むので、一個の定数の $B$ だけから全形状を決められない。

:::example[log ft の計算]
無次元 $f$ と秒単位の全半減期 $t_0$、分岐 $b$ を使う。

計算には次の値を用いる。
$$
f=120,\qquad t_{0}=10\,\mathrm{s},\qquad b=\frac{2}{5}
$$

この条件で求める量を ${\log_{10}ft}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\log_{10}ft}&=\frac{\ln\left(\frac{f\cdot t_{0}}{b\cdot 1\,\mathrm s}\right)}{\ln 10}\\&=\frac{\ln{\left(\frac{\left(10\right) \cdot \left(120\right)}{\left(\frac{2}{5}\right)} \right)}}{\ln{\left(10 \right)}}\\&\simeq 3.47712.\end{aligned}
$$

分岐が弱い方の比較半減期は長い。

見落としやすい小さな分岐と比較する。次の値を変え、ほかは固定する。
$$
b=\frac{1}{10}
$$

このときの値を ${\log_{10}ft}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\log_{10}ft}_{\mathrm{new}}&=\frac{\ln{\left(\frac{\left(10\right) \cdot \left(120\right)}{\left(\frac{1}{10}\right)} \right)}}{\ln{\left(10 \right)}}\\&\simeq 4.07918,\\{\log_{10}ft}_{\mathrm{new}}-{\log_{10}ft}&\simeq 0.60206.\end{aligned}
$$
:::

:::example[非相対論の小端点の積分]
運動エネルギー $T$ で電子の状態密度を $T^{1/2}$ と置く。$\int_0^Q T^{1/2}(Q-T)^2\dd T=(16/105)Q^{7/2}$ となる。

計算には次の値を用いる。
$$
Q=\frac{1}{50}\,\mathrm{MeV}
$$

この条件で求める量を ${I}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{I}&=\frac{16 \cdot Q^{\frac{7}{2}}}{105}\\&=\frac{16 \cdot \left(\frac{1}{50}\right)^{\frac{7}{2}}}{105}\\&\simeq 1.72399\times10^{-7}\,\mathrm{MeV}^{7/2}.\end{aligned}
$$

Coulomb 補正を外した積分であり、実際の低速電子では引力の影響が大きい。

電子静止質量より十分小さい範囲で端点を二倍にする。次の値を変え、ほかは固定する。
$$
Q=\frac{1}{25}\,\mathrm{MeV}
$$

このときの値を ${I}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{I}_{\mathrm{new}}&=\frac{16 \cdot \left(\frac{1}{25}\right)^{\frac{7}{2}}}{105}\\&\simeq 1.95048\times10^{-6}\,\mathrm{MeV}^{7/2},\\{I}_{\mathrm{new}}-{I}&\simeq 1.77808\times10^{-6}\,\mathrm{MeV}^{7/2}.\end{aligned}
$$

**一次近似の検証。** 位相空間がエネルギーに敏感なので、核の窓の誤差が全率に拡大されうる。

${Q}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${I}$ の一階微分を求めると、
$$
\frac{\mathrm d{I}}{\mathrm dQ}=\frac{8 \cdot Q^{\frac{5}{2}}}{15}.
$$

この導関数を変更前の点で評価して $h=\frac{1}{50}$ を掛ける。ここで $Q$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq 6.03398\times10^{-7}\,\mathrm{MeV}^{7/2},\\\Delta-\Delta_{\mathrm{lin}}&\simeq 1.17468\times10^{-6}\,\mathrm{MeV}^{7/2}.\end{aligned}
$$

この変化では一次近似の誤差が大きいため、増分の計算に元の式を使う。
:::

:::example[長波長条件を比較]
代表的なレプトン運動量を $p_c$ のエネルギーに換算し、$R$ での $pR/\hbar$ を求める。

計算には次の値を用いる。
$$
p_{c}=2\,\mathrm{MeV},\qquad R=5\,\mathrm{fm},\qquad h_{c}=\frac{197327}{1000}\,\mathrm{MeV\,fm}
$$

この条件で求める量を ${x}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{x}&=\frac{R \cdot p_{c}}{h_{c}}\\&=\frac{\left(2\right) \cdot \left(5\right)}{\left(\frac{197327}{1000}\right)}\\&\simeq 0.0506773.\end{aligned}
$$

選択則で定数が消える場合、小さくても次の次数が主な振幅になる。

高い端点の例と比較する。次の値を変え、ほかは固定する。
$$
p_{c}=6\,\mathrm{MeV}
$$

このときの値を ${x}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{x}_{\mathrm{new}}&=\frac{\left(5\right) \cdot \left(6\right)}{\left(\frac{197327}{1000}\right)}\\&\simeq 0.152032,\\{x}_{\mathrm{new}}-{x}&\simeq 0.101355.\end{aligned}
$$
:::

:::example[一階追加の強度の尺度]
同じ無次元の核行列要素を仮定し、振幅に $x=pR/\hbar$ が一個多いとする。

計算には次の値を用いる。
$$
x=\frac{1}{20}
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=x^{2}\\&=\left(\frac{1}{20}\right)^{2}\\&\simeq 0.0025.\end{aligned}
$$

これは強度の目安であり、核の実際の行列要素が同じという仮定が必要である。

長波長パラメータを二倍にする。次の値を変え、ほかは固定する。
$$
x=\frac{1}{10}
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=\left(\frac{1}{10}\right)^{2}\\&\simeq 0.01,\\{R}_{\mathrm{new}}-{R}&\simeq 0.0075.\end{aligned}
$$
:::

### Coulomb 補正を外した許容積分の閉形式

$F=1$ とし、無次元の全電子エネルギーを $W$、端点を $W_0>1$ とする。単一の質量を含んだ $p_W=\sqrt{W^2-1}$ は残し、核反跳とニュートリノ質量は捨てる。$(W_0-W)^2$ を展開して $f_0=W_0^2 I_1-2W_0 I_2+I_3$ と分ける。ここで $I_j=\int_1^{W_0} W^j\sqrt{W^2-1}\dd W$ である。$p_0=\sqrt{W_0^2-1}$ と置くと、$u=W^2-1$ の置換で $I_1=p_0^3/3$、$I_3=p_0^5/5+p_0^3/3$ を得る。$I_2$ は $W=\cosh y$ とすると積分関数が $\cosh^2y\sinh^2y=(\cosh4y-1)/8$ になり、
$$
I_2=\frac{W_0 p_0(2W_0^2-1)-\ln(W_0+p_0)}{8}
$$
である。三項を合わせると
$$
f_0=\frac{p_0(2W_0^4-9W_0^2-8)+15W_0\ln(W_0+p_0)}{60}\label{eq-nuc-beta-free-phase-space}
$$
となる。$W_0\gg1$ では主項が $W_0^5/30$ である。端点が $W_0=1$ に近いと分子の二項が強く相殺する。数値積分なら $W=1+(W_0-1)x$ として正の積分関数を $0\le x\le1$ で積分すると安定である。この閉形式は $F=1$ で核の位相空間を見積もる基準になる。実測された $ft$ を計算する段階では、電子の電荷と終核の電荷に応じた補正を入れて積分する。

:::example[電子質量を残す積分の比較]
$W_0=2$ では $p_0=\sqrt3$ であり、:eqref[eq-nuc-beta-free-phase-space] は
$$
f_0=\frac{-12\sqrt3+30\ln(2+\sqrt3)}{60}\simeq0.3121
$$
である。無質量近似で $W_0^5/30$ を使うと $1.0667$ となり、ここでは大きく過大評価する。この差には積分領域の下端を 0 から 1 に変えた効果と、積分関数の $p_W$ が $W$ より小さい効果の両方が入る。一方、$W_0=10$ なら閉形式は約 $3173.5$ で、無質量近似の $3333.3$ に近づく。両方の比を比べれば、同じ核行列要素でも窓が小さい場合には電子の質量を保つ必要があることが分かる。
:::

## Fermi 型遷移と Gamow–Teller 型遷移

許容近似で核子を非相対論的に扱うと、ベクトル電流の時間成分はアイソスピンを変える密度、軸性電流の空間成分はスピン密度を与える。$\beta^-$ に対し
$$
O_F=\sum_i t_-(i),\qquad
\vb O_{\rm GT}=\sum_i\vb\sigma(i)t_-(i)
\label{eq-nuc-beta-operators}
$$
を用いる。$t_-$ は :ref[def-nuc-isospin] の符号規約で中性子を陽子へ変える。回転階数は順に 0 と 1、どちらも空間反転で偶である。従って Fermi 型は $\Delta J=0$、:term[Gamow–Teller 型遷移]{en="Gamow–Teller transition" yomi="がもふてらーがたせんい"} は $|J_i-J_f|\le1\le J_i+J_f$ である。後者では $0\to0$ が除かれる。角度を積分し初期射影を平均した許容率には $g_V^2|M_F|^2+g_A^2|M_{\rm GT}|^2$ が入る。ここで行列要素の二乗には各終射影の和と初期射影の平均を含める。

アイソスピンが厳密なら $O_F=T_-$ は異なる $T$ を結合しない。同じ $T$ のアナログ状態への強度は角運動量の代数により $T(T+1)-T_z(T_z-1)$ である。$\vb O_{\rm GT}$ は全スピンと全アイソスピンの積には一般に等しくない。完全な終状態を和した強度 $S_\mp$ の差は閉包で
$$
S_- - S_+ = \sum_{a=x,y,z}\langle[O_+^a,O_-^a]\rangle=3(N-Z)
$$
となる。異なる粒子の演算子は可換であり、同粒子では $(\sigma_i^a)^2=1$ と $[t_+,t_-]=2t_z$ を使った。有限の励起範囲や切り詰めた模型空間の強度だけでは、この閉包の条件を満たさない。

:::example[純 Fermi 型]
$0^+\to0^+$ の許容遷移には階数 0 だけが寄与する。$T=1,T_z=1\to0$ のアナログ遷移では Fermi 強度が $1(2)-1(0)=2$ である。
:::

このような $0^+\to0^+$ のアナログ遷移では、アイソスピン対称性の極限で Fermi 強度が代数だけで決まる。したがって位相空間を除いた量 $ft$ は、放射補正などを省けば、同じ $T=1$ のアナログ遷移に共通となる。この比較ができるのは $\lambda\propto f g_V^2 B_F$ かつ $t=\ln 2/\lambda$ だからである。実測値を精密に比較するには、荷電の異なる核の波動関数が Coulomb 相互作用などでずれる効果と、弱い遷移の放射補正をそれぞれ評価する。

:::example[純 Gamow–Teller 型の候補]
$1^+\to0^+$ では階数 1 が必要で、パリティは同じである。従って許容項では Gamow–Teller 型だけが候補となる。核構造に応じて振幅がさらに抑えられることはある。
:::

### 許容遷移でのスピン和

許容近似の遷移振幅を $g_V L_0 M_F+g_A\sum_i L_i M_i$ と書く。$M_F$ は階数 0、$M_i$ は階数 1 の三成分である。初期核が無偏極で終状態の磁気射影をすべて足すと、$\sum_{M_i,M_f}M_F^*M_i$ は回転に対して不変でなければならない空間ベクトルになる。全ての回転で不変なベクトルはゼロなので、この平均の交差項は消える。二つの階数 1 の成分の積を平均すると、回転不変な二階テンソルに比例し、対称部分は $\delta_{ij}$ 倍になる。従って角度を積分した許容の全率には、二つの強度の規格化された和が現れる。偏極した核では偏極軸を使ってベクトルが構成できるため、この消去を異方的な分布に適用してはならない。

:::example[混合遷移の全率]
同じ位相空間の二つの遷移で、$g_V^2|M_F|^2$ と $g_A^2|M_{\rm GT}|^2$ が、一方で $(2,0)$、もう一方で $(1,3)$ とする。無偏極の許容全率の比は $2:4=1:2$ である。核行列要素の振幅を足してから二乗する計算では、このスピン和の条件を再現できない。
:::

### Ikeda の和則

$O_{\pm,k}=\sum_i\sigma_{ik}t_{i\pm}$ とし、初期の磁気射影で平均した全強度を $S_\pm$ とする。全終状態の完全性で $S_- -S_+=\sum_k\langle[O_{+,k},O_{-,k}]\rangle$ となる。異なる核子の演算子は可換であり、同じ核子では $\sigma_{ik}^2=1$ と $[t_{i+},t_{i-}]=2t_{iz}$ を使う。三成分を合計して $6\sum_i\langle t_{iz}\rangle=3(N-Z)$ である。この :term[Ikeda の和則]{en="Ikeda sum rule" yomi="いけだのわそく"} は、完全な核子空間での一体のスピン・アイソスピン演算子の強度差を与える。限られたエネルギー窓への GT 強度と置き換えてはならない。荷電方向の符号は :ref[def-nuc-isospin] に従う。

:::example[和則の差]
中性子・陽子数を $N,Z$ とすると $S_- -S_+$ を求められる。

計算には次の値を用いる。
$$
N=24,\qquad Z=20
$$

この条件で求める量を ${D}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{D}&=3 \cdot \left(N - Z\right)\\&=- 3 \cdot \left(20\right) + 3 \cdot \left(24\right)\\&=12.\end{aligned}
$$

それぞれの総強度や個々の準位への強度は差だけからは得られない。

中性子過剰度を増す。次の値を変え、ほかは固定する。
$$
N=28
$$

このときの値を ${D}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{D}_{\mathrm{new}}&=- 3 \cdot \left(20\right) + 3 \cdot \left(28\right)\\&=24,\\{D}_{\mathrm{new}}-{D}&=12.\end{aligned}
$$
:::

:::example[Fermi の類似状態]
同じ $T$ の下降遷移で $T_z=m$ とする。

計算には次の値を用いる。
$$
T=1,\qquad m=1
$$

この条件で求める量を ${B_F}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{B_F}&=T \cdot \left(T + 1\right) - m \cdot \left(m - 1\right)\\&=- \left(1\right) \cdot \left(\left(1\right) - 1\right) + \left(1\right) \cdot \left(\left(1\right) + 1\right)\\&=2.\end{aligned}
$$

最も低い射影は下げ演算子で消える。

同じ多重項の下端からの下降を比べる。次の値を変え、ほかは固定する。
$$
m=-1
$$

このときの値を ${B_F}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{B_F}_{\mathrm{new}}&=- \left(-1\right) \cdot \left(\left(-1\right) - 1\right) + \left(1\right) \cdot \left(\left(1\right) + 1\right)\\&=0,\\{B_F}_{\mathrm{new}}-{B_F}&=-2.\end{aligned}
$$
:::

:::example[有効軸性結合の影響]
GT の振幅にだけ係数 $q$ を掛ける模型をとる。

計算には次の値を用いる。
$$
q=1,\qquad B=4
$$

この条件で求める量を ${B_{\rm eff}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{B_{\rm eff}}&=B \cdot q^{2}\\&=\left(1\right)^{2} \cdot \left(4\right)\\&=4.\end{aligned}
$$

同じ模型空間での比較である。演算子の縮約を採用すると裸の強度の和則と同じ数をそのまま照合できない。

有効演算子を使う場合と比べる。次の値を変え、ほかは固定する。
$$
q=\frac{3}{4}
$$

このときの値を ${B_{\rm eff}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{B_{\rm eff}}_{\mathrm{new}}&=\left(4\right) \cdot \left(\frac{3}{4}\right)^{2}\\&\simeq 2.25,\\{B_{\rm eff}}_{\mathrm{new}}-{B_{\rm eff}}&\simeq -1.75.\end{aligned}
$$

**一次近似の検証。** 振幅の縮約因子の不確かさを強度へ移すとき、二乗を一次で扱った場合の差も見積もる。

${q}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${B_{\rm eff}}$ の一階微分を求めると、
$$
\frac{\mathrm d{B_{\rm eff}}}{\mathrm dq}=2 \cdot B \cdot q.
$$

この導関数を変更前の点で評価して $h=- \frac{1}{4}$ を掛ける。ここで $q$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&=-2,\\\Delta-\Delta_{\mathrm{lin}}&\simeq 0.25.\end{aligned}
$$

残差が示すのは、指定した二点間で生じる線形化の誤差である。基準点からさらに離れた値について同じ精度を保証するものではない。
:::

:::example[混合遷移の軸成分の比]
角度積分した許容率が $B_F+g_A^2 B_G$ に比例するとし、$g_V=1$ とする。

計算には次の値を用いる。
$$
g_{A}=\frac{127}{100},\qquad B_{F}=1,\qquad B_{G}=1
$$

この条件で求める量を ${b_A}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{b_A}&=\frac{B_{G} \cdot g_{A}^{2}}{B_{F} + B_{G} \cdot g_{A}^{2}}\\&=\frac{\left(1\right) \cdot \left(\frac{127}{100}\right)^{2}}{\left(1\right) \cdot \left(\frac{127}{100}\right)^{2} + \left(1\right)}\\&\simeq 0.617283.\end{aligned}
$$

角相関を測るときは、率の各寄与をこの比だけで記述できない。

GT 成分の異なる遷移と比べる。次の値を変え、ほかは固定する。
$$
B_{G}=3
$$

このときの値を ${b_A}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{b_A}_{\mathrm{new}}&=\frac{\left(3\right) \cdot \left(\frac{127}{100}\right)^{2}}{\left(1\right) + \left(3\right) \cdot \left(\frac{127}{100}\right)^{2}}\\&\simeq 0.828729,\\{b_A}_{\mathrm{new}}-{b_A}&\simeq 0.211446.\end{aligned}
$$
:::

## パリティ非保存と Wu の実験

核を偏極させて $\langle\vb J\rangle$ を揃える。角運動量は軸ベクトル、電子の運動量 $\vb p$ は極ベクトルなので、鏡映に相当する空間反転で $\langle\vb J\rangle\cdot\vb p$ が符号を変える。電子の分布が $W(\theta)=W_0(1+a\cos\theta)$ と表される場合、$\theta$ は偏極軸と電子運動量のなす角である。パリティが保存されるなら $\theta$ と $\pi-\theta$ の分布が等しく、$a=0$ となる。

Wu らは低温で偏極した ${}^{60}\mathrm{Co}$ の $\beta$ 電子を測定し、核スピンと反平行な側への放出が多いことを示した。後続の $\gamma$ 線の角度分布を測ることで核の配向を独立に確かめ、試料が暖まって配向が失われると $\beta$ 線の非対称性も消えることを調べた。理論では弱い荷電電流の左カイラルな結合がパリティを破る。パリティは左右のカイラリティを入れ替えるため、左成分だけを選ぶ電流は空間反転で元の式に戻らない。質量を持つレプトンでは、スピンの運動方向への射影であるヘリシティとカイラリティを区別する。

:::example[非対称性の読み方]
角度を鋭く $0,\pi$ に限る理想計測で $(W(0)-W(\pi))/(W(0)+W(\pi))=a$ である。例えば同じ効率で計数が 80 と 120 なら $a=-0.2$ である。
:::

:::example[無偏極の対照]
偏極軸がランダムなら同じ検出方向に対する $\cos\theta$ の平均がゼロになる。弱い相互作用に非対称な項が存在しても、無偏極の全角度計数ではこの項を直接測れない。
:::

### 軸性ベクトルと角相関の符号

空間反転作用で位置と運動量の符号が変わり、角運動量 $\vb J=\vb r\times\vb p+\vb S$ は符号が変わらない。従って $\vb J\cdot\vb p_e$ は反転で負になる。放出分布に $1+a\cos\theta$ が含まれるとき、反転対称なら $a=0$ が必要である。$\theta$ は核の偏極方向に対する電子角とする。前後半球の収率は $\cos\theta=x$ で $\int_0^1(1+ax)\dd x=1+a/2$ と $\int_{-1}^0(1+ax)\dd x=1-a/2$ に比例するので、規格化非対称度は $a/2$ である。Wu の実験では核の配向を準備してこの種の相関を調べ、偏極の大小や反転を手掛かりに装置の前後差と区別した。

:::example[偏極の低下]
分布の一次係数を $a=A P(v/c)$ とし、半球非対称度を求める。

計算には次の値を用いる。
$$
A=-1,\qquad P=\frac{4}{5},\qquad v_{c}=\frac{3}{4}
$$

この条件で求める量を ${\mathcal A}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\mathcal A}&=\frac{A \cdot P \cdot v_{c}}{2}\\&=\frac{\left(-1\right) \cdot \left(\frac{3}{4}\right) \cdot \left(\frac{4}{5}\right)}{2}\\&\simeq -0.3.\end{aligned}
$$

偏極が失われると同じ弱い相互作用でも観測の非対称度は減る。

偏極を半分にする。次の値を変え、ほかは固定する。
$$
P=\frac{2}{5}
$$

このときの値を ${\mathcal A}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\mathcal A}_{\mathrm{new}}&=\frac{\left(-1\right) \cdot \left(\frac{2}{5}\right) \cdot \left(\frac{3}{4}\right)}{2}\\&\simeq -0.15,\\{\mathcal A}_{\mathrm{new}}-{\mathcal A}&\simeq 0.15.\end{aligned}
$$
:::

:::example[前後の半球の収率比]
一次係数 $a$ として前半球を後半球で割る。

計算には次の値を用いる。
$$
a=- \frac{2}{5}
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=\frac{\frac{a}{2} + 1}{1 - \frac{a}{2}}\\&=\frac{\frac{\left(- \frac{2}{5}\right)}{2} + 1}{1 - \frac{\left(- \frac{2}{5}\right)}{2}}\\&\simeq 0.666667.\end{aligned}
$$

受容角と効率が両側で同じという仮定を使っている。

偏極の符号を反転する。次の値を変え、ほかは固定する。
$$
a=\frac{2}{5}
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=\frac{\frac{\left(\frac{2}{5}\right)}{2} + 1}{1 - \frac{\left(\frac{2}{5}\right)}{2}}\\&\simeq 1.5,\\{R}_{\mathrm{new}}-{R}&\simeq 0.833333.\end{aligned}
$$

**一次近似の検証。** 二つの半球の差を比で表すと、係数との関係は非線形になる。

${a}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${R}$ の一階微分を求めると、
$$
\frac{\mathrm d{R}}{\mathrm da}=\frac{4}{\left(a - 2\right)^{2}}.
$$

この導関数を変更前の点で評価して $h=\frac{4}{5}$ を掛ける。ここで $a$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq 0.555556,\\\Delta-\Delta_{\mathrm{lin}}&\simeq 0.277778.\end{aligned}
$$

残差が示すのは、指定した二点間で生じる線形化の誤差である。基準点からさらに離れた値について同じ精度を保証するものではない。
:::

:::example[有限な受容角の平均]
偏極軸側の円錐で $x=\cos\theta$ の下限を $u$ とし、等方な基準に対する平均倍率を求める。

計算には次の値を用いる。
$$
a=- \frac{2}{5},\qquad u=0
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=\frac{a \cdot \left(u + 1\right)}{2} + 1\\&=\frac{\left(- \frac{2}{5}\right) \cdot \left(\left(0\right) + 1\right)}{2} + 1\\&\simeq 0.8.\end{aligned}
$$

小さい前方受容では平均の $\cos\theta$ が大きい。

半球から狭い円錐へ移す。次の値を変え、ほかは固定する。
$$
u=\frac{4}{5}
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=\frac{\left(- \frac{2}{5}\right) \cdot \left(\left(\frac{4}{5}\right) + 1\right)}{2} + 1\\&\simeq 0.64,\\{R}_{\mathrm{new}}-{R}&\simeq -0.16.\end{aligned}
$$
:::

:::example[無偏極との二測定]
前後の効率が同じとし、固定した電子エネルギーで $a=A P\,v_c$ とする。測った半球非対称度 $r$ から $A$ を推定する。

計算には次の値を用いる。
$$
r=- \frac{3}{20},\qquad P=\frac{3}{5},\qquad v_{c}=\frac{3}{4}
$$

この条件で求める量を ${A}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{A}&=\frac{2 \cdot r}{P \cdot v_{c}}\\&=\frac{2 \cdot \left(- \frac{3}{20}\right)}{\left(\frac{3}{4}\right) \cdot \left(\frac{3}{5}\right)}\\&\simeq -0.666667.\end{aligned}
$$

小さな偏極では同じ統計誤差がより大きな係数の誤差になる。

非対称度の推定値を変える。次の値を変え、ほかは固定する。
$$
r=- \frac{1}{10}
$$

このときの値を ${A}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{A}_{\mathrm{new}}&=\frac{2 \cdot \left(- \frac{1}{10}\right)}{\left(\frac{3}{4}\right) \cdot \left(\frac{3}{5}\right)}\\&\simeq -0.444444,\\{A}_{\mathrm{new}}-{A}&\simeq 0.222222.\end{aligned}
$$
:::

### 偏極の反転で固定の効率を分離する

前後の検出効率を $\epsilon_F,\epsilon_B$ とし、偏極反転では大きさだけを保って符号を変える。半球で積分した物理的な非対称度を $b$ とする。正の偏極での計数の期待値は $C_F^+=I_+\epsilon_F(1+b), C_B^+=I_+\epsilon_B(1-b)$ と書ける。負の偏極では $C_F^-=I_-\epsilon_F(1-b), C_B^-=I_-\epsilon_B(1+b)$ である。よって交差する積の比を作ると
$$
R=\sqrt{\frac{C_F^+C_B^-}{C_B^+C_F^-}}=\frac{1+b}{1-b},\qquad b=\frac{R-1}{R+1}.
$$
固定された前後の効率差だけでなく、二測定の供給 $I_+,I_-$ の違いもこの模型で相殺する。実際の計数値でこれを使う場合、計数の分布の平均が比の式と等しいとは限らない。低計数で分母がゼロになる観測や、偏極反転に連動して効率が変わる測定には、率と効率を同時に当てはめる扱いが必要である。この例で比を取っているのは半球の計数であり、微分分布の係数を得るには受容角と偏極率をさらに使う。

## ニュートリノの性質と二重 β 崩壊

### 質量と混合

荷電電流で作る電子型ニュートリノは、質量固有状態の重ね合わせである。二成分の真空模型で $\ket{\nu_e}=\cos\theta\ket{\nu_1}+\sin\theta\ket{\nu_2}$ とする。相対論的近似で各成分は距離 $L$ を進む間に相対位相 $\Delta m^2 c^3 L/(2\hbar E)$ を得る。これは共通運動量で $E_i\simeq pc+m_i^2c^3/(2p)$ と展開し、伝播時間を $L/c$ とした結果である。元の状態に射影して絶対値を二乗すると
$$
P_{ee}=1-\sin^2(2\theta)
\sin^2\left(\frac{\Delta m^2 c^3 L}{4\hbar E}\right)
$$
となる。振動は質量二乗差を測るが、端点の形は絶対質量に感度を持つ。物質中でのコヒーレント散乱や三成分の一般化は [ニュートリノ](/physics/standard-model/neutrinos/) を参照する。

### 二回の核電荷の変化

一回の $\beta^-$ が閉じていても、同じ $A$ で $Z+2$ の中性原子が軽ければ、電子束縛の補正を省いた質量収支で二重 $\beta^-$ が開き得る。二電子と二反ニュートリノを出す $2\nu\beta\beta$ は、通常の弱い相互作用の二次過程である。弱い Hamiltonian を $H_W$、自由な中間状態の Hamiltonian を $H_0$ と書く。遷移演算子の二次項 $H_W(E_i-H_0+i0)^{-1}H_W$ に完全系を挿入すると、最初に電子と反ニュートリノの対を一つ出す順序には

$$
\sum_m\frac{\langle f;\mathrm{leptons}|H_W|m;e,\bar\nu\rangle
\langle m;e,\bar\nu|H_W|i\rangle}
{E_i-E_m-E_e-E_{\bar\nu}+i0}
$$

が現れる。他の放出順序も加え、同種レプトンを反対称化する。分母には中間核のエネルギー $E_m$ に加えて、最初に放出したレプトンの全エネルギーも入る。一回の崩壊が閉じていれば、この中間状態を実際に作ることはできず、エネルギー保存は始状態と全終状態に課される。

:term[ニュートリノを伴わない二重ベータ崩壊]{en="neutrinoless double beta decay" yomi="にゅーとりのをともなわないにじゅうべーたほうかい"} では、終状態が娘核と二電子だけであれば全レプトン数は 2 増える。軽い Majorana ニュートリノの交換で進む模型では振幅が $\sum_i U_{ei}^2m_i$ に比例する。二つの左カイラル頂点の間で質量に比例する項が必要となり、位相により成分が打ち消し合い得る。この質量因子の由来は [ニュートリノを放出しない二重 β 崩壊](/physics/standard-model/neutrinos/) で扱う。核内仮想運動量や多体行列要素の計算は本章の許容近似を超える。

:::example[電子エネルギーの和]
反跳を省くと $2\nu$ 過程では二電子の運動エネルギーの和は $Q$ 未満に連続的に分布する。仮に $0\nu$ 過程があれば、同じチャネルでは和が $Q$ に集中する。単一電子のエネルギーは $0\nu$ でも配分によって変わる。
:::

:::example[位相の相殺]
二成分の仮想模型で $U_{e1}^2m_1=a$、$U_{e2}^2m_2=-a$ ならこの機構の有効質量はゼロである。振動確率に入る $|U_{ei}|^2$ にはこの符号が現れない。
:::

### 二重崩壊の和エネルギーの形

$2\nu$ 過程を、反跳なし、四つのレプトンの質量なし、核行列要素が一定という模型で考える。電子エネルギー和を $K$ とする。電子側の相空間の積分は $\int_0^K E^2(K-E)^2\dd E=K^5/30$ である。一方、二つのニュートリノの和は $Q-K$ なので、その積分は $(Q-K)^5/30$ になる。したがって和スペクトルは
$$
\frac{\dd\lambda_{2\nu}}{\dd K}\propto K^5(Q-K)^5,\qquad 0\le K\le Q
$$
である。実際の電子は質量と電場の影響を受けるので低いエネルギー側はこの形からずれる。また分母に現れる中間準位のエネルギーを一定とした仮定も変更し得る。一方、端点ではほとんど全エネルギーを電子が持ち、二ニュートリノの位相空間が強い抑制を生む。エネルギー分解能が有限なら、この連続成分も端点付近の測定に入る。

:::example[端点の窓幅]
端点近傍で $K^5$ を $Q^5$ に近似し、最後の幅 $\Delta\ll Q$ を積分すると $Q^5\Delta^6/6$ に比例する。端点の近くで窓を半分にすると、この理想模型の連続事象数は $1/64$ になる。線状の成分を測る場合は分解能に応じて、同じ窓に入る線の割合も考える。
:::

### 質量を含む端点と二重崩壊

ニュートリノのエネルギーを $\epsilon$ とする。$p_\nu^2\dd p_\nu=\epsilon\sqrt{\epsilon^2-m_\nu^2c^4}\,\dd\epsilon/c^3$ であり、物理領域は $\epsilon\ge m_\nu c^2$ である。各質量固有状態は検出しない終状態でも互いに直交するので、電子のスペクトルは $|U_{ei}|^2$ で重み付けした和になる。端点への寄与を分解できず、しかも $\epsilon$ が各質量より十分大きい領域では平方根を展開できる。最初の補正は $-\sum_i |U_{ei}|^2m_i^2c^4/2$ である。これが有効な運動学的質量 $m_\beta^2=\sum_i|U_{ei}|^2m_i^2$ の定義につながる。二重 $\beta$ の無ニュートリノ過程を軽い Majorana 質量で解釈する場合に現れるのは、位相を含む別の結合 $\abs{\sum_i U_{ei}^2m_i}$ である。

:::example[端点近くの質量の効果]
エネルギー不足 $\epsilon$ が質量エネルギー $m_c$ 以上の場合、質量ゼロに対するニュートリノ位相空間の比を求める。

計算には次の値を用いる。
$$
m_{c}=\frac{2}{5}\,\mathrm{eV},\qquad eps=2\,\mathrm{eV}
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=\sqrt{1 - \left(\frac{m_{c}}{eps}\right)^{2}}\\&=\sqrt{1 - \frac{\left(\frac{2}{5}\right)^{2}}{\left(2\right)^{2}}}\\&\simeq 0.979796.\end{aligned}
$$

端点に近いほど比の落ち方が大きい。

質量と同程度の不足分を比較する。次の値を変え、ほかは固定する。
$$
eps=\frac{1}{2}\,\mathrm{eV}
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=\sqrt{1 - \frac{\left(\frac{2}{5}\right)^{2}}{\left(\frac{1}{2}\right)^{2}}}\\&\simeq 0.6,\\{R}_{\mathrm{new}}-{R}&\simeq -0.379796.\end{aligned}
$$

**一次近似の検証。** 平方根のしきいに近づくと感度が増す。端点から遠い点で得た一次式を使った場合の誤差を確認する。

${eps}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${R}$ の一階微分を求めると、
$$
\frac{\mathrm d{R}}{\mathrm deps}=\frac{m_{c}^{2}}{eps^{3} \cdot \sqrt{1 - \frac{m_{c}^{2}}{eps^{2}}}}.
$$

この導関数を変更前の点で評価して $h=- \frac{3}{2}$ を掛ける。ここで $eps$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq -0.0306186,\\\Delta-\Delta_{\mathrm{lin}}&\simeq -0.349177.\end{aligned}
$$

この変化では一次近似の誤差が大きいため、増分の計算に元の式を使う。
:::

:::example[二つの質量の非干渉和]
質量エネルギーを $m_1,m_2$、第二の重みを $w$ とし、有効質量を求める。

計算には次の値を用いる。
$$
m_{1}=\frac{1}{100}\,\mathrm{eV},\qquad m_{2}=\frac{1}{20}\,\mathrm{eV},\qquad w=\frac{1}{20}
$$

この条件で求める量を ${m_\beta c^2}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{m_\beta c^2}&=\sqrt{m_{1}^{2} \cdot \left(1 - w\right) + m_{2}^{2} \cdot w}\\&=\sqrt{\left(\frac{1}{100}\right)^{2} \cdot \left(1 - \left(\frac{1}{20}\right)\right) + \left(\frac{1}{20}\right)^{3}}\\&\simeq 0.0148324\,\mathrm{eV}.\end{aligned}
$$

終状態が直交するので相対位相は入らない。

重い成分の混合を増す。次の値を変え、ほかは固定する。
$$
w=\frac{1}{10}
$$

このときの値を ${m_\beta c^2}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{m_\beta c^2}_{\mathrm{new}}&=\sqrt{\left(\frac{1}{100}\right)^{2} \cdot \left(1 - \left(\frac{1}{10}\right)\right) + \left(\frac{1}{10}\right) \cdot \left(\frac{1}{20}\right)^{2}}\\&\simeq 0.0184391\,\mathrm{eV},\\{m_\beta c^2}_{\mathrm{new}}-{m_\beta c^2}&\simeq 0.00360669\,\mathrm{eV}.\end{aligned}
$$
:::

:::example[コヒーレントな質量の打ち消し]
軽い粒子交換の二つの実振幅 $a,b$ を想定する。相対位相が $0$ または $\pi$ なら質量因子は $\abs{a+b}$ である。

計算には次の値を用いる。
$$
a=\frac{1}{40}\,\mathrm{eV},\qquad b=\frac{1}{50}\,\mathrm{eV}
$$

この条件で求める量を ${m_{\beta\beta}c^2}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{m_{\beta\beta}c^2}&=\left|{a + b}\right|\\&=\left|{\left(\frac{1}{40}\right) + \left(\frac{1}{50}\right)}\right|\\&\simeq 0.045\,\mathrm{eV}.\end{aligned}
$$

直接の運動学的質量と違って相対位相で小さくできる。

第二の位相を反転する。次の値を変え、ほかは固定する。
$$
b=- \frac{1}{50}\,\mathrm{eV}
$$

このときの値を ${m_{\beta\beta}c^2}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{m_{\beta\beta}c^2}_{\mathrm{new}}&=\left|{\left(- \frac{1}{50}\right) + \left(\frac{1}{40}\right)}\right|\\&\simeq 0.005\,\mathrm{eV},\\{m_{\beta\beta}c^2}_{\mathrm{new}}-{m_{\beta\beta}c^2}&\simeq -0.04\,\mathrm{eV}.\end{aligned}
$$
:::

:::example[二つのニュートリノの残余エネルギー]
電子の総エネルギーから残る窓を $q$ とする。二ニュートリノの等方かつ質量ゼロの状態数を $\int_0^q x^2(q-x)^2\dd x=q^5/30$ とする。

計算には次の値を用いる。
$$
q=\frac{3}{5}\,\mathrm{MeV}
$$

この条件で求める量を ${I}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{I}&=\frac{q^{5}}{30}\\&=\frac{\left(\frac{3}{5}\right)^{5}}{30}\\&\simeq 0.002592\,\mathrm{MeV}^5.\end{aligned}
$$

これは電子側の他の因子を固定した比較である。端点に寄ると連続成分が強く抑えられる。

残余の窓を半分にする。次の値を変え、ほかは固定する。
$$
q=\frac{3}{10}\,\mathrm{MeV}
$$

このときの値を ${I}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{I}_{\mathrm{new}}&=\frac{\left(\frac{3}{10}\right)^{5}}{30}\\&\simeq 8.1\times10^{-5}\,\mathrm{MeV}^5,\\{I}_{\mathrm{new}}-{I}&\simeq -0.002511\,\mathrm{MeV}^5.\end{aligned}
$$
:::

## 演習問題

::::exercise{level=basic}
原子質量差のエネルギーが $2\,\mathrm{MeV}$ の $\beta^+$ で、基底状態への $Q$ を求めよ。
:::solution
:eqref[eq-nuc-beta-q] より $2m_ec^2\simeq1.022\,\mathrm{MeV}$ を引く。従って $Q\simeq0.978\,\mathrm{MeV}$ である。
:::
::::

::::exercise{level=basic}
全半減期 $20\,\mathrm s$、分岐 $0.25$、$f=25$ の遷移の $ft$ を求めよ。
:::solution
部分半減期は $20/0.25=80\,\mathrm s$ である。$ft=25\cdot80=2000\,\mathrm s$ となる。
:::
::::

::::exercise{level=basic}
$T=2,T_z=1$ のアナログ $\beta^-$ 遷移の Fermi 強度を求めよ。
:::solution
$T_-$ は $T_z=0$ へ移す。強度は $T(T+1)-T_z(T_z-1)=6-0=6$ である。
:::
::::

::::exercise{level=basic}
$2^+\to2^-$ と $0^+\to1^+$ のうち許容 Gamow–Teller 型が候補になるものを選べ。
:::solution
許容演算子は空間反転で偶なので、前者の振幅は消える。後者は同パリティで $0\otimes1$ に $1$ を含むので候補になる。
:::
::::

::::exercise{level=advanced}
$m_e=0,F=1$ の理想模型の $\dd\lambda/\dd E$ が最大になる $E$ を求めよ。
:::solution
$E^2(E_0-E)^2$ は $[E(E_0-E)]^2$ である。$0\le E\le E_0$ では括弧が非負で、平方完成すると $E_0^2/4-(E-E_0/2)^2$ なので $E=E_0/2$ で最大になる。
:::
::::

::::exercise{level=advanced}
$W(\theta)\propto1+a\cos\theta$ を前方半球と後方半球で積分し、全計数に対する差の比を求めよ。
:::solution
方位角積分は共通である。$\mu=\cos\theta$ により前方は $\int_0^1(1+a\mu)\dd\mu=1+a/2$、後方は $\int_{-1}^0(1+a\mu)\dd\mu=1-a/2$ になる。差の比は $a/2$ である。
:::
::::

::::exercise{level=advanced}
質量 $m_\nu>0$ の一成分で電子スペクトルの端点直下のべき依存を求めよ。$E$ の電子側の因子は端点で有限かつ非ゼロとする。
:::solution
$\epsilon=E_0-m_\nu c^2-E$、$M=m_\nu c^2$ とするとニュートリノ因子は $(M+\epsilon)\sqrt{2M\epsilon+\epsilon^2}$ である。$\epsilon\downarrow0$ で先頭は $M\sqrt{2M}\epsilon^{1/2}$ となる。質量ゼロのときの二次の端点とは異なるべきである。
:::
::::

::::exercise{level=advanced}
$\theta=\pi/4$ の二成分振動模型で生存確率が最初にゼロになる $L>0$ を求めよ。
:::solution
混合の係数 $\sin^2(2\theta)=1$ である。正の最小距離では振動項の位相の絶対値が $\pi/2$ なので $L=2\pi\hbar E/(|\Delta m^2|c^3)$ となる。$\Delta m^2=0$ なら振動せず、この距離は存在しない。
:::
::::

## 確認テスト

::::quiz{answer="1.022" tol="0.002"}
娘の電子状態による補正を省く。$\beta^+$ が開くために基底状態間の中性原子質量差に必要な最低エネルギーを MeV で答えよ。
:::explanation
終状態の電子数の変更と陽電子放出のため $2m_ec^2$ が必要で、電子の静止エネルギー $0.511\,\mathrm{MeV}$ の 2 倍となる。等号では位相空間が消える。
:::
::::

::::quiz
許容の $0^+\to0^+$ に寄与し得る演算子を選べ。

- [x] Fermi 型の階数 0
- [ ] Gamow–Teller 型の階数 1
- [ ] パリティ奇の階数 2

:::explanation
角運動量ゼロ同士で結合できる階数は 0 だけである。Fermi 項は階数 0 で同パリティを結合する。
:::
::::

::::quiz{answer="3" tol="1e-6"}
許容の質量ゼロ近似で、$E_0-E$ が 3 倍になると Kurie 縦軸は何倍か。
:::explanation
位相空間を補正した :eqref[eq-nuc-kurie] の縦軸は $\sqrt C(E_0-E)$ であり、横軸端点からの距離に比例するので 3 倍である。
:::
::::

::::quiz
偏極軸と電子運動量の内積は空間反転でどうなるか。

- [ ] 変わらない。
- [x] 符号が逆になる。
- [ ] 常にゼロになる。

:::explanation
角運動量は $\vb r\times\vb p$ と同じく軸ベクトルで反転しない。運動量だけが反転するので内積の符号が逆になる。
:::
::::

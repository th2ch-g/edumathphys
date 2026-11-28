---
title: 中間子・ハドロンと原子核
description: pion 交換、クォーク模型、ハイパー核、QCD 物質を核の観測量に結び付ける。
notation:
  - symbol: 'B_{\rm bar},S,Q_{\rm el}'
    meaning: バリオン数、ストレンジネス、電荷
  - symbol: 'B_\Lambda'
    meaning: Λ の分離エネルギー
  - symbol: 'v_n'
    meaning: 重イオン衝突の方位角分布の Fourier 係数
---

## この章の目標

- pion の量子数と核力の長距離成分を結び付ける。
- クォークの電荷と風味からハドロンを分類する。
- ハイパー核の束縛と弱い崩壊を区別する。
- QCD の自由度と重イオン衝突で測る量の関係を説明する。

## π 中間子とその性質

核子と核子の間の長距離相互作用は最も軽い交換粒子に敏感である。:term[π 中間子]{en="pion" yomi="ぱいちゅうかんし"} はスピン 0、負パリティのハドロンで、$\pi^+,\pi^0,\pi^-$ は近似的なアイソスピン三重項である。質量エネルギーは中性が約 $135\,\mathrm{MeV}$、荷電が約 $140\,\mathrm{MeV}$ である。核子と同じ符号規約では $\pi^-$ がアイソスピン第三成分 $+1$ を持つ。粒子の内部のクォークと反クォークは後で区別する。

核内で仮想的に交換される pion と、外部へ実粒子として生成される pion は異なる。前者は [核力](/physics/nuclear-physics/nuclear-force/) で扱った静的な伝播因子 $1/(q^2+m_\pi^2 c^2/\hbar^2)$ を与える。後者を作るには終状態の静止質量と運動エネルギーが必要である。$\pi^\pm$ は主に弱い相互作用で $\mu^\pm$ とニュートリノへ崩壊し、$\pi^0$ は主に電磁相互作用で二光子へ崩壊する。

価クォーク模型では pion は相対軌道角運動量 $L=0$、スピン一重項の中間子である。フェルミ粒子と反粒子の固有パリティの積が $-1$ なので、軌道の因子を掛けると $P=(-1)^{L+1}=-1$ を得る。中性 pion では荷電共役の固有値も指定でき、$C=(-1)^{L+S}=+1$ である。これはクォークと反クォークの交換で空間部分に $(-1)^L$、スピン部分に $(-1)^{S+1}$、生成演算子の順序に $-1$ が入るためである。

:term[G パリティ]{en="G parity" yomi="じーぱりてぃ"} は荷電共役にアイソスピンの第二軸まわりの $\pi$ 回転を組み合わせた変換である。アイソスピンが良い近似では、三重項の中性成分 $T_z=0$ はこの回転で $(-1)^T=-1$ を受ける。したがって pion の $G$ は $-1$ であり、同じアイソスピン多重項の荷電成分にもこの値を割り当てる。電磁相互作用と質量差によるアイソスピンの破れを含む過程では、この強い相互作用の近似的な選択則に補正が入る。

:::example[静止した中性 pion の崩壊]
二光子の運動量は反対向きで大きさが等しく、各エネルギーは $m_{\pi^0}c^2/2\simeq67.5\,\mathrm{MeV}$ である。一光子だけでは運動量ゼロと正の全エネルギーを同時に満たせない。
:::

:::example[実粒子の生成しきい値]
$p+p\to p+p+\pi^0$ で実験室の一方の陽子が静止している場合、:eqref[eq-nuc-reaction-threshold] と同じ不変量を使う。終状態全体が重心系で静止する極限から $K_{\rm th}=2m_{\pi^0}c^2+m_{\pi^0}^2c^2/(2m_p)$ を得る。$m_pc^2=938\,\mathrm{MeV}$ とすると約 $280\,\mathrm{MeV}$ であり、生成する粒子の静止エネルギーだけでは足りない。
:::

### 荷電 pion のレプトン崩壊

質量の小さい電子の方が終状態の運動空間は広い。それでも $\pi^\pm$ ではミューオンの枝が強い。これを説明するため、この段では $\hbar=c=1$ とし、例えば $\pi^-\to\ell^-+\bar\nu_\ell$ を扱う。荷電弱電流のうち、スピン 0 でパリティが負の pion と真空を結ぶ軸性電流の行列要素は、回転と Lorentz 共変性から四元運動量 $P^\mu$ に比例する。その係数を pion の崩壊定数という。レプトンの行列要素と縮約すると $P_\mu\bar u_\ell\gamma^\mu(1-\gamma_5)v_\nu$ が現れる。

Feynman スラッシュを $\not{p}=\gamma^\mu p_\mu$ と定義する。$P=p_\ell+p_\nu$ と分けると $\bar u_\ell\not{p}_\ell=m_\ell\bar u_\ell$ である。$\not{p}_\nu(1-\gamma_5)=(1+\gamma_5)\not{p}_\nu$ と反交換して $\not{p}_\nu v_\nu=0$ を使うと、ニュートリノ側の項は消える。従って振幅は $m_\ell$ に比例する。スピン和で残る内積は $p_\ell\cdot p_\nu=(m_\pi^2-m_\ell^2)/2$ に比例し、二体の運動量の大きさも $(m_\pi^2-m_\ell^2)/(2m_\pi)$ である。位相空間まで合わせた木レベルの率は $m_\ell^2(1-m_\ell^2/m_\pi^2)^2$ に比例する。放射補正はこの近似では省いた。

:::example[ヘリシティで抑えられた枝]
$m_e=0.511$、$m_\mu=105.66$、$m_\pi=139.57$ を MeV 単位で使う。電子の枝とミューオンの枝の率の比は上の式から約 $1.28\times10^{-4}$ である。電子枝の運動量は大きくても、振幅の質量の係数が二乗で率に入る。
:::

### 同じアイソスピンをもつ散乱チャネル

pion の $T=1$ と核子の $t=1/2$ を合成すると、全アイソスピン $3/2$ と $1/2$ が生じる。例えば本章の第三成分の規約で $\pi^- n$ は $T_z=+3/2$ なので $T=3/2$ のみである。これに下げ演算子を一回作用させると、規格化した $T=3/2,T_z=1/2$ の状態は $(\sqrt2\ket{\pi^0 n}+\ket{\pi^-p})/\sqrt3$ となる。直交する $T=1/2,T_z=1/2$ は $(\ket{\pi^0 n}-\sqrt2\ket{\pi^-p})/\sqrt3$ と取れる。強い相互作用を荷電独立とすれば、この基底で散乱振幅は $f_{3/2},f_{1/2}$ に対角である。$\pi^-p$ の弾性振幅は $(f_{3/2}+2f_{1/2})/3$、電荷交換 $\pi^-p\to\pi^0 n$ は $\sqrt2(f_{3/2}-f_{1/2})/3$ となる。この計算は同じ運動学と部分波で比べたもので、質量差や Coulomb 相互作用を省いた。

:::example[共鳴の現れ方]
ある部分波で $f_{1/2}$ を無視できるとし、チャネルの速度も等しいとする。$\pi^-p$ の弾性強度と荷電交換強度は、それぞれ純粋な $T=3/2$ の強度の $1/9$ と $2/9$ である。$f_{1/2}$ が有限の場合は振幅を混ぜてから二乗するので、その強度を別々に足すことはできない。
:::

例えば $\Delta(1232)$ は $T=3/2$、$J^\pi=3/2^+$ の共鳴である。$\pi N$ の固有パリティの積は $-1$ なので、正パリティには奇数の相対軌道角運動量が要る。核子のスピン $1/2$ と合成して $J=3/2$ になる奇数の階数は $L=1$ である。したがって純粋な $T=3/2$ の $\pi^-n$ や $\pi^+p$ の弾性チャネルでは、対応する $p$ 波にこの共鳴が現れる。

### 擬スカラー中間子の二体崩壊

この節の崩壊の計算では $\hbar=c=1$、内積の計量は $(+,-,-,-)$ とする。中性 pion の二光子振幅の Lorentz 構造は $\epsilon^{\mu\nu\rho\sigma}\epsilon_{1\mu}\epsilon_{2\nu}k_{1\rho}k_{2\sigma}$ であり、各偏極を同じ光子の運動量で置き換えると反対称性でゼロになる。従って縦方向のゲージ成分が消える。係数が軸性カイラル異常から得られる理由は [量子異常](/physics/standard-model/anomalies/) に譲る。静止 pion では二光子が互いに逆向きで、それぞれ $m_\pi/2$ のエネルギーを持つ。運動する pion では Lorentz 変換により二線のエネルギーは一般に異なるが、二つを合わせた不変質量は変化しない。

荷電 pion では、軸性電流の真空との行列要素が pion の四元運動量 $P$ に比例する。Dirac 方程式で $P$ を二つのレプトンに分けると振幅が荷電レプトン質量に比例する。位相空間だけなら電子が軽い分だけ有利だが、この質量因子の二乗による抑制の方が大きい。

:::example[荷電 pion の終状態の運動量]
質量ゼロのニュートリノと質量 $m_l$ のレプトンへの二体崩壊で運動量を求める。

計算には次の値を用いる。
$$
m_{\pi}=\frac{13957}{100}\,\mathrm{MeV},\qquad m_{l}=\frac{511}{1000}\,\mathrm{MeV}
$$

この条件で求める量を ${p}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{p}&=\frac{- m_{l}^{2} + m_{\pi}^{2}}{2 \cdot m_{\pi}}\\&=\frac{\left(\frac{13957}{100}\right)^{2} - \left(\frac{511}{1000}\right)^{2}}{2 \cdot \left(\frac{13957}{100}\right)}\\&\simeq 69.7841\,\mathrm{MeV}.\end{aligned}
$$

電子の運動量が大きくてもレプトン質量に比例する振幅が小さいため分岐は少ない。

ミューオンの枝と比較する。次の値を変え、ほかは固定する。
$$
m_{l}=\frac{5283}{50}\,\mathrm{MeV}
$$

このときの値を ${p}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{p}_{\mathrm{new}}&=\frac{\left(\frac{13957}{100}\right)^{2} - \left(\frac{5283}{50}\right)^{2}}{2 \cdot \left(\frac{13957}{100}\right)}\\&\simeq 29.7906\,\mathrm{MeV},\\{p}_{\mathrm{new}}-{p}&\simeq -39.9935\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[質量抑制と位相空間の積]
共通定数を外した木レベルの荷電レプトンごとの幅を得る。

計算には次の値を用いる。
$$
m_{\pi}=\frac{13957}{100}\,\mathrm{MeV},\qquad m_{l}=\frac{511}{1000}\,\mathrm{MeV}
$$

この条件で求める量を ${W}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{W}&=m_{l}^{2} \cdot \left(1 - \frac{m_{l}^{2}}{m_{\pi}^{2}}\right)^{2}\\&=\left(\frac{511}{1000}\right)^{2} \cdot \left(1 - \frac{\left(\frac{511}{1000}\right)^{2}}{\left(\frac{13957}{100}\right)^{2}}\right)^{2}\\&\simeq 0.261114\,\mathrm{MeV}^2.\end{aligned}
$$

これは幅そのものではなく幅に比例する質量因子である。電磁放射は含まない。

ミューオンでの因子を比較する。次の値を変え、ほかは固定する。
$$
m_{l}=\frac{5283}{50}\,\mathrm{MeV}
$$

このときの値を ${W}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{W}_{\mathrm{new}}&=\left(\frac{5283}{50}\right)^{2} \cdot \left(1 - \frac{\left(\frac{5283}{50}\right)^{2}}{\left(\frac{13957}{100}\right)^{2}}\right)^{2}\\&\simeq 2034.49\,\mathrm{MeV}^2,\\{W}_{\mathrm{new}}-{W}&\simeq 2034.23\,\mathrm{MeV}^2.\end{aligned}
$$
:::

:::example[二光子の不変質量]
二光子のエネルギー $E_1,E_2$ と角 $\theta$ で、二つの四元運動量の和の二乗を取る。

計算には次の値を用いる。
$$
E_{1}=70\,\mathrm{MeV},\qquad E_{2}=65\,\mathrm{MeV},\qquad \theta=\pi\,\mathrm{rad}
$$

この条件で求める量を ${m_{\gamma\gamma}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{m_{\gamma\gamma}}&=\sqrt{2 \cdot E_{1} \cdot E_{2} \cdot \left(1 - \cos{\left(\theta \right)}\right)}\\&=\sqrt{2} \cdot \sqrt{\left(65\right) \cdot \left(70\right) \cdot \left(1 - \cos{\left(\left(\pi\right) \right)}\right)}\\&\simeq 134.907\,\mathrm{MeV}.\end{aligned}
$$

エネルギーの和だけでは親の不変質量にならない。

開き角を小さくする。次の値を変え、ほかは固定する。
$$
\theta=\frac{2 \cdot \pi}{3}\,\mathrm{rad}
$$

このときの値を ${m_{\gamma\gamma}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{m_{\gamma\gamma}}_{\mathrm{new}}&=\sqrt{2} \cdot \sqrt{\left(65\right) \cdot \left(70\right) \cdot \left(1 - \cos{\left(\left(\frac{2 \cdot \pi}{3}\right) \right)}\right)}\\&\simeq 116.833\,\mathrm{MeV},\\{m_{\gamma\gamma}}_{\mathrm{new}}-{m_{\gamma\gamma}}&\simeq -18.0742\,\mathrm{MeV}.\end{aligned}
$$

**一次近似の検証。** 開き角の測定誤差に対する不変質量の感度を調べる。背中合わせでは角度についての一次誤差が消える。

${\theta}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${m_{\gamma\gamma}}$ の一階微分を求めると、
$$
\frac{\mathrm d{m_{\gamma\gamma}}}{\mathrm d\theta}=- \frac{\sqrt{2} \cdot \sqrt{- E_{1} \cdot E_{2} \cdot \left(\cos{\left(\theta \right)} - 1\right)} \cdot \sin{\left(\theta \right)}}{2 \cdot \cos{\left(\theta \right)} - 2}.
$$

この導関数を変更前の点で評価して $h=- \frac{\pi}{3}$ を掛ける。ここで $\theta$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&=0\,\mathrm{MeV},\\\Delta-\Delta_{\mathrm{lin}}&\simeq -18.0742\,\mathrm{MeV}.\end{aligned}
$$

この変化では一次近似の誤差が大きいため、増分の計算に元の式を使う。
:::

:::example[中性 pion 交換の波長]
核力の単位に戻し、軸としての Compton 長を $\hbar/(mc)$ で求める。

計算には次の値を用いる。
$$
h_{c}=\frac{197327}{1000}\,\mathrm{MeV\,fm},\qquad m_{c}=\frac{6749}{50}\,\mathrm{MeV}
$$

この条件で求める量を ${L}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{L}&=\frac{h_{c}}{m_{c}}\\&=\frac{\left(\frac{197327}{1000}\right)}{\left(\frac{6749}{50}\right)}\\&\simeq 1.4619\,\mathrm{fm}.\end{aligned}
$$

質量だけの違いからは交換の結合係数の違いを決められない。

荷電成分の長さと比べる。次の値を変え、ほかは固定する。
$$
m_{c}=\frac{13957}{100}\,\mathrm{MeV}
$$

このときの値を ${L}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{L}_{\mathrm{new}}&=\frac{\left(\frac{197327}{1000}\right)}{\left(\frac{13957}{100}\right)}\\&\simeq 1.41382\,\mathrm{fm},\\{L}_{\mathrm{new}}-{L}&\simeq -0.048077\,\mathrm{fm}.\end{aligned}
$$
:::

:::example[ブースト後の光子端点]
親の速さを $c$ の $\beta$ 倍とし、静止系で進行方向に放出の光子のエネルギーを求める。

計算には次の値を用いる。
$$
m_{\pi}=\frac{6749}{50}\,\mathrm{MeV},\qquad \beta=0
$$

この条件で求める量を ${E_{\max}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{E_{\max}}&=\frac{m_{\pi} \cdot \left(\beta + 1\right)}{2 \cdot \sqrt{1 - \beta^{2}}}\\&=\frac{\left(\frac{6749}{50}\right) \cdot \left(\left(0\right) + 1\right)}{2 \cdot \sqrt{1 - \left(0\right)^{2}}}\\&\simeq 67.49\,\mathrm{MeV}.\end{aligned}
$$

逆方向の光子はエネルギーが低くなる。

親を運動させる。次の値を変え、ほかは固定する。
$$
\beta=\frac{3}{5}
$$

このときの値を ${E_{\max}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{E_{\max}}_{\mathrm{new}}&=\frac{\left(\frac{6749}{50}\right) \cdot \left(\left(\frac{3}{5}\right) + 1\right)}{2 \cdot \sqrt{1 - \left(\frac{3}{5}\right)^{2}}}\\&\simeq 134.98\,\mathrm{MeV},\\{E_{\max}}_{\mathrm{new}}-{E_{\max}}&\simeq 67.49\,\mathrm{MeV}.\end{aligned}
$$
:::

### 電荷とアイソスピンによる $\pi N$ 振幅の組合せ

ここでも :ref[def-nuc-isospin] の $T_z(p)=-1/2$、$T_z(\pi^+)=-1$ を使う。$\ket{\pi^+p}$ は最小射影 $T_z=-3/2$ であり、$T=3/2$ のみに属する。全上げ演算子を一回作用させると $\sqrt{2}\ket{\pi^0p}+\ket{\pi^+n}$ で、規格化は $\sqrt3$ である。符号を一つ選んだ直交状態は $\ket{T=1/2,T_z=-1/2}=(\ket{\pi^0p}-\sqrt{2}\ket{\pi^+n})/\sqrt3$ となる。したがって荷電独立な同じ部分波の散乱演算子を $a_{3/2},a_{1/2}$ で表すと
$$
\begin{aligned}
a(\pi^+p\to\pi^+p)&=a_{3/2},\\
a(\pi^+n\to\pi^+n)&=(a_{3/2}+2a_{1/2})/3,\\
a(\pi^+n\to\pi^0p)&=\sqrt2\,(a_{3/2}-a_{1/2})/3.
\end{aligned}
$$
これはチャネルの基底変換である。二つの固有振幅を先に二乗して足すと、同じ終状態への干渉を失う。例えば $a_{3/2}=a_{1/2}$ のとき電荷交換は消えるが、弾性は残る。現実の荷電チャネルを比べるときには Coulomb 力と pion の質量差も位相空間に効く。

## ハドロンのクォーク模型

クォークとグルーオンが作る色一重項の粒子や共鳴を :term[ハドロン]{en="hadron" yomi="はどろん"} と呼ぶ。価クォークによる最小構成では中間子が $q\bar q$、バリオンが $qqq$ である。この分類は価自由度の模型であり、波動関数には gluon や海のクォークも寄与する。より多くの価成分をもつハドロン候補もある。

$u,d,s$ クォークの電荷をそれぞれ $2e/3,-e/3,-e/3$、バリオン数を $1/3$ とする。ストレンジネスの規約では $s$ が $S=-1$、$u,d$ はゼロであり、反粒子ではすべての加法量の符号を反転する。従って $p=uud$、$n=udd$、$\Lambda=uds$ であり、$\pi^+=u\bar d$、$\pi^-=d\bar u$ である。荷電対称性の極限で $\pi^0$ の風味部分は $(u\bar u-d\bar d)/\sqrt2$ となる。

正味の価クォーク数を $n_u,n_d,n_s$ とすれば、この章の規約では $T_z=(n_d-n_u)/2$ である。:term[ハイパーチャージ]{en="hypercharge" yomi="はいぱーちゃーじ"} をここでは軽い三風味の $Y=B_{\rm bar}+S$ とする。$B_{\rm bar}=(n_u+n_d+n_s)/3$ と $S=-n_s$ を入れると $Y=(n_u+n_d-2n_s)/3$ となるので、各クォークの電荷の和は

$$
\frac{Q_{\rm el}}e=-T_z+\frac{Y}{2}
$$

である。Gell-Mann–Nishijima の関係を陽子の第三成分が正の規約で書く場合は、右辺の $-T_z$ をその規約の第三成分で置き換える。例えば $\Lambda$ と $\Sigma^0$ はどちらも $uds$ で $Y=T_z=0$ だが、$\Lambda$ は $T=0$、$\Sigma^0$ は $T=1$ である。$u,d$ の二体のアイソスピンに一重項と三重項の両方があり、$s$ のアイソスピンは 0 なので、同じ価クォーク数でこの違いが生じる。

### 色とスピン

クォークは三成分の色を持つ。$q\bar q$ の色一重項は $\sum_a\ket{a\bar a}/\sqrt3$ であり、色の回転 $U\otimes U^*$ で係数は $\sum_a U_{ba}U_{ca}^*=\delta_{bc}$ となって不変である。三クォークの色一重項は $\epsilon_{abc}\ket{abc}/\sqrt6$ である。$\SU(3)$ の回転でその係数は $\det U=1$ を掛けられるだけなので不変である。空間が対称な基底軌道にあれば、色部分の反対称性と合わせ、スピンと風味を組み合わせた部分に対称性が要る。スピン・風味の多重項全体の構成は [クォーク模型](/physics/standard-model/quark-model/) を参照する。

:::example[同じクォークでも異なる準位]
$\Delta^{++}=uuu$ では風味が対称であり、空間も対称ならスピンを完全対称な $3/2$ にする。すべて上向きのスピンから下げ演算子で得る四状態がこの対称多重項を作る。電荷は $+2e$、バリオン数は 1 である。
:::

:::example[反クォークのストレンジネス]
$K^+=u\bar s$ は $Q_{\rm el}/e=2/3+1/3=1$、$B_{\rm bar}=1/3-1/3=0$、$S=+1$ である。$\Lambda$ を伴う生成で両者の $S$ が打ち消し合うことがある。
:::

### 反クォークの量子数とスピン依存の質量

反クォークでは電荷、バリオン数、ストレンジネスがクォークと逆符号である。従って $u\bar s$ は電荷 $+1$、バリオン数 0、ストレンジネス $+1$ である。一方 $uds$ は電荷 0、バリオン数 1、ストレンジネス $-1$ であり、同じ組成でもスピンやアイソスピンの結合によって別の状態がある。価クォーク模型で中間子の色一重項は $\sum_i\ket{i\bar i}/\sqrt3$ である。$U\otimes U^*$ を作用させると $\sum_iU_{ai}U_{bi}^*=\delta_{ab}$ なので不変である。バリオンでは $\epsilon_{ijk}\ket{ijk}/\sqrt6$ が $U^{\otimes3}$ で $\det U=1$ 倍されるため不変である。$L=0$ の同種三クォークでは色以外の空間・スピン・味の積が対称になる。

簡単な超微細相互作用を $C\vb s_1\cdot\vb s_2$ とし、スピンは $\hbar$ で割ったものを使う。二スピンの合成が $S=0,1$ なら $\vb s_1\cdot\vb s_2=[S(S+1)-3/2]/2$ は $-3/4,1/4$ を取る。単純な接触模型では係数が逆質量積に比例するため、重いクォークを含むと分裂が小さくなる傾向を表せる。波動関数の原点密度を固定した場合の議論であり、軽い pion の質量をカイラル対称性なしで精密に再現するものではない。

:::example[中間子とバリオンの電荷]
価の $u,d,s$ の数から、反クォークを負符号で数えた数を $n_u,n_d,n_s$ とする。

計算には次の値を用いる。
$$
n_{u}=1,\qquad n_{d}=0,\qquad n_{s}=-1
$$

この条件で求める量を ${Q/e}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{Q/e}&=\frac{- n_{d} - n_{s} + 2 \cdot n_{u}}{3}\\&=- \frac{\left(-1\right)}{3} - \frac{\left(0\right)}{3} + \frac{2 \cdot \left(1\right)}{3}\\&=1.\end{aligned}
$$

初期の組は $u\bar s$、変更後は $uds$ である。

反クォークを含む中間子と三クォークを比べる。次の値を変え、ほかは固定する。
$$
n_{d}=1,\qquad n_{s}=1
$$

このときの値を ${Q/e}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{Q/e}_{\mathrm{new}}&=0\\&=0,\\{Q/e}_{\mathrm{new}}-{Q/e}&=-1.\end{aligned}
$$
:::

:::example[ハイパーチャージ]
三つの軽い味だけを扱うとき $Y=B+S$ である。$S=-n_s$、$B=(n_u+n_d+n_s)/3$ から求める。

計算には次の値を用いる。
$$
n_{u}=1,\qquad n_{d}=1,\qquad n_{s}=1
$$

この条件で求める量を ${Y}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{Y}&=- n_{s} + \frac{n_{d} + n_{s} + n_{u}}{3}\\&=0\\&=0.\end{aligned}
$$

同じバリオン数でもストレンジネスの値によって電荷多重項の関係が変わる。

ストレンジクォークを一つ増やす。次の値を変え、ほかは固定する。
$$
n_{d}=0,\qquad n_{s}=2
$$

このときの値を ${Y}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{Y}_{\mathrm{new}}&=\frac{\left(0\right)}{3} + \frac{\left(1\right)}{3} - \frac{2 \cdot \left(2\right)}{3}\\&=-1,\\{Y}_{\mathrm{new}}-{Y}&=-1.\end{aligned}
$$
:::

:::example[スピンの分裂を測る]
中心値を $M_0$、スピン係数を $C$ とする。三重項の値を計算する。

計算には次の値を用いる。
$$
M_{0}=900\,\mathrm{MeV},\qquad C=100\,\mathrm{MeV}
$$

この条件で求める量を ${M_1}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{M_1}&=\frac{C}{4} + M_{0}\\&=\frac{\left(100\right)}{4} + \left(900\right)\\&=925\,\mathrm{MeV}.\end{aligned}
$$

一重項は $M_0-3C/4$ で下がるため、二準位の差は $C$ になる。

スピン依存を強める。次の値を変え、ほかは固定する。
$$
C=200\,\mathrm{MeV}
$$

このときの値を ${M_1}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{M_1}_{\mathrm{new}}&=\frac{\left(200\right)}{4} + \left(900\right)\\&=950\,\mathrm{MeV},\\{M_1}_{\mathrm{new}}-{M_1}&=25\,\mathrm{MeV}.\end{aligned}
$$

**一次近似の検証。** 超微細結合の変化が、模型での三重項の位置にどのように入るかを調べる。

${C}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${M_1}$ の一階微分を求めると、
$$
\frac{\mathrm d{M_1}}{\mathrm dC}=\frac{1}{4}.
$$

この導関数を変更前の点で評価して $h=100$ を掛ける。ここで $C$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&=25\,\mathrm{MeV},\\\Delta-\Delta_{\mathrm{lin}}&=0\,\mathrm{MeV}.\end{aligned}
$$

入力のこの変化について差は一次の項だけで得られる。従って基準点での傾きで増分を求められる。
:::

:::example[三つの対の内積の和]
三つのスピン $1/2$ の合成スピンを $J$ とする。各対の和を全スピン二乗から抽出する。

計算には次の値を用いる。
$$
J=\frac{1}{2}
$$

この条件で求める量を ${S_\Sigma}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{S_\Sigma}&=\frac{J \cdot \left(J + 1\right) - \frac{9}{4}}{2}\\&=\frac{\left(\frac{1}{2}\right) \cdot \left(\left(\frac{1}{2}\right) + 1\right)}{2} - \frac{9}{8}\\&\simeq -0.75.\end{aligned}
$$

すべての対に同じスピン係数を仮定すると、この差が超微細分裂を決める。

全スピンがそろった状態と比べる。次の値を変え、ほかは固定する。
$$
J=\frac{3}{2}
$$

このときの値を ${S_\Sigma}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{S_\Sigma}_{\mathrm{new}}&=\frac{\left(\frac{3}{2}\right) \cdot \left(\left(\frac{3}{2}\right) + 1\right)}{2} - \frac{9}{8}\\&\simeq 0.75,\\{S_\Sigma}_{\mathrm{new}}-{S_\Sigma}&\simeq 1.5.\end{aligned}
$$
:::

:::example[スピン平均した中心値]
三重項と一重項の質量を $M_1,M_0$ とし、縮退度で重みをつける。

計算には次の値を用いる。
$$
M_{1}=940\,\mathrm{MeV},\qquad M_{0}=800\,\mathrm{MeV}
$$

この条件で求める量を ${\bar M}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\bar M}&=\frac{M_{0} + 3 \cdot M_{1}}{4}\\&=\frac{\left(800\right)}{4} + \frac{3 \cdot \left(940\right)}{4}\\&=905\,\mathrm{MeV}.\end{aligned}
$$

スピン係数の一次寄与が重み付きでは消える。

三重項の観測値を変える。次の値を変え、ほかは固定する。
$$
M_{1}=960\,\mathrm{MeV}
$$

このときの値を ${\bar M}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\bar M}_{\mathrm{new}}&=\frac{\left(800\right)}{4} + \frac{3 \cdot \left(960\right)}{4}\\&=920\,\mathrm{MeV},\\{\bar M}_{\mathrm{new}}-{\bar M}&=15\,\mathrm{MeV}.\end{aligned}
$$
:::

## ストレンジネスとハイパー核

:term[ハイパー核]{en="hypernucleus" yomi="はいぱーかく"} は核子とともにハイペロンを含む束縛系である。例えば $\Lambda$ は電荷ゼロ、ストレンジネス $-1$、スピン $1/2$ のバリオンである。核子と異なる粒子種なので、核子が既に占有した軌道と似た空間軌道に入ることを核子の Pauli 排他は妨げない。ハイペロンの分離エネルギーは、核子だけの芯の質量を $M_{\rm core}$ として
$$
B_\Lambda=(M_{\rm core}+m_\Lambda-M_{\rm hyper})c^2 .
\label{eq-nuc-hyperon-separation}
$$
通常の核子分離エネルギーと同じく、同じ質量規約を用いる。芯の励起を指定して測る場合は $M_{\rm core}$ をその状態に合わせる。

強い相互作用はクォーク風味の正味数を保存する。$S=0$ の入口から $\Lambda$ を作るとき、別の生成物が $S=+1$ を持てば強い相互作用で生成できる。弱い相互作用では正味の $S$ が変わる。核の中の $\Lambda$ は自由な中間子放出型に加え、近傍の核子と $\Lambda N\to NN$ の非中間子弱崩壊を起こし得る。核子終状態の Pauli ブロッキングが自由崩壊と異なる分岐比を生む。

標的の中性子を $\Lambda$ に変換する素過程として、$K^- n\to\pi^-\Lambda$ と $\pi^+n\to K^+\Lambda$ がある。前者では入口の $K^-$ が $S=-1$ を運び、後者では出口の $K^+$ が $\Lambda$ の $S=-1$ を打ち消す。他の粒子が核から出ないチャネルでは、残留核は中性子一つの代わりに $\Lambda$ 一つを含むため、バリオン数は標的と等しい。標的と出射中間子の四元運動量で欠損質量を求めると、生成したハイパー核の準位を調べられる。

:::example[随伴生成]
$p+p\to p+\Lambda+K^+$ では両辺の $B_{\rm bar}=2$、$Q_{\rm el}/e=2$、$S=0$ が一致する。ストレンジネスは生成物ごとにはゼロでなく、和が保存される。
:::

:::example[分離しきい値]
芯と自由な $\Lambda$ の合計静止エネルギーがハイパー核基底状態より $8\,\mathrm{MeV}$ 高いとする。基底状態から $6\,\mathrm{MeV}$ 励起しても自由な $\Lambda$ を出すチャネルは閉じる。弱い崩壊が開いているかは別の質量収支で決まる。
:::

### 分離エネルギーと芯の励起の測定

標的と入射ビームの総四元運動量を $P^\mu$、外へ出るメソンの四元運動量を $p^\mu$ とすると、残る系の質量は $(P-p)^2=M_{\rm miss}^2c^2$ で求まる。既知の芯の質量を使えば、ピークを :eqref[eq-nuc-hyperon-separation] の分離エネルギーに変換できる。ピークの位置にはハイペロンの軌道だけでなく、芯の粒子・穴励起が加わるため、複数の状態を同じ軌道に割り当てる際には量子数や遷移強度を合わせて調べる。分離エネルギーが大きな値は、同じ芯に対して残留系の質量が小さいことを意味する。

:::example[同じ芯への二つの状態]
同じ基底状態の芯を参照して、ハイパー核の二ピークが $B_\Lambda=12$ と $8\,\mathrm{MeV}$ に対応するとする。前者の質量が小さく、後者は前者より $4\,\mathrm{MeV}$ 高い。芯の基準が異なるピークを同様に引く場合は、芯の励起差を先に入れる。
:::

### ハイパー核のしきい値と分光

ハイパー核の $\Lambda$ は核子とは別の種類のフェルミ粒子であり、核子がすでに詰めた軌道に対応する空間形状を占められる。通常の $\Lambda$ は電荷ゼロなので Coulomb の障壁を持たないが、核子との残留相互作用はある。芯と $\Lambda$ の自由なしきい値を $M_{\rm core}+m_\Lambda$ と定め、束縛状態との差を $B_\Lambda$ とする。芯に $E_{\rm core}$ の励起があるときは、対応するしきい値がその分だけ上がる。二重項の間隔は芯の角運動量 $\vb I$ と $\Lambda$ の $\vb s$ が結合する模型では、$\vb I\cdot\vb s$ の二つの固有値の差になる。測った間隔から有効なスピン依存を制限できる。自由な弱崩壊と比べる場合、崩壊で生じた核子に対する Pauli 阻止を含めて評価する必要がある。

:::example[欠損質量の抽出]
四元運動量の差を $P_{\rm miss}$ とする。エネルギー成分を $E$、空間成分の大きさを $p$ として不変質量を求める。

計算には次の値を用いる。
$$
E=12000\,\mathrm{MeV},\qquad p=100\,\mathrm{MeV}
$$

この条件で求める量を ${M_{\rm miss}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{M_{\rm miss}}&=\sqrt{E^{2} - p^{2}}\\&=\sqrt{- \left(100\right)^{2} + \left(12000\right)^{2}}\\&\simeq 11999.6\,\mathrm{MeV}.\end{aligned}
$$

束縛を決めるのは静止質量に相当する値であり、反応の前後のエネルギー差だけではない。

反跳の大きい測定と比べる。次の値を変え、ほかは固定する。
$$
p=300\,\mathrm{MeV}
$$

このときの値を ${M_{\rm miss}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{M_{\rm miss}}_{\mathrm{new}}&=\sqrt{\left(12000\right)^{2} - \left(300\right)^{2}}\\&\simeq 11996.2\,\mathrm{MeV},\\{M_{\rm miss}}_{\mathrm{new}}-{M_{\rm miss}}&\simeq -3.33391\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[芯の励起を基準に戻す]
基底芯と自由 $\Lambda$ のしきい値より下に $B_0$ の束縛があり、同じ状態を芯の励起 $E_c$ をしきい値にして表す。

計算には次の値を用いる。
$$
B_{0}=7\,\mathrm{MeV},\qquad E_{c}=0\,\mathrm{MeV}
$$

この条件で求める量を ${B_\Lambda}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{B_\Lambda}&=B_{0} + E_{c}\\&=\left(0\right) + \left(7\right)\\&=7\,\mathrm{MeV}.\end{aligned}
$$

状態の質量を動かさず、参照するしきい値だけを変更している。

励起芯のしきい値を採る。次の値を変え、ほかは固定する。
$$
E_{c}=2\,\mathrm{MeV}
$$

このときの値を ${B_\Lambda}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{B_\Lambda}_{\mathrm{new}}&=\left(2\right) + \left(7\right)\\&=9\,\mathrm{MeV},\\{B_\Lambda}_{\mathrm{new}}-{B_\Lambda}&=2\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[芯とのスピン結合]
$\Lambda$ のスピンは $1/2$ とし、$a\vb I\cdot\vb s$ による二重項の間隔を求める。

計算には次の値を用いる。
$$
a=80\,\mathrm{keV},\qquad I_{c}=1
$$

この条件で求める量を ${D}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{D}&=\frac{a \cdot \left(2 \cdot I_{c} + 1\right)}{2}\\&=\frac{\left(80\right) \cdot \left(2 \cdot \left(1\right) + 1\right)}{2}\\&=120\,\mathrm{keV}.\end{aligned}
$$

$a>0$ では整列側が高い。芯がスピンゼロなら二重項は存在しない。

芯のスピンが大きい場合と比較する。次の値を変え、ほかは固定する。
$$
I_{c}=2
$$

このときの値を ${D}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{D}_{\mathrm{new}}&=\frac{\left(80\right) \cdot \left(2 \cdot \left(2\right) + 1\right)}{2}\\&=200\,\mathrm{keV},\\{D}_{\mathrm{new}}-{D}&=80\,\mathrm{keV}.\end{aligned}
$$
:::

:::example[媒質によるメソン放出の抑制]
二種の部分率を $\lambda_m,\lambda_n$ とし、前者だけを係数 $f$ で抑える。メソンの分岐比を求める。

計算には次の値を用いる。
$$
f=1,\qquad l_{m}=2\,\mathrm{ns}^{-1},\qquad l_{n}=1\,\mathrm{ns}^{-1}
$$

この条件で求める量を ${b}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{b}&=\frac{f \cdot l_{m}}{f \cdot l_{m} + l_{n}}\\&=\frac{\left(1\right) \cdot \left(2\right)}{\left(1\right) \cdot \left(2\right) + \left(1\right)}\\&\simeq 0.666667.\end{aligned}
$$

分岐比は抑制係数そのものでも全寿命でもない。

Pauli 阻止を模した率と比べる。次の値を変え、ほかは固定する。
$$
f=\frac{1}{4}
$$

このときの値を ${b}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{b}_{\mathrm{new}}&=\frac{\left(2\right) \cdot \left(\frac{1}{4}\right)}{\left(1\right) + \left(2\right) \cdot \left(\frac{1}{4}\right)}\\&\simeq 0.333333,\\{b}_{\mathrm{new}}-{b}&\simeq -0.333333.\end{aligned}
$$

**一次近似の検証。** 同じ部分率を維持した比較でも、分岐の相対変化は抑制係数の相対変化とは等しくない。

${f}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${b}$ の一階微分を求めると、
$$
\frac{\mathrm d{b}}{\mathrm df}=\frac{l_{m} \cdot l_{n}}{\left(f \cdot l_{m} + l_{n}\right)^{2}}.
$$

この導関数を変更前の点で評価して $h=- \frac{3}{4}$ を掛ける。ここで $f$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq -0.166667,\\\Delta-\Delta_{\mathrm{lin}}&\simeq -0.166667.\end{aligned}
$$

この変化では一次近似の誤差が大きいため、増分の計算に元の式を使う。
:::

:::example[相関した二 Λ の尺度]
単一の分離エネルギーを $B_1$、二つを同時に分けるエネルギーを $B_2$ とする。単純加算からの差を取る。

計算には次の値を用いる。
$$
B_{1}=3\,\mathrm{MeV},\qquad B_{2}=\frac{34}{5}\,\mathrm{MeV}
$$

この条件で求める量を ${\Delta B}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\Delta B}&=- 2 \cdot B_{1} + B_{2}\\&=- 2 \cdot \left(3\right) + \left(\frac{34}{5}\right)\\&\simeq 0.8\,\mathrm{MeV}.\end{aligned}
$$

芯の再配列もこの差に含まれうる。

ほぼ加算的な模型を比べる。次の値を変え、ほかは固定する。
$$
B_{2}=\frac{61}{10}\,\mathrm{MeV}
$$

このときの値を ${\Delta B}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\Delta B}_{\mathrm{new}}&=- 2 \cdot \left(3\right) + \left(\frac{61}{10}\right)\\&\simeq 0.1\,\mathrm{MeV},\\{\Delta B}_{\mathrm{new}}-{\Delta B}&\simeq -0.7\,\mathrm{MeV}.\end{aligned}
$$
:::

## 量子色力学と原子核

:term[量子色力学]{en="quantum chromodynamics" yomi="りょうししきりきがく"} は色の局所 $\SU(3)$ ゲージ理論である。ゲージ場の自己相互作用によって、高い運動量尺度で結合が弱まり、低い尺度では結合が強くなる。摂動の係数と閉じ込めの議論は [量子色力学](/physics/standard-model/qcd/) を参照する。摂動が使える高運動量では quark と gluon を自由度に選び、低運動量では色一重項の核子や pion を選ぶ。低エネルギー核力を直接低次の gluon 交換だけで求めることはできない。

核の質量エネルギーは nucleon の静止質量に大きく支配される。陽子質量には quark の質量項だけでなく、QCD の相互作用と場の運動のエネルギーが含まれる。核の結合の $B/A$ はその静止エネルギーより二桁小さい。QCD に基づく低エネルギー理論では近似的なカイラル対称性に基づき、軽い pion の交換と接触項を運動量の次数で整理する。短距離の係数を観測や格子 QCD と比較して定める。次数が上がると多 pion 交換や三体相互作用も現れる。格子での多核子状態の抽出には相関関数の信号と有限体積の扱いが必要であり、本章の価クォーク模型から核結合の数値を予測することはできない。

軽い二風味のクォーク質量をゼロとした極限では、左右の風味を独立に回す対称性がある。真空は両者を同時に回す対角部分だけを保つため、破れた三つの生成子に pion の三成分が対応する。小さいクォーク質量を戻すと対称性は明示的にも破れ、pion は小さな質量を持つ。対称性の破れと有効作用の構成は [カイラル対称性](/physics/standard-model/chiral-symmetry/) で扱う。核力の有効理論に pion を残す理由は、核子の典型的運動量が pion の質量尺度に近いとき、交換の伝播因子を小さな運動量の多項式として少数項で近似できないためである。逆に $\hbar q\ll m_\pi c$ なら伝播因子を展開でき、その先頭は運動量によらない接触項に吸収できる。

:::example[尺度の比較]
典型的な核子の静止エネルギーを $939\,\mathrm{MeV}$、結合を一核子当たり $8\,\mathrm{MeV}$ とすると比は約 $0.0085$ である。核の質量を比較して反応の $Q$ を測るときは、静止エネルギーに対する小さい差を取り出している。
:::

:::example[交換範囲の選別]
遠距離で二種の交換の動径因子が $e^{-m_1cr/\hbar}/r$ と $e^{-m_2cr/\hbar}/r$ で $m_2>m_1$ なら、重い方と軽い方の比は $e^{-(m_2-m_1)cr/\hbar}$ でゼロに近づく。短距離の形を変えても、軽い pion が長距離で重要である理由が分かる。
:::

### 尺度分離と格子の相関関数

摂動での一次の繰り込み群の流れを $\dd\alpha_s/\dd\ln Q=-\beta_0\alpha_s^2/(2\pi)$ とする。変数を分離して $1/\alpha_s(Q)-1/\alpha_s(Q_0)=\beta_0\ln(Q/Q_0)/(2\pi)$ となる。$\beta_0=11-2n_f/3>0$ なら短距離側で結合が減る。重い味のしきいを越えると有効な $n_f$ が変わるので、同じ係数をどの尺度にも延ばさない。核の低い運動量では、この摂動解を使って閉じ込めや結合を証明できない。

この段では $\hbar=c=1$ とする。有限の空間体積と無限の Euclid 時間を仮定し、$C(t)=\langle0|Oe^{-(H-E_{\rm vac})t}O^\dagger|0\rangle-|\langle0|O|0\rangle|^2$ を考える。完全系を入れると、真空以外の状態について $C(t)=\sum_n|\langle 0|O|n\rangle|^2e^{-(E_n-E_{\rm vac})t}$ となる。演算子と重なりを持つ最低状態のエネルギーを真空から測って $E_*$ とする。次の重なりを持つ状態とのギャップが $d>0$ なら、$C(t)=A e^{-E_*t}[1+O(e^{-dt})]$ と書ける。これを二時点の比に入れて対数を取れば、固定した $\Delta t>0$ での実効エネルギー $\ln[C(t)/C(t+\Delta t)]/\Delta t$ は $E_*+O(e^{-dt})$ である。最低状態との重なりが小さい場合は、励起の相対寄与を抑えるためにより大きい $t$ が要る。実際の格子計算では、そこで相関関数の統計誤差も確かめる。

:::example[一次の走る結合]
同じ味数の範囲で、基準の結合を $\alpha_0$ とする。尺度比 $r=Q/Q_0$ から計算する。

計算には次の値を用いる。
$$
a_{0}=\frac{1}{5},\qquad n_{f}=4,\qquad r=1
$$

この条件で求める量を ${\alpha_s}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\alpha_s}&=1 \cdot \frac{1}{\frac{\left(11 - \frac{2 \cdot n_{f}}{3}\right) \cdot \ln{\left(r \right)}}{2 \cdot \pi} + 1 \cdot \frac{1}{a_{0}}}\\&=\frac{1}{\frac{1}{2 \cdot \pi} \cdot \left(11 - \frac{2 \cdot \left(4\right)}{3}\right) \cdot \ln{\left(\left(1\right) \right)} + \frac{1}{\left(\frac{1}{5}\right)}}\\&\simeq 0.2.\end{aligned}
$$

どちらの尺度も摂動が有効であると仮定した結果である。

短距離側に移す。次の値を変え、ほかは固定する。
$$
r=2
$$

このときの値を ${\alpha_s}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\alpha_s}_{\mathrm{new}}&=\frac{1}{\frac{1}{2 \cdot \pi} \cdot \left(11 - \frac{2 \cdot \left(4\right)}{3}\right) \cdot \ln{\left(\left(2\right) \right)} + \frac{1}{\left(\frac{1}{5}\right)}}\\&\simeq 0.168938,\\{\alpha_s}_{\mathrm{new}}-{\alpha_s}&\simeq -0.0310615.\end{aligned}
$$

**一次近似の検証。** スケールの比を先に一次の式で扱う場合と、対数のまま解いた場合を比較する。

${r}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${\alpha_s}$ の一階微分を求めると、
$$
\frac{\mathrm d{\alpha_s}}{\mathrm dr}=\frac{6 \cdot \pi \cdot a_{0}^{2} \cdot \left(2 \cdot n_{f} - 33\right)}{r \cdot \left(a_{0} \cdot \left(2 \cdot n_{f} - 33\right) \cdot \ln{\left(r \right)} - 6 \cdot \pi\right)^{2}}.
$$

この導関数を変更前の点で評価して $h=1$ を掛ける。ここで $r$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq -0.0530516,\\\Delta-\Delta_{\mathrm{lin}}&\simeq 0.0219901.\end{aligned}
$$

この変化では一次近似の誤差が大きいため、増分の計算に元の式を使う。
:::

:::example[有効理論の展開比]
典型的運動量 $Q$ と破綻尺度 $L$ をとる。二次の相対補正の目安を計算する。

計算には次の値を用いる。
$$
Q=100\,\mathrm{MeV},\qquad L=600\,\mathrm{MeV}
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=\left(\frac{Q}{L}\right)^{2}\\&=\frac{\left(100\right)^{2}}{\left(600\right)^{2}}\\&\simeq 0.0277778.\end{aligned}
$$

無次元係数が自然な大きさという推定である。係数を較正した実際の不確かさとは区別する。

低運動量から離す。次の値を変え、ほかは固定する。
$$
Q=200\,\mathrm{MeV}
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=\frac{\left(200\right)^{2}}{\left(600\right)^{2}}\\&\simeq 0.111111,\\{R}_{\mathrm{new}}-{R}&\simeq 0.0833333.\end{aligned}
$$
:::

:::example[第一励起の相対寄与]
基底との重なりの二乗比を $r$、ギャップと時間の積を $u$ とする。

計算には次の値を用いる。
$$
r=\frac{1}{2},\qquad u=1
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=r \cdot e^{- u}\\&=\left(\frac{1}{2}\right) \cdot e^{- \left(1\right)}\\&\simeq 0.18394.\end{aligned}
$$

比が小さければ実効エネルギーの平坦化が期待できる。

長い Euclid 時間へ移す。次の値を変え、ほかは固定する。
$$
u=3
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=\left(\frac{1}{2}\right) \cdot e^{- \left(3\right)}\\&\simeq 0.0248935,\\{R}_{\mathrm{new}}-{R}&\simeq -0.159046.\end{aligned}
$$
:::

:::example[二時点からの実効エネルギー]
相関比 $r=C(t)/C(t+\Delta t)$ と、エネルギーの逆数で表した $\Delta t$ を使う。

計算には次の値を用いる。
$$
r=2,\qquad dt=\frac{1}{100}\,\mathrm{MeV}^{-1}
$$

この条件で求める量を ${E_{\rm eff}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{E_{\rm eff}}&=\frac{\ln{\left(r \right)}}{dt}\\&=\frac{\ln{\left(\left(2\right) \right)}}{\left(\frac{1}{100}\right)}\\&\simeq 69.3147\,\mathrm{MeV}.\end{aligned}
$$

励起混入があれば、この値は $t$ に依存する。

相関比の推定を変える。次の値を変え、ほかは固定する。
$$
r=\frac{11}{5}
$$

このときの値を ${E_{\rm eff}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{E_{\rm eff}}_{\mathrm{new}}&=\frac{\ln{\left(\left(\frac{11}{5}\right) \right)}}{\left(\frac{1}{100}\right)}\\&\simeq 78.8457\,\mathrm{MeV},\\{E_{\rm eff}}_{\mathrm{new}}-{E_{\rm eff}}&\simeq 9.53102\,\mathrm{MeV}.\end{aligned}
$$

**一次近似の検証。** 相関比の不確かさが実効エネルギーに伝わる様子を、この二点推定の式で確認する。

${r}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${E_{\rm eff}}$ の一階微分を求めると、
$$
\frac{\mathrm d{E_{\rm eff}}}{\mathrm dr}=\frac{1}{dt \cdot r}.
$$

この導関数を変更前の点で評価して $h=\frac{1}{5}$ を掛ける。ここで $r$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&=10\,\mathrm{MeV},\\\Delta-\Delta_{\mathrm{lin}}&\simeq -0.468982\,\mathrm{MeV}.\end{aligned}
$$

残差が示すのは、指定した二点間で生じる線形化の誤差である。基準点からさらに離れた値について同じ精度を保証するものではない。
:::

:::example[有限体積の指数を比べる]
最も軽い交換粒子の逆 Compton 長 $\lambda$ と周期箱の長さ $L$ から $e^{-\lambda L}$ を求める。

計算には次の値を用いる。
$$
lam=\frac{7}{10}\,\mathrm{fm}^{-1},\qquad L=4\,\mathrm{fm}
$$

この条件で求める量を ${f}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{f}&=e^{L \cdot \left(- lam\right)}\\&=e^{- \left(4\right) \cdot \left(\frac{7}{10}\right)}\\&\simeq 0.0608101.\end{aligned}
$$

これは一種類の長距離補正の指数部分である。二粒子の散乱状態では別の体積依存もある。

箱を大きくした場合を比べる。次の値を変え、ほかは固定する。
$$
L=6\,\mathrm{fm}
$$

このときの値を ${f}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{f}_{\mathrm{new}}&=e^{- \left(6\right) \cdot \left(\frac{7}{10}\right)}\\&\simeq 0.0149956,\\{f}_{\mathrm{new}}-{f}&\simeq -0.0458145.\end{aligned}
$$
:::

## クォーク・グルーオン・プラズマと重イオン衝突

高温の QCD 物質では、ハドロンの自由度だけで熱平衡を記述することが難しくなる。:term[クォーク・グルーオン・プラズマ]{en="quark–gluon plasma" yomi="くぉーくぐるーおんぷらずま"} は quark と gluon の自由度が熱力学や輸送に寄与する高温物質である。バリオン化学ポテンシャルが小さい領域では、実際のクォーク質量の QCD は温度を上げると滑らかに変化する。有限密度での相構造全体はこの結果から決まらない。

### 衝突から粒子分布へ

超相対論的重イオン衝突では、小さな体積に短時間でエネルギーが入り、生成した物質が膨張する。初期の空間的な形、圧力勾配による流れ、粒子への転換、後続の散乱を経て、検出器ではハドロンや光子を測る。衝突の中心性は重なりの大きさに関連する。中心性の実験的な区分には粒子多重度などの観測量を使う。単一事象の衝突径数がそのまま測れるわけではない。

運動量の方位角 $\phi$ と基準面 $\Psi_n$ に対し、分布の調和成分を
$$
v_n=\left\langle\cos[n(\phi-\Psi_n)]\right\rangle
\label{eq-nuc-flow}
$$
で定める。基準面が共通で、反射対称な分布なら $2\pi$ で規格化した密度は $[1+2\sum_n v_n\cos n(\phi-\Psi_n)]/(2\pi)$ と書ける。$v_2$ は楕円流を表す。圧力勾配による膨張が初期の形の異方性を運動量に写すため、初期状態や輸送係数の情報を含む。高運動量の粒子が少なくなる現象や重いクォークの運動も媒質との相互作用に感度を持つ。ただ一つの観測量から全過程が一意に決まるわけではない。

粒子収量の比を状態の熱力学に結び付けるときは、化学ポテンシャル、保存電荷、崩壊から流れ込む収量を区別する。高温相から検出されるハドロンへ移る理論には流体方程式とハドロン化の仮定が必要である。ここでは Fourier 係数を測る部分を扱い、その逆問題を一つの比例式で解くことはしない。

:::example[直交性での係数抽出]
$p(\phi)=[1+0.2\cos2\phi]/(2\pi)$ なら $v_2=0.1$ である。実際に $\cos2\phi$ を掛けると定数項の積分はゼロ、$\cos^2 2\phi$ の一周期平均は $1/2$ になる。
:::

:::example[熱自由度の尺度]
相互作用を無視した質量ゼロの gluon では、色 8 通り、横偏極 2 通りで 16 自由度である。相互作用のない二風味の quark と反 quark は色 3、スピン 2、風味 2、粒子と反粒子 2 で 24 自由度である。熱エネルギー密度の Bose 粒子単位では Fermi 積分が $7/8$ を与えるため、合計は $16+21=37$ である。$7/8$ は $(e^x+1)^{-1}=(e^x-1)^{-1}-2(e^{2x}-1)^{-1}$ を $x^3\dd x$ で積分して得る。この理想極限の値には相互作用と質量の補正が入る。
:::

### 高運動量の粒子と核の重なり

重イオンの高運動量成分の収量を、同じエネルギーの陽子衝突の収量と比較する。ある運動量区間で原子核衝突一事象当たりの収量を $Y_{AA}$、陽子衝突の収量を $Y_{pp}$ とする。独立な二体衝突の平均数 $\langle N_{\rm coll}\rangle$ で規格化した比 $R_{AA}=Y_{AA}/[\langle N_{\rm coll}\rangle Y_{pp}]$ は、各衝突が同じ確率で独立に粒子を作る基準では 1 である。$R_{AA}<1$ が観測される領域では、媒質中のエネルギー損失や入口の核構造の効果を比較する。$N_{\rm coll}$ は測定された粒子数そのものではなく、中心性区分に対する核の重なりの模型から得る。

:::example[独立衝突の基準との比]
ある区間で $Y_{AA}=4$、$Y_{pp}=0.1$、$\langle N_{\rm coll}\rangle=100$ とすれば $R_{AA}=0.4$ である。分母は 10 であり、陽子衝突の単純な収量より核衝突の収量が大きくても、独立衝突基準に対しては抑制されている。
:::

### 理想熱分布の自由度と流れの測定

ここでは $\hbar=c=k_{\rm B}=1$ とする。質量ゼロの気体では一自由度のエネルギー積分は $\int_0^\infty p^3\dd p/(e^{p/T}\mp1)$ に比例する。$p=Tx$ に変えると $T^4$ が外に出る。Bose の積分の分布を幾何級数に展開すると $6\sum_{n\ge1}n^{-4}$、Fermi は偶数項を引き直して $(1-2^{1-4})$ 倍である。グルーオンに色 8 と偏極 2、質量ゼロのクォークに味 $N_f$ と色 3、スピン 2、粒子・反粒子 2 を掛けると無次元の有効自由度は $16+(7/8)12N_f$ である。これは相互作用を外した高温極限の基準である。カイラル・クロスオーバーの近くでは有限質量と相互作用があり、格子の熱力学量からずれを調べる。

方位分布を $f(\phi)=[1+2v_2\cos 2(\phi-\Psi_2)+\cdots]/(2\pi)$ と書く。余弦の直交積分より $v_2=\langle\cos2(\phi-\Psi_2)\rangle$ である。単一の二次成分だけを残す模型は全方位で非負になるために $\abs{v_2}\le1/2$ を要する。実験の $\Psi_2$ には分解能があり、また非集団的な粒子対相関を落とすため異なる擬ラピディティ領域を比較する。流れだけでなく高い横運動量の粒子の相対減少や複数粒子の相関を合わせて、輸送や初期幾何への依存を制限する。

:::example[熱的自由度]
理想気体の有効自由度で、すべて質量ゼロと仮定した味数を変える。

計算には次の値を用いる。
$$
N_{f}=2
$$

この条件で求める量を ${g_{\rm eff}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{g_{\rm eff}}&=7 \cdot 12 \cdot N_{f} \cdot \frac{1}{8} + 16\\&=\frac{21 \cdot \left(2\right)}{2} + 16\\&=37.\end{aligned}
$$

実際の熱励起での重い味の寄与は連続的に変化するため、この例の整数切替は理想化である。

三味の理想極限と比べる。次の値を変え、ほかは固定する。
$$
N_{f}=3
$$

このときの値を ${g_{\rm eff}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{g_{\rm eff}}_{\mathrm{new}}&=\frac{21 \cdot \left(3\right)}{2} + 16\\&\simeq 47.5,\\{g_{\rm eff}}_{\mathrm{new}}-{g_{\rm eff}}&\simeq 10.5.\end{aligned}
$$
:::

:::example[理想気体の温度依存]
自由度を固定して温度比 $r$ だけでエネルギー密度の比を計算する。

計算には次の値を用いる。
$$
r=\frac{11}{10}
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=r^{4}\\&=\left(\frac{11}{10}\right)^{4}\\&\simeq 1.4641.\end{aligned}
$$

温度の小さな相対変化でも密度の変化は大きい。

温度増分を大きくする。次の値を変え、ほかは固定する。
$$
r=\frac{6}{5}
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=\left(\frac{6}{5}\right)^{4}\\&\simeq 2.0736,\\{R}_{\mathrm{new}}-{R}&\simeq 0.6095.\end{aligned}
$$
:::

:::example[面内と面外の方位比]
二次の方位係数だけを残し、基準角と直角の分布密度の比を取る。

計算には次の値を用いる。
$$
v_{2}=\frac{1}{20}
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=\frac{2 \cdot v_{2} + 1}{1 - 2 \cdot v_{2}}\\&=\frac{2 \cdot \left(\frac{1}{20}\right) + 1}{1 - 2 \cdot \left(\frac{1}{20}\right)}\\&\simeq 1.22222.\end{aligned}
$$

正の $v_2$ なら基準角で分布が多い。

変調の強い事象群と比べる。次の値を変え、ほかは固定する。
$$
v_{2}=\frac{3}{25}
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=\frac{2 \cdot \left(\frac{3}{25}\right) + 1}{1 - 2 \cdot \left(\frac{3}{25}\right)}\\&\simeq 1.63158,\\{R}_{\mathrm{new}}-{R}&\simeq 0.409357.\end{aligned}
$$

**一次近似の検証。** 流れが小さい領域からの外挿で、面内・面外の比がどれだけ変わるか調べる。

${v_{2}}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${R}$ の一階微分を求めると、
$$
\frac{\mathrm d{R}}{\mathrm dv_{2}}=\frac{4}{\left(2 \cdot v_{2} - 1\right)^{2}}.
$$

この導関数を変更前の点で評価して $h=\frac{7}{100}$ を掛ける。ここで $v_{2}$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq 0.345679,\\\Delta-\Delta_{\mathrm{lin}}&\simeq 0.0636777.\end{aligned}
$$

残差が示すのは、指定した二点間で生じる線形化の誤差である。基準点からさらに離れた値について同じ精度を保証するものではない。
:::

:::example[二粒子の係数]
独立放出で共通の事象面を仮定し、二つの粒子群の二次の係数を $v_a,v_b$ とする。方位差の二次モーメントは積になる。

計算には次の値を用いる。
$$
v_{a}=\frac{2}{25},\qquad v_{b}=\frac{1}{10}
$$

この条件で求める量を ${V_{2\Delta}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{V_{2\Delta}}&=v_{a} \cdot v_{b}\\&=\left(\frac{1}{10}\right) \cdot \left(\frac{2}{25}\right)\\&\simeq 0.008.\end{aligned}
$$

ジェットや崩壊相関がある場合には独立性が破れる。

第二群の流れを増す。次の値を変え、ほかは固定する。
$$
v_{b}=\frac{1}{5}
$$

このときの値を ${V_{2\Delta}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{V_{2\Delta}}_{\mathrm{new}}&=\left(\frac{1}{5}\right) \cdot \left(\frac{2}{25}\right)\\&\simeq 0.016,\\{V_{2\Delta}}_{\mathrm{new}}-{V_{2\Delta}}&\simeq 0.008.\end{aligned}
$$
:::

:::example[核修飾因子]
核衝突の収量 $Y$ を、二体衝突数 $N$ と一衝突当たりの参照収量 $y$ の積で割る。

計算には次の値を用いる。
$$
Y=\frac{4}{5},\qquad N=400,\qquad y=\frac{1}{500}
$$

この条件で求める量を ${R_{AA}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R_{AA}}&=\frac{Y}{N \cdot y}\\&=\frac{\left(\frac{4}{5}\right)}{\left(400\right) \cdot \left(\frac{1}{500}\right)}\\&=1.\end{aligned}
$$

比較は同じ運動量幅などの受容条件で行う。変化の解釈には初期状態の効果も検討する。

参照スペクトルの違いを比較する。次の値を変え、ほかは固定する。
$$
y=\frac{1}{250}
$$

このときの値を ${R_{AA}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R_{AA}}_{\mathrm{new}}&=\frac{\left(\frac{4}{5}\right)}{\left(400\right) \cdot \left(\frac{1}{250}\right)}\\&\simeq 0.5,\\{R_{AA}}_{\mathrm{new}}-{R_{AA}}&\simeq -0.5.\end{aligned}
$$
:::

### 縦方向の熱力学的膨張

ここでは $c=1$ とし、速度が縦位置と $v_z=z/t$ で結ばれる理想流体を考える。時空ラピディティ $\eta_s$ を $t=\tau\cosh\eta_s$、$z=\tau\sinh\eta_s$ によって定める。固有時間 $\tau=\sqrt{t^2-z^2}$ が一定の曲面上では、縦方向の微小固有長は $\sqrt{(\dd z)^2-(\dd t)^2}=\tau\dd\eta_s$ である。したがって横面積が固定の流体要素では $V\propto\tau$ となる。エネルギー密度 $\epsilon$ と圧力 $p$ を使い、仕事の収支 $\dd(\epsilon V)=-p\,\dd V$ に代入すると $\dd\epsilon/\dd\tau=-(\epsilon+p)/\tau$ になる。質量ゼロの理想気体では $p=\epsilon/3$ なので、$\epsilon\propto\tau^{-4/3}$、さらに $\epsilon\propto T^4$ から $T\propto\tau^{-1/3}$ である。横方向の流れ、粘性、クロスオーバー近くの状態方程式を加えると、温度の時間変化も変わる。

## 演習問題

::::exercise{level=basic}
$d\bar u$ の電荷とバリオン数を求めよ。
:::solution
電荷は $-e/3-2e/3=-e$ である。クォークと反クォークのバリオン数は $1/3,-1/3$ なので合計はゼロである。
:::
::::

::::exercise{level=basic}
$sud$ のバリオンについて $Q_{\rm el}/e$、$B_{\rm bar}$、$S$ を求めよ。
:::solution
電荷は $-1/3+2/3-1/3=0$、バリオン数は $3/3=1$ である。$s$ が一つあるので $S=-1$ である。これらの加法量だけでは異なるスピン・風味状態は区別できない。
:::
::::

::::exercise{level=basic}
$M_{\rm hyper}c^2=5000\,\mathrm{MeV}$、$(M_{\rm core}+m_\Lambda)c^2=5012\,\mathrm{MeV}$ とする。$B_\Lambda$ を求めよ。
:::solution
:eqref[eq-nuc-hyperon-separation] の通り分離後の静止エネルギーから初期の値を引き、$12\,\mathrm{MeV}$ を得る。
:::
::::

::::exercise{level=basic}
$p(\phi)=1/(2\pi)$ で $n>0$ の $v_n$ を求めよ。
:::solution
任意の固定された基準面について $\cos[n(\phi-\Psi_n)]$ は一周期で平均ゼロである。従ってすべての $v_n$ がゼロとなる。
:::
::::

::::exercise{level=advanced}
静止した $\pi^+$ の二体崩壊 $\pi^+\to\mu^++\nu_\mu$ でニュートリノ質量を無視する。ニュートリノエネルギーを両静止質量で表せ。
:::solution
静止系の運動量は大きさが共通で $p=E_\nu/c$ である。$m_\pi c^2-E_\nu$ がミューオン全エネルギーなので、これを二乗して $m_\mu^2c^4+E_\nu^2$ と等置する。$E_\nu^2$ を消して $E_\nu=(m_\pi^2-m_\mu^2)c^2/(2m_\pi)$ を得る。
:::
::::

::::exercise{level=advanced}
$p(\phi)=(1+2v_2\cos2\phi)/(2\pi)$ に他の成分がない場合、非負性から $v_2$ の範囲を求めよ。
:::solution
余弦が $-1$ と $+1$ の両方を取るので、分子の最小値は $1-2|v_2|$ である。従って必要十分条件は $|v_2|\le1/2$ である。他の成分と合わせた分布では同じ条件を単独で適用できない。
:::
::::

::::exercise{level=advanced}
三クォークの色状態 $\epsilon_{abc}\ket{abc}/\sqrt6$ の規格化と、一対の色を交換したときの固有値を示せ。
:::solution
$\epsilon_{abc}$ は色の全順列 6 個に対してだけ $\pm1$ である。色積基底は正規直交なのでノルム二乗は $6/6=1$ となる。二つの添字の交換は $\epsilon$ を反転させ、交換固有値は $-1$ である。
:::
::::

::::exercise{level=advanced}
自由な質量ゼロ gluon に三種類の味の自由な質量ゼロ quark と反 quark を加えた場合、熱エネルギー密度の自由度の重みを求めよ。
:::solution
gluons は 16、quark と反 quark は $3\cdot2\cdot3\cdot2=36$ である。後者の Fermi 積分を $7/8$ で重み付けし、$16+(7/8)36=47.5$ を得る。
:::
::::

## 確認テスト

::::quiz
強い相互作用による随伴生成で保存するものを選べ。

- [x] 全ストレンジネス
- [ ] 各粒子に割り当てたストレンジネス
- [ ] ハドロンの個数

:::explanation
強い相互作用は各 quark 風味の正味数を保つため、全ストレンジネスが保存する。相殺する quark と反 quark は対で生成でき、個々のハドロンや個数は変わる。
:::
::::

::::quiz{answer="67.5" tol="0.01"}
質量エネルギー $135\,\mathrm{MeV}$ の静止粒子が二光子だけへ崩壊した。片方の光子のエネルギーを MeV で答えよ。
:::explanation
反対向きで運動量の絶対値は等しく、質量ゼロの光子のエネルギーは等しい。従って $135/2=67.5\,\mathrm{MeV}$ である。
:::
::::

::::quiz
ハイパー核の $\Lambda$ と核子を区別する理由を選べ。

- [x] 粒子種が異なり、核子の同一粒子排他の対象にならない。
- [ ] $\Lambda$ のバリオン数がゼロである。
- [ ] $\Lambda$ は電荷正の中間子である。

:::explanation
$\Lambda$ は $uds$ のバリオンで電荷はゼロ、バリオン数は 1 である。核子とは異なる粒子種なので、同種核子の占有済み軌道に似た空間状態へ束縛され得る。
:::
::::

::::quiz{answer="0.2" tol="1e-6"}
$p(\phi)=[1+0.4\cos2\phi]/(2\pi)$ から基準面ゼロの $v_2$ を求めよ。
:::explanation
:eqref[eq-nuc-flow] で $\cos2\phi$ を掛けて積分する。定数項は消え、$\cos^2 2\phi$ の平均が $1/2$ なので $0.4/2=0.2$ となる。
:::
::::

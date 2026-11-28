---
title: α 崩壊と γ 遷移
description: α 放出のトンネル透過と、励起核の電磁遷移の強度・寿命を扱う。
notation:
  - symbol: 'Q_\alpha'
    meaning: α 崩壊で解放されるエネルギー
  - symbol: 'B(XL)'
    meaning: 種類 X、階数 L の換算遷移確率
  - symbol: '\alpha_{\rm IC}'
    meaning: 内部転換係数
  - symbol: 'f_{\rm LM}'
    meaning: 無反跳分率
---

## この章の目標

- 二体崩壊のエネルギー配分と $\alpha$ 崩壊の寿命を求める。
- Coulomb 障壁の WKB 積分から Geiger–Nuttall 則を説明する。
- 多重極の選択則と強度を使って準位の崩壊を判定する。
- 内部転換、核異性体、無反跳共鳴の仕組みを区別する。

## α 崩壊の系統性と Geiger–Nuttall 則

$\alpha$ 粒子は ${}^4\mathrm{He}$ の核であり、$J^\pi=0^+$ である。親核を静止させ、娘核と $\alpha$ の二体崩壊を考える。崩壊エネルギーは核質量により $Q_\alpha=(M_i-M_d-m_\alpha)c^2$ である。測定した中性原子質量を使う場合は初期と最終で電子静止質量が相殺し、電子束縛の差が残る。娘核が励起されるチャネルでは、その励起エネルギーを基底状態への $Q_\alpha$ から引く。

運動量が大きさ $p$ で反対向きなら、非相対論近似で $Q_\alpha=p^2/(2m_\alpha)+p^2/(2M_d)$ となる。したがって $\alpha$ の運動エネルギーは $Q_\alpha M_d/(M_d+m_\alpha)$ である。重い核では $Q_\alpha>0$ でも、このエネルギーで Coulomb 障壁を透過する必要がある。:term[Geiger–Nuttall 則]{en="Geiger–Nuttall law" yomi="がいがーぬったるそく"} は、娘核の電荷が等しい同位体鎖などで $\log_{10}(T_{1/2}/1\,\mathrm s)$ が $1/\sqrt{Q_\alpha}$ の一次式で近似できるという経験則である。ここで $T_{1/2}$ は対象の分岐の部分半減期である。娘核の電荷や核構造が変われば、直線の傾きや切片も変わり得る。

:::example[反跳への配分]
$M_d/m_\alpha=50$、$Q_\alpha=5.1\,\mathrm{MeV}$ なら $\alpha$ に $5.0\,\mathrm{MeV}$、娘核に $0.1\,\mathrm{MeV}$ が渡る。粒子線のエネルギーから $Q_\alpha$ を見積もるにはこの差を加える。
:::

:::example[励起状態への分岐]
基底状態への $Q_\alpha$ が $5\,\mathrm{MeV}$、娘核の励起が $0.5\,\mathrm{MeV}$ なら、この分岐の $Q$ は $4.5\,\mathrm{MeV}$ である。次節の透過率では、この小さな変化が指数に入る。
:::

### エネルギーと寿命の対数を測る

親が静止し、娘核の励起を $E_f$ とする。与えられた基底状態間の $Q_0$ を用いると、運動エネルギーの和は $Q_f=Q_0-E_f$ になる。娘核の励起を含む質量を $M_d^*=M_d+E_f/c^2$ とする。親と $\alpha$ 粒子の四元運動量を $P,p$ とし、計量を $(+,-,-,-)$ と取れば、$(P-p)^2=M_d^{*2}c^2$ である。親の静止系で $P\cdot p=M_PE_\alpha$ を代入すると、$E_\alpha=(M_P^2+m_\alpha^2-M_d^{*2})c^2/(2M_P)$ となる。$M_P=M_d^*+m_\alpha+Q_f/c^2$ を入れて展開すると、運動エネルギー $E_\alpha-m_\alpha c^2$ は $Q_f$ の一次で $Q_f M_d^*/(M_d^*+m_\alpha)$ となる。したがって粒子の線エネルギーに娘核の反跳エネルギーを加えると、質量表から得られる崩壊エネルギーに対応する。

Geiger–Nuttall の図で同位体群を比べるときは、全半減期から目的の分岐の部分半減期を求める。$\alpha$ 崩壊の娘に複数の準位があり、領域ごとに軌道角運動量が違うと、同じ $Q$ に換算しても障壁と形成確率に差が残る。

:::example[α 粒子の線エネルギー]
核の質量比を $r=M_d^*/m_\alpha$ とし、反跳の非相対論近似で $Q_f$ の分配を求める。

計算には次の値を用いる。
$$
r=50,\qquad Q_{f}=\frac{26}{5}\,\mathrm{MeV}
$$

この条件で求める量を ${T_\alpha}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{T_\alpha}&=\frac{Q_{f} \cdot r}{r + 1}\\&=\frac{\left(50\right) \cdot \left(\frac{26}{5}\right)}{\left(50\right) + 1}\\&\simeq 5.09804\,\mathrm{MeV}.\end{aligned}
$$

励起を伴う娘核への線は、運動に使える $Q_f$ が小さいため低い。

娘核の励起で窓が狭い分岐を見る。次の値を変え、ほかは固定する。
$$
Q_{f}=\frac{24}{5}\,\mathrm{MeV}
$$

このときの値を ${T_\alpha}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{T_\alpha}_{\mathrm{new}}&=\frac{\left(50\right) \cdot \left(\frac{24}{5}\right)}{\left(50\right) + 1}\\&\simeq 4.70588\,\mathrm{MeV},\\{T_\alpha}_{\mathrm{new}}-{T_\alpha}&\simeq -0.392157\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[部分半減期]
全半減期 $T$ と $\alpha$ 線の分岐確率 $b$ を使う。率は $b\ln2/T$ である。

計算には次の値を用いる。
$$
T=12\,\mathrm{s},\qquad b=\frac{4}{5}
$$

この条件で求める量を ${T_\alpha}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{T_\alpha}&=\frac{T}{b}\\&=\frac{\left(12\right)}{\left(\frac{4}{5}\right)}\\&=15\,\mathrm{s}.\end{aligned}
$$

部分半減期は全半減期より長い。分岐比で掛けるのは率である。

弱い線の部分半減期を比べる。次の値を変え、ほかは固定する。
$$
b=\frac{3}{10}
$$

このときの値を ${T_\alpha}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{T_\alpha}_{\mathrm{new}}&=\frac{\left(12\right)}{\left(\frac{3}{10}\right)}\\&=40\,\mathrm{s},\\{T_\alpha}_{\mathrm{new}}-{T_\alpha}&=25\,\mathrm{s}.\end{aligned}
$$
:::

:::example[同位体列での傾き]
局所的な当てはめ $y=a/\sqrt Q+b$ で、二つの $Q$ での常用対数の差を求める。

計算には次の値を用いる。
$$
a=120\,\mathrm{MeV}^{1/2},\qquad Q_{1}=5\,\mathrm{MeV},\qquad Q_{2}=\frac{101}{20}\,\mathrm{MeV}
$$

この条件で求める量を ${\Delta y}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\Delta y}&=a \cdot \left(1 \cdot \frac{1}{\sqrt{Q_{2}}} - 1 \cdot \frac{1}{\sqrt{Q_{1}}}\right)\\&=\left(120\right) \cdot \left(\frac{1}{\sqrt{\left(\frac{101}{20}\right)}} - \frac{1}{\sqrt{\left(5\right)}}\right)\\&\simeq -0.266332.\end{aligned}
$$

対数差の 1 は半減期が十倍違うことに対応する。同じ娘電荷の範囲を仮定している。

質量差が少し大きい例と比べる。次の値を変え、ほかは固定する。
$$
Q_{2}=\frac{53}{10}\,\mathrm{MeV}
$$

このときの値を ${\Delta y}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\Delta y}_{\mathrm{new}}&=\left(120\right) \cdot \left(\frac{1}{\sqrt{\left(\frac{53}{10}\right)}} - \frac{1}{\sqrt{\left(5\right)}}\right)\\&\simeq -1.54096,\\{\Delta y}_{\mathrm{new}}-{\Delta y}&\simeq -1.27463.\end{aligned}
$$

**一次近似の検証。** 質量表から得る $Q$ の誤差を半減期の常用対数に移すため、終側の窓を微分する。

${Q_{2}}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${\Delta y}$ の一階微分を求めると、
$$
\frac{\mathrm d{\Delta y}}{\mathrm dQ_{2}}=- \frac{a}{2 \cdot Q_{2}^{\frac{3}{2}}}.
$$

この導関数を変更前の点で評価して $h=\frac{1}{4}$ を掛ける。ここで $Q_{2}$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq -1.32176,\\\Delta-\Delta_{\mathrm{lin}}&\simeq 0.0471348.\end{aligned}
$$

残差が示すのは、指定した二点間で生じる線形化の誤差である。基準点からさらに離れた値について同じ精度を保証するものではない。
:::

:::example[反跳を独立に測る]
線の運動エネルギーを $T$、試料の質量比を $r$ とする。反跳の分を戻した窓を求める。

計算には次の値を用いる。
$$
T=5\,\mathrm{MeV},\qquad r=50
$$

この条件で求める量を ${Q_f}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{Q_f}&=T \cdot \left(1 + 1 \cdot \frac{1}{r}\right)\\&=\left(5\right) \cdot \left(1 + \frac{1}{\left(50\right)}\right)\\&\simeq 5.1\,\mathrm{MeV}.\end{aligned}
$$

軽い娘核ほど反跳が大きい。同じ核から出る二線の分離を測る場合にも補正は必要である。

娘が軽い別の模型と比べる。次の値を変え、ほかは固定する。
$$
r=20
$$

このときの値を ${Q_f}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{Q_f}_{\mathrm{new}}&=\left(5\right) \cdot \left(1 + \frac{1}{\left(20\right)}\right)\\&\simeq 5.25\,\mathrm{MeV},\\{Q_f}_{\mathrm{new}}-{Q_f}&\simeq 0.15\,\mathrm{MeV}.\end{aligned}
$$
:::

## Gamow のトンネル理論

核力が届く内側の半径を $R$、外側の Coulomb ポテンシャルを $C/r$、$C=2Z_de^2/(4\pi\varepsilon_0)$ とする。換算質量 $\mu$ を使い、まず $\ell=0$ を考える。$Q=Q_\alpha$、外側転回点 $b=C/Q$ が $R$ より大きいとき、:term[透過率]{en="penetrability" yomi="とうかりつ"} の WKB 指数は
$$
P\simeq e^{-2G},\qquad
G=\frac{\sqrt{2\mu}}{\hbar}\int_R^b\sqrt{C/r-Q}\dd r .
$$
WKB の包絡線と流束は :eqref[eq-qm-wkb-amplitude-flux] を参照する。$r=b\sin^2\theta$ と置くと平方根は $\sqrt Q\cot\theta$ であり、被積分関数は $2b\sqrt Q\cos^2\theta\dd\theta$ になる。$x=R/b$、$v=\sqrt{2Q/\mu}$ とすれば
$$
G=\frac{2C}{\hbar v}
\left[\arccos\sqrt x-\sqrt{x(1-x)}\right].
\label{eq-nuc-gamow}
$$
$x\ll1$ で角括弧が $\pi/2-2\sqrt x+O(x^{3/2})$ である。$2G$ の主項は $2\pi C/(\hbar v)$、次項は $Q$ によらない。$\alpha$ を形成して障壁へ達する率を $P_\alpha\nu$ とすれば、$\lambda\simeq P_\alpha\nu e^{-2G}$ である。$P_\alpha$ は核内の構造を反映する形成確率、$\nu$ は衝突頻度であり、これらが緩やかに変わるなら :eqref[eq-nuc-decay-law] から $\log T_{1/2}$ の主項が $Q^{-1/2}$ となる。

角運動量保存から $|J_i-J_d|\le\ell\le J_i+J_d$、パリティ保存から $\pi_i=\pi_d(-1)^\ell$ が必要である。非ゼロの $\ell$ では外部ポテンシャルに $\hbar^2\ell(\ell+1)/(2\mu r^2)$ が加わり、禁止帯の幅と積分が変わる。このため $Q$ だけで分岐比を決められない。

:::example[指数の感度]
二状態で $2G$ の差が $\ln100$、その他の因子が共通なら、寿命比は 100 である。基底状態の小さな質量差でも、透過率が指数関数なので寿命の差が大きくなる。
:::

:::example[最低軌道角運動量]
$2^+\to0^+$ では角運動量から $\ell=2$ に限られ、パリティ条件も満たす。$0^-\to0^+$ では角運動量から $\ell=0$ だがパリティは変えられないので、この二体チャネルは閉じる。
:::

### 障壁の傾きと形成率

同じ Coulomb 係数 $C$ と換算質量 $\mu$ で、内側の半径 $R$ も固定してみる。:eqref[eq-nuc-gamow] の $x\ll1$ の近似で $\log_{10}T_{1/2}$ の主項は
$$
\frac{2\pi C}{\hbar\ln10}\sqrt{\frac{\mu}{2}}\,Q^{-1/2}
$$
である。小さい変化 $\Delta Q$ に対する主項の変化はこれに $-\Delta Q/(2Q)$ を掛けた値になる。実験で同じ核の二つの娘準位への強度を比べる場合は、各部分崩壊定数 $\lambda_i=b_i\ln2/T_{1/2}$ を使う。模型で $\nu_iP_i$ を推定したら、各分岐の $\lambda_i e^{2G_i}$ の比較が試行頻度と形成確率の積の違いを反映する。ただし $G_i$ の推定における半径や遠心項の違いも同じ指標に入る。

:::example[異なる分岐の比較]
ある二枝で部分率が等しく、枝 2 の $2G$ が枝 1 より $\ln10$ 大きいとする。模型では枝 2 の $P_\alpha\nu$ が枝 1 の 10 倍でなければ等しい率にならない。逆に形成率を共通と仮定した場合、枝 2 は枝 1 の $1/10$ の率となる。
:::

### 作用の差と角運動量の障壁

遠心項込みの有効ポテンシャルを $V_L(r)=C/r+\hbar^2L(L+1)/(2\mu r^2)$ とする。$C>0$ である。外側の転回点は $Q r^2-Cr-D=0$ の正根、$D=\hbar^2L(L+1)/(2\mu)$ である。$L$ を増やすと区間内部の $\sqrt{V_L-Q}$ と区間の上限がともに増える。ゆえに同じ形成確率・同じ $Q$ での WKB 透過率は下がる。角運動量とパリティ保存は $|I_i-I_f|\le L\le I_i+I_f$ および $\pi_i=\pi_f(-1)^L$ を課す。$\alpha$ はスピン 0 なので、これを満たす最小の $L$ が一次の障壁の目安になる。WKB の接続条件については :eqref[eq-qm-airy] を参照する。

:::example[外部転回点]
Coulomb 定数を $C$、遠心項の定数を $D$ とする。二次方程式の正根を用いる。

計算には次の値を用いる。
$$
C=100\,\mathrm{MeV\,fm},\qquad D=0\,\mathrm{MeV\,fm^2},\qquad Q=5\,\mathrm{MeV}
$$

この条件で求める量を ${b}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{b}&=\frac{C + \sqrt{C^{2} + 4 \cdot D \cdot Q}}{2 \cdot Q}\\&=\frac{\left(100\right) + \sqrt{4 \cdot \left(0\right) \cdot \left(5\right) + \left(100\right)^{2}}}{2 \cdot \left(5\right)}\\&=20\,\mathrm{fm}.\end{aligned}
$$

転回点は幾何学的な核半径とは別の長さである。

遠心項のある部分波と比較する。次の値を変え、ほかは固定する。
$$
D=120\,\mathrm{MeV\,fm^2}
$$

このときの値を ${b}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{b}_{\mathrm{new}}&=\frac{\left(100\right) + \sqrt{\left(100\right)^{2} + 4 \cdot \left(120\right) \cdot \left(5\right)}}{2 \cdot \left(5\right)}\\&\simeq 21.1355\,\mathrm{fm},\\{b}_{\mathrm{new}}-{b}&\simeq 1.13553\,\mathrm{fm}.\end{aligned}
$$
:::

:::example[Coulomb 作用の有限半径補正]
$x=R/b$ の係数 $\arccos\sqrt x-\sqrt{x(1-x)}$ を求める。

計算には次の値を用いる。
$$
x=\frac{1}{10}
$$

この条件で求める量を ${f}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{f}&=- \sqrt{x \cdot \left(1 - x\right)} + \arccos{\left(\sqrt{x} \right)}\\&=- \sqrt{\left(\frac{1}{10}\right) \cdot \left(1 - \left(\frac{1}{10}\right)\right)} + \arccos{\left(\sqrt{\left(\frac{1}{10}\right)} \right)}\\&\simeq 0.949046.\end{aligned}
$$

$x\to0$ で $\pi/2$ に近づくが、有限半径での差は指数に入る。

内側転回点が外にある模型と比べる。次の値を変え、ほかは固定する。
$$
x=\frac{3}{10}
$$

このときの値を ${f}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{f}_{\mathrm{new}}&=- \sqrt{\left(\frac{3}{10}\right) \cdot \left(1 - \left(\frac{3}{10}\right)\right)} + \arccos{\left(\sqrt{\left(\frac{3}{10}\right)} \right)}\\&\simeq 0.532899,\\{f}_{\mathrm{new}}-{f}&\simeq -0.416147.\end{aligned}
$$
:::

:::example[指数と寿命の増幅]
二つの模型の透過率を比べる。無次元作用 $G_1,G_2$ を使い、前因子を固定する。

計算には次の値を用いる。
$$
G_{1}=20,\qquad G_{2}=\frac{41}{2}
$$

この条件で求める量を ${T_2/T_1}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{T_2/T_1}&=e^{2 \cdot \left(- G_{1} + G_{2}\right)}\\&=e^{- 2 \cdot \left(20\right) + 2 \cdot \left(\frac{41}{2}\right)}\\&\simeq 2.71828.\end{aligned}
$$

作用差は寿命比の対数に比例する。大きな絶対値を個別に丸めてから比を取る必要はない。

作用の差が大きい場合と比べる。次の値を変え、ほかは固定する。
$$
G_{2}=22
$$

このときの値を ${T_2/T_1}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{T_2/T_1}_{\mathrm{new}}&=e^{- 2 \cdot \left(20\right) + 2 \cdot \left(22\right)}\\&\simeq 54.5982,\\{T_2/T_1}_{\mathrm{new}}-{T_2/T_1}&\simeq 51.8799.\end{aligned}
$$

**一次近似の検証。** 大きな作用差に一次近似を使えるかどうかを、寿命比そのものの微分から確かめる。

${G_{2}}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${T_2/T_1}$ の一階微分を求めると、
$$
\frac{\mathrm d{T_2/T_1}}{\mathrm dG_{2}}=2 \cdot e^{- 2 \cdot G_{1} + 2 \cdot G_{2}}.
$$

この導関数を変更前の点で評価して $h=\frac{3}{2}$ を掛ける。ここで $G_{2}$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq 8.15485,\\\Delta-\Delta_{\mathrm{lin}}&\simeq 43.725.\end{aligned}
$$

この変化では一次近似の誤差が大きいため、増分の計算に元の式を使う。
:::

:::example[形成率の抽出]
測定率 $\lambda$ を、一回の透過率 $P$ と試行数 $\nu$ の積で割る。

計算には次の値を用いる。
$$
lam=\frac{1}{100}\,\mathrm{s}^{-1},\qquad \nu=1000000000000000000000\,\mathrm{s}^{-1},\qquad P=\frac{1}{10000000000000000000000}
$$

この条件で求める量を ${P_f}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{P_f}&=\frac{lam}{P \cdot \nu}\\&=\frac{\left(\frac{1}{100}\right)}{\left(1000000000000000000000\right) \cdot \left(\frac{1}{10000000000000000000000}\right)}\\&\simeq 0.1.\end{aligned}
$$

これは指定した障壁と内部速度の模型に条件づけた形成率である。

障壁模型が異なる場合を比べる。次の値を変え、ほかは固定する。
$$
P=\frac{1}{2000000000000000000000}
$$

このときの値を ${P_f}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{P_f}_{\mathrm{new}}&=\frac{\left(\frac{1}{100}\right)}{\left(1000000000000000000000\right) \cdot \left(\frac{1}{2000000000000000000000}\right)}\\&\simeq 0.02,\\{P_f}_{\mathrm{new}}-{P_f}&\simeq -0.08.\end{aligned}
$$
:::

## 多重極放射と遷移確率

核の励起エネルギーが光子に移る場合、光子の角運動量を多重極で分類する。波数 $k=E_\gamma/(\hbar c)$ に対して $kR\ll1$ なら、電荷分布の低次のモーメントが主要項になる。電気多重極演算子は $\mathcal M(E L,\mu)=\int r^L Y_{L\mu}(\hat{\vb r})\rho(\vb r)\dd^3r$ である。磁気多重極は電流や磁化を用い、反転に対する符号は電気と逆になる。

:::definition[換算遷移確率]{#def-nuc-reduced-transition}
:term[換算遷移確率]{en="reduced transition probability" yomi="かんさんせんいかくりつ"} を
$$
B(XL;J_i\to J_f)=\frac{|\langle J_f\|\mathcal M(XL)\|J_i\rangle|^2}{2J_i+1}
$$
とする。$X=E,M$ である。磁気射影を和し初期射影を平均すると、Clebsch–Gordan 係数の直交性でこの形になる。
:::

長波長での電気遷移率は、$B(EL)$ を SI の電荷と長さで与える場合、
$$
\lambda_\gamma(EL)=
\frac{8\pi(L+1)}{L[(2L+1)!!]^2}\frac{k^{2L+1}}{\hbar}
 \frac{B(EL)}{4\pi\varepsilon_0}.
\label{eq-nuc-electric-rate}
$$
ここで $n!!$ は $n,n-2,\ldots$ の積である。光子の位相空間が $k^2$、モードの規格化が $1/k$、電気多重極振幅の長波長依存が $k^L$ を与えるので、全体は $k^{2L+1}$ である。角度の規格化を含む係数はこの節の末尾で導く。磁気率も対応する単位で同じ $k^{2L+1}$ 依存をもつ。電気率で $B/(4\pi\varepsilon_0)$ はエネルギーと長さの $2L+1$ 乗の積であり、式は $\mathrm{s}^{-1}$ になる。

光子が出ると核が反跳する。初期核の静止質量を $M_i$、終状態を $M_f$ とし、エネルギーと運動量保存から $M_i c^2=E_\gamma+\sqrt{M_f^2c^4+E_\gamma^2}$ を得る。従って $E_\gamma=(M_i^2-M_f^2)c^2/(2M_i)$ であり、準位差より約 $E_\gamma^2/(2M_fc^2)$ 小さい。

:::example[周波数による率の変化]
$B(E2)$ が共通で $E_\gamma$ が 2 倍になると率は $2^5=32$ 倍になる。結果として寿命は他の枝がなければ $1/32$ になる。
:::

:::example[換算強度の逆向き]
Hermite テンソルでは順逆の換算行列要素の絶対値が等しい。$J_i=2,J_f=0$ の下向きの $B$ は、上向きの $B$ の $1/5$ である。分母の違いは磁気射影についての平均から生じる。
:::

### 電気双極子の係数を確かめる

長波長の電気双極子相互作用を $H_{\rm int}=-\vb d\cdot\vb E$ とする。$\vb d=\int\vb r\rho(\vb r)\dd^3r$ である。体積 $\mathcal V$ に規格化した角周波数 $\omega_k$ の光子一個に対する電場の係数は $\sqrt{\hbar\omega_k/(2\varepsilon_0\mathcal V)}$ である。放射ゲージでは $\vb E=-\partial_t\vb A$ なので、:eqref[eq-aqm-vector-potential] の一光子係数に $\omega_k$ を掛ければ得られる。磁気射影を指定した核行列要素を $\vb d_{fi}$ とする。二つの横偏極の和は $|\vb d_{fi}|^2-|\hat{\vb k}\cdot\vb d_{fi}|^2$ であり、全球面の積分は $(8\pi/3)|\vb d_{fi}|^2$ である。これは $\int\hat k_i\hat k_j\dd\Omega=(4\pi/3)\delta_{ij}$ から分かる。

:ref[law-qm-golden] で、状態密度 $\mathcal V k^2\dd k\dd\Omega/(2\pi)^3$ とデルタ関数 $\delta(E_{\rm tr}-\hbar ck)$ を掛ける。動径積分は $k^2\omega_k/(\hbar c)$ を与え、体積は規格化因子と消える。従って率は $k^3|\vb d_{fi}|^2/(3\pi\varepsilon_0\hbar)$ となる。初期射影の平均と終射影の和を取ると、$\vb d$ の球成分と $\mathcal M(E1,\mu)$ の比は $\sqrt{4\pi/3}$ なので、$|\vb d_{fi}|^2$ の平均は $(4\pi/3)B(E1)$ になる。こうして :eqref[eq-nuc-electric-rate] の $L=1$ を再現する。

:::example[一つの線の寿命]
$B(E1)=10^{-3}e^2\mathrm{fm}^2$、$E_\gamma=1\,\mathrm{MeV}$ とする。$e^2/(4\pi\varepsilon_0)\simeq1.440\,\mathrm{MeV\,fm}$、$\hbar\simeq6.582\times10^{-22}\,\mathrm{MeV\,s}$ を使うと、$\lambda_\gamma(E1)\simeq1.59\times10^{12}\,\mathrm s^{-1}$ である。この枝だけなら平均寿命は約 $0.63\,\mathrm{ps}$ になる。
:::

### 多重極を長波長で展開する

角度部分を分離した平面波では、動径関数として球 Bessel 関数 $j_L(kr)$ が現れる。方程式 $x^2j''+2xj'+[x^2-L(L+1)]j=0$ に $j=x^L\sum_{n\ge0}a_nx^{2n}$ を入れると $a_{n+1}=-a_n/[2(n+1)(2L+2n+3)]$ である。$L=0$ の $j_0=\sin x/x$ と梯子 $j_{L+1}=Lj_L/x-j_L'$ を使えば $a_0=1/(2L+1)!!$ となる。この動径の最初の項が長波長近似に現れ、次の相対補正の大きさは $(kr)^2/[2(2L+3)]$ で見積もれる。

換算遷移確率の順逆では $B(i\to f)$ と $B(f\to i)$ は一般に違う。Hermite な球面テンソルでは、射影ごとの絶対値二乗を全て足した値は随伴性で同じであり、初期射影で割る数が違うことから $(2I_i+1)B(i\to f)=(2I_f+1)B(f\to i)$ が得られる。

:::example[長波長パラメータ]
光子の波数を $E/h_c$ とし、半径 $R$ に対する $kR$ を調べる。

計算には次の値を用いる。
$$
E=1\,\mathrm{MeV},\qquad R=6\,\mathrm{fm},\qquad h_{c}=\frac{197327}{1000}\,\mathrm{MeV\,fm}
$$

この条件で求める量を ${kR}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{kR}&=\frac{E \cdot R}{h_{c}}\\&=\frac{\left(1\right) \cdot \left(6\right)}{\left(\frac{197327}{1000}\right)}\\&\simeq 0.0304064.\end{aligned}
$$

この無次元数で動径の高次項を省けるか判断する。

高い光子エネルギーを比べる。次の値を変え、ほかは固定する。
$$
E=5\,\mathrm{MeV}
$$

このときの値を ${kR}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{kR}_{\mathrm{new}}&=\frac{\left(5\right) \cdot \left(6\right)}{\left(\frac{197327}{1000}\right)}\\&\simeq 0.152032,\\{kR}_{\mathrm{new}}-{kR}&\simeq 0.121626.\end{aligned}
$$
:::

:::example[同じ E2 強度での率の差]
強度が等しい二本の $E2$ でエネルギーを $E_1,E_2$ とする。

計算には次の値を用いる。
$$
E_{1}=\frac{1}{2}\,\mathrm{MeV},\qquad E_{2}=\frac{11}{20}\,\mathrm{MeV}
$$

この条件で求める量を ${\lambda_2/\lambda_1}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\lambda_2/\lambda_1}&=\left(\frac{E_{2}}{E_{1}}\right)^{5}\\&=\frac{\left(\frac{11}{20}\right)^{5}}{\left(\frac{1}{2}\right)^{5}}\\&\simeq 1.61051.\end{aligned}
$$

寿命を比べるには逆数の比を使う。

エネルギー差を増やす。次の値を変え、ほかは固定する。
$$
E_{2}=\frac{4}{5}\,\mathrm{MeV}
$$

このときの値を ${\lambda_2/\lambda_1}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\lambda_2/\lambda_1}_{\mathrm{new}}&=\frac{\left(\frac{4}{5}\right)^{5}}{\left(\frac{1}{2}\right)^{5}}\\&\simeq 10.4858,\\{\lambda_2/\lambda_1}_{\mathrm{new}}-{\lambda_2/\lambda_1}&\simeq 8.87525.\end{aligned}
$$
:::

:::example[吸収と放出の規格化]
基底スピンをゼロとして、励起 $J$ から基底への $B_\downarrow$ から逆向きの $B$ を求める。

計算には次の値を用いる。
$$
J=2,\qquad B_{d}=5
$$

この条件で求める量を ${B_\uparrow}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{B_\uparrow}&=B_{d} \cdot \left(2 \cdot J + 1\right)\\&=\left(5\right) \cdot \left(2 \cdot \left(2\right) + 1\right)\\&=25.\end{aligned}
$$

強度の単位はそれぞれの多重極に対応した共通の単位で与えている。

別の多重極で縮退度の効果を比べる。次の値を変え、ほかは固定する。
$$
J=3
$$

このときの値を ${B_\uparrow}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{B_\uparrow}_{\mathrm{new}}&=\left(5\right) \cdot \left(2 \cdot \left(3\right) + 1\right)\\&=35,\\{B_\uparrow}_{\mathrm{new}}-{B_\uparrow}&=10.\end{aligned}
$$
:::

:::example[E1 の部分率]
$B(E1)=b\,e^2\mathrm{fm}^2$、$a_e=e^2/(4\pi\epsilon_0)$ を使い、一般式の $L=1$ の係数 $16\pi/9$ を用いる。

計算には次の値を用いる。
$$
a_{e}=\frac{36}{25}\,\mathrm{MeV\,fm},\qquad b=\frac{1}{1000},\qquad E=1\,\mathrm{MeV},\qquad h_{c}=\frac{197327}{1000}\,\mathrm{MeV\,fm},\qquad \hbar=\frac{164553}{250000000000000000000000000}\,\mathrm{MeV\,s}
$$

この条件で求める量を ${\lambda}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\lambda}&=\frac{16 \cdot \pi \cdot a_{e} \cdot b \cdot \left(\frac{E}{h_{c}}\right)^{3}}{9 \cdot \hbar}\\&=\frac{16 \cdot \pi \cdot \left(1\right)^{3} \cdot \left(\frac{1}{1000}\right) \cdot \left(\frac{36}{25}\right)}{9 \cdot \left(\frac{164553}{250000000000000000000000000}\right) \cdot \left(\frac{197327}{1000}\right)^{3}}\\&\simeq 1.59025\times10^{12}\,\mathrm{s}^{-1}.\end{aligned}
$$

比較は内部転換なしの光子率で行う。

線エネルギーを二倍にする。次の値を変え、ほかは固定する。
$$
E=2\,\mathrm{MeV}
$$

このときの値を ${\lambda}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\lambda}_{\mathrm{new}}&=\frac{16 \cdot \pi \cdot \left(2\right)^{3} \cdot \left(\frac{1}{1000}\right) \cdot \left(\frac{36}{25}\right)}{9 \cdot \left(\frac{164553}{250000000000000000000000000}\right) \cdot \left(\frac{197327}{1000}\right)^{3}}\\&\simeq 1.2722\times10^{13}\,\mathrm{s}^{-1},\\{\lambda}_{\mathrm{new}}-{\lambda}&\simeq 1.11317\times10^{13}\,\mathrm{s}^{-1}.\end{aligned}
$$

**一次近似の検証。** 光子率の一階の相対感度はエネルギーの三乗に由来する。有限差での誤差も求める。

${E}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${\lambda}$ の一階微分を求めると、
$$
\frac{\mathrm d{\lambda}}{\mathrm dE}=\frac{16 \cdot \pi \cdot E^{2} \cdot a_{e} \cdot b}{3 \cdot h_{c}^{3} \cdot \hbar}.
$$

この導関数を変更前の点で評価して $h=1$ を掛ける。ここで $E$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq 4.77074\times10^{12}\,\mathrm{s}^{-1},\\\Delta-\Delta_{\mathrm{lin}}&\simeq 6.36099\times10^{12}\,\mathrm{s}^{-1}.\end{aligned}
$$

この変化では一次近似の誤差が大きいため、増分の計算に元の式を使う。
:::

### 横波の展開から電気多重極率を導く

電流の核行列要素を $\vb J(\vb r)$、電荷の行列要素を $\rho(\vb r)$ とし、調和時間依存を $e^{-i\omega t}$ とする。量子化した平面波のベクトルポテンシャルの一光子係数は $\sqrt{\hbar/(2\varepsilon_0\mathcal V\omega)}$ である。相互作用 $-\int\vb J\cdot\vb A\,\dd^3r$ を黄金律に代入し、光子の $k$ の積分に $\delta(E_i-E_f-\hbar ck)$ を使う。係数を整理すると、射影を指定した率は
$$
\lambda=\frac{k}{8\pi^2\varepsilon_0\hbar c^2}
\int\sum_{\alpha=1,2}\abs{\vb\epsilon_\alpha^*\cdot\widetilde{\vb J}(\vb k)}^2\dd\Omega_k,\qquad
\widetilde{\vb J}(\vb k)=\int e^{-i\vb k\cdot\vb r}\vb J(\vb r)\dd^3r.
$$
従って必要なのは Fourier 電流の横成分の球面ノルムである。

単位球面の角運動量を無次元 $\vb L=-i\vb r\times\nabla$ で表し、$\vb X_{LM}=\vb L Y_{LM}/\sqrt{L(L+1)}$ を取る。球面 Laplacian の部分積分からこの族のノルムは 1 である。$\hat{\vb k}\times\vb X_{LM}$ も同じノルムを持ち、前者との内積は勾配と回転勾配の内積の球面積分でゼロになる。横空間への展開のうち $\hat{\vb k}\times\vb X_{LM}$ 方向の係数が電気成分である。$\vb X_{LM}$ の角積分は、スカラー平面波展開を $\vb L$ で微分し、球面上で部分積分すると $4\pi j_L(kr)\vb X_{LM}(\hat{\vb r})$ になる。ただし全体の $i$ の位相は率に寄与しない。$\hat{\vb k}$ との外積は Fourier 核の $\vb r$ に関する回転を $k$ で割る操作に置き換えられる。この展開では、電気成分は :term[ベクトル球面調和関数]{en="vector spherical harmonic" yomi="べくとるきゅうめんちょうわかんすう"} の回転を含む。従って各成分の係数の大きさは
$$
4\pi\abs{\int\vb J\cdot
\frac{\nabla\times[j_L(kr)\vb X_{LM}^*(\hat{\vb r})]}{k}\dd^3r}
$$
である。

$S_{LM}=r^LY_{LM}$ は調和な斉次多項式である。$\nabla\times(\vb r\times\nabla S_{LM})=-(L+1)\nabla S_{LM}$ は、$\nabla^2 S_{LM}=0$ と $(\vb r\cdot\nabla)\nabla S_{LM}=(L-1)\nabla S_{LM}$ を外積の微分公式に入れれば得られる。そこで $j_L\simeq (kr)^L/(2L+1)!!$ を入れると、回転を含む因子は位相を除き $\sqrt{(L+1)/L}\,k^{L-1}\nabla S_{LM}^*/(2L+1)!!$ となる。連続の式 $\nabla\cdot\vb J=i\omega\rho$ を使い部分積分すると $\int\vb J\cdot\nabla S_{LM}^*=-i\omega\int\rho S_{LM}^*$ である。境界項は局在した電流で消える。従って長波長での電気成分の二乗は
$$
(4\pi)^2\frac{L+1}{L[(2L+1)!!]^2}
k^{2L-2}\omega^2\abs{\int\rho S_{LM}^*\,\dd^3r}^2.
$$
$\omega=ck$ とし、磁気射影の和と平均を取ると :eqref[eq-nuc-electric-rate] の係数に一致する。$L=0$ では横調和関数が存在しないため、この計算は $L\ge1$ に限る。動径展開の次の項は $(kR)^2$ 高次だが、先頭の核行列要素が抑制された遷移では率の相対誤差をこの比だけから見積もれない。有限波長の電流演算子や、電荷密度との連続の式を満たす多体電流の扱いは :cite[blatt-weisskopf] を参照する。

### 磁気多重極の電流演算子

磁気成分は同じ横波展開の $\vb X_{LM}$ 方向である。Fourier 電流をこの方向へ射影すると、位相を除いた係数は $4\pi\int j_L(kr)\vb X_{LM}^*\cdot\vb J\,\dd^3r$ である。$S_{LM}=r^L Y_{LM}$ と書き、$r^L\vb L Y_{LM}=i(\nabla S_{LM})\times\vb r$ を用いる。長波長で積分の大きさは $k^L|\int(\nabla S_{LM}^*)\cdot(\vb r\times\vb J)\,\dd^3r|/\bigl((2L+1)!!\sqrt{L(L+1)}\bigr)$ となる。したがって磁気演算子を SI で

$$
\mathcal M(ML,M)=\frac1{L+1}\int
(\nabla S_{LM})\cdot(\vb r\times\vb J)\,\dd^3r
$$

と定めると、磁気率は :eqref[eq-nuc-electric-rate] の $B(EL)/(4\pi\varepsilon_0)$ を $\mu_0B(ML)/(4\pi)$ に置き換えた形になる。これは電流の係数を二乗すると $(L+1)/L$ が現れ、残る $1/(\varepsilon_0c^2)$ が $\mu_0$ に等しいためである。演算子の次元は磁気モーメントと長さの $L-1$ 乗の積である。

:::example[磁気双極子の規格化]
$L=1$ では $\nabla S_{1M}$ は球成分を選ぶ一定のベクトルである。磁気モーメント $\vb\mu=\tfrac12\int\vb r\times\vb J\,\dd^3r$ に対し、$\mathcal M(M1,M)=\sqrt{3/(4\pi)}\,\mu_M$ となる。電子電荷の大きさを $e$、陽子質量を $m_p$ とすると、磁気強度の単位に使う核磁子は $\mu_N=e\hbar/(2m_p)$ である。
:::

## Weisskopf 単位と選択則

### 単一粒子の基準

核内に一様に分布した単一粒子を考え、半径を $R=1.2A^{1/3}\,\mathrm{fm}$ とする。動径モーメントは $\int_0^R r^{L+2}\dd r/\int_0^R r^2\dd r=3R^L/(L+3)$ である。電気多重極の単一粒子の基準をこのモーメントの二乗から
$$
B_{\rm W}(EL)=\frac{1}{4\pi}\left(\frac{3}{L+3}\right)^2e^2R^{2L}
$$
と置く。これは :term[Weisskopf 単位]{en="Weisskopf unit" yomi="わいすこっふたんい"} である。磁気の場合の標準的な基準は $B_{\rm W}(ML)=10(3/(L+3))^2 R^{2L-2}\mu_N^2/\pi$ である。基準は特定の核の一粒子波動関数を解いた結果ではない。比較には電気と磁気でそれぞれ適切な単位を使う。

### 幾何学的な制限

初期角運動量は終状態と光子の合成に含まれる必要がある。従って $|J_i-J_f|\le L\le J_i+J_f$、実光子では $L\ge1$ である。電気では $\pi_i\pi_f=(-1)^L$、磁気では $(-1)^{L+1}$ である。前者は $Y_{L\mu}$ の反転から、後者は電流を含む磁気演算子の軸性から従う。条件は行列要素がゼロである場合を除くためのものであり、許された項の大きさまでは決めない。

:::example[混合多重極]
$3/2^+\to1/2^+$ では $L=1,2$ が候補であり、同パリティなので $M1,E2$ が許される。実験の角度相関で両振幅の混合比を調べられる。
:::

### 角度分布と偏光

初期核の磁気副準位を等しく重み付けするなら、終状態の射影と光子の偏光を足した分布は回転によって変わらない。したがって $\dd\lambda/\dd\Omega=\lambda/(4\pi)$ である。一方、準位の占有が軸について偏っていれば、光子角度と偏光にその情報が残る。前の遷移の光子を特定方向に観測して中間準位を作る場合にも、選んだ方向が次の光子の角度分布の基準になる。これが角度相関を用いて準位のスピンや多重極を調べられる理由である。

:::example[磁気副準位を指定した電気双極子放射]
$J_i=1,M_i=0\to J_f=0,M_f=0$ の $E1$ を考える。角運動量の射影が保存されるため双極子の球成分は $M=0$ のみであり、$\vb d_{fi}$ は量子化軸に平行である。放出方向が軸となす角を $\theta$ とすると、横偏極への射影の二乗は $|\vb d_{fi}|^2\sin^2\theta$ である。$\int\sin^2\theta\,\dd\Omega=8\pi/3$ を用いれば、方向の規格化密度は $(3/8\pi)\sin^2\theta$ となる。軸方向には放出されず、垂直方向では偏光は量子化軸に平行である。
:::

:::example[四重極の基準]
$A=125$ なら $R=6\,\mathrm{fm}$ で、$B_{\rm W}(E2)=0.36\cdot6^4/(4\pi)\,e^2\mathrm{fm}^4\simeq37.1\,e^2\mathrm{fm}^4$ である。実測値がこの 20 倍なら 20 W.u. と表す。
:::

### 単一粒子の尺度と選択則の併用

Weisskopf の電気単位では、半径 $R$ の内部に広がる一粒子の規格化動径分布を $3r^2/R^3$ で近似する。この平均で $\langle r^L\rangle=3R^L/(L+3)$ となり、角度の慣用規格化 $1/(4\pi)$ と掛けて $B_{\rm W}(EL)$ を作る。これは核子一個の尺度の定義である。複数の核子が同位相に寄与すると行列要素は振幅で増えるから、集団回転で単位を大きく超えることがある。同じ標準値以下でも、選択則で厳密にゼロになる場合と近似対称性で小さくなる場合を区別する。

:::example[E2 の Weisskopf 単位]
半径 $R=r_0A^{1/3}$ を用い、$B_{\rm W}(E2)/e^2$ を求める。

計算には次の値を用いる。
$$
r_{0}=\frac{6}{5}\,\mathrm{fm},\qquad A=64
$$

この条件で求める量を ${B_{\rm W}/e^2}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{B_{\rm W}/e^2}&=\frac{9 \cdot A^{\frac{4}{3}} \cdot r_{0}^{4}}{100 \cdot \pi}\\&=9 \cdot \frac{1}{100 \cdot \pi} \cdot \left(64\right)^{\frac{4}{3}} \cdot \left(\frac{6}{5}\right)^{4}\\&\simeq 15.2075\,\mathrm{fm}^4.\end{aligned}
$$

核種により単位の大きさが異なる。物理単位の強度は変えず単位だけ変えた比較もできる。

重い核の標準値と比較する。次の値を変え、ほかは固定する。
$$
A=125
$$

このときの値を ${B_{\rm W}/e^2}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{B_{\rm W}/e^2}_{\mathrm{new}}&=9 \cdot \frac{1}{100 \cdot \pi} \cdot \left(125\right)^{\frac{4}{3}} \cdot \left(\frac{6}{5}\right)^{4}\\&\simeq 37.1277\,\mathrm{fm}^4,\\{B_{\rm W}/e^2}_{\mathrm{new}}-{B_{\rm W}/e^2}&\simeq 21.9202\,\mathrm{fm}^4.\end{aligned}
$$

**一次近似の検証。** 滑らかな質量数依存の式を補間し、二核種の違いを近傍の傾きで説明できるか調べる。

${A}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${B_{\rm W}/e^2}$ の一階微分を求めると、
$$
\frac{\mathrm d{B_{\rm W}/e^2}}{\mathrm dA}=\frac{3 \cdot \sqrt[3]{A} \cdot r_{0}^{4}}{25 \cdot \pi}.
$$

この導関数を変更前の点で評価して $h=61$ を掛ける。ここで $A$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq 19.3262\,\mathrm{fm}^4,\\\Delta-\Delta_{\mathrm{lin}}&\simeq 2.59399\,\mathrm{fm}^4.\end{aligned}
$$

残差が示すのは、指定した二点間で生じる線形化の誤差である。基準点からさらに離れた値について同じ精度を保証するものではない。
:::

:::example[集団性を示す無次元数]
実測強度 $b$ を同じ多重極の単位 $b_W$ で割る。

計算には次の値を用いる。
$$
b=250\,\mathrm{e^2\,fm^4},\qquad b_{W}=\frac{25}{2}\,\mathrm{e^2\,fm^4}
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=\frac{b}{b_{W}}\\&=\frac{\left(250\right)}{\left(\frac{25}{2}\right)}\\&=20.\end{aligned}
$$

同じ大きさの実測強度でも質量の異なる核では比が変わる。

標準値を二倍にする。次の値を変え、ほかは固定する。
$$
b_{W}=25\,\mathrm{e^2\,fm^4}
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=\frac{\left(250\right)}{\left(25\right)}\\&=10,\\{R}_{\mathrm{new}}-{R}&=-10.\end{aligned}
$$
:::

:::example[多重極の下限]
スピン差から三角条件の下限を求める。ここでは $J_i>J_f$ とする。

計算には次の値を用いる。
$$
J_{i}=\frac{7}{2},\qquad J_{f}=\frac{3}{2}
$$

この条件で求める量を ${L_{\min}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{L_{\min}}&=- J_{f} + J_{i}\\&=- \left(\frac{3}{2}\right) + \left(\frac{7}{2}\right)\\&=2.\end{aligned}
$$

電気か磁気かはパリティから選ぶ。下限より高い階数も上限まで許されうる。

終状態のスピンを下げる。次の値を変え、ほかは固定する。
$$
J_{f}=\frac{1}{2}
$$

このときの値を ${L_{\min}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{L_{\min}}_{\mathrm{new}}&=- \left(\frac{1}{2}\right) + \left(\frac{7}{2}\right)\\&=3,\\{L_{\min}}_{\mathrm{new}}-{L_{\min}}&=1.\end{aligned}
$$
:::

:::example[固定強度での阻害の評価]
率の Weisskopf 推定を $\lambda_W$、測定した光子率を $\lambda$ とする。遷移エネルギーが同じなら逆比が阻害係数になる。

計算には次の値を用いる。
$$
lam_{W}=100000000\,\mathrm{s}^{-1},\qquad lam=1000000\,\mathrm{s}^{-1}
$$

この条件で求める量を ${F_W}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{F_W}&=\frac{lam_{W}}{lam}\\&=\frac{\left(100000000\right)}{\left(1000000\right)}\\&=100.\end{aligned}
$$

内部転換を含む全率を使うとこの比較は変わる。

光子率がさらに低い線と比べる。次の値を変え、ほかは固定する。
$$
lam=100000\,\mathrm{s}^{-1}
$$

このときの値を ${F_W}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{F_W}_{\mathrm{new}}&=\frac{\left(100000000\right)}{\left(100000\right)}\\&=1000,\\{F_W}_{\mathrm{new}}-{F_W}&=900.\end{aligned}
$$
:::

## 内部転換と核異性体

励起核と原子電子の電磁相互作用によって束縛電子が直接放出される過程を :term[内部転換]{en="internal conversion" yomi="ないぶてんかん"} という。殻 $i$ の束縛エネルギーを $B_i$ とすると、電子の運動エネルギーは反跳を省き $E_{\rm tr}-B_i$ である。空孔の緩和では特性 X 線や Auger 電子が生じる。光子が外で電子に吸収される場合とは発生過程が異なる。

一つの遷移で $\alpha_{\rm IC}=\lambda_e/\lambda_\gamma$ と置けば、光子放出と内部転換だけが開いているとき $\lambda_{\rm tr}=\lambda_\gamma(1+\alpha_{\rm IC})$ である。$0^+\to0^+$ では実光子一個を放出できないため、電子と結合する電気単極子 $E0$ が重要になる。準位差が $2m_ec^2$ を超えると、反跳を除いた運動学では内部対生成も開く。これが寄与する場合は $\lambda_{\rm tr}$ に対生成率を加える。$\lambda_\gamma=0$ のときは通常の転換係数の比は定義しない。

$E0$ 演算子の先頭が電荷の総和だけなら、異なる準位間の振幅は直交性で消える。球対称な電荷の Fourier 成分を展開すると、$j_0(qr)=1-q^2r^2/6+O(q^4r^4)$ であり、最初に残るのは $\langle f|\int r^2\rho(\vb r)\,\dd^3r|i\rangle$ である。したがって $E0$ の強度は、電荷分布の二乗半径が異なる配置の混合に感度を持つ。

高い $L$、小さな準位差、固有構造の違いは遷移率を小さくし、:term[核異性体]{en="nuclear isomer" yomi="かくいせいたい"} を作る。軸対称核で固有軸射影 $K$ がよい近似なら、固有座標の多重極の成分は $|\mu|\le L$ なので $|\Delta K|>L$ の振幅は消える。$K$ の混合があれば弱い遷移が起こり得る。

:::example[転換電子の線]
準位差 $100\,\mathrm{keV}$、束縛エネルギー $10\,\mathrm{keV}$ の殻なら電子線は約 $90\,\mathrm{keV}$ である。複数の殻から出ると束縛エネルギーの差だけ線が分かれる。
:::

:::example[寿命と分岐]
$\lambda_\gamma=100\,\mathrm{s}^{-1}$、$\alpha_{\rm IC}=4$ の一遷移なら全率は $500\,\mathrm{s}^{-1}$ である。光子の分岐は $1/5$、平均寿命は $2\,\mathrm{ms}$ となる。
:::

### 線の強度から率を区別する

同じ初期準位から二つの出口 $i=1,2$ へ電磁遷移する。内部転換係数を $\alpha_i$、検出器効率を $\epsilon_i$、計数を $C_i$ とする。同じ初期準位が供給された回数を $N_{\rm init}$ とする。各線の計数は $N_{\rm init}\epsilon_i\lambda_{\gamma,i}/\lambda_{\rm all}$ に比例する。このとき検出効率補正後の光子数の比は $\lambda_{\gamma,1}/\lambda_{\gamma,2}$ である。電子まで含めた二枝の率の比は、さらに $(1+\alpha_1)/(1+\alpha_2)$ を掛ける。別の枝や核の粒子放出があっても、同じ初期準位から出る二枝の比では $\lambda_{\rm all}$ が消える。

:::example[転換で変わる分岐]
二つの光子線の効率補正済み計数が等しく、$\alpha_1=3,\alpha_2=0$ なら、全遷移率の比は 4 である。枝がこの二つだけなら枝 1 の割合は $4/5$ となる。光子の比をそのまま使うと全遷移の割合を誤る。
:::

### 内部転換の線と競合経路

核から原子軌道の電子へエネルギーが渡ると、軌道束縛 $B_i$ を差し引いた運動エネルギーが残る。反跳を無視すれば $T_i=E_{\rm tr}-B_i$ であり、異なる殻の線の差は束縛の差である。転換係数は電気的な遷移の空間分布と電子の軌道波動関数に依存し、光子が一度放たれてから吸収される二段階過程ではない。転換後の原子の空孔は X 線や Auger 電子で緩和し、これらの線は核遷移光子と起源が異なる。$E0$ では単一光子の放出がないため $\alpha_i=\lambda_i/\lambda_\gamma$ という定義を使えず、電子や対生成の率そのものを扱う。

:::example[転換線の運動エネルギー]
同じ遷移エネルギーに対し、電子軌道の正の束縛を $B_i$ とする。

計算には次の値を用いる。
$$
E=150\,\mathrm{keV},\qquad B_{i}=35\,\mathrm{keV}
$$

この条件で求める量を ${T_i}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{T_i}&=- B_{i} + E\\&=\left(150\right) - \left(35\right)\\&=115\,\mathrm{keV}.\end{aligned}
$$

殻の内側の電子は束縛が深いので線が低い。

外側の殻と比較する。次の値を変え、ほかは固定する。
$$
B_{i}=7\,\mathrm{keV}
$$

このときの値を ${T_i}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{T_i}_{\mathrm{new}}&=\left(150\right) - \left(7\right)\\&=143\,\mathrm{keV},\\{T_i}_{\mathrm{new}}-{T_i}&=28\,\mathrm{keV}.\end{aligned}
$$
:::

:::example[内部転換を含む寿命]
光子だけの平均寿命を $t_\gamma$ とし、総転換係数を $a$ とする。競合する粒子崩壊はない。

計算には次の値を用いる。
$$
t_{g}=200\,\mathrm{ns},\qquad a=\frac{1}{2}
$$

この条件で求める量を ${\tau}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\tau}&=\frac{t_{g}}{a + 1}\\&=\frac{\left(200\right)}{\left(\frac{1}{2}\right) + 1}\\&\simeq 133.333\,\mathrm{ns}.\end{aligned}
$$

大きな転換係数でも光子部分率そのものは変えずに全率を増やす。

電子チャネルが強い場合を比較する。次の値を変え、ほかは固定する。
$$
a=4
$$

このときの値を ${\tau}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\tau}_{\mathrm{new}}&=\frac{\left(200\right)}{\left(4\right) + 1}\\&=40\,\mathrm{ns},\\{\tau}_{\mathrm{new}}-{\tau}&\simeq -93.3333\,\mathrm{ns}.\end{aligned}
$$

**一次近似の検証。** 転換係数の推定の誤差を全寿命に伝える。

${a}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${\tau}$ の一階微分を求めると、
$$
\frac{\mathrm d{\tau}}{\mathrm da}=- \frac{t_{g}}{\left(a + 1\right)^{2}}.
$$

この導関数を変更前の点で評価して $h=\frac{7}{2}$ を掛ける。ここで $a$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq -311.111\,\mathrm{ns},\\\Delta-\Delta_{\mathrm{lin}}&\simeq 217.778\,\mathrm{ns}.\end{aligned}
$$

この変化では一次近似の誤差が大きいため、増分の計算に元の式を使う。
:::

:::example[殻ごとの部分率]
総転換係数を $a$、一つの殻の係数を $a_i$ とし、全遷移率からその部分を分ける。

計算には次の値を用いる。
$$
a=2,\qquad a_{i}=\frac{3}{5}
$$

この条件で求める量を ${b_i}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{b_i}&=\frac{a_{i}}{a + 1}\\&=\frac{\left(\frac{3}{5}\right)}{\left(2\right) + 1}\\&\simeq 0.2.\end{aligned}
$$

分母には光子の 1 を含む。転換事象に条件づけると分母は $a$ になる。

同じ総係数で殻の寄与を変える。次の値を変え、ほかは固定する。
$$
a_{i}=\frac{6}{5}
$$

このときの値を ${b_i}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{b_i}_{\mathrm{new}}&=\frac{\left(\frac{6}{5}\right)}{\left(2\right) + 1}\\&\simeq 0.4,\\{b_i}_{\mathrm{new}}-{b_i}&\simeq 0.2.\end{aligned}
$$
:::

:::example[効率補正つきの転換係数]
電子と光子の計数を $C_e,C_g$、各効率を $\epsilon_e,\epsilon_g$ とする。

計算には次の値を用いる。
$$
C_{e}=600,\qquad C_{g}=1000,\qquad eps_{e}=\frac{1}{5},\qquad eps_{g}=\frac{2}{5}
$$

この条件で求める量を ${a}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{a}&=\frac{C_{e} \cdot eps_{g}}{C_{g} \cdot eps_{e}}\\&=\frac{\left(600\right) \cdot \left(\frac{2}{5}\right)}{\left(1000\right) \cdot \left(\frac{1}{5}\right)}\\&\simeq 1.2.\end{aligned}
$$

同じ線の組でも受光・受電子条件を変えると生の計数比の解釈が変わる。

電子効率の評価を変える。次の値を変え、ほかは固定する。
$$
eps_{e}=\frac{3}{10}
$$

このときの値を ${a}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{a}_{\mathrm{new}}&=\frac{\left(600\right) \cdot \left(\frac{2}{5}\right)}{\left(1000\right) \cdot \left(\frac{3}{10}\right)}\\&\simeq 0.8,\\{a}_{\mathrm{new}}-{a}&\simeq -0.4.\end{aligned}
$$
:::

## Mössbauer 効果

自由核では放出光子のエネルギーが反跳で減り、吸収には反跳分の余分なエネルギーが要る。核の遷移が狭い場合、この二つの線が離れて共鳴吸収しにくい。固体に束縛された核では、格子の振動状態を変えない成分があり、運動量移行を結晶全体で受け持てる。これが :term[Mössbauer 効果]{en="Mössbauer effect" yomi="めすばうあーこうか"} である。

光子の方向に沿う核の変位を $x$ とすると、巨視的な調和結晶の弾性線には長時間で残る変位相関が寄与する。調和格子の Gaussian 変位の分散を $\langle x^2\rangle$ とすれば、平衡平均は $\langle e^{ikx}\rangle=\exp[-k^2\langle x^2\rangle/2]$ である。連続的な格子振動スペクトルで時間に依存する相関が減衰すると、弾性線の無反跳分率はその平均の絶対値二乗
$$
f_{\rm LM}=\exp[-k^2\langle x^2\rangle]
$$
となる。線幅は平均寿命 $\tau$ の状態の振幅 $e^{-t/(2\tau)}$ を Fourier 変換して確率を取ると Lorentz 型で、半値全幅 $\Gamma=\hbar/\tau$ である。線源を速度 $u\ll c$ で動かすと一次の Doppler 変化は $\Delta E/E\simeq u/c$ であり、線の位置や超微細分裂を掃引して測れる。

:::example[自由核の反跳]
$E_\gamma=14.4\,\mathrm{keV}$、$M_fc^2\simeq57\cdot931.5\,\mathrm{MeV}$ とすると、自由核の反跳エネルギーは約 $1.95\,\mathrm{meV}$ である。放出と吸収の線の差はこれの約 2 倍になる。
:::

:::example[Doppler 掃引]
$14.4\,\mathrm{keV}$ の線を $1\,\mathrm{mm\,s^{-1}}$ で掃引すると移動の大きさは約 $4.8\times10^{-8}\,\mathrm{eV}$ である。単位を先に揃え、$u/c$ が約 $3.3\times10^{-12}$ であることを使う。
:::

### 格子の温度と線の幅

一方向の変位を質量 $M$、角周波数 $\omega_0$ の一振動子で表す。$x=\sqrt{\hbar/(2M\omega_0)}(a+a^\dagger)$ なので、占有数 $n$ が決まった状態では $\langle x^2\rangle=\hbar(2n+1)/(2M\omega_0)$ である。熱平衡では確率が $e^{-n\hbar\omega_0/(k_{\rm B}T)}$ に比例し、等比級数の微分で $\langle n\rangle=1/[e^{\hbar\omega_0/(k_{\rm B}T)}-1]$ となる。従って
$$
\langle x^2\rangle=\frac{\hbar}{2M\omega_0}
\coth\left(\frac{\hbar\omega_0}{2k_{\rm B}T}\right).
$$
低温では零点の寄与が残り、高温では $k_{\rm B}T/(M\omega_0^2)$ に近づく。実際の調和結晶では $x$ の各正常モードの変位分散を、対象核への寄与に応じて足すので、分散の温度依存は一周波数の場合と異なる。

等方的な単原子 Debye 模型では、変位の一方向に重み付けした振動数分布は $3\omega^2/\omega_D^3$（$0<\omega<\omega_D$）である。この規格化は三偏極の寄与を方向平均した結果である。各振動子の分散をこの分布で積分し、$\Theta_D=\hbar\omega_D/k_{\rm B}$、$E_R=(\hbar k)^2/(2M)$ と置く。$\coth(u/2)=1+2/(e^u-1)$ の零点項と熱励起項を分けると

$$
-\ln f_{\rm LM}=\frac{6E_R}{k_{\rm B}\Theta_D}
\left[\frac14+
\left(\frac{T}{\Theta_D}\right)^2
\int_0^{\Theta_D/T}\frac{u}{e^u-1}\dd u\right].
$$

零温でも指数の $3E_R/(2k_{\rm B}\Theta_D)$ が残る。高温では積分が $\Theta_D/T$ を先頭項とするため、指数は $6E_R T/(k_{\rm B}\Theta_D^2)$ に近づく。ここでの温度依存は三次元の巨視的な結晶の近似であり、離散的な少数の振動子で熱平均した同一準位への遷移確率には、そのまま適用できない。

:::example[自然幅と自由核の線の差]
平均寿命 $140\,\mathrm{ns}$ なら $\Gamma=\hbar/\tau\simeq4.7\times10^{-9}\,\mathrm{eV}$ である。自由核では先に求めた meV の反跳の方が何桁も大きいので、放出と吸収の線の中心は自然幅内には重ならない。格子を利用する意義はこの二つの尺度を比較すると分かる。
:::

### 自由反跳、無反跳成分、速度走査

自由核で光子エネルギーを $E$ とすると、反跳の運動エネルギーは $E^2/(2Mc^2)$ である。放出線は核準位差より下、吸収に必要なエネルギーは上にずれるので、自由な二核の線の不一致はその二倍になる。固体中で格子の内部状態を変えない成分では、運動量は試料全体が受け持つ。核の変位 $x$ の Gaussian 分布の特性関数は $\langle e^{ikx}\rangle=e^{-k^2\langle x^2\rangle/2}$ なので、強度はその絶対値二乗である。二つの固体を使う共鳴ではこの分率がそれぞれ掛かる。速度走査は $v\ll c$ での移動 $Ev/c$ を用いて微小な線の違いを調べる。

:::example[自由核での反跳]
光子のエネルギーと核の静止エネルギーをともに eV で用いる。

計算には次の値を用いる。
$$
E=14400\,\mathrm{eV},\qquad m_{c}=53000000000\,\mathrm{eV}
$$

この条件で求める量を ${E_R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{E_R}&=\frac{E^{2}}{2 \cdot m_{c}}\\&=\frac{\left(14400\right)^{2}}{2 \cdot \left(53000000000\right)}\\&\simeq 0.00195623\,\mathrm{eV}.\end{aligned}
$$

この量を核準位の自然幅と比較すると自由な共鳴が難しいか判定できる。

高エネルギーの線へ移す。次の値を変え、ほかは固定する。
$$
E=30000\,\mathrm{eV}
$$

このときの値を ${E_R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{E_R}_{\mathrm{new}}&=\frac{\left(30000\right)^{2}}{2 \cdot \left(53000000000\right)}\\&\simeq 0.00849057\,\mathrm{eV},\\{E_R}_{\mathrm{new}}-{E_R}&\simeq 0.00653434\,\mathrm{eV}.\end{aligned}
$$

**一次近似の検証。** 反跳の二次依存が、近傍の一次補正ではどの程度再現できるかを確認する。

${E}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${E_R}$ の一階微分を求めると、
$$
\frac{\mathrm d{E_R}}{\mathrm dE}=\frac{E}{m_{c}}.
$$

この導関数を変更前の点で評価して $h=15600$ を掛ける。ここで $E$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq 0.00423849\,\mathrm{eV},\\\Delta-\Delta_{\mathrm{lin}}&\simeq 0.00229585\,\mathrm{eV}.\end{aligned}
$$

残差が示すのは、指定した二点間で生じる線形化の誤差である。基準点からさらに離れた値について同じ精度を保証するものではない。
:::

:::example[無反跳強度の比較]
一次元の無次元分散 $u=k^2\langle x^2\rangle$ について $e^{-u}$ を求める。

計算には次の値を用いる。
$$
u=\frac{3}{10}
$$

この条件で求める量を ${f}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{f}&=e^{- u}\\&=e^{- \left(\frac{3}{10}\right)}\\&\simeq 0.740818.\end{aligned}
$$

同じ光子エネルギーなら格子の変位が大きいほど弾性成分が減る。

格子振動の大きい条件と比べる。次の値を変え、ほかは固定する。
$$
u=\frac{9}{10}
$$

このときの値を ${f}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{f}_{\mathrm{new}}&=e^{- \left(\frac{9}{10}\right)}\\&\simeq 0.40657,\\{f}_{\mathrm{new}}-{f}&\simeq -0.334249.\end{aligned}
$$
:::

:::example[発光体と吸収体の積]
発光側と吸収側の無反跳分率を $f_s,f_a$ とする。幾何学と吸収断面積を固定して無反跳事象の係数を得る。

計算には次の値を用いる。
$$
f_{s}=\frac{4}{5},\qquad f_{a}=\frac{3}{5}
$$

この条件で求める量を ${f}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{f}&=f_{a} \cdot f_{s}\\&=\left(\frac{3}{5}\right) \cdot \left(\frac{4}{5}\right)\\&\simeq 0.48.\end{aligned}
$$

共鳴線強度の絶対値には線形状なども必要である。

吸収体を変える。次の値を変え、ほかは固定する。
$$
f_{a}=\frac{9}{10}
$$

このときの値を ${f}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{f}_{\mathrm{new}}&=\left(\frac{4}{5}\right) \cdot \left(\frac{9}{10}\right)\\&\simeq 0.72,\\{f}_{\mathrm{new}}-{f}&\simeq 0.24.\end{aligned}
$$
:::

:::example[速度のエネルギー換算]
$E=14.4\,\mathrm{keV}$ の線を速度 $v$ で走査する。$v,c$ をともに m/s で与える。

計算には次の値を用いる。
$$
E=14400\,\mathrm{eV},\qquad v=\frac{1}{1000}\,\mathrm{m/s},\qquad c=299792458\,\mathrm{m/s}
$$

この条件で求める量を ${\Delta E}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\Delta E}&=\frac{E \cdot v}{c}\\&=\frac{\left(14400\right) \cdot \left(\frac{1}{1000}\right)}{\left(299792458\right)}\\&\simeq 4.80332\times10^{-8}\,\mathrm{eV}.\end{aligned}
$$

速度の符号に応じてエネルギー移動も反転する。

走査範囲を三倍にする。次の値を変え、ほかは固定する。
$$
v=\frac{3}{1000}\,\mathrm{m/s}
$$

このときの値を ${\Delta E}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\Delta E}_{\mathrm{new}}&=\frac{\left(14400\right) \cdot \left(\frac{3}{1000}\right)}{\left(299792458\right)}\\&\simeq 1.441\times10^{-7}\,\mathrm{eV},\\{\Delta E}_{\mathrm{new}}-{\Delta E}&\simeq 9.60665\times10^{-8}\,\mathrm{eV}.\end{aligned}
$$
:::

### 自然幅と速度で表した線幅

初期準位を平均寿命 $\tau$ で崩壊する状態とし、終状態の幅をゼロとする。初期状態の複素振幅は $e^{-i\omega_0t-t/(2\tau)}$ である。$t\ge0$ の Fourier 積分は $[1/(2\tau)-i(\omega-\omega_0)]^{-1}$ に比例するので、強度の半値点は $\omega-\omega_0=\pm1/(2\tau)$、全幅は $\Gamma=\hbar/\tau$ となる。発光と吸収の両方が Lorentz 形の場合は、発光エネルギーについての重なり積分で二つの幅が足される。これは規格化した Lorentz 形の Fourier 変換が $e^{-\Gamma |t|/(2\hbar)}$ であり、畳み込みでは積になるためである。狭い無反跳線を速度で測定すると、実験での幅にはこの自然幅の和に加えて、発光体・吸収体の厚みや場のばらつきも関係する。

:::example[速度の幅への変換]
初期の幅が $\Gamma_s=4.7\,\mathrm{neV}$、吸収側も $\Gamma_a=4.7\,\mathrm{neV}$ とする。$14.4\,\mathrm{keV}$ の線で二つの幅の和を速度差に換算すると、$c(\Gamma_s+\Gamma_a)/E\simeq0.196\,\mathrm{mm/s}$ になる。放出側だけの自然幅を速度に換算した値はその半分である。核種の幅を比較する場合、どちらの定義で表示したかを確認する。
:::

### 超微細相互作用による線の変化

核の電荷半径が励起前後で異なると、核内に入り込む電子の静電エネルギーの補正も異なる。その補正の差が線源と吸収体で異なると線の中心がずれる。これを :term[異性体シフト]{en="isomer shift" yomi="いせいたいしふと"} という。四重極モーメントを持つ準位では電場勾配が副準位のエネルギーを変え、内部磁場では核の磁気モーメントが副準位を分ける。

:::example[磁気分裂した線の本数]
量子化軸を静磁場にとり、スピン $J$ の磁気副準位に $-g_J\mu_NBM$ のエネルギーを加える。${}^{57}\mathrm{Fe}$ の $J_i=3/2\to J_f=1/2$ の線で双極子成分だけを考えると、各 $M_f=\pm1/2$ に対し $M_i=M_f-1,M_f,M_f+1$ の三つが許される。したがって一般の $g_i,g_f$ では六つの線の位置が生じる。各線の強度は観測方向や偏光で変わり、磁場に対する特定の方向では一部の線が消える。
:::

## 演習問題

::::exercise{level=basic}
$\alpha$ 崩壊で $M_d/m_\alpha=40$、$Q_\alpha=4.1\,\mathrm{MeV}$ のとき、娘核の運動エネルギーを求めよ。
:::solution
運動量の大きさが等しいのでエネルギー比は質量の逆比である。娘核には全量の $1/41$、すなわち $0.1\,\mathrm{MeV}$ が渡る。
:::
::::

::::exercise{level=basic}
$2^+\to0^+$ の単一光子遷移の多重極を求めよ。
:::solution
三角条件の上下限がともに 2 なので $L=2$ である。パリティ積は $+1$ であり、電気四重極 $E2$ が条件を満たす。
:::
::::

::::exercise{level=basic}
転換係数が 2 の一遷移で全率が $300\,\mathrm{s}^{-1}$ のとき、光子率を求めよ。
:::solution
全率は $(1+2)\lambda_\gamma$ なので $\lambda_\gamma=100\,\mathrm{s}^{-1}$ である。残りの $200\,\mathrm{s}^{-1}$ は転換電子に対応する。
:::
::::

::::exercise{level=basic}
$k^2\langle x^2\rangle=\ln2$ の格子で無反跳分率を求めよ。
:::solution
Gaussian の振幅の二乗が強度であり、$f_{\rm LM}=e^{-\ln2}=1/2$ となる。
:::
::::

::::exercise{level=advanced}
:eqref[eq-nuc-gamow] で $R/b\to0$ のとき、$2G$ の主項と次項を求めよ。
:::solution
$\arccos\sqrt x=\pi/2-\sqrt x+O(x^{3/2})$、$\sqrt{x(1-x)}=\sqrt x+O(x^{3/2})$ を引く。従って $2G=2\pi C/(\hbar v)-8C\sqrt{RQ/C}/(\hbar v)+O(Cx^{3/2}/(\hbar v))$ となる。$v=\sqrt{2Q/\mu}$ を代入すれば次項は $-4\sqrt{2\mu C R}/\hbar$ であり、$Q$ が消える。
:::
::::

::::exercise{level=advanced}
$E1$ の $B$ を固定し、光子エネルギーを 3 倍にした場合の半減期比を求めよ。他のチャネルはないとする。
:::solution
:eqref[eq-nuc-electric-rate] で $L=1$ なので率は $3^3=27$ 倍になる。半減期 $\ln2/\lambda$ は従って $1/27$ 倍である。
:::
::::

::::exercise{level=advanced}
$1^-\to1^+$ で許される単一光子の電気・磁気多重極を全て列挙せよ。
:::solution
角運動量で $0\le L\le2$、光子の条件で $L=1,2$ である。パリティは変わるので奇数 $L$ の電気と偶数 $L$ の磁気が許され、$E1,M2$ を得る。
:::
::::

::::exercise{level=advanced}
初期と終状態の質量差を $\Delta M=M_i-M_f>0$ とする。静止核の光子エネルギーが $\Delta M c^2$ より小さい分を厳密に求めよ。
:::solution
$M_f=M_i-\Delta M$ を $E_\gamma=(M_i^2-M_f^2)c^2/(2M_i)$ に代入する。$E_\gamma=\Delta M c^2-(\Delta M)^2c^2/(2M_i)$ なので、差は $(\Delta M)^2c^2/(2M_i)$ である。
:::
::::

## 確認テスト

::::quiz
$Q_\alpha>0$ の崩壊が長寿命になり得る理由を選べ。

- [ ] 正の $Q$ ではエネルギーが保存されない。
- [x] 生成物が外へ出る際の障壁透過が指数的に小さくなる。
- [ ] 反跳が運動エネルギーを必ず全て消す。

:::explanation
正の $Q$ なら二体の最終状態はエネルギー的に許される。核外の Coulomb ポテンシャルが相対運動エネルギーを超える区間では WKB 透過率が $e^{-2G}$ となり、十分小さくなり得る。
:::
::::

::::quiz{answer="32" tol="1e-6"}
$B(E2)$ を固定したとき $E_\gamma$ を 2 倍にすると率は何倍か。
:::explanation
:eqref[eq-nuc-electric-rate] のべきが $2L+1=5$ なので $2^5=32$ 倍である。
:::
::::

::::quiz
$0^+\to0^+$ の電磁遷移について正しいものを選べ。

- [ ] $E1$ 光子を一つ放出できる。
- [ ] $E0$ 光子を一つ放出できる。
- [x] 内部転換で進み得る。

:::explanation
二つの $J=0$ の間には角運動量階数 $0$ しか結合できない。実光子一個は $L\ge1$ なので放出できないが、単極子の電荷演算子は電子と結合し得る。
:::
::::

::::quiz{answer="0.25" tol="1e-6"}
一遷移の内部転換係数が 3 のとき、光子が出る分岐割合を求めよ。
:::explanation
電子率は光子率の 3 倍で、全率は 4 倍である。光子を選ぶ確率は $1/4=0.25$ である。
:::
::::

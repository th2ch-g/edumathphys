---
title: 核子多体系と平均場理論
description: Fermi 気体から核物質の飽和を考え、平均場と密度汎関数から状態方程式を求める。
notation:
  - symbol: '\rho,\rho_0'
    meaning: 核子数密度と対称核物質の飽和密度
  - symbol: 'k_{Fq}'
    meaning: 核子種 q の Fermi 波数
  - symbol: 'S(\rho)'
    meaning: 密度を固定した対称エネルギー
  - symbol: 'K_0'
    meaning: 対称核物質の非圧縮率
---

## この章の目標

- 核子を Fermi 球に詰めて運動エネルギーを求める。
- 核物質の飽和と圧力をエネルギー曲線から読む。
- Slater 行列式を変分して Hartree–Fock 方程式を導く。
- 密度汎関数と平均場から対称エネルギーを計算する。

## Fermi 気体模型

核の内部を粒子の出入りのない体積 $V=L^3$ の箱で近似する。境界を周期的に取り、温度をゼロとする。波数は $2\pi\vb n/L$、$\vb n\in\mathbb Z^3$ であり、一状態の波数体積は $(2\pi)^3/V$ である。各核子種のスピン縮退度は 2 である。粒子を低い $\varepsilon_k=\hbar^2 k^2/(2m)$ から占有させると、核子種 $q=n,p$ の密度と運動エネルギー密度は
$$
\begin{aligned}
\rho_q &=\frac{2}{(2\pi)^3}4\pi\int_0^{k_{Fq}} k^2\dd k=\frac{k_{Fq}^3}{3\pi^2},\\
\mathcal T_q &=\frac{2}{(2\pi)^3}4\pi\int_0^{k_{Fq}}\frac{\hbar^2k^2}{2m}k^2\dd k
=\frac35\rho_q \frac{\hbar^2k_{Fq}^2}{2m}.
\end{aligned}
\label{eq-nuc-fermi-gas}
$$
球内で一様に数えるため、運動エネルギーの平均は境界の値の $3/5$ になる。対称核物質では $\rho_n=\rho_p=\rho/2$、$\rho=\rho_n+\rho_p$ である。Coulomb 相互作用を除き、陽子と中性子を等質量 $m$ とする。

:::example[飽和密度における運動]
$\rho=0.16\,\mathrm{fm}^{-3}$ で $k_F=(3\pi^2\rho/2)^{1/3}\simeq1.33\,\mathrm{fm}^{-1}$ となる。$mc^2\simeq939\,\mathrm{MeV}$ なら境界のエネルギーは約 $37\,\mathrm{MeV}$、平均は約 $22\,\mathrm{MeV}$ である。$\hbar k_F/(mc)\simeq0.28$ は非相対論近似の目安になる。
:::

:::example[純中性子物質]
同じ $\rho$ を中性子だけで満たすと $k_{Fn}$ は対称物質の $2^{1/3}$ 倍である。運動エネルギーの平均は $2^{2/3}$ 倍になる。陽子へ一部を振り分けると Pauli 原理による高運動量の占有を減らせる。
:::

### 有限核の密度に戻す

一様な対称核物質の飽和密度 $\rho_0$ を有限核の内部で使うなら、$A=(4\pi/3)R^3\rho_0$ から $r_0=(3/(4\pi\rho_0))^{1/3}$ を得る。$\rho_0=0.16\,\mathrm{fm}^{-3}$ では $r_0\simeq1.14\,\mathrm{fm}$ である。有限核では表面近くで密度が下がるため、一様球の等価半径と測定した半径を同じ定義で比較する必要がある。体積に比例するエネルギーも表面、非対称、Coulomb によって変わる。

核内の深い位置で平均場 $U$ が一定とし、真空中の一核子のエネルギーをゼロとして $\epsilon_k=\hbar^2k^2/(2m)+U$ と置く。軌道と場を凍結して最も弱く束縛された中性子を取り去る模型では、最上部の軌道のエネルギー $E_F+U$ の負が分離エネルギーを与える。例えば $E_F=37\,\mathrm{MeV}$、$S_n=8\,\mathrm{MeV}$ なら $U\simeq-45\,\mathrm{MeV}$ と見積もれる。この値と一核子当たり全結合の差には、占有軌道の運動エネルギーの平均に加え、相互作用を固有値和で二重に数える効果がある。

### 状態密度、化学ポテンシャル、有限温度

$k$ 空間の一状態が張る体積は $(2\pi)^3/V$ である。スピン・アイソスピンの縮退度を $g$ とすれば、運動量の方向を積分した状態密度は $Vg k^2\dd k/(2\pi^2)$ である。$\epsilon=\hbar^2k^2/(2m)$ と置いて $D(\epsilon)=\dd N/\dd\epsilon$ を得る。$\dd N/\dd k$ に $m/(\hbar^2 k)$ を掛けると $D(\epsilon)\propto\sqrt\epsilon$、また $D(E_F)=3N/(2E_F)$ である。零温度の化学ポテンシャルは、粒子を一個増やす際の限界エネルギーである。$E_{\rm kin}\propto N^{5/3}V^{-2/3}$ から $(\partial_N E_{\rm kin})_V=(5/3) E_{\rm kin}/N=E_F$ となる。

有限温度で $f(\epsilon)=[1+\exp((\epsilon-\mu)/k_{\rm B}T)]^{-1}$ とする。零温度との差を $u=(\epsilon-\mu)/(k_{\rm B}T)$ で積分すると、分布の差は $u$ に関して奇関数である。滑らかな $D$ を $\mu$ の周りで展開すると、最初の項は $2(k_{\rm B}T)^2D'(\mu)\int_0^\infty u/(1+e^u)\dd u$ である。$k_{\rm B}T\ll E_F$ のとき下限の延長の誤差は指数的に小さい。積分の値は幾何級数を積分した $\sum_{j\ge1}(-1)^{j-1}/j^2=\pi^2/12$ である。よって粒子数一定では $\mu=E_F-\pi^2(k_{\rm B}T)^2/(12E_F)+O(T^4)$ となる。

:::example[縮退度による波数の違い]
同じ総密度で縮退度 $g$ を変える。状態数を球内で積分する。

量の定義に従い、共通の単位で次の値を用いる。
$$
\rho=\frac{4}{25}\,\mathrm{fm}^{-3},\qquad g=4
$$

この条件で求める量を ${k_F}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{k_F}&=\left(\frac{6 \cdot \pi^{2} \cdot \rho}{g}\right)^{1 \cdot \frac{1}{3}}\\&=\left(\frac{6 \cdot \pi^{2} \cdot \left(\frac{4}{25}\right)}{\left(4\right)}\right)^{1 \cdot \frac{1}{3}}\\&\simeq 1.33302\,\mathrm{fm}^{-1}.\end{aligned}
$$

対称系から同種粒子だけの系に移すと、同じ密度で Fermi 球が大きくなる。

スピンの縮退だけを残す。そこで次の値に変更し、ほかの入力を固定する。
$$
g=2
$$

このときの値を ${k_F}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{k_F}_{\mathrm{new}}&=\left(\frac{6 \cdot \pi^{2} \cdot \left(\frac{4}{25}\right)}{\left(2\right)}\right)^{1 \cdot \frac{1}{3}}\\&\simeq 1.6795\,\mathrm{fm}^{-1},\\{k_F}_{\mathrm{new}}-{k_F}&\simeq 0.34648\,\mathrm{fm}^{-1}.\end{aligned}
$$

**一次近似の検証。** 縮退度は離散的だが、三乗根を実数に延ばした感度と、厳密に縮退度を変更した値を比較できる。

${g}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${k_F}$ の一階微分を求めると、
$$
\frac{\mathrm d{k_F}}{\mathrm dg}=- \frac{\sqrt[3]{6} \cdot \pi^{\frac{2}{3}} \cdot \sqrt[3]{\frac{\rho}{g}}}{3 \cdot g}.
$$

この導関数を変更前の点で評価して $h=-2$ を掛ける。ここで $g$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq 0.22217\,\mathrm{fm}^{-1},\\\Delta-\Delta_{\mathrm{lin}}&\simeq 0.12431\,\mathrm{fm}^{-1}.\end{aligned}
$$

残差が示すのは、指定した二点間で生じる線形化の誤差である。基準点からさらに離れた値について同じ精度を保証するものではない。
:::

:::example[Fermi 面の状態密度]
一体準位の連続近似で $D(E_F)$ を粒子数 $N$ と $E_F$ で求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
N=100,\qquad E_{F}=36\,\mathrm{MeV}
$$

この条件で求める量を ${D_F}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{D_F}&=\frac{3 \cdot N}{2 \cdot E_{F}}\\&=\frac{3 \cdot \left(100\right)}{2 \cdot \left(36\right)}\\&\simeq 4.16667\,\mathrm{MeV}^{-1}.\end{aligned}
$$

殻の細かい構造を平均した値である。小さな有限核の準位をこの密度で離散的に予測することはできない。

より高い密度の系で比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
E_{F}=48\,\mathrm{MeV}
$$

このときの値を ${D_F}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{D_F}_{\mathrm{new}}&=\frac{3 \cdot \left(100\right)}{2 \cdot \left(48\right)}\\&\simeq 3.125\,\mathrm{MeV}^{-1},\\{D_F}_{\mathrm{new}}-{D_F}&\simeq -1.04167\,\mathrm{MeV}^{-1}.\end{aligned}
$$
:::

:::example[低温での化学ポテンシャル]
温度のエネルギー換算を $\vartheta$ とする。二次の低温補正を使う。

量の定義に従い、共通の単位で次の値を用いる。
$$
E_{F}=36\,\mathrm{MeV},\qquad \theta=1\,\mathrm{MeV}
$$

この条件で求める量を ${\mu}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\mu}&=E_{F} - \frac{\pi^{2} \cdot \theta^{2}}{12 \cdot E_{F}}\\&=- \frac{\pi^{2} \cdot \left(1\right)^{2}}{12 \cdot \left(36\right)} + \left(36\right)\\&\simeq 35.9772\,\mathrm{MeV}.\end{aligned}
$$

粒子数を固定すると Fermi 面近くの温度による変化は化学ポテンシャルの低下を伴う。

低温域内で温度を上げる。そこで次の値に変更し、ほかの入力を固定する。
$$
\theta=4\,\mathrm{MeV}
$$

このときの値を ${\mu}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\mu}_{\mathrm{new}}&=\left(36\right) - \frac{\pi^{2} \cdot \left(4\right)^{2}}{12 \cdot \left(36\right)}\\&\simeq 35.6345\,\mathrm{MeV},\\{\mu}_{\mathrm{new}}-{\mu}&\simeq -0.342695\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[非相対論の誤差]
運動量のエネルギー換算を $p_c$、静止エネルギーを $m_c$ とする。平方根の展開から四次補正は負である。

量の定義に従い、共通の単位で次の値を用いる。
$$
p_{c}=200\,\mathrm{MeV},\qquad m_{c}=939\,\mathrm{MeV}
$$

この条件で求める量を ${\Delta E}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\Delta E}&=\frac{\left(-1\right) \cdot p_{c}^{4}}{8 \cdot m_{c}^{3}}\\&=- \frac{\left(200\right)^{4}}{8 \cdot \left(939\right)^{3}}\\&\simeq -0.241565\,\mathrm{MeV}.\end{aligned}
$$

この項と二次の運動エネルギーを比べれば非相対論近似の誤差を見積もれる。

運動量を上げて比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
p_{c}=260\,\mathrm{MeV}
$$

このときの値を ${\Delta E}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\Delta E}_{\mathrm{new}}&=- \frac{\left(260\right)^{4}}{8 \cdot \left(939\right)^{3}}\\&\simeq -0.689933\,\mathrm{MeV},\\{\Delta E}_{\mathrm{new}}-{\Delta E}&\simeq -0.448368\,\mathrm{MeV}.\end{aligned}
$$
:::

## 核物質の飽和性

表面効果のない無限核物質の一核子当たり内部エネルギーを $e(\rho)$ とする。静止エネルギーは基準から引く。粒子数 $A$ を固定して $E=Ae(A/V)$ を体積で微分すると
$$
P=-\left.\frac{\partial E}{\partial V}\right|_A
=\rho^2 e'(\rho).
\label{eq-nuc-matter-pressure}
$$
対称核物質では $e$ が $\rho_0\simeq0.16\,\mathrm{fm}^{-3}$ 付近でおよそ $-16\,\mathrm{MeV}$ の極小を取る。この性質を :term[飽和]{en="saturation" yomi="ほうわ"} という。$P(\rho_0)=0$ で、圧縮すると密度が増し正の圧力が生じる。無限の荷電物質には通常の Coulomb 項をそのまま入れない。

曲率を $K_0=9\rho_0^2 e''(\rho_0)$ と定義する。$x=(\rho-\rho_0)/\rho_0$ に対して $e=e_0+K_0x^2/18+O(x^3)$ である。圧力の一次応答は $P\simeq(K_0/9)(\rho-\rho_0)$ である。自由 Fermi 気体の $e\propto\rho^{2/3}$ は単調であるため、飽和の説明には相互作用が必要である。中程度の密度で引力が勝ち、高い密度で運動エネルギーや斥力の寄与が増す。

:::example[圧縮のコスト]
$K_0=240\,\mathrm{MeV}$ の二次近似で $x=0.1$ とすると、粒子当たりの増分は $240(0.1)^2/18\simeq0.133\,\mathrm{MeV}$ である。同じ変位に対する圧力の一次近似は $0.427\,\mathrm{MeV\,fm}^{-3}$ である。
:::

:::example[密度依存模型]
$x=\rho/\rho_0$ として $e=a x^{2/3}+b x+c x^2$ と置く。$e(1)=-16$、$\partial_x e(1)=0$ を要求すると $c=16+a/3$、$b=-32-4a/3$ である。$a>0$ なら $K_0=9(-2a/9+2c)=288+4a$ となる。飽和点を合わせても曲率は独立に合わせられない例である。
:::

### 飽和を満たす局所的な模型

$u=\rho/\rho_0$ とし、$e=a u^{2/3}+b u+c u^2$ とする。この形は Fermi 運動項、密度に比例する引力、より高次の斥力を含む。$e(1)=e_0$ と $\partial_u e(1)=0$ の二式を解けば $c=a/3-e_0$、$b=2e_0-4a/3$ である。曲率は $K_0=9\partial_u^2e(1)=18c-2a$ なので、この簡単な模型では $e_0,a$ を指定すると圧縮率を独立に調節できない。多くの係数をもつ汎関数が必要な理由を、飽和の位置と曲率を区別することで理解できる。

:::example[模型の引力係数]
飽和条件の解から $b$ を求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
e_{0}=-16\,\mathrm{MeV},\qquad a=\frac{45}{2}\,\mathrm{MeV}
$$

この条件で求める量を ${b}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{b}&=- \frac{4 \cdot a}{3} + 2 \cdot e_{0}\\&=2 \cdot \left(-16\right) - \frac{4 \cdot \left(\frac{45}{2}\right)}{3}\\&=-62\,\mathrm{MeV}.\end{aligned}
$$

束縛は $e_0<0$ なので、運動エネルギーと反対符号の係数が必要になる。

弱い飽和結合を仮定する。そこで次の値に変更し、ほかの入力を固定する。
$$
e_{0}=-15\,\mathrm{MeV}
$$

このときの値を ${b}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{b}_{\mathrm{new}}&=2 \cdot \left(-15\right) - \frac{4 \cdot \left(\frac{45}{2}\right)}{3}\\&=-60\,\mathrm{MeV},\\{b}_{\mathrm{new}}-{b}&=2\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[模型の曲率の過大評価]
同じ解を用い $c$ を消去する。$K_0=4a-18e_0$ である。

量の定義に従い、共通の単位で次の値を用いる。
$$
a=\frac{45}{2}\,\mathrm{MeV},\qquad e_{0}=-16\,\mathrm{MeV}
$$

この条件で求める量を ${K_0}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{K_0}&=4 \cdot a - 18 \cdot e_{0}\\&=- 18 \cdot \left(-16\right) + 4 \cdot \left(\frac{45}{2}\right)\\&=378\,\mathrm{MeV}.\end{aligned}
$$

この曲率を経験的に推定される尺度と比べれば、飽和点を合わせるだけでは不十分と分かる。

飽和エネルギーの差に対する敏感さを確認する。そこで次の値に変更し、ほかの入力を固定する。
$$
e_{0}=-15\,\mathrm{MeV}
$$

このときの値を ${K_0}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{K_0}_{\mathrm{new}}&=- 18 \cdot \left(-15\right) + 4 \cdot \left(\frac{45}{2}\right)\\&=360\,\mathrm{MeV},\\{K_0}_{\mathrm{new}}-{K_0}&=-18\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[圧縮時の一次の圧力]
飽和近くで $P\simeq K_0(\rho-\rho_0)/9$ とする。

量の定義に従い、共通の単位で次の値を用いる。
$$
K_{0}=240\,\mathrm{MeV},\qquad d=\frac{1}{125}\,\mathrm{fm}^{-3}
$$

この条件で求める量を ${P}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{P}&=\frac{K_{0} \cdot d}{9}\\&=\frac{\left(240\right) \cdot \left(\frac{1}{125}\right)}{9}\\&\simeq 0.213333\,\mathrm{MeV\,fm}^{-3}.\end{aligned}
$$

$\rho$ が飽和より小さい側で負の圧力になる。これはモデルの均一な枝についての結果である。

同じ幅だけ希薄な側と比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
d=- \frac{1}{125}\,\mathrm{fm}^{-3}
$$

このときの値を ${P}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{P}_{\mathrm{new}}&=\frac{\left(- \frac{1}{125}\right) \cdot \left(240\right)}{9}\\&\simeq -0.213333\,\mathrm{MeV\,fm}^{-3},\\{P}_{\mathrm{new}}-{P}&\simeq -0.426667\,\mathrm{MeV\,fm}^{-3}.\end{aligned}
$$
:::

:::example[二次の圧縮エネルギー]
粒子当たりの飽和からの増分を求める。$u$ は密度比である。

量の定義に従い、共通の単位で次の値を用いる。
$$
K_{0}=240\,\mathrm{MeV},\qquad u=\frac{21}{20}
$$

この条件で求める量を ${\Delta e}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\Delta e}&=\frac{K_{0} \cdot \left(u - 1\right)^{2}}{18}\\&=\frac{\left(240\right) \cdot \left(\left(\frac{21}{20}\right) - 1\right)^{2}}{18}\\&\simeq 0.0333333\,\mathrm{MeV}.\end{aligned}
$$

ずれを二倍にすると二次の増分は四倍になる。

より大きい圧縮で比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
u=\frac{11}{10}
$$

このときの値を ${\Delta e}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\Delta e}_{\mathrm{new}}&=\frac{\left(240\right) \cdot \left(\left(\frac{11}{10}\right) - 1\right)^{2}}{18}\\&\simeq 0.133333\,\mathrm{MeV},\\{\Delta e}_{\mathrm{new}}-{\Delta e}&\simeq 0.1\,\mathrm{MeV}.\end{aligned}
$$

**一次近似の検証。** 飽和点を離れた状態の弾性応答を線形化し、二次式との違いを確認する。

${u}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${\Delta e}$ の一階微分を求めると、
$$
\frac{\mathrm d{\Delta e}}{\mathrm du}=\frac{K_{0} \cdot \left(u - 1\right)}{9}.
$$

この導関数を変更前の点で評価して $h=\frac{1}{20}$ を掛ける。ここで $u$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq 0.0666667\,\mathrm{MeV},\\\Delta-\Delta_{\mathrm{lin}}&\simeq 0.0333333\,\mathrm{MeV}.\end{aligned}
$$

残差が示すのは、指定した二点間で生じる線形化の誤差である。基準点からさらに離れた値について同じ精度を保証するものではない。
:::

### 圧力の微分と音速

一様な冷たい物質の全エネルギー密度を $\epsilon=\rho[mc^2+e(\rho)]$ とする。組成一定の断熱小振動での音速は $c_s^2=c^2\dd P/\dd\epsilon$ である。圧力の式 :eqref[eq-nuc-matter-pressure] を微分すると $\dd P/\dd\rho=2\rho e'+\rho^2e''$ になる。$\rho_0$ で $e'=0$ を用いれば $\dd P/\dd\rho=K_0/9$、また $\dd\epsilon/\dd\rho=mc^2+e(\rho_0)$ である。従って
$$
c_s^2/c^2=K_0/[9(mc^2+e(\rho_0))]
$$
となる。$\rho_0$ で圧力がゼロでも、応答を決めるのは圧力の値でなく密度微分である。極小近くでは $K_0>0$ が音の実周波数を支える。高密度側へこの二次展開をそのまま延ばすと、相対論的に許される音速の上限を超える場合があり、飽和近くでの当てはめの外を独立に調べる必要がある。

## Hartree–Fock 理論

### 変分する集合

[第二量子化](/physics/many-body/second-quantization-review/) の記法を用いる。正規直交軌道 $\phi_i$ を占有した Slater 行列式について、$t$ を一体演算子、$v$ を対称な二体相互作用とする。二粒子の密度は占有空間への射影の反対称化なので、期待値は
$$
E=\sum_i\bra i t\ket i+\frac12\sum_{ij}
  \left(\bra{ij}v\ket{ij}-\bra{ij}v\ket{ji}\right).
$$
交換項は行列式で粒子を置換したときの負号から生じる。$i=j$ は直接項と交換項が相殺する。軌道を独立に変分し、拘束 $\braket{\phi_j}{\phi_i}=\delta_{ji}$ に行列乗数を付ける。相互作用の二つの占有添字から同じ変分が現れるため係数 $1/2$ が消え、
$$
\begin{aligned}
h\ket{\phi_i}&=\sum_j\ket{\phi_j}\Lambda_{ji},\\
h_{ab}&=t_{ab}+\sum_{j\ {\rm occ}}
(\bra{aj}v\ket{bj}-\bra{aj}v\ket{jb})
\end{aligned}
\label{eq-nuc-hartree-fock}
$$
を得る。$\Lambda$ は Hermite 行列であり、占有軌道内のユニタリ変換で対角化できる。これが :term[Hartree–Fock 方程式]{en="Hartree–Fock equation" yomi="はーとりーふぉっくほうていしき"} である。$h$ は選んだ占有軌道自身に依存するので、方程式を反復して密度と整合させる。

### 一体密度を使う変分

$\rho_{ab}=\sum_i\braket{a}{\phi_i}\braket{\phi_i}{b}$ とすると $\rho^2=\rho$、$\Tr\rho=A$ である。$U=e^X$、$X^\dagger=-X$ による軌道回転で $\delta\rho=[X,\rho]$ となり、$\delta E=\Tr(h[X,\rho])=\Tr([\rho,h]X)$ である。任意の $X$ に対して停留する条件は $[\rho,h]=0$ となる。停留条件を満たしていても最小とは限らないため、異なる初期占有を用いた比較や安定性の検査が必要になる。

:::example[平均場固有値の和]
$\varepsilon_i=\bra i h\ket i$ を足すと、二体エネルギーを 2 回数える。運動エネルギー合計 10、二体エネルギー $-6\,\mathrm{MeV}$ なら固有値和は $-2\,\mathrm{MeV}$、全エネルギーは $4\,\mathrm{MeV}$ である。
:::

:::example[凍結軌道での除去]
占有軌道 $k$ だけを除き他を固定すると、残留系と元の系の全エネルギー差は $-\varepsilon_k$ になる。実際、二体和から添字 $k$ を含む二列が消え、自己項はゼロである。軌道を再最適化した除去エネルギーには緩和が加わる。
:::

### 反対称二体密度の成分

占有軌道の作成演算子で作った状態 $\ket\Phi$ に対して
$\langle a_a^\dagger a_b^\dagger a_d a_c\rangle=\rho_{ca}\rho_{db}-\rho_{da}\rho_{cb}$
である。確認するには $\rho$ が対角の軌道基底に移る。消滅演算子二つが取り除いた軌道を作成演算子が戻す場合だけ期待値が残る。直接に戻す順序では $+1$、交差する順序では交換の負号が付き、両軌道が占有でなければゼロである。同じ添字は反交換でゼロである。基底変換の係数を各添字に掛ければ一般の式に戻る。この式を $H=\sum_{ab}t_{ab}a_a^\dagger a_b+\tfrac12\sum_{abcd}v_{ab,cd}a_a^\dagger a_b^\dagger a_d a_c$ に代入すると、エネルギーの直接項と交換項が得られる。

### 自己無撞着解の安定性

Hartree–Fock の変分原理は同じ粒子数の Slater 行列式の範囲で働く。$\phi_i\to\cos\theta\,\phi_i+\sin\theta\,\phi_a$ と回せば直交制約を保てる。相互作用がなければ $E(\theta)-E(0)=(t_{aa}-t_{ii})\sin^2\theta+2\Re t_{ai}\sin\theta\cos\theta$ である。$t_{ai}=0$ により一次の変化が消えても、$t_{aa}<t_{ii}$ であれば二次の変化が負になる。相互作用がある場合には $h$ が回転とともに変わるため、二次変分に密度応答の項も入る。粒子穴方向のこの二次変分が [集団運動](/physics/nuclear-physics/collective-model/) の小振動に現れる。

:::example[停留点と最小点]
一粒子 Hamiltonian が $\operatorname{diag}(1,3)\,\mathrm{MeV}$ で上の軌道だけを占有しても交換子 $[h,\rho]$ はゼロである。しかし下の軌道を振幅 $\sin\theta$ で混ぜると、エネルギー差は $-2\sin^2\theta\,\mathrm{MeV}$ である。従ってこの停留状態は変分の最小値ではない。
:::

### 占有軌道の回転と最低解の判定

占有射影子を $\rho=\sum_{i\in{\rm occ}}\ket i\bra i$ とする。任意の微小な反 Hermite 行列 $X$ で指数共役すると、一次で $\delta\rho=[X,\rho]$ となり、射影性と階数が保たれる。$\delta E=\operatorname{Tr}(h\delta\rho)$ の循環性から $\delta E=\operatorname{Tr}([\rho,h]X)$ である。任意の $X$ に対して停留するには $[\rho,h]=0$ が必要十分である。実際、零でなければ $X=[\rho,h]^\dagger$ を選んだ変分は $\operatorname{Tr}(AA^\dagger)>0$ となる。この条件で占有と非占有の間の行列要素は消える。一方で曲率が負なら停留点は最小ではない。

:::example[一つの占有軌道の回転]
二準位の一体ハミルトニアンが対角で、$\phi=\cos\theta\ket1+\sin\theta\ket2$ とする。第一準位を基準にした変分エネルギーを計算する。

量の定義に従い、共通の単位で次の値を用いる。
$$
d=3\,\mathrm{MeV},\qquad \theta=\frac{1}{10}\,\mathrm{rad}
$$

この条件で求める量を ${\Delta E}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\Delta E}&=d \cdot \sin^{2}{\left(\theta \right)}\\&=\left(3\right) \cdot \sin^{2}{\left(\left(\frac{1}{10}\right) \right)}\\&\simeq 0.0299001\,\mathrm{MeV}.\end{aligned}
$$

小角では上昇は二次である。ギャップ $d$ が負なら同じ停留配置は不安定になる。

軌道回転の角度を広げる。そこで次の値に変更し、ほかの入力を固定する。
$$
\theta=\frac{1}{5}\,\mathrm{rad}
$$

このときの値を ${\Delta E}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\Delta E}_{\mathrm{new}}&=\left(3\right) \cdot \sin^{2}{\left(\left(\frac{1}{5}\right) \right)}\\&\simeq 0.118409\,\mathrm{MeV},\\{\Delta E}_{\mathrm{new}}-{\Delta E}&\simeq 0.0885084\,\mathrm{MeV}.\end{aligned}
$$

**一次近似の検証。** 最小点近くの軌道の許容誤差がエネルギーにどう効くかを、有限な角度変化で検討する。

${\theta}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${\Delta E}$ の一階微分を求めると、
$$
\frac{\mathrm d{\Delta E}}{\mathrm d\theta}=d \cdot \sin{\left(2 \cdot \theta \right)}.
$$

この導関数を変更前の点で評価して $h=\frac{1}{10}$ を掛ける。ここで $\theta$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq 0.0596008\,\mathrm{MeV},\\\Delta-\Delta_{\mathrm{lin}}&\simeq 0.0289076\,\mathrm{MeV}.\end{aligned}
$$

残差が示すのは、指定した二点間で生じる線形化の誤差である。基準点からさらに離れた値について同じ精度を保証するものではない。
:::

:::example[二体エネルギーの二重計数]
対角化した一体の固有値和を $S$、一体運動の期待値和を $T$ とする。密度非依存の二体力では $E=(S+T)/2$ である。

量の定義に従い、共通の単位で次の値を用いる。
$$
S=-300\,\mathrm{MeV},\qquad T=140\,\mathrm{MeV}
$$

この条件で求める量を ${E}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{E}&=\frac{S + T}{2}\\&=\frac{\left(-300\right)}{2} + \frac{\left(140\right)}{2}\\&=-80\,\mathrm{MeV}.\end{aligned}
$$

$S$ は各対を両粒子について一回ずつ数える。多体の全エネルギーは和そのものにならない。

運動項を変えた解と比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
T=160\,\mathrm{MeV}
$$

このときの値を ${E}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{E}_{\mathrm{new}}&=\frac{\left(-300\right)}{2} + \frac{\left(160\right)}{2}\\&=-70\,\mathrm{MeV},\\{E}_{\mathrm{new}}-{E}&=10\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[直接項と交換項]
同種粒子の二つの占有軌道で、ある対の直接行列要素を $D$、交換行列要素を $X$ とする。

量の定義に従い、共通の単位で次の値を用いる。
$$
D=-3\,\mathrm{MeV},\qquad X=-1\,\mathrm{MeV}
$$

この条件で求める量を ${V_{\rm AS}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{V_{\rm AS}}&=D - X\\&=- \left(-1\right) + \left(-3\right)\\&=-2\,\mathrm{MeV}.\end{aligned}
$$

交換の符号は反対称化に由来する。引力の行列要素同士なら差が引力を弱める場合がある。

交換の大きい軌道を比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
X=-2\,\mathrm{MeV}
$$

このときの値を ${V_{\rm AS}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{V_{\rm AS}}_{\mathrm{new}}&=- \left(-2\right) + \left(-3\right)\\&=-1\,\mathrm{MeV},\\{V_{\rm AS}}_{\mathrm{new}}-{V_{\rm AS}}&=1\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[粒子穴の一次応答]
非相互作用の小振幅極限で、ギャップ $d$ の軌道に摂動 $f$ が結合する。下側の射影から上側への振幅の大きさは $f/d$ である。

量の定義に従い、共通の単位で次の値を用いる。
$$
f=\frac{1}{5}\,\mathrm{MeV},\qquad d=4\,\mathrm{MeV}
$$

この条件で求める量を ${P}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{P}&=\left(\frac{f}{d}\right)^{2}\\&=\frac{\left(\frac{1}{5}\right)^{2}}{\left(4\right)^{2}}\\&\simeq 0.0025.\end{aligned}
$$

自己無撞着な応答には密度変化で誘起された場も必要である。

近い空軌道の場合と比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
d=2\,\mathrm{MeV}
$$

このときの値を ${P}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{P}_{\mathrm{new}}&=\frac{\left(\frac{1}{5}\right)^{2}}{\left(2\right)^{2}}\\&\simeq 0.01,\\{P}_{\mathrm{new}}-{P}&\simeq 0.0075.\end{aligned}
$$
:::

## Skyrme 力と核密度汎関数

有限範囲の力を低運動量で展開すると、接触相互作用と座標の勾配を含む項が得られる。Skyrme 型模型はこの構成を局所密度、運動量分布の二次モーメント、スピン軌道密度などで表す。密度依存の項は多体効果を取り込むためにも使う。密度汎関数では、変分の対象を占有軌道から構成した各密度としてエネルギーを指定する。

### 接触型の擬ポテンシャル

$d=\delta^3(\vb r_1-\vb r_2)$、$\vb R=(\vb r_1+\vb r_2)/2$ とする。相対波数演算子 $\vb k=(\grad_1-\grad_2)/(2i)$ は右の波動関数に、$\vb k'$ はその Hermite 共役として左に作用する。$P_\sigma=(1+\vb\sigma_1\cdot\vb\sigma_2)/2$ はスピン交換である。この記法で中心・スピン軌道の代表的な :term[Skyrme 力]{en="Skyrme interaction" yomi="すきるむりょく"} は
$$
\begin{aligned}
v_{\rm Sk}={}&t_0(1+x_0P_\sigma)d
+\frac{t_1}{2}(1+x_1P_\sigma)(\vb k'^2d+d\vb k^2)\\
&+t_2(1+x_2P_\sigma)\vb k'\cdot d\vb k
+\frac{t_3}{6}(1+x_3P_\sigma)\rho(\vb R)^\alpha d\\
&+iW_0(\vb\sigma_1+\vb\sigma_2)\cdot(\vb k'\times d\vb k).
\end{aligned}
$$
$t_i,x_i,W_0,\alpha$ は模型の係数である。密度依存項を含めると、変分では $\rho^\alpha$ の微分がさらに必要になる。実際の密度汎関数ではこの擬ポテンシャルに由来する係数の関係を緩めることもある。

均一でスピンを偏極させない対称物質では、等確率の二スピン・二アイソスピンを平均する。空間交換は散乱終状態の相対波数の符号を反転し、$t_0,t_1,t_3$ は偶、$t_2$ は奇である。従って内部自由度の交換部分 $P_\sigma P_\tau$ を用い、スピン・アイソスピン平均に $\langle P_\sigma\rangle=\langle P_\tau\rangle=1/2$、$\langle P_\sigma P_\tau\rangle=1/4$ を入れると、偶項の重みは $3/4$、奇項の重みは $5/4+x_2$ である。平面波で $\langle|\vb k_1-\vb k_2|^2/4\rangle=3k_F^2/10$ となり、スピン軌道項は外積が消える。各対を一回だけ数えると一核子当たりのエネルギーは
$$
e(\rho,0)=\frac35 E_F+\frac38t_0\rho
+\frac1{16}t_3\rho^{\alpha+1}
+\frac3{80}\big[3t_1+(5+4x_2)t_2\big]k_F^2\rho.
$$
有限核でも同じ反対称化を行うが、波動関数の勾配から表面とスピン軌道の項が残る。以下ではその変分が見えるよう、局所汎関数の一部を取り出す。

### 局所密度で書いたエネルギー

中心で時間反転対称な場合の一部を抜き出す。$\rho=\sum_i|\phi_i|^2$、$\tau=\sum_i|\grad\phi_i|^2$ とし、
$$
E=\int\left[\frac{\hbar^2}{2m}\tau+
C_0\rho^2+C_3\rho^{\alpha+2}+C_\tau\rho\tau+C_\nabla|\grad\rho|^2\right]\dd^3r .
$$
$\phi_i^*$ で変分すると $\tau$ の変分は $\sum_i\grad(\delta\phi_i^*)\cdot\grad\phi_i$ である。境界項が消える条件で部分積分を行う。さらに $\delta|\grad\rho|^2=2\grad\rho\cdot\grad\delta\rho$ を積分すると
$$
h=-\grad\cdot\left(\frac{\hbar^2}{2m}+C_\tau\rho\right)\grad+
2C_0\rho+(\alpha+2)C_3\rho^{\alpha+1}+C_\tau\tau-2C_\nabla\nabla^2\rho
$$
となる。$\rho^{\alpha+2}$ の微分により、密度依存結合を固定した場合には得られない項が現れる。これを再配列の寄与という。Skyrme 汎関数の全項と対称性の拘束は、時間奇の密度や回転を含むため :cite[ring-schuck] を参照する。

:::example[実効質量]
一様密度で $C_\tau\rho=\hbar^2/(4m)$ とする。勾配の係数を $\hbar^2/(2m^*)$ と書くと $m^*=2m/3$ となる。一定の $k$ に対する運動エネルギー分散の係数は $3/2$ 倍になる。
:::

:::example[表面の勾配]
$\rho(z)=\bar\rho+\epsilon\cos kz$ とすると、周期平均した勾配エネルギー密度は $C_\nabla\epsilon^2 k^2/2$ である。$C_\nabla>0$ なら短波長の密度変調ほどコストが大きい。
:::

### 局所汎関数の変分と有効質量

$\int C_\nabla|\nabla\rho|^2\dd^3r$ の変分は $2C_\nabla\int\nabla\rho\cdot\nabla\delta\rho$ である。境界での変分をゼロにして部分積分すると、局所場に $-2C_\nabla\nabla^2\rho$ が加わる。一方、$\tau=\sum_i|\nabla\phi_i|^2$ に掛かる $B=\hbar^2/(2m)+C_\tau\rho$ は、軌道の微分を含むので $-\nabla\cdot(B\nabla)$ になる。係数 $B$ の空間微分を落とすと、一般には Hermite 性を保てない。運動項と電流項を同時に扱う場合、波動関数を $e^{i\chi}\phi_i$ に変えると $\vb j\mapsto\vb j+\rho\nabla\chi$、$\tau\mapsto\tau+2\vb j\cdot\nabla\chi+\rho|\nabla\chi|^2$ なので $\rho\tau-\vb j^2$ は不変である。

:::example[局所の有効質量]
均一な領域で $A_0=\hbar^2/(2m)$ とし、係数 $A_0+C_\tau\rho$ を用いる。

量の定義に従い、共通の単位で次の値を用いる。
$$
A_{0}=\frac{4147}{200}\,\mathrm{MeV\,fm^2},\qquad C_{t}=50\,\mathrm{MeV\,fm^5},\qquad \rho=\frac{2}{25}\,\mathrm{fm}^{-3}
$$

この条件で求める量を ${m^*/m}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{m^*/m}&=\frac{A_{0}}{A_{0} + C_{t} \cdot \rho}\\&=\frac{\left(\frac{4147}{200}\right)}{\left(50\right) \cdot \left(\frac{2}{25}\right) + \left(\frac{4147}{200}\right)}\\&\simeq 0.838286.\end{aligned}
$$

$C_\tau>0$ なら有効質量は裸質量より小さい。

飽和付近の密度と比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
\rho=\frac{4}{25}\,\mathrm{fm}^{-3}
$$

このときの値を ${m^*/m}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{m^*/m}_{\mathrm{new}}&=\frac{\left(\frac{4147}{200}\right)}{\left(50\right) \cdot \left(\frac{4}{25}\right) + \left(\frac{4147}{200}\right)}\\&\simeq 0.721594,\\{m^*/m}_{\mathrm{new}}-{m^*/m}&\simeq -0.116692.\end{aligned}
$$

**一次近似の検証。** 密度変化に対する局所応答を求めると、同じ $C_\tau$ でも傾きが初期密度で違うことが分かる。

${\rho}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${m^*/m}$ の一階微分を求めると、
$$
\frac{\mathrm d{m^*/m}}{\mathrm d\rho}=- \frac{A_{0} \cdot C_{t}}{\left(A_{0} + C_{t} \cdot \rho\right)^{2}}.
$$

この導関数を変更前の点で評価して $h=\frac{2}{25}$ を掛ける。ここで $\rho$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq -0.135563,\\\Delta-\Delta_{\mathrm{lin}}&\simeq 0.0188707.\end{aligned}
$$

残差が示すのは、指定した二点間で生じる線形化の誤差である。基準点からさらに離れた値について同じ精度を保証するものではない。
:::

:::example[密度依存項の再配列]
エネルギー密度 $C\rho^{\alpha+2}$ の一階微分を計算する。

量の定義に従い、共通の単位で次の値を用いる。
$$
C=150\,\mathrm{MeV\,fm^6},\qquad a=1,\qquad \rho=\frac{2}{25}\,\mathrm{fm}^{-3}
$$

この条件で求める量を ${U}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{U}&=C \cdot \rho^{a + 1} \cdot \left(a + 2\right)\\&=\left(150\right) \cdot \left(\frac{2}{25}\right)^{\left(1\right) + 1} \cdot \left(\left(1\right) + 2\right)\\&\simeq 2.88\,\mathrm{MeV}.\end{aligned}
$$

密度依存係数を固定して変分すると $\alpha$ に対応する係数を落とす。

密度を二倍にする。そこで次の値に変更し、ほかの入力を固定する。
$$
\rho=\frac{4}{25}\,\mathrm{fm}^{-3}
$$

このときの値を ${U}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{U}_{\mathrm{new}}&=\left(150\right) \cdot \left(\frac{4}{25}\right)^{\left(1\right) + 1} \cdot \left(\left(1\right) + 2\right)\\&\simeq 11.52\,\mathrm{MeV},\\{U}_{\mathrm{new}}-{U}&\simeq 8.64\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[密度波への表面エネルギー]
振幅 $a$ の密度波 $\rho=\rho_0+a\cos kx$ を周期境界で平均する。$\abs{\nabla \rho}^2$ の平均は $(ka)^2/2$ である。

量の定義に従い、共通の単位で次の値を用いる。
$$
C=40\,\mathrm{MeV\,fm^5},\qquad k=\frac{1}{2}\,\mathrm{fm}^{-1},\qquad a=\frac{1}{50}\,\mathrm{fm}^{-3}
$$

この条件で求める量を ${\Delta\mathcal E}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\Delta\mathcal E}&=\frac{C \cdot a^{2} \cdot k^{2}}{2}\\&=\frac{\left(40\right) \cdot \left(\frac{1}{2}\right)^{2} \cdot \left(\frac{1}{50}\right)^{2}}{2}\\&\simeq 0.002\,\mathrm{MeV\,fm}^{-3}.\end{aligned}
$$

同じ振幅で短い波を作るほど表面項のエネルギーが大きい。

波長を半分にする。そこで次の値に変更し、ほかの入力を固定する。
$$
k=1\,\mathrm{fm}^{-1}
$$

このときの値を ${\Delta\mathcal E}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\Delta\mathcal E}_{\mathrm{new}}&=\frac{\left(1\right)^{2} \cdot \left(40\right) \cdot \left(\frac{1}{50}\right)^{2}}{2}\\&\simeq 0.008\,\mathrm{MeV\,fm}^{-3},\\{\Delta\mathcal E}_{\mathrm{new}}-{\Delta\mathcal E}&\simeq 0.006\,\mathrm{MeV\,fm}^{-3}.\end{aligned}
$$
:::

### Coulomb 交換の局所近似

均一でスピン無偏極な陽子気体を考える。$a_e=e^2/(4\pi\varepsilon_0)$ とし、Coulomb 核の Fourier 成分を $4\pi a_e/q^2$ とする。ここでは交換項だけを計算する。二体和の係数 $1/2$ と、同じスピンの交換を二種類について数える因子 2 が相殺するので、単位体積当たりの交換エネルギーは
$$
\mathcal E_x=-\frac{4\pi a_e}{(2\pi)^6}
\int_{|\vb k|,|\vb k'|<k_{Fp}}\frac{\dd^3k\,\dd^3k'}{|\vb k-\vb k'|^2}
$$
である。$R=k_{Fp}$ とおく。差 $\vb q=\vb k-\vb k'$ を固定したとき、残る積分は半径 $R$ の二球を距離 $q$ だけずらした重なり体積 $W(q)$ である。重なりは高さ $R-q/2$ の球冠二つなので、断面積を積分すると $W(q)=2\pi\int_{q/2}^R(R^2-z^2)\dd z=\pi(4R+q)(2R-q)^2/12$ となる。$q>2R$ ではゼロである。$q$ 空間の体積要素 $4\pi q^2\dd q$ が分母を相殺し、$\int_0^{2R}W(q)\dd q=\pi R^4$ を使うと $\mathcal E_x=-a_e R^4/(4\pi^3)$ である。

陽子数密度 $\rho_p=R^3/(3\pi^2)$ を用い、核内の各位置を均一気体で近似すると $\mathcal E_x(\vb r)\simeq-(3/4)a_e(3/\pi)^{1/3}\rho_p(\vb r)^{4/3}$ となる。これが Coulomb 交換に対する Slater の局所近似であり、密度勾配の効果を省いている。

:::example[Coulomb Slater 交換]
局所近似で一陽子当たりの交換値を求める。$a_e=e^2/(4\pi\epsilon_0)$ とする。

量の定義に従い、共通の単位で次の値を用いる。
$$
a_{e}=\frac{36}{25}\,\mathrm{MeV\,fm},\qquad \rho_{p}=\frac{1}{25}\,\mathrm{fm}^{-3}
$$

この条件で求める量を ${e_x}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{e_x}&=- \frac{3 \cdot a_{e} \cdot \left(\frac{3 \cdot \rho_{p}}{\pi}\right)^{1 \cdot \frac{1}{3}}}{4}\\&=- \frac{3 \cdot \left(\frac{36}{25}\right) \cdot \left(\frac{3 \cdot \left(\frac{1}{25}\right)}{\pi}\right)^{\frac{1}{3}}}{4}\\&\simeq -0.36372\,\mathrm{MeV}.\end{aligned}
$$

交換項は負で、総 Coulomb エネルギーの直接寄与とは区別する。

陽子密度を二倍にする。そこで次の値に変更し、ほかの入力を固定する。
$$
\rho_{p}=\frac{2}{25}\,\mathrm{fm}^{-3}
$$

このときの値を ${e_x}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{e_x}_{\mathrm{new}}&=- \frac{3 \cdot \left(\frac{36}{25}\right) \cdot \left(\frac{3 \cdot \left(\frac{2}{25}\right)}{\pi}\right)^{\frac{1}{3}}}{4}\\&\simeq -0.458259\,\mathrm{MeV},\\{e_x}_{\mathrm{new}}-{e_x}&\simeq -0.0945386\,\mathrm{MeV}.\end{aligned}
$$
:::

## 相対論的平均場理論の概観

相対論的模型では波動関数は Dirac スピノルであり、スカラー場 $S$ と時間成分のベクトル場 $V$ を区別する。一様静止物質で有効質量 $m_D^*c^2=mc^2+S$ とすると、一粒子 Hamiltonian は $H_D=c\vb\alpha\cdot\vb p+\beta m_D^*c^2+V$ となる。Dirac 行列の反交換関係で $(H_D-V)^2=c^2p^2+(m_D^*c^2)^2$ となるため、正エネルギー枝は
$$
\varepsilon(p)=V+\sqrt{c^2p^2+(mc^2+S)^2}.
$$
$S,V$ は核子のスカラー密度と数密度を源として自己無撞着に求める。例えば自然単位の静的線形模型で $( -\nabla^2+m_s^2)S=-g_s^2\rho_s$、$(-\nabla^2+m_v^2)V=g_v^2\rho$ と書けば、一様物質の場は $S=-g_s^2\rho_s/m_s^2$、$V=g_v^2\rho/m_v^2$ になる。有限核では場の空間変化からスピン軌道項が現れる。負エネルギー状態や自己相互作用を含む模型の定式化は場の理論を必要とする。

:::example[引力と斥力の差]
$S=-350\,\mathrm{MeV},V=300\,\mathrm{MeV}$ で $mc^2=939\,\mathrm{MeV}$ なら、$p=0$ の正エネルギーの移動は $-50\,\mathrm{MeV}$ である。一方で Dirac 質量に相当するエネルギーは $589\,\mathrm{MeV}$ になる。
:::

:::example[低運動量展開]
$m_D^*>0$ とし $p\ll m_D^*c$ なら、$\varepsilon=V+m_D^*c^2+p^2/(2m_D^*)+O(p^4)$ である。分散の曲率には $V$ が直接入らず、静止点の移動には両方の場が入る。
:::

### スカラー場とベクトル場でのスピン軌道項

球対称な $S,V$ で Dirac スピノルを上成分 $\phi$ と下成分 $\chi$ に分ける。全エネルギーを $\varepsilon$ とすると、下側の式は $\chi=D^{-1}c\vb\sigma\cdot\vb p\,\phi$、$D=\varepsilon-V+mc^2+S$ となる。これを上側へ代入すると $c^2(\vb\sigma\cdot\vb p)D^{-1}(\vb\sigma\cdot\vb p)\phi=(\varepsilon-V-mc^2-S)\phi$ である。$D^{-1}=g(r)$ と書けば、$p_i g=g p_i-i\hbar\partial_i g$ と Pauli 行列の積から
$$
(\vb\sigma\cdot\vb p)g(\vb\sigma\cdot\vb p)
=g\vb p^2-i\hbar\grad g\cdot\vb p
+\hbar\vb\sigma\cdot(\grad g\times\vb p)
$$
となる。最後の項は $\hbar g'(r)\vb\sigma\cdot\vb L/r$ である。$g'=(V'-S')/D^2$ なので、スピン軌道項の符号と強さには場の差の勾配が関わる。全エネルギーの低運動量での移動には $S+V$ が入ったことと比較する。中心項の和が小さくても差の勾配は大きくなり得る。核の表面で $S$ が増え $V$ が減る場合はスピン軌道の係数が負となり、:eqref[eq-nuc-shell-spin-orbit] の典型的な準位順に対応する。

:::example[一定の場と表面]
$S$ と $V$ が定数なら $D^{-1}$ が運動量演算子と可換であり、交差積の項は消える。同じ内側の場をもつ二つの模型でも、表面の変化する幅を変えれば $g'$ の分布が変わる。核の束縛を合わせることとスピン軌道分裂を合わせることには別の制約がある。
:::

### Dirac 質量と Fermi 面の速度

静的で一様な場では $H=c\boldsymbol\alpha\cdot\vb p+\beta(mc^2+S)+V$ である。$\boldsymbol\alpha,\beta$ の反交換関係を使って $(H-V)^2$ を計算すると、交差項がすべて消えて $p^2c^2+(mc^2+S)^2$ となる。正エネルギー枝を微分すると $v/c=pc/(\varepsilon-V)$ である。準粒子の速度から定義する :term[Landau 質量]{en="Landau mass" yomi="らんだうしつりょう"} を $m_L=p_F/v_F$ とすると $m_Lc^2=\sqrt{(p_Fc)^2+(mc^2+S)^2}$ となる。従って $m_D= m+S/c^2$ は Dirac 方程式の質量で、Fermi 面の速度により測る $m_L$ とは異なる。

:::example[静止点のエネルギー移動]
$m_D>0$ とし、正枝のゼロ運動量で自由核子からのずれを求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
S=-350\,\mathrm{MeV},\qquad V=300\,\mathrm{MeV}
$$

この条件で求める量を ${\Delta\epsilon_0}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\Delta\epsilon_0}&=S + V\\&=\left(-350\right) + \left(300\right)\\&=-50\,\mathrm{MeV}.\end{aligned}
$$

二つの場がそれぞれ大きくてもゼロ運動量での移動は相殺で小さい。

ベクトル場だけを変える。そこで次の値に変更し、ほかの入力を固定する。
$$
V=310\,\mathrm{MeV}
$$

このときの値を ${\Delta\epsilon_0}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\Delta\epsilon_0}_{\mathrm{new}}&=\left(-350\right) + \left(310\right)\\&=-40\,\mathrm{MeV},\\{\Delta\epsilon_0}_{\mathrm{new}}-{\Delta\epsilon_0}&=10\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[速度を決める有効質量]
Fermi 運動量のエネルギー換算 $p_c$ と Dirac 質量 $m_c$ から Landau 質量を求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
p_{c}=260\,\mathrm{MeV},\qquad m_{c}=600\,\mathrm{MeV}
$$

この条件で求める量を ${m_Lc^2}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{m_Lc^2}&=\sqrt{m_{c}^{2} + p_{c}^{2}}\\&=\sqrt{\left(260\right)^{2} + \left(600\right)^{2}}\\&\simeq 653.911\,\mathrm{MeV}.\end{aligned}
$$

一定のベクトル場は速度の微分から消える。

Dirac 質量が大きい場合を比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
m_{c}=650\,\mathrm{MeV}
$$

このときの値を ${m_Lc^2}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{m_Lc^2}_{\mathrm{new}}&=\sqrt{\left(260\right)^{2} + \left(650\right)^{2}}\\&\simeq 700.071\,\mathrm{MeV},\\{m_Lc^2}_{\mathrm{new}}-{m_Lc^2}&\simeq 46.1601\,\mathrm{MeV}.\end{aligned}
$$

**一次近似の検証。** Dirac 質量の誤差が運動学に伝わる割合を、Fermi 運動量を固定して求める。

${m_{c}}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${m_Lc^2}$ の一階微分を求めると、
$$
\frac{\mathrm d{m_Lc^2}}{\mathrm dm_{c}}=\frac{m_{c}}{\sqrt{m_{c}^{2} + p_{c}^{2}}}.
$$

この導関数を変更前の点で評価して $h=50$ を掛ける。ここで $m_{c}$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq 45.8778\,\mathrm{MeV},\\\Delta-\Delta_{\mathrm{lin}}&\simeq 0.282335\,\mathrm{MeV}.\end{aligned}
$$

残差が示すのは、指定した二点間で生じる線形化の誤差である。基準点からさらに離れた値について同じ精度を保証するものではない。
:::

:::example[正枝の運動エネルギー]
同じ場の静止点からの差を取り、相対論的分散を使う。

量の定義に従い、共通の単位で次の値を用いる。
$$
p_{c}=200\,\mathrm{MeV},\qquad m_{c}=600\,\mathrm{MeV}
$$

この条件で求める量を ${T}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{T}&=- m_{c} + \sqrt{m_{c}^{2} + p_{c}^{2}}\\&=- \left(600\right) + \sqrt{\left(200\right)^{2} + \left(600\right)^{2}}\\&\simeq 32.4555\,\mathrm{MeV}.\end{aligned}
$$

速度が大きくなると二次近似との差が増す。

高い運動量を計算する。そこで次の値に変更し、ほかの入力を固定する。
$$
p_{c}=400\,\mathrm{MeV}
$$

このときの値を ${T}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{T}_{\mathrm{new}}&=- \left(600\right) + \sqrt{\left(400\right)^{2} + \left(600\right)^{2}}\\&\simeq 121.11\,\mathrm{MeV},\\{T}_{\mathrm{new}}-{T}&\simeq 88.6547\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[Lorentz スカラー密度の重み]
正エネルギー状態の運動量で $\bar uu/u^\dagger u=m_c/\sqrt{p_c^2+m_c^2}$ になる。

量の定義に従い、共通の単位で次の値を用いる。
$$
m_{c}=600\,\mathrm{MeV},\qquad p_{c}=0\,\mathrm{MeV}
$$

この条件で求める量を ${w}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{w}&=\frac{m_{c}}{\sqrt{m_{c}^{2} + p_{c}^{2}}}\\&=\frac{\left(600\right)}{\sqrt{\left(0\right)^{2} + \left(600\right)^{2}}}\\&=1.\end{aligned}
$$

ベクトル数密度の重み 1 とは運動中に違いが生じる。

静止と Fermi 面付近を比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
p_{c}=300\,\mathrm{MeV}
$$

このときの値を ${w}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{w}_{\mathrm{new}}&=\frac{\left(600\right)}{\sqrt{\left(300\right)^{2} + \left(600\right)^{2}}}\\&\simeq 0.894427,\\{w}_{\mathrm{new}}-{w}&\simeq -0.105573.\end{aligned}
$$
:::

### 一様なスカラー場とベクトル場の自己無撞着条件

この段に限り $\hbar=c=1$ とする。スカラー場 $\phi$ と時間成分だけのベクトル場 $V_0$ を導入し、Dirac 質量を $m_D=m-g_s\phi$、正枝の一体エネルギーを $\sqrt{k^2+m_D^2}+g_v V_0$ とする。近似として Dirac 海と空間勾配を捨て、二つの場の自由質量を $m_s,m_v$ とする。一様な物質のエネルギー密度を、ベクトル場の一体エネルギーを含めて
$$
\mathcal E=\frac{g}{2\pi^2}\int_0^{k_F} k^2\sqrt{k^2+m_D^2}\dd k
+\frac{m_s^2}{2}\phi^2+g_vV_0\rho-\frac{m_v^2}{2}V_0^2
$$
と置く。$V_0$ は制約を課す時間成分であり、この表式をどの場についても最小化するという意味ではない。密度を固定した変分は
$$
m_s^2\phi=g_s\rho_s,\quad
\rho_s=\frac{g}{2\pi^2}\int_0^{k_F}\frac{m_Dk^2}{\sqrt{k^2+m_D^2}}\dd k,\quad
m_v^2V_0=g_v\rho
$$
を与える。後者を消去したベクトル寄与は $g_v^2\rho^2/(2m_v^2)>0$ である。スカラー場の源 $\rho_s$ には運動する核子の重み $m_D/\sqrt{k^2+m_D^2}$ が入るため、数密度と同じようには増えない。これが引力と斥力の異なる密度依存を生む。実際の有限核では勾配、等ベクトル場や非線形な項がこの簡略模型に加わる。

## 対称エネルギーと核物質の状態方程式

Coulomb を除き荷電対称性を仮定する。非対称度 $\delta=(\rho_n-\rho_p)/\rho$ を反転してもエネルギーは同じなので、$\delta=0$ 周辺で十分滑らかなら
$$
e(\rho,\delta)=e(\rho,0)+S(\rho)\delta^2+O(\delta^4),
\qquad S(\rho)=\tfrac12\partial_\delta^2 e(\rho,0).
\label{eq-nuc-symmetry-energy}
$$
これを :term[対称エネルギー]{en="symmetry energy" yomi="たいしょうえねるぎー"} と呼ぶ。:eqref[eq-nuc-fermi-gas] の種ごとの和を展開すると運動寄与 $S_{\rm kin}=E_F/3$ を得る。$S$ の値には相互作用が加わる。飽和密度での傾き $L=3\rho_0 S'(\rho_0)$ は圧力と関係する。一次の密度微分をとると、対称物質の圧力がゼロである $\rho_0$ では $P(\rho_0,\delta)=\rho_0 L\delta^2/3+O(\delta^4)$ である。

粒子種ごとの化学ポテンシャルを $\mu_q=\partial(\rho e)/\partial\rho_q$ とする。中性子を増やし同数の陽子を減らす変化では $\rho$ が固定され、$\dd\delta=2\dd\rho_n/\rho$ なので $\mu_n-\mu_p=2\partial_\delta e=4S\delta+O(\delta^3)$ となる。静止質量差を入れると $(m_n-m_p)c^2$ が加わる。状態方程式は $P$ やエネルギー密度を $\rho$、温度、組成の関数として与える。中性子星で組成を求めるには、別途、荷電中性と弱い相互作用の平衡を用いる。

:::example[運動寄与]
$E_F=37\,\mathrm{MeV}$ なら $S_{\rm kin}\simeq12.3\,\mathrm{MeV}$ である。模型の全係数が $32\,\mathrm{MeV}$ なら残りの約 $19.7\,\mathrm{MeV}$ が相互作用分である。
:::

:::example[傾きから圧力へ]
$\rho_0=0.16\,\mathrm{fm}^{-3}$、$L=60\,\mathrm{MeV}$、$\delta=1/2$ では二次近似の圧力は $0.8\,\mathrm{MeV\,fm}^{-3}$ である。$L$ は対称物質の曲率 $K_0$ とは別の微分を測っている。
:::

### 密度と圧縮率の独立な情報

単純な状態方程式 $e(\rho,\delta)=e_0+(K_0/18)x^2+J\delta^2+(L/3)x\delta^2$ を考える。ここで $x=(\rho-\rho_0)/\rho_0$ である。$\delta$ を固定して極小の条件を課すと $x_{\rm sat}=-3L\delta^2/K_0$ が得られる。従って $L>0$ の模型では、不均衡を増すとこの近似の飽和密度が下がる。導出には $K_0>0$ と $|x|,|\delta|\ll1$ が必要である。中性子だけの系まで式を伸ばして極小を推定するには、捨てた高次項と相の変化を調べる必要がある。

:::example[飽和点の移動]
$K_0=240\,\mathrm{MeV},L=60\,\mathrm{MeV},\delta=0.2$ では $x_{\rm sat}=-0.03$ である。$\rho_0=0.16\,\mathrm{fm}^{-3}$ とすると模型内の密度変化は $-0.0048\,\mathrm{fm}^{-3}$ である。
:::

### 物質の局所安定性

冷たい一成分物質で全エネルギー密度を $\mathcal E=\rho(mc^2+e)$ とする。化学ポテンシャルは $\mu=\partial_\rho\mathcal E=mc^2+e+\rho e'$ である。:eqref[eq-nuc-matter-pressure] と比較すると $P=\rho\mu-\mathcal E$ になる。微分して $\partial_\rho P=\rho\,\partial_\rho\mu$ であるから、圧力が密度に対して低下する領域では一様物質に微小な密度差を作ると成長し得る。安定な単一相で組成を固定した小さい断熱波の音速は $c_s^2=c^2\,\partial P/\partial\mathcal E$ である。相対論的な密度に模型を外挿するときは $c_s\le c$ も点検する。飽和点に合わせた非相対論的な多項式が、高密度でこの条件を保証するわけではない。

:::example[局所曲率の比較]
飽和点で $e_0=-16\,\mathrm{MeV}$、$mc^2=939\,\mathrm{MeV}$、$K_0=240\,\mathrm{MeV}$ とする。$P'=K_0/9$、$\mathcal E'=923\,\mathrm{MeV}$ なので、共通の等温・断熱極限で $c_s/c\simeq\sqrt{240/(9\cdot923)}\simeq0.17$ となる。これは対称核物質の飽和点での値であり、中性子星の全領域で一定ではない。
:::

### 非対称物質の化学ポテンシャル

電磁力と核子質量差を無視し、$\rho=\rho_n+\rho_p$、$\delta=(\rho_n-\rho_p)/\rho$ とする。$e=e_s(\rho)+S(\rho)\delta^2$ なら $\mathcal E=\rho e$ である。$\partial_{\rho_n}\delta=(1-\delta)/\rho$、$\partial_{\rho_p}\delta=-(1+\delta)/\rho$ なので、化学ポテンシャル差は $\mu_n-\mu_p=4S\delta$ となる。$\partial_\rho(\rho e)$ の共通部分は差では消える。圧力は固定組成で $P=\rho^2\partial_\rho e$ である。反応で組成が変わる場合は、反応平衡と中性条件を課してから状態方程式の経路を求める。

:::example[中性子と陽子の化学ポテンシャル差]
二次近似で $S$ と非対称度 $\delta$ を使う。

量の定義に従い、共通の単位で次の値を用いる。
$$
S=32\,\mathrm{MeV},\qquad \delta=\frac{1}{4}
$$

この条件で求める量を ${\Delta\mu}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\Delta\mu}&=4 \cdot S \cdot \delta\\&=4 \cdot \left(32\right) \cdot \left(\frac{1}{4}\right)\\&=32\,\mathrm{MeV}.\end{aligned}
$$

陽子の少ない物質ほど、同じ密度なら差は大きい。

より中性子過剰な組成に変える。そこで次の値に変更し、ほかの入力を固定する。
$$
\delta=\frac{3}{4}
$$

このときの値を ${\Delta\mu}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\Delta\mu}_{\mathrm{new}}&=4 \cdot \left(32\right) \cdot \left(\frac{3}{4}\right)\\&=96\,\mathrm{MeV},\\{\Delta\mu}_{\mathrm{new}}-{\Delta\mu}&=64\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[自由気体の対称エネルギー]
対称系の $E_F$ を使う。$(1+\delta)^{5/3}+(1-\delta)^{5/3}$ の二次係数から $S_{\rm kin}=E_F/3$ である。

量の定義に従い、共通の単位で次の値を用いる。
$$
E_{F}=36\,\mathrm{MeV}
$$

この条件で求める量を ${S_{\rm kin}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{S_{\rm kin}}&=\frac{E_{F}}{3}\\&=\frac{\left(36\right)}{3}\\&=12\,\mathrm{MeV}.\end{aligned}
$$

経験的な全対称エネルギーには相互作用の成分も加わる。

高密度側の運動成分と比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
E_{F}=48\,\mathrm{MeV}
$$

このときの値を ${S_{\rm kin}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{S_{\rm kin}}_{\mathrm{new}}&=\frac{\left(48\right)}{3}\\&=16\,\mathrm{MeV},\\{S_{\rm kin}}_{\mathrm{new}}-{S_{\rm kin}}&=4\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[非対称圧力の傾き依存]
飽和点で $L=3\rho_0S'(\rho_0)$ とおく。二次の非対称圧力だけを計算する。

量の定義に従い、共通の単位で次の値を用いる。
$$
\rho_{0}=\frac{4}{25}\,\mathrm{fm}^{-3},\qquad L=60\,\mathrm{MeV},\qquad \delta=\frac{4}{5}
$$

この条件で求める量を ${P_s}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{P_s}&=\frac{L \cdot \delta^{2} \cdot \rho_{0}}{3}\\&=\frac{\left(60\right) \cdot \left(\frac{4}{25}\right) \cdot \left(\frac{4}{5}\right)^{2}}{3}\\&\simeq 2.048\,\mathrm{MeV\,fm}^{-3}.\end{aligned}
$$

同じ $J=S(\rho_0)$ でも傾きが違えば圧力は変わる。

大きな傾きの模型と比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
L=90\,\mathrm{MeV}
$$

このときの値を ${P_s}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{P_s}_{\mathrm{new}}&=\frac{\left(90\right) \cdot \left(\frac{4}{25}\right) \cdot \left(\frac{4}{5}\right)^{2}}{3}\\&\simeq 3.072\,\mathrm{MeV\,fm}^{-3},\\{P_s}_{\mathrm{new}}-{P_s}&\simeq 1.024\,\mathrm{MeV\,fm}^{-3}.\end{aligned}
$$
:::

:::example[高次の組成効果]
エネルギーに $S_4\delta^4$ を加える。化学ポテンシャル差への補正は $2\partial_\delta e$ なので $8S_4\delta^3$ である。

量の定義に従い、共通の単位で次の値を用いる。
$$
S_{4}=2\,\mathrm{MeV},\qquad \delta=\frac{1}{4}
$$

この条件で求める量を ${\Delta\mu_4}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\Delta\mu_4}&=8 \cdot S_{4} \cdot \delta^{3}\\&=8 \cdot \left(2\right) \cdot \left(\frac{1}{4}\right)^{3}\\&\simeq 0.25\,\mathrm{MeV}.\end{aligned}
$$

ほぼ対称な系で小さな寄与でも、大きい非対称度では精度に効く。

同じ四次係数で組成の変化を見る。そこで次の値に変更し、ほかの入力を固定する。
$$
\delta=\frac{3}{4}
$$

このときの値を ${\Delta\mu_4}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\Delta\mu_4}_{\mathrm{new}}&=8 \cdot \left(2\right) \cdot \left(\frac{3}{4}\right)^{3}\\&\simeq 6.75\,\mathrm{MeV},\\{\Delta\mu_4}_{\mathrm{new}}-{\Delta\mu_4}&\simeq 6.5\,\mathrm{MeV}.\end{aligned}
$$

**一次近似の検証。** 四次の項そのものではなく、その化学ポテンシャル差に対する非対称度の感度を調べる。

${\delta}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${\Delta\mu_4}$ の一階微分を求めると、
$$
\frac{\mathrm d{\Delta\mu_4}}{\mathrm d\delta}=24 \cdot S_{4} \cdot \delta^{2}.
$$

この導関数を変更前の点で評価して $h=\frac{1}{2}$ を掛ける。ここで $\delta$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq 1.5\,\mathrm{MeV},\\\Delta-\Delta_{\mathrm{lin}}&=5\,\mathrm{MeV}.\end{aligned}
$$

省いた分が全変化に比べて大きい。従って増分の大きさを予測するには、少なくとも二次の項を調べるか、元の式を使う必要がある。
:::

## 演習問題

::::exercise{level=basic}
対称 Fermi 気体で密度を 8 倍にしたとき $k_F$ と $E_F$ は何倍か。
:::solution
:eqref[eq-nuc-fermi-gas] から $k_F\propto\rho^{1/3}$ なので 2 倍である。二乗に比例する $E_F$ は 4 倍である。
:::
::::

::::exercise{level=basic}
$e=e_0+K_0[(\rho-\rho_0)/\rho_0]^2/18$ で $\rho=\rho_0$ の圧力を求めよ。
:::solution
微分は $K_0(\rho-\rho_0)/(9\rho_0^2)$ なので $\rho=\rho_0$ でゼロである。:eqref[eq-nuc-matter-pressure] から $P=0$ となる。
:::
::::

::::exercise{level=basic}
一軌道だけを占める Slater 行列式で二体エネルギーがゼロとなることを示せ。
:::solution
二体和に残る添字は $i=j$ だけであり、直接項 $\bra{ii}v\ket{ii}$ と交換項が同じで打ち消す。一粒子には対がない。
:::
::::

::::exercise{level=basic}
$\rho_n=0.12$、$\rho_p=0.04\,\mathrm{fm}^{-3}$、$S=32\,\mathrm{MeV}$ のとき二次の対称エネルギー寄与を求めよ。
:::solution
$\rho=0.16$、$\delta=0.08/0.16=1/2$ なので一核子当たり $S\delta^2=8\,\mathrm{MeV}$ である。
:::
::::

::::exercise{level=advanced}
$e=a\rho^{2/3}$ の気体で $P$ をエネルギー密度 $\rho e$ で表せ。
:::solution
$e'=(2/3)a\rho^{-1/3}$ なので $P=(2/3)a\rho^{5/3}=(2/3)\rho e$ である。この式には粒子の静止エネルギーを含めていない。
:::
::::

::::exercise{level=advanced}
いずれも一体基底である占有 $i$ と非占有 $a$ の混合による一次変分から、$h_{ai}=0$ を導け。
:::solution
$\phi_i\to\phi_i+\epsilon\phi_a$ は一次まで規格化を変えない。密度の変分は $\epsilon\ket a\bra i+\epsilon^*\ket i\bra a$ であり、$\delta E=2\Re(\epsilon^* h_{ai})$ である。任意の実部と虚部に対しゼロとなるには $h_{ai}=0$ が必要である。
:::
::::

::::exercise{level=advanced}
一様物質で $S(\rho)=J(\rho/\rho_0)^\gamma$ と仮定する。$L$ と $\rho_0$ での二次の非対称圧力を求めよ。
:::solution
$S'(\rho_0)=\gamma J/\rho_0$ なので $L=3\gamma J$ である。$P=\rho^2\partial_\rho e$ から、$\rho_0$ では対称成分がゼロとなり、非対称圧力は $\gamma J\rho_0\delta^2$ である。
:::
::::

::::exercise{level=advanced}
$p=0$ の正エネルギーを保ったまま $S$ を $\Delta S$ 変えたとき $V$ と分散の二次係数はどう変わるか。$m_D^*>0$ とする。
:::solution
$\varepsilon(0)=V+mc^2+S$ なので $V$ は $-\Delta S$ 変える。二次係数は $1/(2m_D^*)$ から $1/[2(m_D^*+\Delta S/c^2)]$ へ移る。静止点が同じでも分散は変化する。
:::
::::

## 確認テスト

::::quiz{answer="22.2" tol="1e-6"}
$E_F=37\,\mathrm{MeV}$ の非相対論 Fermi 気体の一粒子運動エネルギーを MeV で求めよ。
:::explanation
球内で $\int k^4\dd k/\int k^2\dd k=(3/5)k_F^2$ であり、平均は $(3/5)37=22.2\,\mathrm{MeV}$ である。
:::
::::

::::quiz
飽和密度での対称核物質について正しいものを選べ。

- [x] 一粒子エネルギーの密度微分がゼロになる。
- [ ] 密度がゼロになる。
- [ ] 非圧縮率も必ずゼロである。

:::explanation
有限な $\rho_0$ でエネルギーが極小になり、一階微分がゼロになる。安定な二次極小では二階微分は正である。
:::
::::

::::quiz{answer="60" tol="1e-6"}
$S=30(\rho/\rho_0)^{2/3}\,\mathrm{MeV}$ の $L$ を MeV で求めよ。
:::explanation
飽和点の微分は $30(2/3)/\rho_0$ である。$3\rho_0$ を掛けて $60\,\mathrm{MeV}$ を得る。
:::
::::

::::quiz
Hartree–Fock で交換項を加える理由を選べ。

- [x] 同種の核子の交換に対して多体状態を反対称にするため。
- [ ] 運動エネルギーを倍に数えるため。
- [ ] 異なる軌道の間に Coulomb 相互作用だけを残すため。

:::explanation
Slater 行列式の反対称化により、二体密度に直接項と負の交換項が現れる。交換項は Coulomb 以外の二体演算子にも現れる。
:::
::::

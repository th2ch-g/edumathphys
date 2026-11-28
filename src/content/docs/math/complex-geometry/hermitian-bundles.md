---
title: Hermite 束と Chern 接続
description: 正則構造と計量から接続を求め、その曲率を Chern 類と正値性に結び付ける
notation:
  - symbol: 'H'
    meaning: 局所正則枠における Hermite 計量の行列
  - symbol: 'F_\nabla'
    meaning: 接続 ∇ の曲率
  - symbol: 'c(E,\nabla)'
    meaning: 曲率で作る全 Chern 形式
  - symbol: '\Lambda_\omega'
    meaning: Kähler 形式との外積の随伴
---

## この章の目標

- 平滑な Hermite 計量を束上に作り、枠の変更でどう表れるかを示す。
- Chern 接続の存在と一意性を局所行列で証明する。
- 曲率の微分形式から第一 Chern 類と高次 Chern 類を計算する。
- 正直線束と Hermite–Einstein 条件を具体例で区別する。

## Hermite 計量

:::definition[Hermite 計量]{#def-cg-hermitian-metric}
複素ベクトル束の各ファイバーに正定値 Hermite 内積を取り、その成分が底点に滑らかに依存するものを :term[Hermite 計量]{en="Hermitian metric" yomi="えるみーとけいりょう"} という。ここでは第一変数に線形とし、列ベクトル $v,w$ に対して $h(v,w)=w^*Hv$ と書く。
:::

局所枠の変更を $e'=eg$ とすると古い係数は $gv'$ なので、$H'=g^*Hg$ になる。任意の滑らかな複素ベクトル束には Hermite 計量が存在する。局所積で標準内積を取り、$1$ の分割で足すと、正定値行列の正係数和は正定値だからである。滑らかな正規直交枠は Gram–Schmidt で作れるが、その枠は一般に正則ではない。

:::example[計量と枠の変更]
直線束の正則枠 $e$ で $\|e\|^2=e^{-|z|^2}$ とする。$e'=e^{2z}e$ では $\|e'\|^2=e^{-|z|^2+4\Re z}$ となる。対数計量の差 $4\Re z$ は調和なので、$\partial\bar\partial$ を取ると消える。
:::

### 商計量と直交補を構成する

$S\subset E$ を正則部分束とする。滑らかな直交分解 $E=S\oplus S^\perp$ から、商束 $Q$ の元は $S^\perp$ にある唯一の代表を持つ。この代表のノルムを商計量とする。正則枠のうち最初の列が $S$ を張るとき、計量行列を $\begin{pmatrix}A&B\\B^*&C\end{pmatrix}$ と分ける。商の成分 $v$ に対して直交条件は $Au+Bv=0$ だから、代表は $(-A^{-1}Bv,v)$ である。ノルムに代入すると商計量は $C-B^*A^{-1}B$ となる。この行列の二次形式は、$v\ne0$ に対して、直交代表のノルムの二乗を与えるので正定値である。

:::example[商で非対角項が引かれる]
$H=\begin{pmatrix}2&1+i\\1-i&4\end{pmatrix}$ で第一座標を $S$ とすると、商係数 $v$ の代表は $(-(1+i)v/2,v)$ である。商計量は $4-|1+i|^2/2=3$ なので、$v=2$ のノルムの二乗は $12$ となる。標準代表 $(0,2)$ のノルムの二乗は $16$ であり、直交代表がこの中で最小になる。
:::

### 行列表示に現れる正値性

:::example[非対角計量の Chern 接続を求める]
正則な枠で計量行列を次のように置く。計算では $\bar z$ を $b$ と書き、$z,b$ を独立に微分する。点 $z=1$ で曲率の係数を求め、計量との両立を確認する。

$$
\begin{aligned}
H&=\left[\begin{matrix}b z + 2 & b\\z & 1\end{matrix}\right] \\
\det H&=2
\end{aligned}
$$

曲率の $(1,1)$ 成分を $C\,\mathrm dz\wedge\mathrm d\bar z$ と書く。計量適合性から正則枠での接続形式は $A\,\mathrm dz=H^{-1}\partial H$ である。$\bar\partial$ の後では外積を入れ替えるため、係数は $-\partial_b A$ になる。

$$
\begin{aligned}
\partial_z H&=\left[\begin{matrix}b & 0\\1 & 0\end{matrix}\right] \\
A&=\left[\begin{matrix}0 & 0\\1 & 0\end{matrix}\right]
\end{aligned}
$$

$$
\begin{aligned}
C&=\left[\begin{matrix}0 & 0\\0 & 0\end{matrix}\right] \\
C(p)&=\left[\begin{matrix}0 & 0\\0 & 0\end{matrix}\right],\qquad \tr C(p)=0
\end{aligned}
$$

実方向での曲率の斜 Hermite 性は、この枠の $(1,1)$ 係数について $C(p)^*H(p)=H(p)C(p)$ となる。実際、両辺は $\left[\begin{matrix}0 & 0\\0 & 0\end{matrix}\right]$ である。跡についても $\tr C=-\partial_z\partial_b\log\det H$ で確かめられる。この計量は行列式が定数で、全曲率も零である。実際、正則な行列を使って定数計量から移したものだからである。
:::

:::example[非対角計量の Chern 接続を求める]
正則な枠で計量行列を次のように置く。計算では $\bar z$ を $b$ と書き、$z,b$ を独立に微分する。点 $z=i$ で曲率の係数を求め、計量との両立を確認する。

$$
\begin{aligned}
H&=\left[\begin{matrix}2 & b\\z & b z + 2\end{matrix}\right] \\
\det H&=b z + 4
\end{aligned}
$$

曲率の $(1,1)$ 成分を $C\,\mathrm dz\wedge\mathrm d\bar z$ と書く。計量適合性から正則枠での接続形式は $A\,\mathrm dz=H^{-1}\partial H$ である。$\bar\partial$ の後では外積を入れ替えるため、係数は $-\partial_b A$ になる。

$$
\begin{aligned}
\partial_z H&=\left[\begin{matrix}0 & 0\\1 & b\end{matrix}\right] \\
A&=\left[\begin{matrix}- \frac{b}{b z + 4} & - \frac{b^{2}}{b z + 4}\\\frac{2}{b z + 4} & \frac{2 b}{b z + 4}\end{matrix}\right]
\end{aligned}
$$

$$
\begin{aligned}
C&=\left[\begin{matrix}\frac{4}{\left(b z + 4\right)^{2}} & \frac{b \left(b z + 8\right)}{\left(b z + 4\right)^{2}}\\\frac{2 z}{\left(b z + 4\right)^{2}} & - \frac{8}{\left(b z + 4\right)^{2}}\end{matrix}\right] \\
C(p)&=\left[\begin{matrix}\frac{4}{25} & - \frac{9 i}{25}\\\frac{2 i}{25} & - \frac{8}{25}\end{matrix}\right],\qquad \tr C(p)=- \frac{4}{25}
\end{aligned}
$$

実方向での曲率の斜 Hermite 性は、この枠の $(1,1)$ 係数について $C(p)^*H(p)=H(p)C(p)$ となる。実際、両辺は $\left[\begin{matrix}\frac{2}{5} & - \frac{2 i}{5}\\\frac{2 i}{5} & - \frac{3}{5}\end{matrix}\right]$ である。跡についても $\tr C=-\partial_z\partial_b\log\det H$ で確かめられる。計量行列の固有値と、曲率の係数行列の固有値を混同しない。計量には正定値性が求められるが、曲率にはそのような制約がない。
:::

## Chern 接続

接続 $\nabla$ は $\nabla(fs)=\dd f\otimes s+f\nabla s$ を満たす一階作用素である。局所枠で $\nabla=\dd+A$ と書く。この規約では曲率は $\nabla^2=F_\nabla=\dd A+A\wedge A$ である。:ref[def-rg-bundle-connection] を参照する。

:::theorem[Chern 接続の存在と一意性]{#thm-cg-chern-connection}
正則 Hermite ベクトル束 $(E,h)$ に対し、計量と両立し $(0,1)$ 成分が $\bar\partial_E$ である接続がただ一つ存在する。これを :term[Chern 接続]{en="Chern connection" yomi="ちぇるんせつぞく"} と呼ぶ。正則枠では $A=H^{-1}\partial H$、$F_\nabla=\bar\partial(H^{-1}\partial H)$ である。
:::

:::proof
正則枠では $\bar\partial_E$ は係数ごとの $\bar\partial$ であり、条件は $A$ に $(0,1)$ 成分がないことである。計量との両立は $\dd H=A^*H+HA$ である。$(1,0)$ 成分から $A=H^{-1}\partial H$ が必要となり、共役転置で $(0,1)$ 成分も成立する。正則な $g$ で枠を変えると $H'=g^*Hg$、$\partial g^*=0$ だから
$$
(H')^{-1}\partial H'=g^{-1}Ag+g^{-1}\partial g.
$$
これは接続の遷移式なので全体に貼り合う。最後に $\partial H^{-1}=-H^{-1}(\partial H)H^{-1}$ を使うと $\partial A=-A\wedge A$ になり、曲率は $\bar\partial A$ で型 $(1,1)$ である。
:::

:::example[非定数計量でも曲率が零になる場合]
自明直線束で $H=e^{z+\bar z}$ なら $A=\dd z$ で $F_\nabla=0$ である。枠 $e'=e^{-z}e$ は正則かつ $\|e'\|=1$ となり、そこで $A'=0$ である。
:::

:::example[直和の対角計算]
$H=\diag(e^{-|z|^2},e^{-2|z|^2})$ なら $A=\diag(-\bar z,-2\bar z)\dd z$ であり、$F=\diag(1,2)\dd z\wedge\dd\bar z$ になる。第二成分の計量を一定倍しても $A$ と $F$ は変わらない。
:::

### 正則枠の中での非対角の計算

自明な階数二の束に
$$
H=\begin{pmatrix}1&\bar z\\ z&1+|z|^2\end{pmatrix}
$$
を入れる。第一主小行列式と行列式がともに $1$ なので正定値である。逆行列は $\left(\begin{smallmatrix}1+|z|^2&-\bar z\\-z&1\end{smallmatrix}\right)$ であり、$\partial H=\left(\begin{smallmatrix}0&0\\1&\bar z\end{smallmatrix}\right)\dd z$ だから
$$
A=\begin{pmatrix}-\bar z&-\bar z^2\\1&\bar z\end{pmatrix}\dd z,\qquad
F=\begin{pmatrix}1&2\bar z\\0&-1\end{pmatrix}\dd z\wedge\dd\bar z
$$
となる。ここで $\bar\partial A$ では $\dd\bar z\wedge\dd z=-\dd z\wedge\dd\bar z$ の符号を使った。跡は零だが曲率自体は零でない。行列式束が平坦でも元のベクトル束の接続が平坦とは限らない。

曲率の行列がこの正則枠で通常の共役転置について斜 Hermite に見えないことは、枠が正規直交でないためである。計量との両立は実接ベクトルで評価した曲率について $F^*H+HF=0$ と書く。ここで $\dd z\wedge\dd\bar z$ を共役すると符号が反転するため、上の係数行列 $C$ には $C^*H=HC$ が必要である。積を取ると両辺は
$$
\begin{pmatrix}1&\bar z\\z&|z|^2-1\end{pmatrix}
$$
に一致する。正則枠での計算とユニタリ枠での計算を比較するときは、どの内積で随伴を取ったかを確認する。

### テンソル積・双対・行列式の接続

$E,G$ の Chern 接続を $D_E,D_G$ とする。テンソル積上で $D(s\otimes t)=D_Es\otimes t+s\otimes D_Gt$ とすると積計量に両立し、$(0,1)$ 部分が正則構造の $\bar\partial$ に等しい。従って一意性により Chern 接続となる。二回作用させると、1 形式を越えるときの符号で交差項が消え、
$$
F_{E\otimes G}=F_E\otimes I+I\otimes F_G
$$
となる。双対では評価の微分が平坦である条件から曲率は $-{}^\top F_E$ になる。行列式では各枠の外積に $D_E$ を作用させ、対角項だけ残るので曲率は $\tr F_E$ である。

:::example[直和に線束を掛ける]
$\mathcal O(2)\oplus\mathcal O(-1)$ に $\mathcal O(3)$ を掛ける。$\omega_{\rm FS}=iF_{\mathcal O(1)}$ に対し元の $iF$ の係数行列は $\diag(2,-1)$ である。掛けた後は $\diag(5,2)$ となる。行列式の係数は $7$ であり、元の $1$ に階数 $2$ のそれぞれへ加えた $3$ を足した値である。
:::

### 接続行列を三種類の正則枠で計算する

:::example[非対角計量の Chern 接続を求める]
正則な枠で計量行列を次のように置く。計算では $\bar z$ を $b$ と書き、$z,b$ を独立に微分する。点 $z=1$ で曲率の係数を求め、計量との両立を確認する。

$$
\begin{aligned}
H&=\left[\begin{matrix}1 & b^{2}\\z^{2} & b^{2} z^{2} + 1\end{matrix}\right] \\
\det H&=1
\end{aligned}
$$

曲率の $(1,1)$ 成分を $C\,\mathrm dz\wedge\mathrm d\bar z$ と書く。計量適合性から正則枠での接続形式は $A\,\mathrm dz=H^{-1}\partial H$ である。$\bar\partial$ の後では外積を入れ替えるため、係数は $-\partial_b A$ になる。

$$
\begin{aligned}
\partial_z H&=\left[\begin{matrix}0 & 0\\2 z & 2 b^{2} z\end{matrix}\right] \\
A&=\left[\begin{matrix}- 2 b^{2} z & - 2 b^{4} z\\2 z & 2 b^{2} z\end{matrix}\right]
\end{aligned}
$$

$$
\begin{aligned}
C&=\left[\begin{matrix}4 b z & 8 b^{3} z\\0 & - 4 b z\end{matrix}\right] \\
C(p)&=\left[\begin{matrix}4 & 8\\0 & -4\end{matrix}\right],\qquad \tr C(p)=0
\end{aligned}
$$

実方向での曲率の斜 Hermite 性は、この枠の $(1,1)$ 係数について $C(p)^*H(p)=H(p)C(p)$ となる。実際、両辺は $\left[\begin{matrix}4 & 4\\4 & 0\end{matrix}\right]$ である。跡についても $\tr C=-\partial_z\partial_b\log\det H$ で確かめられる。非対角項の微分を保持すると、跡だけからは得られない曲率が見える。
:::

:::example[非対角計量の Chern 接続を求める]
正則な枠で計量行列を次のように置く。計算では $\bar z$ を $b$ と書き、$z,b$ を独立に微分する。点 $z=0$ で曲率の係数を求め、計量との両立を確認する。

$$
\begin{aligned}
H&=\left[\begin{matrix}b z + 1 & b\\z & b z + 2\end{matrix}\right] \\
\det H&=b^{2} z^{2} + 2 b z + 2
\end{aligned}
$$

曲率の $(1,1)$ 成分を $C\,\mathrm dz\wedge\mathrm d\bar z$ と書く。計量適合性から正則枠での接続形式は $A\,\mathrm dz=H^{-1}\partial H$ である。$\bar\partial$ の後では外積を入れ替えるため、係数は $-\partial_b A$ になる。

$$
\begin{aligned}
\partial_z H&=\left[\begin{matrix}b & 0\\1 & b\end{matrix}\right] \\
A&=\left[\begin{matrix}\frac{b \left(b z + 1\right)}{b^{2} z^{2} + 2 b z + 2} & - \frac{b^{2}}{b^{2} z^{2} + 2 b z + 2}\\\frac{1}{b^{2} z^{2} + 2 b z + 2} & \frac{b \left(b z + 1\right)}{b^{2} z^{2} + 2 b z + 2}\end{matrix}\right]
\end{aligned}
$$

$$
\begin{aligned}
C&=\left[\begin{matrix}- \frac{b^{2} z^{2} + 4 b z + 2}{\left(b^{2} z^{2} + 2 b z + 2\right)^{2}} & \frac{2 b \left(b z + 2\right)}{\left(b^{2} z^{2} + 2 b z + 2\right)^{2}}\\\frac{2 z \left(b z + 1\right)}{\left(b^{2} z^{2} + 2 b z + 2\right)^{2}} & - \frac{b^{2} z^{2} + 4 b z + 2}{\left(b^{2} z^{2} + 2 b z + 2\right)^{2}}\end{matrix}\right] \\
C(p)&=\left[\begin{matrix}- \frac{1}{2} & 0\\0 & - \frac{1}{2}\end{matrix}\right],\qquad \tr C(p)=-1
\end{aligned}
$$

実方向での曲率の斜 Hermite 性は、この枠の $(1,1)$ 係数について $C(p)^*H(p)=H(p)C(p)$ となる。実際、両辺は $\left[\begin{matrix}- \frac{1}{2} & 0\\0 & -1\end{matrix}\right]$ である。跡についても $\tr C=-\partial_z\partial_b\log\det H$ で確かめられる。原点で $H(0)$ が対角でも、一次偏微分には非対角成分が残る。
:::

:::example[非対角計量の Chern 接続を求める]
正則な枠で計量行列を次のように置く。計算では $\bar z$ を $b$ と書き、$z,b$ を独立に微分する。点 $z=1$ で曲率の係数を求め、計量との両立を確認する。

$$
\begin{aligned}
H&=\left[\begin{matrix}1 & b\\z & 2 b z + 1\end{matrix}\right] \\
\det H&=b z + 1
\end{aligned}
$$

曲率の $(1,1)$ 成分を $C\,\mathrm dz\wedge\mathrm d\bar z$ と書く。計量適合性から正則枠での接続形式は $A\,\mathrm dz=H^{-1}\partial H$ である。$\bar\partial$ の後では外積を入れ替えるため、係数は $-\partial_b A$ になる。

$$
\begin{aligned}
\partial_z H&=\left[\begin{matrix}0 & 0\\1 & 2 b\end{matrix}\right] \\
A&=\left[\begin{matrix}- \frac{b}{b z + 1} & - \frac{2 b^{2}}{b z + 1}\\\frac{1}{b z + 1} & \frac{2 b}{b z + 1}\end{matrix}\right]
\end{aligned}
$$

$$
\begin{aligned}
C&=\left[\begin{matrix}\frac{1}{\left(b z + 1\right)^{2}} & \frac{2 b \left(b z + 2\right)}{\left(b z + 1\right)^{2}}\\\frac{z}{\left(b z + 1\right)^{2}} & - \frac{2}{\left(b z + 1\right)^{2}}\end{matrix}\right] \\
C(p)&=\left[\begin{matrix}\frac{1}{4} & \frac{3}{2}\\\frac{1}{4} & - \frac{1}{2}\end{matrix}\right],\qquad \tr C(p)=- \frac{1}{4}
\end{aligned}
$$

実方向での曲率の斜 Hermite 性は、この枠の $(1,1)$ 係数について $C(p)^*H(p)=H(p)C(p)$ となる。実際、両辺は $\left[\begin{matrix}\frac{1}{2} & 1\\1 & 0\end{matrix}\right]$ である。跡についても $\tr C=-\partial_z\partial_b\log\det H$ で確かめられる。逆行列を微分するため、行列式の変化も曲率に寄与する。
:::

## 曲率と第一 Chern 類

直線束で $h=\|e\|^2=e^{-\varphi}$ と書くと、
$$
A=\partial\log h=-\partial\varphi,\qquad
F=\bar\partial\partial\log h=\partial\bar\partial\varphi.
\label{eq-cg-line-curvature}
$$
枠を $e'=ge$ に変えると $\varphi'=\varphi-\log|g|^2$ であり、局所対数があるので $\partial\bar\partial\log|g|^2=0$ となる。

:::proposition[第一 Chern 形式]{#prop-cg-first-chern-form}
$E$ の第一 Chern 類の複素 de Rham 表示は $\frac{i}{2\pi}\tr F_\nabla$ である。この形式は実かつ閉であり、計量によらず同じ類を表す。
:::

:::proof
直線束の枠 $e_i$ に対し $e_j=e_i g_{ij}$、$g_{ij}=\exp(2\pi i u_{ij})$ と取ると、$A_j-A_i=2\pi i\,\dd u_{ij}$ である。$a_i=iA_i/(2\pi)$ とすると $\delta a=-\dd u$、$\dd a_i=iF/(2\pi)$ となる。Čech–de Rham 二重複体の全微分を $D=\delta+(-1)^p\dd$（$p$ は Čech 次数）とする。$D(a-u)$ の成分を取れば $iF/(2\pi)-\delta u$ である。従って曲率と :ref[prop-cg-exp-sequence] の整数コサイクルが同じ類を表す。計量との両立により単位枠で $A$ は純虚数値だから $iF$ は実であり、$\dd^2=0$ から閉である。

一般の束で行列式の枠に誘導する接続行列は $\tr A$ である。$\tr(A\wedge A)=0$ は $A_{ij}\wedge A_{ji}$ と逆順の項が相殺するためである。従ってその曲率は $\tr F_\nabla$、第一 Chern 類は $\det E$ の第一 Chern 類に等しい。二計量の違いに由来する接続の差は大域 $1$ 形式であり、跡を取った曲率の差はその外微分なので類は変わらない。
:::

:::example[曲率を積分して次数を出す]
$\mathcal O_{\mathbb P^1}(1)$ の計量を $h=(1+|z|^2)^{-1}$ とすると $iF=i\,\dd z\wedge\dd\bar z/(1+|z|^2)^2$ である。$i\dd z\wedge\dd\bar z=2r\dd r\wedge\dd\theta$ なので
$$
\int_{\mathbb P^1}\frac{iF}{2\pi}=\int_0^\infty\frac{2r}{(1+r^2)^2}\dd r=1.
$$
双対束の計量は逆数なので積分は $-1$ になる。
:::

### 曲線上で零点の個数と一致することを確かめる

コンパクト連結 Riemann 面 $X$ 上の直線束 $L$ を取り、その非零有理型切断 $s$ の零点と極の小円板を取り除けば、$s$ 自身を正則枠にできる。$s=z^m u(z)e$、$u(0)\ne0$ の近くでは、接続形式は $A_s=A_e+m\dd z/z+\dd u/u$ である。従って小円の正向き積分は極限で $2\pi i m$ になる。穴を除いた面で $F=\dd A_s$ なので、Stokes の定理で境界の向きに注意して
$$
\int_X\frac{iF}{2\pi}
=\frac{i}{2\pi}\lim_{\epsilon\downarrow0}
 \left(-\sum_p\int_{|z|=\epsilon}A_s\right)
=\sum_p \operatorname{ord}_p s
$$
を得る。よって Chern 形式の積分は因子の次数に一致する。この証明は、計量を使った滑らかな曲率の積分が、特定の切断の零と極を数える理由を示す。

:::example[零と極を別に持つ切断]
球面上の $L=\mathcal O(2[\infty])$ の有理型切断として $s=1/(z-1)$ を取る。関数としては $z=1$ に単純極、無限遠に単純零点がある。束の無限遠の枠は $w^{-2}$ だから切断の係数は $w^2\cdot w/(1-w)=w^3/(1-w)$ であり、無限遠での切断の位数は $3$ になる。従って積分は $3-1=2$ である。
:::

### 曲率と零因子を結ぶ分布の公式

$h$ を $L$ の Hermite 計量とし、$s$ の零因子が滑らかな超曲面 $D_j$ の和である場合を考える。多重度を $m_j$ とし、各点の近くで $s=(u\prod_j f_j^{m_j})e$、$u$ は正則単元、$\mathrm df_j$ は $D_j$ 上で非零と書けるとする。滑らかな成分の積分を $[D_j]$ で表すと、分布の意味で
$$
\frac{i}{2\pi}\partial\bar\partial\log\|s\|_h^2
=\sum_j m_j[D_j]-c_1(L,h)
$$
となる。$s$ の係数に極を許す場合も、同じ局所形で負の $m_j$ を許せば成り立つ。

一変数で $\Delta_{\rm eucl}\log|z|^2=4\pi\delta_0$ である。実際、滑らかな台付き試験関数を掛け、$|z|>\varepsilon$ で Green の公式を使う。内側境界の法微分は $-2/\varepsilon$ であり、境界長の因子 $\varepsilon$ と試験関数の円周平均から $4\pi$ 倍の原点評価が残る。試験関数の法微分を含む項は $O(\varepsilon|\log\varepsilon|)$ で消える。さらに $i\partial\bar\partial v=(\Delta_{\rm eucl}v/2)\mathrm dx\wedge\mathrm dy$ だから、$i\partial\bar\partial\log|z|^2/(2\pi)$ が原点の積分を与える。滑らかな $D_j$ では座標の一つを $f_j$ に取り、残りの変数で積分すると同じ式を各法方向に適用できる。$\log|u|^2$ の二重微分は零で、$\log h$ の二重微分は $-F_L$ である。これらを足すと上の式になる。

### 計量を変えるときの曲率の差

同じ直線束の二つの計量を $h$、$h'=e^{-u}h$ とする。枠変換では両計量に同じ $|g|^2$ が掛かるので $u$ は大域実関数である。接続形式の差は $-\partial u$ であり、曲率の差は $\partial\bar\partial u$ である。従って $iF$ の変化は完全形式であり、コンパクト曲線では積分が同じである。一方、局所での正負は変わりうるので、次数が正という条件と与えられた計量の曲率が正という条件は区別する。

:::example[正の次数でも計量によって負の点が出る]
$\mathbb P^1$ で $h=(1+|z|^2)^{-1}$ とし、$u=-2|z|^2/(1+|z|^2)$ とする。$u$ は無限遠座標 $w$ では $-2/(1+|w|^2)$ なので滑らかである。$0$ で $u_{z\bar z}=-2$、元の曲率係数は $1$ だから、新しい曲率係数は $-1$ である。積分はなお $2\pi$ であり、全体で負になったわけではない。
:::

### 第一 Chern 形式を跡で照合する

:::example[非対角計量の Chern 接続を求める]
正則な枠で計量行列を次のように置く。計算では $\bar z$ を $b$ と書き、$z,b$ を独立に微分する。点 $z=0$ で曲率の係数を求め、計量との両立を確認する。

$$
\begin{aligned}
H&=\left[\begin{matrix}b z + 2 & b\\z & 2\end{matrix}\right] \\
\det H&=b z + 4
\end{aligned}
$$

曲率の $(1,1)$ 成分を $C\,\mathrm dz\wedge\mathrm d\bar z$ と書く。計量適合性から正則枠での接続形式は $A\,\mathrm dz=H^{-1}\partial H$ である。$\bar\partial$ の後では外積を入れ替えるため、係数は $-\partial_b A$ になる。

$$
\begin{aligned}
\partial_z H&=\left[\begin{matrix}b & 0\\1 & 0\end{matrix}\right] \\
A&=\left[\begin{matrix}\frac{b}{b z + 4} & 0\\\frac{2}{b z + 4} & 0\end{matrix}\right]
\end{aligned}
$$

$$
\begin{aligned}
C&=\left[\begin{matrix}- \frac{4}{\left(b z + 4\right)^{2}} & 0\\\frac{2 z}{\left(b z + 4\right)^{2}} & 0\end{matrix}\right] \\
C(p)&=\left[\begin{matrix}- \frac{1}{4} & 0\\0 & 0\end{matrix}\right],\qquad \tr C(p)=- \frac{1}{4}
\end{aligned}
$$

実方向での曲率の斜 Hermite 性は、この枠の $(1,1)$ 係数について $C(p)^*H(p)=H(p)C(p)$ となる。実際、両辺は $\left[\begin{matrix}- \frac{1}{2} & 0\\0 & 0\end{matrix}\right]$ である。跡についても $\tr C=-\partial_z\partial_b\log\det H$ で確かめられる。同じ跡は行列式束の計量 $\det H$ からも求まる。局所の二次形式の値と大域の第一類は別の対象である。
:::

:::example[非対角計量の Chern 接続を求める]
正則な枠で計量行列を次のように置く。計算では $\bar z$ を $b$ と書き、$z,b$ を独立に微分する。点 $z=1$ で曲率の係数を求め、計量との両立を確認する。

$$
\begin{aligned}
H&=\left[\begin{matrix}2 & b\\z & b^{2} z^{2} + 2\end{matrix}\right] \\
\det H&=2 b^{2} z^{2} - b z + 4
\end{aligned}
$$

曲率の $(1,1)$ 成分を $C\,\mathrm dz\wedge\mathrm d\bar z$ と書く。計量適合性から正則枠での接続形式は $A\,\mathrm dz=H^{-1}\partial H$ である。$\bar\partial$ の後では外積を入れ替えるため、係数は $-\partial_b A$ になる。

$$
\begin{aligned}
\partial_z H&=\left[\begin{matrix}0 & 0\\1 & 2 b^{2} z\end{matrix}\right] \\
A&=\left[\begin{matrix}- \frac{b}{2 b^{2} z^{2} - b z + 4} & - \frac{2 b^{3} z}{2 b^{2} z^{2} - b z + 4}\\\frac{2}{2 b^{2} z^{2} - b z + 4} & \frac{4 b^{2} z}{2 b^{2} z^{2} - b z + 4}\end{matrix}\right]
\end{aligned}
$$

$$
\begin{aligned}
C&=\left[\begin{matrix}- \frac{2 \left(b^{2} z^{2} - 2\right)}{\left(2 b^{2} z^{2} - b z + 4\right)^{2}} & \frac{4 b^{2} z \left(b^{2} z^{2} - b z + 6\right)}{\left(2 b^{2} z^{2} - b z + 4\right)^{2}}\\\frac{2 z \left(4 b z - 1\right)}{\left(2 b^{2} z^{2} - b z + 4\right)^{2}} & \frac{4 b z \left(b z - 8\right)}{\left(2 b^{2} z^{2} - b z + 4\right)^{2}}\end{matrix}\right] \\
C(p)&=\left[\begin{matrix}\frac{2}{25} & \frac{24}{25}\\\frac{6}{25} & - \frac{28}{25}\end{matrix}\right],\qquad \tr C(p)=- \frac{26}{25}
\end{aligned}
$$

実方向での曲率の斜 Hermite 性は、この枠の $(1,1)$ 係数について $C(p)^*H(p)=H(p)C(p)$ となる。実際、両辺は $\left[\begin{matrix}\frac{2}{5} & \frac{4}{5}\\\frac{4}{5} & - \frac{12}{5}\end{matrix}\right]$ である。跡についても $\tr C=-\partial_z\partial_b\log\det H$ で確かめられる。原点から外れた点では、行列式の一次微分どうしの積も対数の二階微分に入る。
:::

### 正則部分束と商束の曲率

:::proposition[第二基本形式による曲率の分解]{#prop-cg-subbundle-curvature}
$S\subset E$ を正則部分束、$Q=E/S$ とし、計量を制限と商で入れる。滑らかな直交分解 $E=S\oplus Q$ で Chern 接続の $S$ から $Q$ への成分を $B$ とすると、$B$ は $(1,0)$ 型であり
$$
F_S=(F_E)_{SS}+B^*\wedge B,\qquad
F_Q=(F_E)_{QQ}+B\wedge B^*
$$
である。ここで $B^*$ は共役転置により型を $(0,1)$ に変えたものをいう。
:::

:::proof
$S$ は正則なので $D''s=\bar\partial_Es$ は $S$ に値を持ち、$Q$ への射影は零である。よって $B$ は $(1,0)$ 型である。計量適合性を直交枠で書くと、逆向きの成分は $-B^*$ となる。対角成分は制限・商計量と両立し、$(0,1)$ 成分がそれぞれの正則構造なので Chern 接続の一意性から $D_S,D_Q$ に等しい。従って接続のブロック表示は
$$
D=\begin{pmatrix}D_S&-B^*\\ B&D_Q\end{pmatrix}.
$$
二回作用させたとき、最初の対角の成分は $F_S-B^*\wedge B$、第二の対角の成分は $F_Q-B\wedge B^*$ である。これを移項して主張を得る。
:::

:::example[平坦束から正の商曲率が生まれる]
$\mathbb P^1$ で $S=\mathcal O(-1)$ を $(1,z)$ が張る標準部分束として $\mathcal O^{\oplus2}$ に含める。標準平坦計量を制限した $S$ の計量は $1+|z|^2$ である。従って $F_S=-(1+|z|^2)^{-2}\mathrm dz\wedge\mathrm d\bar z$ である。商の曲率の跡は逆符号になる。例えば $z=1$ では $S$ の係数は $-1/4$、商の係数は $1/4$ である。元の束の曲率が零でも、部分束と商束に非零な曲率が生じる。
:::

## Chern 類の微分形式による表示

曲率は行列値の $2$ 形式である。偶数次数の成分は外積の順序を交換できるので行列式を定義できる。
$$
c(E,\nabla)=\det\bigl(I+\tfrac{i}{2\pi}F_\nabla\bigr)=1+c_1+c_2+\cdots.
$$
枠変換で $F$ は共役されるため形式は大域的である。Chern–Weil の同定により、次数 $2k$ 成分は $k$ 番目の Chern 類を表す。:ref[thm-dtop-chern-weil-class] を参照する。ここでは閉性と計量独立性を具体的に確認する。

:::proposition[閉性と推移公式]{#prop-cg-curvature-transgression}
不変斉次多項式 $P$ の次数を $k$ とし、同じ記号でその対称多重線形化を表す。接続 $A_t$ で曲率を $F_t$ とすると $P(F_t,\ldots,F_t)$ は閉であり、
$$
\frac{\dd}{\dd t}P(F_t,\ldots,F_t)
=k\,\dd P(\dot A_t,F_t,\ldots,F_t)
$$
である。
:::

:::proof
$\dd F+A\wedge F-F\wedge A=0$ は $F=\dd A+A^2$ を展開すると従う。これを Bianchi 恒等式という。部分ごとに微分した式で、フレームの接続の寄与は不変性を微分した $\sum_jP(B_1,\ldots,[A,B_j],\ldots)=0$ により、外積の次数符号を付けて相殺する。従って $P$ の外微分は、各引数に共変外微分を施した和である。曲率だけの項では全て Bianchi 恒等式で零になる。$A_t$ を微分すると $\dot F_t=\dd\dot A_t+A_t\wedge\dot A_t+\dot A_t\wedge A_t$ となり、これは $\dot A_t$ の共変外微分である。よって $\partial_t P$ に代入し、残りの引数の共変微分が消えることを使えば式を得る。接続のアフィン補間を $t$ で積分すると、二つの接続から得られる $2k$ 形式の差は完全である。
:::

$B=iF/(2\pi)$ とすると行列式の二次項は $c_2=\frac12((\tr B)^2-\tr B^2)$ である。実際、二つの異なる行、列を選ぶ主小行列式を足すとこの差の半分になる。直和の接続を選ぶと行列がブロック対角になるから $c(E\oplus E')=c(E)c(E')$ である。また、束の短完全列は滑らかな計量で直交補を取り分裂できるので、同じ積公式が類について成り立つ。正則な分裂は必要ない。

:::example[直線二本の直和]
$c_1(L)=H$ のとき、$E=L\oplus L^{\otimes2}$ は $c_1(E)=3H,c_2(E)=2H^2$ を持つ。局所曲率行列に $B,2B$ が現れるので行列式は $(1+B)(1+2B)=1+3B+2B^2$ になる。
:::

### Chern 指標との関係

$B=iF/(2\pi)$ とすると、Chern 指標の形式は $\tr(e^B)$ である。指数級数の係数を読めば $0$ 次は $\rank E$、二次は $c_1$、四次は $\frac12(c_1^2-2c_2)$ である。曲率の不変多項式であるため :ref[prop-cg-curvature-transgression] によって各項は閉じ、類は接続によらない。直和では跡が加算され、テンソル積では曲率が $F_E\otimes I+I\otimes F_{E'}$ になる。この二項は可換なので指数を分けられ、跡を取ると Chern 指標は積になる。従って Chern 類は直和に対して積、Chern 指標は直和に対して和になる。Riemann–Roch で Chern 指標を用いると Euler 標数の加法性と両立する。

:::example[四次成分だけを比べる]
$E=L\oplus L^2$ で $c_1(L)=H$ なら $\operatorname{ch}(E)=e^H+e^{2H}$ であり、四次成分は $(1+4)H^2/2=5H^2/2$ である。一方 $c_1=3H,c_2=2H^2$ から $\frac12(c_1^2-2c_2)=\frac12(9-4)H^2$ となり一致する。
:::

### 分解する束で Chern 類を読み取る

直和接続の曲率はブロック対角になる。従って $\det(I+iF/(2\pi))$ はブロックごとの積に分かれ、直和に対する全 Chern 類の積公式が得られる。正則短完全列があっても Hermite 計量による直交分解は滑らかにできるため、接続の独立性から同じコホモロジーの公式が成り立つ。局所の形式そのものを比べると、直交分解が接続に保たれない場合には完全形式の差がありうる。

:::example[跡が零でも第二類が残る]
$\mathbb P^2$ 上で $E=\mathcal O(2)\oplus\mathcal O(-2)$ とする。$H=c_1(\mathcal O(1))$ に対し $c(E)=(1+2H)(1-2H)=1-4H^2$ である。従って $c_1=0$ だが $\int c_2=-4$ である。一方 $\operatorname{ch}_2=(c_1^2-2c_2)/2=4H^2$ となる。第二 Chern 類と Chern 指標の二次項は同じものではない。
:::

### 全 Chern 類の係数で差を確認する

:::example[完全交差の標準束と Chern 数を計算する]
$X\subset\mathbb P^{4}$ を次数 $4$ の横断的な滑らかな完全交差とし、$H=c_1(\mathcal O_X(1))$ とする。法束は $\bigoplus_a\mathcal O_X(d_a)$ であり、次元は $3$ である。

接束の完全列から $c(T_X)c(N_{X/\mathbb P^N})=c(T_{\mathbb P^N}|_X)$ を得る。Euler 列により右辺は $(1+H)^{N+1}$ である。$c(T_X)=1+c_1H+\cdots+c_nH^n$ と書いて係数を比べる。

$$
\begin{aligned}
c(N)&=4 H + 1 \\
c_{1}+(4)&=\binom{5}{1},\qquad c_{1}=1 \\
c_{2}+(4)&=\binom{5}{2},\qquad c_{2}=6 \\
c_{3}+(24)&=\binom{5}{3},\qquad c_{3}=-14 \\
c(T_X)&=- 14 H^{3} + 6 H^{2} + H + 1
\end{aligned}
$$

法束の行列式を取った随伴公式は $K_X=\mathcal O_X(-1)$ である。これは $c(T_X)$ の一次の係数の符号を反転した値にも等しい。さらに $X$ を一般の $3$ 個の超平面で切り、交点数を積で数えると $\int_XH^3=4$ となる。

$$
\begin{aligned}
\int_X c_n(T_X)&=-14\cdot4=-56
\end{aligned}
$$

最高 Chern 類は実接束の Euler 類であるため、ここで得た整数は位相的 Euler 標数にも等しい。同じ次元でも各 Chern 次数で符号が異なる。$\operatorname{ch}_2=(c_1^2-2c_2)/2$ には独自の係数が入る。
:::

:::example[完全交差の標準束と Chern 数を計算する]
$X\subset\mathbb P^{5}$ を次数 $3,3$ の横断的な滑らかな完全交差とし、$H=c_1(\mathcal O_X(1))$ とする。法束は $\bigoplus_a\mathcal O_X(d_a)$ であり、次元は $3$ である。

接束の完全列から $c(T_X)c(N_{X/\mathbb P^N})=c(T_{\mathbb P^N}|_X)$ を得る。Euler 列により右辺は $(1+H)^{N+1}$ である。$c(T_X)=1+c_1H+\cdots+c_nH^n$ と書いて係数を比べる。

$$
\begin{aligned}
c(N)&=\left(3 H + 1\right)^{2} \\
c_{1}+(6)&=\binom{6}{1},\qquad c_{1}=0 \\
c_{2}+(9)&=\binom{6}{2},\qquad c_{2}=6 \\
c_{3}+(36)&=\binom{6}{3},\qquad c_{3}=-16 \\
c(T_X)&=- 16 H^{3} + 6 H^{2} + 1
\end{aligned}
$$

法束の行列式を取った随伴公式は $K_X=\mathcal O_X(0)$ である。これは $c(T_X)$ の一次の係数の符号を反転した値にも等しい。さらに $X$ を一般の $3$ 個の超平面で切り、交点数を積で数えると $\int_XH^3=9$ となる。

$$
\begin{aligned}
\int_X c_n(T_X)&=-16\cdot9=-144
\end{aligned}
$$

最高 Chern 類は実接束の Euler 類であるため、ここで得た整数は位相的 Euler 標数にも等しい。接束の全 Chern 類は法束の全 Chern 類を逆数にして掛ける。そのため法束の二次項を無視することはできない。
:::

:::example[完全交差の標準束と Chern 数を計算する]
$X\subset\mathbb P^{6}$ を次数 $2,4$ の横断的な滑らかな完全交差とし、$H=c_1(\mathcal O_X(1))$ とする。法束は $\bigoplus_a\mathcal O_X(d_a)$ であり、次元は $4$ である。

接束の完全列から $c(T_X)c(N_{X/\mathbb P^N})=c(T_{\mathbb P^N}|_X)$ を得る。Euler 列により右辺は $(1+H)^{N+1}$ である。$c(T_X)=1+c_1H+\cdots+c_nH^n$ と書いて係数を比べる。

$$
\begin{aligned}
c(N)&=\left(2 H + 1\right) \left(4 H + 1\right) \\
c_{1}+(6)&=\binom{7}{1},\qquad c_{1}=1 \\
c_{2}+(14)&=\binom{7}{2},\qquad c_{2}=7 \\
c_{3}+(50)&=\binom{7}{3},\qquad c_{3}=-15 \\
c_{4}+(-34)&=\binom{7}{4},\qquad c_{4}=69 \\
c(T_X)&=69 H^{4} - 15 H^{3} + 7 H^{2} + H + 1
\end{aligned}
$$

法束の行列式を取った随伴公式は $K_X=\mathcal O_X(-1)$ である。これは $c(T_X)$ の一次の係数の符号を反転した値にも等しい。さらに $X$ を一般の $4$ 個の超平面で切り、交点数を積で数えると $\int_XH^4=8$ となる。

$$
\begin{aligned}
\int_X c_n(T_X)&=69\cdot8=552
\end{aligned}
$$

最高 Chern 類は実接束の Euler 類であるため、ここで得た整数は位相的 Euler 標数にも等しい。曲面までにしかない項と、四次元まで残る項を同じ漸化式で比較できる。
:::

### 二次の Chern 形式の規格化

$B=iF/(2\pi)$ とすると、$B$ の成分の形式次数は二であり、互いに交換する。従って行列式の二次の係数は
$$
c_2(E,h)=\sum_{j<k}(B_{jj}B_{kk}-B_{jk}B_{kj})
=\tfrac12\bigl((\tr B)^2-\tr(B^2)\bigr)
$$
である。実際、$\tr(B^2)$ の和には対角の二乗と $j\ne k$ の積が含まれる。前者は $(\tr B)^2$ との差で相殺し、後者は各対が二回現れるため $1/2$ が付く。Chern 指標の二次項は $\tr(B^2)/2$ なので、同じ入力から得ても一致するとは限らない。

:::example[三本の直線束の二次類を二通りで数える]
$X=\mathbb P^1\times\mathbb P^1$ とし、規格化した超平面類を $a,b$ とする。$E=\mathcal O(1,2)\oplus\mathcal O(2,-1)\oplus\mathcal O(-3,1)$ を取る。三つの第一類の和は $2b$ で、二乗は $0$ である。各第一類の二乗は $4ab,-4ab,-6ab$ なので、$\tr(B^2)=-6ab$ である。従って $c_2=3ab$ で、Chern 指標の二次項は $-3ab$ となる。積の方法でも、二本ごとの積は $3ab,-5ab,5ab$ であり、その和が $3ab$ になる。
:::

## 正直線束

:::definition[正直線束]{#def-cg-positive-line}
正則直線束が :term[正直線束]{en="positive line bundle" yomi="せいちょくせんそく"} であるとは、$iF$ が正定値の実 $(1,1)$ 形式となる Hermite 計量が存在することである。$iF=i\sum Q_{j\bar k}\dd z_j\wedge\dd\bar z_k$ に対し、正定値性は Hermite 行列 $Q$ の正定値性をいう。
:::

:eqref[eq-cg-line-curvature] により局所重み $\varphi$ の複素 Hessian が $Q$ である。局所枠の変更は $\varphi$ に局所正則関数の実部を加えるだけなので、その Hessian を変えない。正直線束のテンソル積は重みを足すため正である。双対束では符号が反転する。

:::example[超平面束は正である]
$\mathcal O(-1)$ の計量は $\C^{n+1}$ の標準内積を制限したもので、枠 $Z/Z_0$ で $1+|z|^2$ になる。双対 $\mathcal O(1)$ の重みは $\varphi=\log(1+|z|^2)$ であり、
$$
Q_{j\bar k}=\frac{(1+|z|^2)\delta_{jk}-\bar z_j z_k}{(1+|z|^2)^2}.
$$
Cauchy–Schwarz により $Q(v,\bar v)\ge |v|^2/(1+|z|^2)^2>0$ である。従って $\mathcal O(1)$ は正である。
:::

:::example[自明束が正になれない場合]
コンパクト連結な正次元の多様体上で、自明束の重みは大域実関数である。その最大点では各複素直線方向で実 Hessian の円周方向の跡が非正なので、複素 Hessian は正定値になれない。従って局所的に $e^{-|z|^2}$ を選べることから、大域的な正値性は従わない。
:::

### 正値性から曲線の次数を制限する

正直線束 $L$ とコンパクト連結 Riemann 面からの非定数正則写像 $f:C\to X$ を取る。引き戻し束の曲率 $f^*(iF_L)$ は各点で非負であり、接写像が非零なところでは正である。接写像が全点で零なら、局所座標成分の導関数が全て零で $f$ は定数になる。従って曲率積分は正となり、$\deg f^*L>0$ である。例えば二次 Veronese 曲線 $f:\mathbb P^1\to\mathbb P^2$ では、$f^*\mathcal O(1)=\mathcal O(2)$ なので積分は正確に $2$ である。

一方、$\dim_\C X\ge2$ で $\pi:\widetilde X\to X$ が一点のブローアップなら、$\pi^*L$ は例外射影空間の各曲線上で次数 $0$ になる。これらの曲線では $\pi$ が定数なので引き戻し束は自明である。従って $\widetilde X$ 上の $\pi^*L$ には、別の計量を選んでも全方向に正な曲率は入らない。これは埋め込み定理の証明で例外因子の方向を補う直線束が必要になる理由である。

### 積上の曲率と正定値性

$X,Y$ 上の正直線束 $L,M$ を取り、積束 $p_X^*L\otimes p_Y^*M$ に積計量を入れる。その曲率は二つの引き戻しの和である。非零な接ベクトル $(v,w)$ に対して、少なくとも一方の成分が非零なのでこの和は正になる。一方 $p_X^*L$ だけの曲率は $Y$ 方向を消す。曲率の正値性を使う際は、引き戻しの微分の核を確認する必要がある。

:::example[射影直線の積で三種類の符号]
$\mathbb P^1\times\mathbb P^1$ の $\mathcal O(a,b)$ に積計量を使う。$iF=a\omega_1+b\omega_2$ となるので $(a,b)=(2,3)$ なら正、$(2,0)$ なら半正で零方向があり、$(2,-1)$ なら不定符号である。最後の束が他の計量で正にならないことは、第二因子の曲線に制限すると次数 $-1$ になることから分かる。
:::

### 正の曲率が与える固有値を局所で読む

:::example[局所重みの複素 Hessian で曲率を調べる]
$\varphi=\bar z_1 \bar z_2 z_{1} z_{2} + 2 \bar z_1 z_{1} + 3 \bar z_2 z_{2}$ とし、点 $p=(0,1)$ で自明直線束の計量 $h=e^{-\varphi}$ の曲率を求める。$\partial_{z_j}$ は $\bar z_k$ を固定して微分するため、計算中は両者を独立な変数として扱う。

$$
\begin{aligned}
G=(\varphi_{j\bar k})&=\left[\begin{matrix}\bar z_2 z_{2} + 2 & \bar z_1 z_{2}\\\bar z_2 z_{1} & \bar z_1 z_{1} + 3\end{matrix}\right] \\
\det G&=2 \bar z_1 z_{1} + 3 \bar z_2 z_{2} + 6 \\
G(p)&=\left[\begin{matrix}3 & 0\\0 & 3\end{matrix}\right],\qquad\det G(p)=9
\end{aligned}
$$

点での先頭小行列式 $G_{1\bar1}(p)=3$ と行列式がともに正なので、$G(p)$ は正定値である。正定値性は開条件であるため、近傍でも $iF=i\sum G_{j\bar k}\mathrm dz_j\wedge\mathrm d\bar z_k$ が正となる。係数の行列式だけでは、曲率が正定値なのか負定値なのかは分からない。この二次元の例では先頭小行列式も確認することで符号を決めた。
:::

:::example[局所重みの複素 Hessian で曲率を調べる]
$\varphi=\bar z_1 z_{1} + 3 \bar z_2 z_{2} + \left(\bar z_1 + \bar z_2\right)^{2} \left(z_{1} + z_{2}\right)^{2}$ とし、点 $p=(1,-1)$ で自明直線束の計量 $h=e^{-\varphi}$ の曲率を求める。$\partial_{z_j}$ は $\bar z_k$ を固定して微分するため、計算中は両者を独立な変数として扱う。

$$
\begin{aligned}
G=(\varphi_{j\bar k})&=\left[\begin{matrix}4 \left(\bar z_1 + \bar z_2\right) \left(z_{1} + z_{2}\right) + 1 & 4 \left(\bar z_1 + \bar z_2\right) \left(z_{1} + z_{2}\right)\\4 \left(\bar z_1 + \bar z_2\right) \left(z_{1} + z_{2}\right) & 4 \left(\bar z_1 + \bar z_2\right) \left(z_{1} + z_{2}\right) + 3\end{matrix}\right] \\
\det G&=16 \bar z_1 z_{1} + 16 \bar z_1 z_{2} + 16 \bar z_2 z_{1} + 16 \bar z_2 z_{2} + 3 \\
G(p)&=\left[\begin{matrix}1 & 0\\0 & 3\end{matrix}\right],\qquad\det G(p)=3
\end{aligned}
$$

点での先頭小行列式 $G_{1\bar1}(p)=1$ と行列式がともに正なので、$G(p)$ は正定値である。正定値性は開条件であるため、近傍でも $iF=i\sum G_{j\bar k}\mathrm dz_j\wedge\mathrm d\bar z_k$ が正となる。係数の行列式だけでは、曲率が正定値なのか負定値なのかは分からない。この二次元の例では先頭小行列式も確認することで符号を決めた。
:::

:::example[局所重みの複素 Hessian で曲率を調べる]
$\varphi=\bar z_1^{2} z_{1}^{2} + 3 \bar z_1 z_{1} + \bar z_2 z_{2} + \left(\bar z_1 - 2 \bar z_2\right) \left(z_{1} - 2 z_{2}\right)$ とし、点 $p=(1,0)$ で自明直線束の計量 $h=e^{-\varphi}$ の曲率を求める。$\partial_{z_j}$ は $\bar z_k$ を固定して微分するため、計算中は両者を独立な変数として扱う。

$$
\begin{aligned}
G=(\varphi_{j\bar k})&=\left[\begin{matrix}4 \left(\bar z_1 z_{1} + 1\right) & -2\\-2 & 5\end{matrix}\right] \\
\det G&=4 \left(5 \bar z_1 z_{1} + 4\right) \\
G(p)&=\left[\begin{matrix}8 & -2\\-2 & 5\end{matrix}\right],\qquad\det G(p)=36
\end{aligned}
$$

点での先頭小行列式 $G_{1\bar1}(p)=8$ と行列式がともに正なので、$G(p)$ は正定値である。正定値性は開条件であるため、近傍でも $iF=i\sum G_{j\bar k}\mathrm dz_j\wedge\mathrm d\bar z_k$ が正となる。係数の行列式だけでは、曲率が正定値なのか負定値なのかは分からない。この二次元の例では先頭小行列式も確認することで符号を決めた。
:::

### 曲線では正次数から正計量を得られる

連結コンパクト曲線 $X$ で $\deg L>0$ とする。面積形式 $\omega$ を一つ固定して、$\lambda=2\pi\deg L/\int_X\omega$ とおく。正則直線束の任意の滑らかな Hermite 計量 $h$ に対し、$\lambda\omega-iF_h$ の積分は零である。そこで $i\partial\bar\partial u=\lambda\omega-iF_h$ を曲線の Poisson 方程式として解く。局所共形座標でこの作用素の係数は実 Laplacian に比例し、定数が核を尽くすため、右辺の平均が零であることが可解条件の全てである。楕円型正則性で滑らかな $u$ を取れる。新計量 $e^{-u}h$ の $iF$ は $\lambda\omega>0$ である。逆に正の曲率なら積分が正なので $\deg L>0$ である。

:::example[面積を変えたときの定数]
曲線で $\deg L=3$、$\int_X\omega=12\pi$ なら、正に調整した曲率の縮約定数は $1/2$ である。$\omega$ を二倍すると定数は $1/4$ となるが、曲率形式の積 $\lambda\omega$ は同じである。ここで縮約の定数は基準計量と束の次数の両方に依存する。
:::

## Hermite–Einstein 計量の概観

Kähler 形式 $\omega$ との外積の計量随伴を $\Lambda_\omega$ とする。Hermite–Einstein 条件は
$$
i\Lambda_\omega F_\nabla=\lambda\,\id_E
$$
である。これは曲率の全成分を定数にする条件ではなく、$\omega$ に沿う縮約を各ファイバーでスカラーにする条件である。正規直交枠で $\alpha\wedge\omega^{n-1}/(n-1)!=(\Lambda_\omega\alpha)\omega^n/n!$ は、対角成分だけが残る外積計算から従う。従って $\operatorname{Vol}_\omega(X)=\int\omega^n/n!$、$\deg_\omega E=\int c_1(E)\wedge\omega^{n-1}/(n-1)!$ と置くと
$$
\lambda=\frac{2\pi\deg_\omega E}{\rank(E)\operatorname{Vol}_\omega(X)}
$$
である。$\mu_\omega(E)=\deg_\omega E/\rank(E)$ を :term[傾き]{en=slope yomi=かたむき} という。$0<
ank\mathcal S<
ank E$ の任意の連接部分層 $\mathcal S\subset\mathcal O(E)$ が $\mu_\omega(\mathcal S)<\mu_\omega(E)$ を満たすとき、$E$ は傾き安定である。$\mathcal S$ の次数は行列式直線束 $(\wedge^{
ank\mathcal S}\mathcal S)^{**}$ で測る。これは滑らかな底上では階数一の反射的層が局所自由であることによる。同じ傾きの安定束の直和を傾き多安定という。

連結コンパクト Kähler 多様体では Hermite–Einstein 計量を持つことは傾き多安定と同値である（Donaldson–Uhlenbeck–Yau の定理）。一般の存在証明では計量の非線形偏微分方程式の評価と、極限が退化する場合の部分層の構成が必要である。証明は :cite[huybrechts-cg] を参照する。

:::example[射影直線上の定数]
$\omega=iF_{\mathcal O(1)}$ を $\mathbb P^1$ の Kähler 形式とすると $\Lambda_\omega\omega=1$ である。従って $\mathcal O(k)$ の誘導計量は $iF=k\omega$、Hermite–Einstein 定数は $k$ である。$k<0$ でもこの条件は成立し、正値性は要求されない。
:::

### 曲率の縮約を一定にする条件

直線束では計量を $h'=e^{-u}h$ と変えるだけで $i\Lambda F_{h'}=i\Lambda F_h+\Lambda(i\partial\bar\partial u)$ となる。右辺の最後の作用素は、非負符号の Laplacian を使うと $-\frac12\Delta_{\mathrm d}u$ である。実関数 $v$ を掛けて積分し、閉じた $\omega^{n-1}$ を用いて部分積分すると、その値は $-\int_X i\partial v\wedge\bar\partial u\wedge\omega^{n-1}/(n-1)!$ である。実部を取り、ユニタリ枠で外積を計算すると $-\frac12\int_X g(\mathrm dv,\mathrm du)\,\omega^n/n!$ となる。従って関数への Laplacian の弱形式と一致する。Hodge 理論により、連結コンパクトな底では平均零の滑らかな関数がこの作用素の像である。従って $u$ を適当に選んで縮約を平均値にできる。二つの解の差は調和関数なので定数であり、計量の定数倍が残る自由度である。

### 正則部分束の傾きが大きくならない理由

正則部分束 $S\subset E$ を滑らかな直交分解 $E=S\oplus S^\perp$ で表す。$S$ の正則性から $S$ から $S^\perp$ への接続の成分 $\beta$ は型 $(1,0)$ である。計量との両立で逆向きの成分は $-\beta^*$ となる。接続の上左ブロックは $S$ の Chern 接続であり、曲率の上左ブロックは $F_S-\beta^*\wedge\beta$ である。従って
$$
\tr(i\Lambda F_S)=\tr((i\Lambda F_E)|_S)-|\beta|^2
$$
となる。最後の符号は、$\beta=\sum_j B_j\theta_j$ と書けば $\beta^*\wedge\beta$ の対角型が $-B_j^*B_j\theta_j\wedge\bar\theta_j$ になることから従う。従って $E$ が Hermite–Einstein なら積分して $\mu_\omega(S)\le\mu_\omega(E)$ を得る。等号なら $\beta=0$ なので直交補も接続で保たれ、正則な直和分解になる。この計算は滑らかな正則部分束に対するものであり、一般の連接部分層には特異点を扱う議論が必要である。

### 曲線の直線束の場合に計量を調整する

コンパクト連結曲線で面積形式 $\omega$ をとり、目標曲率を $\lambda\omega$、$\lambda=2\pi\deg L/\int_X\omega$ とする。既存の実曲率 $iF_h$ との差 $\lambda\omega-iF_h$ は積分零である。局所では $i\partial\bar\partial u=(\Delta_{\rm eucl}u/2)\mathrm dA$ となるので、必要な $u$ は共形計量の Poisson 方程式で得られる。右辺の平均が零であるため、閉曲面の Laplace 作用素の核である定数と直交する。楕円型 Hodge 理論の逆作用素で滑らかな平均零解を取り、$h'=e^{-u}h$ とすればよい。従って曲線の直線束にはこの条件を満たす計量があり、二つの計量の対数差は調和だから定数である。

:::example[直和で平均傾きが一致しない]
$\mathbb P^1$ の $\omega=\omega_{\rm FS}$ で面積を $2\pi$ とする。$E=\mathcal O(3)\oplus\mathcal O(1)$ の平均傾きは $2$ だが、第一の部分束は傾き $3$ を持つ。Hermite–Einstein 条件なら部分束の傾きは $2$ 以下でなければならないため、この束には条件を満たす計量がない。対して $\mathcal O(2)\oplus\mathcal O(2)$ では積計量の $i\Lambda F$ が $2 I$ となる。
:::

### 縮約の条件から大域例を設計する

:::example[積の計量で対角曲率の縮約を計算する]
複素次元 $2$ の一点で正則余接枠 $\eta_j$ を取り、$a_j=i\eta_j\wedge\bar\eta_j$ とする。$a_j^2=0$ である。最高次形式 $\beta$ の $a_1\wedge\cdots\wedge a_n$ に関する係数を $[\beta]_{\mathrm{top}}$ と書く。基本形式には次の正係数の線形結合を使う。

$$
\begin{aligned}
\omega&=2 a_{1} + 5 a_{2}
\end{aligned}
$$

外積作用素 $L=\omega\wedge$ の $n-2$ 乗は、$a_j$ の張る実空間を $n-1$ 個の添字の積が張る空間に写す。上位側の基底を昇順の添字の積とし、$a_i$ に $\omega^{n-2}$ を掛けると、同じ添字に二度出会う項が消える。残りの各添字集合の係数には並べ替えの数 $(n-2)!$ が掛かる。

$$
\begin{aligned}
[L^{0}]&=\left[\begin{matrix}1 & 0\\0 & 1\end{matrix}\right] \\
\det[L^{0}]&=1
\end{aligned}
$$

最高次の係数 $Q(a,b)=[a\wedge b\wedge\omega^{n-2}]_{\mathrm{top}}$ は、行と列が異なるときに残りの添字の重みの積から求まる。$Q(a,\omega)=0$ となる空間を $P$ とし、基底を列に並べた行列を $N$ とする。

$$
\begin{aligned}
[Q]&=\left[\begin{matrix}0 & 1\\1 & 0\end{matrix}\right],\qquad Q(\omega,\omega)=20 \\
N&=\left[\begin{matrix}- \frac{2}{5}\\1\end{matrix}\right] \\
Q|_P&=\left[\begin{matrix}- \frac{4}{5}\end{matrix}\right]
\end{aligned}
$$

$Q(a,\omega)=0$ は、$a\wedge\omega^{n-1}$ の最高次係数が零になる条件である。従って $a$ を直線束の $iF_L$ の値とみなすと、$P$ はこの点で曲率の縮約が零になる全ての対角成分を与える。行列 $Q|_P=N^\top QN$ は、この条件下で非零な対角形式の二乗に $\omega^{n-2}$ を掛けた形式の係数が負になることを表す。
:::

:::example[積の計量で対角曲率の縮約を計算する]
複素次元 $3$ の一点で正則余接枠 $\eta_j$ を取り、$a_j=i\eta_j\wedge\bar\eta_j$ とする。$a_j^2=0$ である。最高次形式 $\beta$ の $a_1\wedge\cdots\wedge a_n$ に関する係数を $[\beta]_{\mathrm{top}}$ と書く。基本形式には次の正係数の線形結合を使う。

$$
\begin{aligned}
\omega&=a_{1} + 2 a_{2} + 4 a_{3}
\end{aligned}
$$

外積作用素 $L=\omega\wedge$ の $n-2$ 乗は、$a_j$ の張る実空間を $n-1$ 個の添字の積が張る空間に写す。上位側の基底を昇順の添字の積とし、$a_i$ に $\omega^{n-2}$ を掛けると、同じ添字に二度出会う項が消える。残りの各添字集合の係数には並べ替えの数 $(n-2)!$ が掛かる。

$$
\begin{aligned}
[L^{1}]&=\left[\begin{matrix}2 & 1 & 0\\4 & 0 & 1\\0 & 4 & 2\end{matrix}\right] \\
\det[L^{1}]&=-16
\end{aligned}
$$

最高次の係数 $Q(a,b)=[a\wedge b\wedge\omega^{n-2}]_{\mathrm{top}}$ は、行と列が異なるときに残りの添字の重みの積から求まる。$Q(a,\omega)=0$ となる空間を $P$ とし、基底を列に並べた行列を $N$ とする。

$$
\begin{aligned}
[Q]&=\left[\begin{matrix}0 & 4 & 2\\4 & 0 & 1\\2 & 1 & 0\end{matrix}\right],\qquad Q(\omega,\omega)=48 \\
N&=\left[\begin{matrix}- \frac{1}{2} & - \frac{1}{4}\\1 & 0\\0 & 1\end{matrix}\right] \\
Q|_P&=\left[\begin{matrix}-4 & -1\\-1 & -1\end{matrix}\right]
\end{aligned}
$$

$Q(a,\omega)=0$ は、$a\wedge\omega^{n-1}$ の最高次係数が零になる条件である。従って $a$ を直線束の $iF_L$ の値とみなすと、$P$ はこの点で曲率の縮約が零になる全ての対角成分を与える。行列 $Q|_P=N^\top QN$ は、この条件下で非零な対角形式の二乗に $\omega^{n-2}$ を掛けた形式の係数が負になることを表す。
:::

:::example[積の直線束でコホモロジーが残る次数]
$X=\prod_j\mathbb P^{n_j}$ とし、外部テンソル積の直線束 $L=\mathcal O(2,3)$ を考える。$a$ のべきでコホモロジー次数を記録し、係数を各空間の次元とする。

一因子では符号がすべて揃った Taylor–Laurent 指数だけが残る。指数が非負の場合は全て $0$ 次に、全て負の場合は因子の最高次に属する。負指数を $-b_i-1$ と書けば、弱い合成の個数を二項係数で求められる。

$$
\begin{aligned}
j=1&\colon\quad \mathbb P^{1},\; d=2,\quad e_i\ge0,\quad \sum_i e_i=2,\quad P_j(a)=3 \\
j=2&\colon\quad \mathbb P^{2},\; d=3,\quad e_i\ge0,\quad \sum_i e_i=3,\quad P_j(a)=10
\end{aligned}
$$

積では各因子の Dolbeault Laplacian の非負固有値を加える。零になるには各因子の固有値が零でなければならず、調和代表の次元と次数は多項式の積で合成される。

$$
\begin{aligned}
P_X(a)&=\prod_j P_j(a)=30 \\
\dim_\C X&=3 \\
\chi(X,L)&=P_X(-1)=30
\end{aligned}
$$

随伴の双対束 $K_X\otimes L^*$ は $\mathcal O(-4,-6)$ である。$d\ge0$ の因子では双対側が最高次になり、$d\le-n-1$ の因子では逆になる。全消滅の範囲はそのままなので、双対は $q$ 次と $3-q$ 次を対応させる。正の成分を持つ外部テンソル積はこの標準積計量に対し縮約が一定である。
:::

### 傾きが等しいときに起きる分解

$E$ が Hermite–Einstein 条件を満たし、$S$ が正則部分束であるとする。:ref[prop-cg-subbundle-curvature] で $B=\sum_j B_j\theta_j$ と書くと、$i\Lambda(B^*\wedge B)=-\sum_j B_j^*B_j$ である。$\mathrm d\bar z_j\wedge\mathrm dz_j$ を逆順にするため負号が現れる。従って
$$
2\pi\deg_\omega S
=\lambda\,\rank(S)\operatorname{Vol}_\omega X-\int_X |B|^2\,\frac{\omega^n}{n!}.
$$
左辺は跡の縮約を積分した値である。$\lambda$ の式を代入して階数で割ると $\mu_\omega(S)\le\mu_\omega(E)$ が得られる。等号なら $B$ は積分が零の連続なノルム二乗を持つので全点で零である。このとき逆向きの成分 $-B^*$ も零で、$D$ は直交分解を保つ。特に $(0,1)$ 部分も両成分を保つから、$E$ は正則に $S\oplus Q$ へ分解し、各成分が同じ定数 $\lambda$ の条件を満たす。この議論は部分束に対するものである。特異な商を持つ部分層まで含む定理では、特異集合を避けた計算に加え、次数とエネルギーの境界項を制御する必要がある。

:::example[等しい傾きの非分裂拡大]
楕円曲線 $X$ では $H^1(\mathcal O_X)\ne0$ である。非零の加法 Čech 類 $c_{ij}$ で遷移を $\left(\begin{smallmatrix}1&c_{ij}\\0&1\end{smallmatrix}\right)$ とすると、三重交差の条件は行列の積でちょうど加法条件となる。第一成分は自明部分束、商も自明で、行列式も $1$ である。直線束としての両者の傾きは零である。もし正則な分裂があれば、商の定数切断 $1$ の持ち上げの第一成分の差が $c_{ij}$ を消すので、類は零になる。従って分裂せず、上の等号条件によって Hermite–Einstein 計量は存在しない。
:::

## 演習問題

::::exercise{level=basic}
自明直線束の計量 $h=e^{-3|z|^2}$ の Chern 接続と曲率を求めよ。
:::solution
$A=\partial\log h=-3\bar z\dd z$ である。$\bar\partial A=-3\dd\bar z\wedge\dd z=3\dd z\wedge\dd\bar z$ が曲率になる。従って $iF$ は正である。
:::
::::

::::exercise{level=basic}
$h=e^{-|z|^2}$ を $e'=e^{az}e$、$a\in\C$ の枠で書いたときの接続を求めよ。
:::solution
$\log h'=-|z|^2+az+\bar a\bar z$ なので $A'=(-\bar z+a)\dd z$ である。元との差は $a\dd z$ で、その外微分は零である。よって曲率は元と同じである。
:::
::::

::::exercise{level=basic}
$E=L\oplus L^{-1}$、$c_1(L)=H$ とする。$c_1(E),c_2(E)$ を求めよ。
:::solution
直和公式により $c(E)=(1+H)(1-H)=1-H^2$ である。よって $c_1=0,c_2=-H^2$ となる。第一 Chern 類が零でも高次は零とは限らない。
:::
::::

::::exercise{level=basic}
$\mathcal O_{\mathbb P^1}(3)$ の Chern 形式の積分を求めよ。
:::solution
テンソル三乗では局所重みが三倍なので曲率は $\mathcal O(1)$ の三倍である。その第一 Chern 形式の積分は $1$ だから求める値は $3$ である。
:::
::::

::::exercise{level=advanced}
$E$ の曲率が零なら正次数の Chern 類の de Rham 像が全て零になることを示せ。
:::solution
$F_\nabla=0$ なら $\det(I+iF_\nabla/(2\pi))=1$ であり、正次数の Chern 形式が全て零である。一般の接続で計算しても類は同じである。結論は de Rham 像についてのもので、整数係数の torsion の有無はこの計算では区別できない。
:::
::::

::::exercise{level=advanced}
正則直線束 $L,M$ の計量のテンソル積について曲率が加法的であることを証明せよ。
:::solution
局所枠 $e,f$ のノルム二乗を $h_L,h_M$ とすると $e\otimes f$ のノルム二乗は積である。実対数は和になるため、$\bar\partial\partial$ を取ると $F_{L\otimes M}=F_L+F_M$ である。各項は枠によらないため大域的にも成立する。
:::
::::

::::exercise{level=advanced}
$X=\mathbb P^1$、$\omega=iF_{\mathcal O(1)}$ として $E=\mathcal O(1)\oplus\mathcal O(3)$ の直和計量が Hermite–Einstein か判定せよ。
:::solution
曲率の縮約は $\diag(1,3)$ であり、スカラー行列でないため満たさない。仮に満たす計量があれば定数は $\deg E/\rank E=(1+3)/2=2$ となる。部分束 $\mathcal O(3)$ の傾き $3$ がそれより大きいので多安定でなく、対応定理から別の計量でも満たせない。
:::
::::

::::exercise{level=advanced}
$h=e^{-(|z_1|^2-|z_2|^2)}$ の直線束計量の曲率は正か。
:::solution
重みの複素 Hessian は $\diag(1,-1)$ である。従って $iF$ は第一座標方向で正、第二座標方向で負となり、この計量の曲率は正定値ではない。底は $\C^2$ なので、同じ自明束でも重みを $|z_1|^2+|z_2|^2$ に選べば正にできる。
:::
::::

::::exercise{level=advanced}
$E=\mathcal O^{\oplus2}$ に標準平坦計量を入れ、$\mathbb P^1$ の各点 $[Z_0:Z_1]$ に対し、$(Z_0^2,Z_1^2)$ が張る直線を部分束 $S$ とする。$S$ とその商の曲率および次数を求めよ。

:::solution
同次代表を定数倍すると二成分に同じ平方が掛かり、二成分は同時に零にならない。従って確かに正則部分束である。有限チャートで枠を $(1,z^2)$ とすると、平方ノルムは $h=1+|z|^4$ である。$\partial_z\partial_{\bar z}\log(1+|z|^4)=4|z|^2/(1+|z|^4)^2$ なので
$$
F_S=-\frac{4|z|^2}{(1+|z|^4)^2}\mathrm dz\wedge\mathrm d\bar z .
$$
とくに原点で曲率は零だが、原点の外では有限チャートの係数は負となる。元の束の行列式接続は平坦であり、部分束と商の対角接続の跡の和は零なので、直線束である商の曲率は $-F_S$ である。

極座標で $\deg S=-\int_0^\infty 8r^3/(1+r^4)^2\,\mathrm dr=-2$ である。$s=r^4$ と置けば積分の絶対値は $2\int_0^\infty(1+s)^{-2}\mathrm ds=2$ となる。従って商の次数は $2$ である。別の説明では $S$ は $[Z_0:Z_1]\mapsto[Z_0^2:Z_1^2]$ で標準直線束を引き戻したものなので $\mathcal O(-2)$ であり、次数も一致する。
:::
::::

## 確認テスト

::::quiz{answer="5" tol="1e-6"}
$E=L^{\otimes2}\oplus L^{\otimes3}$、$c_1(L)=H$ のとき $c_1(E)$ の $H$ の係数を求めよ。
:::explanation
直和公式で一次項は $2H+3H=5H$ である。積 $6H^2$ は第二 Chern 類に現れる。
:::
::::

::::quiz
正則枠で Chern 接続の接続行列を選べ。

- [x] $H^{-1}\partial H$
- [ ] $H^{-1}\bar\partial H$
- [ ] $H^{-1}\dd H$

:::explanation
正則枠で $(0,1)$ 接続行列は零であり、計量との両立の $(1,0)$ 成分は $\partial H=HA$ である。従って $A=H^{-1}\partial H$ となる。
:::
::::

::::quiz{answer="-2" tol="1e-6"}
$\mathcal O_{\mathbb P^1}(-2)$ の第一 Chern 形式の積分を求めよ。
:::explanation
$\mathcal O(1)$ の積分は $1$ で、双対では曲率の符号が逆になる。$\mathcal O(-2)$ は双対の二乗なので積分は $-2$ である。
:::
::::

::::quiz
直線束の重み $\varphi=-\log h$ から計算する正値性の条件を選べ。

- [ ] $\varphi$ の値が全て正である。
- [x] 複素 Hessian $(\partial_j\partial_{\bar k}\varphi)$ が正定値である。
- [ ] $\varphi$ の一次導関数が全て零である。

:::explanation
曲率は $\partial\bar\partial\varphi$ であり、$iF$ の係数行列がこの複素 Hessian である。重みに定数を加えても条件は変わらない。
:::
::::

::::quiz{answer="3" tol="1e-6"}
$\mathbb P^1\times\mathbb P^1$ の $E=\mathcal O(1,2)\oplus\mathcal O(2,-1)\oplus\mathcal O(-3,1)$ に対して $\int c_2(E)$ を求めよ。

:::explanation
各直線束の第一類を $a+2b,2a-b,-3a+b$ と書く。$a^2=b^2=0$ により二本ごとの積は $3ab,-5ab,5ab$ である。これらを足すと $3ab$ となり、$\int ab=1$ から答えは $3$ である。
:::
::::

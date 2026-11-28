---
title: 集団運動
description: 核の変形、回転、表面振動、巨大共鳴を多核子の集団励起として記述する。
notation:
  - symbol: '\alpha_{\lambda\mu}'
    meaning: 核表面の多重極変形係数
  - symbol: '\mathcal I'
    meaning: 慣性能率
  - symbol: 'K'
    meaning: 固有座標の対称軸への角運動量の射影
  - symbol: 'm_1'
    meaning: エネルギー重み付き遷移強度和
---

## この章の目標

- 多重極係数で核の固有形状を表す。
- 回転バンドと表面振動のエネルギーを導く。
- 粒子穴の集団運動を乱雑位相近似で記述する。
- 巨大共鳴と高スピンでのバンド交差を説明する。

## 原子核の変形と Nilsson 模型

実験室の $J=0$ 状態の平均密度は回転対称である。一方、回転させる前の固有座標では変形した密度を選ぶと準位の関係を記述できる。小変形の表面を
$$
R(\theta,\phi)=R_0\left[1+\sum_{\lambda\ge2,\mu}\alpha_{\lambda\mu}Y_{\lambda\mu}(\theta,\phi)\right]
$$
と表す。表面を実数にする条件は $\alpha_{\lambda,-\mu}=(-1)^\mu\alpha_{\lambda\mu}^*$ である。$\lambda=0$ の成分は体積の調整に、$\lambda=1$ は一次では中心の並進に対応する。球面平均を取れば $\lambda\ge2$ の変形で体積は一次まで一定だが、二次では単極子の補正が必要である。軸対称の $\alpha_{20}=\beta_2$ が正なら極が伸び、負なら扁平になる。

:term[Nilsson 模型]{en="Nilsson model" yomi="にるそんもけい"} では、変形した振動子にスピン軌道項と軌道角運動量の二乗に依存する項を加えて単粒子準位を求める。軸対称振動子を
$$
H_{\rm osc}=\frac{\vb p^2}{2m}
+\frac m2\left[\omega_\perp^2(x^2+y^2)+\omega_z^2 z^2\right]
$$
と置くと波動関数が変数分離し、$E=\hbar\omega_\perp(n_x+n_y+1)+\hbar\omega_z(n_z+1/2)$ となる。体積を一定にするための規約の一つは $\omega_\perp^2\omega_z=\omega_0^3$ である。三次元の $j$ は一般に保存されなくなるが、軸回転に対称なら全角運動量の軸射影とパリティで単粒子を分類できる。実際の準位には [殻模型](/physics/nuclear-physics/shell-model/) のスピン軌道項が大きく効く。

球形振動子の主殻 $N$ に射影する模型では、無次元の $\vb l=\vb L/\hbar,\vb s=\vb S/\hbar$ を用いて補正を
$$
-\kappa_N\hbar\omega_0\left[
2\vb l\cdot\vb s+\mu_N\left(\vb l^2-\frac{N(N+3)}2\right)\right]
$$
と書く。ここで $\kappa_N,\mu_N$ は無次元の模型係数である。$\vb l^2$ から引いた量は、$\ell=N,N-2,\ldots\ge0$ の固有値 $\ell(\ell+1)$ を縮退度 $2\ell+1$ で平均した値である。実際、重み付き和は $N(N+1)(N+2)(N+3)/4$ で、重みの和は $(N+1)(N+2)/2$ である。こうして殻の平均位置を保ったまま、高 $\ell$ と低 $\ell$ の相対位置を補正する。

強い軸変形で漸近的に成分を同定できるときは、$\Omega[Nn_z\Lambda]$ を軌道のラベルに使う。$\Lambda$ は軌道角運動量の軸射影、$\Omega=\Lambda+\Sigma$ はスピン射影 $\Sigma=\pm1/2$ を含む全射影である。$\Lambda,n_z$ は一般には近似的な量子数である。同じ $\Omega$ とパリティの準位は混ざり得るが、異なる厳密な対称性を持つ場合は混ざらない。変形準位の分類と回転バンドとの関係は :cite[bohr-mottelson] を参照する。

### 四重極変形の二つの変数

四重極の主軸を選ぶと、実係数を $\alpha_{20}=\beta_2\cos\gamma$、$\alpha_{22}=\alpha_{2,-2}=\beta_2\sin\gamma/\sqrt2$ と表せる。$\beta_2\ge0$ は四重極変形の大きさ、$\gamma$ は三軸性を表す角度である。この規約では軸の入れ替えで $0\le\gamma\le\pi/3$ を取れ、両端がそれぞれ伸長型と扁平型、内点が三軸型である。軸を固定した前の表記では伸長と扁平を $\alpha_{20}$ の符号で分けていた。

小さい軸変形 $\alpha_{20}$ と均一な電荷密度 $\rho_e=3Ze/(4\pi R_0^3)+O(\alpha_{20}^2)$ なら、固有四重極モーメントは $Q_0=\sqrt{16\pi/5}\int r^2Y_{20}\rho_e\dd^3r$ である。動径を表面まで積分して $R^5/5$ とし、$R=R_0(1+\alpha_{20}Y_{20})$ を一次まで展開する。球の項は角度平均で消え、$\int Y_{20}^2\dd\Omega=1$ から
$$
Q_0=\frac{3}{\sqrt{5\pi}}Ze R_0^2\alpha_{20}+O(\alpha_{20}^2)
$$
となる。変形の大きさを遷移から求めるときは、この近似と次節の回転模型の両方を用いている。

:::example[非球形での分裂]
三振動子の一量子励起は、$z$ 方向が $\hbar\omega_z$、横が $\hbar\omega_\perp$ の増分を持つ。$\omega_z<\omega_\perp$ なら $z$ 方向の節を持つ励起が低くなる。球形では三方向が再び縮退する。
:::

:::example[体積一定の小変形]
$\omega_z=\omega_0e^{-2d}$、$\omega_\perp=\omega_0e^d$ とすれば積は常に $\omega_0^3$ である。$d>0$ で $z$ 方向の振動子長は球形の $e^d$ 倍に伸びる。
:::

### 体積保存と変形振動子の量子数

軸対称振動子を $\omega_x=\omega_y=\omega_\perp,\omega_z$ で定義する。等ポテンシャル面の半軸は各振動数に反比例するので、面の体積を変えない条件は $\omega_\perp^2\omega_z=\omega_0^3$ である。例えば $\omega_\perp=\omega_0e^\epsilon$、$\omega_z=\omega_0e^{-2\epsilon}$ とおくと任意の $\epsilon$ でこれが成立する。$\epsilon>0$ では $z$ 軸が長い。直交基底のエネルギーを微分すると、球形からの傾きは $\hbar\omega_0(n_x+n_y-2n_z)$ である。従って同じ主殻の中でも $n_z$ の大きい軌道は伸長によって下がる。スピン軌道力を結合した場合に残る厳密な量子数は軸射影とパリティであり、振動子の $n_z$ は一般に近似的なラベルとなる。

:::example[変形した振動数の比]
体積保存の指数表示を使う。周方向と軸方向の振動数の比は $e^{3\epsilon}$ である。

量の定義に従い、共通の単位で次の値を用いる。
$$
eps=\frac{1}{20}
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=e^{3 \cdot eps}\\&=e^{3 \cdot \left(\frac{1}{20}\right)}\\&\simeq 1.16183.\end{aligned}
$$

同じ符号では振動数比が大きいほど等ポテンシャル面が伸びる。

変形の大きい設定と比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
eps=\frac{1}{10}
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=e^{3 \cdot \left(\frac{1}{10}\right)}\\&\simeq 1.34986,\\{R}_{\mathrm{new}}-{R}&\simeq 0.188025.\end{aligned}
$$
:::

:::example[同じ主殻の傾き]
主殻を $N$ とすると $n_x+n_y=N-n_z$ である。球形からのエネルギーの一次変化を求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
W=10\,\mathrm{MeV},\qquad eps=\frac{3}{50},\qquad N=2,\qquad n_{z}=0
$$

この条件で求める量を ${\Delta E}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\Delta E}&=W \cdot eps \cdot \left(N - 3 \cdot n_{z}\right)\\&=\left(10\right) \cdot \left(\frac{3}{50}\right) \cdot \left(- 3 \cdot \left(0\right) + \left(2\right)\right)\\&\simeq 1.2\,\mathrm{MeV}.\end{aligned}
$$

同じ主殻の中で符号が反対になる軌道がある。この分裂は変形を支持する占有配置を選別する。

軸方向の節の多い軌道と比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
n_{z}=2
$$

このときの値を ${\Delta E}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\Delta E}_{\mathrm{new}}&=- 2 \cdot \left(10\right) \cdot \left(2\right) \cdot \left(\frac{3}{50}\right)\\&\simeq -2.4\,\mathrm{MeV},\\{\Delta E}_{\mathrm{new}}-{\Delta E}&\simeq -3.6\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[伸長側の等ポテンシャル面]
横半軸を $a$、縦半軸を $c$ とする。$a^2c=R_0^3$ から縦の長さを求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
R_{0}=5\,\mathrm{fm},\qquad a=\frac{49}{10}\,\mathrm{fm}
$$

この条件で求める量を ${c}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{c}&=\frac{R_{0}^{3}}{a^{2}}\\&=\frac{\left(5\right)^{3}}{\left(\frac{49}{10}\right)^{2}}\\&\simeq 5.20616\,\mathrm{fm}.\end{aligned}
$$

これは体積と形を分けるための幾何学である。量子密度の輪郭全体を一つの楕円体で表す近似には限界がある。

横半径を小さくする。そこで次の値に変更し、ほかの入力を固定する。
$$
a=\frac{9}{2}\,\mathrm{fm}
$$

このときの値を ${c}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{c}_{\mathrm{new}}&=\frac{\left(5\right)^{3}}{\left(\frac{9}{2}\right)^{2}}\\&\simeq 6.17284\,\mathrm{fm},\\{c}_{\mathrm{new}}-{c}&\simeq 0.966675\,\mathrm{fm}.\end{aligned}
$$

**一次近似の検証。** 横方向の半径の誤差を、体積を固定したときの縦方向の半径へ伝える。

${a}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${c}$ の一階微分を求めると、
$$
\frac{\mathrm d{c}}{\mathrm da}=- \frac{2 \cdot R_{0}^{3}}{a^{3}}.
$$

この導関数を変更前の点で評価して $h=- \frac{2}{5}$ を掛ける。ここで $a$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq 0.849986\,\mathrm{fm},\\\Delta-\Delta_{\mathrm{lin}}&\simeq 0.116689\,\mathrm{fm}.\end{aligned}
$$

残差が示すのは、指定した二点間で生じる線形化の誤差である。基準点からさらに離れた値について同じ精度を保証するものではない。
:::

:::example[軸方向の振動量子数が異なる準位]
純振動子の同じ $N$ で、$n_z$ の差を $d_n$ とする。エネルギー差は $W d_n(e^{-2\epsilon}-e^\epsilon)$ となる。

量の定義に従い、共通の単位で次の値を用いる。
$$
W=10\,\mathrm{MeV},\qquad d_{n}=1,\qquad eps=- \frac{3}{50}
$$

この条件で求める量を ${\Delta E}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\Delta E}&=W \cdot d_{n} \cdot \left(- e^{eps} + e^{- 2 \cdot eps}\right)\\&=\left(10\right) \cdot \left(1\right) \cdot \left(- e^{\left(- \frac{3}{50}\right)} + e^{- 2 \cdot \left(- \frac{3}{50}\right)}\right)\\&\simeq 1.85732\,\mathrm{MeV}.\end{aligned}
$$

扁平と伸長で順位が逆になる。

変形の符号を反転する。そこで次の値に変更し、ほかの入力を固定する。
$$
eps=\frac{3}{50}
$$

このときの値を ${\Delta E}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\Delta E}_{\mathrm{new}}&=\left(10\right) \cdot \left(1\right) \cdot \left(- e^{\left(\frac{3}{50}\right)} + e^{- 2 \cdot \left(\frac{3}{50}\right)}\right)\\&\simeq -1.74916\,\mathrm{MeV},\\{\Delta E}_{\mathrm{new}}-{\Delta E}&\simeq -3.60648\,\mathrm{MeV}.\end{aligned}
$$
:::

## 回転バンドと慣性能率

軸対称な変形の向きを Euler 角で量子化する。固有状態の軸射影を $\hbar K$、それに垂直な二方向の慣性能率を $\mathcal I$ とする。回転の横成分の運動エネルギーは $(\vb J^2-J_3^2)/(2\mathcal I)$ なので、固定 $K$ の系列内では
$$
E(J)-E(K)=\frac{\hbar^2}{2\mathcal I}[J(J+1)-K(K+1)] .
\label{eq-nuc-rotor}
$$
$K=0$ の偶偶核基底バンドで、垂直な軸の周りの $\pi$ 回転に対して固有状態が不変とする。この回転は固有軸の向きを反転し、回転波動関数の $K=0$ 成分には $(-1)^J$ を掛ける。従って不変な系列は $J=0,2,4,\ldots$ である。$\mathcal I$ は [対相関](/physics/nuclear-physics/pairing/) の影響を受け、剛体値とは一般に異なる。

軸対称な固有四重極モーメントを $Q_0$ とし、固有テンソルの軸成分を $\sqrt{5/(16\pi)}Q_0$ とする。実験室の成分はこれに回転行列 $D^2_{\mu0}$ を掛ける。回転波動関数にもう一つ $D^J_{MK}$ があるので、積を結合角運動量に展開すると、空間射影と固有軸射影の二つの Clebsch–Gordan 係数が現れる。空間射影を和し、初期射影数で平均すると、同一 $K=0$ バンドで
$$
B(E2;J_i\to J_f)=\frac{5Q_0^2}{16\pi}
\big|\langle J_i0,20|J_f0\rangle\big|^2
\label{eq-nuc-rotor-e2}
$$
となる。$B(E2)$ の単位は $e^2\,\mathrm{fm}^4$ である。

:::example[準位比]
$K=0$ で $E(2)=90\,\mathrm{keV}$ なら $E(4)=90\cdot20/6=300\,\mathrm{keV}$ である。このとき $\mathcal I/\hbar^2=3/(0.09\,\mathrm{MeV})\simeq33.3\,\mathrm{MeV}^{-1}$ である。
:::

:::example[上向きと下向き]
$\langle 00,20|20\rangle=1$、$\langle20,20|00\rangle=1/\sqrt5$ なので、$B(E2;0\to2)$ は逆向きの 5 倍である。行列要素の絶対値が共通でも初期状態の射影数が異なる。
:::

### 静的な四重極と遷移の四重極

軸対称の $K=0$ 回転波動関数に対して、実験室で磁気射影 $M=J$ の四重極を求める。固有成分は $Q_0$ なので演算子は $Q_0D^2_{00}$ である。回転行列の積の直交積分から、対角行列要素は
$$
Q_s=Q_0\langle J0,20|J0\rangle\langle JJ,20|JJ\rangle
=-\frac{J}{2J+3}Q_0
$$
となる。ここで $J\ge1$ とした。係数を確かめるには、階数 2 の軸成分を同じ $J$ 内で $3J_z^2-\vb J^2$ に比例させる。両方の Clebsch–Gordan 係数の分母は $\sqrt{J(J+1)(2J-1)(2J+3)}$、分子はそれぞれ $-J(J+1)$ と $J(2J-1)$ であり、積は上の結果になる。$J=0$ では角運動量の三角条件から $Q_s=0$ である。

:::example[固有モーメントからの二つの予測]
$K=0$ の変形核で $Q_0=100\,e\,\mathrm{fm}^2$ とする。$2^+$ の分光学的モーメントは $-200/7\simeq-28.6\,e\,\mathrm{fm}^2$ である。一方 :eqref[eq-nuc-rotor-e2] で $2^+\to0^+$ の強度は $Q_0^2/(16\pi)\simeq199\,e^2\mathrm{fm}^4$ となる。固有の正の四重極は実験室での静的な負の値と両立する。
:::

### 対称こまと遷移強度の比

対称こまの運動エネルギーは $\{J(J+1)-K^2\}\hbar^2/(2\mathcal I_\perp)+K^2\hbar^2/(2\mathcal I_\parallel)$ である。固定 $K$ の差分では後者と前者の $K^2$ は定数として消える。電気四重極の実験室成分は軸対称なら $D^2_{\mu0}$ に比例する。$D^{J_i}_{M_iK}D^2_{\mu0}$ を二つの角運動量の直積表現とみなし、合成基底に展開すると係数は空間固定軸の合成と体固定軸の合成の二つの Clebsch–Gordan 係数の積になる。$D$ の直交積分で終状態を選び、実験室射影を和と平均で消すと $B(E2)$ は体固定軸の係数の二乗に比例する。この導出のため、同じバンド内の比では固有四重極モーメントが消える。

:::example[回転遷移の間隔]
回転定数 $a=\hbar^2/(2\mathcal I)$ とする。$J\to J-2$ の差を取る。

量の定義に従い、共通の単位で次の値を用いる。
$$
a=\frac{1}{50}\,\mathrm{MeV},\qquad J=4
$$

この条件で求める量を ${E_\gamma}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{E_\gamma}&=a \cdot \left(4 \cdot J - 2\right)\\&=\left(\frac{1}{50}\right) \cdot \left(4 \cdot \left(4\right) - 2\right)\\&\simeq 0.28\,\mathrm{MeV}.\end{aligned}
$$

回転定数が固定なら線間隔は一定の増分を示す。

より高い回転準位と比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
J=8
$$

このときの値を ${E_\gamma}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{E_\gamma}_{\mathrm{new}}&=\left(\frac{1}{50}\right) \cdot \left(4 \cdot \left(8\right) - 2\right)\\&\simeq 0.6\,\mathrm{MeV},\\{E_\gamma}_{\mathrm{new}}-{E_\gamma}&\simeq 0.32\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[バンド内の四重極行列要素]
$K=0$ で下降遷移の結合係数の二乗は $3J(J-1)/[2(2J-1)(2J+1)]$ である。$2\to0$ の $1/5$ に対する比を求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
J=4
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=\frac{15 \cdot J \cdot \left(J - 1\right)}{2 \cdot \left(2 \cdot J + 1\right) \cdot \left(2 \cdot J - 1\right)}\\&=\frac{15 \cdot \left(4\right) \cdot \left(\left(4\right) - 1\right)}{2 \cdot \left(2 \cdot \left(4\right) - 1\right) \cdot \left(2 \cdot \left(4\right) + 1\right)}\\&\simeq 1.42857.\end{aligned}
$$

バンド内で $Q_0$ が同一と仮定している。遷移エネルギーの五乗は部分寿命を比べるときだけ追加する。

さらに一段上の下降遷移と比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
J=6
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=\frac{15 \cdot \left(6\right) \cdot \left(\left(6\right) - 1\right)}{2 \cdot \left(2 \cdot \left(6\right) - 1\right) \cdot \left(2 \cdot \left(6\right) + 1\right)}\\&\simeq 1.57343,\\{R}_{\mathrm{new}}-{R}&\simeq 0.144855.\end{aligned}
$$

**一次近似の検証。** 整数スピン間の比較とは別に、係数の有理式を連続化して高い $J$ の変化の尺度を見積もる。

${J}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${R}$ の一階微分を求めると、
$$
\frac{\mathrm d{R}}{\mathrm dJ}=\frac{15 \cdot \left(4 \cdot J^{2} - 2 \cdot J + 1\right)}{2 \cdot \left(16 \cdot J^{4} - 8 \cdot J^{2} + 1\right)}.
$$

この導関数を変更前の点で評価して $h=2$ を掛ける。ここで $J$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq 0.21542,\\\Delta-\Delta_{\mathrm{lin}}&\simeq -0.0705644.\end{aligned}
$$

省いた分が全変化に比べて大きい。従って増分の大きさを予測するには、少なくとも二次の項を調べるか、元の式を使う必要がある。
:::

:::example[剛体球からの期待値]
等質量の $A$ 核子の球で古典的には $\mathcal I=2AmR^2/5$ である。回転定数を求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
h_{c}=\frac{197327}{1000}\,\mathrm{MeV\,fm},\qquad A=150,\qquad m_{c}=939\,\mathrm{MeV},\qquad R=6\,\mathrm{fm}
$$

この条件で求める量を ${a}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{a}&=\frac{5 \cdot h_{c}^{2}}{4 \cdot A \cdot R^{2} \cdot m_{c}}\\&=\frac{5 \cdot \left(\frac{197327}{1000}\right)^{2}}{4 \cdot \left(150\right) \cdot \left(6\right)^{2} \cdot \left(939\right)}\\&\simeq 0.00959895\,\mathrm{MeV}.\end{aligned}
$$

これは質量の分布だけから作った基準である。核の有効慣性能率は対相関や変形にも依存する。

半径の大きい模型を比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
R=\frac{13}{2}\,\mathrm{fm}
$$

このときの値を ${a}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{a}_{\mathrm{new}}&=\frac{5 \cdot \left(\frac{197327}{1000}\right)^{2}}{4 \cdot \left(150\right) \cdot \left(939\right) \cdot \left(\frac{13}{2}\right)^{2}}\\&\simeq 0.00817899\,\mathrm{MeV},\\{a}_{\mathrm{new}}-{a}&\simeq -0.00141996\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[固定軸の分光学的四重極]
$K=0$ の回転模型では $Q_s/Q_0=-J/(2J+3)$ となる。

量の定義に従い、共通の単位で次の値を用いる。
$$
J=2
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=\frac{\left(-1\right) \cdot J}{2 \cdot J + 3}\\&=- \frac{\left(2\right)}{2 \cdot \left(2\right) + 3}\\&\simeq -0.285714.\end{aligned}
$$

伸長した正の固有モーメントでも、実験室の伸張状態の符号が負になる。

高いスピンを比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
J=6
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=- \frac{\left(6\right)}{2 \cdot \left(6\right) + 3}\\&\simeq -0.4,\\{R}_{\mathrm{new}}-{R}&\simeq -0.114286.\end{aligned}
$$
:::

## 表面振動（四重極振動と八重極振動）

球形に安定な核では、小さい $\alpha_{\lambda\mu}$ に復元力が働く。実基底で同じ $\lambda$ の $2\lambda+1$ 成分を $q_i$ とすると、回転不変な二次の Hamiltonian は
$$
H_\lambda=\sum_{i=1}^{2\lambda+1}\left(\frac{p_i^2}{2B_\lambda}+\frac12 C_\lambda q_i^2\right)
$$
と書ける。$B_\lambda,C_\lambda>0$ とすれば、運動方程式は $B_\lambda\ddot q_i=-C_\lambda q_i$ である。量子化した一振動子の増分は $\hbar\omega_\lambda$、$\omega_\lambda=\sqrt{C_\lambda/B_\lambda}$ である。フォノンの作成演算子が階数 $\lambda$、パリティ $(-1)^\lambda$ で変換するので、基底 $0^+$ から一つ作ると四重極で $2^+$、八重極で $3^-$ を得る。

同じ種のフォノンを二つ作ると、二粒子の結合による交換位相 $(-1)^{2\lambda-J}$ が $+1$ でなければならない。従って四重極の二フォノンは $J=0,2,4$ で正パリティである。調和近似では同じエネルギー $2\hbar\omega_2$ に現れ、非調和性と別の配位との混合が分裂を生む。

:::example[振動と回転の区別]
一フォノンが $1\,\mathrm{MeV}$ なら、調和模型の二フォノンは $2\,\mathrm{MeV}$ である。回転バンドで最初の $2^+$ が $1\,\mathrm{MeV}$ なら $4^+$ は $10/3\,\mathrm{MeV}$ である。準位比に加え遷移強度も使って比較する。
:::

:::example[選択則]
線形の変位 $q_i\propto b_i+b_i^\dagger$ はフォノン数を一つ変える。そのため二フォノンから基底状態への直接の線形振幅はゼロであり、二次の演算子や混合で生じる。
:::

### 非圧縮な液滴の振動数

表面張力を $\sigma_s$、質量密度を $\rho_m$ とした非圧縮液滴で、規格化 $\int Y^2\dd\Omega=1$ の実球面調和関数 $Y$ を選ぶ。表面変位を $R_0qY$ とし、液体内部の渦なし速度を $\grad\phi$ とする。非圧縮性は $\nabla^2\phi=0$、運動学的な境界条件は $\partial_r\phi|_{R_0}=R_0\dot qY$ である。原点で正則な解は $\phi=R_0^2\dot q(r/R_0)^\lambda Y/\lambda$ である。Green の公式で内部の運動エネルギーを表面積分にすると
$$
T=\frac{\rho_m}2\int|\grad\phi|^2\dd^3r
=\frac{\rho_m}2\int_{\partial V}\phi\partial_r\phi\,\dd S
=\frac{\rho_m R_0^5}{2\lambda}\dot q^2 .
$$
体積を保つ一様な二次補正を入れて面積を展開すると、面積の増分は $R_0^2[-\int(qY)^2+\tfrac12\int|\grad_\Omega(qY)|^2]$ である。$\nabla_\Omega^2Y=-\lambda(\lambda+1)Y$ より、$C_\lambda=\sigma_s R_0^2(\lambda-1)(\lambda+2)$ となる。従って $\omega_\lambda^2=\sigma_s\lambda(\lambda-1)(\lambda+2)/(\rho_mR_0^3)$ である。Coulomb、圧縮性、単粒子励起を捨てた模型なので、励起準位の値そのものには補正が要る。

:::example[四重極と八重極の尺度]
この液滴では $\omega_3^2/\omega_2^2=30/8=15/4$ となる。なお $\lambda=1$ では復元力がゼロであり、液滴全体の並進に対応する。
:::

### 遷移で数えるフォノン

調和振動の一成分の変位を $q=\sqrt{\hbar/(2B\omega)}(b+b^\dagger)$ とする。規格化した占有基底で $\langle n-1|q|n\rangle=\sqrt{\hbar n/(2B\omega)}$ であり、下向きの強度は $n$ に比例する。球対称な四重極振動では五つの成分を持ち、二フォノンの多重項の一つを選ぶと角運動量結合の係数が入る。例えば二フォノンの $J=4,M=4$ は、最大の磁気射影 $\mu=2$ のフォノンを二つ占有した状態だけから成る。$4,M=4$ から $2,M=2$ への行列要素は $\sqrt2$ の係数を持つ。下降の演算子成分は $\mu=-2$ であり、$2\to0$ の伸び切った射影からの遷移と同様に、許された終射影と演算子成分を足した値を、初期射影一つあたりの強度として用いられる。従ってこの模型で $B(E2;4\to2)/B(E2;2\to0)=2$ である。

:::example[同じ強度比を回転と比較する]
$K=0$ 剛体回転では $\langle40,20|20\rangle^2=2/7$、$\langle20,20|00\rangle^2=1/5$ で、:eqref[eq-nuc-rotor-e2] の比は $10/7$ である。調和振動の 2 と比べると、準位比の他に遷移強度も模型の比較に使える。
:::

### 表面変位の量子化と非調和性

一つの実正準座標 $a$ に対し $H=p_a^2/(2B)+Ca^2/2$ とする。$a=\sqrt{\hbar/(2B\omega)}(b^\dagger+b)$、$p_a=i\sqrt{\hbar B\omega/2}(b^\dagger-b)$ を代入する。$[b,b^\dagger]=1$ で交差項の二乗が打ち消し合い $H=\hbar\omega(b^\dagger b+1/2)$ になる。従って座標の $n\to n-1$ の行列要素二乗は基底遷移の $n$ 倍になる。四重極・八重極の各成分に同じ振動数を与えると回転不変性が保たれ、パリティは一量子ごとの $(-1)^\lambda$ の積である。二量子の分裂や基底への直接遷移は、非調和なポテンシャルや演算子の非線形項を調べる手段になる。

:::example[四重極と八重極の表面張力]
非圧縮液滴で $\omega_\lambda^2\propto\lambda(\lambda-1)(\lambda+2)$ として、$\omega_2$ を基準に比を計算する。

量の定義に従い、共通の単位で次の値を用いる。
$$
l=2
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=\sqrt{\frac{l \cdot \left(l + 2\right) \cdot \left(l - 1\right)}{8}}\\&=\frac{\sqrt{2} \cdot \sqrt{\left(2\right) \cdot \left(\left(2\right) - 1\right) \cdot \left(\left(2\right) + 2\right)}}{4}\\&=1.\end{aligned}
$$

質量密度と表面張力を同じと仮定した比較である。有限核では殻構造と電荷による補正がある。

八重極の比を得る。そこで次の値に変更し、ほかの入力を固定する。
$$
l=3
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=\frac{\sqrt{2} \cdot \sqrt{\left(3\right) \cdot \left(\left(3\right) - 1\right) \cdot \left(\left(3\right) + 2\right)}}{4}\\&\simeq 1.93649,\\{R}_{\mathrm{new}}-{R}&\simeq 0.936492.\end{aligned}
$$
:::

:::example[零点変位]
エネルギー間隔 $W=\hbar\omega$ と剛性 $C$ を使い $\langle a^2\rangle=W/(2C)$ とする。

量の定義に従い、共通の単位で次の値を用いる。
$$
W=1\,\mathrm{MeV},\qquad C=30\,\mathrm{MeV}
$$

この条件で求める量を ${\sigma_a}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\sigma_a}&=\sqrt{\frac{W}{2 \cdot C}}\\&=\frac{\sqrt{2} \cdot \sqrt{\frac{\left(1\right)}{\left(30\right)}}}{2}\\&\simeq 0.129099.\end{aligned}
$$

座標 $a$ は無次元の表面変位である。ここでは観測した間隔を固定した比較である。

剛性を大きくする。そこで次の値に変更し、ほかの入力を固定する。
$$
C=60\,\mathrm{MeV}
$$

このときの値を ${\sigma_a}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\sigma_a}_{\mathrm{new}}&=\frac{\sqrt{2} \cdot \sqrt{\frac{\left(1\right)}{\left(60\right)}}}{2}\\&\simeq 0.0912871,\\{\sigma_a}_{\mathrm{new}}-{\sigma_a}&\simeq -0.0378124.\end{aligned}
$$

**一次近似の検証。** 同じ間隔を保つ模型間で、剛性の推定誤差が固有座標の揺らぎにどう効くか調べる。

${C}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${\sigma_a}$ の一階微分を求めると、
$$
\frac{\mathrm d{\sigma_a}}{\mathrm dC}=- \frac{\sqrt{2} \cdot \sqrt{\frac{W}{C}}}{4 \cdot C}.
$$

この導関数を変更前の点で評価して $h=30$ を掛ける。ここで $C$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq -0.0645497,\\\Delta-\Delta_{\mathrm{lin}}&\simeq 0.0267374.\end{aligned}
$$

省いた分が全変化に比べて大きい。従って増分の大きさを予測するには、少なくとも二次の項を調べるか、元の式を使う必要がある。
:::

:::example[多フォノンの下降強度]
遷移演算子を座標の一次まで取る。最初の $B_1$ と初期の量子数 $n$ を使う一成分の計算である。

量の定義に従い、共通の単位で次の値を用いる。
$$
n=1,\qquad B_{1}=50\,\mathrm{e^2\,fm^4}
$$

この条件で求める量を ${B}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{B}&=B_{1} \cdot n\\&=\left(1\right) \cdot \left(50\right)\\&=50\,\mathrm{e^2\,fm^4}.\end{aligned}
$$

角運動量の異なる終成分がある場合は分岐に合成係数が入る。

三量子からの下降を比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
n=3
$$

このときの値を ${B}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{B}_{\mathrm{new}}&=\left(3\right) \cdot \left(50\right)\\&=150\,\mathrm{e^2\,fm^4},\\{B}_{\mathrm{new}}-{B}&=100\,\mathrm{e^2\,fm^4}.\end{aligned}
$$
:::

:::example[非調和間隔]
有効模型 $E_n=nW+g n(n-1)/2$ で二量子のずれを、$E_2/E_1$ で見る。

量の定義に従い、共通の単位で次の値を用いる。
$$
W=\frac{4}{5}\,\mathrm{MeV},\qquad g=0\,\mathrm{MeV}
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=2 + \frac{g}{W}\\&=\frac{\left(0\right)}{\left(\frac{4}{5}\right)} + 2\\&=2.\end{aligned}
$$

振動模型で比 2 からの変化は有効な量子間相互作用の尺度になる。

弱い斥力を加える。そこで次の値に変更し、ほかの入力を固定する。
$$
g=\frac{3}{20}\,\mathrm{MeV}
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=\frac{\left(\frac{3}{20}\right)}{\left(\frac{4}{5}\right)} + 2\\&\simeq 2.1875,\\{R}_{\mathrm{new}}-{R}&\simeq 0.1875.\end{aligned}
$$
:::

## 乱雑位相近似と集団励起

### 粒子穴振幅

:term[乱雑位相近似]{en="random-phase approximation" yomi="らんざついそうきんじ"} は静的な平均場の周りの小振動を扱う。占有軌道 $h$ と非占有軌道 $p$ の対 $b_{ph}^\dagger=a_p^\dagger a_h$ を用いる。反交換関係で $[a_h^\dagger a_p,a_{p'}^\dagger a_{h'}]=\delta_{pp'}a_h^\dagger a_{h'}-\delta_{hh'}a_{p'}^\dagger a_p$ であり、Slater 行列式上の期待値は $\delta_{pp'}\delta_{hh'}$ になる。この期待値で交換子を置き換え、対を近似的なボソンとして扱う。

実数の対称行列 $A,B$ を使い、二次 Hamiltonian を $b^\dagger A b+\tfrac12(b^\dagger B b^\dagger+b B b)$ とする。$Q^\dagger=\sum_i(X_ib_i^\dagger-Y_ib_i)$ と書いて $[H,Q^\dagger]=\hbar\omega Q^\dagger$ の作成・消滅成分を比べると
$$
\begin{pmatrix}A&B\\-B&-A\end{pmatrix}\binom X Y
=\hbar\omega\binom X Y,\qquad X^{\mathsf T}X-Y^{\mathsf T}Y=1
\label{eq-nuc-rpa}
$$
となる。後式は $[Q,Q^\dagger]=1$ から従う。$B=0$ では励起側の振幅 $X$ だけで解ける。$B$ は後退振幅 $Y$ を生み、復元力を変える。

:::example[一自由度の安定性]
$A=5$、$B=3\,\mathrm{MeV}$ の固有値は $\pm\sqrt{A^2-B^2}=\pm4\,\mathrm{MeV}$ である。$Q$ の正規化から正エネルギー枝は $X^2=(A/\hbar\omega+1)/2=9/8$、$Y^2=1/8$ になる。$|B|>A>0$ なら虚の周波数が出て、展開点の不安定性を示す。
:::

:::example[同相の粒子穴]
$A=\epsilon I+g\vb f\vb f^{\mathsf T}$、$B=0$、$\vb f=(1,1)^{\mathsf T}$ とする。$\vb f$ に平行な状態は $\epsilon+2g$、直交状態は $\epsilon$ である。外場も $\vb f$ に沿うなら、強度は平行な状態にだけ集まる。
:::

### 小振動の二つの行列

:eqref[eq-nuc-rpa] の上段と下段を足し引きすると
$$
(A+B)(X+Y)=\hbar\omega(X-Y),\qquad
(A-B)(X-Y)=\hbar\omega(X+Y)
$$
となる。$A-B>0$ とすれば、$X+Y=(A-B)^{1/2}z$ と置いて
$$
(A-B)^{1/2}(A+B)(A-B)^{1/2}z=(\hbar\omega)^2 z
$$
を得る。左辺は実対称であり、さらに $A+B>0$ なら正定値である。従って全ての $\omega^2$ は正である。$A-B$ に零固有値があれば、この変換に必要な逆平方根は定義できない。自己無撞着な平均場で並進や回転の対称性を自発的に破る場合、対応する方向の復元力がゼロになる。物理的な密度振動を調べるには、この対称変換に沿う零モードを分ける。

:::example[復元力の消失]
一自由度で $A=B>0$ に近づけると $A-B\to0$、$A+B\to2A$ で $\omega\to0$ になる。$A$ 単独の粒子穴エネルギーは有限でも、小振動の周波数はゼロへ近づく。
:::

### RPA の規格化と安定性

一成分で実数 $A,B$ とし、$A>|B|$ とする。$Q^\dagger=Xb^\dagger-Yb$ を用いると交換関係は $[Q,Q^\dagger]=X^2-Y^2$ である。固有値方程式の行列式は $\Omega^2=A^2-B^2$、正枝を $\Omega>0$ として比は $Y/X=-B/(A+\Omega)$ である。規格化条件から $X^2=(A/\Omega+1)/2$、$Y^2=(A/\Omega-1)/2$、$XY=-B/(2\Omega)$ が従う。従って $B\to A$ で正の固有値はゼロに近づき、後ろ向き振幅の大きさが増す。もとの平均場が不安定なら小振幅近似の振動数を実励起として解釈できない。

:::example[後ろ向き振幅]
一成分の正枝の $Y^2$ を求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
A=5\,\mathrm{MeV},\qquad B=0\,\mathrm{MeV}
$$

この条件で求める量を ${Y^2}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{Y^2}&=\frac{\frac{A}{\sqrt{A^{2} - B^{2}}} - 1}{2}\\&=\frac{\left(5\right)}{2 \cdot \sqrt{- \left(0\right)^{2} + \left(5\right)^{2}}} - \frac{1}{2}\\&=0.\end{aligned}
$$

粒子穴の前向き近似ではこの成分はゼロである。

結合を入れる。そこで次の値に変更し、ほかの入力を固定する。
$$
B=4\,\mathrm{MeV}
$$

このときの値を ${Y^2}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{Y^2}_{\mathrm{new}}&=\frac{\left(5\right)}{2 \cdot \sqrt{- \left(4\right)^{2} + \left(5\right)^{2}}} - \frac{1}{2}\\&\simeq 0.333333,\\{Y^2}_{\mathrm{new}}-{Y^2}&\simeq 0.333333.\end{aligned}
$$
:::

:::example[座標型演算子の強度]
$F=f(b^\dagger+b)$ とし、$[Q,F]=f(X+Y)$ を使う。

量の定義に従い、共通の単位で次の値を用いる。
$$
f=2,\qquad A=5\,\mathrm{MeV},\qquad B=-3\,\mathrm{MeV}
$$

この条件で求める量を ${B_F}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{B_F}&=\frac{f^{2} \cdot \left(A - B\right)}{\sqrt{A^{2} - B^{2}}}\\&=\frac{\left(2\right)^{2} \cdot \left(- \left(-3\right) + \left(5\right)\right)}{\sqrt{- \left(-3\right)^{2} + \left(5\right)^{2}}}\\&=8.\end{aligned}
$$

固有値は $B$ の符号を区別しないが、この演算子に対する強度は異なる。

結合の符号を反転する。そこで次の値に変更し、ほかの入力を固定する。
$$
B=3\,\mathrm{MeV}
$$

このときの値を ${B_F}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{B_F}_{\mathrm{new}}&=\frac{\left(2\right)^{2} \cdot \left(- \left(3\right) + \left(5\right)\right)}{\sqrt{- \left(3\right)^{2} + \left(5\right)^{2}}}\\&=2,\\{B_F}_{\mathrm{new}}-{B_F}&=-6.\end{aligned}
$$
:::

:::example[安定限界に近づく振動]
座標と運動量の係数の積から正枝を求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
A=8\,\mathrm{MeV},\qquad B=4\,\mathrm{MeV}
$$

この条件で求める量を ${\Omega}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\Omega}&=\sqrt{\left(A - B\right) \cdot \left(A + B\right)}\\&=\sqrt{\left(- \left(4\right) + \left(8\right)\right) \cdot \left(\left(4\right) + \left(8\right)\right)}\\&\simeq 6.9282\,\mathrm{MeV}.\end{aligned}
$$

$A-B$ が小さくなると固有値が急に下がる。

不安定境界の手前まで移す。そこで次の値に変更し、ほかの入力を固定する。
$$
B=\frac{39}{5}\,\mathrm{MeV}
$$

このときの値を ${\Omega}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\Omega}_{\mathrm{new}}&=\sqrt{\left(\left(8\right) - \left(\frac{39}{5}\right)\right) \cdot \left(\left(8\right) + \left(\frac{39}{5}\right)\right)}\\&\simeq 1.77764\,\mathrm{MeV},\\{\Omega}_{\mathrm{new}}-{\Omega}&\simeq -5.15056\,\mathrm{MeV}.\end{aligned}
$$

**一次近似の検証。** 不安定性の近くでは線形化の誤差が増える。十分離れた点での傾きを境界近くまで延ばすとどうなるかを比べる。

${B}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${\Omega}$ の一階微分を求めると、
$$
\frac{\mathrm d{\Omega}}{\mathrm dB}=- \frac{B}{\sqrt{A^{2} - B^{2}}}.
$$

この導関数を変更前の点で評価して $h=\frac{19}{5}$ を掛ける。ここで $B$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq -2.19393\,\mathrm{MeV},\\\Delta-\Delta_{\mathrm{lin}}&\simeq -2.95663\,\mathrm{MeV}.\end{aligned}
$$

省いた分が全変化に比べて大きい。従って増分の大きさを予測するには、少なくとも二次の項を調べるか、元の式を使う必要がある。
:::

:::example[静的な応答]
外場 $-\lambda F$ に対し、励起エネルギー $\Omega$ の一状態が強度 $b=|\langle n|F|0\rangle|^2$ を持つとする。二次のエネルギー変化は $-\lambda^2 b/\Omega$ である。これを $-\chi\lambda^2/2$ と書いたときの静的感受率 $\chi$ を求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
b=4,\qquad \Omega=8\,\mathrm{MeV}
$$

この条件で求める量を ${\chi}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\chi}&=\frac{2 \cdot b}{\Omega}\\&=\frac{2 \cdot \left(4\right)}{\left(8\right)}\\&=1\,\mathrm{MeV}^{-1}.\end{aligned}
$$

同じ強度で低い励起ほど静的な変形に寄与する。

軟化した状態と比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
\Omega=2\,\mathrm{MeV}
$$

このときの値を ${\chi}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\chi}_{\mathrm{new}}&=\frac{2 \cdot \left(4\right)}{\left(2\right)}\\&=4\,\mathrm{MeV}^{-1},\\{\chi}_{\mathrm{new}}-{\chi}&=3\,\mathrm{MeV}^{-1}.\end{aligned}
$$
:::

### 時間依存の Hartree–Fock から粒子穴方程式へ

Slater 状態の一体密度 $\rho(t)$ は、各軌道に $i\hbar\dot\phi_i=h[\rho]\phi_i$ を課すと $i\hbar\dot\rho=[h[\rho],\rho]$ に従う。実際、$\rho=\sum_i\ket{\phi_i}\bra{\phi_i}$ の微分で右と左の作用が逆符号になる。停留解 $\rho_0$ の周りで $\rho=\rho_0+\delta\rho$、$h=h_0+\delta h$ とすれば、一次の式は
$$
i\hbar\delta\dot\rho=[h_0,\delta\rho]+[\delta h,\rho_0],\qquad
\delta h = \left.\frac{\delta h}{\delta\rho}\right|_0\delta\rho
$$
である。射影条件 $\rho^2=\rho$ の一次は $\rho_0\delta\rho+\delta\rho\rho_0=\delta\rho$ であり、占有・占有と非占有・非占有のブロックはゼロになる。そこで $h_0$ と $\rho_0$ の共通固有基底で $[\delta h,\rho_0]_{ph}=\delta h_{ph}$、$[\delta h,\rho_0]_{hp}=-\delta h_{hp}$ を用いる。密度の振動を
$$
\delta\rho_{ph}(t)=X_{ph}e^{-i\omega t}+Y_{ph}^*e^{i\omega t},\qquad
\delta\rho_{hp}(t)=Y_{ph}e^{-i\omega t}+X_{ph}^*e^{i\omega t}
$$
と置けば Hermite 性が保たれる。$e^{-i\omega t}$ の係数を比べると、上段に $(\epsilon_p-\epsilon_h)X_{ph}+\delta h_{ph}[X,Y]$ が現れ、下段に $-(\epsilon_p-\epsilon_h)Y_{ph}-\delta h_{hp}[X,Y]$ が現れる。実の基底で二階のエネルギー変分から交換対称な行列を組めば :eqref[eq-nuc-rpa] の二つのブロック行列となる。誘起場をゼロにした場合は独立した粒子穴準位へ戻る。残留相互作用は平均場と同じエネルギーを微分したものを用いると、場の対称性を線形化の段階でも保てる。

並進のようにエネルギーを保つ連続変換を $U(s)$ とする。$\rho_s=U(s)\rho_0U(s)^\dagger$ が停留解の族なら $[h[\rho_s],\rho_s]=0$ を $s$ で微分し、小振幅方程式の右辺が $\delta\rho=(\partial_s\rho_s)_{s=0}$ に対してゼロとなる。そのため平均場が対称性を破って $\delta\rho\ne0$ であっても、これは励起エネルギーを持つ内部の振動とは区別される。有限の基底や一部の残留項を捨てた計算では対応する振動数が零からずれることがあり、自己無撞着性の診断に使える。

### 二体力から RPA の行列を求める

密度に依存しない二体相互作用では、反対称化した行列要素を $\bar v_{ab,cd}=\langle ab|v|cd\rangle-\langle ab|v|dc\rangle$ とすれば、:eqref[eq-nuc-hartree-fock] の変分は $\delta h_{ab}=\sum_{cd}\bar v_{ac,bd}\delta\rho_{dc}$ である。上の振動を入れると、$X_{p'h'}$ には $c=h',d=p'$、$Y_{p'h'}$ には $c=p',d=h'$ の項が掛かる。実基底では従って
$$
A_{ph,p'h'}=(\epsilon_p-\epsilon_h)\delta_{pp'}\delta_{hh'}+
\bar v_{ph',hp'},\qquad B_{ph,p'h'}=\bar v_{pp',hh'}
$$
となる。Hermite 性と反対称化により両行列は対称である。密度依存の汎関数を使う場合は、結合係数の密度微分も $\delta h$ に含める。

## 巨大共鳴

多くの粒子穴励起が同相に組み合わさり、外場への強い応答を狭いエネルギー域に作る。これを :term[巨大共鳴]{en="giant resonance" yomi="きょだいきょうめい"} という。等スカラー単極子では密度が膨張と圧縮を繰り返し、核物質の非圧縮率と関係する。等ベクトル双極子では中性子と陽子の分布が逆向きにずれ、全体の並進から区別する。電気双極子の並進不変な $z$ 成分は $D=e[(N/A)\sum_p z_p-(Z/A)\sum_n z_n]$ である。

Hermite 演算子 $F$ に対し、完全なエネルギー固有基底を挿入すると
$$
m_1=\sum_n(E_n-E_0)|\langle n|F|0\rangle|^2
=\frac12\bra0[F,[H,F]]\ket0
\label{eq-nuc-energy-weighted-sum}
$$
となる。実際 $[F,[H,F]]=2FHF-F^2H-HF^2$ である。等質量核子で $D$ と相互作用が可換な模型では、運動項の $[z,[p_z^2/(2m),z]]=\hbar^2/m$ を足して $m_1(D)=e^2\hbar^2 NZ/(2mA)$ となる。速度依存や交換電流の効果がある場合は修正が生じる。

### 二流体模型と双極子和則

陽子と中性子をそれぞれ剛な流体として、$z$ 方向の中心を $x_p,x_n$ とする。全重心を静止させ、相対変位を $x=x_p-x_n$ とすれば $x_p=Nx/A$、$x_n=-Zx/A$ である。従って運動エネルギー $m(Z\dot x_p^2+N\dot x_n^2)/2$ は $\mu\dot x^2/2$、$\mu=mNZ/A$ となる。復元力を $Cx$ とする模型では $\omega=\sqrt{C/\mu}$ であり、双極子演算子は $D=e(NZ/A)x$ になる。

調和振動子で $\langle1|x|0\rangle=\sqrt{\hbar/(2\mu\omega)}$ だから、この一状態のエネルギー重み付き強度は $e^2(NZ/A)^2\hbar^2/(2\mu)$ となり、上の一成分の双極子和則に一致する。これは模型内で強度が一つのモードに集まる結果である。実際には連続状態への粒子放出や複雑な配位への混合で幅が生じる。軸変形があれば軸に平行な変位と垂直な変位で復元力が異なり、二つの共鳴成分が現れ得る。

:::example[強度重心]
強度 2 と 1 の二つのピークが $10$ と $16\,\mathrm{MeV}$ にあるなら、全強度 $m_0=3$、エネルギー重み付きの和 $m_1=36\,\mathrm{MeV}$ に強度単位を掛けた値である。重心 $m_1/m_0$ は $12\,\mathrm{MeV}$ である。
:::

:::example[並進を取り除く]
$z_i\mapsto z_i+a$ とすると $D$ の変化は $ea[(N/A)Z-(Z/A)N]=0$ である。陽子だけの $\sum_p ez_p$ では $eZa$ が残り、内部励起と重心移動を区別できない。
:::

### 静的な分極率と低エネルギー強度

外場に $H'=H-\lambda F$ を加えて、非縮退基底を一次の摂動論で変える。無摂動の基底を $\ket0$ とすると、一次の状態の増分は $\lambda\sum_{n\ne0}\ket n\mel nF0/(E_n-E_0)$ である。従って $\langle F\rangle$ の変化は $2\lambda\sum_{n\ne0}|\mel nF0|^2/(E_n-E_0)$ になる。この和を逆エネルギー重み付き和 $m_{-1}$ と呼ぶ。:eqref[eq-nuc-energy-weighted-sum] の $m_1$ は高い準位に大きな重みをつけるが、静的分極率 $2m_{-1}$ は低エネルギー強度にも敏感である。

モーメント $m_k=\sum_n E_{n0}^k s_n$、$E_{n0}>0$、$s_n=|\mel nF0|^2$ と書く。Cauchy–Schwarz の不等式を成分 $\sqrt{E_{n0}s_n}$ と $\sqrt{s_n/E_{n0}}$ に使うと $m_0^2\le m_1m_{-1}$ である。従って $m_0>0$ なら $\sqrt{m_1/m_{-1}}\le m_1/m_0$ となる。等号は、強度がある全ての状態の $E_{n0}$ が同じ場合に限られる。

:::example[少量の低励起成分]
$10\,\mathrm{MeV}$ に強度 10 のピークだけがあると $m_1=100$、$m_{-1}=1$ である。$1\,\mathrm{MeV}$ に強度 1 を加えると $m_1$ は 101 へ少し増える一方、$m_{-1}$ は 2 に倍増する。和則を比較すると、どの観測が強度の低エネルギー部分を反映するかを判断できる。
:::

### 和則からみる重心と広がり

強度分布を正の重み $B_i$ の和とし、$m_k=\sum_i E_i^k B_i$ とする。平均 $\bar E=m_1/m_0$ と分散 $m_2/m_0-\bar E^2$ は、重みを $B_i/m_0$ に規格化して二乗を展開すれば得られる。逆エネルギー重みを使う推定 $\sqrt{m_1/m_{-1}}$ は $m_1/m_0$ 以下である。実際 $m_1m_{-1}\ge m_0^2$ は $\sqrt{B_i E_i}$ と $\sqrt{B_i/E_i}$ の Cauchy–Schwarz 不等式である。等号は非零強度の全エネルギーが同一の場合であり、低励起に少量の強度が移ると比の差が目立つ。

:::example[二ピークの重心]
低いエネルギー $E_l$ の割合を $p$ とする。残りは $E_h$ に集中させる。

量の定義に従い、共通の単位で次の値を用いる。
$$
p=\frac{1}{10},\qquad E_{l}=5\,\mathrm{MeV},\qquad E_{h}=15\,\mathrm{MeV}
$$

この条件で求める量を ${\bar E}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\bar E}&=E_{h} \cdot \left(1 - p\right) + E_{l} \cdot p\\&=\left(15\right) \cdot \left(1 - \left(\frac{1}{10}\right)\right) + \left(5\right) \cdot \left(\frac{1}{10}\right)\\&=14\,\mathrm{MeV}.\end{aligned}
$$

低線に強度が移ると重心が低下する。

低いピークの重みを増す。そこで次の値に変更し、ほかの入力を固定する。
$$
p=\frac{3}{10}
$$

このときの値を ${\bar E}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\bar E}_{\mathrm{new}}&=\left(15\right) \cdot \left(1 - \left(\frac{3}{10}\right)\right) + \left(5\right) \cdot \left(\frac{3}{10}\right)\\&=12\,\mathrm{MeV},\\{\bar E}_{\mathrm{new}}-{\bar E}&=-2\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[逆エネルギーの重み]
同じ二ピークでエネルギーの逆数の平均を求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
p=\frac{1}{10},\qquad E_{l}=5\,\mathrm{MeV},\qquad E_{h}=15\,\mathrm{MeV}
$$

この条件で求める量を ${m_{-1}/m_0}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{m_{-1}/m_0}&=\frac{p}{E_{l}} + \frac{1 - p}{E_{h}}\\&=\frac{\left(\frac{1}{10}\right)}{\left(5\right)} + \frac{1 - \left(\frac{1}{10}\right)}{\left(15\right)}\\&\simeq 0.08\,\mathrm{MeV}^{-1}.\end{aligned}
$$

小さな強度でも低いエネルギーにあれば応答を大きくする。

前例と同じ重みの変化を比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
p=\frac{3}{10}
$$

このときの値を ${m_{-1}/m_0}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{m_{-1}/m_0}_{\mathrm{new}}&=\frac{\left(\frac{3}{10}\right)}{\left(5\right)} + \frac{1 - \left(\frac{3}{10}\right)}{\left(15\right)}\\&\simeq 0.106667\,\mathrm{MeV}^{-1},\\{m_{-1}/m_0}_{\mathrm{new}}-{m_{-1}/m_0}&\simeq 0.0266667\,\mathrm{MeV}^{-1}.\end{aligned}
$$

**一次近似の検証。** 低い線へわずかな強度が移るときの感度を直接計算する。

${p}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${m_{-1}/m_0}$ の一階微分を求めると、
$$
\frac{\mathrm d{m_{-1}/m_0}}{\mathrm dp}=\frac{E_{h} - E_{l}}{E_{h} \cdot E_{l}}.
$$

この導関数を変更前の点で評価して $h=\frac{1}{5}$ を掛ける。ここで $p$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq 0.0266667\,\mathrm{MeV}^{-1},\\\Delta-\Delta_{\mathrm{lin}}&=0\,\mathrm{MeV}^{-1}.\end{aligned}
$$

入力のこの変化について差は一次の項だけで得られる。従って基準点での傾きで増分を求められる。
:::

:::example[並進を消した双極子和則]
核子の質量を $m$ とし、電荷を $e$ で割った双極子の一 Cartesian 成分を用いる。各核子の係数は陽子が $N/A$、中性子が $-Z/A$ なので二乗和は $NZ/A$ である。

量の定義に従い、共通の単位で次の値を用いる。
$$
h_{c}=\frac{197327}{1000}\,\mathrm{MeV\,fm},\qquad m_{c}=939\,\mathrm{MeV},\qquad Z=20,\qquad N=20
$$

この条件で求める量を ${m_1}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{m_1}&=\frac{N \cdot Z \cdot h_{c}^{2}}{2 \cdot m_{c} \cdot \left(N + Z\right)}\\&=\frac{\left(20\right) \cdot \left(\frac{197327}{1000}\right)^{2}}{4 \cdot \left(939\right)}\\&\simeq 207.337\,\mathrm{MeV\,fm^2}.\end{aligned}
$$

座標のみの相互作用を仮定する。速度依存の相互作用はこの運動項の和則を変えうる。

同位体で比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
N=28
$$

このときの値を ${m_1}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{m_1}_{\mathrm{new}}&=\frac{\left(20\right) \cdot \left(28\right) \cdot \left(\frac{197327}{1000}\right)^{2}}{2 \cdot \left(939\right) \cdot \left(\left(20\right) + \left(28\right)\right)}\\&\simeq 241.894\,\mathrm{MeV\,fm^2},\\{m_1}_{\mathrm{new}}-{m_1}&\simeq 34.5562\,\mathrm{MeV\,fm^2}.\end{aligned}
$$
:::

:::example[寿命幅の換算]
単一の指数減衰の Lorentz 形の全幅を $\Gamma$ とする。平均寿命は $\hbar/\Gamma$ である。

量の定義に従い、共通の単位で次の値を用いる。
$$
\hbar=\frac{164553}{250000000000000000000000000}\,\mathrm{MeV\,s},\qquad \Gamma=3\,\mathrm{MeV}
$$

この条件で求める量を ${\tau}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\tau}&=\frac{\hbar}{\Gamma}\\&=\frac{\left(\frac{164553}{250000000000000000000000000}\right)}{\left(3\right)}\\&\simeq 2.19404\times10^{-22}\,\mathrm{s}.\end{aligned}
$$

分布の広がりすべてが単一寿命による場合の換算である。複数状態への分裂は区別する。

幅の大きい単一線と比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
\Gamma=6\,\mathrm{MeV}
$$

このときの値を ${\tau}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\tau}_{\mathrm{new}}&=\frac{\left(\frac{164553}{250000000000000000000000000}\right)}{\left(6\right)}\\&\simeq 1.09702\times10^{-22}\,\mathrm{s},\\{\tau}_{\mathrm{new}}-{\tau}&\simeq -1.09702\times10^{-22}\,\mathrm{s}.\end{aligned}
$$
:::

### 圧縮の振動を密度の曲率に結び付ける

一様な基底密度の球の各点を $\vb r\mapsto(1+a)\vb r$ と伸縮する微小な単極子振動を考える。粒子数を一定にすると密度比は $(1+a)^{-3}=1-3a+O(a^2)$ である。無限物質での非圧縮率 $K_0$ の定義を一様球へ適用すると、一核子当たりの静的なエネルギー上昇は $K_0(3a)^2/18=K_0 a^2/2$ となる。流れの速度は $\dot a\,\vb r$ なので運動エネルギーは $Am\langle r^2\rangle\dot a^2/2$ である。従って角振動数は $\omega=\sqrt{K_0/(m\langle r^2\rangle)}$ になる。有限核の表面、電荷、非対称度を含む場合は、復元力を一様物質の $A K_0$ に置く近似を見直す。測った単極子共鳴のエネルギーから $m\langle r^2\rangle\omega^2$ を復元した量は、有限核全体の圧縮に対する応答であり、この補正を経て無限物質の状態方程式へつなげられる。

## 高スピン状態とバックベンディング

回転周波数を上げると、低周波数では対を作っていた二核子が回転軸に角運動量を揃えた配置が有利になることがある。回転座標系でのエネルギー $E'=E-\omega\langle J_x\rangle$ は整列した配置ほど早く下がる。二つの配置が交差や混合を起こす領域では、実験から求めた慣性能率が急変する。慣性能率を回転周波数に対して描いたときの曲がりを :term[バックベンディング]{en="backbending" yomi="ばっくべんでぃんぐ"} と呼ぶ。

バンド内の光子エネルギーを $E_\gamma(J)=E(J)-E(J-2)$ とすると、剛体回転模型では :eqref[eq-nuc-rotor] から $E_\gamma(J)=(2J-1)\hbar^2/\mathcal I$ である。$\mathcal I^{(1)}=\hbar^2(2J-1)/E_\gamma(J)$ として読み取れる。周波数の増分への応答をみる動的な値は $\mathcal I^{(2)}\simeq4\hbar^2/[E_\gamma(J+2)-E_\gamma(J)]$ である。これは $\hbar\omega\simeq E_\gamma/2$、$\Delta J=2$ の差分である。

:::example[一定の慣性能率]
$\mathcal I=20\hbar^2/\mathrm{MeV}$ なら隣り合う $\gamma$ 線の差は $4\hbar^2/\mathcal I=0.2\,\mathrm{MeV}$ で一定である。
:::

:::example[二配置の交差]
回転座標系で配置の差を $\Delta E'=\Delta E_0-\omega\Delta J$ とし、$\Delta E_0=2\,\mathrm{MeV}$、$\Delta J=8\hbar$ とする。混合しない交差は $\hbar\omega=0.25\,\mathrm{MeV}$ にある。非ゼロの結合があれば固有値の差は交差近傍で有限になる。
:::

### 回転座標系での二つの枝

$\hbar\omega=w$ をエネルギー単位、整列を $i$ を無次元とする。未混合二帯の回転座標系エネルギーを $e_g(w),e_a(w)$ とし、差を $d_0-w i$ と線形化する。$i>0$ なら交点は $w_c=d_0/i$ にある。一定の混合 $v$ を入れると差は $\sqrt{(d_0-w i)^2+4v^2}$ となり、交差は避けられる。下枝の $w$ 微分から整列を求めると、交点付近で上の帯の寄与が連続的に混ざる。近傍の幅は $|d_0-w i|\lesssim 2|v|$ で見積もれる。混合が小さいと急なバックベンディングになり、大きければ滑らかな変化になる。

:::example[回転周波数の交点]
未混合の静止差 $d_0$ と整列の差 $i$ を使う。

量の定義に従い、共通の単位で次の値を用いる。
$$
d_{0}=2\,\mathrm{MeV},\qquad i=8
$$

この条件で求める量を ${w_c}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{w_c}&=\frac{d_{0}}{i}\\&=\frac{\left(2\right)}{\left(8\right)}\\&\simeq 0.25\,\mathrm{MeV}.\end{aligned}
$$

整列が大きい枝ほど低い回転周波数で有利になる。

より大きな整列と比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
i=10
$$

このときの値を ${w_c}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{w_c}_{\mathrm{new}}&=\frac{\left(2\right)}{\left(10\right)}\\&\simeq 0.2\,\mathrm{MeV},\\{w_c}_{\mathrm{new}}-{w_c}&\simeq -0.05\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[混合による最小分離]
未混合の交点での最小準位差を求める。混合行列要素の絶対値を $v$ とする。

量の定義に従い、共通の単位で次の値を用いる。
$$
v=\frac{1}{20}\,\mathrm{MeV}
$$

この条件で求める量を ${D}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{D}&=2 \cdot v\\&=2 \cdot \left(\frac{1}{20}\right)\\&\simeq 0.1\,\mathrm{MeV}.\end{aligned}
$$

ここでは $v$ を正に取った。符号が負の場合も間隔は絶対値で決まる。

強い混合と比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
v=\frac{1}{5}\,\mathrm{MeV}
$$

このときの値を ${D}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{D}_{\mathrm{new}}&=2 \cdot \left(\frac{1}{5}\right)\\&\simeq 0.4\,\mathrm{MeV},\\{D}_{\mathrm{new}}-{D}&\simeq 0.3\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[避けられた交差の周波数幅]
差の絶対値が $2v$ 以下になる区間の長さを求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
v=\frac{1}{20}\,\mathrm{MeV},\qquad i=8
$$

この条件で求める量を ${\Delta w}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\Delta w}&=\frac{4 \cdot v}{i}\\&=\frac{4 \cdot \left(\frac{1}{20}\right)}{\left(8\right)}\\&\simeq 0.025\,\mathrm{MeV}.\end{aligned}
$$

これは混合を分解するための幅の目安である。実験での周波数は遷移エネルギーの差分から近似する。

前例と同じ混合の増加を使う。そこで次の値に変更し、ほかの入力を固定する。
$$
v=\frac{1}{5}\,\mathrm{MeV}
$$

このときの値を ${\Delta w}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\Delta w}_{\mathrm{new}}&=\frac{4 \cdot \left(\frac{1}{5}\right)}{\left(8\right)}\\&\simeq 0.1\,\mathrm{MeV},\\{\Delta w}_{\mathrm{new}}-{\Delta w}&\simeq 0.075\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[動的慣性能率の敏感さ]
二本の隣接する下降遷移エネルギーの差 $d$ を用いて、有限差分の $\mathcal I^{(2)}$ を求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
d=\frac{3}{25}\,\mathrm{MeV}
$$

この条件で求める量を ${\mathcal I^{(2)}/\hbar^2}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\mathcal I^{(2)}/\hbar^2}&=\frac{4}{d}\\&=\frac{4}{\left(\frac{3}{25}\right)}\\&\simeq 33.3333\,\mathrm{MeV}^{-1}.\end{aligned}
$$

隣接線が接近すると推定値が増え、測定誤差の影響も増す。

隣接差の小さい区間と比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
d=\frac{1}{25}\,\mathrm{MeV}
$$

このときの値を ${\mathcal I^{(2)}/\hbar^2}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\mathcal I^{(2)}/\hbar^2}_{\mathrm{new}}&=\frac{4}{\left(\frac{1}{25}\right)}\\&=100\,\mathrm{MeV}^{-1},\\{\mathcal I^{(2)}/\hbar^2}_{\mathrm{new}}-{\mathcal I^{(2)}/\hbar^2}&\simeq 66.6667\,\mathrm{MeV}^{-1}.\end{aligned}
$$

**一次近似の検証。** 隣接する線の差の誤差を評価するとき、差の逆数を線形化できる範囲を見極める必要がある。

${d}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${\mathcal I^{(2)}/\hbar^2}$ の一階微分を求めると、
$$
\frac{\mathrm d{\mathcal I^{(2)}/\hbar^2}}{\mathrm dd}=- \frac{4}{d^{2}}.
$$

この導関数を変更前の点で評価して $h=- \frac{2}{25}$ を掛ける。ここで $d$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq 22.2222\,\mathrm{MeV}^{-1},\\\Delta-\Delta_{\mathrm{lin}}&\simeq 44.4444\,\mathrm{MeV}^{-1}.\end{aligned}
$$

省いた分が全変化に比べて大きい。従って増分の大きさを予測するには、少なくとも二次の項を調べるか、元の式を使う必要がある。
:::

### 混合した帯の整列を定量化する

回転座標系の下側の固有値を $e_-=(e_g+e_a-D)/2$、$D=\sqrt{d(w)^2+4v^2}$ とする。$d=e_a-e_g=d_0-iw$、$i_g=-\partial_we_g$ と置く。$v$ が定数なら
$$
-\partial_we_-=i_g+\frac i2\left(1-\frac{d}{D}\right),\qquad
-\partial_w^2e_-=-\partial_w^2e_g+\frac{2i^2v^2}{D^3}
$$
である。第一式の追加分は未混合の整列帯の重みに $i$ を掛けたものになり、第二式は交差点で最も大きい。$v>0$ なら追加の曲率の最大は $i^2/(4v)$ である。このため弱い混合では整列が急に進む。帯の勾配が異なり、混合が有限なので、この近傍を :term[避けられた交差]{en="avoided crossing" yomi="さけられたこうさ"} と呼ぶ。実験の離散的なスピンから有限差分で慣性能率を出すと、この曲率の最大を直接標本化できない場合もある。変形や複数の整列軌道が同時に変わると一つの $i,v$ による説明は限られる。

## 演習問題

::::exercise{level=basic}
$K=0$ 回転バンドの $E(6)/E(2)$ を求めよ。
:::solution
:eqref[eq-nuc-rotor] により比は $6\cdot7/(2\cdot3)=7$ である。
:::
::::

::::exercise{level=basic}
八重極の一フォノンと、同じフォノン二つのパリティを求めよ。
:::solution
各フォノンは $(-1)^3=-1$ を持つので、一フォノンは負、二フォノンは積が正である。
:::
::::

::::exercise{level=basic}
$\hbar\omega=0.5\,\mathrm{MeV}$ の調和振動で、基底から三フォノン励起までのエネルギー差はいくらか。
:::solution
各量子が同じ $\hbar\omega$ を加える。零点エネルギーは差で消え、$3\cdot0.5=1.5\,\mathrm{MeV}$ である。
:::
::::

::::exercise{level=basic}
一成分の RPA で $A=13,B=5\,\mathrm{MeV}$ の正の固有エネルギーを求めよ。
:::solution
固有方程式は $\lambda^2=A^2-B^2=169-25=144\,\mathrm{MeV}^2$ なので、正の根は $12\,\mathrm{MeV}$ である。
:::
::::

::::exercise{level=advanced}
球面平均ゼロの実関数 $u$ について、$R=R_0(1+u+c)$ とする。体積を二次まで一定にする一様補正 $c=O(u^2)$ を求めよ。
:::solution
$V=\tfrac13\int R^3\dd\Omega$ を展開すると二次の差は $R_0^3\int(c+u^2)\dd\Omega$ になる。一次の $\int u$ はゼロなので、$c=-\int u^2\dd\Omega/(4\pi)$ と選べばよい。
:::
::::

::::exercise{level=advanced}
同じ八重極フォノンを二つ作ったとき許される $J$ をすべて求めよ。
:::solution
二つの $j=3$ の合成で $0\le J\le6$ を取る。交換位相は $(-1)^{6-J}$ であり、ボソンの交換対称性から偶数だけが残る。従って $J=0,2,4,6$ である。
:::
::::

::::exercise{level=advanced}
$F=\sum_i f_i z_i$、共通質量 $m$、$[V,F]=0$ でエネルギー重み付き和を求めよ。
:::solution
異なる粒子の演算子は可換なので $[F,[H,F]]=\sum_i f_i^2[z_i,[p_i^2/(2m),z_i]]=\hbar^2\sum_i f_i^2/m$ である。:eqref[eq-nuc-energy-weighted-sum] によって $m_1=\hbar^2\sum_i f_i^2/(2m)$ を得る。
:::
::::

::::exercise{level=advanced}
$E_\gamma(J)=0.30\,\mathrm{MeV}$、$E_\gamma(J+2)=0.40\,\mathrm{MeV}$ の区間の動的慣性能率を求めよ。
:::solution
差分 $0.10\,\mathrm{MeV}$ が $\Delta(\hbar\omega)$ の 2 倍である。よって $\mathcal I^{(2)}=4\hbar^2/(0.10\,\mathrm{MeV})=40\hbar^2/\mathrm{MeV}$ となる。
:::
::::

## 確認テスト

::::quiz{answer="3.333333" tol="1e-5"}
剛体回転の $K=0$ バンドで $E(4)/E(2)$ を求めよ。
:::explanation
$J(J+1)$ がそれぞれ 20 と 6 なので、比は $10/3$ である。
:::
::::

::::quiz
軸対称な変形の下で単粒子 Hamiltonian が保存する量を選べ。

- [x] 対称軸への全角運動量の射影
- [ ] 一般の方向への全角運動量の射影
- [ ] 任意の三次元回転に関する球対称性

:::explanation
軸の周りの回転にだけ対称性が残るため、その生成子である軸射影が保存される。別の方向に回すと変形ポテンシャルが変わる。
:::
::::

::::quiz{answer="4" tol="1e-6"}
二つの四重極フォノンからなる状態の最大 $J$ を答えよ。
:::explanation
角運動量の合成の上限は $2+2=4$ であり、交換位相も $(-1)^{4-4}=1$ なので許される。
:::
::::

::::quiz
巨大双極子共鳴の内部運動として適切なものを選べ。

- [ ] 核全体の等速並進
- [x] 陽子分布と中性子分布の逆向きの振動
- [ ] 電荷密度の一様な静的増加

:::explanation
双極子の等ベクトル演算子は陽子と中性子に逆符号の重みをつける。重みの和がゼロなので一様な並進を除き、相対運動に応答する。
:::
::::

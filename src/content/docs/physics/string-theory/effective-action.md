---
title: "弦の低エネルギー有効作用"
description: "世界面の Weyl 不変性から時空の有効方程式を求め、ディラトンと温度の効果を調べる。"
notation:
  - symbol: "G_{\\mu\\nu},B_{\\mu\\nu},\\Phi"
    meaning: "時空計量、二形式、ディラトン"
  - symbol: "H"
    meaning: "三形式 \\dd B"
  - symbol: "G_E"
    meaning: "Einstein フレームの計量"
  - symbol: "\\beta_H"
    meaning: "Hagedorn 逆温度"
---

## この章の目標

- 背景場の世界面への結合を区別する。
- 計量の一ループの繰り込みと低エネルギー方程式を結び付ける。
- Weyl 変換で Einstein フレームへ移る。
- 線形ディラトンの中心電荷と弦の状態密度を求める。

## 背景場中の弦と非線形 σ 模型

平坦な世界面理論の質量ゼロの頂点を背景として作用に加える。Euclid 世界面で
$$
S=\frac1{4\pi\alpha'}\int\sqrt\gamma\,\dd^2\sigma\,
\bigl(\gamma^{ab}G_{\mu\nu}\partial_a X^\mu\partial_b X^\nu
+i\varepsilon^{ab}B_{\mu\nu}\partial_aX^\mu\partial_bX^\nu
+\alpha'R^{(2)}\Phi\bigr).
\label{eq-string-sigma-model}
$$
$\varepsilon^{ab}$ は $\varepsilon^{12}=1/\sqrt\gamma$ とした。$G,B,\Phi$ は $X$ に依存するので $X$ の作用は一般に二次式にならない。これを :term[非線形シグマ模型]{en="nonlinear sigma model" yomi="ひせんけいしぐまもけい"} と呼ぶ。$B\mapsto B+\dd\Lambda$ の差は Stokes の定理で境界積分になり、閉世界面では消える。例えば $H_{123}=h$ が定数となる局所表示に $B=hX^1\dd X^2\wedge\dd X^3$ がある。定数 $\Phi_0$ の項は Euler 標数に比例し、同じ種数での運動方程式を変えない。

### 二形式のゲージ変換と標的空間の座標変換

閉世界面 $\Sigma$ の作用の二形式項は引き戻しを使って $i(2\pi\alpha')^{-1}\int_\Sigma X^\ast B$ と書ける。$B\mapsto B+\dd\Lambda$ による変化は $i(2\pi\alpha')^{-1}\int_{\partial\Sigma}X^\ast\Lambda$ である。閉世界面ではゼロであり、境界があれば端の一形式の結合を同時に変換する。一方、標的座標の変換では $G$、$B$、$\Phi$ がそれぞれテンソル、二形式、スカラーとして変わる。作用全体は $X$ と背景の引き戻しを用いるため、この変換に不変である。

:::example[定数二形式が与える位相]
二次元トーラスの角 $\theta^1,\theta^2$ の周期を $2\pi$ とし、$B=\alpha'b\,\dd\theta^1\wedge\dd\theta^2$ とする。写像が各周期を $m,n$ 回巻けば $\int X^\ast B=(2\pi)^2\alpha' bmn$ なので、Euclid の経路積分の因子は $e^{-2\pi i bmn}$ である。$H=\dd B=0$ でも位相は一般に残り、$b$ の整数シフトだけは全巻き付きで見えない。
:::

:::example[局所座標で自由場に近づける]
二次元の標的計量 $\dd r^2+r^2\dd\theta^2$ は $r>0$ で Cartesian 座標 $x=r\cos\theta,y=r\sin\theta$ に移すと $\dd x^2+\dd y^2$ になる。例えば点 $r=2$ のまわりで角方向の長さの一次変化は $2\delta\theta$ である。極座標の $r$ 依存な係数だけを見て物理的な曲率があるとはいえない。
:::

### 質量ゼロの頂点を背景として読む

平坦な背景を $G_{\mu\nu}=\eta_{\mu\nu}+h_{\mu\nu}(X)$、$B_{\mu\nu}=b_{\mu\nu}(X)$ と変形する。線形の変形を Fourier 分解すると、各 $k$ の頂点は $\epsilon_{\mu\nu}\partial X^\mu\bar\partial X^\nu e^{ik\cdot X}$ の形である。正則自由場の OPE から指数の重みは $\alpha' k^2/4$ であり、微分から 1 が加わる。よって $k^2=0$ なら $(1,1)$ の重みを取る。$T(z)$ との OPE の三重極には $k^\mu\epsilon_{\mu\nu}$ が掛かるので、通常の一次の物理偏極には横条件が要る。左と右の微分を一つずつ含むため、$\epsilon$ の対称部は計量変形、反対称部は二形式変形に写る。トレースの物理的なスカラーを曲がった世界面で追うには、計量のトレースだけでなく $\Phi R^{(2)}$ の結合も必要である。

標的座標の無限小変換 $\delta X^\mu=-\xi^\mu(X)$ で背景を同時に動かす。一次で $\delta G_{\mu\nu}=\partial_\mu\xi_\nu+\partial_\nu\xi_\mu$ である。同様に $\delta B_{\mu\nu}=\partial_\mu\Lambda_\nu-\partial_\nu\Lambda_\mu$ の平面波は、$B$ の積分で全微分になる。従って $\epsilon$ に $k_\mu$ に比例する成分を加える自由度には標的のゲージ変換が対応する。これは線形の頂点で見える性質である。有限の変形では背景そのものが $X$ に依存するため、頂点を一回挿入して重みを数えるだけでは Weyl 不変性を保証できない。その判定が次節の $\bar\beta$ である。

:::example[平面波で微分展開の尺度を評価する]
時空のある空間方向に $h_{12}(X)=\epsilon\cos(kX^3)$ を置く。一次の Ricci テンソルでは $\partial_3^2 h_{12}$ が現れるので、その大きさは $|\epsilon| k^2$ に比例する。$k\ell_s=1/10$ なら一対の追加微分による相対因子は $\alpha' k^2=1/100$ である。場の振幅を小さくする条件 $|\epsilon|\ll1$ と、追加微分を小さくする条件 $|k|\ell_s\ll1$ は別である。
:::

### 世界面の計量を変分する

閉じた世界面で $\gamma_{ab}\mapsto e^{2\omega}\gamma_{ab}$ とする。二次元では $\sqrt\gamma\,\gamma^{ab}$ は不変であり、$\sqrt\gamma\,\varepsilon^{ab}$ も向きを固定すれば不変である。一方、曲率について
$$
\delta_\omega(\sqrt\gamma R^{(2)})=-2\sqrt\gamma\,\Delta_\gamma\omega
$$
である。従ってディラトン項の変化は $-(2\pi)^{-1}\int\sqrt\gamma\,\omega\,\Delta_\gamma\Phi$ になる。部分積分したので世界面が閉じていることを使った。開いた世界面では Euler 標数の式に合わせて境界の測地曲率を $\frac1{2\pi}\int_{\partial\Sigma} k_g \Phi$ として加える。一定 $\Phi_0$ に対して両者の和は $\Phi_0\chi$ であり、局所的な背景変分を生じない。$\Phi$ が $X$ に依存すれば、$\Delta_\gamma \Phi(X)$ が応力テンソルの改善を与える。

曲率との結合を $X$ で変分すると、曲がった世界面での方程式にも $R^{(2)}\nabla_\mu\Phi$ が入る。$B$ は $H=\dd B$ を通して寄与する。標的座標に沿う共変な世界面 Laplacian を $\mathcal D^a\partial_a X^\mu$ と書けば、
$$
\mathcal D^a\partial_a X^\mu
=\frac i2\varepsilon^{ab}H^\mu{}_{\nu\rho}\partial_a X^\nu\partial_b X^\rho
+\frac{\alpha'}2 R^{(2)}\nabla^\mu\Phi
$$
である。運動項の変分の係数 $-2$、二形式項の変分の係数 $+i$、曲率結合の変分の $+\alpha'$ を共通の $1/(4\pi\alpha')$ で割って得る。右辺の $i$ は Euclid 世界面の二形式の位相に由来する。Lorentz 世界面の古典運動方程式に戻すと、向きに対応する実の二形式結合になる。

:::example[定数ディラトンと線形ディラトンを境界で比較する]
平坦な長方形を周期的に貼ってトーラスにすると $\chi=0$ であり、定数 $\Phi_0$ の作用はゼロである。円板では面内の曲率をゼロに選んでも、境界の測地曲率の積分は $2\pi$ なので作用は $\Phi_0$ になる。線形の $\Phi$ では $\Delta_\gamma\Phi=V_\mu\Delta_\gamma X^\mu$ の式があるため、トポロジーの重みだけで応力テンソルの局所的な変化を判定できない。
:::

### $H$ の閉性と大域的な流束

$H$ は局所的な二形式ポテンシャルから $\dd B$ と書けるが、閉三サイクル上の流束を持つときには、同じ $B$ を標的の全体で取れない。閉世界面の位相をある三次元鎖 $W$ への延長で $\exp[-i(2\pi\alpha')^{-1}\int_W H]$ と定義する場合を考える。二通りの延長の差は閉三サイクル $C$ である。位相が延長の選び方に依存しないためには
$$
\frac1{(2\pi)^2\alpha'}\int_C H\in\Z
$$
が必要である。一般の写像で延長が存在するとは限らないので、その場合には局所ポテンシャルの貼り合わせを使う。この積分条件は、局所だけで $\dd H=0$ とする条件より強い。例えば面積形式に比例する二形式で局所的に $H=0$ でも、世界面に残るホロノミーがある。一方三サイクルの $H$ が非零なら、全体で一つの滑らかな $B$ を選ぶと Stokes の定理に反する。

:::example[三球での流束の次数]
半径 $L$ の $S^3$ で $H=h\,\mathrm{vol}_{S^3}$ なら $\int H=2\pi^2L^3h$ である。量子化された整数を $k$ とすれば $h=2k\alpha'/L^3$ である。一ループの局所方程式から $|h|=2/L$ が従うので、両者を合わせると最低次数で $L^2=|k|\alpha'$ となる。$|k|\gg1$ のとき曲率は弦長より十分緩やかであり、流束の密度も $|h|\ell_s=2/\sqrt{|k|}$ と小さい。
:::

## 世界面の β 関数と Einstein 方程式

### 背景のまわりの展開

$B=0$、$\Phi$ 一定から始める。$X=\exp_{\bar X}\xi$ とし、$\xi$ を $\bar X$ での接ベクトルとする。曲率の規約は球面の Ricci 曲率が正とする。第一変分をもう一度とり、共変微分の交換を使うと二次作用は $(4\pi\alpha')^{-1}\int(|D\xi|^2-\inner{R(\xi,\partial_a\bar X)\partial^a\bar X}{\xi})$ になる。境界項は閉世界面で落とした。薄い運動量の殻の自由縮約は $2\pi\alpha'\delta^{\mu\nu}\int_{\Lambda/b<|k|<\Lambda}\dd^2k/[(2\pi)^2 k^2]=\alpha'\log b\,\delta^{\mu\nu}$ である。二次の曲率項を縮約するとトレースが $R_{\mu\nu}\partial_a\bar X^\mu\partial^a\bar X^\nu$ となる。低い尺度に向けた $G$ の変化は $-\alpha'R_{\mu\nu}\log b$ なので $\mu\dv{G}{\mu}=\alpha'R$ である。

$H=\dd B$ と $\Phi$ もある場合、世界面の Weyl 応答の係数は、$\alpha'$ の一次までで
$$
\begin{aligned}
\bar\beta^G_{\mu\nu}/\alpha'&=R_{\mu\nu}-\tfrac14H_{\mu\rho\sigma}H_\nu{}^{\rho\sigma}+2\nabla_\mu\nabla_\nu\Phi,\\
\bar\beta^B_{\mu\nu}/\alpha'&=-\tfrac12\nabla^\rho H_{\rho\mu\nu}+\nabla^\rho\Phi H_{\rho\mu\nu}.
\end{aligned}
\label{eq-string-weyl-beta}
$$
ここでバーは Weyl 改善項を含むことを示す。二形式による二次作用と曲がった世界面での改善を、この節の下で導く。背景場の計算には :cite[polchinski-1] の規約を用いる。$H=0$ かつ $\Phi$ 一定では $R_{\mu\nu}=0$ となる。例えば曲率半径 $L$ の $d$ 次元球面だけを因子にした背景は $R_{ij}=(d-1)G_{ij}/L^2$ を持ち、その因子単独では条件を満たさない。近似には $\alpha'/L^2\ll1$ を要する。

### 一ループで現れる曲率の検算

背景写像 $\bar X$ に沿う接ベクトル $\xi$ を使い、$X=\exp_{\bar X}\xi$ とする。$B=0$ の計量項の二次変分は、Jacobi 場の二次形式で
$\int[(D_a\xi)^2-R_{\mu\rho\nu\sigma}\xi^\mu\xi^\nu\partial_a\bar X^\rho\partial^a\bar X^\sigma]/(4\pi\alpha')$
である。ここで球面に正の Ricci 曲率を割り当てる。実際、測地線に沿う変分と世界面方向の共変微分を入れ替えると $[D_s,D_a]\xi=R(\partial_s X,\partial_a X)\xi$ が生じ、$\partial_s^2 X=0$ を使って二次式に戻すと曲率項が得られる。高運動量の局所縮約は

$$
2\pi\alpha'G^{\mu\nu}\int_{\mu<|k|<\Lambda}\frac{\dd^2k}{(2\pi)^2k^2}
=\alpha'G^{\mu\nu}\log(\Lambda/\mu).
$$

従って曲率項は $-\alpha'R_{\rho\sigma}\log(\Lambda/\mu)$ だけ計量を変える。裸の計量を固定し、低いカットオフ $\mu$ における計量を微分すれば $\mu\partial_\mu G_{\rho\sigma}=\alpha'R_{\rho\sigma}+O(\alpha'^2)$ である。背景の長さ $L$ が $\ell_s$ に近いと高次の項を捨てることはできない。

:::example[球面座標での流れ]
半径 $L$ の $S^2$ で $G=L^2\hat G$、$\hat G$ が単位球面計量なら $R_{ij}=\hat G_{ij}$ である。従って $\mu\partial_\mu L^2=\alpha'$、$L^2(\mu)=L^2(\mu_0)+\alpha'\log(\mu/\mu_0)$ となる。$\log(\mu/\mu_0)=-2$ なら $L^2$ の変化は $-2\alpha'$ である。この結果が使えるのは全区間で $L^2\gg\alpha'$ のときである。
:::

:::example[二形式の線形化]
平坦背景で $\Phi$ 一定、$B_{12}=b\cos(kX^3)$、その他ゼロとする。$\partial^\rho H_{\rho12}=-k^2 B_{12}$ なので :eqref[eq-string-weyl-beta] の線形項は $\beta^B_{12}=(\alpha' k^2/2)B_{12}$ である。一定の $B$ ではこの項は消え、コンパクト時空での大域的な位相はこの局所検査とは別に残る。
:::

### 捩率付き接続で計量と二形式をまとめる

$H$ が全反対称であるとして、$\Gamma^{+\rho}_{\mu\nu}=\Gamma^\rho_{\mu\nu}+\tfrac12 H^\rho{}_{\mu\nu}$ と定める。ここでは $\nabla_\mu v^\rho=\partial_\mu v^\rho+\Gamma^\rho_{\mu\nu}v^\nu$、$R_{\sigma\nu}=R^\rho{}_{\sigma\rho\nu}$ とする。$K^\rho{}_{\mu\nu}=\tfrac12H^\rho{}_{\mu\nu}$ を接続の差とすれば、曲率の定義に $\Gamma+K$ を代入し、

$$
R^+_{\sigma\nu}=R_{\sigma\nu}+\nabla_\rho K^\rho{}_{\nu\sigma}
-\nabla_\nu K^\rho{}_{\rho\sigma}
+K^\rho{}_{\rho\lambda}K^\lambda{}_{\nu\sigma}
-K^\rho{}_{\nu\lambda}K^\lambda{}_{\rho\sigma}
$$

を得る。縮約 $K^\rho{}_{\rho\sigma}$ はゼロである。最後の積は添字の巡回置換により $-\tfrac14H_{\sigma\rho\lambda}H_\nu{}^{\rho\lambda}$、最初の微分は $-\tfrac12\nabla^\rho H_{\rho\sigma\nu}$ である。従って $R^+_{\mu\nu}$ の対称部と反対称部は、ディラトンが一定の :eqref[eq-string-weyl-beta] の括弧をそれぞれ与える。捩率を定義しただけで量子補正が求まるわけではないが、両方の係数の幾何学的な組合せをこの式で照合できる。

:::example[平坦な三つの空間座標に一定の $H$ を入れる]
直交座標で $H_{ijk}=h\epsilon_{ijk}$、$i,j,k\in\{1,2,3\}$ とする。$H_{i\rho\sigma}H_j{}^{\rho\sigma}=2h^2\delta_{ij}$ であるから、平坦で $\Phi$ が一定なら $\bar\beta^G_{ij}/\alpha'=-h^2\delta_{ij}/2$ になる。$\nabla^\rho H_{\rho ij}=0$ で二形式側の係数が消えても、計量側はゼロにならない。
:::

:::example[三次元の球面での流束と曲率]
半径 $L$ の三球を正規直交基底で見て $H=h\,\mathrm{vol}_3$ とする。$R_{ij}=2\delta_{ij}/L^2$、体積形式は平行なので $\nabla H=0$ である。計量と二形式の条件は $2/L^2-h^2/2=0$、すなわち $|h|=2/L$ でそろう。ディラトンのベータ関数と H の量子化はさらに課す。球面の曲率が弦長の尺度に近い場合、高次の補正をこの一ループの式だけで判定できない。
:::

### 二形式を含む二次作用を共変に展開する

背景写像 $\bar X$ のまわりを $X_s=\exp_{\bar X}(s\xi)$ とし、$e_a^\mu=\partial_a\bar X^\mu$ を定義する。すると $D_s\xi=0$、$D_s\partial_a X_s=D_a\xi$ である。局所的な高運動量の計算では標的の正規直交枠を使い、Lorentz の時間成分は Wick 回転で縮約を定義して最後にテンソル式へ戻す。運動項の二次部分は $|D\xi|^2-\langle R(\xi,e_a)e^a,\xi\rangle$ である。二形式の第一変分は、引き戻しに対する Cartan の公式 $\mathcal L_\xi B=\iota_\xi H+\dd(\iota_\xi B)$ から
$$
\frac{i}{4\pi\alpha'}\int\sqrt\gamma\,\varepsilon^{ab}
 H_{\lambda\mu\nu}\xi^\lambda e_a^\mu e_b^\nu
$$
となる。全微分は閉世界面で消える。これを $s$ でさらに微分し、Taylor の $1/2$ を掛けると、
$$
\frac{i}{4\pi\alpha'}\int\sqrt\gamma\,\varepsilon^{ab}
\left[\frac12\nabla_\kappa H_{\lambda\mu\nu}
\xi^\kappa\xi^\lambda e_a^\mu e_b^\nu+
H_{\lambda\mu\nu}\xi^\lambda D_a\xi^\mu e_b^\nu\right]
$$
を得る。$D_s\xi=0$ なので $\xi$ の差分は不要であり、二本の $e_a$ の変分は $\varepsilon$ と $H$ の二つの反対称性から同じ項を作る。

ここで $(K_a)_{\mu\lambda}=\frac i2\varepsilon_{ab}H_{\lambda\mu\nu}e^{b\nu}$ と置く。これを使うと $\xi D\xi$ の項は $2\langle D_a\xi,K^a\xi\rangle$ である。従って二次の運動項を $|(D+K)\xi|^2$ にまとめると、残るポテンシャルに $-\langle K_a\xi,K^a\xi\rangle$ が加わる。ここで二次形式は Euclid 経路積分の複素対称な形式であり、$i$ を共役しない。$K_a$ の成分を縮約すると、ポテンシャルの跡は
$$
\frac14H_{\mu\rho\sigma}H_\nu{}^{\rho\sigma} e_a^\mu e^{a\nu}
$$
である。実際 $i^2=-1$ と $\varepsilon_{ab}\varepsilon^{ac}=\delta_b{}^c$ を使い、前のマイナスを掛ければ正の $1/4$ が残る。従って $H^2$ の係数は、二形式に一次の混合を共変微分に組み替えた後の二次項から得られる。

### 局所的な縮約で繰り込み群の符号を決める

平らな世界面で二次作用を $-\mathcal D^2+U$ の形にすると、二次元の対数発散で二微分の背景に依存する項は $U$ の跡である。これは中心で接続をゼロにした正規枠の熱核の最初の係数から分かる。定数 $U$ の場合は $\exp[-s(k^2+U)]=e^{-sk^2}(1-sU+\cdots)$ なので、その係数は $-\Tr U$ である。接続の曲率は世界面の反対称な二添字を持つ。平坦な Laplace 型の運動演算子自体にはその二添字と縮約できる反対称テンソルがなく、曲率を一つだけ含む跡の項は出ない。次の次数には曲率二乗が現れる。曲がった世界面の熱核にはさらに $R^{(2)}/6$ があり、それは曲率結合に戻る。

従って $U$ を $\langle\xi^\mu\xi^\nu\rangle=\alpha'G^{\mu\nu}\log(\Lambda/\mu)$ で縮約すると、
$$
\begin{aligned}
\delta G_{\mu\nu}&=\alpha'\log(\Lambda/\mu)
\bigl[-R_{\mu\nu}+\tfrac14 H_{\mu\rho\sigma}H_\nu{}^{\rho\sigma}\bigr],\\
\delta B_{\mu\nu}&=\alpha'\log(\Lambda/\mu)\,
\tfrac12 \nabla^\rho H_{\rho\mu\nu}
\end{aligned}
$$
となる。二行目には元の作用の $i\varepsilon^{ab}$ を既に取り出した。$\Lambda$ を固定して $\mu$ で微分すれば、両行で対数の符号が反転する。従って $\beta^G_{\mu\nu}=\alpha'[R_{\mu\nu}-H_{\mu\rho\sigma}H_\nu{}^{\rho\sigma}/4]$、$\beta^B_{\mu\nu}=-\alpha'\nabla^\rho H_{\rho\mu\nu}/2$ である。ここでは局所的に一定に近い背景を使い、曲率、$H^2$、$\nabla H$ に対して次の微分次数を省いた。

次に Weyl 変分の $-\Delta_\gamma\Phi/(2\pi)$ を用いる。合成関数の Laplacian は $\Delta_\gamma\Phi=\nabla_\mu\nabla_\nu\Phi\,e_a^\mu e^{a\nu}+\nabla_\mu\Phi\,\mathcal D^a\partial_a X^\mu$ である。後者に $X$ の古典方程式を使うと、曲率の項以外の係数は
$$
-\frac1{4\pi\alpha'}\left[
2\alpha'\nabla_\mu\nabla_\nu\Phi\,e_a^\mu e^{a\nu}
+i\alpha'\varepsilon^{ab}\nabla^\rho\Phi\,H_{\rho\mu\nu}e_a^\mu e_b^\nu
\right]
$$
となる。Weyl 伸縮は物理的な繰り込みスケールの逆変化に対応するので、$\delta_\omega S_{\rm local}$ の係数は $-\bar\beta/(4\pi\alpha')$ の規約で読む。従って通常の $\beta$ にこれらの改善を足して :eqref[eq-string-weyl-beta] を得る。古典方程式を使って移した項は場の再定義に対応する。このため通常の繰り込み群の $\beta$ の零点と、改善を含めた Weyl 条件を同じ式として扱わない。

:::example[平坦な横波の Weyl 条件]
$h_{12}=\epsilon\cos(kX^3)$、$B=0$、$\Phi$ 一定とし、時空の第 3 方向を空間的とする。線形の $R_{12}$ は $-\Box h_{12}/2=k^2h_{12}/2$ である。従って $\bar\beta^G_{12}=\alpha' k^2h_{12}/2$ となる。$X^3$ を光的位相 $k\cdot X$ に置き換え、かつ $k^2=0$ と横条件を保てば、同じ式で $\bar\beta^G_{12}=0$ となる。例えば振幅が小さいだけでは、空間に固定された非零波数の波を共形な背景にはできない。
:::

:::example[三球の流束とディラトンの条件を分ける]
三球で $h=2/L$ とすると $\bar\beta^G=\bar\beta^B=0$ である。しかし $H^2=6h^2=24/L^2$ なので、一定ディラトンで $D=26$ なら $\bar\beta^\Phi=-\alpha'/L^2$ が残る。別の平坦で空間的な方向に $V^2=1/L^2$ の線形ディラトンを入れると $\alpha'V^2$ で相殺できる。この方向と三球は直積であり、$H$ にその方向の成分はないため、改善 $V^\rho H_{\rho\mu\nu}$ はゼロである。これは二微分の式の解であり、有限のレベルにおける正確な中心電荷の一致は別に確認する。
:::

### 標的の座標変換と有効方程式の整合

計量、二形式、ディラトンの Weyl 条件は独立なテンソルの形を持つが、任意の場についてまったく無関係に選べるわけではない。一般座標変換に不変な作用 $S$ の変分を
$$
\delta S=\int\sqrt{-G}\,
(E^{\mu\nu}\delta G_{\mu\nu}+J^{\mu\nu}\delta B_{\mu\nu}+K_\Phi\delta\Phi)
$$
と書く。$J^{\mu\nu}$ は反対称であり、二形式のゲージ変換だけを入れると $\nabla_\mu J^{\mu\nu}=0$ を得る。座標変換では $\delta G_{\mu\nu}=2\nabla_{(\mu}\xi_{\nu)}$、$\delta B=\iota_\xi H+\dd(\iota_\xi B)$、$\delta\Phi=\xi^\rho\partial_\rho\Phi$ である。ゲージ変換の項を除き、計量の微分を部分積分すると
$$
-2\nabla_\mu E^\mu{}_\rho+J^{\mu\nu}H_{\rho\mu\nu}
+K_\Phi\,\partial_\rho\Phi=0
$$
となる。これは場が方程式を満たす前に成り立つ恒等式である。従って計量と二形式の変分式をゼロとした場合には、$K_\Phi\partial_\rho\Phi=0$ である。$\Phi$ が一定ならこの恒等式だけでディラトン式を要求できない。このことは、平坦かつ一定 $\Phi$ の背景で $D=26$ の条件を計量から導けないこととも一致する。

### Weyl スカラーと残りの方程式の関係

$E^G_{\mu\nu}=\bar\beta^G_{\mu\nu}/\alpha'$、$E^B_{\mu\nu}=\bar\beta^B_{\mu\nu}/\alpha'$ と書く。また $f=e^{-2\Phi}$、$\mathcal C=2(D-26)/(3\alpha')$ とし、
$$
C=R+4\Box\Phi-4(\nabla\Phi)^2-\frac1{12}H^2-\mathcal C
$$
を定める。最低微分次数の時空作用を変分すると、変数を $G_{\mu\nu}$ とした計量の係数は $-f(E^{G,\mu\nu}-G^{\mu\nu}C/2)$ であり、$B$ と $\Phi$ の係数はそれぞれ $-fE^{B,\mu\nu}$、$-2fC$ になる。前に求めた座標変換の恒等式へこれらを入れると
$$
\partial_\rho C
=2 f^{-1}\nabla_\mu(fE^{G,\mu}{}_\rho)
-H_{\rho\mu\nu}E^{B,\mu\nu}
$$
である。$\partial_\rho f=-2f\partial_\rho\Phi$ を代入すると、$C\partial_\rho\Phi$ の二項は相殺する。従って $E^G=E^B=0$ なら各連結成分で $C$ は一定である。しかし一定値がゼロであることはこの式からは分からない。ゴーストまで含む世界面の中心電荷の条件が、その定数をゼロに定める。

$\bar\beta^\Phi$ としてよく用いる量は $-\alpha'(C-G^{\mu\nu}E^G_{\mu\nu})/4$ である。これにより $\frac{D-26}{6}+\alpha'[(\nabla\Phi)^2-\Box\Phi/2-H^2/24]$ を得る。$C$ と $\bar\beta^\Phi$ のどちらを使っても、三つの Weyl 条件を解く集合は同じである。二つのテンソルの式だけから局所的に任意の $C$ を消去することはできない。例えば平坦かつ定数の背景なら $E^G=E^B=0$ であるが、$C=-\mathcal C$ が残る。

:::example[平坦な不足次元での検算]
$D=24$ かつ $H=0$ とすると、$\Phi$ が一定の場合は $C=4/(3\alpha')$ である。空間方向に勾配を入れて $\alpha'V^2=1/3$ にすると、$-4V^2$ がこれを打ち消す。勾配が一定なので $E^G$ と $E^B$ の零は保たれる。$V^2$ の符号は時空の Lorentz 計量で判定している。
:::

### 座標の取り方と場の再定義を区別する

二微分の結果を越えると、有効作用の局所項の表示には場の再定義の自由度がある。例えば $G_{\mu\nu}$ を $G_{\mu\nu}+a\alpha'R_{\mu\nu}$ に置き換える。二微分の作用 $S_0$ にこれを入れると
$$
S_0[G+a\alpha'R]=S_0[G]+a\alpha'\int
\frac{\delta S_0}{\delta G_{\mu\nu}}R_{\mu\nu}+O(\alpha'^2)
$$
である。元の変分は二微分なので追加項は四微分であり、既に最低次数の方程式を満たす背景ではこの一次の追加項は消える。従ってあるスキームの四微分の項と別のスキームの項を比較するとき、Ricci テンソルに比例する全ての係数を物理的な違いと読むことはできない。一方オンシェルの散乱振幅や、再定義を含めて解を対応づけたときの物理量は保たれる。世界面の低エネルギー条件 $\alpha'/L^2\ll1$ が必要なのは、このような再定義も微分の次数ごとに扱っているためである。

### 光的な曲率で局所条件を検査する

例として Brinkmann 座標で
$$
\dd s^2=-2\dd u\dd v+K(u,x)\dd u^2+\sum_i(\dd x^i)^2
$$
を取る。$G^{uu}=0$、$G^{uv}=-1$、$G^{vv}=-K$ である。接続の公式を入れると、曲率を作る独立な成分は $R_{uiuj}=-\frac12\partial_i\partial_j K$ で、$R_{uu}=-\frac12\sum_i\partial_i^2K$ となる。実際 $\Gamma^i_{uu}=-\partial_i K/2$ であり、$R^i{}_{uju}$ に空間微分を取ればこの項が出る。$v$ 方向の微分も下付き $v$ を持つ曲率もない。$B=0$、$\Phi$ 一定、$D=26$ では $K$ が横で調和なら二微分の Weyl 条件を満たす。

例えば $K=A(u)[(x^1)^2-(x^2)^2]$ は横の Laplacian がゼロであるが、$R_{u1u1}=-A$、$R_{u2u2}=A$ は一般に非零である。曲率スカラーや曲率ノルムを作るためには下付き $u$ を縮約しなければならないが、$G^{uu}=0$ で、$G^{uv}$ と組になる曲率の $v$ 成分もない。このため全ての曲率スカラーが消え得る。二微分近似の検査において、曲率の不変スカラーがゼロであることだけから、任意の観測者の潮汐力までゼロとは言えない。この平面波での高次の局所補正の検査にも同じ添字の数え上げが使える。二つ以上の曲率因子は合計で四つ以上の下付き $u$ を持つから、二階テンソルにしても残りの $u$ を縮約できない。平面波で横微分がない曲率の共変微分は、同じ数え上げにさらに $u$ を加える。これは共変な局所微分展開の議論であり、時空の境界条件や種数補正の判定は別に必要である。

## 弦フレームと Einstein フレーム

臨界次元の NS–NS 作用は最低微分次数で
$$
S_{\rm eff}=\frac1{2\kappa_0^2}\int\dd^D x\,\sqrt{-G}\,e^{-2\Phi}
\left(R[G]+4|\nabla\Phi|^2-\frac1{12}H_{\mu\nu\rho}H^{\mu\nu\rho}\right).
\label{eq-string-nsns-action}
$$
$B$ を変分すると $\delta H_{\mu\nu\rho}=3\nabla_{[\mu}\delta B_{\nu\rho]}$ より $\nabla_\mu(e^{-2\Phi}H^{\mu\nu\rho})=0$ となる。$\Phi$ を変分して部分積分すると $R+4\Box\Phi-4|\nabla\Phi|^2-H^2/12=0$ である。計量変分で $e^{-2\Phi}R$ の微分項 $(G_{\mu\nu}\Box-\nabla_\mu\nabla_\nu)e^{-2\Phi}$ を含め、上のスカラー式でトレースに比例する部分を落とすと :eqref[eq-string-weyl-beta] の計量方程式になる。

### ディラトンと中心電荷の不足

世界面の曲率に掛かる Weyl 応答を $\bar\beta^\Phi$ と書く。平坦な $D$ 個のボソンとゴーストでは $c_{\rm tot}=D-26$ なので、定数項を $(D-26)/6$ と規格化する。非臨界次元では :eqref[eq-string-nsns-action] の括弧に $-2(D-26)/(3\alpha')$ が加わる。その係数は後の節の線形ディラトンで照合できる。実際、線形背景の正確な中心電荷を 26 にする勾配は $V^2=(26-D)/(6\alpha')$ となる。スカラーの方程式にこれを入れると $-4V^2-2(D-26)/(3\alpha')=0$ になる。

この定数を含むスカラー方程式から計量方程式のトレース $R-H^2/4+2\Box\Phi=0$ を引き、$-4$ で割る。すると

$$
\bar\beta^\Phi=\frac{D-26}{6}
+\alpha'\left(|\nabla\Phi|^2-\frac12\Box\Phi-\frac1{24}H^2\right)=0
$$

となる。ここでは計量ベータ関数に比例する項を移し替えた規約を用いた。例えば $H=0$ で線形な背景では、空間的な勾配の二乗が $D<26$ の負の定数を相殺する。$R_{\mu\nu}=0$ だけでは定数項への条件が得られないので、三つの Weyl 条件を合わせて使う。

### Einstein フレームへの変換

$D>2$ で $G=e^{4(\Phi-\Phi_0)/(D-2)}G_E$ とする。これを :term[Einstein フレーム]{en="Einstein frame" yomi="あいんしゅたいんふれーむ"} という。$G=e^{2\omega}G_E$ の曲率は $e^{-2\omega}[R_E-2(D-1)\Box_E\omega-(D-1)(D-2)|\nabla\omega|^2]$ である。従って体積と曲率と $e^{-2\Phi}$ の積の前置因子は定数となる。全微分を落として元の $4|\nabla\Phi|^2$ を合わせると係数は $-4/(D-2)$ となる。$H^2$ は逆計量三個を含むため残る因子は $e^{-8(\Phi-\Phi_0)/(D-2)}$ である。例えば $D=10$ で運動項は $-|\nabla\Phi|^2/2$、三形式の因子は $e^{-(\Phi-\Phi_0)}$ になる。

### 二形式の作用とソースの規約

$g_{\rm s}=e^{4\Phi/(D-2)}g_{\rm E}$ とし、定数ディラトンの差は座標のスケールに含める。$H$ は三形式なので三つの逆計量が入り、$\sqrt{-g_{\rm s}}e^{-2\Phi}H_{\rm s}^2=\sqrt{-g_{\rm E}}e^{-8\Phi/(D-2)}H_{\rm E}^2$ となる。Einstein フレームでは $B$ に関する作用の変分は

$$
\delta S_H
=-\frac{1}{2\kappa^2}\int e^{-8\Phi/(D-2)}\dd(\delta B)\wedge *_{\rm E}H.
$$

境界で $\delta B=0$ とすれば、部分積分により $\dd(e^{-8\Phi/(D-2)}*_{\rm E}H)=0$ を得る。$H=\dd B$ からの $\dd H=0$ は変分前から成り立つ。電気的な弦源があると前者の右辺が、磁気的な源があると後者が変わる。

:::example[十次元における規格化]
$D=10$ では場の運動項は $-(1/2)(\partial\Phi)^2$、$H$ の因子は $e^{-\Phi}$ である。$\Phi=\log2$ で同じ Einstein 計量と同じ三形式成分を比較すれば、三形式の運動項は $\Phi=0$ の半分である。作用全体の係数 $1/(2\kappa^2)$ はここでは共通に取った。
:::

:::example[二種類の測地線]
定数でない $\Phi$ のもとで、点粒子の弦フレーム固有長は $\dd s_{\rm s}=e^{2\Phi/(D-2)}\dd s_{\rm E}$ である。Einstein フレームに移すと固定質量の世界線は位置依存の質量 $m e^{2\Phi/(D-2)}$ を持つ。従って変分ではこの因子の勾配も力に入る。無質量の光的世界線の方向は共形変換で保たれる。
:::

### 作用の三つの変分をそろえる

$f=e^{-2\Phi}$ と書く。$G^{\mu\nu}$ による変分では $\delta\sqrt{-G}=-\frac12\sqrt{-G}G_{\mu\nu}\delta G^{\mu\nu}$ であり、$H$ の三つの逆計量から $3H_{\mu\rho\sigma}H_\nu{}^{\rho\sigma}\delta G^{\mu\nu}$ が出る。曲率を部分積分すると $f(R_{\mu\nu}-\frac12G_{\mu\nu}R)+(G_{\mu\nu}\Box-\nabla_\mu\nabla_\nu)f$ である。$f^{-1}\nabla_\mu\nabla_\nu f=-2\nabla_\mu\nabla_\nu\Phi+4\nabla_\mu\Phi\nabla_\nu\Phi$ だから、ディラトンの運動項の $4\nabla_\mu\Phi\nabla_\nu\Phi$ と後者が相殺する。残る計量変分は
$$
R_{\mu\nu}-\frac14 H_{\mu\rho\sigma}H_\nu{}^{\rho\sigma}
+2\nabla_\mu\nabla_\nu\Phi
-\frac12 G_{\mu\nu}\left(R+4\Box\Phi-4(\nabla\Phi)^2-\frac1{12}H^2\right)
$$
である。括弧は $\Phi$ の変分を $-2f$ で割った式である。従って両方の変分をゼロにすれば、前の三項だけの Weyl 条件に書き換えられる。$B$ の方程式 $\nabla_\mu(fH^{\mu\nu\rho})=0$ は、$f^{-1}\partial_\mu f=-2\partial_\mu\Phi$ から $\bar\beta^B=0$ と等価である。

中心電荷不足の定数を作用の括弧に $-\mathcal C$ として加えた場合、計量変分の最後の括弧にも $-\mathcal C$ が加わる。従って $\Phi$ の式を使ってトレースに比例する項を取り除いた後の計量条件は変わらない。$H=0$、$G=\eta$、$\Phi=V\cdot X$ ならスカラー方程式は $-4V^2-\mathcal C=0$ である。次節の自由場の厳密な OPE は $D+6\alpha'V^2=26$ を与えるので、$\mathcal C=2(D-26)/(3\alpha')$ である。ボソンの中心電荷の規約であるため、内部 CFT や世界面フェルミオンがある場合は同じ $D-26$ を無条件に代入しない。

### Einstein フレームの作用を全項で書く

$\varphi=\Phi-\Phi_0$ とし、$G=e^{2\omega}G_E$、$\omega=2\varphi/(D-2)$ とする。接続の差は
$$
\Gamma^\rho_{\mu\nu}[G]-\Gamma^\rho_{\mu\nu}[G_E]
=\delta^\rho_\mu\partial_\nu\omega+\delta^\rho_\nu\partial_\mu\omega
-G_{E,\mu\nu}G_E^{\rho\sigma}\partial_\sigma\omega
$$
である。接続の差を Ricci 曲率に入れ、逆計量 $e^{-2\omega}G_E^{-1}$ で縮約すれば、前に用いた曲率の Weyl 変換を得る。体積の $e^{D\omega}$ と曲率の $e^{-2\omega}$ および $e^{-2\Phi}$ を掛けると $e^{-2\Phi_0}$ が残る。この一定の値は $\kappa^2=\kappa_0^2 e^{2\Phi_0}$ に吸収する。従って臨界次元での作用は
$$
S_E=\frac1{2\kappa^2}\int\dd^D x\sqrt{-G_E}
\left[
R_E-\frac4{D-2}(\nabla_E\varphi)^2
-\frac1{12}e^{-8\varphi/(D-2)}H_E^2
\right],
\label{eq-string-einstein-nsns}
$$
である。ここでは $\Box_E\omega$ を境界項として落とせる条件を仮定した。非臨界の定数 $-\mathcal C$ には曲率の $e^{-2\omega}$ がないので、Einstein フレームでは $-\mathcal C e^{4\varphi/(D-2)}$ となる。これが定数のポテンシャルではなく指数関数になるのは、体積とディラトンの共通の因子が曲率とは異なるためである。$D=2$ では曲率の前の指数をこの変換で消せないため、$D>2$ が必要である。

正準な実スカラー $\varphi_c$ の運動項を全作用の中で $-(\partial\varphi_c)^2/2$ とするなら、
$$
\varphi_c=\frac{2}{\kappa\sqrt{D-2}}\,\varphi
$$
である。例えば $D=10$ では $\varphi_c=\varphi/(\sqrt2\kappa)$ になる。括弧の中だけで運動項を $-1/2$ とする規約と、全作用でスカラーを正準化する規約を区別する。張力に結び付くディラトンの線形源を計算する場合には、場の二点関数の規格化も同時に合わせる。

:::example[弦の線形ディラトン結合を Einstein フレームで読む]
世界体積次元が 2 の基本弦の面積要素は $G=e^{2\omega}G_E$ で $e^{2\omega}$ 倍になる。従って固定した弦フレーム張力 $1/(2\pi\alpha')$ は、Einstein フレームの作用では $T_{\rm F1}e^{4\varphi/(D-2)}$ を伴う。$D=10$ で $\varphi=\log4$ なら係数は $4^{1/2}=2$ であり、同じ Einstein フレームの面積に対して作用の絶対値は二倍になる。
:::

### 基本弦と磁気源を変分から区別する

三形式の運動項を微分形式で $S_H=-\frac1{4\kappa^2}\int f_E H\wedge *H$ と書く。ここで $f_E=e^{-8\varphi/(D-2)}$ であり、標的時空には境界がないか、変分を境界でゼロとする。向き付きの基本弦 $\Sigma_2$ の結合を $+q\int_{\Sigma_2} B$ と取る。Poincaré 双対の分布形式 $J_{D-2}$ を $\int\delta B\wedge J_{D-2}=\int_{\Sigma_2}\delta B$ によって定義する。すると $\dd(\delta B\wedge f_E*H)=\dd\delta B\wedge f_E*H+\delta B\wedge \dd(f_E*H)$ だから、
$$
\dd(f_E*H)=-2\kappa^2 q J_{D-2}
$$
となる。この符号は源の向きと $+q\int B$ の選び方に対応する。$\Sigma_2$ が閉じていれば $\dd J_{D-2}=0$ である。これは $B\mapsto B+\dd\Lambda$ に対する源のゲージ不変性、あるいは方程式の $\dd^2=0$ からも分かる。一方、磁気源は $\dd H$ に加わる。電気源の方程式を書き換えても、それだけで Bianchi 恒等式に右辺は出ない。

例えば十次元の直線状の基本弦の横空間は八次元である。十分遠くで平坦、かつ $\Phi$ がほぼ一定なら電位の横 Laplace 方程式の解は $r^{-6}$ に比例する。それを微分した場は $r^{-7}$ であり、半径 $r$ の七球の面積の $r^7$ と相殺する。従ってフラックスの積分は半径に依存しない。この近似は遠方の線形化であり、源に近い領域の曲率や結合をこの冪だけから判定できない。

### 円周の半径と下の次元のディラトン

空間的な周期 $2\pi$ の座標 $\theta$ に対し、簡単のため混合成分と二形式をゼロとし、
$$
\dd s_{D}^2=\dd s_{D-1}^2+\ell_s^2 e^{2\sigma(x)}\dd\theta^2
$$
と書く。$\sigma$ や $\Phi$ は $\theta$ に依存しないとする。体積の因子は $\ell_s e^\sigma$ であるため、$\phi=\Phi-\sigma/2$ と置くと $\sqrt{-G_D}e^{-2\Phi}=\ell_s\sqrt{-G_{D-1}}e^{-2\phi}$ である。接続の追加成分は $\Gamma^\theta_{\mu\theta}=\partial_\mu\sigma$、$\Gamma^\mu_{\theta\theta}=-G_{\theta\theta}\partial^\mu\sigma$ である。これを Ricci テンソルに入れると
$$
R_D=R_{D-1}-2\Box\sigma-2(\nabla\sigma)^2
$$
となる。例えば $R_{\mu\nu}^{(D)}$ の差は $-\nabla_\mu\nabla_\nu\sigma-\partial_\mu\sigma\partial_\nu\sigma$、縮約した円周成分 $G^{\theta\theta}R_{\theta\theta}$ は $-\Box\sigma-(\nabla\sigma)^2$ であり、両方を足してこの式になる。

元の $4(\nabla\Phi)^2=4(\nabla\phi)^2+4\nabla\phi\cdot\nabla\sigma+(\nabla\sigma)^2$ を入れる。$\int e^{-2\phi}(-2\Box\sigma)=-4\int e^{-2\phi}\nabla\phi\cdot\nabla\sigma$ なので混合項が落ちる。従って下の次元の弦フレームでは括弧が
$$
R_{D-1}+4(\nabla\phi)^2-(\nabla\sigma)^2
$$
となる。周期方向の積分から全体には $2\pi\ell_s$ が掛かる。円周の T 双対性は $\sigma'=-\sigma$、$\Phi'=\Phi-\sigma$ を与えるから、$\phi'=\phi$ である。この部分の作用は $\sigma$ の微分の二乗にしか依存しないため、縮約後にも同じ形である。計量に混合成分があれば KK ベクトルが入り、二形式の混合成分との交換も合わせて扱う。

:::example[半径の期待値が結合に入る]
初期値 $\sigma=\log2$、$\Phi=\log g_s$ とする。双対後は $\sigma'=-\log2$、$\Phi'=\log(g_s/2)$ である。両方とも $\phi=\log g_s-\frac12\log2$ を与える。従って十次元の $g_s$ は変わるが、円周縮約後の重力結合の前置因子を決める値は同じである。
:::

## 線形ディラトン背景

平坦な $G$ と $H=0$ で $\Phi=V_\mu X^\mu$ とする。曲率結合を計量で変分すると、正則応力テンソルは $T=-:\partial X\cdot\partial X:/\alpha'+V\cdot\partial^2X$ である。:eqref[eq-string-x-ope] を $z,w$ で各二回微分すると $\partial^2 X^\mu(z)\partial^2X^\nu(w)\sim3\alpha'\eta^{\mu\nu}/(z-w)^4$ になる。混合項には恒等演算子の四重極がないので $c=D+6\alpha' V^2$ である。例えば $D=20$ なら空間的 $V$ に対して $\alpha'V^2=1$ が条件となる。$D>26$ なら形式上は時間的な $V$ が必要である。局所の Weyl 条件を満たしても $g_s(X)=e^{V\cdot X}$ が大きい領域では種数展開を制御できない。

### 線形背景での重みと零モード

$\Phi(X)=V_\mu X^\mu$ のとき、正則応力テンソルの改善は $V\cdot\partial^2X$ である。これと $V_k=:e^{ikX}:$ の縮約は $\partial^2X^\mu(z)V_k(w)\sim i\alpha'k^\mu V_k(w)/[2(z-w)^2]$ なので、重みは $\alpha'(k^2+2iV\cdot k)/4$ となる。改善同士の縮約は $\partial^2X^\mu(z)\partial^2X^\nu(w)\sim3\alpha'\eta^{\mu\nu}/(z-w)^4$ であり、中心項は $6\alpha' V^2$ だけずれる。さらに曲率結合の零モードは $e^{-\chi V\cdot x}$ なので、相関の零モード積分を解析接続で扱うと $\sum_i k_i+i\chi V=0$ という背景電荷の条件になる。平坦な世界面で局所縮約が従来と同じでも、曲率の積分によってこの条件は変わる。

:::example[不足した一つの中心電荷]
$D=25$ かつ空間的な $V$ を取ると、$c=26$ には $|V|=1/\sqrt{6\alpha'}$ が必要である。$X$ を勾配方向に $\sqrt{6\alpha'}\log2$ だけ進めると、局所 $g_s(X)=e^{\Phi(X)}$ は 2 倍になる。中心電荷が一定でもループ近似の精度は位置に依存する。
:::

:::example[重みを実にする運動量表示]
$k=p-iV$ とおけば $k^2+2iV\cdot k=p^2+V^2$ である。例えば $\alpha'V^2=1/6$、$p^2=0$ なら $h=1/24$ である。普通の実運動量 $k$ を代入した複素の重みと、背景に合わせて平面波を再定義したものを区別する。
:::

### 改善項は局所の Weyl 変分から出る

応力テンソルの規約を $T_{ab}=-4\pi\,(\sqrt\gamma)^{-1}\delta S/\delta\gamma^{ab}$ として、平らな共形座標で $T_{zz}=-:\partial X\cdot\partial X:/\alpha'$ となるように複素座標を用いる。$\int\sqrt\gamma R^{(2)}\Phi$ の計量変分では二次元の Einstein テンソルがゼロであり、$(\gamma_{ab}\Delta_\gamma-\nabla_a\nabla_b)\Phi$ だけが残る。従って平らな世界面の $zz$ 成分では $+\partial^2\Phi=V\cdot\partial^2X$ が加わる。作用上で $R^{(2)}=0$ を先に代入して変分を消すと、この項は得られない。

OPE の残りの特異項も比較できる。一組の自由座標において $T_0(z)\partial^2X(w)$ は、$T_0(z)\partial X(w)$ を $w$ で微分して
$$
\frac{2\partial X(w)}{(z-w)^3}
+\frac{2\partial^2 X(w)}{(z-w)^2}+\frac{\partial^3 X(w)}{z-w}
$$
となる。逆向きの $\partial^2X(z)T_0(w)$ は三重極が $-2\partial X(w)/(z-w)^3$ であり、これらは相殺する。改善同士では恒等演算子の四重極だけがあり、係数は $3\alpha'V^2$ である。従って $T(z)T(w)$ の二重極と単極には改善を含む $2T(w)$、$\partial T(w)$ が現れる。中心項だけを変えた式が、応力テンソル全体の OPE と整合することをこの計算で確認する。

### 背景電荷と二点関数の共役運動量

球面では $\chi=2$ なので、零モード積分の重みは $\exp[-2V\cdot x+i(k+k')\cdot x]$ である。非コンパクトな線形ディラトンの方向では、実数の運動量を入れたこの積分は通常の Fourier 積分として収束しない。解析接続した背景電荷の規則を使うか、波動関数の内積に測度 $e^{-2\Phi}$ を入れる。二点関数の零モードの条件は $k'=-k-2iV$ であり、平坦背景の $k'=-k$ からずれる。$k=p-iV$ と置くと $e^{ikX}=e^{V\cdot X}e^{ipX}$ である。二つの波動から来る $e^{2V\cdot X}$ が内積の測度で消えるので、対応する $p$ の内積は通常の連続な平面波のものになる。これはデルタ関数規格化の波であって、個々の実 $p$ の波を二乗可積分と主張しているわけではない。

さらに自由なレベル $N$ の正則頂点に対して
$$
L_0-1=\frac{\alpha'}4(p^2+V^2)+N-1
$$
である。非コンパクト方向に時間を残した空間的 $V$ を取れば、$p$ の分散関係の質量二乗は $4(N-1)/\alpha'+V^2$ となる。例えば $N=1$ の自由な波では $V^2$ が質量二乗のずれとして入る。背景は並進に不変な結合定数を持たないため、これは線形化した波の分散関係である。全相互作用の S 行列に平坦背景と同じ保存則を課せるという意味ではない。

:::example[二つの運動量表示で重みを照合する]
球面の二点で $k=p-iV$、$k'=-p-iV$ を取る。和は $-2iV$ で背景電荷の条件を満たす。各重みには同じ $p^2+V^2$ が現れ、零モードでも二点が共役になる。種数 1 なら $\chi=0$ なので、同じ相関に使う背景電荷条件が変わる。局所 OPE の改善は種数に依存せず、零モードの積分だけが種数を区別する。
:::

### 空間的、光的、時間的な勾配

Lorentz 内積 $V^2=-(V^0)^2+\sum_i(V^i)^2$ である。$D<26$ なら $V^2>0$ の空間的勾配が中心電荷を補う。$D=26$ では $V^2=0$ が必要だが、$V=0$ に限らず光的な $V$ もこの条件を満たす。例えば光円錐座標で $\Phi=vX^+$ とすれば、$G^{++}=0$ なので $V^2=0$ である。それでも $g_s=e^{vX^+}$ は一定でない。$D>26$ では $V^2<0$ が必要で、時間的な勾配になる。この形式的な中心電荷の相殺から、時間依存背景での真空や散乱状態の選び方は決まらない。

例えば $D=23$ の空間的背景では $\alpha'V^2=1/2$ である。勾配に沿って距離 $\sqrt{2\alpha'}$ だけ進むと $\Phi$ は 1 増え、$g_s$ は $e$ 倍になる。$G$ の曲率はゼロであり、線形な自由場の OPE は勾配の全ての位置で同じである。それでも高種数の抑制は同じではない。弱結合側と強結合側を分けるには $\alpha'$ による微分展開と $g_s$ によるループ展開を独立に調べる。

### 三球の共形場理論で曲率近似を照合する

量子化した流束を持つ三球には、群 $SU(2)$ 上の Wess–Zumino–Witten 模型という共形な実現がある。この厳密な対応には、量子化した Wess–Zumino 項と保存カレントの正則性を示す必要があるため :cite[difrancesco-cft] を参照する。ここではそのレベル $k>0$ のカレント代数を使い、中心電荷の差を直接計算する。$[J_m^a,J_n^b]=i\epsilon^{abc}J_{m+n}^c+(k/2)m\delta^{ab}\delta_{m+n,0}$ とし、真空は $m\ge0$ で消える。カレントを重み 1 にする Sugawara の応力テンソルは、真空の状態として
$$
L_{-2}\ket0=\frac1{k+2}\sum_a J_{-1}^aJ_{-1}^a\ket0
$$
を与える。係数は $J_1^b$ を作用させて調べられる。交換子の中心項は $kJ_{-1}^b\ket0$ となり、残る $J_0$ をもう一つの $J_{-1}$ と交換すると $\epsilon$ の二つの積から $2J_{-1}^b\ket0$ が出る。従って $[J_1^b,L_{-2}]\ket0=J_{-1}^b\ket0$ となる係数は $1/(k+2)$ である。正規順序した二次のカレントから構成した Virasoro の $L_{-2}\ket0$ のノルムは、その式に $J_1^b$ を二回使って $3k/[2(k+2)]$ である。$[L_2,L_{-2}]\ket0=(c/2)\ket0$ と比較して $c=3k/(k+2)$ となる。

残りを 23 個の平坦方向とし、そのうち一つの空間方向に線形ディラトンを入れる。臨界条件は $23+3k/(k+2)+6\alpha'V^2=26$ なので $\alpha'V^2=1/(k+2)$ である。例えば $k=2$ なら $c=3/2$、$\alpha'V^2=1/4$ である。$k\gg1$ では $1/(k+2)=1/k+O(k^{-2})$ なので、三球の $L^2=k\alpha'$ を二微分の流束の式で求めたときの $\alpha'V^2\simeq\alpha'/L^2$ を再現する。有限の $k$ で半径をどう定義するかはスキームに依存するが、この中心電荷の比較は CFT の量である。球面の幾何学的な近似を使って計算した値が、常に有限 $k$ でも厳密になるわけではない。

## 弦の熱力学と Hagedorn 温度

振動子の状態数 $d_N$ を :eqref[eq-string-open-mass] の臨界開弦で求める。$q=e^{-u}$ として $\log\sum_Nd_N q^N=24\sum_{m\ge1}[m(e^{um}-1)]^{-1}$ である。$u$ を掛けた各項は $24/m^2$ 以下であり $u\to0^+$ でその値へ収束する。よって主項は $4\pi^2/u$ になる。逆 Laplace 変換の鞍点の指数は $uN+4\pi^2/u$ であり、$u=2\pi/\sqrt N$ で $4\pi\sqrt N$ となる。漸近展開の厳密な係数の評価は :cite[polchinski-1] を参照する。閉弦では $d_N^2$、$M\simeq2\sqrt{N/\alpha'}$ なので、いずれも指数部分は $\exp(4\pi\ell_s M)$ となる。

一弦の熱的トレースは、質量積分の指数に $-(\beta-4\pi\ell_s)M$ を含む。従って $\beta<\beta_H=4\pi\ell_s$ では高い質量の和が発散する。この指数の境界を :term[Hagedorn 温度]{en="Hagedorn temperature" yomi="はげどるんおんど"} $T_H=1/\beta_H$ と呼ぶ。境界そのものでの収束には冪の前因子と空間次元も必要である。熱的な Euclid 時間円で巻き付き一つのボソン閉弦は $M^2=(\beta/(2\pi\alpha'))^2-4/\alpha'$ を持ち、$\beta_H$ でゼロになる。ボソン弦では元から巻き付きゼロのタキオンがあるため、ここでは高質量の成長を独立に調べた。

### 一弦の熱トレースを質量密度に直す

非コンパクトな空間を $d$ 次元とする。$M\beta\gg1$ では $E(\vb p)=M+\vb p^2/(2M)+O(p^4/M^3)$ であり、積分を支配する運動量は $p\sim\sqrt{M/\beta}$ だから高次の指数の補正は $O((M\beta)^{-1})$ である。従って $\int\dd^dp/(2\pi)^d\,e^{-\beta E}\sim (M/(2\pi\beta))^{d/2}e^{-\beta M}$ となる。質量当たりの内部状態密度が $M^{-a}e^{\beta_H M}$ に比例するなら、$\beta=\beta_H$ での一弦の跡の高質量側は $\int^\infty M^{d/2-a}\dd M$ で判定する。収束条件は $a>1+d/2$ である。エネルギーの期待値にはさらに $M$ が一つ入るため、跡が有限でもエネルギーは発散しうる。

:::example[境界で跡とエネルギーを分ける]
$d=3$、模型の $a=3$ では高質量の積分は跡で $\int^\infty M^{-3/2}\dd M$、エネルギーで $\int^\infty M^{-1/2}\dd M$ となる。下限 $M_0>0$ から積分すると前者は $2/\sqrt{M_0}$、後者は上限の平方根に比例して発散する。
:::

:::example[熱巻き付きの質量の符号]
$\beta=6\pi\ell_s$ なら $w=1$ のボソン閉弦の最低質量は $\alpha'M^2=9-4=5$ である。$\beta=2\pi\ell_s$ なら $\alpha'M^2=1-4=-3$ となる。温度を上げると Euclid 時間円が短くなり、熱円を巻く弦の伸びのエネルギーが減る。
:::

### 状態密度の冪をモジュラー変換から求める

ボソン弦一側の振動子和を $F(u)=\prod_{n\ge1}(1-e^{-un})^{-24}$ とする。$t=u/(2\pi)$ と置くと $F(u)=e^{-u}\eta(it)^{-24}$ である。$\eta(it)=t^{-1/2}\eta(i/t)$ と $t\to0^+$ での $\eta(i/t)\sim \exp[-\pi/(12t)]$ を使うと
$$
F(u)\sim \left(\frac u{2\pi}\right)^{12}\exp\left(\frac{4\pi^2}{u}-u\right).
$$
従って $d_N$ の逆変換の主指数は $(N-1)u+4\pi^2/u$ である。$N\gg1$ では $N-1$ と $N$ の差は先頭の冪を変えない。$A=4\pi^2$ として $u_*=\sqrt{A/N}$ を取り、虚軸方向のずれを $iv$ と置く。主指数の二階微分は $2A/u_*^3=2N^{3/2}/\sqrt A$ なので、$v$ の Gaussian 積分は $N^{-3/4}$ に比例する。$u_*^{12}$ の $N^{-6}$ と合わせると
$$
d_N\sim C\,N^{-27/4}e^{4\pi\sqrt N}
$$
になる。定数 $C$ の値はここでは不要である。厳密な漸近の誤差評価では単位円の他の有理点の近傍も扱う。そこではモジュラー変換した尖点での指数の $1/c^2$ が $c\ge2$ となって主項より小さい。小弧の寄与を一様に足す厳密な方法には分割関数の円法が要るため、ここでは主指数と冪を与える鞍点まで計算した。状態密度の漸近評価については :cite[gsw-1] を参照する。

開弦では $N\simeq\alpha'M^2$、閉弦では整合した両側の準位が $N\simeq\alpha'M^2/4$ である。質量あたりの状態密度にする際はどちらも $\dd N/\dd M\propto M$ を乗じる。従って非コンパクトな平坦背景の一弦の内部状態密度は
$$
\rho_{\rm open}(M)\propto M^{-25/2}e^{4\pi\ell_s M},
\qquad
\rho_{\rm closed}(M)\propto M^{-26}e^{4\pi\ell_s M}
$$
となる。ここでの冪は内部の準位密度であり、連続な重心運動量はまだ積分していない。また、不安定なボソン弦の低準位のタキオンを形式的に除いて、高質量の裾だけを比較している。

:::example[同じ指数で異なる熱の跡]
臨界ボソン弦の非コンパクトな空間次元は $d=25$ である。空間充填開弦の高質量の熱トレースでは $\rho$ の冪と運動量積分の冪 $M^{25/2}$ が相殺し、$\beta=\beta_H$ で $\int^\infty \dd M$ が残る。一方閉弦の冪は $25/2-26=-27/2$ である。この場合、高質量側の跡は収束し、エネルギーの $k$ 次モーメントは $k-27/2<-1$、すなわち整数 $k\le12$ まで収束する。境界の温度の値が同じでも、自由エネルギーの特異性は同じにならない。
:::

### 一弦から多弦へ移るときの統計

相互作用を無視した気体で、一弦の固有状態のエネルギーを $E_j>0$ とする。ボソンは各固有状態を $n_j=0,1,\ldots$ 個占有できるので、分配関数の因子は $(1-e^{-\beta E_j})^{-1}$ である。フェルミオンは $n_j=0,1$ なので $1+e^{-\beta E_j}$ である。対数を展開して一弦の跡 $Z_{\rm B}(\beta)=\sum_{j\in{\rm B}}e^{-\beta E_j}$ と $Z_{\rm F}$ を使えば
$$
\log Z_{\rm gas}=\sum_{m\ge1}\frac1m
\left[Z_{\rm B}(m\beta)+(-1)^{m-1}Z_{\rm F}(m\beta)\right].
$$
この展開は各一弦の跡が収束し、ゼロエネルギーのボソンの占有に赤外の正則化がある領域で用いる。例えば連続の非コンパクトな空間で質量ゼロの場を扱うときは、体積当たりの量を取る。$m=1$ の項だけが $\beta=\beta_H$ で同じ高質量の境界に達する。$m\ge2$ では $m\beta>\beta_H$ なので、高質量の発散がこの境界を動かすことはない。これは重力や弦の相互作用の補正を含まない判定である。

:::example[フェルミオンの熱トレースの符号]
同じエネルギー $E$ を持つ一つのボソンと一つのフェルミオンを考え、$x=e^{-\beta E}$ とする。$\log Z=\log(1+x)-\log(1-x)=2(x+x^3/3+x^5/5+\cdots)$ である。偶数の $m$ の寄与は消えるが、$m=1$ では両状態を加える。従って零温度の超対称な符号付き真空図の相殺は、熱の状態数の指数を消さない。
:::

### Euclid 時間を巻く長管と熱特異性

ボソン閉弦の熱円の半径は $R_\beta=\beta/(2\pi)$ である。円周の質量式で $n=0,w=\pm1,N_L=N_R=0$ を取り、非コンパクトな空間の質量を $m_{\rm th}^2=(R_\beta/\alpha')^2-4/\alpha'$ とする。$\beta=\beta_H+\delta$ で展開すると
$$
m_{\rm th}^2=\frac{\beta_H}{2\pi^2\alpha'^2}\delta
+\frac{\delta^2}{4\pi^2\alpha'^2}
$$
である。$\beta>\beta_H$ では正、$\beta<\beta_H$ の近くで負になる。一ループの長管のその二つの巻き付きの寄与は、運動量の連続次元を $d$ として
$$
\int_{s_0}^\infty\frac{\dd s}{s}\,(4\pi s)^{-d/2}e^{-s m_{\rm th}^2}
$$
に比例する。$m_{\rm th}^2\to0^+$ において下限を変えると、$m_{\rm th}^2$ について解析的な項だけが変わる。上端が作る非解析項は、$d/2$ が整数でなければ $(m_{\rm th}^2)^{d/2}$、正の整数なら $(m_{\rm th}^2)^{d/2}\log m_{\rm th}^2$ である。非整数の場合は $u=sm_{\rm th}^2$ として、下端で積分できない Taylor 項を有限個引けば $\Gamma(-d/2)$ の項が残る。整数の場合にはその Gamma 関数の極と解析的な引き算が合わさって対数を残す。例えば $d=25$ では $\delta^{25/2}$ となり、十二回微分まで有限で十三回目で発散する。これは上の状態密度からのモーメント判定と一致する。

この長管には $w=0$ のボソン弦のタキオンは含めていない。$m_{\rm th}^2$ が負になることは、もとの場の理論で温度を上げたときの正則な相転移を自動的に与えるわけではない。熱円の近くの背景、弦の相互作用、重力による不安定性を調べるには自由な一ループを越える必要がある。Type II では零温度の NS 真空の射影を [超弦理論](/physics/string-theory/superstrings/) で導き、熱円のスピン構造を [Type I 弦とヘテロティック弦](/physics/string-theory/heterotic-type-i/) の一ループの方法で区別する。

## 演習問題

::::exercise[二形式の強度]{#exr-string-effective-action-1 level=basic}
$B=bX^1\dd X^2\wedge\dd X^3$ の $H_{123}$ を求めよ。

:::solution
外微分は $X^1$ の微分だけに作用して $H=b\dd X^1\wedge\dd X^2\wedge\dd X^3$ となるので $H_{123}=b$ である。
:::
::::

::::exercise[計量の補正]{#exr-string-effective-action-2 level=basic}
$H=0$ かつ $\Phi$ 一定とする。三次元球面で $L^2=100\alpha'$ の場合に、最低次数の計量ベータ関数と計量の係数比を求めよ。

:::solution
$R_{ij}=2G_{ij}/L^2$ より $\bar\beta^G_{ij}=2\alpha'G_{ij}/L^2$ である。比は $1/50$ となる。
:::
::::

::::exercise[次元とフレーム]{#exr-string-effective-action-3 level=basic}
$D=6$ で $G=e^{a(\Phi-\Phi_0)}G_E$ の $a$ を求めよ。

:::solution
曲率と体積の組は $e^{(D-2)a(\Phi-\Phi_0)/2}$ なので $e^{-2\Phi}$ を相殺するには $(D-2)a/2=2$ である。従って $a=1$ である。
:::
::::

::::exercise[中心電荷]{#exr-string-effective-action-4 level=basic}
$D=24$、$\alpha'V^2=1/3$ の物質中心電荷を求めよ。

:::solution
$24+6(1/3)=26$ である。ゴーストの $-26$ と相殺する。
:::
::::

::::exercise[ディラトンの変分]{#exr-string-effective-action-5 level=advanced}
:eqref[eq-string-nsns-action] のディラトン式を求めよ。

:::solution
前置因子の変分は $-2(R+4|\nabla\Phi|^2-H^2/12)$、勾配の変分は部分積分して $-8(\Box\Phi-2|\nabla\Phi|^2)$ となる。和を $-2$ で割ると $R+4\Box\Phi-4|\nabla\Phi|^2-H^2/12=0$ となる。
:::
::::

::::exercise[三形式の Weyl 重み]{#exr-string-effective-action-6 level=advanced}
$D=6$ で Einstein フレームの $H_E^2$ の前置因子を求めよ。

:::solution
$G=e^{\Phi-\Phi_0}G_E$ なので体積とディラトンで $e^{\Phi-\Phi_0}$、逆計量三個で $e^{-3(\Phi-\Phi_0)}$ となる。従って $e^{-2(\Phi-\Phi_0)}$ である。
:::
::::

::::exercise[Fourier の殻]{#exr-string-effective-action-7 level=advanced}
運動量の大きさが $\Lambda/2$ から $\Lambda$ の縮約を求めよ。

:::solution
角積分をすると $2\pi\alpha'\int_{\Lambda/2}^{\Lambda} k\dd k/(2\pi k^2)$ であり、$\alpha'\log2$ となる。
:::
::::

::::exercise[熱的な巻き付き]{#exr-string-effective-action-8 level=advanced}
ボソン閉弦で $\beta=2\beta_H$ のとき、巻き付き一つ、運動量と振動子励起がゼロの状態の $\alpha'M^2$ を求めよ。

:::solution
$\beta=8\pi\ell_s$ より円半径は $4\ell_s$ である。巻き付きの寄与 16 と真空の $-4$ を足すと 12 である。
:::
::::

::::exercise[境界のエネルギーモーメント]{#exr-string-effective-action-9 level=advanced}
非コンパクト空間が $d=4$ で、一弦の高質量密度を $M^{-a}e^{\beta_H M}$ と仮定する。$\beta=\beta_H$ でエネルギーの二次モーメントが有限となる条件を求め、$a=5$ の場合を判定せよ。

:::solution
動径運動量でなく全運動量の Gaussian から $M^{d/2}=M^2$ が掛かる。$M\beta_H\gg1$ では典型的な運動エネルギーは $M$ より小さく、$E^2=(M+O(1))^2$ の先頭は $M^2$ である。従って二次モーメントでの冪は $M^{4-a}$ であり、$4-a<-1$ すなわち $a>5$ が収束条件である。$a=5$ では跡は $\int M^{-3}\dd M$、一次は $\int M^{-2}\dd M$ でともに高質量側で有限だが、二次は $\int M^{-1}\dd M$ で対数発散する。
:::
::::

::::exercise[流束の量子化と中心電荷の差]{#exr-string-effective-action-10 level=advanced}
三球の $H$ の量子数を $k=6$ とする。二微分の背景方程式から得る半径 $L/\ell_s$ と、残り 23 方向の一つに必要な勾配の $\alpha'V^2$ を求めよ。また三球をレベル 6 の $SU(2)$ 共形場理論として扱ったときの値を比較せよ。

:::solution
$|h|=2/L$ と $\int_{S^3}H/(4\pi^2\alpha')=hL^3/(2\alpha')=6$ から $L^2=6\alpha'$ である。この近似では $H^2=24/L^2$ なので臨界の $\bar\beta^\Phi=0$ に必要な勾配は $V^2=1/L^2$、すなわち $\alpha'V^2=1/6$ である。共形場理論では三球側の中心電荷が $3\cdot6/(6+2)=9/4$ である。従って $23+9/4+6\alpha'V^2=26$ から $\alpha'V^2=1/8$ になる。両値の差は $1/24$ であり、$k=6$ の有限レベルで二微分の幾何学だけから中心電荷を代用できないことが分かる。
:::
::::

## 確認テスト

::::quiz{answer="26" tol="1e-6"}
$D=14$、$\alpha' V^2=2$ の $c$ を求めよ。

:::explanation
平坦なボソンから 14、改善項の自己縮約から $6\cdot2=12$ であり、和は 26 である。
:::
::::

::::quiz
一定のディラトンが与えるものを選べ。

- [x] Euler 標数に応じた結合の重み
- [ ] 同じ種数での $X$ の質量項
- [ ] 世界面並進の破れ

:::explanation
作用は $\Phi_0\chi$ なので、世界面の位相によって決まる。局所的な $X$ の変分はゼロである。
:::
::::

::::quiz{answer="4" tol="1e-6"}
ボソン弦の $\beta_H/(\pi\ell_s)$ を求めよ。

:::explanation
状態密度の指数は $4\pi\ell_s M$ である。逆温度の境界はその係数であり、$\pi\ell_s$ で割ると 4 になる。
:::
::::

::::quiz
有効作用の微分展開に必要な条件を選べ。

- [x] $\alpha'/L^2\ll1$
- [ ] $\alpha'/L^2\gg1$
- [ ] 任意の曲率半径で二微分項が厳密

:::explanation
曲率の桁は $L^{-2}$ であり、微分の高い補正は $\alpha'$ と追加の曲率や微分を含む。上の比が小さいと抑制される。
:::
::::

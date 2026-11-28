---
title: CP 対称性の破れと K・B 中間子
description: "中性中間子の有効 Hamiltonian を解き、混合・崩壊・干渉に現れる CP の破れを区別する。"
notation:
  - symbol: 'M_{12},\Gamma_{12}'
    meaning: "中性中間子の分散的・吸収的混合"
  - symbol: 'p,q'
    meaning: "中性中間子の固有状態の混合係数"
  - symbol: '\lambda_f'
    meaning: "混合と共通終状態への崩壊の振幅比"
  - symbol: '\bar\theta'
    meaning: "QCD に残る物理的な CP 位相"
---

## この章の目標

- 中性中間子の固有状態と時間発展を求める。
- K 中間子の二種類の CP の破れを分離する。
- B 中間子の時間依存非対称性を導く。
- 重いクォークの展開と強い CP 問題の仮定を理解する。

以下では自然単位系 $\hbar=c=1$ と計量 $g=\diag(1,-1,-1,-1)$ を用いる。

## 中性 K 中間子の混合と CP の破れ

### 非 Hermite の二状態問題

粒子と反粒子をまとめて $i\,\dd\psi/\dd t=(M-i\Gamma/2)\psi$ とする。$M,\Gamma$ は Hermite 行列であり、崩壊チャネルを除いた状態のノルムの減少が $\Gamma$ を与える。CPT を仮定すれば対角要素は共通 $H_0$ である。非対角を $H_{12}=M_{12}-i\Gamma_{12}/2$、$H_{21}=M_{12}^*-i\Gamma_{12}^*/2$ とする。固有値方程式は $(\mu-H_0)^2=H_{12}H_{21}$ であり、固有ベクトルは $p\ket{K^0}\pm q\ket{\bar K^0}$ で
$$
\mu_\pm=H_0\pm\sqrt{H_{12}H_{21}},\qquad (q/p)^2=H_{21}/H_{12}
\label{eq-sm-neutral-mixing}
$$
となる。ここでは $H_{12}H_{21}\ne0$ を仮定し、固有値の枝と $q/p$ の符号を同じ固有ベクトル方程式から選ぶ。両方が零なら $H=H_0 I$ であり、一方だけが零なら固有ベクトルが一つなので二つの指数への分解は使えない。$\mu=m-i\Gamma/2$ と分解すれば各指数の振動と減衰が分かる。二状態近似は崩壊チャネルを消去して遅い変化だけを残したものである。

:term[CP 変換]{en="charge conjugation and parity transformation" yomi="しーぴーへんかん"} は荷電共役と空間反転の積である。規約を $\mathrm{CP}\ket{K^0}=\ket{\bar K^0}$ と選ぶ。CP が保存すれば $H_{12}=H_{21}$ であり、$\ket{K_{1,2}}=(\ket{K^0}\pm\ket{\bar K^0})/\sqrt2$ は CP の固有値 $\pm1$ を持つ。二つの π 中間子の $J=0$ 状態は CP が正なので $K_1$ が主に速く崩壊し、$K_2$ は位相空間の小さい三 π 中間子などを通る。実際の短寿命 $K_S$、長寿命 $K_L$ では $p\ne q$ を許す。$K_L$ の小さい二 π 中間子 分岐が CP の破れを調べる手掛かりとなる。

### 吸収部分と共通終状態

初期の二状態を射影 $P$、崩壊を含む連続状態を $Q=1-P$ とする。全 Hamiltonian の固有方程式で $Q$ 部分を形式的に解くと $H_{\rm eff}(E)=PHP+PHQ(E-QHQ+i0)^{-1}QHP$ である。エネルギー分母を主値と $-i\pi\delta(E-QHQ)$ に分けると、前者は分散的な補正、後者が $\Gamma_{ij}=2\pi\sum_f\langle i|H|f\rangle\langle f|H|j\rangle\delta(E-E_f)$ を与える。この行列は任意の係数 $c_i$ について $\sum_{ij}c_i^*\Gamma_{ij}c_j$ が振幅の絶対値二乗の和になるため半正定値である。非対角の吸収部分には、二つの状態が同じ終状態に崩壊できることが必要である。分母を共通質量の近くで評価してゆっくり変化する部分を残した結果が、章の初めの二状態方程式である。

### 固有状態が直交しないことを幅から調べる

非 Hermite な $H=M-i\Gamma/2$ の右固有状態を $|S\rangle,|L\rangle$ とすると、一般にその二つは直交しない。$\langle S|H|L\rangle=\mu_L\langle S|L\rangle$ と $\langle S|H^\dagger|L\rangle=\mu_S^*\langle S|L\rangle$ の差を取れば
$$
\left[\frac{\Gamma_S+\Gamma_L}{2}+i(m_L-m_S)\right]\langle S|L\rangle
=\langle S|\Gamma|L\rangle
$$
となる。右辺は共通終状態についての $A_S^*A_L$ の位相空間積分である。従って同じような質量の二状態の重なりは、共通崩壊チャンネルの振幅の重なりを測ることでも制約できる。CPT が成り立つ場合の $S,L=p|K^0\rangle\pm q|\bar K^0\rangle$ では $\langle S|L\rangle=|p|^2-|q|^2$ である。半幅の和と質量差の両方を持つ複素係数のため、右辺の位相も重要である。

なお生産時の $K^0$ は $S$ だけでも $L$ だけでもない。ある時刻の二 π 中間子の振幅には、速い $S$ の指数と、幅の狭い $L$ の小さな崩壊振幅の積が足される。十分遅い時刻には前者が抑えられ、長寿命成分だけを近似的に選べる。時間を分けずに同じ二 π 中間子 終状態を観測する場合は、二振幅の干渉まで考える必要がある。

### 物質で短寿命成分を再生する

物質中での前方散乱は、ストレンジネスを交換しない近似では $K^0,\bar K^0$ の対角に異なる複素エネルギー $\chi,\bar\chi$ を加える。実部は同位相の前方振幅、虚部は非弾性や散乱によるビームの減衰を含む。CP の破れを無視した $\ket{K_S}=(\ket{K^0}+\ket{\bar K^0})/\sqrt2$ と $\ket{K_L}=(\ket{K^0}-\ket{\bar K^0})/\sqrt2$ を使うと、この項は
$$
\tfrac12(\chi+\bar\chi)I+
\tfrac12(\chi-\bar\chi)\begin{pmatrix}0&1\\1&0\end{pmatrix}
$$
となる。従って $\chi\ne\bar\chi$ の物質を通ると、純粋な $\ket{K_L}$ からも短寿命の成分を作れる。これは媒質が粒子と反粒子に同じ条件を与えない効果であり、真空の Hamiltonian の CP の破れを仮定しない。

薄い領域を通る固有時間を $\delta t$ とする。真空での固有値差と媒質項の積がともに $\delta t$ で小さいとき、共通位相 $e^{-i\mu_L\delta t}$ を分離すると、残る時間発展は $I-i(H_{\rm eff}-\mu_L I)\delta t+\order{\delta t^2}$ である。$\ket{K_L}$ 入射の短寿命の係数は $-i(\chi-\bar\chi)\delta t/2$ となる。領域の後に真空を時間 $t$ だけ進むと、二パイオンの振幅は、再生係数を $\mathcal R$ として $A_S[\mathcal R e^{-i\mu_S t}+\eta e^{-i\mu_L t}]$ になる。干渉の位相は $m_L-m_S$ の振動を持ち、短寿命の指数が小さくなる時刻では $\eta$ の寄与だけになりやすい。

:::example[再生の振幅と率の桁]
$(\chi-\bar\chi)\delta t=0.02$ の実数を仮定すると、出口での短寿命の振幅は $-0.01i$、生成確率は先頭で $10^{-4}$ である。$\eta=0.002i$ を同じ規約で加えると、出口直後の二π振幅を $A_S$ で割ったものは $-0.008i$、その二乗は $6.4\times10^{-5}$ になる。二つの確率を別々に足した $1.04\times10^{-4}$ との差が干渉である。
:::

### 計算例

:::example[寿命を比較する]
二固有状態の幅 $G_S,G_L$ を同じ単位で与え、寿命の比 $R=\tau_L/\tau_S$ を求める。入力値を次のようにとる。

$$
G_{S}=600,\qquad G_{L}=1
$$

与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
R&=\frac{G_{S}}{G_{L}}=600
\end{aligned}
$$

寿命と幅は逆数の関係にある。
:::

:::example[純粋に分散的な混合]
$\Gamma_{12}=0$ とし、$M_{12}$ の大きさ $m$ を MeV で与え、質量差 $D$ を求める。入力値を次のようにとる。

$$
m=\frac{1}{1000}
$$

与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
D&=2 m=\frac{1}{500}\simeq 0.002\quad[\mathrm{MeV}]
\end{aligned}
$$

K 系で吸収的な混合を実際に落とせるという仮定ではない。
:::

:::example[固有状態の重なり]
規格化された $p$ の二乗の大きさ $P$ を与え、CPT を仮定した二固有状態の重なり $D$ を求める。入力値を次のようにとる。

$$
P=\frac{251}{500}
$$

与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
D&=2 P - 1=\frac{1}{250}\simeq 0.004
\end{aligned}
$$

共通の位相規約で $q$ の大きさの二乗に $1-P$ を使った。
:::

## ε と ε′ パラメータ

### 二つの π 中間子の終状態を比較する

$K_S\propto K_1+\epsilon_m K_2$、$K_L\propto K_2+\epsilon_m K_1$ と書くと $\epsilon_m=(p-q)/(p+q)$ である。崩壊振幅をアイソスピン $I=0,2$ で $A_I e^{i\delta_I}$ と書き、$\delta_I$ は強い散乱位相とする。CPT の下で反粒子には $A_I^*e^{i\delta_I}$ が現れる。和の状態は $\sqrt2\Re A_I e^{i\delta_I}$、差の状態は $i\sqrt2\Im A_I e^{i\delta_I}$ を与える。従って一次の共通係数は $\epsilon=\epsilon_m+i\Im A_0/\Re A_0$ である。$\omega=\Re A_2/\Re A_0$ と置く。Clebsch–Gordan 係数から荷電対は $A_0e^{i\delta_0}+A_2e^{i\delta_2}/\sqrt2$、中性対は $A_0e^{i\delta_0}-\sqrt2A_2e^{i\delta_2}$ に比例する。分母を $\omega$ の一次まで展開すると
$$
\begin{aligned}
\eta_{+-}&=\epsilon+\epsilon',&\eta_{00}&=\epsilon-2\epsilon',\\
\epsilon'&=\frac{i\omega}{\sqrt2}e^{i(\delta_2-\delta_0)}
\left(\frac{\Im A_2}{\Re A_2}-\frac{\Im A_0}{\Re A_0}\right).
\end{aligned}
$$
$\eta$ は $K_L$ 振幅を $K_S$ 振幅で割った値である。この表示ではアイソスピン破れ、CP 破れの高次、$\omega$ の二次以上を無視した。$\eta$ の比を二乗すると $|\eta_{00}/\eta_{+-}|^2\simeq1-6\Re(\epsilon'/\epsilon)$ である。$\epsilon'$ が非零なら共通の混合だけでは二チャンネルを説明できない。強い行列要素を含む定量的な標準模型予測には非摂動 QCD と弱 Hamiltonian の繰り込みが要るため、:cite[donoghue-dsm] を参照する。

混合だけに由来する一次の $\epsilon_m$ も有効 Hamiltonian から求められる。CP 基底の非対角要素 $\langle K_1|H|K_2\rangle=(H_{21}-H_{12})/2$ に冒頭の定義を代入すると
$$
\epsilon_m=\frac{-i\Im M_{12}-\tfrac12\Im\Gamma_{12}}
{(m_L-m_S)+i(\Gamma_S-\Gamma_L)/2}
$$
となる。これは $K_L=K_2+\epsilon_m K_1$ の第一成分の式を一次で解いた値である。$A_0$ を実にした規約で $\Gamma_{12}$ をほぼ二パイオンの $I=0$ 成分で飽和させると、吸収部の虚部を落とせる。このとき $\epsilon\simeq-\Im M_{12}\,e^{i\phi_K}\sin\phi_K/(m_L-m_S)$ であり、$\tan\phi_K=2(m_L-m_S)/(\Gamma_S-\Gamma_L)$ である。分子の負号は $M_{12}=\langle K^0|M|\bar K^0\rangle$ の順序による。別の位相規約では崩壊側の $i\Im A_0/\Re A_0$ も合わせる。

### 混合の位相と崩壊の位相を移す

フレーバー状態を $|K^0\rangle\mapsto e^{i\xi}|K^0\rangle$、$|\bar K^0\rangle\mapsto e^{-i\xi}|\bar K^0\rangle$ とすると、混合比 $q/p$ は $e^{2i\xi}$ を受け、同じ終状態への振幅比 $\bar A_f/A_f$ は $e^{-2i\xi}$ を受ける。従って積は保たれる。$\epsilon_m$ 自体を絶対的な観測量とみなせないのはこのためである。$I=0$ 振幅を実にする規約を選ぶと、一次の $\epsilon$ における $i\Im A_0/\Re A_0$ は消え、その寄与は混合の位相の定義へ移る。二つの終状態の比 $\eta_{00}/\eta_{+-}$ には共通位相の変更が打ち消され、直接の破れに感度が残る。

$\epsilon'/\epsilon$ を純虚数と仮定した場合、二重比の一次変化は消える。しかし振幅の位相差の測定には虚部が残り、二重比だけで複素比全体を決めることはできない。有限の比に対しては $|1-2z|^2/|1+z|^2$ を使い、$z=\epsilon'/\epsilon$ の二次以上を制御する。

### 二 π 中間子の係数をアイソスピンから作る

二つの $I=1$ の π 中間子で、$I_3=0$、交換に対称な基底を $|c\rangle=(|+-\rangle+|-+\rangle)/\sqrt2$ と $|n\rangle=|00\rangle$ とする。一重項の位相を $|0,0\rangle=(\sqrt2|c\rangle-|n\rangle)/\sqrt3$、五重項を $|2,0\rangle=(|c\rangle+\sqrt2|n\rangle)/\sqrt3$ と取る。前者が $I_+$ によって消えること、両者が規格直交することは昇降係数から直接確認できる。従って $c$ の振幅は $\sqrt{2/3}A_0+\sqrt{1/3}A_2$、$n$ では $-\sqrt{1/3}A_0+\sqrt{2/3}A_2$ である。各チャンネルで共通の全体係数をくくると $\epsilon'$ の導出で使った $1/\sqrt2$ と $-\sqrt2$ を得る。振幅を幅へ移すときは識別しない中性二粒子の位相空間の規格化と、電荷の順序を採用した規約をそろえる。ここでの $\eta$ は各チャンネル内で同じ外部状態同士を割るため、全体係数が消える。

アイソスピン零の振幅が大きいことは $|\omega|\ll1$ に表れる。直接の CP 破れは弱位相差だけでなく、二つのアイソスピン成分の相対重みを含む。$\eta_{+-}$ と $\eta_{00}$ の差を測るときは、小さな $\eta$ の絶対値の比を取るため背景の寄与を適切に制御する必要がある。一つのモードの分岐比が小さいことだけから $\epsilon'$ の値は決まらない。

### 共役の半レプトン終状態で混合の大きさを調べる

荷電カレントの木の過程で $\Delta S=\Delta Q$ を仮定すると、$K^0=d\bar s$ は正の荷電レプトン、$\bar K^0=s\bar d$ は負の荷電レプトンを含む半レプトン終状態に対応する。両方の共役した直接振幅の大きさが等しいと仮定し、$K_L=p|K^0\rangle-q|\bar K^0\rangle$ の半レプトンの正負を比べると、差を和で割った量は
$$
\delta_L=\frac{|p|^2-|q|^2}{|p|^2+|q|^2}
=\frac{2\Re\epsilon_m}{1+|\epsilon_m|^2}
$$
である。二番目の式は $p\propto1+\epsilon_m$、$q\propto1-\epsilon_m$ から従う。二 π 中間子で求める量には各終状態の弱位相が入り、半レプトンの大きさの比較は混合の絶対値への別の情報を与える。$\Delta S\ne\Delta Q$ の振幅や直接の半レプトン率の CP 差も許す分析では、この抽出の前提を変更する必要がある。

### 二重比で一次の直接 CP を取り出す

同じ終状態 $f$ について、部分幅比 $\Gamma(K_L\to f)/\Gamma(K_S\to f)$ を作る。固有状態の微小な質量差による位相空間の違いを無視すれば、これは $|\eta_f|^2$ である。各試料の分岐比だけを用いる場合は全幅の比も掛ける。そこで二つの $f$ についてこの商をさらに割れば全幅の比が消える。$\eta_{00}\simeq\epsilon-2\epsilon'$、$\eta_{+-}\simeq\epsilon+\epsilon'$ を入れて
$$
D=\frac{|\eta_{00}|^2}{|\eta_{+-}|^2}\simeq
\frac{|1-2r|^2}{|1+r|^2},\qquad r=\epsilon'/\epsilon
$$
となる。このパラメータ化の範囲で分子は $1-4\Re r+4|r|^2$、分母は $1+2\Re r+|r|^2$ である。従って $D=1-6\Re r+\order{|r|^2}$ であり、共通する間接 CP の大きさが消え、比の実部を測れる。逆順の二重比を使えば一次の符号が変わる。

:::example[複素位相の違いを区別する]
$r=0.002$ なら $D=(0.996/1.002)^2\simeq0.98806$、一次では $0.988$ である。$r=0.002i$ なら一次では 1 で、比の模型を二次まで保持した値は $(1+0.000016)/(1+0.000004)\simeq1.000012$ である。二重比が 1 に近いことだけでは、比の虚部が小さいとは結論できない。
:::

:::example[同じ全幅の換算を二度行わない]
観測した分岐比の商が $B_L^{00}/B_S^{00}=0.003$、$B_L^{+-}/B_S^{+-}=0.00303$ なら二重比は $100/101$ である。部分幅の換算でどちらにも同じ $\Gamma_L/\Gamma_S$ が掛かるので、二重比ではその数値を追加しない。受理率や放射を含む終状態の定義はチャンネルで違い得るため、それらの補正は先に行う。
:::

### 計算例

:::example[二重比から直接の破れへ]
幅の二重比 $R$ から一次近似の $\Re(\epsilon'/\epsilon)=r$ を求める。入力値を次のようにとる。

$$
R=\frac{99}{100}
$$

与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
r&=\frac{1}{6} - \frac{R}{6}=\frac{1}{600}\simeq 0.00166667
\end{aligned}
$$

幅の比を用いるため、振幅の比をさらに二乗した結果である。
:::

:::example[二つの複素振幅比の差]
$\epsilon'$ の実部を $e$ とし、$\Re(\eta_{+-}-\eta_{00})=d$ を求める。入力値を次のようにとる。

$$
e=\frac{1}{250000}
$$

与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
d&=3 e=\frac{3}{250000}\simeq 1.2\times10^{-5}
\end{aligned}
$$

共通の混合寄与は差を取ると消える。
:::

:::example[二重比の二次効果]
$z=x+iy$ を指定し、二重比 $R$ を展開せずに求める。入力値を次のようにとる。

$$
x=\frac{1}{500},\qquad y=\frac{1}{1000}
$$

与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
R&=\frac{4 y^{2} + \left(1 - 2 x\right)^{2}}{y^{2} + \left(x + 1\right)^{2}}=\frac{198404}{200801}\simeq 0.988063
\end{aligned}
$$

一次の実部が零になる場合でも二次の変化は残る。
:::

:::example[半レプトンの規格化差]
$\epsilon_m=a+ib$ を与え、共役の直接振幅の大きさが同じと仮定して $\delta_L=d$ を求める。入力値を次のようにとる。

$$
a=\frac{1}{500},\qquad b=\frac{1}{500}
$$

与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
d&=\frac{2 a}{a^{2} + b^{2} + 1}=\frac{500}{125001}\simeq 0.00399997
\end{aligned}
$$

混合の位相の全てをこの一つの率から決めることはできない。
:::

## B 中間子の混合と振動

### 混合確率の時間依存

規約を固有状態 $\ket{B_L}=p\ket B+q\ket{\bar B}$、$\ket{B_H}=p\ket B-q\ket{\bar B}$ とする。逆変換により $\ket{B(t)}=g_+(t)\ket B+(q/p)g_-(t)\ket{\bar B}$、$g_\pm=(e^{-i\mu_Lt}\pm e^{-i\mu_Ht})/2$ となる。$\Delta m=m_H-m_L$、$\Delta\Gamma=\Gamma_L-\Gamma_H$ とし、$|q/p|=1$、$\Delta\Gamma=0$ の場合、
$$
P(B\to B;t)=\tfrac12e^{-\Gamma t}(1+\cos\Delta m t),\quad
P(B\to\bar B;t)=\tfrac12e^{-\Gamma t}(1-\cos\Delta m t).
$$
二式の和が生存率である。これを同じ幅で時間積分した混合率は $\chi=x^2/[2(1+x^2)]$、$x=\Delta m/\Gamma$ である。積分で $\int e^{-\Gamma t}\cos(\Delta mt)\dd t=\Gamma/(\Gamma^2+\Delta m^2)$ を使った。有限の幅差で $y=\Delta\Gamma/(2\Gamma)$ とすると、和が $e^{-\Gamma t}\cosh(y\Gamma t)$ になり、同じ計算で $\chi=(x^2+y^2)/[2(1+x^2)]$ となる。半レプトン終状態の電荷で崩壊時フレーバーを判定する場合には、生産時のフレーバーを別にタグ付けする。

### 質量差と幅差を同じ行列から求める

$\Delta m$ と $\Delta\Gamma$ の定義を上と同じにする。式 :eqref[eq-sm-neutral-mixing] の二固有値の差を二乗して、実部と虚部を分けると
$$
(\Delta m)^2-\tfrac14(\Delta\Gamma)^2=4|M_{12}|^2-|\Gamma_{12}|^2,\qquad
\Delta m\,\Delta\Gamma=-4\Re(M_{12}\Gamma_{12}^*)
$$
となる。虚部の符号は $\Delta\Gamma=\Gamma_L-\Gamma_H$ による。$|\Gamma_{12}/M_{12}|\ll1$ なら $\Delta m\simeq2|M_{12}|$ であり、幅差は相対位相の余弦を含む。混合比の大きさの 1 からのずれには相対位相の正弦が現れる。従って質量の振動周期だけでは $|q/p|$ を決められない。

またフレーバーを特定の電荷で判定する半レプトンの終状態で、二つの混合経路の率を比べる。対応する直接の半レプトン振幅の大きさが等しければ、$B\to\bar B$ と逆方向の率はそれぞれ $|q/p|^2|g_-|^2$、$|p/q|^2|g_-|^2$ に比例する。従って時間に依存する共通因子が比で消えるので、吸収部の相対位相を別に調べられる。

### 対で生産する B 中間子を使う

電子陽電子の共鳴から中性 $B\bar B$ を反対方向に生成し、二中間子のフレーバー状態が反対称である場合を考える。一つの向きと逆の向きに付くラベルを 1,2 とし、初期状態を $(|B\rangle_1|\bar B\rangle_2-|\bar B\rangle_1|B\rangle_2)/\sqrt2$ とする。二つに同じ一粒子の時間発展行列 $U$ を掛けると、その反対称二階テンソルに $\det U$ が掛かる。従って同時刻でどちらかのフレーバーを測れば、残りはその反対に条件付けられる。初めの崩壊で片側をタグとし、二番目との崩壊時刻差 $\Delta t$ を使えば、一粒子の $g_\pm(\Delta t)$ による振動を調べられる。

$\Delta\Gamma=0$ とし、後で同じ電荷の半レプトンに崩壊した事象の割合を考えると、タグ直後からのフレーバー反転が対応する。同符号率と異符号率の差の時間依存は理想的には $\cos(\Delta m\Delta t)$ に比例する。タグが誤り確率 $w$ を持てば正しい分類の寄与と交換された寄与を引いて係数 $1-2w$ が掛かる。時刻差の分解能にも畳み込みがあり、早い振動を時間平均してしまう可能性がある。

### 小さい吸収混合を展開する

$z=\Gamma_{12}/M_{12}$ とすると、式 :eqref[eq-sm-neutral-mixing] から
$$
(q/p)^2=\frac{M_{12}^*}{M_{12}}\frac{1-iz^*/2}{1-iz/2}
=\frac{M_{12}^*}{M_{12}}[1-\Im z+\order{|z|^2}]
$$
となる。従って $|q/p|\simeq1-\Im z/2$ である。混合した二方向の率を、$\bar B\to B$ を正にして差を取れば $a_{\rm sl}=(1-|q/p|^4)/(1+|q/p|^4)\simeq\Im z$ である。$|q/p|$ と $\Delta\Gamma$ を同時に調べると、$\Gamma_{12}/M_{12}$ の虚部と実部にそれぞれ感度がある。後者は上の $\Delta m\Delta\Gamma=-4\Re(M_{12}\Gamma_{12}^*)$ を $\Delta m\simeq2|M_{12}|$ で割ると得られる。$M_{12}$ が非常に小さい場合にはこの展開を使わず、両非対角の固有値方程式を解く。

短距離の箱図から $B_q$ の混合を計算する場合、$B_q=\bar b q$ の規約で $M_{12}=\langle B_q|\mathcal H_{\Delta B=2}|\bar B_q\rangle/(2m_{B_q})$ と書く。トップ支配の近似では、演算子 $O_q=[\bar q\gamma_\mu(1-\gamma_5)b][\bar q\gamma^\mu(1-\gamma_5)b]$ の係数は $G_F^2 M_W^2 (V_{tq}^*V_{tb})^2 S_0/(16\pi^2)$ である。ここで $S_0=S_0(m_t^2/M_W^2)$ は式 :eqref[eq-sm-gim-inami-lim] であり、低尺度への強い発展を掛ける前の係数である。崩壊定数を $\langle0|\bar q\gamma^\mu\gamma_5 b|\bar B_q(p)\rangle=if_{B_q}p^\mu$ と取る。四クォーク行列要素の大きさを $8m_{B_q}^2 f_{B_q}^2 B_{B_q}/3$ としてバッグパラメータ $B_{B_q}$ を定義する。すると $\Delta m_q\simeq2|M_{12}|$ は $m_{B_q} f_{B_q}^2 B_{B_q} |V_{tq}^*V_{tb}|^2$ に比例する。

$B_d$ と $B_s$ の比ではトップの関数と共通の強い発展が相殺し、$\Delta m_d/\Delta m_s\simeq (m_{B_d}/m_{B_s})\xi_B^{-2}|V_{td}/V_{ts}|^2$ となる。ここで $\xi_B=f_{B_s}\sqrt{B_{B_s}}/(f_{B_d}\sqrt{B_{B_d}})$ である。従って CKM の比を取り出すには、二つの行列要素を同じ低エネルギー規約で計算して $\xi_B$ を決める。

### 時刻積分と有限な時間分解能

$|q/p|=1$、幅差を無視し、生成時にフレーバーが分かっている試料を考える。崩壊時間で全事象を積分した混合割合は
$$
\chi=\Gamma\int_0^\infty e^{-\Gamma t}\sin^2(\Delta m\,t/2)\,\dd t
=\frac{x^2}{2(1+x^2)},\qquad x=\Delta m/\Gamma. \label{eq-sm-cp-integrated-mixing}
$$
余弦の積分に $\Re[1/(\Gamma-i\Delta m)]$ を用いればよい。$x\ll1$ では $\chi\simeq x^2/2$、$x\gg1$ では $\chi\to1/2$ である。後者は瞬間の振動が消えたのではなく、半周期よりはるかに長い範囲で平均している。

一方、測る時刻 $t_{\rm obs}=t+\delta$ に平均零の Gaussian の幅 $\sigma_t$ があるとする。分布の包絡がゆっくり変化し、端点の影響を無視できる範囲では、振動因子 $e^{i\Delta m t}$ の畳み込みに $\langle e^{-i\Delta m\delta}\rangle=e^{-(\Delta m\sigma_t)^2/2}$ が掛かる。Gaussian の指数を平方完成すればこの特性関数が出る。時刻零の近くや、包絡を無視できない場合は確率密度全体を畳み込む。

:::example[速い振動を積分と時刻分布で比べる]
$x=3$ では $\chi=9/20=0.45$ である。$\Delta m\sigma_t=1$ なら時刻分布の振幅は $e^{-1/2}\simeq0.607$ 倍になる。タグの誤り率が $w=0.2$ で時刻誤差と独立なら、さらに $1-2w=0.6$ が掛かり、振幅の係数は約 $0.364$ になる。積分率と振幅の減少は異なる観測量である。
:::

### 計算例

:::example[時間積分した混合]
大きさの比 $x=\Delta m/\Gamma$ を与え、幅差を無視した $\chi$ を得る。入力値を次のようにとる。

$$
x=\frac{77}{100}
$$

与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
\chi&=\frac{x^{2}}{2 x^{2} + 2}=\frac{5929}{31858}\simeq 0.186107
\end{aligned}
$$

非常に速い振動を時刻を区別せず測ると混合率は半分に近い。
:::

:::example[時間分解した混合]
位相 $\phi=\Delta mt$ と生存因子 $S=e^{-\Gamma t}$ を与え、混合確率 $P$ を求める。入力値を次のようにとる。

$$
\phi=\pi,\qquad S=\frac{1}{2}
$$

与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
P&=\frac{S \left(1 - \cos{\left(\phi \right)}\right)}{2}=\frac{1}{2}
\end{aligned}
$$

生存した事象内での条件付き確率なら $S$ で割る。
:::

:::example[二方向の混合の比]
混合係数比の大きさ $r=\abs{q/p}$ を与え、$B\to\bar B$ と逆方向の率比 $R$ を求める。入力値を次のようにとる。

$$
r=\frac{101}{100}
$$

与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
R&=r^{4}=\frac{104060401}{100000000}\simeq 1.0406
\end{aligned}
$$

直接の半レプトン振幅の絶対値が等しい近似を置いた。
:::

:::example[タグの希釈]
タグ誤り確率 $w$ と真の振動コントラスト $A$ から、観測する振幅の係数 $D$ を求める。入力値を次のようにとる。

$$
w=\frac{1}{5},\qquad A=\frac{9}{10}
$$

与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
D&=A \left(1 - 2 w\right)=\frac{27}{50}
\end{aligned}
$$

時刻の分解能と背景の効果はここでは落とした。
:::

:::example[混合だけの半レプトン非対称]
$r=\abs{q/p}$ を与え、正の向きを本文のように取った $a$ を求める。入力値を次のようにとる。

$$
r=\frac{999}{1000}
$$

与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
a&=\frac{1 - r^{4}}{r^{4} + 1}=\frac{3994003999}{1996005996001}\simeq 0.002001
\end{aligned}
$$

$\Im z$ の展開を使わず二方向の率を比べた値である。
:::

## B 中間子崩壊における CP の破れ

### 振幅比と非対称性

共通終状態 $f$ への振幅を $A_f,\bar A_f$ とする。:term[混合と崩壊の干渉]{en="interference between mixing and decay" yomi="こんごうとほうかいのかんしょう"} は $A_f(g_++\lambda_f g_-)$、$\lambda_f=(q/p)(\bar A_f/A_f)$ に現れる。上の規約で幅差を零とすると $g_-/g_+=i\tan(\Delta mt/2)$ である。$|q/p|=1$ を仮定し、$\bar B$ から始める率から $B$ からの率を引いて和で割る。交差項 $2\Re[i\lambda_f]\sin(\Delta mt/2)\cos(\Delta mt/2)$ を展開すると
$$
A_{\rm CP}(t)=S_f\sin\Delta mt-C_f\cos\Delta mt,\quad
S_f=\frac{2\Im\lambda_f}{1+|\lambda_f|^2},\quad C_f=\frac{1-|\lambda_f|^2}{1+|\lambda_f|^2}.
\label{eq-sm-cp-time}
$$
$\bar A_f/A_f$ にほぼ一種類の弱位相しかない CP 固有終状態では、共通の強位相が消えて $C_f$ は小さくなる。$B_d\to J/\psi K_S$ では標準的な CKM 規約で $\lambda_f\simeq-e^{-2i\beta}$、従って $S_f\simeq\sin2\beta$ である。中性 K の混合位相も含めた結果である。

直接の崩壊を二経路 $a_1e^{i(\delta_1+\phi_1)}+a_2e^{i(\delta_2+\phi_2)}$ とすれば、反粒子では弱位相 $\phi_i$ のみが反転する。率の差は、順を反粒子から粒子として $4a_1a_2\sin(\delta_1-\delta_2)\sin(\phi_1-\phi_2)$ となる。従って $C_f$ と $S_f$ が同じ位相組合せを測る必要はない。ペンギン振幅が競合する系では強位相と弱位相の両方を制御して解釈する。

### 幅差で見える振幅の実部

上と同じ $g_\pm$ では $g_+^*g_-=\tfrac12e^{-\Gamma t}[-\sinh(\Delta\Gamma t/2)+i\sin(\Delta mt)]$ である。これは $|e_L|^2-|e_H|^2$ と、二つの指数の干渉の差を計算した結果である。$|g_\pm|^2=\tfrac12e^{-\Gamma t}[\cosh(\Delta\Gamma t/2)\pm\cos(\Delta mt)]$ も合わせると、$B$ から始まる率に $\Re\lambda_f$ に比例する $-\sinh(\Delta\Gamma t/2)$ の項が出る。$|q/p|=1$ で二つの初期フレーバーを同数ずつ足すと正弦・余弦は消え、$\cosh$ とこの $\sinh$ が残る。したがって幅差の大きい系では生産時のフレーバーを測らなくても、寿命の形から $\Re\lambda_f$ を調べられる。$\Delta\Gamma=0$ の式 :eqref[eq-sm-cp-time] にはこの情報は独立には現れない。

例えば $|q/p|=1$ なら、タグなしの率は $e^{-\Gamma t}[\cosh(y\Gamma t)+A_{\Delta\Gamma}\sinh(y\Gamma t)]$ に比例し、$y=\Delta\Gamma/(2\Gamma)$、$A_{\Delta\Gamma}=-2\Re\lambda_f/(1+|\lambda_f|^2)$ である。$\cosh,\sinh$ を指数に戻すと、幅 $\Gamma_H$ と $\Gamma_L$ の成分の重みは順に $(1+A_{\Delta\Gamma})/2$ と $(1-A_{\Delta\Gamma})/2$ になる。この形で平均崩壊時刻を定義すれば、各成分の規格化積分は $1/\Gamma_{H,L}$、時刻付きの積分は $1/\Gamma_{H,L}^2$ である。例えば $\lambda_f=+1$ なら軽い固有状態だけが寄与して平均は $1/\Gamma_L$ となる。また定義から $C_f^2+S_f^2+A_{\Delta\Gamma}^2=1$ なので、時間分布から独立に取り出す係数の整合性を検査できる。

### 直接 CP のための二つの位相差

直接振幅を $A_1[1+r e^{i(\delta+\phi)}]$、共役を $A_1'[1+r e^{i(\delta-\phi)}]$ と書く。$|A_1'|=|A_1|$ とすれば、共役率から元の率を引いた非対称性は
$$
A_{\rm dir}=\frac{2r\sin\delta\sin\phi}{1+r^2+2r\cos\delta\cos\phi}
$$
である。分子は二つの余弦の差、分母は和から得た。$r\ll1$ では $2r\sin\delta\sin\phi$ の次数であるが、二つが同程度で分母が小さくなる場合は一次では扱えない。強位相は終状態の再散乱や吸収部分によって生じ得る。従って複素 CKM だけを使って未知の強位相を零としてしまうと、直接の率の差を予測できなくなる。一方時間依存の干渉は混合の振幅をもう一経路として使うため、一種類の直接崩壊振幅だけでも非対称性を作れる。

### 有限時刻から積分した CP の係数

幅差を無視し、最初の $B$ と $\bar B$ の数が等しいとする。式 :eqref[eq-sm-cp-time] を時刻全体で積分する場合、率の差を積分してから率の和で割る。$\int_0^\infty e^{-\Gamma t}\sin(\Delta mt)\dd t=\Delta m/(\Gamma^2+\Delta m^2)$、余弦の積分は $\Gamma/(\Gamma^2+\Delta m^2)$ である。無偏りの受入れなら $x=\Delta m/\Gamma$ を用いて $(S_f x-C_f)/(1+x^2)$ が得られる。$x\gg1$ では振動の正負が平均され、小さな時刻差を分解する測定よりも CP の干渉への感度が減る。二つの初期フレーバーで異なる時間の効率があると、同じ規格化を使えない。共通の時間切り取りでも正弦と余弦への重みは変わるため、単に時間の中央値で式を評価するとは異なる。

### 小さな第二振幅が混合位相をずらす

CP 固有状態の固有値を $\eta_f$ とし、支配的な振幅の弱位相を $\phi_1$ とする。強い位相差 $\delta$、弱い位相差 $\phi$、大きさの比 $r$ を用いて、$A=A_1 e^{i\phi_1}(1+r e^{i(\delta+\phi)})$、$\bar A=\eta_f A_1 e^{-i\phi_1}(1+r e^{i(\delta-\phi)})$ と置く。CP 共役で強い相互作用の再散乱位相は同じであり、弱い位相だけを反転した。すると $|q/p|=1$ の近似で
$$
\lambda_f=\lambda_0[1-2ir e^{i\delta}\sin\phi+\order{r^2}],\qquad
\lambda_0=\eta_f(q/p)e^{-2i\phi_1}
$$
となる。従って $\arg(\lambda_f/\lambda_0)=-2r\cos\delta\sin\phi+\order{r^2}$、$|\lambda_f|^2=1+4r\sin\delta\sin\phi+\order{r^2}$ である。$C=(1-|\lambda_f|^2)/(1+|\lambda_f|^2)$ なら $C=-2r\sin\delta\sin\phi+\order{r^2}$ である。一方 $S=2\Im\lambda_f/(1+|\lambda_f|^2)$ では一次の絶対値の変化が打ち消され、$\arg\lambda_f$ の変化が残る。

:::example[直接非対称が小さくても角がずれる]
$r=0.02$、$\phi=\pi/2$、$\delta=0$ とすると $|\lambda_f|=1$ が厳密に成り立つ。位相のずれは $-2\arctan(0.02)\simeq-0.040\,\mathrm{rad}$ である。逆に $\delta=\pi/2$ なら比は実の $(1+r)/(1-r)$ となり、位相はずれないが $C=-2r/(1+r^2)\simeq-0.040$ である。同じ大きさの第二振幅で異なる効果が出るので、時間分布の二成分を分けて測る。
:::

### 計算例

:::example[干渉係数の計算]
$\lambda_f=a+ib$ として $S$ を求める。入力値を次のようにとる。

$$
a=\frac{4}{5},\qquad b=\frac{3}{5}
$$

与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
S&=\frac{2 b}{a^{2} + b^{2} + 1}=\frac{3}{5}
\end{aligned}
$$

$\lambda_f$ の符号規約は固有状態の表記と合わせる。
:::

:::example[直接非対称性の係数]
$\lambda_f$ の絶対値を $r$ とし、余弦項の係数 $C$ を求める。入力値を次のようにとる。

$$
r=1
$$

与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
C&=\frac{1 - r^{2}}{r^{2} + 1}=0
\end{aligned}
$$

この式では混合の絶対値を 1 と仮定している。
:::

:::example[実部の規格化]
$\lambda_f=a+ib$ の実部と虚部を与え、$\sinh$ の規格化係数 $D$ を求める。入力値を次のようにとる。

$$
a=1,\qquad b=0
$$

与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
D&=- \frac{2 a}{a^{2} + b^{2} + 1}=-1
\end{aligned}
$$

幅差の符号を交換する規約では係数の符号を合わせる。
:::

:::example[二経路の最大位相差の例]
強位相差と弱位相差の両方を $\pi/2$ とし、強度比 $r$ から率の非対称性 $A$ を計算する。入力値を次のようにとる。

$$
r=\frac{1}{10}
$$

与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
A&=\frac{2 r}{r^{2} + 1}=\frac{20}{101}\simeq 0.19802
\end{aligned}
$$

分母の干渉が消えるが、分子の干渉は残る。
:::

:::example[時間を積分した非対称]
$x=\Delta m/\Gamma$ と係数 $S,C$ を与え、一様な時間の効率での非対称 $A$ を求める。入力値を次のようにとる。

$$
x=\frac{4}{5},\qquad S=\frac{7}{10},\qquad C=\frac{1}{20}
$$

与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
A&=\frac{- C + S x}{x^{2} + 1}=\frac{51}{164}\simeq 0.310976
\end{aligned}
$$

幅差は零で $\abs{q/p}=1$ とした。
:::

## 重いクォーク有効理論の概観

### 重い質量を運動量から分ける

:term[重いクォーク有効理論]{en="heavy quark effective theory" yomi="おもいくぉーくゆうこうりろん"} ではハドロン内の $p_Q=m_Qv+k$ に $v^2=1$、$|k|\sim\Lambda_{\rm QCD}\ll m_Q$ を仮定する。$Q=e^{-im_Q v\cdot x}(h_v+H_v)$ と書き、$\not v$ の固有値 $+1,-1$ に射影する。Dirac 方程式の小成分は $(2m_Q+iv\cdot D)H_v=i\not D_\perp h_v$ を満たす。これを作用に代入すると
$$
\mathcal L=\bar h_v iv\cdot D h_v+
\bar h_v i\not D_\perp(2m_Q+iv\cdot D)^{-1}i\not D_\perp h_v
$$
となる。$1/(2m_Q)$ に展開した二次項は共変微分の反交換部分から運動演算子、交換部分から色磁気演算子を与える。先頭項に $m_Q$ もスピン行列もないため、同じ速度の重いフレーバーとスピンを変えても軽い場の作用が不変である。従って重いスピンの二重項が縮退し、同じ軽い状態の超微細分裂は $1/m_Q$ で減る。ゼロ反跳で軽い状態が同じなら、その重なりが 1 なので先頭の Isgur–Wise 関数は $\xi(1)=1$ である。短距離係数と有限質量補正を含む崩壊の解析は :cite[donoghue-dsm] を参照する。

### ゼロ反跳の運動学

重いメソン $H(v)\to H'(v')\ell\nu$ で速度積を $w=v\cdot v'$ とする。$q=Mv-mv'$ なので $w=(M^2+m^2-q^2)/(2Mm)$ である。未来向きの単位速度二つの積は $w\ge1$ であり、等しいとき $q^2=(M-m)^2$ になる。これがゼロ反跳である。質量零のレプトンでは低い端が $q^2=0$ で、$1\le w\le(M^2+m^2)/(2Mm)$ となる。$\xi(1)=1$ は行列要素の規格化であるが、崩壊率には娘粒子の位相空間が掛かるため、端点でそのまま率が最大になるという意味ではない。有限の質量では形状因子を端近くから外挿し、短距離係数を同じ規約で入れる。

ゼロ反跳では、軽い自由度をそろえた基底で先頭のカレントが重いフレーバーまたはスピンを変える対称性生成子になる。規格化した基底状態への一次の状態補正は元の状態に直交するため、生成子を挟んだ重なりの変化は二次から始まる。さらに方程式 $iv\cdot D h_v=0$ と $P_+\gamma_\perp P_+=0$ を使うと、ゼロ反跳で規格化される成分への一次の微分カレントも消える。この保護は Luke の定理と呼ばれ、例えば $B\to D$ の $h_+(1)$ と $B\to D^*$ の $h_{A_1}(1)$ に適用される。ここで $h_+$ はベクトル電流の $(v+v')^\mu$ の係数、$h_{A_1}$ は軸性電流の $D^*$ の偏極ベクトルの係数を規格化した形状因子である。短距離の照合係数やゼロ反跳から離れた形状には別の補正が残る。演算子基底を含む一般の扱いは :cite[donoghue-dsm] を参照する。

### 有効理論の次の次数でスピンを戻す

$D_\perp^\mu=D^\mu-v^\mu(v\cdot D)$ とし、$i\not D_\perp i\not D_\perp=(iD_\perp)^2+(g_s/2)\sigma^{\mu\nu}G_{\mu\nu}^{\perp}$ の関係を使う。対称な $\gamma^\mu\gamma^\nu$ の部分が二次の運動項、反対称部分が $[D^\mu,D^\nu]=-ig_sG^{\mu\nu}$ から色磁気項を作る。重い静止系で後者はスピンと色磁場を結ぶ。二つの演算子は短距離照合で係数が変わり得るが、どちらも作用では $1/(2m_Q)$ の次数が最初である。運動項は同じ軽い状態で重いスピンの向きを変えても不変であり、スピン二重項の分裂を直接作るのは磁気項である。

軽い自由度の角運動量を $j_\ell$ とすれば、$\vb S_Q\cdot\vb j_\ell=[J(J+1)-3/4-j_\ell(j_\ell+1)]/2$ である。$j_\ell>0$ の二重項で $J=j_\ell\pm1/2$ を代入すると $j_\ell/2$ と $-(j_\ell+1)/2$ になる。多重度はそれぞれ $2j_\ell+2$ と $2j_\ell$ なので、二つをその重みで足せば消える。従ってスピン平均質量から重い質量と軽い束縛を推定する際には、一次の色磁気分裂を取り除ける。

### ゼロ反跳近傍で運動学と形状因子を分ける

擬スカラーから擬スカラーへの半レプトン崩壊で、軽いレプトンの質量を落とす。重いクォーク極限のベクトル電流は $\sqrt{M_BM_D}(v+v')^\mu\xi(w)$ と書ける。$w=v\cdot v'$、$\xi$ は同じ軽い自由度の重なりであり、先頭では $\xi(1)=1$ である。$J^\mu$ のうちレプトン対の全運動量 $q^\mu$ に平行な成分は、質量零のレプトン電流の保存で落ちる。残りを通常の形状因子 $f_+(q^2)(p+p')^\mu$ に書き換えるには、$v+v'=p/M_B+p'/M_D$ を $(p+p')$ と $(p-p')$ に分ければよい。従って $f_+=(M_B+M_D)\xi(w)/(2\sqrt{M_BM_D})$ である。

親の静止系では $|\vb p_D|=M_D\sqrt{w^2-1}$ である。$q$ の静止系で二レプトンを積分したテンソルは、電流の保存と回転不変性により $q^\mu q^\nu-q^2g^{\mu\nu}$ に比例する。$J$ と縮約した運動学因子は $|f_+|^2[(p+p')\cdot q]^2-q^2|f_+|^2(p+p')^2=4M_B^2|\vb p_D|^2|f_+|^2$ である。二体位相空間がさらに $|\vb p_D|$ を一個与え、$|\dd q^2/\dd w|=2M_BM_D$ は一定である。従って $\dd\Gamma/\dd w$ は $(w^2-1)^{3/2}|\xi(w)|^2$ に比例する。

:::example[規格化点で率が零になる場合]
$\xi$ を一定にした模型で $w=1.10$ と $1.05$ の率の比は $(0.21/0.1025)^{3/2}\simeq2.93$ である。形式上もっとも確かな形状因子の規格化は $w=1$ だが、そこで位相空間は閉じる。近傍の有限な $w$ から傾きを含めて外挿する必要がある。実際の重い質量では短距離照合と $1/m_Q$ の補正も掛かる。
:::

### 計算例

:::example[重い質量展開の小さい量]
残差尺度 $L$ と質量 $m$ を GeV で与えて $r=L/m$ を求める。入力値を次のようにとる。

$$
L=\frac{2}{5},\qquad m=\frac{24}{5}
$$

与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
r&=\frac{L}{m}=\frac{1}{12}
\end{aligned}
$$

質量が小さい系では高次の補正が相対的に大きい。
:::

:::example[超微細分裂を比較する]
同じ軽い場を仮定し、重い質量 $m_1,m_2$ から第 2 対第 1 の分裂比 $r$ を最低次で求める。入力値を次のようにとる。

$$
m_{1}=\frac{3}{2},\qquad m_{2}=\frac{9}{2}
$$

与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
r&=\frac{m_{1}}{m_{2}}=\frac{1}{3}
\end{aligned}
$$

短距離係数の走りはこの単純な比に入れていない。
:::

:::example[速度積を計算する]
親と娘の質量 $M,m$ を GeV、$q_2=q^2$ を $\mathrm{GeV}^2$ で与え、$w$ を求める。入力値を次のようにとる。

$$
M=\frac{132}{25},\qquad m=\frac{187}{100},\qquad q_{2}=8
$$

与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
w&=\frac{M^{2} + m^{2} - q_{2}}{2 M m}=\frac{233753}{197472}\simeq 1.18373
\end{aligned}
$$

$q_2$ の上端で速度積は 1 になる。
:::

:::example[軽い角運動量と分裂]
$j_\ell=j>0$ とモデルの一対の係数 $C$ を MeV で与え、二重項の質量差 $D$ を求める。入力値を次のようにとる。

$$
j=\frac{1}{2},\qquad C=140
$$

与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
D&=\frac{C \left(2 j + 1\right)}{2}=140\quad[\mathrm{MeV}]
\end{aligned}
$$

$C$ が重い質量の逆数に比例する範囲に限って比較する。
:::

## 強い CP 問題とアクシオン

### 物理的な位相と動的な緩和

QCD は $\theta g_s^2G\widetilde G/(32\pi^2)$ を含められる。局所では全微分だが、異なるトポロジーの場を積分すると相対位相を生む。軸性回転の測度は [カイラル対称性](/physics/standard-model/chiral-symmetry/) の一重項発散に従い $\theta$ を移す。同じ回転で複素質量行列の行列式の位相も移るため、規約をそろえると $\bar\theta=\theta+\arg\det(M_uM_d)$ が不変である。$G\widetilde G$ は色電場と色磁場の内積に比例し、パリティで符号が反転する。低エネルギーの中性子の電気双極子の小ささは、自然に大きくなり得るこの位相を強く制限する。弱カレントの CKM 位相とは別の量である。

:term[アクシオン]{en="axion" yomi="あくしおん"} $a$ は Peccei–Quinn 対称性の破れにより生じ、結合を $(\bar\theta+a/f_a)q_{\rm top}$ にする。真空エネルギー $E(\vartheta)$ の最小を $\vartheta=0$ とし、$\chi=E''(0)>0$ を取ると、最小で $\langle a\rangle/f_a=-\bar\theta$ となり、励起の質量二乗は $\chi/f_a^2$ である。$\chi$ はここでは軽いクォークを含む QCD の感受率である。他の相互作用がこの対称性を明示的に破る場合には最小がずれる。この解決がどの程度成り立つかは追加の破れの大きさに依存する。

### 小さい位相に対する真空の曲率

二つの軽いクォークが正の質量を持つとする。異常により $\theta$ を軽い質量の位相に移し、低エネルギーで $\phi_u+\phi_d=\theta$ の制約を使う。一つの枝の小さい位相について $-F^2B_0[m_u\cos\phi_u+m_d\cos\phi_d]$ の二次項は $\Sigma(m_u\phi_u^2+m_d\phi_d^2)/2$、$\Sigma=F^2B_0$ である。制約付き最小化で $m_u\phi_u=m_d\phi_d$ となるので、最小エネルギーは $\theta^2\Sigma/[2(1/m_u+1/m_d)]$ だけ増える。従って最低次の感受率は $\chi=\Sigma m_u m_d/(m_u+m_d)$ である。一つの質量を零に近付ければ曲率は消える。これは質量零のカイラル回転で $\theta$ を除けることと一致する。実際の CP 位相の制約をこの式だけから得ることはできず、中性子の電気双極子などに対する行列要素を要する。

### 有限の強い位相で真空のエネルギーを最小化する

二フレーバーの最低次で一枝の $\phi_d=\theta-\phi_u$ とし、$\phi_u$ の角度を変える。質量の和は
$$
m_u\cos\phi_u+m_d\cos(\theta-\phi_u)
=(m_u+m_d\cos\theta)\cos\phi_u+m_d\sin\theta\sin\phi_u .
$$
右辺は二成分のベクトルと単位ベクトルの内積なので、最大値は $\sqrt{m_u^2+m_d^2+2m_um_d\cos\theta}$ である。従って凝縮の係数 $\Sigma$ を掛けた負が真空エネルギーの角依存を与える。小さい $\theta$ で二次まで展開すると、前節で制約付き最小化を用いた $\chi=\Sigma m_um_d/(m_u+m_d)$ に一致する。クォーク質量が等しく $\theta\simeq\pi$ に近い場合は根が非常に小さくなるため、$\theta=0$ の曲率だけで全周期を表すことはできない。特に高次の有効項や複数の枝を考える必要がある。アクシオンの小さい振動質量は最小での曲率で足りるが、場が最小から大きく離れた初期条件の時間発展には角度全体のポテンシャルが必要である。

### アクシオンの小振動と混合の低エネルギー展開

正の二クォーク質量の先頭の有効理論で、有効角を $\vartheta=\bar\theta+a/f_a$ とする。純ゲージの重い位相を消去して $\phi_u+\phi_d=\vartheta$ と近似できる場合、小さい位相の真空エネルギーは $\Sigma(m_u\phi_u^2+m_d\phi_d^2)/2$ である。$\phi_d=\vartheta-\phi_u$ を入れて微分すると $\phi_u=m_d\vartheta/(m_u+m_d)$、$\phi_d=m_u\vartheta/(m_u+m_d)$ となる。従って曲率は $\chi=\Sigma\,m_um_d/(m_u+m_d)$ であり、正準運動項 $\tfrac12(\partial a)^2$ の場の質量は $m_a^2=\chi/f_a^2$ である。$F^2m_\pi^2=(m_u+m_d)\Sigma$ を代入すると
$$
m_a=\frac{F m_\pi}{f_a}\frac{\sqrt z}{1+z},\qquad z=m_u/m_d .
$$
全ての質量と減衰定数は同じ自然単位系で代入する。この式は軽いアクシオンで $F/f_a$ が小さく、π 中間子とアクシオンの運動項の混合が高次である場合の先頭である。

:::example[大きい対称性の尺度と小さい質量]
$F=0.092\,\mathrm{GeV}$、$m_\pi=0.135\,\mathrm{GeV}$、$z=1/2$、$f_a=10^{12}\,\mathrm{GeV}$ なら $m_a\simeq5.85\times10^{-15}\,\mathrm{GeV}=5.85\,\mathrm{\mu eV}$ である。$f_a$ を十倍に上げると質量は十分の一になる。一方、曲率を支えるハドロンの尺度は変えない。
:::

:::example[一つの軽い質量を零へ近づける]
$m_u/m_d=z\to0^+$ とすると $m_a\propto\sqrt z$ となる。この極限では $\phi_u\simeq\vartheta$ が角をほぼすべて吸収し、質量によるエネルギーの代価が小さくなる。$\theta$ の効果が消える極限と、両方の物理的な軽い質量が有限である場合のアクシオンによる緩和を区別できる。
:::

### 計算例

:::example[最小からの場の変位]
$f_a$ を GeV、位相を $\theta$ として、同じ枝の最小の場の値 $a$ を求める。入力値を次のようにとる。

$$
f_{a}=1000,\qquad \theta=\frac{1}{10}
$$

与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
a&=- f_{a} \theta=-100\quad[\mathrm{GeV}]
\end{aligned}
$$

角の周期を考慮する場合は等価な枝を選べる。
:::

:::example[感受率と質量]
感受率の四乗根 $X$ と $f$ を GeV で与え、最低次のアクシオン質量 $m$ を求める。入力値を次のようにとる。

$$
X=\frac{3}{40},\qquad f=1000
$$

与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
m&=\frac{X^{2}}{f}=\frac{9}{1600000}\simeq 5.625\times10^{-6}\quad[\mathrm{GeV}]
\end{aligned}
$$

これは QCD が与える真空曲率だけを残した値である。
:::

:::example[感受率の低エネルギー模型]
$\Sigma$ を $\mathrm{GeV}^3$、$m_1,m_2$ を GeV で与え、$\chi$ を求める。入力値を次のようにとる。

$$
\Sigma=\frac{1}{50},\qquad m_{1}=\frac{1}{500},\qquad m_{2}=\frac{1}{200}
$$

与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
\chi&=\frac{\Sigma m_{1} m_{2}}{m_{1} + m_{2}}=\frac{1}{35000}\simeq 2.85714\times10^{-5}\quad[\mathrm{GeV^4}]
\end{aligned}
$$

この感受率には軽い動的クォークを含めている。
:::

:::example[有限位相でのエネルギー差]
$m_u,m_d$ を GeV、$\Sigma$ を $\mathrm{GeV}^3$ として、$\theta=\pi$ と零の真空エネルギー密度の差 $D$ を求める。入力値を次のようにとる。

$$
m_{u}=\frac{1}{500},\qquad m_{d}=\frac{1}{200},\qquad \Sigma=\frac{1}{50}
$$

与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
D&=\Sigma \left(m_{d} + m_{u} - \left|{m_{d} - m_{u}}\right|\right)=\frac{1}{12500}\simeq 8\times10^{-5}\quad[\mathrm{GeV^4}]
\end{aligned}
$$

最低次の二種の質量項だけを使った。
:::

## 演習問題

::::exercise[寿命を比較する]{level=basic}
二固有状態の幅 $G_S,G_L$ を同じ単位で与え、寿命の比 $R=\tau_L/\tau_S$ を求める。入力値を次のようにとる。

$$
G_{S}=500,\qquad G_{L}=1
$$

:::solution
与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
R&=\frac{G_{S}}{G_{L}}=500
\end{aligned}
$$

寿命と幅は逆数の関係にある。
:::
::::

::::exercise[純粋に分散的な混合]{level=advanced}
$\Gamma_{12}=0$ とし、$M_{12}$ の大きさ $m$ を MeV で与え、質量差 $D$ を求める。入力値を次のようにとる。

$$
m=\frac{1}{500}
$$

:::solution
与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
D&=2 m=\frac{1}{250}\simeq 0.004\quad[\mathrm{MeV}]
\end{aligned}
$$

K 系で吸収的な混合を実際に落とせるという仮定ではない。
:::
::::

::::exercise[固有状態の重なり]{level=advanced}
規格化された $p$ の二乗の大きさ $P$ を与え、CPT を仮定した二固有状態の重なり $D$ を求める。入力値を次のようにとる。

$$
P=\frac{51}{100}
$$

:::solution
与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
D&=2 P - 1=\frac{1}{50}
\end{aligned}
$$

共通の位相規約で $q$ の大きさの二乗に $1-P$ を使った。
:::
::::

::::exercise{level=advanced}
吸収行列について $\abs{\Gamma_{12}}^2\le\Gamma_{11}\Gamma_{22}$ を示せ。

:::solution
各初期状態から終状態への振幅を、エネルギー殻上の測度 $2\pi\delta(E-E_f)$ で積分する内積空間のベクトルとみなす。$\Gamma_{12}$ は二ベクトルの内積、対角成分はノルム二乗なので Cauchy–Schwarz から従う。同じ結論は二次元の正半定値行列の行列式が非負であることでも得られる。
:::
::::

::::exercise[二重比から直接の破れへ]{level=basic}
幅の二重比 $R$ から一次近似の $\Re(\epsilon'/\epsilon)=r$ を求める。入力値を次のようにとる。

$$
R=\frac{247}{250}
$$

:::solution
与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
r&=\frac{1}{6} - \frac{R}{6}=\frac{1}{500}\simeq 0.002
\end{aligned}
$$

幅の比を用いるため、振幅の比をさらに二乗した結果である。
:::
::::

::::exercise[二重比の二次効果]{level=advanced}
$z=x+iy$ を指定し、二重比 $R$ を展開せずに求める。入力値を次のようにとる。

$$
x=0,\qquad y=\frac{1}{100}
$$

:::solution
与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
R&=\frac{4 y^{2} + \left(1 - 2 x\right)^{2}}{y^{2} + \left(x + 1\right)^{2}}=\frac{10004}{10001}\simeq 1.0003
\end{aligned}
$$

一次の実部が零になる場合でも二次の変化は残る。
:::
::::

::::exercise[半レプトンの規格化差]{level=advanced}
$\epsilon_m=a+ib$ を与え、共役の直接振幅の大きさが同じと仮定して $\delta_L=d$ を求める。入力値を次のようにとる。

$$
a=\frac{3}{1000},\qquad b=\frac{1}{250}
$$

:::solution
与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
d&=\frac{2 a}{a^{2} + b^{2} + 1}=\frac{240}{40001}\simeq 0.00599985
\end{aligned}
$$

混合の位相の全てをこの一つの率から決めることはできない。
:::
::::

::::exercise[時間積分した混合]{level=basic}
大きさの比 $x=\Delta m/\Gamma$ を与え、幅差を無視した $\chi$ を得る。入力値を次のようにとる。

$$
x=2
$$

:::solution
与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
\chi&=\frac{x^{2}}{2 x^{2} + 2}=\frac{2}{5}
\end{aligned}
$$

非常に速い振動を時刻を区別せず測ると混合率は半分に近い。
:::
::::

::::exercise[時間分解した混合]{level=advanced}
位相 $\phi=\Delta mt$ と生存因子 $S=e^{-\Gamma t}$ を与え、混合確率 $P$ を求める。入力値を次のようにとる。

$$
\phi=\frac{\pi}{2},\qquad S=\frac{4}{5}
$$

:::solution
与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
P&=\frac{S \left(1 - \cos{\left(\phi \right)}\right)}{2}=\frac{2}{5}
\end{aligned}
$$

生存した事象内での条件付き確率なら $S$ で割る。
:::
::::

::::exercise[二方向の混合の比]{level=advanced}
混合係数比の大きさ $r=\abs{q/p}$ を与え、$B\to\bar B$ と逆方向の率比 $R$ を求める。入力値を次のようにとる。

$$
r=\frac{51}{50}
$$

:::solution
与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
R&=r^{4}=\frac{6765201}{6250000}\simeq 1.08243
\end{aligned}
$$

直接の半レプトン振幅の絶対値が等しい近似を置いた。
:::
::::

::::exercise[タグの希釈]{level=advanced}
タグ誤り確率 $w$ と真の振動コントラスト $A$ から、観測する振幅の係数 $D$ を求める。入力値を次のようにとる。

$$
w=\frac{1}{10},\qquad A=\frac{4}{5}
$$

:::solution
与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
D&=A \left(1 - 2 w\right)=\frac{16}{25}
\end{aligned}
$$

時刻の分解能と背景の効果はここでは落とした。
:::
::::

::::exercise[混合だけの半レプトン非対称]{level=advanced}
$r=\abs{q/p}$ を与え、正の向きを本文のように取った $a$ を求める。入力値を次のようにとる。

$$
r=\frac{201}{200}
$$

:::solution
与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
a&=\frac{1 - r^{4}}{r^{4} + 1}=- \frac{32240801}{3232240801}\simeq -0.00997475
\end{aligned}
$$

$\Im z$ の展開を使わず二方向の率を比べた値である。
:::
::::

::::exercise[干渉係数の計算]{level=basic}
$\lambda_f=a+ib$ として $S$ を求める。入力値を次のようにとる。

$$
a=0,\qquad b=\frac{1}{2}
$$

:::solution
与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
S&=\frac{2 b}{a^{2} + b^{2} + 1}=\frac{4}{5}
\end{aligned}
$$

$\lambda_f$ の符号規約は固有状態の表記と合わせる。
:::
::::

::::exercise[直接非対称性の係数]{level=advanced}
$\lambda_f$ の絶対値を $r$ とし、余弦項の係数 $C$ を求める。入力値を次のようにとる。

$$
r=\frac{1}{2}
$$

:::solution
与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
C&=\frac{1 - r^{2}}{r^{2} + 1}=\frac{3}{5}
\end{aligned}
$$

この式では混合の絶対値を 1 と仮定している。
:::
::::

::::exercise[実部の規格化]{level=advanced}
$\lambda_f=a+ib$ の実部と虚部を与え、$\sinh$ の規格化係数 $D$ を求める。入力値を次のようにとる。

$$
a=\frac{4}{5},\qquad b=\frac{3}{5}
$$

:::solution
与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
D&=- \frac{2 a}{a^{2} + b^{2} + 1}=- \frac{4}{5}
\end{aligned}
$$

幅差の符号を交換する規約では係数の符号を合わせる。
:::
::::

::::exercise[二経路の最大位相差の例]{level=advanced}
強位相差と弱位相差の両方を $\pi/2$ とし、強度比 $r$ から率の非対称性 $A$ を計算する。入力値を次のようにとる。

$$
r=\frac{1}{5}
$$

:::solution
与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
A&=\frac{2 r}{r^{2} + 1}=\frac{5}{13}
\end{aligned}
$$

分母の干渉が消えるが、分子の干渉は残る。
:::
::::

::::exercise[時間を積分した非対称]{level=advanced}
$x=\Delta m/\Gamma$ と係数 $S,C$ を与え、一様な時間の効率での非対称 $A$ を求める。入力値を次のようにとる。

$$
x=2,\qquad S=\frac{3}{5},\qquad C=\frac{1}{10}
$$

:::solution
与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
A&=\frac{- C + S x}{x^{2} + 1}=\frac{11}{50}
\end{aligned}
$$

幅差は零で $\abs{q/p}=1$ とした。
:::
::::

::::exercise{level=advanced}
$\abs{q/p}=1$ のとき $S_f^2+C_f^2\le1$ を示し、等号条件を求めよ。

:::solution
$\lambda_f=a+ib$ と置く。共通分母 $(1+a^2+b^2)^2$ から分子 $4b^2+(1-a^2-b^2)^2$ を引くと $4a^2$ になる。分母は正であるから不等式が成り立ち、等号は $\Re\lambda_f=0$ のときである。
:::
::::

::::exercise[重い質量展開の小さい量]{level=basic}
残差尺度 $L$ と質量 $m$ を GeV で与えて $r=L/m$ を求める。入力値を次のようにとる。

$$
L=\frac{2}{5},\qquad m=\frac{8}{5}
$$

:::solution
与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
r&=\frac{L}{m}=\frac{1}{4}
\end{aligned}
$$

質量が小さい系では高次の補正が相対的に大きい。
:::
::::

::::exercise[超微細分裂を比較する]{level=advanced}
同じ軽い場を仮定し、重い質量 $m_1,m_2$ から第 2 対第 1 の分裂比 $r$ を最低次で求める。入力値を次のようにとる。

$$
m_{1}=2,\qquad m_{2}=5
$$

:::solution
与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
r&=\frac{m_{1}}{m_{2}}=\frac{2}{5}
\end{aligned}
$$

短距離係数の走りはこの単純な比に入れていない。
:::
::::

::::exercise[速度積を計算する]{level=advanced}
親と娘の質量 $M,m$ を GeV、$q_2=q^2$ を $\mathrm{GeV}^2$ で与え、$w$ を求める。入力値を次のようにとる。

$$
M=5,\qquad m=2,\qquad q_{2}=9
$$

:::solution
与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
w&=\frac{M^{2} + m^{2} - q_{2}}{2 M m}=1
\end{aligned}
$$

$q_2$ の上端で速度積は 1 になる。
:::
::::

::::exercise[軽い角運動量と分裂]{level=advanced}
$j_\ell=j>0$ とモデルの一対の係数 $C$ を MeV で与え、二重項の質量差 $D$ を求める。入力値を次のようにとる。

$$
j=\frac{3}{2},\qquad C=50
$$

:::solution
与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
D&=\frac{C \left(2 j + 1\right)}{2}=100\quad[\mathrm{MeV}]
\end{aligned}
$$

$C$ が重い質量の逆数に比例する範囲に限って比較する。
:::
::::

::::exercise[最小からの場の変位]{level=basic}
$f_a$ を GeV、位相を $\theta$ として、同じ枝の最小の場の値 $a$ を求める。入力値を次のようにとる。

$$
f_{a}=2000,\qquad \theta=- \frac{1}{5}
$$

:::solution
与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
a&=- f_{a} \theta=400\quad[\mathrm{GeV}]
\end{aligned}
$$

角の周期を考慮する場合は等価な枝を選べる。
:::
::::

::::exercise[感受率の低エネルギー模型]{level=advanced}
$\Sigma$ を $\mathrm{GeV}^3$、$m_1,m_2$ を GeV で与え、$\chi$ を求める。入力値を次のようにとる。

$$
\Sigma=\frac{1}{40},\qquad m_{1}=\frac{3}{1000},\qquad m_{2}=\frac{3}{500}
$$

:::solution
与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
\chi&=\frac{\Sigma m_{1} m_{2}}{m_{1} + m_{2}}=\frac{1}{20000}\simeq 5\times10^{-5}\quad[\mathrm{GeV^4}]
\end{aligned}
$$

この感受率には軽い動的クォークを含めている。
:::
::::

::::exercise[有限位相でのエネルギー差]{level=advanced}
$m_u,m_d$ を GeV、$\Sigma$ を $\mathrm{GeV}^3$ として、$\theta=\pi$ と零の真空エネルギー密度の差 $D$ を求める。入力値を次のようにとる。

$$
m_{u}=\frac{3}{1000},\qquad m_{d}=\frac{3}{500},\qquad \Sigma=\frac{1}{40}
$$

:::solution
与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
D&=\Sigma \left(m_{d} + m_{u} - \left|{m_{d} - m_{u}}\right|\right)=\frac{3}{20000}\simeq 0.00015\quad[\mathrm{GeV^4}]
\end{aligned}
$$

最低次の二種の質量項だけを使った。
:::
::::

::::exercise[希釈を含む混合信号の振動数を逆算する]{level=advanced}
幅差が零で $|q/p|=1$ の模型で、タグの誤り確率を $w=1/10$ とする。タグと崩壊のフレーバーが違う事象の時間積分割合が $3/10$ であった。タグ誤りが真の混合と独立なとき、$\Delta m/\Gamma$ の絶対値を求めよ。

:::solution
真の混合割合を $\chi$ とすると、タグが正しく混合が起きた場合と、タグが誤り混合が起きなかった場合を足して $\chi_{\rm obs}=(1-w)\chi+w(1-\chi)=w+(1-2w)\chi$ である。従って $\chi=(3/10-1/10)/(1-2/10)=1/4$ となる。式 :eqref[eq-sm-cp-integrated-mixing] に入れれば $x^2=2\chi/(1-2\chi)=1$、よって $|\Delta m/\Gamma|=1$ である。この積分では $\Delta m$ の符号が余弦の中で消えるため、符号を決めるには別の位相情報が要る。
:::
::::

## 確認テスト

::::quiz
同じ終状態に混合の有無の二経路で到達できるとき、どの量が両経路を比較するか。

- [x] $\lambda_f=(q/p)(\bar A_f/A_f)$。
- [ ] $A_f+\bar A_f$ の位相だけ。
- [ ] $q$ の絶対的な位相。

:::explanation
初めの粒子から直接崩壊する振幅は $A_f$ に、反粒子へ混合して崩壊する振幅は $(q/p)\bar A_f$ に比例する。比を取ると二経路の相対位相となり、フレーバー状態の再位相化で不変である。
:::
::::

::::quiz
二つの崩壊振幅で弱位相差だけが非零のとき、直接の率の CP 非対称性はどうなるか。

- [ ] 最大となる。
- [x] 強位相差が零なら率の差は零である。
- [ ] 混合がなくても必ず半分になる。

:::explanation
二経路の干渉で粒子と反粒子には $\cos(\Delta\delta\pm\Delta\phi)$ が現れる。差は $\sin\Delta\delta\sin\Delta\phi$ に比例するので強位相差が必要である。
:::
::::

::::quiz{answer="9e-06" tol="9e-11"}
$\epsilon'$ の実部を $e$ とし、$\Re(\eta_{+-}-\eta_{00})=d$ を求める。入力値を次のようにとる。

$$
e=\frac{3}{1000000}
$$

最後に求める量の数値を答えよ。

:::explanation
与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
d&=3 e=\frac{9}{1000000}\simeq 9\times10^{-6}
\end{aligned}
$$

共通の混合寄与は差を取ると消える。
:::
::::

::::quiz{answer="6.4e-07" tol="6.4e-12"}
感受率の四乗根 $X$ と $f$ を GeV で与え、最低次のアクシオン質量 $m$ を求める。入力値を次のようにとる。

$$
X=\frac{2}{25},\qquad f=10000
$$

最後に求める量の数値を答えよ。

:::explanation
与えた値を代入し、中間量から順に求める。

$$
\begin{aligned}
m&=\frac{X^{2}}{f}=\frac{1}{1562500}\simeq 6.4\times10^{-7}\quad[\mathrm{GeV}]
\end{aligned}
$$

これは QCD が与える真空曲率だけを残した値である。
:::
::::

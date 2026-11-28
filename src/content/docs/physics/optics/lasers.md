---
title: レーザーの原理
description: 熱平衡と遷移率から誘導放出を導き、反転分布、利得飽和、発振とパルス化を速度方程式で解析する。
notation:
  - symbol: 'A_{21},B_{21},B_{12}'
    meaning: Einstein の遷移係数
  - symbol: 'N_1,N_2'
    meaning: レーザー遷移の下準位と上準位の数密度
  - symbol: '\sigma(\nu)'
    meaning: 遷移の断面積
  - symbol: 'g'
    meaning: パワーの利得係数
  - symbol: '\tau_{\mathrm{ph}}'
    meaning: 共振器の光子寿命
---

## この章の目標

- 熱平衡から Einstein 係数の関係を導く。
- 線形状と反転分布を用いて増幅と飽和の条件を求める。
- ポンプと共振器損失から発振しきい値を計算する。
- Q スイッチとモード同期のパルス形成を区別する。

## 黒体放射と Einstein の A 係数・B 係数

### 光の状態密度と熱平衡
真空の一辺 $L$ の周期箱で、波数格子の間隔は $2\pi/L$ である。周波数 $\nu$ に対応する球殻を二偏光で数えると、体積および周波数当たりのモード数は $2(4\pi k^2)/(2\pi)^3\,\dd k/\dd\nu=8\pi\nu^2/c^3$。各モードに光子が $m\ge0$ 個入る Boltzmann の重みは $e^{-mh\nu/(k_BT)}$ に比例する。幾何級数を規格化し、その比で微分して平均個数 $1/(e^{h\nu/(k_BT)}-1)$ を得る。これに光子エネルギーを掛けた黒体の分光エネルギー密度は
$$
\rho_\nu(T)=\frac{8\pi h\nu^3/c^3}{e^{h\nu/(k_BT)}-1}.
$$
ここでは零点エネルギーを熱放射に含めない。光子の占有数は :ref[def-stm-photons-phonons] でも扱う。

二準位間隔を $h\nu$、縮退度を $d_1,d_2$ とする。黒体場が線幅内でほぼ一定なら、下準位の一粒子は $\rho_\nu B_{12}$ の率で吸収し、上準位の一粒子は $\rho_\nu B_{21}$ の率で誘導放出、$A_{21}$ の率で自発放出する。これが :term[Einstein 係数]{en="Einstein coefficients" yomi="あいんしゅたいんけいすう"} の定義である。$A_{21}$ は $\mathrm{s^{-1}}$、各 $B$ は $\rho_\nu$ と掛けて $\mathrm{s^{-1}}$ になる単位を持つ。熱平衡で $N_2/N_1=(d_2/d_1)e^{-h\nu/(k_BT)}$、上下の流量が等しいので
$$
\rho_\nu=\frac{A_{21}}{(d_1/d_2)e^{h\nu/(k_BT)}B_{12}-B_{21}}.
$$
全ての $T$ で黒体の式と等しいことから
$$
d_1B_{12}=d_2B_{21},\qquad A_{21}=\frac{8\pi h\nu^3}{c^3}B_{21}.
\label{eq-optics-einstein}
$$
真空モードと電気双極子遷移の微視的導出は :eqref[eq-qm-spontaneous] にある。

:::example[縮退を取り違えない]
$d_2=3,d_1=1$ なら $B_{12}=3B_{21}$ である。高温で $N_2/N_1\to3$ となっても、吸収と誘導放出が同じになるだけで正味の増幅は生じない。
:::
:::example[古典極限]
$h\nu\ll k_BT$ では分母は $h\nu/(k_BT)$ に近づき、$\rho_\nu\simeq(8\pi\nu^2/c^3)k_BT$。一モード当たりの熱エネルギーは $k_BT$ になる。
:::

### 真空と熱光で同じ遷移を数える
式 :eqref[eq-optics-einstein] は $\rho_\nu$ が周波数 Hz 当たりである約束に依存する。角周波数で光を数えるなら $\rho_\omega=\rho_\nu/(2\pi)$ なので、同じ物理的な誘導率を得る $B^{(\omega)}=2\pi B^{(\nu)}$ である。不変なのは $B\rho$ であり、係数そのものの比較には密度の規約が必要である。真空中の光の占有数 $\bar m=(e^{h\nu/(k_BT)}-1)^{-1}$ を使うと、平衡放射での誘導放出率は $A_{21}\bar m$。自発を加えた率は $A_{21}(1+\bar m)$ となる。

### 双極子の断面積に変換する
一方向の平面波の照度は真空で $I=c\rho$ である。ここでは $B$ の定義が、偏光と方向を平均した黒体場に対するものであることに注意する。固定の双極子と固定の入射偏光の遷移は $|\vb d\cdot\vb e|^2$ に比例する（[時間に依存する摂動論の電気双極子遷移](/physics/quantum-mechanics/time-dependent-perturbation/)）。双極子方向を等方に平均して、入射偏光の選択に依存しない粒子集合を考える。Hz 当たりの面積 1 の線形状 $L(\nu)$ を使い、狭帯域の一方向照度 $I$ での放出率を $B_{21}(I/c)L(\nu)$ とすると、断面積の定義から $\sigma_e=h\nu B_{21}L(\nu)/c$。式 :eqref[eq-optics-einstein] を入れれば
$$
\sigma_e(\nu)=\frac{\lambda_0^2}{8\pi}A_{21}L(\nu)
$$
となる。狭い線なので前置因子の波長や周波数は中心の値とする。線形状は $\mathrm{Hz^{-1}}$、$A$ は $\mathrm{s^{-1}}$ なので断面積は $\mathrm{m^2}$。固定の双極子では偏光への射影を用い、誘電体では媒質中のモード密度とエネルギー密度に合わせて換算する。

:::example[線幅が同じでない媒質]
同じ中心周波数と $A_{21}$ を持つ Lorentz 線の全半値幅を二倍にすると、積分した断面積は同じでピークが半分になる。狭帯域光のしきい値に効くのは、発振する周波数での断面積である。
:::

### 熱放射のモード占有と自然放出の比

温度 $T$ で周波数 $\nu$ の一モードの光子数は $m=0,1,\ldots$ で、熱平衡の重みが $q^m$、$q=e^{-h\nu/(k_BT)}$ になる。等比和を規格化し、$q\,\dd/\dd q$ を作用させれば平均占有数は $q/(1-q)$。真空中で二偏光を含む単位体積・単位周波数当たりのモード数 $8\pi\nu^2/c^3$ を掛ければ Planck 分布を得る。狭い線がこの熱放射に接するなら、:eqref[eq-optics-einstein] により一個の励起原子に対する誘導放出率と自然放出率の比が平均占有数に等しい。熱平衡の基底状態についても同時に吸収が起こるため、誘導放出が大きい帯域を選んだだけでは光は増幅しない。

:::example[室温のマイクロ波]
低い光子エネルギーで誘導と自然放出の比を確認する。

入力は遷移の周波数 $\nu\simeq 1\times 10^{10}\,\mathrm{Hz}$、放射の温度 $T=300\,\mathrm K$、自然放出率 $A_{21}=1\,\mathrm{s^{-1}}$ である。

光子エネルギーを $k_BT$ で割る。

$$
\begin{aligned}
x&=\frac{h\nu}{k_B T},\\
x&\simeq \frac{132521403 \cdot \left(1\times 10^{10}\right)}{2761298000000000000 \cdot \left(300\right)}\simeq 0.0015997.
\end{aligned}
$$

一個を追加するための熱平衡の重みを求める。

$$
\begin{aligned}
q&=e^{- x},\\
q&\simeq e^{- \left(0.001599747691\right)}\simeq 0.9984.
\end{aligned}
$$

可能な占有数を等比和で平均する。

$$
\begin{aligned}
\bar n&=\frac{q}{1 - q},\\
\bar n&\simeq \frac{\left(0.9984015312\right)}{1 - \left(0.9984015312\right)}\simeq 624.6.
\end{aligned}
$$

Einstein 係数の比を用いて熱場による放出を求める。

$$
\begin{aligned}
W_{\rm stim}&=A_{21} \bar n,\\
W_{\rm stim}/(\mathrm{s^{-1}})&\simeq \left(1\right) \cdot \left(624.598707\right)\simeq 624.6.
\end{aligned}
$$

自然放出と誘導放出の二過程を足す。

$$
\begin{aligned}
f&=\frac{W_{\rm stim}}{A_{21} + W_{\rm stim}},\\
f&\simeq \frac{\left(624.598707\right)}{\left(1\right) + \left(624.598707\right)}\simeq 0.9984.
\end{aligned}
$$

熱占有数は一より大きい。
:::

:::example[室温の可視光]
温度を変えず遷移エネルギーを光学周波数へ上げる。

入力は遷移の周波数 $\nu\simeq 5\times 10^{14}\,\mathrm{Hz}$、放射の温度 $T=300\,\mathrm K$、自然放出率 $A_{21}=100000000\,\mathrm{s^{-1}}$ である。

光子エネルギーを $k_BT$ で割る。

$$
\begin{aligned}
x&=\frac{h\nu}{k_B T},\\
x&\simeq \frac{132521403 \cdot \left(5\times 10^{14}\right)}{2761298000000000000 \cdot \left(300\right)}\simeq 79.987.
\end{aligned}
$$

一個を追加するための熱平衡の重みを求める。

$$
\begin{aligned}
q&=e^{- x},\\
q&\simeq e^{- \left(79.98738456\right)}\simeq 1.8278\times 10^{-35}.
\end{aligned}
$$

可能な占有数を等比和で平均する。

$$
\begin{aligned}
\bar n&=\frac{q}{1 - q},\\
\bar n&\simeq \frac{\left(1.827764616\times 10^{-35}\right)}{1 - \left(1.827764616\times 10^{-35}\right)}\simeq 1.8278\times 10^{-35}.
\end{aligned}
$$

Einstein 係数の比を用いて熱場による放出を求める。

$$
\begin{aligned}
W_{\rm stim}&=A_{21} \bar n,\\
W_{\rm stim}/(\mathrm{s^{-1}})&\simeq \left(100000000\right) \cdot \left(1.827764616\times 10^{-35}\right)\simeq 1.8278\times 10^{-27}.
\end{aligned}
$$

自然放出と誘導放出の二過程を足す。

$$
\begin{aligned}
f&=\frac{W_{\rm stim}}{A_{21} + W_{\rm stim}},\\
f&\simeq \frac{\left(1.827764616\times 10^{-27}\right)}{\left(100000000\right) + \left(1.827764616\times 10^{-27}\right)}\simeq 1.8278\times 10^{-35}.
\end{aligned}
$$

室温の熱放射による可視域の誘導放出は抑えられる。
:::

:::example[高温の放射場]
可視光の周波数を保ち、温度を上げる。

入力は遷移の周波数 $\nu\simeq 5\times 10^{14}\,\mathrm{Hz}$、放射の温度 $T=6000\,\mathrm K$、自然放出率 $A_{21}=100000000\,\mathrm{s^{-1}}$ である。

光子エネルギーを $k_BT$ で割る。

$$
\begin{aligned}
x&=\frac{h\nu}{k_B T},\\
x&\simeq \frac{132521403 \cdot \left(5\times 10^{14}\right)}{2761298000000000000 \cdot \left(6000\right)}\simeq 3.9994.
\end{aligned}
$$

一個を追加するための熱平衡の重みを求める。

$$
\begin{aligned}
q&=e^{- x},\\
q&\simeq e^{- \left(3.999369228\right)}\simeq 0.018327.
\end{aligned}
$$

可能な占有数を等比和で平均する。

$$
\begin{aligned}
\bar n&=\frac{q}{1 - q},\\
\bar n&\simeq \frac{\left(0.01832719553\right)}{1 - \left(0.01832719553\right)}\simeq 0.018669.
\end{aligned}
$$

Einstein 係数の比を用いて熱場による放出を求める。

$$
\begin{aligned}
W_{\rm stim}&=A_{21} \bar n,\\
W_{\rm stim}/(\mathrm{s^{-1}})&\simeq \left(100000000\right) \cdot \left(0.0186693524\right)\simeq 1.8669\times 10^{6}.
\end{aligned}
$$

自然放出と誘導放出の二過程を足す。

$$
\begin{aligned}
f&=\frac{W_{\rm stim}}{A_{21} + W_{\rm stim}},\\
f&\simeq \frac{\left(1.86693524\times 10^{6}\right)}{\left(100000000\right) + \left(1.86693524\times 10^{6}\right)}\simeq 0.018327.
\end{aligned}
$$

誘導の比は室温より増すが、放出だけから正味利得は判定できない。
:::

## 誘導放出と反転分布

### 増幅に必要な占有数
:term[誘導放出]{en="stimulated emission" yomi="ゆうどうほうしゅつ"} では入射場と同じモードの占有数が増える。量子化した場で光子数 $m$ の吸収行列要素には $\sqrt m$、放出には $\sqrt{m+1}$ が掛かる。これは消滅・生成演算子を数状態に作用させた係数である。二乗すると放出の $m+1$ は誘導の $m$ と真空からの 1 に分かれる。この因子は占有数の変化を表す。入射が数状態なら絶対位相は定まらない。入射が位相を持つ場であれば、媒質の誘導分極がその場と一定の位相関係を持つため、伝播方向と偏光を保った増幅が生じる。遷移率を行列要素の二乗から作る一般式は :ref[law-qm-golden] を参照する。

誘導放出率から吸収率を引いた差は $B_{21}[N_2-(d_2/d_1)N_1]$ に比例する。角括弧が正になることを :term[反転分布]{en="population inversion" yomi="はんてんぶんぷ"} と呼ぶ。正温度の平衡では $N_2/d_2<N_1/d_1$ なので、選んだ遷移で反転を維持するには外部からのポンプが要る。閉じた二準位系を同じ吸収・放出率の共鳴光だけで駆動する等縮退の定常モデルでは、上の粒子数は全粒子の半分を超えない。実際、 $N_2'=W(N-2N_2)-N_2/T_1$ なので定常解は $N_2/N=W/(2W+1/T_1)<1/2$ である。

:::example[上の粒子数だけでは決まらない]
$d_2/d_1=3$、$N_2=2N_1$ では差は $-N_1$ であり、全粒子数では上準位が多くても正味は吸収である。
:::
:::example[強い共鳴光の極限]
等縮退閉二準位で $WT_1=100$ なら上準位の割合は $100/201\simeq0.498$ となる。強くして吸収を飽和させても定常な反転にはならない。
:::

### 過渡反転と定常反転
二準位でも短い共鳴パルスで過渡的に上準位を占有させることは可能である。緩和を無視した共鳴相互作用を回転座標で $H=(\hbar\Omega_R/2)\begin{pmatrix}0&1\\1&0\end{pmatrix}$ とすると、初期の下準位から $N_2/N=\sin^2(\Omega_R t/2)$ を得る。これは行列の二乗が単位行列の定数倍になることから指数級数の偶奇を分けた解である。パルス面積が $\pi$ になると全粒子が上に移る。先の定常速度模型は位相を捨てて上下の率で記述するため、この短時間の振動は表せない。共鳴を保ち駆動位相が一定なら、時間ごとの Hamiltonian は同じ行列の定数倍で可換なので、$\Omega_R t$ をパルス面積 $\int\Omega_R(t)\dd t$ に置き換えられる。

:::example[弱い反転を損失と比べる]
等縮退で $N_2=0.55N$、$N_1=0.45N$ なら反転差は $0.10N$。$\sigma=10^{-20}\,\mathrm{m^2}$、$N=10^{22}\,\mathrm{m^{-3}}$ では利得は $10\,\mathrm{m^{-1}}$ である。分布損失が $12\,\mathrm{m^{-1}}$ なら反転しているが進行方向のパワーは減る。
:::

### 縮退のある反転を占有分率で判定する

全粒子数密度 $N=N_1+N_2$、上準位分率 $f=N_2/N$ とする。同一線形状で Einstein 関係 $B_{12}=(d_2/d_1)B_{21}$ を使うと、誘導での正味の光子数の増加率は $N B_{21}[(1+d_2/d_1)f-d_2/d_1]$ に比例する。従って光が強くなる条件は $f>d_2/(d_1+d_2)$。占有数を縮退の数で割った値を比べた結果である。熱平衡の $f/(1-f)=(d_2/d_1)e^{-h\nu/(k_BT)}$ は正温度でこの条件を満たさない。光によって二準位を対称に励起し、強度だけを増すと、吸収と誘導放出が釣り合う占有分率に近づく。

:::example[縮退のない反転]
上準位が全体の六割で両準位の縮退が等しい。

入力は下準位の縮退度 $d_{1}=1$、上準位の縮退度 $d_{2}=1$、上準位の占有分率 $f=\frac{3}{5}$、総粒子密度 $N\simeq 1\times 10^{22}\,\mathrm{m^{-3}}$、ピークでの誘導放出断面積 $\sigma\simeq 1\times 10^{-20}\,\mathrm{m^2}$ である。

同じ線形状の Einstein 係数の関係を使う。

$$
\begin{aligned}
r&=\frac{d_{2}}{d_{1}},\\
r&= \frac{\left(1\right)}{\left(1\right)}=1.
\end{aligned}
$$

上と下の占有を縮退一個当たりで等しくする。

$$
\begin{aligned}
f_{\rm tr}&=\frac{d_{2}}{d_{1} + d_{2}},\\
f_{\rm tr}&= \frac{\left(1\right)}{\left(1\right) + \left(1\right)}=\frac{1}{2}.
\end{aligned}
$$

全体の密度を入射時の分率で配分する。

$$
\begin{aligned}
N_2&=N f,\\
N_2/(\mathrm{m^{-3}})&\simeq \left(1\times 10^{22}\right) \cdot \left(\frac{3}{5}\right)\simeq 6\times 10^{21}.
\end{aligned}
$$

下準位からの吸収を断面積比で重み付けして引く。

$$
\begin{aligned}
\Delta N&=N_2 - r \left(N - N_2\right),\\
\Delta N/(\mathrm{m^{-3}})&\simeq \left(6\times 10^{21}\right) - \left(1\right) \cdot \left(\left(1\times 10^{22}\right) - \left(6\times 10^{21}\right)\right)\simeq 2\times 10^{21}.
\end{aligned}
$$

同じ線での誘導放出断面積を掛ける。

$$
\begin{aligned}
g&=\Delta N \sigma,\\
g/(\mathrm{m^{-1}})&\simeq \left(2\times 10^{21}\right) \cdot \left(1\times 10^{-20}\right)=20.
\end{aligned}
$$

正味の誘導効果は増幅である。
:::

:::example[上準位の縮退を増す]
同じ粒子分率を上に置くが、上の縮退を二倍にする。

入力は下準位の縮退度 $d_{1}=1$、上準位の縮退度 $d_{2}=2$、上準位の占有分率 $f=\frac{3}{5}$、総粒子密度 $N\simeq 1\times 10^{22}\,\mathrm{m^{-3}}$、ピークでの誘導放出断面積 $\sigma\simeq 1\times 10^{-20}\,\mathrm{m^2}$ である。

同じ線形状の Einstein 係数の関係を使う。

$$
\begin{aligned}
r&=\frac{d_{2}}{d_{1}},\\
r&= \frac{\left(2\right)}{\left(1\right)}=2.
\end{aligned}
$$

上と下の占有を縮退一個当たりで等しくする。

$$
\begin{aligned}
f_{\rm tr}&=\frac{d_{2}}{d_{1} + d_{2}},\\
f_{\rm tr}&= \frac{\left(2\right)}{\left(1\right) + \left(2\right)}=\frac{2}{3}.
\end{aligned}
$$

全体の密度を入射時の分率で配分する。

$$
\begin{aligned}
N_2&=N f,\\
N_2/(\mathrm{m^{-3}})&\simeq \left(1\times 10^{22}\right) \cdot \left(\frac{3}{5}\right)\simeq 6\times 10^{21}.
\end{aligned}
$$

下準位からの吸収を断面積比で重み付けして引く。

$$
\begin{aligned}
\Delta N&=N_2 - r \left(N - N_2\right),\\
\Delta N/(\mathrm{m^{-3}})&\simeq \left(6\times 10^{21}\right) - \left(2\right) \cdot \left(\left(1\times 10^{22}\right) - \left(6\times 10^{21}\right)\right)\simeq -2\times 10^{21}.
\end{aligned}
$$

同じ線での誘導放出断面積を掛ける。

$$
\begin{aligned}
g&=\Delta N \sigma,\\
g/(\mathrm{m^{-1}})&\simeq \left(-2\times 10^{21}\right) \cdot \left(1\times 10^{-20}\right)=-20.
\end{aligned}
$$

上準位に過半数がいても吸収が勝る。
:::

:::example[縮退を含めて釣り合わせる]
上準位の分率を透明の条件にする。

入力は下準位の縮退度 $d_{1}=1$、上準位の縮退度 $d_{2}=2$、上準位の占有分率 $f=\frac{2}{3}$、総粒子密度 $N\simeq 1\times 10^{22}\,\mathrm{m^{-3}}$、ピークでの誘導放出断面積 $\sigma\simeq 1\times 10^{-20}\,\mathrm{m^2}$ である。

同じ線形状の Einstein 係数の関係を使う。

$$
\begin{aligned}
r&=\frac{d_{2}}{d_{1}},\\
r&= \frac{\left(2\right)}{\left(1\right)}=2.
\end{aligned}
$$

上と下の占有を縮退一個当たりで等しくする。

$$
\begin{aligned}
f_{\rm tr}&=\frac{d_{2}}{d_{1} + d_{2}},\\
f_{\rm tr}&= \frac{\left(2\right)}{\left(1\right) + \left(2\right)}=\frac{2}{3}.
\end{aligned}
$$

全体の密度を入射時の分率で配分する。

$$
\begin{aligned}
N_2&=N f,\\
N_2/(\mathrm{m^{-3}})&\simeq \left(1\times 10^{22}\right) \cdot \left(\frac{2}{3}\right)\simeq 6.6667\times 10^{21}.
\end{aligned}
$$

下準位からの吸収を断面積比で重み付けして引く。

$$
\begin{aligned}
\Delta N&=N_2 - r \left(N - N_2\right),\\
\Delta N/(\mathrm{m^{-3}})&\simeq \left(1\times 10^{22}\right)\left[\frac23-2\left(1-\frac23\right)\right]=0.
\end{aligned}
$$

同じ線での誘導放出断面積を掛ける。

$$
\begin{aligned}
g&=\Delta N \sigma,\\
g/(\mathrm{m^{-1}})&\simeq \left(0\right) \cdot \left(1\times 10^{-20}\right)=0.
\end{aligned}
$$

誘導による増減が釣り合う。自然放出の率はゼロではない。
:::

## スペクトル線の広がり（均一広がり・不均一広がり）

### 同種の粒子に共通の線幅
単一の原子の偏極振幅が $e^{-t/T_2}e^{-2\pi i\nu_0t}$ で減衰するとする。片側 Fourier 積分の実部は $T_2^{-1}/[T_2^{-2}+(2\pi(\nu-\nu_0))^2]$。面積を 1 に規格化した Lorentz 形で、その周波数全半値幅は $\Delta\nu_{\rm L}=1/(\pi T_2)$ である。全ての原子が同じ線形状を持つこの機構を :term[均一広がり]{en="homogeneous broadening" yomi="きんいつひろがり"} と呼ぶ。寿命 $T_1$ だけが原因なら、上準位振幅は確率の平方根 $e^{-t/(2T_1)}$ で減るため、下準位が安定な場合の $T_2$ は $2T_1$。衝突で位相が乱れると $T_2$ はさらに短くなる。

### 原子ごとに中心が異なる場合
一方向の速度 $v$ が平均ゼロ、分散 $k_BT/m$ の Maxwell 分布を持つ気体を考える。非相対論近似で離調は $\nu_0 v/c$ に比例するため、周波数も分散 $\sigma_\nu^2=(\nu_0/c)^2 k_BT/m$ の Gaussian に従う。指数が $\ln2$ だけ下がる二点の間隔が $\Delta\nu_{\mathrm D}=2\sqrt{2\ln2}\,\sigma_\nu$ である。これが :term[不均一広がり]{en="inhomogeneous broadening" yomi="ふきんいつひろがり"} の一例である。個別に Lorentz 線を持てば、全体の線は中心分布と個別の線の畳み込みになり、Voigt 形を得る。狭帯域な光で飽和を起こすと、均一線では全原子の反転が影響を受けるが、広い不均一線では主に共鳴する部分が減少する。

:::example[コヒーレンスと幅]
$T_2=10\,\mathrm{ns}$ の均一線では $\Delta\nu_{\rm L}\simeq31.8\,\mathrm{MHz}$。$T_2$ が二倍なら線幅は半分になる。
:::
:::example[温度の依存性]
Doppler 幅だけを比べると、温度を四倍にすれば標準偏差が二倍になる。質量を四倍にして温度を固定すれば幅は半分になる。
:::

### 線形状の面積とピーク
面積を 1 とした Lorentz 線は、周波数の片側半値幅 $\gamma_\nu$ で $L(\Delta\nu)=\gamma_\nu/[\pi(\gamma_\nu^2+\Delta\nu^2)]$ と書ける。$\int_{-\infty}^\infty\dd x/(1+x^2)=\pi$ より規格化され、ピークは $1/(\pi\gamma_\nu)$。同じ遷移強度を保ったまま衝突で線幅が二倍になると、断面積のピークは半分になる。不均一な速度分布の全てを同時に励起しない測定では、速度群による反転差も積分しなければならない。

### 飽和で均一線が広がる
共鳴近くで誘導率が $W(\Delta)=W_0/[1+(\Delta T_2)^2]$ とする。$\Delta$ は角周波数の離調である。閉二準位の定常反転差の絶対値は、光なしの値を $1+2W(\Delta)T_1$ で割った値になる。応答の吸収断面積にも同じ Lorentz 分母が掛かるので、吸収の強度に対する比例係数は $[1+2W_0T_1+(\Delta T_2)^2]^{-1}$ に比例する。従って低強度の半値幅に比べ、強い駆動では幅が $\sqrt{1+2W_0T_1}$ 倍になる。均一な減衰時間が変わったと仮定する必要はなく、中心部が先に飽和して相対的に平坦になる。

:::example[自然幅から位相緩和の寄与を分ける]
安定した下準位と寿命 $T_1=50\,\mathrm{ns}$ の上準位の組の自然な全半値幅は $1/(2\pi T_1)\simeq3.18\,\mathrm{MHz}$。同じ粒子の均一線が $31.8\,\mathrm{MHz}$ であれば $T_2$ は約 $10\,\mathrm{ns}$ であり、寿命のみの $100\,\mathrm{ns}$ より短い。
:::

### Gaussian な Doppler 線の幅から速度分布を読む

原子の視線方向速度 $v$ に対し $|v|\ll c$ なら周波数のずれは $\delta\nu\simeq\nu_0 v/c$。温度 $T$、粒子質量 $M$ で速度の分布が $\exp[-Mv^2/(2k_BT)]$ に比例するなら、変数変換で周波数分布の標準偏差は $(\nu_0/c)\sqrt{k_BT/M}$。半値全幅はこの値の $2\sqrt{2\ln2}$ 倍である。共通の均一な Lorentz 幅が別にある場合、同じ原子の均一線形状を速度分布に沿って平行移動して足すため、測定線は二つの畳み込みになる。いずれかが他方より十分広く、中心付近を測る場合に限り太い方の単独幅で近似できる。

:::example[暖かい原子気体]
可視域で軽めの原子の一次元速度を評価する。

入力は中心周波数 $\nu_0\simeq 5\times 10^{14}\,\mathrm{Hz}$、温度 $T=300\,\mathrm K$、一粒子の質量 $M\simeq 3.3\times 10^{-26}\,\mathrm{kg}$、均一線の半値全幅 $\Delta\nu_H=10000000\,\mathrm{Hz}$ である。

速度の一成分に Maxwell 分布を使う。

$$
\begin{aligned}
\sigma_v&=\sqrt{\frac{k_B T}{M}},\\
\sigma_v/(\mathrm{m/s})&\simeq \frac{\sqrt{13806490} \cdot \sqrt{\frac{\left(300\right)}{\left(3.3\times 10^{-26}\right)}}}{1000000000000000}\simeq 354.28.
\end{aligned}
$$

低速での周波数変換を施す。

$$
\begin{aligned}
\sigma_\nu&=\frac{\sigma_v\nu_0}{c},\\
\sigma_\nu/(\mathrm{Hz})&\simeq \frac{\left(354.2789091\right) \cdot \left(5\times 10^{14}\right)}{299792458}\simeq 5.9087\times 10^{8}.
\end{aligned}
$$

Gaussian の半値点の間隔を求める。

$$
\begin{aligned}
\Delta\nu_D&=2 \sqrt{2} \sigma_\nu \sqrt{\log{\left(2 \right)}},\\
\Delta\nu_D/(\mathrm{Hz})&\simeq 2 \cdot \sqrt{2} \cdot \left(5.908736189\times 10^{8}\right) \cdot \sqrt{\log{\left(2 \right)}}\simeq 1.3914\times 10^{9}.
\end{aligned}
$$

単一線形状を採用できる範囲かを比較する。

$$
\begin{aligned}
g&=\frac{\Delta\nu_D}{\Delta\nu_H},\\
g&\simeq \frac{\left(1.391401042\times 10^{9}\right)}{\left(10000000\right)}\simeq 139.14.
\end{aligned}
$$

線の中心付近は大きい Doppler 幅に支配される。
:::

:::example[重い原子へ換える]
同じ光学周波数で質量を四倍にした模型を比較する。

入力は中心周波数 $\nu_0\simeq 5\times 10^{14}\,\mathrm{Hz}$、温度 $T=300\,\mathrm K$、一粒子の質量 $M\simeq 1.32\times 10^{-25}\,\mathrm{kg}$、均一線の半値全幅 $\Delta\nu_H=10000000\,\mathrm{Hz}$ である。

速度の一成分に Maxwell 分布を使う。

$$
\begin{aligned}
\sigma_v&=\sqrt{\frac{k_B T}{M}},\\
\sigma_v/(\mathrm{m/s})&\simeq \frac{\sqrt{13806490} \cdot \sqrt{\frac{\left(300\right)}{\left(1.32\times 10^{-25}\right)}}}{1000000000000000}\simeq 177.14.
\end{aligned}
$$

低速での周波数変換を施す。

$$
\begin{aligned}
\sigma_\nu&=\frac{\sigma_v\nu_0}{c},\\
\sigma_\nu/(\mathrm{Hz})&\simeq \frac{\left(177.1394546\right) \cdot \left(5\times 10^{14}\right)}{299792458}\simeq 2.9544\times 10^{8}.
\end{aligned}
$$

Gaussian の半値点の間隔を求める。

$$
\begin{aligned}
\Delta\nu_D&=2 \sqrt{2} \sigma_\nu \sqrt{\log{\left(2 \right)}},\\
\Delta\nu_D/(\mathrm{Hz})&\simeq 2 \cdot \sqrt{2} \cdot \left(2.954368094\times 10^{8}\right) \cdot \sqrt{\log{\left(2 \right)}}\simeq 6.957\times 10^{8}.
\end{aligned}
$$

単一線形状を採用できる範囲かを比較する。

$$
\begin{aligned}
g&=\frac{\Delta\nu_D}{\Delta\nu_H},\\
g&\simeq \frac{\left(6.957005209\times 10^{8}\right)}{\left(10000000\right)}\simeq 69.57.
\end{aligned}
$$

温度を保ったまま幅が半分になる。
:::

:::example[十分冷却した原子]
均一幅を固定して運動を抑える。

入力は中心周波数 $\nu_0\simeq 5\times 10^{14}\,\mathrm{Hz}$、温度 $T\simeq 1\times 10^{-4}\,\mathrm K$、一粒子の質量 $M\simeq 3.3\times 10^{-26}\,\mathrm{kg}$、均一線の半値全幅 $\Delta\nu_H=10000000\,\mathrm{Hz}$ である。

速度の一成分に Maxwell 分布を使う。

$$
\begin{aligned}
\sigma_v&=\sqrt{\frac{k_B T}{M}},\\
\sigma_v/(\mathrm{m/s})&\simeq \frac{\sqrt{13806490} \cdot \sqrt{\frac{\left(1\times 10^{-4}\right)}{\left(3.3\times 10^{-26}\right)}}}{1000000000000000}\simeq 0.20454.
\end{aligned}
$$

低速での周波数変換を施す。

$$
\begin{aligned}
\sigma_\nu&=\frac{\sigma_v\nu_0}{c},\\
\sigma_\nu/(\mathrm{Hz})&\simeq \frac{\left(0.2045430236\right) \cdot \left(5\times 10^{14}\right)}{299792458}\simeq 3.4114\times 10^{5}.
\end{aligned}
$$

Gaussian の半値点の間隔を求める。

$$
\begin{aligned}
\Delta\nu_D&=2 \sqrt{2} \sigma_\nu \sqrt{\log{\left(2 \right)}},\\
\Delta\nu_D/(\mathrm{Hz})&\simeq 2 \cdot \sqrt{2} \cdot \left(3.411410429\times 10^{5}\right) \cdot \sqrt{\log{\left(2 \right)}}\simeq 8.0333\times 10^{5}.
\end{aligned}
$$

単一線形状を採用できる範囲かを比較する。

$$
\begin{aligned}
g&=\frac{\Delta\nu_D}{\Delta\nu_H},\\
g&\simeq \frac{\left(8.03325766\times 10^{5}\right)}{\left(10000000\right)}\simeq 0.080333.
\end{aligned}
$$

Doppler 幅が小さくなるため、均一線を無視できない。
:::

## 光増幅と利得飽和

### 断面積と強度の伝播
ほぼ一方向の単色場の光子束を $\Phi_{\gamma}=I/(h\nu)$ とする。断面積 $\sigma_e,\sigma_a$ を、粒子一つ当たりの誘導放出率 $\sigma_e\Phi_\gamma$、吸収率 $\sigma_a\Phi_\gamma$ によって定義する。層 $\dd z$ で増減するパワーから
$$
\dv{I}{z}=(g-\alpha)I,\qquad g=\sigma_e N_2-\sigma_a N_1
\label{eq-optics-laser-gain}
$$
を得る。$\sigma$ は $\mathrm{m^2}$、$g,\alpha$ は $\mathrm{m^{-1}}$。弱信号で $g$ が空間的に一定なら、長さ $L$ のパワー利得は $e^{(g-\alpha)L}$ である。自発放出はこの式では入力または別の源項として扱う。

### 緩和と定常飽和
等縮退で $\sigma_e=\sigma_a=\sigma$ とし、ポンプが光なしで反転差 $D_0$ を保つとする。差 $D=N_2-N_1$ は一回の誘導遷移で二粒子分だけ変わるので、単純な均一緩和模型は
$$
\dot D=(D_0-D)/T_1-2\sigma I D/(h\nu).
$$
定常値は $D=D_0/(1+I/I_s)$、$I_s=h\nu/(2\sigma T_1)$。これを :term[利得飽和]{en="gain saturation" yomi="りとくほうわ"} という。下準位がすぐ空く四準位の近似では反転が上準位だけで決まり、一回の放出の変化は一粒子分なので、この式の係数 2 は外れる。損失ゼロで $g=g_0/(1+I/I_s)$ を入れて変数分離すると
$$
\ln(I_{\rm out}/I_{\rm in})+(I_{\rm out}-I_{\rm in})/I_s=g_0 L .
$$
強度一定という近似を用いなくても増幅の圧縮を求められる。

:::example[小信号利得]
$g_0=2\,\mathrm{m^{-1}}$、$\alpha=0.5\,\mathrm{m^{-1}}$、$L=1\,\mathrm m$ なら弱信号のパワーは $e^{1.5}\simeq4.48$ 倍になる。
:::
:::example[半分の利得]
局所強度 $I=I_s$ で $D=D_0/2$。同じ媒質でも入射から出射まで $I$ が成長するので、一般には入口での利得を全長へ掛けられない。
:::

### 飽和で利得が減る間のパルス増幅
入力パルスが $T_1$ より十分短く、位相緩和時間 $T_2$ より十分長い四準位模型を考える。下準位は速く空き、パルス中のポンプ、自然崩壊、群速度分散、伝播損失を無視する。パルスの帯域内では断面積 $\sigma$ を一定とみなす。反転 $N(z,t)$ は $\partial_t N=-\sigma IN/(h\nu)$、パワーの伝播は遅延時刻で $\partial_z I=\sigma NI$ とする。積算フルエンスを $\mathcal E(z,t)=\int_{-\infty}^t I(z,t')\dd t'$ とし、$\mathcal E_s=h\nu/\sigma$ と置く。時間方程式から $N=N_0(z)e^{-\mathcal E/\mathcal E_s}$。これを強度の $z$ 微分に代入して時間積分すると
$$
\partial_z\mathcal E=\sigma N_0(z)\mathcal E_s(1-e^{-\mathcal E/\mathcal E_s}).
$$
$x=e^{\mathcal E/\mathcal E_s}-1$ で $x_z=\sigma N_0 x$。従って全パルスの入出力フルエンスは $\mathcal E_{\rm out}=\mathcal E_s\ln[1+G_0(e^{\mathcal E_{\rm in}/\mathcal E_s}-1)]$、$G_0=e^{\int\sigma N_0\dd z}$ となる。連続照射での飽和強度とは異なり、短パルスでは面積当たりに取り出したエネルギーが効く。

:::example[短パルスのフルエンス利得]
$G_0=3$、$\mathcal E_{\rm in}=\mathcal E_s\ln 2$ なら $\mathcal E_{\rm out}=\mathcal E_s\ln4=2\mathcal E_{\rm in}$。弱信号の 3 倍より利得が下がる。入射パルスの後半ほど上準位の占有が減っているためである。
:::

### 不均一線の一部が飽和する

原子の中心周波数を $\nu_a$、弱信号での反転分布を $D_0(\nu_a)$、各原子群の断面積を $\sigma(\nu-\nu_a)$ とする。ポンプと緩和の条件を固定し、原子群同士の交換を無視する。周波数 $\nu_p$ の強い場が群 $\nu_a$ に与える誘導率は $W_p(\nu_a)=I_p\sigma(\nu_p-\nu_a)/(h\nu_p)$ である。等縮退の反転差の定常式を群ごとに解けば、弱い探針光の利得は
$$
g_{\rm probe}(\nu)=\int \sigma(\nu-\nu_a)
 \frac{D_0(\nu_a)}{1+2T_1W_p(\nu_a)}\dd\nu_a .
$$
広い不均一線では、飽和した原子群に応答する探針光の周波数の近くに利得のくぼみができる。この :term[スペクトルホールバーニング]{en="spectral hole burning" yomi="すぺくとるほーるばーにんぐ"} が、同じ共振器で周波数の離れたモードが発振できる理由の一つである。原子群が一つなら分母は積分の外に出て、線全体の利得が同じ反転の低下で抑えられる。

### 一様な飽和増幅器の定常出力を反転する

無損失で定常の一次元式 $\dd I/\dd z=g_0 I/(1+I/I_s)$ に対し、$x=I/I_s$ とすると $(x^{-1}+1)\dd x=g_0\dd z$。従って $\ln x_{\rm out}+x_{\rm out}=\ln x_{\rm in}+x_{\rm in}+g_0L$。左辺は $x>0$ で単調増加し、解は一つである。指数を取った式 $x_{\rm out}e^{x_{\rm out}}=x_{\rm in}e^{x_{\rm in}+g_0L}$ を反転する関数が Lambert の $W$ 関数の正実枝である。小信号の指数利得を固定したまま入射を増やすと、この反転式では強度の倍率が下がる。連続波の定常飽和を仮定しており、短パルスのフルエンス模型には置き換えない。

:::example[弱い種が増幅される]
出力でも飽和強度を越えないほど小さい入力を選ぶ。

入力は飽和強度で測った入力 $x_{\rm in}=\frac{1}{1000}$、非飽和の小信号利得係数 $g_0=20\,\mathrm{m^{-1}}$、増幅長 $L=\frac{1}{10}\,\mathrm m$ である。

指数利得の引数を求める。

$$
\begin{aligned}
b&=L g_0,\\
b&= \left(\frac{1}{10}\right) \cdot \left(20\right)=2.
\end{aligned}
$$

比較用に飽和を除いた輸送式を積分する。

$$
\begin{aligned}
x_{\rm lin}&=x_{\rm in} e^{b},\\
x_{\rm lin}&= \left(\frac{1}{1000}\right) \cdot e^{\left(2\right)}\simeq 0.0073891.
\end{aligned}
$$

有限長で途中の強度による局所利得の低下まで含める。

$$
\begin{aligned}
x_{\rm out}&=W\left(x_{\rm in} e^{b + x_{\rm in}}\right),\\
x_{\rm out}&= W\left(\left(\frac{1}{1000}\right) \cdot e^{\left(2\right) + \left(\frac{1}{1000}\right)}\right)\simeq 0.0073423.
\end{aligned}
$$

同じ入力で出力を割る。

$$
\begin{aligned}
R&=\frac{x_{\rm out}}{x_{\rm in}},\\
R&\simeq \frac{\left(0.007342340489\right)}{\left(\frac{1}{1000}\right)}\simeq 7.3423.
\end{aligned}
$$

小信号近似がどこまで使えるかを判定する。

$$
\begin{aligned}
D&=\frac{x_{\rm out}}{x_{\rm lin}},\\
D&\simeq \frac{\left(0.007342340489\right)}{\left(0.007389056099\right)}\simeq 0.99368.
\end{aligned}
$$

非飽和の予測に近い。
:::

:::example[増幅途中で飽和する]
同じ媒質で入力を強くする。

入力は飽和強度で測った入力 $x_{\rm in}=\frac{1}{5}$、非飽和の小信号利得係数 $g_0=20\,\mathrm{m^{-1}}$、増幅長 $L=\frac{1}{10}\,\mathrm m$ である。

指数利得の引数を求める。

$$
\begin{aligned}
b&=L g_0,\\
b&= \left(\frac{1}{10}\right) \cdot \left(20\right)=2.
\end{aligned}
$$

比較用に飽和を除いた輸送式を積分する。

$$
\begin{aligned}
x_{\rm lin}&=x_{\rm in} e^{b},\\
x_{\rm lin}&= \left(\frac{1}{5}\right) \cdot e^{\left(2\right)}\simeq 1.4778.
\end{aligned}
$$

有限長で途中の強度による局所利得の低下まで含める。

$$
\begin{aligned}
x_{\rm out}&=W\left(x_{\rm in} e^{b + x_{\rm in}}\right),\\
x_{\rm out}&= W\left(\left(\frac{1}{5}\right) \cdot e^{\left(2\right) + \left(\frac{1}{5}\right)}\right)\simeq 0.8061.
\end{aligned}
$$

同じ入力で出力を割る。

$$
\begin{aligned}
R&=\frac{x_{\rm out}}{x_{\rm in}},\\
R&\simeq \frac{\left(0.8061042678\right)}{\left(\frac{1}{5}\right)}\simeq 4.0305.
\end{aligned}
$$

小信号近似がどこまで使えるかを判定する。

$$
\begin{aligned}
D&=\frac{x_{\rm out}}{x_{\rm lin}},\\
D&\simeq \frac{\left(0.8061042678\right)}{\left(1.47781122\right)}\simeq 0.54547.
\end{aligned}
$$

出力近くでは小信号利得のまま進まない。
:::

:::example[入射時から強く飽和]
入力を飽和強度より高くする。

入力は飽和強度で測った入力 $x_{\rm in}=5$、非飽和の小信号利得係数 $g_0=20\,\mathrm{m^{-1}}$、増幅長 $L=\frac{1}{10}\,\mathrm m$ である。

指数利得の引数を求める。

$$
\begin{aligned}
b&=L g_0,\\
b&= \left(\frac{1}{10}\right) \cdot \left(20\right)=2.
\end{aligned}
$$

比較用に飽和を除いた輸送式を積分する。

$$
\begin{aligned}
x_{\rm lin}&=x_{\rm in} e^{b},\\
x_{\rm lin}&= \left(5\right) \cdot e^{\left(2\right)}\simeq 36.945.
\end{aligned}
$$

有限長で途中の強度による局所利得の低下まで含める。

$$
\begin{aligned}
x_{\rm out}&=W\left(x_{\rm in} e^{b + x_{\rm in}}\right),\\
x_{\rm out}&= W\left(\left(5\right) \cdot e^{\left(2\right) + \left(5\right)}\right)\simeq 6.7064.
\end{aligned}
$$

同じ入力で出力を割る。

$$
\begin{aligned}
R&=\frac{x_{\rm out}}{x_{\rm in}},\\
R&\simeq \frac{\left(6.706378782\right)}{\left(5\right)}\simeq 1.3413.
\end{aligned}
$$

小信号近似がどこまで使えるかを判定する。

$$
\begin{aligned}
D&=\frac{x_{\rm out}}{x_{\rm lin}},\\
D&\simeq \frac{\left(6.706378782\right)}{\left(36.94528049\right)}\simeq 0.18152.
\end{aligned}
$$

強度の増分が小信号の指数則より小さくなる。
:::

## 3 準位・4 準位レーザーと速度方程式

### 三準位系の粒子収支
基底 1 からポンプ準位 3 に率 $W$ で励起し、3 は速く準安定な上準位 2 へ落ちるとする。ポンプ帯を空とみなし、等縮退、誘導放出がまだ弱いとき、$\dot N_2=W(N-N_2)-N_2/T_1$。従って $N_2=W T_1 N/(1+W T_1)$。下レーザー準位が基底なので反転には $N_2>N/2$、すなわち $WT_1>1$ が必要である。共振器の有限損失を補うにはさらに大きい反転が要る。

### 四準位系の粒子収支
基底 $0$ から高いポンプ準位への励起は速く上レーザー準位へ移る。下レーザー準位 1 が寿命 $T_{10}\ll T_1$ で基底へ落ちるため、弱信号で $N_1\simeq N_2T_{10}/T_1\ll N_2$（ここでは上準位の自然崩壊が 1 に行くとする）。有効ポンプ供給を一定 $R_p$ とするモデルの速度方程式は
$$
\dot N_2=R_p-N_2/T_1-v_g\sigma N_2 S,\qquad
\dot S=(v_g\Gamma_{\rm ov}\sigma N_2-\tau_{\rm ph}^{-1})S
$$
である。$S$ はモード体積 $V_m$ 当たりの光子数、$N_2$ は活性体積 $V_a$ 当たりの上準位粒子数で、$v_g$ は群速度である。場と反転をそれぞれの体積内で一様とみなすと、重なり係数は $\Gamma_{\rm ov}=V_a/V_m$。二式の誘導項にそれぞれ $V_a,V_m$ を掛ければ、消える上準位粒子数と増える光子数が一致する。ここでは光子の式の自発放出を省いている。強い駆動では下準位が十分空であるか、誘導項も含む式で確かめる必要がある。

:::example[三準位の反転]
$WT_1=3$ なら $N_2=3N/4$。反転差は $N/2$ である。
:::
:::example[空になりやすい下準位]
$T_{10}/T_1=10^{-3}$ で弱信号なら、下準位は上準位の $10^{-3}$ 程度である。三準位のように全粒子の半数を励起する条件は要らない。
:::

### ポンプ遷移が有限速度で戻る効果
三準位モデルでポンプ上準位 3 の寿命を $T_{32}$ として式に残すと、上の弱信号の収支は $\dot N_3=W N_1-N_3/T_{32}$、$\dot N_2=N_3/T_{32}-N_2/T_1$。定常で $N_3=WT_{32}N_1$、$N_2=WT_1N_1$ である。従って全数で割った $N_2$ の割合は $WT_1/[1+W(T_1+T_{32})]$。ポンプ上準位も含めると反転の条件は $N_2>N_1$ のままだが、$N_2>N/2$ と同値ではなくなる。後者には $N_3$ を無視する仮定が要る。

### 発振枝を横切る局所安定性
自発放出を省いた四準位式で、非発振解は $S=0,N_2=R_p T_1$ である。光子密度に対する微小摂動の成長率は $v_g\Gamma_{\rm ov}\sigma R_p T_1-1/\tau_{\rm ph}$。ポンプがしきい値を過ぎると、この非発振解は任意に小さな正の光子密度に対して不安定になる。一方、$S>0$ の枝では線形化行列のトレースは $-(1/T_1+v_g\sigma S)<0$、行列式は正であるため、二固有値は両方負の実部を持つ。二つの直線の出力特性は、この局所安定性が交換する点でつながる。実際には自発放出があるため下側でも光出力が厳密にゼロではなく、しきい値近くの曲線の角が丸くなる。

:::example[有限ポンプ準位の占有]
$WT_1=2,WT_{32}=1$ なら $N_1/N=1/4,N_2/N=1/2,N_3/N=1/4$。上準位は全数の半分だが、下準位との差 $N/4$ による反転がすでに存在する。
:::

### 四準位の下準位の抜けを定常収支に入れる

上準位へ粒子が単位体積・単位時間当たり $P$ 注入され、上準位寿命が $\tau_u$、レーザー下準位の排出時間が $\tau_l$ とする。上の自然減衰はすべて下準位へ入り、誘導遷移を無視できる小信号とする。$\dot N_u=P-N_u/\tau_u$、$\dot N_l=N_u/\tau_u-N_l/\tau_l$ の定常状態は $N_u=P\tau_u$、$N_l=P\tau_l$。この注入は貯蔵準位の枯渇を無視した形である。等縮退で、正味反転を得るには $\tau_l<\tau_u$。速く下準位を空けるという四準位の利点は、同じポンプで下準位への自然供給をこの時間の比だけ抑えることにある。発振が始まった後は誘導項を上下の収支に加える。

:::example[下準位を速く空ける]
上準位よりはるかに速い緩和を持つ場合を取る。

入力は上準位への注入率 $P\simeq 1\times 10^{25}\,\mathrm{m^{-3}s^{-1}}$、上準位の寿命 $\tau_u\simeq 1\times 10^{-4}\,\mathrm s$、下準位の排出時間 $\tau_l\simeq 1\times 10^{-8}\,\mathrm s$、誘導遷移の断面積 $\sigma\simeq 1\times 10^{-20}\,\mathrm{m^2}$ である。

供給と自然減衰を釣り合わせる。

$$
\begin{aligned}
N_u&=P \tau_u,\\
N_u/(\mathrm{m^{-3}})&\simeq \left(1\times 10^{25}\right) \cdot \left(1\times 10^{-4}\right)\simeq 1\times 10^{21}.
\end{aligned}
$$

上からの供給と下からの排出を釣り合わせる。

$$
\begin{aligned}
N_l&=P \tau_l,\\
N_l/(\mathrm{m^{-3}})&\simeq \left(1\times 10^{25}\right) \cdot \left(1\times 10^{-8}\right)\simeq 1\times 10^{17}.
\end{aligned}
$$

等縮退で吸収に入る下準位密度を引く。

$$
\begin{aligned}
\Delta N&=- N_l + N_u,\\
\Delta N/(\mathrm{m^{-3}})&\simeq - \left(1\times 10^{17}\right) + \left(1\times 10^{21}\right)\simeq 9.999\times 10^{20}.
\end{aligned}
$$

四準位を空の下準位で近似するときの誤差を求める。

$$
\begin{aligned}
R&=\frac{N_l}{N_u},\\
R&\simeq \frac{\left(1\times 10^{17}\right)}{\left(1\times 10^{21}\right)}\simeq 1\times 10^{-4}.
\end{aligned}
$$

線の中心で正味の反転を光の利得へ換算する。

$$
\begin{aligned}
g&=\Delta N \sigma,\\
g/(\mathrm{m^{-1}})&\simeq \left(9.999\times 10^{20}\right) \cdot \left(1\times 10^{-20}\right)=\frac{9999}{1000}.
\end{aligned}
$$

下準位の吸収は上準位に対して小さい。
:::

:::example[下準位の緩和が遅い]
注入は同じで下準位の寿命を伸ばす。

入力は上準位への注入率 $P\simeq 1\times 10^{25}\,\mathrm{m^{-3}s^{-1}}$、上準位の寿命 $\tau_u\simeq 1\times 10^{-4}\,\mathrm s$、下準位の排出時間 $\tau_l\simeq 5\times 10^{-5}\,\mathrm s$、誘導遷移の断面積 $\sigma\simeq 1\times 10^{-20}\,\mathrm{m^2}$ である。

供給と自然減衰を釣り合わせる。

$$
\begin{aligned}
N_u&=P \tau_u,\\
N_u/(\mathrm{m^{-3}})&\simeq \left(1\times 10^{25}\right) \cdot \left(1\times 10^{-4}\right)\simeq 1\times 10^{21}.
\end{aligned}
$$

上からの供給と下からの排出を釣り合わせる。

$$
\begin{aligned}
N_l&=P \tau_l,\\
N_l/(\mathrm{m^{-3}})&\simeq \left(1\times 10^{25}\right) \cdot \left(5\times 10^{-5}\right)\simeq 5\times 10^{20}.
\end{aligned}
$$

等縮退で吸収に入る下準位密度を引く。

$$
\begin{aligned}
\Delta N&=- N_l + N_u,\\
\Delta N/(\mathrm{m^{-3}})&\simeq - \left(5\times 10^{20}\right) + \left(1\times 10^{21}\right)\simeq 5\times 10^{20}.
\end{aligned}
$$

四準位を空の下準位で近似するときの誤差を求める。

$$
\begin{aligned}
R&=\frac{N_l}{N_u},\\
R&\simeq \frac{\left(5\times 10^{20}\right)}{\left(1\times 10^{21}\right)}=\frac{1}{2}.
\end{aligned}
$$

線の中心で正味の反転を光の利得へ換算する。

$$
\begin{aligned}
g&=\Delta N \sigma,\\
g/(\mathrm{m^{-1}})&\simeq \left(5\times 10^{20}\right) \cdot \left(1\times 10^{-20}\right)=5.
\end{aligned}
$$

上の占有は同じでも正味利得は半分へ下がる。
:::

:::example[下準位に滞留する]
さらに下準位の排出を遅くする。

入力は上準位への注入率 $P\simeq 1\times 10^{25}\,\mathrm{m^{-3}s^{-1}}$、上準位の寿命 $\tau_u\simeq 1\times 10^{-4}\,\mathrm s$、下準位の排出時間 $\tau_l=\frac{1}{5000}\,\mathrm s$、誘導遷移の断面積 $\sigma\simeq 1\times 10^{-20}\,\mathrm{m^2}$ である。

供給と自然減衰を釣り合わせる。

$$
\begin{aligned}
N_u&=P \tau_u,\\
N_u/(\mathrm{m^{-3}})&\simeq \left(1\times 10^{25}\right) \cdot \left(1\times 10^{-4}\right)\simeq 1\times 10^{21}.
\end{aligned}
$$

上からの供給と下からの排出を釣り合わせる。

$$
\begin{aligned}
N_l&=P \tau_l,\\
N_l/(\mathrm{m^{-3}})&\simeq \left(1\times 10^{25}\right) \cdot \left(\frac{1}{5000}\right)\simeq 2\times 10^{21}.
\end{aligned}
$$

等縮退で吸収に入る下準位密度を引く。

$$
\begin{aligned}
\Delta N&=- N_l + N_u,\\
\Delta N/(\mathrm{m^{-3}})&\simeq - \left(2\times 10^{21}\right) + \left(1\times 10^{21}\right)\simeq -1\times 10^{21}.
\end{aligned}
$$

四準位を空の下準位で近似するときの誤差を求める。

$$
\begin{aligned}
R&=\frac{N_l}{N_u},\\
R&\simeq \frac{\left(2\times 10^{21}\right)}{\left(1\times 10^{21}\right)}=2.
\end{aligned}
$$

線の中心で正味の反転を光の利得へ換算する。

$$
\begin{aligned}
g&=\Delta N \sigma,\\
g/(\mathrm{m^{-1}})&\simeq \left(-1\times 10^{21}\right) \cdot \left(1\times 10^{-20}\right)=-10.
\end{aligned}
$$

上下密度の差の符号が逆になり、吸収を示す。
:::

## 発振条件と出力特性

### 一往復で減衰と利得を比べる
活性長 $\ell$ の直線共振器を考える。鏡のパワー反射率を $R_1,R_2$、共振器全長を $L$、分布損失を $\alpha$ とすると、小信号の一往復の生存率は $R_1R_2e^{2g\ell-2\alpha L}$。発振の :term[しきい値]{en="laser threshold" yomi="しきいち"} はこれが 1 になる条件で、
$$
g_{\rm th}=\frac{\alpha L}{\ell}+\frac{1}{2\ell}\ln\frac{1}{R_1R_2}.
\label{eq-optics-laser-threshold}
$$
位相も往復で $2\pi$ の整数倍になる必要がある（:eqref[eq-optics-resonator-modes]）。小信号で複数モードがしきい値を超える場合、利得飽和による競合を調べる。単一モードの定常四準位模型で $S>0$ なら $N_2=N_{\mathrm{th}}=1/(v_g\Gamma_{\rm ov}\sigma\tau_{\rm ph})$ に固定される。第一式から
$$
S=\frac{R_p-N_{\mathrm{th}}/T_1}{v_g\sigma N_{\mathrm{th}}}.
$$
このため出力結合とモード体積が一定なら、出力パワーはしきい値を超えたポンプに比例する。ポンプの強さだけでは出力効率は決まらず、内部損失と取り出し損失の割合も必要である。

### 定常解近くの揺らぎ
$v_g\sigma=b$、定常値を $N_*,S_*>0$ とすると、線形化した変動は $\delta\dot N=-(1/T_1+bS_*)\delta N-bN_*\delta S$、$\delta\dot S=b\Gamma_{\rm ov}S_*\delta N$。消去して $\delta\ddot S+(1/T_1+bS_*)\delta\dot S+b^2\Gamma_{\rm ov}N_*S_*\delta S=0$ を得る。固有値が複素共役となる条件 $b^2\Gamma_{\rm ov}N_*S_*>(1/T_1+bS_*)^2/4$ で緩和振動が現れる。角周波数はこの二量の差の平方根、振幅の減衰率は $(1/T_1+bS_*)/2$ である。共振器寿命の定義は :eqref[eq-optics-photon-lifetime] を使う。

:::example[鏡のしきい値]
$R_1=1,R_2=0.90,\ell=0.10\,\mathrm m,\alpha=0$ なら $g_{\rm th}=\ln(10/9)/0.20\simeq0.527\,\mathrm{m^{-1}}$。
:::
:::example[ポンプ余裕]
$R_p$ が $2R_{\rm th}$ から $3R_{\rm th}$ に増えると、定常の $S$ は二倍になる。反転 $N_2$ 自体はこの模型では増えない。
:::

### 出力結合としきい値の兼ね合い
小損失では共振器内光子の消滅率を内部 $1/\tau_{\mathrm{int}}$ と出力 $1/\tau_{\mathrm{out}}$ の和で書ける。定常出力へ出る光子の割合は $(1/\tau_{\mathrm{out}})/(1/\tau_{\rm ph})$ となる。出力鏡の透過を増やすとこの割合は上がる一方、$\tau_{\rm ph}$ が短くなって維持に要る反転 $N_{\rm th}$ も増える。例えば出力結合を二倍にしたときの出力を計算するには、ポンプ余裕と取り出し割合の両方を更新する必要がある。線形化した速度方程式を使う共振器の雑音と動的安定性の扱いは、縦モード間の競合や揺動源を含むため :cite[siegman-lasers] を参照する。

:::example[増幅長と鏡の取り分]
小信号 $g=2\,\mathrm{m^{-1}}$、活性長 $0.20\,\mathrm m$、内部の一往復対数損失が 0.10 とする。一往復の材料利得の対数は 0.80。鏡の反射率の積が $e^{-0.70}$ なら、残る対数の和がゼロでちょうどしきい値である。積をさらに下げるなら利得も増やす必要がある。
:::

### 損失と鏡反射率から必要なポンプ余裕を求める

長さ $L$ の一様活性部を一往復に二度通り、小信号のパワー利得係数を $g_0$、内部損失係数を $\alpha$ とする。一往復の倍率は $R_1R_2 e^{2(g_0-\alpha)L}$。一を越える必要条件は $g_0>g_{\rm th}=\alpha+\ln(1/(R_1R_2))/(2L)$。この比較は小信号発振の開始を判定するものであり、単一モードの定常発振では飽和した利得がしきい値に下がる。局所利得を一様な $g(P)=g_0/(1+P/P_s)$ と近似できれば、定常の内部パワーは $P_s(g_0/g_{\rm th}-1)$ である。定在波の場所ごとに飽和が違う場合は、この代数式の $P$ を実際のモード強度に対応させるモデルが要る。

:::example[鏡での損失が小さい]
高反射率の入力鏡と部分透過の出力鏡を使う。

入力は第一鏡の反射率 $R_1=\frac{999}{1000}$、第二鏡の反射率 $R_2=\frac{19}{20}$、活性部の長さ $L=\frac{1}{20}\,\mathrm m$、内部損失係数 $\alpha=\frac{1}{5}\,\mathrm{m^{-1}}$、小信号利得係数 $g_0=2\,\mathrm{m^{-1}}$ である。

鏡での一往復損失を伝播係数へ換算する。

$$
\begin{aligned}
g_{\rm mir}&=- \frac{\log{\left(R_1 R_2 \right)}}{2 L},\\
g_{\rm mir}/(\mathrm{m^{-1}})&= - \frac{\log{\left(\left(\frac{999}{1000}\right) \cdot \left(\frac{19}{20}\right) \right)}}{2 \cdot \left(\frac{1}{20}\right)}\simeq 0.52294.
\end{aligned}
$$

活性部内の損失も補う。

$$
\begin{aligned}
g_{\rm th}&=\alpha + g_{\rm mir},\\
g_{\rm th}/(\mathrm{m^{-1}})&\simeq \left(\frac{1}{5}\right) + \left(0.5229379472\right)\simeq 0.72294.
\end{aligned}
$$

発振初期の光が小さいときの伝播を評価する。

$$
\begin{aligned}
G&=R_1 R_2 e^{L \left(- 2 \alpha + 2 g_0\right)},\\
G&= \left(\frac{999}{1000}\right) \cdot \left(\frac{19}{20}\right) \cdot e^{\left(\frac{1}{20}\right) \cdot \left(- 2 \cdot \left(\frac{1}{5}\right) + 2 \cdot \left(2\right)\right)}\simeq 1.1362.
\end{aligned}
$$

正なら飽和定常解のパワーを $P_s$ で割った値にもなる。

$$
\begin{aligned}
d&=\frac{g_0}{g_{\rm th}} - 1,\\
d&\simeq \frac{\left(2\right)}{\left(0.7229379472\right)} - 1\simeq 1.7665.
\end{aligned}
$$

小信号で一往復すると光が増え、飽和なしには定常にならない。
:::

:::example[出力鏡を大きく開く]
利得を保ち、出力鏡の反射率を下げる。

入力は第一鏡の反射率 $R_1=\frac{999}{1000}$、第二鏡の反射率 $R_2=\frac{4}{5}$、活性部の長さ $L=\frac{1}{20}\,\mathrm m$、内部損失係数 $\alpha=\frac{1}{5}\,\mathrm{m^{-1}}$、小信号利得係数 $g_0=2\,\mathrm{m^{-1}}$ である。

鏡での一往復損失を伝播係数へ換算する。

$$
\begin{aligned}
g_{\rm mir}&=- \frac{\log{\left(R_1 R_2 \right)}}{2 L},\\
g_{\rm mir}/(\mathrm{m^{-1}})&= - \frac{\log{\left(\left(\frac{999}{1000}\right) \cdot \left(\frac{4}{5}\right) \right)}}{2 \cdot \left(\frac{1}{20}\right)}\simeq 2.2414.
\end{aligned}
$$

活性部内の損失も補う。

$$
\begin{aligned}
g_{\rm th}&=\alpha + g_{\rm mir},\\
g_{\rm th}/(\mathrm{m^{-1}})&\simeq \left(\frac{1}{5}\right) + \left(2.241440516\right)\simeq 2.4414.
\end{aligned}
$$

発振初期の光が小さいときの伝播を評価する。

$$
\begin{aligned}
G&=R_1 R_2 e^{L \left(- 2 \alpha + 2 g_0\right)},\\
G&= \left(\frac{999}{1000}\right) \cdot \left(\frac{4}{5}\right) \cdot e^{\left(\frac{1}{20}\right) \cdot \left(- 2 \cdot \left(\frac{1}{5}\right) + 2 \cdot \left(2\right)\right)}\simeq 0.95682.
\end{aligned}
$$

正なら飽和定常解のパワーを $P_s$ で割った値にもなる。

$$
\begin{aligned}
d&=\frac{g_0}{g_{\rm th}} - 1,\\
d&\simeq \frac{\left(2\right)}{\left(2.441440516\right)} - 1\simeq -0.18081.
\end{aligned}
$$

取り出し割合を増やしたが一往復では減衰し、小信号発振は始まらない。
:::

:::example[内部損失も変わる]
出力鏡の透過を元に戻して内部損失を増やす。

入力は第一鏡の反射率 $R_1=\frac{999}{1000}$、第二鏡の反射率 $R_2=\frac{19}{20}$、活性部の長さ $L=\frac{1}{20}\,\mathrm m$、内部損失係数 $\alpha=1\,\mathrm{m^{-1}}$、小信号利得係数 $g_0=2\,\mathrm{m^{-1}}$ である。

鏡での一往復損失を伝播係数へ換算する。

$$
\begin{aligned}
g_{\rm mir}&=- \frac{\log{\left(R_1 R_2 \right)}}{2 L},\\
g_{\rm mir}/(\mathrm{m^{-1}})&= - \frac{\log{\left(\left(\frac{999}{1000}\right) \cdot \left(\frac{19}{20}\right) \right)}}{2 \cdot \left(\frac{1}{20}\right)}\simeq 0.52294.
\end{aligned}
$$

活性部内の損失も補う。

$$
\begin{aligned}
g_{\rm th}&=\alpha + g_{\rm mir},\\
g_{\rm th}/(\mathrm{m^{-1}})&\simeq \left(1\right) + \left(0.5229379472\right)\simeq 1.5229.
\end{aligned}
$$

発振初期の光が小さいときの伝播を評価する。

$$
\begin{aligned}
G&=R_1 R_2 e^{L \left(- 2 \alpha + 2 g_0\right)},\\
G&= \left(\frac{999}{1000}\right) \cdot \left(\frac{19}{20}\right) \cdot e^{\left(\frac{1}{20}\right) \cdot \left(- 2 \cdot \left(1\right) + 2 \cdot \left(2\right)\right)}\simeq 1.0489.
\end{aligned}
$$

正なら飽和定常解のパワーを $P_s$ で割った値にもなる。

$$
\begin{aligned}
d&=\frac{g_0}{g_{\rm th}} - 1,\\
d&\simeq \frac{\left(2\right)}{\left(1.522937947\right)} - 1\simeq 0.31325.
\end{aligned}
$$

発振できる場合も、飽和までの利得の余裕は減る。
:::

### 利得媒質が発振周波数を引き寄せる

均一広がりの一遷移と一共振器モードを考える。共振器の振幅を $a$、媒質分極の応答を $p$ とし、位相と規格化を選んで $\dot a=-(\kappa_c+i\omega_c)a+g p$、$\dot p=-(\gamma_\perp+i\omega_a)p+g'Da$ とする。$\kappa_c=1/(2\tau_{\rm ph})$ は場の減衰率、$\gamma_\perp=1/T_2$ は分極の減衰率、$D$ は反転である。$\omega_c,\omega_a$ はそれぞれ受動共振と原子の中心角周波数であり、正の実定数 $g,g'$ に比例して反転が小信号利得を生む。定常な正弦波 $a,p\propto e^{-i\omega t}$ を代入して $p$ を消去すると
$$
[\kappa_c+i(\omega_c-\omega)][\gamma_\perp+i(\omega_a-\omega)]=gg'D.
$$
右辺が実数なので、虚部をゼロとして
$$
\omega=\frac{\gamma_\perp\omega_c+\kappa_c\omega_a}{\gamma_\perp+\kappa_c}
$$
を得る。共振器が狭ければ $\kappa_c\ll\gamma_\perp$ なので発振は受動共振の近くにあり、原子線が狭ければ原子の中心へ寄る。これは :term[周波数引き込み]{en="frequency pulling" yomi="しゅうはすうひきこみ"} の一遷移模型である。実部はこの周波数におけるしきい反転を決める。

### 自発放出で発振位相が拡散する

しきい値より十分上の単一モードで、共振器内の平均光子数を $\bar n\gg1$ とする。光子数規格化の振幅は $\sqrt{\bar n}e^{i\phi}$。発振モードへの自然放出を、毎秒 $R_{\rm sp}$ 個の割合で大きさ 1 の独立なランダム位相の振幅を加える模型で記述する。一事象での位相の微小変化は、発振振幅に直交する成分を $\sqrt{\bar n}$ で割った値である。位相が一様なら直交成分の二乗平均は $1/2$ であり、時間差 $|\tau|$ の位相差の分散は $R_{\rm sp}|\tau|/(2\bar n)$ になる。

事象が多数ある位相拡散の近似では、平均 $\langle e^{i[\phi(t+\tau)-\phi(t)]}\rangle$ は Gaussian の特性関数より $\exp[-R_{\rm sp}|\tau|/(4\bar n)]$ となる。[コヒーレンス](/physics/optics/coherence/) の指数相関の Fourier 変換から、光の周波数全半値幅は $R_{\rm sp}/(4\pi\bar n)$ である。完全反転の理想模型では、発振モードへの誘導放出率と自然放出率の比が光子数 $\bar n$ なので、定常の光子収支から $R_{\rm sp}\simeq\tau_{\rm ph}^{-1}$ となる。この近似が :term[Schawlow–Townes 線幅]{en="Schawlow–Townes linewidth" yomi="しょうろうたうんずせんぷく"} の基本形である。振幅から位相への結合や装置の振動による周波数変動はここでは省いている。この発振線幅と、受動共振器の透過線幅 $1/(2\pi\tau_{\rm ph})$ は、異なる実験で測る量である。

### モード間で異なる利得の飽和

四準位の局所定常反転を $N(\boldsymbol r)=N_0(\boldsymbol r)/[1+\sum_j S_j F_j(\boldsymbol r)]$ とする。$S_j$ はモード $j$ のパワーを表す無次元量、$F_j\ge0$ は局所の飽和への寄与である。例えば定在波の腹で $F_j$ は大きく、節でゼロとなる。モード $i$ の利得率を、その場分布を重みにした $\int W_i N\dd^3r$ とすれば、$S_j F_j$ が小さい場合は
$$
\dot S_i=S_i\left(p_i-\sum_j b_{ij}S_j\right),\qquad
p_i=\int W_iN_0\dd^3r-\tau_i^{-1},\qquad b_{ij}=\int W_iN_0F_j\dd^3r .
$$
ここでは反転が強度変化に追従でき、モードのうなりを時間平均できる範囲を考える。モード 1 だけが発振した枝では $S_1=p_1/b_{11}>0$。モード 2 の微小な種の増加率は $p_2-b_{21}p_1/b_{11}$ である。従って小信号で $p_2>0$ であっても、既に発振したモードの飽和で成長できない場合がある。二つの定在波の腹が異なる位置にあり、モード 2 がまだ反転の残る場所を使えば、この抑制が弱くなる。

## Q スイッチとモード同期

### 反転を蓄積して一気に放出する
:term[Q スイッチ]{en="Q switching" yomi="きゅーすいっち"} では初めは共振器を高損失にして $\tau_{\rm ph}^{-1}$ を上げ、光子が増えない間に反転をためる。損失を素早く下げると $v_g\Gamma_{\rm ov}\sigma N_2-\tau_{\rm ph}^{-1}>0$ となって光子が増え、誘導放出が反転を消費する。増幅率が再び負になるとパルスは減衰する。素子の切替え時間と蓄積された反転の緩和時間が異なることを利用する。

### 縦モード間の位相をそろえる
:term[モード同期]{en="mode locking" yomi="もーどどうき"} は、等間隔 $\Delta\nu$ の縦モードに一定の相対位相を維持させる。各振幅が同じで $M$ 個を同相で足すと、搬送波を除いた振幅は $\sum_{j=0}^{M-1}e^{-2\pi i j\Delta\nu t}$。幾何級数から強度の最大は周期 $1/\Delta\nu$ ごとに繰り返し、各峰から最初の零点まで $1/(M\Delta\nu)$ である。実際の滑らかなスペクトルでは形状が変わるが、時間幅が全帯域の逆数程度になる。能動同期では往復頻度の変調、受動同期ではパルスの強い部分を通しやすい可飽和吸収などで同期状態を選ぶ。個々のモードの周波数がそろっていても、位相がランダムに変わると再現性のある短パルスは得られない。

:::example[発生間隔と幅]
$\Delta\nu=100\,\mathrm{MHz}$、$M=1000$ の理想列なら、間隔は $10\,\mathrm{ns}$、中心から第一零点まで約 $10\,\mathrm{ps}$。
:::
:::example[同じピーク強度でも異なる動作]
Q スイッチでは主に反転の蓄積と損失切替えの周期が、パルス一個のエネルギーと繰返し間隔を左右する。一方、モード同期の繰返しは共振器の往復時間で決まり、同期した帯域を増すと個々のパルスが狭くなる。
:::

### Q スイッチの反転消費
ポンプと上準位の自然減衰をパルスの短時間だけ無視する。重なりを 1、$b=v_g\sigma$、損失率を $\ell_c=1/\tau_{\rm ph}$ とすると $\dot N=-bNS$、$\dot S=(bN-\ell_c)S$ である。$S>0$ で割り算して $\dd S/\dd N=-1+N_{\rm th}/N$、$N_{\rm th}=\ell_c/b$。初期値 $N_i,S_i$ から
$$
S=S_i+N_i-N+N_{\rm th}\ln(N/N_i).
$$
$N_i>N_{\rm th}$ かつ $S_i>0$ なら、$N$ はパルス中に下がり、$\dot S$ が正から負へ変わる $N=N_{\rm th}$ でピークになる。種を $S_i=0$ と厳密に置くと、この決定論的な式からはパルスが始まらない。発振を始める瞬間の反転と、放出する途中で利得が損失に等しくなる瞬間の反転は異なる。

:::example[切替え後に生じるピーク]
Q スイッチを開いたときの反転を $N_i=2N_{\rm th}$、初期光子密度を無視できるほど小さいとすると、ピークの光子密度は $(1-\ln2)N_{\rm th}\simeq0.307N_{\rm th}$。余分な反転 $N_{\rm th}$ の全量がピークに蓄積するわけではなく、光子が成長する間にも共振器から漏れる。
:::

### 等振幅のモード同期で時間の零点を数える

角周波数が $\omega_0+j\Delta\omega$、$j=0,\ldots,M-1$ の縦モードを同じ振幅・初期位相で足す。等比和を取れば、強度をピークで割った値は $\{\sin(M\Delta\omega t/2)/[M\sin(\Delta\omega t/2)]\}^2$。繰返し周期は $T=2\pi/\Delta\omega$ で、隣の零点は $t=T/M$。一周期で積分すると、異なるモードの交差項が一周期でゼロになるので、平均強度は単独モードの $M$ 倍、ピークは $M^2$ 倍である。強度の幅を求める場合は零点間隔と半値全幅とを区別する。位相がモード番号の二次式を含むと、単純な等比和にはならずピークも変化する。

:::example[二十モードを同期]
平均パワーを固定して同期した線から周期と時間幅を求める。

入力は同期する縦モード数 $M=20$、縦モード間隔 $\Delta\nu=100000000\,\mathrm{Hz}$、時間平均パワー $P_{\rm av}=\frac{1}{100}\,\mathrm W$ である。

縦モードのうなりが一周する時刻を求める。

$$
\begin{aligned}
T&=\frac{1}{\Delta\nu},\\
T/(\mathrm s)&= \frac{1}{\left(100000000\right)}\simeq 1\times 10^{-8}.
\end{aligned}
$$

等比和の分子が最初に消える時刻を取る。

$$
\begin{aligned}
t_0&=\frac{T}{M},\\
t_0/(\mathrm s)&\simeq \frac{\left(1\times 10^{-8}\right)}{\left(20\right)}\simeq 5\times 10^{-10}.
\end{aligned}
$$

中心の両側で同じ零点が現れる。

$$
\begin{aligned}
\Delta t_{\rm zero}&=2 t_0,\\
\Delta t_{\rm zero}/(\mathrm s)&\simeq 2 \cdot \left(5\times 10^{-10}\right)\simeq 1\times 10^{-9}.
\end{aligned}
$$

同相時に振幅を足すことと周期平均の強度を比較する。

$$
\begin{aligned}
P_{\rm pk}&=M P_{\rm av},\\
P_{\rm pk}/(\mathrm W)&= \left(20\right) \cdot \left(\frac{1}{100}\right)=\frac{1}{5}.
\end{aligned}
$$

平均パワーと繰返し周期を掛ける。

$$
\begin{aligned}
E_{\rm period}&=P_{\rm av} T,\\
E_{\rm period}/(\mathrm J)&\simeq \left(\frac{1}{100}\right) \cdot \left(1\times 10^{-8}\right)\simeq 1\times 10^{-10}.
\end{aligned}
$$

繰返しはモード間隔が、零点幅は全体の帯域が決める。
:::

:::example[モード数を増やす]
同じ共振器で同期する線を五倍にする。

入力は同期する縦モード数 $M=100$、縦モード間隔 $\Delta\nu=100000000\,\mathrm{Hz}$、時間平均パワー $P_{\rm av}=\frac{1}{100}\,\mathrm W$ である。

縦モードのうなりが一周する時刻を求める。

$$
\begin{aligned}
T&=\frac{1}{\Delta\nu},\\
T/(\mathrm s)&= \frac{1}{\left(100000000\right)}\simeq 1\times 10^{-8}.
\end{aligned}
$$

等比和の分子が最初に消える時刻を取る。

$$
\begin{aligned}
t_0&=\frac{T}{M},\\
t_0/(\mathrm s)&\simeq \frac{\left(1\times 10^{-8}\right)}{\left(100\right)}\simeq 1\times 10^{-10}.
\end{aligned}
$$

中心の両側で同じ零点が現れる。

$$
\begin{aligned}
\Delta t_{\rm zero}&=2 t_0,\\
\Delta t_{\rm zero}/(\mathrm s)&\simeq 2 \cdot \left(1\times 10^{-10}\right)\simeq 2\times 10^{-10}.
\end{aligned}
$$

同相時に振幅を足すことと周期平均の強度を比較する。

$$
\begin{aligned}
P_{\rm pk}&=M P_{\rm av},\\
P_{\rm pk}/(\mathrm W)&= \left(100\right) \cdot \left(\frac{1}{100}\right)=1.
\end{aligned}
$$

平均パワーと繰返し周期を掛ける。

$$
\begin{aligned}
E_{\rm period}&=P_{\rm av} T,\\
E_{\rm period}/(\mathrm J)&\simeq \left(\frac{1}{100}\right) \cdot \left(1\times 10^{-8}\right)\simeq 1\times 10^{-10}.
\end{aligned}
$$

繰返しと一周期のエネルギーは変わらず、ピークは上がり、時間幅は狭くなる。
:::

:::example[縦間隔を拡げる]
同期本数を維持して間隔の大きい共振器を比べる。

入力は同期する縦モード数 $M=100$、縦モード間隔 $\Delta\nu=200000000\,\mathrm{Hz}$、時間平均パワー $P_{\rm av}=\frac{1}{100}\,\mathrm W$ である。

縦モードのうなりが一周する時刻を求める。

$$
\begin{aligned}
T&=\frac{1}{\Delta\nu},\\
T/(\mathrm s)&= \frac{1}{\left(200000000\right)}\simeq 5\times 10^{-9}.
\end{aligned}
$$

等比和の分子が最初に消える時刻を取る。

$$
\begin{aligned}
t_0&=\frac{T}{M},\\
t_0/(\mathrm s)&\simeq \frac{\left(5\times 10^{-9}\right)}{\left(100\right)}\simeq 5\times 10^{-11}.
\end{aligned}
$$

中心の両側で同じ零点が現れる。

$$
\begin{aligned}
\Delta t_{\rm zero}&=2 t_0,\\
\Delta t_{\rm zero}/(\mathrm s)&\simeq 2 \cdot \left(5\times 10^{-11}\right)\simeq 1\times 10^{-10}.
\end{aligned}
$$

同相時に振幅を足すことと周期平均の強度を比較する。

$$
\begin{aligned}
P_{\rm pk}&=M P_{\rm av},\\
P_{\rm pk}/(\mathrm W)&= \left(100\right) \cdot \left(\frac{1}{100}\right)=1.
\end{aligned}
$$

平均パワーと繰返し周期を掛ける。

$$
\begin{aligned}
E_{\rm period}&=P_{\rm av} T,\\
E_{\rm period}/(\mathrm J)&\simeq \left(\frac{1}{100}\right) \cdot \left(5\times 10^{-9}\right)\simeq 5\times 10^{-11}.
\end{aligned}
$$

平均パワー当たりのピークは同じで、繰返し周期が短くなる。
:::

## 代表的なレーザー

### 利得媒質から装置を比較する
レーザーの波長は、ポンプ源の波長だけでなく、利得媒質の遷移と共振器の条件で決まる。He–Ne では放電で He を励起し、衝突で Ne へ励起を渡す。よく用いられる赤色は $632.8\,\mathrm{nm}$ 付近で、気体速度による Doppler 幅を持つ。固体の Nd:YAG は Nd イオンの反転で $1064\,\mathrm{nm}$ 付近を発振し、四準位模型によって低い下準位占有を説明できる。非線形結晶で周波数を二倍にすると $532\,\mathrm{nm}$ になるが、別の発振遷移ではない。

半導体レーザーでは電子と正孔を電気的または光学的に注入し、バンド間の吸収より誘導再結合が優勢になるようにする。電子と正孔の密度、光の閉じ込め、共振器を一体に扱う。利得が周波数で広がるので単一周波数動作には共振器や格子による選択が要る。希土類添加ファイバーではポンプと信号が長い区間で重なるため、局所利得・飽和・損失を長さ方向に積分する。細い導波路では高い強度を長く維持できるが、光非線形性も増え得る。

:::example[波長から一光子のエネルギー]
$h c\simeq1240\,\mathrm{eV\,nm}$ を使うと、$1064\,\mathrm{nm}$ の光子は約 $1.17\,\mathrm{eV}$、その第二高調波は約 $2.33\,\mathrm{eV}$ である。
:::
:::example[電気的な注入]
電子一個の注入ごとに信号光子一個が外へ出る理想限界を考える。電流 $I_{\rm e}$ に対して放出個数は $I_{\rm e}/e$ 毎秒なので、出力は $(h\nu/e)I_{\rm e}$。実際の素子では励起の損失、競合再結合、内部光損失のため同じ電流でもこれより小さくなる。
:::

:::example[気体の振動遷移と半導体の帯間遷移]
同じパワーで $10.6\,\mathrm{\mu m}$ の $\mathrm{CO_2}$ レーザーと $1.06\,\mathrm{\mu m}$ の近赤外レーザーを比べると、前者は一光子のエネルギーが約十分の一なので光子束が約 10 倍である。$\mathrm{CO_2}$ では分子振動と回転の準位が発振波長群を決める。一方、半導体では電子と正孔が占有する帯の状態が遷移の分布を決める。
:::

### 半導体の準 Fermi 準位
価電子帯と伝導帯の状態の占有率を $f_v,f_c$ とする。吸収には下の電子と上の空きが必要なので $f_v(1-f_c)$、誘導放出には上の電子と下の空きが必要なので $f_c(1-f_v)$ に比例する。差は $f_c-f_v$ である。電子と正孔が同じ温度に緩和し、それぞれの準 Fermi 準位を $F_c,F_v$ で表す場合、Fermi 分布の単調性から $f_c>f_v$ は $E_c-F_c<E_v-F_v$ と同値。従ってエネルギー差 $h\nu=E_c-E_v$ の遷移での反転条件は $F_c-F_v>h\nu$ になる。これは不純物イオンの有限準位での数密度の反転に対応するが、半導体では許される波数とエネルギーの状態を積分して利得を作る。

### 励起光子とレーザー光子のエネルギー収支

光ポンプを一光子ずつ吸収し、吸収一個からレーザー光子を高々一個出す理想模型を考える。ポンプとレーザーの真空波長を $\lambda_p,\lambda_s$ とすると、吸収パワーに対する光出力の上限は $E_s/E_p=\lambda_p/\lambda_s$。これは一光子のエネルギーが $hc/\lambda$ であるためである。量子欠損のエネルギーは、熱や蛍光などレーザー出力以外に移る。ここでの上限はこの光子数の変換仮定に依存し、電気注入や一吸収から複数光子を得る過程にはそのまま適用しない。吸収での損失だけでなく、励起移送や鏡の取り出しの効率を掛けると実際の光出力が低くなる。

:::example[出力波長に近い直接励起]
ポンプと出力の波長が近い希土類レーザーの尺度を考える。

入力はポンプ波長 $\lambda_p\simeq 9.76\times 10^{-7}\,\mathrm m$、レーザー波長 $\lambda_s\simeq 1.03\times 10^{-6}\,\mathrm m$、吸収ポンプパワー $P_{\rm abs}=10\,\mathrm W$ である。

吸収する一光子のエネルギーを求める。

$$
\begin{aligned}
E_p&=\frac{hc}{\lambda_p},\\
E_p/(\mathrm J)&\simeq \frac{1986445857}{10000000000000000000000000000000000 \cdot \left(9.76\times 10^{-7}\right)}\simeq 2.0353\times 10^{-19}.
\end{aligned}
$$

同じ真空波長の換算を出力に適用する。

$$
\begin{aligned}
E_s&=\frac{hc}{\lambda_s},\\
E_s/(\mathrm J)&\simeq \frac{1986445857}{10000000000000000000000000000000000 \cdot \left(1.03\times 10^{-6}\right)}\simeq 1.9286\times 10^{-19}.
\end{aligned}
$$

出力光子一個のエネルギーを入射一個で割る。

$$
\begin{aligned}
\eta_{\max}&=\frac{E_s}{E_p},\\
\eta_{\max}&\simeq \frac{\left(1.928588211\times 10^{-19}\right)}{\left(2.035292886\times 10^{-19}\right)}=\frac{488}{515}.
\end{aligned}
$$

吸収したポンプパワーに光子エネルギーの比を掛ける。

$$
\begin{aligned}
P_{s,\max}&=P_{\rm abs} \eta_{\max},\\
P_{s,\max}/(\mathrm W)&= \left(10\right) \cdot \left(\frac{488}{515}\right)=\frac{976}{103}.
\end{aligned}
$$

光出力に回らない量子欠損を差し引きで求める。

$$
\begin{aligned}
P_{\rm rem}&=P_{\rm abs} - P_{s,\max},\\
P_{\rm rem}/(\mathrm W)&= \left(10\right) - \left(\frac{976}{103}\right)=\frac{54}{103}.
\end{aligned}
$$

一対一模型での量子欠損が小さい。
:::

:::example[より短い波長で励起]
同じ出力波長を短いポンプ波長で作る。

入力はポンプ波長 $\lambda_p\simeq 8.08\times 10^{-7}\,\mathrm m$、レーザー波長 $\lambda_s\simeq 1.064\times 10^{-6}\,\mathrm m$、吸収ポンプパワー $P_{\rm abs}=10\,\mathrm W$ である。

吸収する一光子のエネルギーを求める。

$$
\begin{aligned}
E_p&=\frac{hc}{\lambda_p},\\
E_p/(\mathrm J)&\simeq \frac{1986445857}{10000000000000000000000000000000000 \cdot \left(8.08\times 10^{-7}\right)}\simeq 2.4585\times 10^{-19}.
\end{aligned}
$$

同じ真空波長の換算を出力に適用する。

$$
\begin{aligned}
E_s&=\frac{hc}{\lambda_s},\\
E_s/(\mathrm J)&\simeq \frac{1986445857}{10000000000000000000000000000000000 \cdot \left(1.064\times 10^{-6}\right)}\simeq 1.867\times 10^{-19}.
\end{aligned}
$$

出力光子一個のエネルギーを入射一個で割る。

$$
\begin{aligned}
\eta_{\max}&=\frac{E_s}{E_p},\\
\eta_{\max}&\simeq \frac{\left(1.866960392\times 10^{-19}\right)}{\left(2.458472595\times 10^{-19}\right)}=\frac{101}{133}.
\end{aligned}
$$

吸収したポンプパワーに光子エネルギーの比を掛ける。

$$
\begin{aligned}
P_{s,\max}&=P_{\rm abs} \eta_{\max},\\
P_{s,\max}/(\mathrm W)&= \left(10\right) \cdot \left(\frac{101}{133}\right)=\frac{1010}{133}.
\end{aligned}
$$

光出力に回らない量子欠損を差し引きで求める。

$$
\begin{aligned}
P_{\rm rem}&=P_{\rm abs} - P_{s,\max},\\
P_{\rm rem}/(\mathrm W)&= \left(10\right) - \left(\frac{1010}{133}\right)=\frac{320}{133}.
\end{aligned}
$$

光子一個を作るたびに出力以外へ回るエネルギーが増える。
:::

:::example[可視ポンプから近赤外へ]
波長の離れた光励起の例を評価する。

入力はポンプ波長 $\lambda_p\simeq 5.32\times 10^{-7}\,\mathrm m$、レーザー波長 $\lambda_s\simeq 8\times 10^{-7}\,\mathrm m$、吸収ポンプパワー $P_{\rm abs}=10\,\mathrm W$ である。

吸収する一光子のエネルギーを求める。

$$
\begin{aligned}
E_p&=\frac{hc}{\lambda_p},\\
E_p/(\mathrm J)&\simeq \frac{1986445857}{10000000000000000000000000000000000 \cdot \left(5.32\times 10^{-7}\right)}\simeq 3.7339\times 10^{-19}.
\end{aligned}
$$

同じ真空波長の換算を出力に適用する。

$$
\begin{aligned}
E_s&=\frac{hc}{\lambda_s},\\
E_s/(\mathrm J)&\simeq \frac{1986445857}{10000000000000000000000000000000000 \cdot \left(8\times 10^{-7}\right)}\simeq 2.4831\times 10^{-19}.
\end{aligned}
$$

出力光子一個のエネルギーを入射一個で割る。

$$
\begin{aligned}
\eta_{\max}&=\frac{E_s}{E_p},\\
\eta_{\max}&\simeq \frac{\left(2.483057321\times 10^{-19}\right)}{\left(3.733920784\times 10^{-19}\right)}=\frac{133}{200}.
\end{aligned}
$$

吸収したポンプパワーに光子エネルギーの比を掛ける。

$$
\begin{aligned}
P_{s,\max}&=P_{\rm abs} \eta_{\max},\\
P_{s,\max}/(\mathrm W)&= \left(10\right) \cdot \left(\frac{133}{200}\right)=\frac{133}{20}.
\end{aligned}
$$

光出力に回らない量子欠損を差し引きで求める。

$$
\begin{aligned}
P_{\rm rem}&=P_{\rm abs} - P_{s,\max},\\
P_{\rm rem}/(\mathrm W)&= \left(10\right) - \left(\frac{133}{20}\right)=\frac{67}{20}.
\end{aligned}
$$

同じ吸収パワーでも量子欠損の寄与がさらに増す。
:::

## 演習問題

::::exercise{level=basic}
$d_1=d_2$ のとき $B_{12}/B_{21}$ を求めよ。

:::solution
:eqref[eq-optics-einstein] で縮退度が等しいなら両係数は等しく、比は 1 になる。
:::
::::

::::exercise{level=basic}
$N_2=3N_1,d_2=4d_1$ なら反転しているか。

:::solution
$N_2-(d_2/d_1)N_1=3N_1-4N_1=-N_1<0$。従って正味の応答は吸収である。
:::
::::

::::exercise{level=basic}
均一線の $T_2$ を半分にすると全半値幅はどうなるか。

:::solution
幅は $1/(\pi T_2)$ なので二倍になる。
:::
::::

::::exercise{level=basic}
$g=3\,\mathrm{m^{-1}},\alpha=1\,\mathrm{m^{-1}},L=0.5\,\mathrm m$ の小信号パワー利得を求めよ。

:::solution
指数が $(3-1)(0.5)=1$ なので $P_{\rm out}/P_{\rm in}=e$ である。
:::
::::

::::exercise{level=advanced}
等縮退の三準位模型で $WT_1=2$ のとき反転差を全粒子数で割った値を求めよ。

:::solution
$N_2/N=2/3$、$N_1/N=1/3$ で、差を全数で割ると $1/3$ である。
:::
::::

::::exercise{level=advanced}
$R_1=R_2=e^{-0.10},\ell=L=0.20\,\mathrm m,\alpha=0.5\,\mathrm{m^{-1}}$ のしきい値利得を求めよ。

:::solution
式 :eqref[eq-optics-laser-threshold] の第一項は $0.5\,\mathrm{m^{-1}}$、第二項は $0.20/(0.40\,\mathrm m)=0.5\,\mathrm{m^{-1}}$。和は $1.0\,\mathrm{m^{-1}}$ である。
:::
::::

::::exercise{level=advanced}
等振幅の縦モード $M=500$ を位相同期し、間隔を $\Delta\nu=200\,\mathrm{MHz}$ とする。パルス列の繰返し周期と、ピークから第一零点までの時間を求めよ。

:::solution
間隔は $1/\Delta\nu=5\,\mathrm{ns}$。零点までの時間は $1/(M\Delta\nu)=10\,\mathrm{ps}$ である。
:::
::::

::::exercise{level=advanced}
飽和増幅で $I_{\rm out}=2I_{\rm in}$ かつ $I_{\rm in}=I_s$ を得るための $g_0 L$ を求めよ。

:::solution
積分した式へ代入すると対数項が $\ln 2$、強度差を $I_s$ で割ると 1 なので $g_0L=\ln 2+1$ である。
:::
::::

::::exercise{level=advanced}
無損失の連続波増幅器で、小信号の光学的利得長 $g_0L$ を測る。入力が飽和強度の半分、出力が飽和強度の二倍なら、非飽和近似で得る値と飽和を含む値を比較せよ。

:::solution
出力と入力の比は $2/(1/2)=4$ なので、非飽和の指数則なら $g_0L=\ln4$ と見積もる。一方、$(x^{-1}+1)\dd x=g_0\dd z$ を $x=1/2$ から $2$ まで積分すると $g_0L=\ln4+2-1/2=\ln4+3/2$。観測された同じ出力を得るには、途中の飽和を補うため非飽和の値より大きい小信号利得が要る。
:::
::::

## 確認テスト

::::quiz
等縮退の正温度平衡で誘導放出と吸収を比較せよ。

- [x] 吸収の率が大きい。
- [ ] 誘導放出の率が大きい。
- [ ] 占有数は常に等しい。

:::explanation
係数は等しい一方、Boltzmann 比 $N_2/N_1=e^{-h\nu/(k_BT)}<1$ である。従って共通の放射密度を掛けた吸収率が大きい。
:::
::::

::::quiz{answer="0.25" tol="1e-6"}
単純な均一飽和模型で $I=3I_s$ の利得は小信号の何倍か。

:::explanation
反転と利得は同じ比例で抑えられ、$1/(1+3)=1/4$ 倍になる。
:::
::::

::::quiz
モード同期に直接必要なものを選べ。

- [x] 縦モード間の再現性のある位相関係
- [ ] レーザー下準位の寿命を共振器往復時間に等しくすること
- [ ] 常に一つの縦モードだけを発振させること

:::explanation
各周波数の波がそろう時刻に強いパルスとなる。複数の縦モードの相対位相がばらばらになると、この強め合いが周期ごとに再現されない。
:::
::::

::::quiz{answer="2" tol="1e-6"}
$R_1=R_2=e^{-0.20},\ell=0.10\,\mathrm m$、内部損失なしの $g_{\rm th}$ を $\mathrm{m^{-1}}$ で答えよ。

:::explanation
鏡二枚での対数損失が 0.40、往復の活性長が $0.20\,\mathrm m$。従って $g_{\rm th}=0.40/0.20=2\,\mathrm{m^{-1}}$ になる。
:::
::::

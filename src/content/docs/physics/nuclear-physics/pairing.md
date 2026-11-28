---
title: 対相関と原子核の超流動性
description: 偶奇質量差、seniority、BCS と HFB を結び、対相関による回転応答の変化を調べる。
notation:
  - symbol: '\Omega'
    meaning: 殻に含まれる時間反転した対の数
  - symbol: 'v'
    meaning: seniority
  - symbol: '\Delta'
    meaning: 対相関のギャップ
  - symbol: '\kappa'
    meaning: 異常密度
---

## この章の目標

- 質量の偶奇差から対相関の尺度を求める。
- 縮退した一殻の対相互作用を厳密に対角化する。
- BCS のギャップ方程式と準粒子エネルギーを導く。
- HFB の変分条件と回転への応答を説明する。

## 対相互作用と偶奇質量差

同じ種類の核子について、球対称な一つの $j$ 殻では時間反転で結ばれた軌道の組を角運動量ゼロに結合できる。一組だけを占める状態では射影がゼロになる。$J=0$ の対を作るには、後で示すように各組を決まった位相で足す。対を作れる核子が多い偶数系は、隣の奇数系より余分に結合する傾向がある。$Z$ を固定し、正の結合エネルギーを $B(N)$ と書く。三点差分を
$$
\Delta^{(3)}(N)=\frac{(-1)^N}{2}\big[2B(N)-B(N-1)-B(N+1)\big]
\label{eq-nuc-pairing-mass}
$$
と定める。:ref[def-nuc-binding] の規約に従う。滑らかな項 $B_{\rm sm}$ に $b(-1)^N$ を足す模型では、交代項が $2b$ を与える。滑らかな項にも曲率があるため、差分の値を準粒子のエネルギーギャップと同一視するには注意する。

::::remark
固定 $Z$ の質量エネルギーを用いる場合は、結合エネルギーと符号が逆である。核子の静止エネルギーなど $N$ の一次式は二階差分で消える。
::::

:::example[偶数核からの差分]
$B(N-1)=90,B(N)=100,B(N+1)=107\,\mathrm{MeV}$ で $N$ が偶数なら、差分は $(200-90-107)/2=1.5\,\mathrm{MeV}$ となる。
:::

:::example[背景の曲率]
$B(N)=a+bN+cN^2$ に偶奇項がなくても $\Delta^{(3)}(N)=(-1)^{N+1}c$ となる。隣接する偶数点と奇数点を比べると背景から受ける符号が異なる。
:::

### 三点・四点・五点の質量差

固定した $Z$ で結合エネルギーの三点差を $\Delta^{(3)}(N)=(-1)^N[2B(N)-B(N-1)-B(N+1)]/2$ とする。この量に隣接する $\Delta^{(3)}(N+1)$ の平均を取ると四点差になる。交代符号によって一定の曲率が消えるが、評価位置は $N+1/2$ である。前後の四点差をさらに平均すると、$N$ を中心とする五点差
$$
\Delta^{(5)}(N)=\frac{(-1)^N}{8}
[6B(N)-4B(N-1)-4B(N+1)+B(N-2)+B(N+2)]\label{eq-nuc-pairing-five-point}
$$
を得る。角括弧は四次の中心差分であり、定数から三次までの多項式でゼロになる。$B=b(-1)^N$ なら角括弧は $16b(-1)^N$ で $\Delta^{(5)}=2b$ になる。殻閉じの近くでは、滑らかな多項式という前提自体を検討する必要がある。

時間反転対の仮定は、二粒子を同じ場所に局在させるという意味ではない。球対称なら二つの角運動量を $J=0$ に結合して回転不変な対を作る。各 $m$ 成分に $(-1)^{j-m}$ を付けると、上げ演算子を一方に作用させた項が、隣接する $m$ で他方に作用させた項と相殺する。射影はゼロなのでこの和は $J=0$ である。奇核では一つの対軌道が単独に占有され、この Pauli 禁制が偶奇差を作る機構になる。

:::example[偶数と奇数の推定値]
背景を $B_{\rm sm}=a+cN^2$、交代項を $b(-1)^N$ とおく。$p=(-1)^N$ を使って三点差を得る。

量の定義に従い、共通の単位で次の値を用いる。
$$
b=\frac{3}{5}\,\mathrm{MeV},\qquad c=- \frac{1}{5}\,\mathrm{MeV},\qquad p=1
$$

この条件で求める量を ${\Delta^{(3)}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\Delta^{(3)}}&=2 \cdot b - c \cdot p\\&=- \left(- \frac{1}{5}\right) \cdot \left(1\right) + 2 \cdot \left(\frac{3}{5}\right)\\&\simeq 1.4\,\mathrm{MeV}.\end{aligned}
$$

隣接平均はこの模型で $2b$ に戻る。

奇数側を中心にした結果と比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
p=-1
$$

このときの値を ${\Delta^{(3)}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\Delta^{(3)}}_{\mathrm{new}}&=- \left(- \frac{1}{5}\right) \cdot \left(-1\right) + 2 \cdot \left(\frac{3}{5}\right)\\&=1\,\mathrm{MeV},\\{\Delta^{(3)}}_{\mathrm{new}}-{\Delta^{(3)}}&\simeq -0.4\,\mathrm{MeV}.\end{aligned}
$$

**一次近似の検証。** 交代符号を実数の係数としたときの傾きが、偶数・奇数の端点差と整合するか調べる。

${p}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${\Delta^{(3)}}$ の一階微分を求めると、
$$
\frac{\mathrm d{\Delta^{(3)}}}{\mathrm dp}=- c.
$$

この導関数を変更前の点で評価して $h=-2$ を掛ける。ここで $p$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq -0.4\,\mathrm{MeV},\\\Delta-\Delta_{\mathrm{lin}}&=0\,\mathrm{MeV}.\end{aligned}
$$

入力のこの変化について差は一次の項だけで得られる。従って基準点での傾きで増分を求められる。
:::

:::example[五点差の残留曲率]
滑らかな背景に $c_4(N-N_0)^4$ があるとする。:eqref[eq-nuc-pairing-five-point] で中心 $N_0$ の五点差に加わる分は $24p c_4/8$ である。

量の定義に従い、共通の単位で次の値を用いる。
$$
p=1,\qquad c_{4}=\frac{1}{100}\,\mathrm{MeV}
$$

この条件で求める量を ${\delta\Delta}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\delta\Delta}&=3 \cdot c_{4} \cdot p\\&=3 \cdot \left(1\right) \cdot \left(\frac{1}{100}\right)\\&\simeq 0.03\,\mathrm{MeV}.\end{aligned}
$$

四次の背景は消えないため、中心が強く非放物線的な場合は結果の解釈に注意する。

四次成分の大きい模型に変える。そこで次の値に変更し、ほかの入力を固定する。
$$
c_{4}=\frac{1}{25}\,\mathrm{MeV}
$$

このときの値を ${\delta\Delta}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\delta\Delta}_{\mathrm{new}}&=3 \cdot \left(1\right) \cdot \left(\frac{1}{25}\right)\\&\simeq 0.12\,\mathrm{MeV},\\{\delta\Delta}_{\mathrm{new}}-{\delta\Delta}&\simeq 0.09\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[対のブロッキング]
縮退した $\Omega$ 個の対のうち $b$ 個を単独粒子が占める。一対だけ残りの軌道に作ると、完全対称な和の相関エネルギーは $-G(\Omega-b)$ である。

量の定義に従い、共通の単位で次の値を用いる。
$$
G=\frac{2}{5}\,\mathrm{MeV},\qquad \Omega=6,\qquad b=0
$$

この条件で求める量を ${E}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{E}&=- G \cdot \left(\Omega - b\right)\\&=- \left(\frac{2}{5}\right) \cdot \left(- \left(0\right) + \left(6\right)\right)\\&\simeq -2.4\,\mathrm{MeV}.\end{aligned}
$$

単独粒子数の違いによる他の一体エネルギーはここに含めていない。

一軌道をブロックする。そこで次の値に変更し、ほかの入力を固定する。
$$
b=1
$$

このときの値を ${E}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{E}_{\mathrm{new}}&=- \left(\frac{2}{5}\right) \cdot \left(- \left(1\right) + \left(6\right)\right)\\&=-2\,\mathrm{MeV},\\{E}_{\mathrm{new}}-{E}&\simeq 0.4\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[零温度の奇偶指標]
偶数系のエネルギーを $E_0$ とし、中央の奇数系のエネルギーが左右の偶数の平均から $d$ 上にあるとする。質量差から取り出せるのはこの超過である。

量の定義に従い、共通の単位で次の値を用いる。
$$
E_{o}=-49\,\mathrm{MeV},\qquad E_{l}=-45\,\mathrm{MeV},\qquad E_{r}=-55\,\mathrm{MeV}
$$

この条件で求める量を ${d}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{d}&=E_{o} - \frac{E_{l} + E_{r}}{2}\\&=- \frac{\left(-45\right)}{2} + \left(-49\right) - \frac{\left(-55\right)}{2}\\&=1\,\mathrm{MeV}.\end{aligned}
$$

逆符号の結合エネルギーで計算する場合は偶奇の符号も整合させる。

奇数系の相関がさらに失われた場合を比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
E_{o}=-48\,\mathrm{MeV}
$$

このときの値を ${d}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{d}_{\mathrm{new}}&=- \frac{\left(-45\right)}{2} + \left(-48\right) - \frac{\left(-55\right)}{2}\\&=2\,\mathrm{MeV},\\{d}_{\mathrm{new}}-{d}&=1\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[殻境界と対の切り分け]
第一近似で偶数中心の三点差がギャップ $d$ の半分と対指標 $p$ の和、奇数中心では同じ対指標 $p$ だけであると仮定する。二つの差から $d$ を推定する。

量の定義に従い、共通の単位で次の値を用いる。
$$
D_{e}=2\,\mathrm{MeV},\qquad D_{o}=\frac{7}{10}\,\mathrm{MeV}
$$

この条件で求める量を ${d}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{d}&=2 \cdot \left(D_{e} - D_{o}\right)\\&=2 \cdot \left(2\right) - 2 \cdot \left(\frac{7}{10}\right)\\&\simeq 2.6\,\mathrm{MeV}.\end{aligned}
$$

二つの中心で同じ対指標を仮定しているため、対強度が変わる系列ではこの逆算に模型誤差がある。

奇数中心の値を変える。そこで次の値に変更し、ほかの入力を固定する。
$$
D_{o}=\frac{11}{10}\,\mathrm{MeV}
$$

このときの値を ${d}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{d}_{\mathrm{new}}&=2 \cdot \left(2\right) - 2 \cdot \left(\frac{11}{10}\right)\\&\simeq 1.8\,\mathrm{MeV},\\{d}_{\mathrm{new}}-{d}&\simeq -0.8\,\mathrm{MeV}.\end{aligned}
$$
:::

## seniority 模型

### 準スピン

単一の縮退した $j$ 殻の $2\Omega=2j+1$ 状態を、時間反転の対 $(k,\bar k)$ に分ける。位相を対の作成演算子に吸収して
$$
P_k^\dagger=a_k^\dagger a_{\bar k}^\dagger,\qquad
S_+=\sum_k P_k^\dagger,\quad S_-=S_+^\dagger,\quad S_z=(\hat n-\Omega)/2
$$
と置く。$\hat n$ は殻内粒子数演算子である。異なる対の演算子は互いに可換である。同じ対について $P_k^\dagger P_k=n_kn_{\bar k}$ と $P_kP_k^\dagger=(1-n_k)(1-n_{\bar k})$ が成り立つので $[S_+,S_-]=2S_z$ となる。また $[\hat n,S_\pm]=\pm2S_\pm$ である。これらは角運動量と同じ無次元の代数をなし、通常の角運動量と区別して準スピンと呼ぶ。

:::definition[seniority]{#def-nuc-seniority}
準スピン $s$ の既約表現の最小粒子数を :term[seniority]{en="seniority" yomi="せにおりてぃ"} と呼び、$v=\Omega-2s$ と書く。$S_-$ で消える最下位状態の粒子数であり、そこから $S_+$ を作用させて対を追加しても $v$ は変わらない。
:::

### 固有値

一定の強さ $G>0$ の模型を $H=\epsilon\hat n-GS_+S_-$ とする。$\epsilon$ は共通の軌道エネルギーである。$S_+S_-=S^2-S_z(S_z-1)$ を使うと、粒子数 $n$ における固有値は
$$
E(n,v)=\epsilon n-\frac G4(n-v)(2\Omega-n-v+2).
\label{eq-nuc-seniority-energy}
$$
ここで $n-v$ は偶数であり、$v\le n\le 2\Omega-v$ を満たす。差 $E(n,v+2)-E(n,v)=G(\Omega-v)$ は許された範囲で正なので、同じ $n$ なら最小の $v$ が最も結合する。一殻の縮退と一定の対散乱強度を外すと、この固有値は一般に成り立たない。

:::example[対を一つ切る]
$\Omega=4,n=4$ とする。$\epsilon n$ を引いた値は $v=0$ で $-6G$、$v=2$ で $-2G$ である。対を一つ切るエネルギーは $4G$ となる。
:::

:::example[粒子と穴]
$n=2$ と $n=2\Omega-2$ ではどちらも $v=0$ から $v=2$ への差は $G\Omega$ である。対エネルギーの絶対値は異なるが、同じ粒子数の中で比較した励起の差は等しい。
:::

### 厳密な偶奇差とギャップの比較

縮退殻で $\epsilon\hat n$ は二階差分から消える。$n=2p$ の $v=0$ の対エネルギーは $-Gp(\Omega-p+1)$ である。前後の奇数核の最小 $v=1$ ではそれぞれ $-G(p-1)(\Omega-p+1)$ と $-Gp(\Omega-p)$ になる。従って偶数中心の差 $[E(n-1)+E(n+1)-2E(n)]/2$ は $G(\Omega+1)/2$ である。奇数中心 $n=2p+1$ については次の偶数核の値が $-G(p+1)(\Omega-p)$ であることを使う。ここで差に $(-1)^n=-1$ を掛けると $G\Omega/2$ になる。全ての値が殻内にある場合に成り立ち、殻の端を跨ぐ差分には別の軌道が必要である。

:::example[縮退殻での差分]
$\Omega=4,G=0.5\,\mathrm{MeV}$ の殻内では偶数中心が $1.25\,\mathrm{MeV}$、奇数中心が $1.0\,\mathrm{MeV}$ となる。半充填で得られる BCS の $\Delta=G\Omega/2$ は奇数中心の値と一致する。偶数中心には一定の相互作用のもとでも $G/2$ のずれが残る。
:::

### 準スピンの表現と励起の間隔

一対の四状態を空、二重占有、単独占有二つに分ける。空と二重占有には $S_z=-1/2,+1/2$ を割り当て、生成と消滅を上下の演算子とする。単独占有では上下の作用がゼロで $S_z=0$ である。異なる対軌道の演算子は二つのフェルミ演算子の積なので、互いの順序を入れ替えると符号が四回変わり可換である。よって殻全体の和が同じ交換関係を満たす。固定粒子数の準スピン射影を $m=(n-\Omega)/2$ とすれば、$v=\Omega-2s$ は $0\le v\le\min(n,2\Omega-n)$ かつ $v\equiv n\pmod2$ を満たす。特定の $j$ でどの角運動量が存在するかは、この数の条件のほかに反対称化を要する。

:::example[一対の集団固有値]
対が一つで全軌道が等しい。行列要素がすべて $-G$ の $\Omega\times\Omega$ 行列を等成分のベクトルに作用させる。

量の定義に従い、共通の単位で次の値を用いる。
$$
G=\frac{3}{10}\,\mathrm{MeV},\qquad \Omega=3
$$

この条件で求める量を ${E_0}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{E_0}&=- G \cdot \Omega\\&=- \left(3\right) \cdot \left(\frac{3}{10}\right)\\&\simeq -0.9\,\mathrm{MeV}.\end{aligned}
$$

等成分に直交するベクトルは成分和がゼロで固有値もゼロである。

対軌道を増やす。そこで次の値に変更し、ほかの入力を固定する。
$$
\Omega=6
$$

このときの値を ${E_0}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{E_0}_{\mathrm{new}}&=- \left(6\right) \cdot \left(\frac{3}{10}\right)\\&\simeq -1.8\,\mathrm{MeV},\\{E_0}_{\mathrm{new}}-{E_0}&\simeq -0.9\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[二個の対を持つ殻]
一定対作用の式で $n=4,v=0$ を代入する。

量の定義に従い、共通の単位で次の値を用いる。
$$
G=\frac{2}{5}\,\mathrm{MeV},\qquad \Omega=4
$$

この条件で求める量を ${E}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{E}&=- 2 \cdot G \cdot \left(\Omega - 1\right)\\&=- 2 \cdot \left(\frac{2}{5}\right) \cdot \left(\left(4\right) - 1\right)\\&\simeq -2.4\,\mathrm{MeV}.\end{aligned}
$$

独立な一対のエネルギーの二倍との違いは Pauli 原理による。

殻の容量を増やす。そこで次の値に変更し、ほかの入力を固定する。
$$
\Omega=6
$$

このときの値を ${E}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{E}_{\mathrm{new}}&=- 2 \cdot \left(\frac{2}{5}\right) \cdot \left(\left(6\right) - 1\right)\\&=-4\,\mathrm{MeV},\\{E}_{\mathrm{new}}-{E}&\simeq -1.6\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[最初に対を破るエネルギー]
同じ偶数 $n$ で $v=0,2$ が許される場合を考える。固有値を引くとギャップは $G\Omega$ となる。

量の定義に従い、共通の単位で次の値を用いる。
$$
G=\frac{2}{5}\,\mathrm{MeV},\qquad \Omega=5
$$

この条件で求める量を ${D}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{D}&=G \cdot \Omega\\&=\left(5\right) \cdot \left(\frac{2}{5}\right)\\&=2\,\mathrm{MeV}.\end{aligned}
$$

粒子数の項が差から消える。同じ殻内での一定対作用に固有の性質である。

結合定数を変更する。そこで次の値に変更し、ほかの入力を固定する。
$$
G=\frac{1}{2}\,\mathrm{MeV}
$$

このときの値を ${D}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{D}_{\mathrm{new}}&=\left(5\right) \cdot \left(\frac{1}{2}\right)\\&\simeq 2.5\,\mathrm{MeV},\\{D}_{\mathrm{new}}-{D}&\simeq 0.5\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[対付加の行列要素]
準スピンの梯子で $S_+$ のノルム二乗は $(s-m)(s+m+1)$ である。

量の定義に従い、共通の単位で次の値を用いる。
$$
s=3,\qquad m=-2
$$

この条件で求める量を ${B_+}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{B_+}&=\left(- m + s\right) \cdot \left(m + s + 1\right)\\&=\left(- \left(-2\right) + \left(3\right)\right) \cdot \left(\left(-2\right) + \left(3\right) + 1\right)\\&=10.\end{aligned}
$$

最高射影 $m=s$ では Pauli 原理で追加できない。

殻の下部と上部を比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
m=2
$$

このときの値を ${B_+}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{B_+}_{\mathrm{new}}&=\left(- \left(2\right) + \left(3\right)\right) \cdot \left(\left(2\right) + \left(3\right) + 1\right)\\&=6,\\{B_+}_{\mathrm{new}}-{B_+}&=-4.\end{aligned}
$$

**一次近似の検証。** 梯子の端の非対称を、射影の周りで一次まで残した場合と比較する。微分は許された射影を結ぶ式の数学的な補間に対して取る。

${m}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${B_+}$ の一階微分を求めると、
$$
\frac{\mathrm d{B_+}}{\mathrm dm}=- 2 \cdot m - 1.
$$

この導関数を変更前の点で評価して $h=4$ を掛ける。ここで $m$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&=12,\\\Delta-\Delta_{\mathrm{lin}}&=-16.\end{aligned}
$$

省いた分が全変化に比べて大きい。従って増分の大きさを予測するには、少なくとも二次の項を調べるか、元の式を使う必要がある。
:::

:::example[粒子・穴でのエネルギーのずれ]
同じ $v$ について $m\mapsto-m$ で $S_+S_-$ の固有値差を取ると、共通の一体エネルギー $\epsilon n$ を除いた対エネルギーの差は $E_{\rm pair}(m)-E_{\rm pair}(-m)=-2Gm$ である。

量の定義に従い、共通の単位で次の値を用いる。
$$
G=\frac{7}{20}\,\mathrm{MeV},\qquad m=1
$$

この条件で求める量を ${\Delta E}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\Delta E}&=- 2 \cdot G \cdot m\\&=- 2 \cdot \left(1\right) \cdot \left(\frac{7}{20}\right)\\&\simeq -0.7\,\mathrm{MeV}.\end{aligned}
$$

同じ $m$ の中の励起エネルギーを比較するとこの線形差は相殺される。

半充填から離す。そこで次の値に変更し、ほかの入力を固定する。
$$
m=2
$$

このときの値を ${\Delta E}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\Delta E}_{\mathrm{new}}&=- 2 \cdot \left(2\right) \cdot \left(\frac{7}{20}\right)\\&\simeq -1.4\,\mathrm{MeV},\\{\Delta E}_{\mathrm{new}}-{\Delta E}&\simeq -0.7\,\mathrm{MeV}.\end{aligned}
$$
:::

## BCS 理論と準粒子

### 状態とギャップ方程式

異なるエネルギー $\epsilon_k$ の軌道を含める。独立に対を空または占有にできる規格化状態は
$$
\ket{\Phi}=\prod_k(u_k+v_kP_k^\dagger)\ket0,\qquad u_k^2+v_k^2=1
$$
である。実係数を選び、単粒子場に吸収した正常項を除いて、対エネルギーを $-G(\sum_k u_kv_k)^2$ とする平均場模型を採る。このとき制約 $\langle\hat n\rangle=2\sum_kv_k^2$ を化学ポテンシャル $\lambda$ で実装する。$\xi_k=\epsilon_k-\lambda$、$\Delta=G\sum_k u_kv_k$ と置き、$u_k=\cos\theta_k,v_k=\sin\theta_k$ として変分する。$2\xi_k\sin2\theta_k-2\Delta\cos2\theta_k=0$ なので、エネルギーを低くする枝で $\cos2\theta_k=\xi_k/E_k$、$\sin2\theta_k=\Delta/E_k$ となる。従って
$$
E_k=\sqrt{\xi_k^2+\Delta^2},\qquad
v_k^2=\tfrac12(1-\xi_k/E_k),\qquad
1=G\sum_k(2E_k)^{-1}
\label{eq-nuc-bcs-gap}
$$
である。最後の式は $\Delta\ne0$ のときに限る。有限系では準位間隔が大きいと非ゼロ解を失う。粒子数を厳密に固定した状態でも、対の散乱による相関は残り得る。

$\alpha_k=u_ka_k-v_ka_{\bar k}^\dagger$、$\alpha_{\bar k}=u_ka_{\bar k}+v_ka_k^\dagger$ は正準反交換関係を満たし、$\ket\Phi$ を消す。この :term[準粒子]{en="quasiparticle" yomi="じゅんりゅうし"} は粒子と穴を混ぜる。二行二列の平均場行列 $\begin{pmatrix}\xi_k&-\Delta\\-\Delta&-\xi_k\end{pmatrix}$ の固有値が $\pm E_k$ であるため、励起の尺度は $E_k$ となる。奇数核では一方の時間反転軌道が既に一粒子を含み、対を作れない。この軌道を対の和から除く操作をブロッキングという。

:::example[Fermi 面の占有]
$\xi_k=0,\Delta>0$ なら $v_k^2=1/2$、$E_k=\Delta$ である。$\xi_k=3\Delta/4$ では $E_k=5\Delta/4$、$v_k^2=1/5$ となる。
:::

:::example[縮退模型の半充填]
$\Omega$ 対のすべてで $\xi_k=0$ とする。半充填で $u_kv_k=1/2$ だから $\Delta=G\Omega/2$ である。二準粒子のエネルギー $2\Delta=G\Omega$ は、:eqref[eq-nuc-seniority-energy] の $v=0$ と $v=2$ の差と一致する。
:::

### 粒子数の射影

BCS の積状態は粒子数固有状態の重ね合わせである。規格化のもとで粒子数 $n$ を射影する演算子を
$$
\mathcal P_n=\frac1{2\pi}\int_0^{2\pi}e^{i\varphi(\hat n-n)}\dd\varphi
$$
と置く。整数固有値 $m$ の成分では積分が $\delta_{mn}$ になるので、$\mathcal P_n^2=\mathcal P_n^\dagger=\mathcal P_n$ である。縮退したすべての対で同じ $u,v$ を選ぶなら、$2p$ 粒子への射影は空でない場合 $(S_+)^p\ket0$ に比例する。実際、積を展開したとき各 $p$ 対の係数が $u^{\Omega-p}v^p$ で共通であり、$(S_+)^p$ では同じ選択を $p!$ 回数える。従って射影状態は準スピンが最大の $v=0$ の厳密解になる。この事実は縮退と一定の対相互作用に依存する。

:::example[二対軌道の半充填]
$\Omega=2,u=v=1/\sqrt2$ では $n=0,2,4$ の確率は $1/4,1/2,1/4$ である。$n=2$ を射影して規格化すると $(P_1^\dagger+P_2^\dagger)\ket0/\sqrt2$ になり、:eqref[eq-nuc-seniority-energy] で対エネルギーは $-2G$ となる。
:::

### 対 Hamiltonian の期待値と正常項

一定強度の対相互作用 $-G\sum_{kl}P_k^\dagger P_l$ に積状態を代入する。$k\ne l$ の積の期待値は $\langle P_k^\dagger\rangle\langle P_l\rangle=u_kv_ku_lv_l$ である。$k=l$ の場合は $P_k^\dagger P_k$ が占有対への射影なので期待値は $v_k^2$ である。すべてを $(\sum_k u_kv_k)^2$ にまとめると対角項を $u_k^2v_k^2$ として数えているため、実際には差 $v_k^2-u_k^2v_k^2=v_k^4$ が残る。従って全期待値は
$$
E=2\sum_k\epsilon_k v_k^2-G\left(\sum_k u_kv_k\right)^2-G\sum_k v_k^4
$$
である。対角補正も変分する場合、$v_k^2$ の微分は $-2Gv_k^2$ を与え、ギャップ式の $\xi_k$ は $\epsilon_k-\lambda-Gv_k^2$ に替わる。本節前半の式は、この正常場の補正を単粒子場に含めて固定する形で導いた。

:::example[一つの軌道しかない場合]
$\Omega=1$ では軌道間の対散乱はない。それでも $u_k^2v_k^2$ と $v_k^4$ を足せば $v_k^2$ となる。相互作用の期待値は $-Gv_k^2$ であり、対が存在する確率にそのエネルギーを掛けた値である。
:::

### BCS の停留条件と数ゆらぎ

各対の空状態と二重占有の混合角を $u_k=\cos(\theta_k/2)$、$v_k=\sin(\theta_k/2)$ と置く。正常場へ吸収した対角寄与を除いた模型の粒子数制約付き期待値は $\sum_k \xi_k(1-\cos\theta_k)-G(\sum_k\sin\theta_k/2)^2$ となる。角度の微分は $\xi_k\sin\theta_k-\Delta\cos\theta_k$、$\Delta=G\sum_k\sin\theta_k/2$ である。位相の揃った安定な枝で $\cos\theta_k=\xi_k/E_k$、$\sin\theta_k=\Delta/E_k$、$E_k=\sqrt{\xi_k^2+\Delta^2}$ となる。$u_kv_k=\Delta/(2E_k)$ を和に戻すと解の自己無撞着条件になる。$\Delta=0$ も元の式の解であるので、$\Delta$ で割った後にこの枝を落とさない。

一定状態密度 $g$ の対軌道が $-\omega_c$ から $\omega_c$ まであるとすると、非零解は $1=Gg\int_0^{\omega_c}\dd\xi/\sqrt{\xi^2+\Delta^2}$ に従う。$\xi=\Delta\sinh y$ と置いて積分すれば $\Delta=\omega_c/\sinh[1/(Gg)]$ であり、弱結合の式 $2\omega_c e^{-1/(Gg)}$ を得る。連続近似で任意の微小引力に解が出るのは $\xi=0$ に近い状態を取り続けられるためで、離散殻の臨界結合とは仮定が違う。

:::example[離散系で平均場のギャップが消える場合]
単粒子エネルギーが $\pm d$、$d>0$ の二つの対軌道を二粒子で占める。正常項を固定した BCS では化学ポテンシャルを中点に置けるので、非零ギャップの条件は $1=G/\sqrt{d^2+\Delta^2}$ となり、$G\le d$ ではその解がない。一方、粒子数を固定し $-G\sum_{kl}P_k^\dagger P_l$ を残した Hamiltonian は、時間反転対を一つ占めた基底で
$$
\begin{pmatrix}-2d-G&-G\\-G&2d-G\end{pmatrix}
$$
である。下側の固有値 $-G-\sqrt{4d^2+G^2}$ は、低い軌道だけを占めた値 $-2d-G$ より任意の $G>0$ で低い。中間の差は弱結合で $-G^2/(4d)+O(G^4/d^3)$ となる。粒子数を固定しても、同じ数を持つ二配位の間で対が散乱するためである。
:::

:::example[占有のにじみ]
一軌道の $\xi,\Delta$ を使って $v^2$ を求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
\xi=\frac{3}{2}\,\mathrm{MeV},\qquad \Delta=1\,\mathrm{MeV}
$$

この条件で求める量を ${v^2}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{v^2}&=\frac{1 - \frac{\xi}{\sqrt{\Delta^{2} + \xi^{2}}}}{2}\\&=- \frac{\left(\frac{3}{2}\right)}{2 \cdot \sqrt{\left(1\right)^{2} + \left(\frac{3}{2}\right)^{2}}} + \frac{1}{2}\\&\simeq 0.0839749.\end{aligned}
$$

Fermi 面をはさんで占有が余事象になる。

対称なエネルギーの軌道と比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
\xi=- \frac{3}{2}\,\mathrm{MeV}
$$

このときの値を ${v^2}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{v^2}_{\mathrm{new}}&=- \frac{\left(- \frac{3}{2}\right)}{2 \cdot \sqrt{\left(- \frac{3}{2}\right)^{2} + \left(1\right)^{2}}} + \frac{1}{2}\\&\simeq 0.916025,\\{v^2}_{\mathrm{new}}-{v^2}&\simeq 0.83205.\end{aligned}
$$
:::

:::example[対振幅の距離依存]
同じ $\Delta$ で、粒子数を変える対の振幅 $uv$ を計算する。

量の定義に従い、共通の単位で次の値を用いる。
$$
\xi=0\,\mathrm{MeV},\qquad \Delta=1\,\mathrm{MeV}
$$

この条件で求める量を ${uv}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{uv}&=\frac{\Delta}{2 \cdot \sqrt{\Delta^{2} + \xi^{2}}}\\&=\frac{\left(1\right)}{2 \cdot \sqrt{\left(0\right)^{2} + \left(1\right)^{2}}}\\&\simeq 0.5.\end{aligned}
$$

深い軌道も高い軌道も $\abs{\xi}$ が大きいと対振幅は小さくなる。

Fermi 面から離れる。そこで次の値に変更し、ほかの入力を固定する。
$$
\xi=3\,\mathrm{MeV}
$$

このときの値を ${uv}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{uv}_{\mathrm{new}}&=\frac{\left(1\right)}{2 \cdot \sqrt{\left(1\right)^{2} + \left(3\right)^{2}}}\\&\simeq 0.158114,\\{uv}_{\mathrm{new}}-{uv}&\simeq -0.341886.\end{aligned}
$$
:::

:::example[有限のカットオフでの連続解]
無次元結合 $g_0=Gg$ とエネルギー幅 $\omega_c$ からギャップを得る。

量の定義に従い、共通の単位で次の値を用いる。
$$
w=5\,\mathrm{MeV},\qquad g_{0}=\frac{1}{4}
$$

この条件で求める量を ${\Delta}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\Delta}&=\frac{w}{\sinh{\left(1 \cdot \frac{1}{g_{0}} \right)}}\\&=\frac{\left(5\right)}{\sinh{\left(\frac{1}{\left(\frac{1}{4}\right)} \right)}}\\&\simeq 0.183218\,\mathrm{MeV}.\end{aligned}
$$

弱結合では、結合の小さな相対変化が大きなギャップ変化になる。

結合を増す。そこで次の値に変更し、ほかの入力を固定する。
$$
g_{0}=\frac{2}{5}
$$

このときの値を ${\Delta}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\Delta}_{\mathrm{new}}&=\frac{\left(5\right)}{\sinh{\left(\frac{1}{\left(\frac{2}{5}\right)} \right)}}\\&\simeq 0.826418\,\mathrm{MeV},\\{\Delta}_{\mathrm{new}}-{\Delta}&\simeq 0.6432\,\mathrm{MeV}.\end{aligned}
$$

**一次近似の検証。** 弱結合での見積もりでは指数の感度が高い。有限の結合差でも局所的な傾きで近似できるか調べる。

${g_{0}}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${\Delta}$ の一階微分を求めると、
$$
\frac{\mathrm d{\Delta}}{\mathrm dg_{0}}=\frac{w \cdot \cosh{\left(\frac{1}{g_{0}} \right)}}{g_{0}^{2} \cdot \sinh^{2}{\left(\frac{1}{g_{0}} \right)}}.
$$

この導関数を変更前の点で評価して $h=\frac{3}{20}$ を掛ける。ここで $g_{0}$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq 0.440018\,\mathrm{MeV},\\\Delta-\Delta_{\mathrm{lin}}&\simeq 0.203183\,\mathrm{MeV}.\end{aligned}
$$

残差が示すのは、指定した二点間で生じる線形化の誤差である。基準点からさらに離れた値について同じ精度を保証するものではない。
:::

:::example[二対軌道の粒子数分散]
対称な $\xi=\pm d$ にある二軌道に同じ $\Delta$ がある。独立な分散 $4u^2v^2$ を二つ足す。

量の定義に従い、共通の単位で次の値を用いる。
$$
\Delta=1\,\mathrm{MeV},\qquad d=\frac{1}{2}\,\mathrm{MeV}
$$

この条件で求める量を ${\sigma_N^2}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\sigma_N^2}&=\frac{2 \cdot \Delta^{2}}{\Delta^{2} + d^{2}}\\&=\frac{2 \cdot \left(1\right)^{2}}{\left(1\right)^{2} + \left(\frac{1}{2}\right)^{2}}\\&\simeq 1.6.\end{aligned}
$$

粒子数固有状態を得るには射影が必要である。

準位差を広げる。そこで次の値に変更し、ほかの入力を固定する。
$$
d=2\,\mathrm{MeV}
$$

このときの値を ${\sigma_N^2}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\sigma_N^2}_{\mathrm{new}}&=\frac{2 \cdot \left(1\right)^{2}}{\left(1\right)^{2} + \left(2\right)^{2}}\\&\simeq 0.4,\\{\sigma_N^2}_{\mathrm{new}}-{\sigma_N^2}&\simeq -1.2.\end{aligned}
$$
:::

:::example[臨界結合近くの離散解]
対称な二軌道で $G>d$ を満たす場合を使う。

量の定義に従い、共通の単位で次の値を用いる。
$$
G=\frac{11}{10}\,\mathrm{MeV},\qquad d=1\,\mathrm{MeV}
$$

この条件で求める量を ${\Delta}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\Delta}&=\sqrt{G^{2} - d^{2}}\\&=\sqrt{- \left(1\right)^{2} + \left(\frac{11}{10}\right)^{2}}\\&\simeq 0.458258\,\mathrm{MeV}.\end{aligned}
$$

$G=d$ の極限で連続的に正常枝へ戻る。

しきい値から離す。そこで次の値に変更し、ほかの入力を固定する。
$$
G=\frac{8}{5}\,\mathrm{MeV}
$$

このときの値を ${\Delta}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\Delta}_{\mathrm{new}}&=\sqrt{- \left(1\right)^{2} + \left(\frac{8}{5}\right)^{2}}\\&\simeq 1.249\,\mathrm{MeV},\\{\Delta}_{\mathrm{new}}-{\Delta}&\simeq 0.790742\,\mathrm{MeV}.\end{aligned}
$$
:::

### 熱励起で変わるギャップ方程式

独立な正エネルギー準粒子の熱平均占有を $f_k=(1+e^{E_k/(k_{\rm B}T)})^{-1}$ とする。元の対の異常平均に Bogoliubov 変換を入れると、準粒子の非対角平均はゼロであり、対の二つの占有を引いた $u_kv_k(1-2f_k)$ が残る。数密度の方は $v_k^2(1-f_k)+u_k^2f_k$ である。従って同じ $\xi_k,\Delta$ を使っても熱励起が対の異常平均を減らす。非零ギャップでの自己無撞着条件は
$$
1=G\sum_k\frac{\tanh(E_k/(2k_{\rm B}T))}{2E_k}.
$$
有限の殻では化学ポテンシャルも粒子数条件と同時に調整する。$\Delta$ がゼロに近づく点を調べるには右辺で $E_k\to\abs{\xi_k}$ とし、$\xi_k=0$ の項は $\tanh x/x\to1$ を使う。この極限は有限である。平均場の解が急に変わる温度を、そのまま有限粒子系の厳密な相転移と解釈する必要はない。有限の模型空間の分配関数は有限個の正の指数の和であり、$T>0$ で解析的だからである。

## Hartree–Fock–Bogoliubov 理論の概観

BCS では通常場の軌道を選び、その時間反転の対だけを変分した。:term[Hartree–Fock–Bogoliubov 理論]{en="Hartree–Fock–Bogoliubov theory" yomi="はーとりーふぉっくぼごりゅーぼふりろん"} では、軌道と対を同時に変分する。密度 $\rho_{ij}=\langle a_j^\dagger a_i\rangle$ と異常密度 $\kappa_{ij}=\langle a_j a_i\rangle$ は $\rho^\dagger=\rho$、$\kappa^{\mathsf T}=-\kappa$ を満たす。純粋な準粒子真空では
$$
\mathcal R=\begin{pmatrix}\rho&\kappa\\-\kappa^*&1-\rho^*\end{pmatrix},
\qquad \mathcal R^2=\mathcal R.
$$
この冪等性は、準粒子基底で一般化密度が占有または空の射影であることから従う。例えば上左ブロックから $\rho-\rho^2=\kappa\kappa^\dagger$ を得る。

エネルギーから平均粒子数の拘束を引き、通常場 $h$ と対場 $\Delta$ をその密度微分として定義する。一般化密度の一次変分に対し $\delta(E-\lambda\Tr\rho)=\tfrac12\Tr(\mathcal H\delta\mathcal R)$ と書く規約では
$$
\mathcal H=
\begin{pmatrix}h-\lambda&\Delta\\-\Delta^*&-(h-\lambda)^*\end{pmatrix},
\qquad [\mathcal H,\mathcal R]=0 .
\label{eq-nuc-hfb}
$$
$\Delta^{\mathsf T}=-\Delta$ なので $\mathcal H$ は Hermite 行列である。粒子・穴を交換して複素共役を取る写像 $\mathcal C$ に対して、$\mathcal C\mathcal R\mathcal C^{-1}=1-\mathcal R$、$\mathcal C\mathcal H\mathcal C^{-1}=-\mathcal H$ である。許される微小変換は $\delta\mathcal R=[X,\mathcal R]$ であり、$X^\dagger=-X$、$\mathcal C X\mathcal C^{-1}=X$ を満たす。すると変分は $\tfrac12\Tr([\mathcal R,\mathcal H]X)$ となる。$[\mathcal R,\mathcal H]$ は反 Hermite で $\mathcal C$ の下で不変なので、その共役を $X$ に選べる。この変分がゼロなら交換子のノルム二乗がゼロであり、停留条件が従う。

$\mathcal H$ と $\mathcal R$ は互いに依存する。化学ポテンシャルを調整して粒子数を合わせながら場と密度を更新し、停留解のエネルギーや安定性を比較する。通常場と対場が共に時間反転した対ごとの二次元ブロックに分かれると BCS の式へ帰着する。連続状態の境界条件と接触型対汎関数の正則化は有限核の座標空間の理論を要するため、ここでは有限基底に限る。準粒子真空の変分原理の詳しい扱いは :cite[ring-schuck] を参照する。

:::example[対場を消す極限]
$\Delta=0$ で正常な基底状態を採ると、$\rho$ は $\lambda$ より低い軌道を占める射影、$\kappa=0$ となる。上左ブロックの交換子は :eqref[eq-nuc-hartree-fock] の停留条件に戻る。
:::

:::example[二つの空間の次元]
通常の単粒子基底が 6 次元なら、$\mathcal H$ は粒子と穴を合わせた 12 次元の行列である。粒子・穴対称性により固有値は正負に対をなし、ゼロがなければ 6 個の正エネルギーを取ればよい。
:::

### 一般化密度と数の平均

$\mathcal R^2=\mathcal R$ の上左成分は $\rho^2-\kappa\kappa^*=\rho$ である。反対称性 $\kappa^\dagger=-\kappa^*$ により $\rho-\rho^2=\kappa\kappa^\dagger$ となる。右辺は半正定値なので、$\rho$ の固有値 $v_i^2$ は $0\le v_i^2\le1$ を満たす。正準基底では、対の二つの軌道の占有が同じとき $\kappa$ の二行二列のブロックの非ゼロ成分の大きさは $\sqrt{v_i^2(1-v_i^2)}$ になる。占有が 0 または 1 ならその行で対振幅は消える。BCS でフェルミ面付近の軌道が重要になる理由を密度からも読める。

粒子と穴を交換し複素共役をとる反ユニタリ写像を $\mathcal C$ とすると、:eqref[eq-nuc-hfb] は $\mathcal C\mathcal H\mathcal C^{-1}=-\mathcal H$ を満たす。これは二つのブロックを交換して直接確認できる。従って $E$ の固有ベクトルを粒子・穴交換すれば $-E$ を得る。この対称性は通常の単粒子エネルギーの正負とは違い、同じ準粒子モードを作成と消滅の二通りで表したために生じる。

:::example[占有密度の範囲]
ある正準対の $v^2=1/4$ なら、対振幅の大きさは $\sqrt{3}/4$ である。$v^2=3/4$ に移してもこの絶対値は同じである。数の平均は各対に対する $2v^2$ で計算するため、この二つの値は同じ粒子数を与えない。
:::

### 密度行列から対場を導く

反対称化した二体行列要素を $\bar v_{ab,cd}=\langle ab|v|cd\rangle-\langle ab|v|dc\rangle$ とする。有限基底で Hamiltonian を $H=\sum_{ab}t_{ab}a_a^\dagger a_b+\frac14\sum_{abcd}\bar v_{ab,cd}a_a^\dagger a_b^\dagger a_d a_c$ と書く。準粒子真空の四演算子の期待値は
$$
\langle a_a^\dagger a_b^\dagger a_d a_c\rangle
=\rho_{ca}\rho_{db}-\rho_{da}\rho_{cb}+\kappa_{ab}^*\kappa_{cd}
$$
である。各 $a$ を準粒子の作成と消滅に展開し、右端の消滅が真空を消すことを使う。消滅を右に動かす際の反交換関係で二演算子ごとの積が残る。中央の並べ替えに伴う負号が第二項の負号であり、同種の演算子同士を結ぶ項が最後の異常密度の積である。従って
$$
E=\Tr(t\rho)+\tfrac12\sum_{abcd}\bar v_{ab,cd}\rho_{ca}\rho_{db}
+\tfrac14\sum_{abcd}\bar v_{ab,cd}\kappa_{ab}^*\kappa_{cd}
$$
となる。通常場は $h_{ac}=t_{ac}+\sum_{bd}\bar v_{ab,cd}\rho_{db}$、対場は $\Delta_{ab}=\frac12\sum_{cd}\bar v_{ab,cd}\kappa_{cd}$ である。後者は反対称化行列要素のため $\Delta_{ba}=-\Delta_{ab}$ となる。:eqref[eq-nuc-hfb] ではこの $h,\Delta$ で通常と異常のブロックを同時に更新する。

:::example[正常場だけの状態]
$\kappa=0$ なら異常項が消え、二体密度の残りは Slater 行列式の直接と交換に戻る。さらに正準軌道の占有が $0$ または $1$ なら、全エネルギーは [Hartree–Fock 理論](/physics/nuclear-physics/mean-field/) の期待値と一致する。
:::

### 正準条件と一般化密度の射影性

準粒子を $\beta_\mu=\sum_i(U_{i\mu}^*a_i+V_{i\mu}^*a_i^\dagger)$ とする。反交換関係を展開すると $\{\beta_\mu,\beta_\nu^\dagger\}=(U^\dagger U+V^\dagger V)_{\mu\nu}$、$\{\beta_\mu,\beta_\nu\}=(U^\dagger V^*+V^\dagger U^*)_{\mu\nu}$ である。従ってこの二つの行列を単位行列と零にすることが正準性の条件になる。変換の全ブロック行列を $\mathcal W$ とすればこの条件でユニタリとなり、真空の一般化密度は $\mathcal W \operatorname{diag}(0,1)\mathcal W^\dagger$ という射影の共役になる。よって $\mathcal R^2=\mathcal R$ である。この性質は純粋な準粒子真空に対応し、有限温度の統計演算子にはそのまま使えない。

:::example[二準位での正と負の固有値]
規格化された二成分のブロックに対角 $\pm\xi$ と非対角 $\Delta$ を入れる。行列式を取って正枝を求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
\xi=2\,\mathrm{MeV},\qquad \Delta=1\,\mathrm{MeV}
$$

この条件で求める量を ${E}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{E}&=\sqrt{\Delta^{2} + \xi^{2}}\\&=\sqrt{\left(1\right)^{2} + \left(2\right)^{2}}\\&\simeq 2.23607\,\mathrm{MeV}.\end{aligned}
$$

負枝は独立な別の粒子ではなく正枝の粒子・穴共役である。

対場の強い場合と比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
\Delta=2\,\mathrm{MeV}
$$

このときの値を ${E}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{E}_{\mathrm{new}}&=\sqrt{2} \cdot \sqrt{\left(2\right)^{2}}\\&\simeq 2.82843\,\mathrm{MeV},\\{E}_{\mathrm{new}}-{E}&\simeq 0.592359\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[異常密度と正常密度]
一正準対で密度を $n=v^2$ とする。射影性から異常密度の大きさは $\sqrt{n(1-n)}$ である。

量の定義に従い、共通の単位で次の値を用いる。
$$
n=\frac{1}{2}
$$

この条件で求める量を ${\kappa}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\kappa}&=\sqrt{n \cdot \left(1 - n\right)}\\&=\sqrt{\left(\frac{1}{2}\right) \cdot \left(1 - \left(\frac{1}{2}\right)\right)}\\&\simeq 0.5.\end{aligned}
$$

満杯か空では異常密度が消える。

ほぼ満杯の軌道を比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
n=\frac{9}{10}
$$

このときの値を ${\kappa}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\kappa}_{\mathrm{new}}&=\sqrt{\left(\frac{9}{10}\right) \cdot \left(1 - \left(\frac{9}{10}\right)\right)}\\&\simeq 0.3,\\{\kappa}_{\mathrm{new}}-{\kappa}&\simeq -0.2.\end{aligned}
$$
:::

:::example[対場の密度依存]
一つの正準対に対する有効行列要素を $g(\rho)=g_0(1-\rho/\rho_0)$ とする。$g_0>0$、$0\le\rho\le\rho_0$ で、無次元の異常密度成分の絶対値を $k$ とし、対場の大きさ $g(\rho)k$ を求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
g_{0}=8\,\mathrm{MeV},\qquad \rho=\frac{1}{25}\,\mathrm{fm}^{-3},\qquad \rho_{0}=\frac{4}{25}\,\mathrm{fm}^{-3},\qquad k=\frac{1}{5}
$$

この条件で求める量を ${\Delta}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\Delta}&=g_{0} \cdot k \cdot \left(1 - \frac{\rho}{\rho_{0}}\right)\\&=\left(8\right) \cdot \left(\frac{1}{5}\right) \cdot \left(- \frac{\left(\frac{1}{25}\right)}{\left(\frac{4}{25}\right)} + 1\right)\\&\simeq 1.2\,\mathrm{MeV}.\end{aligned}
$$

実際の自己無撞着計算では異常密度も同時に変わる。

異常密度を固定して結合だけの効果を取り出す。そこで次の値に変更し、ほかの入力を固定する。
$$
\rho=\frac{3}{25}\,\mathrm{fm}^{-3}
$$

このときの値を ${\Delta}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\Delta}_{\mathrm{new}}&=\left(8\right) \cdot \left(\frac{1}{5}\right) \cdot \left(- \frac{\left(\frac{3}{25}\right)}{\left(\frac{4}{25}\right)} + 1\right)\\&\simeq 0.4\,\mathrm{MeV},\\{\Delta}_{\mathrm{new}}-{\Delta}&\simeq -0.8\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[固定粒子数への射影確率]
二つの独立した正準対がともに確率 $p$ で占有される。二粒子を含む項の重みはどちらか一方だけが占有される確率である。

量の定義に従い、共通の単位で次の値を用いる。
$$
p=\frac{1}{2}
$$

この条件で求める量を ${P_2}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{P_2}&=2 \cdot p \cdot \left(1 - p\right)\\&=2 \cdot \left(\frac{1}{2}\right) \cdot \left(1 - \left(\frac{1}{2}\right)\right)\\&\simeq 0.5.\end{aligned}
$$

射影した状態はこの重みの平方根で割り、改めて規格化する。

平均占有が大きい場合と比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
p=\frac{4}{5}
$$

このときの値を ${P_2}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{P_2}_{\mathrm{new}}&=2 \cdot \left(\frac{4}{5}\right) \cdot \left(1 - \left(\frac{4}{5}\right)\right)\\&\simeq 0.32,\\{P_2}_{\mathrm{new}}-{P_2}&\simeq -0.18.\end{aligned}
$$

**一次近似の検証。** 半充填の周りでは第一導関数がゼロである。占有を動かしたときの射影確率を一次近似で評価してみる。

${p}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${P_2}$ の一階微分を求めると、
$$
\frac{\mathrm d{P_2}}{\mathrm dp}=2 - 4 \cdot p.
$$

この導関数を変更前の点で評価して $h=\frac{3}{10}$ を掛ける。ここで $p$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&=0,\\\Delta-\Delta_{\mathrm{lin}}&\simeq -0.18.\end{aligned}
$$

省いた分が全変化に比べて大きい。従って増分の大きさを予測するには、少なくとも二次の項を調べるか、元の式を使う必要がある。
:::

:::example[正準対の数の期待値]
二軌道の占有を $p_1,p_2$ とする。粒子数は一対につき 2 である。

量の定義に従い、共通の単位で次の値を用いる。
$$
p_{1}=\frac{7}{10},\qquad p_{2}=\frac{3}{10}
$$

この条件で求める量を ${\langle N\rangle}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\langle N\rangle}&=2 \cdot \left(p_{1} + p_{2}\right)\\&=2 \cdot \left(\frac{3}{10}\right) + 2 \cdot \left(\frac{7}{10}\right)\\&=2.\end{aligned}
$$

制約付き変分ではこの平均を目標値に合わせるように化学ポテンシャルを動かす。

二番目の占有だけを増やす。そこで次の値に変更し、ほかの入力を固定する。
$$
p_{2}=\frac{2}{5}
$$

このときの値を ${\langle N\rangle}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\langle N\rangle}_{\mathrm{new}}&=2 \cdot \left(\frac{2}{5}\right) + 2 \cdot \left(\frac{7}{10}\right)\\&\simeq 2.2,\\{\langle N\rangle}_{\mathrm{new}}-{\langle N\rangle}&\simeq 0.2.\end{aligned}
$$
:::

### 粒子・穴共役が固有値の対を作る

:eqref[eq-nuc-hfb] の上に並ぶ粒子成分と下に並ぶ穴成分の交換を $\mathcal C=\begin{pmatrix}0&1\\1&0\end{pmatrix}\mathcal K$ で表す。ここで $\mathcal K$ は成分の複素共役であり、各 $1$ は一体空間の単位行列である。$\mathcal H$ は対角に $h-\lambda$ と $-h^*+\lambda$、非対角に $\Delta$ と $-\Delta^*$ を持つ。複素共役してブロックを入れ替えれば $\mathcal C\mathcal H\mathcal C^{-1}=-\mathcal H$ である。従って $\mathcal H z=Ez$、$E$ が実数なら $\mathcal H(\mathcal C z)=-E(\mathcal C z)$ になる。この二つの固有ベクトルは同じ核子自由度の生成と消滅に対応する。$E\ne0$ なら Hermite 性により互いに直交するが、$E=0$ では共役が同じ固有空間へ戻るから単純に正と負のエネルギーには分けられない。零近傍の準粒子を含む系では占有する状態をどう選んだかを明示する。

## 対相関と慣性能率

核が角速度 $\omega$ で回る座標系では $H_\omega=H-\omega J_x$ となる。静止状態は非縮退で $J_x$ の期待値がゼロとする。二次の摂動論でエネルギーの変化は $-\omega^2\sum_{\nu\ne0}|\langle\nu|J_x|0\rangle|^2/(E_\nu-E_0)$ である。$E_\omega=E_0-\mathcal I\omega^2/2+\cdots$ と較べると応答の慣性能率 $\mathcal I$ が定まる。有限の回転慣性を求めるときは、以下のように固有軸を持つ変形した平均場を展開点に選ぶ。

変形した準粒子平均場を固定し、二準粒子だけを中間状態とする近似では
$$
\mathcal I_{\rm IB}=2\sum_{i<j}
\frac{|\langle ij|J_x|\Phi\rangle|^2}{E_i+E_j},\qquad
\ket{ij}=\alpha_i^\dagger\alpha_j^\dagger\ket\Phi.
$$
対相関は分母を大きくし、時間反転に奇な $J_x$ の典型的な行列要素に $u_i v_j-v_i u_j$ を生む。占有が近い軌道の寄与は小さくなる。平均場自体の回転による変化まで含めると誘起場が加わる。このため固定場の値をそのまま実測の慣性能率と同一視しない。

一体行列 $J$ が表す $J_x$ を $\hat J_x=\sum_{ab}J_{ab}a_a^\dagger a_b$ と書く。逆 Bogoliubov 変換 $a=U\beta+V^*\beta^\dagger$ を代入すると、$i<j$ の二準粒子作成の係数は $(U^\dagger J V^*)_{ij}-(U^\dagger J V^*)_{ji}$ である。二項目の負号は作成演算子を同じ順にそろえるときに付く。BCS の時間反転対では、$i,\bar j$ に対してこれは $u_i v_j J_{ij}+u_j v_i J_{\bar j\bar i}$ となる。時間反転に奇で Hermite な $J_x$ は $J_{\bar j\bar i}=-J_{ji}^*=-J_{ij}$ を満たすため、差のコヒーレンス因子が得られる。対場がゼロで準位が化学ポテンシャルをまたぐときは一方の $v$ が 1、他方が 0 となり、応答は粒子・穴遷移による Inglis の式に戻る。

:::example[一つの中間状態]
$|\langle ij|J_x|\Phi\rangle|=\hbar$、$E_i+E_j=2\,\mathrm{MeV}$ なら、その項は $1\,\hbar^2/\mathrm{MeV}$ である。分母だけを 2 倍にすれば寄与は半分になる。
:::

:::example[コヒーレンス因子]
$u_i=v_i=u_j=v_j=1/\sqrt2$ の場合、$u_i v_j-v_i u_j=0$ となる。遷移を抑えるのは分母の有限ギャップだけではなく、粒子成分と穴成分の打ち消しもある。
:::

### 固定した軸と回転の応答

連続な回転対称性を保つ厳密な有限核の Hamiltonian では $[H,J_x]=0$ である。その場合、非縮退の $J=0$ 固有状態に $J_x$ を作用させるとゼロとなり、固定した非縮退状態の二次応答式は有限の慣性を与えない。慣性を得るクランキング近似は、固有方向をもつ変形した平均場を作り、そこから回転への応答を調べる。この慣性を使って回転バンドを量子化するとき、回転で向きの違う状態を集める。球状核の厳密な $J=0$ 固有状態と変形した固有状態の区別が必要である。

### 線形応答と二準粒子近似

$H_\omega=H-\omega J_x$ の一次の波動関数変化を求めると、励起 $\ket{\nu}$ の係数は $\omega\langle\nu|J_x|0\rangle/(E_\nu-E_0)$ である。この一次誤差を交差項として $\langle J_x\rangle$ に入れると慣性能率は $2\sum_{\nu\ne0}|\langle\nu|J_x|0\rangle|^2/(E_\nu-E_0)$ となる。非縮退の対称な真空を仮定している。平均場と対場の変化を固定して二準粒子状態で置き換えたものが Inglis–Belyaev の近似である。一体演算子 $J_x$ を準粒子で書くと、真空からの励起では生成演算子二つの項だけが残る。分母は二つの準粒子エネルギーの和になり、分子には対混合の係数が入る。誘起された場も含む応答は同じ近似では計算できず、残留相互作用を含む応答方程式が必要である。

:::example[一つの励起が支配する近似]
行列要素を $\hbar j$、励起エネルギーを $D$ とする。一項だけの応答を計算する。

量の定義に従い、共通の単位で次の値を用いる。
$$
j=2,\qquad D=2\,\mathrm{MeV}
$$

この条件で求める量を ${\mathcal I/\hbar^2}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\mathcal I/\hbar^2}&=\frac{2 \cdot j^{2}}{D}\\&=2 \cdot \left(2\right)\\&=4\,\mathrm{MeV}^{-1}.\end{aligned}
$$

励起が高くなると同じ行列要素の静的応答は小さくなる。

分母を二倍にする。そこで次の値に変更し、ほかの入力を固定する。
$$
D=4\,\mathrm{MeV}
$$

このときの値を ${\mathcal I/\hbar^2}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\mathcal I/\hbar^2}_{\mathrm{new}}&=\frac{2 \cdot \left(2\right)^{2}}{\left(4\right)}\\&=2\,\mathrm{MeV}^{-1},\\{\mathcal I/\hbar^2}_{\mathrm{new}}-{\mathcal I/\hbar^2}&=-2\,\mathrm{MeV}^{-1}.\end{aligned}
$$
:::

:::example[対ギャップで変わる二準粒子の分母]
二軌道の $\xi$ を $\pm d$、対ギャップを $\Delta$ とする。二準粒子のエネルギーを足す。

量の定義に従い、共通の単位で次の値を用いる。
$$
d=1\,\mathrm{MeV},\qquad \Delta=0\,\mathrm{MeV}
$$

この条件で求める量を ${D}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{D}&=2 \cdot \sqrt{\Delta^{2} + d^{2}}\\&=2 \cdot \sqrt{\left(0\right)^{2} + \left(1\right)^{2}}\\&=2\,\mathrm{MeV}.\end{aligned}
$$

実際の慣性能率では分子の対混合も同時に変化する。

対場を入れない基準と比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
\Delta=1\,\mathrm{MeV}
$$

このときの値を ${D}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{D}_{\mathrm{new}}&=2 \cdot \sqrt{2} \cdot \sqrt{\left(1\right)^{2}}\\&\simeq 2.82843\,\mathrm{MeV},\\{D}_{\mathrm{new}}-{D}&\simeq 0.828427\,\mathrm{MeV}.\end{aligned}
$$

**一次近似の検証。** 正常状態の付近では分母が対場の符号反転に対して偶である。そのため一次近似で見落とす変化を調べる。

${\Delta}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${D}$ の一階微分を求めると、
$$
\frac{\mathrm d{D}}{\mathrm d\Delta}=\frac{2 \cdot \Delta}{\sqrt{\Delta^{2} + d^{2}}}.
$$

この導関数を変更前の点で評価して $h=1$ を掛ける。ここで $\Delta$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&=0\,\mathrm{MeV},\\\Delta-\Delta_{\mathrm{lin}}&\simeq 0.828427\,\mathrm{MeV}.\end{aligned}
$$

省いた分が全変化に比べて大きい。従って増分の大きさを予測するには、少なくとも二次の項を調べるか、元の式を使う必要がある。
:::

:::example[等しい軌道の打ち消し]
時間反転に奇な一体演算子のある対の係数を $u_1v_2-v_1u_2$ とする。混合角 $\theta_i$ ならその二乗は $\sin^2((\theta_2-\theta_1)/2)$ になる。

量の定義に従い、共通の単位で次の値を用いる。
$$
t_{1}=\frac{\pi}{3}\,\mathrm{rad},\qquad t_{2}=\frac{\pi}{3}\,\mathrm{rad}
$$

この条件で求める量を ${w}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{w}&=\sin^{2}{\left(\frac{- t_{1} + t_{2}}{2} \right)}\\&=0\\&=0.\end{aligned}
$$

二軌道の占有のにじみが同じならこの項はゼロである。

異なる混合角を比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
t_{2}=\frac{2 \cdot \pi}{3}\,\mathrm{rad}
$$

このときの値を ${w}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{w}_{\mathrm{new}}&=\sin^{2}{\left(\frac{\left(\frac{2 \cdot \pi}{3}\right)}{2} - \frac{\left(\frac{\pi}{3}\right)}{2} \right)}\\&\simeq 0.25,\\{w}_{\mathrm{new}}-{w}&\simeq 0.25.\end{aligned}
$$
:::

:::example[回転準位からの慣性能率]
剛体に近い $K=0$ 回転で $E(2)-E(0)=3\hbar^2/\mathcal I$ とする。

量の定義に従い、共通の単位で次の値を用いる。
$$
E_{2}=\frac{3}{20}\,\mathrm{MeV}
$$

この条件で求める量を ${\mathcal I/\hbar^2}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\mathcal I/\hbar^2}&=\frac{3}{E_{2}}\\&=\frac{3}{\left(\frac{3}{20}\right)}\\&=20\,\mathrm{MeV}^{-1}.\end{aligned}
$$

回転帯の同定をしたうえでの推定である。球形振動の $2^+$ にはこの換算を直接使えない。

より低い初励起の回転帯と比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
E_{2}=\frac{3}{40}\,\mathrm{MeV}
$$

このときの値を ${\mathcal I/\hbar^2}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\mathcal I/\hbar^2}_{\mathrm{new}}&=\frac{3}{\left(\frac{3}{40}\right)}\\&=40\,\mathrm{MeV}^{-1},\\{\mathcal I/\hbar^2}_{\mathrm{new}}-{\mathcal I/\hbar^2}&=20\,\mathrm{MeV}^{-1}.\end{aligned}
$$

**一次近似の検証。** 励起エネルギーの誤差が慣性能率にどう伝わるか、同じ帯の基準で確認する。

${E_{2}}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${\mathcal I/\hbar^2}$ の一階微分を求めると、
$$
\frac{\mathrm d{\mathcal I/\hbar^2}}{\mathrm dE_{2}}=- \frac{3}{E_{2}^{2}}.
$$

この導関数を変更前の点で評価して $h=- \frac{3}{40}$ を掛ける。ここで $E_{2}$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&=10\,\mathrm{MeV}^{-1},\\\Delta-\Delta_{\mathrm{lin}}&=10\,\mathrm{MeV}^{-1}.\end{aligned}
$$

省いた分が全変化に比べて大きい。従って増分の大きさを予測するには、少なくとも二次の項を調べるか、元の式を使う必要がある。
:::

:::example[応答の加算]
二励起の $J_x$ 行列要素を $\hbar j_1,\hbar j_2$ とする。同じ初期状態から別々の終状態なので二乗を足す。

量の定義に従い、共通の単位で次の値を用いる。
$$
j_{1}=1,\qquad j_{2}=2,\qquad D_{1}=\frac{3}{2}\,\mathrm{MeV},\qquad D_{2}=4\,\mathrm{MeV}
$$

この条件で求める量を ${\mathcal I/\hbar^2}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\mathcal I/\hbar^2}&=\frac{2 \cdot j_{2}^{2}}{D_{2}} + \frac{2 \cdot j_{1}^{2}}{D_{1}}\\&=\frac{2 \cdot \left(1\right)^{2}}{\left(\frac{3}{2}\right)} + \frac{2 \cdot \left(2\right)^{2}}{\left(4\right)}\\&\simeq 3.33333\,\mathrm{MeV}^{-1}.\end{aligned}
$$

強い遷移でも励起エネルギーが高ければ応答の寄与が抑えられる。

第二励起を高くする。そこで次の値に変更し、ほかの入力を固定する。
$$
D_{2}=8\,\mathrm{MeV}
$$

このときの値を ${\mathcal I/\hbar^2}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\mathcal I/\hbar^2}_{\mathrm{new}}&=\frac{2 \cdot \left(1\right)^{2}}{\left(\frac{3}{2}\right)} + \frac{2 \cdot \left(2\right)^{2}}{\left(8\right)}\\&\simeq 2.33333\,\mathrm{MeV}^{-1},\\{\mathcal I/\hbar^2}_{\mathrm{new}}-{\mathcal I/\hbar^2}&=-1\,\mathrm{MeV}^{-1}.\end{aligned}
$$
:::

## 演習問題

::::exercise{level=basic}
偶数 $N$ の周辺で $B(N-1)=50$、$B(N)=58$、$B(N+1)=63\,\mathrm{MeV}$ とする。三点差分を求めよ。
:::solution
:eqref[eq-nuc-pairing-mass] により $\Delta^{(3)}=(116-113)/2=1.5\,\mathrm{MeV}$ である。
:::
::::

::::exercise{level=basic}
$\Omega=3,n=2$ の $v=0,2$ の相関エネルギーを求めよ。
:::solution
:eqref[eq-nuc-seniority-energy] により $v=0$ では $-G(2)(6)/4=-3G$、$v=2$ では $n-v=0$ なのでゼロである。
:::
::::

::::exercise{level=basic}
$\xi=0$、$\Delta=1.2\,\mathrm{MeV}$ の軌道の準粒子エネルギーと占有確率を求めよ。
:::solution
:eqref[eq-nuc-bcs-gap] より $E=|\Delta|=1.2\,\mathrm{MeV}$、$v^2=(1-0/E)/2=1/2$ である。
:::
::::

::::exercise{level=basic}
$\Omega=4$ の半充填縮退模型で $G=0.5\,\mathrm{MeV}$ とする。BCS のギャップを求めよ。
:::solution
各軌道の $uv=1/2$ なので $\Delta=G\Omega/2=1\,\mathrm{MeV}$ である。
:::
::::

::::exercise{level=advanced}
$S_z=(n-\Omega)/2$ と $s=(\Omega-v)/2$ から、$S_+S_-$ の固有値を因数分解せよ。
:::solution
角運動量の代数から $s(s+1)-S_z(S_z-1)=(s+S_z)(s-S_z+1)$ である。代入して第 1 因子は $(n-v)/2$、第 2 因子は $(2\Omega-n-v+2)/2$ となる。
:::
::::

::::exercise{level=advanced}
二つの対軌道で $\xi=\pm\epsilon$、$\epsilon>0$ とする。非ゼロの BCS 解が現れる $G$ の条件を求めよ。
:::solution
両方の $E$ が等しいので :eqref[eq-nuc-bcs-gap] は $1=G/\sqrt{\epsilon^2+\Delta^2}$ になる。$\Delta=\sqrt{G^2-\epsilon^2}>0$ であるための条件は $G>\epsilon$ である。化学ポテンシャルを真ん中に置くと二軌道の $v^2$ の和が 1 なので平均粒子数は 2 である。
:::
::::

::::exercise{level=advanced}
BCS 状態の粒子数分散を計算せよ。
:::solution
各対の粒子数は確率 $u_k^2$ で 0、確率 $v_k^2$ で 2 である。従って分散は $4v_k^2-(2v_k^2)^2=4u_k^2v_k^2$ となる。状態は対ごとの積であり、対をまたぐ共分散はゼロなので全分散は $4\sum_k u_k^2v_k^2$ である。
:::
::::

::::exercise{level=advanced}
$H_\omega=H-\omega J_x$ で $E_\omega=E_0-\mathcal I\omega^2/2+O(\omega^3)$ とする。$J_x$ の一次応答を求めよ。
:::solution
規格化固有状態を $\omega$ で変分すると、状態の微分を含む二項が規格化の微分で消え、$\dd E_\omega/\dd\omega=-\langle J_x\rangle_\omega$ だけが残る。従って $\langle J_x\rangle_\omega=\mathcal I\omega+O(\omega^2)$ となる。
:::
::::

## 確認テスト

::::quiz
一殻の一定対相互作用で $S_+$ を一回作用させるとどうなるか。結果がゼロでない場合を考える。

- [x] 粒子数は 2 増え、seniority は変わらない。
- [ ] 粒子数と seniority がともに 2 増える。
- [ ] 粒子数は一定で seniority が 2 減る。

:::explanation
$[\hat n,S_+]=2S_+$ により粒子数は 2 増える。$S_+$ は同じ準スピン表現内の上げ演算子なので $s$ は固定され、$v=\Omega-2s$ も固定される。
:::
::::

::::quiz{answer="5" tol="1e-6"}
$\xi=3\,\mathrm{MeV}$、$\Delta=4\,\mathrm{MeV}$ の準粒子エネルギーを MeV で求めよ。
:::explanation
準粒子行列の正の固有値は $\sqrt{3^2+4^2}=5\,\mathrm{MeV}$ である。
:::
::::

::::quiz
ブロッキングについて正しい文を選べ。

- [x] 単独で占有された軌道は、その同じ時間反転対の対散乱に使えない。
- [ ] 全軌道の占有確率をゼロにする。
- [ ] 一軌道を占有すると残りのすべての対が消える。

:::explanation
Pauli 原理により、単独占有軌道にもう一度同じ状態の核子を作れない。その対を和から除くが、別の対の散乱は残る。
:::
::::

::::quiz{answer="0.5" tol="1e-6"}
$\xi=0$ かつ $\Delta>0$ の軌道の $v^2$ を求めよ。
:::explanation
:eqref[eq-nuc-bcs-gap] で $\xi/E=0$ なので $v^2=1/2$ である。
:::
::::

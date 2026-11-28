---
title: 弾性体の基礎方程式
description: 線形弾性の構成則と平衡式を導き、対称性を使った境界値問題を解く。
notation:
  - symbol: 'C_{ijkl}'
    meaning: 線形弾性定数テンソル
  - symbol: '\lambda,\mu'
    meaning: Lamé 定数
  - symbol: 'E,\nu,K'
    meaning: Young 率、Poisson 比、体積弾性率
  - symbol: '\eps'
    meaning: 微小ひずみ
---

## この章の目標

- 線形弾性定数の対称性と等方性の意味を理解する。
- 弾性率を変換し、拘束による応力の違いを計算する。
- 平衡の境界値問題を変位とエネルギーの両方から解く。
- ひずみから変位を再構成できる条件を導く。

## Hooke の法則と弾性定数テンソル

バネを自然長の近くで伸ばすと、復元力の一次項が変位に比例する。連続体では伸びやせん断を含む六つの成分があり、一つの成分の変形が複数の応力を生む。ここでは応力ゼロの参照状態をとり、変位勾配が十分小さく、温度は一定であるとする。成分はデカルト座標で表し、単位は SI である。

:::law[線形弾性則]{#law-cont-hooke}
可逆な弾性応答の一次近似は

$$
\sigma_{ij}=C_{ijkl}\eps_{kl},
\qquad
\eps_{kl}=(\partial_k u_l+\partial_lu_k)/2
\label{eq-cont-hooke}
$$

である。:term[弾性定数]{en="elastic stiffness" yomi="だんせいていすう"} $C_{ijkl}$ の単位は $\mathrm{Pa}$ である。
:::

応力とひずみの対称性から $C_{ijkl}=C_{jikl}=C_{ijlk}$ とできる。前半は任意の対称ひずみに対して応力が対称であるため、後半はひずみとの縮約に寄与する部分を選ぶためである。ひずみエネルギー $W$ の微分が $\mathrm dW=\sigma:\mathrm d\eps$ なら混合二階微分の交換から $C_{ijkl}=C_{klij}$ も成り立つ。従って六次元の対称行列として独立定数は高々 $6\cdot7/2=21$ 個である。

Voigt 表記の一つではひずみの列を $(\eps_{11},\eps_{22},\eps_{33},2\eps_{23},2\eps_{31},2\eps_{12})$、応力の列を $(\sigma_{11},\sigma_{22},\sigma_{33},\sigma_{23},\sigma_{31},\sigma_{12})$ とする。列の通常の内積がテンソルの二重内積に一致する。せん断に係数 $2$ があるため、Voigt 行列に入れる係数とテンソルの縮約に現れる項数を区別する。

:::example[せん断試験]
$\eps_{12}=\eps_{21}=10^{-3}$ で、$C_{1212}=C_{1221}=30\,\mathrm{GPa}$、その他のひずみはゼロとする。$\sigma_{12}=60\,\mathrm{MPa}$ である。工学せん断ひずみ $2\times10^{-3}$ に $30\,\mathrm{GPa}$ を掛けても同じ値になる。
:::

結晶軸に沿う立方対称な材料を考える。座標軸まわりの半回転で符号が変わる成分はゼロである。軸交換に対して同じ値を取るため、非零成分は $C_{1111},C_{1122},C_{1212}$ の三種類に分かれる。したがって $[100]$ 方向と $[111]$ 方向で弾性応答が異なる場合がある。

### 結晶軸を用いた立方晶の法則

ひずみエネルギーが二次形式で、座標軸が立方体の辺に一致するとする。各軸の符号反転は同時に二つずつ行っても結晶を変えない。そのため、垂直成分とせん断成分の混合積、および異なるせん断成分同士の積は消える。軸を入れ替えてもエネルギーは変わらないから、残った項の係数は三種類である。$\eps_{ij}=\eps_{ji}$ とし、垂直成分の自乗の係数を $C_{11}/2$、異なる垂直成分の積の係数を $C_{12}$、せん断成分の自乗の係数を $2C_{44}$ とする。例えば $\sigma_{11}=C_{11}\eps_{11}+C_{12}(\eps_{22}+\eps_{33})$、$\sigma_{12}=2C_{44}\eps_{12}$ となる。対称ひずみの片側成分だけを変えるときは、仕事に $\sigma_{12}\delta\eps_{12}+\sigma_{21}\delta\eps_{21}$ の両方が入る。

垂直成分に関する行列は $(C_{11}-C_{12})I+C_{12}\vb 1\vb 1^T$。方向 $\vb 1=(1,1,1)^T$ では固有値 $C_{11}+2C_{12}$、その直交補空間では $C_{11}-C_{12}$ である。これら二つと $C_{44}$ が正ならひずみエネルギーは正定値であり、どれかが非正なら対応するひずみで正定値性が破れる。独立したせん断試験から求まる $C_{44}$ は、等方条件では $(C_{11}-C_{12})/2$ と一致する。

:::example[等方性と整合する立方晶定数]
$C_{11}=200 \mathrm{GPa}$、$C_{12}=100 \mathrm{GPa}$、$C_{44}=50 \mathrm{GPa}$ とする。以下では$K$ が等方体積変形に対する係数、$G_{d}$ が垂直偏差ひずみに対するせん断係数、$Z$ がせん断係数の比である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
K&=\frac{C_{11}}{3} + \frac{2 C_{12}}{3}=\frac{400}{3}\quad[\mathrm{GPa}]\\
G_{d}&=\frac{C_{11}}{2} - \frac{C_{12}}{2}=50\quad[\mathrm{GPa}]\\
Z&=\frac{C_{44}}{\frac{C_{11}}{2} - \frac{C_{12}}{2}}=1
\end{aligned}
$$

比は１であり、この場合に限り垂直偏差ひずみとせん断を同じ係数で記述できる。
:::

:::example[せん断方向によって変わる弾性率]
$C_{11}=200 \mathrm{GPa}$、$C_{12}=140 \mathrm{GPa}$、$C_{44}=60 \mathrm{GPa}$ とする。以下では$K$ が等方体積変形に対する係数、$G_{d}$ が垂直偏差ひずみに対するせん断係数、$Z$ がせん断係数の比である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
K&=\frac{C_{11}}{3} + \frac{2 C_{12}}{3}=160\quad[\mathrm{GPa}]\\
G_{d}&=\frac{C_{11}}{2} - \frac{C_{12}}{2}=30\quad[\mathrm{GPa}]\\
Z&=\frac{C_{44}}{\frac{C_{11}}{2} - \frac{C_{12}}{2}}=2
\end{aligned}
$$

体積変化を伴わなくても、座標軸に沿う伸縮と座標面のせん断では剛性が異なる。比 $Z$ はその差を示す。
:::

:::example[立方晶の異方性と結合係数]
$C_{11}=200 \mathrm{GPa}$、$C_{44}=50 \mathrm{GPa}$ に固定し、別の軸との結合 $C_{12}$ を $C_{12}=100 \mathrm{GPa}$ の近くで変える。せん断係数の比 $Z$ の局所的な変化を求める。固定した値を先に代入すると、残る一変数の関数は次のようになる。以下の微分では数値を先の単位系で表す。

$$
\begin{aligned}
Z(C_{12})&=\frac{50}{100 - \frac{C_{12}}{2}},\\
\frac{\mathrm dZ}{\mathrm dC_{12}}&=\frac{100}{\left(C_{12} - 200\right)^{2}}.
\end{aligned}
$$

基準値で $Z=1$、微分は $\frac{1}{100}$ である。したがって一次の展開は、基準値からの相対変化を $C_{12}\mapsto C_{12}(1+\delta)$ と書くと

$$
\frac{\Delta Z}{Z}
=\left(\frac{C_{12}}{Z}
\frac{\mathrm dZ}{\mathrm dC_{12}}\right)\delta+\order{\delta^2}
=1\,\delta+\order{\delta^2}
$$

となる。例えば入力を $1\%$ 増やすとき、一次式で得る出力の変化率は $1\%$ である。元の式に直接代入した変化率は $\frac{100}{99}\%$ なので、この比較で一次化の誤差を点検できる。

逆に、出力の測定値が基準値から $1\%$ 上がったとき、適用した模型が正しく、ほかの入力が固定されていれば、入力の変化率は一次で $1\%$ と推定される。$C_{12}$ が $C_{11}$ に近づくと垂直偏差ひずみの剛性が小さくなる。正定値条件とその余裕を併せて確認する。
:::

## 等方弾性体と Lamé 定数

:::proposition[等方弾性の二定数]{#prop-cont-isotropic-stress}
回転で構成式が変わらない均質な線形弾性体では、二つの :term[Lamé 定数]{en="Lamé constants" yomi="らめていすう"} $\lambda,\mu$ によって

$$
\sigma=2\mu\eps+\lambda(\tr\eps)I
=2\mu\eps_{\mathrm d}+K(\tr\eps)I,\qquad K=\lambda+\frac23\mu,
\label{eq-cont-isotropic}
$$

と表せる。ここで $\eps_{\mathrm d}=\eps-\tr(\eps)I/3$ である。
:::

:::proof
ひずみの実対称性により固有基底を選べる。対角ひずみは座標軸周りの半回転で保たれるため、応力の各非対角成分はその符号を反転する半回転によってゼロと分かる。対角成分は三つのひずみ固有値の線形関数であり、軸の置換で同じ係数を持つ。軸交換に必要なら一軸の符号を反転して固有回転にでき、対角ひずみへの作用は変わらない。従って $\sigma_{ii}=a\eps_{ii}+b\sum_{j\ne i}\eps_{jj}$ で、$2\mu=a-b$、$\lambda=b$ とおけば所要の式になる。基底を戻しても $I,\tr\eps$ は不変なので成立する。
:::

偏差ひずみは一次で体積を変えないため、$\mu$ は形を変える剛性を表す。純粋な膨張では偏差成分が消え、$K$ が体積変化に対する割合を決める。

:::example[等方ひずみと等容ひずみ]
$\lambda=40\,\mathrm{GPa}$、$\mu=30\,\mathrm{GPa}$ なら $K=60\,\mathrm{GPa}$ である。$\eps=10^{-3}I$ では $\sigma=180 I\,\mathrm{MPa}$ となる。$\eps=10^{-3}\diag(1,-1,0)$ では $\sigma=60\diag(1,-1,0)\,\mathrm{MPa}$ で、平均応力は消える。
:::

### 偏差ひずみを投影してエネルギーを分ける

$\theta=\tr\eps$、$e=\eps-\theta I/3$ とすると、$I:e=\tr e=0$ なので $\eps:\eps=e:e+\theta^2/3$。等方弾性の二次エネルギーは $W=\mu e:e+K\theta^2/2$。これを微分するか、Hooke 則の跡と偏差を取れば $\sigma=2\mu e+K\theta I$ を得る。体積変化を固定してせん断する場合に Lamé の第一定数が現れないのは、この直交分解のためである。

例えば $\eps_{11}=\theta/3+d$、$\eps_{22}=\theta/3-d$、$\eps_{33}=\theta/3$ とすれば $e:e=2d^2$。体積応力 $K\theta$ に $2\mu d,-2\mu d,0$ を加えた値が三主応力である。$K$ を非常に大きくして有限の平均応力を保つには $\theta$ が $1/K$ に比例して小さくならなければならない。ここから、非圧縮近似では体積条件を課して圧力を別の未知量にする理由が分かる。

:::example[体積ひずみと偏差ひずみの同時負荷]
$K=100 \mathrm{GPa}$、$\mu=40 \mathrm{GPa}$、$\theta=\frac{3}{1000}$、$d=\frac{1}{1000}$ とする。以下では$p_{m}$ が平均垂直応力、$W_{v}$ が体積変形のエネルギー、$W_{s}$ が形状変化のエネルギー、$Q$ が形状変化と体積変化のエネルギー比である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
p_{m}&=K \theta=\frac{3}{10}\quad[\mathrm{GPa}]\\
W_{v}&=500 K \theta^{2}=\frac{9}{20}\quad[\mathrm{MJ\,m^{-3}}]\\
W_{s}&=2000 d^{2} \mu=\frac{2}{25}\quad[\mathrm{MJ\,m^{-3}}]\\
Q&=\frac{4 d^{2} \mu}{K \theta^{2}}=\frac{8}{45}
\end{aligned}
$$

体積エネルギーとせん断エネルギーは交差項を持たない。応力の平均に現れるのは体積ひずみだけである。
:::

:::example[平均応力を一定に保ったまま硬くする]
$K=200 \mathrm{GPa}$、$\mu=40 \mathrm{GPa}$、$\theta=\frac{3}{2000}$、$d=\frac{1}{1000}$ とする。以下では$p_{m}$ が平均垂直応力、$W_{v}$ が体積変形のエネルギー、$W_{s}$ が形状変化のエネルギー、$Q$ が形状変化と体積変化のエネルギー比である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
p_{m}&=K \theta=\frac{3}{10}\quad[\mathrm{GPa}]\\
W_{v}&=500 K \theta^{2}=\frac{9}{40}\quad[\mathrm{MJ\,m^{-3}}]\\
W_{s}&=2000 d^{2} \mu=\frac{2}{25}\quad[\mathrm{MJ\,m^{-3}}]\\
Q&=\frac{4 d^{2} \mu}{K \theta^{2}}=\frac{16}{45}
\end{aligned}
$$

平均応力は上の例と同じである。体積エネルギーは半分になり、形状変化の分が占める割合が増える。
:::

:::example[小さい体積ひずみの測定精度]
$K=100 \mathrm{GPa}$、$\mu=40 \mathrm{GPa}$、$d=\frac{1}{1000}$ に固定し、体積ひずみ $\theta$ を $\theta=\frac{3}{1000}$ の近くで変える。形状変化と体積変化のエネルギー比 $Q$ の局所的な変化を求める。固定した値を先に代入すると、残る一変数の関数は次のようになる。以下の微分では数値を先の単位系で表す。

$$
\begin{aligned}
Q(\theta)&=\frac{1}{625000 \theta^{2}},\\
\frac{\mathrm dQ}{\mathrm d\theta}&=- \frac{1}{312500 \theta^{3}}.
\end{aligned}
$$

基準値で $Q=\frac{8}{45}$、微分は $- \frac{3200}{27}$ である。したがって一次の展開は、基準値からの相対変化を $\theta\mapsto \theta(1+\delta)$ と書くと

$$
\frac{\Delta Q}{Q}
=\left(\frac{\theta}{Q}
\frac{\mathrm dQ}{\mathrm d\theta}\right)\delta+\order{\delta^2}
=-2\,\delta+\order{\delta^2}
$$

となる。例えば入力を $1\%$ 増やすとき、一次式で得る出力の変化率は $-2\%$ である。元の式に直接代入した変化率は $- \frac{20100}{10201}\%$ なので、この比較で一次化の誤差を点検できる。

逆に、出力の測定値が基準値から $1\%$ 上がったとき、適用した模型が正しく、ほかの入力が固定されていれば、入力の変化率は一次で $- \frac{1}{2}\%$ と推定される。せん断が同じとき体積エネルギーの相対変化は、体積ひずみの相対変化の２倍である。比ではその符号が逆になる。
:::

## Young 率・Poisson 比・体積弾性率

:::definition[一軸応力試験]{#def-cont-elastic-moduli}
一軸応力 $\sigma_{11}=S$、他はゼロのとき、:term[Young 率]{en="Young modulus" yomi="やんぐりつ"} を $E=S/\eps_{11}$、:term[Poisson 比]{en="Poisson ratio" yomi="ぽあそんひ"} を $\nu=-\eps_{22}/\eps_{11}$ と定める。$\mu$ はせん断弾性率、$K$ は平均引張応力を体積ひずみで割った :term[体積弾性率]{en="bulk modulus" yomi="たいせきだんせいりつ"} である。
:::

:ref[prop-cont-isotropic-stress] の跡をとって $\tr\eps=\tr\sigma/(3K)$ を得る。偏差部分は $\eps_{\mathrm d}=\sigma_{\mathrm d}/(2\mu)$ である。結合すると

$$
\eps=\frac{1+\nu}{E}\sigma-\frac{\nu}{E}\tr(\sigma)I,\quad
E=\frac{\mu(3\lambda+2\mu)}{\lambda+\mu},\quad
\nu=\frac{\lambda}{2(\lambda+\mu)}
\label{eq-cont-compliance}
$$

となる。実際、一軸応力の各成分は $\eps_{11}=S(1/(3\mu)+1/(9K))$、$\eps_{22}=S(-1/(6\mu)+1/(9K))$ であり、定義に代入すると二定数の値が出る。逆に $\mu=E/[2(1+\nu)]$、$K=E/[3(1-2\nu)]$、$\lambda=K-2\mu/3$ である。$\mu,K>0$ の範囲は $E>0$、$-1<\nu<1/2$ に対応する。

### 面内の近似

薄板の面に垂直な応力が小さい平面応力では、$\sigma_{33}=\sigma_{13}=\sigma_{23}=0$ である。:eqref[eq-cont-compliance] より $\eps_{33}=-\nu(\sigma_{11}+\sigma_{22})/E$ なので、厚さは変わりうる。面内二成分を逆に解くと $\sigma_{11}=E(\eps_{11}+\nu\eps_{22})/(1-\nu^2)$、$\sigma_{22}=E(\eps_{22}+\nu\eps_{11})/(1-\nu^2)$ である。$\vb u=(u_x(x,y),u_y(x,y),0)$ の平面ひずみでは $\eps_{33}=\eps_{13}=\eps_{23}=0$ なので $\sigma_{33}=\nu(\sigma_{11}+\sigma_{22})$ となる。

:::example[横拘束]
$E=200\,\mathrm{GPa}$、$\nu=0.25$ の材料を $\eps_{11}=10^{-3}$ だけ伸ばす。側面自由なら $\sigma_{11}=200\,\mathrm{MPa}$ である。$\eps_{22}=\eps_{33}=0$ に拘束すると $\mu=80\,\mathrm{GPa}$、$\lambda=80\,\mathrm{GPa}$ なので $\sigma_{11}=240\,\mathrm{MPa}$、横応力は $80\,\mathrm{MPa}$ となる。
:::

### 一軸応力試験と一軸ひずみ試験をつなぐ

等方体の一軸応力 $s$ からは軸ひずみ $s/E$、横ひずみ $-\nu s/E$ を得る。横ひずみを両方ゼロにするには、新たに等しい横応力 $q$ を加える。横ひずみの式は $(q-\nu(q+s))/E=0$ だから $q=\nu s/(1-\nu)$。これを軸ひずみ $(s-2\nu q)/E$ に代入すると、拘束下の見かけの弾性率は $E(1-\nu)/[(1+\nu)(1-2\nu)]$ となる。側面自由の試験と違い、横方向の拘束具が反力を負担する。

$E$ を正の値に固定して $\nu\to1/2$ とすると、この係数が発散する。一軸ひずみなら体積ひずみがそのまま軸ひずみとなるため、ほぼ非圧縮の物体は側方にも動けないと軸方向にも伸びにくい。$\nu<0$ なら軸引張の際の自由膨張は横方向でも正になり、拘束に必要な横応力は圧縮になる。

:::example[横収縮を止めるために必要な力]
$E=70 \mathrm{GPa}$、$\nu=\frac{1}{3}$、$s=21 \mathrm{MPa}$ とする。以下では$e_{f}$ が横自由の軸ひずみ、$q$ が横拘束のための反力応力、$e_{c}$ が横拘束中の軸ひずみ、$M$ が横拘束下の見かけの弾性率である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
e_{f}&=\frac{s}{1000 E}=\frac{3}{10000}\\
q&=\frac{\nu s}{1 - \nu}=\frac{21}{2}\quad[\mathrm{MPa}]\\
e_{c}&=\frac{- \frac{2 \nu^{2} s}{1 - \nu} + s}{1000 E}=\frac{1}{5000}\\
M&=\frac{E s}{- \frac{2 \nu^{2} s}{1 - \nu} + s}=105\quad[\mathrm{GPa}]
\end{aligned}
$$

側面の反力は引張である。両側の横収縮を打ち消すので、同じ軸応力でも軸ひずみは横自由のときより小さい。
:::

:::example[負の Poisson 比と横圧縮]
$E=70 \mathrm{GPa}$、$\nu=- \frac{1}{4}$、$s=21 \mathrm{MPa}$ とする。以下では$e_{f}$ が横自由の軸ひずみ、$q$ が横拘束のための反力応力、$e_{c}$ が横拘束中の軸ひずみ、$M$ が横拘束下の見かけの弾性率である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
e_{f}&=\frac{s}{1000 E}=\frac{3}{10000}\\
q&=\frac{\nu s}{1 - \nu}=- \frac{21}{5}\quad[\mathrm{MPa}]\\
e_{c}&=\frac{- \frac{2 \nu^{2} s}{1 - \nu} + s}{1000 E}=\frac{27}{100000}\\
M&=\frac{E s}{- \frac{2 \nu^{2} s}{1 - \nu} + s}=\frac{700}{9}\quad[\mathrm{GPa}]
\end{aligned}
$$

横反力は圧縮に変わる。横拘束により軸ひずみの大きさが減ることは、Poisson 比が正のときと共通である。
:::

:::example[横拘束の試験が Poisson 比に敏感な範囲]
$E=70 \mathrm{GPa}$、$s=21 \mathrm{MPa}$ に固定し、Poisson 比 $\nu$ を $\nu=\frac{1}{3}$ の近くで変える。横拘束下の見かけの弾性率 $M$ の局所的な変化を求める。固定した値を先に代入すると、残る一変数の関数は次のようになる。以下の微分では数値を先の単位系で表す。

$$
\begin{aligned}
M(\nu)&=\frac{1470}{- \frac{42 \nu^{2}}{1 - \nu} + 21},\\
\frac{\mathrm dM}{\mathrm d\nu}&=- \frac{140 \nu \left(\nu - 2\right)}{\left(\nu + 1\right)^{2} \left(2 \nu - 1\right)^{2}}.
\end{aligned}
$$

基準値で $M=105$、微分は $\frac{1575}{4}$ である。したがって一次の展開は、基準値からの相対変化を $\nu\mapsto \nu(1+\delta)$ と書くと

$$
\frac{\Delta M}{M}
=\left(\frac{\nu}{M}
\frac{\mathrm dM}{\mathrm d\nu}\right)\delta+\order{\delta^2}
=\frac{5}{4}\,\delta+\order{\delta^2}
$$

となる。例えば入力を $1\%$ 増やすとき、一次式で得る出力の変化率は $\frac{5}{4}\%$ である。元の式に直接代入した変化率は $\frac{25100}{19649}\%$ なので、この比較で一次化の誤差を点検できる。

逆に、出力の測定値が基準値から $1\%$ 上がったとき、適用した模型が正しく、ほかの入力が固定されていれば、入力の変化率は一次で $\frac{4}{5}\%$ と推定される。見かけの係数は $\nu=0$ で停留する。その近くではこの係数一つから横ひずみの符号を決めにくい。
:::

## 弾性エネルギーと熱力学

ゆっくり負荷してひずみを $0$ から $\eps$ に沿って直線的に増やす。単位体積の仕事は $\int_0^1\sigma(s\eps):\eps\,\mathrm ds$ である。一般の線形弾性では下式の第 1 の等号を得る。等方体ではさらに第 2 の形に分解できる。

$$
W=\frac12 C_{ijkl}\eps_{ij}\eps_{kl}
=\mu\eps_{\mathrm d}:\eps_{\mathrm d}+\frac K2(\tr\eps)^2 .
\label{eq-cont-elastic-energy}
$$

弾性定数に主要対称性があれば $\mathrm dW=C_{ijkl}\eps_{kl}\mathrm d\eps_{ij}$ となる。従って経路を曲げても仕事は $W$ の差である。等方体で $W>0$ がすべての非零ひずみで成り立つ必要十分条件は $\mu>0,K>0$ である。必要性は跡ゼロと純膨張を別々に代入し、十分性は二項の非負性と直交分解から従う。

:::example[ひずみエネルギー密度]
$\mu=30\,\mathrm{GPa}$ で工学せん断ひずみ $\gamma=2\times10^{-3}$ を加える。$\eps_{12}=\eps_{21}=\gamma/2$ だから $W=\mu\gamma^2/2=6\times10^4\,\mathrm{J\,m^{-3}}$ である。
:::

温度変化 $\vartheta=T-T_0$ を許すとき、体積当たり自由エネルギーを $W_0(T)+\frac12(\eps-\alpha\vartheta I):C:(\eps-\alpha\vartheta I)$ とする。$\alpha$ は線膨張係数で単位は $\mathrm{K^{-1}}$ である。:eqref[eq-cont-dissipation] で可逆過程の散逸をゼロとし、密度を参照値で扱う一次の熱弾性では、応力はひずみでの偏微分、体積当たりエントロピーは温度での偏微分の負である。この模型では応力自由な温度変化で $\eps=\alpha\vartheta I$ が得られる。完全拘束では $\sigma=-3K\alpha\vartheta I$ なので加熱が圧縮応力を生む。混合微分の交換により、温度による応力変化とひずみによるエントロピー変化は符号を除き同じ係数を持つ。

### 温度を変える速さによって現れる弾性率

体積ひずみ $\theta$ と温度変化 $\vartheta$ に対する体積自由エネルギーを $\Psi=K_T\theta^2/2-\beta\theta\vartheta-c\vartheta^2/(2T_0)$ とする。$c$ は基準温度 $T_0$ における単位体積あたりの定ひずみ熱容量、$\beta$ は熱応力係数である。平均応力は $\partial_\theta\Psi=K_T\theta-\beta\vartheta$、基準からのエントロピー密度の変化は $-\partial_\vartheta\Psi=\beta\theta+c\vartheta/T_0$。等エントロピーなら $\vartheta=-\beta T_0\theta/c$ なので、平均応力の勾配は $K_S=K_T+\beta^2T_0/c$ となる。熱が出入りできない変形では、$\beta>0$ の材料では圧縮によって温度が上昇する。$\beta$ の符号によらず、温度変化の寄与は同じ体積ひずみに対する応力の絶対値を大きくする。

$K_T>0$、$c>0$ とする。このとき $K_S\ge K_T$ であり、温度と体積の結合がゼロの場合に等号である。温度の二次係数が負なのは、温度を自然変数に持つ自由エネルギーの曲率が熱容量の負号を含むためである。ひずみを所定の温度で増やす際の安定性は、ひずみに関する曲率 $K_T$ によって決まる。

:::example[断熱圧縮で上昇する温度]
$K_{T}=80 \mathrm{GPa}$、$\beta=3 \mathrm{MPa\,K^{-1}}$、$c=3 \mathrm{MJ\,m^{-3}\,K^{-1}}$、$T_{0}=300 \mathrm K$、$\theta=- \frac{1}{1000}$ とする。以下では$\Delta_{T}$ が断熱温度変化、$K_{S}$ が等エントロピー体積弾性率、$D$ が断熱と等温の平均応力の差である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
\Delta_{T}&=- \frac{T_{0} \beta \theta}{c}=\frac{3}{10}\quad[\mathrm K]\\
K_{S}&=K_{T} + \frac{T_{0} \beta^{2}}{1000 c}=\frac{809}{10}\quad[\mathrm{GPa}]\\
D&=\frac{T_{0} \beta^{2} \theta}{c}=- \frac{9}{10}\quad[\mathrm{MPa}]
\end{aligned}
$$

平均応力の差は負である。断熱圧縮に必要な外圧が、同じ体積変化をさせる等温圧縮より大きい。
:::

:::example[結合が弱い材料での補正]
$K_{T}=80 \mathrm{GPa}$、$\beta=1 \mathrm{MPa\,K^{-1}}$、$c=2 \mathrm{MJ\,m^{-3}\,K^{-1}}$、$T_{0}=600 \mathrm K$、$\theta=- \frac{1}{1000}$ とする。以下では$\Delta_{T}$ が断熱温度変化、$K_{S}$ が等エントロピー体積弾性率、$D$ が断熱と等温の平均応力の差である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
\Delta_{T}&=- \frac{T_{0} \beta \theta}{c}=\frac{3}{10}\quad[\mathrm K]\\
K_{S}&=K_{T} + \frac{T_{0} \beta^{2}}{1000 c}=\frac{803}{10}\quad[\mathrm{GPa}]\\
D&=\frac{T_{0} \beta^{2} \theta}{c}=- \frac{3}{10}\quad[\mathrm{MPa}]
\end{aligned}
$$

初期温度を高くしても、熱応力係数の二乗が効くため弾性率の差は先の例より小さい。
:::

:::example[熱応力係数が圧縮仕事に与える影響]
$K_{T}=80 \mathrm{GPa}$、$c=3 \mathrm{MJ\,m^{-3}\,K^{-1}}$、$T_{0}=300 \mathrm K$、$\theta=- \frac{1}{1000}$ に固定し、熱応力係数 $\beta$ を $\beta=3 \mathrm{MPa\,K^{-1}}$ の近くで変える。断熱と等温の平均応力の差 $D$ の局所的な変化を求める。固定した値を先に代入すると、残る一変数の関数は次のようになる。以下の微分では数値を先の単位系で表す。

$$
\begin{aligned}
D(\beta)&=- \frac{\beta^{2}}{10},\\
\frac{\mathrm dD}{\mathrm d\beta}&=- \frac{\beta}{5}.
\end{aligned}
$$

基準値で $D=- \frac{9}{10}$、微分は $- \frac{3}{5}$ である。したがって一次の展開は、基準値からの相対変化を $\beta\mapsto \beta(1+\delta)$ と書くと

$$
\frac{\Delta D}{D}
=\left(\frac{\beta}{D}
\frac{\mathrm dD}{\mathrm d\beta}\right)\delta+\order{\delta^2}
=2\,\delta+\order{\delta^2}
$$

となる。例えば入力を $1\%$ 増やすとき、一次式で得る出力の変化率は $2\%$ である。元の式に直接代入した変化率は $\frac{201}{100}\%$ なので、この比較で一次化の誤差を点検できる。

逆に、出力の測定値が基準値から $1\%$ 上がったとき、適用した模型が正しく、ほかの入力が固定されていれば、入力の変化率は一次で $\frac{1}{2}\%$ と推定される。断熱効果だけを求めるには、等温状態と断熱状態の応力の差を取る。この差の測定では、両試験で同じ体積変化を与える。
:::

## Navier の方程式

:::proposition[変位の方程式]{#prop-cont-navier}
等方均質な線形弾性体で二つの Lamé 定数が一定なら、微小運動の式は

$$
\rho\,\partial_t^2\vb u=\mu\Delta\vb u+(\lambda+\mu)\grad(\divg\vb u)+\rho\vb b
\label{eq-cont-navier}
$$

となる。静的平衡では左辺がゼロである。
:::

:::proof
:ref[law-cont-momentum] の加速度を参照座標における $\partial_t^2\vb u$ で一次近似する。$\sigma_{ij}=\lambda\delta_{ij}\partial_ku_k+\mu(\partial_i u_j+\partial_j u_i)$ を $\partial_j\sigma_{ij}$ に代入する。定数を微分の外へ出し、混合微分を入れ替えると式を得る。
:::

境界を変位指定部分 $\Gamma_D$ と面力指定部分 $\Gamma_N$ に分ける。後者には $\sigma\vb n=\bar{\vb t}$ を課す。静的平衡はポテンシャル

$$
\Pi[\vb u]=\int_B W\,\mathrm dv-\int_B\rho\vb b\cdot\vb u\,\mathrm dv
-\int_{\Gamma_N}\bar{\vb t}\cdot\vb u\,\mathrm da
$$

を変位指定を満たす範囲で停留させる条件でもある。変分 $\vb h$ は $\Gamma_D$ でゼロであり、$\delta W=\sigma:\grad\vb h$ である。部分積分により $\delta\Pi=-\int_B(\divg\sigma+\rho\vb b)\cdot\vb h+\int_{\Gamma_N}(\sigma\vb n-\bar{\vb t})\cdot\vb h$ となる。内部支持の $\vb h$ と境界に任意の値を持つ $\vb h$ を順に選べば平衡式と自然境界条件が出る。

二つの平衡解の差 $\vb w$ は斉次境界条件を満たす。差の式を $\vb w$ に掛けて部分積分すると $\int_B 2W(\eps(\vb w))=0$ になる。$W$ が正定値ならひずみゼロであり、:eqref[eq-cont-displacement-hessian] より差は微小剛体運動に限られる。境界がこの運動を拘束すれば一意である。純粋な面力指定では、並進と微小回転を仮想変位に入れると外力の全合力と全トルクがゼロであることが必要と分かる。

:::example[体積力と境界拘束]
板 $0<x<L$ が横に拘束されて $\vb u=(u(x),0,0)$、単位質量当たりの力 $b_1$ が一定とする。$u(0)=0$、$\sigma_{11}(L)=0$ より $(\lambda+2\mu)u''=-\rho b_1$、$u'(L)=0$ である。従って $u=\rho b_1(Lx-x^2/2)/(\lambda+2\mu)$ となる。自由端でも体積力があれば変位はゼロにならない。
:::

### 二次元平衡を Airy の応力関数で解く

体積力なしの平面平衡では $\partial_x\sigma_{xx}+\partial_y\sigma_{xy}=0$、$\partial_x\sigma_{xy}+\partial_y\sigma_{yy}=0$ である。単連結領域では各行が発散ゼロなので、$f,g$ を用いて $(\sigma_{xx},\sigma_{xy})=(f_y,-f_x)$、$(\sigma_{xy},\sigma_{yy})=(g_y,-g_x)$ と書ける。応力対称性が $f_x+g_y=0$ を要求するので、さらに $f=\Phi_y,g=-\Phi_x$ と表せる。従って :term[Airy の応力関数]{en="Airy stress function" yomi="えありーのおうりょくかんすう"} $\Phi$ により

$$
\sigma_{xx}=\Phi_{yy},\qquad \sigma_{yy}=\Phi_{xx},\qquad
\sigma_{xy}=-\Phi_{xy}
\label{eq-cont-airy-stress}
$$

とできる。これらは候補の応力について平衡を自動で保証する。

平面応力では $\eps_{xx}=(\Phi_{yy}-\nu\Phi_{xx})/E$、$\eps_{yy}=(\Phi_{xx}-\nu\Phi_{yy})/E$、$\eps_{xy}=-(1+\nu)\Phi_{xy}/E$。ひずみの定義に代入して混合微分を交換すると、変位から得られるひずみは $\partial_{yy}\eps_{xx}+\partial_{xx}\eps_{yy}-2\partial_{xy}\eps_{xy}=0$ を満たす。応力関数でこの式を評価すると、混合微分の係数が $-2\nu+2(1+\nu)=2$ となるので $\Delta^2\Phi=0$ である。平面ひずみの場合は $\sigma_{zz}=\nu(\sigma_{xx}+\sigma_{yy})$ を入れ、面内の応力–ひずみ式を $E/(1-\nu^2)$ と $\nu/(1-\nu)$ の平面応力型に書き直せるため、同じ双調和式となる。ただし厚さ方向の応力と変位は両問題で異なる。

:::example[遠方で引っ張られた無限板の円孔]
半径 $a$ の円孔を持つ平面応力板を $x$ 方向に遠方応力 $S$ で引く。極座標で無孔板の応力関数 $Sy^2/2=Sr^2(1-\cos 2\theta)/4$ に、遠方応力を変えない軸対称項 $\log r$ と二次角モードの項 $1,r^{-2}$ を加える。二次角モードの Laplacian は $\Delta(r^m\cos 2\theta)=(m^2-4)r^{m-2}\cos 2\theta$ なので、$m=2,0,-2$ は双調和である。$\log r$ も孔外で調和である。極座標に変換した応力は $\sigma_{rr}=\Phi_r/r+\Phi_{\theta\theta}/r^2$、$\sigma_{r\theta}=-\partial_r(\Phi_\theta/r)$、$\sigma_{\theta\theta}=\Phi_{rr}$。これを $r=a$ の自由条件 $\sigma_{rr}=\sigma_{r\theta}=0$ に入れると

$$
\Phi=\frac S4\left[r^2-2a^2\log(r/a)-
(r^2-2a^2+a^4/r^2)\cos 2\theta\right]
$$

となる。実際、$r=a$ で軸対称項の半径微分がゼロ、角モードの括弧とその半径微分もゼロとなる。周応力だけが残り $\sigma_{\theta\theta}(a,\theta)=S(1-2\cos 2\theta)$。孔の上下端では $3S$、左右端では $-S$ であり、遠方の一軸荷重でも局所には圧縮領域と応力集中が生じる。

孔の外は単連結でないので、変位が一価になることも確認する。平面応力則から得た半径ひずみを $r$ で積分し、周ひずみとせん断ひずみに合わせると

$$
\begin{aligned}
u_r&=\frac S{2E}\left[(1-\nu)r+\frac{(1+\nu)a^2}{r}
+\left((1+\nu)r+\frac{4a^2}{r}-\frac{(1+\nu)a^4}{r^3}\right)\cos2\theta\right],\\
u_\theta&=-\frac S{2E}\left[(1+\nu)r+\frac{2(1-\nu)a^2}{r}+\frac{(1+\nu)a^4}{r^3}\right]\sin2\theta .
\end{aligned}
$$

$\eps_{rr}=\partial_r u_r$、$\eps_{\theta\theta}=(u_r+\partial_\theta u_\theta)/r$、$2\eps_{r\theta}=\partial_r u_\theta-u_\theta/r+(\partial_\theta u_r)/r$ を計算すれば、上の応力に対する :eqref[eq-cont-compliance] を満たす。この変位は $\theta$ に関して $2\pi$ 周期であり、一周した後も同じ値に戻る。
:::

### 応力成分で書いた三次元の適合式

均質等方体で体積力がなく、応力が十分滑らかとする。:eqref[eq-cont-displacement-hessian] で二つの微分添字を縮約すると $\Delta u_i=2\partial_j\eps_{ij}-\partial_i\tr\eps$。両辺の対称勾配を取って

$$
\Delta\eps_{ij}=\partial_i\partial_k\eps_{jk}
+\partial_j\partial_k\eps_{ik}-\partial_i\partial_j\tr\eps
$$

を得る。ここで :eqref[eq-cont-compliance] と $\partial_j\sigma_{ij}=0$ を使うと、$\Delta\sigma_{ij}$ は $\tr\sigma$ の二階微分で表せる。跡を取って $\Delta\tr\sigma=0$ を得た上で、その項を消去すると

$$
\Delta\sigma_{ij}+\frac{1}{1+\nu}\partial_i\partial_j\tr\sigma=0
$$

となる。これは体積力のない場合の :term[Beltrami–Michell の適合式]{en="Beltrami–Michell compatibility equations" yomi="べるとらみみっちぇるのてきごうしき"} である。応力を未知量にする場合は平衡式と組み合わせる。穴のある領域では、局所式を満たす候補についても変位の一価性を確かめる。

### 弾性の相反仕事

同じ材料と領域で、二つの荷重がそれぞれ平衡変位 $\vb u^{(1)},\vb u^{(2)}$ を与えるとする。どちらも変位指定部でゼロである。第 1 の平衡式に第 2 の変位を掛けて部分積分すると

$$
\int_B \rho\vb b^{(1)}\cdot\vb u^{(2)}\,\mathrm dv
+\int_{\Gamma_N}\vb t^{(1)}\cdot\vb u^{(2)}\,\mathrm da
=\int_B\eps^{(2)}:C:\eps^{(1)}\,\mathrm dv .
$$

主要対称性 $C_{ijkl}=C_{klij}$ で右辺は上添字を交換しても同じである。この等式を Betti の相反定理という。例えば物体の部位 $A$ を押したときに $B$ が力の方向へどれだけ動くかを求めるには、$B$ を単位荷重で押したときの $A$ の変位を代わりに測れる。力が一点で特異になるときは、面積の小さい荷重分布で得た量を極限として扱う。

反対に、線形の応力応答で主要対称性が失われると、変位を増やす順番で仕事が変わる。たとえば二成分 $q_1,q_2$ のひずみを用い、対応する応力を $f_1=a q_1+b q_2$、$f_2=c q_1+d q_2$ とする。原点から $q_1=A$、次に $q_2=B$ まで進む仕事は $aA^2/2+cAB+dB^2/2$。順番を交換すれば交差項は $bAB$ になる。差が任意の $A,B$ で消えるには $b=c$ が必要である。相反性は可逆なひずみエネルギーの存在と結び付く。

### 仮想仕事と混合境界条件

固定境界 $\Gamma_D$ で $\vb u$ を、残りの $\Gamma_N$ で面力 $\vb t$ を与える。十分滑らかな静的平衡解について、$\Gamma_D$ でゼロの試験変位 $\vb w$ を用い、運動方程式に $\vb w$ を掛ける。$\partial_j(\sigma_{ij}w_i)=w_i\partial_j\sigma_{ij}+\sigma_{ij}\partial_jw_i$ と応力の対称性から

$$
\int_B\eps(\vb w):\sigma(\vb u)\,\dd V
=\int_B\rho\vb b\cdot\vb w\,\dd V+\int_{\Gamma_N}\vb t\cdot\vb w\,\dd A.
$$

逆にこの等式をすべての試験変位について仮定し、まず内部だけに台を持つ試験変位を取れば内部の平衡式が従う。次に $\Gamma_N$ の小領域で任意の境界値を持つ試験変位を取ると、指定面力の条件が従う。剛体並進と回転が試験変位として許される全面力境界では、左辺がゼロとなるため合力・合モーメントがゼロでなければならない。

一次元の模型で左固定、右端力 $P$、長さ方向に一様な体積力による分布荷重 $q$ を考える。断面力は $N'+q=0$、$N(L)=P$ から $N=P+q(L-x)$。$u'=N/(EA)$ を積分すれば $u(L)=[PL+qL^2/2]/EA$。同じ $N$ からエネルギー $U=\int N^2/(2EA)\,\dd x$ を求め、$P$ で微分すると $u(L)$ が出る。この関係は、追加の端荷重に対する仮想仕事を表す。

:::example[端と全長に引張荷重が働く棒]
$P=10 \mathrm N$、$q=4 \mathrm{N\,m^{-1}}$、$L=2 \mathrm m$、$A_{E}=20000 \mathrm N$ とする。以下では$N_{0}$ が固定端断面力、$U$ がひずみエネルギー、$u_{L}$ が自由端の変位である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
N_{0}&=L q + P=18\quad[\mathrm N]\\
U&=\frac{\frac{L^{3} q^{2}}{3} + L^{2} P q + L P^{2}}{2 A_{E}}=\frac{151}{15000}\quad[\mathrm J]\\
u_{L}&=\frac{\frac{L^{2} q}{2} + L P}{A_{E}}=\frac{7}{5000}\quad[\mathrm m]
\end{aligned}
$$

固定端で断面力が最大になる。支持反力の向きはその断面力による合力と逆である。
:::

:::example[一部の荷重が引張を緩める棒]
$P=20 \mathrm N$、$q=-4 \mathrm{N\,m^{-1}}$、$L=2 \mathrm m$、$A_{E}=20000 \mathrm N$ とする。以下では$N_{0}$ が固定端断面力、$U$ がひずみエネルギー、$u_{L}$ が自由端の変位である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
N_{0}&=L q + P=12\quad[\mathrm N]\\
U&=\frac{\frac{L^{3} q^{2}}{3} + L^{2} P q + L P^{2}}{2 A_{E}}=\frac{49}{3750}\quad[\mathrm J]\\
u_{L}&=\frac{\frac{L^{2} q}{2} + L P}{A_{E}}=\frac{1}{625}\quad[\mathrm m]
\end{aligned}
$$

逆向きの分布荷重で固定端側の断面力が減る。エネルギーは断面力を二乗するため、合荷重だけを使う計算とは一致しない。
:::

:::example[端変位を荷重の検定に使う]
$q=4 \mathrm{N\,m^{-1}}$、$L=2 \mathrm m$、$A_{E}=20000 \mathrm N$ に固定し、端荷重 $P$ を $P=10 \mathrm N$ の近くで変える。自由端の変位 $u_{L}$ の局所的な変化を求める。固定した値を先に代入すると、残る一変数の関数は次のようになる。以下の微分では数値を先の単位系で表す。

$$
\begin{aligned}
u_{L}(P)&=\frac{P}{10000} + \frac{1}{2500},\\
\frac{\mathrm du_{L}}{\mathrm dP}&=\frac{1}{10000}.
\end{aligned}
$$

基準値で $u_{L}=\frac{7}{5000}$、微分は $\frac{1}{10000}$ である。したがって一次の展開は、基準値からの相対変化を $P\mapsto P(1+\delta)$ と書くと

$$
\frac{\Delta u_{L}}{u_{L}}
=\left(\frac{P}{u_{L}}
\frac{\mathrm du_{L}}{\mathrm dP}\right)\delta+\order{\delta^2}
=\frac{5}{7}\,\delta+\order{\delta^2}
$$

となる。例えば入力を $1\%$ 増やすとき、一次式で得る出力の変化率は $\frac{5}{7}\%$ である。元の式に直接代入した変化率は $\frac{5}{7}\%$ なので、この比較で一次化の誤差を点検できる。

逆に、出力の測定値が基準値から $1\%$ 上がったとき、適用した模型が正しく、ほかの入力が固定されていれば、入力の変化率は一次で $\frac{7}{5}\%$ と推定される。比例係数の絶対値は伸びコンプライアンス $L/(EA)$ である。分布荷重が残っているので端変位は一般に原点を通らない。
:::

### 平衡解の一意性と全ポテンシャルの最小化

静止解 $\vb u$ に、指定変位境界で消える変位 $\vb w$ を加える。弾性テンソルの主要対称性によりエネルギーの増分は $\int_B\eps(\vb w):C\eps(\vb u)+\tfrac12\eps(\vb w):C\eps(\vb w)\,\dd V$。外力を所定の体積力と面力とした全ポテンシャルでは、前者の線形項が仮想仕事と相殺する。したがって指定条件を満たす全変位で全ポテンシャルは平衡解の値以上になる。逆に最小なら、任意の $\vb w$ に対し $\vb u+\eta\vb w$ の $\eta=0$ での一次係数がゼロであるため仮想仕事式を得る。

:::theorem[線形弾性の平衡解の一意性]{#thm-cont-static-uniqueness}
連結領域で弾性エネルギーが各点で正定値であり、体積力と面力が変位に依存しないとする。同じ変位・面力境界条件を満たす二つの滑らかな平衡解は、微小剛体運動を除いて一致する。指定変位境界に平面であっても開いた面片が含まれるなら、両解は一致する。
:::

:::proof
二解の差 $\vb w$ は同次境界条件と体積力ゼロの平衡式を満たす。その解自身を試験変位に取ると、$\int_B\eps(\vb w):C\eps(\vb w)\,\dd V=0$。正定値性と連続性から $\eps(\vb w)=0$ である。:eqref[eq-cont-displacement-hessian] で右辺がゼロになるので $\vb w=\vb a+B\vb x$、$B^T=-B$。指定変位面片上の互いに独立な二つの差ベクトルに $B$ を作用させるとゼロになる。三次元の非零反対称行列の核は一次元なので $B=0$、そして面上一点の値から $\vb a=0$ である。
:::

### 断面を離散化したときの厳密解との比較

一次元棒の変位を区間 $0<x<L$ で二次多項式に制限する。左端を固定し、右に集中力 $P$、全長に一様な分布荷重 $q$ を作用させる。$u=ax+bx^2$ とすると平衡式は $2EA b=-q$、端の力は $EA(a+2bL)=P$。曲率が一定なので、この形で厳密解を再現できる。線形の試験変位しか許さないなら $u=ax$。仮想変位を $x$ として $\int_0^LEA a\,\dd x=PL+qL^2/2$。よって自由端の変位は偶然にも厳密解と等しいが、内部の伸びの分布は異なる。

エネルギー誤差を求めるには勾配の差 $[q(L/2-x)]/(EA)$ を二乗して $\int_0^L EA/2$ の重みで積分する。$\int_0^L(L/2-x)^2\,\dd x=L^3/12$ から差は $q^2L^3/(24EA)$ となる。自由端値が合っても全ポテンシャルはこの量だけ大きい。分布荷重がなければ差はゼロになり、一次変位でも厳密解を再現する。

## 平衡問題の例（一様圧縮・球殻・円柱のねじり）

### 静水圧による一様圧縮

境界全体で $\vb t=-p\vb n$ を与え、体積力はないとする。定数応力 $-pI$ は発散ゼロで境界条件を満たす。対応する変位は $\vb u=-p\vb x/(3K)$ である。例えば $K=60\,\mathrm{GPa}$、$p=6\,\mathrm{MPa}$ なら体積変化率は $-10^{-4}$、各長さの変化率は $-10^{-4}/3$ である。

### 厚肉球殻

内半径 $a$、外半径 $b$ の球殻で $\vb u=f(r)\vb x/r$ と置く。$n_i=x_i/r$ について $\partial_j n_i=(\delta_{ij}-n_in_j)/r$ なので、主ひずみは $f',f/r,f/r$ である。応力は $\sigma=s_r\vb n\vb n^T+s_t(I-\vb n\vb n^T)$ と書け、微分すると発散の半径成分は $s_r'+2(s_r-s_t)/r$ である。構成式を入れると $(\lambda+2\mu)(f''+2f'/r-2f/r^2)=0$。$r^m$ を代入して $(m-1)(m+2)=0$ だから

$$
f=Ar+B/r^2,\quad
s_r=3KA-4\mu B/r^3,\quad s_t=3KA+2\mu B/r^3 .
$$

内外の圧力を $p_i,p_o$ とし、$s_r(a)=-p_i$、$s_r(b)=-p_o$ とする。引き算と代入により

$$
4\mu B=\frac{a^3b^3(p_i-p_o)}{b^3-a^3},\qquad
3KA=\frac{a^3p_i-b^3p_o}{b^3-a^3}.
$$

$a=1\,\mathrm m$、$b=2\,\mathrm m$、$p_o=0$ なら $s_t(a)=5p_i/7$、$s_t(b)=3p_i/14$ である。内側の表面は半径方向に圧縮されながら周方向に伸ばされる。

### 円柱のねじり

半径 $a$ の円柱で $\vb u=\vartheta z(-y,x,0)$ を試す。$\vartheta$ は単位長さ当たりのねじれである。非零応力は $\sigma_{xz}=-\mu\vartheta y$、$\sigma_{yz}=\mu\vartheta x$ と対称成分であり、発散はゼロである。側面の法線 $(x,y,0)/a$ では軸方向面力 $(-\mu\vartheta yx+\mu\vartheta xy)/a$ が消える。端面トルクは

$$
M=\int (x\sigma_{yz}-y\sigma_{xz})\,\mathrm dA
=\mu\vartheta\int r^2\,\mathrm dA=\mu\vartheta\pi a^4/2 .
\label{eq-cont-circular-torsion}
$$

最後は $\int_0^a 2\pi r^3\,\mathrm dr$ を評価した。$a=0.02\,\mathrm m$、$\mu=80\,\mathrm{GPa}$、$\vartheta=0.01\,\mathrm{m^{-1}}$ なら $M=64\pi\,\mathrm{N\,m}$ である。

### 無限弾性体の集中力

均質等方で静的な無限空間に、原点へ合力 $\vb f$ を加えた理想化を考える。力の密度を $\vb f\delta(\vb x)$ と書き、$r=|\vb x|$ とする。静的 Navier 式を $-\mu\Delta\vb u-(\lambda+\mu)\grad\divg\vb u=\vb f\delta$ の形で使う。スカラー関数 $G=1/(4\pi r)$ は $-\Delta G=\delta$ を満たす。原点外で調和であり、原点を囲む球の法線微分の負の積分が 1 であるためである。$H=-r/(8\pi)$ は $-\Delta H=G$ を満たす。原点外で $\Delta r=2/r$ であり、原点を囲む球の $H$ の微分の積分は面積に比例してゼロへ行くので、後者には球中心のデルタ項がない。

解の候補を $\vb u=G\vb f/\mu+b\grad(\vb f\cdot\grad H)$ と置く。第一項の Laplacian が集中力そのものを作る。残る勾配成分は、$\divg\vb u=(1/\mu-b)\vb f\cdot\grad G$ と $-\Delta H=G$ を使うと係数が $\mu b-(\lambda+\mu)(1/\mu-b)$ となる。ゼロにするには $b=(\lambda+\mu)/[\mu(\lambda+2\mu)]$ である。$H$ の二階微分は $-(\delta_{ij}/r-x_ix_j/r^3)/(8\pi)$ なので、Poisson 比に書き直して

$$
u_i(\vb x)=\frac{1}{16\pi\mu(1-\nu)}
\left[(3-4\nu)\frac{\delta_{ij}}{r}+\frac{x_ix_j}{r^3}\right]f_j
$$

を得る。これは Kelvin の基本解である。有限の空間分布の力に対しては右辺の核を重ね合わせると、境界のない問題の変位が得られる。境界のある問題では境界条件を満たす斉次解をさらに足す。

:::example[力の正面と側方における変位]
$\nu=1/4$、$\vb f=f\vb e_1$ とする。同じ半径で $\vb x=r\vb e_1$ では力方向の変位に係数 3 が掛かる。$\vb x=r\vb e_2$ では係数は 2 である。従って変位は両方向で力と平行だが正面の値は側方の $3/2$ 倍になる。一般の斜めの位置では、$x_i x_j f_j$ の項により力と平行でない成分も現れる。
:::

### 厚肉円筒の内圧問題

一様な厚肉管の端から十分遠い範囲で、軸方向ひずみをゼロとする平面ひずみを仮定する。半径変位 $u(r)$ のひずみは $u'$ と $u/r$。$\sigma_{rr}=(\lambda+2\mu)u'+\lambda u/r$、$\sigma_{\theta\theta}=\lambda u'+(\lambda+2\mu)u/r$ を平衡式 $\sigma_{rr}'+(\sigma_{rr}-\sigma_{\theta\theta})/r=0$ に代入すると $u''+u'/r-u/r^2=0$ になる。$u=r^m$ を代入すると $m^2-1=0$ なので $u=Cr+D/r$。これを応力式に入れれば $\sigma_{rr}=A-B/r^2$、$\sigma_{\theta\theta}=A+B/r^2$ となる。$A=2(\lambda+\mu)C$、$B=2\mu D$ である。

内半径 $a$ の圧力 $p_i$、外半径 $b$ の圧力 $p_o$ に対し $\sigma_{rr}(a)=-p_i$、$\sigma_{rr}(b)=-p_o$ とおく。二式を引くと $B=(p_i-p_o)a^2b^2/(b^2-a^2)$、定数項は $A=(p_i a^2-p_o b^2)/(b^2-a^2)$ である。周応力が最大になるのは $B>0$ で内面である。薄い極限 $b=a+h$ では $b^2-a^2=h(2a+h)$ より、外圧ゼロの内面周応力は $p_i a/h$ に近づく。球殻の結果とは曲率方向の数が異なる。

:::example[外圧を受けない厚肉管]
$p_{i}=10 \mathrm{MPa}$、$p_{o}=0 \mathrm{MPa}$、$a=\frac{1}{10} \mathrm m$、$b=\frac{1}{5} \mathrm m$ とする。以下では$A$ が一定応力項、$B$ が半径依存項の係数、$s_{o}$ が外面周応力、$s_{i}$ が内面周応力である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
A&=\frac{a^{2} p_{i} - b^{2} p_{o}}{- a^{2} + b^{2}}=\frac{10}{3}\quad[\mathrm{MPa}]\\
B&=\frac{a^{2} b^{2} \left(p_{i} - p_{o}\right)}{- a^{2} + b^{2}}=\frac{2}{15}\quad[\mathrm{MPa\,m^2}]\\
s_{o}&=\frac{a^{2} \left(p_{i} - p_{o}\right)}{- a^{2} + b^{2}} + \frac{a^{2} p_{i} - b^{2} p_{o}}{- a^{2} + b^{2}}=\frac{20}{3}\quad[\mathrm{MPa}]\\
s_{i}&=\frac{b^{2} \left(p_{i} - p_{o}\right)}{- a^{2} + b^{2}} + \frac{a^{2} p_{i} - b^{2} p_{o}}{- a^{2} + b^{2}}=\frac{50}{3}\quad[\mathrm{MPa}]
\end{aligned}
$$

周応力は内側の方が大きい。半径応力と周応力の双方を比較して、どの成分を強度評価に使うか決める。
:::

:::example[外圧を加えた管]
$p_{i}=10 \mathrm{MPa}$、$p_{o}=4 \mathrm{MPa}$、$a=\frac{1}{10} \mathrm m$、$b=\frac{3}{20} \mathrm m$ とする。以下では$A$ が一定応力項、$B$ が半径依存項の係数、$s_{o}$ が外面周応力、$s_{i}$ が内面周応力である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
A&=\frac{a^{2} p_{i} - b^{2} p_{o}}{- a^{2} + b^{2}}=\frac{4}{5}\quad[\mathrm{MPa}]\\
B&=\frac{a^{2} b^{2} \left(p_{i} - p_{o}\right)}{- a^{2} + b^{2}}=\frac{27}{250}\quad[\mathrm{MPa\,m^2}]\\
s_{o}&=\frac{a^{2} \left(p_{i} - p_{o}\right)}{- a^{2} + b^{2}} + \frac{a^{2} p_{i} - b^{2} p_{o}}{- a^{2} + b^{2}}=\frac{28}{5}\quad[\mathrm{MPa}]\\
s_{i}&=\frac{b^{2} \left(p_{i} - p_{o}\right)}{- a^{2} + b^{2}} + \frac{a^{2} p_{i} - b^{2} p_{o}}{- a^{2} + b^{2}}=\frac{58}{5}\quad[\mathrm{MPa}]
\end{aligned}
$$

内外の差圧が小さくなると勾配が緩くなる。内外に同じ圧力を加えても円筒内の垂直応力の基準値は移動する。
:::

:::example[厚さを変更したときの最大周応力]
$p_{i}=10 \mathrm{MPa}$、$p_{o}=0 \mathrm{MPa}$、$a=\frac{1}{10} \mathrm m$ に固定し、外半径 $b$ を $b=\frac{1}{5} \mathrm m$ の近くで変える。内面周応力 $s_{i}$ の局所的な変化を求める。固定した値を先に代入すると、残る一変数の関数は次のようになる。以下の微分では数値を先の単位系で表す。

$$
\begin{aligned}
s_{i}(b)&=\frac{10 b^{2}}{b^{2} - \frac{1}{100}} + \frac{1}{10 \left(b^{2} - \frac{1}{100}\right)},\\
\frac{\mathrm ds_{i}}{\mathrm db}&=- \frac{4000 b}{\left(10 b - 1\right)^{2} \left(10 b + 1\right)^{2}}.
\end{aligned}
$$

基準値で $s_{i}=\frac{50}{3}$、微分は $- \frac{800}{9}$ である。したがって一次の展開は、基準値からの相対変化を $b\mapsto b(1+\delta)$ と書くと

$$
\frac{\Delta s_{i}}{s_{i}}
=\left(\frac{b}{s_{i}}
\frac{\mathrm ds_{i}}{\mathrm db}\right)\delta+\order{\delta^2}
=- \frac{16}{15}\,\delta+\order{\delta^2}
$$

となる。例えば入力を $1\%$ 増やすとき、一次式で得る出力の変化率は $- \frac{16}{15}\%$ である。元の式に直接代入した変化率は $- \frac{2680}{2567}\%$ なので、この比較で一次化の誤差を点検できる。

逆に、出力の測定値が基準値から $1\%$ 上がったとき、適用した模型が正しく、ほかの入力が固定されていれば、入力の変化率は一次で $- \frac{15}{16}\%$ と推定される。内半径と圧力を固定している。同じ厚さの差でも外半径と内半径の比が近いときほど応力の相対変化が大きい。
:::

## 適合条件と Saint-Venant の原理

ひずみの六成分は三成分の変位から作るため、任意に与えた六つの関数が変位を持つとは限らない。

:::theorem[微小ひずみの適合条件]{#thm-cont-compatibility}
単連結な開領域で対称な $C^2$ 場 $\eps$ が微小ひずみとなる必要十分条件は、すべての添字について

$$
\partial_l\partial_j\eps_{ik}+\partial_k\partial_i\eps_{jl}
-\partial_k\partial_j\eps_{il}-\partial_l\partial_i\eps_{jk}=0
\label{eq-cont-compatibility}
$$

である。
:::

:::proof
$A_{ijk}=\partial_k\eps_{ij}+\partial_j\eps_{ik}-\partial_i\eps_{jk}$ と置く。必要性は :eqref[eq-cont-displacement-hessian] から $A_{ijk}=\partial_k\partial_j u_i$ となり $\partial_l A_{ijk}=\partial_k A_{ijl}$ が従うことである。この差を展開すると :eqref[eq-cont-compatibility] になる。

逆に等式を仮定すると、各 $i,j$ で $k$ 成分が $A_{ijk}$ の場は回転ゼロである。単連結性と Stokes の定理により線積分は経路によらず、$\partial_k H_{ij}=A_{ijk}$ を満たす $H_{ij}$ が得られる。$A_{ijk}=A_{ikj}$ なので $H_{ij}$ は各 $i$ について同様に回転ゼロであり、$H_{ij}=\partial_j u_i$ と積分できる。$\partial_k[(H_{ij}+H_{ji})/2-\eps_{ij}]=0$ なので差は定数対称行列である。この行列を勾配に持つ線形変位を $\vb u$ から引けば所望のひずみを得る。
:::

平面の独立条件は $\partial_y^2\eps_{xx}+\partial_x^2\eps_{yy}=2\partial_x\partial_y\eps_{xy}$ である。例えば $\eps_{xx}=cy^2$、他はゼロなら左辺は $2c$ なので $c\ne0$ で不適合である。これに対して $\eps_{xy}=cxy$ を加えると適合し、実際 $u_x=cxy^2,u_y=0$ がひずみを再現する。

### 荷重の細部が遠方で減衰する例

:term[Saint-Venant の原理]{en="Saint-Venant principle" yomi="さんぶなんのげんり"} は、合力と合モーメントを等しくした局所荷重の差が、負荷領域から離れるにつれて小さくなるという弾性体の性質を指す。断面形状と境界条件が減衰率を決める。一般の有限領域で一律の減衰率を与える主張ではない。半無限の一様な柱では、断面形状と側面条件を指定して断面の弾性固有値を求めることで、端からの減衰率を評価する。この節では反平面変形に限定し、断面の固有モードを陽に求める。

半無限帯 $x>0,0<y<h$ が $z$ 方向に一様で、$u_z=w(x,y)$ のみを持つとする。$\Delta w=0$、側面で $w_y=0$、端面で $-\mu w_x=g_0\cos(2\pi y/h)$ とする。$k=2\pi/h$ と書けば $w=(g_0/\mu k)e^{-kx}\cos(ky)$ は全条件を満たす。端面の積分はゼロであり、$\int_0^h y\cos(2\pi y/h)\,\mathrm dy=0$ なのでモーメントもゼロである。この端面荷重の応力は距離 $h$ で $e^{-2\pi}$ 倍になる。これに対して一定面力に対する $w=-g_0x/\mu$ では応力は全長にわたって一定である。

### 温度場が自由膨張の適合を失う条件

平面内の自由熱ひずみを $\eps_{xx}=\eps_{yy}=\alpha T(x,y)$、$\eps_{xy}=0$ とする。二次元の適合式 $\partial_{yy}\eps_{xx}+\partial_{xx}\eps_{yy}-2\partial_{xy}\eps_{xy}=0$ を用いると、$\alpha$ が定数なら $\Delta T=0$ が必要である。単連結領域で滑らかな場合はこの条件で十分である。実際、調和関数 $\alpha T$ の調和共役 $g$ を取り、正則関数 $\alpha T+ig$ の原始関数を $u_x+iu_y$ とする。Cauchy–Riemann の関係から $\partial_xu_x=\partial_yu_y=\alpha T$、$\partial_yu_x+\partial_xu_y=0$ が得られる。この構成は二次元の自由熱ひずみを扱い、厚さ方向の変位の適合は別の問題である。

薄い円板の平面応力で $T=T_0+A r^2$ を考える。軸対称変位 $u(r)$ を仮定し、$\sigma_{rr}=E[u'+\nu u/r-(1+\nu)\alpha T]/(1-\nu^2)$ と周方向の同様の式を平衡式へ入れると $u''+u'/r-u/r^2=(1+\nu)\alpha T'$。$r^3$ に左辺の作用素を作用させると $8r$ となるため、中心で有限な解は $u=Cr+(1+\nu)\alpha Ar^3/4$。縁 $r=a$ で $\sigma_{rr}=0$ を課して $C=\alpha T_0+(1-\nu)\alpha Aa^2/4$ を得る。代入すると $\sigma_{rr}=E\alpha A(a^2-r^2)/4$、$\sigma_{\theta\theta}=E\alpha A(a^2-3r^2)/4$ である。

:::example[高温の縁が拘束を受ける円板]
$E=70 \mathrm{GPa}$、$\alpha=\frac{1}{50000} \mathrm{K^{-1}}$、$A=2000 \mathrm{K\,m^{-2}}$、$a=\frac{1}{10} \mathrm m$ とする。以下では$T_{D}$ が縁と中心の温度差、$s_{0}$ が中心の等二軸応力、$s_{a}$ が外縁の周応力である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
T_{D}&=A a^{2}=20\quad[\mathrm K]\\
s_{0}&=250 A E a^{2} \alpha=7\quad[\mathrm{MPa}]\\
s_{a}&=- 500 A E a^{2} \alpha=-14\quad[\mathrm{MPa}]
\end{aligned}
$$

縁は半径方向に無荷重でも周方向に圧縮される。周方向の自由膨張が低温の内側によって制約されるからである。
:::

:::example[中心の方が高温の円板]
$E=70 \mathrm{GPa}$、$\alpha=\frac{1}{50000} \mathrm{K^{-1}}$、$A=-2000 \mathrm{K\,m^{-2}}$、$a=\frac{1}{10} \mathrm m$ とする。以下では$T_{D}$ が縁と中心の温度差、$s_{0}$ が中心の等二軸応力、$s_{a}$ が外縁の周応力である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
T_{D}&=A a^{2}=-20\quad[\mathrm K]\\
s_{0}&=250 A E a^{2} \alpha=-7\quad[\mathrm{MPa}]\\
s_{a}&=- 500 A E a^{2} \alpha=14\quad[\mathrm{MPa}]
\end{aligned}
$$

温度の二次係数の符号を変えると応力が全体で反転する。温度に定数を加えるだけなら応力は変わらない。
:::

:::example[半径を変えて熱応力を見積もる]
$E=70 \mathrm{GPa}$、$\alpha=\frac{1}{50000} \mathrm{K^{-1}}$、$A=2000 \mathrm{K\,m^{-2}}$ に固定し、円板半径 $a$ を $a=\frac{1}{10} \mathrm m$ の近くで変える。外縁の周応力 $s_{a}$ の局所的な変化を求める。固定した値を先に代入すると、残る一変数の関数は次のようになる。以下の微分では数値を先の単位系で表す。

$$
\begin{aligned}
s_{a}(a)&=- 1400 a^{2},\\
\frac{\mathrm ds_{a}}{\mathrm da}&=- 2800 a.
\end{aligned}
$$

基準値で $s_{a}=-14$、微分は $-280$ である。したがって一次の展開は、基準値からの相対変化を $a\mapsto a(1+\delta)$ と書くと

$$
\frac{\Delta s_{a}}{s_{a}}
=\left(\frac{a}{s_{a}}
\frac{\mathrm ds_{a}}{\mathrm da}\right)\delta+\order{\delta^2}
=2\,\delta+\order{\delta^2}
$$

となる。例えば入力を $1\%$ 増やすとき、一次式で得る出力の変化率は $2\%$ である。元の式に直接代入した変化率は $\frac{201}{100}\%$ なので、この比較で一次化の誤差を点検できる。

逆に、出力の測定値が基準値から $1\%$ 上がったとき、適用した模型が正しく、ほかの入力が固定されていれば、入力の変化率は一次で $\frac{1}{2}\%$ と推定される。温度の二次係数を固定する。縁と中心の温度差まで固定するなら応力は半径に依存しないため、熱境界条件に応じてどちらを固定するか決める。
:::

### 穴のある領域で線積分を点検する

単連結性がないと、適合の局所微分条件だけでは閉路の後に同じ変位へ戻ると保証できない。例えば $x^2+y^2>0$ の領域で反平面変位を局所的に $u_z=b\theta/(2\pi)$、$u_x=u_y=0$ とおく。$b$ は長さを単位とする定数、$\theta$ は偏角の局所枝である。枝の違いは定数なので、全領域で

$$
2\eps_{zx}=-\frac{b y}{2\pi(x^2+y^2)},\qquad
2\eps_{zy}=\frac{b x}{2\pi(x^2+y^2)}
$$

となる。他の独立成分はゼロである。局所的には微分できる変位を持つので、微分適合式はどこでも成立する。

同じひずみを持つ別の変位を十分小さな開領域で比べると、その差は微小剛体運動である。その勾配は局所的に定数反対称行列である。局所枝同士の変位差は定数なので、局所枝の勾配は大域的に定義される。したがって、もし大域的一価の変位が存在すれば、その勾配と局所枝の勾配の差は連結領域で一つの定数反対称行列となる。定数行列の閉路積分はゼロだが、局所枝の $u_z$ の勾配を一周積分した値は $\int_0^{2\pi}b\,\dd\theta/(2\pi)=b$ である。$b\ne0$ なら矛盾する。これは線状の :term[らせん転位]{en="screw dislocation" yomi="らせんてんい"} の芯を除いた変位とひずみの例である。

:::example[同じ閉路の変位積分と局所ひずみ]
$b=2.5\times10^{-10}\,\mathrm m$、$\mu=40\,\mathrm{GPa}$ とする。半径 $r=1.0\times10^{-8}\,\mathrm m$ の円でせん断ひずみの工学表示は $b/(2\pi r)\simeq3.98\times10^{-3}$。せん断応力は $\mu b/(2\pi r)\simeq159\,\mathrm{MPa}$。閉路距離 $2\pi r$ を掛けると、どの $r$ でも閉路積分は $b$ に戻る。芯に近づくと局所ひずみの方は増大し、線形弾性を使える下限の距離を考える必要がある。
:::

:::example[芯を除外したひずみエネルギー]
円環 $a<r<R$ で単位長さあたりのエネルギーを求める。せん断ひずみを $\gamma=b/(2\pi r)$ と書くと $W=\mu\gamma^2/2$ なので $2\pi\int_a^R W r\,\dd r=\mu b^2\log(R/a)/(4\pi)$。例えば $R/a=100$ から $10000$ に増やすと同じ $a$ でエネルギーは２倍になる。各対数区間が同じエネルギーを持ち、芯近くのひずみの最大値だけでは総エネルギーを見積もれない。
:::

周回で生じる不整合は、局所測定の分解能を高めても解消しない。領域を切り開いて変位を積分するなら、その切断の両側のずれを境界データとして保持することで、同じひずみ場を再現できる。

## 演習問題

::::exercise{level=basic}
$\lambda=20\,\mathrm{GPa}$、$\mu=30\,\mathrm{GPa}$ の $K,E,\nu$ を求めよ。

:::solution
$K=20+20=40\,\mathrm{GPa}$、$E=30(60+60)/(20+30)=72\,\mathrm{GPa}$、$\nu=20/[2(20+30)]=0.2$ である。
:::
::::

::::exercise{level=basic}
$K=50\,\mathrm{GPa}$ の固体に $p=5\,\mathrm{MPa}$ の一様圧を掛ける。体積変化率を求めよ。

:::solution
$\tr\sigma=-3p$ と $\tr\sigma=3K\tr\eps$ より $\tr\eps=-p/K=-10^{-4}$ である。行列式の一次項がこの跡なので体積が約 $0.01\%$ 減る。
:::
::::

::::exercise{level=basic}
$E=100\,\mathrm{GPa}$、$\nu=1/4$ の板が平面応力状態にあり、$\eps_{xx}=10^{-3}$、$\eps_{yy}=0$ とする。$\sigma_{xx},\sigma_{yy}$ を求めよ。

:::solution
$\sigma_{xx}=E\,10^{-3}/(1-1/16)=320/3\,\mathrm{MPa}$、$\sigma_{yy}=\nu\sigma_{xx}=80/3\,\mathrm{MPa}$ である。
:::
::::

::::exercise{level=basic}
$\alpha=10^{-5}\,\mathrm{K^{-1}}$、$K=50\,\mathrm{GPa}$ の固体を完全拘束して $20\,\mathrm K$ 加熱する。応力は何か。

:::solution
自由膨張ひずみは $2\times10^{-4}I$ で、実際のひずみはゼロである。したがって $\sigma=-3K\alpha\vartheta I=-30I\,\mathrm{MPa}$ となる。
:::
::::

::::exercise{level=advanced}
外圧ゼロの球殻で $b=a+h$、$h/a\ll1$ とする。内面周応力の支配項を求めよ。

:::solution
球殻解を代入して $s_t(a)=p_i(a^3+b^3/2)/(b^3-a^3)$ である。分子は $3p_i a^3/2+\order{p_i a^2h}$、分母は $3a^2h+\order{ah^2}$ なので $s_t(a)\sim p_i a/(2h)$ である。
:::
::::

::::exercise{level=advanced}
半径 $a$ の円柱の長さを $L$、ねじれ角を $\varphi=\vartheta L$ とする。ひずみエネルギーが $M\varphi/2$ であることを確かめよ。

:::solution
非零ひずみは $\eps_{xz}=-\vartheta y/2,\eps_{yz}=\vartheta x/2$ と対称成分である。$\tr\eps=0$、$W=\mu\vartheta^2r^2/2$。体積積分は $\mu\vartheta^2L\pi a^4/4$ で、$M=\mu\vartheta\pi a^4/2$ を代入すると $M\varphi/2$ になる。
:::
::::

::::exercise{level=advanced}
$\eps_{xx}=2xy,\eps_{yy}=0,\eps_{xy}=x^2/2$ を実現する平面変位を求めよ。

:::solution
$u_x=x^2 y$、$u_y=0$ とすると $\partial_xu_x=2xy$、$\partial_yu_y=0$、$(\partial_yu_x+\partial_xu_y)/2=x^2/2$ である。任意の微小剛体運動を加えても同じひずみを持つ。
:::
::::

::::exercise{level=advanced}
等方体の伸長波に現れる係数 $\lambda+2\mu$ が、$W$ の正定値性から正であることを示せ。一軸ひずみにおける仕事とも比較せよ。

:::solution
$\lambda+2\mu=K+4\mu/3>0$ である。一軸ひずみ $\eps=\xi\vb e_1\vb e_1^T$ のエネルギーは $W=(\lambda+2\mu)\xi^2/2$ となる。その微分は拘束された一軸ひずみで必要な $\sigma_{11}$ であり、側面自由の試験で定義した $E$ と一般に異なる。
:::
::::

## 確認テスト

::::quiz
平面応力でゼロとおく厚さ方向の成分を選べ。

- [x] $\sigma_{33}$
- [ ] $\eps_{33}$
- [ ] $\sigma_{11}+\sigma_{22}$

:::explanation
自由な薄板面では厚さ方向の面力が消える。一般に $\eps_{33}=-\nu(\sigma_{11}+\sigma_{22})/E$ は残り、面内に等二軸応力があっても平面応力の仮定と両立する。
:::
::::

::::quiz{answer="2" tol="1e-6"}
$\lambda=0$ の等方性体で $E/\mu$ はいくらか。

:::explanation
$\nu=0$ となり $E=2\mu(1+\nu)=2\mu$ である。
:::
::::

::::quiz{answer="4" tol="1e-6"}
$K=2\,\mathrm{GPa}$、$\mu=1.5\,\mathrm{GPa}$ の $\lambda+2\mu$ を $\mathrm{GPa}$ で答えよ。

:::explanation
$\lambda+2\mu=K+4\mu/3=2+2=4\,\mathrm{GPa}$ である。
:::
::::

::::quiz
変位に微小剛体並進を加えたときに変わる量を選べ。

- [x] 変位
- [ ] 微小ひずみ
- [ ] 線形弾性応力

:::explanation
並進ベクトルは一定なので勾配がゼロである。そのため変位だけが変わり、ひずみやそれから計算される応力には加算項がない。
:::
::::

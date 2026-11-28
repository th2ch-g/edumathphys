---
title: MHD 平衡と磁場閉じ込め
description: 静的な磁気力と圧力のつりあいから、円柱とトーラスの磁気面を求める。
notation:
  - symbol: '\psi'
    meaning: ポロイダル磁束を 2π で割った量
  - symbol: 'F'
    meaning: トロイダル磁場関数
  - symbol: 'q'
    meaning: 安全係数
  - symbol: '\beta'
    meaning: 熱圧力と磁気圧の比
---

## この章の目標

- 静的 MHD が許す圧力分布を調べる。
- ピンチと軸対称トーラスの平衡式を導く。
- 磁力線の巻き方と磁気シアを計算する。

## 力のつりあいと磁気面

静止場では :eqref[eq-plasma-ideal-mhd] が $\nabla p=\vb j\times\vb B$ に縮約する。内積をとると $\vb B\cdot\nabla p=\vb j\cdot\nabla p=0$ である。正則な等圧面には磁場と電流が接する。

:::definition[磁気面]{#def-plasma-magnetic-surface}
磁場が全点で接する滑らかな曲面を :term[磁気面]{en="magnetic surface" yomi="じきめん"} と呼ぶ。磁気面を $\psi=\text{一定}$ と表せば $\vb B\cdot\nabla\psi=0$ である。
:::

一本の磁力線が磁気面に稠密で、圧力が連続なら、その線で一定の圧力は面全体でも一定である。磁力線が閉じている場合は、磁場に沿う式だけで面上の一様性を結論できない。

:::example[直線磁場と全圧]
$\vb B=B_y(x)\hat{\vb y}$ では曲率がゼロなので $p+B_y^2/(2\mu_0)=C$ である。磁場が二倍になるとその磁気圧は四倍となり、その差を熱圧力が負担する。
:::

:::example[力のない磁場]
$\vb B=B_0(\sin kz,\cos kz,0)$ では $\nabla\times\vb B=k\vb B$ である。従って電流が流れていても磁気力は消え、圧力は一定にできる。
:::

### 無力場の電流

圧力勾配がゼロで $B\ne0$ の領域では、電流は磁場と同じ直線上を向き、 $\nabla\times\vb B=\alpha(\vb x)\vb B$ と書ける。このような場を :term[無力場]{en="force-free magnetic field" yomi="むりょくば"} と呼ぶ。両辺の発散を取ると $\vb B\cdot\nabla\alpha=0$ である。従って $\alpha$ は磁力線に沿って一定である。$\alpha$ が領域全体で一定なら、もう一度回転を取って $\nabla^2\vb B+\alpha^2\vb B=0$ となる。真空場では $\alpha=0$ であり、各成分は Laplace 方程式を満たす。

例えば $\vb B=B_0(\sin kz,\cos kz,0)$ なら $\alpha=k$ で、$k$ の符号を変えると電流が磁場と平行か反平行かが変わる。どちらでも $\vb j\times\vb B$ はゼロである。有限領域にこの解を使うときは、境界に要求した法線磁場と合うかを確認する。

### 圧力と電流の成分

$\nabla p$ が決まっていれば、横電流は $\vb j_\perp=\vb B\times\nabla p/B^2$ である。これは $\vb B\times(\vb j\times\vb B)=B^2\vb j-\vb B(\vb B\cdot\vb j)$ から従う。$\vb j=\vb j_\perp+\sigma\vb B$ を $\nabla\cdot\vb j=0$ に入れると $\vb B\cdot\nabla\sigma=-\nabla\cdot\vb j_\perp$ である。閉じた磁力線に沿って $\sigma$ を一価に保つには、線上で $(\nabla\cdot\vb j_\perp)/B$ の積分がゼロでなければならない。この制約は局所の力の釣り合いからだけでは見えない。

:::example[直線場での電流]
$\vb B=B_y(x)\hat{\vb y}$、$p=p(x)$ では $\vb j_\perp=-p'\hat{\vb z}/B_y$ で発散はゼロである。Ampère 則の $j_z=B_y'/\mu_0$ と一致するのは $p+B_y^2/(2\mu_0)$ が一定の場合である。
:::

### 有限領域の磁気ビリアル関係

重力や外部の体積力がないとし、静的な熱・磁気応力を $\mathsf T=\vb B\vb B/\mu_0-(p+B^2/(2\mu_0))I$ と置く。平衡は $\partial_jT_{ij}=0$ である。座標 $x_i$ を掛けると $\partial_j(x_iT_{ij})=T_{ii}=-3p-B^2/(2\mu_0)$ なので、体積で積分して
$$
\int_V \left(3p+\frac{B^2}{2\mu_0}\right)\dd V
=-\int_{\partial V}\vb x\cdot\mathsf T\hat{\vb n}\,\dd S
$$
となる。内部に正の圧力または磁場があれば左辺は正であるため、境界で外側と交換する応力が必要である。全空間で表面積分が無限遠で消えるほど急に場が減衰し、$p\ge0$ で他に力がなければ、平衡は自明な $p=B=0$ に限られる。トーラスのコイル、周辺の真空磁場、導体壁を扱うときに、どこまでを積分領域に含めたかで表面の意味が変わる。

:::example[一定の外部圧力]
境界応力が $-p_w\hat{\vb n}$ なら右辺は $p_w\int_{\partial V}\vb x\cdot\hat{\vb n}\dd S=3p_w V$ である。従って磁気エネルギーを含む体積平均は $3p_w$ に等しい。内部熱圧だけを $p_w$ と比較しても、磁場の応力が寄与する。
:::

### 静止する直線場における全圧力の復元

$\vb B=B_z(x)\hat{\vb z}$ では $J_y=-B_z'/\mu_0$、$p'=-B_zB_z'/\mu_0$ となる。積分定数を境界の全圧から定めれば $p=p_a+(B_a^2-B_z^2)/(2\mu_0)$ である。これは磁場分布を自由に選んだときに必要な熱圧力を与えるが、負圧を予言する入力は物理的平衡にならない。逆に与えた熱圧力から磁場を解く場合、二乗だけが決まるので、磁場の向きや電流を求めるには境界条件も指定する。磁場がゼロになる点ではその符号を任意に切り替えると微分可能性が失われ得る。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $b$ | 磁場 $B_z/B_a$ | 1 |
| $g$ | 磁場勾配 $L B_z'/B_a$ | 1 |
| $\beta_a$ | 境界熱圧を $B_a^2/(2\mu_0)$ で割った値 | 1 |

規格化磁気圧を $p_B/p_{B,a}$ とする。符号付き磁場を二乗する。

$$
p_B/p_{B,a}=b^{2}
$$

局所熱圧を $p/p_{B,a}$ とする。全圧を固定する。

$$
p/p_{B,a}=\beta_{a} - p_B/p_{B,a} + 1
$$

電流の規格化値を $\mu_0 L J_y/B_a$ とする。curl の向きで負符号が付く。

$$
\mu_0 L J_y/B_a=- g
$$

必要な熱圧勾配を $Lp'/p_{B,a}$ とする。磁気圧勾配と逆である。

$$
Lp'/p_{B,a}=- 2 b g
$$

全圧力を $(p+p_B)/p_{B,a}$ とする。$x$ に依存しない。

$$
(p+p_B)/p_{B,a}=p/p_{B,a} + p_B/p_{B,a}
$$

局所の熱圧と磁気圧の比を $\beta$ とする。$B_z\ne0$ に限る。

$$
\beta=\frac{p/p_{B,a}}{p_B/p_{B,a}}
$$

:::example[高圧の内部]
外側へ磁場が増えるときの熱圧を求める。

入力値を次のように置く。

$$
\begin{aligned}
b&=\frac{1}{2} \\g&=\frac{1}{2} \\\beta_a&=\frac{1}{10}
\end{aligned}
$$

これらを順に代入する。

規格化磁気圧は

$$
p_B/p_{B,a}=\left(\frac{1}{2}\right)^{2}=\frac{1}{4}
$$

である。

局所熱圧は

$$
p/p_{B,a}=\left(\frac{1}{10}\right) - \left(\frac{1}{4}\right) + 1=\frac{17}{20}
$$

である。

電流の規格化値は

$$
\mu_0 L J_y/B_a=- \left(\frac{1}{2}\right)=- \frac{1}{2}
$$

である。

必要な熱圧勾配は

$$
Lp'/p_{B,a}=- 2 \left(\frac{1}{2}\right)^{2}=- \frac{1}{2}
$$

である。

全圧力は

$$
(p+p_B)/p_{B,a}=\left(\frac{17}{20}\right) + \left(\frac{1}{4}\right)=\frac{11}{10}
$$

である。

局所の熱圧と磁気圧の比は

$$
\beta=\frac{\left(\frac{17}{20}\right)}{\left(\frac{1}{4}\right)}=\frac{17}{5}
$$

である。

内部の熱圧は境界より大きい。
:::

:::example[磁場向きの反転]
磁場を全体に反転させて力を比較する。

入力値を次のように置く。

$$
\begin{aligned}
b&=- \frac{1}{2} \\g&=- \frac{1}{2} \\\beta_a&=\frac{1}{10}
\end{aligned}
$$

これらを順に代入する。

規格化磁気圧は

$$
p_B/p_{B,a}=\left(- \frac{1}{2}\right)^{2}=\frac{1}{4}
$$

である。

局所熱圧は

$$
p/p_{B,a}=\left(\frac{1}{10}\right) - \left(\frac{1}{4}\right) + 1=\frac{17}{20}
$$

である。

電流の規格化値は

$$
\mu_0 L J_y/B_a=- \left(- \frac{1}{2}\right)=\frac{1}{2}
$$

である。

必要な熱圧勾配は

$$
Lp'/p_{B,a}=- 2 \left(- \frac{1}{2}\right)^{2}=- \frac{1}{2}
$$

である。

全圧力は

$$
(p+p_B)/p_{B,a}=\left(\frac{17}{20}\right) + \left(\frac{1}{4}\right)=\frac{11}{10}
$$

である。

局所の熱圧と磁気圧の比は

$$
\beta=\frac{\left(\frac{17}{20}\right)}{\left(\frac{1}{4}\right)}=\frac{17}{5}
$$

である。

電流も反転するので力は同じである。
:::

:::example[磁場が強すぎる分布]
境界圧の値で平衡が可能か調べる。

入力値を次のように置く。

$$
\begin{aligned}
b&=2 \\g&=1 \\\beta_a&=\frac{1}{10}
\end{aligned}
$$

これらを順に代入する。

規格化磁気圧は

$$
p_B/p_{B,a}=\left(2\right)^{2}=4
$$

である。

局所熱圧は

$$
p/p_{B,a}=- \left(4\right) + \left(\frac{1}{10}\right) + 1=- \frac{29}{10}
$$

である。

電流の規格化値は

$$
\mu_0 L J_y/B_a=- \left(1\right)=-1
$$

である。

必要な熱圧勾配は

$$
Lp'/p_{B,a}=- 2 \left(1\right) \left(2\right)=-4
$$

である。

全圧力は

$$
(p+p_B)/p_{B,a}=\left(- \frac{29}{10}\right) + \left(4\right)=\frac{11}{10}
$$

である。

局所の熱圧と磁気圧の比は

$$
\beta=\frac{\left(- \frac{29}{10}\right)}{\left(4\right)}=- \frac{29}{40}
$$

である。

熱圧が負になるため、指定した磁場と境界条件は両立しない。
:::

:::example[勾配のない区間]
熱圧が有限でも電流が必要か調べる。

入力値を次のように置く。

$$
\begin{aligned}
b&=\frac{3}{4} \\g&=0 \\\beta_a&=\frac{1}{4}
\end{aligned}
$$

これらを順に代入する。

規格化磁気圧は

$$
p_B/p_{B,a}=\left(\frac{3}{4}\right)^{2}=\frac{9}{16}
$$

である。

局所熱圧は

$$
p/p_{B,a}=\left(\frac{1}{4}\right) - \left(\frac{9}{16}\right) + 1=\frac{11}{16}
$$

である。

電流の規格化値は

$$
\mu_0 L J_y/B_a=- \left(0\right)=0
$$

である。

必要な熱圧勾配は

$$
Lp'/p_{B,a}=- 2 \left(0\right) \left(\frac{3}{4}\right)=0
$$

である。

全圧力は

$$
(p+p_B)/p_{B,a}=\left(\frac{11}{16}\right) + \left(\frac{9}{16}\right)=\frac{5}{4}
$$

である。

局所の熱圧と磁気圧の比は

$$
\beta=\frac{\left(\frac{11}{16}\right)}{\left(\frac{9}{16}\right)}=\frac{11}{9}
$$

である。

一様な区間では全圧のつりあいに体積電流を要しない。
:::


## ピンチ平衡（θ ピンチ・Z ピンチと Bennett 関係）

円筒の軸対称平衡を $\vb B=B_\theta(r)\hat{\vb\theta}+B_z(r)\hat{\vb z}$ とする。張力は方位角成分から生じる。式 :eqref[eq-plasma-magnetic-stress] の動径成分は
$$
\dv{}{r}\left(p+\frac{B_\theta^2+B_z^2}{2\mu_0}\right)=-\frac{B_\theta^2}{\mu_0r}.
\label{eq-plasma-pinch-balance}
$$
$\theta$ ピンチでは $B_\theta=0$ とし、外部軸磁場と方位角電流が全圧を支える。Z ピンチでは軸電流が $B_\theta$ を作り、その張力が圧縮に寄与する。

### Bennett 関係

純 Z ピンチで半径 $a$ の外部熱圧を $p(a)$ とする。式 :eqref[eq-plasma-pinch-balance] に $r^2$ を掛けると $r^2 p'=-(r^2 B_\theta^2)'/(2\mu_0)$ となる。軸で場は有界とし、積分して部分積分すると $2\int_0^a r[p-p(a)]\dd r=a^2 B_\theta(a)^2/(2\mu_0)$ である。Ampère 則で $B_\theta(a)=\mu_0I/(2\pi a)$ だから
$$
\int_{r<a}[p-p(a)]\,\dd A=\frac{\mu_0 I^2}{8\pi}.
\label{eq-plasma-bennett}
$$
周辺圧力がゼロで一価で準中性かつ種ごとの温度が一様なら左辺は一単位長さ当たりのイオン数 $N_\ell$ と $T_i+T_e$ の積である。これを :term[Bennett 関係]{en="Bennett relation" yomi="べねっとかんけい"} という。

:::example[一様電流]
$B_\theta=br$ なら $p(r)-p(a)=b^2(a^2-r^2)/\mu_0$ である。圧力差の断面平均は軸上の半分である。
:::

:::example[一様軸磁場の効果]
一様な軸方向の磁場を上の Z ピンチに足しても、動径力は変わらない。ただし、軸方向に変化する横変位は軸磁場を曲げるので、安定性を決めるエネルギーが変わる。
:::

### 外部磁気圧の境界条件

プラズマから真空へ接する境界で $B_n=0$、表面が静止しているとする。全応力の法線成分の連続条件は $p_{\rm in}+B_{\rm in}^2/(2\mu_0)=B_{\rm out}^2/(2\mu_0)$ である。接線磁場に跳びがあれば、薄い表面層に表面電流 $\vb K=\hat{\vb n}\times(\vb B_{\rm out}-\vb B_{\rm in})/\mu_0$ がある。表面電流と全圧の条件を合わせて境界を解く。

:::example[磁束排除の極限]
$\theta$ ピンチで内部の軸磁場を十分小さいとすれば、境界圧力は外部磁気圧 $B_{\rm out}^2/(2\mu_0)$ に近づく。外部磁場を二倍にすればこの支圧能力は四倍となる。境界層の厚さ自体は静的な跳びの条件からは分からない。
:::

### べき型 Z ピンチの圧力と Bennett 積分

半径 $a$ の円柱で $B_\theta=B_a(r/a)^m$、$m>0$ とし、$p(a)=p_a$ を与える。動径平衡により $p'=-(m+1)B_a^2(r/a)^{2m-1}/(\mu_0a)$、積分して $p-p_a=(m+1)B_a^2[1-(r/a)^{2m}]/(2m\mu_0)$ である。断面平均 $\langle(r/a)^{2m}\rangle=1/(m+1)$ により $\langle p\rangle-p_a=B_a^2/(2\mu_0)$ となる。電流は $I=2\pi a B_a/\mu_0$ だから、断面積を掛けると $\int(p-p_a)\dd A=\mu_0 I^2/(8\pi)$ になる。軸で有限な電流密度を求めるなら $m\ge1$ が必要である。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $m$ | 磁場の動径指数 | 1 |
| $x$ | 動径 $r/a$ | 1 |
| $p_a/p_{B,a}$ | 辺での熱圧と辺の磁気圧の比 | 1 |

局所磁場を $B_\theta/B_a$ とする。各例で同じ辺の磁場を基準にする。

$$
B_\theta/B_a=x^{m}
$$

軸電流密度を $\mu_0 aJ_z/B_a$ とする。円筒の Ampère の法則を使う。

$$
\mu_0 aJ_z/B_a=x^{m - 1} \left(m + 1\right)
$$

局所熱圧を $p/p_{B,a}$ とする。張力と磁気圧をともに含めて積分する。

$$
p/p_{B,a}=p_a/p_{B,a} + \frac{\left(1 - x^{2 m}\right) \left(m + 1\right)}{m}
$$

中心圧力を $p(0)/p_{B,a}$ とする。$m>0$ では軸で $B_\theta$ が消える。

$$
p(0)/p_{B,a}=p_a/p_{B,a} + \frac{m + 1}{m}
$$

断面平均熱圧を $\langle p\rangle/p_{B,a}$ とする。面積の重み $2x$ を使う。

$$
\langle p\rangle/p_{B,a}=p_a/p_{B,a} + 1
$$

内側に含まれる電流の割合を $I(r)/I(a)$ とする。円周積分では半径の因子も入る。

$$
I(r)/I(a)=(B_\theta/B_a) x
$$

:::example[一様電流の断面]
軸から半径半分での電流と圧力を求める。

入力値を次のように置く。

$$
\begin{aligned}
m&=1 \\x&=\frac{1}{2} \\p_a/p_{B,a}&=0
\end{aligned}
$$

これらを順に代入する。

局所磁場は

$$
B_\theta/B_a=\left(\frac{1}{2}\right)^{\left(1\right)}=\frac{1}{2}
$$

である。

軸電流密度は

$$
\mu_0 aJ_z/B_a=\left(\frac{1}{2}\right)^{\left(1\right) - 1} \left(\left(1\right) + 1\right)=2
$$

である。

局所熱圧は

$$
p/p_{B,a}=\left(0\right) + \frac{\left(1 - \left(\frac{1}{2}\right)^{2 \left(1\right)}\right) \left(\left(1\right) + 1\right)}{\left(1\right)}=\frac{3}{2}
$$

である。

中心圧力は

$$
p(0)/p_{B,a}=\left(0\right) + \frac{\left(1\right) + 1}{\left(1\right)}=2
$$

である。

断面平均熱圧は

$$
\langle p\rangle/p_{B,a}=\left(0\right) + 1=1
$$

である。

内側に含まれる電流の割合は

$$
I(r)/I(a)=\left(\frac{1}{2}\right)^{2}=\frac{1}{4}
$$

である。

内側の電流は面積分率と一致する。
:::

:::example[外側に偏った電流]
同じ全電流で圧力分布を比べる。

入力値を次のように置く。

$$
\begin{aligned}
m&=3 \\x&=\frac{1}{2} \\p_a/p_{B,a}&=0
\end{aligned}
$$

これらを順に代入する。

局所磁場は

$$
B_\theta/B_a=\left(\frac{1}{2}\right)^{\left(3\right)}=\frac{1}{8}
$$

である。

軸電流密度は

$$
\mu_0 aJ_z/B_a=\left(\frac{1}{2}\right)^{\left(3\right) - 1} \left(\left(3\right) + 1\right)=1
$$

である。

局所熱圧は

$$
p/p_{B,a}=\left(0\right) + \frac{\left(1 - \left(\frac{1}{2}\right)^{2 \left(3\right)}\right) \left(\left(3\right) + 1\right)}{\left(3\right)}=\frac{21}{16}
$$

である。

中心圧力は

$$
p(0)/p_{B,a}=\left(0\right) + \frac{\left(3\right) + 1}{\left(3\right)}=\frac{4}{3}
$$

である。

断面平均熱圧は

$$
\langle p\rangle/p_{B,a}=\left(0\right) + 1=1
$$

である。

内側に含まれる電流の割合は

$$
I(r)/I(a)=\left(\frac{1}{2}\right) \left(\frac{1}{8}\right)=\frac{1}{16}
$$

である。

中心圧は低くなるが断面平均は同じである。
:::

:::example[有限の外側圧力]
外圧を一定だけ上げて電流が変わるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
m&=1 \\x&=\frac{1}{2} \\p_a/p_{B,a}&=1
\end{aligned}
$$

これらを順に代入する。

局所磁場は

$$
B_\theta/B_a=\left(\frac{1}{2}\right)^{\left(1\right)}=\frac{1}{2}
$$

である。

軸電流密度は

$$
\mu_0 aJ_z/B_a=\left(\frac{1}{2}\right)^{\left(1\right) - 1} \left(\left(1\right) + 1\right)=2
$$

である。

局所熱圧は

$$
p/p_{B,a}=\left(1\right) + \frac{\left(1 - \left(\frac{1}{2}\right)^{2 \left(1\right)}\right) \left(\left(1\right) + 1\right)}{\left(1\right)}=\frac{5}{2}
$$

である。

中心圧力は

$$
p(0)/p_{B,a}=\left(1\right) + \frac{\left(1\right) + 1}{\left(1\right)}=3
$$

である。

断面平均熱圧は

$$
\langle p\rangle/p_{B,a}=\left(1\right) + 1=2
$$

である。

内側に含まれる電流の割合は

$$
I(r)/I(a)=\left(\frac{1}{2}\right)^{2}=\frac{1}{4}
$$

である。

圧力全体が押し上げられるが、勾配と磁場は変わらない。
:::

:::example[円柱の辺]
境界と断面平均の違いを確認する。

入力値を次のように置く。

$$
\begin{aligned}
m&=3 \\x&=1 \\p_a/p_{B,a}&=\frac{1}{2}
\end{aligned}
$$

これらを順に代入する。

局所磁場は

$$
B_\theta/B_a=\left(1\right)^{\left(3\right)}=1
$$

である。

軸電流密度は

$$
\mu_0 aJ_z/B_a=\left(1\right)^{\left(3\right) - 1} \left(\left(3\right) + 1\right)=4
$$

である。

局所熱圧は

$$
p/p_{B,a}=\left(\frac{1}{2}\right) + \frac{\left(1 - \left(1\right)^{2 \left(3\right)}\right) \left(\left(3\right) + 1\right)}{\left(3\right)}=\frac{1}{2}
$$

である。

中心圧力は

$$
p(0)/p_{B,a}=\left(\frac{1}{2}\right) + \frac{\left(3\right) + 1}{\left(3\right)}=\frac{11}{6}
$$

である。

断面平均熱圧は

$$
\langle p\rangle/p_{B,a}=\left(\frac{1}{2}\right) + 1=\frac{3}{2}
$$

である。

内側に含まれる電流の割合は

$$
I(r)/I(a)=\left(1\right)^{2}=1
$$

である。

辺では設定した外圧と一致し、平均にはピンチによる増分が残る。
:::


## Grad–Shafranov 方程式

トーラスの軸対称平衡を円筒座標 $(R,\varphi,Z)$ で記述する。$\partial_\varphi=0$、$\nabla\varphi=\hat{\vb\varphi}/R$ とし、規約を
$$
\vb B=\nabla\psi\times\nabla\varphi+F\nabla\varphi,\qquad
\Delta^*\psi=R\partial_R(R^{-1}\partial_R\psi)+\partial_Z^2\psi
$$
と取る。$\nabla\cdot\vb B=0$ から$R>0$ の単連結なポロイダル領域ではこの磁束表示が可能である。穴がある領域では、閉曲線に沿う磁束関数の増分がゼロであることも要る。実際 $\partial_R(R B_R)+\partial_Z(R B_Z)=0$ なので、一価の流れ関数で $R B_R=-\partial_Z\psi$、$R B_Z=\partial_R\psi$ と表せる。

方位角方向の力がゼロなので $\nabla F$ と $\nabla\psi$ が平行となり、連結した正則面で $F=F(\psi)$ である。$\vb B_p\cdot\nabla p=0$ から、連結した正則面では圧力も $p=p(\psi)$ となる。curl を計算すると $\mu_0\vb j=\nabla F\times\nabla\varphi-\Delta^*\psi\,\nabla\varphi$ である。外積を取れば $\vb j\times\vb B=-(FF'+\Delta^*\psi)\nabla\psi/(\mu_0R^2)$ が出る。従って
$$
\Delta^*\psi=-\mu_0R^2 p'(\psi)-F F'(\psi).
\label{eq-plasma-grad-shafranov}
$$
これを :term[Grad–Shafranov 方程式]{en="Grad–Shafranov equation" yomi="ぐらーどしゃふらのふほうていしき"} と呼ぶ。磁気軸では $\nabla\psi=0$ なので導出でそこを割らず、正則な周囲の極限で式を延長する。

:::example[多項式解]
$p'=p_1$、$FF'=g_1$ が定数なら $\psi=-\mu_0p_1R^4/8-g_1 Z^2/2$ は特解である。$\Delta^*R^4=8R^2$、$\Delta^*Z^2=2$ だからである。実際の面は境界を満たす同次解を足して求める。
:::

:::example[真空垂直磁場]
$\psi=B_0R^2/2$、$F=0$ なら $\Delta^*\psi=0$ であり $\vb B=B_0\hat{\vb z}$ である。
:::

### 円形断面に近い多項式

磁気軸を $R=R_0,Z=0$ 付近に置き、$\psi=\alpha(R^2-R_0^2)^2+\beta_0 Z^2$ とする。$\Delta^*\psi=8\alpha R^2+2\beta_0$ なので、$p'=-8\alpha/\mu_0$、$FF'=-2\beta_0$ のとき厳密解である。$R-R_0=x$ が小さければ $\psi\simeq4\alpha R_0^2x^2+\beta_0 Z^2$ となり、$\beta_0=4\alpha R_0^2$ なら軸近くで断面は円となる。$\alpha,\beta_0>0$ なら外向きに $\psi$ が増え、圧力は下がる。$F^2$ が全領域で正となるように積分定数と境界を選ぶ。

:::example[縦に長い断面]
$\beta_0=\alpha R_0^2$ なら、小さい $\psi$ の等値面の $Z$ 方向半径は $R$ 方向の約二倍である。これは多項式の軸近傍の形状であり、遠くの面では $(R^2-R_0^2)^2$ の非対称性が効く。
:::

### Grad–Shafranov 演算子の多項式検査

円柱座標で軸対称とし $\vb B_p=(\grad\psi\times\hat{\vb\varphi})/R$ とする。$F=RB_\varphi$、$\Delta_*\psi=\psi_{RR}-R^{-1}\psi_R+\psi_{ZZ}$ とすると $\Delta_*\psi=-\mu_0R^2p'(\psi)-FF'(\psi)$ である。$\psi=aR^4+bR^2Z^2+cR^2+d Z^2$ に作用させれば $\Delta_*\psi=(8a+2b)R^2+2d$ を得る。従って $p'=-(8a+2b)/\mu_0$、$F^2=F_0^2-4d\psi$ なら強制項と一致する。$cR^2$ は演算子の核に属し、垂直磁場を調整する。ここでは長さと磁束を単位で割り、数値は無次元で計算する。$F^2\ge0$ と、選んだ等磁束面が装置の境界に合うかは別に確認する。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $a$ | $R^4$ の係数 | 1 |
| $b$ | $R^2Z^2$ の係数 | 1 |
| $c$ | $R^2$ の係数 | 1 |
| $d$ | $Z^2$ の係数 | 1 |
| $R$ | 評価点の大半径 | 1 |
| $Z$ | 評価点の高さ | 1 |

磁束関数を $\psi$ とする。この関数の等値面で磁場が接する。

$$
\psi=R^{4} a + R^{2} Z^{2} b + R^{2} c + Z^{2} d
$$

演算子の値を $\Delta_*\psi$ とする。$cR^2$ の寄与は消える。

$$
\Delta_*\psi=R^{2} \left(8 a + 2 b\right) + 2 d
$$

磁場の半径成分を $B_R$ とする。規約により $-\partial_Z\psi/R$ である。

$$
B_R=\frac{- 2 R^{2} Z b - 2 Z d}{R}
$$

磁場の垂直成分を $B_Z$ とする。$\partial_R\psi/R$ を使う。

$$
B_Z=4 R^{2} a + 2 Z^{2} b + 2 c
$$

圧力関数の傾きを $\mu_0 p'$ とする。$R^2$ の係数を照合する。

$$
\mu_0 p'=- 8 a - 2 b
$$

トロイダル電流密度を $\mu_0 J_\varphi$ とする。$\vb B_p$ の curl から得る。

$$
\mu_0 J_\varphi=- \frac{\Delta^*\psi}{R}
$$

:::example[垂直一様場だけを足す]
磁束を変えても体積電流が生じないか調べる。

入力値を次のように置く。

$$
\begin{aligned}
a&=0 \\b&=0 \\c&=1 \\d&=0 \\R&=2 \\Z&=1
\end{aligned}
$$

これらを順に代入する。

磁束関数は

$$
\psi=\left(0\right) \left(1\right)^{2} \left(2\right)^{2} + \left(0\right) \left(1\right)^{2} + \left(0\right) \left(2\right)^{4} + \left(1\right) \left(2\right)^{2}=4
$$

である。

演算子の値は

$$
\Delta_*\psi=10 \left(0\right) \left(2\right)^{2} + 2 \left(0\right)=0
$$

である。

磁場の半径成分は

$$
B_R=\frac{- 2 \left(0\right) \left(1\right) \left(2\right)^{2} - 2 \left(0\right) \left(1\right)}{\left(2\right)}=0
$$

である。

磁場の垂直成分は

$$
B_Z=2 \left(0\right) \left(1\right)^{2} + 4 \left(0\right) \left(2\right)^{2} + 2 \left(1\right)=2
$$

である。

圧力関数の傾きは

$$
\mu_0 p'=- 10 \left(0\right)=0
$$

である。

トロイダル電流密度は

$$
\mu_0 J_\varphi=- \frac{\left(0\right)}{\left(2\right)}=0
$$

である。

$R^2$ の項は $B_Z$ を一様に変える。
:::

:::example[圧力勾配をもつ多項式]
圧力勾配とポロイダル磁場を計算する。

入力値を次のように置く。

$$
\begin{aligned}
a&=1 \\b&=1 \\c&=-2 \\d&=0 \\R&=1 \\Z&=1
\end{aligned}
$$

これらを順に代入する。

磁束関数は

$$
\psi=\left(-2\right) \left(1\right)^{2} + \left(0\right) \left(1\right)^{2} + 2 \left(1\right)^{5}=0
$$

である。

演算子の値は

$$
\Delta_*\psi=2 \left(0\right) + 10 \left(1\right)^{3}=10
$$

である。

磁場の半径成分は

$$
B_R=\frac{- 2 \left(0\right) \left(1\right) - 2 \left(1\right)^{4}}{\left(1\right)}=-2
$$

である。

磁場の垂直成分は

$$
B_Z=2 \left(-2\right) + 6 \left(1\right)^{3}=2
$$

である。

圧力関数の傾きは

$$
\mu_0 p'=- 10 \left(1\right)=-10
$$

である。

トロイダル電流密度は

$$
\mu_0 J_\varphi=- \frac{\left(10\right)}{\left(1\right)}=-10
$$

である。

$R^2$ の源を圧力関数の傾きが支える。
:::

:::example[トロイダル磁場関数の源]
圧力が一定でもトロイダル電流があるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
a&=0 \\b&=0 \\c&=0 \\d&=1 \\R&=2 \\Z&=1
\end{aligned}
$$

これらを順に代入する。

磁束関数は

$$
\psi=\left(0\right) \left(1\right)^{2} \left(2\right)^{2} + \left(0\right) \left(2\right)^{4} + \left(0\right) \left(2\right)^{2} + \left(1\right)^{3}=1
$$

である。

演算子の値は

$$
\Delta_*\psi=10 \left(0\right) \left(2\right)^{2} + 2 \left(1\right)=2
$$

である。

磁場の半径成分は

$$
B_R=\frac{- 2 \left(0\right) \left(1\right) \left(2\right)^{2} - 2 \left(1\right)^{2}}{\left(2\right)}=-1
$$

である。

磁場の垂直成分は

$$
B_Z=2 \left(0\right) \left(1\right)^{2} + 4 \left(0\right) \left(2\right)^{2} + 2 \left(0\right)=0
$$

である。

圧力関数の傾きは

$$
\mu_0 p'=- 10 \left(0\right)=0
$$

である。

トロイダル電流密度は

$$
\mu_0 J_\varphi=-1
$$

である。

一定の $\Delta_*\psi$ は $FF'$ によって支えられる。
:::

:::example[軸点近くの設計]
ポロイダル磁場が消える点と電流を比べる。

入力値を次のように置く。

$$
\begin{aligned}
a&=1 \\b&=1 \\c&=-2 \\d&=0 \\R&=1 \\Z&=0
\end{aligned}
$$

これらを順に代入する。

磁束関数は

$$
\psi=\left(-2\right) \left(1\right)^{2} + \left(0\right)^{3} + \left(0\right)^{2} \left(1\right)^{3} + \left(1\right)^{5}=-1
$$

である。

演算子の値は

$$
\Delta_*\psi=2 \left(0\right) + 10 \left(1\right)^{3}=10
$$

である。

磁場の半径成分は

$$
B_R=\frac{- 2 \left(0\right)^{2} - 2 \left(0\right) \left(1\right)^{3}}{\left(1\right)}=0
$$

である。

磁場の垂直成分は

$$
B_Z=2 \left(-2\right) + 2 \left(0\right)^{2} \left(1\right) + 4 \left(1\right)^{3}=0
$$

である。

圧力関数の傾きは

$$
\mu_0 p'=- 10 \left(1\right)=-10
$$

である。

トロイダル電流密度は

$$
\mu_0 J_\varphi=- \frac{\left(10\right)}{\left(1\right)}=-10
$$

である。

磁場の微分が有限なので軸点でも電流はゼロと限らない。
:::


### 磁気軸を持つ Solov'ev 型多項式

真空領域の方程式と境界条件を指定する前に、磁気軸の近くで閉じた磁気面を作れるかを調べる。$R_0>0$、$A,C>0$ とし
$$
\psi(R,Z)=A(R^2-R_0^2)^2+C R^2Z^2
$$
を選ぶ。$\psi_R=4AR(R^2-R_0^2)+2CRZ^2$、$\psi_Z=2CR^2Z$ で、$(R_0,0)$ が停留点である。Hessian の同点での成分は $\psi_{RR}=8AR_0^2$、$\psi_{ZZ}=2CR_0^2$、$\psi_{RZ}=0$ なので非退化な極小となる。軸に十分近い等値線は二次の楕円 $4AR_0^2(R-R_0)^2+CR_0^2 Z^2=\psi$ で近似できる。$\Delta_*\psi=(8A+2C)R^2$ だから $F$ 一定、$p'=-(8A+2C)/\mu_0$ で一致する。圧力を釣り合わせるには辺の磁束値 $\psi_a$ で熱圧力をゼロと置くと $p=(8A+2C)(\psi_a-\psi)/\mu_0$ と選べる。従って濃い中心を持つ圧力平衡を構成できる。

同じ等値線で上下幅を厳密に測る場合には $Z^2=[\psi_a-A(R^2-R_0^2)^2]/(CR^2)$ を $R$ で最大化する。$R=R_0$ で測る高さは最大値と一般には一致しない。以下では軸を通る水平・垂直線との交差を比較する。この模型は内部の境界値問題であり、辺での真空磁場の接続や外部コイル電流を定めるものではない。

:::remark[軸近傍の面と境界]
$R_0$ は曲面の磁気軸である。等値面の内側と外側の端から計算した幾何学的な中心は $R_0$ と一致しない。実際、平方根は正の引数で狭義凹関数だから、$D>0$ なら $(\sqrt{R_0^2+D}+\sqrt{R_0^2-D})/2<R_0$ である。これはこの固定係数の解における幾何学的な差である。外部の真空領域にある磁束まで定めるときは、境界で法線磁場を連続にし、接線磁場の飛びを許すなら、それに対応する表面電流を指定する。
:::

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $R_0$ | 軸の大半径を長さ単位で割った量 | 1 |
| $A$ | $R$ 方向の磁束係数 | 1 |
| $C$ | $Z$ 方向の磁束係数 | 1 |
| $\psi_a$ | 評価する磁気面の磁束値 | 1 |

磁気面の水平端での $\abs{R^2-R_0^2}$ を $D$ とする。$Z=0$ を代入する。

$$
D=\sqrt{\frac{\psi_{a}}{A}}
$$

外側の水平端を $R_+$ とする。正の半径の解を選ぶ。

$$
R_+=\sqrt{D + R_{0}^{2}}
$$

内側の水平端を $R_-$ とする。$D<R_0^2$ の磁気面を取る。

$$
R_-=\sqrt{- D + R_{0}^{2}}
$$

水平半幅を $h_R$ とする。水平端を軸からの距離で単純に置き換えない。

$$
h_R=\frac{R_{+}}{2} - \frac{R_{-}}{2}
$$

軸の大半径上での垂直半幅を $h_Z$ とする。$R=R_0$ を代入する。

$$
h_Z=\frac{\sqrt{\frac{\psi_{a}}{C}}}{R_{0}}
$$

軸近傍楕円の伸長度を $\kappa_{\rm lin}$ とする。二次係数の比の平方根である。

$$
\kappa_{\rm lin}=2 \sqrt{\frac{A}{C}}
$$

:::example[軸に近い円形]
二次で伸長度一の入力の半幅を調べる。

入力値を次のように置く。

$$
\begin{aligned}
R_0&=3 \\A&=1 \\C&=4 \\\psi_a&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

磁気面の水平端での$\abs{R^2-R_0^2}=D$ は

$$
D=\sqrt{\frac{\left(\frac{1}{100}\right)}{\left(1\right)}}=\frac{1}{10}
$$

である。

外側の水平端は

$$
R_+=\sqrt{\left(3\right)^{2} + \left(\frac{1}{10}\right)}\simeq 3.0166
$$

である。

内側の水平端は

$$
R_-=\sqrt{\left(3\right)^{2} - \left(\frac{1}{10}\right)}\simeq 2.9833
$$

である。

水平半幅は

$$
h_R=- \frac{\left(\frac{\sqrt{890}}{10}\right)}{2} + \frac{\left(\frac{\sqrt{910}}{10}\right)}{2}\simeq 0.016667
$$

である。

軸の大半径上での垂直半幅は

$$
h_Z=\frac{\sqrt{\frac{\left(\frac{1}{100}\right)}{\left(4\right)}}}{\left(3\right)}=\frac{1}{60}
$$

である。

軸近傍楕円の伸長度は

$$
\kappa_{\rm lin}=2 \sqrt{\frac{\left(1\right)}{\left(4\right)}}=1
$$

である。

水平と垂直の半幅がほぼそろう。
:::

:::example[軸近くの縦長面]
磁束係数を変えて伸長度を求める。

入力値を次のように置く。

$$
\begin{aligned}
R_0&=3 \\A&=1 \\C&=1 \\\psi_a&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

磁気面の水平端での$\abs{R^2-R_0^2}=D$ は

$$
D=\sqrt{\frac{\left(\frac{1}{100}\right)}{\left(1\right)}}=\frac{1}{10}
$$

である。

外側の水平端は

$$
R_+=\sqrt{\left(3\right)^{2} + \left(\frac{1}{10}\right)}\simeq 3.0166
$$

である。

内側の水平端は

$$
R_-=\sqrt{\left(3\right)^{2} - \left(\frac{1}{10}\right)}\simeq 2.9833
$$

である。

水平半幅は

$$
h_R=- \frac{\left(\frac{\sqrt{890}}{10}\right)}{2} + \frac{\left(\frac{\sqrt{910}}{10}\right)}{2}\simeq 0.016667
$$

である。

軸の大半径上での垂直半幅は

$$
h_Z=\frac{\sqrt{\frac{\left(\frac{1}{100}\right)}{\left(1\right)}}}{\left(3\right)}=\frac{1}{30}
$$

である。

軸近傍楕円の伸長度は

$$
\kappa_{\rm lin}=2
$$

である。

軸近傍の二次の伸長度は二になる。
:::

:::example[軸から遠い磁気面]
幅が大きいときの内外端の距離を調べる。

入力値を次のように置く。

$$
\begin{aligned}
R_0&=3 \\A&=1 \\C&=1 \\\psi_a&=16
\end{aligned}
$$

これらを順に代入する。

磁気面の水平端での$\abs{R^2-R_0^2}=D$ は

$$
D=\sqrt{\frac{\left(16\right)}{\left(1\right)}}=4
$$

である。

外側の水平端は

$$
R_+=\sqrt{\left(3\right)^{2} + \left(4\right)}\simeq 3.6056
$$

である。

内側の水平端は

$$
R_-=\sqrt{\left(3\right)^{2} - \left(4\right)}\simeq 2.2361
$$

である。

水平半幅は

$$
h_R=\frac{\left(\sqrt{13}\right)}{2} - \frac{\left(\sqrt{5}\right)}{2}\simeq 0.68474
$$

である。

軸の大半径上での垂直半幅は

$$
h_Z=\frac{\sqrt{\frac{\left(16\right)}{\left(1\right)}}}{\left(3\right)}=\frac{4}{3}
$$

である。

軸近傍楕円の伸長度は

$$
\kappa_{\rm lin}=2
$$

である。

軸の両側に等しい距離ではなく、円形近似の誤差が増える。
:::

:::example[横長の磁気面]
垂直の曲率を大きくして断面がどう変わるか計算する。

入力値を次のように置く。

$$
\begin{aligned}
R_0&=3 \\A&=1 \\C&=16 \\\psi_a&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

磁気面の水平端での$\abs{R^2-R_0^2}=D$ は

$$
D=\sqrt{\frac{\left(\frac{1}{100}\right)}{\left(1\right)}}=\frac{1}{10}
$$

である。

外側の水平端は

$$
R_+=\sqrt{\left(3\right)^{2} + \left(\frac{1}{10}\right)}\simeq 3.0166
$$

である。

内側の水平端は

$$
R_-=\sqrt{\left(3\right)^{2} - \left(\frac{1}{10}\right)}\simeq 2.9833
$$

である。

水平半幅は

$$
h_R=- \frac{\left(\frac{\sqrt{890}}{10}\right)}{2} + \frac{\left(\frac{\sqrt{910}}{10}\right)}{2}\simeq 0.016667
$$

である。

軸の大半径上での垂直半幅は

$$
h_Z=\frac{\sqrt{\frac{\left(\frac{1}{100}\right)}{\left(16\right)}}}{\left(3\right)}=\frac{1}{120}
$$

である。

軸近傍楕円の伸長度は

$$
\kappa_{\rm lin}=2 \sqrt{\frac{\left(1\right)}{\left(16\right)}}=\frac{1}{2}
$$

である。

同じ磁束に達する垂直距離が短くなる。
:::


## 安全係数と磁気シア

一周のポロイダル移動の間に磁力線がトロイダルに回る角を $2\pi q$ とし、$q$ を :term[安全係数]{en="safety factor" yomi="あんぜんけいすう"} と呼ぶ。トロイダル角の増分は $\dd\varphi=B_\varphi\dd l_p/(R B_p)$ なので
$$
q(\psi)=\frac1{2\pi}\oint\frac{B_\varphi}{R B_p}\dd l_p.
\label{eq-plasma-safety-factor}
$$
積分では $\dd l_p>0$ をポロイダル磁場に沿う弧長とする。この規約で $q$ の符号はトロイダル磁場の向きで決まる。$\psi$ が外向きに増える入れ子状の面では、磁気面間の薄い環の幅は $\dd\psi/|\nabla\psi|$ である。従ってトロイダル磁束 $\Phi_t=\int B_\varphi\,\dd R\dd Z$ の磁気面微分は $\dd\Phi_t/\dd\psi=\oint B_\varphi\dd l_p/|\nabla\psi|=2\pi q$ になる。大アスペクト比でほぼ円形なら $q\simeq rB_\varphi/(R_0B_\theta)$ であり、局所 :term[磁気シア]{en="magnetic shear" yomi="じきしあ"} は $\hat s=(r/q)\dd q/\dd r$ である。位相 $m\theta-n\varphi$ が磁力線に沿って一定となる面は $q=m/n$ である。

:::example[軸電流が一様な円柱模型]
$B_\theta\propto r$、$B_\varphi$ 一定なら $q$ は定数、$\hat s=0$ である。$B_\theta\propto 1/r$ の電流のない外部では $q\propto r^2$、$\hat s=2$ となる。
:::

:::example[有理面]
$q=3/2$ では磁力線は二回のポロイダル周回で三回のトロイダル周回を終える。$m=3,n=2$ の摂動はこの線上で位相が変わらない。
:::

:::example[電流を測って辺の安全係数を出す]
大アスペクト比の円形断面で $B_\theta(a)\simeq\mu_0 I_p/(2\pi a)$ なら $q(a)\simeq 2\pi a^2 B_\varphi/(\mu_0 R_0 I_p)$ である。半径と外部磁場を保ち、電流を二倍にすると辺の $q$ は半分となる。磁気シアには内部の電流分布も必要であり、全電流だけでは決まらない。
:::

### 局所電流と平均電流から磁気シアを得る

大アスペクト比円筒で $q=rB_z/(R_0B_\theta)=2\pi B_z r^2/(\mu_0R_0 I)$、$B_z$ を一定とする。$\hat s=\dd\ln|q|/\dd\ln r=2-rI'/I=2(1-J_z/\bar J_z)$ である。$\bar J_z=I/(\pi r^2)$ と $I'=2\pi rJ_z$ を使った。$J_z=J_0[1-\alpha(r/a)^2]$ では $I=\pi J_0 r^2[1-\alpha(r/a)^2/2]$ なので、電流が外側ほど低いと正のシアになる。安全係数は局所電流密度だけでは決まらず、内側の総電流を要する。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $x$ | 小半径を辺半径で割った $r/a$ | 1 |
| $\alpha$ | 電流密度の曲率係数 | 1 |
| $q_0$ | 軸での安全係数 | 1 |

局所の電流密度比を $J_z/J_0$ とする。二次の電流分布を使う。

$$
J_z/J_0=- \alpha x^{2} + 1
$$

平均電流密度比を $\bar J_z/J_0$ とする。円板積分の重みを含める。

$$
\bar J_z/J_0=- \frac{\alpha x^{2}}{2} + 1
$$

安全係数を $q(r)$ とする。総電流から磁場のピッチを求める。

$$
q(r)=\frac{q_{0}}{(\bar J_z/J_0)}
$$

磁気シアを $\hat s$ とする。局所と平均の電流差を測る。

$$
\hat s=- \frac{2 (J_z/J_0)}{(\bar J_z/J_0)} + 2
$$

電流の内側割合を $I(r)/I(a)$ とする。同じ分布の全電流で規格化する。

$$
I(r)/I(a)=\frac{(\bar J_z/J_0) x^{2}}{1 - \frac{\alpha}{2}}
$$

回転変換を $\iota$ とする。この軸対称規約では安全係数の逆数である。

$$
\iota=\frac{1}{q(r)}
$$

:::example[一様電流]
半径を変えても安全係数にシアが生じない条件を確認する。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{1}{2} \\\alpha&=0 \\q_0&=2
\end{aligned}
$$

これらを順に代入する。

局所の電流密度比は

$$
J_z/J_0=- \left(0\right) \left(\frac{1}{2}\right)^{2} + 1=1
$$

である。

平均電流密度比は

$$
\bar J_z/J_0=- \frac{\left(0\right) \left(\frac{1}{2}\right)^{2}}{2} + 1=1
$$

である。

安全係数は

$$
q(r)=\frac{\left(2\right)}{\left(1\right)}=2
$$

である。

磁気シアは

$$
\hat s=0
$$

である。

電流の内側割合は

$$
I(r)/I(a)=\frac{\left(1\right) \left(\frac{1}{2}\right)^{2}}{1 - \frac{\left(0\right)}{2}}=\frac{1}{4}
$$

である。

回転変換は

$$
\iota=\frac{1}{\left(2\right)}=\frac{1}{2}
$$

である。

$I\propto r^2$ のため $q$ は一定である。
:::

:::example[辺で消える電流密度]
辺の局所電流がゼロでも総電流が流れることを確かめる。

入力値を次のように置く。

$$
\begin{aligned}
x&=1 \\\alpha&=1 \\q_0&=2
\end{aligned}
$$

これらを順に代入する。

局所の電流密度比は

$$
J_z/J_0=1 - \left(1\right)^{3}=0
$$

である。

平均電流密度比は

$$
\bar J_z/J_0=1 - \frac{\left(1\right)^{3}}{2}=\frac{1}{2}
$$

である。

安全係数は

$$
q(r)=\frac{\left(2\right)}{\left(\frac{1}{2}\right)}=4
$$

である。

磁気シアは

$$
\hat s=- \frac{2 \left(0\right)}{\left(\frac{1}{2}\right)} + 2=2
$$

である。

電流の内側割合は

$$
I(r)/I(a)=\frac{\left(1\right)^{2} \left(\frac{1}{2}\right)}{1 - \frac{\left(1\right)}{2}}=1
$$

である。

回転変換は

$$
\iota=\frac{1}{\left(4\right)}=\frac{1}{4}
$$

である。

$q$ は有限でシアは二になる。
:::

:::example[中間半径の有理近傍]
内側電流の割合と $q$ を求める。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{1}{2} \\\alpha&=1 \\q_0&=2
\end{aligned}
$$

これらを順に代入する。

局所の電流密度比は

$$
J_z/J_0=- \left(1\right) \left(\frac{1}{2}\right)^{2} + 1=\frac{3}{4}
$$

である。

平均電流密度比は

$$
\bar J_z/J_0=- \frac{\left(1\right) \left(\frac{1}{2}\right)^{2}}{2} + 1=\frac{7}{8}
$$

である。

安全係数は

$$
q(r)=\frac{\left(2\right)}{\left(\frac{7}{8}\right)}=\frac{16}{7}
$$

である。

磁気シアは

$$
\hat s=- \frac{2 \left(\frac{3}{4}\right)}{\left(\frac{7}{8}\right)} + 2=\frac{2}{7}
$$

である。

電流の内側割合は

$$
I(r)/I(a)=\frac{\left(\frac{1}{2}\right)^{2} \left(\frac{7}{8}\right)}{1 - \frac{\left(1\right)}{2}}=\frac{7}{16}
$$

である。

回転変換は

$$
\iota=\frac{1}{\left(\frac{16}{7}\right)}=\frac{7}{16}
$$

である。

辺と比べて軸の値に近く、シアも小さい。
:::

:::example[外側ほど電流が高い場合]
電流勾配の反転でシアの符号がどう変わるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
x&=1 \\\alpha&=-1 \\q_0&=2
\end{aligned}
$$

これらを順に代入する。

局所の電流密度比は

$$
J_z/J_0=- \left(-1\right) \left(1\right)^{2} + 1=2
$$

である。

平均電流密度比は

$$
\bar J_z/J_0=- \frac{\left(-1\right) \left(1\right)^{2}}{2} + 1=\frac{3}{2}
$$

である。

安全係数は

$$
q(r)=\frac{\left(2\right)}{\left(\frac{3}{2}\right)}=\frac{4}{3}
$$

である。

磁気シアは

$$
\hat s=- \frac{2 \left(2\right)}{\left(\frac{3}{2}\right)} + 2=- \frac{2}{3}
$$

である。

電流の内側割合は

$$
I(r)/I(a)=\frac{\left(1\right)^{2} \left(\frac{3}{2}\right)}{1 - \frac{\left(-1\right)}{2}}=1
$$

である。

回転変換は

$$
\iota=\frac{1}{\left(\frac{4}{3}\right)}=\frac{3}{4}
$$

である。

局所が平均を上回るので負シアになる。
:::


## トロイダル平衡とトカマク

トカマクでは外部コイルが主にトロイダル磁場を作り、プラズマのトロイダル電流がポロイダル磁場を作る。真空の $B_\varphi=F/R$ は外側で弱くなる。磁場の曲率と強度勾配は粒子を径方向にずらしうるが、磁力線が上下を回れば粒子は異なるドリフト方向を経験する。平衡全体はこの局所描像では決まらず、式 :eqref[eq-plasma-grad-shafranov] に圧力、電流関数、境界磁束を与えて解く。

### 圧力と磁場の尺度

:term[プラズマベータ]{en="plasma beta" yomi="ぷらずまべーた"} は $\beta=2\mu_0p/B^2$ である。トロイダル磁場に対する $\beta$ とポロイダル磁場に対する $\beta_p$ は異なる。$B_\theta/B_\varphi\ll1$ なら低い $\beta$ でも $\beta_p$ は大きくなりうる。磁気軸の外向き変位を定量化する近似では、圧力と電流の分布、外部磁場を指定して式 :eqref[eq-plasma-grad-shafranov] をアスペクト比で展開する。

:::example[辺の安全係数]
$R_0=3\,\mathrm m$、$a=1\,\mathrm m$、$B_\varphi=3\,\mathrm T$、$B_\theta(a)=0.5\,\mathrm T$ では円形近似で $q(a)=2$ である。
:::

:::example[二つのベータ]
$B_\theta/B_\varphi=0.1$ なら、同じ局所圧力で $\beta_p/\beta_\varphi=100$ である。
:::

:::example[電流とピッチの尺度]
$a/R_0=1/4$、$q(a)=2$ なら $B_\theta/B_\varphi=1/8$ である。同じ局所圧力でも $\beta_p\simeq64\beta_\varphi$ となる。力の平衡や不安定性の計算で、どちらを用いた数値かを確かめる。
:::

### 円形境界と磁気軸のずれ

磁気軸の位置は、外から与えた境界の中心と一致するとは限らない。例として $x=R-R_0$、$r=\sqrt{x^2+Z^2}$ とし、$r=a$ に $\psi=0$ を与える。$a/R_0\ll1$ で $p'$ と $FF'$ を定数とする。方程式の右辺を $S_0+S_1 x+\order{\mu_0|p'|a^2}$ と展開すると $S_0=-\mu_0R_0^2p'-FF'$、$S_1=-2\mu_0R_0p'$ である。$S_0\ne0$ とし、$S_1$ が $S_0/R_0$ 程度の場合を考える。

最低次は極断面の Laplace 作用素 $\nabla_p^2=\partial_x^2+\partial_Z^2$ に対して $\psi_0=S_0(r^2-a^2)/4$ となる。次の次数では $-R^{-1}\partial_R$ の項を右辺に移し、
$$
\nabla_p^2\psi_1=(S_1+S_0/(2R_0))x
$$
を得る。$\nabla_p^2[x(r^2-a^2)]=8x$ であり、境界でこの積は消える。従って $\psi_1=(S_1+S_0/(2R_0))x(r^2-a^2)/8$ である。軸で $\partial_x(\psi_0+\psi_1)=0$ を満たす点は、先頭の次数で
$$
x_{
m axis}=\frac{a^2}{4}\left(\frac{S_1}{S_0}+\frac{1}{2R_0}\right)
$$
となる。例えば $FF'=0$、$p'<0$ の場合は $S_1/S_0=2/R_0$ なので、磁気軸の変位は外向きの $5a^2/(8R_0)$ である。これは $p'$、$FF'$ と円形境界を指定した模型での値であり、実機で軸の位置を求めるには電流分布と外部磁場の境界条件が要る。

### 円形断面の有限アスペクト比補正

円形のポロイダル断面を $R=R_0+r\cos\theta$ と書く。特定の面での幾何学的比較として $B_\theta$ を一定に置き、$B_\varphi=F/R$ とすると
$$
q=\frac{rF}{2\pi B_\theta}\int_0^{2\pi}\frac{\dd\theta}{(R_0+r\cos\theta)^2}
=\frac{rF R_0}{B_\theta(R_0^2-r^2)^{3/2}} .
$$
積分は $\tan(\theta/2)$ で $\int\dd\theta/(R_0+r\cos\theta)=2\pi/\sqrt{R_0^2-r^2}$ を求め、それを $R_0$ で微分して得る。$F=B_\varphi(R_0)R_0$ とすれば大アスペクト比の式に対する倍率は $(1-(r/R_0)^2)^{-3/2}$ である。ただし、隣り合う円形面にわたって $B_r=0$ かつ $B_\theta$ 一定という仮定を課すと、発散は $(rR)^{-1}\partial_\theta(RB_\theta)$ となり、一般にゼロにならない。発散ゼロを満たす同心円モデルでは $RB_\theta$ を面上で一定に取るため、補正の係数が変わる。次の角積分でこの違いを計算する。

:::example[内外の角度の重み]
$r/R_0=1/3$ なら倍率は $(9/8)^{3/2}\simeq1.19$ である。$R$ の小さい側では同じ $B_\theta$ でも $\varphi$ の変化が大きくなり、積分への寄与が増える。
:::

### トロイダル磁力線の角積分

円形近似の磁気面 $R=R_0+r\cos\theta$ で、トロイダル磁場を $B_\varphi=B_0R_0/R$、ポロイダル磁場を $B_\theta=B_{p0}R_0/R$ とする。線の接線成分の比から $\dd\varphi/\dd\theta=rB_\varphi/(RB_\theta)$ である。$\int_0^{2\pi}\dd\theta/(1+\epsilon\cos\theta)=2\pi/\sqrt{1-\epsilon^2}$ は半角置換 $t=\tan(\theta/2)$ による $2\int_{\R}\dd t/[(1+\epsilon)+(1-\epsilon)t^2]$ で得られる。局所のピッチは内側が大きいが、平均の補正は一次では消える。この同心円模型は $\nabla\cdot\vb B=0$ を満たす。ただし、圧力平衡も満たすかは Grad–Shafranov 方程式に別途代入して調べる。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $\epsilon$ | 半径比 $r/R_0<1$ | 1 |
| $q_{\rm cyl}$ | 円柱近似での安全係数 $rB_0/(R_0B_{p0})$ | 1 |

内側における局所のピッチを $q_{\rm in}$ とする。$R$ の小さい位置に対応する。

$$
q_{\rm in}=\frac{q_{\rm cyl}}{1 - \epsilon}
$$

外側における局所のピッチを $q_{\rm out}$ とする。同じ小半径で $R$ が最大になる。

$$
q_{\rm out}=\frac{q_{\rm cyl}}{\epsilon + 1}
$$

角度平均の安全係数を $q$ とする。ポロイダル一周当たりのトロイダル周回数である。

$$
q=\frac{q_{\rm cyl}}{\sqrt{1 - \epsilon^{2}}}
$$

円柱値に対する補正を $q/q_{\rm cyl}-1$ とする。平均で一次の項は打ち消す。

$$
q/q_{\rm cyl}-1=\frac{q}{q_{\rm cyl}} - 1
$$

内外の局所ピッチ比を $q_{\rm in}/q_{\rm out}$ とする。局所の変化には一次項が残る。

$$
q_{\rm in}/q_{\rm out}=\frac{q_{\rm in}}{q_{\rm out}}
$$

一周当たりのトロイダル角を $\Delta\varphi$ とする。周回数を角に変換する。

$$
\Delta\varphi=2 \pi q
$$

:::example[大アスペクト比]
局所差と平均補正を比べる。

入力値を次のように置く。

$$
\begin{aligned}
\epsilon&=\frac{1}{10} \\q_{\rm cyl}&=3
\end{aligned}
$$

これらを順に代入する。

内側における局所のピッチは

$$
q_{\rm in}=\frac{\left(3\right)}{1 - \left(\frac{1}{10}\right)}=\frac{10}{3}
$$

である。

外側における局所のピッチは

$$
q_{\rm out}=\frac{\left(3\right)}{\left(\frac{1}{10}\right) + 1}=\frac{30}{11}
$$

である。

角度平均の安全係数は

$$
q=\frac{\left(3\right)}{\sqrt{1 - \left(\frac{1}{10}\right)^{2}}}\simeq 3.0151
$$

である。

円柱値に対する補正は

$$
q/q_{\rm cyl}-1=-1 + \frac{\left(\frac{10 \sqrt{11}}{11}\right)}{\left(3\right)}\simeq 0.0050378
$$

である。

内外の局所ピッチ比は

$$
q_{\rm in}/q_{\rm out}=\frac{\left(\frac{10}{3}\right)}{\left(\frac{30}{11}\right)}=\frac{11}{9}
$$

である。

一周当たりのトロイダル角は

$$
\Delta\varphi=2 \pi \left(\frac{10 \sqrt{11}}{11}\right)\simeq 18.945
$$

である。

局所差は一割の尺度だが平均補正は一パーセント未満である。
:::

:::example[やや太い円環]
円柱近似から平均がどれだけ変わるか求める。

入力値を次のように置く。

$$
\begin{aligned}
\epsilon&=\frac{1}{3} \\q_{\rm cyl}&=3
\end{aligned}
$$

これらを順に代入する。

内側における局所のピッチは

$$
q_{\rm in}=\frac{\left(3\right)}{1 - \left(\frac{1}{3}\right)}=\frac{9}{2}
$$

である。

外側における局所のピッチは

$$
q_{\rm out}=\frac{\left(3\right)}{\left(\frac{1}{3}\right) + 1}=\frac{9}{4}
$$

である。

角度平均の安全係数は

$$
q=\frac{\left(3\right)}{\sqrt{1 - \left(\frac{1}{3}\right)^{2}}}\simeq 3.182
$$

である。

円柱値に対する補正は

$$
q/q_{\rm cyl}-1=-1 + \frac{\left(\frac{9 \sqrt{2}}{4}\right)}{\left(3\right)}\simeq 0.06066
$$

である。

内外の局所ピッチ比は

$$
q_{\rm in}/q_{\rm out}=\frac{\left(\frac{9}{2}\right)}{\left(\frac{9}{4}\right)}=2
$$

である。

一周当たりのトロイダル角は

$$
\Delta\varphi=2 \pi \left(\frac{9 \sqrt{2}}{4}\right)\simeq 19.993
$$

である。

内側と外側の角度の重みの差が残る。
:::

:::example[小アスペクト比の形式値]
幾何補正が大きい場合を計算する。

入力値を次のように置く。

$$
\begin{aligned}
\epsilon&=\frac{3}{5} \\q_{\rm cyl}&=3
\end{aligned}
$$

これらを順に代入する。

内側における局所のピッチは

$$
q_{\rm in}=\frac{\left(3\right)}{1 - \left(\frac{3}{5}\right)}=\frac{15}{2}
$$

である。

外側における局所のピッチは

$$
q_{\rm out}=\frac{\left(3\right)}{\left(\frac{3}{5}\right) + 1}=\frac{15}{8}
$$

である。

角度平均の安全係数は

$$
q=\frac{\left(3\right)}{\sqrt{1 - \left(\frac{3}{5}\right)^{2}}}=\frac{15}{4}
$$

である。

円柱値に対する補正は

$$
q/q_{\rm cyl}-1=-1 + \frac{\left(\frac{15}{4}\right)}{\left(3\right)}=\frac{1}{4}
$$

である。

内外の局所ピッチ比は

$$
q_{\rm in}/q_{\rm out}=\frac{\left(\frac{15}{2}\right)}{\left(\frac{15}{8}\right)}=4
$$

である。

一周当たりのトロイダル角は

$$
\Delta\varphi=2 \pi \left(\frac{15}{4}\right)\simeq 23.562
$$

である。

円形の仮定での値であり、実際の平衡には形状を解く必要がある。
:::

:::example[同じ幾何でピッチ変更]
ポロイダル場を半分にしたときの平均を求める。

入力値を次のように置く。

$$
\begin{aligned}
\epsilon&=\frac{1}{3} \\q_{\rm cyl}&=6
\end{aligned}
$$

これらを順に代入する。

内側における局所のピッチは

$$
q_{\rm in}=\frac{\left(6\right)}{1 - \left(\frac{1}{3}\right)}=9
$$

である。

外側における局所のピッチは

$$
q_{\rm out}=\frac{\left(6\right)}{\left(\frac{1}{3}\right) + 1}=\frac{9}{2}
$$

である。

角度平均の安全係数は

$$
q=\frac{\left(6\right)}{\sqrt{1 - \left(\frac{1}{3}\right)^{2}}}\simeq 6.364
$$

である。

円柱値に対する補正は

$$
q/q_{\rm cyl}-1=-1 + \frac{\left(\frac{9 \sqrt{2}}{2}\right)}{\left(6\right)}\simeq 0.06066
$$

である。

内外の局所ピッチ比は

$$
q_{\rm in}/q_{\rm out}=\frac{\left(9\right)}{\left(\frac{9}{2}\right)}=2
$$

である。

一周当たりのトロイダル角は

$$
\Delta\varphi=2 \pi \left(\frac{9 \sqrt{2}}{2}\right)\simeq 39.986
$$

である。

安全係数は二倍となるが内外比は変わらない。
:::


## ステラレータとヘリカル系

外部コイルの三次元配置でも磁力線にねじれを与えられる。これを使う閉じ込め装置を :term[ステラレータ]{en="stellarator" yomi="すてられーた"} と呼ぶ。軸対称の $q$ と対応する量として :term[回転変換]{en="rotational transform" yomi="かいてんへんかん"} $\iota=\lim_{N\to\infty}\Delta\theta/(2\pi N)$ を定める。ここで $N$ はトロイダル周回数である。軸対称の直線磁力線角では $\iota=1/q$ である。

### 磁気面と三次元構造

磁場のトロイダル成分が同符号で、同じ断面を毎周横切る場合に交点を描く。滑らかな入れ子状の磁気面に載った交点は、その面と断面の交線上にある。近い有理面で共鳴する磁場成分は島を作ることがあり、断面に複数の閉曲線の組が現れる。このように、同じ断面で交点を順に取る操作は :term[Poincaré 写像]{en="Poincaré map" yomi="ぽあんかれしゃぞう"} である。写像が保存する面積は幾何学的な断面積そのものではなく、断面法線方向の磁束 $B_n\dd S$ である。実際、二つの断面を同じ細い磁束管の側面でつなぐと、側面を通る磁束はゼロであり、$\nabla\cdot\vb B=0$ から端面の磁束は等しい。保存写像と共鳴の一般的な扱いは :cite[goldston-rutherford] を参照する。以下では孤立した島の幅を求める。

ステラレータでは外部コイルで回転変換を与えるため、トカマクと異なり、変換の主成分をトロイダルプラズマ電流に頼らずに作れる。連続したらせん巻線を使うヘリオトロンなどのヘリカル系も、断面が変化する複数のコイルを使う装置も、この点は共通である。有限圧力では圧力勾配を支える電流が流れるので、真空磁気面とプラズマ平衡には一般に差がある。粒子の軌道には磁気面の形に加えて磁場強度の分布も効く。強度の変化は、粒子が磁力線に沿って通過するか、反射して往復するかを決める。

:::example[断面交点]
理想化して一周ごとに $\theta$ が $2\pi(2/5)$ 進むなら、幾何学的な角度は $0,4\pi/5,8\pi/5,2\pi/5,6\pi/5$ を巡り、五周目で戻る。
:::

:::example[圧力の制約]
磁力線が一つの体積領域で稠密なら、連続な $p$ が線に沿って一定であることからその体積全体で同じ値になる。この領域の内部で有限の滑らかな圧力勾配を指定することはできない。
:::

### 有理面近傍の島の断面模型

磁力線の交点を局所的に共鳴位相 $\theta_h$ と磁束からのずれ $X$ で表す。場に沿う変数を時間とし、単一共鳴の Hamiltonian を $H=a X^2/2-\epsilon\cos\theta_h$、$a,\epsilon>0$ と近似する。これは回転数のずれを $aX$、共鳴成分を正弦波とした模型であり、$\dot\theta_h=aX,\ \dot X=-\epsilon\sin\theta_h$ を与える。例えば $\theta_h=m\theta-n\varphi$ とすれば、固定した $\varphi$ の断面で同じ安定位相が $\theta$ の一周に $|m|$ 回現れる。単一の有理共鳴が作る島の個数はこの位相の反復数に対応する。鞍点 $\theta_h=\pi,X=0$ の $H=\epsilon$ が分離曲線である。$\theta_h=0$ での半幅は $2\sqrt{\epsilon/a}$ となる。複数の共鳴が重なると単一の保存量 $H$ では断面全体を表せない。

:::example[シアによる幅の変化]
共鳴成分を保ち $a$ を四倍にすれば、島の半幅は半分である。$\epsilon$ を四倍に増やす場合は幅が二倍となる。
:::

### 磁力線の共鳴島を振り子で評価する

磁力線をトロイダル角 $\varphi$ に従って追い、共鳴位相 $\alpha=m\theta-n\varphi$ と磁束からの変位 $x$ を正準変数に取る。共鳴近傍で $\dot\alpha=g x$、$\dot x=-\epsilon\sin\alpha$ と線形シアと単一高調波を残すと $H=gx^2/2-\epsilon\cos\alpha$ が保存する。$g,\epsilon>0$ なら中心は $(0,0)$、鞍点は $(\pi,0)$ である。鞍点の値 $H=\epsilon$ を中心位相で評価して半幅 $2\sqrt{\epsilon/g}$ を得る。中心の近くでは $\ddot\alpha=-g\epsilon\alpha$ なのでトロイダル角当たりの往復率は $\sqrt{g\epsilon}$ である。複数の高調波が重なる場合には各島を独立に扱えない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $g$ | 規格化した磁気シア | 1 |
| $\epsilon$ | 共鳴高調波の振幅 | 1 |
| $x$ | 磁束方向の規格化距離 | 1 |
| $\alpha$ | 共鳴位相 | 1 |

島の磁束半幅を $w_{1/2}$ とする。分離軌道の中心位相で評価する。

$$
w_{1/2}=2 \sqrt{\frac{\epsilon}{g}}
$$

島の全幅を $w_{\rm full}$ とする。対称な両側を含める。

$$
w_{\rm full}=2 w_{1/2}
$$

小振幅の往復率を $\omega_b$ とする。中心を線形化して求める。

$$
\omega_b=\sqrt{\epsilon g}
$$

規格化 Hamiltonianを $H/\epsilon$ とする。一と比較して閉じた軌道か判定する。

$$
H/\epsilon=- \cos{\left(\alpha \right)} + \frac{g x^{2}}{2 \epsilon}
$$

位相の進み率を $\dot\alpha$ とする。シアにより共鳴面から離れると変わる。

$$
\dot\alpha=g x
$$

磁束方向への移動率を $\dot x$ とする。位相によって符号が変わる。

$$
\dot x=- \epsilon \sin{\left(\alpha \right)}
$$

:::example[分離軌道上の点]
半幅の端のエネルギーを確かめる。

入力値を次のように置く。

$$
\begin{aligned}
g&=1 \\\epsilon&=\frac{1}{100} \\x&=\frac{1}{5} \\\alpha&=0
\end{aligned}
$$

これらを順に代入する。

島の磁束半幅は

$$
w_{1/2}=2 \sqrt{\frac{\left(\frac{1}{100}\right)}{\left(1\right)}}=\frac{1}{5}
$$

である。

島の全幅は

$$
w_{\rm full}=2 \left(\frac{1}{5}\right)=\frac{2}{5}
$$

である。

小振幅の往復率は

$$
\omega_b=\sqrt{\left(1\right) \left(\frac{1}{100}\right)}=\frac{1}{10}
$$

である。

規格化 Hamiltonianは

$$
H/\epsilon=\frac{\left(1\right) \left(\frac{1}{5}\right)^{2}}{2 \left(\frac{1}{100}\right)} - \cos{\left(\left(0\right) \right)}=1
$$

である。

位相の進み率は

$$
\dot\alpha=\left(1\right) \left(\frac{1}{5}\right)=\frac{1}{5}
$$

である。

磁束方向への移動率は

$$
\dot x=- \left(\frac{1}{100}\right) \sin{\left(\left(0\right) \right)}=0
$$

である。

$H/\epsilon=1$ であり、中心周りの閉曲線との境界である。
:::

:::example[島中心の近く]
閉じた軌道かどうか調べる。

入力値を次のように置く。

$$
\begin{aligned}
g&=1 \\\epsilon&=\frac{1}{100} \\x&=\frac{1}{10} \\\alpha&=\frac{\pi}{3}
\end{aligned}
$$

これらを順に代入する。

島の磁束半幅は

$$
w_{1/2}=2 \sqrt{\frac{\left(\frac{1}{100}\right)}{\left(1\right)}}=\frac{1}{5}
$$

である。

島の全幅は

$$
w_{\rm full}=2 \left(\frac{1}{5}\right)=\frac{2}{5}
$$

である。

小振幅の往復率は

$$
\omega_b=\sqrt{\left(1\right) \left(\frac{1}{100}\right)}=\frac{1}{10}
$$

である。

規格化 Hamiltonianは

$$
H/\epsilon=\frac{\left(1\right) \left(\frac{1}{10}\right)^{2}}{2 \left(\frac{1}{100}\right)} - \cos{\left(\left(\frac{\pi}{3}\right) \right)}=0
$$

である。

位相の進み率は

$$
\dot\alpha=\left(1\right) \left(\frac{1}{10}\right)=\frac{1}{10}
$$

である。

磁束方向への移動率は

$$
\dot x=- \left(\frac{1}{100}\right) \sin{\left(\left(\frac{\pi}{3}\right) \right)}\simeq -0.0086603
$$

である。

エネルギーは鞍点より小さく、位相は往復する。
:::

:::example[強いシア]
同じ点と摂動を保って島の被覆を求める。

入力値を次のように置く。

$$
\begin{aligned}
g&=4 \\\epsilon&=\frac{1}{100} \\x&=\frac{1}{5} \\\alpha&=0
\end{aligned}
$$

これらを順に代入する。

島の磁束半幅は

$$
w_{1/2}=2 \sqrt{\frac{\left(\frac{1}{100}\right)}{\left(4\right)}}=\frac{1}{10}
$$

である。

島の全幅は

$$
w_{\rm full}=2 \left(\frac{1}{10}\right)=\frac{1}{5}
$$

である。

小振幅の往復率は

$$
\omega_b=\sqrt{\left(4\right) \left(\frac{1}{100}\right)}=\frac{1}{5}
$$

である。

規格化 Hamiltonianは

$$
H/\epsilon=\frac{\left(4\right) \left(\frac{1}{5}\right)^{2}}{2 \left(\frac{1}{100}\right)} - \cos{\left(\left(0\right) \right)}=7
$$

である。

位相の進み率は

$$
\dot\alpha=\left(4\right) \left(\frac{1}{5}\right)=\frac{4}{5}
$$

である。

磁束方向への移動率は

$$
\dot x=- \left(\frac{1}{100}\right) \sin{\left(\left(0\right) \right)}=0
$$

である。

島幅は半分に縮み、局所の往復率は倍になる。与えた点は島の外にある。
:::

:::example[弱い共鳴成分]
摂動振幅を百倍弱めたときに幅を比較する。

入力値を次のように置く。

$$
\begin{aligned}
g&=1 \\\epsilon&=\frac{1}{10000} \\x&=\frac{1}{5} \\\alpha&=\frac{\pi}{2}
\end{aligned}
$$

これらを順に代入する。

島の磁束半幅は

$$
w_{1/2}=2 \sqrt{\frac{\left(\frac{1}{10000}\right)}{\left(1\right)}}=\frac{1}{50}
$$

である。

島の全幅は

$$
w_{\rm full}=2 \left(\frac{1}{50}\right)=\frac{1}{25}
$$

である。

小振幅の往復率は

$$
\omega_b=\sqrt{\left(1\right) \left(\frac{1}{10000}\right)}=\frac{1}{100}
$$

である。

規格化 Hamiltonianは

$$
H/\epsilon=\frac{\left(1\right) \left(\frac{1}{5}\right)^{2}}{2 \left(\frac{1}{10000}\right)} - \cos{\left(\left(\frac{\pi}{2}\right) \right)}=200
$$

である。

位相の進み率は

$$
\dot\alpha=\left(1\right) \left(\frac{1}{5}\right)=\frac{1}{5}
$$

である。

磁束方向への移動率は

$$
\dot x=- \left(\frac{1}{10000}\right) \sin{\left(\left(\frac{\pi}{2}\right) \right)}=- \frac{1}{10000}
$$

である。

幅は十倍狭まり、線形応答と島幅は同じ次数ではない。
:::


## 演習問題

::::exercise{level=basic}
直線磁場で磁気圧が $3\,\mathrm{Pa}$ 増すとき熱圧力はどう変わるか。

:::solution
全圧一定なので熱圧力は $3\,\mathrm{Pa}$ 減る。
:::
::::

::::exercise{level=basic}
$q=4$ に対応する $\iota$ を求めよ。

:::solution
トロイダル一周あたりのポロイダル回転は $1/q=1/4$ である。
:::
::::

::::exercise{level=basic}
同一半径で Z ピンチの電流を二倍にすると積分圧力差は何倍か。

:::solution
式 :eqref[eq-plasma-bennett] の右辺は $I^2$ に比例するので四倍である。
:::
::::

::::exercise{level=basic}
$q(r)=1+2(r/a)^2$ の辺のシアを求めよ。

:::solution
$q(a)=3$、$a q'(a)=4$ から $\hat s(a)=4/3$ である。
:::
::::

::::exercise{level=advanced}
純 $\theta$ ピンチで $B_z=B_0(r/a)^2$ と置き、$p(r)-p(a)$ を求めよ。

:::solution
張力がゼロなので $p+B_z^2/(2\mu_0)$ は一定である。辺との差は $B_0^2[1-(r/a)^4]/(2\mu_0)$ となる。
:::
::::

::::exercise{level=advanced}
$p'=0$、$FF'=g$ が一定のとき $\psi$ の $R$ だけに依存する特解を求めよ。

:::solution
$\partial_R(R^{-1}\psi')=-g/R$ なので $\psi'=-gR\ln(R/R_*)$ とできる。さらに積分して $\psi=-gR^2\ln(R/R_*)/2+gR^2/4$ を得る。
:::
::::

::::exercise{level=advanced}
円柱近似で軸電流密度が $j_z\propto r^2$ のとき $q(r)$ とシアのべき依存性を求めよ。

:::solution
包む電流は $\int_0^r j_z(r')2\pi r'\dd r'\propto r^4$ で、$B_\theta\propto r^3$ である。軸磁場一定なら $q\propto r^{-2}$、$r>0$ で $\hat s=-2$ である。
:::
::::

::::exercise{level=advanced}
純 Z ピンチの $B_\theta=br$ について Bennett 関係の両辺を直接比較せよ。

:::solution
圧力差の積分は $2\pi b^2\int_0^a r(a^2-r^2)\dd r/\mu_0=\pi b^2a^4/(2\mu_0)$ である。$I=2\pi b a^2/\mu_0$ を右辺に入れると同じ値になる。
:::
::::

## 確認テスト

::::quiz
静的平衡で必ず成り立つものを選べ。

- [x] 磁力線に沿って圧力は一定である。
- [ ] 電流密度は常にゼロである。
- [ ] 磁場強度は必ず空間一定である。

:::explanation
$\nabla p=\vb j\times\vb B$ と $\vb B$ の内積はゼロである。従って $\vb B\cdot\nabla p=0$ である。
:::
::::

::::quiz{answer="2" tol="1e-6"}
$r/R_0=1/4,\ B_\varphi/B_\theta=8$ の円形近似で $q$ を答えよ。

:::explanation
積 $q\simeq (r/R_0)(B_\varphi/B_\theta)=8/4=2$ である。
:::
::::

::::quiz
純 Z ピンチの曲率力の向きはどれか。

- [x] 軸へ向かう。
- [ ] 軸から遠ざかる。
- [ ] 円筒の周方向を向く。

:::explanation
方位角磁力線の曲率は $-\hat{\vb r}/r$ であるため、張力は動径内向きである。
:::
::::

::::quiz{answer="3" tol="1e-6"}
$q\propto r^3$ の $\hat s$ を答えよ。

:::explanation
$\hat s=\dd\ln|q|/\dd\ln r$ なのでべき指数そのものの $3$ である。
:::
::::

---
title: "有限群の表現と物理"
description: "有限対称群で状態を分類し、分子振動、選択則、結晶場、時間反転による縮退を導く。"
notation:
  - symbol: 'D(g)'
    meaning: "群元 g の表現行列"
  - symbol: '\chi(g)'
    meaning: "表現の指標"
  - symbol: 'm_\alpha'
    meaning: "既約型 α の重複度"
  - symbol: '\Theta'
    meaning: "時間反転を表す反ユニタリ作用素"
---

## この章の目標

- 対称操作を群と表現で記述する。
- 指標の直交関係から既約成分を取り出す。
- 点群で振動モードと遷移の可否を判定する。
- 結晶場と時間反転で生じる縮退の違いを計算する。

[表現の指標](/math/representation-theory/characters/) を、分子や結晶の対称操作に適用する。状態の縮退と、状態間の遷移を許す条件は別の計算で求める。

## 対称性と群

### 対称操作を合成する

有限集合 $G$ に結合的な積があり、単位元と各元の逆元が存在するとき :term[有限群]{en="finite group" yomi="ゆうげんぐん"} という。例えば正三角形を自分自身へ移す平面内の回転 $r$ と鏡映 $s$ は $r^3=s^2=e$、$srs=r^{-1}$ を満たす。$sr=r^{-1}s$ によって任意の積を $r^a$ または $r^as$（$a=0,1,2$）に直せる。三頂点への作用が六通り異なるので群の位数は $6$ である。平面内の操作を法線方向を固定する三次元の操作に延長した群は $C_{3v}$ である。元 $hgh^{-1}$ を $g$ の共役と呼び、この群の共役類は $\{e\},\{r,r^2\},\{s,rs,r^2s\}$ である。鏡映軸は回転で互いに移り、回転の向きは鏡映で反転するためである。

Hamiltonian $H$ の線形な空間対称操作は、定義域を保つユニタリ写像 $D(g)$ で $D(g)HD(g)^{-1}=H$ を満たすものとする。$Hv=Ev$ なら $H D(g)v=E D(g)v$ なので、エネルギーの固有空間を群の作用で分類できる。半整数スピンに対する固有回転を扱う場合には、点群の表現に代えて [連続群と物理](/physics/mathematical-physics/lie-groups-physics/) で述べる二重被覆のもとでの点群の逆像を使う。

### 生成元の関係から操作を整理する

正 $m$ 角形の回転生成元を $R$、鏡映を $S$ とすると $R^m=S^2=I$ であり、行列積を計算すると $SRS=R^{-1}$ になる。どんな積もこの関係で $R^j$ または $R^jS$ に整理できる。前者は行列式 $+1$、後者は $-1$ で区別でき、同じ側では $0\le j<m$ の角だけが異なる。従って $m\ge3$ のとき、これらが相異なる $2m$ 個の操作である。

:::example[平面の二面体群（1）]
平面で正多角形の回転と一つの鏡映を行列で表す。 三角形の操作を取る。

$$
m=3
$$

一つ進む回転角を選ぶ。

$$
\begin{aligned}
R\big|_{m=3} &= \begin{pmatrix}- \frac{1}{2} & - \frac{\sqrt{3}}{2}\\\frac{\sqrt{3}}{2} & - \frac{1}{2}\end{pmatrix}
\end{aligned}
$$

第 2 軸を反転して、回転との積を作る。

$$
\begin{aligned}
RS\big|_{m=3} &= \begin{pmatrix}- \frac{1}{2} & \frac{\sqrt{3}}{2}\\\frac{\sqrt{3}}{2} & \frac{1}{2}\end{pmatrix}
\end{aligned}
$$

鏡映を左右に掛けて回転を逆にする。

$$
\begin{aligned}
SRS\big|_{m=3} &= \begin{pmatrix}- \frac{1}{2} & \frac{\sqrt{3}}{2}\\- \frac{\sqrt{3}}{2} & - \frac{1}{2}\end{pmatrix}
\end{aligned}
$$

行列式の符号を含めると群の位数が決まる。

$$
\begin{aligned}
|D_m|\big|_{m=3} &= 6
\end{aligned}
$$
:::

:::example[平面の二面体群（2）]
平面で正多角形の回転と一つの鏡映を行列で表す。 正方形の操作を取る。

$$
m=4
$$

一つ進む回転角を選ぶ。

$$
\begin{aligned}
R\big|_{m=4} &= \begin{pmatrix}0 & -1\\1 & 0\end{pmatrix}
\end{aligned}
$$

第 2 軸を反転して、回転との積を作る。

$$
\begin{aligned}
RS\big|_{m=4} &= \begin{pmatrix}0 & 1\\1 & 0\end{pmatrix}
\end{aligned}
$$

鏡映を左右に掛けて回転を逆にする。

$$
\begin{aligned}
SRS\big|_{m=4} &= \begin{pmatrix}0 & 1\\-1 & 0\end{pmatrix}
\end{aligned}
$$

行列式の符号を含めると群の位数が決まる。

$$
\begin{aligned}
|D_m|\big|_{m=4} &= 8
\end{aligned}
$$
:::

:::example[平面の二面体群（3）]
平面で正多角形の回転と一つの鏡映を行列で表す。 正六角形の操作を取る。

$$
m=6
$$

一つ進む回転角を選ぶ。

$$
\begin{aligned}
R\big|_{m=6} &= \begin{pmatrix}\frac{1}{2} & - \frac{\sqrt{3}}{2}\\\frac{\sqrt{3}}{2} & \frac{1}{2}\end{pmatrix}
\end{aligned}
$$

第 2 軸を反転して、回転との積を作る。

$$
\begin{aligned}
RS\big|_{m=6} &= \begin{pmatrix}\frac{1}{2} & \frac{\sqrt{3}}{2}\\\frac{\sqrt{3}}{2} & - \frac{1}{2}\end{pmatrix}
\end{aligned}
$$

鏡映を左右に掛けて回転を逆にする。

$$
\begin{aligned}
SRS\big|_{m=6} &= \begin{pmatrix}\frac{1}{2} & \frac{\sqrt{3}}{2}\\- \frac{\sqrt{3}}{2} & \frac{1}{2}\end{pmatrix}
\end{aligned}
$$

行列式の符号を含めると群の位数が決まる。

$$
\begin{aligned}
|D_m|\big|_{m=6} &= 12
\end{aligned}
$$
:::

:::example[平面の二面体群（4）]
平面で正多角形の回転と一つの鏡映を行列で表す。 正八角形の操作を取る。

$$
m=8
$$

一つ進む回転角を選ぶ。

$$
\begin{aligned}
R\big|_{m=8} &= \begin{pmatrix}\frac{\sqrt{2}}{2} & - \frac{\sqrt{2}}{2}\\\frac{\sqrt{2}}{2} & \frac{\sqrt{2}}{2}\end{pmatrix}
\end{aligned}
$$

第 2 軸を反転して、回転との積を作る。

$$
\begin{aligned}
RS\big|_{m=8} &= \begin{pmatrix}\frac{\sqrt{2}}{2} & \frac{\sqrt{2}}{2}\\\frac{\sqrt{2}}{2} & - \frac{\sqrt{2}}{2}\end{pmatrix}
\end{aligned}
$$

鏡映を左右に掛けて回転を逆にする。

$$
\begin{aligned}
SRS\big|_{m=8} &= \begin{pmatrix}\frac{\sqrt{2}}{2} & \frac{\sqrt{2}}{2}\\- \frac{\sqrt{2}}{2} & \frac{\sqrt{2}}{2}\end{pmatrix}
\end{aligned}
$$

行列式の符号を含めると群の位数が決まる。

$$
\begin{aligned}
|D_m|\big|_{m=8} &= 16
\end{aligned}
$$
:::

### 群が保つ軌道と固定部分

群 $G$ が集合 $X$ に作用するとき、$x$ の軌道を $Gx$、固定部分群を $G_x=\{g:gx=x\}$ とする。$g x=h x$ と $h^{-1}g\in G_x$ は同値なので、左剰余類 $gG_x$ を $gx$ に送る写像が一対一対応になる。したがって有限の場合に $|Gx|=|G|/|G_x|$ を得る。位置の軌道と、作用を受けても変わらない方向を分けると、変位表現の計算に移りやすい。

正三角形の頂点が有向の一周に並ぶとし、$r=(123)$、$s=(23)$ とする。頂点 $1$ には $r^0,r,r^2$ が別の像を与える一方、$s$ は位置を保つ。頂点を変位のない印として扱う場合には固定条件はこれで尽きる。頂点 $1$ に接する平面内ベクトルを追加すると、$s$ の鏡映軸に沿うベクトルは引き続き固定されるが、その法線方向の非零ベクトルは符号が変わる。後者では位置と方向の組の固定部分群は自明になり、軌道は六点になる。


:::example[三角形の置換行列と平均]
基底の三成分を三頂点の変位とする。全群を平均すると、すべての頂点で同じ値をもつ変位だけが残る。

座標成分を列に並べる。この順序に合わせ、線形写像の行列の列を選ぶ。

$$
\begin{aligned}v&=\begin{pmatrix}3\\-1\\1\end{pmatrix}\\ A&=\begin{pmatrix}\frac{1}{3} & \frac{1}{3} & \frac{1}{3}\\\frac{1}{3} & \frac{1}{3} & \frac{1}{3}\\\frac{1}{3} & \frac{1}{3} & \frac{1}{3}\end{pmatrix}\end{aligned}
$$

行列の各行と入力の内積を計算する。例えば第一行では $\begin{pmatrix}\frac{1}{3} & \frac{1}{3} & \frac{1}{3}\end{pmatrix}$ と入力列を掛ける。

$$
\begin{aligned}(Av)_1&=1\\Av&=\begin{pmatrix}1\\1\\1\end{pmatrix}\end{aligned}
$$

核の次元は列数から階数を引いて得られる。像の次元は階数である。

$$
\begin{aligned}\operatorname{rank}A&=1\\\dim\ker A&=2\end{aligned}
$$

正方行列の場合は、特性多項式を求めて固有値を調べられる。ここでは $tI-A$ の行列式を取る。

$$
\begin{aligned}\operatorname{tr}A&=1\\ \det A&=0\\ \det(tI-A)&=t^{3} - t^{2}\end{aligned}
$$

行列要素を得るため、出力から指定した共変ベクトルで成分を取り出す。

$$
\begin{aligned}\ell&=\begin{pmatrix}1 & 1 & 1\end{pmatrix}\\ \ell Av&=3\end{aligned}
$$

次に別の方向を入力にする。出力の差は入力の差に $A$ を掛けたものである。

$$
\begin{aligned}w&=\begin{pmatrix}2\\-2\\0\end{pmatrix}\\ Aw&=\begin{pmatrix}0\\0\\0\end{pmatrix}\\ Aw-Av&=\begin{pmatrix}-1\\-1\\-1\end{pmatrix}\end{aligned}
$$
:::

## 表現と指標による既約分解

### 指標の直交性から成分を数える

:term[表現]{en="representation" yomi="ひょうげん"} は $D(e)=I$、$D(gh)=D(g)D(h)$ を満たす可逆な線形写像の割当である。有限群では任意の正定値内積を $\sum_g\inner{D(g)u}{D(g)v}/|G|$ で平均すると不変になる。よって不変部分空間の直交補も不変であり、次元について帰納的に既約空間へ分解できる。既約とは非自明な不変部分空間がないことをいう。:term[指標]{en="character" yomi="しひょう"} は $\chi_D(g)=\tr D(g)$ であり、共役類上で一定で、直和に対して加法的である。

:::proposition[指標と重複度]{#prop-mp-character-multiplicity}
複素既約表現 $\alpha$ が $D$ に含まれる重複度は
$$
m_\alpha=\frac1{|G|}\sum_{g\in G}\overline{\chi_\alpha(g)}\,\chi_D(g).
\label{eq-mp-character-multiplicity}
$$
:::

:::proof
二表現の間の線形写像 $A$ を $\Pi(A)=|G|^{-1}\sum_gD(g)A D_\alpha(g)^{-1}$ で平均すると、任意の群元に対し可換な写像に射影する。既約表現間の非零可換写像は、核と像の不変性から同型である。また同じ既約空間の可換写像は、複素数上では固有値を一つ持ち、その固有空間が不変なのでスカラーである。従って射影の階数は $m_\alpha$ である。行列単位で計算すると写像 $A\mapsto BAC$ のトレースは $\tr B\,\tr C$ であり、$D_\alpha(g)^{-1}$ の指標はユニタリ性により $\overline{\chi_\alpha(g)}$ である。よって平均のトレースを取ると式を得る。
:::

$C_{3v}$ の自明表現を $A_1$、鏡映だけに負号をつける一次元表現を $A_2$、平面上の表現を $E$ とする。

| 既約型 | $e$ | $r,r^2$ | $s,rs,r^2s$ |
| --- | ---: | ---: | ---: |
| $A_1$ | 1 | 1 | 1 |
| $A_2$ | 1 | 1 | -1 |
| $E$ | 2 | -1 | 0 |

$E$ は回転の二つの固有直線が鏡映で交換されるので既約である。任意の複素既約表現は正則表現の中に次元の回数だけ現れる。実際、正則表現の指標は恒等元で $|G|$、他で零なので :eqref[eq-mp-character-multiplicity] を使える。表の次元の二乗和が $1+1+4=6$ で既に全次元を使うので、他の既約型はない。例えば三頂点の置換表現の指標は固定点数 $(3,0,1)$ で、表との内積から $A_1+E$ に分かれる。

### 既約型の状態だけを残す射影

指標の式を成分の抽出にも使える。$d_\alpha=\chi_\alpha(e)$ とし、
$$
P_\alpha=\frac{d_\alpha}{|G|}\sum_g\overline{\chi_\alpha(g)}D(g)
$$
と置く。係数が共役類上で一定なので $P_\alpha$ は群と交換する。各既約型 $\beta$ 上ではスカラーで、そのトレースは指標の直交性により $d_\alpha\delta_{\alpha\beta}$ である。従って $\beta=\alpha$ で恒等、それ以外で零になる。例えば三頂点の置換表現の $A_1$ 射影は各成分を三成分の平均 $(x_1+x_2+x_3)/3$ に置き換える。$E$ 射影はその平均を引く。$(2,1,0)$ は $(1,1,1)+(1,0,-1)$ に分かれる。

### 表現の積と指標での重複度

回転と鏡映を持つ $C_{3v}$ の既約指標を $A_1=(1,1,1),A_2=(1,1,-1),E=(2,-1,0)$ とする。テンソル積の指標は指標どうしの積である。これは行列の直積の対角成分が元の対角成分の積になるためである。各共役類の大きさ $1,2,3$ を掛けて内積を取れば、同じ既約型の成分を数えられる。

:::example[二次元表現のテンソル冪（1）]
$E^{\otimes m}$ の重複度を指標から求める。

$$
m=1
$$

各共役類で指標の積を取る。

$$
\begin{aligned}
\chi\big|_{m=1} &= \left( 2, \  -1, \  0\right)
\end{aligned}
$$

重み付き内積で二つの一次元表現を数える。

$$
\begin{aligned}
n(A_1)=n(A_2)\big|_{m=1} &= 0
\end{aligned}
$$

二次元表現の指標を掛けると回転類の符号が変わる。

$$
\begin{aligned}
n(E)\big|_{m=1} &= 1
\end{aligned}
$$

次元に戻して合計が元の直積の次元になる。

$$
\begin{aligned}
\dim E^{\otimes m}\big|_{m=1} &= 2
\end{aligned}
$$
:::

:::example[二次元表現のテンソル冪（2）]
$E^{\otimes m}$ の重複度を指標から求める。

$$
m=2
$$

各共役類で指標の積を取る。

$$
\begin{aligned}
\chi\big|_{m=2} &= \left( 4, \  1, \  0\right)
\end{aligned}
$$

重み付き内積で二つの一次元表現を数える。

$$
\begin{aligned}
n(A_1)=n(A_2)\big|_{m=2} &= 1
\end{aligned}
$$

二次元表現の指標を掛けると回転類の符号が変わる。

$$
\begin{aligned}
n(E)\big|_{m=2} &= 1
\end{aligned}
$$

次元に戻して合計が元の直積の次元になる。

$$
\begin{aligned}
\dim E^{\otimes m}\big|_{m=2} &= 4
\end{aligned}
$$
:::

:::example[二次元表現のテンソル冪（3）]
$E^{\otimes m}$ の重複度を指標から求める。

$$
m=3
$$

各共役類で指標の積を取る。

$$
\begin{aligned}
\chi\big|_{m=3} &= \left( 8, \  -1, \  0\right)
\end{aligned}
$$

重み付き内積で二つの一次元表現を数える。

$$
\begin{aligned}
n(A_1)=n(A_2)\big|_{m=3} &= 1
\end{aligned}
$$

二次元表現の指標を掛けると回転類の符号が変わる。

$$
\begin{aligned}
n(E)\big|_{m=3} &= 3
\end{aligned}
$$

次元に戻して合計が元の直積の次元になる。

$$
\begin{aligned}
\dim E^{\otimes m}\big|_{m=3} &= 8
\end{aligned}
$$
:::

:::example[二次元表現のテンソル冪（4）]
$E^{\otimes m}$ の重複度を指標から求める。

$$
m=4
$$

各共役類で指標の積を取る。

$$
\begin{aligned}
\chi\big|_{m=4} &= \left( 16, \  1, \  0\right)
\end{aligned}
$$

重み付き内積で二つの一次元表現を数える。

$$
\begin{aligned}
n(A_1)=n(A_2)\big|_{m=4} &= 3
\end{aligned}
$$

二次元表現の指標を掛けると回転類の符号が変わる。

$$
\begin{aligned}
n(E)\big|_{m=4} &= 5
\end{aligned}
$$

次元に戻して合計が元の直積の次元になる。

$$
\begin{aligned}
\dim E^{\otimes m}\big|_{m=4} &= 16
\end{aligned}
$$
:::

:::example[対称二階テンソルの指標]
平面の $E$ 型変位を二つ掛けて対称化する。固有値 $\lambda_i$ の組 $\lambda_i\lambda_j$ は $i\le j$ の範囲で数える。したがって指標は $(\chi(g)^2+\chi(g^2))/2$ である。

共役類の大きさと対象の指標を同じ順で並べる。

$$
\begin{aligned}
(n_C)_C &= \left( 1, \  2, \  3\right) \\
(\chi(C))_C &= \left( 3, \  0, \  1\right)
\end{aligned}
$$

A_1 について、各共役類からの寄与を並べて足す。

$$
\begin{aligned}
\chi_{A_1} &= \left( 1, \  1, \  1\right) \\
(n_C\chi(C)\overline{\chi_{A_1}(C)})_C &= \left( 3, \  0, \  3\right) \\
n_{A_1} &= 1
\end{aligned}
$$

A_2 について、各共役類からの寄与を並べて足す。

$$
\begin{aligned}
\chi_{A_2} &= \left( 1, \  1, \  -1\right) \\
(n_C\chi(C)\overline{\chi_{A_2}(C)})_C &= \left( 3, \  0, \  -3\right) \\
n_{A_2} &= 0
\end{aligned}
$$

E について、各共役類からの寄与を並べて足す。

$$
\begin{aligned}
\chi_{E} &= \left( 2, \  -1, \  0\right) \\
(n_C\chi(C)\overline{\chi_{E}(C)})_C &= \left( 6, \  0, \  0\right) \\
n_{E} &= 1
\end{aligned}
$$

各重複度に既約表現の次元を掛け、対象全体の次元を回収する。

$$
\begin{aligned}
\sum_\alpha n_\alpha d_\alpha &= 3 \\
|G| &= 6
\end{aligned}
$$

対称二階表現には不変なトレースと、二つのトレース零成分がある。残る一次元の交代部分は面の向きに対応する。
:::

## 点群と結晶点群

### 格子を保つ対称性

原点を固定する直交操作の有限群を :term[点群]{en="point group" yomi="てんぐん"} と呼ぶ。結晶点群では格子の基底で操作が整数行列となる。三次元の固有回転は固有値 $1,e^{i\theta},e^{-i\theta}$ を持つため $\tr=1+2\cos\theta$ が整数である。従って $2\cos\theta$ は $-2,-1,0,1,2$ のいずれかで、回転の位数は $1,2,3,4,6$ に限られる。例えば五回軸の $2\cos(2\pi/5)=(\sqrt5-1)/2$ は整数でない。一方、正五角形分子の点群には格子保存の仮定がないので五回軸を持てる。鏡映は一つの座標の符号を変え、反転は三座標すべての符号を変える。両者の行列式はともに $-1$ だが、固定点集合は異なる。

### 結晶系ごとの点群

Hermann–Mauguin の記法では $n$ は対応する軸の $n$ 回回転、$\bar n$ は回転と中心反転の合成、$m$ は鏡映、$n/m$ はその軸に垂直な鏡面を伴うことを示す。記号の後続の位置は系によって決まる別方向の対称性を指定する。例えば $4mm$ の主軸は四回軸で、縦の座標面と対角面の両方が鏡面である。以下は三次元の点群型を国際記号の短い形で並べたものである。

| 結晶系 | 点群型 | 個数 |
| --- | --- | ---: |
| 三斜 | $1$、$\bar1$ | 2 |
| 単斜 | $2$、$m$、$2/m$ | 3 |
| 斜方 | $222$、$mm2$、$mmm$ | 3 |
| 正方 | $4$、$\bar4$、$4/m$、$422$、$4mm$、$\bar42m$、$4/mmm$ | 7 |
| 三方 | $3$、$\bar3$、$32$、$3m$、$\bar3m$ | 5 |
| 六方 | $6$、$\bar6$、$6/m$、$622$、$6mm$、$\bar62m$、$6/mmm$ | 7 |
| 立方 | $23$、$m\bar3$、$432$、$\bar43m$、$m\bar3m$ | 5 |

例えば $3m$ は上の $C_{3v}$、$mm2$ は非直線三原子分子の例で用いる $C_{2v}$、$m\bar3m$ は結晶場の例で用いる $O_h$ である。純回転だけの $432$ は立方体を回して重ねる $24$ 個の操作からなる。これは座標軸の置換が $3!$、符号の組は行列式 $+1$ を満たす $2^2$ 通りであるためである。中心反転を掛けたものを合わせると $m\bar3m$ の $48$ 個になる。

### なぜ三次元では 32 種類になるのか

まず有限回転群 $H\subset SO(3)$ の位数を $N>1$ とする。非恒等元の回転軸が単位球面と交わる点を集め、各軌道の代表の固定部分群の位数を $m_1,\ldots,m_s$ とする。固定部分群は一つの軸の回転からなる巡回群である。非恒等回転とその固定点の組を数えると、一つの回転が二点を固定するので
$$
2(N-1)=\sum_{i=1}^s\frac{N}{m_i}(m_i-1)
$$
である。各 $m_i\ge2$ なので $s$ は $2$ または $3$ である。$s=2$ なら $m_1=m_2=N$ となり、共通軸を持つ巡回群である。$s=3$ なら、昇順に並べた $m_i$ に対して $\sum_i1/m_i=1+2/N>1$ である。従って候補は $(2,2,n)$、$(2,3,3)$、$(2,3,4)$、$(2,3,5)$ であり、最後の型は格子を保てない。

$(2,2,n)$ では $N=2n$ で、位数 $n$ の軸の固定部分群は指数 $2$ で正規である。残る元は軸の二端を交換するので、軸に垂直な半回転と巡回群からなる $D_n$ となる。$n=2$ では三つの半回転の軸は互いに垂直である。$(2,3,3)$ では三回軸の端の一軌道に四点がある。各点の三回回転が他の三点を巡回し、群が四点に推移的に作用するので、全ての辺長は等しい。従って群は正四面体の $12$ 個の回転 $T$ である。$(2,3,4)$ では四回軸の端の軌道に六点がある。一点を固定する四回回転の下で、同じ軸の反対側の一点と四点の巡回に分かれる。六点の集合は群の推移性により対蹠点を全て含む。残る四点は回転で同じ高さにあり、対蹠点も同じ高さにあるので、その高さは零である。従って四点は軸に垂直な大円で正方形を作る。これは正八面体の回転群 $O$ であり、その位数は $24$ である。

従って格子を保つ回転群の型は $C_1,C_2,C_3,C_4,C_6,D_2,D_3,D_4,D_6,T,O$ の $11$ 個である。$C_n,D_n$ は正方形、三角形または長方形の格子を軸方向に積層して実現でき、$T,O$ は立方格子に含まれる。$O(3)$ の任意の元は固有回転と中心反転の積に一意に分かれる。中心反転は任意の格子を保つので、一般の点群から中心反転の因子を除いて得る回転群も、この $11$ 個のいずれかである。

中心反転を含む点群は $H\times\{I,-I\}$ であり、これで $11$ 型を得る。中心反転を含まなければ各回転に一つずつ符号が付くので、群は準同型 $H\to\{\pm1\}$ のグラフである。自明な符号では固有回転群そのものの $11$ 型を得る。非自明な符号は偶数位数の巡回群三つに一つずつある。$D_n=\langle r,s\rangle$ では $r^n=s^2=e$、$srs=r^{-1}$ に従い、奇数 $n$ なら $(r,s)$ の符号は $(+,-)$ の一通りである。$D_2$ では三つの軸を入れ替えられるため一通りである。$D_4,D_6$ では主軸の半刻みの回転による共役が $s$ を $rs$ に移すため、$r$ の符号に応じて二通りずつとなる。

四面体群 $T$ は頂点の偶置換群と同型である。回転は四面体の向きを保つため偶置換になり、作用は忠実で位数も $12$ である。偶置換は三巡回で生成されるから、$T$ は非自明な符号を持たない。立方体群 $O$ は四本の空間対角線を置換する。この置換が恒等なら各対角線には $+1$ か $-1$ が掛かる。異なる線は互いに直交しないため全て同じ符号になり、固有回転では $+1$ だけが可能である。従ってこの置換作用は忠実で、位数が $24$ なので全置換群を得る。その互換は全て共役で群を生成するため、非自明な符号は置換の符号一つである。合計は $11+11+(3+1+1+2+2+1)=32$ となる。

### 結晶格子で整数行列を作る

格子基底で回転を表したとき、各基底ベクトルの像は整数係数で書ける。回転の逆も格子を保つため、整数行列の行列式は $\pm1$ である。平面で向きを保つ回転のトレース は $2\cos(2\pi/m)$ であり、これが整数であることから $m=1,2,3,4,6$ が候補になる。結晶で許される候補には、実際に正方格子や三角格子で構成できる例がある。

:::example[基底が作る格子と回転の行列（1）]
第二基底を角度 $2\pi/m$ だけ傾けて、回転が格子上でどう見えるか調べる。 三回対称を調べる。

$$
m=3
$$

列が格子基底である行列を作る。

$$
\begin{aligned}
B\big|_{m=3} &= \begin{pmatrix}1 & - \frac{1}{2}\\0 & \frac{\sqrt{3}}{2}\end{pmatrix}
\end{aligned}
$$

二つの列を回転した結果を座標表示する。

$$
\begin{aligned}
RB\big|_{m=3} &= \begin{pmatrix}- \frac{1}{2} & - \frac{1}{2}\\\frac{\sqrt{3}}{2} & - \frac{\sqrt{3}}{2}\end{pmatrix}
\end{aligned}
$$

基底へ戻すと回転が格子を保つか判定できる。

$$
\begin{aligned}
B^{-1}RB\big|_{m=3} &= \begin{pmatrix}0 & -1\\1 & -1\end{pmatrix}
\end{aligned}
$$
:::

:::example[基底が作る格子と回転の行列（2）]
第二基底を角度 $2\pi/m$ だけ傾けて、回転が格子上でどう見えるか調べる。 四回対称を調べる。

$$
m=4
$$

列が格子基底である行列を作る。

$$
\begin{aligned}
B\big|_{m=4} &= \begin{pmatrix}1 & 0\\0 & 1\end{pmatrix}
\end{aligned}
$$

二つの列を回転した結果を座標表示する。

$$
\begin{aligned}
RB\big|_{m=4} &= \begin{pmatrix}0 & -1\\1 & 0\end{pmatrix}
\end{aligned}
$$

基底へ戻すと回転が格子を保つか判定できる。

$$
\begin{aligned}
B^{-1}RB\big|_{m=4} &= \begin{pmatrix}0 & -1\\1 & 0\end{pmatrix}
\end{aligned}
$$
:::

:::example[基底が作る格子と回転の行列（3）]
第二基底を角度 $2\pi/m$ だけ傾けて、回転が格子上でどう見えるか調べる。 五回回転を試す。

$$
m=5
$$

列が格子基底である行列を作る。

$$
\begin{aligned}
B\big|_{m=5} &= \begin{pmatrix}1 & - \frac{1}{4} + \frac{\sqrt{5}}{4}\\0 & \frac{\sqrt{2 \sqrt{5} + 10}}{4}\end{pmatrix}
\end{aligned}
$$

二つの列を回転した結果を座標表示する。

$$
\begin{aligned}
RB\big|_{m=5} &= \begin{pmatrix}- \frac{1}{4} + \frac{\sqrt{5}}{4} & - \frac{\sqrt{5}}{4} - \frac{1}{4}\\\frac{\sqrt{2 \sqrt{5} + 10}}{4} & \frac{\sqrt{10 - 2 \sqrt{5}}}{4}\end{pmatrix}
\end{aligned}
$$

基底へ戻すと回転が格子を保つか判定できる。

$$
\begin{aligned}
B^{-1}RB\big|_{m=5} &= \begin{pmatrix}0 & -1\\1 & - \frac{1}{2} + \frac{\sqrt{5}}{2}\end{pmatrix}
\end{aligned}
$$
:::

:::example[基底が作る格子と回転の行列（4）]
第二基底を角度 $2\pi/m$ だけ傾けて、回転が格子上でどう見えるか調べる。 六回対称を調べる。

$$
m=6
$$

列が格子基底である行列を作る。

$$
\begin{aligned}
B\big|_{m=6} &= \begin{pmatrix}1 & \frac{1}{2}\\0 & \frac{\sqrt{3}}{2}\end{pmatrix}
\end{aligned}
$$

二つの列を回転した結果を座標表示する。

$$
\begin{aligned}
RB\big|_{m=6} &= \begin{pmatrix}\frac{1}{2} & - \frac{1}{2}\\\frac{\sqrt{3}}{2} & \frac{\sqrt{3}}{2}\end{pmatrix}
\end{aligned}
$$

基底へ戻すと回転が格子を保つか判定できる。

$$
\begin{aligned}
B^{-1}RB\big|_{m=6} &= \begin{pmatrix}0 & -1\\1 & 1\end{pmatrix}
\end{aligned}
$$
:::

### 立方体の固有回転を分類する

立方体を座標軸に置くと、固有回転は行列式が $1$ の符号付き置換行列で尽くされる。三本の軸を置換する方法は $3!$ 通りで、それぞれに行列式の条件を満たす符号の組が四つある。従って位数は $24$ である。回転の軸を調べると、面中心を結ぶ三軸について四分の一回転が六つ、半回転が三つ、対頂点を結ぶ四軸について三分の一回転が八つ、対辺の中点を結ぶ六軸について半回転が六つある。これと恒等を足して $24$ となる。

半回転の二類は同じ回転角を持つが、立方体の対称操作で面中心を辺の中点へ移せないため互いに共役でない。一方で同じ種類の軸は符号付き置換によって互いに移せる。回転方向の違いも軸を反転させる別の固有回転で結ばれるので、上の五種類が共役類になる。中心反転 $-I$ を加えると、それはすべての回転と可換であり、回転部分と反転の直積になる。このとき各既約型に偶・奇をそれぞれ付けられる。


:::example[立方体の空間対角線まわりの操作]
三つの正の座標軸を順に入れ替える。対角線方向は固定されるが、一般のベクトルは三回操作して初めて戻る。

座標成分を列に並べる。この順序に合わせ、線形写像の行列の列を選ぶ。

$$
\begin{aligned}v&=\begin{pmatrix}1\\2\\3\end{pmatrix}\\ A&=\begin{pmatrix}0 & 0 & 1\\1 & 0 & 0\\0 & 1 & 0\end{pmatrix}\end{aligned}
$$

行列の各行と入力の内積を計算する。例えば第一行では $\begin{pmatrix}0 & 0 & 1\end{pmatrix}$ と入力列を掛ける。

$$
\begin{aligned}(Av)_1&=3\\Av&=\begin{pmatrix}3\\1\\2\end{pmatrix}\end{aligned}
$$

核の次元は列数から階数を引いて得られる。像の次元は階数である。

$$
\begin{aligned}\operatorname{rank}A&=3\\\dim\ker A&=0\end{aligned}
$$

正方行列の場合は、特性多項式を求めて固有値を調べられる。ここでは $tI-A$ の行列式を取る。

$$
\begin{aligned}\operatorname{tr}A&=0\\ \det A&=1\\ \det(tI-A)&=t^{3} - 1\end{aligned}
$$

行列要素を得るため、出力から指定した共変ベクトルで成分を取り出す。

$$
\begin{aligned}\ell&=\begin{pmatrix}1 & 1 & 1\end{pmatrix}\\ \ell Av&=6\end{aligned}
$$

次に別の方向を入力にする。出力の差は入力の差に $A$ を掛けたものである。

$$
\begin{aligned}w&=\begin{pmatrix}1\\1\\1\end{pmatrix}\\ Aw&=\begin{pmatrix}1\\1\\1\end{pmatrix}\\ Aw-Av&=\begin{pmatrix}-2\\0\\-1\end{pmatrix}\end{aligned}
$$
:::

## 分子振動の基準モードの分類

### 原子の置換と変位を同時に数える

平衡にある $N$ 個の原子の微小変位を並べ、質量重み付き座標 $q_a=\sqrt{m_a}\,\delta r_a$ を取る。対称操作は等しい質量の原子を置換し、変位の三成分にも直交行列 $O_g$ を掛ける。表現のトレースへ寄与するのは固定された原子の対角ブロックだけなので $\chi_{3N}(g)=N_{\rm fixed}(g)\tr O_g$ である。安定な平衡点で調和近似を取ると力定数行列は質量重み付き座標で実対称で対称操作と交換する。各固有値は $\omega^2$（$\mathrm{s^{-2}}$）であり、分子全体の並進・回転を取り除くと振動の分類になる。

交換関係は調和近似の二次エネルギー $V_2(q)=q^{\mathsf T}Kq/2$ から従う。変位の表現を $D(g)$ とすると $V_2(D(g)q)=V_2(q)$ なので、対称行列について $D(g)^{\mathsf T}KD(g)=K$ が全ての $g$ で成り立つ。さらに $D(g)$ が直交だから $KD(g)=D(g)K$ である。運動エネルギーは質量重み付き座標で $\sum_a|\dot q_a|^2/2$ なので、正規直交固有基底で振動方程式は $\ddot Q_i+\omega_i^2 Q_i=0$ に分かれる。

自由な分子では全体の並進と回転の接ベクトルは $K$ の零空間に入る。実際、対称操作で移した平衡点も平衡点なので、その経路に沿ってポテンシャルの勾配を微分すると、Hessian が接ベクトルを零に送る。重心を原点とすれば並進と回転の接空間は、$\sum_a m_a \vb r_a=0$ により直交する。非直線分子では三本の回転方向が独立で $3N-6$ モード、直線分子では分子軸まわりの回転が全原子を固定するので $3N-5$ モードとなる。内部変形にも零固有値がある場合、モード数にはその零方向を含める。

:::example[三回軸のある四原子分子]
中心一原子と三つの等価な周辺原子が非平面の $C_{3v}$ 分子を作るとする。恒等元では $12$、回転では固定数 $1$ だが $\tr O_r=0$ のため $0$、鏡映では固定数 $2$、$\tr O_s=1$ なので $2$ である。並進は極性ベクトルなので $A_1+E$、回転は軸性ベクトルなので鏡映に行列式を掛け $A_2+E$ となる。差の振動指標は $(6,0,2)$ である。式 :eqref[eq-mp-character-multiplicity] から $m_{A_1}=(6+3\cdot2)/6=2$、$m_{A_2}=0$、$m_E=(12+0)/6=2$ であり、振動は $2A_1+2E$ である。非直線分子の $3N-6=6$ と一致する。
:::

:::example[水型の三原子非直線分子]
分子面を $xz$、二回軸を $z$ とする。$C_{2v}$ の操作を $(e,C_2,\sigma_{xz},\sigma_{yz})$ の順に取る。恒等で原子は三つ固定、$C_2$ と $\sigma_{yz}$ で中心だけ固定、$\sigma_{xz}$ で三つ固定である。各直交行列のトレースを掛けると三次元変位全体の指標は $(9,-1,3,1)$ である。並進の指標は $(3,-1,1,1)$、軸性の回転は $(3,-1,-1,-1)$ で、差は $(3,1,3,1)$ になる。$A_1=(1,1,1,1)$ と $x$ 成分の表現 $(1,-1,1,-1)$ に対する内積は順に $2,1$ である。従って二つの対称な振動と一つの $x$ と同型の振動を得る。三つとも双極子の $z$ または $x$ 成分と同じ型なので、基底が自明型なら対称性による赤外禁制は無い。
:::

### 固定点から変位表現を数える

ある対称操作が原子を移す場合、その原子に対応する三つの行列の対角成分は零である。固定した原子だけが空間ベクトルのトレース を寄与する。例えば水平面内の配置に対し、面の鏡映では全原子が固定され、各原子の変位三成分のうち二つが偶、一つが奇になる。剛体並進と回転を除くときは、並進は極性、回転は軸性ベクトルとして変換させる。

:::example[平面反射を持つ非直線分子（1）]
原子数を $m+2$ とし、面を保つ鏡映の偶奇で内部振動を分ける。

$$
m=1
$$

全変位の二つの重複度を数える。

$$
\begin{aligned}
(n_+,n_-)_{\rm all}\big|_{m=1} &= \left( 6, \  3\right)
\end{aligned}
$$

並進は面内二方向が偶で、回転は面に直交する軸一方向が偶である。

$$
\begin{aligned}
(n_+,n_-)_{\rm rigid} &= \left( 3, \  3\right)
\end{aligned}
$$

合計六自由度を除いた値を得る。

$$
\begin{aligned}
(n_+,n_-)_{\rm vib}\big|_{m=1} &= \left( 3, \  0\right)
\end{aligned}
$$
:::

:::example[平面反射を持つ非直線分子（2）]
原子数を $m+2$ とし、面を保つ鏡映の偶奇で内部振動を分ける。

$$
m=2
$$

全変位の二つの重複度を数える。

$$
\begin{aligned}
(n_+,n_-)_{\rm all}\big|_{m=2} &= \left( 8, \  4\right)
\end{aligned}
$$

並進は面内二方向が偶で、回転は面に直交する軸一方向が偶である。

$$
\begin{aligned}
(n_+,n_-)_{\rm rigid} &= \left( 3, \  3\right)
\end{aligned}
$$

合計六自由度を除いた値を得る。

$$
\begin{aligned}
(n_+,n_-)_{\rm vib}\big|_{m=2} &= \left( 5, \  1\right)
\end{aligned}
$$
:::

:::example[平面反射を持つ非直線分子（3）]
原子数を $m+2$ とし、面を保つ鏡映の偶奇で内部振動を分ける。

$$
m=3
$$

全変位の二つの重複度を数える。

$$
\begin{aligned}
(n_+,n_-)_{\rm all}\big|_{m=3} &= \left( 10, \  5\right)
\end{aligned}
$$

並進は面内二方向が偶で、回転は面に直交する軸一方向が偶である。

$$
\begin{aligned}
(n_+,n_-)_{\rm rigid} &= \left( 3, \  3\right)
\end{aligned}
$$

合計六自由度を除いた値を得る。

$$
\begin{aligned}
(n_+,n_-)_{\rm vib}\big|_{m=3} &= \left( 7, \  2\right)
\end{aligned}
$$
:::

:::example[平面反射を持つ非直線分子（4）]
原子数を $m+2$ とし、面を保つ鏡映の偶奇で内部振動を分ける。

$$
m=4
$$

全変位の二つの重複度を数える。

$$
\begin{aligned}
(n_+,n_-)_{\rm all}\big|_{m=4} &= \left( 12, \  6\right)
\end{aligned}
$$

並進は面内二方向が偶で、回転は面に直交する軸一方向が偶である。

$$
\begin{aligned}
(n_+,n_-)_{\rm rigid} &= \left( 3, \  3\right)
\end{aligned}
$$

合計六自由度を除いた値を得る。

$$
\begin{aligned}
(n_+,n_-)_{\rm vib}\big|_{m=4} &= \left( 9, \  3\right)
\end{aligned}
$$
:::

:::example[平面の正三角形分子の内部振動]
三頂点だけに原子があり、点群のうち $C_{3v}$ の操作だけを使う。全変位の指標は $(9,0,1)$、並進と回転の指標はそれぞれ $(3,0,1)$ と $(3,0,-1)$ である。

共役類の大きさと対象の指標を同じ順で並べる。

$$
\begin{aligned}
(n_C)_C &= \left( 1, \  2, \  3\right) \\
(\chi(C))_C &= \left( 3, \  0, \  1\right)
\end{aligned}
$$

A_1 について、各共役類からの寄与を並べて足す。

$$
\begin{aligned}
\chi_{A_1} &= \left( 1, \  1, \  1\right) \\
(n_C\chi(C)\overline{\chi_{A_1}(C)})_C &= \left( 3, \  0, \  3\right) \\
n_{A_1} &= 1
\end{aligned}
$$

A_2 について、各共役類からの寄与を並べて足す。

$$
\begin{aligned}
\chi_{A_2} &= \left( 1, \  1, \  -1\right) \\
(n_C\chi(C)\overline{\chi_{A_2}(C)})_C &= \left( 3, \  0, \  -3\right) \\
n_{A_2} &= 0
\end{aligned}
$$

E について、各共役類からの寄与を並べて足す。

$$
\begin{aligned}
\chi_{E} &= \left( 2, \  -1, \  0\right) \\
(n_C\chi(C)\overline{\chi_{E}(C)})_C &= \left( 6, \  0, \  0\right) \\
n_{E} &= 1
\end{aligned}
$$

各重複度に既約表現の次元を掛け、対象全体の次元を回収する。

$$
\begin{aligned}
\sum_\alpha n_\alpha d_\alpha &= 3 \\
|G| &= 6
\end{aligned}
$$

全点群 $D_{3h}$ では、$C_{3v}$ に水平鏡映 $\sigma_h$ を加える。この鏡映は $C_{3v}$ の全操作と可換である。平面内三原子の全変位は水平鏡映に対して偶が $6$、奇が $3$ である。並進は偶 $2$、奇 $1$、回転は軸性のため偶 $1$、奇 $2$ となる。従って内部振動はすべて偶であり、$D_{3h}$ では $A_1'\oplus E'$ になる。
:::

## 対称性による選択則と行列要素

### 行列要素を群平均する

作用素成分の空間が表現 $\Gamma_O$ に従うとする。$\langle v_f,O v_i\rangle_{\rm qm}$ は Dirac の規約で第 2 引数に線形なので、初期状態・作用素・終状態について $\Gamma_f^*\otimes\Gamma_O\otimes\Gamma_i$ に対する不変な線形汎関数を与える。群平均 $\Pi=|G|^{-1}\sum_gD(g)$ は不変部分空間への射影である。実際どの元を掛けても和が不変で、不変なベクトル上では恒等である。従って自明表現が含まれなければ行列要素は零になる。含まれても動力学的に相殺される可能性があり、対称性だけから値が非零とは言えない。

例えば $C_{3v}$ の電気双極子は極性ベクトルなので $z$ が $A_1$、$(x,y)$ が $E$ である。$A_1$ から $A_1$ への遷移は $z$ 偏光で許され、平面偏光では零である。$A_1$ から $E$ では逆になる。$E\otimes E$ の指標は $(4,1,0)$ であり、表との内積から $A_1+A_2+E$ となる。この場合、単に二重縮退と知るだけでは偏光の条件を区別できない。

調和振動の基本遷移では、平衡点での双極子の一次項 $(\partial\boldsymbol\mu/\partial Q_i)_0Q_i$ を調べる。定数項は振動の基底状態と第一励起状態の直交性により寄与せず、$\mel{1}{Q_i}{0}=\sqrt{\hbar/(2\omega_i)}$ である。これは質量重み付き正規座標を生成・消滅演算子の和に書けば従う。従って赤外吸収で必要なのは、そのモード対称型と双極子成分とが同型で、その成分の微分が零でないことである。Raman 散乱の分極率にも同じ一次展開を使う。反転中心があれば、赤外では奇のモード、Raman では偶のモードを調べるため、一つの基本モードは両方で同時には許容されない。

### Raman の選択則では二つの添字を使う

線形な分極応答の分極率は電場と誘起双極子を結ぶ二階テンソルである。非磁性で対称な分極率を仮定すると、その成分の表現は極性ベクトルの対称二乗である。$C_{3v}$ では $z^2,x^2+y^2$ が二つの $A_1$、$(xz,yz)$ と $(x^2-y^2,2xy)$ が二つの $E$ を与える。後ろの組は平面の回転角に対して二倍の角で動くが、$120^\circ$ では指標は $2\cos240^\circ=-1$、鏡映では零であるから同じ型 $E$ である。次元は $2+2+2=6$ となる。従って自明な基底からの Raman 遷移では $A_1,E$ が対称性で許され、$A_2$ はこの対称テンソルでは許されない。例えば $zz$ の偏光組合せは $z^2$ と同じ $A_1$ なので、$E$ の成分を励起するには別の組合せが要る。

### 巡回群での選択則

$C_n$ の一つの生成元に対する状態の固有値を $e^{2\pi i j/n}$ とする。演算子が共役で $e^{2\pi i q/n}$ 倍されるなら、行列要素には $e^{2\pi i(j_i+q-j_f)/n}$ が掛かる。群平均の指数和はこの位相が $1$ の場合だけ残る。従って $j_f-j_i=q$ を $n$ を法として満たす必要がある。この条件は対称性からの必要条件なので、実際の行列要素が非零かは動力学と状態に依存する。

:::example[四回軸での円偏光の符号（1）]
初期添字を $m$、演算子成分を添字 $+1,-1$ とし、終状態を調べる。

$$
m=1
$$

添字を $0,1,2,3$ に戻す。

$$
\begin{aligned}
j_i\big|_{m=1} &= 1
\end{aligned}
$$

正の円偏光は添字を一つ上げる。

$$
\begin{aligned}
j_f^{(+)}\big|_{m=1} &= 2
\end{aligned}
$$

負の円偏光は添字を一つ下げる。

$$
\begin{aligned}
j_f^{(-)}\big|_{m=1} &= 0
\end{aligned}
$$
:::

:::example[四回軸での円偏光の符号（2）]
初期添字を $m$、演算子成分を添字 $+1,-1$ とし、終状態を調べる。

$$
m=2
$$

添字を $0,1,2,3$ に戻す。

$$
\begin{aligned}
j_i\big|_{m=2} &= 2
\end{aligned}
$$

正の円偏光は添字を一つ上げる。

$$
\begin{aligned}
j_f^{(+)}\big|_{m=2} &= 3
\end{aligned}
$$

負の円偏光は添字を一つ下げる。

$$
\begin{aligned}
j_f^{(-)}\big|_{m=2} &= 1
\end{aligned}
$$
:::

:::example[四回軸での円偏光の符号（3）]
初期添字を $m$、演算子成分を添字 $+1,-1$ とし、終状態を調べる。

$$
m=3
$$

添字を $0,1,2,3$ に戻す。

$$
\begin{aligned}
j_i\big|_{m=3} &= 3
\end{aligned}
$$

正の円偏光は添字を一つ上げる。

$$
\begin{aligned}
j_f^{(+)}\big|_{m=3} &= 0
\end{aligned}
$$

負の円偏光は添字を一つ下げる。

$$
\begin{aligned}
j_f^{(-)}\big|_{m=3} &= 2
\end{aligned}
$$
:::

:::example[四回軸での円偏光の符号（4）]
初期添字を $m$、演算子成分を添字 $+1,-1$ とし、終状態を調べる。

$$
m=4
$$

添字を $0,1,2,3$ に戻す。

$$
\begin{aligned}
j_i\big|_{m=4} &= 0
\end{aligned}
$$

正の円偏光は添字を一つ上げる。

$$
\begin{aligned}
j_f^{(+)}\big|_{m=4} &= 1
\end{aligned}
$$

負の円偏光は添字を一つ下げる。

$$
\begin{aligned}
j_f^{(-)}\big|_{m=4} &= 3
\end{aligned}
$$
:::

:::example[二重項間での偏光の選択]
初期と終状態が $E$ であるとき、垂直方向の双極子は $A_1$ であり、面内の双極子は $E$ である。まず $E^*\otimes E$ を分解して演算子の対称型を判定する。

共役類の大きさと対象の指標を同じ順で並べる。

$$
\begin{aligned}
(n_C)_C &= \left( 1, \  2, \  3\right) \\
(\chi(C))_C &= \left( 4, \  1, \  0\right)
\end{aligned}
$$

A_1 について、各共役類からの寄与を並べて足す。

$$
\begin{aligned}
\chi_{A_1} &= \left( 1, \  1, \  1\right) \\
(n_C\chi(C)\overline{\chi_{A_1}(C)})_C &= \left( 4, \  2, \  0\right) \\
n_{A_1} &= 1
\end{aligned}
$$

A_2 について、各共役類からの寄与を並べて足す。

$$
\begin{aligned}
\chi_{A_2} &= \left( 1, \  1, \  -1\right) \\
(n_C\chi(C)\overline{\chi_{A_2}(C)})_C &= \left( 4, \  2, \  0\right) \\
n_{A_2} &= 1
\end{aligned}
$$

E について、各共役類からの寄与を並べて足す。

$$
\begin{aligned}
\chi_{E} &= \left( 2, \  -1, \  0\right) \\
(n_C\chi(C)\overline{\chi_{E}(C)})_C &= \left( 8, \  -2, \  0\right) \\
n_{E} &= 1
\end{aligned}
$$

各重複度に既約表現の次元を掛け、対象全体の次元を回収する。

$$
\begin{aligned}
\sum_\alpha n_\alpha d_\alpha &= 4 \\
|G| &= 6
\end{aligned}
$$

$A_1$ と $E$ のいずれも含まれるので、両偏光の遷移が対称性だけで排除されるわけではない。個々の成分には別に行列要素の計算が要る。
:::

### 反転と多重極の偶奇

中心反転を含む群では、一次元の符号 $\pi_g=\pm1$ によって偶奇を表す。初期準位の偶奇を $\pi_i$、終準位を $\pi_f$、演算子の偶奇を $\pi_T$ とすると
$$
\langle f|T|i\rangle=\pi_f\pi_T\pi_i\langle f|T|i\rangle
$$
となる。これは反転のユニタリ性を使って左右に同じ反転を挿入した結果である。従って電気双極子 $\vb r$ は奇であり、初状態と終状態の偶奇が異なる必要がある。磁気双極子の $\vb L=\vb r\times\vb p$ は二つの極性ベクトルの積なので偶であり、同じ偶奇を要求する。電気四重極も座標の二次で偶である。同じ偶奇条件を満たすからといって磁気双極子と四重極が同じ全選択則を持つわけではない。回転でのテンソルの階数が違うためである。
:::example[軌道の偶奇による行列要素の制限]
偶準位二つと奇準位一つを取り、奇な Hermite 演算子を作る。偶準位どうしのブロックと奇準位どうしの成分は零である。

座標成分を列に並べる。この順序に合わせ、線形写像の行列の列を選ぶ。

$$
\begin{aligned}v&=\begin{pmatrix}1\\2\\0\end{pmatrix}\\ A&=\begin{pmatrix}0 & 0 & 2\\0 & 0 & 1\\2 & 1 & 0\end{pmatrix}\end{aligned}
$$

行列の各行と入力の内積を計算する。例えば第一行では $\begin{pmatrix}0 & 0 & 2\end{pmatrix}$ と入力列を掛ける。

$$
\begin{aligned}(Av)_1&=0\\Av&=\begin{pmatrix}0\\0\\4\end{pmatrix}\end{aligned}
$$

核の次元は列数から階数を引いて得られる。像の次元は階数である。

$$
\begin{aligned}\operatorname{rank}A&=2\\\dim\ker A&=1\end{aligned}
$$

正方行列の場合は、特性多項式を求めて固有値を調べられる。ここでは $tI-A$ の行列式を取る。

$$
\begin{aligned}\operatorname{tr}A&=0\\ \det A&=0\\ \det(tI-A)&=t^{3} - 5 t\end{aligned}
$$

行列要素を得るため、出力から指定した共変ベクトルで成分を取り出す。

$$
\begin{aligned}\ell&=\begin{pmatrix}0 & 0 & 1\end{pmatrix}\\ \ell Av&=4\end{aligned}
$$

次に別の方向を入力にする。出力の差は入力の差に $A$ を掛けたものである。

$$
\begin{aligned}w&=\begin{pmatrix}0\\0\\1\end{pmatrix}\\ Aw&=\begin{pmatrix}2\\1\\0\end{pmatrix}\\ Aw-Av&=\begin{pmatrix}2\\1\\-4\end{pmatrix}\end{aligned}
$$
:::

## 結晶場分裂と縮退

### 五つの軌道が二組に分かれる

立方体の向きを保つ対称操作を含む八面体環境では、二次調和多項式の空間は対角型 $2z^2-x^2-y^2,x^2-y^2$ と非対角型 $xy,yz,zx$ に分かれる。符号付きの座標置換はそれぞれの空間を保つ。前者は三座標の平方の和が零の係数平面で、三巡回の二つの複素固有線が互換で交換されるため既約である。後者では二つの座標を反転する三種類の対称操作で同時固有線を分けられ、座標置換が三本を推移的に移す。従ってどの真の非零不変空間も存在しない。反転では二次式は偶なので、これらは $E_g$ と $T_{2g}$ である。

結晶場の作用素が群と交換すると、既約部分ごとにスカラーであることは指標の証明に使った可換写像の性質から従う。従って一般に二重と三重へ分裂する。重心を零、間隔を $\Delta=E(E_g)-E(T_{2g})$ と置くと $2E(E_g)+3E(T_{2g})=0$ なのでエネルギーは $3\Delta/5,-2\Delta/5$ である。例えば $\Delta=1\,\mathrm{eV}$ なら $0.6\,\mathrm{eV}$ と $-0.4\,\mathrm{eV}$ になる。どちらが高いかは環境の作用素が決め、群論はこの間隔の符号も大きさも固定しない。

:::example[正方対称への低下で残る縮退]
八面体を $z$ 方向だけに伸縮させ、$D_{4h}$ 対称性を保つとする。$2z^2-x^2-y^2$ は全操作で不変なので $A_{1g}$、$x^2-y^2$ は四分の一回転で奇、座標面鏡映で偶なので $B_{1g}$ である。従って $O_h$ の $E_g$ の二重縮退は一般に解ける。$xy$ は四分の一回転と座標面鏡映の両方で奇なので $B_{2g}$ であり、$(xz,yz)$ は四分の一回転で互いに移り、鏡映で回転の二つの固有線が交換されるため $D_{4h}$ の二次元既約型 $E_g$ になる。従って $O_h$ の $T_{2g}$ は一重項と二重項へ分かれ、残る二重縮退は新しい対称性が保つ。
:::

### 同じ既約型が二つ以上ある場合

ある既約型 $V_\alpha$ が $m$ 回現れる空間を $\C^m\otimes V_\alpha$ と書く。群と交換する作用素を各コピー間のブロックに分けると、各ブロックは既約間の可換写像なのでスカラー倍の恒等である。そのため $H$ の制限は $A\otimes I_{d_\alpha}$ になる。もし $H$ が自己共役なら $A$ は Hermitian である。$A$ を対角化すると、各固有値に少なくとも $d_\alpha$ 個の状態がある。例えば $C_{3v}$ で $E$ が二回現れ、$A=\begin{pmatrix}2&1\\1&2\end{pmatrix}$ なら固有値 $1,3$ がそれぞれ二重になる。$E$ が二回あるだけでは四重縮退は保証されない。

### 同じ対称型の混成

群と可換な Hamiltonian を既約成分に分ける。互いに同型でない成分間のブロックは Schur の補題で零になる。同型の成分が二つ以上あれば、その重複空間では行列を対角化する必要がある。たとえば二次元既約表現を二つ持つなら、Hamiltonian は $I_2\otimes A$ と書ける。$A$ の固有値が異なっても、それぞれには二重縮退が残る。

:::example[同じ二重項を混ぜる行列（1）]
二つの二重項の非混成準位を $m$ と $-m$、混成を $1$ とする。

$$
m=1
$$

対称型の内部では恒等行列である。重複空間の行列を書く。

$$
\begin{aligned}
A\big|_{m=1} &= \begin{pmatrix}1 & 1\\1 & -1\end{pmatrix}
\end{aligned}
$$

二乗を取ると対角成分だけが残る。

$$
\begin{aligned}
A^2\big|_{m=1} &= \begin{pmatrix}2 & 0\\0 & 2\end{pmatrix}
\end{aligned}
$$

二つの固有値に同じ次元の対称型が付く。

$$
\begin{aligned}
(E_-,E_+)\big|_{m=1} &= \left( - \sqrt{2}, \  \sqrt{2}\right)
\end{aligned}
$$
:::

:::example[同じ二重項を混ぜる行列（2）]
二つの二重項の非混成準位を $m$ と $-m$、混成を $1$ とする。

$$
m=2
$$

対称型の内部では恒等行列である。重複空間の行列を書く。

$$
\begin{aligned}
A\big|_{m=2} &= \begin{pmatrix}2 & 1\\1 & -2\end{pmatrix}
\end{aligned}
$$

二乗を取ると対角成分だけが残る。

$$
\begin{aligned}
A^2\big|_{m=2} &= \begin{pmatrix}5 & 0\\0 & 5\end{pmatrix}
\end{aligned}
$$

二つの固有値に同じ次元の対称型が付く。

$$
\begin{aligned}
(E_-,E_+)\big|_{m=2} &= \left( - \sqrt{5}, \  \sqrt{5}\right)
\end{aligned}
$$
:::

:::example[同じ二重項を混ぜる行列（3）]
二つの二重項の非混成準位を $m$ と $-m$、混成を $1$ とする。

$$
m=3
$$

対称型の内部では恒等行列である。重複空間の行列を書く。

$$
\begin{aligned}
A\big|_{m=3} &= \begin{pmatrix}3 & 1\\1 & -3\end{pmatrix}
\end{aligned}
$$

二乗を取ると対角成分だけが残る。

$$
\begin{aligned}
A^2\big|_{m=3} &= \begin{pmatrix}10 & 0\\0 & 10\end{pmatrix}
\end{aligned}
$$

二つの固有値に同じ次元の対称型が付く。

$$
\begin{aligned}
(E_-,E_+)\big|_{m=3} &= \left( - \sqrt{10}, \  \sqrt{10}\right)
\end{aligned}
$$
:::

:::example[同じ二重項を混ぜる行列（4）]
二つの二重項の非混成準位を $m$ と $-m$、混成を $1$ とする。

$$
m=4
$$

対称型の内部では恒等行列である。重複空間の行列を書く。

$$
\begin{aligned}
A\big|_{m=4} &= \begin{pmatrix}4 & 1\\1 & -4\end{pmatrix}
\end{aligned}
$$

二乗を取ると対角成分だけが残る。

$$
\begin{aligned}
A^2\big|_{m=4} &= \begin{pmatrix}17 & 0\\0 & 17\end{pmatrix}
\end{aligned}
$$

二つの固有値に同じ次元の対称型が付く。

$$
\begin{aligned}
(E_-,E_+)\big|_{m=4} &= \left( - \sqrt{17}, \  \sqrt{17}\right)
\end{aligned}
$$
:::

:::example[立方体の二次調和多項式]
固有回転の共役類を $(I,8C_3,3C_2,6C_4,6C_2^{\prime})$ の順に取る。二次調和多項式は空間ベクトルの対称二乗からトレース一成分を除いたものである。二類の半回転は空間ベクトル表現では指標がともに $-1$ である。

共役類の大きさと対象の指標を同じ順で並べる。

$$
\begin{aligned}
(n_C)_C &= \left( 1, \  8, \  3, \  6, \  6\right) \\
(\chi(C))_C &= \left( 5, \  -1, \  1, \  -1, \  1\right)
\end{aligned}
$$

A_1 について、各共役類からの寄与を並べて足す。

$$
\begin{aligned}
\chi_{A_1} &= \left( 1, \  1, \  1, \  1, \  1\right) \\
(n_C\chi(C)\overline{\chi_{A_1}(C)})_C &= \left( 5, \  -8, \  3, \  -6, \  6\right) \\
n_{A_1} &= 0
\end{aligned}
$$

A_2 について、各共役類からの寄与を並べて足す。

$$
\begin{aligned}
\chi_{A_2} &= \left( 1, \  1, \  1, \  -1, \  -1\right) \\
(n_C\chi(C)\overline{\chi_{A_2}(C)})_C &= \left( 5, \  -8, \  3, \  6, \  -6\right) \\
n_{A_2} &= 0
\end{aligned}
$$

E について、各共役類からの寄与を並べて足す。

$$
\begin{aligned}
\chi_{E} &= \left( 2, \  -1, \  2, \  0, \  0\right) \\
(n_C\chi(C)\overline{\chi_{E}(C)})_C &= \left( 10, \  8, \  6, \  0, \  0\right) \\
n_{E} &= 1
\end{aligned}
$$

T_1 について、各共役類からの寄与を並べて足す。

$$
\begin{aligned}
\chi_{T_1} &= \left( 3, \  0, \  -1, \  1, \  -1\right) \\
(n_C\chi(C)\overline{\chi_{T_1}(C)})_C &= \left( 15, \  0, \  -3, \  -6, \  -6\right) \\
n_{T_1} &= 0
\end{aligned}
$$

T_2 について、各共役類からの寄与を並べて足す。

$$
\begin{aligned}
\chi_{T_2} &= \left( 3, \  0, \  -1, \  -1, \  1\right) \\
(n_C\chi(C)\overline{\chi_{T_2}(C)})_C &= \left( 15, \  0, \  -3, \  6, \  6\right) \\
n_{T_2} &= 1
\end{aligned}
$$

各重複度に既約表現の次元を掛け、対象全体の次元を回収する。

$$
\begin{aligned}
\sum_\alpha n_\alpha d_\alpha &= 5 \\
|G| &= 24
\end{aligned}
$$

二重項には $2z^2-x^2-y^2,x^2-y^2$、三重項には $xy,yz,zx$ を取れる。中心反転を含める場合はすべて偶になる。
:::

## 時間反転対称性と Kramers 縮退

### 複素共役を伴う対称性

反ユニタリ $\Theta$ は反線形で $\inner{\Theta u}{\Theta v}=\overline{\inner{u}{v}}$ を満たす。ここでは第 1 引数を線形とする。時間反転は運動量の符号を変えるため、量子力学の交換関係の $i$ も共役する必要がある。複素共役 $K$ と Pauli 行列 $\sigma_y$ を用いたスピン $1/2$ の $\Theta=-i\sigma_y K$ では $-i\sigma_y$ は実行列で、その二乗が $-I$ である。

:::theorem[Kramers 縮退]{#thm-mp-kramers}
自己共役 $H$ が定義域を保つ反ユニタリ $\Theta$ と交換し、$\Theta^2=-I$ なら、有限次元の各エネルギー固有空間は偶数次元である。
:::

:::proof
固有値は実なので $Hv=Ev$ から $H\Theta v=E\Theta v$ となる。$\inner{\Theta v}{\Theta^2v}=\overline{\inner{v}{\Theta v}}$ である一方、左辺は $-\inner{\Theta v}{v}=-\overline{\inner{v}{\Theta v}}$ である。従って $v$ と $\Theta v$ は直交する。この二次元部分空間は $\Theta$ で保たれ、その直交補も反ユニタリ性で保たれる。有限次元なら直交補で同じ組を取り続けられ、次元は偶数になる。
:::

有限次元としなくても、任意の非零固有ベクトルに直交する相方が存在することまでは同じ証明で従う。例えば $\Theta=K$ だけで $\Theta^2=+I$ なら、実対角 $H=\diag(1,2)$ は交換しても固有値はどちらも単純である。

また、自己共役な摂動 $V$ も $\Theta$ と交換するなら、Kramers 対における対角要素は等しく、非対角要素は零である。第一引数を線形とする内積で、$\inner{V\Theta v}{\Theta v}=\overline{\inner{Vv}{v}}=\inner{Vv}{v}$ である。また $\inner{\Theta v}{\Theta(V\Theta v)}=-\inner{\Theta v}{Vv}=-\overline{\inner{Vv}{\Theta v}}$ である一方、反ユニタリ性から左辺は $\overline{\inner{v}{V\Theta v}}=\overline{\inner{Vv}{\Theta v}}$ である。従って一次の摂動行列は一つの対の上でスカラーになる。

### 時間反転の二乗と奇な項

スピン $1/2$ で $\Theta=-i\sigma_yK$ とすると、$U=-i\sigma_y$ は実行列なので $\Theta^2=U\overline U=U^2=-I$ である。直接掛け算して $\Theta\sigma_j\Theta^{-1}=-\sigma_j$ を得る。二次元だけに限れば任意の Hermite 行列は $aI+b_j\sigma_j$ と実係数で書けるので、時間反転不変なら $b_j=0$ である。従って単一の二重項では、Pauli 行列に比例する磁場の項を加えれば時間反転の仮定が失われる。

:::example[スピン二重項に磁場を加える（1）]
無次元の磁場を $m$ として $H=m\sigma_z$ の準位を取る。

$$
m=1
$$

磁場を固定して Hamiltonian を書く。

$$
\begin{aligned}
H\big|_{m=1} &= \begin{pmatrix}1 & 0\\0 & -1\end{pmatrix}
\end{aligned}
$$

時間反転で磁気モーメントの符号が反転する。

$$
\begin{aligned}
\Theta H\Theta^{-1}\big|_{m=1} &= \begin{pmatrix}-1 & 0\\0 & 1\end{pmatrix}
\end{aligned}
$$

固有値の差は磁場の強さに比例する。

$$
\begin{aligned}
E_+-E_-\big|_{m=1} &= 2
\end{aligned}
$$
:::

:::example[スピン二重項に磁場を加える（2）]
無次元の磁場を $m$ として $H=m\sigma_z$ の準位を取る。

$$
m=2
$$

磁場を固定して Hamiltonian を書く。

$$
\begin{aligned}
H\big|_{m=2} &= \begin{pmatrix}2 & 0\\0 & -2\end{pmatrix}
\end{aligned}
$$

時間反転で磁気モーメントの符号が反転する。

$$
\begin{aligned}
\Theta H\Theta^{-1}\big|_{m=2} &= \begin{pmatrix}-2 & 0\\0 & 2\end{pmatrix}
\end{aligned}
$$

固有値の差は磁場の強さに比例する。

$$
\begin{aligned}
E_+-E_-\big|_{m=2} &= 4
\end{aligned}
$$
:::

:::example[スピン二重項に磁場を加える（3）]
無次元の磁場を $m$ として $H=m\sigma_z$ の準位を取る。

$$
m=3
$$

磁場を固定して Hamiltonian を書く。

$$
\begin{aligned}
H\big|_{m=3} &= \begin{pmatrix}3 & 0\\0 & -3\end{pmatrix}
\end{aligned}
$$

時間反転で磁気モーメントの符号が反転する。

$$
\begin{aligned}
\Theta H\Theta^{-1}\big|_{m=3} &= \begin{pmatrix}-3 & 0\\0 & 3\end{pmatrix}
\end{aligned}
$$

固有値の差は磁場の強さに比例する。

$$
\begin{aligned}
E_+-E_-\big|_{m=3} &= 6
\end{aligned}
$$
:::

:::example[スピン二重項に磁場を加える（4）]
無次元の磁場を $m$ として $H=m\sigma_z$ の準位を取る。

$$
m=4
$$

磁場を固定して Hamiltonian を書く。

$$
\begin{aligned}
H\big|_{m=4} &= \begin{pmatrix}4 & 0\\0 & -4\end{pmatrix}
\end{aligned}
$$

時間反転で磁気モーメントの符号が反転する。

$$
\begin{aligned}
\Theta H\Theta^{-1}\big|_{m=4} &= \begin{pmatrix}-4 & 0\\0 & 4\end{pmatrix}
\end{aligned}
$$

固有値の差は磁場の強さに比例する。

$$
\begin{aligned}
E_+-E_-\big|_{m=4} &= 8
\end{aligned}
$$
:::

### 奇数個のスピンと時間反転

複数の独立したスピン空間を直積にし、共通の積基底で複素共役を取る。時間反転を $\Theta_N=U^{\otimes N}K$、$U=-i\sigma_y$ とすれば、この基底で $U$ は実なので $\Theta_N^2=(U^2)^{\otimes N}=(-1)^N I$ となる。交換対称性などにより物理状態をこの中の部分空間に制限する場合も、その部分空間を $\Theta_N$ が保てば同じ二乗が作用する。

$N=2$ の場合について、$\Theta_2\ket{\uparrow\uparrow}=\ket{\downarrow\downarrow}$、$\Theta_2\ket{\downarrow\downarrow}=\ket{\uparrow\uparrow}$ である。さらに $\ket{\uparrow\downarrow}$ と $\ket{\downarrow\uparrow}$ は互いの負号付きの像になる。従って反対称の一重項は $\Theta_2$ で固定される。例えば $H=J\,\vb S_1\cdot\vb S_2/\hbar^2$、$J>0$ の一重項は $-3J/4$、三重項は $J/4$ であり、一重項の非縮退は時間反転と両立する。エネルギー $J$ の単位を $\mathrm{eV}$ に取るなら、$J=0.20\,\mathrm{eV}$ のとき一重項は $-0.15\,\mathrm{eV}$、三重項は $0.05\,\mathrm{eV}$ である。

この模型で $\vb S_1\cdot\vb S_2$ が不変なのは二本のスピンが同時に符号を変えるためである。一方、外から固定した磁場に対する $S_{1z}+S_{2z}$ は符号を変える。磁場を反転した別の Hamiltonian との間に関係があっても、同じ磁場での縮退は時間反転だけでは保証されない。


:::example[二粒子の時間反転の実行列部分]
基底を上上、上下、下上、下下の順に取る。一重項の係数はここでは規格化せず $(0,1,-1,0)$ とする。

座標成分を列に並べる。この順序に合わせ、線形写像の行列の列を選ぶ。

$$
\begin{aligned}v&=\begin{pmatrix}0\\1\\-1\\0\end{pmatrix}\\ A&=\begin{pmatrix}0 & 0 & 0 & 1\\0 & 0 & -1 & 0\\0 & -1 & 0 & 0\\1 & 0 & 0 & 0\end{pmatrix}\end{aligned}
$$

行列の各行と入力の内積を計算する。例えば第一行では $\begin{pmatrix}0 & 0 & 0 & 1\end{pmatrix}$ と入力列を掛ける。

$$
\begin{aligned}(Av)_1&=0\\Av&=\begin{pmatrix}0\\1\\-1\\0\end{pmatrix}\end{aligned}
$$

核の次元は列数から階数を引いて得られる。像の次元は階数である。

$$
\begin{aligned}\operatorname{rank}A&=4\\\dim\ker A&=0\end{aligned}
$$

正方行列の場合は、特性多項式を求めて固有値を調べられる。ここでは $tI-A$ の行列式を取る。

$$
\begin{aligned}\operatorname{tr}A&=0\\ \det A&=1\\ \det(tI-A)&=t^{4} - 2 t^{2} + 1\end{aligned}
$$

行列要素を得るため、出力から指定した共変ベクトルで成分を取り出す。

$$
\begin{aligned}\ell&=\begin{pmatrix}0 & 1 & -1 & 0\end{pmatrix}\\ \ell Av&=2\end{aligned}
$$

次に別の方向を入力にする。出力の差は入力の差に $A$ を掛けたものである。

$$
\begin{aligned}w&=\begin{pmatrix}1\\0\\0\\0\end{pmatrix}\\ Aw&=\begin{pmatrix}0\\0\\0\\1\end{pmatrix}\\ Aw-Av&=\begin{pmatrix}0\\-1\\1\\1\end{pmatrix}\end{aligned}
$$
:::

## 演習問題

::::exercise{level=basic}
正三角形の回転 $r$ と鏡映 $s$ について $(rs)^2$ を求めよ。
:::solution
$s r=r^{-1}s$ だから $(rs)^2=r(sr)s=r r^{-1}s^2=e$ である。
:::
::::

::::exercise{level=basic}
三頂点の表現の $E$ の重複度を計算せよ。
:::solution
指標 $(3,0,1)$ と $(2,-1,0)$ の内積は $(3\cdot2+2\cdot0+3\cdot0)/6=1$ である。
:::
::::

::::exercise{level=basic}
$60^\circ$ の回転行列の三次元トレースと位数を答えよ。
:::solution
$1+2\cos60^\circ=2$、六回で一周するので位数は $6$ である。
:::
::::

::::exercise{level=basic}
$C_{3v}$ の極性ベクトルの鏡映の指標を求めよ。
:::solution
鏡面内の二方向は不変、面に垂直な方向は反転するので $1+1-1=1$ である。
:::
::::

::::exercise{level=advanced}
$A_2\otimes E$ を分解せよ。
:::solution
指標は点ごとの積で $(1,1,-1)\cdot(2,-1,0)=(2,-1,0)$ となる。表を使えば $E$ が一つ、それ以外は零である。
:::
::::

::::exercise{level=advanced}
$A_1$ から $A_2$ への電気双極子遷移を $C_{3v}$ で判定せよ。
:::solution
双極子表現は $A_1+E$ であり、全積は $A_2\otimes(A_1+E)\otimes A_1=A_2+E$ である。自明表現が無いので全偏光で零である。
:::
::::

::::exercise{level=advanced}
八面体の $E_g$ と $T_{2g}$ の重心差の間隔が $2\,\mathrm{eV}$ のとき、重心を原点とする二準位を求めよ。
:::solution
$E_2-E_3=2\,\mathrm{eV}$、$2E_2+3E_3=0$ である。この二元連立を解くと $E_2=1.2\,\mathrm{eV}$、$E_3=-0.8\,\mathrm{eV}$ である。
:::
::::

::::exercise{level=advanced}
$\Theta=-i\sigma_yK$ を $(1,0)^\top$ に二回作用させよ。
:::solution
行列 $-i\sigma_y=\begin{pmatrix}0&-1\\1&0\end{pmatrix}$ は実である。共役は二つの実基底を変えない。一回目は $(0,1)^\top$、二回目は $(-1,0)^\top$ となり、初期状態と一回目が直交する。
:::
::::

## 確認テスト

::::quiz{answer="6" tol="1e-9"}
$C_{3v}$ の位数を答えよ。

:::explanation
回転 $e,r,r^2$ と鏡映 $s,rs,r^2s$ があり、頂点への作用で互いに区別できるので $6$ である。
:::
::::

::::quiz
電気双極子遷移を対称性が禁じる条件はどれか。

- [x] 状態と作用素のテンソル積に自明表現が無い
- [ ] どちらかの状態が一次元表現である
- [ ] 作用素がベクトルの成分である

:::explanation
行列要素は群不変であり、不変な汎関数は自明表現でだけ値を持つ。従って自明成分が無ければ零である。
:::
::::

::::quiz{answer="2" tol="1e-9"}
三回軸を持つ本文の四原子分子の振動表現で $E$ の重複度はいくつか。

:::explanation
振動指標 $(6,0,2)$ と $E$ の $(2,-1,0)$ の重み付き内積は $12/6=2$ である。
:::
::::

::::quiz
Kramers の定理に使う反ユニタリ作用素の二乗はどれか。

- [x] $-I$
- [ ] $+I$
- [ ] 零作用素

:::explanation
$\Theta^2=-I$ の負符号によって $\inner{v}{\Theta v}$ が自分自身の負に等しくなる。そのため二つの同じエネルギーの状態が直交する。
:::
::::

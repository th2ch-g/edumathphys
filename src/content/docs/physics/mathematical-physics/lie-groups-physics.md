---
title: "連続群と物理"
description: "連続対称性を生成子で扱い、回転・内部対称性・Lorentz 変換の表現を物理量に対応させる。"
notation:
  - symbol: 'J_i'
    meaning: "角運動量の自己共役生成子"
  - symbol: 'j,m'
    meaning: "スピンとその z 成分の無次元固有値"
  - symbol: 'C^{JM}_{j_1m_1,j_2m_2}'
    meaning: "Clebsch–Gordan 係数"
  - symbol: '(p,q)'
    meaning: "SU(3) の最高ウェイトの Dynkin ラベル"
---

## この章の目標

- 行列 Lie 群を単位元で微分して生成子を求める。
- SO(3) と SU(2) の被覆関係から整数・半整数スピンを区別する。
- 角運動量を合成し、テンソル演算子の行列要素を絞る。
- SU(3) と Lorentz 群の表現を内部状態と時空の量に対応させる。

[SU(2) と SO(3)](/math/lie-groups/su2-so3/) の関係を出発点に、角運動量と有限次元表現を計算する。量子力学の状態の行列要素には ket に線形な Dirac の規約を用いる。

## Lie 群と Lie 代数の入門

### 単位元の近くで対称操作を調べる

滑らかな多様体でもある群で、積と逆が滑らかなものを :term[リー群]{en="Lie group" yomi="りーぐん"} という。行列 Lie 群の単位元での接空間は :term[リー代数]{en="Lie algebra" yomi="りーだいすう"} であり、行列の括弧 $[X,Y]=XY-YX$ が入る。接ベクトル $X$ は、左移動で定めたベクトル場 $gX$ の流れを単位元から解くと $g(t)=e^{tX}$ を与える。二つの指数の群交換子を二次まで展開すると $e^{tX}e^{tY}e^{-tX}e^{-tY}=I+t^2[X,Y]+O(t^3)$ になり、括弧の意味が分かる。閉性は例えば $e^{tX}Ye^{-tX}$ が各 $t$ で接空間に属し、その導関数が $[X,Y]$ であることから従う。Jacobi 恒等式は三つの二重交換子を展開すると各積が一度ずつ逆符号で現れて消える。

$R^\top R=I$ の微分では反対称行列が、$U^\dagger U=I,\det U=1$ の微分では反 Hermitian かつトレース零の行列が得られる。それぞれ指数を取れば対応する群に入るため、条件は接空間をちょうど特徴付ける。ユニタリ表現で軸の角度 $\theta$ に対して $U=e^{-i\theta J/\hbar}$ と置くと $J$ が自己共役になる。例えばスピン $1/2$ では $J_z=\hbar\diag(1,-1)/2$ であり、回転で二成分の位相が逆向きに進む。

### 微小操作の交換を測る

行列 Lie 群の単位元を通る曲線の微分を $X,Y$ とする。二方向に短く進み、逆順の正負を掛け合わせると $e^{tX}e^{tY}e^{-tX}e^{-tY}=I+t^2[X,Y]+O(t^3)$ である。これは四つの Taylor 展開の一次項が消え、二次項から $XY-YX$ が残るためである。従って Lie 代数の括弧は二つの微小変換の順序の差を表す。

積の対数を取ったときにも交換子が現れる。$e^{tX}e^{tY}$ の二次項は $X^2/2+XY+Y^2/2$ である。一方、指数の中身を $t(X+Y)+t^2C$ とすると、二次項は $(X+Y)^2/2+C$ になる。比較すると $C=[X,Y]/2$ であり、これは :ref[thm-lie-bch] の最初の補正項である。また滑らかな群表現 $D$ の単位元での微分を $\dd D$ とすれば、$D(e^{tX})=e^{t\,\dd D(X)}$ である。群交換子の二次項を $D$ の前後で比較すると $\dd D([X,Y])=[\dd D(X),\dd D(Y)]$ が従う。従って生成子を行列で定めるときは、まずこの交換関係を確認する。

:::example[Pauli 基底の交換子（1）]
反 Hermite な生成元を二つ選び、第二次数に残る項を調べる。

$$
m=1
$$

一方の生成元の強さを変える。

$$
\begin{aligned}
X\big|_{m=1} &= \begin{pmatrix}0 & - \frac{i}{2}\\- \frac{i}{2} & 0\end{pmatrix}
\end{aligned}
$$

他方は固定する。

$$
\begin{aligned}
Y &= \begin{pmatrix}0 & - \frac{1}{2}\\\frac{1}{2} & 0\end{pmatrix}
\end{aligned}
$$

行列の積の差を取る。

$$
\begin{aligned}
[X,Y]\big|_{m=1} &= \begin{pmatrix}- \frac{i}{2} & 0\\0 & \frac{i}{2}\end{pmatrix}
\end{aligned}
$$

trace は零になる。

$$
\begin{aligned}
\operatorname{tr}[X,Y] &= 0
\end{aligned}
$$
:::

:::example[Pauli 基底の交換子（2）]
反 Hermite な生成元を二つ選び、第二次数に残る項を調べる。

$$
m=2
$$

一方の生成元の強さを変える。

$$
\begin{aligned}
X\big|_{m=2} &= \begin{pmatrix}0 & - i\\- i & 0\end{pmatrix}
\end{aligned}
$$

他方は固定する。

$$
\begin{aligned}
Y &= \begin{pmatrix}0 & - \frac{1}{2}\\\frac{1}{2} & 0\end{pmatrix}
\end{aligned}
$$

行列の積の差を取る。

$$
\begin{aligned}
[X,Y]\big|_{m=2} &= \begin{pmatrix}- i & 0\\0 & i\end{pmatrix}
\end{aligned}
$$

trace は零になる。

$$
\begin{aligned}
\operatorname{tr}[X,Y] &= 0
\end{aligned}
$$
:::

:::example[Pauli 基底の交換子（3）]
反 Hermite な生成元を二つ選び、第二次数に残る項を調べる。

$$
m=3
$$

一方の生成元の強さを変える。

$$
\begin{aligned}
X\big|_{m=3} &= \begin{pmatrix}0 & - \frac{3 i}{2}\\- \frac{3 i}{2} & 0\end{pmatrix}
\end{aligned}
$$

他方は固定する。

$$
\begin{aligned}
Y &= \begin{pmatrix}0 & - \frac{1}{2}\\\frac{1}{2} & 0\end{pmatrix}
\end{aligned}
$$

行列の積の差を取る。

$$
\begin{aligned}
[X,Y]\big|_{m=3} &= \begin{pmatrix}- \frac{3 i}{2} & 0\\0 & \frac{3 i}{2}\end{pmatrix}
\end{aligned}
$$

trace は零になる。

$$
\begin{aligned}
\operatorname{tr}[X,Y] &= 0
\end{aligned}
$$
:::

:::example[Pauli 基底の交換子（4）]
反 Hermite な生成元を二つ選び、第二次数に残る項を調べる。

$$
m=4
$$

一方の生成元の強さを変える。

$$
\begin{aligned}
X\big|_{m=4} &= \begin{pmatrix}0 & - 2 i\\- 2 i & 0\end{pmatrix}
\end{aligned}
$$

他方は固定する。

$$
\begin{aligned}
Y &= \begin{pmatrix}0 & - \frac{1}{2}\\\frac{1}{2} & 0\end{pmatrix}
\end{aligned}
$$

行列の積の差を取る。

$$
\begin{aligned}
[X,Y]\big|_{m=4} &= \begin{pmatrix}- 2 i & 0\\0 & 2 i\end{pmatrix}
\end{aligned}
$$

trace は零になる。

$$
\begin{aligned}
\operatorname{tr}[X,Y] &= 0
\end{aligned}
$$
:::

### 指数写像の局所性

行列群で定数の生成元 $X$ を取ると、$U'(t)=XU(t)$、$U(0)=I$ の解は行列指数である。指数級数を項別微分して方程式を満たすことが分かり、積 $e^{sX}e^{tX}=e^{(s+t)X}$ は同じ冪をまとめる二項公式で得られる。逆に微分可能な一径数部分群 $g(s+t)=g(s)g(t)$ で $s$ を $0$ で微分すると、同じ初期値問題を得る。このため一径数部分群は生成元から一意に決まる。

例えば $\SU(2)$ の生成元 $X$ が $X^*=-X$ と $\tr X=0$ を満たすなら、$(e^{tX})^*e^{tX}=I$ である。行列式が $1$ になることは、$X$ の二つの固有値が互いに反対であることから従う。逆に $\SU(2)$ の曲線を単位元で微分すれば、ユニタリ性と行列式の制限からこの二条件を得る。指数写像は単位元の近くでは微分が恒等であるため局所座標を与えるが、群全体では異なる生成元の指数が一致する。例えば回転角の周期による一致は局所座標に含めて考えられない。


:::example[実回転の接行列]
三次元で第 3 軸のまわりの等角速度を選ぶ。接行列は第 3 成分を変えず、平面内では位置と直交する速度を与える。

座標成分を列に並べる。この順序に合わせ、線形写像の行列の列を選ぶ。

$$
\begin{aligned}v&=\begin{pmatrix}2\\1\\3\end{pmatrix}\\ A&=\begin{pmatrix}0 & -2 & 0\\2 & 0 & 0\\0 & 0 & 0\end{pmatrix}\end{aligned}
$$

行列の各行と入力の内積を計算する。例えば第一行では $\begin{pmatrix}0 & -2 & 0\end{pmatrix}$ と入力列を掛ける。

$$
\begin{aligned}(Av)_1&=-2\\Av&=\begin{pmatrix}-2\\4\\0\end{pmatrix}\end{aligned}
$$

核の次元は列数から階数を引いて得られる。像の次元は階数である。

$$
\begin{aligned}\operatorname{rank}A&=2\\\dim\ker A&=1\end{aligned}
$$

正方行列の場合は、特性多項式を求めて固有値を調べられる。ここでは $tI-A$ の行列式を取る。

$$
\begin{aligned}\operatorname{tr}A&=0\\ \det A&=0\\ \det(tI-A)&=t^{3} + 4 t\end{aligned}
$$

行列要素を得るため、出力から指定した共変ベクトルで成分を取り出す。

$$
\begin{aligned}\ell&=\begin{pmatrix}2 & 1 & 3\end{pmatrix}\\ \ell Av&=0\end{aligned}
$$

次に別の方向を入力にする。出力の差は入力の差に $A$ を掛けたものである。

$$
\begin{aligned}w&=\begin{pmatrix}0\\0\\1\end{pmatrix}\\ Aw&=\begin{pmatrix}0\\0\\0\end{pmatrix}\\ Aw-Av&=\begin{pmatrix}2\\-4\\0\end{pmatrix}\end{aligned}
$$
:::

## 回転群 SO(3) と SU(2)

### 二成分から三次元の回転を作る

Pauli 行列 $\sigma_1,\sigma_2,\sigma_3$ を
$$
\sigma_1=\begin{pmatrix}0&1\\1&0\end{pmatrix},\quad
\sigma_2=\begin{pmatrix}0&-i\\i&0\end{pmatrix},\quad
\sigma_3=\begin{pmatrix}1&0\\0&-1\end{pmatrix}
$$
とする。実ベクトル $x$ を $X=x_i\sigma_i$ に対応させると $\det X=-|x|^2$ である。$U\in\SU(2)$ で $UXU^\dagger$ は Hermitian かつトレース零で、行列式が保存される。従って $\R^3$ に直交変換を誘導する。$U=\begin{pmatrix}a&b\\-\bar b&\bar a\end{pmatrix}$、$|a|^2+|b|^2=1$ の形から $\SU(2)$ は連結なので、誘導する変換は行列式 $+1$ である。全 $\sigma_i$ と交換する $U$ は対角で $\sigma_1$ とも交換するからスカラー、従って核は $\{\pm I\}$ である。$\sigma_i\sigma_j=\delta_{ij}I+i[ijk]\sigma_k$ をすべて直接掛け算で確認できる。この式によって $U=\cos(\theta/2)I-i\sin(\theta/2)\vb n\cdot\sigma$ は軸 $\vb n$ の角 $\theta$ の回転を誘導する。任意の三次元の固有回転には固有軸があるので写像は全射である。各回転に二つの元が対応する。

例えば $z$ 軸の $2\pi$ 回転は $U=-I$、$4\pi$ で $I$ になる。スピン $1/2$ のベクトルは $2\pi$ で符号を変えるが、共役で記述した三次元ベクトルは変わらない。座標ベクトルとスピノルの変換はこの二重被覆で結び付く。

### 共役作用で三次元の角度を読む

$X=\vb x\cdot\vb\sigma$ を Hermite かつ trace 零の行列とする。$X\mapsto UXU^\dagger$ は $\frac12\operatorname{tr}(X^2)$ を保つ。$U=\diag(e^{-i\theta/2},e^{i\theta/2})$ では上三角成分 $x-iy$ に $e^{-i\theta}$ が掛かるから、三次元では $x+iy$ が $e^{i\theta}$ 倍になる。この作用では $U$ と $-U$ が同じ回転を与える。

:::example[半角から回転を組み立てる（1）]
角度 $m\pi/2$ の回転を二成分と三成分で比べる。

$$
m=1
$$

スピノルの行列に半角を入れる。

$$
\begin{aligned}
U\big|_{m=1} &= \begin{pmatrix}e^{- \frac{i \pi}{4}} & 0\\0 & e^{\frac{i \pi}{4}}\end{pmatrix}
\end{aligned}
$$

三成分の回転では座標角がそのまま入る。

$$
\begin{aligned}
R\big|_{m=1} &= \begin{pmatrix}0 & -1 & 0\\1 & 0 & 0\\0 & 0 & 1\end{pmatrix}
\end{aligned}
$$

同じ操作を二度行ったときのスピノルの行列を求める。

$$
\begin{aligned}
U^2\big|_{m=1} &= \begin{pmatrix}- i & 0\\0 & i\end{pmatrix}
\end{aligned}
$$
:::

:::example[半角から回転を組み立てる（2）]
角度 $m\pi/2$ の回転を二成分と三成分で比べる。

$$
m=2
$$

スピノルの行列に半角を入れる。

$$
\begin{aligned}
U\big|_{m=2} &= \begin{pmatrix}- i & 0\\0 & i\end{pmatrix}
\end{aligned}
$$

三成分の回転では座標角がそのまま入る。

$$
\begin{aligned}
R\big|_{m=2} &= \begin{pmatrix}-1 & 0 & 0\\0 & -1 & 0\\0 & 0 & 1\end{pmatrix}
\end{aligned}
$$

同じ操作を二度行ったときのスピノルの行列を求める。

$$
\begin{aligned}
U^2\big|_{m=2} &= \begin{pmatrix}-1 & 0\\0 & -1\end{pmatrix}
\end{aligned}
$$
:::

:::example[半角から回転を組み立てる（3）]
角度 $m\pi/2$ の回転を二成分と三成分で比べる。

$$
m=3
$$

スピノルの行列に半角を入れる。

$$
\begin{aligned}
U\big|_{m=3} &= \begin{pmatrix}e^{- \frac{3 i \pi}{4}} & 0\\0 & e^{\frac{3 i \pi}{4}}\end{pmatrix}
\end{aligned}
$$

三成分の回転では座標角がそのまま入る。

$$
\begin{aligned}
R\big|_{m=3} &= \begin{pmatrix}0 & 1 & 0\\-1 & 0 & 0\\0 & 0 & 1\end{pmatrix}
\end{aligned}
$$

同じ操作を二度行ったときのスピノルの行列を求める。

$$
\begin{aligned}
U^2\big|_{m=3} &= \begin{pmatrix}i & 0\\0 & - i\end{pmatrix}
\end{aligned}
$$
:::

:::example[半角から回転を組み立てる（4）]
角度 $m\pi/2$ の回転を二成分と三成分で比べる。

$$
m=4
$$

スピノルの行列に半角を入れる。

$$
\begin{aligned}
U\big|_{m=4} &= \begin{pmatrix}-1 & 0\\0 & -1\end{pmatrix}
\end{aligned}
$$

三成分の回転では座標角がそのまま入る。

$$
\begin{aligned}
R\big|_{m=4} &= \begin{pmatrix}1 & 0 & 0\\0 & 1 & 0\\0 & 0 & 1\end{pmatrix}
\end{aligned}
$$

同じ操作を二度行ったときのスピノルの行列を求める。

$$
\begin{aligned}
U^2\big|_{m=4} &= \begin{pmatrix}1 & 0\\0 & 1\end{pmatrix}
\end{aligned}
$$
:::

### 回転の被覆と四元数の積

任意の $U\in\SU(2)$ は $U=a I-i\vb b\cdot\vb\sigma$ と書ける。$a,\vb b$ は実で $a^2+|\vb b|^2=1$ である。実際、第一列を $z,w$ とすると、第二列はユニタリ性と行列式から $-\overline w,\overline z$ となる。その四実成分を実部と虚部に分けるとこの形になる。Pauli 行列の積 $\sigma_i\sigma_j=\delta_{ij}I+i\epsilon_{ijk}\sigma_k$ を使うと、二つの組の積は
$$
(a,\vb b)(c,\vb d)=(ac-\vb b\cdot\vb d,\ a\vb d+c\vb b+\vb b\times\vb d)
$$
である。従って単位四元数として回転を合成できる。$\theta$ 回転に対し $a=\cos(\theta/2),\vb b=\vb n\sin(\theta/2)$ を選べば、主軸を $\vb n$ とする回転になる。共役作用が恒等ならすべての Pauli 行列と可換である。$\sigma_z$ と可換な行列は対角であり、さらに $\sigma_x$ と可換なら対角成分が等しいので $\pm I$ に限る。どの三次元回転も軸角で表せるため、この共役写像は全射で核は $\{\pm I\}$ である。


:::example[四半回転を順に合成する]
第 1 軸の四半回転のあと第 3 軸の四半回転を行う。合成行列の固定方向は対角線であり、二つの操作の順番を交換すると別の回転になる。

座標成分を列に並べる。この順序に合わせ、線形写像の行列の列を選ぶ。

$$
\begin{aligned}v&=\begin{pmatrix}1\\0\\0\end{pmatrix}\\ A&=\begin{pmatrix}0 & 0 & 1\\1 & 0 & 0\\0 & 1 & 0\end{pmatrix}\end{aligned}
$$

行列の各行と入力の内積を計算する。例えば第一行では $\begin{pmatrix}0 & 0 & 1\end{pmatrix}$ と入力列を掛ける。

$$
\begin{aligned}(Av)_1&=0\\Av&=\begin{pmatrix}0\\1\\0\end{pmatrix}\end{aligned}
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
\begin{aligned}\ell&=\begin{pmatrix}1 & 1 & 1\end{pmatrix}\\ \ell Av&=1\end{aligned}
$$

次に別の方向を入力にする。出力の差は入力の差に $A$ を掛けたものである。

$$
\begin{aligned}w&=\begin{pmatrix}1\\1\\1\end{pmatrix}\\ Aw&=\begin{pmatrix}1\\1\\1\end{pmatrix}\\ Aw-Av&=\begin{pmatrix}1\\0\\1\end{pmatrix}\end{aligned}
$$
:::

## SO(3) の既約表現と球面調和関数

### 昇降演算子と球面のモード

$L_i=J_i/\hbar$、$L_\pm=L_x\pm iL_y$ とする。回転の生成子から $[L_i,L_j]=i[ijk]L_k$ であり、$[L_z,L_\pm]=\pm L_\pm$、$[L_+,L_-]=2L_z$ である。Casimir $L^2=\sum L_i^2$ は括弧を積に展開すると各 $L_j$ と交換する。既約な有限次元ユニタリ表現で $L_z$ の最大固有値を $j$ とすれば、その固有ベクトルは $L_+$ に消される。$L^2=L_-L_++L_z(L_z+1)$ より値は $j(j+1)$ である。$L_-$ のノルムは $\|L_-|j,m\rangle\|^2=j(j+1)-m(m-1)=(j+m)(j-m+1)$ であり、有限次元では $m=-j$ で止まる。従って $2j$ は非負整数で、次元 $2j+1$ になる。この鎖は全生成子で閉じるため既約性で全空間である。逆に任意の $2j=n$ は $\C^2$ の $n$ 次対称冪で、最高ウェイトからの鎖が次元 $n+1$ を尽くすことにより実現する。$-I$ は対称冪で $(-1)^n$ なので $\SO(3)$ へ降りるのは整数 $j$ である。

球面では回転の微分が $L_z=-i\partial_\phi$、$L^2=-\Delta_{S^2}$ になる。後者は $\vb L=-i\vb r\times\grad$ を二乗し、交代記号を縮約すると $L^2=-r^2\Delta+(r\partial_r)^2+r\partial_r$ となるためである。角度だけの関数では動径微分は零である。従って :ref[def-mp-spherical-harmonic] の次数 $\ell$ の空間は $j=\ell$ で、最高ウェイト $Y_{\ell\ell}$ から降ろすと全ての $m$ が現れる。例えば $\ell=1$ は三成分のベクトル表現である。

### 昇降演算子の積と系列の端

$J_\pm=J_x\pm iJ_y$ とすれば $J_+J_-=J^2-J_z(J_z-\hbar I)$ である。これは積を展開して $[J_x,J_y]=i\hbar J_z$ を使えばよい。従って規格化された状態に $J_-$ を作用させた二乗ノルムは $\hbar^2(j+m)(j-m+1)$ になる。無次元の $L_-$ には前の $\hbar^2$ が付かない。系列の下端ではいずれも零である。球の回転では角度の周期が $2\pi$ なので $m$ が整数になり、最高ウェイトも整数になる。

:::example[最高ウェイトから二段下げる（1）]
スピン $j=m$ とし、最高状態から一つずつ降りる。

$$
m=1
$$

第一段階の係数の二乗を求める。

$$
\begin{aligned}
\|L_-\ket{j,j}\|^2\big|_{m=1} &= 2
\end{aligned}
$$

規格化してから第二段を作用させると別の係数が出る。

$$
\begin{aligned}
\|L_-\ket{j,j-1}\|^2\big|_{m=1} &= 2
\end{aligned}
$$

規格化を途中で行わない二段の作用では積になる。

$$
\begin{aligned}
\|L_-^2\ket{j,j}\|^2\big|_{m=1} &= 4
\end{aligned}
$$
:::

:::example[最高ウェイトから二段下げる（2）]
スピン $j=m$ とし、最高状態から一つずつ降りる。

$$
m=2
$$

第一段階の係数の二乗を求める。

$$
\begin{aligned}
\|L_-\ket{j,j}\|^2\big|_{m=2} &= 4
\end{aligned}
$$

規格化してから第二段を作用させると別の係数が出る。

$$
\begin{aligned}
\|L_-\ket{j,j-1}\|^2\big|_{m=2} &= 6
\end{aligned}
$$

規格化を途中で行わない二段の作用では積になる。

$$
\begin{aligned}
\|L_-^2\ket{j,j}\|^2\big|_{m=2} &= 24
\end{aligned}
$$
:::

:::example[最高ウェイトから二段下げる（3）]
スピン $j=m$ とし、最高状態から一つずつ降りる。

$$
m=3
$$

第一段階の係数の二乗を求める。

$$
\begin{aligned}
\|L_-\ket{j,j}\|^2\big|_{m=3} &= 6
\end{aligned}
$$

規格化してから第二段を作用させると別の係数が出る。

$$
\begin{aligned}
\|L_-\ket{j,j-1}\|^2\big|_{m=3} &= 10
\end{aligned}
$$

規格化を途中で行わない二段の作用では積になる。

$$
\begin{aligned}
\|L_-^2\ket{j,j}\|^2\big|_{m=3} &= 60
\end{aligned}
$$
:::

:::example[最高ウェイトから二段下げる（4）]
スピン $j=m$ とし、最高状態から一つずつ降りる。

$$
m=4
$$

第一段階の係数の二乗を求める。

$$
\begin{aligned}
\|L_-\ket{j,j}\|^2\big|_{m=4} &= 8
\end{aligned}
$$

規格化してから第二段を作用させると別の係数が出る。

$$
\begin{aligned}
\|L_-\ket{j,j-1}\|^2\big|_{m=4} &= 14
\end{aligned}
$$

規格化を途中で行わない二段の作用では積になる。

$$
\begin{aligned}
\|L_-^2\ket{j,j}\|^2\big|_{m=4} &= 112
\end{aligned}
$$
:::

### 調和多項式で既約成分を作る

三変数の斉次多項式 $P_\ell$ に対して、$P_\ell=H_\ell+r^2Q_{\ell-2}$ と分ける。ここで $H_\ell$ は調和である。分解を帰納法で作るには、次数 $\ell-2$ の多項式に対しすでに調和成分への分解を持つと仮定する。$\Delta(r^{2a}H_j)=2a(2j+2a+1)r^{2a-2}H_j$ は積の微分と $x\cdot\nabla H_j=jH_j$ で確かめられ、$a>0$ なら係数が非零である。従って $\Delta P_\ell$ を低い次数で分解し、それぞれに $r^2$ を掛けて係数で割れば、$r^2Q_{\ell-2}$ の Laplacian を合わせられる。唯一性も、同じ係数が零でないことと低い次数での分解の唯一性から従う。この分解により $\dim H_\ell=\binom{\ell+2}{2}-\binom{\ell}{2}=2\ell+1$ になる。

例えば二次の $z^2$ は $z^2-r^2/3+r^2/3$ と分かれる。多項式 $z^2-r^2/3$ は $\Delta z^2=2$ と $\Delta r^2=6$ により調和である。球面へ制限すれば $x=\cos\theta$ として $x^2-1/3=(2/3)P_2(x)$ になる。続いて三次の $z^3$ では $\Delta z^3=6z$、$\Delta(r^2z)=10z$ なので $z^3-(3/5)r^2z$ が調和である。その制限は $(2/5)P_3(x)$ である。空間の多項式をそのまま球面の単独モードと見なすと、上の低次成分を落とすことになる。

## 表現のテンソル積と Clebsch–Gordan 係数

### 直積を全角運動量で分ける

二つの空間のテンソル積では $L_i=L_i^{(1)}\otimes I+I\otimes L_i^{(2)}$ であり、異なる因子は交換する。従って $M=m_1+m_2$ である。$j_1\ge j_2$ とすると $M\ge0$ の重みの多重度は、$m_1+m_2=M$ を満たす格子点数である。$M$ を最大 $j_1+j_2$ から $j_1-j_2$ まで一つずつ下げると多重度が一つずつ増え、以後一定である。各既約成分は範囲内の各 $M$ を一度ずつ持つので、隣接多重度の差を取ると
$$
V_{j_1}\otimes V_{j_2}=\bigoplus_{J=|j_1-j_2|}^{j_1+j_2} V_J
\label{eq-mp-spin-addition}
$$
を得る。$J$ の刻みは $1$ である。内積が不変なので全既約成分へ直交分解できる。$C^{JM}_{j_1m_1,j_2m_2}$ は積基底に対する $|J,M\rangle$ の展開係数で、最高ウェイトの位相を一つ定め、全体の $L_-$ で他の状態を得る。

:::example[二つの半整数スピン]
最高状態 $|\uparrow\uparrow\rangle$ を一度降ろすと $|\downarrow\uparrow\rangle+|\uparrow\downarrow\rangle$ であり、$\sqrt2$ で割ると $J=1,M=0$ になる。その直交補の $M=0$ 成分 $(|\uparrow\downarrow\rangle-|\downarrow\uparrow\rangle)/\sqrt2$ は $L_\pm$ の両方で零で $J=0$ になる。他の三重項の成分は $|\uparrow\uparrow\rangle,|\downarrow\downarrow\rangle$ である。
:::

:::example[スピン 1 と 1/2 の合成係数]
最高状態 $|1,1\rangle|\uparrow\rangle$ は $|3/2,3/2\rangle$ である。全 $L_-$ を作用させると $\sqrt2\,|1,0\rangle|\uparrow\rangle+|1,1\rangle|\downarrow\rangle$ になる。ノルムが $\sqrt3$ なので、$|3/2,1/2\rangle$ の係数はそれぞれ $\sqrt{2/3},1/\sqrt3$ である。直交する最高ウェイトを
$$
|1/2,1/2\rangle
=\sqrt{2/3}|1,1\rangle|\downarrow\rangle-\frac1{\sqrt3}|1,0\rangle|\uparrow\rangle
$$
と選ぶ。この二項に $L_+$ を掛けると、$|1,1\rangle|\uparrow\rangle$ の係数が $\sqrt{2/3}$ と $-\sqrt2/\sqrt3$ になって消える。従って第二の状態がスピン $1/2$ の最高ウェイトであることを位相込みで確かめられる。
:::

### 最高状態を一段下げて合成係数を出す

$j>0$ と $1/2$ の積で最高状態は $\ket{j,j}\ket{+}$ である。全降下演算子を一度作用させると、第一因子を降ろした項に $\sqrt{2j}$、第二因子を降ろした項に $1$ が掛かる。二項は直交するので全体のノルムは $\sqrt{2j+1}$ である。同じ磁気量子数においてその直交補空間を取ると、もう一つの既約成分 $J=j-1/2$ の最高状態になる。$j=0$ では第一因子を降ろせず、直交補に当たる成分は現れない。

一般の $J=j_1+j_2-n$ に対し、$J$ の最高状態を
$$
\ket{J,J}=\sum_{a=0}^{n}c_a
\ket{j_1,j_1-a}\otimes\ket{j_2,j_2-n+a}
$$
と書く。$0\le n\le 2\min(j_1,j_2)$ なので全添字が許される。$L_+$ がこの状態を消す条件を、第一因子を $j_1-a$ に固定して読むと
$$
c_{a+1}\sqrt{(a+1)(2j_1-a)}
+c_a\sqrt{(n-a)(2j_2-n+a+1)}=0
\qquad (0\le a<n)
$$
となる。従って漸化式と $\sum_a|c_a|^2=1$ で最高状態が求まる。本章では Condon–Shortley の位相に合わせ $c_0>0$ とし、そこから $L_-$ で降ろす際は正の規格化係数を使う。この直交基底変換は実であり、各 $M$ のブロックで係数行列は直交行列である。従って $C$ の行どうしと列どうしはそれぞれ正規直交する。

### Wigner の $3j$ 記号と結合順序

同じ規約で Wigner の $3j$ 記号は
$$
\begin{pmatrix}j_1&j_2&J\\ m_1&m_2&-M\end{pmatrix}
=\frac{(-1)^{j_1-j_2+M}}{\sqrt{2J+1}}C^{JM}_{j_1m_1,j_2m_2}
$$
と関係する。特に下面の三成分の和が零であることが必要である。三つのスピンを合成するときには、最初の二つを結合してから第三のスピンを加える基底と、後ろの二つを先に結合する基底がある。固定 $J,M$ における変換行列は両基底の内積であり、例えば
$$
\sum_{m_1,m_2,m_3}
C^{j_{12},m_1+m_2}_{j_1m_1,j_2m_2}
C^{JM}_{j_{12},m_1+m_2,j_3m_3}
C^{j_{23},m_2+m_3}_{j_2m_2,j_3m_3}
C^{JM}_{j_1m_1,j_{23},m_2+m_3}
$$
で計算できる。この式は二つの基底ベクトルを同じ積基底に展開した結果である。全ての $M$ が全 $L_-$ で同じ規格化により作られるので、この内積は $M$ に依存しない。和の値を次の形に書くと、Wigner の $6j$ 記号の規約が定まる。
$$
(-1)^{j_1+j_2+j_3+J}\sqrt{(2j_{12}+1)(2j_{23}+1)}
\begin{Bmatrix}j_1&j_2&j_{12}\\ j_3&J&j_{23}\end{Bmatrix}.
$$
例えば三つとも $1/2$、全スピン $J=M=1/2$、中間の二つのスピン $j_{12}=j_{23}=0$ とする。最初の結合では $(|\uparrow\downarrow\uparrow\rangle-|\downarrow\uparrow\uparrow\rangle)/\sqrt2$、後の結合では $(|\uparrow\uparrow\downarrow\rangle-|\uparrow\downarrow\uparrow\rangle)/\sqrt2$ なので、内積は $-1/2$ である。この場合は位相も平方根も $1$ であり、$6j$ 記号は $-1/2$ となる。

:::example[一段下の合成係数（1）]
基底の順を $(\ket{j,j-1}\ket{+},\ket{j,j}\ket{-})$ とする。

$$
m=1
$$

高い方の全スピンの成分を規格化する。

$$
\begin{aligned}
v_+\big|_{m=1} &= \begin{pmatrix}\frac{\sqrt{6}}{3}\\\frac{\sqrt{3}}{3}\end{pmatrix}
\end{aligned}
$$

低い方の全スピンには直交する列を選ぶ。ここでの全体位相は左の成分を負にする規約である。

$$
\begin{aligned}
v_-\big|_{m=1} &= \begin{pmatrix}- \frac{\sqrt{3}}{3}\\\frac{\sqrt{6}}{3}\end{pmatrix}
\end{aligned}
$$

例えば高い方で第二因子が下向きの確率は第二係数の二乗である。

$$
\begin{aligned}
P_-\big|_{m=1} &= \frac{1}{3}
\end{aligned}
$$
:::

:::example[一段下の合成係数（2）]
基底の順を $(\ket{j,j-1}\ket{+},\ket{j,j}\ket{-})$ とする。

$$
m=2
$$

高い方の全スピンの成分を規格化する。

$$
\begin{aligned}
v_+\big|_{m=2} &= \begin{pmatrix}\frac{2 \sqrt{5}}{5}\\\frac{\sqrt{5}}{5}\end{pmatrix}
\end{aligned}
$$

低い方の全スピンには直交する列を選ぶ。ここでの全体位相は左の成分を負にする規約である。

$$
\begin{aligned}
v_-\big|_{m=2} &= \begin{pmatrix}- \frac{\sqrt{5}}{5}\\\frac{2 \sqrt{5}}{5}\end{pmatrix}
\end{aligned}
$$

例えば高い方で第二因子が下向きの確率は第二係数の二乗である。

$$
\begin{aligned}
P_-\big|_{m=2} &= \frac{1}{5}
\end{aligned}
$$
:::

:::example[一段下の合成係数（3）]
基底の順を $(\ket{j,j-1}\ket{+},\ket{j,j}\ket{-})$ とする。

$$
m=3
$$

高い方の全スピンの成分を規格化する。

$$
\begin{aligned}
v_+\big|_{m=3} &= \begin{pmatrix}\frac{\sqrt{42}}{7}\\\frac{\sqrt{7}}{7}\end{pmatrix}
\end{aligned}
$$

低い方の全スピンには直交する列を選ぶ。ここでの全体位相は左の成分を負にする規約である。

$$
\begin{aligned}
v_-\big|_{m=3} &= \begin{pmatrix}- \frac{\sqrt{7}}{7}\\\frac{\sqrt{42}}{7}\end{pmatrix}
\end{aligned}
$$

例えば高い方で第二因子が下向きの確率は第二係数の二乗である。

$$
\begin{aligned}
P_-\big|_{m=3} &= \frac{1}{7}
\end{aligned}
$$
:::

:::example[一段下の合成係数（4）]
基底の順を $(\ket{j,j-1}\ket{+},\ket{j,j}\ket{-})$ とする。

$$
m=4
$$

高い方の全スピンの成分を規格化する。

$$
\begin{aligned}
v_+\big|_{m=4} &= \begin{pmatrix}\frac{2 \sqrt{2}}{3}\\\frac{1}{3}\end{pmatrix}
\end{aligned}
$$

低い方の全スピンには直交する列を選ぶ。ここでの全体位相は左の成分を負にする規約である。

$$
\begin{aligned}
v_-\big|_{m=4} &= \begin{pmatrix}- \frac{1}{3}\\\frac{2 \sqrt{2}}{3}\end{pmatrix}
\end{aligned}
$$

例えば高い方で第二因子が下向きの確率は第二係数の二乗である。

$$
\begin{aligned}
P_-\big|_{m=4} &= \frac{1}{9}
\end{aligned}
$$
:::

:::example[合成スピン 2 の磁気量子数 0 の成分]
積の二つのスピンをそれぞれ $1$ と $1$ にする。全演算子の磁気量子数は $0$ に固定する。

基底を、磁気量子数の和を満たす組だけに絞る。

$$
\begin{aligned}
(m_1,m_2)_{1} &= \left( 1, \  -1\right) \\
(m_1,m_2)_{2} &= \left( 0, \  0\right) \\
(m_1,m_2)_{3} &= \left( -1, \  1\right)
\end{aligned}
$$

同じ磁気量子数の中で、上の全スピンの系列を最高状態から降ろして作る。各段で規格化し、低い系列はその直交補を取る。係数を Condon–Shortley の位相で並べる。

$$
\begin{aligned}
C_{1} &= \frac{\sqrt{6}}{6} \\
C_{2} &= \frac{\sqrt{6}}{3} \\
C_{3} &= \frac{\sqrt{6}}{6}
\end{aligned}
$$

最後に積基底が正規直交であることから、成分の確率を足して確かめる。

$$
\begin{aligned}
|C_{1}|^2 &= \frac{1}{6} \\
|C_{2}|^2 &= \frac{2}{3} \\
|C_{3}|^2 &= \frac{1}{6} \\
\sum_i |C_i|^2 &= 1
\end{aligned}
$$
:::


:::example[合成スピン 0 の磁気量子数 0 の成分]
積の二つのスピンをそれぞれ $1$ と $1$ にする。全演算子の磁気量子数は $0$ に固定する。

基底を、磁気量子数の和を満たす組だけに絞る。

$$
\begin{aligned}
(m_1,m_2)_{1} &= \left( 1, \  -1\right) \\
(m_1,m_2)_{2} &= \left( 0, \  0\right) \\
(m_1,m_2)_{3} &= \left( -1, \  1\right)
\end{aligned}
$$

同じ磁気量子数の中で、上の全スピンの系列を最高状態から降ろして作る。各段で規格化し、低い系列はその直交補を取る。係数を Condon–Shortley の位相で並べる。

$$
\begin{aligned}
C_{1} &= \frac{\sqrt{3}}{3} \\
C_{2} &= - \frac{\sqrt{3}}{3} \\
C_{3} &= \frac{\sqrt{3}}{3}
\end{aligned}
$$

最後に積基底が正規直交であることから、成分の確率を足して確かめる。

$$
\begin{aligned}
|C_{1}|^2 &= \frac{1}{3} \\
|C_{2}|^2 &= \frac{1}{3} \\
|C_{3}|^2 &= \frac{1}{3} \\
\sum_i |C_i|^2 &= 1
\end{aligned}
$$
:::


:::example[合成スピン 3/2 の磁気量子数 1/2 の成分]
積の二つのスピンをそれぞれ $2$ と $\frac{1}{2}$ にする。全演算子の磁気量子数は $\frac{1}{2}$ に固定する。

基底を、磁気量子数の和を満たす組だけに絞る。

$$
\begin{aligned}
(m_1,m_2)_{1} &= \left( 1, \  - \frac{1}{2}\right) \\
(m_1,m_2)_{2} &= \left( 0, \  \frac{1}{2}\right)
\end{aligned}
$$

同じ磁気量子数の中で、上の全スピンの系列を最高状態から降ろして作る。各段で規格化し、低い系列はその直交補を取る。係数を Condon–Shortley の位相で並べる。

$$
\begin{aligned}
C_{1} &= \frac{\sqrt{15}}{5} \\
C_{2} &= - \frac{\sqrt{10}}{5}
\end{aligned}
$$

最後に積基底が正規直交であることから、成分の確率を足して確かめる。

$$
\begin{aligned}
|C_{1}|^2 &= \frac{3}{5} \\
|C_{2}|^2 &= \frac{2}{5} \\
\sum_i |C_i|^2 &= 1
\end{aligned}
$$
:::

## テンソル演算子と Wigner–Eckart の定理

### 方向依存を一つの係数へまとめる

$T^{(k)}_q$（$q=-k,\ldots,k$）が共役作用で既約型 $V_k$ として変換されるとき、これらを :term[球面テンソル演算子]{en="spherical tensor operator" yomi="きゅうめんてんそるえんざんし"} という。回転操作を微分すると $[L_z,T_q]=qT_q$、$[L_\pm,T_q]=\sqrt{(k\mp q)(k\pm q+1)}T_{q\pm1}$ となる。$|j,m\rangle$ に作用させた $T_q$ を $V_{j'}$ へ射影する写像 $V_j\otimes V_k\to V_{j'}$ は群と可換である。式 :eqref[eq-mp-spin-addition] で各 $j'$ は高々一度現れ、既約間の可換写像は同型ならスカラー倍である。従って、位相を共通に選べば
$$
\langle j',m'|T^{(k)}_q|j,m\rangle
=\frac{\langle j'\|T^{(k)}\|j\rangle}{\sqrt{2j'+1}}\,
 C^{j'm'}_{jm,kq}
\label{eq-mp-wigner-eckart}
$$
となる。これは Wigner–Eckart の定理で、式は縮約行列要素の規格化も定める。従って $m'=m+q$ と三角条件が必要である。例えば $j=0$ から階数 $1$ の演算子を用いると $j'=1$ のみが残る。空間反転が対称で演算子が奇なら、同じ反転固有値の状態間では行列要素がさらに消える。回転の合成則に反転の条件を追加する必要がある。

### スカラーとベクトルの選択則

階数 $0$ の演算子は群作用と可換であり、既約表現内では Schur の補題により恒等の倍数になる。階数 $1$ の場合、$[J_z,T_q]=\hbar qT_q$ なので $q$ の成分は磁気量子数を $q$ だけずらす。同じ $j>0$ の一つのコピーに射影したベクトル演算子では、$V_j\otimes V_1$ に $V_j$ が一度出るため、行列要素は $J_q$ のものに比例する。$j=0$ のコピー内では階数 $1$ の行列要素は零である。初状態・終状態に同じ $j$ のコピーが複数ある場合には、コピーごとに縮約行列要素が変わり得る。

:::example[角運動量自体をテンソル演算子に取る（1）]
角運動量を無次元の $L=J/\hbar$ で表す。初期状態を $\ket{j,j-1}$、スピンを $j=m$ とする。

$$
m=1
$$

零成分の期待値は初期の磁気量子数である。

$$
\begin{aligned}
\mel{j,j-1}{L_0}{j,j-1}\big|_{m=1} &= 0
\end{aligned}
$$

球成分 $L_{+1}=-L_+/\sqrt2$ は最高状態へ移す。

$$
\begin{aligned}
\mel{j,j}{L_{+1}}{j,j-1}\big|_{m=1} &= -1
\end{aligned}
$$

もう一方の球成分 $L_{-1}=L_-/\sqrt2$ の係数も求める。

$$
\begin{aligned}
\mel{j,j-2}{L_{-1}}{j,j-1}\big|_{m=1} &= 1
\end{aligned}
$$
:::

:::example[角運動量自体をテンソル演算子に取る（2）]
角運動量を無次元の $L=J/\hbar$ で表す。初期状態を $\ket{j,j-1}$、スピンを $j=m$ とする。

$$
m=2
$$

零成分の期待値は初期の磁気量子数である。

$$
\begin{aligned}
\mel{j,j-1}{L_0}{j,j-1}\big|_{m=2} &= 1
\end{aligned}
$$

球成分 $L_{+1}=-L_+/\sqrt2$ は最高状態へ移す。

$$
\begin{aligned}
\mel{j,j}{L_{+1}}{j,j-1}\big|_{m=2} &= - \sqrt{2}
\end{aligned}
$$

もう一方の球成分 $L_{-1}=L_-/\sqrt2$ の係数も求める。

$$
\begin{aligned}
\mel{j,j-2}{L_{-1}}{j,j-1}\big|_{m=2} &= \sqrt{3}
\end{aligned}
$$
:::

:::example[角運動量自体をテンソル演算子に取る（3）]
角運動量を無次元の $L=J/\hbar$ で表す。初期状態を $\ket{j,j-1}$、スピンを $j=m$ とする。

$$
m=3
$$

零成分の期待値は初期の磁気量子数である。

$$
\begin{aligned}
\mel{j,j-1}{L_0}{j,j-1}\big|_{m=3} &= 2
\end{aligned}
$$

球成分 $L_{+1}=-L_+/\sqrt2$ は最高状態へ移す。

$$
\begin{aligned}
\mel{j,j}{L_{+1}}{j,j-1}\big|_{m=3} &= - \sqrt{3}
\end{aligned}
$$

もう一方の球成分 $L_{-1}=L_-/\sqrt2$ の係数も求める。

$$
\begin{aligned}
\mel{j,j-2}{L_{-1}}{j,j-1}\big|_{m=3} &= \sqrt{5}
\end{aligned}
$$
:::

:::example[角運動量自体をテンソル演算子に取る（4）]
角運動量を無次元の $L=J/\hbar$ で表す。初期状態を $\ket{j,j-1}$、スピンを $j=m$ とする。

$$
m=4
$$

零成分の期待値は初期の磁気量子数である。

$$
\begin{aligned}
\mel{j,j-1}{L_0}{j,j-1}\big|_{m=4} &= 3
\end{aligned}
$$

球成分 $L_{+1}=-L_+/\sqrt2$ は最高状態へ移す。

$$
\begin{aligned}
\mel{j,j}{L_{+1}}{j,j-1}\big|_{m=4} &= -2
\end{aligned}
$$

もう一方の球成分 $L_{-1}=L_-/\sqrt2$ の係数も求める。

$$
\begin{aligned}
\mel{j,j-2}{L_{-1}}{j,j-1}\big|_{m=4} &= \sqrt{7}
\end{aligned}
$$
:::

### 縮約行列要素の規約

:eqref[eq-mp-wigner-eckart] の規約では、同じ $j_i,k,j_f$ に対して、Clebsch–Gordan 係数を取り除いた比例定数が全磁気量子数に共通である。係数を $3j$ 記号で書く流儀では、行と列の順序で位相が変わり、さらに $\sqrt{2j_f+1}$ が比例定数の定義に入る場合がある。例えば本文の形で $T_q=J_q$ とし、$j_i=j_f=j$ を取る。最高状態で $q=0$ の行列要素は $\hbar j$ である。対応する合成係数を $\sqrt{j/(j+1)}$ とする位相を選べば、比例定数は $\hbar\sqrt{j(j+1)}$ である。式における縮約行列要素そのものには、さらに $\sqrt{2j+1}$ を掛ける。

$j=1$ ではこの比例定数は $\sqrt2\hbar$、縮約行列要素は $\sqrt6\hbar$ になる。初期 $m=0$ から終状態 $m=1$ へ移る $J_{+1}=-J_+/\sqrt2$ の行列要素は $-\hbar$ であり、合成係数は $-1/\sqrt2$ である。初期 $m=0$ から $m=-1$ の $J_{-1}=J_-/\sqrt2$ では逆に $\hbar$ になり、係数は $1/\sqrt2$ になる。一方 $q=0$ で $m=0$ のままでは行列要素が零である。よって三角条件と磁気量子数の条件を満たしても、対応する成分の係数が零なら遷移は消える。


:::example[階数一の上昇成分]
スピン $1$ の基底順を $m=1,0,-1$ とし、行列は $\hbar$ で割った $J_{+1}$ とする。

座標成分を列に並べる。この順序に合わせ、線形写像の行列の列を選ぶ。

$$
\begin{aligned}v&=\begin{pmatrix}0\\1\\0\end{pmatrix}\\ A&=\begin{pmatrix}0 & -1 & 0\\0 & 0 & -1\\0 & 0 & 0\end{pmatrix}\end{aligned}
$$

行列の各行と入力の内積を計算する。例えば第一行では $\begin{pmatrix}0 & -1 & 0\end{pmatrix}$ と入力列を掛ける。

$$
\begin{aligned}(Av)_1&=-1\\Av&=\begin{pmatrix}-1\\0\\0\end{pmatrix}\end{aligned}
$$

核の次元は列数から階数を引いて得られる。像の次元は階数である。

$$
\begin{aligned}\operatorname{rank}A&=2\\\dim\ker A&=1\end{aligned}
$$

正方行列の場合は、特性多項式を求めて固有値を調べられる。ここでは $tI-A$ の行列式を取る。

$$
\begin{aligned}\operatorname{tr}A&=0\\ \det A&=0\\ \det(tI-A)&=t^{3}\end{aligned}
$$

行列要素を得るため、出力から指定した共変ベクトルで成分を取り出す。

$$
\begin{aligned}\ell&=\begin{pmatrix}1 & 0 & 0\end{pmatrix}\\ \ell Av&=-1\end{aligned}
$$

次に別の方向を入力にする。出力の差は入力の差に $A$ を掛けたものである。

$$
\begin{aligned}w&=\begin{pmatrix}0\\0\\1\end{pmatrix}\\ Aw&=\begin{pmatrix}0\\-1\\0\end{pmatrix}\\ Aw-Av&=\begin{pmatrix}1\\-1\\0\end{pmatrix}\end{aligned}
$$
:::

## SU(3) と Young 図形

### 三成分を対称化し、縮約を除く

$\SU(3)$ の基本表現を $V=\C^3$ とする。反対称積 $\bigwedge^2V$ は体積形式 $v_1\wedge v_2\wedge v_3$ とのペアリングで $V^*$ と同型である。演算子の行列式が $1$ なのでこのペアリングを保つ。対称部分は次元 $6$ で、従って $3\otimes3=6\oplus\bar3$ である。$V\otimes V^*$ は線形作用素全体で、共役の作用を受ける。トレース部分とトレース零部分に分けると $3\otimes\bar3=1\oplus8$ となる。トレース零部分の既約性は、複素化した接代数 $\mathfrak{sl}_3$ に対する随伴不変部分がイデアルであることから従う。実際、対角行列の括弧で根の異なる非対角成分を分離できる。非零の非対角成分を一つ得れば、行列単位との括弧で全非対角成分と対角差を生成できる。対角成分だけの非零元でも、トレース零かつ非スカラーだからある行列単位との括弧で非対角成分を得る。

:term[ヤング図形]{en="Young diagram" yomi="やんぐずけい"} は添字の対称化を示す。行の置換で対称化し、その結果を列の置換で反対称化する。長さが四以上の列は三次元では零になる。長さ三の列は行列式表現を与え、$\SU(3)$ 上で自明なので列ごと除ける。残る行長 $p+q,q,0$ の隣接差が Dynkin ラベル $(p,q)$ である。実際、第 $i$ 行の全ての箱を $e_i$ で埋め、行と列の操作を適用すると、元のテンソルが非零係数で残る。上昇作用 $E_{i,i+1}$ を先に適用した項は、同じ $e_i$ が列に二つあるため反対称化で消える。従ってこのテンソルは最高ウェイトを持ち、二つの対角差が掛ける係数は行長の隣接差になる。[古典群の表現](/math/lie-groups/classical-representations/) ではこの構成を一般の階数で扱う。

この章では $p$ 個の $V$ 添字と $q$ 個の $V^*$ 添字を各々対称化して縮約を除く。後で示すとおり、その空間が既約表現 $(p,q)$ である。次元は
$$
\dim(p,q)=\binom{p+2}{2}\binom{q+2}{2}
-\binom{p+1}{2}\binom{q+1}{2}
=\frac{(p+1)(q+1)(p+q+2)}2
$$
である。縮約は対を一つずつ落とす全射である。多項式表示でその随伴は、定数の規格化を除けば $\sum_i z_i w_i$ を掛ける操作であり、多項式環に零因子がないので単射である。$p=0$ または $q=0$ では縮約項は無い。例えば $(1,1)$ は $8$、$(3,0)$ は $10$ である。

### 重みの座標から三つの状態を区別する

基本表現で対角の生成子を $T_3=\diag(1,-1,0)/2$、$Y=\diag(1,1,-2)/3$ とする。同時固有値 $(T_3,Y)$ は $(1/2,1/3),(-1/2,1/3),(0,-2/3)$ である。双対表現では符号が反転する。行列単位 $E_{ij}$ は第 $j$ 基底を第 $i$ 基底へ移すので、括弧による重みは第 $i$ の重みから第 $j$ の重みを引いたものである。例えば $E_{12}$ は $(1,0)$、$E_{23}$ は $(-1/2,1)$ になる。随伴表現では六つの非対角行列が非零の根を与え、二つの独立なトレース零対角行列が零ウェイトを与える。従って重みから数えても次元は八である。

### 対称二階の六成分が既約であること

$\operatorname{Sym}^2V$ の基底を $e_i e_j$（$i\le j$）とする。対角な $\SU(3)$ の元はこれらをそれぞれ別の指標で対角化する。同じ部分空間で二つの指標が等しいには、指数の差が $(1,1,1)$ の整数倍でなければならないが、指数の総和の差が零なので別の基底間では起こらない。従って不変部分空間は対角作用で射影して、一つの基底ベクトルを含む。例えば $e_2e_3$ なら $E_{12}$ で $e_1e_3$、$E_{13}$ で $e_1^2$ を得る。任意の単項式からこのように $e_1^2$ まで上げられ、逆に $E_{21},E_{31}$ で六つ全部を得られる。複素線形の不変部分空間は複素化した Lie 代数でも保たれるため、これらの行列単位を使える。従って六次元の対称部分は既約である。

### 対称テンソルと縮約の次元

$V=\C^3$ の対称 $p$ 階テンソルには、三つの成分が何回現れるかの組が対応する。非負の三数の和が $p$ になる組は $\binom{p+2}{2}$ 個である。これと双対の対称 $q$ 階を掛けると積の次元になる。縮約は上と下の一組を取り除く。例えば $q=1$ では $V^*$ の成分を多項式 $f_j$ に持たせると、作用は $\sum_j\partial_j f_j$ であり、任意の次数 $p-1$ のモノミアルに対して一変数の原始関数を作れば全射である。

### 縮約の核の既約性

$V$ と $V^*$ の基底をそれぞれ変数 $z_i,w_i$ で表し、二重次数 $(p,q)$ の多項式を考える。縮約は定数倍を除き $D=\sum_i\partial_{z_i}\partial_{w_i}$ である。$s=\sum_i z_i w_i$ を掛ける操作は、単項式の二乗ノルムを多重添字の階乗の積とする内積では $D$ の随伴になる。$s$ を掛ける操作は単射なので、$p,q>0$ のとき $D$ は $(p-1,q-1)$ の空間への全射である。

最高ウェイトのベクトルは $E_{12},E_{23},E_{13}$ に消される。ここで $E_{ij}$ は微分作用素 $z_i\partial_{z_j}-w_j\partial_{w_i}$ である。これらに消される多項式は、三つの生成子の指数作用で不変である。$z_1w_3\ne0$ の点で、$E_{12}$ の作用で $z_2$、$E_{23}$ の作用で $w_2$、最後に $E_{13}$ の作用で $z_3$ を順に零にできる。後の操作は先に零にした変数を保つ。この間 $z_1,w_3,s$ は不変であり、最後の $w_1$ は $s/z_1$ である。従って二重次数 $(p,q)$ の不変多項式は $z_1^{p-k}w_3^{q-k}s^k$ の線形結合である。評価時の $w_1$ の次数から $0\le k\le q$ となり、$z_1$ の負冪を含む項は多項式ではあり得ない。実際、最大の $k>p$ の項だけが $z_1\to0$ で最も強い極を持ち、$w_3,s$ を非零のまま保てる。

$D$ と $s$ の交換子は $3+\sum_i(z_i\partial_{z_i}+w_i\partial_{w_i})$ である。これを $s^k$ の各因子に用いると
$$
D(s^k z_1^{p-k}w_3^{q-k})
= k(p+q-k+2)s^{k-1}z_1^{p-k}w_3^{q-k}
$$
となる。右辺は $k\ge1$ で線形独立な項である。従って $D$ の核の最高ウェイトの直線は $z_1^p w_3^q$ だけであり、単純根に対応する対角差 $E_{11}-E_{22}$ と $E_{22}-E_{33}$ はこれにそれぞれ $p,q$ を掛ける。各非零既約成分は、対角ウェイトのうち $\diag(1,0,-1)$ の値が最大のものを取れば、上昇作用で消えるベクトルを持つ。ユニタリ表現は直交分解できるので、最高直線が一つであれば核は既約である。

任意の有限次元既約表現でも、群の平均でユニタリ内積を選び、同じ方法で最高ベクトル $v$ を選ぶ。最初の二座標と後の二座標に作用する $\SU(2)$ へそれぞれ制限すると、$v$ は対応する上昇作用で消える。従ってそれぞれの対角差の固有値 $p,q$ は、$\SU(2)$ の無次元 $L_z$ の最高固有値を二倍した非負整数である。

同じ $p,q$ の二つの既約表現は同型になる。これを示すため、Lie 代数の積を低下作用、対角差、上昇作用の順に並べる。隣り合う二つを交換するたびに交換子が増えるが、その語の長さは一つ短い。長さと順序の逆転数についての帰納法により、どの積もこの順の積の線形結合で書ける。上昇作用は $v$ を消し、対角差は $p,q$ を掛けるので、低下作用の積で作るベクトルは $v$ が生成する表現全体を張る。ここで既約性と群の連結性を使った。

ユニタリ内積では低下作用 $E_{ji}$ の随伴は $E_{ij}$ である。低下作用の積 $P,Q$ からできる内積 $\braket{Pv}{Qv}$ は $\braket{v}{P^\dagger Qv}$ に移せる。ここで $P^\dagger Q$ を先ほどの順に並べると、上昇作用を含む項は消え、低下作用だけが残る項は $\diag(1,0,-1)$ の固有値が低くなるため $v$ と直交する。従って内積は $p,q$ と $\|v\|$ だけで決まる。$\|v\|=1$ と規格化すれば、二つの表現の低下作用の積の間の対応は Gram 行列を保ち、同じ線形関係を保つ。対応は並べ替えの式によって全生成子と可換なので同型である。これで $(p,q)$ の核が全既約型を尽くす。[最高ウェイト理論](/math/lie-groups/highest-weight/) では、一般の半単純 Lie 代数に対して分類する。

:::example[下付き添字が一つある既約成分（1）]
上付きの対称階数を $m$ とし、縮約の像を除く。

$$
m=1
$$

縮約前の次元は上の対称空間と三次元を掛ける。

$$
\begin{aligned}
d_{\rm all}\big|_{m=1} &= 9
\end{aligned}
$$

像の次元は一段下の対称空間に等しい。

$$
\begin{aligned}
d_{\rm image}\big|_{m=1} &= 1
\end{aligned}
$$

核の次元を引き算する。

$$
\begin{aligned}
d_{\rm ker}\big|_{m=1} &= 8
\end{aligned}
$$
:::

:::example[下付き添字が一つある既約成分（2）]
上付きの対称階数を $m$ とし、縮約の像を除く。

$$
m=2
$$

縮約前の次元は上の対称空間と三次元を掛ける。

$$
\begin{aligned}
d_{\rm all}\big|_{m=2} &= 18
\end{aligned}
$$

像の次元は一段下の対称空間に等しい。

$$
\begin{aligned}
d_{\rm image}\big|_{m=2} &= 3
\end{aligned}
$$

核の次元を引き算する。

$$
\begin{aligned}
d_{\rm ker}\big|_{m=2} &= 15
\end{aligned}
$$
:::

:::example[下付き添字が一つある既約成分（3）]
上付きの対称階数を $m$ とし、縮約の像を除く。

$$
m=3
$$

縮約前の次元は上の対称空間と三次元を掛ける。

$$
\begin{aligned}
d_{\rm all}\big|_{m=3} &= 30
\end{aligned}
$$

像の次元は一段下の対称空間に等しい。

$$
\begin{aligned}
d_{\rm image}\big|_{m=3} &= 6
\end{aligned}
$$

核の次元を引き算する。

$$
\begin{aligned}
d_{\rm ker}\big|_{m=3} &= 24
\end{aligned}
$$
:::

:::example[下付き添字が一つある既約成分（4）]
上付きの対称階数を $m$ とし、縮約の像を除く。

$$
m=4
$$

縮約前の次元は上の対称空間と三次元を掛ける。

$$
\begin{aligned}
d_{\rm all}\big|_{m=4} &= 45
\end{aligned}
$$

像の次元は一段下の対称空間に等しい。

$$
\begin{aligned}
d_{\rm image}\big|_{m=4} &= 10
\end{aligned}
$$

核の次元を引き算する。

$$
\begin{aligned}
d_{\rm ker}\big|_{m=4} &= 35
\end{aligned}
$$
:::

### 二箱と三箱の表現を具体化する

$V^{\otimes2}$ を交換で対称な空間と反対称な空間に分ける。前者には三つの同じ添字と三つの異なる添字の組があり、次元は $6$ である。後者には相異なる添字だけが入り、次元は $3$ である。交代三形式 $\epsilon$ が $\SU(3)$ で不変なので、反対称二階テンソルを $\epsilon$ で縮約すると $V^*$ が得られる。この同型を使うと $3\otimes3=6\oplus\bar3$ になる。さらに $V\otimes V^*$ は行列の空間であり、共役でトレースを保つ。従って $3\otimes\bar3=1\oplus8$ であり、八次元の部分にはトレース零の行列を取る。

三つを掛ける場合、反対称部分からは $\bar3\otimes3=1\oplus8$ が得られる。最初の二添字が対称な $T^{abd}$ に対し $M^a{}_c=\sum_{b,d}T^{abd}\epsilon_{bdc}$ と置く。$a,b$ が対称なので $\tr M=0$ である。$M=0$ と $b,d$ の対称性は同値であり、$a,b$ の対称性と合わせれば全対称になる。従ってこの写像の核は次元 $10$ の $\operatorname{Sym}^3V$ であり、像は次元 $18-10=8$ のトレース零部分を尽くす。したがって $3\otimes3\otimes3=10\oplus8\oplus8\oplus1$ になり、次元 $27$ が回収される。二つの八次元成分に同じ固有値が割り当てられるかは群だけでは決まらない。群と可換な作用素でも、重複度二の空間の上では任意の $2\times2$ 行列を持てるからである。

表現 $V$ と双対 $V^*$ に中心元 $\omega I$、$\omega=e^{2\pi i/3}$ を作用させると、それぞれ $\omega$ と $\omega^{-1}$ が掛かる。したがって $p$ 個の上付き添字と $q$ 個の下付き添字を持つテンソルでは、中心の位相は $\omega^{p-q}$ である。縮約で添字が一つずつ消えても、この位相は変わらない。例えば $3\otimes3$ の $6$ と $\bar3$ にはいずれも $\omega^2$ が掛かり、三重積のすべての成分には $1$ が掛かる。この性質から、中心の作用が非自明な積には群不変なスカラー成分が現れないと分かる。

## Lorentz 群とその有限次元表現

### 時空の符号を保つ変換

計量符号を $(-,+,+,+)$ とし、$x^0=ct$ を他の座標と同じ長さの単位で測る。$X=x^0I+\sum x^i\sigma_i$ は Hermitian で、$\det X=(x^0)^2-|\vb x|^2$ である。$A\in\SL(2,\C)$ に対し $X\mapsto AXA^\dagger$ はこの行列式を保つ。正定値な $X$ は正定値なままなので、未来向き時間的なベクトルを保つ。$\SL(2,\C)$ は極分解 $A=UP$ で $P$ が正定値、$\det P=1$、$U\in\SU(2)$ であり、$P^s$ と $U$ の経路で単位元につながる。核では $X=I$ を代入すると $A$ がユニタリであり、前の核の計算で $\pm I$ に限られる。$A=e^{\eta\sigma_3/2}$ は $x^0,x^3$ を $\cosh\eta,\sinh\eta$ の二次行列で混ぜる。任意の未来向き単位時間ベクトルはこのブーストを回転で挟んだ操作で得られ、元の時間軸を固定する Lorentz 変換は空間回転である。従ってこの写像は空間の向きと時間の向きを保つ連結成分 $\mathrm{SO}^+(1,3)$ への二重被覆である。

この節の $J_i,K_i$ は無次元とする。交換関係は $[J_i,J_j]=i\epsilon_{ijk}J_k$、$[J_i,K_j]=i\epsilon_{ijk}K_k$、$[K_i,K_j]=-i\epsilon_{ijk}J_k$ である。回転に $\sigma_i/2$、ブーストに $i\sigma_i/2$ を用い、実 Lie 群の指数を $e^{-i\theta J}$、$e^{-i\eta K}$ と書けば Pauli 行列の積で確認できる。複素化すると $A_i=(J_i-iK_i)/2$ と $B_i=(J_i+iK_i)/2$ は互いに交換し、それぞれ回転の括弧を満たす。実際、符号の異なる組の交換子では $[J,J]+[K,K]$ が零で、混合項も相殺する。同じ符号の組では $[J,J]-[K,K]$ と二つの混合項が残り、元と同じ符号の結合を与える。

有限次元複素既約型は :ref[prop-lie-lorentz-finite-reps] に従い、二つのスピン $(j_L,j_R)$ で表す。二つの因子をそれぞれ既約に選んだテンソル積なので、次元は $(2j_L+1)(2j_R+1)$ である。回転の生成子は $J_i=A_i+B_i$ だから、回転部分の中心元 $-I$ は二因子でそれぞれ $(-1)^{2j_L},(-1)^{2j_R}$ を掛ける。従って $j_L+j_R$ が整数のとき Lorentz 群へ降りる。

二つの Weyl スピノルは $(1/2,0),(0,1/2)$、四元ベクトルは $(1/2,1/2)$ である。上の生成子の規約では、行列 $A\in\SL(2,\C)$ がそのまま作用する二成分表現において $J+iK=0$ なので、これは $(1/2,0)$ である。複素共役表現ではもう一方だけが非零になる。四元ベクトルを表す行列 $X$ の二添字にはこれら二つが作用する。Dirac スピノルは二つの Weyl 表現の直和である。空間反転は回転を変えずブーストの符号を反転させるので、複素化した二つの因子を交換する。従って単独の Weyl 表現に空間反転を線形に加えることはできず、両方を含む直和なら、それらを交換する作用を定められる。

非自明な連結 Lorentz 群の有限次元ユニタリ表現は存在しない。実際、もしユニタリなら $J_i,K_i$ は Hermitian である。交換関係とトレースの巡回性より
$$
-i\tr(J_3^2)
=\tr([K_1,K_2]J_3)
=\tr(K_1[K_2,J_3])
=i\tr(K_1^2)
$$
となる。二つの二乗のトレースはともに非負なので $J_3=K_1=0$ である。添字を巡回して全生成子が零になり、連結性から表現は自明になる。例えば Weyl 表現ではブースト $e^{\eta\sigma_3/2}$ の固有値が $e^{\pm\eta/2}$ であり、$\eta\ne0$ では正定値内積をどう選んでもユニタリにならない。

### 二形式の Lorentz 表現

四次元の時空二形式には $6$ 成分がある。計量符号 $(-,+,+,+)$ の Hodge 双対は二形式に二回作用すると $-1$ になるので、複素化した二形式は $\star$ の固有値 $+i,-i$ によって三成分ずつに分かれる。向きと計量を保つ Lorentz 変換は $\star$ と可換なので、二つの空間はそれぞれ不変である。実際、四元ベクトル $V_{1/2}\otimes V_{1/2}$ の外積は、一方の二添字を対称、もう一方を反対称とする二通りからなる。二つの交換符号の積が負になるのはその二通りだけであり、次元 $3+3$ で全てを尽くす。式 :eqref[eq-mp-spin-addition] により対称部分はスピン $1$、反対称部分はスピン $0$ だから型は $(1,0)\oplus(0,1)$ である。実の電磁場では複素共役で両空間が交換され、二つの複素三成分を独立に選べるわけではない。

### ブーストと光円錐の保存

Minkowski の符号を $(-,+,+,+)$ とする。$z$ 方向のブーストは $(x^0,z)$ の二成分に $\cosh\eta,\sinh\eta$ を並べた行列で作用する。光円錐成分 $x^0+z,x^0-z$ には $e^\eta,e^{-\eta}$ が掛かる。二つの積は変わらないので、横成分を固定したまま計量の二次形式を保つ。二つのブーストを同じ方向に行うと指数の引数が足される。

:::example[光円錐の伸縮で速さを求める（1）]
ラピディティを $\eta=\log(m+1)$ とする。

$$
m=1
$$

光円錐の二方向の倍率を求める。

$$
\begin{aligned}
(e^\eta,e^{-\eta})\big|_{m=1} &= \left( 2, \  \frac{1}{2}\right)
\end{aligned}
$$

時空成分での行列を二つの倍率の和と差で作る。

$$
\begin{aligned}
\Lambda_{0z}\big|_{m=1} &= \begin{pmatrix}\frac{5}{4} & \frac{3}{4}\\\frac{3}{4} & \frac{5}{4}\end{pmatrix}
\end{aligned}
$$

物理的な速さは $c$ を掛けた $\tanh\eta$ である。

$$
\begin{aligned}
v/c\big|_{m=1} &= \frac{3}{5}
\end{aligned}
$$
:::

:::example[光円錐の伸縮で速さを求める（2）]
ラピディティを $\eta=\log(m+1)$ とする。

$$
m=2
$$

光円錐の二方向の倍率を求める。

$$
\begin{aligned}
(e^\eta,e^{-\eta})\big|_{m=2} &= \left( 3, \  \frac{1}{3}\right)
\end{aligned}
$$

時空成分での行列を二つの倍率の和と差で作る。

$$
\begin{aligned}
\Lambda_{0z}\big|_{m=2} &= \begin{pmatrix}\frac{5}{3} & \frac{4}{3}\\\frac{4}{3} & \frac{5}{3}\end{pmatrix}
\end{aligned}
$$

物理的な速さは $c$ を掛けた $\tanh\eta$ である。

$$
\begin{aligned}
v/c\big|_{m=2} &= \frac{4}{5}
\end{aligned}
$$
:::

:::example[光円錐の伸縮で速さを求める（3）]
ラピディティを $\eta=\log(m+1)$ とする。

$$
m=3
$$

光円錐の二方向の倍率を求める。

$$
\begin{aligned}
(e^\eta,e^{-\eta})\big|_{m=3} &= \left( 4, \  \frac{1}{4}\right)
\end{aligned}
$$

時空成分での行列を二つの倍率の和と差で作る。

$$
\begin{aligned}
\Lambda_{0z}\big|_{m=3} &= \begin{pmatrix}\frac{17}{8} & \frac{15}{8}\\\frac{15}{8} & \frac{17}{8}\end{pmatrix}
\end{aligned}
$$

物理的な速さは $c$ を掛けた $\tanh\eta$ である。

$$
\begin{aligned}
v/c\big|_{m=3} &= \frac{15}{17}
\end{aligned}
$$
:::

:::example[光円錐の伸縮で速さを求める（4）]
ラピディティを $\eta=\log(m+1)$ とする。

$$
m=4
$$

光円錐の二方向の倍率を求める。

$$
\begin{aligned}
(e^\eta,e^{-\eta})\big|_{m=4} &= \left( 5, \  \frac{1}{5}\right)
\end{aligned}
$$

時空成分での行列を二つの倍率の和と差で作る。

$$
\begin{aligned}
\Lambda_{0z}\big|_{m=4} &= \begin{pmatrix}\frac{13}{5} & \frac{12}{5}\\\frac{12}{5} & \frac{13}{5}\end{pmatrix}
\end{aligned}
$$

物理的な速さは $c$ を掛けた $\tanh\eta$ である。

$$
\begin{aligned}
v/c\big|_{m=4} &= \frac{12}{13}
\end{aligned}
$$
:::

### 実の二形式の複素分解

四次元の向きと Lorentz 計量を固定すると、二形式の Hodge 双対は $*^2=-1$ を満たす。この符号は直交した正規共基底で二度作用させると、向きの交換の符号に時間成分の負号が一つ加わることから出る。複素化した二形式に対し $P_\pm=(1\mp i*)/2$ と置けば $*P_\pm=\pm i P_\pm$、$P_++P_-=1$、$P_+P_-=0$ となる。向きと計量を保つ変換は $*$ と可換なので、それぞれの像が表現空間になる。実二形式では二空間に入る成分が複素共役で結ばれる。

例えば $*\alpha=\beta$、$*\beta=-\alpha$ の基底を一対取る。すると $\alpha-i\beta$ には $+i$、$\alpha+i\beta$ には $-i$ が掛かる。平行な電気成分と磁気成分が $\alpha,\beta$ の係数 $E,B$ であるなら、複素化した係数は $(E+iB)/2$ と $(E-iB)/2$ になる。$E=2,B=1$ の場合も足し戻して二つの実係数が回収される。電気と磁気を SI で比較する場合には、時間座標を $ct$ とするなどして二形式の成分の単位を合わせてから係数を比べる。


:::example[時空平面での標準ブースト]
光円錐成分の倍率を $3$ と $1/3$ にする。行ベクトルは光円錐成分のうち時間と空間の和を選び出す。

座標成分を列に並べる。この順序に合わせ、線形写像の行列の列を選ぶ。

$$
\begin{aligned}v&=\begin{pmatrix}2\\1\end{pmatrix}\\ A&=\begin{pmatrix}\frac{5}{3} & \frac{4}{3}\\\frac{4}{3} & \frac{5}{3}\end{pmatrix}\end{aligned}
$$

行列の各行と入力の内積を計算する。例えば第一行では $\begin{pmatrix}\frac{5}{3} & \frac{4}{3}\end{pmatrix}$ と入力列を掛ける。

$$
\begin{aligned}(Av)_1&=\frac{14}{3}\\Av&=\begin{pmatrix}\frac{14}{3}\\\frac{13}{3}\end{pmatrix}\end{aligned}
$$

核の次元は列数から階数を引いて得られる。像の次元は階数である。

$$
\begin{aligned}\operatorname{rank}A&=2\\\dim\ker A&=0\end{aligned}
$$

正方行列の場合は、特性多項式を求めて固有値を調べられる。ここでは $tI-A$ の行列式を取る。

$$
\begin{aligned}\operatorname{tr}A&=\frac{10}{3}\\ \det A&=1\\ \det(tI-A)&=t^{2} - \frac{10 t}{3} + 1\end{aligned}
$$

行列要素を得るため、出力から指定した共変ベクトルで成分を取り出す。

$$
\begin{aligned}\ell&=\begin{pmatrix}1 & 1\end{pmatrix}\\ \ell Av&=9\end{aligned}
$$

次に別の方向を入力にする。出力の差は入力の差に $A$ を掛けたものである。

$$
\begin{aligned}w&=\begin{pmatrix}1\\-1\end{pmatrix}\\ Aw&=\begin{pmatrix}\frac{1}{3}\\- \frac{1}{3}\end{pmatrix}\\ Aw-Av&=\begin{pmatrix}- \frac{13}{3}\\- \frac{14}{3}\end{pmatrix}\end{aligned}
$$
:::

## 演習問題

::::exercise{level=basic}
$\SU(2)$ の $z$ 回転に $\theta=\pi$ を入れよ。
:::solution
$U=\diag(e^{-i\pi/2},e^{i\pi/2})=\diag(-i,i)$ である。共役作用で $x,y$ の符号が反転し、$z$ は不変である。
:::
::::

::::exercise{level=basic}
$j=3/2$ の $m$ と次元を求めよ。
:::solution
昇降の鎖は $-3/2,-1/2,1/2,3/2$ の四つを持つ。従って次元は $4$ である。
:::
::::

::::exercise{level=basic}
$1\otimes1/2$ のスピン分解を求めよ。
:::solution
最大は $3/2$、最小は $1/2$ なので $V_{3/2}\oplus V_{1/2}$ である。次元は $4+2=6$ で $3\cdot2$ と一致する。
:::
::::

::::exercise{level=basic}
$\SU(3)$ の $(2,0)$ の次元を求めよ。
:::solution
基本表現の対称二階テンソルなので三つの対角と三つの非対角を持ち、次元 $6$ である。次元式でも $3\cdot1\cdot4/2=6$ である。
:::
::::

::::exercise{level=advanced}
二つのスピン $1/2$ の一重項が全 $L_z,L_\pm$ で零になることを直接確認せよ。
:::solution
二つの積状態とも全 $m=0$ で $L_z$ は零である。$L_+$ は両方を $|\uparrow\uparrow\rangle$ に、$L_-$ は両方を $|\downarrow\downarrow\rangle$ に同じ係数で移す。差でキャンセルする。
:::
::::

::::exercise{level=advanced}
$j=2,m=2$ から階数 $k=1$ の $q=1$ 成分が許す終状態を判定せよ。
:::solution
式 :eqref[eq-mp-wigner-eckart] から $m'=3$、三角条件では $j'=1,2,3$ である。$|m'|\le j'$ を満たすのは $j'=3,m'=3$ だけである。
:::
::::

::::exercise{level=advanced}
$\SU(3)$ で $(2,1)$ の次元を縮約から求めよ。
:::solution
対称二階の $V$ と一個の $V^*$ の積は $6\cdot3=18$ 次元である。縮約は一個の $V$ への全射なので三次元を引き、$15$ となる。
:::
::::

::::exercise{level=advanced}
$A=e^{\eta\sigma_3/2}$ の $X=x^0I+x^3\sigma_3$ への作用を計算せよ。
:::solution
対角成分 $x^0+x^3$、$x^0-x^3$ はそれぞれ $e^\eta,e^{-\eta}$ 倍される。半和と半差を取ると $x'^0=x^0\cosh\eta+x^3\sinh\eta$、$x'^3=x^0\sinh\eta+x^3\cosh\eta$ である。二乗差は $\cosh^2\eta-\sinh^2\eta=1$ により不変である。
:::
::::

## 確認テスト

::::quiz{answer="5" tol="1e-9"}
$j=2$ の回転既約表現の次元を答えよ。

:::explanation
$m=-2,-1,0,1,2$ で $2j+1=5$ である。
:::
::::

::::quiz
$\SO(3)$ の通常の表現に降りるのはどのスピンか。

- [x] 整数スピン
- [ ] 半整数で整数でないスピン
- [ ] すべてのスピン

:::explanation
核の $-I$ が $(-1)^{2j}$ と作用する。これが $+1$ のときのみ被覆の二点を同一視できる。
:::
::::

::::quiz{answer="8" tol="1e-9"}
$\SU(3)$ の随伴表現の次元を答えよ。

:::explanation
複素 $3\times3$ 行列の九成分からトレース一条件を引いて $8$ である。
:::
::::

::::quiz
時間成分 $x^0=ct$ を含む Lorentz 四元ベクトルの型はどれか。

- [ ] $(1/2,0)$
- [x] $(1/2,1/2)$
- [ ] $(1,0)$

:::explanation
$X\mapsto AXA^\dagger$ で一つの行列添字が基本、他方が共役の二成分表現で変わる。次元も $2\cdot2=4$ である。
:::
::::

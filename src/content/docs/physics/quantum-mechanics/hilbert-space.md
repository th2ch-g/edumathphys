---
title: Hilbert 空間と Dirac 記法
description: 重ね合わせを内積空間で表し、固有値による分解を証明してから、位置と運動量を連続変数で表す。
notation:
  - symbol: '\mathcal H'
    meaning: 複素 Hilbert 空間
  - symbol: 'A^\dagger'
    meaning: 随伴演算子
  - symbol: 'P_a'
    meaning: 固有値 a に属する直交射影
---

## この章の目標

- 内積の線形性と共役を使い分ける。
- 有限次元のスペクトル分解と同時対角化を証明する。
- 離散基底の和と一般化固有状態の積分を区別する。

## 状態空間としての Hilbert 空間

二つの経路を通った粒子を同じ検出器で数えると、各経路を単独で開いた場合の計数の和からずれることがある。
複素数の振幅なら、和の絶対値の二乗に干渉項が現れる。
各検出器に届く振幅を成分として並べると、重ね合わせを複素ベクトルの和で表せる。

:::definition[状態の内積]{#def-qm-hilbert}
複素ベクトル空間に、共役対称性 $\braket{\phi}{\psi}^*=\braket{\psi}{\phi}$、正定値性 $\braket{\psi}{\psi}>0$（$\psi\ne0$）、および第2引数についての線形性を満たす内積を定める。
ノルム $\norm{\psi}=\sqrt{\braket{\psi}{\psi}}$ に関して完備なものを :term[Hilbert 空間]{en="Hilbert space" yomi="ひるべるとくうかん"} という。
完備とは、任意の Cauchy 列が空間内の元へノルム収束することである。
:::

この科目の内積は**第2引数（ケット）について線形、第1引数（ブラ）について共役線形**である。
数学科目の第1引数を線形とする規約とは逆である。
例えば $\braket{i\phi}{\psi}=-i\braket{\phi}{\psi}$ となる。

:::example[二つの振幅]
$\C^2$ では $\braket{u}{v}=\sum_j u_j^*v_j$ とする。
$u=(1,i)^\top/\sqrt2$ と $v=(1,-i)^\top/\sqrt2$ はいずれも単位ベクトルであり、内積は $(1+(-i)^2)/2=0$ となる。
:::

:::example[極限としての波束]
$L^2(\R,\dd x)$ の内積は $\int\phi^*(x)\psi(x)\dd x$ である。
測度零の集合上だけ異なる関数を同じ元とする。
区間 $[0,L]$ で $1/\sqrt L$、外で零の波束は規格化される。
端点の値をどのように決めてもノルムは変わらない。
:::

### 直交成分と不等式

$u\ne0$ に対し $v_\perp=v-u\braket{u}{v}/\norm{u}^2$ と置くと $\braket{u}{v_\perp}=0$ である。
直交する二項のノルム二乗は加法的なので

$$
0\le\norm{v_\perp}^2=\norm v^2-\frac{|\braket{u}{v}|^2}{\norm u^2}.
$$

従って $|\braket{u}{v}|\le\norm u\,\norm v$ であり、$u\ne0$ の場合の等号は $v_\perp=0$ と同値である。
この Cauchy–Schwarz 不等式は測定の分散を比較するときに用いる。

### 有限個の成分から無限列へ

可算個のモードに振幅を割り振るときの模型は $\ell^2$ である。
例えば $c_n=\sqrt3\,2^{-n}$（$n=1,2,\ldots$）なら $\sum_n|c_n|^2=1$。
第 $N$ 項までの部分和状態と元の状態の距離の二乗は $4^{-N}$ である。
各成分が収束することとノルム収束は異なる。
列 $e_n=(0,\ldots,0,1,0,\ldots)$ で成分を一つ固定すればやがて零になるが、$\norm{e_n-e_m}=\sqrt2$（$n\ne m$）なので Cauchy 列ではない。

:::example[測定を有限個の基底に限る]
任意の正規直交族に $P_N=\sum_{j=1}^N\ketbra j j$ を対応させると $\norm\psi^2=\norm{P_N\psi}^2+\norm{(I-P_N)\psi}^2$。
これが $\sum_{j=1}^N|\braket j\psi|^2\le\norm\psi^2$ という Bessel の不等式である。
測っていない方向への残りの確率を右の差として読み取れる。
:::

### 完備化を実験の近似と結び付ける

一つの波形を有限個のモードの和として記録する装置を考える。読めるモードの個数を増やすたびに別の状態を準備するのでは、同じ波形を近似していることにはならない。第 $N$ 段階の近似 $\psi_N$ が任意の精度で互いに近づき、その極限が状態として扱えることが必要になる。例えば有限個の成分だけが非零の複素数列は内積をもつが、それだけでは Hilbert 空間にならない。$v_N=\sum_{j=1}^N2^{-j}\ket j$ と置くと、$M>N$ に対し

$$
\|v_M-v_N\|^2=\sum_{j=N+1}^M4^{-j}
=\frac{4^{-N}}3(1-4^{-(M-N)})\le\frac{4^{-N}}3 .
$$

従ってこの列は Cauchy 列である。もし有限個しか非零成分をもたない $v$ がその極限なら、各成分に対して $|\braket j{v_N-v}|\le\|v_N-v\|\to0$ より $\braket j v=2^{-j}$ でなければならない。どの $j$ でも零にならないので矛盾する。この欠落を補うのが $\ell^2$ への完備化である。

強度と位相の違いを幾何学的に見るには、ノルムだけから内積を復元する方法がある。差 $D(u,v)=\|u+v\|^2-\|u-v\|^2$ は $4\Re\braket u v$ に等しい。また $\braket u{iv}=i\braket u v$ なので $D(u,iv)=-4\Im\braket u v$。二つを合わせれば

$$
\braket u v=\tfrac14\{D(u,v)-iD(u,iv)\}.
$$

第 1 引数で線形な教科書の公式をそのまま使うと虚部だけが逆になる。干渉の強度を測るときには、基準の経路に $+\pi/2$ の位相を足したか、信号の経路に足したかを書いておくと符号を照合できる。


### 閉部分空間への最良近似

有限個のモードだけを残す場合には、直交射影を外積の有限和で計算できる。一般の閉部分空間 $M$ への最良近似も存在する。これを見るため、$d=\inf_{w\in M}\|v-w\|$ と置き、$\|v-w_n\|\to d$ となる列をとる。平行四辺形の等式と $(w_n+w_k)/2\in M$ より

$$
\|w_n-w_k\|^2
=2\|v-w_n\|^2+2\|v-w_k\|^2-4\|v-(w_n+w_k)/2\|^2
\le2\|v-w_n\|^2+2\|v-w_k\|^2-4d^2.
$$

従って $(w_n)$ は Cauchy 列であり、完備性と閉性から $w\in M$ に収束する。このとき $\|v-w\|=d$ である。任意の $u\in M$ に対し $\|v-w-tu\|^2$ は実数 $t=0$ で最小となるから、$\Re\braket{u}{v-w}=0$。$iu$ にも同じことを適用すると虚部も零となる。したがって $v-w\perp M$ である。この直交分解は、二通りの分解の差が $M$ と $M^\perp$ の両方に属することから一意である。よって $P_Mv=w$ と置けば、線形性は分解の一意性から、$\|P_Mv\|\le\|v\|$ は直交性から従う。

### 有界な読みと状態の近似誤差

同じ状態を二つの計算で近似し、その規格化した結果を $\psi,\phi$ とする。装置の読みを表す演算子 $A$ が有界であるとは、どんな入力 $v$ に対しても $\|Av\|\le M\|v\|$ を満たす有限な定数 $M$ があることをいう。このような $M$ の下限を $\|A\|$ とする。期待値の差に中間項を一つ加えて引くと

$$
\mel\psi A\psi-\mel\phi A\phi=
\langle\psi-\phi|A\psi\rangle+\langle\phi|A(\psi-\phi)\rangle .
$$

各項に Cauchy–Schwarz を使えば、差の絶対値は $2\|A\|\|\psi-\phi\|$ 以下である。従ってノルム近似の精度から、有界な読みのずれを一様に保証できる。射影では $\|A\|\le1$ なので、確率だけを問う場合に特に便利である。

全直線での位置のように有界でない演算子では事情が違う。離散模型で $A\ket n=n^2\ket n$ とし、$\psi_N=\sqrt{1-N^{-2}}\ket0+N^{-1}\ket N$（$N\ge2$）を用意する。$\|\psi_N-\ket0\|^2=2-2\sqrt{1-N^{-2}}\to0$ だが $\mel{\psi_N}A{\psi_N}=1$ であり、$\mel0A0=0$ に収束しない。滅多に出ない高い測定値が平均に有限の寄与を残す例である。エネルギーの近似を評価するとき、規格化の誤差だけを見る危険が分かる。


### 有限分解能で波形の最良近似を作る

幅 $2L$ の検出面を横切る一成分の波形を考える。位置を $L=2.0\,\mathrm{mm}$ で割ると $-1<x<1$。装置は、明るさが一定のモード、左右で符号が異なるモード、さらに二次・三次の曲がりを持つモードを分離できる。候補 $1,x,x^2,x^3$ は独立だが互いには直交しない。

積分変数を $x$ とし、この計算の内積を $\langle f|g\rangle=\int_{-1}^{1}f^*(x)g(x)w(x)\,\dd x$ と書く。重みは $w(x)=1$ である。順に次の候補を入れて独立な方向に分ける。

$$
(v_0,\ldots,v_{3})=\left(1,\;x,\;x^{2},\;x^{3}\right).
$$

**積分で基底を作る。** 第 $k$ 段では、既に長さ一にした $u_j$ に沿う成分を候補から引き去る。引いた後の残りを $r_k$ とすると

$$
r_k=v_k-\sum_{j<k}u_j\braket{u_j}{v_k},\qquad u_k=\frac{r_k}{\sqrt{\braket{r_k}{r_k}}}.
$$



$$
\begin{aligned}r_0&=1,\qquad \|r_0\|^2=2,\\ u_0(x)&=\frac{\sqrt{2}}{2}.\end{aligned}
$$



$$
\begin{aligned}\braket{u_0}{v_1}&=0,\\r_1&=x,\qquad \|r_1\|^2=\frac{2}{3},\\ u_1(x)&=\frac{\sqrt{6} x}{2}.\end{aligned}
$$



$$
\begin{aligned}(\braket{u_0}{v_2},\ldots,\braket{u_1}{v_2})&=\left(\frac{\sqrt{2}}{3},\,0\right),\\r_2&=\frac{3 x^{2} - 1}{3},\qquad \|r_2\|^2=\frac{8}{45},\\ u_2(x)&=\frac{\sqrt{10} \left(3 x^{2} - 1\right)}{4}.\end{aligned}
$$



$$
\begin{aligned}(\braket{u_0}{v_3},\ldots,\braket{u_2}{v_3})&=\left(0,\,\frac{\sqrt{6}}{5},\,0\right),\\r_3&=\frac{x \left(5 x^{2} - 3\right)}{5},\qquad \|r_3\|^2=\frac{8}{175},\\ u_3(x)&=\frac{\sqrt{14} x \left(5 x^{2} - 3\right)}{4}.\end{aligned}
$$

分母に現れた正数は、各段で本当に新しい方向が加わったことを保証する。得られた各関数に演算子 $D$ を作用させる。作用の結果に対して左から順に内積を取ると、列ごとに行列が埋まる。 ここで $D$ は位置を $L$ で割った演算子であり、$bD$ の実際の単位はメートルである。

| 列 $k$ | $D u_k(x)$ | 各基底への射影係数 |
|---|---|---|
| $0$ | $\frac{\sqrt{2} x}{2}$ | $\begin{pmatrix}0\\\frac{\sqrt{3}}{3}\\0\\0\end{pmatrix}$ |
| $1$ | $\frac{\sqrt{6} x^{2}}{2}$ | $\begin{pmatrix}\frac{\sqrt{3}}{3}\\0\\\frac{2 \sqrt{15}}{15}\\0\end{pmatrix}$ |
| $2$ | $\frac{\sqrt{10} x \left(3 x^{2} - 1\right)}{4}$ | $\begin{pmatrix}0\\\frac{2 \sqrt{15}}{15}\\0\\\frac{3 \sqrt{35}}{35}\end{pmatrix}$ |
| $3$ | $\frac{\sqrt{14} x^{2} \left(5 x^{2} - 3\right)}{4}$ | $\begin{pmatrix}0\\0\\\frac{3 \sqrt{35}}{35}\\0\end{pmatrix}$ |



$$
D_{\rm sub}=(\mel{u_j}D{u_k})_{jk}=\begin{pmatrix}0 & \frac{\sqrt{3}}{3} & 0 & 0\\\frac{\sqrt{3}}{3} & 0 & \frac{2 \sqrt{15}}{15} & 0\\0 & \frac{2 \sqrt{15}}{15} & 0 & \frac{3 \sqrt{35}}{35}\\0 & 0 & \frac{3 \sqrt{35}}{35} & 0\end{pmatrix},\qquad b=0.002\,\mathrm m.
$$

次に、有限個の成分を残したことによる影響を確認する。$D u_k$ から表に含まれる部分を引き、その残りのノルムを調べれば、計算空間の外側へ出る成分を見落とさずに済む。

| 入力 $u_k$ | 射影で落ちる関数 | 落ちる成分のノルムの二乗 |
|---|---|---|
| 0 | $0$ | $0$ |
| 1 | $0$ | $0$ |
| 2 | $0$ | $0$ |
| 3 | $\frac{\sqrt{14} \left(35 x^{4} - 30 x^{2} + 3\right)}{28}$ | $\frac{16}{63}$ |

この基底を使って別の準備 $g(x)=\frac{3 \sqrt{170} \left(x^{4} + 1\right)}{68}$ を調べる。成分は各 $u_j$ と $g$ の内積なので

$$
c=\begin{pmatrix}\frac{9 \sqrt{85}}{85}\\0\\\frac{6 \sqrt{17}}{119}\\0\end{pmatrix},\qquad \sum_j|c_j|^2=\frac{4149}{4165},\qquad \|g-\sum_jc_ju_j\|^2=\frac{16}{4165}.
$$

最後の等式は直交性による。落ちた成分がある場合、残った列を長さ一に直した量は、この部分空間内にいるという条件を付けて作った状態である。

任意の近似 $v$ を同じ部分空間に選ぶと $g-v=(g-Pg)+(Pg-v)$。両項の内積が零なので誤差の二乗は和になり、$v=Pg$ が一意の最良近似である。


## ブラとケット

:term[ケット]{en="ket" yomi="けっと"} $\ket\psi$ はベクトルそのものを表す。
対応する :term[ブラ]{en="bra" yomi="ぶら"} $\bra\psi$ は $\ket v\mapsto\braket{\psi}{v}$ という連続線形汎関数である。
対応 $\ket\psi\mapsto\bra\psi$ は共役線形であり、基底で列をブラに直す操作は共役転置である。

有限次元では :ref[thm-linalg-riesz-finite] に対応して、任意の線形汎関数 $f$ はこの形である。
正規直交基底 $\ket j$ で $u=\sum_j f(\ket j)^*\ket j$ と置けば、展開の線形性から $\braket u v=f(v)$ が従う。
二つの代表の差が全ベクトルに直交すれば、その差自身のノルムも零なので代表は一意である。

無限次元でも、連続線形汎関数は同じ形で表せる。$f=0$ なら代表は零である。$f\ne0$ なら核 $M=\{v:f(v)=0\}$ は閉部分空間であり、$f(v_0)\ne0$ となる $v_0$ がある。$u=v_0-P_Mv_0$ とすれば $u\perp M$ かつ $f(u)=f(v_0)\ne0$ である。任意の $v$ に対し $v-u f(v)/f(u)\in M$ なので

$$
f(v)=\frac{f(u)}{\|u\|^2}\braket{u}{v}
=\left\langle \frac{f(u)^*}{\|u\|^2}u\,\middle|\,v\right\rangle .
$$

この代表も有限次元の場合と同じ理由で一意である。これが第 2 引数を線形とした Riesz の表現定理である。連続性を使ったのは核の閉性を保証するためである。

### 基底と恒等演算子

:ref[def-linalg-basis] の基底を直交化するには、順に $w_k=v_k-\sum_{j<k}\ket{e_j}\braket{e_j}{v_k}$ を作り、$e_k=w_k/\norm{w_k}$ とする。
引いた各成分との内積が零となり、元の基底の独立性から $w_k\ne0$ である。
従って途中までの張る空間を保ったまま正規直交基底を得る。

$$
\begin{aligned}
\ket\psi&=\sum_j c_j\ket j, \\
c_j&=\braket j\psi, \\
I&=\sum_j\ketbra j j, \\
\norm{\psi}^2&=\sum_j|c_j|^2.
\end{aligned}
$$

最後の式は内積を展開して直交性を使えば出る。
無限個の離散基底では、有限和への射影がノルムで $\psi$ に近づくことが完全性である。
級数に任意個の項を許す点で代数的な基底と異なる。

逆に、正規直交族の全成分と直交するベクトルが零だけなら、ノルム収束が従う。Bessel の不等式から $\sum_j |c_j|^2\le\|\psi\|^2$ であり、部分和どうしの距離の二乗は対応する係数の二乗和である。従って部分和は Cauchy 列で、完備性によりある $v$ に収束する。内積の連続性から $\braket{j}{\psi-v}=0$ となるので、仮定より $v=\psi$。このとき $\sum_j\ketbra{j}{j}=I$ は各ベクトルへの作用がノルム収束する等式である。

:::example[基底を回す]
$\ket\pm=(\ket0\pm\ket1)/\sqrt2$ で $\ket\psi=(\ket0+i\ket1)/\sqrt2$ を展開すると係数は $(1\pm i)/2$ である。
それぞれの絶対値二乗は $1/2$ で総和は $1$ になる。
:::

:::example[外積の作用]
$\ketbra0{1}$ は $a\ket0+b\ket1$ を $b\ket0$ に写す。
左右の順を逆にした $\ketbra1{0}$ は $a\ket1$ を返す。
従って内積が数を与えるのに対し、外積は向きを持つ写像である。
:::

### 受動的な基底変換

実験装置も状態も固定して、記述に使う正規直交基底だけを $\ket{e_j}$ から $\ket{f_k}$ に替える。
$S_{jk}=\braket{e_j}{f_k}$ とすれば、列 $k$ は新基底の古い座標である。
完全性から $(S^\dagger S)_{kl}=\sum_j\braket{f_k}{e_j}\braket{e_j}{f_l}=\delta_{kl}$。
新しい座標は $c'=S^\dagger c$、行列は $A'=S^\dagger AS$ となるので $c'^\dagger A'c'=c^\dagger Ac$ である。
状態を実際に回す演算と区別するには、どのケットの成分を並べたかを明示するとよい。
数学の基底変換は :ref[thm-linalg-change-of-basis] に対応する。

:::example[位相を含む基底]
$f_1=(e_1+ie_2)/\sqrt2,\ f_2=(ie_1+e_2)/\sqrt2$ とする。
元の $e_1$ の新座標は $(1,-i)^\top/\sqrt2$。
二列にこの係数を掛けて足すと $e_2$ 成分は消え、確かに元のベクトルに戻る。
:::

### テンソル積の基底と一成分の添字

どちらの経路を通ったかと、どちらの偏光で通ったかを同じ光子について記述する。経路のケットを $\ket{a},\ket b$、偏光のケットを $\ket H,\ket V$ とし、それぞれ正規直交すると仮定する。一つずつ指定した組は $\ket{aH}=\ket a\otimes\ket H$ のように書く。二つの振幅が独立に掛かるため、積ケットの内積は $\braket{aH}{bV}=\braket a b\braket H V$ である。四つの組は正規直交し、それらへの射影の和は $(P_a+P_b)\otimes(P_H+P_V)=I$ になる。

一般の状態は四つの係数 $c_{ij}$ の和であり、いつでも一つの積に戻せるわけではない。積なら $c_{ij}=u_i v_j$ で行列の階数が一になる。逆に係数行列の階数が一なら、非零の列を $u$ と選び、ほかの列との比例係数を $v_j$ とすれば積に分けられる。特に $\ket\Psi=(\ket{aH}+i\ket{bV})/\sqrt2$ の係数行列の行列式は $i/2$ で零でない。経路だけの装置 $P_a\otimes I$ を作用させると残るのは $\ket{aH}/\sqrt2$ である。偏光が二通りあるという理由で確率を二倍することはできない。


### 非直交な二つの偏光を座標軸にする

偏光の準備装置から水平偏光と円偏光を取り出せるが、垂直偏光用の出力がない場合を考える。二つの出力は一次独立なので状態の展開には使える。しかし各出力への射影は直交する二結果の測定を構成しない。

直交する基準の座標系で、各候補状態を列に並べた行列を $C$ とする。問い合わせたいケットは右に示した $v$ である。

$$
C=\begin{pmatrix}1 & \frac{\sqrt{2}}{2}\\0 & \frac{\sqrt{2} i}{2}\end{pmatrix},\qquad v=\begin{pmatrix}\frac{\sqrt{2}}{2}\\- \frac{\sqrt{2}}{2}\end{pmatrix}.
$$

候補同士の内積をまとめた :term[グラム行列]{en="Gram matrix" yomi="ぐらむぎょうれつ"} と、候補に対する振幅の列を計算する。

$$
G=C^\dagger C=\begin{pmatrix}1 & \frac{\sqrt{2}}{2}\\\frac{\sqrt{2}}{2} & 1\end{pmatrix},\qquad y=C^\dagger v=\begin{pmatrix}\frac{\sqrt{2}}{2}\\\frac{1}{2} + \frac{i}{2}\end{pmatrix}.
$$

このとき $v=Cc$ の両辺に $C^\dagger$ を掛けると $y=Gc$。候補は互いに直交していないため、$y$ は探している展開係数 $c$ と一般には違う。$d^\dagger Gd=\|Cd\|^2$ であり、候補が一次独立なら $d\ne0$ で正である。したがって $G$ は逆行列をもつ。

$$
G^{-1}=\begin{pmatrix}2 & - \sqrt{2}\\- \sqrt{2} & 2\end{pmatrix},\qquad c=G^{-1}y=\begin{pmatrix}\frac{\sqrt{2} \left(1 - i\right)}{2}\\i\end{pmatrix},\qquad Cc=\begin{pmatrix}\frac{\sqrt{2}}{2}\\- \frac{\sqrt{2}}{2}\end{pmatrix}.
$$

この係数を内積一回で与える別の基底を作れる。列を $\tilde C=CG^{-1}$ とすると $\tilde C^\dagger C=I$。従って $c_i=\langle\tilde c_i|v\rangle$、$\sum_i|c_i^{\rm basis}\rangle\langle\tilde c_i|=I$ となる。ここで $|c_i^{\rm basis}\rangle$ は $C$ の第 $i$ 列である。元の候補の外積だけを足した演算子と比較しよう。

$$
\tilde C=\begin{pmatrix}1 & 0\\- i & \sqrt{2} i\end{pmatrix},\qquad CC^\dagger=\begin{pmatrix}\frac{3}{2} & - \frac{i}{2}\\\frac{i}{2} & \frac{1}{2}\end{pmatrix},\qquad C\tilde C^\dagger=\begin{pmatrix}1 & 0\\0 & 1\end{pmatrix}.
$$

次に、直交する座標で物理量が $bM$、$b=1.0546\times10^{-34}\,\mathrm{J\,s}$ と与えられたとする。非直交な座標で使う二つの行列を実際に求める。

$$
M=\begin{pmatrix}0 & - i\\i & 0\end{pmatrix},\qquad D=C^{-1}MC=\begin{pmatrix}-1 & 0\\\sqrt{2} & 1\end{pmatrix},\qquad B=C^\dagger MC=\begin{pmatrix}0 & \frac{\sqrt{2}}{2}\\\frac{\sqrt{2}}{2} & 1\end{pmatrix}.
$$

前者は座標に左から掛けて作用後の座標を返す。後者は候補状態どうしの**行列要素**である。関係 $B=GD$ から、エルミート性の座標表示は $D^\dagger G=GD$ になる。規格化と期待値を実際に戻すと

$$
c^\dagger Gc=1,\qquad c^\dagger Bc=0,\qquad b\,v^\dagger Mv\simeq 0\,\mathrm{J\,s}.
$$

となる。ここで物理量は伝播軸まわりの光子のスピン成分である。水平・垂直の座標で直線偏光は両方のヘリシティを等しい重みで含むため、平均が零になる。


### 斜めの二軸で展開係数を読む双対基底

二つのフィルターが作る規格化した状態を基底に使いたい。二方向が独立なら平面の位置は定まる。しかし両者が直交しないと、一方への射影には他方に沿う係数も入る。展開係数を読み分けるブラを作る。

線形独立な列を $v_1,\ldots,v_d$ として $V=(v_1,\ldots,v_d)$、$G=V^\dagger V$ と置く。任意の非零列 $b$ に対し $b^\dagger Gb=\|Vb\|^2>0$ だから $G$ は可逆である。$\ket\psi=Vc$ に左から $V^\dagger$ を掛けると、測った内積の列は $V^\dagger\psi=Gc$。したがって $c=G^{-1}V^\dagger\psi$。$W=VG^{-1}$ の第 $j$ 列を $w_j$ とすれば、$W^\dagger V=I$ なので


$$
\langle w_i|v_j\rangle=\delta_{ij},\qquad P_{\operatorname{span} V}=\sum_i\ket{v_i}\bra{w_i}=V(V^\dagger V)^{-1}V^\dagger.\label{eq-qm-dual-basis}
$$

右辺の随伴は自分自身であり、二度掛けると中央の $V^\dagger V$ を一つ消せる。さらに各 $v_i$ を変えず、すべての $v_i$ に直交する列を消す。従って基底が全空間を張らない場合も直交射影として使える。個々の $\ket{v_i}\bra{w_i}$ は一般には自己共役ではなく、各項を一つずつ測定確率と解釈することはできない。

```tikz caption="斜交基底。双対ベクトルは相手側の基底を消す向きに取る。"
\begin{tikzpicture}[>=Stealth,x=1cm,y=1cm]
\draw[->] (-0.4,0)--(3.3,0); \draw[->] (0,-1.8)--(0,2.2);
\draw[thick,->] (0,0)--(2.4,0) node[below] {$v_1$};
\draw[thick,->] (0,0)--(1.7,1.7) node[right] {$v_2$};
\draw[dashed,->] (0,0)--(1.4,-1.4) node[right] {$w_1$};
\draw[dashed,->] (0,0)--(0,1.9) node[left] {$w_2$};
\end{tikzpicture}
```


## 線形演算子とエルミート共役

:ref[def-linalg-linear-map] のようにベクトルの和と複素スカラー倍を保つ写像を用いる。
正規直交基底で $A_{ij}=\mel i A j$ と置くと、積の間に $I$ を入れて $(AB)_{ij}=\sum_k A_{ik}B_{kj}$ が従う。

:::definition[随伴と直交射影]{#def-qm-adjoint}
有限次元で $\braket{\phi}{A\psi}=\braket{A^\dagger\phi}{\psi}$ をすべての $\phi,\psi$ について満たす $A^\dagger$ を :term[随伴演算子]{en="adjoint operator" yomi="ずいはんえんざんし"} という。
$A^\dagger=A$ ならエルミート、$U^\dagger U=UU^\dagger=I$ ならユニタリ、$P^2=P=P^\dagger$ なら直交射影である。
:::

基底の各対に定義を適用すると $(A^\dagger)_{ij}=A_{ji}^*$ だから存在と一意性が従う。
また $\braket u{ABv}=\braket{A^\dagger u}{Bv}=\braket{B^\dagger A^\dagger u}v$ より $(AB)^\dagger=B^\dagger A^\dagger$ である。
ユニタリ演算子が内積を保つことも定義から従う。

:::example[エルミートでない積]
$X=\begin{pmatrix}0&1\\1&0\end{pmatrix}$、$Z=\diag(1,-1)$ とする。
両方ともエルミートだが $XZ=-ZX$ だから $(XZ)^\dagger=-XZ$ である。
$iXZ$ はエルミートになる。
:::

:::example[射影が残す成分]
$P=\ketbra++$ の行列は全成分が $1/2$ である。
$P\ket0=(\ket0+\ket1)/2$、そのノルム二乗は $1/2$。
$\ket-$ に作用させれば零となる。
:::

### 随伴と平均値の実数性

有限次元で $A=A^\dagger$ なら $\mel\psi A\psi$ はその共役と等しく実数である。
逆にすべての $\psi$ について実数なら、$D=(A-A^\dagger)/(2i)$ はエルミートで $\mel\psi D\psi=0$。
$\psi=e_j$ で対角成分が零、$\psi=e_j+e_k$ で $\Re D_{jk}=0$、$\psi=e_j+ie_k$ で $\Im D_{jk}=0$ となる。
従って $D=0$ である。
固有値がすべて実数であることだけでは、この結論は出ない。

:::example[実固有値をもつ非エルミート行列]
$A=\begin{pmatrix}1&1\\0&2\end{pmatrix}$ は固有値が $1,2$ だが、$\psi=(1,i)^\top/\sqrt2$ で $\mel\psi A\psi=(3+i)/2$。
平均を実数として解釈するには状態空間の内積と随伴が必要である。
:::

### 平行移動の随伴と生成子

並進 $T(d)$ が局在した波束を右へ $d$ だけ動かすとする。固定した検出点 $x$ の新しい振幅はもとの $x-d$ の振幅なので、$(T(d)\psi)(x)=\psi(x-d)$ である。全直線での積分の変数を $y=x-d$ と置けば $\langle T(d)\phi|T(d)\psi\rangle=\langle\phi|\psi\rangle$。さらに $T(-d)T(d)=I$ だから $T(d)^\dagger=T(-d)$ である。滑らかで急減少する関数について $T(d)=I-d\partial_x+\cdots=I-idP/\hbar+\cdots$ を比較すると $P=-i\hbar\partial_x$ が得られる。一般の $L^2$ 関数には微分可能性を仮定できないが、有限の $d$ の並進そのものは定義できる。

位置を測る装置との順序も検査できる。$(XT(d)\psi)(x)=x\psi(x-d)$ と $(T(d)X\psi)(x)=(x-d)\psi(x-d)$ の差は $d(T(d)\psi)(x)$ である。$d$ で一次まで比較すれば $\comm X P=i\hbar I$。同じ関係を無限次元行列の成分として扱う場合、どの積も作用できる状態を指定する必要がある。


### 並進と運動量の移動を有限量で比べる

実数 $p$ に対して $(B(p)\psi)(x)=e^{ipx/\hbar}\psi(x)$ と置く。この操作は位置の密度を保ち、Fourier 変換では運動量の振幅を $p$ だけ移動する。$T(d)$ を先に作用させると $B(p)T(d)\psi=e^{ipx/\hbar}\psi(x-d)$、逆順なら指数の $x$ が $x-d$ に変わる。従って全ての $L^2$ 関数に対して

$$
B(p)T(d)=e^{ipd/\hbar}T(d)B(p)
$$

となる。これが有限の移動で表した Weyl の交換関係である。両辺を微分して微分演算子の交換関係を取り出す段階で、共通領域が必要になる。逆向きに、稠密な領域での交換子だけから有限の関係が従うとは限らない。

一組の正準変数では、Weyl の関係を満たす強連続な二つの一径数ユニタリ群が共通の非自明な閉不変部分空間を持たなければ、全直線の位置・運動量の表示とユニタリ同値である。これが Stone–von Neumann の一意性定理である。証明には連続なユニタリ群の既約表現の構成を要するため、定理とその条件については :cite[shimizu-qm] を参照されたい。

### 等長写像が逆を持つための条件

有限次元の正方行列では $V^\dagger V=I$ ならば、$Vv=0$ が $\|v\|=0$ を意味するため単射である。同じ有限次元の間の単射は全射でもあるので $V^{-1}=V^\dagger$、従って $VV^\dagger=I$。この最後の段階に有限次元性が入っている。

半無限の基底 $\ket0,\ket1,\ldots$ に対して $V\ket n=\ket{n+1}$ と置く。どの状態のノルムも保たれるが、$\ket0$ を出力する状態はない。随伴は $V^\dagger\ket0=0$、$V^\dagger\ket{n+1}=\ket n$ と求まる。よって $V^\dagger V=I$、$VV^\dagger=I-\ketbra00$ である。こうした操作を時間発展として採用すると、現在の全状態に対して一意な過去を与えられない。可逆な閉じた系の発展には全射性も要求する理由である。

### 波束を動かす操作と座標を取り換える操作

座標の原点を右へ $a$ 移すと、新しい座標は $x'=x-a$ である。同じ検出点の波動関数は $\psi'(x')=\psi(x'+a)$。これに対し、実験台に対して波束を右へ $a$ 移したときは $\psi_{\rm moved}(x)=\psi(x-a)$ であった。両者の符号が違うのは、同じ $x=0$ で何を比べているかが異なるからである。密度の重心を計算して確かめられる。変数 $y=x-a$ により、実際に動かした後の重心は $\int (y+a)|\psi(y)|^2\dd y=x_0+a$。原点を変えただけの場合の新座標での重心は $x_0-a$ になる。


### 位相板が表すユニタリ演算子を長さで校正する

二つの光路の一方に平行な透明板を挿入する。反射と吸収を無視し、屈折率が波長帯内で一定であるとする。真空中で同じ距離を進む経路との差を位相に直す。



| 設定する量 | 記号 | 計算に用いる値 |
|---|---|---|
| 真空波長 | $\lambda$ | $6.328\times10^{-7}\,\mathrm{m}$ |
| 屈折率 | $n$ | $1.5$ |
| 挿入する板の厚さ | $d$ | $3.164\times10^{-7}\,\mathrm{m}$ |

**光路差。** 光路長の差は板の内部だけに生じる。


$$
\begin{aligned}\delta&=d \left(n - 1\right),\\ \frac{\delta}{\mathrm{m}}&\simeq (3.164\times10^{-7}) \left(1.5 - 1\right)\simeq 1.582\times10^{-7}.\end{aligned}
$$

**相対位相。** 真空中で一波長に対応する位相が二πなので、余分な光路を波長で割る。


$$
\begin{aligned}\phi&=\frac{2 \pi \delta}{\lambda},\\ \frac{\phi}{\mathrm{rad}}&\simeq \frac{2 \pi (1.582\times10^{-7})}{(6.328\times10^{-7})}\simeq 1.5708.\end{aligned}
$$

**対称モードへの射影確率。** 等振幅を準備すると出口の対称成分は二経路の振幅の和の半分である。


$$
\begin{aligned}p_+&=\cos^{2}{\left(\frac{\phi}{2} \right)},\\ p_+&\simeq \cos^{2}{\left(\frac{1.5707963}{2} \right)}\simeq 0.5.\end{aligned}
$$

**反対称モードへの射影確率。** 反対称な検出器では位相因子と一の差が振幅に入る。両確率の和が一になる。


$$
\begin{aligned}p_-&=\sin^{2}{\left(\frac{\phi}{2} \right)},\\ p_-&\simeq \sin^{2}{\left(\frac{1.5707963}{2} \right)}\simeq 0.5.\end{aligned}
$$

**位相の遅れに相当する時間。** 波長帯内で屈折率が一定なら同じ余分な光路から時間差も決まる。


$$
\begin{aligned}\tau&=\frac{\delta}{c},\\ \frac{\tau}{\mathrm{s}}&\simeq \frac{(1.582\times10^{-7})}{(2.9979246\times10^{8})}\simeq 5.277\times10^{-16}.\end{aligned}
$$

厚さを変え、どの出力が明るくなるかを比べる。


| $d\;[\mathrm{m}]$ | $\phi\;[\mathrm{rad}]$ | $p_+$ | $p_-$ |
|---|---|---|---|
| $1.582\times10^{-7}$ | $0.7854$ | $0.85355$ | $0.14645$ |
| $3.164\times10^{-7}$ | $1.5708$ | $0.5$ | $0.5$ |
| $6.328\times10^{-7}$ | $3.1416$ | $0$ | $1$ |

演算子の行列は二経路の基底で $\operatorname{diag}(1,e^{i\phi})$ である。絶対値一の対角成分ゆえに随伴と逆が一致する。吸収があれば失われる光の状態を除いた二経路だけではこの性質は成り立たない。


### 光の二経路で強度差と経路間の位相を読む

二経路の強度差と、両経路を干渉させて読んだ出力の差を比較する。固定した一個の光子について、どちらの読みにも出力ポートに応じて正負のパルスを割り当てる。

一個の検出イベントを符号つき電荷一個として換算する。これは光子が電荷をもつという仮定ではなく、検出器の出力の校正である。

各物理量を $Q_z=a_0 A$、$Q_x=b_0 B$ と分け、$a_0=1.6022\times10^{-19}\,\mathrm C$、$b_0=1.6022\times10^{-19}\,\mathrm C$ とする。計算に必要な無次元の列と行列は

$$
\ket\psi=\begin{pmatrix}\frac{2 \sqrt{5}}{5}\\\frac{\sqrt{5} i}{5}\end{pmatrix},\qquad A=\begin{pmatrix}1 & 0\\0 & -1\end{pmatrix},\qquad B=\begin{pmatrix}0 & 1\\1 & 0\end{pmatrix}.
$$

である。まず作用後の列を作り、元の状態に沿う成分を測る。

$$
\begin{aligned}A\ket\psi&=\begin{pmatrix}\frac{2 \sqrt{5}}{5}\\- \frac{\sqrt{5} i}{5}\end{pmatrix},&\quad \mel\psi A\psi&=\frac{3}{5},\\ B\ket\psi&=\begin{pmatrix}\frac{\sqrt{5} i}{5}\\\frac{2 \sqrt{5}}{5}\end{pmatrix},& \mel\psi B\psi&=0.\end{aligned}
$$

平均をスカラー倍した元の列を引けば、どちらも $\psi$ と直交する。これらの残りを $u,w$ とし、対角に分散、非対角に複素相関を含む小さなグラム行列を作る。

$$
u=\begin{pmatrix}\frac{4 \sqrt{5}}{25}\\- \frac{8 \sqrt{5} i}{25}\end{pmatrix},\quad w=\begin{pmatrix}\frac{\sqrt{5} i}{5}\\\frac{2 \sqrt{5}}{5}\end{pmatrix},\quad \begin{pmatrix}\braket u u&\braket u w\\\braket w u&\braket w w\end{pmatrix}=\begin{pmatrix}\frac{16}{25} & \frac{4 i}{5}\\- \frac{4 i}{5} & 1\end{pmatrix}.
$$

実部は対称化した共分散、虚部は交換子の平均に対応する。積の順序を入れ替えれば相関の虚部だけが反転するため、次の二つを個別に照合できる。

$$
\frac{\langle AB+BA\rangle}{2}-\langle A\rangle\langle B\rangle=0,\qquad \frac{\langle AB-BA\rangle}{2i}=\frac{4}{5}.
$$



| 読み | 単位目盛による表示 | SI での値 |
|---|---|---|
| 第一の標準偏差 | $\frac{4}{5}\,a_0$ | $1.2817\times10^{-19}\,\mathrm C$ |
| 第二の標準偏差 | $1\,b_0$ | $1.6022\times10^{-19}\,\mathrm C$ |

二つの幅の積が取り得る下限を比較するには、無次元の行列式を計算すればよい。

$$
\|u\|^2\|w\|^2-(\Re\braket u w)^2-(\Im\braket u w)^2=\frac{16}{25}-(0)^2-(\frac{4}{5})^2=0.
$$

これは $\|\alpha u+\beta w\|^2\ge0$ に対応する条件である。この例では二成分の空間で両方の残差が同じ一次元の直交補空間に入る。従って行列式が零になることは、実験の二通りの読みがどちらも確定することを意味しない。


## エルミート演算子の固有値問題とスペクトル分解

この節では有限次元に限り、:ref[def-qm-adjoint] のエルミート演算子を扱う。
観測値ごとに状態を分けたいとき、ある一組の固有ベクトルを選ぶよりも、同じ値をもつベクトル全体への射影を求める方が便利である。

:::theorem[有限次元のスペクトル分解]{#thm-qm-spectral}
有限次元複素内積空間で $A=A^\dagger$ なら固有値は実数であり、正規直交固有基底が存在する。
異なる固有値 $a$ ごとの:term[直交射影]{en="orthogonal projection" yomi="ちょっこうしゃえい"}を $P_a$ とすると $A=\sum_a aP_a$、$P_aP_b=\delta_{ab}P_a$、$\sum_a P_a=I$ である。
:::

:::proof
次元が正なら特性多項式には複素根があり、対応する単位固有ベクトル $u$ がある。
$a=\braket u{Au}=\braket{Au}u=a^*$ より固有値は実数である。
$w\perp u$ なら $\braket u{Aw}=\braket{Au}w=0$ だから $u^\perp$ は不変である。
この空間への制限もエルミートであり、次元が一つ減る。
次元零を出発点とする帰納法で、その正規直交固有基底に $u$ を加えればよい。

各固有値に属する基底ベクトルの外積を足したものが $P_a$ である。
各基底ベクトルに両辺を作用させれば上記三つの等式が得られる。
固有空間とその直交補空間への作用が決まるので、射影はその中の基底の選び方によらない。
:::

:::theorem[同時対角化]{#thm-qm-common-basis}
有限次元で二つのエルミート演算子 $A,B$ が共通の正規直交固有基底をもつ必要十分条件は $AB=BA$ である。
:::

:::proof
共通基底があれば両者は対角行列なので可換である。
逆に可換なら $Av=av$ より $A(Bv)=aBv$ となる。
従って各 $A$ 固有空間を $B$ が保つ。
その中で $B$ の制限に :ref[thm-qm-spectral] を適用して得た基底をすべて合わせればよい。
:::

:::example[射影による関数計算]
$X^2=I$ の例では射影は $P_\pm=(I\pm X)/2$ である。
$f(X)=f(1)P_++f(-1)P_-$ と置くと、例えば $e^{itX}=I\cos t+iX\sin t$ となる。
これは偶数次と奇数次に分けたべき級数とも一致する。
:::

:::example[一つの縮退部分を分ける]
$A=\diag(2,2,3)$ に対し、$B$ が最初の二成分を交換し第三成分を固定すると $AB=BA$ である。
第一・第二基底ベクトルを $e_1,e_2$ とし、代わりに $(e_1\pm e_2)/\sqrt2$ をとると、$A$ の固有値を保ったまま $B$ を対角化できる。
:::

### 縮退した固有空間の関数計算

射影の分解から、固有値の集合上の複素数値関数 $f$ に対し $f(A)=\sum_a f(a)P_a$ と定められる。随伴は各係数を共役にするので、$f$ が固有値上で実数値なら $f(A)$ はエルミートである。
二つの関数では $f(A)g(A)=\sum_a f(a)g(a)P_a$。
特に $A\ge0$ なら $B=\sum_a\sqrt a\,P_a$ が正の平方根である。
一意性を見るには、正の $C$ が $C^2=A$ を満たせば $C$ が $A$ と可換であることを使う。
同時固有基底で $C$ の各固有値は非負で二乗が $a$、従ってすべて $\sqrt a$ である。
この平方根を使えば混合状態の不確定性関係を純粋状態と同じ形で証明できる。
数学科目での一般の正規行列に対する分解は :ref[thm-linalg-spectral] である。

:::example[二つの粗い測定を組み合わせる]
$A=\diag(0,0,1,1)$ と $B=\diag(0,1,0,1)$ はそれぞれの固有値が二重に縮退している。
各射影の積をとると四つの一次元射影に分かれる。
同じ $A,B$ の組を指定してなお多次元の空間が残るなら、この二量だけでは状態を一意に指定できない。
:::

### 単位的な軸の交換と反ユニタリ写像

基底を固定して成分をすべて複素共役にする写像を $K$ と書く。$\ket\psi=\sum_jc_j\ket j$ に対し $K\ket\psi=\sum_j c_j^*\ket j$ であり、$K(i\ket\psi)=-iK\ket\psi$ だから線形ではない。一方 $\braket{K\phi}{K\psi}=\braket{\phi}{\psi}^*$ であり、遷移確率は保存される。この性質をもつ共役線形な全単射を :term[反ユニタリ演算子]{en="antiunitary operator" yomi="はんゆにたりえんざんし"} という。たとえば位置表示で $K(-i\hbar\partial_x)K^{-1}=+i\hbar\partial_x$、$KXK^{-1}=X$ である。スピンのない粒子で時間を逆向きにたどる操作に必要な運動量の符号反転が、この二式に含まれる。

:::example[固有値の組と複素共役]
$\psi=(1,i)^\mathsf T/\sqrt2$、$Y=\left(\begin{smallmatrix}0&-i\\i&0\end{smallmatrix}\right)$ とすれば $Y\psi=\psi$。$K\psi=(1,-i)^\mathsf T/\sqrt2$ は $Y$ の固有値 $-1$ に属する。これは $K Y K^{-1}=-Y$ となるためである。行列が実の成分だけをもつ場合の共役の計算と区別する。
:::


### 逆演算子から固有空間を取り出す

固有ベクトルの位相を揃えることなく固有空間を追跡したい場合には、固有値の外側で逆演算子を調べる方法がある。有限次元で $z$ を $A$ のどの固有値とも異なる複素数とし、$R(z)=(zI-A)^{-1}$ を :term[レゾルベント]{en="resolvent" yomi="れぞるべんと"} と呼ぶ。各固有方向では単に数で割るので

$$
R(z)=\sum_a\frac{P_a}{z-a},\qquad \lim_{z\to a}(z-a)R(z)=P_a .
$$

逆行列の差にも正負の検算がある。$R(z)[(wI-A)-(zI-A)]R(w)$ の両端で逆を消去すれば、$R(z)-R(w)=(w-z)R(z)R(w)$。従って $z$ で微分すると $R'(z)=-R(z)^2$ である。固有値で逆行列が壊れることと、各固有空間の情報が極に残ることは両立する。

例えば二つの同等なサイトの間を $t>0$ で結ぶ $A=t\left(\begin{smallmatrix}0&1\\1&0\end{smallmatrix}\right)$ なら、逆行列は $(z^2-t^2)^{-1}\left(\begin{smallmatrix}z&t\\t&z\end{smallmatrix}\right)$。$z=t$ で $z-t$ を掛けて極限を取ると全成分が $1/2$ の対称状態の射影を得る。$z=-t$ では対角が $1/2$、非対角が $-1/2$ の反対称状態の射影を得る。二つを足すと $I$ である。


### 三つの局在モードと縮退した二つの通り道

三つの同等なサイトの間を同じ実振幅で結ぶ一粒子を考える。サイトの順に基底を並べ、オンサイトエネルギーを零とする。どのサイトに置いたかと、どのエネルギーを持つかは別の質問である。

以下の行列は指定した順の正規直交基底で書く。測定量を $H=bM$ とし、単位を担う係数を $b=1.6022\times10^{-22}\,\mathrm J$ とする。サイト間の結合は $1.00\,\mathrm{meV}$ である。


$$
M=\begin{pmatrix}0 & 1 & 1\\1 & 0 & 1\\1 & 1 & 0\end{pmatrix},\qquad \ket{\psi}=\begin{pmatrix}1\\0\\0\end{pmatrix}.
$$

**読み取り値を決める。** $u$ を無次元の候補値とする。$M-uI$ が逆を持つと固有方程式の解は零だけになる。そこで特性多項式を計算して因数を分ける。

$$
\det(uI-M)=\begin{vmatrix}u & -1 & -1\\-1 & u & -1\\-1 & -1 & u\end{vmatrix}=\left(u - 2\right) \left(u + 1\right)^{2}.
$$

各根に対応する空間を取り出すには、$u=\lambda$ で一、ほかの根で零になる多項式を使う。異なる根だけを積に入れて $P_\lambda=\prod_{\mu\ne\lambda}(M-\mu I)/(\lambda-\mu)$ とすると、どの固有ベクトルにもその補間値が掛かる。したがって $P_\lambda$ は行列全体の次数にかかわらず、固有空間への直交射影になる。

$$
P_{-1}=\begin{pmatrix}\frac{2}{3} & - \frac{1}{3} & - \frac{1}{3}\\- \frac{1}{3} & \frac{2}{3} & - \frac{1}{3}\\- \frac{1}{3} & - \frac{1}{3} & \frac{2}{3}\end{pmatrix},\qquad P_{-1}\ket{\psi}=\begin{pmatrix}\frac{2}{3}\\- \frac{1}{3}\\- \frac{1}{3}\end{pmatrix}.
$$



$$
P_{2}=\begin{pmatrix}\frac{1}{3} & \frac{1}{3} & \frac{1}{3}\\\frac{1}{3} & \frac{1}{3} & \frac{1}{3}\\\frac{1}{3} & \frac{1}{3} & \frac{1}{3}\end{pmatrix},\qquad P_{2}\ket{\psi}=\begin{pmatrix}\frac{1}{3}\\\frac{1}{3}\\\frac{1}{3}\end{pmatrix}.
$$

**測定の分布を求める。** 射影した列を先に長さ一に直すと、各結果の出やすさという情報が消えてしまう。ここではそのまま列の長さの二乗を計算する。

| 根 $\lambda$ | 結果 $b\lambda$ | 部分空間の次元 | 準備した状態での確率 |
|---|---|---|---|
| $-1$ | $-1.6022\times10^{-22}\,\mathrm J$ | 2 | $\frac{2}{3}$ |
| $2$ | $3.2044\times10^{-22}\,\mathrm J$ | 1 | $\frac{1}{3}$ |

確率の総和が一になることに加え、もとの基底で計算したモーメントを照合できる。

$$
\begin{aligned} M\ket{\psi}&=\begin{pmatrix}0\\1\\1\end{pmatrix},\\\frac{\langle H\rangle}{b}&=\braket{\psi}{M\psi}=0,\\ \frac{\langle H^2\rangle}{b^2}&=\braket{M\psi}{M\psi}=2,\qquad \frac{(\Delta H)^2}{b^2}=2.\end{aligned}
$$


**続けて別の分析器を通す。** 最後は二番目のサイトで粒子を検出する。 最後に通す状態を $\ket f=\begin{pmatrix}0\\1\\0\end{pmatrix}$ とする。第一段で $\lambda$ を得て第二段を通る同時確率は $|\mel f{P_\lambda}\psi|^2$ である。第一段が起きたという条件を付けるときはこれを $p_\lambda$ で割る。

| 第一段 | 合成振幅 | 同時確率 | 第一段で選別した後の透過率 |
|---|---|---|---|
| $-1$ | $- \frac{1}{3}$ | $\frac{1}{9}$ | $\frac{1}{6}$ |
| $2$ | $\frac{1}{3}$ | $\frac{1}{9}$ | $\frac{1}{3}$ |

第一段を測って値を読み捨てると、最終確率は同時確率の列の和 $\frac{2}{9}$ である。第一段を取り除いたときには $|\braket f\psi|^2=0$。分岐を記録する物理操作の有無を指定すれば、二つの計算を区別できる。

縮退した空間の二本の軸を勝手に測定したことにすると、その内部で保存されていた振幅どうしの干渉も消える。どの装置が読んだかを射影の階数とともに指定する理由がここにある。


## 連続スペクトルと位置・運動量表示

### Fourier 変換としての基底変換

直線上で位置 $X$ は $(X\psi)(x)=x\psi(x)$ と作用する。
$X\psi=x_0\psi$ を満たす通常の関数は一点以外で零になり、$L^2$ では零ベクトルである。
そこで $\braket x\psi=\psi(x)$ は滑らかな試験関数への評価として読み、位置ケットの重なりを :term[デルタ規格化]{en="delta normalization" yomi="でるたきかくか"} で表す。

$$
\begin{aligned}
\braket{x}{x'}&=\delta(x-x'), \\
\int\ketbra x x\,\dd x&=I, \\
\braket x p&=\frac{e^{ipx/\hbar}}{\sqrt{2\pi\hbar}}.
\end{aligned}
$$

最後の核は $-i\hbar\partial_x\braket x p=p\braket x p$ を解き、$\int e^{i(p'-p)x/\hbar}\dd x=2\pi\hbar\delta(p'-p)$ で係数を定めたものである。
その結果

$$
\begin{aligned}
\tilde\psi(p)&=\frac1{\sqrt{2\pi\hbar}}\int e^{-ipx/\hbar}\psi(x)\dd x, \\
\psi(x)&=\frac1{\sqrt{2\pi\hbar}}\int e^{ipx/\hbar}\tilde\psi(p)\dd p.
\end{aligned}
$$

順変換と逆変換を代入するとデルタが元の関数を返す。
位置密度と運動量密度の単位はそれぞれ $\mathrm{m^{-1}}$ と $(\mathrm{kg\,m\,s^{-1}})^{-1}$ である。

:::example[箱での規格化から極限へ]
周期箱では $p_n=2\pi\hbar n/L$、波動関数は $e^{ip_nx/\hbar}/\sqrt L$ である。
$\Delta p=2\pi\hbar/L$ としてこれを $\sqrt{\Delta p}$ で割ると連続核を得る。
確率の和が積分に移るので、固定の一点の確率と密度は区別する。
:::

### 定義域も演算子の一部である

非有界演算子では作用式に加え、作用後も空間内にあるような稠密な定義域 $D(A)$ を指定する。
随伴の定義域は $\psi\mapsto\braket\phi{A\psi}$ が状態のノルムで有界になる $\phi$ の全体である。
ここで $A\subset A^\dagger$ は $D(A)\subset D(A^\dagger)$ かつ $D(A)$ 上で両者の作用が一致するという意味である。
この条件を対称性と呼び、定義域も含め $A=A^\dagger$ が成り立つとき自己共役という。
稠密に定義された自己共役演算子 $A$ には、実数直線の Borel 集合を直交射影に写す一意な射影値測度 $P$ が対応する。$P(\R)=I$、$P(S)P(T)=P(S\cap T)$ であり、互いに素な集合の可算和は、各ベクトルへの作用がノルム収束する意味で射影の和に写る。各ベクトル $\psi$ に対して $\mu_\psi(S)=\|P(S)\psi\|^2$ と置くと、$D(A)=\{\psi:\int a^2\,\dd\mu_\psi(a)<\infty\}$ である。この定義域上で $A=\int a\,\dd P(a)$ と書く。より一般に、複素 Borel 関数 $f$ に対し $\int |f(a)|^2\,\dd\mu_\psi(a)<\infty$ の範囲で $f(A)=\int f(a)\,\dd P(a)$ が定まる。
一般の場合の証明には作用素と測度の収束を扱う関数解析が必要になるため、証明は :cite[reed-simon-1] を参照されたい。以下では位置と運動量の射影を具体的に構成する。
連続基底の扱いは :cite[dirac-principles] を参照されたい。

:::example[区間の運動量]
$[0,L]$ で絶対連続かつ導関数が二乗可積分な関数のうち、両端で零となるものを $P=-i\hbar\partial_x$ の定義域とする。随伴の定義を、まず区間内部に台をもつ滑らかな $\psi$ に適用すると、随伴の定義域に属する $\phi$ は二乗可積分な弱微分をもち、絶対連続な代表を選べる。そこで部分積分すると境界項は $-i\hbar[\phi^*\psi]_0^L$ である。
$\psi(0)=\psi(L)=0$ なら項は消えるが、随伴側には端点条件が要らず、対称でも自己共役ではない。
絶対連続で導関数が二乗可積分な関数に $\psi(L)=e^{i\theta}\psi(0)$ を課せば、境界項が任意の $\psi$ で消える条件は $\phi$ にも同じ端点条件を課す。
この場合は自己共役である。
:::

### スペクトル射影を通常の関数で扱う

位置をちょうど一点に固定する代わりに幅をもつ区間 $S$ を測れば、$(P(S)\psi)(x)=\chi_S(x)\psi(x)$ はノルムが高々一の有界作用素である。
有限測定ではこの射影が直接使える。
たとえば $S=[0,L/3]$、$[0,L]$ 上で一様な波束なら $\norm{P(S)\psi}^2=1/3$。
区間を狭めていくと確率は零へ近づくので、連続変数の一点と離散変数の一点の確率の違いが見える。

全直線での標準的な自己共役な位置演算子の定義域は $D(X)=\{\psi\in L^2:x\psi\in L^2\}$、運動量の定義域は $D(P)=\{\psi\in L^2:\psi'\in L^2\}$ で、微分は弱微分である。
位置が自己共役であることは随伴の定義から確認できる。$\phi\in D(X^\dagger)$ なら $\braket{\phi}{X\psi}=\braket{\eta}{\psi}$ を満たす $\eta\in L^2$ がある。任意の有界区間に台をもつ $\psi$ を入れると $\eta=x\phi$ がその区間上で従うので、全直線で $x\phi\in L^2$ である。逆に $x\phi\in L^2$ なら $\eta=x\phi$ を選べる。Fourier 変換 $F$ は Plancherel の定理により $L^2$ 上のユニタリ写像である。部分積分で $F(-i\hbar\psi')=pF\psi$ が得られ、弱微分にも同じ等式が成り立つ。従って $P$ は $p$ を掛ける自己共役演算子と Fourier 変換で結ばれ、運動量の射影は $F^{-1}\chi_S F$ で与えられる。

位置演算子の $X-\lambda I$ は、$\lambda$ が非実数なら $(x-\lambda)^{-1}$ を掛ける有界な逆をもつ。実数 $\lambda$ の場合は、$\lambda$ の幅 $2\epsilon$ の近傍に台をもつ単位ベクトル $\psi_\epsilon$ に対して $\|(X-\lambda I)\psi_\epsilon\|\le\epsilon$ である。したがって逆は有界になり得ない。一方、$X-\lambda I$ の核は零であり、値域の直交補空間も随伴の核に等しいため零である。値域は稠密なので、実数 $\lambda$ はすべて連続スペクトルに属する。運動量にも Fourier 変換を介して同じ結論が従う。
例えば矩形波束は前者に属するが、端に跳びがあるため後者に属さない。
正規化できることだけから全物理量の有限な二乗平均を仮定してはいけない。

:::example[波束から平面波への極限]
$p_0$ を中心とする Gaussian の運動量密度を幅 $\epsilon$ とする。
任意の有界連続な検出感度 $g$ について $\int g(p)|\tilde\psi_\epsilon(p)|^2\dd p\to g(p_0)$ である。
これは確率分布の収束であり、ノルム一の状態が Hilbert 空間内の平面波へ収束するという主張ではない。
:::

### 連続スペクトルの変数を運動量からエネルギーへ変える

自由粒子の運動量とエネルギーには $E=p^2/(2m)$ の関係がある。$E>0$ では $p=\pm\sqrt{2mE}$ の二つの根がある。根の番号を $\sigma=\pm$ と書いて $\ket{E,\sigma}=\sqrt{m/|p|}\ket{p=\sigma\sqrt{2mE}}$ と置く。$\dd p/\dd E=m/p$ の絶対値の平方根を掛けるのはケットの**振幅**を変換しているためである。実際、二つの枝を別々に積分すると

$$
\sum_{\sigma=\pm}\int_0^\infty\ketbra{E,\sigma}{E,\sigma}\,\dd E
=\int_{-\infty}^\infty\ketbra p p\,\dd p=I .
$$

従って $E$ だけを記録する密度は $w(E)=(m/|p|)\{|\tilde\psi(p)|^2+|\tilde\psi(-p)|^2\}$。二項に干渉項がないのは、それぞれが直交する枝の確率だからである。$p=0$ でヤコビアンが発散しても、$\tilde\psi$ が原点で連続で有界なら $E^{-1/2}$ の特異性は積分可能である。エネルギー零の一点に有限な確率が生じるという解釈は誤りである。


### Gaussian の Fourier 変換を微分方程式から求める

$\psi(x)=C\exp[-(x-d)^2/(4s^2)+ip_0 x/\hbar]$ を全直線で規格化する。絶対値二乗の積分から $C=(2\pi s^2)^{-1/4}$ と選べる。$y=x-d$、$k=(p-p_0)/\hbar$ とおき、$F(k)=\int e^{-y^2/(4s^2)}e^{-iky}\dd y$ とする。無限遠で境界項が消えるので
$$
F'(k)=2is^2\int \partial_y(e^{-y^2/(4s^2)}) e^{-iky}\dd y
=-2s^2 k F(k).
$$
この一次方程式は $F(k)=F(0)e^{-s^2k^2}$ を与える。$F(0)=2s\sqrt\pi$ を入れ、変数の平行移動で外に出した位相を戻すと
$$
\tilde\psi(p)=\left(\frac{2s^2}{\pi\hbar^2}\right)^{1/4}
e^{-i(p-p_0)d/\hbar}e^{-s^2(p-p_0)^2/\hbar^2}.
$$
従って位置を $d$ だけ変えても運動量密度は変わらず、位置表示で傾き $p_0/\hbar$ の位相を付けると運動量分布が $p_0$ だけ動く。運動量の分散は $\hbar^2/(4s^2)$。位置の幅 $s=0.50\,\mathrm{nm}$ の波束なら $\Delta p=1.05\times10^{-25}\,\mathrm{kg\,m\,s^{-1}}$ となる。両分布の中心も幅も、状態の同じ規格化から得られたものである。

:::remark[デルタ関数の次元]
正規化条件 $\int\delta(p-p')\dd p=1$ より、このデルタ関数は運動量の逆数の次元をもつ。従って $\delta(p-p')$ と $\delta(k-k')$ は同じ数値にならない。$p=\hbar k$ なら $\delta(p-p')=\hbar^{-1}\delta(k-k')$ である。箱規格化から連続規格化に移る際、一つのモードの波動関数を $\sqrt{\Delta p}$ で割ったのは、この密度の変換に振幅の二乗が対応するからである。
:::


### 波長を絞ったビームの運動量表示

速度選別器を通った右向きの粒子から、区間 $0<p<b$ の運動量だけを取り出す。運動量表示の振幅を窓内で $p(b-p)$ に比例させ、両端で零、外側で零になるように準備する。

無次元化は $x=p/b$、$b=\hbar/(1.00\,\mathrm{nm})$ である。運動量を長さで割って無次元化することはできない。

以下では無次元の読みを $x$、その一単位に対応する実際の目盛を $b=1.0546\times10^{-25}\,\mathrm{kg\,m\,s^{-1}}$ とする。まず全領域を積分したときに確率が一になるよう係数を決める。形だけを表す関数を $w$ とすると

$$
\begin{aligned}w(x)&=x^{2} \left(x - 1\right)^{2},\qquad 0<x<1,\\ Z&=\int_{0}^{1}w(x)\,\dd x=\frac{1}{30},\\ f(x)&=\frac{w(x)}{Z}=30 x^{2} \left(x - 1\right)^{2}.\end{aligned}
$$

目盛を有限幅に区切ると、区間に入る確率は曲線下の面積で決まる。実際の幅に対する密度は $f(x)/b$ である。したがって幅を変換してから積分しても、無次元の区間の両端で原始関数を引いても同じ値を得る。ここで使う原始関数は

$$
F(x)=x^{3} \left(6 x^{2} - 15 x + 10\right),\qquad F^\prime(x)=f(x).
$$

各区間に入る割合を求める。境界ちょうどにある点は連続分布の確率を持たないため、隣り合う区間に二重に数える心配はない。

| 無次元の区間 | $F(右端)-F(左端)$ | 確率の概数 |
|---|---|---|
| $0\ \text{–}\ \frac{1}{4}$ | $\frac{53}{512}$ | 0.1035 |
| $\frac{1}{4}\ \text{–}\ \frac{3}{4}$ | $\frac{203}{256}$ | 0.793 |
| $\frac{3}{4}\ \text{–}\ 1$ | $\frac{53}{512}$ | 0.1035 |

次に、このうち $\frac{1}{4}<x<\frac{3}{4}$ だけを受ける検出器を置く。ゲートを通過した事象を分母にして平均すると

$$
\begin{aligned}p_{\rm gate}&=\frac{203}{256},\\\int_{\frac{1}{4}}^{\frac{3}{4}}x f(x)\,\dd x&=\frac{203}{512},\\ \overline{x}_{\rm gate}&=\frac{1}{2}\simeq 0.5.\end{aligned}
$$

となる。これを元の単位に戻した条件付き平均は $5.2729\times10^{-26}\,\mathrm{kg\,m\,s^{-1}}$ である。一方、入射した集団全体では二つの積分

$$
\begin{aligned}\langle x\rangle&=\int_{0}^{1}x f(x)\,\dd x=\frac{1}{2},\\\langle x^2\rangle&=\int_{0}^{1}x^2 f(x)\,\dd x=\frac{2}{7}\end{aligned}
$$

を使える。差 $\langle x^2\rangle-\langle x\rangle^2=\frac{1}{28}$ の平方根に $b$ を掛ければ、運動量の標準偏差は $1.993\times10^{-26}\,\mathrm{kg\,m\,s^{-1}}$ となる。

ゲート前後の分布はどちらも正規化できる。元ビームのうち何割を捨てたかを知りたい実験では、条件付き密度とは別にゲートの確率を保存しておく。


### 周期箱に閉じ込めた有限幅の波束

長さ $L$ の環を $-1/2<x<1/2$ に写す。全周の半分だけ一様な状態から、読める運動量に制限がある装置でどこまで波形を復元できるかを調べる。

周期一の無次元座標 $x$ を用いる。実際の一周期は $L=4\times10^{-9}\,\mathrm m$、正規直交する基底は $u_n=e^{2\pi i n x}$ である。隣接する整数の差が零以外なら、積分に含まれる指数関数がちょうど整数周期を回り、内積は零になる。この例では $- \frac{1}{4}<x<\frac{1}{4}$ の波形を $v(x)=1$ とし、一周期の残る部分は零とする。

$$
Z=\int_{- \frac{1}{4}}^{\frac{1}{4}}|v(x)|^2\,\dd x=\frac{1}{2},\qquad c_n=\frac{1}{\sqrt Z}\int_{- \frac{1}{4}}^{\frac{1}{4}}e^{-2\pi i nx}v(x)\,\dd x.
$$

零モードの係数は $c_0=\frac{\sqrt{2}}{2}$ である。他のモードは相対位相が一周する回数によって相殺の程度が違う。この環の一目盛は $2\pi\hbar/L$ である。$n\ne0$ の積分は両端の指数の差を $-2\pi i n$ で割ることで得られ、$\sin(n\pi/2)/n$ に比例する。



| モード $n$ | 複素振幅 $c_n$ | $\vert c_n\vert ^2$ | 重みの概数 |
|---|---|---|---|
| -4 | $0$ | $0$ | $0$ |
| -3 | $- \frac{\sqrt{2}}{3 \pi}$ | $\frac{2}{9 \pi^{2}}$ | $0.022516$ |
| -2 | $0$ | $0$ | $0$ |
| -1 | $\frac{\sqrt{2}}{\pi}$ | $\frac{2}{\pi^{2}}$ | $0.20264$ |
| 0 | $\frac{\sqrt{2}}{2}$ | $\frac{1}{2}$ | $0.5$ |
| 1 | $\frac{\sqrt{2}}{\pi}$ | $\frac{2}{\pi^{2}}$ | $0.20264$ |
| 2 | $0$ | $0$ | $0$ |
| 3 | $- \frac{\sqrt{2}}{3 \pi}$ | $\frac{2}{9 \pi^{2}}$ | $0.022516$ |
| 4 | $0$ | $0$ | $0$ |

ここに挙げた有限個のモードだけを記録する。射影後のノルム二乗は、列から選んだ確率を互いの位相と無関係に足すことで得られる。一方、同じモードから位置波形を再構成する場合には複素係数を保持して足さねばならない。

$$
\begin{aligned}\psi_{\rm kept}(x)&=\sum_{n\ {\rm kept}}c_n e^{2\pi i nx},\\ \|\psi_{\rm kept}\|^2&=\frac{80 + 9 \pi^{2}}{18 \pi^{2}}\simeq 0.95032,\\ \|\psi-\psi_{\rm kept}\|^2&=1-\sum_{n\ {\rm kept}}|c_n|^2\simeq 0.049684.\end{aligned}
$$

奇数モードの振幅は $1/n$ で減る。従って確率の和は収束するが、運動量の二乗を掛けた和の各奇数項は零に近づかず、二乗平均は発散する。

境界で波形が跳ぶ状態は $L^2$ に属しても運動量の定義域にはない、という区別を実際の係数列から読める。


### 同じ周期箱で波数・運動量・エネルギーの目盛を作る

細い環に閉じ込めた電子を考え、横方向の励起を凍結する。磁束のない環を一周した波動関数には位相のずれがない。電子のある一方向のモードに注目する。



| 設定する量 | 記号 | 計算に用いる値 |
|---|---|---|
| 環の周長 | $L$ | $1\times10^{-7}\,\mathrm{m}$ |
| 巻き数 | $n$ | $3$ |

**波数。** 一周で位相が整数倍の二πだけ進む条件を周長で割る。


$$
\begin{aligned}k_n&=\frac{2 \pi n}{L},\\ \frac{k_n}{\mathrm{m^{-1}}}&\simeq \frac{2 \pi 3}{(1\times10^{-7})}\simeq 1.885\times10^{8}.\end{aligned}
$$

**運動量固有値。** 波数演算子の固有値にディラック定数を掛けて機械的運動量の単位に直す。


$$
\begin{aligned}p_n&=\hbar k_n,\\ \frac{p_n}{\mathrm{kg\,m\,s^{-1}}}&\simeq (1.0545718\times10^{-34}) (1.8849556\times10^{8})\simeq 1.9878\times10^{-26}.\end{aligned}
$$

**運動エネルギー。** 環の接線方向にだけ自由に動く電子の分散関係を使う。


$$
\begin{aligned}E_n&=\frac{p_n^{2}}{2 m_e},\\ \frac{E_n}{\mathrm{J}}&\simeq \frac{\left((1.987821\times10^{-26})\right)^{2}}{2 (9.1093837\times10^{-31})}\simeq 2.1689\times10^{-22}.\end{aligned}
$$

**運動量の刻み。** 隣り合う巻き数の差が一であることから求める。


$$
\begin{aligned}\Delta p&=\frac{h}{L},\\ \frac{\Delta p}{\mathrm{kg\,m\,s^{-1}}}&\simeq \frac{(6.6260701\times10^{-34})}{(1\times10^{-7})}\simeq 6.6261\times10^{-27}.\end{aligned}
$$

**隣の正向き準位までの差。** 二乗してから差を取る。エネルギーの刻みは運動量の刻みと異なり一定ではない。


$$
\begin{aligned}\delta E&=\frac{- p_n^{2} + \left(\Delta p + p_n\right)^{2}}{2 m_e},\\ \frac{\delta E}{\mathrm{J}}&\simeq \frac{- \left((1.987821\times10^{-26})\right)^{2} + \left((6.6260701\times10^{-27}) + (1.987821\times10^{-26})\right)^{2}}{2 (9.1093837\times10^{-31})}\simeq 1.6869\times10^{-22}.\end{aligned}
$$

環を固定しても巻き数を変えると隣接準位の間隔は動く。


| $n$ | $p_n\;[\mathrm{kg\,m\,s^{-1}}]$ | $E_n\;[\mathrm{J}]$ | $\delta E\;[\mathrm{J}]$ |
|---|---|---|---|
| $1$ | $6.6261\times10^{-27}$ | $2.4099\times10^{-23}$ | $7.2296\times10^{-23}$ |
| $2$ | $1.3252\times10^{-26}$ | $9.6395\times10^{-23}$ | $1.2049\times10^{-22}$ |
| $3$ | $1.9878\times10^{-26}$ | $2.1689\times10^{-22}$ | $1.6869\times10^{-22}$ |
| $4$ | $2.6504\times10^{-26}$ | $3.8558\times10^{-22}$ | $2.1689\times10^{-22}$ |

逆向きの巻き数は同じエネルギーをもつ。連続極限でエネルギーから運動量に戻すときに二つの枝を数えた理由が、有限の箱ですでに見えている。


## 演習問題

::::exercise{level=basic}
$\ket u=(2,i)^\top/\sqrt5$、$\ket v=(i,1)^\top/\sqrt2$ の内積を順序を変えて計算せよ。

:::solution
$\braket u v=(2i-i)/\sqrt{10}=i/\sqrt{10}$。
共役対称性から $\braket v u=-i/\sqrt{10}$ である。
:::
::::

::::exercise{level=basic}
$v=(1,2)^\top$ の $u=(1,1)^\top$ に沿う成分と直交成分を求めよ。

:::solution
$\braket u v/\braket u u=3/2$ なので平行成分は $(3/2,3/2)^\top$。
差は $(-1/2,1/2)^\top$ であり $u$ との内積は零である。
:::
::::

::::exercise{level=basic}
$A=\ketbra0{1}$ の $A^\dagger A$ と $AA^\dagger$ を求めよ。

:::solution
随伴は $\ketbra1{0}$。
中央の内積をとると順に $\ketbra1{1}$ と $\ketbra0{0}$ であり、それぞれ入力側と出力側の部分空間への射影である。
:::
::::

::::exercise{level=basic}
$A=\begin{pmatrix}2&1\\1&2\end{pmatrix}$ のスペクトル分解を書け。

:::solution
$(1,\pm1)^\top/\sqrt2$ に作用すると固有値は順に $3,1$。
従って $A=3\ketbra++ +\ketbra--$ である。
:::
::::

::::exercise{level=advanced}
異なる固有値に属するエルミート演算子の固有ベクトルが直交することを、帰納法を用いず示せ。

:::solution
$Au=au,Av=bv$ なら $b\braket u v=\braket u{Av}=\braket{Au}v=a\braket u v$。
実固有値の差は零でないので $\braket u v=0$ である。
:::
::::

::::exercise{level=advanced}
有限次元で $P^2=P=P^\dagger$ とする。
$\Img P$ と $\Ker P$ が直交補空間であることを示せ。

:::solution
任意の $v$ は $Pv+(I-P)v$ に分かれ、後者は $P(I-P)=0$ から核に属する。
$w\in\Ker P$ なら $\braket{Pu}w=\braket u{Pw}=0$。
したがってこの分解は直交直和である。
:::
::::

::::exercise{level=advanced}
有限次元のエルミート演算子 $A$ について、固有値 $a$ の射影を $A$ の多項式で表せ。

:::solution
$f_a(t)=\prod_{b\ne a}(t-b)/(a-b)$ と置く。
各相異なる固有値で $f_a(b)=\delta_{ab}$ なので固有基底全体で $f_a(A)=P_a$ が成り立つ。
:::
::::

::::exercise{level=advanced}
区間の位相つき端点条件から運動量の許容値を求めよ。

:::solution
固有方程式を解くと $\psi(x)=Ce^{ipx/\hbar}$。
端点から $e^{ipL/\hbar}=e^{i\theta}$、従って $p=\hbar(2\pi n+\theta)/L$、$n\in\Z$ となる。
規格化は $|C|=L^{-1/2}$ である。
:::
::::

::::exercise{level=advanced}
$\psi(x)=Ae^{-\kappa|x|}$（$\kappa>0$）を規格化し、運動量振幅とその次元を求めよ。
正負の運動量を測る確率も求めよ。

:::solution
位相を選び $A>0$ とする。
規格化積分は $2A^2/(2\kappa)=1$ なので $A=\sqrt\kappa$。
偶関数であるため Fourier 積分は余弦変換にできる。

$$
\tilde\psi(p)=\frac{2\sqrt\kappa}{\sqrt{2\pi\hbar}}\Re\int_0^\infty e^{-(\kappa+ip/\hbar)x}\dd x=\frac{2\kappa^{3/2}}{\sqrt{2\pi\hbar}\,[\kappa^2+(p/\hbar)^2]} .
$$

確率 $\int|\tilde\psi(p)|^2\dd p$ が無次元なので振幅の単位は運動量の逆平方根である。
密度は偶関数で一点の確率は零だから正負は各 $1/2$。
さらに $p^{-4}$ の裾により運動量の二乗平均は有限だが四乗平均は発散する。
:::
::::

::::exercise{level=advanced}
$A=\diag(1,1,-1)$ と可換なエルミート行列 $B$ の最も一般的な形を書け。
$B_{12}=i$、対角成分がすべて零の例で共通固有基底を求めよ。

:::solution
$A$ の固有値を $a_i$ とすれば $[A,B]_{ij}=(a_i-a_j)B_{ij}$。
異なる固有値を結ぶ成分が零であれば十分であり、それ以外にはエルミート性だけを課す。
従って

$$
B=\begin{pmatrix}r&z&0\\z^*&s&0\\0&0&t\end{pmatrix},\qquad r,s,t\in\R,\quad z\in\C .
$$

指定された例の二次元ブロックでは $(1,-i)^\top/\sqrt2$ が固有値 $+1$、$(1,i)^\top/\sqrt2$ が $-1$ を与える。
三成分目を零としてこの二本を延長し、$(0,0,1)^\top$ を加える。
最後のベクトルで $A=-1,B=0$、前二者では $A=1$ である。
:::
::::

::::exercise{level=advanced}
非零ケット $u$ への射影を $P$ とし、任意の規格化した $v$ に対し $\|v-e^{i\theta}u/\|u\|\|$ を最小にする位相を求めよ。射影後の確率との関係も導け。

:::solution
まず $w=u/\|u\|$ とする。二乗を展開すると $2-2\Re(e^{-i\theta}\braket w v)$。$\braket w v=re^{i\alpha}\ne0$ なら $\theta=\alpha$ で実部が最大となり、最小距離の二乗は $2(1-r)$ である。直交する場合はどの位相でも同じである。射影の確率は $\|Pv\|^2=r^2$ なので、全体位相を無視した距離はこの確率だけで決まる。
:::
::::


::::exercise{level=advanced}
次の計算には :eqref[eq-qm-dual-basis] を使える。$v_1=(1,0)^\mathsf T$、$v_2=(1,1)^\mathsf T/\sqrt2$ とする。双対基底を求め、$\psi=(0,i)^\mathsf T$ の展開を復元せよ。係数の絶対値二乗を足すとどうなるか。

:::solution
$V$ は上三角で、$V^{-1}=\left(\begin{smallmatrix}1&-1\\0&\sqrt2\end{smallmatrix}\right)$。双対ブラはこの二行であり、双対ケットは $w_1=(1,-1)^\mathsf T,w_2=(0,\sqrt2)^\mathsf T$。ケットについて線形なので $c_1=-i,c_2=i\sqrt2$。戻すと $-iv_1+i\sqrt2v_2=(0,i)^\mathsf T$ である。係数の絶対値二乗の和は三だが、実際のノルム二乗は $c^\dagger Gc=1$。斜交座標の係数を射影確率として足せない例である。
:::
::::


## 確認テスト

::::quiz
$\ket v$ を $i\ket v$ にすると対応するブラに掛かる係数はどれか。

- [ ] $i$
- [x] $-i$
- [ ] $1$

:::explanation
ブラへの対応は共役線形だから $i^*=-i$ を掛ける。
:::
::::

::::quiz{answer="2" tol="1e-6"}
$\diag(2,2,5)$ の固有値 $2$ への射影の階数はいくつか。

:::explanation
最初の二つの座標を保つ射影 $\diag(1,1,0)$ の像は二次元である。
:::
::::

::::quiz{answer="0.64" tol="1e-6"}
$\ket\psi=(3\ket0+4\ket1)/5$ の $\ket1$ 成分の絶対値二乗はいくつか。

:::explanation
直交性から係数は $4/5$ なので絶対値二乗は $16/25=0.64$ である。
:::
::::

::::quiz
全直線の位置ケットをどう解釈するか。

- [ ] ノルムが一の局在関数
- [x] 試験関数の位置での値を取り出す一般化固有状態
- [ ] 通常の有限次元ベクトル

:::explanation
一点だけに台をもつ二乗可積分関数は零元である。
デルタを積分の中で作用させて位置の振幅を得る。
:::
::::

::::quiz{answer="0.5" tol="0.0001"}
偏光の基底候補を水平と円偏光の二つにとる。グラム行列の行列式はいくらか。

:::explanation
二つの列のノルムは一、内積の絶対値は $1/\sqrt2$ であるから、行列式は $1-1/2=1/2$ となる。正でも一とは限らない。
:::
::::

---
title: 量子ビットと量子状態
description: 量子ビットの状態を Bloch ベクトルで表し、複合系の縮約、純粋化、一般化測定を有限次元で導く。
notation:
  - symbol: '\ket{\psi}'
    meaning: 状態を表す規格化ベクトル
  - symbol: 'X,Y,Z'
    meaning: Pauli 行列
  - symbol: '\rho'
    meaning: 密度演算子
  - symbol: '\Tr_B'
    meaning: 系 B の部分トレース
  - symbol: 'E_a'
    meaning: 測定結果 a の効果
---

## この章の目標

- 相対位相と大域位相を区別し、Bloch 球上の点から測定確率を求める。
- テンソル積で複合系の状態を表し、積状態かどうかを判定する。
- 密度行列と縮約状態を計算し、補助系の最小次元を求める。
- 一般化測定の確率と状態更新を区別し、射影測定による実装を構成する。

## 量子ビットと Bloch 球

### 二準位系と位相

古典ビットを読むと、値は $0$ または $1$ である。一方、二準位の量子系では、読み出し装置がどの方向を測るかによって確率が変わる。この差を記述するには、二つの実確率に加えて振幅の相対位相が必要である。

:::definition[量子ビット]{#def-qi-qubit}
:term[量子ビット]{en="qubit" yomi="りょうしびっと"}は、状態空間が $\C^2$ の量子系である。標準正規直交基底 $\ket0=(1,0)^\top$、$\ket1=(0,1)^\top$ を計算基底という。単位ベクトル $\ket\psi=\alpha\ket0+\beta\ket1$ は純粋状態を表し、$\abs{\alpha}^2+\abs{\beta}^2=1$ である。
:::

Dirac 記法の内積は第 2 引数に線形であり、$\bra\psi$ は $\ket\psi$ の共役転置である。[内積空間](/math/linear-algebra/inner-product-spaces/) の内積と照合するときは引数の順序に注意する。:ref[post-qm-state] と :ref[post-qm-born] により、$\ket\psi$ と $e^{i\gamma}\ket\psi$ はすべての射影測定で同じ確率を与える。実際、振幅に掛かった $e^{i\gamma}$ は絶対値の二乗を取ると消える。

逆に、単位ベクトル $\ket\phi$ が $\ket\psi$ とすべての射影測定で同じ確率を与えるなら、$\ketbra\psi\psi$ の結果確率は 1 である。Cauchy–Schwarz の等号条件より $\ket\phi=e^{i\gamma}\ket\psi$ となる。この共通の位相因子を大域位相という。

:::example[相対位相を読む]
$\ket{\pm}=(\ket0\pm\ket1)/\sqrt2$ を計算基底で測ると、どちらの状態でも $0,1$ が等確率で出る。一方、$\braket{+}{-}=0$ であり、基底 $(\ket+,\ket-)$ で測ると、$\ket+$ と $\ket-$ は確率 1 で区別できる。振幅の符号は相対位相なので測定に現れる。
:::

### Pauli 行列と座標

:term[Pauli 行列]{en="Pauli matrices" yomi="ぱうりぎょうれつ"}を次で定める。

$$
X=\begin{pmatrix}0&1\\1&0\end{pmatrix},\qquad
Y=\begin{pmatrix}0&-i\\i&0\end{pmatrix},\qquad
Z=\begin{pmatrix}1&0\\0&-1\end{pmatrix}.

$$

これらは自己共役であり、各自の二乗は $I$ である。行列積を取ると $XY=iZ$、$YZ=iX$、$ZX=iY$ となり、逆順の積は随伴を取ると符号が反転する。$\boldsymbol\sigma=(X,Y,Z)$ とおけば、二つの実ベクトル $\boldsymbol a,\boldsymbol b$ について、対角項と反対称項を分けて次式を得る。

$$

(\boldsymbol a\cdot\boldsymbol\sigma)(\boldsymbol b\cdot\boldsymbol\sigma)
=(\boldsymbol a\cdot\boldsymbol b)I+i(\boldsymbol a\times\boldsymbol b)\cdot\boldsymbol\sigma.
\label{eq-qi-pauli-product}

$$

:::proposition[Bloch 球の座標]{#prop-qi-bloch-pure}
任意の純粋状態は、位相を除いて $\cos(\theta/2)\ket0+e^{i\phi}\sin(\theta/2)\ket1$ と書ける。$0\le\theta\le\pi$、$0\le\phi<2\pi$ とすると、Pauli 行列の期待値は $\boldsymbol r=(\sin\theta\cos\phi,\sin\theta\sin\phi,\cos\theta)$ である。
:::

:::proof
$\alpha\ne0$ の場合は $e^{-i\arg\alpha}$ を全体に掛ける。これにより第 1 成分は非負となり、$\cos(\theta/2)=\abs\alpha$ とおけば、第 2 成分の絶対値は規格化より $\sin(\theta/2)$ となる。$\alpha=0$ の場合は $\theta=\pi$ と取れる。期待値は行列積から順に $2\Re(\bar\alpha\beta)$、$2\Im(\bar\alpha\beta)$、$\abs\alpha^2-\abs\beta^2$ である。半角公式を代入すれば主張を得る。極では $\phi$ を変えても同じ点を表す。
:::

:term[Bloch 球]{en="Bloch sphere" yomi="ぶろっほきゅう"}の単位球面が純粋状態に対応する。単位ベクトル $\boldsymbol n$ に対し、式 :eqref[eq-qi-pauli-product] により、$(I\pm\boldsymbol n\cdot\boldsymbol\sigma)/2$ は直交する射影で、和は $I$ である。この射影測定の確率は $(1\pm\boldsymbol n\cdot\boldsymbol r)/2$ となる。球上の方向が直交する場合は両結果が等確率で出る。状態ベクトルの内積が零になるのは、球上の方向が反対を向く場合である。

:::example[虚数成分を測る]
$(\ket0+i\ket1)/\sqrt2$ について $\bar\alpha\beta=i/2$ であるから $\boldsymbol r=(0,1,0)$ となる。$Y$ の結果 $+1$ は確実に得られる。一方、$X$ の結果 $+1$ と $Z$ の結果 $+1$ はそれぞれ確率 $1/2$ である。
:::

### 測定軸、射影、分散

同じ実験を複数回行うと、平均値だけでなく測定結果の散らばりも見える。自己共役演算子 $L=\sum_j\ell_j P_j$ に対し、平均は $\langle L\rangle_\psi=\sum_j\ell_j\mel{\psi}{P_j}{\psi}=\mel{\psi}{L}{\psi}$、分散は $\operatorname{Var}_\psi L=\mel{\psi}{L^2}{\psi}-\langle L\rangle_\psi^2$ である。ここで $P_j$ は異なる固有値 $\ell_j$ ごとの射影とし、縮退があっても測定ラベルは固有値とする。分散は $\norm{(L-\langle L\rangle_\psi I)\ket\psi}^2$ に等しいので、零になるのは $\ket\psi$ が $L$ の固有ベクトルであるときに限る。

:::proposition[二つの測定の分散]{#prop-qi-pauli-uncertainty}
単位ベクトル $\ket\psi$ と二つの自己共役演算子 $A,B$ について
$$
\operatorname{Var}_\psi A\,\operatorname{Var}_\psi B
\ge
\frac14\abs{\mel{\psi}{[A,B]}{\psi}}^2
+\frac14\abs{\mel{\psi}{\{\Delta A,\Delta B\}}{\psi}}^2
$$
である。ここで $\Delta A=A-\langle A\rangle_\psi I$、$\Delta B=B-\langle B\rangle_\psi I$ であり、$\{D,E\}=DE+ED$ とする。
:::

:::proof
$u=\Delta A\ket\psi$、$v=\Delta B\ket\psi$ とおく。内積の実部は $\mel{\psi}{\{\Delta A,\Delta B\}}{\psi}/2$ であり、虚部は $\mel{\psi}{[A,B]}{\psi}/(2i)$ である。これは $\braket{u}{v}=\mel{\psi}{\Delta A\Delta B}{\psi}$ と、その複素共役が積の逆順に対応することから分かる。Cauchy–Schwarz の $\norm u^2\norm v^2\ge\abs{\braket{u}{v}}^2$ に両成分を代入する。等号は $u,v$ が線形従属である場合に成立する。二次元の純粋状態では $u,v$ はともに $\psi$ と直交するため、常に等号になる。
:::

:::example[斜めの二つの測定軸]
$\ket\psi=(\sqrt3\ket0+e^{i\pi/3}\ket1)/2$ とし、$Z$ と $N=(X+Z)/\sqrt2$ を比較する。Bloch ベクトルは $(\sqrt3/4,3/4,1/2)$ である。したがって $Z$ の $+1$ 確率は $3/4$、$N$ の $+1$ 確率は $\{1+(\sqrt3+2)/(4\sqrt2)\}/2$ である。二つの結果はどちらも $\pm1$ なので、二乗の平均は 1 になる。分散はそれぞれ $3/4$ と $1-(\sqrt3+2)^2/32$ である。

可換子は $[Z,N]=\sqrt2 iY$ なので、可換子から来る下界は $9/32$ である。反可換子の項には $\langle ZN+NZ\rangle-2\langle Z\rangle\langle N\rangle=(6-\sqrt3)/(4\sqrt2)$ を代入する。二乗して 4 で割った値を $9/32$ に足せば、二つの分散の積 $3(25-4\sqrt3)/128$ が得られる。
:::

:::example[逆向きと直交する軸]
$\ket\phi=(3\ket0+4i\ket1)/5$ では $\boldsymbol r=(0,24/25,-7/25)$ である。$\boldsymbol n=(0,7/25,24/25)$ を選ぶと $\boldsymbol n\cdot\boldsymbol r=0$ になる。したがって $N=(7Y+24Z)/25$ のどちらの結果も確率 $1/2$ で出て、分散は 1 である。

一方、軸を $-\boldsymbol r$ に取ると、$+1$ の確率は 0 である。対応する $+1$ の固有ベクトルは $(4\ket0-3i\ket1)/5$ と取れる。二つの ket の内積を直接取ると $12/25+(-4i)(-3i)/25=0$ になる。測定軸の直交と、測定固有ベクトルの直交は、このように異なる操作で調べる。
:::

### Bloch ベクトルで測定基底を再構成する

単位ベクトル $\boldsymbol n=(\sin\theta\cos\phi,\sin\theta\sin\phi,\cos\theta)$ に対する固有基底を
$$
\begin{aligned}
\ket{n_+}&=\cos(\theta/2)\ket0+e^{i\phi}\sin(\theta/2)\ket1,\\
\ket{n_-}&=-e^{-i\phi}\sin(\theta/2)\ket0+\cos(\theta/2)\ket1
\end{aligned}
$$
と取る。この二つのノルムは 1、内積は 0 であるから、第一の射影が $(I+\boldsymbol n\cdot\boldsymbol\sigma)/2$ なら第二の射影はその補射影である。第一の外積の対角成分は $(1\pm\cos\theta)/2$、上の非対角成分は $e^{-i\phi}\sin\theta/2$ である。したがってこの基底で測ることと $\boldsymbol n$ の向きに Pauli 測定を行うことは等価である。

:::example[準備軸と測定軸の間の角度]
$\ket{n_+}$ を準備し、$\boldsymbol m$ の方向を測るとき、二軸の角を $\gamma$ とすれば確率は $(1+\cos\gamma)/2=\cos^2(\gamma/2)$ である。例えば $\boldsymbol n=(3/5,0,4/5)$、$\boldsymbol m=(0,3/5,4/5)$ では内積が $16/25$ なので $+1$ 確率は $41/50$ となる。準備の固有ベクトルと測定の固有ベクトルをそれぞれ $(3,1)^\top/\sqrt{10}$ と $(3,i)^\top/\sqrt{10}$ に取っても、重なりの二乗は $\abs{9-i}^2/100=41/50$ である。半角が現れるのは ket のパラメータを角 $\theta$ でなく $\theta/2$ に取るためである。
:::

### 幾何位相が測定に現れる条件

三つの正規化ベクトル $\ket a,\ket b,\ket c$ について $B=\braket ab\braket bc\braket ca$ を考える。ベクトルごとに位相を変えても、ある ket に掛けた位相因子は同じベクトルが bra として現れるときに共役となるので、積は変わらない。$B=\Tr(P_aP_bP_c)$ であり、Pauli 展開から
$$
B=\frac{1+\boldsymbol a\cdot\boldsymbol b+\boldsymbol b\cdot\boldsymbol c+\boldsymbol c\cdot\boldsymbol a+i(\boldsymbol a\times\boldsymbol b)\cdot\boldsymbol c}{4}
$$
を得る。ただし右辺の太字は各射影の Bloch ベクトルである。単独の Pauli 行列は跡零であり、二つの積は内積の成分、三つの積は式 :eqref[eq-qi-pauli-product] の反対称項が跡を持つためである。従って順序を逆にすると $B$ は複素共役になる。三つの測定を順に行って成功した回数だけを数える場合は、位相を読み出す干渉経路がないので、この虚部は出現しない。

:::example[三つの軸の順序を干渉で比較する]
$\ket a=\ket0,\ket b=\ket+,\ket c=(\ket0+i\ket1)/\sqrt2$ とする。三軸は $\boldsymbol e_z,\boldsymbol e_x,\boldsymbol e_y$ なので $B=(1+i)/4$ である。直接、$\braket ab=1/\sqrt2$、$\braket bc=(1+i)/2$、$\braket ca=1/\sqrt2$ を掛けても同じ値になる。この積の位相は $\pi/4$ である。$P_b$ と $P_c$ の順を交換して $\mel a{P_cP_b}a$ を測ると $B^*=(1-i)/4$ になる。

複素振幅を比較するには、制御ビットを $\ket+$、データを $\ket a$ に準備する。演算子 $K=\ketbra00\otimes I+\ketbra11\otimes P_bP_c$ は、各射影がベクトルの長さを増やさないので $K^\dagger K\le I$ を満たす。従って二結果の測定演算子を $K,\sqrt{I-K^\dagger K}$ に取れる。これは結果に対応する演算子を ket に作用させ、そのノルムの二乗を結果確率とする測定である。第一結果を得た後にデータを $\ket a$ へ射影すると、制御の非規格化 ket は $(\ket0+B\ket1)/\sqrt2$ になる。その規格化後の $X$ の期待値は $2\Re B/(1+|B|^2)=4/9$、$Y$ は $2\Im B/(1+|B|^2)=4/9$ である。順を交換すると $Y$ は $-4/9$ になる。この測定ではどちらの制御値から第一結果が出たかを記録しないので、二つの振幅が制御ビットで干渉する。
:::

## 多量子ビット系とテンソル積

### 計算基底と係数行列

複数のレジスタを一つの量子系として扱うには、局所基底のすべての組を基底として持つ空間が必要である。:ref[post-qm-composition] に従って、$A,B$ の合成を $\mathcal H_A\otimes\mathcal H_B$ とする。:ref[thm-linalg-tensor-basis] により、$\ket i_A\otimes\ket j_B=\ket{ij}$ は基底をなす。内積は $\braket{u\otimes v}{x\otimes y}=\braket ux\braket vy$ を満たす。

$n$ 量子ビットの場合、計算基底を $\ket{x_1\cdots x_n}$ と書き、$x_1$ を最上位ビットとする。行列を作るときは辞書順に並べる。この順序を変えると同じ演算の行列表現も変わる。

:::proposition[積ベクトルの判定]{#prop-qi-product-rank}
零でない $\ket\psi=\sum_{ij}C_{ij}\ket{ij}$ が積ベクトルであるための必要十分条件は $\rank C=1$ である。
:::

:::proof
$\ket\psi=\ket a\otimes\ket b$ ならば $C_{ij}=a_i b_j$ なので、列空間は $\ket a$ の張る直線となる。逆に階数が 1 ならば、ある非零列 $a$ を基準に各列を $b_j a$ と書ける。したがって元のベクトルは $a\otimes b$ である。単位ベクトルの場合は $\norm a\norm b=1$ なので各因子を規格化できる。
:::

:::example[二つのビットが一致する重ね合わせ]
$\ket{\Phi^+}=(\ket{00}+\ket{11})/\sqrt2$ の係数行列は $I/\sqrt2$ であり、階数は 2 である。したがって単一量子ビットの積には分解できない。$Z$ 基底を両側で測ると $00,11$ が各 $1/2$ であり、$X$ 基底を両側で測ると $\ket{\Phi^+}=(\ket{++}+\ket{--})/\sqrt2$ より $++,--$ が各 $1/2$ である。
:::

$\ket{\Phi^\pm}=(\ket{00}\pm\ket{11})/\sqrt2$、$\ket{\Psi^\pm}=(\ket{01}\pm\ket{10})/\sqrt2$ の四本を :term[Bell 基底]{en="Bell basis" yomi="べるきてい"}という。$\Phi$ と $\Psi$ は計算基底での支持が互いに素である。同じ支持を持つ二本は内積が $(1-1)/2=0$ になる。従って四本は二量子ビット空間の正規直交基底である。

:::example[積状態でも一様な計算基底分布を作れる]
$\ket+\otimes\ket-=(\ket{00}-\ket{01}+\ket{10}-\ket{11})/2$ では全結果が確率 $1/4$ である。係数行列の 2 行は等しく、階数 1 の条件を満たす。したがって、計算基底での確率の形だけから積状態であるかを判定することはできない。
:::

演算 $L\otimes M$ は積ベクトルに $L\ket a\otimes M\ket b$ を対応させる。積基底は空間を張るので、この指定で線形演算全体が決まる。積基底に作用させれば $(L\otimes M)(P\otimes Q)=LP\otimes MQ$ と $(L\otimes M)^\dagger=L^\dagger\otimes M^\dagger$ が従う。

### 部分系の順序と局所演算

テンソル積の結合則は、$\ket{(ij)k}$ と $\ket{i(jk)}$ を対応させる一意なユニタリで与えられる。この対応は直交基底から直交基底への写像だからである。各因子を交換する基本操作は $\mathrm{SWAP}\ket{i}_A\ket{j}_B=\ket{j}_B\ket{i}_A$ と書く。$A,B$ の次元が一致する場合、同じ順序のレジスタの状態を入れ替えるゲートともみなせる。それぞれの周辺演算子を追跡するとき、空間の並べ替えとその空間における状態の変化を区別する必要がある。

:::proposition[両側の係数行列の変換]{#prop-qi-tensor-coefficient-transform}
$\ket\psi=\sum_{ij} C_{ij}\ket{ij}$ に $L\otimes M$ を作用させた後の係数行列は $LCM^\top$ である。$L,M$ が可逆ならば、変換前後で係数行列の階数は等しい。
:::

:::proof
$(a,b)$ 成分を取ると、入力 $\ket{ij}$ の寄与は $L_{ai}M_{bj}C_{ij}$ である。この積を $i,j$ で足すと $(LCM^\top)_{ab}$ となる。両側の因子が可逆なら、行空間と列空間の次元は変わらない。ここで右から掛けるのは転置である。ket の係数を展開すると複素共役は入らない。
:::

:::example[位相を含む局所ゲート]
$\ket\psi=(\ket{00}+2i\ket{11})/\sqrt5$ に対し、第 2 系に $S=\diag(1,i)$ を掛ける。$C=\diag(1,2i)/\sqrt5$ なので $CS^\top=\diag(1,-2)/\sqrt5$ である。$X\otimes X$ は $\ket{00}$ と $\ket{11}$ を交換する。したがって期待値は変換前に $2\Re(2i)/5=0$、変換後に $2\Re(-2)/5=-4/5$ となる。一方、計算基底での確率 $(1/5,0,0,4/5)$ は位相の変化を受けない。
:::

:::example[三量子ビットの二種類の分け方]
$\ket\psi=(\ket{000}+\ket{001}+\ket{110}+\ket{111})/2$ を考える。$AB|C$ で見ると、$AB=00,11$ の両行は $[1/2,1/2]$ で、他の行は零である。したがって係数行列の階数は 1 であり、$(\ket{00}+\ket{11})_{AB}/\sqrt2\otimes\ket+_C$ と因数分解できる。これに対し $A|BC$ では、$A=0$ の行は $[1/2,1/2,0,0]$、$A=1$ の行は $[0,0,1/2,1/2]$ であり、階数は 2 である。積かどうかを尋ねるには、先にどの部分系の間を分けるか指定する。
:::

:::example[見かけの一様性と条件付き位相]
$\ket{\psi}=(\ket{00}+\ket{01}+\ket{10}-\ket{11})/2$ を第 1 ビットの計算基底で分解すると $\ket0\ket+/\sqrt2+\ket1\ket-/\sqrt2$ となる。第 1 ビットを $Z$ で測ってから第 2 ビットを $X$ で測れば、$(0,+1)$ と $(1,-1)$ だけが各 $1/2$ で出る。両方を $Z$ で測る場合は全結果が $1/4$ である。条件付きの第 2 系の二つの状態が互いに直交しているため、係数行列 $[1,1;1,-1]/2$ の階数は 2 である。
:::

### 局所基底での相関の展開

各量子ビットで $I,X,Y,Z$ を選んでテンソル積を取った演算子を Pauli 文字列と呼ぶ。文字列 $P,Q$ のトレース内積は、各因子のトレース内積の積である。したがって $\Tr(PQ)=2^n\delta_{PQ}$ となり、$4^n$ 個の文字列は Hermite 演算子の実基底である。任意の Hermite 演算子 $M$ は
$$
M=2^{-n}\sum_P\Tr(PM)P
$$
と展開できる。$\Tr(PM)$ は実数である。実際、共役は $\Tr(MP)$ になり、巡回性で元に戻る。二体以上の文字列では、各側の測定結果の積を取って平均を得る。単独の測定結果を平均して掛ける手順と一致するのは、その観測量について相関がない場合である。

:::example[三体相関の符号]
$\ket{\psi_\phi}=(\ket{000}+e^{i\phi}\ket{111})/\sqrt2$ では $XXX$ が二項を交換するため $\langle XXX\rangle=\cos\phi$ である。$XXY$ の作用では $\ket{000}$ に $i$、$\ket{111}$ に $-i$ が掛かって交換され、期待値は $\sin\phi$ になる。$XYY$ では交換に $-1$ が掛かり、期待値は $-\cos\phi$ である。一方、一体の $X$ は $\ket{000},\ket{111}$ の張る空間の外へ移すので期待値は零である。例えば $\phi=\pi/6$ で平均は順に $\sqrt3/2,1/2,-\sqrt3/2$ であり、一体の期待値だけからこの三体位相は分からない。
:::

### 交換演算で二つの入力の重なりを測る

等次元の二系で $\mathrm{SWAP}=\sum_{ij}\ketbra ij\otimes\ketbra ji$ である。行列単位をトレースすると、任意の $L,M$ に対して $\Tr[\mathrm{SWAP}(L\otimes M)]=\Tr(LM)$ になる。二つの独立な入力 $\rho,\sigma$ で期待値を取ると $\Tr(\rho\sigma)$ を得る。特に同一状態の二つの独立準備では $\Tr\rho^2$ が求まる。$\mathrm{SWAP}$ は自己共役で二乗が $I$ なので、固有値 $\pm1$ の射影は $(I\pm\mathrm{SWAP})/2$ である。二つの状態ベクトルを入れる場合、$-1$ の枝は反対称な結合 $\ket a\ket b-\ket b\ket a$ に比例する。

:::example[二コピーの純度検査と、相関のある入力]
$\rho=(I+Z/2)/2$ を独立に二つ準備する。$\Tr\rho^2=5/8$ なので、反対称な結果の確率は $(1-5/8)/2=3/16$ になる。直接計算すると、反対称部分は $\Psi^-=(\ket{01}-\ket{10})/\sqrt2$ のみであり、$\rho\otimes\rho$ の $01,10$ の確率は各 $3/16$ だから一致する。一方、同じ縮約をもつ相関状態 $(3/4)\ketbra{00}{00}+(1/4)\ketbra{11}{11}$ を入れると、反対称部分の確率は零である。純度検査で求めているのは一入力だけの性質だが、二コピーの独立性を前提としている。

この検査を補助ビットで行う場合は、補助を $\ket+$ にして制御 SWAP を施し、$X$ を測る。データが純粋 $\ket\chi$ なら全状態は $(\ket0\ket\chi+\ket1\mathrm{SWAP}\ket\chi)/\sqrt2$ であり、$X$ の平均は $\Re\mel\chi{\mathrm{SWAP}}\chi$ である。混合入力では分解について平均すると同じ公式になる。二つのデータを破壊せず結果を記録できるが、結果ごとのデータは対称または反対称部分へ射影される。
:::

## 密度行列と混合状態

### 平均確率を表す演算子

同じ実験を繰り返すとき、準備装置が $\ket{\psi_j}$ を確率 $p_j$ で選び、選んだラベルを記録しない場合を考える。このような確率と状態の組の族を :term[アンサンブル]{en="ensemble" yomi="あんさんぶる"}という。測定結果の確率は各純粋状態の Born 確率の平均である。従って、後から統計を計算するために必要な情報は次に定義する演算子で表せる。

:::definition[密度演算子]{#def-qi-density}
有限次元空間上で $\rho\ge0$ かつ $\Tr\rho=1$ を満たす演算子を :term[密度演算子]{en="density operator" yomi="みつどえんざんし"}と呼ぶ。その基底表示を密度行列という。階数 1 の状態を純粋状態、それ以外を :term[混合状態]{en="mixed state" yomi="こんごうじょうたい"}という。
:::

:::proposition[アンサンブルと状態]{#prop-qi-ensemble}
密度演算子は単位ベクトルの射影の凸結合 $\rho=\sum_j p_j\ketbra{\psi_j}{\psi_j}$ として表される。射影 $P$ の結果確率は $\Tr(P\rho)$ であり、すべての階数 1 射影の結果確率は $\rho$ を一意に定める。
:::

:::proof
凸結合は、任意の $v$ に対して二次形式が $\sum_jp_j\abs{\braket{\psi_j}{v}}^2\ge0$ であり、トレースは $\sum_jp_j=1$ なので密度演算子である。逆に :ref[thm-linalg-spectral] で対角化すれば、その固有値は非負で和が 1 なので求める凸結合になる。確率の式は $\Tr(P\ketbra{\psi_j}{\psi_j})=\mel{\psi_j}{P}{\psi_j}$ と期待値の平均による。二つの密度演算子が同じ確率を与えるとき、差の正規化固有ベクトルへの射影で測れば、対応する固有値は零である。すべての固有値が零なので差は零である。
:::

### Bloch 球の内部と純度

:::proposition[混合状態の Bloch 表示]{#prop-qi-bloch-mixed}
任意の 1 量子ビットの密度行列は $\rho=(I+\boldsymbol r\cdot\boldsymbol\sigma)/2$ と一意に書ける。$\boldsymbol r=(\Tr\rho X,\Tr\rho Y,\Tr\rho Z)$ であり、固有値は $(1\pm\norm{\boldsymbol r})/2$ である。よって状態全体は閉単位球と一対一に対応する。
:::

:::proof
$I,X,Y,Z$ はトレース内積について互いに直交し、それぞれのノルムの二乗は 2 である。Hermite 行列全体の実次元 4 を張るので係数は一意であり、$\Tr\rho=1$ から $I$ の係数は $1/2$ である。式 :eqref[eq-qi-pauli-product] から $(\boldsymbol r\cdot\boldsymbol\sigma)^2=\norm{\boldsymbol r}^2 I$ であり、トレースが零なので固有値は $\pm\norm{\boldsymbol r}$ である。非負性は $\norm{\boldsymbol r}\le1$ と同値である。
:::

$\Tr\rho^2$ を :term[純度]{en="purity" yomi="じゅんど"}と呼ぶ。$d$ 次元で固有値を $\lambda_j$ とすれば $1/d\le\sum_j\lambda_j^2\le1$ である。下限は Cauchy–Schwarz、上限は $0\le\lambda_j\le1$ から出る。上限を取るのは固有値の一つが 1 の場合に限り、下限を取るのは固有値がすべて $1/d$ の場合に限る。量子ビットでは純度は $(1+\norm{\boldsymbol r}^2)/2$ である。

:::example[同じ密度行列を作る二種類の準備]
$\ket0,\ket1$ を半々で準備すると状態は $I/2$ である。$\ket+,\ket-$ を半々で準備しても、射影の和が $I$ なので同じ状態になる。一方、$\ket+$ だけを準備すると $\rho=(I+X)/2$ となり、$X$ の結果は常に $+1$ である。
:::

:::example[偏った混合]
$\rho=\frac34\ketbra00+\frac14\ketbra11$ は $\boldsymbol r=(0,0,1/2)$ を持ち、純度は $9/16+1/16=5/8$ である。$Z$ 測定の確率は $3/4,1/4$ であるが、$X$ 測定の確率は $1/2,1/2$ となる。
:::

### 凸集合の端点と三つの純度指標

密度行列を「準備装置が古典乱数で選べる状態」の集合として考えると、行列の階数から、さらに分解できるかどうかを決められる。凸集合の点 $x$ が、異なる二点 $y,z$ を $0<t<1$ で混ぜた $ty+(1-t)z$ として書けないとき、$x$ を端点という。

:::proposition[状態集合の端点]{#prop-qi-density-extreme}
密度演算子の集合の端点は、階数 1 の射影に一致する。
:::

:::proof
$\rho=\ketbra\psi\psi=t\sigma+(1-t)\tau$ とする。$\psi$ に直交する任意の $v$ に対し、$\mel v\rho v=0$ は二つの非負項の正の重み付き和である。よって $\mel v\sigma v=\mel v\tau v=0$ である。正行列 $D$ に対し $\mel vDv=\norm{\sqrt D v}^2=0$ なら $Dv=0$ だから、$\sigma,\tau$ の支持は $\psi$ の直線に含まれる。両者のトレースは 1 なので $\sigma=\tau=\rho$ である。

逆に正固有値 $\lambda_1,\lambda_2$ を二つ持つなら、$0<\eps<\min(\lambda_1,\lambda_2)$ を選び、対応する二つの射影を $P_1,P_2$ とおく。$\rho\pm\eps(P_1-P_2)$ は共に密度演算子であり、二つは異なり、その中点は $\rho$ である。
:::

:::example[量子ビットの純度と行列式]
Bloch ベクトル $(1/3,2/3,2/3)$ の長さは 1 なので、その固有値は 1 と 0 である。一方、長さを半分にした状態は固有値 $3/4,1/4$ を持ち、$\Tr\rho^2=5/8$、$\det\rho=3/16$ である。一般の $2\times2$ 状態に対して $(\Tr\rho)^2-\Tr\rho^2=2\det\rho$ だから、行列式から純度を再計算できる。対角成分が $(2/3,1/3)$ の場合でも、非対角成分が零なら純度は $5/9$ になる。計算基底の頻度だけを揃えても、行列の純度は揃わない。
:::

:::example[同じ純度で異なる階数]
三準位系の $\rho=\diag(1/2,1/2,0)$ と $\sigma=\diag(2/3,1/6,1/6)$ を比較する。純度はいずれも $1/2$ である。$\rho$ は第 3 基底状態への射影で結果が決して出ないが、$\sigma$ はその結果が確率 $1/6$ で出る。したがって純度から階数を決定できるのは、純粋性を判定する上限や、次元と等分布の関係を使う場合などに限られる。
:::

### 非直交アンサンブルの固有値を小さい行列から求める

アンサンブル $(p_j,\psi_j)$ から、列を $\sqrt{p_j}\ket{\psi_j}$ とする行列 $W$ を作ると、$\rho=WW^\dagger$ である。$G=W^\dagger W$ を重み付き Gram 行列と呼ぶ。$Gv=\lambda v$ で $\lambda>0$ なら $Wv$ は零でなく、$\rho Wv=WGv=\lambda Wv$ である。逆に $W^\dagger$ を適用すれば $\rho$ の正固有空間から $G$ の固有空間へ移る。二つの写像の合成は $\lambda I$ だから、正固有値は重複度まで同じである。

:::example[三次元にある二本の信号]
三準位系で $\ket0$ を確率 $1/3$、$(\ket0+2\ket2)/\sqrt5$ を確率 $2/3$ で用意する。内積は $1/\sqrt5$ なので、
$$
G=
\begin{pmatrix}
1/3 & \sqrt{2/5}/3\\
\sqrt{2/5}/3 & 2/3
\end{pmatrix}.
$$
トレースは 1、行列式は $8/45$ である。したがって $\rho$ の固有値は $0,(1+\sqrt{13/45})/2,(1-\sqrt{13/45})/2$ である。純度は $1-16/45=29/45$ で、二つの確率だけから計算した $5/9$ より大きい。内積の二乗 $1/5$ に交差項の係数 $2(1/3)(2/3)$ を掛けると差 $4/45$ を得る。
:::

### 正性を行列要素から検査する

任意の密度行列 $\rho$ に対し $|\rho_{ij}|^2\le\rho_{ii}\rho_{jj}$ が成り立つ。$\rho_{ij}$ は $\sqrt\rho\ket i$ と $\sqrt\rho\ket j$ の内積なので Cauchy–Schwarz で出る。対角成分が零なら、同じ行と列はすべて零である。二次元の場合は、全対角成分の非負性とこの不等式が正性の十分条件でもある。固有値の和は非負で積は非負なので、負の固有値があるなら二つとも負になるか、もう一方が零になるが、どちらも和の条件に反する。高次元では二項ずつの条件だけでは足りない。

:::example[二項ずつ正でも全体では正でない]
三準位で $\rho_{ii}=1/3$、$i\ne j$ で $\rho_{ij}=-1/4$ とする。各二次主小行列の固有値は $1/3\pm1/4=7/12,1/12$ であり、ともに正である。ところが $\ket s=(\ket0+\ket1+\ket2)/\sqrt3$ には行和 $1/3-1/2=-1/6$ が固有値として対応する。行列を $(7/12)I-(1/4)\mathbf1\mathbf1^\top$ と書けば、$\ket s$ に直交する二次元では固有値は $7/12$ である。従ってこの擬似状態は各二次元の部分空間での基底検査に通っても、$\ketbra ss$ の測定確率を負にするため物理状態ではない。

同じ対角成分で非対角を $-t$ と置くと、$t\ge0$ では全体の正性が $t\le1/6$ となる。境界の $t=1/6$ では固有値 $(1/2,1/2,0)$ であり、$\ket s$ を除いた二次元の一様状態である。正性を回復するため一様状態を元の行列に確率 $q$ で混ぜるなら、非対角は $-(1-q)/4$ となるので $q\ge1/3$ が必要である。
:::

## 部分トレースと縮約密度行列

### 局所統計の取り出し方

$A,B$ を合わせた状態を知っていても、手元に $A$ しかなければ $A$ の測定統計だけを計算したい。このとき $B$ で何も選別しないので、$A$ の効果には $B$ の恒等演算子を掛ける。

:::definition[部分トレース]{#def-qi-partial-trace}
$T$ を $A\otimes B$ 上の演算子とする。$B$ の正規直交基底 $(\ket j)$ を用いて :term[部分トレース]{en="partial trace" yomi="ぶぶんとれーす"}を

$$

\Tr_B T=\sum_j(I_A\otimes\bra j)T(I_A\otimes\ket j)
\label{eq-qi-partial-trace}
$$

と定義する。$\rho_A=\Tr_B\rho_{AB}$ を $A$ の :term[縮約密度行列]{en="reduced density matrix" yomi="しゅくやくみつどぎょうれつ"}と呼ぶ。
:::

:::proposition[部分トレースの特徴付け]{#prop-qi-partial-trace}
任意の $A$ 上の演算子 $L$ について $\Tr(L\Tr_B T)=\Tr((L\otimes I_B)T)$ であり、この性質が $\Tr_B T$ を一意に定める。特に部分トレースは基底によらず、密度演算子を密度演算子に写す。
:::

:::proof
積基底でトレースを展開すると、左辺は $\sum_{ij}\mel{ij}{(L\otimes I_B)T}{ij}$ となって右辺に等しい。任意の $L$ に対して $\Tr(LD)=0$ となる差 $D$ には $L=D^\dagger$ を選べば $\sum_{ij}\abs{D_{ij}}^2=0$ が従い、一意性を得る。$\rho_{AB}\ge0$ なら各項の二次形式は $\mel{v j}{\rho_{AB}}{v j}\ge0$ である。$L=I_A$ と取ればトレースも保存する。
:::

:::example[積基底の非対角成分]
$\Tr_B\ketbra{ij}{kl}=\delta_{jl}\ketbra ik$ である。したがって $\Tr_B\ketbra{\Phi^+}{\Phi^+}=I/2$ となる。全体が純粋でも局所状態は混合になりうるので、局所の混合を直ちに準備装置の古典的な乱択へ帰着することはできない。
:::

:::example[同じ周辺状態と異なる相関]
$\rho_{AB}=(\ketbra{00}{00}+\ketbra{11}{11})/2$ と Bell 状態はともに $A,B$ の縮約が $I/2$ である。ところが $\Tr[(X\otimes X)\rho_{AB}]=0$ であり、Bell 状態では $(X\otimes X)\ket{\Phi^+}=\ket{\Phi^+}$ より期待値が 1 である。局所測定の記録を照合すれば両者を区別できる。
:::

$B$ だけにユニタリ $U$ を作用させても $A$ は変化しない。実際、$I\otimes U$ と $L\otimes I$ が可換であり、:ref[prop-qi-partial-trace] の右辺はトレースの巡回性により変わらない。$B$ の測定結果で条件付ければ $A$ の状態は変化しうる。この場合、条件付けに用いる結果を $A$ 側に別途伝える必要がある。

### 条件付き状態と無条件状態

$B$ の各測定結果 $b$ には、結果確率を与える正演算子 $E_b$ が対応する。$p_b=\Tr(E_b\rho_B)$ であり、$\sum_b E_b=I_B$ とする。このような族を POVM と呼ぶ。以下では、測定後の状態まで指定するため、結果 $b$ に演算子 $M_{b,j}$ の族を対応させる。添字 $j$ は装置が区別して記録しない枝である。このとき条件 $b$ における $A$ の非規格化状態は
$$
\tau_A(b)=\Tr_B\sum_j (I\otimes M_{b,j})\rho_{AB}(I\otimes M_{b,j}^\dagger)
$$
である。$\Tr_A\tau_A(b)$ が条件の確率であり、その値が正なら規格化する。任意の $L$ で $\Tr[L\tau_A(b)]=\Tr[(L\otimes E_b)\rho_{AB}]$ となる。これは :ref[prop-qi-partial-trace] とトレースの巡回性から従う。ただし $E_b=\sum_jM_{b,j}^\dagger M_{b,j}$ である。$B$ の測定後出力を捨てるので、$\tau_A(b)$ を求めるには $M_{b,j}$ の個別の情報を省ける。また $\sum_b E_b=I$ が成り立てば $\sum_b\tau_A(b)=\rho_A$ である。条件付きの局所状態と、記録を捨てた周辺状態が異なる理由はこの式に表れている。

:::example[偏ったもつれから条件付き状態を作る]
$\ket\psi=(2\ket{00}+\ket{11})/\sqrt5$、$B$ の効果を $E_1=(I+X/2)/2$、$E_2=I-E_1$ とする。$\tau_A(b)$ の $(i,j)$ 成分は $\lambda_i^{1/2}\lambda_j^{1/2}(E_b)_{ji}$ であり、$(\lambda_0,\lambda_1)=(4/5,1/5)$ である。これを二つの効果に代入すると
$$
\tau_A(1)=\frac1{10}\begin{pmatrix}4&1\\1&1\end{pmatrix},
\quad
\tau_A(2)=\frac1{10}\begin{pmatrix}4&-1\\-1&1\end{pmatrix}
$$
となる。どちらも確率 $1/2$ であり、条件付き状態の $X$ の期待値は $\pm2/5$ である。平均では非対角成分が消える。$B$ の結果を知らずに $A$ を測る限り、$X$ の結果は半々である。
:::

:::example[縮約の順序と相関]
$\ket\psi=(\ket{000}+2\ket{011}+2i\ket{101})/3$ を考える。$BC$ の基底ラベル $00,11,01$ は互いに異なるので、$\rho_A=\diag(5/9,4/9)$ になる。先に $C$ を捨てると、
$$
\rho_{AB}=\tfrac19\ketbra{00}{00}
+\tfrac49(\ket{01}+i\ket{10})(\bra{01}-i\bra{10})
$$
である。この式で $\Tr_B$ を取っても同じ $\rho_A$ が得られる。途中の交差項により $\langle X\otimes Y\rangle_{AB}=-8/9$ となる。その符号は $X\otimes Y$ が $\ket{01}$ に $-i\ket{10}$ を与えることから決まる。先に $B$ まで捨てると、この相関は復元できない。
:::

### 周辺状態の適合性

別々に指定した縮約が同じ全体系から得られるかどうかには制約がある。例えば $\rho_A$ が純粋なら $\rho_{AB}=\rho_A\otimes\rho_B$ である。これを示すため、$\rho_A=\ketbra aa$、$v\perp a$ とする。$B$ の基底 $(\ket j)$ に沿って $\sum_j\mel{v j}{\rho_{AB}}{v j}=0$ である。各項は非負なので全て零となり、$\rho_{AB}$ は $v\otimes\ket j$ を消す。したがって支持は $\operatorname{span}\{a\}\otimes\mathcal H_B$ に含まれ、$\rho_{AB}$ は $\ketbra aa\otimes T$ の形となる。$A$ の部分トレースを取ると $T=\rho_B$ である。この議論は複合系 $AB$ 自体を一つの部分系に取っても成り立つ。

:::example[二つの Bell 縮約を同時に指定できるか]
$\rho_{AB}=\ketbra{\Phi^+}{\Phi^+}$ とする。このとき任意の拡張 $\rho_{ABC}$ は $\rho_{AB}\otimes\rho_C$ と書けるので、$\rho_{AC}=I_A/2\otimes\rho_C$ である。仮に $\rho_C=I/2$ なら $\rho_{AC}=I/4$ であり、$\langle X_A X_C\rangle=0$ である。$AC$ も $\Phi^+$ にする指定なら相関は 1 になるので両立しない。ここでは複合系 $AB$ の純粋性から積分解を得た。
:::

### 異なる全体が同じ局所測定を再現する

二体系で両方の縮約を指定しても、一般に全体は決まらない。$\rho_A\otimes\rho_B$ は常に一つの解である。また両者が正定値なら小さい $\eps$ に対して、跡零の Hermite 演算子 $L,M$ を用いた $\rho_A\otimes\rho_B+\eps L\otimes M$ も解になる。部分トレースで追加項は零であり、最小固有値 $\lambda_{\min}(\rho_A)\lambda_{\min}(\rho_B)$ より $\abs\eps\norm{L\otimes M}_\infty$ が小さければ、全体は正である。和の跡は変わらない。従って全体の未知パラメータには、局所記録をどれだけ増やしても決定できない方向がある。

:::example[局所状態を保ったまま相関を動かす範囲]
$\rho_A=\diag(3/4,1/4),\rho_B=\diag(2/3,1/3)$ とし、$\rho_{AB}=\rho_A\otimes\rho_B+t Z\otimes Z$ とする。対角成分は $1/2+t,1/4-t,1/6-t,1/12+t$ である。各行を足すと $3/4,1/4$、各列を足すと $2/3,1/3$ となり、$t$ は見えない。正性の範囲は $-1/12\le t\le1/6$ である。積の $ZZ$ 期待値は $(1/2)(1/3)=1/6$ だが、この全体状態の値は $1/6+4t$ となり、$-1/6$ から $5/6$ まで変わる。例えば $t=1/6$ では全体のラベル $10$ が消えるが、各ビットの確率は依然正である。この族は全て同じ基底に対角なので、相関の変化を実現するのに純粋な重ね合わせは要らない。
:::

## 純粋化

### 補助系の構成

混合状態を含む式でも、対象と補助系を合わせて純粋状態にできれば、ベクトルに対する議論が使える。部分トレースが直交したラベル間の交差項を消すことが、その構成に役立つ。

:::definition[純粋化]{#def-qi-purification}
$\Tr_R\ketbra{\Psi}{\Psi}=\rho_A$ を満たす単位ベクトル $\ket\Psi\in\mathcal H_A\otimes\mathcal H_R$ を $\rho_A$ の :term[純粋化]{en="purification" yomi="じゅんすいか"}と呼ぶ。$R$ を参照系という。
:::

:::theorem[純粋化の存在と自由度]{#thm-qi-purification}
$\rho_A=\sum_{j=1}^r\lambda_j\ketbra{u_j}{u_j}$ を正固有値に限ったスペクトル分解とする。純粋化は $\dim R\ge r$ のとき存在する。任意の純粋化は $\sum_j\sqrt{\lambda_j}\ket{u_j}\ket{v_j}$ と書け、$v_j$ は $R$ の正規直交系である。共通の参照空間上の任意の二つの純粋化は $R$ 上のユニタリで結ばれる。
:::

:::proof
正規直交系 $(\ket{v_j})$ を選び、定理のベクトルを作る。外積を取り式 :eqref[eq-qi-partial-trace] を適用すると、$\braket{v_k}{v_j}=\delta_{jk}$ により所望の $\rho_A$ が得られる。逆に、任意の純粋化を $\sum_j\ket{u_j}\ket{w_j}$ と展開する。ここでは零固有値の固有ベクトルも含める。その縮約の $(j,k)$ 成分は $\braket{w_k}{w_j}=\lambda_j\delta_{jk}$ なので、零固有値の成分は消え、$\ket{v_j}=\ket{w_j}/\sqrt{\lambda_j}$ は正規直交系となる。従って $\dim R\ge r$ が必要である。二つの正規直交系をそれぞれ基底に延長し、対応する基底を互いに移せば所望のユニタリを得る。
:::

:::example[偏った混合の参照系]
$\sqrt{3/4}\ket{00}+\sqrt{1/4}\ket{11}$ は $\diag(3/4,1/4)$ の純粋化である。第 2 量子ビットに $X$ を作用させた $\sqrt{3/4}\ket{01}+\sqrt{1/4}\ket{10}$ も同じ縮約を与える。階数 2 であるから参照系に使う量子ビットを省くことはできない。
:::

:::example[古典的なラベルを測る]
$\ket\Psi=(\ket0\ket0+\ket+\ket1)/\sqrt2$ の第 2 系を計算基底で測ると、第 1 系に $\ket0$ と $\ket+$ が各 $1/2$ で準備される。測定前の $A$ の密度行列は $\frac12(\ketbra00+\ketbra++)$ であり、Bloch ベクトルは $(1/2,0,1/2)$ である。アンサンブルの二つの状態が直交していなくても、ラベルを参照系の直交基底に記録できる。
:::


### アンサンブルを参照系の測定で選ぶ

:::proposition[純粋化とアンサンブルの関係]{#prop-qi-purification-ensemble}
$\rho=\sum_{j=1}^r\lambda_j\ketbra{u_j}{u_j}$ を正固有値に限った分解とし、$r$ 次元の参照系で $\ket\Psi=\sum_j\sqrt{\lambda_j}\ket{u_j}\ket j$ を固定する。任意の有限純粋アンサンブル $\rho=\sum_a p_a\ketbra{\psi_a}{\psi_a}$ は、この純粋化の参照系に対する適切な測定から得られる。
:::

:::proof
$p_a>0$ の項だけを取る。$\psi_a$ は $\rho$ の支持に属する。実際、核の任意のベクトルへの二次形式は、各非負項の和が零になる。参照系で $\braket{w_a}{j}=\sqrt{p_a}\braket{u_j}{\psi_a}/\sqrt{\lambda_j}$ と定める。アンサンブルの等式の $(i,j)$ 成分を比べれば $\sum_a\ketbra{w_a}{w_a}=I_R$ である。従ってこれらは効果であり、$(I\otimes\bra{w_a})\ket\Psi=\sqrt{p_a}\ket{\psi_a}$ である。参照系の測定には、例えば固定した単位ベクトル $\ket v$ を用いて測定演算子 $M_a=\ketbra v{w_a}$ を選べる。このとき主系に所望の非規格化状態が残る。より大きな参照系に対しては、定理 :ref[thm-qi-purification] の正規直交系が張る部分空間で同じ効果を使い、その補空間の射影を確率零の結果として加えればよい。
:::

:::example[同じ純粋化から異なるアンサンブルを得る]
$\Phi^+$ の参照ビットを $Z$ で測ると $\ket0,\ket1$ が半々で得られ、$X$ で測ると $\ket+,\ket-$ が半々で得られる。参照ビットを $Y$ で測ると、固有値 $+1$ に対して $(\ket0-i\ket1)/\sqrt2$ が出る。標準純粋化で $\bra{y_+}$ を掛けるので、第 1 系の係数は測定ベクトルの複素共役になる。この符号は実基底の例だけでは現れない。
:::

### ベクトル化と参照系の転置

$A$ の基底を固定し、同じ次元の参照系を取って $\ket\Omega=\sum_j\ket{jj}$ とおく。これはノルム $\sqrt{\dim A}$ の非規格化ベクトルである。行列 $K$ に $\operatorname{vec}(K)=(K\otimes I)\ket\Omega$ を対応させる。$P_K$ をこのベクトルの非規格化外積とおくと、
$$
\begin{aligned}
\braket{\operatorname{vec}(K)}{\operatorname{vec}(L)}&=\Tr(K^\dagger L),\\
\Tr_R P_K&=KK^\dagger,\qquad
\Tr_A P_K=(K^\dagger K)^\top,\\
(I\otimes M)\operatorname{vec}(K)&=\operatorname{vec}(KM^\top)
\end{aligned}
$$
となる。内積は $\sum_{ij}\overline{K_{ij}}L_{ij}$ である。二つの部分トレースは、それぞれ $\sum_j K_{ij}\overline{K_{i'j}}$、$\sum_i K_{ij}\overline{K_{ij'}}$ を成分として持つ。最後の式は $(i,j)$ 成分を取ると $\sum_l K_{il}M_{jl}$ であることから従う。特に $K=\sqrt\rho$ を取れば規格化した純粋化が得られ、参照系は $\rho^\top$ を持つ。

:::example[複素成分のある標準純粋化]
$\rho=(I+3Y/5)/2$ を考える。$Y$ の固有値 $\pm1$ に対する $\rho$ の固有値は $4/5,1/5$ である。射影 $(I\pm Y)/2$ を用いると、正の平方根は $\sqrt\rho=(3I+Y)/(2\sqrt5)$ と書ける。したがって規格化した純粋化は
$$
(3\ket{00}-i\ket{01}+i\ket{10}+3\ket{11})/(2\sqrt5)
$$
である。第 1 系の $Y$ の期待値は $3/5$ だが、第 2 系では $-3/5$ となる。実際、$Y^\top=-Y$ だからである。第 2 系で $+1$ を得る確率は $1/5$ であり、第 1 系の対応する条件付き状態は $Y$ の $-1$ 固有状態になる。
:::

:::example[最小次元と共通の拡張]
三準位の $\rho=\diag(1/3,0,2/3)$ には、$\ket\Psi=(\ket{0,0}+\sqrt2\ket{2,1})/\sqrt3$ が二次元参照系での純粋化となる。三次元参照系を許すなら、$\ket\Phi=(\ket{0,2}+\sqrt2\ket{2,0})/\sqrt3$ と取ることもできる。先の二次元系を $\ket0,\ket1$ の張る空間に埋め込む。三次元系上の基底の置換 $0\mapsto2,1\mapsto0,2\mapsto1$ を参照系に作用させると $\Phi$ が得られる。小さい参照系から大きい参照系への写像は等長であり、余った方向を含む共通の空間を指定すればユニタリで表せる。
:::

### 異なる個数の準備ラベルを同じ状態に写す

列行列 $W$ が $WW^\dagger=\rho$ を満たすとする。$\rho$ の支持への逆平方根を $\rho^{-1/2}$、支持外で零と定める。$V=\rho^{-1/2}W$ とおくと $VV^\dagger=P$ であり、$P$ は支持射影である。また $\rho$ の核方向には $WW^\dagger$ の二次形式が零なので、$W$ の各列は支持に入る。よって $W=\sqrt\rho V$ である。逆に $VV^\dagger=P$ なら $WW^\dagger=\rho$ が従う。この因数分解により、正規直交する参照基底をどの準備ラベルに合わせるかという自由度と、参照系で非射影測定を選ぶ自由度を一つの式で扱える。

:::example[二準位参照系から三つの非直交信号]
$\omega=e^{2\pi i/3}$ とおき、$\ket{u_j}=(\ket0+\omega^j\ket1)/\sqrt2$、$j=0,1,2$ とする。$\sum_j\omega^j=0$ なので $\frac13\sum_j\ketbra{u_j}{u_j}=I/2$ である。$\Phi^+$ の参照系を効果 $E_j=(2/3)\ketbra{\overline{u_j}}{\overline{u_j}}$ で測ると、無条件の和は $I$ になり、各分枝で $A$ の非規格化状態は $E_j^\top/2=\ketbra{u_j}{u_j}/3$ である。参照系は二次元だが、結果数が三つの一般化測定を用いることで、三状態を確率 $1/3$ ずつ準備できる。
:::

### 参照系の測定を有限な行列として解く

標準純粋化の行列を $K=\sqrt\rho$ とすると、参照側の効果 $E_a$ に対応する主系の非規格化状態は $K E_a^\top K^\dagger$ である。基底成分で、$(i,j)$ 成分は $\sum_{kl}K_{ik}\overline{K_{jl}}(E_a)_{lk}$ になるからである。従って $\rho$ が支持で可逆なら、任意の正な分解 $\sum_a\tau_a=\rho$ に対し、その支持上で $E_a=(\rho^{-1/2}\tau_a\rho^{-1/2})^\top$ と逆算できる。転置は正性と和を保つので POVM になる。この形は、純粋なアンサンブルだけでなく、条件付き状態 $\tau_a/\Tr\tau_a$ が混合である場合にも使える。

:::example[混合した二つの条件付き状態]
$\rho=\diag(2/3,1/3)$ とし、$\tau_\pm=\rho/2\pm X/12$ を目標にする。各 $\tau_\pm$ の対角成分は $1/3,1/6$、行列式は $1/18-1/144=7/144$ なので正で、跡は $1/2$ である。参照側で必要な効果の非対角成分は $\pm (1/12)/\sqrt{(2/3)(1/3)}=\pm1/(4\sqrt2)$ である。従って $E_\pm=I/2\pm X/(4\sqrt2)$ となる。各効果は全階数であり、どちらの分枝でも主系は混合のままである。規格化した後は Bloch ベクトル $(\pm1/3,0,1/3)$ を持つ。平方根により $K$ で挟む操作は、参照側と同じ長さの Bloch ベクトルをそのまま転写する操作ではない。
:::

## 一般化測定（POVM）と Naimark の拡張定理

### 結果確率と状態更新

量子測定の基本は :ref[def-aqm-povm]{full} で導入した。この節では出力の統計を指定した後で、測定後の状態をどこまで選べるかを調べる。例えば同じ確率でラベルを返す装置でも、入力を保持する場合と所定の状態に準備し直す場合とでは、次の測定の統計が変わる。

:::definition[効果と POVM]{#def-qi-povm}
$0\le E\le I$ を満たす演算子を :term[効果]{en="effect" yomi="こうか"}という。正演算子の有限族 $(E_a)$ が $\sum_aE_a=I$ を満たすとき、これを :term[正作用素値測度]{en="positive operator-valued measure" yomi="せいさようそちそくど"}（POVM）という。状態 $\rho$ における結果 $a$ の確率は $p_a=\Tr(E_a\rho)$ である。
:::

非負性は $\Tr(E_a\rho)=\Tr(\sqrt{E_a}\rho\sqrt{E_a})\ge0$ から、規格化は $\sum_aE_a=I$ から従う。射影測定は各 $E_a$ が互いに直交する射影である場合を指す。

結果 $a$ に対応する測定演算子を $M_{a,\alpha}$ とすると、出力の非規格化状態は $\sum_\alpha M_{a,\alpha}\rho M_{a,\alpha}^\dagger$ である。トレースの巡回性より $E_a=\sum_\alpha M_{a,\alpha}^\dagger M_{a,\alpha}$ となる。条件付き状態を得るには $p_a>0$ のときに $p_a$ で割る。例えば $(\sqrt{E_a})_a$ と $(U_a\sqrt{E_a})_a$ は同じ POVM を持つが、結果に応じたユニタリ $U_a$ は測定後の状態を変える。

:::theorem[Naimark の拡張定理]{#thm-qi-naimark}
有限個の結果を持つ POVM $(E_a)$ に対し、等長写像 $V$ と直交射影の分解 $(P_a)$ が存在し、$E_a=V^\dagger P_a V$ が成り立つ。
:::

:::proof
結果を記録する空間 $\mathcal K$ に正規直交基底 $(\ket a)$ を取り、$V=\sum_a\sqrt{E_a}\otimes\ket a$ とする。すると $V^\dagger V=\sum_aE_a=I$ である。$P_a=I\otimes\ketbra aa$ とすれば $V^\dagger P_aV=E_a$ となる。物理的には、任意に固定した $\ket0_{\mathcal K}$ に対して、部分空間 $\mathcal H\otimes\ket0$ 上で $v\otimes\ket0\mapsto Vv$ とする。定義域と像の正規直交基底を全空間の基底に延長すれば、対応する基底を移すユニタリ $U$ が得られる。補助系を $\ket0$ に準備し、$U$ を適用してから補助系を結果基底で測れば、確率は $\Tr(E_a\rho)$ であり、系に残る非規格化状態は $\sqrt{E_a}\rho\sqrt{E_a}$ となる。
:::

:::example[三結果を持つ量子ビット測定]
$\omega=e^{2\pi i/3}$ とし、$\ket{t_a}=(\ket0+\omega^a\ket1)/\sqrt2$、$a=0,1,2$ とおく。$\sum_a\omega^a=0$ より $\sum_a\ketbra{t_a}{t_a}=3I/2$ なので、$E_a=(2/3)\ketbra{t_a}{t_a}$ は POVM となる。入力 $\ket{t_0}$ では $p_0=2/3$ で、$\abs{1+\omega}=\abs{1+\omega^2}=1$ より $p_1=p_2=1/6$ である。各効果の非零固有値は $2/3$ なので、効果自体は射影ではない。
:::

:::example[確率的に情報を得ない測定]
$E_0=I/4,E_1=3I/4$ とする。結果分布は任意の $\rho$ に対して $(1/4,3/4)$ である。測定演算子を $M_0=I/2,M_1=(\sqrt3/2)X$ と選ぶと、結果 0 では状態は保たれ、結果 1 では $X\rho X$ になる。結果の統計だけでは、このビット反転の有無は分からない。
:::

### 状態更新も指定した測定の実装

一つの効果 $E_a$ に複数の測定演算子 $M_{a,\alpha}$ が対応する場合も、装置の補助空間でまとめて実装できる。補助の直交基底を $\ket{a,\alpha}$ とし、

$$
V\ket\psi=\sum_{a,\alpha}M_{a,\alpha}\ket\psi\otimes\ket{a,\alpha}
$$

とおく。$\sum_{a,\alpha}M_{a,\alpha}^\dagger M_{a,\alpha}=I$ なので $V^\dagger V=I$ である。補助を固定した入力部分空間から $V$ の像への等長写像を、Naimark の証明と同じ基底延長によりユニタリにできる。補助系で $\sum_\alpha\ketbra{a,\alpha}{a,\alpha}$ を測り、補助を部分トレースすると $\alpha$ の交差項が消え、系に $\sum_\alpha M_{a,\alpha}\rho M_{a,\alpha}^\dagger$ が残る。したがって、結果を粗く記録しても、補助内に残った細かいラベルは状態更新に影響する。

:::example[同じ一結果から異なる後状態を得る]
二つの演算子を $M_{0,0}=\ketbra00$、$M_{0,1}=\ketbra11$ として、装置は常に結果 $a=0$ を返す。このとき $E_0=I$ である。入力 $(\ket0+\ket1)/\sqrt2$ に $V$ を掛けると $(\ket0\ket{0,0}+\ket1\ket{0,1})/\sqrt2$ となり、系の縮約は $I/2$ になる。同じ $E_0=I$ を演算子 $M_0=I$ 一つで実装すれば、入力は純粋なまま残る。
:::

### 拡張空間の最小次元

:::proposition[有限 POVM の最小拡張]{#prop-qi-naimark-min}
$E_a=V^\dagger P_a V$ を満たす等長 $V$ と射影分解 $(P_a)$ の拡張空間の最小次元は $\sum_a\rank E_a$ である。
:::

:::proof
任意の拡張で $E_a=(P_aV)^\dagger(P_aV)$ なので、$\rank E_a\le\rank P_a$ である。後者を足せば空間の次元となる。逆に空間を $\bigoplus_a \operatorname{supp}E_a$ と取り、$Vv=(\sqrt{E_a}v)_a$ とする。$\sum_a E_a=I$ より等長であり、各成分への射影の圧縮は $E_a$ となる。
:::

:::example[三結果測定の補助次元]
前の三結果測定では各効果の階数が 1 なので、全拡張空間の最小次元は 3 である。一方、元の量子ビットをそのまま保持するテンソル積実装で補助レジスタに三つのラベルを記録すると、全次元は $2\cdot3=6$ になる。両者は拡張空間に要求した形が異なる。
:::

### 二値の非鋭い測定と後状態

単位ベクトル $\boldsymbol n$ と強さ $\eta\in[0,1]$ を選び、$E_\pm=(I\pm\eta\boldsymbol n\cdot\boldsymbol\sigma)/2$ とおく。この測定の結果の平均は $\eta\boldsymbol n\cdot\boldsymbol r$ になる。$\eta=1$ では射影測定であり、$\eta=0$ では状態とは独立なコインである。平方根を測定演算子に用いれば、測定方向に並んだ固有成分は測定前後で追跡しやすい。

:::proposition[二値測定が無条件に残すコヒーレンス]{#prop-qi-unsharp-backaction}
$\boldsymbol n$ 方向の固有基底において、測定を無視した更新 $\sum_s\sqrt{E_s}\rho\sqrt{E_s}$ は対角成分を保ち、非対角成分を $\sqrt{1-\eta^2}$ 倍する。
:::

:::proof
$E_+$ の対角成分を $a=(1+\eta)/2,b=(1-\eta)/2$ とおく。$E_-$ では $a,b$ が入れ替わる。各測定演算子はその成分の平方根を対角に持つので、対角成分の倍率は $a+b=1$ である。上の非対角成分の倍率は $\sqrt{ab}+\sqrt{ba}=2\sqrt{ab}=\sqrt{1-\eta^2}$ となる。下の非対角成分でも同じ倍率である。
:::

:::example[同じ結果確率、異なる再測定]
$E_+=\diag(4/5,1/5)$ と $E_-=I-E_+$ で $\ket+$ を測る。平方根による分枝では $+$ が確率 $1/2$ で出て、後状態は $(2\ket0+\ket1)/\sqrt5$ である。直後に同じ測定を行うと $+$ の条件付き確率は $(4/5)^2+(1/5)^2=17/25$ になる。これに対し、$+$ の分枝でその後に $X$ を掛ける装置は第一結果の確率を保つが、第二結果の条件付き確率は $2(4/5)(1/5)=8/25$ になる。したがって一度の結果だけを調べて装置の擾乱を定めることはできない。
:::

:::example[補助ビットへの実装]
$E_+=\diag(9/10,1/10)$ とする。補助ビットを $\ket0_R$ に準備し、データの基底値 $x$ によって条件付けた回転を補助系へ掛ける。$x=0$ なら $\ket0_R$ を $(3\ket0_R+\ket1_R)/\sqrt{10}$、$x=1$ なら $(\ket0_R+3\ket1_R)/\sqrt{10}$ に移す。各写像は規格化ベクトルを単位行列の一列に取る $2\times2$ 回転で実装できる。得られた $\ket0_R,\ket1_R$ のラベルを読むと、演算子はそれぞれ $\diag(3,1)/\sqrt{10}$ と $\diag(1,3)/\sqrt{10}$ になる。$\ket+$ を入力にして結果を捨てれば、非対角成分は $1/2$ から $3/10$ に減る。したがって $X$ を測ったときの $+1$ 確率は $1/2+3/10=4/5$ である。
:::

### 粗視化、測定設定の乱択、共同測定

測定結果 $a$ を $q(b|a)$ という条件付き確率で $b$ に変換すると、新しい効果は $F_b=\sum_a q(b|a)E_a$ となる。各 $F_b$ は正であり、$\sum_b q(b|a)=1$ から総和は $I$ である。逆に、効果を分解して細かいラベルを作る場合、その細分は確率だけでなく、どの測定演算子を保持したかにも依存する。例えば $I$ という一結果の効果は $\ketbra00+\ketbra11$ に分解できるが、この分解を測定し結果を捨てると非対角項を消せる。一結果の装置が単に何もしないなら非対角項は残る。

:::proposition[射影測定と共同測定]{#prop-qi-sharp-joint}
射影測定 $(P_a)$ と POVM $(E_b)$ を同じ POVM の周辺として実現できる必要十分条件は、すべての $a,b$ で $P_aE_b=E_bP_a$ が成り立つことである。
:::

:::proof
共同効果を $G_{ab}$ とすると $0\le G_{ab}\le P_a$ である。従って $P_a$ の核の上で $G_{ab}$ の二次形式は零であり、$\sqrt{G_{ab}}$ が核を消す。よって $G_{ab}=P_aG_{ab}P_a$ である。直交する射影 $P_{a'}$ の支持を使うと、$P_a E_b=P_a\sum_{a'}G_{a'b}=G_{ab}=E_bP_a$ が従う。逆に可換なら $G_{ab}=P_aE_b=P_aE_bP_a\ge0$ とおける。$a$ について足すと $E_b$、$b$ について足すと $P_a$ になる。
:::

:::example[一方を弱めるだけでは両立しない場合]
$Z$ の射影測定と $E_+=(I+\eta X)/2$ を含む二値測定を考える。$[(I+Z)/2,E_+]=i\eta Y/2$ なので、共同測定があるのは $\eta=0$ の場合である。両方の測定を弱める場合は射影が消えるため、この可換性の条件はそのまま適用できない。
:::

:::proposition[直交する非鋭い二値測定の両立条件]{#prop-qi-joint-unsharp}
$E_s=(I+s\eta X)/2$ と $F_t=(I+t\mu Z)/2$、$s,t\in\{+1,-1\}$ は、$\eta^2+\mu^2\le1$ のとき、かつそのときに限り一つの POVM の周辺として実現できる。
:::

:::proof
十分性については $G_{st}=(I+s\eta X+t\mu Z)/4$ とおく。固有値は $(1\pm\sqrt{\eta^2+\mu^2})/4$ なので正であり、$s$ または $t$ を足せば対応する周辺になる。必要性については、任意の共同 POVM を $G_{st}=(g_{st}I+\boldsymbol v_{st}\cdot\boldsymbol\sigma)/4$ と書く。正性は $g_{st}\ge\norm{\boldsymbol v_{st}}$ を含意する。周辺条件を解くと $g_{st}=1+st c$、$\boldsymbol v_{st}=s\eta\boldsymbol e_x+t\mu\boldsymbol e_z+st\boldsymbol w$ と書ける。$++$ と $--$ のベクトルの差は $2(\eta\boldsymbol e_x+\mu\boldsymbol e_z)$ なので、三角不等式と正性から $\sqrt{\eta^2+\mu^2}\le1+c$ である。$+-$ と $-+$ についても同じ議論により $\sqrt{\eta^2+\mu^2}\le1-c$ である。両式を足せば必要条件を得る。
:::

:::example[両立の境界で結果を照合する]
$\eta=3/5$、$\mu=4/5$ では上記の $G_{st}$ の各効果は階数 1 である。入力を $\ket0$ とすれば $p(s,t)=(1+4t/5)/4$ なので、$t=+1$ の結果二つが各 $9/20$、$t=-1$ の二つが各 $1/20$ で出る。$s$ を捨てれば $Z$ の粗い測定、$t$ を捨てれば $X$ の粗い測定に一致する。他方、$\eta=\mu=4/5$ では二乗和が $32/25$ に達するため、二つの指定を同じ POVM で再現することはできない。
:::

### 射影測定で実現できる剛性

同じ空間上で二つの POVM を選び、古典乱数で混ぜても新たな POVM になる。座標ごとの凸集合の端点に関心があるとき、射影測定 $(P_j)$ は端点である。実際、$P_j=tE_j+(1-t)F_j$、$0<t<1$ なら、核の二次形式から $E_j,F_j$ は $P_j$ の支持に限られる。各 $E_j$ の和が $I$ だから両側を $P_j$ で挟んで $E_j=P_j$、同様に $F_j=P_j$ である。一方、非射影の POVM も端点になりうる。階数 1 の場合は、各正効果の支持上で効果の係数しか変えられないので、効果の線形独立性が端点であることを保証する。逆に実線形従属なら、従属関係を各係数の小さな増減に使って二つの異なる POVM に分解できる。

:::example[三つの階数 1 効果が線形独立になる]
赤道面の三方向 $(1,0,0),(-1/2,\sqrt3/2,0),(-1/2,-\sqrt3/2,0)$ で $E_j=(I+\boldsymbol n_j\cdot\boldsymbol\sigma)/3$ を取る。$\sum_j d_j E_j=0$ なら $Y$ の係数から $d_2=d_3$、$X$ の係数から $d_1=(d_2+d_3)/2$、定数の係数から和が 0 である。従って全係数は 0 となる。この POVM は量子ビット上の射影測定を単に古典的に選ぶことで再現できない。実際、各選択後の三ラベルの効果が射影または零なら、端点性によりそれぞれが元の POVM に等しい必要があるが、元の固有値は $2/3,0$ である。追加空間を使った射影測定なら、Naimark の構成が適用できる。
:::

## 演習問題

::::exercise{level=basic}
$\ket\psi=(\ket0+(1+i)\ket1)/\sqrt3$ の Bloch ベクトルを求めよ。

:::solution
$2\bar\alpha\beta=2(1+i)/3$ なので、$r_x=r_y=2/3$、$r_z=1/3-2/3=-1/3$ である。二乗和は $1$ となる。
:::
::::

::::exercise{level=basic}
$\ket0\otimes(\ket0+i\ket1)/\sqrt2$ を辞書順の座標で書き、第 2 系の $Y$ の期待値を求めよ。

:::solution
座標は $(1,i,0,0)^\top/\sqrt2$ である。第 2 因子は $Y$ の固有値 $+1$ に属するので $\mel{\psi}{I\otimes Y}{\psi}=1$ である。
:::
::::

::::exercise{level=basic}
$\frac14\ketbra+++\frac34\ketbra--$ の固有値と純度を求めよ。

:::solution
$\ket\pm$ は正規直交基底なので、固有値は $1/4,3/4$ である。純度は $(1/4)^2+(3/4)^2=5/8$ である。
:::
::::

::::exercise{level=basic}
$\ket\Psi=(\ket{00}+\ket{01}+\ket{11})/\sqrt3$ の第 1 系の縮約を求めよ。

:::solution
第 2 系が $0$ の成分は $\ket0/\sqrt3$、$1$ の成分は $(\ket0+\ket1)/\sqrt3$ である。これらの外積を足すと $\frac13\begin{pmatrix}2&1\\1&1\end{pmatrix}$ となる。
:::
::::

::::exercise{level=advanced}
$\rho=(I+X/2)/2$ を最小次元の空間で純粋化せよ。

:::solution
$\ket+$ と $\ket-$ が固有ベクトルで、固有値はそれぞれ $3/4,1/4$ である。したがって $\frac{\sqrt3}2\ket+\ket0+\frac12\ket-\ket1$ を取ればよい。正固有値が 2 個あるため、:ref[thm-qi-purification] より参照系の最小次元は 2 である。
:::
::::

::::exercise{level=advanced}
純粋状態が密度行列全体の凸集合の端点であり、混合状態が端点でないことを証明せよ。

:::solution
$\ketbra\psi\psi=p\rho+(1-p)\sigma$、$0<p<1$ とする。$\psi$ の直交補の任意の $v$ について左辺の二次形式は零である。右辺の項は非負なので $\norm{\sqrt\rho v}^2=\norm{\sqrt\sigma v}^2=0$ となる。よって両者は直交補を消し、トレース 1 であるから $\ketbra\psi\psi$ に等しい。混合状態に正固有値が二つあれば、その二つにそれぞれ $\pm\eps$ を加え、符号を交換した二つの密度行列の平均にできる。$0<\eps$ を両固有値より小さく取れば両行列は相異なる。
:::
::::

::::exercise{level=advanced}
$E_\pm=(I\pm Z/2)/2$ の測定を $\ket+$ に施す。$M_\pm=\sqrt{E_\pm}$ を選んだとき、結果 $+$ に条件付けた状態を求めよ。

:::solution
$E_+=\diag(3/4,1/4)$ なので、$M_+\ket+=(\sqrt3\ket0+\ket1)/(2\sqrt2)$ である。ノルムの二乗は $1/2$ であり、規格化後のベクトルは $(\sqrt3\ket0+\ket1)/2$ となる。
:::
::::

::::exercise{level=advanced}
$B$ の演算を $M_a$、$\sum_aM_a^\dagger M_a=I_B$ とする。結果を捨てると $A$ の縮約が変わらないことを示せ。

:::solution
任意の $L_A$ に対して結果を足した測定後の期待値は $\sum_a\Tr[(L_A\otimes M_a^\dagger M_a)\rho_{AB}]$ である。完備性により $\Tr[(L_A\otimes I_B)\rho_{AB}]$ となるので、:ref[prop-qi-partial-trace] の一意性により縮約が等しい。
:::
::::

## 確認テスト

::::quiz
$\ket+$ と $\ket-$ について正しいものを選べ。

- [x] $X$ 測定で確実に区別できる。
- [ ] $Z$ 測定で確実に区別できる。
- [ ] 二つは大域位相だけが異なる。

:::explanation
両者は $X$ の固有値 $+1$ と $-1$ に属する。$Z$ 測定では両者の結果が等確率で出る。また内積が零なので大域位相による同一視はできない。
:::
::::

::::quiz{answer="0.625" tol="1e-6"}
$\boldsymbol r=(0,0,1/2)$ の状態の純度を求めよ。

:::explanation
固有値は $3/4$ と $1/4$ であり、二乗して足すと $5/8=0.625$ となる。
:::
::::

::::quiz{answer="2" tol="1e-6"}
$\ket{\Phi^+}$ を縮約して得られる 1 量子ビット状態の階数を求めよ。

:::explanation
異なる第 2 系ラベルの交差項は部分トレースで消え、第 1 系には正固有値が $1/2$ ずつ残る。よって階数は 2 である。
:::
::::

::::quiz
POVM が指定する情報はどれか。

- [x] 任意の入力に対する結果の確率。
- [ ] 任意の入力に対する条件付き状態。
- [ ] 補助系の一意な次元。

:::explanation
$E_a$ から $\Tr(E_a\rho)$ が決まる。測定演算子を左からユニタリで変えても効果は変わらず、条件付き状態は変わりうる。また補助系は拡張して使えるので次元は一意でない。
:::
::::

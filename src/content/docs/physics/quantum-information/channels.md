---
title: 量子チャネル
description: 開放系の状態変換を完全正値性で特徴付け、Kraus、Stinespring、Choi の表示を相互に変換する。
notation:
  - symbol: '\Phi'
    meaning: 量子チャネル
  - symbol: 'K_a'
    meaning: Kraus 演算子
  - symbol: 'J(\Phi)'
    meaning: 非規格化した Choi 行列
  - symbol: '\Phi^c'
    meaning: 相補チャネル
---

## この章の目標

- 正性と完全正値性を区別し、転置の反例を計算する。
- 量子チャネルを Kraus、等長、Choi の各表示で表す。
- 代表的な雑音の作用を Bloch 座標で計算する。
- ユニタル写像と固有値の多数化を結び付ける。

## 完全正値写像と CPTP 写像

以下では有限次元の量子系を扱い、Dirac 内積は ket に線形とする。系の一部に作用する変換は、他の系と相関した入力に対しても正の状態を返さなければならない。[開放系](/physics/advanced-quantum-mechanics/open-systems/) の力学を抽象化すると、変換に求める条件は次のように書ける。一般の作用素環での定式化は [完全正値写像](/math/operator-algebras/cp-maps/) で扱う。

:::definition[完全正値性]{#def-qi-cptp}
線形写像 $\Phi:\operatorname{End}(A)\to\operatorname{End}(B)$ が任意の $n\ge1$ で $\id_n\otimes\Phi$ の正性を満たすとき、:term[完全正値]{en="completely positive" yomi="かんぜんせいち"}という。さらにすべての入力 $M$ で $\Tr\Phi(M)=\Tr M$ ならトレース保存である。両条件を満たす写像を :term[量子チャネル]{en="quantum channel" yomi="りょうしちゃねる"}という。
:::

確率混合を準備ラベルによらず変換するため、状態上では凸結合を保つことを要求する。正演算子は非負数と状態の積であり、自己共役演算子は正演算子の差、任意の演算子は二つの自己共役演算子の複素結合だから、状態上の変換が線形拡張を定める。

:::example[正だが完全正でない演算]
転置は固有値を保つので正である。ところが Bell 射影の部分転置は交換演算子の半分であり、反対称部分の固有値が $-1/2$ である。従って補助系を合わせた入力に物理的な状態を返せない。
:::

### Schrödinger 描像と Heisenberg 描像

入力の状態を出力へ写す $\Phi$ に対し、測定効果を入力へ引き戻す $\Phi^\dagger$ を $\Tr[\Phi(T)E]=\Tr[T\Phi^\dagger(E)]$ で定める。非退化なトレース対があるので一意である。$\Phi$ が正なら $\mel v{\Phi^\dagger(E)}v=\Tr(\Phi(\ketbra vv)E)\ge0$ が任意の $E\ge0$ で成り立ち、随伴も正である。この議論を $\id_R\otimes\Phi$ に適用すると完全正値性が随伴にも移る。また、$\Tr\Phi(T)=\Tr T$ が任意の $T$ で成り立つことと、$\Phi^\dagger(I_{\mathrm{out}})=I_{\mathrm{in}}$ は同値である。出力次元と入力次元が異なる場合にも、恒等演算子の下付き文字を残せば区別できる。

:::example[測定統計だけを運ぶ写像]
効果 $E_0=(I+Z/2)/2$、$E_1=I-E_0$ を古典ラベルへ書く写像を $\Phi(T)=\sum_{j=0}^1\Tr(E_jT)\ketbra jj$ とする。$\sqrt{E_j}$ の各行を $\ket j$ に送る Kraus 演算子で実装できる。出力の $X$ の随伴像は零、$Z$ の像は $E_0-E_1=Z/2$ である。例えば入力 $\ket0$ から得る古典分布は $(3/4,1/4)$ なので、$Z$ の出力期待値は $1/2$ に縮む。一方、出力の非対角測定は常に平均零である。
:::

### 条件付き操作と規格化

正の入力のトレースを増やさない完全正写像を、選別された一つの操作として使う。$\Phi^\dagger(I)\le I$ であることと、このトレース条件は同値である。実際、任意の正入力 $T$ に対し $\Tr[T(I-\Phi^\dagger(I))]\ge0$ を調べればよく、必要性には階数 1 の $T$ を使える。複数の操作 $\mathcal I_a$ の和がチャネルなら、結果 $a$ を古典レジスタへ記録できる。付随する POVM の効果は $E_a=\mathcal I_a^\dagger(I_B)$ であり、和は $I_A$ になる。例えば同じ効果を持つ測定でも、$K_a=\sqrt{E_a}$ と $K'_a=U_a\sqrt{E_a}$ では後状態が異なることがある。この操作の族を :term[量子インストルメント]{en="quantum instrument" yomi="りょうしいんすとるめんと"}という。結果 $a$ の確率は $\Tr\mathcal I_a(\rho)$、正の場合の条件付き状態は $\mathcal I_a(\rho)/\Tr\mathcal I_a(\rho)$ である。非規格化演算は線形だが、条件付き状態は規格化の分母が入力に依存する。

### Bloch ベクトルと参照付き入力で正性を判定する

:term[Choi 行列]{en="Choi matrix" yomi="ちょいぎょうれつ"}は、入力の各行列単位の出力を一つのブロック行列に収める。入力と同次元の $R$ を選び、$\ket\Omega=\sum_i\ket i_R\ket i_A$ とする。$J(\Phi)=(\id_R\otimes\Phi)(\ketbra\Omega\Omega)$ は非規格化である。 完全正値性と $J(\Phi)\ge0$ の同値性は後の Choi–Jamiołkowski 同型の節で証明する。

転置と一様状態への置換を混ぜて $\Phi_c(T)=cT^\top+(1-c)\Tr(T)I/2$ とおく。$c$ は実数とする。Bloch ベクトルを $(cr_x,-cr_y,cr_z)$ に変えるので、正値性の必要十分条件は $|c|\le1$ である。$c<0$ では式は二写像の確率混合を意味しない。Choi 行列は $c F+(1-c)I_4/2$ となる。$F$ は交換演算子であり、対称部分で $+1$、反対称部分で $-1$ なので完全正値性には $-1\le c\le1/3$ が必要十分である。各方向の縮小だけを確認しても、参照系付きの正性まで保証できない。

:::example[転置に雑音を足しても負方向が残る]

係数を $c=\frac{1}{2}$ とし、単一の $\ket0$ と参照付きの最大もつれ入力を用意する。後者の出力は $J/2$ になる。

$$
\begin{aligned}
\Phi_c(\ketbra00)&=\left(\begin{matrix}\frac{3}{4} & 0\\0 & \frac{1}{4}\end{matrix}\right)\\
J&=c F+(1-c)I_4/2=\left(\begin{matrix}\frac{3}{4} & 0 & 0 & 0\\0 & \frac{1}{4} & \frac{1}{2} & 0\\0 & \frac{1}{2} & \frac{1}{4} & 0\\0 & 0 & 0 & \frac{3}{4}\end{matrix}\right)\\
\Tr_B J&=\left(\begin{matrix}1 & 0\\0 & 1\end{matrix}\right)\\
\lambda_{\mathrm{sym}}&=(1+c)/2=\frac{3}{4}\\
\lambda_{\mathrm{anti}}&=(1-3c)/2=- \frac{1}{4}
\end{aligned}
$$

反対称部分の負固有値があるので、参照を付けた実験を許す量子チャネルにはできない。$Z$ の固有状態を入れて正の出力が得られても判定には足りない。
:::

:::example[完全正値性の境界]

係数を $c=\frac{1}{3}$ とし、単一の $\ket0$ と参照付きの最大もつれ入力を用意する。後者の出力は $J/2$ になる。

$$
\begin{aligned}
\Phi_c(\ketbra00)&=\left(\begin{matrix}\frac{2}{3} & 0\\0 & \frac{1}{3}\end{matrix}\right)\\
J&=c F+(1-c)I_4/2=\left(\begin{matrix}\frac{2}{3} & 0 & 0 & 0\\0 & \frac{1}{3} & \frac{1}{3} & 0\\0 & \frac{1}{3} & \frac{1}{3} & 0\\0 & 0 & 0 & \frac{2}{3}\end{matrix}\right)\\
\Tr_B J&=\left(\begin{matrix}1 & 0\\0 & 1\end{matrix}\right)\\
\lambda_{\mathrm{sym}}&=(1+c)/2=\frac{2}{3}\\
\lambda_{\mathrm{anti}}&=(1-3c)/2=0
\end{aligned}
$$

反対称方向の重みがちょうど零である。Choi 階数は 3 なので環境を二次元に制限した実装はできない。
:::

:::example[転置の係数に負号を付ける]

係数を $c=- \frac{1}{2}$ とし、単一の $\ket0$ と参照付きの最大もつれ入力を用意する。後者の出力は $J/2$ になる。

$$
\begin{aligned}
\Phi_c(\ketbra00)&=\left(\begin{matrix}\frac{1}{4} & 0\\0 & \frac{3}{4}\end{matrix}\right)\\
J&=c F+(1-c)I_4/2=\left(\begin{matrix}\frac{1}{4} & 0 & 0 & 0\\0 & \frac{3}{4} & - \frac{1}{2} & 0\\0 & - \frac{1}{2} & \frac{3}{4} & 0\\0 & 0 & 0 & \frac{1}{4}\end{matrix}\right)\\
\Tr_B J&=\left(\begin{matrix}1 & 0\\0 & 1\end{matrix}\right)\\
\lambda_{\mathrm{sym}}&=(1+c)/2=\frac{1}{4}\\
\lambda_{\mathrm{anti}}&=(1-3c)/2=\frac{5}{4}
\end{aligned}
$$

この場合は対称部分の重みが小さくなるが、全て非負である。Choi 行列が正なので、これらの重みから Kraus 演算子を構成できる。
:::

### k 正値性と補助系の次元

$k$ 次元の恒等写像をテンソル積しても正であることを $k$ 正値性という。$k$ 正なら $j\le k$ でも正である。$j$ 次元系を $k$ 次元へ等長に埋め込み、写像の出力を同じ支持へ圧縮すればよい。入力の次元を $d$ とすると、$d$ 正値性は完全正値性と同値である。実際、任意次元の参照系付きの純粋入力には Schmidt 階数 $d$ 以下の展開がある。参照系をその Schmidt 支持から $d$ 次元へ等長に埋め込み、その中での正性を使えば、元の入力の出力は正である。参照次元が $d$ より大きい場合には、必要な参照基底だけを選んで埋め込む。混合入力は純粋射影の非負係数付き和なので、線形性で全部の正入力を扱える。

:::example[三準位で二正だが完全正でない写像]
$d\ge2$ で $\mathcal R_c(T)=c\Tr(T)I-T$、$c\ge0$ とする。規格化した純粋な参照付き入力 $\psi$ の出力は $c\,\rho_R\otimes I-\ketbra\psi\psi$ である。$\rho_R$ の支持で逆平方根を両側から掛けると $cI-\ketbra vv$ となる。ここで $\ket v=\sum_{j=1}^r \ket{u_j}\ket{w_j}$ は元の Schmidt 展開から全重みの平方根を消したベクトルであり、$\norm v^2=r$ である。従って正性の条件は $c\ge r$ である。$d=3,c=2$ なら、参照が二次元以下で常に正だが三次元最大もつれに対して負になる。跡は $(cd-1)\Tr T=5\Tr T$ だから、写像全体を 5 で割ればトレース保存のまま反例となる。正の出力が得られる補助次元の範囲を指定しないと、完全正値性の試験にならない。
:::

### ユニタルな完全正値写像の Schwarz 不等式

$\Phi(I)=I$ とする。任意の行列 $A$ に対し
$$
\begin{pmatrix} A^\dagger A&A^\dagger\\ A&I\end{pmatrix}
=\begin{pmatrix}A^\dagger\\ I\end{pmatrix}\begin{pmatrix}A&I\end{pmatrix}
\ge0.
$$
完全正値性から各ブロックに $\Phi$ を掛けても正である。正写像は自己共役行列を自己共役行列に移すので、任意の $A$ を Hermite 部と反 Hermite 部に分ければ $\Phi(A^\dagger)=\Phi(A)^\dagger$ である。従って変換後の二次形式に列 $(v,-\Phi(A)v)$ を入れると、
$$
\Phi(A^\dagger A)\ge\Phi(A)^\dagger\Phi(A)
$$
となる。この証明では $2$ 正値性で足りる。例えば $A$ がユニタリなら $\Phi(A)^\dagger\Phi(A)\le I$ である。ユニタル正写像一般についてこの形の Schwarz 不等式を使うには仮定を確認する。転置 $T$ と $A=\ketbra01$ に対して、左辺は $\ketbra11$、右辺は $\ketbra00$ なので反例となる。

:::example[位相緩和で不等式の差を測る]
$\Phi(T)=(1-p)T+p ZTZ$、$p=1/4$ とする。$A=\ketbra01$ では $\Phi(A)=A/2$、$\Phi(A^\dagger A)=\ketbra11$ なので、差は $(3/4)\ketbra11$ である。一方、$A=Z$ なら $\Phi(Z)=Z$ が保たれ、差は零である。この差は、入力の積を先に計算してから作用させる手順と、個々に作用させて積を取る手順で何が失われるかを表す。全成分を保つユニタリ共役なら、任意の $A$ で差は零になる。
:::

### Schwarz の等号から積の保存を導く

ユニタル完全正写像を $\Phi(M)=V^\dagger(M\otimes I)V$ とする。これは $\Phi^\dagger$ がトレース保存であることと、後で証明する等長表示を随伴へ移した式である。$P=VV^\dagger$ と書くと、Schwarz の差は
$$
\Phi(A^\dagger A)-\Phi(A)^\dagger\Phi(A)
=V^\dagger(A^\dagger\otimes I)(I-P)(A\otimes I)V
$$
である。従って等号なら $(I-P)(A\otimes I)V=0$ である。$A^\dagger$ でも等号なら、同様に $V^\dagger(A\otimes I)(I-P)=0$ を得る。この二つの零条件を使い積の途中に $P$ を挿入すると、任意の $M$ に対し $\Phi(AM)=\Phi(A)\Phi(M)$ と $\Phi(MA)=\Phi(M)\Phi(A)$ となる。このような $A$ の領域を :term[乗法領域]{en="multiplicative domain" yomi="じょうほうりょういき"}という。

:::example[計算基底の測定で保存する積]
非対角成分を消す写像では、対角な全ての $A$ が二つの等号を満たす。逆に $A=\begin{pmatrix}a&b\\c&d\end{pmatrix}$ に対する Schwarz の差は $\diag(|c|^2,|b|^2)$ である。これが零になるには $b=c=0$ が必要なので、このチャネルの乗法領域は対角行列全体である。
:::

## Kraus 表示

:::proposition[Kraus 演算子による構成]{#prop-qi-kraus}
$\Phi(M)=\sum_a K_a M K_a^\dagger$ の形の写像は完全正である。トレース保存は $\sum_aK_a^\dagger K_a=I_A$ と同値である。
:::

:::proof
拡張写像の式では各 $K_a$ を $I_n\otimes K_a$ に置き換える。正演算子の合同変換は正なので和も正である。巡回性により、すべての $M$ で $\Tr(M(\sum_aK_a^\dagger K_a-I_A))=0$ となることがトレース保存条件である。差の随伴を $M$ に選べば必要性も従う。
:::

この形を :term[Kraus 表示]{en="Kraus representation" yomi="くらうすひょうじ"}と呼び、各 $K_a$ は入力から出力への :term[Kraus 演算子]{en="Kraus operator" yomi="くらうすえんざんし"}である。例えば $K_0=\sqrt{1-p}\,I$、$K_1=\sqrt p\,X$ では確率 $p$ でビット反転する。$K'_b=\sum_a U_{ba}K_a$ とユニタリで混ぜても、和では $\sum_bU_{ba}\overline{U_{bc}}=\delta_{ac}$ によって写像は変わらない。個別の枝を測定する装置の選択にはこの自由度が現れる。

### Kraus 演算子の自由度
$\operatorname{vec}K=\sum_i\ket i\otimes K\ket i$ とし、Kraus 演算子をベクトル化した列を行列 $L$ にまとめる。そのとき非規格化 Choi 行列は $J=LL^\dagger$ である。$J$ の支持で $J^{-1/2}$ を定めれば、$L=\sqrt J\,W$ で $WW^\dagger$ は支持射影になる。これにより、長さが同じ二つの Kraus リストでは $W$ の行の正規直交系をそれぞれ基底に延長し、一方のリストを他方へユニタリ混合できる。長さが違う場合は零演算子を付け足して同じ長さにできる。支持の次元 $r=\rank J$ は必要な列数の下限であり、$J$ のスペクトル分解で $r$ 個の列を取れば達成する。

:::example[位相反転とコヒーレンスの消去]
$K_0=I/\sqrt2,K_1=Z/\sqrt2$ のリストに二次元の Hadamard の係数を掛けると、$L_0=(I+Z)/2$、$L_1=(I-Z)/2$ になる。すなわち公平な位相反転と、計算基底測定の結果を捨てる操作は同じチャネルである。一方、補助系の測定ラベルは前者で入力に依存せず、後者で $\rho$ の対角成分を返す。Kraus 表示の自由度はデータだけを保持した場合の同値性である。
:::

### チャネルの端点を判定する
最小 Kraus リストを $(K_j)_{j=1}^r$ とする。$\{K_i^\dagger K_j\}_{i,j}$ が複素線形独立であれば CPTP の凸集合の端点であり、従属なら端点でない。この主張を示すため、$J$ の任意の凸分解 $tJ_1+(1-t)J_2$ を取る。正性により $J_\alpha$ は $J$ の支持に含まれる。Kraus ベクトルを列に持つ列階数 $r$ の行列 $L$ で $J_\alpha=L M_\alpha L^\dagger$ と一意に書ける。$M_\alpha$ は正で、$tM_1+(1-t)M_2=I$ である。トレース保存は $\sum_{ij}(M_\alpha)_{ij}K_j^\dagger K_i=I$ と同値なので、独立なら $M_\alpha=I$ に限る。従属なら $\sum_{ij}D_{ij}K_j^\dagger K_i=0$ となる非零 $D$ を取る。この線形写像は随伴と交換するので、$D$ の Hermite 部と反 Hermite 部から非零の Hermite 核元を選べる。$I\pm\eps D\ge0$ とでき、対応する二つの $J$ は異なるチャネルを表す。従ってその中点は端点でない。

### 行列要素で Kraus 演算子の正規化を検査する

以下では $K_i^\dagger K_i$ の和を計算してから、特定の入力 $\rho$ に対し各 $\tau_i=K_i\rho K_i^\dagger$ を求める。入力は各測定の結果で条件付ける前の状態である。従って $p_i=\Tr\tau_i$ の総和は 1 になり、無条件の結果は各非規格化行列の和になる。

:::example[位相ではなくビットを反転する環境]

入力を次の密度行列とし、演算子番号を環境で読める実装を仮定する。

$$
\begin{aligned}
\rho&=\left(\begin{matrix}\frac{9}{10} & - \frac{3 i}{10}\\\frac{3 i}{10} & \frac{1}{10}\end{matrix}\right)\\
K_0&=\left(\begin{matrix}\frac{\sqrt{3}}{2} & 0\\0 & \frac{\sqrt{3}}{2}\end{matrix}\right)\\
K_1&=\left(\begin{matrix}0 & \frac{1}{2}\\\frac{1}{2} & 0\end{matrix}\right)\\
K_0^\dagger K_0&=\left(\begin{matrix}\frac{3}{4} & 0\\0 & \frac{3}{4}\end{matrix}\right)\\
K_1^\dagger K_1&=\left(\begin{matrix}\frac{1}{4} & 0\\0 & \frac{1}{4}\end{matrix}\right)\\
\tau_0&=K_0\rho K_0^\dagger=\left(\begin{matrix}\frac{27}{40} & - \frac{9 i}{40}\\\frac{9 i}{40} & \frac{3}{40}\end{matrix}\right)\\
\tau_1&=K_1\rho K_1^\dagger=\left(\begin{matrix}\frac{1}{40} & \frac{3 i}{40}\\- \frac{3 i}{40} & \frac{9}{40}\end{matrix}\right)\\
(p_0,p_1)&=(\frac{3}{4},\frac{1}{4})\\
\Phi(\rho)&=\tau_0+\tau_1=\left(\begin{matrix}\frac{7}{10} & - \frac{3 i}{20}\\\frac{3 i}{20} & \frac{3}{10}\end{matrix}\right)
\end{aligned}
$$

反転枝だけを選ぶと $Y,Z$ 成分の符号は共に変わるが、無条件では二つの期待値が半分に縮む。跡を保存することは、各枝でノルムを保存することとは異なる。
:::

:::example[確率的な選択と入力に依存する人口測定]

入力を次の密度行列とし、演算子番号を環境で読める実装を仮定する。

$$
\begin{aligned}
\rho&=\left(\begin{matrix}\frac{9}{10} & \frac{3}{10}\\\frac{3}{10} & \frac{1}{10}\end{matrix}\right)\\
K_0&=\left(\begin{matrix}\frac{\sqrt{6}}{3} & 0\\0 & \frac{\sqrt{6}}{6}\end{matrix}\right)\\
K_1&=\left(\begin{matrix}\frac{\sqrt{3}}{3} & 0\\0 & \frac{\sqrt{30}}{6}\end{matrix}\right)\\
K_0^\dagger K_0&=\left(\begin{matrix}\frac{2}{3} & 0\\0 & \frac{1}{6}\end{matrix}\right)\\
K_1^\dagger K_1&=\left(\begin{matrix}\frac{1}{3} & 0\\0 & \frac{5}{6}\end{matrix}\right)\\
\tau_0&=K_0\rho K_0^\dagger=\left(\begin{matrix}\frac{3}{5} & \frac{1}{10}\\\frac{1}{10} & \frac{1}{60}\end{matrix}\right)\\
\tau_1&=K_1\rho K_1^\dagger=\left(\begin{matrix}\frac{3}{10} & \frac{\sqrt{10}}{20}\\\frac{\sqrt{10}}{20} & \frac{1}{12}\end{matrix}\right)\\
(p_0,p_1)&=(\frac{37}{60},\frac{23}{60})\\
\Phi(\rho)&=\tau_0+\tau_1=\left(\begin{matrix}\frac{9}{10} & \frac{1}{10} + \frac{\sqrt{10}}{20}\\\frac{1}{10} + \frac{\sqrt{10}}{20} & \frac{1}{10}\end{matrix}\right)
\end{aligned}
$$

各効果の対角の違いにより、結果の確率は入力人口に依存する。一方、全分枝を足すと対角人口は保存され、非対角成分だけが減衰する。
:::

:::example[固定した準備へ戻す二つの演算子]

入力を次の密度行列とし、演算子番号を環境で読める実装を仮定する。

$$
\begin{aligned}
\rho&=\left(\begin{matrix}\frac{4}{5} & - \frac{2 i}{5}\\\frac{2 i}{5} & \frac{1}{5}\end{matrix}\right)\\
K_0&=\left(\begin{matrix}1 & 0\\0 & 0\end{matrix}\right)\\
K_1&=\left(\begin{matrix}0 & 1\\0 & 0\end{matrix}\right)\\
K_0^\dagger K_0&=\left(\begin{matrix}1 & 0\\0 & 0\end{matrix}\right)\\
K_1^\dagger K_1&=\left(\begin{matrix}0 & 0\\0 & 1\end{matrix}\right)\\
\tau_0&=K_0\rho K_0^\dagger=\left(\begin{matrix}\frac{4}{5} & 0\\0 & 0\end{matrix}\right)\\
\tau_1&=K_1\rho K_1^\dagger=\left(\begin{matrix}\frac{1}{5} & 0\\0 & 0\end{matrix}\right)\\
(p_0,p_1)&=(\frac{4}{5},\frac{1}{5})\\
\Phi(\rho)&=\tau_0+\tau_1=\left(\begin{matrix}1 & 0\\0 & 0\end{matrix}\right)
\end{aligned}
$$

どの枝でも後状態が同じであり、出力の純度が上がる。Kraus が複数あることだけから、任意の入力が混合出力になると判断することはできない。
:::

### 跡を減らす枝を一つの測定へ組み込む

演算子 $L_j$ の和 $E=\sum_jL_j^\dagger L_j$ が $I$ 以下なら、成功枝は $T\mapsto\sum_j L_jT L_j^\dagger$ と書ける。失敗枝に $L_f=\sqrt{I-E}$ を加えると、同一の出力空間の場合には完備になる。出力空間が違う場合は失敗用の直交フラグを別に設け、$L_f$ の各行をフラグへ移す演算子で測定の失敗を実装できる。成功枝の非規格化状態をそのトレースで割ると条件付き状態になり、成功と失敗の両方を合わせるまではトレース保存ではない。

:::example[振幅減衰の雑音を見つける場合と見失う場合]
$K_0=\diag(1,1/2)$、$K_1=(\sqrt3/2)\ketbra01$ とし、$\ket+$ を入れる。ジャンプなしの非規格化 ket は $(\ket0+\ket1/2)/\sqrt2$、確率は $5/8$ である。ジャンプ枝は $\sqrt{3/8}\ket0$ で確率 $3/8$ である。ジャンプを検出して捨てれば、残りの状態は $(2\ket0+\ket1)/\sqrt5$ となり $Z$ の期待値は $3/5$ である。環境のラベルを読まず全枝を足せば、平均は $(5/8)(3/5)+(3/8)=3/4$ になる。検出器で記録された条件付き状態と、無条件のチャネルの出力では一般に同じ観測値を持たない。
:::

:::example[端点としての振幅減衰]
$0<\gamma<1$ の振幅減衰の最小 Kraus リスト $K_0,K_1$ を使う。$K_0^\dagger K_1=\sqrt\gamma\ketbra01$、$K_1^\dagger K_0=\sqrt\gamma\ketbra10$ である。対角の二積は $\diag(1,1-\gamma)$ と $\diag(0,\gamma)$ で、第一成分と第二成分から独立性が分かる。従って四積は線形独立であり、このチャネルは CPTP の凸集合の端点である。一方、$I$ と $Z$ を確率混合する位相反転チャネルでは $K_0^\dagger K_0$ と $K_1^\dagger K_1$ がともに $I$ に比例する。出力に同じ数の Kraus 演算子が必要であっても、凸混合へ分解できるかどうかは異なる。
:::

## Stinespring 表現

:::proposition[有限次元の等長表示]{#prop-qi-stinespring}
有限次元の量子チャネルは $\Phi(M)=\Tr_E(VMV^\dagger)$ と書ける。ここで $V:A\to B\otimes E$ は等長写像である。逆にこの式で定めた写像は量子チャネルである。
:::

:::proof
完全正値写像が Kraus 表示を持つことは次節の Choi 判定で示す。Kraus 表示を使い $V=\sum_aK_a\otimes\ket a_E$ と置くと、$V^\dagger V=\sum_aK_a^\dagger K_a=I$ となる。部分トレースは直交ラベル間の交差項を消すので所望の写像を得る。逆に任意の等長写像を環境基底で $K_a=(I_B\otimes\bra a)V$ と展開すれば、Kraus 表示と完備性が出る。
:::

等長写像と部分トレースで表す方法を :term[Stinespring 表現]{en="Stinespring representation" yomi="すたいんすぷりんぐひょうげん"}と呼ぶ。$E$ を残した $\Phi^c(M)=\Tr_B VMV^\dagger$ を :term[相補チャネル]{en="complementary channel" yomi="そうほちゃねる"}という。例えば $V\ket x=\ket x\ket x$ では両側で同じ対角成分が得られる。重ね合わせの非対角成分は片側の状態から消えるが、全体の相関に残る。

### 相補チャネルと環境へ移る情報

等長実装 $V=\sum_i K_i\otimes\ket i_E$ に対し、相補チャネルの $(i,j)$ 成分は $\Tr(K_i T K_j^\dagger)$ である。最小 Kraus 表示のユニタリ自由度は、環境出力をユニタリ共役する自由度に対応する。非最小実装ではさらに等長埋め込みを許す。このため環境出力の正固有値は、どの最小表現を取っても同じである。

:::example[完全消去による環境への移送]
$V\ket j=\ket0_B\ket j_E$ なら $\Phi(T)=\Tr(T)\ketbra00$、$\Phi^c(T)=T$ である。出力 $B$ の測定が入力の情報を得ない場合でも、環境の量子状態は入力を保存しうる。一方、$V\ket j=\ket j_B\ket0_E$ なら $\Phi$ は恒等で $\Phi^c(T)=\Tr(T)\ketbra00$ になる。二つの実装は $B$ と $E$ の役割を交換したものである。
:::

:::example[位相緩和の環境の二つの状態]
$V\ket0=\ket0(\sqrt{1-p}\ket0+\sqrt p\ket1)$、$V\ket1=\ket1(\sqrt{1-p}\ket0-\sqrt p\ket1)$ とする。環境ベクトルの内積は $1-2p$ であり、受信側の非対角成分の倍率も $1-2p$ である。$p=1/4$ なら環境が等確率の二入力を識別する最適成功率は $(1+\sqrt3/2)/2$ になる。$p=1/2$ で環境の二状態は直交し、受信側のコヒーレンスは零になる。
:::

### 受け取る系と環境に残る系を比較する

二つの Kraus 演算子を環境の二基底に記録すると、純粋な入力 $\psi$ に対して係数行列は $C=(K_0\psi,K_1\psi)$ になる。等長性は $\Tr(CC^\dagger)=1$ で読み取れる。受信出力は $CC^\dagger$、環境出力は $(C^\dagger C)^\top$ であり、同じ正固有値を持つ。以下では同じ物理実装に対して入力を選び、純粋出力と混合出力が入れ替わる様子まで計算する。

:::example[基底入力が環境に何も記録しない場合]

減衰率を $3/4$ に固定し、環境を $\ket0$ に準備する。各欄で $B$ は受信系、$E$ は補助系である。

$$
\begin{aligned}
\psi&=\left(\begin{matrix}1\\0\end{matrix}\right)\\
C&=(K_0\psi,K_1\psi)=\left(\begin{matrix}1 & 0\\0 & 0\end{matrix}\right)\\
\rho_B&=CC^\dagger=\left(\begin{matrix}1 & 0\\0 & 0\end{matrix}\right)\\
\rho_E&=(C^\dagger C)^\top=\left(\begin{matrix}1 & 0\\0 & 0\end{matrix}\right)\\
\det \rho_B=\det \rho_E&=0\\
\lambda(\rho_B)=\lambda(\rho_E)&=(0,1)
\end{aligned}
$$

環境も受信側も純粋である。入力の $\ket0$ は変化せず、ジャンプラベルの確率も零である。
:::

:::example[励起状態を送って環境のジャンプを読む]

減衰率を $3/4$ に固定し、環境を $\ket0$ に準備する。各欄で $B$ は受信系、$E$ は補助系である。

$$
\begin{aligned}
\psi&=\left(\begin{matrix}0\\1\end{matrix}\right)\\
C&=(K_0\psi,K_1\psi)=\left(\begin{matrix}0 & \frac{\sqrt{3}}{2}\\\frac{1}{2} & 0\end{matrix}\right)\\
\rho_B&=CC^\dagger=\left(\begin{matrix}\frac{3}{4} & 0\\0 & \frac{1}{4}\end{matrix}\right)\\
\rho_E&=(C^\dagger C)^\top=\left(\begin{matrix}\frac{1}{4} & 0\\0 & \frac{3}{4}\end{matrix}\right)\\
\det \rho_B=\det \rho_E&=\frac{3}{16}\\
\lambda(\rho_B)=\lambda(\rho_E)&=(\frac{1}{4},\frac{3}{4})
\end{aligned}
$$

二つの系の状態は共に混合だが、環境が励起なら受信側は基底であり、その逆も成り立つ。各側の局所確率を二つ掛けるとこの反相関は出ない。
:::

:::example[重ね合わせが環境の非対角成分にも現れる]

減衰率を $3/4$ に固定し、環境を $\ket0$ に準備する。各欄で $B$ は受信系、$E$ は補助系である。

$$
\begin{aligned}
\psi&=\left(\begin{matrix}\frac{3}{5}\\\frac{4 i}{5}\end{matrix}\right)\\
C&=(K_0\psi,K_1\psi)=\left(\begin{matrix}\frac{3}{5} & \frac{2 \sqrt{3} i}{5}\\\frac{2 i}{5} & 0\end{matrix}\right)\\
\rho_B&=CC^\dagger=\left(\begin{matrix}\frac{21}{25} & - \frac{6 i}{25}\\\frac{6 i}{25} & \frac{4}{25}\end{matrix}\right)\\
\rho_E&=(C^\dagger C)^\top=\left(\begin{matrix}\frac{13}{25} & - \frac{6 \sqrt{3} i}{25}\\\frac{6 \sqrt{3} i}{25} & \frac{12}{25}\end{matrix}\right)\\
\det \rho_B=\det \rho_E&=\frac{48}{625}\\
\lambda(\rho_B)=\lambda(\rho_E)&=(\frac{1}{2} - \frac{\sqrt{433}}{50},\frac{\sqrt{433}}{50} + \frac{1}{2})
\end{aligned}
$$

環境は単なる古典ジャンプ記録ではない。環境を別の基底で測ると、非対角成分にも応じた確率を得る。受信側との同じ固有値は、行列要素が一致することではなく、全出力が純粋なことから従う。
:::

### Stinespring 実装を一つのユニタリに埋め込む

入出力空間の次元が等しい場合、等長写像 $V:A\to A\otimes E$ の各基底像は正規直交系である。入力を $\ket j\ket0_E$ の部分空間に埋め込み、そこを $V\ket j$ へ移すようにユニタリの列を指定する。残りの列を正規直交基底に延長して $U$ を得る。異なる次元の入出力なら、大きい共通空間に両方の列を埋め込んだ上で同じ構成を行う。環境の最小次元はチャネルの Choi 階数だが、ユニタリ全体の残りの列には自由度がある。この自由度は、補助入力を $\ket0$ に固定しない実験では観測に現れうる。

:::example[減衰の角度を交換ゲートで実現する]
$\gamma=\sin^2\theta$ とし、$\ket{00},\ket{11}$ を固定し、$\ket{10},\ket{01}$ の張る平面で
$$
\begin{aligned}
U\ket{10}&=\cos\theta\ket{10}+\sin\theta\ket{01},\\
U\ket{01}&=-\sin\theta\ket{10}+\cos\theta\ket{01}
\end{aligned}
$$
とする。この二列は直交しているので全体はユニタリである。補助ビットを $\ket0$ にすると振幅減衰になり、$\theta=\pi/6$ なら励起人口は $3/4$ 倍、非対角成分は $\sqrt3/2$ 倍である。補助を $\ket1$ で用意した場合は、反対に $\ket0$ から $\ket1$ へ人口を移す。この実装ではチャネルの性質は回路だけでなく、環境の初期条件を含めて決まる。
:::

## Choi–Jamiołkowski 同型

Choi 行列を使うと、完全正値性を一つの有限行列で判定できる。

:::theorem[Choi の判定条件]{#thm-qi-choi}
$\Phi$ の完全正値性、$J(\Phi)\ge0$、Kraus 表示の存在は同値である。トレース保存は $\Tr_BJ(\Phi)=I_R$ と同値である。同じ基底の同一視で $\Tr_R J(\Phi)=\Phi(I_A)$ であるから、出力側の縮約はユニタル性を判定する。完全正値写像の最少 Kraus 演算子数は $\rank J(\Phi)$ である。トレース保存も満たす場合、等長表示の最小環境次元は同じ階数である。
:::

:::proof
完全正値性から $J\ge0$ は定義である。逆に $J\ge0$ ならスペクトル分解で $J=\sum_a\ketbra{w_a}{w_a}$ とし、固有値の平方根を $w_a$ に含める。$w_a=\sum_i\ket i\otimes K_a\ket i$ で $K_a$ を定める。$ij$ ブロックを比較すると $\Phi(\ketbra ij)=\sum_aK_a\ketbra ijK_a^\dagger$ である。線形に拡張して Kraus 表示を得る。:ref[prop-qi-kraus] により完全正値性へ戻る。

$\Tr_BJ$ の $ij$ 成分は $\Tr\Phi(\ketbra ij)$ なので、行列単位でのトレース保存がこの縮約条件に等しい。$k$ 個の Kraus 演算子による表示は $k$ 個の階数 1 の外積として $J$ を表すので $\rank J\le k$ である。上のスペクトル構成が等号を実現する。環境については :ref[prop-qi-stinespring] の両方向の構成から従う。
:::

入力 $M$ への写像の作用は $\Phi(M)=\Tr_R[(M^\top\otimes I_B)J(\Phi)]$ で復元できる。実際、$\Tr(M^\top\ketbra ij)=M_{ij}$ が対応するブロックを選ぶ。異なる最小 Kraus 表示はユニタリで関係する。対応する列行列 $W,W'$ が $WW^\dagger=W'W'^\dagger=J$ を満たすので、$J$ の支持上の逆平方根を左から掛けると同じ支持空間の正規直交基底が得られる。

### Choi 行列から実験可能な量を読む

入力を固定した純粋状態 $\ket\psi$ とすると、$(\bra{\bar\psi}\otimes I)J(\ket{\bar\psi}\otimes I)=\Phi(\ketbra\psi\psi)$ である。Choi 行列を $\sum_{ij}\ketbra ij\otimes\Phi(\ketbra ij)$ と展開し、両側の縮約で得る係数が $\psi_i\bar\psi_j$ であることを確認する。また二写像を合成した行列には $J_{\Psi\circ\Phi}=(\id\otimes\Psi)(J_\Phi)$ を使える。正規化した最大もつれ入力から実験的に状態 $J/d$ を推定した場合、チャネルの行列へ戻すときに入力次元 $d$ を掛ける必要がある。

:::example[振幅減衰の Choi 行列を逆に解く]
入力基底を左に置き、減衰確率を $\gamma=3/4$ とする。$J$ は $\ket{00}+\frac12\ket{11}$ の非規格化外積と、$\frac34\ketbra{10}{10}$ の和である。正固有値は $5/4,3/4$ で、跡は 2 である。各列ベクトルを $\operatorname{vec}K$ から戻すと $K_0=\diag(1,1/2)$、$K_1=(\sqrt3/2)\ketbra01$ になる。$K_0^\dagger K_0+K_1^\dagger K_1=I$ であるため、出力を部分トレースした Choi 行列も $I$ になる。二つの正固有値があるので一演算子の実装は不可能である。
:::

### 測定して作り直すチャネル

チャネルが $\Phi(T)=\sum_a\Tr(E_aT)\sigma_a$ なら、正規化 Choi 状態は $\sum_a E_a^\top\otimes\sigma_a/d$ であり分離可能である。逆に $J/d=\sum_a p_a\tau_a\otimes\sigma_a$ が状態の分離可能分解なら、$\Tr_B J=I$ より効果 $E_a=dp_a\tau_a^\top$ の和は $I$ である。Choi 逆写像を適用すれば測定と準備の形を回復する。この形を持つチャネルを :term[もつれ破壊チャネル]{en="entanglement-breaking channel" yomi="もつれはかいちゃねる"}と呼ぶ。参照系付き入力へ作用させると、結果は $\sum_a\Tr_A[(I_R\otimes E_a)\rho_{RA}]\otimes\sigma_a$ である。第一因子は $\sqrt{E_a}$ を両側から作用させて部分トレースした結果に等しく、非規格化正状態なので出力は分離可能である。逆に任意の参照系とのもつれを消すチャネルは、最大もつれ入力に対しても消すので、Choi 状態が分離可能になる。

### Choi 行列の固有値を変換の組として読む

非規格化 Choi 行列の正固有値を $\lambda_j$、対応する正規直交固有ベクトルを $v_j$ とする。ベクトル $\sqrt{\lambda_j}v_j$ の参照ラベルごとの出力成分を列に戻せば Kraus 演算子が得られる。縮退した固有値では別の正規直交基底を選んでもよい。同じ結果が $K_j$ のユニタリ混合によっても得られるのは、Choi 行列が各列の外積の和だからである。

:::example[三つの Pauli 成分の Choi 像]

入出力の基底を固定し、第一レジスタを入力の参照ラベルとし、計算では $\Tr_A$ を参照側のトレースに用いる。次の Choi 行列から入出力の条件を確認する。

$$
\begin{aligned}
J&=\left(\begin{matrix}\frac{1}{2} & 0 & 0 & \frac{1}{2}\\0 & \frac{1}{2} & 0 & 0\\0 & 0 & \frac{1}{2} & 0\\\frac{1}{2} & 0 & 0 & \frac{1}{2}\end{matrix}\right)\\
\lambda(J)&=(1,\frac{1}{2},\frac{1}{2},0)\\
\Tr_B J&=\left(\begin{matrix}1 & 0\\0 & 1\end{matrix}\right)\\
\rho_{\mathrm{in}}&=\left(\begin{matrix}\frac{4}{5} & 0\\0 & \frac{1}{5}\end{matrix}\right)\\
\Phi(\rho_{\mathrm{in}})&=\Tr_A[(\rho_{\mathrm{in}}^\top\otimes I)J]=\left(\begin{matrix}\frac{1}{2} & 0\\0 & \frac{1}{2}\end{matrix}\right)
\end{aligned}
$$

三つの非零固有値は直交する Pauli ベクトルに対応し、第 2 と第 3 は縮退する。縮退部分の二列に相対位相を付けて混ぜてもチャネルは同じである。ただし、環境をその基底で測る場合は読める記録が変わりうる。
:::

:::example[測定して非直交状態を準備する Choi 像]

入出力の基底を固定し、第一レジスタを入力の参照ラベルとし、計算では $\Tr_A$ を参照側のトレースに用いる。次の Choi 行列から入出力の条件を確認する。

$$
\begin{aligned}
J&=\left(\begin{matrix}1 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & \frac{9}{25} & \frac{12}{25}\\0 & 0 & \frac{12}{25} & \frac{16}{25}\end{matrix}\right)\\
\lambda(J)&=(1,1,0,0)\\
\Tr_B J&=\left(\begin{matrix}1 & 0\\0 & 1\end{matrix}\right)\\
\rho_{\mathrm{in}}&=\left(\begin{matrix}\frac{4}{5} & - \frac{2 i}{5}\\\frac{2 i}{5} & \frac{1}{5}\end{matrix}\right)\\
\Phi(\rho_{\mathrm{in}})&=\Tr_A[(\rho_{\mathrm{in}}^\top\otimes I)J]=\left(\begin{matrix}\frac{109}{125} & \frac{12}{125}\\\frac{12}{125} & \frac{16}{125}\end{matrix}\right)
\end{aligned}
$$

準備する二つの状態は非直交だが、Choi ベクトルの参照ラベルが異なるため二列は直交している。入力の非対角成分は最初の測定で失われ、出力の非対角成分は準備した状態に由来する。
:::

:::example[三準位の入力を二ラベルへ圧縮する Choi 像]

入出力の基底を固定し、第一レジスタを入力の参照ラベルとし、計算では $\Tr_A$ を参照側のトレースに用いる。次の Choi 行列から入出力の条件を確認する。

$$
\begin{aligned}
J&=\left(\begin{matrix}1 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 1 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 1\end{matrix}\right)\\
\lambda(J)&=(1,1,1,0,0,0)\\
\Tr_B J&=\left(\begin{matrix}1 & 0 & 0\\0 & 1 & 0\\0 & 0 & 1\end{matrix}\right)\\
\rho_{\mathrm{in}}&=\left(\begin{matrix}\frac{1}{2} & 0 & 0\\0 & \frac{1}{3} & 0\\0 & 0 & \frac{1}{6}\end{matrix}\right)\\
\Phi(\rho_{\mathrm{in}})&=\Tr_A[(\rho_{\mathrm{in}}^\top\otimes I)J]=\left(\begin{matrix}\frac{1}{2} & 0\\0 & \frac{1}{2}\end{matrix}\right)
\end{aligned}
$$

全体の跡は出力次元 2 ではなく入力次元 3 になる。入力の第 2 と第 3 の人口を足したものが出力の第 2 人口であり、正規化には入力側の次元を使う。
:::

### 過程トモグラフィーと Pauli 伝達行列

Choi 行列は、参照側と出力側の行列の積を測定することで推定できる。例えば $J/d$ の $\Tr[(A^\top\otimes B)(J/d)]$ は $\Tr[B\Phi(A)]/d$ に等しい。行列単位の展開で第一因子の跡は $\Tr(A^\top\ketbra ij)=A_{ij}$ となるためである。$A$ が状態である必要はなく、Hermite な基底行列を使ってその係数をまとめて求められる。二量子ビットの Choi 状態で Pauli 積を測ると、見ているのは状態の再構成と同時に、元の変換の Pauli 伝達係数である。第一因子に $Y$ が現れる場合に転置の符号を保存する必要がある。

:::example[Pauli 伝達行列から位相緩和を復元する]
一量子ビットで $J/2=(II+cXX-cYY+ZZ)/4$ が得られたとする。Choi 逆写像から $\Phi(I)=I,\Phi(X)=cX,\Phi(Y)=cY,\Phi(Z)=Z$ を得る。$J/2$ の固有基底は Bell 基底で、$\Phi^+,\Phi^-$ に固有値 $(1+c)/2,(1-c)/2$、残りは 0 である。したがって $|c|\le1$ のとき物理的であり、位相反転確率は $(1-c)/2$ である。例として $c=-1/2$ なら確率 $3/4$ で位相を変える。伝達係数の負符号は物理性の破れではなく、平均した Bloch 成分の向きの反転を表す。
:::

:::example[Choi 行列のトレースだけでは足りない]
$J=2\ketbra{00}{00}$ を候補とすると、正で跡は 2 である。対応する写像は $\Phi(T)=2\mel0T0\ketbra00$ で、入力 $\ket0$ の跡を 2 にする。実際 $\Tr_B J=2\ketbra00$ で $I$ ではない。二次元入力すべての跡を保存するには、$J$ の全トレースを 2 に合わせるだけでなく参照側の行列を $I$ に合わせる必要がある。例えば $J'= \ketbra{00}{00}+\ketbra{10}{10}$ に直せば $\Tr_B J'=I$ となり、固定純粋状態への置換チャネルになる。
:::

## 代表的なチャネル（脱分極・位相緩和・振幅減衰）

### 脱分極と位相緩和

$\Phi_p(\rho)=(1-p)\rho+pI/2$、$0\le p\le1$ を:term[脱分極チャネル]{en="depolarizing channel" yomi="だつぶんきょくちゃねる"}とし、$p$ を置換確率とする。状態の Bloch ベクトルは $(1-p)\boldsymbol r$ となる。Pauli 共役を平均すると $I/2$ になるので、$I$ を確率 $1-3p/4$、$X,Y,Z$ を各 $p/4$ で作用させる実装と等しい。例えば $p=1/2$ で $\ket0$ の出力は $\diag(3/4,1/4)$ である。

位相反転チャネル $\Phi_q(\rho)=(1-q)\rho+qZ\rho Z$ は非対角成分を $1-2q$ 倍し、対角成分は保つ。$q=1/2$ で完全に位相が緩和され、$\ket+$ は $I/2$ になる。$q>1/2$ では非対角成分の符号も反転するので、減衰率そのものと反転確率を区別する。

### 高次元の脱分極と完全正値性の範囲

$d\ge2$ とし、$\Phi_a(T)=aT+(1-a)\Tr(T)I/d$ を実数 $a$ について考える。Choi 行列は $a\ketbra\Omega\Omega+(1-a)I_{d^2}/d$ である。$\ket\Omega$ の方向の固有値は $(1+(d^2-1)a)/d$、直交補では $(1-a)/d$ なので、完全正値性は $-1/(d^2-1)\le a\le1$ と同値である。トレース保存は全ての $a$ で成り立つ。一系だけの正性では、純粋入力の固有値 $[1+(d-1)a]/d,(1-a)/d$ から $-1/(d-1)\le a\le1$ で十分である。純粋射影の正の線形結合が全正演算子を尽くすので、この範囲は必要でもある。

任意の $d$ 次元で、巡回シフト $S\ket j=\ket{j+1\bmod d}$ と位相 $Z\ket j=e^{2\pi ij/d}\ket j$ を用いる。$S^uZ^v$ の一様共役では、まず $v$ の平均で非対角成分が消え、次の $u$ の平均で対角成分が一様になる。従って $\Phi_a$ は $I$ を重み $[1+(d^2-1)a]/d^2$、残る $d^2-1$ 個の共役を各 $(1-a)/d^2$ で適用する混合である。この表示では、$a<0$ でも完全正値性の範囲内なら重みが全て非負になる。

:::example[量子ビットで軸が反転する境界]
$d=2,a=-1/3$ なら恒等を使わず $X,Y,Z$ を各 $1/3$ で適用する。入力 $\ket0$ の出力は $\diag(1/3,2/3)$ であり、Bloch ベクトルの向きは反転する。二値の直交基底を入れ替える処理を全ての純粋入力について厳密に行う $a=-1$ は、単独系の正性は満たすが Choi 行列に負固有値がある。
:::

### 振幅減衰

$K_0=\diag(1,\sqrt{1-\gamma})$、$K_1=\sqrt{\gamma}\ketbra01$ とする。:term[振幅減衰チャネル]{en="amplitude damping channel" yomi="しんぷくげんすいちゃねる"}では $\gamma$ が励起の落ちる確率となる。積を展開すると

$$
\begin{pmatrix}a&b\\b^*&d\end{pmatrix}\longmapsto
\begin{pmatrix}a+\gamma d&\sqrt{1-\gamma}b\\\sqrt{1-\gamma}b^*&(1-\gamma)d\end{pmatrix}.

$$

従って Bloch ベクトルは $(\sqrt{1-\gamma}r_x,\sqrt{1-\gamma}r_y,(1-\gamma)r_z+\gamma)$ へ写る。例えば $\gamma=1/4$ で $\ket1$ の出力は $\diag(1/4,3/4)$ である。零温度の指数緩和模型では $\gamma=1-e^{-t/T_1}$ と置き、$t$ と $T_1$ を同じ時間単位で表す。

### 時間発展としての振幅減衰と位相緩和

零温度で再励起がない二準位の緩和模型では、励起確率 $p_1(t)$ を $\dot p_1=-\Gamma_1 p_1$ とする。従って $\gamma(t)=1-e^{-\Gamma_1 t}$ を減衰確率に取る。Kraus 表示の非対角成分は $e^{-\Gamma_1 t/2}$ 倍になり、純粋な位相雑音を独立に重ねて $e^{-\Gamma_\phi t}$ を掛ければ $T_2^{-1}=(2T_1)^{-1}+\Gamma_\phi$ である。ここで $T_1=1/\Gamma_1$ は人口差の減衰時間、$T_2$ はコヒーレンスの減衰時間であり、率の単位は $\mathrm{s^{-1}}$ とする。$\Gamma_\phi\ge0$ の模型では $T_2\le2T_1$ である。

:::example[二つの緩和率の推定]
$T_1=40\,\mu\mathrm s,T_2=50\,\mu\mathrm s$ なら、$\Gamma_\phi=1/(50\,\mu\mathrm s)-1/(80\,\mu\mathrm s)=7500\,\mathrm{s^{-1}}$ である。$20\,\mu\mathrm s$ 後に励起人口は $e^{-1/2}$ 倍、コヒーレンスは $e^{-2/5}$ 倍である。振幅減衰だけならコヒーレンスは $e^{-1/4}$ 倍であり、差として純粋位相緩和の $e^{-3/20}$ が掛かる。
:::

:::example[減衰を連続して作用させる]
$\gamma_1=1/4$、$\gamma_2=1/3$ の二回の振幅減衰では、励起生存確率は $(3/4)(2/3)=1/2$ なので一回の $\gamma=1/2$ と等価である。非対角成分の倍率も $\sqrt{3/4}\sqrt{2/3}=1/\sqrt2$ となって一致する。二回の Kraus リストをそのまま合成すると四項になるが、積 $K_1^{(2)}K_1^{(1)}$ は零で、残る二つのジャンプ演算子は同じ行列に比例する。したがって最小表示は二項に戻る。
:::

### 同じ人口差の入力に種々の模型を掛ける

入力 $\boldsymbol r=(3/5,0,4/5)$ は純粋である。$Z$ 基底の人口と $X$ の平均を両方読めば、人口緩和と位相緩和の違いが確認できる。ここでは各チャネルの作用を Kraus 演算子の和から計算し、出力の行列式と純度も別に求める。

:::example[等方的な縮小]

脱分極チャネルの置換確率を $p=1/2$ とする。

上で定めた入力に以下の変換を掛ける。出力の行列と平均の両方を計算して、古典人口とコヒーレンスの変化を読む。

$$
\begin{aligned}
\rho_{\mathrm{in}}&=\left(\begin{matrix}\frac{9}{10} & \frac{3}{10}\\\frac{3}{10} & \frac{1}{10}\end{matrix}\right)\\
\Phi(\rho_{\mathrm{in}})&=\sum_jK_j\rho_{\mathrm{in}} K_j^\dagger=\left(\begin{matrix}\frac{7}{10} & \frac{3}{20}\\\frac{3}{20} & \frac{3}{10}\end{matrix}\right)\\
(r_x,r_y,r_z)_{\mathrm{out}}&=(\frac{3}{10},0,\frac{2}{5})\\
\det\Phi(\rho_{\mathrm{in}})&=\frac{3}{16}\\
\Tr[\Phi(\rho_{\mathrm{in}})^2]&=(1+\|\boldsymbol r_{\mathrm{out}}\|^2)/2=\frac{5}{8}
\end{aligned}
$$

両成分は同じ割合で縮むため、平均の軸の方向は変わらない。これに対し、人口緩和や位相緩和で特定の軸だけに変化があると方向も動く。
:::

:::example[位相緩和で人口を保持する]

位相反転確率を $q=1/4$ とする。

上で定めた入力に以下の変換を掛ける。出力の行列と平均の両方を計算して、古典人口とコヒーレンスの変化を読む。

$$
\begin{aligned}
\rho_{\mathrm{in}}&=\left(\begin{matrix}\frac{9}{10} & \frac{3}{10}\\\frac{3}{10} & \frac{1}{10}\end{matrix}\right)\\
\Phi(\rho_{\mathrm{in}})&=\sum_jK_j\rho_{\mathrm{in}} K_j^\dagger=\left(\begin{matrix}\frac{9}{10} & \frac{3}{20}\\\frac{3}{20} & \frac{1}{10}\end{matrix}\right)\\
(r_x,r_y,r_z)_{\mathrm{out}}&=(\frac{3}{10},0,\frac{4}{5})\\
\det\Phi(\rho_{\mathrm{in}})&=\frac{27}{400}\\
\Tr[\Phi(\rho_{\mathrm{in}})^2]&=(1+\|\boldsymbol r_{\mathrm{out}}\|^2)/2=\frac{173}{200}
\end{aligned}
$$

入力と同じ $Z$ の平均でも、状態全体は異なる。$X$ が半分になるので異なる固有値を持つ状態として識別できる。
:::

:::example[零温度の振幅減衰で基底へ傾ける]

振幅減衰率を $\gamma=3/4$ とする。

上で定めた入力に以下の変換を掛ける。出力の行列と平均の両方を計算して、古典人口とコヒーレンスの変化を読む。

$$
\begin{aligned}
\rho_{\mathrm{in}}&=\left(\begin{matrix}\frac{9}{10} & \frac{3}{10}\\\frac{3}{10} & \frac{1}{10}\end{matrix}\right)\\
\Phi(\rho_{\mathrm{in}})&=\sum_jK_j\rho_{\mathrm{in}} K_j^\dagger=\left(\begin{matrix}\frac{39}{40} & \frac{3}{20}\\\frac{3}{20} & \frac{1}{40}\end{matrix}\right)\\
(r_x,r_y,r_z)_{\mathrm{out}}&=(\frac{3}{10},0,\frac{19}{20})\\
\det\Phi(\rho_{\mathrm{in}})&=\frac{3}{1600}\\
\Tr[\Phi(\rho_{\mathrm{in}})^2]&=(1+\|\boldsymbol r_{\mathrm{out}}\|^2)/2=\frac{797}{800}
\end{aligned}
$$

励起人口は元の四分の一になり、非対角成分は半分になる。$Z$ の平均は正方向に上がるので、等方的な縮小とは区別できる。
:::

### 有限温度で平衡人口が零にならない模型

環境ビットが $\ket0$ の確率 $q$、$\ket1$ の確率 $1-q$ で上の交換実装に入るとする。環境が $\ket0$ の枝では下向きの減衰、$\ket1$ の枝では上向きの励起となる。混合した変換は励起人口を
$$
p_1\longmapsto (1-\gamma)p_1+(1-q)\gamma
$$
とし、非対角成分を $\sqrt{1-\gamma}$ 倍する。$\gamma>0$ なら固定点は $\diag(q,1-q)$ である。$q=1/2$ なら固定点は一様でユニタルであり、$q=1$ で零温度模型に戻る。連続時間で上下の遷移率を $\Gamma_\uparrow,\Gamma_\downarrow$ に取る人口方程式
$$
\frac{\mathrm d p_1}{\mathrm d t}=\Gamma_\uparrow(1-p_1)-\Gamma_\downarrow p_1
$$
の解も、$\gamma=1-e^{-(\Gamma_\uparrow+\Gamma_\downarrow)t}$、$1-q=\Gamma_\uparrow/(\Gamma_\uparrow+\Gamma_\downarrow)$ とすれば同じ形になる。率は $\mathrm{s^{-1}}$、時間は $\mathrm s$ である。これは独立な位相減衰を追加する前の二準位の緩和模型である。

:::example[有限温度で励起人口が途中で下げ止まる]
$\gamma=1/2,q=9/10$ とし、初期励起人口を 1 とする。1 回目は $(1/2)+(1/10)(1/2)=11/20$、2 回目は $13/40$ になる。一般の $m$ 回後は $p_1=(1/10)+(9/10)2^{-m}$ である。初期状態が $\ket+$ なら、$m$ 回後の非対角成分は $2^{-m/2}/2$ である。従って集団の緩和先は混合し、位相の平均は零になる。環境を同じ温度で独立に準備する仮定を外すと、各段で同じ係数を掛けられるかを再確認する必要がある。
:::

### 位相の確率過程と Lindblad 生成子

幅 $\mathrm dt$ で確率 $\kappa\,\mathrm dt$ の $Z$ 共役を行う模型では、微分方程式は $\dot\rho=\kappa(Z\rho Z-\rho)$ となる。$Z$ の対角成分は不変で、非対角成分は $\dot\rho_{01}=-2\kappa\rho_{01}$ である。従って時間 $t$ 後の位相反転確率は $(1-e^{-2\kappa t})/2$ である。零温度の振幅減衰では $\sigma_-=\ketbra01$ を使い、短時間の Kraus 演算子を $K_0=I-(\Gamma\,\mathrm dt/2)\sigma_+ \sigma_-+O((\mathrm dt)^2)$、$K_1=\sqrt{\Gamma\,\mathrm dt}\sigma_-+O((\mathrm dt)^{3/2})$ と展開できる。$\sigma_+=\sigma_-^\dagger$ である。一次項を取れば
$$
\dot\rho=\Gamma\left(\sigma_-\rho\sigma_+-\tfrac12\{\sigma_+\sigma_-,\rho\}\right)
$$
を得る。跡は第一項と反可換子の跡が一致して保存される。人口は $e^{-\Gamma t}$、コヒーレンスは $e^{-\Gamma t/2}$ を掛けるので、有限時間の Kraus 公式に戻る。ここで扱った生成子は明示的な模型から得たもので、あらゆる開放系の時間発展が同じ形の定数率で書けるわけではない。

## ユニタルチャネルと多数化

$\Phi(I)=I$ の写像を :term[ユニタル]{en="unital" yomi="ゆにたる"}という。同次元トレース保存写像なら一様状態を保つ。実ベクトル $x,y$ を降順に並べ、各上位 $k$ 個の和が $x$ で $y$ 以下、全和が等しいとき $x\prec y$ と書き、:term[多数化]{en="majorization" yomi="たすうか"}の関係にあり、$y$ が $x$ を多数化するという。

:::proposition[ユニタル写像と固有値]{#prop-qi-unital-majorization}
同次元のユニタルチャネルは $\lambda(\Phi(\rho))\prec\lambda(\rho)$ を満たす。
:::

:::proof
$\rho,\Phi(\rho)$ の固有基底をそれぞれ $(u_j),(v_i)$ とし、$B_{ij}=\mel{v_i}{\Phi(\ketbra{u_j}{u_j})}{v_i}$ とする。完全正値性で非負、トレース保存で各列和 1、ユニタル性で各行和 1 となる。出力の固有値は $B$ と入力固有値の積である。上位 $k$ 行を足すと、列ごとの重みは $[0,1]$ にあり重み総和は $k$ である。降順の入力固有値に対する加重和は、先頭 $k$ 個へ重み 1 を集中させた値を超えない。全和はトレース保存から等しい。
:::

ユニタリ共役の確率混合はユニタルである。例えば位相反転では入力の固有値が $(1,0)$ の場合でも、$\ket+$ は完全緩和で $(1/2,1/2)$ になる。一方、振幅減衰は $\Phi_\gamma(I)=\diag(1+\gamma,1-\gamma)$ なので $\gamma>0$ でユニタルではない。


### 量子ビットではユニタル性から確率混合が従う

量子ビットのトレース保存かつ Hermite 保存な写像は実数のアフィン変換 $\boldsymbol r\mapsto M\boldsymbol r+\boldsymbol t$ で表せる。本節ではユニタルなので $\boldsymbol t=0$ である。実特異値分解の直交因子で行列式が $-1$ なら、対角因子に符号を移して、両側を回転行列かつ中央を符号付き対角 $\diag(\lambda_1,\lambda_2,\lambda_3)$ にできる。:ref[prop-qi-bloch-rotation] により両側の回転はユニタリ共役である。

対角の場合は Pauli 共役の線形結合として

$$

p_0=\tfrac14(1+\lambda_1+\lambda_2+\lambda_3),\qquad
p_j=\tfrac14(1+\lambda_j-\lambda_k-\lambda_l)
$$

と書ける。ここで $\{j,k,l\}=\{1,2,3\}$ である。例えば $X$ の共役は $X$ 成分を固定し $Y,Z$ 成分を反転するので、係数を代入すると指定の対角行列になる。この表示の Choi 行列は四つの直交する Pauli ベクトルの外積の和である。したがって完全正値性は各 $p_j\ge0$ と同値であり、総和が 1 なので確率混合が得られる。両側の回転を戻せば、元の写像もユニタリ共役の確率混合になる。

:::example[球を球へ送るだけでは足りない]
$(\lambda_1,\lambda_2,\lambda_3)=(1,-1,1)$ は転置である。Bloch 球を保つので正だが、式の $Y$ 共役の係数は $(1-1-1-1)/4=-1/2$ となる。これは Choi 行列の負固有値と対応する。一方 $(1/2,1/2,0)$ では $p_0=1/2,p_1=p_2=1/4,p_3=0$ となり物理的に実現できる。
:::

### 縦方向の平行移動を含む量子ビットチャネル

Bloch 座標の変換が $(r_x,r_y,r_z)\mapsto(\lambda_1 r_x,\lambda_2 r_y,\lambda_3 r_z+t)$ の形なら、非規格化 Choi 行列は
$$
J=\frac12\begin{pmatrix}
1+t+\lambda_3 &0&0&\lambda_1+\lambda_2\\
0&1-t-\lambda_3&\lambda_1-\lambda_2&0\\
0&\lambda_1-\lambda_2&1+t-\lambda_3&0\\
\lambda_1+\lambda_2&0&0&1-t+\lambda_3
\end{pmatrix}.
$$
実際、$E_{00}=(I+Z)/2$、$E_{11}=(I-Z)/2$ が二つの対角ブロックに対応し、$E_{01}=(X+iY)/2$ が右上ブロックに対応する。偶奇の二つの支持に分けると、完全正値性は二つの $2\times2$ ブロックの対角成分と行列式が非負であることと同値である。つまり $1\pm\lambda_3\ge |t|$ および $(1\pm\lambda_3)^2-t^2\ge(\lambda_1\pm\lambda_2)^2$ が要る。零温度の振幅減衰では $\lambda_1=\lambda_2=\sqrt{1-\gamma}$、$\lambda_3=1-\gamma,t=\gamma$ で、両方の行列式で得た二つの不等式は、ともに等号になる。

### 三次元では同じ結論にならない

:::example[反対称空間に支持を持つチャネル]
$d=3$ で $\Phi(M)=(\Tr(M)I-M^\top)/2$ とする。交換行列を $F$ とすると $J=(I-F)/2$ であり、反対称部分空間への射影なので正である。また、$\Tr\Phi(M)=\Tr(M)$、$\Phi(I)=I$ である。仮にユニタリ共役の凸結合で表せるなら、各ユニタリのベクトル化は $J$ の支持に入る。これは核で二次形式を取ると各項が非負で和が零だからである。ベクトル化が反対称である行列は $U^\top=-U$ を満たす。ところが奇数次元では $\det U=\det(-U)=-\det U$ なので特異であり、ユニタリではありえない。従ってユニタル性だけでユニタリ混合が保証されるのは次元に依存する。
:::

### 出力の固有値と混合性

:::lemma[多数化の凸包表示]{#lem-qi-majorization-convex}
同じ次元の確率ベクトル $a,b$ について、$a\prec b$ と $a$ が $b$ の置換の凸結合であることは同値である。
:::

:::proof
まず $a\prec b$ とし、両ベクトルを降順とする。仮に $a$ がその凸包の外にあれば、凸包への最近点から $a$ を隔てる線形汎関数 $w$ を取れる。逆順の成分を交換するたびに内積が増えるので、$w\cdot a\le w^\downarrow\cdot a$、かつ置換した $b$ との最大内積は $w^\downarrow\cdot b$ である。さらに差は
$$
w^\downarrow\cdot(b-a)=\sum_{k=1}^{d-1}(w^\downarrow_k-w^\downarrow_{k+1})\sum_{j=1}^k(b_j-a_j)\ge0
$$
となる。ここに多数化の部分和条件を使った。従って分離は不可能である。すると対応する置換行列の凸結合を $D$ として $a=Db$ と書ける。逆に $D$ が二重確率なら、ユニタル写像の証明と同じ行和・列和の議論で $Db\prec b$ が成り立つ。置換行列の凸結合は二重確率なので、これで逆向きも示された。
:::

この補題より、$a\prec b$ ならある二重確率行列 $D$ で $a=Db$ と書ける。従って凹関数 $f$ に対し $\sum_i f(a_i)\ge\sum_i f(b_i)$ である。各行で Jensen を適用し、列和が 1 なので右辺を整理できる。これを $f(t)=-t\log_2t$ に適用するとユニタルチャネルは von Neumann エントロピーを減らさず、$f(t)=-t^2$ なら純度を上げないと分かる。一方でユニタル性は、任意の入力が必ずより混合になるという厳密不等式を意味しない。

:::example[純度が保存される入力と縮む入力]
確率 $p=1/4$ の位相反転では $Z$ は保存し、$X,Y$ は半分に縮む。入力 $\ket0$ は固定され、純度 1 を保つ。入力 $\ket+$ の出力は $(I+X/2)/2$ となり固有値 $3/4,1/4$、純度 $5/8$ である。これに対し振幅減衰で入力 $I/2$ を $\gamma=3/4$ で緩和させると $\diag(7/8,1/8)$ となり、純度が $1/2$ から $25/32$ へ上がる。入力のエントロピーが減るかどうかは、トレース保存性だけでは判断できない。
:::

### 人口の二重確率遷移と実装を組み立てる

対角入力を扱うとき、置換ユニタリの確率分布が人口の古典通信路を決める。次の例では $e_i$ の列をどの位置へ移すかで $P$ を定め、二重確率行列 $D$ と固有値の累積和の両方を計算する。$D$ を知っていても、対角以外に同じ作用を持つ量子写像が一意に決まるわけではない。

:::example[三準位の循環だけで第一人口を移す]

複数の置換を独立の古典ラベルで選ぶ。行列 $D$ は各入力基底がどの出力へ移るかを確率として記録したものである。

$$
\begin{aligned}
D&=\left(\begin{matrix}\frac{1}{2} & 0 & \frac{1}{2}\\\frac{1}{2} & \frac{1}{2} & 0\\0 & \frac{1}{2} & \frac{1}{2}\end{matrix}\right)\\
x&=(\frac{3}{5},\frac{3}{10},\frac{1}{10})\\
y&=Dx=(\frac{7}{20},\frac{9}{20},\frac{1}{5})\\
L(x)&=(\frac{3}{5},\frac{9}{10},1)\\
L(y)&=(\frac{9}{20},\frac{4}{5},1)
\end{aligned}
$$

第二人口が最大になり、入力の基底順で見た第一成分と順序が変わる。多数化を調べるときは、出力を降順に並べ直してから累積和を取る。
:::

:::example[一つの人口だけを固定する交換]

複数の置換を独立の古典ラベルで選ぶ。行列 $D$ は各入力基底がどの出力へ移るかを確率として記録したものである。

$$
\begin{aligned}
D&=\left(\begin{matrix}\frac{1}{4} & 0 & \frac{3}{4}\\0 & 1 & 0\\\frac{3}{4} & 0 & \frac{1}{4}\end{matrix}\right)\\
x&=(\frac{1}{2},\frac{1}{3},\frac{1}{6})\\
y&=Dx=(\frac{1}{4},\frac{1}{3},\frac{5}{12})\\
L(x)&=(\frac{1}{2},\frac{5}{6},1)\\
L(y)&=(\frac{5}{12},\frac{3}{4},1)
\end{aligned}
$$

第二の人口は入力と一致するが、別の二成分には平均化がかかる。最大値の減少だけでなく、上位二つの和も比較すると多数化の全条件を確認できる。
:::

:::example[全置換でなく循環の平均で一様にする]

複数の置換を独立の古典ラベルで選ぶ。行列 $D$ は各入力基底がどの出力へ移るかを確率として記録したものである。

$$
\begin{aligned}
D&=\left(\begin{matrix}\frac{1}{3} & \frac{1}{3} & \frac{1}{3}\\\frac{1}{3} & \frac{1}{3} & \frac{1}{3}\\\frac{1}{3} & \frac{1}{3} & \frac{1}{3}\end{matrix}\right)\\
x&=(\frac{3}{4},\frac{1}{4},0)\\
y&=Dx=(\frac{1}{3},\frac{1}{3},\frac{1}{3})\\
L(x)&=(\frac{3}{4},1,1)\\
L(y)&=(\frac{1}{3},\frac{2}{3},1)
\end{aligned}
$$

零の人口も他の位置の重みを受け取って正になる。対角入力が一様になるという結論から、非対角成分まで消えるとは言えない。たとえば三成分が等しい純粋ベクトルは循環置換で不変である。
:::

### 二状態の順序を量子チャネルで実装する

$\lambda(\sigma)\prec\lambda(\rho)$ であるとする。多数化の凸包表示から、$\lambda(\sigma)$ は $\lambda(\rho)$ の置換ベクトルの凸結合として書ける。各置換を $\rho$ の固有基底での置換ユニタリに置き換える。混合の出力はその固有基底で $\lambda(\sigma)$ を対角に持つので、最後に $\sigma$ の固有基底までユニタリを掛ければ $\sigma$ になる。従って一組の $\rho,\sigma$ について、あるユニタリ混合が $\rho$ を $\sigma$ へ送ることと、あるユニタルチャネルがそうすることと、$\lambda(\sigma)\prec\lambda(\rho)$ は同値である。また全状態を多数化の向きに送るチャネルは、一様入力を固定するのでユニタルである。一方、既に固定されたユニタル写像が全入力に対して単一のリストで確率混合表示を持つかどうかは別の問いである。

:::example[三次元の二状態だけを合わせる]
$\rho=\diag(3/5,3/10,1/10)$ と $\sigma=\diag(9/20,9/20,1/10)$ を取る。$\sigma$ の最上位は $9/20\le3/5$、上位二つの和はともに $9/10$ なので多数化条件を満たす。第一と第二の基底を交換するユニタリ $U$ と恒等を半々で混ぜれば、その二人口はどちらも $(3/5+3/10)/2=9/20$ になり、第三は変わらない。逆方向では最上位和を増やす必要があるので、ユニタル写像で実現できない。非ユニタルな置換チャネルなら $\rho$ を準備し直すことで逆方向を実現できる。
:::

### 二重確率行列の置換分解

行列 $D=(D_{ij})$ の行和と列和が 1、成分が非負とする。正成分だけを辺とする二部グラフを作る。任意の行集合 $S$ から出る重みの和は $|S|$ であり、隣接する列集合 $N(S)$ の全重みは $|N(S)|$ なので $|N(S)|\ge|S|$ である。この条件で完全マッチングが存在することを、最大マッチングの増加経路で示す。未対応の行があるとき、そこから行から列へは未使用辺、列から行へは対応辺をたどる。未対応の列へ届けば経路上の使用と未使用を入れ替えて対応数を増やせる。届かない場合、到達可能な列にはすべて相手の行があり、これらの行と出発点を合わせた集合は到達した列より一つ多い。またその行集合の隣列は全て到達した列に入る。対応辺も当該行へ入る際に辿っているからである。従って $|N(S)|\ge |S|$ に反する。

選んだ置換行列 $P$ に対し、対応辺での最小成分を $s>0$ とおく。$s=1$ なら $D=P$ である。そうでなければ $(D-sP)/(1-s)$ は再び二重確率であり、少なくとも一つの正成分が零になる。同じ操作は有限回で終わる。従って $D$ は置換行列の凸結合である。この手順は状態の人口の配分を構成する際に使える。非対角成分への作用を合わせたい場合には、同じ変換がコヒーレンスも所望の形へ送るかを別に計算する。

## 演習問題

::::exercise{level=basic}
恒等チャネルの Choi 行列の階数を求めよ。

:::solution
$J=\ketbra\Omega\Omega$ なので非零固有値は一つであり、階数 1 である。
:::
::::

::::exercise{level=basic}
$\Phi(M)=\Tr(M)I_2/2$ の Choi 行列を求めよ。

:::solution
行列単位は $\delta_{ij}I_2/2$ に写り、ブロックを合わせると $I_4/2$ である。
:::
::::

::::exercise{level=basic}
位相反転の $q=1$ で $\ket+$ はどうなるか。

:::solution
$Z$ が確率 1 で働くので、符号の異なる $\ket-$ に移る。
:::
::::

::::exercise{level=basic}
$\gamma=1/2$ の振幅減衰を $\ket1$ に施した出力を求めよ。

:::solution
対角成分は $\gamma,1-\gamma$ なので $I/2$ である。
:::
::::

::::exercise{level=advanced}
二つの Kraus 表示 $(K_a),(L_b)$ を持つ連続作用の Kraus 表示を求めよ。

:::solution
$\sum_bL_b(\sum_aK_a M K_a^\dagger)L_b^\dagger$ を展開して $(L_bK_a)_{b,a}$ を得る。両完備性を順に使うと合成の完備性も成立する。
:::
::::

::::exercise{level=advanced}
$\gamma=1/4$ の振幅減衰の Choi 行列の正固有値を求めよ。

:::solution
二つの Kraus 演算子はトレース内積で直交するので、それらのベクトル化の二乗ノルムが固有値である。$K_0$ から $1+3/4=7/4$、$K_1$ から $1/4$ を得る。
:::
::::

::::exercise{level=advanced}
$\lambda=(0.7,0.2,0.1)$ は $(0.6,0.3,0.1)$ を多数化するか。

:::solution
上位 1 個は $0.7\ge0.6$、上位 2 個はともに $0.9$、全和はともに 1 である。従って多数化する。
:::
::::

::::exercise{level=advanced}
振幅減衰の相補チャネルが振幅減衰率 $1-\gamma$ のチャネルに一致することを示せ。

:::solution
$V\ket0=\ket{00}$、$V\ket1=\sqrt{1-\gamma}\ket{10}+\sqrt\gamma\ket{01}$ と置く。$B$ を捨てると励起人口は $\gamma$ 倍、非対角成分は $\sqrt\gamma$ 倍となるので、環境出力への減衰率は $1-\gamma$ である。
:::
::::

::::exercise{level=advanced}
一量子ビットの写像 $\Phi$ が $\Phi(I)=I$、$\Phi(X)=X/2$、$\Phi(Y)=Y/2$、$\Phi(Z)=0$ を満たすとする。完全正値性を確かめ、最小数の Kraus 演算子を求めよ。また $\ket+$ と $\ket0$ の出力の純度を比較せよ。

:::solution
任意の Pauli 混合 $\sum_{j=0}^3p_jP_jT P_j^\dagger$ では、各因子と同じ Pauli はその項で固有値 $+1$、反可換な二つは $-1$ になる。そのため Bloch の三倍率は順に $p_0+p_1-p_2-p_3$、$p_0-p_1+p_2-p_3$、$p_0-p_1-p_2+p_3$ である。総和 1 を加えた四元連立方程式を解くと $(p_0,p_1,p_2,p_3)=(1/2,1/4,1/4,0)$ である。三つの正重みを使い、$I/\sqrt2,X/2,Y/2$ が Kraus 演算子になる。三つのベクトル化した演算子はトレース内積で直交するので、Choi 行列は階数 3 である。従って三つより少ない個数では実現できない。

$\ket+$ は $X$ の正固有状態であり、出力の Bloch ベクトルは $(1/2,0,0)$ になる。固有値は $3/4,1/4$ なので純度は $5/8$ である。$\ket0$ の Bloch ベクトルは $Z$ 軸であり、出力は $I/2$ になる。純度は $1/2$ である。入力の純度が同じでも、チャネルの軸の選び方により結果は異なる。
:::
::::

::::exercise{level=advanced}
位相反転 $\mathcal P_p(T)=(1-p)T+pZTZ$ を二回作用させる。$p=1/4$ のときの相補チャネルで、二つの計算基底入力から得られる純粋環境状態の重なりを求めよ。環境を一つの二次元系にまとめられる理由も説明せよ。

:::solution
各回の環境ベクトルを $e_j=\sqrt{1-p}\ket0+(-1)^j\sqrt p\ket1$ とすると、入力 $\ket j$ は変わらず環境へ $e_j$ が残る。二回に独立な環境を用いれば $e_j\otimes e_j$ が残り、重なりは $\braket{e_0}{e_1}^2=(1-2p)^2=1/4$ になる。両ベクトルの張る空間は二次元なので、環境をその支持へ制限しても内積と全入力への作用を保てる。

同じ結論を受信側から求めると、二反転で元に戻るため奇数反転の確率は $2p(1-p)=3/8$ である。したがって単一環境を使う等長表示の二つのベクトルの重なりは $1-2(3/8)=1/4$ であり、同じ Gram 行列を持つ。二つの環境表示の間で張る空間を対応させれば等長になる。
:::
::::

## 確認テスト

::::quiz
Choi 行列が正であるとき保証される性質を選べ。

- [x] 完全正値性
- [ ] トレース保存性

:::explanation
正固有値の分解から Kraus 表示を作れる。トレース保存には別途 $\Tr_BJ=I_R$ が要る。
:::
::::

::::quiz{answer="4" tol="1e-6"}
量子ビットの完全混合への置換チャネルに必要な最小環境次元はいくつか。

:::explanation
$J=I_4/2$ の階数が 4 なので、:ref[thm-qi-choi] により最小環境次元は 4 である。
:::
::::

::::quiz{answer="0.75" tol="1e-6"}
脱分極チャネルで置換確率が $1/2$ のとき、入力 $\ket0$ を計算基底で測って 0 を得る出力確率を求めよ。

:::explanation
保持分岐から $1/2$、置換分岐から $(1/2)(1/2)=1/4$ を得るので合計 $3/4$ である。
:::
::::

::::quiz
すべての正写像が完全正値であるか。

- [ ] はい。
- [x] いいえ。

:::explanation
転置は一系では非負固有値を保つ。Bell 射影に片側転置を施すと負固有値が生じるので完全正値ではない。
:::
::::

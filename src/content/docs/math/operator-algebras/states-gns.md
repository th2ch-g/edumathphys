---
title: 状態と GNS 表現
description: 正値汎関数から巡回表現を構成し、純粋性、忠実表現、状態空間の凸構造を証明する。
notation:
  - symbol: 'S(A)'
    meaning: A の状態空間
  - symbol: 'N_\varphi'
    meaning: GNS 構成で零長となる左イデアル
  - symbol: '(\pi_\varphi,H_\varphi,\xi_\varphi)'
    meaning: 状態の GNS 表現と巡回単位ベクトル
---

## この章の目標

- 正値汎関数と状態を定義し、有界性を示す。
- GNS 構成を計算し、純粋状態と既約性の関係を証明する。
- 状態を使って抽象 C*-環を作用素として忠実に表す。
- 状態空間の極点を使って正元のノルムを検出する。

## 正値線形汎関数と状態

観測量を元 $a$ で表すとき、期待値には和とスカラー倍の規則があり、正元では非負の値を取る。規格化には単位元の値または汎関数のノルムを使う。

:::definition[正値汎関数と状態]{#def-oa-state}
複素線形汎関数 $\varphi:A\to\C$ が各 $a$ について $\varphi(a^*a)\ge0$ を満たすとき :term[正値汎関数]{en="positive linear functional" yomi="せいちせんけいはんかんすう"} という。ノルムが $1$ の正値汎関数を :term[状態]{en="state" yomi="じょうたい"} という。
:::

自己共役元は二つの正元の差なので、その値は実数である。自己共役成分に分けると $\varphi(a^*)=\overline{\varphi(a)}$ となる。正値汎関数は非単位的な場合も自動的に有界である。正縮小元上で上に有界でないなら、正縮小元 $a_n$ を $\varphi(a_n)\ge n2^n$ となるように選べる。正錐の閉性から $\sum_{n\ge1}2^{-n}a_n$ はどの $2^{-n}a_n$ よりも大きいので、その値がどの $n$ 以上にもなり矛盾する。自己共役元の正負分解と一般元の自己共役成分への分解により、有界性が正縮小元から単位球に広がる。

$\varphi(a^*a)>0$ がすべての $a\ne0$ に対して成り立つ状態を :term[忠実状態]{en="faithful state" yomi="ちゅうじつじょうたい"} という。これはすべての非零正元に正の値を取ることと同値である。また $\varphi(ab)=\varphi(ba)$ がすべての $a,b$ について成り立つ状態を :term[トレース状態]{en="tracial state" yomi="とれーすじょうたい"} という。

:::example[行列の状態]
$D=\diag(1/4,3/4)$ とし、$\varphi(a)=\tr(Da)$ とする。$\varphi(a^*a)=\norm{ae_1}^2/4+3\norm{ae_2}^2/4$ である。$a=E_{11}$ の期待値は $1/4$、$a=E_{12}$ の期待値は $0$ となる。一般の $M_n$ の正値汎関数は一意に $\tr(Da)$ と書ける。行列単位で $D_{ji}=\varphi(E_{ij})$ と定め、階数一射影上の値を調べると $D\ge0$ である。状態条件は $\tr D=1$ である。
:::

### 正値汎関数の密度行列を復元する

$M_d$ の任意の線形汎関数 $\varphi$ に対し、$D_{ij}=\varphi(E_{ji})$ と置く。任意の $a=\sum a_{ij}E_{ij}$ に対し $\varphi(a)=\tr(Da)$ である。$\varphi$ が正なら $\varphi(xx^*)=x^*Dx\ge0$ がすべての列 $x$ に対して成り立つ。偏極で $D=D^*$ が従い、固有ベクトルを $x$ に取ると各固有値が非負である。逆に $D\ge0$ なら $\tr(Da^*a)=\tr((aD^{1/2})^*aD^{1/2})\ge0$ である。従って状態となる条件は $D\ge0$ と $\tr D=1$ である。

:::example[位相を観測する状態]
単位列 $v=(1,i)^{\mathsf T}/\sqrt2$ でベクトル状態を作る。$D=vv^*=\frac12\begin{pmatrix}1&-i\\i&1\end{pmatrix}$ であり、$\varphi(E_{12})=D_{21}=i/2$ となる。行列 $\diag(1/2,1/2)$ に替えても $\varphi(E_{11})=\varphi(E_{22})=1/2$ は同じだが、$\varphi(E_{12})$ は零になる。対角観測だけでは二つの状態を区別できないので、非対角成分を含む $-iE_{12}+iE_{21}$ を測ると、それぞれ $1$ と $0$ が得られる。
:::

:::example[全点で正の重みと忠実性]
$\varphi((a_n))=\sum_{n\ge1}2^{-n}a_n$ は $c_0$ の正値汎関数である。最初の $N$ 座標で $1$ を取る射影の値は $1-2^{-N}$ なので $\norm\varphi=1$ である。さらに $\varphi(a^*a)=0$ なら、非負の和の各項が零でありすべての $a_n=0$ となる。従ってこの状態は忠実である。どれか一つの重みが零ならその座標の射影を用いて忠実でないことが分かる。
:::

:::example[一座標に集中する状態]
列 $v$ が定める規格化ベクトル状態を $\varphi$ とする。密度行列を使って自己共役行列 $h$ の期待値と分散を計算し、ベクトルによる式と照合する。

まず $D=vv^*/(v^*v)$ を作る。トレースは $1$、階数は一である。$h$ と積を取ると非対角成分も期待値に寄与することが分かる。

$$
\begin{aligned}
v&=\begin{pmatrix}1\\0\end{pmatrix}\\
D&=\begin{pmatrix}1 & 0\\0 & 0\end{pmatrix}\\
h&=\begin{pmatrix}2 & 1\\1 & 3\end{pmatrix}
\end{aligned}
$$

$\varphi(h)=\tr(Dh)$ で対角和を取り、$h^2$ とその期待値にも同じことをする。$\varphi(h^2)$ と $\varphi(h)^2$ は別に計算する。

$$
\begin{aligned}
Dh&=\begin{pmatrix}2 & 1\\0 & 0\end{pmatrix}\\
\varphi(h)&=2\\
h^2&=\begin{pmatrix}5 & 5\\5 & 10\end{pmatrix}\\
\varphi(h^2)&=5
\end{aligned}
$$

差が分散である。平均を引いた行列を $v$ に掛け、そのベクトルの二乗ノルムを $v^*v$ で割ってもよい。

$$
\begin{aligned}
(h-\varphi(h)I)v&=\begin{pmatrix}0\\1\end{pmatrix}\\
\operatorname{Var}_\varphi(h)&=1
\end{aligned}
$$

分散が零である場合に限り、状態ベクトルは $h$ の固有ベクトルになる。これは $h$ の全固有値が同じであることを要求しない。
:::

:::example[複素位相の相関を含む状態]
列 $v$ が定める規格化ベクトル状態を $\varphi$ とする。密度行列を使って自己共役行列 $h$ の期待値と分散を計算し、ベクトルによる式と照合する。

まず $D=vv^*/(v^*v)$ を作る。トレースは $1$、階数は一である。$h$ と積を取ると非対角成分も期待値に寄与することが分かる。

$$
\begin{aligned}
v&=\begin{pmatrix}1\\i\end{pmatrix}\\
D&=\begin{pmatrix}\frac{1}{2} & - \frac{i}{2}\\\frac{i}{2} & \frac{1}{2}\end{pmatrix}\\
h&=\begin{pmatrix}1 & - i\\i & 2\end{pmatrix}
\end{aligned}
$$

$\varphi(h)=\tr(Dh)$ で対角和を取り、$h^2$ とその期待値にも同じことをする。$\varphi(h^2)$ と $\varphi(h)^2$ は別に計算する。

$$
\begin{aligned}
Dh&=\begin{pmatrix}1 & - \frac{3 i}{2}\\i & \frac{3}{2}\end{pmatrix}\\
\varphi(h)&=\frac{5}{2}\\
h^2&=\begin{pmatrix}2 & - 3 i\\3 i & 5\end{pmatrix}\\
\varphi(h^2)&=\frac{13}{2}
\end{aligned}
$$

差が分散である。平均を引いた行列を $v$ に掛け、そのベクトルの二乗ノルムを $v^*v$ で割ってもよい。

$$
\begin{aligned}
(h-\varphi(h)I)v&=\begin{pmatrix}- \frac{1}{2}\\\frac{i}{2}\end{pmatrix}\\
\operatorname{Var}_\varphi(h)&=\frac{1}{4}
\end{aligned}
$$

分散が零である場合に限り、状態ベクトルは $h$ の固有ベクトルになる。これは $h$ の全固有値が同じであることを要求しない。
:::

:::example[符号の異なる観測値]
列 $v$ が定める規格化ベクトル状態を $\varphi$ とする。密度行列を使って自己共役行列 $h$ の期待値と分散を計算し、ベクトルによる式と照合する。

まず $D=vv^*/(v^*v)$ を作る。トレースは $1$、階数は一である。$h$ と積を取ると非対角成分も期待値に寄与することが分かる。

$$
\begin{aligned}
v&=\begin{pmatrix}1\\2\end{pmatrix}\\
D&=\begin{pmatrix}\frac{1}{5} & \frac{2}{5}\\\frac{2}{5} & \frac{4}{5}\end{pmatrix}\\
h&=\begin{pmatrix}2 & 3\\3 & -1\end{pmatrix}
\end{aligned}
$$

$\varphi(h)=\tr(Dh)$ で対角和を取り、$h^2$ とその期待値にも同じことをする。$\varphi(h^2)$ と $\varphi(h)^2$ は別に計算する。

$$
\begin{aligned}
Dh&=\begin{pmatrix}\frac{8}{5} & \frac{1}{5}\\\frac{16}{5} & \frac{2}{5}\end{pmatrix}\\
\varphi(h)&=2\\
h^2&=\begin{pmatrix}13 & 3\\3 & 10\end{pmatrix}\\
\varphi(h^2)&=13
\end{aligned}
$$

差が分散である。平均を引いた行列を $v$ に掛け、そのベクトルの二乗ノルムを $v^*v$ で割ってもよい。

$$
\begin{aligned}
(h-\varphi(h)I)v&=\begin{pmatrix}6\\-3\end{pmatrix}\\
\operatorname{Var}_\varphi(h)&=9
\end{aligned}
$$

分散が零である場合に限り、状態ベクトルは $h$ の固有ベクトルになる。これは $h$ の全固有値が同じであることを要求しない。
:::

:::example[非対角成分だけを測る観測量]
列 $v$ が定める規格化ベクトル状態を $\varphi$ とする。密度行列を使って自己共役行列 $h$ の期待値と分散を計算し、ベクトルによる式と照合する。

まず $D=vv^*/(v^*v)$ を作る。トレースは $1$、階数は一である。$h$ と積を取ると非対角成分も期待値に寄与することが分かる。

$$
\begin{aligned}
v&=\begin{pmatrix}2\\1\end{pmatrix}\\
D&=\begin{pmatrix}\frac{4}{5} & \frac{2}{5}\\\frac{2}{5} & \frac{1}{5}\end{pmatrix}\\
h&=\begin{pmatrix}0 & 1\\1 & 0\end{pmatrix}
\end{aligned}
$$

$\varphi(h)=\tr(Dh)$ で対角和を取り、$h^2$ とその期待値にも同じことをする。$\varphi(h^2)$ と $\varphi(h)^2$ は別に計算する。

$$
\begin{aligned}
Dh&=\begin{pmatrix}\frac{2}{5} & \frac{4}{5}\\\frac{1}{5} & \frac{2}{5}\end{pmatrix}\\
\varphi(h)&=\frac{4}{5}\\
h^2&=\begin{pmatrix}1 & 0\\0 & 1\end{pmatrix}\\
\varphi(h^2)&=1
\end{aligned}
$$

差が分散である。平均を引いた行列を $v$ に掛け、そのベクトルの二乗ノルムを $v^*v$ で割ってもよい。

$$
\begin{aligned}
(h-\varphi(h)I)v&=\begin{pmatrix}- \frac{3}{5}\\\frac{6}{5}\end{pmatrix}\\
\operatorname{Var}_\varphi(h)&=\frac{9}{25}
\end{aligned}
$$

分散が零である場合に限り、状態ベクトルは $h$ の固有ベクトルになる。これは $h$ の全固有値が同じであることを要求しない。
:::

## 正値汎関数に対する Cauchy–Schwarz の不等式

:::theorem[Cauchy–Schwarz と規格化]{#thm-oa-positive-cs}
正値汎関数について
$$
\abs{\varphi(b^*a)}^2\le\varphi(a^*a)\varphi(b^*b).
$$
単位的なら $\norm\varphi=\varphi(1)$ である。非単位的な場合も、正縮小近似単位元 $(e_\lambda)$ について $\varphi(e_\lambda),\varphi(e_\lambda^2)\to\norm\varphi$ となる。
:::
:::proof
$\varphi((a+zb)^*(a+zb))$ を展開する。$\varphi(b^*b)>0$ なら $z=-\varphi(b^*a)/\varphi(b^*b)$ を入れて不等式を得る。$\varphi(b^*b)=0$ のときは $z$ の位相と大きさを変えて負値を防ぐには $\varphi(b^*a)=0$ が必要である。単位的なら $b=1$ と $a^*a\le\norm a^2 1$ により $\abs{\varphi(a)}\le\norm a\varphi(1)$ であり、$a=1$ で等号である。

非単位的な場合は $\abs{\varphi(ae_\lambda)}^2\le\norm a^2\norm\varphi\,\varphi(e_\lambda^2)$ である。左辺は $\abs{\varphi(a)}^2$ に収束する。$\norm a\le1$ の上限から $\liminf_\lambda\varphi(e_\lambda^2)\ge\norm\varphi$ が従い、$\varphi(e_\lambda^2)\le\varphi(e_\lambda)\le\norm\varphi$ で結論を得る。$\varphi=0$ の場合は式が直ちに成立する。
:::

従って非単位的 $A$ の状態は $\widetilde\varphi(a+\lambda1)=\varphi(a)+\lambda$ で単位化に状態として延長できる。実際、$(a+\lambda1)^*(a+\lambda1)$ の延長値は $\varphi((a+\lambda e_\nu)^*(a+\lambda e_\nu))$ の極限なので非負である。規格化から延長は一意となる。

:::example[零長の元]
$M_2$ の $\varphi(a)=a_{11}$ では $\varphi(a^*a)=\abs{a_{11}}^2+\abs{a_{21}}^2$ である。$E_{12}$ の長さは零だが $\varphi(E_{12}E_{12}^*)=1$ である。この差は、以下の零空間が一般に右イデアルにならないことに現れる。
:::

### 等号と分散

Cauchy–Schwarz の差を調べると、商空間で二つのベクトルが比例するかどうかが分かる。$\varphi(b^*b)>0$ のとき、$c=\varphi(b^*a)/\varphi(b^*b)$ に対して
$$
\varphi((a-cb)^*(a-cb))
=\varphi(a^*a)-\frac{|\varphi(b^*a)|^2}{\varphi(b^*b)}
$$
である。従って等号は $a-cb$ が零長であることと同値である。状態 $\varphi$ と単位的な環で $b=1$ を取れば、左辺は $a$ の平均 $\varphi(a)$ を引いた分散である。自己共役 $a$ で分散零とは、状態が $a$ の一つの測定値に集中することを意味する。

:::example[非可換な環の中で測定値が確定する]
$\varphi(a)=a_{11}$、$h=\diag(2,-1)$ とすると $\varphi(h^2)=4$、$\varphi(h)^2=4$ で分散は零である。$k=E_{12}+E_{21}$ に替えると $k^2=I$、$\varphi(k)=0$ なので分散は $1$ である。二つの成分を等確率で取る状態 $\psi=\tr/2$ では $\psi(h)=1/2$、$\psi(h^2)=5/2$ なので $h$ の分散も $9/4$ となる。分散は元と状態の両方に依存する。
:::

:::example[状態を変えずに測定の余りを評価する]
状態 $\varphi$ に対し $p$ を射影とし、$\varphi(p)=1-\delta$ とする。$\norm a\le1$ なら
$$
\begin{aligned}
|\varphi(a)-\varphi(pap)|
&\le|\varphi((1-p)a)|+|\varphi(pa(1-p))|\\
&\le\sqrt{\delta}+\sqrt{\delta}
\end{aligned}
$$
である。最初の項の Cauchy–Schwarz では $\varphi(a^*a)\le1$ を使う。二つ目では $\varphi(pa(1-p))=\overline{\varphi((1-p)a^*p)}$ と形を直し、$\varphi(pa a^*p)\le1$ を使う。例えば $\varphi(p)=99/100$ の場合には誤差は高々 $1/5$ である。$\delta=0$ なら圧縮によって状態の全ての値が保たれる。
:::

:::example[左零イデアルによって等号が出る]
$\varphi(x)=\tr(Dx)$ に対する Cauchy–Schwarz の左辺と右辺を比べる。等号が出た場合、比例しているのが行列自身か商空間の類かを区別する。

自己積と混合積を作る。混合積の順序は $\varphi(b^*a)$ で、第 1 引数を線形とする規約に合う。

$$
\begin{aligned}
D&=\begin{pmatrix}1 & 0\\0 & 0\end{pmatrix}\\
a&=\begin{pmatrix}1 & 3\\2 & 4\end{pmatrix}\\
b&=\begin{pmatrix}2 & 0\\4 & 1\end{pmatrix}\\
b^*a&=\begin{pmatrix}10 & 22\\2 & 4\end{pmatrix}
\end{aligned}
$$

対角密度では各列の内積を列の重みで足すことになる。$\varphi(b^*b)$ は正なので、比例成分の係数は分母で割って求められる。

$$
\begin{aligned}
\varphi(a^*a)&=5\\
\varphi(b^*b)&=20\\
\varphi(b^*a)&=10\\
c=\varphi(b^*a)/\varphi(b^*b)&=\frac{1}{2}
\end{aligned}
$$

$a-cb$ の零長性を調べる。差の行列が非零でも、正の重みが付いた列をすべて消していれば商では零になる。

$$
\begin{aligned}
a-cb&=\begin{pmatrix}0 & 3\\0 & \frac{7}{2}\end{pmatrix}\\
\varphi((a-cb)^*(a-cb))&=0
\end{aligned}
$$

残差の二乗が正なら不等号は真である。零なら $[a]=c[b]$ であり、$D$ が可逆な場合は元の行列の比例も従う。
:::

:::example[忠実な状態で等号が破れる]
$\varphi(x)=\tr(Dx)$ に対する Cauchy–Schwarz の左辺と右辺を比べる。等号が出た場合、比例しているのが行列自身か商空間の類かを区別する。

自己積と混合積を作る。混合積の順序は $\varphi(b^*a)$ で、第 1 引数を線形とする規約に合う。

$$
\begin{aligned}
D&=\begin{pmatrix}\frac{1}{2} & 0\\0 & \frac{1}{2}\end{pmatrix}\\
a&=\begin{pmatrix}1 & 2\\3 & 4\end{pmatrix}\\
b&=\begin{pmatrix}2 & 0\\1 & -1\end{pmatrix}\\
b^*a&=\begin{pmatrix}5 & 8\\-3 & -4\end{pmatrix}
\end{aligned}
$$

対角密度では各列の内積を列の重みで足すことになる。$\varphi(b^*b)$ は正なので、比例成分の係数は分母で割って求められる。

$$
\begin{aligned}
\varphi(a^*a)&=15\\
\varphi(b^*b)&=3\\
\varphi(b^*a)&=\frac{1}{2}\\
c=\varphi(b^*a)/\varphi(b^*b)&=\frac{1}{6}
\end{aligned}
$$

$a-cb$ の零長性を調べる。差の行列が非零でも、正の重みが付いた列をすべて消していれば商では零になる。

$$
\begin{aligned}
a-cb&=\begin{pmatrix}\frac{2}{3} & 2\\\frac{17}{6} & \frac{25}{6}\end{pmatrix}\\
\varphi((a-cb)^*(a-cb))&=\frac{179}{12}
\end{aligned}
$$

残差の二乗が正なら不等号は真である。零なら $[a]=c[b]$ であり、$D$ が可逆な場合は元の行列の比例も従う。
:::

:::example[二つの方向が直交する]
$\varphi(x)=\tr(Dx)$ に対する Cauchy–Schwarz の左辺と右辺を比べる。等号が出た場合、比例しているのが行列自身か商空間の類かを区別する。

自己積と混合積を作る。混合積の順序は $\varphi(b^*a)$ で、第 1 引数を線形とする規約に合う。

$$
\begin{aligned}
D&=\begin{pmatrix}\frac{2}{3} & 0\\0 & \frac{1}{3}\end{pmatrix}\\
a&=\begin{pmatrix}1 & 0\\0 & 0\end{pmatrix}\\
b&=\begin{pmatrix}0 & 1\\1 & 0\end{pmatrix}\\
b^*a&=\begin{pmatrix}0 & 0\\1 & 0\end{pmatrix}
\end{aligned}
$$

対角密度では各列の内積を列の重みで足すことになる。$\varphi(b^*b)$ は正なので、比例成分の係数は分母で割って求められる。

$$
\begin{aligned}
\varphi(a^*a)&=\frac{2}{3}\\
\varphi(b^*b)&=1\\
\varphi(b^*a)&=0\\
c=\varphi(b^*a)/\varphi(b^*b)&=0
\end{aligned}
$$

$a-cb$ の零長性を調べる。差の行列が非零でも、正の重みが付いた列をすべて消していれば商では零になる。

$$
\begin{aligned}
a-cb&=\begin{pmatrix}1 & 0\\0 & 0\end{pmatrix}\\
\varphi((a-cb)^*(a-cb))&=\frac{2}{3}
\end{aligned}
$$

残差の二乗が正なら不等号は真である。零なら $[a]=c[b]$ であり、$D$ が可逆な場合は元の行列の比例も従う。
:::

:::example[複素定数で比例する]
$\varphi(x)=\tr(Dx)$ に対する Cauchy–Schwarz の左辺と右辺を比べる。等号が出た場合、比例しているのが行列自身か商空間の類かを区別する。

自己積と混合積を作る。混合積の順序は $\varphi(b^*a)$ で、第 1 引数を線形とする規約に合う。

$$
\begin{aligned}
D&=\begin{pmatrix}\frac{1}{3} & 0\\0 & \frac{2}{3}\end{pmatrix}\\
a&=\begin{pmatrix}1 & 0\\0 & i\end{pmatrix}\\
b&=\begin{pmatrix}i & 0\\0 & -1\end{pmatrix}\\
b^*a&=\begin{pmatrix}- i & 0\\0 & - i\end{pmatrix}
\end{aligned}
$$

対角密度では各列の内積を列の重みで足すことになる。$\varphi(b^*b)$ は正なので、比例成分の係数は分母で割って求められる。

$$
\begin{aligned}
\varphi(a^*a)&=1\\
\varphi(b^*b)&=1\\
\varphi(b^*a)&=- i\\
c=\varphi(b^*a)/\varphi(b^*b)&=- i
\end{aligned}
$$

$a-cb$ の零長性を調べる。差の行列が非零でも、正の重みが付いた列をすべて消していれば商では零になる。

$$
\begin{aligned}
a-cb&=\begin{pmatrix}0 & 0\\0 & 0\end{pmatrix}\\
\varphi((a-cb)^*(a-cb))&=0
\end{aligned}
$$

残差の二乗が正なら不等号は真である。零なら $[a]=c[b]$ であり、$D$ が可逆な場合は元の行列の比例も従う。
:::

## GNS 構成

C*-環 $A$ から $B(H)$ への *-準同型を :term[表現]{en="representation" yomi="ひょうげん"} という。$\pi(A)$ が作用して得られるベクトルの線形包が $H$ に稠密であるとき :term[非退化表現]{en="nondegenerate representation" yomi="ひたいかひょうげん"} と呼ぶ。$A$ が単位的なら、非退化性は $\pi(1)=1_H$ と同値である。実際、$\pi(1)$ はすべての $\pi(a)\eta$ を固定する射影であり、線形包が稠密なら $1_H$ になる。逆向きは $\pi(1)H=H$ から従う。

:::theorem[GNS 構成]{#thm-oa-gns}
状態 $\varphi$ には、表現 $\pi_\varphi:A\to B(H_\varphi)$ と単位ベクトル $\xi_\varphi$ が対応し、
$$
\varphi(a)=\inner{\pi_\varphi(a)\xi_\varphi}{\xi_\varphi},\qquad
\overline{\pi_\varphi(A)\xi_\varphi}=H_\varphi
$$
となる。この組は指定したベクトルを保つユニタリ同値を除いて一意である。この表現を :term[GNS 表現]{en="GNS representation" yomi="じーえぬえすひょうげん"}、稠密性を満たすベクトルを :term[巡回ベクトル]{en="cyclic vector" yomi="じゅんかいべくとる"} という。
:::

:::proof
単位的な場合、$N_\varphi=\{a:\varphi(a^*a)=0\}$ とする。Cauchy–Schwarz により $N_\varphi$ は半内積 $\varphi(b^*a)$ の零空間であり、線形部分空間である。$\varphi(a^*c^*ca)\le\norm c^2\varphi(a^*a)$ なので左イデアルである。また写像 $a\mapsto\varphi(a^*a)$ はノルム連続なので $N_\varphi$ は閉じている。商の内積を $\inner{[a]}{[b]}=\varphi(b^*a)$ とし完備化する。左乗法 $\pi_\varphi(c)[a]=[ca]$ は先の不等式で有界に延び、積を保つ。$\varphi(b^*ca)=\varphi((c^*b)^*a)$ により対合も保つ。$\xi_\varphi=[1]$ に各 $a$ を作用させると $[a]$ となる。

非単位的なら単位化の延長状態から構成する。$\norm{[1]-[e_\lambda]}^2=1-2\varphi(e_\lambda)+\varphi(e_\lambda^2)\to0$ なので、$[A]$ は単位化から得た空間にも稠密である。一意性は二つの巡回表現の間で $\pi(a)\xi\mapsto\pi'(a)\xi'$ と定める。内積がどちらも $\varphi(b^*a)$ であり、等長写像が稠密像をもつためユニタリへ延びる。さらに左乗法と可換であり、非単位的には近似単位元を使って巡回ベクトルも移す。
:::

:::example[行列と測度での商]
$\varphi(a)=a_{11}$ では $[a]\mapsto ae_1$ が等長であり、$H_\varphi=\C^2$、$\pi_\varphi(a)=a$ である。$C([0,1])$ で $\varphi(f)=\tfrac14 f(0)+\tfrac34 f(1)$ なら、$[f]\mapsto(f(0)/2,\sqrt3 f(1)/2)$ が $\C^2$ へのユニタリである。作用は二点での対角乗算になる。
:::

### GNS 核と左零イデアルの違い

:::proposition[GNS 表現の核]
$\Ker\pi_\varphi$ は $N_\varphi=\{a:\varphi(a^*a)=0\}$ に含まれる最大の閉両側イデアルである。$\varphi$ がトレースなら $N_\varphi=\Ker\pi_\varphi$ である。
:::

:::proof
$\pi_\varphi(a)=0$ はすべての $b$ について $[ab]=0$ と同値であり、巡回ベクトルに作用させれば $a\in N_\varphi$ である。$J\subset N_\varphi$ が両側イデアルなら $a\in J$ に対して $ab\in J$ だから $\pi_\varphi(a)=0$ になる。トレースの場合は
$$
\varphi(b^*a^*ab)
=\varphi(abb^*a^*)
\le\norm b^2\varphi(aa^*)=\norm b^2\varphi(a^*a)
$$
によって右からの積も零長である。
:::

:::example[表現が忠実でも状態に零長元がある]
$M_2$ で $\varphi(a)=a_{11}$ なら $N_\varphi$ は第 1 列が零の行列である。例えば $E_{12}\in N_\varphi$ だが $E_{12}E_{21}=E_{11}$ はそこに属さない。一方、$[a]\mapsto ae_1$ による同一視で $\pi_\varphi(b)$ は $\C^2$ 上の $b$ 自身なので、核は零である。巡回ベクトルを消すだけでは表現上の零作用素にならないことが分かる。
:::

:::example[一列だけを残す表現]
有限行列の状態 $\varphi_D$ の GNS 構成を、商 $[a]\mapsto aD^{1/2}$ として具体化する。この同一視は $\tr((bD^{1/2})^*aD^{1/2})=\varphi_D(b^*a)$ を満たす。

各列は密度の固有値の平方根だけスケールされる。$D$ の零成分に対応する列は消えるので、その列が零長イデアルの自由度になる。

$$
\begin{aligned}
D&=\begin{pmatrix}1 & 0 & 0\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix}\\
a&=\begin{pmatrix}1 & 2 & 0\\0 & 1 & 2\\3 & 0 & 1\end{pmatrix}\\
aD^{1/2}&=\begin{pmatrix}1 & 0 & 0\\0 & 0 & 0\\3 & 0 & 0\end{pmatrix}
\end{aligned}
$$

左掛け作用では $E_{ij}$ が第 $j$ 行を第 $i$ 行へ移し、他の行を消す。その操作は列をスケールする右掛けと交換する。

$$
\begin{aligned}
E_{21}aD^{1/2}&=\begin{pmatrix}0 & 0 & 0\\1 & 0 & 0\\0 & 0 & 0\end{pmatrix}\\
\norm{[a]}^2&=10
\end{aligned}
$$

正の重みの数を $r$ とする。各残った列は三次元で、その列を個別に指定できる。したがって GNS 空間の次元は $3r$ になる。巡回ベクトルは $D^{1/2}$ で、左掛けにより各列が自由に動く。

$$
\begin{aligned}
r&=1\\
\dim H_\varphi&=3\\
\xi_\varphi&=\begin{pmatrix}1 & 0 & 0\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix}
\end{aligned}
$$

零長イデアルが非零でも、少なくとも一列が残っていれば左掛けの表現は忠実である。残った列に任意のベクトルを入れて試せるからである。
:::

:::example[中間の列を捨てる表現]
有限行列の状態 $\varphi_D$ の GNS 構成を、商 $[a]\mapsto aD^{1/2}$ として具体化する。この同一視は $\tr((bD^{1/2})^*aD^{1/2})=\varphi_D(b^*a)$ を満たす。

各列は密度の固有値の平方根だけスケールされる。$D$ の零成分に対応する列は消えるので、その列が零長イデアルの自由度になる。

$$
\begin{aligned}
D&=\begin{pmatrix}\frac{2}{3} & 0 & 0\\0 & 0 & 0\\0 & 0 & \frac{1}{3}\end{pmatrix}\\
a&=\begin{pmatrix}1 & 0 & 2\\3 & 1 & 1\\2 & 4 & 0\end{pmatrix}\\
aD^{1/2}&=\begin{pmatrix}\frac{\sqrt{6}}{3} & 0 & \frac{2 \sqrt{3}}{3}\\\sqrt{6} & 0 & \frac{\sqrt{3}}{3}\\\frac{2 \sqrt{6}}{3} & 0 & 0\end{pmatrix}
\end{aligned}
$$

左掛け作用では $E_{ij}$ が第 $j$ 行を第 $i$ 行へ移し、他の行を消す。その操作は列をスケールする右掛けと交換する。

$$
\begin{aligned}
E_{23}aD^{1/2}&=\begin{pmatrix}0 & 0 & 0\\\frac{2 \sqrt{6}}{3} & 0 & 0\\0 & 0 & 0\end{pmatrix}\\
\norm{[a]}^2&=11
\end{aligned}
$$

正の重みの数を $r$ とする。各残った列は三次元で、その列を個別に指定できる。したがって GNS 空間の次元は $3r$ になる。巡回ベクトルは $D^{1/2}$ で、左掛けにより各列が自由に動く。

$$
\begin{aligned}
r&=2\\
\dim H_\varphi&=6\\
\xi_\varphi&=\begin{pmatrix}\frac{\sqrt{6}}{3} & 0 & 0\\0 & 0 & 0\\0 & 0 & \frac{\sqrt{3}}{3}\end{pmatrix}
\end{aligned}
$$

零長イデアルが非零でも、少なくとも一列が残っていれば左掛けの表現は忠実である。残った列に任意のベクトルを入れて試せるからである。
:::

:::example[二番目の列から巡回する表現]
有限行列の状態 $\varphi_D$ の GNS 構成を、商 $[a]\mapsto aD^{1/2}$ として具体化する。この同一視は $\tr((bD^{1/2})^*aD^{1/2})=\varphi_D(b^*a)$ を満たす。

各列は密度の固有値の平方根だけスケールされる。$D$ の零成分に対応する列は消えるので、その列が零長イデアルの自由度になる。

$$
\begin{aligned}
D&=\begin{pmatrix}0 & 0 & 0\\0 & 1 & 0\\0 & 0 & 0\end{pmatrix}\\
a&=\begin{pmatrix}1 & 0 & 0\\0 & i & 2\\1 & 1 & 1\end{pmatrix}\\
aD^{1/2}&=\begin{pmatrix}0 & 0 & 0\\0 & i & 0\\0 & 1 & 0\end{pmatrix}
\end{aligned}
$$

左掛け作用では $E_{ij}$ が第 $j$ 行を第 $i$ 行へ移し、他の行を消す。その操作は列をスケールする右掛けと交換する。

$$
\begin{aligned}
E_{12}aD^{1/2}&=\begin{pmatrix}0 & i & 0\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix}\\
\norm{[a]}^2&=2
\end{aligned}
$$

正の重みの数を $r$ とする。各残った列は三次元で、その列を個別に指定できる。したがって GNS 空間の次元は $3r$ になる。巡回ベクトルは $D^{1/2}$ で、左掛けにより各列が自由に動く。

$$
\begin{aligned}
r&=1\\
\dim H_\varphi&=3\\
\xi_\varphi&=\begin{pmatrix}0 & 0 & 0\\0 & 1 & 0\\0 & 0 & 0\end{pmatrix}
\end{aligned}
$$

零長イデアルが非零でも、少なくとも一列が残っていれば左掛けの表現は忠実である。残った列に任意のベクトルを入れて試せるからである。
:::

:::example[忠実状態による三重の行列作用]
有限行列の状態 $\varphi_D$ の GNS 構成を、商 $[a]\mapsto aD^{1/2}$ として具体化する。この同一視は $\tr((bD^{1/2})^*aD^{1/2})=\varphi_D(b^*a)$ を満たす。

各列は密度の固有値の平方根だけスケールされる。$D$ の零成分に対応する列は消えるので、その列が零長イデアルの自由度になる。

$$
\begin{aligned}
D&=\begin{pmatrix}\frac{1}{4} & 0 & 0\\0 & \frac{1}{4} & 0\\0 & 0 & \frac{1}{2}\end{pmatrix}\\
a&=\begin{pmatrix}1 & 0 & 2\\1 & 2 & 1\\3 & 1 & 0\end{pmatrix}\\
aD^{1/2}&=\begin{pmatrix}\frac{1}{2} & 0 & \sqrt{2}\\\frac{1}{2} & 1 & \frac{\sqrt{2}}{2}\\\frac{3}{2} & \frac{1}{2} & 0\end{pmatrix}
\end{aligned}
$$

左掛け作用では $E_{ij}$ が第 $j$ 行を第 $i$ 行へ移し、他の行を消す。その操作は列をスケールする右掛けと交換する。

$$
\begin{aligned}
E_{32}aD^{1/2}&=\begin{pmatrix}0 & 0 & 0\\0 & 0 & 0\\\frac{1}{2} & 1 & \frac{\sqrt{2}}{2}\end{pmatrix}\\
\norm{[a]}^2&=\frac{13}{2}
\end{aligned}
$$

正の重みの数を $r$ とする。各残った列は三次元で、その列を個別に指定できる。したがって GNS 空間の次元は $3r$ になる。巡回ベクトルは $D^{1/2}$ で、左掛けにより各列が自由に動く。

$$
\begin{aligned}
r&=3\\
\dim H_\varphi&=9\\
\xi_\varphi&=\begin{pmatrix}\frac{1}{2} & 0 & 0\\0 & \frac{1}{2} & 0\\0 & 0 & \frac{\sqrt{2}}{2}\end{pmatrix}
\end{aligned}
$$

零長イデアルが非零でも、少なくとも一列が残っていれば左掛けの表現は忠実である。残った列に任意のベクトルを入れて試せるからである。
:::

## 純粋状態と既約表現

極点かどうかは、期待値を二つの測定規則に分けられるかで決まる。GNS 空間では、この分解は表現と可換な正作用素による重み付けとして現れる。

:::theorem[純粋性と既約性]{#thm-oa-pure-gns}
状態空間の極点を :term[純粋状態]{en="pure state" yomi="じゅんすいじょうたい"} という。$\varphi$ が純粋であることは、$\pi_\varphi$ に非自明な閉不変部分空間がないことと同値である。そのような表現を :term[既約表現]{en="irreducible representation" yomi="きやくひょうげん"} という。
:::

:::proof
$0\le\psi\le\varphi$ を正値汎関数の順序とする。$B([a],[b])=\psi(b^*a)$ は Cauchy–Schwarz により $\abs{B([a],[b])}\le\norm{[a]}\norm{[b]}$ を満たし、$0\le B(\eta,\eta)\le\norm\eta^2$ である。Riesz 表現で一意な $0\le T\le1$ があり $B(\eta,\zeta)=\inner{T\eta}{\zeta}$ と書ける。$B([ca],[b])=B([a],[c^*b])$ から $T$ は $\pi_\varphi(A)$ と可換である。単位的なら $[1]$、非単位的なら $[e_\lambda]\to\xi_\varphi$ を使うと $\psi(a)=\inner{T\pi_\varphi(a)\xi_\varphi}{\xi_\varphi}$ になる。逆に可換な正縮小作用素 $T$ はこの式で $0\le\psi\le\varphi$ を与える。二次形式に戻せるので対応は単射である。

$\varphi$ が極点という条件は、$0\le\psi\le\varphi$ がすべて $\varphi$ のスカラー倍であることと同値である。$\psi$ と $\varphi-\psi$ が非零なら、それぞれのノルムで割って二状態の凸結合とする。係数の和が $1$ になることは近似単位元で評価して分かる。従って純粋性は可換子の正縮小作用素がスカラーだけであることに等しい。不変部分空間は対合保存性により還元部分空間で、そこへの射影が可換子に入る。逆に自己共役な非スカラーの可換子があれば、その異なる二つのスペクトル点の近くを台とする非零な連続関数 $f,g$ で $fg=0$ とできる。$\overline{\Img f(T)}$ は非零で、$\Img g(T)$ と直交し、かつ不変である。これで同値性が従う。
:::

### 純粋状態を延長する

:::proposition[純粋状態の純粋延長]
同じ単位元を持つ単位的 C*-環 $B\subset A$ に対し、$B$ の純粋状態は $A$ の純粋状態に延長できる。
:::

:::proof
延長全体 $F$ は :ref[lem-oa-state-extension] によって空でない。$F$ は $S(A)$ の中で各 $b\in B$ の値を固定した弱 * 閉集合であり、凸かつコンパクトである。$\psi\in F$ と $\psi=t\psi_1+(1-t)\psi_2$（$0<t<1$）があれば、制限先での純粋性により $\psi_i|_B=\psi|_B$ となり $\psi_i\in F$ である。この性質を持つ凸部分集合を :term[面]{en="face" yomi="めん"} という。Krein–Milman の定理による $F$ の極点は、面の性質から $S(A)$ 全体でも極点であり、求める延長になる。
:::

:::example[対角評価は延長の非対角成分を固定する]
$B\subset M_2$ を対角部分環とし、第 1 成分の評価を延長する。密度行列の対角は $(1,0)$ になる。$\det D\ge0$ から非対角成分 $z$ は $-|z|^2\ge0$ を満たし $z=0$ である。従って延長は $E_{11}$ のベクトル状態に一意である。これに対し等確率の状態を延長すると、対角は $(1/2,1/2)$、非対角成分は $|z|\le1/2$ を満たす円板を動く。境界では行列式が零で階数一となり、環全体で純粋な延長になる。
:::

:::example[中心の密度を二通りに分ける]
対角密度 $D$ を純粋状態の凸結合に分ける。同じ密度でも異なる純粋状態の分解があるか、階数一射影の式で確かめる。

$D$ の対角を $p,q$ と置くと $D=pE_{11}+qE_{22}$ が一つの分解である。各行列単位は単位ベクトル状態を表す。

$$
\begin{aligned}
D&=\begin{pmatrix}\frac{1}{2} & 0\\0 & \frac{1}{2}\end{pmatrix}\\
(p,q)&=\left( \frac{1}{2}, \  \frac{1}{2}\right)
\end{aligned}
$$

次に $v=(\sqrt p,\sqrt q)^{\mathsf T}$、$w=(\sqrt p,-\sqrt q)^{\mathsf T}$ を取る。両方とも長さ $1$ であり、外積の対角は同じ、非対角は反対符号である。

$$
\begin{aligned}
vv^*&=\begin{pmatrix}\frac{1}{2} & \frac{1}{2}\\\frac{1}{2} & \frac{1}{2}\end{pmatrix}\\
ww^*&=\begin{pmatrix}\frac{1}{2} & - \frac{1}{2}\\- \frac{1}{2} & \frac{1}{2}\end{pmatrix}\\
\tfrac12(vv^*+ww^*)&=\begin{pmatrix}\frac{1}{2} & 0\\0 & \frac{1}{2}\end{pmatrix}
\end{aligned}
$$

二つの分解の違いを $E_{12}+E_{21}$ の期待値で読む。純粋状態ごとの値は次のように違っても平均は同じである。

$$
\begin{aligned}
\varphi_v(E_{12}+E_{21})&=1\\
\varphi_w(E_{12}+E_{21})&=-1\\
\varphi_D(E_{12}+E_{21})&=0
\end{aligned}
$$

$p,q$ がともに正なら二つ目の分解は一つ目とは異なる。どちらかが零なら $vv^*=ww^*$ となり、同じ純粋状態に退化する。
:::

:::example[対角の一方に重みが偏る]
対角密度 $D$ を純粋状態の凸結合に分ける。同じ密度でも異なる純粋状態の分解があるか、階数一射影の式で確かめる。

$D$ の対角を $p,q$ と置くと $D=pE_{11}+qE_{22}$ が一つの分解である。各行列単位は単位ベクトル状態を表す。

$$
\begin{aligned}
D&=\begin{pmatrix}\frac{4}{5} & 0\\0 & \frac{1}{5}\end{pmatrix}\\
(p,q)&=\left( \frac{4}{5}, \  \frac{1}{5}\right)
\end{aligned}
$$

次に $v=(\sqrt p,\sqrt q)^{\mathsf T}$、$w=(\sqrt p,-\sqrt q)^{\mathsf T}$ を取る。両方とも長さ $1$ であり、外積の対角は同じ、非対角は反対符号である。

$$
\begin{aligned}
vv^*&=\begin{pmatrix}\frac{4}{5} & \frac{2}{5}\\\frac{2}{5} & \frac{1}{5}\end{pmatrix}\\
ww^*&=\begin{pmatrix}\frac{4}{5} & - \frac{2}{5}\\- \frac{2}{5} & \frac{1}{5}\end{pmatrix}\\
\tfrac12(vv^*+ww^*)&=\begin{pmatrix}\frac{4}{5} & 0\\0 & \frac{1}{5}\end{pmatrix}
\end{aligned}
$$

二つの分解の違いを $E_{12}+E_{21}$ の期待値で読む。純粋状態ごとの値は次のように違っても平均は同じである。

$$
\begin{aligned}
\varphi_v(E_{12}+E_{21})&=\frac{4}{5}\\
\varphi_w(E_{12}+E_{21})&=- \frac{4}{5}\\
\varphi_D(E_{12}+E_{21})&=0
\end{aligned}
$$

$p,q$ がともに正なら二つ目の分解は一つ目とは異なる。どちらかが零なら $vv^*=ww^*$ となり、同じ純粋状態に退化する。
:::

:::example[均等から少しずれた場合]
対角密度 $D$ を純粋状態の凸結合に分ける。同じ密度でも異なる純粋状態の分解があるか、階数一射影の式で確かめる。

$D$ の対角を $p,q$ と置くと $D=pE_{11}+qE_{22}$ が一つの分解である。各行列単位は単位ベクトル状態を表す。

$$
\begin{aligned}
D&=\begin{pmatrix}\frac{3}{5} & 0\\0 & \frac{2}{5}\end{pmatrix}\\
(p,q)&=\left( \frac{3}{5}, \  \frac{2}{5}\right)
\end{aligned}
$$

次に $v=(\sqrt p,\sqrt q)^{\mathsf T}$、$w=(\sqrt p,-\sqrt q)^{\mathsf T}$ を取る。両方とも長さ $1$ であり、外積の対角は同じ、非対角は反対符号である。

$$
\begin{aligned}
vv^*&=\begin{pmatrix}\frac{3}{5} & \frac{\sqrt{6}}{5}\\\frac{\sqrt{6}}{5} & \frac{2}{5}\end{pmatrix}\\
ww^*&=\begin{pmatrix}\frac{3}{5} & - \frac{\sqrt{6}}{5}\\- \frac{\sqrt{6}}{5} & \frac{2}{5}\end{pmatrix}\\
\tfrac12(vv^*+ww^*)&=\begin{pmatrix}\frac{3}{5} & 0\\0 & \frac{2}{5}\end{pmatrix}
\end{aligned}
$$

二つの分解の違いを $E_{12}+E_{21}$ の期待値で読む。純粋状態ごとの値は次のように違っても平均は同じである。

$$
\begin{aligned}
\varphi_v(E_{12}+E_{21})&=\frac{2 \sqrt{6}}{5}\\
\varphi_w(E_{12}+E_{21})&=- \frac{2 \sqrt{6}}{5}\\
\varphi_D(E_{12}+E_{21})&=0
\end{aligned}
$$

$p,q$ がともに正なら二つ目の分解は一つ目とは異なる。どちらかが零なら $vv^*=ww^*$ となり、同じ純粋状態に退化する。
:::

:::example[境界で分解が退化する]
対角密度 $D$ を純粋状態の凸結合に分ける。同じ密度でも異なる純粋状態の分解があるか、階数一射影の式で確かめる。

$D$ の対角を $p,q$ と置くと $D=pE_{11}+qE_{22}$ が一つの分解である。各行列単位は単位ベクトル状態を表す。

$$
\begin{aligned}
D&=\begin{pmatrix}1 & 0\\0 & 0\end{pmatrix}\\
(p,q)&=\left( 1, \  0\right)
\end{aligned}
$$

次に $v=(\sqrt p,\sqrt q)^{\mathsf T}$、$w=(\sqrt p,-\sqrt q)^{\mathsf T}$ を取る。両方とも長さ $1$ であり、外積の対角は同じ、非対角は反対符号である。

$$
\begin{aligned}
vv^*&=\begin{pmatrix}1 & 0\\0 & 0\end{pmatrix}\\
ww^*&=\begin{pmatrix}1 & 0\\0 & 0\end{pmatrix}\\
\tfrac12(vv^*+ww^*)&=\begin{pmatrix}1 & 0\\0 & 0\end{pmatrix}
\end{aligned}
$$

二つの分解の違いを $E_{12}+E_{21}$ の期待値で読む。純粋状態ごとの値は次のように違っても平均は同じである。

$$
\begin{aligned}
\varphi_v(E_{12}+E_{21})&=0\\
\varphi_w(E_{12}+E_{21})&=0\\
\varphi_D(E_{12}+E_{21})&=0
\end{aligned}
$$

$p,q$ がともに正なら二つ目の分解は一つ目とは異なる。どちらかが零なら $vv^*=ww^*$ となり、同じ純粋状態に退化する。
:::

## Gelfand–Naimark の表現定理

:::lemma[状態の延長とノルムの達成]{#lem-oa-state-extension}
単位元を共有する部分環の状態は全環へ状態として延長できる。非零な正元 $h$ には $\varphi(h)=\norm h$ となる状態が存在する。
:::
:::proof
:ref[thm-fa-hahn-banach] でノルムを保って延長すれば $\norm f=f(1)=1$ である。自己共役 $k$ について $\abs{f(1+itk)}^2\le\norm{1+itk}^2\le1+t^2\norm k^2$ である。$t$ の両符号で一次項を調べると $f(k)$ は実数である。$0\le k\le1$ なら $\norm{1-k}\le1$ だから $1-f(k)\le1$、従って $f(k)\ge0$ である。これで延長は正値となる。$C^*(1,h)$ のスペクトル点 $\norm h$ での評価を延長すればノルムが達成される。非単位的な場合は単位化から制限し、$h$ でノルム $1$ を確認できる。
:::

:::theorem[忠実表現]{#thm-oa-faithful-gn}
任意の C*-環 $A$ はある $B(H)$ のノルム閉 *-部分環と等長 *-同型である。
:::
:::proof
全状態の GNS 表現を Hilbert 直和に取り $\pi=\bigoplus_\varphi\pi_\varphi$ とする。各表現は縮小的なので直和は有界で、演算を保つ。$a\ne0$ に対し :ref[lem-oa-state-extension] で $\varphi(a^*a)=\norm a^2$ と取ると $\norm{\pi_\varphi(a)\xi_\varphi}^2=\norm a^2$ である。よって $\pi$ は等長となり、その像は完備で閉じる。零環は零作用素だけで表す。
:::

### 可分環の忠実状態

:::proposition[可分 C*-環の忠実状態]
零でない可分 C*-環は忠実な状態を持つ。
:::

:::proof
正縮小元全体の可算稠密族 $(h_n)$ を取る。各 $h_n$ について $\varphi_n(h_n)=\norm{h_n}$ となる状態を選ぶ。$h_n=0$ なら任意の状態でよい。$\varphi=\sum_{n\ge1}2^{-n}\varphi_n$ は双対ノルムで収束する正値汎関数である。単位的なら単位元の値が $1$ である。非単位的なら近似単位元で、有限部分和の評価をその重みの和に任意に近づけられる。全体の評価は有限部分和以上で、ノルムは $1$ 以下なので、やはりノルムが $1$ である。$a$ が非零の正縮小元なら $\norm{a-h_n}<\norm a/3$ の $n$ を取ると
$$
\varphi_n(a)\ge\norm{h_n}-\norm{a-h_n}>\norm a/3
$$
である。従って $\varphi(a)\ge2^{-n}\varphi_n(a)>0$ である。
:::

:::example[忠実状態が存在しない場合]
$X$ を非可算離散空間とする。$c_0(X)$ の各座標の単位射影を $p_x$ と書く。状態 $\varphi$ の値 $\varphi(p_x)\ge\eps$ となる点は高々 $1/\eps$ 個である。実際、有限個の相異なる座標射影の和はノルム $1$ の射影であり、その値の和が $1$ 以下である。従って $\varphi(p_x)>0$ となる点は可算個しかなく、忠実状態は存在しない。それでも全点の評価を直和にすれば $\ell^2(X)$ 上の忠実表現が得られる。任意の非零関数はどこかの座標で零でないからである。
:::

:::example[一つのイデアルを完全に消す状態]
$A=\bigoplus_jM_{d_j}$ 上で各成分の規格化トレース $\tau_j$ を非負の重みで足した状態を考える。この GNS 表現で見える成分を調べ、別の状態を足したときに忠実性が回復する条件も考える。

各 $a_j$ は指定された固有値の対角行列とする。状態の値は成分内で平均し、成分間の重みを掛けてから足す。

| 成分 | 行列サイズ | $a_j$ の固有値 | 重み | 規格化トレース |
| --- | --- | --- | --- | --- |
| $1$ | $2$ | $\left( 1, \  2\right)$ | $1$ | $\frac{3}{2}$ |
| $2$ | $3$ | $\left( -3, \  1, \  4\right)$ | $0$ | $\frac{2}{3}$ |

正の重みの各成分ではトレースは忠実で、零重みの成分は全体が GNS の零イデアルに入る。正の成分の GNS 空間は $M_{d_j}$ 自身であるため、次元はその二乗の和である。

$$
\begin{aligned}
\varphi(a)&=\frac{3}{2}\\
\dim H_\varphi&=4
\end{aligned}
$$

別の状態の GNS 表現を直和したとき、直和の核は二つの核の共通部分である。各成分に少なくとも一つ正の重みが付けば共通部分は零になる。そこで消したどの成分にも正の重みを置く状態を一つ追加すれば忠実表現を得る。

状態の評価は重みを変えると値が変わるが、ここでの表現の核は重みが零かどうかだけで決まる。
:::

:::example[直和の全成分を検出する状態]
$A=\bigoplus_jM_{d_j}$ 上で各成分の規格化トレース $\tau_j$ を非負の重みで足した状態を考える。この GNS 表現で見える成分を調べ、別の状態を足したときに忠実性が回復する条件も考える。

各 $a_j$ は指定された固有値の対角行列とする。状態の値は成分内で平均し、成分間の重みを掛けてから足す。

| 成分 | 行列サイズ | $a_j$ の固有値 | 重み | 規格化トレース |
| --- | --- | --- | --- | --- |
| $1$ | $2$ | $\left( 1, \  -1\right)$ | $\frac{1}{3}$ | $0$ |
| $2$ | $1$ | $\left( 4,\right)$ | $\frac{2}{3}$ | $4$ |

正の重みの各成分ではトレースは忠実で、零重みの成分は全体が GNS の零イデアルに入る。正の成分の GNS 空間は $M_{d_j}$ 自身であるため、次元はその二乗の和である。

$$
\begin{aligned}
\varphi(a)&=\frac{8}{3}\\
\dim H_\varphi&=5
\end{aligned}
$$

別の状態の GNS 表現を直和したとき、直和の核は二つの核の共通部分である。各成分に少なくとも一つ正の重みが付けば共通部分は零になる。そこで消したどの成分にも正の重みを置く状態を一つ追加すれば忠実表現を得る。

状態の評価は重みを変えると値が変わるが、ここでの表現の核は重みが零かどうかだけで決まる。
:::

:::example[中間の成分だけ消す状態]
$A=\bigoplus_jM_{d_j}$ 上で各成分の規格化トレース $\tau_j$ を非負の重みで足した状態を考える。この GNS 表現で見える成分を調べ、別の状態を足したときに忠実性が回復する条件も考える。

各 $a_j$ は指定された固有値の対角行列とする。状態の値は成分内で平均し、成分間の重みを掛けてから足す。

| 成分 | 行列サイズ | $a_j$ の固有値 | 重み | 規格化トレース |
| --- | --- | --- | --- | --- |
| $1$ | $1$ | $\left( 5,\right)$ | $\frac{1}{3}$ | $5$ |
| $2$ | $2$ | $\left( -2, \  10\right)$ | $0$ | $4$ |
| $3$ | $1$ | $\left( 3,\right)$ | $\frac{2}{3}$ | $3$ |

正の重みの各成分ではトレースは忠実で、零重みの成分は全体が GNS の零イデアルに入る。正の成分の GNS 空間は $M_{d_j}$ 自身であるため、次元はその二乗の和である。

$$
\begin{aligned}
\varphi(a)&=\frac{11}{3}\\
\dim H_\varphi&=2
\end{aligned}
$$

別の状態の GNS 表現を直和したとき、直和の核は二つの核の共通部分である。各成分に少なくとも一つ正の重みが付けば共通部分は零になる。そこで消したどの成分にも正の重みを置く状態を一つ追加すれば忠実表現を得る。

状態の評価は重みを変えると値が変わるが、ここでの表現の核は重みが零かどうかだけで決まる。
:::

:::example[左端の点が見えない可換環]
$A=\bigoplus_jM_{d_j}$ 上で各成分の規格化トレース $\tau_j$ を非負の重みで足した状態を考える。この GNS 表現で見える成分を調べ、別の状態を足したときに忠実性が回復する条件も考える。

各 $a_j$ は指定された固有値の対角行列とする。状態の値は成分内で平均し、成分間の重みを掛けてから足す。

| 成分 | 行列サイズ | $a_j$ の固有値 | 重み | 規格化トレース |
| --- | --- | --- | --- | --- |
| $1$ | $1$ | $\left( 8,\right)$ | $0$ | $8$ |
| $2$ | $1$ | $\left( 3,\right)$ | $\frac{1}{3}$ | $3$ |
| $3$ | $1$ | $\left( -4,\right)$ | $\frac{2}{3}$ | $-4$ |

正の重みの各成分ではトレースは忠実で、零重みの成分は全体が GNS の零イデアルに入る。正の成分の GNS 空間は $M_{d_j}$ 自身であるため、次元はその二乗の和である。

$$
\begin{aligned}
\varphi(a)&=- \frac{5}{3}\\
\dim H_\varphi&=2
\end{aligned}
$$

別の状態の GNS 表現を直和したとき、直和の核は二つの核の共通部分である。各成分に少なくとも一つ正の重みが付けば共通部分は零になる。そこで消したどの成分にも正の重みを置く状態を一つ追加すれば忠実表現を得る。

状態の評価は重みを変えると値が変わるが、ここでの表現の核は重みが零かどうかだけで決まる。
:::

## 状態空間と Krein–Milman の定理

:::proposition[極点による状態と正元の検出]{#prop-oa-state-compactness}
非零単位的 C*-環の状態空間は弱 * コンパクト凸で、純粋状態の凸包の弱 * 閉包である。各非零正元のノルムは純粋状態で達成される。
:::
:::proof
:ref[thm-fa-alaoglu] で双対単位球は弱 * コンパクトである。$\varphi(1)=1$ と $\varphi(a^*a)\ge0$ は弱 * 閉条件なので状態空間はコンパクトな凸集合である。:ref[thm-fa-krein-milman] を適用する。$h\ge0$ に対し $\varphi(h)=\norm h$ を満たす集合は :ref[lem-oa-state-extension] により非空で、連続アフィン関数の最大値集合だから閉面である。このコンパクト凸集合の極点は全状態空間でも極点となる。非単位的な場合はこのノルム達成面を単位化で作る。$h\ne0$ に対しその面の状態は制限してノルム $1$ であり、もし制限に分解があれば一意な状態延長で分解が戻るので純粋性も保たれる。
:::

:::example[非単位的な状態空間]
$c_0$ の座標評価 $\delta_n$ は弱 * で零に収束する。従って状態空間は弱 * コンパクトでない。一方、単位化の状態空間で同じ列の極限は無限遠点の評価である。
:::

:::corollary[純粋状態による順序の判定]{#cor-oa-pure-state-order}
自己共役元 $h\in A$ が正であることは、すべての純粋状態 $\varphi$ に対して $\varphi(h)\ge0$ となることと同値である。
:::

:::proof
正元なら定義からすべての状態で非負である。$h$ が正でなければ、負部分 $h_-$ は非零である。:ref[prop-oa-state-compactness] により $\varphi(h_-)=\norm{h_-}=r>0$ となる純粋状態を選ぶ。単位化で $h_+h_-=0$ と $0\le h_-\le r1$ を使うと、$h_+\le\norm{h_+}(1-h_-/r)$ となる。実際この不等式はスペクトルの各点で、$h_+>0$ なら $h_-=0$、$h_+=0$ なら右辺が非負であることから従う。この状態は右辺を零に送るので $\varphi(h_+)=0$ であり、$\varphi(h)=-r<0$ となる。
:::

:::example[関数環の状態と測度]
局所コンパクト Hausdorff 空間 $X$ に対し、:ref[thm-meas-rmk]{text="Riesz 表現定理"} によれば $C_0(X)$ の正値汎関数は有限正則 Borel 測度と一対一に対応し、ノルムは全質量である。従って状態は確率測度に対応する。点質量 $\delta_x$ は純粋である。実際、その凸分解の両成分は $X\setminus\{x\}$ に質量を持てず、規格化からともに $\delta_x$ である。逆に点質量でない正則確率測度なら、質量が $0$ と $1$ の間の Borel 集合が取れる。すべての Borel 集合の質量が $0$ か $1$ なら、内正則性で質量 $1$ のコンパクト集合を選べる。その中の質量 $1$ の閉集合族は有限交差性を持つので、共通点 $x$ がある。$x$ の任意の開近傍の補集合は質量 $1$ を持てないため、その近傍の質量は $1$ である。よって外正則性から $\mu(\{x\})=1$ となる。質量が中間の集合を取れた場合は、その集合と補集合に測度を制限し、それぞれ規格化すれば非自明な凸分解になる。

測度 $\mu$ に対応する GNS 空間は $L^2(X,\mu)$ である。$[f]\mapsto f$ は $\varphi(g^*f)=\int f\bar g\,\dd\mu$ によって等長で、正則性から $C_c(X)$ が $L^2$ に稠密である。具体的には可測集合の指示関数を、内側のコンパクト集合と外側の開集合で挟み、その間で値を変える連続切断関数を取って $L^2$ 近似する。単関数による近似で稠密性が出る。表現は乗算、巡回単位ベクトルは定数 $1$ である。
:::

:::proposition[有界汎関数の正値汎関数による分解]{#prop-oa-functional-positive-span}
C*-環の自己共役な有界線形汎関数 $f$ は $f=f_+-f_-$、$\norm f=\norm{f_+}+\norm{f_-}$ となる正値汎関数 $f_+,f_-$ を一意に持つ。これを :term[Jordan 分解]{en="Jordan decomposition" yomi="じょるだんぶんかい"} という。一般の有界線形汎関数は四つの正値汎関数の複素線形結合である。
:::

:::proof
自己共役な汎関数 $f$ とは $f(a^*)=\overline{f(a)}$ を満たすものである。$f$ のノルムは自己共役単位球で測れる。任意の $a$ について $\abs{f(a)}=\Re f(za)$ となる $\abs z=1$ を選び、$za$ の自己共役部分を取ると、ノルムは増えず $f$ の値が $\abs{f(a)}$ となる。次に $a=a^*$ を状態空間の実関数 $\widehat a(\varphi)=\varphi(a)$ に送る。スペクトルの最大絶対値を持つ端点で評価する状態を :ref[lem-oa-state-extension] と同じ方法で取れるので $\sup_\varphi\abs{\widehat a(\varphi)}=\norm a$ である。よって $f$ をこの実関数空間の部分空間上の汎関数とみなし、実 Hahn–Banach で全ての連続実関数上にノルムを保って延長する。

:ref[thm-meas-rmk]{text="Riesz 表現定理"} により延長は有限符号付き正則測度で表される。その正負分解で二つの正値汎関数を作り、自己共役部分に制限した後で複素化する。各ノルムは単位元での値であり、その和は全変動ノルム、従って $\norm f$ に等しい。非単位的な場合は $f$ を外部単位化へノルムを保って延長し、その自己共役部分を取る。上の構成で得た二つの正値汎関数を $A$ に制限するとノルム和は $\norm f$ 以下であり、逆の不等式は三角不等式から従う。

一意性のため、$f(a_n)\to\norm f$ となる自己共役縮小元 $a_n\in A$ を選ぶ。$p_n=(1+a_n)/2$ は単位化で正縮小元である。$f=g-h$ が所定のノルム等式を満たすとし、$g,h$ を正値かつノルムを保って単位化へ延長して計算する。すると
$$
2\bigl(g(1-p_n)+h(p_n)\bigr)=\norm f-f(a_n)\longrightarrow0
$$
である。Cauchy–Schwarz と $(1-p_n)^2\le1-p_n$ により $g(p_n x p_n)\to g(x)$ となる。同様に $\abs{h(p_nxp_n)}\le\norm x\,h(p_n^2)\le\norm x\,h(p_n)$ なので $h(p_nxp_n)\to0$ である。従って $g(x)=\lim_n f(p_nxp_n)$ が成り立つ。$x\in A$ なら $p_nxp_n$ も $A$ にあり、この式で $g$ が一意に決まる。$h=g-f$ も一意である。

一般の $f$ は $f^\sharp(a)=\overline{f(a^*)}$ を使って $f=(f+f^\sharp)/2+i(f-f^\sharp)/(2i)$ と分ける。二つの自己共役な汎関数にそれぞれ上の分解を適用すればよい。
:::

:::example[正負の成分を持つ行列汎関数]
$D=\diag(2,-3)$ から $f(a)=\tr(Da)$ を作る。この場合、正値汎関数は $2a_{11}$ と $3a_{22}$ となり、ノルムの和は $5$ である。$\diag(1,-1)$ を入れると値は $5$ であるから、実際に差のノルムも $5$ である。二つの正値成分によって $f(1)=-1$ となるので、$f(1)$ は一般の有界汎関数のノルムを決めない。
:::

### 状態空間の切断と凸分解

:::example[二次行列環の密度行列を幾何的に読む]
三つの自己共役行列 $\sigma_1=E_{12}+E_{21}$、$\sigma_2=-iE_{12}+iE_{21}$、$\sigma_3=E_{11}-E_{22}$ を用いる。自己共役でトレース $1$ の行列は $D=(I+\sum_j r_j\sigma_j)/2$ と一意に書ける。$\sigma_j^2=I$、$\sigma_i\sigma_j+\sigma_j\sigma_i=0$（$i\ne j$）なので $(\sum_jr_j\sigma_j)^2=|r|^2I$ であり、トレース零から固有値は $\pm|r|$ である。従って $D$ の正性は $|r|\le1$ と同値で、純粋性は $|r|=1$ に対応する。例えば $r=(0,0,1/3)$ では $D=\diag(2/3,1/3)$ であり、対角の二つの純粋状態の凸結合として表される。中心 $r=0$ は任意の互いに反対向きの二点の平均になるため、純粋状態への分解は一般に一意でない。
:::

:::example[赤道上の純粋状態]
三つの自己共役行列に対する期待値 $r=(r_1,r_2,r_3)$ が指定された状態を復元し、純粋性を判定する。

$h=\sum_jr_j\sigma_j$ と置き $D=(I+h)/2$ とする。全ての $j$ で $\tr(D\sigma_j)=r_j$ になる。$\tr(\sigma_i\sigma_j)=2\delta_{ij}$ だからである。

$$
\begin{aligned}
r&=\left( \frac{3}{5}, \  \frac{4}{5}, \  0\right)\\
h&=\begin{pmatrix}0 & \frac{3}{5} - \frac{4 i}{5}\\\frac{3}{5} + \frac{4 i}{5} & 0\end{pmatrix}\\
D&=\begin{pmatrix}\frac{1}{2} & \frac{3}{10} - \frac{2 i}{5}\\\frac{3}{10} + \frac{2 i}{5} & \frac{1}{2}\end{pmatrix}
\end{aligned}
$$

$h^2=|r|^2I$ なので、二つの射影 $p=(I+h/|r|)/2$ と $q=I-p$ を使うと $D$ を固有値ごとに分けられる。

$$
\begin{aligned}
|r|&=1\\
p&=\begin{pmatrix}\frac{1}{2} & \frac{3}{10} - \frac{2 i}{5}\\\frac{3}{10} + \frac{2 i}{5} & \frac{1}{2}\end{pmatrix}\\
\sigma(D)&=\left\{0, 1\right\}
\end{aligned}
$$

確率は $(1+|r|)/2$ と $(1-|r|)/2$ である。$\tr(D^2)$ はその二乗和になる。これは $|r|=1$ のときだけ $1$ であり、そのとき $D$ は階数一である。

$$
\begin{aligned}
\tr(D^2)&=\frac{1}{2} + 2 \left(\frac{3}{10} - \frac{2 i}{5}\right) \left(\frac{3}{10} + \frac{2 i}{5}\right)
\end{aligned}
$$

ノルム $|r|$ はどの方向へ偏っているかと独立に純粋性を決める。方向を変えても、同じ長さなら密度行列の固有値は同じになる。
:::

:::example[軸上の混合状態]
三つの自己共役行列に対する期待値 $r=(r_1,r_2,r_3)$ が指定された状態を復元し、純粋性を判定する。

$h=\sum_jr_j\sigma_j$ と置き $D=(I+h)/2$ とする。全ての $j$ で $\tr(D\sigma_j)=r_j$ になる。$\tr(\sigma_i\sigma_j)=2\delta_{ij}$ だからである。

$$
\begin{aligned}
r&=\left( 0, \  0, \  - \frac{1}{2}\right)\\
h&=\begin{pmatrix}- \frac{1}{2} & 0\\0 & \frac{1}{2}\end{pmatrix}\\
D&=\begin{pmatrix}\frac{1}{4} & 0\\0 & \frac{3}{4}\end{pmatrix}
\end{aligned}
$$

$h^2=|r|^2I$ なので、二つの射影 $p=(I+h/|r|)/2$ と $q=I-p$ を使うと $D$ を固有値ごとに分けられる。

$$
\begin{aligned}
|r|&=\frac{1}{2}\\
p&=\begin{pmatrix}0 & 0\\0 & 1\end{pmatrix}\\
\sigma(D)&=\left\{\frac{1}{4}, \frac{3}{4}\right\}
\end{aligned}
$$

確率は $(1+|r|)/2$ と $(1-|r|)/2$ である。$\tr(D^2)$ はその二乗和になる。これは $|r|=1$ のときだけ $1$ であり、そのとき $D$ は階数一である。

$$
\begin{aligned}
\tr(D^2)&=\frac{5}{8}
\end{aligned}
$$

ノルム $|r|$ はどの方向へ偏っているかと独立に純粋性を決める。方向を変えても、同じ長さなら密度行列の固有値は同じになる。
:::

:::example[三方向に成分を持つ純粋状態]
三つの自己共役行列に対する期待値 $r=(r_1,r_2,r_3)$ が指定された状態を復元し、純粋性を判定する。

$h=\sum_jr_j\sigma_j$ と置き $D=(I+h)/2$ とする。全ての $j$ で $\tr(D\sigma_j)=r_j$ になる。$\tr(\sigma_i\sigma_j)=2\delta_{ij}$ だからである。

$$
\begin{aligned}
r&=\left( \frac{1}{3}, \  \frac{2}{3}, \  \frac{2}{3}\right)\\
h&=\begin{pmatrix}\frac{2}{3} & \frac{1}{3} - \frac{2 i}{3}\\\frac{1}{3} + \frac{2 i}{3} & - \frac{2}{3}\end{pmatrix}\\
D&=\begin{pmatrix}\frac{5}{6} & \frac{1}{6} - \frac{i}{3}\\\frac{1}{6} + \frac{i}{3} & \frac{1}{6}\end{pmatrix}
\end{aligned}
$$

$h^2=|r|^2I$ なので、二つの射影 $p=(I+h/|r|)/2$ と $q=I-p$ を使うと $D$ を固有値ごとに分けられる。

$$
\begin{aligned}
|r|&=1\\
p&=\begin{pmatrix}\frac{5}{6} & \frac{1}{6} - \frac{i}{3}\\\frac{1}{6} + \frac{i}{3} & \frac{1}{6}\end{pmatrix}\\
\sigma(D)&=\left\{0, 1\right\}
\end{aligned}
$$

確率は $(1+|r|)/2$ と $(1-|r|)/2$ である。$\tr(D^2)$ はその二乗和になる。これは $|r|=1$ のときだけ $1$ であり、そのとき $D$ は階数一である。

$$
\begin{aligned}
\tr(D^2)&=2 \left(\frac{1}{6} - \frac{i}{3}\right) \left(\frac{1}{6} + \frac{i}{3}\right) + \frac{13}{18}
\end{aligned}
$$

ノルム $|r|$ はどの方向へ偏っているかと独立に純粋性を決める。方向を変えても、同じ長さなら密度行列の固有値は同じになる。
:::

:::example[三方向への弱い偏り]
三つの自己共役行列に対する期待値 $r=(r_1,r_2,r_3)$ が指定された状態を復元し、純粋性を判定する。

$h=\sum_jr_j\sigma_j$ と置き $D=(I+h)/2$ とする。全ての $j$ で $\tr(D\sigma_j)=r_j$ になる。$\tr(\sigma_i\sigma_j)=2\delta_{ij}$ だからである。

$$
\begin{aligned}
r&=\left( \frac{1}{3}, \  \frac{1}{3}, \  \frac{1}{3}\right)\\
h&=\begin{pmatrix}\frac{1}{3} & \frac{1}{3} - \frac{i}{3}\\\frac{1}{3} + \frac{i}{3} & - \frac{1}{3}\end{pmatrix}\\
D&=\begin{pmatrix}\frac{2}{3} & \frac{1}{6} - \frac{i}{6}\\\frac{1}{6} + \frac{i}{6} & \frac{1}{3}\end{pmatrix}
\end{aligned}
$$

$h^2=|r|^2I$ なので、二つの射影 $p=(I+h/|r|)/2$ と $q=I-p$ を使うと $D$ を固有値ごとに分けられる。

$$
\begin{aligned}
|r|&=\frac{\sqrt{3}}{3}\\
p&=\begin{pmatrix}\frac{\sqrt{3}}{6} + \frac{1}{2} & \frac{\sqrt{3} \left(\frac{1}{3} - \frac{i}{3}\right)}{2}\\\frac{\sqrt{3} \left(\frac{1}{3} + \frac{i}{3}\right)}{2} & \frac{1}{2} - \frac{\sqrt{3}}{6}\end{pmatrix}\\
\sigma(D)&=\left\{\frac{1}{2} - \frac{\sqrt{3}}{6}, \frac{\sqrt{3}}{6} + \frac{1}{2}\right\}
\end{aligned}
$$

確率は $(1+|r|)/2$ と $(1-|r|)/2$ である。$\tr(D^2)$ はその二乗和になる。これは $|r|=1$ のときだけ $1$ であり、そのとき $D$ は階数一である。

$$
\begin{aligned}
\tr(D^2)&=2 \left(\frac{1}{6} - \frac{i}{6}\right) \left(\frac{1}{6} + \frac{i}{6}\right) + \frac{5}{9}
\end{aligned}
$$

ノルム $|r|$ はどの方向へ偏っているかと独立に純粋性を決める。方向を変えても、同じ長さなら密度行列の固有値は同じになる。
:::

## 演習問題

::::exercise{level=basic}
$D=\diag(2/5,3/5)$ による状態で $a=\begin{pmatrix}1&2\\0&3\end{pmatrix}$ の $\varphi(a)$ と $\varphi(a^*a)$ を求めよ。
:::solution
対角成分の重み付き和から $\varphi(a)=2/5+9/5=11/5$ である。$a$ の列の長さの二乗は $1,13$ なので $\varphi(a^*a)=2/5+39/5=41/5$ である。
:::
::::

::::exercise{level=basic}
$\varphi(a)=a_{22}$ の $M_2$ での GNS 零空間と Hilbert 空間を求めよ。
:::solution
$\varphi(a^*a)=\norm{ae_2}^2$ なので零空間は第 2 列が零の行列全体である。商の内積を保つ写像 $[a]\mapsto ae_2$ は任意のベクトルを第 2 列に置けるので $\C^2$ への全射である。
:::
::::

::::exercise{level=basic}
$c_0$ で $a=(1,2,0,\ldots)$ のノルムを達成する状態を一つ挙げよ。
:::solution
第 2 座標の評価は正値かつノルム $1$ であり、$a$ の値 $2$ は上限ノルムと一致する。
:::
::::

::::exercise{level=basic}
$\varphi(a)=a_{11}$ に対し $a=E_{21}$、$b=2E_{21}$ で Cauchy–Schwarz の等号を確認せよ。
:::solution
$a^*a=E_{11}$、$b^*b=4E_{11}$、$b^*a=2E_{11}$ なので両辺は $4$ である。GNS 空間で $[b]=2[a]$ となることが等号の理由である。
:::
::::

::::exercise{level=advanced}
単位的 $A$ で $\varphi$ が状態、$u$ がユニタリ、$\varphi(u)=1$ なら $\pi_\varphi(u)\xi_\varphi=\xi_\varphi$ を示せ。
:::solution
$\norm{[u]-[1]}^2=\varphi(2-u-u^*)=2-1-1=0$ なので二つのベクトルは等しい。
:::
::::

::::exercise{level=advanced}
$M_n$ の状態 $\tr(Da)$ が純粋である条件を求めよ。
:::solution
$D$ は正でトレースが $1$ である。階数が二以上なら固有分解から二つの異なる密度行列の凸結合になる。階数一のとき $D=tD_1+(1-t)D_2$、$0<t<1$ とすると、$\Ker D$ 上の二次形式がともに零であり、正平方根から $\Ker D\subset\Ker D_i$ である。$D_i$ の像は同じ一次元に含まれ、トレース $1$ により $D_i=D$ となる。
:::
::::

::::exercise{level=advanced}
状態 $\varphi$ が $A_+\setminus\{0\}$ ですべて正値を取るなら $\pi_\varphi$ が忠実であることを示せ。
:::solution
$\pi_\varphi(a)=0$ なら $\varphi(a^*a)=\norm{\pi_\varphi(a)\xi_\varphi}^2=0$ である。仮定から $a^*a=0$、C* 条件から $a=0$ となる。逆は一般に成り立たず、$M_2$ の第 1 座標ベクトル状態は $\pi_\varphi=\id$ だが $E_{22}$ を零にする。
:::
::::

::::exercise{level=advanced}
$C(X)$ の二点 $x\ne y$ について $(\delta_x+\delta_y)/2$ の GNS 表現が既約でないことを示せ。
:::solution
評価で $[f]\mapsto(f(x),f(y))/\sqrt2$ と移すと内積が保たれる。二点を分離する連続関数によりこれは $\C^2$ へ全射である。乗算は対角作用素なので第 1 座標が非自明な閉不変部分空間である。
:::
::::

::::exercise{level=advanced}
単位的な $A$ の状態 $\varphi$ に対し、$\varphi(p)=1$ を満たす射影 $p$ があるとする。すべての $a$ について $\varphi(a)=\varphi(pap)$ を示せ。
:::solution
$\varphi(1-p)=0$ なので、Cauchy–Schwarz を $(1-p)^*(1-p)$ に使えば $\varphi((1-p)a)=0$ と $\varphi(a(1-p))=0$ が出る。ここで後者は前者を $a^*$ に適用して共役を取っても得られる。$a=pap+(1-p)a+pa(1-p)$ と分けると、後ろ二項の値が零なので主張が出る。
:::
::::

## 確認テスト

::::quiz{answer="0.5" tol="1e-6"}
$M_2$ の規格化トレースで階数一射影の値を求めよ。
:::explanation
射影の固有値は $1,0$ なのでトレースは $1$ である。単位元のトレース $2$ で割って $1/2$ になる。
:::
::::

::::quiz
GNS 構成で零長となる元は一般にどの構造をもつか。

- [x] 左イデアル
- [ ] 必ず右イデアル
- [ ] 必ず非零

:::explanation
$\varphi(a^*c^*ca)\le\norm c^2\varphi(a^*a)$ により左乗法で閉じる。$\varphi(a)=a_{11}$ では $E_{12}$ が零長だが $E_{12}E_{21}=E_{11}$ は零長でない。規格化トレースの GNS 零空間は零だけである。
:::
::::

::::quiz{answer="4" tol="1e-6"}
$M_2$ における $D=I_2/2$ の状態の GNS 空間の複素次元を答えよ。
:::explanation
$\varphi(a^*a)=\tfrac12\sum_{i,j}\abs{a_{ij}}^2$ なので零空間が零である。商は $M_2$ 自身となり次元は $4$ である。
:::
::::

::::quiz
単位的 C*-環の状態が純粋である条件を選べ。

- [x] 状態空間の極点である。
- [x] GNS 表現が既約である。
- [ ] GNS 表現がいつも一次元である。

:::explanation
極点は定義であり、支配された正値汎関数を可換子に対応させると既約性と同値となる。$M_2$ のベクトル状態はその通常の二次元既約表現を与えるので一次元である必要はない。
:::
::::

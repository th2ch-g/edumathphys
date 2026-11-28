---
title: Hilbert 空間上の有界作用素
description: 共役、自己共役性と正値性、平方根と極分解、強圧形式に対する変分問題を扱う。
notation:
  - symbol: 'T^*'
    meaning: T の Hilbert 空間での共役作用素
  - symbol: 'A\ge0'
    meaning: A は正作用素
  - symbol: '\abs T'
    meaning: T* T の正平方根
---

## この章の目標

- 内積から作用素の共役を構成し、そのノルムと値域の関係を証明する。
- 自己共役、ユニタリ、正規、正作用素の条件を区別する。
- 正平方根と極分解をスペクトル定理に頼らず構成する。
- 強圧形式から等式と不等式の変分問題を解く。

## 共役作用素

行列の共役転置は、内積の左側にある作用を右側に移す。一般の Hilbert 空間では基底を選ぶ代わりに内積による汎関数の表現を使う。

:::theorem[共役作用素]{#thm-fa-adjoint}
$T\in\mathcal B(H,K)$ に対し、全ての $x\in H,y\in K$ で
$$
\inner{Tx}{y}=\inner{x}{T^*y}
$$
を満たす一意の $T^*\in\mathcal B(K,H)$ がある。これを :term[共役作用素]{en="adjoint operator" yomi="きょうやくさようそ"} という。共役は $\norm{T^*}=\norm T$、$(T^*)^*=T$、$(ST)^*=T^*S^*$、$(aT+bU)^*=\conj aT^*+\conj bU^*$ を満たす。
:::

:::proof
固定した $y$ に対する $x\mapsto\inner{Tx}{y}$ のノルムは $\norm T\norm y$ 以下である。:ref[thm-fa-riesz-representation] により表現ベクトル $T^*y$ が一意に定まり、そのノルムは上の評価を満たす。$y$ と第 2 引数の共役線形性を比較して $T^*$ の線形性を得る。定義式を共役すれば二重共役の等式が得られ、二重共役にも同じノルム上界を適用すればノルムは等しい。積則は $\inner{STx}{y}=\inner{Tx}{S^*y}=\inner{x}{T^*S^*y}$ で従い、和の法則も定義式の両辺を展開して得る。
:::

:::corollary[核と値域の直交関係]{#cor-fa-adjoint-ranges}
$\Ker T^*=(\Img T)^\perp$、$\overline{\Img T^*}=(\Ker T)^\perp$ であり、$\norm{T^*T}=\norm T^2$ である。
:::

:::proof
$y\perp\Img T$ は全ての $x$ について $\inner{x}{T^*y}=0$ であることと同値である。$T$ と $T^*$ を交換して直交補を取ると第二式を得る。第三式の上界は合成の評価であり、下界は $\norm{Tx}^2=\inner{x}{T^*Tx}\le\norm{T^*T}\norm x^2$ で従う。
:::

:::example[列をずらす作用]
$\ell^2(\N)$ の右シフト $S(x_1,x_2,\ldots)=(0,x_1,x_2,\ldots)$ に対し、内積の和の添字を一つずらせば $S^*(y_1,y_2,\ldots)=(y_2,y_3,\ldots)$ である。従って $S^*S=I$、$SS^*=I-P_{\spn e_1}$ である。
:::

:::example[積分の順序を替える]
$Vf(t)=\int_0^t f(s)\,ds$ を $L^2[0,1]$ で考える。Cauchy–Schwarz により $\abs{Vf(t)}^2\le t\norm f_2^2$ なので有界である。$f(s)\conj{g(t)}$ は単位正方形上で可積分であり、積分順序を交換すると $V^*g(s)=\int_s^1g(t)\,dt$ となる。例えば $V1=t$ だが $V^*1=1-t$ である。
:::

:::example[長方形行列の共役と左逆を区別する]
標準内積を入れた $\C^3$ で、列行列 $V=\begin{pmatrix}1 & 0\\i & 3\\0 & 1\end{pmatrix}$ の列 $v_j$ が張る部分空間を $M$ とする。$x=\begin{pmatrix}2\\1 - i\\1\end{pmatrix}$ を $M$ の成分 $p$ と直交成分 $r$ に分ける。係数で $p=Vc$ と書くと、各列との直交条件は $V^*(x-Vc)=0$ になる。列が一次独立なので $\inner{Gd}{d}=\norm{Vd}^2>0$ が $d\ne0$ に対して成り立ち、$G=V^*V$ は可逆である。

$$
\begin{aligned}
G=V^*V&=\begin{pmatrix}2 & - 3 i\\3 i & 10\end{pmatrix},\\
b=V^*x&=\begin{pmatrix}2 - i \left(1 - i\right)\\4 - 3 i\end{pmatrix},\\
c=G^{-1}b&=\begin{pmatrix}\frac{20}{11} - \frac{10 i \left(1 - i\right)}{11} + \frac{3 i \left(4 - 3 i\right)}{11}\\\frac{8}{11} - \frac{6 i}{11} - \frac{3 i \left(2 - i \left(1 - i\right)\right)}{11}\end{pmatrix},\\
(p,r)&=\bigl(\begin{pmatrix}\frac{20}{11} - \frac{10 i \left(1 - i\right)}{11} + \frac{3 i \left(4 - 3 i\right)}{11}\\\frac{24}{11} - \frac{18 i}{11} - \frac{9 i \left(2 - i \left(1 - i\right)\right)}{11} + i \left(\frac{20}{11} - \frac{10 i \left(1 - i\right)}{11} + \frac{3 i \left(4 - 3 i\right)}{11}\right)\\\frac{8}{11} - \frac{6 i}{11} - \frac{3 i \left(2 - i \left(1 - i\right)\right)}{11}\end{pmatrix},\begin{pmatrix}\frac{2}{11} - \frac{3 i \left(4 - 3 i\right)}{11} + \frac{10 i \left(1 - i\right)}{11}\\- \frac{13}{11} - i \left(\frac{20}{11} - \frac{10 i \left(1 - i\right)}{11} + \frac{3 i \left(4 - 3 i\right)}{11}\right) + \frac{7 i}{11} + \frac{9 i \left(2 - i \left(1 - i\right)\right)}{11}\\\frac{3}{11} + \frac{3 i \left(2 - i \left(1 - i\right)\right)}{11} + \frac{6 i}{11}\end{pmatrix}\bigr).
\end{aligned}
$$

ここで残差が全ての列に直交することは $V^*r=b-Gc=0$ から確認できる。また $\norm x^2=\norm p^2+\norm r^2$ により計算結果を検算できる。

| 成分 | ノルムの平方 | 計算 |
| --- | --- | --- |
| 入力 | $7$ | $x^*x$ |
| 射影 | $\frac{64}{11}$ | $c^*Gc$ |
| 残差 | $\frac{13}{11}$ | $r^*r$ |

同じ部分空間を直交する方向で表すと、連立方程式の代わりに各方向について内積を 1 回ずつ取ればよい。

| $j$ | 直交化後の $q_j$ | 平方ノルム | $\inner{x}{q_j}/\norm{q_j}^2$ |
| --- | --- | --- | --- |
| 1 | $\begin{pmatrix}1\\i\\0\end{pmatrix}$ | $2$ | $\frac{1}{2} - \frac{i}{2}$ |
| 2 | $\begin{pmatrix}\frac{3 i}{2}\\\frac{3}{2}\\1\end{pmatrix}$ | $\frac{11}{2}$ | $\frac{5}{11} - \frac{9 i}{11}$ |

列写像の共役は $V^*y$ であり、各成分は $\inner{y}{v_j}$ である。$V^*V$ が単位行列なら $V^*$ は値域上で逆写像として働く。この例で実際に値域上の逆として働くのは $G^{-1}V^*$ である。$G^{-1}V^*V=I$ であるから左逆になり、$V G^{-1}V^*$ は値域を固定する。
:::

:::example[一つの成分を別の方向へ送る共役]
標準内積を入れた $\C^3$ で $a=(1,i,0)$、$b=(2,1,i)$ とする。$Tu=\inner{u}{b}a$ は、入力を一つの数にしてから方向 $a$ を付ける写像である。行列は $T=ab^*=\begin{pmatrix}2 & 1 & - i\\2 i & i & 1\\0 & 0 & 0\end{pmatrix}$ である。例えば入力 $x=(1,2,- i)$ では $\inner{x}{b}=3$ である。$T^*T$ には $a^*a=\norm a^2$、$T^2$ には $b^*a=c$ が現れる。

$$
\begin{aligned}
Tx&=\begin{pmatrix}3\\3 i\\0\end{pmatrix},\\
T^*T&=\begin{pmatrix}8 & 4 & - 4 i\\4 & 2 & - 2 i\\4 i & 2 i & 2\end{pmatrix},\\
T^2&=\begin{pmatrix}4 + 2 i & 2 + i & 1 - 2 i\\-2 + 4 i & -1 + 2 i & 2 + i\\0 & 0 & 0\end{pmatrix}.
\end{aligned}
$$

$T^2=cT$ の係数 $c=\inner{a}{b}=2 + i$ は、出力に再び作用した際に得られる数である。従って零以外の固有値があれば $c$ だけである。実際、$Tu=\lambda u$ で $\lambda\ne0$ なら $u\in\spn a$ であり、$Ta=ca$ である。一方 $b^\perp$ は全て核に入る。

任意の $u,v$ に対して $\inner{Tu}{v}=\inner{u}{b}\inner{a}{v}=\inner{u}{\inner{v}{a}b}$ である。従って $T^*v=\inner{v}{a}b$ となる。$T^*T$ は $b$ の方向に $\norm a^2\norm b^2=12$ を掛け、$b^\perp$ を消す。一方 $TT^*$ でこの役割を果たす方向は $a$ である。$a,b$ が一次独立なら二つの正作用素の核が異なり、$T$ は正規ではない。
:::

:::proposition[積分核の共役と計算可能な上界]{#prop-fa-kernel-adjoint}
測度空間の積上の $k\in L^2(Y\times X)$ から $(Tf)(y)=\int_X k(y,x)f(x)\,d\mu(x)$ を定めるとき、$\norm T\le\norm k_2$ である。ここでは $X,Y$ は $\sigma$ 有限とする。共役は $(T^*g)(x)=\int_Y\conj{k(y,x)}g(y)\,d\nu(y)$ で与えられる。
:::

:::proof
Fubini によりほとんど全ての $y$ に対し切断 $k(y,\cdot)$ が $L^2$ に属する。Cauchy–Schwarz から $\abs{Tf(y)}^2\le\int_X\abs k^2\,d\mu\,\norm f_2^2$ であり、これを $y$ で積分して上界を得る。$\int\abs{k(y,x)f(x)\overline{g(y)}}\le\norm k_2\norm f_2\norm g_2$ なので内積に Fubini を使える。積分順序を交換し、積分の共役を第 2 引数へ移すと所望の式を得る。例えば $k(y,x)=y+ix$、$X=Y=[0,1]$ なら $T1(y)=y+i/2$、$T^*1(x)=1/2-ix$ であり、核の変数交換と複素共役の両方が現れる。
:::

## 自己共役作用素とユニタリ作用素

:::definition[自己共役とユニタリ]{#def-fa-selfadjoint-unitary}
$A\in\mathcal B(H)$ は $A=A^*$ のとき :term[自己共役]{en="self-adjoint" yomi="じこきょうやく"} である。$U:H\to K$ が有界線形で $U^*U=I_H$、$UU^*=I_K$ なら :term[ユニタリ]{en="unitary" yomi="ゆにたり"} である。
:::

:::proposition[二次形式と内積保存]{#prop-fa-selfadjoint-unitary}
複素 Hilbert 空間では $A$ が自己共役であることと全ての $\inner{Ax}{x}$ が実であることは同値である。$U$ がユニタリであることと全射線形等長写像であることは同値である。
:::

:::proof
自己共役なら $\inner{Ax}{x}=\inner{x}{Ax}=\conj{\inner{Ax}{x}}$ である。逆に $D=A-A^*$ とすれば $\inner{Dx}{x}=0$ である。これを $x+y$ と $x+iy$ に適用すると、$\inner{Dx}{y}+\inner{Dy}{x}=0$ と $-i\inner{Dx}{y}+i\inner{Dy}{x}=0$ が得られるので $D=0$ である。

$U^*U=I_H$ ならノルムを保ち、$UU^*=I_K$ なら全射である。逆に線形等長写像は偏極により内積を保つので $U^*U=I_H$ である。さらに全射なら $UU^*Ux=Ux$ で全ての $K$ の元が尽くされる。
:::

:::example[反射と回転]
$\R^2$ の $A(x,y)=(y,x)$ は自己共役でユニタリでもある。これに対して $U(x,y)=(-y,x)$ は直交変換であるが、$U^*=-U$ なので自己共役でない。実数体では $\inner{Ux}{x}$ が常に実でも自己共役性の判定にはならない。
:::

:::example[四つの位相をもつユニタリ作用素]
$\C^4$ の標準座標で $D=\diag(1,-1,i,- i)$、$x=(1,2,i,1 - i)$ とする。作用は各成分を $\lambda_j$ 倍するので、計算を座標ごとに進められる。

| $j$ | $\lambda_j$ | $x_j$ | $(Dx)_j$ | $\abs{(Dx)_j}^2$ |
| --- | --- | --- | --- | --- |
| 1 | $1$ | $1$ | $1$ | $1$ |
| 2 | $-1$ | $2$ | $-2$ | $4$ |
| 3 | $i$ | $i$ | $-1$ | $1$ |
| 4 | $- i$ | $1 - i$ | $- i \left(1 - i\right)$ | $2$ |

これを足すと $\norm{Dx}^2=8$ である。一般の $x$ に対して $\norm{Dx}^2\le 1\norm x^2$ であり、係数の絶対値が最大の座標を一つ選べば等号を得る。そのため $\norm D=1$ である。複数の座標で最大値を取る場合、最大値に対応する部分空間のどの単位ベクトルでも同じノルムに達する。

共役は $D^*=\begin{pmatrix}1 & 0 & 0 & 0\\0 & -1 & 0 & 0\\0 & 0 & - i & 0\\0 & 0 & 0 & i\end{pmatrix}$ である。実際、$\inner{Dx}{v}=\sum_j \lambda_jx_j\conj{v_j}$ は $\inner{x}{D^*v}$ と等しい。$D^*D=DD^*=\diag(1,1,1,1)$ なので $D$ は正規である。$\lambda_j$ が全て実なら自己共役であり、その上で全てが非負なら正である。また全てが絶対値 $1$ なら $D^*D=DD^*=I$ なのでユニタリになる。この例では定義の行列等式を成分ごとの条件で判定できる。

この例の極分解は後の「正作用素の平方根と極分解」で、一般の極分解を構成した後に確認する。
:::

:::proposition[自己共役部分と反自己共役部分]{#prop-fa-cartesian-operator}
複素 Hilbert 空間上の各有界作用素 $T$ は、自己共役な $A,B$ を用いて一意に $T=A+iB$ と書ける。$A=(T+T^*)/2$、$B=(T-T^*)/(2i)$ で、$\norm A,\norm B\le\norm T$ である。
:::

:::proof
共役の計算で $A=A^*, B=B^*$ が分かる。二つを合わせると $T$ に戻る。任意の自己共役な表示からは $T^*=A-iB$ と書けるので、和と差によって $A,B$ が決まる。ノルム上界は三角不等式と $\norm{T^*}=\norm T$ による。例えば $T=\left(\begin{smallmatrix}0&2\\0&0\end{smallmatrix}\right)$ では $A=\left(\begin{smallmatrix}0&1\\1&0\end{smallmatrix}\right)$、$B=\left(\begin{smallmatrix}0&-i\\i&0\end{smallmatrix}\right)$ であり、各部分はノルム $1$ だが $T$ はノルム $2$ である。二次形式を取り、その実部または虚部を抽出すると $A,B$ の二次形式になる。
:::

:::example[片側シフトと両側シフトで全射性を調べる]
$\ell^2(\Z)$ で $(Ux)_j=x_{j-1}$ とすると、添字集合が整数なので、和の変数を $k=j-1$ に変えるだけで $\norm{Ux}=\norm x$ を得る。任意の $y$ に対し $x_j=y_{j+1}$ は同じノルムを持ち $Ux=y$ を満たす。従って $U^*y=(y_{j+1})_j=U^{-1}y$ である。$x=e_{-1}+2e_0+3i e_1$ の場合は $Ux=e_0+2e_1+3i e_2$ であり、平方ノルムはどちらも $14$ になる。

$\ell^2(\N)$ の右シフトでは、同じ計算でノルム保存を示せる一方、像は第 1 座標が零である。従って逆方向のシフトは左逆であるが、右から合成すると第 1 成分を消す射影になる。これは $U^*U=I$ と $UU^*=I$ を別々に確認する理由である。また $U+U^*$ は両側シフトでも片側シフトでも自己共役である。両側の場合に上の $x$ へ作用させると $e_{-2}+2e_{-1}+(1+3i)e_0+2e_1+3i e_2$ になる。
:::

:::proposition[ユニタリな共役による作用素の移し替え]{#prop-fa-unitary-equivalence}
$U:H\to K$ がユニタリなら、写像 $T\mapsto UTU^*$ は作用素ノルム、積、共役を保ち、$\mathcal B(H)$ から $\mathcal B(K)$ への全射である。自己共役性、正規性、正値性も互いに移る。
:::

:::proof
$k=Ux$ と置けば $\norm k=\norm x$、$\norm{UTU^*k}=\norm{Tx}$ であり、単位球の全てを走るのでノルムが一致する。$U^*U=I$ により積を保ち、共役により $(UTU^*)^*=UT^*U^*$ である。逆は $S\mapsto U^*SU$ である。正値性は $\inner{UTU^*k}{k}=\inner{Tx}{x}$ から保たれ、自己共役性と正規性は積と共役の等式だから保たれる。
:::

:::example[非等長な座標変換では直交性が保存されない例]
$P=\begin{pmatrix}1 & 1 & 0\\0 & 1 & 2\\0 & 0 & 1\end{pmatrix}$ は可逆であり、$\Lambda=\diag(1,3,-1)$ を取って $A=P\Lambda P^{-1}=\begin{pmatrix}1 & 2 & -4\\0 & 3 & -8\\0 & 0 & -1\end{pmatrix}$ とする。行列 $P$ の列はそれぞれの固有値の固有ベクトルである。入力 $x=\begin{pmatrix}1\\2\\3\end{pmatrix}$ をこれらに分解するには $Pc=x$ を解いて $c=P^{-1}x$ を得る。標準内積で各列が正規直交しなければ、各列との内積だけで $c$ は求まらない。

| $j$ | 固有値 | 固有ベクトル | 展開係数 | 入力のその成分 |
| --- | --- | --- | --- | --- |
| 1 | $1$ | $\begin{pmatrix}1\\0\\0\end{pmatrix}$ | $5$ | $\begin{pmatrix}5\\0\\0\end{pmatrix}$ |
| 2 | $3$ | $\begin{pmatrix}1\\1\\0\end{pmatrix}$ | $-4$ | $\begin{pmatrix}-4\\-4\\0\end{pmatrix}$ |
| 3 | $-1$ | $\begin{pmatrix}0\\2\\1\end{pmatrix}$ | $3$ | $\begin{pmatrix}0\\6\\3\end{pmatrix}$ |

この表で入力の各成分を足せば $x$ に戻る。各成分に対応する固有値を掛けて足すと $Ax=\begin{pmatrix}-7\\-18\\-3\end{pmatrix}$ である。従って固有座標の操作と標準座標での行列積が一致する。$\Lambda$ が対角であることは作用の仕組みを簡単にするが、$P$ によるノルムの変化は残る。

最初の固有値 $\lambda_1=1$ だけを残す座標射影は $P_1=P[:,1](P^{-1})[1,:]$ と書ける。ここで添字 $1$ は第 1列と第 1行の意味である。具体的には
$$
P_1=\begin{pmatrix}1 & -1 & 2\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix},\qquad P_1x=\begin{pmatrix}5\\0\\0\end{pmatrix},\qquad P_1-P_1^*=\begin{pmatrix}0 & -1 & 2\\1 & 0 & 0\\-2 & 0 & 0\end{pmatrix}.
$$
$P_1^2=P_1$ であり $AP_1=P_1A=\lambda_1P_1$ である。実際、変換後の対角座標では最初の成分のみを残すだけである。$P_1^*$ と一致するかどうかは、標準内積で像と核が直交するかを調べて決まる。

$\lambda_1$ が他と異なる場合、補間多項式 $p_1(z)=\prod_{j>1}(z-\lambda_j)/(\lambda_1-\lambda_j)$ でも $p_1(A)=P_1$ となる。これは重複のない固有値の成分に $1$、他に $0$ を指定した関数を作用させたものになっている。有限次元でも斜交射影を直交射影と混同すると $P_1 x$ の値が変わりうる。
:::

## 正規作用素と正作用素

:::definition[正規と正]{#def-fa-normal-positive}
$T^*T=TT^*$ を満たす $T\in\mathcal B(H)$ を :term[正規作用素]{en="normal operator" yomi="せいきさようそ"} という。自己共役な $A$ が全ての $x$ に対して $\inner{Ax}{x}\ge0$ を満たすとき :term[正作用素]{en="positive operator" yomi="せいさようそ"} といい $A\ge0$ と書く。$A\le B$ は $B-A\ge0$ を意味する。
:::

:::proposition[正規性とノルムの計算]{#prop-fa-positive-order}
実または複素 Hilbert 空間で、$T$ の正規性は $\norm{Tx}=\norm{T^*x}$ が全ての $x$ に対して成り立つことと同値である。自己共役 $A$ について $\norm A=\sup_{\norm x=1}\abs{\inner{Ax}{x}}$ が成り立つ。
:::

:::proof
平方ノルムの差は $\inner{(T^*T-TT^*)x}{x}$ である。自己共役作用素 $D=T^*T-TT^*$ の二次形式が零なら、$x+y$ と $x-y$ での差から $\Re\inner{Dx}{y}=0$ を得る。複素の場合は $y$ を $iy$ に替えて虚部も零と分かり、実の場合は実部だけで $D=0$ となる。ノルム表示で右辺を $m$ とおく。単位ベクトル $x,y$ に対して位相を選び $\inner{Ax}{y}$ を非負実数としてよい。実空間では符号を選ぶ。自己共役性から
$$
4\abs{\inner{Ax}{y}}
=\inner{A(x+y)}{x+y}-\inner{A(x-y)}{x-y}
\le m(\norm{x+y}^2+\norm{x-y}^2)=4m.
$$
ここで $y=Ax/\norm{Ax}$ を使えば $\norm A\le m$ であり、逆向きは内積の評価である。零空間では両辺を零とする。
:::

:::example[正規性と正値性の違い]
$\C^2$ の $\diag(1,i)$ は正規でユニタリだが正ではない。$\diag(2,3)$ は正で $\norm A=3$ である。$T^*T$ は常に正である一方、右シフトでは $T^*T\ne TT^*$ なので、任意の $T$ 自身が正規とは限らない。
:::

:::example[二次形式で上下の作用素順序を調べる]
$\C^4$ の標準座標で $D=\diag(0,2,7,7)$、$x=(1,2,1 + i,1)$ とする。作用は各成分を $\lambda_j$ 倍するので、計算を座標ごとに進められる。

| $j$ | $\lambda_j$ | $x_j$ | $(Dx)_j$ | $\abs{(Dx)_j}^2$ |
| --- | --- | --- | --- | --- |
| 1 | $0$ | $1$ | $0$ | $0$ |
| 2 | $2$ | $2$ | $4$ | $16$ |
| 3 | $7$ | $1 + i$ | $7 + 7 i$ | $98$ |
| 4 | $7$ | $1$ | $7$ | $49$ |

これを足すと $\norm{Dx}^2=163$ である。一般の $x$ に対して $\norm{Dx}^2\le 49\norm x^2$ であり、係数の絶対値が最大の座標を一つ選べば等号を得る。そのため $\norm D=7$ である。複数の座標で最大値を取る場合、最大値に対応する部分空間のどの単位ベクトルでも同じノルムに達する。

共役は $D^*=\begin{pmatrix}0 & 0 & 0 & 0\\0 & 2 & 0 & 0\\0 & 0 & 7 & 0\\0 & 0 & 0 & 7\end{pmatrix}$ である。実際、$\inner{Dx}{v}=\sum_j \lambda_jx_j\conj{v_j}$ は $\inner{x}{D^*v}$ と等しい。$D^*D=DD^*=\diag(0,4,49,49)$ なので $D$ は正規である。$\lambda_j$ が全て実なら自己共役であり、その上で全てが非負なら正である。また全てが絶対値 $1$ なら $D^*D=DD^*=I$ なのでユニタリになる。この例では定義の行列等式を成分ごとの条件で判定できる。

極分解を求めるときは、絶対値 $H=\abs D=\begin{pmatrix}0 & 0 & 0 & 0\\0 & 2 & 0 & 0\\0 & 0 & 7 & 0\\0 & 0 & 0 & 7\end{pmatrix}$ と $U=\begin{pmatrix}0 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{pmatrix}$ を取る。$H\ge0$、$H^2=D^*D$ である。成分ごとに $UH=D$ となり、$U^*U=UU^*=\begin{pmatrix}0 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{pmatrix}$ は $D$ によって消されない座標への射影である。$U$ はその部分空間上で各座標の位相だけを変える。例えば零の成分に任意の絶対値 $1$ の値を入れても積 $UH$ は変わらないが、$U$ を核上で零にする条件によって極分解の $U$ が定まる。
:::

:::proposition[正作用素による半内積と順序]{#prop-fa-positive-form-cauchy-schwarz}
$A\ge0$ なら $\abs{\inner{Ax}{y}}^2\le\inner{Ax}{x}\inner{Ay}{y}$ である。$0\le A\le B$ なら $\norm A\le\norm B$ であり、任意の有界 $T$ に対し $T^*AT\le T^*BT$ である。
:::

:::proof
$a=\inner{Ay}{y}>0$ なら $0\le\inner{A(x-\lambda y)}{x-\lambda y}$ に $\lambda=\inner{Ax}{y}/a$ を代入する。$a=0$ で $\inner{Ax}{y}\ne0$ とすると、$\lambda$ の位相と絶対値を選んで一次項を負で任意に大きくできるため非負性に反する。従って $a=0$ でも不等式は成り立つ。ノルムの大小は :ref[prop-fa-positive-order] で単位ベクトル上の二次形式を比較すれば従う。最後に $\inner{T^*(B-A)Tx}{x}=\inner{(B-A)Tx}{Tx}\ge0$ である。例えば $A=\diag(1,0)$、$B=\diag(2,1)$、$T(x,y)=(x+y,y)$ なら差の二次形式は $\abs{x+y}^2+\abs y^2$ となる。
:::

:::example[正作用素の積が正になる条件]
$A=\left(\begin{smallmatrix}2&0\\0&1\end{smallmatrix}\right)$、$B=\left(\begin{smallmatrix}2&1\\1&2\end{smallmatrix}\right)$ とする。$\inner{Bx}{x}=\abs{x_1+x_2}^2+\abs{x_1}^2+\abs{x_2}^2$ なので $B\ge I$ であり、$A\ge I$ である。積は $AB=\left(\begin{smallmatrix}4&2\\1&2\end{smallmatrix}\right)$ なので自己共役性を失う。実際に $x=(1,i)$ では $\inner{ABx}{x}=6+i$ である。正作用素の積に正値性を求めるときには積の共役の順序を調べる必要がある。

一方 $ABA=\left(\begin{smallmatrix}8&2\\2&2\end{smallmatrix}\right)$ では、$\inner{ABA x}{x}=\inner{BAx}{Ax}\ge\norm{Ax}^2\ge\norm x^2$ である。従って挟む形では可換性を仮定せずに順序を保てる。もし $A,B$ が可換なら、次節で構成する $A^{1/2}$ も $B$ と可換になり、$AB=A^{1/2}BA^{1/2}\ge0$ である。必要性は $(AB)^*=BA$ と自己共役性から従う。
:::

:::example[作用素と共役の核から正規性を判定する]
標準内積を入れた $\C^3$ で $a=(2,1,0)$、$b=(1,2,3)$ とする。$Tu=\inner{u}{b}a$ は、入力を一つの数にしてから方向 $a$ を付ける写像である。行列は $T=ab^*=\begin{pmatrix}2 & 4 & 6\\1 & 2 & 3\\0 & 0 & 0\end{pmatrix}$ である。例えば入力 $x=(2,-1,1)$ では $\inner{x}{b}=3$ である。$T^*T$ には $a^*a=\norm a^2$、$T^2$ には $b^*a=c$ が現れる。

$$
\begin{aligned}
Tx&=\begin{pmatrix}6\\3\\0\end{pmatrix},\\
T^*T&=\begin{pmatrix}5 & 10 & 15\\10 & 20 & 30\\15 & 30 & 45\end{pmatrix},\\
T^2&=\begin{pmatrix}8 & 16 & 24\\4 & 8 & 12\\0 & 0 & 0\end{pmatrix}.
\end{aligned}
$$

$T^2=cT$ の係数 $c=\inner{a}{b}=4$ は、出力に再び作用した際に得られる数である。従って零以外の固有値があれば $c$ だけである。実際、$Tu=\lambda u$ で $\lambda\ne0$ なら $u\in\spn a$ であり、$Ta=ca$ である。一方 $b^\perp$ は全て核に入る。

任意の $u,v$ に対して $\inner{Tu}{v}=\inner{u}{b}\inner{a}{v}=\inner{u}{\inner{v}{a}b}$ である。従って $T^*v=\inner{v}{a}b$ となる。$T^*T$ は $b$ の方向に $\norm a^2\norm b^2=70$ を掛け、$b^\perp$ を消す。一方 $TT^*$ でこの役割を果たす方向は $a$ である。$a,b$ が一次独立なら二つの正作用素の核が異なり、$T$ は正規ではない。
:::

:::corollary[正規性と実部・虚部の可換性]{#cor-fa-cartesian-normal}
$T=A+iB$、$A,B$ は自己共役とする。$T$ が正規であることと $AB=BA$ は同値である。正規なら $\Ker T=\Ker T^*$ であり、この核の直交補は両作用素で不変である。
:::

:::proof
$T^*T=(A-iB)(A+iB)=A^2+B^2+i(AB-BA)$、$TT^*=A^2+B^2-i(AB-BA)$ である。差が零であることと可換性が同値になる。正規性なら $\norm{Tx}=\norm{T^*x}$ なので核は一致する。$x$ が核に直交し $y$ が核に属すれば $\inner{Tx}{y}=\inner{x}{T^*y}=0$、$\inner{T^*x}{y}=\inner{x}{Ty}=0$ である。従って直交補に制限しても両作用素が閉じる。例えば $A=\diag(1,-1)$ と座標交換 $B$ は可換でないため $A+iB$ は正規でない。
:::

## 直交射影

:::proposition[直交射影の特徴付け]{#prop-fa-orthogonal-projection}
$P\in\mathcal B(H)$ が閉部分空間への直交射影であることと $P^2=P=P^*$ であることは同値である。直交射影 $P,Q$ について $PQ=QP$ なら $PQ$ は二つの値域の共通部分への直交射影である。また $PQ$ が直交射影であれば $P,Q$ は可換である。
:::

:::proof
直交分解から $P^2=P$ と $\inner{Px}{y}=\inner{Px}{Py}=\inner{x}{Py}$ を得る。逆に $M=\Img P=\Ker(I-P)$ は閉じており、$m=Py$ に対して $\inner{(I-P)x}{m}=\inner{P(I-P)x}{y}=0$ なので $P=P_M$ である。可換なら $(PQ)^2=PQ$ と $(PQ)^*=PQ$ が成り立つ。逆に $PQ$ が直交射影なら $PQ=(PQ)^*=QP$ である。$PQx$ は両射影で不変であり、両射影で不変な元は $PQ$ でも不変であるから、値域は共通部分である。
:::

:::example[斜交射影]
$\R^2$ で $P(x,y)=(x+2y,0)$ は $P^2=P$ だが、$P^*(x,y)=(x,2x)$ である。$e_1$ と $(1,2)/\sqrt5$ への直交射影の積を計算すると
$$
P_{e_1}P_{(1,2)}=\frac15\begin{pmatrix}1&2\\0&0\end{pmatrix}.
$$
二つの直線の共通部分は零である。積の値域が共通部分になる議論には可換性が必要である。
:::

:::example[射影の行列を正規直交基底なしで計算する]
標準内積を入れた $\K^3$ で、列行列 $V=\begin{pmatrix}1 & 1\\0 & 2\\2 & 1\end{pmatrix}$ の列 $v_j$ が張る部分空間を $M$ とする。$x=\begin{pmatrix}2\\-1\\0\end{pmatrix}$ を $M$ の成分 $p$ と直交成分 $r$ に分ける。係数で $p=Vc$ と書くと、各列との直交条件は $V^*(x-Vc)=0$ になる。列が一次独立なので $\inner{Gd}{d}=\norm{Vd}^2>0$ が $d\ne0$ に対して成り立ち、$G=V^*V$ は可逆である。

$$
\begin{aligned}
G=V^*V&=\begin{pmatrix}5 & 3\\3 & 6\end{pmatrix},\\
b=V^*x&=\begin{pmatrix}2\\0\end{pmatrix},\\
c=G^{-1}b&=\begin{pmatrix}\frac{4}{7}\\- \frac{2}{7}\end{pmatrix},\\
(p,r)&=\bigl(\begin{pmatrix}\frac{2}{7}\\- \frac{4}{7}\\\frac{6}{7}\end{pmatrix},\begin{pmatrix}\frac{12}{7}\\- \frac{3}{7}\\- \frac{6}{7}\end{pmatrix}\bigr).
\end{aligned}
$$

ここで残差が全ての列に直交することは $V^*r=b-Gc=0$ から確認できる。また $\norm x^2=\norm p^2+\norm r^2$ により計算結果を検算できる。

| 成分 | ノルムの平方 | 計算 |
| --- | --- | --- |
| 入力 | $5$ | $x^*x$ |
| 射影 | $\frac{8}{7}$ | $c^*Gc$ |
| 残差 | $\frac{27}{7}$ | $r^*r$ |

同じ部分空間を直交する方向で表すと、連立方程式の代わりに各方向について内積を 1 回ずつ取ればよい。

| $j$ | 直交化後の $q_j$ | 平方ノルム | $\inner{x}{q_j}/\norm{q_j}^2$ |
| --- | --- | --- | --- |
| 1 | $\begin{pmatrix}1\\0\\2\end{pmatrix}$ | $5$ | $\frac{2}{5}$ |
| 2 | $\begin{pmatrix}\frac{2}{5}\\2\\- \frac{1}{5}\end{pmatrix}$ | $\frac{21}{5}$ | $- \frac{2}{7}$ |

$P=V(V^*V)^{-1}V^*=\begin{pmatrix}\frac{5}{21} & \frac{4}{21} & \frac{8}{21}\\\frac{4}{21} & \frac{20}{21} & - \frac{2}{21}\\\frac{8}{21} & - \frac{2}{21} & \frac{17}{21}\end{pmatrix}$ が行列表示である。$V^*V=G$ を中央に挟めば $P^2=P$ であり、$G=G^*$ を使うと $P^*=P$ である。$V^*r=0$ は $\Ker P$ への所属を直接確認している。$P$ と $I-P$ は和が $I$、積が零である。従って $P$ のノルムは $1$ であり、値域上の非零ベクトルによって達成する。
:::

:::proposition[閉部分空間の包含を射影で判定する]{#prop-fa-projection-order}
直交射影 $P,Q$ について、$\Img P\subset\Img Q$、$QP=P$、$P\le Q$ は同値である。このとき $Q-P$ は $\Img Q\cap(\Img P)^\perp$ への直交射影である。
:::

:::proof
包含は $Q$ が $P$ の像を固定することと同値である。$QP=P$ のとき共役を取って $PQ=P$ を得る。従って $Q-P$ は自己共役で、その平方は $Q-2P+P=Q-P$ になる。直交射影は $\inner{(Q-P)x}{x}=\norm{(Q-P)x}^2$ を満たすので $P\le Q$ である。

逆に $x=Px$ に対して $\norm x^2\le\inner{Qx}{x}=\norm{Qx}^2\le\norm x^2$ であり、直交分解の等号から $Qx=x$ となる。像に属するベクトルは $Q$ で固定され $P$ で消される。逆に両条件を満たす $x$ では $(Q-P)x=x$ なので像の式が従う。
:::

:::example[二本の直線への射影の積を反復する]
$\R^2$ で $P$ を $e_1$ 方向、$Q$ を $v=(3,4)/5$ 方向への直交射影とする。$P Qx=\inner{x}{v}(3/5)e_1$ なので $\norm{PQ}=3/5$ である。二度目以降は $e_1$ に対して $PQe_1=(9/25)e_1$ となるため
$$
(PQ)^n=(9/25)^{n-1}PQ\qquad(n\ge1).
$$
例えば $x=v$ の一回目の像は $(3/5)e_1$、二回目は $(27/125)e_1$ である。従って積は射影とは異なるが、そのべきは二直線の共通部分 $\{0\}$ への射影へノルムで近づく。

$P,Q$ を $\R^3$ へ拡張し、第三の方向 $e_3$ を両方の像に加えると、積は $e_3$ を固定する。直交補 $\spn(e_1,e_2)$ 上では前の計算がそのまま成り立つ。従って $(PQ)^n-P_{e_3}$ のノルムは $(3/5)(9/25)^{n-1}$ である。二つの作用素が可換でなくても、共通部分を固定することと、その直交補で縮むことを別に調べられる。
:::

:::proposition[斜交射影のノルムと補空間の角度]{#prop-fa-oblique-projection-norm}
Hilbert 空間が閉部分空間の代数的直和 $H=M\oplus N$ で、$M\ne\{0\}$ とする。成分射影 $P(m+n)=m$ が有界なら
$$
\norm P=\frac{1}{\inf_{\norm m=1,\,m\in M}\operatorname{dist}(m,N)}
$$
であり、$\norm P=1$ となるのは $N=M^\perp$ の場合に限る。また分母が正なら $P$ は有界である。
:::

:::proof
任意の $m\ne0$ に対して $\norm{P(m+n)}/\norm{m+n}=\norm m/\norm{m+n}$ の $n$ に関する上限は $\norm m/\operatorname{dist}(m,N)$ である。斉次性で $m$ を単位にし、さらに上限を取れば式が従う。有界性から分母は正となり、逆に分母が正なら全ての $m,n$ で比が一様に抑えられる。$\norm P=1$ なら任意の $m\in M$、$n\in N$ について $\norm m\le\norm{m+tn}$ が全 $t\in\K$ で成り立つ。平方を展開して実数と純虚数の方向を使うと $\inner{m}{n}=0$ であり、直和から $N=M^\perp$ となる。逆は Pythagoras の等式で従う。
:::

閉じた補空間から成分射影の有界性を導く一般の議論は、[Baire の定理の帰結](/math/functional-analysis/baire-theorems/) の「補空間と有界射影」で扱う。

:::example[直線の角度から射影の誤差増幅を求める]
$\R^2$ で $M=\spn e_1$、$N=\spn(3,4)$ とする。$e_1$ の $N$ への直交射影は $3(3,4)/25$ なので、距離の平方は $1-9/25=16/25$ である。従って $M$ 成分を取り出す斜交射影のノルムは $5/4$ である。成分の式で確認すると $P(x,y)=(x-3y/4,0)$ となり、Cauchy–Schwarz で同じ値が上界になる。$(4,-3)/5$ は単位ベクトルで、$P(4,-3)/5=(5/4,0)$ となり等号に達する。

$N$ を $\spn(1,\eps)$ に替えると、距離は $\abs\eps/\sqrt{1+\eps^2}$ である。$\eps\to0$ で二方向が近づくと、射影のノルムは $\sqrt{1+\eps^2}/\abs\eps$ に発散する。和の分解自体が一意でも、その分解を計算したときに入力の誤差が増幅される場合がある。
:::

## 正作用素の平方根と極分解

:::theorem[正平方根]{#thm-fa-positive-root}
任意の有界正作用素 $A$ に、一意な正作用素 $R$ が存在して $R^2=A$ となる。$A$ と可換な有界作用素は $R$ と可換である。
:::

:::proof
$A=0$ なら $\norm{Rx}^2=\inner{R^2x}{x}=0$ で一意である。$A\ne0$ なら $B=I-A/\norm A$ とおく。$0\le B\le I$ で、:ref[prop-fa-positive-order] より $\norm B\le1$ である。二項級数を $1-\sqrt{1-t}=\sum_{n\ge1}c_n t^n$ と書く。係数は $c_1=1/2$、$c_{n+1}=c_n(2n-1)/(2n+2)>0$ である。$0\le t<1$ の等式で $t\uparrow1$ として単調収束を使うと $\sum c_n=1$ となる。従って $R=\sqrt{\norm A}(I-\sum c_n B^n)$ はノルムで定まる。$B^n$ は自己共役でノルムが $1$ 以下なので $I-B^n\ge0$ であり、$R\ge0$ となる。積級数は絶対収束するので、二項級数を二乗した係数恒等式から $R^2=A$ である。級数表示から可換性を得る。

別の正平方根 $Q$ は $A=Q^2$ と可換なので $R$ とも可換である。$(R+Q)(R-Q)=0$ である。正作用素 $C$ について $\inner{Cx}{x}=0$ なら、$x+ty$ に対する非負性の一次項を実数と純虚数の $t$ で調べて $\inner{Cx}{y}=0$ を得る。よって $\Ker(R+Q)=\Ker R\cap\Ker Q$ である。$z=(R-Q)x$ はこの核に入り、その一方で $\norm z^2=\inner{x}{(R-Q)z}=0$ である。従って $R=Q$ である。
:::

:::theorem[極分解]{#thm-fa-polar}
$T\in\mathcal B(H,K)$ とし $\abs T=\sqrt{T^*T}$ とおく。一意な $U$ が存在し、$T=U\abs T$、$\Ker T$ 上で $U=0$、かつ $U$ は $(\Ker T)^\perp$ から $\overline{\Img T}$ への全射等長写像となる。この $U$ を :term[部分等長作用素]{en="partial isometry" yomi="ぶぶんとうちょうさようそ"} という。
:::

:::proof
$\norm{\abs T x}^2=\norm{Tx}^2$ より $\abs T x\mapsto Tx$ は良定義な線形等長写像である。稠密部分空間の等長写像として連続的に閉包へ延長する。その像は完備で閉じており、$\Img T$ を含むので $\overline{\Img T}$ である。$\overline{\Img\abs T}=(\Ker\abs T)^\perp=(\Ker T)^\perp$ なので、残る直交成分で零と定める。稠密部分空間と核での値が指定されているため一意である。
:::

:::example[符号と大きさ]
四つの位相を持つ先のユニタリ作用素 $D=\diag(1,-1,i,-i)$ では $D^*D=I$ なので $\abs D=I$、$U=D$ である。一方、$T=\diag(-4,0,3i)$ のとき $\abs T=\diag(4,0,3)$、$U=\diag(-1,0,i)$ である。$U$ は核の方向を除いて長さを保つ。右シフト $S$ では $\abs S=I$、$U=S$ であり、全空間で長さを保っても値域は真部分空間になりうる。
:::

:::proposition[正規作用素の極分解と可換性]{#prop-fa-normal-polar}
正規作用素 $T\in\mathcal B(H)$ の極分解 $T=U\abs T$ では、$U$ と $\abs T$ は可換である。初期空間と終空間はともに $(\Ker T)^\perp$ であり、そこでは $U$ はユニタリである。
:::

:::proof
正規性から $T$ と $T^*T$ は可換である。正平方根の構成から $T$ は $\abs T$ とも可換で、$(U\abs T-\abs T U)\abs T=T\abs T-\abs T T=0$ である。従って差は $\Img\abs T$ の閉包で零である。$\Ker\abs T=\Ker T$ 上では $U$ と $\abs T$ がともに零なので、全空間で可換となる。また $\Ker T=\Ker T^*$ だから $\overline{\Img T}=(\Ker T^*)^\perp=(\Ker T)^\perp$ であり、極分解での等長写像がこの空間のユニタリ作用素となる。
:::

### 極分解から共役の極分解を求める

:::proposition[共役と絶対値]{#prop-fa-polar-adjoint}
$T=U\abs T$ を :ref[thm-fa-polar] の極分解とする。このとき $\abs{T^*}=U\abs T U^*$、$T^*=U^*\abs{T^*}$ であり、後者は $T^*$ の極分解である。
:::

:::proof
$U^*U$ は初期空間 $(\Ker T)^\perp$ への直交射影で、その射影は $\abs T$ の左右に掛けても $\abs T$ を変えない。$U\abs T U^*$ は正であり、その二乗は $U\abs T^2U^*=TT^*$ である。従って正平方根の一意性から第一式が従う。これに左から $U^*$ を掛けると $\abs T U^*=T^*$ になる。$U^*$ は $\overline{\Img T}=(\Ker T^*)^\perp$ で等長であり、補空間で零である。従って :ref[thm-fa-polar] の特徴付けを満たす。
:::

:::example[二方向の絶対値]
$T(x,y)=(3y,0)$ なら $\abs T$ は第 2 成分の方向に固有値 $3$ を持つ。共役 $T^*(x,y)=(0,3x)$ では、$\abs{T^*}$ の固有値 $3$ の方向は第 1 成分へ移る。二つの正平方根が一致しないことと $T^*T\ne TT^*$ は同じ情報である。
:::

:::example[核のある行列の極分解を方向で求める]
標準内積を入れた $\C^3$ で $a=(2,1,i)$、$b=(1,- i,1)$ とする。$Tu=\inner{u}{b}a$ は、入力を一つの数にしてから方向 $a$ を付ける写像である。行列は $T=ab^*=\begin{pmatrix}2 & 2 i & 2\\1 & i & 1\\i & -1 & i\end{pmatrix}$ である。例えば入力 $x=(2,1,0)$ では $\inner{x}{b}=2 + i$ である。$T^*T$ には $a^*a=\norm a^2$、$T^2$ には $b^*a=c$ が現れる。

$$
\begin{aligned}
Tx&=\begin{pmatrix}4 + 2 i\\2 + i\\-1 + 2 i\end{pmatrix},\\
T^*T&=\begin{pmatrix}6 & 6 i & 6\\- 6 i & 6 & - 6 i\\6 & 6 i & 6\end{pmatrix},\\
T^2&=\begin{pmatrix}4 + 4 i & -4 + 4 i & 4 + 4 i\\2 + 2 i & -2 + 2 i & 2 + 2 i\\-2 + 2 i & -2 - 2 i & -2 + 2 i\end{pmatrix}.
\end{aligned}
$$

$T^2=cT$ の係数 $c=\inner{a}{b}=2 + 2 i$ は、出力に再び作用した際に得られる数である。従って零以外の固有値があれば $c$ だけである。実際、$Tu=\lambda u$ で $\lambda\ne0$ なら $u\in\spn a$ であり、$Ta=ca$ である。一方 $b^\perp$ は全て核に入る。

初期空間は $\spn b$、終空間は $\spn a$ である。$P_b=bb^*/\norm b^2$ に対して $\abs T=\norm a\norm b\,P_b$ となる。実際にこの二乗は $T^*T$ であり、直交射影の非負倍なので正である。$U=(a/\norm a)(b/\norm b)^*$ と置くと $U(b/\norm b)=a/\norm a$、$U|_{b^\perp}=0$ である。

従って $U\abs T=T$、$U^*U=P_b$、$UU^*=P_a$ となる。この例で $T$ の絶対値を取ると、大きさ $3 \sqrt{2}$ と入力側の方向が残る。出力側の方向を復元するには $U$ を続けて掛ける。
:::

:::example[固有値の位相と長さを分ける]
$\C^4$ の標準座標で $D=\diag(-2,0,3 i,1)$、$x=(1,1 - i,2,i)$ とする。作用は各成分を $\lambda_j$ 倍するので、計算を座標ごとに進められる。

| $j$ | $\lambda_j$ | $x_j$ | $(Dx)_j$ | $\abs{(Dx)_j}^2$ |
| --- | --- | --- | --- | --- |
| 1 | $-2$ | $1$ | $-2$ | $4$ |
| 2 | $0$ | $1 - i$ | $0$ | $0$ |
| 3 | $3 i$ | $2$ | $6 i$ | $36$ |
| 4 | $1$ | $i$ | $i$ | $1$ |

これを足すと $\norm{Dx}^2=41$ である。一般の $x$ に対して $\norm{Dx}^2\le 9\norm x^2$ であり、係数の絶対値が最大の座標を一つ選べば等号を得る。そのため $\norm D=3$ である。複数の座標で最大値を取る場合、最大値に対応する部分空間のどの単位ベクトルでも同じノルムに達する。

共役は $D^*=\begin{pmatrix}-2 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & - 3 i & 0\\0 & 0 & 0 & 1\end{pmatrix}$ である。実際、$\inner{Dx}{v}=\sum_j \lambda_jx_j\conj{v_j}$ は $\inner{x}{D^*v}$ と等しい。$D^*D=DD^*=\diag(4,0,9,1)$ なので $D$ は正規である。$\lambda_j$ が全て実なら自己共役であり、その上で全てが非負なら正である。また全てが絶対値 $1$ なら $D^*D=DD^*=I$ なのでユニタリになる。この例では定義の行列等式を成分ごとの条件で判定できる。

極分解を求めるときは、絶対値 $H=\abs D=\begin{pmatrix}2 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 3 & 0\\0 & 0 & 0 & 1\end{pmatrix}$ と $U=\begin{pmatrix}-1 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & i & 0\\0 & 0 & 0 & 1\end{pmatrix}$ を取る。$H\ge0$、$H^2=D^*D$ である。成分ごとに $UH=D$ となり、$U^*U=UU^*=\begin{pmatrix}1 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{pmatrix}$ は $D$ によって消されない座標への射影である。$U$ はその部分空間上で各座標の位相だけを変える。例えば零の成分に任意の絶対値 $1$ の値を入れても積 $UH$ は変わらないが、$U$ を核上で零にする条件によって極分解の $U$ が定まる。
:::

:::proposition[極分解と閉値域]{#prop-fa-polar-closed-range}
$T\in\mathcal B(H,K)$ の値域が閉じていることと、$\abs T$ の値域が閉じていることは同値である。閉じている場合、$\abs T$ は $(\Ker T)^\perp$ から同じ空間への有界な線形全単射である。
:::

:::proof
$T=U\abs T$ の $U$ は $(\Ker T)^\perp$ から $\overline{\Img T}$ への等長全射である。$\Img\abs T$ はその定義域側に含まれ、像は $\Img T$ である。等長写像は閉性を保存するため同値性を得る。$\abs T$ の値域の閉包は $(\Ker T)^\perp$ であり、閉性のもとで一致する。そこに制限した作用素は核が零で全射である。有界全単射の逆の有界性は、[有界逆写像定理](/math/functional-analysis/baire-theorems/) で扱う。例えば $T=\diag(1/n)$ では初期空間は全 $\ell^2$ だが、$(1/n)_n$ は値域に入らないため逆を全空間に定められない。
:::

:::corollary[正作用素の平方根の下界と逆]{#cor-fa-positive-root-lower-bound}
有界正作用素 $A$ が $A\ge cI$、$c>0$ を満たすなら $A^{1/2}\ge\sqrt c I$、$A^{-1}\ge0$、$\norm{A^{-1}}\le1/c$ である。
:::

:::proof
$R=A^{1/2}$ とする。$\norm{Rx}^2=\inner{Ax}{x}\ge c\norm x^2$ のため、$R$ は単射かつ閉値域である。自己共役性から値域は稠密でもあり、$R^{-1}$ を持つ。同じ議論で $R+\sqrt c I$ は可逆となる。$(R+\sqrt c I)^{-1}$ は正である。実際 $x=(R+\sqrt c I)y$ に対し逆の二次形式は $\inner{y}{(R+\sqrt c I)y}\ge0$ である。従って可換な正作用素の積 $(A-cI)(R+\sqrt c I)^{-1}=R-\sqrt c I$ は正である。$A^{-1}=R^{-2}=(R^{-1})^*R^{-1}$ も正で、ノルムは $\norm{R^{-1}}^2\le1/c$ である。
:::

:::example[複素係数の部分等長作用素を構成する]
標準内積を入れた $\C^4$ で $a=(1,- i,2,0)$、$b=(2,1,0,-1)$ とする。$Tu=\inner{u}{b}a$ は、入力を一つの数にしてから方向 $a$ を付ける写像である。行列は $T=ab^*=\begin{pmatrix}2 & 1 & 0 & -1\\- 2 i & - i & 0 & i\\4 & 2 & 0 & -2\\0 & 0 & 0 & 0\end{pmatrix}$ である。例えば入力 $x=(1,2,i,3)$ では $\inner{x}{b}=1$ である。$T^*T$ には $a^*a=\norm a^2$、$T^2$ には $b^*a=c$ が現れる。

$$
\begin{aligned}
Tx&=\begin{pmatrix}1\\- i\\2\\0\end{pmatrix},\\
T^*T&=\begin{pmatrix}24 & 12 & 0 & -12\\12 & 6 & 0 & -6\\0 & 0 & 0 & 0\\-12 & -6 & 0 & 6\end{pmatrix},\\
T^2&=\begin{pmatrix}4 - 2 i & 2 - i & 0 & -2 + i\\-2 - 4 i & -1 - 2 i & 0 & 1 + 2 i\\8 - 4 i & 4 - 2 i & 0 & -4 + 2 i\\0 & 0 & 0 & 0\end{pmatrix}.
\end{aligned}
$$

$T^2=cT$ の係数 $c=\inner{a}{b}=2 - i$ は、出力に再び作用した際に得られる数である。従って零以外の固有値があれば $c$ だけである。実際、$Tu=\lambda u$ で $\lambda\ne0$ なら $u\in\spn a$ であり、$Ta=ca$ である。一方 $b^\perp$ は全て核に入る。

初期空間は $\spn b$、終空間は $\spn a$ である。$P_b=bb^*/\norm b^2$ に対して $\abs T=\norm a\norm b\,P_b$ となる。実際にこの二乗は $T^*T$ であり、直交射影の非負倍なので正である。$U=(a/\norm a)(b/\norm b)^*$ と置くと $U(b/\norm b)=a/\norm a$、$U|_{b^\perp}=0$ である。従って $U\abs T=T$、$U^*U=P_b$、$UU^*=P_a$ となる。この例で $T$ の絶対値を取ると、大きさ $6$ と入力側の方向が残る。出力側の方向を復元するには $U$ を続けて掛ける。
:::

## Lax–Milgram の定理と Stampacchia の定理

変分問題では形式 $a(u,v)$ を、未知の $u$ を作用素で動かしてから $v$ と内積を取る形にする。強圧性は、その作用がどの方向でも零に近づきすぎないことを保証する。

:::theorem[Stampacchia の定理]{#thm-fa-stampacchia}
$C$ を Hilbert 空間 $H$ の空でない閉凸集合、$a$ を実双線形または第 1 引数について線形な複素半双線形形式とする。
$$
\abs{a(u,v)}\le M\norm u\norm v,\qquad
\Re a(u,u)\ge\alpha\norm u^2\quad(\alpha>0)
$$
を仮定し、後者を :term[強圧性]{en="coercivity" yomi="きょうあつせい"} と呼ぶ。任意の有界な共役線形汎関数 $F$（実の場合は線形）に対して
$$
\Re a(u,v-u)\ge\Re F(v-u)\qquad(v\in C)
$$
を満たす $u\in C$ が一意に存在する。
:::

:::proof
Riesz の定理を共役した形により $a(u,v)=\inner{Au}{v}$、$F(v)=\inner{f}{v}$ と書ける。一意性から $A$ は線形で $\norm A\le M$ である。$P_C$ の特徴付け :eqref[eq-fa-variational-projection] により変分不等式は $u=P_C(u-\tau(Au-f))$ と同値である。零空間なら結論は直ちに成立する。正数 $0<\tau<2\alpha/M^2$ に対し、
$$
\norm{(I-\tau A)w}^2
=\norm w^2-2\tau\Re\inner{Aw}{w}+\tau^2\norm{Aw}^2
\le(1-2\tau\alpha+\tau^2M^2)\norm w^2.
$$
従って $P_C$ の非拡大性と合わせて右辺は $C$ 上の縮小写像である。$C$ は完備なので :ref[thm-fa-contraction] を適用できる。
:::

:::corollary[Lax–Milgram の定理]{#cor-fa-lax-milgram}
同じ形式 $a$ と汎関数 $F$ について、$a(u,v)=F(v)$ が全ての $v\in H$ で成り立つ $u$ が一意に存在し、$\norm u\le\norm F/\alpha$ である。
:::

:::proof
$C=H$ とすると、$v=u+w$ と $v=u-w$ を比較して $\Re(a(u,w)-F(w))=0$ となる。複素の場合は $w$ を $iw$ に置き換えると虚部も零となる。$v=u$ で $\alpha\norm u^2\le\norm F\norm u$ を得る。
:::

:::example[非対称な一次方程式と制約]
$A=\begin{pmatrix}2&1\\-1&2\end{pmatrix}$ なら 実空間 $\R^2$ では $\inner{Au}{u}=2\norm u^2$ である。$f=(1,0)$ に対する無制約解は $(2/5,1/5)$ となる。一方 $C=\{(x,0):x\ge0\}$ に制限すると $u=(1/2,0)$ であり、$Au-f=(0,-1/2)$ が全ての $v-u$ と直交するので変分不等式を満たす。
:::

:::example[境界条件から変分問題を作る]
実 Hilbert 空間 $H=\{u:u(t)=\int_0^t g(s)\,ds,\ g\in L^2[0,1]\}$ に $\inner{u}{v}_H=\int_0^1u'v'$ を入れる。ここで $u'=g$ はほとんど至る所での導関数であり、$g$ は積分の基本定理で一意である。$u\mapsto g$ が $L^2$ への等長同型なので $H$ は完備である。$\abs{u(t)}\le\sqrt t\,\norm u_H$ より $\norm u_2^2\le\norm u_H^2/2$ である。そこで $a(u,v)=\int_0^1(u'v'+uv)$、$F(v)=\int_0^1v$ とすれば $M=3/2$、$\alpha=1$ を使え、$F$ も有界である。Lax–Milgram により $a(u,v)=F(v)$ の解が一意に存在する。

この例では $u(t)=1-\cosh(1-t)/\cosh1$ が解になる。実際、$u(0)=0$、$u'(1)=0$、$-u''+u=1$ が成り立つ。積分表示を持つ任意の $v\in H$ に対し部分積分すると、端点の寄与は $\left[u'v\right]_0^1=0$ であり、$a(u,v)=\int_0^1(-u''+u)v=F(v)$ を得る。$u(0)=0$ は空間の定義に含まれ、$u'(1)=0$ は部分積分で境界項を消す条件である。
:::

:::example[変分方程式で残差のモーメントを消す]
区間 $[0,2]$ の Lebesgue 測度で内積を取る。関数 $f(t)=t^{3} + t$ の $M=\spn\{1,t^{2} + t,t^{2}\}$ への直交成分を計算する。ここでは多項式を係数だけで比較せず、積分で内積を求める。例えば基底 $v_j$ に対して $G_{kj}=\int_a^b v_j\overline{v_k}$ と置く。成分の積は多項式なので、$\int_a^b t^d\,dt=(b^{d+1}-a^{d+1})/(d+1)$ を項ごとに使える。

$$
\begin{aligned}
G&=\begin{pmatrix}2 & \frac{14}{3} & \frac{8}{3}\\\frac{14}{3} & \frac{256}{15} & \frac{52}{5}\\\frac{8}{3} & \frac{52}{5} & \frac{32}{5}\end{pmatrix},\\
b=(\inner{f}{v_k})_k&=\begin{pmatrix}6\\\frac{356}{15}\\\frac{44}{3}\end{pmatrix},\\
c=G^{-1}b&=\begin{pmatrix}\frac{2}{5}\\- \frac{7}{5}\\\frac{22}{5}\end{pmatrix},\\
p(t),\qquad r(t)&=3 t^{2} - \frac{7 t}{5} + \frac{2}{5},\qquad t^{3} - 3 t^{2} + \frac{12 t}{5} - \frac{2}{5}.
\end{aligned}
$$

$r=f-p$ の各モーメントは $\int_a^b r\overline{v_j}=0$ となる。従って $p$ が射影であり、$\norm{f-m}^2=\norm r^2+\norm{p-m}^2$ を全ての $m\in M$ に対して得る。積分で確かめても $\norm f^2=\frac{3544}{105}$、$\norm p^2=\frac{2528}{75}$、$\norm r^2=\frac{8}{175}$ となり、三つの値がこの関係を満たす。

別の計算法は、同じ基底を直交化して内積の表を作ることである。第 4 列は正規化前の $q_j$ を掛けるための係数で、第 5 列はその方向までを残した場合の誤差平方である。

| 段階 | 直交方向 $q_j$ | $\norm{q_j}^2$ | $\inner{f}{q_j}/\norm{q_j}^2$ | 残差の平方 |
| --- | --- | --- | --- | --- |
| 1 | $1$ | $2$ | $3$ | $\frac{1654}{105}$ |
| 2 | $t^{2} + t - \frac{7}{3}$ | $\frac{278}{45}$ | $\frac{219}{139}$ | $\frac{6088}{14595}$ |
| 3 | $\frac{45 t^{2}}{139} - \frac{94 t}{139} + \frac{34}{139}$ | $\frac{8}{417}$ | $\frac{22}{5}$ | $\frac{8}{175}$ |

残差の直交性を変分問題として読むと、$p\in M$ を未知とする条件は全ての $v\in M$ について $\inner{p}{v}=\inner{f}{v}$ である。有界性と強圧性の定数は共に $1$ であり、解のノルムは $\norm f$ 以下である。もし許される候補を $M$ 内の閉凸集合へ狭めれば、一般には等式を全方向で課せない。その場合は許される変位 $v-p$ ごとの実部の不等式を使う。
:::

:::proposition[強圧形式による Galerkin 近似の誤差]{#prop-fa-galerkin-error}
$a$ を :ref[cor-fa-lax-milgram] の形式、$u$ を解とする。任意の閉部分空間 $H_0$ に対し、$a(u_0,v)=F(v)$ が全ての $v\in H_0$ で成り立つ $u_0\in H_0$ が一意に存在し、
$$
\norm{u-u_0}\le(M/\alpha)\inf_{w\in H_0}\norm{u-w}
$$
となる。これを :term[Galerkin 近似]{en="Galerkin approximation" yomi="がれるきんきんじ"} と呼ぶ。
:::

:::proof
$H_0$ は Hilbert 空間であり、制限しても同じ有界性と強圧性が成り立つので定理を適用できる。差 $e=u-u_0$ は全ての $v\in H_0$ に対し $a(e,v)=0$ となる。任意の $w\in H_0$ について $w-u_0\in H_0$ だから $a(e,e)=a(e,u-w)+a(e,w-u_0)=a(e,u-w)$ である。従って $\alpha\norm e^2\le M\norm e\norm{u-w}$ であり、$e=0$ の場合を含めて結論を得る。$a$ が内積なら射影の計算になるが、証明は形式の対称性を使っていない。
:::

:::example[非対称な強圧形式で近似空間を変える]
実空間 $\R^2$ で $A=\left(\begin{smallmatrix}3&2\\-2&3\end{smallmatrix}\right)$、$f=(1,2)$ とする。$a(u,v)=\inner{Au}{v}$ とすれば $\inner{Au}{u}=3\norm u^2$ であり、$A^{\mathsf T}A=13I$ なので $\alpha=3,\ M=\sqrt{13}$ を選べる。解は $u=(-1,8)/13$ である。実際に $A(-1,8)=(13,26)$ になる。

$H_0=\spn e_2$ では $u_0=(0,2/3)$ である。近似誤差は $u-u_0=(-1,-2/3)/13$ で、ノルムは $\sqrt{13}/39$ となる。最良の直交近似は $v=(0,8/13)$、その誤差は $1/13$ である。従って Galerkin の誤差評価の係数 $M/\alpha=\sqrt{13}/3$ がこの組では等号になる。さらに $C=\{(0,s):s\le1/2\}$ に制約すると、解は $u_C=(0,1/2)$ になる。$Au_C-f=(0,-1/2)$ なので、許される $v-u_C$ の第 2 成分 $s-1/2\le0$ と内積を取ると非負になる。
:::

:::proposition[対称な強圧形式のエネルギー最小化]{#prop-fa-energy-minimization}
実 Hilbert 空間の連続対称双線形形式 $a$ が強圧的であり、$F$ が有界線形であるとする。空でない閉凸集合 $C$ 上で $J(v)=a(v,v)/2-F(v)$ の一意な最小点は Stampacchia の不等式の解である。
:::

:::proof
$u,v\in C$ に対し $J(v)-J(u)=a(u,v-u)-F(v-u)+a(v-u,v-u)/2$ である。変分不等式の解では差は $(\alpha/2)\norm{v-u}^2$ 以上となるので一意な最小点である。逆に最小点 $u$ に対し $u+t(v-u)$ を比較候補に取る。式を $t>0$ で割って $t\downarrow0$ とすると変分不等式になる。存在は :ref[thm-fa-stampacchia] で保証される。非対称形式では対角 $a(v,v)$ を微分すると対称部分だけが現れるため、この最小化には形式の対称性を使っている。
:::

:::example[制約付きのエネルギーを二次式から解く]
$\R^2$ で $A=\left(\begin{smallmatrix}2&1\\1&2\end{smallmatrix}\right)$、$F(v)=v_1-2v_2$、$C=\{v_2\ge0\}$ とする。行列の二次形式は $v_1^2+v_2^2+(v_1+v_2)^2$ なので強圧定数 $1$ を使える。無制約解は $A^{-1}(1,-2)=(4/3,-5/3)$ であり制約に入らない。制約面 $v_2=0$ で第 1 成分を最適化すると $2v_1=1$ だから候補は $u=(1/2,0)$ である。$Au-f=(0,5/2)$ なので、任意の $v\in C$ に対し $\inner{Au-f}{v-u}=(5/2)v_2\ge0$ である。この検算で全ての許容方向を扱えており、最小値は $J(u)=-1/4$ となる。

一方、無制約解を Euclid ノルムで $C$ に射影した点は $(4/3,0)$ であり、同じ値にならない。これは形式 $a$ が座標の交差項を含むためである。$J(v)=J(A^{-1}f)+a(v-A^{-1}f,v-A^{-1}f)/2$ と平方完成すれば、正しい最近点は $a$ から得る距離で決まると分かる。
:::

:::example[射影反復の刻み幅と誤差を決める]
実 Hilbert 空間における強圧形式の制約付き問題を、$u_{n+1}=P_C(u_n-\tau(Au_n-f))$ で解く。$\alpha=2$、$M=3$ の評価が使えるなら、縮小率の平方の上界は $1-4\tau+9\tau^2$ である。これを $9(\tau-2/9)^2+5/9$ と平方完成すると、$\tau=2/9$ で最小となり、縮小率 $q=\sqrt5/3$ が使える。初期誤差の評価 $q^n\norm{u_1-u_0}/(1-q)$ に加え、途中の $u_n$ からの残差 $\norm{u_{n+1}-u_n}$ でも誤差を評価できる。この残差を $1-q$ で割った値が $u_n$ と解の距離の上界になる。

例えば $A=\diag(2,3)$、$C=\{(x,y):x,y\ge0\}$、$f=(1,-3)$ なら、射影は各成分を零との最大値に置き換える。$u_0=0$ からは第 2 成分が常に零で、第 1 成分 $x_n$ は $x_{n+1}=(5/9)x_n+2/9$ を満たす。従って $x_n=(1-(5/9)^n)/2$ であり、解は $(1/2,0)$ になる。$Au-f=(0,3)$ なので変分条件を直接検算できる。この例の真の誤差は $(5/9)^n/2$ で、すべての強圧形式に共通の上の率より小さい。対角行列とこの制約では各成分を独立に追えるためである。
:::

## 演習問題

::::exercise{level=basic}
$\begin{pmatrix}1&i\\2&0\end{pmatrix}$ の共役を求めよ。
:::solution
行列成分 $a_{ij}=\inner{Ae_j}{e_i}$ を用いれば $(A^*)_{ij}=\conj{a_{ji}}$ であり、$\begin{pmatrix}1&2\\-i&0\end{pmatrix}$ となる。
:::
::::

::::exercise{level=basic}
$A=\begin{pmatrix}2&1\\1&2\end{pmatrix}$ が正であることと、$\norm A=3$ を示せ。
:::solution
二次形式は $\norm x^2+\abs{x_1+x_2}^2$ であり非負である。$\abs{x_1+x_2}^2\le2\norm x^2$ より上界は $3$、$(1,1)/\sqrt2$ で等号になる。
:::
::::

::::exercise{level=basic}
$A=\diag(9,4)$ の正平方根を求めよ。
:::solution
$\diag(3,2)$ は正で二乗が $A$ になる。正平方根の一意性によりこれが答えである。
:::
::::

::::exercise{level=basic}
$U=2P-I$ が直交射影 $P$ に対して自己共役かつユニタリであることを示せ。
:::solution
$P^*=P$ より $U^*=U$ である。$U^2=4P^2-4P+I=I$ なので $U^*U=UU^*=I$ である。
:::
::::

::::exercise{level=advanced}
$T=\begin{pmatrix}0&3\\0&0\end{pmatrix}$ の極分解を求めよ。
:::solution
$T^*T=\diag(0,9)$、$\abs T=\diag(0,3)$ である。$Ue_1=0$、$Ue_2=e_1$ とおくと $U\abs T=T$ で、$U$ は $\spn e_2$ で等長である。一意性からこれが極分解である。
:::
::::

::::exercise{level=advanced}
可換な正作用素 $A,B$ の積が正であることを示せ。
:::solution
$B$ は $\sqrt A$ と可換なので $AB=\sqrt A B\sqrt A$ である。右辺は自己共役で二次形式が $\inner{B\sqrt A x}{\sqrt A x}\ge0$ となる。
:::
::::

::::exercise{level=advanced}
$A=\diag(2,4)$、$f=(-1,3)$、$C=\{(x,y):x,y\ge0\}$ とする。対応する変分不等式を解け。
:::solution
$u=(0,3/4)$ なら $Au-f=(1,0)$ である。任意の $v=(x,y)\in C$ に対し $\inner{Au-f}{v-u}=x\ge0$ なので条件を満たす。強圧定数 $2$ によって解は一意である。
:::
::::

::::exercise{level=advanced}
強圧形式に対応する $A$ について、Lax–Milgram の定理を変分不等式を経由せず示せ。
:::solution
$\norm{Au}\ge\alpha\norm u$ なので単射で、$Au_n$ が収束すれば $u_n$ は Cauchy 列になる。完備性と連続性から $\Img A$ は閉じる。また $\Re\inner{A^*u}{u}=\Re\inner{Au}{u}\ge\alpha\norm u^2$ なので $\Ker A^*=0$ である。:ref[cor-fa-adjoint-ranges] により $\Img A$ は稠密でもあるので全射である。$Au=f$ を解けば所望の等式になる。
:::
::::

::::exercise{level=advanced}
$L^2[0,2]$ で $M=\spn\{1,t^2,t^3\}$、$F(v)=\inner{t-t^4}{v}$ とする。形式 $a(u,v)=\inner{u}{v}$ を用いて $M$ 上の変分方程式を解き、元の関数との差を平方ノルムで求めよ。
:::solution
**変分方程式を係数に直す**

区間 $[0,2]$ の Lebesgue 測度で内積を取る。関数 $f(t)=- t^{4} + t$ の $M=\spn\{1,t^{2},t^{3}\}$ への直交成分を計算する。ここでは多項式を係数だけで比較せず、積分で内積を求める。例えば基底 $v_j$ に対して $G_{kj}=\int_a^b v_j\overline{v_k}$ と置く。成分の積は多項式なので、$\int_a^b t^d\,dt=(b^{d+1}-a^{d+1})/(d+1)$ を項ごとに使える。

$$
\begin{aligned}
G&=\begin{pmatrix}2 & \frac{8}{3} & 4\\\frac{8}{3} & \frac{32}{5} & \frac{32}{3}\\4 & \frac{32}{3} & \frac{128}{7}\end{pmatrix},\\
b=(\inner{f}{v_k})_k&=\begin{pmatrix}- \frac{22}{5}\\- \frac{100}{7}\\- \frac{128}{5}\end{pmatrix},\\
c=G^{-1}b&=\begin{pmatrix}- \frac{1}{35}\\\frac{207}{56}\\- \frac{71}{20}\end{pmatrix},\\
p(t),\qquad r(t)&=- \frac{71 t^{3}}{20} + \frac{207 t^{2}}{56} - \frac{1}{35},\qquad - t^{4} + \frac{71 t^{3}}{20} - \frac{207 t^{2}}{56} + t + \frac{1}{35}.
\end{aligned}
$$

$r=f-p$ の各モーメントは $\int_a^b r\overline{v_j}=0$ となる。従って $p$ が射影であり、$\norm{f-m}^2=\norm r^2+\norm{p-m}^2$ を全ての $m\in M$ に対して得る。積分で確かめても $\norm f^2=\frac{344}{9}$、$\norm p^2=\frac{93589}{2450}$、$\norm r^2=\frac{499}{22050}$ となり、三つの値がこの関係を満たす。

別の計算法は、同じ基底を直交化して内積の表を作ることである。第 4 列は正規化前の $q_j$ を掛けるための係数で、第 5 列はその方向までを残した場合の誤差平方である。

| 段階 | 直交方向 $q_j$ | $\norm{q_j}^2$ | $\inner{f}{q_j}/\norm{q_j}^2$ | 残差の平方 |
| --- | --- | --- | --- | --- |
| 1 | $1$ | $2$ | $- \frac{11}{5}$ | $\frac{6422}{225}$ |
| 2 | $t^{2} - \frac{4}{3}$ | $\frac{128}{45}$ | $- \frac{663}{224}$ | $\frac{319579}{88200}$ |
| 3 | $t^{3} - \frac{15 t^{2}}{8} + \frac{1}{2}$ | $\frac{2}{7}$ | $- \frac{71}{20}$ | $\frac{499}{22050}$ |

残差の直交性を変分問題として読むと、$p\in M$ を未知とする条件は全ての $v\in M$ について $\inner{p}{v}=\inner{f}{v}$ である。有界性と強圧性の定数は共に $1$ であり、解のノルムは $\norm f$ 以下である。もし許される候補を $M$ 内の閉凸集合へ狭めれば、一般には等式を全方向で課せない。その場合は許される変位 $v-p$ ごとの実部の不等式を使う。
:::
::::

::::exercise{level=advanced}
$\C^4$ で $Tu=\inner{u}{b}a$、$a=(2,i,1,0),b=(0,2,-i,1)$ とする。入力 $x=(1,-1,2,i)$ の像、$T^*T$、正平方根と部分等長作用素を求めよ。
:::solution
**階数一の極分解**

標準内積を入れた $\C^4$ で $a=(2,i,1,0)$、$b=(0,2,- i,1)$ とする。$Tu=\inner{u}{b}a$ は、入力を一つの数にしてから方向 $a$ を付ける写像である。行列は $T=ab^*=\begin{pmatrix}0 & 4 & 2 i & 2\\0 & 2 i & -1 & i\\0 & 2 & i & 1\\0 & 0 & 0 & 0\end{pmatrix}$ である。例えば入力 $x=(1,-1,2,i)$ では $\inner{x}{b}=-2 + 3 i$ である。$T^*T$ には $a^*a=\norm a^2$、$T^2$ には $b^*a=c$ が現れる。

$$
\begin{aligned}
Tx&=\begin{pmatrix}-4 + 6 i\\-3 - 2 i\\-2 + 3 i\\0\end{pmatrix},\\
T^*T&=\begin{pmatrix}0 & 0 & 0 & 0\\0 & 24 & 12 i & 12\\0 & - 12 i & 6 & - 6 i\\0 & 12 & 6 i & 6\end{pmatrix},\\
T^2&=\begin{pmatrix}0 & 12 i & -6 & 6 i\\0 & -6 & - 3 i & -3\\0 & 6 i & -3 & 3 i\\0 & 0 & 0 & 0\end{pmatrix}.
\end{aligned}
$$

$T^2=cT$ の係数 $c=\inner{a}{b}=3 i$ は、出力に再び作用した際に得られる数である。従って零以外の固有値があれば $c$ だけである。実際、$Tu=\lambda u$ で $\lambda\ne0$ なら $u\in\spn a$ であり、$Ta=ca$ である。一方 $b^\perp$ は全て核に入る。

初期空間は $\spn b$、終空間は $\spn a$ である。$P_b=bb^*/\norm b^2$ に対して $\abs T=\norm a\norm b\,P_b$ となる。実際にこの二乗は $T^*T$ であり、直交射影の非負倍なので正である。$U=(a/\norm a)(b/\norm b)^*$ と置くと $U(b/\norm b)=a/\norm a$、$U|_{b^\perp}=0$ である。従って $U\abs T=T$、$U^*U=P_b$、$UU^*=P_a$ となる。この例で $T$ の絶対値を取ると、大きさ $6$ と入力側の方向が残る。出力側の方向を復元するには $U$ を続けて掛ける。
:::
::::

## 確認テスト

::::quiz
どの作用素でも正になるものを選べ。

- [x] $T^*T$
- [ ] $T+T^*$
- [ ] $T^2$

:::explanation
$T^*T$ は自己共役で、二次形式は $\norm{Tx}^2$ である。$T=-I$ なら $T+T^*=-2I$、複素空間の $T=iI$ なら $T^2=-I$ なので他は常に正とは限らない。
:::
::::

::::quiz{answer="4" tol="1e-8"}
$A=\diag(4,16)$ に対し $\norm{\sqrt A}$ を求めよ。
:::explanation
正平方根は $\diag(2,4)$ である。平方ノルムは $4\abs{x_1}^2+16\abs{x_2}^2$ で、単位ベクトル $e_2$ で最大 $16$ なのでノルムは $4$ である。
:::
::::

::::quiz
右シフト $S$ の値域の直交補を選べ。

- [x] $\spn e_1$
- [ ] $\{0\}$
- [ ] $e_1^\perp$

:::explanation
値域は第 1 成分が零である列全体なので、その直交補は第 1 成分だけで支えられる列全体である。これは $\Ker S^*$ とも一致する。
:::
::::

::::quiz{answer="0.25" tol="1e-8"}
実数上で $a(u,v)=4uv$、$F(v)=v$ とした Lax–Milgram の解を求めよ。
:::explanation
全ての $v$ で $4uv=v$ が必要である。$v=1$ で $u=1/4$ を得て、代入すればどの $v$ でも成立する。
:::
::::

::::quiz{answer="5" tol="1e-8"}
$\R^2$ で $Tu=\inner{u}{(1,2)}(2,1)$ の作用素ノルムを求めよ。
:::explanation
Cauchy–Schwarz により $\norm{Tu}\le\norm{(1,2)}\norm{(2,1)}\norm u$ である。$u=(1,2)/\sqrt5$ なら等号が成立するため、ノルムは $\sqrt5\sqrt5=5$ である。
:::
::::

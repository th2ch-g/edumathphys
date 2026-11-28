---
title: 複素構造の変形と Calabi–Yau 多様体
description: 複素構造の無限小変形と障害を記述し、Ricci 平坦計量と具体的な射影多様体を調べる
notation:
  - symbol: 'T_X'
    meaning: X の正則接束
  - symbol: '\kappa'
    meaning: Kodaira–Spencer 写像
  - symbol: '\mu'
    meaning: 概複素構造を表す Beltrami 微分
  - symbol: '\omega_\varphi'
    meaning: 実ポテンシャル φ で変えた Kähler 形式
---

## この章の目標

- 族の遷移関数を微分して Kodaira–Spencer 類を求める。
- Maurer–Cartan 方程式から一次変形と二次障害を取り出す。
- Ricci 形式を指定する問題を複素 Monge–Ampère 方程式に直す。
- Yau の定理の仮定と結論を説明する。
- 四次 K3 曲面と五次三次元超曲面の標準束と Hodge 数を計算する。

## 複素構造の変形と Kodaira–Spencer 写像

同じ滑らかな多様体でも、正則座標の貼り方を変えると異なる複素構造が得られる。例えばトーラスでは格子の一つの生成元を動かせる。近くのファイバーを比較するときには、座標の取り替えによる変化を差し引く必要がある。

:::definition[複素多様体の変形]{#def-cg-deformation}
コンパクト複素多様体 $X$ の :term[変形]{en="deformation" yomi="へんけい"} は、複素多様体の間の固有正則沈め込み $\pi:\mathcal X\to S$ と、指定した点 $0\in S$ でのファイバーの同型 $\mathcal X_0\cong X$ からなる。ここではまず滑らかな基底 $S$ を考える。一般の変形理論では非被約な解析空間も基底に許す。
:::

相対座標を $(z_i,t)$ とし、重なりで $z_i=g_{ij}(z_j,t)$ と書く。$v\in T_0S$ をパラメータ方向に微分すると $\partial_v g_{ij}(z_j,0)$ が得られる。これを $i$ の座標で接ベクトルとみなして $X$ 上の重なりに置く。

:::proposition[Kodaira–Spencer 写像]{#prop-cg-kodaira-spencer}
上のベクトル場の族は $\mathcal O(T_X)$ の $1$ 次 Čech コサイクルを定め、相対座標の変更は類を変えない。得られる線形写像 $\kappa:T_0S\to H^1(X,T_X)$ を :term[Kodaira–Spencer 写像]{en="Kodaira–Spencer map" yomi="こだいらすぺんさーしゃぞう"} と呼ぶ。
:::

:::proof
$g_{ik}=g_{ij}\circ g_{jk}$ を $v$ で微分すると、$\partial_v g_{ik}=\partial_v g_{ij}+(\dd_zg_{ij})\partial_v g_{jk}$ である。右の第二項の Jacobian は $j$ 座標のベクトルを $i$ 座標に直しているため、ベクトル場として読むと Čech の条件になる。座標変更 $z_i'=h_i(z_i,t)$ のパラメータ微分は、両端での $\partial_v h_i$ の差を加える。従って変化はコバウンダリとなる。

接束の式では、$0\to T_X\to T_{\mathcal X}|_X\to \mathcal O_X\otimes T_0S\to0$ の連結準同型である。実際、相対座標で $z_i$ を固定して $v$ を持ち上げ、二つの座標で持ち上げを比較すると上のコサイクルを得る。内部での符号は $j$ での持ち上げから $i$ での持ち上げを引く規約に対応する。
:::

任意の一次変形はこの群で分類される。一次とは基底を $\C[\epsilon]/(\epsilon^2)$ で見ることをいう。座標遷移を $g_{ij}(z)+\epsilon c_{ij}(z)$ とすると三重交差の条件はちょうど上のコサイクル条件である。$\epsilon^2=0$ のため高次の補正は生じず、座標変更 $z_i\mapsto z_i+\epsilon b_i(z_i)$ はコバウンダリを加える。従って商は $H^1(X,T_X)$ となる。一次の類が存在することは、実際のパラメータ領域上でその変形が存在することまでは意味しない。

:::example[楕円曲線の格子を動かす]
$\Im\tau>0$ とし、実座標 $(x,y)\bmod\Z^2$ を $z=x+\tau y$ と書く。$\tau$ を $\tau+t$ に変えると $z_t=z+t(z-\bar z)/(\tau-\bar\tau)$ となる。新しい $T^{0,1}$ を $\partial_{\bar z}+\mu_t\partial_z$ で表すと、$z_t$ を消す条件は $\mu_t=t/(\tau-\bar\tau+t)$ である。従って Beltrami 微分の一次係数は $(\tau-\bar\tau)^{-1}\dd\bar z\otimes\partial_z$ になる。:ref[prop-cg-kodaira-spencer] の Čech 規約では、対応する Dolbeault 代表には負号が付く。実際、各チャートで $y$ の局所持ち上げを $y_i$ とすると、遷移のパラメータ微分は $(y_i-y_j)\partial_z$ である。これを $\delta u$ と書くには $u_i=-y_i\partial_z$ を取る。Čech–Dolbeault 二重複体での比較は $c$ を $-\bar\partial u$ に移すので、代表は $-(\tau-\bar\tau)^{-1}\dd\bar z\otimes\partial_z$ となる。非零な定数形式はトーラスの Dolbeault コホモロジーで非零なので、この方向は座標変更だけでは消えない。
:::

### 部分多様体の方程式を動かす場合

滑らかな超曲面 $Y=\{F=0\}$ について、$F+tG=0$ で変形を作る。$G|_Y$ は法直線束 $N_{Y/X}$ の切断とみなせる。実際、$F$ と $G$ は同じ束の切断であり、$\dd F|_Y$ が接束の法剰余をその束に移す。変形の各局所座標で中央ファイバーの点を動かす速度 $v_i$ は $\dd F(v_i)+G=0$ を満たす。従って $v_i$ は $-G|_Y$ の局所的な法束からの持ち上げである。重なりの差 $v_j-v_i$ は $T_Y$ の切断で、:ref[prop-cg-kodaira-spencer] と同じ座標比較を与える。よって方程式を動かした変形の類は $-G|_Y$ を
$$
0\longrightarrow T_Y\longrightarrow T_X|_Y\longrightarrow N_{Y/X}\longrightarrow0
$$
の連結準同型で送ったものである。大域ベクトル場で法方向の変化を実現できるなら、連結類は零となる。その場合は無限小に座標を動かしただけである。

:::example[斉次方程式の座標変更]
$F=\sum_j Z_j^d$ とし、座標を $Z_i\mapsto Z_i+t\sum_j a_{ij}Z_j$ と変える。一次で追加される式は $d\sum_{ij}a_{ij}Z_i^{d-1}Z_j$ である。これは接束の大域切断の像として零の Kodaira–Spencer 類を与える。例えば $a_{01}=1$ だけが非零なら追加項は $dZ_0^{d-1}Z_1$ である。
:::

### 滑らかなファイバーを同一視する

固有な沈め込みでは、基底の小さい座標球上で、ファイバーを滑らかに同一視できる。全空間に滑らかな Riemann 計量を取り、垂直接束の直交補を水平部分束とする。微分 $\mathrm d\pi$ は各点の水平部分を基底の接空間へ同型に移すので、基底の曲線を水平に持ち上げられる。小さい座標球の中心から各点へ伸びる線分を使う。持ち上げは閉じた小球の逆像に留まり、固有性でこの逆像はコンパクトである。このため常微分方程式の解は線分の終点まで続く。逆向きの線分で持ち上げると逆写像が得られ、解のパラメータへの滑らかな依存性から滑らかな積表示になる。これは :term[Ehresmann の局所自明化]{en="Ehresmann local trivialization" yomi="えーれすまんのきょくしょじめいか"} である。

近傍の各ファイバーが Kähler であれば、ここで得た滑らかな積表示で Betti 数が一定であるため、Hodge 数も局所的に一定になる。実際、相対接束に滑らかな計量を選ぶと、Dolbeault Laplacian は固定した滑らかな多様体上で係数が滑らかに変わる楕円型作用素の族になる。調和空間の次元は上半連続である。もし近づくパラメータで中央の次元より大きい調和空間があれば、各調和空間から中央の次元より一本多い正規直交列を取れる。パラメータを小さくすれば同じ楕円型評価定数を使えるので、Rellich の定理から各ベクトルに $L^2$ 極限がある。作用素の係数と内積の収束から、極限も中央で調和かつ正規直交となり矛盾する。各全次数で Hodge 数の和は一定の Betti 数であるため、どれか一つだけが減ることもない。

### 底変換で一次方向を比較する

$\pi:\mathcal X\to S$ の族を正則写像 $f:(B,0)\to(S,0)$ で引き戻す。局所遷移は $g_{ij}(z,f(b))$ となるので、方向 $v$ での微分は $\partial_{\mathrm df_0(v)}g_{ij}$ である。従って引き戻した族の Kodaira–Spencer 写像は $\kappa_\pi\circ\mathrm df_0$ である。ある方向の一次類が零でも、その方向に沿う族が全て同型とは限らない。$f$ が原点で高い次数を持てば、一次成分だけを見たときに変化が消えるためである。

:::example[格子の二次の変化]
$E_t=\C/(\Z+(i+t^2)\Z)$ を $t=0$ の近くで考える。$\tau=i$ と $s=t^2$ を楕円曲線の式に入れると $\mu_t=t^2/(2i+t^2)$ である。従って一次類は零で、二次の Beltrami 係数は $1/(2i)$ である。例えば小さな正実数 $t$ で $i+t^2$ は $i$ の $\mathrm{SL}_2(\Z)$ 軌道に属さない。軌道で虚部が $1$ となるには $c^2+d^2=1$ が必要で、その場合は実部が整数となるためである。
:::

### 方程式の微分と法方向の持ち上げ

:::example[正則方程式の独立性と接方向を計算する]
$F:\C^3\to\C^2$ の成分を次とする。点 $p=(1,1,2)$ におけるファイバー $F=F(p)$ の局所座標を調べる。

$$
\begin{aligned}
F(x,y,z)&=\left[\begin{matrix}- x^{2} - y^{2} + z\\x y - 1\end{matrix}\right] \\
F(p)&=\left[\begin{matrix}0\\0\end{matrix}\right]
\end{aligned}
$$

正則写像の微分は複素 Jacobian である。反正則変数を別の列に追加する必要はない。点での評価と、その核を線形代数で求める。

$$
\begin{aligned}
\mathrm d F&=\left[\begin{matrix}- 2 x & - 2 y & 1\\y & x & 0\end{matrix}\right] \\
\mathrm dF(p)&=\left[\begin{matrix}-2 & -2 & 1\\1 & 1 & 0\end{matrix}\right] \\
\rank\,\mathrm dF(p)&=2 \\
\Ker\,\mathrm dF(p)&=\operatorname{span}(\text{次の列}),\quad \left[\begin{matrix}-1\\1\\0\end{matrix}\right]
\end{aligned}
$$

写像の成分数と階数が等しいため、独立な小行列式に対応する変数を他の変数の正則関数として解ける。従ってファイバーの複素次元は $3-2=1$ である。パラメータを微分して得られる接ベクトルは行列の核に入り、次元が同じなので核を尽くす。

方程式の右辺を動かす場合、同じ行列で接速度 $v$ の像を指定する。任意の二つの持ち上げの差はここで求めた核に属する。
:::

:::example[正則方程式の独立性と接方向を計算する]
$F:\C^3\to\C^2$ の成分を次とする。点 $p=(2,8,16)$ におけるファイバー $F=F(p)$ の局所座標を調べる。

$$
\begin{aligned}
F(x,y,z)&=\left[\begin{matrix}- x y + z\\- x^{3} + y\end{matrix}\right] \\
F(p)&=\left[\begin{matrix}0\\0\end{matrix}\right]
\end{aligned}
$$

正則写像の微分は複素 Jacobian である。反正則変数を別の列に追加する必要はない。点での評価と、その核を線形代数で求める。

$$
\begin{aligned}
\mathrm d F&=\left[\begin{matrix}- y & - x & 1\\- 3 x^{2} & 1 & 0\end{matrix}\right] \\
\mathrm dF(p)&=\left[\begin{matrix}-8 & -2 & 1\\-12 & 1 & 0\end{matrix}\right] \\
\rank\,\mathrm dF(p)&=2 \\
\Ker\,\mathrm dF(p)&=\operatorname{span}(\text{次の列}),\quad \left[\begin{matrix}\frac{1}{32}\\\frac{3}{8}\\1\end{matrix}\right]
\end{aligned}
$$

写像の成分数と階数が等しいため、独立な小行列式に対応する変数を他の変数の正則関数として解ける。従ってファイバーの複素次元は $3-2=1$ である。パラメータを微分して得られる接ベクトルは行列の核に入り、次元が同じなので核を尽くす。

第一、第二式を同時に動かせるが、接ベクトルは両方の微分を同時に消す一方向だけである。
:::

:::example[完全交差の標準束と Chern 数を計算する]
$X\subset\mathbb P^{4}$ を次数 $4$ の横断的な滑らかな完全交差とし、$H=c_1(\mathcal O_X(1))$ とする。法束は $\bigoplus_a\mathcal O_X(d_a)$ であり、次元は $3$ である。

接束の完全列から $c(T_X)c(N_{X/\mathbb P^N})=c(T_{\mathbb P^N}|_X)$ を得る。Euler 列により右辺は $(1+H)^{N+1}$ である。$c(T_X)=1+c_1H+\cdots+c_nH^n$ と書いて係数を比べる。

$$
\begin{aligned}
c(N)&=4 H + 1 \\
c_{1}+(4)&=\binom{5}{1},\qquad c_{1}=1 \\
c_{2}+(4)&=\binom{5}{2},\qquad c_{2}=6 \\
c_{3}+(24)&=\binom{5}{3},\qquad c_{3}=-14 \\
c(T_X)&=- 14 H^{3} + 6 H^{2} + H + 1
\end{aligned}
$$

法束の行列式を取った随伴公式は $K_X=\mathcal O_X(-1)$ である。これは $c(T_X)$ の一次の係数の符号を反転した値にも等しい。さらに $X$ を一般の $3$ 個の超平面で切り、交点数を積で数えると $\int_XH^3=4$ となる。

$$
\begin{aligned}
\int_X c_n(T_X)&=-14\cdot4=-56
\end{aligned}
$$

最高 Chern 類は実接束の Euler 類であるため、ここで得た整数は位相的 Euler 標数にも等しい。法束の次数を知れば、どの次数の多項式がこの超曲面の埋め込み変形を与えるかが分かる。
:::

## 障害と倉西族の概観

固定した滑らかな多様体のもとの分解 $T_\C X=T^{1,0}\oplus T^{0,1}$ を用いる。それに十分近い概複素構造の新 $T^{0,1}$ は $\mu:T^{0,1}\to T^{1,0}$ のグラフである。$\mu$ を :term[Beltrami 微分]{en="Beltrami differential" yomi="べるとらみびぶん"} と呼ぶ。局所座標で $L_j=\partial_{\bar j}+\sum_k\mu_j^k\partial_k$ と置く。

:::proposition[Maurer–Cartan 方程式]{#prop-cg-maurer-cartan}
グラフの概複素構造が可積分である条件は
$$
\bar\partial\mu+\tfrac12[\mu,\mu]=0
\label{eq-cg-maurer-cartan}
$$
である。ここで括弧は $(j,\ell,k)$ 成分で $2\sum_a(\mu_j^a\partial_a\mu_\ell^k-\mu_\ell^a\partial_a\mu_j^k)$ と定める。
:::

:::proof
$[L_j,L_\ell]$ は $\partial_{\bar a}$ の項を持たず、$\partial_k$ の係数は
$$
\partial_{\bar j}\mu_\ell^k-\partial_{\bar\ell}\mu_j^k
+\sum_a(\mu_j^a\partial_a\mu_\ell^k-\mu_\ell^a\partial_a\mu_j^k)
$$
である。グラフと元の $T^{1,0}$ の交わりは零なので、括弧がグラフに属するにはこの係数が全て零でなければならない。逆に消えれば局所枠同士が閉じ、Leibniz 則で全切断の括弧も閉じる。:ref[thm-cg-newlander-nirenberg] から可積分性を得る。
:::

$\mu=t\mu_1+t^2\mu_2+\cdots$ を代入すると $\bar\partial\mu_1=0$、$\bar\partial\mu_2=-[\mu_1,\mu_1]/2$ になる。従って $[\mu_1,\mu_1]/2$ の $H^2(X,T_X)$ における類が二次障害である。括弧と $\bar\partial$ は次数付き Leibniz 則を満たし、局所式を $\bar z$ で微分して確認できる。これにより $\mu_1$ が閉なら括弧も閉である。また閉じた $\mu_1$ に $\bar\partial\xi$ を加えると、括弧の変化は $2[\mu_1,\bar\partial\xi]+[\bar\partial\xi,\bar\partial\xi]$ であり、各項は $\bar\partial$ の像なので類は変わらない。

:term[倉西族]{en="Kuranishi family" yomi="くらにしぞく"} は、コンパクト複素多様体 $X$ に対する半普遍変形である。その基底の芽は $H^1(X,T_X)$ の原点近傍における正則写像 $k:H^1(X,T_X)\to H^2(X,T_X)$ の零点の芽として記述できる。$k(0)=0$、$\dd k(0)=0$ であり、二次項は上の障害である。半普遍とは変形が局所的にこの族から底変換で得られ、接空間での Kodaira–Spencer 写像が同型になることである。底変換の写像そのものの一意性は一般にはない。

証明では $\bar\partial^*\mu=0$ の条件を課して座標変更の自由度を除き、楕円型の Green 作用素で方程式の非調和成分を解く。調和成分が有限次元の方程式 $k=0$ になる。収束と任意の族の半普遍性には Banach 空間での解析と正則性が要るため、残差が調和成分だけで決まることは下で示す。正則な族への移行と半普遍性の証明は :cite[huybrechts-cg] を参照する。特に $H^2(X,T_X)=0$ ならこの方程式の終域が零であり、基底の芽は滑らかである。

:::example[曲線の変形の次元]
曲線では $H^2(T_X)=0$ なので倉西基底は滑らかである。種数 $g\ge2$ では $T_X=K_X^{-1}$ だから双対性で $h^1(T_X)=h^0(K_X^2)$ である。$K_X^2$ の次数は $4g-4>2g-2$ なので Riemann–Roch から $h^0=3g-3$ となる。$g=2$ では次元 $3$ である。トーラスでは $T_X$ は自明なので次元は $1$ になる。
:::

### 三次まで展開して障害の順序を確かめる

対称な次数一の括弧を使うと、$\mu=t\mu_1+t^2\mu_2+t^3\mu_3$ の三次の係数は $\bar\partial\mu_3+[\mu_1,\mu_2]$ である。二次まで解けている場合、次数付き Leibniz 則によりその括弧の $\bar\partial$ は、符号を除いて $[\mu_1,[\mu_1,\mu_1]]/2$ である。これは次数付き Jacobi 恒等式に $\mu_1$ を三度入れると零になる。従って三次にもコホモロジー類が生じる。ただし $\mu_2$ の選択を変えるとこの類は変わりうる。この段階では、一次の類だけを使って二次障害と同じ方法で一意の三次類が得られるわけではない。

:::example[平坦トーラスの定数方向では括弧が消える]
$\C^2/\Lambda$ の定数行列 $\mu=\begin{pmatrix}a&b\\c&d\end{pmatrix}$ を使い、$L_j=\partial_{\bar z_j}+\sum_k\mu_j^k\partial_{z_k}$ とする。全係数が定数なので括弧は零であり、十分小さければグラフは共役のグラフと直和になる。$\zeta_k=z_k-\sum_j\mu_j^k\bar z_j$ は $L_j$ に消され、実線形変換として小変形なので逆を持つ。従って $\zeta$ は新しい正則座標であり、格子は同じ線形変換の像に変わる。
:::

### 標準束を保つ方程式と障害空間

:::example[完全交差の標準束と Chern 数を計算する]
$X\subset\mathbb P^{5}$ を次数 $2,4$ の横断的な滑らかな完全交差とし、$H=c_1(\mathcal O_X(1))$ とする。法束は $\bigoplus_a\mathcal O_X(d_a)$ であり、次元は $3$ である。

接束の完全列から $c(T_X)c(N_{X/\mathbb P^N})=c(T_{\mathbb P^N}|_X)$ を得る。Euler 列により右辺は $(1+H)^{N+1}$ である。$c(T_X)=1+c_1H+\cdots+c_nH^n$ と書いて係数を比べる。

$$
\begin{aligned}
c(N)&=\left(2 H + 1\right) \left(4 H + 1\right) \\
c_{1}+(6)&=\binom{6}{1},\qquad c_{1}=0 \\
c_{2}+(8)&=\binom{6}{2},\qquad c_{2}=7 \\
c_{3}+(42)&=\binom{6}{3},\qquad c_{3}=-22 \\
c(T_X)&=- 22 H^{3} + 7 H^{2} + 1
\end{aligned}
$$

法束の行列式を取った随伴公式は $K_X=\mathcal O_X(0)$ である。これは $c(T_X)$ の一次の係数の符号を反転した値にも等しい。さらに $X$ を一般の $3$ 個の超平面で切り、交点数を積で数えると $\int_XH^3=8$ となる。

$$
\begin{aligned}
\int_X c_n(T_X)&=-22\cdot8=-176
\end{aligned}
$$

最高 Chern 類は実接束の Euler 類であるため、ここで得た整数は位相的 Euler 標数にも等しい。随伴で標準束が自明になる場合には、接束値形式を最高次形式で収縮できる。
:::

:::example[完全交差の標準束と Chern 数を計算する]
$X\subset\mathbb P^{6}$ を次数 $2,2,3$ の横断的な滑らかな完全交差とし、$H=c_1(\mathcal O_X(1))$ とする。法束は $\bigoplus_a\mathcal O_X(d_a)$ であり、次元は $3$ である。

接束の完全列から $c(T_X)c(N_{X/\mathbb P^N})=c(T_{\mathbb P^N}|_X)$ を得る。Euler 列により右辺は $(1+H)^{N+1}$ である。$c(T_X)=1+c_1H+\cdots+c_nH^n$ と書いて係数を比べる。

$$
\begin{aligned}
c(N)&=\left(2 H + 1\right)^{2} \left(3 H + 1\right) \\
c_{1}+(7)&=\binom{7}{1},\qquad c_{1}=0 \\
c_{2}+(16)&=\binom{7}{2},\qquad c_{2}=5 \\
c_{3}+(47)&=\binom{7}{3},\qquad c_{3}=-12 \\
c(T_X)&=- 12 H^{3} + 5 H^{2} + 1
\end{aligned}
$$

法束の行列式を取った随伴公式は $K_X=\mathcal O_X(0)$ である。これは $c(T_X)$ の一次の係数の符号を反転した値にも等しい。さらに $X$ を一般の $3$ 個の超平面で切り、交点数を積で数えると $\int_XH^3=12$ となる。

$$
\begin{aligned}
\int_X c_n(T_X)&=-12\cdot12=-144
\end{aligned}
$$

最高 Chern 類は実接束の Euler 類であるため、ここで得た整数は位相的 Euler 標数にも等しい。この三次元例でも標準束が自明であり、同じ障害消滅定理の仮定を満たす。
:::

:::example[積の直線束でコホモロジーが残る次数]
$X=\prod_j\mathbb P^{n_j}$ とし、外部テンソル積の直線束 $L=\mathcal O(1)$ を考える。$a$ のべきでコホモロジー次数を記録し、係数を各空間の次元とする。

一因子では符号がすべて揃った Taylor–Laurent 指数だけが残る。指数が非負の場合は全て $0$ 次に、全て負の場合は因子の最高次に属する。負指数を $-b_i-1$ と書けば、弱い合成の個数を二項係数で求められる。

$$
\begin{aligned}
j=1&\colon\quad \mathbb P^{4},\; d=1,\quad e_i\ge0,\quad \sum_i e_i=1,\quad P_j(a)=5
\end{aligned}
$$

積では各因子の Dolbeault Laplacian の非負固有値を加える。零になるには各因子の固有値が零でなければならず、調和代表の次元と次数は多項式の積で合成される。

$$
\begin{aligned}
P_X(a)&=\prod_j P_j(a)=5 \\
\dim_\C X&=4 \\
\chi(X,L)&=P_X(-1)=5
\end{aligned}
$$

随伴の双対束 $K_X\otimes L^*$ は $\mathcal O(-6)$ である。$d\ge0$ の因子では双対側が最高次になり、$d\le-n-1$ の因子では逆になる。全消滅の範囲はそのままなので、双対は $q$ 次と $4-q$ 次を対応させる。射影空間の Euler 列の中央に現れる $\mathcal O(1)$ は高次コホモロジーを持たない。超曲面に制限した後は法束の列も見る必要がある。
:::

### 調和成分を残して方程式を解く

滑らかな Hermite 計量を選び、$D=\bar\partial$、調和射影を $P$、Green 作用素を $G$ とする。$D$ と $G$ は可換であり、$1=P+DD^*G+D^*DG$ である。ゲージ $D^*\mu=0$ のもと、調和成分 $t=P\mu$ を指定して
$$
\mu=t-\tfrac12D^*G[\mu,\mu]
$$
を先に解く。実次元を $d$ とし、$s>d/2+2$ の Sobolev ノルムを使う。括弧は係数の積に一階微分を一つ含むため、Sobolev の積の評価から $H^s\times H^s\to H^{s-1}$ の有界な双線形写像になる。$G$ は楕円性で二階分を回復し、$D^*$ は一階なので、$D^*G[\cdot,\cdot]$ は $H^s$ に戻る。原点で $\mu\mapsto\mu+\frac12D^*G[\mu,\mu]$ の微分は恒等である。従って Banach 空間の正則逆写像定理により、小さい $t$ に一意な小さい解 $\mu(t)$ が正則に依存する。実際、この式を反復すると二次の係数は $-D^*G[t,t]/2$ になる。$D^*$ を作用させれば $D^*\mu=0$、$P$ を作用させれば $P\mu=t$ となる。

残差を $F=D\mu+[\mu,\mu]/2$、$k=P[\mu,\mu]/2$ と置く。上の式を $D$ で微分し、恒等分解を使うと
$$
F=k+\tfrac12D^*GD[\mu,\mu].
$$
次数付き Leibniz 則と Jacobi 恒等式から $D[\mu,\mu]=2[D\mu,\mu]=2[F,\mu]$ である。ここで $[[\mu,\mu],\mu]=0$ を使った。もし $k=0$ なら $F=D^*G[F,\mu]$ となる。$F$ を $H^{s-1}$ で測ると、積と一階微分の評価で $\|[F,\mu]\|_{H^{s-2}}\le C\|F\|_{H^{s-1}}\|\mu\|_{H^s}$ である。$D^*G$ で一階回復するため、$\mu$ を小さく選び $C'\|\mu\|_{H^s}<1$ とすれば $F=0$ が従う。逆に $F=0$ なら $P D\mu=0$ なので $k=0$ となる。従ってこのゲージで残る方程式は有限次元の $k(t)=0$ であり、その二次項は $P[t,t]/2$ である。族の滑らかさと任意の変形をこのゲージに移すことは、上に挙げた倉西の定理の解析を必要とする。

:::example[二次の純粋な障害の読み方]
調和基底を $\eta_1,\eta_2$、$t=a\eta_1+b\eta_2$ とする。二次の項は
$$
\tfrac12a^2P[\eta_1,\eta_1]+abP[\eta_1,\eta_2]
+\tfrac12b^2P[\eta_2,\eta_2]
$$
である。たとえば最初と最後の項が零で、中央の係数が非零なら、座標軸の二方向は二次まで障害を持たないが、方向 $(1,1)$ は二次に障害を持つ。各方向を独立に検査してから線形結合を作っても、障害は線形にはならない。
:::

## Calabi–Yau 多様体と Calabi 予想

この章では、連結コンパクト Kähler 多様体で $K_X$ が正則に自明なものを :term[Calabi–Yau 多様体]{en="Calabi–Yau manifold" yomi="からびやうたようたい"} と呼ぶ。文献によって単連結性や中間次数の正則微分の消滅を追加する。$K_X$ の自明性は零にならない正則最高次形式 $\Omega$ の存在と同値である。一方、以下の Ricci 平坦計量の存在には実係数で $c_1(T_X)=0$ であれば足りる。

Kähler 形式 $\omega_\varphi=\omega+i\partial\bar\partial\varphi$ の Ricci 形式は :eqref[eq-cg-ricci-form] から
$$
\rho_{\omega_\varphi}=\rho_\omega-i\partial\bar\partial
     \log\frac{\omega_\varphi^n}{\omega^n}
$$
となる。指定した実閉 $(1,1)$ 形式 $\rho'$ が $[\rho']=2\pi c_1(T_X)$ を満たせば、:ref[thm-cg-ddbar] から $\rho_\omega-\rho'=i\partial\bar\partial f$ と書ける。定数を足して $\int_X e^f\omega^n=\int_X\omega^n$ と規格化する。そこで
$$
(\omega+i\partial\bar\partial\varphi)^n=e^f\omega^n,\qquad
\omega+i\partial\bar\partial\varphi>0
\label{eq-cg-monge-ampere}
$$
を解けば $\rho_{\omega_\varphi}=\rho'$ になる。逆に Ricci 形式が一致すれば対数体積比と $f$ の差は $\partial\bar\partial$ で消え、コンパクト性と最大値原理で定数である。体積積分でその定数を決めると方程式に戻る。Kähler 類で体積が一定なのは、差が完全で各形式が閉じているため Stokes の定理から従う。

:::example[自明な標準束が与える体積形式]
$\Omega$ が零にならない正則最高次形式なら $i^{n^2}\Omega\wedge\bar\Omega$ は正の体積形式である。座標の正則係数が非零なので、その対数密度の $\partial\bar\partial$ は零である。よってこの体積形式の定数倍を Kähler 計量の体積にできれば、計量の Ricci 形式は零になる。
:::

### Ricci 平坦性と標準束の平行切断

Ricci 形式は接束の $i\tr F$ であり、標準束 $K_X=\det T_X^*$ の曲率は負の跡である。従って Ricci 平坦では標準束の Chern 接続が平坦である。しかし平坦な直線束には非自明なホロノミーがありうるため、平坦性だけから正則な自明化は従わない。一方 $K_X$ が正則に自明で $\Omega$ が零を持たない正則切断なら、曲率零の条件は $\partial\bar\partial\log\|\Omega\|^2=0$ である。コンパクト連結ならこの実関数は定数で、接続形式 $\partial\log\|\Omega\|^2$ も零である。従って $\Omega$ は平行である。

:::example[平坦な商で最高次形式が降りない]
$E=\C/(\Z+i\Z)$ とし、$E^2$ に $\sigma(z,w)=(z+1/2,-w)$ を作用させる。$\sigma^2=\id$ であり、第一座標の平行移動が固定点を持たないので商は滑らかである。定数係数の積計量は不変であり、商の Ricci 形式は零である。他方、$\sigma^*(\mathrm dz\wedge\mathrm dw)=-\mathrm dz\wedge\mathrm dw$ なので非零な正則最高次形式は商へ降りない。実際、$E^2$ 上の最高次形式の係数はトーラス上の関数で定数だからこの一方向しかない。二乗は不変であり、商上の標準束は非自明な二捩れとなる。
:::

### 標準束の次数から Ricci 類を読む

:::example[完全交差の標準束と Chern 数を計算する]
$X\subset\mathbb P^{3}$ を次数 $3$ の横断的な滑らかな完全交差とし、$H=c_1(\mathcal O_X(1))$ とする。法束は $\bigoplus_a\mathcal O_X(d_a)$ であり、次元は $2$ である。

接束の完全列から $c(T_X)c(N_{X/\mathbb P^N})=c(T_{\mathbb P^N}|_X)$ を得る。Euler 列により右辺は $(1+H)^{N+1}$ である。$c(T_X)=1+c_1H+\cdots+c_nH^n$ と書いて係数を比べる。

$$
\begin{aligned}
c(N)&=3 H + 1 \\
c_{1}+(3)&=\binom{4}{1},\qquad c_{1}=1 \\
c_{2}+(3)&=\binom{4}{2},\qquad c_{2}=3 \\
c(T_X)&=3 H^{2} + H + 1
\end{aligned}
$$

法束の行列式を取った随伴公式は $K_X=\mathcal O_X(-1)$ である。これは $c(T_X)$ の一次の係数の符号を反転した値にも等しい。さらに $X$ を一般の $2$ 個の超平面で切り、交点数を積で数えると $\int_XH^2=3$ となる。

$$
\begin{aligned}
\int_X c_n(T_X)&=3\cdot3=9
\end{aligned}
$$

最高 Chern 類は実接束の Euler 類であるため、ここで得た整数は位相的 Euler 標数にも等しい。この曲面では標準束が負であり、実第一類も零でないので Ricci 平坦性の前提を満たさない。
:::

:::example[完全交差の標準束と Chern 数を計算する]
$X\subset\mathbb P^{3}$ を次数 $6$ の横断的な滑らかな完全交差とし、$H=c_1(\mathcal O_X(1))$ とする。法束は $\bigoplus_a\mathcal O_X(d_a)$ であり、次元は $2$ である。

接束の完全列から $c(T_X)c(N_{X/\mathbb P^N})=c(T_{\mathbb P^N}|_X)$ を得る。Euler 列により右辺は $(1+H)^{N+1}$ である。$c(T_X)=1+c_1H+\cdots+c_nH^n$ と書いて係数を比べる。

$$
\begin{aligned}
c(N)&=6 H + 1 \\
c_{1}+(6)&=\binom{4}{1},\qquad c_{1}=-2 \\
c_{2}+(-12)&=\binom{4}{2},\qquad c_{2}=18 \\
c(T_X)&=18 H^{2} - 2 H + 1
\end{aligned}
$$

法束の行列式を取った随伴公式は $K_X=\mathcal O_X(2)$ である。これは $c(T_X)$ の一次の係数の符号を反転した値にも等しい。さらに $X$ を一般の $2$ 個の超平面で切り、交点数を積で数えると $\int_XH^2=6$ となる。

$$
\begin{aligned}
\int_X c_n(T_X)&=18\cdot6=108
\end{aligned}
$$

最高 Chern 類は実接束の Euler 類であるため、ここで得た整数は位相的 Euler 標数にも等しい。こちらでは標準束が正であり、第一類の符号は三次曲面と逆になる。
:::

:::example[完全交差の標準束と Chern 数を計算する]
$X\subset\mathbb P^{5}$ を次数 $3,3$ の横断的な滑らかな完全交差とし、$H=c_1(\mathcal O_X(1))$ とする。法束は $\bigoplus_a\mathcal O_X(d_a)$ であり、次元は $3$ である。

接束の完全列から $c(T_X)c(N_{X/\mathbb P^N})=c(T_{\mathbb P^N}|_X)$ を得る。Euler 列により右辺は $(1+H)^{N+1}$ である。$c(T_X)=1+c_1H+\cdots+c_nH^n$ と書いて係数を比べる。

$$
\begin{aligned}
c(N)&=\left(3 H + 1\right)^{2} \\
c_{1}+(6)&=\binom{6}{1},\qquad c_{1}=0 \\
c_{2}+(9)&=\binom{6}{2},\qquad c_{2}=6 \\
c_{3}+(36)&=\binom{6}{3},\qquad c_{3}=-16 \\
c(T_X)&=- 16 H^{3} + 6 H^{2} + 1
\end{aligned}
$$

法束の行列式を取った随伴公式は $K_X=\mathcal O_X(0)$ である。これは $c(T_X)$ の一次の係数の符号を反転した値にも等しい。さらに $X$ を一般の $3$ 個の超平面で切り、交点数を積で数えると $\int_XH^3=9$ となる。

$$
\begin{aligned}
\int_X c_n(T_X)&=-16\cdot9=-144
\end{aligned}
$$

最高 Chern 類は実接束の Euler 類であるため、ここで得た整数は位相的 Euler 標数にも等しい。この二方程式では標準束が自明であり、滑らかな射影多様体なので Kähler 計量が存在する。
:::

### Ricci 形式が計量の選択に依存する範囲

曲率の代表を所定の $2\pi c_1(T_X)$ に固定することと、Kähler 類そのものを固定することは別の条件である。例えば Ricci 平坦計量 $g$ を $c>0$ 倍すると、行列式は $c^n$ 倍されるため Ricci 形式は変わらず零である。しかし Kähler 類は $c[\omega]$ になる。$[\omega]\ne0$ なので $c\ne1$ では別の類であり、各類で一意であるという主張には反しない。

:::example[二次元トーラスに異なる平坦類を入れる]
$X=(\C/(\Z+i\Z))^2$ で $\omega_{a,b}=i(a\mathrm dz_1\wedge\mathrm d\bar z_1+b\mathrm dz_2\wedge\mathrm d\bar z_2)$、$a,b>0$ とする。行列式は $ab$ で定数なので全て Ricci 平坦である。二つの因子曲線で積分した値は $2a,2b$ であるため、類の一致には $a,b$ がそれぞれ一致することが必要である。体積は $\int\omega_{a,b}^2/2=4ab$ である。例えば $(a,b)=(1,4)$ と $(2,2)$ は体積がともに $16$ でも、同じ Kähler 類ではない。
:::

## Yau の定理の主張

:::theorem[Calabi 予想に対する Yau の定理]{#thm-cg-yau}
$X$ を連結コンパクト Kähler 多様体、$\omega$ を Kähler 形式とする。滑らかな実関数 $f$ が $\int_X e^f\omega^n=\int_X\omega^n$ を満たすなら、:eqref[eq-cg-monge-ampere] を満たす滑らかな実関数 $\varphi$ が存在し、定数の差を除いて一意である。従って、各 Kähler 類には指定した第一 Chern 類の代表を Ricci 形式に持つ Kähler 計量がただ一つ存在する。
:::

:::proof
一意性を示す。解 $\varphi,\psi$ に対し $u=\varphi-\psi$、$T=\sum_{j=0}^{n-1}\omega_\varphi^j\wedge\omega_\psi^{n-1-j}$ と取る。体積の差を因数分解すると $i\partial\bar\partial u\wedge T=0$ である。$u$ を掛けて積分し、$T$ が閉じていることを使って部分積分すると $\int_X i\partial u\wedge\bar\partial u\wedge T=0$ である。各項は非負で、特に $\omega_\varphi^{n-1}$ の項は $\partial u\ne0$ の点で正である。従って $u$ は定数になる。

一般の存在証明は連続法で非線形楕円型方程式を解く。例えば $e^{tf+c_t}$ の体積を規格化して $t$ を $0$ から動かす。開性では線形化した Laplacian を定数の直交補で逆にする。閉性にはポテンシャル、計量の固有値、高階微分の事前評価を一様に得る必要がある。これらの評価は本章の範囲を超えるため、一般の存在は定理として用いる。証明の構成と文献への案内は :cite[huybrechts-cg] を参照する。前節の Ricci 形式との対応により最後の主張が従う。
:::

:::example[平坦トーラスの一意性]
トーラスの定数係数 Kähler 計量は既に Ricci 平坦である。同じ Kähler 類で別の Ricci 平坦計量があるなら :ref[thm-cg-ddbar] で実ポテンシャルを持つ。その体積比は定数かつ全体積が同じなので $1$ である。従って一意性から元の平坦計量と一致する。
:::

### 曲線の場合は線形方程式になる

$n=1$ では :eqref[eq-cg-monge-ampere] の左辺に高次の外積がなく、
$$
i\partial\bar\partial\varphi=(e^f-1)\omega
$$
となる。右辺は実二次形式で積分が零なので、連結なコンパクト向き付け可能曲面では完全である。Kähler の $\partial\bar\partial$ 補題により右辺は $i\partial\bar\partial$ の実関数の像となる。得られた新形式は $e^f\omega$ と等しく、$e^f>0$ のため正である。従って曲線については非線形の事前評価を使わずに存在まで分かる。複素次元が二以上では $(i\partial\bar\partial\varphi)^2$ などの積が現れ、単純な線形逆作用素では解けない。

### 体積の写像を線形化する

$\omega_\varphi=i\sum g_{\varphi,j\bar k}\mathrm dz_j\wedge\mathrm d\bar z_k$ に対して、対数体積比の $\varphi$ 方向の変分は
$$
\left.\frac{\mathrm d}{\mathrm dt}\right|_{0}
\log\frac{\omega_{\varphi+tu}^n}{\omega_\varphi^n}
=\tr(g_\varphi^{-1}(u_{j\bar k}))=\Lambda_{\omega_\varphi}(i\partial\bar\partial u).
$$
行列式を $\det(g_\varphi)\det(I+tg_\varphi^{-1}(u_{j\bar k}))$ に分け、$t$ の一次係数が跡になることを使った。正の Hodge Laplacian の規約では関数上で後者は $-\frac12\Delta_{\mathrm d}u$ である。従って微分の核は定数、像は $\omega_\varphi^n$ に関する平均零の関数である。連続法ではこの線形化を使って解をパラメータ方向に動かす。

:::example[平坦トーラスに非定数な体積を指定する]
$\C/(\Z+i\Z)$ で $\omega=i\mathrm dz\wedge\mathrm d\bar z$、$e^f=1+\epsilon\cos(2\pi x)$、$|\epsilon|<1$ とする。平均は $1$ である。$u=-\epsilon\cos(2\pi x)/\pi^2$ とすれば $u_{z\bar z}=(u_{xx}+u_{yy})/4=\epsilon\cos(2\pi x)$ なので、新しい Kähler 形式は $e^f\omega$ となる。正値性は $|\epsilon|<1$ で保たれ、$u$ は平均零の一意な解である。
:::

### 連続法の方程式を直接解ける例

:::example[平坦トーラスの体積指定を積の形で解く]
複素次元 $2$ の正方トーラスで、座標 $z_j=x_j+iy_j$ の両実周期を $1$ とする。初期の計量行列を $I$ とし、体積比として次の積を指定する。

$$
\begin{aligned}
e^f&=\left(\frac{\cos{\left(2 \pi x_{1} \right)}}{2} + 1\right) \left(\frac{\cos{\left(4 \pi x_{2} \right)}}{3} + 1\right)
\end{aligned}
$$

各因子は正で、周期積分では余弦の平均が零なので因子ごとの平均は $1$ である。従って積の全平均も $1$ となる。ポテンシャルは各一変数の関数を足した形に取る。混合偏微分が消えるため、行列式が指定の積になるには対角成分だけを調整すればよい。

$$
\begin{aligned}
\varphi&=- \frac{\cos{\left(2 \pi x_{1} \right)}}{2 \pi^{2}} - \frac{\cos{\left(4 \pi x_{2} \right)}}{12 \pi^{2}} \\
g_\varphi&=\left[\begin{matrix}\frac{\cos{\left(2 \pi x_{1} \right)}}{2} + 1 & 0\\0 & \frac{\cos{\left(4 \pi x_{2} \right)}}{3} + 1\end{matrix}\right] \\
\det(g_\varphi)&=\left(\frac{\cos{\left(2 \pi x_{1} \right)}}{2} + 1\right) \left(\frac{\cos{\left(4 \pi x_{2} \right)}}{3} + 1\right)
\end{aligned}
$$

各方向について $u_{z_j\bar z_j}=u_{x_jx_j}/4$ を計算すれば対角成分を再現する。従って新しい複素 Hessian は正定値で、同じ Kähler 類にある。さらに対数行列式の二階微分から、新しい Ricci 形式の係数は次の対角行列になる。

$$
\begin{aligned}
R_{1\overline{1}}&=\frac{\pi^{2} \left(\sin^{2}{\left(2 \pi x_{1} \right)} + \cos^{2}{\left(2 \pi x_{1} \right)} + 2 \cos{\left(2 \pi x_{1} \right)}\right)}{\left(\cos{\left(2 \pi x_{1} \right)} + 2\right)^{2}} \\
R_{2\overline{2}}&=\frac{4 \pi^{2} \left(\sin^{2}{\left(4 \pi x_{2} \right)} + \cos^{2}{\left(4 \pi x_{2} \right)} + 3 \cos{\left(4 \pi x_{2} \right)}\right)}{\left(\cos{\left(4 \pi x_{2} \right)} + 3\right)^{2}}
\end{aligned}
$$

新しい体積が旧体積と同じ平均を持つことと、Ricci が点ごとに零であることは異なる。各 $\lambda_j$ の正値性から計量は定義されるが、一般の点で Ricci は零にならない。右辺が積の形を持つため、一般の非線形方程式を解く事前評価を使わずに構成できる。
:::

:::example[平坦トーラスの体積指定を積の形で解く]
複素次元 $3$ の正方トーラスで、座標 $z_j=x_j+iy_j$ の両実周期を $1$ とする。初期の計量行列を $I$ とし、体積比として次の積を指定する。

$$
\begin{aligned}
e^f&=\left(1 - \frac{\cos{\left(2 \pi x_{1} \right)}}{4}\right) \left(\frac{\cos{\left(2 \pi x_{2} \right)}}{5} + 1\right) \left(\frac{\cos{\left(4 \pi x_{3} \right)}}{3} + 1\right)
\end{aligned}
$$

各因子は正で、周期積分では余弦の平均が零なので因子ごとの平均は $1$ である。従って積の全平均も $1$ となる。ポテンシャルは各一変数の関数を足した形に取る。混合偏微分が消えるため、行列式が指定の積になるには対角成分だけを調整すればよい。

$$
\begin{aligned}
\varphi&=\frac{\cos{\left(2 \pi x_{1} \right)}}{4 \pi^{2}} - \frac{\cos{\left(2 \pi x_{2} \right)}}{5 \pi^{2}} - \frac{\cos{\left(4 \pi x_{3} \right)}}{12 \pi^{2}} \\
g_\varphi&=\left[\begin{matrix}1 - \frac{\cos{\left(2 \pi x_{1} \right)}}{4} & 0 & 0\\0 & \frac{\cos{\left(2 \pi x_{2} \right)}}{5} + 1 & 0\\0 & 0 & \frac{\cos{\left(4 \pi x_{3} \right)}}{3} + 1\end{matrix}\right] \\
\det(g_\varphi)&=\left(1 - \frac{\cos{\left(2 \pi x_{1} \right)}}{4}\right) \left(\frac{\cos{\left(2 \pi x_{2} \right)}}{5} + 1\right) \left(\frac{\cos{\left(4 \pi x_{3} \right)}}{3} + 1\right)
\end{aligned}
$$

各方向について $u_{z_j\bar z_j}=u_{x_jx_j}/4$ を計算すれば対角成分を再現する。従って新しい複素 Hessian は正定値で、同じ Kähler 類にある。さらに対数行列式の二階微分から、新しい Ricci 形式の係数は次の対角行列になる。

$$
\begin{aligned}
R_{1\overline{1}}&=\frac{\pi^{2} \left(\sin^{2}{\left(2 \pi x_{1} \right)} + \cos^{2}{\left(2 \pi x_{1} \right)} - 4 \cos{\left(2 \pi x_{1} \right)}\right)}{\left(\cos{\left(2 \pi x_{1} \right)} - 4\right)^{2}} \\
R_{2\overline{2}}&=\frac{\pi^{2} \left(\sin^{2}{\left(2 \pi x_{2} \right)} + \cos^{2}{\left(2 \pi x_{2} \right)} + 5 \cos{\left(2 \pi x_{2} \right)}\right)}{\left(\cos{\left(2 \pi x_{2} \right)} + 5\right)^{2}} \\
R_{3\overline{3}}&=\frac{4 \pi^{2} \left(\sin^{2}{\left(4 \pi x_{3} \right)} + \cos^{2}{\left(4 \pi x_{3} \right)} + 3 \cos{\left(4 \pi x_{3} \right)}\right)}{\left(\cos{\left(4 \pi x_{3} \right)} + 3\right)^{2}}
\end{aligned}
$$

新しい体積が旧体積と同じ平均を持つことと、Ricci が点ごとに零であることは異なる。各 $\lambda_j$ の正値性から計量は定義されるが、一般の点で Ricci は零にならない。最高次の外積では三因子が掛け合わされ、一次近似だけではこの体積を再現しない。
:::

:::example[平坦トーラスの体積指定を積の形で解く]
複素次元 $2$ の正方トーラスで、座標 $z_j=x_j+iy_j$ の両実周期を $1$ とする。初期の計量行列を $I$ とし、体積比として次の積を指定する。

$$
\begin{aligned}
e^f&=\left(1 - \frac{2 \cos{\left(6 \pi x_{2} \right)}}{5}\right) \left(\frac{2 \cos{\left(4 \pi x_{1} \right)}}{3} + 1\right)
\end{aligned}
$$

各因子は正で、周期積分では余弦の平均が零なので因子ごとの平均は $1$ である。従って積の全平均も $1$ となる。ポテンシャルは各一変数の関数を足した形に取る。混合偏微分が消えるため、行列式が指定の積になるには対角成分だけを調整すればよい。

$$
\begin{aligned}
\varphi&=- \frac{\cos{\left(4 \pi x_{1} \right)}}{6 \pi^{2}} + \frac{2 \cos{\left(6 \pi x_{2} \right)}}{45 \pi^{2}} \\
g_\varphi&=\left[\begin{matrix}\frac{2 \cos{\left(4 \pi x_{1} \right)}}{3} + 1 & 0\\0 & 1 - \frac{2 \cos{\left(6 \pi x_{2} \right)}}{5}\end{matrix}\right] \\
\det(g_\varphi)&=\left(1 - \frac{2 \cos{\left(6 \pi x_{2} \right)}}{5}\right) \left(\frac{2 \cos{\left(4 \pi x_{1} \right)}}{3} + 1\right)
\end{aligned}
$$

各方向について $u_{z_j\bar z_j}=u_{x_jx_j}/4$ を計算すれば対角成分を再現する。従って新しい複素 Hessian は正定値で、同じ Kähler 類にある。さらに対数行列式の二階微分から、新しい Ricci 形式の係数は次の対角行列になる。

$$
\begin{aligned}
R_{1\overline{1}}&=\frac{8 \pi^{2} \left(2 \sin^{2}{\left(4 \pi x_{1} \right)} + 2 \cos^{2}{\left(4 \pi x_{1} \right)} + 3 \cos{\left(4 \pi x_{1} \right)}\right)}{\left(2 \cos{\left(4 \pi x_{1} \right)} + 3\right)^{2}} \\
R_{2\overline{2}}&=\frac{18 \pi^{2} \left(2 \sin^{2}{\left(6 \pi x_{2} \right)} + 2 \cos^{2}{\left(6 \pi x_{2} \right)} - 5 \cos{\left(6 \pi x_{2} \right)}\right)}{\left(2 \cos{\left(6 \pi x_{2} \right)} - 5\right)^{2}}
\end{aligned}
$$

新しい体積が旧体積と同じ平均を持つことと、Ricci が点ごとに零であることは異なる。各 $\lambda_j$ の正値性から計量は定義されるが、一般の点で Ricci は零にならない。ここでは各方向を別々に解いた後の積が、正確に指定した体積になる。
:::

## 例：K3 曲面と五次超曲面

:::definition[K3 曲面]{#def-cg-k3}
ここでは連結コンパクト Kähler 曲面で $K_X\cong\mathcal O_X$ かつ $H^1(X,\mathcal O_X)=0$ のものを :term[K3 曲面]{en="K3 surface" yomi="けーすりーきょくめん"} と呼ぶ。単連結性と標準束の自明性による定義との同値には曲面の分類の結果を使う。詳細は :cite[huybrechts-cg] を参照する。
:::

### 四次曲面

$X=\{F_4=0\}\subset\mathbb P^3$ を滑らかな四次曲面とする。例えば Fermat 式 $\sum_{i=0}^3 Z_i^4=0$ は全偏微分が同時に零となる射影点を持たない。:ref[thm-cg-adjunction] で $K_X=\mathcal O_X(4-4)$ となり自明である。四次式は四次 Veronese 埋め込みでの超平面切断である。従って :ref[thm-cg-weak-lefschetz] から $X$ は連結かつ $H^1(X,\Z)=0$ となり、Hodge 分解で $H^1(X,\mathcal O_X)=0$ となる。また同章のホモトピー版により $\pi_1(X)\cong\pi_1(\mathbb P^3)=0$ であり、この例は単連結性による定義も満たす。

$H=c_1(\mathcal O_X(1))$ と置く。法束は $\mathcal O_X(4)$ であり、接束の完全列と射影空間の Euler 列から
$$
c(T_X)=\frac{(1+H)^4}{1+4H}=1+6H^2
$$
となる。$\int_XH^2=4$ だから $\int_Xc_2=24$ である。最高 Chern 類は実接束の Euler 類であるため Euler 標数は $24$ である。[Chern 類と Euler 類](/math/differential-topology/chern-pontryagin/) を参照する。$b_1=b_3=0$、$b_0=b_4=1$ から $b_2=22$、$h^{2,0}=h^{0,2}=1$ を引いて $h^{1,1}=20$ となる。$\Omega$ との収縮で $T_X\cong\Omega_X^1$、従って $h^1(T_X)=20$ である。Serre 双対性で $h^2(T_X)=h^0(\Omega_X^1)=0$ だから、倉西基底の次元は $20$ になる。

### 五次超曲面

$Y=\{F_5=0\}\subset\mathbb P^4$ を滑らかな五次超曲面とする。随伴公式は $K_Y=\mathcal O_Y(5-5)$ を与える。五次 Veronese 埋め込みで超平面切断として扱うと、$b_1=0,b_2=1$ である。超平面類は $(1,1)$ なので $h^{1,1}=1$、$h^{2,0}=h^{0,2}=0$ となる。一方 $K_Y$ が自明なので $h^{3,0}=1$ である。法束が $\mathcal O_Y(5)$ のため
$$
c(T_Y)=\frac{(1+H)^5}{1+5H}=1+10H^2-40H^3
$$
となる。$\int_YH^3=5$ だから Euler 標数は $-200$ である。Poincaré 双対性で $b_4=b_2=1,b_5=b_1=0,b_6=1$ なので、$-200=4-b_3$ より $b_3=204$ となる。Hodge 対称性で $b_3=2+2h^{2,1}$、従って $h^{2,1}=101$ である。

:::example[多項式の係数から見える変形]
五変数の五次斉次単項式は $\binom94=126$ 個である。方程式の定数倍を除くと射影係数空間は次元 $125$ であり、射影座標変更の群 $\mathrm{PGL}_5$ は次元 $24$ である。安定化部分群が離散なら軌道で引かれる次元は $24$ なので、係数から得る局所的な次元の差は $101$ である。Fermat 式ではその無限小安定化部分群は、微分 $\sum_{i,j}5a_{ij}Z_i^4Z_j$ が元の式の定数倍になる行列で求まる。非対角項はすべて零で、対角成分は等しくなるので射影群では無限小安定化は零である。
:::

$Y$ の非消失最高次形式 $\Omega$ との収縮は $T_Y\cong\Omega_Y^2$ を与え、$h^1(T_Y)=h^{2,1}=101$ である。障害空間は $H^2(T_Y)\cong H^{2,2}(Y)$ で一次元であるが、次の議論により障害写像は恒等的に零になる。

### 標準束が自明な場合の障害の消滅

:::theorem[Bogomolov–Tian–Todorov の定理]
連結コンパクト Kähler 多様体で標準束が正則に自明なら、倉西基底は滑らかであり、その次元は $h^1(X,T_X)$ である。
:::

:::proof
非消失正則最高次形式 $\Omega$ を固定し、収縮による同型 $I:A^{0,q}(T_X)\to A^{n-1,q}$ を取る。ここでは反正則形式を左に置き、$I(\alpha\otimes v)=\alpha\wedge\iota_v\Omega$ とする。この順序で $I$ は $\bar\partial$ と可換である。Hodge 分解により、各一次類は $\partial I(\mu_1)=0=\bar\partial I(\mu_1)$ を満たす代表を持つ。

$\partial I(\mu)=\partial I(\nu)=0$ を満たす次数一の元に対し、$I([\mu,\nu])$ は $\partial$ の像に入る。これを見るため、局所で $\mu=\sum_j\mathrm d\bar z_j\otimes v_j$、$\nu=\sum_k\mathrm d\bar z_k\otimes w_k$ と書く。各 $v_j,w_k$ は $\Omega$ に関して発散零である。実際、$\partial I(\mu)=0$ を反正則添字ごとに読むと $\partial\iota_{v_j}\Omega=0$ になる。Cartan の公式を正則方向の微分 $\partial$ で用いると、
$$
\iota_{[v_j,w_k]}\Omega=\partial(\iota_{v_j}\iota_{w_k}\Omega)
$$
である。右辺の収縮を $\mathrm d\bar z_j\wedge\mathrm d\bar z_k$ と掛けて和を取る。この和は二重収縮で大域的に定義され、反正則次数が二なので $\partial$ を反正則形式の前に移す符号は正である。これで括弧項の完全性が従う。

$H^1(T_X)$ の基底を表す調和代表の線形結合から出発し、パラメータの総次数 $m$ ごとに Maurer–Cartan 方程式を解く。次数 $m$ で解くべき式は
$$
\bar\partial I(\mu_m)=I\left(-\tfrac12\sum_{a+b=m}[\mu_a,\mu_b]\right).
$$
右辺は帰納法で $\partial$ の像に入る。また前の次数の方程式と Jacobi 恒等式により $\bar\partial$ で閉じている。従って :ref[thm-cg-ddbar] により $\partial\bar\partial\gamma$ と書ける。$I(\mu_m)=-\partial\gamma$ と取れば、方程式と次の帰納に必要な $\partial I(\mu_m)=0$ を同時に満たす。

この形式解を $D^*\mu=0$ の条件に移す。次数 $m$ で、正則接束に値を持つ滑らかなベクトル場 $\xi_m$ をパラメータの $m$ 次斉次式として取り、その形式的な流れで座標を変える。代数的には、複素化した滑らかな関数環に導分 $\xi_m$ の指数を作用させ、その逆で導分を共役する。各有限次数では指数が有限和であり、共役は Lie 括弧を保つので可積分性を保つ。また $[\partial_{\bar j},\xi_m]=\partial_{\bar j}\xi_m$ である。従って向きを選べば、$m$ 次以下では次数 $m$ の成分にだけ $D\xi_m$ が加わる。$D^*D\xi_m=-D^*\mu_m$ を解けばよい。右辺は $\Ker D$ と直交し、零次の Green 作用素 $G$ を適用すれば解が得られる。これを各次数で行うと、一次のコホモロジー類を保ったまま形式的に $D^*\mu=0$ となる。

この条件と Maurer–Cartan 方程式から、$\mu=P\mu-D^*G[\mu,\mu]/2$ となる。これは倉西の構成で得た解析的な固定点方程式と同じであり、$t=P\mu$ を指定すれば、次数ごとに係数が一意に決まる。変数を $t=P\mu$ に替える操作の一次係数は可逆である。最初に全一次類を基底から選び、補正では厳密形式だけを加えたためである。よってこの変数変換は形式冪級数環上で逆を持ち、残差の調和成分 $k(t)$ の全 Taylor 係数は零になる。倉西の方程式 $k$ 自体は正則なので、実際に近傍で恒等的に零である。従って基底は滑らかである。
:::

### 曲面と三次元の方程式による変形を比較する

五次超曲面 $Y$ では $H^1(T_{\mathbb P^4}|_Y)=0$ である。実際、制限した Euler 列でこの群の前に $H^1(\mathcal O_Y(1))^{\oplus5}$、後ろに $H^2(\mathcal O_Y)$ がある。前者は $K_Y$ が自明で $\mathcal O_Y(1)$ が正であるため Kodaira 消滅で零であり、後者は $h^{0,2}=0$ により零である。従って法束の列で $H^0(\mathcal O_Y(5))\to H^1(T_Y)$ は全射となる。さらに $\mathcal O_Y(5)$ の切断は五次斉次多項式の制限で尽くされる。これは $0\to\mathcal O_{\mathbb P^4}\xrightarrow{F_5}\mathcal O_{\mathbb P^4}(5)\to\mathcal O_Y(5)\to0$ と $H^1(\mathbb P^4,\mathcal O)=0$ から従う。従って全ての一次変形が方程式の係数から得られる。核は座標変更の像であり、$h^0(T_Y)=h^{2,0}=0$ と Euler 列からその次元は $24$ である。ここで $H^0(Y,\mathcal O_Y(1))$ は $\mathbb P^4$ の線形形式の制限で、次元は $5$ である。

一方、四次 K3 曲面 $X$ では Euler 列と Kodaira 消滅から $H^1(T_{\mathbb P^3}|_X)\cong H^2(\mathcal O_X)\cong\C$ となる。$H^1(\mathcal O_X(4))=0$ なので法束の列から $H^1(T_X)\to H^1(T_{\mathbb P^3}|_X)$ は全射である。従って方程式の係数から来る一次変形は $H^1(T_X)$ の余次元一の部分で、次元は $19$ である。四変数の四次単項式 $\binom73=35$ 個から定数倍を一つと射影座標変更の次元 $15$ を引いても $19$ になる。従って全ての K3 の一次変形が、固定した射影空間の四次方程式を動かすことで得られるわけではない。

### 四次・五次の係数を Jacobian 環で数える

Fermat 超曲面 $F=\sum_{j=0}^{N}Z_j^d$ では、斉次 $d$ 次の変形を射影座標変更と方程式の定数倍で割る計算が明示的になる。座標変更の微分は $Z_i^{d-1}Z_j$ の線形結合であり、元の $F$ 自身もその線形結合である。従って商の基底は各指数が高々 $d-2$ の $d$ 次単項式である。$d\ge3$ では次数 $d$ の単項式が異なる二添字で指数 $d-1$ 以上を持てないので、除く単項式の数は重複なしに $(N+1)^2$ 個となる。

:::example[四次と五次の違いを基底から読む]
四変数四次では $\binom73-16=19$、五変数五次では $\binom94-25=101$ となる。四次の $Z_0^2Z_1^2$ は指数がいずれも $2$ 以下なので残るが、$Z_0^3Z_1$ は座標変更で消える。五次の $Z_0Z_1Z_2Z_3Z_4$ も残り、すべての変数を同時に動かす方向の一例になる。四次で抽象変形が一つ多いのは、射影空間内の変形が超平面束の類を保つことに対応する。
:::

:::example[全 Chern 類から別の hypersurface を検算する]
$\mathbb P^5$ の滑らかな六次超曲面では随伴公式から $K_X$ は自明である。$\dim_\C X=4$ なので
$$
\frac{(1+H)^6}{1+6H}=1+15H^2-70H^3+435H^4
$$
まで取る。$\int H^4=6$ より Euler 標数は $2610$ である。曲面や三次元と同じ Chern 類の方法が通用するが、Euler 標数だけでは中間次数の各 Hodge 数まで決まらない。
:::

### Chern 数と変形例を比較する

:::example[完全交差の標準束と Chern 数を計算する]
$X\subset\mathbb P^{4}$ を次数 $2,3$ の横断的な滑らかな完全交差とし、$H=c_1(\mathcal O_X(1))$ とする。法束は $\bigoplus_a\mathcal O_X(d_a)$ であり、次元は $2$ である。

接束の完全列から $c(T_X)c(N_{X/\mathbb P^N})=c(T_{\mathbb P^N}|_X)$ を得る。Euler 列により右辺は $(1+H)^{N+1}$ である。$c(T_X)=1+c_1H+\cdots+c_nH^n$ と書いて係数を比べる。

$$
\begin{aligned}
c(N)&=\left(2 H + 1\right) \left(3 H + 1\right) \\
c_{1}+(5)&=\binom{5}{1},\qquad c_{1}=0 \\
c_{2}+(6)&=\binom{5}{2},\qquad c_{2}=4 \\
c(T_X)&=4 H^{2} + 1
\end{aligned}
$$

法束の行列式を取った随伴公式は $K_X=\mathcal O_X(0)$ である。これは $c(T_X)$ の一次の係数の符号を反転した値にも等しい。さらに $X$ を一般の $2$ 個の超平面で切り、交点数を積で数えると $\int_XH^2=6$ となる。

$$
\begin{aligned}
\int_X c_n(T_X)&=4\cdot6=24
\end{aligned}
$$

最高 Chern 類は実接束の Euler 類であるため、ここで得た整数は位相的 Euler 標数にも等しい。標準束が自明な曲面の Euler 標数は $24$ であり、四次曲面と一致する。
:::

:::example[完全交差の標準束と Chern 数を計算する]
$X\subset\mathbb P^{5}$ を次数 $2,2,2$ の横断的な滑らかな完全交差とし、$H=c_1(\mathcal O_X(1))$ とする。法束は $\bigoplus_a\mathcal O_X(d_a)$ であり、次元は $2$ である。

接束の完全列から $c(T_X)c(N_{X/\mathbb P^N})=c(T_{\mathbb P^N}|_X)$ を得る。Euler 列により右辺は $(1+H)^{N+1}$ である。$c(T_X)=1+c_1H+\cdots+c_nH^n$ と書いて係数を比べる。

$$
\begin{aligned}
c(N)&=\left(2 H + 1\right)^{3} \\
c_{1}+(6)&=\binom{6}{1},\qquad c_{1}=0 \\
c_{2}+(12)&=\binom{6}{2},\qquad c_{2}=3 \\
c(T_X)&=3 H^{2} + 1
\end{aligned}
$$

法束の行列式を取った随伴公式は $K_X=\mathcal O_X(0)$ である。これは $c(T_X)$ の一次の係数の符号を反転した値にも等しい。さらに $X$ を一般の $2$ 個の超平面で切り、交点数を積で数えると $\int_XH^2=8$ となる。

$$
\begin{aligned}
\int_X c_n(T_X)&=3\cdot8=24
\end{aligned}
$$

最高 Chern 類は実接束の Euler 類であるため、ここで得た整数は位相的 Euler 標数にも等しい。三つの二次曲面で切った曲面でも同じ Euler 標数になる。
:::

:::example[完全交差の標準束と Chern 数を計算する]
$X\subset\mathbb P^{6}$ を次数 $2,2,3$ の横断的な滑らかな完全交差とし、$H=c_1(\mathcal O_X(1))$ とする。法束は $\bigoplus_a\mathcal O_X(d_a)$ であり、次元は $3$ である。

接束の完全列から $c(T_X)c(N_{X/\mathbb P^N})=c(T_{\mathbb P^N}|_X)$ を得る。Euler 列により右辺は $(1+H)^{N+1}$ である。$c(T_X)=1+c_1H+\cdots+c_nH^n$ と書いて係数を比べる。

$$
\begin{aligned}
c(N)&=\left(2 H + 1\right)^{2} \left(3 H + 1\right) \\
c_{1}+(7)&=\binom{7}{1},\qquad c_{1}=0 \\
c_{2}+(16)&=\binom{7}{2},\qquad c_{2}=5 \\
c_{3}+(47)&=\binom{7}{3},\qquad c_{3}=-12 \\
c(T_X)&=- 12 H^{3} + 5 H^{2} + 1
\end{aligned}
$$

法束の行列式を取った随伴公式は $K_X=\mathcal O_X(0)$ である。これは $c(T_X)$ の一次の係数の符号を反転した値にも等しい。さらに $X$ を一般の $3$ 個の超平面で切り、交点数を積で数えると $\int_XH^3=12$ となる。

$$
\begin{aligned}
\int_X c_n(T_X)&=-12\cdot12=-144
\end{aligned}
$$

最高 Chern 類は実接束の Euler 類であるため、ここで得た整数は位相的 Euler 標数にも等しい。こちらは三次元であり、標準束が自明でも Euler 標数は曲面の値と異なる。
:::

:::example[完全交差の標準束と Chern 数を計算する]
$X\subset\mathbb P^{7}$ を次数 $2,2,2,2$ の横断的な滑らかな完全交差とし、$H=c_1(\mathcal O_X(1))$ とする。法束は $\bigoplus_a\mathcal O_X(d_a)$ であり、次元は $3$ である。

接束の完全列から $c(T_X)c(N_{X/\mathbb P^N})=c(T_{\mathbb P^N}|_X)$ を得る。Euler 列により右辺は $(1+H)^{N+1}$ である。$c(T_X)=1+c_1H+\cdots+c_nH^n$ と書いて係数を比べる。

$$
\begin{aligned}
c(N)&=\left(2 H + 1\right)^{4} \\
c_{1}+(8)&=\binom{8}{1},\qquad c_{1}=0 \\
c_{2}+(24)&=\binom{8}{2},\qquad c_{2}=4 \\
c_{3}+(64)&=\binom{8}{3},\qquad c_{3}=-8 \\
c(T_X)&=- 8 H^{3} + 4 H^{2} + 1
\end{aligned}
$$

法束の行列式を取った随伴公式は $K_X=\mathcal O_X(0)$ である。これは $c(T_X)$ の一次の係数の符号を反転した値にも等しい。さらに $X$ を一般の $3$ 個の超平面で切り、交点数を積で数えると $\int_XH^3=16$ となる。

$$
\begin{aligned}
\int_X c_n(T_X)&=-8\cdot16=-128
\end{aligned}
$$

最高 Chern 類は実接束の Euler 類であるため、ここで得た整数は位相的 Euler 標数にも等しい。四つの正の法直線束の積から三次の Chern 類を計算する。
:::

### 四次曲面の偏極条件を一次で読む

標準束の非消失切断 $\Omega$ で収縮すると、一次変形の空間 $H^1(T_X)$ は $H^{1,1}(X)$ と同一視される。超平面類 $H$ が変形後も $(1,1)$ 型に留まる条件の一次は、対応するクラスが $H$ と直交することである。実際、$H$ が $(2,0)$ 形式の変形 $\Omega_t$ と直交するためには $\int_XH\wedge\Omega_t=0$ でなければならない。この積分を微分すると、$\Omega_t$ の変分のうち $(1,1)$ 成分だけが積分へ寄与する。新しい $(1,0)$ 形式 $\mathrm dz_k-\sum_j\mu_j^k\mathrm d\bar z_j$ を最高外積に代入し、$\mu$ の一次項を取る。反正則形式を左に置く収縮の規約では、変分の $(1,1)$ 成分は $-\mu\lrcorner\Omega$ になる。ここで $\mu$ は変形方向の Beltrami 微分の一次係数を表す。

$H$ 自身との積分は $\int_XH^2=4$ で非零なので、直交条件は一次方程式一本である。そして $h^{1,1}=20$ であるため直交部分の次元は $19$ となる。射影空間内で四次方程式を変えても、固定した $\mathcal O(1)$ の制限として超平面束は保たれる。従ってその一次方向はこの直交部分に入る。先ほど法束の完全列で求めた次元も $19$ なので、二つの部分空間は一致する。

## 演習問題

::::exercise{level=basic}
種数 $4$ の Riemann 面の倉西基底の次元を求めよ。
:::solution
曲線の接束の $H^2$ は零なので基底は滑らかである。Serre 双対性で $h^1(T_X)=h^0(K_X^2)$、次数 $4g-4=12$ は $2g-2=6$ より大きい。Riemann–Roch で次元は $12+1-4=9$ である。
:::
::::

::::exercise{level=basic}
複素次元 $n$ の射影空間内の滑らかな次数 $n+1$ の超曲面の標準束を求めよ。
:::solution
標準束は随伴公式で $\mathcal O_X((n+1)-(n+1))=\mathcal O_X$ である。従って大域で消えない正則最高次形式をもつ。
:::
::::

::::exercise{level=basic}
五次超曲面の $c_2(T_Y)$ の $H^2$ の係数を求めよ。
:::solution
$(1+H)^5=1+5H+10H^2+\cdots$、$(1+5H)^{-1}=1-5H+25H^2+\cdots$ である。二次係数は $10-25+25=10$ なので $c_2(T_Y)=10H^2$ である。
:::
::::

::::exercise{level=basic}
$g=1$ の曲線で $H^1(X,T_X)$ の次元を求めよ。
:::solution
$K_X$ が自明なので $T_X$ も自明である。従って $h^1(T_X)=h^1(\mathcal O_X)=g=1$ である。格子比の微小変化がその非零な方向を与える。
:::
::::

::::exercise{level=advanced}
$\mu=t\mu_1+t^2\mu_2+O(t^3)$ が可積分性方程式を満たすとき、二次係数に必要な条件を書け。
:::solution
$\bar\partial$ は線形で、括弧の最低の積は $t^2[\mu_1,\mu_1]$ である。従って一次係数から $\bar\partial\mu_1=0$、二次係数から $\bar\partial\mu_2=-[\mu_1,\mu_1]/2$ を得る。後者を解くには右辺のコホモロジー類が零でなければならない。逆に右辺の類が零なら定義により二次までは解ける。
:::
::::

::::exercise{level=advanced}
既に Kähler なコンパクト $X$ で $c_1(T_X)=0$ が実コホモロジーで成り立つとき、Ricci 平坦計量が各 Kähler 類に一意にあることを導け。
:::solution
Ricci 形式の類は $2\pi c_1(T_X)$ なので、目標の実 $(1,1)$ 形式に $0$ を取れる。$\partial\bar\partial$ 補題で元の Ricci 形式を $i\partial\bar\partial f$ と表し、$f$ の定数を体積条件に合わせる。Yau の定理が同じ類の計量を与え、Ricci 形式との差の公式からその計量は平坦な Ricci 形式を持つ。一意性も同じ定理に従う。
:::
::::

::::exercise{level=advanced}
四次 K3 曲面の $h^{1,1}$ を最高 Chern 類から求めよ。
:::solution
$(1+H)^4/(1+4H)$ の二次係数は $6-16+16=6$ である。$\int H^2=4$ なので Euler 標数は $24$、$b_1=b_3=0$ だから $b_2=22$ となる。標準束が自明で連結なので $h^{2,0}=h^{0,2}=1$ であり、差を取って $h^{1,1}=20$ を得る。
:::
::::

::::exercise{level=advanced}
五次超曲面の $b_3$ と $h^{2,1}$ を求めよ。
:::solution
Euler 標数は $-40\int H^3=-40\cdot5=-200$ である。超平面定理と Poincaré 双対性で $b_0,b_2,b_4,b_6$ は $1$、$b_1,b_5$ は $0$ である。従って $-200=4-b_3$ で $b_3=204$ である。$h^{3,0}=h^{0,3}=1$、$h^{2,1}=h^{1,2}$ より $h^{2,1}=(204-2)/2=101$ となる。
:::
::::

::::exercise{level=advanced}
四次の Fermat 曲面 $X$ で $Z_0^2Z_1^2$ と $Z_0^3Z_1$ を方程式に一次に加える。どちらが零でない埋め込み変形を与えるか。さらに係数空間から抽象変形への写像の余次元を求めよ。

:::solution
$F=\sum_{j=0}^3Z_j^4$ に対し、無限小線形座標変更は $\sum_{i,j}4a_{ij}Z_i^3Z_j$ を与える。これら十六個の単項式は互いに異なり、その線形結合の中に $F$ 自身も含まれる。このため方程式を定数倍する自由度もこの部分空間に含められる。$Z_0^3Z_1$ は $a_{01}=1/4$ によって作れるので零類である。$Z_0^2Z_1^2$ には三以上の指数がなく、十六個の単項式の線形結合では表せない。従ってこの商で非零である。

確かにこの商が抽象変形の一部として単射になることを確認する。制限した Euler 列より $H^0(T_{\mathbb P^3}|_X)$ は一次切断四組から定数方向を除いた空間である。$H^1(\mathcal O_X)=0$ が全射を保証する。$H^0(\mathcal O_X(4))$ は四次式を $F$ の倍数で割った空間であり、法束の列の連結写像の核は $H^0(T_{\mathbb P^3}|_X)$ の像である。従って先ほどの商の次元 $\binom73-16=19$ がそのまま $H^1(T_X)$ の中での次元になる。後者は $\Omega$ による収縮と Hodge 数の計算から次元 $20$ なので、像の余次元は $1$ である。
:::
::::

## 確認テスト

::::quiz{answer="101" tol="1e-6"}
滑らかな五次三次元超曲面の $h^{2,1}$ を求めよ。
:::explanation
$c_3=-40H^3$ で $\int H^3=5$ なので Euler 標数は $-200$ である。低次数と双対性から $b_3=204$ であり、$h^{3,0}=h^{0,3}=1$ を引いて残りを二等分すると $101$ になる。
:::
::::

::::quiz
Kodaira–Spencer 写像の値域を選べ。

- [ ] $H^0(X,K_X)$
- [x] $H^1(X,T_X)$
- [ ] $H^2(X,\Z)$

:::explanation
遷移のパラメータ微分は重なり上の正則接ベクトル場を与える。三重交差の条件は $1$ 次コサイクル条件であり、座標変更はコバウンダリとなるため $H^1(X,T_X)$ に値を持つ。
:::
::::

::::quiz{answer="24" tol="1e-6"}
滑らかな四次 K3 曲面の位相的 Euler 標数を求めよ。
:::explanation
接束の全 Chern 類は $(1+H)^4/(1+4H)$ なので最高 Chern 類は $6H^2$ である。$\int H^2=4$ より Euler 標数は $24$ になる。
:::
::::

::::quiz
Yau の定理で方程式 $(\omega+i\partial\bar\partial\varphi)^n=e^f\omega^n$ に必要な規格化を選べ。

- [x] $\int_Xe^f\omega^n=\int_X\omega^n$
- [ ] $\int_X f\omega^n=0$
- [ ] $f$ を負の定数とする。

:::explanation
左辺の計量は $\omega$ と同じコホモロジー類なので体積が等しい。従って積分すると最初の条件が必要であり、Yau の定理はこの条件の下での存在と一意性を保証する。
:::
::::

::::quiz{answer="19" tol="1e-6"}
四次 K3 曲面で超平面束の類を保つ一次変形の次元を求めよ。

:::explanation
収縮により $H^1(T_X)$ は次元 $20$ の $H^{1,1}(X)$ である。超平面類との交点が零であることを課す。超平面類自身の二乗は $4\ne0$ なのでこれは独立な方程式一本であり、次元は $20-1=19$ となる。
:::
::::

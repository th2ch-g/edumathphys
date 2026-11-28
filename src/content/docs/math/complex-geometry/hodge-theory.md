---
title: コンパクト Kähler 多様体の Hodge 理論
description: 調和代表から Hodge 分解、強 Lefschetz 定理と交点形式の符号を導く
notation:
  - symbol: '\mathcal H^{p,q}(X)'
    meaning: 型 (p,q) の調和形式の空間
  - symbol: 'h^{p,q}'
    meaning: 型 (p,q) の Dolbeault コホモロジーの次元
  - symbol: 'P^k'
    meaning: k 次原始コホモロジー
  - symbol: 'G'
    meaning: 調和形式の直交補上での Laplacian の逆
---

## この章の目標

- 楕円型 Laplacian の調和代表を型に分ける。
- $\partial\bar\partial$ 補題から、完全形式を二重微分で表す。
- Lefschetz 作用をコホモロジーに移して次元の制約を得る。
- 実 $(1,1)$ 類の交点形式の符号と、有理 Hodge 類の問題を説明する。

以下、$X$ は複素次元 $n$ のコンパクト連結 Kähler 多様体とする。

## Hodge 分解定理

閉形式の類に最も扱いやすい代表を選ぶには、閉性に加えて随伴でも消えることを課す。$D=\dd,\bar\partial,\partial$ に対し $\Delta_D=DD^*+D^*D$ と書く。部分積分で
$$
\inner{\Delta_D\alpha}{\alpha}=\|D\alpha\|^2+\|D^*\alpha\|^2
$$
なので、$\Ker\Delta_D=\Ker D\cap\Ker D^*$ である。$\Delta_D$ は楕円型で自己随伴であり、コンパクト性から滑らかな有限次元の核と、その直交補上での逆を持つ。ここで使う楕円型 Fredholm 理論と正則性は [Hodge 理論](/math/riemannian-geometry/hodge-bochner/) と [楕円型正則性](/math/pde/elliptic-regularity/) を参照する。

:::theorem[調和代表と Hodge 分解]{#thm-cg-hodge-decomposition}
各 de Rham 類は一意な調和代表を持ち、
$$
H^k_{\mathrm{dR}}(X,\C)=
\bigoplus_{p+q=k} H^{p,q}(X),\qquad
H^{p,q}(X)\cong H^q(X,\Omega_X^p)
$$
である。$H^{p,q}(X)$ は型 $(p,q)$ の調和形式が表す de Rham 類である。
:::

:::proof
$P$ を $\Ker\Delta_D$ への直交射影とし、$\Delta_D$ の逆を核上零へ延長して Green 作用素 $G_D$ とする。$\alpha=P\alpha+DD^*G_D\alpha+D^*DG_D\alpha$ である。$D^2=0$ から $D$ は $\Delta_D$ と可換である。また調和成分は $D$ で消え、$D\beta$ は調和形式と直交する。従って $\Delta_D$ を作用させて比較し、核の直交補上での一意性を使えば $DG_D=G_DD$ である。$D\alpha=0$ なら最後の項は零なので $P\alpha$ が同じコホモロジー類を表す。調和で $D$ の像に入る元は、$D^*$ で消えるため自分自身と直交し零となる。このため代表は一意である。

Kähler 恒等式 :ref[thm-cg-kahler-identities] から $\Delta_\dd=2\Delta_{\bar\partial}$ で、演算子は型を保つ。よって調和形式も型に分けられ、各型の調和空間は $\bar\partial$ コホモロジーと同一視できる。最後に :ref[thm-cg-global-dolbeault] を適用する。
:::

:::example[平坦トーラスの調和形式]
$\C^n/\Lambda$ では係数が定数の枠 $\dd z_j,\dd\bar z_j$ を使える。Laplacian は係数ごとの平坦 Laplacian なので、調和形式の係数はコンパクト面上の調和関数として定数になる。任意の定数係数形式は調和であり、従って $h^{p,q}=\binom np\binom nq$ である。
:::

:::example[射影平面の二次類]
$\mathbb P^2$ の $b_2$ は CW 分解で $1$ である。$\omega_{\mathrm{FS}}$ は閉じており、$\omega_{\mathrm{FS}}^2$ の積分が正だから $[\omega_{\mathrm{FS}}]\ne0$ である。これは $(1,1)$ 型の閉形式であり、調和射影もその型を保つ。従って $h^{1,1}=1,h^{2,0}=h^{0,2}=0$ となる。
:::

### 調和代表がノルムを最小にする

$\alpha$ を閉形式、$P\alpha$ を調和代表とする。同じ類の別の代表は $P\alpha+\mathrm d\beta$ であり、$\mathrm d^*P\alpha=0$ から両項は直交する。従ってノルムの二乗は $\|P\alpha\|^2+\|\mathrm d\beta\|^2$ であり、最小値を実現する代表は調和代表ただ一つである。$D=\bar\partial_E$ にも同じ議論が使える。これにより差の完全形式をどのように選んでも、調和部分だけは変わらない理由が分かる。

:::example[トーラスで定数項だけ残る]
$X=\C/(\Z+i\Z)$ で $\alpha=3\mathrm dx-2\mathrm dy+\mathrm d(\sin2\pi x)$ とする。最後の項は完全であり、前二項は定数係数なので調和である。水平、垂直の周期は $3,-2$ である。Hermite 規格化 $g=2(\mathrm dx^2+\mathrm dy^2)$ では、定数部分の $L^2$ ノルムの二乗は $3^2+(-2)^2=13$ となり、完全部分のノルムの二乗は $2\pi^2$ である。
:::

### 調和代表の定数係数を抽出する

:::example[平均と非零 Fourier モードを分離する]
$X=\C/(\Z+5i\Z)$ 上で $\alpha=a\,\mathrm d\bar z$ とし、周期座標を $z=x+iy$ とする。係数と微分の規約は次である。

$$
\begin{aligned}
a(x,y)&=2 e^{2 i \pi \left(x + \frac{2 y}{5}\right)} - 3 + e^{- \frac{2 i \pi y}{5}} \\
\partial_{\bar z}&=\tfrac12(\partial_x+i\partial_y)
\end{aligned}
$$

周期積分では全微分の平均が零になる。したがって定数モードは $\bar\partial$ の方程式の障害である。他方、指数関数 $e^{2\pi i(mx+ny/L)}$ の二つの微分で $2\pi im,2\pi in/L$ が掛かるので、固有値は $\pi(im-n/L)$ になる。

$$
\begin{aligned}
(1,2) &\colon\quad \lambda=\pi \left(- \frac{2}{5} + i\right),\quad a=2,\quad a/\lambda=\frac{10 \left(-2 - 5 i\right)}{29 \pi} \\
(0,-1) &\colon\quad \lambda=\frac{\pi}{5},\quad a=1,\quad a/\lambda=\frac{5}{\pi} \\
a_{00}&=-3
\end{aligned}
$$

定数モードを除いて各係数を固有値で割る。これで得られる関数は平均零で、微分を取り直すと次を満たす。

$$
\begin{aligned}
u&=\frac{50 \left(- \frac{2}{5} - i\right) e^{2 i \pi \left(x + \frac{2 y}{5}\right)}}{29 \pi} + \frac{5 e^{- \frac{2 i \pi y}{5}}}{\pi} \\
u_{\bar z}&=a--3
\end{aligned}
$$

全てのモードの周期が元のトーラスと一致するため $u$ は大域的である。計量が平坦なので調和な $(0,1)$ 形式は定数係数であり、この場合は $-3\,\mathrm d\bar z$ である。
:::

:::example[平均と非零 Fourier モードを分離する]
$X=\C/(\Z+2i\Z)$ 上で $\alpha=a\,\mathrm d\bar z$ とし、周期座標を $z=x+iy$ とする。係数と微分の規約は次である。

$$
\begin{aligned}
a(x,y)&=- 2 e^{2 i \pi \left(- x + \frac{y}{2}\right)} + e^{2 i \pi \left(x - \frac{3 y}{2}\right)} + 8 \\
\partial_{\bar z}&=\tfrac12(\partial_x+i\partial_y)
\end{aligned}
$$

周期積分では全微分の平均が零になる。したがって定数モードは $\bar\partial$ の方程式の障害である。他方、指数関数 $e^{2\pi i(mx+ny/L)}$ の二つの微分で $2\pi im,2\pi in/L$ が掛かるので、固有値は $\pi(im-n/L)$ になる。

$$
\begin{aligned}
(1,-3) &\colon\quad \lambda=\pi \left(\frac{3}{2} + i\right),\quad a=1,\quad a/\lambda=\frac{2 \left(3 - 2 i\right)}{13 \pi} \\
(-1,1) &\colon\quad \lambda=\pi \left(- \frac{1}{2} - i\right),\quad a=-2,\quad a/\lambda=\frac{4 \left(1 - 2 i\right)}{5 \pi} \\
a_{00}&=8
\end{aligned}
$$

定数モードを除いて各係数を固有値で割る。これで得られる関数は平均零で、微分を取り直すと次を満たす。

$$
\begin{aligned}
u&=- \frac{8 \left(- \frac{1}{2} + i\right) e^{2 i \pi \left(- x + \frac{y}{2}\right)}}{5 \pi} + \frac{4 \left(\frac{3}{2} - i\right) e^{2 i \pi \left(x - \frac{3 y}{2}\right)}}{13 \pi} \\
u_{\bar z}&=a-8
\end{aligned}
$$

全てのモードの周期が元のトーラスと一致するため $u$ は大域的である。従って Hodge 射影が残す係数は八である。
:::

### 積の各型の調和代表

:::example[積の直線束でコホモロジーが残る次数]
$X=\prod_j\mathbb P^{n_j}$ とし、外部テンソル積の直線束 $L=\mathcal O(-4,3)$ を考える。$a$ のべきでコホモロジー次数を記録し、係数を各空間の次元とする。

一因子では符号がすべて揃った Taylor–Laurent 指数だけが残る。指数が非負の場合は全て $0$ 次に、全て負の場合は因子の最高次に属する。負指数を $-b_i-1$ と書けば、弱い合成の個数を二項係数で求められる。

$$
\begin{aligned}
j=1&\colon\quad \mathbb P^{2},\; d=-4,\quad e_i<0,\quad \sum_i e_i=-4,\quad P_j(a)=3 a^{2} \\
j=2&\colon\quad \mathbb P^{1},\; d=3,\quad e_i\ge0,\quad \sum_i e_i=3,\quad P_j(a)=4
\end{aligned}
$$

積では各因子の Dolbeault Laplacian の非負固有値を加える。零になるには各因子の固有値が零でなければならず、調和代表の次元と次数は多項式の積で合成される。

$$
\begin{aligned}
P_X(a)&=\prod_j P_j(a)=12 a^{2} \\
\dim_\C X&=3 \\
\chi(X,L)&=P_X(-1)=12
\end{aligned}
$$

随伴の双対束 $K_X\otimes L^*$ は $\mathcal O(1,-5)$ である。$d\ge0$ の因子では双対側が最高次になり、$d\le-n-1$ の因子では逆になる。全消滅の範囲はそのままなので、双対は $q$ 次と $3-q$ 次を対応させる。ここでは係数束が非自明なので、調和代表の次数を自明束の Hodge 数から直接読むことはできない。
:::

:::example[積の直線束でコホモロジーが残る次数]
$X=\prod_j\mathbb P^{n_j}$ とし、外部テンソル積の直線束 $L=\mathcal O(3,-4,1)$ を考える。$a$ のべきでコホモロジー次数を記録し、係数を各空間の次元とする。

一因子では符号がすべて揃った Taylor–Laurent 指数だけが残る。指数が非負の場合は全て $0$ 次に、全て負の場合は因子の最高次に属する。負指数を $-b_i-1$ と書けば、弱い合成の個数を二項係数で求められる。

$$
\begin{aligned}
j=1&\colon\quad \mathbb P^{1},\; d=3,\quad e_i\ge0,\quad \sum_i e_i=3,\quad P_j(a)=4 \\
j=2&\colon\quad \mathbb P^{2},\; d=-4,\quad e_i<0,\quad \sum_i e_i=-4,\quad P_j(a)=3 a^{2} \\
j=3&\colon\quad \mathbb P^{1},\; d=1,\quad e_i\ge0,\quad \sum_i e_i=1,\quad P_j(a)=2
\end{aligned}
$$

積では各因子の Dolbeault Laplacian の非負固有値を加える。零になるには各因子の固有値が零でなければならず、調和代表の次元と次数は多項式の積で合成される。

$$
\begin{aligned}
P_X(a)&=\prod_j P_j(a)=24 a^{2} \\
\dim_\C X&=4 \\
\chi(X,L)&=P_X(-1)=24
\end{aligned}
$$

随伴の双対束 $K_X\otimes L^*$ は $\mathcal O(-5,1,-3)$ である。$d\ge0$ の因子では双対側が最高次になり、$d\le-n-1$ の因子では逆になる。全消滅の範囲はそのままなので、双対は $q$ 次と $4-q$ 次を対応させる。因子の調和形式の外積を使うと、全体のコホモロジーの型が決まる。
:::

### 調和分解を使うときの計量依存性

$X$ の同じ複素構造に二つの Kähler 計量を入れたとする。de Rham の調和代表は計量によって変わるが、型 $(p,q)$ の類が張る部分は同じである。実際、一方の計量で得た調和代表は純型でかつ $\mathrm d$ 閉である。他方の計量の調和射影も型を保ち、$\mathrm d$ 閉な代表を同じ de Rham 類の調和形式に移すため、同じ類は他方でも純型である。逆向きも同じ議論である。従って Hodge 分解は調和代表で作れるが、定義された部分空間は選んだ Kähler 計量に依存しない。

:::example[面積の密度が変わると調和代表も変わる]
連結コンパクト曲線で体積形式 $\omega$ を取る。任意の実二次形式は $f\omega$ と書け、その星は $f$ である。調和である条件は $\mathrm df=0$ なので、積分 $1$ の類の調和代表は $\omega/\int_X\omega$ である。別の計量の面積形式が $\omega'=g\omega$、$g$ が正の非定数なら、新しい代表は $g\omega/\int_Xg\omega$ であり点ごとの値は変わる。一方どちらも最高次で積分 $1$ の類を表すので、同じ実 $(1,1)$ 類である。
:::

## ∂∂̄ 補題

:::theorem[$\partial\bar\partial$ 補題]{#thm-cg-ddbar}
純型の形式 $\alpha$ が $\partial,\bar\partial$ の両方で閉じているとする。$\alpha$ が $\dd,\partial,\bar\partial$ のいずれかの像に入るなら、$\alpha=\partial\bar\partial\beta$ と書ける。逆にこの表示は三つの像全てに入ることを意味する。
:::

:::proof
三つの Laplacian の調和空間は同じであるため、いずれかの像に入る $\alpha$ はその核に直交する。$G=G_{\bar\partial}$ とする。Kähler 恒等式と $\partial^2=\bar\partial^2=0$ により $\partial,\bar\partial$ はこの Laplacian と可換で、従って $G$ とも可換である。調和成分が零で $\bar\partial G\alpha=0$ なので
$$
\alpha=\bar\partial\bar\partial^*G\alpha
=-i\bar\partial[\Lambda,\partial]G\alpha
=-i\partial\bar\partial\Lambda G\alpha .
$$
最後では $\partial G\alpha=0$ と $\bar\partial\partial=-\partial\bar\partial$ を使った。逆向きはこれらの交換関係と $\dd=\partial+\bar\partial$ から成り立つ。例えば $\partial\bar\partial\beta=\dd(\bar\partial\beta)$ である。
:::

:::example[同じ Kähler 類の二形式]
二つの実 Kähler 形式の差 $\alpha$ が完全なら $\alpha=\partial\bar\partial f$ である。実性から $\alpha=-\partial\bar\partial\bar f$ でもあるため、平均を取ると $\alpha=i\partial\bar\partial\Im f$ となる。従って同じ類の計量は大域実ポテンシャルによって比較できる。
:::

:::example[型の境界にある場合]
完全な正則 $1$ 形式 $\alpha$ は $(1,0)$ 型で両微分に閉じている。補題の原始形式の型は $(0,-1)$ となってその空間が零なので $\alpha=0$ となる。実際、正則微分は調和で、完全な調和類は零である。
:::

### 同じ型の三つのコホモロジー

純型の形式で $\partial,\bar\partial$ の両方に閉じたものを $\partial\bar\partial$ の像で割った空間を :term[Bott–Chern コホモロジー]{en=Bott–Chern cohomology yomi=ぼっとちぇるんこほもろじー} と呼ぶ。各類を Dolbeault 類に送る写像がある。コンパクト Kähler の場合、任意の Dolbeault 類の調和代表は $\partial$ にも閉じているのでこの写像は全射である。核に入る代表は $\bar\partial$ の像であるため :ref[thm-cg-ddbar] で $\partial\bar\partial$ の像になり、単射でもある。同様に、閉代表の de Rham 類へ送る写像は同型 $H^{p,q}_{\mathrm{BC}}\cong H^{p,q}(X)$ を与える。このとき de Rham 側の終域は全次数の全空間でなく、同じ型の部分である。

:::example[コンパクト性を外した場合]
$\C$ 上の $\dd z$ は $\partial$ と $\bar\partial$ でともに閉じ、$\dd z=\dd(z)$ なので完全である。それでも $\partial\bar\partial$ の像には入らない。原始形式の型に負の成分が要るからである。従って局所的な Dolbeault の補題からコンパクト Kähler の二重微分補題を代用することはできない。
:::

### 実ポテンシャルの一意性と積分ゼロ

連結コンパクト Kähler 多様体で $i\partial\bar\partial u=0$、$u$ は実とする。局所座標で縮約を取ると Laplacian に比例するので $\Delta_{\mathrm d}u=0$ である。従って $\|\mathrm du\|^2=0$ となり $u$ は定数である。従って実 $(1,1)$ 完全形式のポテンシャルは平均を零とすれば一意である。また $\alpha=i\partial\bar\partial u$ なら $\int_X\alpha\wedge\omega^{n-1}=0$ である。$n=1$ では積分が唯一の二次コホモロジー成分なので逆も成り立つが、高次元ではこの積分だけでは完全性を判定できない。

:::example[積分が零でも類が残る曲面]
$X=\mathbb P^1\times\mathbb P^1$ で $\omega=\omega_1+\omega_2$、$\alpha=\omega_1-\omega_2$ とする。$\alpha\wedge\omega=\omega_1^2-\omega_2^2=0$ であるが、第一因子の曲線で $\int\alpha=2\pi$ となる。従って $\alpha$ は完全でなく、大域ポテンシャルを持たない。
:::

### 同じ類の代表のポテンシャル

:::example[平坦トーラスの体積指定を積の形で解く]
複素次元 $2$ の正方トーラスで、座標 $z_j=x_j+iy_j$ の両実周期を $1$ とする。初期の計量行列を $I$ とし、体積比として次の積を指定する。

$$
\begin{aligned}
e^f&=\left(1 - \frac{\cos{\left(4 \pi x_{2} \right)}}{3}\right) \left(\frac{\cos{\left(2 \pi x_{1} \right)}}{3} + 1\right)
\end{aligned}
$$

各因子は正で、周期積分では余弦の平均が零なので因子ごとの平均は $1$ である。従って積の全平均も $1$ となる。ポテンシャルは各一変数の関数を足した形に取る。混合偏微分が消えるため、行列式が指定の積になるには対角成分だけを調整すればよい。

$$
\begin{aligned}
\varphi&=- \frac{\cos{\left(2 \pi x_{1} \right)}}{3 \pi^{2}} + \frac{\cos{\left(4 \pi x_{2} \right)}}{12 \pi^{2}} \\
g_\varphi&=\left[\begin{matrix}\frac{\cos{\left(2 \pi x_{1} \right)}}{3} + 1 & 0\\0 & 1 - \frac{\cos{\left(4 \pi x_{2} \right)}}{3}\end{matrix}\right] \\
\det(g_\varphi)&=\left(1 - \frac{\cos{\left(4 \pi x_{2} \right)}}{3}\right) \left(\frac{\cos{\left(2 \pi x_{1} \right)}}{3} + 1\right)
\end{aligned}
$$

各方向について $u_{z_j\bar z_j}=u_{x_jx_j}/4$ を計算すれば対角成分を再現する。従って新しい複素 Hessian は正定値で、同じ Kähler 類にある。さらに対数行列式の二階微分から、新しい Ricci 形式の係数は次の対角行列になる。

$$
\begin{aligned}
R_{1\overline{1}}&=\frac{\pi^{2} \left(\sin^{2}{\left(2 \pi x_{1} \right)} + \cos^{2}{\left(2 \pi x_{1} \right)} + 3 \cos{\left(2 \pi x_{1} \right)}\right)}{\left(\cos{\left(2 \pi x_{1} \right)} + 3\right)^{2}} \\
R_{2\overline{2}}&=\frac{4 \pi^{2} \left(\sin^{2}{\left(4 \pi x_{2} \right)} + \cos^{2}{\left(4 \pi x_{2} \right)} - 3 \cos{\left(4 \pi x_{2} \right)}\right)}{\left(\cos{\left(4 \pi x_{2} \right)} - 3\right)^{2}}
\end{aligned}
$$

新しい体積が旧体積と同じ平均を持つことと、Ricci が点ごとに零であることは異なる。各 $\lambda_j$ の正値性から計量は定義されるが、一般の点で Ricci は零にならない。差は大域の $i\partial\bar\partial\varphi$ であり、従って任意の閉な二次元サイクルにおいて差の積分は零になる。
:::

:::example[平坦トーラスの体積指定を積の形で解く]
複素次元 $3$ の正方トーラスで、座標 $z_j=x_j+iy_j$ の両実周期を $1$ とする。初期の計量行列を $I$ とし、体積比として次の積を指定する。

$$
\begin{aligned}
e^f&=\left(1 - \frac{\cos{\left(2 \pi x_{3} \right)}}{4}\right) \left(\frac{\cos{\left(4 \pi x_{1} \right)}}{2} + 1\right) \left(\frac{\cos{\left(2 \pi x_{2} \right)}}{5} + 1\right)
\end{aligned}
$$

各因子は正で、周期積分では余弦の平均が零なので因子ごとの平均は $1$ である。従って積の全平均も $1$ となる。ポテンシャルは各一変数の関数を足した形に取る。混合偏微分が消えるため、行列式が指定の積になるには対角成分だけを調整すればよい。

$$
\begin{aligned}
\varphi&=- \frac{\cos{\left(4 \pi x_{1} \right)}}{8 \pi^{2}} - \frac{\cos{\left(2 \pi x_{2} \right)}}{5 \pi^{2}} + \frac{\cos{\left(2 \pi x_{3} \right)}}{4 \pi^{2}} \\
g_\varphi&=\left[\begin{matrix}\frac{\cos{\left(4 \pi x_{1} \right)}}{2} + 1 & 0 & 0\\0 & \frac{\cos{\left(2 \pi x_{2} \right)}}{5} + 1 & 0\\0 & 0 & 1 - \frac{\cos{\left(2 \pi x_{3} \right)}}{4}\end{matrix}\right] \\
\det(g_\varphi)&=\left(1 - \frac{\cos{\left(2 \pi x_{3} \right)}}{4}\right) \left(\frac{\cos{\left(4 \pi x_{1} \right)}}{2} + 1\right) \left(\frac{\cos{\left(2 \pi x_{2} \right)}}{5} + 1\right)
\end{aligned}
$$

各方向について $u_{z_j\bar z_j}=u_{x_jx_j}/4$ を計算すれば対角成分を再現する。従って新しい複素 Hessian は正定値で、同じ Kähler 類にある。さらに対数行列式の二階微分から、新しい Ricci 形式の係数は次の対角行列になる。

$$
\begin{aligned}
R_{1\overline{1}}&=\frac{4 \pi^{2} \left(\sin^{2}{\left(4 \pi x_{1} \right)} + \cos^{2}{\left(4 \pi x_{1} \right)} + 2 \cos{\left(4 \pi x_{1} \right)}\right)}{\left(\cos{\left(4 \pi x_{1} \right)} + 2\right)^{2}} \\
R_{2\overline{2}}&=\frac{\pi^{2} \left(\sin^{2}{\left(2 \pi x_{2} \right)} + \cos^{2}{\left(2 \pi x_{2} \right)} + 5 \cos{\left(2 \pi x_{2} \right)}\right)}{\left(\cos{\left(2 \pi x_{2} \right)} + 5\right)^{2}} \\
R_{3\overline{3}}&=\frac{\pi^{2} \left(\sin^{2}{\left(2 \pi x_{3} \right)} + \cos^{2}{\left(2 \pi x_{3} \right)} - 4 \cos{\left(2 \pi x_{3} \right)}\right)}{\left(\cos{\left(2 \pi x_{3} \right)} - 4\right)^{2}}
\end{aligned}
$$

新しい体積が旧体積と同じ平均を持つことと、Ricci が点ごとに零であることは異なる。各 $\lambda_j$ の正値性から計量は定義されるが、一般の点で Ricci は零にならない。積の体積が変動しても、同じ類であるため全体積は変わらない。
:::

### Aeppli コホモロジーと Frölicher フィルトレーション

型 $(p,q)$ の :term[Aeppli コホモロジー]{en="Aeppli cohomology" yomi="あえぷりこほもろじー"} は $\Ker(\partial\bar\partial)$ を $\Img\partial+\Img\bar\partial$ で割った空間である。Bott–Chern から Aeppli への自然な写像がある。$X$ がコンパクト Kähler のとき、これは同型である。まず閉代表 $\alpha$ が Aeppli 側で零なら、$\alpha=\partial u+\bar\partial v$ なので共通の調和空間に直交する。:ref[thm-cg-hodge-decomposition] の証明から $\alpha$ は $\bar\partial$ の像に入る。従って :ref[thm-cg-ddbar] により Bott–Chern 側でも零になる。

次に $\partial\bar\partial\alpha=0$ とする。$\bar\partial\alpha$ は $\partial,\bar\partial$ で閉じ、$\bar\partial$ の像に入るため $\bar\partial\alpha=\partial\bar\partial\beta$ と書ける。$\alpha$ に $\partial\beta$ を足すと $\bar\partial$ で閉じた代表になる。この代表から $\bar\partial$ の像を引いて Dolbeault 調和代表を得ると、それは同時に $\partial$ でも閉じている。従って最初の写像は全射である。

de Rham 複体を、正則次数 $p$ が一定値以上の形式でフィルター付けする。付随する段階で最初に取る微分は $\bar\partial$ なので、スペクトル系列の第一頁は $E_1^{p,q}=H^q(X,\Omega_X^p)$ である。コンパクト Kähler の場合、各 $E_1$ の元は $\mathrm d$ でも閉じた純型の調和代表を持つ。この代表は全複体で既に閉じているので、それ以降の全ての微分は零になる。これが :term[Frölicher スペクトル系列]{en="Frölicher spectral sequence" yomi="ふれりひゃーすぺくとるけいれつ"} の第一頁での退化である。同じフィルトレーションをコホモロジーへ移すと、調和射影は型を保つため、正則次数が指定値以上の Hodge 成分の直和を得る。

## Hodge 対称性と Hodge 数

$h^{p,q}=\dim_\C H^{p,q}(X)$ を :term[Hodge 数]{en="Hodge number" yomi="ほっじすう"} と呼ぶ。$\Delta_\dd$ は実作用素なので、複素共役は調和形式の型を $(p,q)$ から $(q,p)$ に移す。従って $h^{p,q}=h^{q,p}$ である。また $H^q(X,\Omega_X^p)$ に :ref[thm-cg-serre-duality] を適用し、外積の完全ペアリングから $(\Omega_X^p)^*\otimes K_X\cong\Omega_X^{n-p}$ と分かるので、$h^{p,q}=h^{n-p,n-q}$ である。各反射は由来が異なり、前者は共役、後者は補次数の双対性である。

:::example[奇数 Betti 数は偶数]
$b_{2m+1}=\sum_{p+q=2m+1}h^{p,q}$ では $p=q$ の項がない。対称性で二つずつ等しい項を組にできるため偶数となる。例えば $S^1\times S^3$ は $b_1=1$ なので、どの複素構造に対しても Kähler 計量を持てない。
:::

:::example[二つの曲線の積]
種数 $g,h$ の曲線の積では $h^{1,0}=g+h$、$h^{2,0}=gh$、$h^{1,1}=2+2gh$ である。積計量で Laplacian は二因子の Laplacian の和になる。固有形式のテンソル積で固有値を足すと、零固有値には両方が調和のものだけが寄与する。各曲線の $h^{0,0}=h^{1,1}=1$、$h^{1,0}=h^{0,1}$ に種数を入れると、例えば $g=1,h=2$ では $h^{1,1}=6$ となる。
:::

### 積多様体の Hodge 数

二つのコンパクト Kähler 多様体に積計量を入れると外積形式について Laplacian は各因子の Laplacian の和になる。符号を確かめるには $\mathrm d(\alpha\wedge\beta)=\mathrm d_X\alpha\wedge\beta+(-1)^{\deg\alpha}\alpha\wedge\mathrm d_Y\beta$ と随伴の同様の式を使うと、混合項が消える。各 Laplacian の正規直交固有基底の積を取れば固有値が加わる。両方非負なので、積で固有値が零となるのは両方が零の場合だけである。従って Hodge 数を多項式 $\sum_{p,q}h^{p,q}s^pt^q$ に並べると、積に対して多項式が積になる。

:::example[二つの曲線の積]
種数 $2$ と種数 $3$ の曲線では、多項式は $1+2s+2t+st$ と $1+3s+3t+st$ である。積の $s^2$ 係数は $6$、$st$ 係数は $1+1+2\cdot3+2\cdot3=14$ である。従って $b_2=6+14+6=26$ であり、実 Künneth の公式 $1+1+(4)(6)=26$ と一致する。
:::

### 積と双対の次数を分ける

:::example[積の直線束でコホモロジーが残る次数]
$X=\prod_j\mathbb P^{n_j}$ とし、外部テンソル積の直線束 $L=\mathcal O(-2,1)$ を考える。$a$ のべきでコホモロジー次数を記録し、係数を各空間の次元とする。

一因子では符号がすべて揃った Taylor–Laurent 指数だけが残る。指数が非負の場合は全て $0$ 次に、全て負の場合は因子の最高次に属する。負指数を $-b_i-1$ と書けば、弱い合成の個数を二項係数で求められる。

$$
\begin{aligned}
j=1&\colon\quad \mathbb P^{1},\; d=-2,\quad e_i<0,\quad \sum_i e_i=-2,\quad P_j(a)=a \\
j=2&\colon\quad \mathbb P^{3},\; d=1,\quad e_i\ge0,\quad \sum_i e_i=1,\quad P_j(a)=4
\end{aligned}
$$

積では各因子の Dolbeault Laplacian の非負固有値を加える。零になるには各因子の固有値が零でなければならず、調和代表の次元と次数は多項式の積で合成される。

$$
\begin{aligned}
P_X(a)&=\prod_j P_j(a)=4 a \\
\dim_\C X&=4 \\
\chi(X,L)&=P_X(-1)=-4
\end{aligned}
$$

随伴の双対束 $K_X\otimes L^*$ は $\mathcal O(0,-5)$ である。$d\ge0$ の因子では双対側が最高次になり、$d\le-n-1$ の因子では逆になる。全消滅の範囲はそのままなので、双対は $q$ 次と $4-q$ 次を対応させる。Serre 双対は係数束も変える。Hodge 対称性が主張する共役による型の交換とは別の操作である。
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

最高 Chern 類は実接束の Euler 類であるため、ここで得た整数は位相的 Euler 標数にも等しい。曲面の $K_X$ は二次式の切断を持つ。$h^{2,0}$ を計算したら、共役から $h^{0,2}$ が同じになり、Euler 標数と比較できる。
:::

:::example[完全交差の標準束と Chern 数を計算する]
$X\subset\mathbb P^{4}$ を次数 $2,4$ の横断的な滑らかな完全交差とし、$H=c_1(\mathcal O_X(1))$ とする。法束は $\bigoplus_a\mathcal O_X(d_a)$ であり、次元は $2$ である。

接束の完全列から $c(T_X)c(N_{X/\mathbb P^N})=c(T_{\mathbb P^N}|_X)$ を得る。Euler 列により右辺は $(1+H)^{N+1}$ である。$c(T_X)=1+c_1H+\cdots+c_nH^n$ と書いて係数を比べる。

$$
\begin{aligned}
c(N)&=\left(2 H + 1\right) \left(4 H + 1\right) \\
c_{1}+(6)&=\binom{5}{1},\qquad c_{1}=-1 \\
c_{2}+(2)&=\binom{5}{2},\qquad c_{2}=8 \\
c(T_X)&=8 H^{2} - H + 1
\end{aligned}
$$

法束の行列式を取った随伴公式は $K_X=\mathcal O_X(1)$ である。これは $c(T_X)$ の一次の係数の符号を反転した値にも等しい。さらに $X$ を一般の $2$ 個の超平面で切り、交点数を積で数えると $\int_XH^2=8$ となる。

$$
\begin{aligned}
\int_X c_n(T_X)&=8\cdot8=64
\end{aligned}
$$

最高 Chern 類は実接束の Euler 類であるため、ここで得た整数は位相的 Euler 標数にも等しい。曲面の $K_X=\mathcal O_X(1)$ には五つの独立な切断がある。従って $h^{2,0}=h^{0,2}=5$ となる。
:::

## 強 Lefschetz 定理

:::theorem[強 Lefschetz 定理]{#thm-cg-hard-lefschetz}
$0\le k\le n$ において $L^{n-k}:H^k(X,\C)\to H^{2n-k}(X,\C)$、$L[\alpha]=[\omega\wedge\alpha]$ は同型である。さらに $H^k=\bigoplus_{j\ge0}L^jP^{k-2j}$ であり、$P^r=\Ker(L^{n-r+1}:H^r\to H^{2n-r+2})$、$0\le r\le n$ である。
:::

:::proof
$\omega$ が閉じているので外積はコホモロジーに降りる。:ref[thm-cg-kahler-identities] に続く交換関係から $L,\Lambda$ は調和空間を保つ。この有限次元空間上で :ref[prop-cg-lefschetz-linear] の計算をそのまま使える。実際、全次数 $k$ の交換子は $(n-k)\id$ であり、随伴関係も残る。各原始列は相補次数で対になり $L^{n-k}$ が同型になる。次数 $r\le n$ の列を $L^{n-r+1}$ で送ると、その次数で始まる原始成分だけが零になる。それより下の次数から来た列では、対応する底から測ったべきが列の長さを超えないため零にならない。従って記載の核が原始部分に等しく、分解が得られる。
:::

:::example[Betti 数の単調性]
$k<n$ で $L[\alpha]=0$ なら $L^{n-k}[\alpha]=0$ なので $[\alpha]=0$ である。従って $k\le n-1$ で $b_k\le b_{k+2}$ である。また $L$ は型を $(1,1)$ だけ上げるので、例えば $h^{2,1}\ge h^{1,0}$ は次元 $n\ge2$ で成り立つ。
:::

### 原始成分を Betti 数から逆算する

$k\le n$ では $L:H^{k-2}\to H^k$ が単射であり、Lefschetz 分解は $H^k=P^k\oplus L H^{k-2}$ を与える。従って $\dim P^k=b_k-b_{k-2}$ である。型を指定すると $\dim P^{p,q}=h^{p,q}-h^{p-1,q-1}$ となる。負の添字の空間は零とする。この差は、全次元から同じ Kähler 類を一回掛けて得られる成分を除いた次元である。

:::example[三つの曲線の積で二次を分ける]
$X=(\mathbb P^1)^3$ とし、各因子の規格化された面積類を $a,b,c$ とする。$\omega=a+b+c$ を使えば
$$
L(a)=ab+ac,\quad L(b)=ab+bc,\quad L(c)=ac+bc.
$$
基底 $(a,b,c)$ と $(ab,ac,bc)$ での行列は $\begin{pmatrix}1&1&0\\1&0&1\\0&1&1\end{pmatrix}$ で行列式が $-2$ である。従って $L:H^2\to H^4$ は有理および複素係数では同型であるが、整数格子では指数 $2$ の像を持つ。$P^2$ は $L^2$ の核であり、係数の和が零の二次元部分空間である。
:::

### Lefschetz 写像の次数を具体化する

:::example[積射影直線の原始部分と交点形式を計算する]
$X=(\mathbb P^1)^{3}$ で各因子の正規化された二次類を $a_j$ とする。$a_j^2=0$、$\int a_1\cdots a_n=1$ である。Kähler 類は次の正係数の線形結合とする。

$$
\begin{aligned}
\omega&=a_{1} + 3 a_{2} + 5 a_{3}
\end{aligned}
$$

強 Lefschetz の $L^{n-2}$ は $H^2$ を $H^{2n-2}$ へ送る。上位側の基底を昇順の添字の積とし、$a_i$ に $\omega^{n-2}$ を掛けると、同じ添字に二度出会う項が消える。残りの各添字集合の係数には並べ替えの数 $(n-2)!$ が掛かる。

$$
\begin{aligned}
[L^{1}]&=\left[\begin{matrix}3 & 1 & 0\\5 & 0 & 1\\0 & 5 & 3\end{matrix}\right] \\
\det[L^{1}]&=-30
\end{aligned}
$$

二次類の交点形式 $Q(a,b)=\int ab\omega^{n-2}$ は、行と列が異なるときに残りの添字の重みの積を数える。原始成分は $\omega$ に直交する核であり、列にその基底を並べると $N$ を得る。

$$
\begin{aligned}
[Q]&=\left[\begin{matrix}0 & 5 & 3\\5 & 0 & 1\\3 & 1 & 0\end{matrix}\right],\qquad Q(\omega,\omega)=90 \\
N&=\left[\begin{matrix}- \frac{1}{3} & - \frac{1}{5}\\1 & 0\\0 & 1\end{matrix}\right] \\
Q|_P&=\left[\begin{matrix}- \frac{10}{3} & -1\\-1 & - \frac{6}{5}\end{matrix}\right]
\end{aligned}
$$

最後の行列の先頭小行列式の符号は順に交替するので負定値である。これで実係数の正方向と原始の負方向を数値でも判定できる。三因子の積における $L$ の行列の行列式を非零な整数として確認する。
:::

:::example[積射影直線の原始部分と交点形式を計算する]
$X=(\mathbb P^1)^{4}$ で各因子の正規化された二次類を $a_j$ とする。$a_j^2=0$、$\int a_1\cdots a_n=1$ である。Kähler 類は次の正係数の線形結合とする。

$$
\begin{aligned}
\omega&=a_{1} + a_{2} + 2 a_{3} + 3 a_{4}
\end{aligned}
$$

強 Lefschetz の $L^{n-2}$ は $H^2$ を $H^{2n-2}$ へ送る。上位側の基底を昇順の添字の積とし、$a_i$ に $\omega^{n-2}$ を掛けると、同じ添字に二度出会う項が消える。残りの各添字集合の係数には並べ替えの数 $(n-2)!$ が掛かる。

$$
\begin{aligned}
[L^{2}]&=\left[\begin{matrix}4 & 4 & 2 & 0\\6 & 6 & 0 & 2\\12 & 0 & 6 & 4\\0 & 12 & 6 & 4\end{matrix}\right] \\
\det[L^{2}]&=-1728
\end{aligned}
$$

二次類の交点形式 $Q(a,b)=\int ab\omega^{n-2}$ は、行と列が異なるときに残りの添字の重みの積を数える。原始成分は $\omega$ に直交する核であり、列にその基底を並べると $N$ を得る。

$$
\begin{aligned}
[Q]&=\left[\begin{matrix}0 & 12 & 6 & 4\\12 & 0 & 6 & 4\\6 & 6 & 0 & 2\\4 & 4 & 2 & 0\end{matrix}\right],\qquad Q(\omega,\omega)=144 \\
N&=\left[\begin{matrix}-1 & - \frac{1}{2} & - \frac{1}{3}\\1 & 0 & 0\\0 & 1 & 0\\0 & 0 & 1\end{matrix}\right] \\
Q|_P&=\left[\begin{matrix}-24 & -6 & -4\\-6 & -6 & -2\\-4 & -2 & - \frac{8}{3}\end{matrix}\right]
\end{aligned}
$$

最後の行列の先頭小行列式の符号は順に交替するので負定値である。これで実係数の正方向と原始の負方向を数値でも判定できる。整数類を持つ場合にも整数格子上の全射性が自動で従うわけではない。
:::

### 原始形式に対する Hodge 星の公式

:::proposition[原始形式の Hodge 星]{#prop-cg-primitive-star}
全次数 $k=p+q\le n$ の原始 $(p,q)$ 形式 $\alpha$ に対して
$$
*\alpha=(-1)^{k(k+1)/2}i^{p-q}
          \frac{L^{n-k}\alpha}{(n-k)!}
$$
である。星作用素は実計量から作り、複素線形に拡張したものを使う。
:::

:::proof
一点でユニタリ枠を取り、各方向の二次形式を $\omega_j=i\theta_j\wedge\bar\theta_j$ と書く。基底の単項式を、$\theta_j$ だけを含む添字集合 $I$ と $\bar\theta_j$ だけを含む添字集合 $J$ によって分ける。これらは互いに素である。その一つを固定し、単独の因子を添字順に掛けたものを $\eta$、残った添字集合を $M$ とする。$s=|I|+|J|$、$m=|M|=n-s$ とすると、その成分は
$$
\eta\wedge\sum_{\substack{A\subset M\\ |A|=\ell}}c_A\omega_A,\qquad
\omega_A=\prod_{a\in A}\omega_a,\qquad k=s+2\ell
$$
と書ける。$L$ と $\Lambda$ は $I,J$ を変えない。$\Lambda$ は $\omega_a$ を一つ除く操作の和なので、原始性から大きさ $\ell-1$ の各 $B\subset M$ について $\sum_{A\supset B}c_A=0$ である。さらに $\Lambda$ のべきを作用させると、大きさ $\ell$ 未満の全ての $B$ で同じ和が零となる。各項が現れる回数は除いた添字の順列の数なので、非零の共通因子を除ける。

二方向を共に含む添字が事前に $\ell$ 個ある。従って $L^{n-k}/(n-k)!$ を作用させた後に、二方向を含む添字集合が $M\setminus C$ となる項の係数は $\sum_{A\cap C=\varnothing}c_A$ である。ここで $|C|=\ell$ とする。包含排除により
$$
\sum_{A\cap C=\varnothing}c_A
  =\sum_{B\subset C}(-1)^{|B|}\sum_{A\supset B}c_A
  =(-1)^\ell c_C
$$
である。最後の等式では $B\subsetneq C$ の項が原始性で消える。

実二次元の各方向で $*1=\omega_j$、$*\omega_j=1$、$*\theta_j=-i\theta_j$、$*\bar\theta_j=i\bar\theta_j$ である。積で対を組むとき、単独因子の余りを各方向へ戻す交換が $s(s-1)/2$ 回起きる。よって星は $c_C\eta\omega_C$ を、係数 $(-1)^{s(s-1)/2}(-i)^{|I|}i^{|J|}$ を掛けて $c_C\eta\omega_{M\setminus C}$ へ送る。$p-q=|I|-|J|$ なので、この係数は $(-1)^{s(s+1)/2}i^{p-q}$ である。上の包含排除の符号 $(-1)^\ell$ と合わせれば、$k=s+2\ell$ より主張の符号になる。各 $I,J$ の成分について成り立つので、和についても成り立つ。
:::

:::example[空き方向を一つ掛ける]
複素三次元で $\alpha=\omega_1-\omega_2$ は縮約が零の原始 $(1,1)$ 形式である。このとき $L\alpha=\omega_1\omega_3-\omega_2\omega_3$、一方 $*\alpha=\omega_2\omega_3-\omega_1\omega_3$ である。$k=2$ の符号 $(-1)^3=-1$ に一致する。複素四次元で $\beta=\theta_1\wedge\theta_2$ を取ると $L^2\beta/2!=\beta\wedge\omega_3\wedge\omega_4$ であり、$k=2$ の負号と $i^2$ の負号が相殺して $*\beta=+\beta\wedge\omega_3\wedge\omega_4$ となる。
:::

## Hodge 指数定理

$H^{1,1}(X,\R)$ を複素共役で固定される $(1,1)$ 類の空間とする。$n\ge2$ で実双線形形式 $Q(a,b)=\int_X a\wedge b\wedge\omega^{n-2}$ を考える。閉代表を完全形式だけ変えても値は変わらないのでコホモロジー上で定義される。

:::theorem[Hodge 指数定理]{#thm-cg-hodge-index}
$Q$ は $H^{1,1}(X,\R)$ 上で符号 $(1,h^{1,1}-1)$ を持つ。$\R[\omega]$ 上で正定値であり、その $Q$ 直交補では負定値である。
:::

:::proof
調和な実 $(1,1)$ 代表 $\alpha$ を取る。$\Lambda\alpha$ は調和な関数なので定数である。ユニタリ枠で $\alpha$ を実固有値 $(\lambda_j)$ を持つ対角形式にすると、$\alpha\wedge\omega^{n-1}=(\sum_j\lambda_j)\omega^n/n$ である。従って $Q(a,[\omega])=0$ なら $\Lambda\alpha=0$ であり、全点で $\sum_j\lambda_j=0$ である。

同じ枠で二乗を展開すると
$$
\alpha^2\wedge\omega^{n-2}
=\frac{(\sum_j\lambda_j)^2-\sum_j\lambda_j^2}{n(n-1)}\omega^n
$$
となる。これは二つの異なる方向を選ぶものだけが外積に残り、各順序の積 $\lambda_i\lambda_j$ が現れるためである。原始の場合、係数は非正で、$\alpha\ne0$ ならある開集合で負となる。よって積分は負である。一方 $Q([\omega],[\omega])=\int\omega^n>0$ であり、全ての類はこの直線と直交補の和に一意に分かれる。
:::

:::example[$\mathbb P^1\times\mathbb P^1$]
各因子からの規格化された面積類を $a,b$ とすると $a^2=b^2=0,\int ab=1$ である。$\omega$ の類を $a+b$ と取れば $(a+b)^2$ の積分は $2$、原始類 $a-b$ の二乗は $-2$ となる。
:::

### 実二次コホモロジー全体の正方向

実二次類は実 $(1,1)$ 部分と、$\alpha+\bar\alpha$、$\alpha$ は $(2,0)$ 型調和、という部分に分かれる。両部分は $Q$ について直交する。実際、$(1,1)$ と $(2,0)$ の積に $\omega^{n-2}$ を掛けると型は $(n+1,n-1)$ になり零である。$(2,0)$ 型の $\alpha=\sum_{i<j}a_{ij}\theta_i\wedge\theta_j$ については、$\alpha\wedge\bar\alpha\wedge\omega^{n-2}$ に同じ添字対だけが残る。外積の順を揃えると各項は $|a_{ij}|^2\omega_i\wedge\omega_j\wedge\omega^{n-2}$ であり、$\omega_i=i\theta_i\wedge\bar\theta_i$ なので係数は正である。$\alpha^2$ は $\omega^{n-2}$ を掛けると型の範囲を超えて零になる。従って $\alpha+\bar\alpha$ の二乗積分は非零な $\alpha$ について正である。この部分の実次元は $2h^{2,0}$ だから、実二次全体の符号は $(1+2h^{2,0},h^{1,1}-1)$ になる。

:::example[同じ交点形式の異なる制限]
曲面で $h^{2,0}=1$、$h^{1,1}=20$ なら、全実二次類の符号は $(3,19)$ である。一方、$(1,1)$ 部分に限ると $(1,19)$ となる。Hodge 指数定理で正方向が一つという記述は、後者の制限についてである。
:::

### 交点形式の負方向を計算する

$n=2$ で $X=\mathbb P^1\times\mathbb P^1$ とし、Kähler 類を $\omega=2a+3b$ と取る。$Q(xa+yb,\omega)=3x+2y$ である。従って原始部分は $2a-3b$ で張られ、その二乗は $-12$ である。Kähler 類自身の二乗は $12$ である。$\omega$ を変えると正方向として選ぶ直線と原始部分は変わるが、交点形式の符号は変わらない。

:::example[二つの因子類で不等式を確認する]
$a^2=b^2=0$、$ab=1$ として $D=a+4b$ を取ると、$D^2=8$、$D\cdot\omega=11$ である。従って $(D\cdot\omega)^2-D^2\omega^2=121-96=25$ となる。実際 $D-(11/12)\omega$ は原始であるから、その二乗は負であり、差の値はその二乗の $-12$ 倍になる。
:::

### 原始方向に交点を制限する

:::example[積射影直線の原始部分と交点形式を計算する]
$X=(\mathbb P^1)^{2}$ で各因子の正規化された二次類を $a_j$ とする。$a_j^2=0$、$\int a_1\cdots a_n=1$ である。Kähler 類は次の正係数の線形結合とする。

$$
\begin{aligned}
\omega&=a_{1} + 5 a_{2}
\end{aligned}
$$

強 Lefschetz の $L^{n-2}$ は $H^2$ を $H^{2n-2}$ へ送る。上位側の基底を昇順の添字の積とし、$a_i$ に $\omega^{n-2}$ を掛けると、同じ添字に二度出会う項が消える。残りの各添字集合の係数には並べ替えの数 $(n-2)!$ が掛かる。

$$
\begin{aligned}
[L^{0}]&=\left[\begin{matrix}1 & 0\\0 & 1\end{matrix}\right] \\
\det[L^{0}]&=1
\end{aligned}
$$

二次類の交点形式 $Q(a,b)=\int ab\omega^{n-2}$ は、行と列が異なるときに残りの添字の重みの積を数える。原始成分は $\omega$ に直交する核であり、列にその基底を並べると $N$ を得る。

$$
\begin{aligned}
[Q]&=\left[\begin{matrix}0 & 1\\1 & 0\end{matrix}\right],\qquad Q(\omega,\omega)=10 \\
N&=\left[\begin{matrix}- \frac{1}{5}\\1\end{matrix}\right] \\
Q|_P&=\left[\begin{matrix}- \frac{2}{5}\end{matrix}\right]
\end{aligned}
$$

最後の行列の先頭小行列式の符号は順に交替するので負定値である。これで実係数の正方向と原始の負方向を数値でも判定できる。曲面では $L^0$ は恒等である。負の方向の情報は写像の行列式ではなく交点形式から出る。
:::

:::example[積射影直線の原始部分と交点形式を計算する]
$X=(\mathbb P^1)^{3}$ で各因子の正規化された二次類を $a_j$ とする。$a_j^2=0$、$\int a_1\cdots a_n=1$ である。Kähler 類は次の正係数の線形結合とする。

$$
\begin{aligned}
\omega&=2 a_{1} + 2 a_{2} + 3 a_{3}
\end{aligned}
$$

強 Lefschetz の $L^{n-2}$ は $H^2$ を $H^{2n-2}$ へ送る。上位側の基底を昇順の添字の積とし、$a_i$ に $\omega^{n-2}$ を掛けると、同じ添字に二度出会う項が消える。残りの各添字集合の係数には並べ替えの数 $(n-2)!$ が掛かる。

$$
\begin{aligned}
[L^{1}]&=\left[\begin{matrix}2 & 2 & 0\\3 & 0 & 2\\0 & 3 & 2\end{matrix}\right] \\
\det[L^{1}]&=-24
\end{aligned}
$$

二次類の交点形式 $Q(a,b)=\int ab\omega^{n-2}$ は、行と列が異なるときに残りの添字の重みの積を数える。原始成分は $\omega$ に直交する核であり、列にその基底を並べると $N$ を得る。

$$
\begin{aligned}
[Q]&=\left[\begin{matrix}0 & 3 & 2\\3 & 0 & 2\\2 & 2 & 0\end{matrix}\right],\qquad Q(\omega,\omega)=72 \\
N&=\left[\begin{matrix}-1 & - \frac{2}{3}\\1 & 0\\0 & 1\end{matrix}\right] \\
Q|_P&=\left[\begin{matrix}-6 & -2\\-2 & - \frac{8}{3}\end{matrix}\right]
\end{aligned}
$$

最後の行列の先頭小行列式の符号は順に交替するので負定値である。これで実係数の正方向と原始の負方向を数値でも判定できる。原始部分の具体的な基底は係数から選べる。二次形式の負定値性はその基底の選び方によらない。
:::

:::example[積射影直線の原始部分と交点形式を計算する]
$X=(\mathbb P^1)^{4}$ で各因子の正規化された二次類を $a_j$ とする。$a_j^2=0$、$\int a_1\cdots a_n=1$ である。Kähler 類は次の正係数の線形結合とする。

$$
\begin{aligned}
\omega&=a_{1} + 2 a_{2} + 4 a_{3} + a_{4}
\end{aligned}
$$

強 Lefschetz の $L^{n-2}$ は $H^2$ を $H^{2n-2}$ へ送る。上位側の基底を昇順の添字の積とし、$a_i$ に $\omega^{n-2}$ を掛けると、同じ添字に二度出会う項が消える。残りの各添字集合の係数には並べ替えの数 $(n-2)!$ が掛かる。

$$
\begin{aligned}
[L^{2}]&=\left[\begin{matrix}16 & 8 & 4 & 0\\4 & 2 & 0 & 4\\8 & 0 & 2 & 8\\0 & 8 & 4 & 16\end{matrix}\right] \\
\det[L^{2}]&=-3072
\end{aligned}
$$

二次類の交点形式 $Q(a,b)=\int ab\omega^{n-2}$ は、行と列が異なるときに残りの添字の重みの積を数える。原始成分は $\omega$ に直交する核であり、列にその基底を並べると $N$ を得る。

$$
\begin{aligned}
[Q]&=\left[\begin{matrix}0 & 8 & 4 & 16\\8 & 0 & 2 & 8\\4 & 2 & 0 & 4\\16 & 8 & 4 & 0\end{matrix}\right],\qquad Q(\omega,\omega)=192 \\
N&=\left[\begin{matrix}- \frac{1}{2} & - \frac{1}{4} & -1\\1 & 0 & 0\\0 & 1 & 0\\0 & 0 & 1\end{matrix}\right] \\
Q|_P&=\left[\begin{matrix}-8 & -2 & -8\\-2 & -2 & -4\\-8 & -4 & -32\end{matrix}\right]
\end{aligned}
$$

最後の行列の先頭小行列式の符号は順に交替するので負定値である。これで実係数の正方向と原始の負方向を数値でも判定できる。原始部分の階数が増えても、符号を小行列式で確認する方法は同じである。
:::

### Hodge–Riemann の符号を型ごとに求める

原始調和な $\alpha\ne0$ が $(p,q)$ 型、$k=p+q\le n$ なら、:ref[prop-cg-primitive-star] に $\bar\alpha$ を代入する。$\int\alpha\wedge *\bar\alpha=\|\alpha\|^2$ を使うと
$$
(-1)^{k(k-1)/2}i^{p-q}
\int_X\alpha\wedge\bar\alpha\wedge\omega^{n-k}
=(n-k)!\;\|\alpha\|^2>0
$$
となる。ここでは $(-1)^k i^{q-p}=i^{p-q}$ を使った。これが原始成分上の Hodge–Riemann の双線形関係の正値性である。型 $(p,q)$ と型 $(r,s)$ の原始類を共役なしに積分すると、$p+r=k$ と $q+s=k$ でなければ積の型が最高型にならない。従って原始部分の符号付き交点形式は実 Hodge 構造を偏極する。Kähler 類が有理類なら、後述の通り原始部分にも有理構造が入り、この形式を有理係数の偏極にできる。$k=2$ では $(1,1)$ に対して負、$(2,0)$ に対して正となり、前段の指数定理の符号が戻る。

:::example[三次の階数と原始の符号]
複素三次元では $(3,0)$ 型の $\alpha$ に対し $\Lambda\alpha=0$ は型から自動である。$k=3$、$p-q=3$ なので前の係数は $(-1)^3i^3=i$ である。従って $i\int\alpha\wedge\bar\alpha>0$ である。一方、原始 $(2,1)$ 型では同じ次数でも係数は $-i$ になる。二つを合わせて扱う場合に、型に応じた位相因子を落とすと正値性を失う。
:::

## Hodge 構造と Hodge 予想の概観

:::definition[純粋 Hodge 構造]{#def-cg-hodge-structure}
有限次元有理ベクトル空間 $V_\Q$ に対して $V_\C=\bigoplus_{p+q=k}V^{p,q}$、$\overline{V^{p,q}}=V^{q,p}$ を与えたものを重み $k$ の :term[純粋 Hodge 構造]{en="pure Hodge structure" yomi="じゅんすいほっじこうぞう"} という。
:::

Hodge 分解と共役対称性から $H^k(X,\Q)$ は重み $k$ の純粋 Hodge 構造を持つ。正則写像の引き戻しは純型の閉代表を同じ型の閉形式に移す。調和射影で型が保たれるので、Hodge 構造の各成分を保つ。型 $(r,r)$ に入り、元が有理係数であるクラスを有理 Hodge 類と呼ぶ。

:::proposition[Lefschetz の $(1,1)$ 定理]{#prop-cg-lefschetz-one-one}
整数二次コホモロジー類が $\operatorname{Pic}(X)$ の $c_1$ の像に入ることと、その複素化が型 $(1,1)$ に属することは同値である。
:::

:::proof
指数完全列では $c_1$ の像は $H^2(X,\Z)\to H^2(X,\mathcal O_X)$ の核である。後者の写像は、de Rham 分解から Dolbeault 分解への $(0,q)$ 型射影で記述できる。実際、外微分の $(0,q+1)$ 成分は $(0,q)$ 成分への $\bar\partial$ であり、定数関数は構造層に入るので分解の射を与える。従って二次の写像は $H^{0,2}$ 成分への射影である。整数類の複素化は実類なので、$(0,2)$ 成分が零なら共役の $(2,0)$ 成分も零になる。逆は明らかであり、torsion の類も複素化で零になるため核に含まれる。
:::

滑らかな複素射影多様体上で、複素余次元 $r$ の代数サイクルのコホモロジー類は型 $(r,r)$ に属する。滑らかな部分多様体では、積分で定義する Poincaré 双対が相補的な型 $(n-r,n-r)$ のみに非零の値を取り得ることから従う。特異な部分の場合も滑らかな部分の積分で閉カレントを定めて同じ結論を得る。特異部分まで含めた基本類の構成にはカレントの理論が必要であり、:cite[griffiths-harris] を参照する。

Hodge 予想は、滑らかな複素射影多様体の $H^{2r}(X,\Q)\cap H^{r,r}(X)$ の元は、余次元 $r$ の代数サイクルの類の有理線形結合であるという主張である。$r=1$ では、上の定理が有理類の整数倍を正則直線束の第一 Chern 類として実現する。その束を代数化して予想を示す段階は、後の [Chow の定理と GAGA](/math/complex-geometry/vanishing-embedding/) で扱う。一般の余次元での主張は未解決である。仮定には射影性が含まれ、係数は有理数である。

:::example[射影空間での予想]
$\mathbb P^n$ の偶数次元セルによる計算では各 $H^{2r}$ が一次元である。余次元 $r$ の線形部分空間の類は超平面類の $r$ 乗であり、補余次元の一般の線形部分空間との横断交点が一点なので非零である。従って各有理 Hodge 類はその有理倍である。
:::

### Hodge フィルトレーションと整数類の制約

複素化 $V_\C$ の Hodge 構造から
$$
F^rV_\C=\bigoplus_{p\ge r}V^{p,k-p}
$$
という減少フィルトレーションを作る。逆にこのフィルトレーションが $V_\C=F^rV_\C\oplus\overline{F^{k-r+1}V_\C}$ を全ての $r$ で満たすなら、$V^{p,q}=F^pV_\C\cap\overline{F^qV_\C}$、$p+q=k$ で分解を復元できる。実際、$F^p$ の各元を $F^{p+1}\oplus\overline{F^q}$ で分けると、後ろの成分は $F^p$ にも入る。従って $F^p=F^{p+1}\oplus(F^p\cap\overline{F^q})$ となる。これを次数の上から繰り返せば直和分解が得られ、共役で添字が交換される。

:::example[二次コホモロジーの次元から読めること]
$h^{2,0}=1,h^{1,1}=20,h^{0,2}=1$ なら $F^2$ は一次元、$F^1$ は二十一次元、$F^0$ は二十二次元である。正則直線束の $c_1$ は実 $(1,1)$ 部分に入るので、相互に独立な整数の第一 Chern 類は高々二十個である。これは上界であり、$H^{1,1}$ の実部分に入る整数格子の階数がいつも二十になるとは限らない。
:::

### 複素三次元までの有理 Hodge 類

複素三次元の滑らかな射影多様体 $X$ で $h=c_1(\mathcal O_X(1))$ を用いる。曲率の正値性により $h$ は実係数で Kähler 類の正定数倍である。強 Lefschetz 定理で $h\cup-:H^2(X,\Q)\to H^4(X,\Q)$ は同型となる。実際、複素化した写像が同型で、その行列は有理係数であるからである。この写像は型 $(1,1)$ を加えるので、有理 $(2,2)$ 類の逆像は有理 $(1,1)$ 類である。分母を払って Lefschetz の $(1,1)$ 定理を適用すれば、逆像は因子類の有理倍として得られる。それに超平面類を掛けると代数サイクルの類になる。従って次元三まででは、余次元一の定理と強 Lefschetz 定理で全ての余次元を扱える。最高次数は一点の類、零次は成分ごとの定数で生成される。

### Hodge フィルトレーションから型を取り出す

重み $k$ で $F^p=\bigoplus_{r\ge p}V^{r,k-r}$ とすれば $V^{p,q}=F^p\cap\bar F^q$ である。実際、右辺の成分には $r\ge p$、$k-r\ge q$ が必要であり、$p+q=k$ から $r=p$ しか残らない。実なクラスの場合は $F^p$ への所属から共役の条件も自動で従う。例えば重み二で、実なクラスが $F^1$ に属するなら $(0,2)$ が零であり、共役により $(2,0)$ も零となる。

:::example[楕円曲線の二つの部分空間]
$E=\C/(\Z+\tau\Z)$ とし、$z=x+\tau y$ と書くと $F^1=\C[\mathrm dx+\tau\mathrm dy]$、$\bar F^1=\C[\mathrm dx+\bar\tau\mathrm dy]$ である。$\Im\tau>0$ なのでこの二直線は異なる。例えば $\tau=2i$ では $F^1$ には有理係数の非零類がない。もし $a,b\in\Q$ が $b/a=2i$ を満たせば $a,b$ が共に零になるからである。
:::

### 整数類と有理類の基底

:::example[積射影直線の原始部分と交点形式を計算する]
$X=(\mathbb P^1)^{3}$ で各因子の正規化された二次類を $a_j$ とする。$a_j^2=0$、$\int a_1\cdots a_n=1$ である。Kähler 類は次の正係数の線形結合とする。

$$
\begin{aligned}
\omega&=2 a_{1} + a_{2} + 5 a_{3}
\end{aligned}
$$

強 Lefschetz の $L^{n-2}$ は $H^2$ を $H^{2n-2}$ へ送る。上位側の基底を昇順の添字の積とし、$a_i$ に $\omega^{n-2}$ を掛けると、同じ添字に二度出会う項が消える。残りの各添字集合の係数には並べ替えの数 $(n-2)!$ が掛かる。

$$
\begin{aligned}
[L^{1}]&=\left[\begin{matrix}1 & 2 & 0\\5 & 0 & 2\\0 & 5 & 1\end{matrix}\right] \\
\det[L^{1}]&=-20
\end{aligned}
$$

二次類の交点形式 $Q(a,b)=\int ab\omega^{n-2}$ は、行と列が異なるときに残りの添字の重みの積を数える。原始成分は $\omega$ に直交する核であり、列にその基底を並べると $N$ を得る。

$$
\begin{aligned}
[Q]&=\left[\begin{matrix}0 & 5 & 1\\5 & 0 & 2\\1 & 2 & 0\end{matrix}\right],\qquad Q(\omega,\omega)=60 \\
N&=\left[\begin{matrix}-2 & - \frac{2}{5}\\1 & 0\\0 & 1\end{matrix}\right] \\
Q|_P&=\left[\begin{matrix}-20 & -2\\-2 & - \frac{4}{5}\end{matrix}\right]
\end{aligned}
$$

最後の行列の先頭小行列式の符号は順に交替するので負定値である。これで実係数の正方向と原始の負方向を数値でも判定できる。全ての二次基底 $a_j$ は各因子の超平面の逆像で、実 $(1,1)$ 型である。この積ではそれらが整数二次群を張る。
:::

:::example[完全交差の標準束と Chern 数を計算する]
$X\subset\mathbb P^{4}$ を次数 $5$ の横断的な滑らかな完全交差とし、$H=c_1(\mathcal O_X(1))$ とする。法束は $\bigoplus_a\mathcal O_X(d_a)$ であり、次元は $3$ である。

接束の完全列から $c(T_X)c(N_{X/\mathbb P^N})=c(T_{\mathbb P^N}|_X)$ を得る。Euler 列により右辺は $(1+H)^{N+1}$ である。$c(T_X)=1+c_1H+\cdots+c_nH^n$ と書いて係数を比べる。

$$
\begin{aligned}
c(N)&=5 H + 1 \\
c_{1}+(5)&=\binom{5}{1},\qquad c_{1}=0 \\
c_{2}+(0)&=\binom{5}{2},\qquad c_{2}=10 \\
c_{3}+(50)&=\binom{5}{3},\qquad c_{3}=-40 \\
c(T_X)&=- 40 H^{3} + 10 H^{2} + 1
\end{aligned}
$$

法束の行列式を取った随伴公式は $K_X=\mathcal O_X(0)$ である。これは $c(T_X)$ の一次の係数の符号を反転した値にも等しい。さらに $X$ を一般の $3$ 個の超平面で切り、交点数を積で数えると $\int_XH^3=5$ となる。

$$
\begin{aligned}
\int_X c_n(T_X)&=-40\cdot5=-200
\end{aligned}
$$

最高 Chern 類は実接束の Euler 類であるため、ここで得た整数は位相的 Euler 標数にも等しい。超平面類は二次の有理 $(1,1)$ 類である。この例で二次群全体を決める方法と、中間三次がどう増えるかは、後の [Lefschetz の超平面定理](/math/complex-geometry/vanishing-embedding/) と [五次超曲面の例](/math/complex-geometry/deformations/) で扱う。
:::

### 偏極が有理係数になる条件

Hodge–Riemann の形式を有理 Hodge 構造の偏極として使うには、積分で定義した双線形形式が有理数値であることも確認する。Kähler 類が有理類であれば、これは有理係数コホモロジーのカップ積と基本類の評価から従う。原始部分も有理線形写像 $L^{n-k+1}$ の核として有理構造を持つ。例えば $\mathbb P^2$ の超平面類を二倍した Kähler 類では、$H^0$ 上の $Q(1,1)=\int(2H)^2=4$ であり、有理数値で正である。正の無理数倍として $\omega=\sqrt[3]{2}\,H$ を使うと同じ形式は正であるが、その値は $2^{2/3}$ であり有理数値という条件は失われる。

## 演習問題

::::exercise{level=basic}
複素二次元の平坦トーラスの $h^{1,1}$ と $b_2$ を求めよ。
:::solution
定数係数形式から $h^{1,1}=\binom21^2=4$ である。$h^{2,0}=h^{0,2}=1$ なので $b_2=1+4+1=6$ になる。
:::
::::

::::exercise{level=basic}
$h^{2,0}=3,h^{1,1}=7$ とする。$b_2$ を求めよ。
:::solution
共役対称性により $h^{0,2}=3$ である。二次 Hodge 分解の三項を足し $b_2=3+7+3=13$ となる。
:::
::::

::::exercise{level=basic}
$n=4,k=2$ で強 Lefschetz 同型の作用を具体的に書け。
:::solution
$n-k=2$ なので $[\alpha]\mapsto[\omega^2\wedge\alpha]$ である。全次数は $2+4=6$ になり、$H^2(X,\C)\cong H^6(X,\C)$ を与える。
:::
::::

::::exercise{level=basic}
$\mathbb P^1\times\mathbb P^1$ で $2a-3b$ の二乗の積分を求めよ。
:::solution
$a^2=b^2=0$、二次形式の外積は可換である。従って $(2a-3b)^2=-12ab$ となり、$\int ab=1$ から積分は $-12$ である。
:::
::::

::::exercise{level=advanced}
$X$ がコンパクト Kähler 多様体で $h^{0,2}=0$ なら、全整数二次類が正則直線束の第一 Chern 類であることを示せ。
:::solution
指数完全列の $H^2(X,\Z)\to H^2(X,\mathcal O_X)$ の終域は Dolbeault 同型で $H^{0,2}$ なので零である。従って核は全体で、前の写像である $c_1$ の像と一致する。
:::
::::

::::exercise{level=advanced}
調和な $(0,1)$ 形式で $\bar\partial$ 完全なものは零であることを示せ。
:::solution
$\alpha=\bar\partial u$ なら、調和性で $\bar\partial^*\alpha=0$ なので $\|\alpha\|^2=\langle\bar\partial u,\alpha\rangle=\langle u,\bar\partial^*\alpha\rangle=0$ である。従って $\alpha$ は全点で零になる。
:::
::::

::::exercise{level=advanced}
$n\ge2$ で、実 $(1,1)$ 類 $a$ が $Q(a,a)>0$ を満たすとき、その $Q$ 直交補が負定値であることを示せ。
:::solution
指数定理の基底で $Q(x,x)=x_0^2-|x'|^2$ と書く。$Q(a,a)>0$ なので $|a'|<|a_0|$ である。直交する $x$ では $a_0x_0=a'\cdot x'$ だから $|x_0|\le|a'|\,|x'|/|a_0|$ となる。$x\ne0$ なら $x'\ne0$ で、この上界が $|x'|$ より小さいため $Q(x,x)<0$ となる。
:::
::::

::::exercise{level=advanced}
$\mathbb P^1\times\mathbb P^1$ で Kähler 類 $2a+b$ に直交する実 $(1,1)$ 類の一次元空間を求めよ。
:::solution
$(xa+yb)(2a+b)$ の積分は $x+2y$ である。従って $x=-2y$ が直交条件で、$2a-b$ が空間を張る。この生成元の二乗の積分は $-4$ となり負定値である。
:::
::::

::::exercise{level=advanced}
$X=(\mathbb P^1)^3$ の超平面類を $a,b,c$、Kähler 類を $\omega=2a+3b+c$ とする。$H^2(X,\R)$ において $a+2b$ を Kähler 類の成分と原始成分に分け、原始成分の二次形式 $Q(P,P)=\int P^2\omega$ を求めよ。

:::solution
$Q$ に関する直交射影の係数は $Q(a+2b,\omega)/Q(\omega,\omega)$ である。$a^2=b^2=c^2=0$ と $\int abc=1$ から
$$
\omega^2=12ab+4ac+6bc,\qquad \int\omega^3=6\cdot2\cdot3=36.
$$
従って $\int(a+2b)\omega^2=6+8=14$ となり、係数は $7/18$ である。原始成分は
$$
P=a+2b-\tfrac7{18}\omega
=\tfrac29a+\tfrac56b-\tfrac7{18}c .
$$
二乗を全部展開する代わりに、直交性から
$$
Q(P,P)=Q(a+2b,a+2b)-\frac{14^2}{36}
=4-\frac{49}{9}=-\frac{13}{9}.
$$
原始成分は非零であり、計算した値は Hodge 指数定理の負定値性と一致する。全体の $Q(a+2b,a+2b)=4$ が正でも、Kähler 成分を除いたクラスでは負になることに注意する。
:::
::::

## 確認テスト

::::quiz{answer="6" tol="1e-6"}
複素二次元の平坦トーラスの $b_2$ を求めよ。
:::explanation
定数係数の $2$ 形式は実座標四つから二つを選んで得られ、全て調和である。従って $b_2=\binom42=6$ である。
:::
::::

::::quiz
コンパクト Kähler 多様体の $H^{p,q}$ と $H^{q,p}$ を対応させる操作を選べ。

- [ ] $\omega$ との外積。
- [x] 複素共役。
- [ ] $\bar\partial$。

:::explanation
実 Laplacian は共役と可換であり、共役は $(p,q)$ 型を $(q,p)$ 型に移す。従って調和空間に制限して対称性を与える。
:::
::::

::::quiz{answer="-2" tol="1e-6"}
$\mathbb P^1\times\mathbb P^1$ の規格化された面積類 $a,b$ について $\int(a-b)^2$ を求めよ。
:::explanation
積の各因子が実二次元なので $a^2=b^2=0$ である。交差項は $-2ab$ で、$\int ab=1$ より答えは $-2$ となる。
:::
::::

::::quiz
$\partial\bar\partial$ 補題を用いる際の正しい前提を選べ。

- [x] コンパクト Kähler 多様体上で、純型の形式が $\partial,\bar\partial$ で閉じ、さらに $\dd$ の像に入る。
- [ ] 任意の複素多様体上で、任意の関数が $\dd$ で閉じている。
- [ ] コンパクトな空間で、形式の係数がすべて実数である。

:::explanation
前者では調和成分が消え、Green 作用素と Kähler 恒等式により $\partial\bar\partial$ の像に入る。閉性と完全性は別の条件である。
:::
::::

::::quiz{answer="26" tol="1e-6"}
種数 $2$ と $3$ の二つのコンパクト Riemann 面の積の第二 Betti 数を求めよ。

:::explanation
各面の Betti 数を $1,2g,1$ とし、Künneth の公式で総次数二の積を取る。二次と零次、一次と一次、零次と二次から $1+4\cdot6+1=26$ を得る。
:::
::::

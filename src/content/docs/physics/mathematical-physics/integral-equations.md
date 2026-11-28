---
title: "積分方程式"
description: "積分方程式の核の構造を使い、有限次元化、反復、可解条件、弱特異核の逆変換を導く。"
notation:
  - symbol: 'K(x,t)'
    meaning: "積分作用素の核"
  - symbol: 'T'
    meaning: "核 K が定める作用素"
  - symbol: 'T^*'
    meaning: "随伴作用素"
  - symbol: 'I^\alpha'
    meaning: "片側の分数積分"
---

## この章の目標

- 積分区間と未知関数の位置から方程式を分類する。
- 分離核や反復級数によって具体的な解を求める。
- 斉次解が存在する場合に非斉次項の可解条件を導く。
- Abel 核と微分方程式の積分化を計算する。

## Fredholm 型と Volterra 型の積分方程式

### 積分区間と未知量

$u(x)-\lambda\int_a^b K(x,t)u(t)\dd t=f(x)$ を第 2 種の :term[フレドホルム積分方程式]{en="Fredholm integral equation" yomi="ふれどほるむせきぶんほうていしき"} と呼ぶ。上限が $x$ の場合は :term[ヴォルテラ積分方程式]{en="Volterra integral equation" yomi="ゔぉるてらせきぶんほうていしき"} である。第 1 種は積分項だけが未知関数を含む。例えば $f(x)=\int_0^x u(t)\dd t$ は第 1 種で $u=f'$、$f(0)=0$ が必要である。微分で誤差が拡大する。実際 $f_n(x)=\sin(nx)/n$ は一様に零に近づくが、$u_n=\cos(nx)$ の一様ノルムは $1$ のままである。一方、$u(x)=1+\int_0^x u(t)\dd t$ は第 2 種で、微分すると $u'=u,u(0)=1$ だから $e^x$ である。上限が現在時刻なら過去の値だけが寄与する。

### 第一種と第二種の逆の違い

第一種方程式 $Tu=f$ で核が滑らかな場合、$T$ は入力の細かい振動を弱める。例えば $Tu(x)=\int_0^x u(t)\dd t$ は速い余弦を振幅が小さい正弦に送る。逆は微分であり、与えた右辺に高周波が混じると増幅される。第二種 $u-\lambda Tu=f$ では、有界作用素 $T$ に対し $I-\lambda T$ が全単射なら、有界逆写像定理により逆も有界である。

この違いはコンパクト性で説明できる。無限次元 Hilbert 空間でコンパクト作用素 $T$ が単射とする。もし $\|u\|\le C\|Tu\|$ が成り立つなら、正規直交列 $(e_n)$ から $\|Te_n-Te_m\|\ge\sqrt2/C$ を得る。これは $(Te_n)$ が収束部分列を持つことと矛盾する。さらに値域が閉じていれば、値域への有界全単射 $T$ は有界な逆を持つ。従って値域は閉じず、近いデータの全てが解を持つとも限らない。

以下の $L^2$ の内積は $\inner{u}{v}=\int u\bar v$ とし、第一引数を線形とする。連続核と連続データについて点ごとの解を求める箇所では、$C([a,b])$ に一様ノルムを用いる。

:::example[積分で振幅が下がる例（1）]
単位区間で $u_m(x)=\cos(2\pi m x)$ を入力とする。

$$
m=1
$$

Volterra 積分を行って原点で零になることを確認する。

$$
\begin{aligned}
Tu_m(x)\big|_{m=1} &= \frac{\sin{\left(2 \pi x \right)}}{2 \pi}
\end{aligned}
$$

入力の二乗ノルムを求める。

$$
\begin{aligned}
\|u_m\|_2^2 &= \frac{1}{2}
\end{aligned}
$$

出力を同じ測度で積分し、逆での増幅が無限に大きくなることを確かめる。

$$
\begin{aligned}
\|Tu_m\|_2^2\big|_{m=1} &= \frac{1}{8 \pi^{2}}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
Tu_m(x)\big|_{x=\frac{1}{3}} &= \frac{\sqrt{3}}{4 \pi} \\
\frac{\partial^1 (Tu_m(x))}{\partial x^1} &= \cos{\left(2 \pi x \right)} \\
\left.\frac{\partial^1 (Tu_m(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - \frac{1}{2}
\end{aligned}
$$
:::

:::example[積分で振幅が下がる例（2）]
単位区間で $u_m(x)=\cos(2\pi m x)$ を入力とする。

$$
m=2
$$

Volterra 積分を行って原点で零になることを確認する。

$$
\begin{aligned}
Tu_m(x)\big|_{m=2} &= \frac{\sin{\left(4 \pi x \right)}}{4 \pi}
\end{aligned}
$$

入力の二乗ノルムを求める。

$$
\begin{aligned}
\|u_m\|_2^2 &= \frac{1}{2}
\end{aligned}
$$

出力を同じ測度で積分し、逆での増幅が無限に大きくなることを確かめる。

$$
\begin{aligned}
\|Tu_m\|_2^2\big|_{m=2} &= \frac{1}{32 \pi^{2}}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
Tu_m(x)\big|_{x=\frac{1}{3}} &= - \frac{\sqrt{3}}{8 \pi} \\
\frac{\partial^1 (Tu_m(x))}{\partial x^1} &= \cos{\left(4 \pi x \right)} \\
\left.\frac{\partial^1 (Tu_m(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - \frac{1}{2}
\end{aligned}
$$
:::

:::example[積分で振幅が下がる例（3）]
単位区間で $u_m(x)=\cos(2\pi m x)$ を入力とする。

$$
m=3
$$

Volterra 積分を行って原点で零になることを確認する。

$$
\begin{aligned}
Tu_m(x)\big|_{m=3} &= \frac{\sin{\left(6 \pi x \right)}}{6 \pi}
\end{aligned}
$$

入力の二乗ノルムを求める。

$$
\begin{aligned}
\|u_m\|_2^2 &= \frac{1}{2}
\end{aligned}
$$

出力を同じ測度で積分し、逆での増幅が無限に大きくなることを確かめる。

$$
\begin{aligned}
\|Tu_m\|_2^2\big|_{m=3} &= \frac{1}{72 \pi^{2}}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
Tu_m(x)\big|_{x=\frac{1}{3}} &= 0 \\
\frac{\partial^1 (Tu_m(x))}{\partial x^1} &= \cos{\left(6 \pi x \right)} \\
\left.\frac{\partial^1 (Tu_m(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= 1
\end{aligned}
$$
:::

:::example[積分で振幅が下がる例（4）]
単位区間で $u_m(x)=\cos(2\pi m x)$ を入力とする。

$$
m=4
$$

Volterra 積分を行って原点で零になることを確認する。

$$
\begin{aligned}
Tu_m(x)\big|_{m=4} &= \frac{\sin{\left(8 \pi x \right)}}{8 \pi}
\end{aligned}
$$

入力の二乗ノルムを求める。

$$
\begin{aligned}
\|u_m\|_2^2 &= \frac{1}{2}
\end{aligned}
$$

出力を同じ測度で積分し、逆での増幅が無限に大きくなることを確かめる。

$$
\begin{aligned}
\|Tu_m\|_2^2\big|_{m=4} &= \frac{1}{128 \pi^{2}}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
Tu_m(x)\big|_{x=\frac{1}{3}} &= \frac{\sqrt{3}}{16 \pi} \\
\frac{\partial^1 (Tu_m(x))}{\partial x^1} &= \cos{\left(8 \pi x \right)} \\
\left.\frac{\partial^1 (Tu_m(x))}{\partial x^1}\right|_{x=\frac{1}{3}} &= - \frac{1}{2}
\end{aligned}
$$
:::

## 分離核の積分方程式

### 核を有限個の積に分ける

$K(x,t)=\sum_{j=1}^n a_j(x)\overline{b_j(t)}$ では $Tu=\sum a_j\inner{u}{b_j}$ である。$c_j=\inner{u}{b_j}$ と置くと $u=f+\lambda\sum a_jc_j$ なので $c_i-\lambda\sum_j\inner{a_j}{b_i}c_j=\inner{f}{b_i}$ となる。逆にこの有限系の解から $u$ を作ると、その $c_i$ が実際の積分係数である。

:::example[定数核と共鳴]
$[0,1]$ で $K=1$ なら $c=\int_0^1u$ は $(1-\lambda)c=\int_0^1 f$ に従う。$\lambda\ne1$ では $u=f+\lambda\int f/(1-\lambda)$ である。$\lambda=1$ では $\int f=0$ が必要十分で、$u=f+c$ と任意の定数を足せる。例えば $f=x$, $\lambda=1/2$ なら $u=x+1/2$ である。
:::

:::example[非対称核で随伴の条件を計算する]
$[0,1]$ で $K(x,t)=x$ とする。$Tu=x\int u$ なので $Tx=x/2$、一方 $T^*u=\int t u(t)\dd t$ なので $T^*1=1/2$ である。$\lambda=2$ での斉次解は $x$ の倍数、随伴の斉次解は定数である。実際、$c=\int u$ を使うと $u=f+2cx$ であり、積分すると $c=\int f+c$ になる。従って必要十分条件は $\int f=0$ である。例えば $f=x-1/2$ は解を持ち、$f=x-2/3$ は $x$ と直交しても解を持たない。この例では元の核との直交性だけでは可解性を判定できない。
:::

### 分離核の応答を一つのスカラーにする

$Tu=p\int q(t)u(t)\dd t$ と書けるなら、$c=\int qu$ と置いて $u=f+\lambda p c$ を使う。両辺に $q$ を掛けて積分すると $(1-\lambda\int qp)c=\int qf$ になる。積分値が零だと $T^2=0$ となるため、任意の $\lambda$ で逆が一次の多項式になる場合もある。

:::example[単位区間で冪を持つ分離核（1）]
$K(x,t)=x^m t$、右辺 $f(x)=1$ とする。

$$
m=1
$$

唯一の非零固有値は二つの因子の積の積分である。

$$
\begin{aligned}
\int_0^1 t^{m+1}\dd t\big|_{m=1} &= \frac{1}{3}
\end{aligned}
$$

源と左因子の内積を求め、スカラー方程式を解く。

$$
\begin{aligned}
c\big|_{m=1} &= - \frac{3}{2 \lambda - 6}
\end{aligned}
$$

積分係数を元の式に戻して解を求める。

$$
\begin{aligned}
f+\lambda x^m c\big|_{m=1} &= \frac{- \frac{3 \lambda x}{2} + \lambda - 3}{\lambda - 3}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
f+\lambda x^m c\big|_{x=\frac{1}{3}} &= \frac{\lambda - 6}{2 \left(\lambda - 3\right)} \\
\frac{\partial^1 (f+\lambda x^m c)}{\partial x^1} &= - \frac{3 \lambda}{2 \lambda - 6} \\
\left.\frac{\partial^1 (f+\lambda x^m c)}{\partial x^1}\right|_{x=\frac{1}{3}} &= - \frac{3 \lambda}{2 \lambda - 6}
\end{aligned}
$$
:::

:::example[単位区間で冪を持つ分離核（2）]
$K(x,t)=x^m t$、右辺 $f(x)=1$ とする。

$$
m=2
$$

唯一の非零固有値は二つの因子の積の積分である。

$$
\begin{aligned}
\int_0^1 t^{m+1}\dd t\big|_{m=2} &= \frac{1}{4}
\end{aligned}
$$

源と左因子の内積を求め、スカラー方程式を解く。

$$
\begin{aligned}
c\big|_{m=2} &= - \frac{2}{\lambda - 4}
\end{aligned}
$$

積分係数を元の式に戻して解を求める。

$$
\begin{aligned}
f+\lambda x^m c\big|_{m=2} &= \frac{- 2 \lambda x^{2} + \lambda - 4}{\lambda - 4}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
f+\lambda x^m c\big|_{x=\frac{1}{3}} &= \frac{7 \lambda - 36}{9 \left(\lambda - 4\right)} \\
\frac{\partial^1 (f+\lambda x^m c)}{\partial x^1} &= - \frac{4 \lambda x}{\lambda - 4} \\
\left.\frac{\partial^1 (f+\lambda x^m c)}{\partial x^1}\right|_{x=\frac{1}{3}} &= - \frac{4 \lambda}{3 \lambda - 12}
\end{aligned}
$$
:::

:::example[単位区間で冪を持つ分離核（3）]
$K(x,t)=x^m t$、右辺 $f(x)=1$ とする。

$$
m=3
$$

唯一の非零固有値は二つの因子の積の積分である。

$$
\begin{aligned}
\int_0^1 t^{m+1}\dd t\big|_{m=3} &= \frac{1}{5}
\end{aligned}
$$

源と左因子の内積を求め、スカラー方程式を解く。

$$
\begin{aligned}
c\big|_{m=3} &= - \frac{5}{2 \lambda - 10}
\end{aligned}
$$

積分係数を元の式に戻して解を求める。

$$
\begin{aligned}
f+\lambda x^m c\big|_{m=3} &= \frac{- \frac{5 \lambda x^{3}}{2} + \lambda - 5}{\lambda - 5}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
f+\lambda x^m c\big|_{x=\frac{1}{3}} &= \frac{49 \lambda - 270}{54 \left(\lambda - 5\right)} \\
\frac{\partial^1 (f+\lambda x^m c)}{\partial x^1} &= - \frac{15 \lambda x^{2}}{2 \lambda - 10} \\
\left.\frac{\partial^1 (f+\lambda x^m c)}{\partial x^1}\right|_{x=\frac{1}{3}} &= - \frac{5 \lambda}{6 \lambda - 30}
\end{aligned}
$$
:::

:::example[単位区間で冪を持つ分離核（4）]
$K(x,t)=x^m t$、右辺 $f(x)=1$ とする。

$$
m=4
$$

唯一の非零固有値は二つの因子の積の積分である。

$$
\begin{aligned}
\int_0^1 t^{m+1}\dd t\big|_{m=4} &= \frac{1}{6}
\end{aligned}
$$

源と左因子の内積を求め、スカラー方程式を解く。

$$
\begin{aligned}
c\big|_{m=4} &= - \frac{3}{\lambda - 6}
\end{aligned}
$$

積分係数を元の式に戻して解を求める。

$$
\begin{aligned}
f+\lambda x^m c\big|_{m=4} &= \frac{- 3 \lambda x^{4} + \lambda - 6}{\lambda - 6}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
f+\lambda x^m c\big|_{x=\frac{1}{3}} &= \frac{2 \left(13 \lambda - 81\right)}{27 \left(\lambda - 6\right)} \\
\frac{\partial^1 (f+\lambda x^m c)}{\partial x^1} &= - \frac{12 \lambda x^{3}}{\lambda - 6} \\
\left.\frac{\partial^1 (f+\lambda x^m c)}{\partial x^1}\right|_{x=\frac{1}{3}} &= - \frac{4 \lambda}{9 \lambda - 54}
\end{aligned}
$$
:::

### 有限階数の Fredholm 行列式

分離核を $K(x,t)=\sum_{j=1}^n a_j(x)b_j(t)$ とする。$A:\C^n\to L^2$ を $Ac=\sum a_j c_j$、$B:L^2\to\C^n$ を $(Bu)_j=\int b_j u$ と書くと核の作用素は $AB$ である。$I-\lambda BA$ が可逆なとき
$$
(I-\lambda AB)^{-1}=I+\lambda A(I-\lambda BA)^{-1}B
$$
になる。左から $I-\lambda AB$ を掛けると、$A-\lambda ABA=A(I-\lambda BA)$ によって追加分が $\lambda AB$ となり相殺する。右から掛けた場合も $B-\lambda BAB=(I-\lambda BA)B$ で相殺する。逆に $\lambda\ne0$ で $(I-\lambda BA)c=0$ かつ $c\ne0$ なら $Ac\ne0$ である。実際 $Ac=0$ なら $c=\lambda BAc=0$ となる。その $Ac$ は $(I-\lambda AB)$ の核に入る。従って可逆性は有限次の行列式 $D(\lambda)=\det(I-\lambda BA)$ の零点で判定できる。

例えば $K(x,t)=1+xt$、区間 $[0,1]$ なら $BA=\begin{pmatrix}1&1/2\\1/2&1/3\end{pmatrix}$ であり、$D(\lambda)=1-4\lambda/3+\lambda^2/12$ となる。$\lambda=1$ では $D=-1/4\ne0$ である。右辺を $f=1$ とすると $Bf=(1,1/2)^{\mathsf T}$ で、$(I-BA)c=Bf$ の解は $c=(-11/3,-2)^{\mathsf T}$ になる。従って $u=1+c_1+xc_2=-8/3-2x$ である。検算では $\int u=-11/3$、$\int tu=-2$ を先に求めて式に戻せばよい。

## Neumann 級数

### 小さい作用素または三角形の反復

任意の Banach 空間で $\|\lambda T\|<1$ なら $\sum_{n=0}^\infty(\lambda T)^n$ はノルム収束する。有限和を $I-\lambda T$ に掛けると $I-(\lambda T)^{N+1}$ なので、極限は逆作用素である。打ち切り誤差は $u$ に対し $\|\lambda T\|^{N+1}\|f\|/(1-\|\lambda T\|)$ 以下である。

Volterra 核が $|K|\le M$ なら、$T^nu(x)$ は $a<t_n<\cdots<t_1<x$ 上の積分である。積分領域は立方体を $n!$ 個の添字順に分けた一つなので、$\|T^n\|_\infty\le M^n(b-a)^n/n!$ である。従って $\lambda$ の大きさによらず級数は収束し、同じ有限積で逆を得る。例えば $K=1$ は $T^n1=(x-a)^n/n!$ なので $u=1+\lambda Tu$ の解は $e^{\lambda(x-a)}$ である。

### 反復の核をまとめてから積分する

積分作用素の $n$ 乗には、反復核 $K_n$ が対応する。$K_1=K$、Fredholm 型では $K_{n+1}(x,t)=\int_a^b K(x,s)K_n(s,t)\dd s$ である。Volterra 型では積分区間を $t\le s\le x$ にする。Fubini の定理で作用素の合成と核の反復が一致する。Volterra 核が $|K|\le M$ を満たせば、帰納法で
$$
|K_n(x,t)|\le M^n\frac{(x-t)^{n-1}}{(n-1)!}
$$
となる。$n$ から $n+1$ へ移る際は $(s-t)^{n-1}$ を $s\in[t,x]$ で積分する。従って :term[レゾルベント核]{en="resolvent kernel" yomi="れぞるべんとかく"} $R_\lambda=\sum_{n\ge1}\lambda^{n-1}K_n$ は有限三角領域で一様収束し、$u=f+\lambda\int_a^x R_\lambda(x,t)f(t)\dd t$ と表せる。

Volterra の $K(x,t)=1$ では $T^n$ の核は $(x-t)^{n-1}/(n-1)!$（$a\le t\le x,\ n\ge1$）である。$t$ を一番内側で固定すると、残りの $n-1$ 変数の順序領域の体積がこの値になる。従って $(I-\lambda T)^{-1}f=f+\lambda\int_a^x e^{\lambda(x-t)}f(t)\dd t$ を得る。例えば $a=0,\lambda=1,f=x$ なら $\int_0^x e^{x-t}t\dd t=e^x-1-x$ だから、解は $e^x-1$ になる。微分すると $u'=1+u$、原点では $u=0$ であり、元の式を満たす。

### Volterra の反復核を積分する

$Tu(x)=\int_0^xu(t)\dd t$ とする。帰納法で $T^j u(x)=\int_0^x (x-t)^{j-1}u(t)\dd t/(j-1)!$ を得る。一段進むと積分順序を交換して、固定した $t$ と $x$ の間で前の多項式を積分する。階乗が出るので、有限区間では $\lambda$ の大きさを問わず逆の級数が収束する。級数を核の内側で足すと、残りの逆の核は $\lambda e^{\lambda(x-t)}$ になる。

:::example[冪入力への逐次積分（1）]
入力 $f(x)=x^m$ を固定し、最初と次の反復を求める。

$$
m=1
$$

最初の積分の原始関数を原点で規格化する。

$$
\begin{aligned}
Tf\big|_{m=1} &= \frac{x^{2}}{2}
\end{aligned}
$$

もう一回積分すると別の分母が加わる。

$$
\begin{aligned}
T^2f\big|_{m=1} &= \frac{x^{3}}{6}
\end{aligned}
$$

三回までの和を取り、原点の値を確かめる。

$$
\begin{aligned}
f+\lambda Tf+\lambda^2T^2f+\lambda^3T^3f\big|_{m=1} &= \frac{x \left(\lambda^{3} x^{3} + 4 \lambda^{2} x^{2} + 12 \lambda x + 24\right)}{24}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
Tf\big|_{x=\frac{1}{3}} &= \frac{1}{18} \\
\frac{\partial^1 (Tf)}{\partial x^1} &= x \\
\left.\frac{\partial^1 (Tf)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{1}{3} \\
\frac{\partial^2 (Tf)}{\partial x^2} &= 1 \\
\left.\frac{\partial^2 (Tf)}{\partial x^2}\right|_{x=\frac{1}{3}} &= 1
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
T^2f\big|_{x=\frac{1}{3}} &= \frac{1}{162} \\
\frac{\partial^1 (T^2f)}{\partial x^1} &= \frac{x^{2}}{2} \\
\left.\frac{\partial^1 (T^2f)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{1}{18} \\
\frac{\partial^2 (T^2f)}{\partial x^2} &= x \\
\left.\frac{\partial^2 (T^2f)}{\partial x^2}\right|_{x=\frac{1}{3}} &= \frac{1}{3}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
f+\lambda Tf+\lambda^2T^2f+\lambda^3T^3f\big|_{x=\frac{1}{3}} &= \frac{\lambda^{3}}{1944} + \frac{\lambda^{2}}{162} + \frac{\lambda}{18} + \frac{1}{3} \\
\frac{\partial^1 (f+\lambda Tf+\lambda^2T^2f+\lambda^3T^3f)}{\partial x^1} &= \frac{\lambda^{3} x^{3}}{6} + \frac{\lambda^{2} x^{2}}{2} + \lambda x + 1 \\
\left.\frac{\partial^1 (f+\lambda Tf+\lambda^2T^2f+\lambda^3T^3f)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{\lambda^{3}}{162} + \frac{\lambda^{2}}{18} + \frac{\lambda}{3} + 1 \\
\frac{\partial^2 (f+\lambda Tf+\lambda^2T^2f+\lambda^3T^3f)}{\partial x^2} &= \lambda \left(\frac{\lambda^{2} x^{2}}{2} + \lambda x + 1\right) \\
\left.\frac{\partial^2 (f+\lambda Tf+\lambda^2T^2f+\lambda^3T^3f)}{\partial x^2}\right|_{x=\frac{1}{3}} &= \frac{\lambda \left(\lambda^{2} + 6 \lambda + 18\right)}{18}
\end{aligned}
$$
:::

:::example[冪入力への逐次積分（2）]
入力 $f(x)=x^m$ を固定し、最初と次の反復を求める。

$$
m=2
$$

最初の積分の原始関数を原点で規格化する。

$$
\begin{aligned}
Tf\big|_{m=2} &= \frac{x^{3}}{3}
\end{aligned}
$$

もう一回積分すると別の分母が加わる。

$$
\begin{aligned}
T^2f\big|_{m=2} &= \frac{x^{4}}{12}
\end{aligned}
$$

三回までの和を取り、原点の値を確かめる。

$$
\begin{aligned}
f+\lambda Tf+\lambda^2T^2f+\lambda^3T^3f\big|_{m=2} &= \frac{x^{2} \left(\lambda^{3} x^{3} + 5 \lambda^{2} x^{2} + 20 \lambda x + 60\right)}{60}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
Tf\big|_{x=\frac{1}{3}} &= \frac{1}{81} \\
\frac{\partial^1 (Tf)}{\partial x^1} &= x^{2} \\
\left.\frac{\partial^1 (Tf)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{1}{9} \\
\frac{\partial^2 (Tf)}{\partial x^2} &= 2 x \\
\left.\frac{\partial^2 (Tf)}{\partial x^2}\right|_{x=\frac{1}{3}} &= \frac{2}{3}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
T^2f\big|_{x=\frac{1}{3}} &= \frac{1}{972} \\
\frac{\partial^1 (T^2f)}{\partial x^1} &= \frac{x^{3}}{3} \\
\left.\frac{\partial^1 (T^2f)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{1}{81} \\
\frac{\partial^2 (T^2f)}{\partial x^2} &= x^{2} \\
\left.\frac{\partial^2 (T^2f)}{\partial x^2}\right|_{x=\frac{1}{3}} &= \frac{1}{9}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
f+\lambda Tf+\lambda^2T^2f+\lambda^3T^3f\big|_{x=\frac{1}{3}} &= \frac{\lambda^{3}}{14580} + \frac{\lambda^{2}}{972} + \frac{\lambda}{81} + \frac{1}{9} \\
\frac{\partial^1 (f+\lambda Tf+\lambda^2T^2f+\lambda^3T^3f)}{\partial x^1} &= \frac{x \left(\lambda^{3} x^{3} + 4 \lambda^{2} x^{2} + 12 \lambda x + 24\right)}{12} \\
\left.\frac{\partial^1 (f+\lambda Tf+\lambda^2T^2f+\lambda^3T^3f)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{\lambda^{3}}{972} + \frac{\lambda^{2}}{81} + \frac{\lambda}{9} + \frac{2}{3} \\
\frac{\partial^2 (f+\lambda Tf+\lambda^2T^2f+\lambda^3T^3f)}{\partial x^2} &= \frac{\lambda^{3} x^{3}}{3} + \lambda^{2} x^{2} + 2 \lambda x + 2 \\
\left.\frac{\partial^2 (f+\lambda Tf+\lambda^2T^2f+\lambda^3T^3f)}{\partial x^2}\right|_{x=\frac{1}{3}} &= \frac{\lambda^{3}}{81} + \frac{\lambda^{2}}{9} + \frac{2 \lambda}{3} + 2
\end{aligned}
$$
:::

:::example[冪入力への逐次積分（3）]
入力 $f(x)=x^m$ を固定し、最初と次の反復を求める。

$$
m=3
$$

最初の積分の原始関数を原点で規格化する。

$$
\begin{aligned}
Tf\big|_{m=3} &= \frac{x^{4}}{4}
\end{aligned}
$$

もう一回積分すると別の分母が加わる。

$$
\begin{aligned}
T^2f\big|_{m=3} &= \frac{x^{5}}{20}
\end{aligned}
$$

三回までの和を取り、原点の値を確かめる。

$$
\begin{aligned}
f+\lambda Tf+\lambda^2T^2f+\lambda^3T^3f\big|_{m=3} &= \frac{x^{3} \left(\lambda^{3} x^{3} + 6 \lambda^{2} x^{2} + 30 \lambda x + 120\right)}{120}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
Tf\big|_{x=\frac{1}{3}} &= \frac{1}{324} \\
\frac{\partial^1 (Tf)}{\partial x^1} &= x^{3} \\
\left.\frac{\partial^1 (Tf)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{1}{27} \\
\frac{\partial^2 (Tf)}{\partial x^2} &= 3 x^{2} \\
\left.\frac{\partial^2 (Tf)}{\partial x^2}\right|_{x=\frac{1}{3}} &= \frac{1}{3}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
T^2f\big|_{x=\frac{1}{3}} &= \frac{1}{4860} \\
\frac{\partial^1 (T^2f)}{\partial x^1} &= \frac{x^{4}}{4} \\
\left.\frac{\partial^1 (T^2f)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{1}{324} \\
\frac{\partial^2 (T^2f)}{\partial x^2} &= x^{3} \\
\left.\frac{\partial^2 (T^2f)}{\partial x^2}\right|_{x=\frac{1}{3}} &= \frac{1}{27}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
f+\lambda Tf+\lambda^2T^2f+\lambda^3T^3f\big|_{x=\frac{1}{3}} &= \frac{\lambda^{3}}{87480} + \frac{\lambda^{2}}{4860} + \frac{\lambda}{324} + \frac{1}{27} \\
\frac{\partial^1 (f+\lambda Tf+\lambda^2T^2f+\lambda^3T^3f)}{\partial x^1} &= x^{2} \left(\frac{\lambda^{3} x^{3}}{20} + \frac{\lambda^{2} x^{2}}{4} + \lambda x + 3\right) \\
\left.\frac{\partial^1 (f+\lambda Tf+\lambda^2T^2f+\lambda^3T^3f)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{\lambda^{3}}{4860} + \frac{\lambda^{2}}{324} + \frac{\lambda}{27} + \frac{1}{3} \\
\frac{\partial^2 (f+\lambda Tf+\lambda^2T^2f+\lambda^3T^3f)}{\partial x^2} &= \frac{x \left(\lambda^{3} x^{3} + 4 \lambda^{2} x^{2} + 12 \lambda x + 24\right)}{4} \\
\left.\frac{\partial^2 (f+\lambda Tf+\lambda^2T^2f+\lambda^3T^3f)}{\partial x^2}\right|_{x=\frac{1}{3}} &= \frac{\lambda^{3}}{324} + \frac{\lambda^{2}}{27} + \frac{\lambda}{3} + 2
\end{aligned}
$$
:::

:::example[冪入力への逐次積分（4）]
入力 $f(x)=x^m$ を固定し、最初と次の反復を求める。

$$
m=4
$$

最初の積分の原始関数を原点で規格化する。

$$
\begin{aligned}
Tf\big|_{m=4} &= \frac{x^{5}}{5}
\end{aligned}
$$

もう一回積分すると別の分母が加わる。

$$
\begin{aligned}
T^2f\big|_{m=4} &= \frac{x^{6}}{30}
\end{aligned}
$$

三回までの和を取り、原点の値を確かめる。

$$
\begin{aligned}
f+\lambda Tf+\lambda^2T^2f+\lambda^3T^3f\big|_{m=4} &= \frac{x^{4} \left(\lambda^{3} x^{3} + 7 \lambda^{2} x^{2} + 42 \lambda x + 210\right)}{210}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
Tf\big|_{x=\frac{1}{3}} &= \frac{1}{1215} \\
\frac{\partial^1 (Tf)}{\partial x^1} &= x^{4} \\
\left.\frac{\partial^1 (Tf)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{1}{81} \\
\frac{\partial^2 (Tf)}{\partial x^2} &= 4 x^{3} \\
\left.\frac{\partial^2 (Tf)}{\partial x^2}\right|_{x=\frac{1}{3}} &= \frac{4}{27}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
T^2f\big|_{x=\frac{1}{3}} &= \frac{1}{21870} \\
\frac{\partial^1 (T^2f)}{\partial x^1} &= \frac{x^{5}}{5} \\
\left.\frac{\partial^1 (T^2f)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{1}{1215} \\
\frac{\partial^2 (T^2f)}{\partial x^2} &= x^{4} \\
\left.\frac{\partial^2 (T^2f)}{\partial x^2}\right|_{x=\frac{1}{3}} &= \frac{1}{81}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
f+\lambda Tf+\lambda^2T^2f+\lambda^3T^3f\big|_{x=\frac{1}{3}} &= \frac{\lambda^{3}}{459270} + \frac{\lambda^{2}}{21870} + \frac{\lambda}{1215} + \frac{1}{81} \\
\frac{\partial^1 (f+\lambda Tf+\lambda^2T^2f+\lambda^3T^3f)}{\partial x^1} &= x^{3} \left(\frac{\lambda^{3} x^{3}}{30} + \frac{\lambda^{2} x^{2}}{5} + \lambda x + 4\right) \\
\left.\frac{\partial^1 (f+\lambda Tf+\lambda^2T^2f+\lambda^3T^3f)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{\lambda^{3}}{21870} + \frac{\lambda^{2}}{1215} + \frac{\lambda}{81} + \frac{4}{27} \\
\frac{\partial^2 (f+\lambda Tf+\lambda^2T^2f+\lambda^3T^3f)}{\partial x^2} &= x^{2} \left(\frac{\lambda^{3} x^{3}}{5} + \lambda^{2} x^{2} + 4 \lambda x + 12\right) \\
\left.\frac{\partial^2 (f+\lambda Tf+\lambda^2T^2f+\lambda^3T^3f)}{\partial x^2}\right|_{x=\frac{1}{3}} &= \frac{\lambda^{3}}{1215} + \frac{\lambda^{2}}{81} + \frac{4 \lambda}{27} + \frac{4}{3}
\end{aligned}
$$
:::

## Fredholm の交代定理

### 逆が無いときにどの源が許されるか

:::theorem[Fredholm の交代定理]{#thm-mp-fredholm-alternative}
Hilbert 空間上のコンパクト作用素 $T$ と複素数 $\lambda$ について $A=I-\lambda T$ とする。$A$ の値域は閉じており、その余次元と $\dim\Ker A$ は有限で等しい。$Au=f$ の可解条件は $f\perp\Ker A^*$ であり、斉次核が零なら全ての $f$ に対して解は一意である。
:::

:::proof
任意の $\eta>0$ に対し、$T$ が単位球を移す集合は有限個の半径 $\eta$ の球で覆える。その中心の張る空間への直交射影を $P$ とすると $\|T-PT\|\le\eta$ である。$F=\lambda PT$ と $E=\lambda(T-PT)$ に分け、$\eta$ を小さくして $\|E\|<1$ とする。$A=(I-E)(I-S)$、$S=(I-E)^{-1}F$ は有限階数である。

$M$ を $S$ の値域を含む有限次元空間とし、$M\oplus M^\perp$ で $I-S$ は $\begin{pmatrix}B&C\\0&I\end{pmatrix}$ と書ける。核は $(m,0)$ で $Bm=0$ となるもの全体である。値域は $(y,z)$ で $y-Cz\in\operatorname{ran} B$ となるもの全体だから閉じ、余次元は $\dim M-\rank B=\dim\Ker B$ である。可逆な $I-E$ を掛けても核次元、閉性、余次元は変わらない。最後に $\operatorname{ran} A$ の直交補空間は随伴の定義から $\Ker A^*$ である。値域の閉性を使うと主張の可解条件を得る。
:::

### 共鳴で零になる一成分

自己共役な階数一核 $Tu=\eta\inner{u}{e}e$ で $\|e\|=1$、$\eta\in\R$ とする。$I-\lambda T$ は $e^\perp$ で恒等、$e$ で値 $1-\lambda\eta$ を掛ける。従って $\lambda\eta=1$ で可解な源は $e^\perp$ に限られる。可解なときは一解に $e$ の任意の倍数を足せる。非対称核では $T^*$ の固有値 $\overline{\lambda}^{-1}$ に属する固有空間への直交性が可解条件になる。

### 連続核の Fredholm 行列式と余因子

連続核 $K$ を有限区間 $[a,b]^2$ に取る。有限階数の場合の行列式に対応する関数と、その余因子に対応する核を次の級数で定める。
$$
\begin{aligned}
D(\lambda)&=\sum_{n=0}^\infty\frac{(-\lambda)^n}{n!}
\int_{[a,b]^n}\det[K(x_i,x_j)]_{i,j=1}^{n}\,\dd x_1\cdots\dd x_n,\\
M_\lambda(x,t)&=\sum_{n=0}^\infty\frac{(-\lambda)^n}{n!}
\int_{[a,b]^n}\det[K(\xi_i,\tau_j)]_{i,j=0}^{n}\,\dd x_1\cdots\dd x_n,
\end{aligned}
$$
ここで $\xi_0=x$、$\tau_0=t$、$\xi_i=\tau_i=x_i$ $(i\ge1)$ である。$n=0$ の項はそれぞれ $1,K(x,t)$ とする。$|K|\le M$ なら $n$ 次の行列式は、列の長さの積による Hadamard の不等式で $M^n n^{n/2}$ 以下である。$n!$ で割った級数は Stirling の公式により、$\lambda$ が有界な範囲で絶対一様収束する。余因子の行列は $n+1$ 次だが、同じ評価で収束する。

$M_\lambda$ の行列式を第 $0$ 行で展開する。第 $0$ 列の項は $K(x,t)D$ になる。他の $n$ 項は、選んだ積分変数を $s$ に付け替え、残った $s$ の行を先頭に移すと、同じ $-K(x,s)$ 倍の余因子になる。この $n$ が $n!$ を $(n-1)!$ に変え、負号と $(-\lambda)$ が $\lambda$ を与える。列で展開する計算も同様に行うと
$$
M_\lambda=K D+\lambda K*M_\lambda
=K D+\lambda M_\lambda*K
$$
を得る。$*$ は積分による核の合成である。従って $D(\lambda)\ne0$ なら $R_\lambda=M_\lambda/D$ は逆の核であり、両側から掛けて $(I-\lambda T)^{-1}=I+\lambda R_\lambda$ と確認できる。最後の式の $R_\lambda$ は核の定める積分作用素を表す。

逆に $D$ の零点 $\lambda_0$ では逆は存在しない。実際、定義を微分して $D'(\lambda)=-\int_a^b M_\lambda(x,x)\dd x$ となる。もし $\lambda_0$ で $L^2$ 上の逆が存在するなら、連続右辺の解は $u=f+\lambda_0 Tu$ から連続である。この逆は有界逆写像定理で $C([a,b])$ 上も有界であり、Neumann 級数によって近傍でも連続に有界である。そこで $D\ne0$ の点では $R_\lambda(\cdot,t)=(I-\lambda T)^{-1}K(\cdot,t)$ が $t$ も含め一様有界なので、$D'/D=-\int R_\lambda(x,x)\dd x$ は有界となる。しかし $D(0)=1$ であり、非恒等零の正則関数の零点では $D'/D$ は極を持つ。従って可逆性の仮定と矛盾する。

:::example[定数核の余因子]
$K=1$ では大きさ $2$ 以上の行列は全て行が等しく、行列式が消える。従って $D=1-\lambda(b-a)$、$M_\lambda=1$ である。$R_\lambda=1/[1-\lambda(b-a)]$ から分離核と同じ解を得る。
:::

:::example[共鳴を引き算で避ける右辺（1）]
$e$ の方向を定数とする。右辺に $x^m$ を選び、定数成分を引いて平均零にする。

$$
m=1
$$

単位区間で平均を計算する。

$$
\begin{aligned}
\langle x^m,1\rangle\big|_{m=1} &= \frac{1}{2}
\end{aligned}
$$

平均を除けば可解条件を満たす。

$$
\begin{aligned}
f_\perp\big|_{m=1} &= x - \frac{1}{2}
\end{aligned}
$$

右辺のノルムを求め、引いた成分のノルムを除く。

$$
\begin{aligned}
\|f_\perp\|^2\big|_{m=1} &= \frac{1}{12}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
f_\perp\big|_{x=\frac{1}{3}} &= - \frac{1}{6} \\
\frac{\partial^1 (f_\perp)}{\partial x^1} &= 1 \\
\left.\frac{\partial^1 (f_\perp)}{\partial x^1}\right|_{x=\frac{1}{3}} &= 1
\end{aligned}
$$
:::

:::example[共鳴を引き算で避ける右辺（2）]
$e$ の方向を定数とする。右辺に $x^m$ を選び、定数成分を引いて平均零にする。

$$
m=2
$$

単位区間で平均を計算する。

$$
\begin{aligned}
\langle x^m,1\rangle\big|_{m=2} &= \frac{1}{3}
\end{aligned}
$$

平均を除けば可解条件を満たす。

$$
\begin{aligned}
f_\perp\big|_{m=2} &= x^{2} - \frac{1}{3}
\end{aligned}
$$

右辺のノルムを求め、引いた成分のノルムを除く。

$$
\begin{aligned}
\|f_\perp\|^2\big|_{m=2} &= \frac{4}{45}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
f_\perp\big|_{x=\frac{1}{3}} &= - \frac{2}{9} \\
\frac{\partial^1 (f_\perp)}{\partial x^1} &= 2 x \\
\left.\frac{\partial^1 (f_\perp)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{2}{3}
\end{aligned}
$$
:::

:::example[共鳴を引き算で避ける右辺（3）]
$e$ の方向を定数とする。右辺に $x^m$ を選び、定数成分を引いて平均零にする。

$$
m=3
$$

単位区間で平均を計算する。

$$
\begin{aligned}
\langle x^m,1\rangle\big|_{m=3} &= \frac{1}{4}
\end{aligned}
$$

平均を除けば可解条件を満たす。

$$
\begin{aligned}
f_\perp\big|_{m=3} &= x^{3} - \frac{1}{4}
\end{aligned}
$$

右辺のノルムを求め、引いた成分のノルムを除く。

$$
\begin{aligned}
\|f_\perp\|^2\big|_{m=3} &= \frac{9}{112}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
f_\perp\big|_{x=\frac{1}{3}} &= - \frac{23}{108} \\
\frac{\partial^1 (f_\perp)}{\partial x^1} &= 3 x^{2} \\
\left.\frac{\partial^1 (f_\perp)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{1}{3}
\end{aligned}
$$
:::

:::example[共鳴を引き算で避ける右辺（4）]
$e$ の方向を定数とする。右辺に $x^m$ を選び、定数成分を引いて平均零にする。

$$
m=4
$$

単位区間で平均を計算する。

$$
\begin{aligned}
\langle x^m,1\rangle\big|_{m=4} &= \frac{1}{5}
\end{aligned}
$$

平均を除けば可解条件を満たす。

$$
\begin{aligned}
f_\perp\big|_{m=4} &= x^{4} - \frac{1}{5}
\end{aligned}
$$

右辺のノルムを求め、引いた成分のノルムを除く。

$$
\begin{aligned}
\|f_\perp\|^2\big|_{m=4} &= \frac{16}{225}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
f_\perp\big|_{x=\frac{1}{3}} &= - \frac{76}{405} \\
\frac{\partial^1 (f_\perp)}{\partial x^1} &= 4 x^{3} \\
\left.\frac{\partial^1 (f_\perp)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{4}{27}
\end{aligned}
$$
:::

## Hilbert–Schmidt 理論

### 二乗可積分な核を近似する

$K\in L^2([a,b]^2)$ なら $Tu(x)=\int K(x,t)u(t)\dd t$ は Cauchy–Schwarz により $\|Tu\|_2\le\|K\|_2\|u\|_2$ を満たす。矩形ごとの定数関数は $L^2$ で稠密なので有限分離核で近似でき、作用素はコンパクトである。正規直交基底 $e_j$ に対して $Te_j(x)=\inner{K(x,\cdot)}{\overline{e_j}}$ だから、共役した基底で Parseval を適用し $x$ で積分すると $\sum_j\|Te_j\|^2=\|K\|_2^2$ となる。これを :term[ヒルベルト・シュミットノルム]{en="Hilbert-Schmidt norm" yomi="ひるべるとしゅみっとのるむ"} の二乗という。Fubini で随伴の核は $\overline{K(t,x)}$ である。ノルム恒等式により、自己共役であることは $K(x,t)=\overline{K(t,x)}$ がほとんど至る所で成り立つことと同値である。

自己共役の場合、[コンパクト作用素](/math/functional-analysis/compact-operators/) のスペクトル定理で、非零固有値 $\mu_j$ の基底 $\phi_j$ を取ると $K=\sum \mu_j\phi_j(x)\overline{\phi_j(t)}$ が $L^2$ で収束する。実際、固有基底に零空間の基底を足した場合の上のノルム恒等式が $\sum\mu_j^2=\|K\|_2^2$ を与える。例えば $K(x,t)=xt$ は $\phi=\sqrt3 x$、$\mu=1/3$ の一項だけである。非斉次項のこの成分は $1-\lambda/3$ で割り、残りはそのままになる。

### 対称でない核を特異値で展開する

$T^*T$ はコンパクトで正である。正の固有値を $\sigma_j^2$、対応する正規直交基底を $v_j$ とし、$u_j=Tv_j/\sigma_j$ と置く。$\inner{Tv_i}{Tv_j}=\inner{v_i}{T^*Tv_j}$ によって $u_j$ は正規直交する。$T^*T$ の零空間は $T$ の零空間と一致するので、コンパクト自己共役作用素の完全性から
$$
Tf=\sum_j\sigma_j\inner{f}{v_j}u_j
$$
である。Hilbert–Schmidt の恒等式から $\sum\sigma_j^2=\|K\|_2^2$ であり、対応する核 $\sum_j\sigma_j u_j(x)\overline{v_j(t)}$ は $L^2$ で収束する。例えば $[0,1]$ の $K(x,t)=x$ では右の $v=1$、左の $u=\sqrt3 x$、$\sigma=1/\sqrt3$ である。固有値 $1/2$ と特異値 $1/\sqrt3$ は異なり、非対称核では二つを区別する必要がある。

### 連続な正値核では一様収束する

連続な Hermite 核 $K$ が正の作用素を定めるとする。正固有値 $\mu_j$ に対して $\phi_j=\mu_j^{-1}T\phi_j$ とすれば $\phi_j$ は連続に取れる。このとき $K$ の固有関数展開は実際には正方形全体で絶対一様収束し、$\int K(x,x)\dd x=\sum_j\mu_j$ が成り立つ。これを :term[マーサーの定理]{en="Mercer's theorem" yomi="まーさーのていり"} という。

有限部分和を引いた残りの核を $R_N$ とする。$R_N$ と $K-R_N$ は正の作用素を定める。点ごとの有限行列 $[R_N(x_i,x_j)]$ も半正定値である。実際、各点の近くで積分が $1$ の非負試験関数を取り、その線形結合に作用素の正性を使い、台を点へ縮めればよい。同様に $[K(x_i,x_j)-R_N(x_i,x_j)]$ も半正定値である。

$C=\max_x K(x,x)$ とすれば $0\le R_N(x,x)\le C$ である。三点での行列の Cauchy–Schwarz 不等式を、$x$ と $x'$ を引いたベクトルと、$y$ のベクトルに適用すると
$$
|R_N(x,y)-R_N(x',y)|^2
\le C[K(x,x)+K(x',x')-2\Re K(x,x')]
$$
を得る。括弧の式は $x'\to x$ で一様に零へ近づく。Hermite 対称性から第二変数についても同じ評価があり、$R_N$ は一様有界かつ同程度連続である。すでに $R_N\to0$ が $L^2$ で分かっているので、Arzelà–Ascoli の定理を適用すると一様にも零へ収束する。さらに固有関数級数の $N$ 項目以降の絶対値の和は、数列の Cauchy–Schwarz によって $\sqrt{R_N(x,x)R_N(y,y)}$ 以下になる。ここでは対角の級数が $K(x,x)-R_N(x,x)$ に収束することを使った。最後に対角の一様収束を積分し、各 $\phi_j$ のノルムが $1$ であることを用いればトレースの式が従う。

:::example[折れ線の核の固有値と対角積分]
$[0,1]$ で $K(x,t)=\min(x,t)$ とする。$K(x,t)=\int_0^1\mathbf1_{s<x}\mathbf1_{s<t}\dd s$ なので作用素は正である。$v=Tf$ と書くと $v(0)=v'(1)=0$、$-v''=f$ である。従って $T\phi=\mu\phi$ の非零固有対は $-\phi''=\mu^{-1}\phi$ と混合端条件を満たし、$\phi_n=\sqrt2\sin((n+1/2)\pi x)$、$\mu_n=1/[(n+1/2)^2\pi^2]$（$n\ge0$）となる。逆にこの二階方程式の解に $T$ を掛ければ、同じ端条件から固有対になる。対角 $K(x,x)=x$ の積分は $1/2$ なので、$\sum_{n\ge0}1/(2n+1)^2=\pi^2/8$ を再現する。
:::

### 核のノルムと階数一の特異値

$K(x,y)=p(x)q(y)$ の $L^2$ ノルムは、積を分けると $\|p\|\|q\|$ になる。同時に Cauchy–Schwarz と入力 $\overline q/\|q\|$ を使えば、作用素ノルムも同じ値である。複数の直交方向へ別々に送る場合は、Hilbert–Schmidt ノルムの二乗が特異値の二乗和、作用素ノルムが最大特異値になる。二つを区別すると、小さい固有値からの寄与を数えられる。

:::example[非対称な核のノルム（1）]
単位区間で $p=x^m,q=x$ を取る。

$$
m=1
$$

二つの二乗ノルムの積を求める。

$$
\begin{aligned}
\|K\|_{\rm HS}^2\big|_{m=1} &= \frac{1}{9}
\end{aligned}
$$

正規化した入力 $q/\|q\|$ への応答を計算する。

$$
\begin{aligned}
T(q/\|q\|)\big|_{m=1} &= \frac{\sqrt{3} x}{3}
\end{aligned}
$$

固有値と特異値は一般に異なる。固有関数 $p$ へ作用させた倍率を求める。

$$
\begin{aligned}
\eta_{\rm eig}\big|_{m=1} &= \frac{1}{3}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
T(q/\|q\|)\big|_{x=\frac{1}{3}} &= \frac{\sqrt{3}}{9} \\
\frac{\partial^1 (T(q/\|q\|))}{\partial x^1} &= \frac{\sqrt{3}}{3} \\
\left.\frac{\partial^1 (T(q/\|q\|))}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{\sqrt{3}}{3}
\end{aligned}
$$
:::

:::example[非対称な核のノルム（2）]
単位区間で $p=x^m,q=x$ を取る。

$$
m=2
$$

二つの二乗ノルムの積を求める。

$$
\begin{aligned}
\|K\|_{\rm HS}^2\big|_{m=2} &= \frac{1}{15}
\end{aligned}
$$

正規化した入力 $q/\|q\|$ への応答を計算する。

$$
\begin{aligned}
T(q/\|q\|)\big|_{m=2} &= \frac{\sqrt{3} x^{2}}{3}
\end{aligned}
$$

固有値と特異値は一般に異なる。固有関数 $p$ へ作用させた倍率を求める。

$$
\begin{aligned}
\eta_{\rm eig}\big|_{m=2} &= \frac{1}{4}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
T(q/\|q\|)\big|_{x=\frac{1}{3}} &= \frac{\sqrt{3}}{27} \\
\frac{\partial^1 (T(q/\|q\|))}{\partial x^1} &= \frac{2 \sqrt{3} x}{3} \\
\left.\frac{\partial^1 (T(q/\|q\|))}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{2 \sqrt{3}}{9}
\end{aligned}
$$
:::

:::example[非対称な核のノルム（3）]
単位区間で $p=x^m,q=x$ を取る。

$$
m=3
$$

二つの二乗ノルムの積を求める。

$$
\begin{aligned}
\|K\|_{\rm HS}^2\big|_{m=3} &= \frac{1}{21}
\end{aligned}
$$

正規化した入力 $q/\|q\|$ への応答を計算する。

$$
\begin{aligned}
T(q/\|q\|)\big|_{m=3} &= \frac{\sqrt{3} x^{3}}{3}
\end{aligned}
$$

固有値と特異値は一般に異なる。固有関数 $p$ へ作用させた倍率を求める。

$$
\begin{aligned}
\eta_{\rm eig}\big|_{m=3} &= \frac{1}{5}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
T(q/\|q\|)\big|_{x=\frac{1}{3}} &= \frac{\sqrt{3}}{81} \\
\frac{\partial^1 (T(q/\|q\|))}{\partial x^1} &= \sqrt{3} x^{2} \\
\left.\frac{\partial^1 (T(q/\|q\|))}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{\sqrt{3}}{9}
\end{aligned}
$$
:::

:::example[非対称な核のノルム（4）]
単位区間で $p=x^m,q=x$ を取る。

$$
m=4
$$

二つの二乗ノルムの積を求める。

$$
\begin{aligned}
\|K\|_{\rm HS}^2\big|_{m=4} &= \frac{1}{27}
\end{aligned}
$$

正規化した入力 $q/\|q\|$ への応答を計算する。

$$
\begin{aligned}
T(q/\|q\|)\big|_{m=4} &= \frac{\sqrt{3} x^{4}}{3}
\end{aligned}
$$

固有値と特異値は一般に異なる。固有関数 $p$ へ作用させた倍率を求める。

$$
\begin{aligned}
\eta_{\rm eig}\big|_{m=4} &= \frac{1}{6}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
T(q/\|q\|)\big|_{x=\frac{1}{3}} &= \frac{\sqrt{3}}{243} \\
\frac{\partial^1 (T(q/\|q\|))}{\partial x^1} &= \frac{4 \sqrt{3} x^{3}}{3} \\
\left.\frac{\partial^1 (T(q/\|q\|))}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{4 \sqrt{3}}{81}
\end{aligned}
$$
:::

## Abel の積分方程式

### 弱い特異性を二度積分する

任意の $\alpha>0$ に対して、有限区間 $[0,b]$ 上の分数積分を $I^\alpha u(x)=\Gamma(\alpha)^{-1}\int_0^x (x-t)^{\alpha-1}u(t)\dd t$ とする。$u\in L^1$ なら、絶対値を入れて積分順序を交換すると $\|I^\alpha u\|_1\le b^\alpha\|u\|_1/\Gamma(1+\alpha)$ を得る。さらに $\beta>0$ に対し、合成の二重積分の内側は
$$
\int_t^x (x-s)^{\beta-1}(s-t)^{\alpha-1}\dd s
=(x-t)^{\alpha+\beta-1}B(\alpha,\beta)
$$
である。絶対値でも $x$ について積分すれば有限なので、Fubini が適用でき $I^\beta I^\alpha=I^{\alpha+\beta}$ が $L^1$ で成り立つ。

$0<\alpha<1$ の Abel 方程式 $I^\alpha u=f$ を $L^1$ で解く。$I^{1-\alpha}$ を適用すれば $I^1u$ となるので、$I^{1-\alpha}f$ が絶対連続で原点値が零の代表を持つことが必要である。この条件のもとでその代表の導関数を $u$ とすれば $I^{1-\alpha}(I^\alpha u-f)=0$ となる。$I^\alpha$ をさらに適用して一階微分すれば $I^\alpha u=f$ がほとんど至る所で成り立つ。同じ単射性から解は一意である。連続解を求める場合は、導関数が連続な代表を持つことが追加条件になる。

:::example[Abel 核で定数を作る]
$I^{1/2}u=1$ に対し $I^{1/2}1=2\sqrt{x/\pi}$ である。微分すれば $u(x)=1/\sqrt{\pi x}$ を得る。この解は原点で発散するが $L^1$ に属する。逆に連続な $u$ なら $|I^{1/2}u(x)|\le2\sqrt{x/\pi}\|u\|_\infty$ なので定数 $1$ にはならない。元の方程式へ特異解を戻すと、$t=xv$ で積分が $B(1/2,1/2)/\pi=1$ になる。
:::

### 半階積分を二度合成する

$I^\alpha I^\beta$ の二重積分で順序を交換し、内側の変数を $s=t+(x-t)v$ と置く。内側の積分は $(x-t)^{\alpha+\beta-1}B(\alpha,\beta)$ となる。外側にある二つの Gamma で割ると $I^{\alpha+\beta}$ の規格化になる。特に $\alpha=\beta=1/2$ では普通の積分に戻る。端の可積分性には $\alpha,\beta>0$ が必要である。

弱特異な Volterra 方程式 $u-\lambda I^\alpha u=f$ でも、合成則から $(I^\alpha)^n=I^{n\alpha}$ となる。連続関数に対して $\|I^{n\alpha}\|_\infty\le b^{n\alpha}/\Gamma(1+n\alpha)$ なので、Stirling の公式から全ての $\lambda$ で Neumann 級数が収束する。例えば $f=1$ のとき $u(x)=\sum_{n\ge0}\lambda^n x^{n\alpha}/\Gamma(1+n\alpha)$ である。$\alpha=1$ なら $e^{\lambda x}$ に戻る。

:::example[半階積分の冪への作用（1）]
入力を $x^m$ とし、原点からの半階積分を求める。

$$
m=1
$$

変数を $t=xv$ として冪の因子と係数を分ける。

$$
\begin{aligned}
I^{1/2}x^m\big|_{m=1} &= \frac{4 x^{\frac{3}{2}}}{3 \sqrt{\pi}}
\end{aligned}
$$

半階をもう一回加えると通常の積分である。

$$
\begin{aligned}
I^{1/2}I^{1/2}x^m\big|_{m=1} &= \frac{x^{2}}{2}
\end{aligned}
$$

最後に通常の微分を行うと入力へ戻る。

$$
\begin{aligned}
\partial_x I^{1/2}I^{1/2}x^m\big|_{m=1} &= x
\end{aligned}
$$

観測点 $x=2$ での値と導関数を計算する。

$$
\begin{aligned}
I^{1/2}x^m\big|_{x=2} &= \frac{8 \sqrt{2}}{3 \sqrt{\pi}} \\
\frac{\partial^1 (I^{1/2}x^m)}{\partial x^1} &= \frac{2 \sqrt{x}}{\sqrt{\pi}} \\
\left.\frac{\partial^1 (I^{1/2}x^m)}{\partial x^1}\right|_{x=2} &= \frac{2 \sqrt{2}}{\sqrt{\pi}} \\
\frac{\partial^2 (I^{1/2}x^m)}{\partial x^2} &= \frac{1}{\sqrt{\pi} \sqrt{x}} \\
\left.\frac{\partial^2 (I^{1/2}x^m)}{\partial x^2}\right|_{x=2} &= \frac{\sqrt{2}}{2 \sqrt{\pi}}
\end{aligned}
$$

観測点 $x=2$ での値と導関数を計算する。

$$
\begin{aligned}
I^{1/2}I^{1/2}x^m\big|_{x=2} &= 2 \\
\frac{\partial^1 (I^{1/2}I^{1/2}x^m)}{\partial x^1} &= x \\
\left.\frac{\partial^1 (I^{1/2}I^{1/2}x^m)}{\partial x^1}\right|_{x=2} &= 2 \\
\frac{\partial^2 (I^{1/2}I^{1/2}x^m)}{\partial x^2} &= 1 \\
\left.\frac{\partial^2 (I^{1/2}I^{1/2}x^m)}{\partial x^2}\right|_{x=2} &= 1
\end{aligned}
$$

観測点 $x=2$ での値と導関数を計算する。

$$
\begin{aligned}
\partial_x I^{1/2}I^{1/2}x^m\big|_{x=2} &= 2 \\
\frac{\partial^1 (\partial_x I^{1/2}I^{1/2}x^m)}{\partial x^1} &= 1 \\
\left.\frac{\partial^1 (\partial_x I^{1/2}I^{1/2}x^m)}{\partial x^1}\right|_{x=2} &= 1 \\
\frac{\partial^2 (\partial_x I^{1/2}I^{1/2}x^m)}{\partial x^2} &= 0 \\
\left.\frac{\partial^2 (\partial_x I^{1/2}I^{1/2}x^m)}{\partial x^2}\right|_{x=2} &= 0
\end{aligned}
$$
:::

:::example[半階積分の冪への作用（2）]
入力を $x^m$ とし、原点からの半階積分を求める。

$$
m=2
$$

変数を $t=xv$ として冪の因子と係数を分ける。

$$
\begin{aligned}
I^{1/2}x^m\big|_{m=2} &= \frac{16 x^{\frac{5}{2}}}{15 \sqrt{\pi}}
\end{aligned}
$$

半階をもう一回加えると通常の積分である。

$$
\begin{aligned}
I^{1/2}I^{1/2}x^m\big|_{m=2} &= \frac{x^{3}}{3}
\end{aligned}
$$

最後に通常の微分を行うと入力へ戻る。

$$
\begin{aligned}
\partial_x I^{1/2}I^{1/2}x^m\big|_{m=2} &= x^{2}
\end{aligned}
$$

観測点 $x=2$ での値と導関数を計算する。

$$
\begin{aligned}
I^{1/2}x^m\big|_{x=2} &= \frac{64 \sqrt{2}}{15 \sqrt{\pi}} \\
\frac{\partial^1 (I^{1/2}x^m)}{\partial x^1} &= \frac{8 x^{\frac{3}{2}}}{3 \sqrt{\pi}} \\
\left.\frac{\partial^1 (I^{1/2}x^m)}{\partial x^1}\right|_{x=2} &= \frac{16 \sqrt{2}}{3 \sqrt{\pi}} \\
\frac{\partial^2 (I^{1/2}x^m)}{\partial x^2} &= \frac{4 \sqrt{x}}{\sqrt{\pi}} \\
\left.\frac{\partial^2 (I^{1/2}x^m)}{\partial x^2}\right|_{x=2} &= \frac{4 \sqrt{2}}{\sqrt{\pi}}
\end{aligned}
$$

観測点 $x=2$ での値と導関数を計算する。

$$
\begin{aligned}
I^{1/2}I^{1/2}x^m\big|_{x=2} &= \frac{8}{3} \\
\frac{\partial^1 (I^{1/2}I^{1/2}x^m)}{\partial x^1} &= x^{2} \\
\left.\frac{\partial^1 (I^{1/2}I^{1/2}x^m)}{\partial x^1}\right|_{x=2} &= 4 \\
\frac{\partial^2 (I^{1/2}I^{1/2}x^m)}{\partial x^2} &= 2 x \\
\left.\frac{\partial^2 (I^{1/2}I^{1/2}x^m)}{\partial x^2}\right|_{x=2} &= 4
\end{aligned}
$$

観測点 $x=2$ での値と導関数を計算する。

$$
\begin{aligned}
\partial_x I^{1/2}I^{1/2}x^m\big|_{x=2} &= 4 \\
\frac{\partial^1 (\partial_x I^{1/2}I^{1/2}x^m)}{\partial x^1} &= 2 x \\
\left.\frac{\partial^1 (\partial_x I^{1/2}I^{1/2}x^m)}{\partial x^1}\right|_{x=2} &= 4 \\
\frac{\partial^2 (\partial_x I^{1/2}I^{1/2}x^m)}{\partial x^2} &= 2 \\
\left.\frac{\partial^2 (\partial_x I^{1/2}I^{1/2}x^m)}{\partial x^2}\right|_{x=2} &= 2
\end{aligned}
$$
:::

:::example[半階積分の冪への作用（3）]
入力を $x^m$ とし、原点からの半階積分を求める。

$$
m=3
$$

変数を $t=xv$ として冪の因子と係数を分ける。

$$
\begin{aligned}
I^{1/2}x^m\big|_{m=3} &= \frac{32 x^{\frac{7}{2}}}{35 \sqrt{\pi}}
\end{aligned}
$$

半階をもう一回加えると通常の積分である。

$$
\begin{aligned}
I^{1/2}I^{1/2}x^m\big|_{m=3} &= \frac{x^{4}}{4}
\end{aligned}
$$

最後に通常の微分を行うと入力へ戻る。

$$
\begin{aligned}
\partial_x I^{1/2}I^{1/2}x^m\big|_{m=3} &= x^{3}
\end{aligned}
$$

観測点 $x=2$ での値と導関数を計算する。

$$
\begin{aligned}
I^{1/2}x^m\big|_{x=2} &= \frac{256 \sqrt{2}}{35 \sqrt{\pi}} \\
\frac{\partial^1 (I^{1/2}x^m)}{\partial x^1} &= \frac{16 x^{\frac{5}{2}}}{5 \sqrt{\pi}} \\
\left.\frac{\partial^1 (I^{1/2}x^m)}{\partial x^1}\right|_{x=2} &= \frac{64 \sqrt{2}}{5 \sqrt{\pi}} \\
\frac{\partial^2 (I^{1/2}x^m)}{\partial x^2} &= \frac{8 x^{\frac{3}{2}}}{\sqrt{\pi}} \\
\left.\frac{\partial^2 (I^{1/2}x^m)}{\partial x^2}\right|_{x=2} &= \frac{16 \sqrt{2}}{\sqrt{\pi}}
\end{aligned}
$$

観測点 $x=2$ での値と導関数を計算する。

$$
\begin{aligned}
I^{1/2}I^{1/2}x^m\big|_{x=2} &= 4 \\
\frac{\partial^1 (I^{1/2}I^{1/2}x^m)}{\partial x^1} &= x^{3} \\
\left.\frac{\partial^1 (I^{1/2}I^{1/2}x^m)}{\partial x^1}\right|_{x=2} &= 8 \\
\frac{\partial^2 (I^{1/2}I^{1/2}x^m)}{\partial x^2} &= 3 x^{2} \\
\left.\frac{\partial^2 (I^{1/2}I^{1/2}x^m)}{\partial x^2}\right|_{x=2} &= 12
\end{aligned}
$$

観測点 $x=2$ での値と導関数を計算する。

$$
\begin{aligned}
\partial_x I^{1/2}I^{1/2}x^m\big|_{x=2} &= 8 \\
\frac{\partial^1 (\partial_x I^{1/2}I^{1/2}x^m)}{\partial x^1} &= 3 x^{2} \\
\left.\frac{\partial^1 (\partial_x I^{1/2}I^{1/2}x^m)}{\partial x^1}\right|_{x=2} &= 12 \\
\frac{\partial^2 (\partial_x I^{1/2}I^{1/2}x^m)}{\partial x^2} &= 6 x \\
\left.\frac{\partial^2 (\partial_x I^{1/2}I^{1/2}x^m)}{\partial x^2}\right|_{x=2} &= 12
\end{aligned}
$$
:::

:::example[半階積分の冪への作用（4）]
入力を $x^m$ とし、原点からの半階積分を求める。

$$
m=4
$$

変数を $t=xv$ として冪の因子と係数を分ける。

$$
\begin{aligned}
I^{1/2}x^m\big|_{m=4} &= \frac{256 x^{\frac{9}{2}}}{315 \sqrt{\pi}}
\end{aligned}
$$

半階をもう一回加えると通常の積分である。

$$
\begin{aligned}
I^{1/2}I^{1/2}x^m\big|_{m=4} &= \frac{x^{5}}{5}
\end{aligned}
$$

最後に通常の微分を行うと入力へ戻る。

$$
\begin{aligned}
\partial_x I^{1/2}I^{1/2}x^m\big|_{m=4} &= x^{4}
\end{aligned}
$$

観測点 $x=2$ での値と導関数を計算する。

$$
\begin{aligned}
I^{1/2}x^m\big|_{x=2} &= \frac{4096 \sqrt{2}}{315 \sqrt{\pi}} \\
\frac{\partial^1 (I^{1/2}x^m)}{\partial x^1} &= \frac{128 x^{\frac{7}{2}}}{35 \sqrt{\pi}} \\
\left.\frac{\partial^1 (I^{1/2}x^m)}{\partial x^1}\right|_{x=2} &= \frac{1024 \sqrt{2}}{35 \sqrt{\pi}} \\
\frac{\partial^2 (I^{1/2}x^m)}{\partial x^2} &= \frac{64 x^{\frac{5}{2}}}{5 \sqrt{\pi}} \\
\left.\frac{\partial^2 (I^{1/2}x^m)}{\partial x^2}\right|_{x=2} &= \frac{256 \sqrt{2}}{5 \sqrt{\pi}}
\end{aligned}
$$

観測点 $x=2$ での値と導関数を計算する。

$$
\begin{aligned}
I^{1/2}I^{1/2}x^m\big|_{x=2} &= \frac{32}{5} \\
\frac{\partial^1 (I^{1/2}I^{1/2}x^m)}{\partial x^1} &= x^{4} \\
\left.\frac{\partial^1 (I^{1/2}I^{1/2}x^m)}{\partial x^1}\right|_{x=2} &= 16 \\
\frac{\partial^2 (I^{1/2}I^{1/2}x^m)}{\partial x^2} &= 4 x^{3} \\
\left.\frac{\partial^2 (I^{1/2}I^{1/2}x^m)}{\partial x^2}\right|_{x=2} &= 32
\end{aligned}
$$

観測点 $x=2$ での値と導関数を計算する。

$$
\begin{aligned}
\partial_x I^{1/2}I^{1/2}x^m\big|_{x=2} &= 16 \\
\frac{\partial^1 (\partial_x I^{1/2}I^{1/2}x^m)}{\partial x^1} &= 4 x^{3} \\
\left.\frac{\partial^1 (\partial_x I^{1/2}I^{1/2}x^m)}{\partial x^1}\right|_{x=2} &= 32 \\
\frac{\partial^2 (\partial_x I^{1/2}I^{1/2}x^m)}{\partial x^2} &= 12 x^{2} \\
\left.\frac{\partial^2 (\partial_x I^{1/2}I^{1/2}x^m)}{\partial x^2}\right|_{x=2} &= 48
\end{aligned}
$$
:::

## 微分方程式から積分方程式への変換

### 初期値を積分の自由項に置く

$y'=F(x,y)$、$y(a)=y_0$ の解は $y(x)=y_0+\int_a^xF(t,y(t))\dd t$ と同値である。連続な右辺を微分すれば逆向きも成り立つ。$F$ が $y$ で Lipschitz なら、十分短い区間で積分写像の Lipschitz 定数は $L(b-a)<1$ であり、反復の差は等比級数で収束する。区間を広げるには解が係数の定義域にとどまることを確認する。二階式 $y''+q(x)y=g(x)$ では
$$
y(x)=y_0+(x-a)v_0+\int_a^x(x-t)(g(t)-q(t)y(t))\dd t
$$
となる。二回積分し三角領域の順序を交換すれば得られる。初期値問題の存在を示すには、初期点を含む長方形で $|F|\le M$、$y$ に関する Lipschitz 定数を $L$ とする。時刻幅 $h$ は長方形の内部に収め、$Mh$ が長方形の $y$ 方向の半幅以下で、$Lh<1$ となるように選ぶ。これにより積分写像は閉球を自分に移し、その上で縮小となる。

一方、境界値問題 $(L_0+V)u=f$ では、境界条件込みで $L_0$ の逆が存在すれば、その Green 核を用いて $u=G_0 f-G_0Vu$ という Fredholm 方程式に移る。$G_0$ が $L_0$ の定義域に写すので、積分方程式から元の問題へも戻れる。初期値問題と違い、共鳴で逆が存在しない場合がある。

:::example[同じ方程式を二端で固定する]
$-y''-\lambda y=f$、$y(0)=y(1)=0$ で、$K(x,t)=\min(x,t)[1-\max(x,t)]$ を基準の逆核に取る。すると $y-\lambda Ty=Tf$ となる。$f=\sin\pi x$ なら $Tf=\pi^{-2}f$ である。$\lambda$ が $(n\pi)^2$ のどれにも等しくなければ $y=(\pi^2-\lambda)^{-1}f$ が一意解である。$\lambda=\pi^2$ なら $Tf$ は斉次固有関数と直交せず、交代定理により解はない。$\lambda=(n\pi)^2$、$n\ge2$ なら上の特解に $\sin(n\pi x)$ の任意の倍数を加えたものが全解である。
:::

### 初期傾斜を積分方程式へ保存する

$y^{\prime\prime}=F(x,y)$ の積分形は $y(x)=y(0)+xy^\prime(0)+\int_0^x(x-s)F(s,y(s))\dd s$ である。積分を二回微分すれば元の式になり、$x=0$ で積分とその一階微分がともに零になる。これで形式的な積分だけでなく、初期値との一対一の対応も示せる。外力を与える場合は、既知関数だけの積分を自由項に含められる。

:::example[振動方程式の最初の反復（1）]
初期値を $(0,1)$、方程式を $y^{\prime\prime}=-m^2 y$ とする。

$$
m=1
$$

最初の自由項は線形である。

$$
\begin{aligned}
y_0 &= x
\end{aligned}
$$

一度反復すると三次の補正が出る。

$$
\begin{aligned}
y_1\big|_{m=1} &= - \frac{x^{3}}{6} + x
\end{aligned}
$$

反復をさらに二回行うと五次と七次が加わる。

$$
\begin{aligned}
y_3\big|_{m=1} &= - \frac{x^{7}}{5040} + \frac{x^{5}}{120} - \frac{x^{3}}{6} + x
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
y_0\big|_{x=\frac{1}{3}} &= \frac{1}{3} \\
\frac{\partial^1 (y_0)}{\partial x^1} &= 1 \\
\left.\frac{\partial^1 (y_0)}{\partial x^1}\right|_{x=\frac{1}{3}} &= 1 \\
\frac{\partial^2 (y_0)}{\partial x^2} &= 0 \\
\left.\frac{\partial^2 (y_0)}{\partial x^2}\right|_{x=\frac{1}{3}} &= 0
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
y_1\big|_{x=\frac{1}{3}} &= \frac{53}{162} \\
\frac{\partial^1 (y_1)}{\partial x^1} &= 1 - \frac{x^{2}}{2} \\
\left.\frac{\partial^1 (y_1)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{17}{18} \\
\frac{\partial^2 (y_1)}{\partial x^2} &= - x \\
\left.\frac{\partial^2 (y_1)}{\partial x^2}\right|_{x=\frac{1}{3}} &= - \frac{1}{3}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
y_3\big|_{x=\frac{1}{3}} &= \frac{3606497}{11022480} \\
\frac{\partial^1 (y_3)}{\partial x^1} &= - \frac{x^{6}}{720} + \frac{x^{4}}{24} - \frac{x^{2}}{2} + 1 \\
\left.\frac{\partial^1 (y_3)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{495989}{524880} \\
\frac{\partial^2 (y_3)}{\partial x^2} &= - \frac{x^{5}}{120} + \frac{x^{3}}{6} - x \\
\left.\frac{\partial^2 (y_3)}{\partial x^2}\right|_{x=\frac{1}{3}} &= - \frac{9541}{29160}
\end{aligned}
$$
:::

:::example[振動方程式の最初の反復（2）]
初期値を $(0,1)$、方程式を $y^{\prime\prime}=-m^2 y$ とする。

$$
m=2
$$

最初の自由項は線形である。

$$
\begin{aligned}
y_0 &= x
\end{aligned}
$$

一度反復すると三次の補正が出る。

$$
\begin{aligned}
y_1\big|_{m=2} &= - \frac{2 x^{3}}{3} + x
\end{aligned}
$$

反復をさらに二回行うと五次と七次が加わる。

$$
\begin{aligned}
y_3\big|_{m=2} &= \frac{x \left(- 4 x^{6} + 42 x^{4} - 210 x^{2} + 315\right)}{315}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
y_0\big|_{x=\frac{1}{3}} &= \frac{1}{3} \\
\frac{\partial^1 (y_0)}{\partial x^1} &= 1 \\
\left.\frac{\partial^1 (y_0)}{\partial x^1}\right|_{x=\frac{1}{3}} &= 1 \\
\frac{\partial^2 (y_0)}{\partial x^2} &= 0 \\
\left.\frac{\partial^2 (y_0)}{\partial x^2}\right|_{x=\frac{1}{3}} &= 0
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
y_1\big|_{x=\frac{1}{3}} &= \frac{25}{81} \\
\frac{\partial^1 (y_1)}{\partial x^1} &= 1 - 2 x^{2} \\
\left.\frac{\partial^1 (y_1)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{7}{9} \\
\frac{\partial^2 (y_1)}{\partial x^2} &= - 4 x \\
\left.\frac{\partial^2 (y_1)}{\partial x^2}\right|_{x=\frac{1}{3}} &= - \frac{4}{3}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
y_3\big|_{x=\frac{1}{3}} &= \frac{212999}{688905} \\
\frac{\partial^1 (y_3)}{\partial x^1} &= - \frac{4 x^{6}}{45} + \frac{2 x^{4}}{3} - 2 x^{2} + 1 \\
\left.\frac{\partial^1 (y_3)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{25781}{32805} \\
\frac{\partial^2 (y_3)}{\partial x^2} &= \frac{4 x \left(- 2 x^{4} + 10 x^{2} - 15\right)}{15} \\
\left.\frac{\partial^2 (y_3)}{\partial x^2}\right|_{x=\frac{1}{3}} &= - \frac{4508}{3645}
\end{aligned}
$$
:::

:::example[振動方程式の最初の反復（3）]
初期値を $(0,1)$、方程式を $y^{\prime\prime}=-m^2 y$ とする。

$$
m=3
$$

最初の自由項は線形である。

$$
\begin{aligned}
y_0 &= x
\end{aligned}
$$

一度反復すると三次の補正が出る。

$$
\begin{aligned}
y_1\big|_{m=3} &= - \frac{3 x^{3}}{2} + x
\end{aligned}
$$

反復をさらに二回行うと五次と七次が加わる。

$$
\begin{aligned}
y_3\big|_{m=3} &= \frac{x \left(- 81 x^{6} + 378 x^{4} - 840 x^{2} + 560\right)}{560}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
y_0\big|_{x=\frac{1}{3}} &= \frac{1}{3} \\
\frac{\partial^1 (y_0)}{\partial x^1} &= 1 \\
\left.\frac{\partial^1 (y_0)}{\partial x^1}\right|_{x=\frac{1}{3}} &= 1 \\
\frac{\partial^2 (y_0)}{\partial x^2} &= 0 \\
\left.\frac{\partial^2 (y_0)}{\partial x^2}\right|_{x=\frac{1}{3}} &= 0
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
y_1\big|_{x=\frac{1}{3}} &= \frac{5}{18} \\
\frac{\partial^1 (y_1)}{\partial x^1} &= 1 - \frac{9 x^{2}}{2} \\
\left.\frac{\partial^1 (y_1)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{1}{2} \\
\frac{\partial^2 (y_1)}{\partial x^2} &= - 9 x \\
\left.\frac{\partial^2 (y_1)}{\partial x^2}\right|_{x=\frac{1}{3}} &= -3
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
y_3\big|_{x=\frac{1}{3}} &= \frac{4241}{15120} \\
\frac{\partial^1 (y_3)}{\partial x^1} &= - \frac{81 x^{6}}{80} + \frac{27 x^{4}}{8} - \frac{9 x^{2}}{2} + 1 \\
\left.\frac{\partial^1 (y_3)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{389}{720} \\
\frac{\partial^2 (y_3)}{\partial x^2} &= \frac{9 x \left(- 27 x^{4} + 60 x^{2} - 40\right)}{40} \\
\left.\frac{\partial^2 (y_3)}{\partial x^2}\right|_{x=\frac{1}{3}} &= - \frac{101}{40}
\end{aligned}
$$
:::

:::example[振動方程式の最初の反復（4）]
初期値を $(0,1)$、方程式を $y^{\prime\prime}=-m^2 y$ とする。

$$
m=4
$$

最初の自由項は線形である。

$$
\begin{aligned}
y_0 &= x
\end{aligned}
$$

一度反復すると三次の補正が出る。

$$
\begin{aligned}
y_1\big|_{m=4} &= - \frac{8 x^{3}}{3} + x
\end{aligned}
$$

反復をさらに二回行うと五次と七次が加わる。

$$
\begin{aligned}
y_3\big|_{m=4} &= \frac{x \left(- 256 x^{6} + 672 x^{4} - 840 x^{2} + 315\right)}{315}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
y_0\big|_{x=\frac{1}{3}} &= \frac{1}{3} \\
\frac{\partial^1 (y_0)}{\partial x^1} &= 1 \\
\left.\frac{\partial^1 (y_0)}{\partial x^1}\right|_{x=\frac{1}{3}} &= 1 \\
\frac{\partial^2 (y_0)}{\partial x^2} &= 0 \\
\left.\frac{\partial^2 (y_0)}{\partial x^2}\right|_{x=\frac{1}{3}} &= 0
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
y_1\big|_{x=\frac{1}{3}} &= \frac{19}{81} \\
\frac{\partial^1 (y_1)}{\partial x^1} &= 1 - 8 x^{2} \\
\left.\frac{\partial^1 (y_1)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{1}{9} \\
\frac{\partial^2 (y_1)}{\partial x^2} &= - 16 x \\
\left.\frac{\partial^2 (y_1)}{\partial x^2}\right|_{x=\frac{1}{3}} &= - \frac{16}{3}
\end{aligned}
$$

観測点 $x=\frac{1}{3}$ での値と導関数を計算する。

$$
\begin{aligned}
y_3\big|_{x=\frac{1}{3}} &= \frac{167387}{688905} \\
\frac{\partial^1 (y_3)}{\partial x^1} &= - \frac{256 x^{6}}{45} + \frac{32 x^{4}}{3} - 8 x^{2} + 1 \\
\left.\frac{\partial^1 (y_3)}{\partial x^1}\right|_{x=\frac{1}{3}} &= \frac{7709}{32805} \\
\frac{\partial^2 (y_3)}{\partial x^2} &= \frac{16 x \left(- 32 x^{4} + 40 x^{2} - 15\right)}{15} \\
\left.\frac{\partial^2 (y_3)}{\partial x^2}\right|_{x=\frac{1}{3}} &= - \frac{14192}{3645}
\end{aligned}
$$
:::

## 演習問題

::::exercise{level=basic}
$K=1$、$[0,1]$、$\lambda=1/2$、$f=x$ の Fredholm 方程式を解け。
:::solution
$c=\int u$ に対して $c=1/2+c/2$ なので $c=1$ である。$u=x+c/2=x+1/2$ となる。
:::
::::

::::exercise{level=basic}
Volterra の $Tu(x)=\int_0^x u(t)\dd t$ で $T^31$ を求めよ。
:::solution
順に $x,x^2/2,x^3/6$ となる。
:::
::::

::::exercise{level=basic}
$K(x,t)=xt$ の唯一の非零固有値を求めよ。
:::solution
像は $x$ の一次元空間である。$Tx=x\int_0^1t^2\dd t=x/3$ なので固有値は $1/3$ である。
:::
::::

::::exercise{level=basic}
$I^{1/2}1$ を積分せよ。
:::solution
$\int_0^x (x-t)^{-1/2}\dd t=2\sqrt x$ であり、$\Gamma(1/2)=\sqrt\pi$ で割って $2\sqrt{x/\pi}$ になる。
:::
::::

::::exercise{level=advanced}
$u-\int_0^1u(t)\dd t=x-1/2$ の全解を求めよ。
:::solution
右辺の積分は零である。$u=x-1/2+c$ を代入すると $\int u=c$ なので任意の $c$ で成り立つ。元の式を移項すると解は必ずこの形になる。
:::
::::

::::exercise{level=advanced}
$\|T\|\le1/4$、$\|f\|=1$、$\lambda=2$ で二次までの Neumann 和の誤差を評価せよ。
:::solution
$\|\lambda T\|\le1/2$ なので、三次以降の和は $(1/2)^3/(1-1/2)=1/4$ 以下である。
:::
::::

::::exercise{level=advanced}
$I^{1/2}u=x$ を解け。
:::solution
$I^{1/2}x=x^{3/2}\Gamma(2)/\Gamma(5/2)=4x^{3/2}/(3\sqrt\pi)$ である。これは $t=xv$ と Beta 積分から得る。微分して $u=2\sqrt{x/\pi}$ となる。合成則により直接の像も $I^11=x$ である。
:::
::::

::::exercise{level=advanced}
$y''=-y,\ y(0)=1,y'(0)=0$ を Volterra 反復で解け。
:::solution
$y=1-Ty$、$Tu=\int_0^x(x-t)u(t)\dd t$ である。$Tx^m=x^{m+2}/[(m+1)(m+2)]$ より $T^n1=x^{2n}/(2n)!$ である。交代級数の和は $\cos x$ となる。
:::
::::

## 確認テスト

::::quiz{answer="3" tol="1e-9"}
$K=xt$ の第 2 種 Fredholm 方程式が特異になる $\lambda$ を求めよ。

:::explanation
唯一の非零固有値は $\int_0^1t^2\dd t=1/3$ であり、$1-\lambda/3=0$ から $\lambda=3$ である。
:::
::::

::::quiz
Volterra 核の反復に階乗が現れる理由を選べ。

- [x] 順序付き積分領域の体積が立方体の $1/n!$ である
- [ ] 核はいつも定数である
- [ ] 反復で区間が一つ増えるたび同じ幅になる

:::explanation
立方体の点は測度零の等値面を除いて $n!$ 通りの並び順に分かれ、どれも座標交換で同じ体積を持つ。
:::
::::

::::quiz{answer="0.25" tol="1e-9"}
$K(x,t)=1/2$ を $[0,1]^2$ に置いたときの Hilbert–Schmidt ノルムの二乗を答えよ。

:::explanation
核の絶対値の二乗は $1/4$、積分領域の面積は $1$ なので $1/4$ である。
:::
::::

::::quiz
$Au=f$ の Fredholm 可解条件で直交すべき空間はどれか。

- [ ] 常に $\Ker A$
- [x] $\Ker A^*$
- [ ] $\operatorname{ran}A$

:::explanation
$\inner{Au}{v}=\inner{u}{A^*v}$ なので $v\in\Ker A^*$ なら零になる。値域が閉じることにより、この直交条件は十分でもある。
:::
::::

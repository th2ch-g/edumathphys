---
title: 楕円型方程式の弱解
description: Sobolev 空間における変分形式から解の存在、エネルギー評価、固有値の変分表示を導く。
notation:
  - symbol: 'V=H_0^1(\Omega)'
    meaning: 零境界値を持つ試験空間
  - symbol: 'V^*=H^{-1}(\Omega)'
    meaning: V の連続双対空間
  - symbol: 'a(u,v)'
    meaning: 作用素に対応する双線形形式
---

## この章の目標

- 微分方程式と境界条件を一階の積分式に直す。
- Lax–Milgram の定理で弱解を作り、データに対する評価を求める。
- コンパクト性から Fredholm の分岐と固有系を導く。
- Rayleigh 商で固有値を評価し、Weyl の主要項を求める。

## 2 階楕円型作用素と弱解

$\Omega\subset\R^n$ は有界開集合とする。$H^1(\Omega)$ は関数と一階弱微分が $L^2$ に入る空間である。$H^1_0(\Omega)$ は $C_c^\infty(\Omega)$ の $H^1$ 閉包である。境界が Lipschitz なら、後者はトレースが $0$ の関数に一致する。弱微分とトレースは [Sobolev 空間](/math/functional-analysis/sobolev-spaces/) を参照する。以下では実数値の関数を扱い、$H^1_0(\Omega)$ の連続双対空間を $H^{-1}(\Omega)$ と書く。

:::definition[一様楕円性と弱形式]{#def-pde-elliptic-weak}
実係数 $A=(a_{ij}),b,c$ が有界可測であり、ある $\lambda>0$ に対し $\xi^\top A(x)\xi\ge\lambda|\xi|^2$ がほとんどすべての $x$ とすべての $\xi\in\R^n$ について成り立つとする。この条件を :term[一様楕円性]{en="uniform ellipticity" yomi="いちようだえんせい"} という。$V=H^1_0(\Omega)$、$f\in V^*$ とし、$L=-\divg(A D\cdot)+b\cdot D+c$ の斉次 Dirichlet 弱解とは、

$$
a(u,v)=\int_\Omega(ADu\cdot Dv+b\cdot Du\,v+cuv)\,\dd x=\langle f,v\rangle
\quad(v\in V)
\label{eq-pde-elliptic-weak}
$$

を満たす $u\in V$ のことである。
:::

古典解には部分積分で弱形式が従う。弱解には二階の古典微分を要求せず、可測係数に対しても定義できる。反対に $A,u$ が十分滑らかなら $v\in C_c^\infty$ から分布の式が得られ、$Lu$ と $f$ が連続なら点ごとの式に戻る。ここで $V^*$ の元は $L^2$ 関数に限らない。例えば $F\in L^2(\Omega;\R^n)$ なら $\langle f,v\rangle=-\int F\cdot Dv$ は $f=\divg F$ というデータである。

:::example[媒質の界面]
$(-1,1)$ で $a=1$（左）、$a=2$（右）、境界値 $u(-1)=0$、$u(1)=3$ とする。$u(x)=2(x+1)$（左）、$u(x)=2+x$（右）とすれば、導関数は跳ぶが $au'=2$ は一定である。従って $\int au'v'=2\int v'=0$ となり、$-(au')'=0$ の弱解である。
:::

:::example[点源の空間]
一次元では $v(0)=\int_{-1}^0 v'$ より $|v(0)|\le\norm{v'}_2$ である。従って $f=\delta_0$ は $H^{-1}(-1,1)$ に属する。$u=(1-|x|)/2$ の傾きは左 $1/2$、右 $-1/2$ なので $\int u'v'=v(0)$ となる。
:::

### 発散形の界面条件とデータの双対表示

分割した二領域で $A$ が滑らかで、共通の滑らかな界面を持つとする。区分的な $C^2$ 解 $u$ が $H^1$ に属するには界面における値の跳躍が零でなければならない。実際、値の跳躍があれば分布勾配には面測度が現れ、$L^2$ の勾配にならない。共通の法線 $\nu$ を第一領域から第二領域へ向ける。双方で部分積分すると、試験関数の係数として $(A_1Du_1-A_2Du_2)\cdot\nu$ が界面に残る。右辺に面測度がなければこれも零である。従って発散形で連続するのは法線勾配そのものではなく法線流束である。

:::example[二層区間の流束と傾き]
$(0,2)$ で $a=1$（$x<1$）、$a=3$（$x>1$）とし、$-(au')'=0$、$u(0)=0,u(2)=4$ を課す。流束を $q=au'$ と置くと $u(2)-u(0)=q(1+1/3)=4$ より $q=3$ である。従って $u=3x$（$x\le1$）、$u=x+2$（$x\ge1$）である。界面の値は $3$ に一致するが、勾配は $3$ から $1$ に下がる。$\int au'\phi'=3[\phi]_0^2=0$ なので全ての零トレース試験関数で弱式を満たす。
:::

$V=H_0^1(\Omega)$ に $\norm{Dv}_2$ を入れる。有界開集合では :ref[thm-fa-poincare] により、このノルムは $H^1$ のノルムと同値である。任意の $f\in V^*$ はある $q\in L^2(\Omega;\R^n)$ で $\langle f,v\rangle=\int q\cdot Dv$ と書ける。実際、Riesz の表現で $\langle f,v\rangle=\int Dw\cdot Dv$ となる $w\in V$ が存在するので $q=Dw$ とすればよい。従って $H^{-1}$ の外力は $-\divg q$ として実現できる。

:::example[非定数係数の流束と弱式]
$(0,2)$ で $A=1+x$、$f=4x$、同次 Dirichlet 条件とする。流束を積分して値の差が零になるように定数を決める。

まず流束 $q=Au_x$ を未知量にする。$q_x=-f$ を積分すると、任意の定数 $q(0)$ について $q(x)=q(0)-\int_0^x f$ である。さらに $q/A$ を積分し、指定された両端の値または平均に合わせて積分定数を決める。以下では計算後の $u$ について関係を逆に確かめる。

$$
\begin{aligned}
u(x)&=x \left(2 - x\right)\\
u_x&=\frac{d}{d x} x \left(2 - x\right)=- 2 \left(x - 1\right)\\
q&=\left(2 - 2 x\right) \left(x + 1\right)=- 2 \left(x - 1\right) \left(x + 1\right)\\
-f&=\frac{d}{d x} \left(- 2 \left(x - 1\right) \left(x + 1\right)\right)=- 4 x
\end{aligned}
$$

次に値の条件と流出量を両端で区別して評価する。区間の外法線は左端で $-1$、右端で $1$ である。このため源項の積分は左側の流束から右側の流束を引いたものになる。流束を係数で割った積分は両端の値の差に等しい。

$$
\begin{aligned}
u(0)&=0,\qquad x=0\\
u(L)&=0,\qquad x=2\\
q(0)&=2,\qquad x=0\\
q(L)&=-6,\qquad x=2\\
\int_0^L f\,\mathrm{d}x&=\int\limits_{0}^{2} 4 x\, dx=8\\
\int_0^L (q/A)\,\mathrm{d}x&=\int\limits_{0}^{2} \left(2 - 2 x\right)\, dx=0
\end{aligned}
$$

弱形式も具体的に比べる。試験関数 $\phi=x(L-x)$ は両端で零であり、積分の間に境界項は出ない。以下の値が一致するのは、任意の零トレース試験関数に対しても部分積分できることの一例である。

$$
\begin{aligned}
\phi&=x \left(2 - x\right)\\
\phi^{\prime}&=\frac{d}{d x} x \left(2 - x\right)=- 2 \left(x - 1\right)\\
\int_0^L q\phi^{\prime}\,\mathrm{d}x&=\int\limits_{0}^{2} \left(- 2 \left(2 - 2 x\right) \left(x - 1\right) \left(x + 1\right)\right)\, dx=\frac{16}{3}\\
\int_0^L f\phi\,\mathrm{d}x&=\int\limits_{0}^{2} 4 x^{2} \left(2 - x\right)\, dx=\frac{16}{3}
\end{aligned}
$$

最後に内部での値、勾配、流束を同じ点で測る。区間全体での勾配の寄与には $A$ の重みが付くので、点での値だけでは全体のエネルギーは分からない。

$$
\begin{aligned}
u&=\frac{3}{4},\qquad x=\frac{3}{2}\\
u_x&=-1,\qquad x=\frac{3}{2}\\
q&=- \frac{5}{2},\qquad x=\frac{3}{2}\\
\frac12\int_0^L A(u_x)^2\,\mathrm{d}x&=\int\limits_{0}^{2} \frac{\left(2 - 2 x\right)^{2} \left(x + 1\right)}{2}\, dx=\frac{8}{3}
\end{aligned}
$$

流束は $q=(1+x)(2-2x)=2-2x^2$ となり、定数項は微分で消える。$A$ を微分する項を落とすと、$f=4x$ と $-Au_{xx}=2+2x$ を混同する。
:::

## Lax–Milgram の定理による存在

$V$ には $\norm u_V=\norm{Du}_2$ を使う。$\Omega$ を幅 $D$ の帯に含め、$v\in C_c^\infty$ を零延長すると、帯の端からの積分で $|v(x)|^2\le D\int|\partial_1v|^2$ となる。幅方向に積分して $\norm v_2\le D\norm{Dv}_2$ を得る。密度により $V$ でも成立するので $\norm{\cdot}_V$ は $H^1$ ノルムと同値である。以下ではこの Poincaré 定数を $C_P$ と書く。

:::theorem[強圧形式の弱解]{#thm-pde-elliptic-coercive}
$a$ が有界で $a(v,v)\ge\alpha\norm v_V^2$（$\alpha>0$）なら、各 $f\in V^*$ に一意な弱解があり、$\norm u_V\le\alpha^{-1}\norm f_{V^*}$ である。
:::

:::proof
[Hilbert 空間の Riesz 表現定理](/math/functional-analysis/hilbert-spaces/) により $(Tu,v)_V=a(u,v)$ となる有界作用素 $T$ を定める。強圧性から $\norm{Tu}_V\ge\alpha\norm u_V$ であり、$T$ は単射で像は閉である。実際 $Tu_j$ が収束すれば $u_j$ は Cauchy であり、$T$ の連続性で極限も像に入る。像に直交する $v$ があれば全ての $u$ で $a(u,v)=0$ となり、$u=v$ として $v=0$ である。従って像は稠密かつ閉なので全空間である。$f$ の Riesz 表現を $Tu$ として解けば存在と一意性が従う。最後に $v=u$ を弱形式に入れ $\alpha\norm u_V^2\le\norm f_{V^*}\norm u_V$ を得る。
:::

例えば $b=0$、$c\ge0$ なら $\alpha=\lambda$ と選べる。一般の有界な $b,c$ では $B=\norm b_\infty$、$C_-=\norm{c_-}_\infty$ とおくと、$B\norm{Du}_2\norm u_2\le\lambda\norm{Du}_2^2/2+B^2\norm u_2^2/(2\lambda)$ なので

$$
a(u,u)\ge\tfrac{\lambda}{2}\norm{Du}_2^2-\gamma_0\norm u_2^2,\qquad
\gamma_0=B^2/(2\lambda)+C_- .
\label{eq-pde-garding}
$$

この形を :term[Gårding の不等式]{en="Gårding inequality" yomi="ごーでぃんぐのふとうしき"} という。$L+\gamma I$ で $\gamma>\gamma_0$ とすれば強圧性を得る。非斉次境界値が $G\in H^1$ のトレースで与えられる場合、未知量を $u=G+w$ と置き、右辺を $f-a(G,\cdot)$ に替えて $w\in V$ を求める。延長 $G$ を変えてもその差は $V$ に属するため、最終的な $u$ は同じである。

:::example[零階項で強圧性を失う場合]
$\Omega=(0,\pi)$、$L=-\partial_{xx}-1$ とすると $a(\sin x,\sin x)=0$ である。従って $\lambda=1$ の主部が正でも、全作用素に強圧性はない。$L+2$ なら形式は $\int(|u'|^2+u^2)$ となり強圧性が回復する。
:::

:::example[低次微分項があっても解ける範囲]
$L=-u''+b u'$ を $(0,1)$ の零端点で考え、$b$ は実定数とする。$\int b u'u=b[u^2]_0^1/2=0$ なので、$b$ の大きさにかかわらず $a(u,u)=\norm{u'}_2^2$ である。:eqref[eq-pde-garding] の一般評価で負の項が出るのは絶対値で評価したためであり、この係数では部分積分でちょうど消える。
:::

### 非対称な作用素での強圧性の役割

強圧形式 $a$ を $V$ の内積に関して $(Tu,v)_V=a(u,v)$ と表す。$|(Tu,u)_V|\le\norm{Tu}_V\norm u_V$ より $\norm{Tu}_V\ge\alpha\norm u_V$ である。従って $Tu_n$ が収束すれば $u_n$ は Cauchy で、像は閉じている。さらに $w$ が像に直交すれば $a(u,w)=0$ が全 $u$ について成り立ち、$u=w$ を入れて $w=0$ を得る。従って像は稠密でもあり、全空間である。対称性を使わずに解が存在するのは、強圧性が作用素とその随伴の両方を制御しているためである。

:::example[定数ドリフトはエネルギーの対角項を変えない]
$(0,1)$ で $a(u,v)=\int(u'v'+6u'v)$ とする。$u\in H_0^1$ で $\int u'u=0$ だから $a(u,u)=\int|u'|^2$。他方 $u=x(1-x)$ なら $Lu=-u''+6u'=8-12x$ であり、ドリフトは実際の解を変える。例えば $f=8-12x$ に対しこの $u$ が一意解で、$\int|u'|^2=1/3$。対角エネルギーに現れないことと、方程式に影響しないことは別である。
:::

:::example[近似空間での強圧定数]
正規直交な二基底上の剛性行列を考え、強圧性の定数を求める。

剛性行列を $K$ とする。対角成分は座標方向での値を与える。固有方程式から直交する方向を取り出すため、特性多項式 $\det(\lambda I-K)$ を計算する。

$$
\begin{aligned}
K&=\left[\begin{matrix}\frac{178}{25} & - \frac{96}{25}\\- \frac{96}{25} & \frac{122}{25}\end{matrix}\right]\\
\det(\lambda I-K)&=\lambda^{2} - 12 \lambda + 20=\left(\lambda - 10\right) \left(\lambda - 2\right)\\
d_1&=2,\quad v_1=\left[\begin{matrix}\frac{3}{5}\\\frac{4}{5}\end{matrix}\right]\\
d_2&=10,\quad v_2=\left[\begin{matrix}- \frac{4}{5}\\\frac{3}{5}\end{matrix}\right]
\end{aligned}
$$

$V$ の列に固有ベクトルを並べる。$V^\top V=I$ と $V^\top KV=D$ を使えば、元の座標 $X$ と固有座標 $Y=V^\top X$ の双方で基底係数を記述できる。二次形式は固有座標にすると交差項を持たない。

$$
\begin{aligned}
V&=\left[\begin{matrix}\frac{3}{5} & - \frac{4}{5}\\\frac{4}{5} & \frac{3}{5}\end{matrix}\right]\\
V^\top KV&=\left[\begin{matrix}2 & 0\\0 & 10\end{matrix}\right]\\
X^\top KX&=\frac{2 \left(89 x_{1}^{2} - 96 x_{1} x_{2} + 61 x_{2}^{2}\right)}{25}\\
\sum_j d_jY_j^2&=2 Y_{1}^{2} + 10 Y_{2}^{2}
\end{aligned}
$$

一つの基底係数で固有係数を計算し、エネルギーの二次形式を求める。直交変換なので二乗ノルムは変わらない。二次形式を二乗ノルムで割ると、固有値を非負の重みで平均した値になる。

$$
\begin{aligned}
X&=\left[\begin{matrix}1\\2\end{matrix}\right]\\
Y&=\left[\begin{matrix}\frac{11}{5}\\\frac{2}{5}\end{matrix}\right]\\
|X|^2&=5\\
X^\top KX&=\frac{282}{25}\\
\frac{X^\top KX}{|X|^2}&=\frac{282}{125}\\
\min_j d_j&=2\\
\max_j d_j&=10
\end{aligned}
$$

下限は固有ベクトル $v_1$ で達し、$|X|=1$ でエネルギーは $2$ 以上である。逆行列のノルムは $1/2$。強圧定数を大きく見積もると、小さい固有方向の誤差を過小評価する。
:::

### 非斉次データと自然境界条件を同じ形式で扱う

境界を、値を固定する部分 $\Gamma_D$ と法線流束を与える部分 $\Gamma_N$ に分ける。Lipschitz 領域 $\Omega$ を連結かつ有界とし、$\Gamma_D$ は正の表面測度を持つとする。$V_D=\{v\in H^1(\Omega):\operatorname{Tr}v=0\text{ on }\Gamma_D\}$ は閉部分空間である。この空間では $\norm v_2\le C\norm{Dv}_2$ が成り立つ。実際、反例列を $\norm{v_j}_2=1$、$\norm{Dv_j}_2\to0$ に規格化し、Rellich の定理から $L^2$ 強極限を取る。極限は勾配が零なので定数である。トレース条件は $H^1$ の弱収束で保たれ、定数が正の測度を持つ境界で零なら定数全体が零である。しかし $L^2$ ノルムは $1$ なので矛盾する。

$-\divg(A Du)=f$ に $\partial_\nu^A u=(A Du)\cdot\nu=g$ を $\Gamma_N$ に与え、値の延長 $G\in H^1$ を取ると、$v=u-G$ の方程式は
$$
\int_\Omega A Dv\cdot D\phi=\langle f,\phi\rangle
+\int_{\Gamma_N}g\,\operatorname{Tr}\phi-\int_\Omega A DG\cdot D\phi
\qquad(\phi\in V_D)
$$
となる。例えば $f\in L^2(\Omega)$ と $g\in L^2(\Gamma_N)$ ならトレース評価により右辺は有界である。左辺の強圧性は楕円性と上の Poincaré 不等式によるので、Lax–Milgram から解が得られる。全境界が Neumann の場合はこの証明で定数方向を除けない。そのとき $\phi=1$ から適合条件が出て、平均零の部分空間で解くことになる。

:::example[区間で Dirichlet と Neumann を一つずつ使う]
$(0,2)$ で $-u''=3$、$u(0)=1$、$u'(2)=2$ とする。$u=-3x^2/2+Ax+1$ とおくと $-6+A=2$ なので $A=8$、従って $u(2)=11$ である。$\phi=x$ を試験関数に取ると、左辺は $\int_0^2u'=10$、右辺は $\int_0^2 3x+2\phi(2)=6+4=10$ である。$\phi(2)$ は零に固定されず、この境界項が自然条件の情報を伝える。
:::

### 純 Neumann 問題の適合条件

連結な有界 Lipschitz 領域 $\Omega$ で、$A$ は有界かつ一様楕円とする。全境界で $(ADu)\cdot\nu=g$ を指定する弱問題を

$$
\int_\Omega ADu\cdot Dv=\ell(v)
\quad(v\in H^1(\Omega)),\qquad
\ell(v)=\int_\Omega fv+\int_{\partial\Omega}g\,\operatorname{Tr}v
$$

で表す。ここでは $f\in L^2(\Omega)$、$g\in L^2(\partial\Omega)$ とする。$v=1$ を入れると $\int_\Omega f+\int_{\partial\Omega}g=0$ が必要である。この条件は十分でもある。平均零の閉部分空間では :ref[thm-fa-poincare] により勾配ノルムが $H^1$ ノルムと同値になり、トレース評価で右辺も有界になる。そこで :ref[thm-pde-elliptic-coercive] と同じ Hilbert 空間の議論で一意な平均零解を得る。任意の $v$ を平均零の部分と定数に分けると、適合条件により定数部分に対しても弱式が成り立つ。他の解との差を試験関数にすると勾配が零になるため、全解は平均零解に定数を加えたものである。

:::example[両端の外向き流束を指定する]
$(0,1)$ で $-u''=2$、外向き流束を $g(0)=-1,\ g(1)=-1$ とする。適合条件は $2-1-1=0$ である。$-u'(0)=-1$ と $u'(1)=-1$ から $u=x(1-x)+C$ を得る。平均零の解では $\int_0^1x(1-x)\,\dd x=1/6$ より $C=-1/6$ である。両端で外向きの法線が反対を向くので、二つの流束が等しくても $u'$ の値は異なる。
:::

## エネルギー評価

:::proposition[データと係数への依存性]{#prop-pde-elliptic-stability}
強圧定数が共通に $\alpha>0$ である二つの形式 $a_1,a_2$ と解 $u_1,u_2$ に対して

$$
\norm{u_1-u_2}_V
\le\alpha^{-1}\left(\norm{f_1-f_2}_{V^*}
+\norm{a_1-a_2}\norm{u_2}_V\right)
\label{eq-pde-elliptic-stability}
$$

である。形式ノルムは $\norm v_V,\norm w_V$ が $1$ 以下の組での絶対値の上限である。
:::

:::proof
$w=u_1-u_2$ とし、$a_1(w,v)=\langle f_1-f_2,v\rangle+(a_2-a_1)(u_2,v)$ を得る。$v=w$ を選び強圧性と右辺の双対評価を適用する。$w=0$ なら結論は自明であり、それ以外では $\norm w_V$ で割る。
:::

$f\in L^2$ なら $\norm f_{V^*}\le C_P\norm f_2$ である。$-\divg(A Du)=\divg F$ なら $\lambda\norm{Du}_2^2\le\norm F_2\norm{Du}_2$ なので $C_P$ は必要ない。右辺の書き方で必要なノルムが変わる。

:::example[エネルギーを直接計算する]
$(0,1)$ で $-u''=2$、両端 $0$ なら $u=x(1-x)$ であり、$\norm{u'}_2^2=\int_0^1(1-2x)^2=1/3$ である。$\int2u=1/3$ なので弱式で $v=u$ を選んだ値と一致する。
:::

:::example[係数摂動に対する解の誤差]
$(0,1)$ で $L_1=-\partial_{xx}$、$L_2=-2\partial_{xx}$、$f=2$ とする。対応する解は $u_1=x(1-x)$、$u_2=x(1-x)/2$ である。よって $\norm{u_1-u_2}_V=1/(2\sqrt3)$ である。式 :eqref[eq-pde-elliptic-stability] は $\alpha=1$、$\norm{a_1-a_2}=1$ で $\norm{u_1-u_2}_V\le\norm{u_2}_V=1/(2\sqrt3)$ を与える。この例では上界が達成される。
:::

### 残差から得る両側の評価

任意の近似 $v\in V$ に対し $r=f-Lv\in V^*$ を定める。厳密解との差 $e=u-v$ は $Le=r$ である。$a$ が有界定数 $M$、強圧定数 $\alpha$ を持てば
$$
\alpha\norm e_V\le \norm r_{V^*}\le M\norm e_V .
$$
左は $\alpha\norm e_V^2\le a(e,e)=\langle r,e\rangle$、右は $r$ の定義から出る。例えば $(0,1)$ で $L=-\partial_{xx}$、$f=2$、$v=\tfrac12x(1-x)$ なら $r=1$ であり、$u-v=\tfrac12x(1-x)$ となる。この作用素は Riesz 表現と一致するため $\norm r_{V^*}=\norm e_V=1/(2\sqrt3)$ である。残差の空間を $L^2$ に置くと、近似の微分が存在しない場合に計算できなくなる。

:::example[境界誤差と内部残差を分ける]
$(0,1)$ で $A=4$、$f=8$、$u(0)=0,u(1)=1$ とする。境界のアフィン延長と原方程式の解を比べる。

まず流束 $q=Au_x$ を未知量にする。$q_x=-f$ を積分すると、任意の定数 $q(0)$ について $q(x)=q(0)-\int_0^x f$ である。さらに $q/A$ を積分し、指定された両端の値または平均に合わせて積分定数を決める。以下では計算後の $u$ について関係を逆に確かめる。

$$
\begin{aligned}
u(x)&=- x^{2} + 2 x\\
u_x&=\frac{d}{d x} \left(- x^{2} + 2 x\right)=- 2 \left(x - 1\right)\\
q&=8 - 8 x=- 8 \left(x - 1\right)\\
-f&=\frac{d}{d x} \left(- 8 \left(x - 1\right)\right)=-8
\end{aligned}
$$

次に値の条件と流出量を両端で区別して評価する。区間の外法線は左端で $-1$、右端で $1$ である。このため源項の積分は左側の流束から右側の流束を引いたものになる。流束を係数で割った積分は両端の値の差に等しい。

$$
\begin{aligned}
u(0)&=0,\qquad x=0\\
u(L)&=1,\qquad x=1\\
q(0)&=8,\qquad x=0\\
q(L)&=0,\qquad x=1\\
\int_0^L f\,\mathrm{d}x&=\int\limits_{0}^{1} 8\, dx=8\\
\int_0^L (q/A)\,\mathrm{d}x&=\int\limits_{0}^{1} \left(2 - 2 x\right)\, dx=1
\end{aligned}
$$

弱形式も具体的に比べる。試験関数 $\phi=x(L-x)$ は両端で零であり、積分の間に境界項は出ない。以下の値が一致するのは、任意の零トレース試験関数に対しても部分積分できることの一例である。

$$
\begin{aligned}
\phi&=x \left(1 - x\right)\\
\phi^{\prime}&=\frac{d}{d x} x \left(1 - x\right)=1 - 2 x\\
\int_0^L q\phi^{\prime}\,\mathrm{d}x&=\int\limits_{0}^{1} \left(- 8 \left(1 - 2 x\right) \left(x - 1\right)\right)\, dx=\frac{4}{3}\\
\int_0^L f\phi\,\mathrm{d}x&=\int\limits_{0}^{1} 8 x \left(1 - x\right)\, dx=\frac{4}{3}
\end{aligned}
$$

最後に内部での値、勾配、流束を同じ点で測る。区間全体での勾配の寄与には $A$ の重みが付くので、点での値だけでは全体のエネルギーは分からない。

$$
\begin{aligned}
u&=\frac{3}{4},\qquad x=\frac{1}{2}\\
u_x&=1,\qquad x=\frac{1}{2}\\
q&=4,\qquad x=\frac{1}{2}\\
\frac12\int_0^L A(u_x)^2\,\mathrm{d}x&=\int\limits_{0}^{1} 2 \left(2 - 2 x\right)^{2}\, dx=\frac{8}{3}
\end{aligned}
$$

近似を $v=x$ とすると端点は正確だが $Lv=0$、従って残差は定数 $8$ である。誤差 $e=x-x^2$ の勾配ノルムは $1/\sqrt3$、残差の双対ノルムは $4/\sqrt3$ であり、強圧定数 $4$ の下界が等号になる。近似の境界値も異なる場合は、その誤差を先に延長で引いてから $H_0^1$ の評価を使う。
:::

:::example[残差が直交しても誤差が零とは限らない]
$(0,1)$ で $A=1+x$、$c=0$、厳密解 $u=x^2(1-x)^2$ となる外力を置く。試験空間は二次多項式と三次多項式の二次元空間である。

剛性行列を作る前に、方程式と試験空間を決める。ここでは全ての基底が両端で消える。厳密解を含まない試験空間に制限すると、近似は割り当てた方向の弱形式を満たすように決まる。厳密解は誤差を測るときに使う。

$$
\begin{aligned}
A&=x + 1\\
u&=x^{2} \left(1 - x\right)^{2}\\
f&=- 2 \left(8 x^{3} - 3 x^{2} - 4 x + 1\right)\\
\psi_1&=x \left(1 - x\right)\\
\psi_2&=x \left(1 - x\right) \left(2 x - 1\right)
\end{aligned}
$$

剛性の成分は $K_{ij}=\int A\psi^{\prime}_j\psi^{\prime}_i+\int c\psi_j\psi_i$、右辺は $F_i=\int f\psi_i$ である。基底を正規直交には取っていないので、内積の質量行列 $M_{ij}=\int\psi_j\psi_i$ も別に求める。

$$
\begin{aligned}
K&=\left[\begin{matrix}\frac{1}{2} & \frac{1}{15}\\\frac{1}{15} & \frac{3}{10}\end{matrix}\right]\\
M&=\left[\begin{matrix}\frac{1}{30} & 0\\0 & \frac{1}{210}\end{matrix}\right]\\
F&=\left[\begin{matrix}\frac{1}{10}\\\frac{1}{210}\end{matrix}\right]\\
\det(K-\lambda M)&=\frac{\lambda^{2}}{6300} - \frac{13 \lambda}{1050} + \frac{131}{900}=\frac{\lambda^{2} - 78 \lambda + 917}{6300}
\end{aligned}
$$

連立方程式 $Kc=F$ を解き、元の関数に戻す。残差 $r=f-Lu_N$ を計算すると一般に零ではないが、選んだ試験関数に対する積分は全て零になる。PDE を点で満たすことと有限次元の弱形式を満たすことを区別する。

$$
\begin{aligned}
c&=\left[\begin{matrix}\frac{187}{917}\\- \frac{27}{917}\end{matrix}\right]\\
u_N&=\frac{2 x \left(x - 1\right) \left(27 x - 107\right)}{917}\\
e&=u-u_N=\frac{x \left(x - 1\right) \left(917 x^{2} - 971 x + 214\right)}{917}\\
r&=- \frac{4 \left(3668 x^{3} - 1497 x^{2} - 1647 x + 539\right)}{917}\\
(\int r\psi_i)_i&=\left[\begin{matrix}0\\0\end{matrix}\right]
\end{aligned}
$$

最後に誤差の二乗ノルムを二通りに計算する。近似空間へのエネルギー直交性により、誤差の勾配エネルギーは全体のエネルギーから近似の分を引いた値になる。ここでは二次形式そのものを用いており、エネルギー汎関数の二分の一は付けていない。

$$
\begin{aligned}
a(u,u)&=\frac{1}{35}\\
a(u_N,u_N)&=\frac{130}{6419}\\
a(e,e)&=\frac{267}{32095}\\
\norm e_2^2&=\int\limits_{0}^{1} \frac{x^{2} \left(x - 1\right)^{2} \left(917 x^{2} - 971 x + 214\right)^{2}}{840889}\, dx=\frac{17057}{264880035}
\end{aligned}
$$

$u$ の四次項は試験空間では作れない。従って二つの試験方向で残差が消えても、どの試験関数でも消えることにはならない。試験空間を増やすと最小となるエネルギー誤差は増えない。
:::

### Galerkin 直交性からの準最良近似

有界・強圧な形式 $a$ の解 $u$ と有限次元部分空間 $V_h$ での解 $u_h$ を考える。部分空間でも同じ強圧性が使えるので、有限次元方程式は一意に可解である。差 $e=u-u_h$ に対し $a(e,v_h)=0$（$v_h\in V_h$）。従って任意の $w_h\in V_h$ について
$$
\alpha\norm e_V^2\le a(e,e)=a(e,u-w_h)
\le M\norm e_V\norm{u-w_h}_V .
$$
$e=0$ なら結論は明らかであり、それ以外は割って $\norm e_V\le (M/\alpha)\inf_{w_h}\norm{u-w_h}_V$ を得る。これを :term[Céa の評価]{en="Céa estimate" yomi="せあのひょうか"} という。対称形式なら $a(e,e)^{1/2}$ が内積ノルムであり、$a(u-w_h,u-w_h)=a(e,e)+a(u_h-w_h,u_h-w_h)$ と分解できるので、エネルギーノルムでの最良近似になる。非対称な場合に同じ等式を使おうとすると $a(u_h-w_h,e)$ は零とは限らないため、係数 $M/\alpha$ が残る。

## 弱解に対する Fredholm の択一定理

:::theorem[楕円型の Fredholm 条件]{#thm-pde-elliptic-fredholm}
:ref[def-pde-elliptic-weak] の条件下で、$L:V\to V^*$ の核と随伴形式の核は有限次元で次元が等しい。$Lu=f$ が可解である必要十分条件は、$a(v,w)=0$ が全ての $v\in V$ で成立する $w$ に対して $\langle f,w\rangle=0$ である。
:::

:::proof
$J:V\to V^*$ を $\langle Ju,v\rangle=(u,v)_{L^2}$ で定める。有界領域で $V$ の元を零延長すると、$V\hookrightarrow L^2$ は Rellich のコンパクト性によりコンパクトである。よって $J$ もコンパクトである。:eqref[eq-pde-garding] により $B=L+\gamma J$ は可逆なので、$L=B(I-\gamma B^{-1}J)$ と因数分解できる。$I$ からコンパクト作用素を引いたものには [Fredholm の択一定理](/math/functional-analysis/compact-operators/) が適用できる。ゆえに $L$ は閉像、有限次元の核、指数 $0$ を持つ。双対の直交関係で像の消去空間は $a(v,w)=0$ を満たす $w$ の集合である。閉像のため、その消去空間を消す $f$ は像に入り、条件は十分でもある。
:::

### コンパクト性を零延長から見る

$H^1_0(\Omega)$ の有界列を全空間へ零延長すると、同じ有界な箱に台を持つ $H^1$ の有界列になる。延長は $C_c^\infty(\Omega)$ の近似から得られ、勾配ノルムを保つ。線分積分と Cauchy–Schwarz により $\norm{u(\cdot+h)-u}_2\le |h|\norm{Du}_2$ である。滑らかな場合に証明して $H^1$ の近似に移せる。格子で箱を一辺 $\ell$ の立方体に分け、各立方体の平均で置き換えると、Jensen の不等式により

$$
\sum_Q\int_Q|u-u_Q|^2\le \sum_Q |Q|^{-1}\int_Q\int_Q|u(x)-u(y)|^2\,\dd y\,\dd x
\le C_n\ell^2\norm{Du}_2^2 .
$$

最後の不等式は $y=x+h$ として、$|h_i|\le\ell$ で平行移動の上界を積分すれば得られる。平均への置換は有限階であり、近似誤差は列に一様に小さくできる。従って対角部分列を取れば $L^2$ で Cauchy になる。この場合のコンパクト性には境界が滑らかであることは要らない。

:::example[共鳴と適合するデータ]
$L=-\partial_{xx}-1$ を $(0,\pi)$ の Dirichlet 条件で用いる。核は $\sin x$ の張る空間である。これは弱方程式から $u''=-u$ が従い、一次元では $u\in H^2$ で ODE を満たすためである。$f=\sin x$ なら $\int f\sin x=\pi/2\ne0$ なので不可解である。$f=\sin2x$ ならこの積分は $0$ で、$u=\sin2x/3+C\sin x$ が全ての解である。
:::

:::example[随伴の核は元の核と異なる]
$(0,\pi)$ の零端点で $L=-\partial_{xx}+2\partial_x-2$ とする。$u=e^x v$ を代入すると $Lu=e^x(-v''-v)$ なので、核は $e^x\sin x$ が張る。随伴は微分項の符号を逆にした作用素で、核は $e^{-x}\sin x$ が張る。従って $f=e^x\sin 2x$ は $\int f(e^{-x}\sin x)=0$ で可解であり、特解は $e^x\sin2x/3$ である。
:::

### 共鳴における核と商空間

$T=I-C$、$C$ はコンパクトとする。核の直交補空間で $\norm{Tu}\ge c\norm u$ が成り立つ。そうでなければ $\norm{u_n}=1$、$u_n\perp\ker T$、$Tu_n\to0$ を取れる。$Cu_n$ の収束部分列を取ると $u_n=Cu_n+Tu_n$ も強収束し、極限はノルム $1$ で核に入ると同時に核に直交する。これは矛盾である。核を足す自由度を取り除けば可解データから解への写像が有界であり、不安定さは共鳴そのものよりも共鳴へ近づくパラメータに現れる。

:::example[隣のモードだけを解く]
$(0,\pi)$ で $L=-\partial_{xx}-9$ とする。$\sin(3x)$ の方向が核であり、$f=5\sin(2x)$ はそれに直交する。$L\sin(2x)=-5\sin(2x)$ なので全解は $-\sin(2x)+c\sin(3x)$ である。$f$ に任意の非零の $\sin(3x)$ 成分を足すと、弱式を $\sin(3x)$ で試験したとき左辺が零なので解がない。この場合の直交条件は微分方程式を解く前に調べられる。
:::

:::example[離散の共鳴条件と直交射影]
三次元近似の剛性行列に零固有値が現れるとする。右辺として与えたベクトルの可解性を判定する。

剛性行列を $K$ とする。対角成分は座標方向での値を与える。固有方程式から直交する方向を取り出すため、特性多項式 $\det(\lambda I-K)$ を計算する。

$$
\begin{aligned}
K&=\left[\begin{matrix}\frac{48}{25} & - \frac{36}{25} & 0\\- \frac{36}{25} & \frac{27}{25} & 0\\0 & 0 & 8\end{matrix}\right]\\
\det(\lambda I-K)&=\lambda^{3} - 11 \lambda^{2} + 24 \lambda=\lambda \left(\lambda - 8\right) \left(\lambda - 3\right)\\
d_1&=0,\quad v_1=\left[\begin{matrix}\frac{3}{5}\\\frac{4}{5}\\0\end{matrix}\right]\\
d_2&=3,\quad v_2=\left[\begin{matrix}- \frac{4}{5}\\\frac{3}{5}\\0\end{matrix}\right]\\
d_3&=8,\quad v_3=\left[\begin{matrix}0\\0\\1\end{matrix}\right]
\end{aligned}
$$

$V$ の列に固有ベクトルを並べる。$V^\top V=I$ と $V^\top KV=D$ を使えば、元の座標 $X$ と固有座標 $Y=V^\top X$ の双方で基底係数を記述できる。二次形式は固有座標にすると交差項を持たない。

$$
\begin{aligned}
V&=\left[\begin{matrix}\frac{3}{5} & - \frac{4}{5} & 0\\\frac{4}{5} & \frac{3}{5} & 0\\0 & 0 & 1\end{matrix}\right]\\
V^\top KV&=\left[\begin{matrix}0 & 0 & 0\\0 & 3 & 0\\0 & 0 & 8\end{matrix}\right]\\
X^\top KX&=\frac{48 x_{1}^{2} - 72 x_{1} x_{2} + 27 x_{2}^{2} + 200 x_{3}^{2}}{25}\\
\sum_j d_jY_j^2&=3 Y_{2}^{2} + 8 Y_{3}^{2}
\end{aligned}
$$

一つの基底係数で固有係数を計算し、エネルギーの二次形式を求める。直交変換なので二乗ノルムは変わらない。二次形式を二乗ノルムで割ると、固有値を非負の重みで平均した値になる。

$$
\begin{aligned}
X&=\left[\begin{matrix}1\\2\\0\end{matrix}\right]\\
Y&=\left[\begin{matrix}\frac{11}{5}\\\frac{2}{5}\\0\end{matrix}\right]\\
|X|^2&=5\\
X^\top KX&=\frac{12}{25}\\
\frac{X^\top KX}{|X|^2}&=\frac{12}{125}\\
\min_j d_j&=0\\
\max_j d_j&=8
\end{aligned}
$$

零固有値に対して右辺が持つ係数 $v_1^\top X$ は非零である。したがってこの右辺は可解でない。$X$ から $(v_1^\top X)v_1$ を引いた右辺なら、残りの固有係数を $3,8$ で割れば解を求められる。
:::

### 対称性がない場合の適合条件を解き分ける

$(0,\pi)$ で $Lu=-u''+2u'-2u$ を零境界で考える。$u=e^xv$ を代入すると $Lu=e^x(-v''-v)$ である。実際 $u'=e^x(v+v')$、$u''=e^x(v+2v'+v'')$ なので一階項が相殺する。従って核は $e^x\sin x$ の張る空間である。一方 $L^*=-\partial_{xx}-2\partial_x-2$ であり、$w=e^{-x}z$ を入れると $L^*w=e^{-x}(-z''-z)$ になる。従って適合条件は $\int_0^\pi f e^{-x}\sin x=0$ である。

:::example[元の核に直交する条件と混同しない]
$f=e^x\sin 2x$ とすると、随伴核との内積は $\int\sin2x\sin x=0$ なので可解である。$-v''-v=\sin2x$ より一つの解は $u=e^x\sin2x/3$。ところが元の核との内積は $\int_0^\pi e^{2x}\sin2x\sin x$ であり、積を $(\cos x-\cos3x)/2$ と直して積分すると $-8(e^{2\pi}+1)/65\ne0$ になる。方程式の可解性を判定する相手は随伴核であり、元の核は解に自由に足せる方向を示す。
:::

## 固有値問題と固有関数

:::theorem[対称形式の固有関数展開]{#thm-pde-elliptic-spectrum}
$a$ が対称で :eqref[eq-pde-garding] を満たすとする。$L$ の固有値は有限重複度で下から有界な列 $\lambda_1\le\lambda_2\le\cdots\to+\infty$ に並び、対応する固有関数は $L^2(\Omega)$ の正規直交基底を作る。
:::

:::proof
$\gamma>\gamma_0$ とし、$f\in L^2$ を弱解 $u=(L+\gamma I)^{-1}f\in V$ へ写し、$L^2$ に埋め込んだ作用素を $K$ とする。強圧評価とコンパクト埋め込みから $K$ はコンパクトである。$a_\gamma=a+\gamma(\cdot,\cdot)_2$ が対称なので $(Kf,g)_2=a_\gamma(Kf,Kg)=(f,Kg)_2$ であり、$(Kf,f)_2=a_\gamma(Kf,Kf)\ge0$ である。また $Kf=0$ なら全 $v\in V$ について $(f,v)_2=0$ で、$V$ の稠密性から $f=0$ である。従って :ref[thm-fa-compact-spectral] が正の固有値 $\mu_j\downarrow0$ と正規直交基底 $e_j$ を与える。固有式は $e_j\in V$ と $a(e_j,v)=(\mu_j^{-1}-\gamma)(e_j,v)_2$ を意味する。よって $\lambda_j=\mu_j^{-1}-\gamma$ である。
:::

### レゾルベントと共鳴近くでの増幅

$\zeta$ が実数で、どの固有値とも一致しないとする。$f\in L^2$ なら $(L-\zeta I)u=f$ の係数は $u_j=f_j/(\lambda_j-\zeta)$ である。$\lambda_j\to\infty$ なので $\sup_j(\lambda_j+\gamma)/|\lambda_j-\zeta|^2<\infty$ であり、この係数列は $V$ で収束する。従って弱式を満たし、スペクトルから一意である。$L^2$ から $L^2$ への逆作用素のノルムは $\sup_j|\lambda_j-\zeta|^{-1}$ である。これは Parseval から上界が従い、最も近い固有値の固有ベクトルを右辺にすれば達成できる。$\zeta$ が固有値に近づくと増幅は固有値までの距離の逆数になる。一方、その固有空間に直交する右辺に限れば、他の固有値との間隔によってこの発散を避けられる。

:::example[同じ大きさの入力に対する増幅の差]
$(0,\pi)$ の $L=-\partial_{xx}$ で $\zeta=4+1/10$ とする。$f=\sin2x$ のとき $u=-10\sin2x$ になるが、$f=\sin x$ のときは $u=-10\sin x/31$ である。二つの右辺の $L^2$ ノルムは等しい。増幅の差は $\zeta$ と $4,1$ との距離が異なるためである。
:::

変分ノルムでも展開できる。$e_j/\sqrt{\lambda_j+\gamma}$ は $(V,a_\gamma)$ で正規直交である。これらに直交する $v$ は $(v,e_j)_2=0$ が全ての $j$ で成立するため $v=0$ である。従って $a_\gamma(v,v)=\sum_j(\lambda_j+\gamma)|(v,e_j)_2|^2$ である。逆に右辺が有限な係数列は $V$ で収束するため、この等式でエネルギー空間の元を特徴付けられる。

:::example[区間の固有値]
$(0,L)$ の Dirichlet 作用素 $-\partial_{xx}+q$、$q$ は定数とする。ODE と両端条件により $\lambda_j=(j\pi/L)^2+q$、$e_j=(2/L)^{1/2}\sin(j\pi x/L)$ となる。例えば $L=\pi$、$q=-2$ なら最初の二つは $-1,2$ である。
:::

:::example[逆作用素の固有値の順序]
$(0,\pi)$ の Dirichlet Laplace 作用素に $3I$ を足して逆を取る。$\sin x,\sin2x,\sin3x$ での逆作用素の固有値は $1/4,1/7,1/12$ である。もとの作用素で小さい固有値が、コンパクトな逆作用素で大きい固有値になる。例えば $f=\sin x+7\sin2x$ の解に $K$ を使うと $Kf=\sin x/4+\sin2x$ である。
:::

### $L^2$ 上の作用素の定義域

$V$ 上の形式と $L^2$ 上の閉作用素を区別する。後者の定義域は $D(L)=\{u\in V:\exists f\in L^2,\ a(u,v)=(f,v)_2\quad\forall v\in V\}$ であり、その $f$ を $Lu$ とする。$f$ は $V$ の稠密性で一意である。固有展開から $f_j=\lambda_j u_j$ であり、逆に $\sum_j(1+\lambda_j^2)|u_j|^2<\infty$ なら $u\in V$ でこの等式を満たす $f$ が存在する。実際 $(\lambda_j+\gamma)/(1+\lambda_j^2)$ は有界なので $u$ の部分和は $V$ で収束し、$Lu$ の部分和も $L^2$ で収束する。$u_m\to u$、$Lu_m\to f$ の $L^2$ 収束があれば各係数で $f_j=\lambda_j u_j$ が従うから作用素は閉じている。

複素数の $\zeta$ が固有値の列に含まれなければ $\sup_j(1+\lambda_j^2)/|\lambda_j-\zeta|^2<\infty$ なので、$u_j=f_j/(\lambda_j-\zeta)$ が領域内の一意解を与える。従って $L^2$ スペクトルは固有値の列そのものである。複素 Hilbert 空間への拡張では内積を第一引数について線形とする。

この作用素は自己共役である。実際、定義域に属する二つの元について対称性が成り立つので $L\subset L^*$ である。逆に $w\in D(L^*)$ なら、固有関数 $e_j$ を試験して $(L^*w,e_j)_2=\lambda_j(w,e_j)_2$ を得る。左辺の二乗和は有限なので $\sum_j\lambda_j^2|(w,e_j)_2|^2<\infty$ であり、上の定義域の特徴付けから $w\in D(L)$ である。そのため同じスペクトル展開を、固有値問題だけでなく時間発展の作用素を定義するときにも使える。

:::example[重複固有値の固有空間に不変な量]
三つの固有方向のうち二方向で同じ固有値を持つモデルを考える。

剛性行列を $K$ とする。対角成分は座標方向での値を与える。固有方程式から直交する方向を取り出すため、特性多項式 $\det(\lambda I-K)$ を計算する。

$$
\begin{aligned}
K&=\left[\begin{matrix}\frac{3994}{845} & - \frac{1728}{845} & \frac{240}{169}\\- \frac{1728}{845} & \frac{2986}{845} & - \frac{180}{169}\\\frac{240}{169} & - \frac{180}{169} & \frac{463}{169}\end{matrix}\right]\\
\det(\lambda I-K)&=\lambda^{3} - 11 \lambda^{2} + 32 \lambda - 28=\left(\lambda - 7\right) \left(\lambda - 2\right)^{2}\\
d_1&=2,\quad v_1=\left[\begin{matrix}\frac{3}{5}\\\frac{4}{5}\\0\end{matrix}\right]\\
d_2&=2,\quad v_2=\left[\begin{matrix}- \frac{4}{13}\\\frac{3}{13}\\\frac{12}{13}\end{matrix}\right]\\
d_3&=7,\quad v_3=\left[\begin{matrix}\frac{48}{65}\\- \frac{36}{65}\\\frac{5}{13}\end{matrix}\right]
\end{aligned}
$$

$V$ の列に固有ベクトルを並べる。$V^\top V=I$ と $V^\top KV=D$ を使えば、元の座標 $X$ と固有座標 $Y=V^\top X$ の双方で基底係数を記述できる。二次形式は固有座標にすると交差項を持たない。

$$
\begin{aligned}
V&=\left[\begin{matrix}\frac{3}{5} & - \frac{4}{13} & \frac{48}{65}\\\frac{4}{5} & \frac{3}{13} & - \frac{36}{65}\\0 & \frac{12}{13} & \frac{5}{13}\end{matrix}\right]\\
V^\top KV&=\left[\begin{matrix}2 & 0 & 0\\0 & 2 & 0\\0 & 0 & 7\end{matrix}\right]\\
X^\top KX&=\frac{3994 x_{1}^{2} - 3456 x_{1} x_{2} + 2400 x_{1} x_{3} + 2986 x_{2}^{2} - 1800 x_{2} x_{3} + 2315 x_{3}^{2}}{845}\\
\sum_j d_jY_j^2&=2 Y_{1}^{2} + 2 Y_{2}^{2} + 7 Y_{3}^{2}
\end{aligned}
$$

一つの基底係数で固有係数を計算し、エネルギーの二次形式を求める。直交変換なので二乗ノルムは変わらない。二次形式を二乗ノルムで割ると、固有値を非負の重みで平均した値になる。

$$
\begin{aligned}
X&=\left[\begin{matrix}1\\2\\3\end{matrix}\right]\\
Y&=\left[\begin{matrix}\frac{11}{5}\\\frac{38}{13}\\\frac{51}{65}\end{matrix}\right]\\
|X|^2&=14\\
X^\top KX&=\frac{26261}{845}\\
\frac{X^\top KX}{|X|^2}&=\frac{26261}{11830}\\
\min_j d_j&=2\\
\max_j d_j&=7
\end{aligned}
$$

固有値 $2$ に属する二本は、その平面内の別の正規直交基底に取り替えられる。個々の係数は変化するが、その二乗和と二本の寄与の和は直交射影のため不変である。固有値の列では $2$ を二回数える。
:::

### 固有空間に射影したときのレゾルベント

固有値 $\lambda_*$ の直交射影を $P_*$、他の固有値との最短距離を $d>0$ とする。$0<|z-\lambda_*|<d/2$ なら
$$
(L-z)^{-1}=(\lambda_*-z)^{-1}P_*+(L-z)^{-1}(I-P_*)
$$
である。第二項の作用素ノルムは $2/d$ 以下である。実際他の固有値 $\lambda_j$ で $|\lambda_j-z|\ge|\lambda_j-\lambda_*|-|\lambda_*-z|\ge d/2$ だからである。従って $z\to\lambda_*$ において $(\lambda_*-z)(L-z)^{-1}\to P_*$ が作用素ノルムで成り立つ。この極限で、共鳴する成分の有無だけでなく、その固有空間への射影が回復する。一方、源が $P_*f=0$ ならこの極を通過しても、$u\perp\ker(L-\lambda_*)$ に限定した解へ収束する。

:::example[同じ固有値を持つ二つの空間形]
正方形 $(0,\pi)^2$ で $-\Delta$ の固有値 $5$ は $e_{12},e_{21}$ の二本に対応する。ここで $e_{mn}$ は正規直交化された $\sin(mx)\sin(ny)$ である。$f=2e_{12}+e_{21}+3e_{11}$ なら $(L-z)^{-1}f=(2e_{12}+e_{21})/(5-z)+3e_{11}/(2-z)$ である。例えば $z=9/2$ なら係数は $(4,2,-6/5)$。$(5-z)$ を掛けて $z\to5$ とすれば $2e_{12}+e_{21}$ が残り、二つの空間形を同時に含む射影が得られる。
:::

## Rayleigh 商と min-max 原理

:::definition[Rayleigh 商]{#def-pde-rayleigh-quotient}
$u\in V\setminus\{0\}$ に対し $R(u)=a(u,u)/\norm u_2^2$ を :term[Rayleigh 商]{en="Rayleigh quotient" yomi="れいりーしょう"} と呼ぶ。
:::

:::theorem[min-max 原理]{#thm-pde-elliptic-minmax}
:ref[thm-pde-elliptic-spectrum] の条件下で

$$
\lambda_k=\min_{\substack{S\subset V\\ \dim S=k}}\ 
\max_{u\in S\setminus\{0\}}R(u).
\label{eq-pde-elliptic-minmax}
$$

:::

:::proof
$S=\mathrm{span}(e_1,\dots,e_k)$ では $R(u)$ が $\lambda_1,\dots,\lambda_k$ の重み付き平均なので最大は $\lambda_k$ である。任意の $k$ 次元 $S$ では、$e_1,\dots,e_{k-1}$ に $L^2$ 直交する $u\ne0$ がある。スペクトル展開でその商は $\lambda_k$ 以上である。両方を合わせると証明が完了する。
:::

既知の低い固有モードを除外して下界を求めるときは、これと双対な式

$$
\lambda_k=\max_{\substack{F\subset L^2\\\dim F=k-1}}
\inf_{\substack{u\in V\setminus\{0\}\\u\perp_{L^2} F}}R(u)
$$

が使える。実際、$F=\mathrm{span}(e_1,\ldots,e_{k-1})$ を取れば下限は $e_k$ で達成される。任意の $F$ に対しては $\mathrm{span}(e_1,\ldots,e_k)$ 内に $F$ と直交する非零ベクトルが存在する。その商は $\lambda_k$ 以下なので、どの $F$ を選んでも下限は $\lambda_k$ を超えない。例えば $(0,\pi)$ で $\int u\sin x=0$ を課せば $\norm{u'}_2^2\ge4\norm u_2^2$ であり、$\sin2x$ が等号を与える。

Dirichlet 領域を広げると、零延長により小領域の試験部分空間を大領域で使えるため、同じ係数なら固有値は上昇しない。また二つの対称形式が点ごとに $a_1(u,u)\le a_2(u,u)$ を満たすなら、それぞれの固有値も同じ順序に並ぶ。

:::example[多項式を試験する]
$(0,1)$ の $-\partial_{xx}$ に対し $v=x(1-x)$ は $\int (v')^2=1/3$、$\int v^2=1/30$ を満たす。最低固有値は $\lambda_1\le10$ であり、正確な値 $\pi^2$ に近い。各試験関数で得られるのは上界である。
:::

:::example[二番目の固有値を二次元部分空間で挟む]
$(0,\pi)$ の $-\partial_{xx}$ を $S=\mathrm{span}(\sin x,\sin3x)$ 上で調べる。この部分空間上の Rayleigh 商の最大は $9$ なので $\lambda_2\le9$ である。$\sin3x$ を $\sin2x$ に替えれば上界は $4$ に下がる。さらに任意の二次元空間で $\sin x$ と直交する方向を選ぶと、残る正弦係数は $j\ge2$ だけなので $4$ が下界となる。
:::

### 比較による固有値の摂動評価

同じ定義域と同じ $L^2$ ノルムで二つの対称形式 $a,b$ を比べ、$|a(v,v)-b(v,v)|\le\delta\norm v_2^2$ とする。Rayleigh 商は各非零 $v$ について高々 $\delta$ だけ異なる。従って任意の $j$ 次元部分空間での最大値の差も高々 $\delta$ であり、部分空間に関する下限を取っても $|\lambda_j(a)-\lambda_j(b)|\le\delta$ である。$a(v,v)\le b(v,v)$ が全 $v$ で成り立てば同じ議論で $\lambda_j(a)\le\lambda_j(b)$ となる。

:::example[有界なポテンシャルでの固有値区間]
$(0,\pi)$ で $L=-\partial_{xx}+2+\cos x$ とする。ポテンシャルは $[1,3]$ にあるので $j^2+1\le\lambda_j\le j^2+3$ である。第 $1$ 固有値は $1+2=3$ 以下にもできる。実際 $v=\sin x$ に対して $\int_0^\pi\cos x\sin^2x\dd x=[\sin^3x/3]_0^\pi=0$ なので Rayleigh 商は $3$ である。定数束縛は全ての固有値に使え、試験関数は個別の上界を鋭くする。
:::

:::example[二番目の Rayleigh 値を与える平面]
三次元の剛性行列について、第2固有値の min-max に使う平面を選ぶ。

剛性行列を $K$ とする。対角成分は座標方向での値を与える。固有方程式から直交する方向を取り出すため、特性多項式 $\det(\lambda I-K)$ を計算する。

$$
\begin{aligned}
K&=\left[\begin{matrix}\frac{6538}{845} & - \frac{3636}{845} & \frac{336}{169}\\- \frac{3636}{845} & \frac{4417}{845} & - \frac{252}{169}\\\frac{336}{169} & - \frac{252}{169} & \frac{1020}{169}\end{matrix}\right]\\
\det(\lambda I-K)&=\lambda^{3} - 19 \lambda^{2} + 94 \lambda - 120=\left(\lambda - 12\right) \left(\lambda - 5\right) \left(\lambda - 2\right)\\
d_1&=2,\quad v_1=\left[\begin{matrix}\frac{3}{5}\\\frac{4}{5}\\0\end{matrix}\right]\\
d_2&=5,\quad v_2=\left[\begin{matrix}- \frac{4}{13}\\\frac{3}{13}\\\frac{12}{13}\end{matrix}\right]\\
d_3&=12,\quad v_3=\left[\begin{matrix}\frac{48}{65}\\- \frac{36}{65}\\\frac{5}{13}\end{matrix}\right]
\end{aligned}
$$

$V$ の列に固有ベクトルを並べる。$V^\top V=I$ と $V^\top KV=D$ を使えば、元の座標 $X$ と固有座標 $Y=V^\top X$ の双方で基底係数を記述できる。二次形式は固有座標にすると交差項を持たない。

$$
\begin{aligned}
V&=\left[\begin{matrix}\frac{3}{5} & - \frac{4}{13} & \frac{48}{65}\\\frac{4}{5} & \frac{3}{13} & - \frac{36}{65}\\0 & \frac{12}{13} & \frac{5}{13}\end{matrix}\right]\\
V^\top KV&=\left[\begin{matrix}2 & 0 & 0\\0 & 5 & 0\\0 & 0 & 12\end{matrix}\right]\\
X^\top KX&=\frac{6538 x_{1}^{2} - 7272 x_{1} x_{2} + 3360 x_{1} x_{3} + 4417 x_{2}^{2} - 2520 x_{2} x_{3} + 5100 x_{3}^{2}}{845}\\
\sum_j d_jY_j^2&=2 Y_{1}^{2} + 5 Y_{2}^{2} + 12 Y_{3}^{2}
\end{aligned}
$$

一つの基底係数で固有係数を計算し、エネルギーの二次形式を求める。直交変換なので二乗ノルムは変わらない。二次形式を二乗ノルムで割ると、固有値を非負の重みで平均した値になる。

$$
\begin{aligned}
X&=\left[\begin{matrix}2\\1\\3\end{matrix}\right]\\
Y&=\left[\begin{matrix}2\\\frac{31}{13}\\\frac{27}{13}\end{matrix}\right]\\
|X|^2&=14\\
X^\top KX&=\frac{14905}{169}\\
\frac{X^\top KX}{|X|^2}&=\frac{14905}{2366}\\
\min_j d_j&=2\\
\max_j d_j&=12
\end{aligned}
$$

$v_1,v_2$ が張る平面上での最大商は $5$ である。どの二次元平面でも $v_1$ に直交する非零ベクトルがあり、その商は $5$ 以上になる。従ってこの最大値を $5$ より小さくすることはできない。
:::

:::example[一つの試験方向に沿うエネルギー差]
$g=x(1-x)$ を用い、$J(ag)=\int_0^1|(ag)_x|^2/2-\int_0^1 3ag$ を最小化する。

まず候補の選択を一変数 $a$ の問題にする。停留条件は $J_a=0$ であり、境界や無限遠も許される場合は停留点と分けて扱う。次の式では一階微分を取って零になる点を求め、二階微分を残しておく。

$$
\begin{aligned}
J(a)&=\frac{a^{2}}{6} - \frac{a}{2}\\
J_a&=\frac{d}{d a} \left(\frac{a^{2}}{6} - \frac{a}{2}\right)=\frac{2 a - 3}{6}\\
J_{aa}&=\frac{d^{2}}{d a^{2}} \left(\frac{a^{2}}{6} - \frac{a}{2}\right)=\frac{1}{3}\\
a_*&=\frac{3}{2}\\
J_a(a_*)&=0,\qquad a=\frac{3}{2}
\end{aligned}
$$

停留値と、その点から $h$ だけ動かした差を計算する。一次項が相殺しているため、展開の先頭から局所の振る舞いが分かる。一方、停留点が全体で最小か最大かを結論するには、許容範囲全体でこの差の符号を調べる必要がある。

$$
\begin{aligned}
J(a_*)&=- \frac{3}{8},\qquad a=\frac{3}{2}\\
J_{aa}(a_*)&=\frac{1}{3}\\
J(a_*+h)-J(a_*)&=- \frac{h}{2} + \frac{\left(h + \frac{3}{2}\right)^{2}}{6} - \frac{3}{8}=\frac{h^{2}}{6}\\
\text{二次係数}&=\frac{1}{6}
\end{aligned}
$$

指定した値で目的関数を比較する。$a$ 自体の値と、そこでの目的関数の値は別の量である。

$$
\begin{aligned}
a_*\text{ の具体値}&=\frac{3}{2}\\
J_*\text{ の具体値}&=- \frac{3}{8}\\
J&=- \frac{1}{3},\qquad a=1\\
J&=- \frac{1}{3},\qquad a=2\\
\text{二つの値の差}&=0\\
\text{停留点まわりの差}&=\frac{h^{2}}{6}
\end{aligned}
$$

最小点の係数 $3/2$ から $u=(3/2)x(1-x)$ を得る。この場合は $-u_{xx}=3$ と一致するので、制限した方向で得た最小点が全体の解にもなる。同じ試験関数で $f=x$ を選ぶと最小化係数は得られても、$-u_{xx}$ は定数なので全体の解にはならない。
:::

:::example[Rayleigh–Ritz に現れる二種類の行列]
$(0,1)$ で $A=1$、$c=0$ とし、$u=x^3(1-x)$ を使う。二つの基底に関する近似を作り、同じ行列から固有値を近似する場合の式も読む。

剛性行列を作る前に、方程式と試験空間を決める。ここでは全ての基底が両端で消える。厳密解を含まない試験空間に制限すると、近似は割り当てた方向の弱形式を満たすように決まる。厳密解は誤差を測るときに使う。

$$
\begin{aligned}
A&=1\\
u&=x^{3} \left(1 - x\right)\\
f&=6 x \left(2 x - 1\right)\\
\psi_1&=x \left(1 - x\right)\\
\psi_2&=x \left(1 - x\right) \left(2 x - 1\right)
\end{aligned}
$$

剛性の成分は $K_{ij}=\int A\psi^{\prime}_j\psi^{\prime}_i+\int c\psi_j\psi_i$、右辺は $F_i=\int f\psi_i$ である。基底を正規直交には取っていないので、内積の質量行列 $M_{ij}=\int\psi_j\psi_i$ も別に求める。

$$
\begin{aligned}
K&=\left[\begin{matrix}\frac{1}{3} & 0\\0 & \frac{1}{5}\end{matrix}\right]\\
M&=\left[\begin{matrix}\frac{1}{30} & 0\\0 & \frac{1}{210}\end{matrix}\right]\\
F&=\left[\begin{matrix}\frac{1}{10}\\\frac{1}{10}\end{matrix}\right]\\
\det(K-\lambda M)&=\frac{\lambda^{2}}{6300} - \frac{13 \lambda}{1575} + \frac{1}{15}=\frac{\left(\lambda - 42\right) \left(\lambda - 10\right)}{6300}
\end{aligned}
$$

連立方程式 $Kc=F$ を解き、元の関数に戻す。残差 $r=f-Lu_N$ を計算すると一般に零ではないが、選んだ試験関数に対する積分は全て零になる。PDE を点で満たすことと有限次元の弱形式を満たすことを区別する。

$$
\begin{aligned}
c&=\left[\begin{matrix}\frac{3}{10}\\\frac{1}{2}\end{matrix}\right]\\
u_N&=- \frac{x \left(x - 1\right) \left(5 x - 1\right)}{5}\\
e&=u-u_N=- \frac{x \left(x - 1\right) \left(5 x^{2} - 5 x + 1\right)}{5}\\
r&=\frac{12 \left(5 x^{2} - 5 x + 1\right)}{5}\\
(\int r\psi_i)_i&=\left[\begin{matrix}0\\0\end{matrix}\right]
\end{aligned}
$$

最後に誤差の二乗ノルムを二通りに計算する。近似空間へのエネルギー直交性により、誤差の勾配エネルギーは全体のエネルギーから近似の分を引いた値になる。ここでは二次形式そのものを用いており、エネルギー汎関数の二分の一は付けていない。

$$
\begin{aligned}
a(u,u)&=\frac{3}{35}\\
a(u_N,u_N)&=\frac{2}{25}\\
a(e,e)&=\frac{1}{175}\\
\norm e_2^2&=\int\limits_{0}^{1} \frac{x^{2} \left(x - 1\right)^{2} \left(5 x^{2} - 5 x + 1\right)^{2}}{25}\, dx=\frac{1}{15750}
\end{aligned}
$$

近似固有値は $Kv=\lambda Mv$ で求める。上の特性式の二つの根は、二次元試験空間内での最小と最大の商である。全空間での min-max はより多くの部分空間を使えるので、試験空間の $j$ 番目の根は真の $j$ 番目の固有値の上界である。
:::

### 両側を用いた近似固有値の確認

対称な $L$ の定義域に属する $v\ne0$ と実数 $\mu$ に対し、固有展開を使うと
$$
\norm{(L-\mu)v}_2^2=\sum_j|\lambda_j-\mu|^2|v_j|^2
\ge \operatorname{dist}(\mu,\sigma(L))^2\norm v_2^2 .
$$
従って残差を $L^2$ で評価できる試行ベクトルなら、その商の近くに固有値があることを保証できる。一方、min-max は指定した番号の上界を与える。残差だけでは番号や重複度まで分からないので、番号を指定するには既知の上下界などを合わせる。例えば区間 $(0,\pi)$ で $v=\sin x+\eps\sin2x$ とすると Rayleigh 商は $(1+4\eps^2)/(1+\eps^2)$ である。$\eps=1/3$ では商は $13/10$、残差と $L^2$ ノルムの比は
$$
\left(\frac{(1-13/10)^2+(4-13/10)^2/9}{1+1/9}\right)^{1/2}=9/10 .
$$
よって固有値が $[2/5,11/5]$ にあり、この場合は第一固有値 $1$ が入る。試行方向の混合に対し商の誤差は $O(\eps^2)$ であるが、二階作用素を使う残差は $O(|\eps|)$ である。

## Weyl の漸近公式（概観）

:::theorem[Dirichlet Laplace 作用素の計数]{#thm-pde-weyl-laplacian}
有界開集合 $\Omega$ の境界が体積 $0$ であれば、$N(\Lambda)=\#\{j:\lambda_j\le\Lambda\}$ は

$$
N(\Lambda)\sim \frac{\omega_n|\Omega|}{(2\pi)^n}\Lambda^{n/2}
\qquad(\Lambda\to\infty)
\label{eq-pde-weyl}
$$

を満たす。ここで $\omega_n=|B_1(0)|$ である。
:::

:::proof
辺長 $h$ の立方体では分離変数により $\lambda=(\pi/h)^2|k|^2$ である。Dirichlet の場合は各 $k_i\ge1$、Neumann では各 $k_i\ge0$ である。直交基底の完全性は区間の正弦または余弦基底を各変数で順に使って従う。半径 $R=h\sqrt\Lambda/\pi$ の第 1 象限の格子点を単位立方体で数えると、半径 $R\pm\sqrt n$ の象限球の体積に挟まれ、座標面の点は $O(R^{n-1}+1)$ である。よって両境界条件の主要項は $\omega_n h^n\Lambda^{n/2}/(2\pi)^n$ である。

$\Omega$ の内部に閉包が含まれる立方体に Dirichlet 条件を置き、各関数を零延長すると、それらの直和は $H^1_0(\Omega)$ の部分空間である。従って min-max からその固有値数の和が $N(\Lambda)$ の下界になる。$\Omega$ の閉包と交わる立方体の和で外側を覆い、$H^1_0(\Omega)$ の元を全空間に零延長して各立方体に制限する。立方体ごとに独立な $H^1$ の直和では境界の一致を要求しないため、Neumann 固有値数の和が上界になる。$h$ を固定して $\Lambda\to\infty$ とし、次に $h\downarrow0$ とする。境界が体積 $0$ なので内外の立方体の体積はともに $|\Omega|$ に収束し、上下界が一致する。
:::

:::example[区間と正方形]
長さ $L$ では $N(\Lambda)=\lfloor L\sqrt\Lambda/\pi\rfloor$ で、$\omega_1=2$ による主要項と一致する。辺長 $\pi$ の正方形では $\lambda=m^2+n^2$ なので、主要項は第 1 象限の面積に対応する $\pi\Lambda/4$ である。
:::

:::example[領域を拡大したときの固有値数]
$\Omega$ を倍率 $s$ で拡大すると、座標変換で二階微分に $s^{-2}$ が掛かる。従って $N_{s\Omega}(\Lambda)=N_\Omega(s^2\Lambda)$ であり、主要係数は $s^n$ 倍となる。例えば二次元で正方形の辺を $2$ 倍にすると、同じ高い閾値までの固有値数の主要項は $4$ 倍である。
:::

### 計数から固有値の成長と熱跡を得る

$N(\Lambda)\sim C\Lambda^{n/2}$ なら $\lambda_j\sim(j/C)^{2/n}$ である。実際、十分大きな $\Lambda$ で $(1-\eps)C\Lambda^{n/2}\le N(\Lambda)\le(1+\eps)C\Lambda^{n/2}$ と挟む。$\lambda_j=\inf\{\Lambda:N(\Lambda)\ge j\}$ に両側の評価を適用し $\eps\downarrow0$ とすればよい。重複度が有限でも一様に有界とは限らないため $N(\lambda_j)=j$ を仮定しない。

非負の固有値について $\sum_j e^{-t\lambda_j}=t\int_0^\infty e^{-ts}N(s)\dd s$ である。各項を $t\int_{\lambda_j}^\infty e^{-ts}\dd s$ と書いて非負の積分を交換すると従う。$r=ts$ と変数変換し、$t^{n/2}N(r/t)\to Cr^{n/2}$ と $N(s)\le C'(1+s^{n/2})$ を用いると、$t^{n/2}\sum_j e^{-t\lambda_j}\to C\Gamma(1+n/2)$ である。従って Dirichlet Laplacian の Weyl 係数では極限は $|\Omega|(4\pi)^{-n/2}$ になる。これは熱核の対角の主要項と同じ体積係数である。

### 格子点の個数を面積で挟む

辺長 $a_i$ の直方体の Dirichlet 固有値は $\sum_i(\pi k_i/a_i)^2$（$k_i\ge1$）である。$N(\Lambda)$ は正象限で軸半径 $R_i=a_i\sqrt\Lambda/\pi$ の楕円体に入る格子点数である。各点 $k$ に単位立方体 $\prod_i[k_i-1,k_i]$ を対応させる。辺が座標軸に平行なので、正象限で各座標を減らせば楕円体に留まる。従って立方体の和は楕円体の体積以下である。逆に、十分大きな $\Lambda$ で $\delta<1$ とし、座標平面上の測度零の点を除く。楕円体ノルム $\|x\|_R=(\sum(x_i/R_i)^2)^{1/2}$ が $1-\delta$ 以下の点を考える。各成分を上へ丸めるとノルムは高々 $\delta=(\sum_iR_i^{-2})^{1/2}$ だけ増えるので、この点は対応する立方体の和に入る。これで
$$
(1-\delta)^n\frac{\omega_n\prod_iR_i}{2^n}
\le N(\Lambda)\le\frac{\omega_n\prod_iR_i}{2^n}
$$
を得る。$\delta=O(\Lambda^{-1/2})$ なので、直方体では誤差が $O(\Lambda^{(n-1)/2})$ 以下である。Neumann では軸上の格子点も数えるが、いずれかの座標が零の集合は一つ低い次元の楕円体に入り、同じ誤差で抑えられる。従って両方に同じ主要項が現れる。この同じ主要項を使うことが小立方体による Dirichlet と Neumann の挟み撃ちに必要である。

:::example[長方形の重複を含めた計数]
辺長が $\pi,2\pi$ なら固有値は $m^2+n^2/4$ である。$\Lambda=5$ では $m=1$ のとき $1\le n\le4$、$m=2$ のとき $1\le n\le2$、$m\ge3$ は不可能なので $N(5)=6$。特に $m=1,n=4$ と $m=2,n=2$ はどちらも固有値 $5$ であり、二つとも数える。無限の計数を面積へ移す前に、重複を値の集合へ縮めると誤る。
:::

### 逆作用素の固有値の和が収束する条件

計数公式から $\lambda_j\sim Cj^{2/n}$ なら、正の作用素 $L$ に対して $\sum_j\lambda_j^{-p}$ が有限であることと $p>n/2$ は同値である。実際、十分大きな $j$ では各項を定数倍の $j^{-2p/n}$ で上下から挟める。$p=n/2$ は調和級数に相当するので発散する。零固有値がある場合は逆をその直交補空間に限定するか、$L+I$ を使う。例えば区間 $(0,\pi)$ の Dirichlet Laplacian の逆では $p=1$ の和は有限だが、二次元の面積を持つ領域では $p=1$ の和が発散する。一方、熱半群の係数 $e^{-t\lambda_j}$ は各 $t>0$ で全ての多項式のべきよりも速く下がるので、同じ次元でも熱跡は有限である。

## 演習問題

::::exercise{level=basic}
$-\partial_{xx}+3$ を $(0,\pi)$ の両端 Dirichlet 条件で用いる。第 2 固有値を求めよ。

:::solution
$\sin2x$ の負の二階微分は $4\sin2x$ である。定数の零階項で固有値は $3$ 上がるので第 2 固有値は $7$ である。
:::
::::

::::exercise{level=basic}
$-\partial_{xx}-4$ の $(0,\pi)$ の Dirichlet 問題に右辺 $\sin2x$ を与えられるか判定せよ。

:::solution
同次核は $\sin2x$ であり、右辺との内積は $\int_0^\pi\sin^2(2x)=\pi/2$ である。Fredholm 条件に反するため解はない。
:::
::::

::::exercise{level=basic}
$u=x(1-x)$ の $(0,1)$ における $-\partial_{xx}$ の Rayleigh 商を求めよ。

:::solution
分子は $\int_0^1(1-2x)^2=1/3$、分母は $\int_0^1x^2(1-x)^2=1/30$ なので商は $10$ である。
:::
::::

::::exercise{level=basic}
$f=\divg F$、$F\in L^2$ のとき $\norm f_{V^*}$ を見積もれ。

:::solution
$\langle f,v\rangle=-\int F\cdot Dv$ なので Cauchy–Schwarz より絶対値は $\norm F_2\norm v_V$ 以下である。単位球上で上限を取って $\norm f_{V^*}\le\norm F_2$ を得る。
:::
::::

::::exercise{level=advanced}
$(0,1)$ で $-u''=\delta_{1/2}$、両端 $0$ を満たす弱解を求めよ。

:::solution
傾きを左 $1/2$、右 $-1/2$ に取ると連続な $u=\min(x,1-x)/2$ になる。積分で $\int u'v'=(v(1/2)-v(0))/2-(v(1)-v(1/2))/2=v(1/2)$ である。強圧性により一意である。
:::
::::

::::exercise{level=advanced}
二つの固有値 $\lambda_1=2,\lambda_2=5$ の正規直交固有関数で $v=2e_1+e_2$ とする。商を求めよ。

:::solution
交差項は固有式と直交性で消える。分子は $4\cdot2+5=13$、分母は $4+1=5$ なので $13/5$ である。
:::
::::

::::exercise{level=advanced}
$\Omega_1\subset\Omega_2$ なら Dirichlet Laplace 固有値が $\lambda_k(\Omega_2)\le\lambda_k(\Omega_1)$ になることを示せ。

:::solution
$\Omega_1$ の試験空間の零延長は $H^1_0(\Omega_2)$ に入り、$L^2$ と勾配ノルムを保つ。$\Omega_2$ の min-max ではその延長による $k$ 次元部分空間も選べるので、小領域の値以下になる。
:::
::::

::::exercise{level=advanced}
$L=-\partial_{xx}-4$、$(0,\pi)$ の Dirichlet 問題で $f=\sin x$ の全解を求めよ。

:::solution
$\sin x$ モードの固有値は $1-4=-3$ なので特解は $-\sin x/3$ である。同次解は $C\sin2x$ のみだから、全解は $-\sin x/3+C\sin2x$ である。
:::
::::

## 確認テスト

::::quiz
弱形式の移動で何が起きるか。

- [x] 二階発散のうち一階を試験関数の微分に移す
- [ ] 未知関数を必ず二度古典微分する
- [ ] 境界条件が全て不要になる

:::explanation
$\int(-\divg(ADu))v=\int ADu\cdot Dv$ であり、試験関数の零境界値で境界項が消える。$u$ の零境界値は $H_0^1$ に入れることで課している。
:::
::::

::::quiz{answer="7" tol="1e-6"}
$(0,\pi)$ の $-\partial_{xx}+3$ の第 2 Dirichlet 固有値を求めよ。

:::explanation
第 2 正弦の負の二階微分は $4$ 倍であり、定数項を加えて $4+3=7$ である。
:::
::::

::::quiz
対称な Dirichlet 形式で同次核がある場合に可解性が要求するものはどれか。

- [x] 右辺が核に直交する
- [ ] 右辺が核の正の倍数である
- [ ] 右辺が常に $0$ である

:::explanation
$w$ が核なら $a(u,w)=a(w,u)=0$ のため $\langle f,w\rangle=0$ が必要である。Fredholm の定理によりこの直交条件は十分でもある。
:::
::::

::::quiz{answer="10" tol="1e-6"}
$(0,1)$ の Dirichlet 作用素 $-\partial_{xx}$ に対する $x(1-x)$ の Rayleigh 商を求めよ。

:::explanation
勾配の二乗積分が $1/3$、関数の二乗積分が $1/30$ なので商は $10$ である。
:::
::::

::::quiz{answer="6" tol="1e-6"}
辺長 $\pi,2\pi$ の長方形の Dirichlet Laplace 作用素について、固有値 $5$ 以下の個数を重複込みで求めよ。

:::explanation
固有値は $m^2+n^2/4$ である。$m=1$ のとき $n\le4$、$m=2$ のとき $n\le2$ であり、それ以上の $m$ は不可能である。したがって $4+2=6$ となる。値 $5$ は二つの独立な固有関数から出るので二つ数える。
:::
::::


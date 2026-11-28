---
title: 非線形偏微分方程式の手法
description: 単調性、固定点、比較、対称性、エネルギーの散逸を用いて非線形偏微分方程式を解析する。
notation:
  - symbol: '\langle A(u),v\rangle'
    meaning: 作用素の値と試験関数の双対積
  - symbol: '\underline u,\overline u'
    meaning: 劣解と優解
  - symbol: '\mathcal E'
    meaning: 勾配流で減少するエネルギー
---

## この章の目標

- 単調作用素の弱極限を同定して解を作る。
- 固定点と優解・劣解の反復で存在を証明する。
- 固有関数による評価と積分恒等式で爆発や非存在を示す。
- 移動平面で正値解の対称性を導き、勾配流の散逸を調べる。

## 単調作用素の方法

線形の場合は有界な作用素が弱収束を保つ。非線形では $u_j\rightharpoonup u$ と $A(u_j)\rightharpoonup \chi$ から直接 $\chi=A(u)$ は従わない。そこで二つの点で作用素の差を取り、双対積の符号から極限を同定する。

:::definition[単調作用素]{#def-pde-monotone-operator}
実 Banach 空間 $V$ の作用素 $A:V\to V^*$ が $\langle A(u)-A(v),u-v\rangle\ge0$ を満たすとき :term[単調作用素]{en="monotone operator" yomi="たんちょうさようそ"} という。$u\ne v$ で厳密な不等号なら狭義単調である。$s\mapsto\langle A(u+sv),w\rangle$ が各 $u,v,w$ について連続なら半連続とする。
:::

:::lemma[Minty の極限同定]{#lem-pde-minty}
$A$ が単調かつ半連続で、$\langle\chi-A(v),u-v\rangle\ge0$ が全 $v\in V$ で成り立つなら $\chi=A(u)$ である。
:::

:::proof
$v=u-sw$、$s>0$ を代入し $s$ で割ると $\langle\chi-A(u-sw),w\rangle\ge0$ となる。$s\downarrow0$ で $\langle\chi-A(u),w\rangle\ge0$ であり、$-w$ も使うと等号になる。任意の $w$ なので結論を得る。
:::

:::theorem[単調性からの存在]{#thm-pde-browder-minty}
$V$ は可分反射的 Banach 空間とし、$A:V\to V^*$ は有界集合を有界集合に写し、単調、半連続で $\langle A(u),u\rangle/\norm u_V\to+\infty$（$\norm u_V\to\infty$）を満たすとする。任意の $f\in V^*$ に対し $A(u)=f$ が可解である。狭義単調なら解は一意である。
:::

:::proof
増加する有限次元部分空間 $V_j$ で稠密な合併を取る。各 $V_j$ への制限 $u\mapsto A(u)|_{V_j}$ は連続である。実際 $u_k\to u$ なら有界性で部分列の像の極限 $\chi$ が存在し、単調性から $v\in V_j$ について $\langle\chi-A(v),u-v\rangle\ge0$ を得る。:ref[lem-pde-minty] を $V_j$ 内で使えば $\chi=A(u)|_{V_j}$ であり、任意の収束部分列の極限が同じなので連続である。$V_j$ 上の Euclid 内積で $A(u)-f$ をベクトル $B(u)$ として表す。有限次元の閉球では :ref[thm-fa-schauder] の不動点定理を使える。これから、球面で $(B(u),u)>0$ なら球内に零点があることを示す。実際、零点がなければ連続写像 $u\mapsto -R B(u)/|B(u)|$ は閉球を球面に写し、不動点では内積が負になり矛盾する。有限次元ではノルムが同値なので Euclid 球の半径を大きくしてこの符号を実現できる。従って $\langle A(u_j),v\rangle=\langle f,v\rangle$（$v\in V_j$）を満たす $u_j$ がある。

$v=u_j$ を入れ強圧条件を使うと $\norm{u_j}_V$ は $j$ によらず有界である。反射性で $u_j\rightharpoonup u$、また $A(u_j)$ も弱収束部分列を持つ。その極限は稠密な合併上で $f$ と等しく、従って $f$ である。任意の $v\in V$ に対して単調性を展開すると

$$
0\le\langle f,u_j\rangle-\langle A(u_j),v\rangle-\langle A(v),u_j-v\rangle.
$$

極限から $\langle f-A(v),u-v\rangle\ge0$ を得て Minty の補題を使う。単調性が狭義なら二つの解の差の双対積が $0$ になる場合は差が $0$ である。
:::

:::example[$p$-Laplace の弱解]
有界開集合 $\Omega$ 上で $V=W^{1,p}_0(\Omega)$、$1<p<\infty$ とし、$\langle A(u),v\rangle=\int |Du|^{p-2}Du\cdot Dv$ とする。Hölder の不等式で $\norm{A(u)}_{V^*}\le\norm{Du}_p^{p-1}$ であり、$\langle A(u),u\rangle=\norm{Du}_p^p$ である。$p$ 乗ノルムの狭義凸性から、その勾配 $|\xi|^{p-2}\xi$ は狭義単調である。これは異なる点での二つの接平面の厳密な不等式を足せば従う。従って零トレースで作用素も狭義単調である。半連続性は $|Du+sDv|^{p-1}\le C_p(|Du|^{p-1}+|Dv|^{p-1})$（$|s|\le1$）と Hölder から従う。よって各 $f\in V^*$ に一意解がある。この場合は :ref[thm-pde-direct-method] からも同じ解を得られるが、単調作用素の定理はポテンシャルを持たない作用素にも使える。
:::

:::example[ポテンシャルを持たない単調作用素]
$\R^2$ で $A(x,y)=(x^3-y,y^3+x)$ とする。二点の差との内積で $(x,y)\mapsto(-y,x)$ の項は $0$ になり、三乗の二項は非負で異なる二点なら少なくとも一つは正である。さらに $\langle A(x,y),(x,y)\rangle=x^4+y^4\ge (x^2+y^2)^2/2$ なので強圧条件を満たす。交差微分が $-1$ と $1$ で異なるため、この写像は $C^2$ スカラーポテンシャルの勾配にはならない。
:::

### 強単調性と右辺への依存性

$A$ が $\langle A(u)-A(v),u-v\rangle\ge c\norm{u-v}_V^p$（$c>0,p>1$）を満たすとする。右辺 $f,g$ に対応する解を $u,v$ とすれば双対ノルムから
$$
\norm{u-v}_V\le c^{-1/(p-1)}\norm{f-g}_{V^*}^{1/(p-1)} .
$$
これは解の存在定理とは別に、二解の差を一回試験するだけで得られる。単調性が狭義であっても、一様に差を下から評価できなければ同じ定数付き評価は出ない。

:::example[三乗で逆写像の Hölder 指数が決まる]
$A(s)=s^3$ では $(u^3-v^3)(u-v)=(u-v)^2(u^2+uv+v^2)\ge(u-v)^4/4$ である。実際 $u^2+uv+v^2=(u-v)^2/4+3(u+v)^2/4$。従って $|u-v|\le 4^{1/3}|f-g|^{1/3}$ であり、$v=-u$ で等号となる。$g=0$ では実際 $u=f^{1/3}$ なので $f\to0$ で Lipschitz 定数は発散する。滑らかな作用素であっても、逆写像の正則性は単調性の強さに依存する。
:::

:::example[強単調性で逆写像のノルムを押さえる]
対称な線形写像を単調作用素と見なし、右辺への依存性に使う下界を求める。

線形化の対称行列を $A$ とする。対角成分は座標方向での値を与える。固有方程式から直交する方向を取り出すため、特性多項式 $\det(\lambda I-A)$ を計算する。

$$
\begin{aligned}
A&=\left[\begin{matrix}\frac{171}{25} & - \frac{72}{25}\\- \frac{72}{25} & \frac{129}{25}\end{matrix}\right]\\
\det(\lambda I-A)&=\lambda^{2} - 12 \lambda + 27=\left(\lambda - 9\right) \left(\lambda - 3\right)\\
d_1&=3,\quad v_1=\left[\begin{matrix}\frac{3}{5}\\\frac{4}{5}\end{matrix}\right]\\
d_2&=9,\quad v_2=\left[\begin{matrix}- \frac{4}{5}\\\frac{3}{5}\end{matrix}\right]
\end{aligned}
$$

$V$ の列に固有ベクトルを並べる。$V^\top V=I$ と $V^\top AV=D$ を使えば、元の座標 $X$ と固有座標 $Y=V^\top X$ の双方で解の差を記述できる。二次形式は固有座標にすると交差項を持たない。

$$
\begin{aligned}
V&=\left[\begin{matrix}\frac{3}{5} & - \frac{4}{5}\\\frac{4}{5} & \frac{3}{5}\end{matrix}\right]\\
V^\top AV&=\left[\begin{matrix}3 & 0\\0 & 9\end{matrix}\right]\\
X^\top AX&=\frac{3 \left(57 x_{1}^{2} - 48 x_{1} x_{2} + 43 x_{2}^{2}\right)}{25}\\
\sum_j d_jY_j^2&=3 Y_{1}^{2} + 9 Y_{2}^{2}
\end{aligned}
$$

一つの解の差で固有係数を計算し、単調性を測る二次形式を求める。直交変換なので二乗ノルムは変わらない。二次形式を二乗ノルムで割ると、固有値を非負の重みで平均した値になる。

$$
\begin{aligned}
X&=\left[\begin{matrix}3\\2\end{matrix}\right]\\
Y&=\left[\begin{matrix}\frac{17}{5}\\- \frac{6}{5}\end{matrix}\right]\\
|X|^2&=13\\
X^\top AX&=\frac{1191}{25}\\
\frac{X^\top AX}{|X|^2}&=\frac{1191}{325}\\
\min_j d_j&=3\\
\max_j d_j&=9
\end{aligned}
$$

$X^\top AX\ge3|X|^2$ なので $|AX|\ge3|X|$。右辺の変化を $f_1-f_2$、解の変化を $X$ とすれば $|X|\le|f_1-f_2|/3$ である。最小固有方向の摂動で等号となる。
:::

:::example[凸なエネルギーの勾配を単調方程式に使う]
$A(a)=a^3+2a$、$f=3$ に対して $A(a)=f$ を解く。汎関数は $J(a)=a^4/4+a^2-3a$ である。

まず候補の選択を一変数 $a$ の問題にする。停留条件は $J_a=0$ であり、境界や無限遠も許される場合は停留点と分けて扱う。次の式では一階微分を取って零になる点を求め、二階微分を残しておく。

$$
\begin{aligned}
J(a)&=\frac{a^{4}}{4} + a^{2} - 3 a\\
J_a&=\frac{d}{d a} \left(\frac{a^{4}}{4} + a^{2} - 3 a\right)=\left(a - 1\right) \left(a^{2} + a + 3\right)\\
J_{aa}&=\frac{d^{2}}{d a^{2}} \left(\frac{a^{4}}{4} + a^{2} - 3 a\right)=3 a^{2} + 2\\
a_*&=1\\
J_a(a_*)&=0,\qquad a=1
\end{aligned}
$$

停留値と、その点から $h$ だけ動かした差を計算する。一次項が相殺しているため、展開の先頭から局所の振る舞いが分かる。一方、停留点が全体で最小か最大かを結論するには、許容範囲全体でこの差の符号を調べる必要がある。

$$
\begin{aligned}
J(a_*)&=- \frac{7}{4},\qquad a=1\\
J_{aa}(a_*)&=5\\
J(a_*+h)-J(a_*)&=- 3 h + \frac{\left(h + 1\right)^{4}}{4} + \left(h + 1\right)^{2} - \frac{5}{4}=\frac{h^{2} \left(h^{2} + 4 h + 10\right)}{4}\\
\text{二次係数}&=\frac{5}{2}
\end{aligned}
$$

指定した値で目的関数を比較する。$a$ 自体の値と、そこでの目的関数の値は別の量である。

$$
\begin{aligned}
a_*\text{ の具体値}&=1\\
J_*\text{ の具体値}&=- \frac{7}{4}\\
J&=0,\qquad a=0\\
J&=2,\qquad a=2\\
\text{二つの値の差}&=2\\
\text{停留点まわりの差}&=\frac{h^{2} \left(h^{2} + 4 h + 10\right)}{4}
\end{aligned}
$$

二階微分は $3a^2+2\ge2$ なので、停留点が全体で唯一の最小点である。$J(a)=J(1)$ のもう一つの根があるかを調べる必要はなく、下からの二次評価で排除できる。二つの右辺で解を比較すれば、逆写像の Lipschitz 定数は $1/2$ 以下である。
:::

### 一次元の $p$-Laplace 方程式を流束から解く

$(0,L)$、$p>1$、零境界で $-(|u'|^{p-2}u')'=f$ とする。$f$ が積分可能なら流束 $q=C-F(x)$、$F(x)=\int_0^x f$ である。写像 $s\mapsto |s|^{p-2}s$ は連続で厳密に増加し、逆は $\Phi(q)=\sgn(q)|q|^{1/(p-1)}$。従って $u(x)=\int_0^x\Phi(C-F(y))\dd y$ であり、もう一端の条件は
$$
G(C)=\int_0^L\Phi(C-F(y))\dd y=0
$$
である。$G$ は支配収束により連続で、厳密単調かつ $C\to\pm\infty$ で $\pm\infty$ に行く。従って $C$ は一つ存在し、$u\in W_0^{1,p}$ である。流束の絶対連続性から部分積分して弱式を得る。逆に弱解の流束は同じ分布微分を持つので定数分しか異ならず、境界条件により同じ $C$ になる。この一次元の構成は、抽象的な存在証明における単調性と強圧性が何を決めていたかを直接示す。

:::example[勾配が零の点における二階微分]
$L=2,p=3,f=2$ とする。$F(x)=2x$ で $C=2$ を選ぶと勾配 $u'=\sgn(1-x)\sqrt{2|1-x|}$ は $x=1$ について反対称なので両端の値が一致する。積分すると
$$
u(x)=\frac{2\sqrt2}{3}\big(1-|1-x|^{3/2}\big).
$$
特に中心値は $2\sqrt2/3$ である。$u''$ の絶対値は中心近くで $1/\sqrt{2|1-x|}$ なので局所 $L^1$ に入るが $L^2$ には入らない。流束は $2-2x$ と滑らかであり、流束を勾配から求める写像の導関数が原点で零になるため、その逆写像の導関数が発散する。
:::

## 不動点定理による存在証明

半線形方程式 $Lu=N(u)$ を解くとき、線形問題を解く写像 $L^{-1}$ があれば、未知量は写像 $T=L^{-1}N$ の不動点になる。:term[不動点]{en="fixed point" yomi="ふどうてん"} とは $T(u)=u$ を満たす点である。

:::example[縮小による一意性]
$-u''=1+\lambda\sin u$、$u(0)=u(1)=0$ とする。$G$ を区間の Dirichlet 逆作用素とし、$T(u)=G(1+\lambda\sin u)$ と置く。$G$ は比較原理で正値性を保ち、$G1=x(1-x)/2$ なので $\norm{Gg}_\infty\le\norm g_\infty/8$ である。$\sin$ は Lipschitz 定数 $1$ を持つので $\norm{T(u)-T(v)}_\infty\le|\lambda|\norm{u-v}_\infty/8$ となる。$|\lambda|<8$ なら :ref[thm-anl-contraction] により唯一の連続不動点があり、右辺の連続性から $C^2$ になる。例えば $\lambda=4$ なら誤差は一反復ごとに最大 $1/2$ 倍である。
:::

:::example[縮小がなくてもコンパクトな写像]
同じ区間で $-u''=1+\tfrac12\sin(64u)$、零境界を考える。右辺は $1/2$ と $3/2$ の間にあるので $T(u)=G(1+\tfrac12\sin(64u))$ は $C([0,1])$ の閉凸集合 $0\le u\le3/16$ をそれ自身に写す。定数関数 $v=\delta>0$ と $0$ に対し、$\norm{T(v)-T(0)}_\infty/\delta=|\sin(64\delta)|/(16\delta)\to4$ なので、この集合で縮小写像ではない。$G$ の解は二階微分が一様に有界であり、両端の値が等しいため一点で一階微分が $0$ になる。従って一階微分も一様に有界で、像は Arzelà–Ascoli により相対コンパクトである。$T$ は一様ノルムで連続なので :ref[thm-fa-schauder] により解がある。非線形写像が確実に値を戻す範囲を見つけることが、この適用の要点である。
:::

### 事前評価から不動点の球を作る

連続コンパクト写像 $T$ が Banach 空間 $X$ で定義されるとする。$u=\theta T(u)$（$0\le\theta\le1$）の全解が $\norm u<R$ に入ると仮定する。半径 $R$ の閉球への放射状の写像 $P_R z=z$（$\norm z\le R$）、$P_Rz=Rz/\norm z$（$\norm z>R$）を使う。$P_RT$ は閉球を自身に写す連続コンパクト写像であり、:ref[thm-fa-schauder] で不動点 $u$ がある。もし $\norm{T(u)}>R$ なら $u=\theta T(u)$、$\norm u=R$ となって仮定に反する。従って $T(u)=u$ である。この証明により、像全体の上界を与える代わりに、候補方程式 $u=\theta T(u)$ の解だけを評価すればよい。

:::example[増大する反応を含む方程式の一様評価]
$(0,1)$ で $-u''=\theta(1-u^3)$、零境界とする。負の最小点では $-u''\le0$ だが $\theta>0$ なら右辺は正であり矛盾する。$1$ を超える最大点では逆に左辺が非負、右辺が負である。従って $0\le u\le1$ であり、$\theta=0$ では零解しかない。$T(u)=G(1-u^3)$ は有界集合を $C^1$ の有界集合へ写すので $C([0,1])$ のコンパクト写像である。例えば $R=2$ として上の方法を使えば $\theta=1$ の解が得られる。$1-u^3$ は減少するので、二解の差を試験すれば一意性も従う。
:::

:::example[固定点反復の線形段階で何を解くか]
$(0,1)$ の $-u_{xx}+u$ を線形部とし、ある反復段階の源が $u=x(1-x)(1+x^2)$ に対応するとする。試験空間を二方向に制限して線形ステップを計算する。

剛性行列を作る前に、方程式と試験空間を決める。ここでは全ての基底が両端で消える。厳密解を含まない試験空間に制限すると、近似は割り当てた方向の弱形式を満たすように決まる。厳密解は誤差を測るときに使う。

$$
\begin{aligned}
A&=1\\
u&=x \left(1 - x\right) \left(x^{2} + 1\right)\\
f&=- x^{4} + x^{3} + 11 x^{2} - 5 x + 2\\
\psi_1&=x \left(1 - x\right)\\
\psi_2&=x \left(1 - x\right) \left(2 x - 1\right)
\end{aligned}
$$

剛性の成分は $K_{ij}=\int A\psi^{\prime}_j\psi^{\prime}_i+\int c\psi_j\psi_i$、右辺は $F_i=\int f\psi_i$ である。基底を正規直交には取っていないので、内積の質量行列 $M_{ij}=\int\psi_j\psi_i$ も別に求める。

$$
\begin{aligned}
K&=\left[\begin{matrix}\frac{11}{30} & 0\\0 & \frac{43}{210}\end{matrix}\right]\\
M&=\left[\begin{matrix}\frac{1}{30} & 0\\0 & \frac{1}{210}\end{matrix}\right]\\
F&=\left[\begin{matrix}\frac{10}{21}\\\frac{43}{420}\end{matrix}\right]\\
\det(K-\lambda M)&=\frac{\lambda^{2}}{6300} - \frac{3 \lambda}{350} + \frac{473}{6300}=\frac{\left(\lambda - 43\right) \left(\lambda - 11\right)}{6300}
\end{aligned}
$$

連立方程式 $Kc=F$ を解き、元の関数に戻す。残差 $r=f-Lu_N$ を計算すると一般に零ではないが、選んだ試験関数に対する積分は全て零になる。PDE を点で満たすことと有限次元の弱形式を満たすことを区別する。

$$
\begin{aligned}
c&=\left[\begin{matrix}\frac{100}{77}\\\frac{1}{2}\end{matrix}\right]\\
u_N&=- \frac{x \left(x - 1\right) \left(154 x + 123\right)}{154}\\
e&=u-u_N=- \frac{x \left(x - 1\right) \left(154 x^{2} - 154 x + 31\right)}{154}\\
r&=- \frac{154 x^{4} - 308 x^{3} - 1663 x^{2} + 1817 x - 370}{154}\\
(\int r\psi_i)_i&=\left[\begin{matrix}0\\0\end{matrix}\right]
\end{aligned}
$$

最後に誤差の二乗ノルムを二通りに計算する。近似空間へのエネルギー直交性により、誤差の勾配エネルギーは全体のエネルギーから近似の分を引いた値になる。ここでは二次形式そのものを用いており、エネルギー汎関数の二分の一は付けていない。

$$
\begin{aligned}
a(u,u)&=\frac{851}{1260}\\
a(u_N,u_N)&=\frac{14437}{21560}\\
a(e,e)&=\frac{1121}{194040}\\
\norm e_2^2&=\int\limits_{0}^{1} \frac{x^{2} \left(x - 1\right)^{2} \left(154 x^{2} - 154 x + 31\right)^{2}}{23716}\, dx=\frac{19}{304920}
\end{aligned}
$$

この一段階の精度は残差のうち試験空間の外にある成分で決まる。固定点の方程式を近似している場合、各段階で同じ $K$ を使えるが、右辺 $F$ は反復ごとに変わる。有限次元での固定点が元の PDE の解へ近づくには、試験空間の稠密性と反復または近似解の一様評価が別に必要である。
:::

### 縮小写像で局所解を構成できる時間

$C_0([0,L])$ を両端で零となる連続関数の空間とし、$S(t)$ をその上の Dirichlet 熱半群とする。$\norm{S(t)}\le1$ とする。$N(0)=0$ で、$N$ が有界区間上で Lipschitz と仮定する。$X=C([0,T];C_0([0,L]))$ に一様ノルムを入れ、積分写像
$$
(\mathcal T u)(t)=S(t)g+\int_0^t S(t-s)N(u(s))\dd s
$$
を考える。$R>\norm g_\infty$ とし、$M_R=\max_{|z|\le R}|N(z)|$、$L_R$ を同じ範囲の Lipschitz 定数とする。$\norm{\mathcal T u}_X\le\norm g_\infty+TM_R$、$\norm{\mathcal T u-\mathcal T v}_X\le TL_R\norm{u-v}_X$ なので、$TM_R\le R-\norm g_\infty$ かつ $TL_R<1$ なら閉球上で唯一の不動点がある。$N(u)$ は境界で零であるため各時刻で $C_0$ に入り、半群の強連続性から積分は時刻に連続である。これは境界を保ったまま反復を定義するために使った仮定である。

二つの局所解が同じ有界集合内にあれば、積分式を引いて短区間で同じ縮小を使い、繰り返すことで重なる区間全体で一致する。さらに $\sup_{t<T_*}\norm{u(t)}_\infty\le R_0$ なら、どの時刻 $t_0<T_*$ から始めても半径 $R=2R_0+1$ の同じ時間幅を使える。半群則により再始動した解は重なりで元の解に一致し、$t_0$ を十分 $T_*$ に近く選ぶとそれを超えられる。よって有限最大時刻には $L^\infty$ ノルムが有界であり続けることはない。

:::example[三乗反応の場合に安全な局所時間を決める]
$N(z)=z^3$、$\norm g_\infty=1$ なら $R=2$ で $M_R=8$、$L_R=12$ である。$T=1/24$ を取ると写像の像は半径 $1+8/24=4/3$ に入り、縮小率は $12/24=1/2$ である。$u^{(0)}=S(t)g$ から始めると最初の差は $\norm{u^{(1)}-u^{(0)}}_X\le T$（$\norm{u^{(0)}}\le1$）である。よって $k$ 回目以後の差の和は $2^{1-k}T$ 以下になる。例えば $k=3$ なら $1/96$ 以下である。ここでの $T$ は縮小を保証する下界であり、真の存在区間の長さと一致する必要はない。
:::

## 優解と劣解の方法

:::definition[優解と劣解]{#def-pde-upper-lower-solutions}
$Lu=f(x,u)$ について全ての非負試験関数に対し $L\underline u\le f(x,\underline u)$ が弱形式で成り立ち、境界のトレースが指定値以下なら、$\underline u$ を :term[劣解]{en="subsolution" yomi="れっかい"} という。二つの不等号を逆にしたものを :term[優解]{en="supersolution" yomi="ゆうかい"} という。
:::

:::theorem[順序を保つ反復]{#thm-pde-monotone-iteration}
有界 Lipschitz 領域で $L=-\divg(A D\cdot)$、$A$ は一様楕円かつ有界とする。斉次 Dirichlet 問題に対して $H^1\cap L^\infty$ に属する順序付きの劣解と優解 $\underline u\le\overline u$ があるとする。$f(x,s)$ は $x$ に可測、$s$ に連続で、両関数の値域を含む区間で一様有界であるとする。さらに $s\mapsto f(x,s)+Ms$ がこの区間で非減少になる $M\ge0$ があるとする。このとき両関数の間に最小の弱解と最大の弱解がある。
:::

:::proof
$(L+M)T(v)=f(x,v)+Mv$ を零境界で解く。:ref[thm-pde-elliptic-coercive] で可解であり、:ref[prop-pde-weak-elliptic-maximum] で $T$ は順序を保つ。劣解の不等式とトレース条件から $\underline u\le T(\underline u)$、優解から $T(\overline u)\le\overline u$ である。従って $u_k=T^k(\underline u)$ は増加し、$\overline u$ 以下である。点ごとの極限 $u$ があり、優収束で $u_k\to u$ と $f(x,u_k)\to f(x,u)$ が $L^2$ で成り立つ。線形エネルギー評価で $u_{k+1}$ は $H^1_0$ でも $T(u)$ に収束するため $u=T(u)$ である。任意の両者の間の解 $w$ は $T(w)=w$ なので、$\underline u\le w$ から各 $u_k\le w$ が従う。従って $u$ は最小である。$\overline u$ から始める減少列で最大解も得られる。
:::

:::example[ロジスティック型の境界値問題]
$(0,\pi)$ で $-u''=\lambda u(1-u)$ とし、$\lambda>1$ とする。$\underline u=\delta\sin x$ では $-\underline u''=\underline u\le\lambda\underline u(1-\underline u)$ が $0<\delta\le1-1/\lambda$ で成り立つ。定数 $\overline u=1$ は優解であり、境界で指定値 $0$ 以上である。例えば $\lambda=2$ なら $\delta=1/2$ を選べる。$f'(s)=\lambda(1-2s)\ge-\lambda$（$0\le s\le1$）なので $M=\lambda$ で反復が使え、正の内部解が得られる。
:::

:::example[非単調な反応をずらす]
$(0,1)$ で $-u''=1-u^2$ を考える。$0$ は劣解、定数 $1$ は優解である。反応 $1-s^2$ は $[0,1]$ で減少するが $1-s^2+2s$ は非減少である。従って $M=2$ を使う。最初の反復は $(-\partial_{xx}+2)u_1=1$、零境界なので $u_1=(1-\cosh(\sqrt2(x-1/2))/\cosh(\sqrt2/2))/2$ である。内部で $u_1>0$ となり、反復の極限も正である。
:::

### 区間での比較と幅の縮小

$-u''=f(x,u)$、区間 $(0,L)$、零境界について、優解・劣解が与える範囲で $f(x,z)-f(x,w)\le K(z-w)$（$z\ge w$）を仮定する。二解の差の正部分を $p$ とすると $\int |p'|^2\le K\int p^2$。従って $K<(\pi/L)^2$ なら Poincaré の不等式から $p=0$ である。二解を交換して一致を得る。単調反復そのものは上下から最小解・最大解を与えるが、両者が一致するにはこのような比較が必要である。

:::example[優解と劣解を多項式で選ぶ]
$(0,1)$ で $f(x,u)=1+u^2$ を考える。$\underline u=0$ は劣解、$\overline u=x(1-x)$ は $-\overline u''=2$、$1+\overline u^2\le1+1/16<2$ なので優解である。$0\le w\le z\le1/4$ で $f(z)-f(w)=(z+w)(z-w)\le (z-w)/2$。従ってこの順序区間で解は唯一である。第一反復は $-u_1''=1$ なので $u_1=x(1-x)/2$、中心で $1/8$。優解の値 $1/4$ とこの値の間に真の中心値が入る。大きい値の解まで一意と結論するには、順序区間に全ての解が入ることを別に示す必要がある。
:::

:::example[優解を探すために線形問題を解く]
$(0,2)$ で $-u_{xx}=1$、両端の値を零とする。得られた $u$ を非線形方程式の比較関数として使う。

まず流束 $q=Au_x$ を未知量にする。$q_x=-f$ を積分すると、任意の定数 $q(0)$ について $q(x)=q(0)-\int_0^x f$ である。さらに $q/A$ を積分し、指定された両端の値または平均に合わせて積分定数を決める。以下では計算後の $u$ について関係を逆に確かめる。

$$
\begin{aligned}
u(x)&=- \frac{x^{2}}{2} + x\\
u_x&=\frac{d}{d x} \left(- \frac{x^{2}}{2} + x\right)=1 - x\\
q&=1 - x\\
-f&=\frac{d}{d x} \left(1 - x\right)=-1
\end{aligned}
$$

次に値の条件と流出量を両端で区別して評価する。区間の外法線は左端で $-1$、右端で $1$ である。このため源項の積分は左側の流束から右側の流束を引いたものになる。流束を係数で割った積分は両端の値の差に等しい。

$$
\begin{aligned}
u(0)&=0,\qquad x=0\\
u(L)&=0,\qquad x=2\\
q(0)&=1,\qquad x=0\\
q(L)&=-1,\qquad x=2\\
\int_0^L f\,\mathrm{d}x&=\int\limits_{0}^{2} 1\, dx=2\\
\int_0^L (q/A)\,\mathrm{d}x&=\int\limits_{0}^{2} \left(1 - x\right)\, dx=0
\end{aligned}
$$

弱形式も具体的に比べる。試験関数 $\phi=x(L-x)$ は両端で零であり、積分の間に境界項は出ない。以下の値が一致するのは、任意の零トレース試験関数に対しても部分積分できることの一例である。

$$
\begin{aligned}
\phi&=x \left(2 - x\right)\\
\phi^{\prime}&=\frac{d}{d x} x \left(2 - x\right)=- 2 \left(x - 1\right)\\
\int_0^L q\phi^{\prime}\,\mathrm{d}x&=\int\limits_{0}^{2} \left(1 - x\right) \left(2 - 2 x\right)\, dx=\frac{4}{3}\\
\int_0^L f\phi\,\mathrm{d}x&=\int\limits_{0}^{2} x \left(2 - x\right)\, dx=\frac{4}{3}
\end{aligned}
$$

最後に内部での値、勾配、流束を同じ点で測る。区間全体での勾配の寄与には $A$ の重みが付くので、点での値だけでは全体のエネルギーは分からない。

$$
\begin{aligned}
u&=\frac{1}{2},\qquad x=1\\
u_x&=0,\qquad x=1\\
q&=0,\qquad x=1\\
\frac12\int_0^L A(u_x)^2\,\mathrm{d}x&=\int\limits_{0}^{2} \frac{\left(1 - x\right)^{2}}{2}\, dx=\frac{1}{3}
\end{aligned}
$$

この解は $0\le u\le1/2$ を満たす。例えば $-\partial_{xx}w=1/2+w^2$ では $0$ が劣解、$u$ が優解である。実際 $1/2+u^2\le3/4<1$ である。線形問題から、反復が入る順序区間を具体的に定められる。
:::

### 上下の反復の差から解の誤差を読む

$T$ が順序を保ち、$\underline u\le T\underline u$、$T\overline u\le\overline u$ とする。$u_k=T^k\underline u$、$v_k=T^k\overline u$ なら、帰納的に $u_k\le u_{k+1}\le v_{k+1}\le v_k$ である。従って真の固定点 $u$ がこの順序区間にある限り、位置ごとの誤差は $0\le u-u_k\le v_k-u_k$ である。さらに線形比較が $\norm{Gf}_\infty\le C_G\norm f_\infty$ を与え、反応がこの区間で Lipschitz 定数 $L_f$ を持つとき、$T=Gf$ なら
$$
\norm{v_k-u_k}_\infty\le (C_GL_f)^k\norm{\overline u-\underline u}_\infty .
$$
ここで $T$ が順序を保つために $f$ の増加性を使い、幅を定量化するために小さな Lipschitz 定数を使った。前者だけでは減少の速度や極限の一致までは分からない。

:::example[非線形反応を加えた区間の二回目の下界]
$(0,1)$ で $-u''=1+u^2$、$\underline u=0$、$\overline u=x(1-x)$ とする。$u_1=x(1-x)/2$ なので、$-u_2''=1+x^2(1-x)^2/4$ である。二回積分し両端を零にすると
$$
u_2=\frac{x(1-x)}2+\frac{x}{240}-\frac{x^4}{48}+\frac{x^5}{40}-\frac{x^6}{120}.
$$
実際、追加項の二階微分は $-x^2/4+x^3/2-x^4/4$ であり、追加項の $x=1$ での和は零である。中心値は $1/8+11/7680$。右辺の Lipschitz 定数は $1/2$、零境界の逆作用素のノルムは $1/8$ なので、二回反復後の幅は初めの最大幅 $1/4$ に $1/16^2$ を掛けた $1/1024$ 以下である。
:::

### 劣解・優解による正値解の厳密な範囲

$(0,\pi)$ で $-u''=\lambda u-u^3$、零境界とする。$\lambda>1$ なら $\underline u=a\sin x$、$0<a\le\sqrt{\lambda-1}$ は劣解である。実際 $-\underline u''=\underline u$、右辺は $\underline u(\lambda-\underline u^2)\ge\underline u$ である。優解 $\overline u=\sqrt\lambda$ は定数であり、右辺は零、境界で非負である。$f'(u)=\lambda-3u^2\ge-2\lambda$ なので、$M=2\lambda$ を足した右辺は単調になり、順序反復が使える。境界値は反復ごとに零であり、優解自身が境界で零である必要はない。従ってこの区間では正値解が得られる。

逆に $\lambda\le1$ で非負非零解があれば、$\sin x$ との弱い内積で $(\lambda-1)\int u\sin x=\int u^3\sin x$ を得る。左は非正、右は正なので矛盾する。従って正値解が存在する範囲はちょうど $\lambda>1$ である。例えば $\lambda=5$ なら $\underline u=2\sin x$、上側の定数は $\sqrt5$ であり、中心で $2\le u(\pi/2)\le\sqrt5$ と挟める。閾値は立方項の大きさだけでなく、零境界と最初の空間モードから決まる。

### 正値解の比を試験する一意性

$(0,L)$ で $-u''=f(u)$、斉次 Dirichlet 条件とする。$f(0)=0$ で局所 Lipschitz、$f(s)/s$ は $s>0$ で厳密に減少と仮定する。$u,v$ は閉区間で $C^2$ な正値解とする。いずれかの端で微分が零なら $u(0)=u'(0)=0$ または右端の同じ条件を持つことになり、ODE の一意性から零解になってしまう。従って両解は端に向かって一次で消え、$u^2/v,v^2/u$ は連続に端で零となり、勾配も有界である。例えば $u^2/v$ の微分は $2(u/v)u'-(u/v)^2v'$ であり、比が有界なので弱い試験関数に使える。

$v$ の方程式を $u^2/v$ で試験し、$u$ の方程式を $u$ で試験する。差を取ると
$$
\int_0^L u^2\left(\frac{f(u)}u-\frac{f(v)}v\right)
=\int_0^L\left(u'-\frac uv v'\right)^2\ge0 .
$$
$u,v$ を交換した式を足せば左辺は $\int (u^2-v^2)(f(u)/u-f(v)/v)$ となる。厳密減少性で被積分関数は $u\ne v$ の点で負であるから、両解は同一でなければならない。この方法では $f$ そのものの単調性を仮定しない。例えば $f(s)=\lambda s-s^3$ では $f(s)/s=\lambda-s^2$ が厳密減少なので、先の劣解・優解で作った正値解は唯一である。零解も同じ境界値を持つため、ここでの一意性は正値解についての結論である。

## 解の爆発と非存在

:term[爆発]{en="blow-up" yomi="ばくはつ"} は、有限時刻までに解の指定されたノルムが有界でいられない現象である。空間平均を取ると一部の拡散が失われるが、適切な重みを選べば拡散の効果を正確に一次項にできる。

:::proposition[固有関数による爆発の判定]{#prop-pde-kaplan-blowup}
滑らかな連結有界領域の斉次 Dirichlet 条件で $u_t-\Delta u=u^p$、$p>1$、$u\ge0$ の古典解を考える。第一固有関数 $\phi_1\ge0$ を $\int\phi_1=1$ と規格化し、固有値を $\lambda_1>0$ とする。$y_0=\int u(0)\phi_1>(2\lambda_1)^{1/(p-1)}$ なら、有界な古典解は

$$
T_*=\frac{2 y_0^{1-p}}{p-1}
\label{eq-pde-kaplan-time}
$$

を越えて存在できない。
:::

:::proof
$y=\int u\phi_1$ とする。部分積分で境界項は消え、$y'=-\lambda_1y+\int u^p\phi_1$ である。Jensen より $y'\ge-\lambda_1y+y^p$ である。初期値が閾値より大きいので最初から $y'\ge y^p/2>0$ であり、初めて閾値に戻ることはできない。よって $(y^{1-p})'\le-(p-1)/2$ であり、積分すると $t=T_*$ で左辺の値が非正になって矛盾する。$y\le\norm u_\infty$ であるから、有界性の破れは $\norm u_\infty$ の上界にも及ぶ。古典的な局所理論で有界解を継続できるクラスでは、この上限までに $\norm u_\infty$ が発散する。
:::

### 局所解と継続できる条件

連続な初期値に対し、$L^\infty$ の有界性から継続できることを、固定境界の一次元の場合に示す。区間 $(0,\pi)$ で端点零の連続関数の空間を $X$ とする。反射による解作用素 $S(t)$ は :ref[thm-pde-heat-maximum] で $\norm{S(t)}_{X\to X}\le1$、奇周期延長と熱核から $S(t)g\to g$ が一様に成り立つ。実数関数 $N(s)=(s_+)^p$ は局所 Lipschitz で $N(0)=0$ である。そこで

$$
\mathcal T u(t)=S(t)g+\int_0^t S(t-s)N(u(s))\,\dd s
\label{eq-pde-semilinear-mild}
$$

を $C([0,T];X)$ に作用させる。$R>\norm g_\infty$ とし、$\norm u_{C(X)}\le R$ の球上では $\norm{\mathcal T u}\le\norm g+TR^p$、二つの像の差は $Tp R^{p-1}\norm{u-v}$ 以下である。$T$ を小さくすれば球内の縮小になる。従って有界な連続解を一意に得る。熱核の非負性と比較により $S(t)$ は非負性を保つので、$g\ge0$ なら反復でも非負の解を得る。

解を始め直す時間幅はその時点の一様ノルムの上界にだけ依存する。従って有限の最大時刻 $T_{\max}$ に向けたある時刻列で $\norm{u(t_j)}_\infty$ が有界なら、共通の長さの局所解を十分遅い $t_j$ から作り、一意性でつなげて $T_{\max}$ を越えられる。これは矛盾なので、有限の最大時刻では一様ノルムが発散する。固有係数に :eqref[eq-pde-semilinear-mild] を適用すれば第一モードの ODE が得られるため、Kaplan の判定はこの緩和解にもそのまま使える。

第一固有関数を非負に取れる理由は、Rayleigh 商の最小化で $\phi$ を $|\phi|$ に置き換えても商が変わらないことである。$|\phi|$ も最小点なので弱固有式を満たす。次の例では固有関数を明示できる。

:::example[爆発時刻の上界]
$(0,\pi)$ で $p=2$、$u(0,x)=4\sin x$ とする。$\phi_1=\sin x/2$、$\lambda_1=1$ なので $y_0=\pi>2$ である。上の評価から時刻 $2/\pi$ を越える有界な解はない。拡散を無視して点ごとに ODE を解くとこの重み付きの上界は得られない。
:::

比較する ODE をそのまま解くと、粗い判定より閾値と上界を改善できる。$z'=-\lambda_1 z+z^p$、$z(0)=y_0>0$ とし、$w=z^{1-p}$ と置く。$w'=(p-1)(\lambda_1 w-1)$ だから

$$
w(t)=\lambda_1^{-1}+
(y_0^{1-p}-\lambda_1^{-1})e^{(p-1)\lambda_1 t}.
$$

$y_0>\lambda_1^{1/(p-1)}$ ならこれは有限時刻で零になる。$y$ に対しても $y>0$ の範囲で $(y^{1-p})'\le(p-1)(\lambda_1 y^{1-p}-1)$ となり、同じ初期値から積分因子で $y^{1-p}\le w$、従って $y\ge z$ である。この比較が有効な間は下から正の $z$ で抑えられるので、$y$ は零へ達しない。これにより爆発時刻の上界

$$
\frac{1}{(p-1)\lambda_1}\log\frac{1}{1-\lambda_1 y_0^{1-p}}
$$

を得る。例えば $p=2,\lambda_1=1,y_0=4$ では $\log(4/3)$ である。これは粗い上界 $1/2$ より小さい。閾値以下ではこの微分不等式自体は爆発を強制しない。小さい初期値の大域存在には、後で示すように比較できる優解を別に作る。

:::proposition[Pohozaev 恒等式と臨界指数]{#prop-pde-pohozaev}
$n\ge3$、滑らかな有界領域 $\Omega$ が原点に関して厳密に星型、すなわち $x\cdot\nu>0$ が境界で成り立つとする。$p\ge(n+2)/(n-2)$ なら $-\Delta u=u^p$、$u|_{\partial\Omega}=0$ は正の $C^2(\overline\Omega)$ 解を持たない。
:::

:::proof
両辺に $x\cdot Du$ を掛ける。右辺は $x\cdot D(u^{p+1}/(p+1))$ なので積分は $-n\int u^{p+1}/(p+1)$ である。左辺は

$$
\int Du\cdot D(x\cdot Du)-\int_{\partial\Omega}\partial_\nu u(x\cdot Du)
=\frac{2-n}{2}\int|Du|^2-\frac12\int_{\partial\Omega}(x\cdot\nu)(\partial_\nu u)^2
$$

となる。ここでは $Du\cdot D(x\cdot Du)=|Du|^2+x\cdot D(|Du|^2/2)$ と、境界で $Du=\nu\partial_\nu u$ を用いた。また $u$ を掛けて積分すれば $\int|Du|^2=\int u^{p+1}$ である。従って

$$
\frac12\int_{\partial\Omega}(x\cdot\nu)(\partial_\nu u)^2
=\left(\frac n{p+1}-\frac{n-2}{2}\right)\int u^{p+1}.
\label{eq-pde-pohozaev}
$$

右辺は非正であるが、$-u$ に Hopf の補題を適用すると $\partial_\nu u<0$ であり、厳密星型の仮定により左辺は正である。矛盾する。
:::

:::example[小さい初期値の大域的上界]
$(0,\pi)$ の $u_t-u_{xx}=u^p$ で $0\le g\le\sin x$ とする。時間一定の $\overline u=\sin x$ は $-\overline u''=\sin x\ge(\sin x)^p$ を満たすため優解である。有界な時間区間で $u^p$ と $\overline u^p$ の差は有界係数と $u-\overline u$ の積であり、時間指数を掛けて係数の符号を調整すれば弱比較を適用できる。従って非負解は $u\le\sin x$ のままである。局所解を一様ノルムで継続できるクラスでは、有限時間爆発は起こらない。
:::

### 小さい初期値と大きい初期値での違い

$u_t-\Delta u=u^p$、$p>1$、非負、斉次 Dirichlet 条件とする。第一固有関数 $\phi_1$ を正で $\norm{\phi_1}_\infty=1$ に規格化し、$a_0^{p-1}<\lambda_1$ とする。$a'=-\lambda_1 a+a^p$、$a(0)=a_0$ の非負解は非増加であり、$a'\le-(\lambda_1-a_0^{p-1})a$ である。$W=a\phi_1$ は $W_t-\Delta W=a^p\phi_1\ge a^p\phi_1^p$ なので優解である。$0\le g\le a_0\phi_1$ なら、局所解を比較でこの下に保つことで $L^\infty$ の発散を排除できる。反対に第一固有関数で重みを付けた平均が大きいと、Kaplan の不等式が有限時間内の発散を与える。この二条件は初期値の異なる数量を測っている。

:::example[二乗反応の爆発時刻の上界]
$\lambda_1=1$、$p=2$、規格化した平均を $J(0)=2$ とする。比較する $z'=z(z-1)$、$z(0)=2$ は $1/z$ を使うと $(1/z)'=1/z-1$ となり、解は $z=1/(1-e^t/2)$ である。従って $t=\log2$ で爆発する。Jensen の不等式による微分不等式から $J\ge z$ なので、元の解の存在時刻もこの上界を超えられない。$J'$ の式に等号を入れたのは比較する補助 ODE であり、元の解の爆発時刻が厳密に $\log2$ とは限らない。
:::

### 方程式に関する Pohozaev の境界項の向き

星形領域では $x\cdot\nu\ge0$ である。$-\Delta u=u^p$、$u=0$ 上で $u>0$ が内部にあるとき、境界で接方向の微分は零だから $Du=(\partial_\nu u)\nu$ である。$P=(x\cdot Du)Du-|Du|^2x/2$ と置くと
$$
\divg P=(x\cdot Du)\Delta u+(1-n/2)|Du|^2,qquad
P\cdot\nu=(x\cdot\nu)(\partial_\nu u)^2/2 .
$$
第一の等式では $x\cdot D(|Du|^2)/2$ が両項から出て相殺する。次に $(x\cdot Du)\Delta u=-x\cdot D(u^{p+1}/(p+1))$ を積分する。境界で $u=0$ なので積分は $n\int u^{p+1}/(p+1)$。また方程式に $u$ を掛けると $\int|Du|^2=\int u^{p+1}$。従って境界項は $(n/(p+1)-(n-2)/2)\int u^{p+1}$ である。$n=3,p=6$ なら係数は $3/7-1/2=-1/14$ で負なので、左の非負性と両立しない。$p=5$ では係数が零で直接は符号が矛盾しないが、球のように $x\cdot\nu>0$ の部分があり Hopf が適用できれば境界積分が正であるため、やはり正値解は存在しない。

境界条件を変えたり、星形性を外したりすると、上の右辺で使った境界の符号や零値が変わる。そのため臨界指数の非存在は、指数だけから任意の領域と境界条件に従う結論ではない。

## 移動平面法と解の対称性

対称領域でも非線形方程式の一意性は一般にない。そこで解を対称移動して同じデータに戻すだけでは、もとの解と等しいとは言えない。その代わり、領域の端に近い小さい部分から比較を始める。

:::theorem[球内の正値解の対称性]{#thm-pde-moving-plane}
$u\in C^2(\overline B_R)$ が $-\Delta u=f(u)$、$u|_{\partial B_R}=0$、$u>0$（内部）を満たすとする。$f$ が $u$ の値域で Lipschitz 連続なら、$u$ は中心からの距離だけに依存し、半径の増加に対して狭義に減少する。
:::

:::proof
$0<\mu<R$ に対し $\Sigma_\mu=B_R\cap\{x_1>\mu\}$、$x^\mu=(2\mu-x_1,x')$、$w_\mu=u(x^\mu)-u(x)$ とする。反射点は球内にある。境界で $w_\mu\ge0$ であり、内部では $-\Delta w_\mu=c_\mu(x)w_\mu$、$|c_\mu|\le M=\operatorname{Lip}f$ である。$w_{\mu,-}=\max(-w_\mu,0)$ を試験すると $\norm{Dw_{\mu,-}}_2^2\le M\norm{w_{\mu,-}}_2^2$ を得る。任意の $q>2$ で $H_0^1(B_R)\hookrightarrow L^q(B_R)$ が成り立つものを固定する。零延長して Hölder と Sobolev を用いると $\norm{w_{\mu,-}}_2^2\le C|\{w_\mu<0\}|^{1-2/q}\norm{Dw_{\mu,-}}_2^2$ である。従って負集合の測度が十分小さければ勾配が $0$ である。$\mu$ が $R$ に近いときは $\Sigma_\mu$ 自体が小さいので $w_\mu\ge0$ から始められる。

この順序が保たれる平面を $0$ まで動かせることを示す。$0<\mu_0<R$ で初めて止まると仮定すると連続性から $w_{\mu_0}\ge0$ である。球面の一部の反射は真に内部に入るため、この差は恒等的に $0$ ではない。$w_{\mu_0}\ge0$ について $\Delta w_{\mu_0}-M w_{\mu_0}=-(c_{\mu_0}+M)w_{\mu_0}\le0$ であり、強最大値原理を $-w_{\mu_0}$ に適用して内部で正である。帽子の中のコンパクト集合 $K$ を選んで補集合の測度を十分小さくし、$K$ 上の正の下限を使うと $\mu<\mu_0$ に少し動かしても $K$ 上で正である。従って負の集合は $\Sigma_\mu\setminus K$ に含まれ、先の小測度評価で空になる。これは停止の仮定に反する。従って $\mu=0$ でも比較が成り立ち、反対側からも動かすと中心面で対称である。任意の方向に同じ議論を行えば放射対称になる。$0<r_1<r_2<R$ に対し $\mu=(r_1+r_2)/2$ の面で $r_2e_1$ を反射すると $r_1e_1$ になる。$\mu>0$ の比較は強最大値原理で厳密なので $u(r_1)>u(r_2)$ である。
:::

:::example[対称性と明示解]
単位球で $f(u)=2n$ のとき $u=1-|x|^2$ は放射対称で狭義に減少する。$\mu>0$ での反射差は $x_1^2-(2\mu-x_1)^2=4\mu(x_1-\mu)$ であり、帽子の内部で正、平面上で $0$ になる。
:::

球の場合は半径方向の微分も $0<r<R$ で厳密に負である。実際、点 $re_1$ を通る $\mu=r$ の反射平面を使う。$w_\mu$ は帽子の内部で正、平面上で零であり、$\Delta w_\mu-M w_\mu\le0$ を満たす。$-w_\mu$ への Hopf の補題で、外法線 $-e_1$ に沿う微分 $\partial_{-e_1}(-w_\mu)$ が正となる。これは $\partial_1 w_\mu=-2\partial_1 u$ に等しいため、$u_r(r)<0$ が従う。

同じ正則性と境界条件を仮定すると、球での証明を一つの方向で対称な滑らかな有界凸領域にも使える。例えば $\{x_1=0\}$ に関して対称とする。固定した $x'$ に対する断面は空か区間 $(-a(x'),a(x'))$ である。$0<\mu<x_1<a(x')$ なら $-a(x')<2\mu-x_1<a(x')$ なので、帽子の反射は内部に入る。帽子は凸で連結であり、$a(x')>\mu$ を満たす断面の右端では反射先が内部に入るため、境界差の一部が厳密に正である。よって停止点での強最大値原理と小測度の比較がともに使える。両側から平面を動かすと、その方向での対称性が従う。

:::example[対称面が中心にないときに反射が保つもの]
球の中心を $(2,0,\dots,0)$ に移すと、中心対称面は $x_1=2$ になる。半径 $1$ の $u=1-(x_1-2)^2-|x'|^2$ に対し平面 $x_1=9/4$ で $x_1=5/2$ を写すと $2$ に移る。元の点と反射点の値の差は反射点を先に取って $1/4$ である。一般の式を使うときも座標の原点を球の中心に合わせる。
:::

### 負部分を排除する小領域評価

反射差 $w$ が $-\Delta w=c(x)w$ を満たし、境界で $w\ge0$ とする。$w_-\in H_0^1(D)$ を試験すると $\int|Dw_-|^2=\int c w_-^2$ である。$n>2$ で負集合を $E$ とすれば Hölder と Sobolev の不等式により
$$
\int_E w_-^2\le |E|^{2/n}\norm{w_-}_{2n/(n-2)}^2
\le C_n |E|^{2/n}\norm{Dw_-}_2^2.
$$
従って $C_n\norm{c_+}_\infty|E|^{2/n}<1$ なら負部分は零である。この評価によって、反射差が既に正のコンパクト集合を確保しておけば、平面を少し動かしても残りの小さい領域で負にならないことが分かる。一次元では区間の長さによる Poincaré の定数で同じ比較を行える。例えば幅 $L=1/2$、$c\le12$ なら第一固有値 $4\pi^2>12$ なので反射差の符号が保たれる。

### 一次元の正値減衰解を対称性から決める

全直線で $-u''+u=u^3$、$u>0$、$u,u'\to0$（$|x|\to\infty$）を考える。$u'$ を掛けると $(-(u')^2/2+u^2/2-u^4/4)'=0$ であり、減衰から定数は零である。従って
$$
(u')^2=u^2-u^4/2.
$$
正で無限遠に零へ行く連続関数は内点で正の最大を取る。その点 $x_0$ で $u'=0$ なので $u(x_0)=\sqrt2$。反射 $u(2x_0-x)$ は同じ ODE を満たし、$x_0$ で値と微分が一致する。右辺が滑らかな ODE の一意性から両者は一致する。さらに最大点から右では $u''(x_0)=\sqrt2-(\sqrt2)^3=-\sqrt2<0$ で下がり始める。$0<u<\sqrt2$ では一階微分は零でないので、右半直線で単調に減少する。

式 $u(x)=\sqrt2\,\mathrm{sech}(x-x_0)$ は $u''=u-u^3$ を満たす。実際 $q=\mathrm{sech}s$ に対し $q'=-q\tanh s$、$q''=q(1-2q^2)$ である。中心の値と微分が一致するため、一意性からすべての正値減衰解はこの形である。この一次元の議論は、対称中心を先に極大点として見つけ、反射後の解が同じ初期値を持つことを使う。多次元では一点の値と勾配だけでは楕円方程式の解は決まらないため、移動平面法では領域全体での差の符号が必要になる。

:::example[半値幅と中心エネルギー]
$u(x_0+s)=u(x_0)/2$ の条件は $\cosh s=2$ で、正の根は $s=\log(2+\sqrt3)$ である。従って左右の半値点の距離は $2\log(2+\sqrt3)$。中心では $u'=0$ かつ $u^2/2-u^4/4=1-1=0$ なので、空間の第一積分は零である。解の高さに応じて $u'$ の大きさが決まるため、半値点で中心と同じ傾きを持つと仮定することはできない。
:::

## 勾配流

:::definition[勾配流]{#def-pde-gradient-flow}
Hilbert 空間 $H$ 上でエネルギー $\mathcal E$ の微分を $(\nabla_H\mathcal E(u),v)_H$ で表せるとき、$u_t=-\nabla_H\mathcal E(u)$ を :term[勾配流]{en="gradient flow" yomi="こうばいりゅう"} という。
:::

連鎖律が使える解では

$$
\frac{\dd}{\dd t}\mathcal E(u(t))=(\nabla_H\mathcal E(u),u_t)_H=-\norm{u_t}_H^2 .
\label{eq-pde-gradient-flow-dissipation}
$$

例えば $\mathcal E=\int(|Du|^2/2+W(u))$、$H=L^2$、斉次 Dirichlet 条件なら $L^2$ 勾配は $-\Delta u+W'(u)$ であり、発展は $u_t-\Delta u+W'(u)=0$ である。部分積分により勾配を導くには、この式が $L^2$ に属する時刻を考えればよい。自由境界のエネルギーなら第一変分の自然条件が同時に課される。

:::proposition[凸エネルギーの縮小と時間離散化]{#prop-pde-gradient-contract}
$\mathcal E$ が凸で、勾配が定まる滑らかな二つの解 $u,v$ があるとき、$\norm{u-v}_H$ は非増加である。また $\mathcal E:H\to(-\infty,+\infty]$ が有限な値を取る点を少なくとも一つ持つ凸関数で、弱下半連続かつ下に有界とする。任意の $u^k\in H$ と時間刻み $\tau>0$ に対し

$$
u^{k+1}\in\operatorname*{arg\,min}_w\left(\mathcal E(w)+\frac{\norm{w-u^k}_H^2}{2\tau}\right)
\label{eq-pde-minimizing-movement}
$$

の最小化元は一意に存在し、その点で微分可能なら $(u^{k+1}-u^k)/\tau=-\nabla_H\mathcal E(u^{k+1})$ である。
:::

:::proof
凸性の二つの接平面の不等式を足すと $(\nabla\mathcal E(u)-\nabla\mathcal E(v),u-v)_H\ge0$ である。よって $\tfrac12(\norm{u-v}_H^2)'$ は非正である。離散式では被最小化汎関数を任意方向に微分して $0$ と置けば表示を得る。また $\mathcal E(u^k)<\infty$ のとき、$w=u^k$ と比較すると $\mathcal E(u^{k+1})+\norm{u^{k+1}-u^k}_H^2/(2\tau)\le \mathcal E(u^k)$ となる。最小化は例えば $\mathcal E$ が有界下、凸、弱下半連続なら直接法で実行でき、移動費用が狭義凸なので一意である。
:::

:::example[熱方程式の陰的時間刻み]
$-\Delta$ の固有値 $\lambda>0$ のモードについて $\mathcal E=\lambda a^2/2$ である。離散式は $a_{k+1}=a_k/(1+\tau\lambda)$ であり、例えば $\lambda=4$、$\tau=1/4$ で振幅は一歩で半分になる。$k\tau=t$ を固定して $\tau\downarrow0$ とすれば $(1+\tau\lambda)^{-k}\to e^{-t\lambda}$ で、熱半群と一致する。
:::

:::example[反応項のエネルギーと平衡点]
空間一定の $a$ に $\mathcal E(a)=(a^2-1)^2/4$ を取る。勾配流は $a'=a-a^3$ で、平衡点は $-1,0,1$ である。$a(0)=1/2$ では初速は $3/8>0$ であり、エネルギー変化率は $-9/64$ である。変化率の負号はエネルギーの非凸性にかかわらず成立するが、凸性から得た距離の縮小は、この例では一般に成り立たない。例えば初期値 $0$ と $1/2$ の間の距離は最初に増える。
:::

### 劣微分と射影による離散化

凸汎関数の :term[劣微分]{en="subdifferential" yomi="れつびぶん"} は $\xi\in\partial\mathcal E(u)$ を、全 $v$ に対する $\mathcal E(v)\ge\mathcal E(u)+(\xi,v-u)_H$ で定義する。近接最小化の点 $w$ では $(u-w)/\tau\in\partial\mathcal E(w)$ である。実際、最小化不等式を $w+s(v-w)$ に適用して凸性で $\mathcal E(w+s(v-w))$ を上から評価し、$s\downarrow0$ とすると劣微分の式になる。逆にその式を二乗項の展開と足せば、任意の $v$ へのエネルギー差は $\norm{v-w}^2/(2\tau)$ 以上である。従ってこれは必要十分条件である。

二つの近接点 $w_1,w_2$ について劣微分の二式を足すと、$\big((u_1-w_1)-(u_2-w_2),w_1-w_2\big)_H\ge0$ となる。従って $\norm{w_1-w_2}^2\le(u_1-u_2,w_1-w_2)_H$。距離が縮小するだけでなく、写像前後の差の角度にも制約がある。

:::example[絶対値エネルギーで零に達する]
$H=\R$、$\mathcal E(x)=|x|$ とする。$x\ne0$ の劣微分は $\operatorname{sgn}x$、零では $[-1,1]$ である。近接条件から $u>\tau$ なら $w=u-\tau$、$u<-\tau$ なら $w=u+\tau$、$|u|\le\tau$ なら $w=0$ である。例えば $u^0=5/4$、$\tau=1/2$ なら $3/4,1/4,0$ と進む。連続時間では $x(t)=(5/4-t)_+$ が包含式 $-x'\in\partial\mathcal E(x)$ を a.e. に満たす。滑らかでないエネルギーでも、零に着いた後に零から離れる必要はない。
:::

### 離散最小化から凸エネルギーの流れを作る

:::theorem[凸汎関数の劣微分による勾配流]{#thm-pde-subgradient-flow}
:ref[prop-pde-gradient-contract] の凸汎関数 $\mathcal E$ に対し、$\mathcal E(u_0)<\infty$ なら任意の $T>0$ で $u\in H^1(0,T;H)$ が存在し、

$$
u(0)=u_0,\qquad -u'(t)\in\partial\mathcal E(u(t))\quad\text{a.e.}
$$

を満たす。二解の距離は非増加なので、この解は一意である。
:::

:::proof
時刻 $T$ を等分する刻み $\tau$ を取り、:eqref[eq-pde-minimizing-movement] から列を作る。時刻 $k\tau$ の値が $u^k$ の区分線形補間を $U_\tau$、区間 $((k-1)\tau,k\tau]$ の値が $u^k$ の階段補間を $V_\tau$ とする。$\xi_\tau=-U_\tau'$ と置けば $\xi_\tau\in\partial\mathcal E(V_\tau)$ である。最小化で前の点を候補に取った不等式を足し合わせると

$$
\int_0^T\norm{\xi_\tau}^2\le 2(\mathcal E(u_0)-\inf\mathcal E)=M^2,\qquad
\sup_{0\le t\le T}\norm{U_\tau(t)-V_\tau(t)}\le M\sqrt\tau .
$$

第二の評価は一歩の移動の二乗が $\tau\int_{(k-1)\tau}^{k\tau}\norm{\xi_\tau}^2$ に等しいことから従う。二つの刻み $\tau,\sigma$ を比べる。劣微分は二つの支持不等式を足すと単調になるので、$\frac12(\norm{U_\tau-U_\sigma}^2)'=-(\xi_\tau-\xi_\sigma,U_\tau-U_\sigma)$ に $V_\tau-V_\sigma$ を加減すると

$$
\tfrac12(\norm{U_\tau-U_\sigma}^2)'
\le (\norm{\xi_\tau}+\norm{\xi_\sigma})M(\sqrt\tau+\sqrt\sigma)
$$

である。初期値が共通なので、積分と Cauchy–Schwarz から $U_\tau$ は $C([0,T];H)$ で Cauchy になる。その極限を $u$ とすると $V_\tau$ も一様に $u$ へ収束する。$\xi_\tau$ の $L^2(H)$ 弱収束部分列を取ると、時間の部分積分により弱極限は $\xi=-u'$ である。

極限で劣微分の式を示す。正の $\gamma$ を固定し、近接写像 $J_\gamma$ を :eqref[eq-pde-minimizing-movement] で定める。$z(t)=J_\gamma(u(t)+\gamma\xi(t))$ とし、$\eta=(u+\gamma\xi-z)/\gamma$ と置く。近接写像は非拡大的なので $z,\eta\in L^2(H)$ であり、各時刻で $\eta\in\partial\mathcal E(z)$ である。単調性から $\int_0^T(\xi_\tau-\eta,V_\tau-z)\ge0$ となる。強収束する $V_\tau$ と弱収束する $\xi_\tau$ を用いると、その極限は $\int_0^T(\xi-\eta,u-z)=-\gamma^{-1}\int_0^T\norm{u-z}^2$ である。従って $u=z$、$\xi=\eta\in\partial\mathcal E(u)$ が a.e. に成り立つ。

二解の差に対して劣微分の単調性を使うと $\frac12(\norm{u-v}^2)'\le0$ である。従って同じ初期値の解は一致する。さらに階段補間のエネルギー不等式を各時刻で下半連続性と弱下半連続性により極限に移すと、$\mathcal E(u(t))+\frac12\int_0^t\norm{u'}^2\le\mathcal E(u_0)$ となる。
:::

:::example[二次の凸勾配流で異なる減衰を測る]
$\mathcal E(X)=X^\top AX/2$ に対して勾配流と時間離散化を比べる。

線形化の対称行列を $A$ とする。対角成分は座標方向での値を与える。固有方程式から直交する方向を取り出すため、特性多項式 $\det(\lambda I-A)$ を計算する。

$$
\begin{aligned}
A&=\left[\begin{matrix}\frac{178}{25} & - \frac{96}{25}\\- \frac{96}{25} & \frac{122}{25}\end{matrix}\right]\\
\det(\lambda I-A)&=\lambda^{2} - 12 \lambda + 20=\left(\lambda - 10\right) \left(\lambda - 2\right)\\
d_1&=2,\quad v_1=\left[\begin{matrix}\frac{3}{5}\\\frac{4}{5}\end{matrix}\right]\\
d_2&=10,\quad v_2=\left[\begin{matrix}- \frac{4}{5}\\\frac{3}{5}\end{matrix}\right]
\end{aligned}
$$

$V$ の列に固有ベクトルを並べる。$V^\top V=I$ と $V^\top AV=D$ を使えば、元の座標 $X$ と固有座標 $Y=V^\top X$ の双方で解の差を記述できる。二次形式は固有座標にすると交差項を持たない。

$$
\begin{aligned}
V&=\left[\begin{matrix}\frac{3}{5} & - \frac{4}{5}\\\frac{4}{5} & \frac{3}{5}\end{matrix}\right]\\
V^\top AV&=\left[\begin{matrix}2 & 0\\0 & 10\end{matrix}\right]\\
X^\top AX&=\frac{2 \left(89 x_{1}^{2} - 96 x_{1} x_{2} + 61 x_{2}^{2}\right)}{25}\\
\sum_j d_jY_j^2&=2 Y_{1}^{2} + 10 Y_{2}^{2}
\end{aligned}
$$

一つの解の差で固有係数を計算し、単調性を測る二次形式を求める。直交変換なので二乗ノルムは変わらない。二次形式を二乗ノルムで割ると、固有値を非負の重みで平均した値になる。

$$
\begin{aligned}
X&=\left[\begin{matrix}1\\2\end{matrix}\right]\\
Y&=\left[\begin{matrix}\frac{11}{5}\\\frac{2}{5}\end{matrix}\right]\\
|X|^2&=5\\
X^\top AX&=\frac{282}{25}\\
\frac{X^\top AX}{|X|^2}&=\frac{282}{125}\\
\min_j d_j&=2\\
\max_j d_j&=10
\end{aligned}
$$

連続時間では固有係数は $e^{-2t},e^{-10t}$ 倍される。陰的刻みでは $\tau=1/2$ とすれば $1/2,1/6$ 倍である。いずれも高い固有値の方向を強く減らし、初期値を二つ変えても距離は増えない。
:::

:::example[非二次エネルギーの陰的な一歩]
$\mathcal E(a)=a^4/4$、初期点 $u^0=2$、刻み $\tau=1$ とする。$J(a)=\mathcal E(a)+|a-u^0|^2/2$ を最小にする。

まず候補の選択を一変数 $a$ の問題にする。停留条件は $J_a=0$ であり、境界や無限遠も許される場合は停留点と分けて扱う。次の式では一階微分を取って零になる点を求め、二階微分を残しておく。

$$
\begin{aligned}
J(a)&=\frac{a^{4}}{4} + \frac{\left(a - 2\right)^{2}}{2}\\
J_a&=\frac{d}{d a} \left(\frac{a^{4}}{4} + \frac{\left(a - 2\right)^{2}}{2}\right)=\left(a - 1\right) \left(a^{2} + a + 2\right)\\
J_{aa}&=\frac{d^{2}}{d a^{2}} \left(\frac{a^{4}}{4} + \frac{\left(a - 2\right)^{2}}{2}\right)=3 a^{2} + 1\\
a_*&=1\\
J_a(a_*)&=0,\qquad a=1
\end{aligned}
$$

停留値と、その点から $h$ だけ動かした差を計算する。一次項が相殺しているため、展開の先頭から局所の振る舞いが分かる。一方、停留点が全体で最小か最大かを結論するには、許容範囲全体でこの差の符号を調べる必要がある。

$$
\begin{aligned}
J(a_*)&=\frac{3}{4},\qquad a=1\\
J_{aa}(a_*)&=4\\
J(a_*+h)-J(a_*)&=\frac{\left(h - 1\right)^{2}}{2} + \frac{\left(h + 1\right)^{4}}{4} - \frac{3}{4}=\frac{h^{2} \left(h^{2} + 4 h + 8\right)}{4}\\
\text{二次係数}&=2
\end{aligned}
$$

指定した値で目的関数を比較する。$a$ 自体の値と、そこでの目的関数の値は別の量である。

$$
\begin{aligned}
a_*\text{ の具体値}&=1\\
J_*\text{ の具体値}&=\frac{3}{4}\\
J&=2,\qquad a=0\\
J&=4,\qquad a=2\\
\text{二つの値の差}&=2\\
\text{停留点まわりの差}&=\frac{h^{2} \left(h^{2} + 4 h + 8\right)}{4}
\end{aligned}
$$

陰的刻みでは $a-2=-a^3$ を満たす $a=1$ になる。これに対して陽的に一歩進むと $2-2^3=-6$ になり、エネルギーは $4$ から $324$ に増える。陰的な最小化の散逸は刻みが大きくても成立するが、陽的更新には同じ原理はない。
:::

### 強凸なエネルギーと指数収束

$\mathcal E$ が微分可能で $\mathcal E(v)\ge\mathcal E(u)+(D\mathcal E(u),v-u)_H+\lambda\norm{v-u}_H^2/2$（$\lambda>0$）を満たすとする。二点での不等式を足すと $\langle D\mathcal E(u)-D\mathcal E(v),u-v\rangle\ge\lambda\norm{u-v}^2$。従って二つの滑らかな流れの差は
$$
\frac12\frac{\dd}{\dd t}\norm{u-v}^2\le-\lambda\norm{u-v}^2,qquad
\norm{u(t)-v(t)}\le e^{-\lambda t}\norm{u(0)-v(0)}
$$
を満たす。最小化元 $v$ があれば $D\mathcal E(v)=0$ なので、定常解への指数収束になる。離散解では $w_1-w_2=(u_1-u_2)-\tau(D\mathcal E(w_1)-D\mathcal E(w_2))$ であり、$w_1-w_2$ と内積を取ると $(1+\tau\lambda)\norm{w_1-w_2}\le\norm{u_1-u_2}$ を得る。連続と離散で異なる係数だが、$k\tau=t$ で刻みを小さくすれば $(1+\tau\lambda)^{-k}\to e^{-\lambda t}$ となる。

:::example[Allen–Cahn 型のエネルギーで拡散が凸性を補う]
$(0,L)$ の零境界で $\mathcal E(u)=\int(|u'|^2/2+(u^2-1)^2/4)$ とする。二階変分は $\int(|v'|^2+(3u^2-1)v^2)$ である。$L=\pi/2$ なら第一 Dirichlet 固有値は $4$ なので、この二次形式は $3\norm v_2^2$ 以上である。線分に沿って二階微分を積分すれば $L^2$ の意味で強凸定数 $3$ が得られる。$u=0$ は停留点なので最小化元であり、滑らかな勾配流は $\norm{u(t)}_2\le e^{-3t}\norm{u(0)}_2$ を満たす。点ごとのポテンシャルは非凸でも、固定端の勾配の項が全体の下界を与える場合がある。
:::

### エネルギーの差を用いる勾配流の評価

凸な $\mathcal E$ の滑らかな勾配流で、任意の固定した許容点 $v$ に対し接平面の不等式を使うと
$$
\frac12\frac{\dd}{\dd t}\norm{u(t)-v}^2
=-(D\mathcal E(u),u-v)\le \mathcal E(v)-\mathcal E(u(t)).
$$
$\mathcal E(u(t))$ は非増加なので、$0$ から $t$ まで積分し、左辺の終端ノルムを捨てると
$$
\mathcal E(u(t))-\mathcal E(v)
\le \frac{\norm{u(0)-v}^2}{2t}.
$$
従って強凸でない場合でも、$v$ を最小化元と選べれば最小値へ近づく速度を評価できる。距離の収束を得るには、エネルギー差で距離を支配できる条件が別に必要である。陰的刻みで同じ結論を得るには、最小化条件の劣微分を $u^{j+1}-v$ と組み合わせ、
$$
\frac{\norm{u^{j+1}-v}^2-\norm{u^j-v}^2+\norm{u^{j+1}-u^j}^2}{2\tau}
\le\mathcal E(v)-\mathcal E(u^{j+1})
$$
とする。左の最後の項は離散化の非負の寄与なので、和を取れば $k\tau$ における同じ上界が出る。例えば $\mathcal E(a)=a^4/4$ では $a_t=-a^3$、$a(0)=1$ から $(a^{-2})'=2$、$a(t)=(1+2t)^{-1/2}$。エネルギーの解そのものは $1/(4(1+2t)^2)$ であり、一般上界 $1/(2t)$ より速く減る。この一般評価は個々のエネルギーに対して鋭い漸近式を与えるものではない。

:::example[拘束集合への近接ステップ]
$H=\R^2$ で $\mathcal E(u)=0$（$|u|\le1$）、$+\infty$（$|u|>1$）とする。このエネルギーの近接点は単位円板への直交射影である。実際、$u\ne0$ に対し $|w-u|^2=|w|^2+|u|^2-2u\cdot w$ は、半径 $|w|$ を固定すると $w$ を $u$ と同じ向きにしたとき最小である。次に $0\le |w|\le1$ で $(|w|-|u|)^2$ を最小にすれば、$|u|>1$ のとき $w=u/|u|$ を得る。例えば $u=(3,4)$ は $w=(3/5,4/5)$ へ写る。$|u|\le1$ なら $w=u$ であり、反復してもそれ以上動かない。有限なエネルギーの初期値から始めれば定常流となるが、近接問題自体は制約外の点からでも定義できる。
:::

## 演習問題

::::exercise{level=basic}
$A(s)=s^3$ の狭義単調性を示せ。

:::solution
$(a^3-b^3)(a-b)=(a-b)^2(a^2+ab+b^2)$ である。第 2 因子は $(a+b/2)^2+3b^2/4$ であり、$a\ne b$ なら正である。
:::
::::

::::exercise{level=basic}
$-u''=1+\lambda\sin u$ の区間 $(0,1)$ の零 Dirichlet 条件で $\lambda=6$ とする。対応する不動点写像の縮小率の上界を求めよ。

:::solution
Dirichlet 逆は $G1=x(1-x)/2$、$\norm G_{\infty\to\infty}\le1/8$ を満たす。$\sin$ の Lipschitz 定数は $1$ なので縮小率は $6/8=3/4$ 以下である。
:::
::::

::::exercise{level=basic}
$u_t-\Delta u=u^2$ の重み付き平均について $y_0=4$、$\lambda_1=1$ のとき固有関数による粗い爆発時刻の上界を求めよ。

:::solution
$4>2\lambda_1$ なので $y'\ge y^2/2$ を使える。$(1/y)'\le-1/2$ を積分し、$1/y_0=1/4$ が消える時刻は $1/2$ である。
:::
::::

::::exercise{level=basic}
$\mathcal E(a)=3a^2/2$ の実数上の勾配流を初期値 $a(0)=2$ から解け。

:::solution
微分は $3a$ なので $a'=-3a$、$a=2e^{-3t}$ である。エネルギーは $6e^{-6t}$ であり、その微分は $-36e^{-6t}=-(a')^2$ になる。
:::
::::

::::exercise{level=advanced}
$(0,\pi)$ で $\lambda=3$ のロジスティック問題に $\underline u=\delta\sin x$ を使える最大の $\delta\in(0,1]$ を求めよ。

:::solution
内部で $-\underline u''=\underline u>0$ なので劣解条件を割ると $1\le3(1-\delta\sin x)$ である。$\sin x$ の最大は $1$ だから $\delta\le2/3$ が必要十分である。
:::
::::

::::exercise{level=advanced}
三次元の厳密星型領域で $-\Delta u=u^5$ の正解を Pohozaev 恒等式で判定せよ。

:::solution
体積項の係数は $3/6-(3-2)/2=0$ である。境界では Hopf の補題により $\partial_\nu u$ は厳密に負で、$x\cdot\nu>0$ だから積分が正になる。従って恒等式に矛盾し、正の古典解は存在しない。
:::
::::

::::exercise{level=advanced}
$u=1-|x|^2$ の球内で $x_1=1/4$ の反射差を点 $x=(1/2,0,\dots,0)$ で求めよ。

:::solution
反射点の第 1 座標は $2\mu-x_1=0$ であり、反射点で $u=1$、元の点で $u=3/4$ なので差は $1/4$ である。一般の反射差 $4\mu(x_1-\mu)$ からも同じ値を得る。
:::
::::

::::exercise{level=advanced}
滑らかな凸エネルギーが全ての $u,v$ で $(\nabla\mathcal E(u)-\nabla\mathcal E(v),u-v)\ge m\norm{u-v}^2$ を満たすとき二解の距離を評価せよ。

:::solution
差の二乗の時間微分は $-2(\nabla\mathcal E(u)-\nabla\mathcal E(v),u-v)\le-2m\norm{u-v}^2$ である。積分因子により $\norm{u(t)-v(t)}\le e^{-mt}\norm{u(0)-v(0)}$ を得る。
:::
::::

## 確認テスト

::::quiz
Minty の補題で $u$ の近傍に選ぶ試験点はどれか。

- [x] $v=u-sw$ とし $s\downarrow0$
- [ ] 常に $v=0$ だけ
- [ ] $v$ を無限遠へ送る

:::explanation
$v=u-sw$ とすると双対積の $u-v$ が $sw$ になり、$s$ で割って半連続性を使える。$w$ と $-w$ の両方で不等式を得るため等号が従う。
:::
::::

::::quiz{answer="0.5" tol="1e-6"}
$\lambda_1=1$、$p=2$、$y_0=4$ の :eqref[eq-pde-kaplan-time] の値を求めよ。

:::explanation
式に代入すると $2\cdot4^{-1}/1=1/2$ である。これは平均が有限でいられる時間の上界である。
:::
::::

::::quiz
移動平面の最初の比較を可能にする条件はどれか。

- [x] 球の端の帽子の測度が十分小さい
- [ ] 方程式の全ての解が一意である
- [ ] 非線形項が常に定数である

:::explanation
反射差の負部分に Sobolev を適用すると、負集合の測度の正のべきがエネルギー評価に掛かる。帽子を小さく取るとその係数を $1$ 未満にでき、負部分が消える。
:::
::::

::::quiz{answer="0.5" tol="1e-6"}
熱方程式の陰的刻みで $\lambda=4$、$\tau=1/4$ の一歩の振幅比を求めよ。

:::explanation
係数方程式は $(a_{k+1}-a_k)/\tau=-\lambda a_{k+1}$ である。比は $(1+\tau\lambda)^{-1}=(1+1)^{-1}=1/2$ である。
:::
::::

::::quiz{answer="1" tol="1e-6"}
単調方程式 $a^3+2a=3$ の解を求めよ。

:::explanation
$a=1$ を代入すると $1+2=3$ になる。左辺は導関数 $3a^2+2>0$ を持つ狭義増加関数なので、それ以外の解はない。
:::
::::


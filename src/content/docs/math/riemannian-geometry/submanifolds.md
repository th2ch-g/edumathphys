---
title: 'Riemann 部分多様体と等長はめ込み'
description: '接続の接成分と法成分を分け、内在曲率、法曲率および体積の第一変分を結ぶ。'
notation:
  - symbol: '\mathrm{II}'
    meaning: '第二基本形式'
  - symbol: 'A_\xi'
    meaning: '法方向 ξ の形作用素'
  - symbol: '\nabla^\perp'
    meaning: '法接続'
  - symbol: 'H'
    meaning: '第二基本形式の平均'
---

## この章の目標

- 等長はめ込みの接・法分解を行う。
- Gauss、Codazzi、Ricci の方程式を導く。
- 全測地性と極小性を判定する。
- 超曲面で主曲率を計算する。

## 第二基本形式と形作用素

計量の引き戻しで接方向の長さを比較できる。次に、接方向に動いたときの接ベクトルの微分が法方向をどれだけ向くかを測る。曲面の場合の :ref[thm-cs-gauss-codazzi]{text="Gauss と Codazzi の方程式"} から、任意の次元と余次元へ議論を広げる。

:::definition[第二基本形式と形作用素]{#def-rg-second-fundamental}
等長はめ込み $F:(M,g)\to(\overline M,\bar g)$ をとり、引き戻した周囲接束を $TM\oplus NM$ と直交分解する。ここで $NM$ は法束である。:term[第二基本形式]{en="second fundamental form" yomi="だいにきほんけいしき"} は $\mathrm{II}(X,Y)=(\bar\nabla_XY)^\perp$ である。法ベクトル $\xi$ に対して $g(A_\xi X,Y)=\bar g(\mathrm{II}(X,Y),\xi)$ によって :term[形作用素]{en="shape operator" yomi="けいさようそ"} を定める。
:::

第 2 引数の $X(f)Y$ は接方向なので法射影で消え、第 1 引数は初めから関数線形である。また $\bar\nabla_XY-\bar\nabla_YX=dF[X,Y]$ は接方向なので $\mathrm{II}$ は対称である。従って $A_\xi$ は自己共役である。等長はめ込みは局所的に埋め込みであるが、自己交差では同じ周囲点に別の接空間が対応し得る。そこで分解は $M$ 上の引き戻し束で行う。

:::example[平面曲線]
単位速の $\gamma(s)$ では $\mathrm{II}(T,T)=\gamma''$ である。半径 $2$ の円で外向き法線を $N=\gamma/2$ とすると $\gamma''=-N/2$ なので $A_NT=-T/2$ となる。
:::

:::example[異なる法方向の情報]
$\R^4$ で $F(x,y)=(x,y,x^2,xy)$ を考える。原点では $g=I$、$e_3,e_4$ は法標構である。各成分の Hessian から $\mathrm{II}(e_1,e_1)=2e_3$、$\mathrm{II}(e_1,e_2)=e_4$、$\mathrm{II}(e_2,e_2)=0$ である。このとき $A_{e_3}=\diag(2,0)$ だが $A_{e_4}$ は対角が零、非対角成分が $1$ である。法方向を指定して初めて一つの形作用素が決まる。
:::

### 周囲加速度から測る法成分

曲線 $\gamma$ を $M$ 内に取り、$\dot\gamma(0)=v$ とする。周囲での速度は $dF(\dot\gamma)$ である。局所座標で $\bar\nabla_{\dot\gamma}dF(\dot\gamma)$ を展開すると、速度成分を微分する項 $\ddot x^i\,dF(\partial_i)$ は接方向である。残りの法部分は $(\dot x^i\dot x^j\bar\nabla_i dF(\partial_j))^\perp$ であり、$t=0$ では $\mathrm{II}(v,v)$ になる。このように、法加速度は曲線の接加速度と無関係で、初速度だけから決まる。異なる二方向の情報は
$$
2\mathrm{II}(u,v)=\mathrm{II}(u+v,u+v)-\mathrm{II}(u,u)-\mathrm{II}(v,v)
$$
で復元できる。従ってまず座標曲線や直線に沿う加速度を計算し、最後に対称双線形形式に戻す手順が使える。

:::example[グラフの法成分を直接求める]
$f:\R^m\supset U\to\R$ のグラフ $F(x)=(x,f(x))$ を考える。$W=\sqrt{1+|\grad f|_{\rm E}^2}$ とすると上向き単位法線は $N=(-\grad_{\rm E} f,1)/W$ である。実際、$F_i=(e_i,f_i)$ と内積を取ると $-f_i+f_i=0$ になる。$F_{ij}=(0,f_{ij})$ の $N$ 成分は $f_{ij}/W$ だから
$$
h_{ij}=f_{ij}/W,\qquad A_N=G^{-1}h,\qquad G=I+df\,df^\top
$$
である。最後では $df$ の成分を列とみなした。$f=x^2+2y^2$ の点 $(1,0)$ では $G=\diag(5,1)$、$W=\sqrt5$、$h=\diag(2,4)/\sqrt5$ なので形作用素の固有値は $2/(5\sqrt5),4/\sqrt5$ である。$h$ の固有値をそのまま読めるのは、接標構で $G=I$ の場合である。
:::

:::example[速度の尺度と単位方向]
半径 $3$ の円を $\gamma(t)=3(\cos2t,\sin2t)$ とたどる。$\dot\gamma$ の長さは $6$、$\ddot\gamma=-12N$ だから $\mathrm{II}(\dot\gamma,\dot\gamma)=-12N$ である。単位方向 $T=\dot\gamma/6$ については $\mathrm{II}(T,T)=-N/3$ となる。曲率を速度と独立に比べるには、法加速度を速度の二乗で割る。
:::

:::example[放物面の臍点の曲率の計算]
グラフ $F(x,y)=(x,y,f(x,y))$、$f=2 x^{2} + 2 y^{2}$ を点 $\left( 0, \  0\right)$ で調べる。

$F_x=(1,0,f_x)$、$F_y=(0,1,f_y)$ だから、計量行列は $G=I+df\,df^\top$ である。点での $df$ の列は $\begin{pmatrix}0\\0\end{pmatrix}$ なので、
$$
\begin{aligned}
G&=\begin{pmatrix}1 & 0\\0 & 1\end{pmatrix}\\
W=\sqrt{1+|df|^2}&=1\\
N=\frac{(-f_x,-f_y,1)}W&=\begin{pmatrix}0\\0\\1\end{pmatrix}
\end{aligned}
$$
となる。$N$ は分子の二乗長が $W^2$ であるため単位で、$N\cdot F_x=N\cdot F_y=0$ のため法線である。座標での面積要素の倍率は $\sqrt{\det G}=W$ である。実際、$df$ に垂直な方向の二乗長は $1$、平行方向は $1+|df|^2$ となる。$df=0$ の場合は両方とも $1$ である。

周囲微分 $F_{ij}$ は鉛直成分 $f_{ij}$ だけを持つ。第二基本形式の $N$ 成分は $F_{ij}\cdot N=f_{ij}/W$ である。従って
$$
\begin{aligned}
\bigl(h_{ij}\bigr)=\frac{\operatorname{Hess}_{\mathrm E} f}W&=\begin{pmatrix}4 & 0\\0 & 4\end{pmatrix}\\
A_N=G^{-1}h&=\begin{pmatrix}4 & 0\\0 & 4\end{pmatrix}
\end{aligned}
$$
となる。$G$ の逆を掛けるのは $g(A_N\partial_i,\partial_j)=h_{ij}$ により、二つの共変添字の一方を出力の添字に直すためである。この行列の特性多項式は $\det(\lambda I-A_N)=\lambda^{2} - 8 \lambda + 16$ であり、主曲率は $4,\quad 4$ である。主方向を単位化するときも、座標の通常のノルムではなく $G$ のノルムを使う。

平均曲率ベクトルは $(\tr A_N/2)N=4N$、Gauss 曲率は $\det A_N=16$ になる。同じ $K$ は $\det h/\det G=\det(\operatorname{Hess}_{\mathrm E}f)/W^4$ でも得られる。二つの $W$ は第二基本形式の法線への射影に、残る二つは Gram 行列式に由来する。法線を $-N$ に変えれば二つの主曲率は反転するので $K$ は変わらない。平均曲率の法成分は反転するが、掛ける法線も逆になるため平均曲率ベクトルは変わらない。

二つの固有値が一致するのでどの単位方向でも法曲率は同じである。
:::

## Gauss の公式と Weingarten の公式

:::proposition[Gauss と Weingarten の公式]{#prop-rg-gauss-weingarten}
接場 $Y$ と法場 $\xi$ に対し
$$
\bar\nabla_XY=\nabla_XY+\mathrm{II}(X,Y),\qquad
\bar\nabla_X\xi=-A_\xi X+\nabla_X^\perp\xi
\label{eq-rg-gauss-weingarten}
$$
である。$\nabla$ は誘導計量の Levi-Civita 接続、$\nabla^\perp$ は法束の計量接続である。
:::

:::proof
接射影は接続の Leibniz 則を満たし、計量両立式を接場の間で評価すると法成分が寄与しない。周囲の捩率が零で、Lie 括弧が接方向なので接射影の捩率も零である。従って :ref[thm-rg-levi-civita] により $\nabla$ となる。$X\bar g(\xi,Y)=0$ を微分すると $\bar g(\bar\nabla_X\xi,Y)=-\bar g(\xi,\mathrm{II}(X,Y))$ である。法射影 $\nabla^\perp$ は Leibniz 則を満たす。法場同士の内積を微分すると接成分は寄与しないため計量と両立する。
:::

:::example[球と円柱の形作用素]
$S^n(a)$ の外向き法線は $N=F/a$ だから $\bar\nabla_XN=X/a$ であり、$A_N=-\id/a$ である。円柱 $F(\theta,z)=(a\cos\theta,a\sin\theta,z)$ では $N=(\cos\theta,\sin\theta,0)$ なので、単位周方向に $A_N=-1/a$、軸方向に $A_N=0$ となる。法線を反転すると両方の作用素の符号が反転する。
:::

:::example[高さの方向を接方向から分ける]
$F(x,y)=(x,y,xy)$ の原点で、$\bar\nabla_{\partial_x}dF(\partial_y)=(0,0,1)$ である。原点の接空間は水平だから、このベクトルの接射影は零、法射影は $(0,0,1)$ である。一方、上向き法場 $\xi=(-y,-x,1)/\sqrt{1+x^2+y^2}$ の $\partial_x$ 微分は原点で $(0,-1,0)$ である。これは $-A_\xi\partial_x$ に等しく法接続の成分は零である。
:::

### 法標構の回転と法接続

二つの正規直交法場 $\nu_1,\nu_2$ を取れる領域を考える。法束が二次元で、$\nabla^\perp\nu_1=\omega\nu_2$ とすると計量両立性により $\nabla^\perp\nu_2=-\omega\nu_1$ である。新しい法標構を
$$
\nu'_1=\cos\phi\,\nu_1+\sin\phi\,\nu_2,\qquad
\nu'_2=-\sin\phi\,\nu_1+\cos\phi\,\nu_2
$$
とおく。Leibniz 則の係数をまとめると $\nabla^\perp\nu'_1=(\omega+d\phi)\nu'_2$ となる。一方、$A_{\nu'_1}=\cos\phi\, A_{\nu_1}+\sin\phi\,A_{\nu_2}$ であり、こちらには $d\phi$ が入らない。Weingarten の公式で接成分を取れば $X(\cos\phi)\nu_1+X(\sin\phi)\nu_2$ が消えるからである。

:::example[平面に沿って法標構だけを回す]
$\R^4$ の平面 $(x,y,0,0)$ で最初の法場を $e_3,e_4$、回転角を $\phi=2x-y$ とする。全ての $\mathrm{II}$ と $A$ は零であるが、$\nabla^\perp_{\partial_x}\nu'_1=2\nu'_2$、$\nabla^\perp_{\partial_y}\nu'_1=-\nu'_2$ となる。曲線 $(x,y)=(t,0)$ 上で一定の周囲ベクトル $e_3$ を法標構で表すと、成分は $(\cos2t,-\sin2t)$ である。この成分を微分した項は $\omega'=2\,dx-dy$ の接続項と相殺するので、$e_3$ は新しい標構でも法平行である。
:::

### 位置関数の Hessian と制限

周囲の滑らかな関数 $f$ を $M$ に制限する。Hessian の定義と Gauss の公式を用いれば
$$
\operatorname{Hess}_M(f\circ F)(X,Y)
=\operatorname{Hess}_{\overline M}f(dFX,dFY)+\bar g(\grad_{\overline M}f,\mathrm{II}(X,Y)).
$$
実際、両方の Hessian で $X(Y(f\circ F))$ は同じであり、引く接続の差が $\mathrm{II}$ である。Euclid 空間の線形関数 $f(z)=a\cdot z$ では最初の項が零になる。例えば球面 $S^m(R)$ で外向き法線 $F/R$ を選ぶと $\mathrm{II}=-gF/R^2$ だから、制限した高さ関数 $u=a\cdot F$ は $\operatorname{Hess}u=-u g/R^2$ を満たす。跡を取れば $\Delta u=-mu/R^2$ が得られる。

:::example[放物柱面の斜面の曲率の計算]
グラフ $F(x,y)=(x,y,f(x,y))$、$f=x^{2}$ を点 $\left( 1, \  0\right)$ で調べる。

$F_x=(1,0,f_x)$、$F_y=(0,1,f_y)$ だから、計量行列は $G=I+df\,df^\top$ である。点での $df$ の列は $\begin{pmatrix}2\\0\end{pmatrix}$ なので、
$$
\begin{aligned}
G&=\begin{pmatrix}5 & 0\\0 & 1\end{pmatrix}\\
W=\sqrt{1+|df|^2}&=\sqrt{5}\\
N=\frac{(-f_x,-f_y,1)}W&=\begin{pmatrix}- \frac{2 \sqrt{5}}{5}\\0\\\frac{\sqrt{5}}{5}\end{pmatrix}
\end{aligned}
$$
となる。$N$ は分子の二乗長が $W^2$ であるため単位で、$N\cdot F_x=N\cdot F_y=0$ のため法線である。座標での面積要素の倍率は $\sqrt{\det G}=W$ である。実際、$df$ に垂直な方向の二乗長は $1$、平行方向は $1+|df|^2$ となる。$df=0$ の場合は両方とも $1$ である。

周囲微分 $F_{ij}$ は鉛直成分 $f_{ij}$ だけを持つ。第二基本形式の $N$ 成分は $F_{ij}\cdot N=f_{ij}/W$ である。従って
$$
\begin{aligned}
\bigl(h_{ij}\bigr)=\frac{\operatorname{Hess}_{\mathrm E} f}W&=\begin{pmatrix}\frac{2 \sqrt{5}}{5} & 0\\0 & 0\end{pmatrix}\\
A_N=G^{-1}h&=\begin{pmatrix}\frac{2 \sqrt{5}}{25} & 0\\0 & 0\end{pmatrix}
\end{aligned}
$$
となる。$G$ の逆を掛けるのは $g(A_N\partial_i,\partial_j)=h_{ij}$ により、二つの共変添字の一方を出力の添字に直すためである。この行列の特性多項式は $\det(\lambda I-A_N)=\lambda^{2} - \frac{2 \sqrt{5} \lambda}{25}$ であり、主曲率は $\frac{2 \sqrt{5}}{25},\quad 0$ である。主方向を単位化するときも、座標の通常のノルムではなく $G$ のノルムを使う。

平均曲率ベクトルは $(\tr A_N/2)N=\frac{\sqrt{5}}{25}N$、Gauss 曲率は $\det A_N=0$ になる。同じ $K$ は $\det h/\det G=\det(\operatorname{Hess}_{\mathrm E}f)/W^4$ でも得られる。二つの $W$ は第二基本形式の法線への射影に、残る二つは Gram 行列式に由来する。法線を $-N$ に変えれば二つの主曲率は反転するので $K$ は変わらない。平均曲率の法成分は反転するが、掛ける法線も逆になるため平均曲率ベクトルは変わらない。

一方の変数に高さが依存しないので、その座標曲線は周囲の直線であり、法加速度は零である。
:::

## Gauss・Codazzi・Ricci の方程式

式 :eqref[eq-rg-gauss-weingarten] を二回適用すると、周囲の曲率の接成分と法成分が得られる。$\nabla_X\mathrm{II}$ には出力側の法接続を用い、二つの入力側には接接続を用いる。法接続の交換子を $R^\perp$ と書く。

:::theorem[曲率の三つの方程式]{#thm-rg-submanifold-equations}
接場 $X,Y,Z,W$ と法場 $\xi,\eta$ について
$$
\begin{aligned}
g(R(X,Y)Z,W)
={}&\bar g(\bar R(X,Y)Z,W)\\
&+\bar g(\mathrm{II}(X,W),\mathrm{II}(Y,Z))
-\bar g(\mathrm{II}(Y,W),\mathrm{II}(X,Z)),\\
(\bar R(X,Y)Z)^\perp
={}&(\nabla_X\mathrm{II})(Y,Z)-(\nabla_Y\mathrm{II})(X,Z),\\
\bar g(R^\perp(X,Y)\xi,\eta)
={}&\bar g(\bar R(X,Y)\xi,\eta)+g([A_\xi,A_\eta]X,Y).
\end{aligned}
\label{eq-rg-submanifold-curvature}
$$
上から Gauss、Codazzi、Ricci の方程式という。
:::

:::proof
$\bar\nabla_X\bar\nabla_YZ$ の接成分は $\nabla_X\nabla_YZ-A_{\mathrm{II}(Y,Z)}X$ である。$X,Y$ を交換した差から括弧の微分を引き、形作用素を第二基本形式に戻すと第 1 式を得る。法成分は $\mathrm{II}(X,\nabla_YZ)+\nabla_X^\perp\mathrm{II}(Y,Z)$ である。差を取り、$[X,Y]=\nabla_XY-\nabla_YX$ を代入すると、二入力の共変微分が揃って第 2 式となる。$\bar\nabla_X\bar\nabla_Y\xi$ の法成分は $-\mathrm{II}(X,A_\xi Y)+\nabla_X^\perp\nabla_Y^\perp\xi$ である。差と括弧の項を取った上で $\eta$ と組にすると、形作用素の二項は $-g(A_\eta X,A_\xi Y)+g(A_\eta Y,A_\xi X)=-g([A_\xi,A_\eta]X,Y)$ である。移項して最後の式を得る。
:::

:::example[同じ内在曲率を持つ異なる形]
半径 $2$ の円柱では単位主方向に形作用素の固有値 $-1/2,0$ が対応する。Gauss の式から $K=(-1/2)\cdot0=0$ である。平面でも $K=0$ だが、円柱では周方向の第二基本形式が非零である。
:::

:::example[余次元二の法曲率]
Euclid 空間で二つの法方向の形作用素が $A_\xi=\diag(1,-1)$、$A_\eta=\begin{pmatrix}0&1\\1&0\end{pmatrix}$ なら、$[A_\xi,A_\eta]=\begin{pmatrix}0&2\\-2&0\end{pmatrix}$ である。従って $\bar g(R^\perp(e_1,e_2)\xi,\eta)=-2$ である。この点では法接続を近傍全体で平坦にできない。
:::

### 逆問題と超曲面の基本定理

Gauss と Codazzi の式は、与えた計量と形作用素が Euclid 超曲面から得られるための条件でもある。

:::theorem[Euclid 超曲面の局所実現]{#thm-rg-local-hypersurface}
$(M^m,g)$ に滑らかな自己共役作用素 $A:TM\to TM$ があり、$h(X,Y)=g(AX,Y)$ とする。$g$ の曲率と $h$ が周囲平坦の場合の Gauss と Codazzi の式を満たすなら、各点の近くに第一基本形式 $g$、指定した単位法線での形作用素 $A$ を持つ $\R^{m+1}$ への等長はめ込みがある。連結な座標近傍上で実現が二つあれば、周囲の剛体運動によって互いに移る。
:::

:::proof
$E=TM\oplus\R$ に直和の計量を入れ、接続を
$$
D_X(Y,a)=(\nabla_XY-aAX,\ X(a)+h(X,Y))
$$
で定める。積を微分すると $-aAX$ と $h(X,Y)$ の交差項が相殺するため計量と両立する。曲率を $(Z,0)$ に評価すると、接成分は $R(X,Y)Z-h(Y,Z)AX+h(X,Z)AY$、法成分は $(\nabla_Xh)(Y,Z)-(\nabla_Yh)(X,Z)$ である。仮定により両者が零である。曲率は反対称作用素なので、残る法方向の作用も零となる。実際、接空間との全内積が零で、法部分は一次元上の反対称作用素だからである。

座標直方体上で平行な正規直交標構を作る。:ref[thm-rg-flat] の証明の座標軸に沿う延長は、ここでも $D$ の曲率が零なので使える。この標構で $E$ を $\R^{m+1}$ と同一視し、$\mathbb R^{m+1}$ 値 1 形式 $\alpha(X)$ を $(X,0)$ の成分とする。この形式の外微分は $D_X(Y,0)-D_Y(X,0)-([X,Y],0)$ の成分である。接部分は捩率零、法部分は $h$ の対称性で消える。よって各成分を積分して $\alpha=dF$ とできる。$\alpha$ は点ごとに等長なので $F$ は等長はめ込みであり、$(0,1)$ の成分を単位法線に取れば $D_X(0,1)=(-AX,0)$ が形作用素を与える。

逆に二つのはめ込みは、それぞれ周囲の通常の微分と $D$ を結ぶ平行自明化を与える。二つの平行標構の変換は微分が零なので連結領域で一定の直交行列である。従って微分 $dF$ はその直交行列だけ異なり、積分定数で平行移動が加わる。
:::

### 空間形の中での縮約

周囲の断面曲率が $c$、部分多様体の次元が $m>0$ とする。ここで $H=m^{-1}\sum_i\mathrm{II}(e_i,e_i)$ は第二基本形式の平均である。曲率作用素は $\bar R(X,Y)Z=c(\bar g(Y,Z)X-\bar g(X,Z)Y)$ なので、三つの入力が接方向なら法成分は零である。接方向 $X,Y$ と法方向 $\xi$ では $\bar R(X,Y)\xi=0$ である。従って Codazzi の式の左辺は零で、Ricci の式には形作用素の交換子だけが残る。

接ベクトル $X,Y$ を単位直交に選べば、Gauss の式は
$$
K_M(X,Y)=c+\inner{\mathrm{II}(X,X)}{\mathrm{II}(Y,Y)}
-\norm{\mathrm{II}(X,Y)}^2
$$
になる。正規直交基底 $e_i$ で $X=e_i$、$Y=e_j$ とし、$i,j$ の順序付きの組で和を取る。$i=j$ の補正は二項が等しいので消える。従って
$$
S_M=m(m-1)c+m^2\norm{H}^2-\norm{\mathrm{II}}^2,\qquad
\norm{\mathrm{II}}^2=\sum_{i,j}\norm{\mathrm{II}(e_i,e_j)}^2
$$
となる。特に Euclid 空間の極小部分多様体では $S_M=-\norm{\mathrm{II}}^2$ である。二次元なら $K_M=S_M/2$ なので、極小曲面で Gauss 曲率が零となる点では第二基本形式全体が零になる。

:::example[平均曲率が零でも法曲率が残る点]
Euclid 空間で $F(x,y)=(x,y,x^2-y^2,2xy)$ を考える。原点で $\mathrm{II}(e_1,e_1)=2e_3$、$\mathrm{II}(e_2,e_2)=-2e_3$、$\mathrm{II}(e_1,e_2)=2e_4$ である。従って $H=0$、$\norm{\mathrm{II}}^2=4+4+2\cdot4=16$、$K=-8$ である。また
$$
[A_{e_3},A_{e_4}]=
\begin{pmatrix}0&8\\-8&0\end{pmatrix}
$$
なので $\inner{R^\perp(e_1,e_2)e_3}{e_4}=-8$ になる。平均曲率が零であることから法接続の平坦性は従わない。
:::

### 余次元が二以上の場合の実現

抽象的な Euclid はめ込みのデータとして、計量付き束 $N\to M$、その計量接続 $D^\perp$、対称な $N$ 値双線形形式 $B$ を用意する。$g(A_\xi X,Y)=\inner{B(X,Y)}{\xi}$ と定め、$TM\oplus N$ に
$$
D_X(Y,\xi)=(\nabla_XY-A_\xi X,\ D^\perp_X\xi+B(X,Y))
$$
を入れる。接成分と法成分の交差項は互いに反対符号なので $D$ は計量と両立する。接入力 $(Z,0)$ への曲率の接成分と法成分は、Gauss と Codazzi の式の残差である。法入力 $(0,\xi)$ への法成分は $R^\perp(X,Y)\xi-B(X,A_\xi Y)+B(Y,A_\xi X)$ で、任意の $\eta$ に内積を取ると Ricci の式の残差である。曲率は反対称作用素なので、接入力への法成分が零なら法入力への接成分も零になる。従って三つの式を全て仮定すると $D$ は平坦である。

平行標構で $D$ を通常の微分に変え、$(X,0)$ を表す $\R^{m+\rank N}$ 値 1 形式を $\alpha$ とする。$d\alpha(X,Y)$ は接部分で捩率、法部分で $B(X,Y)-B(Y,X)$ に対応するので零である。閉 1 形式を各成分で積分すると $dF=\alpha$ となり、$\alpha$ の等長性から $F$ は等長はめ込みである。$N$ の平行自明化での像は法空間全体と一致するから、$B,D^\perp$ も構成したはめ込みの第二基本形式と法接続である。一意性は二つの平行標構が連結領域で一定の直交変換だけ異なることから従う。従って超曲面の場合の平坦接続による証明は同じ構造を持つが、余次元が大きいと Ricci の条件を別に確認する必要がある。

### 周囲の曲率が非零の場合の基本定理

周囲の断面曲率を定数 $c\ne0$ とする。上と同じ抽象データ $(g,N,D^\perp,B)$ が、この周囲曲率に対する Gauss、Codazzi、Ricci の方程式を満たせば、定曲率模型への等長はめ込みが局所的に存在する。周囲の等長変換を除いて一意であり、$M$ が単連結なら大域的に実現できる。ここで模型は $c>0$ の場合は半径 $1/\sqrt c$ の球面、$c<0$ の場合は半径尺度 $1/\sqrt{-c}$ の双曲面である。

実際、$E=TM\oplus N\oplus\R$ の最後の成分に双線形形式 $ab/c$ を入れ、接続を
$$
D_X(Y,\xi,a)=
\bigl(\nabla_XY-A_\xi X+aX,\ D^\perp_X\xi+B(X,Y),\
X(a)-c\,g(X,Y)\bigr)
$$
と定める。新たに加えた交差項は、内積の微分で $a\,g(X,Y)$ と $-c\,g(X,Y)(a/c)$ が相殺するため、$D$ は $E$ の双線形形式を保つ。接入力 $(Z,0,0)$ への曲率の接成分には、Euclid の場合の式に加え $-c(g(Y,Z)X-g(X,Z)Y)$ が現れる。それは最後の成分 $-c\,g(Y,Z)$ を微分する際に、接成分へ移す項である。従って接成分は Gauss の式により零となる。法成分には最後の直線からの寄与がなく、Codazzi の式で零となる。法入力への法成分も Ricci の式で零である。

さらに $s=(0,0,1)$ とおくと $D_Xs=(X,0,0)$ であり、$R^D(X,Y)s=0$ である。実際、接成分は $\nabla_XY-\nabla_YX-[X,Y]=0$、法成分は $B(X,Y)-B(Y,X)=0$、最後の成分は $-c(g(X,Y)-g(Y,X))=0$ となる。計量との両立により曲率作用素は反対称なので、残る成分もすべて零である。平行標構で $s$ の成分を $F$ とすれば、周囲空間の形式で $F$ の二乗長は $1/c$、$dF(X)$ は $(X,0,0)$ の像である。よって $F$ は該当する球面または双曲面への等長はめ込みとなる。双曲面では標構の向きを選んで所定の連結成分に置く。$D$ の接・法成分から、指定した第二基本形式と法接続が得られる。平行標構の定数の変更は模型の等長変換を与えるので、一意性も従う。単連結の場合は、[曲率の章の平坦接続](/math/riemannian-geometry/curvature/) と同じくホモトピーに沿う局所平行移動の合成で大域平行標構を取れる。$c=0$ の場合も、この標構で得た閉 1 形式を経路積分し、単連結性で経路に依存しないことを使えば大域実現が得られる。

:::example[双曲放物面の鞍点の曲率の計算]
グラフ $F(x,y)=(x,y,f(x,y))$、$f=2 x^{2} - 2 y^{2}$ を点 $\left( 0, \  0\right)$ で調べる。

$F_x=(1,0,f_x)$、$F_y=(0,1,f_y)$ だから、計量行列は $G=I+df\,df^\top$ である。点での $df$ の列は $\begin{pmatrix}0\\0\end{pmatrix}$ なので、
$$
\begin{aligned}
G&=\begin{pmatrix}1 & 0\\0 & 1\end{pmatrix}\\
W=\sqrt{1+|df|^2}&=1\\
N=\frac{(-f_x,-f_y,1)}W&=\begin{pmatrix}0\\0\\1\end{pmatrix}
\end{aligned}
$$
となる。$N$ は分子の二乗長が $W^2$ であるため単位で、$N\cdot F_x=N\cdot F_y=0$ のため法線である。座標での面積要素の倍率は $\sqrt{\det G}=W$ である。実際、$df$ に垂直な方向の二乗長は $1$、平行方向は $1+|df|^2$ となる。$df=0$ の場合は両方とも $1$ である。

周囲微分 $F_{ij}$ は鉛直成分 $f_{ij}$ だけを持つ。第二基本形式の $N$ 成分は $F_{ij}\cdot N=f_{ij}/W$ である。従って
$$
\begin{aligned}
\bigl(h_{ij}\bigr)=\frac{\operatorname{Hess}_{\mathrm E} f}W&=\begin{pmatrix}4 & 0\\0 & -4\end{pmatrix}\\
A_N=G^{-1}h&=\begin{pmatrix}4 & 0\\0 & -4\end{pmatrix}
\end{aligned}
$$
となる。$G$ の逆を掛けるのは $g(A_N\partial_i,\partial_j)=h_{ij}$ により、二つの共変添字の一方を出力の添字に直すためである。この行列の特性多項式は $\det(\lambda I-A_N)=\lambda^{2} - 16$ であり、主曲率は $4,\quad -4$ である。主方向を単位化するときも、座標の通常のノルムではなく $G$ のノルムを使う。

平均曲率ベクトルは $(\tr A_N/2)N=0$、Gauss 曲率は $\det A_N=-16$ になる。同じ $K$ は $\det h/\det G=\det(\operatorname{Hess}_{\mathrm E}f)/W^4$ でも得られる。二つの $W$ は第二基本形式の法線への射影に、残る二つは Gram 行列式に由来する。法線を $-N$ に変えれば二つの主曲率は反転するので $K$ は変わらない。平均曲率の法成分は反転するが、掛ける法線も逆になるため平均曲率ベクトルは変わらない。

平均曲率は零であるが、二つの主曲率はともに非零である。
:::

## 全測地的部分多様体

:::definition[全測地性]{#def-rg-totally-geodesic}
内在的なすべての測地線が周囲の測地線となる等長はめ込みを :term[全測地的]{en="totally geodesic" yomi="ぜんそくちてき"} という。
:::

式 :eqref[eq-rg-gauss-weingarten] より、$\mathrm{II}=0$ なら全測地的である。逆に任意の $v$ を初速度とする内在測地線を取ると、周囲加速度が零だから $\mathrm{II}(v,v)=0$ である。対称双線形性と偏極から $\mathrm{II}=0$ となる。また初期速度が接する周囲測地線は、内在測地線との一意性により局所的に像に留まる。完備性を仮定しなければ、この記述を全時刻へ延長することはできない。

:::example[反射の固定点]
単位球面の赤道は周囲の線形反射の固定点成分である。固定点において赤道に接するベクトルは反射の微分で固定される。この初期値の測地線は反射後も同じ初期値を持つので赤道に留まる。例えば赤道の長さは $2\pi$ で、周囲の大円の長さと一致する。
:::

:::example[積のスライス]
$M_1\times\{p_2\}$ は積接続が第 1 因子だけを微分するので全測地的である。$M$ の対角線 $p\mapsto(p,p)$ の誘導計量は $2g$ であり、測地線 $(\gamma,\gamma)$ が積の測地線になるため対角線も全測地的である。
:::

### 等長対合の固定点の局所構造

等長写像 $\sigma$ が $\sigma^2=\id$ を満たし $p$ を固定するとする。その微分 $L=d\sigma_p$ は直交し $L^2=I$ なので、$T_pM$ は $+1,-1$ 固有空間に直交分解する。等長写像は Levi-Civita 接続を保つため、測地線 ODE の一意性により $\sigma(\exp_p v)=\exp_p(Lv)$ である。小さい正規球で $\exp_p$ が単射なら、その中の固定点は正確に $+1$ 固有空間の小球の像である。これにより固定点集合が局所的に滑らかな部分多様体となり、接空間が $+1$ 側になる。任意のその接初速度から出る測地線は $\sigma$ で不変だから、固定点集合の各成分は全測地的である。

:::example[入れ替えの固定点]
$M\times M$ に積計量を入れ、$\sigma(x,y)=(y,x)$ を考える。固定点は対角線で、$+1$ 固有空間は $(v,v)$、$-1$ 側は $(v,-v)$ である。例えば各因子が単位球面なら、対角線の誘導計量は $2g$ だからその断面曲率は $1/2$ となる。周囲の対角方向の単位直交対 $(e_1,e_1)/\sqrt2,(e_2,e_2)/\sqrt2$ に対し、積の曲率の分子は各因子で $1/4$ ずつであり、同じ値を得る。
:::

:::example[全測地的でも周囲に近道がある]
平坦トーラス $\R^2/\Z^2$ に $t\mapsto[(2t,t)]$ を $t\in\R/\Z$ からはめ込む。これは埋め込みである。実際、二点の差 $t-s$ が同じ商の点を与えるなら第 2 座標から $t-s\in\Z$ となる。周囲では直線の射影なので全測地的で、誘導円の周長は $\sqrt5$ である。$t=0,1/2$ の内在距離は $\sqrt5/2$ だが、周囲では $[(1,1/2)]=[(0,1/2)]$ を縦に結べば距離 $1/2$ である。測地線方程式を保つことから、大域的な距離保存までは従わない。
:::

## 平均曲率と極小部分多様体

:::definition[平均曲率]{#def-rg-mean-curvature}
$m=\dim M>0$ とする。:term[平均曲率ベクトル]{en="mean curvature vector" yomi="へいきんきょくりつべくとる"} を $H=m^{-1}\sum_i\mathrm{II}(e_i,e_i)$ とする。$H=0$ のとき :term[極小部分多様体]{en="minimal submanifold" yomi="きょくしょうぶぶんたようたい"} という。
:::

:::proposition[体積の第一変分]{#prop-rg-area-first-variation}
はめ込みの変分 $F_s$ があるコンパクト集合の外で $F_0=F$ と一致するとし、$V=\partial_s F_s|_0$ とする。変分の台を含む相対コンパクトな領域を取ると、その領域の体積の微分は $-m\int_M\bar g(H,V^\perp)\mu_g$ である。台の外では密度が変化しないので、微分の値は領域の選び方によらない。従って全てのそのような変分で停留することは $H=0$ と同値である。
:::

:::proof
誘導計量行列の微分は $\dot g_{ij}=\bar g(\bar\nabla_iV,dF\partial_j)+\bar g(dF\partial_i,\bar\nabla_jV)$ である。混合微分を交換できるのは周囲の捩率が零だからである。体積密度の微分は $\frac12\tr_g(\dot g)\mu_g$ なので、正規直交標構では $\sum_i\bar g(\bar\nabla_iV,e_i)\mu_g$ となる。接成分は $\divg V^\top$、法成分は :eqref[eq-rg-gauss-weingarten] より $-m\bar g(H,V^\perp)$ である。コンパクト台の発散積分は :ref[prop-rg-green] の証明と同様に零である。逆に $H\ne0$ の近傍を選び $V=\chi H$ とすると、周囲の指数写像で変分を作れ、非負の $\chi$ を適切に選べば積分は負になる。
:::

Euclid 空間のはめ込み $F$ の成分に $\Delta_g$ を作用させると $\Delta_g F=mH$ となる。Hessian の各成分は周囲の通常の微分から $dF(\nabla_XY)$ を引いたものなので、トレースは $\mathrm{II}$ の和である。このため正の次元の閉じたコンパクト多様体は Euclid 空間に極小にはめ込めない。各成分が調和であり、:ref[prop-rg-green] により各連結成分上で一定になってしまうからである。

:::example[Clifford トーラス]
$S^3(1)$ に $F(\theta,\varphi)=2^{-1/2}(\cos\theta,\sin\theta,\cos\varphi,\sin\varphi)$ を入れる。球内の単位法線は前半と後半の符号を反転関係にした $N=2^{-1/2}(\cos\theta,\sin\theta,-\cos\varphi,-\sin\varphi)$ である。$dN$ は二つの単位接方向でそれぞれ $+1,-1$ 倍なので $A_N$ の固有値は $-1,1$ である。従って $H=0$ だが $\mathrm{II}\ne0$ である。Gauss の式から内在曲率は $1-1=0$ となる。
:::

:::example[カテノイドの極小性]
$F(u,v)=(a\cosh u\cos v,a\cosh u\sin v,au)$、$a>0$ とする。$F_u,F_v$ は直交し、両方の二乗長が $a^2\cosh^2u$ である。従って :eqref[eq-rg-laplace-conformal] の二次元の式で $\Delta_gF=(F_{uu}+F_{vv})/(a^2\cosh^2u)$ となる。最初の二成分では $u$ と $v$ の二階微分が反対で、第三成分ではともに零だから $\Delta_gF=0$ である。従って $H=0$ であり、円筒に似た形でも平均曲率が消えることが確かめられる。
:::

### グラフの極小曲面方程式

Euclid グラフの面積は $\int_U W\,dx$、$W=\sqrt{1+\sum_i f_i^2}$ である。逆計量は
$$
G^{-1}=I-\frac{df\,df^\top}{1+|df|^2}
$$
であり、両行列の積を取ると階数一の項の係数が $1-1/(1+|df|^2)-|df|^2/(1+|df|^2)=0$ となる。従って平均曲率の法成分は
$$
m\,\inner{H}{N}
=\frac1W\left(\sum_i f_{ii}-\frac{\sum_{i,j}f_i f_j f_{ij}}{1+|df|^2}\right)
=\sum_i\partial_i(f_i/W).
$$
最後の等式では $\partial_i W=W^{-1}\sum_j f_j f_{ij}$ を用いた。同じ式はグラフを $f+s\phi$ に変えた面積の微分からも得られる。$\phi$ をコンパクト台に取れば微分は $\int\sum_i f_i\partial_i\phi/W=-\int\phi\sum_i\partial_i(f_i/W)$ である。この値が任意の $\phi$ に対し零であることと $H=0$ は同値である。

:::example[一次元のグラフ]
$m=1$ なら方程式は $(f'/\sqrt{1+(f')^2})'=f''/(1+(f')^2)^{3/2}=0$ である。よって極小なグラフは区間ごとに直線である。例えば $f(x)=x^2$ では上向きの平均曲率は $2/(1+4x^2)^{3/2}$ であり、原点で $2$、$x=1$ で $2/(5\sqrt5)$ となる。分母は同じ座標方向が誘導計量で伸びることと、法線が傾くことの両方から出る。
:::

### 調和的な座標関数とコンパクト性

高さ関数の Hessian の式を各定数ベクトルに適用し、跡を取ると
$$
\Delta_g F=mH
$$
である。ここで左辺は Euclid 座標ごとに Laplacian を取ったベクトルである。従って極小性と全ての座標関数の調和性は同値である。正次元のコンパクトで境界のない $M$ が Euclid 空間へ極小にはめ込まれると仮定すると、:ref[prop-rg-green] により各座標は連結成分ごとに定数になる。すると $dF=0$ となり、はめ込みに反する。

球面 $S^{N-1}(a)$ に極小にはめ込んだ場合は、周囲を Euclid 空間に変えると余分に $-gF/a^2$ が第二基本形式に加わる。従って $\Delta_g F=-mF/a^2$ になる。逆にこの式が成り立てば、球面内の平均曲率は Euclid 平均曲率から $-F/a^2$ を引いたもので零となる。Euclid 座標関数を Laplacian の固有関数として調べることで、球面内の極小性を判定できる。

:::example[積の円から固有関数を作る]
$F(\theta,\varphi)=(\cos\theta,\sin\theta,\cos\varphi,\sin\varphi)/\sqrt2$ は単位 $S^3$ 上にある。誘導計量は $(d\theta^2+d\varphi^2)/2$ なので $\Delta_g=2(\partial_\theta^2+\partial_\varphi^2)$ である。各成分ではどちらか一方の偏微分だけが効き $\Delta_g F=-2F$ となる。ここで $m=2,a=1$ であるから、球面の中で $H=0$ である。Euclid 空間での平均曲率は $-F$ であり、同じはめ込みでも周囲の取り方で平均曲率が変わる。
:::

:::example[平均が零になる鞍形の点の曲率の計算]
グラフ $F(x,y)=(x,y,f(x,y))$、$f=x^{2} - \frac{y^{2}}{5}$ を点 $\left( 1, \  0\right)$ で調べる。

$F_x=(1,0,f_x)$、$F_y=(0,1,f_y)$ だから、計量行列は $G=I+df\,df^\top$ である。点での $df$ の列は $\begin{pmatrix}2\\0\end{pmatrix}$ なので、
$$
\begin{aligned}
G&=\begin{pmatrix}5 & 0\\0 & 1\end{pmatrix}\\
W=\sqrt{1+|df|^2}&=\sqrt{5}\\
N=\frac{(-f_x,-f_y,1)}W&=\begin{pmatrix}- \frac{2 \sqrt{5}}{5}\\0\\\frac{\sqrt{5}}{5}\end{pmatrix}
\end{aligned}
$$
となる。$N$ は分子の二乗長が $W^2$ であるため単位で、$N\cdot F_x=N\cdot F_y=0$ のため法線である。座標での面積要素の倍率は $\sqrt{\det G}=W$ である。実際、$df$ に垂直な方向の二乗長は $1$、平行方向は $1+|df|^2$ となる。$df=0$ の場合は両方とも $1$ である。

周囲微分 $F_{ij}$ は鉛直成分 $f_{ij}$ だけを持つ。第二基本形式の $N$ 成分は $F_{ij}\cdot N=f_{ij}/W$ である。従って
$$
\begin{aligned}
\bigl(h_{ij}\bigr)=\frac{\operatorname{Hess}_{\mathrm E} f}W&=\begin{pmatrix}\frac{2 \sqrt{5}}{5} & 0\\0 & - \frac{2 \sqrt{5}}{25}\end{pmatrix}\\
A_N=G^{-1}h&=\begin{pmatrix}\frac{2 \sqrt{5}}{25} & 0\\0 & - \frac{2 \sqrt{5}}{25}\end{pmatrix}
\end{aligned}
$$
となる。$G$ の逆を掛けるのは $g(A_N\partial_i,\partial_j)=h_{ij}$ により、二つの共変添字の一方を出力の添字に直すためである。この行列の特性多項式は $\det(\lambda I-A_N)=\lambda^{2} - \frac{4}{125}$ であり、主曲率は $\frac{2 \sqrt{5}}{25},\quad - \frac{2 \sqrt{5}}{25}$ である。主方向を単位化するときも、座標の通常のノルムではなく $G$ のノルムを使う。

平均曲率ベクトルは $(\tr A_N/2)N=0$、Gauss 曲率は $\det A_N=- \frac{4}{125}$ になる。同じ $K$ は $\det h/\det G=\det(\operatorname{Hess}_{\mathrm E}f)/W^4$ でも得られる。二つの $W$ は第二基本形式の法線への射影に、残る二つは Gram 行列式に由来する。法線を $-N$ に変えれば二つの主曲率は反転するので $K$ は変わらない。平均曲率の法成分は反転するが、掛ける法線も逆になるため平均曲率ベクトルは変わらない。

この点で平均曲率は零である。指定点では二階偏微分の通常の跡は $2-2/5=8/5$ であり、誘導計量で形作用素を求めてから取る跡は零である。
:::

### 極小超曲面の面積の第二変分

平均曲率が零であれば面積は停留する。そこから面積が局所最小かを判定するには、法方向への変分の二次項が必要である。以下、二側的、すなわち大域単位法線 $N$ を持つ極小超曲面を考える。コンパクト台の滑らかな関数 $u$ で、各点を初速度 $uN$ の周囲測地線に沿って動かす。

:::proposition[極小超曲面の法変分]{#prop-rg-minimal-second-variation}
周囲の Ricci 曲率を $\overline{\mathrm{Ric}}$、形作用素を $A$ とすると、上の変分での面積の第二変分は
$$
\int_M\left(\norm{\grad u}^2-
\bigl(\norm A^2+\overline{\mathrm{Ric}}(N,N)\bigr)u^2\right)\mu_g .
\label{eq-rg-minimal-index}
$$
ここで $\norm A^2$ は正規直交基底での形作用素の成分の二乗和である。この積分がすべての $u$ に対して非負のとき、超曲面は安定という。
:::

:::proof
変分を $F_s$ とし、$S=\partial_sF_s$、局所座標ベクトルの像を $X_i=\partial_iF_s$ とおく。$s=0$ の点で $X_i$ を正規直交に選べる。捩率零により $\overline D_sX_i=\overline D_iS$ である。また $\overline D_sS=0$ だから、曲率の交換子を用いて
$$
\overline D_s^2X_i=\overline D_s\overline D_iS
=\bar R(S,X_i)S.
$$
$s=0$ で $S=uN$ であるので、Weingarten の式から $\overline D_sX_i=u_iN-uAX_i$ である。従って $G_{ij}(s)=\bar g(X_i,X_j)$ の一階微分は $G'_{ij}=-2u h_{ij}$、二階微分の跡は
$$
\tr G''=2\sum_i\bigl(
\bar g(\bar R(uN,X_i)uN,X_i)+\norm{u_iN-uAX_i}^2\bigr)
=2\bigl(\norm{\grad u}^2+u^2\norm A^2-u^2\overline{\mathrm{Ric}}(N,N)\bigr).
$$
曲率項の負号は最初の添字対の交換による。法接続で $\overline D_i N$ の法成分が零なのは、法束が一次元で $N$ の長さが一定だからである。

密度の対数は $\frac12\log\det G$ である。行列式の微分と $(G^{-1})'=-G^{-1}G'G^{-1}$ から、$G(0)=I$ で密度の二階微分を密度自身で割ると
$$
\tfrac12\tr G''-\tfrac12\tr((G')^2)+\tfrac14(\tr G')^2
$$
となる。ここで $\tr G'=-2u\tr A=0$ は極小性による。また $\tr((G')^2)=4u^2\norm A^2$ である。従って $\frac12\tr G''$ が持つ $u^2\norm A^2$ から $2u^2\norm A^2$ を引き、式 :eqref[eq-rg-minimal-index] を得る。変分の台がコンパクトなので、密度の二階微分を領域全体で積分できる。周囲測地線で定めた変分は、この値を求めるための一つの変分である。停留点では、同じ初速度を持つ変分で二階微分を計算しても、加速度の違いには面積の一階微分が掛かるため値は等しい。
:::

:::example[球面の赤道の面積が減る方向]
単位 $S^{m+1}$ の赤道 $S^m$ は全測地的で $A=0$ だから極小である。一方、周囲の $\overline{\mathrm{Ric}}(N,N)=m$ である。$u=1$ を使うと第二変分は $-m\operatorname{Vol}(S^m)$ となる。実際、赤道を距離 $s$ だけ法に運ぶと、その誘導計量は $\cos^2s$ 倍になり、面積は $\operatorname{Vol}(S^m)\cos^m s$ である。初回の微分は原点で零だが、二階微分は上の負の値になる。周囲に境界を置かず全体の面積を比べると、赤道は安定でない。これに対し Euclid 空間の超平面は $A=0$、周囲 Ricci も零であり、第二変分は $\int|\grad u|^2$ だから安定である。二つの例はともに全測地的であるが、周囲曲率の項が異なる。
:::

:::example[カテノイドの広い領域で面積を減らす]
標準カテノイドを $F(t,\theta)=(\cosh t\cos\theta,\cosh t\sin\theta,t)$ と書く。計量は $\cosh^2t(dt^2+d\theta^2)$ で、Gauss 曲率は共形計量の公式から $-\cosh^{-4}t$ である。実際、共形指数は $\log\cosh t$、その二階微分は $\cosh^{-2}t$ である。極小なので主曲率の和が零であり、積がこの $K$ である。従って $\norm A^2=-2K=2\cosh^{-4}t$ である。$t$ だけの関数 $u$ に対し式 :eqref[eq-rg-minimal-index] は
$$
2\pi\int_\R\left((u')^2-2\operatorname{sech}^2t\,u^2\right)dt
$$
になる。まず $u(t)=\operatorname{sech}t$ と取ると $u''=(1-2\operatorname{sech}^2t)u$ である。境界項は指数的減衰により零なので、積分の中身を部分積分して $-\int_\R u^2dt=-2$ となる。従って値は $-4\pi$ である。

この $u$ 自身はコンパクト台を持たない。そこで $[-L,L]$ で $1$、$[-L-1,L+1]$ の外で零、導関数が $L$ によらず有界な滑らかな切断関数 $\chi_L$ を取り、$u_L=\chi_Lu$ を代入する。$(u_L')^2$ は $(|u'|+C|u|)^2$ で押さえられ、これは可積分である。ポテンシャル項も $2\operatorname{sech}^2t\,u^2$ で押さえられるので、支配収束により第二変分は $-4\pi$ に収束する。従って十分大きい有限の $L$ で負となり、コンパクトに支えられた変分だけでも不安定性が分かる。
:::

### 任意の余次元での安定性

極小はめ込みのコンパクト台を持つ法変分場を $\xi$ とする。各点を $\xi$ を初速度とする周囲の測地線に沿って動かすと、面積の第二変分は
$$
\int_M\left(
\sum_i\norm{\nabla^\perp_{e_i}\xi}^2
-\sum_i\bar g(\bar R(\xi,e_i)e_i,\xi)
-\sum_{i,j}\bar g(\mathrm{II}(e_i,e_j),\xi)^2
\right)\mu_g
$$
である。ここで $e_i$ は正規直交接標構である。超曲面の場合と同じく、座標微分したベクトル $X_i$ に対して $\overline D_s^2 X_i=\bar R(\xi,e_i)\xi$ である。一方、一階微分は合成した Weingarten の公式により $\overline D_s X_i=-A_\xi e_i+\nabla^\perp_{e_i}\xi$ である。従って $\frac12\tr G''$ には曲率項と、この直交和の二乗長が入る。$G'$ は $-2g(A_\xi e_i,e_j)$ なので、行列式の二階微分で $2\sum_i\norm{A_\xi e_i}^2$ を引く。$\tr G'=-2m\,\bar g(H,\xi)=0$ だから跡の二乗は消える。残る形作用素の項は $-\sum_i\norm{A_\xi e_i}^2$ であり、第二基本形式で書けば表示した式になる。従って大域的な単位法線が選べない超曲面でも、法束のコンパクト台を持つ切断を試験場にして安定性を定義できる。

:::example[平面を任意の法方向へ動かす]
Euclid 空間内の平面部分では曲率項と第二基本形式がともに零である。定数直交法標構で $\xi=\sum_\alpha u_\alpha N_\alpha$ と書けば、第二変分は $\sum_\alpha\int|\grad u_\alpha|^2$ である。例えば平面全体でコンパクト台を持つ非零の $\xi$ に対しては正である。積分が零なら全成分が定数であり、コンパクト台の条件により各定数は零となるからである。
:::

### 極小グラフを閉形式で比較する

極小グラフ $F(x)=(x,f(x))$ を計算するとき、第一変分の方程式は $\sum_i\partial_i(f_i/W)=0$ であった。$U\times\R$ 全体に上向き単位ベクトル場 $N(x,z)=(-df,1)/W$ を延長する。成分は $z$ に依存しないため、周囲の発散は $-\sum_i\partial_i(f_i/W)=0$ である。Euclid の体積形式 $\omega$ を $N$ に挿入した $m$ 形式 $\beta=\iota_N\omega$ は $d\beta=(\divg N)\omega=0$ を満たす。任意の向き付き直交単位の $m$ 枠に対する $\beta$ の絶対値は $1$ 以下である。実際、$m$ 枠と一緒に作る単位法線を $\nu$ とすれば、評価は選んだ向きのもとで $\inner N\nu$ である。グラフ上では $N$ が法線なので、向きを合わせれば $\beta$ は誘導体積形式と一致する。

グラフのコンパクトな部分を、境界を固定した滑らかな変形で別の向き付きはめ込みへ移す。その変形が $U\times\R$ 内に留まるなら、変形された面に沿う $\int\beta$ は変わらない。実際、変形のパラメータ区間と領域の積に $\beta$ を引き戻して Stokes を使うと、側面は境界が固定されているため階数が低く積分が零であり、内部の $d\beta$ も零である。グラフでは $\int\beta$ が面積そのもの、変形後では点ごとの評価が面積密度以下である。従って元のグラフの面積はその変形の面積以下である。ここでは停留から最小性を直ちに結論したのではなく、閉形式の積分によって面積を下から評価した。グラフの領域外に出る場合には、同じ閉形式をそこで定義できるかを確認しなければならない。

## 超曲面と主曲率

:::definition[主曲率]{#def-rg-principal-curvatures}
超曲面で単位法線 $N$ を局所的に選ぶ。自己共役な $A_N$ の固有値 $\kappa_i$ を :term[主曲率]{en="principal curvature" yomi="しゅきょくりつ"}、単位固有ベクトルを主方向という。
:::

主方向 $i\ne j$ の断面で $K=\bar K+\kappa_i\kappa_j$ であり、$H=(\sum_i\kappa_i/m)N$ である。$N$ を $-N$ にすると $\kappa_i$ は逆符号だが、積および $H$ は不変である。任意の接単位ベクトル $X=\sum_i x_ie_i$ に対して法曲率 $g(A_NX,X)=\sum_i\kappa_ix_i^2$ なので、最大と最小をとる方向は対応する固有空間である。

:::example[グラフの主曲率]
$z=x^2-2y^2$ の原点では勾配が零なので $g=I$、上向き単位法線は $(0,0,1)$ である。第二基本形式の成分は二階偏微分 $\diag(2,-4)$ である。従って主曲率は $2,-4$、Gauss 曲率は $-8$、平均曲率ベクトルは $-N$ である。
:::

### グラフの極小方程式

$f:\Omega\subset\R^m\to\R$ のグラフで $G=I+\nabla f\,(\nabla f)^\top$ だから面積密度は $W=\sqrt{1+|\nabla f|^2}$ である。$f+s\psi$ についてこの積分を微分すると $\int\nabla f\cdot\nabla\psi/W$ となり、部分積分で $-\int\psi\,\divg_{\mathrm E}(\nabla f/W)$ となる。垂直変分の法成分は $\psi/W$ で、面積密度の $W$ と相殺する。従って $m(H\cdot N)=\divg_{\mathrm E}(\nabla f/W)$ であり、この発散が零なら極小である。例えば $f=ax+by+c$ は勾配が定数なので右辺が零である。

:::example[法曲率の範囲と方向]
原点で $z=x^2-2y^2$ の上向き法線を用いる。$45$ 度の方向 $X=(e_1+e_2)/\sqrt2$ の法曲率は $2/2-4/2=-1$ である。法曲率零の方向は $2\cos^2\theta-4\sin^2\theta=0$ から $\tan^2\theta=1/2$ を満たす。従って最小と最大の間で零になる方向があるが、固有値そのものには零がない。
:::

### 平行超曲面の形作用素

Euclid 空間の向き付けた超曲面で $F_s=F+sN$ を考える。$dF_s(X)=dF((I-sA_N)X)$ である。従って $1-s\kappa_i\ne0$ が全ての主曲率で成立する間ははめ込みであり、$N$ は同じ単位法線として使える。Weingarten の式を $F_s$ で書くと
$$
-dF_s(A_sX)=D_XN=-dF(A_NX)
$$
だから $A_s=(I-sA_N)^{-1}A_N$ となる。$X$ を元の主方向に取れば、新しい主曲率は $\kappa_i/(1-s\kappa_i)$ である。また誘導体積の倍率は $\prod_i|1-s\kappa_i|$ である。実際、$X$ を元の正規直交主標構に選ぶと、各ベクトルの像は元の方向に $1-s\kappa_i$ を掛けたものになる。

:::example[円柱を法方向へ移す]
半径 $2$ の円柱の外向き法線について、主曲率は $-1/2,0$ である。距離 $s=1$ だけ外へ移ると、周方向の固有値は $(-1/2)/(1+1/2)=-1/3$、軸方向は $0$ のままである。面積倍率は $3/2$ となり、円周だけが半径に比例して伸びることと一致する。$s=-2$ では周方向の微分が零になり、像は軸へ潰れる。この点では形作用素の式の分母も零になる。
:::

### 臍点と Codazzi の条件

超曲面で $A=\lambda I$ となる点を :term[臍点]{en="umbilic point" yomi="せいてん"} という。Euclid 空間の連結超曲面で $m\ge2$ とし、すべての点が臍点とする。$\lambda=(\tr A)/m$ なので滑らかである。Codazzi の式は
$$
(X\lambda)Y=(Y\lambda)X
$$
となる。任意の単位 $X$ に直交する単位 $Y$ を選び、両辺を $Y$ と組にすると $X\lambda=0$ である。従って $\lambda$ は一定である。$\lambda=0$ なら $D_X N=0$ で $N$ は一定、$D_X\inner{F}{N}=0$ だから像は超平面に含まれる。$\lambda\ne0$ なら $D_X(F+N/\lambda)=0$ なので中心 $F+N/\lambda$ は一定であり、像は半径 $1/|\lambda|$ の球面に含まれる。いずれもはめ込みの次元が像を含む模型の次元と同じなので、局所的にその開部分となる。単位法線が局所的にしか選べないときもこの結論は局所的に成り立つ。一次元では全点が自動的に臍点で、Codazzi に使う直交二方向がないため同じ結論にならない。

:::example[混合項のある鞍面の曲率の計算]
グラフ $F(x,y)=(x,y,f(x,y))$、$f=2 x y$ を点 $\left( 0, \  0\right)$ で調べる。

$F_x=(1,0,f_x)$、$F_y=(0,1,f_y)$ だから、計量行列は $G=I+df\,df^\top$ である。点での $df$ の列は $\begin{pmatrix}0\\0\end{pmatrix}$ なので、
$$
\begin{aligned}
G&=\begin{pmatrix}1 & 0\\0 & 1\end{pmatrix}\\
W=\sqrt{1+|df|^2}&=1\\
N=\frac{(-f_x,-f_y,1)}W&=\begin{pmatrix}0\\0\\1\end{pmatrix}
\end{aligned}
$$
となる。$N$ は分子の二乗長が $W^2$ であるため単位で、$N\cdot F_x=N\cdot F_y=0$ のため法線である。座標での面積要素の倍率は $\sqrt{\det G}=W$ である。実際、$df$ に垂直な方向の二乗長は $1$、平行方向は $1+|df|^2$ となる。$df=0$ の場合は両方とも $1$ である。

周囲微分 $F_{ij}$ は鉛直成分 $f_{ij}$ だけを持つ。第二基本形式の $N$ 成分は $F_{ij}\cdot N=f_{ij}/W$ である。従って
$$
\begin{aligned}
\bigl(h_{ij}\bigr)=\frac{\operatorname{Hess}_{\mathrm E} f}W&=\begin{pmatrix}0 & 2\\2 & 0\end{pmatrix}\\
A_N=G^{-1}h&=\begin{pmatrix}0 & 2\\2 & 0\end{pmatrix}
\end{aligned}
$$
となる。$G$ の逆を掛けるのは $g(A_N\partial_i,\partial_j)=h_{ij}$ により、二つの共変添字の一方を出力の添字に直すためである。この行列の特性多項式は $\det(\lambda I-A_N)=\lambda^{2} - 4$ であり、主曲率は $-2,\quad 2$ である。主方向を単位化するときも、座標の通常のノルムではなく $G$ のノルムを使う。

平均曲率ベクトルは $(\tr A_N/2)N=0$、Gauss 曲率は $\det A_N=-4$ になる。同じ $K$ は $\det h/\det G=\det(\operatorname{Hess}_{\mathrm E}f)/W^4$ でも得られる。二つの $W$ は第二基本形式の法線への射影に、残る二つは Gram 行列式に由来する。法線を $-N$ に変えれば二つの主曲率は反転するので $K$ は変わらない。平均曲率の法成分は反転するが、掛ける法線も逆になるため平均曲率ベクトルは変わらない。

両座標軸での法曲率は零だが、主方向は座標軸を $\pi/4$ だけ回した方向である。
:::

:::example[だ円放物面の非水平な点の曲率の計算]
グラフ $F(x,y)=(x,y,f(x,y))$、$f=x^{2} + 2 y^{2}$ を点 $\left( 1, \  0\right)$ で調べる。

$F_x=(1,0,f_x)$、$F_y=(0,1,f_y)$ だから、計量行列は $G=I+df\,df^\top$ である。点での $df$ の列は $\begin{pmatrix}2\\0\end{pmatrix}$ なので、
$$
\begin{aligned}
G&=\begin{pmatrix}5 & 0\\0 & 1\end{pmatrix}\\
W=\sqrt{1+|df|^2}&=\sqrt{5}\\
N=\frac{(-f_x,-f_y,1)}W&=\begin{pmatrix}- \frac{2 \sqrt{5}}{5}\\0\\\frac{\sqrt{5}}{5}\end{pmatrix}
\end{aligned}
$$
となる。$N$ は分子の二乗長が $W^2$ であるため単位で、$N\cdot F_x=N\cdot F_y=0$ のため法線である。座標での面積要素の倍率は $\sqrt{\det G}=W$ である。実際、$df$ に垂直な方向の二乗長は $1$、平行方向は $1+|df|^2$ となる。$df=0$ の場合は両方とも $1$ である。

周囲微分 $F_{ij}$ は鉛直成分 $f_{ij}$ だけを持つ。第二基本形式の $N$ 成分は $F_{ij}\cdot N=f_{ij}/W$ である。従って
$$
\begin{aligned}
\bigl(h_{ij}\bigr)=\frac{\operatorname{Hess}_{\mathrm E} f}W&=\begin{pmatrix}\frac{2 \sqrt{5}}{5} & 0\\0 & \frac{4 \sqrt{5}}{5}\end{pmatrix}\\
A_N=G^{-1}h&=\begin{pmatrix}\frac{2 \sqrt{5}}{25} & 0\\0 & \frac{4 \sqrt{5}}{5}\end{pmatrix}
\end{aligned}
$$
となる。$G$ の逆を掛けるのは $g(A_N\partial_i,\partial_j)=h_{ij}$ により、二つの共変添字の一方を出力の添字に直すためである。この行列の特性多項式は $\det(\lambda I-A_N)=\lambda^{2} - \frac{22 \sqrt{5} \lambda}{25} + \frac{8}{25}$ であり、主曲率は $\frac{2 \sqrt{5}}{25},\quad \frac{4 \sqrt{5}}{5}$ である。主方向を単位化するときも、座標の通常のノルムではなく $G$ のノルムを使う。

平均曲率ベクトルは $(\tr A_N/2)N=\frac{11 \sqrt{5}}{25}N$、Gauss 曲率は $\det A_N=\frac{8}{25}$ になる。同じ $K$ は $\det h/\det G=\det(\operatorname{Hess}_{\mathrm E}f)/W^4$ でも得られる。二つの $W$ は第二基本形式の法線への射影に、残る二つは Gram 行列式に由来する。法線を $-N$ に変えれば二つの主曲率は反転するので $K$ は変わらない。平均曲率の法成分は反転するが、掛ける法線も逆になるため平均曲率ベクトルは変わらない。

座標 Hessian の比と主曲率の比は異なる。最初の接方向の二乗長が大きいためである。
:::

### 法指数写像と平行な面の正則性

埋め込まれたコンパクトな二側的 Euclid 超曲面 $F:M\to\R^{m+1}$ に対し、$\Phi(p,s)=F(p)+sN(p)$ とする。$s=0$ では微分は $dF(v)+aN$ で接空間と法空間の直交和だから同型である。従って各 $p$ の近くで $\Phi$ は微分同相となる。コンパクト性を使って一様な幅を取れることを示す。もし任意に小さい幅で異なる二入力が同じ像を持つなら、$(p_j,s_j)$、$(q_j,t_j)$ を $s_j,t_j\to0$ で取り、底点の収束部分列を選ぶ。像の一致から極限で $F(p)=F(q)$ なので、埋め込みの単射性により $p=q$ である。すると二入力は最後には同じ局所逆の近傍にあり、単射性に反する。微分の非特異性も同様に、零切断のコンパクト性と連続性で一様な幅に確保できる。

この幅を縮めれば、像内の点 $\Phi(p,s)$ の超曲面までの距離は $|s|$ である。実際、コンパクト性から最も近い $q$ があり、$|F(q)-\Phi(p,s)|^2$ を $q$ で微分すると差は法線に平行である。候補の距離は $|s|$ 以下なので、この点は同じ幅内の別の法座標表示 $\Phi(q,t)$ を与える。単射性により $p=q,s=t$ である。従って符号付き距離関数はこの管状近傍で $s$ に等しい。微分 $d\Phi(v,a)=dF((I-sA)v)+aN$ により $\grad s=N$ である。平行面の接方向では $\operatorname{Hess}s(X,Y)=g(\bar\nabla_XN,Y)=-g(A_sX,Y)$ である。従って符号付き距離の Hessian を使えば、平行超曲面の曲率を座標に依存せず読める。

:::example[球面の内外で符号付き距離を測る]
半径 $a$ の球面の外向き法線に対し $s(x)=|x|-a$ である。球面近くの点で $\grad s=x/|x|$ を微分すると、動径方向では零、球面に接する方向では $1/|x|$ 倍である。従って平行面での形作用素は $-1/(a+s)$ 倍である。前に求めた $\kappa/(1-s\kappa)$ に $\kappa=-1/a$ を入れても同じ値になる。例えば $a=3$、$s=-1$ では主曲率は $-1/2$ である。$s=-a$ になると全法線が原点で集まり、最も近い点の一意性と滑らかな距離座標が同時に失われる。
:::

### 球面内の二球の積

$m=p+q$ とし、$a,b>0,\ a^2+b^2=1$ のもとで $S^p(a)\times S^q(b)$ を $S^{m+1}(1)$ にはめ込む。点を $(au,bv)$、$|u|=|v|=1$ と書くと、球面内の単位法線は $N=(bu,-av)$ である。実際、位置との内積は $ab-ab=0$ で球面に接し、両因子の接ベクトルと直交し、二乗長は $a^2+b^2=1$ である。第 1 因子の接ベクトルを $(X,0)$ とすれば $D_XN=(bX/a,0)$、第 2 因子では $D_YN=(0,-aY/b)$ である。ともに接方向なので、球面の接続でもこの微分は同じである。従って主曲率は
$$
\kappa_1=-b/a\quad(p\text{ 重}),\qquad
\kappa_2=a/b\quad(q\text{ 重})
$$
である。平均曲率が零である条件は $-pb/a+qa/b=0$ であり、$a^2=p/(p+q)$、$b^2=q/(p+q)$ が必要十分である。Gauss の式で混合平面は $1+\kappa_1\kappa_2=0$、次元が二以上の因子について、第 1 因子内の平面では $1+b^2/a^2=1/a^2$、第 2 因子内の平面では $1+a^2/b^2=1/b^2$ である。これは積計量の内在曲率と一致する。

:::example[異なる次元の因子のバランス]
$p=1,q=3$ とすると極小条件は $a=1/2,\ b=\sqrt3/2$ である。主曲率は $-\sqrt3$ が一つ、$1/\sqrt3$ が三つで、和は零になる。平方和は $3+3(1/3)=4$ である。単位五次元球面で $\overline{\mathrm{Ric}}(N,N)=4$ だから、式 :eqref[eq-rg-minimal-index] の定数法変分では面積の二階微分は $-8\operatorname{Vol}(M)$ となる。第 1 因子は一次元なので純粋な第 1 因子断面は存在しない。一方、第 2 因子の断面は $4/3$、混合断面は零である。平均曲率の跡で方向の重複度を落とすと、等しい半径を誤って選んでしまう。
:::

## 演習問題

::::exercise{level=basic}
半径 $3$ の円柱で外向き法線を取り、二つの主曲率を求めよ。

:::solution
周角の変化に対して $dN=dF/3$、軸方向では $dN=0$ である。$A_N=-dN$ を接空間と同一視すると $-1/3,0$ である。
:::
::::

::::exercise{level=basic}
半径 $2$ の球面の平均曲率ベクトルを求めよ。

:::solution
外向き法線で各主曲率は $-1/2$ なので平均は $-1/2$ であり、$H=-N/2$ となる。
:::
::::

::::exercise{level=basic}
Euclid 空間の超曲面で直交主方向の主曲率が $2,3$ なら断面曲率はいくつか。

:::solution
:eqref[eq-rg-submanifold-curvature] で周囲の曲率と交差第二基本形式がともに零だから、$K=2\cdot3=6$ である。
:::
::::

::::exercise{level=basic}
グラフ $z=3x^2+y^2$ の原点の Gauss 曲率を求めよ。

:::solution
原点で誘導計量は単位行列で法線は垂直なので、形作用素は Hessian $\diag(6,2)$ に等しい。従って $K=12$ である。
:::
::::

::::exercise{level=advanced}
凸な上半球を $z=\sqrt{a^2-x^2-y^2}$ と表示する。平均曲率のグラフ方程式を直接確認せよ。

:::solution
$f_i=-x_i/f$ なので $W=\sqrt{1+(x^2+y^2)/f^2}=a/f$ である。従って $\nabla f/W=-(x,y)/a$ で、発散は $-2/a$ である。平均曲率はその半分の $-1/a$ となる。
:::
::::

::::exercise{level=advanced}
Euclid 超曲面が連結で $A_N=\lambda\,\id$、次元 $m\ge2$ とする。局所的に $\lambda$ が定数であることを示せ。

:::solution
単位法線の法接続は零であり、Codazzi の式から $X(\lambda)Y-Y(\lambda)X=0$ となる。任意の単位 $X$ に直交する単位 $Y$ を選び $Y$ と内積を取れば $X(\lambda)=0$ である。すべての $X$ で成立するので $\lambda$ は局所定数である。
:::
::::

::::exercise{level=advanced}
前問で $\lambda$ が非零定数なら、像がある球面に含まれることを示せ。

:::solution
Weingarten の式で $dN(X)=-\lambda dF(X)$ となる。従って $F+N/\lambda$ の微分は零で、連結領域上では定ベクトル $c$ に等しい。$\norm{F-c}=1/|\lambda|$ なので球面の方程式を満たす。
:::
::::

::::exercise{level=advanced}
$S^3(1)$ の $S^1(a)\times S^1(b)$、$a^2+b^2=1$、$a,b>0$ が極小となる条件を求めよ。

:::solution
$F=(au,bv)$、$N=(bu,-av)$ と取ると、$dN$ は二因子上で $b/a,-a/b$ 倍である。従って $\tr A_N=-b/a+a/b=(a^2-b^2)/(ab)$ である。これが零となるのは $a=b=1/\sqrt2$ のときに限る。
:::
::::

::::exercise{level=advanced}
$S^{p+q+1}(1)$ の中の二球の積を $F(t,u,v)=(\cos t\,u,\sin t\,v)$、$0<t<\pi/2$ と表す。ここで $u\in S^p(1)$、$v\in S^q(1)$、$p,q\ge1$ である。$t$ を固定した超曲面の面積を求め、面積が停留する $t$ と、その点での面積の二階微分を求めよ。式 :eqref[eq-rg-minimal-index] の定数変分とも比較せよ。
:::solution
因子の計量はそれぞれ $\cos^2t$ 倍と $\sin^2t$ 倍なので、面積は
$$
V(t)=V_pV_q\,\cos^p t\,\sin^q t
$$
である。ここで $V_j$ は単位 $j$ 次元球面の体積である。従って $V'/V=-p\tan t+q\cot t$ となる。停留条件は $q\cos^2t=p\sin^2t$ であり、$\cos^2t=p/(p+q)$ である。区間内にその値をとる $t$ はただ一つある。この点で
$$
V''/V=(\log V)''=-p\sec^2t-q\csc^2t=-2(p+q).
$$
一方、$\partial_tF$ は単位法線で、$t$ を動かす曲線は周囲球面の大円である。主曲率の二乗和は $p\tan^2t+q\cot^2t$ であり、停留点では $q+p$ になる。周囲の法方向 Ricci 値も $p+q$ なので、式 :eqref[eq-rg-minimal-index] に $u=1$ を入れると $-2(p+q)V(t)$ となり直接の微分結果と一致する。従ってこの停留は、半径だけを変える変分の族の中ですでに局所最大である。
:::
::::

## 確認テスト

::::quiz
任意の点で第二基本形式が零の場合に成り立つものを選べ。

- [x] はめ込みは全測地的である。
- [x] はめ込みは極小である。
- [ ] 周囲の曲率は必ず零である。

:::explanation
接測地線の周囲加速度は第二基本形式なので零である。またトレースも零である。球面の赤道はこの例だが周囲の球面の曲率は正である。
:::
::::

::::quiz{answer="-1" tol="1e-6"}
外向き法線に関する単位球面の任意の主曲率はいくつか。

:::explanation
$N=F$、$dN=dF$ より $A_N=-\id$ であり、固有値は $-1$ である。
:::
::::

::::quiz{answer="0" tol="1e-6"}
半径 $2$ の Euclid 円柱の Gauss 曲率はいくつか。

:::explanation
主曲率は $-1/2$ と $0$ なので、Gauss の式で積は $0$ である。
:::
::::

::::quiz
正しい体積変分の判定を選べ。

- [x] すべてのコンパクト台の変分で体積が停留するなら $H=0$ である。
- [ ] $H=0$ は全測地性と同値である。

:::explanation
体積の第一変分に $V=\chi H$ を用いると最初が従う。Clifford トーラスでは主曲率が相殺するが個々は非零なので後者は誤りである。
:::
::::

::::quiz{answer="12" tol="1e-6"}
Euclid 空間のグラフ $z=3x^2+y^2$ の原点で、誘導計量の Gauss 曲率を求めよ。
:::explanation
原点では $df=0$ なので誘導計量の行列は $I$ であり、上向き単位法線は $(0,0,1)$ である。従って形作用素の行列は $f$ の Hessian $\diag(6,2)$ になる。周囲の断面曲率は零だから、Gauss の式で主曲率の積 $6\cdot2=12$ を得る。法線を下向きに選ぶと両主曲率が反転するが、積は同じである。
:::
::::

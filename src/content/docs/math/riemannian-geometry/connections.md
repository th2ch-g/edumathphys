---
title: 'アフィン接続と Levi-Civita 接続'
description: 'ベクトルの微分を接続で定め、計量に適合する捩率零の接続の一意性と計算法を証明する。'
notation:
  - symbol: '\nabla_XY'
    meaning: '方向 X に沿う Y の共変微分'
  - symbol: 'D_t'
    meaning: '曲線に沿う共変微分'
  - symbol: 'T'
    meaning: '捩率'
  - symbol: '\Gamma^k_{ij}'
    meaning: '座標基底に関する接続係数'
---

## この章の目標

- 接続の Leibniz 則から座標表示を得る。
- 平行移動を線形常微分方程式として構成する。
- Koszul の公式から Levi-Civita 接続を構成する。
- テンソルの共変微分を添字と標構で計算する。

## アフィン接続と共変微分

別の点に属する接ベクトルは異なるベクトル空間の元である。したがって、差分を取るには点と点の間でベクトルを比較する規則が要る。座標を固定すれば成分の微分ができるが、標構そのものが変わればその微分も加わる。

:::definition[アフィン接続]{#def-rg-affine-connection}
ベクトル場の組にベクトル場を対応させる実双線形写像 $\nabla$ が $\nabla_{fX}Y=f\nabla_XY$ と $\nabla_X(fY)=X(f)Y+f\nabla_XY$ を満たすとき、:term[アフィン接続]{en="affine connection" yomi="あふぃんせつぞく"} という。$\nabla_XY$ を :term[共変微分]{en="covariant derivative" yomi="きょうへんびぶん"} という。
:::

この規則は局所的である。例えば $Y$ が $p$ の近傍で零ならば、そこで $1$ となるコンパクト台の関数 $\chi$ を選んで $\chi Y=0$ と書ける。Leibniz 則を $p$ で評価すると $\nabla_XY(p)=0$ である。第 1 引数にも同じ局所化を使える。よって局所標構 $e_i$ に対して、$\nabla_{e_i}e_j=\Gamma^k_{ij}e_k$ とおけば $\nabla_XY=(X(Y^k)+X^iY^j\Gamma^k_{ij})e_k$ である。逆に係数を任意に指定したこの式は局所接続を与える。1 の分割で局所接続を加重平均すると Leibniz 則の $X(f)Y$ の係数は重みの和 $1$ となるから、大域接続も存在する。

### 接続の差

$\nabla-\nabla'$ は両引数で $C^\infty(M)$ 線形になる。第 2 引数での微分項が相殺するからである。そのため差は $(1,2)$ テンソルであり、逆に任意のその型のテンソル $A$ を $\nabla$ に加えれば接続である。

:::example[移動標構]
平面の通常の微分 $D$ について、$e_1=(\cos x,\sin x)$、$e_2=(-\sin x,\cos x)$ とおく。$D_{\partial_x}e_1=e_2$、$D_{\partial_x}e_2=-e_1$ である。たとえ成分が定数でも、$e_1$ 自体は平行でない。
:::

:::example[接続の自由度]
$\R$ で $\nabla_{\partial_x}\partial_x=2\partial_x$ とする。$\nabla_{x\partial_x}(x^2\partial_x)=(2x^2+2x^3)\partial_x$ である。通常の微分との差は最後の $2x^3\partial_x$ だけである。
:::

### 局所係数が定める接続

接続の第二引数で係数の微分が生じるのは、その引数を曲線に沿って変化させるためである。一方、第一引数はどの方向へ微分するかを選ぶので、点での値だけで決まる。この区別を用いると、接続を局所係数で作れる。

:::proposition[局所性と任意の係数からの構成]{#prop-rg-connection-local}
局所標構 $(E_i)$ を固定する。任意の滑らかな係数 $A^k_{ij}$ は、その座標域で $\nabla_{E_i}E_j=A^k_{ij}E_k$ を満たす接続を一意に定める。また二つの接続の差は $(1,2)$ テンソルであり、接続に任意の $(1,2)$ テンソルを加えれば別の接続となる。
:::

:::proof
$X=X^iE_i$、$Y=Y^jE_j$ に公理を適用すると
$$
\nabla_XY=\sum_k\left(X(Y^k)+X^iY^j A^k_{ij}\right)E_k
$$
が強制される。右辺を定義とすると、$X$ に関する関数線形性は $X(Y^k)$ が点での $X$ の値に線形であることから従い、$Y$ に関する Leibniz 則は $X(fY^k)=(Xf)Y^k+fX(Y^k)$ から従う。もし最初に接続を大域的な場に対して定義した場合も、局所場に切断関数を掛けて扱える。実際、$Y$ が $p$ の近傍で零なら $p$ の近くで $1$ になる関数 $\chi$ で $\chi Y=0$ とできる。Leibniz 則により $p$ で $\nabla_XY=0$ である。$X$ の方は関数線形性で同じ結論になる。従って局所場の延長方法は値に影響しない。

差 $S(X,Y)=\widehat\nabla_XY-\nabla_XY$ を取ると、第一引数はもともと関数線形で、第二引数でも $X(f)Y$ が消えて関数線形になる。局所係数の差は滑らかなので $S$ はテンソルである。逆に $\nabla+S$ は同じ Leibniz 則を満たす。
:::

:::example[非座標標構の接続行列]
Euclid 平面で $E_1=\cos x\,\partial_x+\sin x\,\partial_y$、$E_2=-\sin x\,\partial_x+\cos x\,\partial_y$ とおく。標準接続で $\nabla_XE_1=dx(X)E_2$、$\nabla_XE_2=-dx(X)E_1$ である。例えば $x=\pi/3$ では $\nabla_{E_1}E_1=\frac12E_2$、$\nabla_{E_2}E_1=-\frac{\sqrt3}{2}E_2$ である。同じ点で直交標構を一つ固定しても、どのように周囲まで延長したかにより微分の値は変わる。
:::

### 接続を貼り合わせる

もし計量を先に選ばずに接続を作るなら、座標域ごとに座標標構を平行と宣言する。従属する局所有限な $1$ の分割 $\rho_\alpha$ を選び $\nabla_XY=\sum_\alpha\rho_\alpha\nabla_X^\alpha Y$ とする。各項は $\rho_\alpha$ の台の外で零に延長できる。第二引数に $fY$ を入れたとき、微分項は $\sum_\alpha\rho_\alpha(Xf)Y=(Xf)Y$ となるので接続が得られる。重みの和が $1$ でなければ、この係数がずれてしまう。さらに局所接続が全て捩率零だった場合、差 $\nabla_XY-\nabla_YX$ は $\sum\rho_\alpha[X,Y]=[X,Y]$ となり、貼り合わせ後も捩率は零である。

## 曲線に沿う共変微分と平行移動

曲線 $\gamma(t)$ に沿う場 $V$ は各 $t$ で $T_{\gamma(t)}M$ に値を持つ。同じ点を複数回通る曲線では、各時刻で違う値を取ってよい。引き戻し束 $\gamma^*TM$ の切断として扱うとこの自由度が明確になる。

:::proposition[曲線に沿う微分]{#prop-rg-curve-derivative}
曲線に沿う微分 $D_t$ は $D_t(fV)=f'V+fD_tV$ と $D_t(Y\circ\gamma)=\nabla_{\dot\gamma}Y$ を満たす唯一の実線形演算であり、座標で
$$
(D_tV)^k=(V^k)'+\Gamma^k_{ij}(\gamma)\dot\gamma^iV^j
\label{eq-rg-curve-derivative}
$$
となる。
:::

:::proof
局所基底を用いて $V=V^ie_i\circ\gamma$ と分解し、両条件を適用すると成分式が強制される。逆に別の基底 $\widetilde e_a=S^i{}_ae_i$ を使うと、$\nabla_{\dot\gamma}\widetilde e_a=(S^i{}_a\circ\gamma)'e_i+S^i{}_a\nabla_{\dot\gamma}e_i$ である。この式と $(S^i{}_a\widetilde V^a)'$ を比較すれば、成分式は基底によらない。よって時刻区間ごとに作った微分は一致して貼り合わさる。
:::

:::definition[平行移動]{#def-rg-parallel-transport}
$D_tV=0$ を満たす場を :term[平行]{en="parallel" yomi="へいこう"} という。初期値 $V(a)=v$ に対して $P_{a,t}v=V(t)$ とおく。この線形写像を :term[平行移動]{en="parallel transport" yomi="へいこういどう"} という。
:::

線形 ODE の存在と一意性を標構ごとに用いる。コンパクトなパラメータ区間は有限個の標構域に分けられ、各区間の係数は有界なので線形 ODE は端まで解ける。解の線形性と一意性より $P_{a,t}$ は線形であり、逆は $P_{t,a}$ である。また $P_{b,c}P_{a,b}=P_{a,c}$ である。再パラメータ化では $D_sV=(dt/ds)D_tV$ となるから、向きの反転で輸送が逆転する。

:::example[指数的な成分変化]
$\R$ の $\Gamma^1_{11}=2$ で $\gamma(t)=t$ とすると $(V^1)'+2V^1=0$ であり、初期値 $3$ に対して $V^1(t)=3e^{-2t}$ である。$t=\log 2$ では $3/4$ になる。
:::

:::example[平面の極座標]
$e_r=(\cos\theta,\sin\theta)$、$e_\theta=(-\sin\theta,\cos\theta)$ とする。円 $\theta=t$ に沿う平行場は Cartesian 成分が一定である。初期値 $(1,0)$ の場は極標構で $\cos t\, e_r-\sin t\, e_\theta$ である。成分が変わってもベクトルは平行である。
:::

### 基本行列と合成

曲線に沿う標構で $D_tV$ の成分を $\dot v+A(t)v$ と書く。可逆行列 $P(t,s)$ を、各列が $\dot v=-A(t)v$ を満たし $P(s,s)=I$ となる解として作る。線形 ODE の一意性により
$$
P(t,r)P(r,s)=P(t,s),\qquad P(s,t)=P(t,s)^{-1}.
$$
実際、最初の等式の両辺を $t$ の関数とみると同じ ODE を満たし $t=r$ で一致する。右辺に $t=s$ を入れれば第二式も得る。複数の標構をまたぐ曲線では、各区間の終端で標構変換を挿入し、これらの線形写像を順に合成する。途中で分割を変更しても最初の式により同じ輸送となる。

非斉次方程式 $\dot v+Av=f$ では $v(t)=P(t,s)u(t)$ とおくと $P(t,s)\dot u=f$ なので
$$
v(t)=P(t,s)v(s)+\int_s^tP(t,r)f(r)\,dr
$$
となる。これを曲線に沿う微分の積分表示として使える。標構を明示しない場合、積分の各 $f(r)$ を時刻 $t$ の同じ接空間まで運んでから足している。

:::example[順序の効く平行移動]
$\R^2$ で標準標構に関する接続行列を $A\,dx+B\,dy$ とし、
$$
A=\begin{pmatrix}0&1\\0&0\end{pmatrix},\qquad
B=\begin{pmatrix}0&0\\1&0\end{pmatrix}
$$
を取る。$A^2=B^2=0$ なので、それぞれ単位長の $x,y$ 方向の輸送は $I-A,I-B$ である。$(0,0)$ から $(1,1)$ まで先に横へ進めば $(I-B)(I-A)=\begin{pmatrix}1&-1\\-1&2\end{pmatrix}$ となり、先に縦へ進めば $(I-A)(I-B)=\begin{pmatrix}2&-1\\-1&1\end{pmatrix}$ となる。両者は出発点と到着点が同じである。この差を微小な長方形で測ると [曲率](/math/riemannian-geometry/curvature/) の曲率作用素に結び付く。
:::

:::example[動く標構でも一定なベクトル]
平面の極座標で正規直交標構を $E_r=\partial_r$、$E_\theta=r^{-1}\partial_\theta$ とする。半径を固定した曲線 $\theta=t$ では $D_tE_r=E_\theta$、$D_tE_\theta=-E_r$ となる。従って初期値 $E_r(0)$ の平行場は $V(t)=\cos t\,E_r(t)-\sin t\,E_\theta(t)$ である。実際、その二成分を微分した $( -\sin t,-\cos t)$ は標構の微分から出る成分 $(\sin t,\cos t)$ と相殺する。$t=\pi/2$ では $V=-E_\theta$ となるが、Euclid 座標では始終 $(1,0)$ である。
:::

## 捩率テンソル

:::definition[捩率]{#def-rg-torsion}
$T(X,Y)=\nabla_XY-\nabla_YX-[X,Y]$ を :term[捩率]{en="torsion" yomi="ねじれりつ"} と呼ぶ。
:::

$[fX,Y]=f[X,Y]-Y(f)X$ なので $T(fX,Y)=fT(X,Y)$ である。交代性から第 2 引数も関数線形となる。座標基底は可換なので $T^k{}_{ij}=\Gamma^k_{ij}-\Gamma^k_{ji}$ である。一般標構では $[e_i,e_j]=c^k{}_{ij}e_k$ を引き、$T^k{}_{ij}=\Gamma^k_{ij}-\Gamma^k_{ji}-c^k{}_{ij}$ となる。この補正がないと標構の非可換性と接続の捩率を混同する。

### 二変数の変分

$\alpha(s,t)$ に沿って $S=\partial_s\alpha$、$V=\partial_t\alpha$ とする。局所座標で $D_sV-D_tS$ の通常の微分部分は混合偏微分が交換するため零になる。残る係数は $(\Gamma^k_{ij}-\Gamma^k_{ji})S^iV^j$ である。従って $D_sV-D_tS=T(S,V)$ であり、捩率が零なら変分の速度と変分場を交換できる。

:::example[非可換標構での捩率]
上半平面で $e_1=y\partial_x,\ e_2=y\partial_y$ とすると $[e_1,e_2]=-e_1$ である。この標構を平行と定めた接続では $T(e_1,e_2)=e_1$ になる。一方、通常の座標微分の捩率は零であり、同じ標構では $D_{e_2}e_1=e_1$ が Lie 括弧を相殺する。
:::

任意の接続に $-\frac12T$ を足すと捩率が零になる。実際、新しい捩率は $T(X,Y)-T(X,Y)/2+T(Y,X)/2=0$ である。ただし、この操作は次節の計量との両立を一般には保たない。

:::example[接続の反対称部分を調整する]
平面で $\nabla_{\partial_x}\partial_y=2\partial_x$ とし、ほかの座標基底の組の共変微分を零とする。$T(\partial_x,\partial_y)=2\partial_x$ なので $\overline\nabla=\nabla-\frac12T$ では $\overline\nabla_{\partial_x}\partial_y=\overline\nabla_{\partial_y}\partial_x=\partial_x$ となる。従って捩率は零だが、$g$ を Euclid 計量としたとき $\partial_xg(\partial_y,\partial_x)=0$ に対し $g(\overline\nabla_{\partial_x}\partial_y,\partial_x)+g(\partial_y,\overline\nabla_{\partial_x}\partial_x)=1$ となる。捩率の補正で計量の微分まで決まるわけではない。
:::

### 共標構で捩率を測る

標構 $E_a$ と双対余標構 $\theta^a$ を用い、$\nabla_XE_b=\omega^a{}_b(X)E_a$ とおく。捩率の成分を 2 形式で表すと
$$
\theta^a(T(X,Y))=(d\theta^a+\omega^a{}_b\wedge\theta^b)(X,Y).
\label{eq-rg-torsion-in-frame}
$$
実際、$Y=\theta^b(Y)E_b$ を微分して $\theta^a(\nabla_XY)=X(\theta^a(Y))+\omega^a{}_b(X)\theta^b(Y)$ となる。$X,Y$ を交換した式と $\theta^a([X,Y])$ を引くと、最初の微分の組は $d\theta^a$ に、残る二項は外積になる。座標の余標構では $d\theta^a=0$ なので下添字の反対称部分だけでよい。一般の標構では $d\theta^a$ の項も要る。

:::example[極座標の計量接続と捩率]
平面の $\theta^1=dr,\ \theta^2=r\,d\varphi$ で $d\theta^2=r^{-1}\theta^1\wedge\theta^2$ である。標構を平行と宣言すれば $\omega=0$ なので $T(E_1,E_2)=E_2/r$ となる。一方、Euclid 接続では $\omega^2{}_1=d\varphi$、$\omega^1{}_2=-d\varphi$ である。従って第 2 成分で $d\varphi\wedge dr=-dr\wedge d\varphi$ が $d\theta^2$ を消し、第 1 成分で $-d\varphi\wedge r\,d\varphi=0$ となる。いずれの接続も計量と両立するが、後者では捩率が零になる。
:::

### 標構を変えたときの接続形式

局所標構を列ベクトルに作用する写像 $E=(E_1,\ldots,E_n)$ と書き、$\nabla E=E\omega$ とする。新しい標構を $\widetilde E=ES$ とすると、Leibniz 則より $\nabla\widetilde E=E(\omega S+dS)$ なので
$$
\widetilde\omega=S^{-1}\omega S+S^{-1}dS,\qquad
\widetilde\theta=S^{-1}\theta .
$$
第 2 式は双対性から従う。$\widetilde\theta$ を微分すると $dS^{-1}=-S^{-1}(dS)S^{-1}$ の項が生じる。これは $\widetilde\omega\wedge\widetilde\theta$ に含まれる $S^{-1}dS\wedge S^{-1}\theta$ と相殺する。その結果、捩率の成分は $d\widetilde\theta+\widetilde\omega\wedge\widetilde\theta=S^{-1}(d\theta+\omega\wedge\theta)$ と変換する。接続形式と違い、捩率の成分には $S$ の微分が残らない。

:::example[直交標構を回して生じる接続形式]
Euclid 平面の定数標構を角 $\phi(x,y)$ だけ回す行列を $S$ とする。元の $\omega$ は零なので、新しい接続形式は $S^{-1}dS=\begin{pmatrix}0&-1\\1&0\end{pmatrix}d\phi$ となる。例えば $\phi=2y$ なら $y$ 方向の輸送の成分は、移動した距離を $\ell$ として角 $-2\ell$ だけ回る。新しい標構も点ごとに正規直交なので、行列は反対称であり、成分の二乗和を保つ。
:::

### 対称化と測地線

接続の捩率を $T$ とし $\nabla^0_XY=\nabla_XY-\frac12T(X,Y)$ とおく。差はテンソルだから接続となり、その捩率は $T(X,Y)-\frac12T(X,Y)+\frac12T(Y,X)=0$ である。また曲線の速度 $V$ に対し $T(V,V)=0$ なので、曲線に沿う加速度 $D_tV$ は同じである。このことから、加速度だけを観察しても一般の接続の捩率を決められない。逆に二つの接続が各点のあらゆる初速度で加速度を等しくするなら、差 $S$ は $S(v,v)=0$ を満たす。$v+w$ で展開して $S(v,w)+S(w,v)=0$ が従うので、加速度は接続差の対称部分だけを測っている。

:::example[一定の三重積が作る捩率]
標準 $\R^3$ で $\nabla_XY=D_XY+bX\times Y$ とする。$b$ は定数である。捩率は $2bX\times Y$ なので $b=2$ で $T(e_1,e_2)=4e_3$ である。一方、曲線の加速度には $\dot\gamma\times\dot\gamma$ が入るので普通の二階微分と一致する。また $\inner{X\times Y}{Z}+\inner{Y}{X\times Z}=0$ なので、後で見る計量両立性も満たす。
:::

## 計量と両立する接続

:::definition[計量接続]{#def-rg-metric-connection}
$X(g(Y,Z))=g(\nabla_XY,Z)+g(Y,\nabla_XZ)$ を満たす接続を :term[計量接続]{en="metric connection" yomi="けいりょうせつぞく"} という。
:::

座標式を曲線に制限すると $(g(V,W))'=g(D_tV,W)+g(V,D_tW)$ になるので、平行移動は内積を保つ。逆に平行移動が内積を保つなら、任意の曲線上で平行基底を取り、その基底の成分の積を微分する。内積の係数は定数なので同じ積の微分公式を得る。任意の $X_p$ を初速度とする曲線を選んで $Y,Z$ の制限に適用すれば、定義式が成り立つ。

局所正規直交標構では接続行列 $\omega^i{}_j(X)=g(e_i,\nabla_Xe_j)$ が反対称となる。実際 $g(e_i,e_j)=\delta_{ij}$ を微分すると $\omega^i{}_j+\omega^j{}_i=0$ である。この条件から逆に、成分を微分する Leibniz 則により計量と両立する。

:::example[計量と両立しても捩率を持つ接続]
Euclid 平面で $J\partial_x=\partial_y,\ J\partial_y=-\partial_x$ とし、$\nabla_XY=D_XY+dx(X)JY$ とする。$J$ は反対称だから積の微分公式を満たす。ところが $T(\partial_x,\partial_y)=-\partial_x$ である。$\gamma(t)=(t,0)$ に沿う平行移動は $V'+JV=0$、すなわち角 $-t$ の回転である。
:::

:::example[一次元での両立条件]
$g=e^{4x}dx^2$ とする。両立条件を $\partial_x$ 三つに用いると $4e^{4x}=2\Gamma^1_{11}e^{4x}$ なので $\Gamma^1_{11}=2$ である。$\gamma(t)=t$ に沿う場 $e^{-2t}\partial_x$ の成分は減るが二乗長は $1$ のままである。
:::

### 平行移動が内積を保つこととの同値性

接続が計量と両立すれば、曲線に沿う任意の二場 $U,V$ に対して
$$
\frac{d}{dt}g(U,V)=g(D_tU,V)+g(U,D_tV)
$$
となる。証明には局所標構で係数の積を微分し、$g(E_i,E_j)$ の微分に両立式を使えばよい。従って平行な二場では右辺が零となる。逆にどの曲線でも平行移動が内積を保存するとする。点 $p$ とベクトル $X_p,U_p,V_p$ を選び、$\dot\gamma(0)=X_p$ を満たす曲線に沿って $U,V$ を平行に延ばす。その微分の左辺が零であることから $(\nabla_{X_p}g)(U_p,V_p)=0$ である。最後の等号は $(\nabla_{X_p}g)(U_p,V_p)$ が点での値で決まるため、曲線に沿う延長を使っても変わらないことによる。

### 接続差の自由度

$g$ に両立する接続 $\nabla$ が与えられたとする。差 $S(X,Y)$ を加えた接続も両立する条件は
$$
g(S(X,Y),Z)+g(Y,S(X,Z))=0
$$
である。これは各 $X$ で作用素 $S_X:Y\mapsto S(X,Y)$ が反対称であるという条件である。従って正規直交標構では接続行列の各成分は $\omega_{ab}=-\omega_{ba}$ を満たす。逆にその関係を満たす行列 1 形式を任意に選んで接続を定めれば、標構の内積は定数なので両立式が成り立つ。非正規直交な標構の計量行列を $G$ とすると関係式は $dG=\omega^\top G+G\omega$ となる。成分で右辺の $(a,b)$ は $g(\nabla E_a,E_b)+g(E_a,\nabla E_b)$ である。

:::example[拡大する標構の係数]
Euclid 直線で $E=e^{2x}\partial_x$ とする。標準接続では $\nabla_{\partial_x}E=2E$ である。係数 $2$ は零でないが、計量成分は $G=e^{4x}$ であり、$dG=4G\,dx=2G\,dx+G\cdot2dx$ を満たす。この標構で $E$ を平行と宣言すると $dG\ne0$ のため両立しない。一次元の反対称行列が零になるのは標構が正規直交しているときである。
:::

:::example[捩率を 3 形式で与える]
$H$ を 3 形式とし、$g(S(X,Y),Z)=H(X,Y,Z)/2$ で $S$ を定める。$Y,Z$ の反対称性から $\nabla+S$ は計量接続である。$\nabla$ が捩率零であれば、新しい捩率を計量で下げたものは $H$ になる。例えば標準 $\R^3$ で $H=4\,dx\wedge dy\wedge dz$ を取ると $S(e_1,e_2)=2e_3$ となり、前節の $b=2$ の接続を再現する。$S(V,V)=0$ なのでこの場合は曲線の加速度も同じである。
:::

## Levi-Civita 接続の存在と一意性（Koszul の公式）

内積の微分に現れる未知量は $g(\nabla_XY,Z)$ である。$X,Y,Z$ を順番に置き換えた三つの式を組み合わせると、捩率零の条件で未知量を一つに絞れる。

:::theorem[Levi-Civita 接続]{#thm-rg-levi-civita}
計量と両立し捩率が零の接続は一つ存在する。これを :term[Levi-Civita 接続]{en="Levi-Civita connection" yomi="れびちびたせつぞく"} という。この接続は
$$
\begin{aligned}
2g(\nabla_XY,Z)={}&Xg(Y,Z)+Yg(Z,X)-Zg(X,Y)\\
&+g([X,Y],Z)-g([Y,Z],X)+g([Z,X],Y)
\end{aligned}
\label{eq-rg-koszul}
$$
を満たす。
:::

:::proof
計量の両立式を巡回置換して最初の二式を加え、三番目を引く。$\nabla_XY-\nabla_YX=[X,Y]$ を用いて、残る三つの差を括弧に置き換えると :eqref[eq-rg-koszul] を得る。正定値性により右辺は $\nabla_XY$ を一意に定める。

逆に右辺を $K(X,Y,Z)$ とする。$Z\to fZ$ で $X(f)g(Y,Z)$ と $Y(f)g(Z,X)$ がそれぞれ $[fZ,X]$ と $[Y,fZ]$ から生じる項に打ち消される。よって $K$ は $Z$ に関数線形であり、計量でベクトルに戻せる。$X\to fX$ では $Y(f)g(Z,X)$ と $-Z(f)g(X,Y)$ が括弧項と打ち消し合う。$Y\to fY$ では $X(f)g(Y,Z)$ が二つ、$Z(f)g(X,Y)$ が相殺して $2X(f)g(Y,Z)$ が残る。従って接続の公理を満たす。さらに $K(X,Y,Z)-K(Y,X,Z)=2g([X,Y],Z)$、$K(X,Y,Z)+K(X,Z,Y)=2Xg(Y,Z)$ は、式の対称項と括弧の交代性から直接従う。この二式が捩率零と両立条件を与える。
:::

:::example[局所等長写像と接続]
局所等長写像で Levi-Civita 接続を運ぶと、内積の微分と Lie 括弧は保存されるので、運んだ接続は計量に両立し捩率が零である。一意性から像側の Levi-Civita 接続と一致する。例えば $(x,y)\mapsto(x,2y)$ は $dx^2+4dy^2$ から Euclid 計量への等長写像なので、元の座標における全接続係数は零である。
:::

定数倍 $c^2g$ では :eqref[eq-rg-koszul] の両辺に同じ倍率がかかるので接続は変わらない。積計量の持ち上げた異なる因子の場の共変微分は零である。式の括弧と交差内積が零であり、一つの因子の係数は他因子に依存しないためである。

:::example[積接続の成分]
$dr^2+4d\theta^2$ の円と直線の積で $X=r\partial_\theta$、$Y=\sin\theta\,\partial_r$ とする。成分係数が定数の積座標では $\nabla_{\partial_\theta}\partial_r=0$ である。従って $\nabla_XY=r\cos\theta\,\partial_r$ となる一方、$\nabla_YX=\sin\theta\,\partial_\theta$ である。その差が $[X,Y]$ と一致するのは第 2 引数の係数だけが微分されるためである。
:::

### はめ込み上での構成

誘導計量のはめ込み $F:M\to(N,h)$ を考える。$h$ の Levi-Civita 接続で $dFY$ を $X$ 方向に微分し、$dF(TM)$ に正射影して $dF$ で戻す。この操作が誘導計量の接続である。実際、正射影は関数線形なので接続の公理を保つ。さらに接ベクトルに内積を取ると法成分が消えるから両立式も保つ。局所座標 $x,y$ を用いると、引き戻し微分の反対称部分
$$
D_i(dF(\partial_j))^\alpha-D_j(dF(\partial_i))^\alpha
=\partial_i\partial_j F^\alpha-\partial_j\partial_i F^\alpha
+\Gamma^\alpha_{\beta\gamma}
(\partial_iF^\beta\partial_jF^\gamma-\partial_jF^\beta\partial_iF^\gamma)
$$
は混合偏微分と $\Gamma^\alpha_{\beta\gamma}$ の対称性によって零である。従って座標場で捩率が零となる。捩率はテンソルなので任意の場で零となり、:ref[thm-rg-levi-civita] の一意性が適用できる。

:::example[球面上で接方向だけを残す]
半径 $a$ の球面を Euclid 空間に入れ、$p$ を位置ベクトルとする。接場 $Y$ について $\inner{Y}{p}=0$ を $X$ 方向に微分すると $\inner{D_XY}{p}=-\inner{X}{Y}$ となる。よって
$$
\nabla_XY=D_XY+\frac{\inner{X}{Y}}{a^2}p.
$$
例えば曲線 $p(t)=(a\cos(t/a),a\sin(t/a),0)$ の速度は単位で、Euclid 加速度は $-p/a^2$ である。式に $X=Y=\dot p$ を入れると接加速度は零となる。周囲空間で加速していても球面上の接続では速度は平行である。
:::

### 接続差から Koszul の一意性を確かめる

二つの計量接続の差を $S$ とし、$A(X,Y,Z)=g(S(X,Y),Z)$ とおく。両立条件から $A(X,Y,Z)=-A(X,Z,Y)$ である。両方とも捩率零なら $A(X,Y,Z)=A(Y,X,Z)$ である。これを順に用いると
$$
A(X,Y,Z)=A(Y,X,Z)=-A(Y,Z,X)=-A(Z,Y,X)
=A(Z,X,Y)=A(X,Z,Y)=-A(X,Y,Z)
$$
なので $A=0$ となる。同じ方法は捩率を指定した計量接続の構成にも使える。実際、Levi-Civita 接続との差に対し、$t(X,Y,Z)=g(T(X,Y),Z)$ と書くと $t(X,Y,Z)=A(X,Y,Z)+A(Y,Z,X)$ である。これを巡回した三式から
$$
2A(X,Y,Z)=t(X,Y,Z)-t(Y,Z,X)+t(Z,X,Y)
$$
が従う。右辺は $Y,Z$ の交換で符号が変わり、$X,Y$ の反対称部分から $t$ が復元されるので、任意の捩率テンソル $T$ を持つ計量接続が一意に作れる。ここで自由なのは $T$ の各点での選択であって、指定した後に接続をさらに変えることはできない。

:::example[四元数群の移動標構（尺度の二乗 1/2）]
単位四元数群の左不変計量で、虚四元数 $i,j,k$ の長さを $1,1,a$ とし、$a^2=\frac{1}{2}$ を選ぶ。正規直交標構を $(E_1,E_2,E_3)=(i,j,k/a)$ とする。左移動で延長した場の括弧は虚四元数の交換子である。実際、群の点 $q$ で $E_u(q)=qu$ とすれば、$[E_u,E_v](q)=(qu)v-(qv)u=q(uv-vu)$ となる。ここで水平とは $E_1,E_2$ の張る面を指す。

非座標標構での Koszul 公式を、$X=E_i,Y=E_j,Z=E_k$ と置いて使う。内積の成分は定数なので $2\inner{\nabla_i E_j}{E_k}=c^k_{ij}-c^i_{jk}+c^j_{ki}$ となる。括弧係数は $[E_1,E_2]=2aE_3$、$[E_2,E_3]=(2/a)E_1$、$[E_3,E_1]=(2/a)E_2$ である。例えば $2\inner{\nabla_1 E_2}{E_3}=2a-2/a+2/a$ なので $\nabla_1 E_2=aE_3$ になる。残りも同じ三つの係数を添字の順に入れ、$E_i$ 方向の接続行列を $C_i$ とすれば
$$
\begin{aligned}
C_1&=\begin{pmatrix}0 & 0 & 0\\0 & 0 & - \frac{\sqrt{2}}{2}\\0 & \frac{\sqrt{2}}{2} & 0\end{pmatrix}\\
C_2&=\begin{pmatrix}0 & 0 & \frac{\sqrt{2}}{2}\\0 & 0 & 0\\- \frac{\sqrt{2}}{2} & 0 & 0\end{pmatrix}\\
C_3&=\begin{pmatrix}0 & - \frac{3 \sqrt{2}}{2} & 0\\\frac{3 \sqrt{2}}{2} & 0 & 0\\0 & 0 & 0\end{pmatrix}
\end{aligned}
$$
となる。各行列は反対称なので、内積の成分が定数である条件と一致する。一方、$C_1E_2-C_2E_1=2aE_3$ で零ではないが、これは $[E_1,E_2]$ と一致するため捩率は零である。座標基底で捩率零を表す $\Gamma^k_{ij}=\Gamma^k_{ji}$ を、この移動標構にそのまま使ってはいけない。

例えば $E_3$ に沿う積分曲線では $D_tV$ の成分は $\dot v+C_3v$ である。$C_3$ は最初の二成分の平面内で角速度 $2/a-a$ の回転を作る。従って初期値 $E_1$ の平行移動の成分は $(\cos((2/a-a)t),-\sin((2/a-a)t),0)$ であり、実際の係数は $\frac{3 \sqrt{2}}{2}$ である。一定の Lie 代数方向に沿う曲線でも、移動標構の他の方向が平行とは限らない。
:::

:::example[四元数群の移動標構（尺度の二乗 3/4）]
単位四元数群の左不変計量で、虚四元数 $i,j,k$ の長さを $1,1,a$ とし、$a^2=\frac{3}{4}$ を選ぶ。正規直交標構を $(E_1,E_2,E_3)=(i,j,k/a)$ とする。左移動で延長した場の括弧は虚四元数の交換子である。実際、群の点 $q$ で $E_u(q)=qu$ とすれば、$[E_u,E_v](q)=(qu)v-(qv)u=q(uv-vu)$ となる。ここで水平とは $E_1,E_2$ の張る面を指す。

非座標標構での Koszul 公式を、$X=E_i,Y=E_j,Z=E_k$ と置いて使う。内積の成分は定数なので $2\inner{\nabla_i E_j}{E_k}=c^k_{ij}-c^i_{jk}+c^j_{ki}$ となる。括弧係数は $[E_1,E_2]=2aE_3$、$[E_2,E_3]=(2/a)E_1$、$[E_3,E_1]=(2/a)E_2$ である。例えば $2\inner{\nabla_1 E_2}{E_3}=2a-2/a+2/a$ なので $\nabla_1 E_2=aE_3$ になる。残りも同じ三つの係数を添字の順に入れ、$E_i$ 方向の接続行列を $C_i$ とすれば
$$
\begin{aligned}
C_1&=\begin{pmatrix}0 & 0 & 0\\0 & 0 & - \frac{\sqrt{3}}{2}\\0 & \frac{\sqrt{3}}{2} & 0\end{pmatrix}\\
C_2&=\begin{pmatrix}0 & 0 & \frac{\sqrt{3}}{2}\\0 & 0 & 0\\- \frac{\sqrt{3}}{2} & 0 & 0\end{pmatrix}\\
C_3&=\begin{pmatrix}0 & - \frac{5 \sqrt{3}}{6} & 0\\\frac{5 \sqrt{3}}{6} & 0 & 0\\0 & 0 & 0\end{pmatrix}
\end{aligned}
$$
となる。各行列は反対称なので、内積の成分が定数である条件と一致する。一方、$C_1E_2-C_2E_1=2aE_3$ で零ではないが、これは $[E_1,E_2]$ と一致するため捩率は零である。座標基底で捩率零を表す $\Gamma^k_{ij}=\Gamma^k_{ji}$ を、この移動標構にそのまま使ってはいけない。

例えば $E_3$ に沿う積分曲線では $D_tV$ の成分は $\dot v+C_3v$ である。$C_3$ は最初の二成分の平面内で角速度 $2/a-a$ の回転を作る。従って初期値 $E_1$ の平行移動の成分は $(\cos((2/a-a)t),-\sin((2/a-a)t),0)$ であり、実際の係数は $\frac{5 \sqrt{3}}{6}$ である。一定の Lie 代数方向に沿う曲線でも、移動標構の他の方向が平行とは限らない。
:::

:::example[四元数群の移動標構（尺度の二乗 3/2）]
単位四元数群の左不変計量で、虚四元数 $i,j,k$ の長さを $1,1,a$ とし、$a^2=\frac{3}{2}$ を選ぶ。正規直交標構を $(E_1,E_2,E_3)=(i,j,k/a)$ とする。左移動で延長した場の括弧は虚四元数の交換子である。実際、群の点 $q$ で $E_u(q)=qu$ とすれば、$[E_u,E_v](q)=(qu)v-(qv)u=q(uv-vu)$ となる。ここで水平とは $E_1,E_2$ の張る面を指す。

非座標標構での Koszul 公式を、$X=E_i,Y=E_j,Z=E_k$ と置いて使う。内積の成分は定数なので $2\inner{\nabla_i E_j}{E_k}=c^k_{ij}-c^i_{jk}+c^j_{ki}$ となる。括弧係数は $[E_1,E_2]=2aE_3$、$[E_2,E_3]=(2/a)E_1$、$[E_3,E_1]=(2/a)E_2$ である。例えば $2\inner{\nabla_1 E_2}{E_3}=2a-2/a+2/a$ なので $\nabla_1 E_2=aE_3$ になる。残りも同じ三つの係数を添字の順に入れ、$E_i$ 方向の接続行列を $C_i$ とすれば
$$
\begin{aligned}
C_1&=\begin{pmatrix}0 & 0 & 0\\0 & 0 & - \frac{\sqrt{6}}{2}\\0 & \frac{\sqrt{6}}{2} & 0\end{pmatrix}\\
C_2&=\begin{pmatrix}0 & 0 & \frac{\sqrt{6}}{2}\\0 & 0 & 0\\- \frac{\sqrt{6}}{2} & 0 & 0\end{pmatrix}\\
C_3&=\begin{pmatrix}0 & - \frac{\sqrt{6}}{6} & 0\\\frac{\sqrt{6}}{6} & 0 & 0\\0 & 0 & 0\end{pmatrix}
\end{aligned}
$$
となる。各行列は反対称なので、内積の成分が定数である条件と一致する。一方、$C_1E_2-C_2E_1=2aE_3$ で零ではないが、これは $[E_1,E_2]$ と一致するため捩率は零である。座標基底で捩率零を表す $\Gamma^k_{ij}=\Gamma^k_{ji}$ を、この移動標構にそのまま使ってはいけない。

例えば $E_3$ に沿う積分曲線では $D_tV$ の成分は $\dot v+C_3v$ である。$C_3$ は最初の二成分の平面内で角速度 $2/a-a$ の回転を作る。従って初期値 $E_1$ の平行移動の成分は $(\cos((2/a-a)t),-\sin((2/a-a)t),0)$ であり、実際の係数は $\frac{\sqrt{6}}{6}$ である。一定の Lie 代数方向に沿う曲線でも、移動標構の他の方向が平行とは限らない。
:::

## Christoffel 記号の座標計算

### 座標と標構

$\partial_i$ を :eqref[eq-rg-koszul] に代入すると括弧がすべて消える。両辺に逆計量を掛ければ
$$
\Gamma^k_{ij}=\tfrac12g^{k\ell}(\partial_i g_{j\ell}+\partial_j g_{i\ell}-\partial_\ell g_{ij})
\label{eq-rg-christoffel}
$$
となる。これを :term[Christoffel 記号]{en="Christoffel symbol" yomi="くりすとっふぇるきごう"} という。新座標 $y$ では、$\nabla_{\partial_{y^a}}\partial_{y^b}$ を $x$ の成分に展開すると $\partial_{y^a}\partial_{y^b}x^k+\Gamma^k_{ij}(\partial_{y^a}x^i)(\partial_{y^b}x^j)$ になる。これに $\partial y^c/\partial x^k$ を掛けて新係数を得る。二階微分が加わるため、接続係数はテンソルの成分とは異なる変換をする。

:::example[回転対称計量]
$g=dr^2+f(r)^2d\theta^2$、$f>0$ とする。唯一の非定数成分の微分は $\partial_r g_{\theta\theta}=2ff'$ なので $\Gamma^r_{\theta\theta}=-ff'$、$\Gamma^\theta_{r\theta}=\Gamma^\theta_{\theta r}=f'/f$ であり、他は零である。平面の $f=r$ では $r=2$ でこれらは $-2,\frac12,\frac12$ となる。
:::

:::example[共形座標]
$\widetilde g=e^{2u}g$ を :eqref[eq-rg-christoffel] に代入する。$\partial_i(e^{2u}g_{j\ell})=e^{2u}(2u_i g_{j\ell}+\partial_i g_{j\ell})$ だから、逆行列の $e^{-2u}$ が消え、差は $\delta_i^k u_j+\delta_j^ku_i-g_{ij}g^{k\ell}u_\ell$ である。従って座標と独立に
$$
\widetilde\nabla_XY-\nabla_XY=du(X)Y+du(Y)X-g(X,Y)\grad u.
\label{eq-rg-conformal-connection}
$$
上半平面では $u=-\log y$ なので $\widetilde\Gamma^y_{xx}=1/y$、$\widetilde\Gamma^y_{yy}=\widetilde\Gamma^x_{xy}=-1/y$ である。
:::

### 点で消す座標

固定した $p$ を座標の原点とし、$x^k=y^k-\Gamma^k_{ij}(p)y^iy^j/2$ とおく。捩率が零なのでこの二階微分は $-\Gamma^k_{ij}(p)$ になる。一次微分は単位行列なので逆関数定理が適用でき、新係数は原点で零となる。計量の両立式からその点で $\partial_a g_{bc}=0$ となる。まず線形座標変換で $g_{ij}(p)=\delta_{ij}$ にできるため、両方を同時に実現できる。近傍全体で係数が零となる条件は [曲率](/math/riemannian-geometry/curvature/) で扱う。

### 交差項のある座標での計算

最初に逆計量を求め、微分が零になる成分を消し、それから上添字を上げると計算の重複を減らせる。下添字が対称であることは捩率零を反映する。変換則で二階微分が加わるため、平坦な計量でも接続記号が現れる。

:::example[座標の曲げ方が作る記号]
Euclid 平面で $X=x$、$Y=y+x^2$ を元の直交座標とする。引き戻し計量は
$$
g=(1+4x^2)\,dx^2+4x\,dx\,dy+dy^2,\qquad
G^{-1}=\begin{pmatrix}1&-2x\\-2x&1+4x^2\end{pmatrix}.
$$
計量偏微分から $\Gamma_{xx,\ell}=\partial_xg_{x\ell}-\frac12\partial_\ell g_{xx}$ を求めると、$\ell=x$ で $4x$、$\ell=y$ で $2$ である。逆行列を掛けると $\Gamma^x_{xx}=4x-4x=0$、$\Gamma^y_{xx}=-8x^2+2(1+4x^2)=2$ となる。残る $\Gamma_{xy,\ell}$ と $\Gamma_{yy,\ell}$ は全て零である。例えば $\Gamma_{xy,x}=\frac12(\partial_xg_{yx}-\partial_xg_{xy})=0$ である。これは $(X,Y)$ の二階微分 $\partial_x^2Y=2$ が唯一残ることと一致する。
:::

:::example[同じ回転対称面の異なる座標]
円錐から頂点を除いた領域の計量 $dr^2+4r^2d\theta^2$ を考える。ここで $r=e^s$ とすると $g=e^{2s}(ds^2+4d\theta^2)$ である。$(r,\theta)$ では $\Gamma^r_{\theta\theta}=-4r$、$\Gamma^\theta_{r\theta}=1/r$ である。一方、$(s,\theta)$ では共形接続の式から $\Gamma^s_{ss}=1$、$\Gamma^s_{\theta\theta}=-4$、$\Gamma^\theta_{s\theta}=1$ となる。$s$ による二階微分 $\partial_s^2r=r$ が $\Gamma^s_{ss}$ に寄与するため、元の $\Gamma^r_{rr}=0$ から単なる線形変換では新係数を得られない。
:::

:::example[円錐の頂点を除いた領域の Christoffel 記号と余標構]
$g=dr^2+f(r)^2d\theta^2$ に $f(r)=2 r$ を入れ、$r=1$ の円を考える。ここでは $r>0$ とする。角を $2\theta$ に替えると局所 Euclid 表示を得るが、角の周期は二倍になる。

まず $\partial_r g_{\theta\theta}=2ff'$、角方向の微分は零である。逆計量は $\diag(1,f^{-2})$ だから、角方向の下添字が二つのとき上添字は $r$、動径と角方向が一つずつなら上添字は $\theta$ だけ残る。従って
$$
\begin{aligned}
\Gamma^r_{\theta\theta}&=-ff'=- 4 r\\
\Gamma^\theta_{r\theta}=\Gamma^\theta_{\theta r}&=f'/f=\frac{1}{r}
\end{aligned}
$$
を得る。指定した円では $f=2$、$f'=2$ で、記号の二つの値は $-4$ と $1$ である。それ以外の記号は、定数成分と $\theta$ 微分しか含まないため零となる。

次に長さをそろえた標構 $E_1=\partial_r,\ E_2=f^{-1}\partial_\theta$ に移す。$\nabla_{\partial_\theta}E_1=(f'/f)\partial_\theta=f'E_2$ であり、$\nabla_{\partial_\theta}E_2=-f'E_1$ である。また $r$ 微分では $f^{-1}$ の微分 $-f'/f^2$ が接続記号の項と相殺するので $\nabla_{\partial_r}E_2=0$ になる。よって角 $\theta=\tau$ で進む場合の成分方程式は
$$
\frac{d}{d\tau}\begin{pmatrix}v_1\\v_2\end{pmatrix}
=2\begin{pmatrix}0&1\\-1&0\end{pmatrix}
\begin{pmatrix}v_1\\v_2\end{pmatrix} .
$$
初期値 $(1,0)$ での解は $\cos(f'(r_0)\tau)E_1-\sin(f'(r_0)\tau)E_2$ である。この式の成分を微分すると、右辺の反対称行列が作用した成分と一致する。円を弧長でたどるなら $\tau=\ell/f(r_0)$ を入れるので係数は $f'/f$ になる。座標成分、正規直交成分、パラメータの尺度を分けると、同じ接続に異なる数が現れる理由が分かる。
:::

:::example[カテノイドのくびれの外側の Christoffel 記号と余標構]
$g=dr^2+f(r)^2d\theta^2$ に $f(r)=\sqrt{r^{2} + 1}$ を入れ、$r=1$ の円を考える。半径は $r=0$ で最小になる。左右の同じ半径の円では $f'$ の符号が異なる。

まず $\partial_r g_{\theta\theta}=2ff'$、角方向の微分は零である。逆計量は $\diag(1,f^{-2})$ だから、角方向の下添字が二つのとき上添字は $r$、動径と角方向が一つずつなら上添字は $\theta$ だけ残る。従って
$$
\begin{aligned}
\Gamma^r_{\theta\theta}&=-ff'=- r\\
\Gamma^\theta_{r\theta}=\Gamma^\theta_{\theta r}&=f'/f=\frac{r}{r^{2} + 1}
\end{aligned}
$$
を得る。指定した円では $f=\sqrt{2}$、$f'=\frac{\sqrt{2}}{2}$ で、記号の二つの値は $-1$ と $\frac{1}{2}$ である。それ以外の記号は、定数成分と $\theta$ 微分しか含まないため零となる。

次に長さをそろえた標構 $E_1=\partial_r,\ E_2=f^{-1}\partial_\theta$ に移す。$\nabla_{\partial_\theta}E_1=(f'/f)\partial_\theta=f'E_2$ であり、$\nabla_{\partial_\theta}E_2=-f'E_1$ である。また $r$ 微分では $f^{-1}$ の微分 $-f'/f^2$ が接続記号の項と相殺するので $\nabla_{\partial_r}E_2=0$ になる。よって角 $\theta=\tau$ で進む場合の成分方程式は
$$
\frac{d}{d\tau}\begin{pmatrix}v_1\\v_2\end{pmatrix}
=\frac{\sqrt{2}}{2}\begin{pmatrix}0&1\\-1&0\end{pmatrix}
\begin{pmatrix}v_1\\v_2\end{pmatrix} .
$$
初期値 $(1,0)$ での解は $\cos(f'(r_0)\tau)E_1-\sin(f'(r_0)\tau)E_2$ である。この式の成分を微分すると、右辺の反対称行列が作用した成分と一致する。円を弧長でたどるなら $\tau=\ell/f(r_0)$ を入れるので係数は $f'/f$ になる。座標成分、正規直交成分、パラメータの尺度を分けると、同じ接続に異なる数が現れる理由が分かる。
:::

:::example[双曲円柱の軸から離れた位置の Christoffel 記号と余標構]
$g=dr^2+f(r)^2d\theta^2$ に $f(r)=\cosh{\left(r \right)}$ を入れ、$r=\log{\left(2 \right)}$ の円を考える。座標 $r$ は閉測地線の軸からの符号付き距離である。角方向を周期で割っても局所計算は同じである。

まず $\partial_r g_{\theta\theta}=2ff'$、角方向の微分は零である。逆計量は $\diag(1,f^{-2})$ だから、角方向の下添字が二つのとき上添字は $r$、動径と角方向が一つずつなら上添字は $\theta$ だけ残る。従って
$$
\begin{aligned}
\Gamma^r_{\theta\theta}&=-ff'=- \frac{\sinh{\left(2 r \right)}}{2}\\
\Gamma^\theta_{r\theta}=\Gamma^\theta_{\theta r}&=f'/f=\tanh{\left(r \right)}
\end{aligned}
$$
を得る。指定した円では $f=\frac{5}{4}$、$f'=\frac{3}{4}$ で、記号の二つの値は $- \frac{15}{16}$ と $\frac{3}{5}$ である。それ以外の記号は、定数成分と $\theta$ 微分しか含まないため零となる。

次に長さをそろえた標構 $E_1=\partial_r,\ E_2=f^{-1}\partial_\theta$ に移す。$\nabla_{\partial_\theta}E_1=(f'/f)\partial_\theta=f'E_2$ であり、$\nabla_{\partial_\theta}E_2=-f'E_1$ である。また $r$ 微分では $f^{-1}$ の微分 $-f'/f^2$ が接続記号の項と相殺するので $\nabla_{\partial_r}E_2=0$ になる。よって角 $\theta=\tau$ で進む場合の成分方程式は
$$
\frac{d}{d\tau}\begin{pmatrix}v_1\\v_2\end{pmatrix}
=\frac{3}{4}\begin{pmatrix}0&1\\-1&0\end{pmatrix}
\begin{pmatrix}v_1\\v_2\end{pmatrix} .
$$
初期値 $(1,0)$ での解は $\cos(f'(r_0)\tau)E_1-\sin(f'(r_0)\tau)E_2$ である。この式の成分を微分すると、右辺の反対称行列が作用した成分と一致する。円を弧長でたどるなら $\tau=\ell/f(r_0)$ を入れるので係数は $f'/f$ になる。座標成分、正規直交成分、パラメータの尺度を分けると、同じ接続に異なる数が現れる理由が分かる。
:::

:::example[回転トーラスの内側の赤道の Christoffel 記号と余標構]
$g=dr^2+f(r)^2d\theta^2$ に $f(r)=\cos{\left(r \right)} + 2$ を入れ、$r=\pi$ の円を考える。内側の赤道では角方向の半径が最小になる。外側の赤道 $r=0$ と比べると二階微分の符号が逆である。

まず $\partial_r g_{\theta\theta}=2ff'$、角方向の微分は零である。逆計量は $\diag(1,f^{-2})$ だから、角方向の下添字が二つのとき上添字は $r$、動径と角方向が一つずつなら上添字は $\theta$ だけ残る。従って
$$
\begin{aligned}
\Gamma^r_{\theta\theta}&=-ff'=\left(\cos{\left(r \right)} + 2\right) \sin{\left(r \right)}\\
\Gamma^\theta_{r\theta}=\Gamma^\theta_{\theta r}&=f'/f=- \frac{\sin{\left(r \right)}}{\cos{\left(r \right)} + 2}
\end{aligned}
$$
を得る。指定した円では $f=1$、$f'=0$ で、記号の二つの値は $0$ と $0$ である。それ以外の記号は、定数成分と $\theta$ 微分しか含まないため零となる。

次に長さをそろえた標構 $E_1=\partial_r,\ E_2=f^{-1}\partial_\theta$ に移す。$\nabla_{\partial_\theta}E_1=(f'/f)\partial_\theta=f'E_2$ であり、$\nabla_{\partial_\theta}E_2=-f'E_1$ である。また $r$ 微分では $f^{-1}$ の微分 $-f'/f^2$ が接続記号の項と相殺するので $\nabla_{\partial_r}E_2=0$ になる。よって角 $\theta=\tau$ で進む場合の成分方程式は
$$
\frac{d}{d\tau}\begin{pmatrix}v_1\\v_2\end{pmatrix}
=0\begin{pmatrix}0&1\\-1&0\end{pmatrix}
\begin{pmatrix}v_1\\v_2\end{pmatrix} .
$$
初期値 $(1,0)$ での解は $\cos(f'(r_0)\tau)E_1-\sin(f'(r_0)\tau)E_2$ である。この式の成分を微分すると、右辺の反対称行列が作用した成分と一致する。円を弧長でたどるなら $\tau=\ell/f(r_0)$ を入れるので係数は $f'/f$ になる。座標成分、正規直交成分、パラメータの尺度を分けると、同じ接続に異なる数が現れる理由が分かる。
:::

:::example[指数的に縮む円柱の端の Christoffel 記号と余標構]
$g=dr^2+f(r)^2d\theta^2$ に $f(r)=e^{- \frac{r}{2}}$ を入れ、$r=2 \log{\left(2 \right)}$ の円を考える。底の弧長を増やしたとき円周が一定の割合で縮む。座標ベクトルを使うと角方向の基底の長さが場所によって異なる。

まず $\partial_r g_{\theta\theta}=2ff'$、角方向の微分は零である。逆計量は $\diag(1,f^{-2})$ だから、角方向の下添字が二つのとき上添字は $r$、動径と角方向が一つずつなら上添字は $\theta$ だけ残る。従って
$$
\begin{aligned}
\Gamma^r_{\theta\theta}&=-ff'=\frac{e^{- r}}{2}\\
\Gamma^\theta_{r\theta}=\Gamma^\theta_{\theta r}&=f'/f=- \frac{1}{2}
\end{aligned}
$$
を得る。指定した円では $f=\frac{1}{2}$、$f'=- \frac{1}{4}$ で、記号の二つの値は $\frac{1}{8}$ と $- \frac{1}{2}$ である。それ以外の記号は、定数成分と $\theta$ 微分しか含まないため零となる。

次に長さをそろえた標構 $E_1=\partial_r,\ E_2=f^{-1}\partial_\theta$ に移す。$\nabla_{\partial_\theta}E_1=(f'/f)\partial_\theta=f'E_2$ であり、$\nabla_{\partial_\theta}E_2=-f'E_1$ である。また $r$ 微分では $f^{-1}$ の微分 $-f'/f^2$ が接続記号の項と相殺するので $\nabla_{\partial_r}E_2=0$ になる。よって角 $\theta=\tau$ で進む場合の成分方程式は
$$
\frac{d}{d\tau}\begin{pmatrix}v_1\\v_2\end{pmatrix}
=- \frac{1}{4}\begin{pmatrix}0&1\\-1&0\end{pmatrix}
\begin{pmatrix}v_1\\v_2\end{pmatrix} .
$$
初期値 $(1,0)$ での解は $\cos(f'(r_0)\tau)E_1-\sin(f'(r_0)\tau)E_2$ である。この式の成分を微分すると、右辺の反対称行列が作用した成分と一致する。円を弧長でたどるなら $\tau=\ell/f(r_0)$ を入れるので係数は $f'/f$ になる。座標成分、正規直交成分、パラメータの尺度を分けると、同じ接続に異なる数が現れる理由が分かる。
:::

## テンソル場の共変微分

### 双対と積

ペアリングの積の微分を要求すると、1 形式 $\alpha$ の共変微分は $(\nabla_X\alpha)(Y)=X(\alpha(Y))-\alpha(\nabla_XY)$ に決まる。右辺では $Y$ の微分項が相殺するので 1 形式である。同様に $(0,k)$ テンソルを引数に評価した式の積の微分から
$$
(\nabla_XB)(Y_1,\dots,Y_k)
=X(B(Y_1,\dots,Y_k))-\sum_{i=1}^kB(Y_1,\dots,\nabla_XY_i,\dots,Y_k)
$$
を得る。一般の型ではベクトルと共変ベクトルのテンソル積に Leibniz 則を適用する。基底の変換に由来する微分は双対基底の微分と相殺するので貼り合わさり、縮約はペアリングであるため微分と可換である。上添字ごとに $+\Gamma$、下添字ごとに $-\Gamma$ が現れる。例えば $(\nabla_iB)^k{}_j=\partial_iB^k{}_j+\Gamma^k_{im}B^m{}_j-\Gamma^m_{ij}B^k{}_m$ である。

:::definition[Hessian]{#def-rg-hessian}
滑らかな関数 $f$ の :term[Hessian]{en="Hessian" yomi="へっしあん"} は $\nabla df$ である。従って $\operatorname{Hess}f(X,Y)=X(Yf)-(\nabla_XY)f=g(\nabla_X\grad f,Y)$ である。最後の式は計量との両立を用いた。
:::

反対称部分は $[X,Y]f-(\nabla_XY-\nabla_YX)f=-T(X,Y)f$ なので Levi-Civita 接続では対称である。さらに $\Delta_g f=\tr_g\operatorname{Hess}f$ である。実際 $\tr(X\mapsto\nabla_XY)=\partial_iY^i+\Gamma^i_{ij}Y^j$ と書け、:eqref[eq-rg-christoffel] を縮約すると $\Gamma^i_{ij}=\frac12g^{ab}\partial_jg_{ab}=\partial_j\log\sqrt{\det G}$ となる。最後は行列式を列ごとに微分した式であり、発散の式 :eqref[eq-rg-laplace-coordinate] と一致する。

:::example[極座標の Hessian]
平面で $f=r^2/2$ とする。$rr$ 成分は $1$、$\theta\theta$ 成分は $-\Gamma^r_{\theta\theta}r=r^2$、交差成分は零である。従って $\operatorname{Hess} f=g$、$\Delta f=2$ である。普通の二階偏微分だけでは $\theta\theta$ 成分を失う。
:::

:::example[上半平面の Hessian]
上半平面で $f=\log y$ とする。$xx$ 成分は $-1/y^2$、$yy$ 成分は $-1/y^2+1/y^2=0$ である。よって $\operatorname{Hess}f=-(g-df\otimes df)$、トレースは $-1$ である。
:::

### 外微分との関係

1 形式 $\alpha$ では
$$
d\alpha(X,Y)=(\nabla_X\alpha)(Y)-(\nabla_Y\alpha)(X)+\alpha(T(X,Y)).
$$
実際、右辺の接続項をまとめると $-\alpha(\nabla_XY-\nabla_YX-T(X,Y))=-\alpha([X,Y])$ となる。一般の $k$ 形式 $\beta$ でも捩率零なら
$$
d\beta(X_0,\ldots,X_k)=\sum_i(-1)^i(\nabla_{X_i}\beta)(X_0,\ldots,\widehat X_i,\ldots,X_k).
$$
両辺はテンソルなので座標場で確認すればよい。外微分の括弧項は座標場の可換性で零である。共変微分側では、$i<j$ ごとに $\nabla_{\partial_i}\partial_j$ と $\nabla_{\partial_j}\partial_i$ を引数に含む項が反対符号で現れる。捩率零ゆえ二つのベクトルが等しく、残るのは係数の偏微分の交代和である。

:::example[閉形式と平行形式]
極座標の $\alpha=dr$ は閉 1 形式である。平面の接続で
$$
(\nabla_{\partial_\theta}dr)(\partial_\theta)=-dr(-r\partial_r)=r
$$
となるので $\alpha$ は平行でない。閉性が強制するのは共変微分の反対称成分が零であることだけである。ここでは $\nabla dr=r\,d\theta^2$ が対称なので、この条件を満たす。
:::

### 共変微分と体積

向き付けた Riemann 多様体の体積形式が Levi-Civita 接続で平行であることを確かめる。正に向いた正規直交標構 $E_i$ 上で、$\omega_g(E_1,\ldots,E_n)=1$ である。共変微分を展開すると、評価値の微分は零で、各引数を $\nabla_XE_i$ に変えた項では $E_i$ 成分だけが残る。その係数は $g(\nabla_XE_i,E_i)=\frac12X g(E_i,E_i)=0$ である。従って $\nabla\omega_g=0$ である。同じ議論で、任意の計量接続は計量が誘導する体積形式を平行に保つ。向きがない場合も局所での形式の符号が変わるだけなので、密度として平行である。

:::example[交差項があっても体積を保つ]
座標変換 $(X,Y)=(x,y+x^2)$ の例では $\det G=1$ なので体積形式は $dx\wedge dy$ である。接続記号は $\Gamma^y_{xx}=2$ だけなので $\nabla_{\partial_x}dx=0$、$\nabla_{\partial_x}dy=-2dx$ である。従って
$$
\nabla_{\partial_x}(dx\wedge dy)=dx\wedge(-2dx)=0.
$$
$\partial_y$ 方向では両方の 1 形式が平行である。個々の共標構が平行でなくても、最高次の外積は平行になり得る。
:::

:::example[線形の共形指数で Hessian と跡を分ける]
平面に $g=e^{2u}(dx^2+dy^2)$、$u=x$ を入れ、$f=x^2+y^2$ を点 $\left( 1, \  0\right)$ で微分する。指数の勾配が一定で、二階微分は零である。ここでは接続による補正が Hessian の成分と跡にどう現れるかを調べる。

共形接続差の式で、上添字 $x,y$ ごとに下添字を行列に並べると
$$
\begin{aligned}
(\Gamma^x_{ij})&=\begin{pmatrix}1 & 0\\0 & -1\end{pmatrix}\\
(\Gamma^y_{ij})&=\begin{pmatrix}0 & 1\\1 & 0\end{pmatrix}
\end{aligned}
$$
である。例えば $\Gamma^x_{xy}=u_y$、$\Gamma^y_{xx}=-u_y$ となる。符号の違いは前者では $\delta_x^x u_y$、後者では $-\delta_{xx}u_y$ が残るためである。指定点では $du$ の成分が $\begin{pmatrix}1\\0\end{pmatrix}$、計量行列が $\begin{pmatrix}e^{2} & 0\\0 & e^{2}\end{pmatrix}$ となる。

普通の二階微分は $2I$ であり、共変 Hessian はこれから $f_x(\Gamma^x_{ij})+f_y(\Gamma^y_{ij})$ を引く。従って指定点での結果は
$$
\begin{aligned}
\operatorname{Hess} f&=\begin{pmatrix}0 & 0\\0 & 4\end{pmatrix}\\
\Delta_g f&=e^{-2u}\tr(\operatorname{Hess} f)=\frac{4}{e^{2}}
\end{aligned}
$$
となる。二次元では二つの接続行列の通常の跡がともに零であり、Hessian の補正項の跡は消える。従って $\Delta_g f$ の式には $du$ が見えない。しかし Hessian の成分を求める場合、その補正まで消えるとは限らない。例えば単位ベクトルの評価には逆計量の倍率も必要であり、$g$ 単位方向 $e^{-u}\partial_x$ に対する値は $0$ である。
:::

:::example[原点で対称な二次式で Hessian と跡を分ける]
平面に $g=e^{2u}(dx^2+dy^2)$、$u=x^{2} + y^{2}$ を入れ、$f=x^2+y^2$ を点 $\left( 0, \  0\right)$ で微分する。原点では一階微分が零であるが、二階微分は零でない。点での記号と曲率を区別する例である。

共形接続差の式で、上添字 $x,y$ ごとに下添字を行列に並べると
$$
\begin{aligned}
(\Gamma^x_{ij})&=\begin{pmatrix}2 x & 2 y\\2 y & - 2 x\end{pmatrix}\\
(\Gamma^y_{ij})&=\begin{pmatrix}- 2 y & 2 x\\2 x & 2 y\end{pmatrix}
\end{aligned}
$$
である。例えば $\Gamma^x_{xy}=u_y$、$\Gamma^y_{xx}=-u_y$ となる。符号の違いは前者では $\delta_x^x u_y$、後者では $-\delta_{xx}u_y$ が残るためである。指定点では $du$ の成分が $\begin{pmatrix}0\\0\end{pmatrix}$、計量行列が $\begin{pmatrix}1 & 0\\0 & 1\end{pmatrix}$ となる。

普通の二階微分は $2I$ であり、共変 Hessian はこれから $f_x(\Gamma^x_{ij})+f_y(\Gamma^y_{ij})$ を引く。従って指定点での結果は
$$
\begin{aligned}
\operatorname{Hess} f&=\begin{pmatrix}2 & 0\\0 & 2\end{pmatrix}\\
\Delta_g f&=e^{-2u}\tr(\operatorname{Hess} f)=4
\end{aligned}
$$
となる。二次元では二つの接続行列の通常の跡がともに零であり、Hessian の補正項の跡は消える。従って $\Delta_g f$ の式には $du$ が見えない。しかし Hessian の成分を求める場合、その補正まで消えるとは限らない。例えば単位ベクトルの評価には逆計量の倍率も必要であり、$g$ 単位方向 $e^{-u}\partial_x$ に対する値は $2$ である。
:::

:::example[調和な二次式で Hessian と跡を分ける]
平面に $g=e^{2u}(dx^2+dy^2)$、$u=x^{2} - y^{2}$ を入れ、$f=x^2+y^2$ を点 $\left( 1, \  1\right)$ で微分する。二方向の二階微分は反対符号なので平面 Laplacian が零になる。計量の係数は定数でない。

共形接続差の式で、上添字 $x,y$ ごとに下添字を行列に並べると
$$
\begin{aligned}
(\Gamma^x_{ij})&=\begin{pmatrix}2 x & - 2 y\\- 2 y & - 2 x\end{pmatrix}\\
(\Gamma^y_{ij})&=\begin{pmatrix}2 y & 2 x\\2 x & - 2 y\end{pmatrix}
\end{aligned}
$$
である。例えば $\Gamma^x_{xy}=u_y$、$\Gamma^y_{xx}=-u_y$ となる。符号の違いは前者では $\delta_x^x u_y$、後者では $-\delta_{xx}u_y$ が残るためである。指定点では $du$ の成分が $\begin{pmatrix}2\\-2\end{pmatrix}$、計量行列が $\begin{pmatrix}1 & 0\\0 & 1\end{pmatrix}$ となる。

普通の二階微分は $2I$ であり、共変 Hessian はこれから $f_x(\Gamma^x_{ij})+f_y(\Gamma^y_{ij})$ を引く。従って指定点での結果は
$$
\begin{aligned}
\operatorname{Hess} f&=\begin{pmatrix}-6 & 0\\0 & 10\end{pmatrix}\\
\Delta_g f&=e^{-2u}\tr(\operatorname{Hess} f)=4
\end{aligned}
$$
となる。二次元では二つの接続行列の通常の跡がともに零であり、Hessian の補正項の跡は消える。従って $\Delta_g f$ の式には $du$ が見えない。しかし Hessian の成分を求める場合、その補正まで消えるとは限らない。例えば単位ベクトルの評価には逆計量の倍率も必要であり、$g$ 単位方向 $e^{-u}\partial_x$ に対する値は $-6$ である。
:::

:::example[混合二次式で Hessian と跡を分ける]
平面に $g=e^{2u}(dx^2+dy^2)$、$u=x y$ を入れ、$f=x^2+y^2$ を点 $\left( 1, \  2\right)$ で微分する。対角の二階微分は零だが、混合二階微分は零でない。Hessian の跡の情報と行列の情報を分ける。

共形接続差の式で、上添字 $x,y$ ごとに下添字を行列に並べると
$$
\begin{aligned}
(\Gamma^x_{ij})&=\begin{pmatrix}y & x\\x & - y\end{pmatrix}\\
(\Gamma^y_{ij})&=\begin{pmatrix}- x & y\\y & x\end{pmatrix}
\end{aligned}
$$
である。例えば $\Gamma^x_{xy}=u_y$、$\Gamma^y_{xx}=-u_y$ となる。符号の違いは前者では $\delta_x^x u_y$、後者では $-\delta_{xx}u_y$ が残るためである。指定点では $du$ の成分が $\begin{pmatrix}2\\1\end{pmatrix}$、計量行列が $\begin{pmatrix}e^{4} & 0\\0 & e^{4}\end{pmatrix}$ となる。

普通の二階微分は $2I$ であり、共変 Hessian はこれから $f_x(\Gamma^x_{ij})+f_y(\Gamma^y_{ij})$ を引く。従って指定点での結果は
$$
\begin{aligned}
\operatorname{Hess} f&=\begin{pmatrix}2 & -10\\-10 & 2\end{pmatrix}\\
\Delta_g f&=e^{-2u}\tr(\operatorname{Hess} f)=\frac{4}{e^{4}}
\end{aligned}
$$
となる。二次元では二つの接続行列の通常の跡がともに零であり、Hessian の補正項の跡は消える。従って $\Delta_g f$ の式には $du$ が見えない。しかし Hessian の成分を求める場合、その補正まで消えるとは限らない。例えば単位ベクトルの評価には逆計量の倍率も必要であり、$g$ 単位方向 $e^{-u}\partial_x$ に対する値は $\frac{2}{e^{4}}$ である。
:::

## 演習問題

::::exercise{level=basic}
$\Gamma^1_{11}=3$ の直線上で $X=x\partial_x,\ Y=x^2\partial_x$ とする。$x=2$ における $\nabla_XY$ を求めよ。

:::solution
$X(x^2)=2x^2$ であり、接続項は $3x^3$ である。$2x^2+3x^3$ に $x=2$ を代入すると $32\partial_x$ である。
:::
::::

::::exercise{level=basic}
$\Gamma^1_{11}=2$ の直線で $\gamma(t)=2t$ に沿って $V(0)=\partial_x$ を平行移動せよ。

:::solution
:eqref[eq-rg-curve-derivative] から $v'+4v=0$ となる。$e^{4t}v$ の導関数が零だから $v(t)=e^{-4t}$ である。
:::
::::

::::exercise{level=basic}
$g=dr^2+e^{2r}d\theta^2$ の非零 Christoffel 記号を求めよ。

:::solution
$f=e^r$ として $f'=f$ を代入すると $\Gamma^r_{\theta\theta}=-e^{2r}$、$\Gamma^\theta_{r\theta}=\Gamma^\theta_{\theta r}=1$ である。他の記号では :eqref[eq-rg-christoffel] の計量微分が全て零になる。
:::
::::

::::exercise{level=basic}
前問の計量で $f=r^2$ の Hessian と Laplacian を求めよ。

:::solution
$rr$ 成分は $2$、$\theta\theta$ 成分は $2re^{2r}$、交差成分は零である。逆計量で縮約すると $2+2r$ を得る。
:::
::::

::::exercise{level=advanced}
二つの計量接続の差 $A(X,Y)$ は固定した $X$ に対し $A_X$ が反対称であることと同値である。この意味を精密に述べ、証明せよ。

:::solution
$\nabla$ を計量接続、$A$ を $(1,2)$ テンソルとすると、$\nabla+A$ が計量接続であるための必要十分条件は $g(A_XY,Z)+g(Y,A_XZ)=0$ である。両立条件の二式を引けば必要性が得られ、逆にこの和が零なら $\nabla+A$ の両立式が $\nabla$ の式に還元される。
:::
::::

::::exercise{level=advanced}
向き付けられた Riemann 多様体の体積形式が Levi-Civita 接続で平行であることを示せ。

:::solution
正規直交標構上では $\omega_g(e_1,\ldots,e_n)=1$ である。共変微分の引数置換項は、置換するベクトルの同じ添字の成分だけが寄与するので $-\sum_i\omega^i{}_i(X)$ を得る。接続行列は反対称で対角が零だから、体積形式の共変微分は零である。
:::
::::

::::exercise{level=advanced}
Euclid 平面で直交座標を $(x,y)$ とする。同じ計量を新座標 $(u,v)=(x,y+x^2)$ で表すとき、$\Gamma^v_{uu}$ を求めよ。

:::solution
逆変換は $x=u,\ y=v-u^2$ で、旧係数は零である。$\partial_u^2y=-2$ だけが二階微分に寄与する。$\partial_yv=1$ なので新係数は $-2$ である。
:::
::::

::::exercise{level=advanced}
$X$ をベクトル場とし、$A(Y)=\nabla_YX$ とおく。$\mathcal L_Xg(Y,Z)=g(A(Y),Z)+g(Y,A(Z))$ を示せ。

:::solution
Lie 微分は $Xg(Y,Z)-g([X,Y],Z)-g(Y,[X,Z])$ である。捩率零から $[X,Y]=\nabla_XY-\nabla_YX$ と書き、計量の積の微分を代入すると $\nabla_X$ を含む項が相殺し、所望の二項が残る。
:::
::::

## 確認テスト

::::quiz
接続の定義に含まれる式を選べ。

- [x] $\nabla_{fX}Y=f\nabla_XY$
- [ ] $\nabla_X(fY)=f\nabla_XY$

:::explanation
第 1 引数は関数線形である。第 2 引数では $X(f)Y$ が加わる。例えば通常の微分で定数場に非定数関数を掛ければこの項が必要になる。
:::
::::

::::quiz{answer="-3" tol="1e-6"}
$g=dr^2+r^2d\theta^2$ の $r=3$ で $\Gamma^r_{\theta\theta}$ を求めよ。

:::explanation
$g^{rr}=1$ なので $\Gamma^r_{\theta\theta}=-\partial_r(r^2)/2=-r=-3$ である。
:::
::::

::::quiz{answer="2" tol="1e-6"}
Euclid 平面で $f=x^2$ の Hessian のトレースはいくらか。

:::explanation
Cartesian 座標で接続係数は零である。二階微分行列は $\diag(2,0)$ なのでトレースは $2$ である。
:::
::::

::::quiz
Levi-Civita 接続の特徴として正しいものをすべて選べ。

- [x] 捩率は零である。
- [x] 平行移動は内積を保つ。
- [ ] 任意の座標で係数は零である。

:::explanation
捩率零と計量との両立が定義の条件である。計量の両立式を平行場に適用すると内積の時間微分は零である。平面でも極座標では $\Gamma^r_{\theta\theta}=-r$ なので最後は誤りである。
:::
::::

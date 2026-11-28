---
title: '曲率'
description: '接続の交換子から曲率を構成し、代数的対称性、曲率の縮約、共形変化を調べる。'
notation:
  - symbol: 'R(X,Y)Z'
    meaning: '曲率作用素'
  - symbol: '\mathcal R(X,Y,Z,W)'
    meaning: 'g(R(X,Y)W,Z)'
  - symbol: 'K'
    meaning: '断面曲率'
  - symbol: '\mathrm{Ric},\,S'
    meaning: 'Ricci 曲率とスカラー曲率'
  - symbol: 'W'
    meaning: 'Weyl テンソル'
---

## この章の目標

- 曲率がテンソルであることと対称性を証明する。
- 平坦性を平行標構から局所等長性につなげる。
- 断面曲率、Ricci 曲率、Weyl テンソルを区別する。
- 球面、双曲平面、積の曲率を計算する。

## Riemann 曲率テンソルとその対称性

同じ場を二方向に微分し、順序を交換する。方向場の括弧そのものが生む差を引くと、接続が作る点ごとの量が残る。

:::definition[曲率テンソル]{#def-rg-curvature}
:term[Riemann 曲率テンソル]{en="Riemann curvature tensor" yomi="りーまんきょくりつてんそる"} を
$$
R(X,Y)Z=\nabla_X\nabla_YZ-\nabla_Y\nabla_XZ-\nabla_{[X,Y]}Z
\label{eq-rg-curvature-convention}
$$
で定める。4 添字の形は $\mathcal R(X,Y,Z,W)=g(R(X,Y)W,Z)$ とする。この規約では球面の $g(R(X,Y)Y,X)$ が正である。
:::

$[fX,Y]=f[X,Y]-Y(f)X$ を代入すると第 1 引数の余計な項が消える。交代性で第 2 引数も関数線形になる。第 3 引数に $fZ$ を代入すると一次微分項は打ち消され、残る $X(Yf)-Y(Xf)-[X,Y]f$ は零である。従って点で決まるテンソルである。

:::proposition[Bianchi 恒等式と対称性]{#prop-rg-curvature-symmetries}
$\mathcal R$ は最初と最後の添字対で交代的で、対を交換しても不変である。さらに $\sum_{\mathrm{cyc}(X,Y,Z)} R(X,Y)Z=0$ および $\sum_{\mathrm{cyc}(X,Y,Z)}(\nabla_X R)(Y,Z)=0$ が成り立つ。
:::

:::proof
最初の交代性は定義から従う。計量の両立によって $XYg(Z,W)-YXg(Z,W)-[X,Y]g(Z,W)$ を展開すると、交差項が相殺し $g(R(X,Y)Z,W)+g(Z,R(X,Y)W)=0$ が残る。

点 $p$ で正規座標を用いると $\Gamma=0$ であり、$R(\partial_i,\partial_j)\partial_\ell$ の成分は $\partial_i\Gamma^m_{j\ell}-\partial_j\Gamma^m_{i\ell}$ になる。下添字の対称性により巡回和は項ごとに相殺する。さらに :eqref[eq-rg-christoffel] を代入する。$g_{ab,cd}=\partial_c\partial_dg_{ab}$ と書けば
$$
2\mathcal R_{ijk\ell}
=g_{jk,i\ell}-g_{j\ell,ik}-g_{ik,j\ell}+g_{i\ell,jk}.
$$
これは $(i,j)$ と $(k,\ell)$ の交換で不変である。最後に同じ座標で $\nabla R$ を計算すると、$\Gamma\Gamma$ の微分は点 $p$ で零なので、微分二回の $\Gamma$ の差だけが残る。その巡回和は混合微分の可換性によって零となる。
:::

:::example[一定曲率型のテンソル]
$R(X,Y)Z=k(g(Y,Z)X-g(X,Z)Y)$ とおく。巡回和では同じ内積が同じベクトルの係数として反対符号で現れる。$k=2$、$X=e_1+e_2$、$Y=e_2$、$Z=e_2$ なら $R(X,Y)Z=2e_1$ である。
:::

:::example[曲率作用素の核]
三次元の定曲率 $k=2$ で $X=e_1,Y=e_2$ を取り、$R(X,Y)$ を $\operatorname{span}(e_1,e_2)$ に制限すると行列は $\begin{pmatrix}0&2\\-2&0\end{pmatrix}$ である。$R(X,Y)e_3=0$ なので、反対称作用素として回転する面と固定する軸が分かれる。$R(X,Y)Y=2e_1$ となるので正の曲率を表す。
:::

:::example[球面の接射影から曲率を求める]
半径 $a$ の球面の位置ベクトルを $p$ とすれば、[接続の章](/math/riemannian-geometry/connections/) の接射影の式は $D_XY=\nabla_XY-g(X,Y)p/a^2$ である。そこで周囲の曲率が零である式 $D_XD_YZ-D_YD_XZ-D_{[X,Y]}Z=0$ を接方向へ射影する。$D_X(\nabla_YZ)$ の接成分は $\nabla_X\nabla_YZ$ である。一方、$D_X(-g(Y,Z)p/a^2)$ では係数を微分した項は法方向で、$D_Xp=X$ の項だけが接方向に残る。従って $R(X,Y)Z=a^{-2}(g(Y,Z)X-g(X,Z)Y)$ となる。この符号では直交単位対に対して $g(R(X,Y)Y,X)=a^{-2}>0$ である。
:::

### 無限小の長方形に沿う平行移動

ある点 $p$ の座標で二方向 $\partial_1,\partial_2$ を選び、辺長 $t$ の座標正方形を $+\partial_1,+\partial_2,-\partial_1,-\partial_2$ の順に回る。初期ベクトルを同じ点まで平行に運ぶ線形写像を $P_t$ とすると、
$$
P_t=I-t^2R(\partial_1,\partial_2)|_p+O(t^3)
$$
となる。符号は曲率の定義と回る向きの組で決まる。

この展開は、$p$ で接続行列を零にする正規座標で示せる。$p$ を原点とし、$\nabla_{\partial_i}Y$ の成分を $\partial_iY+A_iY$ と書く。各辺上で $A_i=O(t)$ だから、平行移動方程式の積分形から成分は初期値に対し $O(t^2)$ しか変わらない。したがって全周で、初期値 $v$ に対する変化は $-(\oint A_i\,dx^i)v+O(t^4)$ である。各行列成分ごとに平面の Stokes の公式を使うと、積分は
$$
\int_0^t\int_0^t(\partial_1A_2-\partial_2A_1)(x^1,x^2)\,dx^1\,dx^2
=t^2(\partial_1A_2-\partial_2A_1)(0)+O(t^3)
$$
となる。座標原点では $A_i=0$ なので曲率の行列は右辺の微分差であり、所望の展開を得る。座標を取り直しても $P_t$ と $R$ は同じ接空間の作用素であり、座標正方形の接する二方向だけが二次の項に現れる。

### 代数的曲率テンソルの自由度

内積空間 $V$ に対し、二つの添字対で交代的で対交換にも対称なテンソルは $\Lambda^2V$ 上の対称双線形形式 $B$ と対応する。対応は $B(u\wedge v,w\wedge z)=\mathcal R(u,v,w,z)$ である。第一 Bianchi 恒等式はさらに、$B$ の完全交代化が零という線形条件を課す。実際、交代化は
$$
\operatorname{Alt}B(u,v,w,z)=
\tfrac13\{B(u\wedge v,w\wedge z)+B(v\wedge w,u\wedge z)+B(w\wedge u,v\wedge z)\}
$$
となる。24 項を二対それぞれの交換と対の交換でまとめると、符号込みで同じ項が 8 個ずつ現れる。4 形式 $\omega$ 自身を $B_\omega$ として入れると、三つの巡回置換は偶置換なので $\operatorname{Alt}B_\omega=\omega$ である。従ってこの線形写像は全射で、第一 Bianchi 恒等式を満たすテンソルの次元は
$$
\frac{m(m+1)}2-\binom n4=\frac{n^2(n^2-1)}{12},\qquad m=\binom n2
$$
である。例えば $n=4$ で対称行列として数えると $21$ 成分だが、交代化で一つの条件が加わり $20$ 成分となる。

:::example[対称な作用素でも第一 Bianchi 恒等式が必要]
向き付けた四次元 Euclid 空間の体積形式を $B_\omega$ とみなすと、対交換は偶置換なので対称である。どの $u,v$ に対しても $B_\omega(u\wedge v,u\wedge v)=0$ であるが、正の基底で $B_\omega(e_1\wedge e_2,e_3\wedge e_4)=1$ である。巡回和は $3$ となり曲率の恒等式を満たさない。断面曲率から全テンソルを決める議論では、この恒等式が自由度を一つ取り除いている。
:::

## 平坦性と曲率

:::theorem[局所平坦性]{#thm-rg-flat}
ある開集合で $R=0$ であることと、その各点の近くに $g=\sum_i(dy^i)^2$ となる座標があることは同値である。
:::

:::proof
そのような座標では $\Gamma=0$ なので $R=0$ である。逆に座標直方体を取り、原点の正規直交基底を順に各座標軸の方向へ平行移動して標構 $E_a$ を作る。最初の $k$ 座標だけが変わる部分で $\nabla_{\partial_i}E_a=0$（$i\le k$）を帰納的に示す。構成より $\nabla_{\partial_k}E_a=0$ である。$i<k$ では $R=0$ と $[\partial_k,\partial_i]=0$ から $\nabla_{\partial_k}\nabla_{\partial_i}E_a=\nabla_{\partial_i}\nabla_{\partial_k}E_a=0$ である。初期面で $\nabla_{\partial_i}E_a=0$ なので ODE の一意性から延長先でも零になる。

この標構は内積を保つ。捩率零より $[E_a,E_b]=0$ であるから、双対共標構 $\theta^a$ について $d\theta^a(E_i,E_j)=-\theta^a([E_i,E_j])=0$ となる。座標直方体上で $y^a(x)=\int_0^1\theta^a_{tx}(x)\,dt$ とおけば、閉性から $\partial_i y^a=\int_0^1\partial_t(t\theta^a_i(tx))\,dt=\theta^a_i(x)$ である。$dy^a$ は基底なので逆関数定理が適用できる。
:::

:::example[局所性と周期]
$dr^2+r^2d\theta^2$ は $(r\cos\theta,r\sin\theta)$ を使って局所的に Euclid 計量となる。平坦トーラス $\R^2/\Z^2$ も局所的に Euclid 的であるが、コンパクトなので $\R^2$ と等長ではない。
:::

:::example[平坦円柱の閉じない座標]
$g=4\,d\theta^2+dz^2$ の円柱で $E_1=\partial_\theta/2,E_2=\partial_z$ は大域平行な正規直交標構である。双対形式 $2d\theta$ は閉だが、周方向を一周すると積分が $4\pi$ であるため大域関数の微分ではない。平坦性の証明が座標直方体を使う理由は、双対共標構の閉性から座標を得る段階にある。
:::

### 局所平坦性と大域平行移動

平坦な計量では局所平行標構の間の遷移行列は定数の直交行列である。実際、二つの平行標構 $E,F$ が $F=EQ$ の関係にあれば $\nabla F=E\,dQ$ となり、$dQ=0$ である。曲線を平坦標構の領域に分割すると、その中では成分は一定なので、輸送を決めるのは境界でのこの定数行列である。閉曲線を少しずつ変形したとき、隣接する小矩形の境界を平坦標構内に取れば、各小矩形の輸送は恒等写像になる。内部の辺は逆向きに二度現れるので、その輸送が相殺する。従って始終点を固定したホモトピーで輸送は不変である。もし開集合が単連結なら、任意の始点からの輸送が経路に依存せず、大域平行標構が得られる。

単連結な平坦多様体でも、得た閉余標構を大域関数に積分した局所等長写像が単射になるとは限らない。例えば $(0,\infty)\times\R$ の $dr^2+r^2d\theta^2$ は単連結で、$(r,\theta)\mapsto(r\cos\theta,r\sin\theta)$ は局所等長写像だが、$\theta$ に $2\pi$ を足しても像が同じになる。動径を $r\downarrow0$ に進めば長さが有限なので、この例は完備ではない。完備性を加えた場合の大域的な結論は [完備性と Hopf–Rinow の定理](/math/riemannian-geometry/completeness/) に進む。

:::example[平坦な閉曲線で反転が起こる]
$\R^2$ に等長変換 $F(x,y)=(x+2,-y)$ を作用させ、生成する群で割る。作用は自由で、任意のコンパクト集合を自身と交わる位置へ移す反復は有限個である。従って商は平坦な滑らかな面になる。$(0,0)$ を $(2,0)$ に結ぶ線分の射影は閉曲線で、被覆上の平行移動は成分を保つ。終端の接空間を始点と同一視するには $dF^{-1}=\diag(1,-1)$ が必要なので、縦方向のベクトルは符号が反転する。この閉曲線の小近傍では平坦でも、経路を一周した結果は恒等写像でない。
:::

## 断面曲率

:::definition[断面曲率]{#def-rg-sectional}
二次元部分空間 $\sigma=\operatorname{span}(X,Y)\subset T_pM$ の :term[断面曲率]{en="sectional curvature" yomi="だんめんきょくりつ"} を
$$
K(\sigma)=\frac{g(R(X,Y)Y,X)}{g(X,X)g(Y,Y)-g(X,Y)^2}
\label{eq-rg-sectional}
$$
とする。
:::

分子は $\mathcal R$ の二つの交代添字対に同じ基底変換を施した値なので、変換行列の行列式の二乗倍になる。分母は Gram 行列式なので同じ倍率で変わる。従って定義は平面だけで決まる。正規直交基底では分母は $1$ である。

:::lemma[断面曲率による決定]{#lem-rg-sectional-determines}
曲率の代数的対称性と第一 Bianchi 恒等式を満たすテンソルは、すべての $K(\sigma)$ によって一意に決まる。
:::

:::proof
差 $A$ に対して $A(x,y,x,y)=0$ と仮定する。$x$ を偏極し対交換を使うと $A(x,y,z,y)=0$ である。次に $y$ を偏極して $A(x,y,z,w)=-A(x,w,z,y)$ を得る。これは第 2 と第 4 引数の交換で反対称という条件である。第 1・第 2 と第 3・第 4 の交換でも反対称なので、これらの互換が生成する全置換について交代的になる。第一 Bianchi 恒等式では 3 項が偶置換になり $3A=0$ となる。
:::

従って全平面の $K=k$ は $R(X,Y)Z=k(g(Y,Z)X-g(X,Z)Y)$ と同値である。各点で等しいという仮定から $k$ の局所定数性が従うかどうかは、次元と縮約 Bianchi 恒等式に関係する。

:::example[非直交基底の計算]
定曲率 $k=3$ の空間で $X=e_1+e_2$、$Y=2e_2$ とすると分母は $2\cdot4-2^2=4$、分子は $12$ である。基底を $e_1,e_2$ に換えると分子は $3$ になるが断面曲率は同じである。
:::

:::example[混合した断面を取る]
$S^2(1)\times\R$ で $e_1,e_2$ を球面、$e_3$ を直線の単位方向とする。$X=e_1$、$Y=(e_2+e_3)/\sqrt2$ なら $X,Y$ は直交単位であり、分子には球面成分 $e_2/\sqrt2$ の二乗が寄与するので $K=1/2$ である。純粋断面での $1$ と混合断面での $0$ は典型値だが、それらの中間値も取る。
:::

### 楔積と曲率作用素

$\Lambda^2T_pM$ に $\inner{u\wedge v}{w\wedge z}=g(u,w)g(v,z)-g(u,z)g(v,w)$ を入れる。正規直交基底から作る $e_i\wedge e_j$（$i<j$）が正規直交するので正定値である。対交換の対称性を用いて $\inner{\mathscr R(u\wedge v)}{w\wedge z}=\mathcal R(u,v,w,z)$ で自己共役作用素 $\mathscr R$ を定める。従って断面曲率は単位の分解可能な 2 ベクトル $u\wedge v$ 上の二次形式の値である。一方、曲率作用素の最小固有値を測るには $\Lambda^2$ の全単位ベクトルを使う。この二種類の制約は四次元以上で異なる。例えば $e_1\wedge e_2+e_3\wedge e_4$ の自身との外積は $2e_1\wedge e_2\wedge e_3\wedge e_4$ で零でない。分解可能なら自身との外積は零なので、このベクトルはどの平面からも得られない。

:::example[四次元の断面曲率と作用素の符号の差]
向き付けた四次元 Euclid 空間の $\Lambda^2$ を Hodge スターの $+1,-1$ 固有空間に分ける。ここではスターは $e_1\wedge e_2$ を $e_3\wedge e_4$ に移すように、向きと外積で定める。各固有空間は三次元であり、正規直交基底に関する実固有値を $+1$ 側で $(-1,2,2)$、$-1$ 側で $(1,1,1)$ とする対称作用素 $A$ を考える。第一 Bianchi 恒等式は $\tr(A\circ *)=0$ と同値である。実際、標準基底 $e_i\wedge e_j$ でこのトレースを計算すると、完全交代化が $(e_1,e_2,e_3,e_4)$ に取る値の $6$ 倍となる。この例では両側の固有値の和が $3$ なので条件を満たす。

分解可能な単位 2 ベクトル $\xi=\xi_++\xi_-$ について、$\xi\wedge\xi=0$ から $\inner{\xi}{*\xi}=0$ となり、$\norm{\xi_+}^2=\norm{\xi_-}^2=1/2$ である。従って $\inner{A\xi}{\xi}\ge-1/2+1/2=0$ であるが、$A$ の最小固有値は $-1$ である。ここでは点ごとの代数を見ている。断面曲率の符号だけで作用素の正値性を置き換えることができないと分かる。
:::

## Ricci 曲率とスカラー曲率

一つの方向 $X$ の周りの断面曲率を合計すると、点の近くの体積の変化を扱える。

:::definition[Ricci 曲率とスカラー曲率]{#def-rg-ricci}
正規直交基底 $e_i$ に対して :term[Ricci 曲率]{en="Ricci curvature" yomi="りっちきょくりつ"} を $\mathrm{Ric}(Y,Z)=\sum_i g(R(e_i,Y)Z,e_i)$、:term[スカラー曲率]{en="scalar curvature" yomi="すからーきょくりつ"} を $S=\sum_i\mathrm{Ric}(e_i,e_i)$ とする。
:::

最初は線形写像のトレース、次は計量での縮約なので基底に依存しない。対交換から $\mathrm{Ric}$ は対称である。単位 $X=e_1$ に対して $\mathrm{Ric}(X,X)=\sum_{j=2}^nK(X,e_j)$、さらに $S=2\sum_{i<j}K(e_i,e_j)$ である。

:::lemma[縮約 Bianchi 恒等式]{#lem-rg-contracted-bianchi}
$\nabla^i\mathrm{Ric}_{ij}=\tfrac12\partial_j S$ が成り立つ。
:::

:::proof
正規直交基底を原点で持つ正規座標で第二恒等式を $\nabla_a\mathcal R_{bcde}+\nabla_b\mathcal R_{cade}+\nabla_c\mathcal R_{abde}=0$ と書く。$a=d$ として和を取ると $\sum_a\nabla_a\mathcal R_{bcae}-\nabla_b\mathrm{Ric}_{ce}+\nabla_c\mathrm{Ric}_{be}=0$ である。ここで $b=e$ として和を取ると、第 1 項が $-\nabla^a\mathrm{Ric}_{ca}$ となるので $-2\nabla^a\mathrm{Ric}_{ca}+\nabla_c S=0$ を得る。
:::

:::example[同じスカラー曲率を持つ異なる幾何]
$S^2(1)\times\R$ では $\mathrm{Ric}$ は正規直交基底で $\diag(1,1,0)$、$S=2$ である。三次元定曲率 $1/3$ の球面でも $S=6/3=2$ だが $\mathrm{Ric}=(2/3)g$ である。積の混合曲率が零であることは、積接続で各因子の係数が独立であることから従う。
:::

:::example[Ricci 曲率を足す方向]
半径 $2$ の三次元球面で、単位方向 $e_1$ に対し $K(e_1,e_2)=K(e_1,e_3)=1/4$ である。従って $\mathrm{Ric}(e_1,e_1)=1/2$ で、$S=3/2$ である。断面曲率を全体で足すときは $(i,j)$ と $(j,i)$ の双方が現れるため、独立な三断面の和 $3/4$ を二倍する。
:::

### 単位球面上の平均

点 $p$ の単位球面の回転不変確率測度を $d\nu$ とする。任意の対称形式 $A$ に対し
$$
\int A(u,u)\,d\nu(u)=\frac{\tr_g A}{n}
$$
である。正規直交座標で $u^iu^j$ の積分は $i\ne j$ なら一方の成分の反射で符号が変わるので零である。対角の積分は軸を交換しても同じで、和は $\int\norm u^2\,d\nu=1$ だからそれぞれ $1/n$ である。$A=\mathrm{Ric}$ に適用すると平均は $S/n$ である。また $u$ を固定し、$u^\perp$ の単位球面で $v\mapsto K(u,v)$ を平均すると $\mathrm{Ric}(u,u)/(n-1)$ となる。なぜなら計量で $J_u(v)=R(v,u)u$ と置くと $K(u,v)=g(J_uv,v)$ で、その跡が $\mathrm{Ric}(u,u)$ だからである。両方を平均した断面曲率は $S/(n(n-1))$ である。

:::example[負の見える方向があっても平均は正になる]
$S^2(1)\times H^2(2)$ で、後者の曲率を $-1/4$ とする。正規直交基底で Ricci の固有値は $(1,1,-1/4,-1/4)$ なので $S=3/2$ である。方向 $u=(e_1+e_3)/\sqrt2$ では $\mathrm{Ric}(u,u)=(1-1/4)/2=3/8$ である。$u^\perp$ の基底に $e_2,e_4,(e_1-e_3)/\sqrt2$ を取ると、対応する断面曲率は $1/2,-1/8,0$ となり和は $3/8$ である。従って微分量をある方向に固定した平均と、全方向での平均には異なる情報が入る。
:::

### 多様体の直積と曲率の縮約

積の接続は成分ごとに作用するので、座標で別の因子の方向に接続係数を微分しても零である。曲率は $\mathcal R=\pi_1^*\mathcal R_1+\pi_2^*\mathcal R_2$ となる。混合断面の消失に加え、$\mathrm{Ric}=\pi_1^*\mathrm{Ric}_1+\pi_2^*\mathrm{Ric}_2$、$S=S_1+S_2$ を得る。縮約に使う正規直交基底を二因子の基底の和として取れば、片方の因子の曲率にもう片方の基底を入れた項は零になるからである。この計算は一定曲率の仮定を必要としない。

### 左不変な標構から断面曲率を求める

[接続の章の四元数群の例](/math/riemannian-geometry/connections/) と同じ正規直交標構 $(E_1,E_2,E_3)=(i,j,k/a)$ を使う。$J_i$ をこの向き付き標構における $E_i$ との外積作用素とする。Koszul の公式で得た接続行列は $C_1=aJ_1$、$C_2=aJ_2$、$C_3=(2/a-a)J_3$ である。括弧係数を $[E_i,E_j]=c^k_{ij}E_k$ とすると、行列の係数が一定なので曲率行列は $[C_i,C_j]-c^k_{ij}C_k$ となる。$[J_i,J_j]=J_k$ は巡回順の $i,j,k$ について、各基底ベクトルに二つの外積を順に作用させれば確認できる。従って
$$
R(E_1,E_2)=(3a^2-4)J_3,\qquad
R(E_2,E_3)=-a^2J_1,\qquad R(E_3,E_1)=-a^2J_2 .
$$
例えば最初の係数は $a^2-2a(2/a-a)=3a^2-4$ である。$J_3E_2=-E_1$ だから $K_{12}=4-3a^2$ となる。他の二つも同じ順序で評価すれば $K_{23}=K_{31}=a^2$ である。各 $R(E_i,E_j)$ はその二方向の平面以外を消すので、曲率作用素は $\Lambda^2$ の対応する基底で対角になる。

:::example[三次元の左不変計量の縮約（尺度の二乗 3/4）]
単位四元数群の左不変計量で、虚四元数 $i,j,k$ の長さを $1,1,a$ とし、$a^2=\frac{3}{4}$ を選ぶ。正規直交標構を $(E_1,E_2,E_3)=(i,j,k/a)$ とする。ここで水平とは $E_1,E_2$ の張る面を指す。

正規直交な二方向に掛かる曲率は、$12$ 断面に $4-3a^2$、$13$ と $23$ の断面に $a^2$ である。曲率作用素はこの標構の二形式の基底で対角だから、任意の断面曲率はこれらの重み付き平均である。例えば $E_1$ と $Y=(E_2+E_3)/\sqrt2$ の断面では、二つの二形式の成分の二乗がともに $1/2$ なので、$K=(4-3a^2+a^2)/2=2-a^2=\frac{5}{4}$ になる。従って
$$
\begin{aligned}
(K_{12},K_{13},K_{23})&=\left( \frac{7}{4}, \  \frac{3}{4}, \  \frac{3}{4}\right)\\
(K_{\min},K_{\max})&=\left( \frac{3}{4}, \  \frac{7}{4}\right)\\
K_{\min}/K_{\max}&=\frac{3}{7}
\end{aligned}
$$
である。絶対曲率を正規化するには計量を $K_{\max}$ 倍する。曲率はその逆数倍になるから、最大が $1$、最小がこの比になる。

Ricci 曲率では、固定した一方向を含む二つの断面を足す。例えば $E_1$ では $K_{12}+K_{13}=4-2a^2$、$E_3$ では $2a^2$ である。曲率の混合成分が消えるので Ricci の非対角項も零になる。よって
$$
\begin{aligned}
\mathrm{Ric}&=\begin{pmatrix}\frac{5}{2} & 0 & 0\\0 & \frac{5}{2} & 0\\0 & 0 & \frac{3}{2}\end{pmatrix}\\
S&=\frac{13}{2}\\
\lambda_{\min}(\mathrm{Ric})&=\frac{3}{2}
\end{aligned}
$$
となる。計量は左移動で不変なので、全点で $\mathrm{Ric}\ge(3/2)g$ となる。このような下界と直径の関係は [曲率と位相](/math/riemannian-geometry/curvature-topology/) で扱う。

最小最大の比が $1/4$ より大きい。すべての点の曲率域が同じなので、点ごとの条件と大域ピンチの条件の両方を確認できる。
:::

:::example[三次元の左不変計量の縮約（尺度の二乗 3/2）]
単位四元数群の左不変計量で、虚四元数 $i,j,k$ の長さを $1,1,a$ とし、$a^2=\frac{3}{2}$ を選ぶ。正規直交標構を $(E_1,E_2,E_3)=(i,j,k/a)$ とする。ここで水平とは $E_1,E_2$ の張る面を指す。

正規直交な二方向に掛かる曲率は、$12$ 断面に $4-3a^2$、$13$ と $23$ の断面に $a^2$ である。曲率作用素はこの標構の二形式の基底で対角だから、任意の断面曲率はこれらの重み付き平均である。例えば $E_1$ と $Y=(E_2+E_3)/\sqrt2$ の断面では、二つの二形式の成分の二乗がともに $1/2$ なので、$K=(4-3a^2+a^2)/2=2-a^2=\frac{1}{2}$ になる。従って
$$
\begin{aligned}
(K_{12},K_{13},K_{23})&=\left( - \frac{1}{2}, \  \frac{3}{2}, \  \frac{3}{2}\right)\\
(K_{\min},K_{\max})&=\left( - \frac{1}{2}, \  \frac{3}{2}\right)\\
K_{\min}/K_{\max}&=- \frac{1}{3}
\end{aligned}
$$
である。絶対曲率を正規化するには計量を $K_{\max}$ 倍する。曲率はその逆数倍になるから、最大が $1$、最小がこの比になる。

Ricci 曲率では、固定した一方向を含む二つの断面を足す。例えば $E_1$ では $K_{12}+K_{13}=4-2a^2$、$E_3$ では $2a^2$ である。曲率の混合成分が消えるので Ricci の非対角項も零になる。よって
$$
\begin{aligned}
\mathrm{Ric}&=\begin{pmatrix}1 & 0 & 0\\0 & 1 & 0\\0 & 0 & 3\end{pmatrix}\\
S&=5\\
\lambda_{\min}(\mathrm{Ric})&=1
\end{aligned}
$$
となる。計量は左移動で不変なので、全点で $\mathrm{Ric}\ge g$ となる。

水平断面は負になるが、Ricci は全方向で正である。これは正の Ricci 曲率を持つ計量に負の断面も現れる例である。
:::

## Einstein 計量

:::definition[Einstein 計量]{#def-rg-einstein}
$\mathrm{Ric}=\lambda g$ を満たす計量を :term[Einstein 計量]{en="Einstein metric" yomi="あいんしゅたいんけいりょう"} という。以下 $\lambda$ は定数とする。次元が $3$ 以上で連結なら、初めは関数 $\lambda$ と仮定しても同じ定義になる。
:::

実際 $\nabla g=0$ により $\divg(\lambda g)=d\lambda$ であり、$S=n\lambda$ である。:ref[lem-rg-contracted-bianchi] から $(n-2)d\lambda=0$ となる。各点で全平面の曲率 $k(p)$ が等しい場合は $\lambda=(n-1)k$ を代入できる。これが $n\ge3$ における Schur の補題である。二次元では $\mathrm{Ric}=Kg$ が恒等的に成り立つのでこの議論で $K$ は一定とならない。

:::example[Einstein でも断面曲率が一定でない積]
$S^2(1)\times S^2(1)$ では $\mathrm{Ric}=g$ である。一つの因子に含まれる断面は $1$、異なる因子から各一本取った断面は $0$ である。一方、$S^2(1)\times S^2(2)$ では各因子の $\lambda$ が $1$ と $1/4$ なので積は Einstein でない。
:::

:::example[次元の異なる因子]
$S^2(a)\times S^3(b)$ の Ricci 固有値は各因子で $1/a^2$ と $2/b^2$ である。従って $b=\sqrt2 a$ のとき Einstein となる。例えば $a=1$、$b=\sqrt2$ で $\lambda=1$、$S=5$ である。第 1 因子内の断面は $1$、第 2 因子内は $1/2$、異なる因子を混ぜた断面は $0$ である。
:::

### 次元に依存する逆向きの主張

定断面曲率 $k$ なら Ricci は $(n-1)kg$ である。三次元での逆を直接確認するには、$\mathrm{Ric}=\lambda g$ を満たす正規直交基底を任意に取る。$K_{12}+K_{13}=\lambda$、$K_{12}+K_{23}=\lambda$、$K_{13}+K_{23}=\lambda$ から $K_{12}=K_{13}=K_{23}=\lambda/2$ である。基底の最初の二方向は任意の平面を張るように選べるので、全断面が同じ曲率を持つ。$\lambda$ は縮約 Bianchi 恒等式から一定である。二次元では $\mathrm{Ric}=Kg$ が自動的に成立するが、Einstein 計量の定義では係数を定数としている。従って次元 $2$ の Einstein 計量は一定 Gauss 曲率の計量である。

:::example[積の尺度をそろえる]
$S^3(2)\times S^5(b)$ の Ricci 固有値は $2/4=1/2$ と $4/b^2$ なので、$b=2\sqrt2$ の場合に Einstein となる。このとき次元は $8$ で $S=4$ となる。しかしそれぞれの因子の断面曲率は $1/4,1/8$ で混合断面は $0$ である。従ってスカラー曲率 $4$ から定曲率模型として $4/(8\cdot7)=1/14$ を当てると、その値はどの純粋断面にも一致しない。
:::

### 発散が零の曲率の組み合わせ

$G=\mathrm{Ric}-\frac12 Sg$ とおくと $\nabla^iG_{ij}=0$ である。実際、$\nabla g=0$ と :ref[lem-rg-contracted-bianchi] を使えば $\nabla^i(Sg_{ij})=\partial_jS$ が Ricci の発散の二倍になる。この対称形式は Einstein テンソルと呼ばれる。$n=2$ では $G$ は恒等的に零となる。$n\ne2$ で $G=cg$ とすると、発散の式から $dc=0$ であり、跡から $(1-n/2)S=nc$ が得られる。従って $S$ は一定で、元の式に戻せば Ricci も計量の定数倍である。ここでは Riemann 計量での縮約だけを使っている。

## Weyl テンソルと共形平坦性

四次元以上では、Ricci 曲率を記録しても全ての断面曲率は復元できない。縮約を含む部分を明示的に取り除く。本章では対称形式 $A$ と計量 $g$ の :term[Kulkarni–Nomizu 積]{en="Kulkarni–Nomizu product" yomi="かるかるにのみずせき"} を $A\mathbin{\diamond}g$ と書き、次式で定める。
$$
(A\mathbin{\diamond} g)_{ijk\ell}=A_{ik}g_{j\ell}+A_{j\ell}g_{ik}-A_{i\ell}g_{jk}-A_{jk}g_{i\ell}
$$
とする。対称性と巡回和を展開すれば曲率と同じ代数的対称性があり、最初と三つ目の縮約は $(n-2)A+(\tr_g A)g$ である。

:::definition[Weyl テンソル]{#def-rg-weyl}
$n\ge3$ とし $P=(\mathrm{Ric}-Sg/(2(n-1)))/(n-2)$ を :term[Schouten テンソル]{en="Schouten tensor" yomi="すはうてんてんそる"} とする。:term[Weyl テンソル]{en="Weyl tensor" yomi="わいるてんそる"} は $W=\mathcal R-P\mathbin{\diamond} g$ である。
:::

$P$ のトレースは $S/(2(n-1))$ なので $W$ の縮約は零になる。$n=3$ では任意の正規直交基底において、縮約の三つの対角成分が $W_{1212}, W_{1313}, W_{2323}$ の二項和を与える。この三つが零なら各成分は零である。非対角成分も、例えば縮約の $12$ 成分に寄与するのは $W_{3132}$ だけなので零である。添字の対称性で全成分を尽くす。

### 共形変化

$U=\grad u$、$H=\nabla du$ とする。:eqref[eq-rg-conformal-connection] の接続差 $A_XY=du(X)Y+du(Y)X-g(X,Y)U$ を曲率の定義に代入すると
$$
\widetilde R(X,Y)Z=R(X,Y)Z+(\nabla_XA)(Y,Z)-(\nabla_YA)(X,Z)+[A_X,A_Y]Z.
$$
$g$ で最後の引数と対にすると、微分項は $-(H\mathbin{\diamond} g)$、二次項は $(du\otimes du-\norm U^2g/2)\mathbin{\diamond} g$ になる。例えば微分項は $(\nabla_XA)(Y,Z)=H(X,Y)Z+H(X,Z)Y-g(Y,Z)\nabla_X U$ を反対称化すれば得られる。二次項では交換子 $[A_X,A_Y]Z$ を展開し、$du(Y)du(Z)X-du(X)du(Z)Y+g(Y,Z)du(X)U-g(X,Z)du(Y)U-\norm U^2(g(Y,Z)X-g(X,Z)Y)$ を得る。従って
$$
e^{-2u}\widetilde{\mathcal R}=\mathcal R+B\mathbin{\diamond} g,\qquad
B=-H+du\otimes du-\tfrac12\norm U^2g .
$$
縮約とその除去の一意性から $\widetilde P=P+B$、$\widetilde W=e^{2u}W$ である。

各点の近傍で正の関数を計量に掛けて平坦にできるとき、計量は :term[局所共形平坦]{en="locally conformally flat" yomi="きょくしょきょうへいへいたん"} であるという。

:::theorem[局所共形平坦性の判定]{#thm-rg-conformal-flatness}
$n\ge4$ では、局所共形平坦性と $W=0$ は同値である。$n=3$ では $W$ は常に零であり、局所共形平坦性は Cotton テンソル $C_{ijk}=\nabla_iP_{jk}-\nabla_jP_{ik}$ が零であることと同値である。
:::

:::proof
共形変化で平坦な計量が得られれば、$\widetilde W=e^{2u}W$ から $W=0$ である。まず $W=0$ の場合に $C$ にどの条件が加わるかを求める。縮約 Bianchi 恒等式と $\mathrm{Ric}=(n-2)P+(\tr_g P)g$ から $\nabla^aP_{ab}=\partial_b(\tr_gP)$ である。第二 Bianchi 恒等式を一回縮約すると
$$
\nabla^a\mathcal R_{bc a e}=\nabla_b\mathrm{Ric}_{ce}-\nabla_c\mathrm{Ric}_{be}.
$$
$\mathcal R=P\mathbin{\diamond} g$ を入れると、左辺は $C_{bce}+\partial_b(\tr_gP)g_{ce}-\partial_c(\tr_gP)g_{be}$ であり、右辺は同じ跡の項と $(n-2)C_{bce}$ の和である。従って $(n-3)C=0$ となる。

次に $u$ の代わりに 1 形式 $v$ を未知量とし、$\widetilde P=0$ を作るための方程式を
$$
\nabla_i v_j=Q_{ij},\qquad
Q_{ij}=P_{ij}+v_i v_j-\tfrac12|v|_g^2g_{ij}
$$
とする。右辺は対称なので、解があれば $dv=0$ となり、局所的に $v=du$ と書ける。この方程式の可積分性を調べるには、右辺を使って $\nabla_i Q_{jk}-\nabla_j Q_{ik}$ を展開する。得られる式は
$$
C_{ijk}+v_jQ_{ik}-v_iQ_{jk}-v^aQ_{ia}g_{jk}+v^aQ_{ja}g_{ik}
=C_{ijk}-(P\mathbin{\diamond} g)_{ijak}v^a
$$
である。最後の等式では、三次の項は $v_jv_iv_k-v_iv_jv_k=0$ と $v^aQ_{ia}=v^aP_{ia}+\frac12|v|^2v_i$ により相殺する。一方、1 形式の共変微分の交換子は $-\mathcal R_{ijak}v^a$ である。従って可積分条件は $C_{ijk}+W_{ijak}v^a=0$ となる。特に三次元で共形平坦なら、この条件を平坦化する $u$ に適用して $C=0$ を得る。

逆に $W=C=0$ とする。$T^*M$ の局所座標 $(x^i,v_j)$ で
$$
D_i=\partial_{x^i}+(\Gamma^a_{ij}v_a+Q_{ij})\partial_{v_j}
$$
を考える。$[D_i,D_j]$ は垂直で、その成分は上の可積分条件の左辺である。実際、点で $\Gamma=0$ の座標を取ると、$\Gamma$ を微分した項が $\mathcal R_{ijak}v^a$、$Q$ の全微分の反対称部分が上の展開式となる。従って分布は包合的であり、:ref[thm-mfd-frobenius]{text="Frobenius の定理"} により局所積分多様体を持つ。$D_i$ の $x$ 成分は座標基底なので、積分多様体は局所的にある 1 形式 $v$ のグラフである。それは $\nabla_i v_j=Q_{ij}$ を満たす。$v$ は閉形式なので局所的に $du$ と書け、$\widetilde P=0$ である。$\widetilde W=e^{2u}W=0$ も成立し、従って $\widetilde{\mathcal R}=0$ となる。:ref[thm-rg-flat] により $\widetilde g$ は局所 Euclid 計量である。$n\ge4$ では $W=0$ から $C=0$ をすでに得た。$n=3$ では $W=0$ が恒等式であるため、必要な仮定は $C=0$ だけである。
:::

二次元では、共形変化の式の $B$ を縮約すると $\widetilde K=e^{-2u}(K-\Delta_g u)$ となる。勾配の二次項のトレースが $|\grad u|^2-2|\grad u|^2/2=0$ で消えるためである。小座標球で $\Delta_g u=K$ の Dirichlet 問題を解き、弱解の存在には [変分法](/math/pde/variational/) の Dirichlet 問題の理論を、内部の滑らかさには [楕円型正則性](/math/pde/elliptic-regularity/) を用いる。そこで $\widetilde g$ は平坦となるので、どの二次元計量も局所共形平坦である。

:::example[Weyl 成分の検出]
$S^2(1)\times S^2(1)$ は $n=4$、$\mathrm{Ric}=g$、$S=4$ なので $P=g/6$ である。混合平面では $W_{ijij}=0-1/3=-1/3$ となる。従ってこの計量は局所共形平坦でない。
:::

:::example[共形変化しても残る成分]
$S^2(1)\times S^2(1)$ で $g$ を $e^{2u}g$ に変える。元の単位混合対を $X,Y$ とすると、新計量で単位の対は $e^{-u}X,e^{-u}Y$ である。共形変換則からこの対での成分は $-e^{-2u}/3$ となる。例えば $u=\log2$ の点では $-1/12$ であり、定数倍で曲率成分が $1/4$ 倍に変わることとも一致する。
:::

:::example[三次元で Cotton テンソルが消えない積]
二次元計量 $h=e^{2x^2}(dx^2+dy^2)$ と実線の積を取る。二次元の共形曲率公式から $K_h=-2e^{-2x^2}$ である。積では $\mathrm{Ric}=K_hh$、$S=2K_h$ なので、三次元の Schouten テンソルは $P=(K_h/2)(h-dz^2)$ である。積接続は $h,dz$ を平行に保つため、$C_{xzz}=-\partial_x K_h/2=-4xe^{-2x^2}$ となる。従って $x\ne0$ の点では、Weyl テンソルが零であっても局所共形平坦でない。
:::

### 二つの定曲率因子の Weyl 曲率

$p,q\ge2$ とし、断面曲率 $k_1,k_2$ の $p,q$ 次元の因子の積を取る。$n=p+q$、$S=p(p-1)k_1+q(q-1)k_2$ とする。Schouten テンソルの第 1、第 2 因子の固有値はそれぞれ
$$
a=\frac{(p-1)k_1-S/(2(n-1))}{n-2},\qquad
b=\frac{(q-1)k_2-S/(2(n-1))}{n-2}.
$$
従って純粋な二断面と混合断面での Weyl 成分は
$$
\begin{aligned}
W_{ijij}^{(1)}&=k_1-2a
=\frac{q(q-1)(k_1+k_2)}{(n-1)(n-2)},\\
W_{\alpha\beta\alpha\beta}^{(2)}&=k_2-2b
=\frac{p(p-1)(k_1+k_2)}{(n-1)(n-2)},\\
W_{i\alpha i\alpha}&=-a-b
=-\frac{(p-1)(q-1)(k_1+k_2)}{(n-1)(n-2)}.
\end{aligned}
$$
例えば最初の式では $k_1(n-1)(n-2)-2(p-1)k_1(n-1)+S$ を展開すると、$k_1$ の係数が $q(q-1)$ になり、$k_2$ の係数も $q(q-1)$ である。残る式は同じ分母で分子を整理すればよい。因子に沿う正規直交基底では $\mathscr R$ も $P$ も対角表示であり、Weyl の曲率作用素も対角だから、以上が全固有値を尽くす。従ってこの積が局所共形平坦である条件は $k_1+k_2=0$ である。

:::example[正負の曲率が Weyl 成分で打ち消し合う]
$S^2(1)\times H^2(1)$ では $k_1=1,k_2=-1$ なので $S=0$ である。$P$ は正の因子で $g/2$、負の因子で $-g/2$ を持つ。$P\mathbin{\diamond}g$ はそれぞれの純粋断面で $1,-1$、混合断面で $0$ となり、$\mathcal R$ と一致する。従って $W=0$ で局所共形平坦であるが、曲率自体は零でない。この例では、Schouten テンソルが残す曲率は共形因子によって局所的に取り除ける。
:::

### Ricci とスカラー曲率の共形変換

先に得た $\widetilde P=P+B$ と $\widetilde g=e^{2u}g$ を縮約する。$\tr_g B=-\Delta_g u-(n-2)\norm{\grad u}^2/2$ なので、
$$
\begin{aligned}
\widetilde{\mathrm{Ric}}
&=\mathrm{Ric}-(n-2)(\operatorname{Hess}u-du\otimes du)
-\bigl(\Delta_g u+(n-2)\norm{\grad u}^2\bigr)g,\\
\widetilde S
&=e^{-2u}\bigl(S-2(n-1)\Delta_g u
-(n-1)(n-2)\norm{\grad u}^2\bigr).
\end{aligned}
$$
実際、曲率の最初と三つ目を縮約すると $B\mathbin{\diamond}g$ は $(n-2)B+(\tr_g B)g$ を与える。$\widetilde{\mathcal R}$ の倍率 $e^{2u}$ と逆計量の倍率 $e^{-2u}$ は Ricci の計算で相殺する。最後にもう一度 $\widetilde g$ で縮約する分だけ、スカラー曲率には $e^{-2u}$ が残る。例えば $g$ を三次元 Euclid 計量、$u=-\log y$ とすると、$\Delta_g u=\norm{\grad u}^2=1/y^2$ だから $\widetilde S=-4-2=-6$ であり、全断面曲率 $-1$ の縮約と一致する。

:::example[線形の共形指数の二次元と三次元の比較]
$u=x$ として $e^{2u}g_{\rm E}$ を考える。三次元で扱う場合は $u$ は第 3 座標 $z$ に依存しないとする。指数の勾配が一定で、二階微分は零である。二次元と三次元で曲率がどう変わるかを比べる。

まず $du$ の接座標成分は $\begin{pmatrix}1\\0\end{pmatrix}$、Euclid Laplacian は $0$ である。二次元の共形 Gauss 曲率では勾配二次項が相殺するので
$$
K_2=-e^{-2u}\Delta_{\rm E}u=0 .
$$
これに対し三次元では最初のスカラー曲率が零、$n-1=2$、$n-2=1$ だから
$$
S_3=e^{-2u}(-4\Delta_{\rm E}u-2|du|_{\rm E}^2)
=- 2 e^{- 2 x} .
$$
ここで $|du|_{\rm E}^2=1$ とした。一つの関数に同じ計量倍率を使っても、跡を取る空間の次元によりどの項が残るかが変わる。

点 $\left( 1, \  0\right)$ では、$u=1$、$K_2=0$、$S_3=- \frac{2}{e^{2}}$ である。2 次元では $K_2$ が恒等的に零なら局所平坦性が従う。3 次元の対応する計量で $S_3$ が一点でも非零なら、その点の曲率テンソルは非零である。いずれの次元でも出発点が平坦なので共形平坦性は定義から成立する。従って共形平坦性、スカラー曲率の消失、曲率テンソルの消失を判定する式はそれぞれ異なる。
:::

:::example[原点で対称な二次式の二次元と三次元の比較]
$u=x^{2} + y^{2}$ として $e^{2u}g_{\rm E}$ を考える。三次元で扱う場合は $u$ は第 3 座標 $z$ に依存しないとする。原点では一階微分が零であるが、二階微分は零でない。点での記号と曲率を区別する例である。

まず $du$ の接座標成分は $\begin{pmatrix}2 x\\2 y\end{pmatrix}$、Euclid Laplacian は $4$ である。二次元の共形 Gauss 曲率では勾配二次項が相殺するので
$$
K_2=-e^{-2u}\Delta_{\rm E}u=- 4 e^{- 2 x^{2} - 2 y^{2}} .
$$
これに対し三次元では最初のスカラー曲率が零、$n-1=2$、$n-2=1$ だから
$$
S_3=e^{-2u}(-4\Delta_{\rm E}u-2|du|_{\rm E}^2)
=8 \left(- x^{2} - y^{2} - 2\right) e^{- 2 x^{2} - 2 y^{2}} .
$$
ここで $|du|_{\rm E}^2=4 x^{2} + 4 y^{2}$ とした。一つの関数に同じ計量倍率を使っても、跡を取る空間の次元によりどの項が残るかが変わる。

点 $\left( 0, \  0\right)$ では、$u=0$、$K_2=-4$、$S_3=-16$ である。2 次元では $K_2$ が恒等的に零なら局所平坦性が従う。3 次元の対応する計量で $S_3$ が一点でも非零なら、その点の曲率テンソルは非零である。いずれの次元でも出発点が平坦なので共形平坦性は定義から成立する。従って共形平坦性、スカラー曲率の消失、曲率テンソルの消失を判定する式はそれぞれ異なる。
:::

:::example[調和な二次式の二次元と三次元の比較]
$u=x^{2} - y^{2}$ として $e^{2u}g_{\rm E}$ を考える。三次元で扱う場合は $u$ は第 3 座標 $z$ に依存しないとする。二方向の二階微分は反対符号なので平面 Laplacian が零になる。計量の係数は定数でない。

まず $du$ の接座標成分は $\begin{pmatrix}2 x\\- 2 y\end{pmatrix}$、Euclid Laplacian は $0$ である。二次元の共形 Gauss 曲率では勾配二次項が相殺するので
$$
K_2=-e^{-2u}\Delta_{\rm E}u=0 .
$$
これに対し三次元では最初のスカラー曲率が零、$n-1=2$、$n-2=1$ だから
$$
S_3=e^{-2u}(-4\Delta_{\rm E}u-2|du|_{\rm E}^2)
=8 \left(- x^{2} - y^{2}\right) e^{- 2 x^{2} + 2 y^{2}} .
$$
ここで $|du|_{\rm E}^2=4 x^{2} + 4 y^{2}$ とした。一つの関数に同じ計量倍率を使っても、跡を取る空間の次元によりどの項が残るかが変わる。

点 $\left( 1, \  1\right)$ では、$u=0$、$K_2=0$、$S_3=-16$ である。2 次元では $K_2$ が恒等的に零なら局所平坦性が従う。3 次元の対応する計量で $S_3$ が一点でも非零なら、その点の曲率テンソルは非零である。いずれの次元でも出発点が平坦なので共形平坦性は定義から成立する。従って共形平坦性、スカラー曲率の消失、曲率テンソルの消失を判定する式はそれぞれ異なる。
:::

:::example[混合二次式の二次元と三次元の比較]
$u=x y$ として $e^{2u}g_{\rm E}$ を考える。三次元で扱う場合は $u$ は第 3 座標 $z$ に依存しないとする。対角の二階微分は零だが、混合二階微分は零でない。Hessian の跡の情報と行列の情報を分ける。

まず $du$ の接座標成分は $\begin{pmatrix}y\\x\end{pmatrix}$、Euclid Laplacian は $0$ である。二次元の共形 Gauss 曲率では勾配二次項が相殺するので
$$
K_2=-e^{-2u}\Delta_{\rm E}u=0 .
$$
これに対し三次元では最初のスカラー曲率が零、$n-1=2$、$n-2=1$ だから
$$
S_3=e^{-2u}(-4\Delta_{\rm E}u-2|du|_{\rm E}^2)
=2 \left(- x^{2} - y^{2}\right) e^{- 2 x y} .
$$
ここで $|du|_{\rm E}^2=x^{2} + y^{2}$ とした。一つの関数に同じ計量倍率を使っても、跡を取る空間の次元によりどの項が残るかが変わる。

点 $\left( 1, \  2\right)$ では、$u=2$、$K_2=0$、$S_3=- \frac{10}{e^{4}}$ である。2 次元では $K_2$ が恒等的に零なら局所平坦性が従う。3 次元の対応する計量で $S_3$ が一点でも非零なら、その点の曲率テンソルは非零である。いずれの次元でも出発点が平坦なので共形平坦性は定義から成立する。従って共形平坦性、スカラー曲率の消失、曲率テンソルの消失を判定する式はそれぞれ異なる。
:::

### 跡を除いた成分の次元

$n\ge3$ では、代数的曲率テンソルの空間から対称 2 形式への Ricci 縮約は全射である。任意の対称形式 $B$ に対して
$$
A=\frac1{n-2}\left(B-\frac{\tr_gB}{2(n-1)}g\right)
$$
を取ると、$A\mathbin{\diamond}g$ の縮約は $(n-2)A+(\tr_g A)g=B$ であるからである。従って Ricci の補空間、すなわち Weyl の空間の次元は
$$
\frac{n^2(n^2-1)}{12}-\frac{n(n+1)}2
=\frac{n(n+1)(n+2)(n-3)}{12}
$$
となる。$n=3$ では零、$n=4$ では $10$、$n=5$ では $35$ である。例えば四次元では曲率全体が $20$ 成分、対称 Ricci が $10$ 成分なので、残り十次元の違いを Weyl が取り出している。この数え方は一点の代数的自由度であり、曲率の微分が満たす第二 Bianchi 恒等式の制約は含めていない。

## 曲率の座標計算と例

座標で $\nabla_i(\Gamma^m_{j\ell}\partial_m)$ を展開し $i,j$ を交換した差を取ると
$$
R^k{}_{\ell ij}=\partial_i\Gamma^k_{j\ell}-\partial_j\Gamma^k_{i\ell}
+\Gamma^k_{im}\Gamma^m_{j\ell}-\Gamma^k_{jm}\Gamma^m_{i\ell}
\label{eq-rg-curvature-coordinate}
$$
を得る。Ricci 成分は $R^i{}_{\ell i j}$ であり、ここでは $i,j$ の曲率方向と $\ell$ の作用先を区別する。

非座標標構では接続の行列 $1$ 形式 $\omega$ を用い、$\Omega(X,Y)$ を同じ標構での $R(X,Y)$ の行列とする。成分列 $v$ に対して $\nabla_X$ は $Xv+\omega(X)v$ と作用する。交換子から $[X,Y]v+\omega([X,Y])v$ を引くと、$v$ の微分項が相殺し、
$$
\Omega=d\omega+\omega\wedge\omega
$$
となる。ここで $(\omega\wedge\omega)(X,Y)=[\omega(X),\omega(Y)]$ である。これは曲率に関する :term[Cartan の構造方程式]{en="Cartan's structure equation" yomi="かるたんのこうぞうほうていしき"} である。各成分を外微分すると $d\Omega=d\omega\wedge\omega-\omega\wedge d\omega$ となる。$\Omega$ の式を代入すれば $d\Omega+\omega\wedge\Omega-\Omega\wedge\omega=0$ である。これを三つの座標場で評価すると、$\Omega$ の出力と作用先に対する接続差が行列の交換子として現れる。方向を表す二つの引数に対する接続項は、巡回和では捩率零により相殺する。従ってこの式は第二 Bianchi 恒等式の標構表示である。

:::example[回転対称計量]
$dr^2+f^2d\theta^2$ では $R^r{}_{\theta r\theta}=\partial_r(-ff')-(-ff')(f'/f)=-ff''$ である。断面曲率は $K=-f''/f$ となる。$f=a\sin(r/a)$、$r$、$a\sinh(r/a)$ はそれぞれ $1/a^2,0,-1/a^2$ を与える。$f=2+\cos r$ の場合は $K=\cos r/(2+\cos r)$ で、$r=0,\pi$ で $1/3,-1$ となる。
:::

:::example[二次元の共形座標]
$g=e^{2u}(dx^2+dy^2)$ では $K=-e^{-2u}(u_{xx}+u_{yy})$ となる。実際、平坦計量との共形変化の $B\mathbin{\diamond} g_{\mathrm E}$ の $1212$ 成分は $\tr B=-\Delta_{\mathrm E}u$ であり、$\mathcal R_{1212}$ の倍率は $e^{2u}$、Gram 行列式は $e^{4u}$ である。$u=-\log y$ では $u_{yy}=1/y^2$ なので $K=-1$ である。$u=x$ では $K=0$ であり、局所写像 $(x,y)\mapsto(e^x\cos y,e^x\sin y)$ も平坦性を示す。
:::

定数倍 $\widetilde g=a^2g$ では $\widetilde R=R$、$\widetilde{\mathrm{Ric}}=\mathrm{Ric}$ である。後者は $Z\mapsto R(Z,X)Y$ のトレースで示せる。断面曲率とスカラー曲率は逆計量または Gram 行列式の倍率により $a^{-2}$ 倍になる。

### 正規座標の二次の係数を復元する

正規座標 $x$ では $g_{ij}(x)x^j=x^i$ である。固定した $v$ に対し、この式の三次の項を取ると $g_{ij,ab}(0)v^jv^av^b=0$ である。これと曲率の式から
$$
g_{ij}(x)=\delta_{ij}-\tfrac13\mathcal R_{ia jb}(p)x^ax^b+O(|x|^3)
$$
を得る。復元の手順を示す。二次部分を $h_{ij}(x)$ とし $h_{ij}(x)x^j=0$ と置く。任意の $x$ で $h_{ij,ab}x^ax^b=2h_{ij}$ である。$h_{ib}x^b=0$ をまず $x^a$、次いで $x^j$ で微分し、微分前の $h$ が二次同次であることを使うと
$$
h_{ib,aj}x^ax^b=-h_{ij},\qquad
h_{ab,ij}x^ax^b=2h_{ij}
$$
となる。第 1 式は一階の条件 $h_{ib,j}x^b=-h_{ij}$ に動径微分 $x^a\partial_a$ を作用させて得る。第 2 式は $h_{ab}x^ax^b=0$ を二階微分し、交差項を同じ一階の条件で整理する。曲率式
$$
2\mathcal R_{ia jb}=
h_{aj,ib}-h_{ab,ij}-h_{ij,ab}+h_{ib,aj}
$$
に $x^ax^b$ を掛けると、右辺は $-h_{ij}-2h_{ij}-2h_{ij}-h_{ij}=-6h_{ij}$ になる。これが二次の係数である。

二次の行列式展開 $\det(I+H)=1+\tr H+O(\norm H^2)$ は、行列式の積の和で $H$ を一つだけ取る項が対角の和になることから従う。平方根を取れば、正規座標での体積密度は
$$
\sqrt{\det G(x)}=1-\tfrac16\mathrm{Ric}_{ab}(p)x^ax^b+O(|x|^3)
$$
となる。半径 $r$ の小正規球で座標積分すると、非対角積分は反射で消え、対角積分はそれぞれ $\omega_n r^{n+2}/(n+2)$ である。実際、和は極座標で $n\omega_n\int_0^r t^{n+1}dt$ で、各成分は等しい。従って
$$
\operatorname{Vol}(B_p(r))=\omega_n r^n
\left(1-\frac{S(p)}{6(n+2)}r^2+O(r^4)\right).
$$
ここで $\omega_n$ は Euclid 単位球の体積である。余りが $O(r^4)$ になるのは、滑らかな体積密度を三次まで展開したとき、三次の斉次項が原点反転で符号を変えて接球上の積分から消えるためである。この公式は $S$ が小さい距離球の最初の体積補正として現れることを示す。

:::example[球面上の小円板の面積]
半径 $2$ の球面では $K=1/4$、$S=1/2$ だから、公式は $\pi r^2(1-r^2/48+O(r^4))$ となる。一方、極座標の体積を積分すると $8\pi(1-\cos(r/2))$ であり、Taylor 展開の最初の二項は $\pi r^2-\pi r^4/48$ である。曲率の正符号が面積を減らす符号として現れる。
:::

:::example[円錐の頂点を除いた領域の座標曲率と断面曲率]
$g=dr^2+f^2d\theta^2$、$f(r)=2 r$ とし、$r=1$ での曲率を求める。ここでは $r>0$ とする。角を $2\theta$ に替えると局所 Euclid 表示を得るが、角の周期は二倍になる。

記号 $\Gamma^r_{\theta\theta}=-ff'$ を $r$ で微分すると $-(f')^2-ff''$ である。一方、曲率式で引かれる積 $\Gamma^r_{\theta\theta}\Gamma^\theta_{r\theta}$ は $-(f')^2$ である。従って二乗の項が相殺し、
$$
\begin{aligned}
R^r_{\theta r\theta}&=-ff''=0\\
g(R(\partial_r,\partial_\theta)\partial_\theta,\partial_r)&=R^r_{\theta r\theta}\\
K&=\frac{R^r_{\theta r\theta}}{f^2}=0
\end{aligned}
$$
となる。指定点で $f=2$、$f''=0$ なので、曲率の座標成分は $0$、断面曲率は $0$ である。座標断面の分母は $g_{rr}g_{\theta\theta}-g_{r\theta}^2=f^2$ だから、これを割った値が長さをそろえた断面の曲率になる。

正規直交標構 $E_1=\partial_r,E_2=f^{-1}\partial_\theta$ では $R(E_1,E_2)E_2=KE_1$ であり、最後の二引数での交代性から $R(E_1,E_2)E_1=-KE_2$ である。従って作用素の行列は
$$
R(E_1,E_2)=\begin{pmatrix}0 & 0\\0 & 0\end{pmatrix}.
$$
二次元の $\mathrm{Ric}$ は両単位方向で同じ断面を一回ずつ数えるので $Kg$ であり、スカラー曲率は $2K=0$ となる。こうして同じ点での座標成分、断面曲率、スカラー曲率を、どの標構で何回縮約するかによって区別できる。
:::

:::example[カテノイドのくびれの外側の座標曲率と断面曲率]
$g=dr^2+f^2d\theta^2$、$f(r)=\sqrt{r^{2} + 1}$ とし、$r=1$ での曲率を求める。半径は $r=0$ で最小になる。左右の同じ半径の円では $f'$ の符号が異なる。

記号 $\Gamma^r_{\theta\theta}=-ff'$ を $r$ で微分すると $-(f')^2-ff''$ である。一方、曲率式で引かれる積 $\Gamma^r_{\theta\theta}\Gamma^\theta_{r\theta}$ は $-(f')^2$ である。従って二乗の項が相殺し、
$$
\begin{aligned}
R^r_{\theta r\theta}&=-ff''=- \frac{1}{r^{2} + 1}\\
g(R(\partial_r,\partial_\theta)\partial_\theta,\partial_r)&=R^r_{\theta r\theta}\\
K&=\frac{R^r_{\theta r\theta}}{f^2}=- \frac{1}{\left(r^{2} + 1\right)^{2}}
\end{aligned}
$$
となる。指定点で $f=\sqrt{2}$、$f''=\frac{\sqrt{2}}{4}$ なので、曲率の座標成分は $- \frac{1}{2}$、断面曲率は $- \frac{1}{4}$ である。座標断面の分母は $g_{rr}g_{\theta\theta}-g_{r\theta}^2=f^2$ だから、これを割った値が長さをそろえた断面の曲率になる。

正規直交標構 $E_1=\partial_r,E_2=f^{-1}\partial_\theta$ では $R(E_1,E_2)E_2=KE_1$ であり、最後の二引数での交代性から $R(E_1,E_2)E_1=-KE_2$ である。従って作用素の行列は
$$
R(E_1,E_2)=\begin{pmatrix}0 & - \frac{1}{4}\\\frac{1}{4} & 0\end{pmatrix}.
$$
二次元の $\mathrm{Ric}$ は両単位方向で同じ断面を一回ずつ数えるので $Kg$ であり、スカラー曲率は $2K=- \frac{1}{2}$ となる。こうして同じ点での座標成分、断面曲率、スカラー曲率を、どの標構で何回縮約するかによって区別できる。
:::

:::example[双曲円柱の軸から離れた位置の座標曲率と断面曲率]
$g=dr^2+f^2d\theta^2$、$f(r)=\cosh{\left(r \right)}$ とし、$r=\log{\left(2 \right)}$ での曲率を求める。座標 $r$ は閉測地線の軸からの符号付き距離である。角方向を周期で割っても局所計算は同じである。

記号 $\Gamma^r_{\theta\theta}=-ff'$ を $r$ で微分すると $-(f')^2-ff''$ である。一方、曲率式で引かれる積 $\Gamma^r_{\theta\theta}\Gamma^\theta_{r\theta}$ は $-(f')^2$ である。従って二乗の項が相殺し、
$$
\begin{aligned}
R^r_{\theta r\theta}&=-ff''=- \cosh^{2}{\left(r \right)}\\
g(R(\partial_r,\partial_\theta)\partial_\theta,\partial_r)&=R^r_{\theta r\theta}\\
K&=\frac{R^r_{\theta r\theta}}{f^2}=-1
\end{aligned}
$$
となる。指定点で $f=\frac{5}{4}$、$f''=\frac{5}{4}$ なので、曲率の座標成分は $- \frac{25}{16}$、断面曲率は $-1$ である。座標断面の分母は $g_{rr}g_{\theta\theta}-g_{r\theta}^2=f^2$ だから、これを割った値が長さをそろえた断面の曲率になる。

正規直交標構 $E_1=\partial_r,E_2=f^{-1}\partial_\theta$ では $R(E_1,E_2)E_2=KE_1$ であり、最後の二引数での交代性から $R(E_1,E_2)E_1=-KE_2$ である。従って作用素の行列は
$$
R(E_1,E_2)=\begin{pmatrix}0 & -1\\1 & 0\end{pmatrix}.
$$
二次元の $\mathrm{Ric}$ は両単位方向で同じ断面を一回ずつ数えるので $Kg$ であり、スカラー曲率は $2K=-2$ となる。こうして同じ点での座標成分、断面曲率、スカラー曲率を、どの標構で何回縮約するかによって区別できる。
:::

:::example[回転トーラスの内側の赤道の座標曲率と断面曲率]
$g=dr^2+f^2d\theta^2$、$f(r)=\cos{\left(r \right)} + 2$ とし、$r=\pi$ での曲率を求める。内側の赤道では角方向の半径が最小になる。外側の赤道 $r=0$ と比べると二階微分の符号が逆である。

記号 $\Gamma^r_{\theta\theta}=-ff'$ を $r$ で微分すると $-(f')^2-ff''$ である。一方、曲率式で引かれる積 $\Gamma^r_{\theta\theta}\Gamma^\theta_{r\theta}$ は $-(f')^2$ である。従って二乗の項が相殺し、
$$
\begin{aligned}
R^r_{\theta r\theta}&=-ff''=\left(\cos{\left(r \right)} + 2\right) \cos{\left(r \right)}\\
g(R(\partial_r,\partial_\theta)\partial_\theta,\partial_r)&=R^r_{\theta r\theta}\\
K&=\frac{R^r_{\theta r\theta}}{f^2}=\frac{\cos{\left(r \right)}}{\cos{\left(r \right)} + 2}
\end{aligned}
$$
となる。指定点で $f=1$、$f''=1$ なので、曲率の座標成分は $-1$、断面曲率は $-1$ である。座標断面の分母は $g_{rr}g_{\theta\theta}-g_{r\theta}^2=f^2$ だから、これを割った値が長さをそろえた断面の曲率になる。

正規直交標構 $E_1=\partial_r,E_2=f^{-1}\partial_\theta$ では $R(E_1,E_2)E_2=KE_1$ であり、最後の二引数での交代性から $R(E_1,E_2)E_1=-KE_2$ である。従って作用素の行列は
$$
R(E_1,E_2)=\begin{pmatrix}0 & -1\\1 & 0\end{pmatrix}.
$$
二次元の $\mathrm{Ric}$ は両単位方向で同じ断面を一回ずつ数えるので $Kg$ であり、スカラー曲率は $2K=-2$ となる。こうして同じ点での座標成分、断面曲率、スカラー曲率を、どの標構で何回縮約するかによって区別できる。
:::

:::example[指数的に縮む円柱の端の座標曲率と断面曲率]
$g=dr^2+f^2d\theta^2$、$f(r)=e^{- \frac{r}{2}}$ とし、$r=2 \log{\left(2 \right)}$ での曲率を求める。底の弧長を増やしたとき円周が一定の割合で縮む。座標ベクトルを使うと角方向の基底の長さが場所によって異なる。

記号 $\Gamma^r_{\theta\theta}=-ff'$ を $r$ で微分すると $-(f')^2-ff''$ である。一方、曲率式で引かれる積 $\Gamma^r_{\theta\theta}\Gamma^\theta_{r\theta}$ は $-(f')^2$ である。従って二乗の項が相殺し、
$$
\begin{aligned}
R^r_{\theta r\theta}&=-ff''=- \frac{e^{- r}}{4}\\
g(R(\partial_r,\partial_\theta)\partial_\theta,\partial_r)&=R^r_{\theta r\theta}\\
K&=\frac{R^r_{\theta r\theta}}{f^2}=- \frac{1}{4}
\end{aligned}
$$
となる。指定点で $f=\frac{1}{2}$、$f''=\frac{1}{8}$ なので、曲率の座標成分は $- \frac{1}{16}$、断面曲率は $- \frac{1}{4}$ である。座標断面の分母は $g_{rr}g_{\theta\theta}-g_{r\theta}^2=f^2$ だから、これを割った値が長さをそろえた断面の曲率になる。

正規直交標構 $E_1=\partial_r,E_2=f^{-1}\partial_\theta$ では $R(E_1,E_2)E_2=KE_1$ であり、最後の二引数での交代性から $R(E_1,E_2)E_1=-KE_2$ である。従って作用素の行列は
$$
R(E_1,E_2)=\begin{pmatrix}0 & - \frac{1}{4}\\\frac{1}{4} & 0\end{pmatrix}.
$$
二次元の $\mathrm{Ric}$ は両単位方向で同じ断面を一回ずつ数えるので $Kg$ であり、スカラー曲率は $2K=- \frac{1}{2}$ となる。こうして同じ点での座標成分、断面曲率、スカラー曲率を、どの標構で何回縮約するかによって区別できる。
:::

### 一次元を底とする歪積の全曲率

$g=dr^2+f(r)^2h$ を考える。$h$ は $m$ 次元多様体 $N$ の計量、$f>0$ である。$N$ の場 $U,V,W$ は $r$ に依存しないよう延長し、$T=\partial_r$ と書く。全ての $[T,U]$ は零である。Koszul の公式で $T$ を一つ含む成分を計算すると
$$
\nabla_TT=0,\qquad
\nabla_TU=\nabla_UT=(f'/f)U,\qquad
\nabla_U V=\nabla^h_U V-ff'h(U,V)T
$$
である。例えば最後の式の動径成分は $-Tg(U,V)/2=-ff'h(U,V)$ であり、接成分では Koszul の各項に共通の倍率 $f^2$ が掛かるので $\nabla^h$ が残る。$\nabla_TU$ の角方向には $Tg(U,V)/2$ が入り、$g=(f^2h)$ で戻すと $f'/f$ である。

これを曲率の定義に使う。$c=f'/f$ とすると $R(U,T)T=-\nabla_T(cU)=-(c'+c^2)U=-(f''/f)U$ である。また
$$
R(T,U)V=-ff''h(U,V)T
$$
となる。実際、$\nabla_T\nabla_UV$ の角成分は $c\nabla_U^hV$、動径成分は $-(ff')'h(U,V)T$ である。$\nabla_U\nabla_T V=c\nabla_U V$ を引くと角成分が消え、動径では $-(ff')'+c ff'=-ff''$ が残る。$R(U,V)T=c(\nabla_UV-\nabla_VU-[U,V])=0$ である。角方向だけの曲率では
$$
R(U,V)W=R^h(U,V)W-(f')^2\bigl(h(V,W)U-h(U,W)V\bigr).
$$
第一項は三つの共変微分を全て角方向に取った部分である。途中の微分が $-ff'h(V,W)T$ を作ると、その次の微分の角成分が $-ff'(f'/f)h(V,W)U$ になり、もう一つの項はこれの $U,V$ の交換である。法成分は $g(R(U,V)W,T)=-g(W,R(U,V)T)=0$ である。

従って動径平面と純粋角平面の断面曲率は
$$
K(T,U)=-f''/f,\qquad
K(U,V)=\frac{K_h(U,V)-(f')^2}{f^2}
$$
である。第二式では $h$ 単位な場を $g$ 単位にすると $f$ で割ることを使った。例えば $h$ が単位球面計量で $f(r)=r+r^3$、$r>0$ とすれば、$r=1$ で $f=2,f'=4,f''=6$ となる。動径平面は $-3$、角平面は $(1-16)/4=-15/4$ である。$m=2$ の場合、正規直交基底で Ricci は動径方向が $-6$、二つの角方向はそれぞれ $-3-15/4=-27/4$、スカラー曲率は $-39/2$ である。同じ点でも、動径の Jacobi 方程式が使う曲率だけからスカラー曲率を求めるには、残りの角平面も数えなければならない。

## 演習問題

::::exercise{level=basic}
$dr^2+\sinh^2r\,d\theta^2$ の曲率を求めよ。

:::solution
$f''=f$ なので $-f''/f=-1$ である。二次元では $\mathrm{Ric}=-g$、スカラー曲率は $-2$ である。
:::
::::

::::exercise{level=basic}
四次元定曲率 $2$ の空間の Ricci 曲率とスカラー曲率を求めよ。

:::solution
単位方向あたり直交する $3$ 個の断面曲率を足すので $\mathrm{Ric}=6g$ である。さらに四方向を足して $S=24$ となる。
:::
::::

::::exercise{level=basic}
$K=3$ の計量を $9$ 倍すると $K$ はいくつか。

:::solution
接続が同じで分子が $9$ 倍、分母が $81$ 倍なので $3/9=1/3$ である。
:::
::::

::::exercise{level=basic}
$e^{2(x^2+y^2)}(dx^2+dy^2)$ の原点での曲率を求めよ。

:::solution
$u_{xx}=u_{yy}=2$、$u(0)=0$ なので $K(0)=-4$ である。
:::
::::

::::exercise{level=advanced}
$dr^2+f(r)^2d\theta^2$ において $f(0)=0,f'(0)=1$、$r>0$ で $K=4$ のときの $f$ を求めよ。

:::solution
$f''+4f=0$ である。一般解 $A\cos(2r)+B\sin(2r)$ に初期値を入れると $A=0$、$2B=1$ なので $f(r)=\frac12\sin(2r)$ である。正の範囲では例えば $0<r<\pi/2$ をとる。
:::
::::

::::exercise{level=advanced}
二次元では $S$ から曲率テンソルが決まることを示せ。

:::solution
すべての二次元断面は接空間自身であり、$S=2K$ である。従って :ref[lem-rg-sectional-determines] によって $R(X,Y)Z=\frac S2(g(Y,Z)X-g(X,Z)Y)$ と決まる。
:::
::::

::::exercise{level=advanced}
$S^2(1)\times S^2(2)$ のスカラー曲率と二種類の純粋断面の Weyl 成分を求めよ。

:::solution
Ricci 固有値は $1,1,1/4,1/4$ で $S=5/2$ である。$n=4$ なので $P$ の固有値は $7/24,7/24,-1/12,-1/12$ となる。よって第 1 因子では $W_{1212}=1-7/12=5/12$、第 2 因子でも $W_{3434}=1/4+1/6=5/12$ である。
:::
::::

::::exercise{level=advanced}
$n=3$ の Einstein 計量は定曲率であることを証明せよ。

:::solution
$\mathrm{Ric}=\lambda g$、$S=3\lambda$ で $P=\lambda g/4$ である。次元 $3$ では $W=0$ だから $\mathcal R=P\mathbin{\diamond} g$ となる。直交単位平面で $K=2(\lambda/4)=\lambda/2$ となり、縮約 Bianchi 恒等式により $\lambda$ は連結成分ごとに一定である。
:::
::::

## 確認テスト

::::quiz
平坦性の判定として正しいものを選べ。

- [x] 近傍で $R=0$ なら局所 Euclid 座標が取れる。
- [ ] 一点で $\Gamma=0$ なら近傍全体で $R=0$ である。

:::explanation
前者は :ref[thm-rg-flat] で証明した。後者では $\Gamma$ の微分が残る。例えば球面でも一点の正規座標では接続係数が零だが曲率は正である。
:::
::::

::::quiz{answer="6" tol="1e-6"}
三次元定曲率 $1$ の空間でスカラー曲率はいくつか。

:::explanation
三種類の独立な座標断面を二回ずつ加えるので $S=2\cdot3=6$ である。
:::
::::

::::quiz{answer="-1" tol="1e-6"}
$dr^2+e^{2r}d\theta^2$ の断面曲率を求めよ。

:::explanation
$f=e^r$ は $f''=f$ を満たすから曲率は $-1$ である。
:::
::::

::::quiz
正しい主張をすべて選べ。

- [x] 四次元以上にも Einstein で断面曲率が一定でない例がある。
- [x] Weyl テンソルの Ricci 縮約は零である。
- [ ] Weyl テンソルが零ならスカラー曲率も零である。

:::explanation
前者の例は $S^2(1)\times S^2(1)$ である。縮約は Schouten テンソルによって取り除いた。球面は $W=0$ かつ $S>0$ なので最後は誤りである。
:::
::::

---
title: "偏微分方程式への応用"
description: "熱、波動、Laplace、Schrödinger の各方程式と Radon 変換を Fourier 法で解く。"
notation:
  - symbol: 'G_{\kappa t}'
    meaning: "全空間の熱核"
  - symbol: 'M_r f'
    meaning: "半径 r の球面上の平均"
  - symbol: 'P_y'
    meaning: "上半平面の Poisson 核"
  - symbol: 'U(t)'
    meaning: "自由 Schrödinger 発展"
  - symbol: '\mathcal R f'
    meaning: "Radon 変換"
---

## この章の目標

- 熱方程式の周波数解を、円周と全空間の核に戻す。
- 波動方程式の球面平均表示を導き、伝播範囲を求める。
- 上半平面の境界値問題と自由 Schrödinger 方程式を解く。
- Fourier 切片定理から CT の逆投影公式を導く。

## 円周上と ℝ 上の熱方程式

一定の熱拡散率 $\kappa>0$ に対し $u_t=\kappa\Delta u$ を考える。一次元では $\Delta u=u_{xx}$ である。位置の単位が m、時間が s なら $\kappa$ は $\mathrm{m^2/s}$ である。円周では半径を長さの単位に選び、角度 $\theta$ に対する係数を $\kappa$ と書く。半径 $R$ で角度を無次元とすると、その係数は物理的な拡散率を $R^2$ で割ったものである。

:::theorem[熱半群と初期値]{#thm-four-heat}
$f\in L^2(\mathbb T)$ に対する熱方程式の解は $t>0$ で

$$
u(t,\theta)=\sum_{k\in\Z}e^{-\kappa k^2t}\widehat f(k)e^{ik\theta}.
$$

全空間 $\R^n$ で $f\in L^2$ なら $\widehat u(t,\xi)=e^{-4\pi^2\kappa t|\xi|^2}\widehat f(\xi)$ であり、

$$
u(t)=G_{\kappa t}*f,\qquad G_s(x)=(4\pi s)^{-n/2}e^{-|x|^2/(4s)}.
\label{eq-four-heat-kernel}
$$

両者は $t\downarrow0$ で $L^2$ 収束し、$C([0,\infty);L^2)\cap C^1((0,\infty);L^2)$ で $\Delta u(t)\in L^2$ を満たす解として一意である。
:::

:::proof
円周の各係数は $c_k'=-\kappa k^2c_k$ を満たす。係数の解を代入すると、任意の正の時刻以後では多項式因子を掛けても級数が一様絶対収束する。このため項別微分できる。実空間では変換で $\Delta$ が $-4\pi^2|\xi|^2$ になる。Gauss 変換と畳み込み公式が核を与える。$\widehat f\in L^2$ に指数因子と任意の多項式を掛けたものは $t>0$ で $L^1$ と $L^2$ に属するため、逆積分を微分できる。

初期値の収束と時刻の $L^2$ 連続性は :ref[thm-four-circle-parseval] または :ref[thm-four-plancherel] に優収束を組み合わせれば従う。一意性は差の周波数を $|k|\le M$ または $|\xi|\le M$ に切断する。切断上で係数の乗算は有界作用素なので、$e^{\kappa tk^2}$ または $e^{4\pi^2\kappa t|\xi|^2}$ を掛けた差の $L^2$ 微分はゼロである。初期値がゼロなので差はゼロであり、$M\to\infty$ で一意性が従う。
:::

時刻 $t$ に初期値を解へ送る作用素を $S(t)$ と書く。乗数の積から $S(t)S(s)=S(t+s)$ であり、$S(0)=I$ と置ける。この作用素族を :term[熱半群]{en="heat semigroup" yomi="ねつはんぐん"} という。

円周の核 $K_t=\sum_k e^{-\kappa tk^2}e^{ik\theta}$ は実空間の熱核の周期化 $2\pi\sum_{m\in\Z}G_{\kappa t}(\theta+2\pi m)$ に等しい。係数を計算すると実空間の変換の周波数 $k/(2\pi)$ での値になり、一意性定理が等式を与える。係数 $2\pi$ は畳み込みに使う測度が $d\theta/(2\pi)$ だからである。両核は正で質量 $1$ であり、初期値の上下界を保存する。

:::example[熱による高周波の消滅]
$\kappa=1$、$f=2+\cos\theta+3\sin(2\theta)$ なら $u=2+e^{-t}\cos\theta+3e^{-4t}\sin(2\theta)$ である。$t=\log 2$ で二つの振幅は $1/2$ と $3/16$ になる。一般の $L^2$ データでも平均を除いた二乗ノルムは $e^{-2\kappa t}$ 倍以下に減る。係数 $k\ne0$ の減衰率に $k^2\ge1$ を使えばよい。
:::

:::example[点から広がる熱]
一次元で初期値を質量 $1$ の $\delta_0$ と読むと、正の時刻では $G_{\kappa t}$ になる。$\int x^2G_{\kappa t}(x)dx=2\kappa t$ であり、$\kappa=10^{-6}\,\mathrm{m^2/s}$、$t=50\,\mathrm{s}$ なら標準偏差は $0.01\,\mathrm m$ である。分散は Gauss 積分を指数の係数で微分して求められる。
:::

### 熱源がある場合と Duhamel の公式

熱源はその時刻に加えた新たな初期値として扱える。その後に残る寄与は、加えた時刻から現在までの熱半群を適用して求められる。それを時刻について積分する表示を :term[Duhamel の公式]{en="Duhamel formula" yomi="でゅあめるのこうしき"} という。

:::proposition[非斉次熱方程式の時間積分]{#prop-four-heat-duhamel}
$A=\kappa\Delta$ とし、$e^{tA}$ は上で得た熱半群とする。$f\in\mathcal S$、$F:[0,\infty)\to\mathcal S$ が Schwartz 空間の位相で連続なら、$u_t=Au+F(t)$、$u(0)=f$ の解は
$$
u(t)=e^{tA}f+\int_0^t e^{(t-s)A}F(s)\,ds
$$
である。円周では滑らかな周期関数と導関数の上限を使えば同じ式が成り立つ。
:::

:::proof
周波数を固定すると一階方程式 $v'=-\omega v+\widehat F$ となり、積分因子 $e^{\omega t}$ で解くと式を得る。空間微分を乗数に移しても各時間区間で急減少する周波数積分で支配されるため、逆変換と時間積分が交換できる。時間微分では積分の上端から $F(t)$ が出て、残りは $A$ を掛けた積分となる。差は斉次の熱方程式を満たすので :ref[thm-four-heat] の一意性を使える。
:::

:::example[一定の周期熱源]
円周で $f=0$、$F(t,\theta)=\cos(2\theta)$、$\kappa=1$ とする。公式から解の係数は $\int_0^t e^{-4(t-s)}ds=(1-e^{-4t})/4$ である。長時間後には $\cos(2\theta)/4$ に近づく。定常解を $\Delta u+F=0$ から求めても同じである。一方、空間定数の熱源では減衰率がゼロなので、時間積分は $t$ に比例する。
:::

### 正則化の尺度と最大値の原理

熱核を微分すると空間の尺度 $(\kappa t)^{1/2}$ が分母に現れる。したがって初期値が積分可能なだけでも、正の時刻で得られる導関数を定量的に制御できる。この評価では、入力の大きさを測る指数と出力を測る指数を別々に選べる。

:::proposition[熱半群の $L^p$–$L^q$ 評価]{#prop-four-heat-smoothing}
$1\le p\le q\le\infty$、$\alpha$ を多重指数、$s=\kappa t>0$ とする。$f\in L^p(\R^n)$ に対し
$$
\|\partial_x^\alpha(G_s*f)\|_q
 \le C_{n,\alpha,p,q}\,s^{-|\alpha|/2-\frac n2(1/p-1/q)}\|f\|_p .
$$
$f$ が実数値で $A\le f\le B$ a.e. なら $A\le (G_s*f)(x)\le B$ である。
:::

:::proof
$G_s(x)=s^{-n/2}G_1(x/\sqrt s)$ なので $\partial^\alpha G_s=s^{-(n+|\alpha|)/2}(\partial^\alpha G_1)(x/\sqrt s)$ である。$\partial^\alpha G_1$ は多項式と Gauss 関数の積であり、全ての $L^r$（$1\le r\le\infty$）に入る。したがって
$$
\|\partial^\alpha G_s\|_r
=s^{-|\alpha|/2-\frac n2(1-1/r)}\|\partial^\alpha G_1\|_r .
$$
$r$ を $1/p+1/r=1+1/q$ で選び、Young の不等式を適用する。微分と畳み込みの交換は、正の $s$ の近傍で核とその微分が $C(1+|x|)^{-n-1}$ で支配されることから Hölder と優収束で従う。最後は $G_s\ge0$、$\int G_s=1$ だから $A G_s\le G_s f(x-\cdot)\le B G_s$ を積分すればよい。
:::

:::example[区間を暖めた直後の温度]
$n=1$、$\kappa=1$、$f=1_{[-1,1]}$ とする。原点では
$$
u(t,0)=\frac1{\sqrt\pi}\int_{-1/(2\sqrt t)}^{1/(2\sqrt t)}e^{-v^2}dv
$$
である。$t=1/4$ なら積分範囲は $[-1,1]$ であり、$2e^{-1}/\sqrt\pi\le u(1/4,0)<1$ である。$x>1$ でも $u(t,x)>0$ となる。積分する区間の長さが正で、核は全空間で正だからである。熱方程式には明瞭な有限速度の境界がなく、後に扱う波動方程式と伝播の仕方が異なる。
:::

:::example[熱量がゼロの場合の速い減衰]
$n=1$、$f=1_{[0,1]}-1_{[-1,0]}$ とする。平均はゼロで、$\int|y||f(y)|dy=1$ である。よって
$$
u(t,x)=\int (G_t(x-y)-G_t(x))f(y)dy,\qquad
\|u(t)\|_\infty\le\|G_t'\|_\infty.
$$
$G_t'(x)=-xG_t(x)/(2t)$ である。$|x|e^{-x^2/(4t)}$ は $|x|=\sqrt{2t}$ で最大なので $\|G_t'\|_\infty=e^{-1/2}/(2\sqrt{2\pi}\,t)$ となる。一般の質量を持つ $L^1$ データより $t^{-1/2}$ 一因子だけ速い。
:::

### 大時間で残る質量と一次モーメント

:::proposition[最初の二つの熱漸近項]{#prop-four-heat-moment-asymptotic}
$f\in L^1(\R^n)$ かつ $\int |x|^2|f(x)|dx<\infty$ とし、$m=\int f$、$m_j=\int x_jf(x)dx$ とする。$1\le q\le\infty$ なら
$$
\left\|G_s*f-mG_s+\sum_{j=1}^n m_j\partial_jG_s\right\|_q
\le C_{n,q}s^{-1-\frac n2(1-1/q)}\int |x|^2|f(x)|dx .
$$
:::

:::proof
$y$ を積分変数として、一変数の Taylor 公式を $\tau\mapsto G_s(x-\tau y)$ に適用すると
$$
G_s(x-y)=G_s(x)-\sum_j y_j\partial_jG_s(x)
+\sum_{j,k}y_jy_k\int_0^1(1-\tau)\partial_j\partial_kG_s(x-\tau y)d\tau .
$$
$m,m_j$ を積分すると定数項と一次項が取り出される。余項の $L^q$ ノルムには Minkowski と平行移動の等長性を適用する。$q=\infty$ でも上限を積分の中に入れられる。$\sum_{j,k}|y_jy_k|\le n|y|^2$ と :ref[prop-four-heat-smoothing] の証明に現れた $r=q$ の核ノルムを使えば主張が従う。
:::

:::example[平行移動した Gauss 核の展開]
$f=G_1(\cdot-2)$ なら $m=1$、$m_1=2$、$\int x^2f(x)dx=6$ である。最後の値は中心座標 $z=x-2$ を用い、$\int z^2G_1(z)dz=2$ と $\int zG_1(z)dz=0$ を合わせた。正確な解は $u=G_{t+1}(\cdot-2)$ であり、上の式で $G_t-2G_t'$ との差の上限ノルムは $O(t^{-3/2})$ である。$2$ の移動を初期値の質量に含めるだけでは一次項が残り、精度は $O(t^{-1})$ となる。
:::

### 有限区間における境界条件の選択

円周で解いた波形を奇関数として延長すれば端の値をゼロにできる。偶関数として延長すれば端での垂直方向の微分をゼロにできる。周期を $2L$ とすれば、反射は二つの端 $0,L$ で同時に使える。

:::proposition[反射による区間の熱核]{#prop-four-heat-interval-images}
長さ $L>0$ の区間で $f\in L^2(0,L)$ とする。時刻 $0$ まで $L^2$ で連続であり、正の時刻で閉区間まで古典的な微分ができる熱解を考える。Dirichlet 条件 $u(t,0)=u(t,L)=0$ なら $K_D$、Neumann 条件 $u_x(t,0)=u_x(t,L)=0$ なら $K_N$ を用いて、唯一の解は $\int_0^L K(t,x,y)f(y)dy$ で表される。ただし $s=\kappa t$、
$$
K_{D,N}(t,x,y)=\sum_{m\in\Z}\bigl[G_s(x-y+2mL)\mp G_s(x+y+2mL)\bigr]
$$
であり、上の符号が $D$、下が $N$ である。
:::

:::proof
$f$ を $(-L,L)$ に奇または偶延長し、さらに $2L$ 周期とする。全空間核との畳み込みの積分を長さ $2L$ の各区間に分け、負の半区間の変数を反射する。すると二つの項の和が得られる。$t$ を正のコンパクトな区間に制限すれば、各導関数を含む核の和は $C(1+|m|)^{B}e^{-c m^2}$ で支配される。従って積分と項別微分を交換でき、方程式を満たす。

熱核の偶性により延長時の奇偶性が保たれる。周期と奇性から $u(t,0)=u(t,L)=0$、周期と偶性から $u_x(t,0)=u_x(t,L)=0$ を得る。初期値への $L^2$ 収束は、$x$ から $\theta=\pi x/L$ に置き換えた :ref[thm-four-heat] から従う。差の $L^2(0,L)$ の二乗を微分して部分積分すると境界項は $u_x\overline u$ であり、いずれの境界条件でもゼロである。従ってその二乗ノルムは非増加であり、初期値がゼロなら常にゼロである。
:::

:::example[断熱端と温度固定端]
$L=\pi$、$\kappa=1$ とする。初期値 $f=1+\cos x$ は偶延長で滑らかであり、Neumann 解は $1+e^{-t}\cos x$ である。端の温度は $1\pm e^{-t}$ だが、その間の熱量 $\int_0^\pi u\,dx=\pi$ は保存される。Dirichlet 条件の例として $f=\sin x$ を取ると、奇延長は $\sin x$ なので $u=e^{-t}\sin x$ となり、熱量は $2e^{-t}$ である。境界の熱流がゼロであるかを確認すると、熱量の時間微分 $\kappa[u_x(t,L)-u_x(t,0)]$ の違いを直接読める。
:::

### 逆向きの熱方程式と高周波

順向きの一意性は、固定時刻で観測したデータから初期値を安定に求められることを意味しない。例えば円周上で $u_j(0,\theta)=\sqrt2\cos(j\theta)$ は規格化測度でノルム $1$ である。正の時刻 $T$ ではノルムが $e^{-\kappa j^2T}$ になる。したがって $u(T)\mapsto u(0)$ の作用素は $L^2$ 上の値域に限っても有界でない。対照的に、周波数を $|j|\le N$ に限れば最大増幅率は $e^{\kappa N^2T}$ であり、$j=\pm N$ で実現する。初期値復元で帯域を制限するのは、ノイズの増幅と捨てる成分とを両方制御するためである。

:::example[増加中の熱源と減衰するモード]
初期値をゼロ、$F(t,\theta)=e^{-4t}\cos2\theta$、$\kappa=1$ とすると、積分因子の中で二つの指数が消えて $u(t,\theta)=te^{-4t}\cos2\theta$ になる。係数を微分すると $e^{-4t}(1-4t)$ で、最大値は $t=1/4$ で $1/(4e)$ である。熱源が単調に弱まっても、短い時刻には供給した寄与が蓄積するため、解の振幅が増える。
:::

## 波動方程式と球面平均

波動方程式 $u_{tt}=c^2\Delta u$ では、各周波数が振動数を決める。$c>0$ は波の速度で、例えば m/s で測る。$f=u(0)$、$g=u_t(0)$ を Schwartz 関数とする。

:::proposition[波動方程式の周波数解]{#prop-four-wave}
$u\in C^2(\R_t;L^2(\R^n_x))$、$\Delta u(t)\in L^2$ の解は一意であり、

$$
\widehat u(t,\xi)=\cos(2\pi ct|\xi|)\widehat f(\xi)
+\frac{\sin(2\pi ct|\xi|)}{2\pi c|\xi|}\widehat g(\xi)
$$

であり、第 2 乗数は原点で $t$ とする。一次元では

$$
u(t,x)=\tfrac12[f(x+ct)+f(x-ct)]+\frac1{2c}\int_{x-ct}^{x+ct}g(y)\,dy
$$

となる。
:::

:::proof
変換後の方程式は二階の定係数常微分方程式である。余弦と正弦の係数を初期条件で決める。時刻微分後も多項式増大以下なので Schwartz データで逆変換と微分を交換できる。一意性は周波数を有界集合へ切断して成立する常微分方程式の一意性から従う。一次元では余弦を指数二つに分けて移動に戻す。区間指示関数 $1_{[-ct,ct]}/(2c)$ の変換が正弦の乗数に一致するので積分項が得られる。
:::

:::theorem[三次元の Kirchhoff 公式と二次元への降下]{#thm-four-wave-means}
三次元で $M_rh(x)=(4\pi)^{-1}\int_{S^2}h(x+r\omega)\,d\sigma(\omega)$ とする。$t>0$ の解は

$$
u(t,x)=\partial_t\bigl(tM_{ct}f(x)\bigr)+tM_{ct}g(x).
$$

二次元では $M_r$ を $(2\pi r)^{-1}\int_{|y|<r}h(x+y)(r^2-|y|^2)^{-1/2}dy$ に置き換える。
:::

:::proof
三次元では周波数の向きを極軸に選ぶと球面平均の乗数は $\frac12\int_{-1}^1 e^{2\pi ir|\xi|s}ds=\sin(2\pi r|\xi|)/(2\pi r|\xi|)$ である。$r=ct$ とし $t$ を掛ければ正弦の項となり、それを $t$ で微分すれば余弦の項となる。従って :ref[prop-four-wave] に一致する。

二次元のデータを第 3 座標に依存しないよう三次元へ延長する。公式はコンパクト球面上の積分なので延長に対して定義できる。各有界な $(t,x)$ の範囲では、そこで使うデータだけを残して外側を滑らかに切断しても公式とその微分は一致する。したがって同じ公式が三次元の方程式を満たし、第 3 座標に依存しないため二次元の方程式を満たす。半径 $r$ の球面を平面の円板に射影すると上下二枚のグラフがあり、各面積要素は $r/\sqrt{r^2-|y|^2}$ 倍である。球面積 $4\pi r^2$ で割り二枚分を足すと二次元の核を得る。
:::

:::example[球面上と球面内の違い]
$g$ が原点中心の半径 $1$ の球内に台を持ち、$f=0,c=2$ とする。三次元の原点で $t>1/2$ なら積分球面は台の外にあり、$u(t,0)=0$ である。二次元では $g\ge0$ が非零なら $t>1/2$ でも円板積分の核は正なので $u(t,0)>0$ である。どちらの公式でも、データが観測点から距離 $ct$ より遠い場合にはその影響は現れない。
:::

波のエネルギー $E(t)=\frac12\int(|u_t|^2+c^2|\nabla u|^2)dx$ は、Schwartz データでは時刻によらない。時間微分で $u_{tt}$ を方程式に置き換え、$\Delta u$ を一度部分積分すると二つの項が相殺する。熱方程式の二乗ノルムの時間微分は逆に $-2\kappa\int|\nabla u|^2dx$ である。波動方程式では運動とひずみのエネルギーが互いに移る。熱方程式では Laplace 作用素が各非零周波数を直接減衰させる。

:::example[二次元の円板平均を定数で確かめる]
二次元の核の質量は
$$
\frac1{2\pi r}\int_{|y|<r}\frac{dy}{\sqrt{r^2-|y|^2}}
=\frac1r\int_0^r\frac{s\,ds}{\sqrt{r^2-s^2}}=1
$$
である。このため $g=1$ が観測点の半径 $ct$ の円板の近傍で成り立ち、$f=0$ なら $u(t,x)=t$ となる。全空間で定数の初期速度を与えなくても、その時刻の値は同じである。
:::

### 円錐上のエネルギーと依存領域

球面平均の明示公式から一次元から三次元までの有限伝播速度が分かる。任意の次元でも、移動する球の中でエネルギーが増えないことを直接示せる。この方法は核を計算できない場合にも使われる。

:::proposition[後方円錐での局所エネルギー]{#prop-four-wave-local-energy}
$u$ は $u_{tt}=c^2\Delta u$ の滑らかな解とし、固定した $x_0\in\R^n$、$T>0$ に対して
$$
E(s)=\frac12\int_{|x-x_0|<c(T-s)}\big(|u_t(s,x)|^2+c^2|\nabla u(s,x)|^2\big)dx
\qquad(0\le s<T)
$$
と置く。$E$ は非増加である。従って初期位置と初期速度が球 $B(x_0,cT)$ 内でゼロなら $u(T,x_0)=0$ である。
:::

:::proof
エネルギー密度を $e=(|u_t|^2+c^2|\nabla u|^2)/2$ とすると、方程式と積の微分から $\partial_t e=c^2\operatorname{div}(\Re(\overline{u_t}\nabla u))$ である。動く球を極座標で微分すると、外側の半径の速度 $-c$ が面積積分に掛かる。発散定理と合わせて
$$
E'(s)=\int_{|x-x_0|=c(T-s)}
\left[c^2\Re(\overline{u_t}\partial_\nu u)-ce\right]d\sigma .
$$
$2|\overline{u_t}\,c\partial_\nu u|\le|u_t|^2+c^2|\partial_\nu u|^2$ と $|\partial_\nu u|\le|\nabla u|$ により被積分関数は非正である。$E(0)=0$ なら円錐の内部全体で $u_t=0$ となる。特に中心線では $u(s,x_0)=u(0,x_0)=0$ が $s<T$ で成り立ち、連続性で頂点でもゼロとなる。
:::

:::example[二つの初期波形が同じ値を与える範囲]
速度 $c=3$ とし、二つの初期位置と速度が $|x|<6$ で一致するとする。差を上の命題に入れると原点で $0\le t\le2$ の解は一致する。同様に $x_0$ を原点から距離 $3$ の点に移すと $B(x_0,3t)$ が元の球に入るのは $t\le1$ である。従って一致区間の保証は短くなる。境界の外から観測点へ最短で到達する距離が評価を決める。
:::

### 球面平均の微分方程式

:::lemma[球面平均の Darboux 方程式]{#lem-four-wave-darboux}
$n\ge2$、$f\in C^\infty(\R^n)$ とする。正規化球面平均 $M_rf(x)=|S^{n-1}|^{-1}\int_{S^{n-1}}f(x+r\omega)d\sigma(\omega)$ は
$$
\partial_r^2 M_rf+\frac{n-1}{r}\partial_r M_rf=\Delta_x M_rf,\qquad
M_0f=f,\quad \partial_r M_0f=0
$$
を満たす。
:::

:::proof
$r$ について微分すると $\partial_r M_r f$ は法線微分の球面平均であり、半径 $r$ の面積測度に換算すると
$$
r^{n-1}\partial_r M_r f
=|S^{n-1}|^{-1}\int_{B(x,r)}\Delta f(z)dz .
$$
ここで発散定理を用いた。もう一度 $r$ で微分すると右辺は $r^{n-1} M_r(\Delta f)$ になる。$x$ による微分は元の積分の内側に入れられるので右辺は $r^{n-1}\Delta_x M_rf$ である。初期値は連続性から得られ、初期微分のゼロは各 $\omega_j$ の平均がゼロであることから従う。
:::

:::example[四次元で球面平均の二次項を調べる]
$f(x)=|x|^2$ なら球面上での一次の交差項 $2r x\cdot\omega$ は平均ゼロなので $M_r f=|x|^2+r^2$ である。$n=4$ で左辺は $2+3(2r)/r=8$、右辺も $\Delta_x|x|^2=8$ となる。球面平均そのものは波動方程式の解ではなく、余分な $r^{-1}\partial_r$ がある。この項の係数は球の面積が $r^{n-1}$ に比例することから来る。
:::

### 奇数次元の球面公式

:::proposition[奇数次元における球面表示]{#prop-four-wave-odd-dimension}
$n=2m+1\ge3$、$m\ge1$、$r=ct$ とする。滑らかな初期速度 $g$ と初期位置 $0$ に対し
$$
u(t,x)=\frac1{c(2m-1)!!}
      \left(\frac1r\partial_r\right)^{m-1}\big(r^{2m-1}M_rg(x)\big)\Big|_{r=ct}.
$$
初期位置 $f$、初期速度 $0$ の解は、右辺で $g$ を $f$ に置き換えた式をさらに $t$ で微分して得られる。
:::

:::proof
$F=M_r g$ とし $W=r^{n-2}F$ と置く。:ref[lem-four-wave-darboux] の左辺に代入して積を微分すると $W''+(3-n)W'/r=\Delta_x W$ である。$D=r^{-1}\partial_r$、$L_b=\partial_r^2+(b/r)\partial_r$ とすると
$$
L_{b+2}D v=D L_b v
$$
である。実際、両辺は $v'''/r+bv''/r^2-bv'/r^3$ に等しい。このため $D^{m-1}W$ には $b=3-n+2(m-1)=0$ の方程式が成り立つ。$r=ct$ へ戻すと $u_{tt}=c^2\Delta_xu$ である。

原点では $F(r,x)$ は $r$ の滑らかな偶関数で $F(0,x)=g(x)$ である。Taylor 公式を十分高い階数まで使うと、$D^{m-1}(r^{2m-1}F)=(2m-1)!!\,rg+O(r^3)$ であり、微分した余項も $O(r^2)$ である。実際 $D r^a=ar^{a-2}$ であり、$r^{2m-1}$ から $r$ に至る積は $(2m-1)(2m-3)\cdots3$ である。これにより $u(0)=0$、$u_t(0)=g$ を得る。$t$ 微分した解の初期値は $g$、初期速度は $0$ である。局所一意性は :ref[prop-four-wave-local-energy] を差に適用して従う。
:::

:::example[五次元で一定の初期速度]
$n=5,c=2$ では $D(r^3M_rg)=3rM_rg+r^2\partial_rM_rg$ である。$g$ が観測点を中心とする半径 $2t$ の球の近くで定数 $6$ なら、その時刻の値は $3r\cdot6/(2\cdot3)=6t$ である。球面近傍で $g=0$ なら平均とその $r$ 微分はともにゼロであり、球の内部に初期速度が残っていても、その時刻での寄与はゼロとなる。
:::

### 有限エネルギーと外力の Duhamel 積分

:::proposition[エネルギーデータと外力]{#prop-four-wave-energy-duhamel}
$A$ を周波数で $2\pi c|\xi|$ を掛ける作用素とし、$C(t)=\cos(tA)$、$S(t)=\sin(tA)/A$ とする。$f\in H^1$、$g\in L^2$、$F\in C([0,T];L^2)$ なら、$u_{tt}=c^2\Delta u+F$ の分布解で $u\in C([0,T];H^1)\cap C^1([0,T];L^2)$ に属するものは一意で、
$$
u(t)=C(t)f+S(t)g+\int_0^tS(t-s)F(s)ds.
$$
また $e(t)=(\|u_t(t)\|_2^2+c^2\|\nabla u(t)\|_2^2)^{1/2}$ に対し $e(t)\le e(0)+\int_0^t\|F(s)\|_2ds$ である。
:::

:::proof
$|S(t)|$ の乗数は $|t|$ 以下、$|A S(t)|$ と $|C(t)|$ は $1$ 以下である。したがって $S(t):L^2\to H^1$ は有界時間区間で一様有界で、その強連続性は優収束から従う。時間微分 $S'(t)=C(t)$ は $L^2$ で行え、積分の上端では $S(0)=0$ なので現れない。二階目は周波数を有限半径以内に切って常微分方程式として確認する。切断を除く極限は $u\in C(H^1)$、$F\in C(L^2)$ で分布の意味で取れる。同じ有界周波数切断を解の差に行うと初期値ゼロの有界作用素常微分方程式となり、各切断で差はゼロである。

周波数ごとの二成分 $(A\widehat u,\widehat u_t)$ は、斉次部分では行列
$$
\begin{pmatrix}\cos(tA)&\sin(tA)\\-\sin(tA)&\cos(tA)\end{pmatrix}
$$
を初期成分 $(A\widehat f,\widehat g)$ に掛けたものである。実直交行列で複素成分の二乗和も保存されるから、Plancherel によりその $L^2$ ノルムは $e(0)$ である。外力の各時刻 $s$ の寄与は、同じ行列を $(0,\widehat F(s))$ に時間 $t-s$ だけ作用させたものである。積分の三角不等式で評価すればよい。
:::

:::example[外力を加え始めた直後の周期波]
円周上で波速 $c=1$、$f=g=0$、外力 $F(t,\theta)=\cos(2\theta)$ とする。周期版の乗数で $S(t)$ はこのモードに $\sin(2t)/2$ を掛ける。時間積分で $(1-\cos(2t))/4$ が得られるので $u(t,\theta)=(1-\cos(2t))\cos(2\theta)/4$ である。$t=\pi/2$ では振幅は $1/2$ だが時間微分はゼロである。力が一定でも時間の効果は単調ではない。
:::

## 上半平面の Laplace 方程式と Poisson 核

半平面の高さ $y$ を増やすと境界データの細かい振動は減る。$x$ について変換すると $y$ について二つの指数解が得られ、上半平面での有界性が成長する解を除く。

:::theorem[上半平面の Dirichlet 問題]{#thm-four-upper-poisson}
有界連続な境界データ $f$ に対し、有界で閉半平面まで連続な調和拡張は一意であり、

$$
u(x,y)=(P_y*f)(x),\qquad P_y(x)=\frac{y}{\pi(x^2+y^2)}
$$

で与えられる。
:::

:::proof
Schwartz データでは水平変換により $u_{yy}$ が $4\pi^2|\xi|^2\widehat u$ に対応する。減衰する解の乗数 $e^{-2\pi y|\xi|}$ を逆変換すると前章の Poisson 核になる。一般の有界連続データに対しては $P_y$ の二階までの微分を局所的な $y>0$ の範囲で $C/(1+x^2)$ で支配できるため、畳み込みを微分できる。$P_y$ は $i/(\pi(x+iy))$ の実部なので調和的である。質量は $\arctan(x/y)/\pi$ を両端で評価して $1$ であり、$y\downarrow0$ で原点外の質量は消える。境界点で座標が動く場合も、局所連続性と遠方の質量の評価に分ければ $u(x,y)\to f(x_0)$（$(x,y)\to(x_0,0)$）となる。

一意性には差 $w$ が有界で実軸上でゼロであることを用いる。複素数 $z=x+iy$ とし、半円板 $|z|<R,y>0$ 上で $v_R(z)=(2/\pi)\arg((R+z)/(R-z))$ を取る。この商は半円板内で第 1 象限に入り、直径上で正の実数、円弧上で正の純虚数になる。従って $v_R$ は調和、非負で、直径で $0$、円弧で $1$ である。端点で $\limsup(|w|-\|w\|_\infty v_R)\le0$ は $w$ の連続性から従う。最大値原理を実部と虚部の正負に適用し、固定 $z$ で $R\to\infty$ とすれば $v_R(z)\to0$ なので $w(z)=0$ である。
:::

:::example[区間の調和測度]
$f=1_{[-1,1]}$ なら

$$
u(x,y)=\pi^{-1}\left(\arctan\frac{1-x}{y}+\arctan\frac{1+x}{y}\right).
$$

原点の真上 $y=1$ では $u=1/2$ である。また $x=1$ に沿って $y\downarrow0$ とすると $u\to1/2$ であり、跳躍の左右の平均が得られる。
:::

たとえば境界波 $f(x)=\cos(2\pi x)$ には $u(x,y)=e^{-2\pi y}\cos(2\pi x)$ が対応する。正負の周波数を二つの指数関数に分け、$P_y$ の変換を $\xi=\pm1$ で評価すればよい。

:::example[点の真上での二つの高さ]
境界値 $f=P_2$ とすると $\widehat f=e^{-4\pi|\xi|}$ である。高さ $3$ の Poisson 核を畳み込むと乗数が $e^{-10\pi|\xi|}$ となり、解は $P_5$ である。従って $x=0$ の値は $1/(5\pi)$、$x=5$ ではその半分の $1/(10\pi)$ である。
:::

### Dirichlet 値から Neumann 値を求める

Poisson 積分は境界の値を半平面へ延ばす。拡張された関数の法線微分を境界で取ると、境界上だけで定義された別の作用素を得る。上半平面の外向き法線は下向きであるため、法線微分は $-\partial_y$ になる。

:::proposition[上半平面の Dirichlet–Neumann 作用素]{#prop-four-half-plane-dirichlet-neumann}
$f\in\mathcal S(\R)$、$u=P_y*f$ とする。境界の法線微分 $\Lambda f=\lim_{y\downarrow0}(-\partial_yu)$ は一様かつ $L^2$ で存在し、
$$
\widehat{\Lambda f}(\xi)=2\pi|\xi|\widehat f(\xi),\qquad
\Lambda f(x)=\frac1\pi\lim_{\eps\downarrow0}
       \int_{|t|>\eps}\frac{f(x)-f(x-t)}{t^2}dt .
$$
:::

:::proof
周波数では $-\partial_yu$ の乗数は $2\pi|\xi|e^{-2\pi y|\xi|}$ である。$|\xi|\widehat f$ は $L^1\cap L^2$ なので、優収束と反転、Plancherel により両方の極限を得る。空間の表示を得るには
$$
\frac{f(x)-P_y*f(x)}y
 =\frac1\pi\int_\R\frac{f(x)-f(x-t)}{t^2+y^2}dt
$$
と書く。$|t|\le1$ では被積分関数を $t$ と $-t$ で平均し、分子を $f(x)-(f(x-t)+f(x+t))/2$ に置く。この分子は $\|f''\|_\infty t^2/2$ 以下である。遠方では $2\|f\|_\infty/t^2$ によって支配できる。従って優収束で主値を得る。一方、差商の乗数 $(1-e^{-2\pi y|\xi|})/y$ は $2\pi|\xi|$ 以下で同じ極限を持つ。両表示は同じ作用素を与える。
:::

:::example[境界に正弦波を置く]
$f(x)=\cos(4\pi x)$ では $u=e^{-4\pi y}\cos(4\pi x)$ なので、境界の法線微分は $4\pi\cos(4\pi x)$ である。$f$ は Schwartz ではないが、核積分で同じ公式を確認できる。$x=0$ で主値は
$$
\frac2\pi\int_0^\infty\frac{1-\cos(4\pi t)}{t^2}dt
=8\int_0^\infty\frac{\sin(4\pi t)}{t}dt=4\pi
$$
である。第 1 等号では部分積分し、端点の境界項がゼロであることを使った。第 2 等号は :ref[prop-four-sine-integral-limit] の変数変換による。正負の周波数は同じ符号の係数を得る。
:::

### 調和拡張が最小にするエネルギー

:::proposition[Dirichlet エネルギーの最小化]{#prop-four-half-plane-dirichlet-principle}
$f\in\mathcal S(\R)$、$u=P_y*f$ とし、$\mathcal E(v)=\int_{y>0}(|v_x|^2+|v_y|^2)dxdy$ とする。$\mathcal E(u)=2\pi\int|\xi||\widehat f(\xi)|^2d\xi$ である。$w$ が閉半平面上で滑らか、コンパクト台を持ち、境界でゼロなら $\mathcal E(u+w)=\mathcal E(u)+\mathcal E(w)$ である。従ってこの比較類で $u$ は一意の最小化関数である。
:::

:::proof
水平と高さの微分はそれぞれ $2\pi i\xi e^{-2\pi y|\xi|}$ と $-2\pi|\xi|e^{-2\pi y|\xi|}$ を $\widehat f$ に掛ける。Plancherel の後、二乗を足して $y$ で積分すると $\xi\ne0$ で $8\pi^2|\xi|^2/(4\pi|\xi|)=2\pi|\xi|$ となる。Tonelli は非負の積分に使える。交差項 $\int\nabla u\cdot\overline{\nabla w}$ は部分積分で $-\int\Delta u\,\bar w$ に等しい。遠方の境界では $w=0$、下の境界でも $w=0$ だからである。$u$ は調和的なので交差項はゼロであり、エネルギーの恒等式が従う。$\mathcal E(w)=0$ なら $w$ の二つの一次微分が連続関数としてゼロであり、境界値から $w=0$ である。
:::

:::example[幅によって変わらない二次元のエネルギー]
$f_a(x)=e^{-\pi x^2/a^2}$（$a>0$）なら $\widehat f_a(\xi)=a e^{-\pi a^2\xi^2}$ である。エネルギーは
$$
4\pi a^2\int_0^\infty \xi e^{-2\pi a^2\xi^2}d\xi=1
$$
となり、$a$ によらない。幅を $a$ 倍にすると拡張の高さ方向の尺度も $a$ 倍になる。面積は $a^2$ 倍、導関数の二乗は $a^{-2}$ 倍となるので、この二次元の尺度変換でも同じ結論を読める。
:::

### 帯領域の二つの境界

半平面では減衰する指数を一つ選べば境界値が定まる。幅が有限の帯では上下の値を指定できるため、二つの指数を残して二本の条件を解く。そのとき原点の周波数だけは指数の根が重なるので、一次関数の補間になる。

:::proposition[帯領域の Dirichlet 問題]{#prop-four-laplace-strip}
$L>0$、$\Omega=\{(x,y):0<y<L\}$ とする。$f_0,f_L\in\mathcal S(\R)$ に対して閉帯まで連続で、境界値がそれぞれ $f_0,f_L$ の有界調和解は一意であり、水平 Fourier 変換で
$$
\widehat u(\xi,y)=
\frac{\sinh((L-y)2\pi|\xi|)}{\sinh(L2\pi|\xi|)}\widehat f_0(\xi)
+\frac{\sinh(y2\pi|\xi|)}{\sinh(L2\pi|\xi|)}\widehat f_L(\xi)
$$
となる。$\xi=0$ での比は順に $(L-y)/L$、$y/L$ とする。
:::

:::proof
各 $\xi$ に対して $v''=(2\pi|\xi|)^2v$ の二つの境界値を満たす式を解くと上の比になる。両比は $[0,L]$ で $0$ 以上 $1$ 以下である。$\xi=0$ 付近では $\sinh$ の奇冪級数で共通の一次因子を約分すると、両比は $\xi$ の滑らかな偶関数になる。従って固定した $y$ で、遠方での指数表示と合わせると右辺は Schwartz 関数である。$y$ の導関数は固定した $y$ の近傍で高々多項式増大である。実際、1 階で現れる $\cosh(a 2\pi|\xi|)/\sinh(L2\pi|\xi|)$ に $2\pi|\xi|$ が掛かり、$\xi$ が小さいときはこの積が有限である。偶数階は方程式を使い、奇数階には同じ 1 階の上界を使う。水平微分も多項式を掛けるだけなので、急減少する $\widehat f_j$ によって逆変換と全ての微分が交換できる。境界値は優収束と反転で得る。

差 $w$ は有界、閉帯まで連続で、上下でゼロである。$0<\lambda<\pi/L$ とし、$v(x,y)=\cosh(\lambda x)\cos(\lambda(y-L/2))$ を置く。$v$ は調和的で正であり、$|x|=R$ での最小値は無限に増える。任意の $\eps>0$ について十分大きい長方形で $\pm w\le\eps v$ を境界で満たし、最大値原理を適用できる。複素値では各位相を掛けた実部について同じことを行う。固定点で $\eps\downarrow0$ とすれば $w=0$ である。
:::

:::example[上の境界が下からの波を弱める割合]
$L=1$、下の波が $\cos(2\pi x)$、上の値はゼロとする。直接代入して確かめると $y=1/2$ で振幅は $\sinh\pi/\sinh(2\pi)=1/(2\cosh\pi)$ である。半平面での対応する減衰率は $e^{-\pi}$ だから、帯の場合の比は $e^{-\pi}/(1+e^{-2\pi})$ である。高さが固定された上の境界の影響が分母に現れる。一方、下と上に定数 $2,4$ を置くと $u=2+2y$ であり、真ん中は $3$ となる。ゼロ周波数の極限は定数の上下条件を反映する。
:::

## Schrödinger 方程式と分散評価

$i u_t=-a\Delta u$（$a>0$）を考える。自由粒子では $a=\hbar/(2m)$ で単位は $\mathrm{m^2/s}$ である。変換後に現れる指数は絶対値 $1$ なので、空間の二乗ノルムを保存する。一方で異なる速度の位相が離れていくため、上限ノルムは小さくなり得る。

:::theorem[自由発展と分散評価]{#thm-four-schrodinger}
$U(t)f=\mathcal F^{-1}(e^{-4\pi^2iat|\xi|^2}\widehat f)$ は $L^2$ のユニタリ群である。$f\in L^1\cap L^2$、$t\ne0$ では

$$
U(t)f(x)=(4\pi iat)^{-n/2}\int e^{i|x-y|^2/(4at)}f(y)dy,\qquad
\|U(t)f\|_\infty\le(4\pi a|t|)^{-n/2}\|f\|_1 .
$$

複素冪は右半平面から $iat$ へ近づけた値を取る。
:::

:::proof
乗数の積で群則が、乗数の絶対値と Plancherel でユニタリ性が従う。時刻については優収束で強連続であり、Schwartz データでは乗数を微分して方程式を得る。

$\eps>0$ で乗数を $e^{-4\pi^2(\eps+iat)|\xi|^2}$ に替える。Gauss 変換を右半平面で正則に延長すると逆変換は $(4\pi(\eps+iat))^{-n/2}e^{-|x|^2/(4(\eps+iat))}$ になる。正則な延長は、正の実数での等式と :ref[thm-cpx-identity] による。$L^2$ では $\eps\downarrow0$ で元の発展に収束する。空間核は全点で絶対値が $(4\pi a|t|)^{-n/2}$ 以下なので、$L^1$ データと積分して各点で極限を取れる。$L^2$ 極限の部分列は a.e. に収束するので、得られた積分表示が同じ解を表す。絶対値を取れば評価が従う。
:::

:::example[自由発展する Gauss 波束]
$f(x)=e^{-\pi x^2}$、$n=1$ なら $\widehat f=e^{-\pi\xi^2}$ なので

$$
U(t)f(x)=(1+4\pi iat)^{-1/2}e^{-\pi x^2/(1+4\pi iat)}.
$$

原点の絶対値は $(1+16\pi^2a^2t^2)^{-1/4}$ である。大きな $|t|$ では $|t|^{-1/2}$ の速さで減り、発展の分散率と一致する。
:::

補間により $2\le r\le\infty$ で $\|U(t)f\|_r\le(4\pi a|t|)^{-n(1/2-1/r)}\|f\|_{r'}$ を得る。実際 :ref[thm-meas-riesz-thorin] で $L^2\to L^2$ と $L^1\to L^\infty$ を指数 $1-2/r$ で補間する。

:::example[波束の広がりとノルムの保存]
$z=1+4\pi iat$ とすると上の Gauss 波束の二乗密度は $|z|^{-1}e^{-2\pi x^2/|z|^2}$ である。$x=|z|u$ の置換で積分は $1/\sqrt2$ に保たれる。例えば $t=1/(4\pi a)$ なら幅は $\sqrt2$ 倍になり、原点での二乗密度は $1/\sqrt2$ 倍になる。上限の低下と $L^2$ 保存が同時に起きることを確認できる。
:::

### 平行移動する包絡と群速度

Schrödinger 発展の位相は周波数の二乗である。中心周波数だけずらすと、二乗を展開した交差項が空間移動に変わる。このため、中心周波数を持つ波束の包絡が進む速さを正確に計算できる。

:::proposition[Galilei 変換]{#prop-four-schrodinger-galilean}
$f\in L^2(\R^n)$、$\eta\in\R^n$ とする。
$$
U(t)(e^{2\pi i x\cdot\eta}f)(x)=
e^{2\pi i x\cdot\eta-4\pi^2iat|\eta|^2}
(U(t)f)(x-4\pi at\eta)
$$
が $L^2$ で成り立つ。
:::

:::proof
変調した $f$ の変換は $\widehat f(\xi-\eta)$ である。逆変換の積分で $\zeta=\xi-\eta$ と置くと
$$
|\xi|^2=|\zeta|^2+2\zeta\cdot\eta+|\eta|^2 .
$$
交差項 $e^{-8\pi^2iat\zeta\cdot\eta}$ は逆変換の指数と合わさって位置 $x-4\pi at\eta$ を与え、$\eta$ だけの項は外に出る。まず Schwartz データで積分を実行し、一般の $L^2$ では移動、変調、$U(t)$ の全てが等長なので近似して等式を得る。
:::

:::example[波数と古典速度の照合]
一次元で $a=\hbar/(2m)$ とする。中心周波数 $\eta$ の平面波の運動量は $p=2\pi\hbar\eta$ であり、包絡の速度は $4\pi a\eta=p/m$ である。例えば $m=2$、$\hbar=1$ の単位を用い、$\eta=1/(2\pi)$ とすると速度は $1/2$ になる。初期包絡 $e^{-\pi x^2}$ の二乗密度は中心 $t/2$ で最大となり、その周りの幅は変調なしの Gauss 波束と同じである。電子を自由粒子とみなす式なら、外力や格子の効果が無視できる範囲でこの速度を読む。
:::

### 分散の位相を除いた大時間極限

:::proposition[拡大座標での Fourier 極限]{#prop-four-schrodinger-far-field}
$f\in L^1\cap L^2$、$t\ne0$ とし、$x=4\pi at\xi$ で観測する。そのとき
$$
(4\pi iat)^{n/2}e^{-i|4\pi at\xi|^2/(4at)} U(t)f(4\pi at\xi)
=\mathcal F\big(e^{i|\cdot|^2/(4at)}f\big)(\xi).
$$
右辺は $|t|\to\infty$ で $\widehat f$ へ一様かつ $L^2_\xi$ で収束する。さらに $\int |y|^2|f(y)|dy<\infty$ なら一様誤差は $(4a|t|)^{-1}\int |y|^2|f(y)|dy$ 以下である。
:::

:::proof
自由伝播の核の位相を $|x-y|^2=|x|^2-2x\cdot y+|y|^2$ と展開する。外に出る位相を左辺で相殺し、交差項に $x=4\pi at\xi$ を代入すると $e^{-2\pi i\xi\cdot y}$ になる。差は $\mathcal F((e^{i|\cdot|^2/(4at)}-1)f)$ である。括弧内は $2$ 以下で各点でゼロに近づくので、優収束で $L^1$ と $L^2$ の誤差が消える。変換の $L^1\to L^\infty$ 評価と Plancherel が二つの結論を与える。二次モーメントの仮定では $|e^{i\theta}-1|\le|\theta|$ を用いる。
:::

:::example[区間から放した波の零点]
$n=1$、$f=1_{[-1,1]}$ とする。拡大座標での極限は $\widehat f(\xi)=\sin(2\pi\xi)/(\pi\xi)$、原点では $2$ である。誤差は $(4a|t|)^{-1}\cdot2/3=1/(6a|t|)$ 以下である。したがって $\xi=1/2$、すなわち $x=2\pi at$ に沿うと、通常の分散の $|t|^{-1/2}$ に加えて $|t|^{-1}$ が掛かる。零点でない方向の振幅は $|t|^{-1/2}$ の主項を持つ。特定の空間方向を調べるには、上限ノルムだけよりもこの位相付きの極限が詳しい。
:::

### 波束の二次モーメント

:::proposition[自由発展の二次モーメント恒等式]{#prop-four-schrodinger-virial}
$f\in\mathcal S(\R^n)$ とし、$u=U(t)f$ と置く。このとき
$$
x_j u(t)=U(t)\bigl(x_j f-2iat\partial_j f\bigr),\qquad
\int |x|^2|u(t,x)|^2dx=\sum_{j=1}^n\|x_j f-2iat\partial_j f\|_2^2 .
$$
特に二次モーメントの時間二階微分は $8a^2\|\nabla f\|_2^2$ である。
:::

:::proof
$\mathcal F(x_j v)=(i/(2\pi))\partial_{\xi_j}\widehat v$ を用いる。発展の乗数を微分した項の係数は $(i/(2\pi))(-8\pi^2iat\xi_j)=4\pi at\xi_j$ であり、$-2iat\partial_j f$ の乗数と同じである。残る項は $U(t)x_jf$ である。二乗和の等式はユニタリ性で従う。右辺は二次式で、その $t^2$ の係数は $4a^2\sum_j\|\partial_j f\|_2^2$ なので二階微分を得る。
:::

:::example[初期チャープで一度絞られる幅]
$n=1$、$f(x)=e^{-\pi x^2}e^{-i b x^2}$、$b>0$ とする。$f'=-(2\pi+2ib)xf$ なので $xf-2iatf'=(1-4abt+4\pi iat)xf$ である。二次モーメントは
$$
\big((1-4abt)^2+16\pi^2a^2t^2\big)\|xf\|_2^2
$$
となる。$t=b/(4a(b^2+\pi^2))$ で最小になり、その倍率は $\pi^2/(b^2+\pi^2)$ である。例えば $b=\pi$ なら最小値は初期値の半分である。分散評価は長時間の減衰を制御するが、初期の位相に応じて一時的に波束が狭くなることもある。
:::

### 時空間積分に対する非端点評価

固定時刻の分散評価には $t=0$ で特異性がある。一方、時間変数にも積分ノルムを取ると、この特異性を時間の分数積分として扱える。

:::proposition[自由 Schrödinger 発展の非端点 Strichartz 評価]{#prop-four-schrodinger-strichartz}
$2<q<\infty$、$2<r<\infty$、$2/q+n/r=n/2$ とする。$f\in L^2(\R^n)$ に対し
$$
\|U(t)f\|_{L^q(\R_t;L^r(\R^n_x))}\le C_{n,a,q,r}\|f\|_2
$$
である。混合ノルムは、まず $x$ について $L^r$ を取り、その結果を $t$ について $L^q$ に入れる。
:::

:::proof
時間と空間にコンパクト台を持つ滑らかな $F(t,x)$ を取り、$A^*F=\int U(-s)F(s)ds$ と置く。ユニタリ性と Fubini で
$$
\|A^*F\|_2^2
=\int\!\!\int \inner{U(t-s)F(s)}{F(t)}\,ds\,dt.
$$
右辺の絶対値を空間の Hölder と分散補間で抑えると
$$
C\int \|F(t)\|_{r'}
      \left(\int |t-s|^{-\beta}\|F(s)\|_{r'}ds\right)dt,\qquad
\beta=n(1/2-1/r)=2/q
$$
以下である。時間方向の分数積分には、:ref[lem-four-fractional-integral] と同じ分割を入力指数 $p=q'$ で用いる。実際 $\alpha=1-\beta$ とし、$h\ge0$ に対して
$$
\int_{|t-s|\le R}|t-s|^{\alpha-1}h(s)ds\le C_\alpha R^\alpha Mh(t),\qquad
\int_{|t-s|>R}|t-s|^{\alpha-1}h(s)ds\le C_{\alpha,p}R^{\alpha-1/p}\|h\|_p
$$
である。第 1 式では半径 $2^{-j}R$ の環に分けて球平均を使い、第 2 式では Hölder を使う。後者の核の冪は無限遠で可積分であり、その条件 $\alpha<1/p$ は $1/p-\alpha=1/q>0$ から従う。$Mh(t)>0$ の点では $R=(\|h\|_p/Mh(t))^p$ と取ると、和は $C\|h\|_p^{\alpha p}(Mh(t))^{1-\alpha p}$ 以下である。$Mh(t)=0$ なら $h=0$ a.e. で自明である。$q(1-\alpha p)=p$ と Hardy–Littlewood 極大関数の $L^p$ 有界性を使うと、時間積分の作用素は $L^p\to L^q$ となる。極大関数の強型評価は :ref[thm-meas-maximal-weak] の弱型評価と $L^\infty$ 評価の補間 :ref[thm-meas-marcinkiewicz] から従う。

この評価を $h(s)=\|F(s)\|_{r'}$ に適用し、時間の Hölder を使って $\|A^*F\|_2^2\le C\|F\|_{L^{q'}L^{r'}}^2$ を得る。$\int\inner{U(t)f}{F(t)}dt=\inner{f}{A^*F}$ に双対ノルムを取ると主張が従う。左辺の関数と双対性で得た $L^q(L^r)$ の元は、滑らかなコンパクト台関数とのペアリングで一致する。
:::

:::example[一次元と二次元での時間指数]
一次元で $r=6$ とすると $n(1/2-1/r)=1/3$ であり $q=6$ となる。時間核は $|t|^{-1/3}$、分数積分の次数は $2/3$ である。二次元で $r=4$ とすると $q=4$、時間核は $|t|^{-1/2}$ となる。いずれも指数は $1$ より小さい。一方、$q=2$ をこの証明にそのまま入れると時間核 $|t|^{-1}$ は原点で可積分でなく、同じ絶対値による評価では扱えない。
:::

## Radon 変換と CT 再構成

平面内の密度を直線に沿って積分すると、方向ごとに一次元の信号を得る。この信号を Fourier 変換すると、元の二次元変換を原点を通る直線上で測った値になる。

:::definition[Radon 変換]{#def-four-radon}
$f\in\mathcal S(\R^2)$、$\theta\in S^1$ に対し

$$
\mathcal Rf(\theta,s)=\int_{\R} f(s\theta+r\theta^\perp)\,dr
$$

を :term[Radon 変換]{en="Radon transform" yomi="らどんへんかん"} という。$\theta^\perp$ は $\theta$ を $\pi/2$ 回転した単位ベクトルである。
:::

:::theorem[Fourier 切片とフィルタ付き逆投影]{#thm-four-radon-inversion}
$\mathcal F_s$ を $s$ だけの変換とすると、

$$
\mathcal F_s\mathcal Rf(\theta,\sigma)=\widehat f(\sigma\theta),\qquad
f(x)=\frac12\int_{S^1}\int_{\R}|\sigma|\mathcal F_s\mathcal Rf(\theta,\sigma)e^{2\pi i\sigma x\cdot\theta}\,d\sigma\,d\theta .
$$

$d\theta$ は全質量 $2\pi$ の弧長測度である。
:::

:::proof
切片の式は Fubini と直交変換 $x=s\theta+r\theta^\perp$ に従い、指数は $e^{-2\pi i\sigma x\cdot\theta}$ になる。二次元の反転積分を $\xi=\rho\theta$（$\rho>0$）に置くと Jacobian は $\rho$ である。$\sigma$ を両符号に拡げると $\xi$ を二回数えるため $1/2$ を掛ける。Schwartz 減衰により全ての積分は絶対収束する。
:::

周波数の重み $|\sigma|$ を掛けた後、$s=x\cdot\theta$ で評価して角度を積分するのが :term[フィルタ付き逆投影]{en="filtered back projection" yomi="ふぃるたつきぎゃくとうえい"} である。透過強度 $I$ が線減弱係数 $\mu$ に対して $dI/dr=-\mu I$ に従う理想モデルでは、対数比 $\log(I_0/I)=\int\mu\,dr$ を測る。$\mu$ の単位は逆長さである。角度と位置を連続に測る仮定の下で上式を使える。有限標本の処理では $|\sigma|$ が高周波のノイズも増幅するため、観測帯域で切る操作が再構成の分解能を変える。

:::example[Gauss 密度の投影]
$f(x)=e^{-\pi|x|^2}$ では $\mathcal Rf(\theta,s)=e^{-\pi s^2}$ であり、方向に依存しない。原点の逆投影は $\pi\int_\R|\sigma|e^{-\pi\sigma^2}d\sigma=1$ である。$\int_0^\infty\sigma e^{-\pi\sigma^2}d\sigma=1/(2\pi)$ を使った。
:::

### 逆投影の合成で現れる空間核

フィルタを掛ける前に全角度の投影を戻すと、元の画像に対する平滑化になる。この作用を空間積分で直接計算すると、フィルタが何を打ち消すかを確認できる。

:::proposition[Radon 変換の正規作用素]{#prop-four-radon-normal}
$\mathcal R^*g(x)=\int_{S^1}g(\theta,x\cdot\theta)d\theta$ とすると、Schwartz 関数について
$$
\mathcal R^*\mathcal Rf(x)=2\int_{\R^2}\frac{f(x-y)}{|y|}\,dy .
$$
右辺は $x$ ごとに絶対収束する。
:::

:::proof
$\mathcal Rf(\theta,x\cdot\theta)=\int_\R f(x+r\theta^\perp)dr$ である。元の定義の積分変数を $x$ の $\theta^\perp$ 成分だけ移した。正負の $r$ は角度全体で同じ寄与を持つため、$r>0$ だけを二倍する。$y=-r\theta^\perp$ の極座標では $dy=r\,dr\,d\theta$ であり、元の積分の $dr\,d\theta$ との比が $1/|y|$ になる。局所ではこの冪が可積分で、遠方では $f$ の急減少で収束する。
:::

:::example[原点での逆投影]
$f(x)=e^{-\pi|x|^2}$ では $\mathcal R f(\theta,0)=1$ なので、フィルタなしの逆投影の原点は $2\pi$ である。空間核の式でも $2\cdot 2\pi\int_0^\infty e^{-\pi r^2}dr=2\pi$ となる。フィルタ付きの場合は :ref[thm-four-radon-inversion] で原点は元の値 $1$ に戻る。
:::

### 高次元の超平面積分

平面で直線に沿って積分したときと同じく、高次元でも法線方向の Fourier 変換が切片になる。ただし極座標の体積要素が変わるため、逆投影に必要な微分の次数も変わる。

:::proposition[超平面 Radon 変換と反転]{#prop-four-radon-hyperplane}
$n\ge2$、$f\in\mathcal S(\R^n)$ とし、$\mathcal R_n f(\theta,s)=\int_{\theta^\perp} f(s\theta+y)dy$ と定める。ここで $dy$ は超平面上の Euclid 測度である。$g=\mathcal R_n f$ とすると
$$
\mathcal F_sg(\theta,\sigma)=\widehat f(\sigma\theta),\qquad
f(x)=\frac12\int_{S^{n-1}}\int_\R
|\sigma|^{n-1}\mathcal F_sg(\theta,\sigma)e^{2\pi i\sigma x\cdot\theta}\,d\sigma\,d\theta .
$$
さらに
$$
\int_{S^{n-1}}\int_{\R}|\sigma|^{n-1}|\mathcal F_s g(\theta,\sigma)|^2\,d\sigma\,d\theta
=2\|f\|_2^2 .
$$
奇数次元 $n=2m+1$ では、反転積分の内側は $(-1)^m(2\pi)^{-2m}\partial_s^{2m}g(\theta,x\cdot\theta)$ に等しい。$d\theta$ は正規化していない球面測度とする。
:::

:::proof
$x=s\theta+y$ は直交分解で、$dx=ds\,dy$ である。$x\cdot\theta=s$ だから、$s$ の変換の指数を超平面積分と合わせれば切片の式になる。反転積分に極座標 $\xi=\rho\theta$、$\rho>0$ を入れると体積要素は $\rho^{n-1}d\rho\,d\theta$ である。$\sigma$ の負側も使うと、$\theta\mapsto-\theta$ で同じ積分となる。従って係数は $1/2$ である。二乗積分も切片の式から同じ極座標積分に戻り、Plancherel でノルム等式を得る。最後に、$n=2m+1$ では $(2\pi i\sigma)^{2m}=(-1)^m(2\pi)^{2m}|\sigma|^{n-1}$ なので、逆変換で整数階の微分に置き換えられる。積分と微分に必要な周波数の減衰は $\widehat f$ の急減少性から従う。
:::

:::example[三次元の Gauss 密度を二階微分で戻す]
$f(x)=e^{-\pi|x|^2}$ なら直交する二方向で Gauss 積分の質量はともに $1$ なので、$g(\theta,s)=e^{-\pi s^2}$ となる。従って $\partial_s^2g(\theta,0)=-2\pi$ であり、原点の再構成は $-\frac1{8\pi^2}\int_{S^2}(-2\pi)\,d\theta=1$ である。二次元の重みは $|\sigma|$ であるため、同じ整数階の微分による置き換えはできない。
:::

### 投影の対称性とモーメント

同じ直線は $(\theta,s)$ と $(-\theta,-s)$ で表される。したがって任意の方向に独立な信号を置いても、画像の投影になるとは限らない。方向を変えても保たれる積分量から、観測データが満たす条件を得られる。

:::proposition[投影の多項式モーメント]{#prop-four-radon-moments}
$f\in\mathcal S(\R^2)$、$k$ は非負整数とし、$g=\mathcal R f$ とする。$g(-\theta,-s)=g(\theta,s)$ であり、
$$
\int_{\R}s^k g(\theta,s)\,ds=\int_{\R^2}(x\cdot\theta)^k f(x)dx
$$
は $\theta$ を変数とする次数 $k$ の斉次多項式を単位円に制限したものである。
:::

:::proof
方向を反転すると $\theta^\perp$ も反転するので、線積分で $r\mapsto-r$ と置けば対称性を得る。モーメントについては、$x=s\theta+r\theta^\perp$ の直交座標で $x\cdot\theta=s$ を使う。Fubini は $\int |x|^k|f(x)|dx<\infty$ から適用できる。$k$ 乗を二項展開すると $\theta_1^j\theta_2^{k-j}$ の係数は $\binom kj\int x_1^jx_2^{k-j}f(x)dx$ である。
:::

:::example[楕円形 Gauss 密度の投影]
$f(x_1,x_2)=e^{-\pi(x_1^2/4+x_2^2)}$ とする。変換は $\widehat f(\xi)=2e^{-\pi(4\xi_1^2+\xi_2^2)}$ なので、切片定理と一次元 Gauss の反転から
$$
g(\theta,s)=\frac{2}{w(\theta)}e^{-\pi s^2/w(\theta)^2},\qquad
w(\theta)=\sqrt{4\theta_1^2+\theta_2^2}.
$$
$\theta=(1,0)$ の原点値は $1$、$\theta=(0,1)$ の原点値は $2$ だが、どちらも総質量は $2$ である。二次モーメントは $2 w(\theta)^2/(2\pi)=w(\theta)^2/\pi$ であり、角度とともに二次多項式として変化する。
:::

:::example[対称性だけを満たす偽の投影]
$g(\theta,s)=(1+\theta_1^2-\theta_2^2)e^{-\pi s^2}$ と置く。$s$ が偶で $\theta$ も同時反転で値が変わらないから、$g(-\theta,-s)=g(\theta,s)$ である。しかし $s$ 積分は $1+\theta_1^2-\theta_2^2$ であり、方向 $(1,0)$ で $2$、$(0,1)$ で $0$ である。$k=0$ の条件に反するので、この関数を全方向の投影に持つ Schwartz 画像は存在しない。
:::

### 放射対称な画像の Abel 型反転

放射対称な画像では方向が一つで足りる。直線内を進むと原点からの距離が $r=\sqrt{s^2+t^2}$ と変わる。$t$ を $r$ に換えると、接点 $r=|s|$ 付近で距離の変化が遅くなることが平方根の核に現れる。

:::proposition[放射対称な投影と反転]{#prop-four-radon-abel}
$f(x)=F(|x|)$ を放射対称な Schwartz 関数とし、$g(s)=\mathcal R f(\theta,s)$ とする。$g$ は方向によらず偶関数であり、
$$
g(s)=2\int_{|s|}^\infty \frac{F(r)r}{\sqrt{r^2-s^2}}dr .
$$
さらに $r>0$ に対して
$$
F(r)=-\frac1\pi\int_r^\infty\frac{g'(s)}{\sqrt{s^2-r^2}}ds .
$$
:::

:::proof
直線の変数 $t$ の正負から二倍し、$r=\sqrt{s^2+t^2}$、$dt=r(r^2-s^2)^{-1/2}dr$ と置けば第一式を得る。$v=s^2$、$w=r^2$、$G(v)=g(\sqrt v)$、$\Phi(w)=F(\sqrt w)$ と置くと
$$
G(v)=\int_v^\infty \frac{\Phi(w)}{\sqrt{w-v}}dw
$$
である。$u>0$ から $\int_u^\infty (v-u)^{-1/2}G(v)dv$ を取り、積分を交換する。絶対値を入れても、内側の
$$
\int_u^w \frac{dv}{\sqrt{(v-u)(w-v)}}=\int_0^1\frac{dt}{\sqrt{t(1-t)}}=\pi
$$
が有限なので Fubini を適用できる。最後の値は $t=\sin^2\alpha$ の置換で得られる。従って二度重ねた積分は $\pi\int_u^\infty\Phi(w)dw$ に等しい。左辺で $v=u+t$ としてから $u$ で微分すると $\int_0^\infty t^{-1/2}G'(u+t)dt$ になる。$u>0$ 固定の近傍で $G'$ は有界であり、無限遠では急減少なので優収束を使える。右辺の微分は $-\pi\Phi(u)$ である。$u=r^2$ とし $v=s^2$ へ戻すと反転式を得る。
:::

:::example[一様円板が作る弦の長さ]
$f=1_{\{|x|<R\}}$ なら直線の円板内での長さは $g(s)=2\sqrt{R^2-s^2}$（$|s|<R$）である。従って $R=2$、$s=1$ なら $g(1)=2\sqrt3$ である。滑らかさの仮定の外でも、この例では反転積分を直接計算できる。$0<r<R$ で $-g'(s)=2s/\sqrt{R^2-s^2}$（$r<s<R$）だから
$$
-\frac1\pi\int_r^\infty\frac{g'(s)}{\sqrt{s^2-r^2}}ds
=\frac1\pi\int_{r^2}^{R^2}\frac{dv}{\sqrt{(R^2-v)(v-r^2)}}=1.
$$
$r>R$ ではゼロである。端点 $R$ で $g'$ は発散するが絶対可積分であり、$g$ 自身は連続であるため、弱微分に点質量は現れない。
:::

### 帯域を制限した再構成と誤差

高周波に掛ける $|\sigma|$ は、逆変換の正しさに必要な Jacobian である。同時に、投影の $L^2$ 誤差から再構成の $L^2$ 誤差へ移るときには重みが現れる。周波数帯域を有限にしておくと、この増幅を定量的に抑えられる。

:::proposition[有限帯域でのノイズ評価]{#prop-four-radon-bandlimited}
$B>0$ とし、反転積分で $|\sigma|\le B$ に制限した結果を $I_B g$ と書く。$f\in L^1\cap L^2(\R^2)$ とし、a.e. に定まる線積分を $g=\mathcal Rf$ とする。このとき $g\in L^2(S^1\times\R)$ であり、$\widehat{I_B g}=1_{\{|\xi|\le B\}}\widehat f$ である。投影に誤差 $h\in L^2(S^1\times\R)$ を加えたとき
$$
\|I_B h\|_{L^2(\R^2)}\le\sqrt{B/2}\,\|h\|_{L^2(S^1\times\R)}.
$$
更に $f\in H^s$、$s>0$ で $\|f\|_{\dot H^s}=\||\xi|^s\widehat f\|_2$ とすれば
$$
\|I_B(g+h)-f\|_2\le B^{-s}\|f\|_{\dot H^s}+\sqrt{B/2}\|h\|_2.
$$
:::

:::proof
$\xi=\rho\theta$、$\rho>0$ とする。正負の $\sigma$ を別に分け、負の項で $\theta\mapsto-\theta$ と置くと、反転積分は
$$
\widehat{I_B h}(\rho\theta)
=1_{\{\rho\le B\}}\tfrac12\bigl(\mathcal F_s h(\theta,\rho)+\mathcal F_s h(-\theta,-\rho)\bigr)
$$
の二次元逆変換である。まず滑らかなデータに対して計算し、一般の $L^2$ データでは以下の評価で定義を延長する。$f\in L^1$ なので直交座標と Fubini により $g(\theta,\cdot)\in L^1$ であり、切片定理と同じ積分計算ができる。$|\sigma|\le1$ では $\widehat f$ の有界性を使い、$|\sigma|>1$ では極座標で
$$
\int_{S^1}\int_{|\sigma|>1}|\widehat f(\sigma\theta)|^2d\sigma\,d\theta
=2\int_{|\xi|>1}|\widehat f(\xi)|^2|\xi|^{-1}d\xi
\le2\|f\|_2^2
$$
となる。従って切片は二乗可積分である。一次元の Plancherel の逆向きと変換の一意性により $g\in L^2(S^1\times\R)$ が従う。より具体的には、$L^2$ で得た逆変換と $L^1$ で定義した $g$ は緩増加超関数として変換が一致するため等しい。上の式の $g$ の二項はともに $\widehat f(\rho\theta)$ である。誤差では $|(a+b)/2|^2\le(|a|^2+|b|^2)/2$ と Jacobian $\rho$ を使って
$$
\|\widehat{I_B h}\|_2^2
\le B/2\int_{S^1}\int_{\R}|\mathcal F_s h(\theta,\sigma)|^2d\sigma\,d\theta
= B/2\|h\|_2^2
$$
を得る。最後に $|\xi|>B$ で $1\le B^{-2s}|\xi|^{2s}$ として切断誤差の二乗積分を抑え、誤差の三角不等式を適用する。
:::

:::example[Gauss 密度を途中まで戻す]
$f(x)=e^{-\pi|x|^2}$ なら原点で $I_B\mathcal R f(0)=2\pi\int_0^B\rho e^{-\pi\rho^2}d\rho=1-e^{-\pi B^2}$ である。$B=1$ では原点の欠損は $e^{-\pi}$ である。一方、誤差ノルムが $\|h\|_2=10^{-2}$ ならその再構成の $L^2$ ノルムは $B=1$ で $10^{-2}/\sqrt2$ 以下である。前者の計算は一点での値、後者は平面全体での二乗平均なので、比較するときは測っている量をそろえる必要がある。
:::

### 観測しない角度があるとき

ある方向を観測しても、その方向と直交する変化を直接測ったことにはならない。切片定理で見ると、観測方向が不足したときに欠落するのは周波数平面内の扇形である。画像がコンパクト台を持つかどうかで、一意性の結論が変わる。

:::proposition[角度制限と一意性]{#prop-four-radon-limited-angles}
$A\subset S^1$ を空でない開弧とする。$f\in C_c^\infty(\R^2)$ について $\mathcal R f(\theta,s)=0$ が全 $\theta\in A$、$s\in\R$ で成り立てば $f=0$ である。一方、$S^1\setminus(\overline A\cup(-\overline A))$ に開弧が残る場合、同じ観測を全てゼロにする非零の Schwartz 画像が存在する。
:::

:::proof
第一の仮定では切片定理により $\widehat f$ が非空の開錐上でゼロになる。その錐内に開長方形 $I_1\times I_2$ を取る。コンパクト台のため、$\widehat f$ の積分で一方の周波数を複素変数にしても、各コンパクトな複素集合上では全ての微分を可積分な関数で支配できる。従って一方の変数を固定した各関数は整関数である。まず $\xi_2\in I_2$ を固定し、第 1 変数に :ref[thm-cpx-identity] を適用すると $\widehat f(\xi_1,\xi_2)=0$ が全実数 $\xi_1$ に拡がる。今度は任意の $\xi_1$ を固定し第 2 変数に同じ原理を使う。これで $\widehat f=0$ が平面全体に拡がり、反転により $f=0$ である。

第二の仮定では、欠落する開錐内で原点から離れた小球を選び、そこに台を持つ非零の $C_c^\infty$ 関数 $v$ を置く。$\widehat f=v$ として反転すれば $f\in\mathcal S$、$f\ne0$ である。観測領域の切片は全て台を外れるのでゼロであり、一次元の反転から各投影もゼロである。実数値の反例が必要なら $v$ を実非負に選び、$v(\xi)+v(-\xi)$ を使う。
:::

:::example[水平に近い角度の欠落領域]
観測する単位方向が $|\theta_2|<1/2$ を満たすとする。この方向とその反対向きにある周波数は $|\xi_2|<|\xi|/2$ に入る。$v$ の台を $B((0,2),1/4)$ に置くと $|\xi_2|\ge7/4$、$|\xi|\le9/4$ なので、比は $7/9>1/2$ である。従って全観測投影がゼロとなる。方向数を増やさずに各直線上の標本数だけを増やしても、この扇形の欠落は埋まらない。
:::

## 演習問題

::::exercise{level=basic}
円周上で $\kappa=2$、初期値 $\cos(3\theta)$ の熱解を求めよ。
:::solution
両周波数 $\pm3$ の乗数は $e^{-2\cdot9t}$ なので解は $e^{-18t}\cos(3\theta)$ である。
:::
::::

::::exercise{level=basic}
一次元の波動方程式で $c=1,f=0,g(x)=e^{-x^2}$ とする。$\partial_tu(t,0)$ を求めよ。
:::solution
$u(t,0)=\frac12\int_{-t}^{t}e^{-y^2}dy$ である。積分の上下端を微分すると $(e^{-t^2}+e^{-(-t)^2})/2=e^{-t^2}$ となる。
:::
::::

::::exercise{level=basic}
上半平面で境界値 $2+\cos(4\pi x)$ の有界な調和拡張を求めよ。
:::solution
定数は質量 $1$ の核で保存される。余弦の周波数は $\pm2$ なので減衰率は $4\pi y$ である。よって $2+e^{-4\pi y}\cos(4\pi x)$ である。
:::
::::

::::exercise{level=basic}
二次元で $f(x)=e^{-2\pi|x|^2}$ の Radon 変換を求めよ。
:::solution
直交座標では $|x|^2=s^2+r^2$ となる。$r$ に関する Gauss 積分は $1/\sqrt2$ なので $\mathcal Rf(\theta,s)=2^{-1/2}e^{-2\pi s^2}$ である。
:::
::::

::::exercise{level=advanced}
円周の熱方程式で平均ゼロの $L^2$ データに対し $\|u(t)\|_2\le e^{-\kappa t}\|f\|_2$ を証明し、等号条件を求めよ。
:::solution
第 $0$ 係数はゼロで、他の二乗減衰率は $e^{-2\kappa k^2t}\le e^{-2\kappa t}$ である。Parseval に入れて平方根を取る。$t>0$ で等号は $|k|\ge2$ の係数が全てゼロの場合に限る。
:::
::::

::::exercise{level=advanced}
波動方程式の周波数解から全エネルギー保存を示せ。
:::solution
固定 $\xi$ で $v=\widehat u$、$\omega=2\pi c|\xi|$ とすると $v''=-\omega^2v$ である。$\frac{d}{dt}(|v'|^2+\omega^2|v|^2)=2\Re(v''\overline{v'}+\omega^2v'\bar v)=0$ である。積分して Plancherel を使えば勾配と時間微分のエネルギーに戻る。
:::
::::

::::exercise{level=advanced}
$t\ne0$ を固定する。Schrödinger の分散評価の定数が鋭いことを示せ。
:::solution
非負で非零の $\phi\in\mathcal D$ を取り $f(y)=e^{-i|y|^2/(4at)}\phi(y)$ とする。$x=0$ では核の位相と相殺し、$|U(t)f(0)|=(4\pi a|t|)^{-n/2}\int\phi$ である。解は連続なので本質上限はこの値以上であり、$\|f\|_1=\int\phi$ で等号が従う。
:::
::::

::::exercise{level=advanced}
$\mathcal Rf=0$ なら Schwartz 関数 $f$ はゼロであることを、逆投影公式を使わずに示せ。
:::solution
切片の式から $\widehat f(\sigma\theta)=0$ である。任意の $\xi\ne0$ を $(\sigma,\theta)=(|\xi|,\xi/|\xi|)$ で表し、原点は $\sigma=0$ とすれば変換は全点でゼロである。Fourier 変換の一意性で $f=0$ となる。
:::
::::

::::exercise{level=advanced}
$n=1$、$\kappa>0$、$f(x)=e^{-\pi x^2}$ とする。熱方程式の解の中心値と二次モーメント $\int x^2u(t,x)dx$ を求めよ。
:::solution
$f=G_{1/(4\pi)}$ なので $u=G_{\kappa t+1/(4\pi)}$ である。中心値は $(1+4\pi\kappa t)^{-1/2}$ となる。質量は $1$ であり、Gauss 核 $G_s$ の二次モーメントは $2s$ だから、ここでは $2\kappa t+1/(2\pi)$ である。
:::
::::

::::exercise{level=advanced}
幅 $L=2$ の帯で $u(x,0)=\sin(\pi x)$、$u(x,2)=\sin(\pi x)$ を満たす有界調和解を求めよ。
:::solution
$u(x,y)=a(y)\sin(\pi x)$ と置けば $a''=\pi^2a$ である。上下の値が同じなので中心の反射に関して対称な $a(y)=\cosh(\pi(y-1))/\cosh\pi$ が両端で $1$ となる。直接代入して調和性と有界性を確認でき、:ref[prop-four-laplace-strip] の最大値原理による一意性はこの有界なデータでも使える。
:::
::::

## 確認テスト

::::quiz{answer="0.25" tol="1e-6"}
$\kappa=1$、初期値 $\cos\theta$ の円周熱解の $\theta=0,t=\log4$ における値を求めよ。
:::explanation
減衰因子は $e^{-t}$ で、空間因子は $1$ である。従って $e^{-\log4}=1/4$ となる。
:::
::::

::::quiz
三次元の波動方程式で $f=0$ のとき、位置 $x$ の時刻 $t$ の値を直接決める $g$ の範囲はどれか。

- [ ] $x$ を中心とする半径 $ct$ の球の外部
- [x] $x$ を中心とする半径 $ct$ の球面
- [ ] 半径 $ct$ の開球内部だけ

:::explanation
:ref[thm-four-wave-means] で $u=tM_{ct}g$ となる。この平均の積分領域は球面である。例えば球の内部だけに台を持つ $g$ はその時刻の値に寄与しない。
:::
::::

::::quiz{answer="0.5" tol="1e-6"}
境界が $1_{[-1,1]}$ の Poisson 積分を $(0,1)$ で評価せよ。
:::explanation
$\pi^{-1}\int_{-1}^{1}(1+s^2)^{-1}ds=(\arctan1-\arctan(-1))/\pi=1/2$ である。
:::
::::

::::quiz
二次元のフィルタ付き逆投影で使う一次元周波数の重みはどれか。

- [ ] $1$
- [x] $|\sigma|$
- [ ] $\sigma^2$

:::explanation
二次元の反転積分を極座標にすると Jacobian は半径である。切片の周波数 $\sigma$ を両符号で扱うと半径は $|\sigma|$ になる。全角度での二重計数は外側の $1/2$ が補う。
:::
::::

---
title: ポテンシャル論と Dirichlet 問題
description: Newton ポテンシャルの Hölder 評価と Perron 法を導き、層ポテンシャルで境界条件を積分方程式に変える。
notation:
  - symbol: '[f]_\alpha'
    meaning: 指数 α の Hölder 半ノルム
  - symbol: 'S\mu'
    meaning: 一重層ポテンシャル
  - symbol: 'D\mu'
    meaning: 二重層ポテンシャル
---

## この章の目標

- 基本解の二階微分に現れる特異性を定数の相殺で扱う。
- 劣調和関数の上限から調和関数を作り、境界値の実現を障壁で示す。
- 層ポテンシャルの跳躍の向きと大きさを求める。
- Dirichlet 問題と Neumann 問題で積分方程式の適合条件を区別する。

## Poisson 方程式と Hölder 連続性

$0<\alpha<1$ に対し $[f]_\alpha=\sup_{x\ne y}|f(x)-f(y)|/|x-y|^\alpha<\infty$ を :term[Hölder 連続]{en="Hölder continuity" yomi="へるだーれんぞく"} という。Newton ポテンシャル $u=\Gamma*f$ を二階微分すると核の大きさが $|z|^{-n}$ になる。球面上で平均が消えることを使えば $f(x-z)$ を $f(x-z)-f(x)$ に置き換えられるので、積分が収束する。

:::proposition[Newton ポテンシャルの二階微分]{#prop-pde-newton-holder}
$f\in C_c^{0,\alpha}(\R^n)$、$n\ge2$ なら $u=\Gamma*f$ は $C^{2,\alpha}_{\rm loc}$ であり、$-\Delta u=f$ を満たす。任意の $R>0$ について、

$$
u_{ij}(x)=\int_{|z|<R}K_{ij}(z)(f(x-z)-f(x))\,\dd z
+\int_{|z|\ge R}K_{ij}(z)f(x-z)\,\dd z-\frac{\delta_{ij}}{n}f(x),
\label{eq-pde-newton-hessian}
$$

ここで $K_{ij}=(n z_i z_j/|z|^2-\delta_{ij})/(\sigma_{n-1}|z|^n)$ である。
:::

:::proof
:ref[prop-pde-laplace-fundamental] の $\partial_i\Gamma=-z_i/(\sigma_{n-1}|z|^n)$ を原点の外で微分すると $K_{ij}$ になる。球面上で $\theta_i\theta_j$ の平均は $\delta_{ij}/n$ だから $K_{ij}$ の環状領域上での積分は $0$ である。$\partial_i\Gamma$ を試験関数に対して部分積分すると、小球境界から $-\delta_{ij}\phi(0)/n$ が残る。よって滑らかな $f$ では :eqref[eq-pde-newton-hessian] が従う。$|z|<r$ の寄与は $C[f]_\alpha r^\alpha$ 以下である。$f$ を平滑化すると Hölder 半ノルムは増えないので、まず $r$ の外で一様収束を用い、次に $r\downarrow0$ とすれば二階微分は局所一様に右辺へ収束する。一階以下は核の局所可積分性で収束するから右辺は古典的二階微分である。跡を取ると $\sum_iK_{ii}=0$ なので $-\Delta u=f$ である。

Hölder 評価では主値部分を $x,y$ で比較する。$d=|x-y|$ とし、両積分を中心 $x$ の球 $B_{2d}(x)$ と外側に分ける。各中心で定数 $f(x),f(y)$ を引くと、内側の部分は $C[f]_\alpha d^\alpha$ 以下である。$y$ を中心にした場合も積分領域に $B_d(y)$ が含まれ $B_{3d}(y)$ に含まれるため、残る定数の主値積分は絶対値が $C$ 以下である。より明示的には両方から $f(x)$ を引いてよい。全空間の核の定数に対する主値を中心ごとの球で切れば $0$ であり、$y$ 側の内側には $|f(y)-f(x)|$ と有界な主値係数の積が加わるだけである。外側の差は $|K(x-z)-K(y-z)|\le Cd|x-z|^{-n-1}$ と $|f(z)-f(x)|\le [f]_\alpha|x-z|^\alpha$ より、$Cd[f]_\alpha\int_{2d}^\infty r^{\alpha-2}\,\dd r\le C_\alpha[f]_\alpha d^\alpha$ である。異なる中心の無限遠切断の差は有界な $f$ に対して $O(d/R)$ で消える。最後の局所項も Hölder なので結論を得る。
:::

:::example[右辺の指数がそのまま残る場合]
一次元で $-u''=|x|^{1/2}$、$u(0)=u'(0)=0$ なら $u=-4|x|^{5/2}/15$ である。$u''(x)=-|x|^{1/2}$ だから、原点では $u''$ は $1/2$ より大きい任意の指数について Hölder 連続でない。
:::

:::example[二次元の相殺]
$n=2$ における $K_{11}$ は極座標で $\cos(2\theta)/(2\pi r^2)$ である。角度積分は $0$ であり、$f$ が近傍で一定なら、その近傍から二階微分に寄与するのは局所項 $-f(x)/2$ のみである。
:::

### Hessian の相殺と Hölder 指数

Newton 核の Hessian の通常の部分を $K$ とし、$|K(z)|\le C|z|^{-n}$、$|DK(z)|\le C|z|^{-n-1}$ とする。$K$ の球面平均は零である。近い二点 $x,x'$ の距離を $d$ とする。二点を中心とする半径 $2d$ の近傍では、$f(y)-f(x)$ または $f(y)-f(x')$ を掛ければ $C[f]_\alpha\int_0^{3d}r^{\alpha-1}\dd r$ 以下である。遠方では核の差に平均値の定理を使い、積を
$$
C[f]_\alpha d\int_{2d}^{\infty}r^{\alpha-2}\dd r
\le \frac{C[f]_\alpha}{1-\alpha}d^\alpha
$$
で抑える。ここで引いた定数の違いは核の平均零を用いる。二点ともに同じ中心の環状領域を選び直すと、動いた切断面の帯は原点側で $2d\le r\le3d$、外側で幅 $d$ であり、そこで残る項も $C[f]_\alpha d^\alpha$ 以下になる。従って $0<\alpha<1$ が指数を保つ範囲である。$\alpha=1$ では遠方積分を同じ方法で有限距離まで積分しても $d\log(1/d)$ が残るので、この計算から Lipschitz の結論は得られない。

### Poisson 方程式の内部 Schauder 評価

$-\Delta u=f$ を $B_{2R}(x_0)$ の分布方程式とし、$u$ は有界、$f\in C^{0,\alpha}(B_{2R}(x_0))$ とする。このとき $u$ は $B_R(x_0)$ で $C^{2,\alpha}$ の代表を持ち、
$$
R^2\norm{D^2u}_{L^\infty(B_R)}
+R^{2+\alpha}[D^2u]_{\alpha;B_R}
\le C_{n,\alpha}\big(\norm{u}_{L^\infty(B_{2R})}
+R^2\norm{f}_{L^\infty(B_{2R})}+R^{2+\alpha}[f]_{\alpha;B_{2R}}\big)
$$
が成り立つ。ここで $[\,\cdot\,]_{\alpha;B}$ は二点を $B$ に限った Hölder 半ノルムである。

まず $R=1$、$x_0=0$ とする。$B_{3/2}$ で $1$、$B_2$ の内側に台を持つ滑らかな切断 $\chi$ を選ぶ。$v=N(\chi f)$ について、核の局所可積分性から $\norm v_{L^\infty(B_{3/2})}\le C\norm f_\infty$ である。:eqref[eq-pde-newton-hessian] では近くの積分に Hölder 半ノルム、遠くの積分に有界性と台の大きさを使い、$\norm{D^2v}_{C^{0,\alpha}(B_1)}\le C\norm f_{C^{0,\alpha}(B_2)}$ を得る。一方 $u-v$ は $B_{3/2}$ で分布調和なので :ref[prop-pde-weyl-harmonic] を適用できる。:ref[prop-pde-harmonic-derivative] で二階と三階の微分を評価すれば、球 $B_1$ は凸なので二階微分の Hölder 半ノルムも評価できる。両方の評価を足し、$x=x_0+Rz$ に戻すと上の式になる。

### Lipschitz な源に対する端点の反例

平面で $u=r^3\log r\cos(3\theta)$、$u(0)=0$ とする。$P_3(x,y)=x^3-3xy^2$ と書けば $u=P_3\log r$ であり、$|D^ku|=O(r^{3-k}(1+|\log r|))$（$k\le2$）なので、原点まで $C^2$ である。$r>0$ で極座標の Laplacian を使うと、$R=r^3\log r$ に対し
$$
R^{\prime\prime}+r^{-1}R^{\prime}-9r^{-2}R=6r .
$$
従って $\Delta u=6r\cos(3\theta)$ であり、原点でも連続に零へ延びる。右辺は半径微分が $6\cos3\theta$、単位接方向の微分が $-18\sin3\theta$ なので、原点を除いた平面で勾配が有界である。線分が原点を通る場合は二つに分けて積分すれば、原点を含む円板で Lipschitz になる。一方、正の $x$ 軸上で
$$
u_{xx}(x,0)=(x^3\log x)^{\prime\prime}=6x\log x+5x .
$$
$u_{xx}(0)=0$ からの差を $x$ で割ると $6\log x+5$ なので有界にならない。従って Poisson 方程式の内部 Schauder 評価で $\alpha=1$ をそのまま含められない。$\chi u$ を原点付近で $1$ の滑らかな切断でコンパクトに延長すると、源もコンパクト台で Lipschitz となり、局所で同じ反例が得られる。実際 $(\chi u)-N(-\Delta(\chi u))$ は分布調和であり、その滑らかな二階微分を足しても上の対数発散は消えない。

### Green 関数で体積源の符号を確認する

二次元の上半平面で $a=(0,1)$、反射点を $a^*=(0,-1)$ とする。関数
$$
U(x,y)=\frac1{2\pi}\log\frac{\sqrt{x^2+(y+1)^2}}{\sqrt{x^2+(y-1)^2}}
$$
は $a$ における基本解と、領域外の $a^*$ における基本解の差である。従って上半平面では $-\Delta U=\delta_a$、境界では両距離が同じなので $U(x,0)=0$。内部では $y>0$ なので分子の距離の方が長く $U>0$ である。例えば $(0,2)$ では距離が $3,1$ で $U(0,2)=\log3/(2\pi)$。法線を $-e_y$ とすると
$$
-\partial_\nu U(x,0)=\partial_y U(x,0)
=\frac1{\pi(1+x^2)}
$$
であり、境界上の積分は $[\arctan x]_{-\infty}^{\infty}/\pi=1$ となる。内点の点源を境界の Poisson 核へ移すと、分布の総量が境界における規格化に現れる。

:::example[点源を平滑な密度へ置き換える]
$0<\eps<1/2$ とし、$B_\eps(a)$ に台を持つ滑らかな非負の密度 $f_\eps$ で $\int f_\eps=1$ とする。各点 $z$ の反射点 $z^*$ を用いた上の差を $f_\eps(z)$ に対して積分する。これは滑らかな体積源 $f_\eps$ を持つ零境界解である。点 $(0,2)$ は密度の台から距離 $1/2$ 以上離れるため、Green 核は $z$ について連続で、積分値は $\eps\downarrow0$ で $U(0,2)$ に収束する。点源に向かう極限でも観測点が特異点から離れていれば、核の連続性を直接使える。
:::

## 劣調和関数と Perron の方法

:::definition[劣調和関数]{#def-pde-subharmonic}
連続関数 $v$ が各 $\overline B\subset\Omega$ で、その球の調和関数 $h$ と $v|_{\partial B}\le h|_{\partial B}$ なら $v\le h$ を満たすとき :term[劣調和関数]{en="subharmonic function" yomi="れつちょうわかんすう"} という。$-v$ が劣調和なら $v$ は優調和である。
:::

Poisson 積分を比較対象に選ぶと $v$ は球面平均以下である。逆に球面での平均不等式があれば $v-h$ にも平均不等式がある。内部最大点の開閉の議論を :ref[thm-pde-harmonic-maximum] と同様に適用して比較が従う。$C^2$ 関数では平均の Taylor 展開で $\Delta v\ge0$ が必要であり、逆に球面平均の微分が非負になるので十分である。二つの劣調和関数の最大は比較の定義から劣調和である。さらに球内で $v$ を境界値 $v|_{\partial B}$ の調和関数に置き換えても劣調和である。これを調和持ち上げと呼ぶ。実際、別の球と調和関数で境界から比較する際、持ち上げ前の $v$ はその球全体で比較対象以下である。二球の共通部分では持ち上げ後の関数も境界全体で比較対象以下であり、調和関数の差に最大値原理を使える。

:::theorem[Perron の包絡線]{#thm-pde-perron}
$\Omega$ は有界領域、$g\in C(\partial\Omega)$ とする。境界で $\limsup v\le g$ を満たす連続な有界劣調和関数全体を $\mathcal S_g$ とする。$U(x)=\sup_{v\in\mathcal S_g}v(x)$ は $\Omega$ で調和である。
:::

:::proof
$\min g$ が候補であり、境界の上極限条件と内部最大値の議論で各候補は $\max g$ 以下である。閉包が $\Omega$ に含まれる球 $B$ と $x_0\in B$ を固定し、$v_j(x_0)\to U(x_0)$ となる候補を取る。$h_0=\min g$ とし、$\max(h_{j-1},v_j)$ を $B$ で調和に持ち上げた関数を $\Omega$ 全体で $h_j$ と書く。持ち上げは元の候補以上であり、境界上極限を変えないので、$h_j\in\mathcal S_g$ は増加する。球内では調和であり、:ref[prop-pde-harmonic-derivative] による等連続性を用いると極限 $h$ への収束は局所一様である。平均値性から $h$ は球内で調和であり、$h\le U$ かつ $h(x_0)=U(x_0)$ となる。別の候補 $v$ を加え $\max(h_j,v)$ を同じ球で持ち上げた極限を $k$ とすると $h\le k\le U$ であり、$x_0$ で $k-h=0$ になる。最大値原理で $k=h$ である。一方、$k\ge v$ なので全候補について $v\le h$ が従う。従って球上で $U=h$ である。
:::

:::example[球での一致]
単位円盤で $g=\cos(2\theta)$ のとき $h=x_1^2-x_2^2$ は候補である。差の比較で全候補 $v\le h$ だから、Perron の上限は $h$ に等しい。特に $(1/2,0)$ の値は $1/4$ である。
:::

:::example[点を除いた円盤]
$\Omega=B_1\setminus\{0\}\subset\R^2$ で、外周の値を $0$、原点で $1$ とする。$v\in\mathcal S_g$ と $\eps>0$ に対し $v\le \eps\log(1/|x|)$ が成り立つ。実際、内側の半径を十分小さくすれば $v$ の上界より右辺が大きく、外周では境界上極限が $0$ 以下であり、環状領域で比較できる。$\eps\downarrow0$ で $v\le0$、定数 $0$ も候補なので $U=0$ である。原点への境界値は実現されない。
:::

### Perron の上下からの比較

下の候補 $v$ が劣調和、上の候補 $w$ が優調和で、各境界点における $\limsup v$ と $\liminf w$ の間に所与の $g$ が入るとする。差 $v-w$ は球面平均以下であり、各境界点での上極限は非正である。$v-w>\eps$ の集合が空でなければその閉包は内部にコンパクトに入り、そこで正の最大値を取る。球面平均の等号により最大点集合は開でもあり、領域成分全体に広がる。境界での上極限と矛盾するので $v\le w$ である。この議論で下側の Perron 包絡線は上側の包絡線以下になる。両方が連続に境界値を実現すれば、差に最大値原理を使って両者が一致する。

:::example[Perron 包絡線を球の調和関数で挟む]
単位円板で $g=2+\cos\theta-2\sin2\theta$ とする。滑らかな境界値に対する明示解を Perron 法の候補として調べる。

円周データを定数と三角関数に分ける。角振動数 $j$ について半径の方程式は $R^{\prime\prime}+r^{-1}R^{\prime}-j^2r^{-2}R=0$ である。元の条件に合う根を選び、境界半径で係数を規格化する。この例での関数は次のとおりである。

$$
\begin{aligned}
u(r,\theta)&=- 2 r^{2} \sin{\left(2 \theta \right)} + r \cos{\left(\theta \right)} + 2\\
g&=- 2 \sin{\left(2 \theta \right)} + \cos{\left(\theta \right)} + 2,\qquad r=1\\
u_r&=\frac{\partial}{\partial r} \left(- 2 r^{2} \sin{\left(2 \theta \right)} + r \cos{\left(\theta \right)} + 2\right)=- 4 r \sin{\left(2 \theta \right)} + \cos{\left(\theta \right)}\\
u_\theta&=\frac{\partial}{\partial \theta} \left(- 2 r^{2} \sin{\left(2 \theta \right)} + r \cos{\left(\theta \right)} + 2\right)=- r \left(4 r \cos{\left(2 \theta \right)} + \sin{\left(\theta \right)}\right)
\end{aligned}
$$

二階微分を別々に計算する。半径方向の二階微分だけを角度方向の二階微分と足しても相殺しない。極座標の Laplacian の一階項も入れると各振動数で相殺する。

$$
\begin{aligned}
u_{rr}&=\frac{\partial^{2}}{\partial r^{2}} \left(- 2 r^{2} \sin{\left(2 \theta \right)} + r \cos{\left(\theta \right)} + 2\right)=- 4 \sin{\left(2 \theta \right)}\\
u_{\theta\theta}&=\frac{\partial^{2}}{\partial \theta^{2}} \left(- 2 r^{2} \sin{\left(2 \theta \right)} + r \cos{\left(\theta \right)} + 2\right)=r \left(16 r \sin{\left(\theta \right)} - 1\right) \cos{\left(\theta \right)}\\
\Delta u&=- 4 \sin{\left(2 \theta \right)} + \frac{- 4 r \sin{\left(2 \theta \right)} + \cos{\left(\theta \right)}}{r} + \frac{8 r \sin{\left(2 \theta \right)} - \cos{\left(\theta \right)}}{r}=0
\end{aligned}
$$

境界における値の平均と法線微分を読む。$j\ne0$ の角方向の平均は零であり、全体の平均は定数項である。またエネルギーは体積積分の方法と Green の境界積分の方法の両方で計算できる。

$$
\begin{aligned}
\text{円周の平均}&=\int\limits_{- \pi}^{\pi} \frac{- 2 \sin{\left(2 \theta \right)} + \cos{\left(\theta \right)} + 2}{2 \pi}\, d\theta=2\\
\partial_\nu u&=- 4 \sin{\left(2 \theta \right)} + \cos{\left(\theta \right)}\\
\int_{-\pi}^{\pi}|Du|^2\,\mathrm{d}\theta&=32 \pi r^{2} + 2 \pi=2 \pi \left(16 r^{2} + 1\right)\\
\int |Du|^2&=\int\limits_{0}^{1} r \left(32 \pi r^{2} + 2 \pi\right)\, dr=9 \pi\\
\int_{\partial\Omega} u\partial_\nu u&=\int\limits_{- \pi}^{\pi} \left(- 4 \sin{\left(2 \theta \right)} + \cos{\left(\theta \right)}\right) \left(- 2 \sin{\left(2 \theta \right)} + \cos{\left(\theta \right)} + 2\right)\, d\theta=9 \pi
\end{aligned}
$$

次に同じ振動数の半径因子を比較し、角度を固定した一点で値と半径方向の傾きを調べる。角度の取り方によって、振幅のあるモードでも点での値が零になることがある。

$$
\begin{aligned}
u&=\frac{\sqrt{2}}{4} + \frac{3}{2}=\frac{\sqrt{2} + 6}{4},\qquad r=\frac{1}{2},\theta=\frac{\pi}{4}\\
u_r&=-2 + \frac{\sqrt{2}}{2}=\frac{-4 + \sqrt{2}}{2},\qquad r=\frac{1}{2},\theta=\frac{\pi}{4}\\
u_\theta&=- \frac{\sqrt{2}}{4},\qquad r=\frac{1}{2},\theta=\frac{\pi}{4}
\end{aligned}
$$

この関数は境界で $g$ に一致し、調和なので下からの候補である。任意の別の候補との差は領域で劣調和、境界上極限が非正である。従って最大値原理でどの候補もこの関数以下になり、包絡線は式と一致する。
:::

## 障壁と正則境界点

:::definition[障壁と正則点]{#def-pde-barrier}
境界点 $\xi$ で $U(x)\to g(\xi)$ が全ての連続な $g$ について成り立つとき、$\xi$ を :term[正則境界点]{en="regular boundary point" yomi="せいそくきょうかいてん"} という。$\xi$ で $0$ に収束し、$\xi$ の各近傍の外では正の下限を持つ有界な正の優調和関数を :term[障壁]{en="barrier" yomi="しょうへき"} という。
:::

:::proposition[障壁による境界値の実現]{#prop-pde-barrier-regularity}
障壁 $b$ が存在すれば $\xi$ は正則である。逆に $\xi$ が正則なら障壁が存在する。
:::

:::proof
$g$ の連続性により $\xi$ の近くの境界では $g(\xi)-\eps\le g\le g(\xi)+\eps$ である。その外では $b$ の下限が正なので、$M$ を大きく取れば $g(\xi)-\eps-Mb$ は劣解の候補、$g(\xi)+\eps+Mb$ は全候補を上から比較する優調和関数になる。これらの境界上極限の比較から領域内で $U$ を挟み、$x\to\xi$、$\eps\downarrow0$ とする。逆を示すため $g(y)=|y-\xi|$ を取る。上限 $U$ は非負な調和関数で、正則性により $\xi$ で $0$ に収束する。各 $r>0$ で $\Omega\setminus B_r(\xi)$ に下限を与えるには $v(y)=A(|y-\xi|^2-r^2/2)$ を使う。$\Delta v=2nA\ge0$ であり、有界性を使って $A>0$ を十分小さくすると $v|_{\partial\Omega}\le |y-\xi|$ になる。よって $U\ge v\ge Ar^2/2$ がこの範囲で成り立つ。内部での正値性は強最大値原理で従う。
:::

:::example[外接球による障壁]
外部にある閉球 $\overline B_R(a)$ が $\overline\Omega$ と点 $\xi$ だけで接するとする。$n\ge3$ では $b=R^{2-n}-|x-a|^{2-n}$、$n=2$ では $b=\log(|x-a|/R)$ が障壁である。いずれも外側で調和であり、距離が $R$ より大きい点で正になる。例えば $\Omega=B_1(0)\subset\R^3$、$\xi=(1,0,0)$、$a=(2,0,0)$、$R=1$ では $b(0)=1/2$ である。
:::

:::example[環状領域の内側]
二次元の $1<r<2$ で、内周の値 $1$、外周の値 $0$ の解は $u=\log(2/r)/\log 2$ である。各内周点には内側の穴の内部から小さい円を一点で接触させられるので、内周も正則である。
:::

### 外接球から作る具体的な障壁

$n\ge3$ で $\overline\Omega$ と外側の球の共通点が $\xi$ のみとし、中心を $a$、半径を $R>0$ とする。$b(x)=R^{2-n}-|x-a|^{2-n}$ は領域で調和で、$\xi$ 以外の閉包において正である。二次元なら $b=\log(|x-a|/R)$ が同じ性質を持つ。$|g-g(\xi)|<\eps$ を満たす境界近傍を選んでから、それを除いた境界上で $\min b>0$ を使う。十分大きい $C$ に対して $g(\xi)-\eps-Cb$ と $g(\xi)+\eps+Cb$ は下と上の候補になり、解を挟む。$x\to\xi$、次に $\eps\to0$ とすれば正則性が得られる。

:::example[球が外に接する場合の近づき方]
三次元で外接球の中心を $a=(2,0,0)$、半径を $1$ とし、接点を $\xi=(1,0,0)$ とする。内部側の半径上 $x=(1-s,0,0)$ では $b=1-(1+s)^{-1}=s/(1+s)$ である。例えば $s=1/4$ では $b=1/5$ であり、接点へ近づくと $b/s\to1$ となる。障壁による挟み込みで一次の制御が現れるが、一般の連続境界データに同じ速度を要求するには $\eps$ に依存する係数 $C$ の制御も必要である。
:::

:::example[二次元の環から対数障壁を作る]
円環 $1<r<3$ で内側を $0$、外側を $1$ とする。円の外側に接したときに用いる障壁を規格化する。

放射方向だけに依存する解を考える。半径 $r$ の球面積が $\sigma_{n-1}r^{n-1}$ であることから、発散を計算すると $-(r^{n-1}u_r)^{\prime}=r^{n-1}f$ となる。これを一度積分し、中心を含む場合には中心に余分な点源がないよう流束の極限を零にする。環状領域では内側の境界条件で積分定数を決める。

$$
\begin{aligned}
u&=\frac{\log{\left(r \right)}}{\log{\left(3 \right)}}\\
u_r&=\frac{d}{d r} \frac{\log{\left(r \right)}}{\log{\left(3 \right)}}=\frac{1}{r \log{\left(3 \right)}}\\
r^{n-1}u_r&=\frac{1}{\log{\left(3 \right)}}\\
-r^{n-1}f&=\frac{d}{d r} \frac{1}{\log{\left(3 \right)}}=0\\
-\Delta u&=0
\end{aligned}
$$

次に値と流束を境界で読む。外側では法線は半径が増える向きである。内側にも境界がある場合は法線が反対向きになる。源項の体積積分を角度方向に積分すると面積係数が現れ、流束との差を確かめられる。

$$
\begin{aligned}
u&=0,\qquad r=1\\
u&=1,\qquad r=3\\
u_r&=\frac{1}{3 \log{\left(3 \right)}},\qquad r=3\\
\sigma_{n-1}&=2 \pi\\
\text{全源量}&=\int\limits_{1}^{3} 0\, dr=0\\
\text{負の全流束}&=0
\end{aligned}
$$

観測半径で値と半径微分を計算する。勾配の二乗を用いてエネルギーも計算できる。観測点で勾配が小さい場合でも、球面積の因子があるため、そこを通る全流束が小さいとは限らない。

$$
\begin{aligned}
u&=\frac{\log{\left(\sqrt{3} \right)}}{\log{\left(3 \right)}}=\frac{1}{2},\qquad r=\sqrt{3}\\
u_r&=\frac{\sqrt{3}}{3 \log{\left(3 \right)}},\qquad r=\sqrt{3}\\
\sigma_{n-1}r^{n-1}u_r&=\frac{2 \pi}{\log{\left(3 \right)}},\qquad r=\sqrt{3}\\
\frac12\int |Du|^2&=\int\limits_{1}^{3} \frac{\pi}{r \log{\left(3 \right)}^{2}}\, dr=\frac{\pi}{\log{\left(3 \right)}}
\end{aligned}
$$

半径 $\sqrt3$ で値は $1/2$ である。$\partial_r u(1)=1/\log3$ なので半径 $1$ に近いとき $u$ は $r-1$ に比例する。三次元以上の $r^{2-n}$ の代わりに対数が現れるのは、$r^{1-n}$ を一度積分した結果である。
:::

### 局所障壁を領域全体へ延ばす

$B_{2r}(\xi)\cap\Omega$ で非負な連続優調和関数 $b$ があり、$b(x)\to0$（$x\to\xi$）とする。さらに $\overline\Omega\cap\partial B_r(\xi)$ 上で連続に正とし、その最小値を $m>0$ とする。さらにどの $0<\rho<r$ に対しても $\Omega\cap(B_r(\xi)\setminus B_\rho(\xi))$ 上で $b$ の下限は正とする。$0<c<m$ を選び、$\Omega\cap B_r(\xi)$ で $\min(b,c)$、外側で定数 $c$ と置く。内側での最小は優調和である。これは負を取った二つの劣調和関数の最大が劣調和であることから従う。接合面では $b>c$ なので、その近傍でも関数は定数 $c$ となり、優調和性が保たれる。従って他の境界点から離れた場所で正の下界を持つ局所障壁があれば、この切り詰めで大域的な比較に使える。

:::example[外側の円に接する局所領域]
二次元で $a=(0,-1)$、接点 $\xi=(0,0)$ とし、領域が $y>0$ にあるとする。$b(x,y)=\log\sqrt{x^2+(y+1)^2}$ は調和で、閉じた上半平面では原点を除けば正である。原点から距離 $r=1/2$ の上半円では $x^2+(y+1)^2=1+1/4+2y\ge5/4$ なので、$b\ge\tfrac12\log(5/4)>0$ である。例えばその半分を $c$ に取って延長できる。境界に沿って接点に近づくと $b(x,0)=\tfrac12\log(1+x^2)$ で二次、垂直なら $b(0,y)=\log(1+y)$ で一次であり、障壁の増え方は方向によって異なる。
:::

### 外部円錐による障壁

境界点 $\xi$ から出る開いた円錐が局所的に領域の外にあるとする。$\xi$ を原点、円錐の軸との角度を $\theta$ とすれば、領域側では $\theta\ge\theta_0>0$ である。$n\ge2$ とし、
$$
\varphi(\theta)=1+\int_{\theta_0}^{\theta}
\frac{\int_s^\pi\sin^{n-2}q\,\dd q}{\sin^{n-2}s}\,\dd s,
\qquad \theta_0\le\theta\le\pi
$$
と置く。$\varphi$ は正で有界である。$\theta\to\pi$ では積分核は $\pi-\theta$ の滑らかな奇関数として延びるので、$\varphi$ は軸の反対側でも滑らかな球面関数である。球面 Laplacian の軸対称部分を使うと $(\sin^{n-2}\theta\,\varphi')'=-\sin^{n-2}\theta$、すなわち $\Delta_{\mathbb S^{n-1}}\varphi=-1$ である。

そこで $b=r^\beta\varphi(\theta)$ とすると
$$
\Delta b=r^{\beta-2}\big(\beta(\beta+n-2)\varphi-1\big)\le0
$$
となるよう十分小さい $\beta>0$ を選べる。$b$ は原点で零へ近づき、どの $\rho>0$ に対しても $r\ge\rho$ では $\rho^\beta$ 以上である。従って局所障壁を切り詰めて大域障壁を得られる。例えば平面の開き角 $3\pi/2$ の内角では外接円を頂点に置けないが、補集合に円錐を取れるので頂点は正則である。

## 一重層ポテンシャルと二重層ポテンシャル

:::definition[層ポテンシャル]{#def-pde-layer-potentials}
滑らかな有界領域 $\Omega$ の境界 $\Sigma$ と外向き法線 $\nu$ に対し

$$
S\mu(x)=\int_\Sigma\Gamma(x-y)\mu(y)\,\dd S_y,\qquad
D\mu(x)=\int_\Sigma\partial_{\nu_y}\Gamma(x-y)\mu(y)\,\dd S_y
$$

を :term[一重層ポテンシャル]{en="single-layer potential" yomi="ひとえそうぽてんしゃる"} と :term[二重層ポテンシャル]{en="double-layer potential" yomi="ふたえそうぽてんしゃる"} という。いずれも境界から離れた点で調和である。
:::

:::proposition[層ポテンシャルの境界極限]{#prop-pde-layer-jumps}
$\Sigma$ は $C^2$、$\mu$ は Hölder 連続とする。境界上の核を直接積分した作用素を $K\mu=D\mu|_{\rm kernel}$、$K'\mu=(\partial_{\nu_x}S\mu)|_{\rm kernel}$ と書く。内側と外側からの法線方向極限をそれぞれ $\mathrm{in},\mathrm{out}$ とすれば

$$
\begin{aligned}
(S\mu)_{\rm in}&=(S\mu)_{\rm out},\\
(D\mu)_{\rm in}&=K\mu-\tfrac12\mu,\qquad (D\mu)_{\rm out}=K\mu+\tfrac12\mu,\\
(\partial_\nu S\mu)_{\rm in}&=K'\mu+\tfrac12\mu,\qquad
(\partial_\nu S\mu)_{\rm out}=K'\mu-\tfrac12\mu .
\end{aligned}
\label{eq-pde-layer-jumps}
$$

:::

:::proof
境界点を原点に取り、外法線を $e_n$ とする。境界は $y=(z,\psi(z))$、$\psi(z)=O(|z|^2)$、$D\psi(z)=O(|z|)$ と書け、評価点は $x=h e_n$ とする。一重層の核は小球内で $C|z|^{2-n}$（二次元では $C(1+|\log |z||)$）以下であり、境界次元は $n-1$ なので優収束で両側の値が一致する。微分核では、表面要素を含めた $D$ の分子は $h-\psi+D\psi\cdot z=h+O(|z|^2)$、$\partial_{\nu_x}S$ の分子は $-(h-\psi)(1+O(|z|^2))$ になる。分母は $\sigma_{n-1}(|z|^2+(h-\psi)^2)^{n/2}$ である。$h$ を分子に持つ平面核 $\pm h/(\sigma_{n-1}(|z|^2+h^2)^{n/2})$ を引く。差の $O(|z|^2)$ の項は $C|z|^{2-n}$ で抑えられる。分母の差に由来する項も、二乗距離の差が $O(|h||z|^2+|z|^4)$ で、両方の距離が $\sqrt{|z|^2+h^2}$ と同程度であることから $C|z|^{2-n}$ 以下である。よって差の極限は境界核の積分である。平面核は符号を除き全平面で積分 $1/2$ を持ち、原点の外で積分が $0$ に収束する。積分値は球の上半面への中心射影の面積要素 $h(|z|^2+h^2)^{-n/2}\dd z$ から得られる。従って残りは $D$ では $\mathrm{sgn}(h)\mu(0)/2$、$\partial_\nu S$ ではその負である。$h<0$ が内側なので表示の符号が従う。
:::

:::example[一定密度の二重層]
$\mu=1$ とすると $D1=-1$ が内側、$D1=0$ が外側である。これは $y$ について発散定理を適用し、$x\in\Omega$ の場合だけ $-\Delta\Gamma=\delta_x$ の寄与があるからである。式 :eqref[eq-pde-layer-jumps] によれば $K1=-1/2$ である。
:::

:::example[球面の一定一重層]
半径 $R$ の三次元球で $\mu=1$ とする。内部は球対称で有界なので定数、外部は球対称かつ無限遠で $0$ になるので $C/r$ である。全量が $4\pi R^2$ なので $C=R^2$ である。連続性から内部は $R$ となる。法線微分の内側は $0$、外側は $-1$ で、差が $1$ である。
:::

### 二重層の法線微分は両側で一致する

滑らかな境界と密度を取り、密度 $\mu$ を $E\in C^2(\overline\Omega)$ に延長する。Green の恒等式を領域内部・外部の点について書けば
$$
D\mu=N((-\Delta E)\mathbf1_\Omega)+S(\partial_\nu E)-\mathbf1_\Omega E
$$
となる。符号は $-\Delta\Gamma=\delta$ と外向き法線による。体積ポテンシャルの一階微分は連続である。実際、源は有界で $D\Gamma$ の小球での積分は $O(r)$ である。従って特異部分を一様に捨てて連続な積分に近似できる。一重層の法線微分の内部跡と外部跡の差は $\partial_\nu E$ であり、最後の項の法線微分の差は $-\partial_\nu E$ なので相殺する。従って二重層は値に跳躍を持ち、法線微分は連続する。ここで法線は両側で同じ $\nu$ を用いる。

:::example[円周の一重層を二つの調和拡張でつなぐ]
単位円周で跡 $g=\cos(3\theta)$ を持つ体積源のないポテンシャルを作る。内側と外側の法線の違いを調べる。

円周データを定数と三角関数に分ける。角振動数 $j$ について半径の方程式は $R^{\prime\prime}+r^{-1}R^{\prime}-j^2r^{-2}R=0$ である。元の条件に合う根を選び、境界半径で係数を規格化する。この例での関数は次のとおりである。

$$
\begin{aligned}
u(r,\theta)&=\frac{\cos{\left(3 \theta \right)}}{r^{3}}\\
g&=\cos{\left(3 \theta \right)},\qquad r=1\\
u_r&=\frac{\partial}{\partial r} \frac{\cos{\left(3 \theta \right)}}{r^{3}}=- \frac{3 \cos{\left(3 \theta \right)}}{r^{4}}\\
u_\theta&=\frac{\partial}{\partial \theta} \frac{\cos{\left(3 \theta \right)}}{r^{3}}=- \frac{3 \sin{\left(3 \theta \right)}}{r^{3}}
\end{aligned}
$$

二階微分を別々に計算する。半径方向の二階微分だけを角度方向の二階微分と足しても相殺しない。極座標の Laplacian の一階項も入れると各振動数で相殺する。

$$
\begin{aligned}
u_{rr}&=\frac{\partial^{2}}{\partial r^{2}} \frac{\cos{\left(3 \theta \right)}}{r^{3}}=\frac{12 \cos{\left(3 \theta \right)}}{r^{5}}\\
u_{\theta\theta}&=\frac{\partial^{2}}{\partial \theta^{2}} \frac{\cos{\left(3 \theta \right)}}{r^{3}}=- \frac{9 \cos{\left(3 \theta \right)}}{r^{3}}\\
\Delta u&=0
\end{aligned}
$$

境界における値の平均と法線微分を読む。$j\ne0$ の角方向の平均は零であり、全体の平均は定数項である。またエネルギーは体積積分の方法と Green の境界積分の方法の両方で計算できる。

$$
\begin{aligned}
\text{円周の平均}&=\int\limits_{- \pi}^{\pi} \frac{\cos{\left(3 \theta \right)}}{2 \pi}\, d\theta=0\\
\partial_\nu u&=3 \cos{\left(3 \theta \right)}\\
\int_{-\pi}^{\pi}|Du|^2\,\mathrm{d}\theta&=\frac{18 \pi}{r^{8}}\\
\int |Du|^2&=\int\limits_{1}^{\infty} \frac{18 \pi}{r^{7}}\, dr=3 \pi\\
\int_{\partial\Omega} u\partial_\nu u&=\int\limits_{- \pi}^{\pi} 3 \cos^{2}{\left(3 \theta \right)}\, d\theta=3 \pi
\end{aligned}
$$

次に同じ振動数の半径因子を比較し、角度を固定した一点で値と半径方向の傾きを調べる。角度の取り方によって、振幅のあるモードでも点での値が零になることがある。

$$
\begin{aligned}
u&=\frac{1}{8},\qquad r=2,\theta=0\\
u_r&=- \frac{3}{16},\qquad r=2,\theta=0\\
u_\theta&=0,\qquad r=2,\theta=0
\end{aligned}
$$

外側の式は $r^{-3}\cos(3\theta)$、内側の式は $r^3\cos(3\theta)$ である。円周上で値は一致し、円板の外向き法線に固定すると内外の法線微分の差は $6\cos(3\theta)$ になる。従って密度 $6\cos(3\theta)$ の一重層に一致する。実際両者の差は全平面で分布調和で、平均零の密度により無限遠で零に行くので Weyl と Liouville から零になる。
:::

### 二次元で定数モードに現れる対数

半径 $R$ の円周上に密度 $\mu=m$ を一様に置く。二次元の基本解は $-(2\pi)^{-1}\log|x-y|$ であり、$r=|x|$ だけに依存する。$r<R$ では $\log|R e^{i\theta}-x|=\log R+\Re\log(1-(x/R)e^{-i\theta})$ と分ける。後者は $|x/R|<1$ の収束べき級数であり、各非零振動数の積分が零になる。従って内部値は $-mR\log R$ である。$r>R$ では $r$ をくくり出して同じ議論を使い、外部値 $-mR\log r$ を得る。境界で値は一致し、同じ半径方向で測った微分の内外差は $0-(-mR/R)=m$ である。

:::example[単位円の境界作用素の核]
$R=1$ では密度 $m=2$ を置いても一重層の境界値は零である。実際内部も恒等的に零だが、外部では $-2\log r$ になる。この境界作用素の零モードは「外部で無限遠に零へ行く」という仮定を満たさない。平均零密度の非零角モードなら外部は減衰し、単位円の境界作用素の固有値は $1/(2j)$（振動数 $j\ge1$）である。従って同じ作用素で定数モードを逆に解こうとすると、他の角モードとは異なる適合条件や定数の補正が必要になる。例えば $R=2$ へ拡大しただけで、密度 $2$ の境界値は $-4\log2$ になり、定数モードの零は消える。これは二次元で基本解に加法定数の自由があることにも対応する。
:::

## 積分方程式による境界値問題の解法

調和関数を $u=D\mu$ とおくと Dirichlet 条件は $(K-\frac12 I)\mu=g$ に帰着する。右辺 $f$ がある場合は、先に $Nf$ を引き、境界データを $g-Nf|_\Sigma$ に変える。調和関数の一意性だけから、表現する密度の単射性は従わない。そこで $K$ と $K'$ の随伴関係を使い、領域の外側の一重層解も調べる。

:::lemma[$L^2$ の境界作用素]{#lem-pde-layer-l2}
$C^2$ 境界上で $K,K'$ は $L^2(\Sigma)$ 上のコンパクト作用素で、互いに随伴である。:eqref[eq-pde-layer-jumps] の極限は $\mu\in L^2(\Sigma)$ に対しても、平行面の値を法線に沿って境界へ引き戻した $L^2(\Sigma)$ で成り立つ。
:::

:::proof
境界上の $K,K'$ の核は $C|x-y|^{2-n}$ 以下である。これは境界の二点間の差と各法線の内積が $O(|x-y|^2)$ であるためである。$|x-y|<\eta$ の核の絶対値の積分は $x,y$ のどちらを固定しても一様に $O(\eta)$ である。核 $k$ が両変数について絶対積分の上界 $M_1,M_2$ を持つとき、重み $|k(x,y)|$ 付きの Cauchy–Schwarz で $\norm{K\mu}_2^2\le M_1M_2\norm\mu_2^2$ である。距離 $\eta$ 以下で零、$2\eta$ 以上で $1$ の連続な切断を核に掛ける。捨てた部分の作用素ノルムは $O(\eta)$ である。残る連続核は、有限個の小領域に分けて各領域で定数近似することで有限階作用素の極限になる。よってコンパクト性が従う。基本解の対称性から $K$ の核の変数を入れ替えると $K'$ の核となるので、Fubini と近似で随伴性が従う。

同じ局所座標で、平行面を境界に引き戻した微分核は $C(|h|/(|z|^2+h^2)^{n/2}+|z|^{2-n})$ 以下である。各変数についての絶対積分は、面素が一様比較できる局所座標で見ると $h$ に依存せず有界である。平滑な密度では :ref[prop-pde-layer-jumps] の証明の誤差評価と連続性の一様性により一様収束する。$L^2$ の密度を平滑なものに近似し、一様な作用素ノルム評価を使うと同じ極限が $L^2$ で従う。一重層の連続極限でも、局所可積分な核の同様の評価を使える。
:::

:::theorem[外部が連結な領域の Dirichlet 積分方程式]{#thm-pde-layer-dirichlet}
$n\ge2$、$\Omega$ は有界で $C^2$ 境界を持ち、$\R^n\setminus\overline\Omega$ が連結とする。このとき $K-\frac12 I$ は $L^2(\Sigma)$ 上で可逆である。
:::

:::proof
[関数解析の Fredholm の択一定理](/math/functional-analysis/compact-operators/) により随伴作用素 $K'-\frac12 I$ の核が零であることを示せばよい。その核に属する $\mu$ を取り、$v=S\mu$ とする。外側からの法線微分は $0$ である。$n\ge3$ では、$\mu\in L^1(\Sigma)$ と基本解の評価から無限遠で $v=O(|x|^{2-n})$、$Dv=O(|x|^{1-n})$ である。

$n=2$ では $M=\int_\Sigma\mu$ とおく。$|x|\to\infty$ で基本解を $y$ について展開すると $v=-(M/2\pi)\log|x|+O(|x|^{-1})$、$\partial_r v=-M/(2\pi|x|)+O(|x|^{-2})$ である。一方、外側平行面と外円の間で発散定理を使い、平行面を境界に近づけると外円の総流束は零になる。これを無限遠へ送ると $M=0$ である。従って $v=O(|x|^{-1})$、$Dv=O(|x|^{-2})$ が得られる。

境界から距離 $h>0$ の外側平行面と半径 $R$ の外球で囲んだ領域で $\int |Dv|^2=\int v\partial_\nu v$ を用いる。平行面の法線は境界の法線と逆向きであり、$h\downarrow0$ で積分は $L^2$ の値と法線微分の積に収束して $0$ になる。無限遠で得た減衰により外球の寄与も消える。従って外部全体で $Dv=0$ であり、連結性と減衰により $v=0$ である。一重層の両側値は等しいため、内部でも境界値が $0$ である。内部の平行面で同じエネルギー式を使うと、値は $L^2$ で $0$ に収束し、法線微分は一様 $L^2$ 有界なので内部でも $Dv=0$ である。従って法線微分の跳躍は $0$ であり、:eqref[eq-pde-layer-jumps] により $\mu=0$ である。
:::

:::example[円の二重層方程式]
単位円で $x,y$ が境界上にあるとき $\nu_y=y$、$(x-y)\cdot y=\cos(\theta_x-\theta_y)-1=-|x-y|^2/2$ である。従って $K\mu=-(4\pi)^{-1}\int_0^{2\pi}\mu(\theta)\,\dd\theta$ は平均の $-1/2$ 倍である。$\mu=m+\mu_0$ と定数部分と平均 $0$ の部分に分解すると $(K-I/2)\mu=-m-\mu_0/2$ である。例えば $g=3+2\cos\theta$ なら $\mu=-3-4\cos\theta$ であり、得られる内部解は $u=3+2x_1$ である。
:::

:::example[穴があると密度に核が生じる]
環状領域で内側の境界にのみ密度 $1$ を置く。環状領域の外法線は穴の外法線の反対である。穴には評価点が入らないため発散定理で、この二重層は環状領域内で $0$ になる。したがって同次 Dirichlet 問題の解が $0$ しかなくても、境界積分方程式は非自明な核を持ち得る。コンパクト作用素の択一性により、可解性を調べるには随伴の核を求めればよい。外部の成分ごとの条件は、次の補正の構成で調べる。
:::

### 外部に複数の成分がある場合の補正

$\Omega$ は有界連結な $C^2$ 領域、境界データは $L^2$ とし、$n\ge2$ とする。外部の有界成分を $O_1,\ldots,O_m$、各内部の一点を $a_j$ とする。内部 Dirichlet の作用素を $A=K-I/2$ とすると、その随伴の核に属する $\mu$ は $S\mu$ の外側法線微分が零である。この関数は各穴で定数 $c_j$、無限遠を含む外部で零になる。各成分で Green のエネルギー恒等式を使えば分かる。$n=2$ では非有界成分の総流束が零であることから、:ref[thm-pde-layer-dirichlet] の証明と同じく対数項が消える。従っていずれの次元でも無限遠の境界積分は零に収束する。跡の条件が $L^2$ の場合は、平行面で恒等式を使って極限を取る。

境界の関数 $\Gamma(\cdot-a_j)$ と $\mu$ の内積は $S\mu(a_j)=c_j$ である。従って $\mu$ が全ての補正関数に直交すれば、$S\mu$ の外側の値は零である。一重層の連続性から内側の境界値も零となり、エネルギー恒等式で内部でも零になる。法線微分の跳躍から $\mu=0$ である。Fredholm 性により $\operatorname{Ran} A$ は閉じて余次元有限なので、これら $m$ 個の補正関数を足せば境界空間全体を張る。実際、和に直交する元はいま示したように零であり、有限次元を加えた和は閉じている。

従って $D\sigma+\sum_j c_j'\Gamma(\cdot-a_j)$ という表現なら全ての境界データを扱える。一方、各穴の境界だけに $1$ を置いた密度は内部の二重層を零にし、これらは独立である。指数零と上の単射から核の次元はちょうど $m$ なので、補正関数の剰余類は像を除いた次元 $m$ の商空間を張るため独立であり、係数 $c_j^{\prime}$ は一意になる。その上で各穴の密度の積分を零に正規化すれば密度も一意にできる。

:::example[積分方程式の零平均モードを逆算する]
単位円周で一重層の境界値に $2\cos2\theta+\sin4\theta$ を指定する。先に外側の調和拡張から密度を求める。

円周データを定数と三角関数に分ける。角振動数 $j$ について半径の方程式は $R^{\prime\prime}+r^{-1}R^{\prime}-j^2r^{-2}R=0$ である。元の条件に合う根を選び、境界半径で係数を規格化する。この例での関数は次のとおりである。

$$
\begin{aligned}
u(r,\theta)&=\frac{2 \cos{\left(2 \theta \right)}}{r^{2}} + \frac{\sin{\left(4 \theta \right)}}{r^{4}}\\
g&=\sin{\left(4 \theta \right)} + 2 \cos{\left(2 \theta \right)},\qquad r=1\\
u_r&=\frac{\partial}{\partial r} \left(\frac{2 \cos{\left(2 \theta \right)}}{r^{2}} + \frac{\sin{\left(4 \theta \right)}}{r^{4}}\right)=- \frac{4 \left(r^{2} \cos{\left(2 \theta \right)} + \sin{\left(4 \theta \right)}\right)}{r^{5}}\\
u_\theta&=\frac{\partial}{\partial \theta} \left(\frac{2 \cos{\left(2 \theta \right)}}{r^{2}} + \frac{\sin{\left(4 \theta \right)}}{r^{4}}\right)=- \frac{4 \left(r^{2} \sin{\left(2 \theta \right)} - \cos{\left(4 \theta \right)}\right)}{r^{4}}
\end{aligned}
$$

二階微分を別々に計算する。半径方向の二階微分だけを角度方向の二階微分と足しても相殺しない。極座標の Laplacian の一階項も入れると各振動数で相殺する。

$$
\begin{aligned}
u_{rr}&=\frac{\partial^{2}}{\partial r^{2}} \left(\frac{2 \cos{\left(2 \theta \right)}}{r^{2}} + \frac{\sin{\left(4 \theta \right)}}{r^{4}}\right)=\frac{4 \left(3 r^{2} \cos{\left(2 \theta \right)} + 5 \sin{\left(4 \theta \right)}\right)}{r^{6}}\\
u_{\theta\theta}&=\frac{\partial^{2}}{\partial \theta^{2}} \left(\frac{2 \cos{\left(2 \theta \right)}}{r^{2}} + \frac{\sin{\left(4 \theta \right)}}{r^{4}}\right)=- \frac{8 \left(r^{2} \cos{\left(2 \theta \right)} + 2 \sin{\left(4 \theta \right)}\right)}{r^{4}}\\
\Delta u&=\frac{- \frac{4 \cos{\left(2 \theta \right)}}{r^{3}} - \frac{4 \sin{\left(4 \theta \right)}}{r^{5}}}{r} - \frac{8 \left(\cos{\left(2 \theta \right)} + \frac{2 \sin{\left(4 \theta \right)}}{r^{2}}\right)}{r^{4}} + \frac{4 \left(3 \cos{\left(2 \theta \right)} + \frac{5 \sin{\left(4 \theta \right)}}{r^{2}}\right)}{r^{4}}=0
\end{aligned}
$$

境界における値の平均と法線微分を読む。$j\ne0$ の角方向の平均は零であり、全体の平均は定数項である。またエネルギーは体積積分の方法と Green の境界積分の方法の両方で計算できる。

$$
\begin{aligned}
\text{円周の平均}&=\int\limits_{- \pi}^{\pi} \frac{\sin{\left(4 \theta \right)} + 2 \cos{\left(2 \theta \right)}}{2 \pi}\, d\theta=0\\
\partial_\nu u&=4 \sin{\left(4 \theta \right)} + 4 \cos{\left(2 \theta \right)}=4 \left(\sin{\left(4 \theta \right)} + \cos{\left(2 \theta \right)}\right)\\
\int_{-\pi}^{\pi}|Du|^2\,\mathrm{d}\theta&=\frac{32 \pi}{r^{6}} + \frac{32 \pi}{r^{10}}=\frac{32 \pi \left(r^{4} + 1\right)}{r^{10}}\\
\int |Du|^2&=\int\limits_{1}^{\infty} r \left(\frac{32 \pi}{r^{6}} + \frac{32 \pi}{r^{10}}\right)\, dr=12 \pi\\
\int_{\partial\Omega} u\partial_\nu u&=\int\limits_{- \pi}^{\pi} \left(\sin{\left(4 \theta \right)} + 2 \cos{\left(2 \theta \right)}\right) \left(4 \sin{\left(4 \theta \right)} + 4 \cos{\left(2 \theta \right)}\right)\, d\theta=12 \pi
\end{aligned}
$$

次に同じ振動数の半径因子を比較し、角度を固定した一点で値と半径方向の傾きを調べる。角度の取り方によって、振幅のあるモードでも点での値が零になることがある。

$$
\begin{aligned}
u&=\frac{1}{16} + \frac{\sqrt{2}}{4}=\frac{1 + 4 \sqrt{2}}{16},\qquad r=2,\theta=\frac{\pi}{8}\\
u_r&=- \frac{\sqrt{2}}{4} - \frac{1}{8}=- \frac{1 + 2 \sqrt{2}}{8},\qquad r=2,\theta=\frac{\pi}{8}\\
u_\theta&=- \frac{\sqrt{2}}{2},\qquad r=2,\theta=\frac{\pi}{8}
\end{aligned}
$$

零平均の角モードでは $j$ 次の境界値 $g_j$ に対し密度は $2jg_j$ である。理由は内部の法線微分が $jg_j$、外部では同じ外法線に沿って $-jg_j$ になるためである。この場合 $\mu=8\cos2\theta+8\sin4\theta$。$j=0$ はこの式では決められず、二次元での全電荷と無限遠の対数項を別に扱う。
:::

### 環状領域で補正項を一つ使う

三次元の環 $1<|x|<2$ を考え、内側で $g=1$、外側で $g=0$ とする。調和解が放射状であることは、回転した解も同じ境界値を持ち、一意性があることから従う。求める解は $u=A+B/r$ であり、連立条件 $A+B=1,A+B/2=0$ を解くと $A=-1$、$B=2$ である。これを二重層と補正の式に入れる。外側の球面に密度 $1$ を置けば、環のすべての点はその球の内側にあるので二重層の値は $-1$ である。内側の球面の密度を零にし、穴の中心を $a=0$ に取れば、補正 $8\pi\Gamma(x-a)=2/r$ が残る。従ってこの一つの補正で両方の境界値を実現できる。

内側の球面だけに密度 $1$ を置いても、環の点はその球の外側にあるため二重層は環で恒等的に零である。従って境界値が同じでも密度は自由に変えられる。内側の球面で密度の平均を零にするという正規化により、上の選び方が一意になる。一方、補正を使わず純粋な二重層でこの $g$ を表すことはできない。実際、回転平均を取っても境界値は変わらず、回転平均した密度は二つの球面上でそれぞれ定数になる。しかし定数密度から得られる環内の値は外側の球面の定数から作る一つの定数であり、内外の異なる値を作れない。この具体的な可解性の不足が、外部の穴ごとに補正を足した理由である。

:::example[同じ境界値の解が持つ流束]
上の $u=-1+2/r$ は $u_r=-2/r^2$ である。外側の面で環の外向き微分は $-1/2$、内側では外向きが半径と逆なので $2$ である。面積を掛けるとそれぞれ $-8\pi,8\pi$ となり、領域内に源がないので和は零である。勾配ノルムの二乗は Green の公式で $\int u\partial_\nu u=8\pi$。体積積分からも $4\pi\int_1^2 (4/r^4)r^2\dd r=8\pi$ を得る。補正を境界値の代数だけで選んでも、この流束とエネルギーを同時に満たす。
:::

## Neumann 問題

:::proposition[適合条件と一意性]{#prop-pde-neumann-compatibility}
$\Omega$ は滑らかな有界連結領域とする。$-\Delta u=f$、$\partial_\nu u=h$ の古典解が存在するためには

$$
\int_\Omega f\,\dd x+\int_{\partial\Omega}h\,\dd S=0
\label{eq-pde-neumann-compatibility}
$$

が必要である。解は加法定数を除いて一意である。
:::

:::proof
方程式を積分して発散定理を使うと条件が必要である。二解の差 $w$ に Green の公式を用いると $\int|Dw|^2=\int_\Sigma w\partial_\nu w-\int_\Omega w\Delta w=0$ であるから、連結性より $w$ は定数である。
:::

弱解の存在も証明できる。平均 $0$ の $H^1(\Omega)$ 上で $a(u,v)=\int Du\cdot Dv$ とし、右辺を $\ell(v)=\int_\Omega fv+\int_\Sigma h v$ とする。$f,h$ がそれぞれ領域と境界で $L^2$ に属すれば、トレース評価で $\ell$ は有界である。平均 $0$ の Poincaré 不等式で $a$ は強圧的であるから、[Hilbert 空間上の有界作用素](/math/functional-analysis/hilbert-operators/) の Lax–Milgram の定理で唯一の平均 $0$ の弱解が存在する。:eqref[eq-pde-neumann-compatibility] により $\ell(1)=0$ であるため、任意の試験関数に定数を足しても弱形式は同じである。一重層で調和部分を表す場合は $(K'+\frac12 I)\mu=h-\partial_\nu Nf$ を解く。左辺の積分は常に $0$ であり、右辺の積分が $0$ になる条件は :eqref[eq-pde-neumann-compatibility] である。

:::proposition[一重層による Neumann 解の構成]{#prop-pde-neumann-layer}
$n\ge3$、$\Omega$ は有界連結な $C^2$ 領域とする。$g\in L^2(\Sigma)$、$\int_\Sigma g=0$ なら $(K'+\frac12 I)\mu=g$ が解ける。
:::

:::proof
$B=K'+\frac12 I$ とする。:ref[lem-pde-layer-l2] により Fredholm 作用素で指数は $0$ である。任意の $\mu$ について $v=S\mu$ は調和であり、内部の平行面で発散定理を使い極限を取ると $\int_\Sigma B\mu=0$ である。従って余次元は少なくとも $1$ である。一方、$\mu\in\ker B$ なら、内部法線微分が零である。平行面を使ったエネルギー評価で内部の $v$ は定数 $c_\mu$ となる。$c_\mu=0$ なら外部の境界値も零である。外部で同じエネルギー計算を行い、無限遠の項を $O(R^{2-n})$ で消すと全ての外部成分で勾配が零になる。よって法線微分の跳躍から $\mu=0$ である。従って $\mu\mapsto c_\mu$ は核から $\R$ への単射であり、核の次元は高々 $1$ である。指数零から余次元はちょうど $1$ であり、像は積分零の閉部分空間に等しい。
:::

二次元では $B=K'+\frac12 I$ の核で、内部の定数に代えて総密度 $M=\int_\Sigma\mu$ を使う。$B\mu=0$ なら $v=S\mu$ は内部で定数 $c$ である。もし $M=0$ なら、外部の非有界成分では $v=O(r^{-1})$、$Dv=O(r^{-2})$ となる。一重層の連続性から全外部成分の有限境界で値は同じ $c$ である。外部を外円で切り、平行面から境界への極限を取ると、外部のエネルギーは $-c\int_\Sigma(\partial_\nu v)_{\rm out}$ に等しい。外円の寄与は減衰する。跳躍式と $B\mu=0$ によりこの積分は $-M=0$ だから、外部では勾配が零になる。非有界成分の減衰から $c=0$ となり、跳躍式により $\mu=0$ である。したがって $M$ は核上で単射である。境界積分が零という像の制約と Fredholm 指数零を合わせると、核と像の余次元はともに一である。従って任意の平均零データを実現でき、密度も平均零と指定すれば一意に選べる。調和解の定数は必要に応じて一重層の値に加える。

:::example[定数の源を法線流束で出す]
三次元単位球で $f=6$、$h=-2$ とする。$6(4\pi/3)-2(4\pi)=0$ なので適合する。$u=-r^2+C$ が解であり、球上の平均が $3/5$ の $r^2$ を引いているので、平均 $0$ の解は $C=3/5$ と選ぶ。
:::

:::example[円の Neumann データ]
単位円で $h=2\cos(2\theta)$ は平均 $0$ である。$u=r^2\cos(2\theta)$ の半径方向微分を $r=1$ で取ると $h$ になり、領域内の平均も $0$ なので、この $u$ は規格化した解である。
:::

### 弱 Neumann 問題の空間

連結な有界 Lipschitz 領域で $-\Delta u=f,\ \partial_\nu u=g$ を考える。$f\in L^2(\Omega)$、$g\in L^2(\partial\Omega)$ とする。試験関数 $1$ による適合条件は $\int_\Omega f+\int_{\partial\Omega}g=0$ である。この条件の下で、平均零の $H^1$ 関数の空間に
$$
\int Du\cdot Dv=\int fv+\int_{\partial\Omega}gv
$$
を課す。トレース評価により右辺は有界である。平均零の空間で Poincaré の不等式を確かめるには、勾配が零に近づくノルム $1$ の列があると仮定し、Rellich のコンパクト性で $L^2$ 極限を取ればよい。極限の勾配は零なので連結性から定数であり、平均零とノルム $1$ に矛盾する。この不等式により勾配の内積が完全なので Riesz の表現で唯一の平均零解を得る。適合条件が定数を消すため試験関数は $H^1$ 全体に拡張できる。

:::example[一次元で法線の符号が異なる二端点]
$(0,2)$ で $-u''=2$、$g(0)=1$、$g(2)=-5$ とする。領域積分は $4$、端点の和は $-4$ なので適合する。左端の外法線は $-1$ なので $u'(0)=-1$、従って $u=-x^2-x+C$ となる。右端で $u'(2)=-5$ を満たす。平均は $-4/3-1+C$ なので $C=7/3$ が平均零の解である。両方の端点で $u'$ をそのまま指定した場合は、適合条件が別になる。
:::

:::example[Neumann 解の定数を平均零で決める]
区間 $(0,3)$ で $-u_{xx}=2x$ とし、左の外向き微分を $-2$、右を $-7$ とする。境界値の代わりに最後に平均零を要求する。

まず流束 $q=Au_x$ を未知量にする。$q_x=-f$ を積分すると、任意の定数 $q(0)$ について $q(x)=q(0)-\int_0^x f$ である。さらに $q/A$ を積分し、指定された両端の値または平均に合わせて積分定数を決める。以下では計算後の $u$ について関係を逆に確かめる。

$$
\begin{aligned}
u(x)&=- \frac{x^{3}}{3} + 2 x - \frac{3}{4}\\
u_x&=\frac{d}{d x} \left(- \frac{x^{3}}{3} + 2 x - \frac{3}{4}\right)=2 - x^{2}\\
q&=2 - x^{2}\\
-f&=\frac{d}{d x} \left(2 - x^{2}\right)=- 2 x
\end{aligned}
$$

次に値の条件と流出量を両端で区別して評価する。区間の外法線は左端で $-1$、右端で $1$ である。このため源項の積分は左側の流束から右側の流束を引いたものになる。流束を係数で割った積分は両端の値の差に等しい。

$$
\begin{aligned}
u(0)&=- \frac{3}{4},\qquad x=0\\
u(L)&=- \frac{15}{4},\qquad x=3\\
q(0)&=2,\qquad x=0\\
q(L)&=-7,\qquad x=3\\
\int_0^L f\,\mathrm{d}x&=\int\limits_{0}^{3} 2 x\, dx=9\\
\int_0^L (q/A)\,\mathrm{d}x&=\int\limits_{0}^{3} \left(2 - x^{2}\right)\, dx=-3
\end{aligned}
$$

弱形式も具体的に比べる。試験関数 $\phi=x(L-x)$ は両端で零であり、積分の間に境界項は出ない。以下の値が一致するのは、任意の零トレース試験関数に対しても部分積分できることの一例である。

$$
\begin{aligned}
\phi&=x \left(3 - x\right)\\
\phi^{\prime}&=\frac{d}{d x} x \left(3 - x\right)=3 - 2 x\\
\int_0^L q\phi^{\prime}\,\mathrm{d}x&=\int\limits_{0}^{3} \left(2 - x^{2}\right) \left(3 - 2 x\right)\, dx=\frac{27}{2}\\
\int_0^L f\phi\,\mathrm{d}x&=\int\limits_{0}^{3} 2 x^{2} \left(3 - x\right)\, dx=\frac{27}{2}
\end{aligned}
$$

最後に内部での値、勾配、流束を同じ点で測る。区間全体での勾配の寄与には $A$ の重みが付くので、点での値だけでは全体のエネルギーは分からない。

$$
\begin{aligned}
u&=\frac{7}{12},\qquad x=2\\
u_x&=-2,\qquad x=2\\
q&=-2,\qquad x=2\\
\frac12\int_0^L A(u_x)^2\,\mathrm{d}x&=\int\limits_{0}^{3} \frac{\left(2 - x^{2}\right)^{2}}{2}\, dx=\frac{123}{10}
\end{aligned}
$$

外向き微分の和は $-2-7=-9$、源項の積分は $9$ であり、適合条件が成立する。$q=u_x=2-x^2$ から $u$ を求め、端点値の差は結果として決まる。積分 $\int_0^3(2x-x^3/3)=9/4$ なので定数は $-3/4$ である。$u$ に定数を加えても流束の検算は変わらない。
:::

:::example[Neumann データから定数を除いて戻す]
半径 $2$ の円板で $\partial_\nu u=2\cos\theta+3\cos(3\theta)$ とする。内部の調和解で円板平均が零になるものを求める。

円周データを定数と三角関数に分ける。角振動数 $j$ について半径の方程式は $R^{\prime\prime}+r^{-1}R^{\prime}-j^2r^{-2}R=0$ である。元の条件に合う根を選び、境界半径で係数を規格化する。この例での関数は次のとおりである。

$$
\begin{aligned}
u(r,\theta)&=\frac{r^{3} \cos{\left(3 \theta \right)}}{4} + 2 r \cos{\left(\theta \right)}\\
g&=4 \cos{\left(\theta \right)} + 2 \cos{\left(3 \theta \right)}=2 \left(2 \cos{\left(\theta \right)} + \cos{\left(3 \theta \right)}\right),\qquad r=2\\
u_r&=\frac{\partial}{\partial r} \left(\frac{r^{3} \cos{\left(3 \theta \right)}}{4} + 2 r \cos{\left(\theta \right)}\right)=\frac{3 r^{2} \cos{\left(3 \theta \right)} + 8 \cos{\left(\theta \right)}}{4}\\
u_\theta&=\frac{\partial}{\partial \theta} \left(\frac{r^{3} \cos{\left(3 \theta \right)}}{4} + 2 r \cos{\left(\theta \right)}\right)=- \frac{r \left(3 r^{2} \sin{\left(3 \theta \right)} + 8 \sin{\left(\theta \right)}\right)}{4}
\end{aligned}
$$

二階微分を別々に計算する。半径方向の二階微分だけを角度方向の二階微分と足しても相殺しない。極座標の Laplacian の一階項も入れると各振動数で相殺する。

$$
\begin{aligned}
u_{rr}&=\frac{\partial^{2}}{\partial r^{2}} \left(\frac{r^{3} \cos{\left(3 \theta \right)}}{4} + 2 r \cos{\left(\theta \right)}\right)=\frac{3 r \cos{\left(3 \theta \right)}}{2}\\
u_{\theta\theta}&=\frac{\partial^{2}}{\partial \theta^{2}} \left(\frac{r^{3} \cos{\left(3 \theta \right)}}{4} + 2 r \cos{\left(\theta \right)}\right)=- \frac{r \left(9 r^{2} \cos{\left(3 \theta \right)} + 8 \cos{\left(\theta \right)}\right)}{4}\\
\Delta u&=\frac{3 r \cos{\left(3 \theta \right)}}{2} + \frac{\frac{3 r^{2} \cos{\left(3 \theta \right)}}{4} + 2 \cos{\left(\theta \right)}}{r} - \frac{\frac{9 r^{2} \cos{\left(3 \theta \right)}}{4} + 2 \cos{\left(\theta \right)}}{r}=0
\end{aligned}
$$

境界における値の平均と法線微分を読む。$j\ne0$ の角方向の平均は零であり、全体の平均は定数項である。またエネルギーは体積積分の方法と Green の境界積分の方法の両方で計算できる。

$$
\begin{aligned}
\text{円周の平均}&=\int\limits_{- \pi}^{\pi} \frac{4 \cos{\left(\theta \right)} + 2 \cos{\left(3 \theta \right)}}{2 \pi}\, d\theta=0\\
\partial_\nu u&=2 \cos{\left(\theta \right)} + 3 \cos{\left(3 \theta \right)}\\
\int_{-\pi}^{\pi}|Du|^2\,\mathrm{d}\theta&=\frac{9 \pi r^{4}}{8} + 8 \pi=\frac{\pi \left(9 r^{4} + 64\right)}{8}\\
\int |Du|^2&=\int\limits_{0}^{2} r \left(\frac{9 \pi r^{4}}{8} + 8 \pi\right)\, dr=28 \pi\\
\int_{\partial\Omega} u\partial_\nu u&=\int\limits_{- \pi}^{\pi} \left(2 \cos{\left(\theta \right)} + 3 \cos{\left(3 \theta \right)}\right) \left(8 \cos{\left(\theta \right)} + 4 \cos{\left(3 \theta \right)}\right)\, d\theta=28 \pi
\end{aligned}
$$

次に同じ振動数の半径因子を比較し、角度を固定した一点で値と半径方向の傾きを調べる。角度の取り方によって、振幅のあるモードでも点での値が零になることがある。

$$
\begin{aligned}
u&=\frac{3}{4},\qquad r=1,\theta=\frac{\pi}{3}\\
u_r&=\frac{1}{4},\qquad r=1,\theta=\frac{\pi}{3}\\
u_\theta&=- \sqrt{3},\qquad r=1,\theta=\frac{\pi}{3}
\end{aligned}
$$

法線微分には各角モードで $j/2$ が掛かる。そこで境界値は $4\cos\theta+2\cos(3\theta)$ になる。すべての同心円で平均零なので円板平均も零である。法線微分の平均も零だから適合し、定数を足しても法線微分は変わらない。
:::

### 外部 Neumann 問題における遠方の流束

$n\ge3$ の外部領域 $\R^n\setminus\overline\Omega$ で調和な $u$ が境界まで $C^1$ であり、無限遠で零になるとする。以下の $R$ は $\overline\Omega\subset B_R$ を満たす。内部領域の外法線を $\nu$ と定めると、外部領域の内側境界での法線は $-\nu$ である。従って $\int_{\partial B_R}\partial_r u=\int_{\partial\Omega}\partial_\nu u$ であり、一般にはこの共通値は零でなくてよい。例えば除いた領域に原点が入るとし、$u=c|x|^{2-n}$ とすれば遠方の値は零でも球面積が $R^{n-1}$ で増えるので、総流束は $(2-n)c\sigma_{n-1}$ で一定である。内部 Neumann 問題の適合条件をそのまま外部へ課すと、この解を除いてしまう。

一方、$u=o(|x|^{2-n})$ かつ $Du=o(|x|^{1-n})$ を追加すれば遠方流束が零になり、適合条件が戻る。例えば三次元で単位球外、境界に同じ向きの法線微分 $g=-2$ を与えれば $u=2/|x|$ が解である。境界での総流束は $-8\pi$、半径 $R$ でも $4\pi R^2(-2/R^2)=-8\pi$。固定した法線 $\nu$ で読んだ値と、外部領域の外法線で読んだ値には負号の差がある。

## 演習問題

::::exercise{level=basic}
$0<\alpha<1$ とし、$f(x)=|x|^\alpha$ に対し $-u''=f$、$u(0)=u'(0)=0$ を解け。

:::solution
正の半直線で二度積分すると $u=-x^{2+\alpha}/((1+\alpha)(2+\alpha))$ である。負の半直線には偶延長でき、原点で両側の一階と二階導関数がそれぞれ $0$ に一致する。
:::
::::

::::exercise{level=basic}
$v=x_1^2$ が劣調和であることを示し、原点中心の半径 $2$ の円での平均を求めよ。

:::solution
$\Delta v=2\ge0$ なので劣調和である。平均は $4$ と $\cos^2\theta$ の平均 $1/2$ の積で $2$ であり、$v(0)=0$ より大きい。
:::
::::

::::exercise{level=basic}
円環 $1<r<4$ で $u(1)=0$、$u(4)=2$ の放射対称調和関数を求めよ。

:::solution
$(ru')'=0$ なので $u=A\log r+B$ である。内周から $B=0$、外周から $A=2/\log4$ を得る。従って $u(2)=1$ である。
:::
::::

::::exercise{level=basic}
単位円盤で $-\Delta u=4$、$\partial_\nu u=1$ が可能か判定せよ。

:::solution
適合条件の左辺は $4\pi+2\pi=6\pi\ne0$ である。従って解は存在しない。$f=4$ のまま $h$ を定数にするなら $h=-2$ が必要である。
:::
::::

::::exercise{level=advanced}
一重層の内外の法線微分差が一定密度で $1$ になることを、半径 $2$ の三次元球で検証せよ。

:::solution
外側は $4/r$、内側は連続性と放射対称調和性から $2$ である。外側の法線微分は $-4/r^2$ なので境界で $-1$、内側は $0$ であり、内側から外側を引くと $1$ になる。
:::
::::

::::exercise{level=advanced}
単位円で $g=1-\sin(3\theta)$ を二重層で実現する密度を求めよ。

:::solution
定数に $(K-I/2)$ を作用させると $-1$ 倍、平均 $0$ の関数には $-1/2$ 倍である。従って $\mu=-1+2\sin(3\theta)$ が解であり、この密度から得る調和関数は境界値の一意性により $1-r^3\sin(3\theta)$ である。
:::
::::

::::exercise{level=advanced}
Perron の構成で二つの候補の最小値を取る操作が許されない例を挙げよ。

:::solution
区間 $(-1,1)$ で $v_1=x$、$v_2=-x$ は調和である。最小値 $-|x|$ は $0$ で値 $0$、左右 $\pm r$ の平均 $-r$ である。平均不等式を破るので劣調和でない。二つの候補をまとめるときは最大を用いる。
:::
::::

::::exercise{level=advanced}
連結な有界 Lipschitz 領域で、同じデータを持つ平均 $0$ の弱 Neumann 解が二つあるとする。差の式を使って一意性を示せ。

:::solution
差 $w$ の弱形式は $\int Dw\cdot Dv=0$ である。$w$ 自身は $H^1$ かつ平均 $0$ なので試験に使える。$\norm{Dw}_2^2=0$ となり、平均 $0$ の Poincaré 不等式から $\norm w_2=0$ である。
:::
::::

## 確認テスト

::::quiz
$C^2$ の劣調和関数の Laplace 作用素の符号はどれか。

- [x] $\Delta v\ge0$
- [ ] $\Delta v<0$
- [ ] 常に $\Delta v=0$

:::explanation
球面平均は中心値以上であり、その差の二次の項が $r^2\Delta v/(2n)$ である。従って $\Delta v\ge0$ であり、例えば $v=|x|^2$ では正である。
:::
::::

::::quiz{answer="1" tol="1e-6"}
円環 $1<r<4$ で内周 $0$、外周 $2$ の放射対称解を $r=2$ で評価せよ。

:::explanation
解は $u=2\log r/\log4$ なので $u(2)=1$ である。
:::
::::

::::quiz
外接球の障壁が存在する点で Perron 包絡線はどう振る舞うか。

- [x] その点の連続境界値に近づく
- [ ] どの境界データに対しても $0$ に近づく
- [ ] 領域内からの極限は常に無限大である

:::explanation
障壁の倍数を境界値の上下に足すと包絡線を挟める。点へ近づくと障壁の項が消え、上下の差を任意に小さくできる。
:::
::::

::::quiz{answer="-2" tol="1e-6"}
二次元単位円盤で $-\Delta u=4$ とする。一定 Neumann データ $h$ の適合値を求めよ。

:::explanation
$\int f=4\pi$、$\int h=2\pi h$ なので適合条件から $h=-2$ である。
:::
::::

::::quiz{answer="3" tol="1e-6"}
単位円板で調和な $u(r,\theta)=B r^2\cos2\theta$ を考える。$B=3/2$ のとき $r=1,\theta=0$ の外向き微分を求めよ。

:::explanation
外法線は半径の増加方向なので $\partial_\nu u=\partial_r u=2Br\cos2\theta$ である。指定点では $r=1$、余弦は $1$、従って $2\cdot3/2=3$ である。角方向の微分を用いると違う量になる。
:::
::::


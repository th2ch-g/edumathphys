---
title: Hamilton–Jacobi 理論
description: 作用を未知関数とする方程式を解き、周期運動を作用と角で記述して緩やかな変化に応用する。
notation:
  - symbol: 'S,W'
    meaning: Hamilton の主関数と特性関数
  - symbol: 'I_j,\theta^j'
    meaning: 一周を 2π で規格化した作用と共役な角
  - symbol: '\alpha_j,\beta^j'
    meaning: 完全積分の定数と共役な定数
  - symbol: '\omega_j'
    meaning: 角変数の時間変化率
---

## この章の目標

- 完全積分から正準方程式の軌道を再構成する。
- 調和振動子と Kepler 問題で変数を分離し、積分定数の意味を調べる。
- 閉曲線の作用から周期と周波数を求める。
- 断熱不変量の誤差と、光線近似の成立条件を説明する。

## Hamilton–Jacobi 方程式

### 運動を静止に変える母関数

時間に依存する正準変換で、新しい Hamiltonian を零にできれば、新しい位置も運動量も初期条件を記憶する定数になる。
第二型母関数を $S(q,\alpha,t)$、新運動量を $\alpha=(\alpha_1,\ldots,\alpha_n)$ と置く。
旧運動量は $p_i=\partial_{q^i}S$ なので、$K=H+S_t=0$ が満たすべき偏微分方程式は

$$
S_t+H(q,\partial_q S,t)=0 .
\label{eq-am-hj}
$$

これを :term[Hamilton–Jacobi 方程式]{en="Hamilton–Jacobi equation" yomi="はみるとんやこびほうていしき"}、$S$ を Hamilton の主関数という。
$\partial_q S$ は $n$ 成分の列である。

:::proposition[完全積分による軌道の復元]{#prop-am-hj-complete}
$C^2$ 級の解 $S(q,\alpha,t)$ が $n$ 個の定数を含み、$\det(S_{q^i\alpha_j})\ne0$ とする。
このような解を :term[完全積分]{en="complete integral" yomi="かんぜんせきぶん"} という。
時刻 $t$ ごとに定数 $\beta^j=S_{\alpha_j}$ から $q$ を解き、$p_i=S_{q^i}$ とおけば、局所的な一般解を得る。
:::

:::proof
初期の $(q,p)$ からは混合 Hessian が正則なので、:ref[thm-anl-implicit-function] により $p=S_q$ を $\alpha$ について解け、次に $\beta=S_\alpha$ が決まる。
逆の復元も同じ行列の転置が正則だから可能である。
$\beta$ を時間微分し、方程式を $\alpha_j$ で微分した式と引き算すると $\sum_i S_{\alpha_jq^i}(\dot q^i-H_{p_i})=0$。
正則性により $\dot q^i=H_{p_i}$。
次に $p_i$ を時間微分すると $\dot p_i=S_{q^it}+\sum_jS_{q^iq^j}H_{p_j}$ であり、方程式を $q^i$ で微分した式からこれは $-H_{q^i}$ である。
:::

### 端点と多価性

運動の両端を少し動かして :eqref[eq-am-phase-action] の停留値を比較する。
時刻も動かすと、固定時刻での経路変分は端点の総変位 $\delta q^i$ から $\dot q^i\delta t$ を引いたものである。
部分積分の境界項に積分区間の伸びを加え、

$$
\delta S_{\rm cl}=\left[\sum_i p_i\delta q^i-H\delta t\right]_{1}^{2}
$$

を得る。
初点を固定したとき、終点での勾配は $p$、時刻微分は $-H$ となり、主関数の単位は $\mathrm{J\,s}$ である。
二点を結ぶ軌道が複数あれば停留作用にも枝がある。
転回点や軌道が集まる焦点では単一の $S(q,t)$ がすべての運動量を与えるとは限らない。

:::example[二つの自由粒子解の意味]
一自由度で $S=q\alpha-\alpha^2t/(2m)$ は完全積分であり、$p=\alpha$、$q=\beta+\alpha t/m$。
初点 $q_0$ を固定して $\beta=q_0$ と置き、$S-\alpha q_0$ から $\alpha=m(q-q_0)/t$ を消去すると $S_{\rm cl}=m(q-q_0)^2/(2t)$ を得る。
前者は初期運動量を、後者は初期位置を引数に取る。
$m=2\,\mathrm{kg}$、$q-q_0=3\,\mathrm m$、$t=1\,\mathrm s$ のとき停留作用は $9\,\mathrm{J\,s}$ である。
:::

### 特性曲線として正準方程式を読む

初期の作用 $S_0(q)$ を与える初期値問題では、初期運動量を $p_i(0)=\partial_i S_0(q(0))$ に取る。
そこから正準方程式を解いた曲線を主関数の**特性曲線**という。
曲線に沿っては

$$
\frac{\dd S}{\dd t}=S_t+\sum_iS_{q^i}\dot q^i
=\sum_i p_iH_{p_i}-H.
$$

従って $S_0$ に軌道上の位相空間作用を足せば終点での解が得られる。
出発位置から時刻 $t$ の位置への写像が局所的に逆に解ける間は、この構成で一価な $S(q,t)$ が得られる。
相空間では軌道が交わらなくても、位置への射影では重なるため、局所条件が必要である。

:::example[自由粒子の集束]
$S_0(q_0)=-mq_0^2/(2t_f)$、$t_f>0$ とする。
初速度は $-q_0/t_f$ だから $q(t)=q_0(1-t/t_f)$。
$t<t_f$ なら $q_0$ を消去でき、

$$
S(q,t)=-\frac{mq^2}{2(t_f-t)}
$$

である。
実際 $S_t=-mq^2/[2(t_f-t)^2]$ と $S_q^2/(2m)=mq^2/[2(t_f-t)^2]$ が消し合う。
$t=t_f$ ではすべての位置が原点へ集まり、同じ終点に多くの運動量が到着する。
自由粒子の相流は正則のままでも、位置で表した主関数のこの枝は破綻する。
:::

### 完全積分の定数を選び直す

$\alpha$ を可逆な関数 $\alpha=\alpha(A)$ で置き換えても混合 Hessian の行列式は $\det S_{q\alpha}\det(\partial\alpha/\partial A)$ なので正則性は保たれる。
さらに定数だけの関数 $g(\alpha)$ を $S$ に足しても $p$ と方程式は変わらず、$\beta^j$ が $\partial_{\alpha_j}g$ だけ移る。
特性関数の積分下限を変えたときの初期位相の違いはこの自由度である。

### 初期位相から出発する特性曲線と焦線

完全積分の定数から個々の軌道を得る問題と、波面全体の初期値から主関数を発展させる問題を区別しよう。自由粒子 $H=p^2/(2m)$ に初期関数 $S_0(u)$ を与えると、初期点 $u$ での運動量は $S_0'(u)$ である。特性曲線上で運動量は一定だから
$$
q=u+tS_0'(u)/m,\qquad S=S_0(u)+t[S_0'(u)]^2/(2m).
\label{eq-am-hj-characteristic-map}
$$
右式は $\dd S/\dd t=p\dot q-H$ を軌道に沿って積分したものと分かる。また固定した $q,t>0$ で $S_0(u)+m(q-u)^2/(2t)$ を $u$ について停留させても同じ式になる。$\partial q/\partial u=1+tS_0''(u)/m\ne0$ なら $u$ を消して滑らかな枝が得られる。合成関数を微分すると、停留条件で $u$ の微分が消え、$S_q=m(q-u)/t$、$S_t=-m(q-u)^2/(2t^2)$。従って求めた枝は確かに Hamilton–Jacobi 方程式を満たす。

特性を位置へ投影する写像の臨界値の集まりを :term[焦線]{en="caustic" yomi="しょうせん"} という。一点に集束する場合もあれば、折り目を境に同じ位置へ投影される特性の本数が変わる場合もある。相空間では各点に運動量が付いているので、一意に進む解同士が状態を共有するという矛盾は生じない。多価になった時刻以降も古典力学では各枝を保持する。別の問題で最小作用の枝だけを選ぶ場合には、その選択原理を付け加える必要がある。

:::example[平行な初期束]
$m=1\,\mathrm{kg}$ とし、位置は m、時刻は s、$S_0$ は $\mathrm{J\,s}$ の数値で $S_0(u)=2 u$ とする。微分して初期運動量を定め、時刻 $t$ での位置と作用を初期点で表すと

$$
p=2,\quad q(u,t)=2 t + u,\quad S(u,t)=2 \left(t + u\right).
$$

$u$ 方向に微分した $\partial_u S$ は $p\,\partial_u q$ に等しい。従って投影が可逆なところでは $S_q=p$ に戻る。$t=2\,\mathrm s,u=1\,\mathrm m$ の特性では $q=5\,\mathrm m,p=2\,\mathrm{kg\,m/s},S=6\,\mathrm{J\,s}$、投影の微分は $1$ である。全特性の速さが等しく、投影の微分は時間によらず一である。初期の位相の傾きは伝播後も一様に保たれる。
:::

:::example[一点に集まる二次の位相]
$m=1\,\mathrm{kg}$ とし、位置は m、時刻は s、$S_0$ は $\mathrm{J\,s}$ の数値で $S_0(u)=- \frac{u^{2}}{2}$ とする。微分して初期運動量を定め、時刻 $t$ での位置と作用を初期点で表すと

$$
p=- u,\quad q(u,t)=- u \left(t - 1\right),\quad S(u,t)=\frac{u^{2} \left(t - 1\right)}{2}.
$$

$u$ 方向に微分した $\partial_u S$ は $p\,\partial_u q$ に等しい。従って投影が可逆なところでは $S_q=p$ に戻る。$t=1\,\mathrm s,u=2\,\mathrm m$ の特性では $q=0\,\mathrm m,p=-2\,\mathrm{kg\,m/s},S=0\,\mathrm{J\,s}$、投影の微分は $0$ である。この時刻には出発点の全区間が原点に投影される。原点に至る全特性の作用が一致しても、到着運動量は一つに決まらない。
:::

:::example[折り目を作る三次の位相]
$m=1\,\mathrm{kg}$ とし、位置は m、時刻は s、$S_0$ は $\mathrm{J\,s}$ の数値で $S_0(u)=- \frac{u^{3}}{3}$ とする。微分して初期運動量を定め、時刻 $t$ での位置と作用を初期点で表すと

$$
p=- u^{2},\quad q(u,t)=- u \left(t u - 1\right),\quad S(u,t)=\frac{u^{3} \left(3 t u - 2\right)}{6}.
$$

$u$ 方向に微分した $\partial_u S$ は $p\,\partial_u q$ に等しい。従って投影が可逆なところでは $S_q=p$ に戻る。$t=1\,\mathrm s,u=\frac{1}{2}\,\mathrm m$ の特性では $q=\frac{1}{4}\,\mathrm m,p=- \frac{1}{4}\,\mathrm{kg\,m/s},S=- \frac{1}{96}\,\mathrm{J\,s}$、投影の微分は $0$ である。投影は二次曲線で、選んだ初期点が極値を与える。終点を極値より少し戻すと二つの異なる傾きをもつ解が現れる。
:::

:::example[広がっていく四次の位相]
$m=1\,\mathrm{kg}$ とし、位置は m、時刻は s、$S_0$ は $\mathrm{J\,s}$ の数値で $S_0(u)=\frac{u^{4}}{4}$ とする。微分して初期運動量を定め、時刻 $t$ での位置と作用を初期点で表すと

$$
p=u^{3},\quad q(u,t)=u \left(t u^{2} + 1\right),\quad S(u,t)=\frac{u^{4} \left(2 t u^{2} + 1\right)}{4}.
$$

$u$ 方向に微分した $\partial_u S$ は $p\,\partial_u q$ に等しい。従って投影が可逆なところでは $S_q=p$ に戻る。$t=2\,\mathrm s,u=1\,\mathrm m$ の特性では $q=3\,\mathrm m,p=1\,\mathrm{kg\,m/s},S=\frac{5}{4}\,\mathrm{J\,s}$、投影の微分は $7$ である。正の時刻にはどの初期点でも投影の微分が正である。自由運動でも非一様な運動量が初期面の形を変える。
:::

図では初期運動量を位置 $u$ に比例する $-mu/t_f$ とした。水平線 $t=\text{一定}$ と特性との交点が、その瞬間の粒子の位置である。$t_f$ に近づくほど同じ初期点間隔に対応する位置の間隔が狭くなる。焦点直後は出発順が逆転するが、衝突による運動量の交換は仮定していない。

```tikz caption="一点へ向かう自由粒子の特性。線の傾きは異なり、焦点を通過した後も各粒子の運動量は保存される。"
\draw[->] (-2.5,0)--(2.65,0) node[right] {$q$};
\draw[->] (0,-.2)--(0,3.4) node[above] {$t$};
\foreach \a in {-2,-1,1,2}{
 \draw[thick] (\a,0)--({-.5*\a},3);
 \fill (\a,0) circle (1.2pt);
}
\fill (0,2) circle (1.8pt) node[right] {$t_f$};
```

## 変数分離

### 自律系と調和振動子

$H$ が自律的なら $S=W(q;\alpha)-Et$ として $H(q,\partial_qW)=E$ を解く。
$W$ を :term[特性関数]{en="characteristic function" yomi="とくせいかんすう"} という。
さらに $W=\sum_i W_i(q^i)$ と書き、各 $W_i'$ が自分の座標と定数だけで決まれば偏微分方程式は一変数の積分に帰着する。

$H=p^2/(2m)+m\omega^2q^2/2$ で $E>0$ とする。
$A=\sqrt{2E/(m\omega^2)}$ と置くと $W'=\pm m\omega\sqrt{A^2-q^2}$。
速度が正の枝を原点から積分すると

$$
W_+(q,E)=\frac{m\omega}{2}\left[q\sqrt{A^2-q^2}+A^2\arcsin(q/A)\right].
$$

$\partial_E W_+=\int_0^q m\,\dd x/p(x)=\omega^{-1}\arcsin(q/A)$。
従って $\beta=S_E$ から $q=A\sin[\omega(t+\beta)]$ となる。
$p<0$ の枝では積分の符号と位相の原点を変えて転回点でつなぐ。
一つの枝で $\arcsin$ を単調に増やし続けることはできない。

### Kepler 問題の分離定数

$V=-k/r$ で球座標 $(r,\vartheta,\phi)$ を選ぶ。
$r>0,\ 0<\vartheta<\pi$ での Hamiltonian は

$$
H=\frac1{2\mu}\left(p_r^2+\frac{p_\vartheta^2}{r^2}+\frac{p_\phi^2}{r^2\sin^2\vartheta}\right)-\frac{k}{r}.
$$

方位が循環するので $W=W_r(r)+W_\vartheta(\vartheta)+L_z\phi$ とする。
方程式に $2\mu r^2$ を掛けると角だけの量と動径だけの量が等しい。
共通値を $\ell^2$ と置けば、

$$
W_\vartheta'=\pm\sqrt{\ell^2-L_z^2/\sin^2\vartheta},\qquad
W_r'=\pm\sqrt{2\mu(E+k/r)-\ell^2/r^2}.
$$

従って $E,\ell,L_z$ が分離定数で、$\ell^2=p_\vartheta^2+p_\phi^2/\sin^2\vartheta$ は角運動量の大きさの二乗である。
各平方根が実になる区間で積分し、端で枝をつなぐ。
軌道面を先に選べば二自由度に減り $W=W_r+\ell\phi$ となる。

平面で $p_r>0$ の枝をとり $S_\ell$ と $S_E$ を求めると

$$
\phi-\beta_\ell=\int^r\frac{\ell\,\dd r'}{r'^2p_r},\qquad
t+\beta_E=\int^r\frac{\mu\,\dd r'}{p_r}.
$$

基準点の定数は $\beta$ に吸収した。
第一積分で $u=1/r$ として平方完成すると $p_r^2=\ell^2[(\mu ke/\ell^2)^2-(u-\mu k/\ell^2)^2]$。
従って角を適切にずらせば $u=(\mu k/\ell^2)(1+e\cos(\phi-\phi_0))$ となり :eqref[eq-am-kepler-conic] が戻る。
第二積分は $\dot r=p_r/\mu$ による飛行時間であり、楕円では $r=a(1-e\cos\mathcal E)$ を入れると Kepler の時刻方程式になる。

:::example[二方向で異なる振動数]
$H=\sum_{j=1}^2[p_j^2/(2m)+m\omega_j^2(q^j)^2/2]$ は $E=E_1+E_2$ として上の振動子の $W_+(q^j,E_j)$ の和で分離する。
定数 $\beta_j$ は二つの初期位相に対応する。
例えば $\omega_2=2\omega_1$ なら両位相は $2\pi/\omega_1$ 後にともに戻る。
:::

### 分離してから時刻を一つに戻す

異方的二次元振動子では各 $W_j$ に別の $E_j$ が入るが、時間は一つしかない。
$S=W_1+W_2-(E_1+E_2)t$ とすると、

$$
\partial_{E_1}W_1=t+\beta_1,\qquad
\partial_{E_2}W_2=t+\beta_2 .
$$

一つの時刻に対する位相のずれが二つの初期定数に保存されている。
各自由度に別々の時刻を割り当ててから解を重ねる必要はない。
$E_2=0$ に近づけると二番目の振幅が消え、その位相を初期位置から読み取れなくなる。
完全積分の非退化性はこのような退化点を除く局所的なものなのである。

:::example[一様重力中の平面運動]
$H=(p_x^2+p_z^2)/(2m)+mgz$ とする。
$S=\alpha x+W_z(z)-Et$ と置くと

$$
W_z'=\pm\sqrt{2m(E-mgz)-\alpha^2}.
$$

上昇の枝で積分すれば $W_z=-[2m(E-mgz)-\alpha^2]^{3/2}/(3m^2g)$（定数を除く）。
$\partial_EW_z=-p_z/(mg)$、$\partial_\alpha W_z=\alpha p_z/(m^2g)$ である。
従って $\beta_E=-p_z/(mg)-t$ は $p_z$ の直線的な減少を、$\beta_\alpha=x+\alpha p_z/(m^2g)$ は $\dot x=\alpha/m$ を与える。
頂点 $p_z=0$ を過ぎたら負の枝に移る。
分離で現れた定数は水平方向の運動量と全エネルギーである。
:::

### 分離定数を決める Stäckel 型の条件

変数分離ができるかどうかを、個々の軌道を求めずに調べる方法がある。定常な自然系で、分離関係を
$$
\sum_{j=1}^n A_{ij}(q^i)\alpha_j=\tfrac12p_i^2+U_i(q^i),\qquad i=1,\ldots,n
$$
の形に書けたとする。各行が一つの座標のみに依存し、$\det A\ne0$ である。これを :term[Stäckel 行列]{en="Stäckel matrix" yomi="しゅてっけるぎょうれつ"} を用いた分離という。線形方程式を解いて $\alpha_1=H$ が元の Hamiltonian に一致すれば、$S=\sum_i\int p_i(q^i,\alpha)\dd q^i-\alpha_1t$ は解である。各 $p_i\ne0$ となる枝では $\partial p_i/\partial\alpha_j=A_{ij}/p_i$ なので混合微分の行列式は $\det A/\prod_i p_i$。よって完全積分が得られる。逆に $\alpha_j(q,p)$ を積分と見たとき、$\partial_{p_i}\alpha_j=(A^{-1})_{ji}p_i$、$\partial_{q^i}\alpha_j=(A^{-1})_{ji}(U_i'-\sum_k A_{ik}'\alpha_k)$。同じ $i$ について二つの積分の勾配が比例するため $\{\alpha_j,\alpha_k\}=0$。これは分離が可積分性につながる具体的な仕組みである。

### 一様な外力を受ける Kepler 系の放物線座標

鉛直に一定の力を加えると球対称性が破れ、球面極座標での角運動量の大きさは保存しない。ところが $z$ 軸を外力に平行に取り、$\xi=r+z,\eta=r-z$ とすると $r=(\xi+\eta)/2,\ z=(\xi-\eta)/2,\ x^2+y^2=\xi\eta$ である。軸以外で $\xi,\eta>0$、方位を $\phi$ とすると
$$
\dd l^2=\frac{\xi+\eta}{4\xi}\dd\xi^2+
\frac{\xi+\eta}{4\eta}\dd\eta^2+\xi\eta\,\dd\phi^2 .
$$
これは円筒半径 $\sqrt{\xi\eta}$ と $z$ を微分して二乗和を作れば交差項が消えるためである。$V=-k/r+Fz$ として、$S=W_\xi+W_\eta+p_\phi\phi-Et$ を代入する。$H=E$ の両辺に $\xi+\eta$ を掛けると
$$
\begin{aligned}
2\xi p_\xi^2/m+p_\phi^2/(2m\xi)+F\xi^2/2-E\xi&=\alpha,\\
2\eta p_\eta^2/m+p_\phi^2/(2m\eta)-F\eta^2/2-E\eta&=2k-\alpha .
\end{aligned}
$$
二つの左辺の和が $2k$ になるので、独立な $\xi,\eta$ に対してそれぞれが定数である。$F$ は符号付きで、実際の一様力は $-F\vb e_z$。$\alpha$ の単位は $\mathrm{J\,m}$ である。$W_\xi=\int p_\xi\dd\xi,\ W_\eta=\int p_\eta\dd\eta$ の平方根の符号は初期速度から選ぶ。各転回点までの積分が分かれば $W_E-t,\ W_\alpha,\ W_{p_\phi}+\phi$ が一定という条件から同一時刻の三座標を復元できる。ここで最後の定数の $W$ は $W_\xi+W_\eta$ を指す。

:::example[赤道面を横切る初期値]
質量 $1\,\mathrm{kg}$ とし SI の数値で $(\xi,\eta)=(1,1)$、$(p_\xi,p_\eta,p_\phi)=(1,-1,1)$、$F=2,k=4$ とする。円筒半径は $\sqrt{\xi\eta}=1\,\mathrm m$、高さは $0\,\mathrm m$。逆計量から運動エネルギーを計算すれば

$$
H=\frac{2(\xi p_\xi^2+\eta p_\eta^2)}{m(\xi+\eta)}+\frac{p_\phi^2}{2m\xi\eta}-\frac{2k}{\xi+\eta}+\frac F2(\xi-\eta)=- \frac{3}{2}\,\mathrm J.
$$

この値を各座標の式に戻すと第一の分離定数は $\alpha=5\,\mathrm{J\,m}$、第二の左辺は $3\,\mathrm{J\,m}$。両者の和は $2k$ になり、$\dot\xi=4\xi p_\xi/[m(\xi+\eta)]=2\,\mathrm{m/s}$、$\dot\eta=-2\,\mathrm{m/s}$ が得られる。初め $z=0$ でも $\xi$ と $\eta$ の速度は逆向きなので、その後も赤道面にとどまるとは限らない。
:::

:::example[片方の座標だけが折り返す点]
質量 $1\,\mathrm{kg}$ とし SI の数値で $(\xi,\eta)=(2,1)$、$(p_\xi,p_\eta,p_\phi)=(0,1,1)$、$F=1,k=3$ とする。円筒半径は $\sqrt{\xi\eta}=\sqrt{2}\,\mathrm m$、高さは $\frac{1}{2}\,\mathrm m$。逆計量から運動エネルギーを計算すれば

$$
H=\frac{2(\xi p_\xi^2+\eta p_\eta^2)}{m(\xi+\eta)}+\frac{p_\phi^2}{2m\xi\eta}-\frac{2k}{\xi+\eta}+\frac F2(\xi-\eta)=- \frac{7}{12}\,\mathrm J.
$$

この値を各座標の式に戻すと第一の分離定数は $\alpha=\frac{41}{12}\,\mathrm{J\,m}$、第二の左辺は $\frac{31}{12}\,\mathrm{J\,m}$。両者の和は $2k$ になり、$\dot\xi=4\xi p_\xi/[m(\xi+\eta)]=0\,\mathrm{m/s}$、$\dot\eta=\frac{4}{3}\,\mathrm{m/s}$ が得られる。第一座標の運動量が零なので、その点そのもので混合微分を評価しない。隣接する正則な枝の積分を接続する。
:::

:::example[方位の運動がない子午面]
質量 $1\,\mathrm{kg}$ とし SI の数値で $(\xi,\eta)=(3,1)$、$(p_\xi,p_\eta,p_\phi)=(1,0,0)$、$F=2,k=5$ とする。円筒半径は $\sqrt{\xi\eta}=\sqrt{3}\,\mathrm m$、高さは $1\,\mathrm m$。逆計量から運動エネルギーを計算すれば

$$
H=\frac{2(\xi p_\xi^2+\eta p_\eta^2)}{m(\xi+\eta)}+\frac{p_\phi^2}{2m\xi\eta}-\frac{2k}{\xi+\eta}+\frac F2(\xi-\eta)=1\,\mathrm J.
$$

この値を各座標の式に戻すと第一の分離定数は $\alpha=12\,\mathrm{J\,m}$、第二の左辺は $-2\,\mathrm{J\,m}$。両者の和は $2k$ になり、$\dot\xi=4\xi p_\xi/[m(\xi+\eta)]=3\,\mathrm{m/s}$、$\dot\eta=0\,\mathrm{m/s}$ が得られる。方位角は一定であり、極付近での方位の不定性と実際の回転運動を区別する。
:::

## 作用角変数

### 閉曲線に沿う作用

周期軌道を一周した積分は $W$ の一つの枝の値に依存しない。
正則な一自由度の周期族で運動の向きに積分し、

$$
I(E)=\frac1{2\pi}\oint p\,\dd q,\qquad
\theta=\frac{\partial W(q,I)}{\partial I}\pmod{2\pi}
\label{eq-am-action-angle}
$$

とする。
$I$ を :term[作用変数]{en="action variable" yomi="さようへんすう"}、$\theta$ を :term[角変数]{en="angle variable" yomi="かくへんすう"} という。
自然な振動系では $2\pi I$ は位相図で囲む面積である。
$H(q,p(q,E))=E$ を微分して $\partial_Ep=1/H_p=1/\dot q$。
折り返しでの端点項は両枝の運動量が等しく相殺するので、

$$
\frac{\dd I}{\dd E}=\frac1{2\pi}\oint\frac{\dd q}{\dot q}=\frac{T}{2\pi},\qquad
\dot\theta=\frac{\dd E}{\dd I}=\frac{2\pi}{T},\quad \dot I=0 .
$$

一周で $W$ は $2\pi I$ 増えるため、$\theta$ はちょうど $2\pi$ 増える。
多自由度では独立な閉路 $\gamma_j$ ごとに $I_j=(2\pi)^{-1}\oint_{\gamma_j}\sum_i p_i\dd q^i$ を定める。
こうした閉路を選べる条件は [Liouville–Arnold の定理](/physics/analytical-mechanics/integrability-chaos/) で調べる。

:::example[振動子の明示的な変換]
楕円の面積から $I=E/\omega$。
角の原点を正の速度で原点を通るときに選べば

$$
q=\sqrt{2I/(m\omega)}\sin\theta,\qquad
p=\sqrt{2m\omega I}\cos\theta .
$$

$\partial_\theta q\,\partial_I p-\partial_I q\,\partial_\theta p=\cos^2\theta+\sin^2\theta=1$ なので :ref[prop-am-symplectic] を満たす。
$H=\omega I$ から $\theta(t)=\theta(0)+\omega t$。
$E=0.60\,\mathrm J,\omega=3\,\mathrm{s^{-1}}$ なら $I=0.20\,\mathrm{J\,s}$ である。
:::

```tikz caption="調和振動子の作用は楕円の面積を 2π で割った値である。上の枝は正の速度で進む。"
\draw[->] (-2.65,0) -- (2.7,0) node[right] {$q$};
\draw[->] (0,-1.72) -- (0,1.9) node[above] {$p$};
\draw[thick] (0,0) ellipse (2.1 and 1.3);
\draw[->,thick] (90:2.1 and 1.3) arc[start angle=90,end angle=62,x radius=2.1,y radius=1.3];
\node[above left] at (0,1.3) {$m\omega A$};
\node[below right] at (2.1,0) {$A$};
\node[below left] at (-2.1,0) {$-A$};
```

### Kepler 楕円の動径作用

平面で $\ell>0$ の束縛運動を考える。
角の作用は $I_\phi=\ell$ であり、$I_r=\pi^{-1}\int_{r_-}^{r_+}p_r\dd r$ は正の枝を用いる。
$\partial_E I_r=T/(2\pi)$ は $\ell$ を固定して成り立つ。
$a=-k/(2E)$ と周期を用いると

$$
\partial_E I_r=\frac{k\sqrt{\mu}}{(-2E)^{3/2}},\qquad
I_r=\frac{\mu k}{\sqrt{-2\mu E}}-\ell .
$$

積分定数は円の $E=-\mu k^2/(2\ell^2)$ で動径振幅が消え、$I_r=0$ となることから決まる。
従って

$$
H=-\frac{\mu k^2}{2(I_r+I_\phi)^2}.
$$

二作用に関する偏微分は等しいので動径と方位の周波数は等しく、近点が移らない。
円極限では動径角自体は定まらない。

### Kepler 問題の残る角作用

軌道面を先に固定しなかった場合、方位の正方向を選んで測る作用と極角の往復作用を区別する。
$\ell>|L_z|>0$ とし、方位一周を運動の向きに取ると $I_\phi=|L_z|$ である。
$b=|L_z|/\ell$、$\sin\vartheta_0=b$ と置くと

$$
I_\vartheta=\frac1\pi\int_{\vartheta_0}^{\pi-\vartheta_0}
\sqrt{\ell^2-\frac{L_z^2}{\sin^2\vartheta}}\,\dd\vartheta .
$$

$|L_z|$ を固定して $\ell$ で微分すると端点は運動量が零なので寄与しない。
$u=\cos\vartheta$ に替えれば

$$
\frac{\partial I_\vartheta}{\partial\ell}
=\frac1\pi\int_{-\sqrt{1-b^2}}^{\sqrt{1-b^2}}
\frac{\dd u}{\sqrt{1-b^2-u^2}}=1 .
$$

$\ell=|L_z|$ では軌道は赤道面内にあり往復幅が消えるから $I_\vartheta=\ell-|L_z|$。
従って三次元でも

$$
H=-\frac{\mu k^2}{2(I_r+I_\vartheta+I_\phi)^2}.
\label{eq-am-kepler-actions}
$$

三作用の特定の和しかエネルギーに出ないことは Kepler 問題の縮退を表す。
符号を含めて $I_\phi=L_z$ と取る局所座標なら同じ式の分母で $I_\phi$ を $|I_\phi|$ にする必要がある。

### 周期から作用を求める逆向きの利用

積分 $\oint p\,\dd q$ が直接は難しくても、周期 $T(E)$ が既知なら $I(E)=I(E_0)+(2\pi)^{-1}\int_{E_0}^E T(\epsilon)\dd\epsilon$ を使える。
その逆に実測した周期のエネルギー依存から、位相図の面積がどう増えるかが分かる。

:::example[四次ポテンシャルの尺度]
$H=p^2/(2m)+gq^4,\ g>0$ とする。
振幅 $A=(E/g)^{1/4}$ と置き $q=Ax$ に替えると、

$$
I=\frac{2\sqrt{2m}}{\pi g^{1/4}}E^{3/4}
\int_0^1\sqrt{1-x^4}\,\dd x .
$$

積分値はエネルギーによらない。
従って $T=2\pi I_E$ は $E^{-1/4}$ に比例し、振幅が二倍になればエネルギーは十六倍、周期は半分となる。
調和振動子で振幅に周期が依存しないのは二次ポテンシャルの特別な性質である。
:::

### 振り子の作用を閉軌道の内側から求める

調和振動子では作用がエネルギーに比例した。
振幅が大きくなるとどう変わるかを、長さ $l$ の平面振り子で調べる。
$C=ml^2$、$h=mgl$ とし、

$$
H=\frac{p^2}{2C}+h(1-\cos q).
$$

往復運動では $0<E<2h$。
$k=\sqrt{E/(2h)}$、最大角 $q_0$ を $\sin(q_0/2)=k$ で定めると、正の枝の運動量と作用は

$$
p=2\sqrt{Ch}\sqrt{k^2-\sin^2(q/2)},\qquad
I=\frac{2}{\pi}\int_0^{q_0}p\,\dd q .
$$

$\sin(q/2)=k\sin s$ と置けば独立変数の範囲は常に $0\le s\le\pi/2$ となり、

$$
I=\frac{8\sqrt{Ch}}{\pi}\int_0^{\pi/2}
\frac{k^2\cos^2s}{\sqrt{1-k^2\sin^2s}}\dd s .
$$

ここで第一種と第二種の完全楕円積分を

$$
K(k)=\int_0^{\pi/2}\frac{\dd s}{\sqrt{1-k^2\sin^2s}},\qquad
\mathcal E(k)=\int_0^{\pi/2}\sqrt{1-k^2\sin^2s}\,\dd s
$$

と定める。
二つの被積分関数を通分すれば、

$$
I=\frac{8\sqrt{Ch}}{\pi}
[\mathcal E(k)-(1-k^2)K(k)] .
$$

作用の積分に現れる平方根は、転回点では零になる。
これに対して周期を同じ変数で表すと逆平方根の積分となり、

$$
T=4\int_0^{q_0}\frac{C\,\dd q}{p}=4\sqrt{C/h}\,K(k)
$$

である。
$k\to0$ なら $K(0)=\pi/2$ より $T\to2\pi\sqrt{C/h}$。
作用の被積分関数を $k^2\cos^2s$ まで残せば $I\simeq2\sqrt{Ch}\,k^2=E/\sqrt{h/C}$ となり、小振幅の振動子とも一致する。

逆に $k\to1$ では $s=\pi/2$ 近くの分母が $\sqrt{1-k^2+k^2(\pi/2-s)^2}$ に近づく。
従って $K$ は対数発散する。
しかし作用の積分では分子の $\cos^2s$ があるので発散せず、

$$
I_{\rm sep}=\frac{8\sqrt{Ch}}{\pi}\int_0^{\pi/2}\cos s\,\dd s
=\frac{8\sqrt{Ch}}{\pi}.
$$

セパラトリクスで壊れるのは作用の値自体よりも、その変化率 $I_E=T/(2\pi)$ と、有限の周波数で進む角の表示である。

:::example[回転側では一周の閉路が変わる]
$E>2h$ で $p>0$ の回転なら、$q$ が $2\pi$ 増えるだけで位相円筒上の初点に戻る。
従って

$$
I_{\rm rot}=\frac1{2\pi}\int_{-\pi}^{\pi}\sqrt{2C[E-h(1-\cos q)]}\,\dd q .
$$

回転側から $E\to2h$ とすると被積分関数は $2\sqrt{Ch}\cos(q/2)$ となるので、$I_{\rm rot}\to4\sqrt{Ch}/\pi$。
これは往復運動側の極限の半分である。
往復軌道は正負両方の運動量の枝をたどり、正向き回転では一方の枝だけをたどるからである。
セパラトリクスを越えて同じ閉路のまま断熱不変量を比較することはできない。
:::

### Morse 振動子の作用と解離極限

結合の伸縮を $V=D(1-e^{-a q})^2$（$D,a>0$）で近似すると、井戸の底は調和的だが $q\to\infty$ で $V\to D$ となる。$0<E<D$ のとき $y=e^{-a q}=1+s\cos u,\ s=\sqrt{E/D}$ と置けば、二つの転回点は $y=1\pm s$。$\sqrt{E-V}=\sqrt E\sin u$ と $\dd q=s\sin u\,\dd u/[a(1+s\cos u)]$ を半周期の積分に入れる。作用を直接積分する代わりに $\dd I/\dd E=T/(2\pi)$ を使うと
$$
T=\frac{\sqrt{2m}}{a\sqrt D}\int_0^\pi\frac{\dd u}{1+s\cos u}
=\frac{\pi\sqrt{2m}}{a\sqrt{D-E}} .
$$
二つ目の等号は $v=\tan(u/2)$ によって積分が $\int_0^\infty 2\dd v/[1+s+(1-s)v^2]$ となるためである。$I(0)=0$ としてエネルギーで積分し、逆に解けば
$$
I(E)=\frac{\sqrt{2m}}a(\sqrt D-\sqrt{D-E}),\qquad
H(I)=a\sqrt{2D/m}\,I-\frac{a^2}{2m}I^2 .
\label{eq-am-morse-action}
$$
作用の上限 $I_d=\sqrt{2mD}/a$ は有限でも、$\omega=H'(I)$ はそこで零、周期は無限大となる。軌道の位相面積と一周するのに要する時間を区別できる例である。

:::example[井戸の底の調和近似との比較]
$m=2\,\mathrm{kg},a=1\,\mathrm{m^{-1}},D=4\,\mathrm J$ とする。井戸の底での角振動数は $a\sqrt{2D/m}=2\,\mathrm{s^{-1}}$。$E=\frac{7}{4}\,\mathrm J$ を与えて、エネルギーから二つの位置と作用を別々に求めると

$$
q_-=-\log(1+\frac{\sqrt{7}}{4})\,\mathrm m,\quad q_+=-\log(1-\frac{\sqrt{7}}{4})\,\mathrm m,\quad I=1\,\mathrm{J\,s}.
$$

反転式はこの単位で $H=2I-I^2/4$ なので、$\omega=2-I/2=\frac{3}{2}\,\mathrm{s^{-1}}$、周期は $\frac{4 \pi}{3}\,\mathrm s$。二次方程式のもう一つの根も代数的には $H=E$ を満たすが、$0<I<I_d=4\,\mathrm{J\,s}$ から外れる。近点の間隔が十分短ければ底の二次近似が使えるが、すでに周波数は底での値より小さい。
:::

:::example[井戸の中程での非対称な振動]
$m=2\,\mathrm{kg},a=1\,\mathrm{m^{-1}},D=4\,\mathrm J$ とする。井戸の底での角振動数は $a\sqrt{2D/m}=2\,\mathrm{s^{-1}}$。$E=3\,\mathrm J$ を与えて、エネルギーから二つの位置と作用を別々に求めると

$$
q_-=-\log(1+\frac{\sqrt{3}}{2})\,\mathrm m,\quad q_+=-\log(1-\frac{\sqrt{3}}{2})\,\mathrm m,\quad I=2\,\mathrm{J\,s}.
$$

反転式はこの単位で $H=2I-I^2/4$ なので、$\omega=2-I/2=1\,\mathrm{s^{-1}}$、周期は $2 \pi\,\mathrm s$。正負の折返し距離が等しくない。伸びる向きに柔らかい結合を対称な二次式に置き換えると、位置の時間平均も変わる。
:::

:::example[解離に近い軌道]
$m=2\,\mathrm{kg},a=1\,\mathrm{m^{-1}},D=4\,\mathrm J$ とする。井戸の底での角振動数は $a\sqrt{2D/m}=2\,\mathrm{s^{-1}}$。$E=\frac{63}{16}\,\mathrm J$ を与えて、エネルギーから二つの位置と作用を別々に求めると

$$
q_-=-\log(1+\frac{3 \sqrt{7}}{8})\,\mathrm m,\quad q_+=-\log(1-\frac{3 \sqrt{7}}{8})\,\mathrm m,\quad I=\frac{7}{2}\,\mathrm{J\,s}.
$$

反転式はこの単位で $H=2I-I^2/4$ なので、$\omega=2-I/2=\frac{1}{4}\,\mathrm{s^{-1}}$、周期は $8 \pi\,\mathrm s$。遠い折返し点で長く滞在するため周期が長くなる。作用が有限に見えることを、断熱変化を続けられる根拠にはできない。
:::

### 楕円の軌道要素を正準な三組にする

Kepler 軌道の幾何学的要素には直接観測しやすい意味があるが、任意の三つの角と長さを組み合わせても正準座標にはならない。束縛楕円で $0<e<1$、傾き $0<i<\pi$、まず $L_z>0$ の座標域を取る。動径、極角、方位の作用を $I_r,I_\theta,I_\phi$ とすると :eqref[eq-am-kepler-actions] から
$$
{\cal L}=I_r+I_\theta+I_\phi=\sqrt{\mu k a},\qquad
G=I_\theta+I_\phi=|\boldsymbol\ell|,\qquad N=I_\phi=L_z .
$$
従って $G={\cal L}\sqrt{1-e^2},\ N=G\cos i$。$L_z<0$ の座標域で正の作用 $I_\phi=|L_z|$ を使うなら $N=-I_\phi$ とする。このとき $I_r=\mathcal L-G,\ I_\theta=G+N,\ I_\phi=-N$ なので、第三の共役角は $\theta_\theta-\theta_\phi$ となる。$N=L_z$ を符号付き作用として最初から取る方法でも、同じ節線の角が得られる。三つとも作用の単位をもち、相互の不等式 $|N|<G<{\cal L}$ には傾斜と離心率が残っている。

もとの角を $\theta_r,\theta_\theta,\theta_\phi$ とする。作用に線形な変換を行うので一次形式を展開すると
$$
I_r\dd\theta_r+I_\theta\dd\theta_\theta+I_\phi\dd\theta_\phi
={\cal L}\dd\theta_r+G\,\dd(\theta_\theta-\theta_r)+N\,\dd(\theta_\phi-\theta_\theta).
$$
従って $\theta_r,\theta_\theta-\theta_r,\theta_\phi-\theta_\theta$ が新しい共役角になる。$H=-\mu k^2/(2{\cal L}^2)$ に $G,N$ は現れないので、第一角だけが平均運動 $n=\mu k^2/{\cal L}^3$ で等速に進む。角の原点を軌道要素に対応させるには、分離した特性関数 $W=\int p_r\dd r+\int p_\theta\dd\theta+N\phi$ の各偏微分を調べる。

$G,N$ を固定して $\partial_{\cal L}W$ を取ると、エネルギーを通じて動径項だけが変わる。$\partial_Ep_r=\mu/p_r$ より $\partial_{\cal L}W=n\int\mu\,\dd r/p_r=n(t-t_p)$ であり、近点で零に選べば平均近点角になる。次に $E,N$ を固定すれば $\partial_Gp_r=-G/(r^2p_r)$ で、その積分は軌道面内で近点から測った角の負 $-\nu$。一方 $\partial_Gp_\theta=G/p_\theta$ だから、$\dot\theta=p_\theta/(\mu r^2)$ を使って積分すると $\int G\,\dd t/(\mu r^2)$ となる。昇交点を基準に取ればこれは昇交点から位置までの面内の角である。従って和が近点引数 $g$ を与える。最後に $\partial_Np_\theta=-N/(p_\theta\sin^2\theta)$ の積分は $-\int\dot\phi\,\dd t$ であり、$N\phi$ の微分を加えると昇交点の経度 $h$ が残る。基準点を近点・昇交点の順に選んだとき、三共役角は平均近点角、近点引数、昇交点経度である。これらの組を :term[Delaunay 変数]{en="Delaunay variables" yomi="どろねーへんすう"} という。$e=0$ では近点引数が、$i=0,\pi$ では昇交点経度が定まらない。$L_z=0$ の軌道では球座標によるこの分離の枝が極を通るため、局所表示をそのまま使えない。

:::example[軌道面だけを傾ける]
$\mu=1\,\mathrm{kg},k=1\,\mathrm{J\,m}$ とし、$a=4\,\mathrm m,e=\frac{3}{5},\cos i=\frac{1}{2}$ を与える。まず $r_\pm=a(1\pm e)$ から二つの距離を得る。角運動量の大きさと射影は

$$
({\cal L},G,N)=\left(2,\frac{8}{5},\frac{4}{5}\right)\,\mathrm{J\,s}.
$$

元の積分に戻すと $(I_r,I_\theta,I_\phi)=(\frac{2}{5},\frac{4}{5},\frac{4}{5})\,\mathrm{J\,s}$。角の速さは $(n,0,0)=(\frac{1}{8},0,0)\,\mathrm{s^{-1}}$、$H=- \frac{1}{8}\,\mathrm J$ である。節線の角を変えても三作用は不変である。同じ作用をもつ別の初期状態は、異なる共役角によって区別する。
:::

:::example[長半径を保ち離心率を増やす]
$\mu=1\,\mathrm{kg},k=1\,\mathrm{J\,m}$ とし、$a=4\,\mathrm m,e=\frac{4}{5},\cos i=\frac{1}{2}$ を与える。まず $r_\pm=a(1\pm e)$ から二つの距離を得る。角運動量の大きさと射影は

$$
({\cal L},G,N)=\left(2,\frac{6}{5},\frac{3}{5}\right)\,\mathrm{J\,s}.
$$

元の積分に戻すと $(I_r,I_\theta,I_\phi)=(\frac{4}{5},\frac{3}{5},\frac{3}{5})\,\mathrm{J\,s}$。角の速さは $(n,0,0)=(\frac{1}{8},0,0)\,\mathrm{s^{-1}}$、$H=- \frac{1}{8}\,\mathrm J$ である。同じエネルギーでも近点距離が縮み、角運動量が小さくなる。近点から近点までの時間は変わらない。
:::

:::example[軌道を相似に拡大する]
$\mu=1\,\mathrm{kg},k=1\,\mathrm{J\,m}$ とし、$a=9\,\mathrm m,e=\frac{3}{5},\cos i=\frac{1}{2}$ を与える。まず $r_\pm=a(1\pm e)$ から二つの距離を得る。角運動量の大きさと射影は

$$
({\cal L},G,N)=\left(3,\frac{12}{5},\frac{6}{5}\right)\,\mathrm{J\,s}.
$$

元の積分に戻すと $(I_r,I_\theta,I_\phi)=(\frac{3}{5},\frac{6}{5},\frac{6}{5})\,\mathrm{J\,s}$。角の速さは $(n,0,0)=(\frac{1}{27},0,0)\,\mathrm{s^{-1}}$、$H=- \frac{1}{18}\,\mathrm J$ である。長さ比の平方根が三作用の共通の倍率、周期の倍率は長さ比の $3/2$ 乗である。
:::

## 断熱不変量

### 一周期の間にどれだけ変わるか

外部パラメータ $\lambda$ をゆっくり変える。
周期 $T$ が変化時間に比べて短く、滑らかな周期軌道の族にとどまるとき、作用は :term[断熱不変量]{en="adiabatic invariant" yomi="だんねつふへんりょう"} になる。
ここで「断熱」は力学的な時間尺度の分離を指す。

各時刻で凍結した $\lambda$ の母関数を $W(q,I;\lambda)$ とすると、時間依存の正準変換から $K=H_0(I,\lambda)+\dot\lambda\,G(\theta,I,\lambda)$、
$G=(\partial_\lambda W)_{q,I}$ となる。
$W$ の一周の増分 $2\pi I$ は固定した $I$ で $\lambda$ に依存しないから $G$ は角について周期的である。
ゆえに $\dot I=-\dot\lambda G_\theta$ の一次の一周期平均は零となる。

より詳しく、無次元化した時刻で $\lambda=\lambda(\epsilon t)$ とし、$\omega=\partial_IH_0$ が零から離れ、
必要な導関数が有界な正則領域にいるとする。
$\dot\theta=\omega+O(\epsilon)$ なので

$$
\frac{\dd}{\dd t}\left(I+\frac{\dot\lambda}{\omega}G\right)
=-\dot\lambda G_\theta+\frac{\dot\lambda}{\omega}G_\theta\omega+O(\epsilon^2)=O(\epsilon^2).
\label{eq-am-adiabatic-action}
$$

$\dot I=O(\epsilon)$、$\ddot\lambda=O(\epsilon^2)$ を余りに用いた。
従って $O(1/\epsilon)$ の時間にも $I$ の変化は $O(\epsilon)$ に収まる。
セパラトリクスへの接近では周期が発散し、この見積もりの分母が小さくなる。

:::example[ゆっくり硬くするばね]
$H=p^2/(2m)+m\omega(t)^2q^2/2$ では $\dot E=m\omega\dot\omega q^2$。
凍結した角表示を入れると

$$
\frac{\dd}{\dd t}(E/\omega)=-\frac{I\dot\omega}{\omega}\cos 2\theta.
$$

振動一周期の正負が相殺される。
$|\dot\omega|/\omega^2\ll1$ で $\omega$ が2倍になれば $E$ はほぼ2倍、振幅は $A=\sqrt{2I/(m\omega)}$ によりほぼ $1/\sqrt2$ 倍になる。
急に変えた場合、例えば $q=0$ の瞬間には初期エネルギーは変わらず、この予測から外れる。
:::

### ゆっくり変えるとは何との比較か

単振り子の長さ $l$ をゆっくり変える。
振幅が小さい範囲で角度 $\alpha$ を座標とすると、凍結した Hamiltonian は

$$
H=\frac{p_\alpha^2}{2ml^2}+\frac{mgl}{2}\alpha^2,\qquad
\omega=\sqrt{g/l}.
$$

振幅を $\alpha_{\max}$ とすれば

$$
I=\frac{m\sqrt{g}}2\,l^{3/2}\alpha_{\max}^2 .
$$

従って $I$ の近似保存から $\alpha_{\max}\propto l^{-3/4}$、横の振れ幅 $l\alpha_{\max}\propto l^{1/4}$ が得られる。
長さを増やすと角度は小さくなるが、長さで測った振幅はむしろ大きくなりうる。
条件は $|\dot l|/(l\omega)\ll1$ で、初めから大角度でセパラトリクス近くに振れる場合は同じ評価を使えない。

:::example[急変と緩変を同じ初期状態で比べる]
振動子を最大変位 $A$ で静止させ、質量を固定して角振動数を急に二倍にする。
$q,p$ は直後に連続で $E'=m(2\omega)^2A^2/2=4E$ だから $I'=2I$。
これを $q=0$ で行えば全エネルギーは運動エネルギーだけなので $E'=E,\ I'=I/2$。
ゆっくり行ったときの $I'\simeq I$ とはどちらも異なる。
急変の結果は開始位相に依存する。
:::

### 断熱変化の仕事率を作用の微分から求める

パラメータ $\lambda$ を固定した周期軌道の作用を $I(E,\lambda)$ とする。
$H(q,p;\lambda)=E$ を、$q,E$ を固定して $\lambda$ で微分すると

$$
\left(\frac{\partial p}{\partial\lambda}\right)_{q,E}
=-\frac{H_\lambda}{H_p}.
$$

往復積分の境界項は二つの枝で相殺されるため、

$$
I_\lambda=-\frac{1}{2\pi}\oint H_\lambda\,\frac{\dd q}{H_p}
=-\frac{T}{2\pi}\langle H_\lambda\rangle_E .
$$

ここで平均は凍結した一周期について取る。
既知の $I_E=T/(2\pi)$ と合わせて、作用を一定とする方向のエネルギー変化は

$$
\left(\frac{\partial E}{\partial\lambda}\right)_I
=-\frac{I_\lambda}{I_E}=\langle H_\lambda\rangle_E .
$$

実際の駆動で $\dot H=H_\lambda\dot\lambda$ だから、外部変化が一周期内でほぼ一定ならこの平均仕事率が得られる。
作用の近似保存と力の一周期平均は、別の原理を仮定せずに導ける。

:::example[ばね定数をゆっくり増やす]
$\lambda=k$、$H=p^2/(2m)+kq^2/2$ なら $H_k=q^2/2$。
振幅 $A$ の周期平均は $\langle q^2\rangle=A^2/2=E/k$ である。
従って $I$ 一定で $\partial E/\partial k=E/(2k)$。
積分すると $E\propto\sqrt k$ が戻る。
$k$ を少しずつ $1.21$ 倍にすると予測は $E'=1.10E$。
ばねが自然長にある一瞬の $H_k$ は零なのに有限のエネルギー増加があるのは、圧縮・伸長の位相も合わせて平均したからである。
:::

### ゆっくり閉路を回った後の角のずれ

作用が初期値に戻っても、パラメータの閉路に沿って角のずれが蓄積することがある。一自由度の正則な周期族を取り、$\lambda=\lambda(\epsilon t)$ を滑らかに変え、閉路の始点と終点で $\dot\lambda=0$ とする。前節の $K=H_0+\sum_j\dot\lambda_jG_j$ から角依存の一次項を除くには、第二型母関数 $\theta J+S_1(\theta,J,t)$ に対して
$$
\partial_\theta S_1=-\frac{\sum_j\dot\lambda_j(G_j-\langle G_j\rangle_\theta)}{\omega(J,\lambda)}
$$
と選ぶ。右辺は角平均が零なので周期的な $S_1=O(\epsilon)$ が存在する。$\partial_t S_1=O(\epsilon^2)$ であり、変換後の Hamiltonian は $H_0(J,\lambda)+\sum_j\dot\lambda_j\langle G_j\rangle_\theta+O(\epsilon^2)$ となる。周波数が零から離れ、必要な導関数が一様に有界なら、二次の振動部分にも同じ操作を適用できる。残る角依存は三次であり、一周の駆動時間 $O(1/\epsilon)$ に平均化した作用の変動は $O(\epsilon^2)$ に収まる。従って初期作用で評価した $\int\omega(J,\lambda)\dd t$ を角から差し引くと、極限では $\partial_J\oint\sum_j\langle G_j\rangle_\theta\dd\lambda_j$ が残る。始点で $\dot\lambda=0$ としたので、一次変換による初期作用のずれはない。

この式を軌道の幾何から計算する。$\theta=W_I$、$p=W_q$ なので
$$
p\,\dd q=I\,\dd\theta+\dd(W-I\theta)-\sum_j G_j\,\dd\lambda_j .
$$
$W$ と $I\theta$ は一周で同じ $2\pi I$ だけ増えるため、その差は周期的である。固定した $I,\theta$ で $\lambda_j$ の係数を比較して角平均を取ると、$\mathcal A_j=\langle p\partial_{\lambda_j}q\rangle_\theta=-\langle G_j\rangle_\theta+\partial_{\lambda_j}\langle W-I\theta\rangle_\theta$。最後の項は、始点と終点で同じ角と作用の規約を使う閉路では積分すると消える。従って :term[Hannay 角]{en="Hannay angle" yomi="はんねいかく"} は
$$
\Delta\theta_{\rm H}=-\partial_I\oint\sum_j\mathcal A_j\,\dd\lambda_j
$$
である。以下の例のように $p\partial_Iq$ の平均が零でない座標でも、この導出は成り立つ。

無次元の模型 $q=c\sqrt{2I}\sin\theta,\ p=\sqrt{2I}\cos\theta/c+bq,\ c>0$ を取る。瞬間の Hamiltonian は $H=\omega[c^2(p-bq)^2+q^2/c^2]/2=\omega I$。平均で $\langle\sin^2\theta\rangle=1/2,\ \langle\sin\theta\cos\theta\rangle=0$ を使えば ${\cal A}_c=bIc,\ {\cal A}_b=0,\ \langle p q_I\rangle=bc^2/2$。長方形の閉路で $b=0$ のまま $c_1$ から $c_2$、次に $b=B$、その値で $c_1$ に戻し、最後に $b=0$ に戻すと、$\Delta\theta_{\rm H}=B(c_2^2-c_1^2)/2$。$c_1=1,c_2=2,B=1/3$ なら動的位相 $\int\omega\dd t$ のほかに $0.50$ ラジアンが残る。向きを逆にするとこの符号だけが反転する。角度は $2\pi$ を法として読む。

### 反射する壁をゆっくり動かしたときの作用

断熱不変性の近似を、瞬間衝突を含む運動で直接確かめよう。左壁を $x=0$ に固定し、右壁を $L(t)$ とする。粒子が右壁に到着した瞬間の長さを $L$、入射速度を $v>0$、壁速度を一定 $u$ とし $v>2|u|$ を仮定し、壁の後退時にはさらに $v>3u$ とする。完全弾性反射の直後に左へ戻る速さは $v'=v-2u$。左壁を経て次に右壁と会うまでの時間 $\Delta t$ は $v'\Delta t=2L+u\Delta t$ を満たす。従って次の右壁衝突直前の長さは $L'=L(v'+u)/(v'-u)$。箱の瞬間作用は面積から $I=mvL/\pi$ であるから
$$
\frac{I'}{I}=\frac{(v-2u)(v-u)}{v(v-3u)}
=1+\frac{2u^2}{v(v-3u)}.
$$
後退する壁に再び追いつくには $v>3u$ が必要である。速度の一次で粒子は冷えるが、壁間隔の一次の増加がこれを補い、作用の衝突間の変化は二次から始まる。伸縮に要する時間に含まれる衝突数は $v/|u|$ の桁であり、その間の累積誤差は一般には一次の桁になる。伸長の端で突然停止する場合には、最後の衝突から停止までの位相も最終速度に影響する。瞬間瞬間のエネルギーを滑らかな関数とみなせない模型でも、尺度を指定すれば断熱変化の意味を定められる。

:::example[伸びる箱]
右壁に到着したとき $L=1\,\mathrm m,\ v=4\,\mathrm{m/s},u=\frac{1}{10}\,\mathrm{m/s}$ とする。衝突により速度の大きさは $\frac{19}{5}\,\mathrm{m/s}$ に変わる。次の到着までの時間は $\frac{20}{37}\,\mathrm s$、その時点で右壁は $L'=\frac{39}{37}\,\mathrm m$ にある。従ってエネルギー比は $(v'/v)^2=\frac{361}{400}$ だが、作用比は $I'/I=\frac{741}{740}$。作用の相対誤差は約 $0.135\%$ である。壁の速度が同じ実験でも入射する粒子の速さが違えば、断熱性を支配する $|u|/v$ の値は違う。
:::

:::example[縮む箱]
右壁に到着したとき $L=1\,\mathrm m,\ v=4\,\mathrm{m/s},u=- \frac{1}{10}\,\mathrm{m/s}$ とする。衝突により速度の大きさは $\frac{21}{5}\,\mathrm{m/s}$ に変わる。次の到着までの時間は $\frac{20}{43}\,\mathrm s$、その時点で右壁は $L'=\frac{41}{43}\,\mathrm m$ にある。従ってエネルギー比は $(v'/v)^2=\frac{441}{400}$ だが、作用比は $I'/I=\frac{861}{860}$。作用の相対誤差は約 $0.116\%$ である。
:::

:::example[遅い粒子についての誤差]
右壁に到着したとき $L=1\,\mathrm m,\ v=1\,\mathrm{m/s},u=\frac{1}{10}\,\mathrm{m/s}$ とする。衝突により速度の大きさは $\frac{4}{5}\,\mathrm{m/s}$ に変わる。次の到着までの時間は $\frac{20}{7}\,\mathrm s$、その時点で右壁は $L'=\frac{9}{7}\,\mathrm m$ にある。従ってエネルギー比は $(v'/v)^2=\frac{16}{25}$ だが、作用比は $I'/I=\frac{36}{35}$。作用の相対誤差は約 $2.857\%$ である。
:::

## 光学と力学の対応

### 等作用面と光線

質量 $m$、等方的な運動エネルギー $p^2/(2m)$、時間不変の $V$ に対する定常方程式は

$$
|\nabla W|=\sqrt{2m(E-V)} .
$$

$E>V$ の領域で運動量は等 $W$ 面の法線方向にあり、運動曲線に沿って $\dd W=p\,\dd l$。
短縮作用の停留は $\delta\int p\,\dd l=0$ と書ける。
一方、等方的で静止した屈折率 $n(\vb r)$ の媒質で光速が $c/n$ なら、Fermat の原理は $\delta\int n\,\dd l=0$。
光路長 $\mathcal L$ の :term[アイコナール方程式]{en="eikonal equation" yomi="あいこなーるほうていしき"} は $|\nabla\mathcal L|=n$ となる。
力学の $p$ と光学の $n$ は光路の形を決める同じ位置に現れるが、単位は異なる。

波長が屈折率の変化尺度より短いと仮定し、スカラー波の Helmholtz 方程式 $\nabla^2u+k_0^2n^2u=0$ に $u=Ae^{ik_0\mathcal L}$ を代入する。
$k_0^2$ の係数は $A(n^2-|\nabla\mathcal L|^2)$、$ik_0$ の係数は $2\nabla A\cdot\nabla\mathcal L+A\nabla^2\mathcal L$。
前者が光線を、後者が振幅の輸送を決める。
焦線付近では振幅の緩変化が破れ、光線だけで回折を記述できない。

:::example[平面境界での屈折]
平らな界面の両側で屈折率が一定なら光路はそれぞれ直線である。
交点を界面に沿って $\delta x$ 動かすと、光路長の変化は $(n_1\sin i-n_2\sin r)\delta x$。
停留条件は Snell の法則となる。
$n_1=1,n_2=1.5,i=30^\circ$ なら $\sin r=1/3$。
力学でも境界に平行な座標が循環し、$p_1\sin i=p_2\sin r$ が保存する。
:::

### 波面の法線から光線の曲率を求める

等方的媒質で光路を弧長 $l$ で測り、単位接ベクトルを $\vb t=\dd\vb r/\dd l$ とする。
任意のパラメータ $s$ を使って光路を $\int n(\vb r)|\vb r'|\dd s$ と書き変分すると

$$
\frac{\dd}{\dd l}(n\vb t)=\nabla n .
$$

接ベクトルと平行な部分は $\dd n/\dd l=\vb t\cdot\nabla n$ の恒等式である。
垂直部分は

$$
n\frac{\dd\vb t}{\dd l}=\nabla n-\vb t(\vb t\cdot\nabla n)
$$

となり、屈折率が速く増える側へ曲がることを示す。
例えば $n=n(z)$ だけに依存する層状媒質では $n\,\dd x/\dd l$ が一定で、法線と光線の角を $i$ として $n\sin i$ が連続的に保存する。
境界の屈折則が滑らかな媒質にもつながる。

### 波の位相と粒子の作用の次元

$W$ は作用の単位をもち、光路長 $\mathcal L$ は長さである。
両者を位相として比べる場合には無次元化する尺度が要る。
単色光の位相は $k_0\mathcal L$、非相対論的粒子の波の位相は $W/\hbar$ である。
後者の対応は定常 Schrödinger 方程式に $\psi=Ae^{iW/\hbar}$ を代入すると直接分かる。
二階微分を積の公式で展開すると

$$
\left(\frac{|\nabla W|^2}{2m}+V-E\right)A
-\frac{i\hbar}{2m}(2\nabla A\cdot\nabla W+A\nabla^2W)
-\frac{\hbar^2}{2m}\nabla^2A=0 .
$$

$A$ を実振幅にとり、$A,W$ が滑らかで位相が振幅に比べて速く変わる領域では最初の括弧が先頭の近似で零となり、次の次数で $\nabla\cdot(A^2\nabla W)=0$ が得られる。
転回点で $p=|\nabla W|$ が零に近づくと位相が速いという仮定が崩れる。
転回点をまたぐ波動関数の接続には波動方程式を直接扱う必要がある。
[量子力学](/physics/quantum-mechanics/)

### 位相速度と力学的な進行時間

定常な特性関数の定義には $E$ が自由なパラメータとして含まれる。
そのため $p\,\dd l$ の積分を $E$ で微分すれば、軌道の形から移動に要する時間も引き出せる。
二端点を固定し、エネルギーごとに短縮作用を停留させる軌道を選ぶとする。
積分経路自体の変化は停留条件で一次では消えるから

$$
\left(\frac{\partial W}{\partial E}\right)_{\rm endpoints}
=\int\frac{\partial}{\partial E}\sqrt{2m(E-V)}\,\dd l
=\int\frac{m}{p}\,\dd l.
$$

$m/p$ は速さの逆数であり、右辺は飛行時間に等しい。
角の分離定数を固定して $W_r$ を微分したときにも同じ理由で時間積分が得られた。

:::example[一様な領域における波と粒子]
$V=V_0$ 一定とし $\psi\propto e^{i(kx-\omega t)}$ の成分を考える。
$p=\hbar k,\ E=\hbar\omega$ とすれば

$$
\omega(k)=\frac{\hbar k^2}{2m}+\frac{V_0}{\hbar}.
$$

狭い波数幅の二成分を足すと包絡線の位相は $(\Delta k)x-(\Delta\omega)t$ であり、移動速度は群速度 $\dd\omega/\dd k=\hbar k/m=p/m$ になる。
一方、個々の山の位相速度は $\omega/k=p/(2m)+V_0/p$。
ポテンシャルの基準を変えると位相速度の式は変わるが、群速度と力学的な飛行時間は変わらない。
光線と力学軌道の対応で、何を追跡しているのかを区別する例である。
:::

### 分離の可否は座標に依存する

平面の中心力では極座標にすると方位が循環するが、非自明な $V(\sqrt{x^2+y^2})$ を直交座標で $W_x(x)+W_y(y)$ に分けられるとは限らない。
実際、加法分離できるなら

$$
2mV(x,y)=2mE-[W_x'(x)]^2-[W_y'(y)]^2
$$

が必要なので $\partial_x\partial_y V=0$ である。
一つの長方形領域で混合微分が零なら $V$ は $V_x(x)+V_y(y)$ に分かれ、逆に二つのエネルギーへ分離できる。
例えば等方振動子は直交座標でも極座標でも分離する。
$V=-k/\sqrt{x^2+y^2}$ では $\partial_x\partial_yV=-3kxy/(x^2+y^2)^{5/2}$ が一般に零でないから、
この直交座標での加法分離は失敗する。
それでも極座標で解けることはすでに示した。
座標選択による分離の失敗を、系が可積分でないという判定に使うことはできない。

### 短波長極限での振幅輸送と転回点

波面の運動だけを追っていると焦線で何が破れるか分かりにくい。一次元の定常波動方程式 $-\hbar^2\psi''/(2m)+V\psi=E\psi$ に $\psi=Ae^{iW/\hbar}$ を代入して階数を分ける。実部と虚部から
$$
(W')^2=2m(E-V)+\hbar^2 A''/A,\qquad (A^2W')'=0
$$
を得る。$A$ が波長に比べてゆっくり変わる領域で最後の曲率補正を落とせば $W'=\pm p$、$A\propto p^{-1/2},\ p=\sqrt{2m(E-V)}$。透過する定常流の確率流束は $A^2 W'/m$ で一定である。近似を点検するには、第一近似から $\hbar^2|A''/A|\ll p^2$ を調べればよい。例えば運動量が非零でも急に変化する場所では条件が崩れる。

単純転回点 $x_*$、$V'(x_*)=F>0$ では $V-E\simeq F(x-x_*)$。$\ell=(\hbar^2/(2mF))^{1/3}$、$u=(x-x_*)/\ell$ と置くと局所方程式は $\psi_{uu}=u\psi$ となる。有限な Airy 解でつなぐ必要があり、$p^{-1/2}$ の発散は波そのものの発散を表さない。二つの単純転回点とその外で減衰する条件を課すと、各端での Airy 接続の位相 $\pi/4$ から、半古典近似の精度で $\int_{x_-}^{x_+}p\,\dd x\simeq(n+1/2)\pi\hbar$ を得る。位相の由来は、減衰解の積分表示 $\operatorname{Ai}(-s)=\pi^{-1}\int_0^\infty\cos(t^3/3-st)\dd t$ の停留点 $t=\sqrt{s}$ で二次の位相を積分する際の $+\pi/4$ である。従って左からの振動解は $\sin(\int_{x_-}^x p\,\dd x/\hbar+\pi/4)$、右からは $\sin(\int_x^{x_+}p\,\dd x/\hbar+\pi/4)$ に比例する。一つの定数倍で全区間で一致する条件が上の式である。硬い壁には別の境界位相があり、この半整数をそのまま使えない。

### トンネル透過を決める虚の作用

$E<V$ では実運動量で記した古典軌道はないが、波動方程式の解はその領域を越えて続く。$\kappa(x)=\sqrt{2m(V-E)}/\hbar$ とおくと、前節で得た位相を虚数に延長した減衰解は、十分滑らかな領域で $\psi(x)\propto\kappa(x)^{-1/2}\exp[-\int\kappa(x)\dd x]$。この表示の有効範囲は $|(\kappa^{-1/2})''|\ll\kappa^2\kappa^{-1/2}$ である。二つの単純転回点に挟まれた厚い障壁なら、振幅が指数で小さくなるので、透過流の主な指数因子は $\exp[-2\int_{x_1}^{x_2}\kappa(x)\dd x]$ となる。係数まで求めるには両側で進行波と反射波を接続する必要がある。

段差のある長方形障壁でこの区別を具体的に確かめられる。外側の波数を $k$、内部の減衰定数を $\kappa$、幅を $d$ とする。左境界から右境界へ $(\psi,\psi')$ を運ぶ行列は内部の二解 $e^{\pm\kappa x}$ から
$$
\begin{pmatrix}\cosh\kappa d&\sinh(\kappa d)/\kappa\\
\kappa\sinh\kappa d&\cosh\kappa d\end{pmatrix}.
$$
右境界で $(t,ikt)$ とし、逆行列を掛けて左の $(1+r,ik(1-r))$ に合わせる。第一成分と第二成分を $ik$ で割ったものの和を取れば $2=t[2\cosh\kappa d+i(\kappa/k-k/\kappa)\sinh\kappa d]$。左右の外部媒質が同じなら流束比は $|t|^2=[1+(k^2+\kappa^2)^2\sinh^2(\kappa d)/(4k^2\kappa^2)]^{-1}$ である。例えば $k=\kappa,\ \kappa d=3$ なら厳密値は $1/\cosh^2 3\simeq0.00987$、厚い障壁の近似は $4e^{-6}\simeq0.00992$。指数の $e^{-6}$ だけを透過確率とすれば係数四を落とす。古典的な実軌道から離れる場合にも、作用が与える量と境界条件が与える量を区別できる。

## 演習問題

::::exercise{level=basic}
一様な力 $F$ のもとで $H=p^2/(2m)-Fq$ とする。
$S=q(\alpha+Ft)+s(t)$ を代入し $s$ と軌道を求めよ。

:::solution
$H+S_t=(\alpha+Ft)^2/(2m)+s'$ なので $s=-(\alpha^2t+\alpha Ft^2+F^2t^3/3)/(2m)$ と選ぶ。
$\beta=S_\alpha=q-\alpha t/m-Ft^2/(2m)$、よって $q=\beta+\alpha t/m+Ft^2/(2m)$、$p=\alpha+Ft$ である。
:::
::::

::::exercise{level=basic}
$E=4\,\mathrm J,\omega=8\,\mathrm{s^{-1}}$ の振動子の作用と周期を求めよ。

:::solution
$I=E/\omega=0.5\,\mathrm{J\,s}$。
$H=\omega I$ より一周に要する時間は $2\pi/\omega=\pi/4\,\mathrm s$ である。
:::
::::

::::exercise{level=basic}
慣性モーメント $C$ の自由回転子で正向き回転の $I$ と $H(I)$ を求めよ。

:::solution
$p_\phi$ は一定だから $I=(2\pi)^{-1}\int_0^{2\pi}p_\phi\dd\phi=p_\phi$。
従って $H=I^2/(2C)$、角速度は $I/C$ である。
:::
::::

::::exercise{level=basic}
調和振動子の角振動数をゆっくり4倍にする。
エネルギーと振幅は何倍か。

:::solution
$E/\omega$ がほぼ一定だからエネルギーは4倍。
$m\omega^2A^2/2=E$ より、質量一定なら振幅は $\sqrt4/4=1/2$ 倍である。
:::
::::

::::exercise{level=advanced}
$S=q\alpha^2-\alpha^4t/(2m)$ は自由粒子の完全積分としてどこで使えるか。

:::solution
$S_t+(S_q)^2/(2m)=0$。
混合微分は $2\alpha$ なので $\alpha\ne0$ で局所的に使える。
$p=\alpha^2$ のため正の運動量だけを覆い、全体では $\alpha$ の符号まで一意には決まらない。
:::
::::

::::exercise{level=advanced}
Kepler 楕円で角運動量を固定して $E$ を円の値から $0$ に近づける。
動径作用と周期の極限を求めよ。

:::solution
$I_r=\mu k/\sqrt{-2\mu E}-\ell$ は発散し、$T=2\pi\partial_E I_r=2\pi k\sqrt\mu/(-2E)^{3/2}$ も発散する。
遠点が無限遠へ移る放物線極限に対応する。
:::
::::

::::exercise{level=advanced}
長さ $d$ の一次元箱で質量 $m$ の粒子が壁に弾性反射する。
$I$ を求め、緩やかな圧縮で $E d^2$ がほぼ一定になる理由を説明せよ。

:::solution
往路と復路で運動量の符号と移動方向がともに変わり $\oint p\,\dd q=2d\sqrt{2mE}$。
従って $I=d\sqrt{2mE}/\pi$。
動く右壁の速度を $u$ とすると衝突後の速度は $2u-v$ で、一衝突で $\Delta E=-2muv+O(u^2)$。
往復時間 $2d/v$ で割れば $\dot E\simeq-muv^2/d=-2E\dot d/d$ となり $Ed^2$ は一次まで一定である。
:::
::::

::::exercise{level=advanced}
一定の $n$ に対し $\mathcal L=n|\vb r|$ がアイコナール方程式を満たす範囲と光線の形を求めよ。

:::solution
$\vb r\ne0$ で $\nabla\mathcal L=n\vb r/|\vb r|$ の長さは $n$。
法線は放射方向なので光線は原点からの半直線である。
原点で距離関数は微分できず、この表示の通常の解の範囲外となる。
:::
::::

::::exercise{level=advanced}
端点 $q_1,q_2$ を時間 $T$ で結ぶ調和振動子の主関数を求めよ。
$\sin\omega T\ne0$ とする。

:::solution
運動方程式上では $\dd(q\dot q)/\dd t=\dot q^2-\omega^2q^2$ なので $S_{\rm cl}=m[q\dot q]_0^T/2$。
$q(t)=q_1\cos\omega t+b\sin\omega t$ に終点を入れると $b=(q_2-q_1\cos\omega T)/\sin\omega T$。
両端の速度を代入して

$$
S_{\rm cl}=\frac{m\omega}{2\sin\omega T}
[(q_1^2+q_2^2)\cos\omega T-2q_1q_2].
$$

$\partial_{q_2}S_{\rm cl}=m\dot q(T)$、$-\partial_{q_1}S_{\rm cl}=m\dot q(0)$ となる。
$T=\pi/\omega$ では $q_2=-q_1$ が初速にかかわらず成立し、与えた端点から初速が一意に戻らない。
:::
::::

::::exercise{level=advanced}
$H=p^2/(2m)+g|q|^s$（$g,s>0$）で作用の $E,g,m$ に関するべき依存を求めよ。
また $g$ をゆっくり変えたときの $E$ の変化を求めよ。

:::solution
転回点は $A=(E/g)^{1/s}$。
往復の対称性を使って

$$
I=\frac{2}{\pi}\int_0^A\sqrt{2m(E-gq^s)}\,\dd q
=C_s\sqrt m\,g^{-1/s}E^{1/2+1/s},\quad
C_s=\frac{2\sqrt2}{\pi}\int_0^1\sqrt{1-x^s}\,\dd x .
$$

質量一定で正則な周期族をゆっくりたどる場合、$I$ を固定すれば $E\propto g^{2/(s+2)}$。
特に $s=2$ で $E\propto\sqrt g$ は振動数に比例する結果と一致する。
原点で力が滑らかでない指数では、前節の滑らかな平均化定理の仮定は別途確認を要する。
:::
::::

::::exercise{level=advanced}
初期位相が $S_0(u)=-mu^2/(2t_f)$ の一次元自由粒子の束を考える。焦点前後の主関数を求め、$t=t_f$ で Hamilton 相流が特異になるかを判定せよ。

:::solution
初期運動量は $-mu/t_f$ なので $q=u(1-t/t_f)$。$t\ne t_f$ なら $u$ を消去でき、実軌道の作用を加えた式は $S=-mq^2/[2(t_f-t)]$。微分すると $S_t= -mq^2/[2(t_f-t)^2]$、$S_q^2/(2m)=mq^2/[2(t_f-t)^2]$ が打ち消し合う。$t=t_f$ ではすべての位置が零となるが運動量は $-mu/t_f$ で異なる。全相空間の流れは常に $(q_0,p_0)\mapsto(q_0+tp_0/m,p_0)$ で可逆であり、初期位相を与えた曲線の位置軸への投影だけが特異になった。
:::
::::

## 確認テスト

::::quiz{answer="3" tol="1e-6"}
$E=6\,\mathrm J,\omega=2\,\mathrm{s^{-1}}$ の調和振動子の作用を $\mathrm{J\,s}$ で答えよ。

:::explanation
位相楕円の面積を $2\pi$ で割って $I=E/\omega=3\,\mathrm{J\,s}$ となる。
:::
::::

::::quiz
主関数 $S$ から旧運動量を得る操作はどれか。

- [x] 座標についての偏微分
- [ ] 時刻についての偏微分
- [ ] 新運動量についての偏微分

:::explanation
第二型の母関数で $p_i=S_{q^i}$ である。
$S_t$ は $-H$、$S_{\alpha_i}$ は共役な積分定数になる。
:::
::::

::::quiz{answer="2" tol="1e-6"}
$H=I^2/(2C), I=6\,\mathrm{J\,s}, C=3\,\mathrm{kg\,m^2}$ の角速度を毎秒で答えよ。

:::explanation
正準方程式より $\dot\theta=H_I=I/C=2\,\mathrm{s^{-1}}$ である。
:::
::::

::::quiz
断熱近似で特に問題となる場所はどこか。

- [ ] 周波数が一定の正則軌道
- [x] 周期が発散するセパラトリクスの近傍
- [ ] 角の原点を変えた点

:::explanation
周期が発散すると外部変化の間に多数回振動する仮定が破れ、$\dot\lambda/\omega$ を小量として扱えない。
:::
::::

::::quiz
振り子を往復側からセパラトリクスに近づけたときの作用と周期はどうなるか。

- [ ] 作用と周期がともに零に近づく
- [x] 作用は有限に近づき、周期は発散する
- [ ] 周期が有限に近づき、作用が発散する

:::explanation
囲まれる相平面の面積は有限である。
一方、上端近くで粒子がほぼ止まるため周期積分が対数発散し、$\partial I/\partial E=T/(2\pi)$ が大きくなる。
:::
::::

::::quiz{answer="1.1" tol="1e-6"}
調和振動子のばね定数をゆっくり $1.21$ 倍にする。
質量一定なら、終状態のエネルギーは初めのおよそ何倍か。

:::explanation
$I=E/\omega$ は断熱的に一定で $\omega=\sqrt{k/m}$。
従ってエネルギー比は $\sqrt{1.21}=1.1$ である。
:::
::::

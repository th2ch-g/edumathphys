---
title: 応力と保存則
description: 面力から応力を構成し、運動量・エネルギーの収支と構成式の制約を導く。
notation:
  - symbol: '\sigma'
    meaning: 引張を正とする Cauchy 応力
  - symbol: '\vb t(\vb n)'
    meaning: 外向き法線 n を持つ面に外側が及ぼす面力密度
  - symbol: '\vb b,\vb q,r'
    meaning: 単位質量当たりの体積力、熱流束、質量比加熱率
  - symbol: 'e,s,T,\psi'
    meaning: 比内部エネルギー、比エントロピー、絶対温度、比自由エネルギー
---

## この章の目標

- 各法線方向の面力を一つのテンソルで表す。
- 物質体の収支を局所的な運動方程式とエネルギー式へ変換する。
- エントロピー不等式と客観性から構成式を制約する。
- 粘弾性の緩和とクリープの違いを計算する。

## Cauchy の応力原理と応力テンソル

物体を仮想的に切ると、切断面を介して隣の部分から力が伝わる。引張られた棒では断面に垂直な力が伝わり、滑る二層の間では接線方向にも力が伝わる。本章は接触力を面積に比例する力で近似し、境界の縁に集中する力や接触偶力を除く。場は滑らかとし、添字はデカルト成分、単位は SI とする。

:::principle[Cauchy の応力原理]{#prin-cont-traction}
切断面を介する接触力は面積積分 $\int_S\vb t(\vb x,t,\vb n)\,\mathrm da$ であり、:term[面力]{en="traction" yomi="めんりょく"} 密度 $\vb t$ は位置、時刻、外向き単位法線だけで決まる。ここで位置に関して連続かつ局所有界とする。
:::

:::theorem[応力の存在]{#thm-cont-cauchy-stress}
質量密度、体積力、加速度が局所有界であり、任意の物質部分体に並進の運動量収支が成り立つとする。このとき面力は

$$
t_i(\vb n)=\sigma_{ij}n_j
\label{eq-cont-traction}
$$

と表せる。$\sigma_{ij}$ は法線が $\vb e_j$ の面に作用する力密度の第 $i$ 成分であり、:term[Cauchy 応力]{en="Cauchy stress" yomi="こーしーおうりょく"} と呼ぶ。
:::

:::proof
まず点を中心に、辺長 $\ell h$ の正方形を法線 $\vb n$ に垂直にとり、厚さ $\ell h^2$ の薄い直方体を作る。ここで $\ell$ は固定した長さ、$h$ は無次元である。二つの正方形面の面積は $\order{h^2}$、側面積は $\order{h^3}$、体積は $\order{h^4}$ である。六面の法線は $h$ によらないので、各面力の局所有界性で側面の寄与を抑えられる。並進の運動量収支を面積で割って $h\to0$ とすれば $\vb t(\vb n)+\vb t(-\vb n)=0$ を得る。

次に傾斜面の法線を $\vb n$ とし、他の三面を座標面に平行にとった相似な四面体を縮める。傾斜面の面積 $A$ に対し座標面は $A|n_j|$、外向き法線は $-\sgn(n_j)\vb e_j$ となる。これは閉曲面の向き付き面積の和がゼロであることからも従う。体積力と慣性項は $\order{h^3}$、面積は $\order{h^2}$ である。収支を $A$ で割って極限を取ると $\vb t(\vb n)-\sum_jn_j\vb t(\vb e_j)=0$ となる。$n_j=0$ の場合は柱状に切り、残りの法線が互いに逆向きの二面を加えればその寄与が相殺され、同じ面積比の議論が使える。基底での値が全成分を指定するので一意である。
:::

:::example[圧力とせん断の重ね合わせ]
$\sigma=\left(\begin{smallmatrix}-100&20&0\\20&-100&0\\0&0&-100\end{smallmatrix}\right)\,\mathrm{kPa}$ とする。$\vb n=\vb e_2$ の面には $(20,-100,0)^T\,\mathrm{kPa}$ が作用する。面積 $10^{-3}\,\mathrm{m^2}$ の平面では接線方向 $20\,\mathrm N$、法線の逆向きに $100\,\mathrm N$ である。法線を反転すると両成分とも符号が変わる。
:::

参照面積当たりの力を表す第 1 Piola–Kirchhoff 応力 $P$ は、:eqref[eq-cont-nanson] に :eqref[eq-cont-traction] を代入して $P=J\sigma F^{-T}$ となる。参照面に作用する力は $P\vb N\,\mathrm dA$ である。成分の片方は現在配置、もう片方は参照配置に属するため、通常 $P$ の転置をとっても同じ行列にはならない。

:::example[同じ力を二つの面積で割る]
$F=\diag(2,1,1/2)$、$\sigma=-pI$ とする。$\vb N=\vb e_3$ の現在面積は初めの 2 倍である。従って $P\vb N=-2p\vb e_3$ であり、参照面積当たりの圧縮力は現在面積当たりの圧縮力の 2 倍である。
:::

### 接触力の局所化と二種類の面積

Cauchy の仮定では、力密度は位置と向き付き単位法線によって決まる。辺の丸みや微小片の全体形状にも依存するモデルには、面積に比例する項のほかに縁の寄与が現れうる。応力原理を使うときは、面積に比例する局所極限を採用している。相対する二面について、平らな側面をもつ厚さ $\delta$ の薄片を使えば側面の面積は $\order{\delta}$、体積力と慣性項も $\order{\delta}$ である。幅を固定して $\delta\to0$ とし、その後に幅を縮めると、上下の力が打ち消し合うので $\vb t(-\vb n)=-\vb t(\vb n)$ を得る。四面体では全辺を同時に縮め、面積の主項を比較する。これら二つの極限で比較する力は異なる。

平面応力の例として $\sigma_{11}=s$、$\sigma_{12}=\sigma_{21}=\tau$、他をゼロとする。単位法線を $\vb n=(\cos\theta,\sin\theta,0)$、正の接線を $\vb m=(-\sin\theta,\cos\theta,0)$ とする。$\sigma\vb n$ にそれぞれを内積して $s_n=s\cos^2\theta+2\tau\sin\theta\cos\theta$、$t_m=\tau\cos2\theta-s\sin2\theta/2$ を得る。面の向きを指定する角度は、面そのものと $x$ 軸の角度とは $\pi/2$ 異なる。

:::example[斜めの断面に働く力]
$s=8 \mathrm{MPa}$、$\tau=2 \mathrm{MPa}$、$\theta=\frac{\pi}{6}\simeq 0.5236$ とする。以下では$s_{n}$ が法線応力、$t_{m}$ が符号付き接線応力、$T$ が面力密度の大きさである。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
s_{n}&=s \cos^{2}{\left(\theta \right)} + 2 \tau \sin{\left(\theta \right)} \cos{\left(\theta \right)}=\sqrt{3} + 6\simeq 7.732\quad[\mathrm{MPa}]\\
t_{m}&=- \frac{s \sin{\left(2 \theta \right)}}{2} + \tau \cos{\left(2 \theta \right)}=1 - 2 \sqrt{3}\simeq -2.4641\quad[\mathrm{MPa}]\\
T&=\sqrt{\left(- \frac{s \sin{\left(2 \theta \right)}}{2} + \tau \cos{\left(2 \theta \right)}\right)^{2} + \left(s \cos^{2}{\left(\theta \right)} + 2 \tau \sin{\left(\theta \right)} \cos{\left(\theta \right)}\right)^{2}}=\sqrt{\left(1 - 2 \sqrt{3}\right)^{2} + \left(\sqrt{3} + 6\right)^{2}}=2 \sqrt{2 \sqrt{3} + 13}\quad[\mathrm{MPa}]
\end{aligned}
$$

接線応力の負号は $\vb m$ と逆向きの力を表す。これを法線応力の符号と読み違えると、同じ面が引張か圧縮かを誤る。
:::

:::example[せん断荷重を逆転させた斜面]
$s=8 \mathrm{MPa}$、$\tau=-2 \mathrm{MPa}$、$\theta=\frac{\pi}{3}\simeq 1.0472$ とする。以下では$s_{n}$ が法線応力、$t_{m}$ が符号付き接線応力、$T$ が面力密度の大きさである。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
s_{n}&=s \cos^{2}{\left(\theta \right)} + 2 \tau \sin{\left(\theta \right)} \cos{\left(\theta \right)}=2 - \sqrt{3}\simeq 0.26795\quad[\mathrm{MPa}]\\
t_{m}&=- \frac{s \sin{\left(2 \theta \right)}}{2} + \tau \cos{\left(2 \theta \right)}=1 - 2 \sqrt{3}\simeq -2.4641\quad[\mathrm{MPa}]\\
T&=\sqrt{\left(- \frac{s \sin{\left(2 \theta \right)}}{2} + \tau \cos{\left(2 \theta \right)}\right)^{2} + \left(s \cos^{2}{\left(\theta \right)} + 2 \tau \sin{\left(\theta \right)} \cos{\left(\theta \right)}\right)^{2}}=\sqrt{\left(2 - \sqrt{3}\right)^{2} + \left(1 - 2 \sqrt{3}\right)^{2}}=2 \sqrt{5 - 2 \sqrt{3}}\quad[\mathrm{MPa}]
\end{aligned}
$$

法線応力にもせん断荷重が寄与する。$\tau$ の符号を変えると、法線の二成分が同符号であるこの面ではその寄与が引張から圧縮に変わる。
:::

:::example[法線の測定角と面力密度の誤差]
$s=8 \mathrm{MPa}$、$\tau=2 \mathrm{MPa}$ に固定し、法線角 $\theta$ を $\theta=\frac{\pi}{6}\simeq 0.5236$ の近くで変える。面力密度の大きさ $T$ の局所的な変化を求める。固定した値を先に代入すると、残る一変数の関数は次のようになる。以下の微分では数値を先の単位系で表す。

$$
\begin{aligned}
T(\theta)&=\sqrt{\left(4 \sin{\left(\theta \right)} \cos{\left(\theta \right)} + 8 \cos^{2}{\left(\theta \right)}\right)^{2} + \left(- 4 \sin{\left(2 \theta \right)} + 2 \cos{\left(2 \theta \right)}\right)^{2}},\\
\frac{\mathrm dT}{\mathrm d\theta}&=- \frac{4 \left(2 \sin^{3}{\left(\theta \right)} \cos{\left(\theta \right)} + 12 \sin^{2}{\left(\theta \right)} \cos^{2}{\left(\theta \right)} + 14 \sin{\left(\theta \right)} \cos^{3}{\left(\theta \right)} - 2 \sin^{2}{\left(2 \theta \right)} - 3 \sin{\left(2 \theta \right)} \cos{\left(2 \theta \right)} - 4 \cos^{4}{\left(\theta \right)} + 2 \cos^{2}{\left(2 \theta \right)}\right)}{\sqrt{4 \sin^{2}{\left(\theta \right)} \cos^{2}{\left(\theta \right)} + 16 \sin{\left(\theta \right)} \cos^{3}{\left(\theta \right)} + 4 \sin^{2}{\left(2 \theta \right)} - 4 \sin{\left(2 \theta \right)} \cos{\left(2 \theta \right)} + 16 \cos^{4}{\left(\theta \right)} + \cos^{2}{\left(2 \theta \right)}}}.
\end{aligned}
$$

基準値で $T=2 \sqrt{2 \sqrt{3} + 13}$、微分は $\frac{4 \left(1 - 2 \sqrt{3}\right)}{\sqrt{2 \sqrt{3} + 13}}$ である。したがって一次の展開は、基準値からの相対変化を $\theta\mapsto \theta(1+\delta)$ と書くと

$$
\frac{\Delta T}{T}
=\left(\frac{\theta}{T}
\frac{\mathrm dT}{\mathrm d\theta}\right)\delta+\order{\delta^2}
=\frac{(25-28\sqrt{3})\pi}{471}\,\delta+\order{\delta^2}
$$

となる。例えば入力を $1\%$ 増やすとき、一次式で得る出力の変化率は $\frac{\pi \left(25 - 28 \sqrt{3}\right)}{471}\simeq -0.15673\%$ である。元の式に直接代入した変化率は $-100 + \frac{100 \sqrt{4 \cos{\left(\frac{49 \pi}{300} \right)} + 8 \sin{\left(\frac{49 \pi}{300} \right)} + 9}}{\sqrt{2 \sqrt{3} + 13}}\%$ なので、この比較で一次化の誤差を点検できる。

逆に、出力の測定値が基準値から $1\%$ 上がったとき、適用した模型が正しく、ほかの入力が固定されていれば、入力の変化率は一次で $\frac{471}{\pi \left(25 - 28 \sqrt{3}\right)}\simeq -6.3804\%$ と推定される。角度はラジアンで測る。面力の大きさを角度から求める際は、法線応力と接線応力をそれぞれ求めてから大きさを合成する。
:::

## 応力テンソルの対称性と主応力

偶力を独立に伝えない物体では、角運動量の収支から $\sigma=\sigma^T$ が従う。導出は次節で与える。この対称性により :ref[thm-linalg-real-spectral] を適用できる。

:::definition[主応力]{#def-cont-principal-stress}
単位法線 $\vb n$ の面に働く法線応力を $\sigma_n=\vb n^T\sigma\vb n$、せん断面力を $(I-\vb n\vb n^T)\sigma\vb n$ とする。せん断面力が消える向きは $\sigma\vb n=\sigma_n\vb n$ を満たす。その向きを主方向、対応する固有値を :term[主応力]{en="principal stress" yomi="しゅおうりょく"} という。
:::

主方向を基底にして $\sigma_1\ge\sigma_2\ge\sigma_3$ と並べると $\sigma_n=\sum_i\sigma_i n_i^2$、$\sum_i n_i^2=1$ なので範囲は $[\sigma_3,\sigma_1]$ である。せん断面力の二乗は $\sum_i\sigma_i^2 n_i^2-\sigma_n^2$ である。各 $\sigma_i$ で $(\sigma_i-\sigma_3)(\sigma_1-\sigma_i)\ge0$ とすると

$$
|\vb t-\sigma_n\vb n|^2
\le(\sigma_1-\sigma_n)(\sigma_n-\sigma_3)
\le(\sigma_1-\sigma_3)^2/4.
$$

第 1 と第 3 主方向に等しい二乗成分を持つ法線が上限を達成する。平面問題では主方向と成す角を $\theta$ とすると $\sigma_n=(\sigma_1+\sigma_2)/2+(\sigma_1-\sigma_2)\cos 2\theta/2$、接線 $(-\sin\theta,\cos\theta)$ の面力は $-(\sigma_1-\sigma_2)\sin 2\theta/2$ である。これは Mohr の円の媒介表示である。

:::example[純せん断の主軸]
$\sigma_{12}=\sigma_{21}=30\,\mathrm{kPa}$ だけが非零なら、$(1,1,0)/\sqrt2$ に $30\,\mathrm{kPa}$、$(1,-1,0)/\sqrt2$ に $-30\,\mathrm{kPa}$、第 3 軸にゼロが対応する。最大せん断は $30\,\mathrm{kPa}$ で座標面に現れる。
:::

平均圧縮応力を $p=-\tr\sigma/3$、偏差応力を $S=\sigma+pI$ とおくと $\tr S=0$ である。$\sigma=-pI$ ならすべての面力が法線方向であり、剪断を支えない静止流体の応力に当たる。運動中の圧縮性流体では、ここで定めた $p$ と熱力学的圧力を区別する。

主軸を計算せずに応力の固有値を求めるには、特性多項式を用いる。対角化した基底で $\det(zI-\sigma)=\prod_i(z-\sigma_i)$ を展開すると、係数は $I_1=\tr\sigma$、$I_2=[(\tr\sigma)^2-\tr(\sigma^2)]/2$、$I_3=\det\sigma$ であり、多項式は $z^3-I_1z^2+I_2z-I_3$ となる。これらは直交基底の取り替えで変わらない。偏差応力については第二主不変量 $I_2(S)=-\tr(S^2)/2$ の負を $J_2=\tfrac12 S:S\ge0$ と書く。例えば純せん断で主応力が $(\tau,0,-\tau)$ なら $J_2=\tau^2$ である。一様な圧力を加えても $S$ が変わらないため、$J_2$ は変わらない。

### Mohr 円と三次元の最大せん断

主軸を用いて $\vb n=(\cos\theta,\sin\theta,0)$ を取る。$\vb t=\sigma\vb n$ の法線成分は $s_n=(s_1+s_2)/2+(s_1-s_2)\cos2\theta/2$ で、接線成分は $t_m=-(s_1-s_2)\sin2\theta/2$ である。従って $(s_n,t_m)$ は中心 $((s_1+s_2)/2,0)$、半径 $(s_1-s_2)/2$ の :term[Mohr の応力円]{en="Mohr circle for stress" yomi="もーるのおうりょくえん"} 上を動く。面を $\theta$ 回したとき円上の角度が $2\theta$ 動くのは、面力を求める法線と射影方向の両方を回すためである。

三次元では $s_1\ge s_2\ge s_3$ と順序を決め、主軸に対する法線の成分を $n_i$ とする。$w_i=n_i^2$ は正またはゼロで総和が $1$。法線応力 $s_n=\sum_i w_i s_i$、接線応力の二乗は $\sum_iw_i s_i^2-s_n^2$ である。各 $s_i$ は $(s_1-s_i)(s_i-s_3)\ge0$ を満たすため、後式は $(s_1-s_n)(s_n-s_3)$ 以下である。この二次式の最大は $(s_1-s_3)^2/4$。中間の軸の成分をゼロとし、残りに $w_1=w_3=1/2$ を取れば両方の不等号が等号になる。

:::example[三主応力から臨界面を探す]
$s_{1}=20 \mathrm{MPa}$、$s_{2}=5 \mathrm{MPa}$、$s_{3}=-4 \mathrm{MPa}$ とする。以下では$m$ が平均垂直応力、$J_{2}$ が偏差応力の第２不変量、$\tau_{M}$ が最大せん断応力である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
m&=\frac{s_{1}}{3} + \frac{s_{2}}{3} + \frac{s_{3}}{3}=7\quad[\mathrm{MPa}]\\
J_{2}&=\frac{\left(- \frac{s_{1}}{3} - \frac{s_{2}}{3} + \frac{2 s_{3}}{3}\right)^{2}}{2} + \frac{\left(- \frac{s_{1}}{3} + \frac{2 s_{2}}{3} - \frac{s_{3}}{3}\right)^{2}}{2} + \frac{\left(\frac{2 s_{1}}{3} - \frac{s_{2}}{3} - \frac{s_{3}}{3}\right)^{2}}{2}=147\quad[\mathrm{MPa^2}]\\
\tau_{M}&=\frac{s_{1}}{2} - \frac{s_{3}}{2}=12\quad[\mathrm{MPa}]
\end{aligned}
$$

最大せん断は第１軸と第３軸に等しい絶対値の方向余弦を持つ面で現れる。この面の法線応力は $8\,\mathrm{MPa}$。
:::

:::example[純せん断に対応する主応力]
$s_{1}=15 \mathrm{MPa}$、$s_{2}=0 \mathrm{MPa}$、$s_{3}=-15 \mathrm{MPa}$ とする。以下では$m$ が平均垂直応力、$J_{2}$ が偏差応力の第２不変量、$\tau_{M}$ が最大せん断応力である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
m&=\frac{s_{1}}{3} + \frac{s_{2}}{3} + \frac{s_{3}}{3}=0\quad[\mathrm{MPa}]\\
J_{2}&=\frac{\left(- \frac{s_{1}}{3} - \frac{s_{2}}{3} + \frac{2 s_{3}}{3}\right)^{2}}{2} + \frac{\left(- \frac{s_{1}}{3} + \frac{2 s_{2}}{3} - \frac{s_{3}}{3}\right)^{2}}{2} + \frac{\left(\frac{2 s_{1}}{3} - \frac{s_{2}}{3} - \frac{s_{3}}{3}\right)^{2}}{2}=225\quad[\mathrm{MPa^2}]\\
\tau_{M}&=\frac{s_{1}}{2} - \frac{s_{3}}{2}=15\quad[\mathrm{MPa}]
\end{aligned}
$$

平均垂直応力はゼロであり、せん断の最大面でも法線応力はゼロである。第２不変量は非零なので、平均垂直応力だけではせん断の大きさを判断できない。
:::

:::example[最大せん断から主方向の荷重変化を読む]
$s_{2}=5 \mathrm{MPa}$、$s_{3}=-4 \mathrm{MPa}$ に固定し、最大主応力 $s_{1}$ を $s_{1}=20 \mathrm{MPa}$ の近くで変える。最大せん断応力 $\tau_{M}$ の局所的な変化を求める。固定した値を先に代入すると、残る一変数の関数は次のようになる。以下の微分では数値を先の単位系で表す。

$$
\begin{aligned}
\tau_{M}(s_{1})&=\frac{s_{1}}{2} + 2,\\
\frac{\mathrm d\tau_{M}}{\mathrm ds_{1}}&=\frac{1}{2}.
\end{aligned}
$$

基準値で $\tau_{M}=12$、微分は $\frac{1}{2}$ である。したがって一次の展開は、基準値からの相対変化を $s_{1}\mapsto s_{1}(1+\delta)$ と書くと

$$
\frac{\Delta \tau_{M}}{\tau_{M}}
=\left(\frac{s_{1}}{\tau_{M}}
\frac{\mathrm d\tau_{M}}{\mathrm ds_{1}}\right)\delta+\order{\delta^2}
=\frac{5}{6}\,\delta+\order{\delta^2}
$$

となる。例えば入力を $1\%$ 増やすとき、一次式で得る出力の変化率は $\frac{5}{6}\%$ である。元の式に直接代入した変化率は $\frac{5}{6}\%$ なので、この比較で一次化の誤差を点検できる。

逆に、出力の測定値が基準値から $1\%$ 上がったとき、適用した模型が正しく、ほかの入力が固定されていれば、入力の変化率は一次で $\frac{6}{5}\%$ と推定される。全主応力に同じ圧力を加える操作では、最大せん断の差が消える。この推定では最小主応力を独立に固定している。
:::

## 運動量・角運動量の保存則と Cauchy の運動方程式

:::law[Cauchy の運動方程式]{#law-cont-momentum}
$\vb b$ を単位質量当たりの体積力とすると、物質部分体 $B_t$ の運動量収支の局所形は

$$
\rho D_tv_i=\partial_j\sigma_{ij}+\rho b_i,
\qquad
\partial_t(\rho v_i)+\partial_j(\rho v_i v_j-\sigma_{ij})=\rho b_i.
\label{eq-cont-momentum}
$$

:::

第 1 式は $\frac{\mathrm d}{\mathrm dt}\int_{B_t}\rho v_i\,\mathrm dv=\int_{B_t}\rho b_i\,\mathrm dv+\int_{\partial B_t}\sigma_{ij}n_j\,\mathrm da$ から従う。:ref[thm-cont-transport] と :ref[law-cont-mass] を使うと左辺は $\int\rho D_tv_i\,\mathrm dv$ となる。右辺の面積分には発散定理を使い、任意部分体に対する等式から積分の中身が等しいと分かる。第 2 式は積の微分で得られ、流出運動量の流束に応力の負号が付く。

:::proposition[角運動量と応力対称性]{#prop-cont-stress-symmetry}
接触偶力と体積偶力がなく、角運動量密度が $\rho\vb x\times\vb v$ であるとする。線運動量収支の下で、角運動量収支は応力の対称性に同値である。
:::

:::proof
角運動量の物質微分で $D_t\vb x=\vb v$ の寄与は $\vb v\times\vb v=0$ で消える。成分 $i$ の表面トルクを発散定理で変換すると

$$
\partial_l(\epsilon_{ijk}x_j\sigma_{kl})
=\epsilon_{ijk}x_j\partial_l\sigma_{kl}+\epsilon_{ijk}\sigma_{kj}.
$$

力のモーメントに相当する第 1 項と体積力は線運動量収支により左辺と同じになる。従って任意部分体で $\int\epsilon_{ijk}\sigma_{kj}\,\mathrm dv=0$ が必要十分である。連続性と任意性から $\epsilon_{ijk}\sigma_{kj}=0$、すなわち三つの反対称成分がゼロになる。
:::

:::example[自重を受ける柱]
一定断面の柱で上端を $z=H$、鉛直上向きを正とし、$\sigma_{zz}(H)=0$ とする。静止では $\partial_z\sigma_{zz}=\rho g$ なので $\sigma_{zz}=-\rho g(H-z)$。$\rho=2500\,\mathrm{kg\,m^{-3}}$、$H=2\,\mathrm m$、$g=9.8\,\mathrm{m\,s^{-2}}$ なら底で $-49\,\mathrm{kPa}$ である。
:::

### 界面を越える面力の釣合い

静止する接合界面をはさむ薄い検査体積を考える。界面から負側を $-$、正側を $+$ とし、法線 $\vb n$ を負側から正側へ向ける。体積力と側面の力が有限で界面自体に表面荷重がなければ、厚さをゼロにして残るのは二つの大面積面で受ける力だけである。この力を加えると $(\sigma^+-\sigma^-)\vb n=0$ になる。一様な引張を受ける完全接合では、材料定数が異なっても界面を横切る面力は等しい。ひずみや面力に現れない応力成分は跳べる。

質量が界面を通過するなら運動量の輸送も残る。界面法線速度を $w$、共有される質量流束を $j=\rho(v_n-w)$ とおく。保存形の運動量式で動く検査体積を取ると、相対流出は $j\vb v-\sigma\vb n$。界面上に集中質量も集中力もなければ $j[\vb v]=[\sigma]\vb n$ となる。圧力だけが働く流体では法線方向の条件が $j[v_n]=-[p]$ に簡約される。したがって流速の跳びを推定するには圧力差だけでなく質量流束も要る。

:::example[流速が半分になる不連続面]
$j=20 \mathrm{kg\,m^{-2}\,s^{-1}}$、$u_{1}=10 \mathrm{m\,s^{-1}}$、$u_{2}=5 \mathrm{m\,s^{-1}}$ とする。以下では$\rho_{1}$ が負側密度、$\rho_{2}$ が正側密度、$\Delta_{p}$ が界面を渡るときの圧力増分である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
\rho_{1}&=\frac{j}{u_{1}}=2\quad[\mathrm{kg\,m^{-3}}]\\
\rho_{2}&=\frac{j}{u_{2}}=4\quad[\mathrm{kg\,m^{-3}}]\\
\Delta_{p}&=j \left(u_{1} - u_{2}\right)=100\quad[\mathrm{Pa}]
\end{aligned}
$$

密度が２倍になるので、単位面積を単位時間に通る質量が等しい。流体の法線方向の運動量が減るため、正側の圧力が高い。
:::

:::example[高密度側への運動量流束]
$j=400 \mathrm{kg\,m^{-2}\,s^{-1}}$、$u_{1}=4 \mathrm{m\,s^{-1}}$、$u_{2}=1 \mathrm{m\,s^{-1}}$ とする。以下では$\rho_{1}$ が負側密度、$\rho_{2}$ が正側密度、$\Delta_{p}$ が界面を渡るときの圧力増分である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
\rho_{1}&=\frac{j}{u_{1}}=100\quad[\mathrm{kg\,m^{-3}}]\\
\rho_{2}&=\frac{j}{u_{2}}=400\quad[\mathrm{kg\,m^{-3}}]\\
\Delta_{p}&=j \left(u_{1} - u_{2}\right)=1200\quad[\mathrm{Pa}]
\end{aligned}
$$

両側の密度は質量収支から導いた。ここで得た数値は質量と運動量の二条件を満たすもので、実在流体の衝撃波として採用するにはエネルギーと構成式の条件も要る。
:::

:::example[下流速度から圧力差を再評価する]
$j=20 \mathrm{kg\,m^{-2}\,s^{-1}}$、$u_{1}=10 \mathrm{m\,s^{-1}}$ に固定し、正側の界面相対速度 $u_{2}$ を $u_{2}=5 \mathrm{m\,s^{-1}}$ の近くで変える。界面を渡るときの圧力増分 $\Delta_{p}$ の局所的な変化を求める。固定した値を先に代入すると、残る一変数の関数は次のようになる。以下の微分では数値を先の単位系で表す。

$$
\begin{aligned}
\Delta_{p}(u_{2})&=200 - 20 u_{2},\\
\frac{\mathrm d\Delta_{p}}{\mathrm du_{2}}&=-20.
\end{aligned}
$$

基準値で $\Delta_{p}=100$、微分は $-20$ である。したがって一次の展開は、基準値からの相対変化を $u_{2}\mapsto u_{2}(1+\delta)$ と書くと

$$
\frac{\Delta \Delta_{p}}{\Delta_{p}}
=\left(\frac{u_{2}}{\Delta_{p}}
\frac{\mathrm d\Delta_{p}}{\mathrm du_{2}}\right)\delta+\order{\delta^2}
=-1\,\delta+\order{\delta^2}
$$

となる。例えば入力を $1\%$ 増やすとき、一次式で得る出力の変化率は $-1\%$ である。元の式に直接代入した変化率は $-1\%$ なので、この比較で一次化の誤差を点検できる。

逆に、出力の測定値が基準値から $1\%$ 上がったとき、適用した模型が正しく、ほかの入力が固定されていれば、入力の変化率は一次で $-1\%$ と推定される。下流速度だけ増えると圧力増分が減る。質量流束を固定した計算なので、同時に正側の密度を調整している。
:::

## エネルギー保存則

応力のする仕事は粒子の運動エネルギーだけでなく、内部エネルギーも変える。各面における力の仕事率は $\vb t\cdot\vb v$ である。熱流束 $\vb q$ の正の法線成分は外向きの熱流出、$r$ は単位質量・単位時間当たりの加熱率と決める。

:::law[内部エネルギー収支]{#law-cont-energy}
総比エネルギー $\mathcal E=e+|\vb v|^2/2$ と内部エネルギーは

$$
\begin{aligned}
\partial_t(\rho\mathcal E)+\partial_j(\rho\mathcal E v_j-\sigma_{ij}v_i+q_j)
 &=\rho\vb b\cdot\vb v+\rho r,\\
\rho D_te&=\sigma_{ij}D_{ij}-\divg\vb q+\rho r
\end{aligned}
\label{eq-cont-energy}
$$

を満たす。
:::

第 1 式では、物質体の総エネルギーの増加を面力の仕事率、体積力の仕事率、正味の熱入力の和とした。移動領域を固定した形に直すには :ref[law-cont-mass] と発散定理を使う。第 2 式を導くには :eqref[eq-cont-momentum] に $v_i$ を掛けて

$$
\rho D_t(|\vb v|^2/2)=\partial_j(\sigma_{ij}v_i)-\sigma_{ij}\partial_j v_i+\rho\vb b\cdot\vb v
$$

を得る。これを総エネルギー式から引くと内部エネルギー式になる。応力の対称性から $\sigma:W=0$ なので $\sigma_{ij}\partial_jv_i=\sigma:D$ である。ここで二重内積 $A:B=A_{ij}B_{ij}$ を用いた。

:::example[圧縮仕事の符号]
$\sigma=-pI$、$p=10^5\,\mathrm{Pa}$、$\divg\vb v=-2\,\mathrm{s^{-1}}$ とする。応力仕事は $2\times10^5\,\mathrm{W\,m^{-3}}$ である。断熱で熱源がなければ内部エネルギーが増える。
:::

:::example[移動する熱流と符号]
板の左から $100\,\mathrm{W\,m^{-2}}$ が入り右から $60\,\mathrm{W\,m^{-2}}$ が出るなら、両側の面積が $0.5\,\mathrm{m^2}$ のとき熱入力は $20\,\mathrm W$ である。左面の法線は右面の逆なので、積分 $-\int\vb q\cdot\vb n\,\mathrm da$ はこの正値を与える。
:::

### 参照配置に戻した収支

参照配置の変位や変形勾配を未知量にする固体の計算では、力を積分する面も参照配置に戻す。$P=J\sigma F^{-T}$ と :eqref[eq-cont-nanson] により、任意の参照部分体に働く表面力は $\int_{\partial B_0}P\vb N\,\mathrm dA$。体積力を $\rho_0\vb b(\chi_t)$ とし、運動量 $\int_{B_0}\rho_0\dot\chi$ を微分する。発散定理で参照座標に局所化して

$$
\rho_0\ddot\chi_i=\partial_{X_A}P_{iA}+\rho_0 b_i(\chi_t)
$$

を得る。これと現在配置の式を比較すれば $\partial_{X_A}P_{iA}=J(\partial_j\sigma_{ij})\circ\chi_t$ である。これは両表示の力の一致を用いた発散の変換である。

現在応力が対称なら、参照体積当たりの機械的仕事は

$$
P:\dot F
=J\sigma_{ij}(F^{-1})_{Aj}\dot F_{iA}
=J\sigma:L=J\sigma:D .
$$

また $S=F^{-1}P$ を第 2 Piola–Kirchhoff 応力とすると $S=JF^{-1}\sigma F^{-T}$ は対称である。$\dot E=\operatorname{sym}(F^T\dot F)$ より $S:\dot E=P:\dot F$ となる。従って $\dot E$ に対応する仕事の係数は $S$ であり、Cauchy 応力と $E$ を直接縮約すると面積と配置の違いを取り逃がす。

:::example[一軸伸長における仕事]
$F=\diag(\lambda,1,1)$、$\sigma=\diag(\sigma_1,0,0)$ とする。$P_{11}=\sigma_1$、$S_{11}=\sigma_1/\lambda$ であり、$E_{11}=(\lambda^2-1)/2$ なので $S_{11}\dot E_{11}=\sigma_1\dot\lambda$。一方 $D_{11}=\dot\lambda/\lambda$、$J=\lambda$ だから $J\sigma:D$ も同じ値である。例えば $\lambda=2$、$\dot\lambda=1\,\mathrm{s^{-1}}$、$\sigma_1=10^5\,\mathrm{Pa}$ なら、仕事率は参照体積当たり $10^5\,\mathrm{W\,m^{-3}}$、現在体積当たり $5\times10^4\,\mathrm{W\,m^{-3}}$ である。
:::

### 支配体積で力の仕事と熱を分ける

$\sigma=-pI+\tau$、$\bar e=e+|\vb v|^2/2$ とおく。局所総エネルギー式で力の仕事を発散にすると、固定した面のエネルギー流出は

$$
\rho\bar e\,\vb v\cdot\vb n-(\sigma\vb n)\cdot\vb v+\vb q\cdot\vb n
=\rho(\bar e+p/\rho)\vb v\cdot\vb n-(\tau\vb n)\cdot\vb v+\vb q\cdot\vb n
$$

である。従って流入口と出口で粘性面力の仕事が無視できるなら、移流と圧力仕事をまとめて比エンタルピー $h=e+p/\rho$ を用いることができる。壁が静止かつ滑りなしであれば $\vb v=0$ なので、壁の面に作用する粘性力が非零でも壁からの仕事はゼロである。内部ではその場合にも粘性散逸が起こる。

単一の入口・出口を持つ定常流で、流量を $\dot m$、装置から流体への軸仕事率を $\dot W$、加熱率を $\dot Q$ とする。外力が重力だけなら $\rho\vb b\cdot\vb v=-\divg(\rho gz\vb v)$ を用いて、$\dot m(\Delta h+\Delta(|\vb v|^2/2)+g\Delta z)=\dot Q+\dot W$。入口出口の値を一様とする。この式の $\dot W>0$ は装置が流体に仕事する場合である。

動く界面の法線速度を $w$、負側から正側へ向く法線を $\vb n$ とする。界面にエネルギーの蓄積も熱源もなく、体積仕事と加熱率が両側で有界なら、厚さをゼロへ縮めた検査体積で総エネルギー収支の体積項が消える。:eqref[eq-cont-energy] の面流束から境界の移動による $\rho\mathcal E w$ を引き、両面の流出を足し合わせると

$$
\bigl[j\mathcal E-(\sigma\vb n)\cdot\vb v+\vb q\cdot\vb n\bigr]=0,
\qquad j=\rho(\vb v\cdot\vb n-w)
$$

を得る。角括弧は正側の値から負側の値を引いた跳びである。例えば静止する非透過界面では $j=0$ であり、両面が静止し界面熱源もなければ法線熱流束が連続する。流体が通過する定常界面で $\sigma=-pI$、熱流束ゼロなら $p\vb v\cdot\vb n=jp/\rho$ となり、$j\ne0$ で総比エンタルピー $e+p/\rho+|\vb v|^2/2$ が両側で等しい。

:::example[軸仕事で加速する流体]
$m_{dot}=2 \mathrm{kg\,s^{-1}}$、$W_{dot}=80 \mathrm W$、$Q_{dot}=20 \mathrm W$、$v_{1}=2 \mathrm{m\,s^{-1}}$、$v_{2}=6 \mathrm{m\,s^{-1}}$ とする。以下では$D_{k}$ が運動エネルギーの増分、$D_{h}$ が同じ高さにおける比エンタルピー増分である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
D_{k}&=- \frac{v_{1}^{2}}{2} + \frac{v_{2}^{2}}{2}=16\quad[\mathrm{J\,kg^{-1}}]\\
D_{h}&=\frac{v_{1}^{2}}{2} - \frac{v_{2}^{2}}{2} + \frac{Q_{dot} + W_{dot}}{m_{dot}}=34\quad[\mathrm{J\,kg^{-1}}]
\end{aligned}
$$

入口と出口の高さを同じとすると、入力の一部が流速を増し、残りがエンタルピーを増す。$\dot W$ を内部エネルギーに全量加えると運動エネルギーの増加を落とす。
:::

:::example[加熱と同時に減速する流体]
$m_{dot}=3 \mathrm{kg\,s^{-1}}$、$W_{dot}=0 \mathrm W$、$Q_{dot}=120 \mathrm W$、$v_{1}=3 \mathrm{m\,s^{-1}}$、$v_{2}=1 \mathrm{m\,s^{-1}}$ とする。以下では$D_{k}$ が運動エネルギーの増分、$D_{h}$ が同じ高さにおける比エンタルピー増分である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
D_{k}&=- \frac{v_{1}^{2}}{2} + \frac{v_{2}^{2}}{2}=-4\quad[\mathrm{J\,kg^{-1}}]\\
D_{h}&=\frac{v_{1}^{2}}{2} - \frac{v_{2}^{2}}{2} + \frac{Q_{dot} + W_{dot}}{m_{dot}}=44\quad[\mathrm{J\,kg^{-1}}]
\end{aligned}
$$

出口で流速が下がるため、比エンタルピーの増加は熱入力を質量流量で割った値より大きい。
:::

:::example[速度計の誤差が熱収支に伝わる場合]
$m_{dot}=2 \mathrm{kg\,s^{-1}}$、$W_{dot}=80 \mathrm W$、$Q_{dot}=20 \mathrm W$、$v_{1}=2 \mathrm{m\,s^{-1}}$ に固定し、出口速さ $v_{2}$ を $v_{2}=6 \mathrm{m\,s^{-1}}$ の近くで変える。同じ高さにおける比エンタルピー増分 $D_{h}$ の局所的な変化を求める。固定した値を先に代入すると、残る一変数の関数は次のようになる。以下の微分では数値を先の単位系で表す。

$$
\begin{aligned}
D_{h}(v_{2})&=52 - \frac{v_{2}^{2}}{2},\\
\frac{\mathrm dD_{h}}{\mathrm dv_{2}}&=- v_{2}.
\end{aligned}
$$

基準値で $D_{h}=34$、微分は $-6$ である。したがって一次の展開は、基準値からの相対変化を $v_{2}\mapsto v_{2}(1+\delta)$ と書くと

$$
\frac{\Delta D_{h}}{D_{h}}
=\left(\frac{v_{2}}{D_{h}}
\frac{\mathrm dD_{h}}{\mathrm dv_{2}}\right)\delta+\order{\delta^2}
=- \frac{18}{17}\,\delta+\order{\delta^2}
$$

となる。例えば入力を $1\%$ 増やすとき、一次式で得る出力の変化率は $- \frac{18}{17}\%$ である。元の式に直接代入した変化率は $- \frac{1809}{1700}\%$ なので、この比較で一次化の誤差を点検できる。

逆に、出力の測定値が基準値から $1\%$ 上がったとき、適用した模型が正しく、ほかの入力が固定されていれば、入力の変化率は一次で $- \frac{17}{18}\%$ と推定される。出口の速さを大きく測るほど、同じ投入仕事から推定される比エンタルピー増分は小さい。
:::

## エントロピー生成と Clausius–Duhem 不等式

単一成分の局所平衡を仮定し、比エントロピー $s$ の熱による流入を $-\vb q\cdot\vb n/T$、体積加熱に伴う流入を $\rho r/T$ とする。:ref[thm-thermo-entropy-balance] の物質部分体への適用により

$$
\dot s_{\mathrm i}=\rho D_ts+\divg(\vb q/T)-\rho r/T\ge0
$$

となる。$\dot s_{\mathrm i}$ の単位は $\mathrm{W\,m^{-3}\,K^{-1}}$ である。比 Helmholtz 自由エネルギーを $\psi=e-Ts$ とおく。

:::proposition[Clausius–Duhem 不等式]{#prop-cont-dissipation}
温度 $T>0$ のとき、エントロピー生成の非負性は

$$
\sigma:D-\rho(D_t\psi+sD_tT)-\frac{\vb q\cdot\grad T}{T}\ge0
\label{eq-cont-dissipation}
$$

と同値である。
:::

:::proof
エントロピー式に $T$ を掛け、$T\divg(\vb q/T)=\divg\vb q-\vb q\cdot\grad T/T$ と展開する。:eqref[eq-cont-energy] から $\divg\vb q-\rho r=\sigma:D-\rho D_te$ を代入する。$D_te-TD_ts=D_t\psi+sD_tT$ なので左辺は $T\dot s_{\mathrm i}$ となり、$T>0$ から同値性が従う。
:::

Fourier 則 $\vb q=-k\grad T$ は $k\ge0$ のとき熱伝導による寄与を $k|\grad T|^2/T$ とする。$k=0.6\,\mathrm{W\,m^{-1}\,K^{-1}}$、$T=300\,\mathrm K$、$|\grad T|=10\,\mathrm{K\,m^{-1}}$ ならこの寄与は $0.2\,\mathrm{W\,m^{-3}}$、エントロピー生成は $0.2/300\,\mathrm{W\,m^{-3}\,K^{-1}}$ である。

単純圧縮物質で $\psi=\psi(\rho,T)$ とする。$D_t\rho=-\rho\tr D$ を使うと

$$
T\dot s_{\mathrm i}=(\sigma+\rho^2\psi_\rho I):D-\rho(\psi_T+s)D_tT-\vb q\cdot\grad T/T.
$$

温度の変化率を独立に任意の符号で選べる構成理論では $s=-\psi_T$ が必要である。可逆応力を $-pI$ として $p=\rho^2\psi_\rho$ を選べば可逆な圧縮仕事は散逸から消える。残りを粘性応力 $\tau$ とすれば機械的散逸は $\tau:D$ である。この分離により、圧縮の際の可逆な温度上昇と摩擦による発熱を区別できる。

### 熱流を含む全体のエントロピー

Fourier 則を仮定して $\vb q=-k\grad T$ とすれば、熱伝導による体積エントロピー生成は $k|\grad T|^2/T^2$ である。例えば厚さ $L$ の静止板で一様な熱伝導率 $k$、内部熱源ゼロ、両面 $T_c<T_h$ を考える。熱流は高温面から低温面へ流れ、流束の大きさは $q=k(T_h-T_c)/L$。局所生成を厚さに沿って積分すると

$$
\int_0^L k \frac{(T')^2}{T^2}\,\dd x
=-k T'\left[\frac{1}{T}\right]_0^L
=q(1/T_c-1/T_h).
$$

これは同じ熱流が高温槽から抜けて低温槽に入るときの両熱槽のエントロピー変化率の和でもある。板内のエントロピーの蓄積率は定常なのでゼロである。生成が正でも、熱とともに境界から運び出される。

圧縮性 Newton 流体では粘性応力を $\tau=2\mu D'+\zeta(\divg\vb v)I$、$D'=D-(\divg\vb v)I/3$ と分ける。$\tr D'=0$ なので $\tau:D=2\mu D':D'+\zeta(\divg\vb v)^2$。温度勾配をゼロにし、純粋なせん断と純粋な膨張を別々に与えれば、任意の勾配に対する非負性から $\mu,\zeta\ge0$ が必要である。この条件と $k\ge0$ があれば、和は任意の勾配について非負である。体積変化が小さい試験だけでは $\zeta$ の値を絞れない。

:::example[温度差が 30 K の定常熱伝導]
$k=2 \mathrm{W\,m^{-1}\,K^{-1}}$、$T_{c}=300 \mathrm K$、$T_{h}=330 \mathrm K$、$L=\frac{1}{10} \mathrm m$ とする。以下では$q$ が熱流束、$s_{i}$ が単位面積当たり生成率である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
q&=\frac{k \left(- T_{c} + T_{h}\right)}{L}=600\quad[\mathrm{W\,m^{-2}}]\\
s_{i}&=\frac{k \left(- T_{c} + T_{h}\right) \left(- \frac{1}{T_{h}} + \frac{1}{T_{c}}\right)}{L}=\frac{2}{11}\quad[\mathrm{W\,m^{-2}\,K^{-1}}]
\end{aligned}
$$

正の生成は、熱力学的に低温側へ熱が流れる向きと一致する。単位面積当たりの生成を求めるには、厚さ方向に局所密度を積分する。
:::

:::example[低温域で同じ熱流量を通す]
$k=2 \mathrm{W\,m^{-1}\,K^{-1}}$、$T_{c}=200 \mathrm K$、$T_{h}=230 \mathrm K$、$L=\frac{1}{10} \mathrm m$ とする。以下では$q$ が熱流束、$s_{i}$ が単位面積当たり生成率である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
q&=\frac{k \left(- T_{c} + T_{h}\right)}{L}=600\quad[\mathrm{W\,m^{-2}}]\\
s_{i}&=\frac{k \left(- T_{c} + T_{h}\right) \left(- \frac{1}{T_{h}} + \frac{1}{T_{c}}\right)}{L}=\frac{9}{23}\quad[\mathrm{W\,m^{-2}\,K^{-1}}]
\end{aligned}
$$

熱流束は同じである。低い絶対温度では逆温度の差が大きいため、単位面積当たりの生成率が増す。
:::

:::example[高温面の温度の変化と生成率]
$k=2 \mathrm{W\,m^{-1}\,K^{-1}}$、$T_{c}=300 \mathrm K$、$L=\frac{1}{10} \mathrm m$ に固定し、高温面の温度 $T_{h}$ を $T_{h}=330 \mathrm K$ の近くで変える。単位面積当たり生成率 $s_{i}$ の局所的な変化を求める。固定した値を先に代入すると、残る一変数の関数は次のようになる。以下の微分では数値を先の単位系で表す。

$$
\begin{aligned}
s_{i}(T_{h})&=20 \left(\frac{1}{300} - \frac{1}{T_{h}}\right) \left(T_{h} - 300\right),\\
\frac{\mathrm ds_{i}}{\mathrm dT_{h}}&=\frac{\left(T_{h} - 300\right) \left(T_{h} + 300\right)}{15 T_{h}^{2}}.
\end{aligned}
$$

基準値で $s_{i}=\frac{2}{11}$、微分は $\frac{7}{605}$ である。したがって一次の展開は、基準値からの相対変化を $T_{h}\mapsto T_{h}(1+\delta)$ と書くと

$$
\frac{\Delta s_{i}}{s_{i}}
=\left(\frac{T_{h}}{s_{i}}
\frac{\mathrm ds_{i}}{\mathrm dT_{h}}\right)\delta+\order{\delta^2}
=21\,\delta+\order{\delta^2}
$$

となる。例えば入力を $1\%$ 増やすとき、一次式で得る出力の変化率は $21\%$ である。元の式に直接代入した変化率は $\frac{2221}{101}\%$ なので、この比較で一次化の誤差を点検できる。

逆に、出力の測定値が基準値から $1\%$ 上がったとき、適用した模型が正しく、ほかの入力が固定されていれば、入力の変化率は一次で $\frac{1}{21}\%$ と推定される。低温面を固定して高温面を上げると、熱流束と逆温度差が同時に増す。小温度差では生成率は温度差の二乗から始まる。
:::

### 非負の機械散逸を圧縮とせん断で別に調べる

偏差成分がゼロの変形 $D=dI$ では $\tau:D=9\zeta d^2$。したがって粘性による仕事を測るには小振幅の膨張・収縮が使える。一方、単純せん断 $v_x=\gamma y$ では $\tr D=0$、二つの対称成分がともに $\gamma/2$ なので $\tau:D=\mu\gamma^2$。同じ運動の中に両成分があれば、交差項は $D':I=0$ のため生じない。

例えば $\mu=0.2\,\mathrm{Pa\,s}$、$\zeta=0.4\,\mathrm{Pa\,s}$ とし、同じ点で $d=-1\,\mathrm{s^{-1}}$ と $\gamma=3\,\mathrm{s^{-1}}$ を重ねる。このとき体積仕事から得る粘性散逸は $3.6\,\mathrm{W\,m^{-3}}$、せん断からは $1.8\,\mathrm{W\,m^{-3}}$ である。圧縮の符号を反転しても散逸の和は同じだが、圧力による仕事 $-p\tr D$ は符号が変わる。この違いを使えば、内部エネルギーの増減と不可逆な熱生成とを区別できる。

## 構成方程式の一般原理（客観性）

力と熱の収支をそろえても、$\sigma,\vb q,e$ が変形や温度にどう応答するかは決まらない。材料の性質を指定する関係を :term[構成方程式]{en="constitutive equation" yomi="こうせいほうていしき"} という。観測者の剛体運動を付加しても材料の応答を同じ法則で記述できる性質が :term[客観性]{en="objectivity" yomi="きゃっかんせい"} である。

:::proposition[重ね合わせた剛体運動]{#prop-cont-objectivity}
$\vb x^*=Q(t)\vb x+\vb c(t)$、$Q^TQ=I$、$\det Q=1$ とする。応力と熱流は $\sigma^*=Q\sigma Q^T$、$\vb q^*=Q\vb q$ と変換する。運動学量は

$$
F^*=QF,\qquad L^*=\dot Q Q^T+QLQ^T,\qquad D^*=QDQ^T
$$

を満たす。
:::

:::proof
力と法線はともに $Q$ を掛けたベクトルに移り、面積は保たれる。従って $\sigma^*Q\vb n=Q\sigma\vb n$ がすべての法線で成り立つ。熱流も各法線への投影を保つためベクトルとして変換する。$\vb x^*$ の参照微分は $F^*=QF$、時刻微分は $\vb v^*=\dot Q\vb x+Q\vb v+\dot{\vb c}$ である。$\partial\vb x/\partial\vb x^*=Q^T$ なので速度勾配の式を得る。さらに $QQ^T=I$ の微分により $\dot Q Q^T$ は反対称であり、対称部分をとれば $D^*$ の式になる。
:::

例えば同一の物体に重ね合わせた剛体回転では、$L$ の反対称部分だけが非零であり $D=0$ となる。従って $L$ に比例する応力は剛体回転の寄与を拾う一方、$D$ に比例する応力は所定のテンソル変換を満たす。客観性は観測者の選び方に関する。結晶の材料軸を物体に対して変えると応答が変わる異方性は、材料軸も一緒に変換すれば客観性と両立する。

テンソルの普通の物質微分は、時間依存する $Q$ の微分項を含む。$A^*=QAQ^T$ を微分し、$W^*=QWQ^T+\dot Q Q^T$ を代入すると、$D_tA-WA+AW$ の変換では余分な二項が相殺される。これが Jaumann 率である。どの客観率を使うべきかは構成模型の仮定で決まる。

### 構成式の客観性を有限弾性で使う

:term[超弾性体]{en="hyperelastic material" yomi="ちょうだんせいたい"} は変形を準静的に変えたとき、参照体積当たりの応力仕事をエネルギー関数 $\mathcal W(F)$ の変化で表せる材料である。等温で内部散逸がなく、任意の変形速度を許すとすれば $P:\dot F=\mathrm d\mathcal W/\mathrm dt$ であり、各成分の係数を比較して $P_{iA}=\partial\mathcal W/\partial F_{iA}$ となる。客観性のため $\mathcal W(QF)=\mathcal W(F)$ が必要である。:ref[thm-cont-polar] の $F=RU$ に $Q=R^T$ を用いれば $\mathcal W(F)=\mathcal W(U)$。$U$ は $C=F^TF$ によって一意に定まるので、客観的エネルギーを $C$ の関数として書ける。逆に $C(QF)=C(F)$ なのでそのようなエネルギーは客観的である。

例えば非圧縮の等方固体の簡単な模型として、制約 $J=1$ の下で $\mathcal W=\mu(\tr C-3)/2$ を選ぶ。$J$ の一次変化は $\delta J=JF^{-T}:\delta F$ であるから、圧力乗数 $p$ を導入すると $P=\mu F-pF^{-T}$。$P=J\sigma F^{-T}$ を戻して $\sigma=\mu FF^T-pI$ となる。$p$ は非圧縮制約のために生じ、境界条件と運動式から決まる。これは Neo-Hooke 模型の一つであり、$\mu$ の単位は $\mathrm{Pa}$ である。

:::example[ゴムの一軸伸長]
変形を $F=\diag(\lambda,\lambda^{-1/2},\lambda^{-1/2})$ とすると $J=1$。横面自由なら $p=\mu\lambda^{-1}$ となり、現在面積当たりの軸応力は $\mu(\lambda^2-\lambda^{-1})$。参照面積当たりはこれを $\lambda$ で割った $\mu(\lambda-\lambda^{-2})$ である。$\lambda=2$ ではそれぞれ $3.5\mu,1.75\mu$。$\lambda=1+\delta$ の一次ではどちらも $3\mu\delta$ となり、微小変形での一軸弾性率 $3\mu$ に一致する。
:::

### 観測者の回転と応力率

$\sigma^*=Q\sigma Q^T$、$L^*=\dot Q Q^T+QLQ^T$ とする。同じ粒子に沿って応力を微分すると $D_t^*\sigma^*=Q(D_t\sigma)Q^T+H\sigma^*-\sigma^*H$、$H=\dot Q Q^T$ である。したがって単なる物質微分は一般に応力と同じ変換則を満たさない。一方、$W^*=H+QWQ^T$ なので

$$
\overset{\circ}{\sigma}=D_t\sigma-W\sigma+\sigma W
$$

とおけば、$H$ を含む項が相殺し $\overset{\circ}{\sigma}{}^*=Q\overset{\circ}{\sigma}Q^T$ となる。これを :term[Jaumann 応力率]{en="Jaumann stress rate" yomi="やうまんおうりょくりつ"} と呼ぶ。$D_t\sigma-L\sigma-\sigma L^T$ も $H^T=-H$ によって同じ変換則を満たす。両者の差は $D\sigma+\sigma D$ である。客観性は応力率の候補を制約するが、その中から材料に合う式を一つに選ぶにはさらに物理的な条件が必要である。

非圧縮の neo-Hooke 体では $\sigma=-p I+\mu B$ とする。一軸伸長比 $\lambda$ のとき、断面が横二方向に同じ比で縮むなら $B=\diag(\lambda^2,\lambda^{-1},\lambda^{-1})$。側面の面力ゼロを課すと $p=\mu\lambda^{-1}$、軸の真応力は $\mu(\lambda^2-\lambda^{-1})$。参照断面積あたりの軸力は、現在断面積が参照の $1/\lambda$ 倍なので $\mu(\lambda-\lambda^{-2})$ となる。

:::example[軸力を現在の断面積で割る]
$\mu=1 \mathrm{MPa}$、$\lambda=\frac{3}{2}$ とする。以下では$p$ が側面条件が決める圧力、$sigma_{11}$ が軸方向の真応力、$P_{11}$ が参照面積当たりの軸力である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
p&=\frac{\mu}{\lambda}=\frac{2}{3}\quad[\mathrm{MPa}]\\
\sigma_{11}&=\mu \left(\lambda^{2} - \frac{1}{\lambda}\right)=\frac{19}{12}\quad[\mathrm{MPa}]\\
P_{11}&=\mu \left(\lambda - \frac{1}{\lambda^{2}}\right)=\frac{19}{18}\quad[\mathrm{MPa}]
\end{aligned}
$$

伸長によって断面が小さくなるので真応力は参照面積あたりの値より大きい。同じ引張荷重がどちらで報告されているかを確かめる必要がある。
:::

:::example[均一な圧縮を仮定した計算]
$\mu=2 \mathrm{MPa}$、$\lambda=\frac{4}{5}$ とする。以下では$p$ が側面条件が決める圧力、$sigma_{11}$ が軸方向の真応力、$P_{11}$ が参照面積当たりの軸力である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
p&=\frac{\mu}{\lambda}=\frac{5}{2}\quad[\mathrm{MPa}]\\
\sigma_{11}&=\mu \left(\lambda^{2} - \frac{1}{\lambda}\right)=- \frac{61}{50}\quad[\mathrm{MPa}]\\
P_{11}&=\mu \left(\lambda - \frac{1}{\lambda^{2}}\right)=- \frac{61}{40}\quad[\mathrm{MPa}]
\end{aligned}
$$

軸応力は両面積の規約で負になる。この計算は一様変形を与えたもので、細長い試験片が座屈して一様変形を失う場合は使えない。
:::

:::example[参照応力から軸伸長比を読む]
$\mu=1 \mathrm{MPa}$ に固定し、軸伸長比 $\lambda$ を $\lambda=\frac{3}{2}$ の近くで変える。参照面積当たりの軸力 $P_{11}$ の局所的な変化を求める。固定した値を先に代入すると、残る一変数の関数は次のようになる。以下の微分では数値を先の単位系で表す。

$$
\begin{aligned}
P_{11}(\lambda)&=\lambda - \frac{1}{\lambda^{2}},\\
\frac{\mathrm dP_{11}}{\mathrm d\lambda}&=\frac{\lambda^{3} + 2}{\lambda^{3}}.
\end{aligned}
$$

基準値で $P_{11}=\frac{19}{18}$、微分は $\frac{43}{27}$ である。したがって一次の展開は、基準値からの相対変化を $\lambda\mapsto \lambda(1+\delta)$ と書くと

$$
\frac{\Delta P_{11}}{P_{11}}
=\left(\frac{\lambda}{P_{11}}
\frac{\mathrm dP_{11}}{\mathrm d\lambda}\right)\delta+\order{\delta^2}
=\frac{43}{19}\,\delta+\order{\delta^2}
$$

となる。例えば入力を $1\%$ 増やすとき、一次式で得る出力の変化率は $\frac{43}{19}\%$ である。元の式に直接代入した変化率は $\frac{436227}{193819}\%$ なので、この比較で一次化の誤差を点検できる。

逆に、出力の測定値が基準値から $1\%$ 上がったとき、適用した模型が正しく、ほかの入力が固定されていれば、入力の変化率は一次で $\frac{19}{43}\%$ と推定される。$\lambda=1$ が無応力状態で、近くでの接線係数は $3\mu$ である。有限伸長では接線係数そのものが変わる。
:::

### 参照面に戻すと運動の回転から独立になる応力

第１ Piola–Kirchhoff 応力は $P=J\sigma F^{-T}$ である。参照配置で両方の添字を測る :term[第２ Piola–Kirchhoff 応力]{en="second Piola–Kirchhoff stress" yomi="だいにぴおらきるひほっふおうりょく"} を $S=F^{-1}P$ と定義する。$\sigma$ が対称なら $S=JF^{-1}\sigma F^{-T}$ も対称である。空間で剛体回転 $Q$ を重ねた場合、$F^*=QF$、$P^*=QP$ だから $S^*=(F^*)^{-1}P^*=F^{-1}P=S$。ここで参照座標は固定している。また $\dot E=(F^T\dot F+\dot F^TF)/2$ に代入すると、$S$ の対称性により

$$
S:\dot E=S:(F^T\dot F)=(FS):\dot F=P:\dot F
$$

を得る。よってこの応力は参照体積当たりの仕事を、有限ひずみ $E$ の変化で表すときに用いる。

:::example[一軸試験で三種類の応力を換算する]
$F=\diag(2,1/\sqrt2,1/\sqrt2)$、$\sigma=\diag(4,0,0)\,\mathrm{MPa}$ なら $J=1$、$F^{-1}$ の軸方向成分は $1/2$ である。右から一度掛けた $P$ の軸成分は $2\,\mathrm{MPa}$、左からもう一度掛けた $S$ の軸成分は $1\,\mathrm{MPa}$。伸長比 $\lambda$ を単位時間当たり $\dot\lambda$ だけ変えると、この時点では $\dot E_{11}=2\dot\lambda$ であり、どちらの参照応力でも参照体積当たりの仕事率は $2\dot\lambda\,\mathrm{MPa}$ となる。
:::

## 粘弾性体のモデル（Maxwell 模型・Kelvin–Voigt 模型）

線形の一軸模型を考え、応力 $\sigma$ は $\mathrm{Pa}$、ひずみ $\eps$ は無次元である。ばねの弾性率 $E>0$、ダッシュポットの粘性率 $\eta>0$ はそれぞれ $\mathrm{Pa}$、$\mathrm{Pa\,s}$ とする。

:::definition[直列と並列の粘弾性模型]{#def-cont-viscoelastic-models}
ばねとダッシュポットを直列にしたものを :term[Maxwell 模型]{en="Maxwell model" yomi="まくすうぇるもけい"}、並列にしたものを :term[Kelvin–Voigt 模型]{en="Kelvin–Voigt model" yomi="けるびんふぉいくともけい"} という。
:::

直列では応力を共通として全ひずみが $\eps_{\mathrm s}+\eps_{\mathrm d}$ である。$\eps_{\mathrm s}=\sigma/E$、$\dot\eps_{\mathrm d}=\sigma/\eta$ より $\dot\sigma+\sigma/\tau=E\dot\eps$、$\tau=\eta/E$ を得る。並列ではひずみを共通として全応力が $\sigma_{\mathrm s}+\sigma_{\mathrm d}=E\eps+\eta\dot\eps$ となる。直列の応力仕事は $\sigma\dot\eps=\mathrm d(\sigma^2/2E)/\mathrm dt+\sigma^2/\eta$、並列では $\sigma\dot\eps=\mathrm d(E\eps^2/2)/\mathrm dt+\eta\dot\eps^2$ である。どちらも弾性貯蔵と非負の散逸に分かれる。

:::example[応力緩和とクリープ]
$E=2\,\mathrm{MPa}$、$\eta=6\,\mathrm{MPa\,s}$ なら $\tau=3\,\mathrm s$ である。直列模型に瞬時に $\eps_0=0.01$ を与え、その後一定にすると初期にばねが全変形を引き受け、$\sigma=E\eps_0e^{-t/\tau}$ となる。$t=3\,\mathrm s$ で $20/e\,\mathrm{kPa}$ である。並列模型に同じ値 $20\,\mathrm{kPa}$ の応力を保持すると $\eta\dot\eps+E\eps=\sigma_0$ より $\eps=(\sigma_0/E)(1-e^{-t/\tau})$ である。零から上昇して $t=3\,\mathrm s$ に $0.01(1-e^{-1})$ に達する。
:::

逆に直列模型に一定応力を印加すると直後に $\sigma_0/E$ 伸び、続いて $\sigma_0 t/\eta$ だけ伸びる。並列模型に瞬時のひずみステップを課す理想化では $\dot\eps$ にデルタ関数が現れるので、初期応力の衝撃を伴う。その後の応力は $E\eps_0$ のままである。有限時間で荷重を立ち上げれば、衝撃は大きいが有限な応力へ置き換わる。

調和ひずみを $\Re(\eps_0 e^{i\omega t})$ で表す。直列の複素弾性率は $E^*=E i\omega\tau/(1+i\omega\tau)$、並列では $E^*=E+i\eta\omega$ である。$E^*=E'+iE''$ に対し周期平均 $\langle\sigma\dot\eps\rangle=\omega E''|\eps_0|^2/2$ は、時間依存が $e^{i\omega t}$ であるとき $E''\ge0$ を要求する。直列の低周波極限は粘性支配、高周波極限は弾性支配になる。

### 緩和時間が複数ある材料

単一の Maxwell 模型だけでは、時刻の異なる複数の緩和を表せない。平衡応力を持つばねと複数の Maxwell 枝を並列に結ぶと、同じひずみを各枝へ与えるので全応力は

$$
\sigma=E_\infty\eps+\sum_{i=1}^N \sigma_i,\qquad
\dot\sigma_i+\sigma_i/\tau_i=E_i\dot\eps,\quad
\tau_i=\eta_i/E_i
$$

となる。材料を $t\to-\infty$ で十分ゆっくり準備し、初期応力の寄与が消える場合、積分因子で $\sigma_i(t)=\int_{-\infty}^t E_i e^{-(t-s)/\tau_i}\dot\eps(s)\,\mathrm ds$ を得る。ステップひずみ $\eps_0$ への応答は $\eps_0[E_\infty+\sum_i E_i e^{-t/\tau_i}]$ である。短い時間では全ばねが応力を負担し、長い時間では平衡ばねだけが残る。

任意の滑らかな履歴で仕事を求めるには枝ごとの仕事を足す。内部エネルギーの機械的貯蔵部分を $W=E_\infty\eps^2/2+\sum_i\sigma_i^2/(2E_i)$ とすれば、構成式から $\sigma\dot\eps-\dot W=\sum_i\sigma_i^2/\eta_i\ge0$ である。従って各枝の正の弾性率と粘性率により全体も受動的である。例えば $E_\infty=1\,\mathrm{MPa}$、$E_1=2\,\mathrm{MPa}$、$E_2=3\,\mathrm{MPa}$ で $\tau_1=1\,\mathrm s$、$\tau_2=100\,\mathrm s$ なら、$\eps_0=0.01$ のステップ直後は $60\,\mathrm{kPa}$、$t=10\,\mathrm s$ では約 $37.15\,\mathrm{kPa}$、十分後は $10\,\mathrm{kPa}$ である。

### 周期負荷で蓄積と散逸を分ける

微小ひずみ $\eps=\Re(\eps_0e^{i\omega t})$ に対し、複素応力振幅を $E^*(\omega)\eps_0$ とする。Kelvin–Voigt 則に代入すると $E^*=E+i\omega\eta$。Maxwell 則 $\dot\eps=\dot\sigma/E+\sigma/\eta$ では $E^*=E\,i\omega t_r/(1+i\omega t_r)$、$t_r=\eta/E$ である。後式に分母の複素共役を掛けると

$$
E'=\frac{E(\omega t_r)^2}{1+(\omega t_r)^2},\qquad
E''=\frac{E\omega t_r}{1+(\omega t_r)^2}.
$$

$E'$ と $E''$ をそれぞれ :term[貯蔵弾性率]{en="storage modulus" yomi="ちょぞうだんせいりつ"}、:term[損失弾性率]{en="loss modulus" yomi="そんしつだんせいりつ"} と呼ぶ。正の実数 $\eps_0$ では $\sigma=\eps_0(E'\cos\omega t-E''\sin\omega t)$。一周期の仕事は $\int_0^{2\pi/\omega}\sigma\dot\eps\,\dd t$ であり、$\sin\omega t\cos\omega t$ の積分がゼロ、$\sin^2\omega t$ の平均が $1/2$ であることから $\pi E''\eps_0^2$ となる。単位は体積あたりのエネルギーである。Maxwell 模型では $r=\omega t_r$ として $\mathrm d[r/(1+r^2)]/\mathrm dr=(1-r^2)/(1+r^2)^2$。損失弾性率だけを測っても、一般には $r<1$ と $r>1$ の枝を区別できない。

:::example[緩和より遅い周期負荷]
$E=1000000 \mathrm{Pa}$、$t_{r}=\frac{1}{2} \mathrm s$、$\omega=1 \mathrm{s^{-1}}$、$\epsilon_{0}=\frac{1}{100}$ とする。以下では$r$ が無次元角周波数、$E_{s}$ が貯蔵弾性率、$E_{l}$ が損失弾性率、$W$ が一周期の散逸である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
r&=\omega t_{r}=\frac{1}{2}\\
E_{s}&=\frac{E \omega^{2} t_{r}^{2}}{\omega^{2} t_{r}^{2} + 1}=200000\quad[\mathrm{Pa}]\\
E_{l}&=\frac{E \omega t_{r}}{\omega^{2} t_{r}^{2} + 1}=400000\quad[\mathrm{Pa}]\\
W&=\frac{\pi E \epsilon_{0}^{2} \omega t_{r}}{\omega^{2} t_{r}^{2} + 1}=40 \pi\simeq 125.66\quad[\mathrm{J\,m^{-3}}]
\end{aligned}
$$

応力の位相ずれには損失の寄与が大きい。$E_l/E_s=1/r$ であるため、一周期の仕事と貯蔵応力の両方からこの領域を判定できる。
:::

:::example[緩和より速い周期負荷]
$E=1000000 \mathrm{Pa}$、$t_{r}=\frac{1}{2} \mathrm s$、$\omega=4 \mathrm{s^{-1}}$、$\epsilon_{0}=\frac{1}{100}$ とする。以下では$r$ が無次元角周波数、$E_{s}$ が貯蔵弾性率、$E_{l}$ が損失弾性率、$W$ が一周期の散逸である。右端の角括弧は式に使った単位を表す。

$$
\begin{aligned}
r&=\omega t_{r}=2\\
E_{s}&=\frac{E \omega^{2} t_{r}^{2}}{\omega^{2} t_{r}^{2} + 1}=800000\quad[\mathrm{Pa}]\\
E_{l}&=\frac{E \omega t_{r}}{\omega^{2} t_{r}^{2} + 1}=400000\quad[\mathrm{Pa}]\\
W&=\frac{\pi E \epsilon_{0}^{2} \omega t_{r}}{\omega^{2} t_{r}^{2} + 1}=40 \pi\simeq 125.66\quad[\mathrm{J\,m^{-3}}]
\end{aligned}
$$

$r$ が先の例の逆数になり損失弾性率は同じである。貯蔵弾性率は増えるので、応力の振幅や位相を併せて測ると区別できる。
:::

:::example[周波数を調整して散逸を測る]
$E=1000000 \mathrm{Pa}$、$t_{r}=\frac{1}{2} \mathrm s$、$\epsilon_{0}=\frac{1}{100}$ に固定し、角周波数 $\omega$ を $\omega=1 \mathrm{s^{-1}}$ の近くで変える。一周期の散逸 $W$ の局所的な変化を求める。固定した値を先に代入すると、残る一変数の関数は次のようになる。以下の微分では数値を先の単位系で表す。

$$
\begin{aligned}
W(\omega)&=\frac{50 \pi \omega}{\frac{\omega^{2}}{4} + 1},\\
\frac{\mathrm dW}{\mathrm d\omega}&=- \frac{200 \pi \left(\omega - 2\right) \left(\omega + 2\right)}{\left(\omega^{2} + 4\right)^{2}}.
\end{aligned}
$$

基準値で $W=40 \pi\simeq 125.66$、微分は $24 \pi\simeq 75.398$ である。したがって一次の展開は、基準値からの相対変化を $\omega\mapsto \omega(1+\delta)$ と書くと

$$
\frac{\Delta W}{W}
=\left(\frac{\omega}{W}
\frac{\mathrm dW}{\mathrm d\omega}\right)\delta+\order{\delta^2}
=\frac{3}{5}\,\delta+\order{\delta^2}
$$

となる。例えば入力を $1\%$ 増やすとき、一次式で得る出力の変化率は $\frac{3}{5}\%$ である。元の式に直接代入した変化率は $\frac{29900}{50201}\%$ なので、この比較で一次化の誤差を点検できる。

逆に、出力の測定値が基準値から $1\%$ 上がったとき、適用した模型が正しく、ほかの入力が固定されていれば、入力の変化率は一次で $\frac{5}{3}\%$ と推定される。ピークより低い周波数を使った予測である。高い周波数では感度が負になるため、同じ推定を使う前に枝を確認する。
:::

## 演習問題

::::exercise{level=basic}
$\sigma=\diag(10,30,-20)\,\mathrm{kPa}$、$\vb n=(1,1,0)/\sqrt2$ とする。法線応力とせん断面力の大きさを求めよ。

:::solution
$\vb t=(10,30,0)/\sqrt2\,\mathrm{kPa}$ なので $\sigma_n=(10+30)/2=20\,\mathrm{kPa}$。差は $(-10,10,0)/\sqrt2\,\mathrm{kPa}$ で大きさは $10\,\mathrm{kPa}$ である。
:::
::::

::::exercise{level=basic}
$\sigma_1=80\,\mathrm{MPa}$、$\sigma_2=20\,\mathrm{MPa}$、$\sigma_3=-40\,\mathrm{MPa}$ の最大せん断応力を求めよ。

:::solution
最大と最小の主応力差の半分なので $(80+40)/2=60\,\mathrm{MPa}$。第 1 と第 3 主方向の二等分法線で実現する。
:::
::::

::::exercise{level=basic}
$\sigma=-pI$、$p=p_0+1000x$ を SI 数値で与え、密度は $1000\,\mathrm{kg\,m^{-3}}$、体積力はないとする。加速度を求めよ。

:::solution
$\partial_j\sigma_{ij}=-\partial_i p$ なので :eqref[eq-cont-momentum] より $D_t\vb v=(-1,0,0)\,\mathrm{m\,s^{-2}}$ である。
:::
::::

::::exercise{level=basic}
$E=4\,\mathrm{MPa}$、$\eta=8\,\mathrm{MPa\,s}$ の Maxwell 模型に $0.02$ のひずみを瞬時に与え、保持する。$t=2\,\mathrm s$ の応力を求めよ。

:::solution
不変ひずみでは $\dot\sigma=-\sigma/\tau$、$\tau=2\,\mathrm s$ である。直後の応力 $80\,\mathrm{kPa}$ から指数緩和し、$t=2\,\mathrm s$ で $80/e\,\mathrm{kPa}$ になる。
:::
::::

::::exercise{level=advanced}
熱流束が $\vb q=-K\grad T$、$K$ が実行列であるとする。熱伝導の寄与が任意の温度勾配で非負となる必要十分条件を導け。

:::solution
$\vb a=\grad T$ と置くと寄与は $\vb a^TK\vb a/T$ である。反対称部分 $A$ では $\vb a^TA\vb a=-\vb a^TA\vb a$ で消える。従って必要十分条件は対称部分 $(K+K^T)/2$ が半正定値であることである。
:::
::::

::::exercise{level=advanced}
$s=s_0+c\log(T/T_0)$、$e=cT$ である熱容量一定の模型について、$\psi_T=-s$ を確かめよ。

:::solution
$\psi=cT-Ts$ を微分すると $c-s-T(c/T)=-s$ である。$T_0>0$ に対する温度比を対数に取るので、対数の引数は無次元である。
:::
::::

::::exercise{level=advanced}
初期応力ゼロの Maxwell 模型で $\dot\eps=a$ を時刻 $0$ から加えたときの応力と長時間極限を求めよ。

:::solution
$\dot\sigma+\sigma/\tau=Ea$ に積分因子 $e^{t/\tau}$ を掛ける。積分して $\sigma=Ea\tau(1-e^{-t/\tau})=\eta a(1-e^{-t/\tau})$。$t\gg\tau$ では $\eta a$ に近づく。
:::
::::

::::exercise{level=advanced}
剛体回転を重ねるとき、対称テンソル $A$ の Jaumann 率の客観性を成分によらず証明せよ。

:::solution
$S=\dot Q Q^T$ とすると $D_tA^*=Q(D_tA)Q^T+SA^*-A^*S$ である。一方 $W^*=QWQ^T+S$ なので $D_tA^*-W^*A^*+A^*W^*=Q(D_tA-WA+AW)Q^T$ となる。観測者の回転に由来する二項が相殺される。
:::
::::

## 確認テスト

::::quiz
応力が対称になる仮定を選べ。

- [x] 角運動量を保存し、独立の偶力が作用しない。
- [ ] 温度が一定である。
- [ ] 速度場が定常である。

:::explanation
角運動量式から反対称応力だけを抽出できるのは、接触・体積偶力がない古典連続体の場合である。温度の一様性や定常性はその証明に必要ない。
:::
::::

::::quiz{answer="100" tol="1e-6"}
$\sigma=-100I\,\mathrm{kPa}$ での圧縮を正に数えた圧力を $\mathrm{kPa}$ で答えよ。

:::explanation
$p=-\tr\sigma/3=-(-300)/3=100\,\mathrm{kPa}$ である。
:::
::::

::::quiz{answer="3" tol="1e-6"}
$E=5\,\mathrm{MPa}$、$\eta=15\,\mathrm{MPa\,s}$ の時定数を秒で答えよ。

:::explanation
ばねとダッシュポットから得た時定数は $\tau=\eta/E=15/5=3\,\mathrm s$ である。
:::
::::

::::quiz
Fourier 則で熱伝導による散逸が任意の温度勾配で非負になる条件を選べ。

- [x] $k\ge0$
- [ ] $k<0$
- [ ] $\grad T$ が時間変化しない。

:::explanation
熱の寄与は $k|\grad T|^2/T$ であり $T>0$ なので、ゼロでない温度勾配も許すと必要十分条件は $k\ge0$ である。
:::
::::

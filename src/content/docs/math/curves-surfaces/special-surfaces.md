---
title: 回転面・線織面・極小曲面
description: 対称性や直線によって作られる曲面の曲率を求め、面積変分と複素関数を通じて極小曲面を構成する。
notation:
  - symbol: 'r(s),z(s)'
    meaning: 回転面の子午線の半径と高さ
  - symbol: 'H'
    meaning: 形作用素の跡の半分
  - symbol: '\Phi'
    meaning: 極小曲面の複素微分
  - symbol: 'X_d'
    meaning: 法線方向に距離 d だけずらした表示
---

## この章の目標

- 回転面と線織面の曲率を、生成する曲線の微分から求める。
- 面積の第一変分から極小曲面の方程式を導く。
- 懸垂面、螺旋面、Enneper 曲面の計量と曲率を比較する。
- Weierstrass 表現で極小曲面を作り、平行曲面や管の特異点を判定する。

## 回転面とその曲率

平面内の曲線を鉛直軸の回りに回すと、一点を通る子午線と平行圏が自然に定まる。多くの場合、二方向が主方向になるため、生成曲線の曲がり方だけで形作用素を求められる。

:::definition[回転面]{#def-cs-revolution}
弧長 $s$ で表した平面曲線 $(r(s),z(s))$ で $r>0$ とし、$X(s,\theta)=(r(s)\cos\theta,r(s)\sin\theta,z(s))$ で得る曲面を :term[回転面]{en="surface of revolution" yomi="かいてんめん"} という。$s$ を変える曲線を子午線、$\theta$ を変える曲線を平行圏という。
:::

:::proposition[回転面の主曲率]{#prop-cs-revolution-curvature}
$r'=\cos\alpha$、$z'=\sin\alpha$ とし、法線を $N=(-z'\cos\theta,-z'\sin\theta,r')$ と取る。すると

$$
\mathrm I=\dd s^2+r^2\dd\theta^2,\quad
k_s=\alpha',\quad k_\theta=\frac{z'}r,\quad
K=-\frac{r''}{r},\quad H=\tfrac12(\alpha'+z'/r).
\label{eq-cs-revolution-curvature}
$$

:::

:::proof
円筒基底を $e_r,e_\theta,e_z$ と書く。$X_s=r'e_r+z'e_z$、$X_\theta=re_\theta$ なので計量は対角で指定の法線は正規化した外積である。$X_{ss}=r''e_r+z''e_z=\alpha'N$、$X_{s\theta}=r'e_\theta$、$X_{\theta\theta}=-re_r$。従って第二基本形式の対角成分は $\alpha'$ と $rz'$、混合成分は零である。二つの形式が対角だから主曲率は成分の商となる。$r''=-\alpha'\sin\alpha=-\alpha'z'$ を使えば $K=-r''/r$ を得る。
:::

$N$ は生成曲線の進行向きに応じて外向きにも内向きにもなる。主曲率と平均曲率を比較する場合は符号を合わせる。$r=0$ の点ではこの座標は壊れるが、それだけでは形の特異点かどうかを判断できない。

:::proposition[回転軸上での滑らかさの確認]{#prop-cs-revolution-axis}
軸の近くで側面が $z=\varphi(r^2)$ と書け、$\varphi$ が $0$ の近くで滑らかなら、側面に軸上の点を加えたものはそこで正則である。上向き法線に対する軸上の二主曲率はともに $2\varphi'(0)$ である。
:::

:::proof
直交座標で $X(x,y)=(x,y,\varphi(x^2+y^2))$ とすれば、最初の二成分が座標なので原点も含めて正則である。原点で勾配は零、Hessian は $2\varphi'(0)I_2$ となる。第一基本形式は単位行列であるため、Hessian の二つの固有値が主曲率になる。同じ結論は座標特異点へ向けて子午線方向と平行圏方向の主曲率を極限に移しても得られる。
:::

例えば半径 $R$ の上半球では $\varphi(\xi)=\sqrt{R^2-\xi}$ だから、軸上の外向き法線で二主曲率は $-1/R$ である。対して $z=r$ の円錐の高さは直交座標で $\sqrt{x^2+y^2}$ となり、原点で微分できない。側面の各母線から原点へ向かう方向が一つの接平面に入らないため、別の表示を選んでもその頂点は正則にならない。

:::example[トーラスの曲率の符号]
$r=3+\cos s$、$z=\sin s$ は単位速さである。:eqref[eq-cs-revolution-curvature] より $K=\cos s/(3+\cos s)$。外周 $s=0$ で $1/4$、内周 $s=\pi$ で $-1/2$、上端と下端では零である。外向き法線は公式の $N$ の反対なので、外向きでの主曲率は $-1$ と $-\cos s/(3+\cos s)$ である。
:::

:::example[円錐の頂点を除く部分]
$r=s/\sqrt2$、$z=s/\sqrt2$、$s>0$ なら子午線は直線で $k_s=0$、$k_\theta=1/s$ である。よって $K=0$ だが $H=1/(2s)$。頂点を除けば平坦な計量を持つが、主曲率の一方は頂点に近づくにつれて発散する。
:::

:::example[放物面の中心から離れた点]
子午線を掃いてできる面を考える。角方向を第二座標に取り、指定した子午線パラメータでの曲率を計算する。子午線が等速とは仮定しない。

座標で指定された点の位置と、速度に相当する二つの接ベクトルを計算する。正則性の計算は、二本が張る平行四辺形の面積が正であることを確かめることになる。

$$
\begin{aligned}
&X(u,v)=\left(u \cos{\left(v \right)},u \sin{\left(v \right)},u^{2}\right),\\
&p=\left(1,0,1\right),\\
&(X_u,X_v)=\left(\left(1,0,2\right),\left(0,1,0\right)\right)
\end{aligned}
$$

これらを材料に内積の表と外積を出す。以下で第二基本形式を計算する際には、外積の向きと指定された法線の向きを合わせる。符号が必要な場面では法線を規格化してから代入する。

$$
\begin{aligned}
&g(p)=\begin{pmatrix}5 & 0\\0 & 1\end{pmatrix},\\
&X_u\times X_v=\left(-2,0,1\right),\\
&\sqrt{\det g}=\sqrt{5},\\
&N=\left(- \frac{2 \sqrt{5}}{5},0,\frac{\sqrt{5}}{5}\right)
\end{aligned}
$$

二階微分をする時点ではまだ座標値を固定しない。その後で二階微分を法線に射影すると、加速度のうち接平面内で向きを変える寄与が消える。

$$
\begin{aligned}X_{uu}&=\left(0,0,2\right),& X_{uv}&=\left(0,1,0\right),& X_{vv}&=\left(-1,0,0\right),\\ h(p)&=\begin{pmatrix}\frac{2 \sqrt{5}}{5} & 0\\0 & \frac{2 \sqrt{5}}{5}\end{pmatrix}, & A(p)&=\begin{pmatrix}\frac{2 \sqrt{5}}{25} & 0\\0 & \frac{2 \sqrt{5}}{5}\end{pmatrix}.\end{aligned}
$$

Gauss 曲率は二つの主曲率の積として行列式で、平均曲率は和の半分として跡で求める。主曲率自体が必要なら、この二つを係数とする二次式を解ける。

$$
\begin{aligned}
&K=\frac{\det h}{\det g}=\frac{4}{25},\\
&H=\frac12\tr A=\frac{6 \sqrt{5}}{25},\\
&\det(zI-A)=\frac{25 z^{2} - 12 \sqrt{5} z + 4}{25}
\end{aligned}
$$

法線の鉛直成分は正で、子午線方向の曲率の方が小さい。平面グラフで同じ点を計算すると、半径方向の接速度の長さが円周方向より大きい理由も分かる。
:::

:::example[放物面の中心から離れた点の曲率を子午線から見る]
同じ回転面について、角度に関する微分を含まない公式も確かめる。回転で長さと曲率が保たれるので、結果は角度に依存しない。

半径関数と高さ関数の速度を並べ、その長さを求める。法線は水平の半径成分が高さ微分の負、鉛直成分が半径微分となる。

$$
(r,r',z')=\left(u,1,2 u\right),\qquad w^2=4 u^{2} + 1
$$

子午線方向の法曲率は平面曲線の符号付き曲率と同じ形である。円周方向では円の法線加速度の大きさが半径の逆数であり、その曲面法線への射影だけが残る。

$$
k_{\mathrm{mer}}=\frac{2}{\left(4 u^{2} + 1\right)^{\frac{3}{2}}},\qquad k_{\mathrm{par}}=\frac{2}{\sqrt{4 u^{2} + 1}}
$$

指定した高さで両者を評価して掛け合わせる。実際に回転面で曲率の符号を見たい場合、子午線がどちらに曲がるかと、その場所の傾きを別々に読むとよい。

$$
k_{\mathrm{mer}}(u_0)=\frac{2 \sqrt{5}}{25},\qquad k_{\mathrm{par}}(u_0)=\frac{2 \sqrt{5}}{5},\qquad K(u_0)=\frac{4}{25}
$$

非等速の座標では、子午線方向の二階微分をそのまま法曲率と読めない。速さの三乗が出るのは、外積に一次微分がもう一つ含まれるためである。
:::

:::example[円錐と同じ半径の平行圏]
子午線を掃いてできる面を考える。角方向を第二座標に取り、指定した子午線パラメータでの曲率を計算する。子午線が等速とは仮定しない。

座標で指定された点の位置と、速度に相当する二つの接ベクトルを計算する。正則性の計算は、二本が張る平行四辺形の面積が正であることを確かめることになる。

$$
\begin{aligned}
&X(u,v)=\left(2 u \cos{\left(v \right)},2 u \sin{\left(v \right)},3 u\right),\\
&p=\left(2,0,3\right),\\
&(X_u,X_v)=\left(\left(2,0,3\right),\left(0,2,0\right)\right)
\end{aligned}
$$

これらを材料に内積の表と外積を出す。以下で第二基本形式を計算する際には、外積の向きと指定された法線の向きを合わせる。符号が必要な場面では法線を規格化してから代入する。

$$
\begin{aligned}
&g(p)=\begin{pmatrix}13 & 0\\0 & 4\end{pmatrix},\\
&X_u\times X_v=\left(-6,0,4\right),\\
&\sqrt{\det g}=2 \sqrt{13},\\
&N=\left(- \frac{3 \sqrt{13}}{13},0,\frac{2 \sqrt{13}}{13}\right)
\end{aligned}
$$

二階微分をする時点ではまだ座標値を固定しない。その後で二階微分を法線に射影すると、加速度のうち接平面内で向きを変える寄与が消える。

$$
\begin{aligned}X_{uu}&=\left(0,0,0\right),& X_{uv}&=\left(0,2,0\right),& X_{vv}&=\left(-2,0,0\right),\\ h(p)&=\begin{pmatrix}0 & 0\\0 & \frac{6 \sqrt{13}}{13}\end{pmatrix}, & A(p)&=\begin{pmatrix}0 & 0\\0 & \frac{3 \sqrt{13}}{26}\end{pmatrix}.\end{aligned}
$$

Gauss 曲率は二つの主曲率の積として行列式で、平均曲率は和の半分として跡で求める。主曲率自体が必要なら、この二つを係数とする二次式を解ける。

$$
\begin{aligned}
&K=\frac{\det h}{\det g}=0,\\
&H=\frac12\tr A=\frac{3 \sqrt{13}}{52},\\
&\det(zI-A)=\frac{z \left(26 z - 3 \sqrt{13}\right)}{26}
\end{aligned}
$$

母線方向は直線なので曲率が消え、面の Gauss 曲率も零である。それでも平行圏の方向には正の法曲率がある。
:::

:::example[円錐と同じ半径の平行圏の曲率を子午線から見る]
同じ回転面について、角度に関する微分を含まない公式も確かめる。回転で長さと曲率が保たれるので、結果は角度に依存しない。

半径関数と高さ関数の速度を並べ、その長さを求める。法線は水平の半径成分が高さ微分の負、鉛直成分が半径微分となる。

$$
(r,r',z')=\left(2 u,2,3\right),\qquad w^2=13
$$

子午線方向の法曲率は平面曲線の符号付き曲率と同じ形である。円周方向では円の法線加速度の大きさが半径の逆数であり、その曲面法線への射影だけが残る。

$$
k_{\mathrm{mer}}=0,\qquad k_{\mathrm{par}}=\frac{3 \sqrt{13}}{26 u}
$$

指定した高さで両者を評価して掛け合わせる。実際に回転面で曲率の符号を見たい場合、子午線がどちらに曲がるかと、その場所の傾きを別々に読むとよい。

$$
k_{\mathrm{mer}}(u_0)=0,\qquad k_{\mathrm{par}}(u_0)=\frac{3 \sqrt{13}}{26},\qquad K(u_0)=0
$$

非等速の座標では、子午線方向の二階微分をそのまま法曲率と読めない。速さの三乗が出るのは、外積に一次微分がもう一つ含まれるためである。
:::

## 線織面と可展面

直線を動かして曲面を作ると、その直線に沿う法曲率は零である。従って二方向に同符号で曲がる点はない。混合成分を計算すると、曲率が零となる場合も判定できる。

:::definition[線織面と可展面]{#def-cs-ruled}
$X(u,v)=c(u)+v\,d(u)$、$\norm d=1$ と表せる正則な曲面を :term[線織面]{en="ruled surface" yomi="せんしょくめん"} という。$v$ を動かす直線を母線という。ここでは $K\equiv0$ の線織面を :term[可展面]{en="developable surface" yomi="かてんめん"} と呼ぶ。
:::

:::theorem[線織面の曲率]{#thm-cs-ruled-curvature}
正則な線織面で

$$
K=-\frac{\det(c',d,d')^2}{\norm{(c'+vd')\times d}^4}\le0.
\label{eq-cs-ruled-curvature}
$$

特に可展であることは $\det(c',d,d')=0$ と同値であり、そのとき単位法線は母線に沿って一定である。
:::

:::proof
$X_{vv}=0$ であり、$X_{uv}=d'$ だから第二基本形式の最後の成分は零、混合成分は $\det(c',d,d')/\norm{X_u\times X_v}$ である。第一基本形式の行列式は $\norm{X_u\times X_v}^2$ なので公式を得る。可展の場合、外積 $(c'+vd')\times d$ を構成する $c'\times d$ と $d'\times d$ は平行である。実際、$d'=0$ なら自明で、そうでなければ $c'\in\spn(d,d')$ である。各母線の正則な連結部分上で外積は零にならず、その単位方向は一定である。
:::

$d'\ne0$ の可展面では $c'=a d+b d'$ と書ける。$C=c-bd$ と置くと $C'=(a-b')d$ なので、$a-b'\ne0$ の区間では母線は $C$ の接線である。$a-b'$ が恒等的に零の区間では $C$ が一点となり、曲面は円錐型である。$d$ が一定なら柱面型である。$a-b'$ の零点を含む場合には、境界を含めてこれらの型に分類できるとまではいえない。

:::example[接線曲面]
単位速さ曲線 $\alpha$ から $X(s,v)=\alpha(s)+vT(s)$ を作ると $X_s=T+v\kappa N$ なので $X_s\times X_v=-v\kappa B$。$\kappa>0$ で $v\ne0$ の範囲では正則、かつ $\det(T,T,\kappa N)=0$ で可展である。$v=0$ で座標が退化する。
:::

:::example[螺旋面の母線]
$X(u,v)=(v\cos u,v\sin u,u)$ では $c'=(0,0,1)$、$d=(\cos u,\sin u,0)$、$d'=(-\sin u,\cos u,0)$ なので行列式は $1$。外積のノルムの二乗は $1+v^2$ であり、$K=-1/(1+v^2)^2$ となる。
:::

### 可展面を平面へ展開する

曲率が零である線織面では、一本の母線の上で接平面が変わらない。従って母線を横切る区間に沿って法線がどう回るかを捉えれば、可展性を判断できる。接線方向を母線に選ぶと曲線からも可展面が得られる。

:::proposition[線織表示を用いた局所展開]{#prop-cs-developable-flattening}
可展な線織表示で $d'$ が非零の区間、または $d$ が一定の区間では、各正則点の近傍を第一基本形式を保って平面に写せる。
:::

:::proof
$d'\ne0$ のときは単位球面上の $d$ を弧長 $u$ で表示し直す。可展性により $c'=a d+b d'$ と書ける。$X_u=a d+(b+v)d'$、$X_v=d$ であり、$d,d'$ は正規直交である。平面では $D(u)=(\cos u,\sin u)$ と置き、$C'=aD+bD'$ を積分し、$Y=C+vD$ と定める。$Y_u,Y_v$ は $D,D'$ に関して $X_u,X_v$ と同じ成分を持つ。従って第一基本形式は一致し、正則点では $Y$ も階数二なので逆関数定理により平面の座標になる。

$d$ が一定のときは $a=\inner{c'}d$、$b=\norm{c'-ad}>0$ と置く。$b>0$ は階数条件による。平面の表示を $Y=(v+\int a\,\dd u,\int b\,\dd u)$ とすれば、座標速度の内積の表は $X$ と同じになる。
:::

例えば $X(s,\theta)=(s\cos\theta/\sqrt2,s\sin\theta/\sqrt2,s/\sqrt2)$ では、第一基本形式は $\dd s^2+\tfrac12 s^2\dd\theta^2$ である。平面の極座標を $(s,\theta/\sqrt2)$ と取ればこの式と一致する。切れ目を入れて $\theta$ を一周分に限ると開く扇形の角度は $\sqrt2\pi$ であり、角度を一周分貼り合わせるかどうかで大域的な形が変わる。

:::example[空間曲線の接線曲面]
単位速さ曲線 $\alpha$ に対し $X(s,v)=\alpha(s)+vT(s)$ とする。$\kappa>0$、$v\ne0$ の領域では $X_s=T+v\kappa N$、$X_v=T$ なので $X_s\times X_v=-v\kappa B$ となる。法線は $v$ の符号が一定なら母線に沿って一定である。従って形作用素は $X_v$ を消し、Gauss 曲率は零である。$v=0$ では二つの接ベクトルが一致し、この表示の階数が一に落ちる。元の曲線が正則であっても、接線曲面までその曲線上で正則になるわけではない。

例えば $(2\cos t,2\sin t,t)$ の弧長を $s=\sqrt5t$ とすると、$s=0$ で $B=(0,-1,2)/\sqrt5$ である。$v>0$ の母線の法線は $(0,1,-2)/\sqrt5$ であり、その母線に接する平面はどの $v$ でも $y-2z=0$ である。$s$ を変えると $B$ が変化するので、曲面全体がその平面に含まれるわけではない。
:::

:::example[母線方向が変わっても可展な場合]
原点を頂点とする円錐を線織面として扱う。頂点を除き、各母線は一つの接平面を保つか確かめる。

座標で指定された点の位置と、速度に相当する二つの接ベクトルを計算する。正則性の計算は、二本が張る平行四辺形の面積が正であることを確かめることになる。

$$
\begin{aligned}
&X(u,v)=\left(v \cos{\left(u \right)},v \sin{\left(u \right)},v\right),\\
&p=\left(2,0,2\right),\\
&(X_u,X_v)=\left(\left(0,2,0\right),\left(1,0,1\right)\right)
\end{aligned}
$$

これらを材料に内積の表と外積を出す。以下で第二基本形式を計算する際には、外積の向きと指定された法線の向きを合わせる。符号が必要な場面では法線を規格化してから代入する。

$$
\begin{aligned}
&g(p)=\begin{pmatrix}4 & 0\\0 & 2\end{pmatrix},\\
&X_u\times X_v=\left(2,0,-2\right),\\
&\sqrt{\det g}=2 \sqrt{2},\\
&N=\left(\frac{\sqrt{2}}{2},0,- \frac{\sqrt{2}}{2}\right)
\end{aligned}
$$

二階微分をする時点ではまだ座標値を固定しない。その後で二階微分を法線に射影すると、加速度のうち接平面内で向きを変える寄与が消える。

$$
\begin{aligned}X_{uu}&=\left(-2,0,0\right),& X_{uv}&=\left(0,1,0\right),& X_{vv}&=\left(0,0,0\right),\\ h(p)&=\begin{pmatrix}- \sqrt{2} & 0\\0 & 0\end{pmatrix}, & A(p)&=\begin{pmatrix}- \frac{\sqrt{2}}{4} & 0\\0 & 0\end{pmatrix}.\end{aligned}
$$

Gauss 曲率は二つの主曲率の積として行列式で、平均曲率は和の半分として跡で求める。主曲率自体が必要なら、この二つを係数とする二次式を解ける。

$$
\begin{aligned}
&K=\frac{\det h}{\det g}=0,\\
&H=\frac12\tr A=- \frac{\sqrt{2}}{8},\\
&\det(zI-A)=\frac{z \left(4 z + \sqrt{2}\right)}{4}
\end{aligned}
$$

形作用素で零の方向は母線である。固定した $u$ で距離 $v$ を変えても、単位法線は変わらない。母線の回転を追うだけでは可展性を否定できない。
:::

:::example[牛角状の鞍面の直線を追う]
二つの座標線が直線である曲面を調べ、母線の近くで接平面が回るかを見る。

座標で指定された点の位置と、速度に相当する二つの接ベクトルを計算する。正則性の計算は、二本が張る平行四辺形の面積が正であることを確かめることになる。

$$
\begin{aligned}
&X(u,v)=\left(u,v,2 u v\right),\\
&p=\left(0,1,0\right),\\
&(X_u,X_v)=\left(\left(1,0,2\right),\left(0,1,0\right)\right)
\end{aligned}
$$

これらを材料に内積の表と外積を出す。以下で第二基本形式を計算する際には、外積の向きと指定された法線の向きを合わせる。符号が必要な場面では法線を規格化してから代入する。

$$
\begin{aligned}
&g(p)=\begin{pmatrix}5 & 0\\0 & 1\end{pmatrix},\\
&X_u\times X_v=\left(-2,0,1\right),\\
&\sqrt{\det g}=\sqrt{5},\\
&N=\left(- \frac{2 \sqrt{5}}{5},0,\frac{\sqrt{5}}{5}\right)
\end{aligned}
$$

二階微分をする時点ではまだ座標値を固定しない。その後で二階微分を法線に射影すると、加速度のうち接平面内で向きを変える寄与が消える。

$$
\begin{aligned}X_{uu}&=\left(0,0,0\right),& X_{uv}&=\left(0,0,2\right),& X_{vv}&=\left(0,0,0\right),\\ h(p)&=\begin{pmatrix}0 & \frac{2 \sqrt{5}}{5}\\\frac{2 \sqrt{5}}{5} & 0\end{pmatrix}, & A(p)&=\begin{pmatrix}0 & \frac{2 \sqrt{5}}{25}\\\frac{2 \sqrt{5}}{5} & 0\end{pmatrix}.\end{aligned}
$$

Gauss 曲率は二つの主曲率の積として行列式で、平均曲率は和の半分として跡で求める。主曲率自体が必要なら、この二つを係数とする二次式を解ける。

$$
\begin{aligned}
&K=\frac{\det h}{\det g}=- \frac{4}{25},\\
&H=\frac12\tr A=0,\\
&\det(zI-A)=\frac{\left(5 z - 2\right) \left(5 z + 2\right)}{25}
\end{aligned}
$$

母線上で $u=0$ のまま $v$ を動かすと、法線の第一成分は $v$ に応じて変わる。行列式が負であり、直線が通ることと可展であることを区別できる。
:::

## 極小曲面と面積の第一変分

膜の境界を固定し、内部を少しだけ動かして面積を比べる。そのとき一階の変化を決めるのは平均曲率である。曲面に沿う変位は座標のずれに相当し、境界を固定すると一階の面積変化を作らない。

:::theorem[面積の第一変分]{#thm-cs-first-variation-area}
向き付けられた曲面上の、滑らかな境界を持つコンパクト領域 $S$ を滑らかなはめ込みの族 $X_\eps$ で変分する。面積はその族が引き戻す面積密度を固定した定義域上で積分したものとする。初期の変位を $V=W+\varphi N$ と分解し、$W$ を接ベクトル場とする。面積の微分は

$$
\left.\dv{A(X_\eps)}{\eps}\right|_{\eps=0}
=\int_S(\divg_S W-2H\varphi)\,\dd A.
\label{eq-cs-area-variation}
$$

ここで $\divg_S W$ は接ベクトル場の接方向微分の跡である。境界を固定する場合は発散の積分は零である。
:::

:::proof
計量成分 $g_{ij}(\eps)=\inner{\partial_iX_\eps}{\partial_jX_\eps}$ の微分は $\dot g_{ij}=\inner{\partial_iV}{X_j}+\inner{X_i}{\partial_jV}$ となる。行列式を微分すると $\partial_\eps\sqrt{\det g}=\tfrac12\sqrt{\det g}\sum_{i,j}g^{ij}\dot g_{ij}$。法線変位では $\inner{\partial_i(\varphi N)}{X_j}=-\varphi h_{ij}$ だから寄与は $-2H\varphi$ となる。接変位では縮約 $\sum g^{ij}\inner{\partial_iW}{X_j}$ が接方向微分の跡である。

座標で $W=W^k X_k$ と書くと、この縮約は $\partial_iW^i+W^k\sum_{i,j}g^{ij}\inner{X_{ik}}{X_j}$ となる。内積を微分して跡を取ると和は $\partial_k\log\sqrt{\det g}$ に等しい。従って $\divg_S W\,\sqrt{\det g}=\partial_i(\sqrt{\det g}\,W^i)$ である。これを積分したものは境界の流束であり、固定境界では零である。複数座標の場合は分割の一を使うと各座標で境界項を合わせられる。分割関数の微分の和が零なので余分な内部項は消える。
:::

:::corollary[極小曲面の条件]{#cor-cs-minimal-condition}
曲面の内部にコンパクト台を持つすべての変分に対して面積の一階微分が零であることと $H=0$ は同値である。面積が無限の場合は、変分の台を含むコンパクト領域上で面積の変化を測る。この条件を満たす曲面を :term[極小曲面]{en="minimal surface" yomi="きょくしょうきょくめん"} という。
:::

:::proof
$H=0$ なら :eqref[eq-cs-area-variation] から微分は零。逆に $H$ が非零の点があれば、そこで符号が一定となる近傍を取り、その中に台を持つ非負で非零の $\varphi$ を選ぶ。法線変分 $X+\eps\varphi N$ に対する面積微分は非零であり、仮定に反する。
:::

定義は停留性を表す。面積が全ての比較曲面より小さいかを決めるには、変分のクラスと境界を指定し、別の議論が要る。

:::proposition[極小曲面の法線変分の第二変分]{#prop-cs-minimal-area-second-variation}
極小曲面のコンパクト領域 $S$ 上で、滑らかな関数 $\varphi$ が境界で零になるとする。$X_\eps=X+\eps\varphi N$ の面積について

$$
A''(0)=\int_S\bigl(\norm{\grad_S\varphi}^2-(k_1^2+k_2^2)\varphi^2\bigr)\,\dd A
$$

が成り立つ。ここで $\grad_S\varphi$ は接平面の内積で $\inner{\grad_S\varphi}{w}=d\varphi(w)$ を満たすベクトルである。
:::

:::proof
一点で正規直交な主基底 $e_1,e_2$ を取る。$\varphi_i=d\varphi(e_i)$ と置くと、変分した写像の微分は $e_i$ を $(1-\eps\varphi k_i)e_i+\eps\varphi_iN$ に送る。二つの像が張る平行四辺形の面積の二乗は

$$
(1-\eps\varphi k_1)^2(1-\eps\varphi k_2)^2
+\eps^2\varphi_1^2(1-\eps\varphi k_2)^2
+\eps^2\varphi_2^2(1-\eps\varphi k_1)^2
$$

となる。外積の三成分が直交するためである。$k_1+k_2=0$ を使って平方根を展開すると、面積倍率は

$$
1+\eps^2\left(\varphi^2 k_1k_2+\tfrac12(\varphi_1^2+\varphi_2^2)\right)+\order{\eps^3}
$$

である。コンパクトな領域で密度は $\eps$ に関して滑らかなので積分の下で微分できる。二階微分を取り、$2k_1k_2=-(k_1^2+k_2^2)$ を代入すれば結論を得る。
:::

:::example[極小グラフの方程式]
$X=(u,v,q)$ の面積は $\int\sqrt{1+q_u^2+q_v^2}\,\dd u\,\dd v$。$q+\eps\psi$ に変分して部分積分すると停留条件は

$$
\partial_u(q_u/W)+\partial_v(q_v/W)=0,\qquad W=\sqrt{1+q_u^2+q_v^2}
$$

である。分母を払うと :eqref[eq-cs-graph-curvature] の $H$ の分子と一致する。固定境界値をもつ別のグラフと比較すると、関数 $p\mapsto\sqrt{1+\norm p^2}$ は厳密凸で、その接平面が与える下界を積分した一次項は方程式で消える。従って極小グラフは同じ領域と境界値のグラフの中で面積を最小にする。
:::

:::example[小球片の法線変分]
半径 $2$ の球を外向きに向き付ければ $H=-1/2$。コンパクト台の外向き変位 $\varphi$ に対して面積の微分は $\int\varphi\,\dd A$ であり、$\varphi\ge0$ が非零なら正である。球の半径全体を増す場合も $\dv{}{R}(4\pi R^2)=8\pi R$ となり、式と一致する。
:::

### グラフの面積汎関数と極小曲面方程式

高さだけを動かせるグラフでは、面積の式から直接停留条件を計算できる。この計算には曲面のパラメータの変化が要らないので、境界で高さを固定した問題に使いやすい。

:::proposition[グラフの極小曲面方程式]
領域 $\Omega$ 上の滑らかなグラフ $z=f(u,v)$ に対し、台がコンパクトな任意の垂直変分で面積の第一変分が零であるための必要十分条件は
$$
(1+f_v^2) f_{uu}-2 f_u f_v f_{uv}+(1+f_u^2)f_{vv}=0
$$
である。
:::

:::proof
グラフの面積は $\int_\Omega\sqrt{1+|\grad f|^2}\,\dd u\,\dd v$ である。$f+\eps\eta$ を代入して $\eps=0$ で微分すると、積分の中は $\inner{\grad f}{\grad\eta}/\sqrt{1+|\grad f|^2}$ となる。台が境界に届かないため部分積分で境界項が消え、$\eta$ には $-\divg(\grad f/\sqrt{1+|\grad f|^2})$ が掛かる。これが各点で零でなければ、符号が一定の小領域に台を持つ非負関数を取って第一変分を零でなくできる。発散を微分して共通分母 $(1+|\grad f|^2)^{3/2}$ を掛けると式を得る。逆向きは同じ積分表示から従う。
:::

:::example[加法的な極小グラフを作る]
$f(u,v)=\log\cos v-\log\cos u$ を $|u|,|v|<\pi/2$ で取る。$f_u=\tan u$、$f_v=-\tan v$、混合微分は零である。二つの純粋二階微分は $\sec^2u,-\sec^2v$ なので、極小曲面方程式の左右の寄与は $\sec^2u\sec^2v-\sec^2u\sec^2v=0$ となる。原点では勾配が零で、主曲率は $1,-1$ である。従って平均曲率は零だが Gauss 曲率は $-1$ であり、平面にはならない。

さらに $p,q\in\R^2$ に対し関数 $L(p)=\sqrt{1+|p|^2}$ の Hessian は $(1+|p|^2)^{-1/2}I-(1+|p|^2)^{-3/2}pp^\top$ である。任意の非零方向への二次形式は Cauchy–Schwarz により正なので、$L(q)\ge L(p)+dL_p(q-p)$ となる。これをグラフの面積に適用し、境界を固定して部分積分すると、上のグラフは同じ境界値を持つグラフの間で面積を最小にする。この比較はグラフという比較対象を指定しているため、停留条件の名称から最小性を推測する議論とは異なる。
:::

:::example[非停留なグラフの曲がりを見つける]
原点の上向きの法線速度を持つ変形で、面積の第一変分の被積分関数を求める。

座標で指定された点の位置と、速度に相当する二つの接ベクトルを計算する。正則性の計算は、二本が張る平行四辺形の面積が正であることを確かめることになる。

$$
\begin{aligned}
&X(u,v)=\left(u,v,u^{2} + 2 v^{2}\right),\\
&p=\left(0,0,0\right),\\
&(X_u,X_v)=\left(\left(1,0,0\right),\left(0,1,0\right)\right)
\end{aligned}
$$

これらを材料に内積の表と外積を出す。以下で第二基本形式を計算する際には、外積の向きと指定された法線の向きを合わせる。符号が必要な場面では法線を規格化してから代入する。

$$
\begin{aligned}
&g(p)=\begin{pmatrix}1 & 0\\0 & 1\end{pmatrix},\\
&X_u\times X_v=\left(0,0,1\right),\\
&\sqrt{\det g}=1,\\
&N=\left(0,0,1\right)
\end{aligned}
$$

二階微分をする時点ではまだ座標値を固定しない。その後で二階微分を法線に射影すると、加速度のうち接平面内で向きを変える寄与が消える。

$$
\begin{aligned}X_{uu}&=\left(0,0,2\right),& X_{uv}&=\left(0,0,0\right),& X_{vv}&=\left(0,0,4\right),\\ h(p)&=\begin{pmatrix}2 & 0\\0 & 4\end{pmatrix}, & A(p)&=\begin{pmatrix}2 & 0\\0 & 4\end{pmatrix}.\end{aligned}
$$

Gauss 曲率は二つの主曲率の積として行列式で、平均曲率は和の半分として跡で求める。主曲率自体が必要なら、この二つを係数とする二次式を解ける。

$$
K=\frac{\det h}{\det g}=8,\qquad H=\frac12\tr A=3,\qquad \det(zI-A)=\left(z - 4\right) \left(z - 2\right)
$$

原点で $-2H=-6$ であり、同符号の近傍に正で台がコンパクトな関数を選べば上方への変形で面積が減る。境界を動かさずに面積を変化させられる。
:::

:::example[平均の打消しと全ての方向の曲がり]
原点で平均曲率が消えることを計算する。近傍のどの点でも消えるかを別に確認する。

座標で指定された点の位置と、速度に相当する二つの接ベクトルを計算する。正則性の計算は、二本が張る平行四辺形の面積が正であることを確かめることになる。

$$
\begin{aligned}
&X(u,v)=\left(u,v,\left(u - v\right) \left(u + v\right)\right),\\
&p=\left(0,0,0\right),\\
&(X_u,X_v)=\left(\left(1,0,0\right),\left(0,1,0\right)\right)
\end{aligned}
$$

これらを材料に内積の表と外積を出す。以下で第二基本形式を計算する際には、外積の向きと指定された法線の向きを合わせる。符号が必要な場面では法線を規格化してから代入する。

$$
\begin{aligned}
&g(p)=\begin{pmatrix}1 & 0\\0 & 1\end{pmatrix},\\
&X_u\times X_v=\left(0,0,1\right),\\
&\sqrt{\det g}=1,\\
&N=\left(0,0,1\right)
\end{aligned}
$$

二階微分をする時点ではまだ座標値を固定しない。その後で二階微分を法線に射影すると、加速度のうち接平面内で向きを変える寄与が消える。

$$
\begin{aligned}X_{uu}&=\left(0,0,2\right),& X_{uv}&=\left(0,0,0\right),& X_{vv}&=\left(0,0,-2\right),\\ h(p)&=\begin{pmatrix}2 & 0\\0 & -2\end{pmatrix}, & A(p)&=\begin{pmatrix}2 & 0\\0 & -2\end{pmatrix}.\end{aligned}
$$

Gauss 曲率は二つの主曲率の積として行列式で、平均曲率は和の半分として跡で求める。主曲率自体が必要なら、この二つを係数とする二次式を解ける。

$$
K=\frac{\det h}{\det g}=-4,\qquad H=\frac12\tr A=0,\qquad \det(zI-A)=\left(z - 2\right) \left(z + 2\right)
$$

例えば座標 $(1,0)$ では同じグラフの平均曲率は $-4/(5\sqrt5)$。従って原点で被積分関数が零でも、曲面全体の停留性は得られない。
:::

## 極小曲面の例（懸垂面・螺旋面・Enneper 曲面）

まず回転面で $H=0$ を解く。続いて回転対称性の代わりに螺旋対称性を持つ曲面を調べる。さらに調和多項式から対称性の異なる例を作れる。

:::example[懸垂面と回転極小面]
半径を高さの関数 $r(z)>0$ として $X=(r(z)\cos u,r(z)\sin u,z)$ と置く。子午線の速さは $\sqrt{1+r_z^2}$、接線角の弧長微分は $-r_{zz}/(1+r_z^2)^{3/2}$、$z'=(1+r_z^2)^{-1/2}$ である。:eqref[eq-cs-revolution-curvature] により $H=0$ は $r r_{zz}=1+r_z^2$ と同値になる。$r/\sqrt{1+r_z^2}$ を微分するとこの方程式で零になり、この正の定数を $a$ とする。従って $r_{zz}=r/a^2$ となり、一般解と $r^2=a^2(1+r_z^2)$ から $r=a\cosh((z-z_0)/a)$ を得る。実際、$r=Ae^{z/a}+Be^{-z/a}$ では $4AB=a^2$ で、$r>0$ から $A,B>0$ である。

この回転面を :term[懸垂面]{en="catenoid" yomi="けんすいめん"} という。$z=z_0+av$ を使うと $\mathrm I=a^2\cosh^2v(\dd u^2+\dd v^2)$。二主曲率は互いに反対で絶対値 $1/(a\cosh^2 v)$ だから $K=-1/(a^2\cosh^4v)$ である。高さが定値の回転極小面は平面の一部である。非定値の場合、上の高さ座標が取れる区間で懸垂面となり、そこで $z'=\pm a/r$ だから有限の $r>0$ で高さが停留する点に延長することはない。
:::

:::example[螺旋面]
$X(u,v)=(v\cos u,v\sin u,au)$、$a>0$ を :term[螺旋面]{en="helicoid" yomi="らせんめん"} という。第一基本形式は $(v^2+a^2)\dd u^2+\dd v^2$。$X_{vv}=0$ であり、$X_{uu}$ は母線方向なので法線成分は零である。従って第二基本形式は混合成分だけからなり、$H=0$ となる。また :eqref[eq-cs-ruled-curvature] より $K=-a^2/(a^2+v^2)^2$。$v=a\sinh w$ と置けば計量は $a^2\cosh^2w(\dd u^2+\dd w^2)$ であり、懸垂面の計量と局所的に一致する。
:::

:::example[Enneper 曲面]
$z=u+iv$ として

$$
X(u,v)=\Re\left(z-\frac{z^3}{3},\,i\left(z+\frac{z^3}{3}\right),\,z^2\right)
$$

を考える。括弧内の複素微分は $\Phi=(1-z^2,i(1+z^2),2z)$ で、$\Phi\cdot\Phi=0$、$\norm{\Phi}_{\C}^2=2(1+|z|^2)^2$。従って実部と虚部は等しい正のノルムを持ち直交し、$X$ ははめ込みで計量は $(1+|z|^2)^2(\dd u^2+\dd v^2)$ である。各成分は正則関数の実部なので調和であり、$X_{uu}+X_{vv}=0$ の法線成分を取れば $H=0$ となる。原点で $X_u=(1,0,0)$、$X_v=(0,-1,0)$ だから右手系の法線は $(0,0,-1)$ で、主曲率は $-2,2$、$K=-4$ である。
:::

:::example[螺旋面の軸上の混合項]
螺旋面が軸上でも正則であることと、その点での曲率を同じ表示で確認する。

座標で指定された点の位置と、速度に相当する二つの接ベクトルを計算する。正則性の計算は、二本が張る平行四辺形の面積が正であることを確かめることになる。

$$
\begin{aligned}
&X(u,v)=\left(u \cos{\left(v \right)},u \sin{\left(v \right)},2 v\right),\\
&p=\left(0,0,0\right),\\
&(X_u,X_v)=\left(\left(1,0,0\right),\left(0,0,2\right)\right)
\end{aligned}
$$

これらを材料に内積の表と外積を出す。以下で第二基本形式を計算する際には、外積の向きと指定された法線の向きを合わせる。符号が必要な場面では法線を規格化してから代入する。

$$
\begin{aligned}
&g(p)=\begin{pmatrix}1 & 0\\0 & 4\end{pmatrix},\\
&X_u\times X_v=\left(0,-2,0\right),\\
&\sqrt{\det g}=2,\\
&N=\left(0,-1,0\right)
\end{aligned}
$$

二階微分をする時点ではまだ座標値を固定しない。その後で二階微分を法線に射影すると、加速度のうち接平面内で向きを変える寄与が消える。

$$
\begin{aligned}X_{uu}&=\left(0,0,0\right),& X_{uv}&=\left(0,1,0\right),& X_{vv}&=\left(0,0,0\right),\\ h(p)&=\begin{pmatrix}0 & -1\\-1 & 0\end{pmatrix}, & A(p)&=\begin{pmatrix}0 & -1\\- \frac{1}{4} & 0\end{pmatrix}.\end{aligned}
$$

Gauss 曲率は二つの主曲率の積として行列式で、平均曲率は和の半分として跡で求める。主曲率自体が必要なら、この二つを係数とする二次式を解ける。

$$
\begin{aligned}
&K=\frac{\det h}{\det g}=- \frac{1}{4},\\
&H=\frac12\tr A=0,\\
&\det(zI-A)=\frac{\left(2 z - 1\right) \left(2 z + 1\right)}{4}
\end{aligned}
$$

二つの対角成分が零に見えるのに曲率が消えないのは、第二基本形式に交差項があるためである。直交する適切な二方向を取ると、符号の逆の主曲率 $1/2,-1/2$ が現れる。
:::

:::example[懸垂面の首を面積密度から検算する]
縮尺 $2$ の懸垂面を高さが第一座標となるように表す。首の点では二本の接方向が何を表すか確認する。

座標で指定された点の位置と、速度に相当する二つの接ベクトルを計算する。正則性の計算は、二本が張る平行四辺形の面積が正であることを確かめることになる。

$$
\begin{aligned}
&X(u,v)=\left(2 \cos{\left(v \right)} \cosh{\left(\frac{u}{2} \right)},2 \sin{\left(v \right)} \cosh{\left(\frac{u}{2} \right)},u\right),\\
&p=\left(2,0,0\right),\\
&(X_u,X_v)=\left(\left(0,0,1\right),\left(0,2,0\right)\right)
\end{aligned}
$$

これらを材料に内積の表と外積を出す。以下で第二基本形式を計算する際には、外積の向きと指定された法線の向きを合わせる。符号が必要な場面では法線を規格化してから代入する。

$$
\begin{aligned}
&g(p)=\begin{pmatrix}1 & 0\\0 & 4\end{pmatrix},\\
&X_u\times X_v=\left(-2,0,0\right),\\
&\sqrt{\det g}=2,\\
&N=\left(-1,0,0\right)
\end{aligned}
$$

二階微分をする時点ではまだ座標値を固定しない。その後で二階微分を法線に射影すると、加速度のうち接平面内で向きを変える寄与が消える。

$$
\begin{aligned}X_{uu}&=\left(\frac{1}{2},0,0\right),& X_{uv}&=\left(0,0,0\right),& X_{vv}&=\left(-2,0,0\right),\\ h(p)&=\begin{pmatrix}- \frac{1}{2} & 0\\0 & 2\end{pmatrix}, & A(p)&=\begin{pmatrix}- \frac{1}{2} & 0\\0 & \frac{1}{2}\end{pmatrix}.\end{aligned}
$$

Gauss 曲率は二つの主曲率の積として行列式で、平均曲率は和の半分として跡で求める。主曲率自体が必要なら、この二つを係数とする二次式を解ける。

$$
\begin{aligned}
&K=\frac{\det h}{\det g}=- \frac{1}{4},\\
&H=\frac12\tr A=0,\\
&\det(zI-A)=\frac{\left(2 z - 1\right) \left(2 z + 1\right)}{4}
\end{aligned}
$$

首の上では一方向で内側、他の方向で外側に曲がる。両主曲率の大きさは $1/2$ であるため平均曲率が消える。座標の順序からここでの法線は内向きである。
:::

:::example[懸垂面と同じ境界を持つ二円板]
半径 $r(z)=\cosh z$ の懸垂面を $-2\le z\le2$ に制限する。二境界は半径 $\cosh 2$ の円である。懸垂面の計量から面積は

$$
2\pi\int_{-2}^2\cosh^2 z\,\dd z=4\pi+\pi\sinh 4
$$

となる。一方、それぞれの境界を張る水平円板の面積の和は $2\pi\cosh^2 2=\pi(1+\cosh4)$ である。懸垂面の面積からこの和を引くと $\pi(3-e^{-4})>0$ となる。従って境界を張る曲面として非連結な曲面も許すと、極小な懸垂面より面積が小さい比較対象がある。
:::

:::example[境界を固定した小さな変形で面積が減る懸垂面]
同じ懸垂面で $\varphi(v)=1-v^2/4$ と置く。これは $v=\pm2$ で零になる。:ref[prop-cs-minimal-area-second-variation] で等温計量 $\cosh^2v(\dd u^2+\dd v^2)$ と主曲率 $\pm1/\cosh^2v$ を用いると

$$
\frac{A''(0)}{2\pi}
=\int_{-2}^2\left((\varphi')^2-\frac{2\varphi^2}{\cosh^2v}\right)\,\dd v.
$$

第一項の積分は $4/3$ である。第二項の絶対値の積分は、$-1\le v\le1$ で $\varphi\ge3/4$ だから $(9/4)\tanh1$ 以上である。$e^2>5$ より $\tanh1>2/3$ であり、この下界は $4/3$ より大きい。よって $A''(0)<0$ となる。第一変分は零なので、小さい非零の $\eps$ で境界を動かさずに面積が減る。
:::

## 等温座標と Weierstrass 表現

前節では計量を座標平面の内積の定数倍に見える形に直した。倍率が点ごとに変わってよいとすると、面積の停留条件を通常の調和方程式として扱える。

:::definition[等温座標]{#def-cs-isothermal}
第一基本形式が $\lambda^2(\dd u^2+\dd v^2)$、$\lambda>0$ となる座標を :term[等温座標]{en="isothermal coordinates" yomi="とうおんざひょう"} という。
:::

### 等温座標の存在

計量の係数を同時に簡単にするには、一つ目の座標をその計量に関して調和に選ぶ。二つ目は、一つ目の微分を余接平面で直角に回して積分すると得られる。以下では楕円型方程式の解の存在を前提にし、その解が座標になることを確かめる。

:::theorem[等温座標の局所存在]{#thm-cs-isothermal-existence}
滑らかな二次元の正定値計量には、任意の点の近くで等温座標が存在する。
:::

:::proof
局所座標の原点で $g_{ij}(0)=\delta_{ij}$ となるよう最初に線形変換する。滑らかな関数 $f$ の星作用素 $*\dd f$ は向きに従って余接平面の正規直交基底 $e^1,e^2$ で $*e^1=e^2,\ *e^2=-e^1$ とする。従って $*\dd f$ は $\dd f$ と直交して同じ長さを持つ。また座標表示で微分すると

$$
d(*\dd f)=\partial_i(\sqrt{\det g}\,g^{ij}\partial_j f)\,\dd x^1\wedge\dd x^2
$$

である。例えば $\dd f=f_i\dd x^i$ に対し、$*\dd f$ は $g^{ij}f_j$ を成分に持つベクトルと面積形式 $\sqrt{\det g}\,\dd x^1\wedge\dd x^2$ の縮約なので、この式が成り立つ。

原点で微分が零でない解を得るため、半径 $\eps$ の円板を単位円板へ拡大して変数 $X=x/\eps$ を取る。拡大後の作用素を $L_\eps$ とすると、二階係数は $\delta_{ij}$ に近づき、一階係数は零に近づく。実際 $g_{ij}(\eps X)$ は滑らかで、$X$ 微分には $\eps$ が掛かる。単位円板上で $L_\eps w=-L_\eps X^1$、境界で $w=0$ を解く。[楕円型方程式の正則性](/math/pde/elliptic-regularity/) で扱う Dirichlet 問題の可解性と :ref[prop-pde-boundary-schauder]をここで使う。右辺の Hölder ノルムは零に収束する。楕円性の下界、係数の Hölder ノルム、円板を固定した境界評価により $\norm{w}_{C^{2,\alpha}}\le C\norm{L_\eps X^1}_{C^{0,\alpha}}$ であり、定数は小さな $\eps$ に一様に選べる。零次項がないため、評価の $\norm w_{C^0}$ は最大値原理と単位円板の二次の障壁で右辺により抑えられる。これらの偏微分方程式の一般結果は :cite[gilbarg-trudinger] を参照する。正則性により $w$ は内部で滑らかである。

よって $U=X^1+w$ は調和であり、原点で $\dd U\ne0$ となる。元の座標へ戻しても微分は非零で、十分小さい円板全体で非零にできる。$d(*\dd U)=0$ なので、円板上で $*\dd U=\dd V$ と書ける。$\dd U,\dd V$ が直交し同じ非零の長さなので、逆関数定理により $(U,V)$ は座標である。またその逆計量は内積 $\inner{\dd U}{\dd V}=0$ とノルムの一致により正のスカラー行列となる。従って計量も正のスカラー行列であり、等温座標を得る。
:::

:::theorem[局所 Weierstrass 表現]{#thm-cs-weierstrass}
単連結な複素平面の領域で、等温座標を持つはめ込み $X$ が極小であることと、$\Phi=X_u-iX_v$ が正則であることは同値である。このとき複素双線形内積について $\Phi\cdot\Phi=0$ である。逆に、どの点でも零ベクトルにならない正則な複素ベクトル関数 $\Phi$ がこの条件を満たすなら、$X=\Re\int\Phi\,\dd z$ は等温座標で表された極小はめ込みである。
:::

:::proof
等温座標では $X_{uu}+X_{vv}$ は法線に平行である。例えば $X_u$ との内積は $\tfrac12\partial_u\norm{X_u}^2+\partial_v\inner{X_v}{X_u}-\tfrac12\partial_u\norm{X_v}^2=0$ であり、$X_v$ についても同様である。法線成分は $h_{11}+h_{22}=2\lambda^2 H$ である。従って $H=0$ と調和性は同値であり、$\partial_{\bar z}\Phi=(X_{uu}+X_{vv})/2$ から正則性との同値を得る。正規直交性を倍率まで含めて書けば $\Phi\cdot\Phi=\norm{X_u}^2-\norm{X_v}^2-2i\inner{X_u}{X_v}=0$ となる。

逆に単連結領域では :ref[thm-cpx-primitive-period] と Cauchy の積分定理により積分が経路に依存しない。その実部を微分すると $X_u=\Re\Phi$、$X_v=-\Im\Phi$。零条件から直交し、長さは等しく、$\Phi\ne0$ により正である。調和性から $H=0$ を得る。
:::

$\Phi_1-i\Phi_2\ne0$ の領域では $f=(\Phi_1-i\Phi_2)/2$、$q=\Phi_3/(2f)$ と置くと

$$
\Phi=f(1-q^2,\,i(1+q^2),\,2q),\qquad
\lambda^2=|f|^2(1+|q|^2)^2.
\label{eq-cs-weierstrass-data}
$$

実際、零条件は $(\Phi_1-i\Phi_2)(\Phi_1+i\Phi_2)=-\Phi_3^2$ となり、これを解けばよい。$\Phi_1-i\Phi_2$ が恒等的に零なら零条件から $\Phi_3=0$ であり、$X$ は水平面の一部となる。恒等的に零でなければ $q$ は有理型関数として扱える。

有限の $q$ を持つ点では正則性の条件は $f\ne0$ である。$q$ が位数 $m$ の極を持つ点では、$fq^2$ が正則かつ非零になるように $f$ がちょうど位数 $2m$ の零を持つことが必要十分である。実際、三成分の正則性には $f,fq,fq^2$ の正則な延長が必要であり、極では $fq^2$ が最も低い次数を持つ。この項まで零になると三成分すべてが消える。例えば $f=z^2,q=1/z$ なら $\Phi=(z^2-1,i(z^2+1),2z)$ となり、原点でも極小はめ込みを与える。単連結でない領域では閉路ごとの積分の実部が零であることが、一価の写像を得る追加条件になる。

:::example[共役な極小曲面]
$a=1$ の懸垂面を $X=(\cosh v\cos u,\cosh v\sin u,v)$ とし、$Y=(-\sinh v\sin u,\sinh v\cos u,-u)$ とする。すると $Y_u=-X_v$、$Y_v=X_u$ であり、$Y$ は螺旋面を表す。$X_\theta=\cos\theta X+\sin\theta Y$ と置くと、一階微分が $(X_u,X_v)$ の二次元回転になり計量は同じである。また両者は調和なので $X_\theta$ は極小である。ここで同じ計量は局所の主張である。$u$ を $2\pi$ 増やすと $X$ は戻るが $Y$ の高さは変わるため、大域的には同じ曲面にならない。
:::

:::example[複素積分の周期が消えない場合]
穴を開けた領域 $\C\setminus\{0\}$ で $f=1/z,q=0$ と取ると、核は $\Phi=(1/z,i/z,0)$ であり、計量は $|z|^{-2}|\dd z|^2$ となる。局所的に $\int\Phi\,\dd z=(\log z,i\log z,0)$ と積分できる。ところが原点の周りを反時計回りに一周すると、積分の増分は $(2\pi i,-2\pi,0)$ で、実部に $-2\pi$ の移動が残る。従ってこの領域全体では一価の $X$ にならない。角度を実数値で選べる領域では $X=(\log|z|,-\arg z,0)$ となる。零条件で保証されるのは等温極小はめ込みの局所構成である。
:::

### 極小な線織面の分類

螺旋面には直線の母線があり、しかも極小である。一般の線織面で平均曲率を計算すると、母線を少し動かしたときの係数に強い制約がかかる。

:::theorem[Catalan の定理]{#thm-cs-catalan}
連結な正則線織面が極小なら、平面または螺旋面の開部分である。像が自己交差し得る局所表示の場合、結論は各局所像に対して成り立つ。
:::

:::proof
まず $d'\ne0$ の区間を取り、単位球面内の $d(u)$ をその弧長で表示する。$e=d'$、$n=d\times e$ と置けば、直交性の微分によりある関数 $\delta$ について $e'=-d+\delta n$、$n'=-\delta e$ となる。$X=c+vd$ で、$c$ に $b(u)d(u)$ を加えても母線は変わらない。$b=-\inner{c'}e$ と選ぶと新しい $c'$ の $e$ 成分が零になるので、$c'=\alpha d+\beta n$ と置ける。

この座標で

$$
X_u=\alpha d+ve+\beta n,\quad X_v=d,\quad X_u\times X_v=\beta e-vn.
$$

従って $G=1,F=\alpha$、ノルムの二乗は $v^2+\beta^2$。$X_{uv}=e$、$X_{vv}=0$ であり、$X_{uu}$ の $e,n$ 成分は順に $\alpha-\beta\delta,\ \beta'+v\delta$ となる。よって :eqref[eq-cs-curvature-invariants] の $H$ の分子に法線の分母を掛けた形は

$$
-\alpha\beta-\delta(\beta^2+v^2)-v\beta'
$$

である。ある $u$ に対し $v$ の開区間でこの多項式が零なら、係数比較で $\delta=0$、$\beta'=0$、$\alpha\beta=0$。$\beta\ne0$ なら $\alpha=0$、$n'=0$、$c'=\beta n$ となる。$d''=-d$ で $d$ は固定平面の単位円を速さ $1$ で回るので、剛体運動とパラメータの原点の選び方で $X=(v\cos u,v\sin u,\beta u)$ を得る。$\beta=0$ なら $n$ との内積が $c,d$ について一定または零なので平面の一部である。

次に $d'=0$ の点を考える。そこでは :eqref[eq-cs-ruled-curvature] から $K=0$。曲面の全域で $K=0$ なら $H=0$ と合わせて二主曲率が零になり、局所法線が一定であるため平面の一部となる。$K<0$ の点が一つでもあれば、その近傍で $d'\ne0$ かつ $\beta\ne0$ だから上で得た螺旋面の開部分がある。剛体運動後、その螺旋面は空間の実解析関数 $F(x,y,z)=x\sin(z/\beta)-y\cos(z/\beta)$ の零水準面である。極小面は等温座標で各成分が調和なので、局所的に実解析的なパラメータ表示を持つ。この性質と解析関数の一致の定理により、$F\circ X$ がある開近傍で恒等的に零なら、重なりを通して連結な曲面の全域で零になる。零水準面はどこでも $(F_x,F_y)\ne(0,0)$ なので正則であり、元の曲面はその開部分に含まれる。
:::

:::example[螺旋面の標構で係数を見る]
$X=(v\cos u,v\sin u,2u)$ では $d=(\cos u,\sin u,0)$、$e=(-\sin u,\cos u,0)$、$n=(0,0,1)$ であり、$\alpha=0,\beta=2,\delta=0$ である。上で得た $H$ の分子は零になる。さらに $v=0$ で $K=-1/4$ であり、二主曲率の絶対値はどちらも $1/2$ である。
:::

:::example[Scherk の極小グラフ]
$|u|,|v|<\pi/2$ で $q=\log(\cos v)-\log(\cos u)$ と置く。$q_u=\tan u$、$q_v=-\tan v$、$q_{uu}=\sec^2u$、$q_{vv}=-\sec^2v$、$q_{uv}=0$ である。:eqref[eq-cs-graph-curvature] に代入すると二つの項が相殺され $H=0$ となる。原点で $q$ の Hessian は $\diag(1,-1)$ だから $K=-1$。座標を一つ固定した切り口は対数曲線になる。この例で極小性を得た理由は、二つの二階微分の寄与が相殺されたことである。
:::

### 等温座標での Gauss 写像を読み取る

Weierstrass データの $q$ は、複素数を球面の法線方向へ対応させる座標に当たる。従って $f$ と $q$ の両方から位置を復元する一方、法線方向だけなら $q$ を読み取ればよい。

:::proposition[Weierstrass データから得る単位法線]
:ref[thm-cs-weierstrass] の規約で、有限の $q$ を取る領域では
$$
N=\frac{(2\Re q,2\Im q,|q|^2-1)}{1+|q|^2}
$$
である。
:::

:::proof
分子を $n$ とおく。$q=x+iy$ に対し、その第一、第二成分を $q,\bar q$ で書き直すと $2x=q+\bar q$、$2iy=q-\bar q$ である。従って $n$ と $(1-q^2,i(1+q^2),2q)$ の双線形内積は零になる。また $\norm n^2=(1+|q|^2)^2$ であるから、これは二つの可能な単位法線の一つである。$f$ を非零の複素数倍で変えても実部と虚部の作る向きは正の行列式で保たれる。$q=0,f=1$ のとき $X_u=(1,0,0),X_v=(0,-1,0)$ なので法線の第三成分は $-1$ である。$q$ を連続に動かしても外積は零にならないから、この符号がすべての有限な $q$ で同じである。
:::

:::example[法線が垂直から水平に移る]
$f=1,q=z$ では原点の法線は $(0,0,-1)$ である。$z=i$ では分子が $(0,2,0)$、分母が $2$ なので法線は $(0,1,0)$。さらに $z=2$ では $(4/5,0,3/5)$ となる。どの点でも計量の倍率 $(1+|z|^2)^2$ は正であり、法線が水平になることは座標の特異性を意味しない。
:::

:::corollary[Weierstrass データによる Gauss 曲率]{#cor-cs-weierstrass-curvature}
$q$ が有限で $f\ne0$ の領域では

$$
K=-\frac{4|q'(z)|^2}{|f(z)|^2(1+|q(z)|^2)^4}.
$$

特に Enneper 曲面では $f=1$、$q=z$ なので $K=-4/(1+|z|^2)^4$ である。
:::

:::proof
単位法線の式は $q$ の平面から単位球面への立体射影の逆写像である。その面積倍率の絶対値は $4/(1+|q|^2)^2$ であり、実部・虚部で微分した二本が直交し、それぞれ長さ $2/(1+|q|^2)$ を持つことから分かる。正則関数 $q$ の面積倍率は $|q'|^2$ である。従ってこれらの積を曲面の面積密度 $\lambda^2=|f|^2(1+|q|^2)^2$ で割ると $|K|$ を得る。$H=0$ だから二主曲率は互いに反対で $K\le0$ であり、負号を付ければ所要の式になる。
:::

:::example[Gauss 写像が一定となる平面の表示]
局所データを $f=1$、$q=1$ と置く。複素座標は $z=u+iv$ である。

まず複素双線形な内積で零となる三成分を作る。この式では複素共役を取らない。一方、計量に使う長さの二乗を計算する場面では共役が入る。

$$
\Phi=\left(0,2 i,2\right),\qquad \Phi\cdot\Phi=0,\qquad \lambda^2=4
$$

三成分を正則に積分し、その実部を取る。実部を取る前に積分することで、元の二つの実微分の関係を一度に満たす。

$$
\int\Phi\,\dd z=\left(0,2 i z,2 z\right),\qquad X(u,v)=\left(0,- 2 v,2 u\right)
$$

できた位置関数を実変数でも微分して確認する。二つの方向の長さは等しく、内積は零である。また調和性により平均曲率が零になる。

$$
g=\left(4\right)(\dd u^2+\dd v^2),\qquad X_{uu}+X_{vv}=0,\qquad \lambda^2\big|_{(u,v)=\left(0,0\right)}=4
$$

積分して得た位置は $X=(0,-2v,2u)$ で、像は平面である。零条件を満たす三成分の中に恒等的に零の成分があっても、ベクトル全体が消えなければよい。
:::

:::example[Gauss 写像に二次の変化を与える]
局所データを $f=1$、$q=z^{2}$ と置く。複素座標は $z=u+iv$ である。

まず複素双線形な内積で零となる三成分を作る。この式では複素共役を取らない。一方、計量に使う長さの二乗を計算する場面では共役が入る。

$$
\begin{aligned}
&\Phi=\left(- \left(z - 1\right) \left(z + 1\right) \left(z^{2} + 1\right),i \left(z^{4} + 1\right),2 z^{2}\right),\\
&\Phi\cdot\Phi=0,\\
&\lambda^2=\left(u^{4} + 2 u^{2} v^{2} + v^{4} + 1\right)^{2}
\end{aligned}
$$

三成分を正則に積分し、その実部を取る。実部を取る前に積分することで、元の二つの実微分の関係を一度に満たす。

$$
\begin{aligned}
&\int\Phi\,\dd z=\left(- \frac{z \left(z^{4} - 5\right)}{5},\frac{i z \left(z^{4} + 5\right)}{5},\frac{2 z^{3}}{3}\right),\\
&X(u,v)=\left(- \frac{u \left(u^{4} - 10 u^{2} v^{2} + 5 v^{4} - 5\right)}{5},- \frac{v \left(5 u^{4} - 10 u^{2} v^{2} + v^{4} + 5\right)}{5},\frac{2 u \left(u^{2} - 3 v^{2}\right)}{3}\right)
\end{aligned}
$$

できた位置関数を実変数でも微分して確認する。二つの方向の長さは等しく、内積は零である。また調和性により平均曲率が零になる。

$$
\begin{aligned}
&g=\left(\left(u^{4} + 2 u^{2} v^{2} + v^{4} + 1\right)^{2}\right)(\dd u^2+\dd v^2),\\
&X_{uu}+X_{vv}=0,\\
&\lambda^2\big|_{(u,v)=\left(1,0\right)}=4
\end{aligned}
$$

原点で $q=0$ だが、計量はそこで正である。Gauss 写像の微分が原点で消えることと、はめ込みの微分が消えることは別の条件である。
:::

:::example[係数に零がある場合の正則域]
局所データを $f=z$、$q=z$ と置く。複素座標は $z=u+iv$ である。

まず複素双線形な内積で零となる三成分を作る。この式では複素共役を取らない。一方、計量に使う長さの二乗を計算する場面では共役が入る。

$$
\begin{aligned}
&\Phi=\left(- z \left(z - 1\right) \left(z + 1\right),i z \left(z^{2} + 1\right),2 z^{2}\right),\\
&\Phi\cdot\Phi=0,\\
&\lambda^2=\left(u^{2} + v^{2}\right) \left(u^{2} + v^{2} + 1\right)^{2}
\end{aligned}
$$

三成分を正則に積分し、その実部を取る。実部を取る前に積分することで、元の二つの実微分の関係を一度に満たす。

$$
\begin{aligned}
&\int\Phi\,\dd z=\left(- \frac{z^{2} \left(z^{2} - 2\right)}{4},\frac{i z^{2} \left(z^{2} + 2\right)}{4},\frac{2 z^{3}}{3}\right),\\
&X(u,v)=\left(- \frac{u^{4} - 6 u^{2} v^{2} - 2 u^{2} + v^{4} + 2 v^{2}}{4},- u v \left(u^{2} - v^{2} + 1\right),\frac{2 u \left(u^{2} - 3 v^{2}\right)}{3}\right)
\end{aligned}
$$

できた位置関数を実変数でも微分して確認する。二つの方向の長さは等しく、内積は零である。また調和性により平均曲率が零になる。

$$
\begin{aligned}
&g=\left(\left(u^{2} + v^{2}\right) \left(u^{2} + v^{2} + 1\right)^{2}\right)(\dd u^2+\dd v^2),\\
&X_{uu}+X_{vv}=0,\\
&\lambda^2\big|_{(u,v)=\left(1,1\right)}=18
\end{aligned}
$$

原点では $f=q=0$ であり、三成分の微分が全て消えて計量も退化する。このデータから正則な極小曲面を得るには原点を除いた領域を使う。
:::

## 平行曲面と管状曲面

法線方向に同じ距離だけ移すと、元の曲面と移した曲面の法線は対応する点で平行になる。一方、その移動によって接方向がつぶれることがあり、その距離は主曲率から分かる。

:::proposition[平行曲面の曲率]{#prop-cs-parallel-surface}
:term[平行曲面]{en="parallel surface" yomi="へいこうきょくめん"} $X_d=X+dN$ は $1-dk_i\ne0$ の点で正則である。その点で同じ法線 $N$ を使うと主曲率は

$$
k_i(d)=\frac{k_i}{1-dk_i}.
\label{eq-cs-parallel-surface}
$$

面積倍率は $|(1-dk_1)(1-dk_2)|$ である。
:::

:::proof
接空間上で $dX_d=\id-dA$ となる。主方向について微分の倍率は $1-dk_i$ だから、階数 $2$ の条件と面積倍率を得る。像の接平面は元と平行で法線として $N$ を取れる。その微分は $dN$ と $(dX_d)^{-1}$ の合成なので、新形作用素は $A(\id-dA)^{-1}$ である。
:::

:::example[球の平行曲面]
外向き球では $k_i=-1/R$ であり、元の法線を運ぶと $k_i(d)=-1/(R+d)$ となる。例えば $R=3,d=-1$ なら半径 $2$ になり、面積倍率は $4/9$ である。$d=-3$ では全点が球心に移る。$d<-R$ では点が中心の反対側へ移るので、運んだ法線は新しい球の内向きであり、半径は $|R+d|$ である。
:::

:::proposition[管状曲面の曲率]{#prop-cs-tube-curvature}
正曲率の単位速さ空間曲線 $\alpha$ の各点の法平面に半径 $a>0$ の円を置き、:term[管状曲面]{en="tube surface" yomi="かんじょうきょくめん"} の表示を

$$
X(s,\theta)=\alpha(s)+a(\cos\theta\,N(s)+\sin\theta\,B(s))
$$

と定める。$1-a\kappa\cos\theta\ne0$ なら正則であり、断面の外向き単位法線 $n=\cos\theta N+\sin\theta B$ に対して主曲率は

$$
\frac{\kappa\cos\theta}{1-a\kappa\cos\theta},\qquad -\frac1a.
\label{eq-cs-tube-curvature}
$$

:::

:::proof
$e_\theta=-\sin\theta N+\cos\theta B$ とすると :eqref[eq-cs-space-frenet] より

$$
X_s=(1-a\kappa\cos\theta)T+a\tau e_\theta,\quad
X_\theta=ae_\theta,\quad n_s=-\kappa\cos\theta\,T+\tau e_\theta,\quad n_\theta=e_\theta.
$$

ここで $X_s-\tau X_\theta=(1-a\kappa\cos\theta)T$ であり、対応する法線の微分は $n_s-\tau n_\theta=-\kappa\cos\theta T$ である。従って $T$ と $e_\theta$ が主方向になり、式を得る。捩率は座標線のずれを作るがこの二つの曲率には直接入らない。
:::

:::example[円を芯とする管]
半径 $3$ の円を芯に $a=1$ とすれば全点で正則であり、内周 $\theta=0$ で主曲率は $1/2,-1$、外周 $\theta=\pi$ では $-1/4,-1$ になる。よって $K$ は内周で $-1/2$、外周で $1/4$。一般の芯では局所正則性があっても、離れた断面が交わらないためには別の条件が要る。
:::

正曲率の仮定は、この表示で Frenet 標構を使うためのものである。曲率が零になる芯にも :ref[prop-cs-bishop-frame] の平行標構を使って管を作れる。$n=\cos\theta\,E_1+\sin\theta\,E_2$、$e_\theta=-\sin\theta\,E_1+\cos\theta\,E_2$、$k=k_1\cos\theta+k_2\sin\theta$ と置く。ここで $k_1,k_2$ は芯の曲率ベクトルの二成分である。この標構では $n_s=-kT$、$n_\theta=e_\theta$ なので $X_s=(1-ak)T,\ X_\theta=ae_\theta$ となる。従って $1-ak\ne0$ が正則条件で、主曲率は $k/(1-ak),-1/a$ である。直線の芯では $k=0$ となり、通常の円柱を得る。

### 平行面の写像がつぶれる距離

形作用素の固有方向で微分すると、平行面への写像は接方向を $1-dk_i$ 倍する。従って一方が零になる距離では、元の曲面に正則表示があっても平行面の表示で階数が落ちる。どの方向がつぶれ、もう一方はどうなるかを調べると、曲面の中心や軸が平行面の極限として現れる理由が分かる。

:::example[一方向を失う円柱と二方向を失う球]
半径 $3$ の外向き円柱では主曲率は $-1/3,0$ であり、$d=-3$ とすると円周方向の倍率は零、母線方向の倍率は $1$ のままである。実際 $X=(3\cos u,3\sin u,v)$ に $-3N$ を加えると $(0,0,v)$ となり、像は軸になる。一方、半径 $3$ の球面では両曲率が $-1/3$ なので、同じ距離で両方向が零になる。中心を原点とすれば $X-3N$ は恒等的に原点である。平行面の曲率公式には分母の $1-dk_i$ があるため、これらの位置で代入して曲率を定義することはできない。
:::

:::example[円を芯とする管の内側]
半径 $4$ の円の周りに半径 $1$ の管を作る。内側の赤道を調べ、芯に垂直な断面の曲率と芯に沿う方向の曲率の符号を比べる。

座標で指定された点の位置と、速度に相当する二つの接ベクトルを計算する。正則性の計算は、二本が張る平行四辺形の面積が正であることを確かめることになる。

$$
\begin{aligned}
&X(u,v)=\left(\left(\cos{\left(u \right)} + 4\right) \cos{\left(v \right)},\left(\cos{\left(u \right)} + 4\right) \sin{\left(v \right)},\sin{\left(u \right)}\right),\\
&p=\left(3,0,0\right),\\
&(X_u,X_v)=\left(\left(0,0,-1\right),\left(0,3,0\right)\right)
\end{aligned}
$$

これらを材料に内積の表と外積を出す。以下で第二基本形式を計算する際には、外積の向きと指定された法線の向きを合わせる。符号が必要な場面では法線を規格化してから代入する。

$$
\begin{aligned}
&g(p)=\begin{pmatrix}1 & 0\\0 & 9\end{pmatrix},\\
&X_u\times X_v=\left(3,0,0\right),\\
&\sqrt{\det g}=3,\\
&N=\left(1,0,0\right)
\end{aligned}
$$

二階微分をする時点ではまだ座標値を固定しない。その後で二階微分を法線に射影すると、加速度のうち接平面内で向きを変える寄与が消える。

$$
\begin{aligned}X_{uu}&=\left(1,0,0\right),& X_{uv}&=\left(0,0,0\right),& X_{vv}&=\left(-3,0,0\right),\\ h(p)&=\begin{pmatrix}1 & 0\\0 & -3\end{pmatrix}, & A(p)&=\begin{pmatrix}1 & 0\\0 & - \frac{1}{3}\end{pmatrix}.\end{aligned}
$$

Gauss 曲率は二つの主曲率の積として行列式で、平均曲率は和の半分として跡で求める。主曲率自体が必要なら、この二つを係数とする二次式を解ける。

$$
\begin{aligned}
&K=\frac{\det h}{\det g}=- \frac{1}{3},\\
&H=\frac12\tr A=\frac{1}{3},\\
&\det(zI-A)=\frac{\left(z - 1\right) \left(3 z + 1\right)}{3}
\end{aligned}
$$

ここでの法線は管の内向きである。法線を反転すると両主曲率が反転するが、内側の赤道が負の Gauss 曲率を持つという結論は変わらない。
:::

:::example[平行面で円柱の半径が縮む]
外向き法線を持つ半径 $3$ の円柱を距離 $-1$ だけ移した後の表示を計算する。

座標で指定された点の位置と、速度に相当する二つの接ベクトルを計算する。正則性の計算は、二本が張る平行四辺形の面積が正であることを確かめることになる。

$$
\begin{aligned}
&X(u,v)=\left(2 \cos{\left(u \right)},2 \sin{\left(u \right)},v\right),\\
&p=\left(2,0,0\right),\\
&(X_u,X_v)=\left(\left(0,2,0\right),\left(0,0,1\right)\right)
\end{aligned}
$$

これらを材料に内積の表と外積を出す。以下で第二基本形式を計算する際には、外積の向きと指定された法線の向きを合わせる。符号が必要な場面では法線を規格化してから代入する。

$$
\begin{aligned}
&g(p)=\begin{pmatrix}4 & 0\\0 & 1\end{pmatrix},\\
&X_u\times X_v=\left(2,0,0\right),\\
&\sqrt{\det g}=2,\\
&N=\left(1,0,0\right)
\end{aligned}
$$

二階微分をする時点ではまだ座標値を固定しない。その後で二階微分を法線に射影すると、加速度のうち接平面内で向きを変える寄与が消える。

$$
\begin{aligned}X_{uu}&=\left(-2,0,0\right),& X_{uv}&=\left(0,0,0\right),& X_{vv}&=\left(0,0,0\right),\\ h(p)&=\begin{pmatrix}-2 & 0\\0 & 0\end{pmatrix}, & A(p)&=\begin{pmatrix}- \frac{1}{2} & 0\\0 & 0\end{pmatrix}.\end{aligned}
$$

Gauss 曲率は二つの主曲率の積として行列式で、平均曲率は和の半分として跡で求める。主曲率自体が必要なら、この二つを係数とする二次式を解ける。

$$
K=\frac{\det h}{\det g}=0,\qquad H=\frac12\tr A=- \frac{1}{4},\qquad \det(zI-A)=\frac{z \left(2 z + 1\right)}{2}
$$

元の円周方向の主曲率 $-1/3$ は $(-1/3)/(1-(-1)(-1/3))=-1/2$ になる。母線方向は零のままであり、表示から直接計算した結果と一致する。
:::

## 演習問題

::::exercise{level=basic}
外側の半径が $6$、内側の半径が $2$ の標準トーラスの両赤道の $K$ を求めよ。

:::solution
子午線の円の中心は回転軸から $4$、円の半径は $2$ である。$\mathrm I$ の弧長変数で $r=4+2\cos(s/2)$ とすれば $r''=-\tfrac12\cos(s/2)$。従って外側は $K=1/12$、内側は $K=-1/4$ となる。
:::
::::

::::exercise{level=basic}
$X(u,v)=(u,v,uv)$ が線織面であることを示し、原点で $K$ を求めよ。

:::solution
$u$ を固定すると $v$ の一次式なので直線である。原点の計量は単位行列、第二基本形式は混合項 $f=1$ だけを持つため行列式は $-1$。よって $K=-1$。
:::
::::

::::exercise{level=basic}
半径 $4$ の円柱を外向きに距離 $1$ だけ移した面の曲率を求めよ。

:::solution
主曲率は $-1/4,0$ から $(-1/4)/(1+1/4)=-1/5,0$ に移る。従って $K=0,H=-1/10$ である。
:::
::::

::::exercise{level=basic}
$a=2$ の懸垂面の最小半径の場所で $K$ と $H$ を求めよ。

:::solution
最小半径は $v=0$ で $a=2$ である。そこで主曲率は符号を除き $1/2$ で互いに反対なので $K=-1/4$、$H=0$。
:::
::::

::::exercise{level=advanced}
$X(u,v)=(u,v,p(u))$ が極小なら平面であることを示せ。

:::solution
極小グラフの方程式に $f(u,v)=p(u)$ を代入すると $\partial_u(p'/\sqrt{1+p'^2})=p''/(1+p'^2)^{3/2}=0$。従って $p''=0$ となり $p$ は一次式である。
:::
::::

::::exercise{level=advanced}
Weierstrass データ $f=2,q=0$ が定める曲面と計量を求めよ。

:::solution
$\Phi=(2,2i,0)$ だから $X=(2u,-2v,0)$ と取れる。二本の座標ベクトルは長さ $2$ で直交し、計量は $4(\dd u^2+\dd v^2)$。像は平面である。
:::
::::

::::exercise{level=advanced}
半径 $5$ の円を芯とする半径 $2$ の管で、$\cos\theta=0$ の点の主曲率を求めよ。

:::solution
芯の曲率は $1/5$ だが :eqref[eq-cs-tube-curvature] の第一の分子は零である。主曲率は $0,-1/2$ であり、従って $K=0,H=-1/4$ となる。
:::
::::

::::exercise{level=advanced}
$H=0$ で $K<0$ の点の平行曲面における平均曲率を、距離 $d$ と $k=\sqrt{-K}$ で表せ。

:::solution
元の二主曲率は $k,-k$ である。平行曲面で和を取ると $k/(1-dk)-k/(1+dk)=2dk^2/(1-d^2k^2)$。従って平均曲率は $dk^2/(1-d^2k^2)$ であり、$|d|=1/k$ では正則性が失われる。
:::
::::

## 確認テスト

::::quiz
線織面の正則点で常に成り立つ条件を選べ。

- [x] $K\le0$
- [ ] $K=0$
- [ ] $H=0$

:::explanation
母線は直線なので第二基本形式の母線方向成分は零であり、行列式は混合成分の二乗に負号を付けたものである。従って $K\le0$。螺旋面では $K<0$、円柱では $H\ne0$ なので後二者は一般には成立しない。
:::
::::

::::quiz{answer="2" tol="1e-8"}
等温座標の点で $\lambda^2=4$ とする。長さ $\norm{X_u}$ を求めよ。

:::explanation
等温計量の対角成分は $\inner{X_u}{X_u}=\lambda^2=4$ である。正の平方根の $2$ が長さである。
:::
::::

::::quiz
等温座標を持つ曲面が極小となる条件を選べ。

- [x] $X_{uu}+X_{vv}=0$
- [ ] $X_{uv}=0$
- [ ] $\lambda$ が一定

:::explanation
等温座標で $X_{uu}+X_{vv}$ の接成分は相殺し、法線成分は $2\lambda^2H$ である。従って極小条件は三成分の調和性に等しい。混合項や計量の倍率が一定である必要はない。
:::
::::

::::quiz{answer="4" tol="1e-8"}
半径 $2$ の外向き球面を距離 $2$ だけ外に移す。新しい半径を求めよ。

:::explanation
$N=X/2$ なので $X+2N=2X$。従って新しい半径は $2\cdot2=4$ である。主曲率の式にも $d=2,k=-1/2$ を代入すれば $-1/4$ が得られる。
:::
::::
